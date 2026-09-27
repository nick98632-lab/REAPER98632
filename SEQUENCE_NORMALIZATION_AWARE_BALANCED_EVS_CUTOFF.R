#!/usr/bin/env Rscript

# =============================================================================
# SEQUENCE NORMALIZATION-AWARE BALANCED EVS CUTOFF
# =============================================================================
# Purpose
#   Derive an EVS cutoff k* separately for each normalization state and each
#   RT/ZT comparison without using downstream significance calls.
#
# Principle
#   EVS should preserve PC1-defined structure while making the Leading Edge
#   (LE) and Remainder (REM) as comparable as possible in abundance and
#   mean-adjusted negative-binomial dispersion.
#
#   For every candidate k:
#     1) select top-k PASs independently in control and treatment by |PC1 loading|;
#     2) LE(k) = union(top-k control, top-k treatment); REM(k) = complement;
#     3) quantify abundance imbalance between LE and REM;
#     4) quantify dispersion-residual imbalance after adjusting dispersion for mean;
#     5) quantify spread imbalance of the dispersion residuals;
#     6) quantify fraction of total PC1-attributable variance captured by LE;
#     7) retain the Pareto frontier: maximize PC1 capture, minimize balance cost;
#     8) choose k* as the frontier point with maximum perpendicular deviation
#        from the endpoint chord (objective, no tuned lambda).
#
# Normalization states
#   RAW            : raw counts for PCA
#   DESEQ2_NORM    : DESeq2 median-of-ratios normalized counts for PCA
#   CPM_LOG1P      : log1p(CPM) for PCA
#   VST            : DESeq2 variance-stabilizing transformation for PCA
#
# Outputs
#   Tables/
#     Table_1_Optimal_Cutoffs.csv
#     Table_2_All_Candidate_k.csv
#     Table_3_Pareto_Frontiers.csv
#     Table_4_Selected_LE_REM_Balance.csv
#     Table_5_PAS_Level_Selected_Membership.csv
#     Normalization_Aware_EVS_Cutoff_Manifest.csv
#   Figures/
#     Figure_1_Cutoff_Summary.pdf/png
#     Figure_2_Balance_vs_k_<METHOD>.pdf/png
#     Figure_3_PC1_vs_Balance_Pareto_<METHOD>.pdf/png
#     Figure_4_Selected_Mean_Dispersion_<METHOD>.pdf/png
#     Figure_5_Selected_LE_REM_Distributions_<METHOD>.pdf/png
#     Figure_6_Cross_Normalization_Comparison.pdf/png
#
# Notes
#   * Downstream HBFSS/TWAS significance is NOT used to choose k*.
#   * DESeq2 raw-count NB dispersions are used only to quantify expected
#     mean-dependent count variance and dispersion residuals.
# =============================================================================

options(stringsAsFactors = FALSE)
set.seed(98632)

required <- c("DESeq2", "SummarizedExperiment", "ggplot2", "dplyr", "tidyr", "gridExtra", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(scales)
})

# -----------------------------------------------------------------------------
# Paths and sample definitions
# -----------------------------------------------------------------------------
repo_candidates <- unique(c(
  Sys.getenv("SEQUENCE_REPO_ROOT", unset = ""),
  getwd(), dirname(getwd()), "/root/REAPER98632"
))
repo_candidates <- repo_candidates[nzchar(repo_candidates)]
repo_root <- NULL
for (p in repo_candidates) {
  if (file.exists(file.path(p, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")) ||
      file.exists(file.path(p, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"))) {
    repo_root <- normalizePath(p, winslash = "/", mustWork = TRUE)
    break
  }
}
if (is.null(repo_root)) repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

count_path <- if (file.exists(file.path(repo_root, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"))) {
  file.path(repo_root, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
} else {
  file.path(repo_root, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
}
if (!file.exists(count_path)) stop("Count file not found: ", count_path)

OUT_ROOT <- Sys.getenv(
  "SEQUENCE_BALANCED_EVS_OUT",
  unset = file.path(repo_root, "exports", "Normalization_Aware_Balanced_EVS")
)
FIG_DIR <- file.path(OUT_ROOT, "Figures")
TAB_DIR <- file.path(OUT_ROOT, "Tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

GROUP_PATTERNS <- c(
  RT0="^R0_", ZT6="^ZT6_", RT2="^R2_", ZT8="^ZT8_",
  RT4="^R4_", ZT10="^ZT10_", RT8="^R8_", ZT14="^ZT14_"
)
COMPARISONS <- list(
  RT0_ZT6=c(control="RT0", treatment="ZT6"),
  RT2_ZT8=c(control="RT2", treatment="ZT8"),
  RT4_ZT10=c(control="RT4", treatment="ZT10"),
  RT8_ZT14=c(control="RT8", treatment="ZT14")
)
METHODS <- c("RAW", "DESEQ2_NORM", "CPM_LOG1P", "VST")

# Minimum number of PASs required on BOTH sides of the split.
# This is not a tuned cutoff; it only prevents degenerate LE/REM groups.
MIN_SIDE_N <- as.integer(Sys.getenv("SEQUENCE_BALANCED_EVS_MIN_SIDE_N", unset = "250"))

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
read_counts <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  idx <- sort(unique(unlist(lapply(GROUP_PATTERNS, function(p) grep(p, names(x))))))
  if (!length(idx)) stop("No sample columns matched GROUP_PATTERNS.")
  ids <- trimws(as.character(x[[1]]))
  blank <- is.na(ids) | ids == ""
  ids[blank] <- paste0("feature_", which(blank))
  ids <- make.unique(ids, sep="__dup_")
  m <- do.call(cbind, lapply(x[, idx, drop=FALSE], function(z) suppressWarnings(as.numeric(as.character(z)))))
  rownames(m) <- ids
  colnames(m) <- names(x)[idx]
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  m <- pmax(m, 0)
  m <- m[rowSums(m) > 0, , drop=FALSE]
  m
}

assign_groups <- function(samples) {
  g <- rep(NA_character_, length(samples)); names(g) <- samples
  for (nm in names(GROUP_PATTERNS)) {
    hit <- grep(GROUP_PATTERNS[[nm]], samples)
    if (any(!is.na(g[hit]))) stop("Sample matched multiple groups.")
    g[hit] <- nm
  }
  if (anyNA(g)) stop("Unassigned samples: ", paste(names(g)[is.na(g)], collapse=", "))
  factor(g, levels = names(GROUP_PATTERNS))
}

log1p_cpm <- function(x) {
  lib <- colSums(x)
  lib[!is.finite(lib) | lib <= 0] <- 1
  log1p(sweep(x, 2, lib/1e6, "/"))
}

safe_sd <- function(x) {
  s <- stats::sd(x, na.rm=TRUE)
  if (!is.finite(s) || s <= 1e-12) 1 else s
}

minmax01 <- function(x) {
  r <- range(x, finite=TRUE)
  if (!all(is.finite(r)) || diff(r) <= .Machine$double.eps) return(rep(0, length(x)))
  (x-r[1])/diff(r)
}

pc1_info <- function(mat) {
  # Samples are observations; PASs are variables/features.
  p <- stats::prcomp(t(mat), center=TRUE, scale.=FALSE, rank.=1)
  load <- p$rotation[,1]
  load[!is.finite(load)] <- 0
  contrib <- (p$sdev[1]^2) * load^2
  ord_desc <- order(abs(load), decreasing=TRUE)
  rank_desc <- integer(length(load)); rank_desc[ord_desc] <- seq_along(ord_desc)
  names(rank_desc) <- names(load)
  list(load=load, contribution=contrib, rank=rank_desc, lambda1=p$sdev[1]^2)
}

fit_pair_nb <- function(raw_pair, condition) {
  # Independent NB fit for the pair. Dispersions are estimated from RAW counts.
  keep <- rowSums(raw_pair) > 0
  raw_pair <- raw_pair[keep,,drop=FALSE]
  condition <- droplevels(condition)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(raw_pair),
    colData = data.frame(condition=condition, row.names=colnames(raw_pair)),
    design = ~ condition
  )
  dds <- tryCatch(DESeq2::estimateSizeFactors(dds),
                  error=function(e) DESeq2::estimateSizeFactors(dds, type="poscounts"))
  dds <- DESeq2::estimateDispersions(dds, quiet=TRUE)
  norm <- DESeq2::counts(dds, normalized=TRUE)
  base_mean <- rowMeans(norm)
  alpha <- DESeq2::dispersions(dds)
  alpha[!is.finite(alpha) | alpha <= 0] <- NA_real_
  trend_fun <- DESeq2::dispersionFunction(dds)
  alpha_trend <- suppressWarnings(trend_fun(base_mean))
  alpha_trend[!is.finite(alpha_trend) | alpha_trend <= 0] <- NA_real_
  eps <- 1e-12
  disp_resid <- log((alpha + eps)/(alpha_trend + eps))
  disp_resid[!is.finite(disp_resid)] <- 0

  # Expected NB variance on normalized-count scale. We use the fitted dispersion
  # and mean; this prevents high means from being misread as excess dispersion.
  nb_var <- base_mean + pmax(alpha, eps) * base_mean^2
  nb_var[!is.finite(nb_var)] <- 0

  list(dds=dds, norm=norm, base_mean=base_mean, alpha=alpha,
       alpha_trend=alpha_trend, disp_resid=disp_resid, nb_var=nb_var)
}

build_method_matrix <- function(method, raw_pair, nb_obj) {
  if (method == "RAW") return(raw_pair)
  if (method == "DESEQ2_NORM") return(nb_obj$norm)
  if (method == "CPM_LOG1P") return(log1p_cpm(raw_pair))
  if (method == "VST") {
    vst_obj <- DESeq2::varianceStabilizingTransformation(nb_obj$dds, blind=TRUE)
    return(SummarizedExperiment::assay(vst_obj))
  }
  stop("Unknown method: ", method)
}

method_abundance_axis <- function(method, method_matrix, nb_obj) {
  # Abundance balance is evaluated on the scale seen by PCA for that method.
  # RAW/DESEQ2_NORM are log1p-transformed only for the balance metric so that
  # extreme counts do not dominate the standardized difference.
  if (method %in% c("RAW", "DESEQ2_NORM")) return(log1p(rowMeans(method_matrix)))
  rowMeans(method_matrix)
}

# Exact candidate-k scan using entry rank = min(rank_control, rank_treatment).
# LE(k) is therefore every PAS whose first entry into either arm top-k occurs by k.
scan_exact_k <- function(method, comparison_name, raw_pair, condition, nb_obj, method_matrix) {
  levs <- levels(condition)
  if (length(levs) != 2L) stop("Pairwise comparison must have two conditions.")
  ctrl <- levs[1]; trt <- levs[2]
  idx_c <- which(condition == ctrl)
  idx_t <- which(condition == trt)

  pc_c <- pc1_info(method_matrix[,idx_c,drop=FALSE])
  pc_t <- pc1_info(method_matrix[,idx_t,drop=FALSE])

  ids <- rownames(method_matrix)
  rc <- pc_c$rank[ids]
  rt <- pc_t$rank[ids]
  entry_rank <- pmin(rc, rt)
  pc_mass <- pc_c$contribution[ids] + pc_t$contribution[ids]
  pc_mass[!is.finite(pc_mass)] <- 0

  abundance <- method_abundance_axis(method, method_matrix, nb_obj)[ids]
  disp_resid <- nb_obj$disp_resid[ids]
  nb_var <- nb_obj$nb_var[ids]
  base_mean <- nb_obj$base_mean[ids]
  alpha <- nb_obj$alpha[ids]
  alpha_trend <- nb_obj$alpha_trend[ids]

  # Sort by first k at which a PAS enters LE. Ties enter together.
  o <- order(entry_rank, ids)
  entry_rank <- entry_rank[o]
  ids_o <- ids[o]
  pc_o <- pc_mass[o]
  ab_o <- abundance[o]
  dr_o <- disp_resid[o]

  # cumulative prefix summaries
  n <- length(ids_o)
  c_pc <- cumsum(pc_o)
  total_pc <- sum(pc_o)
  if (!is.finite(total_pc) || total_pc <= 0) stop("No finite PC1 contribution.")

  c_ab <- cumsum(ab_o); c_ab2 <- cumsum(ab_o^2)
  c_dr <- cumsum(dr_o); c_dr2 <- cumsum(dr_o^2)
  total_ab <- sum(ab_o); total_ab2 <- sum(ab_o^2)
  total_dr <- sum(dr_o); total_dr2 <- sum(dr_o^2)

  # Every integer k is considered, but only k values where LE membership changes
  # need evaluation. This is exact because membership is constant between ties.
  candidate_k <- sort(unique(entry_rank))
  candidate_k <- candidate_k[candidate_k >= 1L]

  rows <- vector("list", length(candidate_k))
  for (ii in seq_along(candidate_k)) {
    k <- candidate_k[ii]
    m <- max(which(entry_rank <= k))
    n_le <- m
    n_rem <- n - m
    if (n_le < MIN_SIDE_N || n_rem < MIN_SIDE_N) next

    # means
    mean_ab_le <- c_ab[m]/n_le
    mean_ab_rem <- (total_ab-c_ab[m])/n_rem
    mean_dr_le <- c_dr[m]/n_le
    mean_dr_rem <- (total_dr-c_dr[m])/n_rem

    # unbiased variances
    var_ab_le <- if (n_le > 1) (c_ab2[m] - c_ab[m]^2/n_le)/(n_le-1) else 0
    rem_ab_sum <- total_ab-c_ab[m]; rem_ab2_sum <- total_ab2-c_ab2[m]
    var_ab_rem <- if (n_rem > 1) (rem_ab2_sum - rem_ab_sum^2/n_rem)/(n_rem-1) else 0

    var_dr_le <- if (n_le > 1) (c_dr2[m] - c_dr[m]^2/n_le)/(n_le-1) else 0
    rem_dr_sum <- total_dr-c_dr[m]; rem_dr2_sum <- total_dr2-c_dr2[m]
    var_dr_rem <- if (n_rem > 1) (rem_dr2_sum - rem_dr_sum^2/n_rem)/(n_rem-1) else 0

    pooled_ab <- sqrt(max(((n_le-1)*var_ab_le + (n_rem-1)*var_ab_rem)/(n-2), 1e-12))
    pooled_dr <- sqrt(max(((n_le-1)*var_dr_le + (n_rem-1)*var_dr_rem)/(n-2), 1e-12))

    abundance_smd <- abs(mean_ab_le-mean_ab_rem)/pooled_ab
    dispersion_smd <- abs(mean_dr_le-mean_dr_rem)/pooled_dr
    dispersion_spread <- abs(log((var_dr_le+1e-12)/(var_dr_rem+1e-12)))
    pc1_capture <- c_pc[m]/total_pc

    rows[[ii]] <- data.frame(
      method=method,
      comparison=comparison_name,
      k=as.integer(k),
      LE_n=n_le,
      REM_n=n_rem,
      PC1_capture=pc1_capture,
      abundance_SMD=abundance_smd,
      dispersion_residual_SMD=dispersion_smd,
      dispersion_residual_log_variance_ratio=dispersion_spread,
      abundance_mean_LE=mean_ab_le,
      abundance_mean_REM=mean_ab_rem,
      dispersion_residual_mean_LE=mean_dr_le,
      dispersion_residual_mean_REM=mean_dr_rem,
      dispersion_residual_var_LE=var_dr_le,
      dispersion_residual_var_REM=var_dr_rem,
      stringsAsFactors=FALSE
    )
  }
  scan <- dplyr::bind_rows(rows)
  if (!nrow(scan)) stop("No valid candidate k remained after MIN_SIDE_N constraint.")

  # Scale each imbalance metric across the exact candidate set so no one metric
  # dominates solely because of units/range. Equal weighting is applied only to
  # the three imbalance dimensions; the final cutoff uses Pareto geometry, not a
  # tuned signal-vs-cost weight.
  scan <- scan %>%
    mutate(
      abundance_cost01 = minmax01(abundance_SMD),
      dispersion_location_cost01 = minmax01(dispersion_residual_SMD),
      dispersion_spread_cost01 = minmax01(dispersion_residual_log_variance_ratio),
      balance_cost = (abundance_cost01 + dispersion_location_cost01 + dispersion_spread_cost01)/3
    )

  # Pareto frontier: maximize PC1_capture, minimize balance_cost.
  tmp <- scan %>% arrange(balance_cost, desc(PC1_capture), desc(k))
  tmp <- tmp[!duplicated(signif(tmp$balance_cost, 12)), , drop=FALSE]
  prior_best <- c(-Inf, head(cummax(tmp$PC1_capture), -1L))
  frontier <- tmp[tmp$PC1_capture > prior_best, , drop=FALSE] %>%
    arrange(balance_cost, PC1_capture, k)
  if (nrow(frontier) < 2L) stop("Pareto frontier has fewer than 2 points.")

  # Endpoint-chord maximum deviation in normalized objective space.
  x <- minmax01(frontier$balance_cost)
  y <- minmax01(frontier$PC1_capture)
  x1 <- x[1]; y1 <- y[1]
  x2 <- x[length(x)]; y2 <- y[length(y)]
  denom <- sqrt((y2-y1)^2 + (x2-x1)^2)
  if (!is.finite(denom) || denom <= .Machine$double.eps) stop("Undefined Pareto endpoint chord.")
  frontier$endpoint_deviation <- abs((y2-y1)*x - (x2-x1)*y + x2*y1 - y2*x1)/denom

  best_dev <- max(frontier$endpoint_deviation, na.rm=TRUE)
  cand <- frontier[abs(frontier$endpoint_deviation-best_dev) < 1e-12,,drop=FALSE]
  # Tie-breaker: higher PC1 capture, then lower balance cost, then smaller k.
  cand <- cand[order(-cand$PC1_capture, cand$balance_cost, cand$k),,drop=FALSE]
  selected <- cand[1,,drop=FALSE]

  kstar <- selected$k[1]
  le <- entry_rank <= kstar
  membership <- data.frame(
    feature_id=ids_o,
    method=method,
    comparison=comparison_name,
    control_rank=rc[ids_o],
    treatment_rank=rt[ids_o],
    entry_rank=entry_rank,
    PC1_mass=pc_o,
    abundance_axis=ab_o,
    baseMean=base_mean[ids_o],
    dispersion=alpha[ids_o],
    dispersion_trend=alpha_trend[ids_o],
    dispersion_residual=dr_o,
    expected_NB_variance=nb_var[ids_o],
    selected_k=kstar,
    dataset=ifelse(le, "Leading_Edge", "Remainder"),
    stringsAsFactors=FALSE
  )

  list(scan=scan, frontier=frontier, selected=selected,
       membership=membership, pc_control=pc_c, pc_treatment=pc_t)
}

save_both <- function(p, stem, width=10, height=7) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p, width=width, height=height, device=cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p, width=width, height=height, dpi=320)
}

# -----------------------------------------------------------------------------
# Load counts and run all normalization x comparison analyses
# -----------------------------------------------------------------------------
message("Reading counts: ", count_path)
counts <- read_counts(count_path)
groups <- assign_groups(colnames(counts))
message("PASs: ", nrow(counts), " | samples: ", ncol(counts))

all_scans <- list(); all_frontiers <- list(); all_selected <- list(); all_membership <- list()

for (cmp_name in names(COMPARISONS)) {
  mp <- COMPARISONS[[cmp_name]]
  keep_cols <- groups %in% c(mp[["control"]], mp[["treatment"]])
  raw_pair <- counts[,keep_cols,drop=FALSE]
  cond_chr <- as.character(groups[keep_cols])
  # Force control first, treatment second.
  condition <- factor(cond_chr, levels=c(mp[["control"]], mp[["treatment"]]))

  # Remove PASs zero across this comparison only.
  keep_rows <- rowSums(raw_pair) > 0
  raw_pair <- raw_pair[keep_rows,,drop=FALSE]

  message("\n", cmp_name, ": fitting NB model on ", nrow(raw_pair), " PASs...")
  nb_obj <- fit_pair_nb(raw_pair, condition)

  for (method in METHODS) {
    message("  ", method, ": scanning exact EVS k...")
    mm <- build_method_matrix(method, raw_pair, nb_obj)
    res <- scan_exact_k(method, cmp_name, raw_pair, condition, nb_obj, mm)
    key <- paste(method, cmp_name, sep="__")
    all_scans[[key]] <- res$scan
    all_frontiers[[key]] <- transform(res$frontier, method=method, comparison=cmp_name)
    sel <- res$selected
    sel$method <- method; sel$comparison <- cmp_name
    all_selected[[key]] <- sel
    all_membership[[key]] <- res$membership
    message("    k* = ", sel$k[1],
            " | PC1 capture = ", sprintf("%.3f", sel$PC1_capture[1]),
            " | balance cost = ", sprintf("%.3f", sel$balance_cost[1]))
  }
}

scan_df <- bind_rows(all_scans)
frontier_df <- bind_rows(all_frontiers)
selected_df <- bind_rows(all_selected) %>%
  select(method, comparison, everything()) %>%
  arrange(method, comparison)
membership_df <- bind_rows(all_membership)

# Selected LE/REM balance summary
balance_df <- selected_df %>%
  transmute(
    method, comparison, selected_k=k, LE_n, REM_n,
    PC1_capture, balance_cost,
    abundance_SMD,
    dispersion_residual_SMD,
    dispersion_residual_log_variance_ratio,
    abundance_mean_LE, abundance_mean_REM,
    dispersion_residual_mean_LE, dispersion_residual_mean_REM,
    dispersion_residual_var_LE, dispersion_residual_var_REM,
    endpoint_deviation
  )

# Manifest format convenient for later pipeline integration.
manifest <- selected_df %>%
  select(method, comparison, k) %>%
  mutate(method = recode(method,
                         RAW="Balanced_RAW",
                         DESEQ2_NORM="Balanced_DESeq2Norm",
                         CPM_LOG1P="Balanced_CPM",
                         VST="Balanced_VST")) %>%
  pivot_wider(names_from=comparison, values_from=k, names_sort=FALSE) %>%
  rename(Cutoff_Method=method)
manifest$Source <- "Normalization-aware balanced EVS: Pareto maximum endpoint-chord deviation, maximizing PC1 capture while minimizing LE/REM abundance and mean-adjusted dispersion imbalance"

write.csv(balance_df, file.path(TAB_DIR, "Table_1_Optimal_Cutoffs.csv"), row.names=FALSE)
write.csv(scan_df, file.path(TAB_DIR, "Table_2_All_Candidate_k.csv"), row.names=FALSE)
write.csv(frontier_df, file.path(TAB_DIR, "Table_3_Pareto_Frontiers.csv"), row.names=FALSE)
write.csv(balance_df, file.path(TAB_DIR, "Table_4_Selected_LE_REM_Balance.csv"), row.names=FALSE)
write.csv(membership_df, file.path(TAB_DIR, "Table_5_PAS_Level_Selected_Membership.csv"), row.names=FALSE)
write.csv(manifest, file.path(TAB_DIR, "Normalization_Aware_EVS_Cutoff_Manifest.csv"), row.names=FALSE)
write.csv(manifest, file.path(OUT_ROOT, "Normalization_Aware_EVS_Cutoff_Manifest.csv"), row.names=FALSE)

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
method_labels <- c(
  RAW="Raw counts",
  DESEQ2_NORM="DESeq2 normalized",
  CPM_LOG1P="log1p(CPM)",
  VST="DESeq2 VST"
)

# Figure 1: selected cutoff summary
p1 <- balance_df %>%
  mutate(method=factor(method, levels=METHODS, labels=method_labels[METHODS])) %>%
  ggplot(aes(x=comparison, y=selected_k, group=method, shape=method)) +
  geom_line(position=position_dodge(width=0.25)) +
  geom_point(size=3, position=position_dodge(width=0.25)) +
  geom_hline(yintercept=5000, linetype=2) +
  labs(x=NULL, y="Selected EVS cutoff k*",
       title="Normalization-aware balanced EVS cutoffs",
       subtitle="Dashed line = fixed 5,000 comparator") +
  theme_bw(base_size=12) + theme(legend.title=element_blank())
save_both(p1, "Figure_1_Cutoff_Summary", 10, 6.5)

# Figure 2: balance metrics vs k, one multi-panel figure per method
for (m in METHODS) {
  d <- scan_df %>% filter(method==m) %>%
    select(method, comparison, k, abundance_cost01, dispersion_location_cost01,
           dispersion_spread_cost01, balance_cost) %>%
    pivot_longer(cols=c(abundance_cost01, dispersion_location_cost01,
                        dispersion_spread_cost01, balance_cost),
                 names_to="metric", values_to="value") %>%
    mutate(metric=recode(metric,
                         abundance_cost01="Mean abundance imbalance",
                         dispersion_location_cost01="Dispersion residual location imbalance",
                         dispersion_spread_cost01="Dispersion residual spread imbalance",
                         balance_cost="Combined balance cost"))
  sel <- balance_df %>% filter(method==m)
  p <- ggplot(d, aes(k, value, linetype=metric)) +
    geom_line(linewidth=0.65) +
    geom_vline(data=sel, aes(xintercept=selected_k), inherit.aes=FALSE, linetype=2) +
    facet_wrap(~comparison, scales="free_x", ncol=2) +
    labs(x="EVS cutoff k", y="Scaled imbalance / cost",
         title=paste0(method_labels[[m]], ": LE vs Remainder balance across k")) +
    theme_bw(base_size=11) + theme(legend.title=element_blank(), legend.position="bottom")
  save_both(p, paste0("Figure_2_Balance_vs_k_", m), 11, 8)
}

# Figure 3: Pareto PC1 capture vs balance cost
for (m in METHODS) {
  d <- frontier_df %>% filter(method==m)
  sel <- balance_df %>% filter(method==m)
  p <- ggplot(d, aes(balance_cost, PC1_capture)) +
    geom_path() + geom_point(aes(size=endpoint_deviation), alpha=0.55) +
    geom_point(data=sel, aes(balance_cost, PC1_capture), inherit.aes=FALSE, size=3.5, shape=8) +
    geom_text(data=sel, aes(balance_cost, PC1_capture, label=paste0("k* = ", selected_k)),
              inherit.aes=FALSE, nudge_y=0.015, size=3.3) +
    facet_wrap(~comparison, scales="free", ncol=2) +
    labs(x="LE/REM balance cost (lower is better)", y="PC1 variance captured by LE",
         title=paste0(method_labels[[m]], ": PC1 capture versus LE/REM balance"),
         subtitle="Star = maximum endpoint-chord deviation on Pareto frontier") +
    theme_bw(base_size=11) + theme(legend.position="none")
  save_both(p, paste0("Figure_3_PC1_vs_Balance_Pareto_", m), 11, 8)
}

# Figure 4: mean-dispersion landscape at selected k
for (m in METHODS) {
  d <- membership_df %>% filter(method==m) %>%
    mutate(dataset=factor(dataset, levels=c("Leading_Edge","Remainder"), labels=c("Leading Edge","Remainder")))
  p <- ggplot(d, aes(baseMean, dispersion, shape=dataset)) +
    geom_point(alpha=0.25, size=0.7) +
    geom_line(aes(y=dispersion_trend, group=comparison), linewidth=0.7, alpha=0.8) +
    scale_x_log10(labels=scales::label_number()) +
    scale_y_log10(labels=scales::label_number()) +
    facet_wrap(~comparison, scales="free", ncol=2) +
    labs(x="DESeq2 normalized mean", y="DESeq2 dispersion",
         title=paste0(method_labels[[m]], ": selected LE/Remainder mean-dispersion landscape"),
         subtitle="Curve = fitted mean-dispersion trend; membership comes from selected k*") +
    theme_bw(base_size=11) + theme(legend.title=element_blank(), legend.position="bottom")
  save_both(p, paste0("Figure_4_Selected_Mean_Dispersion_", m), 11, 8)
}

# Figure 5: selected LE/REM distributions in abundance and dispersion residual
for (m in METHODS) {
  d <- membership_df %>% filter(method==m) %>%
    mutate(dataset=factor(dataset, levels=c("Leading_Edge","Remainder"), labels=c("Leading Edge","Remainder")))
  d_long <- bind_rows(
    d %>% transmute(comparison, dataset, variable="Abundance axis", value=abundance_axis),
    d %>% transmute(comparison, dataset, variable="Mean-adjusted dispersion residual", value=dispersion_residual)
  )
  p <- ggplot(d_long, aes(value, linetype=dataset)) +
    geom_density(linewidth=0.8, adjust=1.1) +
    facet_grid(variable ~ comparison, scales="free") +
    labs(x=NULL, y="Density",
         title=paste0(method_labels[[m]], ": LE and Remainder distributions at selected k*")) +
    theme_bw(base_size=10.5) + theme(legend.title=element_blank(), legend.position="bottom")
  save_both(p, paste0("Figure_5_Selected_LE_REM_Distributions_", m), 12, 7.5)
}

# Figure 6: cross-normalization comparison of selected performance
cross <- balance_df %>%
  mutate(method=factor(method, levels=METHODS, labels=method_labels[METHODS])) %>%
  select(method, comparison, selected_k, PC1_capture, balance_cost,
         abundance_SMD, dispersion_residual_SMD,
         dispersion_residual_log_variance_ratio) %>%
  pivot_longer(cols=c(selected_k, PC1_capture, balance_cost, abundance_SMD,
                      dispersion_residual_SMD, dispersion_residual_log_variance_ratio),
               names_to="metric", values_to="value") %>%
  mutate(metric=recode(metric,
                       selected_k="Selected k*",
                       PC1_capture="PC1 capture",
                       balance_cost="Combined balance cost",
                       abundance_SMD="Abundance SMD",
                       dispersion_residual_SMD="Dispersion-residual SMD",
                       dispersion_residual_log_variance_ratio="Dispersion residual log variance ratio"))

p6 <- ggplot(cross, aes(method, value, shape=method)) +
  geom_point(size=2.7) +
  facet_grid(metric ~ comparison, scales="free_y") +
  labs(x=NULL, y=NULL, title="Cross-normalization comparison of balanced EVS solutions") +
  theme_bw(base_size=9.5) +
  theme(axis.text.x=element_text(angle=35, hjust=1), legend.position="none")
save_both(p6, "Figure_6_Cross_Normalization_Comparison", 13, 11)

# -----------------------------------------------------------------------------
# Manuscript-ready methods text
# -----------------------------------------------------------------------------
methods_text <- c(
  "# Normalization-aware balanced eigenvector-splitting cutoff",
  "",
  "For each RT/ZT comparison, the EVS cutoff was derived independently under four preprocessing states: raw counts, DESeq2 median-of-ratios normalized counts, log1p counts-per-million, and DESeq2 variance-stabilized counts. Within each condition, PCA was performed on the transposed PAS-by-sample matrix without feature scaling, and PASs were ranked by the absolute PC1 loading. At a candidate cutoff k, the Leading Edge was defined as the union of the condition-specific top-k PAS sets and the Remainder as all other PASs.",
  "",
  "The cutoff criterion was designed to preserve PC1-defined structure while minimizing systematic abundance and dispersion differences between the Leading Edge and Remainder. For every PAS, DESeq2 negative-binomial dispersions were estimated from the raw counts within the corresponding pairwise comparison. Mean-dependent dispersion was removed by comparing the fitted PAS dispersion with the DESeq2 mean-dispersion trend, using the log dispersion residual log(alpha_i / alpha_trend(mu_i)). Thus, higher variance expected solely from higher mean abundance was not treated as excess biological or technical variability.",
  "",
  "For every candidate k at which Leading-Edge membership changed, three scale-free imbalance quantities were calculated between the Leading Edge and Remainder: the standardized mean difference in abundance, the standardized mean difference in mean-adjusted dispersion residuals, and the absolute log ratio of dispersion-residual variances. Each imbalance quantity was min-max normalized over candidate k values and their mean defined the composite LE/Remainder balance cost. PC1 benefit was defined as the fraction of the combined control- and treatment-arm PC1-attributable variance, lambda1*v_i1^2, captured by the Leading Edge.",
  "",
  "Candidate cutoffs were evaluated on the two-objective plane of PC1 benefit versus LE/Remainder balance cost. Nondominated solutions formed a Pareto frontier. The final k* was selected as the Pareto point with the maximum perpendicular distance from the chord joining the two frontier endpoints. This avoided specifying an arbitrary signal-versus-balance weight and locked the cutoff without reference to downstream DESeq2, HBFSS, or TWAS significance results.",
  "",
  "The resulting normalization-specific k* values can subsequently be entered into the downstream SEQUENCE pipeline. Downstream significance counts are therefore used as validation of the cutoff rule rather than as the criterion used to select the cutoff."
)
writeLines(methods_text, file.path(OUT_ROOT, "METHODS_Normalization_Aware_Balanced_EVS.md"))

# Audit
writeLines(c(
  paste0("Count file: ", normalizePath(count_path, winslash="/", mustWork=TRUE)),
  paste0("PASs in experiment-wide matrix: ", nrow(counts)),
  paste0("Samples: ", ncol(counts)),
  paste0("MIN_SIDE_N: ", MIN_SIDE_N),
  paste0("Methods: ", paste(METHODS, collapse=", ")),
  "Selection rule: Pareto maximize PC1 capture / minimize LE-REM balance cost; maximum endpoint-chord deviation",
  "Balance components: abundance SMD, dispersion-residual SMD, dispersion-residual log variance ratio",
  "Downstream significance is not used in cutoff selection."
), file.path(OUT_ROOT, "Normalization_Aware_Balanced_EVS_Audit.txt"))

message("\n============================================================")
message("Normalization-aware balanced EVS cutoff analysis complete")
message("Output: ", normalizePath(OUT_ROOT, winslash="/", mustWork=TRUE))
message("============================================================")
