#!/usr/bin/env Rscript

# =============================================================================
# SEQUENCE: PC1 / NEGATIVE-BINOMIAL RESIDUAL CUTOFF -- NO PERMUTATION
# Actual-data version for WTTS-Seq
# =============================================================================
# PURPOSE
#   Derive a data-driven EVS cutoff k* for each RT/ZT comparison from the real
#   WTTS-Seq raw-count matrix without permutation tests, Storey pi0, downstream
#   DE p-values, HBFSS, or TWAS.
#
# PRINCIPLE
#   1) DESeq2 median-of-ratios normalization and MAP NB dispersions are fitted
#      from the raw count matrix.
#   2) For each RT/ZT comparison, a pooled NB Pearson-residual PCA provides an
#      abundance-adjusted measure of PC1-attributable structure per PAS.
#   3) Actual EVS is performed independently in RT and ZT on DESeq2-normalized
#      counts. For candidate k:
#         LE(k)  = TopK_RT(k) union TopK_ZT(k)
#         REM(k) = all retained PASs \ LE(k)
#   4) Each candidate k is evaluated by:
#         BENEFIT: fraction of pooled NB-residual PC1 mass captured by LE(k)
#         COST:    LE/REM imbalance in mean-adjusted dispersion location and
#                  dispersion-residual spread. Abundance differences are
#                  reported descriptively but are NOT penalized, because LE is
#                  expected to have a higher mean and therefore higher raw
#                  variance under the NB mean-variance relationship.
#   5) The nondominated Pareto frontier is retained (maximize PC1 capture,
#      minimize mean-adjusted dispersion balance cost), and k* is the frontier point with maximum
#      perpendicular deviation from the endpoint chord.
#
# IMPORTANT
#   * k* is the scalar cutoff; EVS determines the actual Leading Edge by the
#     union of the two arm-specific top-k* sets.
#   * No permutation is used anywhere in this script.
#   * No downstream significance result is used to choose k*.
#   * This objective is designed to preserve structured PC1 signal while avoiding
#     excess mean-adjusted dispersion imbalance; it does not guarantee the maximum
#     number of downstream significant sites and must be validated downstream.
#
# FIGURES (two multi-panel figures, 180 mm wide, PDF + 600 dpi PNG)
#   Figure_1_Cutoff_Selection        a PC1 capture vs k | b dispersion-balance cost vs k | c Pareto frontier
#   Figure_2_EVS_Partition_at_kstar  a mean-dispersion landscape | b abundance | c dispersion residual
# =============================================================================

options(stringsAsFactors = FALSE)
set.seed(98632)

required <- c("DESeq2", "ggplot2", "patchwork", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(patchwork)
})

# =============================================================================
# SETTINGS
# =============================================================================

env_num <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) as.numeric(v) else default
}

MIN_COUNT      <- env_num("PC1_NB_MIN_COUNT", 10)
MIN_SIDE_ABS   <- as.integer(env_num("PC1_NB_MIN_SIDE_ABS", 100))
MIN_SIDE_FRAC  <- env_num("PC1_NB_MIN_SIDE_FRAC", 0.02)
FIXED_K        <- 5000L

GROUP_PATTERNS <- c(
  RT0 = "^R0_", ZT6 = "^ZT6_", RT2 = "^R2_", ZT8 = "^ZT8_",
  RT4 = "^R4_", ZT10 = "^ZT10_", RT8 = "^R8_", ZT14 = "^ZT14_"
)

COMPARISONS <- list(
  RT0_ZT6  = c(rt = "RT0", zt = "ZT6"),
  RT2_ZT8  = c(rt = "RT2", zt = "ZT8"),
  RT4_ZT10 = c(rt = "RT4", zt = "ZT10"),
  RT8_ZT14 = c(rt = "RT8", zt = "ZT14")
)

# =============================================================================
# PATHS
# =============================================================================

COUNT_FILE <- "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
root_candidates <- unique(c(Sys.getenv("SEQUENCE_REPO_ROOT", unset = ""),
                            getwd(), dirname(getwd()), "/root/REAPER98632"))
root_candidates <- root_candidates[nzchar(root_candidates)]

repo_root <- NULL
for (cand in root_candidates) {
  if (file.exists(file.path(cand, COUNT_FILE)) ||
      file.exists(file.path(cand, "data", COUNT_FILE))) {
    repo_root <- normalizePath(cand, winslash = "/", mustWork = TRUE)
    break
  }
}
if (is.null(repo_root)) repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

count_path <- c(file.path(repo_root, COUNT_FILE), file.path(repo_root, "data", COUNT_FILE))
count_path <- count_path[file.exists(count_path)][1]
if (is.na(count_path)) stop("Raw-count file not found under ", repo_root)

OUT_ROOT <- Sys.getenv(
  "SEQUENCE_PC1_NB_OUT",
  unset = file.path(repo_root, "exports", "PC1_NB_Residual_Cutoff_NoPermutation")
)
FIG_DIR <- file.path(OUT_ROOT, "Figures")
TAB_DIR <- file.path(OUT_ROOT, "Tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

message("Repository: ", repo_root)
message("Counts:     ", count_path)
message("Output:     ", OUT_ROOT)
message("Permutation: NONE")

# =============================================================================
# INPUT
# =============================================================================

read_counts <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  sample_idx <- sort(unique(unlist(lapply(GROUP_PATTERNS, function(p) grep(p, names(x))))))
  if (!length(sample_idx)) stop("No sample columns matched GROUP_PATTERNS.")

  id_candidates <- c("OrigID", "feature_id", "FeatureID", "PAS", "pas_id",
                     "GeneID", "gene_id", "id")
  id_hit <- id_candidates[id_candidates %in% names(x)]
  feature_col <- if (length(id_hit)) id_hit[1] else names(x)[1]

  feature_id <- trimws(as.character(x[[feature_col]]))
  blank <- is.na(feature_id) | feature_id == ""
  feature_id[blank] <- paste0("feature_", which(blank))

  m <- do.call(cbind, lapply(x[, sample_idx, drop = FALSE], function(z) {
    suppressWarnings(as.numeric(as.character(z)))
  }))
  colnames(m) <- names(x)[sample_idx]

  n_bad <- sum(!is.finite(m))
  if (n_bad > 0) warning(n_bad, " missing/non-numeric count entries were set to 0.")
  n_neg <- sum(m < 0, na.rm = TRUE)
  if (n_neg > 0) warning(n_neg, " negative count entries were set to 0.")
  m[!is.finite(m)] <- 0
  m <- pmax(m, 0)

  if (anyDuplicated(feature_id)) {
    m <- rowsum(m, group = feature_id, reorder = FALSE)
  } else {
    rownames(m) <- feature_id
  }

  m[rowSums(m) > 0, , drop = FALSE]
}

assign_groups <- function(samples) {
  g <- rep(NA_character_, length(samples))
  names(g) <- samples
  for (nm in names(GROUP_PATTERNS)) {
    idx <- grep(GROUP_PATTERNS[[nm]], samples)
    if (any(!is.na(g[idx]))) stop("A sample matched more than one group.")
    g[idx] <- nm
  }
  if (anyNA(g)) stop("Unassigned samples: ", paste(names(g)[is.na(g)], collapse = ", "))
  factor(g, levels = names(GROUP_PATTERNS))
}

counts_all <- read_counts(count_path)
groups <- assign_groups(colnames(counts_all))
group_n <- table(groups)
if (any(group_n < 2)) {
  stop("Every group needs >=2 samples; found: ",
       paste(names(group_n), group_n, sep = "=", collapse = ", "))
}

min_group_n <- min(group_n)
keep_rows <- rowSums(counts_all >= MIN_COUNT) >= min_group_n
count_mat <- counts_all[keep_rows, , drop = FALSE]
message("Nonzero PASs: ", nrow(counts_all), "; after count filter: ", nrow(count_mat))

# =============================================================================
# DESeq2 MEDIAN-OF-RATIOS NORMALIZATION + NB DISPERSION
# =============================================================================

dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = data.frame(group = groups, row.names = colnames(count_mat)),
  design = ~ group
)
dds <- tryCatch(
  estimateSizeFactors(dds),
  error = function(e) estimateSizeFactors(dds, type = "poscounts")
)
dds <- tryCatch(
  estimateDispersions(dds, quiet = TRUE),
  error = function(e) {
    message("Default dispersion fit failed (", conditionMessage(e), "); using fitType='mean'.")
    estimateDispersions(dds, fitType = "mean", quiet = TRUE)
  }
)

y_raw <- counts(dds, normalized = FALSE)
y_norm <- counts(dds, normalized = TRUE)
sf_all <- sizeFactors(dds)
alpha <- dispersions(dds)
mc <- mcols(dds)
base_mean <- mc$baseMean
disp_trend <- mc$dispFit

disp_ok <- is.finite(alpha) & alpha > 0 & is.finite(base_mean) & base_mean > 0 &
           is.finite(disp_trend) & disp_trend > 0
if (sum(disp_ok) < 100) stop("Too few PASs with finite positive DESeq2 dispersion information.")

y_raw <- y_raw[disp_ok, , drop = FALSE]
y_norm <- y_norm[disp_ok, , drop = FALSE]
alpha <- alpha[disp_ok]
base_mean <- base_mean[disp_ok]
disp_trend <- disp_trend[disp_ok]
feature_ids <- rownames(y_raw)
disp_resid <- log(alpha / disp_trend)
abundance_axis <- log1p(base_mean)

message("PASs with usable dispersion information: ", length(feature_ids))

# =============================================================================
# HELPERS
# =============================================================================

minmax01 <- function(x) {
  out <- rep(0, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  lo <- min(x[ok]); hi <- max(x[ok])
  if (!is.finite(hi - lo) || hi <= lo) return(out)
  out[ok] <- (x[ok] - lo) / (hi - lo)
  out
}

pooled_residual_pc1 <- function(y, sf, alpha) {
  mu <- rowMeans(sweep(y, 2, sf, "/"))
  m <- outer(mu, sf)
  z <- (y - m) / sqrt(m + alpha * m^2)
  zc <- z - rowMeans(z)
  ss <- rowSums(zc^2)
  usable <- mu > 0 & is.finite(ss) & ss > 1e-12

  zcu <- zc[usable, , drop = FALSE]
  e <- eigen(crossprod(zcu), symmetric = TRUE)
  ev <- pmax(e$values, 0)
  u1 <- e$vectors[, 1]
  proj <- as.vector(zcu %*% u1)

  pc1_mass <- rep(0, nrow(y))
  cos2 <- rep(NA_real_, nrow(y))
  resid_var <- rep(NA_real_, nrow(y))
  pc1_mass[usable] <- proj^2 / (ncol(zcu) - 1)
  cos2[usable] <- proj^2 / rowSums(zcu^2)
  resid_var[usable] <- rowSums(zcu^2) / (ncol(zcu) - 1)

  list(
    pc1_mass = pc1_mass,
    cos2 = cos2,
    resid_var = resid_var,
    mu = mu,
    usable = usable,
    pc1_var_frac = if (sum(ev) > 0) ev[1] / sum(ev) else NA_real_,
    scores = setNames(u1 * sqrt(ev[1]), colnames(zcu))
  )
}

# Actual EVS rank: highest |PC1 loading| gets rank 1.
evs_pc1_rank <- function(x) {
  x <- as.matrix(x)
  rv <- apply(x, 1, var)
  usable <- is.finite(rv) & rv > 1e-12

  rk <- rep(NA_integer_, nrow(x))
  loading <- rep(NA_real_, nrow(x))
  pc1_mass <- rep(0, nrow(x))

  if (sum(usable) < 2L) return(list(rank = rk, loading = loading, pc1_mass = pc1_mass))

  pca <- prcomp(t(x[usable, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 1)
  ld <- as.numeric(pca$rotation[, 1])
  mass <- (pca$sdev[1]^2) * ld^2
  ord <- order(abs(ld), mass, decreasing = TRUE)
  idx <- which(usable)

  rk[idx[ord]] <- seq_along(ord)
  loading[idx] <- ld
  pc1_mass[idx] <- mass
  list(rank = rk, loading = loading, pc1_mass = pc1_mass)
}

# Efficient candidate-k scan. Union membership is entry_rank <= k.
scan_candidates <- function(cn, rt_rank, zt_rank, residual_pc1_mass,
                            abundance_axis, disp_resid) {
  entry_rank <- pmin(rt_rank, zt_rank, na.rm = TRUE)
  entry_rank[!is.finite(entry_rank)] <- NA_real_

  valid_metric <- is.finite(abundance_axis) & is.finite(disp_resid) &
                  is.finite(residual_pc1_mass) & residual_pc1_mass >= 0
  valid_entry <- is.finite(entry_rank)
  valid <- valid_metric & valid_entry
  if (sum(valid) < 200L) stop(cn, ": too few PASs for candidate scan.")

  er <- entry_rank[valid]
  ab <- abundance_axis[valid]
  dr <- disp_resid[valid]
  pm <- residual_pc1_mass[valid]

  o <- order(er, seq_along(er))
  er <- er[o]; ab <- ab[o]; dr <- dr[o]; pm <- pm[o]
  n <- length(er)
  min_side_n <- max(MIN_SIDE_ABS, ceiling(MIN_SIDE_FRAC * n))

  # Candidate k only where membership actually changes.
  end_idx <- c(which(diff(er) != 0), n)
  k_vals <- er[end_idx]
  le_n <- end_idx
  rem_n <- n - le_n
  keep <- le_n >= min_side_n & rem_n >= min_side_n & is.finite(k_vals)
  end_idx <- end_idx[keep]
  k_vals <- as.integer(k_vals[keep])
  le_n <- le_n[keep]
  rem_n <- rem_n[keep]
  if (length(k_vals) < 3L) stop(cn, ": fewer than 3 valid candidate cutoffs after side-size constraint.")

  c_ab <- cumsum(ab); c_ab2 <- cumsum(ab^2)
  c_dr <- cumsum(dr); c_dr2 <- cumsum(dr^2)
  c_pm <- cumsum(pm)
  T_ab <- sum(ab); T_ab2 <- sum(ab^2)
  T_dr <- sum(dr); T_dr2 <- sum(dr^2)
  T_pm <- sum(pm)
  if (!is.finite(T_pm) || T_pm <= 0) stop(cn, ": residual PC1 mass is zero/undefined.")

  mean_var_from_sums <- function(sumx, sumx2, nn) {
    mn <- sumx / nn
    vv <- ifelse(nn > 1, pmax((sumx2 - (sumx^2)/nn) / (nn - 1), 0), NA_real_)
    list(mean = mn, var = vv)
  }

  le_ab <- mean_var_from_sums(c_ab[end_idx], c_ab2[end_idx], le_n)
  rem_ab <- mean_var_from_sums(T_ab - c_ab[end_idx], T_ab2 - c_ab2[end_idx], rem_n)
  le_dr <- mean_var_from_sums(c_dr[end_idx], c_dr2[end_idx], le_n)
  rem_dr <- mean_var_from_sums(T_dr - c_dr[end_idx], T_dr2 - c_dr2[end_idx], rem_n)

  pooled_ab <- sqrt(pmax(((le_n - 1) * le_ab$var + (rem_n - 1) * rem_ab$var) / (n - 2), 1e-12))
  pooled_dr <- sqrt(pmax(((le_n - 1) * le_dr$var + (rem_n - 1) * rem_dr$var) / (n - 2), 1e-12))

  abundance_smd <- abs(le_ab$mean - rem_ab$mean) / pooled_ab
  dispersion_smd <- abs(le_dr$mean - rem_dr$mean) / pooled_dr
  dispersion_spread <- abs(log((le_dr$var + 1e-12) / (rem_dr$var + 1e-12)))
  pc1_capture <- c_pm[end_idx] / T_pm

  scan <- data.frame(
    comparison = cn,
    k = k_vals,
    LE_n = le_n,
    REM_n = rem_n,
    PC1_capture = pc1_capture,
    abundance_SMD = abundance_smd,
    dispersion_residual_SMD = dispersion_smd,
    dispersion_residual_log_variance_ratio = dispersion_spread,
    abundance_mean_LE = le_ab$mean,
    abundance_mean_REM = rem_ab$mean,
    dispersion_residual_mean_LE = le_dr$mean,
    dispersion_residual_mean_REM = rem_dr$mean,
    dispersion_residual_var_LE = le_dr$var,
    dispersion_residual_var_REM = rem_dr$var,
    stringsAsFactors = FALSE
  )

  scan$abundance_cost01 <- minmax01(scan$abundance_SMD)
  scan$dispersion_location_cost01 <- minmax01(scan$dispersion_residual_SMD)
  scan$dispersion_spread_cost01 <- minmax01(scan$dispersion_residual_log_variance_ratio)
  # FINAL OBJECTIVE: abundance is descriptive only. A higher LE mean is expected
  # and should naturally imply higher raw variance under Var(Y)=mu+alpha*mu^2.
  # Penalize only mean-adjusted dispersion location and spread differences.
  scan$balance_cost <- rowMeans(scan[, c("dispersion_location_cost01",
                                         "dispersion_spread_cost01")])

  # Pareto frontier: minimize balance_cost, maximize PC1_capture.
  tmp <- scan[order(scan$balance_cost, -scan$PC1_capture, -scan$k), , drop = FALSE]
  tmp <- tmp[!duplicated(signif(tmp$balance_cost, 12)), , drop = FALSE]
  prior_best <- c(-Inf, head(cummax(tmp$PC1_capture), -1L))
  frontier <- tmp[tmp$PC1_capture > prior_best, , drop = FALSE]
  frontier <- frontier[order(frontier$balance_cost, frontier$PC1_capture, frontier$k), , drop = FALSE]
  if (nrow(frontier) < 2L) stop(cn, ": Pareto frontier has fewer than 2 points.")

  x <- minmax01(frontier$balance_cost)
  y <- minmax01(frontier$PC1_capture)
  x1 <- x[1]; y1 <- y[1]
  x2 <- x[length(x)]; y2 <- y[length(y)]
  denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
  if (!is.finite(denom) || denom <= .Machine$double.eps) {
    # Rare degenerate frontier: choose lowest cost, then highest capture.
    frontier$endpoint_deviation <- 0
    ord_best <- order(frontier$balance_cost, -frontier$PC1_capture, frontier$k)
    selected <- frontier[ord_best[1], , drop = FALSE]
  } else {
    frontier$endpoint_deviation <- abs((y2 - y1) * x - (x2 - x1) * y + x2*y1 - y2*x1) / denom
    best_dev <- max(frontier$endpoint_deviation, na.rm = TRUE)
    cand <- frontier[abs(frontier$endpoint_deviation - best_dev) < 1e-12, , drop = FALSE]
    cand <- cand[order(-cand$PC1_capture, cand$balance_cost, cand$k), , drop = FALSE]
    selected <- cand[1, , drop = FALSE]
  }

  list(scan = scan, frontier = frontier, selected = selected, entry_rank = entry_rank)
}

save_fig <- function(p, stem, w_mm = 180, h_mm = 95) {
  pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p,
         width = w_mm, height = h_mm, units = "mm", device = pdf_dev)
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p,
         width = w_mm, height = h_mm, units = "mm", dpi = 600, bg = "white")
}

# =============================================================================
# RUN EACH COMPARISON
# =============================================================================

all_scan <- list()
all_frontier <- list()
all_selected <- list()
all_membership <- list()
all_feature_scores <- list()

for (cn in names(COMPARISONS)) {
  mp <- COMPARISONS[[cn]]
  pair_idx <- which(groups %in% mp)
  rt_idx <- which(groups == mp[["rt"]])
  zt_idx <- which(groups == mp[["zt"]])

  message("\n", cn, ": ", length(pair_idx), " samples (",
          mp[["rt"]], " n=", length(rt_idx), "; ", mp[["zt"]], " n=", length(zt_idx), ")")

  # Abundance-adjusted PC1 structure from pooled NB Pearson residuals.
  r <- pooled_residual_pc1(y_raw[, pair_idx, drop = FALSE], sf_all[pair_idx], alpha)

  # Actual EVS ranks are arm-specific and based on DESeq2 median-of-ratios normalized counts.
  er_rt <- evs_pc1_rank(y_norm[, rt_idx, drop = FALSE])
  er_zt <- evs_pc1_rank(y_norm[, zt_idx, drop = FALSE])

  opt <- scan_candidates(
    cn = cn,
    rt_rank = er_rt$rank,
    zt_rank = er_zt$rank,
    residual_pc1_mass = r$pc1_mass,
    abundance_axis = abundance_axis,
    disp_resid = disp_resid
  )

  kstar <- as.integer(opt$selected$k[1])
  message("  k* = ", kstar,
          " | LE = ", opt$selected$LE_n[1],
          " | REM = ", opt$selected$REM_n[1],
          " | residual PC1 capture = ", sprintf("%.4f", opt$selected$PC1_capture[1]),
          " | balance cost = ", sprintf("%.4f", opt$selected$balance_cost[1]))

  # Raw-EVS ranking retained as a diagnostic comparator only.
  raw_rt <- evs_pc1_rank(y_raw[, rt_idx, drop = FALSE])
  raw_zt <- evs_pc1_rank(y_raw[, zt_idx, drop = FALSE])

  norm_le <- is.finite(opt$entry_rank) & opt$entry_rank <= kstar
  raw_entry <- pmin(raw_rt$rank, raw_zt$rank, na.rm = TRUE)
  raw_entry[!is.finite(raw_entry)] <- NA_real_
  raw_le <- is.finite(raw_entry) & raw_entry <= kstar

  all_scan[[cn]] <- opt$scan
  all_frontier[[cn]] <- opt$frontier
  sel <- opt$selected
  sel$PC1_var_frac_pooled_residual <- r$pc1_var_frac
  sel$max_sample_share_PC1 <- max(r$scores^2) / sum(r$scores^2)
  sel$NormEVS_LE_union_all_features <- sum(norm_le)
  sel$NormEVS_REM_all_features <- length(feature_ids) - sum(norm_le)
  sel$RawEVS_LE_union_all_features <- sum(raw_le)
  sel$RawEVS_REM_all_features <- length(feature_ids) - sum(raw_le)
  all_selected[[cn]] <- sel

  all_membership[[cn]] <- data.frame(
    feature_id = feature_ids,
    comparison = cn,
    control_rank = er_rt$rank,
    treatment_rank = er_zt$rank,
    entry_rank = opt$entry_rank,
    control_loading = er_rt$loading,
    treatment_loading = er_zt$loading,
    pooled_residual_PC1_mass = r$pc1_mass,
    pooled_residual_cos2 = r$cos2,
    pooled_residual_variance = r$resid_var,
    baseMean = base_mean,
    dispersion = alpha,
    dispersion_trend = disp_trend,
    dispersion_residual = disp_resid,
    selected_k = kstar,
    dataset = ifelse(norm_le, "Leading_Edge", "Remainder"),
    raw_entry_rank = raw_entry,
    raw_dataset = ifelse(raw_le, "Leading_Edge", "Remainder"),
    stringsAsFactors = FALSE
  )

  all_feature_scores[[cn]] <- data.frame(
    feature_id = feature_ids,
    comparison = cn,
    pooled_residual_PC1_mass = r$pc1_mass,
    pooled_residual_cos2 = r$cos2,
    pooled_residual_variance = r$resid_var,
    normalized_EVS_control_rank = er_rt$rank,
    normalized_EVS_treatment_rank = er_zt$rank,
    normalized_EVS_entry_rank = opt$entry_rank,
    stringsAsFactors = FALSE
  )
}

scan_df <- do.call(rbind, all_scan)
frontier_df <- do.call(rbind, all_frontier)
selected_df <- do.call(rbind, all_selected)
membership_df <- do.call(rbind, all_membership)
feature_scores_df <- do.call(rbind, all_feature_scores)
rownames(scan_df) <- rownames(frontier_df) <- rownames(selected_df) <-
  rownames(membership_df) <- rownames(feature_scores_df) <- NULL

# =============================================================================
# OUTPUT TABLES
# =============================================================================

summary_tab <- selected_df[, c(
  "comparison", "k", "LE_n", "REM_n", "PC1_capture", "balance_cost",
  "abundance_SMD", "dispersion_residual_SMD",
  "dispersion_residual_log_variance_ratio",
  "PC1_var_frac_pooled_residual", "max_sample_share_PC1",
  "NormEVS_LE_union_all_features", "NormEVS_REM_all_features",
  "RawEVS_LE_union_all_features", "RawEVS_REM_all_features",
  "endpoint_deviation"
)]
names(summary_tab)[names(summary_tab) == "k"] <- "k_star"

manifest <- data.frame(
  Cutoff_Method = "PC1_NB_Residual_NoPermutation",
  RT0_ZT6 = summary_tab$k_star[summary_tab$comparison == "RT0_ZT6"],
  RT2_ZT8 = summary_tab$k_star[summary_tab$comparison == "RT2_ZT8"],
  RT4_ZT10 = summary_tab$k_star[summary_tab$comparison == "RT4_ZT10"],
  RT8_ZT14 = summary_tab$k_star[summary_tab$comparison == "RT8_ZT14"],
  Source = paste(
    "No-permutation cutoff: actual DESeq2-normalized EVS LE/REM split; Pareto knee",
    "maximizes pooled NB Pearson-residual PC1 mass capture while minimizing LE/REM",
    "mean-adjusted dispersion location and spread imbalance; abundance reported descriptively only"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(summary_tab, file.path(TAB_DIR, "Table_1_Selected_Cutoffs.csv"), row.names = FALSE)
write.csv(scan_df, file.path(TAB_DIR, "Table_2_All_Candidate_k.csv"), row.names = FALSE)
write.csv(frontier_df, file.path(TAB_DIR, "Table_3_Pareto_Frontiers.csv"), row.names = FALSE)
write.csv(membership_df, file.path(TAB_DIR, "Table_4_PAS_Membership.csv"), row.names = FALSE)
write.csv(feature_scores_df, file.path(TAB_DIR, "Table_S1_PC1_Residual_Feature_Scores.csv"), row.names = FALSE)
write.csv(manifest, file.path(OUT_ROOT, "PC1_NB_Residual_NoPermutation_Manifest.csv"), row.names = FALSE)

# =============================================================================
# FIGURES -- two multi-panel manuscript figures (180 mm, 8 pt sans)
# =============================================================================
#   Figure 1  Cutoff selection
#             a  Residual PC1 mass captured by the Leading Edge vs k
#             b  dispersion-balance cost vs k (abundance shown as diagnostic only)
#             c  Pareto frontier with the selected k*
#   Figure 2  Actual EVS partition at k*
#             a  Mean-dispersion landscape (Leading Edge over Remainder)
#             b  Abundance distribution, LE vs REM
#             c  Dispersion-residual distribution, LE vs REM
#   Figure data use exactly the PASs evaluated in the candidate scan, so the
#   plotted distributions match the SMDs reported in Table 1.

COL_LE  <- "#0072B2"
COL_REM <- "#9A9A9A"
COL_SET <- c("Leading Edge" = COL_LE, "Remainder" = COL_REM)
COL_COST <- c("Dispersion balance cost"       = "black",
              "Abundance difference (diagnostic only)" = "#E69F00",
              "Dispersion location imbalance" = "#009E73",
              "Dispersion spread imbalance"   = "#CC79A7")

comp_levels <- names(COMPARISONS)
comp_labels <- setNames(sub("_", " / ", comp_levels), comp_levels)
facet_lab   <- labeller(comparison = comp_labels)

fmt     <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
lab_log <- scales::trans_format("log10", scales::math_format(10^.x))
lab_pct <- scales::label_percent(accuracy = 0.01)
k_axis  <- scale_x_continuous(labels = scales::label_comma(), breaks = scales::breaks_extended(n = 3))

theme_ms <- theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.3),
    axis.text        = element_text(colour = "black", size = 6.5),
    axis.title       = element_text(size = 7.5),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 7.5),
    panel.spacing.x  = unit(4, "mm"),
    legend.position  = "none",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 7),
    legend.key.width = unit(6, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.margin    = margin(0, 0, 0, 0),
    plot.title       = element_text(size = 8, face = "bold", margin = margin(b = 3)),
    plot.tag         = element_text(size = 10, face = "bold"),
    plot.margin      = margin(3, 4, 3, 3)
  )

# Rasterise dense point layers when ggrastr is available (keeps PDFs small).
point_layer <- function(...) {
  if (requireNamespace("ggrastr", quietly = TRUE)) ggrastr::geom_point_rast(..., raster.dpi = 600)
  else geom_point(...)
}

as_comp <- function(d) { d$comparison <- factor(d$comparison, levels = comp_levels); d }
sel_df      <- as_comp(summary_tab)
scan_plot   <- as_comp(scan_df)
front_plot  <- as_comp(frontier_df)

# ---- Figure 1: cutoff selection ---------------------------------------------
kstar_lines <- list(
  geom_vline(data = sel_df, aes(xintercept = k_star), colour = COL_LE,
             linetype = "dashed", linewidth = 0.35),
  geom_vline(xintercept = FIXED_K, colour = "grey45", linetype = "dotted", linewidth = 0.3)
)

f1a <- ggplot(scan_plot, aes(k, PC1_capture)) +
  kstar_lines +
  geom_line(linewidth = 0.45) +
  geom_text(data = sel_df, aes(x = k_star, y = -Inf, label = paste0("k* = ", fmt(k_star))),
            inherit.aes = FALSE, hjust = 1.08, vjust = -0.7, size = 2.2, colour = COL_LE) +
  facet_wrap(~ comparison, nrow = 1, scales = "free", labeller = facet_lab) +
  k_axis +
  scale_y_continuous(labels = lab_pct) +
  labs(x = NULL, y = "Residual PC1 mass\ncaptured by LE",
       title = "PC1 capture across candidate cutoffs (dashed: k*; dotted: k = 5,000)") +
  theme_ms

cost_long <- rbind(
  data.frame(comparison = scan_plot$comparison, k = scan_plot$k,
             component = "Abundance difference (diagnostic only)", value = scan_plot$abundance_cost01),
  data.frame(comparison = scan_plot$comparison, k = scan_plot$k,
             component = "Dispersion location imbalance", value = scan_plot$dispersion_location_cost01),
  data.frame(comparison = scan_plot$comparison, k = scan_plot$k,
             component = "Dispersion spread imbalance", value = scan_plot$dispersion_spread_cost01),
  data.frame(comparison = scan_plot$comparison, k = scan_plot$k,
             component = "Dispersion balance cost", value = scan_plot$balance_cost)
)
cost_long$component <- factor(cost_long$component, levels = names(COL_COST))
is_comb <- cost_long$component == "Dispersion balance cost"

f1b <- ggplot(cost_long, aes(k, value, colour = component)) +
  kstar_lines +
  geom_line(data = cost_long[!is_comb, ], linewidth = 0.3, alpha = 0.85) +
  geom_line(data = cost_long[is_comb, ], linewidth = 0.6) +
  facet_wrap(~ comparison, nrow = 1, scales = "free_x", labeller = facet_lab) +
  k_axis +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_colour_manual(values = COL_COST, breaks = names(COL_COST)) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "EVS cutoff k (per condition)", y = "Scaled imbalance\n(0 = best k)",
       title = "LE/REM balance cost across candidate cutoffs") +
  theme_ms + theme(legend.position = "bottom")

f1c <- ggplot(front_plot, aes(balance_cost, PC1_capture)) +
  geom_path(linewidth = 0.4, colour = "grey35") +
  geom_point(size = 0.8, colour = "grey35") +
  geom_point(data = sel_df, aes(balance_cost, PC1_capture), inherit.aes = FALSE,
             shape = 21, size = 2.6, stroke = 0.5, fill = COL_LE, colour = "white") +
  geom_text(data = sel_df, aes(balance_cost, PC1_capture, label = paste0("k* = ", fmt(k_star))),
            inherit.aes = FALSE, vjust = -1.2, size = 2.2, colour = COL_LE) +
  facet_wrap(~ comparison, nrow = 1, scales = "free", labeller = facet_lab) +
  scale_x_continuous(breaks = scales::breaks_extended(n = 3)) +
  scale_y_continuous(labels = lab_pct, expand = expansion(mult = c(0.08, 0.2))) +
  labs(x = "LE/REM balance cost (lower is better)", y = "Residual PC1 mass\ncaptured by LE",
       title = "Pareto frontier; k* = maximum deviation from the endpoint chord") +
  theme_ms

fig1 <- (f1a / f1b / f1c) + plot_annotation(tag_levels = "a")
save_fig(fig1, "Figure_1_Cutoff_Selection", 180, 170)

# ---- Figure 2: actual EVS partition at k* -----------------------------------
part_df <- membership_df[is.finite(membership_df$entry_rank) &
                         is.finite(membership_df$dispersion_residual) &
                         is.finite(membership_df$pooled_residual_PC1_mass) &
                         membership_df$pooled_residual_PC1_mass >= 0, ]
part_df <- as_comp(part_df)
part_df$set <- factor(ifelse(part_df$dataset == "Leading_Edge", "Leading Edge", "Remainder"),
                      levels = c("Remainder", "Leading Edge"))
part_df <- part_df[order(part_df$set), ]   # Remainder drawn first, Leading Edge on top

trend_df <- data.frame(baseMean = base_mean, dispersion_trend = disp_trend)
trend_df <- trend_df[order(trend_df$baseMean), ]

n_lab <- data.frame(comparison = sel_df$comparison,
                    label = sprintf("LE %s\nREM %s", fmt(sel_df$LE_n), fmt(sel_df$REM_n)))

f2a <- ggplot(part_df, aes(baseMean, dispersion, colour = set)) +
  point_layer(size = 0.2, alpha = 0.3, stroke = 0) +
  geom_line(data = trend_df, aes(baseMean, dispersion_trend), inherit.aes = FALSE,
            colour = "black", linewidth = 0.45) +
  geom_text(data = n_lab, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.2, size = 2.2, lineheight = 0.95) +
  facet_wrap(~ comparison, nrow = 1, labeller = facet_lab) +
  scale_x_log10(labels = lab_log) +
  scale_y_log10(labels = lab_log) +
  scale_colour_manual(values = COL_SET) +
  labs(x = "DESeq2 normalized mean", y = "DESeq2 MAP dispersion",
       title = "Mean-dispersion landscape at k* (blue: Leading Edge; grey: Remainder; line: fitted trend)") +
  theme_ms

smd_lab <- function(txt) data.frame(comparison = sel_df$comparison, label = txt)
dens_layers <- list(
  geom_density(aes(colour = set, linetype = set), linewidth = 0.55, adjust = 1.1),
  facet_wrap(~ comparison, nrow = 1, labeller = facet_lab),
  scale_colour_manual(values = COL_SET, breaks = c("Leading Edge", "Remainder")),
  scale_linetype_manual(values = c("Leading Edge" = "solid", "Remainder" = "22"),
                        breaks = c("Leading Edge", "Remainder")),
  theme_ms
)

f2b <- ggplot(part_df, aes(log1p(baseMean))) +
  dens_layers +
  geom_text(data = smd_lab(sprintf("|SMD| = %.2f", sel_df$abundance_SMD)),
            aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.3, size = 2.2) +
  labs(x = "log(1 + DESeq2 baseMean)", y = "Density",
       title = "Abundance: Leading Edge vs Remainder") +
  theme(strip.text = element_blank())

xl <- quantile(part_df$dispersion_residual, c(0.005, 0.995), na.rm = TRUE)
f2c <- ggplot(part_df, aes(dispersion_residual)) +
  dens_layers +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_text(data = smd_lab(sprintf("|SMD| = %.2f\n|log VR| = %.2f", sel_df$dispersion_residual_SMD,
                                   sel_df$dispersion_residual_log_variance_ratio)),
            aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.2, size = 2.2, lineheight = 0.95) +
  coord_cartesian(xlim = xl) +
  labs(x = "Dispersion residual log(alpha / fitted trend); central 99% shown", y = "Density",
       title = "Mean-adjusted dispersion: Leading Edge vs Remainder") +
  theme(strip.text = element_blank(), legend.position = "bottom")

fig2 <- (f2a / f2b / f2c) + plot_layout(heights = c(1.25, 1, 1)) +
  plot_annotation(tag_levels = "a")
save_fig(fig2, "Figure_2_EVS_Partition_at_kstar", 180, 170)

# =============================================================================
# METHODS + AUDIT
# =============================================================================

fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
sel_txt <- paste(sprintf(
  "%s: k*=%s, LE=%s, REM=%s, PC1 capture=%.3f, balance cost=%.3f",
  summary_tab$comparison, fmt(summary_tab$k_star),
  fmt(summary_tab$NormEVS_LE_union_all_features),
  fmt(summary_tab$NormEVS_REM_all_features),
  summary_tab$PC1_capture, summary_tab$balance_cost
), collapse = "; ")

methods_text <- c(
  "Data-derived EVS cutoff without permutation",
  "",
  paste0(
    "PASs with at least ", MIN_COUNT, " reads in at least ", min_group_n,
    " samples were retained. DESeq2 median-of-ratios size factors and MAP negative-binomial dispersions were estimated from the full raw-count matrix with design ~ group. ",
    "For each RT/ZT comparison, pooled negative-binomial Pearson residuals were calculated as z_ij=(y_ij-s_j*mu_i)/sqrt(s_j*mu_i+alpha_i*(s_j*mu_i)^2), using the pooled comparison mean mu_i. ",
    "PCA of this residual matrix yielded a PAS-specific PC1-attributable residual variance P_i=(z_i'u_1)^2/(n-1), thereby measuring coordinated PC1 structure after accounting for expected mean-dependent count variance. ",
    "Actual eigenvector splitting was then performed independently within the RT and ZT arms on DESeq2 median-of-ratios normalized counts, ranking PASs by decreasing absolute PC1 loading. For each candidate k, the Leading Edge was the union of the two arm-specific top-k sets and the Remainder was its complement. ",
    "Candidate cutoffs were evaluated by the fraction of pooled NB-residual PC1 mass captured by the actual Leading Edge and by LE-versus-Remainder imbalance in the mean-adjusted dispersion residual log(alpha_MAP/alpha_trend) and the variance of that dispersion residual. These two dispersion imbalance measures were min-max normalized over candidate k and averaged to form a unitless balance cost. LE-versus-Remainder abundance difference in log(1+baseMean) was retained as a descriptive diagnostic only and was not penalized, because a higher Leading-Edge mean is expected to produce higher raw variance under the negative-binomial mean-variance relationship. ",
    "The final k* was chosen on the nondominated Pareto frontier maximizing residual PC1 capture while minimizing balance cost, using maximum perpendicular deviation from the normalized endpoint chord. No permutation, Storey pi0 estimate, downstream p-value, HBFSS score, or TWAS result was used to select the cutoff."
  ),
  "",
  paste0("Selected solutions: ", sel_txt, "."),
  "",
  paste0("Minimum side-size constraint: max(", MIN_SIDE_ABS, ", ", MIN_SIDE_FRAC,
         " x evaluable PASs). This is a stability constraint only and does not select k*."),
  "",
  "Sample counts:",
  paste(names(group_n), group_n, sep = "=", collapse = ", ")
)
writeLines(methods_text, file.path(OUT_ROOT, "Methods_and_Audit.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT_ROOT, "sessionInfo.txt"))

message("\nSelected cutoffs:")
print(summary_tab[, c("comparison", "k_star", "NormEVS_LE_union_all_features",
                      "NormEVS_REM_all_features", "PC1_capture", "balance_cost")],
      row.names = FALSE, digits = 4)
message("\nFigures: ", FIG_DIR)
message("Tables:  ", TAB_DIR)
message("Methods: ", file.path(OUT_ROOT, "Methods_and_Audit.txt"))
