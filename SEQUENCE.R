
log_step <- function(...) {
  message(sprintf('[%s] %s', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), paste0(..., collapse = '')))
  flush.console()
}

# =============================================================================
# SEQUENCE 10 ACTIVE CLEAN REWRITE
# Current method only:
#   - EVS ranking
#   - local Fourier treatment/control/combined reference curves
#   - last stable crossing before divergence as reference cutoff
#   - raw-Wald bracket search between treatment/control backup ranks
#   - final selected cutoff from minimum total distortion
#   - final split applied to leading-edge and remainder
#   - pure Strimmer fdrtool + hc.thresh
#   - one manuscript volcano per comparison
#   - cleaned tables and figures
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(fdrtool)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(gridExtra)
  library(grid)
  library(scales)
  library(grDevices)
  library(S4Vectors)
})

# -----------------------------------------------------------------------------
# settings
# -----------------------------------------------------------------------------

repo_dir <- getwd()
in_dir <- file.path(repo_dir, "data")
out_root <- file.path(repo_dir, "exports")
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(out_root, paste0("sequence_run_", stamp))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

count_file <- file.path(in_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")

alpha_level <- 0.20
lfc_boundary <- 1.0
hc_upper <- 0.99
lfc_shrink_type <- "apeglm"

fourier_step <- 0.01
fourier_window_frac <- 0.12
fourier_harmonics <- 2L
crossing_min_pct <- 0.05
crossing_max_pct <- 0.95
crossing_stability_n <- 3L
crossing_post_sep_n <- 4L
crossing_post_sep_frac <- 0.75

raw_wald_scan <- TRUE
raw_wald_coarse_step <- 25L
raw_wald_fine_radius <- 50L
raw_wald_min_n <- 200L
raw_wald_w_mean <- 1.00
raw_wald_w_sd <- 1.00
raw_wald_w_skew <- 0.50

fixed_top_n <- 5000L
use_fixed_top_n <- FALSE

figure_dpi <- 320
base_theme_size <- 10
hist_bins <- 60

pal <- list(
  threshold = "#8C2D04",
  hbfss = "#E67E22",
  strong = "#C0392B",
  overlap = "#7D3C98",
  weak = "#4A90E2",
  grey = "#969696",
  control = "#4D4D4D",
  treatment = "#1F78B4"
)

# -----------------------------------------------------------------------------
# metadata
# -----------------------------------------------------------------------------

meta_all <- data.frame(
  id = c(
    "R0_1","R0_2","R0_3","R0_4","R0_5","ZT6_1","ZT6_2","ZT6_3","ZT6_4","ZT6_5",
    "R2_1","R2_2","R2_3","R2_4","R2_5","ZT8_1","ZT8_2","ZT8_3","ZT8_4","ZT8_5",
    "R4_1","R4_2","R4_3","R4_4","R4_5","ZT10_1","ZT10_2","ZT10_3","ZT10_4","ZT10_5",
    "R8_1","R8_2","R8_3","R8_4","R8_5","ZT14_1","ZT14_2","ZT14_3","ZT14_4","ZT14_5"
  ),
  condition = c(
    rep("treatment",5), rep("control",5),
    rep("treatment",5), rep("control",5),
    rep("treatment",5), rep("control",5),
    rep("treatment",5), rep("control",5)
  ),
  stringsAsFactors = FALSE
)
rownames(meta_all) <- meta_all$id
meta_all$condition <- factor(meta_all$condition, levels = c("control", "treatment"))
levels(meta_all$condition) <- c("untrt", "trt")

cmp_tbl <- data.frame(
  cmp = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  trt = c("R0", "R2", "R4", "R8"),
  ctrl = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# helpers
# -----------------------------------------------------------------------------

assert_cols <- function(df, cols, obj = "data frame") {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop(sprintf("%s missing columns: %s", obj, paste(miss, collapse = ", ")), call. = FALSE)
}

save_csv <- function(x, path) write.csv(x, path, row.names = FALSE)

save_plot <- function(p, path, width = 12, height = 8, dpi = figure_dpi) {
  tryCatch(
    ggsave(path, p, width = width, height = height, dpi = dpi, units = "in", limitsize = FALSE, bg = "white"),
    error = function(e) warning(sprintf("plot export failed for %s: %s", basename(path), conditionMessage(e)))
  )
}

safe_log10 <- function(x, eps = 1e-300) log10(pmax(x, eps))
safe_nlog10 <- function(x, eps = 1e-300) -log10(pmax(x, eps))

clip_p <- function(x, eps = 1e-300) {
  x <- as.numeric(x)
  x[is.infinite(x) & x > 0] <- 1 - 1e-12
  x[is.infinite(x) & x < 0] <- eps
  x[is.finite(x)] <- pmin(pmax(x[is.finite(x)], eps), 1 - 1e-12)
  x
}

fmt <- function(x, digits = 3) {
  if (!is.finite(x) || is.na(x)) return("NA")
  formatC(x, digits = digits, format = "fg", flag = "#")
}

theme_seq <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = base_theme_size - 1),
      plot.caption = element_text(size = base_theme_size - 3, colour = "grey30"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88")
    )
}

pretty_ds <- function(x) {
  switch(x,
         raw = "Original dataset",
         lead = "Leading-edge dataset",
         rem = "Remainder dataset",
         x)
}

safe_skew <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) < 5) return(NA_real_)
  sx <- stats::sd(x)
  if (!is.finite(sx) || is.na(sx) || sx == 0) return(NA_real_)
  mean(((x - mean(x)) / sx)^3)
}


roll_mean_dense <- function(x, k = 151L) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0L) return(numeric(0))
  k <- max(3L, as.integer(k))
  if (k %% 2L == 0L) k <- k + 1L
  h <- (k - 1L) / 2L
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    lo <- max(1L, i - h)
    hi <- min(n, i + h)
    vals <- x[lo:hi]
    vals <- vals[is.finite(vals) & !is.na(vals)]
    if (length(vals)) out[i] <- mean(vals)
  }
  out
}

build_nb_rank_profiles <- function(full_results, shared_tbl, prof_k = 151L) {
  res <- as.data.frame(full_results, stringsAsFactors = FALSE)
  assert_cols(res, c("feature_id","baseMean","dispGeneEst"), "full_results")
  assert_cols(shared_tbl, c("feature_id","combined_rank"), "shared_tbl")
  out <- left_join(
    shared_tbl[, c("feature_id","combined_rank"), drop = FALSE],
    res[, c("feature_id","baseMean","dispGeneEst"), drop = FALSE],
    by = "feature_id"
  )
  out <- out[order(out$combined_rank), , drop = FALSE]
  out$mu <- as.numeric(out$baseMean)
  out$alpha <- as.numeric(out$dispGeneEst)
  out <- out[
    is.finite(out$mu) & !is.na(out$mu) & out$mu > 0 &
      is.finite(out$alpha) & !is.na(out$alpha) & out$alpha >= 0,
    ,
    drop = FALSE
  ]
  out$rank <- out$combined_rank
  out$iod <- 1 + out$alpha * out$mu
  out$cv2 <- (1 / out$mu) + out$alpha
  out$iod_prof <- roll_mean_dense(out$iod, prof_k)
  out$cv2_prof <- roll_mean_dense(out$cv2, prof_k)
  out$delta_prof <- out$iod_prof - out$cv2_prof
  out
}

plot_cutoff_explanation_panel <- function(cmp_name, prof_df, cutoff_rank) {
  df <- as.data.frame(prof_df, stringsAsFactors = FALSE)
  req <- c("rank","iod_prof","cv2_prof","delta_prof")
  assert_cols(df, req, "prof_df")
  df <- df[
    is.finite(df$rank) & !is.na(df$rank) &
      is.finite(df$iod_prof) & !is.na(df$iod_prof) &
      is.finite(df$cv2_prof) & !is.na(df$cv2_prof) &
      is.finite(df$delta_prof) & !is.na(df$delta_prof),
    ,
    drop = FALSE
  ]
  y_rng <- range(c(df$iod_prof, df$cv2_prof), na.rm = TRUE)
  p1 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = iod_prof, colour = "IOD(r)"), linewidth = 0.9) +
    geom_line(aes(y = cv2_prof, colour = "CV²(r)"), linewidth = 0.9) +
    geom_vline(xintercept = cutoff_rank, linetype = "dashed", linewidth = 0.8, colour = pal$treatment) +
    scale_colour_manual(values = c("IOD(r)" = pal$treatment, "CV²(r)" = pal$control)) +
    coord_cartesian(ylim = y_rng) +
    labs(title = "A. Rank-linked regime profiles", x = "Rank", y = "amplitude", colour = NULL) +
    theme_seq()

  p2 <- ggplot(df, aes(rank, iod_prof)) +
    geom_line(linewidth = 0.9, colour = pal$treatment) +
    geom_vline(xintercept = cutoff_rank, linetype = "dashed", linewidth = 0.8, colour = pal$treatment) +
    labs(title = "B. IOD as a function of rank", x = "Rank", y = "IOD(r)",
         subtitle = "IOD(r) = Var(Y[r]) / mu[r] = 1 + alpha[r] * mu[r]") +
    theme_seq()

  p3 <- ggplot(df, aes(rank, cv2_prof)) +
    geom_line(linewidth = 0.9, colour = pal$treatment) +
    geom_vline(xintercept = cutoff_rank, linetype = "dashed", linewidth = 0.8, colour = pal$treatment) +
    labs(title = "C. CV² as a function of rank", x = "Rank", y = "CV²(r)",
         subtitle = "CV²(r) = Var(Y[r]) / mu[r]^2 = 1 / mu[r] + alpha[r]") +
    theme_seq()

  p4 <- ggplot(df, aes(rank, delta_prof)) +
    geom_line(linewidth = 0.9, colour = pal$treatment) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, colour = pal$control) +
    geom_vline(xintercept = cutoff_rank, linetype = "dashed", linewidth = 0.8, colour = pal$treatment) +
    labs(title = "D. Crossing condition", x = "Rank", y = "IOD(r) - CV²(r)",
         subtitle = "Δ(r) = IOD(r) - CV²(r); crossing occurs where Δ(r) = 0") +
    theme_seq()

  arrangeGrob(
    p1, p2, p3, p4, ncol = 2,
    top = textGrob(paste0(cmp_name, " | cutoff-explanation panel"), gp = gpar(fontface = "bold", cex = 1.04)),
    bottom = textGrob("Dashed vertical line = projected cutoff across all panels.", gp = gpar(cex = 0.9))
  )
}

# -----------------------------------------------------------------------------
# data
# -----------------------------------------------------------------------------

if (!file.exists(count_file)) stop(sprintf("count file not found: %s", count_file), call. = FALSE)

wt <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
wt$OrigID <- as.character(wt$OrigID)
wt$Symbol <- as.character(wt$Symbol)
assert_cols(wt, c("OrigID", "Symbol"), "count file")
assert_cols(wt, meta_all$id, "count file")
wt <- wt[!is.na(wt$OrigID) & !is.na(wt$Symbol), , drop = FALSE]
wt <- wt[rowSums(is.na(wt[, meta_all$id, drop = FALSE])) == 0, , drop = FALSE]
rownames(wt) <- wt$OrigID

annot <- unique(wt[, c("OrigID", "Symbol"), drop = FALSE])
names(annot) <- c("feature_id", "gene_symbol")
annot$feature_id <- as.character(annot$feature_id)
annot$gene_symbol <- as.character(annot$gene_symbol)
annot <- annot %>%
  mutate(gene_symbol = if_else(is.na(gene_symbol), "", trimws(gene_symbol))) %>%
  arrange(feature_id, desc(gene_symbol != ""), gene_symbol) %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  mutate(gene_symbol = na_if(gene_symbol, ""))

# -----------------------------------------------------------------------------
# DESeq2 + Strimmer
# -----------------------------------------------------------------------------

get_coef <- function(dds) {
  rn <- resultsNames(dds)
  idx <- grep("^condition_", rn)
  if (!length(idx)) stop("condition coefficient not found", call. = FALSE)
  rn[idx[1]]
}

run_fdrtool_strict <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  stat_vec <- unname(stat_vec)

  if (length(stat_vec) < 5) {
    stop(sprintf("[%s] too few finite Wald statistics for fdrtool", dataset_name), call. = FALSE)
  }
  sdx <- stats::sd(stat_vec, na.rm = TRUE)
  if (!is.finite(sdx) || is.na(sdx) || sdx == 0) {
    stop(sprintf("[%s] Wald statistics have zero/invalid spread for fdrtool", dataset_name), call. = FALSE)
  }

  tries <- list(
    list(cutoff.method = "fndr", pct0 = 0.75),
    list(cutoff.method = "pct0", pct0 = 0.75),
    list(cutoff.method = "locfdr", pct0 = 0.75)
  )

  last_err <- NULL
  for (tt in tries) {
    fit <- tryCatch(
      fdrtool(stat_vec, statistic = "normal", plot = FALSE, verbose = FALSE,
              cutoff.method = tt$cutoff.method, pct0 = tt$pct0),
      error = function(e) { last_err <<- conditionMessage(e); NULL }
    )
    if (is.null(fit)) next
    if (!all(c("pval", "qval", "lfdr") %in% names(fit))) next
    if (length(fit$pval) != length(stat_vec) || length(fit$qval) != length(stat_vec) || length(fit$lfdr) != length(stat_vec)) next
    fit$pval <- clip_p(fit$pval)
    fit$qval <- clip_p(fit$qval)
    fit$lfdr <- as.numeric(fit$lfdr)
    return(fit)
  }
  stop(sprintf("[%s] fdrtool failed under strict Strimmer calls. Last error: %s", dataset_name,
               ifelse(is.null(last_err), "unknown", last_err)), call. = FALSE)
}

hc_thresh_strict <- function(emp_p, dataset_name) {
  x <- sort(clip_p(emp_p), na.last = NA, decreasing = FALSE)
  if (length(x) < 5) return(NA_real_)
  out <- suppressWarnings(tryCatch(fdrtool::hc.thresh(as.vector(x)),
                                   error = function(e) { message(sprintf("[%s] hc.thresh failed: %s", dataset_name, conditionMessage(e))); NA_real_ }))
  out <- as.numeric(out[1])
  if (!is.finite(out) || is.na(out) || out <= 0 || out >= 1) return(NA_real_)
  out
}

run_core <- function(count_mat, coldata, dataset_name, annot_df) {
  dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, betaPrior = FALSE)

  res <- results(dds, contrast = c("condition", "trt", "untrt"), alpha = alpha_level, independentFiltering = FALSE)
  res_ga <- results(dds, contrast = c("condition", "trt", "untrt"), lfcThreshold = lfc_boundary,
                    altHypothesis = "greaterAbs", independentFiltering = FALSE)
  res_la <- results(dds, contrast = c("condition", "trt", "untrt"), lfcThreshold = lfc_boundary,
                    altHypothesis = "lessAbs", independentFiltering = FALSE)

  res_df <- as.data.frame(res, stringsAsFactors = FALSE)
  res_df$feature_id <- rownames(res_df)

  stat_ok <- is.finite(res_df$stat) & !is.na(res_df$stat)
  fit <- run_fdrtool_strict(res_df$stat[stat_ok], dataset_name)

  res_df$empirical_p <- NA_real_
  res_df$empirical_q <- NA_real_
  res_df$lfdr <- NA_real_
  res_df$empirical_p[stat_ok] <- as.numeric(fit$pval)
  res_df$empirical_q[stat_ok] <- as.numeric(fit$qval)
  res_df$lfdr[stat_ok] <- as.numeric(fit$lfdr)

  coef_name <- get_coef(dds)
  shr <- if (identical(lfc_shrink_type, "normal")) {
    suppressMessages(suppressWarnings(lfcShrink(dds, coef = coef_name, type = "normal")))
  } else {
    suppressMessages(suppressWarnings(lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)))
  }
  shr_df <- as.data.frame(shr, stringsAsFactors = FALSE)
  shr_df$feature_id <- rownames(shr_df)
  res_df <- left_join(res_df, shr_df[, c("feature_id", "log2FoldChange")], by = "feature_id", suffix = c("", "_shr"))
  names(res_df)[names(res_df) == "log2FoldChange_shr"] <- "lfc_shrunk"

  hc <- hc_thresh_strict(res_df$empirical_p, dataset_name)
  emp_floor <- ifelse(is.na(res_df$empirical_p), NA_real_, pmax(res_df$empirical_p, 1e-300))
  res_df$HBFSS <- abs(res_df$lfc_shrunk * log10(emp_floor))

  if (is.na(hc)) {
    hbfss_thr <- NA_real_
    res_df$hc_empirical_pass <- FALSE
    res_df$HBFSS_significant <- FALSE
  } else {
    hbfss_thr <- abs(log10(hc)) * lfc_boundary
    res_df$hc_empirical_pass <- !is.na(res_df$empirical_p) & (res_df$empirical_p <= hc)
    res_df$HBFSS_significant <- ifelse(is.na(res_df$HBFSS), FALSE, res_df$HBFSS >= hbfss_thr)
  }

  ga_df <- as.data.frame(res_ga, stringsAsFactors = FALSE); ga_df$feature_id <- rownames(ga_df)
  la_df <- as.data.frame(res_la, stringsAsFactors = FALSE); la_df$feature_id <- rownames(la_df)
  res_df <- left_join(res_df, ga_df[, c("feature_id", "padj")], by = "feature_id", suffix = c("", "_ga"))
  res_df <- left_join(res_df, la_df[, c("feature_id", "padj")], by = "feature_id", suffix = c("", "_la"))
  names(res_df)[names(res_df) == "padj_ga"] <- "resGA_padj"
  names(res_df)[names(res_df) == "padj_la"] <- "resLA_padj"

  res_df$deseq2_strong_call <- !is.na(res_df$resGA_padj) & (res_df$resGA_padj < alpha_level) &
    !is.na(res_df$lfc_shrunk) & (abs(res_df$lfc_shrunk) >= lfc_boundary)

  res_df$deseq2_weak_call <- !is.na(res_df$resLA_padj) & (res_df$resLA_padj < alpha_level) &
    !is.na(res_df$lfc_shrunk) & (abs(res_df$lfc_shrunk) < lfc_boundary) &
    res_df$hc_empirical_pass

  res_df$standard_significant <- res_df$deseq2_strong_call | res_df$deseq2_weak_call
  res_df$HBFSS_only_call <- res_df$HBFSS_significant & !res_df$standard_significant
  res_df$overlap_call <- res_df$HBFSS_significant & res_df$deseq2_strong_call

  res_df$regulation_direction <- ifelse(is.na(res_df$lfc_shrunk), NA_character_,
                                        ifelse(res_df$lfc_shrunk > 0, "upregulated",
                                               ifelse(res_df$lfc_shrunk < 0, "downregulated", "no_change")))

  res_df$neglog10_wald <- safe_nlog10(res_df$pvalue)
  res_df$neglog10_emp <- safe_nlog10(res_df$empirical_p)

  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts$feature_id <- rownames(norm_counts)

  mm <- as.data.frame(mcols(dds), stringsAsFactors = FALSE)
  mm$feature_id <- rownames(mm)
  keep_mm <- intersect(c("feature_id", "dispGeneEst", "dispFit", "dispersion", "baseMean"), names(mm))
  mm <- mm[, keep_mm, drop = FALSE]

  out <- res_df %>%
    left_join(annot_df, by = "feature_id") %>%
    left_join(norm_counts, by = "feature_id") %>%
    left_join(mm, by = "feature_id")

  out$dataset_name <- dataset_name
  out$hc_p_threshold_dataset <- hc
  out$hbfss_threshold_dataset <- hbfss_thr

  keep <- c("dataset_name","feature_id","gene_symbol","baseMean","stat","lfc_shrunk","pvalue","padj",
            "empirical_p","empirical_q","lfdr","resGA_padj","resLA_padj","deseq2_strong_call",
            "deseq2_weak_call","standard_significant","HBFSS_significant","HBFSS_only_call",
            "overlap_call","hc_p_threshold_dataset","hbfss_threshold_dataset","HBFSS","regulation_direction",
            "neglog10_wald","neglog10_emp","dispGeneEst","dispFit","dispersion")
  out <- out[, c(intersect(keep, names(out)), setdiff(names(out), keep)), drop = FALSE]

  list(dds = dds, results = out, hc = hc, hbfss_thr = hbfss_thr)
}

# -----------------------------------------------------------------------------
# EVS rank + Fourier
# -----------------------------------------------------------------------------

prep_cmp <- function(cmp, trt_prefix, ctrl_prefix) {
  keep_ids <- grepl(paste0("^", trt_prefix, "_"), meta_all$id) | grepl(paste0("^", ctrl_prefix, "_"), meta_all$id)
  meta_sub <- meta_all[keep_ids, , drop = FALSE]
  coldata <- meta_sub[, "condition", drop = FALSE]
  sample_ids <- rownames(meta_sub)
  count_mat <- as.matrix(wt[, sample_ids, drop = FALSE])
  list(cmp = cmp, count_mat = count_mat, coldata = coldata)
}

cond_feature_metrics <- function(count_sub) {
  cd <- S4Vectors::DataFrame(row.names = colnames(count_sub))
  dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(count_sub)), colData = cd, design = ~ 1)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds, quiet = TRUE)
  dds <- estimateDispersionsFit(dds, quiet = TRUE)
  md <- as.data.frame(SummarizedExperiment::mcols(dds), stringsAsFactors = FALSE)
  md$feature_id <- rownames(md)
  md[, intersect(c("feature_id","baseMean","dispGeneEst","dispFit","dispersion"), names(md)), drop = FALSE]
}

pc1_tbl <- function(value_df, sample_names, metrics = NULL) {
  x <- as.matrix(value_df[, sample_names, drop = FALSE])
  pca <- prcomp(t(x), scale. = FALSE, rank. = min(5L, ncol(x)))
  la <- abs(pca$rotation[, 1])
  out <- data.frame(feature_id = names(la), pc1_loading = unname(pca$rotation[,1]), pc1_abs = unname(la), stringsAsFactors = FALSE)
  if (!is.null(metrics) && nrow(metrics)) out <- left_join(out, metrics, by = "feature_id")
  out <- out[order(out$pc1_abs, decreasing = TRUE), , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  list(pca = pca, tbl = out)
}

fourier_design <- function(x, harmonics = fourier_harmonics) {
  x <- as.numeric(x)
  x01 <- (x - min(x)) / max(1e-12, (max(x) - min(x)))
  out <- data.frame(x = x01)
  for (k in seq_len(harmonics)) {
    out[[paste0("sin_", k)]] <- sin(2*pi*k*x01)
    out[[paste0("cos_", k)]] <- cos(2*pi*k*x01)
  }
  out
}

fit_local_fourier <- function(rank_vec, y_vec, harmonics = fourier_harmonics) {
  rank_vec <- as.numeric(rank_vec)
  y_vec <- as.numeric(y_vec)
  keep <- is.finite(rank_vec) & !is.na(rank_vec) & is.finite(y_vec) & !is.na(y_vec)
  rank_vec <- rank_vec[keep]
  y_vec <- y_vec[keep]
  if (length(rank_vec) < (2 * harmonics + 5L)) return(NULL)
  y_sd <- suppressWarnings(stats::sd(y_vec, na.rm = TRUE))
  if (!is.finite(y_sd) || is.na(y_sd) || y_sd == 0) return(NULL)
  dd <- fourier_design(rank_vec, harmonics)
  dd$y <- y_vec
  rhs <- paste(setdiff(names(dd), "y"), collapse = " + ")
  fm <- as.formula(paste("y ~", rhs))
  fit <- tryCatch(stats::lm(fm, data = dd), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  fy <- as.numeric(predict(fit, newdata = dd))
  amp <- 0.5 * (max(fy, na.rm = TRUE) - min(fy, na.rm = TRUE))
  pk <- which.max(fy); tr <- which.min(fy)
  list(local_amplitude = amp,
       center_rank = mean(c(rank_vec[pk], rank_vec[tr])),
       residual_sd = sd(dd$y - fy, na.rm = TRUE))
}

rank_metric_tbl <- function(tbl) {
  need <- c("feature_id","rank","pc1_abs","baseMean","dispGeneEst")
  names(tbl)[names(tbl) == "pc1_abs"] <- "pc1_abs"
  if (!all(need %in% names(tbl))) return(NULL)
  df <- tbl[, need, drop = FALSE]
  df <- df[is.finite(df$baseMean) & !is.na(df$baseMean) & df$baseMean > 0 &
             is.finite(df$dispGeneEst) & !is.na(df$dispGeneEst) & df$dispGeneEst > 0, , drop = FALSE]
  if (nrow(df) < 25) return(NULL)
  df <- df[order(df$rank), , drop = FALSE]
  mu <- pmax(df$baseMean, 1e-12)
  alpha <- pmax(df$dispGeneEst, 1e-12)
  df$iod <- 1 + alpha * mu
  df$cv2 <- (1 / mu) + alpha
  df$log_iod <- log10(df$iod)
  df$log_cv2 <- log10(df$cv2)
  df
}

build_windows <- function(n_total) {
  pct <- seq(fourier_step, 1, by = fourier_step)
  ctr <- pmax(1L, pmin(n_total, round(pct * n_total)))
  max_rank <- floor(n_total * 1.00)
  keep <- ctr >= 1L & ctr <= max_rank
  pct <- pct[keep]; ctr <- ctr[keep]
  half <- max(5L, round((fourier_window_frac * n_total) / 2))
  data.frame(percentile = pct, center_rank = ctr,
             lo_rank = pmax(1L, ctr - half),
             hi_rank = pmin(n_total, ctr + half),
             stringsAsFactors = FALSE)
}

summarize_window <- function(metric_df, lo_rank, hi_rank) {
  sub <- metric_df[metric_df$rank >= lo_rank & metric_df$rank <= hi_rank, , drop = FALSE]
  if (nrow(sub) < 9L) return(data.frame(iod_amp = NA_real_, cv2_amp = NA_real_,
                                        iod_center = NA_real_, cv2_center = NA_real_, stringsAsFactors = FALSE))
  iod_fit <- fit_local_fourier(sub$rank, sub$log_iod)
  cv2_fit <- fit_local_fourier(sub$rank, sub$log_cv2)
  if (is.null(iod_fit) || is.null(cv2_fit)) return(data.frame(iod_amp = NA_real_, cv2_amp = NA_real_,
                                                              iod_center = NA_real_, cv2_center = NA_real_, stringsAsFactors = FALSE))
  data.frame(iod_amp = iod_fit$local_amplitude,
             cv2_amp = cv2_fit$local_amplitude,
             iod_center = iod_fit$center_rank,
             cv2_center = cv2_fit$center_rank,
             stringsAsFactors = FALSE)
}

wave_map <- function(tbl) {
  metric_df <- rank_metric_tbl(tbl)
  if (is.null(metric_df) || !nrow(metric_df)) return(NULL)
  win <- build_windows(nrow(metric_df))
  rows <- lapply(seq_len(nrow(win)), function(i) cbind(win[i, , drop = FALSE], summarize_window(metric_df, win$lo_rank[i], win$hi_rank[i])))
  wm <- bind_rows(rows)
  iod_s <- if (all(is.na(wm$iod_amp))) rep(NA_real_, nrow(wm)) else scales::rescale(wm$iod_amp, to = c(0,1), from = range(wm$iod_amp, na.rm = TRUE))
  cv2_s <- if (all(is.na(wm$cv2_amp))) rep(NA_real_, nrow(wm)) else scales::rescale(wm$cv2_amp, to = c(0,1), from = range(wm$cv2_amp, na.rm = TRUE))
  wm$center_dist <- abs(wm$iod_center - wm$cv2_center)
  wm$center_agree <- 1 / (1 + wm$center_dist)
  wm$score <- iod_s + cv2_s + wm$center_agree
  list(metric_df = metric_df, wave = wm)
}

combine_waves <- function(trt_wave, ctrl_wave) {
  if (is.null(trt_wave) || is.null(ctrl_wave)) return(NULL)
  a <- trt_wave$wave; b <- ctrl_wave$wave
  if (is.null(a) || is.null(b) || !nrow(a) || !nrow(b)) return(NULL)
  keep <- c("percentile","center_rank","iod_amp","cv2_amp","iod_center","cv2_center","score")
  a2 <- a[, keep, drop = FALSE]; b2 <- b[, keep, drop = FALSE]
  names(a2)[names(a2) != "percentile"] <- paste0(names(a2)[names(a2) != "percentile"], "_trt")
  names(b2)[names(b2) != "percentile"] <- paste0(names(b2)[names(b2) != "percentile"], "_ctrl")
  out <- inner_join(a2, b2, by = "percentile")
  out$combined_center_rank <- round((out$center_rank_trt + out$center_rank_ctrl)/2)
  out$combined_iod_amp <- out$iod_amp_trt + out$iod_amp_ctrl
  out$combined_cv2_amp <- out$cv2_amp_trt + out$cv2_amp_ctrl
  out$combined_center_dist <- abs(((out$iod_center_trt + out$iod_center_ctrl)/2) - ((out$cv2_center_trt + out$cv2_center_ctrl)/2))
  out$combined_center_agree <- 1 / (1 + out$combined_center_dist)
  iod_s <- if (all(is.na(out$combined_iod_amp))) rep(NA_real_, nrow(out)) else scales::rescale(out$combined_iod_amp, to = c(0,1), from = range(out$combined_iod_amp, na.rm = TRUE))
  cv2_s <- if (all(is.na(out$combined_cv2_amp))) rep(NA_real_, nrow(out)) else scales::rescale(out$combined_cv2_amp, to = c(0,1), from = range(out$combined_cv2_amp, na.rm = TRUE))
  out$combined_score <- iod_s + cv2_s + out$combined_center_agree
  out$diff <- out$combined_iod_amp - out$combined_cv2_amp
  out
}

find_crossings <- function(df) {
  df <- df[order(df$percentile), , drop = FALSE]
  df <- df[is.finite(df$percentile) & is.finite(df$diff) &
             df$percentile >= crossing_min_pct & df$percentile <= crossing_max_pct, , drop = FALSE]
  if (nrow(df) < 2L) return(data.frame())
  rows <- list(); kk <- 1L
  for (i in seq_len(nrow(df)-1L)) {
    y1 <- df$diff[i]; y2 <- df$diff[i+1L]
    crossed <- (y1 == 0) || (y2 == 0) || ((y1 > 0) && (y2 < 0)) || ((y1 < 0) && (y2 > 0))
    if (!crossed) next
    x1 <- df$percentile[i]; x2 <- df$percentile[i+1L]
    r1 <- df$combined_center_rank[i]; r2 <- df$combined_center_rank[i+1L]
    if (isTRUE(all.equal(y1, y2))) {
      cp <- mean(c(x1,x2)); cr <- round(mean(c(r1,r2)))
    } else {
      cp <- x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
      cr <- round(r1 + (0 - y1) * (r2 - r1) / (y2 - y1))
    }
    rows[[kk]] <- data.frame(
      crossing_id = paste0("crossing_", kk),
      idx_left = i, idx_right = i+1L,
      percentile_left = x1, percentile_right = x2,
      rank_left = r1, rank_right = r2,
      diff_left = y1, diff_right = y2,
      crossing_percentile = cp, crossing_rank = as.integer(cr),
      stringsAsFactors = FALSE
    )
    kk <- kk + 1L
  }
  if (!length(rows)) return(data.frame())
  bind_rows(rows)
}

label_stable <- function(diff_df, crossing_tbl) {
  if (is.null(crossing_tbl) || !nrow(crossing_tbl)) return(data.frame())
  out <- crossing_tbl
  out$left_pos_frac <- NA_real_; out$right_neg_frac <- NA_real_; out$stable <- FALSE
  for (i in seq_len(nrow(out))) {
    il <- out$idx_left[i]; ir <- out$idx_right[i]
    left_idx <- seq.int(max(1L, il - crossing_stability_n + 1L), il, by = 1L)
    right_idx <- seq.int(ir, min(nrow(diff_df), ir + crossing_stability_n - 1L), by = 1L)
    lp <- mean(diff_df$diff[left_idx] > 0, na.rm = TRUE)
    rn <- mean(diff_df$diff[right_idx] < 0, na.rm = TRUE)
    out$left_pos_frac[i] <- lp
    out$right_neg_frac[i] <- rn
    out$stable[i] <- isTRUE(is.finite(lp) && is.finite(rn) && lp >= 0.67 && rn >= 0.67)
  }
  out
}

last_stable_before_div <- function(diff_df, crossing_tbl) {
  st <- crossing_tbl[crossing_tbl$stable, , drop = FALSE]
  if (!nrow(st)) return(NULL)
  st <- st[order(st$crossing_percentile), , drop = FALSE]
  keep <- rep(FALSE, nrow(st))
  for (i in seq_len(nrow(st))) {
    ir <- st$idx_right[i]
    idx <- seq.int(ir, min(nrow(diff_df), ir + crossing_post_sep_n - 1L), by = 1L)
    neg_frac <- mean(diff_df$diff[idx] < 0, na.rm = TRUE)
    keep[i] <- isTRUE(is.finite(neg_frac) && neg_frac >= crossing_post_sep_frac)
  }
  st2 <- st[keep, , drop = FALSE]
  if (!nrow(st2)) return(st[nrow(st), , drop = FALSE])
  st2[nrow(st2), , drop = FALSE]
}

backup_crossing <- function(wave_obj, label = "group") {
  if (is.null(wave_obj) || is.null(wave_obj$wave) || !nrow(wave_obj$wave)) {
    return(list(selected = NULL, reason = "wave_missing", table = data.frame(), diff_df = data.frame()))
  }
  df <- wave_obj$wave
  diff_df <- data.frame(percentile = df$percentile,
                        combined_center_rank = round(df$center_rank),
                        diff = df$iod_amp - df$cv2_amp, stringsAsFactors = FALSE)
  tab <- find_crossings(diff_df)
  tab <- label_stable(diff_df, tab)
  if (!nrow(tab)) return(list(selected = NULL, reason = paste0(label, "_no_crossing"), table = tab, diff_df = diff_df))
  st <- tab[tab$stable, , drop = FALSE]
  if (nrow(st)) {
    st <- st[order(st$crossing_percentile), , drop = FALSE]
    sel <- st[nrow(st), , drop = FALSE]
    list(selected = sel, reason = paste0(label, "_last_stable"), table = tab, diff_df = diff_df)
  } else {
    tab <- tab[order(tab$crossing_percentile), , drop = FALSE]
    list(selected = tab[nrow(tab), , drop = FALSE], reason = paste0(label, "_last_crossing_fallback"), table = tab, diff_df = diff_df)
  }
}

# -----------------------------------------------------------------------------
# raw-Wald bracket search
# -----------------------------------------------------------------------------

wald_dist_metrics <- function(z) {
  z <- as.numeric(z)
  z <- z[is.finite(z) & !is.na(z)]
  if (length(z) < 5) return(data.frame(n = length(z), mean_abs = Inf, sd_dev = Inf, skew_abs = Inf, distortion = Inf))
  mu <- mean(z); sdv <- sd(z); sk <- safe_skew(z)
  mean_abs <- abs(mu); sd_dev <- abs(sdv - 1); skew_abs <- abs(sk)
  distortion <- raw_wald_w_mean*mean_abs + raw_wald_w_sd*sd_dev + raw_wald_w_skew*skew_abs
  data.frame(n = length(z), mean_abs = mean_abs, sd_dev = sd_dev, skew_abs = skew_abs, distortion = distortion)
}

build_raw_wald_tbl <- function(full_results, shared_tbl) {
  res <- as.data.frame(full_results, stringsAsFactors = FALSE)
  assert_cols(res, c("feature_id","stat"), "full_results")
  assert_cols(shared_tbl, c("feature_id","combined_rank"), "shared_tbl")
  out <- left_join(shared_tbl[, c("feature_id","combined_rank"), drop = FALSE],
                   res[, c("feature_id","stat"), drop = FALSE], by = "feature_id")
  names(out)[names(out) == "stat"] <- "wald_stat"
  out <- out[order(out$combined_rank), , drop = FALSE]
  out <- out[is.finite(out$wald_stat) & !is.na(out$wald_stat), , drop = FALSE]
  rownames(out) <- NULL
  out
}

score_raw_split <- function(rank_index, ordered_wald_df) {
  n_total <- nrow(ordered_wald_df)
  rank_index <- as.integer(rank_index)
  if (!is.finite(rank_index) || is.na(rank_index)) return(NULL)
  if (rank_index < raw_wald_min_n) return(NULL)
  if ((n_total - rank_index) < raw_wald_min_n) return(NULL)

  z_lead <- ordered_wald_df$wald_stat[seq_len(rank_index)]
  z_rem <- ordered_wald_df$wald_stat[(rank_index + 1L):n_total]
  m_lead <- wald_dist_metrics(z_lead)
  m_rem <- wald_dist_metrics(z_rem)
  if (!is.finite(m_lead$distortion[1]) || !is.finite(m_rem$distortion[1])) return(NULL)
  total_distortion <- m_lead$distortion[1] + m_rem$distortion[1]
  data.frame(
    rank_index = rank_index,
    n_lead = rank_index,
    n_rem = n_total - rank_index,
    lead_mean_abs = m_lead$mean_abs[1],
    lead_sd_dev = m_lead$sd_dev[1],
    lead_skew_abs = m_lead$skew_abs[1],
    lead_distortion = m_lead$distortion[1],
    rem_mean_abs = m_rem$mean_abs[1],
    rem_sd_dev = m_rem$sd_dev[1],
    rem_skew_abs = m_rem$skew_abs[1],
    rem_distortion = m_rem$distortion[1],
    total_distortion = total_distortion,
    split_score = -total_distortion,
    stringsAsFactors = FALSE
  )
}

run_raw_scan <- function(full_results, shared_tbl, lower_rank, upper_rank) {
  ow <- build_raw_wald_tbl(full_results, shared_tbl)
  n_total <- nrow(ow)
  lo <- max(raw_wald_min_n, as.integer(lower_rank))
  hi <- min(n_total - raw_wald_min_n, as.integer(upper_rank))
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) stop("invalid raw-Wald search interval", call. = FALSE)

  coarse <- seq.int(lo, hi, by = max(1L, as.integer(raw_wald_coarse_step)))
  if (tail(coarse,1) != hi) coarse <- c(coarse, hi)
  coarse_tbl <- bind_rows(Filter(Negate(is.null), lapply(coarse, score_raw_split, ordered_wald_df = ow)))
  if (!nrow(coarse_tbl)) stop("no valid coarse candidates", call. = FALSE)
  best_coarse <- coarse_tbl$rank_index[which.max(coarse_tbl$split_score)]

  fine_lo <- max(lo, best_coarse - as.integer(raw_wald_fine_radius))
  fine_hi <- min(hi, best_coarse + as.integer(raw_wald_fine_radius))
  fine_tbl <- bind_rows(Filter(Negate(is.null), lapply(seq.int(fine_lo, fine_hi, by = 1L), score_raw_split, ordered_wald_df = ow)))
  if (!nrow(fine_tbl)) stop("no valid fine candidates", call. = FALSE)
  best <- fine_tbl[which.max(fine_tbl$split_score), , drop = FALSE]

  list(ordered_wald = ow, coarse = coarse_tbl, fine = fine_tbl, selected_rank = as.integer(best$rank_index[1]), selected_row = best)
}

# -----------------------------------------------------------------------------
# volcano
# -----------------------------------------------------------------------------

vol_cols <- c(
  "Neither" = "grey70",
  "DC2 weak effect" = pal$weak,
  "DC2 strong effect" = pal$strong,
  "HBFSS only" = pal$hbfss,
  "Overlap" = pal$overlap
)
vol_shapes <- c(
  "Neither" = 21,
  "DC2 weak effect" = 22,
  "DC2 strong effect" = 24,
  "HBFSS only" = 23,
  "Overlap" = 25
)

build_volcano_classes <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  visible <- is.finite(df$neglog10_emp) & !is.na(df$neglog10_emp)
  hc_thr <- as.numeric(df$hc_p_threshold_dataset[1])
  if (is.finite(hc_thr) && !is.na(hc_thr) && hc_thr > 0 && hc_thr < 1) {
    visible <- visible & (df$neglog10_emp >= -log10(hc_thr))
  }

  strong <- !is.na(df$deseq2_strong_call) & df$deseq2_strong_call
  weak <- !is.na(df$deseq2_weak_call) & df$deseq2_weak_call
  hbfss_only <- !is.na(df$HBFSS_only_call) & df$HBFSS_only_call
  overlap <- !is.na(df$overlap_call) & df$overlap_call

  df$plot_class <- dplyr::case_when(
    overlap & visible ~ "Overlap",
    hbfss_only & visible ~ "HBFSS only",
    strong & visible ~ "DC2 strong effect",
    weak & visible ~ "DC2 weak effect",
    TRUE ~ "Neither"
  )
  df$plot_class <- factor(df$plot_class, levels = c("Neither","DC2 weak effect","DC2 strong effect","HBFSS only","Overlap"))
  df$gene_symbol_plot <- ifelse(!is.na(df$gene_symbol) & grepl("^[A-Za-z0-9._-]+$", trimws(df$gene_symbol)), trimws(df$gene_symbol), NA_character_)
  df
}

compute_vol_lims <- function(df) {
  x <- as.numeric(df$lfc_shrunk); y <- as.numeric(df$neglog10_emp)
  x <- x[is.finite(x) & !is.na(x)]; y <- y[is.finite(y) & !is.na(y)]
  if (!length(x)) x <- c(-2,2)
  if (!length(y)) y <- c(0,5)
  xq <- quantile(abs(x), probs = 0.995, na.rm = TRUE, type = 7)
  xlim <- max(lfc_boundary + 0.5, as.numeric(xq)); xlim <- min(xlim, max(abs(x), na.rm = TRUE)); xlim <- max(2.5, xlim)
  yq <- quantile(y, probs = 0.995, na.rm = TRUE, type = 7)
  ylim <- max(3, as.numeric(yq) * 1.08)
  list(x = c(-xlim, xlim), y = c(0, ylim))
}

pick_labels <- function(df, n_labels = 10) {
  df <- df[!is.na(df$gene_symbol_plot) & df$plot_class %in% c("Overlap","DC2 strong effect","DC2 weak effect","HBFSS only"), , drop = FALSE]
  if (!nrow(df)) return(df[0, , drop = FALSE])
  df$prio <- case_when(
    df$plot_class == "Overlap" ~ 1L,
    df$plot_class == "DC2 strong effect" ~ 2L,
    df$plot_class == "HBFSS only" ~ 3L,
    df$plot_class == "DC2 weak effect" ~ 4L,
    TRUE ~ 9L
  )
  ord <- order(df$prio, -df$neglog10_emp, -abs(df$lfc_shrunk), na.last = TRUE)
  df <- df[ord, , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol_plot), , drop = FALSE]
  df[seq_len(min(n_labels, nrow(df))), , drop = FALSE]
}

add_hbfss_guides <- function(p, df, lims) {
  thr <- as.numeric(df$hbfss_threshold_dataset[1])
  hc <- as.numeric(df$hc_p_threshold_dataset[1])
  if (is.finite(thr) && !is.na(thr) && thr > 0) {
    x_abs <- seq(from = max(0.08, min(lfc_boundary, max(abs(lims$x)))), to = max(abs(lims$x)), length.out = 400)
    y <- thr / x_abs
    bd <- data.frame(lfc_shrunk = c(-rev(x_abs), x_abs), neglog10_emp = c(rev(y), y))
    p <- p + geom_path(data = bd, aes(x = lfc_shrunk, y = neglog10_emp), inherit.aes = FALSE,
                       linetype = "22", linewidth = 0.55, colour = pal$hbfss, alpha = 0.80)
  }
  if (is.finite(hc) && !is.na(hc) && hc > 0 && hc < 1) {
    p <- p + geom_hline(yintercept = -log10(hc), linetype = "22", linewidth = 0.55, colour = pal$threshold, alpha = 0.80)
  }
  txt <- paste0("HC p = ", ifelse(is.finite(hc), fmt(hc, 4), "NA"),
                "\nHBFSS = ", ifelse(is.finite(thr), fmt(thr, 4), "NA"), " / |LFC|")
  p + annotate("label", x = lims$x[1] * 0.96, y = lims$y[2] * 0.96, label = txt,
               hjust = 0, vjust = 1, fill = "white", colour = "grey25", size = 2.3, label.size = 0.12)
}

plot_volcano <- function(df, panel_title, show_legend = FALSE) {
  df <- build_volcano_classes(df)
  labs_df <- pick_labels(df, n_labels = 10)
  lims <- compute_vol_lims(df)

  p <- ggplot(df, aes(lfc_shrunk, neglog10_emp)) +
    geom_point(aes(fill = plot_class, color = plot_class, shape = plot_class),
               alpha = 0.82, size = 2.0, stroke = 0.3, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_boundary, lfc_boundary), linetype = 'dashed', linewidth = 0.55, colour = pal$threshold) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = 'grey45') +
    scale_fill_manual(values = vol_cols, drop = FALSE, name = NULL) +
    scale_color_manual(values = vol_cols, drop = FALSE, name = NULL) +
    scale_shape_manual(values = vol_shapes, drop = FALSE, name = NULL) +
    coord_cartesian(xlim = lims$x, ylim = lims$y, clip = 'off') +
    labs(title = panel_title,
         subtitle = NULL,
         x = 'Shrunken log2 fold change',
         y = expression(-log[10](p[empirical])),
         caption = paste0(
           'Standard=', sum(df$standard_significant, na.rm = TRUE),
           ' | Weak=', sum(df$deseq2_weak_call, na.rm = TRUE),
           ' | Strong=', sum(df$deseq2_strong_call, na.rm = TRUE),
           ' | HBFSS=', sum(df$HBFSS_significant, na.rm = TRUE),
           ' | Overlap=', sum(df$overlap_call, na.rm = TRUE)
         )) +
    theme_seq() +
    theme(legend.position = if (show_legend) 'bottom' else 'none') +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                               override.aes = list(size = 3.5, alpha = 1,
                                                   shape = unname(vol_shapes),
                                                   fill = unname(vol_cols), colour = unname(vol_cols))),
           color = 'none', shape = 'none')
  p <- add_hbfss_guides(p, df, lims)
  if (nrow(labs_df)) {
    p <- p + geom_text_repel(data = labs_df, aes(label = gene_symbol_plot), size = 1.85,
                             seed = 1, max.overlaps = 18, force = 1.0, force_pull = 0.4,
                             box.padding = 0.25, point.padding = 0.14,
                             min.segment.length = 0, segment.alpha = 0.5, segment.size = 0.18)
  }
  p
}

get_legend_grob <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == 'guide-box')
  if (!length(idx)) return(NULL)
  g$grobs[[idx[1]]]
}

plot_volcano_panel <- function(cmp, raw_df, lead_df, rem_df) {
  p_raw <- plot_volcano(raw_df, paste0(cmp, ' | Original dataset'), show_legend = FALSE)
  p_lead <- plot_volcano(lead_df, paste0(cmp, ' | Leading-edge dataset'), show_legend = FALSE)
  p_rem <- plot_volcano(rem_df, paste0(cmp, ' | Remainder dataset'), show_legend = FALSE)
  leg <- get_legend_grob(plot_volcano(raw_df, paste0(cmp, ' | Original dataset'), show_legend = TRUE))
  arrangeGrob(p_raw, p_lead, p_rem, leg,
              layout_matrix = rbind(c(1,2,3), c(4,4,4)),
              heights = c(12, 1.8),
              top = textGrob(paste0(cmp, ' | HBFSS volcano panel'), gp = gpar(fontface = 'bold', cex = 1.04)))
}

# -----------------------------------------------------------------------------
# figures

# -----------------------------------------------------------------------------

plot_wave_single <- function(wave_obj, cmp, group_label, backup = NULL) {
  if (is.null(wave_obj) || is.null(wave_obj$wave) || !nrow(wave_obj$wave)) return(NULL)
  df <- wave_obj$wave
  line_col <- if (group_label == "Treatment") pal$treatment else pal$control
  p <- ggplot(df, aes(percentile)) +
    geom_line(aes(y = iod_amp, color = "IOD"), linewidth = 0.9) +
    geom_line(aes(y = cv2_amp, color = "CV2"), linewidth = 0.9) +
    scale_color_manual(values = c("IOD" = pal$treatment, "CV2" = pal$control)) +
    labs(title = paste0(group_label, " | local regime lines"),
         x = "Percentile center", y = "Local amplitude", color = NULL) +
    theme_seq()
  if (!is.null(backup) && nrow(backup)) {
    ymax <- max(c(df$iod_amp, df$cv2_amp), na.rm = TRUE)
    p <- p +
      geom_vline(xintercept = backup$crossing_percentile[1], linetype = "dotted", linewidth = 0.9, colour = line_col) +
      annotate("label", x = backup$crossing_percentile[1], y = ymax,
               label = paste0(group_label, " backup\np = ", fmt(backup$crossing_percentile[1], 4),
                              "\nrank = ", backup$crossing_rank[1]),
               fill = "white", colour = line_col, size = 2.7, label.size = 0.15, vjust = -0.5)
  }
  p
}

plot_combined_wave <- function(cmb, cmp, trt_bk = NULL, ctrl_bk = NULL, crossing = NULL, final_rank = NA_integer_) {
  if (is.null(cmb) || !nrow(cmb)) return(NULL)
  p <- ggplot(cmb, aes(percentile)) +
    geom_line(aes(y = combined_iod_amp, color = "Composite IOD"), linewidth = 0.95) +
    geom_line(aes(y = combined_cv2_amp, color = "Composite CV2"), linewidth = 0.95) +
    scale_color_manual(values = c("Composite IOD" = pal$treatment, "Composite CV2" = pal$control)) +
    labs(title = paste0(cmp, " | combined local Fourier wave"),
         x = "Percentile center", y = "Composite local amplitude", color = NULL) +
    theme_seq()
  ymax <- max(c(cmb$combined_iod_amp, cmb$combined_cv2_amp), na.rm = TRUE)
  if (!is.null(crossing) && nrow(crossing)) {
    p <- p + geom_vline(xintercept = crossing$crossing_percentile[1], linetype = "dashed", linewidth = 0.95, colour = pal$threshold) +
      annotate("label", x = crossing$crossing_percentile[1], y = ymax,
               label = paste0("Crossing ref\np = ", fmt(crossing$crossing_percentile[1], 4),
                              "\nrank = ", crossing$crossing_rank[1]),
               fill = "white", colour = pal$threshold, size = 2.8, label.size = 0.15, vjust = -0.55)
  }
  if (!is.null(trt_bk) && nrow(trt_bk)) {
    p <- p + geom_vline(xintercept = trt_bk$crossing_percentile[1], linetype = "dotted", linewidth = 0.85, colour = pal$treatment) +
      annotate("label", x = trt_bk$crossing_percentile[1], y = ymax*0.78,
               label = paste0("Treatment backup\np = ", fmt(trt_bk$crossing_percentile[1], 4),
                              "\nrank = ", trt_bk$crossing_rank[1]),
               fill = "white", colour = pal$treatment, size = 2.6, label.size = 0.15)
  }
  if (!is.null(ctrl_bk) && nrow(ctrl_bk)) {
    p <- p + geom_vline(xintercept = ctrl_bk$crossing_percentile[1], linetype = "dotted", linewidth = 0.85, colour = pal$control) +
      annotate("label", x = ctrl_bk$crossing_percentile[1], y = ymax*0.58,
               label = paste0("Control backup\np = ", fmt(ctrl_bk$crossing_percentile[1], 4),
                              "\nrank = ", ctrl_bk$crossing_rank[1]),
               fill = "white", colour = pal$control, size = 2.6, label.size = 0.15)
  }
  if (is.finite(final_rank) && !is.na(final_rank)) {
    # map final rank to nearest percentile in combined table
    idx <- which.min(abs(cmb$combined_center_rank - final_rank))
    p <- p + geom_vline(xintercept = cmb$percentile[idx], linetype = "solid", linewidth = 1.05, colour = "black") +
      annotate("label", x = cmb$percentile[idx], y = ymax*0.38,
               label = paste0("Final raw-Wald cutoff\nrank = ", final_rank),
               fill = "white", colour = "black", size = 2.7, label.size = 0.15)
  }
  p
}

plot_score_curve <- function(cmb, cmp) {
  if (is.null(cmb) || !nrow(cmb)) return(NULL)
  ggplot(cmb, aes(percentile, combined_score)) +
    geom_col(width = 0.008, fill = pal$threshold) +
    labs(title = paste0(cmp, " | descriptive score profile"),
         x = "Percentile center", y = "Score") +
    theme_seq()
}

plot_fourier_panel <- function(cmp, trt_wave, ctrl_wave, cmb, trt_bk, ctrl_bk, crossing, final_rank) {
  p1 <- plot_wave_single(trt_wave, cmp, "Treatment", trt_bk)
  p2 <- plot_wave_single(ctrl_wave, cmp, "Control", ctrl_bk)
  p3 <- plot_combined_wave(cmb, cmp, trt_bk, ctrl_bk, crossing, final_rank)
  p4 <- plot_score_curve(cmb, cmp)
  arrangeGrob(p1, p2, p3, p4, ncol = 2,
              top = textGrob(paste0(cmp, " | Fourier summary panel"), gp = gpar(fontface = "bold", cex = 1.04)))
}

fmt_num <- function(x, d = 3) {
  if (!is.finite(x) || is.na(x)) return("NA")
  formatC(x, digits = d, format = "fg", flag = "#")
}

rank_to_pct <- function(rank_index, n_total) {
  if (!is.finite(rank_index) || !is.finite(n_total) || n_total <= 0) return(NA_real_)
  rank_index / n_total
}

wald_sum <- function(z) {
  z <- as.numeric(z)
  z <- z[is.finite(z) & !is.na(z)]
  if (length(z) < 3) {
    return(data.frame(n = length(z), mean = NA_real_, sd = NA_real_, skew = NA_real_, stringsAsFactors = FALSE))
  }
  data.frame(
    n = length(z),
    mean = mean(z, na.rm = TRUE),
    sd = stats::sd(z, na.rm = TRUE),
    skew = safe_skew(z),
    stringsAsFactors = FALSE
  )
}

plot_wald_hist <- function(z, title_txt, subtitle_txt = NULL, fill_col = pal$grey) {
  z <- as.numeric(z)
  z <- z[is.finite(z) & !is.na(z)]
  df <- data.frame(wald = z, stringsAsFactors = FALSE)
  sm <- wald_sum(z)
  stat_txt <- paste0(
    "n = ", sm$n[1],
    " | mean = ", fmt_num(sm$mean[1]),
    " | sd = ", fmt_num(sm$sd[1]),
    " | skew = ", fmt_num(sm$skew[1])
  )
  ggplot(df, aes(wald)) +
    geom_histogram(bins = hist_bins, fill = fill_col, color = "white", alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, colour = pal$threshold) +
    labs(
      title = title_txt,
      subtitle = paste(c(subtitle_txt, stat_txt), collapse = "
"),
      x = "Wald statistic",
      y = "Count"
    ) +
    theme_seq() +
    theme(
      plot.title = element_text(size = base_theme_size, hjust = 0.5),
      plot.subtitle = element_text(size = base_theme_size - 1, hjust = 0.5, lineheight = 1.0)
    )
}

plot_raw_wald_search <- function(scan, cmp, trt_rank = NA_integer_, ctrl_rank = NA_integer_) {
  if (is.null(scan) || is.null(scan$fine) || !nrow(scan$fine)) return(NULL)
  fine <- as.data.frame(scan$fine, stringsAsFactors = FALSE)
  coarse <- if (!is.null(scan$coarse)) as.data.frame(scan$coarse, stringsAsFactors = FALSE) else data.frame()
  sel <- as.data.frame(scan$selected_row, stringsAsFactors = FALSE)
  sel_rank <- as.integer(scan$selected_rank)
  n_total <- max(fine$n_leading_edge + fine$n_remainder, na.rm = TRUE)
  sel_pct <- rank_to_pct(sel_rank, n_total)
  trt_pct <- rank_to_pct(trt_rank, n_total)
  ctrl_pct <- rank_to_pct(ctrl_rank, n_total)
  ann_txt <- paste0(
    "Final selected cutoff
",
    "rank = ", sel_rank,
    " | pct = ", fmt_num(sel_pct, 4),
    "
distortion = ", fmt_num(sel$total_distortion[1]),
    "
lead n = ", sel$n_leading_edge[1],
    " | rem n = ", sel$n_remainder[1]
  )
  p <- ggplot(fine, aes(rank_index, total_distortion)) +
    geom_line(linewidth = 0.85, colour = "grey35") +
    geom_point(size = 1.0, colour = "grey35")
  if (nrow(coarse)) {
    p <- p + geom_point(
      data = coarse,
      aes(rank_index, total_distortion),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.6,
      stroke = 0.25,
      fill = "white",
      colour = pal$threshold
    )
  }
  ymax <- max(fine$total_distortion, na.rm = TRUE)
  if (is.finite(trt_rank) && !is.na(trt_rank)) {
    p <- p +
      geom_vline(xintercept = trt_rank, linetype = "dotted", linewidth = 0.9, colour = pal$treatment) +
      annotate(
        "label",
        x = trt_rank,
        y = ymax * 0.98,
        label = paste0("Treatment backup
rank = ", trt_rank, "
pct = ", fmt_num(trt_pct, 4)),
        fill = "white",
        colour = pal$treatment,
        size = 2.6,
        label.size = 0.15,
        vjust = 1
      )
  }
  if (is.finite(ctrl_rank) && !is.na(ctrl_rank)) {
    p <- p +
      geom_vline(xintercept = ctrl_rank, linetype = "dotted", linewidth = 0.9, colour = pal$control) +
      annotate(
        "label",
        x = ctrl_rank,
        y = ymax * 0.88,
        label = paste0("Control backup
rank = ", ctrl_rank, "
pct = ", fmt_num(ctrl_pct, 4)),
        fill = "white",
        colour = pal$control,
        size = 2.6,
        label.size = 0.15,
        vjust = 1
      )
  }
  p +
    geom_vline(xintercept = sel_rank, linetype = "solid", linewidth = 1.1, colour = pal$threshold) +
    geom_point(
      data = sel,
      aes(rank_index, total_distortion),
      inherit.aes = FALSE,
      shape = 24,
      size = 2.6,
      stroke = 0.35,
      fill = "white",
      colour = pal$threshold
    ) +
    annotate(
      "label",
      x = sel_rank,
      y = sel$total_distortion[1],
      label = ann_txt,
      fill = "white",
      colour = pal$threshold,
      size = 2.7,
      label.size = 0.15,
      hjust = 1,
      vjust = 1
    ) +
    labs(
      title = paste0(cmp, " | Bracketed raw-Wald cutoff search"),
      subtitle = paste0(
        "Search performed only between the treatment and control backup cutoffs.
",
        "Total distortion = |mean| + |sd - 1| + 0.5*|skew|. The final selected cutoff minimizes total distortion."
      ),
      x = "Candidate cutoff rank",
      y = "Total distortion"
    ) +
    theme_seq()
}

plot_raw_wald_panel <- function(cmp, scan, trt_rank = NA_integer_, ctrl_rank = NA_integer_) {
  if (is.null(scan) || is.null(scan$ordered_wald) || !nrow(scan$ordered_wald)) return(NULL)
  ow <- as.data.frame(scan$ordered_wald, stringsAsFactors = FALSE)
  sel_rank <- as.integer(scan$selected_rank)
  sel_pct <- rank_to_pct(sel_rank, nrow(ow))
  z_raw <- ow$wald_stat
  z_lead <- ow$wald_stat[seq_len(sel_rank)]
  z_rem <- ow$wald_stat[(sel_rank + 1L):nrow(ow)]
  p1 <- plot_wald_hist(
    z_raw,
    paste0(cmp, " | Raw dataset"),
    "Full unsplit Wald-statistic distribution used as the baseline reference before EVS splitting."
  )
  p2 <- plot_raw_wald_search(scan, cmp, trt_rank, ctrl_rank)
  p3 <- plot_wald_hist(
    z_lead,
    paste0(cmp, " | Leading-edge subset"),
    paste0("Features ranked above the final selected cutoff.
Selected rank = ", sel_rank,
           " | selected percentile = ", fmt_num(sel_pct, 4)),
    pal$treatment
  )
  p4 <- plot_wald_hist(
    z_rem,
    paste0(cmp, " | Remainder subset"),
    paste0("Features ranked below the final selected cutoff.
Selected rank = ", sel_rank,
           " | selected percentile = ", fmt_num(sel_pct, 4)),
    pal$control
  )
  arrangeGrob(
    p1, p2, p3, p4,
    ncol = 2,
    top = textGrob(paste0(cmp, " | Raw-Wald cutoff selection panel"), gp = gpar(fontface = "bold", cex = 1.06)),
    bottom = textGrob(
      paste0(
        "Top left: raw unsplit Wald-statistic distribution. ",
        "Top right: bracketed cutoff scan between treatment and control backup cutoffs. ",
        "Bottom left/right: Wald-statistic distributions for the leading-edge and remainder subsets created by the final selected cutoff."
      ),
      gp = gpar(cex = 0.86)
    )
  )
}

# -----------------------------------------------------------------------------
# pipeline per comparison
# -----------------------------------------------------------------------------

build_split <- function(count_mat, coldata, cmp_name) {
  dds0 <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ condition)
  dds0 <- dds0[rowSums(counts(dds0)) > 0, ]
  dds0 <- estimateSizeFactors(dds0)
  norm_mat <- as.data.frame(counts(dds0, normalized = TRUE))
  raw_mat <- as.data.frame(count_mat[rownames(dds0), , drop = FALSE])

  samp <- colnames(count_mat)
  trt_ids <- samp[coldata$condition == "trt"]
  ctrl_ids <- samp[coldata$condition == "untrt"]

  trt_metrics <- cond_feature_metrics(count_mat[, trt_ids, drop = FALSE])
  ctrl_metrics <- cond_feature_metrics(count_mat[, ctrl_ids, drop = FALSE])

  trt_fit <- pc1_tbl(norm_mat, trt_ids, trt_metrics)
  ctrl_fit <- pc1_tbl(norm_mat, ctrl_ids, ctrl_metrics)
  trt_fit_raw <- pc1_tbl(raw_mat, trt_ids, trt_metrics)
  ctrl_fit_raw <- pc1_tbl(raw_mat, ctrl_ids, ctrl_metrics)

  shared <- full_join(
    trt_fit$tbl[, c("feature_id","pc1_loading","pc1_abs","rank"), drop = FALSE] %>% rename(pc1_loading_trt = pc1_loading, pc1_abs_trt = pc1_abs, rank_trt = rank),
    ctrl_fit$tbl[, c("feature_id","pc1_loading","pc1_abs","rank"), drop = FALSE] %>% rename(pc1_loading_ctrl = pc1_loading, pc1_abs_ctrl = pc1_abs, rank_ctrl = rank),
    by = "feature_id"
  )
  shared$pc1_abs_trt[is.na(shared$pc1_abs_trt)] <- 0
  shared$pc1_abs_ctrl[is.na(shared$pc1_abs_ctrl)] <- 0
  shared$combined_loading <- pmax(shared$pc1_abs_trt, shared$pc1_abs_ctrl, na.rm = TRUE)
  shared <- shared[order(-shared$combined_loading, shared$feature_id), , drop = FALSE]
  shared$combined_rank <- seq_len(nrow(shared))
  rownames(shared) <- NULL

  trt_wave <- wave_map(trt_fit$tbl)
  ctrl_wave <- wave_map(ctrl_fit$tbl)
  cmb <- combine_waves(trt_wave, ctrl_wave)
  cross_tab <- find_crossings(cmb)
  cross_tab <- label_stable(cmb, cross_tab)
  sel_cross <- if (nrow(cross_tab)) last_stable_before_div(cmb, cross_tab) else NULL

  trt_bk <- backup_crossing(trt_wave, "treatment")
  ctrl_bk <- backup_crossing(ctrl_wave, "control")

  final_rank <- if (!is.null(sel_cross) && nrow(sel_cross)) as.integer(sel_cross$crossing_rank[1]) else if (use_fixed_top_n) min(fixed_top_n, nrow(shared)) else min(fixed_top_n, nrow(shared))
  final_method <- if (!is.null(sel_cross) && nrow(sel_cross)) "last_stable_crossing_before_divergence" else "fixed_top_n"

  raw_scan <- NULL
  full_raw_fit <- run_core(count_mat, coldata, paste0(cmp_name, "_raw_wald_reference"), annot)
  if (isTRUE(raw_wald_scan)) {
    trt_rank <- if (!is.null(trt_bk$selected) && nrow(trt_bk$selected)) as.integer(trt_bk$selected$crossing_rank[1]) else final_rank
    ctrl_rank <- if (!is.null(ctrl_bk$selected) && nrow(ctrl_bk$selected)) as.integer(ctrl_bk$selected$crossing_rank[1]) else final_rank
    bracket_lo <- min(trt_rank, ctrl_rank, final_rank, na.rm = TRUE)
    bracket_hi <- max(trt_rank, ctrl_rank, final_rank, na.rm = TRUE)
    pad <- max(10L, round(0.05 * abs(bracket_hi - bracket_lo)))
    bracket_lo <- max(1L, bracket_lo - pad)
    bracket_hi <- min(nrow(shared), bracket_hi + pad)
    raw_scan <- run_raw_scan(full_raw_fit$results, shared, bracket_lo, bracket_hi)
    final_rank <- as.integer(raw_scan$selected_rank)
    final_method <- "raw_wald_bracket_total_distortion"
  }

  lead_ids <- as.character(shared$feature_id[shared$combined_rank <= final_rank])
  all_ids <- as.character(shared$feature_id)
  rem_ids <- setdiff(all_ids, lead_ids)
  if (!length(lead_ids) || !length(rem_ids)) stop(sprintf("[%s] invalid split: empty leading-edge or remainder", cmp_name), call. = FALSE)

  list(
    trt_fit = trt_fit, ctrl_fit = ctrl_fit,
    trt_fit_raw = trt_fit_raw, ctrl_fit_raw = ctrl_fit_raw,
    trt_wave = trt_wave, ctrl_wave = ctrl_wave, cmb = cmb,
    cross_tab = cross_tab, sel_cross = sel_cross,
    trt_bk = trt_bk, ctrl_bk = ctrl_bk,
    shared = shared,
    raw_ref = full_raw_fit,
    raw_scan = raw_scan,
    final_rank = final_rank,
    final_method = final_method,
    raw_dataset = count_mat,
    lead_dataset = count_mat[lead_ids, , drop = FALSE],
    rem_dataset = count_mat[rem_ids, , drop = FALSE]
  )
}

run_cmp <- function(cmp_name, count_mat, coldata) {
  cmp_dir <- file.path(out_dir, cmp_name)
  fig_dir <- file.path(cmp_dir, "figures")
  tab_dir <- file.path(cmp_dir, "tables")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

  log_step(cmp_name, 'comparison loaded')
  sp <- build_split(count_mat, coldata, cmp_name)

  trt_rank <- if (!is.null(sp$trt_bk$selected) && nrow(sp$trt_bk$selected)) as.integer(sp$trt_bk$selected$crossing_rank[1]) else NA_integer_
  ctrl_rank <- if (!is.null(sp$ctrl_bk$selected) && nrow(sp$ctrl_bk$selected)) as.integer(sp$ctrl_bk$selected$crossing_rank[1]) else NA_integer_
  cross_rank <- if (!is.null(sp$sel_cross) && nrow(sp$sel_cross)) as.integer(sp$sel_cross$crossing_rank[1]) else NA_integer_

  cutoff_summary <- data.frame(
    comparison_name = cmp_name,
    final_selected_rank = sp$final_rank,
    method = sp$final_method,
    crossing_rank_reference = cross_rank,
    treatment_backup_rank = trt_rank,
    control_backup_rank = ctrl_rank,
    raw_wald_selected_rank = if (!is.null(sp$raw_scan)) sp$raw_scan$selected_rank else NA_integer_,
    cutoff_quantile = 1 - (sp$final_rank / nrow(sp$shared)),
    stringsAsFactors = FALSE
  )
  save_csv(cutoff_summary, file.path(tab_dir, paste0(cmp_name, "_cutoff_summary.csv")))

  if (!is.null(sp$raw_scan)) {
    save_csv(sp$raw_scan$coarse, file.path(tab_dir, paste0(cmp_name, "_raw_wald_scan_coarse.csv")))
    save_csv(sp$raw_scan$fine, file.path(tab_dir, paste0(cmp_name, "_raw_wald_scan_fine.csv")))
    save_csv(sp$raw_scan$selected_row, file.path(tab_dir, paste0(cmp_name, "_raw_wald_scan_selected.csv")))
  }

  log_step(cmp_name, 'split built')
  log_step(cmp_name, 'cutoff summary | treatment backup = ', trt_rank, ' | control backup = ', ctrl_rank, ' | crossing ref = ', cross_rank, ' | final rank = ', sp$final_rank, ' | method = ', sp$final_method)

  raw_fit <- run_core(sp$raw_dataset, coldata, paste0(cmp_name, '_raw'), annot)
  log_step(cmp_name, 'original dataset analysis complete')
  lead_fit <- run_core(sp$lead_dataset, coldata, paste0(cmp_name, '_lead'), annot)
  log_step(cmp_name, 'leading-edge dataset analysis complete')
  rem_fit <- run_core(sp$rem_dataset, coldata, paste0(cmp_name, '_rem'), annot)
  log_step(cmp_name, 'remainder dataset analysis complete')

  keep_main <- c("dataset_name","feature_id","gene_symbol","baseMean","stat","lfc_shrunk","pvalue","padj","empirical_p",
                 "empirical_q","lfdr","resGA_padj","resLA_padj","deseq2_strong_call","deseq2_weak_call",
                 "standard_significant","HBFSS_significant","HBFSS_only_call","overlap_call",
                 "hc_p_threshold_dataset","hbfss_threshold_dataset","HBFSS","regulation_direction")
  save_csv(raw_fit$results[, intersect(keep_main, names(raw_fit$results)), drop = FALSE],
           file.path(tab_dir, paste0(cmp_name, "_raw_results_main.csv")))
  save_csv(lead_fit$results[, intersect(keep_main, names(lead_fit$results)), drop = FALSE],
           file.path(tab_dir, paste0(cmp_name, "_lead_results_main.csv")))
  save_csv(rem_fit$results[, intersect(keep_main, names(rem_fit$results)), drop = FALSE],
           file.path(tab_dir, paste0(cmp_name, "_rem_results_main.csv")))

  fig_fourier <- plot_fourier_panel(
    cmp_name, sp$trt_wave, sp$ctrl_wave, sp$cmb,
    sp$trt_bk$selected, sp$ctrl_bk$selected, sp$sel_cross, sp$final_rank
  )
  ggsave(file.path(fig_dir, paste0(cmp_name, '_fourier_summary_panel.png')),
         fig_fourier, width = 16, height = 10, dpi = figure_dpi, units = 'in', limitsize = FALSE, bg = 'white')
  log_step(cmp_name, 'Fourier summary panel exported')

  fig_rawscan <- plot_raw_wald_panel(cmp_name, sp$raw_scan, trt_rank, ctrl_rank)
  if (!is.null(fig_rawscan)) {
    ggsave(file.path(fig_dir, paste0(cmp_name, '_raw_wald_cutoff_selection_panel.png')),
           fig_rawscan, width = 16, height = 10, dpi = figure_dpi, units = 'in', limitsize = FALSE, bg = 'white')
    log_step(cmp_name, 'raw-Wald cutoff panel exported')
  }

  fig_vol <- plot_volcano_panel(cmp_name, raw_fit$results, lead_fit$results, rem_fit$results)
  ggsave(file.path(fig_dir, paste0(cmp_name, '_volcano_panel.png')), fig_vol, width = 18, height = 7.8, dpi = figure_dpi, units = 'in', limitsize = FALSE, bg = 'white')
  log_step(cmp_name, 'volcano panel exported')

  data.frame(
    comparison_name = cmp_name,
    final_selected_rank = sp$final_rank,
    method = sp$final_method,
    n_raw = nrow(raw_fit$results),
    n_lead = nrow(lead_fit$results),
    n_rem = nrow(rem_fit$results),
    n_raw_weak = sum(raw_fit$results$deseq2_weak_call, na.rm = TRUE),
    n_raw_strong = sum(raw_fit$results$deseq2_strong_call, na.rm = TRUE),
    n_raw_hbfss = sum(raw_fit$results$HBFSS_significant, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------

all_sum <- list()
fail <- list()

for (i in seq_len(nrow(cmp_tbl))) {
  cmp_name <- cmp_tbl$cmp[i]
  message("\n=====================================================")
  message("Running comparison: ", cmp_name)
  message("=====================================================")
  x <- prep_cmp(cmp_name, cmp_tbl$trt[i], cmp_tbl$ctrl[i])
  out <- tryCatch(
    run_cmp(x$cmp, x$count_mat, x$coldata),
    error = function(e) {
      fail[[cmp_name]] <<- data.frame(comparison_name = cmp_name, error_message = conditionMessage(e), stringsAsFactors = FALSE)
      NULL
    }
  )
  if (!is.null(out) && nrow(out)) all_sum[[cmp_name]] <- out
}

if (length(all_sum)) save_csv(bind_rows(all_sum), file.path(out_dir, "all_comparisons_summary.csv"))
if (length(fail)) save_csv(bind_rows(fail), file.path(out_dir, "failed_comparisons.csv"))

cat("\n=====================================================\n")
cat("Pipeline complete.\n")
cat("Output directory:\n")
cat(normalizePath(out_dir), "\n")
cat("=====================================================\n\n")

session_info_txt <- capture.output(sessionInfo())
writeLines(session_info_txt, file.path(out_dir, "sessionInfo.txt"))
saveRDS(sessionInfo(), file.path(out_dir, "sessionInfo.rds"))
