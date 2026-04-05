# =============================================================================
# SEQUENCE STAGE 1: NEGATIVE-BINOMIAL MEAN-VARIANCE REGIME ANALYSIS
# -----------------------------------------------------------------------------
# Anchored to the working WTTS import / metadata scaffold from SEQUENCE 8_fixed_v6.R
# to avoid the repeated bootstrap failures from earlier rewrites.
#
# Core NB model:
#   Var(X) = mu + alpha * mu^2
#   IOD    = Var(X) / mu   = 1 + alpha * mu
#   CV^2   = Var(X) / mu^2 = 1 / mu + alpha
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(gridExtra)
})

# =============================================================================
# USER SETTINGS
# =============================================================================

repo_dir <- getwd()
input_dir <- file.path(repo_dir, "data")
output_root <- file.path(repo_dir, "exports")
analysis_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(output_root, paste0("nb_regime_analysis_", analysis_stamp))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
percentile_step <- 0.01
min_bin_n <- 5L

# =============================================================================
# MINIMAL EMBEDDED METADATA (copied from working script scaffold)
# =============================================================================

meta_all <- data.frame(
  id = c(
    "R0_1","R0_2","R0_3","R0_4","R0_5","ZT6_1","ZT6_2","ZT6_3","ZT6_4","ZT6_5",
    "R2_1","R2_2","R2_3","R2_4","R2_5","ZT8_1","ZT8_2","ZT8_3","ZT8_4","ZT8_5",
    "R4_1","R4_2","R4_3","R4_4","R4_5","ZT10_1","ZT10_2","ZT10_3","ZT10_4","ZT10_5",
    "R8_1","R8_2","R8_3","R8_4","R8_5","ZT14_1","ZT14_2","ZT14_3","ZT14_4","ZT14_5"
  ),
  condition = c(
    "treatment","treatment","treatment","treatment","treatment",
    "control","control","control","control","control",
    "treatment","treatment","treatment","treatment","treatment",
    "control","control","control","control","control",
    "treatment","treatment","treatment","treatment","treatment",
    "control","control","control","control","control",
    "treatment","treatment","treatment","treatment","treatment",
    "control","control","control","control","control"
  ),
  stringsAsFactors = FALSE
)
rownames(meta_all) <- meta_all$id
meta_all$condition <- factor(meta_all$condition, levels = c("control", "treatment"))
levels(meta_all$condition) <- c("untrt", "trt")

comparison_map <- list(
  RT0_ZT6   = list(group1_prefix = "R0", group2_prefix = "ZT6"),
  RT2_ZT8   = list(group1_prefix = "R2", group2_prefix = "ZT8"),
  RT4_ZT10  = list(group1_prefix = "R4", group2_prefix = "ZT10"),
  RT8_ZT14  = list(group1_prefix = "R8", group2_prefix = "ZT14")
)

# =============================================================================
# HELPERS
# =============================================================================

assert_required_columns <- function(df, cols, object_name = "data frame") {
  missing_cols <- setdiff(cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "%s is missing required columns: %s",
        object_name,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_var <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::var(x)
}

rescale01 <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(out)
  rng <- range(x[ok], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    out[ok] <- 0
    return(out)
  }
  out[ok] <- (x[ok] - rng[1]) / (rng[2] - rng[1])
  out
}

find_last_zero_crossing <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  if (length(x) < 2L) return(NULL)
  idx <- integer(0)
  for (i in seq_len(length(y) - 1L)) {
    y1 <- y[i]
    y2 <- y[i + 1L]
    if (isTRUE(all.equal(y1, 0))) idx <- c(idx, i)
    else if (isTRUE(all.equal(y2, 0))) idx <- c(idx, i)
    else if (is.finite(y1) && is.finite(y2) && y1 * y2 < 0) idx <- c(idx, i)
  }
  if (!length(idx)) return(NULL)
  i <- max(idx)
  x1 <- x[i]; x2 <- x[i + 1L]
  y1 <- y[i]; y2 <- y[i + 1L]
  if (isTRUE(all.equal(y1, 0))) x_cross <- x1
  else if (isTRUE(all.equal(y2, 0))) x_cross <- x2
  else x_cross <- x1 - y1 * (x2 - x1) / (y2 - y1)
  list(index_left = i, crossing_x = x_cross)
}

make_percentile_bins <- function(n, step = 0.01) {
  probs <- seq(step, 1, by = step)
  if (tail(probs, 1) < 1) probs <- c(probs, 1)
  centers <- unique(pmax(1L, pmin(n, round(probs * n))))
  data.frame(percentile = probs[seq_along(centers)], rank_index = centers, stringsAsFactors = FALSE)
}

compute_group_pc1_loadings <- function(norm_counts_mat, group_cols) {
  mat <- norm_counts_mat[, group_cols, drop = FALSE]
  mat <- log2(mat + 1)
  mat <- t(mat)
  if (nrow(mat) < 2L) return(rep(NA_real_, ncol(mat)))
  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  abs(pca$rotation[, 1L])
}

build_evs_table <- function(norm_counts_mat, ctrl_cols, trt_cols, feature_ids) {
  ctrl_load <- compute_group_pc1_loadings(norm_counts_mat, ctrl_cols)
  trt_load  <- compute_group_pc1_loadings(norm_counts_mat, trt_cols)
  out <- data.frame(
    feature_id = feature_ids,
    ctrl_abs_loading = ctrl_load,
    trt_abs_loading = trt_load,
    stringsAsFactors = FALSE
  )
  out$rank_ctrl <- rank(-out$ctrl_abs_loading, ties.method = "first")
  out$rank_trt  <- rank(-out$trt_abs_loading, ties.method = "first")
  out$combined_rank <- pmin(out$rank_ctrl, out$rank_trt, na.rm = TRUE)
  out <- out[order(out$combined_rank, out$rank_ctrl, out$rank_trt), , drop = FALSE]
  rownames(out) <- NULL
  out
}

compute_feature_metrics <- function(norm_counts_mat, ctrl_cols, trt_cols, dispersions, feature_ids) {
  ctrl_mu  <- apply(norm_counts_mat[, ctrl_cols, drop = FALSE], 1L, safe_mean)
  trt_mu   <- apply(norm_counts_mat[, trt_cols,  drop = FALSE], 1L, safe_mean)
  ctrl_var <- apply(norm_counts_mat[, ctrl_cols, drop = FALSE], 1L, safe_var)
  trt_var  <- apply(norm_counts_mat[, trt_cols,  drop = FALSE], 1L, safe_var)

  out <- data.frame(
    feature_id = feature_ids,
    alpha = as.numeric(dispersions),
    mu_ctrl = ctrl_mu,
    mu_trt = trt_mu,
    var_ctrl = ctrl_var,
    var_trt = trt_var,
    stringsAsFactors = FALSE
  )

  out$iod_emp_ctrl <- ifelse(out$mu_ctrl > 0, out$var_ctrl / out$mu_ctrl, NA_real_)
  out$iod_emp_trt  <- ifelse(out$mu_trt  > 0, out$var_trt  / out$mu_trt,  NA_real_)
  out$cv2_emp_ctrl <- ifelse(out$mu_ctrl > 0, out$var_ctrl / (out$mu_ctrl ^ 2), NA_real_)
  out$cv2_emp_trt  <- ifelse(out$mu_trt  > 0, out$var_trt  / (out$mu_trt  ^ 2), NA_real_)

  out$iod_nb_ctrl <- ifelse(out$mu_ctrl > 0 & is.finite(out$alpha), 1 + out$alpha * out$mu_ctrl, NA_real_)
  out$iod_nb_trt  <- ifelse(out$mu_trt  > 0 & is.finite(out$alpha), 1 + out$alpha * out$mu_trt,  NA_real_)
  out$cv2_nb_ctrl <- ifelse(out$mu_ctrl > 0 & is.finite(out$alpha), (1 / out$mu_ctrl) + out$alpha, NA_real_)
  out$cv2_nb_trt  <- ifelse(out$mu_trt  > 0 & is.finite(out$alpha), (1 / out$mu_trt)  + out$alpha, NA_real_)
  out
}

summarize_by_percentile <- function(df, value_cols, percentile_step = 0.01, min_bin_n = 5L) {
  n <- nrow(df)
  bins <- make_percentile_bins(n, percentile_step)
  rows <- vector("list", nrow(bins))
  for (i in seq_len(nrow(bins))) {
    lo <- if (i == 1L) 1L else bins$rank_index[i - 1L] + 1L
    hi <- bins$rank_index[i]
    chunk <- df[lo:hi, , drop = FALSE]
    row <- data.frame(
      percentile = bins$percentile[i],
      rank_index = bins$rank_index[i],
      lo_rank = lo,
      hi_rank = hi,
      n_bin = nrow(chunk),
      stringsAsFactors = FALSE
    )
    for (nm in value_cols) {
      vals <- chunk[[nm]]
      vals <- vals[is.finite(vals)]
      row[[nm]] <- if (length(vals) >= min_bin_n) mean(vals) else NA_real_
    }
    rows[[i]] <- row
  }
  dplyr::bind_rows(rows)
}

make_group_plots <- function(sum_df, title_prefix, out_file) {
  long_emp <- dplyr::bind_rows(
    data.frame(percentile = sum_df$percentile, value = sum_df$iod_emp_scaled, metric = "IOD empirical"),
    data.frame(percentile = sum_df$percentile, value = sum_df$cv2_emp_scaled, metric = "CV² empirical"),
    data.frame(percentile = sum_df$percentile, value = sum_df$iod_nb_scaled, metric = "IOD NB"),
    data.frame(percentile = sum_df$percentile, value = sum_df$cv2_nb_scaled, metric = "CV² NB")
  )

  diff_df <- data.frame(
    percentile = sum_df$percentile,
    empirical_difference = sum_df$iod_emp_scaled - sum_df$cv2_emp_scaled,
    theoretical_difference = sum_df$iod_nb_scaled - sum_df$cv2_nb_scaled,
    stringsAsFactors = FALSE
  )

  crossing_emp <- find_last_zero_crossing(diff_df$percentile, diff_df$empirical_difference)
  crossing_nb  <- find_last_zero_crossing(diff_df$percentile, diff_df$theoretical_difference)

  p1 <- ggplot(long_emp, aes(percentile, value, color = metric)) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0(title_prefix, ": empirical and NB trajectories"), x = "EVS percentile", y = "Scaled trajectory (0 to 1)") +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")

  p2 <- ggplot(diff_df, aes(percentile, empirical_difference)) +
    geom_hline(yintercept = 0, linetype = 2) + geom_line(linewidth = 0.9) +
    labs(title = paste0(title_prefix, ": empirical difference curve"), x = "EVS percentile", y = "IOD - CV²") +
    theme_bw(base_size = 10)
  if (!is.null(crossing_emp)) p2 <- p2 + geom_vline(xintercept = crossing_emp$crossing_x, linetype = 2)

  p3 <- ggplot(diff_df, aes(percentile, theoretical_difference)) +
    geom_hline(yintercept = 0, linetype = 2) + geom_line(linewidth = 0.9) +
    labs(title = paste0(title_prefix, ": NB-theoretical difference curve"), x = "EVS percentile", y = "IOD - CV²") +
    theme_bw(base_size = 10)
  if (!is.null(crossing_nb)) p3 <- p3 + geom_vline(xintercept = crossing_nb$crossing_x, linetype = 2)

  png(out_file, width = 2000, height = 1600, res = 200)
  gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
  dev.off()

  list(crossing_emp = crossing_emp, crossing_nb = crossing_nb)
}

# =============================================================================
# IMPORT WTTS MASTER COUNT FILE (copied from working script scaffold)
# =============================================================================

if (!file.exists(count_file)) stop(sprintf("Count file not found: %s", count_file), call. = FALSE)

WTTS_Seq <- read.csv(
  count_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
WTTS_Seq <- as.data.frame(WTTS_Seq, stringsAsFactors = FALSE)

# Keep the same required columns as the working scaffold.
assert_required_columns(WTTS_Seq, c("OrigID", "Symbol"), object_name = "WTTS count file")
assert_required_columns(WTTS_Seq, meta_all$id, object_name = "WTTS count file sample columns")

WTTS_Seq$OrigID <- as.character(WTTS_Seq$OrigID)
WTTS_Seq$Symbol <- as.character(WTTS_Seq$Symbol)
WTTS_Seq <- WTTS_Seq[!is.na(WTTS_Seq$OrigID) & !is.na(WTTS_Seq$Symbol), , drop = FALSE]

# Remove rows with missing sample counts exactly like the working scaffold.
sample_na <- rowSums(is.na(WTTS_Seq[, meta_all$id, drop = FALSE])) > 0
WTTS_Seq  <- WTTS_Seq[!sample_na, , drop = FALSE]

# Coerce only the known sample columns to numeric.
for (sid in meta_all$id) {
  WTTS_Seq[[sid]] <- as.numeric(WTTS_Seq[[sid]])
}

bad_after_numeric <- rowSums(!is.finite(as.matrix(WTTS_Seq[, meta_all$id, drop = FALSE]))) > 0
WTTS_Seq <- WTTS_Seq[!bad_after_numeric, , drop = FALSE]

rownames(WTTS_Seq) <- WTTS_Seq$OrigID
feature_ids <- rownames(WTTS_Seq)
count_mat <- as.matrix(WTTS_Seq[, meta_all$id, drop = FALSE])
storage.mode(count_mat) <- "integer"

# =============================================================================
# MAIN LOOP
# =============================================================================

summary_rows <- list()

for (cmp_name in names(comparison_map)) {
  message("Processing ", cmp_name, "...")
  group1_prefix <- comparison_map[[cmp_name]]$group1_prefix
  group2_prefix <- comparison_map[[cmp_name]]$group2_prefix

  keep_ids <- grepl(paste0("^", group1_prefix, "_"), meta_all$id) |
    grepl(paste0("^", group2_prefix, "_"), meta_all$id)

  meta_sub <- meta_all[keep_ids, , drop = FALSE]
  sample_ids <- rownames(meta_sub)
  coldata <- meta_sub[, c("condition"), drop = FALSE]

  missing_samples <- setdiff(sample_ids, colnames(WTTS_Seq))
  if (length(missing_samples) > 0) {
    stop(sprintf("Missing samples in WTTS file for %s: %s", cmp_name, paste(missing_samples, collapse = ", ")), call. = FALSE)
  }

  cmp_counts <- count_mat[, sample_ids, drop = FALSE]
  ctrl_cols <- rownames(meta_sub)[meta_sub$condition == "untrt"]
  trt_cols  <- rownames(meta_sub)[meta_sub$condition == "trt"]

  dds <- DESeqDataSetFromMatrix(
    countData = round(cmp_counts),
    colData = coldata,
    design = ~ condition
  )
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, quiet = TRUE)

  norm_counts <- counts(dds, normalized = TRUE)
  dispersions <- mcols(dds)$dispGeneEst
  if (is.null(dispersions) || all(!is.finite(dispersions))) {
    dispersions <- mcols(dds)$dispersion
  }

  evs_tbl <- build_evs_table(norm_counts, ctrl_cols, trt_cols, rownames(norm_counts))
  metrics_tbl <- compute_feature_metrics(norm_counts, ctrl_cols, trt_cols, dispersions, rownames(norm_counts))
  full_tbl <- dplyr::left_join(evs_tbl, metrics_tbl, by = "feature_id") %>% dplyr::arrange(combined_rank)

  utils::write.csv(full_tbl, file.path(out_dir, paste0(cmp_name, "_feature_metrics.csv")), row.names = FALSE)

  ctrl_rank_tbl <- full_tbl %>% dplyr::arrange(rank_ctrl) %>%
    dplyr::transmute(feature_id, rank = rank_ctrl, iod_emp = iod_emp_ctrl, cv2_emp = cv2_emp_ctrl, iod_nb = iod_nb_ctrl, cv2_nb = cv2_nb_ctrl)
  trt_rank_tbl  <- full_tbl %>% dplyr::arrange(rank_trt) %>%
    dplyr::transmute(feature_id, rank = rank_trt, iod_emp = iod_emp_trt, cv2_emp = cv2_emp_trt, iod_nb = iod_nb_trt, cv2_nb = cv2_nb_trt)

  ctrl_sum <- summarize_by_percentile(ctrl_rank_tbl, c("iod_emp", "cv2_emp", "iod_nb", "cv2_nb"), percentile_step, min_bin_n)
  ctrl_sum$iod_emp_scaled <- rescale01(ctrl_sum$iod_emp)
  ctrl_sum$cv2_emp_scaled <- rescale01(ctrl_sum$cv2_emp)
  ctrl_sum$iod_nb_scaled  <- rescale01(ctrl_sum$iod_nb)
  ctrl_sum$cv2_nb_scaled  <- rescale01(ctrl_sum$cv2_nb)

  trt_sum <- summarize_by_percentile(trt_rank_tbl, c("iod_emp", "cv2_emp", "iod_nb", "cv2_nb"), percentile_step, min_bin_n)
  trt_sum$iod_emp_scaled <- rescale01(trt_sum$iod_emp)
  trt_sum$cv2_emp_scaled <- rescale01(trt_sum$cv2_emp)
  trt_sum$iod_nb_scaled  <- rescale01(trt_sum$iod_nb)
  trt_sum$cv2_nb_scaled  <- rescale01(trt_sum$cv2_nb)

  utils::write.csv(ctrl_sum, file.path(out_dir, paste0(cmp_name, "_control_percentile_summary.csv")), row.names = FALSE)
  utils::write.csv(trt_sum,  file.path(out_dir, paste0(cmp_name, "_treatment_percentile_summary.csv")), row.names = FALSE)

  ctrl_info <- make_group_plots(ctrl_sum, paste0(cmp_name, " control"), file.path(out_dir, paste0(cmp_name, "_control_nb_regime.png")))
  trt_info  <- make_group_plots(trt_sum,  paste0(cmp_name, " treatment"), file.path(out_dir, paste0(cmp_name, "_treatment_nb_regime.png")))

  combined_sum <- dplyr::left_join(
    ctrl_sum[, c("percentile", "rank_index", "iod_emp_scaled", "cv2_emp_scaled", "iod_nb_scaled", "cv2_nb_scaled")],
    trt_sum[, c("percentile", "rank_index", "iod_emp_scaled", "cv2_emp_scaled", "iod_nb_scaled", "cv2_nb_scaled")],
    by = c("percentile", "rank_index"),
    suffix = c("_ctrl", "_trt")
  )

  combined_sum$iod_emp_combined <- combined_sum$iod_emp_scaled_ctrl + combined_sum$iod_emp_scaled_trt
  combined_sum$cv2_emp_combined <- combined_sum$cv2_emp_scaled_ctrl + combined_sum$cv2_emp_scaled_trt
  combined_sum$iod_nb_combined  <- combined_sum$iod_nb_scaled_ctrl + combined_sum$iod_nb_scaled_trt
  combined_sum$cv2_nb_combined  <- combined_sum$cv2_nb_scaled_ctrl + combined_sum$cv2_nb_scaled_trt
  combined_sum$empirical_difference   <- combined_sum$iod_emp_combined - combined_sum$cv2_emp_combined
  combined_sum$theoretical_difference <- combined_sum$iod_nb_combined - combined_sum$cv2_nb_combined

  utils::write.csv(combined_sum, file.path(out_dir, paste0(cmp_name, "_combined_percentile_summary.csv")), row.names = FALSE)

  cross_emp <- find_last_zero_crossing(combined_sum$percentile, combined_sum$empirical_difference)
  cross_nb  <- find_last_zero_crossing(combined_sum$percentile, combined_sum$theoretical_difference)

  crossing_tbl <- data.frame(
    comparison = cmp_name,
    empirical_crossing_percentile = if (is.null(cross_emp)) NA_real_ else cross_emp$crossing_x,
    theoretical_crossing_percentile = if (is.null(cross_nb)) NA_real_ else cross_nb$crossing_x,
    stringsAsFactors = FALSE
  )
  utils::write.csv(crossing_tbl, file.path(out_dir, paste0(cmp_name, "_crossings.csv")), row.names = FALSE)

  p_comb_1 <- ggplot(combined_sum, aes(percentile)) +
    geom_line(aes(y = iod_emp_combined, color = "IOD empirical"), linewidth = 0.9) +
    geom_line(aes(y = cv2_emp_combined, color = "CV² empirical"), linewidth = 0.9) +
    geom_line(aes(y = iod_nb_combined, color = "IOD NB"), linewidth = 0.9, linetype = 2) +
    geom_line(aes(y = cv2_nb_combined, color = "CV² NB"), linewidth = 0.9, linetype = 2) +
    labs(title = paste0(cmp_name, ": combined percentile trajectories"), x = "EVS percentile", y = "Scaled combined trajectory") +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")

  p_comb_2 <- ggplot(combined_sum, aes(percentile, empirical_difference)) +
    geom_hline(yintercept = 0, linetype = 2) + geom_line(linewidth = 0.9) +
    labs(title = paste0(cmp_name, ": empirical combined difference"), x = "EVS percentile", y = "IOD - CV²") +
    theme_bw(base_size = 10)
  if (!is.null(cross_emp)) p_comb_2 <- p_comb_2 + geom_vline(xintercept = cross_emp$crossing_x, linetype = 2)

  p_comb_3 <- ggplot(combined_sum, aes(percentile, theoretical_difference)) +
    geom_hline(yintercept = 0, linetype = 2) + geom_line(linewidth = 0.9) +
    labs(title = paste0(cmp_name, ": NB-theoretical combined difference"), x = "EVS percentile", y = "IOD - CV²") +
    theme_bw(base_size = 10)
  if (!is.null(cross_nb)) p_comb_3 <- p_comb_3 + geom_vline(xintercept = cross_nb$crossing_x, linetype = 2)

  png(file.path(out_dir, paste0(cmp_name, "_combined_nb_regime.png")), width = 2000, height = 1600, res = 200)
  gridExtra::grid.arrange(p_comb_1, p_comb_2, p_comb_3, ncol = 1)
  dev.off()

  summary_rows[[cmp_name]] <- data.frame(
    comparison = cmp_name,
    n_features = nrow(full_tbl),
    empirical_control_crossing = if (is.null(ctrl_info$crossing_emp)) NA_real_ else ctrl_info$crossing_emp$crossing_x,
    empirical_treatment_crossing = if (is.null(trt_info$crossing_emp)) NA_real_ else trt_info$crossing_emp$crossing_x,
    empirical_combined_crossing = if (is.null(cross_emp)) NA_real_ else cross_emp$crossing_x,
    theoretical_combined_crossing = if (is.null(cross_nb)) NA_real_ else cross_nb$crossing_x,
    stringsAsFactors = FALSE
  )
}

summary_tbl <- dplyr::bind_rows(summary_rows)
utils::write.csv(summary_tbl, file.path(out_dir, "nb_regime_summary.csv"), row.names = FALSE)
message("Done. Outputs written to: ", out_dir)
