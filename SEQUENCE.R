# =============================================================================
# SEQUENCE STAGE 1: NEGATIVE-BINOMIAL MEAN-VARIANCE REGIME ANALYSIS
# -----------------------------------------------------------------------------
# This script replaces the wave-based cutoff logic with a negative-binomial
# interpretation of the EVS-ranked axis.
#
# Core model:
#   Var(X) = mu + alpha * mu^2
#
# Therefore:
#   IOD  = Var(X) / mu   = 1 + alpha * mu
#   CV^2 = Var(X) / mu^2 = 1 / mu + alpha
#
# The goal is to examine whether the EVS-ranked axis reveals a transition in
# the relative dominance of these two NB-derived variance normalizations.
#
# This script:
#   1. reads the raw count matrix
#   2. builds EVS ranks for each comparison from treatment and control PC1
#      absolute loadings
#   3. estimates DESeq2 size factors and feature-wise dispersions
#   4. computes empirical and NB-theoretical IOD and CV^2 for each feature
#   5. summarizes trajectories over percentile bins
#   6. computes the empirical and theoretical difference curves
#   7. identifies the last crossing before divergence
#   8. exports tables and figures
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
  library(gridExtra)
})

# =============================================================================
# USER SETTINGS
# =============================================================================

counts_file <- "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
out_dir <- "exports/nb_regime_analysis"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
counts_file <- resolve_counts_file(counts_file)
message(sprintf('Using count file: %s', counts_file))

comparison_map <- list(
  RT0_ZT6   = list(ctrl = "^R0_", trt = "^ZT6_"),
  RT2_ZT8   = list(ctrl = "^R2_", trt = "^ZT8_"),
  RT4_ZT10  = list(ctrl = "^R4_", trt = "^ZT10_"),
  RT8_ZT14  = list(ctrl = "^R8_", trt = "^ZT14_")
)

percentile_step <- 0.01
min_bin_n <- 5L

# =============================================================================
# HELPERS
# =============================================================================


resolve_counts_file <- function(path_hint) {
  candidates <- c(
    path_hint,
    file.path('data', path_hint),
    file.path('.', path_hint),
    '/root/REAPER98632/WTTS-Seq_2022.2_DE_raw_read_numbers.csv',
    '/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv'
  )
  candidates <- unique(candidates)
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) return(existing[[1]])

  found <- Sys.glob('/root/REAPER98632/**/WTTS-Seq_2022.2_DE_raw_read_numbers.csv')
  found <- found[file.exists(found)]
  if (length(found)) return(found[[1]])

  stop(sprintf('Count file not found. Tried: %s', paste(candidates, collapse=', ')))
}


safe_var <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::var(x)
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
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
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) return(NULL)

  crossings <- integer(0)
  for (i in seq_len(length(y) - 1L)) {
    yi <- y[i]
    yj <- y[i + 1L]
    if (!is.finite(yi) || !is.finite(yj)) next
    if (yi == 0 || yi * yj < 0) crossings <- c(crossings, i)
  }

  if (!length(crossings)) return(NULL)

  i <- max(crossings)
  x1 <- x[i]; x2 <- x[i + 1L]
  y1 <- y[i]; y2 <- y[i + 1L]

  if (isTRUE(all.equal(y1, 0))) {
    x_cross <- x1
  } else if (isTRUE(all.equal(y2, 0))) {
    x_cross <- x2
  } else {
    x_cross <- x1 - y1 * (x2 - x1) / (y2 - y1)
  }

  list(index_left = i, crossing_x = x_cross, y_left = y1, y_right = y2)
}

make_percentile_bins <- function(n, step = 0.01) {
  probs <- seq(step, 1, by = step)
  if (tail(probs, 1) < 1) probs <- c(probs, 1)
  centers <- unique(pmax(1L, pmin(n, round(probs * n))))
  tibble(percentile = probs[seq_along(centers)], rank_index = centers)
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

  tibble(
    feature_id = feature_ids,
    ctrl_abs_loading = ctrl_load,
    trt_abs_loading  = trt_load
  ) %>%
    mutate(
      rank_ctrl = rank(-ctrl_abs_loading, ties.method = "first"),
      rank_trt  = rank(-trt_abs_loading, ties.method = "first"),
      combined_rank = pmin(rank_ctrl, rank_trt, na.rm = TRUE)
    ) %>%
    arrange(combined_rank, rank_ctrl, rank_trt)
}

compute_feature_metrics <- function(norm_counts_mat, ctrl_cols, trt_cols, dispersions, feature_ids) {
  ctrl_mu <- apply(norm_counts_mat[, ctrl_cols, drop = FALSE], 1L, safe_mean)
  trt_mu  <- apply(norm_counts_mat[, trt_cols, drop = FALSE], 1L, safe_mean)
  ctrl_var <- apply(norm_counts_mat[, ctrl_cols, drop = FALSE], 1L, safe_var)
  trt_var  <- apply(norm_counts_mat[, trt_cols, drop = FALSE], 1L, safe_var)

  tibble(
    feature_id = feature_ids,
    alpha = as.numeric(dispersions),
    mu_ctrl = ctrl_mu,
    mu_trt  = trt_mu,
    var_ctrl = ctrl_var,
    var_trt  = trt_var
  ) %>%
    mutate(
      iod_emp_ctrl = ifelse(mu_ctrl > 0, var_ctrl / mu_ctrl, NA_real_),
      iod_emp_trt  = ifelse(mu_trt  > 0, var_trt  / mu_trt,  NA_real_),
      cv2_emp_ctrl = ifelse(mu_ctrl > 0, var_ctrl / (mu_ctrl ^ 2), NA_real_),
      cv2_emp_trt  = ifelse(mu_trt  > 0, var_trt  / (mu_trt  ^ 2), NA_real_),
      iod_nb_ctrl  = ifelse(mu_ctrl > 0 & is.finite(alpha), 1 + alpha * mu_ctrl, NA_real_),
      iod_nb_trt   = ifelse(mu_trt  > 0 & is.finite(alpha), 1 + alpha * mu_trt,  NA_real_),
      cv2_nb_ctrl  = ifelse(mu_ctrl > 0 & is.finite(alpha), (1 / mu_ctrl) + alpha, NA_real_),
      cv2_nb_trt   = ifelse(mu_trt  > 0 & is.finite(alpha), (1 / mu_trt)  + alpha, NA_real_)
    )
}

summarize_by_percentile <- function(df, value_cols, percentile_step = 0.01, min_bin_n = 5L) {
  n <- nrow(df)
  bins <- make_percentile_bins(n, percentile_step)
  out <- vector("list", nrow(bins))

  for (i in seq_len(nrow(bins))) {
    lo <- if (i == 1L) 1L else bins$rank_index[i - 1L] + 1L
    hi <- bins$rank_index[i]
    chunk <- df[lo:hi, , drop = FALSE]

    row <- tibble(
      percentile = bins$percentile[i],
      rank_index = bins$rank_index[i],
      lo_rank = lo,
      hi_rank = hi,
      n_bin = nrow(chunk)
    )

    for (nm in value_cols) {
      vals <- chunk[[nm]]
      vals <- vals[is.finite(vals)]
      row[[nm]] <- if (length(vals) >= min_bin_n) mean(vals) else NA_real_
    }

    out[[i]] <- row
  }

  bind_rows(out)
}

make_group_plots <- function(sum_df, title_prefix, out_file) {
  long_emp <- bind_rows(
    sum_df %>% transmute(percentile, value = iod_emp_scaled, metric = "IOD empirical"),
    sum_df %>% transmute(percentile, value = cv2_emp_scaled, metric = "CV² empirical"),
    sum_df %>% transmute(percentile, value = iod_nb_scaled, metric = "IOD NB"),
    sum_df %>% transmute(percentile, value = cv2_nb_scaled, metric = "CV² NB")
  )

  diff_df <- sum_df %>%
    transmute(
      percentile,
      empirical_difference = iod_emp_scaled - cv2_emp_scaled,
      theoretical_difference = iod_nb_scaled - cv2_nb_scaled
    )

  crossing_emp <- find_last_zero_crossing(diff_df$percentile, diff_df$empirical_difference)
  crossing_nb  <- find_last_zero_crossing(diff_df$percentile, diff_df$theoretical_difference)

  p1 <- ggplot(long_emp, aes(percentile, value, color = metric)) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0(title_prefix, ": empirical and NB trajectories"), x = "EVS percentile", y = "Scaled trajectory (0 to 1)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p2 <- ggplot(diff_df, aes(percentile, empirical_difference)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0(title_prefix, ": empirical difference curve"), x = "EVS percentile", y = "IOD - CV²") +
    theme_bw(base_size = 10)
  if (!is.null(crossing_emp)) p2 <- p2 + geom_vline(xintercept = crossing_emp$crossing_x, linetype = 2)

  p3 <- ggplot(diff_df, aes(percentile, theoretical_difference)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0(title_prefix, ": NB-theoretical difference curve"), x = "EVS percentile", y = "IOD - CV²") +
    theme_bw(base_size = 10)
  if (!is.null(crossing_nb)) p3 <- p3 + geom_vline(xintercept = crossing_nb$crossing_x, linetype = 2)

  png(out_file, width = 2000, height = 1600, res = 200)
  grid.arrange(p1, p2, p3, ncol = 1)
  dev.off()

  list(crossing_emp = crossing_emp, crossing_nb = crossing_nb)
}

# =============================================================================
# DATA IMPORT
# =============================================================================

message("Reading raw count matrix...")
raw_df <- read_csv(counts_file, show_col_types = FALSE)
feature_id_col <- names(raw_df)[1]
feature_ids <- raw_df[[feature_id_col]]
count_cols <- setdiff(names(raw_df), feature_id_col)
count_mat <- as.matrix(raw_df[, count_cols, drop = FALSE])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- feature_ids

# =============================================================================
# MAIN LOOP
# =============================================================================

summary_rows <- list()

for (cmp_name in names(comparison_map)) {
  message("Processing ", cmp_name, "...")

  ctrl_pat <- comparison_map[[cmp_name]]$ctrl
  trt_pat  <- comparison_map[[cmp_name]]$trt

  ctrl_cols <- grep(ctrl_pat, colnames(count_mat), value = TRUE)
  trt_cols  <- grep(trt_pat,  colnames(count_mat), value = TRUE)
  if (!length(ctrl_cols) || !length(trt_cols)) next

  cmp_cols <- c(ctrl_cols, trt_cols)
  cmp_counts <- count_mat[, cmp_cols, drop = FALSE]
  col_data <- data.frame(row.names = cmp_cols, condition = factor(c(rep("ctrl", length(ctrl_cols)), rep("trt", length(trt_cols)))))

  dds <- DESeqDataSetFromMatrix(countData = round(cmp_counts), colData = col_data, design = ~ condition)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, quiet = TRUE)

  norm_counts <- counts(dds, normalized = TRUE)
  dispersions <- mcols(dds)$dispersion

  evs_tbl <- build_evs_table(norm_counts, ctrl_cols, trt_cols, rownames(norm_counts))
  metrics_tbl <- compute_feature_metrics(norm_counts, ctrl_cols, trt_cols, dispersions, rownames(norm_counts))
  full_tbl <- evs_tbl %>% left_join(metrics_tbl, by = "feature_id") %>% arrange(combined_rank)

  write_csv(full_tbl, file.path(out_dir, paste0(cmp_name, "_feature_metrics.csv")))

  ctrl_rank_tbl <- full_tbl %>% arrange(rank_ctrl) %>% transmute(feature_id, rank = rank_ctrl, iod_emp = iod_emp_ctrl, cv2_emp = cv2_emp_ctrl, iod_nb = iod_nb_ctrl, cv2_nb = cv2_nb_ctrl)
  trt_rank_tbl  <- full_tbl %>% arrange(rank_trt)  %>% transmute(feature_id, rank = rank_trt,  iod_emp = iod_emp_trt,  cv2_emp = cv2_emp_trt,  iod_nb = iod_nb_trt,  cv2_nb = cv2_nb_trt)

  ctrl_sum <- summarize_by_percentile(ctrl_rank_tbl, c("iod_emp", "cv2_emp", "iod_nb", "cv2_nb"), percentile_step, min_bin_n) %>%
    mutate(iod_emp_scaled = rescale01(iod_emp), cv2_emp_scaled = rescale01(cv2_emp), iod_nb_scaled = rescale01(iod_nb), cv2_nb_scaled = rescale01(cv2_nb))

  trt_sum <- summarize_by_percentile(trt_rank_tbl, c("iod_emp", "cv2_emp", "iod_nb", "cv2_nb"), percentile_step, min_bin_n) %>%
    mutate(iod_emp_scaled = rescale01(iod_emp), cv2_emp_scaled = rescale01(cv2_emp), iod_nb_scaled = rescale01(iod_nb), cv2_nb_scaled = rescale01(cv2_nb))

  write_csv(ctrl_sum, file.path(out_dir, paste0(cmp_name, "_control_percentile_summary.csv")))
  write_csv(trt_sum,  file.path(out_dir, paste0(cmp_name, "_treatment_percentile_summary.csv")))

  ctrl_info <- make_group_plots(ctrl_sum, paste0(cmp_name, " control"), file.path(out_dir, paste0(cmp_name, "_control_nb_regime.png")))
  trt_info  <- make_group_plots(trt_sum,  paste0(cmp_name, " treatment"), file.path(out_dir, paste0(cmp_name, "_treatment_nb_regime.png")))

  combined_sum <- ctrl_sum %>%
    select(percentile, rank_index, iod_emp_scaled, cv2_emp_scaled, iod_nb_scaled, cv2_nb_scaled) %>%
    rename_with(~ paste0(.x, "_ctrl"), -c(percentile, rank_index)) %>%
    left_join(
      trt_sum %>%
        select(percentile, rank_index, iod_emp_scaled, cv2_emp_scaled, iod_nb_scaled, cv2_nb_scaled) %>%
        rename_with(~ paste0(.x, "_trt"), -c(percentile, rank_index)),
      by = c("percentile", "rank_index")
    ) %>%
    mutate(
      iod_emp_combined = iod_emp_scaled_ctrl + iod_emp_scaled_trt,
      cv2_emp_combined = cv2_emp_scaled_ctrl + cv2_emp_scaled_trt,
      iod_nb_combined  = iod_nb_scaled_ctrl + iod_nb_scaled_trt,
      cv2_nb_combined  = cv2_nb_scaled_ctrl + cv2_nb_scaled_trt,
      empirical_difference = iod_emp_combined - cv2_emp_combined,
      theoretical_difference = iod_nb_combined - cv2_nb_combined
    )

  write_csv(combined_sum, file.path(out_dir, paste0(cmp_name, "_combined_percentile_summary.csv")))

  cross_emp <- find_last_zero_crossing(combined_sum$percentile, combined_sum$empirical_difference)
  cross_nb  <- find_last_zero_crossing(combined_sum$percentile, combined_sum$theoretical_difference)

  crossing_tbl <- tibble(
    comparison = cmp_name,
    empirical_crossing_percentile = if (is.null(cross_emp)) NA_real_ else cross_emp$crossing_x,
    theoretical_crossing_percentile = if (is.null(cross_nb)) NA_real_ else cross_nb$crossing_x
  )
  write_csv(crossing_tbl, file.path(out_dir, paste0(cmp_name, "_crossings.csv")))

  p_comb_1 <- ggplot(combined_sum, aes(percentile)) +
    geom_line(aes(y = iod_emp_combined, color = "IOD empirical"), linewidth = 0.9) +
    geom_line(aes(y = cv2_emp_combined, color = "CV² empirical"), linewidth = 0.9) +
    geom_line(aes(y = iod_nb_combined, color = "IOD NB"), linewidth = 0.9, linetype = 2) +
    geom_line(aes(y = cv2_nb_combined, color = "CV² NB"), linewidth = 0.9, linetype = 2) +
    labs(title = paste0(cmp_name, ": combined percentile trajectories"), x = "EVS percentile", y = "Scaled combined trajectory") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p_comb_2 <- ggplot(combined_sum, aes(percentile, empirical_difference)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0(cmp_name, ": empirical combined difference"), x = "EVS percentile", y = "IOD - CV²") +
    theme_bw(base_size = 10)
  if (!is.null(cross_emp)) p_comb_2 <- p_comb_2 + geom_vline(xintercept = cross_emp$crossing_x, linetype = 2)

  p_comb_3 <- ggplot(combined_sum, aes(percentile, theoretical_difference)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0(cmp_name, ": NB-theoretical combined difference"), x = "EVS percentile", y = "IOD - CV²") +
    theme_bw(base_size = 10)
  if (!is.null(cross_nb)) p_comb_3 <- p_comb_3 + geom_vline(xintercept = cross_nb$crossing_x, linetype = 2)

  png(file.path(out_dir, paste0(cmp_name, "_combined_nb_regime.png")), width = 2000, height = 1600, res = 200)
  grid.arrange(p_comb_1, p_comb_2, p_comb_3, ncol = 1)
  dev.off()

  summary_rows[[cmp_name]] <- tibble(
    comparison = cmp_name,
    n_features = nrow(full_tbl),
    empirical_control_crossing = if (is.null(ctrl_info$crossing_emp)) NA_real_ else ctrl_info$crossing_emp$crossing_x,
    empirical_treatment_crossing = if (is.null(trt_info$crossing_emp)) NA_real_ else trt_info$crossing_emp$crossing_x,
    empirical_combined_crossing = if (is.null(cross_emp)) NA_real_ else cross_emp$crossing_x,
    theoretical_combined_crossing = if (is.null(cross_nb)) NA_real_ else cross_nb$crossing_x
  )
}

summary_tbl <- bind_rows(summary_rows)
write_csv(summary_tbl, file.path(out_dir, "nb_regime_summary.csv"))
message("Done. Outputs written to: ", out_dir)
