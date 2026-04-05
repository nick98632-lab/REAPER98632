# =============================================================================
# SEQUENCE STAGE 1: NB CUTOFF PLAYBOOK
# -----------------------------------------------------------------------------
# Purpose
#   Build a less arbitrary cutoff by combining:
#     1. EVS-ranked features
#     2. feature-level mean and variance from DESeq2-normalized counts
#     3. three DESeq2 dispersion layers:
#          - gene-wise dispersion
#          - fitted/trend dispersion
#          - final/MAP dispersion
#     4. median summaries within EVS percentile bins
#
# Core model
#   Var(X) = mu + alpha * mu^2
#   IOD    = Var(X) / mu   = 1 + alpha * mu
#   CV^2   = Var(X) / mu^2 = 1 / mu + alpha
#
# Strategy
#   For each comparison:
#     - build EVS ranks from treatment and control PC1 absolute loadings
#     - compute feature-level mu, variance, empirical IOD, empirical CV^2
#     - compute NB-based IOD and CV^2 using three dispersion versions
#     - summarize all feature-level quantities within percentile bins using medians
#     - build difference curves: IOD - CV^2
#     - compute cutoff candidates from empirical, gene-wise, fitted, and final curves
#     - use final/MAP as the primary cutoff and the min-max span as a robustness band
#     - export everything into a comparison-specific cutoff folder
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(SummarizedExperiment)
})

# =============================================================================
# USER SETTINGS
# =============================================================================

repo_dir <- getwd()
input_dir <- file.path(repo_dir, "data")
count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
out_root <- file.path(repo_dir, "exports", "cutoff_folders")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  treatment_prefix = c("R0", "R2", "R4", "R8"),
  control_prefix = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

meta_all <- data.frame(
  id = c(
    "R0_1","R0_2","R0_3","R0_4","R0_5","ZT6_1","ZT6_2","ZT6_3","ZT6_4","ZT6_5",
    "R2_1","R2_2","R2_3","R2_4","R2_5","ZT8_1","ZT8_2","ZT8_3","ZT8_4","ZT8_5",
    "R4_1","R4_2","R4_3","R4_4","R4_5","ZT10_1","ZT10_2","ZT10_3","ZT10_4","ZT10_5",
    "R8_1","R8_2","R8_3","R8_4","R8_5","ZT14_1","ZT14_2","ZT14_3","ZT14_4","ZT14_5"
  ),
  condition = c(
    rep("trt", 5), rep("untrt", 5),
    rep("trt", 5), rep("untrt", 5),
    rep("trt", 5), rep("untrt", 5),
    rep("trt", 5), rep("untrt", 5)
  ),
  stringsAsFactors = FALSE
)
rownames(meta_all) <- meta_all$id
meta_all$condition <- factor(meta_all$condition, levels = c("untrt", "trt"))

percentile_step <- 0.01
min_bin_n <- 5L

# =============================================================================
# HELPERS
# =============================================================================

resolve_counts_file <- function(path_hint) {
  candidates <- unique(c(
    path_hint,
    file.path(getwd(), path_hint),
    file.path(getwd(), "data", basename(path_hint)),
    "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
    "/root/REAPER98632/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
  ))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0L) return(existing[[1]])

  found <- list.files(
    path = "/root/REAPER98632",
    pattern = "WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
    recursive = TRUE,
    full.names = TRUE
  )
  found <- found[file.exists(found)]
  if (length(found) > 0L) return(found[[1]])

  stop(paste0("Count file not found. Tried: ", paste(candidates, collapse = ", ")), call. = FALSE)
}

detect_feature_id_column <- function(df) {
  candidates <- c("OrigID", "feature_id", "FeatureID", "PAS", "pas_id", "GeneID", "gene_id", "id")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) return(hit[1])
  names(df)[1]
}

detect_gene_symbol_column <- function(df) {
  candidates <- c("Symbol", "symbol", "GeneSymbol", "gene_symbol", "gene", "Gene")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) return(hit[1])
  NULL
}

read_count_matrix <- function(path, meta_ids) {
  raw_df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  feature_col <- detect_feature_id_column(raw_df)
  symbol_col <- detect_gene_symbol_column(raw_df)
  sample_cols <- intersect(meta_ids, names(raw_df))

  if (length(sample_cols) == 0L) stop("No count columns matched metadata sample IDs.", call. = FALSE)

  annot_df <- data.frame(feature_id = as.character(raw_df[[feature_col]]), stringsAsFactors = FALSE)
  annot_df$gene_symbol <- if (!is.null(symbol_col)) as.character(raw_df[[symbol_col]]) else annot_df$feature_id
  keep <- !is.na(annot_df$feature_id) & nzchar(annot_df$feature_id)
  annot_df <- annot_df[keep, , drop = FALSE]

  count_df <- raw_df[keep, sample_cols, drop = FALSE]
  count_mat <- as.matrix(count_df)
  mode(count_mat) <- "numeric"
  rownames(count_mat) <- annot_df$feature_id
  colnames(count_mat) <- sample_cols

  finite_rows <- rowSums(!is.finite(count_mat)) == 0L
  count_mat <- count_mat[finite_rows, , drop = FALSE]
  annot_df <- annot_df[finite_rows, , drop = FALSE]

  if (anyDuplicated(rownames(count_mat))) {
    count_mat <- rowsum(count_mat, rownames(count_mat), reorder = FALSE)
    annot_df <- annot_df[match(rownames(count_mat), annot_df$feature_id), , drop = FALSE]
  }

  list(count_matrix = count_mat, annot_df = annot_df)
}

subset_comparison <- function(count_matrix, comparison_row, meta_all) {
  trt_ids <- meta_all$id[grepl(paste0("^", comparison_row$treatment_prefix, "_"), meta_all$id)]
  ctrl_ids <- meta_all$id[grepl(paste0("^", comparison_row$control_prefix, "_"), meta_all$id)]
  keep_ids <- c(trt_ids, ctrl_ids)

  missing_ids <- setdiff(keep_ids, colnames(count_matrix))
  if (length(missing_ids) > 0L) {
    stop(paste0("Missing samples for comparison ", comparison_row$comparison_name, ": ", paste(missing_ids, collapse = ", ")), call. = FALSE)
  }

  list(
    count_matrix = count_matrix[, keep_ids, drop = FALSE],
    coldata = meta_all[keep_ids, , drop = FALSE],
    trt_ids = trt_ids,
    ctrl_ids = ctrl_ids
  )
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

safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

rescale01 <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(out)
  rng <- range(x[ok], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || abs(rng[2] - rng[1]) < .Machine$double.eps) {
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
    trt_abs_loading = trt_load
  ) %>%
    mutate(
      rank_ctrl = rank(-ctrl_abs_loading, ties.method = "first"),
      rank_trt = rank(-trt_abs_loading, ties.method = "first"),
      combined_rank = pmin(rank_ctrl, rank_trt, na.rm = TRUE)
    ) %>%
    arrange(combined_rank, rank_ctrl, rank_trt)
}

extract_dispersion_columns <- function(dds) {
  disp_df <- as.data.frame(SummarizedExperiment::mcols(dds))
  alpha_gw <- if ("dispGeneEst" %in% names(disp_df)) disp_df$dispGeneEst else if ("dispersion" %in% names(disp_df)) disp_df$dispersion else rep(NA_real_, nrow(disp_df))
  alpha_fit <- if ("dispFit" %in% names(disp_df)) disp_df$dispFit else rep(NA_real_, nrow(disp_df))
  alpha_final <- if ("dispersion" %in% names(disp_df)) disp_df$dispersion else alpha_gw

  tibble(
    feature_id = rownames(disp_df),
    alpha_gene_wise = as.numeric(alpha_gw),
    alpha_fitted = as.numeric(alpha_fit),
    alpha_final = as.numeric(alpha_final)
  )
}

compute_feature_metrics <- function(norm_counts_mat, ctrl_cols, trt_cols, disp_tbl) {
  ctrl_mu <- apply(norm_counts_mat[, ctrl_cols, drop = FALSE], 1L, safe_mean)
  trt_mu  <- apply(norm_counts_mat[, trt_cols, drop = FALSE], 1L, safe_mean)
  ctrl_var <- apply(norm_counts_mat[, ctrl_cols, drop = FALSE], 1L, safe_var)
  trt_var  <- apply(norm_counts_mat[, trt_cols, drop = FALSE], 1L, safe_var)

  tibble(
    feature_id = rownames(norm_counts_mat),
    mu_ctrl = ctrl_mu,
    mu_trt = trt_mu,
    var_ctrl = ctrl_var,
    var_trt = trt_var
  ) %>%
    left_join(disp_tbl, by = "feature_id") %>%
    mutate(
      iod_emp_ctrl = ifelse(mu_ctrl > 0, var_ctrl / mu_ctrl, NA_real_),
      iod_emp_trt  = ifelse(mu_trt  > 0, var_trt / mu_trt, NA_real_),
      cv2_emp_ctrl = ifelse(mu_ctrl > 0, var_ctrl / (mu_ctrl ^ 2), NA_real_),
      cv2_emp_trt  = ifelse(mu_trt  > 0, var_trt / (mu_trt ^ 2), NA_real_),

      iod_nb_gw_ctrl = ifelse(mu_ctrl > 0 & is.finite(alpha_gene_wise), 1 + alpha_gene_wise * mu_ctrl, NA_real_),
      iod_nb_gw_trt  = ifelse(mu_trt  > 0 & is.finite(alpha_gene_wise), 1 + alpha_gene_wise * mu_trt,  NA_real_),
      cv2_nb_gw_ctrl = ifelse(mu_ctrl > 0 & is.finite(alpha_gene_wise), (1 / mu_ctrl) + alpha_gene_wise, NA_real_),
      cv2_nb_gw_trt  = ifelse(mu_trt  > 0 & is.finite(alpha_gene_wise), (1 / mu_trt)  + alpha_gene_wise, NA_real_),

      iod_nb_fit_ctrl = ifelse(mu_ctrl > 0 & is.finite(alpha_fitted), 1 + alpha_fitted * mu_ctrl, NA_real_),
      iod_nb_fit_trt  = ifelse(mu_trt  > 0 & is.finite(alpha_fitted), 1 + alpha_fitted * mu_trt,  NA_real_),
      cv2_nb_fit_ctrl = ifelse(mu_ctrl > 0 & is.finite(alpha_fitted), (1 / mu_ctrl) + alpha_fitted, NA_real_),
      cv2_nb_fit_trt  = ifelse(mu_trt  > 0 & is.finite(alpha_fitted), (1 / mu_trt)  + alpha_fitted, NA_real_),

      iod_nb_final_ctrl = ifelse(mu_ctrl > 0 & is.finite(alpha_final), 1 + alpha_final * mu_ctrl, NA_real_),
      iod_nb_final_trt  = ifelse(mu_trt  > 0 & is.finite(alpha_final), 1 + alpha_final * mu_trt,  NA_real_),
      cv2_nb_final_ctrl = ifelse(mu_ctrl > 0 & is.finite(alpha_final), (1 / mu_ctrl) + alpha_final, NA_real_),
      cv2_nb_final_trt  = ifelse(mu_trt  > 0 & is.finite(alpha_final), (1 / mu_trt)  + alpha_final, NA_real_)
    )
}

summarize_by_percentile_median <- function(df, value_cols, percentile_step = 0.01, min_bin_n = 5L) {
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
      row[[nm]] <- if (length(vals) >= min_bin_n) stats::median(vals) else NA_real_
    }

    out[[i]] <- row
  }

  bind_rows(out)
}

add_group_differences <- function(sum_df) {
  sum_df %>%
    mutate(
      diff_emp = if ("diff_emp" %in% names(.)) diff_emp else iod_emp - cv2_emp,
      diff_nb_gw = if ("diff_nb_gw" %in% names(.)) diff_nb_gw else iod_nb_gw - cv2_nb_gw,
      diff_nb_fit = if ("diff_nb_fit" %in% names(.)) diff_nb_fit else iod_nb_fit - cv2_nb_fit,
      diff_nb_final = if ("diff_nb_final" %in% names(.)) diff_nb_final else iod_nb_final - cv2_nb_final
    )
}

build_group_plot <- function(sum_df, title_prefix, out_file) {
  curve_long <- bind_rows(
    transmute(sum_df, percentile, value = iod_emp, curve = "IOD empirical"),
    transmute(sum_df, percentile, value = cv2_emp, curve = "CV² empirical"),
    transmute(sum_df, percentile, value = iod_nb_gw, curve = "IOD NB gene-wise"),
    transmute(sum_df, percentile, value = cv2_nb_gw, curve = "CV² NB gene-wise"),
    transmute(sum_df, percentile, value = iod_nb_fit, curve = "IOD NB fitted"),
    transmute(sum_df, percentile, value = cv2_nb_fit, curve = "CV² NB fitted"),
    transmute(sum_df, percentile, value = iod_nb_final, curve = "IOD NB final"),
    transmute(sum_df, percentile, value = cv2_nb_final, curve = "CV² NB final")
  )

  mean_var_long <- bind_rows(
    transmute(sum_df, percentile, value = mu, quantity = "Median mean of normalized counts"),
    transmute(sum_df, percentile, value = variance, quantity = "Median variance of normalized counts")
  )

  diff_long <- bind_rows(
    transmute(sum_df, percentile, value = diff_emp, curve = "Difference empirical"),
    transmute(sum_df, percentile, value = diff_nb_gw, curve = "Difference NB gene-wise"),
    transmute(sum_df, percentile, value = diff_nb_fit, curve = "Difference NB fitted"),
    transmute(sum_df, percentile, value = diff_nb_final, curve = "Difference NB final")
  )

  crossings <- tibble(
    method = c("empirical", "nb_gene_wise", "nb_fitted", "nb_final"),
    crossing_percentile = c(
      if (is.null(find_last_zero_crossing(sum_df$percentile, sum_df$diff_emp))) NA_real_ else find_last_zero_crossing(sum_df$percentile, sum_df$diff_emp)$crossing_x,
      if (is.null(find_last_zero_crossing(sum_df$percentile, sum_df$diff_nb_gw))) NA_real_ else find_last_zero_crossing(sum_df$percentile, sum_df$diff_nb_gw)$crossing_x,
      if (is.null(find_last_zero_crossing(sum_df$percentile, sum_df$diff_nb_fit))) NA_real_ else find_last_zero_crossing(sum_df$percentile, sum_df$diff_nb_fit)$crossing_x,
      if (is.null(find_last_zero_crossing(sum_df$percentile, sum_df$diff_nb_final))) NA_real_ else find_last_zero_crossing(sum_df$percentile, sum_df$diff_nb_final)$crossing_x
    )
  )

  primary_cutoff <- crossings$crossing_percentile[crossings$method == "nb_final"]
  support_range <- range(crossings$crossing_percentile[is.finite(crossings$crossing_percentile)], na.rm = TRUE)
  if (length(support_range) == 0L || !all(is.finite(support_range))) support_range <- c(NA_real_, NA_real_)

  p1 <- ggplot(mean_var_long, aes(percentile, value, color = quantity)) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0(title_prefix, ": median mean and variance by EVS percentile bin"), x = "EVS percentile", y = "Median bin summary") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p2 <- ggplot(curve_long, aes(percentile, value, color = curve)) +
    geom_line(linewidth = 0.8) +
    labs(title = paste0(title_prefix, ": median IOD and CV² trajectories"), x = "EVS percentile", y = "Median bin summary") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p3 <- ggplot(diff_long, aes(percentile, value, color = curve)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0(title_prefix, ": difference curves (IOD - CV²)"), x = "EVS percentile", y = "Difference") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
  if (all(is.finite(support_range))) {
    p3 <- p3 + annotate("rect", xmin = support_range[1], xmax = support_range[2], ymin = -Inf, ymax = Inf, alpha = 0.08)
  }
  if (length(primary_cutoff) == 1L && is.finite(primary_cutoff)) {
    p3 <- p3 + geom_vline(xintercept = primary_cutoff, linetype = 2, linewidth = 0.8)
  }

  png(out_file, width = 2200, height = 2000, res = 200)
  grid.arrange(p1, p2, p3, ncol = 1)
  dev.off()

  list(crossings = crossings, primary_cutoff = primary_cutoff, support_range = support_range)
}

# =============================================================================
# DATA IMPORT
# =============================================================================

count_file <- resolve_counts_file(count_file)
message("Reading raw count matrix from: ", count_file)
loaded <- read_count_matrix(count_file, meta_all$id)
count_mat <- loaded$count_matrix
annot_df <- loaded$annot_df

# =============================================================================
# MAIN LOOP
# =============================================================================

overall_summary_rows <- list()

for (i in seq_len(nrow(comparison_table))) {
  comparison_row <- comparison_table[i, , drop = FALSE]
  cmp_name <- comparison_row$comparison_name[[1]]
  message("Processing ", cmp_name, "...")

  cmp_dir <- file.path(out_root, paste0(cmp_name, "_cutoff_folder"))
  dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

  comp <- subset_comparison(count_mat, comparison_row, meta_all)
  cmp_counts <- comp$count_matrix
  col_data <- comp$coldata

  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(cmp_counts)),
    colData = col_data,
    design = ~ condition
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, quiet = TRUE)

  norm_counts <- counts(dds, normalized = TRUE)
  disp_tbl <- extract_dispersion_columns(dds)
  evs_tbl <- build_evs_table(norm_counts, comp$ctrl_ids, comp$trt_ids, rownames(norm_counts))
  metrics_tbl <- compute_feature_metrics(norm_counts, comp$ctrl_ids, comp$trt_ids, disp_tbl)
  full_tbl <- left_join(evs_tbl, metrics_tbl, by = "feature_id") %>%
    left_join(annot_df, by = c("feature_id" = "feature_id")) %>%
    arrange(combined_rank)

  utils::write.csv(full_tbl, file.path(cmp_dir, paste0(cmp_name, "_feature_level_metrics.csv")), row.names = FALSE)

  ctrl_rank_tbl <- full_tbl %>%
    arrange(rank_ctrl) %>%
    transmute(
      feature_id, gene_symbol, rank = rank_ctrl,
      mu = mu_ctrl,
      variance = var_ctrl,
      iod_emp = iod_emp_ctrl,
      cv2_emp = cv2_emp_ctrl,
      iod_nb_gw = iod_nb_gw_ctrl,
      cv2_nb_gw = cv2_nb_gw_ctrl,
      iod_nb_fit = iod_nb_fit_ctrl,
      cv2_nb_fit = cv2_nb_fit_ctrl,
      iod_nb_final = iod_nb_final_ctrl,
      cv2_nb_final = cv2_nb_final_ctrl,
      alpha_gene_wise, alpha_fitted, alpha_final
    )

  trt_rank_tbl <- full_tbl %>%
    arrange(rank_trt) %>%
    transmute(
      feature_id, gene_symbol, rank = rank_trt,
      mu = mu_trt,
      variance = var_trt,
      iod_emp = iod_emp_trt,
      cv2_emp = cv2_emp_trt,
      iod_nb_gw = iod_nb_gw_trt,
      cv2_nb_gw = cv2_nb_gw_trt,
      iod_nb_fit = iod_nb_fit_trt,
      cv2_nb_fit = cv2_nb_fit_trt,
      iod_nb_final = iod_nb_final_trt,
      cv2_nb_final = cv2_nb_final_trt,
      alpha_gene_wise, alpha_fitted, alpha_final
    )

  ctrl_rank_tbl <- ctrl_rank_tbl %>%
    mutate(
      diff_emp = iod_emp - cv2_emp,
      diff_nb_gw = iod_nb_gw - cv2_nb_gw,
      diff_nb_fit = iod_nb_fit - cv2_nb_fit,
      diff_nb_final = iod_nb_final - cv2_nb_final
    )

  trt_rank_tbl <- trt_rank_tbl %>%
    mutate(
      diff_emp = iod_emp - cv2_emp,
      diff_nb_gw = iod_nb_gw - cv2_nb_gw,
      diff_nb_fit = iod_nb_fit - cv2_nb_fit,
      diff_nb_final = iod_nb_final - cv2_nb_final
    )

  combined_rank_tbl <- full_tbl %>%
    arrange(combined_rank) %>%
    transmute(
      feature_id, gene_symbol, rank = combined_rank,
      mu = mu_ctrl + mu_trt,
      variance = var_ctrl + var_trt,
      iod_emp = iod_emp_ctrl + iod_emp_trt,
      cv2_emp = cv2_emp_ctrl + cv2_emp_trt,
      iod_nb_gw = iod_nb_gw_ctrl + iod_nb_gw_trt,
      cv2_nb_gw = cv2_nb_gw_ctrl + cv2_nb_gw_trt,
      iod_nb_fit = iod_nb_fit_ctrl + iod_nb_fit_trt,
      cv2_nb_fit = cv2_nb_fit_ctrl + cv2_nb_fit_trt,
      iod_nb_final = iod_nb_final_ctrl + iod_nb_final_trt,
      cv2_nb_final = cv2_nb_final_ctrl + cv2_nb_final_trt,
      alpha_gene_wise = alpha_gene_wise,
      alpha_fitted = alpha_fitted,
      alpha_final = alpha_final,
      diff_emp = (iod_emp_ctrl - cv2_emp_ctrl) + (iod_emp_trt - cv2_emp_trt),
      diff_nb_gw = (iod_nb_gw_ctrl - cv2_nb_gw_ctrl) + (iod_nb_gw_trt - cv2_nb_gw_trt),
      diff_nb_fit = (iod_nb_fit_ctrl - cv2_nb_fit_ctrl) + (iod_nb_fit_trt - cv2_nb_fit_trt),
      diff_nb_final = (iod_nb_final_ctrl - cv2_nb_final_ctrl) + (iod_nb_final_trt - cv2_nb_final_trt)
    )

  value_cols <- c(
    "mu", "variance",
    "iod_emp", "cv2_emp",
    "iod_nb_gw", "cv2_nb_gw",
    "iod_nb_fit", "cv2_nb_fit",
    "iod_nb_final", "cv2_nb_final",
    "alpha_gene_wise", "alpha_fitted", "alpha_final",
    "diff_emp", "diff_nb_gw", "diff_nb_fit", "diff_nb_final"
  )

  ctrl_sum <- summarize_by_percentile_median(ctrl_rank_tbl, value_cols, percentile_step, min_bin_n) %>% add_group_differences()
  trt_sum  <- summarize_by_percentile_median(trt_rank_tbl, value_cols, percentile_step, min_bin_n) %>% add_group_differences()
  combined_sum <- summarize_by_percentile_median(combined_rank_tbl, value_cols, percentile_step, min_bin_n) %>% add_group_differences()

  utils::write.csv(ctrl_sum, file.path(cmp_dir, paste0(cmp_name, "_control_percentile_median_summary.csv")), row.names = FALSE)
  utils::write.csv(trt_sum, file.path(cmp_dir, paste0(cmp_name, "_treatment_percentile_median_summary.csv")), row.names = FALSE)
  utils::write.csv(combined_rank_tbl, file.path(cmp_dir, paste0(cmp_name, "_combined_feature_level_metrics.csv")), row.names = FALSE)
  utils::write.csv(combined_sum, file.path(cmp_dir, paste0(cmp_name, "_combined_percentile_median_summary.csv")), row.names = FALSE)

  ctrl_plot_info <- build_group_plot(ctrl_sum, paste0(cmp_name, " control"), file.path(cmp_dir, paste0(cmp_name, "_control_cutoff_playbook.png")))
  trt_plot_info <- build_group_plot(trt_sum, paste0(cmp_name, " treatment"), file.path(cmp_dir, paste0(cmp_name, "_treatment_cutoff_playbook.png")))
  combined_plot_info <- build_group_plot(combined_sum, paste0(cmp_name, " combined"), file.path(cmp_dir, paste0(cmp_name, "_combined_cutoff_playbook.png")))

  crossing_tbl <- bind_rows(
    ctrl_plot_info$crossings %>% mutate(dataset = "control"),
    trt_plot_info$crossings %>% mutate(dataset = "treatment"),
    combined_plot_info$crossings %>% mutate(dataset = "combined")
  ) %>%
    select(dataset, method, crossing_percentile)
  utils::write.csv(crossing_tbl, file.path(cmp_dir, paste0(cmp_name, "_cutoff_candidates.csv")), row.names = FALSE)

  combined_cross <- combined_plot_info$crossings
  primary_cutoff <- combined_cross$crossing_percentile[combined_cross$method == "nb_final"]
  support_vals <- combined_cross$crossing_percentile[is.finite(combined_cross$crossing_percentile)]
  support_min <- if (length(support_vals)) min(support_vals) else NA_real_
  support_max <- if (length(support_vals)) max(support_vals) else NA_real_

  playbook_summary <- tibble(
    comparison = cmp_name,
    n_features = nrow(full_tbl),
    primary_cutoff_method = "nb_final",
    primary_cutoff_percentile = if (length(primary_cutoff) == 1L) primary_cutoff else NA_real_,
    support_range_min = support_min,
    support_range_max = support_max,
    empirical_cutoff = combined_cross$crossing_percentile[combined_cross$method == "empirical"],
    gene_wise_cutoff = combined_cross$crossing_percentile[combined_cross$method == "nb_gene_wise"],
    fitted_cutoff = combined_cross$crossing_percentile[combined_cross$method == "nb_fitted"],
    final_cutoff = combined_cross$crossing_percentile[combined_cross$method == "nb_final"]
  )
  utils::write.csv(playbook_summary, file.path(cmp_dir, paste0(cmp_name, "_playbook_summary.csv")), row.names = FALSE)

  overall_summary_rows[[cmp_name]] <- playbook_summary
}

overall_summary_tbl <- bind_rows(overall_summary_rows)
utils::write.csv(overall_summary_tbl, file.path(out_root, "overall_cutoff_playbook_summary.csv"), row.names = FALSE)
message("Done. Outputs written to: ", out_root)
