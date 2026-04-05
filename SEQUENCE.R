# =============================================================================
# SEQUENCE STAGE 1: NEGATIVE-BINOMIAL MEAN-VARIANCE REGIME ANALYSIS
# PRETTY OUTPUT VERSION
# -----------------------------------------------------------------------------
# Keeps the current working analysis logic but improves output structure and
# figure labeling.
#
# Main changes:
#   1. one folder per comparison under exports/nb_regime_analysis_pretty
#   2. cleaner plot themes and legends
#   3. explicit final cutoff labels on the difference panels
#   4. final cutoff rank reported as percentile and approximate EVS rank
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(gridExtra)
})

# =============================================================================
# USER SETTINGS
# =============================================================================

repo_dir <- getwd()
input_dir <- file.path(repo_dir, "data")
count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
out_root <- file.path(repo_dir, "exports", "nb_regime_analysis_pretty")
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

  stop(
    paste0("Count file not found. Tried: ", paste(candidates, collapse = ", ")),
    call. = FALSE
  )
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
  if (!file.exists(path)) stop(paste0("Count file not found: ", path), call. = FALSE)

  raw_df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  feature_col <- detect_feature_id_column(raw_df)
  symbol_col <- detect_gene_symbol_column(raw_df)
  sample_cols <- intersect(meta_ids, names(raw_df))
  if (length(sample_cols) == 0L) stop("No count columns matched the metadata sample IDs.", call. = FALSE)

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
    stop(
      paste0("Missing samples for comparison ", comparison_row$comparison_name, ": ", paste(missing_ids, collapse = ", ")),
      call. = FALSE
    )
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
  dplyr::tibble(percentile = probs[seq_along(centers)], rank_index = centers)
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

  dplyr::tibble(
    feature_id = feature_ids,
    ctrl_abs_loading = ctrl_load,
    trt_abs_loading  = trt_load
  ) %>%
    dplyr::mutate(
      rank_ctrl = rank(-ctrl_abs_loading, ties.method = "first"),
      rank_trt  = rank(-trt_abs_loading, ties.method = "first"),
      combined_rank = pmin(rank_ctrl, rank_trt, na.rm = TRUE)
    ) %>%
    dplyr::arrange(combined_rank, rank_ctrl, rank_trt)
}

compute_feature_metrics <- function(norm_counts_mat, ctrl_cols, trt_cols, dispersions, feature_ids) {
  ctrl_mu  <- apply(norm_counts_mat[, ctrl_cols, drop = FALSE], 1L, safe_mean)
  trt_mu   <- apply(norm_counts_mat[, trt_cols, drop = FALSE], 1L, safe_mean)
  ctrl_var <- apply(norm_counts_mat[, ctrl_cols, drop = FALSE], 1L, safe_var)
  trt_var  <- apply(norm_counts_mat[, trt_cols, drop = FALSE], 1L, safe_var)

  dplyr::tibble(
    feature_id = feature_ids,
    alpha = as.numeric(dispersions),
    mu_ctrl = ctrl_mu,
    mu_trt = trt_mu,
    var_ctrl = ctrl_var,
    var_trt = trt_var
  ) %>%
    dplyr::mutate(
      iod_emp_ctrl = ifelse(mu_ctrl > 0, var_ctrl / mu_ctrl, NA_real_),
      iod_emp_trt  = ifelse(mu_trt  > 0, var_trt / mu_trt, NA_real_),
      cv2_emp_ctrl = ifelse(mu_ctrl > 0, var_ctrl / (mu_ctrl ^ 2), NA_real_),
      cv2_emp_trt  = ifelse(mu_trt  > 0, var_trt / (mu_trt ^ 2), NA_real_),
      iod_nb_ctrl  = ifelse(mu_ctrl > 0 & is.finite(alpha), 1 + alpha * mu_ctrl, NA_real_),
      iod_nb_trt   = ifelse(mu_trt  > 0 & is.finite(alpha), 1 + alpha * mu_trt, NA_real_),
      cv2_nb_ctrl  = ifelse(mu_ctrl > 0 & is.finite(alpha), (1 / mu_ctrl) + alpha, NA_real_),
      cv2_nb_trt   = ifelse(mu_trt  > 0 & is.finite(alpha), (1 / mu_trt) + alpha, NA_real_)
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

    row <- dplyr::tibble(
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
  dplyr::bind_rows(out)
}

pretty_theme <- function() {
  theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

add_cutoff_label <- function(plot_obj, crossing, n_features, color = "black") {
  if (is.null(crossing) || !is.finite(crossing$crossing_x)) return(plot_obj)
  cutoff_rank <- max(1L, min(n_features, round(crossing$crossing_x * n_features)))
  label_txt <- paste0("Final cutoff\nPercentile = ", sprintf("%.3f", crossing$crossing_x),
                      "\nRank = ", format(cutoff_rank, big.mark = ","))
  plot_obj +
    geom_vline(xintercept = crossing$crossing_x, linetype = 2, linewidth = 0.8) +
    annotate("label", x = crossing$crossing_x, y = Inf, vjust = 1.1,
             label = label_txt, size = 3, color = color)
}

make_group_plots <- function(sum_df, title_prefix, out_file, n_features) {
  long_df <- dplyr::bind_rows(
    dplyr::transmute(sum_df, percentile, value = iod_emp_scaled, metric = "IOD empirical"),
    dplyr::transmute(sum_df, percentile, value = cv2_emp_scaled, metric = "CV² empirical"),
    dplyr::transmute(sum_df, percentile, value = iod_nb_scaled, metric = "IOD NB"),
    dplyr::transmute(sum_df, percentile, value = cv2_nb_scaled, metric = "CV² NB")
  )

  diff_df <- dplyr::transmute(
    sum_df,
    percentile,
    empirical_difference = iod_emp_scaled - cv2_emp_scaled,
    theoretical_difference = iod_nb_scaled - cv2_nb_scaled
  )

  crossing_emp <- find_last_zero_crossing(diff_df$percentile, diff_df$empirical_difference)
  crossing_nb <- find_last_zero_crossing(diff_df$percentile, diff_df$theoretical_difference)

  p1 <- ggplot(long_df, aes(percentile, value, color = metric)) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0(title_prefix, ": empirical and NB trajectories"),
         x = "EVS percentile", y = "Scaled trajectory (0 to 1)") +
    pretty_theme()

  p2 <- ggplot(diff_df, aes(percentile, empirical_difference)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0(title_prefix, ": empirical difference curve"),
         x = "EVS percentile", y = "IOD - CV²") +
    pretty_theme()
  p2 <- add_cutoff_label(p2, crossing_emp, n_features)

  p3 <- ggplot(diff_df, aes(percentile, theoretical_difference)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0(title_prefix, ": NB-theoretical difference curve"),
         x = "EVS percentile", y = "IOD - CV²") +
    pretty_theme()
  p3 <- add_cutoff_label(p3, crossing_nb, n_features)

  png(out_file, width = 2200, height = 1800, res = 220)
  gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
  dev.off()

  list(crossing_emp = crossing_emp, crossing_nb = crossing_nb)
}

make_combined_plots <- function(combined_sum, cmp_name, out_file, n_features) {
  cross_emp <- find_last_zero_crossing(combined_sum$percentile, combined_sum$diff_emp)
  cross_nb <- find_last_zero_crossing(combined_sum$percentile, combined_sum$diff_nb)

  p1 <- ggplot(combined_sum, aes(percentile)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6) +
    geom_line(aes(y = diff_emp, color = "Empirical difference"), linewidth = 0.9) +
    geom_line(aes(y = diff_nb, color = "NB difference"), linewidth = 0.9, linetype = 2) +
    labs(title = paste0(cmp_name, ": combined feature-level difference curves"),
         x = "EVS percentile", y = "IOD - CV²") +
    pretty_theme()
  p1 <- add_cutoff_label(p1, cross_emp, n_features)

  png(out_file, width = 2200, height = 1200, res = 220)
  print(p1)
  dev.off()

  list(crossing_emp = cross_emp, crossing_nb = cross_nb)
}

# =============================================================================
# DATA IMPORT
# =============================================================================

count_file <- resolve_counts_file(count_file)
message("Reading raw count matrix from: ", count_file)
loaded <- read_count_matrix(count_file, meta_all$id)
count_mat <- loaded$count_matrix
feature_ids <- rownames(count_mat)

# =============================================================================
# MAIN LOOP
# =============================================================================

summary_rows <- list()

for (i in seq_len(nrow(comparison_table))) {
  comparison_row <- comparison_table[i, , drop = FALSE]
  cmp_name <- comparison_row$comparison_name[[1]]
  message("Processing ", cmp_name, "...")

  cmp_dir <- file.path(out_root, cmp_name)
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
  disp_df <- as.data.frame(SummarizedExperiment::mcols(dds))
  dispersions <- if ("dispGeneEst" %in% names(disp_df)) disp_df$dispGeneEst else disp_df$dispersion

  evs_tbl <- build_evs_table(norm_counts, comp$ctrl_ids, comp$trt_ids, rownames(norm_counts))
  metrics_tbl <- compute_feature_metrics(norm_counts, comp$ctrl_ids, comp$trt_ids, dispersions, rownames(norm_counts))
  full_tbl <- dplyr::left_join(evs_tbl, metrics_tbl, by = "feature_id") %>%
    dplyr::arrange(combined_rank)

  utils::write.csv(full_tbl, file.path(cmp_dir, paste0(cmp_name, "_feature_metrics.csv")), row.names = FALSE)

  ctrl_rank_tbl <- full_tbl %>%
    dplyr::arrange(rank_ctrl) %>%
    dplyr::transmute(feature_id, rank = rank_ctrl, iod_emp = iod_emp_ctrl, cv2_emp = cv2_emp_ctrl, iod_nb = iod_nb_ctrl, cv2_nb = cv2_nb_ctrl)
  trt_rank_tbl <- full_tbl %>%
    dplyr::arrange(rank_trt) %>%
    dplyr::transmute(feature_id, rank = rank_trt, iod_emp = iod_emp_trt, cv2_emp = cv2_emp_trt, iod_nb = iod_nb_trt, cv2_nb = cv2_nb_trt)

  ctrl_sum <- summarize_by_percentile(ctrl_rank_tbl, c("iod_emp", "cv2_emp", "iod_nb", "cv2_nb"), percentile_step, min_bin_n) %>%
    dplyr::mutate(
      iod_emp_scaled = rescale01(iod_emp),
      cv2_emp_scaled = rescale01(cv2_emp),
      iod_nb_scaled = rescale01(iod_nb),
      cv2_nb_scaled = rescale01(cv2_nb)
    )

  trt_sum <- summarize_by_percentile(trt_rank_tbl, c("iod_emp", "cv2_emp", "iod_nb", "cv2_nb"), percentile_step, min_bin_n) %>%
    dplyr::mutate(
      iod_emp_scaled = rescale01(iod_emp),
      cv2_emp_scaled = rescale01(cv2_emp),
      iod_nb_scaled = rescale01(iod_nb),
      cv2_nb_scaled = rescale01(cv2_nb)
    )

  utils::write.csv(ctrl_sum, file.path(cmp_dir, paste0(cmp_name, "_control_percentile_summary.csv")), row.names = FALSE)
  utils::write.csv(trt_sum, file.path(cmp_dir, paste0(cmp_name, "_treatment_percentile_summary.csv")), row.names = FALSE)

  ctrl_info <- make_group_plots(ctrl_sum, paste0(cmp_name, " control"), file.path(cmp_dir, paste0(cmp_name, "_control_nb_regime.png")), nrow(ctrl_rank_tbl))
  trt_info <- make_group_plots(trt_sum, paste0(cmp_name, " treatment"), file.path(cmp_dir, paste0(cmp_name, "_treatment_nb_regime.png")), nrow(trt_rank_tbl))

  combined_rank_tbl <- full_tbl %>%
    dplyr::arrange(combined_rank) %>%
    dplyr::transmute(
      feature_id,
      rank = combined_rank,
      diff_emp = (iod_emp_ctrl - cv2_emp_ctrl) + (iod_emp_trt - cv2_emp_trt),
      diff_nb = (iod_nb_ctrl - cv2_nb_ctrl) + (iod_nb_trt - cv2_nb_trt)
    )

  combined_sum <- summarize_by_percentile(
    combined_rank_tbl,
    c("diff_emp", "diff_nb"),
    percentile_step,
    min_bin_n
  )

  utils::write.csv(combined_rank_tbl, file.path(cmp_dir, paste0(cmp_name, "_combined_feature_level_differences.csv")), row.names = FALSE)
  utils::write.csv(combined_sum, file.path(cmp_dir, paste0(cmp_name, "_combined_percentile_summary.csv")), row.names = FALSE)

  combined_info <- make_combined_plots(combined_sum, cmp_name, file.path(cmp_dir, paste0(cmp_name, "_combined_nb_regime.png")), nrow(combined_rank_tbl))

  crossing_tbl <- dplyr::tibble(
    comparison = cmp_name,
    empirical_control_crossing = if (is.null(ctrl_info$crossing_emp)) NA_real_ else ctrl_info$crossing_emp$crossing_x,
    empirical_treatment_crossing = if (is.null(trt_info$crossing_emp)) NA_real_ else trt_info$crossing_emp$crossing_x,
    empirical_combined_crossing = if (is.null(combined_info$crossing_emp)) NA_real_ else combined_info$crossing_emp$crossing_x,
    theoretical_combined_crossing = if (is.null(combined_info$crossing_nb)) NA_real_ else combined_info$crossing_nb$crossing_x
  )
  utils::write.csv(crossing_tbl, file.path(cmp_dir, paste0(cmp_name, "_crossings.csv")), row.names = FALSE)

  summary_rows[[cmp_name]] <- crossing_tbl
}

summary_tbl <- dplyr::bind_rows(summary_rows)
utils::write.csv(summary_tbl, file.path(out_root, "nb_regime_summary.csv"), row.names = FALSE)
message("Done. Outputs written to: ", out_root)
