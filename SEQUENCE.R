# =============================================================================
# SEQUENCE STAGE 1: EMPIRICAL + GENE-WISE NB LEADING-EDGE CUTOFF ANALYSIS
# -----------------------------------------------------------------------------
# Method summary
#   1. Build EVS ranks separately for treatment and control from PC1 absolute
#      loadings using the observed count matrix for the comparison.
#   2. Compute empirical feature-level mean and variance directly from the
#      actual counts (not DESeq2-normalized counts).
#   3. Extract DESeq2 gene-wise dispersion estimates.
#   4. For each feature compute:
#         empirical IOD   = variance / mean
#         empirical CV^2  = variance / mean^2
#         NB IOD          = 1 + alpha_gene_wise * mean
#         NB CV^2         = 1 / mean + alpha_gene_wise
#   5. Summarize feature-level quantities within EVS percentile bins using the
#      median.
#   6. Define the leading-edge cutoff from the leading-edge side:
#         the first meaningful transition point where IOD and CV^2 come
#         together and then begin to diverge with IOD increasing and CV^2
#         decreasing.
#   7. Build one folder per comparison with:
#         - feature-level tables
#         - control percentile summaries
#         - treatment percentile summaries
#         - cutoff table
#         - a labeled 3-panel figure:
#             control trajectories
#             treatment trajectories
#             control+treatment difference curves with cutoff range
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

# =============================================================================
# USER SETTINGS
# =============================================================================

repo_dir <- getwd()
input_dir <- file.path(repo_dir, "data")
count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
out_root <- file.path(repo_dir, "exports", "cutoff_folders_empirical_gene_wise")
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
forward_window_bins <- 3L

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
  raw_df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  feature_col <- detect_feature_id_column(raw_df)
  symbol_col <- detect_gene_symbol_column(raw_df)
  sample_cols <- intersect(meta_ids, names(raw_df))

  if (length(sample_cols) == 0L) {
    stop("No count columns matched metadata sample IDs.", call. = FALSE)
  }

  annot_df <- data.frame(
    feature_id = as.character(raw_df[[feature_col]]),
    stringsAsFactors = FALSE
  )
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

safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

make_percentile_bins <- function(n, step = 0.01) {
  probs <- seq(step, 1, by = step)
  if (tail(probs, 1) < 1) probs <- c(probs, 1)
  centers <- unique(pmax(1L, pmin(n, round(probs * n))))
  tibble(percentile = probs[seq_along(centers)], rank_index = centers)
}

approx_rank_from_percentile <- function(percentile_value, sum_df) {
  ok <- is.finite(sum_df$percentile) & is.finite(sum_df$rank_index)
  if (!any(ok) || !is.finite(percentile_value)) return(NA_real_)
  approx(
    x = sum_df$percentile[ok],
    y = sum_df$rank_index[ok],
    xout = percentile_value,
    rule = 2
  )$y
}

compute_group_pc1_loadings <- function(count_mat, group_cols) {
  mat <- count_mat[, group_cols, drop = FALSE]
  mat <- log2(mat + 1)
  mat <- t(mat)
  if (nrow(mat) < 2L) return(rep(NA_real_, ncol(mat)))
  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  abs(pca$rotation[, 1L])
}

build_evs_table <- function(count_mat, ctrl_cols, trt_cols, feature_ids) {
  ctrl_load <- compute_group_pc1_loadings(count_mat, ctrl_cols)
  trt_load  <- compute_group_pc1_loadings(count_mat, trt_cols)

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

compute_feature_metrics <- function(raw_count_mat, ctrl_cols, trt_cols, alpha_gene_wise, feature_ids) {
  ctrl_mu  <- apply(raw_count_mat[, ctrl_cols, drop = FALSE], 1L, safe_mean)
  trt_mu   <- apply(raw_count_mat[, trt_cols, drop = FALSE], 1L, safe_mean)
  ctrl_var <- apply(raw_count_mat[, ctrl_cols, drop = FALSE], 1L, safe_var)
  trt_var  <- apply(raw_count_mat[, trt_cols, drop = FALSE], 1L, safe_var)

  tibble(
    feature_id = feature_ids,
    alpha_gene_wise = as.numeric(alpha_gene_wise),
    mu_ctrl = ctrl_mu,
    mu_trt = trt_mu,
    var_ctrl = ctrl_var,
    var_trt = trt_var
  ) %>%
    mutate(
      iod_emp_ctrl = ifelse(mu_ctrl > 0, var_ctrl / mu_ctrl, NA_real_),
      iod_emp_trt  = ifelse(mu_trt  > 0, var_trt / mu_trt, NA_real_),
      cv2_emp_ctrl = ifelse(mu_ctrl > 0, var_ctrl / (mu_ctrl ^ 2), NA_real_),
      cv2_emp_trt  = ifelse(mu_trt  > 0, var_trt / (mu_trt ^ 2), NA_real_),

      iod_nb_ctrl  = ifelse(mu_ctrl > 0 & is.finite(alpha_gene_wise), 1 + alpha_gene_wise * mu_ctrl, NA_real_),
      iod_nb_trt   = ifelse(mu_trt  > 0 & is.finite(alpha_gene_wise), 1 + alpha_gene_wise * mu_trt, NA_real_),
      cv2_nb_ctrl  = ifelse(mu_ctrl > 0 & is.finite(alpha_gene_wise), (1 / mu_ctrl) + alpha_gene_wise, NA_real_),
      cv2_nb_trt   = ifelse(mu_trt  > 0 & is.finite(alpha_gene_wise), (1 / mu_trt) + alpha_gene_wise, NA_real_)
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

find_leading_edge_cutoff <- function(sum_df, iod_col, cv2_col, x_col = "percentile", rank_col = "rank_index", forward_bins = 3L) {
  x <- as.numeric(sum_df[[x_col]])
  rank_index <- as.numeric(sum_df[[rank_col]])
  iod <- as.numeric(sum_df[[iod_col]])
  cv2 <- as.numeric(sum_df[[cv2_col]])
  diff <- iod - cv2
  gap <- abs(diff)

  ok <- is.finite(x) & is.finite(rank_index) & is.finite(iod) & is.finite(cv2) & is.finite(diff)
  x <- x[ok]; rank_index <- rank_index[ok]; iod <- iod[ok]; cv2 <- cv2[ok]; diff <- diff[ok]; gap <- gap[ok]

  if (length(x) < 4L) {
    return(list(
      cutoff_percentile = NA_real_,
      cutoff_rank = NA_real_,
      mode = "insufficient_data"
    ))
  }

  # Primary rule: earliest sign change from the leading-edge side followed by
  # widening separation with IOD increasing and CV^2 decreasing.
  for (i in seq_len(length(diff) - 1L)) {
    sign_change <- (diff[i] == 0) || (diff[i + 1L] == 0) || (diff[i] * diff[i + 1L] < 0)
    if (!sign_change) next

    end_idx <- min(length(diff), i + forward_bins + 1L)
    idx <- i:end_idx
    if (length(idx) < 3L) next

    iod_forward <- mean(diff(iod[idx]), na.rm = TRUE)
    cv2_forward <- mean(diff(cv2[idx]), na.rm = TRUE)
    gap_forward <- mean(diff(gap[idx]), na.rm = TRUE)

    if (is.finite(iod_forward) && is.finite(cv2_forward) && is.finite(gap_forward) &&
        iod_forward > 0 && cv2_forward < 0 && gap_forward > 0) {

      x1 <- x[i]; x2 <- x[i + 1L]
      y1 <- diff[i]; y2 <- diff[i + 1L]
      if (isTRUE(all.equal(y1, 0))) {
        cutoff_x <- x1
      } else if (isTRUE(all.equal(y2, 0))) {
        cutoff_x <- x2
      } else {
        cutoff_x <- x1 - y1 * (x2 - x1) / (y2 - y1)
      }
      cutoff_rank <- approx(x = x, y = rank_index, xout = cutoff_x, rule = 2)$y

      return(list(
        cutoff_percentile = cutoff_x,
        cutoff_rank = cutoff_rank,
        mode = "leading_edge_crossing"
      ))
    }
  }

  # Fallback: earliest local minimum in the gap followed by the same divergence.
  for (i in 2:(length(gap) - 1L)) {
    local_min <- is.finite(gap[i - 1L]) && is.finite(gap[i]) && is.finite(gap[i + 1L]) &&
      gap[i] <= gap[i - 1L] && gap[i] <= gap[i + 1L]
    if (!local_min) next

    end_idx <- min(length(gap), i + forward_bins)
    idx <- i:end_idx
    if (length(idx) < 3L) next

    iod_forward <- mean(diff(iod[idx]), na.rm = TRUE)
    cv2_forward <- mean(diff(cv2[idx]), na.rm = TRUE)
    gap_forward <- mean(diff(gap[idx]), na.rm = TRUE)

    if (is.finite(iod_forward) && is.finite(cv2_forward) && is.finite(gap_forward) &&
        iod_forward > 0 && cv2_forward < 0 && gap_forward > 0) {
      return(list(
        cutoff_percentile = x[i],
        cutoff_rank = rank_index[i],
        mode = "leading_edge_closest_approach"
      ))
    }
  }

  # Last fallback: closest approach to balance from the leading-edge side.
  i <- which.min(gap)
  list(
    cutoff_percentile = x[i],
    cutoff_rank = rank_index[i],
    mode = "closest_approach_fallback"
  )
}

build_dataset_plot <- function(sum_df, dataset_label, out_file) {
  cutoff_emp <- find_leading_edge_cutoff(sum_df, "iod_emp", "cv2_emp", forward_bins = forward_window_bins)
  cutoff_nb  <- find_leading_edge_cutoff(sum_df, "iod_nb", "cv2_nb", forward_bins = forward_window_bins)

  curve_long <- bind_rows(
    transmute(sum_df, percentile, value = iod_emp, curve = "IOD empirical"),
    transmute(sum_df, percentile, value = cv2_emp, curve = "CV² empirical"),
    transmute(sum_df, percentile, value = iod_nb, curve = "IOD NB gene-wise"),
    transmute(sum_df, percentile, value = cv2_nb, curve = "CV² NB gene-wise")
  )

  diff_long <- bind_rows(
    transmute(sum_df, percentile, value = iod_emp - cv2_emp, curve = "Empirical difference"),
    transmute(sum_df, percentile, value = iod_nb - cv2_nb, curve = "NB gene-wise difference")
  )

  p1 <- ggplot(curve_long, aes(percentile, value, color = curve, linetype = curve)) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = c(
      "IOD empirical" = "#1b5e20",
      "CV² empirical" = "#b71c1c",
      "IOD NB gene-wise" = "#2e7d32",
      "CV² NB gene-wise" = "#c62828"
    )) +
    scale_linetype_manual(values = c(
      "IOD empirical" = "solid",
      "CV² empirical" = "solid",
      "IOD NB gene-wise" = "dashed",
      "CV² NB gene-wise" = "dashed"
    )) +
    labs(
      title = paste0(dataset_label, ": median bin trajectories"),
      subtitle = "Median of feature-level empirical means/variances from actual counts; gene-wise DESeq2 dispersion for NB curves",
      x = "EVS percentile",
      y = "Median bin value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  if (is.finite(cutoff_emp$cutoff_percentile)) {
    p1 <- p1 + geom_vline(xintercept = cutoff_emp$cutoff_percentile, color = "#0d47a1", linetype = 2, linewidth = 0.9)
    p1 <- p1 + annotate(
      "label",
      x = cutoff_emp$cutoff_percentile,
      y = max(curve_long$value, na.rm = TRUE),
      label = paste0(
        "Empirical cutoff\nmode: ", cutoff_emp$mode,
        "\npercentile: ", sprintf("%.4f", cutoff_emp$cutoff_percentile),
        "\nleading-edge rank: ", sprintf("%.0f", cutoff_emp$cutoff_rank)
      ),
      hjust = 0, vjust = 1, size = 3
    )
  }

  p2 <- ggplot(diff_long, aes(percentile, value, color = curve, linetype = curve)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = c(
      "Empirical difference" = "#1565c0",
      "NB gene-wise difference" = "#6a1b9a"
    )) +
    scale_linetype_manual(values = c(
      "Empirical difference" = "solid",
      "NB gene-wise difference" = "dashed"
    )) +
    labs(
      title = paste0(dataset_label, ": difference curves (IOD - CV²)"),
      subtitle = "Cutoff defined from the leading-edge side where IOD begins increasing and CV² begins decreasing",
      x = "EVS percentile",
      y = "Difference"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  if (is.finite(cutoff_emp$cutoff_percentile)) {
    p2 <- p2 + geom_vline(xintercept = cutoff_emp$cutoff_percentile, color = "#0d47a1", linetype = 2, linewidth = 0.9)
  }
  if (is.finite(cutoff_nb$cutoff_percentile)) {
    p2 <- p2 + geom_vline(xintercept = cutoff_nb$cutoff_percentile, color = "#4a148c", linetype = 3, linewidth = 0.9)
  }

  png(out_file, width = 2200, height = 1800, res = 200)
  grid.arrange(p1, p2, ncol = 1)
  dev.off()

  tibble(
    dataset = dataset_label,
    empirical_cutoff_mode = cutoff_emp$mode,
    empirical_cutoff_percentile = cutoff_emp$cutoff_percentile,
    empirical_cutoff_rank = cutoff_emp$cutoff_rank,
    nb_cutoff_mode = cutoff_nb$mode,
    nb_cutoff_percentile = cutoff_nb$cutoff_percentile,
    nb_cutoff_rank = cutoff_nb$cutoff_rank
  )
}

build_final_range_panel <- function(ctrl_sum, trt_sum, ctrl_cut_tbl, trt_cut_tbl, cmp_name, out_file) {
  ctrl_diff <- ctrl_sum %>% transmute(percentile, value = iod_emp - cv2_emp, curve = "Control empirical difference")
  trt_diff  <- trt_sum  %>% transmute(percentile, value = iod_emp - cv2_emp, curve = "Treatment empirical difference")
  plot_df <- bind_rows(ctrl_diff, trt_diff)

  cut_vals <- c(ctrl_cut_tbl$empirical_cutoff_percentile[1], trt_cut_tbl$empirical_cutoff_percentile[1])
  cut_vals <- cut_vals[is.finite(cut_vals)]
  range_min <- if (length(cut_vals)) min(cut_vals) else NA_real_
  range_max <- if (length(cut_vals)) max(cut_vals) else NA_real_

  p <- ggplot(plot_df, aes(percentile, value, color = curve)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = c(
      "Control empirical difference" = "#00695c",
      "Treatment empirical difference" = "#ef6c00"
    )) +
    labs(
      title = paste0(cmp_name, ": treatment/control cutoff range"),
      subtitle = "Range taken from treatment and control empirical leading-edge cutoffs",
      x = "EVS percentile",
      y = "IOD - CV²"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")

  if (is.finite(range_min) && is.finite(range_max)) {
    p <- p + annotate("rect", xmin = range_min, xmax = range_max, ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "#90caf9")
    p <- p + geom_vline(xintercept = range_min, color = "#0d47a1", linetype = 2, linewidth = 0.9)
    p <- p + geom_vline(xintercept = range_max, color = "#0d47a1", linetype = 2, linewidth = 0.9)
    mid_x <- mean(c(range_min, range_max))
    max_y <- max(plot_df$value, na.rm = TRUE)
    p <- p + annotate(
      "label",
      x = mid_x,
      y = max_y,
      label = paste0(
        "Leading-edge cutoff range\n",
        "percentile: ", sprintf("%.4f", range_min), " to ", sprintf("%.4f", range_max), "\n",
        "control rank: ", sprintf("%.0f", ctrl_cut_tbl$empirical_cutoff_rank[1]), "\n",
        "treatment rank: ", sprintf("%.0f", trt_cut_tbl$empirical_cutoff_rank[1])
      ),
      size = 3.2, vjust = 1, hjust = 0.5
    )
  }

  png(out_file, width = 2200, height = 1200, res = 200)
  print(p)
  dev.off()

  tibble(
    comparison = cmp_name,
    range_percentile_min = range_min,
    range_percentile_max = range_max,
    control_cutoff_rank = ctrl_cut_tbl$empirical_cutoff_rank[1],
    treatment_cutoff_rank = trt_cut_tbl$empirical_cutoff_rank[1]
  )
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

summary_rows <- list()

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

  disp_df <- as.data.frame(mcols(dds))
  alpha_gene_wise <- if ("dispGeneEst" %in% names(disp_df)) disp_df$dispGeneEst else disp_df$dispersion

  evs_tbl <- build_evs_table(as.matrix(cmp_counts), comp$ctrl_ids, comp$trt_ids, rownames(cmp_counts))
  metrics_tbl <- compute_feature_metrics(as.matrix(cmp_counts), comp$ctrl_ids, comp$trt_ids, alpha_gene_wise, rownames(cmp_counts))

  full_tbl <- left_join(evs_tbl, metrics_tbl, by = "feature_id") %>%
    left_join(annot_df, by = c("feature_id" = "feature_id")) %>%
    arrange(combined_rank)

  write.csv(full_tbl, file.path(cmp_dir, paste0(cmp_name, "_feature_level_metrics.csv")), row.names = FALSE)

  ctrl_rank_tbl <- full_tbl %>%
    arrange(rank_ctrl) %>%
    transmute(
      feature_id, gene_symbol, rank = rank_ctrl,
      mu = mu_ctrl,
      variance = var_ctrl,
      iod_emp = iod_emp_ctrl,
      cv2_emp = cv2_emp_ctrl,
      iod_nb = iod_nb_ctrl,
      cv2_nb = cv2_nb_ctrl
    )

  trt_rank_tbl <- full_tbl %>%
    arrange(rank_trt) %>%
    transmute(
      feature_id, gene_symbol, rank = rank_trt,
      mu = mu_trt,
      variance = var_trt,
      iod_emp = iod_emp_trt,
      cv2_emp = cv2_emp_trt,
      iod_nb = iod_nb_trt,
      cv2_nb = cv2_nb_trt
    )

  value_cols <- c("mu", "variance", "iod_emp", "cv2_emp", "iod_nb", "cv2_nb")

  ctrl_sum <- summarize_by_percentile_median(ctrl_rank_tbl, value_cols, percentile_step, min_bin_n)
  trt_sum  <- summarize_by_percentile_median(trt_rank_tbl, value_cols, percentile_step, min_bin_n)

  write.csv(ctrl_sum, file.path(cmp_dir, paste0(cmp_name, "_control_percentile_median_summary.csv")), row.names = FALSE)
  write.csv(trt_sum, file.path(cmp_dir, paste0(cmp_name, "_treatment_percentile_median_summary.csv")), row.names = FALSE)

  ctrl_cut_tbl <- build_dataset_plot(
    ctrl_sum,
    paste0(cmp_name, " control"),
    file.path(cmp_dir, paste0(cmp_name, "_control_cutoff_panel.png"))
  )

  trt_cut_tbl <- build_dataset_plot(
    trt_sum,
    paste0(cmp_name, " treatment"),
    file.path(cmp_dir, paste0(cmp_name, "_treatment_cutoff_panel.png"))
  )

  range_tbl <- build_final_range_panel(
    ctrl_sum,
    trt_sum,
    ctrl_cut_tbl,
    trt_cut_tbl,
    cmp_name,
    file.path(cmp_dir, paste0(cmp_name, "_final_cutoff_range_panel.png"))
  )

  cutoff_tbl <- bind_rows(
    ctrl_cut_tbl %>% mutate(comparison = cmp_name),
    trt_cut_tbl %>% mutate(comparison = cmp_name)
  )
  write.csv(cutoff_tbl, file.path(cmp_dir, paste0(cmp_name, "_cutoff_values.csv")), row.names = FALSE)
  write.csv(range_tbl, file.path(cmp_dir, paste0(cmp_name, "_cutoff_range.csv")), row.names = FALSE)

  summary_rows[[cmp_name]] <- tibble(
    comparison = cmp_name,
    control_cutoff_percentile = ctrl_cut_tbl$empirical_cutoff_percentile[1],
    control_cutoff_rank = ctrl_cut_tbl$empirical_cutoff_rank[1],
    treatment_cutoff_percentile = trt_cut_tbl$empirical_cutoff_percentile[1],
    treatment_cutoff_rank = trt_cut_tbl$empirical_cutoff_rank[1],
    range_percentile_min = range_tbl$range_percentile_min[1],
    range_percentile_max = range_tbl$range_percentile_max[1]
  )
}

summary_tbl <- bind_rows(summary_rows)
write.csv(summary_tbl, file.path(out_root, "overall_cutoff_summary.csv"), row.names = FALSE)
message("Done. Outputs written to: ", out_root)
