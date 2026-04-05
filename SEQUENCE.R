# =============================================================================
# SEQUENCE STAGE 1: EMPIRICAL + GENE-WISE NB LEADING-EDGE CUTOFF ANALYSIS
# -----------------------------------------------------------------------------
# Purpose
#   Build treatment- and control-specific leading-edge cutoffs using:
#     1. EVS ranks from PC1 absolute loadings
#     2. empirical mean and variance from the actual count data
#     3. DESeq2 gene-wise dispersion estimates only
#     4. median summaries within EVS percentile bins
#
# Core quantities
#   empirical IOD   = variance / mean
#   empirical CV^2  = variance / mean^2
#   NB IOD          = 1 + alpha_gene_wise * mean
#   NB CV^2         = 1 / mean + alpha_gene_wise
#
# Cutoff definition
#   For treatment and control separately, the cutoff is defined on the
#   leading-edge side as the first meaningful near-balance point at which
#   IOD and CV^2 come together and then begin to diverge with:
#     - IOD increasing
#     - CV^2 decreasing
#     - the gap |IOD - CV^2| beginning to widen
#
# Outputs
#   For each comparison folder:
#     - feature-level metrics
#     - treatment percentile-median summary
#     - control percentile-median summary
#     - cutoff table with treatment/control/range
#     - a final 3-panel figure:
#         1) treatment trajectories with labeled cutoff
#         2) control trajectories with labeled cutoff
#         3) treatment vs control difference curves with shaded cutoff range
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(tibble)
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
forward_window_bins <- 3L
search_fraction_limit <- 0.80

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
    stop("No count columns matched the metadata sample IDs.", call. = FALSE)
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
      paste0(
        "Missing samples for comparison ",
        comparison_row$comparison_name,
        ": ",
        paste(missing_ids, collapse = ", ")
      ),
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

compute_feature_metrics <- function(raw_count_mat, ctrl_cols, trt_cols, alpha_gene_wise_named, feature_ids) {
  raw_count_mat <- raw_count_mat[feature_ids, , drop = FALSE]
  alpha_tbl <- tibble(
    feature_id = names(alpha_gene_wise_named),
    alpha_gene_wise = as.numeric(alpha_gene_wise_named)
  )

  ctrl_mu  <- apply(raw_count_mat[, ctrl_cols, drop = FALSE], 1L, safe_mean)
  trt_mu   <- apply(raw_count_mat[, trt_cols, drop = FALSE], 1L, safe_mean)
  ctrl_var <- apply(raw_count_mat[, ctrl_cols, drop = FALSE], 1L, safe_var)
  trt_var  <- apply(raw_count_mat[, trt_cols, drop = FALSE], 1L, safe_var)

  tibble(
    feature_id = feature_ids,
    mu_ctrl = ctrl_mu,
    mu_trt = trt_mu,
    var_ctrl = ctrl_var,
    var_trt = trt_var
  ) %>%
    left_join(alpha_tbl, by = "feature_id") %>%
    mutate(
      iod_emp_ctrl = ifelse(mu_ctrl > 0, var_ctrl / mu_ctrl, NA_real_),
      iod_emp_trt  = ifelse(mu_trt  > 0, var_trt / mu_trt, NA_real_),
      cv2_emp_ctrl = ifelse(mu_ctrl > 0, var_ctrl / (mu_ctrl ^ 2), NA_real_),
      cv2_emp_trt  = ifelse(mu_trt  > 0, var_trt / (mu_trt ^ 2), NA_real_),
      iod_nb_ctrl  = ifelse(mu_ctrl > 0 & is.finite(alpha_gene_wise), 1 + alpha_gene_wise * mu_ctrl, NA_real_),
      iod_nb_trt   = ifelse(mu_trt  > 0 & is.finite(alpha_gene_wise), 1 + alpha_gene_wise * mu_trt, NA_real_),
      cv2_nb_ctrl  = ifelse(mu_ctrl > 0 & is.finite(alpha_gene_wise), (1 / mu_ctrl) + alpha_gene_wise, NA_real_),
      cv2_nb_trt   = ifelse(mu_trt  > 0 & is.finite(alpha_gene_wise), (1 / mu_trt) + alpha_gene_wise, NA_real_),
      diff_emp_ctrl = iod_emp_ctrl - cv2_emp_ctrl,
      diff_emp_trt  = iod_emp_trt - cv2_emp_trt,
      diff_nb_ctrl  = iod_nb_ctrl - cv2_nb_ctrl,
      diff_nb_trt   = iod_nb_trt - cv2_nb_trt
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

find_leading_edge_cutoff <- function(sum_df, iod_col, cv2_col, forward_window_bins = 3L, search_fraction_limit = 0.80) {
  keep <- is.finite(sum_df$percentile) & is.finite(sum_df[[iod_col]]) & is.finite(sum_df[[cv2_col]])
  df <- sum_df[keep, c("percentile", "rank_index", iod_col, cv2_col), drop = FALSE]
  names(df) <- c("percentile", "rank_index", "iod", "cv2")
  n <- nrow(df)
  if (n < (forward_window_bins + 2L)) {
    idx <- which.min(abs(df$iod - df$cv2))
    return(list(
      cutoff_percentile = df$percentile[idx],
      cutoff_rank = df$rank_index[idx],
      cutoff_method = "minimum_gap_fallback"
    ))
  }

  max_i <- max(1L, min(n - forward_window_bins, floor(n * search_fraction_limit)))
  df$gap <- abs(df$iod - df$cv2)

  candidates <- vector("list", max_i)
  k <- 0L
  for (i in seq_len(max_i)) {
    iod0 <- df$iod[i]
    cv20 <- df$cv2[i]
    iod_f <- df$iod[min(n, i + forward_window_bins)]
    cv2_f <- df$cv2[min(n, i + forward_window_bins)]
    gap0 <- df$gap[i]
    gap_f <- df$gap[min(n, i + forward_window_bins)]

    if (!all(is.finite(c(iod0, cv20, iod_f, cv2_f, gap0, gap_f)))) next

    if ((iod_f > iod0) && (cv2_f < cv20) && (gap_f > gap0)) {
      k <- k + 1L
      candidates[[k]] <- tibble(
        idx = i,
        percentile = df$percentile[i],
        rank_index = df$rank_index[i],
        gap = gap0,
        forward_gap = gap_f,
        forward_iod = iod_f - iod0,
        forward_cv2 = cv2_f - cv20
      )
    }
  }

  if (k > 0L) {
    cand_df <- bind_rows(candidates[seq_len(k)]) %>%
      arrange(gap, percentile)
    best <- cand_df[1, , drop = FALSE]
    return(list(
      cutoff_percentile = best$percentile[[1]],
      cutoff_rank = best$rank_index[[1]],
      cutoff_method = "leading_edge_divergence"
    ))
  }

  idx <- which.min(df$gap[seq_len(max_i)])
  list(
    cutoff_percentile = df$percentile[idx],
    cutoff_rank = df$rank_index[idx],
    cutoff_method = "minimum_gap_fallback"
  )
}

build_single_group_plot <- function(sum_df, group_label, iod_emp_col, cv2_emp_col, iod_nb_col, cv2_nb_col, cutoff_info) {
  curve_df <- bind_rows(
    transmute(sum_df, percentile, value = .data[[iod_emp_col]], curve = "IOD empirical"),
    transmute(sum_df, percentile, value = .data[[cv2_emp_col]], curve = "CV² empirical"),
    transmute(sum_df, percentile, value = .data[[iod_nb_col]], curve = "IOD NB gene-wise"),
    transmute(sum_df, percentile, value = .data[[cv2_nb_col]], curve = "CV² NB gene-wise")
  )

  label_text <- paste0(
    group_label,
    " cutoff\nmethod = ", cutoff_info$cutoff_method,
    "\npercentile = ", sprintf("%.3f", cutoff_info$cutoff_percentile),
    "\nrank = ", round(cutoff_info$cutoff_rank)
  )

  ggplot(curve_df, aes(percentile, value, color = curve)) +
    geom_line(linewidth = 0.95) +
    geom_vline(xintercept = cutoff_info$cutoff_percentile, linetype = 2, linewidth = 0.9) +
    annotate(
      "label",
      x = cutoff_info$cutoff_percentile,
      y = max(curve_df$value, na.rm = TRUE),
      label = label_text,
      hjust = -0.02,
      vjust = 1,
      size = 3.1,
      label.size = 0.25
    ) +
    labs(
      title = paste0(group_label, ": median trajectories of empirical and gene-wise NB quantities"),
      subtitle = "Bins summarize feature-level values by the median. Empirical mean and variance come from the actual count data.",
      x = "EVS percentile (leading edge to remainder)",
      y = "Median bin value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
}

build_range_panel <- function(trt_sum, ctrl_sum, trt_cutoff, ctrl_cutoff, trt_diff_col, ctrl_diff_col) {
  diff_df <- bind_rows(
    transmute(trt_sum, percentile, value = .data[[trt_diff_col]], group = "Treatment difference (IOD - CV²)"),
    transmute(ctrl_sum, percentile, value = .data[[ctrl_diff_col]], group = "Control difference (IOD - CV²)")
  )

  range_min <- min(c(trt_cutoff$cutoff_percentile, ctrl_cutoff$cutoff_percentile), na.rm = TRUE)
  range_max <- max(c(trt_cutoff$cutoff_percentile, ctrl_cutoff$cutoff_percentile), na.rm = TRUE)

  ggplot(diff_df, aes(percentile, value, color = group)) +
    annotate("rect", xmin = range_min, xmax = range_max, ymin = -Inf, ymax = Inf, alpha = 0.08) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.95) +
    geom_vline(xintercept = trt_cutoff$cutoff_percentile, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = ctrl_cutoff$cutoff_percentile, linetype = 3, linewidth = 0.8) +
    annotate(
      "label",
      x = mean(c(range_min, range_max)),
      y = max(diff_df$value, na.rm = TRUE),
      label = paste0(
        "Cutoff range\n",
        "treatment rank = ", round(trt_cutoff$cutoff_rank),
        "\ncontrol rank = ", round(ctrl_cutoff$cutoff_rank),
        "\nrange = [", round(min(c(trt_cutoff$cutoff_rank, ctrl_cutoff$cutoff_rank))),
        ", ", round(max(c(trt_cutoff$cutoff_rank, ctrl_cutoff$cutoff_rank))), "]"
      ),
      vjust = 1,
      size = 3.1,
      label.size = 0.25
    ) +
    labs(
      title = "Treatment-control cutoff range from the leading-edge side",
      subtitle = "The shaded band spans the treatment and control leading-edge cutoffs. Curves are median summaries of feature-level differences.",
      x = "EVS percentile (leading edge to remainder)",
      y = "Median difference"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
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

  disp_df <- as.data.frame(SummarizedExperiment::mcols(dds))
  alpha_gene_wise <- if ("dispGeneEst" %in% names(disp_df)) disp_df$dispGeneEst else disp_df$dispersion
  names(alpha_gene_wise) <- rownames(disp_df)

  raw_counts_aligned <- as.matrix(counts(dds, normalized = FALSE))
  evs_tbl <- build_evs_table(raw_counts_aligned, comp$ctrl_ids, comp$trt_ids, rownames(raw_counts_aligned))
  metrics_tbl <- compute_feature_metrics(raw_counts_aligned, comp$ctrl_ids, comp$trt_ids, alpha_gene_wise, rownames(raw_counts_aligned))

  full_tbl <- left_join(evs_tbl, metrics_tbl, by = "feature_id") %>%
    left_join(annot_df, by = c("feature_id" = "feature_id")) %>%
    arrange(combined_rank)

  utils::write.csv(full_tbl, file.path(cmp_dir, paste0(cmp_name, "_feature_level_metrics.csv")), row.names = FALSE)

  trt_rank_tbl <- full_tbl %>%
    arrange(rank_trt) %>%
    transmute(
      feature_id,
      gene_symbol,
      rank = rank_trt,
      mu = mu_trt,
      variance = var_trt,
      iod_emp = iod_emp_trt,
      cv2_emp = cv2_emp_trt,
      iod_nb = iod_nb_trt,
      cv2_nb = cv2_nb_trt,
      diff_emp = diff_emp_trt,
      diff_nb = diff_nb_trt,
      alpha_gene_wise = alpha_gene_wise
    )

  ctrl_rank_tbl <- full_tbl %>%
    arrange(rank_ctrl) %>%
    transmute(
      feature_id,
      gene_symbol,
      rank = rank_ctrl,
      mu = mu_ctrl,
      variance = var_ctrl,
      iod_emp = iod_emp_ctrl,
      cv2_emp = cv2_emp_ctrl,
      iod_nb = iod_nb_ctrl,
      cv2_nb = cv2_nb_ctrl,
      diff_emp = diff_emp_ctrl,
      diff_nb = diff_nb_ctrl,
      alpha_gene_wise = alpha_gene_wise
    )

  value_cols <- c("mu", "variance", "iod_emp", "cv2_emp", "iod_nb", "cv2_nb", "diff_emp", "diff_nb", "alpha_gene_wise")

  trt_sum <- summarize_by_percentile_median(trt_rank_tbl, value_cols, percentile_step, min_bin_n)
  ctrl_sum <- summarize_by_percentile_median(ctrl_rank_tbl, value_cols, percentile_step, min_bin_n)

  utils::write.csv(trt_sum, file.path(cmp_dir, paste0(cmp_name, "_treatment_percentile_median_summary.csv")), row.names = FALSE)
  utils::write.csv(ctrl_sum, file.path(cmp_dir, paste0(cmp_name, "_control_percentile_median_summary.csv")), row.names = FALSE)

  trt_cutoff <- find_leading_edge_cutoff(trt_sum, "iod_emp", "cv2_emp", forward_window_bins, search_fraction_limit)
  ctrl_cutoff <- find_leading_edge_cutoff(ctrl_sum, "iod_emp", "cv2_emp", forward_window_bins, search_fraction_limit)

  cutoff_tbl <- tibble(
    comparison = cmp_name,
    treatment_cutoff_percentile = trt_cutoff$cutoff_percentile,
    treatment_cutoff_rank = trt_cutoff$cutoff_rank,
    treatment_cutoff_method = trt_cutoff$cutoff_method,
    control_cutoff_percentile = ctrl_cutoff$cutoff_percentile,
    control_cutoff_rank = ctrl_cutoff$cutoff_rank,
    control_cutoff_method = ctrl_cutoff$cutoff_method,
    cutoff_range_rank_min = min(c(trt_cutoff$cutoff_rank, ctrl_cutoff$cutoff_rank), na.rm = TRUE),
    cutoff_range_rank_max = max(c(trt_cutoff$cutoff_rank, ctrl_cutoff$cutoff_rank), na.rm = TRUE)
  )
  utils::write.csv(cutoff_tbl, file.path(cmp_dir, paste0(cmp_name, "_cutoff_table.csv")), row.names = FALSE)

  p_trt <- build_single_group_plot(trt_sum, paste0(cmp_name, " treatment"), "iod_emp", "cv2_emp", "iod_nb", "cv2_nb", trt_cutoff)
  p_ctrl <- build_single_group_plot(ctrl_sum, paste0(cmp_name, " control"), "iod_emp", "cv2_emp", "iod_nb", "cv2_nb", ctrl_cutoff)
  p_range <- build_range_panel(trt_sum, ctrl_sum, trt_cutoff, ctrl_cutoff, "diff_emp", "diff_emp")

  png(file.path(cmp_dir, paste0(cmp_name, "_final_cutoff_panel.png")), width = 2200, height = 2200, res = 200)
  grid.arrange(p_trt, p_ctrl, p_range, ncol = 1)
  dev.off()

  overall_summary_rows[[cmp_name]] <- tibble(
    comparison = cmp_name,
    n_features = nrow(full_tbl),
    treatment_cutoff_rank = trt_cutoff$cutoff_rank,
    control_cutoff_rank = ctrl_cutoff$cutoff_rank,
    cutoff_range_rank_min = min(c(trt_cutoff$cutoff_rank, ctrl_cutoff$cutoff_rank), na.rm = TRUE),
    cutoff_range_rank_max = max(c(trt_cutoff$cutoff_rank, ctrl_cutoff$cutoff_rank), na.rm = TRUE)
  )
}

overall_summary_tbl <- bind_rows(overall_summary_rows)
utils::write.csv(overall_summary_tbl, file.path(out_root, "overall_cutoff_summary.csv"), row.names = FALSE)
message("Done. Outputs written to: ", out_root)
