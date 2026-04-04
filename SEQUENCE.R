# =============================================================================
# SEQUENCE STAGE 1 REWRITE
# EVS + FULL-AXIS FOURIER + UNION-BASED LEADING EDGE
# =============================================================================
# What this script does
# 1. Reads the raw count matrix.
# 2. Builds normalized counts with DESeq2 for each comparison.
# 3. Builds treatment and control EVS loading tables from PC1 absolute loadings.
# 4. Computes full-axis Fourier fits for IOD and CV2 within treatment and control.
# 5. Combines the treatment and control fitted curves.
# 6. Selects the last crossing before divergence.
# 7. Projects that rank cutoff back to treatment and control EVS ranks.
# 8. Defines the leading edge as the union of treatment and control members with
#    ranks less than or equal to the shared cutoff rank.
# 9. Defines the remainder as everything outside that union.
# 10. Exports only Stage 1 EVS/Fourier tables and figures.
#
# What this script does not do
# - no quantile cutoff logic
# - no stability screen for crossings
# - no higher criticism
# - no HBFSS
# - no shrinkage workflow
# - no volcano plots
# - no final DESeq2 significance calling
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)
  library(scales)
  library(S4Vectors)
})

# =============================================================================
# USER SETTINGS
# =============================================================================

repo_dir <- getwd()
input_dir <- file.path(repo_dir, "data")
output_root <- file.path(repo_dir, "exports")
analysis_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(output_root, paste0("sequence_stage1_evs_fourier_", analysis_stamp))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

auto_count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
count_file <- if (file.exists(auto_count_file)) auto_count_file else auto_count_file

figure_dpi <- 320
base_theme_size <- 10
fourier_harmonics <- 2L
crossing_min_percentile <- 0.05
crossing_max_percentile <- 0.95

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

plot_palette <- list(
  treatment = "#1F78B4",
  control = "#4D4D4D",
  threshold = "#8C2D04",
  histogram = "#969696",
  background = "white"
)

# =============================================================================
# HELPERS
# =============================================================================

assert_required_columns <- function(df, required_cols, object_name = "object") {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      paste0(object_name, " is missing required columns: ", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
}

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88")
    )
}

save_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE)
}

save_plot <- function(plot_obj, path, width = 12, height = 8) {
  ggplot2::ggsave(path, plot = plot_obj, width = width, height = height, dpi = figure_dpi, bg = "white")
}

save_grob <- function(grob_obj, path, width = 14, height = 10) {
  ggplot2::ggsave(path, plot = grob_obj, width = width, height = height, dpi = figure_dpi, bg = "white")
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

read_sequence_count_matrix <- function(path, meta_ids) {
  if (!file.exists(path)) {
    stop(paste0("Count file not found: ", path), call. = FALSE)
  }

  raw_df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  feature_id_col <- detect_feature_id_column(raw_df)
  gene_symbol_col <- detect_gene_symbol_column(raw_df)

  sample_cols <- intersect(meta_ids, names(raw_df))
  if (length(sample_cols) == 0L) {
    stop("No sample columns matching metadata IDs were found in the count file.", call. = FALSE)
  }

  annot_df <- data.frame(
    feature_id = as.character(raw_df[[feature_id_col]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(gene_symbol_col)) {
    annot_df$gene_symbol <- as.character(raw_df[[gene_symbol_col]])
  } else {
    annot_df$gene_symbol <- annot_df$feature_id
  }

  keep <- !is.na(annot_df$feature_id) & nzchar(annot_df$feature_id)
  annot_df <- annot_df[keep, , drop = FALSE]
  count_df <- raw_df[keep, sample_cols, drop = FALSE]

  count_mat <- as.matrix(count_df)
  mode(count_mat) <- "numeric"
  rownames(count_mat) <- annot_df$feature_id
  colnames(count_mat) <- sample_cols

  dup <- duplicated(rownames(count_mat))
  if (any(dup)) {
    count_mat <- rowsum(count_mat, group = rownames(count_mat), reorder = FALSE)
    annot_df <- annot_df[match(rownames(count_mat), annot_df$feature_id), , drop = FALSE]
  }

  list(count_matrix = count_mat, annot_df = annot_df)
}

subset_for_comparison <- function(count_matrix, comparison_row, meta_all) {
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

  coldata <- meta_all[keep_ids, , drop = FALSE]
  count_sub <- count_matrix[, keep_ids, drop = FALSE]
  list(count_matrix = count_sub, coldata = coldata, trt_ids = trt_ids, ctrl_ids = ctrl_ids)
}

compute_group_feature_metrics <- function(count_submatrix) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_submatrix)),
    colData = S4Vectors::DataFrame(row.names = colnames(count_submatrix)),
    design = ~ 1
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds, quiet = TRUE)

  md <- as.data.frame(SummarizedExperiment::mcols(dds), stringsAsFactors = FALSE)
  md$feature_id <- rownames(md)
  keep_cols <- intersect(c("feature_id", "baseMean", "dispGeneEst"), names(md))
  md[, keep_cols, drop = FALSE]
}

compute_normalized_counts <- function(count_matrix, coldata) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_matrix)),
    colData = coldata,
    design = ~ condition
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  counts(dds, normalized = TRUE)
}

compute_pc1_loading_table <- function(norm_counts, sample_ids, feature_metrics, annot_df) {
  x <- log2(as.matrix(norm_counts[, sample_ids, drop = FALSE]) + 1)
  pca_fit <- prcomp(t(x), scale. = FALSE, rank. = 2)
  loadings <- pca_fit$rotation[, 1]

  loading_tbl <- data.frame(
    feature_id = names(loadings),
    pc1_loading = as.numeric(loadings),
    pc1_loading_abs = abs(as.numeric(loadings)),
    stringsAsFactors = FALSE
  )
  loading_tbl <- loading_tbl %>%
    left_join(feature_metrics, by = "feature_id") %>%
    left_join(annot_df, by = "feature_id") %>%
    arrange(desc(pc1_loading_abs), feature_id)
  loading_tbl$rank <- seq_len(nrow(loading_tbl))

  list(pca_fit = pca_fit, loading_table = loading_tbl)
}

build_local_fourier_design <- function(rank_vec, n_harmonics = fourier_harmonics) {
  x <- as.numeric(rank_vec)
  x01 <- (x - min(x)) / max(1e-12, max(x) - min(x))
  design_df <- data.frame(x = x01)
  for (k in seq_len(n_harmonics)) {
    design_df[[paste0("sin_", k)]] <- sin(2 * pi * k * x01)
    design_df[[paste0("cos_", k)]] <- cos(2 * pi * k * x01)
  }
  design_df
}

fit_fourier_curve <- function(rank_vec, y_vec, n_harmonics = fourier_harmonics) {
  keep <- is.finite(rank_vec) & is.finite(y_vec)
  rank_vec <- as.numeric(rank_vec[keep])
  y_vec <- as.numeric(y_vec[keep])

  if (length(rank_vec) < (2L * n_harmonics + 5L)) {
    stop("Not enough finite points to fit the Fourier model.", call. = FALSE)
  }

  design_df <- build_local_fourier_design(rank_vec, n_harmonics = n_harmonics)
  design_df$y <- y_vec
  rhs <- paste(setdiff(names(design_df), "y"), collapse = " + ")
  fit <- stats::lm(stats::as.formula(paste("y ~", rhs)), data = design_df)
  fitted_vals <- as.numeric(stats::predict(fit, newdata = design_df))

  list(
    fit = fit,
    fitted_y = fitted_vals,
    residual_y = y_vec - fitted_vals
  )
}

build_group_fourier_map <- function(loading_tbl) {
  assert_required_columns(loading_tbl, c("feature_id", "rank", "baseMean", "dispGeneEst"), "loading_tbl")

  df <- loading_tbl %>%
    filter(is.finite(baseMean), baseMean > 0, is.finite(dispGeneEst), dispGeneEst > 0) %>%
    arrange(rank)

  if (nrow(df) < 25L) {
    stop("Too few analyzable features after filtering for baseMean and dispGeneEst.", call. = FALSE)
  }

  mu <- pmax(df$baseMean, 1e-12)
  alpha <- pmax(df$dispGeneEst, 1e-12)

  df$iod_nb <- 1 + alpha * mu
  df$cv2_nb <- (1 / mu) + alpha
  df$log_iod_nb <- log10(df$iod_nb)
  df$log_cv2_nb <- log10(df$cv2_nb)

  iod_fit <- fit_fourier_curve(df$rank, df$log_iod_nb)
  cv2_fit <- fit_fourier_curve(df$rank, df$log_cv2_nb)

  wave_map <- data.frame(
    feature_id = df$feature_id,
    rank = df$rank,
    percentile = df$rank / nrow(df),
    iod_fitted = iod_fit$fitted_y,
    cv2_fitted = cv2_fit$fitted_y,
    stringsAsFactors = FALSE
  )

  list(metric_df = df, wave_map = wave_map, iod_fit = iod_fit, cv2_fit = cv2_fit)
}

build_shared_combined_loading_table <- function(trt_loading_tbl, ctrl_loading_tbl) {
  trt_use <- trt_loading_tbl[, c("feature_id", "gene_symbol", "pc1_loading", "pc1_loading_abs", "rank"), drop = FALSE]
  ctrl_use <- ctrl_loading_tbl[, c("feature_id", "gene_symbol", "pc1_loading", "pc1_loading_abs", "rank"), drop = FALSE]

  names(trt_use) <- c("feature_id", "gene_symbol_trt", "pc1_loading_trt", "pc1_loading_abs_trt", "rank_trt")
  names(ctrl_use) <- c("feature_id", "gene_symbol_ctrl", "pc1_loading_ctrl", "pc1_loading_abs_ctrl", "rank_ctrl")

  merged <- full_join(trt_use, ctrl_use, by = "feature_id")
  merged$gene_symbol <- ifelse(
    !is.na(merged$gene_symbol_trt) & nzchar(merged$gene_symbol_trt),
    merged$gene_symbol_trt,
    merged$gene_symbol_ctrl
  )
  merged <- merged %>%
    mutate(
      pc1_loading_abs_trt = ifelse(is.na(pc1_loading_abs_trt), 0, pc1_loading_abs_trt),
      pc1_loading_abs_ctrl = ifelse(is.na(pc1_loading_abs_ctrl), 0, pc1_loading_abs_ctrl),
      combined_center_rank = round(rowMeans(cbind(rank_trt, rank_ctrl), na.rm = TRUE))
    ) %>%
    arrange(combined_center_rank, feature_id)

  merged
}

combine_treatment_control_fourier_maps <- function(trt_wave_obj, ctrl_wave_obj) {
  trt_map <- trt_wave_obj$wave_map[, c("feature_id", "rank", "percentile", "iod_fitted", "cv2_fitted")]
  ctrl_map <- ctrl_wave_obj$wave_map[, c("feature_id", "rank", "percentile", "iod_fitted", "cv2_fitted")]

  names(trt_map) <- c("feature_id", "rank_trt", "percentile_trt", "iod_fitted_trt", "cv2_fitted_trt")
  names(ctrl_map) <- c("feature_id", "rank_ctrl", "percentile_ctrl", "iod_fitted_ctrl", "cv2_fitted_ctrl")

  out <- inner_join(trt_map, ctrl_map, by = "feature_id") %>%
    mutate(
      percentile = rowMeans(cbind(percentile_trt, percentile_ctrl), na.rm = TRUE),
      combined_center_rank = round(rowMeans(cbind(rank_trt, rank_ctrl), na.rm = TRUE)),
      composite_iod = iod_fitted_trt + iod_fitted_ctrl,
      composite_cv2 = cv2_fitted_trt + cv2_fitted_ctrl,
      regime_difference = composite_iod - composite_cv2
    ) %>%
    arrange(percentile, combined_center_rank, feature_id)

  out
}

find_crossings <- function(combined_wave_df) {
  df <- combined_wave_df %>%
    filter(is.finite(percentile), is.finite(regime_difference)) %>%
    filter(percentile >= crossing_min_percentile, percentile <= crossing_max_percentile) %>%
    arrange(percentile)

  if (nrow(df) < 2L) return(data.frame())

  out <- list()
  idx <- 1L
  for (i in seq_len(nrow(df) - 1L)) {
    y1 <- df$regime_difference[i]
    y2 <- df$regime_difference[i + 1L]
    x1 <- df$percentile[i]
    x2 <- df$percentile[i + 1L]
    r1 <- df$combined_center_rank[i]
    r2 <- df$combined_center_rank[i + 1L]

    crossed <- (y1 == 0) || (y2 == 0) || ((y1 > 0) && (y2 < 0)) || ((y1 < 0) && (y2 > 0))
    if (!crossed) next

    if (identical(y1, y2) || isTRUE(all.equal(y1, y2))) {
      crossing_percentile <- mean(c(x1, x2))
      crossing_rank <- round(mean(c(r1, r2)))
    } else {
      crossing_percentile <- x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
      crossing_rank <- round(r1 + (0 - y1) * (r2 - r1) / (y2 - y1))
    }

    out[[idx]] <- data.frame(
      crossing_id = paste0("crossing_", idx),
      percentile_left = x1,
      percentile_right = x2,
      crossing_percentile = crossing_percentile,
      rank_left = r1,
      rank_right = r2,
      crossing_rank = crossing_rank,
      regime_difference_left = y1,
      regime_difference_right = y2,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }

  if (!length(out)) return(data.frame())
  bind_rows(out) %>% arrange(crossing_percentile)
}

select_last_crossing_before_divergence <- function(crossing_tbl) {
  if (is.null(crossing_tbl) || !nrow(crossing_tbl)) {
    stop("No crossings were found within the requested percentile window.", call. = FALSE)
  }
  crossing_tbl %>%
    arrange(desc(crossing_percentile), desc(crossing_rank)) %>%
    slice(1)
}

build_eigenvector_split <- function(shared_combined_tbl, rank_cutoff, count_matrix, norm_counts_df, annot_df) {
  leading_edge_flag <- (
    (!is.na(shared_combined_tbl$rank_trt) & shared_combined_tbl$rank_trt <= rank_cutoff) |
    (!is.na(shared_combined_tbl$rank_ctrl) & shared_combined_tbl$rank_ctrl <= rank_cutoff)
  )

  membership_tbl <- shared_combined_tbl %>%
    mutate(
      in_leading_edge_union = leading_edge_flag,
      in_remainder = !leading_edge_flag
    )

  leading_edge_ids <- as.character(membership_tbl$feature_id[membership_tbl$in_leading_edge_union])
  remainder_ids <- as.character(membership_tbl$feature_id[membership_tbl$in_remainder])

  if (!length(leading_edge_ids)) stop("Leading edge is empty after union-based split.", call. = FALSE)
  if (!length(remainder_ids)) stop("Remainder is empty after union-based split.", call. = FALSE)

  count_feature_ids <- rownames(count_matrix)
  leading_edge_count_matrix <- count_matrix[count_feature_ids %in% leading_edge_ids, , drop = FALSE]
  remainder_count_matrix <- count_matrix[count_feature_ids %in% remainder_ids, , drop = FALSE]

  norm_core <- as.data.frame(norm_counts_df, stringsAsFactors = FALSE)
  norm_core$feature_id <- rownames(norm_core)
  leading_edge_norm_df <- norm_core[norm_core$feature_id %in% leading_edge_ids, , drop = FALSE]
  remainder_norm_df <- norm_core[norm_core$feature_id %in% remainder_ids, , drop = FALSE]

  if (!is.null(annot_df) && nrow(annot_df)) {
    leading_edge_annot_df <- annot_df[annot_df$feature_id %in% leading_edge_ids, , drop = FALSE]
    remainder_annot_df <- annot_df[annot_df$feature_id %in% remainder_ids, , drop = FALSE]
  } else {
    leading_edge_annot_df <- data.frame(feature_id = leading_edge_ids, stringsAsFactors = FALSE)
    remainder_annot_df <- data.frame(feature_id = remainder_ids, stringsAsFactors = FALSE)
  }

  list(
    membership_tbl = membership_tbl,
    leading_edge_ids = leading_edge_ids,
    remainder_ids = remainder_ids,
    leading_edge_count_matrix = leading_edge_count_matrix,
    remainder_count_matrix = remainder_count_matrix,
    leading_edge_norm_df = leading_edge_norm_df,
    remainder_norm_df = remainder_norm_df,
    leading_edge_annot_df = leading_edge_annot_df,
    remainder_annot_df = remainder_annot_df,
    leading_edge_union_n = length(leading_edge_ids),
    remainder_union_n = length(remainder_ids)
  )
}

resolve_stage1_cutoff <- function(trt_loading_tbl, ctrl_loading_tbl, count_matrix, norm_counts_df, annot_df) {
  trt_wave_obj <- build_group_fourier_map(trt_loading_tbl)
  ctrl_wave_obj <- build_group_fourier_map(ctrl_loading_tbl)
  shared_combined_tbl <- build_shared_combined_loading_table(trt_loading_tbl, ctrl_loading_tbl)
  combined_wave_df <- combine_treatment_control_fourier_maps(trt_wave_obj, ctrl_wave_obj)
  crossing_tbl <- find_crossings(combined_wave_df)
  selected_crossing <- select_last_crossing_before_divergence(crossing_tbl)

  rank_cutoff <- as.integer(selected_crossing$crossing_rank[1])
  rank_cutoff <- max(1L, rank_cutoff)

  evs_split <- build_eigenvector_split(
    shared_combined_tbl = shared_combined_tbl,
    rank_cutoff = rank_cutoff,
    count_matrix = count_matrix,
    norm_counts_df = norm_counts_df,
    annot_df = annot_df
  )

  list(
    trt_wave_obj = trt_wave_obj,
    ctrl_wave_obj = ctrl_wave_obj,
    shared_combined_tbl = shared_combined_tbl,
    combined_wave_df = combined_wave_df,
    crossing_tbl = crossing_tbl,
    selected_crossing = selected_crossing,
    rank_cutoff = rank_cutoff,
    evs_split = evs_split,
    leading_edge_ids = evs_split$leading_edge_ids,
    remainder_ids = evs_split$remainder_ids,
    n_total = nrow(shared_combined_tbl),
    leading_edge_union_n = evs_split$leading_edge_union_n,
    remainder_union_n = evs_split$remainder_union_n,
    selected_reason = "last_crossing_before_divergence"
  )
}

# =============================================================================
# PLOTS
# =============================================================================

plot_group_fourier_map <- function(wave_obj, group_label, comparison_name) {
  df <- wave_obj$wave_map
  ggplot(df, aes(percentile)) +
    geom_line(aes(y = iod_fitted, color = "IOD"), linewidth = 0.9) +
    geom_line(aes(y = cv2_fitted, color = "CV²"), linewidth = 0.9) +
    scale_color_manual(values = c("IOD" = plot_palette$treatment, "CV²" = plot_palette$control)) +
    labs(
      title = paste0(comparison_name, " | ", group_label, " full-axis Fourier fit"),
      x = "Rank percentile",
      y = "Fitted log-scale value"
    ) +
    manuscript_theme()
}

plot_composite_overlap <- function(stage1_obj, comparison_name) {
  df <- stage1_obj$combined_wave_df
  sc <- stage1_obj$selected_crossing
  ymax <- max(c(df$composite_iod, df$composite_cv2), na.rm = TRUE)

  ggplot(df, aes(percentile)) +
    geom_line(aes(y = composite_iod, color = "Composite IOD"), linewidth = 0.9) +
    geom_line(aes(y = composite_cv2, color = "Composite CV²"), linewidth = 0.9) +
    geom_vline(xintercept = sc$crossing_percentile[1], linetype = "dashed", linewidth = 0.9, colour = plot_palette$threshold) +
    annotate(
      "label",
      x = sc$crossing_percentile[1],
      y = ymax,
      label = paste0(
        "Last crossing before divergence\n",
        "rank cutoff = ", stage1_obj$rank_cutoff, "\n",
        "leading edge n = ", stage1_obj$leading_edge_union_n, "\n",
        "remainder n = ", stage1_obj$remainder_union_n
      ),
      fill = "white",
      colour = plot_palette$threshold,
      size = 3,
      label.size = 0.15,
      vjust = -0.5
    ) +
    scale_color_manual(values = c("Composite IOD" = plot_palette$treatment, "Composite CV²" = plot_palette$control)) +
    labs(
      title = paste0(comparison_name, " | treatment-control composite overlap"),
      x = "Rank percentile",
      y = "Composite fitted value"
    ) +
    manuscript_theme()
}

plot_difference_curve <- function(stage1_obj, comparison_name) {
  df <- stage1_obj$combined_wave_df
  sc <- stage1_obj$selected_crossing

  ggplot(df, aes(percentile, regime_difference)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.5) +
    geom_line(colour = plot_palette$threshold, linewidth = 0.9) +
    geom_vline(xintercept = sc$crossing_percentile[1], linetype = "dashed", linewidth = 0.9, colour = plot_palette$threshold) +
    geom_point(data = data.frame(percentile = sc$crossing_percentile[1], regime_difference = 0), aes(x = percentile, y = regime_difference), inherit.aes = FALSE, size = 2) +
    annotate(
      "label",
      x = sc$crossing_percentile[1],
      y = 0,
      label = paste0("rank cutoff = ", stage1_obj$rank_cutoff),
      fill = "white",
      colour = plot_palette$threshold,
      size = 3,
      label.size = 0.15,
      vjust = -0.8
    ) +
    labs(
      title = paste0(comparison_name, " | IOD minus CV² crossing curve"),
      x = "Rank percentile",
      y = "IOD − CV²"
    ) +
    manuscript_theme()
}

plot_loading_rank_curve <- function(loading_tbl, rank_cutoff, group_label, comparison_name) {
  ggplot(loading_tbl, aes(rank, pc1_loading_abs)) +
    geom_line(colour = plot_palette$histogram, linewidth = 0.5) +
    geom_vline(xintercept = rank_cutoff, linetype = "dashed", linewidth = 0.9, colour = plot_palette$threshold) +
    labs(
      title = paste0(comparison_name, " | ", group_label, " EVS loading rank curve"),
      x = "EVS rank",
      y = "Absolute PC1 loading"
    ) +
    manuscript_theme()
}

plot_loading_histogram <- function(loading_tbl, group_label, comparison_name) {
  ggplot(loading_tbl, aes(pc1_loading_abs)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, colour = "white") +
    labs(
      title = paste0(comparison_name, " | ", group_label, " EVS loading histogram"),
      x = "Absolute PC1 loading",
      y = "Count"
    ) +
    manuscript_theme()
}

build_summary_panel <- function(stage1_obj, comparison_name) {
  p1 <- plot_group_fourier_map(stage1_obj$trt_wave_obj, "treatment", comparison_name)
  p2 <- plot_group_fourier_map(stage1_obj$ctrl_wave_obj, "control", comparison_name)
  p3 <- plot_composite_overlap(stage1_obj, comparison_name)
  p4 <- plot_difference_curve(stage1_obj, comparison_name)

  gridExtra::arrangeGrob(
    grobs = list(p1, p2, p3, p4),
    ncol = 2,
    top = grid::textGrob(
      paste0(comparison_name, " | treatment, control, overlap, and crossing panels"),
      gp = grid::gpar(fontface = "bold", cex = 1.05)
    )
  )
}

build_loading_panel <- function(trt_loading_tbl, ctrl_loading_tbl, rank_cutoff, comparison_name) {
  p1 <- plot_loading_rank_curve(trt_loading_tbl, rank_cutoff, "treatment", comparison_name)
  p2 <- plot_loading_rank_curve(ctrl_loading_tbl, rank_cutoff, "control", comparison_name)
  p3 <- plot_loading_histogram(trt_loading_tbl, "treatment", comparison_name)
  p4 <- plot_loading_histogram(ctrl_loading_tbl, "control", comparison_name)

  gridExtra::arrangeGrob(
    grobs = list(p1, p2, p3, p4),
    ncol = 2,
    top = grid::textGrob(
      paste0(comparison_name, " | EVS rank and histogram panels"),
      gp = grid::gpar(fontface = "bold", cex = 1.05)
    )
  )
}

# =============================================================================
# TABLE BUILDERS
# =============================================================================

build_stage1_summary_table <- function(stage1_obj, comparison_name) {
  sc <- stage1_obj$selected_crossing
  data.frame(
    comparison_name = comparison_name,
    selected_reason = stage1_obj$selected_reason,
    crossing_id = sc$crossing_id[1],
    crossing_percentile = sc$crossing_percentile[1],
    crossing_rank = sc$crossing_rank[1],
    rank_cutoff = stage1_obj$rank_cutoff,
    n_total = stage1_obj$n_total,
    leading_edge_union_n = stage1_obj$leading_edge_union_n,
    remainder_union_n = stage1_obj$remainder_union_n,
    stringsAsFactors = FALSE
  )
}

build_split_membership_table <- function(stage1_obj) {
  stage1_obj$evs_split$membership_tbl
}

# =============================================================================
# MAIN COMPARISON RUNNER
# =============================================================================

run_one_comparison <- function(comparison_row, count_matrix, annot_df) {
  comparison_name <- comparison_row$comparison_name
  cat("\n--- Running ", comparison_name, " ---\n", sep = "")

  comp <- subset_for_comparison(count_matrix, comparison_row, meta_all)
  norm_counts <- compute_normalized_counts(comp$count_matrix, comp$coldata)

  trt_metrics <- compute_group_feature_metrics(comp$count_matrix[, comp$trt_ids, drop = FALSE])
  ctrl_metrics <- compute_group_feature_metrics(comp$count_matrix[, comp$ctrl_ids, drop = FALSE])

  norm_counts_df <- as.data.frame(norm_counts, stringsAsFactors = FALSE)
  norm_counts_df$feature_id <- rownames(norm_counts_df)
  norm_counts_df <- norm_counts_df %>% left_join(annot_df, by = "feature_id")
  rownames(norm_counts_df) <- norm_counts_df$feature_id

  trt_loading <- compute_pc1_loading_table(norm_counts_df, comp$trt_ids, trt_metrics, annot_df)
  ctrl_loading <- compute_pc1_loading_table(norm_counts_df, comp$ctrl_ids, ctrl_metrics, annot_df)

  stage1_obj <- resolve_stage1_cutoff(trt_loading$loading_table, ctrl_loading$loading_table, comp$count_matrix, norm_counts_df, annot_df)

  comparison_dir <- file.path(output_dir, comparison_name)
  tab_dir <- file.path(comparison_dir, "tables")
  fig_dir <- file.path(comparison_dir, "figures")
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  save_csv(trt_loading$loading_table, file.path(tab_dir, paste0(comparison_name, "_treatment_loading_table.csv")))
  save_csv(ctrl_loading$loading_table, file.path(tab_dir, paste0(comparison_name, "_control_loading_table.csv")))
  save_csv(stage1_obj$trt_wave_obj$wave_map, file.path(tab_dir, paste0(comparison_name, "_treatment_full_axis_fourier_map.csv")))
  save_csv(stage1_obj$ctrl_wave_obj$wave_map, file.path(tab_dir, paste0(comparison_name, "_control_full_axis_fourier_map.csv")))
  save_csv(stage1_obj$combined_wave_df, file.path(tab_dir, paste0(comparison_name, "_combined_full_axis_fourier_map.csv")))
  save_csv(stage1_obj$crossing_tbl, file.path(tab_dir, paste0(comparison_name, "_crossing_table.csv")))
  save_csv(build_stage1_summary_table(stage1_obj, comparison_name), file.path(tab_dir, paste0(comparison_name, "_stage1_summary.csv")))
  save_csv(build_split_membership_table(stage1_obj), file.path(tab_dir, paste0(comparison_name, "_split_membership_table.csv")))
  save_csv(data.frame(feature_id = stage1_obj$leading_edge_ids, stringsAsFactors = FALSE), file.path(tab_dir, paste0(comparison_name, "_leading_edge_ids.csv")))
  save_csv(data.frame(feature_id = stage1_obj$remainder_ids, stringsAsFactors = FALSE), file.path(tab_dir, paste0(comparison_name, "_remainder_ids.csv")))
  save_csv(cbind(data.frame(feature_id = rownames(stage1_obj$evs_split$leading_edge_count_matrix), stringsAsFactors = FALSE), as.data.frame(stage1_obj$evs_split$leading_edge_count_matrix, check.names = FALSE)), file.path(tab_dir, paste0(comparison_name, "_leading_edge_raw_counts.csv")))
  save_csv(cbind(data.frame(feature_id = rownames(stage1_obj$evs_split$remainder_count_matrix), stringsAsFactors = FALSE), as.data.frame(stage1_obj$evs_split$remainder_count_matrix, check.names = FALSE)), file.path(tab_dir, paste0(comparison_name, "_remainder_raw_counts.csv")))
  save_csv(stage1_obj$evs_split$leading_edge_norm_df, file.path(tab_dir, paste0(comparison_name, "_leading_edge_normalized_counts.csv")))
  save_csv(stage1_obj$evs_split$remainder_norm_df, file.path(tab_dir, paste0(comparison_name, "_remainder_normalized_counts.csv")))
  save_csv(stage1_obj$evs_split$leading_edge_annot_df, file.path(tab_dir, paste0(comparison_name, "_leading_edge_annotation.csv")))
  save_csv(stage1_obj$evs_split$remainder_annot_df, file.path(tab_dir, paste0(comparison_name, "_remainder_annotation.csv")))

  save_plot(plot_group_fourier_map(stage1_obj$trt_wave_obj, "treatment", comparison_name), file.path(fig_dir, paste0(comparison_name, "_treatment_full_axis_fourier_fit.png")))
  save_plot(plot_group_fourier_map(stage1_obj$ctrl_wave_obj, "control", comparison_name), file.path(fig_dir, paste0(comparison_name, "_control_full_axis_fourier_fit.png")))
  save_plot(plot_composite_overlap(stage1_obj, comparison_name), file.path(fig_dir, paste0(comparison_name, "_treatment_control_composite_overlap.png")))
  save_plot(plot_difference_curve(stage1_obj, comparison_name), file.path(fig_dir, paste0(comparison_name, "_iod_minus_cv2_crossing_curve.png")))
  save_grob(build_summary_panel(stage1_obj, comparison_name), file.path(fig_dir, paste0(comparison_name, "_fourier_summary_panel.png")))
  save_grob(build_loading_panel(trt_loading$loading_table, ctrl_loading$loading_table, stage1_obj$rank_cutoff, comparison_name), file.path(fig_dir, paste0(comparison_name, "_evs_loading_panel.png")))

  list(
    comparison_name = comparison_name,
    summary = build_stage1_summary_table(stage1_obj, comparison_name)
  )
}

# =============================================================================
# DRIVER
# =============================================================================

main <- function() {
  loaded <- read_sequence_count_matrix(count_file, meta_all$id)
  count_matrix <- loaded$count_matrix
  annot_df <- loaded$annot_df

  results <- lapply(seq_len(nrow(comparison_table)), function(i) {
    run_one_comparison(comparison_table[i, , drop = FALSE], count_matrix, annot_df)
  })

  all_summary <- bind_rows(lapply(results, function(x) x$summary))
  save_csv(all_summary, file.path(output_dir, "all_comparisons_stage1_summary.csv"))
  cat("\nCompleted Stage 1 EVS/Fourier rewrite. Output directory:\n", output_dir, "\n", sep = "")
}

main()
