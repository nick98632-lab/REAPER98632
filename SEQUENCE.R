
#!/usr/bin/env Rscript

# =============================================================================
# SEQUENCE.R
# Final manuscript pipeline for EVS + DESeq2 + empirical-null HC + HBFSS
#
# PURPOSE
# This script implements the manuscript analysis pipeline for four circadian
# pairwise comparisons:
#   RT0_ZT6, RT2_ZT8, RT4_ZT10, RT8_ZT14
#
# The unit of analysis is the PAS / OrigID feature, not the gene. OrigID and
# Symbol are therefore carried separately throughout the pipeline.
#
# CORE METHODS
# 1. Import the WTTS count matrix and embedded sample metadata.
# 2. For each comparison, perform eigenvector splitting (EVS):
#    - estimate DESeq2 size factors on the full comparison matrix
#    - compute condition-specific PCA on DESeq2-normalized counts
#    - rank features by absolute PC1 loading within each condition
#    - define the leading edge as the union of the top N features in either arm
#    - define the remainder as all features not in that union
# 3. Run DESeq2 independently on Raw, Lead, and Rem datasets.
# 4. Recalibrate the DESeq2 Wald statistic distribution with fdrtool's
#    empirical-null normal model.
# 5. Compute the dataset-specific higher-criticism (HC) threshold on empirical
#    p-values, and derive the dataset-specific HBFSS boundary:
#        HBFSS_threshold = abs(log10(HC_p_threshold)) * lfc_boundary
# 6. Compute manuscript plotting classes using HC-gated precedence:
#        Background < Weak CNH < Strong CNH < Standard < HBFSS
#    No feature below the HC threshold may be colored as Weak CNH, Strong CNH,
#    Standard, or HBFSS.
# 7. Export manuscript figures and audit tables with short GitHub-friendly names.
#
# VOLCANO CLASS DEFINITIONS
# Weak CNH:
#   resLA_padj < alpha_level
#   |lfc_shrunk| < lfc_boundary
#   empirical_p <= hc_p_threshold
#
# Strong CNH:
#   resGA_padj < alpha_level
#   |lfc_shrunk| >= lfc_boundary
#   empirical_p <= hc_p_threshold
#
# Standard:
#   padj < alpha_level
#   |lfc_shrunk| >= lfc_boundary
#   empirical_p <= hc_p_threshold
#
# HBFSS:
#   HBFSS >= hbfss_threshold
#   empirical_p <= hc_p_threshold
#
# Final plotting precedence:
#   Background -> Weak CNH -> Strong CNH -> Standard -> HBFSS
#
# OUTPUT DESIGN
# - one clean volcano per dataset
# - one shared-legend compare volcano panel per comparison
# - one shared-legend compare dispersion panel per comparison
# - EVS support plots for PCA and PC1-loading rank
# - short filenames beginning with Figure_ or Table_
# - a PlotLogic table for each dataset to audit class assignments feature-wise
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
})

options(stringsAsFactors = FALSE)

# =============================================================================
# USER SETTINGS
# =============================================================================
count_file <- file.path("data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
output_dir <- file.path("exports", "manuscript_final_clean")

run_twas_overlap <- FALSE
twas_file <- "3aTWAS_genes_of_11_brain_disorders.csv"

alpha_level <- 0.10
lfc_boundary <- 1.0
top_n_target <- 5000L

figure_dpi <- 320
base_theme_size <- 10
n_top_labels_standard <- 12
n_top_labels_hbfss <- 12

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PALETTE, LABELS, SHORT NAMES
# =============================================================================
plot_palette <- list(
  background   = "#BDBDBD",
  weak_cnh     = "#2C7FB8",
  strong_cnh   = "#0B6E4F",
  standard     = "#33A02C",
  hbfss        = "#6A3D9A",
  threshold    = "#A65628",
  histogram    = "#969696",
  control      = "#4D4D4D",
  treatment    = "#1F78B4"
)

dataset_short <- c(
  raw_dataset = "Raw",
  leading_edge_dataset = "Lead",
  remainder_dataset = "Rem"
)

plot_class_levels <- c("Background", "Weak CNH", "Strong CNH", "Standard", "HBFSS")
plot_class_colors <- c(
  "Background" = plot_palette$background,
  "Weak CNH"   = plot_palette$weak_cnh,
  "Strong CNH" = plot_palette$strong_cnh,
  "Standard"   = plot_palette$standard,
  "HBFSS"      = plot_palette$hbfss
)
plot_class_shapes <- c(
  "Background" = 16,
  "Weak CNH"   = 16,
  "Strong CNH" = 17,
  "Standard"   = 15,
  "HBFSS"      = 18
)

fig_file <- function(dir, cmp, ds_key, tag) {
  file.path(dir, paste0("Figure_", cmp, "_", unname(dataset_short[ds_key]), "_", tag, ".png"))
}
tab_file <- function(dir, cmp, ds_key, tag) {
  file.path(dir, paste0("Table_", cmp, "_", unname(dataset_short[ds_key]), "_", tag, ".csv"))
}

# =============================================================================
# EMBEDDED SAMPLE METADATA
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
  )
)
rownames(meta_all) <- meta_all$id
meta_all$condition <- factor(meta_all$condition, levels = c("control", "treatment"))
levels(meta_all$condition) <- c("untrt", "trt")

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix   = c("R0", "R2", "R4", "R8"),
  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14")
)

# =============================================================================
# HELPERS
# =============================================================================
assert_required_columns <- function(df, required_cols, object_name = "data frame") {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste0("Missing required columns in ", object_name, ": ", paste(missing_cols, collapse = ", ")))
  }
}

safe_log10 <- function(x, pseudocount = 1e-12) {
  log10(pmax(x, pseudocount))
}

safe_neglog10 <- function(x, pseudocount = 1e-12) {
  -log10(pmax(x, pseudocount))
}

clip_probabilities <- function(x, eps = 1e-300) {
  x <- unname(as.numeric(x))
  if (!length(x)) return(numeric(0))
  bad <- !is.finite(x) | is.na(x)
  x[bad] <- NA_real_
  good <- !is.na(x)
  x[good] <- pmin(pmax(x[good], eps), 1 - 1e-12)
  x
}

compact_title <- function(x, width = 56) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

compact_caption <- function(x, width = 110) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

save_csv <- function(df, path) {
  write.csv(df, file = path, row.names = FALSE)
}

save_plot <- function(p, path, width = 7, height = 6, dpi = figure_dpi, bg = "white") {
  ggsave(filename = path, plot = p, width = width, height = height,
         dpi = dpi, units = "in", bg = bg, limitsize = FALSE)
}

save_grob <- function(g, path, width = 13, height = 6, dpi = figure_dpi, bg = "white") {
  ggsave(filename = path, plot = g, width = width, height = height,
         dpi = dpi, units = "in", bg = bg, limitsize = FALSE)
}

extract_legend_grob <- function(p) {
  gt <- ggplotGrob(p)
  idx <- which(vapply(gt$grobs, function(x) x$name, character(1)) == "guide-box")
  if (!length(idx)) return(NULL)
  gt$grobs[[idx[1]]]
}

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_theme_size + 1, hjust = 0.5),
      plot.subtitle = element_text(size = base_theme_size - 2, hjust = 0.5),
      plot.caption = element_text(size = base_theme_size - 1, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      plot.margin = margin(12, 16, 12, 16)
    )
}

condition_shapes <- c("untrt" = 21, "trt" = 24)
condition_fills  <- c("untrt" = plot_palette$control, "trt" = plot_palette$treatment)
condition_labels <- c("untrt" = "Control", "trt" = "Treatment")

pretty_group_label <- function(group_label) {
  switch(group_label,
    trt = "Treatment",
    untrt = "Control",
    treatment = "Treatment",
    control = "Control",
    group_label
  )
}

make_design_formula <- function(coldata) ~ condition

get_condition_coef <- function(dds) {
  rn <- resultsNames(dds)
  idx <- grep("^condition_", rn)
  if (!length(idx)) stop("Could not identify condition coefficient.")
  rn[idx[1]]
}

resolve_top_n_cutoff <- function(sorted_values_desc, top_n = top_n_target) {
  n_total <- length(sorted_values_desc)
  if (n_total == 0) stop("resolve_top_n_cutoff() received an empty vector.")
  top_n_actual <- min(max(1L, as.integer(top_n)), n_total)
  cutoff_value <- sorted_values_desc[top_n_actual]
  cutoff_quantile <- 1 - (top_n_actual / n_total)
  list(top_n_actual = top_n_actual, cutoff_value = cutoff_value,
       cutoff_quantile = cutoff_quantile, n_total = n_total)
}

run_empirical_null_fdrtool <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  if (length(stat_vec) < 5) {
    stop(sprintf("[%s] fewer than 5 finite Wald statistics available for fdrtool.", dataset_name))
  }

  fit <- tryCatch(
    fdrtool(stat_vec, statistic = "normal", plot = FALSE, verbose = FALSE,
            cutoff.method = "fndr", pct0 = 0.75),
    error = function(e1) {
      message(sprintf("[%s] primary fdrtool call failed: %s", dataset_name, conditionMessage(e1)))
      fdrtool(as.vector(stat_vec), statistic = "normal", plot = FALSE, verbose = FALSE,
              cutoff.method = "pct0", pct0 = 0.75)
    }
  )

  fit$pval <- clip_probabilities(fit$pval)
  fit$qval <- clip_probabilities(fit$qval)
  fit$lfdr <- as.numeric(fit$lfdr)
  fit
}

safe_hc_thresh <- function(empirical_p, dataset_name) {
  sorted_empirical_p <- sort(clip_probabilities(empirical_p), na.last = NA, decreasing = FALSE)
  if (length(sorted_empirical_p) < 5) return(NA_real_)
  out <- suppressWarnings(
    tryCatch(fdrtool::hc.thresh(as.vector(sorted_empirical_p)),
             error = function(e) {
               message(sprintf("[%s] hc.thresh failed: %s", dataset_name, conditionMessage(e)))
               NA_real_
             })
  )
  out <- as.numeric(out[1])
  if (!is.finite(out) || is.na(out) || out <= 0 || out >= 1) return(NA_real_)
  out
}

classify_effect_strength <- function(res_strong_padj, res_weak_padj, alpha = alpha_level) {
  out <- rep("intermediate", length(res_strong_padj))
  out[!is.na(res_weak_padj) & res_weak_padj < alpha] <- "weak_effect"
  out[!is.na(res_strong_padj) & res_strong_padj < alpha] <- "strong_effect"
  out
}

make_hbfss_curve_df <- function(xmax_abs, hbfss_threshold, ymax_plot, x_min_abs = 0.08, n = 800) {
  if (!is.finite(hbfss_threshold) || is.na(hbfss_threshold) || hbfss_threshold <= 0) return(NULL)
  x_left <- seq(-xmax_abs, -x_min_abs, length.out = n %/% 2)
  x_right <- seq(x_min_abs, xmax_abs, length.out = n %/% 2)
  x <- c(x_left, x_right)
  y <- hbfss_threshold / abs(x)
  keep <- is.finite(y) & !is.na(y) & y <= ymax_plot
  data.frame(x = x[keep], y = y[keep])
}

# =============================================================================
# IMPORT COUNT MATRIX
# =============================================================================
WTTS_Seq <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
WTTS_Seq <- as.data.frame(WTTS_Seq, stringsAsFactors = FALSE)
assert_required_columns(WTTS_Seq, c("OrigID", "Symbol"), "WTTS count file")
assert_required_columns(WTTS_Seq, meta_all$id, "WTTS count file sample columns")

WTTS_Seq$OrigID <- as.character(WTTS_Seq$OrigID)
WTTS_Seq$Symbol <- as.character(WTTS_Seq$Symbol)
WTTS_Seq <- WTTS_Seq[!is.na(WTTS_Seq$OrigID) & !is.na(WTTS_Seq$Symbol), , drop = FALSE]
sample_na <- rowSums(is.na(WTTS_Seq[, meta_all$id, drop = FALSE])) > 0
WTTS_Seq <- WTTS_Seq[!sample_na, , drop = FALSE]
rownames(WTTS_Seq) <- WTTS_Seq$OrigID

OrigID_Symbol <- unique(WTTS_Seq[, c("OrigID", "Symbol"), drop = FALSE])
colnames(OrigID_Symbol) <- c("feature_id", "gene_symbol")
OrigID_Symbol <- OrigID_Symbol %>%
  dplyr::mutate(
    feature_id = as.character(feature_id),
    gene_symbol = dplyr::if_else(is.na(gene_symbol), "", trimws(as.character(gene_symbol)))
  ) %>%
  dplyr::arrange(feature_id, dplyr::desc(gene_symbol != ""), gene_symbol) %>%
  dplyr::distinct(feature_id, .keep_all = TRUE) %>%
  dplyr::mutate(gene_symbol = dplyr::na_if(gene_symbol, ""))

if (run_twas_overlap) {
  TWAS_Seq <- read.csv(twas_file, header = TRUE, stringsAsFactors = FALSE)
  TWAS_data <- TWAS_Seq[, c(1, 4), drop = FALSE]
  colnames(TWAS_data) <- c("source_id", "gene_symbol")
}

# =============================================================================
# SECTION 2 - EVS PREPARATION
# =============================================================================
prepare_comparison_data <- function(comparison_name, group1_prefix, group2_prefix, WTTS_Seq, meta_all) {
  keep_ids <- grepl(paste0("^", group1_prefix, "_"), meta_all$id) |
    grepl(paste0("^", group2_prefix, "_"), meta_all$id)
  meta_sub <- meta_all[keep_ids, , drop = FALSE]
  coldata <- meta_sub[, "condition", drop = FALSE]
  sample_ids <- rownames(meta_sub)

  missing_samples <- setdiff(sample_ids, colnames(WTTS_Seq))
  if (length(missing_samples) > 0) {
    stop(paste("Missing samples in WTTS file for", comparison_name, ":",
               paste(missing_samples, collapse = ", ")))
  }

  count_sub <- WTTS_Seq[, sample_ids, drop = FALSE]
  count_mat <- as.matrix(count_sub)
  storage.mode(count_mat) <- "numeric"
  rownames(count_mat) <- WTTS_Seq$OrigID
  stopifnot(all(colnames(count_mat) == rownames(coldata)))

  list(comparison_name = comparison_name, count_matrix = count_mat, coldata = coldata)
}

compute_pc1_loading_table <- function(value_df, sample_names, top_n = top_n_target,
                                      preprocessing_label = "Normalized prior to eigenvector splitting") {
  x <- as.matrix(value_df[, sample_names, drop = FALSE])
  pca_fit <- prcomp(t(x), scale. = FALSE, rank. = 2)
  loading_abs <- abs(pca_fit$rotation[, 1])
  loading_tbl <- data.frame(feature_id = names(loading_abs),
                            pc1_loading_abs = unname(loading_abs))
  loading_tbl <- loading_tbl[order(loading_tbl$pc1_loading_abs, decreasing = TRUE), , drop = FALSE]
  loading_tbl$rank <- seq_len(nrow(loading_tbl))

  cutoff_info <- resolve_top_n_cutoff(loading_tbl$pc1_loading_abs, top_n = top_n)
  cutoff <- cutoff_info$cutoff_value

  loading_tbl$split_class <- ifelse(
    loading_tbl$pc1_loading_abs >= cutoff,
    "high_loading",
    "background_loading"
  )

  list(
    pca_fit = pca_fit,
    loading_table = loading_tbl,
    cutoff = cutoff,
    top_n_used = cutoff_info$top_n_actual,
    cutoff_quantile = cutoff_info$cutoff_quantile,
    preprocessing_label = preprocessing_label
  )
}

build_eigenvector_split <- function(count_matrix, coldata) {
  design_formula <- make_design_formula(coldata)

  dds_init <- DESeqDataSetFromMatrix(countData = round(count_matrix), colData = coldata, design = design_formula)
  dds_init <- dds_init[rowSums(counts(dds_init)) > 0, ]
  dds_init <- estimateSizeFactors(dds_init)

  norm_counts_init <- as.data.frame(counts(dds_init, normalized = TRUE))
  raw_counts_init <- as.data.frame(count_matrix)
  sf <- sizeFactors(dds_init)

  sample_ids <- colnames(count_matrix)
  trt_ids <- sample_ids[coldata$condition == "trt"]
  untrt_ids <- sample_ids[coldata$condition == "untrt"]

  fit_trt <- compute_pc1_loading_table(norm_counts_init, trt_ids, preprocessing_label = "Normalized prior to eigenvector splitting")
  fit_untrt <- compute_pc1_loading_table(norm_counts_init, untrt_ids, preprocessing_label = "Normalized prior to eigenvector splitting")

  trt_high <- as.character(subset(fit_trt$loading_table, split_class == "high_loading")$feature_id)
  untrt_high <- as.character(subset(fit_untrt$loading_table, split_class == "high_loading")$feature_id)

  leading_edge_ids <- union(trt_high, untrt_high)
  remainder_ids <- setdiff(rownames(count_matrix), leading_edge_ids)

  if (!length(leading_edge_ids)) stop("Leading-edge dataset is empty.")

  list(
    fit_trt = fit_trt,
    fit_untrt = fit_untrt,
    size_factors = sf,
    normalized_counts = norm_counts_init,
    raw_counts = raw_counts_init,
    raw_dataset = count_matrix,
    leading_edge_dataset = count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset = count_matrix[remainder_ids, , drop = FALSE]
  )
}

# =============================================================================
# EVS SUPPORT PLOTS
# =============================================================================
plot_pca_scatter <- function(pca_fit, dataset_label, group_label) {
  pca_var <- pca_fit$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)

  cond_key <- ifelse(group_label %in% c("control", "untrt"), "untrt", "trt")
  pca_df <- data.frame(
    Sample = rownames(pca_fit$x),
    PC1 = pca_fit$x[, 1],
    PC2 = pca_fit$x[, 2],
    Condition = cond_key
  )

  ggplot(pca_df, aes(PC1, PC2, label = Sample, shape = Condition, fill = Condition)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey70") +
    geom_point(size = 3.0, colour = "white", stroke = 0.55) +
    geom_text_repel(size = 2.0, max.overlaps = 8, force = 1.0,
                    box.padding = 0.22, point.padding = 0.10, min.segment.length = 0) +
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition") +
    scale_fill_manual(values = condition_fills, labels = condition_labels, name = "Condition") +
    labs(
      title = compact_title(paste(dataset_label, pretty_group_label(group_label), "PCA")),
      subtitle = paste0("PC1 = ", pca_var_per[1], "%; PC2 = ", pca_var_per[2], "%"),
      x = paste0("PC1 (", pca_var_per[1], "%)"),
      y = paste0("PC2 (", pca_var_per[2], "%)")
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()
}

plot_pc1_loading_rank <- function(loading_tbl, cutoff, dataset_label, group_label, top_n_used) {
  ggplot(loading_tbl, aes(rank, pc1_loading_abs)) +
    geom_line(linewidth = 0.4, color = "grey35") +
    geom_hline(yintercept = cutoff, color = "red", linewidth = 0.9) +
    labs(
      title = compact_title(paste(dataset_label, pretty_group_label(group_label), "PC1 loading rank")),
      subtitle = paste0("Top ", top_n_used, " features retained in EVS"),
      x = "Ranked PAS feature",
      y = "Absolute PC1 loading"
    ) +
    manuscript_theme()
}

# =============================================================================
# SECTION 3 - CORE ANALYSIS
# =============================================================================
run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
  design_formula <- make_design_formula(coldata)

  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = coldata,
    design = design_formula
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, betaPrior = FALSE)

  res <- results(dds, contrast = c("condition", "trt", "untrt"), alpha = alpha_level)
  res_strong <- results(dds, contrast = c("condition", "trt", "untrt"),
                        lfcThreshold = lfc_boundary, altHypothesis = "greaterAbs")
  res_weak <- results(dds, contrast = c("condition", "trt", "untrt"),
                      lfcThreshold = lfc_boundary, altHypothesis = "lessAbs")

  res_all_df <- as.data.frame(res)
  res_all_df$feature_id <- as.character(rownames(res_all_df))

  valid_stat <- is.finite(res_all_df$stat) & !is.na(res_all_df$stat)
  stat_vec <- as.numeric(res_all_df$stat[valid_stat])
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  if (length(stat_vec) < 5) stop(sprintf("[%s] fewer than 5 finite Wald statistics available.", dataset_name))

  message(sprintf("[%s] Mean Wald stat sent to fdrtool: %.4f", dataset_name, mean(stat_vec, na.rm = TRUE)))
  fdr_fit <- run_empirical_null_fdrtool(stat_vec, dataset_name = dataset_name)

  res_df <- res_all_df
  res_df$empirical_p <- NA_real_
  res_df$empirical_q <- NA_real_
  res_df$lfdr <- NA_real_
  res_df$empirical_p[valid_stat] <- as.numeric(fdr_fit$pval)
  res_df$empirical_q[valid_stat] <- as.numeric(fdr_fit$qval)
  res_df$lfdr[valid_stat] <- as.numeric(fdr_fit$lfdr)

  res_df$empirical_bh <- NA_real_
  valid_empirical <- is.finite(res_df$empirical_p) & !is.na(res_df$empirical_p)
  if (any(valid_empirical)) {
    res_df$empirical_bh[valid_empirical] <- p.adjust(res_df$empirical_p[valid_empirical], method = "BH")
  }

  res_df$pval <- res_df$empirical_p
  res_df$padjc <- res_df$empirical_bh
  res_df$qval <- res_df$empirical_q

  coef_name <- get_condition_coef(dds)
  shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)
  shr_df <- as.data.frame(shr)
  shr_df$feature_id <- as.character(rownames(shr_df))

  res_df <- dplyr::left_join(
    res_df,
    shr_df[, c("feature_id", "log2FoldChange")],
    by = "feature_id",
    suffix = c("", "_shrunk")
  )
  colnames(res_df)[colnames(res_df) == "log2FoldChange_shrunk"] <- "lfc_shrunk"

  hc_p_threshold_dataset <- safe_hc_thresh(res_df$empirical_p, dataset_name = dataset_name)
  res_df$gene_empirical_pvalue <- res_df$empirical_p
  res_df$HBFSS <- abs(res_df$lfc_shrunk * log10(pmax(res_df$gene_empirical_pvalue, 1e-300)))

  if (is.na(hc_p_threshold_dataset) || hc_p_threshold_dataset <= 0 || hc_p_threshold_dataset >= 1) {
    hbfss_threshold_dataset <- NA_real_
    res_df$HBFSS_core_pass <- FALSE
    res_df$HBFSS_significant <- FALSE
  } else {
    hbfss_threshold_dataset <- abs(log10(hc_p_threshold_dataset)) * lfc_boundary
    res_df$HBFSS_core_pass <- res_df$HBFSS >= hbfss_threshold_dataset
    res_df$HBFSS_significant <- res_df$HBFSS_core_pass
  }

  res_df$regulation_direction <- ifelse(
    is.na(res_df$lfc_shrunk), NA_character_,
    ifelse(res_df$lfc_shrunk > 0, "upregulated",
           ifelse(res_df$lfc_shrunk < 0, "downregulated", "no_change"))
  )

  res_df$raw_lfc_pass <- !is.na(res_df$log2FoldChange) & (abs(res_df$log2FoldChange) >= lfc_boundary)
  res_df$shrunk_lfc_pass <- !is.na(res_df$lfc_shrunk) & (abs(res_df$lfc_shrunk) >= lfc_boundary)
  res_df$standard_significant <- !is.na(res_df$padj) &
    (res_df$padj < alpha_level) &
    res_df$raw_lfc_pass &
    res_df$shrunk_lfc_pass

  res_strong_df <- as.data.frame(res_strong)
  res_strong_df$feature_id <- as.character(rownames(res_strong_df))
  res_weak_df <- as.data.frame(res_weak)
  res_weak_df$feature_id <- as.character(rownames(res_weak_df))

  res_df <- dplyr::left_join(res_df, res_strong_df[, c("feature_id", "padj")], by = "feature_id", suffix = c("", "_strong"))
  res_df <- dplyr::left_join(res_df, res_weak_df[, c("feature_id", "padj")], by = "feature_id", suffix = c("", "_weak"))
  colnames(res_df)[colnames(res_df) == "padj_strong"] <- "resGA_padj"
  colnames(res_df)[colnames(res_df) == "padj_weak"] <- "resLA_padj"

  res_df$effect_class <- classify_effect_strength(res_df$resGA_padj, res_df$resLA_padj, alpha = alpha_level)

  hbfss_native_call <- (res_df$HBFSS_core_pass & (is.na(res_df$resLA_padj) | !(res_df$resLA_padj < 0.2))) |
    (!is.na(res_df$resGA_padj) & (res_df$resGA_padj < alpha_level))
  res_df$HBFSS_significant <- hbfss_native_call & res_df$raw_lfc_pass

  base_mean_vec <- res_df$baseMean[!is.na(res_df$baseMean)]

  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts$feature_id <- as.character(rownames(norm_counts))

  mm <- as.data.frame(mcols(dds))
  mm$feature_id <- as.character(rownames(mm))
  disp_cols_available <- intersect(c("feature_id", "dispGeneEst", "dispFit", "dispersion", "baseMean"), colnames(mm))
  disp_df <- mm[, disp_cols_available, drop = FALSE]

  annot_df$feature_id <- as.character(annot_df$feature_id)
  annot_df$gene_symbol <- as.character(annot_df$gene_symbol)

  final_df <- res_df %>%
    dplyr::left_join(annot_df, by = "feature_id") %>%
    dplyr::left_join(norm_counts, by = "feature_id") %>%
    dplyr::left_join(disp_df, by = "feature_id")

  if (!"baseMean" %in% colnames(final_df) && "baseMean.x" %in% colnames(final_df)) {
    final_df$baseMean <- final_df$baseMean.x
  }
  if (!"baseMean" %in% colnames(final_df) && "baseMean.y" %in% colnames(final_df)) {
    final_df$baseMean <- final_df$baseMean.y
  }

  final_df$neglog10_padj <- safe_neglog10(final_df$padj)
  final_df$neglog10_empirical_p <- safe_neglog10(final_df$empirical_p)
  final_df$neglog10_padjc <- safe_neglog10(final_df$empirical_bh)
  final_df$dataset_name <- dataset_name
  final_df$hc_p_threshold_dataset <- hc_p_threshold_dataset
  final_df$hbfss_threshold_dataset <- hbfss_threshold_dataset

  list(
    dds = dds,
    results = final_df,
    base_mean_vec = base_mean_vec,
    hc_p_threshold = hc_p_threshold_dataset,
    hbfss_threshold = hbfss_threshold_dataset
  )
}

assign_final_plot_classes <- function(df, alpha_level = 0.10, lfc_boundary = 1.0) {
  df <- as.data.frame(df)

  hc_thr <- suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1]))
  hbfss_thr <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))

  df$hc_pass <- !is.na(df$empirical_p) &
    is.finite(df$empirical_p) &
    is.finite(hc_thr) &
    !is.na(hc_thr) &
    hc_thr > 0 & hc_thr < 1 &
    (df$empirical_p <= hc_thr)

  df$weak_cnh_pass <- df$hc_pass &
    !is.na(df$resLA_padj) & (df$resLA_padj < alpha_level) &
    !is.na(df$lfc_shrunk) & (abs(df$lfc_shrunk) < lfc_boundary)

  df$strong_cnh_pass <- df$hc_pass &
    !is.na(df$resGA_padj) & (df$resGA_padj < alpha_level) &
    !is.na(df$lfc_shrunk) & (abs(df$lfc_shrunk) >= lfc_boundary)

  df$standard_pass <- df$hc_pass &
    !is.na(df$padj) & (df$padj < alpha_level) &
    !is.na(df$lfc_shrunk) & (abs(df$lfc_shrunk) >= lfc_boundary)

  df$hbfss_pass <- df$hc_pass &
    is.finite(hbfss_thr) & !is.na(hbfss_thr) &
    !is.na(df$HBFSS) & is.finite(df$HBFSS) &
    (df$HBFSS >= hbfss_thr)

  df$plot_class <- "Background"
  df$plot_class[df$weak_cnh_pass] <- "Weak CNH"
  df$plot_class[df$strong_cnh_pass] <- "Strong CNH"
  df$plot_class[df$standard_pass] <- "Standard"
  df$plot_class[df$hbfss_pass] <- "HBFSS"
  df$plot_class <- factor(df$plot_class, levels = plot_class_levels)

  df$overlap_standard_hbfss <- df$standard_pass & df$hbfss_pass
  df
}

select_volcano_labels <- function(df, n_labels = 12) {
  df <- as.data.frame(df)
  df <- df[df$plot_class != "Background", , drop = FALSE]
  df <- df[!is.na(df$gene_symbol) & nzchar(trimws(df$gene_symbol)), , drop = FALSE]
  if (!nrow(df)) return(df[0, , drop = FALSE])

  df$label_priority <- dplyr::case_when(
    df$plot_class == "HBFSS" ~ 1,
    df$plot_class == "Standard" ~ 2,
    df$plot_class == "Strong CNH" ~ 3,
    df$plot_class == "Weak CNH" ~ 4,
    TRUE ~ 5
  )

  ord <- order(df$label_priority, -df$neglog10_empirical_p, -abs(df$lfc_shrunk), na.last = TRUE)
  df <- df[ord, , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol), , drop = FALSE]
  df[seq_len(min(n_labels, nrow(df))), , drop = FALSE]
}

build_single_volcano <- function(df, dataset_label, label_n = 12, show_legend = TRUE) {
  df <- assign_final_plot_classes(df = df, alpha_level = alpha_level, lfc_boundary = lfc_boundary)

  hc_thr <- suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1]))
  hbfss_thr <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))

  finite_x <- df$lfc_shrunk[is.finite(df$lfc_shrunk) & !is.na(df$lfc_shrunk)]
  finite_y <- df$neglog10_empirical_p[is.finite(df$neglog10_empirical_p) & !is.na(df$neglog10_empirical_p)]
  x_max <- if (length(finite_x)) max(2.5, ceiling(max(abs(finite_x), na.rm = TRUE) * 1.05)) else 3
  y_max <- if (length(finite_y)) max(4, ceiling(max(finite_y, na.rm = TRUE) * 1.08)) else 4
  hc_y <- if (is.finite(hc_thr) && !is.na(hc_thr) && hc_thr > 0 && hc_thr < 1) -log10(hc_thr) else NA_real_

  hbfss_curve_df <- make_hbfss_curve_df(xmax_abs = x_max, hbfss_threshold = hbfss_thr, ymax_plot = y_max, x_min_abs = 0.08)
  label_df <- select_volcano_labels(df, n_labels = label_n)

  caption_text <- paste0(
    "Weak=", sum(df$weak_cnh_pass, na.rm = TRUE),
    "  Strong=", sum(df$strong_cnh_pass, na.rm = TRUE),
    "  Std=", sum(df$standard_pass, na.rm = TRUE),
    "  HBFSS=", sum(df$hbfss_pass, na.rm = TRUE),
    "  Overlap=", sum(df$overlap_standard_hbfss, na.rm = TRUE)
  )

  p <- ggplot(df, aes(x = lfc_shrunk, y = neglog10_empirical_p)) +
    geom_point(aes(color = plot_class, shape = plot_class),
               size = 1.1, alpha = 0.86, stroke = 0.20, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_boundary, lfc_boundary),
               linewidth = 0.45, linetype = "dashed", color = plot_palette$threshold) +
    geom_vline(xintercept = 0, linewidth = 0.28, color = "grey55") +
    scale_color_manual(values = plot_class_colors, breaks = plot_class_levels, drop = FALSE, name = "Class") +
    scale_shape_manual(values = plot_class_shapes, breaks = plot_class_levels, drop = FALSE, name = "Class") +
    labs(
      title = dataset_label,
      subtitle = "x = shrunken log2FC   ·   y = -log10(empirical p)",
      x = "Shrunken log2FC",
      y = expression(-log[10]("Empirical p")),
      caption = caption_text
    ) +
    coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, y_max), clip = "off") +
    manuscript_theme() +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      plot.caption = element_text(hjust = 0.5)
    )

  if (is.finite(hc_y) && !is.na(hc_y)) {
    p <- p +
      geom_hline(yintercept = hc_y, linewidth = 0.45, linetype = "dotted", color = plot_palette$threshold) +
      annotate("text", x = -x_max * 0.88, y = min(y_max - 0.20, hc_y + 0.16),
               label = paste0("HC=", signif(hc_thr, 3)), color = plot_palette$threshold,
               size = 2.6, hjust = 0)
  }

  if (!is.null(hbfss_curve_df) && nrow(hbfss_curve_df) > 0) {
    p <- p +
      geom_line(data = hbfss_curve_df, aes(x = x, y = y),
                inherit.aes = FALSE, linewidth = 0.75, color = plot_palette$hbfss) +
      annotate("text", x = x_max * 0.56, y = min(y_max - 0.20, max(0.8, hbfss_thr + 0.20)),
               label = paste0("HBFSS=", signif(hbfss_thr, 3)), color = plot_palette$hbfss,
               size = 2.6, hjust = 0)
  }

  if (nrow(label_df) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(label = gene_symbol),
        size = 1.9,
        seed = 1,
        max.overlaps = 20,
        box.padding = 0.22,
        point.padding = 0.10,
        min.segment.length = 0,
        segment.alpha = 0.55,
        segment.size = 0.20
      )
  }

  p
}

build_compare_volcano_panel <- function(raw_df, lead_df, rem_df, comparison_name) {
  p_raw <- build_single_volcano(raw_df, "Raw", label_n = 7, show_legend = TRUE)
  p_lead <- build_single_volcano(lead_df, "Lead", label_n = 7, show_legend = FALSE)
  p_rem <- build_single_volcano(rem_df, "Rem", label_n = 7, show_legend = FALSE)

  legend_grob <- extract_legend_grob(p_raw)
  p_raw <- p_raw + theme(legend.position = "none")

  top_row <- arrangeGrob(
    grobs = list(p_raw, p_lead, p_rem),
    ncol = 3,
    top = textGrob(paste0(comparison_name, " volcano"), gp = gpar(fontface = "bold", cex = 1.03))
  )

  if (is.null(legend_grob)) return(top_row)

  arrangeGrob(
    grobs = list(top_row, legend_grob),
    ncol = 1,
    heights = unit.c(unit(1, "npc") - unit(0.52, "in"), unit(0.52, "in"))
  )
}

build_single_dispersion <- function(df, dataset_label, show_legend = TRUE) {
  req <- c("baseMean", "dispersion")
  if (!all(req %in% colnames(df))) {
    stop("build_single_dispersion() requires columns: baseMean, dispersion")
  }

  plot_df <- assign_final_plot_classes(df)
  plot_df <- plot_df[
    is.finite(plot_df$baseMean) & !is.na(plot_df$baseMean) &
      is.finite(plot_df$dispersion) & !is.na(plot_df$dispersion),
    , drop = FALSE
  ]

  ggplot(plot_df, aes(x = baseMean, y = dispersion)) +
    geom_point(aes(color = plot_class, shape = plot_class),
               size = 0.8, alpha = 0.72, stroke = 0.18, na.rm = TRUE) +
    scale_x_log10(labels = scales::label_number(accuracy = 0.1)) +
    scale_y_log10(labels = scales::label_number(accuracy = 0.01)) +
    scale_color_manual(values = plot_class_colors, breaks = plot_class_levels, drop = FALSE, name = "Class") +
    scale_shape_manual(values = plot_class_shapes, breaks = plot_class_levels, drop = FALSE, name = "Class") +
    labs(title = dataset_label, subtitle = "Final dispersion vs mean", x = "baseMean", y = "Dispersion") +
    manuscript_theme() +
    theme(legend.position = if (show_legend) "bottom" else "none")
}

build_compare_dispersion_panel <- function(raw_df, lead_df, rem_df, comparison_name) {
  p_raw <- build_single_dispersion(raw_df, "Raw", show_legend = TRUE)
  p_lead <- build_single_dispersion(lead_df, "Lead", show_legend = FALSE)
  p_rem <- build_single_dispersion(rem_df, "Rem", show_legend = FALSE)

  legend_grob <- extract_legend_grob(p_raw)
  p_raw <- p_raw + theme(legend.position = "none")

  top_row <- arrangeGrob(
    grobs = list(p_raw, p_lead, p_rem),
    ncol = 3,
    top = textGrob(paste0(comparison_name, " dispersion"), gp = gpar(fontface = "bold", cex = 1.03))
  )

  if (is.null(legend_grob)) return(top_row)

  arrangeGrob(
    grobs = list(top_row, legend_grob),
    ncol = 1,
    heights = unit.c(unit(1, "npc") - unit(0.52, "in"), unit(0.52, "in"))
  )
}

# =============================================================================
# SECTION 4 - RUN ONE COMPARISON
# =============================================================================
run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  cmp_dir <- file.path(output_dir, comparison_name)
  tab_dir <- file.path(cmp_dir, "tables")
  fig_raw_dir <- file.path(cmp_dir, "fig_raw")
  fig_le_dir <- file.path(cmp_dir, "fig_lead")
  fig_rem_dir <- file.path(cmp_dir, "fig_rem")
  for (d in c(cmp_dir, tab_dir, fig_raw_dir, fig_le_dir, fig_rem_dir)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  evs <- build_eigenvector_split(count_matrix, coldata)

  save_plot(plot_pca_scatter(evs$fit_trt$pca_fit, comparison_name, "treatment"),
            file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_PCA.png")), width = 8, height = 6)
  save_plot(plot_pca_scatter(evs$fit_untrt$pca_fit, comparison_name, "control"),
            file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_PCA.png")), width = 8, height = 6)
  save_plot(plot_pc1_loading_rank(evs$fit_trt$loading_table, evs$fit_trt$cutoff, comparison_name, "treatment", evs$fit_trt$top_n_used),
            file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_Rank.png")), width = 8, height = 6)
  save_plot(plot_pc1_loading_rank(evs$fit_untrt$loading_table, evs$fit_untrt$cutoff, comparison_name, "control", evs$fit_untrt$top_n_used),
            file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_Rank.png")), width = 8, height = 6)

  dataset_list <- list(
    raw_dataset = evs$raw_dataset,
    leading_edge_dataset = evs$leading_edge_dataset,
    remainder_dataset = evs$remainder_dataset
  )
  dataset_fig_dirs <- list(
    raw_dataset = fig_raw_dir,
    leading_edge_dataset = fig_le_dir,
    remainder_dataset = fig_rem_dir
  )

  analysis_results <- list()

  for (nm in names(dataset_list)) {
    fig_subdir <- dataset_fig_dirs[[nm]]
    full_dataset_name <- paste(comparison_name, nm, sep = "_")

    fit <- run_core_analysis(
      count_mat = dataset_list[[nm]],
      coldata = coldata,
      dataset_name = full_dataset_name,
      annot_df = annot_df
    )

    df <- assign_final_plot_classes(fit$results, alpha_level = alpha_level, lfc_boundary = lfc_boundary)

    save_csv(df, tab_file(tab_dir, comparison_name, nm, "Results"))
    save_csv(subset(df, standard_pass), tab_file(tab_dir, comparison_name, nm, "Std"))
    save_csv(subset(df, hbfss_pass), tab_file(tab_dir, comparison_name, nm, "HBFSS"))
    save_csv(
      df %>% dplyr::select(
        feature_id, gene_symbol, lfc_shrunk, empirical_p, padj, HBFSS,
        hc_pass, weak_cnh_pass, strong_cnh_pass, standard_pass, hbfss_pass,
        plot_class, overlap_standard_hbfss, baseMean, dispersion
      ),
      tab_file(tab_dir, comparison_name, nm, "PlotLogic")
    )

    summary_row <- data.frame(
      comparison_name = comparison_name,
      dataset_name = full_dataset_name,
      n_features = nrow(df),
      hc_p_threshold = fit$hc_p_threshold,
      hbfss_threshold = fit$hbfss_threshold,
      n_weak = sum(df$weak_cnh_pass, na.rm = TRUE),
      n_strong = sum(df$strong_cnh_pass, na.rm = TRUE),
      n_standard = sum(df$standard_pass, na.rm = TRUE),
      n_hbfss = sum(df$hbfss_pass, na.rm = TRUE),
      n_overlap = sum(df$overlap_standard_hbfss, na.rm = TRUE)
    )
    save_csv(summary_row, tab_file(tab_dir, comparison_name, nm, "Summary"))

    save_plot(
      build_single_volcano(df, unname(dataset_short[nm]), label_n = 12, show_legend = TRUE),
      fig_file(fig_subdir, comparison_name, nm, "Volcano"),
      width = 7.2,
      height = 6.0
    )

    if ("baseMean" %in% colnames(df) && "dispersion" %in% colnames(df)) {
      save_plot(
        build_single_dispersion(df, unname(dataset_short[nm]), show_legend = TRUE),
        fig_file(fig_subdir, comparison_name, nm, "Disp"),
        width = 7.2,
        height = 6.0
      )
    }

    analysis_results[[nm]] <- list(results = df, summary = summary_row)
  }

  save_grob(
    build_compare_volcano_panel(
      raw_df = analysis_results[["raw_dataset"]]$results,
      lead_df = analysis_results[["leading_edge_dataset"]]$results,
      rem_df = analysis_results[["remainder_dataset"]]$results,
      comparison_name = comparison_name
    ),
    file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_Volcano.png")),
    width = 13.0,
    height = 5.8
  )

  if (all(vapply(analysis_results, function(x) all(c("baseMean", "dispersion") %in% colnames(x$results)), logical(1)))) {
    save_grob(
      build_compare_dispersion_panel(
        raw_df = analysis_results[["raw_dataset"]]$results,
        lead_df = analysis_results[["leading_edge_dataset"]]$results,
        rem_df = analysis_results[["remainder_dataset"]]$results,
        comparison_name = comparison_name
      ),
      file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_Disp.png")),
      width = 13.0,
      height = 5.8
    )
  }

  dplyr::bind_rows(lapply(analysis_results, `[[`, "summary"))
}

# =============================================================================
# BATCH EXECUTION
# =============================================================================
comparison_inputs <- lapply(seq_len(nrow(comparison_table)), function(i) {
  prepare_comparison_data(
    comparison_name = comparison_table$comparison_name[i],
    group1_prefix = comparison_table$group1_prefix[i],
    group2_prefix = comparison_table$group2_prefix[i],
    WTTS_Seq = WTTS_Seq,
    meta_all = meta_all
  )
})
names(comparison_inputs) <- comparison_table$comparison_name

all_summaries_list <- list()
failed_comparisons <- list()

for (cmp in names(comparison_inputs)) {
  message("\n=====================================================")
  message("Running comparison: ", cmp)
  message("=====================================================")

  input_obj <- comparison_inputs[[cmp]]

  out <- tryCatch(
    run_full_comparison_pipeline(
      comparison_name = input_obj$comparison_name,
      count_matrix = input_obj$count_matrix,
      coldata = input_obj$coldata,
      annot_df = OrigID_Symbol
    ),
    error = function(e) {
      failed_comparisons[[cmp]] <<- data.frame(
        comparison_name = cmp,
        error_message = conditionMessage(e)
      )
      NULL
    }
  )

  if (!is.null(out) && nrow(out) > 0) {
    all_summaries_list[[cmp]] <- out
  }
}

all_summaries <- if (length(all_summaries_list) > 0) dplyr::bind_rows(all_summaries_list) else data.frame()

if (nrow(all_summaries) > 0) {
  save_csv(all_summaries, file.path(output_dir, "Table_Overall_Summary.csv"))
}

if (length(failed_comparisons) > 0) {
  failed_df <- dplyr::bind_rows(failed_comparisons)
  save_csv(failed_df, file.path(output_dir, "Table_Failed.csv"))
}

cat("\n=====================================================\n")
cat("Pipeline complete.\n")
cat("Output directory:\n")
cat(normalizePath(output_dir), "\n")
cat("=====================================================\n\n")

if (nrow(all_summaries) > 0) {
  print(all_summaries)
} else {
  message("No comparison summaries were written. Check Table_Failed.csv for the exact error message.")
}
