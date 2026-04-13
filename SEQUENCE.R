#!/usr/bin/env Rscript

# =============================================================================
# FINAL MANUSCRIPT PIPELINE
# EVS + DESeq2 + empirical-null + HBFSS
#
# PURPOSE
# This script is the final manuscript-ready analysis pipeline for the WTTS-Seq
# PAS-level count matrix. It performs:
#
# 1. Pairwise comparison setup for four circadian contrasts
# 2. Eigenvector splitting (EVS) into original, leading-edge, and remainder sets
# 3. DESeq2 modeling with design = ~ condition
# 4. Empirical-null recalibration of DESeq2 Wald statistics using fdrtool
# 5. Higher-criticism thresholding on empirical-null p-values
# 6. HBFSS scoring using apeglm-shrunken effect sizes and empirical-null p-values
# 7. Manuscript figures, cross-dataset comparison figures, and export tables
#
# DESIGN PRINCIPLES
# - Unit of observation is PAS / feature_id (OrigID), not gene
# - Symbol is retained only as annotation
# - Every comparison is modeled as ~ condition
# - No alternate design is used in this final manuscript version
# - File names are short, sortable, and GitHub-friendly
# - Plotting code guards against NA/invalid HC thresholds
# - Volcano legends use actual plotted shapes and colors
#
# METHODS SUMMARY
# A. EVS stage
#    Features are ranked by absolute PC1 loading within treatment and control
#    separately. The top-N features in either condition define the leading-edge
#    union. Features low-loading in both conditions form the remainder.
#
# B. DESeq2 stage
#    Each dataset (original, leading-edge, remainder) is fit independently with
#    DESeq2 using design = ~ condition and contrast = trt / untrt.
#
# C. Empirical-null stage
#    Only finite DESeq2 Wald statistics are sent to fdrtool to estimate the
#    empirical-null p-value, q-value, and local FDR.
#
# D. Higher-criticism stage
#    hc.thresh(sort(empirical_p)) is applied directly to the empirical-null
#    p-values. This returns a dataset-specific empirical p threshold.
#
# E. HBFSS stage
#    HBFSS = | lfc_shrunk * log10(empirical_p) |
#    with threshold:
#      hbfss_threshold_dataset = abs(log10(hc_p_threshold_dataset)) * lfc_boundary
#
# F. Manuscript interpretation
#    - standard_significant = classical DESeq2 positive call
#    - HBFSS_significant    = empirical-null / hybrid positive call
#    - effect_class         = strong / weak / intermediate effect class
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
# SECTION 1 OF 4
# SETUP, INPUTS, METADATA, GLOBAL HELPERS
# =============================================================================

# -----------------------------------------------------------------------------
# USER INPUT
# -----------------------------------------------------------------------------
# Uses repo-local data/ path so the script runs from the repository root.
count_file <- file.path("data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")

# Optional TWAS overlap
run_twas_overlap <- FALSE
twas_file <- "3aTWAS_genes_of_11_brain_disorders.csv"

# DESeq2 / manuscript settings
alpha_level <- 0.10
lfc_boundary <- 1.0
top_n_target <- 5000L

# Plot settings
figure_dpi <- 320
base_theme_size <- 10
n_top_labels_standard <- 20
n_top_labels_hbfss <- 20

# Output root
output_dir <- file.path("exports", "manuscript_final_clean")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# SHORT EXPORT NAMES
# -----------------------------------------------------------------------------
dataset_short <- c(
  raw_dataset = "Raw",
  leading_edge_dataset = "Lead",
  remainder_dataset = "Rem"
)

fig_file <- function(dir, cmp, ds_key, tag) {
  file.path(dir, paste0("Figure_", cmp, "_", unname(dataset_short[ds_key]), "_", tag, ".png"))
}

tab_file <- function(dir, cmp, ds_key, tag) {
  file.path(dir, paste0("Table_", cmp, "_", unname(dataset_short[ds_key]), "_", tag, ".csv"))
}

# -----------------------------------------------------------------------------
# CONSISTENT MANUSCRIPT PALETTE
# -----------------------------------------------------------------------------
plot_palette <- list(
  background   = "#BDBDBD",
  threshold    = "#8C2D04",
  hbfss        = "#E66101",
  deseq2       = "#B2182B",
  empirical    = "#5E3C99",
  strong       = "#1F78B4",
  weak         = "#33A02C",
  intermediate = "#7F7F7F",
  overlap      = "#000000",
  histogram    = "#969696",
  control      = "#4D4D4D",
  treatment    = "#1F78B4"
)

# -----------------------------------------------------------------------------
# MINIMAL EMBEDDED SAMPLE METADATA
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# FOUR PAIRWISE COMPARISONS
# -----------------------------------------------------------------------------
comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix   = c("R0", "R2", "R4", "R8"),
  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# BASIC HELPERS
# -----------------------------------------------------------------------------
assert_required_columns <- function(df, required_cols, object_name = "data frame") {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns in ", object_name, ": ",
        paste(missing_cols, collapse = ", ")
      )
    )
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

compact_caption <- function(x, width = 118) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

save_csv <- function(df, path) {
  write.csv(df, file = path, row.names = FALSE)
}

save_plot <- function(p, path, width = 12.2, height = 9.1, dpi = figure_dpi, bg = "white") {
  ggsave(
    filename = path,
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = bg,
    limitsize = FALSE
  )
}

save_grob <- function(g, path, width = 15.8, height = 9.3, dpi = figure_dpi, bg = "white") {
  ggsave(
    filename = path,
    plot = g,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = bg,
    limitsize = FALSE
  )
}

# -----------------------------------------------------------------------------
# LABEL HELPERS
# -----------------------------------------------------------------------------
dataset_key_order <- c("raw_dataset", "leading_edge_dataset", "remainder_dataset")
dataset_key_labels <- c(
  raw_dataset = "Original dataset",
  leading_edge_dataset = "Leading-edge dataset",
  remainder_dataset = "Remainder dataset"
)

pretty_dataset_type <- function(dataset_key) {
  switch(
    dataset_key,
    raw_dataset = "Original dataset",
    leading_edge_dataset = "Leading-edge dataset",
    remainder_dataset = "Remainder dataset",
    dataset_key
  )
}

pretty_dataset_label <- function(dataset_name) {
  parts <- strsplit(dataset_name, "_", fixed = TRUE)[[1]]
  if (length(parts) < 4) return(dataset_name)
  comparison_name <- paste(parts[1], parts[2], sep = "_")
  dataset_key <- paste(parts[3:length(parts)], collapse = "_")
  paste(comparison_name, pretty_dataset_type(dataset_key), sep = " | ")
}

pretty_group_label <- function(group_label) {
  switch(
    group_label,
    trt = "Treatment",
    untrt = "Control",
    treatment = "Treatment",
    control = "Control",
    group_label
  )
}

preprocessing_panel_label <- function(x) {
  x <- as.character(x)[1]
  if (is.na(x) || !nzchar(x)) return("Preprocessing not specified")
  dplyr::case_when(
    x %in% c("DESeq2-normalized counts", "normalized", "Normalized prior to eigenvector splitting") ~
      "Normalized prior to eigenvector splitting",
    x %in% c("raw counts without DESeq2 normalization", "raw_counts", "Eigenvector splitting without prior normalization") ~
      "Eigenvector splitting without prior normalization",
    TRUE ~ x
  )
}

# -----------------------------------------------------------------------------
# SHARED THEME
# -----------------------------------------------------------------------------
manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title        = element_text(face = "bold", size = base_theme_size + 1, hjust = 0.5, lineheight = 1.00, margin = margin(b = 5)),
      plot.subtitle     = element_text(size = base_theme_size - 1, hjust = 0.5, lineheight = 1.00, margin = margin(b = 7)),
      plot.caption      = element_text(size = base_theme_size - 3, hjust = 0.5, colour = "grey30", lineheight = 0.98, margin = margin(t = 8)),
      axis.title        = element_text(face = "bold"),
      axis.text         = element_text(colour = "black"),
      legend.title      = element_text(face = "bold"),
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.text       = element_text(size = base_theme_size - 1),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(linewidth = 0.25, colour = "grey88"),
      plot.margin       = margin(t = 14, r = 22, b = 18, l = 18)
    )
}

plot_expand_xy <- function() {
  list(
    scale_x_continuous(expand = expansion(mult = c(0.09, 0.18))),
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.24)))
  )
}

# -----------------------------------------------------------------------------
# MODEL / STATISTICS HELPERS
# -----------------------------------------------------------------------------
make_design_formula <- function(coldata) {
  ~ condition
}

get_condition_coef <- function(dds) {
  rn <- resultsNames(dds)
  idx <- grep("^condition_", rn)
  if (length(idx) == 0) stop("Could not identify condition coefficient in resultsNames(dds).")
  if (length(idx) > 1) message("Multiple condition coefficients found; using: ", rn[idx[1]])
  rn[idx[1]]
}

classify_effect_strength <- function(res_strong_padj, res_weak_padj, alpha = alpha_level) {
  out <- rep("intermediate", length(res_strong_padj))
  out[!is.na(res_weak_padj) & res_weak_padj < alpha] <- "weak_effect"
  out[!is.na(res_strong_padj) & res_strong_padj < alpha] <- "strong_effect"
  out
}

run_empirical_null_fdrtool <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  stat_vec <- unname(stat_vec)

  if (length(stat_vec) < 5) {
    stop(sprintf("[%s] Fewer than 5 finite Wald statistics were available for fdrtool.", dataset_name))
  }

  fit <- tryCatch(
    fdrtool(
      stat_vec,
      statistic = "normal",
      plot = FALSE,
      verbose = FALSE,
      cutoff.method = "fndr",
      pct0 = 0.75
    ),
    error = function(e1) {
      message(sprintf("[%s] Primary fdrtool call failed: %s", dataset_name, conditionMessage(e1)))
      tryCatch(
        fdrtool(
          as.vector(stat_vec),
          statistic = "normal",
          plot = FALSE,
          verbose = FALSE,
          cutoff.method = "pct0",
          pct0 = 0.75
        ),
        error = function(e2) {
          stop(sprintf("[%s] fdrtool failed after retry: %s", dataset_name, conditionMessage(e2)))
        }
      )
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
    tryCatch(
      fdrtool::hc.thresh(as.vector(sorted_empirical_p)),
      error = function(e) {
        message(sprintf("[%s] hc.thresh failed: %s", dataset_name, conditionMessage(e)))
        NA_real_
      }
    )
  )

  out <- as.numeric(out[1])
  if (!is.finite(out) || is.na(out) || out <= 0 || out >= 1) return(NA_real_)
  out
}

resolve_top_n_cutoff <- function(sorted_values_desc, top_n = top_n_target) {
  n_total <- length(sorted_values_desc)
  if (n_total == 0) stop("resolve_top_n_cutoff() received an empty vector.")
  top_n_actual <- min(max(1L, as.integer(top_n)), n_total)
  cutoff_value <- sorted_values_desc[top_n_actual]
  cutoff_quantile <- 1 - (top_n_actual / n_total)
  list(
    top_n_actual = top_n_actual,
    cutoff_value = cutoff_value,
    cutoff_quantile = cutoff_quantile,
    n_total = n_total
  )
}

clean_gene_set <- function(x) {
  unique(tolower(trimws(x[!is.na(x) & x != ""])))
}

# -----------------------------------------------------------------------------
# INPUT COUNT MATRIX
# -----------------------------------------------------------------------------
if (!file.exists(count_file)) {
  stop(
    paste0(
      "Count file not found: ", count_file,
      "\nRun this script from the repository root, or correct count_file."
    )
  )
}

WTTS_Seq <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
WTTS_Seq <- as.data.frame(WTTS_Seq, stringsAsFactors = FALSE)
WTTS_Seq$OrigID <- as.character(WTTS_Seq$OrigID)
WTTS_Seq$Symbol <- as.character(WTTS_Seq$Symbol)

assert_required_columns(WTTS_Seq, c("OrigID", "Symbol"), object_name = "WTTS count file")
assert_required_columns(WTTS_Seq, meta_all$id, object_name = "WTTS count file sample columns")

WTTS_Seq <- WTTS_Seq[!is.na(WTTS_Seq$OrigID) & !is.na(WTTS_Seq$Symbol), , drop = FALSE]
sample_na <- rowSums(is.na(WTTS_Seq[, meta_all$id, drop = FALSE])) > 0
WTTS_Seq <- WTTS_Seq[!sample_na, , drop = FALSE]
rownames(WTTS_Seq) <- WTTS_Seq$OrigID

OrigID_Symbol <- unique(WTTS_Seq[, c("OrigID", "Symbol"), drop = FALSE])
colnames(OrigID_Symbol) <- c("feature_id", "gene_symbol")
OrigID_Symbol$feature_id <- as.character(OrigID_Symbol$feature_id)
OrigID_Symbol$gene_symbol <- as.character(OrigID_Symbol$gene_symbol)
OrigID_Symbol <- OrigID_Symbol %>%
  dplyr::mutate(gene_symbol = dplyr::if_else(is.na(gene_symbol), "", trimws(gene_symbol))) %>%
  dplyr::arrange(feature_id, dplyr::desc(gene_symbol != ""), gene_symbol) %>%
  dplyr::distinct(feature_id, .keep_all = TRUE) %>%
  dplyr::mutate(gene_symbol = dplyr::na_if(gene_symbol, ""))

if (run_twas_overlap) {
  TWAS_Seq <- read.csv(twas_file, header = TRUE, stringsAsFactors = FALSE)
  TWAS_Seq <- as.data.frame(TWAS_Seq)
  if (ncol(TWAS_Seq) < 4) stop("TWAS file must contain at least 4 columns.")
  TWAS_data <- TWAS_Seq[, c(1, 4), drop = FALSE]
  colnames(TWAS_data) <- c("source_id", "gene_symbol")
}

# =============================================================================
# SECTION 2 OF 4
# COMPARISON PREPARATION, EVS, PCA / EVS SUPPORT FIGURES
# =============================================================================

prepare_comparison_data <- function(comparison_name, group1_prefix, group2_prefix, WTTS_Seq, meta_all) {
  keep_ids <- grepl(paste0("^", group1_prefix, "_"), meta_all$id) |
    grepl(paste0("^", group2_prefix, "_"), meta_all$id)

  meta_sub <- meta_all[keep_ids, , drop = FALSE]
  coldata <- meta_sub[, c("condition"), drop = FALSE]
  sample_ids <- rownames(meta_sub)

  missing_samples <- setdiff(sample_ids, colnames(WTTS_Seq))
  if (length(missing_samples) > 0) {
    stop(
      paste(
        "Missing samples in WTTS file for", comparison_name, ":",
        paste(missing_samples, collapse = ", ")
      )
    )
  }

  count_sub <- WTTS_Seq[, sample_ids, drop = FALSE]
  stopifnot(all(colnames(count_sub) == rownames(coldata)))

  list(
    comparison_name = comparison_name,
    count_matrix = as.matrix(count_sub),
    coldata = coldata
  )
}

compute_pc1_loading_table <- function(value_df, sample_names, top_n = top_n_target,
                                      preprocessing_label = "Normalized prior to eigenvector splitting") {
  x <- as.matrix(value_df[, sample_names, drop = FALSE])
  pca_fit <- prcomp(t(x), scale. = FALSE, rank. = 2)

  loading_abs <- abs(pca_fit$rotation[, 1])
  loading_tbl <- data.frame(
    feature_id = names(loading_abs),
    pc1_loading_abs = unname(loading_abs),
    stringsAsFactors = FALSE
  )
  loading_tbl <- loading_tbl[order(loading_tbl$pc1_loading_abs, decreasing = TRUE), ]
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

  dds_init <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = design_formula
  )
  dds_init <- dds_init[rowSums(counts(dds_init)) > 0, ]
  dds_init <- estimateSizeFactors(dds_init)

  norm_counts_init <- as.data.frame(counts(dds_init, normalized = TRUE))
  raw_counts_init <- as.data.frame(count_matrix)
  sf <- sizeFactors(dds_init)

  sample_ids <- colnames(count_matrix)
  trt_ids <- sample_ids[coldata$condition == "trt"]
  untrt_ids <- sample_ids[coldata$condition == "untrt"]

  fit_trt <- compute_pc1_loading_table(
    norm_counts_init, trt_ids, top_n = top_n_target,
    preprocessing_label = "Normalized prior to eigenvector splitting"
  )
  fit_untrt <- compute_pc1_loading_table(
    norm_counts_init, untrt_ids, top_n = top_n_target,
    preprocessing_label = "Normalized prior to eigenvector splitting"
  )

  fit_trt_raw <- compute_pc1_loading_table(
    raw_counts_init, trt_ids, top_n = top_n_target,
    preprocessing_label = "Eigenvector splitting without prior normalization"
  )
  fit_untrt_raw <- compute_pc1_loading_table(
    raw_counts_init, untrt_ids, top_n = top_n_target,
    preprocessing_label = "Eigenvector splitting without prior normalization"
  )

  trt_high <- as.character(subset(fit_trt$loading_table, split_class == "high_loading")$feature_id)
  untrt_high <- as.character(subset(fit_untrt$loading_table, split_class == "high_loading")$feature_id)

  leading_edge_ids <- union(trt_high, untrt_high)
  remainder_ids <- setdiff(rownames(count_matrix), leading_edge_ids)

  if (length(leading_edge_ids) == 0) {
    stop("Leading-edge dataset is empty. Check sample mapping or top_n_target.")
  }

  list(
    fit_trt = fit_trt,
    fit_untrt = fit_untrt,
    fit_trt_raw = fit_trt_raw,
    fit_untrt_raw = fit_untrt_raw,
    leading_edge_ids = leading_edge_ids,
    remainder_ids = remainder_ids,
    size_factors = sf,
    normalized_counts = norm_counts_init,
    raw_counts = raw_counts_init,
    raw_dataset = count_matrix,
    leading_edge_dataset = count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset = count_matrix[remainder_ids, , drop = FALSE]
  )
}

condition_shapes <- c("untrt" = 21, "trt" = 24)
condition_fills  <- c("untrt" = plot_palette$control, "trt" = plot_palette$treatment)
condition_labels <- c("untrt" = "Control", "trt" = "Treatment")

plot_pca_scatter <- function(pca_fit, dataset_label, group_label,
                             preprocessing_label = "Normalized prior to eigenvector splitting") {
  pca_var <- pca_fit$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)

  cond_key <- ifelse(group_label %in% c("control", "untrt"), "untrt", "trt")
  pca_df <- data.frame(
    Sample = rownames(pca_fit$x),
    PC1 = pca_fit$x[, 1],
    PC2 = pca_fit$x[, 2],
    Condition = cond_key,
    stringsAsFactors = FALSE
  )

  ggplot(pca_df, aes(PC1, PC2, label = Sample, shape = Condition, fill = Condition)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey70") +
    geom_point(size = 3.0, colour = "white", stroke = 0.55) +
    geom_text_repel(size = 2.0, max.overlaps = 8, force = 1.0, box.padding = 0.22, point.padding = 0.10, min.segment.length = 0) +
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition") +
    scale_fill_manual(values = condition_fills, labels = condition_labels, name = "Condition") +
    labs(
      title = compact_title(paste(dataset_label, "|", pretty_group_label(group_label), "PCA"), width = 44),
      subtitle = paste0(
        preprocessing_panel_label(preprocessing_label),
        " · Samples colored by condition. PC1 = ",
        pca_var_per[1], "%; PC2 = ", pca_var_per[2],
        "%. This panel compares preprocessing before eigenvector splitting."
      ),
      x = paste0("PC1 (", pca_var_per[1], "%)"),
      y = paste0("PC2 (", pca_var_per[2], "%)")
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    guides(
      fill = guide_legend(override.aes = list(size = 3.0, shape = 21, colour = "white")),
      shape = guide_legend(override.aes = list(size = 3.0, fill = "grey70", colour = "white"))
    )
}

plot_pc1_loading_rank <- function(loading_tbl, cutoff, dataset_label, group_label,
                                  top_n_used = top_n_target, cutoff_quantile = NA_real_,
                                  preprocessing_label = "Normalized prior to eigenvector splitting") {
  quantile_label <- if (is.finite(cutoff_quantile)) {
    paste0("Upper-tail quantile for this rank cutoff = ", signif(cutoff_quantile, 4))
  } else {
    "Corresponding upper-tail quantile = NA"
  }

  ggplot(loading_tbl, aes(rank, pc1_loading_abs)) +
    geom_line(linewidth = 0.4, color = "grey35") +
    geom_hline(yintercept = cutoff, color = "red", linewidth = 0.9) +
    annotate(
      "text",
      x = max(loading_tbl$rank) * 0.72,
      y = cutoff,
      label = paste0(
        "Absolute PC1 loading cutoff for top ", top_n_used,
        " sites = ", signif(cutoff, 4), "\n", quantile_label
      ),
      color = "red", vjust = -0.8, size = 4
    ) +
    labs(
      title = compact_title(paste(dataset_label, "|", pretty_group_label(group_label), "PC1 loading rank"), width = 52),
      subtitle = paste(
        "Absolute PC1 loading ranked within condition; the red line marks the leading-edge cutoff used for eigenvector splitting.",
        preprocessing_panel_label(preprocessing_label),
        "The reported upper-tail quantile shows how extreme that cutoff is within the ranked loading distribution."
      ),
      x = "Ranked PAS feature",
      y = "Absolute PC1 loading"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()
}

plot_eigenvector_histograms <- function(loading_tbl, cutoff, dataset_label, group_label,
                                        preprocessing_label = "Normalized prior to eigenvector splitting") {
  full_vals <- loading_tbl$pc1_loading_abs
  lead_vals <- loading_tbl$pc1_loading_abs[loading_tbl$pc1_loading_abs >= cutoff]
  rem_vals  <- loading_tbl$pc1_loading_abs[loading_tbl$pc1_loading_abs < cutoff]
  cutoff_tx <- -log(pmax(cutoff * 100, 1e-6))

  p1 <- ggplot(data.frame(x = -log(pmax(full_vals * 100, 1e-6))), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    geom_vline(xintercept = cutoff_tx, color = "red", linewidth = 0.9) +
    labs(title = "All ranked sites", x = "-log(|PC1 loading| × 100)", y = "Count") +
    manuscript_theme()

  p2 <- ggplot(data.frame(x = -log(pmax(rem_vals * 100, 1e-6))), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    geom_vline(xintercept = cutoff_tx, color = "red", linewidth = 0.9) +
    labs(title = "Remainder sites", x = "-log(|PC1 loading| × 100)", y = "Count") +
    manuscript_theme()

  p3 <- ggplot(data.frame(x = -log(pmax(lead_vals * 100, 1e-6))), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    geom_vline(xintercept = cutoff_tx, color = "red", linewidth = 0.9) +
    labs(title = "Leading-edge candidate sites", x = "-log(|PC1 loading| × 100)", y = "Count") +
    manuscript_theme()

  arrangeGrob(
    p1, p2, p3, ncol = 3,
    top = textGrob(
      paste(dataset_label, "|", pretty_group_label(group_label), "eigenvector loading distributions"),
      gp = gpar(fontface = "bold", cex = 1.15)
    )
  )
}

compute_mean_expression_table <- function(raw_counts, coldata) {
  sample_ids <- colnames(raw_counts)
  trt_ids <- sample_ids[coldata$condition == "trt"]
  untrt_ids <- sample_ids[coldata$condition == "untrt"]

  data.frame(
    feature_id = as.character(rownames(raw_counts)),
    mean_trt = rowMeans(raw_counts[, trt_ids, drop = FALSE]),
    mean_untrt = rowMeans(raw_counts[, untrt_ids, drop = FALSE]),
    mean_all = rowMeans(raw_counts),
    stringsAsFactors = FALSE
  )
}

plot_mean_histogram_panel <- function(df_means, base_mean_vec, dataset_name) {
  p1 <- ggplot(data.frame(x = safe_log10(df_means$mean_trt + 1)), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    labs(title = paste(dataset_name, "Treatment mean"), x = "log10(mean + 1)", y = "Count") +
    manuscript_theme()

  p2 <- ggplot(data.frame(x = safe_log10(df_means$mean_untrt + 1)), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    labs(title = paste(dataset_name, "Control mean"), x = "log10(mean + 1)", y = "Count") +
    manuscript_theme()

  p3 <- ggplot(data.frame(x = safe_log10(df_means$mean_all + 1)), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    labs(title = paste(dataset_name, "Pooled mean"), x = "log10(mean + 1)", y = "Count") +
    manuscript_theme()

  p4 <- ggplot(data.frame(x = safe_log10(base_mean_vec + 1)), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    labs(title = paste(dataset_name, "DESeq2 baseMean"), x = "log10(baseMean + 1)", y = "Count") +
    manuscript_theme()

  arrangeGrob(
    p1, p2, p3, p4, ncol = 2,
    top = textGrob(paste(dataset_name, "Mean-Expression Histograms"),
                   gp = gpar(fontface = "bold", cex = 1.2))
  )
}

# =============================================================================
# SECTION 3 OF 4
# DESEQ2, EMPIRICAL-NULL, HBFSS, VOLCANOES, DISPERSION
# =============================================================================

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
  design_formula <- make_design_formula(coldata)

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = design_formula
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, betaPrior = FALSE)

  res <- results(dds, contrast = c("condition", "trt", "untrt"), alpha = alpha_level)
  res_strong <- results(dds, contrast = c("condition", "trt", "untrt"), lfcThreshold = lfc_boundary, altHypothesis = "greaterAbs")
  res_weak <- results(dds, contrast = c("condition", "trt", "untrt"), lfcThreshold = lfc_boundary, altHypothesis = "lessAbs")

  res_all_df <- as.data.frame(res)
  res_all_df$feature_id <- as.character(rownames(res_all_df))

  valid_stat <- is.finite(res_all_df$stat) & !is.na(res_all_df$stat)
  stat_vec <- as.numeric(res_all_df$stat[valid_stat])
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]

  if (length(stat_vec) < 5) {
    stop(sprintf("[%s] Fewer than 5 finite Wald statistics were available for fdrtool.", dataset_name))
  }

  stat_mean <- mean(stat_vec, na.rm = TRUE)
  message(sprintf("[%s] Mean Wald stat sent to fdrtool: %.4f", dataset_name, stat_mean))

  fdr_fit <- run_empirical_null_fdrtool(stat_vec, dataset_name = dataset_name)

  res_df <- res_all_df
  n_valid <- sum(valid_stat)
  if (length(fdr_fit$pval) != n_valid || length(fdr_fit$qval) != n_valid || length(fdr_fit$lfdr) != n_valid) {
    stop(sprintf(
      "[%s] fdrtool output length mismatch: n_valid=%d, length(pval)=%d, length(qval)=%d, length(lfdr)=%d.",
      dataset_name, n_valid, length(fdr_fit$pval), length(fdr_fit$qval), length(fdr_fit$lfdr)
    ))
  }

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
  colnames(res_df)[colnames(res_df) == "padj_strong"] <- "padj_strong_effect"
  colnames(res_df)[colnames(res_df) == "padj_weak"] <- "padj_weak_effect"

  res_df$effect_class <- classify_effect_strength(
    res_df$padj_strong_effect,
    res_df$padj_weak_effect,
    alpha = alpha_level
  )

  res_df$resGA_padj <- res_df$padj_strong_effect
  res_df$resLA_padj <- res_df$padj_weak_effect

  hbfss_native_call <- (
    (res_df$HBFSS_core_pass & (is.na(res_df$resLA_padj) | !(res_df$resLA_padj < 0.2))) |
      (!is.na(res_df$resGA_padj) & (res_df$resGA_padj < alpha_level))
  )
  res_df$HBFSS_significant <- hbfss_native_call & res_df$raw_lfc_pass

  base_mean_vec <- res_df$baseMean[!is.na(res_df$baseMean)]

  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts$feature_id <- as.character(rownames(norm_counts))

  mm <- as.data.frame(mcols(dds))
  mm$feature_id <- as.character(rownames(mm))

  disp_cols_available <- intersect(
    c("feature_id", "dispGeneEst", "dispFit", "dispersion", "dispIter", "baseMean", "dispOutlier"),
    colnames(mm)
  )
  disp_df <- mm[, disp_cols_available, drop = FALSE]

  if ("dispersion" %in% colnames(disp_df) && "baseMean" %in% colnames(disp_df)) {
    disp_df$IOD <- disp_df$dispersion * disp_df$baseMean
  } else {
    disp_df$IOD <- NA_real_
  }

  if ("baseMean" %in% colnames(disp_df)) {
    disp_df <- disp_df[, setdiff(colnames(disp_df), "baseMean"), drop = FALSE]
  }

  annot_df$feature_id <- as.character(annot_df$feature_id)
  annot_df$gene_symbol <- as.character(annot_df$gene_symbol)
  annot_df <- annot_df %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(is.na(gene_symbol), "", trimws(gene_symbol))) %>%
    dplyr::arrange(feature_id, dplyr::desc(gene_symbol != ""), gene_symbol) %>%
    dplyr::distinct(feature_id, .keep_all = TRUE) %>%
    dplyr::mutate(gene_symbol = dplyr::na_if(gene_symbol, ""))

  res_df$feature_id <- as.character(res_df$feature_id)
  norm_counts$feature_id <- as.character(norm_counts$feature_id)
  norm_counts <- norm_counts[!duplicated(norm_counts$feature_id), , drop = FALSE]
  disp_df$feature_id <- as.character(disp_df$feature_id)
  disp_df <- disp_df[!duplicated(disp_df$feature_id), , drop = FALSE]

  final_df <- res_df %>%
    dplyr::left_join(annot_df, by = "feature_id") %>%
    dplyr::left_join(norm_counts, by = "feature_id") %>%
    dplyr::left_join(disp_df, by = "feature_id")

  final_df$neglog10_padj <- safe_neglog10(final_df$padj)
  final_df$neglog10_empirical_p <- safe_neglog10(final_df$empirical_p)
  final_df$neglog10_padjc <- safe_neglog10(final_df$empirical_bh)
  final_df$dataset_name <- dataset_name
  final_df$hc_p_threshold_dataset <- hc_p_threshold_dataset
  final_df$hbfss_threshold_dataset <- hbfss_threshold_dataset
  final_df$Apeglm_L2FC <- final_df$lfc_shrunk
  final_df$pi_valueE <- final_df$HBFSS
  final_df$`empirical p-value / pval` <- final_df$empirical_p
  final_df$`empirical q-value / qval` <- final_df$empirical_q
  final_df$`empirical BH / padjc` <- final_df$empirical_bh

  preferred_cols <- c(
    "dataset_name", "feature_id", "gene_symbol", "baseMean",
    "lfc_shrunk", "Apeglm_L2FC", "regulation_direction",
    "gene_empirical_pvalue", "pval", "padjc", "qval", "lfdr",
    "pvalue", "padj", "empirical_p", "empirical_q", "empirical_bh",
    "HBFSS", "pi_valueE", "hc_p_threshold_dataset", "hbfss_threshold_dataset",
    "resLA_padj", "resGA_padj",
    "standard_significant", "HBFSS_significant", "effect_class"
  )

  final_df <- final_df[, c(intersect(preferred_cols, names(final_df)), setdiff(names(final_df), preferred_cols)), drop = FALSE]

  list(
    dds = dds,
    results = final_df,
    base_mean_vec = base_mean_vec,
    hc_p_threshold = hc_p_threshold_dataset,
    hbfss_threshold = hbfss_threshold_dataset
  )
}

# -----------------------------------------------------------------------------
# VOLCANO CLASS ENCODING
# -----------------------------------------------------------------------------
build_reviewer_volcano_classes <- function(df) {
  df <- as.data.frame(df)

  df$effect_color_class <- dplyr::case_when(
    df$effect_class == "weak_effect" ~ "Weak effect",
    df$effect_class == "strong_effect" ~ "Strong effect",
    TRUE ~ "Intermediate effect"
  )
  df$effect_color_class <- factor(
    df$effect_color_class,
    levels = c("Strong effect", "Intermediate effect", "Weak effect")
  )

  df$method_call_class <- dplyr::case_when(
    !is.na(df$standard_significant) & df$standard_significant &
      !is.na(df$HBFSS_significant) & df$HBFSS_significant ~ "Strong effect + HBFSS",
    !is.na(df$standard_significant) & df$standard_significant ~ "Strong-effect only",
    !is.na(df$HBFSS_significant) & df$HBFSS_significant ~ "HBFSS candidate only",
    TRUE ~ "Empirical-null / neither"
  )
  df$method_call_class <- factor(
    df$method_call_class,
    levels = c("Empirical-null / neither", "Strong-effect only", "HBFSS candidate only", "Strong effect + HBFSS")
  )

  df$has_valid_gene_symbol <- !is.na(df$gene_symbol) & grepl("[A-Za-z0-9]", trimws(df$gene_symbol))
  df$gene_symbol_plot <- ifelse(df$has_valid_gene_symbol, trimws(df$gene_symbol), NA_character_)
  df$is_overlap <- !is.na(df$method_call_class) & as.character(df$method_call_class) == "Strong effect + HBFSS"
  df
}

effect_volcano_colors <- c(
  "Weak effect" = plot_palette$weak,
  "Intermediate effect" = plot_palette$intermediate,
  "Strong effect" = plot_palette$strong
)

method_call_colors <- c(
  "Empirical-null / neither" = "grey60",
  "Strong-effect only" = plot_palette$deseq2,
  "HBFSS candidate only" = plot_palette$hbfss,
  "Strong effect + HBFSS" = "#1E8449"
)

method_call_shapes <- c(
  "Empirical-null / neither" = 16,
  "Strong-effect only" = 17,
  "HBFSS candidate only" = 15,
  "Strong effect + HBFSS" = 18
)

method_fill_colors <- method_call_colors

volcano_guides <- function() {
  guides(
    fill = "none",
    color = guide_legend(
      order = 1,
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        size = 3.2,
        alpha = 1,
        stroke = 0.65,
        shape = unname(method_call_shapes)
      )
    ),
    shape = guide_legend(
      order = 1,
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        size = 3.2,
        alpha = 1,
        stroke = 0.65,
        color = unname(method_call_colors)
      )
    )
  )
}

volcano_label_layer <- function(lab_df) {
  if (!nrow(lab_df)) return(NULL)
  ggrepel::geom_text_repel(
    data = lab_df,
    aes(label = gene_symbol_plot),
    size = 1.8,
    seed = 1,
    max.overlaps = 20,
    force = 1.15,
    force_pull = 0.4,
    box.padding = 0.34,
    point.padding = 0.16,
    min.segment.length = 0,
    segment.alpha = 0.55,
    segment.size = 0.20,
    bg.color = "white",
    bg.r = 0.08
  )
}

select_volcano_labels <- function(df, y_col, n_labels = 20) {
  df <- as.data.frame(df)
  if (!nrow(df)) return(df[0, , drop = FALSE])

  df <- df[df$has_valid_gene_symbol, , drop = FALSE]
  if (!nrow(df)) return(df[0, , drop = FALSE])

  df$label_priority <- dplyr::case_when(
    df$method_call_class == "Strong effect + HBFSS" ~ 1,
    df$effect_color_class == "Strong effect" ~ 2,
    df$method_call_class == "Strong-effect only" ~ 3,
    df$method_call_class == "HBFSS candidate only" ~ 4,
    df$effect_color_class == "Intermediate effect" ~ 5,
    TRUE ~ 6
  )

  metric <- suppressWarnings(as.numeric(df[[y_col]]))
  metric[!is.finite(metric)] <- -Inf
  ord <- order(df$label_priority, -metric, -abs(df$lfc_shrunk), na.last = TRUE)
  df <- df[ord, , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol_plot), , drop = FALSE]
  df[seq_len(min(n_labels, nrow(df))), , drop = FALSE]
}

add_overlap_outline_gg <- function(plot_obj, df, x_col, y_col, size = 2.5, stroke = 0.8) {
  overlap_df <- df[df$is_overlap & is.finite(df[[x_col]]) & is.finite(df[[y_col]]), , drop = FALSE]
  if (!nrow(overlap_df)) return(plot_obj)

  plot_obj +
    geom_point(
      data = overlap_df,
      aes_string(x = x_col, y = y_col),
      inherit.aes = FALSE,
      shape = 21,
      size = size,
      stroke = stroke,
      fill = NA,
      color = plot_palette$overlap
    )
}

reviewer_volcano_subtitle <- function(y_label_text) {
  paste0(
    "Color + shape = interpretive tier. ",
    "Strong-effect only = classical DESeq2 significance; ",
    "HBFSS candidate only = intermediate alternative-hypothesis region; ",
    "Strong effect + HBFSS = agreement between both layers. ",
    y_label_text
  )
}

volcano_count_caption <- function(df) {
  method_counts <- table(factor(df$method_call_class, levels = levels(df$method_call_class)))
  effect_counts <- table(factor(df$effect_color_class, levels = levels(df$effect_color_class)))
  paste0(
    "Empirical-null / neither: n=", method_counts["Empirical-null / neither"],
    "  |  Strong-effect only: n=", method_counts["Strong-effect only"],
    "  |  HBFSS candidate only: n=", method_counts["HBFSS candidate only"],
    "  |  Strong effect + HBFSS: n=", method_counts["Strong effect + HBFSS"],
    "  ||  Strong effect: n=", effect_counts["Strong effect"],
    "  |  Intermediate effect: n=", effect_counts["Intermediate effect"],
    "  |  Weak effect: n=", effect_counts["Weak effect"]
  )
}

plot_empirical_p_histogram <- function(df, dataset_name, hc_p_threshold) {
  subtitle_text <- if (is.na(hc_p_threshold) || hc_p_threshold <= 0 || hc_p_threshold >= 1) {
    "Empirical-null p-values; HC threshold unavailable"
  } else {
    paste0("Empirical-null p-values with HC p-threshold = ", signif(hc_p_threshold, 4), " using direct hc.thresh(sort(empirical_p))")
  }

  p <- ggplot(df, aes(empirical_p)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, color = "white") +
    labs(
      title = compact_title(paste(compact_title(pretty_dataset_label(dataset_name)), "| Empirical-null p-value distribution")),
      subtitle = subtitle_text,
      x = "Empirical-null p-value",
      y = "Count"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()

  if (!is.na(hc_p_threshold) && is.finite(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < 1) {
    p <- p +
      geom_vline(xintercept = hc_p_threshold, color = plot_palette$hbfss, linewidth = 1) +
      annotate(
        "label",
        x = hc_p_threshold,
        y = Inf,
        label = paste0("HC threshold = ", signif(hc_p_threshold, 4)),
        vjust = 1.8,
        hjust = -0.02,
        fill = "white",
        label.size = 0.15,
        color = plot_palette$hbfss,
        size = 3.6
      )
  }

  p
}

plot_hc_profile <- function(df, dataset_name, hc_p_threshold) {
  pvals <- sort(df$empirical_p[is.finite(df$empirical_p) & !is.na(df$empirical_p) & df$empirical_p > 0 & df$empirical_p < 1])
  if (length(pvals) < 5) return(NULL)

  n <- length(pvals)
  i <- seq_len(n)
  v <- (i / n) * (1 - (i / n)) / n
  v[v == 0] <- min(v[v > 0])
  hc_score <- abs((i / n) - pvals) / sqrt(v)
  hc_df <- data.frame(rank = i, empirical_p = pvals, hc_score = hc_score)
  idx_peak <- which.max(hc_df$hc_score)
  peak_df <- hc_df[idx_peak, , drop = FALSE]

  p <- ggplot(hc_df, aes(empirical_p, hc_score)) +
    geom_line(linewidth = 0.5, colour = "grey30") +
    geom_point(data = peak_df, aes(empirical_p, hc_score), inherit.aes = FALSE, size = 2.4, colour = plot_palette$hbfss) +
    labs(
      title = compact_title(paste(compact_title(pretty_dataset_label(dataset_name)), "| Higher-criticism profile")),
      subtitle = paste0("Higher criticism computed from sorted empirical-null p-values with fdrtool::hc.thresh; returned threshold = ", signif(hc_p_threshold, 4), "."),
      x = "Sorted empirical p-value",
      y = "Higher-criticism score"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()

  if (!is.na(hc_p_threshold) && is.finite(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < 1) {
    p <- p + geom_vline(xintercept = hc_p_threshold, color = plot_palette$hbfss, linewidth = 1)
  }

  p
}

plot_standard_volcano <- function(df, dataset_name) {
  df <- build_reviewer_volcano_classes(df)
  lab_df <- select_volcano_labels(df, y_col = "neglog10_padj", n_labels = n_top_labels_standard)
  count_caption <- volcano_count_caption(df)

  p <- ggplot(df, aes(lfc_shrunk, neglog10_padj, color = method_call_class, shape = method_call_class)) +
    geom_point(alpha = 0.84, size = 1.6, stroke = 0.45) +
    geom_vline(xintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", linewidth = 0.6, colour = plot_palette$threshold) +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.45, colour = "grey45") +
    geom_hline(yintercept = -log10(alpha_level), linetype = "dashed", linewidth = 0.6, colour = plot_palette$threshold) +
    scale_color_manual(values = method_call_colors, drop = FALSE, name = "Interpretive tier") +
    scale_shape_manual(values = method_call_shapes, drop = FALSE, name = "Interpretive tier") +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "| Classical strong-effect volcano"), width = 42),
      subtitle = reviewer_volcano_subtitle("Standard volcano: y = -log10(adjusted p-value). Dashed lines mark the classical strong-effect decision boundary."),
      x = "Shrunken log2 fold change",
      y = expression(-log[10](padj)),
      caption = compact_caption(paste0(
        count_caption,
        "\nInterpretive tier encodes empirical-null vs intermediate HBFSS candidate vs classical strong effect  |  Dashed lines = ±",
        lfc_boundary, " and adjusted p-value ", alpha_level
      ))
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    plot_expand_xy() +
    volcano_guides() +
    theme(legend.position = "bottom", legend.box = "vertical", plot.margin = margin(18, 26, 22, 18))

  p <- add_overlap_outline_gg(p, df, "lfc_shrunk", "neglog10_padj", size = 2.15, stroke = 0.7)
  if (nrow(lab_df) > 0) p <- p + volcano_label_layer(lab_df)
  p
}

plot_hbfss_volcano <- function(df, dataset_name) {
  df <- build_reviewer_volcano_classes(df)
  lab_df <- select_volcano_labels(df, y_col = "neglog10_empirical_p", n_labels = n_top_labels_hbfss)
  count_caption <- volcano_count_caption(df)

  p <- ggplot(df, aes(lfc_shrunk, neglog10_empirical_p, color = method_call_class, shape = method_call_class)) +
    geom_point(alpha = 0.84, size = 1.6, stroke = 0.45) +
    scale_color_manual(values = method_call_colors, drop = FALSE, name = "Interpretive tier") +
    scale_shape_manual(values = method_call_shapes, drop = FALSE, name = "Interpretive tier") +
    geom_vline(xintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", linewidth = 0.6, colour = plot_palette$threshold) +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.45, colour = "grey45") +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "| HBFSS intermediate-candidate volcano"), width = 42),
      subtitle = reviewer_volcano_subtitle("HBFSS volcano: y = -log10(empirical p-value). This panel emphasizes the intermediate candidate region."),
      x = "Shrunken log2 fold change",
      y = expression(-log[10]("Empirical p-value")),
      caption = compact_caption(paste0(
        count_caption,
        "\nHBFSS = |shrunken log2FC × log10(empirical p)|  |  Empirical p-values derived from the empirical-null model; color + shape = interpretive tier"
      ))
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    plot_expand_xy() +
    volcano_guides() +
    theme(legend.position = "bottom", legend.box = "vertical", plot.margin = margin(18, 26, 22, 18))

  p <- add_overlap_outline_gg(p, df, "lfc_shrunk", "neglog10_empirical_p", size = 2.15, stroke = 0.7)
  if (nrow(lab_df) > 0) p <- p + volcano_label_layer(lab_df)
  p
}

plot_hbfss_distribution <- function(df, dataset_name, hbfss_threshold) {
  subtitle_text <- if (is.na(hbfss_threshold)) {
    "Observed HBFSS distribution; no usable HC-derived cutoff was returned."
  } else {
    paste0("HBFSS distribution; dataset-specific cutoff = ", signif(abs(hbfss_threshold), 4), ".")
  }

  p <- ggplot(df, aes(HBFSS)) +
    geom_histogram(bins = 70, fill = plot_palette$histogram, color = "white") +
    labs(
      title = compact_title(paste(compact_title(pretty_dataset_label(dataset_name)), "| HBFSS distribution")),
      subtitle = subtitle_text,
      x = "HBFSS = |Apeglm log2FC × log10(empirical p)|",
      y = "Count",
      caption = compact_caption("Orange line = dataset-specific HBFSS cutoff.")
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(plot.margin = margin(14, 18, 18, 14))

  if (!is.na(hbfss_threshold) && is.finite(hbfss_threshold)) {
    p <- p + geom_vline(xintercept = abs(hbfss_threshold), color = plot_palette$hbfss, linewidth = 0.9)
  }

  p
}

plot_composite_volcano <- function(df, dataset_name) {
  df <- build_reviewer_volcano_classes(df)
  lab_df <- select_volcano_labels(df, y_col = "neglog10_padjc", n_labels = n_top_labels_standard)
  count_caption <- volcano_count_caption(df)

  p <- ggplot(df, aes(lfc_shrunk, neglog10_padjc, color = method_call_class, shape = method_call_class)) +
    geom_point(alpha = 0.84, size = 1.55, stroke = 0.45) +
    scale_color_manual(values = method_call_colors, drop = FALSE, name = "Interpretive tier") +
    scale_shape_manual(values = method_call_shapes, drop = FALSE, name = "Interpretive tier") +
    geom_vline(xintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", linewidth = 0.6, colour = plot_palette$threshold) +
    geom_hline(yintercept = -log10(alpha_level), linetype = "dashed", linewidth = 0.6, colour = plot_palette$threshold) +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.45, colour = "grey45") +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "| Empirical-null adjusted volcano"), width = 46),
      subtitle = "Color + shape = method call; y = -log10(BH-adjusted empirical p-value).",
      x = "Shrunken log2 fold change (Apeglm)",
      y = expression(-log[10](padjc)),
      caption = paste0(
        count_caption,
        "\nColor + shape = method call  |  Dashed lines = ±", lfc_boundary,
        " and BH cutoff ", alpha_level
      )
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    plot_expand_xy() +
    volcano_guides() +
    theme(legend.position = "bottom", legend.box = "vertical", plot.margin = margin(18, 26, 22, 18))

  p <- add_overlap_outline_gg(p, df, "lfc_shrunk", "neglog10_padjc", size = 2.1, stroke = 0.68)
  if (nrow(lab_df) > 0) p <- p + volcano_label_layer(lab_df)
  p
}

plot_publication_volcano_panel <- function(df, dataset_name) {
  build <- build_reviewer_volcano_classes(df)
  top_df <- select_volcano_labels(build, y_col = "neglog10_empirical_p", n_labels = 20)

  hc_line <- NA_real_
  if ("hc_p_threshold_dataset" %in% colnames(build)) {
    hc_raw <- suppressWarnings(as.numeric(build$hc_p_threshold_dataset[1]))
    if (is.finite(hc_raw) && !is.na(hc_raw) && hc_raw > 0 && hc_raw < 1) {
      hc_line <- safe_neglog10(hc_raw)
    }
  }

  p <- ggplot(build, aes(lfc_shrunk, neglog10_empirical_p, color = method_call_class, shape = method_call_class)) +
    geom_vline(xintercept = c(-lfc_boundary, lfc_boundary), linewidth = 0.35, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.30, linetype = "solid", color = "grey55") +
    geom_point(alpha = 0.82, size = 1.5, stroke = 0.45) +
    scale_color_manual(values = method_call_colors, drop = FALSE, name = "Interpretive tier") +
    scale_shape_manual(values = method_call_shapes, drop = FALSE, name = "Interpretive tier") +
    labs(
      title = compact_title(pretty_dataset_label(dataset_name), width = 42),
      subtitle = "Color + shape = method call. Horizontal line = HC p-threshold when available.",
      x = "Shrunken log2 fold change",
      y = "-log10(empirical p-value)",
      caption = compact_caption(
        paste0(
          volcano_count_caption(build),
          "\nInterpretive tiers: empirical-null / neither, HBFSS-only candidate, classical strong-effect only, or overlap."
        )
      )
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    plot_expand_xy() +
    volcano_guides() +
    theme(legend.position = "bottom", legend.box = "vertical", plot.margin = margin(18, 26, 22, 18))

  if (is.finite(hc_line) && !is.na(hc_line)) {
    p <- p + geom_hline(yintercept = hc_line, linewidth = 0.35, linetype = "dashed", color = plot_palette$threshold)
  }

  p <- add_overlap_outline_gg(p, build, "lfc_shrunk", "neglog10_empirical_p", size = 2.05, stroke = 0.65)
  if (nrow(top_df) > 0) p <- p + volcano_label_layer(top_df)
  p
}

save_publication_volcano <- function(df, dataset_name, outfile) {
  p <- plot_publication_volcano_panel(df, dataset_name)
  save_plot(p, outfile)
  invisible(outfile)
}

plot_shrunken_ma <- function(df, dataset_name) {
  ggplot(df, aes(safe_log10(baseMean + 1), lfc_shrunk)) +
    geom_point(alpha = 0.6, size = 1.1, color = "grey40") +
    geom_hline(yintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(
      title = compact_title(paste(compact_title(pretty_dataset_label(dataset_name)), "| MA plot")),
      subtitle = "Shrunken effect size versus baseMean",
      x = expression(log[10]("baseMean + 1")),
      y = "Shrunken log2 fold change"
    ) +
    manuscript_theme()
}
plot_shrunken_MA <- plot_shrunken_ma

plot_dispersion_cloud <- function(dds, dataset_name, fig_subdir, cmp_short, ds_key) {
  outfile <- fig_file(fig_subdir, cmp_short, ds_key, "DispEst")
  png(outfile, width = 2400, height = 1900, res = figure_dpi)
  op <- par(no.readonly = TRUE)
  on.exit({par(op); dev.off()}, add = TRUE)
  par(mar = c(5.2, 5.2, 4.6, 2.2), mgp = c(2.8, 0.9, 0), cex.main = 1.2, cex.lab = 1.05)
  plotDispEsts(dds, main = paste(compact_title(pretty_dataset_label(dataset_name)), "Dispersion estimates"))
}

plot_dispersion_panel_for_dataset <- function(df, dataset_name) {
  df <- build_reviewer_volcano_classes(df)
  ggplot(df, aes(baseMean, dispersion, color = method_call_class, shape = method_call_class)) +
    geom_point(alpha = 0.55, size = 1.3, stroke = 0.35) +
    scale_x_log10(labels = label_number(accuracy = 0.1)) +
    scale_y_log10(labels = label_number(accuracy = 0.1)) +
    scale_color_manual(values = method_call_colors, drop = FALSE, name = "Interpretive tier") +
    scale_shape_manual(values = method_call_shapes, drop = FALSE, name = "Interpretive tier") +
    labs(
      title = compact_title(pretty_dataset_label(dataset_name), width = 42),
      subtitle = "Final dispersion vs mean count; color + shape = method call.",
      x = "baseMean (log10 scale)",
      y = "Final dispersion (log10 scale)"
    ) +
    manuscript_theme() +
    volcano_guides()
}

plot_dispersion_relationships <- function(df, dataset_name) {
  df <- build_reviewer_volcano_classes(df)
  plots <- list()

  if (all(c("baseMean", "dispersion") %in% names(df))) {
    plots[[length(plots) + 1]] <- ggplot(df, aes(baseMean, dispersion, fill = method_call_class, shape = method_call_class, color = effect_color_class)) +
      geom_point(alpha = 0.58, size = 1.05, stroke = 0.30) +
      scale_x_log10(labels = comma_format()) +
      scale_y_log10() +
      scale_fill_manual(values = method_fill_colors, drop = FALSE, name = "Interpretive tier") +
      scale_color_manual(values = effect_volcano_colors, drop = FALSE, name = "Effect size") +
      scale_shape_manual(values = method_call_shapes, drop = FALSE, name = "Interpretive tier") +
      labs(title = "Final dispersion vs mean", x = "baseMean", y = "Final dispersion") +
      manuscript_theme() + volcano_guides()
  }

  if (all(c("baseMean", "dispFit") %in% names(df))) {
    plots[[length(plots) + 1]] <- ggplot(df, aes(baseMean, dispFit, fill = method_call_class, shape = method_call_class, color = effect_color_class)) +
      geom_point(alpha = 0.58, size = 1.05, stroke = 0.30) +
      scale_x_log10(labels = comma_format()) +
      scale_y_log10() +
      scale_fill_manual(values = method_fill_colors, drop = FALSE, name = "Interpretive tier") +
      scale_color_manual(values = effect_volcano_colors, drop = FALSE, name = "Effect size") +
      scale_shape_manual(values = method_call_shapes, drop = FALSE, name = "Interpretive tier") +
      labs(title = "Trend dispersion vs mean", x = "baseMean", y = "Trend dispersion") +
      manuscript_theme() + volcano_guides()
  }

  if (all(c("baseMean", "dispGeneEst") %in% names(df))) {
    plots[[length(plots) + 1]] <- ggplot(df, aes(baseMean, dispGeneEst, fill = method_call_class, shape = method_call_class, color = effect_color_class)) +
      geom_point(alpha = 0.58, size = 1.05, stroke = 0.30) +
      scale_x_log10(labels = comma_format()) +
      scale_y_log10() +
      scale_fill_manual(values = method_fill_colors, drop = FALSE, name = "Interpretive tier") +
      scale_color_manual(values = effect_volcano_colors, drop = FALSE, name = "Effect size") +
      scale_shape_manual(values = method_call_shapes, drop = FALSE, name = "Interpretive tier") +
      labs(title = "Gene-wise dispersion vs mean", x = "baseMean", y = "Gene-wise dispersion") +
      manuscript_theme() + volcano_guides()
  }

  if (length(plots) == 0) return(NULL)

  do.call(
    arrangeGrob,
    c(plots, list(
      ncol = min(2, length(plots)),
      top = textGrob(paste(compact_title(pretty_dataset_label(dataset_name)), "Dispersion relationships"),
                     gp = gpar(fontface = "bold", cex = 1.2))
    ))
  )
}

dispersion_residual_section <- function(df, dataset_name, fig_subdir, tab_dir, cmp_short, ds_key) {
  if (!all(c("dispGeneEst", "dispFit", "dispersion") %in% colnames(df))) return(NULL)

  out <- data.frame(
    feature_id = df$feature_id,
    dispGeneEst = df$dispGeneEst,
    dispFit = df$dispFit,
    dispersion = df$dispersion,
    residual_fit_minus_final = df$dispFit - df$dispersion,
    residual_fit_minus_gene = df$dispFit - df$dispGeneEst,
    residual_final_minus_gene = df$dispersion - df$dispGeneEst,
    stringsAsFactors = FALSE
  )

  p1 <- ggplot(out, aes(residual_fit_minus_final)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, color = "white") +
    labs(title = "dispFit - final dispersion", x = "Residual", y = "Count") +
    manuscript_theme()

  p2 <- ggplot(out, aes(residual_fit_minus_gene)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, color = "white") +
    labs(title = "dispFit - gene-wise dispersion", x = "Residual", y = "Count") +
    manuscript_theme()

  p3 <- ggplot(out, aes(residual_final_minus_gene)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, color = "white") +
    labs(title = "final dispersion - gene-wise dispersion", x = "Residual", y = "Count") +
    manuscript_theme()

  panel <- arrangeGrob(
    p1, p2, p3, ncol = 3,
    top = textGrob(paste(compact_title(pretty_dataset_label(dataset_name)), "Dispersion residuals"),
                   gp = gpar(fontface = "bold", cex = 1.2))
  )

  save_grob(panel, fig_file(fig_subdir, cmp_short, ds_key, "DispResid"), width = 16.8, height = 5.8)
  save_csv(out, tab_file(tab_dir, cmp_short, ds_key, "DispResid"))
  out
}

# -----------------------------------------------------------------------------
# CROSS-DATASET COMPARISON PANELS
# -----------------------------------------------------------------------------
plot_empirical_histogram_for_panel <- function(df, dataset_name, hc_p_threshold) {
  p <- ggplot(df, aes(empirical_p)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, color = "white") +
    labs(
      title = compact_title(pretty_dataset_label(dataset_name), width = 42),
      subtitle = if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold) &&
                     hc_p_threshold > 0 && hc_p_threshold < 1) {
        paste0("Empirical-null p-values; HC threshold = ", signif(hc_p_threshold, 4))
      } else {
        "Empirical-null p-values; HC threshold unavailable"
      },
      x = "Empirical-null p-value",
      y = "Count"
    ) +
    manuscript_theme()

  if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < 1) {
    p <- p + geom_vline(xintercept = hc_p_threshold, color = plot_palette$threshold, linewidth = 0.9)
  }

  p
}

compute_dataset_pca_plot <- function(count_df, coldata, dataset_name, preprocessing = c("normalized", "raw_counts")) {
  preprocessing <- match.arg(preprocessing)
  design_formula <- make_design_formula(coldata)

  if (preprocessing == "normalized") {
    dds <- DESeqDataSetFromMatrix(countData = count_df, colData = coldata, design = design_formula)
    dds <- dds[rowSums(counts(dds)) > 0, ]
    dds <- estimateSizeFactors(dds)
    x <- counts(dds, normalized = TRUE)
    preprocessing_label <- "Normalized prior to eigenvector splitting"
  } else {
    x <- as.matrix(count_df)
    preprocessing_label <- "Eigenvector splitting without prior normalization"
  }

  pca_fit <- prcomp(t(x), scale. = FALSE, rank. = 2)
  pca_var <- pca_fit$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)

  pca_df <- data.frame(
    Sample = rownames(pca_fit$x),
    PC1 = pca_fit$x[, 1],
    PC2 = pca_fit$x[, 2],
    Condition = as.character(coldata[rownames(pca_fit$x), "condition"]),
    stringsAsFactors = FALSE
  )

  ggplot(pca_df, aes(PC1, PC2, label = Sample, shape = Condition, fill = Condition)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey70") +
    geom_point(size = 3, colour = "white", stroke = 0.55) +
    geom_text_repel(size = 2.0, max.overlaps = 8, force = 1.0, box.padding = 0.22, point.padding = 0.10, min.segment.length = 0) +
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition") +
    scale_fill_manual(values = condition_fills, labels = condition_labels, name = "Condition") +
    labs(
      title = compact_title(pretty_dataset_label(dataset_name)),
      subtitle = paste0(
        preprocessing_label,
        " · Samples colored by condition. PC1 = ",
        pca_var_per[1], "%; PC2 = ", pca_var_per[2],
        "%. This panel is for the methods comparison of preprocessing before eigenvector splitting."
      ),
      x = paste0("PC1 (", pca_var_per[1], "% variance)"),
      y = paste0("PC2 (", pca_var_per[2], "% variance)")
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    guides(
      fill = guide_legend(override.aes = list(size = 3.0, shape = 21, colour = "white")),
      shape = guide_legend(override.aes = list(size = 3.0, fill = "grey70", colour = "white"))
    )
}

save_cross_dataset_comparison_panels <- function(comparison_name, analysis_results, cmp_dir, dataset_list, coldata) {
  keys_present <- base::intersect(as.character(dataset_key_order), as.character(names(analysis_results)))
  if (length(keys_present) == 0) return(invisible(NULL))

  hist_grobs <- lapply(keys_present, function(k) {
    plot_empirical_histogram_for_panel(
      analysis_results[[k]]$results,
      analysis_results[[k]]$summary$dataset_name[1],
      analysis_results[[k]]$summary$hc_p_threshold[1]
    )
  })
  hist_panel <- do.call(arrangeGrob, c(hist_grobs, list(
    ncol = length(hist_grobs),
    top = textGrob(paste(comparison_name, "| Empirical-null p-value histograms"),
                   gp = gpar(fontface = "bold", cex = 1.15))
  )))
  save_grob(hist_panel, file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_EmpHist.png")),
            width = 6.4 * length(hist_grobs), height = 5.9)

  disp_grobs <- lapply(keys_present, function(k) {
    plot_dispersion_panel_for_dataset(
      analysis_results[[k]]$results,
      analysis_results[[k]]$summary$dataset_name[1]
    )
  })
  disp_panel <- do.call(arrangeGrob, c(disp_grobs, list(
    ncol = length(disp_grobs),
    top = textGrob(paste(comparison_name, "| Dispersion comparison across original, leading-edge, and remainder"),
                   gp = gpar(fontface = "bold", cex = 1.15))
  )))
  save_grob(disp_panel, file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_Disp.png")),
            width = 6.4 * length(disp_grobs), height = 5.9)

  std_volcano_grobs <- lapply(keys_present, function(k) {
    plot_standard_volcano(analysis_results[[k]]$results, analysis_results[[k]]$summary$dataset_name[1])
  })
  std_volcano_panel <- do.call(arrangeGrob, c(std_volcano_grobs, list(
    ncol = length(std_volcano_grobs),
    top = textGrob(paste(comparison_name, "| Classical strong-effect volcano comparison across original, leading-edge, and remainder"),
                   gp = gpar(fontface = "bold", cex = 1.10))
  )))
  save_grob(std_volcano_panel, file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_VolcanoStd.png")),
            width = 6.8 * length(std_volcano_grobs), height = 6.1)

  hbfss_volcano_grobs <- lapply(keys_present, function(k) {
    plot_hbfss_volcano(analysis_results[[k]]$results, analysis_results[[k]]$summary$dataset_name[1])
  })
  hbfss_volcano_panel <- do.call(arrangeGrob, c(hbfss_volcano_grobs, list(
    ncol = length(hbfss_volcano_grobs),
    top = textGrob(paste(comparison_name, "| HBFSS volcano comparison across original, leading-edge, and remainder"),
                   gp = gpar(fontface = "bold", cex = 1.10))
  )))
  save_grob(hbfss_volcano_panel, file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_VolcanoHBFSS.png")),
            width = 6.8 * length(hbfss_volcano_grobs), height = 6.1)

  publication_volcano_grobs <- lapply(keys_present, function(k) {
    plot_publication_volcano_panel(analysis_results[[k]]$results, analysis_results[[k]]$summary$dataset_name[1])
  })
  publication_volcano_panel <- do.call(arrangeGrob, c(publication_volcano_grobs, list(
    ncol = length(publication_volcano_grobs),
    top = textGrob(paste(comparison_name, "| Integrated publication volcano comparison across original, leading-edge, and remainder"),
                   gp = gpar(fontface = "bold", cex = 1.10))
  )))
  save_grob(publication_volcano_panel, file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_VolcanoPub.png")),
            width = 6.8 * length(publication_volcano_grobs), height = 6.1)

  pca_norm_grobs <- lapply(keys_present, function(k) {
    compute_dataset_pca_plot(dataset_list[[k]], coldata, analysis_results[[k]]$summary$dataset_name[1], preprocessing = "normalized")
  })
  pca_raw_grobs <- lapply(keys_present, function(k) {
    compute_dataset_pca_plot(dataset_list[[k]], coldata, analysis_results[[k]]$summary$dataset_name[1], preprocessing = "raw_counts")
  })

  pca_panel <- arrangeGrob(
    grobs = c(pca_norm_grobs, pca_raw_grobs),
    ncol = length(keys_present),
    top = textGrob(
      paste(
        comparison_name,
        "| PCA comparison across original, leading-edge, and remainder datasets.",
        "Top row: DESeq2 size-factor normalization before eigenvector splitting.",
        "Bottom row: raw counts without prior normalization before eigenvector splitting."
      ),
      gp = gpar(fontface = "bold", cex = 1.08)
    )
  )
  save_grob(pca_panel, file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_PCA.png")),
            width = 6.2 * length(keys_present), height = 10.5)

  summary_table <- dplyr::bind_rows(lapply(keys_present, function(k) {
    sm <- analysis_results[[k]]$summary
    data.frame(
      comparison_name = comparison_name,
      dataset_key = k,
      dataset_label = unname(dataset_key_labels[k]),
      dataset_name = sm$dataset_name[1],
      n_features = sm$n_features[1],
      n_standard_significant = sm$n_standard_significant[1],
      n_HBFSS_significant = sm$n_HBFSS_significant[1],
      n_overlap = sum(analysis_results[[k]]$results$standard_significant & analysis_results[[k]]$results$HBFSS_significant, na.rm = TRUE),
      n_strong_effect = sm$n_strong_effect[1],
      n_weak_effect = sm$n_weak_effect[1],
      hc_p_threshold = sm$hc_p_threshold[1],
      hbfss_threshold = sm$hbfss_threshold[1],
      top_n_target = top_n_target,
      stringsAsFactors = FALSE
    )
  }))
  save_csv(summary_table, file.path(cmp_dir, paste0("Table_", comparison_name, "_Compare_Summary.csv")))
  invisible(NULL)
}

# =============================================================================
# SECTION 4 OF 4
# BATCH EXECUTION, EXPORTS, FINAL SUMMARIES
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

  save_plot(plot_pca_scatter(evs$fit_trt$pca_fit, comparison_name, "treatment", evs$fit_trt$preprocessing_label),
            file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_PCA.png")), width = 9, height = 7)
  save_plot(plot_pca_scatter(evs$fit_untrt$pca_fit, comparison_name, "control", evs$fit_untrt$preprocessing_label),
            file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_PCA.png")), width = 9, height = 7)

  save_plot(plot_pc1_loading_rank(
    evs$fit_trt$loading_table, evs$fit_trt$cutoff,
    comparison_name, "treatment",
    top_n_used = evs$fit_trt$top_n_used,
    cutoff_quantile = evs$fit_trt$cutoff_quantile,
    preprocessing_label = evs$fit_trt$preprocessing_label
  ), file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_Rank.png")), width = 10, height = 7)

  save_plot(plot_pc1_loading_rank(
    evs$fit_untrt$loading_table, evs$fit_untrt$cutoff,
    comparison_name, "control",
    top_n_used = evs$fit_untrt$top_n_used,
    cutoff_quantile = evs$fit_untrt$cutoff_quantile,
    preprocessing_label = evs$fit_untrt$preprocessing_label
  ), file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_Rank.png")), width = 10, height = 7)

  trt_hist_grob <- plot_eigenvector_histograms(
    evs$fit_trt$loading_table, evs$fit_trt$cutoff, comparison_name, "treatment", evs$fit_trt$preprocessing_label
  )
  save_grob(trt_hist_grob, file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_Hist.png")), width = 16, height = 5)

  ctrl_hist_grob <- plot_eigenvector_histograms(
    evs$fit_untrt$loading_table, evs$fit_untrt$cutoff, comparison_name, "control", evs$fit_untrt$preprocessing_label
  )
  save_grob(ctrl_hist_grob, file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_Hist.png")), width = 16, height = 5)

  save_plot(plot_pca_scatter(evs$fit_trt_raw$pca_fit, comparison_name, "treatment", evs$fit_trt_raw$preprocessing_label),
            file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_PCA_Raw.png")), width = 9, height = 7)
  save_plot(plot_pca_scatter(evs$fit_untrt_raw$pca_fit, comparison_name, "control", evs$fit_untrt_raw$preprocessing_label),
            file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_PCA_Raw.png")), width = 9, height = 7)

  save_plot(plot_pc1_loading_rank(
    evs$fit_trt_raw$loading_table, evs$fit_trt_raw$cutoff,
    comparison_name, "treatment",
    top_n_used = evs$fit_trt_raw$top_n_used,
    cutoff_quantile = evs$fit_trt_raw$cutoff_quantile,
    preprocessing_label = evs$fit_trt_raw$preprocessing_label
  ), file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_Rank_Raw.png")), width = 10, height = 7)

  save_plot(plot_pc1_loading_rank(
    evs$fit_untrt_raw$loading_table, evs$fit_untrt_raw$cutoff,
    comparison_name, "control",
    top_n_used = evs$fit_untrt_raw$top_n_used,
    cutoff_quantile = evs$fit_untrt_raw$cutoff_quantile,
    preprocessing_label = evs$fit_untrt_raw$preprocessing_label
  ), file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_Rank_Raw.png")), width = 10, height = 7)

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
    full_dataset_name <- paste(comparison_name, nm, sep = "_")
    fig_subdir <- dataset_fig_dirs[[nm]]

    fit <- run_core_analysis(
      count_mat = dataset_list[[nm]],
      coldata = coldata,
      dataset_name = full_dataset_name,
      annot_df = annot_df
    )

    df <- fit$results

    save_csv(df, tab_file(tab_dir, comparison_name, nm, "Results"))
    save_csv(subset(df, standard_significant), tab_file(tab_dir, comparison_name, nm, "StdSig"))
    save_csv(subset(df, HBFSS_significant), tab_file(tab_dir, comparison_name, nm, "HBFSSSig"))
    save_csv(subset(df, effect_class == "strong_effect"), tab_file(tab_dir, comparison_name, nm, "Strong"))
    save_csv(subset(df, effect_class == "weak_effect"), tab_file(tab_dir, comparison_name, nm, "Weak"))

    summary_row <- data.frame(
      comparison_name = comparison_name,
      dataset_name = full_dataset_name,
      n_features = nrow(df),
      hc_p_threshold = fit$hc_p_threshold,
      hbfss_threshold = fit$hbfss_threshold,
      n_standard_significant = sum(df$standard_significant, na.rm = TRUE),
      n_HBFSS_significant = sum(df$HBFSS_significant, na.rm = TRUE),
      n_overlap_significant = sum(df$standard_significant & df$HBFSS_significant, na.rm = TRUE),
      n_strong_effect = sum(df$effect_class == "strong_effect", na.rm = TRUE),
      n_weak_effect = sum(df$effect_class == "weak_effect", na.rm = TRUE),
      top_n_target = top_n_target,
      trt_cutoff_quantile = evs$fit_trt$cutoff_quantile,
      ctrl_cutoff_quantile = evs$fit_untrt$cutoff_quantile,
      stringsAsFactors = FALSE
    )
    save_csv(summary_row, tab_file(tab_dir, comparison_name, nm, "Summary"))

    mean_df <- compute_mean_expression_table(dataset_list[[nm]], coldata)
    mean_panel <- plot_mean_histogram_panel(mean_df, fit$base_mean_vec, full_dataset_name)
    save_grob(mean_panel, fig_file(fig_subdir, comparison_name, nm, "MeanHist"), width = 14, height = 10)

    p_hist <- plot_empirical_p_histogram(df, full_dataset_name, fit$hc_p_threshold)
    p_std  <- plot_standard_volcano(df, full_dataset_name)
    p_hbf  <- plot_hbfss_volcano(df, full_dataset_name)
    p_hbfd <- plot_hbfss_distribution(df, full_dataset_name, fit$hbfss_threshold)
    p_hc   <- plot_hc_profile(df, full_dataset_name, fit$hc_p_threshold)
    p_comp <- plot_composite_volcano(df, full_dataset_name)
    p_ma   <- plot_shrunken_ma(df, full_dataset_name)

    save_plot(p_hist, fig_file(fig_subdir, comparison_name, nm, "EmpHist"))
    save_plot(p_std, fig_file(fig_subdir, comparison_name, nm, "VolcanoStd"))
    save_plot(p_hbf, fig_file(fig_subdir, comparison_name, nm, "VolcanoHBFSS"))
    save_plot(p_hbfd, fig_file(fig_subdir, comparison_name, nm, "HBFSSDist"))
    if (!is.null(p_hc)) save_plot(p_hc, fig_file(fig_subdir, comparison_name, nm, "HCProfile"))
    save_plot(p_comp, fig_file(fig_subdir, comparison_name, nm, "VolcanoEmpAdj"))
    save_publication_volcano(df, full_dataset_name, fig_file(fig_subdir, comparison_name, nm, "VolcanoPub"))
    save_plot(p_ma, fig_file(fig_subdir, comparison_name, nm, "MA"))

    plot_dispersion_cloud(fit$dds, full_dataset_name, fig_subdir, comparison_name, nm)

    disp_panel <- plot_dispersion_relationships(df, full_dataset_name)
    if (!is.null(disp_panel)) {
      save_grob(disp_panel, fig_file(fig_subdir, comparison_name, nm, "DispRel"), width = 14, height = 8)
    }

    dispersion_residual_section(df, full_dataset_name, fig_subdir, tab_dir, comparison_name, nm)

    summary_plot_list <- list(p_hist, p_std, p_hbf, p_hbfd, p_comp, p_ma)
    if (!is.null(p_hc)) summary_plot_list <- c(summary_plot_list, list(p_hc))

    summary_panel <- do.call(
      arrangeGrob,
      c(summary_plot_list, list(
        ncol = 2,
        top = textGrob(
          paste(pretty_dataset_label(full_dataset_name), "Summary panel"),
          gp = gpar(fontface = "bold", cex = 1.2)
        )
      ))
    )
    save_grob(summary_panel, fig_file(fig_subdir, comparison_name, nm, "Summary"), width = 16, height = 18)

    analysis_results[[nm]] <- list(
      dds = fit$dds,
      results = df,
      summary = summary_row,
      dataset_mat = dataset_list[[nm]],
      fig_subdir = fig_subdir
    )
  }

  tryCatch({
    save_cross_dataset_comparison_panels(comparison_name, analysis_results, cmp_dir, dataset_list, coldata)
  }, error = function(e) {
    warning(paste0(comparison_name, ": cross-dataset panel generation failed: ", conditionMessage(e)))
  })

  if (run_twas_overlap) {
    twas_genes <- clean_gene_set(TWAS_data$gene_symbol)

    get_twas_overlap <- function(result_df, dataset_nm, cmp_name, out_dir, mode = c("union", "standard_only", "HBFSS_only")) {
      mode <- match.arg(mode)

      sig_df <- switch(
        mode,
        union = subset(result_df, standard_significant | HBFSS_significant),
        standard_only = subset(result_df, standard_significant),
        HBFSS_only = subset(result_df, HBFSS_significant)
      )

      sig_df$gene_symbol_clean <- tolower(trimws(sig_df$gene_symbol))
      overlap_df <- subset(sig_df, gene_symbol_clean %in% twas_genes)

      summary_df <- data.frame(
        comparison_name = cmp_name,
        dataset_name = dataset_nm,
        selection_mode = mode,
        n_selected_features = nrow(sig_df),
        n_overlap_features = nrow(overlap_df),
        n_overlap_genes = length(unique(overlap_df$gene_symbol_clean)),
        stringsAsFactors = FALSE
      )

      mode_tag <- switch(
        mode,
        union = "Union",
        standard_only = "StdOnly",
        HBFSS_only = "HBFSSOnly"
      )

      write.csv(overlap_df, file.path(out_dir, paste0("Table_", dataset_nm, "_TWAS_", mode_tag, ".csv")), row.names = FALSE)
      write.csv(summary_df, file.path(out_dir, paste0("Table_", dataset_nm, "_TWAS_", mode_tag, "_Summary.csv")), row.names = FALSE)
      summary_df
    }

    twas_summaries <- dplyr::bind_rows(lapply(names(analysis_results), function(nm) {
      full_nm <- paste(comparison_name, nm, sep = "_")
      dplyr::bind_rows(
        get_twas_overlap(analysis_results[[nm]]$results, full_nm, comparison_name, tab_dir, "union"),
        get_twas_overlap(analysis_results[[nm]]$results, full_nm, comparison_name, tab_dir, "standard_only"),
        get_twas_overlap(analysis_results[[nm]]$results, full_nm, comparison_name, tab_dir, "HBFSS_only")
      )
    }))
    save_csv(twas_summaries, file.path(tab_dir, paste0("Table_", comparison_name, "_TWAS_Summary.csv")))
  }

  sm_list <- lapply(analysis_results, `[[`, "summary")
  sm_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, sm_list)
  if (length(sm_list) == 0) return(data.frame())
  dplyr::bind_rows(sm_list)
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
        error_message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      NULL
    }
  )

  if (!is.null(out) && nrow(out) > 0) {
    all_summaries_list[[cmp]] <- out
  }
}

all_summaries <- if (length(all_summaries_list) > 0) {
  dplyr::bind_rows(all_summaries_list)
} else {
  data.frame()
}

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
