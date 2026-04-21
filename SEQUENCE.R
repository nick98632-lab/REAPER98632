#!/usr/bin/env Rscript

# =============================================================================

# FINAL MANUSCRIPT PIPELINE

# EVS + DESeq2 + empirical-null HC + HBFSS

#

# This script is a full rewrite from scratch.

#

# Manuscript intent

# -----------------

# 1. Run four fixed pairwise comparisons:

#      RT0_ZT6, RT2_ZT8, RT4_ZT10, RT8_ZT14

#

# 2. Within each comparison, define three datasets:

#      Raw dataset

#      Leading-edge dataset

#      Remainder dataset

#

# 3. The leading-edge split is based on condition-specific absolute PC1

#    loadings computed after DESeq2 size-factor normalization. The leading edge

#    is the union of the top-N PAS features from treatment and control.

#

# 4. Within each dataset, fit DESeq2 with design = ~ condition.

#

# 5. Use finite DESeq2 Wald statistics for empirical-null recalibration with

#    fdrtool.

#

# 6. Compute the higher-criticism threshold from sorted empirical-null

#    p-values using hc.thresh(sort(empirical_p)).

#

# 7. Define HBFSS as:

#        HBFSS = abs(lfc_shrunk * -log10(empirical_p))

#

#    Define the dataset-specific HBFSS threshold as:

#        hbfss_threshold = abs(log10(hc_p_threshold)) * lfc_boundary

#

# 8. Volcano plotting rule:

#    HC is the hard floor for all plotted classes.

#    Nothing below HC may be colored Weak CNH, Strong CNH, Standard, or HBFSS.

#

# 9. The only plotted volcano classes are:

#        Background

#        Weak CNH

#        Strong CNH

#        Standard

#        HBFSS

#

#    Overlap is counted and written in figure text, but overlap is not a

#    separate displayed plotting class.

#

# 10. Final mutually exclusive plotting precedence is:

#        HBFSS > Standard > Strong CNH > Weak CNH > Background

#

# 11. Compare figures use exactly three panels:

#        Raw / Lead / Rem

#    with one shared legend only.

#

# 12. Export names are kept short and GitHub-sortable:

#        Figure_<comparison>_<tag>.png

#        Table_<comparison>_<tag>.csv

#

# Notes

# -----

# - This script is written to be runnable as a standalone manuscript script.

# - If the required packages are installed and the WTTS count file exists in one

#   of the candidate locations below, the script should run end to end.

# - If a comparison fails, the failure is captured in Table_Failed.csv.

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

  library(gtable)

})

options(stringsAsFactors = FALSE)

# =============================================================================

# USER SETTINGS

# =============================================================================

count_file_candidates <- c(

  "WTTS-Seq_2022.2_DE_raw_read_numbers.csv",

  file.path("data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"),

  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

)

run_twas_overlap <- FALSE

twas_file <- "3aTWAS_genes_of_11_brain_disorders.csv"

alpha_level <- 0.10

lfc_boundary <- 1.0

top_n_target <- 5000L

figure_dpi <- 320

base_theme_size <- 9

single_volcano_label_n <- 14

compare_volcano_label_n <- 8

output_dir <- file.path("exports", "manuscript_final_clean")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================

# PLOTTING PALETTE

# =============================================================================

plot_palette <- list(

  background = "#BDBDBD",

  weak_cnh   = "#2C7FB8",

  strong_cnh = "#D95F02",

  standard   = "#33A02C",

  hbfss      = "#6A3D9A",

  hc_line    = "#A65628",

  lfc_line   = "#A65628",

  control    = "#4D4D4D",

  treatment  = "#1F78B4"

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

condition_shapes <- c("untrt" = 21, "trt" = 24)

condition_fills  <- c("untrt" = plot_palette$control, "trt" = plot_palette$treatment)

condition_labels <- c("untrt" = "Control", "trt" = "Treatment")

dataset_short <- c(

  raw_dataset = "Raw",

  leading_edge_dataset = "Lead",

  remainder_dataset = "Rem"

)

dataset_pretty <- c(

  raw_dataset = "Raw",

  leading_edge_dataset = "Lead",

  remainder_dataset = "Rem"

)

# =============================================================================

# METADATA

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

comparison_table <- data.frame(

  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),

  group1_prefix   = c("R0", "R2", "R4", "R8"),

  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),

  stringsAsFactors = FALSE

)

# =============================================================================

# GENERAL HELPERS

# =============================================================================

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

find_existing_file <- function(candidates) {

  for (p in candidates) {

    if (file.exists(p)) return(normalizePath(p))

  }

  stop(

    "Could not find WTTS count file. Checked: ",

    paste(candidates, collapse = " | ")

  )

}

safe_log10 <- function(x, pseudocount = 1e-300) {

  log10(pmax(as.numeric(x), pseudocount))

}

safe_neglog10 <- function(x, pseudocount = 1e-300) {

  -log10(pmax(as.numeric(x), pseudocount))

}

clip_probabilities <- function(x, eps = 1e-300) {

  x <- as.numeric(x)

  bad <- !is.finite(x) | is.na(x)

  x[bad] <- NA_real_

  good <- !is.na(x)

  x[good] <- pmin(pmax(x[good], eps), 1 - 1e-12)

  x

}

compact_title <- function(x, width = 48) {

  paste(strwrap(as.character(x), width = width), collapse = "\n")

}

manuscript_theme <- function() {

  theme_bw(base_size = base_theme_size) +

    theme(

      plot.title       = element_text(face = "bold", hjust = 0.5, size = base_theme_size + 1),

      plot.subtitle    = element_text(hjust = 0.5, size = base_theme_size - 1),

      plot.caption     = element_text(size = base_theme_size - 2, colour = "grey25"),

      axis.title       = element_text(face = "bold"),

      axis.text        = element_text(colour = "black"),

      legend.position  = "bottom",

      legend.title     = element_text(face = "bold"),

      panel.grid.minor = element_blank(),

      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),

      plot.margin      = margin(t = 8, r = 8, b = 8, l = 8)

    )

}

save_plot <- function(p, path, width = 7.0, height = 5.2, dpi = figure_dpi) {

  ggsave(

    filename = path,

    plot = p,

    width = width,

    height = height,

    dpi = dpi,

    units = "in",

    bg = "white",

    limitsize = FALSE

  )

}

save_grob <- function(g, path, width = 12.5, height = 5.8, dpi = figure_dpi) {

  ggsave(

    filename = path,

    plot = g,

    width = width,

    height = height,

    dpi = dpi,

    units = "in",

    bg = "white",

    limitsize = FALSE

  )

}

save_csv <- function(df, path) {

  write.csv(df, path, row.names = FALSE)

}

fig_file <- function(dir, cmp, tag) {

  file.path(dir, paste0("Figure_", cmp, "_", tag, ".png"))

}

tab_file <- function(dir, cmp, tag) {

  file.path(dir, paste0("Table_", cmp, "_", tag, ".csv"))

}

get_condition_coef <- function(dds) {

  rn <- resultsNames(dds)

  idx <- grep("^condition_", rn)

  if (length(idx) == 0) stop("Could not identify condition coefficient in resultsNames(dds).")

  rn[idx[1]]

}

get_legend_grob <- function(plot_obj) {

  g <- ggplotGrob(plot_obj)

  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")

  if (length(idx) == 0) return(NULL)

  g$grobs[[idx[1]]]

}

# =============================================================================

# IMPORT WTTS MASTER COUNT FILE

# =============================================================================

count_file <- find_existing_file(count_file_candidates)

WTTS_Seq <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

WTTS_Seq <- as.data.frame(WTTS_Seq, stringsAsFactors = FALSE)

assert_required_columns(WTTS_Seq, c("OrigID", "Symbol"), object_name = "WTTS count file")

assert_required_columns(WTTS_Seq, meta_all$id, object_name = "WTTS count file sample columns")

WTTS_Seq$OrigID <- as.character(WTTS_Seq$OrigID)

WTTS_Seq$Symbol <- as.character(WTTS_Seq$Symbol)

WTTS_Seq <- WTTS_Seq[

  !is.na(WTTS_Seq$OrigID) &

  !is.na(WTTS_Seq$Symbol),

  ,

  drop = FALSE

]

sample_na <- rowSums(is.na(WTTS_Seq[, meta_all$id, drop = FALSE])) > 0

WTTS_Seq <- WTTS_Seq[!sample_na, , drop = FALSE]

rownames(WTTS_Seq) <- WTTS_Seq$OrigID

OrigID_Symbol <- unique(WTTS_Seq[, c("OrigID", "Symbol"), drop = FALSE])

colnames(OrigID_Symbol) <- c("feature_id", "gene_symbol")

OrigID_Symbol$feature_id <- as.character(OrigID_Symbol$feature_id)

OrigID_Symbol$gene_symbol <- as.character(OrigID_Symbol$gene_symbol)

OrigID_Symbol <- OrigID_Symbol %>%

  mutate(gene_symbol = if_else(is.na(gene_symbol), "", trimws(gene_symbol))) %>%

  arrange(feature_id, desc(gene_symbol != ""), gene_symbol) %>%

  distinct(feature_id, .keep_all = TRUE) %>%

  mutate(gene_symbol = na_if(gene_symbol, ""))

if (run_twas_overlap) {

  TWAS_Seq <- read.csv(twas_file, header = TRUE, stringsAsFactors = FALSE)

  TWAS_Seq <- as.data.frame(TWAS_Seq)

  if (ncol(TWAS_Seq) < 4) stop("TWAS file must contain at least 4 columns.")

  TWAS_data <- TWAS_Seq[, c(1, 4), drop = FALSE]

  colnames(TWAS_data) <- c("source_id", "gene_symbol")

}

# =============================================================================

# COMPARISON PREPARATION

# =============================================================================

prepare_comparison_data <- function(comparison_name, group1_prefix, group2_prefix, WTTS_Seq, meta_all) {

  keep_ids <- grepl(paste0("^", group1_prefix, "_"), meta_all$id) |

    grepl(paste0("^", group2_prefix, "_"), meta_all$id)

  meta_sub <- meta_all[keep_ids, , drop = FALSE]

  coldata <- meta_sub[, "condition", drop = FALSE]

  sample_ids <- rownames(meta_sub)

  missing_samples <- setdiff(sample_ids, colnames(WTTS_Seq))

  if (length(missing_samples) > 0) {

    stop(

      paste0(

        "Missing samples in WTTS count file for ", comparison_name, ": ",

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

# =============================================================================

# EVS

# =============================================================================

resolve_top_n_cutoff <- function(sorted_values_desc, top_n = top_n_target) {

  n_total <- length(sorted_values_desc)

  if (n_total == 0) stop("resolve_top_n_cutoff() received an empty vector.")

  top_n_actual <- min(max(1L, as.integer(top_n)), n_total)

  cutoff_value <- sorted_values_desc[top_n_actual]

  cutoff_quantile <- 1 - (top_n_actual / n_total)

  list(

    top_n_actual = top_n_actual,

    cutoff_value = cutoff_value,

    cutoff_quantile = cutoff_quantile

  )

}

compute_pc1_loading_table <- function(value_df, sample_names, top_n = top_n_target) {

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

    cutoff_quantile = cutoff_info$cutoff_quantile

  )

}

build_eigenvector_split <- function(count_matrix, coldata) {

  dds_init <- DESeqDataSetFromMatrix(

    countData = round(count_matrix),

    colData = coldata,

    design = ~ condition

  )

  dds_init <- dds_init[rowSums(counts(dds_init)) > 0, ]

  dds_init <- estimateSizeFactors(dds_init)

  norm_counts_init <- as.data.frame(counts(dds_init, normalized = TRUE))

  raw_counts_init <- as.data.frame(count_matrix)

  sample_ids <- colnames(count_matrix)

  trt_ids <- sample_ids[coldata$condition == "trt"]

  untrt_ids <- sample_ids[coldata$condition == "untrt"]

  fit_trt <- compute_pc1_loading_table(norm_counts_init, trt_ids, top_n = top_n_target)

  fit_untrt <- compute_pc1_loading_table(norm_counts_init, untrt_ids, top_n = top_n_target)

  trt_high <- as.character(subset(fit_trt$loading_table, split_class == "high_loading")$feature_id)

  untrt_high <- as.character(subset(fit_untrt$loading_table, split_class == "high_loading")$feature_id)

  leading_edge_ids <- union(trt_high, untrt_high)

  remainder_ids <- setdiff(rownames(count_matrix), leading_edge_ids)

  if (length(leading_edge_ids) == 0) {

    stop("Leading-edge dataset is empty.")

  }

  list(

    fit_trt = fit_trt,

    fit_untrt = fit_untrt,

    normalized_counts = norm_counts_init,

    raw_counts = raw_counts_init,

    raw_dataset = count_matrix,

    leading_edge_dataset = count_matrix[leading_edge_ids, , drop = FALSE],

    remainder_dataset = count_matrix[remainder_ids, , drop = FALSE]

  )

}

# =============================================================================

# PCA SUPPORT FIGURE

# =============================================================================

plot_pca_scatter <- function(pca_fit, comparison_name, group_label) {

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

      title = compact_title(paste0(comparison_name, " ", ifelse(cond_key == "trt", "treatment", "control"), " PCA")),

      subtitle = paste0("PC1=", pca_var_per[1], "%  PC2=", pca_var_per[2], "%"),

      x = "PC1",

      y = "PC2"

    ) +

    manuscript_theme()

}

# =============================================================================

# EMPIRICAL-NULL AND HC

# =============================================================================

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

      tryCatch(

        fdrtool(

          stat_vec,

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

# =============================================================================

# DESEQ2 CORE ANALYSIS

# =============================================================================

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {

  dds <- DESeqDataSetFromMatrix(

    countData = round(count_mat),

    colData = coldata,

    design = ~ condition

  )

  dds <- dds[rowSums(counts(dds)) > 0, ]

  dds <- DESeq(dds, betaPrior = FALSE, quiet = TRUE)

  res <- results(dds, contrast = c("condition", "trt", "untrt"), alpha = alpha_level)

  res_strong <- results(

    dds,

    contrast = c("condition", "trt", "untrt"),

    lfcThreshold = lfc_boundary,

    altHypothesis = "greaterAbs"

  )

  res_weak <- results(

    dds,

    contrast = c("condition", "trt", "untrt"),

    lfcThreshold = lfc_boundary,

    altHypothesis = "lessAbs"

  )

  res_all_df <- as.data.frame(res)

  res_all_df$feature_id <- as.character(rownames(res_all_df))

  valid_stat <- is.finite(res_all_df$stat) & !is.na(res_all_df$stat)

  stat_vec <- as.numeric(res_all_df$stat[valid_stat])

  if (length(stat_vec) < 5) {

    stop(sprintf("[%s] Fewer than 5 finite Wald statistics were available for fdrtool.", dataset_name))

  }

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

  coef_name <- get_condition_coef(dds)

  shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)

  shr_df <- as.data.frame(shr)

  shr_df$feature_id <- as.character(rownames(shr_df))

  res_df <- left_join(

    res_df,

    shr_df[, c("feature_id", "log2FoldChange")],

    by = "feature_id",

    suffix = c("", "_shrunk")

  )

  colnames(res_df)[colnames(res_df) == "log2FoldChange_shrunk"] <- "lfc_shrunk"

  res_strong_df <- as.data.frame(res_strong)

  res_strong_df$feature_id <- as.character(rownames(res_strong_df))

  res_weak_df <- as.data.frame(res_weak)

  res_weak_df$feature_id <- as.character(rownames(res_weak_df))

  res_df <- left_join(res_df, res_strong_df[, c("feature_id", "padj")], by = "feature_id", suffix = c("", "_strong"))

  res_df <- left_join(res_df, res_weak_df[, c("feature_id", "padj")], by = "feature_id", suffix = c("", "_weak"))

  colnames(res_df)[colnames(res_df) == "padj_strong"] <- "resGA_padj"

  colnames(res_df)[colnames(res_df) == "padj_weak"] <- "resLA_padj"

  hc_p_threshold_dataset <- safe_hc_thresh(res_df$empirical_p, dataset_name = dataset_name)

  hbfss_threshold_dataset <- if (

    is.finite(hc_p_threshold_dataset) &&

    !is.na(hc_p_threshold_dataset) &&

    hc_p_threshold_dataset > 0 &&

    hc_p_threshold_dataset < 1

  ) {

    abs(log10(hc_p_threshold_dataset)) * lfc_boundary

  } else {

    NA_real_

  }

  res_df$HBFSS <- abs(res_df$lfc_shrunk * safe_neglog10(res_df$empirical_p))

  # ---------------------------------------------------------------------------

  # FINAL MANUSCRIPT LOGIC

  #

  # HC is the hard global coloring floor.

  # A point must satisfy empirical_p <= hc_p_threshold_dataset before it can

  # be plotted as Weak CNH, Strong CNH, Standard, or HBFSS.

  #

  # Weak CNH:

  #   - pass HC

  #   - lessAbs adjusted p-value < alpha

  #   - |lfc_shrunk| < lfc_boundary

  #

  # Strong CNH:

  #   - pass HC

  #   - greaterAbs adjusted p-value < alpha

  #   - |lfc_shrunk| >= lfc_boundary

  #

  # Standard:

  #   - pass HC

  #   - native DESeq2 padj < alpha

  #   - |lfc_shrunk| >= lfc_boundary

  #

  # HBFSS:

  #   - pass HC

  #   - HBFSS >= hbfss_threshold_dataset

  #

  # Final display precedence:

  #   HBFSS > Standard > Strong CNH > Weak CNH > Background

  # ---------------------------------------------------------------------------

  res_df$hc_pass <- is.finite(hc_p_threshold_dataset) &

    !is.na(res_df$empirical_p) &

    (res_df$empirical_p <= hc_p_threshold_dataset)

  res_df$weak_cnh_pass <- res_df$hc_pass &

    !is.na(res_df$resLA_padj) &

    (res_df$resLA_padj < alpha_level) &

    !is.na(res_df$lfc_shrunk) &

    (abs(res_df$lfc_shrunk) < lfc_boundary)

  res_df$strong_cnh_pass <- res_df$hc_pass &

    !is.na(res_df$resGA_padj) &

    (res_df$resGA_padj < alpha_level) &

    !is.na(res_df$lfc_shrunk) &

    (abs(res_df$lfc_shrunk) >= lfc_boundary)

  res_df$standard_pass <- res_df$hc_pass &

    !is.na(res_df$padj) &

    (res_df$padj < alpha_level) &

    !is.na(res_df$lfc_shrunk) &

    (abs(res_df$lfc_shrunk) >= lfc_boundary)

  res_df$hbfss_pass <- res_df$hc_pass &

    is.finite(hbfss_threshold_dataset) &

    !is.na(res_df$HBFSS) &

    (res_df$HBFSS >= hbfss_threshold_dataset)

  res_df$plot_class <- "Background"

  res_df$plot_class[res_df$weak_cnh_pass] <- "Weak CNH"

  res_df$plot_class[res_df$strong_cnh_pass] <- "Strong CNH"

  res_df$plot_class[res_df$standard_pass] <- "Standard"

  res_df$plot_class[res_df$hbfss_pass] <- "HBFSS"

  res_df$plot_class <- factor(res_df$plot_class, levels = plot_class_levels)

  res_df$overlap_standard_hbfss <- res_df$standard_pass & res_df$hbfss_pass

  res_df$feature_id <- as.character(res_df$feature_id)

  annot_df$feature_id <- as.character(annot_df$feature_id)

  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))

  norm_counts$feature_id <- as.character(rownames(norm_counts))

  mm <- as.data.frame(mcols(dds))

  mm$feature_id <- as.character(rownames(mm))

  disp_cols_available <- intersect(

    c("feature_id", "dispGeneEst", "dispFit", "dispersion", "dispIter", "baseMean", "dispOutlier"),

    colnames(mm)

  )

  disp_df <- mm[, disp_cols_available, drop = FALSE]

  if ("baseMean" %in% colnames(disp_df) && "dispersion" %in% colnames(disp_df)) {

    disp_df$IOD <- disp_df$baseMean * disp_df$dispersion

  }

  final_df <- res_df %>%

    left_join(annot_df, by = "feature_id") %>%

    left_join(norm_counts, by = "feature_id") %>%

    left_join(disp_df, by = "feature_id")

  final_df$neglog10_empirical_p <- safe_neglog10(final_df$empirical_p)

  final_df$neglog10_padj <- safe_neglog10(final_df$padj)

  final_df$dataset_name <- dataset_name

  final_df$hc_p_threshold_dataset <- hc_p_threshold_dataset

  final_df$hbfss_threshold_dataset <- hbfss_threshold_dataset

  final_df$label_ok <- !is.na(final_df$gene_symbol) & nzchar(trimws(final_df$gene_symbol))

  list(

    dds = dds,

    results = final_df,

    hc_p_threshold = hc_p_threshold_dataset,

    hbfss_threshold = hbfss_threshold_dataset

  )

}

# =============================================================================

# VOLCANO FIGURE HELPERS

# =============================================================================

select_volcano_labels <- function(df, n_labels = 12) {

  df <- as.data.frame(df)

  df <- df[df$label_ok, , drop = FALSE]

  if (!nrow(df)) return(df[0, , drop = FALSE])

  df$label_priority <- ifelse(df$plot_class == "HBFSS", 1,

                       ifelse(df$plot_class == "Standard", 2,

                       ifelse(df$plot_class == "Strong CNH", 3,

                       ifelse(df$plot_class == "Weak CNH", 4, 5))))

  ord <- order(df$label_priority, -df$neglog10_empirical_p, -abs(df$lfc_shrunk), na.last = TRUE)

  df <- df[ord, , drop = FALSE]

  df <- df[!duplicated(df$gene_symbol), , drop = FALSE]

  df[seq_len(min(n_labels, nrow(df))), , drop = FALSE]

}

make_hbfss_curve_df <- function(xlim_max, hbfss_threshold) {

  if (!is.finite(hbfss_threshold) || is.na(hbfss_threshold) || hbfss_threshold <= 0) return(NULL)

  x_left <- seq(-xlim_max, -0.08, length.out = 400)

  x_right <- seq(0.08, xlim_max, length.out = 400)

  x <- c(x_left, x_right)

  y <- hbfss_threshold / abs(x)

  data.frame(x = x, y = y)

}

volcano_counts_text <- function(df) {

  paste0(

    "Weak=", sum(df$weak_cnh_pass, na.rm = TRUE),

    "  Strong=", sum(df$strong_cnh_pass, na.rm = TRUE),

    "  Std=", sum(df$standard_pass, na.rm = TRUE),

    "  HBFSS=", sum(df$hbfss_pass, na.rm = TRUE),

    "  Overlap=", sum(df$overlap_standard_hbfss, na.rm = TRUE)

  )

}

build_single_volcano <- function(df, cmp_name, ds_key, show_legend = TRUE, label_n = single_volcano_label_n) {

  ds_name <- dataset_pretty[[ds_key]]

  hc_thr <- unique(df$hc_p_threshold_dataset)[1]

  hbfss_thr <- unique(df$hbfss_threshold_dataset)[1]

  finite_lfc <- df$lfc_shrunk[is.finite(df$lfc_shrunk)]

  finite_y <- df$neglog10_empirical_p[is.finite(df$neglog10_empirical_p)]

  xlim_max <- if (length(finite_lfc)) max(3, ceiling(max(abs(finite_lfc)))) else 3

  ylim_max <- if (length(finite_y)) max(4, ceiling(max(finite_y, na.rm = TRUE))) else 4

  if (is.finite(hbfss_thr) && !is.na(hbfss_thr)) ylim_max <- max(ylim_max, ceiling(hbfss_thr + 2))

  curve_df <- make_hbfss_curve_df(xlim_max, hbfss_thr)

  lab_df <- select_volcano_labels(df[df$plot_class != "Background", , drop = FALSE], n_labels = label_n)

  p <- ggplot(df, aes(lfc_shrunk, neglog10_empirical_p)) +

    geom_point(

      aes(color = plot_class, shape = plot_class),

      size = 1.15,

      alpha = 0.85,

      stroke = 0.25

    ) +

    geom_vline(xintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", linewidth = 0.45, colour = plot_palette$lfc_line) +

    geom_vline(xintercept = 0, linewidth = 0.30, colour = "grey55") +

    scale_color_manual(values = plot_class_colors, drop = FALSE, name = "Class") +

    scale_shape_manual(values = plot_class_shapes, drop = FALSE, name = "Class") +

    labs(

      title = ds_name,

      subtitle = "x = shrunken log2FC   ·   y = -log10(empirical p)",

      x = "Shrunken log2FC",

      y = expression(-log[10]("Empirical p"))

    ) +

    coord_cartesian(xlim = c(-xlim_max, xlim_max), ylim = c(0, ylim_max), clip = "off") +

    manuscript_theme() +

    theme(

      legend.position = if (show_legend) "bottom" else "none",

      plot.title = element_text(size = base_theme_size + 1),

      plot.subtitle = element_text(size = base_theme_size - 2),

      axis.title = element_text(size = base_theme_size),

      axis.text = element_text(size = base_theme_size - 1)

    )

  if (is.finite(hc_thr) && !is.na(hc_thr) && hc_thr > 0 && hc_thr < 1) {

    hc_y <- safe_neglog10(hc_thr)

    p <- p +

      geom_hline(yintercept = hc_y, linetype = "dotted", linewidth = 0.50, colour = plot_palette$hc_line) +

      annotate(

        "text",

        x = -xlim_max * 0.82,

        y = min(ylim_max - 0.3, hc_y + 0.18),

        label = paste0("HC=", signif(hc_thr, 3)),

        colour = plot_palette$hc_line,

        size = 2.8,

        hjust = 0

      )

  }

  if (!is.null(curve_df)) {

    p <- p +

      geom_line(data = curve_df, aes(x, y), inherit.aes = FALSE, colour = plot_palette$hbfss, linewidth = 0.75) +

      annotate(

        "text",

        x = xlim_max * 0.58,

        y = min(ylim_max - 0.3, max(0.8, hbfss_thr + 0.25)),

        label = paste0("HBFSS=", signif(hbfss_thr, 3)),

        colour = plot_palette$hbfss,

        size = 2.8,

        hjust = 0

      )

  }

  p <- p +

    annotate(

      "text",

      x = 0,

      y = ylim_max * 0.97,

      label = volcano_counts_text(df),

      size = 2.8,

      hjust = 0.5,

      vjust = 1

    )

  if (nrow(lab_df) > 0) {

    p <- p +

      ggrepel::geom_text_repel(

        data = lab_df,

        aes(label = gene_symbol),

        size = 2.0,

        seed = 1,

        max.overlaps = 20,

        box.padding = 0.22,

        point.padding = 0.10,

        min.segment.length = 0

      )

  }

  p

}

build_compare_volcano <- function(res_list, cmp_name) {

  p_raw <- build_single_volcano(res_list$raw_dataset, cmp_name, "raw_dataset", show_legend = TRUE, label_n = compare_volcano_label_n)

  p_lead <- build_single_volcano(res_list$leading_edge_dataset, cmp_name, "leading_edge_dataset", show_legend = FALSE, label_n = compare_volcano_label_n)

  p_rem <- build_single_volcano(res_list$remainder_dataset, cmp_name, "remainder_dataset", show_legend = FALSE, label_n = compare_volcano_label_n)

  legend_grob <- get_legend_grob(p_raw)

  p_raw <- p_raw + theme(legend.position = "none")

  row_grob <- arrangeGrob(

    grobs = list(p_raw, p_lead, p_rem),

    ncol = 3,

    top = textGrob(

      paste0(cmp_name, " volcano"),

      gp = gpar(fontface = "bold", cex = 1.05)

    )

  )

  if (is.null(legend_grob)) {

    arrangeGrob(row_grob, ncol = 1)

  } else {

    arrangeGrob(

      row_grob,

      legend_grob,

      ncol = 1,

      heights = unit.c(unit(1, "npc") - unit(0.55, "in"), unit(0.55, "in"))

    )

  }

}

# =============================================================================

# DISPERSION FIGURES

# =============================================================================

build_single_dispersion <- function(df, ds_key, show_legend = TRUE) {

  ds_name <- dataset_pretty[[ds_key]]

  p <- ggplot(df, aes(baseMean, dispersion)) +

    geom_point(aes(color = plot_class, shape = plot_class), size = 0.8, alpha = 0.70, stroke = 0.20) +

    scale_x_log10(labels = label_number(accuracy = 0.1)) +

    scale_y_log10(labels = label_number(accuracy = 0.01)) +

    scale_color_manual(values = plot_class_colors, drop = FALSE, name = "Class") +

    scale_shape_manual(values = plot_class_shapes, drop = FALSE, name = "Class") +

    labs(

      title = ds_name,

      subtitle = "Final dispersion vs mean",

      x = "baseMean",

      y = "Dispersion"

    ) +

    manuscript_theme() +

    theme(

      legend.position = if (show_legend) "bottom" else "none",

      plot.title = element_text(size = base_theme_size + 1),

      plot.subtitle = element_text(size = base_theme_size - 2),

      axis.title = element_text(size = base_theme_size),

      axis.text = element_text(size = base_theme_size - 1)

    )

  p

}

build_compare_dispersion <- function(res_list, cmp_name) {

  p_raw <- build_single_dispersion(res_list$raw_dataset, "raw_dataset", show_legend = TRUE)

  p_lead <- build_single_dispersion(res_list$leading_edge_dataset, "leading_edge_dataset", show_legend = FALSE)

  p_rem <- build_single_dispersion(res_list$remainder_dataset, "remainder_dataset", show_legend = FALSE)

  legend_grob <- get_legend_grob(p_raw)

  p_raw <- p_raw + theme(legend.position = "none")

  row_grob <- arrangeGrob(

    grobs = list(p_raw, p_lead, p_rem),

    ncol = 3,

    top = textGrob(

      paste0(cmp_name, " dispersion"),

      gp = gpar(fontface = "bold", cex = 1.05)

    )

  )

  if (is.null(legend_grob)) {

    arrangeGrob(row_grob, ncol = 1)

  } else {

    arrangeGrob(

      row_grob,

      legend_grob,

      ncol = 1,

      heights = unit.c(unit(1, "npc") - unit(0.55, "in"), unit(0.55, "in"))

    )

  }

}

# =============================================================================

# RUN ONE FULL COMPARISON

# =============================================================================

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {

  cmp_dir <- file.path(output_dir, comparison_name)

  dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

  evs <- build_eigenvector_split(count_matrix, coldata)

  save_plot(

    plot_pca_scatter(evs$fit_trt$pca_fit, comparison_name, "treatment"),

    fig_file(cmp_dir, comparison_name, "EVS_Trt_PCA"),

    width = 7.2, height = 5.0

  )

  save_plot(

    plot_pca_scatter(evs$fit_untrt$pca_fit, comparison_name, "control"),

    fig_file(cmp_dir, comparison_name, "EVS_Ctrl_PCA"),

    width = 7.2, height = 5.0

  )

  dataset_list <- list(

    raw_dataset = evs$raw_dataset,

    leading_edge_dataset = evs$leading_edge_dataset,

    remainder_dataset = evs$remainder_dataset

  )

  analysis_results <- list()

  summary_rows <- list()

  for (ds_key in names(dataset_list)) {

    res <- run_core_analysis(

      count_mat = dataset_list[[ds_key]],

      coldata = coldata,

      dataset_name = paste0(comparison_name, "_", ds_key),

      annot_df = annot_df

    )

    df <- res$results

    analysis_results[[ds_key]] <- df

    ds_short <- dataset_short[[ds_key]]

    save_csv(df, tab_file(cmp_dir, comparison_name, paste0(ds_short, "_Results")))

    save_csv(

      df %>%

        select(

          feature_id, gene_symbol, lfc_shrunk, empirical_p, padj, HBFSS,

          hc_pass, weak_cnh_pass, strong_cnh_pass, standard_pass, hbfss_pass,

          plot_class, overlap_standard_hbfss

        ),

      tab_file(cmp_dir, comparison_name, paste0(ds_short, "_PlotLogic"))

    )

    save_plot(

      build_single_volcano(df, comparison_name, ds_key, show_legend = TRUE, label_n = single_volcano_label_n),

      fig_file(cmp_dir, comparison_name, paste0(ds_short, "_Volcano")),

      width = 7.6, height = 6.2

    )

    summary_rows[[ds_key]] <- data.frame(

      comparison_name = comparison_name,

      dataset_key = ds_key,

      dataset_label = dataset_pretty[[ds_key]],

      n_features = nrow(df),

      hc = res$hc_p_threshold,

      hbfss = res$hbfss_threshold,

      weak_n = sum(df$weak_cnh_pass, na.rm = TRUE),

      strong_n = sum(df$strong_cnh_pass, na.rm = TRUE),

      standard_n = sum(df$standard_pass, na.rm = TRUE),

      hbfss_n = sum(df$hbfss_pass, na.rm = TRUE),

      overlap_n = sum(df$overlap_standard_hbfss, na.rm = TRUE),

      stringsAsFactors = FALSE

    )

  }

  compare_volcano_grob <- build_compare_volcano(analysis_results, comparison_name)

  save_grob(

    compare_volcano_grob,

    fig_file(cmp_dir, comparison_name, "Compare_Volcano"),

    width = 13.0,

    height = 5.6

  )

  compare_disp_grob <- build_compare_dispersion(analysis_results, comparison_name)

  save_grob(

    compare_disp_grob,

    fig_file(cmp_dir, comparison_name, "Compare_Disp"),

    width = 13.0,

    height = 5.6

  )

  comparison_summary <- bind_rows(summary_rows)

  save_csv(comparison_summary, tab_file(cmp_dir, comparison_name, "Summary"))

  comparison_summary

}

# =============================================================================

# BATCH EXECUTION

# =============================================================================

message("Using count file: ", count_file)

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

all_summaries <- list()

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

    all_summaries[[cmp]] <- out

  }

}

if (length(all_summaries) > 0) {

  save_csv(bind_rows(all_summaries), file.path(output_dir, "Table_Overall_Summary.csv"))

}

if (length(failed_comparisons) > 0) {

  save_csv(bind_rows(failed_comparisons), file.path(output_dir, "Table_Failed.csv"))

}

cat("\n=====================================================\n")

cat("Pipeline complete.\n")

cat("Output directory:\n")

cat(normalizePath(output_dir), "\n")

cat("=====================================================\n\n")
