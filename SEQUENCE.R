Save this as SEQUENCE.R.

#!/usr/bin/env Rscript

# =============================================================================
# MANUSCRIPT PIPELINE
# EVS + DESeq2 + empirical-null calibration + CNH/HBFSS volcano panels
# -----------------------------------------------------------------------------
# Purpose
# This script runs a complete manuscript-ready analysis for the WTTS PAS matrix.
# It performs four pairwise circadian comparisons, builds an eigenvector split
# (EVS), runs DESeq2 on the original / leading-edge / remainder datasets,
# calibrates Wald statistics with an empirical-null model via fdrtool, and
# exports concise manuscript figures and tables directly into the Git-tracked
# repository tree under exports/manuscript_final_clean.
#
# Core manuscript logic
# 1. Ranking / EVS stage
#    - DESeq2 size-factor normalization is used before PCA-based ranking.
#    - Absolute PC1 loadings are computed separately within treatment and
#      control.
#    - The top_n_target PAS sites from either condition define the union-based
#      leading edge.
#    - The remainder contains PAS sites that are not high-loading in either arm.
#
# 2. Differential-expression stage
#    - Each dataset (Raw / Lead / Rem) is fit as its own DESeq2 object with
#      design = ~ condition.
#    - Standard DESeq2, weak composite-null hypothesis (lessAbs), and strong
#      composite-null hypothesis (greaterAbs) results are all retained.
#
# 3. Empirical-null / HBFSS stage
#    - Only finite DESeq2 Wald statistics are sent to fdrtool.
#    - fdrtool returns empirical-null p-values, q-values, and local FDR.
#    - Higher criticism is applied to the sorted empirical-null p-values.
#    - The dataset-specific HBFSS boundary is
#          abs(log10(HC_p_threshold)) * lfc_boundary
#    - Gene-specific HBFSS is
#          abs(shrunken_log2FC * log10(empirical_p))
#
# 4. Volcano figure stage
#    - Each volcano panel plots exactly these four colored classes:
#         Weak CNH
#         Strong CNH
#         Standard
#         HBFSS
#      Background is shown in grey.
#    - Standard/HBFSS overlap counts are written in panel text rather than as
#      separate legend classes.
#    - Each 3-panel comparison figure has one shared legend only.
#    - Each panel shows:
#         * ±lfc_boundary lines
#         * HC threshold line on the y-axis
#         * HBFSS boundary curve
#         * panel-specific HC and HBFSS values in short text
#
# Export policy
# - All exported filenames are short and begin with Figure_ or Table_.
# - All outputs are written under exports/manuscript_final_clean.
# - No Rplots.pdf is allowed to accumulate.
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
  library(gtable)
})

options(stringsAsFactors = FALSE)
unlink("Rplots.pdf", force = TRUE)

# =============================================================================
# USER SETTINGS
# =============================================================================

count_file_candidates <- c(
  "data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
)

count_file <- count_file_candidates[file.exists(count_file_candidates)][1]
if (is.na(count_file) || !nzchar(count_file)) {
  stop("Could not find WTTS-Seq_2022.2_DE_raw_read_numbers.csv in repo/data or working directory.")
}

output_dir <- file.path("exports", "manuscript_final_clean")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

alpha_level     <- 0.10
lfc_boundary    <- 1.0
top_n_target    <- 5000L
figure_dpi      <- 320
base_theme_size <- 9
n_top_labels    <- 14

plot_palette <- list(
  background = "#BDBDBD",
  weak_cnh   = "#1F78B4",
  strong_cnh = "#D95F02",
  standard   = "#33A02C",
  hbfss      = "#6A3D9A",
  threshold  = "#8C2D04",
  control    = "#4D4D4D",
  treatment  = "#1F78B4",
  hist       = "#969696"
)

condition_shapes <- c("untrt" = 21, "trt" = 24)
condition_fills  <- c("untrt" = plot_palette$control, "trt" = plot_palette$treatment)
condition_labels <- c("untrt" = "Control", "trt" = "Treatment")

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix   = c("R0", "R2", "R4", "R8"),
  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),
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

# =============================================================================
# FILE-NAME HELPERS
# =============================================================================

dataset_short <- c(
  raw_dataset = "Raw",
  leading_edge_dataset = "Lead",
  remainder_dataset = "Rem"
)

fig_file <- function(dir, cmp, tag, ds_key = NULL) {
  nm <- if (is.null(ds_key)) {
    paste0("Figure_", cmp, "_", tag, ".png")
  } else {
    paste0("Figure_", cmp, "_", unname(dataset_short[ds_key]), "_", tag, ".png")
  }
  file.path(dir, nm)
}

tab_file <- function(dir, cmp, tag, ds_key = NULL) {
  nm <- if (is.null(ds_key)) {
    paste0("Table_", cmp, "_", tag, ".csv")
  } else {
    paste0("Table_", cmp, "_", unname(dataset_short[ds_key]), "_", tag, ".csv")
  }
  file.path(dir, nm)
}

# =============================================================================
# GENERAL HELPERS
# =============================================================================

assert_required_columns <- function(df, required_cols, object_name = "data frame") {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in %s: %s", object_name, paste(missing_cols, collapse = ", ")))
  }
}

safe_log10 <- function(x, pseudocount = 1e-12) log10(pmax(x, pseudocount))
safe_neglog10 <- function(x, pseudocount = 1e-12) -log10(pmax(x, pseudocount))

clip_probabilities <- function(x, eps = 1e-300) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  ok <- !is.na(x)
  x[ok] <- pmin(pmax(x[ok], eps), 1 - 1e-12)
  x
}

compact_title <- function(x, width = 42) paste(strwrap(as.character(x), width = width), collapse = "\n")
compact_caption <- function(x, width = 110) paste(strwrap(as.character(x), width = width), collapse = "\n")

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = base_theme_size + 1),
      plot.subtitle   = element_text(hjust = 0.5, size = base_theme_size - 1),
      plot.caption    = element_text(hjust = 0.5, size = base_theme_size - 2, colour = "grey30"),
      axis.title      = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      legend.position = "bottom",
      legend.title    = element_text(face = "bold"),
      legend.text     = element_text(size = base_theme_size - 1),
      plot.margin     = margin(10, 14, 12, 14)
    )
}

save_plot <- function(p, path, width = 10, height = 7, dpi = figure_dpi) {
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

save_grob <- function(g, path, width = 12, height = 8, dpi = figure_dpi) {
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
  write.csv(df, file = path, row.names = FALSE)
}

get_condition_coef <- function(dds) {
  rn <- resultsNames(dds)
  idx <- grep("^condition_", rn)
  if (length(idx) == 0) stop("Could not identify condition coefficient in resultsNames(dds).")
  rn[idx[1]]
}

extract_legend <- function(plot_obj) {
  gt <- ggplotGrob(plot_obj + theme(legend.position = "bottom"))
  guide_index <- which(vapply(gt$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(guide_index) == 0) return(NULL)
  gt$grobs[[guide_index[1]]]
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

pretty_dataset_type <- function(dataset_key) {
  switch(
    dataset_key,
    raw_dataset = "Raw",
    leading_edge_dataset = "Lead",
    remainder_dataset = "Rem",
    dataset_key
  )
}

pretty_dataset_label <- function(dataset_name) {
  parts <- strsplit(dataset_name, "_", fixed = TRUE)[[1]]
  if (length(parts) < 4) return(dataset_name)
  comparison_name <- paste(parts[1], parts[2], sep = "_")
  dataset_key <- paste(parts[3:length(parts)], collapse = "_")
  paste(comparison_name, pretty_dataset_type(dataset_key), sep = " ")
}

# =============================================================================
# DATA IMPORT
# =============================================================================

WTTS_Seq <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
WTTS_Seq <- as.data.frame(WTTS_Seq, stringsAsFactors = FALSE)

assert_required_columns(WTTS_Seq, c("OrigID", "Symbol"), "WTTS count file")
assert_required_columns(WTTS_Seq, meta_all$id, "WTTS count file sample columns")

WTTS_Seq$OrigID <- as.character(WTTS_Seq$OrigID)
WTTS_Seq$Symbol <- as.character(WTTS_Seq$Symbol)

WTTS_Seq <- WTTS_Seq[
  !is.na(WTTS_Seq$OrigID) & !is.na(WTTS_Seq$Symbol),
  ,
  drop = FALSE
]
WTTS_Seq <- WTTS_Seq[rowSums(is.na(WTTS_Seq[, meta_all$id, drop = FALSE])) == 0, , drop = FALSE]
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

# =============================================================================
# EVS PREPARATION
# =============================================================================

prepare_comparison_data <- function(comparison_name, group1_prefix, group2_prefix, WTTS_Seq, meta_all) {
  keep_ids <- grepl(paste0("^", group1_prefix, "_"), meta_all$id) |
              grepl(paste0("^", group2_prefix, "_"), meta_all$id)
  meta_sub <- meta_all[keep_ids, , drop = FALSE]
  coldata  <- meta_sub[, c("condition"), drop = FALSE]
  sample_ids <- rownames(meta_sub)

  missing_samples <- setdiff(sample_ids, colnames(WTTS_Seq))
  if (length(missing_samples) > 0) {
    stop(sprintf("Missing samples for %s: %s", comparison_name, paste(missing_samples, collapse = ", ")))
  }

  count_sub <- WTTS_Seq[, sample_ids, drop = FALSE]
  stopifnot(all(colnames(count_sub) == rownames(coldata)))

  list(
    comparison_name = comparison_name,
    count_matrix = as.matrix(count_sub),
    coldata = coldata
  )
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
    cutoff_quantile = cutoff_quantile
  )
}

compute_pc1_loading_table <- function(value_df, sample_names, top_n = top_n_target, preprocessing_label = "Norm") {
  x <- as.matrix(value_df[, sample_names, drop = FALSE])
  pca_fit <- prcomp(t(x), scale. = FALSE, rank. = 2)
  loading_abs <- abs(pca_fit$rotation[, 1])

  loading_tbl <- data.frame(
    feature_id = names(loading_abs),
    pc1_loading_abs = unname(loading_abs),
    stringsAsFactors = FALSE
  )
  loading_tbl <- loading_tbl[order(loading_tbl$pc1_loading_abs, decreasing = TRUE), , drop = FALSE]
  loading_tbl$rank <- seq_len(nrow(loading_tbl))

  cutoff_info <- resolve_top_n_cutoff(loading_tbl$pc1_loading_abs, top_n = top_n)
  cutoff <- cutoff_info$cutoff_value
  loading_tbl$split_class <- ifelse(loading_tbl$pc1_loading_abs >= cutoff, "high_loading", "background_loading")

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
  dds_init <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)
  dds_init <- dds_init[rowSums(counts(dds_init)) > 0, ]
  dds_init <- estimateSizeFactors(dds_init)

  norm_counts_init <- as.data.frame(counts(dds_init, normalized = TRUE))
  raw_counts_init  <- as.data.frame(count_matrix)

  sample_ids <- colnames(count_matrix)
  trt_ids   <- sample_ids[coldata$condition == "trt"]
  untrt_ids <- sample_ids[coldata$condition == "untrt"]

  fit_trt     <- compute_pc1_loading_table(norm_counts_init, trt_ids, preprocessing_label = "Norm")
  fit_untrt   <- compute_pc1_loading_table(norm_counts_init, untrt_ids, preprocessing_label = "Norm")
  fit_trt_raw   <- compute_pc1_loading_table(raw_counts_init, trt_ids, preprocessing_label = "Raw")
  fit_untrt_raw <- compute_pc1_loading_table(raw_counts_init, untrt_ids, preprocessing_label = "Raw")

  trt_high   <- as.character(subset(fit_trt$loading_table, split_class == "high_loading")$feature_id)
  untrt_high <- as.character(subset(fit_untrt$loading_table, split_class == "high_loading")$feature_id)

  leading_edge_ids <- union(trt_high, untrt_high)
  remainder_ids    <- setdiff(rownames(count_matrix), leading_edge_ids)

  if (length(leading_edge_ids) == 0) {
    stop("Leading-edge dataset is empty.")
  }

  list(
    fit_trt = fit_trt,
    fit_untrt = fit_untrt,
    fit_trt_raw = fit_trt_raw,
    fit_untrt_raw = fit_untrt_raw,
    leading_edge_ids = leading_edge_ids,
    remainder_ids = remainder_ids,
    raw_dataset = count_matrix,
    leading_edge_dataset = count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset = count_matrix[remainder_ids, , drop = FALSE]
  )
}

plot_pca_scatter <- function(pca_fit, dataset_label, group_label, preprocessing_label = "Norm") {
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
    geom_point(size = 3, colour = "white", stroke = 0.55) +
    geom_text_repel(size = 2.0, max.overlaps = 8, force = 1, box.padding = 0.22, point.padding = 0.10) +
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition") +
    scale_fill_manual(values = condition_fills, labels = condition_labels, name = "Condition") +
    labs(
      title = compact_title(paste(dataset_label, pretty_group_label(group_label), "PCA"), 40),
      subtitle = paste0(preprocessing_label, " preprocessing; PC1=", pca_var_per[1], "%, PC2=", pca_var_per[2], "%"),
      x = paste0("PC1 (", pca_var_per[1], "%)"),
      y = paste0("PC2 (", pca_var_per[2], "%)")
    ) +
    manuscript_theme()
}

plot_pc1_loading_rank <- function(loading_tbl, cutoff, dataset_label, group_label, top_n_used, cutoff_quantile, preprocessing_label = "Norm") {
  ggplot(loading_tbl, aes(rank, pc1_loading_abs)) +
    geom_line(linewidth = 0.4, color = "grey35") +
    geom_hline(yintercept = cutoff, color = plot_palette$threshold, linewidth = 0.8) +
    annotate(
      "text",
      x = max(loading_tbl$rank) * 0.72,
      y = cutoff,
      label = paste0("topN=", top_n_used, "  q=", signif(cutoff_quantile, 4)),
      color = plot_palette$threshold,
      vjust = -0.8,
      size = 3.1
    ) +
    labs(
      title = compact_title(paste(dataset_label, pretty_group_label(group_label), "PC1 rank"), 42),
      subtitle = paste0(preprocessing_label, " preprocessing"),
      x = "Ranked PAS",
      y = "|PC1 loading|"
    ) +
    manuscript_theme()
}

# =============================================================================
# CORE ANALYSIS
# =============================================================================

run_empirical_null_fdrtool <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
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
      fdrtool(
        stat_vec,
        statistic = "normal",
        plot = FALSE,
        verbose = FALSE,
        cutoff.method = "pct0",
        pct0 = 0.75
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
      error = function(e) NA_real_
    )
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

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
  dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, betaPrior = FALSE)

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

  res_df <- as.data.frame(res)
  res_df$feature_id <- rownames(res_df)

  valid_stat <- is.finite(res_df$stat) & !is.na(res_df$stat)
  stat_vec <- as.numeric(res_df$stat[valid_stat])
  message(sprintf("[%s] Mean Wald stat sent to fdrtool: %.4f", dataset_name, mean(stat_vec, na.rm = TRUE)))

  fdr_fit <- run_empirical_null_fdrtool(stat_vec, dataset_name)
  res_df$empirical_p <- NA_real_
  res_df$empirical_q <- NA_real_
  res_df$lfdr        <- NA_real_
  res_df$empirical_p[valid_stat] <- as.numeric(fdr_fit$pval)
  res_df$empirical_q[valid_stat] <- as.numeric(fdr_fit$qval)
  res_df$lfdr[valid_stat]        <- as.numeric(fdr_fit$lfdr)
  res_df$empirical_bh <- NA_real_

  ok_emp <- is.finite(res_df$empirical_p) & !is.na(res_df$empirical_p)
  if (any(ok_emp)) {
    res_df$empirical_bh[ok_emp] <- p.adjust(res_df$empirical_p[ok_emp], method = "BH")
  }

  coef_name <- get_condition_coef(dds)
  shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)
  shr_df <- as.data.frame(shr)
  shr_df$feature_id <- rownames(shr_df)

  res_df <- left_join(
    res_df,
    shr_df[, c("feature_id", "log2FoldChange")],
    by = "feature_id",
    suffix = c("", "_shr")
  )
  colnames(res_df)[colnames(res_df) == "log2FoldChange_shr"] <- "lfc_shrunk"

  res_strong_df <- as.data.frame(res_strong)
  res_strong_df$feature_id <- rownames(res_strong_df)

  res_weak_df <- as.data.frame(res_weak)
  res_weak_df$feature_id <- rownames(res_weak_df)

  res_df <- left_join(res_df, res_strong_df[, c("feature_id", "padj")], by = "feature_id", suffix = c("", "_strong"))
  res_df <- left_join(res_df, res_weak_df[, c("feature_id", "padj")], by = "feature_id", suffix = c("", "_weak"))
  colnames(res_df)[colnames(res_df) == "padj_strong"] <- "padj_strong_effect"
  colnames(res_df)[colnames(res_df) == "padj_weak"]   <- "padj_weak_effect"

  res_df$effect_class <- classify_effect_strength(
    res_df$padj_strong_effect,
    res_df$padj_weak_effect
  )
  res_df$resGA_padj <- res_df$padj_strong_effect
  res_df$resLA_padj <- res_df$padj_weak_effect

  hc_p_threshold_dataset <- safe_hc_thresh(res_df$empirical_p, dataset_name)
  res_df$gene_empirical_pvalue <- res_df$empirical_p
  res_df$HBFSS <- abs(res_df$lfc_shrunk * log10(pmax(res_df$gene_empirical_pvalue, 1e-300)))

  if (is.finite(hc_p_threshold_dataset) && !is.na(hc_p_threshold_dataset) &&
      hc_p_threshold_dataset > 0 && hc_p_threshold_dataset < 1) {
    hbfss_threshold_dataset <- abs(log10(hc_p_threshold_dataset)) * lfc_boundary
  } else {
    hbfss_threshold_dataset <- NA_real_
  }

  res_df$raw_lfc_pass <- !is.na(res_df$log2FoldChange) & abs(res_df$log2FoldChange) >= lfc_boundary
  res_df$shrunk_lfc_pass <- !is.na(res_df$lfc_shrunk) & abs(res_df$lfc_shrunk) >= lfc_boundary

  res_df$standard_significant <- !is.na(res_df$padj) &
    res_df$padj < alpha_level &
    res_df$raw_lfc_pass &
    res_df$shrunk_lfc_pass

  res_df$HBFSS_significant <- FALSE
  if (!is.na(hbfss_threshold_dataset)) {
    hbfss_native_call <-
      (res_df$HBFSS >= hbfss_threshold_dataset & (is.na(res_df$resLA_padj) | !(res_df$resLA_padj < 0.2))) |
      (!is.na(res_df$resGA_padj) & (res_df$resGA_padj < alpha_level))
    res_df$HBFSS_significant <- hbfss_native_call & res_df$raw_lfc_pass
  }

  res_df$regulation_direction <- ifelse(
    is.na(res_df$lfc_shrunk),
    NA_character_,
    ifelse(res_df$lfc_shrunk > 0, "upregulated",
           ifelse(res_df$lfc_shrunk < 0, "downregulated", "no_change"))
  )

  base_mean_vec <- res_df$baseMean[!is.na(res_df$baseMean)]

  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts$feature_id <- rownames(norm_counts)

  mm <- as.data.frame(mcols(dds))
  mm$feature_id <- rownames(mm)
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
    mutate(gene_symbol = if_else(is.na(gene_symbol), "", trimws(gene_symbol))) %>%
    arrange(feature_id, desc(gene_symbol != ""), gene_symbol) %>%
    distinct(feature_id, .keep_all = TRUE) %>%
    mutate(gene_symbol = na_if(gene_symbol, ""))

  final_df <- res_df %>%
    left_join(annot_df, by = "feature_id") %>%
    left_join(norm_counts[!duplicated(norm_counts$feature_id), , drop = FALSE], by = "feature_id") %>%
    left_join(disp_df[!duplicated(disp_df$feature_id), , drop = FALSE], by = "feature_id")

  final_df$neglog10_padj        <- safe_neglog10(final_df$padj)
  final_df$neglog10_empirical_p <- safe_neglog10(final_df$empirical_p)
  final_df$neglog10_padjc       <- safe_neglog10(final_df$empirical_bh)

  final_df$dataset_name            <- dataset_name
  final_df$hc_p_threshold_dataset  <- hc_p_threshold_dataset
  final_df$hbfss_threshold_dataset <- hbfss_threshold_dataset
  final_df$Apeglm_L2FC             <- final_df$lfc_shrunk
  final_df$pi_valueE               <- final_df$HBFSS
  final_df$pval                    <- final_df$empirical_p
  final_df$padjc                   <- final_df$empirical_bh
  final_df$qval                    <- final_df$empirical_q

  list(
    dds = dds,
    results = final_df,
    base_mean_vec = base_mean_vec,
    hc_p_threshold = hc_p_threshold_dataset,
    hbfss_threshold = hbfss_threshold_dataset
  )
}

# =============================================================================
# VOLCANO CLASSING / PLOTTING
# =============================================================================

volcano_legend_levels <- c("Weak CNH", "Strong CNH", "Standard", "HBFSS")
volcano_legend_colors <- c(
  "Weak CNH" = plot_palette$weak_cnh,
  "Strong CNH" = plot_palette$strong_cnh,
  "Standard" = plot_palette$standard,
  "HBFSS" = plot_palette$hbfss
)
volcano_legend_shapes <- c(
  "Weak CNH" = 16,
  "Strong CNH" = 17,
  "Standard" = 15,
  "HBFSS" = 18
)

build_volcano_classes <- function(df) {
  df <- as.data.frame(df)
  df$plot_class <- "Background"

  weak_only <- !is.na(df$resLA_padj) & df$resLA_padj < alpha_level &
    !(!is.na(df$resGA_padj) & df$resGA_padj < alpha_level) &
    !df$standard_significant &
    !df$HBFSS_significant

  strong_only <- !is.na(df$resGA_padj) & df$resGA_padj < alpha_level &
    !df$standard_significant &
    !df$HBFSS_significant

  std_only <- df$standard_significant & !df$HBFSS_significant
  hbfss_only <- df$HBFSS_significant & !df$standard_significant

  df$plot_class[weak_only]   <- "Weak CNH"
  df$plot_class[strong_only] <- "Strong CNH"
  df$plot_class[std_only]    <- "Standard"
  df$plot_class[hbfss_only]  <- "HBFSS"
  df$plot_class <- factor(df$plot_class, levels = c("Background", volcano_legend_levels))

  df$gene_symbol_plot <- ifelse(
    !is.na(df$gene_symbol) & grepl("[A-Za-z0-9]", trimws(df$gene_symbol)),
    trimws(df$gene_symbol),
    NA_character_
  )

  df
}

select_volcano_labels <- function(df, n_labels = n_top_labels) {
  df <- df[!is.na(df$gene_symbol_plot), , drop = FALSE]
  if (!nrow(df)) return(df[0, , drop = FALSE])

  priority <- match(as.character(df$plot_class), volcano_legend_levels)
  priority[is.na(priority)] <- 99

  ord <- order(priority, -df$neglog10_empirical_p, -abs(df$lfc_shrunk), na.last = TRUE)
  df <- df[ord, , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol_plot), , drop = FALSE]
  df[seq_len(min(n_labels, nrow(df))), , drop = FALSE]
}

build_hbfss_curve <- function(xvals, hc_p_threshold, lfc_boundary) {
  if (!is.finite(hc_p_threshold) || is.na(hc_p_threshold) || hc_p_threshold <= 0 || hc_p_threshold >= 1) {
    return(NULL)
  }
  hbfss_thr <- abs(log10(hc_p_threshold)) * lfc_boundary
  x_abs <- pmax(abs(xvals), 1e-6)
  y <- hbfss_thr / x_abs
  data.frame(x = xvals, y = y)
}

panel_overlap_text <- function(df) {
  std_n <- sum(df$standard_significant, na.rm = TRUE)
  hbf_n <- sum(df$HBFSS_significant, na.rm = TRUE)
  ov_n  <- sum(df$standard_significant & df$HBFSS_significant, na.rm = TRUE)
  paste0("Std=", std_n, "  HBFSS=", hbf_n, "  Overlap=", ov_n)
}

plot_integrated_volcano <- function(df, dataset_name) {
  df <- build_volcano_classes(df)
  lab_df <- select_volcano_labels(df)

  hc_p <- unique(df$hc_p_threshold_dataset)[1]
  hbfss_thr <- unique(df$hbfss_threshold_dataset)[1]
  hc_y <- if (is.finite(hc_p) && !is.na(hc_p) && hc_p > 0 && hc_p < 1) safe_neglog10(hc_p) else NA_real_

  xlim_max <- max(abs(df$lfc_shrunk), na.rm = TRUE)
  if (!is.finite(xlim_max) || is.na(xlim_max)) xlim_max <- 2
  xlim_max <- max(1.5, min(8, xlim_max * 1.08))
  x_grid <- seq(-xlim_max, xlim_max, length.out = 600)
  curve_df <- build_hbfss_curve(x_grid, hc_p, lfc_boundary)

  y_max <- max(df$neglog10_empirical_p, na.rm = TRUE)
  if (!is.finite(y_max) || is.na(y_max)) y_max <- 6
  y_max <- max(y_max, if (is.finite(hc_y)) hc_y + 0.5 else 0) + 0.8

  p <- ggplot(df, aes(lfc_shrunk, neglog10_empirical_p)) +
    geom_point(color = plot_palette$background, alpha = 0.70, size = 0.95) +
    geom_point(
      data = subset(df, plot_class != "Background"),
      aes(color = plot_class, shape = plot_class),
      alpha = 0.92,
      size = 1.35,
      stroke = 0.22
    ) +
    geom_vline(xintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", linewidth = 0.40, colour = plot_palette$threshold) +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.25, colour = "grey55") +
    scale_color_manual(values = volcano_legend_colors, breaks = volcano_legend_levels, drop = FALSE, name = "Class") +
    scale_shape_manual(values = volcano_legend_shapes, breaks = volcano_legend_levels, drop = FALSE, name = "Class") +
    labs(
      title = compact_title(gsub(" \\| ", " ", pretty_dataset_label(dataset_name)), 22),
      subtitle = "x = shrunken log2FC  ·  y = -log10(empirical p)",
      x = "Shrunken log2FC",
      y = expression(-log[10]("Empirical p"))
    ) +
    coord_cartesian(xlim = c(-xlim_max, xlim_max), ylim = c(0, y_max), clip = "off") +
    manuscript_theme() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = base_theme_size - 1),
      legend.text = element_text(size = base_theme_size - 2),
      plot.title = element_text(size = base_theme_size),
      plot.subtitle = element_text(size = base_theme_size - 2),
      axis.title = element_text(size = base_theme_size - 1),
      axis.text = element_text(size = base_theme_size - 2)
    )

  if (!is.na(hc_y) && is.finite(hc_y)) {
    p <- p +
      geom_hline(yintercept = hc_y, linetype = "dotted", linewidth = 0.40, colour = plot_palette$threshold) +
      annotate(
        "text",
        x = -xlim_max * 0.95,
        y = hc_y,
        label = paste0("HC=", signif(hc_p, 3)),
        hjust = 0,
        vjust = -0.35,
        size = 2.25,
        colour = plot_palette$threshold
      )
  }

  if (!is.null(curve_df) && is.finite(hbfss_thr) && !is.na(hbfss_thr)) {
    curve_df <- subset(curve_df, is.finite(y) & y <= y_max)
    p <- p +
      geom_line(data = curve_df, aes(x, y), inherit.aes = FALSE, linewidth = 0.45, colour = plot_palette$hbfss) +
      annotate(
        "text",
        x = xlim_max * 0.95,
        y = min(y_max - 0.25, max(0.4, hbfss_thr / max(lfc_boundary, 1e-6))),
        label = paste0("HBFSS=", signif(hbfss_thr, 3)),
        hjust = 1,
        vjust = -0.25,
        size = 2.25,
        colour = plot_palette$hbfss
      )
  }

  p <- p + annotate("text", x = 0, y = y_max * 0.96, label = panel_overlap_text(df), size = 2.25)

  if (nrow(lab_df) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = lab_df,
        aes(label = gene_symbol_plot, color = plot_class),
        size = 1.55,
        max.overlaps = 18,
        seed = 1,
        box.padding = 0.22,
        point.padding = 0.10,
        min.segment.length = 0,
        segment.alpha = 0.45,
        show.legend = FALSE
      )
  }

  p
}

plot_integrated_volcano_compare <- function(analysis_results, comparison_name) {
  key_order <- c("raw_dataset", "leading_edge_dataset", "remainder_dataset")
  key_order <- key_order[key_order %in% names(analysis_results)]

  grobs <- lapply(key_order, function(k) {
    plot_integrated_volcano(
      analysis_results[[k]]$results,
      analysis_results[[k]]$summary$dataset_name[1]
    ) + theme(legend.position = "none")
  })

  legend_grob <- extract_legend(
    plot_integrated_volcano(
      analysis_results[[key_order[1]]]$results,
      analysis_results[[key_order[1]]]$summary$dataset_name[1]
    )
  )

  title_grob <- textGrob(
    paste0(comparison_name, " integrated volcano"),
    gp = gpar(fontface = "bold", cex = 0.95)
  )

  arrangeGrob(
    grobs = c(grobs, list(legend_grob)),
    layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
    top = title_grob,
    heights = unit.c(unit(1, "null"), unit(0.14, "null"))
  )
}

# =============================================================================
# OTHER PLOTS
# =============================================================================

plot_empirical_p_histogram <- function(df, dataset_name, hc_p_threshold) {
  subtitle_text <- if (is.na(hc_p_threshold) || !is.finite(hc_p_threshold) || hc_p_threshold <= 0 || hc_p_threshold >= 1) {
    "Empirical-null p-values"
  } else {
    paste0("Empirical-null p-values; HC=", signif(hc_p_threshold, 4))
  }

  p <- ggplot(df, aes(empirical_p)) +
    geom_histogram(bins = 60, fill = plot_palette$hist, color = "white") +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "Empirical p"), 38),
      subtitle = subtitle_text,
      x = "Empirical p",
      y = "Count"
    ) +
    manuscript_theme()

  if (!is.na(hc_p_threshold) && is.finite(hc_p_threshold) &&
      hc_p_threshold > 0 && hc_p_threshold < 1) {
    p <- p + geom_vline(xintercept = hc_p_threshold, color = plot_palette$hbfss, linewidth = 0.8)
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
    geom_line(linewidth = 0.45, colour = "grey30") +
    geom_point(data = peak_df, aes(empirical_p, hc_score), inherit.aes = FALSE, size = 2.0, colour = plot_palette$hbfss) +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "HC profile"), 38),
      subtitle = if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold)) paste0("HC threshold=", signif(hc_p_threshold, 4)) else "HC threshold unavailable",
      x = "Sorted empirical p",
      y = "HC score"
    ) +
    manuscript_theme()

  if (!is.na(hc_p_threshold) && is.finite(hc_p_threshold) &&
      hc_p_threshold > 0 && hc_p_threshold < 1) {
    p <- p + geom_vline(xintercept = hc_p_threshold, color = plot_palette$hbfss, linewidth = 0.8)
  }

  p
}

plot_shrunken_ma <- function(df, dataset_name) {
  ggplot(df, aes(safe_log10(baseMean + 1), lfc_shrunk)) +
    geom_point(alpha = 0.55, size = 1.0, color = "grey40") +
    geom_hline(yintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", color = plot_palette$threshold, linewidth = 0.6) +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "MA"), 38),
      subtitle = "Shrunken effect vs baseMean",
      x = expression(log[10]("baseMean + 1")),
      y = "Shrunken log2FC"
    ) +
    manuscript_theme()
}

plot_dispersion_panel_for_dataset <- function(df, dataset_name) {
  df <- build_volcano_classes(df)

  ggplot(df, aes(baseMean, dispersion, color = plot_class, shape = plot_class)) +
    geom_point(alpha = 0.55, size = 1.15, stroke = 0.2) +
    scale_x_log10(labels = label_number(accuracy = 0.1)) +
    scale_y_log10(labels = label_number(accuracy = 0.1)) +
    scale_color_manual(values = c("Background" = plot_palette$background, volcano_legend_colors), drop = FALSE, name = "Class") +
    scale_shape_manual(values = c("Background" = 16, volcano_legend_shapes), drop = FALSE, name = "Class") +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "Dispersion"), 38),
      subtitle = "Final dispersion vs mean",
      x = "baseMean",
      y = "Dispersion"
    ) +
    manuscript_theme() +
    theme(legend.position = "bottom")
}

plot_mean_histogram_panel <- function(df_means, base_mean_vec, dataset_name) {
  p1 <- ggplot(data.frame(x = safe_log10(df_means$mean_trt + 1)), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$hist, color = "white") +
    labs(title = "Trt mean", x = "log10(mean+1)", y = "Count") +
    manuscript_theme()

  p2 <- ggplot(data.frame(x = safe_log10(df_means$mean_untrt + 1)), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$hist, color = "white") +
    labs(title = "Ctrl mean", x = "log10(mean+1)", y = "Count") +
    manuscript_theme()

  p3 <- ggplot(data.frame(x = safe_log10(df_means$mean_all + 1)), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$hist, color = "white") +
    labs(title = "All mean", x = "log10(mean+1)", y = "Count") +
    manuscript_theme()

  p4 <- ggplot(data.frame(x = safe_log10(base_mean_vec + 1)), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$hist, color = "white") +
    labs(title = "baseMean", x = "log10(baseMean+1)", y = "Count") +
    manuscript_theme()

  arrangeGrob(
    p1, p2, p3, p4, ncol = 2,
    top = textGrob(
      paste(pretty_dataset_label(dataset_name), "Mean histograms"),
      gp = gpar(fontface = "bold", cex = 1.0)
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

# =============================================================================
# COMPARISON PANELS
# =============================================================================

save_cross_dataset_comparison_panels <- function(comparison_name, analysis_results, cmp_dir, dataset_list, coldata) {
  keys_present <- intersect(c("raw_dataset", "leading_edge_dataset", "remainder_dataset"), names(analysis_results))
  if (!length(keys_present)) return(invisible(NULL))

  volcano_panel <- plot_integrated_volcano_compare(analysis_results, comparison_name)
  save_grob(
    volcano_panel,
    file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_Volcano.png")),
    width = 13.5,
    height = 6.4
  )

  disp_grobs <- lapply(keys_present, function(k) {
    plot_dispersion_panel_for_dataset(
      analysis_results[[k]]$results,
      analysis_results[[k]]$summary$dataset_name[1]
    ) + theme(legend.position = "none")
  })
  disp_legend <- extract_legend(
    plot_dispersion_panel_for_dataset(
      analysis_results[[keys_present[1]]]$results,
      analysis_results[[keys_present[1]]]$summary$dataset_name[1]
    )
  )
  disp_panel <- arrangeGrob(
    grobs = c(disp_grobs, list(disp_legend)),
    layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
    top = textGrob(paste0(comparison_name, " dispersion"), gp = gpar(fontface = "bold", cex = 0.95)),
    heights = unit.c(unit(1, "null"), unit(0.14, "null"))
  )
  save_grob(
    disp_panel,
    file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_Disp.png")),
    width = 13.5,
    height = 6.4
  )

  hist_grobs <- lapply(keys_present, function(k) {
    plot_empirical_p_histogram(
      analysis_results[[k]]$results,
      analysis_results[[k]]$summary$dataset_name[1],
      analysis_results[[k]]$summary$hc_p_threshold[1]
    )
  })
  hist_panel <- do.call(
    arrangeGrob,
    c(hist_grobs, list(
      ncol = length(hist_grobs),
      top = textGrob(paste0(comparison_name, " empirical p"), gp = gpar(fontface = "bold", cex = 0.95))
    ))
  )
  save_grob(
    hist_panel,
    file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_EmpP.png")),
    width = 12.5,
    height = 4.8
  )

  summary_table <- bind_rows(lapply(keys_present, function(k) {
    sm <- analysis_results[[k]]$summary
    data.frame(
      comparison_name = comparison_name,
      dataset_key = k,
      dataset_label = unname(dataset_short[k]),
      dataset_name = sm$dataset_name[1],
      n_features = sm$n_features[1],
      n_standard = sm$n_standard_significant[1],
      n_hbfss = sm$n_HBFSS_significant[1],
      n_overlap = sm$n_overlap_significant[1],
      n_weak_cnh = sm$n_weak_cnh[1],
      n_strong_cnh = sm$n_strong_cnh[1],
      hc_p_threshold = sm$hc_p_threshold[1],
      hbfss_threshold = sm$hbfss_threshold[1],
      stringsAsFactors = FALSE
    )
  }))
  save_csv(summary_table, file.path(cmp_dir, paste0("Table_", comparison_name, "_Compare_Summary.csv")))
}

# =============================================================================
# FULL COMPARISON PIPELINE
# =============================================================================

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  cmp_dir     <- file.path(output_dir, comparison_name)
  tab_dir     <- file.path(cmp_dir, "tables")
  fig_raw_dir <- file.path(cmp_dir, "fig_raw")
  fig_le_dir  <- file.path(cmp_dir, "fig_lead")
  fig_rem_dir <- file.path(cmp_dir, "fig_rem")

  for (d in c(cmp_dir, tab_dir, fig_raw_dir, fig_le_dir, fig_rem_dir)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  evs <- build_eigenvector_split(count_matrix, coldata)

  save_plot(
    plot_pca_scatter(evs$fit_trt$pca_fit, comparison_name, "treatment", "Norm"),
    file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_PCA.png")),
    width = 8.6, height = 6.2
  )
  save_plot(
    plot_pca_scatter(evs$fit_untrt$pca_fit, comparison_name, "control", "Norm"),
    file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_PCA.png")),
    width = 8.6, height = 6.2
  )
  save_plot(
    plot_pc1_loading_rank(
      evs$fit_trt$loading_table,
      evs$fit_trt$cutoff,
      comparison_name,
      "treatment",
      evs$fit_trt$top_n_used,
      evs$fit_trt$cutoff_quantile,
      "Norm"
    ),
    file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_Rank.png")),
    width = 9.4, height = 6.3
  )
  save_plot(
    plot_pc1_loading_rank(
      evs$fit_untrt$loading_table,
      evs$fit_untrt$cutoff,
      comparison_name,
      "control",
      evs$fit_untrt$top_n_used,
      evs$fit_untrt$cutoff_quantile,
      "Norm"
    ),
    file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_Rank.png")),
    width = 9.4, height = 6.3
  )

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

    df <- fit$results

    summary_row <- data.frame(
      comparison_name = comparison_name,
      dataset_name = full_dataset_name,
      n_features = nrow(df),
      hc_p_threshold = fit$hc_p_threshold,
      hbfss_threshold = fit$hbfss_threshold,
      n_standard_significant = sum(df$standard_significant, na.rm = TRUE),
      n_HBFSS_significant = sum(df$HBFSS_significant, na.rm = TRUE),
      n_overlap_significant = sum(df$standard_significant & df$HBFSS_significant, na.rm = TRUE),
      n_weak_cnh = sum(!is.na(df$resLA_padj) & df$resLA_padj < alpha_level, na.rm = TRUE),
      n_strong_cnh = sum(!is.na(df$resGA_padj) & df$resGA_padj < alpha_level, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    save_csv(df, tab_file(tab_dir, comparison_name, "Results", nm))
    save_csv(summary_row, tab_file(tab_dir, comparison_name, "Summary", nm))
    save_csv(subset(df, standard_significant), tab_file(tab_dir, comparison_name, "Std", nm))
    save_csv(subset(df, HBFSS_significant), tab_file(tab_dir, comparison_name, "HBFSS", nm))

    mean_df <- compute_mean_expression_table(dataset_list[[nm]], coldata)
    mean_panel <- plot_mean_histogram_panel(mean_df, fit$base_mean_vec, full_dataset_name)
    save_grob(mean_panel, fig_file(fig_subdir, comparison_name, "MeanHist", nm), width = 12, height = 8)

    p_hist <- plot_empirical_p_histogram(df, full_dataset_name, fit$hc_p_threshold)
    p_hc   <- plot_hc_profile(df, full_dataset_name, fit$hc_p_threshold)
    p_ma   <- plot_shrunken_ma(df, full_dataset_name)
    p_vol  <- plot_integrated_volcano(df, full_dataset_name)

    save_plot(p_hist, fig_file(fig_subdir, comparison_name, "EmpP", nm), width = 8.6, height = 6.0)
    if (!is.null(p_hc)) {
      save_plot(p_hc, fig_file(fig_subdir, comparison_name, "HC", nm), width = 8.6, height = 6.0)
    }
    save_plot(p_ma, fig_file(fig_subdir, comparison_name, "MA", nm), width = 8.6, height = 6.0)
    save_plot(p_vol, fig_file(fig_subdir, comparison_name, "Volcano", nm), width = 8.6, height = 6.4)

    analysis_results[[nm]] <- list(
      dds = fit$dds,
      results = df,
      summary = summary_row,
      dataset_mat = dataset_list[[nm]],
      fig_subdir = fig_subdir
    )
  }

  save_cross_dataset_comparison_panels(comparison_name, analysis_results, cmp_dir, dataset_list, coldata)
  bind_rows(lapply(analysis_results, `[[`, "summary"))
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
  bind_rows(all_summaries_list)
} else {
  data.frame()
}

if (nrow(all_summaries) > 0) {
  save_csv(all_summaries, file.path(output_dir, "Table_Overall_Summary.csv"))
}

if (length(failed_comparisons) > 0) {
  save_csv(bind_rows(failed_comparisons), file.path(output_dir, "Table_Failed.csv"))
}

unlink("Rplots.pdf", force = TRUE)

cat("\n=====================================================\n")
cat("Pipeline complete.\n")
cat("Using count file:\n")
cat(normalizePath(count_file), "\n")
cat("Output directory:\n")
cat(normalizePath(output_dir), "\n")
cat("=====================================================\n\n")

if (nrow(all_summaries) > 0) {
  print(all_summaries)
} else {
  message("No comparison summaries were written. Check Table_Failed.csv for the exact error message.")
}
