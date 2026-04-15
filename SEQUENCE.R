#!/usr/bin/env Rscript

# =============================================================================
# SEQUENCE.R
# FINAL MANUSCRIPT PIPELINE
#
# PURPOSE
# This script builds the final manuscript-ready EVS + DESeq2 + empirical-null +
# HBFSS analysis pipeline for the WTTS PAS dataset and exports only the figure
# and table sets that are useful for the manuscript.
#
# CORE MANUSCRIPT CLAIM
# 1. Eigenvector splitting (EVS) separates the ranked PAS universe into an
#    original dataset, a leading-edge dataset, and a remainder dataset.
# 2. The leading-edge dataset enriches the structured signal-bearing PAS set.
# 3. Standard DESeq2, composite-null weak/strong effect testing, and HBFSS are
#    shown together on a single integrated volcano panel so that the manuscript
#    can compare classical and hybrid signal calls directly.
#
# DESIGN PRINCIPLES FOR THIS FINAL VERSION
# - Short export names so GitHub folders are easy to scan.
# - One shared legend per multi-panel figure.
# - Short panel titles and short subtitles.
# - No repeated legends inside subplot grids.
# - Only manuscript-useful figures are exported by default.
# - Methods are fully explained in comments, not crowded into plotting space.
#
# DIFFERENTIAL EXPRESSION DESIGN
# All DESeq2 models use:
#   design = ~ condition
# with condition coded as:
#   untrt = control
#   trt   = treatment
#
# UNIT OF OBSERVATION
# The feature is a PAS / OrigID, not a gene. Gene symbols are annotation only.
#
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(fdrtool)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(grid)
  library(gtable)
  library(scales)
  library(grDevices)
})

options(stringsAsFactors = FALSE)
# Avoid accidental default-device files such as Rplots.pdf during batch export.
options(device = function(...) grDevices::png(filename = tempfile(fileext = ".png"), width = 960, height = 720, res = 120))

# =============================================================================
# SECTION 1
# SETTINGS, PATHS, METADATA, AND GLOBAL STYLING
# =============================================================================

# -----------------------------------------------------------------------------
# INPUT FILE DISCOVERY
#
# The script first tries the repository-standard data path, then the current
# working directory, then an absolute fallback path. This fixes the prior
# "file not found" failure when the script was run from the repository root.
# -----------------------------------------------------------------------------
count_file_candidates <- c(
  "data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  "WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
)

count_file <- count_file_candidates[file.exists(count_file_candidates)][1]
if (is.na(count_file) || !nzchar(count_file)) {
  stop(
    "Could not find WTTS count file. Tried:\n",
    paste(" -", count_file_candidates, collapse = "\n")
  )
}

# Optional TWAS overlap
run_twas_overlap <- FALSE
twas_file <- "3aTWAS_genes_of_11_brain_disorders.csv"

# Global significance settings
alpha_level  <- 0.10
lfc_boundary <- 1.0
top_n_target <- 5000L

# Plot settings
figure_dpi            <- 320
base_theme_size       <- 10
n_top_labels_integrated <- 18

# If TRUE, additional diagnostic figures are exported beyond the core manuscript
# figure set. Default FALSE keeps the output manuscript-focused.
export_optional_diagnostics <- FALSE

# Output root
output_dir <- "exports/manuscript_final_clean"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# SHORT EXPORT NAMING
#
# Every file begins with Figure_ or Table_. Dataset abbreviations are short.
# -----------------------------------------------------------------------------
dataset_short <- c(
  raw_dataset          = "Raw",
  leading_edge_dataset = "Lead",
  remainder_dataset    = "Rem"
)

fig_file <- function(dir, cmp, ds_key, tag) {
  file.path(dir, paste0("Figure_", cmp, "_", unname(dataset_short[ds_key]), "_", tag, ".png"))
}

tab_file <- function(dir, cmp, ds_key, tag) {
  file.path(dir, paste0("Table_", cmp, "_", unname(dataset_short[ds_key]), "_", tag, ".csv"))
}

# -----------------------------------------------------------------------------
# FIXED PALETTE
#
# Weak effect is blue as requested.
# Strong composite-null effect is red.
# Standard DESeq2 is dark green.
# HBFSS-only is orange.
# Overlap classes get distinct colors and shapes so they do not disappear.
# -----------------------------------------------------------------------------
plot_palette <- list(
  background = "#BDBDBD",
  weak       = "#1F78B4",
  strong     = "#B2182B",
  standard   = "#1B9E77",
  hbfss      = "#E66101",
  std_hbfss  = "#4D9221",
  weak_hbfss = "#6A3D9A",
  strong_hbfss = "#D95F02",
  threshold  = "#8C2D04",
  histogram  = "#969696",
  control    = "#4D4D4D",
  treatment  = "#1F78B4"
)

condition_shapes <- c("untrt" = 21, "trt" = 24)
condition_fills  <- c("untrt" = plot_palette$control, "trt" = plot_palette$treatment)
condition_labels <- c("untrt" = "Control", "trt" = "Treatment")

integrated_class_levels <- c(
  "Weak CNH",
  "Strong CNH",
  "Standard",
  "HBFSS"
)

integrated_class_colors <- c(
  "Weak CNH"   = plot_palette$weak,
  "Strong CNH" = plot_palette$strong,
  "Standard"   = plot_palette$standard,
  "HBFSS"      = plot_palette$hbfss
)

integrated_class_shapes <- c(
  "Weak CNH"   = 16,
  "Strong CNH" = 17,
  "Standard"   = 15,
  "HBFSS"      = 18
)

# -----------------------------------------------------------------------------
# MINIMAL EMBEDDED METADATA
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

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix   = c("R0", "R2", "R4", "R8"),
  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

# =============================================================================
# SECTION 2
# HELPERS AND GLOBAL UTILITY FUNCTIONS
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

compact_title <- function(x, width = 44) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

compact_caption <- function(x, width = 110) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

save_csv <- function(df, path) {
  write.csv(df, file = path, row.names = FALSE)
}

save_plot <- function(p, path, width = 11.5, height = 7.8, dpi = figure_dpi, bg = "white") {
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

save_grob <- function(g, path, width = 14, height = 8.5, dpi = figure_dpi, bg = "white") {
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

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_theme_size + 1, hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_theme_size - 1, hjust = 0.5, margin = margin(b = 6)),
      plot.caption     = element_text(size = base_theme_size - 3, hjust = 0.5, colour = "grey30", margin = margin(t = 6)),
      axis.title       = element_text(face = "bold"),
      axis.text        = element_text(colour = "black"),
      legend.title     = element_text(face = "bold"),
      legend.position  = "bottom",
      legend.box       = "vertical",
      legend.margin    = margin(t = 2, r = 2, b = 2, l = 2),
      legend.spacing.x = unit(4, "pt"),
      legend.spacing.y = unit(2, "pt"),
      legend.text      = element_text(size = base_theme_size - 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      plot.margin      = margin(t = 10, r = 14, b = 12, l = 12)
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
  cmp <- paste(parts[1], parts[2], sep = "_")
  ds_key <- paste(parts[3:length(parts)], collapse = "_")
  paste(cmp, pretty_dataset_type(ds_key))
}

pretty_group_label <- function(group_label) {
  switch(
    group_label,
    trt = "Trt",
    untrt = "Ctrl",
    treatment = "Trt",
    control = "Ctrl",
    group_label
  )
}

make_design_formula <- function(coldata) {
  ~ condition
}

get_condition_coef <- function(dds) {
  rn <- resultsNames(dds)
  idx <- grep("^condition_", rn)
  if (length(idx) == 0) stop("Could not identify condition coefficient in resultsNames(dds).")
  rn[idx[1]]
}

clean_gene_set <- function(x) {
  unique(tolower(trimws(x[!is.na(x) & x != ""])))
}

resolve_top_n_cutoff <- function(sorted_values_desc, top_n = top_n_target) {
  n_total <- length(sorted_values_desc)
  if (n_total == 0) stop("resolve_top_n_cutoff() received an empty vector.")
  top_n_actual <- min(max(1L, as.integer(top_n)), n_total)
  cutoff_value <- sorted_values_desc[top_n_actual]
  cutoff_quantile <- 1 - (top_n_actual / n_total)
  list(
    top_n_actual    = top_n_actual,
    cutoff_value    = cutoff_value,
    cutoff_quantile = cutoff_quantile,
    n_total         = n_total
  )
}

# -----------------------------------------------------------------------------
# EMPIRICAL-NULL FIT
#
# The DESeq2 Wald statistics are recalibrated using fdrtool under a normal
# empirical-null model. Only finite Wald statistics are passed to fdrtool.
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# LEGEND EXTRACTION
#
# Multi-panel figures use one shared legend only.
# -----------------------------------------------------------------------------
extract_shared_legend <- function(plot_obj) {
  gt <- ggplotGrob(plot_obj)
  guide_index <- which(vapply(gt$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(guide_index) == 0) return(NULL)
  gt$grobs[[guide_index[1]]]
}

remove_legend <- function(p) {
  p + theme(legend.position = "none")
}

assemble_panel_with_shared_legend <- function(plot_list, ncol, legend_grob, top = NULL) {
  plot_grob <- arrangeGrob(grobs = plot_list, ncol = ncol, top = top)
  if (is.null(legend_grob)) return(plot_grob)
  arrangeGrob(
    plot_grob,
    legend_grob,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - unit(1.2, "in"), unit(1.2, "in"))
  )
}

# =============================================================================
# SECTION 3
# DATA IMPORT AND ANNOTATION PREPARATION
# =============================================================================

message("Using count file: ", normalizePath(count_file))

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
# SECTION 4
# COMPARISON PREPARATION AND EVS
# =============================================================================

prepare_comparison_data <- function(comparison_name, group1_prefix, group2_prefix, WTTS_Seq, meta_all) {
  keep_ids <- grepl(paste0("^", group1_prefix, "_"), meta_all$id) |
    grepl(paste0("^", group2_prefix, "_"), meta_all$id)

  meta_sub   <- meta_all[keep_ids, , drop = FALSE]
  coldata    <- meta_sub[, c("condition"), drop = FALSE]
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

# -----------------------------------------------------------------------------
# EVS STEP
#
# PCA is run separately within each condition. Features are ranked by absolute
# PC1 loading. The leading-edge dataset is the union of the top-N features from
# treatment and control. The remainder dataset contains all other PAS features.
# -----------------------------------------------------------------------------
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

  dds_init <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = design_formula
  )
  dds_init <- dds_init[rowSums(counts(dds_init)) > 0, ]
  dds_init <- estimateSizeFactors(dds_init)

  norm_counts_init <- as.data.frame(counts(dds_init, normalized = TRUE))
  raw_counts_init  <- as.data.frame(count_matrix)

  sample_ids <- colnames(count_matrix)
  trt_ids    <- sample_ids[coldata$condition == "trt"]
  untrt_ids  <- sample_ids[coldata$condition == "untrt"]

  fit_trt <- compute_pc1_loading_table(
    norm_counts_init, trt_ids,
    top_n = top_n_target,
    preprocessing_label = "Normalized prior to eigenvector splitting"
  )
  fit_untrt <- compute_pc1_loading_table(
    norm_counts_init, untrt_ids,
    top_n = top_n_target,
    preprocessing_label = "Normalized prior to eigenvector splitting"
  )

  fit_trt_raw <- compute_pc1_loading_table(
    raw_counts_init, trt_ids,
    top_n = top_n_target,
    preprocessing_label = "Eigenvector splitting without prior normalization"
  )
  fit_untrt_raw <- compute_pc1_loading_table(
    raw_counts_init, untrt_ids,
    top_n = top_n_target,
    preprocessing_label = "Eigenvector splitting without prior normalization"
  )

  trt_high   <- as.character(subset(fit_trt$loading_table, split_class == "high_loading")$feature_id)
  untrt_high <- as.character(subset(fit_untrt$loading_table, split_class == "high_loading")$feature_id)

  leading_edge_ids <- union(trt_high, untrt_high)
  remainder_ids    <- setdiff(rownames(count_matrix), leading_edge_ids)

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
    raw_dataset = count_matrix,
    leading_edge_dataset = count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset = count_matrix[remainder_ids, , drop = FALSE]
  )
}

# =============================================================================
# SECTION 5
# CORE DESEQ2 + EMPIRICAL NULL + HBFSS ANALYSIS
# =============================================================================

# -----------------------------------------------------------------------------
# HBFSS / HYBRID LOGIC
#
# Definitions:
#   mu        = empirical mean across samples
#   variance  = empirical variance across samples
#   log_nb2   = log(1 + variance - mu)
#   nb_gap    = log(1 + variance - mu) - log(1 + mu)
#   alpha_mu  = log(1 + alpha*mu), where alpha = max((variance - mu)/mu^2, 0)
#
# Why these quantities are used:
# - variance - mu is the extra-Poisson variance term.
# - nb_gap compares extra-Poisson variance against mean-linked signal.
# - alpha*mu reflects the NB2 variance identity variance = mu + alpha*mu^2.
#
# These quantities are used descriptively and comparatively. They are not a
# formal likelihood-ratio test between NB1 and NB2.
# -----------------------------------------------------------------------------
classify_effect_strength <- function(res_strong_padj, res_weak_padj, alpha = alpha_level) {
  out <- rep("intermediate", length(res_strong_padj))
  out[!is.na(res_weak_padj) & res_weak_padj < alpha] <- "weak_effect"
  out[!is.na(res_strong_padj) & res_strong_padj < alpha] <- "strong_effect"
  out
}

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

  message(sprintf("[%s] Mean Wald stat sent to fdrtool: %.4f", dataset_name, mean(stat_vec, na.rm = TRUE)))

  fdr_fit <- run_empirical_null_fdrtool(stat_vec, dataset_name = dataset_name)

  res_df <- res_all_df
  n_valid <- sum(valid_stat)

  if (length(fdr_fit$pval) != n_valid ||
      length(fdr_fit$qval) != n_valid ||
      length(fdr_fit$lfdr) != n_valid) {
    stop(sprintf("[%s] fdrtool output length mismatch.", dataset_name))
  }

  res_df$empirical_p <- NA_real_
  res_df$empirical_q <- NA_real_
  res_df$lfdr        <- NA_real_

  res_df$empirical_p[valid_stat] <- as.numeric(fdr_fit$pval)
  res_df$empirical_q[valid_stat] <- as.numeric(fdr_fit$qval)
  res_df$lfdr[valid_stat]        <- as.numeric(fdr_fit$lfdr)

  res_df$empirical_bh <- NA_real_
  valid_empirical <- is.finite(res_df$empirical_p) & !is.na(res_df$empirical_p)
  if (any(valid_empirical)) {
    res_df$empirical_bh[valid_empirical] <- p.adjust(res_df$empirical_p[valid_empirical], method = "BH")
  }

  # Legacy aliases
  res_df$pval  <- res_df$empirical_p
  res_df$padjc <- res_df$empirical_bh
  res_df$qval  <- res_df$empirical_q

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

  hc_p_threshold_dataset <- safe_hc_thresh(res_df$empirical_p, dataset_name = dataset_name)
  res_df$gene_empirical_pvalue <- res_df$empirical_p
  res_df$HBFSS <- abs(res_df$lfc_shrunk * log10(pmax(res_df$gene_empirical_pvalue, 1e-300)))

  if (is.na(hc_p_threshold_dataset) || hc_p_threshold_dataset <= 0 || hc_p_threshold_dataset >= 1) {
    hbfss_threshold_dataset  <- NA_real_
    res_df$HBFSS_core_pass   <- FALSE
    res_df$HBFSS_significant <- FALSE
  } else {
    hbfss_threshold_dataset  <- abs(log10(hc_p_threshold_dataset)) * lfc_boundary
    res_df$HBFSS_core_pass   <- res_df$HBFSS >= hbfss_threshold_dataset
    res_df$HBFSS_significant <- res_df$HBFSS_core_pass
  }

  res_df$regulation_direction <- ifelse(
    is.na(res_df$lfc_shrunk), NA_character_,
    ifelse(res_df$lfc_shrunk > 0, "upregulated",
           ifelse(res_df$lfc_shrunk < 0, "downregulated", "no_change"))
  )

  res_df$raw_lfc_pass <- !is.na(res_df$log2FoldChange) &
    abs(res_df$log2FoldChange) >= lfc_boundary

  res_df$shrunk_lfc_pass <- !is.na(res_df$lfc_shrunk) &
    abs(res_df$lfc_shrunk) >= lfc_boundary

  res_df$standard_significant <- !is.na(res_df$padj) &
    (res_df$padj < alpha_level) &
    res_df$raw_lfc_pass &
    res_df$shrunk_lfc_pass

  res_strong_df <- as.data.frame(res_strong)
  res_strong_df$feature_id <- as.character(rownames(res_strong_df))
  res_weak_df <- as.data.frame(res_weak)
  res_weak_df$feature_id <- as.character(rownames(res_weak_df))

  res_df <- left_join(res_df, res_strong_df[, c("feature_id", "padj")], by = "feature_id", suffix = c("", "_strong"))
  res_df <- left_join(res_df, res_weak_df[, c("feature_id", "padj")], by = "feature_id", suffix = c("", "_weak"))

  colnames(res_df)[colnames(res_df) == "padj_strong"] <- "padj_strong_effect"
  colnames(res_df)[colnames(res_df) == "padj_weak"]   <- "padj_weak_effect"

  res_df$effect_class <- classify_effect_strength(
    res_df$padj_strong_effect,
    res_df$padj_weak_effect,
    alpha = alpha_level
  )

  res_df$resGA_padj <- res_df$padj_strong_effect
  res_df$resLA_padj <- res_df$padj_weak_effect

  hbfss_native_call <-
    (res_df$HBFSS_core_pass & (is.na(res_df$resLA_padj) | !(res_df$resLA_padj < 0.2))) |
    (!is.na(res_df$resGA_padj) & (res_df$resGA_padj < alpha_level))
  res_df$HBFSS_significant <- hbfss_native_call & res_df$raw_lfc_pass

  # Feature-level NB2-related diagnostics for summary tables
  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts$feature_id <- as.character(rownames(norm_counts))

  mu <- rowMeans(norm_counts[, colnames(count_mat), drop = FALSE], na.rm = TRUE)
  variance <- apply(norm_counts[, colnames(count_mat), drop = FALSE], 1, var, na.rm = TRUE)
  variance <- pmax(variance, 0)
  mu <- pmax(mu, 0)

  nb2_var_minus_mu <- pmax(variance - mu, 0)
  alpha_hat <- rep(0, length(mu))
  pos <- mu > 0
  alpha_hat[pos] <- pmax((variance[pos] - mu[pos]) / (mu[pos]^2), 0)
  alpha_mu <- alpha_hat * mu

  nb_diag <- data.frame(
    feature_id = as.character(rownames(norm_counts)),
    mu = mu,
    empirical_variance = variance,
    log_nb1 = log1p(mu),
    log_nb2 = log1p(nb2_var_minus_mu),
    nb_gap = log1p(nb2_var_minus_mu) - log1p(mu),
    log_alpha_mu = log1p(alpha_mu),
    stringsAsFactors = FALSE
  )

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

  annot_df <- annot_df %>%
    mutate(
      feature_id = as.character(feature_id),
      gene_symbol = as.character(gene_symbol),
      gene_symbol = if_else(is.na(gene_symbol), "", trimws(gene_symbol))
    ) %>%
    arrange(feature_id, desc(gene_symbol != ""), gene_symbol) %>%
    distinct(feature_id, .keep_all = TRUE) %>%
    mutate(gene_symbol = na_if(gene_symbol, ""))

  norm_counts <- norm_counts[!duplicated(norm_counts$feature_id), , drop = FALSE]
  disp_df <- disp_df[!duplicated(disp_df$feature_id), , drop = FALSE]
  nb_diag <- nb_diag[!duplicated(nb_diag$feature_id), , drop = FALSE]

  final_df <- res_df %>%
    left_join(annot_df, by = "feature_id") %>%
    left_join(nb_diag, by = "feature_id") %>%
    left_join(norm_counts, by = "feature_id") %>%
    left_join(disp_df, by = "feature_id")

  final_df$neglog10_padj        <- safe_neglog10(final_df$padj)
  final_df$neglog10_empirical_p <- safe_neglog10(final_df$empirical_p)
  final_df$neglog10_padjc       <- safe_neglog10(final_df$empirical_bh)

  final_df$dataset_name <- dataset_name
  final_df$hc_p_threshold_dataset <- hc_p_threshold_dataset
  final_df$hbfss_threshold_dataset <- hbfss_threshold_dataset
  final_df$Apeglm_L2FC <- final_df$lfc_shrunk
  final_df$pi_valueE <- final_df$HBFSS

  # Plotting booleans for the manuscript volcano panels.
  # These are intentionally not collapsed into mutually-exclusive classes.
  # Weak CNH, Strong CNH, Standard, and HBFSS are each plotted directly.
  final_df$is_weak_cnh <- !is.na(final_df$effect_class) & final_df$effect_class == "weak_effect"
  final_df$is_strong_cnh <- !is.na(final_df$effect_class) & final_df$effect_class == "strong_effect"
  final_df$is_standard <- !is.na(final_df$standard_significant) & final_df$standard_significant
  final_df$is_hbfss <- !is.na(final_df$HBFSS_significant) & final_df$HBFSS_significant

  final_df$has_valid_gene_symbol <- !is.na(final_df$gene_symbol) & grepl("[A-Za-z0-9]", trimws(final_df$gene_symbol))
  final_df$gene_symbol_plot <- ifelse(final_df$has_valid_gene_symbol, trimws(final_df$gene_symbol), NA_character_)

  preferred_cols <- c(
    "dataset_name", "feature_id", "gene_symbol", "baseMean",
    "lfc_shrunk", "Apeglm_L2FC", "regulation_direction",
    "pvalue", "padj", "empirical_p", "empirical_q", "empirical_bh",
    "pval", "qval", "padjc", "lfdr",
    "HBFSS", "pi_valueE", "hc_p_threshold_dataset", "hbfss_threshold_dataset",
    "resLA_padj", "resGA_padj", "standard_significant", "HBFSS_significant",
    "effect_class", "is_weak_cnh", "is_strong_cnh", "is_standard", "is_hbfss",
    "mu", "empirical_variance", "log_nb1", "log_nb2", "nb_gap", "log_alpha_mu"
  )

  final_df <- final_df[, c(intersect(preferred_cols, names(final_df)), setdiff(names(final_df), preferred_cols)), drop = FALSE]

  list(
    dds = dds,
    results = final_df,
    base_mean_vec = res_df$baseMean[!is.na(res_df$baseMean)],
    hc_p_threshold = hc_p_threshold_dataset,
    hbfss_threshold = hbfss_threshold_dataset
  )
}

# =============================================================================
# SECTION 6
# LABEL SELECTION AND CORE PLOT BUILDERS
# =============================================================================

select_integrated_labels <- function(df, y_col = "neglog10_empirical_p", n_labels = n_top_labels_integrated) {
  if (!nrow(df)) return(df[0, , drop = FALSE])
  df <- df[df$has_valid_gene_symbol, , drop = FALSE]
  if (!nrow(df)) return(df[0, , drop = FALSE])

  df$label_priority <- dplyr::case_when(
    df$is_hbfss & df$is_strong_cnh ~ 1,
    df$is_hbfss & df$is_standard ~ 2,
    df$is_hbfss & df$is_weak_cnh ~ 3,
    df$is_strong_cnh ~ 4,
    df$is_standard ~ 5,
    df$is_weak_cnh ~ 6,
    df$is_hbfss ~ 7,
    TRUE ~ 8
  )

  yv <- suppressWarnings(as.numeric(df[[y_col]]))
  yv[!is.finite(yv)] <- -Inf
  ord <- order(df$label_priority, -yv, -abs(df$lfc_shrunk), na.last = TRUE)
  df <- df[ord, , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol_plot), , drop = FALSE]
  df[seq_len(min(n_labels, nrow(df))), , drop = FALSE]
}

make_hbfss_curve_df <- function(x_range, hbfss_threshold, n = 500) {
  if (!is.finite(hbfss_threshold) || is.na(hbfss_threshold) || hbfss_threshold <= 0) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }
  eps <- 0.05
  left_x <- seq(x_range[1], min(-eps, x_range[2]), length.out = ceiling(n/2))
  right_x <- seq(max(eps, x_range[1]), x_range[2], length.out = ceiling(n/2))
  out <- rbind(
    data.frame(x = left_x, y = hbfss_threshold / abs(left_x)),
    data.frame(x = right_x, y = hbfss_threshold / abs(right_x))
  )
  out[is.finite(out$y), , drop = FALSE]
}

make_volcano_overlap_label <- function(df) {
  paste0(
    "Weak=", sum(df$is_weak_cnh, na.rm = TRUE),
    "  Strong=", sum(df$is_strong_cnh, na.rm = TRUE),
    "  Std=", sum(df$is_standard, na.rm = TRUE),
    "  HBFSS=", sum(df$is_hbfss, na.rm = TRUE),
    "\nW∩H=", sum(df$is_weak_cnh & df$is_hbfss, na.rm = TRUE),
    "  S∩H=", sum(df$is_strong_cnh & df$is_hbfss, na.rm = TRUE),
    "  Std∩H=", sum(df$is_standard & df$is_hbfss, na.rm = TRUE)
  )
}

plot_integrated_volcano <- function(df, panel_title) {
  lab_df <- select_integrated_labels(df, y_col = "neglog10_empirical_p", n_labels = n_top_labels_integrated)

  x_ok <- is.finite(df$lfc_shrunk)
  y_ok <- is.finite(df$neglog10_empirical_p)
  x_range <- range(df$lfc_shrunk[x_ok], na.rm = TRUE)
  y_range <- range(df$neglog10_empirical_p[y_ok], na.rm = TRUE)

  hc_p <- suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1]))
  hc_y <- if (is.finite(hc_p) && !is.na(hc_p) && hc_p > 0 && hc_p < 1) safe_neglog10(hc_p) else NA_real_

  hbfss_thr <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))
  hbfss_curve_df <- make_hbfss_curve_df(x_range, hbfss_thr)
  if (nrow(hbfss_curve_df) > 0) {
    hbfss_curve_df <- hbfss_curve_df[hbfss_curve_df$y <= (y_range[2] * 1.05), , drop = FALSE]
  }

  p <- ggplot() +
    geom_vline(
      xintercept = c(-lfc_boundary, lfc_boundary),
      linetype = "dashed", linewidth = 0.45, colour = plot_palette$threshold
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid", linewidth = 0.30, colour = "grey50"
    ) +
    geom_point(
      data = df[!(df$is_weak_cnh | df$is_strong_cnh | df$is_standard | df$is_hbfss), , drop = FALSE],
      aes(lfc_shrunk, neglog10_empirical_p),
      inherit.aes = FALSE,
      color = "grey75", shape = 16, size = 1.0, alpha = 0.55, stroke = 0
    ) +
    geom_point(
      data = df[df$is_weak_cnh, , drop = FALSE],
      aes(lfc_shrunk, neglog10_empirical_p, color = "Weak CNH", shape = "Weak CNH"),
      inherit.aes = FALSE, alpha = 0.90, size = 1.35, stroke = 0.25
    ) +
    geom_point(
      data = df[df$is_strong_cnh, , drop = FALSE],
      aes(lfc_shrunk, neglog10_empirical_p, color = "Strong CNH", shape = "Strong CNH"),
      inherit.aes = FALSE, alpha = 0.90, size = 1.35, stroke = 0.25
    ) +
    geom_point(
      data = df[df$is_standard, , drop = FALSE],
      aes(lfc_shrunk, neglog10_empirical_p, color = "Standard", shape = "Standard"),
      inherit.aes = FALSE, alpha = 0.92, size = 1.35, stroke = 0.25
    ) +
    geom_point(
      data = df[df$is_hbfss, , drop = FALSE],
      aes(lfc_shrunk, neglog10_empirical_p, color = "HBFSS", shape = "HBFSS"),
      inherit.aes = FALSE, alpha = 0.95, size = 1.55, stroke = 0.30
    ) +
    scale_color_manual(values = integrated_class_colors, breaks = integrated_class_levels, drop = FALSE, name = "Class") +
    scale_shape_manual(values = integrated_class_shapes, breaks = integrated_class_levels, drop = FALSE, name = "Class") +
    labs(
      title = panel_title,
      subtitle = "x = shrunken log2FC; y = -log10(empirical p)",
      x = "Shrunken log2FC",
      y = expression(-log[10]("empirical p"))
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(
      plot.title = element_text(size = base_theme_size, face = "bold"),
      plot.subtitle = element_text(size = base_theme_size - 2),
      legend.position = "bottom"
    )

  if (is.finite(hc_y) && !is.na(hc_y)) {
    p <- p +
      geom_hline(yintercept = hc_y, linetype = "dotted", linewidth = 0.45, colour = plot_palette$threshold) +
      annotate("text", x = x_range[1] + 0.04 * diff(x_range), y = hc_y, label = paste0("HC=", signif(hc_p, 3)), hjust = 0, vjust = -0.35, size = 2.4, colour = plot_palette$threshold)
  }

  if (nrow(hbfss_curve_df) > 0) {
    p <- p +
      geom_line(data = hbfss_curve_df, aes(x, y), inherit.aes = FALSE, colour = plot_palette$hbfss, linewidth = 0.45) +
      annotate("text", x = x_range[2] - 0.22 * diff(x_range), y = max(hbfss_curve_df$y, na.rm = TRUE), label = paste0("HBFSS=", signif(hbfss_thr, 3)), hjust = 0, vjust = -0.3, size = 2.4, colour = plot_palette$hbfss)
  }

  p <- p + annotate(
    "label",
    x = x_range[1] + 0.04 * diff(x_range),
    y = y_range[2] - 0.05 * diff(y_range),
    label = make_volcano_overlap_label(df),
    hjust = 0, vjust = 1, size = 2.2, label.size = 0.15,
    fill = adjustcolor("white", alpha.f = 0.84)
  )

  if (nrow(lab_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      aes(lfc_shrunk, neglog10_empirical_p, label = gene_symbol_plot),
      inherit.aes = FALSE,
      size = 1.8,
      seed = 1,
      max.overlaps = 20,
      force = 1.0,
      force_pull = 0.4,
      box.padding = 0.32,
      point.padding = 0.14,
      min.segment.length = 0,
      segment.alpha = 0.50,
      segment.size = 0.18
    )
  }

  p
}

plot_dispersion_panel_for_dataset <- function(df, panel_title) {
  ggplot(df, aes(baseMean, dispersion, color = plot_class, shape = plot_class)) +
    geom_point(alpha = 0.55, size = 1.15, stroke = 0.30) +
    scale_x_log10(labels = label_number(accuracy = 0.1)) +
    scale_y_log10(labels = label_number(accuracy = 0.1)) +
    scale_color_manual(
      values = integrated_class_colors,
      breaks = integrated_class_levels,
      drop = FALSE,
      name = "Class"
    ) +
    scale_shape_manual(
      values = integrated_class_shapes,
      breaks = integrated_class_levels,
      drop = FALSE,
      name = "Class"
    ) +
    labs(
      title = panel_title,
      subtitle = "Final dispersion vs mean",
      x = "baseMean",
      y = "Dispersion"
    ) +
    manuscript_theme() +
    theme(
      plot.title = element_text(size = base_theme_size, face = "bold"),
      plot.subtitle = element_text(size = base_theme_size - 2)
    )
}

plot_empirical_p_histogram <- function(df, dataset_name, hc_p_threshold) {
  subtitle_text <- if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < 1) {
    paste0("HC = ", signif(hc_p_threshold, 4))
  } else {
    "HC unavailable"
  }

  p <- ggplot(df, aes(empirical_p)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, color = "white") +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "Empirical p")),
      subtitle = subtitle_text,
      x = "Empirical p",
      y = "Count"
    ) +
    manuscript_theme()

  if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < 1) {
    p <- p + geom_vline(xintercept = hc_p_threshold, color = plot_palette$hbfss, linewidth = 0.9)
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
    geom_point(data = peak_df, aes(empirical_p, hc_score), inherit.aes = FALSE, size = 2.2, colour = plot_palette$hbfss) +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "HC profile")),
      subtitle = "Higher-criticism over sorted empirical p",
      x = "Sorted empirical p",
      y = "HC score"
    ) +
    manuscript_theme()

  if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < 1) {
    p <- p + geom_vline(xintercept = hc_p_threshold, color = plot_palette$hbfss, linewidth = 0.9)
  }

  p
}

plot_shrunken_ma <- function(df, dataset_name) {
  ggplot(df, aes(safe_log10(baseMean + 1), lfc_shrunk, color = plot_class, shape = plot_class)) +
    geom_point(alpha = 0.65, size = 1.2, stroke = 0.25) +
    geom_hline(yintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", color = plot_palette$threshold, linewidth = 0.6) +
    scale_color_manual(values = integrated_class_colors, breaks = integrated_class_levels, drop = FALSE, name = "Class") +
    scale_shape_manual(values = integrated_class_shapes, breaks = integrated_class_levels, drop = FALSE, name = "Class") +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "MA")),
      subtitle = "Shrunken effect vs mean",
      x = expression(log[10]("baseMean + 1")),
      y = "Shrunken log2FC"
    ) +
    manuscript_theme()
}

plot_pca_scatter <- function(pca_fit, dataset_label, group_label, preprocessing_label) {
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
    geom_point(size = 2.8, colour = "white", stroke = 0.50) +
    geom_text_repel(size = 1.9, max.overlaps = 8, force = 1.0, box.padding = 0.20, point.padding = 0.10, min.segment.length = 0) +
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition") +
    scale_fill_manual(values = condition_fills, labels = condition_labels, name = "Condition") +
    labs(
      title = compact_title(paste(dataset_label, pretty_group_label(group_label), "PCA")),
      subtitle = paste0(preprocessing_label, " · PC1=", pca_var_per[1], "%; PC2=", pca_var_per[2], "%"),
      x = paste0("PC1 (", pca_var_per[1], "%)"),
      y = paste0("PC2 (", pca_var_per[2], "%)")
    ) +
    manuscript_theme()
}

plot_pc1_loading_rank <- function(loading_tbl, cutoff, dataset_label, group_label, top_n_used, cutoff_quantile, preprocessing_label) {
  quantile_label <- if (is.finite(cutoff_quantile)) {
    paste0("q=", signif(cutoff_quantile, 4))
  } else {
    "q=NA"
  }

  ggplot(loading_tbl, aes(rank, pc1_loading_abs)) +
    geom_line(linewidth = 0.4, color = "grey35") +
    geom_hline(yintercept = cutoff, color = plot_palette$threshold, linewidth = 0.9) +
    annotate(
      "text",
      x = max(loading_tbl$rank) * 0.75,
      y = cutoff,
      label = paste0("top ", top_n_used, " cutoff=", signif(cutoff, 4), "\n", quantile_label),
      color = plot_palette$threshold,
      vjust = -0.8,
      size = 3.1
    ) +
    labs(
      title = compact_title(paste(dataset_label, pretty_group_label(group_label), "EVS rank")),
      subtitle = preprocessing_label,
      x = "Rank",
      y = "|PC1 loading|"
    ) +
    manuscript_theme()
}

plot_eigenvector_histograms <- function(loading_tbl, cutoff, dataset_label, group_label, preprocessing_label) {
  full_vals <- loading_tbl$pc1_loading_abs
  lead_vals <- loading_tbl$pc1_loading_abs[loading_tbl$pc1_loading_abs >= cutoff]
  rem_vals  <- loading_tbl$pc1_loading_abs[loading_tbl$pc1_loading_abs < cutoff]
  cutoff_tx <- -log(pmax(cutoff * 100, 1e-6))

  p1 <- ggplot(data.frame(x = -log(pmax(full_vals * 100, 1e-6))), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    geom_vline(xintercept = cutoff_tx, color = plot_palette$threshold, linewidth = 0.9) +
    labs(title = "All", x = "-log(|PC1|×100)", y = "Count") +
    manuscript_theme()

  p2 <- ggplot(data.frame(x = -log(pmax(rem_vals * 100, 1e-6))), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    geom_vline(xintercept = cutoff_tx, color = plot_palette$threshold, linewidth = 0.9) +
    labs(title = "Rem", x = "-log(|PC1|×100)", y = "Count") +
    manuscript_theme()

  p3 <- ggplot(data.frame(x = -log(pmax(lead_vals * 100, 1e-6))), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    geom_vline(xintercept = cutoff_tx, color = plot_palette$threshold, linewidth = 0.9) +
    labs(title = "Lead", x = "-log(|PC1|×100)", y = "Count") +
    manuscript_theme()

  arrangeGrob(
    p1, p2, p3,
    ncol = 3,
    top = textGrob(
      paste(dataset_label, pretty_group_label(group_label), "EVS loadings"),
      gp = gpar(fontface = "bold", cex = 1.05)
    )
  )
}

compute_mean_expression_table <- function(raw_counts, coldata) {
  sample_ids <- colnames(raw_counts)
  trt_ids    <- sample_ids[coldata$condition == "trt"]
  untrt_ids  <- sample_ids[coldata$condition == "untrt"]

  data.frame(
    feature_id = as.character(rownames(raw_counts)),
    mean_trt   = rowMeans(raw_counts[, trt_ids, drop = FALSE]),
    mean_untrt = rowMeans(raw_counts[, untrt_ids, drop = FALSE]),
    mean_all   = rowMeans(raw_counts),
    stringsAsFactors = FALSE
  )
}

plot_mean_histogram_panel <- function(df_means, base_mean_vec, dataset_name) {
  p1 <- ggplot(data.frame(x = safe_log10(df_means$mean_trt + 1)), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    labs(title = "Trt mean", x = "log10(mean+1)", y = "Count") +
    manuscript_theme()

  p2 <- ggplot(data.frame(x = safe_log10(df_means$mean_untrt + 1)), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    labs(title = "Ctrl mean", x = "log10(mean+1)", y = "Count") +
    manuscript_theme()

  p3 <- ggplot(data.frame(x = safe_log10(df_means$mean_all + 1)), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    labs(title = "All mean", x = "log10(mean+1)", y = "Count") +
    manuscript_theme()

  p4 <- ggplot(data.frame(x = safe_log10(base_mean_vec + 1)), aes(x)) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    labs(title = "baseMean", x = "log10(baseMean+1)", y = "Count") +
    manuscript_theme()

  arrangeGrob(
    p1, p2, p3, p4,
    ncol = 2,
    top = textGrob(
      paste(pretty_dataset_label(dataset_name), "mean histograms"),
      gp = gpar(fontface = "bold", cex = 1.05)
    )
  )
}

plot_dispersion_cloud <- function(dds, dataset_name, fig_subdir, cmp_short, ds_key) {
  outfile <- fig_file(fig_subdir, cmp_short, ds_key, "DispEst")
  png(outfile, width = 2400, height = 1900, res = figure_dpi)
  op <- par(no.readonly = TRUE)
  on.exit({par(op); dev.off()}, add = TRUE)
  par(mar = c(5.2, 5.2, 4.6, 2.2), mgp = c(2.8, 0.9, 0), cex.main = 1.0, cex.lab = 1.0)
  plotDispEsts(dds, main = paste(pretty_dataset_label(dataset_name), "dispersion"))
}

plot_dispersion_relationships <- function(df, dataset_name) {
  plots <- list()

  if (all(c("baseMean", "dispersion") %in% names(df))) {
    plots[[length(plots) + 1]] <- ggplot(df, aes(baseMean, dispersion, color = plot_class, shape = plot_class)) +
      geom_point(alpha = 0.58, size = 1.05, stroke = 0.30) +
      scale_x_log10(labels = comma_format()) +
      scale_y_log10() +
      scale_color_manual(values = integrated_class_colors, breaks = integrated_class_levels, drop = FALSE, name = "Class") +
      scale_shape_manual(values = integrated_class_shapes, breaks = integrated_class_levels, drop = FALSE, name = "Class") +
      labs(title = "Final", x = "baseMean", y = "Dispersion") +
      manuscript_theme()
  }

  if (all(c("baseMean", "dispFit") %in% names(df))) {
    plots[[length(plots) + 1]] <- ggplot(df, aes(baseMean, dispFit, color = plot_class, shape = plot_class)) +
      geom_point(alpha = 0.58, size = 1.05, stroke = 0.30) +
      scale_x_log10(labels = comma_format()) +
      scale_y_log10() +
      scale_color_manual(values = integrated_class_colors, breaks = integrated_class_levels, drop = FALSE, name = "Class") +
      scale_shape_manual(values = integrated_class_shapes, breaks = integrated_class_levels, drop = FALSE, name = "Class") +
      labs(title = "Trend", x = "baseMean", y = "dispFit") +
      manuscript_theme()
  }

  if (all(c("baseMean", "dispGeneEst") %in% names(df))) {
    plots[[length(plots) + 1]] <- ggplot(df, aes(baseMean, dispGeneEst, color = plot_class, shape = plot_class)) +
      geom_point(alpha = 0.58, size = 1.05, stroke = 0.30) +
      scale_x_log10(labels = comma_format()) +
      scale_y_log10() +
      scale_color_manual(values = integrated_class_colors, breaks = integrated_class_levels, drop = FALSE, name = "Class") +
      scale_shape_manual(values = integrated_class_shapes, breaks = integrated_class_levels, drop = FALSE, name = "Class") +
      labs(title = "Gene-wise", x = "baseMean", y = "dispGeneEst") +
      manuscript_theme()
  }

  if (length(plots) == 0) return(NULL)

  legend <- extract_shared_legend(plots[[1]])
  plots_noleg <- lapply(plots, remove_legend)

  assemble_panel_with_shared_legend(
    plot_list = plots_noleg,
    ncol = min(3, length(plots_noleg)),
    legend_grob = legend,
    top = textGrob(
      paste(pretty_dataset_label(dataset_name), "dispersion"),
      gp = gpar(fontface = "bold", cex = 1.05)
    )
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
    labs(title = "fit-final", x = "Residual", y = "Count") +
    manuscript_theme()

  p2 <- ggplot(out, aes(residual_fit_minus_gene)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, color = "white") +
    labs(title = "fit-gene", x = "Residual", y = "Count") +
    manuscript_theme()

  p3 <- ggplot(out, aes(residual_final_minus_gene)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, color = "white") +
    labs(title = "final-gene", x = "Residual", y = "Count") +
    manuscript_theme()

  panel <- arrangeGrob(
    p1, p2, p3, ncol = 3,
    top = textGrob(
      paste(pretty_dataset_label(dataset_name), "disp residuals"),
      gp = gpar(fontface = "bold", cex = 1.05)
    )
  )

  save_grob(panel, fig_file(fig_subdir, cmp_short, ds_key, "DispResid"), width = 15.5, height = 5.4)
  save_csv(out, tab_file(tab_dir, cmp_short, ds_key, "DispResid"))
  out
}

# =============================================================================
# SECTION 7
# CROSS-DATASET MANUSCRIPT FIGURES
# =============================================================================

save_cross_dataset_comparison_panels <- function(comparison_name, analysis_results, cmp_dir, dataset_list, coldata) {
  keys_present <- intersect(names(dataset_short), names(analysis_results))
  if (length(keys_present) == 0) return(invisible(NULL))

  # ---------------------------------------------------------------------------
  # Integrated volcano panel: Raw / Lead / Rem with one shared legend
  # ---------------------------------------------------------------------------
  volcano_plots <- lapply(keys_present, function(k) {
    plot_integrated_volcano(
      analysis_results[[k]]$results,
      panel_title = paste(comparison_name, dataset_short[[k]])
    )
  })
  volcano_legend <- extract_shared_legend(volcano_plots[[1]])
  volcano_panel <- assemble_panel_with_shared_legend(
    plot_list = lapply(volcano_plots, remove_legend),
    ncol = length(volcano_plots),
    legend_grob = volcano_legend,
    top = textGrob(
      paste(comparison_name, "integrated volcano"),
      gp = gpar(fontface = "bold", cex = 1.10)
    )
  )
  save_grob(
    volcano_panel,
    file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_Volcano.png")),
    width = 6.2 * length(volcano_plots),
    height = 7.2
  )

  # ---------------------------------------------------------------------------
  # Dispersion panel: Raw / Lead / Rem with one shared legend
  # ---------------------------------------------------------------------------
  disp_plots <- lapply(keys_present, function(k) {
    plot_dispersion_panel_for_dataset(
      analysis_results[[k]]$results,
      panel_title = paste(comparison_name, dataset_short[[k]])
    )
  })
  disp_legend <- extract_shared_legend(disp_plots[[1]])
  disp_panel <- assemble_panel_with_shared_legend(
    plot_list = lapply(disp_plots, remove_legend),
    ncol = length(disp_plots),
    legend_grob = disp_legend,
    top = textGrob(
      paste(comparison_name, "dispersion"),
      gp = gpar(fontface = "bold", cex = 1.10)
    )
  )
  save_grob(
    disp_panel,
    file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_Disp.png")),
    width = 6.0 * length(disp_plots),
    height = 7.0
  )

  # Optional diagnostics only
  if (isTRUE(export_optional_diagnostics)) {
    hist_grobs <- lapply(keys_present, function(k) {
      plot_empirical_p_histogram_for_panel <- plot_empirical_p_histogram(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1],
        analysis_results[[k]]$summary$hc_p_threshold[1]
      )
      plot_empirical_p_histogram_for_panel
    })

    hist_panel <- do.call(
      arrangeGrob,
      c(hist_grobs, list(
        ncol = length(hist_grobs),
        top = textGrob(
          paste(comparison_name, "empirical p"),
          gp = gpar(fontface = "bold", cex = 1.05)
        )
      ))
    )
    save_grob(hist_panel, file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_EmpP.png")),
              width = 5.6 * length(hist_grobs), height = 5.4)

    pca_norm_grobs <- lapply(keys_present, function(k) {
      compute_dataset_pca_plot(
        dataset_list[[k]],
        coldata,
        analysis_results[[k]]$summary$dataset_name[1],
        preprocessing = "normalized"
      )
    })
    pca_raw_grobs <- lapply(keys_present, function(k) {
      compute_dataset_pca_plot(
        dataset_list[[k]],
        coldata,
        analysis_results[[k]]$summary$dataset_name[1],
        preprocessing = "raw_counts"
      )
    })

    pca_panel <- arrangeGrob(
      grobs = c(pca_norm_grobs, pca_raw_grobs),
      ncol = length(keys_present),
      top = textGrob(
        paste(comparison_name, "PCA normalized vs raw"),
        gp = gpar(fontface = "bold", cex = 1.05)
      )
    )
    save_grob(pca_panel, file.path(cmp_dir, paste0("Figure_", comparison_name, "_Compare_PCA.png")),
              width = 5.8 * length(keys_present), height = 9.8)
  }

  summary_table <- bind_rows(lapply(keys_present, function(k) {
    sm <- analysis_results[[k]]$summary
    data.frame(
      comparison_name = comparison_name,
      dataset_key = k,
      dataset_label = dataset_short[[k]],
      dataset_name = sm$dataset_name[1],
      n_features = sm$n_features[1],
      n_standard_significant = sm$n_standard_significant[1],
      n_HBFSS_significant = sm$n_HBFSS_significant[1],
      n_overlap = sm$n_overlap_significant[1],
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
# SECTION 8
# FULL PIPELINE FOR ONE COMPARISON
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

  # ---------------------------------------------------------------------------
  # EVS
  # ---------------------------------------------------------------------------
  evs <- build_eigenvector_split(count_matrix, coldata)

  if (isTRUE(export_optional_diagnostics)) {
    save_plot(
      plot_pca_scatter(evs$fit_trt$pca_fit, comparison_name, "treatment", evs$fit_trt$preprocessing_label),
      file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_PCA.png")),
      width = 8.6, height = 6.8
    )
    save_plot(
      plot_pca_scatter(evs$fit_untrt$pca_fit, comparison_name, "control", evs$fit_untrt$preprocessing_label),
      file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_PCA.png")),
      width = 8.6, height = 6.8
    )

    save_plot(
      plot_pc1_loading_rank(
        evs$fit_trt$loading_table, evs$fit_trt$cutoff, comparison_name, "treatment",
        top_n_used = evs$fit_trt$top_n_used,
        cutoff_quantile = evs$fit_trt$cutoff_quantile,
        preprocessing_label = evs$fit_trt$preprocessing_label
      ),
      file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_Rank.png")),
      width = 9.2, height = 6.8
    )

    save_plot(
      plot_pc1_loading_rank(
        evs$fit_untrt$loading_table, evs$fit_untrt$cutoff, comparison_name, "control",
        top_n_used = evs$fit_untrt$top_n_used,
        cutoff_quantile = evs$fit_untrt$cutoff_quantile,
        preprocessing_label = evs$fit_untrt$preprocessing_label
      ),
      file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_Rank.png")),
      width = 9.2, height = 6.8
    )

    save_grob(
      plot_eigenvector_histograms(
        evs$fit_trt$loading_table, evs$fit_trt$cutoff, comparison_name, "treatment", evs$fit_trt$preprocessing_label
      ),
      file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_Hist.png")),
      width = 14.5, height = 4.8
    )

    save_grob(
      plot_eigenvector_histograms(
        evs$fit_untrt$loading_table, evs$fit_untrt$cutoff, comparison_name, "control", evs$fit_untrt$preprocessing_label
      ),
      file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_Hist.png")),
      width = 14.5, height = 4.8
    )
  }

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

    # Tables
    save_csv(df, tab_file(tab_dir, comparison_name, nm, "Results"))
    save_csv(subset(df, standard_significant), tab_file(tab_dir, comparison_name, nm, "Std"))
    save_csv(subset(df, HBFSS_significant), tab_file(tab_dir, comparison_name, nm, "HBFSS"))
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

    # Mean histograms
    if (isTRUE(export_optional_diagnostics)) {
      mean_df <- compute_mean_expression_table(dataset_list[[nm]], coldata)
      mean_panel <- plot_mean_histogram_panel(mean_df, fit$base_mean_vec, full_dataset_name)
      save_grob(mean_panel, fig_file(fig_subdir, comparison_name, nm, "MeanHist"), width = 12.5, height = 9.2)
    }

    # Single-dataset figures
    save_plot(
      plot_integrated_volcano(df, panel_title = paste(comparison_name, dataset_short[[nm]])),
      fig_file(fig_subdir, comparison_name, nm, "Volcano"),
      width = 7.2, height = 6.6
    )

    save_plot(
      plot_dispersion_panel_for_dataset(df, panel_title = paste(comparison_name, dataset_short[[nm]])),
      fig_file(fig_subdir, comparison_name, nm, "Disp"),
      width = 7.2, height = 6.6
    )

    if (isTRUE(export_optional_diagnostics)) {
      save_plot(plot_empirical_p_histogram(df, full_dataset_name, fit$hc_p_threshold),
                fig_file(fig_subdir, comparison_name, nm, "EmpP"),
                width = 7.0, height = 5.8)

      p_hc <- plot_hc_profile(df, full_dataset_name, fit$hc_p_threshold)
      if (!is.null(p_hc)) {
        save_plot(p_hc, fig_file(fig_subdir, comparison_name, nm, "HC"), width = 7.0, height = 5.8)
      }

      save_plot(plot_shrunken_ma(df, full_dataset_name),
                fig_file(fig_subdir, comparison_name, nm, "MA"),
                width = 7.0, height = 5.8)
    }

    # Dispersion diagnostics
    if (isTRUE(export_optional_diagnostics)) {
      plot_dispersion_cloud(fit$dds, full_dataset_name, fig_subdir, comparison_name, nm)

      disp_panel <- plot_dispersion_relationships(df, full_dataset_name)
      if (!is.null(disp_panel)) {
        save_grob(disp_panel, fig_file(fig_subdir, comparison_name, nm, "DispRel"), width = 13.2, height = 7.2)
      }

      dispersion_residual_section(df, full_dataset_name, fig_subdir, tab_dir, comparison_name, nm)
    }

    analysis_results[[nm]] <- list(
      dds = fit$dds,
      results = df,
      summary = summary_row,
      dataset_mat = dataset_list[[nm]],
      fig_subdir = fig_subdir
    )
  }

  save_cross_dataset_comparison_panels(comparison_name, analysis_results, cmp_dir, dataset_list, coldata)

  # Optional TWAS overlap
  if (run_twas_overlap) {
    twas_genes <- clean_gene_set(TWAS_data$gene_symbol)

    get_twas_overlap <- function(result_df, dataset_nm, cmp_name, out_dir,
                                 mode = c("union", "standard_only", "HBFSS_only")) {
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

      save_csv(overlap_df, file.path(out_dir, paste0("Table_", dataset_nm, "_TWAS_", mode_tag, ".csv")))
      save_csv(summary_df, file.path(out_dir, paste0("Table_", dataset_nm, "_TWAS_", mode_tag, "_Summary.csv")))

      summary_df
    }

    twas_summaries <- bind_rows(lapply(names(analysis_results), function(nm) {
      full_nm <- paste(comparison_name, nm, sep = "_")
      bind_rows(
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
  bind_rows(sm_list)
}

# =============================================================================
# SECTION 9
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
  failed_df <- bind_rows(failed_comparisons)
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
