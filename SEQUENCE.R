#!/usr/bin/env Rscript

# =============================================================================
# FINAL MANUSCRIPT SCRIPT
# EVS + DESeq2 + empirical-null HC + HBFSS
#
# This version is intentionally rewritten to do the following cleanly:
#
# 1. Read one WTTS master count file from either:
#      - data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv
#      - WTTS-Seq_2022.2_DE_raw_read_numbers.csv
#
# 2. Embed the minimal metadata needed for the four fixed pairwise comparisons:
#      RT0_ZT6, RT2_ZT8, RT4_ZT10, RT8_ZT14
#
# 3. Build the EVS split using DESeq2-normalized counts prior to PCA.
#    The split is defined by the union of the top-N absolute PC1 loading PAS
#    features from treatment and control separately.
#
# 4. Run DESeq2 independently on:
#      - the original dataset
#      - the leading-edge dataset
#      - the remainder dataset
#
# 5. Compute:
#      - native DESeq2 Wald p-values and BH-adjusted p-values
#      - empirical-null p-values, q-values, and local FDR from fdrtool
#      - higher-criticism threshold on sorted empirical-null p-values
#      - apeglm-shrunken log2 fold changes
#      - HBFSS = | lfc_shrunk * log10(empirical_p) |
#
# 6. Produce short, GitHub-friendly export names inside:
#      exports/manuscript_final/
#
# 7. Produce manuscript-focused volcano plots with only these plotted classes:
#      - Weak CNH
#      - Strong CNH
#      - Standard
#      - HBFSS
#
#    IMPORTANT CLASS RULE:
#    No point below the dataset-specific HC threshold is colored.
#    Anything below HC remains Background.
#
# 8. Use one legend per panel.
#
# 9. Keep titles, subtitles, labels, and captions compact enough to fit.
#
# -----------------------------------------------------------------------------
# METHODS SUMMARY
#
# EVS step:
# The EVS split is based on DESeq2-normalized counts. Within each comparison,
# treatment and control are analyzed separately by PCA on the normalized count
# matrix. Absolute PC1 loadings are ranked within each condition. The union of
# the top-N PAS features from the two condition-specific ranked loading tables
# defines the leading-edge dataset. Features not in that union define the
# remainder dataset.
#
# DESeq2 step:
# Each dataset (Original, Lead, Rem) is analyzed as its own DESeq2 object with
# design = ~ condition. No donor blocking term is used in this final version.
#
# Composite null hypotheses:
#   Weak CNH   : lessAbs test significant, |lfc_shrunk| < lfc_boundary
#   Strong CNH : greaterAbs test significant, |lfc_shrunk| >= lfc_boundary
#
# Classical DESeq2 plotted class:
#   Standard   : native padj < alpha, |raw LFC| >= lfc_boundary,
#                |lfc_shrunk| >= lfc_boundary
#
# Empirical-null and HBFSS:
# Finite DESeq2 Wald statistics are passed to fdrtool under the normal
# empirical-null model. The dataset-specific higher-criticism threshold is
# obtained from hc.thresh(sort(empirical_p)).
#
# Let:
#   hc_p_threshold_dataset = HC threshold on empirical-null p-values
#   hbfss_threshold_dataset = abs(log10(hc_p_threshold_dataset)) * lfc_boundary
#
# Then:
#   HBFSS = abs(lfc_shrunk * log10(empirical_p))
#
# Volcano interpretation:
# x-axis = shrunken log2FC
# y-axis = -log10(empirical_p)
#
# The horizontal HC line marks the minimum empirical evidence floor used for
# coloring. No point below that line is assigned Weak CNH, Strong CNH, Standard,
# or HBFSS. The HBFSS curve is:
#
#   y = hbfss_threshold_dataset / |x|
#
# so HBFSS-positive points lie above the curve and above the HC line.
#
# Class plotting precedence:
#   HBFSS > Standard > Strong CNH > Weak CNH > Background
#
# Overlap counts are written in compact subtitle/caption text rather than using
# extra legend entries.
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
# SETTINGS
# =============================================================================

count_file_candidates <- c(
  "data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  "WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  "/root/REAPER98632/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
)

run_twas_overlap <- FALSE
twas_file_candidates <- c(
  "data/3aTWAS_genes_of_11_brain_disorders.csv",
  "3aTWAS_genes_of_11_brain_disorders.csv"
)

alpha_level   <- 0.10
lfc_boundary  <- 1.0
top_n_target  <- 5000L

figure_dpi            <- 320
base_theme_size       <- 10
n_top_labels_standard <- 12L
n_top_labels_hbfss    <- 12L

output_dir <- "exports/manuscript_final"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plot_palette <- list(
  background   = "#BDBDBD",
  weak_cnh     = "#1F78B4",
  strong_cnh   = "#E66101",
  standard     = "#33A02C",
  hbfss        = "#6A3D9A",
  threshold    = "#8C2D04",
  histogram    = "#969696",
  control      = "#4D4D4D",
  treatment    = "#1F78B4"
)

dataset_short <- c(
  raw_dataset = "Raw",
  leading_edge_dataset = "Lead",
  remainder_dataset = "Rem"
)

dataset_label <- c(
  raw_dataset = "Original",
  leading_edge_dataset = "Lead",
  remainder_dataset = "Rem"
)

class_levels <- c("Background", "Weak CNH", "Strong CNH", "Standard", "HBFSS")
class_colors <- c(
  "Background" = plot_palette$background,
  "Weak CNH"   = plot_palette$weak_cnh,
  "Strong CNH" = plot_palette$strong_cnh,
  "Standard"   = plot_palette$standard,
  "HBFSS"      = plot_palette$hbfss
)
class_shapes <- c(
  "Background" = 16,
  "Weak CNH"   = 16,
  "Strong CNH" = 17,
  "Standard"   = 15,
  "HBFSS"      = 18
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

# =============================================================================
# EMBEDDED METADATA
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

# =============================================================================
# FILE HELPERS
# =============================================================================

detect_first_existing_file <- function(candidates) {
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop(
      "None of the expected input files were found.\nChecked:\n",
      paste(" -", candidates, collapse = "\n")
    )
  }
  found[1]
}

fig_file <- function(dir, cmp, ds_key, tag) {
  file.path(dir, paste0("Figure_", cmp, "_", unname(dataset_short[ds_key]), "_", tag, ".png"))
}

tab_file <- function(dir, cmp, ds_key, tag) {
  file.path(dir, paste0("Table_", cmp, "_", unname(dataset_short[ds_key]), "_", tag, ".csv"))
}

save_csv <- function(df, path) {
  write.csv(df, file = path, row.names = FALSE)
}

save_plot <- function(p, path, width = 10, height = 7, dpi = figure_dpi, bg = "white") {
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

save_grob <- function(g, path, width = 12, height = 8, dpi = figure_dpi, bg = "white") {
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

# =============================================================================
# GENERAL HELPERS
# =============================================================================

assert_required_columns <- function(df, required_cols, object_name = "data frame") {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", object_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

safe_log10 <- function(x, pseudocount = 1e-12) {
  log10(pmax(x, pseudocount))
}

safe_neglog10 <- function(x, pseudocount = 1e-300) {
  -log10(pmax(x, pseudocount))
}

compact_title <- function(x, width = 44) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

compact_caption <- function(x, width = 110) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

clip_probabilities <- function(x, eps = 1e-300) {
  x <- as.numeric(x)
  bad <- !is.finite(x) | is.na(x)
  x[bad] <- NA_real_
  ok <- !is.na(x)
  x[ok] <- pmin(pmax(x[ok], eps), 1 - 1e-12)
  x
}

pretty_dataset_label <- function(dataset_name) {
  parts <- strsplit(dataset_name, "_", fixed = TRUE)[[1]]
  if (length(parts) < 4) return(dataset_name)
  cmp <- paste(parts[1], parts[2], sep = "_")
  ds  <- paste(parts[3:length(parts)], collapse = "_")
  paste(cmp, unname(dataset_label[ds]))
}

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = base_theme_size + 1),
      plot.subtitle    = element_text(hjust = 0.5, size = base_theme_size - 1),
      plot.caption     = element_text(hjust = 0.5, size = base_theme_size - 2, colour = "grey25"),
      axis.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold"),
      legend.box       = "horizontal",
      plot.margin      = margin(10, 14, 10, 14)
    )
}

get_condition_coef <- function(dds) {
  rn <- resultsNames(dds)
  idx <- grep("^condition_", rn)
  if (length(idx) == 0) stop("Could not identify condition coefficient in resultsNames(dds).")
  rn[idx[1]]
}

classify_effect_strength <- function(res_strong_padj, res_weak_padj, alpha = alpha_level) {
  out <- rep("intermediate", length(res_strong_padj))
  out[!is.na(res_weak_padj)   & res_weak_padj   < alpha] <- "weak_effect"
  out[!is.na(res_strong_padj) & res_strong_padj < alpha] <- "strong_effect"
  out
}

resolve_top_n_cutoff <- function(sorted_values_desc, top_n = top_n_target) {
  n_total <- length(sorted_values_desc)
  if (n_total == 0L) stop("resolve_top_n_cutoff() received an empty vector.")
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

extract_legend_grob <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx[1]]]
}

panel_with_shared_legend <- function(plot_list, legend_grob, ncol = 3, top_text = NULL) {
  plots_no_legend <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
  panel <- arrangeGrob(grobs = plots_no_legend, ncol = ncol)

  pieces <- list()
  if (!is.null(top_text)) {
    pieces[[length(pieces) + 1L]] <- textGrob(
      top_text,
      gp = gpar(fontface = "bold", cex = 1.05)
    )
  }
  pieces[[length(pieces) + 1L]] <- panel
  if (!is.null(legend_grob)) {
    pieces[[length(pieces) + 1L]] <- legend_grob
  }

  heights <- c()
  if (!is.null(top_text)) heights <- c(0.08)
  heights <- c(heights, 1)
  if (!is.null(legend_grob)) heights <- c(heights, 0.12)

  do.call(arrangeGrob, c(pieces, ncol = 1, heights = heights))
}

# =============================================================================
# INPUT
# =============================================================================

count_file <- detect_first_existing_file(count_file_candidates)
message("Using count file: ", normalizePath(count_file))

WTTS_Seq <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
WTTS_Seq <- as.data.frame(WTTS_Seq, stringsAsFactors = FALSE)

assert_required_columns(WTTS_Seq, c("OrigID", "Symbol"), "WTTS count file")
assert_required_columns(WTTS_Seq, meta_all$id, "WTTS count file sample columns")

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
  twas_file <- detect_first_existing_file(twas_file_candidates)
  TWAS_Seq <- read.csv(twas_file, header = TRUE, stringsAsFactors = FALSE)
  if (ncol(TWAS_Seq) < 4) stop("TWAS file must contain at least 4 columns.")
  TWAS_data <- TWAS_Seq[, c(1, 4), drop = FALSE]
  colnames(TWAS_data) <- c("source_id", "gene_symbol")
}

# =============================================================================
# COMPARISON PREP
# =============================================================================

prepare_comparison_data <- function(comparison_name, group1_prefix, group2_prefix, WTTS_Seq, meta_all) {
  keep_ids <- grepl(paste0("^", group1_prefix, "_"), meta_all$id) |
    grepl(paste0("^", group2_prefix, "_"), meta_all$id)

  meta_sub <- meta_all[keep_ids, , drop = FALSE]
  coldata  <- meta_sub[, c("condition"), drop = FALSE]
  sample_ids <- rownames(meta_sub)

  missing_samples <- setdiff(sample_ids, colnames(WTTS_Seq))
  if (length(missing_samples) > 0) {
    stop(
      "Missing samples in WTTS file for ", comparison_name, ": ",
      paste(missing_samples, collapse = ", ")
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
    pca_fit         = pca_fit,
    loading_table   = loading_tbl,
    cutoff          = cutoff,
    top_n_used      = cutoff_info$top_n_actual,
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
  raw_counts_init  <- as.data.frame(count_matrix)

  sample_ids <- colnames(count_matrix)
  trt_ids    <- sample_ids[coldata$condition == "trt"]
  untrt_ids  <- sample_ids[coldata$condition == "untrt"]

  fit_trt <- compute_pc1_loading_table(norm_counts_init, trt_ids, top_n = top_n_target)
  fit_untrt <- compute_pc1_loading_table(norm_counts_init, untrt_ids, top_n = top_n_target)

  trt_high   <- as.character(subset(fit_trt$loading_table, split_class == "high_loading")$feature_id)
  untrt_high <- as.character(subset(fit_untrt$loading_table, split_class == "high_loading")$feature_id)

  leading_edge_ids <- union(trt_high, untrt_high)
  remainder_ids <- setdiff(rownames(count_matrix), leading_edge_ids)

  if (length(leading_edge_ids) == 0L) {
    stop("Leading-edge dataset is empty. Check sample mapping or top_n_target.")
  }

  list(
    fit_trt              = fit_trt,
    fit_untrt            = fit_untrt,
    normalized_counts    = norm_counts_init,
    raw_counts           = raw_counts_init,
    raw_dataset          = count_matrix,
    leading_edge_dataset = count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset    = count_matrix[remainder_ids, , drop = FALSE]
  )
}

# =============================================================================
# EMPIRICAL NULL / HC
# =============================================================================

run_empirical_null_fdrtool <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  stat_vec <- unname(stat_vec)

  if (length(stat_vec) < 5L) {
    stop("[", dataset_name, "] Fewer than 5 finite Wald statistics were available for fdrtool.")
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
      message("[", dataset_name, "] Primary fdrtool call failed: ", conditionMessage(e1))
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
          stop("[", dataset_name, "] fdrtool failed after retry: ", conditionMessage(e2))
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
  if (length(sorted_empirical_p) < 5L) return(NA_real_)

  out <- suppressWarnings(
    tryCatch(
      fdrtool::hc.thresh(as.vector(sorted_empirical_p)),
      error = function(e) {
        message("[", dataset_name, "] hc.thresh failed: ", conditionMessage(e))
        NA_real_
      }
    )
  )

  out <- as.numeric(out[1])
  if (!is.finite(out) || is.na(out) || out <= 0 || out >= 1) return(NA_real_)
  out
}

# =============================================================================
# CORE ANALYSIS
# =============================================================================

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = coldata,
    design = ~ condition
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
  res_all_df$feature_id <- rownames(res_all_df)

  valid_stat <- is.finite(res_all_df$stat) & !is.na(res_all_df$stat)
  stat_vec <- as.numeric(res_all_df$stat[valid_stat])

  if (length(stat_vec) < 5L) {
    stop("[", dataset_name, "] Fewer than 5 finite Wald statistics were available for fdrtool.")
  }

  message(sprintf("[%s] Mean Wald stat sent to fdrtool: %.4f", dataset_name, mean(stat_vec, na.rm = TRUE)))
  fdr_fit <- run_empirical_null_fdrtool(stat_vec, dataset_name)

  res_df <- res_all_df
  n_valid <- sum(valid_stat)
  if (length(fdr_fit$pval) != n_valid) {
    stop("[", dataset_name, "] fdrtool output length mismatch.")
  }

  res_df$empirical_p  <- NA_real_
  res_df$empirical_q  <- NA_real_
  res_df$lfdr         <- NA_real_
  res_df$empirical_p[valid_stat] <- as.numeric(fdr_fit$pval)
  res_df$empirical_q[valid_stat] <- as.numeric(fdr_fit$qval)
  res_df$lfdr[valid_stat]        <- as.numeric(fdr_fit$lfdr)

  res_df$empirical_bh <- NA_real_
  ok_emp <- is.finite(res_df$empirical_p) & !is.na(res_df$empirical_p)
  if (any(ok_emp)) {
    res_df$empirical_bh[ok_emp] <- p.adjust(res_df$empirical_p[ok_emp], method = "BH")
  }

  res_df$pval  <- res_df$empirical_p
  res_df$qval  <- res_df$empirical_q
  res_df$padjc <- res_df$empirical_bh

  coef_name <- get_condition_coef(dds)
  shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)
  shr_df <- as.data.frame(shr)
  shr_df$feature_id <- rownames(shr_df)

  res_df <- left_join(
    res_df,
    shr_df[, c("feature_id", "log2FoldChange")],
    by = "feature_id",
    suffix = c("", "_shrunk")
  )
  colnames(res_df)[colnames(res_df) == "log2FoldChange_shrunk"] <- "lfc_shrunk"

  hc_p_threshold_dataset <- safe_hc_thresh(res_df$empirical_p, dataset_name)
  res_df$gene_empirical_pvalue <- res_df$empirical_p
  res_df$HBFSS <- abs(res_df$lfc_shrunk * log10(pmax(res_df$gene_empirical_pvalue, 1e-300)))

  if (is.na(hc_p_threshold_dataset) || hc_p_threshold_dataset <= 0 || hc_p_threshold_dataset >= 1) {
    hbfss_threshold_dataset  <- NA_real_
    res_df$HBFSS_core_pass   <- FALSE
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

  res_df$raw_lfc_pass <- !is.na(res_df$log2FoldChange) &
    (abs(res_df$log2FoldChange) >= lfc_boundary)

  res_df$shrunk_lfc_pass <- !is.na(res_df$lfc_shrunk) &
    (abs(res_df$lfc_shrunk) >= lfc_boundary)

  res_df$standard_significant <- !is.na(res_df$padj) &
    (res_df$padj < alpha_level) &
    res_df$raw_lfc_pass &
    res_df$shrunk_lfc_pass

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
    res_df$padj_weak_effect,
    alpha = alpha_level
  )

  res_df$resGA_padj <- res_df$padj_strong_effect
  res_df$resLA_padj <- res_df$padj_weak_effect

  hbfss_native_call <- (
    res_df$HBFSS_core_pass &
      (is.na(res_df$resLA_padj) | !(res_df$resLA_padj < 0.2))
  ) | (
    !is.na(res_df$resGA_padj) &
      (res_df$resGA_padj < alpha_level)
  )
  res_df$HBFSS_significant <- hbfss_native_call & res_df$raw_lfc_pass

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
  if ("dispersion" %in% names(disp_df) && "baseMean" %in% names(disp_df)) {
    disp_df$IOD <- disp_df$dispersion * disp_df$baseMean
  } else {
    disp_df$IOD <- NA_real_
  }
  if ("baseMean" %in% names(disp_df)) {
    disp_df <- disp_df[, setdiff(names(disp_df), "baseMean"), drop = FALSE]
  }

  annot_df$feature_id  <- as.character(annot_df$feature_id)
  annot_df$gene_symbol <- as.character(annot_df$gene_symbol)
  annot_df <- annot_df %>%
    mutate(gene_symbol = if_else(is.na(gene_symbol), "", trimws(gene_symbol))) %>%
    arrange(feature_id, desc(gene_symbol != ""), gene_symbol) %>%
    distinct(feature_id, .keep_all = TRUE) %>%
    mutate(gene_symbol = na_if(gene_symbol, ""))

  norm_counts <- norm_counts[!duplicated(norm_counts$feature_id), , drop = FALSE]
  disp_df <- disp_df[!duplicated(disp_df$feature_id), , drop = FALSE]

  final_df <- res_df %>%
    left_join(annot_df, by = "feature_id") %>%
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

  preferred_cols <- c(
    "dataset_name", "feature_id", "gene_symbol", "baseMean",
    "log2FoldChange", "lfc_shrunk", "Apeglm_L2FC",
    "pvalue", "padj", "empirical_p", "empirical_q", "empirical_bh",
    "pval", "qval", "padjc", "lfdr",
    "HBFSS", "pi_valueE",
    "hc_p_threshold_dataset", "hbfss_threshold_dataset",
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

# =============================================================================
# CLASS LAYER BUILDER
# =============================================================================

build_method_layers <- function(df, alpha = alpha_level, lfc_cut = lfc_boundary) {
  out <- as.data.frame(df)

  hc_p <- suppressWarnings(as.numeric(out$hc_p_threshold_dataset[1]))
  hc_y <- if (is.finite(hc_p) && !is.na(hc_p) && hc_p > 0 && hc_p < 1) {
    safe_neglog10(hc_p)
  } else {
    NA_real_
  }

  out$above_hc <- if (is.finite(hc_y) && !is.na(hc_y)) {
    !is.na(out$neglog10_empirical_p) & (out$neglog10_empirical_p >= hc_y)
  } else {
    FALSE
  }

  out$call_weak_cnh <- out$above_hc &
    !is.na(out$resLA_padj) &
    (out$resLA_padj < alpha) &
    !is.na(out$lfc_shrunk) &
    (abs(out$lfc_shrunk) < lfc_cut)

  out$call_strong_cnh <- out$above_hc &
    !is.na(out$resGA_padj) &
    (out$resGA_padj < alpha) &
    !is.na(out$lfc_shrunk) &
    (abs(out$lfc_shrunk) >= lfc_cut)

  out$call_standard <- out$above_hc &
    !is.na(out$padj) &
    (out$padj < alpha) &
    !is.na(out$log2FoldChange) &
    (abs(out$log2FoldChange) >= lfc_cut) &
    !is.na(out$lfc_shrunk) &
    (abs(out$lfc_shrunk) >= lfc_cut)

  out$call_hbfss <- out$above_hc &
    !is.na(out$HBFSS_significant) &
    out$HBFSS_significant

  out$class <- "Background"
  out$class[out$call_weak_cnh]   <- "Weak CNH"
  out$class[out$call_strong_cnh] <- "Strong CNH"
  out$class[out$call_standard]   <- "Standard"
  out$class[out$call_hbfss]      <- "HBFSS"
  out$class <- factor(out$class, levels = class_levels)

  out$gene_symbol_plot <- ifelse(
    !is.na(out$gene_symbol) & grepl("[A-Za-z0-9]", trimws(out$gene_symbol)),
    trimws(out$gene_symbol),
    NA_character_
  )

  out$hc_y <- hc_y
  out$overlap_std_hbfss <- out$call_standard & out$call_hbfss

  out
}

select_volcano_labels <- function(df, y_col, n_labels = 12L) {
  lab_df <- df[
    !is.na(df$gene_symbol_plot) &
      df$class %in% c("Standard", "HBFSS"),
    ,
    drop = FALSE
  ]
  if (!nrow(lab_df)) return(lab_df)

  metric <- suppressWarnings(as.numeric(lab_df[[y_col]]))
  metric[!is.finite(metric)] <- -Inf
  ord <- order(-metric, -abs(lab_df$lfc_shrunk), na.last = TRUE)
  lab_df <- lab_df[ord, , drop = FALSE]
  lab_df <- lab_df[!duplicated(lab_df$gene_symbol_plot), , drop = FALSE]
  lab_df[seq_len(min(n_labels, nrow(lab_df))), , drop = FALSE]
}

volcano_label_layer <- function(lab_df) {
  if (!nrow(lab_df)) return(NULL)
  geom_text_repel(
    data = lab_df,
    aes(label = gene_symbol_plot),
    size = 1.7,
    seed = 1,
    max.overlaps = 20,
    force = 1.0,
    box.padding = 0.26,
    point.padding = 0.12,
    min.segment.length = 0,
    segment.alpha = 0.55,
    segment.size = 0.20
  )
}

# =============================================================================
# VOLCANO PLOTS
# =============================================================================

build_hbfss_curve_df <- function(df, x_limit = NULL) {
  thr <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))
  if (!is.finite(thr) || is.na(thr) || thr <= 0) return(NULL)

  x_max <- if (is.null(x_limit)) {
    max(1.25, max(abs(df$lfc_shrunk), na.rm = TRUE), lfc_boundary * 1.5)
  } else {
    x_limit
  }

  x_left  <- seq(-x_max, -0.05, length.out = 300)
  x_right <- seq(0.05, x_max, length.out = 300)
  x <- c(x_left, x_right)
  y <- thr / abs(x)

  data.frame(x = x, y = y, stringsAsFactors = FALSE)
}

volcano_caption_text <- function(df) {
  paste0(
    "Weak=", sum(df$class == "Weak CNH", na.rm = TRUE),
    "  Strong=", sum(df$class == "Strong CNH", na.rm = TRUE),
    "  Std=", sum(df$call_standard, na.rm = TRUE),
    "  HBFSS=", sum(df$call_hbfss, na.rm = TRUE),
    "  Overlap=", sum(df$overlap_std_hbfss, na.rm = TRUE)
  )
}

plot_volcano_single <- function(df, dataset_name) {
  d <- build_method_layers(df)
  lab_df <- select_volcano_labels(d, "neglog10_empirical_p", n_labels = n_top_labels_hbfss)

  hc_y <- d$hc_y[1]
  hc_p <- suppressWarnings(as.numeric(d$hc_p_threshold_dataset[1]))
  hbfss_thr <- suppressWarnings(as.numeric(d$hbfss_threshold_dataset[1]))

  x_lim <- max(1.5, quantile(abs(d$lfc_shrunk), 0.995, na.rm = TRUE))
  x_lim <- min(max(x_lim, 2.5), 8)
  curve_df <- build_hbfss_curve_df(d, x_limit = x_lim)

  subtitle_txt <- paste0(
    "HC=", ifelse(is.finite(hc_p) && !is.na(hc_p), signif(hc_p, 3), "NA"),
    "  HBFSS=", ifelse(is.finite(hbfss_thr) && !is.na(hbfss_thr), signif(hbfss_thr, 3), "NA"),
    "  Std=", sum(d$call_standard, na.rm = TRUE),
    "  HBFSS=", sum(d$call_hbfss, na.rm = TRUE),
    "  Overlap=", sum(d$overlap_std_hbfss, na.rm = TRUE)
  )

  p <- ggplot(d, aes(lfc_shrunk, neglog10_empirical_p)) +
    geom_point(aes(color = class, shape = class), alpha = 0.80, size = 1.2, stroke = 0.35) +
    scale_color_manual(values = class_colors, drop = FALSE, limits = class_levels, name = "Class") +
    scale_shape_manual(values = class_shapes, drop = FALSE, limits = class_levels, name = "Class") +
    geom_vline(xintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", linewidth = 0.5, colour = plot_palette$threshold) +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.30, colour = "grey50") +
    labs(
      title = compact_title(pretty_dataset_label(dataset_name), 28),
      subtitle = subtitle_txt,
      x = "Shrunken log2FC",
      y = expression(-log[10]("Empirical p")),
      caption = volcano_caption_text(d)
    ) +
    coord_cartesian(xlim = c(-x_lim, x_lim), clip = "off") +
    manuscript_theme()

  if (is.finite(hc_y) && !is.na(hc_y)) {
    p <- p +
      geom_hline(yintercept = hc_y, linetype = "dotted", linewidth = 0.55, colour = plot_palette$threshold) +
      annotate("text", x = -x_lim * 0.85, y = hc_y, label = paste0("HC=", signif(hc_p, 3)),
               hjust = 0, vjust = -0.3, size = 2.8, colour = plot_palette$threshold)
  }

  if (!is.null(curve_df)) {
    p <- p +
      geom_line(data = curve_df, aes(x, y), inherit.aes = FALSE, linewidth = 0.7, colour = plot_palette$hbfss) +
      annotate("text", x = x_lim * 0.55, y = max(curve_df$y[curve_df$x > 0 & curve_df$x < x_lim * 0.7], na.rm = TRUE),
               label = paste0("HBFSS=", signif(hbfss_thr, 3)),
               hjust = 0, vjust = -0.1, size = 2.8, colour = plot_palette$hbfss)
  }

  if (nrow(lab_df) > 0) {
    p <- p + volcano_label_layer(lab_df)
  }

  p
}

plot_volcano_compare_panel <- function(plot_list, comparison_name) {
  legend_plot <- plot_list[[1]] + theme(legend.position = "bottom")
  legend_grob <- extract_legend_grob(legend_plot)
  panel_with_shared_legend(
    plot_list,
    legend_grob,
    ncol = 3,
    top_text = paste(comparison_name, "volcano")
  )
}

# =============================================================================
# DISPERSION PLOTS
# =============================================================================

plot_dispersion_single <- function(df, dataset_name) {
  d <- build_method_layers(df)

  ggplot(d, aes(baseMean, dispersion)) +
    geom_point(aes(color = class, shape = class), alpha = 0.65, size = 1.0, stroke = 0.30) +
    scale_color_manual(values = class_colors, drop = FALSE, limits = class_levels, name = "Class") +
    scale_shape_manual(values = class_shapes, drop = FALSE, limits = class_levels, name = "Class") +
    scale_x_log10(labels = comma_format()) +
    scale_y_log10(labels = label_number(accuracy = 0.1)) +
    labs(
      title = compact_title(pretty_dataset_label(dataset_name), 28),
      subtitle = "Final dispersion vs mean",
      x = "baseMean",
      y = "Dispersion"
    ) +
    manuscript_theme()
}

plot_dispersion_compare_panel <- function(plot_list, comparison_name) {
  legend_plot <- plot_list[[1]] + theme(legend.position = "bottom")
  legend_grob <- extract_legend_grob(legend_plot)
  panel_with_shared_legend(
    plot_list,
    legend_grob,
    ncol = 3,
    top_text = paste(comparison_name, "dispersion")
  )
}

# =============================================================================
# SUPPORT PLOTS
# =============================================================================

plot_empirical_p_histogram <- function(df, dataset_name, hc_p_threshold) {
  p <- ggplot(df, aes(empirical_p)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, color = "white") +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "Empirical p"), 32),
      subtitle = if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold)) {
        paste0("HC=", signif(hc_p_threshold, 3))
      } else {
        "HC unavailable"
      },
      x = "Empirical p",
      y = "Count"
    ) +
    manuscript_theme()

  if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < 1) {
    p <- p + geom_vline(xintercept = hc_p_threshold, color = plot_palette$threshold, linewidth = 0.8)
  }
  p
}

plot_hc_profile <- function(df, dataset_name, hc_p_threshold) {
  pvals <- sort(df$empirical_p[is.finite(df$empirical_p) & !is.na(df$empirical_p) & df$empirical_p > 0 & df$empirical_p < 1])
  if (length(pvals) < 5L) return(NULL)

  n <- length(pvals)
  i <- seq_len(n)
  v <- (i / n) * (1 - (i / n)) / n
  v[v == 0] <- min(v[v > 0])
  hc_score <- abs((i / n) - pvals) / sqrt(v)
  hc_df <- data.frame(rank = i, empirical_p = pvals, hc_score = hc_score, stringsAsFactors = FALSE)
  idx_peak <- which.max(hc_df$hc_score)
  peak_df <- hc_df[idx_peak, , drop = FALSE]

  p <- ggplot(hc_df, aes(empirical_p, hc_score)) +
    geom_line(linewidth = 0.45, colour = "grey30") +
    geom_point(data = peak_df, aes(empirical_p, hc_score), inherit.aes = FALSE, size = 2.2, colour = plot_palette$hbfss) +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "HC"), 32),
      subtitle = if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold)) {
        paste0("HC threshold=", signif(hc_p_threshold, 3))
      } else {
        "HC unavailable"
      },
      x = "Empirical p",
      y = "HC score"
    ) +
    manuscript_theme()

  if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < 1) {
    p <- p + geom_vline(xintercept = hc_p_threshold, color = plot_palette$threshold, linewidth = 0.8)
  }
  p
}

plot_hbfss_distribution <- function(df, dataset_name, hbfss_threshold) {
  p <- ggplot(df, aes(HBFSS)) +
    geom_histogram(bins = 70, fill = plot_palette$histogram, color = "white") +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "HBFSS"), 32),
      subtitle = if (is.finite(hbfss_threshold) && !is.na(hbfss_threshold)) {
        paste0("Cutoff=", signif(hbfss_threshold, 3))
      } else {
        "Cutoff unavailable"
      },
      x = "HBFSS",
      y = "Count"
    ) +
    manuscript_theme()

  if (is.finite(hbfss_threshold) && !is.na(hbfss_threshold)) {
    p <- p + geom_vline(xintercept = abs(hbfss_threshold), color = plot_palette$hbfss, linewidth = 0.8)
  }
  p
}

plot_shrunken_ma <- function(df, dataset_name) {
  ggplot(df, aes(safe_log10(baseMean + 1), lfc_shrunk)) +
    geom_point(alpha = 0.60, size = 1.0, color = "grey40") +
    geom_hline(yintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", linewidth = 0.6, color = plot_palette$threshold) +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "MA"), 32),
      subtitle = "Shrunken effect vs baseMean",
      x = expression(log[10]("baseMean + 1")),
      y = "Shrunken log2FC"
    ) +
    manuscript_theme()
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
    p1, p2, p3, p4, ncol = 2,
    top = textGrob(paste(pretty_dataset_label(dataset_name), "mean hist"), gp = gpar(fontface = "bold", cex = 1.0))
  )
}

plot_pca_scatter <- function(pca_fit, comparison_name, group_label) {
  pca_var <- pca_fit$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)
  cond_key <- ifelse(group_label == "control", "untrt", "trt")

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
      title = compact_title(paste(comparison_name, group_label, "PCA"), 30),
      subtitle = paste0("PC1=", pca_var_per[1], "%  PC2=", pca_var_per[2], "%"),
      x = paste0("PC1 (", pca_var_per[1], "%)"),
      y = paste0("PC2 (", pca_var_per[2], "%)")
    ) +
    manuscript_theme()
}

plot_pc1_loading_rank <- function(loading_tbl, cutoff, comparison_name, group_label, top_n_used, cutoff_quantile) {
  ggplot(loading_tbl, aes(rank, pc1_loading_abs)) +
    geom_line(linewidth = 0.4, color = "grey35") +
    geom_hline(yintercept = cutoff, color = plot_palette$threshold, linewidth = 0.8) +
    labs(
      title = compact_title(paste(comparison_name, group_label, "PC1 rank"), 30),
      subtitle = paste0("TopN=", top_n_used, "  Quantile=", signif(cutoff_quantile, 3)),
      x = "Ranked PAS",
      y = "|PC1 loading|"
    ) +
    manuscript_theme()
}

compute_mean_expression_table <- function(raw_counts, coldata) {
  sample_ids <- colnames(raw_counts)
  trt_ids    <- sample_ids[coldata$condition == "trt"]
  untrt_ids  <- sample_ids[coldata$condition == "untrt"]

  data.frame(
    feature_id = as.character(rownames(raw_counts)),
    mean_trt = rowMeans(raw_counts[, trt_ids, drop = FALSE]),
    mean_untrt = rowMeans(raw_counts[, untrt_ids, drop = FALSE]),
    mean_all = rowMeans(raw_counts),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# TWAS
# =============================================================================

clean_gene_set <- function(x) {
  unique(tolower(trimws(x[!is.na(x) & x != ""])))
}

# =============================================================================
# FULL PER-COMPARISON PIPELINE
# =============================================================================

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  cmp_dir     <- file.path(output_dir, comparison_name)
  tab_dir     <- file.path(cmp_dir, "tables")
  fig_raw_dir <- file.path(cmp_dir, "fig_raw")
  fig_le_dir  <- file.path(cmp_dir, "fig_lead")
  fig_rem_dir <- file.path(cmp_dir, "fig_rem")
  fig_cmp_dir <- file.path(cmp_dir, "fig_compare")

  for (d in c(cmp_dir, tab_dir, fig_raw_dir, fig_le_dir, fig_rem_dir, fig_cmp_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  evs <- build_eigenvector_split(count_matrix, coldata)

  save_plot(plot_pca_scatter(evs$fit_trt$pca_fit, comparison_name, "treatment"),
            file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_PCA.png")), 8.5, 6.5)
  save_plot(plot_pca_scatter(evs$fit_untrt$pca_fit, comparison_name, "control"),
            file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_PCA.png")), 8.5, 6.5)

  save_plot(
    plot_pc1_loading_rank(
      evs$fit_trt$loading_table, evs$fit_trt$cutoff, comparison_name, "treatment",
      evs$fit_trt$top_n_used, evs$fit_trt$cutoff_quantile
    ),
    file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Trt_Rank.png")),
    9, 6.5
  )
  save_plot(
    plot_pc1_loading_rank(
      evs$fit_untrt$loading_table, evs$fit_untrt$cutoff, comparison_name, "control",
      evs$fit_untrt$top_n_used, evs$fit_untrt$cutoff_quantile
    ),
    file.path(cmp_dir, paste0("Figure_", comparison_name, "_EVS_Ctrl_Rank.png")),
    9, 6.5
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
    save_csv(subset(df, standard_significant), tab_file(tab_dir, comparison_name, nm, "Std"))
    save_csv(subset(df, HBFSS_significant), tab_file(tab_dir, comparison_name, nm, "HBFSS"))
    save_csv(subset(df, effect_class == "strong_effect"), tab_file(tab_dir, comparison_name, nm, "StrongCNH"))
    save_csv(subset(df, effect_class == "weak_effect"), tab_file(tab_dir, comparison_name, nm, "WeakCNH"))

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
    save_grob(mean_panel, fig_file(fig_subdir, comparison_name, nm, "MeanHist"), 11, 8)

    p_emp  <- plot_empirical_p_histogram(df, full_dataset_name, fit$hc_p_threshold)
    p_hc   <- plot_hc_profile(df, full_dataset_name, fit$hc_p_threshold)
    p_hbdf <- plot_hbfss_distribution(df, full_dataset_name, fit$hbfss_threshold)
    p_ma   <- plot_shrunken_ma(df, full_dataset_name)
    p_vol  <- plot_volcano_single(df, full_dataset_name)
    p_disp <- plot_dispersion_single(df, full_dataset_name)

    save_plot(p_emp,  fig_file(fig_subdir, comparison_name, nm, "EmpHist"), 7.5, 6)
    if (!is.null(p_hc))   save_plot(p_hc,   fig_file(fig_subdir, comparison_name, nm, "HC"), 7.5, 6)
    save_plot(p_hbdf, fig_file(fig_subdir, comparison_name, nm, "HBFSSDist"), 7.5, 6)
    save_plot(p_ma,   fig_file(fig_subdir, comparison_name, nm, "MA"), 7.5, 6)
    save_plot(p_vol,  fig_file(fig_subdir, comparison_name, nm, "Volcano"), 7.8, 6.4)
    save_plot(p_disp, fig_file(fig_subdir, comparison_name, nm, "Disp"), 7.5, 6)

    analysis_results[[nm]] <- list(
      dds = fit$dds,
      results = df,
      summary = summary_row,
      fig_subdir = fig_subdir
    )
  }

  # Compare panels with one legend each
  volc_panel <- plot_volcano_compare_panel(
    list(
      plot_volcano_single(analysis_results$raw_dataset$results, analysis_results$raw_dataset$summary$dataset_name[1]),
      plot_volcano_single(analysis_results$leading_edge_dataset$results, analysis_results$leading_edge_dataset$summary$dataset_name[1]),
      plot_volcano_single(analysis_results$remainder_dataset$results, analysis_results$remainder_dataset$summary$dataset_name[1])
    ),
    comparison_name
  )
  save_grob(volc_panel, file.path(fig_cmp_dir, paste0("Figure_", comparison_name, "_Compare_Volcano.png")), 15.5, 6.8)

  disp_panel <- plot_dispersion_compare_panel(
    list(
      plot_dispersion_single(analysis_results$raw_dataset$results, analysis_results$raw_dataset$summary$dataset_name[1]),
      plot_dispersion_single(analysis_results$leading_edge_dataset$results, analysis_results$leading_edge_dataset$summary$dataset_name[1]),
      plot_dispersion_single(analysis_results$remainder_dataset$results, analysis_results$remainder_dataset$summary$dataset_name[1])
    ),
    comparison_name
  )
  save_grob(disp_panel, file.path(fig_cmp_dir, paste0("Figure_", comparison_name, "_Compare_Disp.png")), 15.5, 6.8)

  compare_summary <- bind_rows(lapply(names(analysis_results), function(k) {
    sm <- analysis_results[[k]]$summary
    data.frame(
      comparison_name = comparison_name,
      dataset_key = k,
      dataset_name = sm$dataset_name[1],
      n_features = sm$n_features[1],
      n_standard_significant = sm$n_standard_significant[1],
      n_HBFSS_significant = sm$n_HBFSS_significant[1],
      n_overlap = sm$n_overlap_significant[1],
      n_strong_effect = sm$n_strong_effect[1],
      n_weak_effect = sm$n_weak_effect[1],
      hc_p_threshold = sm$hc_p_threshold[1],
      hbfss_threshold = sm$hbfss_threshold[1],
      stringsAsFactors = FALSE
    )
  }))
  save_csv(compare_summary, file.path(cmp_dir, paste0("Table_", comparison_name, "_Compare_Summary.csv")))

  if (run_twas_overlap) {
    twas_genes <- clean_gene_set(TWAS_data$gene_symbol)

    get_twas_overlap <- function(result_df, dataset_nm, out_dir, mode = c("union", "standard_only", "HBFSS_only")) {
      mode <- match.arg(mode)
      sig_df <- switch(
        mode,
        union         = subset(result_df, standard_significant | HBFSS_significant),
        standard_only = subset(result_df, standard_significant),
        HBFSS_only    = subset(result_df, HBFSS_significant)
      )
      sig_df$gene_symbol_clean <- tolower(trimws(sig_df$gene_symbol))
      overlap_df <- subset(sig_df, gene_symbol_clean %in% twas_genes)
      summary_df <- data.frame(
        dataset_name = dataset_nm,
        selection_mode = mode,
        n_selected_features = nrow(sig_df),
        n_overlap_features = nrow(overlap_df),
        n_overlap_genes = length(unique(overlap_df$gene_symbol_clean)),
        stringsAsFactors = FALSE
      )
      save_csv(overlap_df, file.path(out_dir, paste0("Table_", dataset_nm, "_TWAS_", mode, ".csv")))
      save_csv(summary_df, file.path(out_dir, paste0("Table_", dataset_nm, "_TWAS_", mode, "_Summary.csv")))
      summary_df
    }

    twas_summaries <- bind_rows(lapply(names(analysis_results), function(nm) {
      full_nm <- paste(comparison_name, nm, sep = "_")
      bind_rows(
        get_twas_overlap(analysis_results[[nm]]$results, full_nm, tab_dir, "union"),
        get_twas_overlap(analysis_results[[nm]]$results, full_nm, tab_dir, "standard_only"),
        get_twas_overlap(analysis_results[[nm]]$results, full_nm, tab_dir, "HBFSS_only")
      )
    }))
    save_csv(twas_summaries, file.path(tab_dir, paste0("Table_", comparison_name, "_TWAS_Summary.csv")))
  }

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
