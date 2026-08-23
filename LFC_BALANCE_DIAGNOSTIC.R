#!/usr/bin/env Rscript

# =============================================================================
# LFC BALANCE DIAGNOSTIC: does EVS partitioning improve apeglm shrinkage
# calibration by making the LFC distribution more balanced?
# =============================================================================
#
# This script is standalone. It does not modify SEQUENCE.R, EMPERICALCUTOFF.R,
# REGIME_DIAGNOSTIC.r, or REMNB1_LEADNB2.R. It reads the same raw count file
# and reuses the same PC1-ranking / Leading-Edge / Remainder conventions as
# REGIME_DIAGNOSTIC.r (log1p-CPM, arm-specific PCA, top-k union), and the
# same per-comparison empirical k* lookup from EMPERICALCUTOFF.R's
# Table_Timepoints.csv, so gene-set membership matches exactly.
#
# CLAIM BEING TESTED
# apeglm fits its shrinkage prior from the empirical distribution of MLE
# log2 fold changes across the genes being tested together. If that
# distribution is skewed or otherwise imbalanced, the fitted prior can be
# miscalibrated (over-shrinking one side, under-shrinking the other). The
# claim under test is that partitioning genes into Leading Edge and
# Remainder before running DESeq2+apeglm -- rather than running apeglm once
# on the unpartitioned Original dataset -- produces a more balanced LFC
# distribution within each partition, and therefore better-calibrated
# shrinkage.
#
# WHAT THIS SCRIPT DOES
# For each comparison, at both the paper-reference k (5000, fixed) and that
# comparison's own empirical k* (from Table_Timepoints.csv):
#   1. Run DESeq2 + apeglm on the Original (unpartitioned, all features)
#      dataset for that comparison.
#   2. Run DESeq2 + apeglm independently on the Leading Edge subset only.
#   3. Run DESeq2 + apeglm independently on the Remainder subset only.
# For each of the three resulting apeglm-shrunk LFC distributions, this
# reports:
#   - skewness (Fisher-Pearson standardized third moment)
#   - an asymptotic z-test of skewness against 0 (skewness / SE, with
#     SE = sqrt(6n(n-1) / ((n-2)(n+1)(n+3))), a standard large-sample
#     approximation), with a two-sided p-value
#   - the proportion of tested genes with LFC > 0 vs LFC < 0, and a
#     two-sided binomial test of that proportion against 0.5
#   - mean and median LFC
# A genuinely more "balanced" partition should show skewness and the
# up/down proportion closer to zero/0.5 than the Original run, with the
# corresponding tests failing to reject balance where Original's test
# rejects it.
#
# IMPORTANT CAVEAT: this is a real, non-circular comparison in one specific
# sense -- Leading Edge and Remainder are each refit independently with
# their own DESeq2 size factors and dispersions (as SEQUENCE.R already does
# for its Lead/Rem views), so nothing about the DESeq2/apeglm fit itself
# was told about the other partition. It is not a claim that the partition
# criterion (PC1 loading) is independent of the LFC distribution -- that is
# the separate question already addressed by the permutation diagnostic in
# REGIME_DIAGNOSTIC.r.
#
# OUTPUTS (written to a directory separate from every other script's output):
#   Table_LFC_Balance.csv
#   Figure_LFC_Balance_<comparison>.png  (one density plot per comparison)
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Settings (kept consistent with REGIME_DIAGNOSTIC.r)
# -----------------------------------------------------------------------------

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT   <- "/root/REAPER98632/exports/lfc_balance_diagnostic"
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

GROUP_PATTERNS <- c(
  RT0  = "^R0_",
  ZT6  = "^ZT6_",
  RT2  = "^R2_",
  ZT8  = "^ZT8_",
  RT4  = "^R4_",
  ZT10 = "^ZT10_",
  RT8  = "^R8_",
  ZT14 = "^ZT14_"
)

COMPARISONS <- list(
  RT0_ZT6  = c(control = "RT0", treatment = "ZT6"),
  RT2_ZT8  = c(control = "RT2", treatment = "ZT8"),
  RT4_ZT10 = c(control = "RT4", treatment = "ZT10"),
  RT8_ZT14 = c(control = "RT8", treatment = "ZT14")
)

CUTOFF_TIMEPOINTS_FILE <- "/root/REAPER98632/exports/pc1_nb_manuscript_final/Table_Timepoints.csv"
PAPER_REFERENCE_K <- 5000L

load_pairwise_k_star <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Could not find ", path, ". Run EMPERICALCUTOFF_MANUSCRIPT_METHODS_ONLY.r ",
      "first so its per-comparison pairwise_weighted_k values exist to read."
    )
  }
  tp <- read.csv(path, stringsAsFactors = FALSE)
  required_cols <- c("comparison", "pairwise_weighted_k")
  missing_cols <- setdiff(required_cols, names(tp))
  if (length(missing_cols) > 0L) {
    stop("Table_Timepoints.csv is missing expected column(s): ", paste(missing_cols, collapse = ", "))
  }
  setNames(as.integer(round(tp$pairwise_weighted_k)), tp$comparison)
}

# -----------------------------------------------------------------------------
# Data loading (mirrors REGIME_DIAGNOSTIC.r conventions)
# -----------------------------------------------------------------------------

read_count_matrix <- function(path, group_patterns) {
  if (!file.exists(path)) stop("Count file does not exist: ", path)
  raw_df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(raw_df) < 1L || ncol(raw_df) < 2L) stop("Count file is empty or malformed.")

  sample_idx <- sort(unique(unlist(
    lapply(group_patterns, function(pattern) grep(pattern, colnames(raw_df)))
  )))
  if (length(sample_idx) == 0L) stop("No sample columns matched GROUP_PATTERNS.")

  feature_ids <- as.character(raw_df[[1]])
  count_mat <- as.matrix(raw_df[, sample_idx, drop = FALSE])
  storage.mode(count_mat) <- "numeric"
  rownames(count_mat) <- feature_ids

  n_na <- sum(is.na(count_mat))
  if (n_na > 0L) {
    message(sprintf("Count matrix contained %d NA cell(s) after loading; treating as 0.", n_na))
    count_mat[is.na(count_mat)] <- 0
  }

  all_zero <- rowSums(count_mat, na.rm = TRUE) == 0
  if (any(all_zero)) {
    message(sprintf("Dropping %d feature(s) with zero counts across every sample.", sum(all_zero)))
    count_mat <- count_mat[!all_zero, , drop = FALSE]
  }

  count_mat
}

normalize_cpm_log1p <- function(count_mat_arm) {
  lib_size <- colSums(count_mat_arm, na.rm = TRUE)
  lib_size[!is.finite(lib_size) | lib_size <= 0] <- 1
  cpm <- sweep(count_mat_arm, 2L, lib_size / 1e6, "/")
  result <- log1p(cpm)
  result[!is.finite(result)] <- 0
  result
}

compute_pc1_rank <- function(rank_matrix_arm) {
  input_mat <- t(rank_matrix_arm)
  if (any(!is.finite(input_mat))) {
    stop("Non-finite values reached compute_pc1_rank after normalization.")
  }
  pca <- stats::prcomp(input_mat, center = TRUE, scale. = FALSE, rank. = 1)
  loading <- pca$rotation[, 1L]
  loading[!is.finite(loading)] <- 0
  abs(loading)
}

top_k_features <- function(abs_loading, k) {
  k <- min(k, length(abs_loading))
  names(sort(abs_loading, decreasing = TRUE))[seq_len(k)]
}

leading_edge_set <- function(count_mat, sample_ids_arm1, sample_ids_arm2, k) {
  mat1 <- normalize_cpm_log1p(count_mat[, sample_ids_arm1, drop = FALSE])
  mat2 <- normalize_cpm_log1p(count_mat[, sample_ids_arm2, drop = FALSE])
  loading1 <- compute_pc1_rank(mat1)
  loading2 <- compute_pc1_rank(mat2)
  top1 <- top_k_features(loading1, k)
  top2 <- top_k_features(loading2, k)
  union(top1, top2)
}

# -----------------------------------------------------------------------------
# DESeq2 + apeglm on an arbitrary gene subset
# -----------------------------------------------------------------------------

run_deseq2_apeglm_lfc <- function(count_mat, control_ids, treatment_ids, feature_subset = NULL) {
  if (!is.null(feature_subset)) {
    keep <- rownames(count_mat) %in% feature_subset
    count_mat <- count_mat[keep, , drop = FALSE]
  }
  if (nrow(count_mat) < 10L) {
    return(NULL)
  }

  all_ids <- c(control_ids, treatment_ids)
  sub_mat <- round(count_mat[, all_ids, drop = FALSE])
  condition <- factor(
    c(rep("untrt", length(control_ids)), rep("trt", length(treatment_ids))),
    levels = c("untrt", "trt")
  )
  col_data <- data.frame(condition = condition, row.names = all_ids)

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = sub_mat, colData = col_data, design = ~condition)
  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
  dds <- tryCatch(
    suppressMessages(DESeq2::DESeq(dds, betaPrior = FALSE)),
    error = function(e) NULL
  )
  if (is.null(dds)) return(NULL)

  coef_name <- "condition_trt_vs_untrt"
  if (!coef_name %in% DESeq2::resultsNames(dds)) {
    coef_name <- DESeq2::resultsNames(dds)[length(DESeq2::resultsNames(dds))]
  }

  shrunk <- tryCatch(
    suppressMessages(DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")),
    error = function(e) NULL
  )
  if (is.null(shrunk)) return(NULL)

  lfc <- as.numeric(shrunk$log2FoldChange)
  lfc <- lfc[is.finite(lfc)]
  lfc
}

# -----------------------------------------------------------------------------
# Balance statistics
# -----------------------------------------------------------------------------

skewness_stat <- function(x) {
  n <- length(x)
  if (n < 3L) return(NA_real_)
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  (sum((x - m)^3) / n) / (s^3)
}

skewness_z_test <- function(x) {
  n <- length(x)
  skew <- skewness_stat(x)
  if (n < 8L || is.na(skew)) {
    return(list(skewness = skew, se = NA_real_, z = NA_real_, p_value = NA_real_))
  }
  se <- sqrt((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3)))
  z <- skew / se
  p <- 2 * stats::pnorm(-abs(z))
  list(skewness = skew, se = se, z = z, p_value = p)
}

updown_balance_test <- function(x) {
  n_up <- sum(x > 0)
  n_down <- sum(x < 0)
  n_total <- n_up + n_down
  if (n_total == 0L) {
    return(list(prop_up = NA_real_, p_value = NA_real_))
  }
  bt <- stats::binom.test(n_up, n_total, p = 0.5)
  list(prop_up = n_up / n_total, p_value = bt$p.value)
}

summarize_lfc_balance <- function(lfc, comparison_name, k, k_type, view_name) {
  if (is.null(lfc) || length(lfc) < 10L) {
    return(data.frame(
      comparison = comparison_name, k = k, k_type = k_type, view = view_name,
      n_genes = length(lfc), mean_lfc = NA_real_, median_lfc = NA_real_,
      skewness = NA_real_, skewness_p = NA_real_,
      prop_up = NA_real_, updown_p = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  skew <- skewness_z_test(lfc)
  ud <- updown_balance_test(lfc)
  data.frame(
    comparison = comparison_name, k = k, k_type = k_type, view = view_name,
    n_genes = length(lfc),
    mean_lfc = mean(lfc), median_lfc = stats::median(lfc),
    skewness = skew$skewness, skewness_p = skew$p_value,
    prop_up = ud$prop_up, updown_p = ud$p_value,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------

message("Reading count matrix...")
count_mat <- read_count_matrix(COUNT_FILE, GROUP_PATTERNS)

message("Reading per-comparison empirical cutoffs from ", CUTOFF_TIMEPOINTS_FILE, " ...")
pairwise_k_star <- load_pairwise_k_star(CUTOFF_TIMEPOINTS_FILE)

result_rows <- list()
lfc_long_rows <- list()

for (comparison_name in names(COMPARISONS)) {
  pair <- COMPARISONS[[comparison_name]]
  control_ids <- grep(GROUP_PATTERNS[[pair[["control"]]]], colnames(count_mat), value = TRUE)
  treatment_ids <- grep(GROUP_PATTERNS[[pair[["treatment"]]]], colnames(count_mat), value = TRUE)

  if (!comparison_name %in% names(pairwise_k_star) || is.na(pairwise_k_star[[comparison_name]])) {
    stop("No pairwise_weighted_k found for comparison '", comparison_name, "' in ", CUTOFF_TIMEPOINTS_FILE)
  }
  k_values <- c(PAPER_REFERENCE_K, pairwise_k_star[[comparison_name]])

  message(sprintf("[%s] control n=%d, treatment n=%d", comparison_name, length(control_ids), length(treatment_ids)))

  # Original (unpartitioned) is the same for both k values, so compute once.
  message("  Running Original (unpartitioned)...")
  lfc_original <- run_deseq2_apeglm_lfc(count_mat, control_ids, treatment_ids, feature_subset = NULL)

  for (k in k_values) {
    k_type <- if (k == PAPER_REFERENCE_K) "paper_reference" else "empirical_k_star"

    le_ids <- leading_edge_set(count_mat, control_ids, treatment_ids, k)
    rem_ids <- setdiff(rownames(count_mat), le_ids)

    message(sprintf("  k=%d (%s): Leading Edge n=%d, Remainder n=%d", k, k_type, length(le_ids), length(rem_ids)))

    lfc_le <- run_deseq2_apeglm_lfc(count_mat, control_ids, treatment_ids, feature_subset = le_ids)
    lfc_rem <- run_deseq2_apeglm_lfc(count_mat, control_ids, treatment_ids, feature_subset = rem_ids)

    result_rows[[length(result_rows) + 1L]] <- summarize_lfc_balance(lfc_original, comparison_name, k, k_type, "Original")
    result_rows[[length(result_rows) + 1L]] <- summarize_lfc_balance(lfc_le, comparison_name, k, k_type, "Leading_Edge")
    result_rows[[length(result_rows) + 1L]] <- summarize_lfc_balance(lfc_rem, comparison_name, k, k_type, "Remainder")

    if (identical(k_type, "empirical_k_star")) {
      if (!is.null(lfc_original)) {
        lfc_long_rows[[length(lfc_long_rows) + 1L]] <- data.frame(
          comparison = comparison_name, view = "Original", lfc = lfc_original
        )
      }
      if (!is.null(lfc_le)) {
        lfc_long_rows[[length(lfc_long_rows) + 1L]] <- data.frame(
          comparison = comparison_name, view = "Leading_Edge", lfc = lfc_le
        )
      }
      if (!is.null(lfc_rem)) {
        lfc_long_rows[[length(lfc_long_rows) + 1L]] <- data.frame(
          comparison = comparison_name, view = "Remainder", lfc = lfc_rem
        )
      }
    }
  }
}

result_table <- dplyr::bind_rows(result_rows)
write.csv(result_table, file.path(OUT_ROOT, "Table_LFC_Balance.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Figure: LFC density by view, one panel per comparison (empirical k* only)
# -----------------------------------------------------------------------------

if (length(lfc_long_rows) > 0L) {
  lfc_long <- dplyr::bind_rows(lfc_long_rows)
  lfc_long$view <- factor(lfc_long$view, levels = c("Original", "Leading_Edge", "Remainder"))

  p <- ggplot(lfc_long, aes(x = lfc, color = view, fill = view)) +
    geom_density(alpha = 0.15, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    facet_wrap(~comparison, scales = "free", ncol = 2) +
    labs(
      title = "apeglm-shrunk LFC distribution: Original vs Leading Edge vs Remainder",
      subtitle = "Each comparison's own empirical k*. A more symmetric, less skewed curve indicates better-balanced shrinkage input.",
      x = "apeglm-shrunk log2 fold change",
      y = "Density",
      color = NULL, fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(OUT_ROOT, "Figure_LFC_Balance_AllComparisons.png"), p, width = 11, height = 8, dpi = 300)
}

message("Diagnostic complete. Outputs written to: ", OUT_ROOT)
message("  Table_LFC_Balance.csv")
message("  Figure_LFC_Balance_AllComparisons.png")

# -----------------------------------------------------------------------------
# Zip archives: everything in one download for figures and for tables.
# Wrapped in tryCatch so that if the zip utility is unavailable, the actual
# diagnostic results above are still kept -- only the packaging step is lost.
# -----------------------------------------------------------------------------

message("Creating zip archives...")

zip_result <- tryCatch(
  {
    figure_files <- list.files(OUT_ROOT, pattern = "\\.png$", full.names = TRUE)
    table_files <- list.files(OUT_ROOT, pattern = "\\.csv$", full.names = TRUE)

    if (length(figure_files) > 0L) {
      figures_zip_path <- file.path(OUT_ROOT, "LFC_Balance_Figures.zip")
      if (file.exists(figures_zip_path)) file.remove(figures_zip_path)
      utils::zip(figures_zip_path, files = figure_files, flags = "-j")
      message("  LFC_Balance_Figures.zip (", length(figure_files), " files)")
    } else {
      message("  No .png files found; skipping LFC_Balance_Figures.zip")
    }

    if (length(table_files) > 0L) {
      tables_zip_path <- file.path(OUT_ROOT, "LFC_Balance_Tables.zip")
      if (file.exists(tables_zip_path)) file.remove(tables_zip_path)
      utils::zip(tables_zip_path, files = table_files, flags = "-j")
      message("  LFC_Balance_Tables.zip (", length(table_files), " files)")
    } else {
      message("  No .csv files found; skipping LFC_Balance_Tables.zip")
    }

    TRUE
  },
  error = function(e) {
    message(
      "Zip archive creation failed (individual files above are still ",
      "intact and usable): ", conditionMessage(e)
    )
    FALSE
  }
)
