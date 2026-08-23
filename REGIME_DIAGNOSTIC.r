#!/usr/bin/env Rscript

# =============================================================================
# REGIME DIAGNOSTIC: EVS LEADING-EDGE / REMAINDER VALIDATION
# =============================================================================
#
# This script is standalone. It does not modify or overwrite SEQUENCE.R or
# EMPERICALCUTOFF.R. It reads the same raw count file and reuses the same
# PC1-ranking convention (log1p-CPM, arm-specific PCA, ranked by ascending
# |PC1 loading|) so results are directly comparable to the existing pipeline.
#
# DEPENDENCY: this script requires EMPERICALCUTOFF.R to have already been
# run at least once, because it reads each comparison's own empirically-
# derived cutoff (pairwise_weighted_k) directly from that script's
# Table_Timepoints.csv output. Run EMPERICALCUTOFF.R first if that file
# does not yet exist.
#
# It answers two questions, run under TWO cutoffs per comparison:
#   - the paper-reference k = 5000, the same fixed value for every comparison
#   - that comparison's own empirical k* (pairwise_weighted_k), which
#     differs by comparison (e.g. RT0_ZT6 and RT4_ZT10 do not share a value)
#
# DIAGNOSTIC 1: Does Remainder behave like an NB1 (linear mean-variance)
# regime and Leading Edge like an NB2 (quadratic mean-variance) regime?
#   Within each set, per-feature mean and variance come from DESeq2's own
#   dispersion-shrinkage pipeline (estimateDispersions), not raw unshrunk
#   per-feature sample variance: Var = mean + dispersion_final * mean^2,
#   using DESeq2's final (empirical-Bayes MAP shrunk) dispersion estimate.
#   This avoids letting a handful of noisy, low-replicate-count features
#   dominate the curve fit, which raw sample variance is prone to with
#   small n per group. Excess-over-Poisson variance (variance minus mean,
#   floored at 0) is regressed against mean (NB1, linear) and against
#   mean^2 (NB2, quadratic) separately for Remainder and for Leading Edge.
#   Both R-squared and AIC are reported for each model in each set; AIC
#   (delta_aic / preferred columns) is the primary comparison, since it is
#   a more rigorous way to judge two non-nested models than raw R-squared
#   on a through-origin fit. This directly tests whether Remainder fits
#   NB1 better and Leading Edge fits NB2 better.
#
# DIAGNOSTIC 2: is the observed NB2-NB1 gap between Leading Edge and
# Remainder bigger than a random partition of the same size would give?
#   Diagnostic A tests whether gene-SET MEMBERSHIP survives relabeling --
#   a question about whether Leading Edge is a gene-intrinsic property.
#   Diagnostic C instead directly tests the specific quantity motivating
#   the two-regime claim: gap = log1p(variance - mean) - log1p(mean), i.e.
#   NB2 minus NB1, using the exact same formula as REMNB1_LEADNB2.R. Using
#   the real (unpermuted) per-gene mean/variance throughout, the real
#   median gap of Leading Edge minus the real median gap of Remainder is
#   compared against a null built by repeatedly drawing a random gene
#   subset of the same size as Leading Edge (independent of PC1 rank) and
#   computing the same statistic. A one-sided permutation p-value reports
#   how often a random partition matches or exceeds the real, PC1-defined
#   partition's gap; a z-like effect size (real minus null mean, divided
#   by null sd) gives a standardized magnitude.
#
# OUTPUTS (written to a directory separate from the existing pipeline):
#   Table_NB_Regime_Fit.csv
#   Table_Gap_Permutation.csv
#   Figure_NB_Regime_Fit.png
#   Figure_NB_Regime_AIC.png
#   Figure_Gap_Permutation.png
#   RegimeDiagnostic_Figures.zip
#   RegimeDiagnostic_Tables.zip
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

options(stringsAsFactors = FALSE)
set.seed(20260823L)

# -----------------------------------------------------------------------------
# Settings (kept consistent with EMPERICALCUTOFF_MANUSCRIPT_METHODS_ONLY.R)
# -----------------------------------------------------------------------------

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <- "/root/REAPER98632/exports/regime_diagnostic"
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
    stop(
      "Table_Timepoints.csv is missing expected column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  setNames(as.integer(round(tp$pairwise_weighted_k)), tp$comparison)
}
N_PERMUTATIONS <- 500L

# -----------------------------------------------------------------------------
# Data loading (mirrors EMPERICALCUTOFF_MANUSCRIPT_METHODS_ONLY.R conventions)
# -----------------------------------------------------------------------------

read_count_matrix <- function(path, group_patterns) {
  if (!file.exists(path)) {
    stop("Count file does not exist: ", path)
  }

  raw_df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  if (nrow(raw_df) < 1L || ncol(raw_df) < 2L) {
    stop("Count file is empty or malformed.")
  }

  sample_idx <- sort(unique(unlist(
    lapply(group_patterns, function(pattern) grep(pattern, colnames(raw_df)))
  )))

  if (length(sample_idx) == 0L) {
    stop("No sample columns matched GROUP_PATTERNS.")
  }

  feature_ids <- as.character(raw_df[[1]])
  count_mat <- as.matrix(raw_df[, sample_idx, drop = FALSE])
  storage.mode(count_mat) <- "numeric"
  rownames(count_mat) <- feature_ids

  # The raw CSV can contain blank cells, non-numeric placeholders (e.g. "NA",
  # "-"), or genuinely missing values, all of which become NA on coercion to
  # numeric above. PCA/svd cannot run with any NA/Inf present anywhere in the
  # matrix, so these are resolved explicitly here rather than failing deep
  # inside prcomp with an opaque "infinite or missing values" error.
  n_na <- sum(is.na(count_mat))
  n_inf <- sum(is.infinite(count_mat))
  if (n_na > 0L || n_inf > 0L) {
    message(sprintf(
      "Count matrix contained %d NA and %d non-finite cell(s) after loading; treating as 0.",
      n_na, n_inf
    ))
    count_mat[is.na(count_mat) | is.infinite(count_mat)] <- 0
  }

  # A feature with zero counts across every sample carries no PC1 information
  # and can still destabilize downstream variance calculations; drop it here
  # rather than silently propagating zeros/NaNs through PCA and DESeq2.
  all_zero <- rowSums(count_mat, na.rm = TRUE) == 0
  if (any(all_zero)) {
    message(sprintf(
      "Dropping %d feature(s) with zero counts across every sample.",
      sum(all_zero)
    ))
    count_mat <- count_mat[!all_zero, , drop = FALSE]
  }

  count_mat
}

normalize_cpm_log1p <- function(count_mat_arm) {
  lib_size <- colSums(count_mat_arm, na.rm = TRUE)
  lib_size[!is.finite(lib_size) | lib_size <= 0] <- 1
  cpm <- sweep(count_mat_arm, 2L, lib_size / 1e6, "/")
  result <- log1p(cpm)
  # Defensive final guard: no NA/Inf should reach prcomp regardless of cause.
  result[!is.finite(result)] <- 0
  result
}

compute_pc1_rank <- function(rank_matrix_arm) {
  input_mat <- t(rank_matrix_arm)
  if (any(!is.finite(input_mat))) {
    stop(
      "Non-finite values reached compute_pc1_rank after normalization; ",
      "this indicates malformed input data upstream of PCA."
    )
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
# DIAGNOSTIC 1 helpers: NB1 (linear) vs NB2 (quadratic) mean-variance fit
# -----------------------------------------------------------------------------

deseq2_shrunk_mean_variance <- function(count_mat, group_labels) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required for the NB regime diagnostic.")
  }
  col_data <- data.frame(group = factor(group_labels), row.names = colnames(count_mat))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = col_data,
    design = ~group
  )
  dds <- tryCatch(
    DESeq2::estimateSizeFactors(dds),
    error = function(e) DESeq2::estimateSizeFactors(dds, type = "poscounts")
  )

  # estimateDispersions fits gene-wise dispersion, the mean-dispersion trend,
  # and the final MAP (empirical-Bayes shrunk) dispersion per gene -- this is
  # the same shrinkage DESeq2 normally applies before testing. Using the
  # final shrunk dispersion here (instead of raw unshrunk per-feature sample
  # variance) avoids letting a handful of noisy, low-replicate-count
  # features dominate the NB1-vs-NB2 curve fit.
  dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)

  base_mean <- rowMeans(DESeq2::counts(dds, normalized = TRUE))
  disp_final <- DESeq2::dispersions(dds)

  # DESeq2 flags a gene as a dispersion outlier when its gene-wise estimate
  # sits too far above the fitted trend (default ~2 residual SD); for those
  # genes, dispersion_final above is NOT shrunk toward the trend -- it
  # reverts to the raw gene-wise MLE. Because the trend is fit on the whole
  # dataset (dominated by the much larger Remainder population), Leading
  # Edge genes -- if they systematically sit higher on mean/variance -- are
  # disproportionately likely to be flagged this way, which would make
  # "dispersion_final" for Leading Edge closer to an unshrunk estimate than
  # the label implies. This flag makes that checkable rather than assumed.
  disp_outlier <- tryCatch(
    S4Vectors::mcols(dds)$dispOutlier,
    error = function(e) rep(NA, nrow(dds))
  )
  if (is.null(disp_outlier)) disp_outlier <- rep(NA, nrow(dds))

  # NB2 parameterization: Var = mu + dispersion * mu^2. This is the variance
  # implied by DESeq2's own final shrunk dispersion estimate, analogous to
  # the "final dispersion estimates" used in the Figure 22 mean-variance fit.
  shrunk_variance <- base_mean + disp_final * base_mean^2

  keep <- is.finite(base_mean) & is.finite(shrunk_variance) & !is.na(disp_final)

  data.frame(
    feature_id = rownames(count_mat)[keep],
    mean = base_mean[keep],
    variance = shrunk_variance[keep],
    dispersion_final = disp_final[keep],
    dispersion_is_outlier = disp_outlier[keep],
    stringsAsFactors = FALSE
  )
}

fit_nb_regimes <- function(mean_var_df, feature_subset) {
  sub <- mean_var_df[mean_var_df$feature_id %in% feature_subset, , drop = FALSE]

  outlier_col <- sub$dispersion_is_outlier
  frac_outlier <- if (length(outlier_col) > 0L && !all(is.na(outlier_col))) {
    mean(outlier_col, na.rm = TRUE)
  } else {
    NA_real_
  }

  sub <- sub[is.finite(sub$mean) & is.finite(sub$variance) & sub$mean > 0, , drop = FALSE]
  sub$excess <- pmax(sub$variance - sub$mean, 0)

  if (nrow(sub) < 10L) {
    return(data.frame(
      model = c("NB1_linear", "NB2_quadratic"),
      r_squared = NA_real_,
      aic = NA_real_,
      delta_aic = NA_real_,
      preferred = NA,
      n_features = nrow(sub),
      frac_dispersion_outlier = frac_outlier
    ))
  }

  nb1_fit <- lm(excess ~ 0 + mean, data = sub)
  nb2_fit <- lm(excess ~ 0 + I(mean^2), data = sub)

  r2 <- function(fit, y) {
    pred <- predict(fit)
    ss_res <- sum((y - pred)^2)
    ss_tot <- sum((y - mean(y))^2)
    if (ss_tot <= 0) return(NA_real_)
    1 - ss_res / ss_tot
  }

  aic1 <- stats::AIC(nb1_fit)
  aic2 <- stats::AIC(nb2_fit)
  best_aic <- min(aic1, aic2)

  data.frame(
    model = c("NB1_linear", "NB2_quadratic"),
    r_squared = c(r2(nb1_fit, sub$excess), r2(nb2_fit, sub$excess)),
    aic = c(aic1, aic2),
    delta_aic = c(aic1 - best_aic, aic2 - best_aic),
    preferred = c(aic1 < aic2, aic2 < aic1),
    n_features = nrow(sub),
    frac_dispersion_outlier = frac_outlier
  )
}

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------

message("Reading count matrix...")
count_mat <- read_count_matrix(COUNT_FILE, GROUP_PATTERNS)

message("Reading per-comparison empirical cutoffs from ", CUTOFF_TIMEPOINTS_FILE, " ...")
pairwise_k_star <- load_pairwise_k_star(CUTOFF_TIMEPOINTS_FILE)

nb_rows <- list()
gap_permutation_rows <- list()
gap_null_long_rows <- list()

for (comparison_name in names(COMPARISONS)) {
  pair <- COMPARISONS[[comparison_name]]
  control_label <- pair[["control"]]
  treatment_label <- pair[["treatment"]]

  control_ids <- grep(GROUP_PATTERNS[[control_label]], colnames(count_mat), value = TRUE)
  treatment_ids <- grep(GROUP_PATTERNS[[treatment_label]], colnames(count_mat), value = TRUE)
  all_ids <- c(control_ids, treatment_ids)
  n_control <- length(control_ids)
  n_treatment <- length(treatment_ids)

  group_labels <- c(rep("control", n_control), rep("treatment", n_treatment))
  names(group_labels) <- all_ids

  if (!comparison_name %in% names(pairwise_k_star) || is.na(pairwise_k_star[[comparison_name]])) {
    stop(
      "No pairwise_weighted_k found for comparison '", comparison_name,
      "' in ", CUTOFF_TIMEPOINTS_FILE
    )
  }
  comparison_k_values <- c(PAPER_REFERENCE_K, pairwise_k_star[[comparison_name]])
  names(comparison_k_values) <- NULL

  message(sprintf(
    "[%s] control n=%d, treatment n=%d, k values = %s (paper reference + comparison-specific empirical k*)",
    comparison_name, n_control, n_treatment, paste(comparison_k_values, collapse = ", ")
  ))

  comparison_mat <- count_mat[, all_ids, drop = FALSE]

  # Computed once per comparison (does not depend on k) and shared by
  # Diagnostic B and the new Diagnostic C below. Previously this was
  # recomputed inside the k loop even though it doesn't depend on k.
  message("  Fitting DESeq2 shrunk dispersion for this comparison...")
  mean_var_df <- deseq2_shrunk_mean_variance(comparison_mat, group_labels)

  # Same NB2 / NB1 / gap formulas as REMNB1_LEADNB2.R's compute_ranked_feature_metrics,
  # computed here genome-wide (every gene in the comparison) rather than
  # only within a local LEFT/RIGHT window at the boundary.
  mean_var_df$nb2 <- log1p(pmax(mean_var_df$variance - mean_var_df$mean, 0))
  mean_var_df$nb1 <- log1p(mean_var_df$mean)
  mean_var_df$gap <- mean_var_df$nb2 - mean_var_df$nb1

  for (k in comparison_k_values) {

    k_type <- if (k == PAPER_REFERENCE_K) "paper_reference" else "empirical_k_star"

    # --- Real (unpermuted) Leading Edge / Remainder ---------------------
    real_le <- leading_edge_set(comparison_mat, control_ids, treatment_ids, k)
    real_remainder <- setdiff(rownames(comparison_mat), real_le)

    # --- DIAGNOSTIC 1: NB1 vs NB2 regime fit -----------------------------
    # Uses DESeq2's final shrunk dispersion (empirical-Bayes MAP estimate),
    # not raw unshrunk per-feature sample variance, so a handful of noisy
    # low-replicate features cannot dominate the curve fit. mean_var_df was
    # computed once per comparison above, since it does not depend on k.
    fit_le <- fit_nb_regimes(mean_var_df, real_le)
    fit_remainder <- fit_nb_regimes(mean_var_df, real_remainder)

    fit_le$feature_set <- "Leading_Edge"
    fit_remainder$feature_set <- "Remainder"
    fit_le$comparison <- comparison_name
    fit_remainder$comparison <- comparison_name
    fit_le$k <- k
    fit_remainder$k <- k
    fit_le$k_type <- k_type
    fit_remainder$k_type <- k_type

    nb_rows[[length(nb_rows) + 1L]] <- fit_le
    nb_rows[[length(nb_rows) + 1L]] <- fit_remainder

    # --- DIAGNOSTIC 2: direct permutation test of the NB2-NB1 gap --------
    # Diagnostic A tests whether gene-SET MEMBERSHIP survives relabeling.
    # This is a different, more direct question: is the actual NB2-NB1 gap
    # between Leading Edge and Remainder bigger than what an arbitrary
    # random partition of the same sizes would produce by chance, using
    # the exact same gap statistic already used throughout REMNB1_LEADNB2.R
    # (gap = log1p(variance - mean) - log1p(mean), i.e. NB2 minus NB1).
    #
    # gene mu/variance come from the real (unpermuted) data throughout --
    # only which genes get labeled "Leading Edge" vs "Remainder" is
    # permuted, drawing a random gene subset of the same size as the real
    # Leading Edge set, uniformly at random from every gene in the
    # comparison, independent of PC1 rank.
    real_gap_le <- mean_var_df$gap[mean_var_df$feature_id %in% real_le]
    real_gap_rem <- mean_var_df$gap[mean_var_df$feature_id %in% real_remainder]
    diff_gap_real <- stats::median(real_gap_le, na.rm = TRUE) - stats::median(real_gap_rem, na.rm = TRUE)

    all_feature_ids <- mean_var_df$feature_id
    n_le <- length(real_le)
    null_diffs <- numeric(N_PERMUTATIONS)

    for (p in seq_len(N_PERMUTATIONS)) {
      random_le_ids <- sample(all_feature_ids, n_le)
      random_rem_ids <- setdiff(all_feature_ids, random_le_ids)
      gap_random_le <- mean_var_df$gap[mean_var_df$feature_id %in% random_le_ids]
      gap_random_rem <- mean_var_df$gap[mean_var_df$feature_id %in% random_rem_ids]
      null_diffs[p] <- stats::median(gap_random_le, na.rm = TRUE) - stats::median(gap_random_rem, na.rm = TRUE)
    }

    # One-sided: is the real gap bigger than a random partition would give?
    # The +1 in numerator and denominator is the standard permutation-test
    # correction (the observed statistic counts as one of its own draws),
    # so the p-value is never reported as exactly 0.
    perm_p_gap <- (1 + sum(null_diffs >= diff_gap_real)) / (N_PERMUTATIONS + 1)

    gap_permutation_rows[[length(gap_permutation_rows) + 1L]] <- data.frame(
      comparison = comparison_name,
      k = k,
      k_type = k_type,
      diff_gap_real = diff_gap_real,
      null_mean = mean(null_diffs, na.rm = TRUE),
      null_sd = stats::sd(null_diffs, na.rm = TRUE),
      effect_size_z = (diff_gap_real - mean(null_diffs, na.rm = TRUE)) / stats::sd(null_diffs, na.rm = TRUE),
      perm_p_one_sided = perm_p_gap,
      n_permutations = N_PERMUTATIONS,
      stringsAsFactors = FALSE
    )

    gap_null_long_rows[[length(gap_null_long_rows) + 1L]] <- data.frame(
      comparison = comparison_name,
      k = k,
      k_type = k_type,
      null_diff = null_diffs,
      stringsAsFactors = FALSE
    )
  }
}

nb_table <- dplyr::bind_rows(nb_rows)
gap_permutation_table <- dplyr::bind_rows(gap_permutation_rows)
gap_null_long_table <- dplyr::bind_rows(gap_null_long_rows)

write.csv(
  nb_table,
  file.path(OUT_ROOT, "Table_NB_Regime_Fit.csv"),
  row.names = FALSE
)
write.csv(
  gap_permutation_table,
  file.path(OUT_ROOT, "Table_Gap_Permutation.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

nb_plot_df <- nb_table %>%
  dplyr::mutate(
    k_type_label = dplyr::recode(
      k_type,
      paper_reference = "Paper reference (k=5000, same for all)",
      empirical_k_star = "Empirical k* (comparison-specific)"
    ),
    comparison_k_label = paste0(comparison, " (k=", k, ")")
  )

p_nb_r2 <- ggplot(nb_plot_df, aes(x = model, y = r_squared, fill = feature_set)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(k_type_label ~ comparison_k_label) +
  labs(
    title = "Diagnostic B: NB1 (linear) vs NB2 (quadratic) mean-variance fit (R-squared)",
    subtitle = "Remainder fitting NB1 better and Leading Edge fitting NB2 better would support the two-regime claim. Each comparison's empirical panel uses its own k* (see facet labels).",
    x = NULL,
    y = expression(R^2),
    fill = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(
  file.path(OUT_ROOT, "Figure_NB_Regime_Fit.png"),
  p_nb_r2, width = 13, height = 7, dpi = 300
)

# AIC is the more rigorous comparison: delta_aic is 0 for the preferred model
# and positive for the other, so a taller bar means a more decisively rejected
# alternative. This is the primary figure for judging the two-regime claim;
# the R-squared figure above is kept for interpretability alongside it.
p_nb_aic <- ggplot(nb_plot_df, aes(x = model, y = delta_aic, fill = feature_set)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(k_type_label ~ comparison_k_label) +
  labs(
    title = "Diagnostic B: NB1 vs NB2 model comparison by AIC",
    subtitle = "Delta AIC = 0 marks the preferred model for that set; taller bars mean the alternative model is more decisively rejected. Remainder preferring NB1 and Leading Edge preferring NB2 would support the two-regime claim.",
    x = NULL,
    y = "Delta AIC (0 = preferred model for that set)",
    fill = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(
  file.path(OUT_ROOT, "Figure_NB_Regime_AIC.png"),
  p_nb_aic, width = 13, height = 7, dpi = 300
)

# Diagnostic C figure: for each comparison/k, the null distribution of
# diff_gap under random gene-set partitions, with the real observed
# Leading-Edge-vs-Remainder gap marked as a vertical line. If the real line
# sits far into the right tail of its null histogram, that is direct,
# interpretable evidence that Leading Edge is more NB2-like than Remainder
# by more than chance would produce for a random partition of that size.
gap_null_plot_df <- gap_null_long_table %>%
  dplyr::mutate(
    k_type_label = dplyr::recode(
      k_type,
      paper_reference = "Paper reference (k=5000)",
      empirical_k_star = "Empirical k*"
    )
  )

gap_real_plot_df <- gap_permutation_table %>%
  dplyr::mutate(
    k_type_label = dplyr::recode(
      k_type,
      paper_reference = "Paper reference (k=5000)",
      empirical_k_star = "Empirical k*"
    )
  )

p_gap_perm <- ggplot(gap_null_plot_df, aes(x = null_diff)) +
  geom_histogram(bins = 40, fill = "grey70", color = "white") +
  geom_vline(
    data = gap_real_plot_df,
    aes(xintercept = diff_gap_real),
    color = "#C0392B", linewidth = 1
  ) +
  geom_text(
    data = gap_real_plot_df,
    aes(x = diff_gap_real, y = Inf, label = paste0("p=", signif(perm_p_one_sided, 3))),
    vjust = 1.3, hjust = -0.05, size = 3, color = "#C0392B"
  ) +
  facet_grid(k_type_label ~ comparison, scales = "free") +
  labs(
    title = "Diagnostic C: observed Leading-Edge-vs-Remainder NB2-NB1 gap vs a random-partition null",
    subtitle = "Grey histogram = gap expected from a random gene subset of the same size, repeated many times. Red line = the real, PC1-defined Leading Edge partition's gap.",
    x = "median(gap | subset) - median(gap | rest)",
    y = "Count across permutations"
  ) +
  theme_bw(base_size = 10)

ggsave(
  file.path(OUT_ROOT, "Figure_Gap_Permutation.png"),
  p_gap_perm, width = 13, height = 7, dpi = 300
)

message("Diagnostic complete. Outputs written to: ", OUT_ROOT)
message("  Table_NB_Regime_Fit.csv")
message("  Table_Gap_Permutation.csv")
message("  Figure_NB_Regime_Fit.png")
message("  Figure_NB_Regime_AIC.png")
message("  Figure_Gap_Permutation.png")

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
      figures_zip_path <- file.path(OUT_ROOT, "RegimeDiagnostic_Figures.zip")
      if (file.exists(figures_zip_path)) file.remove(figures_zip_path)
      utils::zip(figures_zip_path, files = figure_files, flags = "-j")
      message("  RegimeDiagnostic_Figures.zip (", length(figure_files), " files)")
    } else {
      message("  No .png files found; skipping RegimeDiagnostic_Figures.zip")
    }

    if (length(table_files) > 0L) {
      tables_zip_path <- file.path(OUT_ROOT, "RegimeDiagnostic_Tables.zip")
      if (file.exists(tables_zip_path)) file.remove(tables_zip_path)
      utils::zip(tables_zip_path, files = table_files, flags = "-j")
      message("  RegimeDiagnostic_Tables.zip (", length(table_files), " files)")
    } else {
      message("  No .csv files found; skipping RegimeDiagnostic_Tables.zip")
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
