#!/usr/bin/env Rscript

# =============================================================================
# REGIME DIAGNOSTIC: EVS LEADING-EDGE / REMAINDER VALIDATION
# =============================================================================
#
# This script is standalone. It does not modify, overwrite, or depend on
# SEQUENCE.R or EMPERICALCUTOFF outputs. It reads the same raw count file
# and reuses the same PC1-ranking convention (log1p-CPM, arm-specific PCA,
# ranked by ascending |PC1 loading|) so results are directly comparable to
# the existing pipeline.
#
# It answers two questions, run under BOTH candidate cutoffs
# (k = 5000 and k = 4077, the empirically-derived global k*):
#
# DIAGNOSTIC A: Is Leading-Edge membership gene-intrinsic, or does it track
# the true condition split?
#   For each comparison, condition labels are randomly permuted many times.
#   Under each permutation, PC1 is recomputed independently within each
#   pseudo-arm and a permuted Leading-Edge set is derived exactly as in the
#   real pipeline. If Leading-Edge membership reflects an intrinsic property
#   of each gene (e.g. a stable dispersion regime), it should be largely
#   insensitive to which samples get called "control" vs "treatment", and
#   permuted Leading-Edge sets should still substantially overlap with the
#   real Leading-Edge set. If Leading-Edge membership is actually tracking
#   the true condition effect, scrambling the labels destroys that signal
#   and permuted Leading-Edge sets should look close to a structureless
#   random baseline.
#
#   Two null references are reported alongside the real-vs-permuted
#   overlap:
#     - label-permutation null: PCA is still run, only labels are shuffled.
#     - fully-random null: no PCA at all; two random k-subsets are drawn
#       directly, as an absolute floor for how much overlap is expected
#       from set size alone.
#
# DIAGNOSTIC B: Does Remainder behave like an NB1 (linear mean-variance)
# regime and Leading Edge like an NB2 (quadratic mean-variance) regime?
#   Within each set, per-feature pooled within-group mean and variance are
#   computed from DESeq2-normalized counts. Excess-over-Poisson variance
#   (variance minus mean, floored at 0) is regressed against mean (NB1,
#   linear) and against mean^2 (NB2, quadratic) separately for Remainder
#   and for Leading Edge. R-squared for each model in each set is reported,
#   directly testing whether Remainder fits NB1 better and Leading Edge
#   fits NB2 better.
#
# OUTPUTS (written to a directory separate from the existing pipeline):
#   Table_Permutation_Stability.csv
#   Table_NB_Regime_Fit.csv
#   Figure_Permutation_Stability.png
#   Figure_NB_Regime_Fit.png
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

CANDIDATE_K <- c(5000L, 4077L)
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

jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0L) return(NA_real_)
  inter / uni
}

chance_jaccard <- function(a_size, b_size, n_total) {
  exp_inter <- (a_size * b_size) / n_total
  exp_union <- a_size + b_size - exp_inter
  if (exp_union <= 0) return(NA_real_)
  exp_inter / exp_union
}

# -----------------------------------------------------------------------------
# DIAGNOSTIC B helpers: NB1 (linear) vs NB2 (quadratic) mean-variance fit
# -----------------------------------------------------------------------------

deseq2_normalize <- function(count_mat, group_labels) {
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
  DESeq2::counts(dds, normalized = TRUE)
}

pooled_mean_variance <- function(normalized_counts, group_labels) {
  groups <- levels(factor(group_labels))
  sse <- rep(0, nrow(normalized_counts))
  sum_mu <- rep(0, nrow(normalized_counts))
  residual_df <- 0L
  n_groups_used <- 0L

  for (g in groups) {
    idx <- which(group_labels == g)
    if (length(idx) < 2L) next
    xg <- normalized_counts[, idx, drop = FALSE]
    mu_g <- rowMeans(xg)
    resid_g <- sweep(xg, 1L, mu_g, "-")
    sse <- sse + rowSums(resid_g^2)
    sum_mu <- sum_mu + mu_g
    residual_df <- residual_df + length(idx) - 1L
    n_groups_used <- n_groups_used + 1L
  }

  V <- sse / max(residual_df, 1L)
  V[!is.finite(V)] <- 0
  V <- pmax(V, 0)
  mu <- sum_mu / max(n_groups_used, 1L)

  data.frame(
    feature_id = rownames(normalized_counts),
    mean = mu,
    variance = V,
    stringsAsFactors = FALSE
  )
}

fit_nb_regimes <- function(mean_var_df, feature_subset) {
  sub <- mean_var_df[mean_var_df$feature_id %in% feature_subset, , drop = FALSE]
  sub <- sub[is.finite(sub$mean) & is.finite(sub$variance) & sub$mean > 0, , drop = FALSE]
  sub$excess <- pmax(sub$variance - sub$mean, 0)

  if (nrow(sub) < 10L) {
    return(data.frame(model = c("NB1_linear", "NB2_quadratic"), r_squared = NA_real_, n_features = nrow(sub)))
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

  data.frame(
    model = c("NB1_linear", "NB2_quadratic"),
    r_squared = c(r2(nb1_fit, sub$excess), r2(nb2_fit, sub$excess)),
    n_features = nrow(sub)
  )
}

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------

message("Reading count matrix...")
count_mat <- read_count_matrix(COUNT_FILE, GROUP_PATTERNS)
n_total_features <- nrow(count_mat)

stability_rows <- list()
nb_rows <- list()

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

  message(sprintf(
    "[%s] control n=%d, treatment n=%d", comparison_name, n_control, n_treatment
  ))

  comparison_mat <- count_mat[, all_ids, drop = FALSE]

  for (k in CANDIDATE_K) {

    # --- Real (unpermuted) Leading Edge / Remainder ---------------------
    real_le <- leading_edge_set(comparison_mat, control_ids, treatment_ids, k)
    real_remainder <- setdiff(rownames(comparison_mat), real_le)

    # --- DIAGNOSTIC A: permutation stability (both sets, symmetric) ------
    # Both Leading Edge and Remainder retain 100% of the data (they are
    # exact complements of each other); this diagnostic tests each side
    # explicitly rather than reporting only Leading Edge and leaving
    # Remainder's stability implied.
    perm_jaccard_le <- numeric(N_PERMUTATIONS)
    perm_jaccard_rem <- numeric(N_PERMUTATIONS)
    random_jaccard_le <- numeric(N_PERMUTATIONS)
    random_jaccard_rem <- numeric(N_PERMUTATIONS)

    for (p in seq_len(N_PERMUTATIONS)) {
      shuffled <- sample(all_ids)
      pseudo_g1 <- shuffled[seq_len(n_control)]
      pseudo_g2 <- shuffled[(n_control + 1L):length(shuffled)]

      perm_le <- leading_edge_set(comparison_mat, pseudo_g1, pseudo_g2, k)
      perm_remainder <- setdiff(rownames(comparison_mat), perm_le)

      perm_jaccard_le[p] <- jaccard(real_le, perm_le)
      perm_jaccard_rem[p] <- jaccard(real_remainder, perm_remainder)

      random_le <- union(
        sample(rownames(comparison_mat), min(k, n_total_features)),
        sample(rownames(comparison_mat), min(k, n_total_features))
      )
      random_remainder <- setdiff(rownames(comparison_mat), random_le)

      random_jaccard_le[p] <- jaccard(real_le, random_le)
      random_jaccard_rem[p] <- jaccard(real_remainder, random_remainder)
    }

    chance_ref_le <- chance_jaccard(length(real_le), length(real_le), n_total_features)
    chance_ref_rem <- chance_jaccard(length(real_remainder), length(real_remainder), n_total_features)

    make_stability_row <- function(feature_set, real_size, perm_j, random_j, chance_ref) {
      data.frame(
        comparison = comparison_name,
        k = k,
        feature_set = feature_set,
        real_set_size = real_size,
        label_permutation_jaccard_mean = mean(perm_j, na.rm = TRUE),
        label_permutation_jaccard_sd = sd(perm_j, na.rm = TRUE),
        fully_random_jaccard_mean = mean(random_j, na.rm = TRUE),
        analytic_chance_jaccard = chance_ref,
        n_permutations = N_PERMUTATIONS,
        interpretation = ifelse(
          mean(perm_j, na.rm = TRUE) > mean(random_j, na.rm = TRUE) * 2,
          "Membership survives relabeling (supports gene-intrinsic regime)",
          "Membership collapses toward chance under relabeling (tracks condition split)"
        ),
        stringsAsFactors = FALSE
      )
    }

    stability_rows[[length(stability_rows) + 1L]] <- make_stability_row(
      "Leading_Edge", length(real_le), perm_jaccard_le, random_jaccard_le, chance_ref_le
    )
    stability_rows[[length(stability_rows) + 1L]] <- make_stability_row(
      "Remainder", length(real_remainder), perm_jaccard_rem, random_jaccard_rem, chance_ref_rem
    )

    # --- DIAGNOSTIC B: NB1 vs NB2 regime fit -----------------------------
    norm_counts <- deseq2_normalize(comparison_mat, group_labels)
    mean_var_df <- pooled_mean_variance(norm_counts, group_labels)

    fit_le <- fit_nb_regimes(mean_var_df, real_le)
    fit_remainder <- fit_nb_regimes(mean_var_df, real_remainder)

    fit_le$feature_set <- "Leading_Edge"
    fit_remainder$feature_set <- "Remainder"
    fit_le$comparison <- comparison_name
    fit_remainder$comparison <- comparison_name
    fit_le$k <- k
    fit_remainder$k <- k

    nb_rows[[length(nb_rows) + 1L]] <- fit_le
    nb_rows[[length(nb_rows) + 1L]] <- fit_remainder
  }
}

stability_table <- dplyr::bind_rows(stability_rows)
nb_table <- dplyr::bind_rows(nb_rows)

write.csv(
  stability_table,
  file.path(OUT_ROOT, "Table_Permutation_Stability.csv"),
  row.names = FALSE
)
write.csv(
  nb_table,
  file.path(OUT_ROOT, "Table_NB_Regime_Fit.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

stability_long <- stability_table %>%
  dplyr::select(
    comparison, k, feature_set,
    label_permutation_jaccard_mean,
    fully_random_jaccard_mean
  ) %>%
  tidyr::pivot_longer(
    cols = c(label_permutation_jaccard_mean, fully_random_jaccard_mean),
    names_to = "null_type",
    values_to = "jaccard"
  ) %>%
  dplyr::mutate(
    null_type = dplyr::recode(
      null_type,
      label_permutation_jaccard_mean = "Label-permutation null",
      fully_random_jaccard_mean = "Fully-random null"
    ),
    k_label = paste0("k=", k)
  )

p_stability <- ggplot(stability_long, aes(x = comparison, y = jaccard, fill = null_type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(k_label ~ feature_set) +
  labs(
    title = "Diagnostic A: membership stability under label permutation",
    subtitle = "Reported symmetrically for Leading Edge and Remainder, since both retain 100% of the data. Higher label-permutation overlap relative to the fully-random floor supports a gene-intrinsic regime for that set.",
    x = NULL,
    y = "Mean Jaccard overlap with the real (unpermuted) set",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(
  file.path(OUT_ROOT, "Figure_Permutation_Stability.png"),
  p_stability, width = 11, height = 7, dpi = 300
)

nb_plot_df <- nb_table %>%
  dplyr::mutate(k_label = paste0("k=", k))

p_nb <- ggplot(nb_plot_df, aes(x = model, y = r_squared, fill = feature_set)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(k_label ~ comparison) +
  labs(
    title = "Diagnostic B: NB1 (linear) vs NB2 (quadratic) mean-variance fit",
    subtitle = "Remainder fitting NB1 better and Leading Edge fitting NB2 better would support the two-regime claim",
    x = NULL,
    y = expression(R^2),
    fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(
  file.path(OUT_ROOT, "Figure_NB_Regime_Fit.png"),
  p_nb, width = 12, height = 6.5, dpi = 300
)

message("Diagnostic complete. Outputs written to: ", OUT_ROOT)
message("  Table_Permutation_Stability.csv")
message("  Table_NB_Regime_Fit.csv")
message("  Figure_Permutation_Stability.png")
message("  Figure_NB_Regime_Fit.png")
