#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# FINAL MANUSCRIPT SCRIPT
# BOOTSTRAP-STABLE, DATA-DERIVED HIGH-PC1 CUTOFF
# =============================================================================
#
# PURPOSE
#
# This script produces:
#
# 1. Main manuscript figures
#    - Rank: log-transformed CPM-style library-size normalization
#    - Metrics: raw counts
#
# 2. DESeq2 supplementary figures
#    - Rank: DESeq2-normalized counts, optionally VST
#    - Metrics: DESeq2-normalized counts
#
# =============================================================================
# CORE SCIENTIFIC LOGIC
# =============================================================================
#
# A. RANKING
#
# Features are ranked within each arm by absolute PC1 loading:
#
#     low |PC1 loading|  ---------------------->  high |PC1 loading|
#
# Therefore, the RIGHT side of the ranked series contains the features with
# the largest absolute PC1 loadings.
#
# =============================================================================
# B. VARIANCE GEOMETRY
# =============================================================================
#
# At every feature:
#
#     empirical variance
#
# is calculated from the metric matrix and transformed:
#
#     log(1 + variance)
#
# A smoothing spline is fitted to this ranked variance trajectory.
#
# The spline's second derivative is then calculated over a dense rank grid.
# Sign changes in the second derivative identify curvature transitions
# (inflection points).
#
# =============================================================================
# C. NO FIXED 5,000-FEATURE REFERENCE
# =============================================================================
#
# The previous fixed 5,000-feature reference has been removed completely.
#
# Instead, candidate minimum RIGHT-tail fractions are prespecified:
#
#     5%, 7.5%, 10%, 12.5%, 15%
#
# For each candidate minimum fraction f:
#
# 1. Find all spline second-derivative zero crossings.
#
# 2. Restrict to crossings that:
#
#    a. leave at least f * N features from Anchor through the right edge;
#
#    b. leave enough features on the left to construct an equally sized
#       matched LEFT block.
#
# 3. Among those valid crossings, select the RIGHTMOST crossing.
#
# Thus, for each fraction:
#
#     Anchor = rightmost curvature transition satisfying the minimum
#              tail-size requirement.
#
# RIGHT = Anchor through rank N.
#
# LEFT  = immediately preceding block having exactly the same number of
#         features as RIGHT.
#
# =============================================================================
# D. BOOTSTRAP STABILITY
# =============================================================================
#
# The minimum fraction is NOT selected using NB2 results.
#
# Instead, samples within each arm are resampled with replacement.
#
# For each bootstrap replicate:
#
#     sample resampling
#          ->
#     PC1 reranking
#          ->
#     variance spline
#          ->
#     second derivative
#          ->
#     zero crossings
#          ->
#     candidate Anchor
#
# Stability is assessed separately for every candidate minimum fraction.
#
# A fraction is considered stable only if ALL criteria are met:
#
# 1. Bootstrap-valid rate >= MIN_BOOTSTRAP_VALID_RATE
#
# 2. Anchor recovery rate >= MIN_ANCHOR_RECOVERY_RATE
#    Recovery means:
#
#       |bootstrap Anchor - full-data Anchor|
#             <= ANCHOR_TOLERANCE_FRACTION * N
#
# 3. Anchor IQR <= MAX_ANCHOR_IQR_FRACTION * N
#
# 4. Median Jaccard overlap between bootstrap RIGHT-feature membership
#    and full-data RIGHT-feature membership >= MIN_MEDIAN_JACCARD
#
# The SMALLEST candidate minimum fraction satisfying all criteria is chosen.
#
# This means the fraction is chosen from geometric reproducibility alone.
#
# =============================================================================
# E. FINAL CUTOFF
# =============================================================================
#
# After the minimum fraction has been selected by bootstrap stability,
# the final Anchor is the Anchor obtained from the FULL dataset using that
# selected minimum fraction.
#
# The bootstrap median Anchor is NOT substituted for the full-data Anchor.
#
# =============================================================================
# F. NB2 CORROBORATION
# =============================================================================
#
# Only AFTER the Anchor has been selected are NB2-related quantities examined.
#
# Let:
#
#     mu = empirical mean
#     variance = empirical variance
#
# 1. NB2
#
#     log(1 + max(variance - mu, 0))
#
# 2. NB2-NB1
#
#     log(1 + max(variance - mu, 0)) - log(1 + mu)
#
# 3. alpha*mu
#
#     alpha = max((variance - mu) / mu^2, 0)
#
#     alpha*mu signal = log(1 + alpha*mu)
#
# These are descriptive corroborative diagnostics.
#
# They are NOT used to select:
#
#     - minimum fraction
#     - Anchor
#     - spline smoothing
#     - bootstrap stability
#
# Therefore, the geometric selection and NB2 corroboration remain separated.
#
# =============================================================================
# G. METHODS-LEVEL VALIDATION
# =============================================================================
#
# For every arm and track, the script verifies:
#
# - selected fraction passed all stability criteria
# - Anchor is a rounded spline d2 zero crossing
# - RIGHT contains at least the selected minimum fraction
# - LEFT and RIGHT have exactly equal size
# - LEFT does not extend below rank 1
# - summary medians exactly match the plotted/sliced regions
#
# If NO candidate fraction is stable:
#
# - the track is marked UNSTABLE
# - stability tables are still written
# - no cutoff is forced
# - no NB2 corroboration figure is produced for that track
#
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <- "/root/REAPER98632/exports/manuscript_bootstrap_stable"


# -----------------------------------------------------------------------------
# Variance spline
# -----------------------------------------------------------------------------

VAR_SPLINE_SPAR <- 0.60


# -----------------------------------------------------------------------------
# Candidate minimum RIGHT-tail fractions
#
# IMPORTANT:
# These must be specified before examining NB2 corroboration.
# -----------------------------------------------------------------------------

MIN_FRACTION_GRID <- c(
  0.050,
  0.075,
  0.100,
  0.125,
  0.150
)


# -----------------------------------------------------------------------------
# Bootstrap settings
#
# 500 is reasonable for manuscript analysis.
# For development/testing, temporarily reduce to 100-200.
# -----------------------------------------------------------------------------

BOOTSTRAP_N <- 500L

BOOTSTRAP_SEED_BASE <- 20260820L


# -----------------------------------------------------------------------------
# Stability thresholds
#
# These should be fixed before examining NB2 results.
# -----------------------------------------------------------------------------

MIN_BOOTSTRAP_VALID_RATE <- 0.90

ANCHOR_TOLERANCE_FRACTION <- 0.03

MIN_ANCHOR_RECOVERY_RATE <- 0.80

MAX_ANCHOR_IQR_FRACTION <- 0.03

MIN_MEDIAN_JACCARD <- 0.80


# -----------------------------------------------------------------------------
# Figure settings
# -----------------------------------------------------------------------------

PNG_WIDTH_IN  <- 14
PNG_HEIGHT_IN <- 10.8
PNG_DPI       <- 260


# -----------------------------------------------------------------------------
# DESeq2 supplementary analysis
# -----------------------------------------------------------------------------

RUN_DESEQ2_SUPPLEMENT <- TRUE

DESEQ2_RANK_METHOD <- "normalized_log1p"
# Options:
#   "normalized_log1p"
#   "vst"


# -----------------------------------------------------------------------------
# Comparisons
# -----------------------------------------------------------------------------

COMPARISONS <- list(
  RT0_ZT6  = list(control = "^R0_", treatment = "^ZT6_"),
  RT2_ZT8  = list(control = "^R2_", treatment = "^ZT8_"),
  RT4_ZT10 = list(control = "^R4_", treatment = "^ZT10_"),
  RT8_ZT14 = list(control = "^R8_", treatment = "^ZT14_")
)


dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# COLORS
# =============================================================================

COL <- list(

  var_curve = "#117A65",

  nb2      = "#1B9E77",
  nb_gap   = "#CC1E8C",
  alpha_mu = "#386CB0",

  left_fill  = "#CBE3F8",
  right_fill = "#DDF2D5",

  anchor = "#000000",

  zero_crossing = "#777777",

  left_pt  = "#5B8FD1",
  right_pt = "#43A047"
)


TRACE_LEVELS <- c(
  "NB2",
  "NB2-NB1",
  "alpha*mu"
)

TRACE_COLORS <- c(
  "NB2"      = COL$nb2,
  "NB2-NB1"  = COL$nb_gap,
  "alpha*mu" = COL$alpha_mu
)


REGION_LEVELS <- c(
  "LEFT",
  "RIGHT"
)

REGION_COLORS <- c(
  "LEFT"  = COL$left_pt,
  "RIGHT" = COL$right_pt
)


# =============================================================================
# BASIC HELPERS
# =============================================================================

first_numeric_col_index <- function(df) {

  idx <- which(
    vapply(
      df,
      is.numeric,
      logical(1)
    )
  )

  if (length(idx) == 0L) {
    return(NA_integer_)
  }

  idx[1L]
}


# -----------------------------------------------------------------------------

read_count_matrix <- function(path) {

  raw_df <- read.csv(
    path,
    check.names = FALSE
  )

  if (
    nrow(raw_df) == 0L ||
    ncol(raw_df) < 2L
  ) {
    stop(
      "Count file is empty or malformed: ",
      path
    )
  }

  first_num <- first_numeric_col_index(raw_df)

  if (is.na(first_num)) {
    stop("No numeric count columns detected.")
  }

  feature_ids <- make.unique(
    as.character(raw_df[[1L]])
  )

  count_df <- raw_df[
    ,
    first_num:ncol(raw_df),
    drop = FALSE
  ]

  count_mat <- as.matrix(count_df)

  storage.mode(count_mat) <- "numeric"

  rownames(count_mat) <- feature_ids

  count_mat[
    !is.finite(count_mat)
  ] <- 0

  count_mat <- pmax(
    count_mat,
    0
  )

  keep <- rowSums(
    count_mat
  ) > 0

  count_mat[
    keep,
    ,
    drop = FALSE
  ]
}


# -----------------------------------------------------------------------------

normalize_cpm_log1p <- function(count_mat_arm) {

  lib_sizes <- colSums(
    count_mat_arm,
    na.rm = TRUE
  )

  lib_sizes[
    !is.finite(lib_sizes) |
      lib_sizes <= 0
  ] <- 1

  cpm <- sweep(
    count_mat_arm,
    2,
    lib_sizes / 1e6,
    "/"
  )

  log1p(cpm)
}


# =============================================================================
# DESEQ2
# =============================================================================

compute_deseq2_matrices <- function(
    count_mat_arm,
    rank_method = "normalized_log1p") {

  if (
    !requireNamespace(
      "DESeq2",
      quietly = TRUE
    )
  ) {
    return(NULL)
  }

  if (
    !requireNamespace(
      "SummarizedExperiment",
      quietly = TRUE
    )
  ) {
    return(NULL)
  }

  col_data <- data.frame(

    row.names = colnames(count_mat_arm),

    intercept = factor(
      rep(
        "one",
        ncol(count_mat_arm)
      )
    )
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(

    countData = round(count_mat_arm),

    colData = col_data,

    design = ~ 1
  )

  dds <- DESeq2::estimateSizeFactors(dds)

  norm_counts <- DESeq2::counts(
    dds,
    normalized = TRUE
  )

  vst_mat <- NULL

  rank_method_used <- rank_method

  if (rank_method == "vst") {

    vst_obj <- tryCatch(

      DESeq2::vst(
        dds,
        blind = TRUE
      ),

      error = function(e) NULL
    )

    if (!is.null(vst_obj)) {

      vst_mat <- SummarizedExperiment::assay(
        vst_obj
      )

    } else {

      rank_method_used <- "normalized_log1p_fallback"
    }
  }

  ranking_matrix <- switch(

    rank_method,

    normalized_log1p =
      log1p(norm_counts),

    vst =
      if (!is.null(vst_mat)) {
        vst_mat
      } else {
        log1p(norm_counts)
      },

    log1p(norm_counts)
  )

  list(

    normalized_counts = norm_counts,

    ranking_matrix = ranking_matrix,

    size_factors = DESeq2::sizeFactors(dds),

    rank_method_used = rank_method_used
  )
}


# =============================================================================
# PC1 RANKING
# =============================================================================

compute_abs_pc1_loadings <- function(norm_mat_arm) {

  if (
    nrow(norm_mat_arm) < 2L ||
    ncol(norm_mat_arm) < 2L
  ) {
    return(NULL)
  }

  pca <- tryCatch(

    stats::prcomp(
      t(norm_mat_arm),
      center = TRUE,
      scale. = FALSE,
      rank. = 1
    ),

    error = function(e) NULL
  )

  if (
    is.null(pca) ||
    is.null(pca$rotation) ||
    ncol(pca$rotation) < 1L
  ) {
    return(NULL)
  }

  out <- abs(
    pca$rotation[, 1L]
  )

  out[
    !is.finite(out)
  ] <- 0

  out
}


# =============================================================================
# VARIANCE GEOMETRY
# =============================================================================

compute_ranked_variance_curve <- function(
    metric_mat_arm,
    rank_order,
    spar = 0.60) {

  empirical_var <- apply(
    metric_mat_arm,
    1L,
    stats::var,
    na.rm = TRUE
  )

  empirical_var[
    !is.finite(empirical_var)
  ] <- 0

  empirical_var <- pmax(
    empirical_var,
    0
  )

  ranked_var <- empirical_var[
    rank_order
  ]

  ranked_log_var <- log1p(
    ranked_var
  )

  ranks <- seq_along(
    rank_order
  )

  spline_fit <- stats::smooth.spline(

    x = ranks,

    y = ranked_log_var,

    spar = spar
  )

  smooth_y <- as.numeric(
    stats::predict(
      spline_fit,
      x = ranks,
      deriv = 0
    )$y
  )

  smooth_d2 <- as.numeric(
    stats::predict(
      spline_fit,
      x = ranks,
      deriv = 2
    )$y
  )

  dense_n <- max(
    5000L,
    length(ranks) * 4L
  )

  dense_x <- seq(
    min(ranks),
    max(ranks),
    length.out = dense_n
  )

  dense_y <- as.numeric(
    stats::predict(
      spline_fit,
      x = dense_x,
      deriv = 0
    )$y
  )

  dense_d2 <- as.numeric(
    stats::predict(
      spline_fit,
      x = dense_x,
      deriv = 2
    )$y
  )

  out <- data.frame(

    rank = ranks,

    empirical_variance =
      ranked_var,

    log1p_empirical_variance =
      ranked_log_var,

    smooth_log1p_empirical_variance =
      smooth_y,

    d2_spline =
      smooth_d2,

    stringsAsFactors = FALSE
  )

  attr(
    out,
    "dense_curve_df"
  ) <- data.frame(

    dense_rank =
      dense_x,

    dense_smooth_log1p_empirical_variance =
      dense_y,

    dense_d2 =
      dense_d2,

    stringsAsFactors = FALSE
  )

  out
}


# -----------------------------------------------------------------------------

find_d2_zero_crossings <- function(dense_df) {

  x <- dense_df$dense_rank

  y <- dense_df$dense_d2

  ok <- (
    is.finite(x) &
      is.finite(y)
  )

  x <- x[ok]

  y <- y[ok]

  if (length(x) < 2L) {

    return(
      data.frame(
        crossing_rank = numeric(0),
        crossing_type = character(0),
        stringsAsFactors = FALSE
      )
    )
  }

  crossings <- numeric(0)

  for (
    i in seq_len(
      length(x) - 1L
    )
  ) {

    a <- y[i]

    b <- y[i + 1L]

    xa <- x[i]

    xb <- x[i + 1L]

    if (
      !is.finite(a) ||
      !is.finite(b)
    ) {
      next
    }

    # Exact zero on left point
    if (
      a == 0 &&
      b != 0
    ) {

      crossings <- c(
        crossings,
        xa
      )

      next
    }

    # Exact zero on right point
    if (
      b == 0 &&
      a != 0
    ) {

      crossings <- c(
        crossings,
        xb
      )

      next
    }

    # True sign change
    if (
      (a < 0 && b > 0) ||
      (a > 0 && b < 0)
    ) {

      frac <- abs(a) /
        (
          abs(a) +
            abs(b)
        )

      xr <- xa +
        frac *
        (
          xb - xa
        )

      crossings <- c(
        crossings,
        xr
      )
    }
  }

  crossings <- crossings[
    is.finite(crossings)
  ]

  if (
    length(crossings) == 0L
  ) {

    return(
      data.frame(
        crossing_rank = numeric(0),
        crossing_type = character(0),
        stringsAsFactors = FALSE
      )
    )
  }

  crossings <- sort(
    unique(
      round(
        crossings,
        8
      )
    )
  )

  data.frame(

    crossing_rank =
      crossings,

    crossing_type =
      "d2_sign_change",

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# COMPLETE GEOMETRY FOR ONE DATA MATRIX
# =============================================================================

compute_geometry <- function(
    rank_matrix,
    metric_matrix,
    spar = VAR_SPLINE_SPAR) {

  if (
    nrow(rank_matrix) !=
      nrow(metric_matrix)
  ) {

    stop(
      "rank_matrix and metric_matrix have different numbers of features."
    )
  }

  if (
    ncol(rank_matrix) !=
      ncol(metric_matrix)
  ) {

    stop(
      "rank_matrix and metric_matrix have different numbers of samples."
    )
  }

  abs_loadings <- compute_abs_pc1_loadings(
    rank_matrix
  )

  if (is.null(abs_loadings)) {

    stop(
      "PC1 calculation failed."
    )
  }

  rank_order <- order(
    abs_loadings,
    decreasing = FALSE
  )

  variance_df <- compute_ranked_variance_curve(

    metric_mat_arm =
      metric_matrix,

    rank_order =
      rank_order,

    spar =
      spar
  )

  dense_curve_df <- attr(
    variance_df,
    "dense_curve_df"
  )

  zero_df <- find_d2_zero_crossings(
    dense_curve_df
  )

  list(

    abs_loadings =
      abs_loadings,

    rank_order =
      rank_order,

    variance_df =
      variance_df,

    dense_curve_df =
      dense_curve_df,

    zero_df =
      zero_df
  )
}


# =============================================================================
# SELECT ANCHOR FOR ONE MINIMUM FRACTION
# =============================================================================

select_anchor_for_min_fraction <- function(
    zero_df,
    total_n,
    min_fraction) {

  if (
    !is.finite(min_fraction) ||
    min_fraction <= 0 ||
    min_fraction >= 0.50
  ) {

    stop(
      "min_fraction must be > 0 and < 0.50."
    )
  }

  minimum_right_n <- ceiling(
    min_fraction *
      total_n
  )

  # Equal-sized LEFT requires:
  #
  # anchor - 1 >= total_n - anchor + 1
  #
  # therefore:
  #
  # anchor >= (total_n + 2) / 2

  minimum_anchor_for_match <- ceiling(
    (
      total_n + 2L
    ) / 2
  )

  # To leave at least minimum_right_n features on the right:
  #
  # total_n - anchor + 1 >= minimum_right_n
  #
  # anchor <= total_n - minimum_right_n + 1

  maximum_anchor_for_tail <- (
    total_n -
      minimum_right_n +
      1L
  )

  if (
    maximum_anchor_for_tail <
      minimum_anchor_for_match
  ) {

    return(
      data.frame(

        min_fraction =
          min_fraction,

        valid =
          FALSE,

        anchor =
          NA_integer_,

        right_n =
          NA_integer_,

        left_n =
          NA_integer_,

        actual_right_fraction =
          NA_real_,

        minimum_right_n =
          minimum_right_n,

        minimum_anchor_for_match =
          minimum_anchor_for_match,

        maximum_anchor_for_tail =
          maximum_anchor_for_tail,

        stringsAsFactors = FALSE
      )
    )
  }

  if (
    nrow(zero_df) == 0L
  ) {

    return(
      data.frame(

        min_fraction =
          min_fraction,

        valid =
          FALSE,

        anchor =
          NA_integer_,

        right_n =
          NA_integer_,

        left_n =
          NA_integer_,

        actual_right_fraction =
          NA_real_,

        minimum_right_n =
          minimum_right_n,

        minimum_anchor_for_match =
          minimum_anchor_for_match,

        maximum_anchor_for_tail =
          maximum_anchor_for_tail,

        stringsAsFactors = FALSE
      )
    )
  }

  crossing_ranks <- sort(
    unique(
      as.integer(
        round(
          zero_df$crossing_rank
        )
      )
    )
  )

  crossing_ranks <- crossing_ranks[
    crossing_ranks >= 1L &
      crossing_ranks <= total_n
  ]

  eligible <- crossing_ranks[
    crossing_ranks >=
      minimum_anchor_for_match &
      crossing_ranks <=
      maximum_anchor_for_tail
  ]

  if (
    length(eligible) == 0L
  ) {

    return(
      data.frame(

        min_fraction =
          min_fraction,

        valid =
          FALSE,

        anchor =
          NA_integer_,

        right_n =
          NA_integer_,

        left_n =
          NA_integer_,

        actual_right_fraction =
          NA_real_,

        minimum_right_n =
          minimum_right_n,

        minimum_anchor_for_match =
          minimum_anchor_for_match,

        maximum_anchor_for_tail =
          maximum_anchor_for_tail,

        stringsAsFactors = FALSE
      )
    )
  }

  # Critical rule:
  #
  # Choose the RIGHTMOST curvature crossing that still leaves
  # the required minimum RIGHT-tail size.

  anchor <- max(
    eligible
  )

  right_n <- (
    total_n -
      anchor +
      1L
  )

  left_n <- right_n

  actual_right_fraction <- (
    right_n /
      total_n
  )

  data.frame(

    min_fraction =
      min_fraction,

    valid =
      TRUE,

    anchor =
      anchor,

    right_n =
      right_n,

    left_n =
      left_n,

    actual_right_fraction =
      actual_right_fraction,

    minimum_right_n =
      minimum_right_n,

    minimum_anchor_for_match =
      minimum_anchor_for_match,

    maximum_anchor_for_tail =
      maximum_anchor_for_tail,

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# JACCARD FEATURE-SET STABILITY
# =============================================================================

jaccard_similarity <- function(
    set_a,
    set_b) {

  set_a <- unique(set_a)

  set_b <- unique(set_b)

  union_n <- length(
    union(
      set_a,
      set_b
    )
  )

  if (
    union_n == 0L
  ) {
    return(NA_real_)
  }

  intersection_n <- length(
    intersect(
      set_a,
      set_b
    )
  )

  intersection_n /
    union_n
}


# =============================================================================
# DETERMINISTIC SEED
# =============================================================================

seed_from_label <- function(
    label,
    base_seed = BOOTSTRAP_SEED_BASE) {

  label_value <- sum(
    utf8ToInt(
      as.character(label)
    )
  )

  out <- (
    base_seed +
      label_value * 1009L
  ) %% 2147483647L

  as.integer(out)
}


# =============================================================================
# BOOTSTRAP STABILITY
# =============================================================================

assess_geometry_stability <- function(
    rank_matrix,
    metric_matrix,
    candidate_fractions = MIN_FRACTION_GRID,
    bootstrap_n = BOOTSTRAP_N,
    spar = VAR_SPLINE_SPAR,
    seed = BOOTSTRAP_SEED_BASE) {

  candidate_fractions <- sort(
    unique(
      as.numeric(
        candidate_fractions
      )
    )
  )

  if (
    any(
      !is.finite(candidate_fractions)
    ) ||
    any(
      candidate_fractions <= 0
    ) ||
    any(
      candidate_fractions >= 0.50
    )
  ) {

    stop(
      "All candidate fractions must lie between 0 and 0.50."
    )
  }

  if (
    bootstrap_n < 1L
  ) {

    stop(
      "bootstrap_n must be >= 1."
    )
  }

  if (
    is.null(
      rownames(rank_matrix)
    )
  ) {

    stop(
      "rank_matrix must have feature row names."
    )
  }

  if (
    !identical(
      rownames(rank_matrix),
      rownames(metric_matrix)
    )
  ) {

    stop(
      "Feature row names do not match between rank and metric matrices."
    )
  }

  total_n <- nrow(
    rank_matrix
  )

  sample_n <- ncol(
    rank_matrix
  )


  # ---------------------------------------------------------------------------
  # Full-data geometry
  # ---------------------------------------------------------------------------

  original_geometry <- compute_geometry(

    rank_matrix =
      rank_matrix,

    metric_matrix =
      metric_matrix,

    spar =
      spar
  )


  original_candidates <- bind_rows(

    lapply(
      candidate_fractions,
      function(f) {

        select_anchor_for_min_fraction(

          zero_df =
            original_geometry$zero_df,

          total_n =
            total_n,

          min_fraction =
            f
        )
      }
    )
  )


  # ---------------------------------------------------------------------------
  # Full-data RIGHT-feature sets for each candidate fraction
  # ---------------------------------------------------------------------------

  original_right_sets <- vector(
    mode = "list",
    length = length(candidate_fractions)
  )

  names(original_right_sets) <- sprintf(
    "%.6f",
    candidate_fractions
  )

  for (
    j in seq_along(candidate_fractions)
  ) {

    candidate_row <- original_candidates[
      j,
      ,
      drop = FALSE
    ]

    key <- sprintf(
      "%.6f",
      candidate_fractions[j]
    )

    if (
      isTRUE(
        candidate_row$valid
      )
    ) {

      anchor <- candidate_row$anchor

      right_order_idx <- original_geometry$rank_order[
        seq.int(
          anchor,
          total_n
        )
      ]

      original_right_sets[[key]] <- rownames(
        rank_matrix
      )[
        right_order_idx
      ]

    } else {

      original_right_sets[[key]] <- character(0)
    }
  }


  # ---------------------------------------------------------------------------
  # Bootstrap
  # ---------------------------------------------------------------------------

  set.seed(seed)

  boot_rows <- vector(
    mode = "list",
    length = bootstrap_n *
      length(candidate_fractions)
  )

  row_counter <- 0L


  add_invalid_rows <- function(
      bootstrap_id,
      reason) {

    out <- vector(
      "list",
      length(candidate_fractions)
    )

    for (
      jj in seq_along(candidate_fractions)
    ) {

      out[[jj]] <- data.frame(

        bootstrap =
          bootstrap_id,

        min_fraction =
          candidate_fractions[jj],

        valid =
          FALSE,

        anchor =
          NA_integer_,

        right_n =
          NA_integer_,

        right_fraction =
          NA_real_,

        jaccard =
          NA_real_,

        invalid_reason =
          reason,

        stringsAsFactors = FALSE
      )
    }

    out
  }


  for (
    b in seq_len(
      bootstrap_n
    )
  ) {

    boot_sample_idx <- sample.int(

      n =
        sample_n,

      size =
        sample_n,

      replace =
        TRUE
    )


    # A bootstrap consisting of only one unique biological sample
    # contains no meaningful between-sample variance information.

    if (
      length(
        unique(
          boot_sample_idx
        )
      ) < 2L
    ) {

      invalid_list <- add_invalid_rows(
        bootstrap_id = b,
        reason = "fewer_than_2_unique_samples"
      )

      for (
        jj in seq_along(invalid_list)
      ) {

        row_counter <- row_counter + 1L

        boot_rows[[row_counter]] <-
          invalid_list[[jj]]
      }

      next
    }


    boot_rank_matrix <- rank_matrix[
      ,
      boot_sample_idx,
      drop = FALSE
    ]

    boot_metric_matrix <- metric_matrix[
      ,
      boot_sample_idx,
      drop = FALSE
    ]


    boot_geometry <- tryCatch(

      compute_geometry(

        rank_matrix =
          boot_rank_matrix,

        metric_matrix =
          boot_metric_matrix,

        spar =
          spar
      ),

      error = function(e) NULL
    )


    if (
      is.null(
        boot_geometry
      )
    ) {

      invalid_list <- add_invalid_rows(
        bootstrap_id = b,
        reason = "geometry_failed"
      )

      for (
        jj in seq_along(invalid_list)
      ) {

        row_counter <- row_counter + 1L

        boot_rows[[row_counter]] <-
          invalid_list[[jj]]
      }

      next
    }


    for (
      j in seq_along(
        candidate_fractions
      )
    ) {

      f <- candidate_fractions[j]

      selection <- select_anchor_for_min_fraction(

        zero_df =
          boot_geometry$zero_df,

        total_n =
          total_n,

        min_fraction =
          f
      )

      row_counter <- row_counter + 1L


      if (
        !isTRUE(
          selection$valid
        )
      ) {

        boot_rows[[row_counter]] <- data.frame(

          bootstrap =
            b,

          min_fraction =
            f,

          valid =
            FALSE,

          anchor =
            NA_integer_,

          right_n =
            NA_integer_,

          right_fraction =
            NA_real_,

          jaccard =
            NA_real_,

          invalid_reason =
            "no_valid_anchor",

          stringsAsFactors = FALSE
        )

        next
      }


      boot_anchor <- selection$anchor

      boot_right_order_idx <- boot_geometry$rank_order[
        seq.int(
          boot_anchor,
          total_n
        )
      ]

      boot_right_features <- rownames(
        rank_matrix
      )[
        boot_right_order_idx
      ]

      key <- sprintf(
        "%.6f",
        f
      )

      original_right_features <-
        original_right_sets[[key]]


      jac <- if (
        length(
          original_right_features
        ) > 0L
      ) {

        jaccard_similarity(
          original_right_features,
          boot_right_features
        )

      } else {

        NA_real_
      }


      boot_rows[[row_counter]] <- data.frame(

        bootstrap =
          b,

        min_fraction =
          f,

        valid =
          TRUE,

        anchor =
          boot_anchor,

        right_n =
          selection$right_n,

        right_fraction =
          selection$actual_right_fraction,

        jaccard =
          jac,

        invalid_reason =
          NA_character_,

        stringsAsFactors = FALSE
      )
    }
  }


  boot_rows <- boot_rows[
    seq_len(
      row_counter
    )
  ]

  bootstrap_df <- bind_rows(
    boot_rows
  )


  # ---------------------------------------------------------------------------
  # Stability summary
  # ---------------------------------------------------------------------------

  stability_rows <- vector(
    mode = "list",
    length = length(candidate_fractions)
  )


  for (
    j in seq_along(
      candidate_fractions
    )
  ) {

    f <- candidate_fractions[j]

    original_row <- original_candidates[
      j,
      ,
      drop = FALSE
    ]

    boot_sub <- bootstrap_df[
      abs(
        bootstrap_df$min_fraction -
          f
      ) < 1e-12,
      ,
      drop = FALSE
    ]

    valid_idx <- which(
      boot_sub$valid
    )

    valid_rate <- mean(
      boot_sub$valid
    )


    if (
      length(valid_idx) > 0L
    ) {

      valid_anchors <- boot_sub$anchor[
        valid_idx
      ]

      valid_right_fraction <- boot_sub$right_fraction[
        valid_idx
      ]

      anchor_median <- median(
        valid_anchors,
        na.rm = TRUE
      )

      anchor_iqr <- stats::IQR(
        valid_anchors,
        na.rm = TRUE
      )

      anchor_iqr_fraction <- (
        anchor_iqr /
          total_n
      )

      right_fraction_median <- median(
        valid_right_fraction,
        na.rm = TRUE
      )

      right_fraction_iqr <- stats::IQR(
        valid_right_fraction,
        na.rm = TRUE
      )

    } else {

      anchor_median <- NA_real_

      anchor_iqr <- NA_real_

      anchor_iqr_fraction <- NA_real_

      right_fraction_median <- NA_real_

      right_fraction_iqr <- NA_real_
    }


    if (
      isTRUE(
        original_row$valid
      )
    ) {

      full_anchor <- original_row$anchor

      tolerance_n <- ceiling(
        ANCHOR_TOLERANCE_FRACTION *
          total_n
      )

      recovered <- rep(
        FALSE,
        nrow(boot_sub)
      )

      recovered[
        boot_sub$valid
      ] <- abs(
        boot_sub$anchor[
          boot_sub$valid
        ] -
          full_anchor
      ) <= tolerance_n

      # Denominator is ALL bootstrap replicates.
      # Failed geometry therefore does not count as recovery.

      anchor_recovery_rate <- mean(
        recovered
      )

    } else {

      full_anchor <- NA_integer_

      tolerance_n <- ceiling(
        ANCHOR_TOLERANCE_FRACTION *
          total_n
      )

      anchor_recovery_rate <- NA_real_
    }


    valid_jaccard <- boot_sub$jaccard[
      is.finite(
        boot_sub$jaccard
      )
    ]

    median_jaccard <- if (
      length(
        valid_jaccard
      ) > 0L
    ) {

      median(
        valid_jaccard
      )

    } else {

      NA_real_
    }


    passes <- (
      isTRUE(
        original_row$valid
      ) &&
        is.finite(valid_rate) &&
        valid_rate >=
          MIN_BOOTSTRAP_VALID_RATE &&
        is.finite(
          anchor_recovery_rate
        ) &&
        anchor_recovery_rate >=
          MIN_ANCHOR_RECOVERY_RATE &&
        is.finite(
          anchor_iqr_fraction
        ) &&
        anchor_iqr_fraction <=
          MAX_ANCHOR_IQR_FRACTION &&
        is.finite(
          median_jaccard
        ) &&
        median_jaccard >=
          MIN_MEDIAN_JACCARD
    )


    stability_rows[[j]] <- data.frame(

      min_fraction =
        f,

      original_valid =
        original_row$valid,

      anchor =
        original_row$anchor,

      right_n =
        original_row$right_n,

      actual_right_fraction =
        original_row$actual_right_fraction,

      bootstrap_n =
        bootstrap_n,

      valid_bootstraps =
        length(valid_idx),

      bootstrap_valid_rate =
        valid_rate,

      anchor_median =
        anchor_median,

      anchor_iqr =
        anchor_iqr,

      anchor_iqr_fraction =
        anchor_iqr_fraction,

      anchor_tolerance_n =
        tolerance_n,

      anchor_recovery_rate =
        anchor_recovery_rate,

      right_fraction_median =
        right_fraction_median,

      right_fraction_iqr =
        right_fraction_iqr,

      median_jaccard =
        median_jaccard,

      pass_stability =
        passes,

      stringsAsFactors = FALSE
    )
  }


  stability_summary <- bind_rows(
    stability_rows
  )


  # ---------------------------------------------------------------------------
  # Select the SMALLEST stable minimum fraction
  # ---------------------------------------------------------------------------

  stable_rows <- stability_summary %>%

    filter(
      pass_stability
    ) %>%

    arrange(
      min_fraction
    )


  selected <- if (
    nrow(stable_rows) > 0L
  ) {

    stable_rows[
      1L,
      ,
      drop = FALSE
    ]

  } else {

    NULL
  }


  list(

    original_geometry =
      original_geometry,

    original_candidates =
      original_candidates,

    bootstrap_df =
      bootstrap_df,

    stability_summary =
      stability_summary,

    selected =
      selected
  )
}


# =============================================================================
# FEATURE-LEVEL NB2 METRICS
# =============================================================================

compute_ranked_feature_metrics <- function(
    metric_mat_arm,
    rank_order) {

  ranked_mat <- metric_mat_arm[
    rank_order,
    ,
    drop = FALSE
  ]

  mu <- rowMeans(
    ranked_mat,
    na.rm = TRUE
  )

  empirical_var <- apply(
    ranked_mat,
    1L,
    stats::var,
    na.rm = TRUE
  )

  mu[
    !is.finite(mu)
  ] <- 0

  empirical_var[
    !is.finite(empirical_var)
  ] <- 0

  mu <- pmax(
    mu,
    0
  )

  empirical_var <- pmax(
    empirical_var,
    0
  )


  # Extra-Poisson variance component

  nb2_var_minus_mu <- pmax(
    empirical_var -
      mu,
    0
  )


  # Method-of-moments alpha estimate under:
  #
  # variance = mu + alpha*mu^2

  alpha_hat <- rep(
    0,
    length(mu)
  )

  pos <- mu > 0

  alpha_hat[pos] <- pmax(

    (
      empirical_var[pos] -
        mu[pos]
    ) /
      (
        mu[pos]^2
      ),

    0
  )

  alpha_mu_val <- (
    alpha_hat *
      mu
  )


  data.frame(

    rank =
      seq_along(rank_order),

    feature_id =
      rownames(metric_mat_arm)[
        rank_order
      ],

    mu =
      mu,

    empirical_variance =
      empirical_var,

    NB2 =
      log1p(
        nb2_var_minus_mu
      ),

    NB2_NB1 =
      log1p(
        nb2_var_minus_mu
      ) -
      log1p(mu),

    alpha_mu =
      log1p(
        alpha_mu_val
      ),

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# MATCHED LEFT / RIGHT SUMMARY
# =============================================================================

summarize_regions <- function(
    feature_df,
    anchor_rank,
    total_n) {

  right_idx <- seq.int(
    anchor_rank,
    total_n
  )

  right_n <- length(
    right_idx
  )

  left_end <- (
    anchor_rank -
      1L
  )

  left_start <- (
    left_end -
      right_n +
      1L
  )

  if (
    left_start < 1L
  ) {

    stop(
      "LEFT block extends below rank 1."
    )
  }

  left_idx <- seq.int(
    left_start,
    left_end
  )

  left_df <- feature_df[
    left_idx,
    ,
    drop = FALSE
  ]

  right_df <- feature_df[
    right_idx,
    ,
    drop = FALSE
  ]


  data.frame(

    left_start =
      left_start,

    left_end =
      left_end,

    right_start =
      anchor_rank,

    right_end =
      total_n,

    left_n =
      nrow(left_df),

    right_n =
      nrow(right_df),


    left_NB2 =
      median(
        left_df$NB2,
        na.rm = TRUE
      ),

    right_NB2 =
      median(
        right_df$NB2,
        na.rm = TRUE
      ),

    diff_NB2 =
      median(
        right_df$NB2,
        na.rm = TRUE
      ) -
      median(
        left_df$NB2,
        na.rm = TRUE
      ),


    left_gap =
      median(
        left_df$NB2_NB1,
        na.rm = TRUE
      ),

    right_gap =
      median(
        right_df$NB2_NB1,
        na.rm = TRUE
      ),

    diff_gap =
      median(
        right_df$NB2_NB1,
        na.rm = TRUE
      ) -
      median(
        left_df$NB2_NB1,
        na.rm = TRUE
      ),


    left_alpha =
      median(
        left_df$alpha_mu,
        na.rm = TRUE
      ),

    right_alpha =
      median(
        right_df$alpha_mu,
        na.rm = TRUE
      ),

    diff_alpha =
      median(
        right_df$alpha_mu,
        na.rm = TRUE
      ) -
      median(
        left_df$alpha_mu,
        na.rm = TRUE
      ),

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# VALIDATION
# =============================================================================

validate_method_level <- function(
    feature_df,
    zero_df,
    anchor,
    selected_fraction,
    selected_stability_row,
    region_summary,
    total_n) {

  # ---------------------------------------------------------------------------
  # Stability must have passed
  # ---------------------------------------------------------------------------

  if (
    !isTRUE(
      selected_stability_row$pass_stability
    )
  ) {

    stop(
      "Validation failed: selected fraction did not pass stability."
    )
  }


  # ---------------------------------------------------------------------------
  # Anchor validity
  # ---------------------------------------------------------------------------

  if (
    !is.finite(anchor) ||
    anchor < 1L ||
    anchor > total_n
  ) {

    stop(
      "Validation failed: invalid Anchor."
    )
  }


  # ---------------------------------------------------------------------------
  # Anchor must correspond to a rounded d2 zero crossing
  # ---------------------------------------------------------------------------

  rounded_crossings <- unique(
    as.integer(
      round(
        zero_df$crossing_rank
      )
    )
  )

  if (
    !(anchor %in% rounded_crossings)
  ) {

    stop(
      "Validation failed: Anchor is not a rounded d2 zero crossing."
    )
  }


  # ---------------------------------------------------------------------------
  # RIGHT must satisfy the selected minimum fraction
  # ---------------------------------------------------------------------------

  right_idx <- seq.int(
    anchor,
    total_n
  )

  right_n <- length(
    right_idx
  )

  actual_fraction <- (
    right_n /
      total_n
  )

  if (
    actual_fraction +
      1e-12 <
      selected_fraction
  ) {

    stop(
      "Validation failed: RIGHT is smaller than selected minimum fraction."
    )
  }


  # ---------------------------------------------------------------------------
  # LEFT must be equal-sized
  # ---------------------------------------------------------------------------

  left_end <- (
    anchor -
      1L
  )

  left_start <- (
    left_end -
      right_n +
      1L
  )

  if (
    left_start < 1L
  ) {

    stop(
      "Validation failed: LEFT block extends below rank 1."
    )
  }

  left_idx <- seq.int(
    left_start,
    left_end
  )

  if (
    length(left_idx) !=
      length(right_idx)
  ) {

    stop(
      "Validation failed: LEFT and RIGHT blocks are not equal size."
    )
  }


  # ---------------------------------------------------------------------------
  # Verify reported medians exactly match sliced regions
  # ---------------------------------------------------------------------------

  left_df <- feature_df[
    left_idx,
    ,
    drop = FALSE
  ]

  right_df <- feature_df[
    right_idx,
    ,
    drop = FALSE
  ]


  expected <- list(

    left_n =
      nrow(left_df),

    right_n =
      nrow(right_df),

    left_NB2 =
      median(
        left_df$NB2,
        na.rm = TRUE
      ),

    right_NB2 =
      median(
        right_df$NB2,
        na.rm = TRUE
      ),

    diff_NB2 =
      median(
        right_df$NB2,
        na.rm = TRUE
      ) -
      median(
        left_df$NB2,
        na.rm = TRUE
      ),

    left_gap =
      median(
        left_df$NB2_NB1,
        na.rm = TRUE
      ),

    right_gap =
      median(
        right_df$NB2_NB1,
        na.rm = TRUE
      ),

    diff_gap =
      median(
        right_df$NB2_NB1,
        na.rm = TRUE
      ) -
      median(
        left_df$NB2_NB1,
        na.rm = TRUE
      ),

    left_alpha =
      median(
        left_df$alpha_mu,
        na.rm = TRUE
      ),

    right_alpha =
      median(
        right_df$alpha_mu,
        na.rm = TRUE
      ),

    diff_alpha =
      median(
        right_df$alpha_mu,
        na.rm = TRUE
      ) -
      median(
        left_df$alpha_mu,
        na.rm = TRUE
      )
  )


  for (
    nm in names(expected)
  ) {

    chk <- all.equal(

      as.numeric(
        region_summary[[nm]][1]
      ),

      as.numeric(
        expected[[nm]]
      ),

      tolerance = 1e-10
    )

    if (
      !isTRUE(chk)
    ) {

      stop(
        "Validation failed for summary field: ",
        nm,
        " | ",
        chk
      )
    }
  }


  invisible(TRUE)
}


# =============================================================================
# FIGURE HELPERS
# =============================================================================

save_three_panel_plot <- function(
    plot_list,
    filename) {

  png(

    filename,

    width =
      PNG_WIDTH_IN,

    height =
      PNG_HEIGHT_IN,

    units =
      "in",

    res =
      PNG_DPI,

    bg =
      "white"
  )

  grid.newpage()

  pushViewport(

    viewport(

      layout = grid.layout(

        nrow = 3,

        ncol = 1,

        heights = unit(
          c(
            1.18,
            1.18,
            0.88
          ),
          "null"
        )
      )
    )
  )

  for (
    i in seq_along(
      plot_list
    )
  ) {

    print(

      plot_list[[i]],

      vp = viewport(

        layout.pos.row =
          i,

        layout.pos.col =
          1
      )
    )
  }

  dev.off()
}


# -----------------------------------------------------------------------------

compute_text_positions <- function(
    total_n,
    anchor) {

  pad <- max(
    50L,
    round(
      total_n *
        0.04
    )
  )

  safe_right_limit <- max(
    200L,
    anchor -
      pad
  )

  left_x <- max(
    5L,
    round(
      total_n *
        0.035
    )
  )

  stat_x <- min(

    max(
      150L,
      round(
        total_n *
          0.26
      )
    ),

    safe_right_limit
  )

  if (
    stat_x <=
      left_x +
      50L
  ) {

    stat_x <- (
      left_x +
        60L
    )
  }

  list(

    left_x =
      left_x,

    stat_x =
      stat_x
  )
}


# =============================================================================
# FIGURE BUILDER
# =============================================================================

build_main_figure <- function(
    comparison_name,
    arm_name,
    variance_df,
    zero_df,
    feature_df,
    anchor,
    selected_stability_row,
    region_summary,
    out_file,
    fig_tag,
    rank_tag,
    metric_tag) {

  total_n <- nrow(
    feature_df
  )

  left_n <- region_summary$left_n

  left_min <- (
    anchor -
      left_n
  )

  left_max <- (
    anchor -
      1L
  )

  right_min <- anchor

  right_max <- total_n


  # ---------------------------------------------------------------------------
  # Zero-crossing display positions
  # ---------------------------------------------------------------------------

  zero_plot_df <- data.frame(

    rank =
      zero_df$crossing_rank,

    stringsAsFactors = FALSE
  )

  if (
    nrow(zero_plot_df) > 0L
  ) {

    zero_plot_df$y <- approx(

      x =
        variance_df$rank,

      y =
        variance_df$smooth_log1p_empirical_variance,

      xout =
        zero_plot_df$rank,

      rule =
        2
    )$y
  }


  anchor_y <- approx(

    x =
      variance_df$rank,

    y =
      variance_df$smooth_log1p_empirical_variance,

    xout =
      anchor,

    rule =
      2
  )$y


  nb_long <- feature_df %>%

    select(
      rank,
      NB2,
      NB2_NB1,
      alpha_mu
    ) %>%

    pivot_longer(

      cols = c(
        NB2,
        NB2_NB1,
        alpha_mu
      ),

      names_to =
        "metric",

      values_to =
        "value"
    ) %>%

    mutate(

      metric = factor(

        metric,

        levels = c(
          "NB2",
          "NB2_NB1",
          "alpha_mu"
        ),

        labels = c(
          "NB2",
          "NB2-NB1",
          "alpha*mu"
        )
      )
    )


  pos <- compute_text_positions(
    total_n,
    anchor
  )


  top_y <- max(
    variance_df$smooth_log1p_empirical_variance,
    na.rm = TRUE
  )

  mid_y <- max(
    nb_long$value,
    na.rm = TRUE
  )


  selected_fraction_pct <- round(
    selected_stability_row$min_fraction *
      100,
    1
  )

  actual_fraction_pct <- round(
    selected_stability_row$actual_right_fraction *
      100,
    1
  )


  box1 <- paste(

    fig_tag,

    rank_tag,

    metric_tag,

    "Blue = LEFT",

    "Green = RIGHT",

    sep = "\n"
  )


  box2 <- paste0(

    "Anchor = ",
    anchor,
    "\n",

    "Selected minimum = ",
    selected_fraction_pct,
    "%\n",

    "Actual RIGHT = ",
    actual_fraction_pct,
    "%\n",

    "LEFT n = ",
    region_summary$left_n,
    "\n",

    "RIGHT n = ",
    region_summary$right_n,
    "\n",

    "Bootstrap recovery = ",
    round(
      selected_stability_row$anchor_recovery_rate,
      3
    ),
    "\n",

    "Median Jaccard = ",
    round(
      selected_stability_row$median_jaccard,
      3
    )
  )


  box3 <- paste(

    "RIGHT = Anchor to end",

    "LEFT = matched block",

    "Higher RIGHT = more NB2-like",

    sep = "\n"
  )


  box4 <- paste0(

    "LEFT NB2 = ",
    round(
      region_summary$left_NB2,
      3
    ),
    "\n",

    "RIGHT NB2 = ",
    round(
      region_summary$right_NB2,
      3
    ),
    "\n",

    "RIGHT-LEFT NB2 = ",
    round(
      region_summary$diff_NB2,
      3
    ),
    "\n",

    "LEFT NB2-NB1 = ",
    round(
      region_summary$left_gap,
      3
    ),
    "\n",

    "RIGHT NB2-NB1 = ",
    round(
      region_summary$right_gap,
      3
    ),
    "\n",

    "RIGHT-LEFT NB2-NB1 = ",
    round(
      region_summary$diff_gap,
      3
    ),
    "\n",

    "LEFT alpha*mu = ",
    round(
      region_summary$left_alpha,
      3
    ),
    "\n",

    "RIGHT alpha*mu = ",
    round(
      region_summary$right_alpha,
      3
    ),
    "\n",

    "RIGHT-LEFT alpha*mu = ",
    round(
      region_summary$diff_alpha,
      3
    )
  )


  # ---------------------------------------------------------------------------
  # PANEL 1: GEOMETRY
  # ---------------------------------------------------------------------------

  p1 <- ggplot(

    variance_df,

    aes(
      rank,
      smooth_log1p_empirical_variance
    )
  ) +

    annotate(

      "rect",

      xmin =
        left_min,

      xmax =
        left_max,

      ymin =
        -Inf,

      ymax =
        Inf,

      fill =
        COL$left_fill,

      alpha =
        0.70
    ) +

    annotate(

      "rect",

      xmin =
        right_min,

      xmax =
        right_max,

      ymin =
        -Inf,

      ymax =
        Inf,

      fill =
        COL$right_fill,

      alpha =
        0.70
    ) +

    geom_line(

      color =
        COL$var_curve,

      linewidth =
        1.0
    ) +

    geom_point(

      data =
        zero_plot_df,

      aes(
        rank,
        y
      ),

      inherit.aes =
        FALSE,

      color =
        COL$zero_crossing,

      size =
        1.25,

      alpha =
        0.65
    ) +

    geom_vline(

      xintercept =
        anchor,

      color =
        COL$anchor,

      linetype =
        "solid",

      linewidth =
        1.0
    ) +

    geom_point(

      aes(
        x =
          anchor,

        y =
          anchor_y
      ),

      inherit.aes =
        FALSE,

      color =
        COL$anchor,

      shape =
        16,

      size =
        3.5
    ) +

    annotate(

      "label",

      x =
        pos$left_x,

      y =
        top_y *
          0.96,

      label =
        box1,

      hjust =
        0,

      vjust =
        1,

      size =
        2.8,

      label.size =
        0.25,

      fill =
        grDevices::adjustcolor(
          "white",
          alpha.f = 0.96
        )
    ) +

    annotate(

      "label",

      x =
        pos$stat_x,

      y =
        top_y *
          0.70,

      label =
        box2,

      hjust =
        0,

      vjust =
        1,

      size =
        2.8,

      label.size =
        0.25,

      fill =
        grDevices::adjustcolor(
          "white",
          alpha.f = 0.96
        )
    ) +

    labs(

      title =
        paste0(
          comparison_name,
          " ",
          arm_name,
          ": geometry"
        ),

      subtitle =
        "Anchor is the bootstrap-stable right-tail curvature cutoff; grey points are d2 zero crossings",

      x =
        "Rank",

      y =
        "Smoothed log(1 + variance)"
    ) +

    theme_bw(
      base_size = 11
    ) +

    theme(

      panel.grid.minor =
        element_blank(),

      legend.position =
        "none"
    )


  # ---------------------------------------------------------------------------
  # PANEL 2: FEATURE-WISE CORROBORATION
  # ---------------------------------------------------------------------------

  p2 <- ggplot() +

    annotate(

      "rect",

      xmin =
        left_min,

      xmax =
        left_max,

      ymin =
        -Inf,

      ymax =
        Inf,

      fill =
        COL$left_fill,

      alpha =
        0.70
    ) +

    annotate(

      "rect",

      xmin =
        right_min,

      xmax =
        right_max,

      ymin =
        -Inf,

      ymax =
        Inf,

      fill =
        COL$right_fill,

      alpha =
        0.70
    ) +

    geom_vline(

      xintercept =
        anchor,

      color =
        COL$anchor,

      linewidth =
        0.8
    ) +

    geom_line(

      data =
        nb_long,

      aes(
        rank,
        value,
        color = metric
      ),

      linewidth =
        0.95
    ) +

    annotate(

      "label",

      x =
        pos$left_x,

      y =
        mid_y *
          0.96,

      label =
        box3,

      hjust =
        0,

      vjust =
        1,

      size =
        2.8,

      label.size =
        0.25,

      fill =
        grDevices::adjustcolor(
          "white",
          alpha.f = 0.96
        )
    ) +

    annotate(

      "label",

      x =
        pos$stat_x,

      y =
        mid_y *
          0.70,

      label =
        box4,

      hjust =
        0,

      vjust =
        1,

      size =
        2.8,

      label.size =
        0.25,

      fill =
        grDevices::adjustcolor(
          "white",
          alpha.f = 0.96
        )
    ) +

    scale_color_manual(

      values =
        TRACE_COLORS,

      breaks =
        TRACE_LEVELS
    ) +

    labs(

      title =
        paste0(
          comparison_name,
          " ",
          arm_name,
          ": corroboration"
        ),

      subtitle =
        "NB2-related quantities are evaluated only after geometric cutoff selection",

      x =
        "Rank",

      y =
        "NB2-related signal",

      color =
        NULL
    ) +

    theme_bw(
      base_size = 11
    ) +

    theme(

      panel.grid.minor =
        element_blank(),

      legend.position =
        "bottom"
    )


  # ---------------------------------------------------------------------------
  # PANEL 3: REGION SUMMARY
  # ---------------------------------------------------------------------------

  summary_df <- data.frame(

    metric = factor(

      c(
        "NB2",
        "NB2-NB1",
        "alpha*mu"
      ),

      levels = rev(
        c(
          "NB2",
          "NB2-NB1",
          "alpha*mu"
        )
      )
    ),

    LEFT = c(

      region_summary$left_NB2,

      region_summary$left_gap,

      region_summary$left_alpha
    ),

    RIGHT = c(

      region_summary$right_NB2,

      region_summary$right_gap,

      region_summary$right_alpha
    ),

    stringsAsFactors = FALSE
  )


  p3 <- ggplot(

    summary_df,

    aes(
      y = metric
    )
  ) +

    geom_segment(

      aes(
        x =
          LEFT,

        xend =
          RIGHT,

        yend =
          metric
      ),

      color =
        "#7A7A7A",

      linewidth =
        0.8
    ) +

    geom_point(

      aes(
        x =
          LEFT,

        color =
          "LEFT"
      ),

      size =
        3.4
    ) +

    geom_point(

      aes(
        x =
          RIGHT,

        color =
          "RIGHT"
      ),

      size =
        3.4
    ) +

    scale_color_manual(

      values =
        REGION_COLORS,

      breaks =
        REGION_LEVELS
    ) +

    labs(

      title =
        paste0(
          comparison_name,
          " ",
          arm_name,
          ": summary"
        ),

      subtitle =
        "Points farther right indicate stronger NB2-related corroboration",

      x =
        "Median",

      y =
        NULL,

      color =
        NULL
    ) +

    theme_bw(
      base_size = 11
    ) +

    theme(

      panel.grid.minor =
        element_blank(),

      legend.position =
        "bottom"
    )


  save_three_panel_plot(

    plot_list =
      list(
        p1,
        p2,
        p3
      ),

    filename =
      out_file
  )
}


# =============================================================================
# ONE ANALYSIS TRACK
# =============================================================================

run_one_track <- function(
    comparison_name,
    arm_name,
    rank_matrix,
    metric_matrix,
    output_dir,
    track,
    fig_tag,
    rank_tag,
    metric_tag,
    metric_name,
    seed) {

  total_n <- nrow(
    rank_matrix
  )


  # ---------------------------------------------------------------------------
  # Output paths
  # ---------------------------------------------------------------------------

  stability_summary_path <- file.path(

    output_dir,

    paste0(
      "Table_StabilitySummary_",
      comparison_name,
      "_",
      arm_name,
      "_",
      track,
      ".csv"
    )
  )


  stability_bootstrap_path <- file.path(

    output_dir,

    paste0(
      "Table_StabilityBootstrap_",
      comparison_name,
      "_",
      arm_name,
      "_",
      track,
      ".csv"
    )
  )


  zero_path <- file.path(

    output_dir,

    paste0(
      "Table_Zero_",
      comparison_name,
      "_",
      arm_name,
      "_",
      track,
      ".csv"
    )
  )


  valid_path <- file.path(

    output_dir,

    paste0(
      "Table_Valid_",
      comparison_name,
      "_",
      arm_name,
      "_",
      track,
      ".csv"
    )
  )


  rank_path <- file.path(

    output_dir,

    paste0(
      "Table_Rank_",
      comparison_name,
      "_",
      arm_name,
      "_",
      track,
      ".csv"
    )
  )


  cut_path <- file.path(

    output_dir,

    paste0(
      "Table_Cutoff_",
      comparison_name,
      "_",
      arm_name,
      "_",
      track,
      ".csv"
    )
  )


  fig_path <- file.path(

    output_dir,

    paste0(
      "Figure_",
      track,
      "_",
      comparison_name,
      "_",
      arm_name,
      ".png"
    )
  )


  # ---------------------------------------------------------------------------
  # Geometry stability
  # ---------------------------------------------------------------------------

  stability_obj <- assess_geometry_stability(

    rank_matrix =
      rank_matrix,

    metric_matrix =
      metric_matrix,

    candidate_fractions =
      MIN_FRACTION_GRID,

    bootstrap_n =
      BOOTSTRAP_N,

    spar =
      VAR_SPLINE_SPAR,

    seed =
      seed
  )


  write.csv(

    stability_obj$stability_summary,

    stability_summary_path,

    row.names = FALSE
  )


  write.csv(

    stability_obj$bootstrap_df,

    stability_bootstrap_path,

    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # If no stable fraction exists, DO NOT force a cutoff.
  # ---------------------------------------------------------------------------

  if (
    is.null(
      stability_obj$selected
    )
  ) {

    write.csv(

      data.frame(

        comp =
          comparison_name,

        arm =
          arm_name,

        track =
          track,

        Status =
          "UNSTABLE",

        Reason =
          "No candidate minimum fraction passed all bootstrap stability criteria.",

        stringsAsFactors = FALSE
      ),

      valid_path,

      row.names = FALSE
    )


    fail_summary <- data.frame(

      comp =
        comparison_name,

      arm =
        arm_name,

      track =
        track,

      rank_method =
        rank_tag,

      metric_matrix =
        metric_name,

      Status =
        "UNSTABLE",

      SelectedMinFraction =
        NA_real_,

      Anchor =
        NA_integer_,

      LeftSize =
        NA_integer_,

      RightSize =
        NA_integer_,

      left_NB2 =
        NA_real_,

      right_NB2 =
        NA_real_,

      diff_NB2 =
        NA_real_,

      left_gap =
        NA_real_,

      right_gap =
        NA_real_,

      diff_gap =
        NA_real_,

      left_alpha =
        NA_real_,

      right_alpha =
        NA_real_,

      diff_alpha =
        NA_real_,

      stringsAsFactors = FALSE
    )


    write.csv(

      fail_summary,

      cut_path,

      row.names = FALSE
    )


    message(

      "[",
      track,
      "] ",
      comparison_name,
      " ",
      arm_name,
      " | UNSTABLE: no candidate fraction passed."
    )


    return(

      list(

        summary =
          fail_summary,

        selected =
          NULL,

        stability_summary =
          stability_obj$stability_summary
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Stable solution
  # ---------------------------------------------------------------------------

  selected_stability <- stability_obj$selected

  selected_fraction <- selected_stability$min_fraction

  anchor <- selected_stability$anchor

  geometry <- stability_obj$original_geometry

  rank_order <- geometry$rank_order

  variance_df <- geometry$variance_df

  zero_df <- geometry$zero_df


  feature_df <- compute_ranked_feature_metrics(

    metric_mat_arm =
      metric_matrix,

    rank_order =
      rank_order
  )


  feature_df$abs_pc1_loading <- geometry$abs_loadings[
    rank_order
  ]


  region_summary <- summarize_regions(

    feature_df =
      feature_df,

    anchor_rank =
      anchor,

    total_n =
      total_n
  )


  # ---------------------------------------------------------------------------
  # Validation
  # ---------------------------------------------------------------------------

  validate_method_level(

    feature_df =
      feature_df,

    zero_df =
      zero_df,

    anchor =
      anchor,

    selected_fraction =
      selected_fraction,

    selected_stability_row =
      selected_stability,

    region_summary =
      region_summary,

    total_n =
      total_n
  )


  # ---------------------------------------------------------------------------
  # Selected cutoff summary
  # ---------------------------------------------------------------------------

  selected_df <- data.frame(

    comp =
      comparison_name,

    arm =
      arm_name,

    track =
      track,

    rank_method =
      rank_tag,

    metric_matrix =
      metric_name,

    Status =
      "PASS",

    SelectedMinFraction =
      selected_fraction,

    Anchor =
      anchor,

    ActualRightFraction =
      selected_stability$actual_right_fraction,

    BootstrapValidRate =
      selected_stability$bootstrap_valid_rate,

    AnchorRecoveryRate =
      selected_stability$anchor_recovery_rate,

    AnchorIQR =
      selected_stability$anchor_iqr,

    AnchorIQRFraction =
      selected_stability$anchor_iqr_fraction,

    MedianJaccard =
      selected_stability$median_jaccard,

    stringsAsFactors = FALSE
  )


  cutoff_summary <- bind_cols(

    selected_df,

    region_summary
  ) %>%

    mutate(

      LeftSize =
        region_summary$left_n,

      RightSize =
        region_summary$right_n,

      Call_NB2 =
        ifelse(
          diff_NB2 > 0,
          "RIGHT",
          "NOT_RIGHT"
        ),

      Call_Gap =
        ifelse(
          diff_gap > 0,
          "RIGHT",
          "NOT_RIGHT"
        ),

      Call_Alpha =
        ifelse(
          diff_alpha > 0,
          "RIGHT",
          "NOT_RIGHT"
        )
    )


  # ---------------------------------------------------------------------------
  # Zero-crossing table
  # ---------------------------------------------------------------------------

  zero_out <- zero_df %>%

    mutate(

      crossing_rank_rounded =
        as.integer(
          round(
            crossing_rank
          )
        ),

      selected_anchor =
        crossing_rank_rounded ==
        anchor
    )


  write.csv(

    zero_out,

    zero_path,

    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Validation table
  # ---------------------------------------------------------------------------

  write.csv(

    data.frame(

      comp =
        comparison_name,

      arm =
        arm_name,

      track =
        track,

      Status =
        "PASS",

      SelectedMinFraction =
        selected_fraction,

      Anchor =
        anchor,

      ActualRightFraction =
        selected_stability$actual_right_fraction,

      FractionRequirementMet =
        selected_stability$actual_right_fraction >=
        selected_fraction,

      BootstrapValidRate =
        selected_stability$bootstrap_valid_rate,

      AnchorRecoveryRate =
        selected_stability$anchor_recovery_rate,

      AnchorIQRFraction =
        selected_stability$anchor_iqr_fraction,

      MedianJaccard =
        selected_stability$median_jaccard,

      LeftN =
        region_summary$left_n,

      RightN =
        region_summary$right_n,

      EqualRegionSize =
        region_summary$left_n ==
        region_summary$right_n,

      AnchorIsZeroCrossing =
        anchor %in%
        as.integer(
          round(
            zero_df$crossing_rank
          )
        ),

      stringsAsFactors = FALSE
    ),

    valid_path,

    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Rank-level table
  # ---------------------------------------------------------------------------

  write.csv(

    feature_df %>%
      left_join(
        variance_df,
        by = "rank"
      ),

    rank_path,

    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Cutoff table
  # ---------------------------------------------------------------------------

  write.csv(

    cutoff_summary,

    cut_path,

    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Figure
  # ---------------------------------------------------------------------------

  build_main_figure(

    comparison_name =
      comparison_name,

    arm_name =
      arm_name,

    variance_df =
      variance_df,

    zero_df =
      zero_df,

    feature_df =
      feature_df,

    anchor =
      anchor,

    selected_stability_row =
      selected_stability,

    region_summary =
      region_summary,

    out_file =
      fig_path,

    fig_tag =
      fig_tag,

    rank_tag =
      rank_tag,

    metric_tag =
      metric_tag
  )


  message(

    "[",
    track,
    "] ",
    comparison_name,
    " ",
    arm_name,

    " | min_fraction=",
    round(
      selected_fraction,
      4
    ),

    " | Anchor=",
    anchor,

    " | RIGHT n=",
    region_summary$right_n,

    " | recovery=",
    round(
      selected_stability$anchor_recovery_rate,
      3
    ),

    " | Jaccard=",
    round(
      selected_stability$median_jaccard,
      3
    )
  )


  list(

    summary =
      cutoff_summary,

    selected =
      selected_df,

    stability_summary =
      stability_obj$stability_summary
  )
}


# =============================================================================
# RUN
# =============================================================================

count_mat <- read_count_matrix(
  COUNT_FILE
)


message(
  "Using count file: ",
  COUNT_FILE
)


message(
  "Count matrix dimensions: ",
  nrow(count_mat),
  " features x ",
  ncol(count_mat),
  " samples"
)


message(
  "Candidate minimum fractions: ",
  paste(
    MIN_FRACTION_GRID,
    collapse = ", "
  )
)


message(
  "Bootstrap replicates per track: ",
  BOOTSTRAP_N
)


overall_rows <- list()

overall_stability_rows <- list()


for (
  comparison_name in names(
    COMPARISONS
  )
) {

  comp_dir <- file.path(
    OUT_ROOT,
    comparison_name
  )

  dir.create(

    comp_dir,

    recursive = TRUE,

    showWarnings = FALSE
  )


  pats <- COMPARISONS[
    [
      comparison_name
    ]
  ]


  for (
    arm_name in c(
      "control",
      "treatment"
    )
  ) {

    sample_idx <- grep(

      pats[[arm_name]],

      colnames(count_mat)
    )


    if (
      length(sample_idx) < 2L
    ) {

      stop(
        "Not enough samples for ",
        comparison_name,
        " ",
        arm_name
      )
    }


    count_mat_arm <- count_mat[
      ,
      sample_idx,
      drop = FALSE
    ]


    # =========================================================================
    # MAIN TRACK
    # =========================================================================

    main_rank_matrix <- normalize_cpm_log1p(
      count_mat_arm
    )


    main_seed <- seed_from_label(

      paste(
        comparison_name,
        arm_name,
        "Main",
        sep = "_"
      )
    )


    main_res <- run_one_track(

      comparison_name =
        comparison_name,

      arm_name =
        arm_name,

      rank_matrix =
        main_rank_matrix,

      metric_matrix =
        count_mat_arm,

      output_dir =
        comp_dir,

      track =
        "Main",

      fig_tag =
        "Main",

      rank_tag =
        "Rank: CPM log1p",

      metric_tag =
        "Metrics: raw counts",

      metric_name =
        "raw_counts",

      seed =
        main_seed
    )


    overall_rows[
      [
        length(
          overall_rows
        ) +
          1L
      ]
    ] <- main_res$summary


    overall_stability_rows[
      [
        length(
          overall_stability_rows
        ) +
          1L
      ]
    ] <- main_res$stability_summary %>%

      mutate(

        comp =
          comparison_name,

        arm =
          arm_name,

        track =
          "Main",

        .before =
          1
      )


    # =========================================================================
    # DESEQ2 SUPPLEMENT
    # =========================================================================

    if (
      RUN_DESEQ2_SUPPLEMENT
    ) {

      deseq2_obj <- compute_deseq2_matrices(

        count_mat_arm,

        rank_method =
          DESEQ2_RANK_METHOD
      )


      if (
        is.null(
          deseq2_obj
        )
      ) {

        message(

          "[DESeq2] skipped: package not available for ",
          comparison_name,
          " ",
          arm_name
        )

      } else {


        deseq2_rank_tag <- switch(

          deseq2_obj$rank_method_used,

          vst =
            "Rank: DESeq2 VST",

          normalized_log1p =
            "Rank: DESeq2 log1p",

          normalized_log1p_fallback =
            "Rank: DESeq2 log1p (VST fallback)",

          "Rank: DESeq2 log1p"
        )


        deseq2_seed <- seed_from_label(

          paste(
            comparison_name,
            arm_name,
            "DESeq2",
            sep = "_"
          )
        )


        deseq2_res <- run_one_track(

          comparison_name =
            comparison_name,

          arm_name =
            arm_name,

          rank_matrix =
            deseq2_obj$ranking_matrix,

          metric_matrix =
            deseq2_obj$normalized_counts,

          output_dir =
            comp_dir,

          track =
            "DESeq2",

          fig_tag =
            "DESeq2 supplement",

          rank_tag =
            deseq2_rank_tag,

          metric_tag =
            "Metrics: DESeq2 normalized",

          metric_name =
            "deseq2_normalized_counts",

          seed =
            deseq2_seed
        )


        write.csv(

          data.frame(

            sample =
              names(
                deseq2_obj$size_factors
              ),

            size_factor =
              as.numeric(
                deseq2_obj$size_factors
              ),

            stringsAsFactors = FALSE
          ),

          file.path(

            comp_dir,

            paste0(
              "Table_SizeFactor_",
              comparison_name,
              "_",
              arm_name,
              ".csv"
            )
          ),

          row.names = FALSE
        )


        overall_rows[
          [
            length(
              overall_rows
            ) +
              1L
          ]
        ] <- deseq2_res$summary


        overall_stability_rows[
          [
            length(
              overall_stability_rows
            ) +
              1L
          ]
        ] <- deseq2_res$stability_summary %>%

          mutate(

            comp =
              comparison_name,

            arm =
              arm_name,

            track =
              "DESeq2",

            .before =
              1
          )
      }
    }
  }
}


# =============================================================================
# OVERALL OUTPUTS
# =============================================================================

overall_summary <- bind_rows(
  overall_rows
)


write.csv(

  overall_summary,

  file.path(
    OUT_ROOT,
    "Table_Overall_Cutoff.csv"
  ),

  row.names = FALSE
)


overall_stability <- bind_rows(
  overall_stability_rows
)


write.csv(

  overall_stability,

  file.path(
    OUT_ROOT,
    "Table_Overall_Stability.csv"
  ),

  row.names = FALSE
)


message(
  "Done. Outputs written to: ",
  OUT_ROOT
)
