#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# MANUSCRIPT ANALYSIS
# BOOTSTRAP-STABLE, DATA-DERIVED HIGH-PC1 CUTOFF
# =============================================================================
#
# METHODS
#
# Within each experimental arm, features are ordered from lowest to highest
# absolute loading on the first principal component (|PC1 loading|). For the
# primary analysis, PC1 is calculated from log1p-transformed counts normalized
# to counts per million (CPM), whereas empirical feature variance is calculated
# from the corresponding raw count matrix. For the supplementary DESeq2
# analysis, PC1 ranking and variance estimation are based on DESeq2-normalized
# quantities, with ranking performed on either log1p-normalized counts or a
# variance-stabilizing transformation (VST), as specified below.
#
# To identify a data-derived boundary in the high-|PC1|-loading tail, empirical
# feature variance is ordered along the PC1 rank axis, transformed as
# log(1 + variance), and represented by a smoothing spline. The second
# derivative of the spline is evaluated over a dense rank grid. Sign changes
# in the second derivative define candidate curvature transitions
# (inflection points).
#
# No fixed feature-count cutoff (e.g., 5,000 features) is used.
#
# Instead, a prespecified grid of minimum high-loading tail fractions is
# evaluated:
#
#     5%, 7.5%, 10%, 12.5%, and 15%.
#
# For a candidate minimum fraction f, an eligible Anchor must:
#
#   1. correspond to a spline second-derivative sign-change crossing;
#   2. leave at least f*N features from the Anchor through the right edge;
#   3. permit construction of an immediately adjacent, equally sized LEFT
#      comparison block.
#
# Among all eligible crossings for a given f, the rightmost crossing is chosen.
# Thus, the procedure identifies the most terminal curvature transition that
# still supports the prespecified minimum tail size and matched comparison.
#
# Cutoff reproducibility is then evaluated by bootstrap resampling of biological
# samples within each arm. For every bootstrap replicate, the complete
# track-specific preprocessing, PC1 ranking, empirical variance trajectory,
# spline fit, second derivative, and Anchor selection are recomputed.
#
# A candidate fraction is classified as geometrically stable only if all of
# the following prespecified criteria are satisfied:
#
#   - >= 90% of bootstrap replicates yield a valid Anchor;
#   - >= 80% of valid replicates recover an Anchor within 3% of the total
#     rank axis from the full-data Anchor;
#   - bootstrap Anchor IQR is <= 3% of the total rank axis;
#   - bootstrap RIGHT-region fraction IQR is <= 3 percentage points;
#   - median Jaccard overlap between bootstrap and full-data RIGHT feature
#     membership is >= 0.80.
#
# The smallest candidate minimum fraction satisfying all stability criteria is
# selected. The final Anchor is the corresponding Anchor obtained from the full
# dataset; the bootstrap median Anchor is not substituted for the full-data
# estimate.
#
# After cutoff selection is complete, the RIGHT region is defined as all
# features from the Anchor through the highest-|PC1|-loading rank. LEFT is the
# immediately preceding block containing exactly the same number of features.
#
# NB2-related quantities are evaluated only after the geometric cutoff has been
# selected:
#
#     NB2 = log(1 + max(variance - mu, 0))
#
#     NB2-NB1 =
#       log(1 + max(variance - mu, 0)) - log(1 + mu)
#
#     alpha = max((variance - mu) / mu^2, 0)
#
#     alpha*mu signal = log(1 + alpha*mu)
#
# where mu and variance are the empirical feature mean and variance.
#
# These quantities are descriptive measures of excess-variance behavior and
# are not formal likelihood-ratio statistics or independent validation tests.
# Importantly, none of the NB2-related quantities is used to select the minimum
# fraction, Anchor, spline smoothing parameter, or bootstrap stability result.
#
# If no candidate fraction satisfies all prespecified stability criteria, the
# analysis track is reported as UNSTABLE and no cutoff or corroborative figure
# is forced.
#
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <- "/root/REAPER98632/exports/manuscript_bootstrap_stable"


# -----------------------------------------------------------------------------
# Spline geometry
# -----------------------------------------------------------------------------

VAR_SPLINE_SPAR <- 0.60

DENSE_GRID_MULTIPLIER <- 2L
DENSE_GRID_MIN <- 5000L


# -----------------------------------------------------------------------------
# Prespecified candidate minimum RIGHT-tail fractions
# -----------------------------------------------------------------------------

MIN_FRACTION_GRID <- c(
  0.050,
  0.075,
  0.100,
  0.125,
  0.150
)


# -----------------------------------------------------------------------------
# Bootstrap
# -----------------------------------------------------------------------------

BOOTSTRAP_N <- 500L
BOOTSTRAP_SEED_BASE <- 20260820L

detected_cores <- suppressWarnings(
  parallel::detectCores(logical = FALSE)
)

if (is.na(detected_cores) || detected_cores < 2L) {
  BOOTSTRAP_CORES <- 1L
} else {
  BOOTSTRAP_CORES <- max(
    1L,
    min(
      4L,
      detected_cores - 1L
    )
  )
}


# -----------------------------------------------------------------------------
# Prespecified stability criteria
# -----------------------------------------------------------------------------

MIN_BOOTSTRAP_VALID_RATE <- 0.90

ANCHOR_TOLERANCE_FRACTION <- 0.03
MIN_ANCHOR_RECOVERY_RATE <- 0.80

MAX_ANCHOR_IQR_FRACTION <- 0.03

MAX_RIGHT_FRACTION_IQR <- 0.03

MIN_MEDIAN_JACCARD <- 0.80


# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

PNG_WIDTH_IN <- 14
PNG_HEIGHT_IN <- 10.8
PNG_DPI <- 260


# -----------------------------------------------------------------------------
# DESeq2 supplementary analysis
# -----------------------------------------------------------------------------

RUN_DESEQ2_SUPPLEMENT <- TRUE

DESEQ2_RANK_METHOD <- "normalized_log1p"
# Allowed:
#   "normalized_log1p"
#   "vst"


# -----------------------------------------------------------------------------
# Experimental comparisons
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

  anchor        = "#000000",
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
# INPUT
# =============================================================================

read_count_matrix <- function(path, comparisons) {

  if (!file.exists(path)) {
    stop("Count file does not exist: ", path)
  }

  raw_df <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(raw_df) == 0L || ncol(raw_df) < 2L) {
    stop("Count file is empty or malformed: ", path)
  }

  # Identify sample columns from the experimental sample-name patterns.
  # This avoids coercing unrelated metadata columns to numeric.

  sample_patterns <- unique(
    unname(
      unlist(
        comparisons,
        recursive = TRUE,
        use.names = FALSE
      )
    )
  )

  sample_idx <- sort(
    unique(
      unlist(
        lapply(
          sample_patterns,
          function(pattern) {
            grep(
              pattern,
              colnames(raw_df)
            )
          }
        )
      )
    )
  )

  if (length(sample_idx) == 0L) {
    stop(
      "No sample columns matched the patterns specified in COMPARISONS."
    )
  }

  if (1L %in% sample_idx) {
    stop(
      "The first column matched a sample pattern. ",
      "The first column is expected to contain feature identifiers."
    )
  }

  feature_ids <- trimws(
    as.character(
      raw_df[[1L]]
    )
  )

  if (
    any(is.na(feature_ids)) ||
    any(feature_ids == "")
  ) {
    stop(
      "Missing or empty feature identifiers were found in column 1."
    )
  }

  feature_ids <- make.unique(feature_ids)

  parsed_columns <- vector(
    mode = "list",
    length = length(sample_idx)
  )

  for (j in seq_along(sample_idx)) {

    col_index <- sample_idx[j]
    col_name <- colnames(raw_df)[col_index]

    original <- raw_df[[col_index]]

    if (is.numeric(original) || is.integer(original)) {

      numeric_values <- as.numeric(original)

    } else {

      numeric_values <- suppressWarnings(
        as.numeric(
          trimws(
            as.character(original)
          )
        )
      )
    }

    if (any(!is.finite(numeric_values))) {

      bad_rows <- which(
        !is.finite(numeric_values)
      )

      preview_rows <- paste(
        head(
          bad_rows,
          10L
        ),
        collapse = ", "
      )

      stop(
        "Sample column '",
        col_name,
        "' contains missing or non-numeric count values. ",
        "Example row indices: ",
        preview_rows,
        ". Counts are not silently replaced by zero."
      )
    }

    if (any(numeric_values < 0)) {
      stop(
        "Negative values detected in sample column: ",
        col_name
      )
    }

    parsed_columns[[j]] <- numeric_values
  }

  count_mat <- do.call(
    cbind,
    parsed_columns
  )

  rownames(count_mat) <- feature_ids

  colnames(count_mat) <- colnames(raw_df)[sample_idx]

  storage.mode(count_mat) <- "double"

  if (any(!is.finite(count_mat))) {
    stop("Non-finite count values remain after parsing.")
  }

  if (any(count_mat < 0)) {
    stop("Negative count values remain after parsing.")
  }

  # This is a raw-read-count analysis. Noninteger values indicate that
  # the wrong columns or an inappropriate input matrix may have been supplied.

  noninteger <- abs(
    count_mat -
      round(count_mat)
  ) > 1e-8

  if (any(noninteger)) {

    bad_col <- which(
      colSums(noninteger) > 0
    )[1L]

    stop(
      "Noninteger values were detected in raw count column '",
      colnames(count_mat)[bad_col],
      "'. Verify that COUNT_FILE contains raw read counts."
    )
  }

  count_mat <- round(count_mat)

  # Remove features with zero counts across every selected sample.

  keep <- rowSums(count_mat) > 0

  if (!any(keep)) {
    stop("All features have zero total counts.")
  }

  count_mat[
    keep,
    ,
    drop = FALSE
  ]
}


# =============================================================================
# NORMALIZATION
# =============================================================================

normalize_cpm_log1p <- function(count_mat_arm) {

  lib_sizes <- colSums(
    count_mat_arm
  )

  if (
    any(!is.finite(lib_sizes)) ||
    any(lib_sizes <= 0)
  ) {
    stop(
      "At least one sample has a nonpositive or invalid library size."
    )
  }

  cpm <- sweep(
    count_mat_arm,
    2L,
    lib_sizes / 1e6,
    "/"
  )

  log1p(cpm)
}


prepare_main_matrices <- function(count_mat_arm) {

  list(
    rank_matrix = normalize_cpm_log1p(
      count_mat_arm
    ),

    metric_matrix = count_mat_arm,

    rank_method_used = "CPM log1p",

    metric_method_used = "raw counts",

    size_factors = NULL
  )
}


prepare_deseq2_matrices <- function(
    count_mat_arm,
    rank_method = DESEQ2_RANK_METHOD) {

  if (
    !requireNamespace(
      "DESeq2",
      quietly = TRUE
    )
  ) {
    stop("DESeq2 is not installed.")
  }

  if (
    !requireNamespace(
      "SummarizedExperiment",
      quietly = TRUE
    )
  ) {
    stop("SummarizedExperiment is not installed.")
  }

  if (
    !(rank_method %in% c("normalized_log1p", "vst"))
  ) {
    stop(
      "DESEQ2_RANK_METHOD must be 'normalized_log1p' or 'vst'."
    )
  }

  if (
    any(
      abs(
        count_mat_arm -
          round(count_mat_arm)
      ) > 1e-8
    )
  ) {
    stop(
      "DESeq2 requires integer raw counts."
    )
  }

  # Bootstrap resampling can repeat the same biological sample.
  # DESeq2 requires unique column names, so duplicate draws are assigned
  # unique computational identifiers without changing their values.

  if (
    is.null(
      colnames(count_mat_arm)
    )
  ) {
    colnames(count_mat_arm) <- paste0(
      "sample_",
      seq_len(
        ncol(count_mat_arm)
      )
    )
  }

  colnames(count_mat_arm) <- make.unique(
    colnames(count_mat_arm)
  )

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

  if (rank_method == "normalized_log1p") {

    ranking_matrix <- log1p(
      norm_counts
    )

    rank_method_used <- "DESeq2 normalized log1p"

  } else {

    vst_obj <- DESeq2::vst(
      dds,
      blind = TRUE
    )

    ranking_matrix <- SummarizedExperiment::assay(
      vst_obj
    )

    rank_method_used <- "DESeq2 VST"
  }

  list(
    rank_matrix = ranking_matrix,

    metric_matrix = norm_counts,

    rank_method_used = rank_method_used,

    metric_method_used = "DESeq2 normalized counts",

    size_factors = DESeq2::sizeFactors(dds)
  )
}


# =============================================================================
# NUMERICAL HELPERS
# =============================================================================

fast_row_variance <- function(mat) {

  n <- ncol(mat)

  if (n < 2L) {
    stop(
      "At least two samples are required to calculate feature variance."
    )
  }

  mu <- rowMeans(mat)

  ss <- rowSums(
    mat * mat
  ) -
    n * mu^2

  # Floating-point cancellation can produce tiny negative values.

  tolerance <- .Machine$double.eps *
    pmax(
      1,
      abs(
        rowSums(
          mat * mat
        )
      )
    )

  ss[
    ss < 0 &
      abs(ss) <= tolerance
  ] <- 0

  variance <- ss /
    (n - 1L)

  variance[
    !is.finite(variance)
  ] <- NA_real_

  variance <- pmax(
    variance,
    0,
    na.rm = FALSE
  )

  variance
}


# =============================================================================
# PC1 RANKING
# =============================================================================

compute_abs_pc1_loadings <- function(norm_mat_arm) {

  if (
    nrow(norm_mat_arm) < 2L ||
    ncol(norm_mat_arm) < 2L
  ) {
    stop(
      "PC1 calculation requires at least two features and two samples."
    )
  }

  # X:
  #   rows    = samples
  #   columns = features

  X <- t(
    norm_mat_arm
  )

  feature_means <- colMeans(X)

  X_centered <- sweep(
    X,
    2L,
    feature_means,
    "-"
  )

  # Because the number of samples is much smaller than the number of
  # features, PC1 is obtained through the sample-space Gram matrix.
  # This is mathematically equivalent to the leading right singular
  # vector used by prcomp() but substantially faster for repeated
  # bootstrap analyses.

  gram <- tcrossprod(
    X_centered
  )

  eig <- eigen(
    gram,
    symmetric = TRUE
  )

  lambda1 <- eig$values[1L]

  scale_reference <- max(
    1,
    max(
      abs(
        eig$values
      ),
      na.rm = TRUE
    )
  )

  if (
    !is.finite(lambda1) ||
    lambda1 <=
      .Machine$double.eps *
      scale_reference
  ) {
    stop(
      "PC1 is undefined because the bootstrap sample matrix has ",
      "insufficient between-sample variation."
    )
  }

  u1 <- eig$vectors[
    ,
    1L
  ]

  singular_value <- sqrt(
    lambda1
  )

  loading <- as.numeric(
    crossprod(
      X_centered,
      u1
    )
  ) /
    singular_value

  names(loading) <- colnames(
    X_centered
  )

  loading[
    !is.finite(loading)
  ] <- 0

  abs(loading)
}


# =============================================================================
# VARIANCE GEOMETRY
# =============================================================================

compute_ranked_variance_curve <- function(
    metric_mat_arm,
    rank_order,
    spar = VAR_SPLINE_SPAR) {

  empirical_var <- fast_row_variance(
    metric_mat_arm
  )

  if (
    any(
      !is.finite(
        empirical_var
      )
    )
  ) {
    stop(
      "Non-finite empirical feature variances were produced."
    )
  }

  ranked_var <- empirical_var[
    rank_order
  ]

  ranked_log_var <- log1p(
    ranked_var
  )

  ranks <- seq_along(
    rank_order
  )

  spline_fit <- tryCatch(
    stats::smooth.spline(
      x = ranks,
      y = ranked_log_var,
      spar = spar
    ),
    error = function(e) {
      NULL
    }
  )

  if (is.null(spline_fit)) {
    stop(
      "Smoothing spline fit failed."
    )
  }

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
    DENSE_GRID_MIN,
    length(ranks) *
      DENSE_GRID_MULTIPLIER
  )

  dense_x <- seq(
    from = min(ranks),
    to = max(ranks),
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

    empirical_variance = ranked_var,

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
    dense_rank = dense_x,

    dense_smooth_log1p_empirical_variance =
      dense_y,

    dense_d2 = dense_d2,

    stringsAsFactors = FALSE
  )

  out
}


find_d2_zero_crossings <- function(dense_df) {

  x <- dense_df$dense_rank
  y <- dense_df$dense_d2

  ok <- is.finite(x) &
    is.finite(y)

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

    if (
      (a < 0 && b > 0) ||
      (a > 0 && b < 0)
    ) {

      fraction_between_points <- abs(a) /
        (
          abs(a) +
            abs(b)
        )

      crossing_rank <- xa +
        fraction_between_points *
        (
          xb - xa
        )

      crossings <- c(
        crossings,
        crossing_rank
      )
    }
  }

  crossings <- crossings[
    is.finite(crossings)
  ]

  if (length(crossings) == 0L) {

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
        digits = 8L
      )
    )
  )

  data.frame(
    crossing_rank = crossings,
    crossing_type = "d2_sign_change",
    stringsAsFactors = FALSE
  )
}


compute_geometry <- function(
    rank_matrix,
    metric_matrix,
    spar = VAR_SPLINE_SPAR) {

  if (
    nrow(rank_matrix) !=
      nrow(metric_matrix)
  ) {
    stop(
      "rank_matrix and metric_matrix contain different feature counts."
    )
  }

  if (
    ncol(rank_matrix) !=
      ncol(metric_matrix)
  ) {
    stop(
      "rank_matrix and metric_matrix contain different sample counts."
    )
  }

  if (
    !identical(
      rownames(rank_matrix),
      rownames(metric_matrix)
    )
  ) {
    stop(
      "Feature identifiers differ between rank_matrix and metric_matrix."
    )
  }

  abs_loadings <- compute_abs_pc1_loadings(
    rank_matrix
  )

  rank_order <- order(
    abs_loadings,
    decreasing = FALSE
  )

  variance_df <- compute_ranked_variance_curve(
    metric_mat_arm = metric_matrix,
    rank_order = rank_order,
    spar = spar
  )

  dense_curve_df <- attr(
    variance_df,
    "dense_curve_df"
  )

  zero_df <- find_d2_zero_crossings(
    dense_curve_df
  )

  list(
    abs_loadings = abs_loadings,
    rank_order = rank_order,
    variance_df = variance_df,
    dense_curve_df = dense_curve_df,
    zero_df = zero_df
  )
}


# =============================================================================
# CANDIDATE ANCHOR SELECTION
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
      "min_fraction must lie strictly between 0 and 0.50."
    )
  }

  minimum_right_n <- ceiling(
    min_fraction *
      total_n
  )

  # Equal-size LEFT comparison requires:
  #
  #   Anchor - 1 >= N - Anchor + 1
  #
  # therefore:
  #
  #   Anchor >= (N + 2) / 2

  minimum_anchor_for_match <- ceiling(
    (
      total_n +
        2L
    ) /
      2
  )

  # Minimum RIGHT size requires:
  #
  #   N - Anchor + 1 >= minimum_right_n
  #
  # therefore:
  #
  #   Anchor <= N - minimum_right_n + 1

  maximum_anchor_for_tail <- total_n -
    minimum_right_n +
    1L

  invalid_result <- function() {

    data.frame(
      min_fraction = min_fraction,
      valid = FALSE,
      anchor = NA_integer_,
      right_n = NA_integer_,
      left_n = NA_integer_,
      actual_right_fraction = NA_real_,
      minimum_right_n = minimum_right_n,
      minimum_anchor_for_match = minimum_anchor_for_match,
      maximum_anchor_for_tail = maximum_anchor_for_tail,
      stringsAsFactors = FALSE
    )
  }

  if (
    maximum_anchor_for_tail <
      minimum_anchor_for_match
  ) {
    return(
      invalid_result()
    )
  }

  if (nrow(zero_df) == 0L) {
    return(
      invalid_result()
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

  if (length(eligible) == 0L) {
    return(
      invalid_result()
    )
  }

  # Primary geometric rule:
  #
  # Choose the most terminal/rightmost curvature transition that still
  # leaves the required minimum tail size and allows a matched LEFT block.

  anchor <- max(
    eligible
  )

  right_n <- total_n -
    anchor +
    1L

  left_n <- right_n

  actual_right_fraction <- right_n /
    total_n

  data.frame(
    min_fraction = min_fraction,
    valid = TRUE,
    anchor = anchor,
    right_n = right_n,
    left_n = left_n,
    actual_right_fraction = actual_right_fraction,
    minimum_right_n = minimum_right_n,
    minimum_anchor_for_match = minimum_anchor_for_match,
    maximum_anchor_for_tail = maximum_anchor_for_tail,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# FEATURE-SET STABILITY
# =============================================================================

jaccard_similarity <- function(set_a, set_b) {

  set_a <- unique(set_a)
  set_b <- unique(set_b)

  union_set <- union(
    set_a,
    set_b
  )

  if (length(union_set) == 0L) {
    return(
      NA_real_
    )
  }

  length(
    intersect(
      set_a,
      set_b
    )
  ) /
    length(
      union_set
    )
}


seed_from_label <- function(
    label,
    base_seed = BOOTSTRAP_SEED_BASE) {

  label_value <- sum(
    utf8ToInt(
      as.character(label)
    )
  )

  seed <- (
    base_seed +
      label_value *
      1009L
  ) %%
    2147483647L

  as.integer(seed)
}


# =============================================================================
# BOOTSTRAP GEOMETRIC STABILITY
# =============================================================================

assess_geometry_stability <- function(
    count_mat_arm,
    prepare_function,
    candidate_fractions = MIN_FRACTION_GRID,
    bootstrap_n = BOOTSTRAP_N,
    spar = VAR_SPLINE_SPAR,
    seed = BOOTSTRAP_SEED_BASE,
    cores = BOOTSTRAP_CORES) {

  candidate_fractions <- sort(
    unique(
      as.numeric(
        candidate_fractions
      )
    )
  )

  if (
    any(
      !is.finite(
        candidate_fractions
      )
    ) ||
    any(
      candidate_fractions <= 0
    ) ||
    any(
      candidate_fractions >= 0.50
    )
  ) {
    stop(
      "All candidate fractions must lie strictly between 0 and 0.50."
    )
  }

  if (bootstrap_n < 1L) {
    stop(
      "bootstrap_n must be >= 1."
    )
  }

  total_n <- nrow(
    count_mat_arm
  )

  sample_n <- ncol(
    count_mat_arm
  )

  if (sample_n < 2L) {
    stop(
      "At least two biological samples are required."
    )
  }


  # ---------------------------------------------------------------------------
  # Full-data preprocessing and geometry
  # ---------------------------------------------------------------------------

  original_prepared <- prepare_function(
    count_mat_arm
  )

  original_geometry <- compute_geometry(
    rank_matrix = original_prepared$rank_matrix,
    metric_matrix = original_prepared$metric_matrix,
    spar = spar
  )

  original_candidates <- bind_rows(
    lapply(
      candidate_fractions,
      function(f) {

        select_anchor_for_min_fraction(
          zero_df = original_geometry$zero_df,
          total_n = total_n,
          min_fraction = f
        )
      }
    )
  )


  # ---------------------------------------------------------------------------
  # Full-data RIGHT feature sets
  # ---------------------------------------------------------------------------

  original_right_sets <- vector(
    mode = "list",
    length = length(
      candidate_fractions
    )
  )

  names(original_right_sets) <- sprintf(
    "%.6f",
    candidate_fractions
  )

  for (
    j in seq_along(
      candidate_fractions
    )
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
        count_mat_arm
      )[
        right_order_idx
      ]

    } else {

      original_right_sets[[key]] <- character(0)
    }
  }


  # ---------------------------------------------------------------------------
  # Generate all bootstrap samples before parallel execution.
  #
  # This makes the bootstrap sample draws deterministic regardless of the
  # number of worker processes.
  # ---------------------------------------------------------------------------

  set.seed(seed)

  bootstrap_indices <- lapply(
    seq_len(
      bootstrap_n
    ),
    function(b) {

      sample.int(
        n = sample_n,
        size = sample_n,
        replace = TRUE
      )
    }
  )


  # ---------------------------------------------------------------------------
  # One bootstrap replicate
  # ---------------------------------------------------------------------------

  bootstrap_worker <- function(b) {

    boot_sample_idx <- bootstrap_indices[[b]]

    make_invalid_rows <- function(reason) {

      bind_rows(
        lapply(
          candidate_fractions,
          function(f) {

            data.frame(
              bootstrap = b,
              min_fraction = f,
              valid = FALSE,
              anchor = NA_integer_,
              right_n = NA_integer_,
              right_fraction = NA_real_,
              jaccard = NA_real_,
              invalid_reason = reason,
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }

    # A resample containing only one unique biological sample cannot
    # provide meaningful between-sample variance or PCA geometry.

    if (
      length(
        unique(
          boot_sample_idx
        )
      ) < 2L
    ) {
      return(
        make_invalid_rows(
          "fewer_than_2_unique_samples"
        )
      )
    }

    boot_counts <- count_mat_arm[
      ,
      boot_sample_idx,
      drop = FALSE
    ]

    # Repeated bootstrap draws are intentionally preserved, but computational
    # column identifiers are made unique for packages such as DESeq2.

    colnames(boot_counts) <- paste0(
      colnames(
        count_mat_arm
      )[
        boot_sample_idx
      ],
      "__bootstrap_",
      b,
      "_draw_",
      seq_along(
        boot_sample_idx
      )
    )

    boot_prepared <- tryCatch(
      prepare_function(
        boot_counts
      ),
      error = function(e) {
        NULL
      }
    )

    if (is.null(boot_prepared)) {
      return(
        make_invalid_rows(
          "preprocessing_failed"
        )
      )
    }

    boot_geometry <- tryCatch(
      compute_geometry(
        rank_matrix = boot_prepared$rank_matrix,
        metric_matrix = boot_prepared$metric_matrix,
        spar = spar
      ),
      error = function(e) {
        NULL
      }
    )

    if (is.null(boot_geometry)) {
      return(
        make_invalid_rows(
          "geometry_failed"
        )
      )
    }

    result_rows <- vector(
      "list",
      length(
        candidate_fractions
      )
    )

    for (
      j in seq_along(
        candidate_fractions
      )
    ) {

      f <- candidate_fractions[j]

      selection <- select_anchor_for_min_fraction(
        zero_df = boot_geometry$zero_df,
        total_n = total_n,
        min_fraction = f
      )

      if (
        !isTRUE(
          selection$valid
        )
      ) {

        result_rows[[j]] <- data.frame(
          bootstrap = b,
          min_fraction = f,
          valid = FALSE,
          anchor = NA_integer_,
          right_n = NA_integer_,
          right_fraction = NA_real_,
          jaccard = NA_real_,
          invalid_reason = "no_valid_anchor",
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
        count_mat_arm
      )[
        boot_right_order_idx
      ]

      key <- sprintf(
        "%.6f",
        f
      )

      original_right_features <- original_right_sets[[key]]

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

      result_rows[[j]] <- data.frame(
        bootstrap = b,
        min_fraction = f,
        valid = TRUE,
        anchor = boot_anchor,
        right_n = selection$right_n,
        right_fraction = selection$actual_right_fraction,
        jaccard = jac,
        invalid_reason = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    bind_rows(
      result_rows
    )
  }


  # ---------------------------------------------------------------------------
  # Run bootstrap
  # ---------------------------------------------------------------------------

  bootstrap_ids <- seq_len(
    bootstrap_n
  )

  if (
    .Platform$OS.type == "unix" &&
    cores > 1L
  ) {

    bootstrap_list <- parallel::mclapply(
      bootstrap_ids,
      bootstrap_worker,
      mc.cores = cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )

  } else {

    bootstrap_list <- lapply(
      bootstrap_ids,
      bootstrap_worker
    )
  }

  bootstrap_df <- bind_rows(
    bootstrap_list
  )


  # ---------------------------------------------------------------------------
  # Summarize stability for each candidate minimum fraction
  # ---------------------------------------------------------------------------

  stability_rows <- vector(
    "list",
    length(
      candidate_fractions
    )
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

    valid_sub <- boot_sub[
      boot_sub$valid,
      ,
      drop = FALSE
    ]

    valid_bootstraps <- nrow(
      valid_sub
    )

    valid_rate <- valid_bootstraps /
      bootstrap_n


    if (valid_bootstraps > 0L) {

      anchor_median <- median(
        valid_sub$anchor,
        na.rm = TRUE
      )

      anchor_iqr <- stats::IQR(
        valid_sub$anchor,
        na.rm = TRUE
      )

      anchor_iqr_fraction <- anchor_iqr /
        total_n

      right_fraction_median <- median(
        valid_sub$right_fraction,
        na.rm = TRUE
      )

      right_fraction_iqr <- stats::IQR(
        valid_sub$right_fraction,
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
      ) &&
      valid_bootstraps > 0L
    ) {

      full_anchor <- original_row$anchor

      tolerance_n <- ceiling(
        ANCHOR_TOLERANCE_FRACTION *
          total_n
      )

      anchor_recovery_rate <- mean(
        abs(
          valid_sub$anchor -
            full_anchor
        ) <=
          tolerance_n
      )

    } else {

      full_anchor <- NA_integer_

      tolerance_n <- ceiling(
        ANCHOR_TOLERANCE_FRACTION *
          total_n
      )

      anchor_recovery_rate <- NA_real_
    }


    finite_jaccard <- valid_sub$jaccard[
      is.finite(
        valid_sub$jaccard
      )
    ]

    median_jaccard <- if (
      length(
        finite_jaccard
      ) > 0L
    ) {

      median(
        finite_jaccard
      )

    } else {

      NA_real_
    }


    pass_valid_rate <- is.finite(
      valid_rate
    ) &&
      valid_rate >=
      MIN_BOOTSTRAP_VALID_RATE


    pass_anchor_recovery <- is.finite(
      anchor_recovery_rate
    ) &&
      anchor_recovery_rate >=
      MIN_ANCHOR_RECOVERY_RATE


    pass_anchor_iqr <- is.finite(
      anchor_iqr_fraction
    ) &&
      anchor_iqr_fraction <=
      MAX_ANCHOR_IQR_FRACTION


    pass_right_fraction_iqr <- is.finite(
      right_fraction_iqr
    ) &&
      right_fraction_iqr <=
      MAX_RIGHT_FRACTION_IQR


    pass_jaccard <- is.finite(
      median_jaccard
    ) &&
      median_jaccard >=
      MIN_MEDIAN_JACCARD


    pass_stability <- isTRUE(
      original_row$valid
    ) &&
      pass_valid_rate &&
      pass_anchor_recovery &&
      pass_anchor_iqr &&
      pass_right_fraction_iqr &&
      pass_jaccard


    stability_rows[[j]] <- data.frame(
      min_fraction = f,

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
        valid_bootstraps,

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

      threshold_valid_rate =
        MIN_BOOTSTRAP_VALID_RATE,

      threshold_anchor_recovery =
        MIN_ANCHOR_RECOVERY_RATE,

      threshold_anchor_iqr_fraction =
        MAX_ANCHOR_IQR_FRACTION,

      threshold_right_fraction_iqr =
        MAX_RIGHT_FRACTION_IQR,

      threshold_median_jaccard =
        MIN_MEDIAN_JACCARD,

      pass_valid_rate =
        pass_valid_rate,

      pass_anchor_recovery =
        pass_anchor_recovery,

      pass_anchor_iqr =
        pass_anchor_iqr,

      pass_right_fraction_iqr =
        pass_right_fraction_iqr,

      pass_jaccard =
        pass_jaccard,

      pass_stability =
        pass_stability,

      stringsAsFactors = FALSE
    )
  }

  stability_summary <- bind_rows(
    stability_rows
  )


  # ---------------------------------------------------------------------------
  # Select the smallest candidate fraction satisfying all criteria
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
    original_prepared = original_prepared,
    original_geometry = original_geometry,
    original_candidates = original_candidates,
    bootstrap_df = bootstrap_df,
    stability_summary = stability_summary,
    selected = selected
  )
}


# =============================================================================
# NB2-RELATED FEATURE METRICS
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
    ranked_mat
  )

  empirical_var <- fast_row_variance(
    ranked_mat
  )

  if (
    any(!is.finite(mu)) ||
    any(!is.finite(empirical_var))
  ) {
    stop(
      "Non-finite feature means or variances were produced."
    )
  }

  mu <- pmax(
    mu,
    0
  )

  empirical_var <- pmax(
    empirical_var,
    0
  )

  extra_variance <- pmax(
    empirical_var -
      mu,
    0
  )

  alpha_hat <- numeric(
    length(mu)
  )

  positive_mu <- mu > 0

  alpha_hat[positive_mu] <- pmax(
    (
      empirical_var[positive_mu] -
        mu[positive_mu]
    ) /
      (
        mu[positive_mu]^2
      ),
    0
  )

  alpha_mu_value <- alpha_hat *
    mu

  data.frame(
    rank = seq_along(
      rank_order
    ),

    feature_id = rownames(
      metric_mat_arm
    )[
      rank_order
    ],

    mu = mu,

    empirical_variance =
      empirical_var,

    NB2 = log1p(
      extra_variance
    ),

    NB2_NB1 =
      log1p(
        extra_variance
      ) -
      log1p(mu),

    alpha_mu = log1p(
      alpha_mu_value
    ),

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# MATCHED LEFT/RIGHT REGIONS
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

  left_end <- anchor_rank -
    1L

  left_start <- left_end -
    right_n +
    1L

  if (left_start < 1L) {
    stop(
      "LEFT block extends below rank 1."
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
      "LEFT and RIGHT blocks are not equal in size."
    )
  }

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
    left_start = left_start,
    left_end = left_end,

    right_start = anchor_rank,
    right_end = total_n,

    left_n = nrow(left_df),
    right_n = nrow(right_df),

    left_NB2 = median(
      left_df$NB2,
      na.rm = TRUE
    ),

    right_NB2 = median(
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

    left_gap = median(
      left_df$NB2_NB1,
      na.rm = TRUE
    ),

    right_gap = median(
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

    left_alpha = median(
      left_df$alpha_mu,
      na.rm = TRUE
    ),

    right_alpha = median(
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
# METHOD-LEVEL VALIDATION
# =============================================================================

validate_method_level <- function(
    feature_df,
    zero_df,
    anchor,
    selected_fraction,
    selected_stability_row,
    region_summary,
    total_n) {

  if (
    !isTRUE(
      selected_stability_row$pass_stability
    )
  ) {
    stop(
      "Validation failed: selected fraction did not pass stability criteria."
    )
  }

  if (
    !is.finite(anchor) ||
    anchor < 1L ||
    anchor > total_n
  ) {
    stop(
      "Validation failed: Anchor is outside the ranked feature series."
    )
  }

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
      "Validation failed: Anchor is not a rounded spline d2 zero crossing."
    )
  }

  right_idx <- seq.int(
    anchor,
    total_n
  )

  right_n <- length(
    right_idx
  )

  actual_right_fraction <- right_n /
    total_n

  if (
    actual_right_fraction +
      1e-12 <
      selected_fraction
  ) {
    stop(
      "Validation failed: RIGHT is smaller than the selected minimum fraction."
    )
  }

  left_end <- anchor -
    1L

  left_start <- left_end -
    right_n +
    1L

  if (left_start < 1L) {
    stop(
      "Validation failed: matched LEFT block extends below rank 1."
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
      "Validation failed: LEFT and RIGHT blocks are not equal in size."
    )
  }

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
    left_n = nrow(left_df),

    right_n = nrow(right_df),

    left_NB2 = median(
      left_df$NB2,
      na.rm = TRUE
    ),

    right_NB2 = median(
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

    left_gap = median(
      left_df$NB2_NB1,
      na.rm = TRUE
    ),

    right_gap = median(
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

    left_alpha = median(
      left_df$alpha_mu,
      na.rm = TRUE
    ),

    right_alpha = median(
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

  for (nm in names(expected)) {

    comparison <- all.equal(
      as.numeric(
        region_summary[[nm]][1L]
      ),
      as.numeric(
        expected[[nm]]
      ),
      tolerance = 1e-10
    )

    if (!isTRUE(comparison)) {
      stop(
        "Validation failed for summary field '",
        nm,
        "': ",
        comparison
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
    width = PNG_WIDTH_IN,
    height = PNG_HEIGHT_IN,
    units = "in",
    res = PNG_DPI,
    bg = "white"
  )

  grid.newpage()

  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 3L,
        ncol = 1L,
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
        layout.pos.row = i,
        layout.pos.col = 1L
      )
    )
  }

  dev.off()
}


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
    stat_x <- left_x +
      60L
  }

  list(
    left_x = left_x,
    stat_x = stat_x
  )
}


# =============================================================================
# THREE-PANEL FIGURE
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

  left_min <- region_summary$left_start
  left_max <- region_summary$left_end

  right_min <- region_summary$right_start
  right_max <- region_summary$right_end


  # ---------------------------------------------------------------------------
  # Curvature-transition coordinates
  # ---------------------------------------------------------------------------

  if (nrow(zero_df) > 0L) {

    zero_y <- approx(
      x = variance_df$rank,
      y = variance_df$smooth_log1p_empirical_variance,
      xout = zero_df$crossing_rank,
      rule = 2
    )$y

    zero_plot_df <- data.frame(
      rank = zero_df$crossing_rank,
      y = zero_y,
      stringsAsFactors = FALSE
    )

  } else {

    zero_plot_df <- data.frame(
      rank = numeric(0),
      y = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  anchor_y <- approx(
    x = variance_df$rank,
    y = variance_df$smooth_log1p_empirical_variance,
    xout = anchor,
    rule = 2
  )$y


  # ---------------------------------------------------------------------------
  # Long-form NB2 metrics
  # ---------------------------------------------------------------------------

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
      names_to = "metric",
      values_to = "value"
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
    total_n = total_n,
    anchor = anchor
  )

  variance_range <- range(
    variance_df$smooth_log1p_empirical_variance,
    finite = TRUE
  )

  variance_span <- diff(
    variance_range
  )

  if (
    !is.finite(variance_span) ||
    variance_span <= 0
  ) {
    variance_span <- 1
  }

  top_y <- variance_range[2L]

  mid_range <- range(
    nb_long$value,
    finite = TRUE
  )

  mid_span <- diff(
    mid_range
  )

  if (
    !is.finite(mid_span) ||
    mid_span <= 0
  ) {
    mid_span <- 1
  }

  mid_top <- mid_range[2L]


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
    "Minimum tail = ",
    selected_fraction_pct,
    "%\n",
    "Observed RIGHT = ",
    actual_fraction_pct,
    "%\n",
    "LEFT n = ",
    region_summary$left_n,
    "\n",
    "RIGHT n = ",
    region_summary$right_n,
    "\n",
    "Recovery = ",
    round(
      selected_stability_row$anchor_recovery_rate,
      3
    ),
    "\n",
    "Jaccard = ",
    round(
      selected_stability_row$median_jaccard,
      3
    )
  )


  box3 <- paste(
    "RIGHT = Anchor to end",
    "LEFT = matched block",
    "NB2 metrics evaluated after selection",
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
  # PANEL 1: DATA-DERIVED GEOMETRY
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
      xmin = left_min,
      xmax = left_max,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$left_fill,
      alpha = 0.70
    ) +

    annotate(
      "rect",
      xmin = right_min,
      xmax = right_max,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$right_fill,
      alpha = 0.70
    ) +

    geom_line(
      color = COL$var_curve,
      linewidth = 1.0
    ) +

    geom_point(
      data = zero_plot_df,
      mapping = aes(
        x = rank,
        y = y
      ),
      inherit.aes = FALSE,
      color = COL$zero_crossing,
      size = 1.2,
      alpha = 0.65
    ) +

    geom_vline(
      xintercept = anchor,
      color = COL$anchor,
      linewidth = 1.0
    ) +

    annotate(
      "point",
      x = anchor,
      y = anchor_y,
      color = COL$anchor,
      shape = 16,
      size = 3.5
    ) +

    annotate(
      "label",
      x = pos$left_x,
      y = top_y -
        0.04 *
        variance_span,
      label = box1,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.25,
      fill = grDevices::adjustcolor(
        "white",
        alpha.f = 0.96
      )
    ) +

    annotate(
      "label",
      x = pos$stat_x,
      y = top_y -
        0.30 *
        variance_span,
      label = box2,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.25,
      fill = grDevices::adjustcolor(
        "white",
        alpha.f = 0.96
      )
    ) +

    labs(
      title = paste0(
        comparison_name,
        " ",
        arm_name,
        ": geometry"
      ),
      subtitle = paste0(
        "Bootstrap-stable curvature cutoff; grey points denote ",
        "spline second-derivative zero crossings"
      ),
      x = "Rank",
      y = "Smoothed log(1 + variance)"
    ) +

    theme_bw(
      base_size = 11
    ) +

    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )


  # ---------------------------------------------------------------------------
  # PANEL 2: NB2-RELATED CORROBORATION
  # ---------------------------------------------------------------------------

  p2 <- ggplot() +

    annotate(
      "rect",
      xmin = left_min,
      xmax = left_max,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$left_fill,
      alpha = 0.70
    ) +

    annotate(
      "rect",
      xmin = right_min,
      xmax = right_max,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$right_fill,
      alpha = 0.70
    ) +

    geom_vline(
      xintercept = anchor,
      color = COL$anchor,
      linewidth = 0.8
    ) +

    geom_line(
      data = nb_long,
      mapping = aes(
        rank,
        value,
        color = metric
      ),
      linewidth = 0.95
    ) +

    annotate(
      "label",
      x = pos$left_x,
      y = mid_top -
        0.04 *
        mid_span,
      label = box3,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.25,
      fill = grDevices::adjustcolor(
        "white",
        alpha.f = 0.96
      )
    ) +

    annotate(
      "label",
      x = pos$stat_x,
      y = mid_top -
        0.30 *
        mid_span,
      label = box4,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.25,
      fill = grDevices::adjustcolor(
        "white",
        alpha.f = 0.96
      )
    ) +

    scale_color_manual(
      values = TRACE_COLORS,
      breaks = TRACE_LEVELS
    ) +

    labs(
      title = paste0(
        comparison_name,
        " ",
        arm_name,
        ": corroboration"
      ),
      subtitle = paste0(
        "NB2-related quantities are evaluated only after ",
        "geometric cutoff selection"
      ),
      x = "Rank",
      y = "NB2-related signal",
      color = NULL
    ) +

    theme_bw(
      base_size = 11
    ) +

    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )


  # ---------------------------------------------------------------------------
  # PANEL 3: MATCHED REGION SUMMARY
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
        x = LEFT,
        xend = RIGHT,
        yend = metric
      ),
      color = "#7A7A7A",
      linewidth = 0.8
    ) +

    geom_point(
      aes(
        x = LEFT,
        color = "LEFT"
      ),
      size = 3.4
    ) +

    geom_point(
      aes(
        x = RIGHT,
        color = "RIGHT"
      ),
      size = 3.4
    ) +

    scale_color_manual(
      values = REGION_COLORS,
      breaks = REGION_LEVELS
    ) +

    labs(
      title = paste0(
        comparison_name,
        " ",
        arm_name,
        ": summary"
      ),
      subtitle = paste0(
        "Matched-region medians; farther right indicates ",
        "stronger NB2-related signal"
      ),
      x = "Median",
      y = NULL,
      color = NULL
    ) +

    theme_bw(
      base_size = 11
    ) +

    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  save_three_panel_plot(
    plot_list = list(
      p1,
      p2,
      p3
    ),
    filename = out_file
  )
}


# =============================================================================
# RUN ONE ANALYSIS TRACK
# =============================================================================

run_one_track <- function(
    comparison_name,
    arm_name,
    count_mat_arm,
    prepare_function,
    output_dir,
    track,
    fig_tag,
    rank_tag,
    metric_tag,
    metric_name,
    seed) {

  total_n <- nrow(
    count_mat_arm
  )

  prefix <- paste0(
    comparison_name,
    "_",
    arm_name,
    "_",
    track
  )

  stability_summary_path <- file.path(
    output_dir,
    paste0(
      "Table_StabilitySummary_",
      prefix,
      ".csv"
    )
  )

  stability_bootstrap_path <- file.path(
    output_dir,
    paste0(
      "Table_StabilityBootstrap_",
      prefix,
      ".csv"
    )
  )

  candidate_path <- file.path(
    output_dir,
    paste0(
      "Table_CandidateAnchors_",
      prefix,
      ".csv"
    )
  )

  zero_path <- file.path(
    output_dir,
    paste0(
      "Table_Zero_",
      prefix,
      ".csv"
    )
  )

  valid_path <- file.path(
    output_dir,
    paste0(
      "Table_Valid_",
      prefix,
      ".csv"
    )
  )

  rank_path <- file.path(
    output_dir,
    paste0(
      "Table_Rank_",
      prefix,
      ".csv"
    )
  )

  cut_path <- file.path(
    output_dir,
    paste0(
      "Table_Cutoff_",
      prefix,
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
  # Bootstrap stability
  # ---------------------------------------------------------------------------

  stability_obj <- assess_geometry_stability(
    count_mat_arm = count_mat_arm,
    prepare_function = prepare_function,
    candidate_fractions = MIN_FRACTION_GRID,
    bootstrap_n = BOOTSTRAP_N,
    spar = VAR_SPLINE_SPAR,
    seed = seed,
    cores = BOOTSTRAP_CORES
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

  write.csv(
    stability_obj$original_candidates,
    candidate_path,
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Full-data zero crossings are written regardless of whether stability passes
  # ---------------------------------------------------------------------------

  zero_out <- stability_obj$original_geometry$zero_df

  if (nrow(zero_out) > 0L) {
    zero_out$crossing_rank_rounded <- as.integer(
      round(
        zero_out$crossing_rank
      )
    )
  } else {
    zero_out$crossing_rank_rounded <- integer(0)
  }

  write.csv(
    zero_out,
    zero_path,
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # No stable solution: stop without forcing a cutoff
  # ---------------------------------------------------------------------------

  if (
    is.null(
      stability_obj$selected
    )
  ) {

    validation_row <- data.frame(
      comp = comparison_name,
      arm = arm_name,
      track = track,
      Status = "UNSTABLE",
      Reason = paste0(
        "No candidate minimum fraction satisfied all ",
        "prespecified bootstrap stability criteria."
      ),
      stringsAsFactors = FALSE
    )

    write.csv(
      validation_row,
      valid_path,
      row.names = FALSE
    )

    fail_summary <- data.frame(
      comp = comparison_name,
      arm = arm_name,
      track = track,
      rank_method = rank_tag,
      metric_matrix = metric_name,
      Status = "UNSTABLE",
      SelectedMinFraction = NA_real_,
      Anchor = NA_integer_,
      ActualRightFraction = NA_real_,
      LeftSize = NA_integer_,
      RightSize = NA_integer_,
      left_NB2 = NA_real_,
      right_NB2 = NA_real_,
      diff_NB2 = NA_real_,
      left_gap = NA_real_,
      right_gap = NA_real_,
      diff_gap = NA_real_,
      left_alpha = NA_real_,
      right_alpha = NA_real_,
      diff_alpha = NA_real_,
      Call_NB2 = NA_character_,
      Call_Gap = NA_character_,
      Call_Alpha = NA_character_,
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
        summary = fail_summary,
        selected = NULL,
        stability_summary = stability_obj$stability_summary,
        prepared = stability_obj$original_prepared
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Stable solution
  # ---------------------------------------------------------------------------

  selected_stability <- stability_obj$selected

  selected_fraction <- as.numeric(
    selected_stability$min_fraction
  )

  anchor <- as.integer(
    selected_stability$anchor
  )

  geometry <- stability_obj$original_geometry
  prepared <- stability_obj$original_prepared

  rank_order <- geometry$rank_order
  variance_df <- geometry$variance_df
  zero_df <- geometry$zero_df

  feature_df <- compute_ranked_feature_metrics(
    metric_mat_arm = prepared$metric_matrix,
    rank_order = rank_order
  )

  feature_df$abs_pc1_loading <- geometry$abs_loadings[
    rank_order
  ]

  region_summary <- summarize_regions(
    feature_df = feature_df,
    anchor_rank = anchor,
    total_n = total_n
  )


  # ---------------------------------------------------------------------------
  # Validate final result
  # ---------------------------------------------------------------------------

  validate_method_level(
    feature_df = feature_df,
    zero_df = zero_df,
    anchor = anchor,
    selected_fraction = selected_fraction,
    selected_stability_row = selected_stability,
    region_summary = region_summary,
    total_n = total_n
  )


  # ---------------------------------------------------------------------------
  # Selected-cutoff metadata
  # ---------------------------------------------------------------------------

  selected_df <- data.frame(
    comp = comparison_name,
    arm = arm_name,
    track = track,

    rank_method = rank_tag,
    metric_matrix = metric_name,

    Status = "PASS",

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

    RightFractionIQR =
      selected_stability$right_fraction_iqr,

    MedianJaccard =
      selected_stability$median_jaccard,

    stringsAsFactors = FALSE
  )

  cutoff_summary <- bind_cols(
    selected_df,
    region_summary
  ) %>%
    mutate(
      LeftSize = left_n,
      RightSize = right_n,

      Call_NB2 = ifelse(
        diff_NB2 > 0,
        "RIGHT",
        "NOT_RIGHT"
      ),

      Call_Gap = ifelse(
        diff_gap > 0,
        "RIGHT",
        "NOT_RIGHT"
      ),

      Call_Alpha = ifelse(
        diff_alpha > 0,
        "RIGHT",
        "NOT_RIGHT"
      )
    )


  # ---------------------------------------------------------------------------
  # Mark selected zero crossing
  # ---------------------------------------------------------------------------

  zero_selected <- zero_df

  if (nrow(zero_selected) > 0L) {

    zero_selected$crossing_rank_rounded <- as.integer(
      round(
        zero_selected$crossing_rank
      )
    )

    zero_selected$selected_anchor <- (
      zero_selected$crossing_rank_rounded ==
        anchor
    )

  } else {

    zero_selected$crossing_rank_rounded <- integer(0)
    zero_selected$selected_anchor <- logical(0)
  }

  write.csv(
    zero_selected,
    zero_path,
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Validation output
  # ---------------------------------------------------------------------------

  validation_row <- data.frame(
    comp = comparison_name,
    arm = arm_name,
    track = track,

    Status = "PASS",

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

    RightFractionIQR =
      selected_stability$right_fraction_iqr,

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
  )

  write.csv(
    validation_row,
    valid_path,
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Rank-level output
  # ---------------------------------------------------------------------------

  rank_output <- feature_df %>%
    left_join(
      variance_df,
      by = "rank",
      suffix = c(
        "_feature_metric",
        "_geometry"
      )
    )

  write.csv(
    rank_output,
    rank_path,
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Cutoff summary
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
    comparison_name = comparison_name,
    arm_name = arm_name,
    variance_df = variance_df,
    zero_df = zero_df,
    feature_df = feature_df,
    anchor = anchor,
    selected_stability_row = selected_stability,
    region_summary = region_summary,
    out_file = fig_path,
    fig_tag = fig_tag,
    rank_tag = rank_tag,
    metric_tag = metric_tag
  )


  message(
    "[",
    track,
    "] ",
    comparison_name,
    " ",
    arm_name,
    " | min_fraction=",
    sprintf(
      "%.3f",
      selected_fraction
    ),
    " | Anchor=",
    anchor,
    " | RIGHT n=",
    region_summary$right_n,
    " | valid=",
    sprintf(
      "%.3f",
      selected_stability$bootstrap_valid_rate
    ),
    " | recovery=",
    sprintf(
      "%.3f",
      selected_stability$anchor_recovery_rate
    ),
    " | anchor_IQR/N=",
    sprintf(
      "%.4f",
      selected_stability$anchor_iqr_fraction
    ),
    " | RIGHT_fraction_IQR=",
    sprintf(
      "%.4f",
      selected_stability$right_fraction_iqr
    ),
    " | Jaccard=",
    sprintf(
      "%.3f",
      selected_stability$median_jaccard
    )
  )


  list(
    summary = cutoff_summary,
    selected = selected_df,
    stability_summary = stability_obj$stability_summary,
    prepared = prepared
  )
}


# =============================================================================
# RUN
# =============================================================================

count_mat <- read_count_matrix(
  path = COUNT_FILE,
  comparisons = COMPARISONS
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

message(
  "Bootstrap worker processes: ",
  BOOTSTRAP_CORES
)

message(
  "Stability criteria:"
)

message(
  "  valid bootstrap rate >= ",
  MIN_BOOTSTRAP_VALID_RATE
)

message(
  "  Anchor recovery >= ",
  MIN_ANCHOR_RECOVERY_RATE,
  " within ",
  ANCHOR_TOLERANCE_FRACTION *
    100,
  "% of rank axis"
)

message(
  "  Anchor IQR/N <= ",
  MAX_ANCHOR_IQR_FRACTION
)

message(
  "  RIGHT-fraction IQR <= ",
  MAX_RIGHT_FRACTION_IQR
)

message(
  "  median Jaccard >= ",
  MIN_MEDIAN_JACCARD
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

  pats <- COMPARISONS[[comparison_name]]


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

    if (length(sample_idx) < 2L) {
      stop(
        "Not enough samples for ",
        comparison_name,
        " ",
        arm_name,
        ". Matched columns: ",
        length(sample_idx)
      )
    }

    count_mat_arm <- count_mat[
      ,
      sample_idx,
      drop = FALSE
    ]

    message(
      "------------------------------------------------------------"
    )

    message(
      comparison_name,
      " ",
      arm_name,
      ": ",
      ncol(count_mat_arm),
      " samples"
    )


    # =========================================================================
    # MAIN MANUSCRIPT TRACK
    # =========================================================================

    main_seed <- seed_from_label(
      paste(
        comparison_name,
        arm_name,
        "Main",
        sep = "_"
      )
    )

    main_prepare_function <- function(x) {
      prepare_main_matrices(x)
    }

    main_res <- run_one_track(
      comparison_name = comparison_name,
      arm_name = arm_name,
      count_mat_arm = count_mat_arm,
      prepare_function = main_prepare_function,
      output_dir = comp_dir,
      track = "Main",
      fig_tag = "Main",
      rank_tag = "Rank: CPM log1p",
      metric_tag = "Metrics: raw counts",
      metric_name = "raw_counts",
      seed = main_seed
    )

    overall_rows[[
      length(overall_rows) +
        1L
    ]] <- main_res$summary

    overall_stability_rows[[
      length(overall_stability_rows) +
        1L
    ]] <- main_res$stability_summary %>%
      mutate(
        comp = comparison_name,
        arm = arm_name,
        track = "Main",
        .before = 1L
      )


    # =========================================================================
    # DESEQ2 SUPPLEMENTARY TRACK
    # =========================================================================

    if (RUN_DESEQ2_SUPPLEMENT) {

      if (
        !requireNamespace(
          "DESeq2",
          quietly = TRUE
        ) ||
        !requireNamespace(
          "SummarizedExperiment",
          quietly = TRUE
        )
      ) {

        message(
          "[DESeq2] skipped for ",
          comparison_name,
          " ",
          arm_name,
          ": required package not available."
        )

      } else {

        deseq2_seed <- seed_from_label(
          paste(
            comparison_name,
            arm_name,
            "DESeq2",
            sep = "_"
          )
        )

        deseq2_prepare_function <- function(x) {

          prepare_deseq2_matrices(
            count_mat_arm = x,
            rank_method = DESEQ2_RANK_METHOD
          )
        }

        deseq2_rank_tag <- if (
          DESEQ2_RANK_METHOD ==
            "vst"
        ) {

          "Rank: DESeq2 VST"

        } else {

          "Rank: DESeq2 log1p"
        }

        deseq2_res <- run_one_track(
          comparison_name = comparison_name,
          arm_name = arm_name,
          count_mat_arm = count_mat_arm,
          prepare_function = deseq2_prepare_function,
          output_dir = comp_dir,
          track = "DESeq2",
          fig_tag = "DESeq2 supplement",
          rank_tag = deseq2_rank_tag,
          metric_tag = "Metrics: DESeq2 normalized",
          metric_name = "deseq2_normalized_counts",
          seed = deseq2_seed
        )


        # ---------------------------------------------------------------------
        # Full-data DESeq2 size factors
        # ---------------------------------------------------------------------

        if (
          !is.null(
            deseq2_res$prepared$size_factors
          )
        ) {

          size_factors <- deseq2_res$prepared$size_factors

          write.csv(
            data.frame(
              sample = names(
                size_factors
              ),
              size_factor = as.numeric(
                size_factors
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
        }


        overall_rows[[
          length(overall_rows) +
            1L
        ]] <- deseq2_res$summary

        overall_stability_rows[[
          length(overall_stability_rows) +
            1L
        ]] <- deseq2_res$stability_summary %>%
          mutate(
            comp = comparison_name,
            arm = arm_name,
            track = "DESeq2",
            .before = 1L
          )
      }
    }
  }
}


# =============================================================================
# OVERALL OUTPUT TABLES
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


# =============================================================================
# FINAL COMPLETION MESSAGE
# =============================================================================

message(
  "============================================================"
)

message(
  "Analysis complete."
)

message(
  "Outputs written to: ",
  OUT_ROOT
)

message(
  "No fixed 5,000-feature reference was used."
)

message(
  "Cutoffs were selected using variance geometry and bootstrap ",
  "reproducibility before NB2-related corroboration."
)

message(
  "Tracks failing all prespecified stability criteria were retained ",
  "as UNSTABLE rather than being assigned a forced cutoff."
)

message(
  "============================================================"
)
