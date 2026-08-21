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
# Within each experimental arm, features are ranked from lowest to highest
# absolute loading on the first principal component (|PC1 loading|).
#
# PRIMARY ANALYSIS
#   Ranking matrix:
#     log1p-transformed counts-per-million (CPM) values.
#
#   Variance / corroboration matrix:
#     raw sequencing counts.
#
# SUPPLEMENTARY ANALYSIS
#   Ranking matrix:
#     log1p-transformed DESeq2-normalized counts, or VST if requested.
#
#   Variance / corroboration matrix:
#     DESeq2-normalized counts.
#
# ---------------------------------------------------------------------------
# DATA-DERIVED GEOMETRIC CUTOFF
# ---------------------------------------------------------------------------
#
# The previous fixed 5,000-feature reference is removed completely.
#
# After ranking by |PC1 loading|, empirical feature variance is calculated
# along the ranked axis and transformed as:
#
#     log(1 + variance)
#
# A smoothing spline is fitted to this trajectory. Its second derivative is
# evaluated over a dense rank grid. Locations at which the second derivative
# changes sign are treated as candidate curvature transitions.
#
# Candidate minimum high-loading tail fractions are prespecified as:
#
#     5%, 7.5%, 10%, 12.5%, 15%
#
# For each candidate fraction f, an eligible Anchor must:
#
#   1. correspond to a spline second-derivative zero crossing;
#   2. leave at least f*N features from Anchor through the highest-loading
#      right edge of the ranked series; and
#   3. permit construction of an immediately adjacent LEFT block containing
#      exactly the same number of features as RIGHT.
#
# Among eligible crossings, the RIGHTMOST crossing is selected.
#
# Thus:
#
#     RIGHT = Anchor through rank N
#
#     LEFT  = immediately preceding equal-sized block
#
# ---------------------------------------------------------------------------
# BOOTSTRAP STABILITY
# ---------------------------------------------------------------------------
#
# Cutoff stability is assessed by bootstrap resampling of samples within each
# experimental arm.
#
# For every bootstrap replicate, the following are recomputed:
#
#     sample resampling
#          ->
#     PC1
#          ->
#     |PC1 loading| ranking
#          ->
#     empirical feature variance
#          ->
#     smoothing spline
#          ->
#     spline second derivative
#          ->
#     zero crossings
#          ->
#     candidate Anchor
#
# For the DESeq2 supplementary track, DESeq2 normalization is estimated from
# the full arm once. Bootstrap resampling is then performed on those
# track-specific normalized matrices. This evaluates cutoff stability
# conditional on the specified normalization procedure while avoiding
# re-estimation of normalization parameters within every bootstrap sample.
#
# A candidate minimum fraction is considered stable only if ALL of the
# following prespecified criteria are satisfied:
#
#   1. Valid Anchor in >= 90% of bootstrap replicates.
#
#   2. Among valid bootstrap replicates, >= 80% recover an Anchor within
#      3% of the total rank axis of the full-data Anchor.
#
#   3. Bootstrap Anchor IQR <= 3% of the total rank axis.
#
#   4. Bootstrap RIGHT-region fraction IQR <= 3 percentage points.
#
#   5. Median Jaccard overlap between bootstrap and full-data RIGHT-feature
#      membership >= 0.80.
#
# The SMALLEST candidate minimum fraction satisfying all criteria is selected.
#
# The final Anchor is the corresponding Anchor obtained from the full dataset.
# The bootstrap median is NOT substituted for the full-data Anchor.
#
# If no candidate fraction satisfies all criteria, that analysis track is
# reported as UNSTABLE. No cutoff is forced.
#
# ---------------------------------------------------------------------------
# NB2-RELATED CORROBORATION
# ---------------------------------------------------------------------------
#
# NB2-related quantities are calculated ONLY AFTER the cutoff has been selected.
#
# Let:
#
#     mu       = empirical feature mean
#     variance = empirical feature variance
#
# Then:
#
#     NB2 =
#       log(1 + max(variance - mu, 0))
#
#     NB2-NB1 =
#       log(1 + max(variance - mu, 0)) - log(1 + mu)
#
#     alpha =
#       max((variance - mu) / mu^2, 0)
#
#     alpha*mu =
#       log(1 + alpha*mu)
#
# These quantities are descriptive diagnostics of excess-variance behavior.
# They are not formal likelihood-ratio statistics.
#
# Critically, NB2, NB2-NB1, and alpha*mu are NOT used to select:
#
#     - the candidate minimum fraction;
#     - the Anchor;
#     - the spline smoothing parameter; or
#     - bootstrap stability.
#
# This maintains separation between geometric cutoff selection and downstream
# NB2-related corroboration.
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
# Candidate minimum high-loading tail fractions
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
# -----------------------------------------------------------------------------

BOOTSTRAP_N <- 500L

BOOTSTRAP_SEED_BASE <- 20260820L


DETECTED_CORES <- suppressWarnings(
  parallel::detectCores(
    logical = FALSE
  )
)

if (
  is.na(DETECTED_CORES) ||
  DETECTED_CORES < 2L
) {

  BOOTSTRAP_CORES <- 1L

} else {

  BOOTSTRAP_CORES <- min(
    4L,
    DETECTED_CORES - 1L
  )
}


# -----------------------------------------------------------------------------
# Prespecified stability thresholds
# -----------------------------------------------------------------------------

MIN_BOOTSTRAP_VALID_RATE <- 0.90

ANCHOR_TOLERANCE_FRACTION <- 0.03

MIN_ANCHOR_RECOVERY_RATE <- 0.80

MAX_ANCHOR_IQR_FRACTION <- 0.03

MAX_RIGHT_FRACTION_IQR <- 0.03

MIN_MEDIAN_JACCARD <- 0.80


# -----------------------------------------------------------------------------
# Figure settings
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
#
#   "normalized_log1p"
#   "vst"


# -----------------------------------------------------------------------------
# Experimental comparisons
# -----------------------------------------------------------------------------

COMPARISONS <- list(

  RT0_ZT6 =
    list(
      control = "^R0_",
      treatment = "^ZT6_"
    ),

  RT2_ZT8 =
    list(
      control = "^R2_",
      treatment = "^ZT8_"
    ),

  RT4_ZT10 =
    list(
      control = "^R4_",
      treatment = "^ZT10_"
    ),

  RT8_ZT14 =
    list(
      control = "^R8_",
      treatment = "^ZT14_"
    )
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

  nb2 = "#1B9E77",

  nb_gap = "#CC1E8C",

  alpha_mu = "#386CB0",

  left_fill = "#CBE3F8",

  right_fill = "#DDF2D5",

  anchor = "#000000",

  zero_crossing = "#777777",

  left_pt = "#5B8FD1",

  right_pt = "#43A047"
)


TRACE_LEVELS <- c(
  "NB2",
  "NB2-NB1",
  "alpha*mu"
)


TRACE_COLORS <- c(

  "NB2" =
    COL$nb2,

  "NB2-NB1" =
    COL$nb_gap,

  "alpha*mu" =
    COL$alpha_mu
)


REGION_LEVELS <- c(
  "LEFT",
  "RIGHT"
)


REGION_COLORS <- c(

  "LEFT" =
    COL$left_pt,

  "RIGHT" =
    COL$right_pt
)


# =============================================================================
# INPUT
# =============================================================================

read_count_matrix <- function(
    path,
    comparisons) {

  raw_df <- read.csv(

    path,

    check.names = FALSE,

    stringsAsFactors = FALSE
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


  # ---------------------------------------------------------------------------
  # Identify actual experimental sample columns using the same sample patterns
  # used by the analysis. This prevents unrelated metadata columns from being
  # included as samples.
  # ---------------------------------------------------------------------------

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


  if (
    length(sample_idx) == 0L
  ) {

    stop(
      "No sample columns matched the patterns specified in COMPARISONS."
    )
  }


  if (
    1L %in% sample_idx
  ) {

    stop(
      "Column 1 matched a sample pattern; column 1 must contain feature IDs."
    )
  }


  # ---------------------------------------------------------------------------
  # Convert only the actual sample columns to numeric.
  #
  # The original working analysis converted non-finite entries to zero. That
  # behavior is preserved here so the revised cutoff method remains compatible
  # with the same input file that was successfully analyzed previously.
  # ---------------------------------------------------------------------------

  count_df <- raw_df[
    ,
    sample_idx,
    drop = FALSE
  ]


  parsed_columns <- lapply(

    count_df,

    function(x) {

      suppressWarnings(

        as.numeric(

          trimws(
            as.character(x)
          )
        )
      )
    }
  )


  count_mat <- do.call(
    cbind,
    parsed_columns
  )


  colnames(count_mat) <- colnames(
    count_df
  )


  storage.mode(count_mat) <- "numeric"


  bad_count_n <- sum(
    !is.finite(count_mat)
  )


  if (
    bad_count_n > 0L
  ) {

    message(

      "Replacing ",

      bad_count_n,

      " non-finite count entries with 0 ",
      "(same behavior as the original working script)."
    )


    count_mat[
      !is.finite(count_mat)
    ] <- 0
  }


  count_mat <- pmax(
    count_mat,
    0
  )


  # ---------------------------------------------------------------------------
  # Feature identifiers
  #
  # The prior script did not require every feature-ID field to be populated.
  # Blank identifiers therefore do not cause the analysis to terminate.
  #
  # Instead, blank identifiers receive deterministic row-based labels. This
  # preserves those rows and provides stable identities for bootstrap Jaccard
  # calculations.
  # ---------------------------------------------------------------------------

  feature_ids <- trimws(

    as.character(
      raw_df[[1L]]
    )
  )


  blank_feature_id <- (
    is.na(feature_ids) |
    feature_ids == ""
  )


  if (
    any(blank_feature_id)
  ) {

    feature_ids[
      blank_feature_id
    ] <- paste0(

      "__feature_row_",

      which(
        blank_feature_id
      )
    )


    message(

      "Assigned deterministic row IDs to ",

      sum(
        blank_feature_id
      ),

      " blank feature identifier(s)."
    )
  }


  feature_ids <- make.unique(

    feature_ids,

    sep = "__dup_"
  )


  rownames(count_mat) <- feature_ids


  # ---------------------------------------------------------------------------
  # Remove features having zero counts in all selected samples.
  # ---------------------------------------------------------------------------

  keep <- rowSums(
    count_mat
  ) > 0


  count_mat <- count_mat[
    keep,
    ,
    drop = FALSE
  ]


  if (
    nrow(count_mat) < 2L
  ) {

    stop(
      "Fewer than two nonzero features remain after input filtering."
    )
  }


  count_mat
}


# =============================================================================
# CPM NORMALIZATION
# =============================================================================

normalize_cpm_log1p <- function(
    count_mat_arm) {

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

    2L,

    lib_sizes / 1e6,

    "/"
  )


  log1p(cpm)
}


# =============================================================================
# DESEQ2 NORMALIZATION
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


  if (
    !(
      rank_method %in%
      c(
        "normalized_log1p",
        "vst"
      )
    )
  ) {

    stop(
      "DESEQ2_RANK_METHOD must be 'normalized_log1p' or 'vst'."
    )
  }


  col_data <- data.frame(

    row.names =
      colnames(
        count_mat_arm
      ),

    intercept =
      factor(

        rep(
          "one",
          ncol(count_mat_arm)
        )
      )
  )


  dds <- DESeq2::DESeqDataSetFromMatrix(

    countData =
      round(
        count_mat_arm
      ),

    colData =
      col_data,

    design =
      ~ 1
  )


  dds <- DESeq2::estimateSizeFactors(
    dds
  )


  norm_counts <- DESeq2::counts(

    dds,

    normalized = TRUE
  )


  if (
    rank_method == "vst"
  ) {

    vst_obj <- tryCatch(

      DESeq2::vst(
        dds,
        blind = TRUE
      ),

      error =
        function(e) {
          NULL
        }
    )


    if (
      is.null(vst_obj)
    ) {

      ranking_matrix <- log1p(
        norm_counts
      )


      rank_method_used <-
        "DESeq2 log1p (VST fallback)"

    } else {

      ranking_matrix <-
        SummarizedExperiment::assay(
          vst_obj
        )


      rank_method_used <-
        "DESeq2 VST"
    }

  } else {

    ranking_matrix <- log1p(
      norm_counts
    )


    rank_method_used <-
      "DESeq2 log1p"
  }


  list(

    normalized_counts =
      norm_counts,

    ranking_matrix =
      ranking_matrix,

    size_factors =
      DESeq2::sizeFactors(
        dds
      ),

    rank_method_used =
      rank_method_used
  )
}


# =============================================================================
# FAST ROW VARIANCE
# =============================================================================

row_variance_fast <- function(
    mat) {

  n <- ncol(
    mat
  )


  if (
    n < 2L
  ) {

    stop(
      "At least two samples are required for variance estimation."
    )
  }


  mu <- rowMeans(
    mat
  )


  ss <- rowSums(
    mat * mat
  ) -
    n *
    mu^2


  # Numerical roundoff can create very small negative sums of squares.

  ss[
    ss < 0 &
    abs(ss) < 1e-8
  ] <- 0


  out <- ss /
    (
      n - 1L
    )


  out[
    !is.finite(out)
  ] <- 0


  pmax(
    out,
    0
  )
}


# =============================================================================
# PC1
# =============================================================================

compute_abs_pc1_loadings <- function(
    norm_mat_arm) {

  if (
    nrow(norm_mat_arm) < 2L ||
    ncol(norm_mat_arm) < 2L
  ) {

    stop(
      "PC1 requires at least two features and two samples."
    )
  }


  # ---------------------------------------------------------------------------
  # Matrix orientation:
  #
  # rows    = samples
  # columns = features
  #
  # PC1 is calculated through the sample-space Gram matrix. Because the number
  # of samples is much smaller than the number of features, this is faster than
  # repeatedly calling prcomp() during bootstrap analysis while yielding the
  # same leading loading vector up to sign.
  # ---------------------------------------------------------------------------

  X <- t(
    norm_mat_arm
  )


  X_centered <- sweep(

    X,

    2L,

    colMeans(X),

    "-"
  )


  gram <- tcrossprod(
    X_centered
  )


  eig <- eigen(

    gram,

    symmetric = TRUE
  )


  lambda1 <- eig$values[
    1L
  ]


  if (
    !is.finite(lambda1) ||
    lambda1 <=
    .Machine$double.eps
  ) {

    stop(
      "PC1 is undefined because between-sample variation is insufficient."
    )
  }


  u1 <- eig$vectors[
    ,
    1L
  ]


  loading <- as.numeric(

    crossprod(
      X_centered,
      u1
    )
  ) /
    sqrt(
      lambda1
    )


  names(loading) <- colnames(
    X_centered
  )


  loading[
    !is.finite(loading)
  ] <- 0


  abs(
    loading
  )
}


# =============================================================================
# RANKED VARIANCE CURVE
# =============================================================================

compute_ranked_variance_curve <- function(
    metric_mat_arm,
    rank_order,
    spar = VAR_SPLINE_SPAR) {

  empirical_var <- row_variance_fast(
    metric_mat_arm
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

    x =
      ranks,

    y =
      ranked_log_var,

    spar =
      spar
  )


  smooth_y <- as.numeric(

    stats::predict(

      spline_fit,

      x =
        ranks,

      deriv =
        0
    )$y
  )


  smooth_d2 <- as.numeric(

    stats::predict(

      spline_fit,

      x =
        ranks,

      deriv =
        2
    )$y
  )


  dense_n <- max(

    DENSE_GRID_MIN,

    length(ranks) *
    DENSE_GRID_MULTIPLIER
  )


  dense_x <- seq(

    min(ranks),

    max(ranks),

    length.out =
      dense_n
  )


  dense_y <- as.numeric(

    stats::predict(

      spline_fit,

      x =
        dense_x,

      deriv =
        0
    )$y
  )


  dense_d2 <- as.numeric(

    stats::predict(

      spline_fit,

      x =
        dense_x,

      deriv =
        2
    )$y
  )


  out <- data.frame(

    rank =
      ranks,

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


# =============================================================================
# SECOND-DERIVATIVE ZERO CROSSINGS
# =============================================================================

find_d2_zero_crossings <- function(
    dense_df) {

  x <- dense_df$dense_rank

  y <- dense_df$dense_d2


  ok <- (
    is.finite(x) &
    is.finite(y)
  )


  x <- x[
    ok
  ]

  y <- y[
    ok
  ]


  if (
    length(x) < 2L
  ) {

    return(

      data.frame(

        crossing_rank =
          numeric(0),

        crossing_type =
          character(0),

        stringsAsFactors = FALSE
      )
    )
  }


  a <- y[
    -length(y)
  ]


  b <- y[
    -1L
  ]


  xa <- x[
    -length(x)
  ]


  xb <- x[
    -1L
  ]


  sign_change_idx <- which(

    (
      a < 0 &
      b > 0
    ) |

    (
      a > 0 &
      b < 0
    )
  )


  crossings <- numeric(0)


  if (
    length(sign_change_idx) > 0L
  ) {

    frac <- abs(
      a[
        sign_change_idx
      ]
    ) /
      (
        abs(
          a[
            sign_change_idx
          ]
        ) +
        abs(
          b[
            sign_change_idx
          ]
        )
      )


    crossings <- xa[
      sign_change_idx
    ] +
      frac *
      (
        xb[
          sign_change_idx
        ] -
        xa[
          sign_change_idx
        ]
      )
  }


  # Exact zeros are uncommon but are retained explicitly.

  exact_zero_idx <- which(
    y == 0
  )


  if (
    length(exact_zero_idx) > 0L
  ) {

    crossings <- c(

      crossings,

      x[
        exact_zero_idx
      ]
    )
  }


  crossings <- crossings[
    is.finite(crossings)
  ]


  crossings <- sort(

    unique(

      round(
        crossings,
        8L
      )
    )
  )


  if (
    length(crossings) == 0L
  ) {

    return(

      data.frame(

        crossing_rank =
          numeric(0),

        crossing_type =
          character(0),

        stringsAsFactors = FALSE
      )
    )
  }


  data.frame(

    crossing_rank =
      crossings,

    crossing_type =
      "d2_sign_change",

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# COMPLETE GEOMETRY
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
      "Rank and metric matrices have different feature counts."
    )
  }


  if (
    ncol(rank_matrix) !=
    ncol(metric_matrix)
  ) {

    stop(
      "Rank and metric matrices have different sample counts."
    )
  }


  if (
    !identical(
      rownames(rank_matrix),
      rownames(metric_matrix)
    )
  ) {

    stop(
      "Feature IDs differ between rank and metric matrices."
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
# CANDIDATE ANCHOR SELECTION
# =============================================================================

select_anchor_for_min_fraction <- function(
    zero_df,
    total_n,
    min_fraction) {

  minimum_right_n <- ceiling(

    min_fraction *
    total_n
  )


  # Equal-sized LEFT block requires:
  #
  # Anchor - 1 >= N - Anchor + 1

  minimum_anchor_for_match <- ceiling(

    (
      total_n +
      2L
    ) /
    2
  )


  # At least minimum_right_n features must remain on the right.

  maximum_anchor_for_tail <- total_n -
    minimum_right_n +
    1L


  invalid_result <- function() {

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

      stringsAsFactors = FALSE
    )
  }


  if (
    nrow(zero_df) == 0L ||
    maximum_anchor_for_tail <
    minimum_anchor_for_match
  ) {

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
      invalid_result()
    )
  }


  # Rightmost eligible curvature transition.

  anchor <- max(
    eligible
  )


  right_n <- total_n -
    anchor +
    1L


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
      right_n,

    actual_right_fraction =
      right_n /
      total_n,

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# JACCARD SIMILARITY
# =============================================================================

jaccard_similarity <- function(
    set_a,
    set_b) {

  set_a <- unique(
    set_a
  )


  set_b <- unique(
    set_b
  )


  union_set <- union(
    set_a,
    set_b
  )


  if (
    length(union_set) == 0L
  ) {

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


# =============================================================================
# DETERMINISTIC RANDOM SEEDS
# =============================================================================

seed_from_label <- function(
    label,
    base_seed = BOOTSTRAP_SEED_BASE) {

  label_value <- sum(

    utf8ToInt(
      as.character(label)
    )
  )


  as.integer(

    (
      base_seed +
      label_value *
      1009L
    ) %%
    2147483647L
  )
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
    seed = BOOTSTRAP_SEED_BASE,
    cores = BOOTSTRAP_CORES) {

  candidate_fractions <- sort(

    unique(

      as.numeric(
        candidate_fractions
      )
    )
  )


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
  # Full-data RIGHT feature sets
  # ---------------------------------------------------------------------------

  original_right_sets <- setNames(

    vector(

      "list",

      length(
        candidate_fractions
      )
    ),

    sprintf(
      "%.6f",
      candidate_fractions
    )
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

      candidate_fractions[
        j
      ]
    )


    if (
      isTRUE(
        candidate_row$valid
      )
    ) {

      right_order_idx <-
        original_geometry$rank_order[

          seq.int(

            candidate_row$anchor,

            total_n
          )
        ]


      original_right_sets[[key]] <-
        rownames(rank_matrix)[
          right_order_idx
        ]

    } else {

      original_right_sets[[key]] <-
        character(0)
    }
  }


  # ---------------------------------------------------------------------------
  # Generate bootstrap sample draws before parallel processing so results are
  # deterministic regardless of worker count.
  # ---------------------------------------------------------------------------

  set.seed(
    seed
  )


  bootstrap_indices <- lapply(

    seq_len(
      bootstrap_n
    ),

    function(i) {

      sample.int(

        sample_n,

        size =
          sample_n,

        replace =
          TRUE
      )
    }
  )


  # ---------------------------------------------------------------------------
  # One bootstrap replicate
  # ---------------------------------------------------------------------------

  bootstrap_worker <- function(
      bootstrap_id) {

    sample_idx <- bootstrap_indices[[
      bootstrap_id
    ]]


    invalid_rows <- function(
        reason) {

      bind_rows(

        lapply(

          candidate_fractions,

          function(f) {

            data.frame(

              bootstrap =
                bootstrap_id,

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
                reason,

              stringsAsFactors = FALSE
            )
          }
        )
      )
    }


    # At least two distinct biological samples are required.

    if (
      length(
        unique(
          sample_idx
        )
      ) < 2L
    ) {

      return(

        invalid_rows(
          "fewer_than_2_unique_samples"
        )
      )
    }


    bootstrap_rank_matrix <- rank_matrix[
      ,
      sample_idx,
      drop = FALSE
    ]


    bootstrap_metric_matrix <- metric_matrix[
      ,
      sample_idx,
      drop = FALSE
    ]


    bootstrap_geometry <- tryCatch(

      compute_geometry(

        rank_matrix =
          bootstrap_rank_matrix,

        metric_matrix =
          bootstrap_metric_matrix,

        spar =
          spar
      ),

      error =
        function(e) {
          NULL
        }
    )


    if (
      is.null(
        bootstrap_geometry
      )
    ) {

      return(

        invalid_rows(
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

      f <- candidate_fractions[
        j
      ]


      selection <- select_anchor_for_min_fraction(

        zero_df =
          bootstrap_geometry$zero_df,

        total_n =
          total_n,

        min_fraction =
          f
      )


      if (
        !isTRUE(
          selection$valid
        )
      ) {

        result_rows[[j]] <- data.frame(

          bootstrap =
            bootstrap_id,

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


      right_order_idx <-
        bootstrap_geometry$rank_order[

          seq.int(

            selection$anchor,

            total_n
          )
        ]


      bootstrap_right_features <-
        rownames(rank_matrix)[
          right_order_idx
        ]


      key <- sprintf(
        "%.6f",
        f
      )


      full_right_features <-
        original_right_sets[[key]]


      if (
        length(
          full_right_features
        ) > 0L
      ) {

        jaccard_value <- jaccard_similarity(

          full_right_features,

          bootstrap_right_features
        )

      } else {

        jaccard_value <- NA_real_
      }


      result_rows[[j]] <- data.frame(

        bootstrap =
          bootstrap_id,

        min_fraction =
          f,

        valid =
          TRUE,

        anchor =
          selection$anchor,

        right_n =
          selection$right_n,

        right_fraction =
          selection$actual_right_fraction,

        jaccard =
          jaccard_value,

        invalid_reason =
          NA_character_,

        stringsAsFactors = FALSE
      )
    }


    bind_rows(
      result_rows
    )
  }


  # ---------------------------------------------------------------------------
  # Execute bootstrap
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

      mc.cores =
        cores,

      mc.preschedule =
        TRUE,

      mc.set.seed =
        FALSE
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
  # Stability summary
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

    f <- candidate_fractions[
      j
    ]


    original_row <- original_candidates[
      j,
      ,
      drop = FALSE
    ]


    bootstrap_subset <- bootstrap_df[

      abs(
        bootstrap_df$min_fraction -
        f
      ) < 1e-12,

      ,
      drop = FALSE
    ]


    valid_subset <- bootstrap_subset[

      bootstrap_subset$valid,

      ,
      drop = FALSE
    ]


    valid_bootstraps <- nrow(
      valid_subset
    )


    bootstrap_valid_rate <-
      valid_bootstraps /
      bootstrap_n


    if (
      valid_bootstraps > 0L
    ) {

      anchor_median <- median(
        valid_subset$anchor
      )


      anchor_iqr <- stats::IQR(
        valid_subset$anchor
      )


      anchor_iqr_fraction <-
        anchor_iqr /
        total_n


      right_fraction_median <- median(
        valid_subset$right_fraction
      )


      right_fraction_iqr <- stats::IQR(
        valid_subset$right_fraction
      )


      valid_jaccard <- valid_subset$jaccard[

        is.finite(
          valid_subset$jaccard
        )
      ]


      if (
        length(valid_jaccard) > 0L
      ) {

        median_jaccard <- median(
          valid_jaccard
        )

      } else {

        median_jaccard <- NA_real_
      }

    } else {

      anchor_median <- NA_real_

      anchor_iqr <- NA_real_

      anchor_iqr_fraction <- NA_real_

      right_fraction_median <- NA_real_

      right_fraction_iqr <- NA_real_

      median_jaccard <- NA_real_
    }


    anchor_tolerance_n <- ceiling(

      ANCHOR_TOLERANCE_FRACTION *
      total_n
    )


    if (
      isTRUE(
        original_row$valid
      ) &&
      valid_bootstraps > 0L
    ) {

      anchor_recovery_rate <- mean(

        abs(

          valid_subset$anchor -
          original_row$anchor
        ) <=
        anchor_tolerance_n
      )

    } else {

      anchor_recovery_rate <- NA_real_
    }


    pass_valid_rate <-
      is.finite(
        bootstrap_valid_rate
      ) &&
      bootstrap_valid_rate >=
      MIN_BOOTSTRAP_VALID_RATE


    pass_anchor_recovery <-
      is.finite(
        anchor_recovery_rate
      ) &&
      anchor_recovery_rate >=
      MIN_ANCHOR_RECOVERY_RATE


    pass_anchor_iqr <-
      is.finite(
        anchor_iqr_fraction
      ) &&
      anchor_iqr_fraction <=
      MAX_ANCHOR_IQR_FRACTION


    pass_right_fraction_iqr <-
      is.finite(
        right_fraction_iqr
      ) &&
      right_fraction_iqr <=
      MAX_RIGHT_FRACTION_IQR


    pass_jaccard <-
      is.finite(
        median_jaccard
      ) &&
      median_jaccard >=
      MIN_MEDIAN_JACCARD


    pass_stability <- (

      isTRUE(
        original_row$valid
      ) &&

      pass_valid_rate &&

      pass_anchor_recovery &&

      pass_anchor_iqr &&

      pass_right_fraction_iqr &&

      pass_jaccard
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
        valid_bootstraps,

      bootstrap_valid_rate =
        bootstrap_valid_rate,

      anchor_median =
        anchor_median,

      anchor_iqr =
        anchor_iqr,

      anchor_iqr_fraction =
        anchor_iqr_fraction,

      anchor_tolerance_n =
        anchor_tolerance_n,

      anchor_recovery_rate =
        anchor_recovery_rate,

      right_fraction_median =
        right_fraction_median,

      right_fraction_iqr =
        right_fraction_iqr,

      median_jaccard =
        median_jaccard,

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


  stable_rows <- stability_summary %>%

    filter(
      pass_stability
    ) %>%

    arrange(
      min_fraction
    )


  if (
    nrow(stable_rows) > 0L
  ) {

    selected <- stable_rows[
      1L,
      ,
      drop = FALSE
    ]

  } else {

    selected <- NULL
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
# FEATURE-LEVEL NB2-RELATED METRICS
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


  empirical_var <- row_variance_fast(
    ranked_mat
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


  extra_variance <- pmax(

    empirical_var -
    mu,

    0
  )


  alpha_hat <- rep(
    0,
    length(mu)
  )


  positive_mu <- mu > 0


  alpha_hat[
    positive_mu
  ] <- pmax(

    (
      empirical_var[
        positive_mu
      ] -
      mu[
        positive_mu
      ]
    ) /
    (
      mu[
        positive_mu
      ]^2
    ),

    0
  )


  alpha_mu_value <-
    alpha_hat *
    mu


  data.frame(

    rank =
      seq_along(
        rank_order
      ),

    feature_id =
      rownames(
        metric_mat_arm
      )[
        rank_order
      ],

    mu =
      mu,

    empirical_variance =
      empirical_var,

    NB2 =
      log1p(
        extra_variance
      ),

    NB2_NB1 =
      log1p(
        extra_variance
      ) -
      log1p(
        mu
      ),

    alpha_mu =
      log1p(
        alpha_mu_value
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


  left_end <- anchor_rank -
    1L


  left_start <- left_end -
    right_n +
    1L


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

  if (
    !isTRUE(
      selected_stability_row$pass_stability
    )
  ) {

    stop(
      "Selected fraction did not pass stability."
    )
  }


  rounded_crossings <- as.integer(

    round(
      zero_df$crossing_rank
    )
  )


  if (
    !(
      anchor %in%
      rounded_crossings
    )
  ) {

    stop(
      "Anchor is not a rounded second-derivative zero crossing."
    )
  }


  right_n <- total_n -
    anchor +
    1L


  if (
    (
      right_n /
      total_n
    ) +
    1e-12 <
    selected_fraction
  ) {

    stop(
      "RIGHT is smaller than the selected minimum fraction."
    )
  }


  left_end <- anchor -
    1L


  left_start <- left_end -
    right_n +
    1L


  if (
    left_start < 1L
  ) {

    stop(
      "Matched LEFT block extends below rank 1."
    )
  }


  left_idx <- seq.int(
    left_start,
    left_end
  )


  right_idx <- seq.int(
    anchor,
    total_n
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


  expected <- list(

    left_n =
      nrow(left_df),

    right_n =
      nrow(right_df),

    left_NB2 =
      median(
        left_df$NB2
      ),

    right_NB2 =
      median(
        right_df$NB2
      ),

    diff_NB2 =
      median(
        right_df$NB2
      ) -
      median(
        left_df$NB2
      ),

    left_gap =
      median(
        left_df$NB2_NB1
      ),

    right_gap =
      median(
        right_df$NB2_NB1
      ),

    diff_gap =
      median(
        right_df$NB2_NB1
      ) -
      median(
        left_df$NB2_NB1
      ),

    left_alpha =
      median(
        left_df$alpha_mu
      ),

    right_alpha =
      median(
        right_df$alpha_mu
      ),

    diff_alpha =
      median(
        right_df$alpha_mu
      ) -
      median(
        left_df$alpha_mu
      )
  )


  for (
    nm in names(
      expected
    )
  ) {

    check_value <- all.equal(

      as.numeric(
        region_summary[[nm]][1L]
      ),

      as.numeric(
        expected[[nm]]
      ),

      tolerance =
        1e-10
    )


    if (
      !isTRUE(
        check_value
      )
    ) {

      stop(

        "Validation failed for ",

        nm,

        ": ",

        check_value
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

      layout =
        grid.layout(

          nrow =
            3L,

          ncol =
            1L,

          heights =
            unit(

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

      vp =
        viewport(

          layout.pos.row =
            i,

          layout.pos.col =
            1L
        )
    )
  }


  dev.off()
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


  left_min <-
    region_summary$left_start


  left_max <-
    region_summary$left_end


  right_min <-
    region_summary$right_start


  right_max <-
    region_summary$right_end


  # ---------------------------------------------------------------------------
  # Zero-crossing locations on variance curve
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

  } else {

    zero_plot_df$y <-
      numeric(0)
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


  # ---------------------------------------------------------------------------
  # NB2-related curves
  # ---------------------------------------------------------------------------

  nb_long <- feature_df %>%

    select(
      rank,
      NB2,
      NB2_NB1,
      alpha_mu
    ) %>%

    pivot_longer(

      cols =
        c(
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

      metric =
        factor(

          metric,

          levels =
            c(
              "NB2",
              "NB2_NB1",
              "alpha_mu"
            ),

          labels =
            c(
              "NB2",
              "NB2-NB1",
              "alpha*mu"
            )
        )
    )


  selected_pct <-
    100 *
    selected_stability_row$min_fraction


  actual_pct <-
    100 *
    selected_stability_row$actual_right_fraction


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

      mapping =
        aes(
          rank,
          y
        ),

      inherit.aes =
        FALSE,

      color =
        COL$zero_crossing,

      size =
        1.1,

      alpha =
        0.65
    ) +

    geom_vline(

      xintercept =
        anchor,

      color =
        COL$anchor,

      linewidth =
        0.9
    ) +

    annotate(

      "point",

      x =
        anchor,

      y =
        anchor_y,

      color =
        COL$anchor,

      size =
        3.2
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
        sprintf(

          "%s | %s | %s | Anchor=%d | minimum=%.1f%% | RIGHT=%.1f%% | recovery=%.2f | Jaccard=%.2f",

          fig_tag,

          rank_tag,

          metric_tag,

          anchor,

          selected_pct,

          actual_pct,

          selected_stability_row$anchor_recovery_rate,

          selected_stability_row$median_jaccard
        ),

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
        element_blank()
    )


  # ---------------------------------------------------------------------------
  # PANEL 2: NB2 CORROBORATION
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

      mapping =
        aes(
          rank,
          value,
          color = metric
        ),

      linewidth =
        0.95
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
        sprintf(

          "Matched LEFT/RIGHT | RIGHT-LEFT: NB2=%.3f, NB2-NB1=%.3f, alpha*mu=%.3f",

          region_summary$diff_NB2,

          region_summary$diff_gap,

          region_summary$diff_alpha
        ),

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
  # PANEL 3: MATCHED REGION SUMMARY
  # ---------------------------------------------------------------------------

  summary_df <- data.frame(

    metric =
      factor(

        c(
          "NB2",
          "NB2-NB1",
          "alpha*mu"
        ),

        levels =
          rev(

            c(
              "NB2",
              "NB2-NB1",
              "alpha*mu"
            )
          )
      ),

    LEFT =
      c(

        region_summary$left_NB2,

        region_summary$left_gap,

        region_summary$left_alpha
      ),

    RIGHT =
      c(

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
          ": matched-region summary"
        ),

      subtitle =
        paste0(

          "LEFT n = ",

          region_summary$left_n,

          "; RIGHT n = ",

          region_summary$right_n
        ),

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


  prefix <- paste(

    comparison_name,

    arm_name,

    track,

    sep = "_"
  )


  message(

    "[",

    track,

    "] ",

    comparison_name,

    " ",

    arm_name,

    ": bootstrapping ",

    BOOTSTRAP_N,

    " replicates..."
  )


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
      seed,

    cores =
      BOOTSTRAP_CORES
  )


  # ---------------------------------------------------------------------------
  # Stability outputs
  # ---------------------------------------------------------------------------

  write.csv(

    stability_obj$stability_summary,

    file.path(

      output_dir,

      paste0(
        "Table_StabilitySummary_",
        prefix,
        ".csv"
      )
    ),

    row.names = FALSE
  )


  write.csv(

    stability_obj$bootstrap_df,

    file.path(

      output_dir,

      paste0(
        "Table_StabilityBootstrap_",
        prefix,
        ".csv"
      )
    ),

    row.names = FALSE
  )


  write.csv(

    stability_obj$original_candidates,

    file.path(

      output_dir,

      paste0(
        "Table_CandidateAnchors_",
        prefix,
        ".csv"
      )
    ),

    row.names = FALSE
  )


  write.csv(

    stability_obj$original_geometry$zero_df,

    file.path(

      output_dir,

      paste0(
        "Table_Zero_",
        prefix,
        ".csv"
      )
    ),

    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # No stable cutoff
  # ---------------------------------------------------------------------------

  if (
    is.null(
      stability_obj$selected
    )
  ) {

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

      ActualRightFraction =
        NA_real_,

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

      Call_NB2 =
        NA_character_,

      Call_Gap =
        NA_character_,

      Call_Alpha =
        NA_character_,

      stringsAsFactors = FALSE
    )


    write.csv(

      fail_summary,

      file.path(

        output_dir,

        paste0(
          "Table_Cutoff_",
          prefix,
          ".csv"
        )
      ),

      row.names = FALSE
    )


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
          paste0(

            "No candidate minimum fraction passed all ",

            "bootstrap stability criteria."
          ),

        stringsAsFactors = FALSE
      ),

      file.path(

        output_dir,

        paste0(
          "Table_Valid_",
          prefix,
          ".csv"
        )
      ),

      row.names = FALSE
    )


    message(

      "[",

      track,

      "] ",

      comparison_name,

      " ",

      arm_name,

      " | UNSTABLE"
    )


    return(

      list(

        summary =
          fail_summary,

        stability_summary =
          stability_obj$stability_summary
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Stable solution
  # ---------------------------------------------------------------------------

  selected <- stability_obj$selected


  anchor <- as.integer(
    selected$anchor
  )


  selected_fraction <- as.numeric(
    selected$min_fraction
  )


  geometry <-
    stability_obj$original_geometry


  feature_df <- compute_ranked_feature_metrics(

    metric_mat_arm =
      metric_matrix,

    rank_order =
      geometry$rank_order
  )


  feature_df$abs_pc1_loading <-
    geometry$abs_loadings[
      geometry$rank_order
    ]


  region_summary <- summarize_regions(

    feature_df =
      feature_df,

    anchor_rank =
      anchor,

    total_n =
      total_n
  )


  validate_method_level(

    feature_df =
      feature_df,

    zero_df =
      geometry$zero_df,

    anchor =
      anchor,

    selected_fraction =
      selected_fraction,

    selected_stability_row =
      selected,

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
      selected$actual_right_fraction,

    BootstrapValidRate =
      selected$bootstrap_valid_rate,

    AnchorRecoveryRate =
      selected$anchor_recovery_rate,

    AnchorIQR =
      selected$anchor_iqr,

    AnchorIQRFraction =
      selected$anchor_iqr_fraction,

    RightFractionIQR =
      selected$right_fraction_iqr,

    MedianJaccard =
      selected$median_jaccard,

    stringsAsFactors = FALSE
  )


  cutoff_summary <- bind_cols(

    selected_df,

    region_summary
  ) %>%

    mutate(

      LeftSize =
        left_n,

      RightSize =
        right_n,

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
  # Zero crossing table with selected Anchor marked
  # ---------------------------------------------------------------------------

  zero_out <- geometry$zero_df


  zero_out$crossing_rank_rounded <- as.integer(

    round(
      zero_out$crossing_rank
    )
  )


  zero_out$selected_anchor <- (
    zero_out$crossing_rank_rounded ==
    anchor
  )


  write.csv(

    zero_out,

    file.path(

      output_dir,

      paste0(
        "Table_Zero_",
        prefix,
        ".csv"
      )
    ),

    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Cutoff output
  # ---------------------------------------------------------------------------

  write.csv(

    cutoff_summary,

    file.path(

      output_dir,

      paste0(
        "Table_Cutoff_",
        prefix,
        ".csv"
      )
    ),

    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Rank-level output
  # ---------------------------------------------------------------------------

  rank_output <- feature_df %>%

    left_join(

      geometry$variance_df,

      by =
        "rank",

      suffix =
        c(
          "_metric",
          "_geometry"
        )
    )


  write.csv(

    rank_output,

    file.path(

      output_dir,

      paste0(
        "Table_Rank_",
        prefix,
        ".csv"
      )
    ),

    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Validation output
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
        selected$actual_right_fraction,

      BootstrapValidRate =
        selected$bootstrap_valid_rate,

      AnchorRecoveryRate =
        selected$anchor_recovery_rate,

      AnchorIQRFraction =
        selected$anchor_iqr_fraction,

      RightFractionIQR =
        selected$right_fraction_iqr,

      MedianJaccard =
        selected$median_jaccard,

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
            geometry$zero_df$crossing_rank
          )
        ),

      stringsAsFactors = FALSE
    ),

    file.path(

      output_dir,

      paste0(
        "Table_Valid_",
        prefix,
        ".csv"
      )
    ),

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
      geometry$variance_df,

    zero_df =
      geometry$zero_df,

    feature_df =
      feature_df,

    anchor =
      anchor,

    selected_stability_row =
      selected,

    region_summary =
      region_summary,

    out_file =
      file.path(

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
      ),

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

    " | fraction=",

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
      selected$bootstrap_valid_rate
    ),

    " | recovery=",

    sprintf(
      "%.3f",
      selected$anchor_recovery_rate
    ),

    " | IQR/N=",

    sprintf(
      "%.4f",
      selected$anchor_iqr_fraction
    ),

    " | Jaccard=",

    sprintf(
      "%.3f",
      selected$median_jaccard
    )
  )


  list(

    summary =
      cutoff_summary,

    stability_summary =
      stability_obj$stability_summary
  )
}


# =============================================================================
# RUN ANALYSIS
# =============================================================================

count_mat <- read_count_matrix(

  COUNT_FILE,

  COMPARISONS
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

  "Bootstrap cores: ",

  BOOTSTRAP_CORES
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


  pats <- COMPARISONS[[
    comparison_name
  ]]


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
    # MAIN MANUSCRIPT TRACK
    # =========================================================================

    main_rank_matrix <- normalize_cpm_log1p(
      count_mat_arm
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
        seed_from_label(

          paste(

            comparison_name,

            arm_name,

            "Main",

            sep = "_"
          )
        )
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

        comp =
          comparison_name,

        arm =
          arm_name,

        track =
          "Main",

        .before =
          1L
      )


    # =========================================================================
    # DESEQ2 SUPPLEMENTARY TRACK
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

          "[DESeq2] skipped: required packages unavailable for ",

          comparison_name,

          " ",

          arm_name
        )

      } else {

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
            paste0(

              "Rank: ",

              deseq2_obj$rank_method_used
            ),

          metric_tag =
            "Metrics: DESeq2 normalized",

          metric_name =
            "deseq2_normalized_counts",

          seed =
            seed_from_label(

              paste(

                comparison_name,

                arm_name,

                "DESeq2",

                sep = "_"
              )
            )
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


        overall_rows[[
          length(overall_rows) +
          1L
        ]] <- deseq2_res$summary


        overall_stability_rows[[
          length(overall_stability_rows) +
          1L
        ]] <- deseq2_res$stability_summary %>%

          mutate(

            comp =
              comparison_name,

            arm =
              arm_name,

            track =
              "DESeq2",

            .before =
              1L
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


overall_stability <- bind_rows(
  overall_stability_rows
)


write.csv(

  overall_summary,

  file.path(
    OUT_ROOT,
    "Table_Overall_Cutoff.csv"
  ),

  row.names = FALSE
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
  "============================================================"
)


message(
  "Analysis complete."
)


message(
  "No fixed 5,000-feature reference was used."
)


message(
  "Cutoff selection used PC1/variance geometry and bootstrap stability only."
)


message(
  "NB2-related quantities were evaluated only after cutoff selection."
)


message(
  "Outputs written to: ",
  OUT_ROOT
)


message(
  "============================================================"
)
