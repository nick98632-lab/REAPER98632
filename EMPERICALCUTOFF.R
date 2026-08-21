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
# DATA-DERIVED PC1 LEADING EDGE WITH BOOTSTRAP-STABLE GEOMETRIC BOUNDARIES
# =============================================================================
#
# METHODS
#
# Within each experimental arm, features are ranked from lowest to highest
# absolute loading on the first principal component (|PC1 loading|). The right
# edge of the ranked series therefore contains the features contributing most
# strongly to PC1 in absolute loading magnitude.
#
# The previous fixed 5,000-feature reference is removed completely. Instead,
# the characteristic size of the PC1 leading edge is estimated directly from
# the concentration of the PC1 loading vector using the inverse participation
# ratio (effective number of contributing features):
#
#     p_i = loading_i^2 / sum_j(loading_j^2)
#
#     N_eff = 1 / sum_i(p_i^2)
#
# N_eff is the number of equally weighted features that would produce the same
# concentration of squared PC1 loading mass. The data-derived reference rank is
# then defined so that N_eff features lie at or to the right of that rank:
#
#     Ref = N - round(N_eff) + 1
#
# Thus, if the intrinsic PC1 loading concentration corresponds to approximately
# 5,000 features, a ~5,000-feature reference emerges from the data without
# explicitly specifying 5,000.
#
# To prevent a trivial terminal fluctuation involving only a very small number
# of features from being interpreted as a leading edge, the effective PC1 tail
# must contain at least 5% of all ranked features. This 5% value is an
# admissibility guardrail only; it does not select or optimize the cutoff.
#
# Empirical feature variance is then ordered along the PC1 rank axis,
# transformed as log(1 + variance), and represented by a smoothing spline. The
# second derivative of the spline is evaluated on a dense rank grid. Sign
# changes in the second derivative identify curvature transitions.
#
# The geometric interval is defined by the nearest curvature transitions that
# flank the data-derived PC1 reference:
#
#     Anchor   = nearest d2 zero crossing immediately LEFT of Ref
#     Terminal = nearest d2 zero crossing immediately RIGHT of Ref
#
# The downstream leading-edge region is:
#
#     RIGHT = Anchor through rank N
#
# and is compared with an immediately adjacent equal-sized block:
#
#     LEFT = matched block immediately left of Anchor
#
# The matched comparison is allowed only when the Anchor leaves enough features
# on the left to construct an equal-sized LEFT block.
#
# -----------------------------------------------------------------------------
# BOOTSTRAP STABILITY
# -----------------------------------------------------------------------------
#
# Samples are resampled with replacement within each experimental arm. For each
# bootstrap replicate, PC1, N_eff, Ref, the ranked variance curve, spline,
# second derivative, Anchor, and Terminal are recomputed.
#
# A track is considered geometrically stable only if all prespecified criteria
# are met:
#
#   - >= 90% of bootstrap replicates yield a valid flanking interval;
#   - >= 80% of valid replicates recover Ref within 3% of the rank axis;
#   - >= 80% of valid replicates recover Anchor within 3% of the rank axis;
#   - >= 80% of valid replicates recover Terminal within 3% of the rank axis;
#   - Ref IQR <= 3% of the rank axis;
#   - Anchor IQR <= 3% of the rank axis;
#   - Terminal IQR <= 3% of the rank axis;
#   - RIGHT-region fraction IQR <= 3 percentage points.
#
# Jaccard overlap of RIGHT-feature membership is reported descriptively but is
# NOT used as a pass/fail criterion. This is deliberate: the inferential target
# is stability of the leading-edge boundary in rank space, not exact identity
# of every feature after PC1 is re-estimated from a small resampled sample set.
#
# -----------------------------------------------------------------------------
# NB2-RELATED CORROBORATION
# -----------------------------------------------------------------------------
#
# NB2-related quantities are calculated only after the geometric cutoff passes
# bootstrap stability. They do not participate in cutoff selection.
#
# Let mu denote the empirical feature mean and variance the empirical variance:
#
#     NB2 = log(1 + max(variance - mu, 0))
#
#     NB2-NB1 = log(1 + max(variance - mu, 0)) - log(1 + mu)
#
#     alpha = max((variance - mu) / mu^2, 0)
#
#     alpha*mu = log(1 + alpha*mu)
#
# These are descriptive excess-variance diagnostics rather than formal
# likelihood-ratio statistics.
#
# If a track fails the prespecified geometric stability criteria, it is written
# as UNSTABLE and no NB2-based cutoff interpretation is forced.
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <- "/root/REAPER98632/exports/manuscript_effective_pc1_bootstrap"

VAR_SPLINE_SPAR <- 0.60

DENSE_GRID_MULTIPLIER <- 2L
DENSE_GRID_MIN <- 5000L

# Guardrail only; this does NOT define the selected leading edge.
MIN_EFFECTIVE_TAIL_FRACTION <- 0.05

BOOTSTRAP_N <- 500L
BOOTSTRAP_SEED_BASE <- 20260820L

MIN_BOOTSTRAP_VALID_RATE <- 0.90
POSITION_TOLERANCE_FRACTION <- 0.03
MIN_POSITION_RECOVERY_RATE <- 0.80
MAX_POSITION_IQR_FRACTION <- 0.03
MAX_RIGHT_FRACTION_IQR <- 0.03

PNG_WIDTH_IN <- 14
PNG_HEIGHT_IN <- 10.8
PNG_DPI <- 260

RUN_DESEQ2_SUPPLEMENT <- TRUE

DESEQ2_RANK_METHOD <- "normalized_log1p"
# Options:
#   "normalized_log1p"
#   "vst"

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

DETECTED_CORES <- suppressWarnings(
  parallel::detectCores(logical = FALSE)
)

if (is.na(DETECTED_CORES) || DETECTED_CORES < 2L) {
  BOOTSTRAP_CORES <- 1L
} else {
  BOOTSTRAP_CORES <- min(
    4L,
    DETECTED_CORES - 1L
  )
}


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
  interval_fill = "#9E9E9E",

  anchor = "#000000",
  ref = "#E69F00",
  terminal = "#D95F02",

  left_pt = "#5B8FD1",
  right_pt = "#43A047"
)

EVENT_LEVELS <- c(
  "Anchor",
  "Ref",
  "Terminal"
)

EVENT_COLORS <- c(
  "Anchor" = COL$anchor,
  "Ref" = COL$ref,
  "Terminal" = COL$terminal
)

EVENT_SHAPES <- c(
  "Anchor" = 16,
  "Ref" = 18,
  "Terminal" = 1
)

EVENT_LTY <- c(
  "Anchor" = "solid",
  "Ref" = "dashed",
  "Terminal" = "dotted"
)

TRACE_LEVELS <- c(
  "NB2",
  "NB2-NB1",
  "alpha*mu"
)

TRACE_COLORS <- c(
  "NB2" = COL$nb2,
  "NB2-NB1" = COL$nb_gap,
  "alpha*mu" = COL$alpha_mu
)

REGION_LEVELS <- c(
  "LEFT",
  "RIGHT"
)

REGION_COLORS <- c(
  "LEFT" = COL$left_pt,
  "RIGHT" = COL$right_pt
)


# =============================================================================
# INPUT
# =============================================================================

read_count_matrix <- function(path, comparisons) {

  if (!file.exists(path)) {
    stop(
      "Count file does not exist: ",
      path
    )
  }

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
  # Use only columns belonging to the prespecified experimental arms.
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

  if (length(sample_idx) == 0L) {
    stop(
      "No sample columns matched COMPARISONS."
    )
  }

  if (1L %in% sample_idx) {
    stop(
      "Column 1 matched a sample pattern; ",
      "column 1 must contain feature identifiers."
    )
  }


  # ---------------------------------------------------------------------------
  # Parse count columns.
  #
  # This preserves the behavior of the original working analysis:
  # malformed/non-finite cells are replaced by zero rather than terminating
  # the complete pipeline.
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

  bad_n <- sum(
    !is.finite(count_mat)
  )

  if (bad_n > 0L) {

    message(
      "Replacing ",
      bad_n,
      " non-finite count entries with 0 ",
      "(original-script behavior)."
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
  # Feature IDs
  #
  # The previous failure occurred because one feature identifier was blank.
  # Blank IDs receive deterministic row-based IDs instead of stopping.
  # ---------------------------------------------------------------------------

  feature_ids <- trimws(
    as.character(
      raw_df[[1L]]
    )
  )

  blank <- (
    is.na(feature_ids) |
    feature_ids == ""
  )

  if (any(blank)) {

    feature_ids[
      blank
    ] <- paste0(
      "__feature_row_",
      which(blank)
    )

    message(
      "Assigned deterministic row IDs to ",
      sum(blank),
      " blank feature identifier(s)."
    )
  }

  feature_ids <- make.unique(
    feature_ids,
    sep = "__dup_"
  )

  rownames(count_mat) <- feature_ids


  # ---------------------------------------------------------------------------
  # Remove all-zero features.
  # ---------------------------------------------------------------------------

  keep <- rowSums(
    count_mat
  ) > 0

  count_mat <- count_mat[
    keep,
    ,
    drop = FALSE
  ]

  if (nrow(count_mat) < 2L) {
    stop(
      "Fewer than two nonzero features remain after filtering."
    )
  }

  count_mat
}


# =============================================================================
# NORMALIZATION
# =============================================================================

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
    2L,
    lib_sizes / 1e6,
    "/"
  )

  log1p(cpm)
}


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
      "DESEQ2_RANK_METHOD must be ",
      "'normalized_log1p' or 'vst'."
    )
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

  dds <- DESeq2::estimateSizeFactors(
    dds
  )

  norm_counts <- DESeq2::counts(
    dds,
    normalized = TRUE
  )

  if (rank_method == "vst") {

    vst_obj <- tryCatch(
      DESeq2::vst(
        dds,
        blind = TRUE
      ),
      error = function(e) {
        NULL
      }
    )

    if (is.null(vst_obj)) {

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
    normalized_counts = norm_counts,

    ranking_matrix = ranking_matrix,

    size_factors =
      DESeq2::sizeFactors(
        dds
      ),

    rank_method_used =
      rank_method_used
  )
}


# =============================================================================
# NUMERICAL HELPERS
# =============================================================================

row_variance_fast <- function(mat) {

  n <- ncol(
    mat
  )

  if (n < 2L) {
    stop(
      "At least two samples are required ",
      "for variance estimation."
    )
  }

  mu <- rowMeans(
    mat
  )

  ss <- rowSums(
    mat * mat
  ) -
    n * mu^2

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

compute_abs_pc1_loadings <- function(norm_mat_arm) {

  if (
    nrow(norm_mat_arm) < 2L ||
    ncol(norm_mat_arm) < 2L
  ) {
    stop(
      "PC1 requires at least two features ",
      "and two samples."
    )
  }

  # Rows = samples
  # Columns = features

  X <- t(
    norm_mat_arm
  )

  X_centered <- sweep(
    X,
    2L,
    colMeans(X),
    "-"
  )

  # Sample-space PCA.
  #
  # This yields the same leading loading vector, up to sign, as the
  # corresponding SVD/prcomp solution while being efficient when the number
  # of samples is much smaller than the number of features.

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
    lambda1 <= .Machine$double.eps
  ) {
    stop(
      "PC1 is undefined because between-sample ",
      "variation is insufficient."
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
# EFFECTIVE PC1 LEADING-EDGE SIZE
# =============================================================================

compute_effective_pc1_reference <- function(
    abs_loadings,
    total_n,
    min_tail_fraction = MIN_EFFECTIVE_TAIL_FRACTION) {

  squared_loading <- abs_loadings^2

  loading_mass <- sum(
    squared_loading
  )

  if (
    !is.finite(loading_mass) ||
    loading_mass <= 0
  ) {

    return(
      list(
        valid = FALSE,
        reason = "invalid_pc1_loading_mass"
      )
    )
  }

  p <- squared_loading /
    loading_mass

  participation_denominator <- sum(
    p^2
  )

  if (
    !is.finite(
      participation_denominator
    ) ||
    participation_denominator <= 0
  ) {

    return(
      list(
        valid = FALSE,
        reason = "invalid_effective_feature_number"
      )
    )
  }

  n_eff <- 1 /
    participation_denominator

  n_eff_n <- as.integer(
    round(
      n_eff
    )
  )

  minimum_tail_n <- ceiling(
    min_tail_fraction *
      total_n
  )


  # ---------------------------------------------------------------------------
  # Guardrail against microscopic terminal edges.
  # ---------------------------------------------------------------------------

  if (
    n_eff_n <
      minimum_tail_n
  ) {

    return(
      list(
        valid = FALSE,
        reason =
          "effective_pc1_tail_below_guardrail",
        n_eff = n_eff,
        n_eff_n = n_eff_n,
        minimum_tail_n = minimum_tail_n
      )
    )
  }


  if (
    n_eff_n >=
      total_n
  ) {

    return(
      list(
        valid = FALSE,
        reason =
          "effective_pc1_tail_spans_entire_rank_axis",
        n_eff = n_eff,
        n_eff_n = n_eff_n,
        minimum_tail_n = minimum_tail_n
      )
    )
  }


  reference_rank <- total_n -
    n_eff_n +
    1L


  list(
    valid = TRUE,
    reason = NA_character_,

    n_eff = n_eff,

    n_eff_n = n_eff_n,

    n_eff_fraction =
      n_eff_n /
      total_n,

    ref =
      as.integer(
        reference_rank
      ),

    minimum_tail_n =
      minimum_tail_n
  )
}


# =============================================================================
# VARIANCE GEOMETRY
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
    DENSE_GRID_MIN,
    length(ranks) *
      DENSE_GRID_MULTIPLIER
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
    dense_rank = dense_x,

    dense_smooth_log1p_empirical_variance =
      dense_y,

    dense_d2 =
      dense_d2,

    stringsAsFactors = FALSE
  )

  out
}


find_d2_zero_crossings <- function(dense_df) {

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

  if (length(x) < 2L) {

    return(
      data.frame(
        crossing_rank = numeric(0),
        crossing_type = character(0),
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
    length(
      sign_change_idx
    ) > 0L
  ) {

    fraction_between_points <-
      abs(
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

    crossings <-
      xa[
        sign_change_idx
      ] +
      fraction_between_points *
      (
        xb[
          sign_change_idx
        ] -
        xa[
          sign_change_idx
        ]
      )
  }


  # Retain exact zero locations as well.

  exact_zero_idx <- which(
    y == 0
  )

  if (
    length(
      exact_zero_idx
    ) > 0L
  ) {

    crossings <- c(
      crossings,
      x[
        exact_zero_idx
      ]
    )
  }

  crossings <- crossings[
    is.finite(
      crossings
    )
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
    length(
      crossings
    ) == 0L
  ) {

    return(
      data.frame(
        crossing_rank = numeric(0),
        crossing_type = character(0),
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
# SELECT GEOMETRIC INTERVAL AROUND DATA-DERIVED REFERENCE
# =============================================================================

select_flanking_interval <- function(
    zero_df,
    ref,
    total_n) {

  invalid_result <- function(reason) {

    list(
      valid = FALSE,

      reason = reason,

      anchor = NA_integer_,

      ref = as.integer(
        ref
      ),

      terminal = NA_integer_,

      right_n = NA_integer_,

      right_fraction = NA_real_
    )
  }


  if (
    nrow(
      zero_df
    ) == 0L
  ) {

    return(
      invalid_result(
        "no_d2_zero_crossings"
      )
    )
  }


  left_candidates <-
    zero_df$crossing_rank[
      zero_df$crossing_rank <
      ref
    ]

  right_candidates <-
    zero_df$crossing_rank[
      zero_df$crossing_rank >
      ref
    ]


  if (
    length(
      left_candidates
    ) == 0L
  ) {

    return(
      invalid_result(
        "no_left_flanking_crossing"
      )
    )
  }


  if (
    length(
      right_candidates
    ) == 0L
  ) {

    return(
      invalid_result(
        "no_right_flanking_crossing"
      )
    )
  }


  anchor <- as.integer(
    round(
      max(
        left_candidates
      )
    )
  )

  terminal <- as.integer(
    round(
      min(
        right_candidates
      )
    )
  )

  ref <- as.integer(
    round(
      ref
    )
  )


  anchor <- max(
    1L,
    anchor
  )

  terminal <- min(
    total_n,
    terminal
  )


  if (
    !(
      anchor <
      ref &&
      ref <
      terminal
    )
  ) {

    return(
      invalid_result(
        "invalid_anchor_ref_terminal_order"
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Ensure equal-sized matched LEFT region can be constructed.
  # ---------------------------------------------------------------------------

  right_n <- total_n -
    anchor +
    1L

  left_start <- anchor -
    right_n


  if (
    left_start < 1L
  ) {

    return(
      invalid_result(
        "matched_left_block_not_possible"
      )
    )
  }


  list(
    valid = TRUE,

    reason = NA_character_,

    anchor = anchor,

    ref = ref,

    terminal = terminal,

    interval_width =
      terminal -
      anchor +
      1L,

    right_n =
      right_n,

    right_fraction =
      right_n /
      total_n,

    left_start =
      left_start,

    left_end =
      anchor -
      1L
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
      "Rank and metric matrices have ",
      "different feature counts."
    )
  }


  if (
    ncol(rank_matrix) !=
    ncol(metric_matrix)
  ) {

    stop(
      "Rank and metric matrices have ",
      "different sample counts."
    )
  }


  if (
    !identical(
      rownames(rank_matrix),
      rownames(metric_matrix)
    )
  ) {

    stop(
      "Feature IDs differ between ",
      "rank and metric matrices."
    )
  }


  total_n <- nrow(
    rank_matrix
  )


  abs_loadings <- compute_abs_pc1_loadings(
    rank_matrix
  )


  reference_info <-
    compute_effective_pc1_reference(
      abs_loadings,
      total_n
    )


  if (
    !isTRUE(
      reference_info$valid
    )
  ) {

    return(
      list(
        valid = FALSE,

        reason =
          reference_info$reason,

        reference_info =
          reference_info
      )
    )
  }


  rank_order <- order(
    abs_loadings,
    decreasing = FALSE
  )


  variance_df <-
    compute_ranked_variance_curve(
      metric_matrix,
      rank_order,
      spar
    )


  dense_curve_df <- attr(
    variance_df,
    "dense_curve_df"
  )


  zero_df <-
    find_d2_zero_crossings(
      dense_curve_df
    )


  interval <-
    select_flanking_interval(
      zero_df,
      reference_info$ref,
      total_n
    )


  if (
    !isTRUE(
      interval$valid
    )
  ) {

    return(
      list(
        valid = FALSE,

        reason =
          interval$reason,

        abs_loadings =
          abs_loadings,

        reference_info =
          reference_info,

        rank_order =
          rank_order,

        variance_df =
          variance_df,

        dense_curve_df =
          dense_curve_df,

        zero_df =
          zero_df,

        interval =
          interval
      )
    )
  }


  list(
    valid = TRUE,

    reason = NA_character_,

    abs_loadings =
      abs_loadings,

    reference_info =
      reference_info,

    rank_order =
      rank_order,

    variance_df =
      variance_df,

    dense_curve_df =
      dense_curve_df,

    zero_df =
      zero_df,

    interval =
      interval
  )
}


# =============================================================================
# BOOTSTRAP STABILITY
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
    length(
      union_set
    ) == 0L
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


assess_geometry_stability <- function(
    rank_matrix,
    metric_matrix,
    bootstrap_n = BOOTSTRAP_N,
    spar = VAR_SPLINE_SPAR,
    seed = BOOTSTRAP_SEED_BASE,
    cores = BOOTSTRAP_CORES) {

  total_n <- nrow(
    rank_matrix
  )

  sample_n <- ncol(
    rank_matrix
  )


  # ---------------------------------------------------------------------------
  # Full-data geometry
  # ---------------------------------------------------------------------------

  full_geometry <- compute_geometry(
    rank_matrix,
    metric_matrix,
    spar
  )


  if (
    !isTRUE(
      full_geometry$valid
    )
  ) {

    return(
      list(
        full_geometry =
          full_geometry,

        bootstrap_df =
          data.frame(),

        stability_summary =
          data.frame(
            Status =
              "UNSTABLE",

            Reason =
              paste0(
                "Full-data geometry invalid: ",
                full_geometry$reason
              ),

            PassStability =
              FALSE,

            stringsAsFactors = FALSE
          )
      )
    )
  }


  full_interval <-
    full_geometry$interval

  full_ref_info <-
    full_geometry$reference_info


  full_right_idx <-
    full_geometry$rank_order[
      seq.int(
        full_interval$anchor,
        total_n
      )
    ]


  full_right_features <-
    rownames(
      rank_matrix
    )[
      full_right_idx
    ]


  # ---------------------------------------------------------------------------
  # Generate bootstrap draws before parallel processing for reproducibility.
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
        size = sample_n,
        replace = TRUE
      )
    }
  )


  # ---------------------------------------------------------------------------
  # Single bootstrap replicate
  # ---------------------------------------------------------------------------

  bootstrap_worker <- function(b) {

    idx <- bootstrap_indices[[b]]


    invalid_row <- function(reason) {

      data.frame(
        bootstrap = b,

        valid = FALSE,

        invalid_reason =
          reason,

        effective_n =
          NA_real_,

        effective_n_rounded =
          NA_integer_,

        effective_fraction =
          NA_real_,

        ref =
          NA_integer_,

        anchor =
          NA_integer_,

        terminal =
          NA_integer_,

        interval_width =
          NA_integer_,

        right_n =
          NA_integer_,

        right_fraction =
          NA_real_,

        jaccard =
          NA_real_,

        stringsAsFactors = FALSE
      )
    }


    if (
      length(
        unique(
          idx
        )
      ) < 2L
    ) {

      return(
        invalid_row(
          "fewer_than_2_unique_samples"
        )
      )
    }


    boot_rank <- rank_matrix[
      ,
      idx,
      drop = FALSE
    ]


    boot_metric <- metric_matrix[
      ,
      idx,
      drop = FALSE
    ]


    g <- tryCatch(
      compute_geometry(
        boot_rank,
        boot_metric,
        spar
      ),
      error = function(e) {
        NULL
      }
    )


    if (
      is.null(g)
    ) {

      return(
        invalid_row(
          "geometry_error"
        )
      )
    }


    if (
      !isTRUE(
        g$valid
      )
    ) {

      return(
        invalid_row(
          g$reason
        )
      )
    }


    right_idx <-
      g$rank_order[
        seq.int(
          g$interval$anchor,
          total_n
        )
      ]


    right_features <-
      rownames(
        rank_matrix
      )[
        right_idx
      ]


    data.frame(
      bootstrap = b,

      valid = TRUE,

      invalid_reason =
        NA_character_,

      effective_n =
        g$reference_info$n_eff,

      effective_n_rounded =
        g$reference_info$n_eff_n,

      effective_fraction =
        g$reference_info$n_eff_fraction,

      ref =
        g$interval$ref,

      anchor =
        g$interval$anchor,

      terminal =
        g$interval$terminal,

      interval_width =
        g$interval$interval_width,

      right_n =
        g$interval$right_n,

      right_fraction =
        g$interval$right_fraction,

      jaccard =
        jaccard_similarity(
          full_right_features,
          right_features
        ),

      stringsAsFactors = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Run bootstrap
  # ---------------------------------------------------------------------------

  ids <- seq_len(
    bootstrap_n
  )


  if (
    .Platform$OS.type == "unix" &&
    cores > 1L
  ) {

    boot_list <- parallel::mclapply(
      ids,
      bootstrap_worker,
      mc.cores = cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )

  } else {

    boot_list <- lapply(
      ids,
      bootstrap_worker
    )
  }


  bootstrap_df <- bind_rows(
    boot_list
  )


  valid_df <- bootstrap_df[
    bootstrap_df$valid,
    ,
    drop = FALSE
  ]


  valid_n <- nrow(
    valid_df
  )


  valid_rate <- valid_n /
    bootstrap_n


  tolerance_n <- ceiling(
    POSITION_TOLERANCE_FRACTION *
      total_n
  )


  safe_median <- function(x) {

    x <- x[
      is.finite(x)
    ]

    if (length(x) == 0L) {
      return(
        NA_real_
      )
    }

    median(x)
  }


  safe_iqr <- function(x) {

    x <- x[
      is.finite(x)
    ]

    if (length(x) == 0L) {
      return(
        NA_real_
      )
    }

    stats::IQR(x)
  }


  if (
    valid_n > 0L
  ) {

    ref_recovery <- mean(
      abs(
        valid_df$ref -
        full_interval$ref
      ) <=
        tolerance_n
    )


    anchor_recovery <- mean(
      abs(
        valid_df$anchor -
        full_interval$anchor
      ) <=
        tolerance_n
    )


    terminal_recovery <- mean(
      abs(
        valid_df$terminal -
        full_interval$terminal
      ) <=
        tolerance_n
    )


    ref_iqr <- safe_iqr(
      valid_df$ref
    )


    anchor_iqr <- safe_iqr(
      valid_df$anchor
    )


    terminal_iqr <- safe_iqr(
      valid_df$terminal
    )


    right_fraction_iqr <- safe_iqr(
      valid_df$right_fraction
    )


    effective_fraction_iqr <- safe_iqr(
      valid_df$effective_fraction
    )


    median_jaccard <- safe_median(
      valid_df$jaccard
    )

  } else {

    ref_recovery <- NA_real_

    anchor_recovery <- NA_real_

    terminal_recovery <- NA_real_

    ref_iqr <- NA_real_

    anchor_iqr <- NA_real_

    terminal_iqr <- NA_real_

    right_fraction_iqr <- NA_real_

    effective_fraction_iqr <- NA_real_

    median_jaccard <- NA_real_
  }


  # ---------------------------------------------------------------------------
  # Prespecified geometric stability criteria
  # ---------------------------------------------------------------------------

  pass_valid <-
    is.finite(
      valid_rate
    ) &&
    valid_rate >=
      MIN_BOOTSTRAP_VALID_RATE


  pass_ref_recovery <-
    is.finite(
      ref_recovery
    ) &&
    ref_recovery >=
      MIN_POSITION_RECOVERY_RATE


  pass_anchor_recovery <-
    is.finite(
      anchor_recovery
    ) &&
    anchor_recovery >=
      MIN_POSITION_RECOVERY_RATE


  pass_terminal_recovery <-
    is.finite(
      terminal_recovery
    ) &&
    terminal_recovery >=
      MIN_POSITION_RECOVERY_RATE


  pass_ref_iqr <-
    is.finite(
      ref_iqr
    ) &&
    (
      ref_iqr /
      total_n
    ) <=
      MAX_POSITION_IQR_FRACTION


  pass_anchor_iqr <-
    is.finite(
      anchor_iqr
    ) &&
    (
      anchor_iqr /
      total_n
    ) <=
      MAX_POSITION_IQR_FRACTION


  pass_terminal_iqr <-
    is.finite(
      terminal_iqr
    ) &&
    (
      terminal_iqr /
      total_n
    ) <=
      MAX_POSITION_IQR_FRACTION


  pass_right_fraction_iqr <-
    is.finite(
      right_fraction_iqr
    ) &&
    right_fraction_iqr <=
      MAX_RIGHT_FRACTION_IQR


  # Jaccard is deliberately NOT included here.

  pass_stability <- all(
    c(
      pass_valid,
      pass_ref_recovery,
      pass_anchor_recovery,
      pass_terminal_recovery,
      pass_ref_iqr,
      pass_anchor_iqr,
      pass_terminal_iqr,
      pass_right_fraction_iqr
    )
  )


  stability_summary <- data.frame(
    Status =
      ifelse(
        pass_stability,
        "PASS",
        "UNSTABLE"
      ),

    Reason =
      ifelse(
        pass_stability,
        NA_character_,
        "One or more prespecified geometric stability criteria failed"
      ),

    FullEffectiveN =
      full_ref_info$n_eff,

    FullEffectiveNRounded =
      full_ref_info$n_eff_n,

    FullEffectiveFraction =
      full_ref_info$n_eff_fraction,

    FullRef =
      full_interval$ref,

    FullAnchor =
      full_interval$anchor,

    FullTerminal =
      full_interval$terminal,

    FullIntervalWidth =
      full_interval$interval_width,

    FullRightN =
      full_interval$right_n,

    FullRightFraction =
      full_interval$right_fraction,

    BootstrapN =
      bootstrap_n,

    ValidBootstraps =
      valid_n,

    BootstrapValidRate =
      valid_rate,

    RefMedian =
      safe_median(
        valid_df$ref
      ),

    RefIQR =
      ref_iqr,

    RefIQRFraction =
      ref_iqr /
      total_n,

    RefRecoveryRate =
      ref_recovery,

    AnchorMedian =
      safe_median(
        valid_df$anchor
      ),

    AnchorIQR =
      anchor_iqr,

    AnchorIQRFraction =
      anchor_iqr /
      total_n,

    AnchorRecoveryRate =
      anchor_recovery,

    TerminalMedian =
      safe_median(
        valid_df$terminal
      ),

    TerminalIQR =
      terminal_iqr,

    TerminalIQRFraction =
      terminal_iqr /
      total_n,

    TerminalRecoveryRate =
      terminal_recovery,

    EffectiveFractionMedian =
      safe_median(
        valid_df$effective_fraction
      ),

    EffectiveFractionIQR =
      effective_fraction_iqr,

    RightFractionMedian =
      safe_median(
        valid_df$right_fraction
      ),

    RightFractionIQR =
      right_fraction_iqr,

    MedianJaccard =
      median_jaccard,

    PassValidRate =
      pass_valid,

    PassRefRecovery =
      pass_ref_recovery,

    PassAnchorRecovery =
      pass_anchor_recovery,

    PassTerminalRecovery =
      pass_terminal_recovery,

    PassRefIQR =
      pass_ref_iqr,

    PassAnchorIQR =
      pass_anchor_iqr,

    PassTerminalIQR =
      pass_terminal_iqr,

    PassRightFractionIQR =
      pass_right_fraction_iqr,

    PassStability =
      pass_stability,

    stringsAsFactors = FALSE
  )


  list(
    full_geometry =
      full_geometry,

    bootstrap_df =
      bootstrap_df,

    stability_summary =
      stability_summary
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
    !is.finite(
      empirical_var
    )
  ] <- 0


  mu <- pmax(
    mu,
    0
  )


  empirical_var <- pmax(
    empirical_var,
    0
  )


  excess_variance <- pmax(
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
        excess_variance
      ),

    NB2_NB1 =
      log1p(
        excess_variance
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
      nrow(
        left_df
      ),

    right_n =
      nrow(
        right_df
      ),


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
    geometry,
    stability_summary,
    region_summary,
    total_n) {

  if (
    !isTRUE(
      geometry$valid
    )
  ) {
    stop(
      "Validation failed: full geometry is invalid."
    )
  }


  if (
    !isTRUE(
      stability_summary$PassStability[
        1L
      ]
    )
  ) {
    stop(
      "Validation failed: geometry did not pass bootstrap stability."
    )
  }


  anchor <- geometry$interval$anchor
  ref <- geometry$interval$ref
  terminal <- geometry$interval$terminal


  if (
    !(
      anchor <
      ref &&
      ref <
      terminal
    )
  ) {
    stop(
      "Validation failed: expected Anchor < Ref < Terminal."
    )
  }


  rounded_crossings <- unique(
    as.integer(
      round(
        geometry$zero_df$crossing_rank
      )
    )
  )


  if (
    !(
      anchor %in%
      rounded_crossings
    )
  ) {
    stop(
      "Validation failed: Anchor is not a d2 crossing."
    )
  }


  if (
    !(
      terminal %in%
      rounded_crossings
    )
  ) {
    stop(
      "Validation failed: Terminal is not a d2 crossing."
    )
  }


  expected_ref <- total_n -
    geometry$reference_info$n_eff_n +
    1L


  if (
    ref !=
    expected_ref
  ) {
    stop(
      "Validation failed: Ref does not match N_eff."
    )
  }


  right_idx <- seq.int(
    anchor,
    total_n
  )


  right_n <- length(
    right_idx
  )


  left_end <- anchor -
    1L


  left_start <- left_end -
    right_n +
    1L


  if (
    left_start < 1L
  ) {
    stop(
      "Validation failed: matched LEFT block extends below rank 1."
    )
  }


  left_idx <- seq.int(
    left_start,
    left_end
  )


  if (
    length(
      left_idx
    ) !=
    length(
      right_idx
    )
  ) {
    stop(
      "Validation failed: LEFT and RIGHT sizes differ."
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
      nrow(
        left_df
      ),

    right_n =
      nrow(
        right_df
      ),

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
    nm in names(
      expected
    )
  ) {

    chk <- all.equal(
      as.numeric(
        region_summary[[nm]][
          1L
        ]
      ),

      as.numeric(
        expected[[nm]]
      ),

      tolerance =
        1e-10
    )


    if (
      !isTRUE(
        chk
      )
    ) {

      stop(
        "Validation failed for ",
        nm,
        ": ",
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
    geometry,
    feature_df,
    region_summary,
    stability_summary,
    out_file,
    fig_tag,
    rank_tag,
    metric_tag) {

  total_n <- nrow(
    feature_df
  )


  anchor <- geometry$interval$anchor
  ref <- geometry$interval$ref
  terminal <- geometry$interval$terminal


  left_min <- region_summary$left_start
  left_max <- region_summary$left_end

  right_min <- region_summary$right_start
  right_max <- region_summary$right_end


  variance_df <- geometry$variance_df


  event_df <- data.frame(
    event =
      factor(
        EVENT_LEVELS,
        levels = EVENT_LEVELS
      ),

    rank =
      c(
        anchor,
        ref,
        terminal
      ),

    y =
      approx(
        x =
          variance_df$rank,

        y =
          variance_df$smooth_log1p_empirical_variance,

        xout =
          c(
            anchor,
            ref,
            terminal
          ),

        rule =
          2
      )$y,

    stringsAsFactors = FALSE
  )


  vline_df <- data.frame(
    event =
      factor(
        EVENT_LEVELS,
        levels = EVENT_LEVELS
      ),

    xint =
      c(
        anchor,
        ref,
        terminal
      ),

    stringsAsFactors = FALSE
  )


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


  pos <- compute_text_positions(
    total_n,
    anchor
  )


  variance_range <- range(
    variance_df$smooth_log1p_empirical_variance,
    finite = TRUE
  )


  variance_span <- diff(
    variance_range
  )


  if (
    !is.finite(
      variance_span
    ) ||
    variance_span <= 0
  ) {
    variance_span <- 1
  }


  metric_range <- range(
    nb_long$value,
    finite = TRUE
  )


  metric_span <- diff(
    metric_range
  )


  if (
    !is.finite(
      metric_span
    ) ||
    metric_span <= 0
  ) {
    metric_span <- 1
  }


  box1 <- paste(
    fig_tag,
    rank_tag,
    metric_tag,
    "Blue = LEFT",
    "Green = RIGHT",
    "Grey = geometric interval",
    sep = "\n"
  )


  box2 <- paste0(
    "Effective PC1 n = ",
    geometry$reference_info$n_eff_n,
    " (",
    round(
      100 *
      geometry$reference_info$n_eff_fraction,
      1
    ),
    "%)\n",

    "Anchor = ",
    anchor,
    "\n",

    "Ref = ",
    ref,
    "\n",

    "Terminal = ",
    terminal,
    "\n",

    "RIGHT n = ",
    region_summary$right_n,
    "\n",

    "Valid bootstrap = ",
    round(
      stability_summary$BootstrapValidRate,
      3
    ),
    "\n",

    "Anchor recovery = ",
    round(
      stability_summary$AnchorRecoveryRate,
      3
    ),
    "\n",

    "Ref recovery = ",
    round(
      stability_summary$RefRecoveryRate,
      3
    ),
    "\n",

    "Terminal recovery = ",
    round(
      stability_summary$TerminalRecoveryRate,
      3
    )
  )


  box3 <- paste(
    "RIGHT = Anchor to end",
    "LEFT = equal matched block",
    "NB2 metrics evaluated after stability",
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
  # PANEL 1
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

    annotate(
      "rect",
      xmin = anchor,
      xmax = terminal,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$interval_fill,
      alpha = 0.18
    ) +

    geom_line(
      color =
        COL$var_curve,
      linewidth =
        1.0
    ) +

    geom_vline(
      data =
        vline_df,

      aes(
        xintercept =
          xint,

        color =
          event,

        linetype =
          event
      ),

      linewidth =
        0.9
    ) +

    geom_point(
      data =
        event_df,

      aes(
        x =
          rank,

        y =
          y,

        color =
          event,

        shape =
          event
      ),

      inherit.aes =
        FALSE,

      size =
        3.4,

      stroke =
        1.0
    ) +

    annotate(
      "label",

      x =
        pos$left_x,

      y =
        variance_range[2L] -
        0.04 *
        variance_span,

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
        variance_range[2L] -
        0.30 *
        variance_span,

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

    scale_color_manual(
      values =
        EVENT_COLORS,

      breaks =
        EVENT_LEVELS
    ) +

    scale_shape_manual(
      values =
        EVENT_SHAPES,

      guide =
        "none"
    ) +

    scale_linetype_manual(
      values =
        EVENT_LTY,

      guide =
        "none"
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
        paste0(
          "Effective PC1 reference flanked by ",
          "bootstrap-stable curvature transitions"
        ),

      x =
        "Rank",

      y =
        "Smoothed log(1 + variance)",

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
  # PANEL 2
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

    annotate(
      "rect",
      xmin = anchor,
      xmax = terminal,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$interval_fill,
      alpha = 0.18
    ) +

    geom_vline(
      data =
        vline_df,

      aes(
        xintercept =
          xint
      ),

      color =
        "grey35",

      linetype =
        "dashed",

      linewidth =
        0.5
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
        metric_range[2L] -
        0.04 *
        metric_span,

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
        metric_range[2L] -
        0.30 *
        metric_span,

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
        paste0(
          "NB2-related quantities are evaluated only ",
          "after geometric stability is established"
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
  # PANEL 3
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
          "Points farther right indicate stronger ",
          "NB2-related corroboration"
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
    list(
      p1,
      p2,
      p3
    ),
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


  stability_obj <-
    assess_geometry_stability(
      rank_matrix =
        rank_matrix,

      metric_matrix =
        metric_matrix,

      bootstrap_n =
        BOOTSTRAP_N,

      spar =
        VAR_SPLINE_SPAR,

      seed =
        seed,

      cores =
        BOOTSTRAP_CORES
    )


  stability_summary <-
    stability_obj$stability_summary


  stability_summary <- cbind(
    data.frame(
      comp =
        comparison_name,

      arm =
        arm_name,

      track =
        track,

      stringsAsFactors = FALSE
    ),

    stability_summary
  )


  write.csv(
    stability_summary,

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


  # ---------------------------------------------------------------------------
  # Full-data geometry invalid
  # ---------------------------------------------------------------------------

  if (
    !isTRUE(
      stability_obj$full_geometry$valid
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

      Reason =
        stability_obj$full_geometry$reason,

      EffectivePC1N =
        NA_real_,

      EffectivePC1NRounded =
        NA_integer_,

      EffectivePC1Fraction =
        NA_real_,

      Anchor =
        NA_integer_,

      Ref =
        NA_integer_,

      Terminal =
        NA_integer_,

      RightSize =
        NA_integer_,

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


    message(
      "[",
      track,
      "] ",
      comparison_name,
      " ",
      arm_name,
      " | UNSTABLE: ",
      stability_obj$full_geometry$reason
    )


    return(
      list(
        summary =
          fail_summary,

        stability_summary =
          stability_summary
      )
    )
  }


  geometry <-
    stability_obj$full_geometry


  write.csv(
    geometry$zero_df,

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
  # Full geometry exists but bootstrap stability fails.
  # ---------------------------------------------------------------------------

  if (
    !isTRUE(
      stability_summary$PassStability[
        1L
      ]
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

      Reason =
        "Bootstrap geometric stability criteria not met",

      EffectivePC1N =
        geometry$reference_info$n_eff,

      EffectivePC1NRounded =
        geometry$reference_info$n_eff_n,

      EffectivePC1Fraction =
        geometry$reference_info$n_eff_fraction,

      Anchor =
        geometry$interval$anchor,

      Ref =
        geometry$interval$ref,

      Terminal =
        geometry$interval$terminal,

      IntervalWidth =
        geometry$interval$interval_width,

      RightSize =
        geometry$interval$right_n,

      RightFraction =
        geometry$interval$right_fraction,

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


    message(
      "[",
      track,
      "] ",
      comparison_name,
      " ",
      arm_name,
      " | UNSTABLE",
      " | Effective n=",
      geometry$reference_info$n_eff_n,
      " | Anchor=",
      geometry$interval$anchor,
      " Ref=",
      geometry$interval$ref,
      " Terminal=",
      geometry$interval$terminal
    )


    return(
      list(
        summary =
          fail_summary,

        stability_summary =
          stability_summary
      )
    )
  }


  # ===========================================================================
  # STABLE GEOMETRY
  #
  # Only now compute NB2-related corroboration.
  # ===========================================================================

  feature_df <-
    compute_ranked_feature_metrics(
      metric_matrix,
      geometry$rank_order
    )


  feature_df$abs_pc1_loading <-
    geometry$abs_loadings[
      geometry$rank_order
    ]


  region_summary <-
    summarize_regions(
      feature_df,
      geometry$interval$anchor,
      total_n
    )


  validate_method_level(
    feature_df =
      feature_df,

    geometry =
      geometry,

    stability_summary =
      stability_summary,

    region_summary =
      region_summary,

    total_n =
      total_n
  )


  rank_output <- feature_df %>%
    left_join(
      geometry$variance_df,
      by = "rank",
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


  cutoff_summary <- data.frame(
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

    Reason =
      NA_character_,

    EffectivePC1N =
      geometry$reference_info$n_eff,

    EffectivePC1NRounded =
      geometry$reference_info$n_eff_n,

    EffectivePC1Fraction =
      geometry$reference_info$n_eff_fraction,

    Anchor =
      geometry$interval$anchor,

    Ref =
      geometry$interval$ref,

    Terminal =
      geometry$interval$terminal,

    IntervalWidth =
      geometry$interval$interval_width,

    RightSize =
      geometry$interval$right_n,

    RightFraction =
      geometry$interval$right_fraction,

    BootstrapValidRate =
      stability_summary$BootstrapValidRate,

    RefRecoveryRate =
      stability_summary$RefRecoveryRate,

    AnchorRecoveryRate =
      stability_summary$AnchorRecoveryRate,

    TerminalRecoveryRate =
      stability_summary$TerminalRecoveryRate,

    RefIQRFraction =
      stability_summary$RefIQRFraction,

    AnchorIQRFraction =
      stability_summary$AnchorIQRFraction,

    TerminalIQRFraction =
      stability_summary$TerminalIQRFraction,

    RightFractionIQR =
      stability_summary$RightFractionIQR,

    MedianJaccard =
      stability_summary$MedianJaccard,

    left_NB2 =
      region_summary$left_NB2,

    right_NB2 =
      region_summary$right_NB2,

    diff_NB2 =
      region_summary$diff_NB2,

    left_gap =
      region_summary$left_gap,

    right_gap =
      region_summary$right_gap,

    diff_gap =
      region_summary$diff_gap,

    left_alpha =
      region_summary$left_alpha,

    right_alpha =
      region_summary$right_alpha,

    diff_alpha =
      region_summary$diff_alpha,

    Call_NB2 =
      ifelse(
        region_summary$diff_NB2 > 0,
        "RIGHT",
        "NOT_RIGHT"
      ),

    Call_Gap =
      ifelse(
        region_summary$diff_gap > 0,
        "RIGHT",
        "NOT_RIGHT"
      ),

    Call_Alpha =
      ifelse(
        region_summary$diff_alpha > 0,
        "RIGHT",
        "NOT_RIGHT"
      ),

    stringsAsFactors = FALSE
  )


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

      Anchor_lt_Ref =
        geometry$interval$anchor <
        geometry$interval$ref,

      Ref_lt_Terminal =
        geometry$interval$ref <
        geometry$interval$terminal,

      EffectiveTailAboveGuardrail =
        geometry$reference_info$n_eff_n >=
        ceiling(
          MIN_EFFECTIVE_TAIL_FRACTION *
          total_n
        ),

      EqualRegionSize =
        region_summary$left_n ==
        region_summary$right_n,

      PassStability =
        stability_summary$PassStability,

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


  build_main_figure(
    comparison_name =
      comparison_name,

    arm_name =
      arm_name,

    geometry =
      geometry,

    feature_df =
      feature_df,

    region_summary =
      region_summary,

    stability_summary =
      stability_summary,

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
    " | PASS",
    " | Effective n=",
    geometry$reference_info$n_eff_n,
    " (",
    sprintf(
      "%.1f",
      100 *
      geometry$reference_info$n_eff_fraction
    ),
    "%)",
    " | Anchor=",
    geometry$interval$anchor,
    " Ref=",
    geometry$interval$ref,
    " Terminal=",
    geometry$interval$terminal,
    " | RIGHT n=",
    geometry$interval$right_n,
    " | Jaccard(descriptive)=",
    sprintf(
      "%.3f",
      stability_summary$MedianJaccard
    )
  )


  list(
    summary =
      cutoff_summary,

    stability_summary =
      stability_summary
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
  "Fixed 5,000-feature cutoff: REMOVED"
)


message(
  "PC1 effective-tail guardrail: ",
  100 *
  MIN_EFFECTIVE_TAIL_FRACTION,
  "%"
)


message(
  "Bootstrap replicates per track: ",
  BOOTSTRAP_N
)


message(
  "Bootstrap cores: ",
  BOOTSTRAP_CORES
)


all_cutoff_rows <- list()

all_stability_rows <- list()


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
      length(
        sample_idx
      ) < 2L
    ) {

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


    # =========================================================================
    # MAIN MANUSCRIPT TRACK
    # =========================================================================

    main_rank_matrix <-
      normalize_cpm_log1p(
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


    all_cutoff_rows[[
      length(
        all_cutoff_rows
      ) +
      1L
    ]] <- main_res$summary


    all_stability_rows[[
      length(
        all_stability_rows
      ) +
      1L
    ]] <- main_res$stability_summary


    # =========================================================================
    # DESEQ2 SUPPLEMENTARY TRACK
    # =========================================================================

    if (
      RUN_DESEQ2_SUPPLEMENT
    ) {

      deseq2_obj <-
        compute_deseq2_matrices(
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
          "[DESeq2] skipped: required package(s) unavailable for ",
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


        all_cutoff_rows[[
          length(
            all_cutoff_rows
          ) +
          1L
        ]] <- deseq2_res$summary


        all_stability_rows[[
          length(
            all_stability_rows
          ) +
          1L
        ]] <- deseq2_res$stability_summary
      }
    }
  }
}


# =============================================================================
# OVERALL OUTPUTS
# =============================================================================

overall_cutoff <- bind_rows(
  all_cutoff_rows
)


overall_stability <- bind_rows(
  all_stability_rows
)


write.csv(
  overall_cutoff,

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
  "Outputs written to: ",
  OUT_ROOT
)


message(
  "The fixed 5,000-feature reference was not used."
)


message(
  paste0(
    "The reference was derived from the effective number ",
    "of PC1-contributing features."
  )
)


message(
  paste0(
    "Jaccard overlap was reported descriptively and was ",
    "not used to force pass/fail stability."
  )
)


message(
  paste0(
    "NB2-related corroboration was evaluated only for ",
    "bootstrap-stable geometric cutoffs."
  )
)


message(
  "============================================================"
)
