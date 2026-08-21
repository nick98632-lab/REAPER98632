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
# DATA-DERIVED TERMINAL TRANSITION DETECTOR
# =============================================================================
#
# PURPOSE
# -------
# This script preserves the original manuscript analysis architecture while
# removing the fixed 5,000-feature reference.
#
# 1. MAIN MANUSCRIPT TRACK
#    - Rank: log1p CPM-style library-size normalization
#    - Geometry / NB diagnostics: raw counts
#
# 2. DESEQ2 SUPPLEMENTARY TRACK
#    - Rank: DESeq2-normalized log1p counts, optionally VST
#    - Geometry / NB diagnostics: DESeq2-normalized counts
#
# SCIENTIFIC LOGIC
# ----------------
# A. Features are ranked within each arm from LOWEST to HIGHEST absolute PC1
#    loading. The high-|PC1|-loading leading edge is therefore at the RIGHT.
#
# B. Along the complete ranked axis, empirical feature-wise variance is
#    transformed as log(1 + variance).
#
# C. The old fixed reference
#
#        Ref = N - 5000 + 1
#
#    is NOT used.
#
# D. A data-derived reference is estimated as a TERMINAL CHANGE POINT in the
#    ranked log-variance trajectory. The detector uses multiscale continuous
#    segmented regression. At each smoothing scale, every mathematically
#    admissible terminal breakpoint is evaluated with a continuous hinge model:
#
#        y(x) = beta0 + beta1*x + gamma*max(0, x-c)
#
#    versus the no-break model:
#
#        y(x) = beta0 + beta1*x
#
#    The change-point model is fit on one FIXED terminal-half domain. Candidate
#    c must leave at least 5% of the full feature count on BOTH sides of the
#    breakpoint within that fixed domain. This is an anti-degeneracy condition
#    only; it does not encode a 5,000-feature target.
#
#    For each candidate, the exact reduction in residual sum of squares from
#    adding the hinge term is computed. If M is the number of observations in
#    the fixed terminal-half domain, the BIC improvement score is:
#
#        BIC_gain = M*log(SSE_linear / SSE_hinge) - log(M)
#
#    Because every candidate is fit to the SAME M observations with the SAME
#    number of parameters, this score cannot win merely by using a larger LEFT
#    or RIGHT region. This specifically avoids the 50/50 artifact produced by
#    the previous t-like / standard-error objective.
#
#    A terminal breakpoint is eligible only when:
#      - gamma > 0 (slope increases after the break),
#      - terminal fitted slope > 0,
#      - the terminal segment has a higher mean log-variance than its matched
#        immediately preceding segment, and
#      - BIC_gain > 0.
#
# E. MULTISCALE CONSENSUS
#    The breakpoint search is repeated at three prespecified smoothing levels.
#    Boundary solutions are rejected. A robust weighted-median breakpoint,
#    weighted by sqrt(BIC_gain), defines the data-derived Ref. At least two
#    smoothing scales must yield valid interior breakpoints, and their reference
#    IQR must be <= 5% of the full rank axis.
#
# F. ORIGINAL GEOMETRIC RULE IS THEN RESTORED
#    The original spline-based second derivative is computed at spar = 0.60.
#    The final custom interval is:
#
#      Anchor   = nearest d2 zero crossing immediately LEFT of data-derived Ref
#      Terminal = nearest d2 zero crossing immediately RIGHT of data-derived Ref
#
#    with:
#
#      Anchor < Ref < Terminal
#
# G. ORIGINAL LEADING-EDGE DEFINITION IS PRESERVED
#
#      RIGHT = Anchor through the terminal rank N
#      LEFT  = equal-sized block immediately preceding Anchor
#
#    Terminal is a geometric boundary around Ref; it does NOT truncate RIGHT.
#
# H. BOOTSTRAP STABILITY
#    Samples are resampled within each arm. Every replicate recomputes:
#      PC1 -> ranking -> variance trajectory -> multiscale terminal Ref ->
#      d2 crossings -> Anchor / Terminal.
#
#    Stability requires:
#      - >= 90% valid bootstrap replicates
#      - >= 80% recovery within 3% of N for Ref, Anchor, and Terminal
#      - Ref, Anchor, and Terminal IQR <= 3% of N
#      - RIGHT-fraction IQR <= 3 percentage points
#
#    Jaccard feature overlap is reported descriptively only.
#
# I. NB2-RELATED CORROBORATION
#    Only after geometry is fixed, LEFT and RIGHT are compared using:
#
#      NB2      = log(1 + max(variance - mu, 0))
#      NB2-NB1  = log(1 + max(variance - mu, 0)) - log(1 + mu)
#      alpha*mu = log(1 + alpha*mu)
#
#      alpha = max((variance - mu) / mu^2, 0)
#
#    These are descriptive moment-based diagnostics, not likelihood-ratio tests.
#    They do NOT participate in terminal-transition detection or bootstrap
#    selection.
# =============================================================================

# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT   <- "/root/REAPER98632/exports/manuscript_terminal_transition_final"

VAR_SPLINE_SPAR <- 0.60

TRANSITION_SPARS <- c(0.50, 0.60, 0.70)
MIN_VALID_TRANSITION_SCALES <- 2L
MAX_MULTISCALE_REF_IQR_FRACTION <- 0.05

MIN_TRANSITION_SEGMENT_FRACTION <- 0.05
MIN_TRANSITION_SEGMENT_ABSOLUTE <- 100L

DENSE_GRID_MULTIPLIER <- 4L
DENSE_GRID_MIN <- 5000L

BOOTSTRAP_N <- 500L
BOOTSTRAP_SEED_BASE <- 20260820L
MIN_BOOTSTRAP_VALID_RATE <- 0.90
POSITION_TOLERANCE_FRACTION <- 0.03
MIN_POSITION_RECOVERY_RATE <- 0.80
MAX_POSITION_IQR_FRACTION <- 0.03
MAX_RIGHT_FRACTION_IQR <- 0.03

PNG_WIDTH_IN  <- 14
PNG_HEIGHT_IN <- 10.8
PNG_DPI       <- 260

RUN_DESEQ2_SUPPLEMENT <- TRUE
DESEQ2_RANK_METHOD <- "normalized_log1p"

COMPARISONS <- list(
  RT0_ZT6  = list(control = "^R0_", treatment = "^ZT6_"),
  RT2_ZT8  = list(control = "^R2_", treatment = "^ZT8_"),
  RT4_ZT10 = list(control = "^R4_", treatment = "^ZT10_"),
  RT8_ZT14 = list(control = "^R8_", treatment = "^ZT14_")
)

DETECTED_CORES <- suppressWarnings(parallel::detectCores(logical = FALSE))

if (is.na(DETECTED_CORES) || DETECTED_CORES < 2L) {
  BOOTSTRAP_CORES <- 1L
} else {
  BOOTSTRAP_CORES <- min(4L, DETECTED_CORES - 1L)
}

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

# =============================================================================
# COLORS
# =============================================================================

COL <- list(
  var_curve     = "#117A65",
  nb2           = "#1B9E77",
  nb_gap        = "#CC1E8C",
  alpha_mu      = "#386CB0",
  left_fill     = "#CBE3F8",
  right_fill    = "#DDF2D5",
  interval_fill = "#9E9E9E",
  anchor         = "#000000",
  ref            = "#E69F00",
  terminal       = "#D95F02",
  left_pt        = "#5B8FD1",
  right_pt       = "#43A047"
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
# INPUT / NORMALIZATION
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


read_count_matrix <- function(path) {

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

  first_num <- first_numeric_col_index(
    raw_df
  )

  if (is.na(first_num)) {
    stop(
      "No numeric count columns detected."
    )
  }

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

  count_df <- raw_df[
    ,
    first_num:ncol(raw_df),
    drop = FALSE
  ]

  count_mat <- as.matrix(
    count_df
  )

  suppressWarnings({
    storage.mode(count_mat) <- "numeric"
  })

  rownames(count_mat) <- feature_ids

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

  keep <- rowSums(
    count_mat
  ) > 0

  count_mat[
    keep,
    ,
    drop = FALSE
  ]
}


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

  dds <- tryCatch(

    DESeq2::estimateSizeFactors(
      dds
    ),

    error =
      function(e) {

        message(
          "DESeq2 default size-factor estimation failed; ",
          "retrying with type='poscounts'."
        )

        DESeq2::estimateSizeFactors(
          dds,
          type = "poscounts"
        )
      }
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
# NUMERICAL HELPERS
# =============================================================================

row_variance_fast <- function(
    mat) {

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
    n *
    mu^2

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


safe_median <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  median(x)
}


safe_iqr <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  stats::IQR(x)
}


weighted_median <- function(
    x,
    w) {

  ok <- (
    is.finite(x) &
    is.finite(w) &
    w > 0
  )

  x <- x[ok]
  w <- w[ok]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  ord <- order(x)

  x <- x[ord]
  w <- w[ord]

  cw <- cumsum(w) /
    sum(w)

  x[
    which(
      cw >= 0.5
    )[1L]
  ]
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


jaccard_similarity <- function(
    a,
    b) {

  a <- unique(a)
  b <- unique(b)

  u <- union(
    a,
    b
  )

  if (length(u) == 0L) {
    return(NA_real_)
  }

  length(
    intersect(
      a,
      b
    )
  ) /
    length(u)
}


suffix_sum <- function(x) {

  rev(
    cumsum(
      rev(x)
    )
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
      "PC1 requires at least two features ",
      "and two samples."
    )
  }

  X <- t(
    norm_mat_arm
  )

  Xc <- sweep(
    X,
    2L,
    colMeans(X),
    "-"
  )

  gram <- tcrossprod(
    Xc
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
      Xc,
      u1
    )
  ) /
    sqrt(
      lambda1
    )

  names(loading) <- colnames(
    Xc
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

  dense_x <- seq(

    min(ranks),

    max(ranks),

    length.out =
      max(
        DENSE_GRID_MIN,
        length(ranks) *
        DENSE_GRID_MULTIPLIER
      )
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


find_d2_zero_crossings <- function(
    dense_df) {

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

  idx <- which(
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

  if (length(idx) > 0L) {

    frac <- abs(
      a[idx]
    ) /
      (
        abs(a[idx]) +
        abs(b[idx])
      )

    crossings <-
      xa[idx] +
      frac *
      (
        xb[idx] -
        xa[idx]
      )
  }

  exact_idx <- which(
    y == 0
  )

  if (length(exact_idx) > 0L) {

    crossings <- c(
      crossings,
      x[exact_idx]
    )
  }

  crossings <- sort(
    unique(
      round(
        crossings[
          is.finite(crossings)
        ],
        8L
      )
    )
  )

  data.frame(

    crossing_rank =
      crossings,

    crossing_type =
      "sign_change",

    stringsAsFactors = FALSE
  )
}

# =============================================================================
# TERMINAL CHANGE-POINT MODEL
# =============================================================================

fit_terminal_hinge_breakpoint <- function(
    y,
    min_segment_fraction = MIN_TRANSITION_SEGMENT_FRACTION,
    min_segment_absolute = MIN_TRANSITION_SEGMENT_ABSOLUTE) {

  y <- as.numeric(y)
  N <- length(y)

  if (
    N < 10L ||
    any(
      !is.finite(y)
    )
  ) {

    return(
      list(
        valid = FALSE,
        reason = "invalid_transition_series"
      )
    )
  }

  # ---------------------------------------------------------------------------
  # Terminal-domain definition.
  #
  # The original matched LEFT/RIGHT definition means an Anchor left of the
  # midpoint cannot support RIGHT=Anchor->N plus an equal-sized preceding
  # LEFT block. Therefore the terminal half is the mathematically relevant
  # change-point domain.
  # ---------------------------------------------------------------------------

  domain_start <- as.integer(
    ceiling(
      (
        N +
        2L
      ) /
        2L
    )
  )

  min_segment_n <- max(

    as.integer(
      min_segment_absolute
    ),

    as.integer(
      ceiling(
        min_segment_fraction *
        N
      )
    ),

    2L
  )

  # ---------------------------------------------------------------------------
  # Every candidate is compared on ONE FIXED terminal domain.
  #
  # Both pre-break and post-break portions must contain at least
  # min_segment_n points.
  # ---------------------------------------------------------------------------

  candidate_min <-
    domain_start +
    min_segment_n -
    1L

  candidate_max <-
    N -
    min_segment_n +
    1L

  if (
    candidate_min >=
    candidate_max
  ) {

    return(
      list(
        valid = FALSE,
        reason = "no_admissible_terminal_breakpoints"
      )
    )
  }

  domain_rank <- seq.int(
    domain_start,
    N
  )

  yd <- y[
    domain_rank
  ]

  M <- length(
    yd
  )

  xd <- (
    seq_len(M) -
    1
  ) /
    (
      M -
      1
    )

  candidate_rank <- seq.int(
    candidate_min,
    candidate_max
  )

  candidate_pos <-
    candidate_rank -
    domain_start +
    1L

  n <- as.numeric(M)

  sx <- sum(xd)
  sxx <- sum(xd * xd)
  sy <- sum(yd)
  sxy <- sum(xd * yd)

  det0 <-
    n *
    sxx -
    sx *
    sx

  if (
    !is.finite(det0) ||
    abs(det0) <=
    .Machine$double.eps
  ) {

    return(
      list(
        valid = FALSE,
        reason = "singular_baseline_linear_model"
      )
    )
  }

  inv00 <- sxx / det0
  inv01 <- -sx / det0
  inv11 <- n / det0

  beta0_intercept <-
    inv00 *
    sy +
    inv01 *
    sxy

  beta0_slope <-
    inv01 *
    sy +
    inv11 *
    sxy

  residual0 <-
    yd -
    beta0_intercept -
    beta0_slope *
    xd

  sse0 <- sum(
    residual0 *
    residual0
  )

  if (
    !is.finite(sse0) ||
    sse0 <=
    .Machine$double.eps
  ) {

    return(
      list(
        valid = FALSE,
        reason = "degenerate_baseline_sse"
      )
    )
  }

  # ---------------------------------------------------------------------------
  # Suffix statistics permit exact evaluation of every candidate in O(N)
  # rather than fitting thousands of separate lm() objects.
  # ---------------------------------------------------------------------------

  sx_suf <- suffix_sum(
    xd
  )

  sxx_suf <- suffix_sum(
    xd *
    xd
  )

  sy_suf <- suffix_sum(
    yd
  )

  sxy_suf <- suffix_sum(
    xd *
    yd
  )

  tail_start_pos <-
    candidate_pos +
    1L

  n_tail_hinge <-
    M -
    candidate_pos

  xc <- xd[
    candidate_pos
  ]

  sx_tail <- sx_suf[
    tail_start_pos
  ]

  sxx_tail <- sxx_suf[
    tail_start_pos
  ]

  sy_tail <- sy_suf[
    tail_start_pos
  ]

  sxy_tail <- sxy_suf[
    tail_start_pos
  ]

  sum_h <-
    sx_tail -
    n_tail_hinge *
    xc

  sum_xh <-
    sxx_tail -
    xc *
    sx_tail

  sum_h2 <-
    sxx_tail -
    2 *
    xc *
    sx_tail +
    n_tail_hinge *
    xc^2

  sum_yh <-
    sxy_tail -
    xc *
    sy_tail

  # ---------------------------------------------------------------------------
  # Frisch-Waugh-Lovell residualization.
  #
  # This gives the exact one-hinge regression improvement while avoiding
  # thousands of repeated model fits.
  # ---------------------------------------------------------------------------

  q0 <-
    inv00 *
    sum_h +
    inv01 *
    sum_xh

  q1 <-
    inv01 *
    sum_h +
    inv11 *
    sum_xh

  hMh <-
    sum_h2 -
    sum_h *
    q0 -
    sum_xh *
    q1

  hMy <-
    sum_yh -
    sum_h *
    beta0_intercept -
    sum_xh *
    beta0_slope

  valid_h <- (
    is.finite(hMh) &
    hMh >
    .Machine$double.eps &
    is.finite(hMy)
  )

  gamma <- rep(
    NA_real_,
    length(candidate_rank)
  )

  gamma[
    valid_h
  ] <- hMy[
    valid_h
  ] /
    hMh[
      valid_h
    ]

  pre_slope <- rep(
    NA_real_,
    length(candidate_rank)
  )

  pre_slope[
    valid_h
  ] <-
    beta0_slope -
    q1[
      valid_h
    ] *
    gamma[
      valid_h
    ]

  post_slope <-
    pre_slope +
    gamma

  sse1 <- rep(
    NA_real_,
    length(candidate_rank)
  )

  sse1[
    valid_h
  ] <-
    sse0 -
    (
      hMy[
        valid_h
      ]^2 /
      hMh[
        valid_h
      ]
    )

  bic_gain <- rep(
    NA_real_,
    length(candidate_rank)
  )

  ok_sse <- (
    valid_h &
    is.finite(sse1) &
    sse1 >
    .Machine$double.eps &
    sse1 <
    sse0
  )

  bic_gain[
    ok_sse
  ] <-
    M *
    log(
      sse0 /
      sse1[
        ok_sse
      ]
    ) -
    log(M)

  # ---------------------------------------------------------------------------
  # Persistent terminal-level check using the original matched-region concept.
  # No NB quantities are used here.
  # ---------------------------------------------------------------------------

  prefix_y <- c(
    0,
    cumsum(y)
  )

  right_n <-
    N -
    candidate_rank +
    1L

  left_start <-
    candidate_rank -
    right_n

  left_end <-
    candidate_rank -
    1L

  left_sum <-
    prefix_y[
      left_end +
      1L
    ] -
    prefix_y[
      left_start
    ]

  right_sum <-
    prefix_y[
      N +
      1L
    ] -
    prefix_y[
      candidate_rank
    ]

  left_mean <-
    left_sum /
    right_n

  right_mean <-
    right_sum /
    right_n

  terminal_minus_left_mean <-
    right_mean -
    left_mean

  eligible <- (
    is.finite(bic_gain) &
    bic_gain > 0 &
    is.finite(gamma) &
    gamma > 0 &
    is.finite(post_slope) &
    post_slope > 0 &
    is.finite(
      terminal_minus_left_mean
    ) &
    terminal_minus_left_mean > 0
  )

  scan_df <- data.frame(

    reference_rank =
      candidate_rank,

    right_n =
      right_n,

    right_fraction =
      right_n /
      N,

    gamma =
      gamma,

    pre_slope =
      pre_slope,

    post_slope =
      post_slope,

    terminal_minus_left_mean =
      terminal_minus_left_mean,

    BIC_gain =
      bic_gain,

    eligible =
      eligible,

    stringsAsFactors = FALSE
  )

  if (!any(eligible)) {

    return(
      list(
        valid = FALSE,
        reason = "no_eligible_terminal_breakpoint",
        scan_df = scan_df
      )
    )
  }

  idx_eligible <- which(
    eligible
  )

  best_idx <- idx_eligible[
    which.max(
      bic_gain[
        idx_eligible
      ]
    )
  ]

  best <- scan_df[
    best_idx,
    ,
    drop = FALSE
  ]

  boundary_hit <-
    best$reference_rank[
      1L
    ] %in%
    c(
      candidate_min,
      candidate_max
    )

  list(

    valid = TRUE,

    reason =
      NA_character_,

    best =
      best,

    scan_df =
      scan_df,

    domain_start =
      domain_start,

    candidate_min =
      candidate_min,

    candidate_max =
      candidate_max,

    boundary_hit =
      boundary_hit
  )
}

# =============================================================================
# MULTISCALE TERMINAL-TRANSITION CONSENSUS
# =============================================================================

detect_terminal_transition <- function(
    ranked_log_variance,
    spars = TRANSITION_SPARS,
    min_valid_scales = MIN_VALID_TRANSITION_SCALES,
    max_iqr_fraction = MAX_MULTISCALE_REF_IQR_FRACTION) {

  y <- as.numeric(
    ranked_log_variance
  )

  N <- length(y)

  x <- seq_len(
    N
  )

  scale_rows <- vector(
    "list",
    length(spars)
  )

  for (
    i in seq_along(spars)
  ) {

    sp <- spars[i]

    one <- tryCatch({

      fit <- stats::smooth.spline(
        x = x,
        y = y,
        spar = sp
      )

      ys <- as.numeric(
        stats::predict(
          fit,
          x = x,
          deriv = 0
        )$y
      )

      bp <- fit_terminal_hinge_breakpoint(
        ys
      )

      if (
        !isTRUE(
          bp$valid
        )
      ) {

        data.frame(

          spar =
            sp,

          valid =
            FALSE,

          reason =
            bp$reason,

          boundary_hit =
            NA,

          reference_rank =
            NA_integer_,

          BIC_gain =
            NA_real_,

          gamma =
            NA_real_,

          pre_slope =
            NA_real_,

          post_slope =
            NA_real_,

          terminal_minus_left_mean =
            NA_real_,

          stringsAsFactors = FALSE
        )

      } else {

        b <- bp$best

        data.frame(

          spar =
            sp,

          valid =
            !isTRUE(
              bp$boundary_hit
            ),

          reason =
            if (
              isTRUE(
                bp$boundary_hit
              )
            ) {
              "boundary_optimum"
            } else {
              NA_character_
            },

          boundary_hit =
            bp$boundary_hit,

          reference_rank =
            as.integer(
              b$reference_rank[
                1L
              ]
            ),

          BIC_gain =
            b$BIC_gain[
              1L
            ],

          gamma =
            b$gamma[
              1L
            ],

          pre_slope =
            b$pre_slope[
              1L
            ],

          post_slope =
            b$post_slope[
              1L
            ],

          terminal_minus_left_mean =
            b$terminal_minus_left_mean[
              1L
            ],

          stringsAsFactors = FALSE
        )
      }

    }, error = function(e) {

      data.frame(

        spar =
          sp,

        valid =
          FALSE,

        reason =
          paste0(
            "scale_error: ",
            conditionMessage(e)
          ),

        boundary_hit =
          NA,

        reference_rank =
          NA_integer_,

        BIC_gain =
          NA_real_,

        gamma =
          NA_real_,

        pre_slope =
          NA_real_,

        post_slope =
          NA_real_,

        terminal_minus_left_mean =
          NA_real_,

        stringsAsFactors = FALSE
      )
    })

    scale_rows[[i]] <- one
  }

  scale_df <- bind_rows(
    scale_rows
  )

  valid_df <- scale_df[
    scale_df$valid,
    ,
    drop = FALSE
  ]

  if (
    nrow(valid_df) <
    min_valid_scales
  ) {

    return(
      list(

        valid =
          FALSE,

        reason =
          "insufficient_valid_multiscale_breakpoints",

        reference_rank =
          NA_integer_,

        scale_df =
          scale_df,

        scale_iqr =
          NA_real_,

        scale_iqr_fraction =
          NA_real_
      )
    )
  }

  refs <- valid_df$reference_rank

  weights <- sqrt(
    pmax(
      valid_df$BIC_gain,
      .Machine$double.eps
    )
  )

  ref_consensus <- weighted_median(
    refs,
    weights
  )

  ref_iqr <- safe_iqr(
    refs
  )

  ref_iqr_fraction <-
    ref_iqr /
    N

  if (
    !is.finite(
      ref_consensus
    )
  ) {

    return(
      list(

        valid =
          FALSE,

        reason =
          "invalid_multiscale_consensus",

        reference_rank =
          NA_integer_,

        scale_df =
          scale_df,

        scale_iqr =
          ref_iqr,

        scale_iqr_fraction =
          ref_iqr_fraction
      )
    )
  }

  if (
    !is.finite(
      ref_iqr_fraction
    ) ||
    ref_iqr_fraction >
    max_iqr_fraction
  ) {

    return(
      list(

        valid =
          FALSE,

        reason =
          "multiscale_reference_disagreement",

        reference_rank =
          as.integer(
            round(
              ref_consensus
            )
          ),

        scale_df =
          scale_df,

        scale_iqr =
          ref_iqr,

        scale_iqr_fraction =
          ref_iqr_fraction
      )
    )
  }

  list(

    valid =
      TRUE,

    reason =
      NA_character_,

    reference_rank =
      as.integer(
        round(
          ref_consensus
        )
      ),

    scale_df =
      scale_df,

    scale_iqr =
      ref_iqr,

    scale_iqr_fraction =
      ref_iqr_fraction,

    median_BIC_gain =
      safe_median(
        valid_df$BIC_gain
      ),

    min_reference =
      min(refs),

    max_reference =
      max(refs),

    valid_scales =
      nrow(valid_df),

    total_scales =
      length(spars)
  )
}

# =============================================================================
# ORIGINAL GEOMETRIC INTERVAL AROUND DATA-DERIVED REFERENCE
# =============================================================================

select_custom_interval <- function(
    zero_df,
    reference_rank,
    total_n) {

  if (
    nrow(zero_df) == 0L
  ) {

    stop(
      "No d2 sign-change crossings found."
    )
  }

  reference_rank <- as.integer(
    round(
      reference_rank
    )
  )

  min_anchor_for_match <- as.integer(
    ceiling(
      (
        total_n +
        2L
      ) /
        2L
    )
  )

  left_candidates <- zero_df$crossing_rank[
    zero_df$crossing_rank <
      reference_rank &
    zero_df$crossing_rank >=
      min_anchor_for_match
  ]

  right_candidates <- zero_df$crossing_rank[
    zero_df$crossing_rank >
      reference_rank
  ]

  if (
    length(left_candidates) == 0L
  ) {

    stop(
      "No admissible left d2 crossing found."
    )
  }

  if (
    length(right_candidates) == 0L
  ) {

    stop(
      "No right d2 crossing found."
    )
  }

  anchor_rank <- as.integer(
    round(
      max(
        left_candidates
      )
    )
  )

  terminal_rank <- as.integer(
    round(
      min(
        right_candidates
      )
    )
  )

  anchor_rank <- max(
    min_anchor_for_match,
    anchor_rank
  )

  terminal_rank <- min(
    total_n,
    terminal_rank
  )

  if (
    !(
      anchor_rank <
      reference_rank
    )
  ) {

    stop(
      "Invalid interval: Anchor must lie left of Ref."
    )
  }

  if (
    !(
      reference_rank <
      terminal_rank
    )
  ) {

    stop(
      "Invalid interval: Ref must lie left of Terminal."
    )
  }

  right_n <-
    total_n -
    anchor_rank +
    1L

  left_start <-
    anchor_rank -
    right_n

  if (
    left_start < 1L
  ) {

    stop(
      "Invalid interval: matched LEFT block ",
      "cannot be constructed."
    )
  }

  list(

    anchor =
      anchor_rank,

    ref =
      reference_rank,

    terminal =
      terminal_rank,

    interval_min =
      anchor_rank,

    interval_max =
      terminal_rank
  )
}

# =============================================================================
# COMPLETE GEOMETRIC SELECTION
# =============================================================================

compute_geometric_selection <- function(
    rank_matrix,
    metric_matrix) {

  if (
    nrow(rank_matrix) !=
    nrow(metric_matrix) ||
    ncol(rank_matrix) !=
    ncol(metric_matrix)
  ) {

    stop(
      "Rank and metric matrices have ",
      "incompatible dimensions."
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

  abs_loadings <- compute_abs_pc1_loadings(
    rank_matrix
  )

  rank_order <- order(
    abs_loadings,
    decreasing = FALSE
  )

  total_n <- length(
    rank_order
  )

  variance_df <- compute_ranked_variance_curve(

    metric_mat_arm =
      metric_matrix,

    rank_order =
      rank_order,

    spar =
      VAR_SPLINE_SPAR
  )

  detector <- detect_terminal_transition(
    variance_df$log1p_empirical_variance
  )

  if (
    !isTRUE(
      detector$valid
    )
  ) {

    return(
      list(

        valid =
          FALSE,

        reason =
          detector$reason,

        abs_loadings =
          abs_loadings,

        rank_order =
          rank_order,

        total_n =
          total_n,

        variance_df =
          variance_df,

        detector =
          detector,

        zero_df =
          data.frame(),

        interval_info =
          NULL
      )
    )
  }

  dense_curve_df <- attr(
    variance_df,
    "dense_curve_df"
  )

  zero_df <- find_d2_zero_crossings(
    dense_curve_df
  )

  interval_info <- tryCatch(

    select_custom_interval(

      zero_df =
        zero_df,

      reference_rank =
        detector$reference_rank,

      total_n =
        total_n
    ),

    error =
      function(e) {
        e
      }
  )

  if (
    inherits(
      interval_info,
      "error"
    )
  ) {

    return(
      list(

        valid =
          FALSE,

        reason =
          paste0(
            "interval_error: ",
            conditionMessage(
              interval_info
            )
          ),

        abs_loadings =
          abs_loadings,

        rank_order =
          rank_order,

        total_n =
          total_n,

        variance_df =
          variance_df,

        detector =
          detector,

        zero_df =
          zero_df,

        interval_info =
          NULL
      )
    )
  }

  list(

    valid =
      TRUE,

    reason =
      NA_character_,

    abs_loadings =
      abs_loadings,

    rank_order =
      rank_order,

    total_n =
      total_n,

    variance_df =
      variance_df,

    detector =
      detector,

    zero_df =
      zero_df,

    interval_info =
      interval_info
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

  nb2_var_minus_mu <- pmax(
    empirical_var -
    mu,
    0
  )

  alpha_hat <- rep(
    0,
    length(mu)
  )

  pos <- mu > 0

  alpha_hat[
    pos
  ] <- pmax(

    (
      empirical_var[
        pos
      ] -
      mu[
        pos
      ]
    ) /
      (
        mu[
          pos
        ]^2
      ),

    0
  )

  alpha_mu_val <-
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
        nb2_var_minus_mu
      ),

    NB2_NB1 =
      log1p(
        nb2_var_minus_mu
      ) -
      log1p(
        mu
      ),

    alpha_mu =
      log1p(
        alpha_mu_val
      ),

    stringsAsFactors = FALSE
  )
}


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

  left_end <-
    anchor_rank -
    1L

  left_start <-
    left_end -
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


validate_method_level <- function(
    feature_df,
    interval_info,
    region_summary,
    total_n) {

  anchor <- interval_info$anchor
  ref <- interval_info$ref
  terminal <- interval_info$terminal

  if (
    !(
      anchor <
      ref &&
      ref <
      terminal
    )
  ) {

    stop(
      "Validation failed: expected ",
      "Anchor < Ref < Terminal."
    )
  }

  right_idx <- seq.int(
    anchor,
    total_n
  )

  right_n <- length(
    right_idx
  )

  left_end <-
    anchor -
    1L

  left_start <-
    left_end -
    right_n +
    1L

  if (
    left_start < 1L
  ) {

    stop(
      "Validation failed: LEFT block ",
      "extends below rank 1."
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
      "Validation failed: LEFT and RIGHT ",
      "blocks do not have equal size."
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
# BOOTSTRAP STABILITY OF COMPLETE GEOMETRY
# =============================================================================

assess_geometry_stability <- function(
    rank_matrix,
    metric_matrix,
    full_geometry,
    bootstrap_n = BOOTSTRAP_N,
    seed = BOOTSTRAP_SEED_BASE,
    cores = BOOTSTRAP_CORES) {

  if (
    !isTRUE(
      full_geometry$valid
    )
  ) {

    stop(
      "Cannot bootstrap invalid full-data geometry."
    )
  }

  N <- nrow(
    rank_matrix
  )

  sample_n <- ncol(
    rank_matrix
  )

  full_interval <-
    full_geometry$interval_info

  full_right_idx <-
    full_geometry$rank_order[
      seq.int(
        full_interval$anchor,
        N
      )
    ]

  full_right_features <-
    rownames(
      rank_matrix
    )[
      full_right_idx
    ]

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

  worker <- function(b) {

    idx <- bootstrap_indices[[
      b
    ]]

    invalid_row <- function(reason) {

      data.frame(

        bootstrap =
          b,

        valid =
          FALSE,

        invalid_reason =
          reason,

        ref =
          NA_integer_,

        anchor =
          NA_integer_,

        terminal =
          NA_integer_,

        right_n =
          NA_integer_,

        right_fraction =
          NA_real_,

        multiscale_ref_iqr_fraction =
          NA_real_,

        jaccard =
          NA_real_,

        stringsAsFactors = FALSE
      )
    }

    if (
      length(
        unique(idx)
      ) < 2L
    ) {

      return(
        invalid_row(
          "fewer_than_2_unique_samples"
        )
      )
    }

    g <- tryCatch(

      compute_geometric_selection(

        rank_matrix[
          ,
          idx,
          drop = FALSE
        ],

        metric_matrix[
          ,
          idx,
          drop = FALSE
        ]
      ),

      error =
        function(e) {
          NULL
        }
    )

    if (is.null(g)) {

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

    interval <-
      g$interval_info

    right_idx <-
      g$rank_order[
        seq.int(
          interval$anchor,
          N
        )
      ]

    right_features <-
      rownames(
        rank_matrix
      )[
        right_idx
      ]

    data.frame(

      bootstrap =
        b,

      valid =
        TRUE,

      invalid_reason =
        NA_character_,

      ref =
        interval$ref,

      anchor =
        interval$anchor,

      terminal =
        interval$terminal,

      right_n =
        N -
        interval$anchor +
        1L,

      right_fraction =
        (
          N -
          interval$anchor +
          1L
        ) /
        N,

      multiscale_ref_iqr_fraction =
        g$detector$scale_iqr_fraction,

      jaccard =
        jaccard_similarity(
          full_right_features,
          right_features
        ),

      stringsAsFactors = FALSE
    )
  }

  ids <- seq_len(
    bootstrap_n
  )

  if (
    .Platform$OS.type == "unix" &&
    cores > 1L
  ) {

    boot_list <- parallel::mclapply(

      ids,

      worker,

      mc.cores =
        cores,

      mc.preschedule =
        TRUE,

      mc.set.seed =
        FALSE
    )

  } else {

    boot_list <- lapply(
      ids,
      worker
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

  valid_rate <-
    valid_n /
    bootstrap_n

  tolerance_n <- ceiling(
    POSITION_TOLERANCE_FRACTION *
    N
  )

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

    ref_q <- as.numeric(
      stats::quantile(
        valid_df$ref,
        c(
          0.025,
          0.50,
          0.975
        ),
        na.rm = TRUE,
        names = FALSE
      )
    )

    anchor_q <- as.numeric(
      stats::quantile(
        valid_df$anchor,
        c(
          0.025,
          0.50,
          0.975
        ),
        na.rm = TRUE,
        names = FALSE
      )
    )

    terminal_q <- as.numeric(
      stats::quantile(
        valid_df$terminal,
        c(
          0.025,
          0.50,
          0.975
        ),
        na.rm = TRUE,
        names = FALSE
      )
    )

    right_n_q <- as.numeric(
      stats::quantile(
        valid_df$right_n,
        c(
          0.025,
          0.50,
          0.975
        ),
        na.rm = TRUE,
        names = FALSE
      )
    )

    median_jaccard <- safe_median(
      valid_df$jaccard
    )

  } else {

    ref_recovery <-
      NA_real_

    anchor_recovery <-
      NA_real_

    terminal_recovery <-
      NA_real_

    ref_iqr <-
      NA_real_

    anchor_iqr <-
      NA_real_

    terminal_iqr <-
      NA_real_

    right_fraction_iqr <-
      NA_real_

    ref_q <-
      rep(
        NA_real_,
        3L
      )

    anchor_q <-
      rep(
        NA_real_,
        3L
      )

    terminal_q <-
      rep(
        NA_real_,
        3L
      )

    right_n_q <-
      rep(
        NA_real_,
        3L
      )

    median_jaccard <-
      NA_real_
  }

  pass_valid <- (
    is.finite(
      valid_rate
    ) &&
    valid_rate >=
      MIN_BOOTSTRAP_VALID_RATE
  )

  pass_ref_recovery <- (
    is.finite(
      ref_recovery
    ) &&
    ref_recovery >=
      MIN_POSITION_RECOVERY_RATE
  )

  pass_anchor_recovery <- (
    is.finite(
      anchor_recovery
    ) &&
    anchor_recovery >=
      MIN_POSITION_RECOVERY_RATE
  )

  pass_terminal_recovery <- (
    is.finite(
      terminal_recovery
    ) &&
    terminal_recovery >=
      MIN_POSITION_RECOVERY_RATE
  )

  pass_ref_iqr <- (
    is.finite(
      ref_iqr
    ) &&
    ref_iqr /
      N <=
      MAX_POSITION_IQR_FRACTION
  )

  pass_anchor_iqr <- (
    is.finite(
      anchor_iqr
    ) &&
    anchor_iqr /
      N <=
      MAX_POSITION_IQR_FRACTION
  )

  pass_terminal_iqr <- (
    is.finite(
      terminal_iqr
    ) &&
    terminal_iqr /
      N <=
      MAX_POSITION_IQR_FRACTION
  )

  pass_right_iqr <- (
    is.finite(
      right_fraction_iqr
    ) &&
    right_fraction_iqr <=
      MAX_RIGHT_FRACTION_IQR
  )

  pass_stability <- all(
    c(
      pass_valid,
      pass_ref_recovery,
      pass_anchor_recovery,
      pass_terminal_recovery,
      pass_ref_iqr,
      pass_anchor_iqr,
      pass_terminal_iqr,
      pass_right_iqr
    )
  )

  stability_summary <- data.frame(

    Status =
      ifelse(
        pass_stability,
        "PASS",
        "UNSTABLE"
      ),

    FullRef =
      full_interval$ref,

    FullAnchor =
      full_interval$anchor,

    FullTerminal =
      full_interval$terminal,

    FullRightN =
      N -
      full_interval$anchor +
      1L,

    FullRightFraction =
      (
        N -
        full_interval$anchor +
        1L
      ) /
      N,

    BootstrapN =
      bootstrap_n,

    ValidBootstraps =
      valid_n,

    BootstrapValidRate =
      valid_rate,

    PositionToleranceN =
      tolerance_n,

    RefRecoveryRate =
      ref_recovery,

    AnchorRecoveryRate =
      anchor_recovery,

    TerminalRecoveryRate =
      terminal_recovery,

    RefIQR =
      ref_iqr,

    RefIQRFraction =
      ref_iqr /
      N,

    AnchorIQR =
      anchor_iqr,

    AnchorIQRFraction =
      anchor_iqr /
      N,

    TerminalIQR =
      terminal_iqr,

    TerminalIQRFraction =
      terminal_iqr /
      N,

    RightFractionIQR =
      right_fraction_iqr,

    RefCI025 =
      ref_q[1L],

    RefMedian =
      ref_q[2L],

    RefCI975 =
      ref_q[3L],

    AnchorCI025 =
      anchor_q[1L],

    AnchorMedian =
      anchor_q[2L],

    AnchorCI975 =
      anchor_q[3L],

    TerminalCI025 =
      terminal_q[1L],

    TerminalMedian =
      terminal_q[2L],

    TerminalCI975 =
      terminal_q[3L],

    RightNCI025 =
      right_n_q[1L],

    RightNMedian =
      right_n_q[2L],

    RightNCI975 =
      right_n_q[3L],

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
      pass_right_iqr,

    PassStability =
      pass_stability,

    stringsAsFactors = FALSE
  )

  list(

    bootstrap_df =
      bootstrap_df,

    stability_summary =
      stability_summary
  )
}

# =============================================================================
# FIGURE HELPERS
# =============================================================================

make_display_event_df <- function(
    event_df,
    total_n) {

  offset_big <- max(
    12L,
    round(
      total_n *
      0.006
    )
  )

  event_df %>%
    mutate(

      rank_display =
        case_when(

          event == "Anchor" ~
            rank -
            offset_big,

          event == "Ref" ~
            rank,

          event == "Terminal" ~
            rank +
            offset_big,

          TRUE ~
            rank
        )
    )
}


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
            3,

          ncol =
            1,

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
            1
        )
    )
  }

  dev.off()
}


compute_text_positions <- function(
    total_n,
    anchor,
    ref,
    terminal) {

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

    stat_x <-
      left_x +
      60L
  }

  list(
    left_x = left_x,
    stat_x = stat_x
  )
}

# =============================================================================
# FIGURE BUILDER
# =============================================================================

build_main_figure <- function(
    comparison_name,
    arm_name,
    variance_df,
    feature_df,
    interval_info,
    region_summary,
    detector,
    stability_summary,
    out_file,
    fig_tag,
    rank_tag,
    metric_tag) {

  total_n <- nrow(
    feature_df
  )

  anchor <-
    interval_info$anchor

  ref <-
    interval_info$ref

  terminal <-
    interval_info$terminal

  left_n <-
    region_summary$left_n

  left_min <-
    anchor -
    left_n

  left_max <-
    anchor -
    1L

  right_min <-
    anchor

  right_max <-
    total_n

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
      c(
        variance_df$smooth_log1p_empirical_variance[
          anchor
        ],
        variance_df$smooth_log1p_empirical_variance[
          ref
        ],
        variance_df$smooth_log1p_empirical_variance[
          terminal
        ]
      ),

    stringsAsFactors = FALSE
  )

  event_display_df <- make_display_event_df(
    event_df,
    total_n
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
    anchor,
    ref,
    terminal
  )

  top_y <- max(
    variance_df$smooth_log1p_empirical_variance,
    na.rm = TRUE
  )

  mid_y <- max(
    nb_long$value,
    na.rm = TRUE
  )

  box1 <- paste(
    fig_tag,
    rank_tag,
    metric_tag,
    "Blue = LEFT",
    "Green = RIGHT",
    "Grey = transition interval",
    sep = "\n"
  )

  box2 <- paste0(

    "Anchor = ",
    anchor,
    "\n",

    "Data-derived Ref = ",
    ref,
    "\n",

    "Terminal = ",
    terminal,
    "\n",

    "Interval = [",
    anchor,
    ", ",
    terminal,
    "]\n",

    "RIGHT n = ",
    region_summary$right_n,
    "\n",

    "Scale Ref IQR/N = ",
    round(
      detector$scale_iqr_fraction,
      4
    ),
    "\n",

    "Bootstrap status = ",
    stability_summary$Status[
      1L
    ]
  )

  box3 <- paste(

    "RIGHT = Anchor to end",

    "LEFT = equal matched block",

    "NB metrics do not select geometry",

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
      color = COL$var_curve,
      linewidth = 1.0
    ) +

    geom_vline(

      data =
        vline_df,

      aes(
        xintercept = xint,
        color = event,
        linetype = event
      ),

      linewidth = 0.9
    ) +

    geom_point(

      data =
        event_display_df,

      aes(
        rank_display,
        y,
        color = event,
        shape = event
      ),

      size = 3.4,
      stroke = 1.0
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

      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.25,

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

      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.25,

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
        EVENT_LEVELS,

      guide =
        guide_legend(

          override.aes =
            list(

              shape =
                unname(
                  EVENT_SHAPES
                ),

              linetype =
                unname(
                  EVENT_LTY
                ),

              linewidth =
                1.0,

              size =
                3.4
            )
        )
    ) +

    scale_shape_manual(
      values = EVENT_SHAPES,
      guide = "none"
    ) +

    scale_linetype_manual(
      values = EVENT_LTY,
      guide = "none"
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
          "Multiscale terminal change point ",
          "-> original Anchor/Ref/Terminal geometry"
        ),

      x =
        "Rank (ascending absolute PC1 loading)",

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
        xintercept = xint
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
        mid_y *
        0.96,

      label =
        box3,

      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.25,

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

      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.25,

      fill =
        grDevices::adjustcolor(
          "white",
          alpha.f = 0.96
        )
    ) +

    scale_color_manual(
      values = TRACE_COLORS,
      breaks = TRACE_LEVELS
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
          "RIGHT is compared directly with the matched LEFT ",
          "block after geometry is fixed"
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
        x = LEFT,
        xend = RIGHT,
        yend = metric
      ),

      color =
        "#7A7A7A",

      linewidth =
        0.8
    ) +

    geom_point(

      aes(
        x = LEFT,
        color = "LEFT"
      ),

      size =
        3.4
    ) +

    geom_point(

      aes(
        x = RIGHT,
        color = "RIGHT"
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
        paste0(

          "Bootstrap Anchor 95% interval: ",

          round(
            stability_summary$AnchorCI025[
              1L
            ]
          ),

          "-",

          round(
            stability_summary$AnchorCI975[
              1L
            ]
          ),

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
    ": detecting data-derived terminal transition..."
  )

  geometry <- compute_geometric_selection(
    rank_matrix,
    metric_matrix
  )

  if (
    !is.null(
      geometry$detector$scale_df
    )
  ) {

    write.csv(

      geometry$detector$scale_df,

      file.path(
        output_dir,
        paste0(
          "Table_TransitionScales_",
          prefix,
          ".csv"
        )
      ),

      row.names = FALSE
    )
  }

  if (
    !isTRUE(
      geometry$valid
    )
  ) {

    fail <- data.frame(

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
        "NO_VALID_GEOMETRY",

      Reason =
        geometry$reason,

      Anchor =
        NA_integer_,

      Ref =
        if (
          !is.null(
            geometry$detector$reference_rank
          )
        ) {
          geometry$detector$reference_rank
        } else {
          NA_integer_
        },

      Terminal =
        NA_integer_,

      LeftSize =
        NA_integer_,

      RightSize =
        NA_integer_,

      stringsAsFactors = FALSE
    )

    write.csv(

      fail,

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
      " | NO_VALID_GEOMETRY | ",
      geometry$reason
    )

    return(
      list(
        summary = fail,
        stability_summary = data.frame()
      )
    )
  }

  interval_info <-
    geometry$interval_info

  total_n <-
    geometry$total_n

  feature_df <- compute_ranked_feature_metrics(
    metric_matrix,
    geometry$rank_order
  )

  feature_df$abs_pc1_loading <-
    geometry$abs_loadings[
      geometry$rank_order
    ]

  region_summary <- summarize_regions(
    feature_df,
    interval_info$anchor,
    total_n
  )

  validate_method_level(

    feature_df =
      feature_df,

    interval_info =
      interval_info,

    region_summary =
      region_summary,

    total_n =
      total_n
  )

  message(

    "[",
    track,
    "] ",
    comparison_name,
    " ",
    arm_name,

    " | full-data Anchor=",
    interval_info$anchor,

    " Ref=",
    interval_info$ref,

    " Terminal=",
    interval_info$terminal,

    " RIGHT n=",
    region_summary$right_n,

    " | bootstrapping ",
    BOOTSTRAP_N,
    " replicates..."
  )

  stability <- assess_geometry_stability(

    rank_matrix =
      rank_matrix,

    metric_matrix =
      metric_matrix,

    full_geometry =
      geometry,

    bootstrap_n =
      BOOTSTRAP_N,

    seed =
      seed,

    cores =
      BOOTSTRAP_CORES
  )

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

    stability$stability_summary
  )

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
      stability_summary$Status[
        1L
      ],

    Anchor =
      interval_info$anchor,

    Ref =
      interval_info$ref,

    Terminal =
      interval_info$terminal,

    IntMin =
      interval_info$interval_min,

    IntMax =
      interval_info$interval_max,

    TransitionValidScales =
      geometry$detector$valid_scales,

    TransitionTotalScales =
      geometry$detector$total_scales,

    TransitionScaleRefIQR =
      geometry$detector$scale_iqr,

    TransitionScaleRefIQRFraction =
      geometry$detector$scale_iqr_fraction,

    TransitionMedianBICGain =
      geometry$detector$median_BIC_gain,

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

      BootstrapValidRate =
        stability_summary$BootstrapValidRate[
          1L
        ],

      RefRecoveryRate =
        stability_summary$RefRecoveryRate[
          1L
        ],

      AnchorRecoveryRate =
        stability_summary$AnchorRecoveryRate[
          1L
        ],

      TerminalRecoveryRate =
        stability_summary$TerminalRecoveryRate[
          1L
        ],

      AnchorCI025 =
        stability_summary$AnchorCI025[
          1L
        ],

      AnchorMedian =
        stability_summary$AnchorMedian[
          1L
        ],

      AnchorCI975 =
        stability_summary$AnchorCI975[
          1L
        ],

      RightNCI025 =
        stability_summary$RightNCI025[
          1L
        ],

      RightNMedian =
        stability_summary$RightNMedian[
          1L
        ],

      RightNCI975 =
        stability_summary$RightNCI975[
          1L
        ],

      MedianJaccard =
        stability_summary$MedianJaccard[
          1L
        ],

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

  boot_path <- file.path(
    output_dir,
    paste0(
      "Table_Bootstrap_",
      prefix,
      ".csv"
    )
  )

  stability_path <- file.path(
    output_dir,
    paste0(
      "Table_StabilitySummary_",
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

  write.csv(
    geometry$zero_df,
    zero_path,
    row.names = FALSE
  )

  write.csv(
    stability$bootstrap_df,
    boot_path,
    row.names = FALSE
  )

  write.csv(
    stability_summary,
    stability_path,
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
        stability_summary$Status[
          1L
        ],

      Anchor_lt_Ref =
        interval_info$anchor <
        interval_info$ref,

      Ref_lt_Terminal =
        interval_info$ref <
        interval_info$terminal,

      LeftN =
        region_summary$left_n,

      RightN =
        region_summary$right_n,

      EqualRegionSize =
        region_summary$left_n ==
        region_summary$right_n,

      PassStability =
        stability_summary$PassStability[
          1L
        ],

      stringsAsFactors = FALSE
    ),

    valid_path,

    row.names = FALSE
  )

  write.csv(

    feature_df %>%
      left_join(
        geometry$variance_df,
        by = "rank"
      ),

    rank_path,

    row.names = FALSE
  )

  write.csv(
    cutoff_summary,
    cut_path,
    row.names = FALSE
  )

  build_main_figure(

    comparison_name =
      comparison_name,

    arm_name =
      arm_name,

    variance_df =
      geometry$variance_df,

    feature_df =
      feature_df,

    interval_info =
      interval_info,

    region_summary =
      region_summary,

    detector =
      geometry$detector,

    stability_summary =
      stability_summary,

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

    " | ",
    stability_summary$Status[
      1L
    ],

    " | Anchor=",
    interval_info$anchor,

    " Ref=",
    interval_info$ref,

    " Terminal=",
    interval_info$terminal,

    " | RIGHT n=",
    region_summary$right_n,

    " | bootstrap Anchor median=",
    round(
      stability_summary$AnchorMedian[
        1L
      ]
    ),

    " | 95% Anchor interval=",
    round(
      stability_summary$AnchorCI025[
        1L
      ]
    ),

    "-",

    round(
      stability_summary$AnchorCI975[
        1L
      ]
    ),

    " | Jaccard(descriptive)=",

    sprintf(
      "%.3f",
      stability_summary$MedianJaccard[
        1L
      ]
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
  "Fixed 5,000-feature reference: REMOVED"
)

message(
  "Terminal detector: multiscale continuous segmented regression"
)

message(
  "Transition smoothing scales: ",
  paste(
    TRANSITION_SPARS,
    collapse = ", "
  )
)

message(
  "Minimum change-point segment guardrail: max(",
  MIN_TRANSITION_SEGMENT_ABSOLUTE,
  " features, ",
  100 *
  MIN_TRANSITION_SEGMENT_FRACTION,
  "% of N)"
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

    if (
      nrow(
        main_res$stability_summary
      ) > 0L
    ) {

      overall_stability_rows[[
        length(overall_stability_rows) +
        1L
      ]] <- main_res$stability_summary
    }

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
          "[DESeq2] skipped: package not available for ",
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
            deseq2_obj$rank_method_used,

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

        if (
          nrow(
            deseq2_res$stability_summary
          ) > 0L
        ) {

          overall_stability_rows[[
            length(overall_stability_rows) +
            1L
          ]] <- deseq2_res$stability_summary
        }
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
  "The terminal reference was estimated from variance geometry only."
)

message(
  "Anchor/Ref/Terminal and RIGHT=Anchor->end preserve the original method."
)

message(
  "NB2-related quantities were evaluated only after geometric selection."
)

message(
  "Outputs written to: ",
  OUT_ROOT
)

message(
  "============================================================"
)
