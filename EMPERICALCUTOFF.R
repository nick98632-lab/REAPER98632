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
# FULL-AXIS SEARCH FOR A STABLE NB1-LIKE -> NB2-LIKE LEADING-EDGE CUTOFF
# =============================================================================
#
# STUDY LOGIC
# -----------
# Features are ranked from LOWEST to HIGHEST absolute PC1 loading.
#
# The terminal RIGHT side therefore represents the high-|PC1|-loading
# leading edge.
#
# The previous fixed 5,000-feature reference is REMOVED.
#
# Instead, every admissible integer cutoff across the ranked axis is examined.
#
# For candidate cutoff c:
#
#     RIGHT(c) = c, ..., N
#
#     LEFT(c)  = the immediately preceding block containing exactly the
#                same number of features as RIGHT(c)
#
# The goal is to identify the cutoff at which the terminal RIGHT population
# shows the strongest reproducible shift toward NB2-like excess-variance
# behavior relative to the matched LEFT population.
#
# ---------------------------------------------------------------------------
# FEATURE-LEVEL DESCRIPTIVE QUANTITIES
# ---------------------------------------------------------------------------
#
# For each feature:
#
#     mu       = empirical mean
#     variance = empirical variance
#
#     extra = max(variance - mu, 0)
#
#     NB2 =
#         log(1 + extra)
#
#     NB2-NB1 =
#         log(1 + extra) - log(1 + mu)
#
#     alpha =
#         max((variance - mu) / mu^2, 0)
#
#     alpha*mu =
#         log(1 + alpha*mu)
#
# NB2-NB1 is used as the PRIMARY cutoff-selection signal because it directly
# measures the excess-variance term relative to the mean term.
#
# For every candidate cutoff, the program calculates:
#
#     mean(NB2-NB1 RIGHT) - mean(NB2-NB1 LEFT)
#
# and standardizes this difference using the estimated standard error.
#
# The resulting standardized separation score is an OPTIMIZATION SCORE.
# It is not interpreted as a formal p-value.
#
# NB2 and alpha*mu must also move in the expected direction:
#
#     RIGHT > LEFT
#
# ---------------------------------------------------------------------------
# FULL-AXIS SEARCH
# ---------------------------------------------------------------------------
#
# The program scans EVERY admissible integer cutoff.
#
# It does NOT:
#
#     - use 5,000 as a reference;
#     - use N_eff;
#     - restrict the search to 5%, 7.5%, 10%, 12.5%, or 15%;
#     - restrict selection to a predetermined location.
#
# A 5% minimum terminal-tail size is retained solely as an anti-degeneracy
# guardrail so that an extreme solution involving only a tiny number of
# terminal sites cannot win.
#
# For approximately 31,756 features:
#
#     5% ~= 1,588 features
#
# Thus this guardrail does NOT encode the previous 5,000-feature result.
#
# The maximum possible RIGHT size is determined automatically by the
# requirement that an equal-sized LEFT block must exist.
#
# ---------------------------------------------------------------------------
# CUTOFF SELECTION
# ---------------------------------------------------------------------------
#
# The standardized NB2-NB1 separation score is calculated at every admissible
# cutoff and smoothed across the cutoff axis.
#
# The final full-data cutoff is the integer rank with the maximum smoothed
# separation score among locations where:
#
#     delta NB2-NB1 > 0
#     delta NB2     > 0
#     delta alpha*mu > 0
#
# Therefore, if the biological transition truly occurs with approximately
# 5,000 features remaining in the leading edge, that scale should emerge
# naturally from the data rather than being supplied to the algorithm.
#
# ---------------------------------------------------------------------------
# BOOTSTRAP STABILITY
# ---------------------------------------------------------------------------
#
# Samples are resampled with replacement within each arm.
#
# EACH bootstrap replicate repeats the complete procedure:
#
#     bootstrap samples
#          ->
#     recompute PC1
#          ->
#     rank all features by |PC1 loading|
#          ->
#     recompute feature means/variances
#          ->
#     recompute NB-related quantities
#          ->
#     scan every admissible cutoff
#          ->
#     select optimal cutoff
#
# The full-data optimum remains the primary point estimate.
#
# Bootstrap output provides:
#
#     - bootstrap median cutoff;
#     - 95% percentile interval for cutoff;
#     - bootstrap median terminal-tail size;
#     - 95% percentile interval for terminal-tail size;
#     - cutoff IQR;
#     - cutoff recovery rate;
#     - descriptive Jaccard overlap.
#
# Stability requires:
#
#     >= 90% valid bootstrap replicates
#
#     >= 80% of valid replicates within 3% of N
#             of the full-data optimum
#
#     cutoff IQR <= 3% of N
#
# Jaccard overlap is descriptive only and does NOT determine PASS/UNSTABLE.
#
# ---------------------------------------------------------------------------
# VARIANCE GEOMETRY
# ---------------------------------------------------------------------------
#
# The original smoothed variance curve and second-derivative zero crossings
# are retained.
#
# They do NOT determine the cutoff.
#
# Instead, they provide geometric corroboration by reporting the nearest
# curvature zero crossing to the empirically optimized NB transition.
#
# ---------------------------------------------------------------------------
# IMPORTANT MANUSCRIPT INTERPRETATION
# ---------------------------------------------------------------------------
#
# Because NB2-related quantities are used to SELECT the cutoff here, those
# same quantities are not independent confirmatory evidence for that cutoff.
#
# This analysis establishes the empirically optimal and bootstrap-stable
# NB1-like -> NB2-like transition along the PC1-ranked axis.
#
# Formal downstream biological inference should be treated separately.
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <-
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <-
  "/root/REAPER98632/exports/manuscript_full_axis_nb_transition"


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


# -----------------------------------------------------------------------------
# Anti-degeneracy tail guardrail.
#
# This does NOT encode a 5,000-feature target.
# -----------------------------------------------------------------------------

MIN_TAIL_FRACTION <- 0.05

MIN_TAIL_ABSOLUTE <- 100L


# -----------------------------------------------------------------------------
# Cutoff-response smoothing.
#
# Prevents one individual rank from winning because of a microscopic
# rank-level fluctuation.
# -----------------------------------------------------------------------------

CUTOFF_SCORE_SPLINE_SPAR <- 0.60


# -----------------------------------------------------------------------------
# Original variance-curve geometry.
# -----------------------------------------------------------------------------

VAR_SPLINE_SPAR <- 0.60

DENSE_GRID_MULTIPLIER <- 2L

DENSE_GRID_MIN <- 5000L


# -----------------------------------------------------------------------------
# Bootstrap.
# -----------------------------------------------------------------------------

BOOTSTRAP_N <- 500L

BOOTSTRAP_SEED_BASE <- 20260820L


MIN_BOOTSTRAP_VALID_RATE <- 0.90

CUTOFF_TOLERANCE_FRACTION <- 0.03

MIN_CUTOFF_RECOVERY_RATE <- 0.80

MAX_CUTOFF_IQR_FRACTION <- 0.03


# -----------------------------------------------------------------------------
# DESeq2 supplementary track.
# -----------------------------------------------------------------------------

RUN_DESEQ2_SUPPLEMENT <- TRUE

DESEQ2_RANK_METHOD <- "normalized_log1p"

# Other allowed option:
#
#     "vst"


# -----------------------------------------------------------------------------
# Figures.
# -----------------------------------------------------------------------------

PNG_WIDTH_IN <- 14

PNG_HEIGHT_IN <- 11

PNG_DPI <- 260


# -----------------------------------------------------------------------------
# Cores.
# -----------------------------------------------------------------------------

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


dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# INPUT
# =============================================================================

read_count_matrix <- function(
    path,
    comparisons) {


  if (
    !file.exists(path)
  ) {

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
  # Identify experimental sample columns from COMPARISONS.
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
      "No sample columns matched the patterns in COMPARISONS."
    )
  }


  if (
    1L %in% sample_idx
  ) {

    stop(
      "Column 1 matched a sample pattern; ",
      "column 1 must contain feature identifiers."
    )
  }


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


  colnames(count_mat) <-
    colnames(count_df)


  storage.mode(count_mat) <-
    "numeric"


  bad_n <- sum(
    !is.finite(count_mat)
  )


  if (
    bad_n > 0L
  ) {

    message(

      "Replacing ",

      bad_n,

      " non-finite count entries with 0 ",
      "(same behavior used by the previous working script)."
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
  # Feature IDs.
  #
  # This preserves compatibility with the file that previously failed because
  # of a blank feature-ID field.
  # ---------------------------------------------------------------------------

  feature_ids <- trimws(

    as.character(
      raw_df[[1L]]
    )
  )


  blank_ids <- (

    is.na(feature_ids) |

    feature_ids == ""
  )


  if (
    any(blank_ids)
  ) {

    feature_ids[
      blank_ids
    ] <- paste0(

      "__feature_row_",

      which(
        blank_ids
      )
    )


    message(

      "Assigned deterministic row IDs to ",

      sum(
        blank_ids
      ),

      " blank feature identifier(s)."
    )
  }


  feature_ids <- make.unique(

    feature_ids,

    sep = "__dup_"
  )


  rownames(count_mat) <-
    feature_ids


  # ---------------------------------------------------------------------------
  # Remove features that contain zero counts in every selected sample.
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
      "Fewer than two nonzero features remain after filtering."
    )
  }


  count_mat
}


# =============================================================================
# NORMALIZATION
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


  log1p(
    cpm
  )
}


compute_deseq2_matrices <- function(
    count_mat_arm,
    rank_method = "normalized_log1p") {


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

    return(
      NULL
    )
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
          "DESeq2 default size factors failed; trying type='poscounts'."
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
# NUMERICAL HELPERS
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


safe_median <- function(
    x) {


  x <- x[
    is.finite(x)
  ]


  if (
    length(x) == 0L
  ) {

    return(
      NA_real_
    )
  }


  median(
    x
  )
}


safe_iqr <- function(
    x) {


  x <- x[
    is.finite(x)
  ]


  if (
    length(x) == 0L
  ) {

    return(
      NA_real_
    )
  }


  stats::IQR(
    x
  )
}


jaccard_similarity <- function(
    a,
    b) {


  a <- unique(
    a
  )


  b <- unique(
    b
  )


  union_set <- union(
    a,
    b
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
      a,
      b
    )
  ) /
    length(
      union_set
    )
}


# =============================================================================
# PC1 RANKING
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
  # Rows    = samples
  # Columns = features
  #
  # Sample-space PCA is used because the number of samples is much smaller
  # than the number of features.
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


  names(loading) <-
    colnames(
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
# FEATURE-LEVEL NB-RELATED QUANTITIES
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
        alpha_hat *
        mu
      ),

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# FAST WINDOW MOMENTS
# =============================================================================

window_moments <- function(
    x,
    starts,
    ends,
    n_vec) {


  # ---------------------------------------------------------------------------
  # Prefix arrays allow every candidate cutoff to be evaluated in O(N)
  # total time rather than recalculating every LEFT/RIGHT region separately.
  #
  # Range [a,b]:
  #
  #     prefix[b + 1] - prefix[a]
  # ---------------------------------------------------------------------------

  prefix <- c(
    0,
    cumsum(x)
  )


  prefix2 <- c(
    0,
    cumsum(
      x * x
    )
  )


  sums <-
    prefix[
      ends + 1L
    ] -
    prefix[
      starts
    ]


  sums2 <-
    prefix2[
      ends + 1L
    ] -
    prefix2[
      starts
    ]


  means <- sums /
    n_vec


  vars <- rep(
    0,
    length(means)
  )


  ok <- n_vec > 1L


  vars[
    ok
  ] <- (

    sums2[
      ok
    ] -
      (
        sums[
          ok
        ]^2 /
        n_vec[
          ok
        ]
      )

  ) /
    (
      n_vec[
        ok
      ] -
      1L
    )


  vars[
    !is.finite(vars)
  ] <- 0


  vars <- pmax(
    vars,
    0
  )


  list(

    mean =
      means,

    var =
      vars
  )
}


# =============================================================================
# FULL-AXIS CUTOFF SCAN
# =============================================================================

scan_all_cutoffs <- function(
    feature_df,
    min_tail_fraction = MIN_TAIL_FRACTION,
    min_tail_absolute = MIN_TAIL_ABSOLUTE,
    spline_spar = CUTOFF_SCORE_SPLINE_SPAR) {


  N <- nrow(
    feature_df
  )


  # ---------------------------------------------------------------------------
  # Minimum allowable terminal region.
  #
  # For N ~= 31,756 this is about 1,588 features.
  #
  # It prevents microscopic tails but does NOT steer the result toward 5,000.
  # ---------------------------------------------------------------------------

  min_tail_n <- max(

    as.integer(
      min_tail_absolute
    ),

    as.integer(

      ceiling(
        min_tail_fraction *
        N
      )
    ),

    2L
  )


  # ---------------------------------------------------------------------------
  # Equal LEFT/RIGHT regions require:
  #
  #     cutoff - 1 >= N - cutoff + 1
  #
  # therefore:
  #
  #     cutoff >= ceil((N + 2) / 2)
  # ---------------------------------------------------------------------------

  min_cutoff <- as.integer(

    ceiling(
      (
        N +
        2L
      ) /
        2
    )
  )


  max_cutoff <- as.integer(

    N -
    min_tail_n +
    1L
  )


  if (
    min_cutoff >
    max_cutoff
  ) {

    return(

      list(

        valid =
          FALSE,

        reason =
          "no_admissible_cutoffs",

        scan_df =
          data.frame(),

        selected =
          NULL
      )
    )
  }


  # ---------------------------------------------------------------------------
  # EVERY admissible integer cutoff.
  # ---------------------------------------------------------------------------

  cutoffs <- seq.int(

    min_cutoff,

    max_cutoff
  )


  right_n <- N -
    cutoffs +
    1L


  left_n <- right_n


  left_start <- cutoffs -
    right_n


  left_end <- cutoffs -
    1L


  right_start <- cutoffs


  right_end <- rep.int(

    N,

    length(
      cutoffs
    )
  )


  pref <- feature_df$NB2_NB1

  nb2 <- feature_df$NB2

  amu <- feature_df$alpha_mu


  if (
    any(
      !is.finite(pref)
    ) ||
    any(
      !is.finite(nb2)
    ) ||
    any(
      !is.finite(amu)
    )
  ) {

    return(

      list(

        valid =
          FALSE,

        reason =
          "nonfinite_nb_metric",

        scan_df =
          data.frame(),

        selected =
          NULL
      )
    )
  }


  # ---------------------------------------------------------------------------
  # LEFT/RIGHT moments for NB2-NB1.
  # ---------------------------------------------------------------------------

  L_pref <- window_moments(

    pref,

    left_start,

    left_end,

    left_n
  )


  R_pref <- window_moments(

    pref,

    right_start,

    right_end,

    right_n
  )


  # ---------------------------------------------------------------------------
  # LEFT/RIGHT moments for NB2.
  # ---------------------------------------------------------------------------

  L_nb2 <- window_moments(

    nb2,

    left_start,

    left_end,

    left_n
  )


  R_nb2 <- window_moments(

    nb2,

    right_start,

    right_end,

    right_n
  )


  # ---------------------------------------------------------------------------
  # LEFT/RIGHT moments for alpha*mu.
  # ---------------------------------------------------------------------------

  L_amu <- window_moments(

    amu,

    left_start,

    left_end,

    left_n
  )


  R_amu <- window_moments(

    amu,

    right_start,

    right_end,

    right_n
  )


  # ---------------------------------------------------------------------------
  # Directional contrasts.
  # ---------------------------------------------------------------------------

  delta_pref <-
    R_pref$mean -
    L_pref$mean


  delta_nb2 <-
    R_nb2$mean -
    L_nb2$mean


  delta_amu <-
    R_amu$mean -
    L_amu$mean


  # ---------------------------------------------------------------------------
  # Standardized NB2-NB1 separation score.
  #
  # This behaves like a standardized two-group separation statistic but is
  # used ONLY for optimization. It is not interpreted as an inferential
  # p-value because features need not be statistically independent.
  # ---------------------------------------------------------------------------

  se_pref <- sqrt(

    (
      L_pref$var /
      left_n
    ) +

    (
      R_pref$var /
      right_n
    )
  )


  se_pref[

    !is.finite(se_pref) |

    se_pref <=
    .Machine$double.eps

  ] <- NA_real_


  raw_score <- delta_pref /
    se_pref


  scan_df <- data.frame(

    cutoff =
      cutoffs,

    tail_n =
      right_n,

    tail_fraction =
      right_n /
      N,

    left_start =
      left_start,

    left_end =
      left_end,

    right_start =
      right_start,

    right_end =
      right_end,


    left_mean_NB2_NB1 =
      L_pref$mean,

    right_mean_NB2_NB1 =
      R_pref$mean,

    delta_NB2_NB1 =
      delta_pref,


    left_mean_NB2 =
      L_nb2$mean,

    right_mean_NB2 =
      R_nb2$mean,

    delta_NB2 =
      delta_nb2,


    left_mean_alpha_mu =
      L_amu$mean,

    right_mean_alpha_mu =
      R_amu$mean,

    delta_alpha_mu =
      delta_amu,


    raw_separation_score =
      raw_score,

    stringsAsFactors = FALSE
  )


  finite_score <- is.finite(
    scan_df$raw_separation_score
  )


  if (
    sum(
      finite_score
    ) < 4L
  ) {

    return(

      list(

        valid =
          FALSE,

        reason =
          "insufficient_finite_cutoff_scores",

        scan_df =
          scan_df,

        selected =
          NULL
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Smooth the RESPONSE to cutoff across the entire axis.
  #
  # The spline is NOT centered on 5,000 or any other predefined tail size.
  # ---------------------------------------------------------------------------

  score_spline <- tryCatch(

    stats::smooth.spline(

      x =
        scan_df$cutoff[
          finite_score
        ],

      y =
        scan_df$raw_separation_score[
          finite_score
        ],

      spar =
        spline_spar
    ),

    error =
      function(e) {
        NULL
      }
  )


  if (
    is.null(
      score_spline
    )
  ) {

    scan_df$smoothed_separation_score <-
      scan_df$raw_separation_score

  } else {

    scan_df$smoothed_separation_score <-
      as.numeric(

        stats::predict(

          score_spline,

          x =
            scan_df$cutoff,

          deriv =
            0
        )$y
      )
  }


  # ---------------------------------------------------------------------------
  # The transition must consistently point toward greater NB2-like behavior
  # on the RIGHT.
  # ---------------------------------------------------------------------------

  eligible <- (

    is.finite(
      scan_df$smoothed_separation_score
    ) &

    scan_df$delta_NB2_NB1 > 0 &

    scan_df$delta_NB2 > 0 &

    scan_df$delta_alpha_mu > 0
  )


  scan_df$eligible <-
    eligible


  scan_df$selected <-
    FALSE


  if (
    !any(
      eligible
    )
  ) {

    return(

      list(

        valid =
          FALSE,

        reason =
          "no_cutoff_with_consistent_right_greater_than_left_direction",

        scan_df =
          scan_df,

        selected =
          NULL
      )
    )
  }


  eligible_idx <- which(
    eligible
  )


  best_idx <- eligible_idx[

    which.max(

      scan_df$smoothed_separation_score[
        eligible_idx
      ]
    )
  ]


  scan_df$selected[
    best_idx
  ] <- TRUE


  list(

    valid =
      TRUE,

    reason =
      NA_character_,

    scan_df =
      scan_df,

    selected =
      scan_df[
        best_idx,
        ,
        drop = FALSE
      ],

    min_tail_n =
      min_tail_n
  )
}


# =============================================================================
# ORIGINAL VARIANCE GEOMETRY
# =============================================================================

find_d2_zero_crossings <- function(
    dense_x,
    dense_d2) {


  ok <- (

    is.finite(
      dense_x
    ) &

    is.finite(
      dense_d2
    )
  )


  x <- dense_x[
    ok
  ]


  y <- dense_d2[
    ok
  ]


  if (
    length(x) < 2L
  ) {

    return(

      data.frame(

        crossing_rank =
          numeric(0),

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


  if (
    length(idx) > 0L
  ) {

    frac <- abs(
      a[
        idx
      ]
    ) /
      (
        abs(
          a[
            idx
          ]
        ) +
        abs(
          b[
            idx
          ]
        )
      )


    crossings <-

      xa[
        idx
      ] +

      frac *

      (
        xb[
          idx
        ] -
        xa[
          idx
        ]
      )
  }


  exact_idx <- which(
    y == 0
  )


  if (
    length(exact_idx) > 0L
  ) {

    crossings <- c(

      crossings,

      x[
        exact_idx
      ]
    )
  }


  crossings <- sort(

    unique(

      round(

        crossings[
          is.finite(
            crossings
          )
        ],

        8L
      )
    )
  )


  data.frame(

    crossing_rank =
      crossings,

    stringsAsFactors = FALSE
  )
}


compute_variance_geometry <- function(
    feature_df,
    spar = VAR_SPLINE_SPAR) {


  ranks <- feature_df$rank


  log_var <- log1p(
    feature_df$empirical_variance
  )


  fit <- stats::smooth.spline(

    x =
      ranks,

    y =
      log_var,

    spar =
      spar
  )


  smooth_y <- as.numeric(

    stats::predict(

      fit,

      x =
        ranks,

      deriv =
        0
    )$y
  )


  smooth_d2 <- as.numeric(

    stats::predict(

      fit,

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

      fit,

      x =
        dense_x,

      deriv =
        0
    )$y
  )


  dense_d2 <- as.numeric(

    stats::predict(

      fit,

      x =
        dense_x,

      deriv =
        2
    )$y
  )


  zero_df <- find_d2_zero_crossings(

    dense_x,

    dense_d2
  )


  list(

    variance_df =
      data.frame(

        rank =
          ranks,

        log1p_empirical_variance =
          log_var,

        smooth_log1p_empirical_variance =
          smooth_y,

        d2_spline =
          smooth_d2,

        stringsAsFactors = FALSE
      ),


    dense_df =
      data.frame(

        rank =
          dense_x,

        smooth_log1p_empirical_variance =
          dense_y,

        d2 =
          dense_d2,

        stringsAsFactors = FALSE
      ),


    zero_df =
      zero_df
  )
}


nearest_d2_crossing <- function(
    zero_df,
    cutoff) {


  if (
    nrow(zero_df) == 0L
  ) {

    return(

      list(

        rank =
          NA_real_,

        distance =
          NA_real_
      )
    )
  }


  idx <- which.min(

    abs(

      zero_df$crossing_rank -
      cutoff
    )
  )


  list(

    rank =
      zero_df$crossing_rank[
        idx
      ],

    distance =
      zero_df$crossing_rank[
        idx
      ] -
      cutoff
  )
}


# =============================================================================
# COMPLETE FULL-DATA ANALYSIS
# =============================================================================

analyze_full_track <- function(
    rank_matrix,
    metric_matrix) {


  if (
    nrow(rank_matrix) !=
    nrow(metric_matrix) ||
    ncol(rank_matrix) !=
    ncol(metric_matrix)
  ) {

    stop(
      "Rank and metric matrices have incompatible dimensions."
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


  feature_df <- compute_ranked_feature_metrics(

    metric_matrix,

    rank_order
  )


  feature_df$abs_pc1_loading <-

    abs_loadings[
      rank_order
    ]


  scan <- scan_all_cutoffs(
    feature_df
  )


  if (
    !isTRUE(
      scan$valid
    )
  ) {

    return(

      list(

        valid =
          FALSE,

        reason =
          scan$reason,

        abs_loadings =
          abs_loadings,

        rank_order =
          rank_order,

        feature_df =
          feature_df,

        scan =
          scan,

        geometry =
          NULL,

        nearest_crossing =
          list(
            rank = NA_real_,
            distance = NA_real_
          )
      )
    )
  }


  geometry <- compute_variance_geometry(
    feature_df
  )


  nearest_crossing <- nearest_d2_crossing(

    geometry$zero_df,

    scan$selected$cutoff[
      1L
    ]
  )


  list(

    valid =
      TRUE,

    reason =
      NA_character_,

    abs_loadings =
      abs_loadings,

    rank_order =
      rank_order,

    feature_df =
      feature_df,

    scan =
      scan,

    geometry =
      geometry,

    nearest_crossing =
      nearest_crossing
  )
}


# =============================================================================
# SELECTED LEFT/RIGHT SUMMARY
# =============================================================================

summarize_selected_cutoff <- function(
    feature_df,
    cutoff) {


  N <- nrow(
    feature_df
  )


  right_n <- N -
    cutoff +
    1L


  left_start <- cutoff -
    right_n


  left_end <- cutoff -
    1L


  if (
    left_start < 1L
  ) {

    stop(
      "Selected cutoff cannot form an equal-sized LEFT block."
    )
  }


  left_df <- feature_df[

    left_start:left_end,

    ,

    drop = FALSE
  ]


  right_df <- feature_df[

    cutoff:N,

    ,

    drop = FALSE
  ]


  data.frame(

    left_start =
      left_start,

    left_end =
      left_end,

    right_start =
      cutoff,

    right_end =
      N,

    left_n =
      nrow(
        left_df
      ),

    right_n =
      nrow(
        right_df
      ),


    left_median_NB2_NB1 =
      median(
        left_df$NB2_NB1,
        na.rm = TRUE
      ),

    right_median_NB2_NB1 =
      median(
        right_df$NB2_NB1,
        na.rm = TRUE
      ),

    delta_median_NB2_NB1 =
      median(
        right_df$NB2_NB1,
        na.rm = TRUE
      ) -
      median(
        left_df$NB2_NB1,
        na.rm = TRUE
      ),


    left_median_NB2 =
      median(
        left_df$NB2,
        na.rm = TRUE
      ),

    right_median_NB2 =
      median(
        right_df$NB2,
        na.rm = TRUE
      ),

    delta_median_NB2 =
      median(
        right_df$NB2,
        na.rm = TRUE
      ) -
      median(
        left_df$NB2,
        na.rm = TRUE
      ),


    left_median_alpha_mu =
      median(
        left_df$alpha_mu,
        na.rm = TRUE
      ),

    right_median_alpha_mu =
      median(
        right_df$alpha_mu,
        na.rm = TRUE
      ),

    delta_median_alpha_mu =
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
# BOOTSTRAP THE ENTIRE SEARCH
# =============================================================================

assess_cutoff_stability <- function(
    rank_matrix,
    metric_matrix,
    full_analysis,
    bootstrap_n = BOOTSTRAP_N,
    seed = BOOTSTRAP_SEED_BASE,
    cores = BOOTSTRAP_CORES) {


  N <- nrow(
    rank_matrix
  )


  sample_n <- ncol(
    rank_matrix
  )


  if (
    !isTRUE(
      full_analysis$valid
    )
  ) {

    stop(
      "Cannot bootstrap an invalid full-data analysis."
    )
  }


  full_cutoff <- as.integer(

    full_analysis$scan$selected$cutoff[
      1L
    ]
  )


  full_right_idx <-

    full_analysis$rank_order[

      full_cutoff:N
    ]


  full_right_features <-

    rownames(
      rank_matrix
    )[
      full_right_idx
    ]


  # ---------------------------------------------------------------------------
  # Generate all bootstrap draws before parallel execution so results remain
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


  worker <- function(
      b) {


    idx <- bootstrap_indices[[
      b
    ]]


    invalid_row <- function(
        reason) {


      data.frame(

        bootstrap =
          b,

        valid =
          FALSE,

        invalid_reason =
          reason,

        cutoff =
          NA_integer_,

        tail_n =
          NA_integer_,

        tail_fraction =
          NA_real_,

        raw_score =
          NA_real_,

        smoothed_score =
          NA_real_,

        delta_NB2_NB1 =
          NA_real_,

        delta_NB2 =
          NA_real_,

        delta_alpha_mu =
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


    ans <- tryCatch(


      {

        abs_loadings <-
          compute_abs_pc1_loadings(
            boot_rank
          )


        rank_order <- order(

          abs_loadings,

          decreasing = FALSE
        )


        feature_df <-
          compute_ranked_feature_metrics(

            boot_metric,

            rank_order
          )


        scan <-
          scan_all_cutoffs(
            feature_df
          )


        if (
          !isTRUE(
            scan$valid
          )
        ) {

          invalid_row(
            scan$reason
          )

        } else {


          selected <-
            scan$selected


          cutoff <- as.integer(

            selected$cutoff[
              1L
            ]
          )


          right_idx <-

            rank_order[

              cutoff:N
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

            cutoff =
              cutoff,

            tail_n =
              as.integer(

                selected$tail_n[
                  1L
                ]
              ),

            tail_fraction =
              selected$tail_fraction[
                1L
              ],

            raw_score =
              selected$raw_separation_score[
                1L
              ],

            smoothed_score =
              selected$smoothed_separation_score[
                1L
              ],

            delta_NB2_NB1 =
              selected$delta_NB2_NB1[
                1L
              ],

            delta_NB2 =
              selected$delta_NB2[
                1L
              ],

            delta_alpha_mu =
              selected$delta_alpha_mu[
                1L
              ],

            jaccard =
              jaccard_similarity(

                full_right_features,

                right_features
              ),

            stringsAsFactors = FALSE
          )
        }
      },


      error =
        function(e) {

          invalid_row(

            paste0(

              "analysis_error: ",

              conditionMessage(e)
            )
          )
        }
    )


    ans
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


  valid_rate <- valid_n /
    bootstrap_n


  tolerance_n <- ceiling(

    CUTOFF_TOLERANCE_FRACTION *
    N
  )


  if (
    valid_n > 0L
  ) {

    cutoff_recovery_rate <- mean(

      abs(

        valid_df$cutoff -
        full_cutoff

      ) <=
        tolerance_n
    )


    cutoff_iqr <- safe_iqr(
      valid_df$cutoff
    )


    cutoff_iqr_fraction <-
      cutoff_iqr /
      N


    cutoff_ci <- as.numeric(

      stats::quantile(

        valid_df$cutoff,

        probs =
          c(
            0.025,
            0.50,
            0.975
          ),

        na.rm =
          TRUE,

        names =
          FALSE,

        type =
          7
      )
    )


    tail_ci <- as.numeric(

      stats::quantile(

        valid_df$tail_n,

        probs =
          c(
            0.025,
            0.50,
            0.975
          ),

        na.rm =
          TRUE,

        names =
          FALSE,

        type =
          7
      )
    )


    median_jaccard <-
      safe_median(
        valid_df$jaccard
      )


    median_score <-
      safe_median(
        valid_df$smoothed_score
      )

  } else {

    cutoff_recovery_rate <-
      NA_real_

    cutoff_iqr <-
      NA_real_

    cutoff_iqr_fraction <-
      NA_real_

    cutoff_ci <-
      c(
        NA_real_,
        NA_real_,
        NA_real_
      )

    tail_ci <-
      c(
        NA_real_,
        NA_real_,
        NA_real_
      )

    median_jaccard <-
      NA_real_

    median_score <-
      NA_real_
  }


  pass_valid <- (

    is.finite(
      valid_rate
    ) &&

    valid_rate >=
      MIN_BOOTSTRAP_VALID_RATE
  )


  pass_recovery <- (

    is.finite(
      cutoff_recovery_rate
    ) &&

    cutoff_recovery_rate >=
      MIN_CUTOFF_RECOVERY_RATE
  )


  pass_iqr <- (

    is.finite(
      cutoff_iqr_fraction
    ) &&

    cutoff_iqr_fraction <=
      MAX_CUTOFF_IQR_FRACTION
  )


  # Jaccard intentionally excluded from PASS/FAIL.

  pass_stability <- (

    pass_valid &&
    pass_recovery &&
    pass_iqr
  )


  stability_summary <- data.frame(

    Status =
      ifelse(
        pass_stability,
        "PASS",
        "UNSTABLE"
      ),


    FullCutoff =
      full_cutoff,

    FullTailN =
      N -
      full_cutoff +
      1L,

    FullTailFraction =
      (
        N -
        full_cutoff +
        1L
      ) /
      N,


    BootstrapN =
      bootstrap_n,

    ValidBootstraps =
      valid_n,

    BootstrapValidRate =
      valid_rate,


    CutoffToleranceN =
      tolerance_n,

    CutoffRecoveryRate =
      cutoff_recovery_rate,

    CutoffIQR =
      cutoff_iqr,

    CutoffIQRFraction =
      cutoff_iqr_fraction,


    BootstrapCutoffCI025 =
      cutoff_ci[
        1L
      ],

    BootstrapCutoffMedian =
      cutoff_ci[
        2L
      ],

    BootstrapCutoffCI975 =
      cutoff_ci[
        3L
      ],


    BootstrapTailNCI025 =
      tail_ci[
        1L
      ],

    BootstrapTailNMedian =
      tail_ci[
        2L
      ],

    BootstrapTailNCI975 =
      tail_ci[
        3L
      ],


    MedianBootstrapScore =
      median_score,

    MedianJaccard =
      median_jaccard,


    PassValidRate =
      pass_valid,

    PassCutoffRecovery =
      pass_recovery,

    PassCutoffIQR =
      pass_iqr,

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
# FIGURES
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
                1.05,
                1.05,
                0.90
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


build_figure <- function(
    comparison_name,
    arm_name,
    track,
    full_analysis,
    selected_summary,
    stability_summary,
    out_file) {


  cutoff <- as.integer(

    full_analysis$scan$selected$cutoff[
      1L
    ]
  )


  N <- nrow(
    full_analysis$feature_df
  )


  left_start <-
    selected_summary$left_start[
      1L
    ]


  left_end <-
    selected_summary$left_end[
      1L
    ]


  geom_df <-
    full_analysis$geometry$variance_df


  zero_df <-
    full_analysis$geometry$zero_df


  # ---------------------------------------------------------------------------
  # PANEL 1: variance geometry.
  # ---------------------------------------------------------------------------

  p1 <- ggplot(

    geom_df,

    aes(
      rank,
      smooth_log1p_empirical_variance
    )
  ) +

    annotate(

      "rect",

      xmin =
        left_start,

      xmax =
        left_end,

      ymin =
        -Inf,

      ymax =
        Inf,

      fill =
        "grey85",

      alpha =
        0.55
    ) +

    annotate(

      "rect",

      xmin =
        cutoff,

      xmax =
        N,

      ymin =
        -Inf,

      ymax =
        Inf,

      fill =
        "grey70",

      alpha =
        0.35
    ) +

    geom_line(
      linewidth = 0.9
    ) +

    geom_vline(

      xintercept =
        cutoff,

      linewidth =
        0.9
    )


  if (
    nrow(zero_df) > 0L
  ) {

    p1 <- p1 +

      geom_vline(

        data =
          zero_df,

        aes(
          xintercept =
            crossing_rank
        ),

        inherit.aes =
          FALSE,

        linewidth =
          0.25,

        linetype =
          "dotted",

        alpha =
          0.5
      )
  }


  p1 <- p1 +

    labs(

      title =
        paste0(

          comparison_name,
          " ",
          arm_name,
          " ",
          track,
          ": variance geometry"
        ),

      subtitle =
        paste0(

          "Selected cutoff = ",
          cutoff,

          "; terminal tail n = ",
          selected_summary$right_n[
            1L
          ],

          "; nearest d2 crossing = ",

          ifelse(

            is.finite(
              full_analysis$nearest_crossing$rank
            ),

            round(
              full_analysis$nearest_crossing$rank,
              1
            ),

            "NA"
          )
        ),

      x =
        "Rank (ascending |PC1 loading|)",

      y =
        "Smoothed log(1 + empirical variance)"
    ) +

    theme_bw(
      base_size = 11
    ) +

    theme(
      panel.grid.minor =
        element_blank()
    )


  # ---------------------------------------------------------------------------
  # PANEL 2: entire cutoff-response curve.
  # ---------------------------------------------------------------------------

  scan_df <-
    full_analysis$scan$scan_df


  p2 <- ggplot(

    scan_df,

    aes(
      cutoff,
      smoothed_separation_score
    )
  ) +

    geom_line(
      linewidth = 0.9
    ) +

    geom_vline(

      xintercept =
        cutoff,

      linewidth =
        0.9
    ) +

    geom_point(

      data =
        scan_df[
          scan_df$selected,
          ,
          drop = FALSE
        ],

      aes(
        cutoff,
        smoothed_separation_score
      ),

      size =
        2.8
    ) +

    labs(

      title =
        "Full-axis NB1-like to NB2-like cutoff scan",

      subtitle =
        paste0(

          "Every admissible cutoff scanned; no 5,000-feature reference supplied. ",

          "Bootstrap 95% cutoff interval: ",

          round(
            stability_summary$BootstrapCutoffCI025[
              1L
            ]
          ),

          "-",

          round(
            stability_summary$BootstrapCutoffCI975[
              1L
            ]
          )
        ),

      x =
        "Candidate cutoff rank",

      y =
        "Smoothed standardized NB2-NB1 separation score"
    ) +

    theme_bw(
      base_size = 11
    ) +

    theme(
      panel.grid.minor =
        element_blank()
    )


  # ---------------------------------------------------------------------------
  # PANEL 3: median LEFT vs RIGHT feature-level quantities.
  # ---------------------------------------------------------------------------

  region_df <- data.frame(

    metric =
      factor(

        c(
          "NB2-NB1",
          "NB2",
          "alpha*mu"
        ),

        levels =
          rev(

            c(
              "NB2-NB1",
              "NB2",
              "alpha*mu"
            )
          )
      ),


    LEFT =
      c(

        selected_summary$left_median_NB2_NB1[
          1L
        ],

        selected_summary$left_median_NB2[
          1L
        ],

        selected_summary$left_median_alpha_mu[
          1L
        ]
      ),


    RIGHT =
      c(

        selected_summary$right_median_NB2_NB1[
          1L
        ],

        selected_summary$right_median_NB2[
          1L
        ],

        selected_summary$right_median_alpha_mu[
          1L
        ]
      )
  )


  p3 <- ggplot(

    region_df,

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

      linewidth =
        0.8
    ) +

    geom_point(

      aes(

        x =
          LEFT,

        shape =
          "LEFT"
      ),

      size =
        3.2
    ) +

    geom_point(

      aes(

        x =
          RIGHT,

        shape =
          "RIGHT"
      ),

      size =
        3.2
    ) +

    scale_shape_manual(

      values =
        c(
          "LEFT" = 16,
          "RIGHT" = 17
        )
    ) +

    labs(

      title =
        "Matched LEFT versus terminal RIGHT region",

      subtitle =
        paste0(

          "Status: ",
          stability_summary$Status[
            1L
          ],

          "; cutoff recovery = ",

          sprintf(

            "%.3f",

            stability_summary$CutoffRecoveryRate[
              1L
            ]
          ),

          "; Jaccard = ",

          sprintf(

            "%.3f",

            stability_summary$MedianJaccard[
              1L
            ]
          )
        ),

      x =
        "Median feature-level diagnostic",

      y =
        NULL,

      shape =
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
# RUN ONE TRACK
# =============================================================================

run_one_track <- function(
    comparison_name,
    arm_name,
    rank_matrix,
    metric_matrix,
    output_dir,
    track,
    rank_method,
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

    ": scanning every admissible cutoff..."
  )


  full_analysis <- analyze_full_track(

    rank_matrix,

    metric_matrix
  )


  # ---------------------------------------------------------------------------
  # No valid full-data cutoff.
  # ---------------------------------------------------------------------------

  if (
    !isTRUE(
      full_analysis$valid
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
        rank_method,

      metric_matrix =
        metric_name,

      Status =
        "NO_VALID_CUTOFF",

      Reason =
        full_analysis$reason,

      stringsAsFactors = FALSE
    )


    write.csv(

      full_analysis$scan$scan_df,

      file.path(

        output_dir,

        paste0(
          "Table_CutoffScan_",
          prefix,
          ".csv"
        )
      ),

      row.names = FALSE
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

      " | NO VALID CUTOFF"
    )


    return(

      list(

        summary =
          fail,

        stability_summary =
          data.frame()
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Write complete full-axis scan.
  # ---------------------------------------------------------------------------

  write.csv(

    full_analysis$scan$scan_df,

    file.path(

      output_dir,

      paste0(
        "Table_CutoffScan_",
        prefix,
        ".csv"
      )
    ),

    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Write variance curvature transitions.
  # ---------------------------------------------------------------------------

  write.csv(

    full_analysis$geometry$zero_df,

    file.path(

      output_dir,

      paste0(
        "Table_ZeroCrossings_",
        prefix,
        ".csv"
      )
    ),

    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Ranked feature table.
  # ---------------------------------------------------------------------------

  write.csv(

    full_analysis$feature_df,

    file.path(

      output_dir,

      paste0(
        "Table_RankedFeatures_",
        prefix,
        ".csv"
      )
    ),

    row.names = FALSE
  )


  selected_cutoff <- as.integer(

    full_analysis$scan$selected$cutoff[
      1L
    ]
  )


  selected_summary <-
    summarize_selected_cutoff(

      full_analysis$feature_df,

      selected_cutoff
    )


  message(

    "[",

    track,

    "] ",

    comparison_name,

    " ",

    arm_name,

    ": full-data optimum cutoff=",

    selected_cutoff,

    ", terminal tail n=",

    selected_summary$right_n[
      1L
    ],

    "; bootstrapping ",

    BOOTSTRAP_N,

    " replicates..."
  )


  # ---------------------------------------------------------------------------
  # Bootstrap complete selection procedure.
  # ---------------------------------------------------------------------------

  stability <- assess_cutoff_stability(

    rank_matrix =
      rank_matrix,

    metric_matrix =
      metric_matrix,

    full_analysis =
      full_analysis,

    bootstrap_n =
      BOOTSTRAP_N,

    seed =
      seed,

    cores =
      BOOTSTRAP_CORES
  )


  write.csv(

    stability$bootstrap_df,

    file.path(

      output_dir,

      paste0(
        "Table_Bootstrap_",
        prefix,
        ".csv"
      )
    ),

    row.names = FALSE
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


  selected_row <-
    full_analysis$scan$selected


  # ---------------------------------------------------------------------------
  # Main cutoff summary.
  # ---------------------------------------------------------------------------

  summary_row <- data.frame(

    comp =
      comparison_name,

    arm =
      arm_name,

    track =
      track,

    rank_method =
      rank_method,

    metric_matrix =
      metric_name,


    Status =
      stability_summary$Status[
        1L
      ],


    Cutoff =
      selected_cutoff,

    TailN =
      selected_summary$right_n[
        1L
      ],

    TailFraction =
      selected_summary$right_n[
        1L
      ] /
      nrow(
        full_analysis$feature_df
      ),


    LeftStart =
      selected_summary$left_start[
        1L
      ],

    LeftEnd =
      selected_summary$left_end[
        1L
      ],

    RightStart =
      selected_summary$right_start[
        1L
      ],

    RightEnd =
      selected_summary$right_end[
        1L
      ],


    FullRawSeparationScore =
      selected_row$raw_separation_score[
        1L
      ],

    FullSmoothedSeparationScore =
      selected_row$smoothed_separation_score[
        1L
      ],


    MeanDelta_NB2_NB1 =
      selected_row$delta_NB2_NB1[
        1L
      ],

    MeanDelta_NB2 =
      selected_row$delta_NB2[
        1L
      ],

    MeanDelta_alpha_mu =
      selected_row$delta_alpha_mu[
        1L
      ],


    MedianDelta_NB2_NB1 =
      selected_summary$delta_median_NB2_NB1[
        1L
      ],

    MedianDelta_NB2 =
      selected_summary$delta_median_NB2[
        1L
      ],

    MedianDelta_alpha_mu =
      selected_summary$delta_median_alpha_mu[
        1L
      ],


    NearestD2Crossing =
      full_analysis$nearest_crossing$rank,

    DistanceToNearestD2 =
      full_analysis$nearest_crossing$distance,


    BootstrapValidRate =
      stability_summary$BootstrapValidRate[
        1L
      ],

    CutoffRecoveryRate =
      stability_summary$CutoffRecoveryRate[
        1L
      ],

    CutoffIQR =
      stability_summary$CutoffIQR[
        1L
      ],

    CutoffIQRFraction =
      stability_summary$CutoffIQRFraction[
        1L
      ],


    BootstrapCutoffCI025 =
      stability_summary$BootstrapCutoffCI025[
        1L
      ],

    BootstrapCutoffMedian =
      stability_summary$BootstrapCutoffMedian[
        1L
      ],

    BootstrapCutoffCI975 =
      stability_summary$BootstrapCutoffCI975[
        1L
      ],


    BootstrapTailNCI025 =
      stability_summary$BootstrapTailNCI025[
        1L
      ],

    BootstrapTailNMedian =
      stability_summary$BootstrapTailNMedian[
        1L
      ],

    BootstrapTailNCI975 =
      stability_summary$BootstrapTailNCI975[
        1L
      ],


    MedianJaccard =
      stability_summary$MedianJaccard[
        1L
      ],

    stringsAsFactors = FALSE
  )


  write.csv(

    summary_row,

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
  # Figure.
  # ---------------------------------------------------------------------------

  build_figure(

    comparison_name =
      comparison_name,

    arm_name =
      arm_name,

    track =
      track,

    full_analysis =
      full_analysis,

    selected_summary =
      selected_summary,

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
      )
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

    " | cutoff=",

    selected_cutoff,

    " | terminal tail n=",

    selected_summary$right_n[
      1L
    ],

    " | bootstrap median cutoff=",

    round(

      stability_summary$BootstrapCutoffMedian[
        1L
      ]
    ),

    " | 95% cutoff interval=",

    round(

      stability_summary$BootstrapCutoffCI025[
        1L
      ]
    ),

    "-",

    round(

      stability_summary$BootstrapCutoffCI975[
        1L
      ]
    ),

    " | 95% tail-size interval=",

    round(

      stability_summary$BootstrapTailNCI025[
        1L
      ]
    ),

    "-",

    round(

      stability_summary$BootstrapTailNCI975[
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
      summary_row,

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
  "Fixed 5,000-feature reference: REMOVED"
)


message(

  "Admissible terminal-tail guardrail: max(",

  MIN_TAIL_ABSOLUTE,

  " features, ",

  100 *
  MIN_TAIL_FRACTION,

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


all_summary_rows <- list()

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
      length(sample_idx) < 2L
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
    #
    # Ranking:
    #     CPM log1p
    #
    # NB metrics:
    #     raw counts
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

      rank_method =
        "CPM log1p",

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


    all_summary_rows[[
      length(
        all_summary_rows
      ) +
      1L
    ]] <- main_res$summary


    if (
      nrow(
        main_res$stability_summary
      ) > 0L
    ) {

      all_stability_rows[[
        length(
          all_stability_rows
        ) +
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

          rank_method =
            deseq2_obj$rank_method_used,

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


        all_summary_rows[[
          length(
            all_summary_rows
          ) +
          1L
        ]] <- deseq2_res$summary


        if (
          nrow(
            deseq2_res$stability_summary
          ) > 0L
        ) {

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
}


# =============================================================================
# OVERALL OUTPUTS
# =============================================================================

overall_summary <- bind_rows(
  all_summary_rows
)


overall_stability <- bind_rows(
  all_stability_rows
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
  "Every admissible integer cutoff was searched for every track."
)


message(
  "No 5,000-feature reference was supplied to the algorithm."
)


message(
  "N_eff was not used."
)


message(
  "The reported terminal-tail size is therefore an empirical result."
)


message(
  "The full-data optimum is the primary cutoff estimate."
)


message(
  "The bootstrap gives uncertainty and stability of that optimum."
)


message(
  "Jaccard is descriptive only."
)


message(
  "Outputs written to: ",
  OUT_ROOT
)


message(
  "============================================================"
)
