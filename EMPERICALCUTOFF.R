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
# FULL-RANK-AXIS SEARCH FOR THE NB1-LIKE -> NB2-LIKE LEADING-EDGE TRANSITION
# =============================================================================
#
# RATIONALE
# ---------
#
# Features are ranked from LOWEST to HIGHEST absolute PC1 loading.
#
# Thus:
#
#     LEFT  = relatively low-|PC1|-loading portion of the ranked series
#     RIGHT = high-|PC1|-loading terminal leading edge
#
# The biological/statistical hypothesis is that the LEFT side is relatively
# NB1-like whereas the high-loading terminal RIGHT side becomes increasingly
# NB2-like because excess variance grows disproportionately toward the
# leading edge.
#
# IMPORTANT:
#
#     The historical fixed 5,000-feature reference is NOT used.
#
# The algorithm searches the rank axis itself to determine where the strongest
# NB1-like -> NB2-like transition occurs.
#
#
# =============================================================================
# FEATURE-LEVEL NB-RELATED DIAGNOSTICS
# =============================================================================
#
# For each ranked feature:
#
#     mu       = empirical mean across samples
#     variance = empirical variance across samples
#
#     excess_variance = max(variance - mu, 0)
#
#     NB2 =
#         log(1 + excess_variance)
#
#     NB2-NB1 =
#         log(1 + excess_variance) - log(1 + mu)
#
#     alpha =
#         max((variance - mu) / mu^2, 0)
#
#     alpha*mu =
#         log(1 + alpha*mu)
#
# NB2-NB1 is the PRIMARY transition variable.
#
# NB2 and alpha*mu are used as directional corroboration.
#
# These are moment-derived/descriptive diagnostics. They are not formal
# likelihood-ratio statistics.
#
#
# =============================================================================
# CUTOFF SEARCH
# =============================================================================
#
# Every admissible integer cutoff c is examined.
#
# For every cutoff, the algorithm constructs the LARGEST POSSIBLE symmetric
# pair of adjacent windows:
#
#     w(c) = min(c - 1, N - c + 1)
#
#     LEFT(c)  = ranks [c - w(c), ..., c - 1]
#     RIGHT(c) = ranks [c, ..., c + w(c) - 1]
#
# Therefore both comparison regions always contain exactly the same number
# of features.
#
# Importantly, when c lies in the terminal half of the rank axis:
#
#     w(c) = N - c + 1
#
# and therefore:
#
#     RIGHT(c) = c ... N
#
# This is exactly the original leading-edge construction:
#
#     RIGHT = cutoff through terminal rank
#
# with an immediately preceding equal-sized LEFT region.
#
#
# =============================================================================
# WHY THIS FIXES THE PREVIOUS 15,879 ARTIFACT
# =============================================================================
#
# The previous code used:
#
#     difference / standard error
#
# Since:
#
#     SE decreases approximately as 1/sqrt(n)
#
# very large LEFT/RIGHT regions received an artificial mathematical advantage.
# That caused the solution to collapse toward the 50/50 midpoint:
#
#     cutoff ~ 15,879
#
# for N ~ 31,756.
#
# THIS VERSION DOES NOT USE STANDARD ERROR TO SELECT THE CUTOFF.
#
# Instead it uses a scale-standardized EFFECT SIZE:
#
#                       mean_RIGHT - mean_LEFT
#     score(c) = -----------------------------------------
#                sqrt((variance_LEFT + variance_RIGHT)/2)
#
# for NB2-NB1.
#
# There is NO sqrt(n) term.
#
# Therefore a cutoff does not become better merely because it contains more
# sites.
#
# The intended behavior is:
#
#     Too far LEFT:
#         RIGHT contains too many NB1-like sites -> contrast diluted.
#
#     Near transition:
#         LEFT relatively NB1-like and RIGHT relatively NB2-like
#         -> contrast maximized.
#
#     Too far RIGHT:
#         LEFT and RIGHT both increasingly lie inside the NB2-like leading edge
#         -> contrast falls again.
#
# This creates an interior maximum at the transition if such a transition
# exists.
#
#
# =============================================================================
# ANTI-DEGENERACY GUARDRAIL
# =============================================================================
#
# The only location restriction is that both local comparison windows must
# contain at least 5% of all features.
#
# For approximately 31,756 features:
#
#     5% ~ 1,588 features
#
# This exists only to prevent an extreme terminal fluctuation involving a tiny
# number of sites (e.g. 20 or 50 sites) from winning.
#
# It does NOT encode the historical 5,000-site result.
#
# Therefore the program can select, for example:
#
#     2,000 sites
#     3,700 sites
#     4,800 sites
#     5,300 sites
#     7,000 sites
#
# if that is where the actual transition is strongest.
#
#
# =============================================================================
# SELECTION CRITERIA
# =============================================================================
#
# A candidate cutoff is eligible only when:
#
#     RIGHT mean NB2-NB1 > LEFT mean NB2-NB1
#     RIGHT mean NB2     > LEFT mean NB2
#     RIGHT mean alpha*mu > LEFT mean alpha*mu
#
# Among eligible cutoffs, the one with the maximum NB2-NB1 standardized
# effect-size score is selected.
#
# No smoothing of the cutoff-response curve is used for selection.
#
# This avoids introducing an additional smoothing bandwidth into the
# optimization problem. Moving the cutoff by one rank already changes the
# windows by only one or a few features, making the response curve naturally
# highly correlated from rank to rank.
#
#
# =============================================================================
# BOUNDARY CHECK
# =============================================================================
#
# If the optimum occurs exactly at either admissible search boundary, it is
# flagged as a BOUNDARY_OPTIMUM rather than interpreted as a validated
# transition.
#
# This protects against a monotonic score curve in which no internal optimum
# was actually found.
#
#
# =============================================================================
# VARIANCE GEOMETRY
# =============================================================================
#
# The original variance-curve analysis is retained independently:
#
#     empirical variance
#         ->
#     log1p variance
#         ->
#     smoothing spline
#         ->
#     second derivative
#         ->
#     curvature zero crossings
#
# The nearest curvature zero crossing to the optimized cutoff is reported as
# geometric corroboration.
#
# It does NOT determine the optimized cutoff.
#
#
# =============================================================================
# BOOTSTRAP
# =============================================================================
#
# Each bootstrap replicate resamples biological samples WITHIN THE ARM and
# repeats the COMPLETE analysis:
#
#     sample resampling
#         ->
#     PC1
#         ->
#     rank all features
#         ->
#     mean/variance
#         ->
#     NB diagnostics
#         ->
#     scan every admissible cutoff
#         ->
#     select optimal transition
#
# The bootstrap therefore estimates stability of the complete selection
# procedure rather than merely resampling a previously selected region.
#
# Stability requirements:
#
#     >= 90% valid interior bootstrap optima
#
#     >= 80% of valid bootstrap optima within 3% of the total rank axis
#             of the full-data optimum
#
#     bootstrap cutoff IQR <= 3% of total N
#
# Jaccard feature overlap is reported descriptively but is NOT used as a
# pass/fail criterion because the inferential target is the rank-space
# transition, not exact feature identity under PC1 re-estimation.
#
#
# =============================================================================
# INTERPRETATION
# =============================================================================
#
# If the selected cutoff is c:
#
#     terminal leading-edge size = N - c + 1
#
# Thus, with N = 31,756:
#
#     c = 26,757  ->  5,000 terminal sites
#
# But 26,757 and 5,000 are NEVER supplied to the algorithm.
#
# If a ~5,000-site leading edge is intrinsic to the data, it should emerge
# from the full-axis search and bootstrap distribution.
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <-
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"


OUT_ROOT <-
  "/root/REAPER98632/exports/manuscript_nb_transition_repaired"


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
# Anti-degeneracy guardrail.
# -----------------------------------------------------------------------------

MIN_WINDOW_FRACTION <- 0.05

MIN_WINDOW_ABSOLUTE <- 100L


# -----------------------------------------------------------------------------
# Variance geometry.
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
# Supplementary DESeq2 analysis.
# -----------------------------------------------------------------------------

RUN_DESEQ2_SUPPLEMENT <- TRUE

DESEQ2_RANK_METHOD <- "normalized_log1p"

# Allowed:
#
#     "normalized_log1p"
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
  # Preserve compatibility with the previously encountered blank feature ID.
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
      which(blank_ids)
    )


    message(
      "Assigned deterministic row IDs to ",
      sum(blank_ids),
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
# BASIC NUMERICAL HELPERS
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


  if (
    length(u) == 0L
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
    length(u)
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
# FEATURE-LEVEL NB DIAGNOSTICS
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
        alpha_hat *
          mu
      ),

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# FAST PREFIX-MOMENT CALCULATIONS
# =============================================================================

window_moments <- function(
    x,
    starts,
    ends,
    n_vec) {


  prefix_sum <- c(
    0,
    cumsum(x)
  )


  prefix_sum2 <- c(
    0,
    cumsum(
      x * x
    )
  )


  sums <-
    prefix_sum[
      ends + 1L
    ] -
      prefix_sum[
        starts
      ]


  sums2 <-
    prefix_sum2[
      ends + 1L
    ] -
      prefix_sum2[
        starts
      ]


  means <- sums /
    n_vec


  vars <- rep(
    NA_real_,
    length(means)
  )


  valid_n <- n_vec > 1L


  vars[
    valid_n
  ] <- (

    sums2[
      valid_n
    ] -
      (
        sums[
          valid_n
        ]^2 /
          n_vec[
            valid_n
          ]
      )

  ) /
    (
      n_vec[
        valid_n
      ] -
        1L
    )


  vars[
    is.finite(vars) &
    vars < 0 &
    abs(vars) < 1e-10
  ] <- 0


  vars[
    !is.finite(vars)
  ] <- NA_real_


  vars <- pmax(
    vars,
    0
  )


  list(
    mean = means,
    var = vars
  )
}


# =============================================================================
# FULL-RANK-AXIS TRANSITION SCAN
# =============================================================================

scan_all_cutoffs <- function(
    feature_df,
    min_window_fraction = MIN_WINDOW_FRACTION,
    min_window_absolute = MIN_WINDOW_ABSOLUTE) {


  N <- nrow(
    feature_df
  )


  min_window_n <- max(

    as.integer(
      min_window_absolute
    ),

    as.integer(
      ceiling(
        min_window_fraction *
          N
      )
    ),

    2L
  )


  # ---------------------------------------------------------------------------
  # Scan every rank where a symmetric LEFT/RIGHT comparison of at least the
  # minimum allowable size exists.
  #
  # This searches essentially the entire rank axis except the protected
  # terminal margins.
  # ---------------------------------------------------------------------------

  min_cutoff <- min_window_n +
    1L


  max_cutoff <- N -
    min_window_n +
    1L


  if (
    min_cutoff >
      max_cutoff
  ) {

    return(
      list(
        valid = FALSE,
        reason = "no_admissible_cutoffs",
        scan_df = data.frame(),
        selected = NULL
      )
    )
  }


  cutoffs <- seq.int(
    min_cutoff,
    max_cutoff
  )


  # ---------------------------------------------------------------------------
  # Largest symmetric window available around every cutoff.
  #
  # For cutoffs in the RIGHT half:
  #
  #     window_n = N - cutoff + 1
  #
  # so RIGHT automatically extends all the way to N.
  # ---------------------------------------------------------------------------

  window_n <- pmin(
    cutoffs - 1L,
    N - cutoffs + 1L
  )


  left_start <- cutoffs -
    window_n


  left_end <- cutoffs -
    1L


  right_start <- cutoffs


  right_end <- cutoffs +
    window_n -
    1L


  # ---------------------------------------------------------------------------
  # Primary and secondary diagnostics.
  # ---------------------------------------------------------------------------

  primary <- feature_df$NB2_NB1
  nb2 <- feature_df$NB2
  alpha_mu <- feature_df$alpha_mu


  if (
    any(!is.finite(primary)) ||
    any(!is.finite(nb2)) ||
    any(!is.finite(alpha_mu))
  ) {

    return(
      list(
        valid = FALSE,
        reason = "nonfinite_nb_metric",
        scan_df = data.frame(),
        selected = NULL
      )
    )
  }


  L_primary <- window_moments(
    primary,
    left_start,
    left_end,
    window_n
  )


  R_primary <- window_moments(
    primary,
    right_start,
    right_end,
    window_n
  )


  L_nb2 <- window_moments(
    nb2,
    left_start,
    left_end,
    window_n
  )


  R_nb2 <- window_moments(
    nb2,
    right_start,
    right_end,
    window_n
  )


  L_alpha <- window_moments(
    alpha_mu,
    left_start,
    left_end,
    window_n
  )


  R_alpha <- window_moments(
    alpha_mu,
    right_start,
    right_end,
    window_n
  )


  # ---------------------------------------------------------------------------
  # Directional differences.
  # ---------------------------------------------------------------------------

  delta_primary <-
    R_primary$mean -
      L_primary$mean


  delta_nb2 <-
    R_nb2$mean -
      L_nb2$mean


  delta_alpha <-
    R_alpha$mean -
      L_alpha$mean


  # ---------------------------------------------------------------------------
  # SCALE-STANDARDIZED EFFECT SIZE.
  #
  # CRITICAL:
  #
  # There is deliberately NO division by standard error and NO sqrt(n).
  #
  # Therefore large windows do not automatically receive larger scores.
  # ---------------------------------------------------------------------------

  pooled_scale <- sqrt(
    (
      L_primary$var +
        R_primary$var
    ) /
      2
  )


  effect_size <- rep(
    NA_real_,
    length(cutoffs)
  )


  scale_ok <- (
    is.finite(pooled_scale) &
      pooled_scale >
      .Machine$double.eps
  )


  effect_size[
    scale_ok
  ] <- delta_primary[
    scale_ok
  ] /
    pooled_scale[
      scale_ok
    ]


  # ---------------------------------------------------------------------------
  # Candidate table.
  # ---------------------------------------------------------------------------

  scan_df <- data.frame(

    cutoff =
      cutoffs,

    window_n =
      window_n,

    window_fraction =
      window_n /
        N,

    final_terminal_tail_n =
      N -
        cutoffs +
        1L,

    final_terminal_tail_fraction =
      (
        N -
          cutoffs +
          1L
      ) /
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
      L_primary$mean,

    right_mean_NB2_NB1 =
      R_primary$mean,

    delta_NB2_NB1 =
      delta_primary,


    left_var_NB2_NB1 =
      L_primary$var,

    right_var_NB2_NB1 =
      R_primary$var,


    left_mean_NB2 =
      L_nb2$mean,

    right_mean_NB2 =
      R_nb2$mean,

    delta_NB2 =
      delta_nb2,


    left_mean_alpha_mu =
      L_alpha$mean,

    right_mean_alpha_mu =
      R_alpha$mean,

    delta_alpha_mu =
      delta_alpha,


    transition_effect_size =
      effect_size,

    stringsAsFactors = FALSE
  )


  # ---------------------------------------------------------------------------
  # Direction must be consistent with the NB1 -> NB2 hypothesis.
  # ---------------------------------------------------------------------------

  eligible <- (

    is.finite(
      scan_df$transition_effect_size
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
    !any(eligible)
  ) {

    return(
      list(
        valid = FALSE,
        reason = "no_directionally_consistent_transition",
        scan_df = scan_df,
        selected = NULL,
        min_window_n = min_window_n
      )
    )
  }


  eligible_idx <- which(
    eligible
  )


  best_idx <- eligible_idx[
    which.max(
      scan_df$transition_effect_size[
        eligible_idx
      ]
    )
  ]


  scan_df$selected[
    best_idx
  ] <- TRUE


  selected <- scan_df[
    best_idx,
    ,
    drop = FALSE
  ]


  # ---------------------------------------------------------------------------
  # Boundary diagnostic.
  #
  # If maximum occurs exactly at the minimum or maximum admissible cutoff,
  # the data did not demonstrate an internal optimum.
  # ---------------------------------------------------------------------------

  boundary_hit <- (
    selected$cutoff[
      1L
    ] ==
      min_cutoff ||

      selected$cutoff[
        1L
      ] ==
      max_cutoff
  )


  selected$boundary_hit <-
    boundary_hit


  list(

    valid = TRUE,

    reason = NA_character_,

    scan_df = scan_df,

    selected = selected,

    min_window_n =
      min_window_n,

    min_cutoff =
      min_cutoff,

    max_cutoff =
      max_cutoff,

    boundary_hit =
      boundary_hit
  )
}


# =============================================================================
# VARIANCE SPLINE / CURVATURE GEOMETRY
# =============================================================================

find_d2_zero_crossings <- function(
    dense_x,
    dense_d2) {


  ok <- (
    is.finite(dense_x) &
      is.finite(dense_d2)
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
        crossing_rank = numeric(0),
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


  spline_fit <- stats::smooth.spline(
    x = ranks,
    y = log_var,
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
        rank = NA_real_,
        distance = NA_real_
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
# COMPLETE FULL-DATA TRACK ANALYSIS
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


  geometry <- compute_variance_geometry(
    feature_df
  )


  if (
    !isTRUE(scan$valid)
  ) {

    return(
      list(

        valid = FALSE,

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
          geometry,

        nearest_crossing =
          list(
            rank = NA_real_,
            distance = NA_real_
          )
      )
    )
  }


  selected_cutoff <- as.integer(
    scan$selected$cutoff[
      1L
    ]
  )


  nearest_crossing <- nearest_d2_crossing(
    geometry$zero_df,
    selected_cutoff
  )


  list(

    valid = TRUE,

    reason = NA_character_,

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
# SELECTED REGION SUMMARY
# =============================================================================

summarize_selected_transition <- function(
    feature_df,
    selected_row) {


  cutoff <- as.integer(
    selected_row$cutoff[
      1L
    ]
  )


  left_start <- as.integer(
    selected_row$left_start[
      1L
    ]
  )


  left_end <- as.integer(
    selected_row$left_end[
      1L
    ]
  )


  right_start <- as.integer(
    selected_row$right_start[
      1L
    ]
  )


  right_end <- as.integer(
    selected_row$right_end[
      1L
    ]
  )


  N <- nrow(
    feature_df
  )


  left_df <- feature_df[
    left_start:left_end,
    ,
    drop = FALSE
  ]


  right_window_df <- feature_df[
    right_start:right_end,
    ,
    drop = FALSE
  ]


  terminal_df <- feature_df[
    cutoff:N,
    ,
    drop = FALSE
  ]


  data.frame(

    cutoff =
      cutoff,

    comparison_window_n =
      nrow(left_df),

    final_terminal_tail_n =
      nrow(terminal_df),

    final_terminal_tail_fraction =
      nrow(terminal_df) /
        N,


    left_start =
      left_start,

    left_end =
      left_end,

    right_window_start =
      right_start,

    right_window_end =
      right_end,


    left_median_NB2_NB1 =
      median(
        left_df$NB2_NB1,
        na.rm = TRUE
      ),

    right_median_NB2_NB1 =
      median(
        right_window_df$NB2_NB1,
        na.rm = TRUE
      ),

    delta_median_NB2_NB1 =
      median(
        right_window_df$NB2_NB1,
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
        right_window_df$NB2,
        na.rm = TRUE
      ),

    delta_median_NB2 =
      median(
        right_window_df$NB2,
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
        right_window_df$alpha_mu,
        na.rm = TRUE
      ),

    delta_median_alpha_mu =
      median(
        right_window_df$alpha_mu,
        na.rm = TRUE
      ) -
        median(
          left_df$alpha_mu,
          na.rm = TRUE
        ),


    terminal_median_NB2_NB1 =
      median(
        terminal_df$NB2_NB1,
        na.rm = TRUE
      ),

    terminal_median_NB2 =
      median(
        terminal_df$NB2,
        na.rm = TRUE
      ),

    terminal_median_alpha_mu =
      median(
        terminal_df$alpha_mu,
        na.rm = TRUE
      ),

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# BOOTSTRAP COMPLETE CUTOFF SEARCH
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
    !isTRUE(full_analysis$valid)
  ) {

    stop(
      "Cannot bootstrap an invalid full-data track."
    )
  }


  full_selected <- full_analysis$scan$selected


  full_cutoff <- as.integer(
    full_selected$cutoff[
      1L
    ]
  )


  full_terminal_idx <- full_analysis$rank_order[
    full_cutoff:N
  ]


  full_terminal_features <- rownames(
    rank_matrix
  )[
    full_terminal_idx
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

        comparison_window_n =
          NA_integer_,

        terminal_tail_n =
          NA_integer_,

        terminal_tail_fraction =
          NA_real_,

        transition_effect_size =
          NA_real_,

        delta_NB2_NB1 =
          NA_real_,

        delta_NB2 =
          NA_real_,

        delta_alpha_mu =
          NA_real_,

        boundary_hit =
          NA,

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


    result <- tryCatch(

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


        scan <- scan_all_cutoffs(
          feature_df
        )


        if (
          !isTRUE(scan$valid)
        ) {

          invalid_row(
            scan$reason
          )

        } else if (
          isTRUE(scan$boundary_hit)
        ) {

          invalid_row(
            "boundary_optimum"
          )

        } else {

          selected <- scan$selected


          cutoff <- as.integer(
            selected$cutoff[
              1L
            ]
          )


          terminal_idx <- rank_order[
            cutoff:N
          ]


          terminal_features <- rownames(
            rank_matrix
          )[
            terminal_idx
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

            comparison_window_n =
              as.integer(
                selected$window_n[
                  1L
                ]
              ),

            terminal_tail_n =
              N -
                cutoff +
                1L,

            terminal_tail_fraction =
              (
                N -
                  cutoff +
                  1L
              ) /
                N,

            transition_effect_size =
              selected$transition_effect_size[
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

            boundary_hit =
              FALSE,

            jaccard =
              jaccard_similarity(
                full_terminal_features,
                terminal_features
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


    result
  }


  bootstrap_ids <- seq_len(
    bootstrap_n
  )


  if (
    .Platform$OS.type == "unix" &&
    cores > 1L
  ) {

    bootstrap_list <- parallel::mclapply(

      bootstrap_ids,

      worker,

      mc.cores = cores,

      mc.preschedule = TRUE,

      mc.set.seed = FALSE
    )

  } else {

    bootstrap_list <- lapply(
      bootstrap_ids,
      worker
    )
  }


  bootstrap_df <- bind_rows(
    bootstrap_list
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

    recovery_rate <- mean(
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


    cutoff_quantiles <- as.numeric(
      stats::quantile(
        valid_df$cutoff,
        probs = c(
          0.025,
          0.50,
          0.975
        ),
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
    )


    tail_quantiles <- as.numeric(
      stats::quantile(
        valid_df$terminal_tail_n,
        probs = c(
          0.025,
          0.50,
          0.975
        ),
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
    )


    median_effect <- safe_median(
      valid_df$transition_effect_size
    )


    median_jaccard <- safe_median(
      valid_df$jaccard
    )

  } else {

    recovery_rate <- NA_real_

    cutoff_iqr <- NA_real_

    cutoff_iqr_fraction <- NA_real_

    cutoff_quantiles <- c(
      NA_real_,
      NA_real_,
      NA_real_
    )

    tail_quantiles <- c(
      NA_real_,
      NA_real_,
      NA_real_
    )

    median_effect <- NA_real_

    median_jaccard <- NA_real_
  }


  pass_full_interior <-
    !isTRUE(
      full_analysis$scan$boundary_hit
    )


  pass_valid_rate <- (
    is.finite(valid_rate) &&
      valid_rate >=
      MIN_BOOTSTRAP_VALID_RATE
  )


  pass_recovery <- (
    is.finite(recovery_rate) &&
      recovery_rate >=
      MIN_CUTOFF_RECOVERY_RATE
  )


  pass_iqr <- (
    is.finite(cutoff_iqr_fraction) &&
      cutoff_iqr_fraction <=
      MAX_CUTOFF_IQR_FRACTION
  )


  pass_stability <- (
    pass_full_interior &&
      pass_valid_rate &&
      pass_recovery &&
      pass_iqr
  )


  if (
    isTRUE(
      full_analysis$scan$boundary_hit
    )
  ) {

    status <- "BOUNDARY_OPTIMUM"

  } else if (
    pass_stability
  ) {

    status <- "PASS"

  } else {

    status <- "UNSTABLE"
  }


  stability_summary <- data.frame(

    Status =
      status,


    FullCutoff =
      full_cutoff,

    FullComparisonWindowN =
      full_selected$window_n[
        1L
      ],

    FullTerminalTailN =
      N -
        full_cutoff +
        1L,

    FullTerminalTailFraction =
      (
        N -
          full_cutoff +
          1L
      ) /
        N,

    FullTransitionEffectSize =
      full_selected$transition_effect_size[
        1L
      ],

    FullBoundaryHit =
      full_analysis$scan$boundary_hit,


    BootstrapN =
      bootstrap_n,

    ValidBootstraps =
      valid_n,

    BootstrapValidRate =
      valid_rate,


    CutoffToleranceN =
      tolerance_n,

    CutoffRecoveryRate =
      recovery_rate,

    CutoffIQR =
      cutoff_iqr,

    CutoffIQRFraction =
      cutoff_iqr_fraction,


    BootstrapCutoffCI025 =
      cutoff_quantiles[
        1L
      ],

    BootstrapCutoffMedian =
      cutoff_quantiles[
        2L
      ],

    BootstrapCutoffCI975 =
      cutoff_quantiles[
        3L
      ],


    BootstrapTailNCI025 =
      tail_quantiles[
        1L
      ],

    BootstrapTailNMedian =
      tail_quantiles[
        2L
      ],

    BootstrapTailNCI975 =
      tail_quantiles[
        3L
      ],


    MedianBootstrapEffectSize =
      median_effect,

    MedianJaccard =
      median_jaccard,


    PassFullInteriorOptimum =
      pass_full_interior,

    PassValidRate =
      pass_valid_rate,

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
      layout =
        grid.layout(
          nrow = 3L,
          ncol = 1L,
          heights =
            unit(
              c(
                1.0,
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
          layout.pos.row = i,
          layout.pos.col = 1L
        )
    )
  }


  dev.off()
}


# =============================================================================
# FIGURE
# =============================================================================

build_figure <- function(
    comparison_name,
    arm_name,
    track,
    full_analysis,
    selected_summary,
    stability_summary,
    out_file) {


  N <- nrow(
    full_analysis$feature_df
  )


  selected <- full_analysis$scan$selected


  cutoff <- as.integer(
    selected$cutoff[
      1L
    ]
  )


  left_start <- as.integer(
    selected$left_start[
      1L
    ]
  )


  left_end <- as.integer(
    selected$left_end[
      1L
    ]
  )


  right_start <- as.integer(
    selected$right_start[
      1L
    ]
  )


  right_end <- as.integer(
    selected$right_end[
      1L
    ]
  )


  # ---------------------------------------------------------------------------
  # PANEL 1: original variance geometry.
  # ---------------------------------------------------------------------------

  variance_df <-
    full_analysis$geometry$variance_df


  zero_df <-
    full_analysis$geometry$zero_df


  p1 <- ggplot(
    variance_df,
    aes(
      rank,
      smooth_log1p_empirical_variance
    )
  ) +

    annotate(
      "rect",
      xmin = left_start,
      xmax = left_end,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey90"
    ) +

    annotate(
      "rect",
      xmin = right_start,
      xmax = right_end,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey75",
      alpha = 0.65
    ) +

    geom_line(
      linewidth = 0.9
    ) +

    geom_vline(
      xintercept = cutoff,
      linewidth = 0.9
    )


  if (
    nrow(zero_df) > 0L
  ) {

    p1 <- p1 +

      geom_vline(
        data = zero_df,
        aes(
          xintercept = crossing_rank
        ),
        inherit.aes = FALSE,
        linetype = "dotted",
        linewidth = 0.25,
        alpha = 0.45
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
          "; terminal leading edge = ",
          N - cutoff + 1L,
          " features; nearest d2 crossing = ",
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
        "Rank (ascending absolute PC1 loading)",

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
  # PANEL 2: complete transition-response curve.
  # ---------------------------------------------------------------------------

  scan_df <-
    full_analysis$scan$scan_df


  selected_df <- scan_df[
    scan_df$selected,
    ,
    drop = FALSE
  ]


  ci_low <-
    stability_summary$BootstrapCutoffCI025[
      1L
    ]


  ci_high <-
    stability_summary$BootstrapCutoffCI975[
      1L
    ]


  p2 <- ggplot(
    scan_df,
    aes(
      cutoff,
      transition_effect_size
    )
  )


  if (
    is.finite(ci_low) &&
    is.finite(ci_high)
  ) {

    p2 <- p2 +

      annotate(
        "rect",
        xmin = ci_low,
        xmax = ci_high,
        ymin = -Inf,
        ymax = Inf,
        alpha = 0.10
      )
  }


  p2 <- p2 +

    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +

    geom_line(
      linewidth = 0.8
    ) +

    geom_vline(
      xintercept = cutoff,
      linewidth = 0.9
    ) +

    geom_point(
      data = selected_df,
      aes(
        cutoff,
        transition_effect_size
      ),
      size = 3
    ) +

    labs(

      title =
        "Full-rank-axis NB1-like to NB2-like transition scan",

      subtitle =
        paste0(
          "Primary score = standardized NB2-NB1 effect size; ",
          "no standard-error or sqrt(n) term; status = ",
          stability_summary$Status[
            1L
          ]
        ),

      x =
        "Candidate cutoff rank",

      y =
        "NB2-NB1 transition effect size"
    ) +

    theme_bw(
      base_size = 11
    ) +

    theme(
      panel.grid.minor =
        element_blank()
    )


  # ---------------------------------------------------------------------------
  # PANEL 3: selected LEFT/RIGHT medians.
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
      ),

    stringsAsFactors = FALSE
  )


  p3 <- ggplot(
    region_df,
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
      linewidth = 0.8
    ) +

    geom_point(
      aes(
        x = LEFT,
        shape = "LEFT"
      ),
      size = 3.2
    ) +

    geom_point(
      aes(
        x = RIGHT,
        shape = "RIGHT"
      ),
      size = 3.2
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
        "Selected matched LEFT versus RIGHT transition",

      subtitle =
        paste0(
          "Comparison-window n = ",
          selected_summary$comparison_window_n[
            1L
          ],
          "; terminal leading-edge n = ",
          selected_summary$final_terminal_tail_n[
            1L
          ],
          "; bootstrap cutoff median = ",
          round(
            stability_summary$BootstrapCutoffMedian[
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
    ": scanning full admissible rank axis..."
  )


  full_analysis <- analyze_full_track(
    rank_matrix,
    metric_matrix
  )


  # ---------------------------------------------------------------------------
  # Always save ranked features and geometry when available.
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
  # No directionally valid full-data transition.
  # ---------------------------------------------------------------------------

  if (
    !isTRUE(full_analysis$valid)
  ) {

    fail_row <- data.frame(

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
        "NO_VALID_TRANSITION",

      Reason =
        full_analysis$reason,

      stringsAsFactors = FALSE
    )


    write.csv(
      fail_row,
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
      " | NO_VALID_TRANSITION"
    )


    return(
      list(
        summary = fail_row,
        stability_summary = data.frame()
      )
    )
  }


  selected <- full_analysis$scan$selected


  cutoff <- as.integer(
    selected$cutoff[
      1L
    ]
  )


  selected_summary <- summarize_selected_transition(
    full_analysis$feature_df,
    selected
  )


  message(
    "[",
    track,
    "] ",
    comparison_name,
    " ",
    arm_name,
    ": full-data optimum cutoff=",
    cutoff,
    ", terminal leading edge n=",
    selected_summary$final_terminal_tail_n[
      1L
    ],
    ", effect size=",
    sprintf(
      "%.4f",
      selected$transition_effect_size[
        1L
      ]
    ),
    "; bootstrapping ",
    BOOTSTRAP_N,
    " replicates..."
  )


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
      cutoff,

    ComparisonWindowN =
      selected$window_n[
        1L
      ],

    TerminalLeadingEdgeN =
      selected_summary$final_terminal_tail_n[
        1L
      ],

    TerminalLeadingEdgeFraction =
      selected_summary$final_terminal_tail_fraction[
        1L
      ],


    TransitionEffectSize =
      selected$transition_effect_size[
        1L
      ],


    MeanDelta_NB2_NB1 =
      selected$delta_NB2_NB1[
        1L
      ],

    MeanDelta_NB2 =
      selected$delta_NB2[
        1L
      ],

    MeanDelta_alpha_mu =
      selected$delta_alpha_mu[
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


    BoundaryHit =
      full_analysis$scan$boundary_hit,


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
    cutoff,
    " | terminal leading edge n=",
    selected_summary$final_terminal_tail_n[
      1L
    ],
    " | effect size=",
    sprintf(
      "%.4f",
      selected$transition_effect_size[
        1L
      ]
    ),
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
    " | 95% terminal-tail interval=",
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
# RUN COMPLETE ANALYSIS
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
  "N_eff reference: REMOVED"
)


message(
  "Standard-error / t-like cutoff score: REMOVED"
)


message(
  "Primary cutoff criterion: NB2-NB1 standardized effect size"
)


message(
  "Minimum comparison window: max(",
  MIN_WINDOW_ABSOLUTE,
  " features, ",
  100 * MIN_WINDOW_FRACTION,
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
    # MAIN TRACK
    #
    # PC1 ranking:
    #     CPM log1p
    #
    # Mean/variance/NB diagnostics:
    #     raw counts
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
      length(all_summary_rows) +
        1L
    ]] <- main_res$summary


    if (
      nrow(
        main_res$stability_summary
      ) > 0L
    ) {

      all_stability_rows[[
        length(all_stability_rows) +
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
        is.null(deseq2_obj)
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
          length(all_summary_rows) +
            1L
        ]] <- deseq2_res$summary


        if (
          nrow(
            deseq2_res$stability_summary
          ) > 0L
        ) {

          all_stability_rows[[
            length(all_stability_rows) +
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
  "All admissible cutoff ranks were searched."
)


message(
  "The historical 5,000-site cutoff was not supplied anywhere."
)


message(
  "The previous sample-size-biased standard-error score was removed."
)


message(
  "The selected cutoff maximizes NB2-NB1 effect-size separation."
)


message(
  "NB2 and alpha*mu must independently move in the expected RIGHT direction."
)


message(
  "Boundary optima are explicitly flagged rather than accepted."
)


message(
  "Bootstrap repeats the entire PC1/ranking/cutoff-selection procedure."
)


message(
  "Outputs written to: ",
  OUT_ROOT
)


message(
  "============================================================"
)
