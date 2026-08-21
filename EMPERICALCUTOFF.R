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
# DESIGN-AWARE POOLED VARIANCE + DATA-DERIVED TERMINAL TRANSITION
# =============================================================================
#
# CORE QUESTION
# -------------
# Where, along the COMPLETE PC1-ranked feature axis, does a persistent
# high-variance terminal regime begin?
#
# Features are ranked independently within each experimental arm from:
#
#     lowest -> highest absolute PC1 loading
#
# Therefore the biologically relevant leading edge lies on the RIGHT.
#
#
# =============================================================================
# KEY IMPROVEMENT OVER THE ORIGINAL METHOD
# =============================================================================
#
# ORIGINAL:
#
#     arm-specific PC1
#         ->
#     arm-specific empirical variance estimated from ~5 samples
#         ->
#     fixed 5,000-feature reference
#         ->
#     nearest d2 crossings
#
#
# THIS VERSION:
#
#     all 40 samples
#         ->
#     DESeq2 normalization with design ~ group
#         ->
#     remove the mean of each of the 8 experimental groups
#         ->
#     pool WITHIN-GROUP residual sum of squares
#         ->
#     pooled feature variance with approximately 32 residual df
#         ->
#     arm-specific PC1 ranking
#         ->
#     multiscale terminal-transition detector
#         ->
#     data-derived Ref
#         ->
#     original d2 Anchor / Terminal geometry
#
#
# Thus treatment/control/time-group mean differences are NOT counted as
# variance, but all groups contribute information to estimating underlying
# within-group variance.
#
#
# =============================================================================
# PRIMARY GEOMETRY VARIANCE
# =============================================================================
#
# Let y_ij denote the DESeq2-normalized count for feature i, sample j.
#
# For each experimental group g:
#
#     group mean = mean(y_ij | j belongs to g)
#
# Residuals are:
#
#     e_ij = y_ij - group_mean_ig
#
# The pooled within-group variance is:
#
#                       sum_g sum_j e_ij^2
#     V_i(pooled) = --------------------------------
#                     sum_g (n_g - 1)
#
# With 40 samples and 8 groups of approximately 5 samples:
#
#     residual df ~= 40 - 8 = 32
#
# rather than approximately 4 df from estimating variance separately inside
# each arm.
#
# This pooled variance is COMMON across arms, but each arm orders that same
# variance vector according to its OWN absolute-PC1-loading ranking.
#
# Therefore each arm retains its own leading-edge geometry.
#
#
# =============================================================================
# WHY PRIMARY GEOMETRY DOES NOT USE NB2 MODEL-IMPLIED VARIANCE
# =============================================================================
#
# DESeq2 also estimates:
#
#     Var(K_ij) = mu_ij + alpha_i * mu_ij^2
#
# using shrunken gene-wise dispersion alpha_i.
#
# Those model-implied variances are calculated below as a sensitivity analysis.
#
# They are NOT used for primary cutoff selection because the downstream
# hypothesis is precisely that RIGHT is more NB2-like than LEFT.
#
# Using an NB2 variance model to select RIGHT and then using NB2 behavior to
# "confirm" RIGHT would make the argument partly circular.
#
# PRIMARY selection therefore uses pooled empirical within-group variance.
#
#
# =============================================================================
# TERMINAL TRANSITION DETECTOR
# =============================================================================
#
# For each arm:
#
# 1. Rank all features by ascending |PC1 loading|.
#
# 2. Reorder pooled within-group variance according to that ranking.
#
# 3. Transform:
#
#        y(r) = log(1 + pooled variance)
#
# 4. Fit smoothing splines at several prespecified smoothing levels.
#
# 5. At each smoothing level, search the mathematically admissible terminal
#    portion of the rank axis for a CONTINUOUS SEGMENTED-REGRESSION breakpoint.
#
#    Null model:
#
#        y(x) = beta0 + beta1*x
#
#    Change-point model:
#
#        y(x) = beta0 + beta1*x + gamma*max(0, x-c)
#
#    where c is the candidate terminal-transition location.
#
# 6. Every candidate breakpoint is evaluated over the SAME fixed terminal-half
#    domain. Therefore candidate locations are directly comparable and cannot
#    win simply because one model contains more observations.
#
# 7. Breakpoint evidence is quantified by:
#
#        BIC_gain =
#            M * log(SSE_linear / SSE_hinge) - log(M)
#
#    where M is constant across all candidate breakpoints.
#
# 8. An eligible terminal transition requires:
#
#        gamma > 0
#
#        terminal fitted slope > 0
#
#        terminal region mean variance >
#        immediately preceding equal-sized region mean variance
#
#        BIC_gain > 0
#
# 9. The detector is repeated at multiple smoothing scales.
#
# 10. A weighted-median consensus breakpoint across scales becomes the
#     DATA-DERIVED REFERENCE, Ref.
#
#
# =============================================================================
# ORIGINAL GEOMETRY IS PRESERVED AFTER Ref IS FOUND
# =============================================================================
#
# Using the original spline at spar = 0.60:
#
#     Anchor =
#         nearest second-derivative zero crossing LEFT of Ref
#
#     Terminal =
#         nearest second-derivative zero crossing RIGHT of Ref
#
# and:
#
#     Anchor < Ref < Terminal
#
#
# Critically:
#
#     RIGHT = Anchor -> final rank N
#
#     LEFT =
#         equal-sized immediately preceding block
#
# Terminal does NOT truncate RIGHT.
#
#
# =============================================================================
# BOOTSTRAP
# =============================================================================
#
# Each bootstrap replicate:
#
#     resamples samples WITHIN EACH of all 8 experimental groups
#         ->
#     recalculates pooled within-group variance using all groups
#         ->
#     resamples the focal arm using the SAME group-specific bootstrap draw
#         ->
#     recalculates PC1
#         ->
#     reranks all features
#         ->
#     reruns terminal-transition detection
#         ->
#     reruns Anchor / Ref / Terminal geometry
#
# Thus variance estimation and arm-specific ranking uncertainty are both
# propagated.
#
#
# Bootstrap stability requires:
#
#     valid geometry >= 90%
#
#     Ref recovery within 3% of N >= 80%
#
#     Anchor recovery within 3% of N >= 80%
#
#     Terminal recovery within 3% of N >= 80%
#
#     Ref IQR / N <= 3%
#
#     Anchor IQR / N <= 3%
#
#     Terminal IQR / N <= 3%
#
#     RIGHT-fraction IQR <= 3 percentage points
#
# Jaccard overlap is descriptive only.
#
#
# =============================================================================
# DOWNSTREAM NB2-RELATED CORROBORATION
# =============================================================================
#
# These quantities are calculated AFTER primary geometric selection:
#
#     NB2 =
#         log(1 + max(variance - mu, 0))
#
#     NB2-NB1 =
#         log(1 + max(variance - mu, 0)) - log(1 + mu)
#
#     alpha =
#         max((variance - mu) / mu^2, 0)
#
#     alpha*mu =
#         log(1 + alpha*mu)
#
# Main track:
#
#     PC1 ranking = CPM log1p
#     NB metrics  = raw counts
#
# Supplement:
#
#     PC1 ranking = DESeq2 normalized log1p or VST
#     NB metrics  = DESeq2 normalized counts
#
#
# =============================================================================
# DESEQ2 MODEL-IMPLIED VARIANCE SENSITIVITY ANALYSIS
# =============================================================================
#
# A single DESeq2 model is also fit to ALL 40 samples:
#
#     design ~ group
#
# The MAP/shrunken dispersion estimate alpha_i is combined with each group's
# fitted normalized mean:
#
#     V_model,ig =
#         mu_ig + alpha_i * mu_ig^2
#
# The same terminal-transition detector is applied to this model-implied
# variance as a SENSITIVITY ANALYSIS ONLY.
#
# Agreement between primary pooled-residual geometry and DESeq2 model-implied
# geometry strengthens the result.
#
# No historical 5,000-site reference is used anywhere.
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <-
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <-
  "/root/REAPER98632/exports/manuscript_pooled_terminal_transition"


# -----------------------------------------------------------------------------
# Explicit experimental groups
# -----------------------------------------------------------------------------

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

  RT0_ZT6 =
    list(
      control = "RT0",
      treatment = "ZT6"
    ),

  RT2_ZT8 =
    list(
      control = "RT2",
      treatment = "ZT8"
    ),

  RT4_ZT10 =
    list(
      control = "RT4",
      treatment = "ZT10"
    ),

  RT8_ZT14 =
    list(
      control = "RT8",
      treatment = "ZT14"
    )
)


# -----------------------------------------------------------------------------
# Original variance spline
# -----------------------------------------------------------------------------

VAR_SPLINE_SPAR <- 0.60

DENSE_GRID_MULTIPLIER <- 4L

DENSE_GRID_MIN <- 5000L


# -----------------------------------------------------------------------------
# Multiscale terminal-transition detector
# -----------------------------------------------------------------------------

TRANSITION_SPARS <- c(
  0.50,
  0.60,
  0.70
)

MIN_VALID_TRANSITION_SCALES <- 2L

MAX_MULTISCALE_REF_IQR_FRACTION <- 0.05


# -----------------------------------------------------------------------------
# Anti-degeneracy guardrail
#
# This does NOT specify the leading-edge size.
# -----------------------------------------------------------------------------

MIN_TRANSITION_SEGMENT_FRACTION <- 0.05

MIN_TRANSITION_SEGMENT_ABSOLUTE <- 100L


# -----------------------------------------------------------------------------
# Bootstrap
# -----------------------------------------------------------------------------

BOOTSTRAP_N <- 500L

BOOTSTRAP_SEED_BASE <- 20260820L

MIN_BOOTSTRAP_VALID_RATE <- 0.90

POSITION_TOLERANCE_FRACTION <- 0.03

MIN_POSITION_RECOVERY_RATE <- 0.80

MAX_POSITION_IQR_FRACTION <- 0.03

MAX_RIGHT_FRACTION_IQR <- 0.03


# -----------------------------------------------------------------------------
# Supplement
# -----------------------------------------------------------------------------

RUN_DESEQ2_SUPPLEMENT <- TRUE

DESEQ2_RANK_METHOD <- "normalized_log1p"

# Allowed:
#
#     "normalized_log1p"
#     "vst"


RUN_MODEL_VARIANCE_SENSITIVITY <- TRUE


# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

PNG_WIDTH_IN <- 14

PNG_HEIGHT_IN <- 10.8

PNG_DPI <- 260


# -----------------------------------------------------------------------------
# Cores
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
# COLORS
# =============================================================================

COL <- list(

  var_curve =
    "#117A65",

  nb2 =
    "#1B9E77",

  nb_gap =
    "#CC1E8C",

  alpha_mu =
    "#386CB0",

  left_fill =
    "#CBE3F8",

  right_fill =
    "#DDF2D5",

  interval_fill =
    "#9E9E9E",

  anchor =
    "#000000",

  ref =
    "#E69F00",

  terminal =
    "#D95F02",

  sensitivity =
    "#7B1FA2",

  left_pt =
    "#5B8FD1",

  right_pt =
    "#43A047"
)


EVENT_LEVELS <- c(
  "Anchor",
  "Ref",
  "Terminal"
)

EVENT_COLORS <- c(

  "Anchor" =
    COL$anchor,

  "Ref" =
    COL$ref,

  "Terminal" =
    COL$terminal
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
    group_patterns) {

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
      "Count file is empty or malformed."
    )
  }


  sample_idx <- sort(
    unique(
      unlist(
        lapply(
          group_patterns,
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
      "No sample columns matched GROUP_PATTERNS."
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


  parsed <- lapply(
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
    parsed
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
      " non-finite count entries with 0."
    )

    count_mat[
      !is.finite(count_mat)
    ] <- 0
  }


  count_mat <- pmax(
    count_mat,
    0
  )


  feature_ids <- trimws(
    as.character(
      raw_df[[1L]]
    )
  )


  blank <- (
    is.na(feature_ids) |
    feature_ids == ""
  )


  if (
    any(blank)
  ) {

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


  rownames(count_mat) <-
    feature_ids


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
      "Fewer than two nonzero features remain."
    )
  }


  count_mat
}


assign_sample_groups <- function(
    sample_names,
    group_patterns) {

  out <- rep(
    NA_character_,
    length(sample_names)
  )


  for (
    group_name in names(
      group_patterns
    )
  ) {

    idx <- grep(
      group_patterns[[group_name]],
      sample_names
    )

    if (
      any(
        !is.na(
          out[idx]
        )
      )
    ) {

      stop(
        "At least one sample matched multiple group patterns."
      )
    }

    out[idx] <- group_name
  }


  if (
    any(
      is.na(out)
    )
  ) {

    stop(
      "Some selected samples could not be assigned to a group: ",
      paste(
        sample_names[
          is.na(out)
        ],
        collapse = ", "
      )
    )
  }


  factor(
    out,
    levels = names(
      group_patterns
    )
  )
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


  log1p(cpm)
}


# =============================================================================
# DESIGN-AWARE POOLED WITHIN-GROUP VARIANCE
# =============================================================================

compute_pooled_within_group_variance <- function(
    normalized_counts,
    group_labels) {

  groups <- levels(
    factor(
      group_labels
    )
  )


  sse <- rep(
    0,
    nrow(normalized_counts)
  )


  residual_df <- 0L


  for (
    g in groups
  ) {

    idx <- which(
      group_labels == g
    )


    if (
      length(idx) < 2L
    ) {

      next
    }


    mat_g <- normalized_counts[
      ,
      idx,
      drop = FALSE
    ]


    mu_g <- rowMeans(
      mat_g
    )


    residuals_g <- sweep(
      mat_g,
      1L,
      mu_g,
      "-"
    )


    sse <- sse +
      rowSums(
        residuals_g *
        residuals_g
      )


    residual_df <-
      residual_df +
      length(idx) -
      1L
  }


  if (
    residual_df < 2L
  ) {

    stop(
      "Insufficient pooled residual degrees of freedom."
    )
  }


  variance <- sse /
    residual_df


  variance[
    !is.finite(variance)
  ] <- 0


  variance <- pmax(
    variance,
    0
  )


  names(variance) <-
    rownames(
      normalized_counts
    )


  list(

    variance =
      variance,

    residual_df =
      residual_df
  )
}


compute_bootstrap_pooled_variance <- function(
    normalized_counts,
    group_labels,
    draw_by_group) {

  sse <- rep(
    0,
    nrow(normalized_counts)
  )


  residual_df <- 0L


  for (
    g in names(
      draw_by_group
    )
  ) {

    idx <- draw_by_group[[g]]


    if (
      length(idx) < 2L
    ) {

      next
    }


    mat_g <- normalized_counts[
      ,
      idx,
      drop = FALSE
    ]


    mu_g <- rowMeans(
      mat_g
    )


    residuals_g <- sweep(
      mat_g,
      1L,
      mu_g,
      "-"
    )


    sse <- sse +
      rowSums(
        residuals_g *
        residuals_g
      )


    residual_df <-
      residual_df +
      length(idx) -
      1L
  }


  if (
    residual_df < 2L
  ) {

    stop(
      "Bootstrap pooled residual df < 2."
    )
  }


  variance <- sse /
    residual_df


  variance[
    !is.finite(variance)
  ] <- 0


  variance <- pmax(
    variance,
    0
  )


  names(variance) <-
    rownames(
      normalized_counts
    )


  variance
}


# =============================================================================
# GLOBAL DESEQ2 MODEL
# =============================================================================

fit_global_deseq2_model <- function(
    count_mat,
    group_labels,
    rank_method = DESEQ2_RANK_METHOD) {

  if (
    !requireNamespace(
      "DESeq2",
      quietly = TRUE
    )
  ) {

    stop(
      "DESeq2 is required for this analysis."
    )
  }


  if (
    !requireNamespace(
      "SummarizedExperiment",
      quietly = TRUE
    )
  ) {

    stop(
      "SummarizedExperiment is required."
    )
  }


  col_data <- data.frame(

    group =
      factor(
        group_labels
      ),

    row.names =
      colnames(
        count_mat
      )
  )


  dds <- DESeq2::DESeqDataSetFromMatrix(

    countData =
      round(
        count_mat
      ),

    colData =
      col_data,

    design =
      ~ group
  )


  dds <- tryCatch(

    DESeq2::estimateSizeFactors(
      dds
    ),

    error =
      function(e) {

        message(
          "Default DESeq2 size factors failed; ",
          "using type='poscounts'."
        )

        DESeq2::estimateSizeFactors(
          dds,
          type = "poscounts"
        )
      }
  )


  normalized_counts <- DESeq2::counts(
    dds,
    normalized = TRUE
  )


  pooled_obj <- compute_pooled_within_group_variance(

    normalized_counts =
      normalized_counts,

    group_labels =
      group_labels
  )


  # ---------------------------------------------------------------------------
  # Fit the full NB model for sensitivity analysis.
  # Try parametric, then local, then mean dispersion trend.
  # ---------------------------------------------------------------------------

  fit_types <- c(
    "parametric",
    "local",
    "mean"
  )


  dds_fit <- NULL

  fit_type_used <- NA_character_


  for (
    ft in fit_types
  ) {

    candidate <- tryCatch(

      DESeq2::DESeq(
        dds,
        fitType = ft,
        quiet = TRUE
      ),

      error =
        function(e) {
          NULL
        }
    )


    if (
      !is.null(candidate)
    ) {

      dds_fit <- candidate

      fit_type_used <- ft

      break
    }
  }


  if (
    is.null(dds_fit)
  ) {

    stop(
      "DESeq2 failed with parametric, local, and mean dispersion fits."
    )
  }


  dispersions <- DESeq2::dispersions(
    dds_fit
  )


  good_disp <- (
    is.finite(dispersions) &
    dispersions >= 0
  )


  if (
    !any(good_disp)
  ) {

    stop(
      "No finite DESeq2 dispersion estimates."
    )
  }


  if (
    any(!good_disp)
  ) {

    dispersions[
      !good_disp
    ] <- median(
      dispersions[
        good_disp
      ],
      na.rm = TRUE
    )
  }


  # ---------------------------------------------------------------------------
  # Fitted normalized means.
  #
  # DESeq2 stores fitted raw means in assay "mu" after model fitting.
  # Divide by size factors to express them on normalized-count scale.
  # ---------------------------------------------------------------------------

  fitted_mu_raw <- tryCatch(

    SummarizedExperiment::assay(
      dds_fit,
      "mu"
    ),

    error =
      function(e) {
        NULL
      }
  )


  size_factors <- DESeq2::sizeFactors(
    dds_fit
  )


  if (
    is.null(fitted_mu_raw)
  ) {

    message(
      "DESeq2 fitted 'mu' assay unavailable; ",
      "using group means of normalized counts for sensitivity means."
    )


    fitted_mu_norm <- normalized_counts

  } else {

    fitted_mu_norm <- sweep(

      fitted_mu_raw,

      2L,

      size_factors,

      "/"
    )
  }


  group_levels <- levels(
    factor(
      group_labels
    )
  )


  group_mu <- matrix(

    NA_real_,

    nrow =
      nrow(count_mat),

    ncol =
      length(group_levels),

    dimnames =
      list(
        rownames(count_mat),
        group_levels
      )
  )


  model_variance <- group_mu


  for (
    g in group_levels
  ) {

    idx <- which(
      group_labels == g
    )


    mu_g <- rowMeans(
      fitted_mu_norm[
        ,
        idx,
        drop = FALSE
      ]
    )


    mu_g[
      !is.finite(mu_g)
    ] <- 0


    mu_g <- pmax(
      mu_g,
      0
    )


    group_mu[
      ,
      g
    ] <- mu_g


    model_variance[
      ,
      g
    ] <- mu_g +
      dispersions *
      mu_g^2
  }


  ranking_matrix <- NULL

  rank_method_used <- NULL


  if (
    rank_method == "vst"
  ) {

    vst_obj <- tryCatch(

      DESeq2::vst(
        dds_fit,
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
        normalized_counts
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
      normalized_counts
    )

    rank_method_used <-
      "DESeq2 log1p"
  }


  list(

    dds =
      dds_fit,

    normalized_counts =
      normalized_counts,

    ranking_matrix =
      ranking_matrix,

    rank_method_used =
      rank_method_used,

    size_factors =
      size_factors,

    dispersions =
      dispersions,

    group_mu =
      group_mu,

    model_variance =
      model_variance,

    pooled_variance =
      pooled_obj$variance,

    pooled_residual_df =
      pooled_obj$residual_df,

    dispersion_fit_type =
      fit_type_used
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
      "At least two samples required for variance."
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

  if (
    length(x) == 0L
  ) {
    return(NA_real_)
  }

  median(x)
}


safe_iqr <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (
    length(x) == 0L
  ) {
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


  if (
    length(x) == 0L
  ) {

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


suffix_sum <- function(x) {

  rev(
    cumsum(
      rev(x)
    )
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
      "PC1 undefined due to insufficient between-sample variation."
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


  names(loading) <-
    colnames(
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
# VARIANCE CURVE FROM PRECOMPUTED VARIANCE VECTOR
# =============================================================================

compute_ranked_variance_curve <- function(
    variance_vector,
    rank_order,
    spar = VAR_SPLINE_SPAR) {

  if (
    length(variance_vector) !=
    length(rank_order)
  ) {

    stop(
      "Variance vector and rank order differ in length."
    )
  }


  empirical_var <- as.numeric(
    variance_vector
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


# =============================================================================
# SECOND-DERIVATIVE CROSSINGS
# =============================================================================

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


  if (
    length(x) < 2L
  ) {

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


  if (
    length(idx) > 0L
  ) {

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


  if (
    length(exact_idx) > 0L
  ) {

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
# TERMINAL SEGMENTED-REGRESSION BREAKPOINT
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
  # The original RIGHT=Anchor->N plus equal-sized LEFT construction means the
  # relevant terminal transition must lie in the right half of the ranking.
  # This is a mathematical consequence of equal region sizes, not a 5,000-site
  # assumption.
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

  sxx <- sum(
    xd *
    xd
  )

  sy <- sum(yd)

  sxy <- sum(
    xd *
    yd
  )


  det0 <-
    n *
    sxx -
    sx^2


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
    residual0^2
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
  # Suffix sums permit exact O(N) evaluation of all hinge locations.
  # ---------------------------------------------------------------------------

  sx_suf <- suffix_sum(
    xd
  )

  sxx_suf <- suffix_sum(
    xd^2
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
  # Original matched LEFT versus terminal RIGHT variance direction.
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


  if (
    !any(eligible)
  ) {

    return(
      list(
        valid = FALSE,
        reason = "no_eligible_terminal_breakpoint",
        scan_df = scan_df
      )
    )
  }


  eligible_idx <- which(
    eligible
  )


  best_idx <- eligible_idx[
    which.max(
      bic_gain[
        eligible_idx
      ]
    )
  ]


  best <- scan_df[
    best_idx,
    ,
    drop = FALSE
  ]


  boundary_hit <- (
    best$reference_rank[
      1L
    ] ==
      candidate_min ||
    best$reference_rank[
      1L
    ] ==
      candidate_max
  )


  list(

    valid =
      !boundary_hit,

    reason =
      if (
        boundary_hit
      ) {
        "boundary_optimum"
      } else {
        NA_character_
      },

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
# MULTISCALE CONSENSUS
# =============================================================================

detect_terminal_transition <- function(
    ranked_log_variance,
    spars = TRANSITION_SPARS) {

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

        x =
          x,

        y =
          y,

        spar =
          sp
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

          reference_rank =
            if (
              !is.null(
                bp$best
              )
            ) {
              as.integer(
                bp$best$reference_rank[
                  1L
                ]
              )
            } else {
              NA_integer_
            },

          BIC_gain =
            if (
              !is.null(
                bp$best
              )
            ) {
              bp$best$BIC_gain[
                1L
              ]
            } else {
              NA_real_
            },

          gamma =
            if (
              !is.null(
                bp$best
              )
            ) {
              bp$best$gamma[
                1L
              ]
            } else {
              NA_real_
            },

          post_slope =
            if (
              !is.null(
                bp$best
              )
            ) {
              bp$best$post_slope[
                1L
              ]
            } else {
              NA_real_
            },

          terminal_minus_left_mean =
            if (
              !is.null(
                bp$best
              )
            ) {
              bp$best$terminal_minus_left_mean[
                1L
              ]
            } else {
              NA_real_
            },

          stringsAsFactors = FALSE
        )

      } else {

        b <- bp$best


        data.frame(

          spar =
            sp,

          valid =
            TRUE,

          reason =
            NA_character_,

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

        reference_rank =
          NA_integer_,

        BIC_gain =
          NA_real_,

        gamma =
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
    MIN_VALID_TRANSITION_SCALES
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


  consensus_ref <- weighted_median(
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
      consensus_ref
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
      MAX_MULTISCALE_REF_IQR_FRACTION
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
              consensus_ref
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
          consensus_ref
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

    valid_scales =
      nrow(valid_df),

    total_scales =
      length(spars),

    min_reference =
      min(refs),

    max_reference =
      max(refs)
  )
}


# =============================================================================
# ORIGINAL ANCHOR / REF / TERMINAL RULE
# =============================================================================

select_custom_interval <- function(
    zero_df,
    reference_rank,
    total_n) {

  if (
    nrow(zero_df) == 0L
  ) {

    stop(
      "No d2 zero crossings found."
    )
  }


  reference_rank <- as.integer(
    round(
      reference_rank
    )
  )


  # ---------------------------------------------------------------------------
  # Anchor must permit RIGHT=Anchor->N plus an equal-sized preceding LEFT.
  # ---------------------------------------------------------------------------

  minimum_anchor <- as.integer(
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
      minimum_anchor
  ]


  right_candidates <- zero_df$crossing_rank[
    zero_df$crossing_rank >
      reference_rank
  ]


  if (
    length(left_candidates) == 0L
  ) {

    stop(
      "No admissible left d2 crossing."
    )
  }


  if (
    length(right_candidates) == 0L
  ) {

    stop(
      "No right d2 crossing."
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


  anchor <- max(
    minimum_anchor,
    anchor
  )


  terminal <- min(
    total_n,
    terminal
  )


  if (
    !(
      anchor <
      reference_rank &&
      reference_rank <
      terminal
    )
  ) {

    stop(
      "Expected Anchor < Ref < Terminal."
    )
  }


  right_n <-
    total_n -
    anchor +
    1L


  left_start <-
    anchor -
    right_n


  if (
    left_start < 1L
  ) {

    stop(
      "Matched LEFT region cannot be constructed."
    )
  }


  list(

    anchor =
      anchor,

    ref =
      reference_rank,

    terminal =
      terminal,

    interval_min =
      anchor,

    interval_max =
      terminal
  )
}


# =============================================================================
# COMPLETE GEOMETRIC SELECTION
# =============================================================================

compute_geometric_selection <- function(
    rank_matrix,
    geometry_variance) {

  if (
    nrow(rank_matrix) !=
    length(geometry_variance)
  ) {

    stop(
      "Ranking matrix and geometry variance differ in feature count."
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

    variance_vector =
      geometry_variance,

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


  dense_df <- attr(
    variance_df,
    "dense_curve_df"
  )


  zero_df <- find_d2_zero_crossings(
    dense_df
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
# DESEQ2 MODEL-VARIANCE SENSITIVITY GEOMETRY
# =============================================================================

compute_sensitivity_geometry <- function(
    rank_matrix,
    model_variance) {

  g <- tryCatch(

    compute_geometric_selection(
      rank_matrix =
        rank_matrix,
      geometry_variance =
        model_variance
    ),

    error =
      function(e) {
        NULL
      }
  )


  if (
    is.null(g) ||
    !isTRUE(
      g$valid
    )
  ) {

    return(
      data.frame(

        SensitivityStatus =
          "FAILED",

        SensitivityAnchor =
          NA_integer_,

        SensitivityRef =
          NA_integer_,

        SensitivityTerminal =
          NA_integer_,

        SensitivityRightN =
          NA_integer_,

        stringsAsFactors = FALSE
      )
    )
  }


  data.frame(

    SensitivityStatus =
      "PASS",

    SensitivityAnchor =
      g$interval_info$anchor,

    SensitivityRef =
      g$interval_info$ref,

    SensitivityTerminal =
      g$interval_info$terminal,

    SensitivityRightN =
      g$total_n -
      g$interval_info$anchor +
      1L,

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# DOWNSTREAM FEATURE-LEVEL NB DIAGNOSTICS
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
# BOOTSTRAP COMPLETE GEOMETRY
# =============================================================================

assess_geometry_stability <- function(
    rank_matrix_arm,
    arm_group,
    global_normalized_counts,
    global_group_labels,
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
      "Cannot bootstrap invalid full geometry."
    )
  }


  N <- nrow(
    rank_matrix_arm
  )


  full_interval <-
    full_geometry$interval_info


  full_right_idx <- full_geometry$rank_order[
    seq.int(
      full_interval$anchor,
      N
    )
  ]


  full_right_features <- rownames(
    rank_matrix_arm
  )[
    full_right_idx
  ]


  group_levels <- levels(
    factor(
      global_group_labels
    )
  )


  group_indices <- setNames(
    lapply(
      group_levels,
      function(g) {

        which(
          global_group_labels == g
        )
      }
    ),
    group_levels
  )


  # ---------------------------------------------------------------------------
  # Verify focal-arm ordering against global matrix.
  # ---------------------------------------------------------------------------

  focal_global_idx <- group_indices[[
    arm_group
  ]]


  focal_names <- colnames(
    global_normalized_counts
  )[
    focal_global_idx
  ]


  if (
    !setequal(
      focal_names,
      colnames(
        rank_matrix_arm
      )
    )
  ) {

    stop(
      "Arm ranking matrix samples do not match global group samples."
    )
  }


  set.seed(
    seed
  )


  bootstrap_draws <- lapply(

    seq_len(
      bootstrap_n
    ),

    function(b) {

      setNames(
        lapply(
          group_levels,
          function(g) {

            idx <- group_indices[[g]]

            sample(
              idx,
              size = length(idx),
              replace = TRUE
            )
          }
        ),
        group_levels
      )
    }
  )


  worker <- function(b) {

    draw <- bootstrap_draws[[
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

        jaccard =
          NA_real_,

        stringsAsFactors = FALSE
      )
    }


    # -------------------------------------------------------------------------
    # Recalculate pooled variance using all 8 bootstrapped groups.
    # -------------------------------------------------------------------------

    pooled_variance_b <- tryCatch(

      compute_bootstrap_pooled_variance(

        normalized_counts =
          global_normalized_counts,

        group_labels =
          global_group_labels,

        draw_by_group =
          draw
      ),

      error =
        function(e) {
          NULL
        }
    )


    if (
      is.null(
        pooled_variance_b
      )
    ) {

      return(
        invalid_row(
          "pooled_variance_failed"
        )
      )
    }


    # -------------------------------------------------------------------------
    # Use exactly the same bootstrap draw for the focal group's PC1 ranking.
    # -------------------------------------------------------------------------

    focal_draw_global <- draw[[
      arm_group
    ]]


    focal_draw_names <- colnames(
      global_normalized_counts
    )[
      focal_draw_global
    ]


    focal_local_idx <- match(
      focal_draw_names,
      colnames(
        rank_matrix_arm
      )
    )


    if (
      any(
        is.na(
          focal_local_idx
        )
      )
    ) {

      return(
        invalid_row(
          "focal_sample_mapping_failed"
        )
      )
    }


    if (
      length(
        unique(
          focal_local_idx
        )
      ) < 2L
    ) {

      return(
        invalid_row(
          "fewer_than_2_unique_focal_samples"
        )
      )
    }


    boot_rank <- rank_matrix_arm[
      ,
      focal_local_idx,
      drop = FALSE
    ]


    g <- tryCatch(

      compute_geometric_selection(

        rank_matrix =
          boot_rank,

        geometry_variance =
          pooled_variance_b
      ),

      error =
        function(e) {
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


    interval <- g$interval_info


    right_idx <- g$rank_order[
      seq.int(
        interval$anchor,
        N
      )
    ]


    right_features <- rownames(
      rank_matrix_arm
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

    ref_recovery <- NA_real_

    anchor_recovery <- NA_real_

    terminal_recovery <- NA_real_

    ref_iqr <- NA_real_

    anchor_iqr <- NA_real_

    terminal_iqr <- NA_real_

    right_fraction_iqr <- NA_real_

    ref_q <- rep(
      NA_real_,
      3L
    )

    anchor_q <- rep(
      NA_real_,
      3L
    )

    terminal_q <- rep(
      NA_real_,
      3L
    )

    right_n_q <- rep(
      NA_real_,
      3L
    )

    median_jaccard <- NA_real_
  }


  pass_valid <- (
    is.finite(valid_rate) &&
    valid_rate >=
      MIN_BOOTSTRAP_VALID_RATE
  )


  pass_ref_recovery <- (
    is.finite(ref_recovery) &&
    ref_recovery >=
      MIN_POSITION_RECOVERY_RATE
  )


  pass_anchor_recovery <- (
    is.finite(anchor_recovery) &&
    anchor_recovery >=
      MIN_POSITION_RECOVERY_RATE
  )


  pass_terminal_recovery <- (
    is.finite(terminal_recovery) &&
    terminal_recovery >=
      MIN_POSITION_RECOVERY_RATE
  )


  pass_ref_iqr <- (
    is.finite(ref_iqr) &&
    ref_iqr /
      N <=
      MAX_POSITION_IQR_FRACTION
  )


  pass_anchor_iqr <- (
    is.finite(anchor_iqr) &&
    anchor_iqr /
      N <=
      MAX_POSITION_IQR_FRACTION
  )


  pass_terminal_iqr <- (
    is.finite(terminal_iqr) &&
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


build_main_figure <- function(
    comparison_name,
    arm_name,
    geometry,
    feature_df,
    region_summary,
    stability_summary,
    sensitivity_summary,
    pooled_residual_df,
    out_file,
    track) {

  total_n <- nrow(
    feature_df
  )


  anchor <- geometry$interval_info$anchor

  ref <- geometry$interval_info$ref

  terminal <- geometry$interval_info$terminal


  left_min <- region_summary$left_start

  left_max <- region_summary$left_end

  right_min <- region_summary$right_start

  right_max <- region_summary$right_end


  variance_df <- geometry$variance_df


  vline_df <- data.frame(

    event =
      factor(
        EVENT_LEVELS,
        levels = EVENT_LEVELS
      ),

    x =
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
        xintercept = x,
        color = event,
        linetype = event
      ),

      linewidth = 0.9
    )


  if (
    RUN_MODEL_VARIANCE_SENSITIVITY &&
    sensitivity_summary$SensitivityStatus[
      1L
    ] == "PASS"
  ) {

    p1 <- p1 +

      geom_vline(

        xintercept =
          sensitivity_summary$SensitivityRef[
            1L
          ],

        color =
          COL$sensitivity,

        linetype =
          "dotdash",

        linewidth =
          0.7
      )
  }


  p1 <- p1 +

    scale_color_manual(
      values =
        EVENT_COLORS,
      breaks =
        EVENT_LEVELS
    ) +

    scale_linetype_manual(
      values =
        EVENT_LTY,
      breaks =
        EVENT_LEVELS
    ) +

    labs(

      title =
        paste0(
          comparison_name,
          " ",
          arm_name,
          " ",
          track,
          ": pooled-variance geometry"
        ),

      subtitle =
        paste0(
          "Geometry variance estimated across all groups; residual df = ",
          pooled_residual_df,
          "; Anchor=",
          anchor,
          ", Ref=",
          ref,
          ", Terminal=",
          terminal,
          ", RIGHT n=",
          region_summary$right_n
        ),

      x =
        "Rank (ascending absolute PC1 loading)",

      y =
        "Smoothed log(1 + pooled within-group variance)",

      color =
        NULL,

      linetype =
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

    geom_vline(
      xintercept = anchor,
      linewidth = 0.8
    ) +

    geom_line(

      data =
        nb_long,

      aes(
        rank,
        value,
        color = metric
      ),

      linewidth = 0.95
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
          ": downstream NB corroboration"
        ),

      subtitle =
        paste0(
          "RIGHT-LEFT medians: NB2=",
          round(
            region_summary$diff_NB2,
            3
          ),
          "; NB2-NB1=",
          round(
            region_summary$diff_gap,
            3
          ),
          "; alpha*mu=",
          round(
            region_summary$diff_alpha,
            3
          )
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
      values =
        REGION_COLORS
    ) +

    labs(

      title =
        "Matched LEFT versus terminal RIGHT",

      subtitle =
        paste0(
          "Bootstrap status=",
          stability_summary$Status[
            1L
          ],
          "; Anchor recovery=",
          round(
            stability_summary$AnchorRecoveryRate[
              1L
            ],
            3
          ),
          "; 95% RIGHT-size interval=",
          round(
            stability_summary$RightNCI025[
              1L
            ]
          ),
          "-",
          round(
            stability_summary$RightNCI975[
              1L
            ]
          )
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
# RUN ONE TRACK
# =============================================================================

run_one_track <- function(
    comparison_name,
    arm_name,
    arm_group,
    rank_matrix,
    metric_matrix,
    pooled_variance,
    sensitivity_variance,
    global_normalized_counts,
    global_group_labels,
    pooled_residual_df,
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
    ": detecting pooled terminal transition..."
  )


  geometry <- compute_geometric_selection(

    rank_matrix =
      rank_matrix,

    geometry_variance =
      pooled_variance
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

      group =
        arm_group,

      track =
        track,

      Status =
        "NO_VALID_GEOMETRY",

      Reason =
        geometry$reason,

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
        stability = data.frame()
      )
    )
  }


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
      geometry$interval_info$anchor,

    total_n =
      geometry$total_n
  )


  # ---------------------------------------------------------------------------
  # DESeq2 model-implied variance sensitivity.
  # ---------------------------------------------------------------------------

  if (
    RUN_MODEL_VARIANCE_SENSITIVITY
  ) {

    sensitivity_summary <- compute_sensitivity_geometry(

      rank_matrix =
        rank_matrix,

      model_variance =
        sensitivity_variance
    )

  } else {

    sensitivity_summary <- data.frame(

      SensitivityStatus =
        "DISABLED",

      SensitivityAnchor =
        NA_integer_,

      SensitivityRef =
        NA_integer_,

      SensitivityTerminal =
        NA_integer_,

      SensitivityRightN =
        NA_integer_,

      stringsAsFactors = FALSE
    )
  }


  write.csv(

    sensitivity_summary,

    file.path(
      output_dir,
      paste0(
        "Table_ModelVarianceSensitivity_",
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
    " | full Anchor=",
    geometry$interval_info$anchor,
    " Ref=",
    geometry$interval_info$ref,
    " Terminal=",
    geometry$interval_info$terminal,
    " RIGHT n=",
    region_summary$right_n,
    " | bootstrapping ",
    BOOTSTRAP_N,
    " replicates..."
  )


  stability <- assess_geometry_stability(

    rank_matrix_arm =
      rank_matrix,

    arm_group =
      arm_group,

    global_normalized_counts =
      global_normalized_counts,

    global_group_labels =
      global_group_labels,

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

      group =
        arm_group,

      track =
        track,

      stringsAsFactors = FALSE
    ),

    stability$stability_summary
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


  rank_output <- feature_df %>%

    left_join(
      geometry$variance_df,
      by = "rank",
      suffix =
        c(
          "_arm",
          "_pooled_geometry"
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


  summary_row <- data.frame(

    comp =
      comparison_name,

    arm =
      arm_name,

    group =
      arm_group,

    track =
      track,

    rank_method =
      rank_method,

    metric_matrix =
      metric_name,

    GeometryVariance =
      "DESeq2-normalized pooled within-group residual variance",

    PooledResidualDF =
      pooled_residual_df,

    Status =
      stability_summary$Status[
        1L
      ],

    Anchor =
      geometry$interval_info$anchor,

    Ref =
      geometry$interval_info$ref,

    Terminal =
      geometry$interval_info$terminal,

    RightSize =
      region_summary$right_n,

    RightFraction =
      region_summary$right_n /
      geometry$total_n,

    TransitionScaleIQR =
      geometry$detector$scale_iqr,

    TransitionScaleIQRFraction =
      geometry$detector$scale_iqr_fraction,

    TransitionMedianBICGain =
      geometry$detector$median_BIC_gain,

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

    SensitivityStatus =
      sensitivity_summary$SensitivityStatus[
        1L
      ],

    SensitivityAnchor =
      sensitivity_summary$SensitivityAnchor[
        1L
      ],

    SensitivityRef =
      sensitivity_summary$SensitivityRef[
        1L
      ],

    SensitivityTerminal =
      sensitivity_summary$SensitivityTerminal[
        1L
      ],

    SensitivityRightN =
      sensitivity_summary$SensitivityRightN[
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

    sensitivity_summary =
      sensitivity_summary,

    pooled_residual_df =
      pooled_residual_df,

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

    track =
      track
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
    geometry$interval_info$anchor,
    " Ref=",
    geometry$interval_info$ref,
    " Terminal=",
    geometry$interval_info$terminal,
    " | RIGHT n=",
    region_summary$right_n,
    " | bootstrap Anchor median=",
    round(
      stability_summary$AnchorMedian[
        1L
      ]
    ),
    " | model-sensitivity Ref=",
    sensitivity_summary$SensitivityRef[
      1L
    ]
  )


  list(

    summary =
      summary_row,

    stability =
      stability_summary
  )
}


# =============================================================================
# RUN
# =============================================================================

count_mat <- read_count_matrix(

  path =
    COUNT_FILE,

  group_patterns =
    GROUP_PATTERNS
)


group_labels <- assign_sample_groups(

  sample_names =
    colnames(count_mat),

  group_patterns =
    GROUP_PATTERNS
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
  "Experimental groups: ",
  paste(
    levels(group_labels),
    collapse = ", "
  )
)


message(
  "Fitting one global design-aware DESeq2 model: design ~ group ..."
)


global_model <- fit_global_deseq2_model(

  count_mat =
    count_mat,

  group_labels =
    group_labels,

  rank_method =
    DESEQ2_RANK_METHOD
)


message(
  "DESeq2 dispersion fit type: ",
  global_model$dispersion_fit_type
)


message(
  "Pooled within-group residual df: ",
  global_model$pooled_residual_df
)


message(
  "Fixed 5,000-feature reference: REMOVED"
)


message(
  "Primary geometry: pooled within-group variance across all groups"
)


message(
  "DESeq2 model-implied NB variance: sensitivity analysis only"
)


message(
  "Bootstrap replicates per track: ",
  BOOTSTRAP_N
)


message(
  "Bootstrap cores: ",
  BOOTSTRAP_CORES
)


# -----------------------------------------------------------------------------
# Save global DESeq2 / pooled-variance information
# -----------------------------------------------------------------------------

write.csv(

  data.frame(

    feature_id =
      rownames(count_mat),

    pooled_within_group_variance =
      global_model$pooled_variance,

    DESeq2_dispersion =
      global_model$dispersions,

    stringsAsFactors = FALSE
  ),

  file.path(
    OUT_ROOT,
    "Table_Global_Variance_Estimates.csv"
  ),

  row.names = FALSE
)


write.csv(

  data.frame(

    sample =
      colnames(count_mat),

    group =
      as.character(
        group_labels
      ),

    size_factor =
      as.numeric(
        global_model$size_factors
      ),

    stringsAsFactors = FALSE
  ),

  file.path(
    OUT_ROOT,
    "Table_Global_Samples_SizeFactors.csv"
  ),

  row.names = FALSE
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


  comp_groups <- COMPARISONS[[
    comparison_name
  ]]


  for (
    arm_name in c(
      "control",
      "treatment"
    )
  ) {

    arm_group <- comp_groups[[
      arm_name
    ]]


    global_idx <- which(
      group_labels ==
      arm_group
    )


    if (
      length(global_idx) < 2L
    ) {

      stop(
        "Not enough samples for group ",
        arm_group
      )
    }


    count_mat_arm <- count_mat[
      ,
      global_idx,
      drop = FALSE
    ]


    normalized_arm <- global_model$normalized_counts[
      ,
      global_idx,
      drop = FALSE
    ]


    model_variance_arm <- global_model$model_variance[
      ,
      arm_group
    ]


    # =========================================================================
    # MAIN TRACK
    #
    # Rank:
    #     arm-specific CPM log1p
    #
    # Geometry:
    #     pooled within-group variance estimated using ALL groups
    #
    # Downstream NB diagnostics:
    #     raw counts from the focal arm
    # =========================================================================

    main_rank_matrix <- normalize_cpm_log1p(
      count_mat_arm
    )


    geometry_seed <- seed_from_label(
      paste(
        comparison_name,
        arm_name,
        "shared_geometry",
        sep = "_"
      )
    )


    main_res <- run_one_track(

      comparison_name =
        comparison_name,

      arm_name =
        arm_name,

      arm_group =
        arm_group,

      rank_matrix =
        main_rank_matrix,

      metric_matrix =
        count_mat_arm,

      pooled_variance =
        global_model$pooled_variance,

      sensitivity_variance =
        model_variance_arm,

      global_normalized_counts =
        global_model$normalized_counts,

      global_group_labels =
        group_labels,

      pooled_residual_df =
        global_model$pooled_residual_df,

      output_dir =
        comp_dir,

      track =
        "Main",

      rank_method =
        "CPM log1p",

      metric_name =
        "raw_counts",

      seed =
        geometry_seed
    )


    overall_rows[[
      length(overall_rows) +
      1L
    ]] <- main_res$summary


    if (
      nrow(
        main_res$stability
      ) > 0L
    ) {

      overall_stability_rows[[
        length(overall_stability_rows) +
        1L
      ]] <- main_res$stability
    }


    # =========================================================================
    # DESEQ2 SUPPLEMENTARY TRACK
    #
    # Rank:
    #     DESeq2 normalized log1p or VST
    #
    # Geometry:
    #     identical pooled design-aware variance estimator
    #
    # Downstream NB diagnostics:
    #     DESeq2-normalized counts
    # =========================================================================

    if (
      RUN_DESEQ2_SUPPLEMENT
    ) {

      deseq2_rank_arm <- global_model$ranking_matrix[
        ,
        global_idx,
        drop = FALSE
      ]


      deseq2_res <- run_one_track(

        comparison_name =
          comparison_name,

        arm_name =
          arm_name,

        arm_group =
          arm_group,

        rank_matrix =
          deseq2_rank_arm,

        metric_matrix =
          normalized_arm,

        pooled_variance =
          global_model$pooled_variance,

        sensitivity_variance =
          model_variance_arm,

        global_normalized_counts =
          global_model$normalized_counts,

        global_group_labels =
          group_labels,

        pooled_residual_df =
          global_model$pooled_residual_df,

        output_dir =
          comp_dir,

        track =
          "DESeq2",

        rank_method =
          global_model$rank_method_used,

        metric_name =
          "deseq2_normalized_counts",

        # Same resampling seed as Main so both tracks see the same
        # biological bootstrap draws.
        seed =
          geometry_seed
      )


      overall_rows[[
        length(overall_rows) +
        1L
      ]] <- deseq2_res$summary


      if (
        nrow(
          deseq2_res$stability
        ) > 0L
      ) {

        overall_stability_rows[[
          length(overall_stability_rows) +
          1L
        ]] <- deseq2_res$stability
      }
    }
  }
}


# =============================================================================
# OVERALL TABLES
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
  "No fixed 5,000-site reference was used."
)


message(
  "Variance geometry was estimated using all experimental groups."
)


message(
  "Between-group mean differences were removed before pooling variance."
)


message(
  "Primary geometry therefore uses substantially more residual df than ",
  "arm-specific empirical variance."
)


message(
  "Arm-specific PC1 ranking was preserved."
)


message(
  "RIGHT remains Anchor -> terminal rank N."
)


message(
  "LEFT remains the immediately preceding equal-sized block."
)


message(
  "DESeq2 model-implied NB variance was used only as sensitivity analysis."
)


message(
  "NB2-related LEFT/RIGHT diagnostics remain downstream of primary geometry."
)


message(
  "Outputs written to: ",
  OUT_ROOT
)


message(
  "============================================================"
)
