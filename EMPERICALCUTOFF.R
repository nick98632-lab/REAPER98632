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
# SHARED PC1–NEGATIVE-BINOMIAL DIVERGENCE BOUNDARIES
# =============================================================================
#
# STUDY AXIS
# ----------
# Within each experimental arm, features are ranked from the smallest to the
# largest contribution to variance represented by principal component 1 (PC1).
#
# For a centered sample-by-feature matrix X with singular value decomposition
#
#     X = U D V^T,
#
# the feature-specific variance contribution of PC1 is
#
#                     d1^2 * v_i1^2
#     P_i =          ---------------- ,
#                         n - 1
#
# where d1 is the first singular value, v_i1 is feature i's PC1 loading, and n
# is the number of samples in that arm. Ranking by P_i is equivalent to ranking
# by |v_i1| within an arm, while P_i additionally retains the variance scale
# implied by the singular value decomposition.
#
# The resulting ranked axis is interpreted as
#
#     REMAINDER  ------------------------------------------>  LEADING EDGE
#     low PC1 contribution                                  high PC1 contribution
#
#
# DESIGN-AWARE SEQUENCING VARIANCE
# --------------------------------
# Counts are normalized once with DESeq2 size factors across the complete
# experiment. For each feature, a pooled within-group variance is then computed
# after removing each experimental group's mean:
#
#                  sum_g sum_{j in g} (y_ij - ybar_ig)^2
#     V_i =        --------------------------------------- .
#                           sum_g (n_g - 1)
#
# Thus all experimental groups contribute to the sequencing-variance estimate,
# while differences between experimental-group means do not inflate V_i.
#
#
# NEGATIVE-BINOMIAL EXCESS-VARIANCE SIGNAL
# ----------------------------------------
# For feature i in experimental group g:
#
#     E_ig = max(V_i - mu_ig, 0),
#
# where mu_ig is the normalized group mean. E_ig is the extra-Poisson variance
# component used to quantify the sequencing variance that exceeds the linear
# Poisson mean term.
#
#
# DIMENSIONLESS PC1–NB DIVERGENCE
# ------------------------------
# PC1 variance contribution and sequencing excess variance are not directly
# commensurate. Each is therefore log-transformed and robustly standardized:
#
#     P*_ig(r) = robustZ[ log(1 + P_ig(r)) ]
#
#     E*_ig(r) = robustZ[ log(1 + E_ig(r)) ]
#
# and their rank-wise divergence is
#
#     D_g(r) = E*_g(r) - P*_g(r).
#
# The shared experimental divergence trajectory is the rank-wise median:
#
#     D_cons(r) = median_g D_g(r).
#
#
# SHARED REMAINDER AND LEADING-EDGE BOUNDARIES
# --------------------------------------------
# A generalized cross-validation smoothing spline is fit to D_cons(r). Its
# first and second derivatives are evaluated on the complete normalized rank
# axis.
#
# A transition episode begins where curvature changes from non-positive to
# positive:
#
#     D_cons''(r):  -  ->  +
#
# and continues until the next positive-to-non-positive curvature crossing (or
# the end of the applicable rank domain). Episode strength is the increase in
# the smoothed divergence trajectory across that acceleration episode.
#
# The REMAINDER_BOUNDARY is the strongest positive-divergence acceleration
# episode whose onset lies in the lower-rank half of the ordered axis.
#
# The LEADING_EDGE_BOUNDARY is the strongest positive-divergence acceleration
# episode whose onset lies in the upper-rank half of the ordered axis.
#
# The same two shared rank boundaries are applied to every arm:
#
#     REMAINDER:
#         ranks < REMAINDER_BOUNDARY
#
#     TRANSITION INTERVAL:
#         REMAINDER_BOUNDARY <= rank < LEADING_EDGE_BOUNDARY
#
#     LEADING EDGE:
#         ranks >= LEADING_EDGE_BOUNDARY
#
#
# BOOTSTRAP UNCERTAINTY
# ---------------------
# Biological samples are resampled with replacement within each experimental
# group. Every bootstrap replicate recomputes:
#
#     pooled within-group variance
#     -> arm-specific PC1 contribution and ranking
#     -> arm-specific dimensionless divergence
#     -> shared consensus divergence
#     -> both derivative-defined boundaries.
#
# The manuscript reports the full-data boundary estimate, bootstrap median,
# and bootstrap interquartile range (25th-75th percentiles).
#
#
# REGIONAL NB CHARACTERIZATION
# ----------------------------
# After the shared boundaries are determined, the three rank regions are
# summarized using:
#
#     NB2_i      = log(1 + E_i)
#
#     NB2-NB1_i  = log(1 + E_i) - log(1 + mu_i)
#
#     alpha_i    = E_i / mu_i^2
#
#     alpha*mu_i = log(1 + alpha_i * mu_i).
#
# These are descriptive moment-based diagnostics of the mean-variance
# relationship and are not likelihood-ratio tests.
#
#
# FIGURES
# -------
# Mathematical equations are printed directly on the consensus and arm-level
# figures so that the statistical construction is explicit in the figure
# itself.
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <-
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <-
  "/root/REAPER98632/exports/manuscript_pc1_nb_divergence_final"


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
  RT0_ZT6  = list(control = "RT0", treatment = "ZT6"),
  RT2_ZT8  = list(control = "RT2", treatment = "ZT8"),
  RT4_ZT10 = list(control = "RT4", treatment = "ZT10"),
  RT8_ZT14 = list(control = "RT8", treatment = "ZT14")
)


# A small endpoint exclusion protects derivative estimation from spline boundary
# behavior. It does not target any biological cutoff location.
ENDPOINT_GUARD_FRACTION <- 0.02


BOOTSTRAP_N <- 500L
BOOTSTRAP_SEED <- 20260820L


DETECTED_CORES <- suppressWarnings(
  parallel::detectCores(logical = FALSE)
)

if (
  is.na(DETECTED_CORES) ||
  DETECTED_CORES < 2L
) {
  BOOTSTRAP_CORES <- 1L
} else {
  BOOTSTRAP_CORES <- min(4L, DETECTED_CORES - 1L)
}


PNG_WIDTH_IN <- 15
PNG_HEIGHT_IN <- 14
PNG_DPI <- 260


dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# COLORS
# =============================================================================

COL <- list(
  pc1          = "#386CB0",
  excess       = "#D95F02",
  divergence   = "#7B1FA2",
  derivative1  = "#1B9E77",
  derivative2  = "#CC1E8C",
  remainder    = "#DCEAF7",
  transition   = "#ECECEC",
  leading      = "#DDF2D5",
  rem_line     = "#2166AC",
  lead_line    = "#1B7837",
  nb2          = "#1B9E77",
  nb_gap       = "#CC1E8C",
  alpha_mu     = "#386CB0"
)


# =============================================================================
# INPUT
# =============================================================================

read_count_matrix <- function(
    path,
    group_patterns) {

  if (!file.exists(path)) {
    stop("Count file does not exist: ", path)
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
    stop("Count file is empty or malformed: ", path)
  }


  sample_idx <- sort(
    unique(
      unlist(
        lapply(
          group_patterns,
          function(pattern) {
            grep(pattern, colnames(raw_df))
          }
        )
      )
    )
  )


  if (length(sample_idx) == 0L) {
    stop("No sample columns matched GROUP_PATTERNS.")
  }


  if (1L %in% sample_idx) {
    stop(
      "Column 1 matched a sample pattern; column 1 must contain feature IDs."
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


  colnames(count_mat) <- colnames(count_df)
  storage.mode(count_mat) <- "numeric"


  nonfinite_n <- sum(
    !is.finite(count_mat)
  )


  if (nonfinite_n > 0L) {
    message(
      "Replacing ",
      nonfinite_n,
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


  keep <- rowSums(
    count_mat
  ) > 0


  count_mat <- count_mat[
    keep,
    ,
    drop = FALSE
  ]


  if (nrow(count_mat) < 2L) {
    stop("Fewer than two nonzero features remain after filtering.")
  }


  count_mat
}


assign_sample_groups <- function(
    sample_names,
    group_patterns) {

  assigned <- rep(
    NA_character_,
    length(sample_names)
  )


  for (
    group_name in names(group_patterns)
  ) {

    idx <- grep(
      group_patterns[[group_name]],
      sample_names
    )


    if (
      length(idx) > 0L &&
      any(!is.na(assigned[idx]))
    ) {
      stop("At least one sample matched more than one group pattern.")
    }


    assigned[idx] <- group_name
  }


  if (any(is.na(assigned))) {
    stop(
      "Unassigned sample columns: ",
      paste(
        sample_names[is.na(assigned)],
        collapse = ", "
      )
    )
  }


  factor(
    assigned,
    levels = names(group_patterns)
  )
}


# =============================================================================
# NORMALIZATION
# =============================================================================

normalize_cpm_log1p <- function(
    count_mat) {

  lib_sizes <- colSums(
    count_mat,
    na.rm = TRUE
  )


  lib_sizes[
    !is.finite(lib_sizes) |
    lib_sizes <= 0
  ] <- 1


  cpm <- sweep(
    count_mat,
    2L,
    lib_sizes / 1e6,
    "/"
  )


  log1p(cpm)
}


normalize_deseq2 <- function(
    count_mat,
    group_labels) {

  if (
    !requireNamespace(
      "DESeq2",
      quietly = TRUE
    )
  ) {
    stop(
      "DESeq2 is required. Install DESeq2 before running this manuscript script."
    )
  }


  col_data <- data.frame(
    group = factor(group_labels),
    row.names = colnames(count_mat)
  )


  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = col_data,
    design = ~ group
  )


  dds <- tryCatch(
    DESeq2::estimateSizeFactors(dds),
    error = function(e) {
      message(
        "Default DESeq2 size-factor estimation failed; using type='poscounts'."
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


  list(
    normalized_counts = normalized_counts,
    size_factors = DESeq2::sizeFactors(dds)
  )
}


# =============================================================================
# BASIC NUMERICAL HELPERS
# =============================================================================

robust_z <- function(x) {

  x <- as.numeric(x)

  finite <- is.finite(x)

  if (!any(finite)) {
    return(rep(0, length(x)))
  }


  med <- median(
    x[finite],
    na.rm = TRUE
  )


  scale_value <- stats::mad(
    x[finite],
    center = med,
    constant = 1.4826,
    na.rm = TRUE
  )


  if (
    !is.finite(scale_value) ||
    scale_value <= .Machine$double.eps
  ) {

    scale_value <- stats::IQR(
      x[finite],
      na.rm = TRUE
    ) / 1.349
  }


  if (
    !is.finite(scale_value) ||
    scale_value <= .Machine$double.eps
  ) {

    scale_value <- stats::sd(
      x[finite],
      na.rm = TRUE
    )
  }


  if (
    !is.finite(scale_value) ||
    scale_value <= .Machine$double.eps
  ) {
    scale_value <- 1
  }


  out <- rep(
    0,
    length(x)
  )


  out[finite] <- (
    x[finite] -
    med
  ) /
    scale_value


  out
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


safe_quantile <- function(
    x,
    p) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  as.numeric(
    stats::quantile(
      x,
      probs = p,
      names = FALSE,
      na.rm = TRUE,
      type = 7
    )
  )
}


interp_at <- function(
    x,
    y,
    xout) {

  as.numeric(
    stats::approx(
      x = x,
      y = y,
      xout = xout,
      rule = 2
    )$y
  )
}


# =============================================================================
# DESIGN-AWARE POOLED WITHIN-GROUP VARIANCE
# =============================================================================

compute_pooled_within_group_variance <- function(
    normalized_counts,
    group_labels) {

  group_levels <- levels(
    factor(group_labels)
  )


  sse <- rep(
    0,
    nrow(normalized_counts)
  )


  residual_df <- 0L


  for (
    g in group_levels
  ) {

    idx <- which(
      group_labels == g
    )


    if (length(idx) < 2L) {
      next
    }


    xg <- normalized_counts[
      ,
      idx,
      drop = FALSE
    ]


    mu_g <- rowMeans(
      xg
    )


    resid_g <- sweep(
      xg,
      1L,
      mu_g,
      "-"
    )


    sse <- sse +
      rowSums(
        resid_g *
        resid_g
      )


    residual_df <- residual_df +
      length(idx) -
      1L
  }


  if (residual_df < 2L) {
    stop("Pooled within-group residual degrees of freedom < 2.")
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


  names(variance) <- rownames(
    normalized_counts
  )


  list(
    variance = variance,
    residual_df = residual_df
  )
}


compute_pooled_variance_from_draws <- function(
    normalized_counts,
    draws_by_group) {

  sse <- rep(
    0,
    nrow(normalized_counts)
  )


  residual_df <- 0L


  for (
    g in names(draws_by_group)
  ) {

    idx <- draws_by_group[[g]]


    if (length(idx) < 2L) {
      next
    }


    xg <- normalized_counts[
      ,
      idx,
      drop = FALSE
    ]


    mu_g <- rowMeans(
      xg
    )


    resid_g <- sweep(
      xg,
      1L,
      mu_g,
      "-"
    )


    sse <- sse +
      rowSums(
        resid_g *
        resid_g
      )


    residual_df <- residual_df +
      length(idx) -
      1L
  }


  if (residual_df < 2L) {
    stop("Bootstrap pooled residual degrees of freedom < 2.")
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


  names(variance) <- rownames(
    normalized_counts
  )


  variance
}


# =============================================================================
# PC1 VARIANCE CONTRIBUTION FROM SVD GEOMETRY
# =============================================================================

compute_pc1_variance_contribution <- function(
    rank_matrix_arm) {

  if (
    nrow(rank_matrix_arm) < 2L ||
    ncol(rank_matrix_arm) < 2L
  ) {
    stop("PC1 requires at least two features and two samples.")
  }


  # Samples x features
  X <- t(
    rank_matrix_arm
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
    lambda1 <= .Machine$double.eps
  ) {
    stop("PC1 is undefined because between-sample variation is insufficient.")
  }


  u1 <- eig$vectors[
    ,
    1L
  ]


  d1 <- sqrt(
    lambda1
  )


  v1 <- as.numeric(
    crossprod(
      Xc,
      u1
    )
  ) /
    d1


  names(v1) <- colnames(
    Xc
  )


  v1[
    !is.finite(v1)
  ] <- 0


  n_samples <- nrow(
    Xc
  )


  pc1_variance_contribution <- (
    lambda1 *
    v1^2
  ) /
    (
      n_samples -
      1L
    )


  total_feature_variance <- colSums(
    Xc^2
  ) /
    (
      n_samples -
      1L
    )


  fraction_feature_variance_pc1 <- rep(
    0,
    length(v1)
  )


  positive_total <- (
    is.finite(total_feature_variance) &
    total_feature_variance > 0
  )


  fraction_feature_variance_pc1[
    positive_total
  ] <- pc1_variance_contribution[
    positive_total
  ] /
    total_feature_variance[
      positive_total
    ]


  total_ss <- sum(
    Xc^2
  )


  pc1_fraction_total <- if (
    is.finite(total_ss) &&
    total_ss > 0
  ) {
    lambda1 / total_ss
  } else {
    NA_real_
  }


  list(
    loading = v1,
    abs_loading = abs(v1),
    singular_value = d1,
    eigenvalue = lambda1,
    contribution = pc1_variance_contribution,
    feature_fraction = fraction_feature_variance_pc1,
    pc1_fraction_total = pc1_fraction_total
  )
}


# =============================================================================
# GROUP-SPECIFIC PC1–NB DIVERGENCE CURVE
# =============================================================================

compute_one_group_curve <- function(
    group_name,
    rank_matrix_arm,
    normalized_counts_arm,
    pooled_variance) {

  pc1 <- compute_pc1_variance_contribution(
    rank_matrix_arm
  )


  mu <- rowMeans(
    normalized_counts_arm
  )


  mu[
    !is.finite(mu)
  ] <- 0


  mu <- pmax(
    mu,
    0
  )


  excess <- pmax(
    pooled_variance -
    mu,
    0
  )


  rank_order <- order(
    pc1$contribution,
    decreasing = FALSE
  )


  P <- pc1$contribution[
    rank_order
  ]


  E <- excess[
    rank_order
  ]


  mu_ranked <- mu[
    rank_order
  ]


  V_ranked <- pooled_variance[
    rank_order
  ]


  P_star <- robust_z(
    log1p(P)
  )


  E_star <- robust_z(
    log1p(E)
  )


  divergence <- E_star -
    P_star


  alpha <- rep(
    0,
    length(mu_ranked)
  )


  positive_mu <- mu_ranked > 0


  alpha[
    positive_mu
  ] <- E[
    positive_mu
  ] /
    (
      mu_ranked[
        positive_mu
      ]^2
    )


  feature_df <- data.frame(
    group = group_name,
    rank = seq_along(rank_order),
    feature_id = rownames(rank_matrix_arm)[rank_order],
    pc1_loading = pc1$loading[rank_order],
    abs_pc1_loading = pc1$abs_loading[rank_order],
    pc1_variance_contribution = P,
    pc1_feature_variance_fraction = pc1$feature_fraction[rank_order],
    group_mean_normalized = mu_ranked,
    pooled_within_group_variance = V_ranked,
    excess_variance = E,
    P_star = P_star,
    E_star = E_star,
    divergence = divergence,
    NB2 = log1p(E),
    NB2_NB1 = log1p(E) - log1p(mu_ranked),
    alpha_mu = log1p(alpha * mu_ranked),
    stringsAsFactors = FALSE
  )


  pc1_summary <- data.frame(
    group = group_name,
    n_samples = ncol(rank_matrix_arm),
    singular_value_1 = pc1$singular_value,
    eigenvalue_1 = pc1$eigenvalue,
    pc1_fraction_total_variance = pc1$pc1_fraction_total,
    stringsAsFactors = FALSE
  )


  list(
    feature_df = feature_df,
    pc1_summary = pc1_summary
  )
}


compute_all_group_curves <- function(
    rank_matrix_all,
    normalized_counts,
    pooled_variance,
    group_labels) {

  groups <- levels(
    factor(group_labels)
  )


  curve_list <- vector(
    "list",
    length(groups)
  )


  names(curve_list) <- groups


  pc1_rows <- vector(
    "list",
    length(groups)
  )


  names(pc1_rows) <- groups


  for (
    g in groups
  ) {

    idx <- which(
      group_labels == g
    )


    one <- compute_one_group_curve(
      group_name = g,
      rank_matrix_arm = rank_matrix_all[
        ,
        idx,
        drop = FALSE
      ],
      normalized_counts_arm = normalized_counts[
        ,
        idx,
        drop = FALSE
      ],
      pooled_variance = pooled_variance
    )


    curve_list[[g]] <- one$feature_df
    pc1_rows[[g]] <- one$pc1_summary
  }


  list(
    curves = curve_list,
    pc1_summary = bind_rows(pc1_rows)
  )
}


# =============================================================================
# SHARED CONSENSUS DIVERGENCE
# =============================================================================

compute_consensus_curve <- function(
    group_curves) {

  groups <- names(
    group_curves
  )


  N_values <- vapply(
    group_curves,
    nrow,
    integer(1)
  )


  if (
    length(
      unique(N_values)
    ) != 1L
  ) {
    stop("All group curves must contain the same number of ranked features.")
  }


  N <- N_values[
    1L
  ]


  P_mat <- do.call(
    cbind,
    lapply(
      group_curves,
      function(df) {
        df$P_star
      }
    )
  )


  E_mat <- do.call(
    cbind,
    lapply(
      group_curves,
      function(df) {
        df$E_star
      }
    )
  )


  D_mat <- do.call(
    cbind,
    lapply(
      group_curves,
      function(df) {
        df$divergence
      }
    )
  )


  colnames(P_mat) <- groups
  colnames(E_mat) <- groups
  colnames(D_mat) <- groups


  data.frame(
    rank = seq_len(N),
    consensus_P_star = apply(
      P_mat,
      1L,
      median,
      na.rm = TRUE
    ),
    consensus_E_star = apply(
      E_mat,
      1L,
      median,
      na.rm = TRUE
    ),
    consensus_divergence = apply(
      D_mat,
      1L,
      median,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# CALCULUS-BASED TWO-BOUNDARY DETECTOR
# =============================================================================

find_curvature_crossings <- function(
    x,
    d2) {

  ok <- (
    is.finite(x) &
    is.finite(d2)
  )


  x <- x[ok]
  d2 <- d2[ok]


  if (length(x) < 2L) {
    return(
      data.frame(
        x = numeric(0),
        direction = character(0),
        stringsAsFactors = FALSE
      )
    )
  }


  a <- d2[
    -length(d2)
  ]


  b <- d2[
    -1L
  ]


  xa <- x[
    -length(x)
  ]


  xb <- x[
    -1L
  ]


  idx_neg_pos <- which(
    a <= 0 &
    b > 0
  )


  idx_pos_neg <- which(
    a >= 0 &
    b < 0
  )


  interpolate_crossing <- function(idx) {

    if (length(idx) == 0L) {
      return(numeric(0))
    }


    denom <- abs(a[idx]) +
      abs(b[idx])


    frac <- rep(
      0.5,
      length(idx)
    )


    positive_denom <- denom > 0


    frac[
      positive_denom
    ] <- abs(
      a[
        idx[
          positive_denom
        ]
      ]
    ) /
      denom[
        positive_denom
      ]


    xa[idx] +
      frac *
      (
        xb[idx] -
        xa[idx]
      )
  }


  out <- bind_rows(
    data.frame(
      x = interpolate_crossing(
        idx_neg_pos
      ),
      direction = "negative_to_positive",
      stringsAsFactors = FALSE
    ),
    data.frame(
      x = interpolate_crossing(
        idx_pos_neg
      ),
      direction = "positive_to_negative",
      stringsAsFactors = FALSE
    )
  )


  out %>%
    filter(is.finite(x)) %>%
    arrange(x)
}


build_acceleration_episodes <- function(
    crossing_df,
    x,
    smooth_y,
    d1,
    endpoint_guard_fraction = ENDPOINT_GUARD_FRACTION) {

  if (nrow(crossing_df) == 0L) {
    return(
      data.frame()
    )
  }


  guard_low <- endpoint_guard_fraction

  guard_high <- 1 -
    endpoint_guard_fraction


  onsets <- crossing_df %>%
    filter(
      direction == "negative_to_positive",
      x >= guard_low,
      x <= guard_high
    )


  ends <- crossing_df %>%
    filter(
      direction == "positive_to_negative"
    )


  if (nrow(onsets) == 0L) {
    return(
      data.frame()
    )
  }


  rows <- vector(
    "list",
    nrow(onsets)
  )


  for (
    i in seq_len(
      nrow(onsets)
    )
  ) {

    onset_x <- onsets$x[
      i
    ]


    next_end <- ends$x[
      ends$x >
        onset_x
    ]


    if (length(next_end) > 0L) {
      end_x <- min(next_end)
    } else {
      end_x <- 1
    }


    onset_y <- interp_at(
      x,
      smooth_y,
      onset_x
    )


    end_y <- interp_at(
      x,
      smooth_y,
      end_x
    )


    x_idx <- which(
      x >= onset_x &
      x <= end_x
    )


    if (length(x_idx) == 0L) {
      peak_d1 <- NA_real_
    } else {
      peak_d1 <- max(
        d1[x_idx],
        na.rm = TRUE
      )
    }


    rise <- end_y -
      onset_y


    rows[[i]] <- data.frame(
      onset_x = onset_x,
      end_x = end_x,
      rise = rise,
      peak_d1 = peak_d1,
      valid = (
        is.finite(rise) &&
        rise > 0 &&
        is.finite(peak_d1) &&
        peak_d1 > 0
      ),
      stringsAsFactors = FALSE
    )
  }


  bind_rows(rows)
}


detect_shared_boundaries <- function(
    consensus_df) {

  N <- nrow(
    consensus_df
  )


  if (N < 20L) {
    stop("Too few ranked features for boundary detection.")
  }


  x <- (
    consensus_df$rank -
    1
  ) /
    (
      N -
      1
    )


  y <- consensus_df$consensus_divergence


  if (any(!is.finite(y))) {
    stop("Consensus divergence contains non-finite values.")
  }


  # GCV selects the smoothing level from the observed consensus trajectory.
  spline_fit <- stats::smooth.spline(
    x = x,
    y = y,
    cv = FALSE
  )


  smooth_y <- as.numeric(
    stats::predict(
      spline_fit,
      x = x,
      deriv = 0
    )$y
  )


  d1 <- as.numeric(
    stats::predict(
      spline_fit,
      x = x,
      deriv = 1
    )$y
  )


  d2 <- as.numeric(
    stats::predict(
      spline_fit,
      x = x,
      deriv = 2
    )$y
  )


  crossings <- find_curvature_crossings(
    x = x,
    d2 = d2
  )


  episodes <- build_acceleration_episodes(
    crossing_df = crossings,
    x = x,
    smooth_y = smooth_y,
    d1 = d1
  )


  if (
    nrow(episodes) == 0L ||
    !any(episodes$valid)
  ) {
    return(
      list(
        valid = FALSE,
        reason = "no_valid_positive_acceleration_episodes",
        curve = data.frame(
          rank = consensus_df$rank,
          x = x,
          smooth_divergence = smooth_y,
          derivative_1 = d1,
          derivative_2 = d2,
          stringsAsFactors = FALSE
        ),
        crossings = crossings,
        episodes = episodes,
        spline = spline_fit
      )
    )
  }


  episodes <- episodes %>%
    mutate(
      onset_rank = 1 +
        onset_x *
        (
          N -
          1
        ),
      end_rank = 1 +
        end_x *
        (
          N -
          1
        ),
      domain = ifelse(
        onset_x <= 0.5,
        "REMAINDER_SIDE",
        "LEADING_EDGE_SIDE"
      )
    )


  remainder_candidates <- episodes %>%
    filter(
      valid,
      domain == "REMAINDER_SIDE"
    )


  leading_candidates <- episodes %>%
    filter(
      valid,
      domain == "LEADING_EDGE_SIDE"
    )


  if (nrow(remainder_candidates) == 0L) {
    return(
      list(
        valid = FALSE,
        reason = "no_valid_remainder_side_transition",
        curve = data.frame(
          rank = consensus_df$rank,
          x = x,
          smooth_divergence = smooth_y,
          derivative_1 = d1,
          derivative_2 = d2,
          stringsAsFactors = FALSE
        ),
        crossings = crossings,
        episodes = episodes,
        spline = spline_fit
      )
    )
  }


  if (nrow(leading_candidates) == 0L) {
    return(
      list(
        valid = FALSE,
        reason = "no_valid_leading_edge_side_transition",
        curve = data.frame(
          rank = consensus_df$rank,
          x = x,
          smooth_divergence = smooth_y,
          derivative_1 = d1,
          derivative_2 = d2,
          stringsAsFactors = FALSE
        ),
        crossings = crossings,
        episodes = episodes,
        spline = spline_fit
      )
    )
  }


  # The strongest sustained increase in smoothed divergence identifies each
  # side's transition episode. Ties are resolved by larger peak slope.
  remainder_candidates <- remainder_candidates %>%
    arrange(
      desc(rise),
      desc(peak_d1),
      onset_rank
    )


  leading_candidates <- leading_candidates %>%
    arrange(
      desc(rise),
      desc(peak_d1),
      onset_rank
    )


  rem <- remainder_candidates[
    1L,
    ,
    drop = FALSE
  ]


  lead <- leading_candidates[
    1L,
    ,
    drop = FALSE
  ]


  remainder_boundary <- as.integer(
    round(
      rem$onset_rank[
        1L
      ]
    )
  )


  leading_edge_boundary <- as.integer(
    round(
      lead$onset_rank[
        1L
      ]
    )
  )


  if (
    remainder_boundary >=
    leading_edge_boundary
  ) {
    return(
      list(
        valid = FALSE,
        reason = "boundary_order_failure",
        curve = data.frame(
          rank = consensus_df$rank,
          x = x,
          smooth_divergence = smooth_y,
          derivative_1 = d1,
          derivative_2 = d2,
          stringsAsFactors = FALSE
        ),
        crossings = crossings,
        episodes = episodes,
        spline = spline_fit
      )
    )
  }


  list(
    valid = TRUE,
    reason = NA_character_,
    remainder_boundary = remainder_boundary,
    leading_edge_boundary = leading_edge_boundary,
    remainder_episode = rem,
    leading_episode = lead,
    curve = data.frame(
      rank = consensus_df$rank,
      x = x,
      smooth_divergence = smooth_y,
      derivative_1 = d1,
      derivative_2 = d2,
      stringsAsFactors = FALSE
    ),
    crossings = crossings,
    episodes = episodes,
    spline = spline_fit,
    spline_spar = spline_fit$spar,
    spline_lambda = spline_fit$lambda,
    spline_df = spline_fit$df
  )
}


# =============================================================================
# REGION ASSIGNMENT AND SUMMARIES
# =============================================================================

assign_regions <- function(
    rank,
    remainder_boundary,
    leading_edge_boundary) {

  ifelse(
    rank < remainder_boundary,
    "REMAINDER",
    ifelse(
      rank < leading_edge_boundary,
      "TRANSITION",
      "LEADING_EDGE"
    )
  )
}


summarize_group_regions <- function(
    feature_df,
    remainder_boundary,
    leading_edge_boundary) {

  x <- feature_df %>%
    mutate(
      region = assign_regions(
        rank,
        remainder_boundary,
        leading_edge_boundary
      )
    )


  region_levels <- c(
    "REMAINDER",
    "TRANSITION",
    "LEADING_EDGE"
  )


  x$region <- factor(
    x$region,
    levels = region_levels
  )


  x %>%
    group_by(
      group,
      region
    ) %>%
    summarise(
      n_features = n(),
      median_P_star = median(
        P_star,
        na.rm = TRUE
      ),
      median_E_star = median(
        E_star,
        na.rm = TRUE
      ),
      median_divergence = median(
        divergence,
        na.rm = TRUE
      ),
      median_NB2 = median(
        NB2,
        na.rm = TRUE
      ),
      median_NB2_NB1 = median(
        NB2_NB1,
        na.rm = TRUE
      ),
      median_alpha_mu = median(
        alpha_mu,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
}


# =============================================================================
# BOOTSTRAP OF THE COMPLETE SHARED-BOUNDARY PROCEDURE
# =============================================================================

bootstrap_shared_boundaries <- function(
    count_mat,
    rank_matrix_all,
    normalized_counts,
    group_labels,
    bootstrap_n = BOOTSTRAP_N,
    seed = BOOTSTRAP_SEED,
    cores = BOOTSTRAP_CORES) {

  groups <- levels(
    factor(group_labels)
  )


  group_indices <- setNames(
    lapply(
      groups,
      function(g) {
        which(
          group_labels == g
        )
      }
    ),
    groups
  )


  set.seed(
    seed
  )


  bootstrap_draws <- lapply(
    seq_len(bootstrap_n),
    function(b) {
      setNames(
        lapply(
          groups,
          function(g) {
            idx <- group_indices[[g]]

            sample(
              idx,
              size = length(idx),
              replace = TRUE
            )
          }
        ),
        groups
      )
    }
  )


  worker <- function(b) {

    draws <- bootstrap_draws[[
      b
    ]]


    invalid_row <- function(reason) {
      data.frame(
        bootstrap = b,
        valid = FALSE,
        reason = reason,
        remainder_boundary = NA_integer_,
        leading_edge_boundary = NA_integer_,
        transition_width = NA_integer_,
        leading_edge_size = NA_integer_,
        spline_spar = NA_real_,
        spline_df = NA_real_,
        stringsAsFactors = FALSE
      )
    }


    pooled_variance_b <- tryCatch(
      compute_pooled_variance_from_draws(
        normalized_counts = normalized_counts,
        draws_by_group = draws
      ),
      error = function(e) {
        NULL
      }
    )


    if (is.null(pooled_variance_b)) {
      return(
        invalid_row(
          "pooled_variance_failed"
        )
      )
    }


    curve_list <- vector(
      "list",
      length(groups)
    )

    names(curve_list) <- groups


    for (
      g in groups
    ) {

      idx <- draws[[g]]


      if (
        length(
          unique(idx)
        ) < 2L
      ) {
        return(
          invalid_row(
            paste0(
              "fewer_than_2_unique_samples_",
              g
            )
          )
        )
      }


      one <- tryCatch(
        compute_one_group_curve(
          group_name = g,
          rank_matrix_arm = rank_matrix_all[
            ,
            idx,
            drop = FALSE
          ],
          normalized_counts_arm = normalized_counts[
            ,
            idx,
            drop = FALSE
          ],
          pooled_variance = pooled_variance_b
        ),
        error = function(e) {
          NULL
        }
      )


      if (is.null(one)) {
        return(
          invalid_row(
            paste0(
              "group_curve_failed_",
              g
            )
          )
        )
      }


      curve_list[[g]] <- one$feature_df
    }


    consensus_b <- tryCatch(
      compute_consensus_curve(
        curve_list
      ),
      error = function(e) {
        NULL
      }
    )


    if (is.null(consensus_b)) {
      return(
        invalid_row(
          "consensus_failed"
        )
      )
    }


    detector_b <- tryCatch(
      detect_shared_boundaries(
        consensus_b
      ),
      error = function(e) {
        NULL
      }
    )


    if (
      is.null(detector_b) ||
      !isTRUE(
        detector_b$valid
      )
    ) {

      reason <- if (
        is.null(detector_b)
      ) {
        "detector_error"
      } else {
        detector_b$reason
      }

      return(
        invalid_row(
          reason
        )
      )
    }


    rem <- detector_b$remainder_boundary

    lead <- detector_b$leading_edge_boundary

    N <- nrow(
      consensus_b
    )


    data.frame(
      bootstrap = b,
      valid = TRUE,
      reason = NA_character_,
      remainder_boundary = rem,
      leading_edge_boundary = lead,
      transition_width = lead - rem,
      leading_edge_size = N - lead + 1L,
      spline_spar = detector_b$spline_spar,
      spline_df = detector_b$spline_df,
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
      mc.cores = cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
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


  if (nrow(valid_df) == 0L) {
    stop("No valid bootstrap boundary estimates were obtained.")
  }


  summary_df <- data.frame(
    BootstrapN = bootstrap_n,
    ValidBootstraps = nrow(valid_df),
    ValidRate = nrow(valid_df) / bootstrap_n,

    RemainderMedian = safe_quantile(
      valid_df$remainder_boundary,
      0.50
    ),
    RemainderQ25 = safe_quantile(
      valid_df$remainder_boundary,
      0.25
    ),
    RemainderQ75 = safe_quantile(
      valid_df$remainder_boundary,
      0.75
    ),

    LeadingMedian = safe_quantile(
      valid_df$leading_edge_boundary,
      0.50
    ),
    LeadingQ25 = safe_quantile(
      valid_df$leading_edge_boundary,
      0.25
    ),
    LeadingQ75 = safe_quantile(
      valid_df$leading_edge_boundary,
      0.75
    ),

    TransitionWidthMedian = safe_quantile(
      valid_df$transition_width,
      0.50
    ),
    TransitionWidthQ25 = safe_quantile(
      valid_df$transition_width,
      0.25
    ),
    TransitionWidthQ75 = safe_quantile(
      valid_df$transition_width,
      0.75
    ),

    LeadingEdgeSizeMedian = safe_quantile(
      valid_df$leading_edge_size,
      0.50
    ),
    LeadingEdgeSizeQ25 = safe_quantile(
      valid_df$leading_edge_size,
      0.25
    ),
    LeadingEdgeSizeQ75 = safe_quantile(
      valid_df$leading_edge_size,
      0.75
    ),

    MedianSplineSpar = safe_median(
      valid_df$spline_spar
    ),

    MedianSplineDF = safe_median(
      valid_df$spline_df
    ),

    stringsAsFactors = FALSE
  )


  list(
    bootstrap_df = bootstrap_df,
    summary_df = summary_df
  )
}


# =============================================================================
# FIGURE LAYOUT HELPERS
# =============================================================================

save_four_panel_plot <- function(
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
        nrow = 4L,
        ncol = 1L,
        heights = unit(
          c(
            1.05,
            1.05,
            0.95,
            0.95
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


save_three_panel_plot <- function(
    plot_list,
    filename) {

  png(
    filename,
    width = PNG_WIDTH_IN,
    height = 11.5,
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
            1.10,
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
      vp = viewport(
        layout.pos.row = i,
        layout.pos.col = 1L
      )
    )
  }


  dev.off()
}


add_region_backgrounds <- function(
    p,
    remainder_boundary,
    leading_edge_boundary,
    N) {

  p +
    annotate(
      "rect",
      xmin = 1,
      xmax = remainder_boundary,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$remainder,
      alpha = 0.55
    ) +
    annotate(
      "rect",
      xmin = remainder_boundary,
      xmax = leading_edge_boundary,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$transition,
      alpha = 0.42
    ) +
    annotate(
      "rect",
      xmin = leading_edge_boundary,
      xmax = N,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$leading,
      alpha = 0.55
    )
}


# =============================================================================
# CONSENSUS MANUSCRIPT FIGURE
# =============================================================================

build_consensus_figure <- function(
    consensus_df,
    detector,
    bootstrap,
    out_file) {

  N <- nrow(
    consensus_df
  )


  rem <- detector$remainder_boundary

  lead <- detector$leading_edge_boundary


  curve_df <- consensus_df %>%
    left_join(
      detector$curve %>%
        select(
          rank,
          smooth_divergence,
          derivative_1,
          derivative_2
        ),
      by = "rank"
    )


  component_long <- curve_df %>%
    select(
      rank,
      consensus_P_star,
      consensus_E_star,
      consensus_divergence
    ) %>%
    pivot_longer(
      cols = c(
        consensus_P_star,
        consensus_E_star,
        consensus_divergence
      ),
      names_to = "quantity",
      values_to = "value"
    ) %>%
    mutate(
      quantity = factor(
        quantity,
        levels = c(
          "consensus_P_star",
          "consensus_E_star",
          "consensus_divergence"
        ),
        labels = c(
          "PC1 variance contribution, P*",
          "NB excess variance, E*",
          "Divergence, D"
        )
      )
    )


  p1 <- ggplot(
    component_long,
    aes(
      rank,
      value,
      color = quantity
    )
  )


  p1 <- add_region_backgrounds(
    p1,
    rem,
    lead,
    N
  )


  p1 <- p1 +
    geom_line(
      linewidth = 0.75
    ) +
    geom_vline(
      xintercept = rem,
      color = COL$rem_line,
      linewidth = 0.8
    ) +
    geom_vline(
      xintercept = lead,
      color = COL$lead_line,
      linewidth = 0.8
    ) +
    scale_color_manual(
      values = c(
        "PC1 variance contribution, P*" = COL$pc1,
        "NB excess variance, E*" = COL$excess,
        "Divergence, D" = COL$divergence
      )
    ) +
    labs(
      title = "Shared PC1–negative-binomial variance divergence",
      subtitle = paste0(
        "Shared remainder boundary = ",
        rem,
        "; shared leading-edge boundary = ",
        lead
      ),
      x = "Rank (low -> high PC1 variance contribution)",
      y = "Robust standardized value",
      color = NULL
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )


  y_top_1 <- max(
    component_long$value,
    na.rm = TRUE
  )

  y_bottom_1 <- min(
    component_long$value,
    na.rm = TRUE
  )

  y_span_1 <- y_top_1 -
    y_bottom_1

  if (
    !is.finite(y_span_1) ||
    y_span_1 <= 0
  ) {
    y_span_1 <- 1
  }


  p1 <- p1 +
    annotate(
      "label",
      x = max(1, round(N * 0.03)),
      y = y_top_1 - 0.05 * y_span_1,
      label = "P[i] == d[1]^2*v[i*1]^2/(n-1)",
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      fill = "white"
    ) +
    annotate(
      "label",
      x = max(1, round(N * 0.03)),
      y = y_top_1 - 0.19 * y_span_1,
      label = "E[i*g] == max(V[i]-mu[i*g],0)",
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      fill = "white"
    ) +
    annotate(
      "label",
      x = max(1, round(N * 0.03)),
      y = y_top_1 - 0.33 * y_span_1,
      label = "D[g](r) == Z[R](log(1+E[g](r)))-Z[R](log(1+P[g](r)))",
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      fill = "white"
    )


  p2 <- ggplot(
    curve_df,
    aes(
      rank,
      smooth_divergence
    )
  )


  p2 <- add_region_backgrounds(
    p2,
    rem,
    lead,
    N
  )


  p2 <- p2 +
    geom_line(
      color = COL$divergence,
      linewidth = 1.0
    ) +
    geom_vline(
      xintercept = rem,
      color = COL$rem_line,
      linewidth = 0.9
    ) +
    geom_vline(
      xintercept = lead,
      color = COL$lead_line,
      linewidth = 0.9
    ) +
    labs(
      title = "Smoothed consensus divergence and shared boundaries",
      subtitle = paste0(
        "Spline selected by generalized cross-validation; effective df = ",
        round(
          detector$spline_df,
          1
        )
      ),
      x = "Rank",
      y = "Smoothed consensus divergence"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank()
    )


  y_top_2 <- max(
    curve_df$smooth_divergence,
    na.rm = TRUE
  )

  y_bottom_2 <- min(
    curve_df$smooth_divergence,
    na.rm = TRUE
  )

  y_span_2 <- y_top_2 -
    y_bottom_2

  if (
    !is.finite(y_span_2) ||
    y_span_2 <= 0
  ) {
    y_span_2 <- 1
  }


  p2 <- p2 +
    annotate(
      "label",
      x = max(1, round(N * 0.03)),
      y = y_top_2 - 0.05 * y_span_2,
      label = "D[cons](r) == median[g](D[g](r))",
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 3.1,
      fill = "white"
    ) +
    annotate(
      "text",
      x = rem,
      y = y_bottom_2 + 0.08 * y_span_2,
      label = "REMAINDER\nBOUNDARY",
      color = COL$rem_line,
      hjust = 0,
      size = 3.0
    ) +
    annotate(
      "text",
      x = lead,
      y = y_bottom_2 + 0.08 * y_span_2,
      label = "LEADING-EDGE\nBOUNDARY",
      color = COL$lead_line,
      hjust = 1,
      size = 3.0
    )


  derivative_df <- curve_df %>%
    transmute(
      rank = rank,
      `D'(r)` = robust_z(
        derivative_1
      ),
      `D''(r)` = robust_z(
        derivative_2
      )
    ) %>%
    pivot_longer(
      cols = c(
        `D'(r)`,
        `D''(r)`
      ),
      names_to = "derivative",
      values_to = "value"
    )


  p3 <- ggplot(
    derivative_df,
    aes(
      rank,
      value,
      color = derivative
    )
  )


  p3 <- add_region_backgrounds(
    p3,
    rem,
    lead,
    N
  )


  p3 <- p3 +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_line(
      linewidth = 0.85
    ) +
    geom_vline(
      xintercept = rem,
      color = COL$rem_line,
      linewidth = 0.8
    ) +
    geom_vline(
      xintercept = lead,
      color = COL$lead_line,
      linewidth = 0.8
    ) +
    scale_color_manual(
      values = c(
        "D'(r)" = COL$derivative1,
        "D''(r)" = COL$derivative2
      )
    ) +
    labs(
      title = "Calculus of the consensus divergence trajectory",
      subtitle = "Boundaries mark the strongest positive-divergence acceleration episodes on the remainder and leading-edge sides",
      x = "Rank",
      y = "Robust standardized derivative",
      color = NULL
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )


  y_top_3 <- max(
    derivative_df$value,
    na.rm = TRUE
  )

  y_bottom_3 <- min(
    derivative_df$value,
    na.rm = TRUE
  )

  y_span_3 <- y_top_3 -
    y_bottom_3

  if (
    !is.finite(y_span_3) ||
    y_span_3 <= 0
  ) {
    y_span_3 <- 1
  }


  p3 <- p3 +
    annotate(
      "label",
      x = max(1, round(N * 0.03)),
      y = y_top_3 - 0.05 * y_span_3,
      label = "D''_cons(r): negative -> positive",
      parse = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      fill = "white"
    )


  valid_boot <- bootstrap$bootstrap_df %>%
    filter(valid)


  boot_long <- bind_rows(
    data.frame(
      boundary = "Remainder boundary",
      rank = valid_boot$remainder_boundary,
      stringsAsFactors = FALSE
    ),
    data.frame(
      boundary = "Leading-edge boundary",
      rank = valid_boot$leading_edge_boundary,
      stringsAsFactors = FALSE
    )
  )


  q_df <- data.frame(
    boundary = c(
      "Remainder boundary",
      "Leading-edge boundary"
    ),
    full = c(
      rem,
      lead
    ),
    q25 = c(
      bootstrap$summary_df$RemainderQ25[
        1L
      ],
      bootstrap$summary_df$LeadingQ25[
        1L
      ]
    ),
    median = c(
      bootstrap$summary_df$RemainderMedian[
        1L
      ],
      bootstrap$summary_df$LeadingMedian[
        1L
      ]
    ),
    q75 = c(
      bootstrap$summary_df$RemainderQ75[
        1L
      ],
      bootstrap$summary_df$LeadingQ75[
        1L
      ]
    ),
    stringsAsFactors = FALSE
  )


  p4 <- ggplot(
    boot_long,
    aes(
      rank
    )
  ) +
    geom_histogram(
      bins = 40,
      alpha = 0.70
    ) +
    geom_vline(
      data = q_df,
      aes(
        xintercept = full
      ),
      linewidth = 0.8
    ) +
    geom_vline(
      data = q_df,
      aes(
        xintercept = median
      ),
      linetype = "dashed",
      linewidth = 0.8
    ) +
    geom_segment(
      data = q_df,
      aes(
        x = q25,
        xend = q75,
        y = Inf,
        yend = Inf
      ),
      inherit.aes = FALSE,
      linewidth = 2.0
    ) +
    facet_wrap(
      ~ boundary,
      nrow = 1,
      scales = "free_y"
    ) +
    labs(
      title = "Bootstrap-supported shared boundary ranges",
      subtitle = paste0(
        "Valid bootstrap replicates = ",
        bootstrap$summary_df$ValidBootstraps[
          1L
        ],
        "/",
        bootstrap$summary_df$BootstrapN[
          1L
        ],
        "; solid = full-data estimate; dashed = bootstrap median; top bar = IQR"
      ),
      x = "Rank",
      y = "Bootstrap count"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank()
    )


  save_four_panel_plot(
    list(
      p1,
      p2,
      p3,
      p4
    ),
    out_file
  )
}


# =============================================================================
# ARM-LEVEL MANUSCRIPT FIGURE
# =============================================================================

build_group_figure <- function(
    feature_df,
    region_summary,
    remainder_boundary,
    leading_edge_boundary,
    out_file,
    title_text) {

  N <- nrow(
    feature_df
  )


  transformed_long <- feature_df %>%
    select(
      rank,
      P_star,
      E_star,
      divergence
    ) %>%
    pivot_longer(
      cols = c(
        P_star,
        E_star,
        divergence
      ),
      names_to = "quantity",
      values_to = "value"
    ) %>%
    mutate(
      quantity = factor(
        quantity,
        levels = c(
          "P_star",
          "E_star",
          "divergence"
        ),
        labels = c(
          "P*",
          "E*",
          "D = E* - P*"
        )
      )
    )


  p1 <- ggplot(
    transformed_long,
    aes(
      rank,
      value,
      color = quantity
    )
  )


  p1 <- add_region_backgrounds(
    p1,
    remainder_boundary,
    leading_edge_boundary,
    N
  )


  p1 <- p1 +
    geom_line(
      linewidth = 0.72
    ) +
    geom_vline(
      xintercept = remainder_boundary,
      color = COL$rem_line,
      linewidth = 0.8
    ) +
    geom_vline(
      xintercept = leading_edge_boundary,
      color = COL$lead_line,
      linewidth = 0.8
    ) +
    scale_color_manual(
      values = c(
        "P*" = COL$pc1,
        "E*" = COL$excess,
        "D = E* - P*" = COL$divergence
      )
    ) +
    labs(
      title = paste0(
        title_text,
        ": PC1–NB divergence"
      ),
      subtitle = paste0(
        "Shared remainder boundary = ",
        remainder_boundary,
        "; shared leading-edge boundary = ",
        leading_edge_boundary
      ),
      x = "Rank (low -> high PC1 variance contribution)",
      y = "Robust standardized value",
      color = NULL
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )


  y_top <- max(
    transformed_long$value,
    na.rm = TRUE
  )

  y_bottom <- min(
    transformed_long$value,
    na.rm = TRUE
  )

  y_span <- y_top -
    y_bottom

  if (
    !is.finite(y_span) ||
    y_span <= 0
  ) {
    y_span <- 1
  }


  p1 <- p1 +
    annotate(
      "label",
      x = max(1, round(N * 0.03)),
      y = y_top - 0.05 * y_span,
      label = "P[i] == d[1]^2*v[i*1]^2/(n-1)",
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      fill = "white"
    ) +
    annotate(
      "label",
      x = max(1, round(N * 0.03)),
      y = y_top - 0.19 * y_span,
      label = "E[i*g] == max(V[i]-mu[i*g],0)",
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      fill = "white"
    ) +
    annotate(
      "label",
      x = max(1, round(N * 0.03)),
      y = y_top - 0.33 * y_span,
      label = "D_g(r) = E*_g(r) - P*_g(r)",
      parse = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      fill = "white"
    )


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


  p2 <- ggplot(
    nb_long,
    aes(
      rank,
      value,
      color = metric
    )
  )


  p2 <- add_region_backgrounds(
    p2,
    remainder_boundary,
    leading_edge_boundary,
    N
  )


  p2 <- p2 +
    geom_line(
      linewidth = 0.78
    ) +
    geom_vline(
      xintercept = remainder_boundary,
      color = COL$rem_line,
      linewidth = 0.8
    ) +
    geom_vline(
      xintercept = leading_edge_boundary,
      color = COL$lead_line,
      linewidth = 0.8
    ) +
    scale_color_manual(
      values = c(
        "NB2" = COL$nb2,
        "NB2-NB1" = COL$nb_gap,
        "alpha*mu" = COL$alpha_mu
      )
    ) +
    labs(
      title = paste0(
        title_text,
        ": regional negative-binomial variance structure"
      ),
      x = "Rank",
      y = "Moment-based NB signal",
      color = NULL
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )


  y_top_2 <- max(
    nb_long$value,
    na.rm = TRUE
  )

  y_bottom_2 <- min(
    nb_long$value,
    na.rm = TRUE
  )

  y_span_2 <- y_top_2 -
    y_bottom_2

  if (
    !is.finite(y_span_2) ||
    y_span_2 <= 0
  ) {
    y_span_2 <- 1
  }


  p2 <- p2 +
    annotate(
      "label",
      x = max(1, round(N * 0.03)),
      y = y_top_2 - 0.05 * y_span_2,
      label = "NB2[i] == log(1+E[i])",
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      fill = "white"
    ) +
    annotate(
      "label",
      x = max(1, round(N * 0.03)),
      y = y_top_2 - 0.19 * y_span_2,
      label = "NB2-NB1_i = log(1+E_i) - log(1+mu_i)",
      parse = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.9,
      fill = "white"
    ) +
    annotate(
      "label",
      x = max(1, round(N * 0.03)),
      y = y_top_2 - 0.33 * y_span_2,
      label = "alpha[i] == E[i]/mu[i]^2",
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      fill = "white"
    )


  summary_long <- region_summary %>%
    select(
      region,
      median_NB2,
      median_NB2_NB1,
      median_alpha_mu
    ) %>%
    pivot_longer(
      cols = c(
        median_NB2,
        median_NB2_NB1,
        median_alpha_mu
      ),
      names_to = "metric",
      values_to = "median_value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c(
          "median_NB2",
          "median_NB2_NB1",
          "median_alpha_mu"
        ),
        labels = c(
          "NB2",
          "NB2-NB1",
          "alpha*mu"
        )
      )
    )


  p3 <- ggplot(
    summary_long,
    aes(
      x = median_value,
      y = metric,
      shape = region
    )
  ) +
    geom_point(
      size = 3.4
    ) +
    facet_wrap(
      ~ region,
      nrow = 1
    ) +
    labs(
      title = paste0(
        title_text,
        ": remainder / transition / leading-edge summary"
      ),
      subtitle = "Shared boundaries are applied at the same rank positions in every experimental arm",
      x = "Regional median",
      y = NULL,
      shape = NULL
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none"
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
# RUN ANALYSIS
# =============================================================================

count_mat <- read_count_matrix(
  COUNT_FILE,
  GROUP_PATTERNS
)


group_labels <- assign_sample_groups(
  sample_names = colnames(count_mat),
  group_patterns = GROUP_PATTERNS
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


# -----------------------------------------------------------------------------
# Normalize sequencing counts for pooled design-aware variance.
# -----------------------------------------------------------------------------

deseq2_norm <- normalize_deseq2(
  count_mat = count_mat,
  group_labels = group_labels
)


normalized_counts <- deseq2_norm$normalized_counts


# -----------------------------------------------------------------------------
# PC1 ranking matrix.
#
# CPM-log1p preserves the main manuscript ranking convention while the SVD
# contribution P_i explicitly quantifies how much variance each feature
# contributes through PC1.
# -----------------------------------------------------------------------------

rank_matrix_all <- normalize_cpm_log1p(
  count_mat
)


# -----------------------------------------------------------------------------
# Full-data pooled within-group sequencing variance.
# -----------------------------------------------------------------------------

pooled_obj <- compute_pooled_within_group_variance(
  normalized_counts = normalized_counts,
  group_labels = group_labels
)


pooled_variance <- pooled_obj$variance


message(
  "Pooled within-group residual degrees of freedom: ",
  pooled_obj$residual_df
)


# -----------------------------------------------------------------------------
# Full-data arm-specific divergence curves.
# -----------------------------------------------------------------------------

all_curves <- compute_all_group_curves(
  rank_matrix_all = rank_matrix_all,
  normalized_counts = normalized_counts,
  pooled_variance = pooled_variance,
  group_labels = group_labels
)


group_curves <- all_curves$curves


write.csv(
  all_curves$pc1_summary,
  file.path(
    OUT_ROOT,
    "Table_PC1_Summary.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# Shared rank-wise consensus divergence.
# -----------------------------------------------------------------------------

consensus_df <- compute_consensus_curve(
  group_curves
)


detector <- detect_shared_boundaries(
  consensus_df
)


if (
  !isTRUE(
    detector$valid
  )
) {
  stop(
    "Full-data shared boundary detector failed: ",
    detector$reason
  )
}


REMAINDER_BOUNDARY <- detector$remainder_boundary

LEADING_EDGE_BOUNDARY <- detector$leading_edge_boundary


message(
  "Shared REMAINDER_BOUNDARY = ",
  REMAINDER_BOUNDARY
)

message(
  "Shared LEADING_EDGE_BOUNDARY = ",
  LEADING_EDGE_BOUNDARY
)

message(
  "Shared transition interval width = ",
  LEADING_EDGE_BOUNDARY -
  REMAINDER_BOUNDARY
)

message(
  "Shared terminal leading-edge size = ",
  nrow(consensus_df) -
  LEADING_EDGE_BOUNDARY +
  1L
)


# -----------------------------------------------------------------------------
# Bootstrap the entire shared-boundary procedure.
# -----------------------------------------------------------------------------

message(
  "Bootstrapping ",
  BOOTSTRAP_N,
  " complete shared-boundary replicates on ",
  BOOTSTRAP_CORES,
  " core(s)..."
)


bootstrap <- bootstrap_shared_boundaries(
  count_mat = count_mat,
  rank_matrix_all = rank_matrix_all,
  normalized_counts = normalized_counts,
  group_labels = group_labels,
  bootstrap_n = BOOTSTRAP_N,
  seed = BOOTSTRAP_SEED,
  cores = BOOTSTRAP_CORES
)


write.csv(
  bootstrap$bootstrap_df,
  file.path(
    OUT_ROOT,
    "Table_Bootstrap_Boundaries.csv"
  ),
  row.names = FALSE
)


boundary_summary <- data.frame(
  RemainderBoundaryFull = REMAINDER_BOUNDARY,
  RemainderBoundaryBootstrapMedian =
    bootstrap$summary_df$RemainderMedian[1L],
  RemainderBoundaryBootstrapQ25 =
    bootstrap$summary_df$RemainderQ25[1L],
  RemainderBoundaryBootstrapQ75 =
    bootstrap$summary_df$RemainderQ75[1L],

  LeadingEdgeBoundaryFull = LEADING_EDGE_BOUNDARY,
  LeadingEdgeBoundaryBootstrapMedian =
    bootstrap$summary_df$LeadingMedian[1L],
  LeadingEdgeBoundaryBootstrapQ25 =
    bootstrap$summary_df$LeadingQ25[1L],
  LeadingEdgeBoundaryBootstrapQ75 =
    bootstrap$summary_df$LeadingQ75[1L],

  TransitionWidthFull =
    LEADING_EDGE_BOUNDARY -
    REMAINDER_BOUNDARY,
  TransitionWidthBootstrapMedian =
    bootstrap$summary_df$TransitionWidthMedian[1L],
  TransitionWidthBootstrapQ25 =
    bootstrap$summary_df$TransitionWidthQ25[1L],
  TransitionWidthBootstrapQ75 =
    bootstrap$summary_df$TransitionWidthQ75[1L],

  LeadingEdgeSizeFull =
    nrow(consensus_df) -
    LEADING_EDGE_BOUNDARY +
    1L,
  LeadingEdgeSizeBootstrapMedian =
    bootstrap$summary_df$LeadingEdgeSizeMedian[1L],
  LeadingEdgeSizeBootstrapQ25 =
    bootstrap$summary_df$LeadingEdgeSizeQ25[1L],
  LeadingEdgeSizeBootstrapQ75 =
    bootstrap$summary_df$LeadingEdgeSizeQ75[1L],

  BootstrapN =
    bootstrap$summary_df$BootstrapN[1L],
  ValidBootstraps =
    bootstrap$summary_df$ValidBootstraps[1L],
  BootstrapValidRate =
    bootstrap$summary_df$ValidRate[1L],

  SplineSparFull =
    detector$spline_spar,
  SplineEffectiveDFFull =
    detector$spline_df,

  stringsAsFactors = FALSE
)


write.csv(
  boundary_summary,
  file.path(
    OUT_ROOT,
    "Table_Shared_Boundaries.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# Save the exact consensus curve and calculus quantities used in the figures.
# -----------------------------------------------------------------------------

consensus_output <- consensus_df %>%
  left_join(
    detector$curve %>%
      select(
        rank,
        smooth_divergence,
        derivative_1,
        derivative_2
      ),
    by = "rank"
  ) %>%
  mutate(
    region = assign_regions(
      rank,
      REMAINDER_BOUNDARY,
      LEADING_EDGE_BOUNDARY
    )
  )


write.csv(
  consensus_output,
  file.path(
    OUT_ROOT,
    "Table_Consensus_Divergence_Curve.csv"
  ),
  row.names = FALSE
)


write.csv(
  detector$crossings,
  file.path(
    OUT_ROOT,
    "Table_Consensus_Curvature_Crossings.csv"
  ),
  row.names = FALSE
)


write.csv(
  detector$episodes,
  file.path(
    OUT_ROOT,
    "Table_Consensus_Acceleration_Episodes.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# Global sample / normalization table.
# -----------------------------------------------------------------------------

write.csv(
  data.frame(
    sample = colnames(count_mat),
    group = as.character(group_labels),
    DESeq2_size_factor = as.numeric(
      deseq2_norm$size_factors
    ),
    stringsAsFactors = FALSE
  ),
  file.path(
    OUT_ROOT,
    "Table_Samples_Normalization.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# Consensus mathematical figure.
# -----------------------------------------------------------------------------

build_consensus_figure(
  consensus_df = consensus_df,
  detector = detector,
  bootstrap = bootstrap,
  out_file = file.path(
    OUT_ROOT,
    "Figure_Consensus_PC1_NB_Divergence.png"
  )
)


# -----------------------------------------------------------------------------
# Apply the SAME shared rank boundaries to every experimental arm.
# Save ranked features, region summaries, and manuscript figures.
# -----------------------------------------------------------------------------

region_rows <- list()


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


  arm_map <- COMPARISONS[[
    comparison_name
  ]]


  for (
    arm_name in names(
      arm_map
    )
  ) {

    group_name <- arm_map[[
      arm_name
    ]]


    feature_df <- group_curves[[
      group_name
    ]] %>%
      mutate(
        region = assign_regions(
          rank,
          REMAINDER_BOUNDARY,
          LEADING_EDGE_BOUNDARY
        )
      )


    region_summary <- summarize_group_regions(
      feature_df = feature_df,
      remainder_boundary = REMAINDER_BOUNDARY,
      leading_edge_boundary = LEADING_EDGE_BOUNDARY
    )


    region_summary <- region_summary %>%
      mutate(
        comparison = comparison_name,
        arm = arm_name,
        .before = 1L
      )


    region_rows[[
      length(region_rows) +
      1L
    ]] <- region_summary


    write.csv(
      feature_df,
      file.path(
        comp_dir,
        paste0(
          "Table_Ranked_",
          comparison_name,
          "_",
          arm_name,
          ".csv"
        )
      ),
      row.names = FALSE
    )


    write.csv(
      region_summary,
      file.path(
        comp_dir,
        paste0(
          "Table_Regions_",
          comparison_name,
          "_",
          arm_name,
          ".csv"
        )
      ),
      row.names = FALSE
    )


    build_group_figure(
      feature_df = feature_df,
      region_summary = region_summary,
      remainder_boundary = REMAINDER_BOUNDARY,
      leading_edge_boundary = LEADING_EDGE_BOUNDARY,
      out_file = file.path(
        comp_dir,
        paste0(
          "Figure_Main_",
          comparison_name,
          "_",
          arm_name,
          ".png"
        )
      ),
      title_text = paste0(
        comparison_name,
        " ",
        arm_name,
        " (",
        group_name,
        ")"
      )
    )
  }
}


overall_regions <- bind_rows(
  region_rows
)


write.csv(
  overall_regions,
  file.path(
    OUT_ROOT,
    "Table_Overall_Regions.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

message(
  "============================================================"
)

message(
  "Shared PC1-NB divergence analysis complete."
)

message(
  "REMAINDER_BOUNDARY: full = ",
  REMAINDER_BOUNDARY,
  "; bootstrap median = ",
  round(
    boundary_summary$RemainderBoundaryBootstrapMedian[
      1L
    ]
  ),
  "; IQR = ",
  round(
    boundary_summary$RemainderBoundaryBootstrapQ25[
      1L
    ]
  ),
  "-",
  round(
    boundary_summary$RemainderBoundaryBootstrapQ75[
      1L
    ]
  )
)

message(
  "LEADING_EDGE_BOUNDARY: full = ",
  LEADING_EDGE_BOUNDARY,
  "; bootstrap median = ",
  round(
    boundary_summary$LeadingEdgeBoundaryBootstrapMedian[
      1L
    ]
  ),
  "; IQR = ",
  round(
    boundary_summary$LeadingEdgeBoundaryBootstrapQ25[
      1L
    ]
  ),
  "-",
  round(
    boundary_summary$LeadingEdgeBoundaryBootstrapQ75[
      1L
    ]
  )
)

message(
  "Leading-edge size: full = ",
  boundary_summary$LeadingEdgeSizeFull[
    1L
  ],
  "; bootstrap median = ",
  round(
    boundary_summary$LeadingEdgeSizeBootstrapMedian[
      1L
    ]
  ),
  "; IQR = ",
  round(
    boundary_summary$LeadingEdgeSizeBootstrapQ25[
      1L
    ]
  ),
  "-",
  round(
    boundary_summary$LeadingEdgeSizeBootstrapQ75[
      1L
    ]
  )
)

message(
  "Bootstrap valid rate = ",
  round(
    boundary_summary$BootstrapValidRate[
      1L
    ],
    3
  )
)

message(
  "Outputs written to: ",
  OUT_ROOT
)

message(
  "============================================================"
)
