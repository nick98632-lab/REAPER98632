#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# PC1–NB DIFFERENTIAL GEOMETRY DIAGNOSTIC
# =============================================================================
#
# PURPOSE
# -------
# This script does NOT choose a cutoff and does NOT bootstrap.
#
# It preserves the original feature ordering:
#
#     ascending |PC1 loading|
#
# and calculates two differential quantities along that ranked axis.
#
#
# 1) PC1 VARIANCE CONTRIBUTION
# ----------------------------
# For the centered sample-by-feature matrix X = U D V^T,
#
#                 d1^2 * v_i1^2
#     P_i =       ------------- ,
#                     n - 1
#
# where d1 is the first singular value and v_i1 is feature i's PC1 loading.
#
# Ranking is performed explicitly by ascending |v_i1|.
#
#
# 2) DESIGN-AWARE SEQUENCING VARIANCE
# -----------------------------------
# DESeq2-normalized counts from all eight experimental groups are used to
# estimate pooled within-group variance:
#
#                 sum_g sum_{j in g} (y_ij - ybar_ig)^2
#     V_i =       --------------------------------------- .
#                          sum_g (n_g - 1)
#
# For feature i in group g:
#
#     mu_ig = mean normalized count
#
#     E_ig = max(V_i - mu_ig, 0)
#
# E_ig is the empirical excess-over-Poisson variance signal.
#
#
# 3) PC1–NB ELASTICITY
# --------------------
# Instead of subtracting standardized P and E values, calculate their local
# log-slope along normalized PC1 rank x:
#
#                  d log(E_g) / dx
#     eta_EP,g =  ------------------ .
#                  d log(P_g) / dx
#
# eta_EP is dimensionless.  Where defined:
#
#     eta_EP = 1   -> E and P change proportionally on the log scale
#     eta_EP > 1   -> NB excess variance changes faster than PC1 contribution
#     eta_EP < 1   -> PC1 contribution changes faster than NB excess variance
#
#
# 4) LOCAL NEGATIVE-BINOMIAL MEAN–VARIANCE EXPONENT
# -------------------------------------------------
# Write the excess variance law locally as
#
#     E = alpha * mu^p .
#
# Taking logs gives
#
#     log(E) = log(alpha) + p log(mu).
#
# Therefore the local effective exponent along the PC1-ranked trajectory is
#
#                d log(E_g) / dx
#     p_NB,g =  ------------------- .
#                d log(mu_g) / dx
#
# When alpha changes slowly over the local rank neighborhood:
#
#     p_NB approximately 1  -> NB1-like mean–variance scaling
#     p_NB approximately 2  -> NB2-like mean–variance scaling
#
# Because a derivative ratio is undefined when its denominator is numerically
# zero, eta_EP and p_NB are reported as NA at those positions.  The numerical
# tolerance is tied only to machine precision and derivative magnitude; it is
# not a biological cutoff.
#
#
# 5) SMOOTHING
# ------------
# log(P), log(E), and log(mu) are each fit as functions of normalized rank
#
#     x = (rank - 1)/(N - 1)
#
# with generalized-cross-validation-selected smoothing splines.  Derivatives
# are taken analytically from those spline fits.
#
# Only strictly positive raw values are used for the exact logarithms.
#
#
# 6) CONSENSUS
# ------------
# eta_EP and p_NB are calculated independently in each of the eight
# experimental groups and then summarized rank-wise across groups by the
# median and interquartile range.
#
# The figures contain the defining equations directly on the panels.
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <-
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <-
  "/root/REAPER98632/exports/PC1_NB_DIFFERENTIAL_GEOMETRY"

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
    stop("Count file is empty or malformed.")
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

  if (length(sample_idx) == 0L) {
    stop("No sample columns matched GROUP_PATTERNS.")
  }

  if (1L %in% sample_idx) {
    stop("Column 1 matched a sample pattern; column 1 must contain feature IDs.")
  }

  count_df <- raw_df[
    ,
    sample_idx,
    drop = FALSE
  ]

  count_mat <- do.call(
    cbind,
    lapply(
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
  )

  colnames(count_mat) <- colnames(
    count_df
  )

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

  lib_size <- colSums(
    count_mat,
    na.rm = TRUE
  )

  lib_size[
    !is.finite(lib_size) |
    lib_size <= 0
  ] <- 1

  cpm <- sweep(
    count_mat,
    2L,
    lib_size / 1e6,
    "/"
  )

  log1p(
    cpm
  )
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
      "DESeq2 is required. Install DESeq2 before running this script."
    )
  }

  col_data <- data.frame(
    group = factor(
      group_labels
    ),
    row.names = colnames(
      count_mat
    )
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(
      count_mat
    ),
    colData = col_data,
    design = ~ group
  )

  dds <- tryCatch(
    DESeq2::estimateSizeFactors(
      dds
    ),
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

  list(
    normalized_counts = DESeq2::counts(
      dds,
      normalized = TRUE
    ),
    size_factors = DESeq2::sizeFactors(
      dds
    )
  )
}


# =============================================================================
# POOLED WITHIN-GROUP VARIANCE
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
    nrow(
      normalized_counts
    )
  )

  residual_df <- 0L

  for (
    g in groups
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

    residuals_g <- sweep(
      xg,
      1L,
      mu_g,
      "-"
    )

    sse <- sse +
      rowSums(
        residuals_g^2
      )

    residual_df <- residual_df +
      length(idx) -
      1L
  }

  if (residual_df < 2L) {
    stop("Pooled residual degrees of freedom < 2.")
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


# =============================================================================
# PC1 / SVD
# =============================================================================

compute_pc1 <- function(
    rank_matrix_arm) {

  if (
    nrow(rank_matrix_arm) < 2L ||
    ncol(rank_matrix_arm) < 2L
  ) {
    stop("PC1 requires at least two features and two samples.")
  }

  # samples x features
  X <- t(
    rank_matrix_arm
  )

  Xc <- sweep(
    X,
    2L,
    colMeans(
      X
    ),
    "-"
  )

  # Efficient sample-space SVD geometry.
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

  # Right singular vector / feature loading:
  #     v1 = X^T u1 / d1
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

  abs_v1 <- abs(
    v1
  )

  n_samples <- nrow(
    Xc
  )

  P <- (
    d1^2 *
    v1^2
  ) /
    (
      n_samples -
      1L
    )

  total_ss <- sum(
    Xc^2
  )

  pc1_fraction_total <- if (
    is.finite(total_ss) &&
    total_ss > 0
  ) {
    d1^2 / total_ss
  } else {
    NA_real_
  }

  list(
    loading = v1,
    abs_loading = abs_v1,
    contribution = P,
    singular_value = d1,
    eigenvalue = d1^2,
    pc1_fraction_total = pc1_fraction_total
  )
}


# =============================================================================
# SPLINE / DIFFERENTIAL HELPERS
# =============================================================================

fit_positive_log_spline <- function(
    x,
    quantity,
    label) {

  keep <- (
    is.finite(x) &
    is.finite(quantity) &
    quantity > 0
  )

  if (sum(keep) < 10L) {
    stop(
      "Too few positive observations to fit log spline for ",
      label,
      "."
    )
  }

  x_fit <- x[
    keep
  ]

  y_fit <- log(
    quantity[
      keep
    ]
  )

  fit <- stats::smooth.spline(
    x = x_fit,
    y = y_fit,
    cv = FALSE
  )

  domain <- range(
    x_fit,
    finite = TRUE
  )

  pred0 <- rep(
    NA_real_,
    length(x)
  )

  pred1 <- rep(
    NA_real_,
    length(x)
  )

  inside <- (
    x >= domain[
      1L
    ] &
    x <= domain[
      2L
    ]
  )

  pred0[
    inside
  ] <- as.numeric(
    stats::predict(
      fit,
      x = x[
        inside
      ],
      deriv = 0
    )$y
  )

  pred1[
    inside
  ] <- as.numeric(
    stats::predict(
      fit,
      x = x[
        inside
      ],
      deriv = 1
    )$y
  )

  list(
    fit = fit,
    log_quantity = pred0,
    derivative = pred1,
    positive = quantity > 0,
    domain = domain
  )
}


safe_derivative_ratio <- function(
    numerator,
    denominator) {

  out <- rep(
    NA_real_,
    length(numerator)
  )

  finite_den <- denominator[
    is.finite(
      denominator
    )
  ]

  if (length(finite_den) == 0L) {
    return(
      list(
        ratio = out,
        tolerance = NA_real_
      )
    )
  }

  den_scale <- max(
    abs(
      finite_den
    ),
    na.rm = TRUE
  )

  if (
    !is.finite(den_scale) ||
    den_scale <= 0
  ) {
    return(
      list(
        ratio = out,
        tolerance = NA_real_
      )
    )
  }

  # Purely numerical guard: sqrt(machine precision) relative to the observed
  # derivative scale.  This prevents division by a derivative indistinguishable
  # from zero without imposing a biological threshold.
  tolerance <- sqrt(
    .Machine$double.eps
  ) *
    den_scale

  valid <- (
    is.finite(numerator) &
    is.finite(denominator) &
    abs(denominator) >
      tolerance
  )

  out[
    valid
  ] <- numerator[
    valid
  ] /
    denominator[
      valid
    ]

  list(
    ratio = out,
    tolerance = tolerance
  )
}


rankwise_median <- function(
    mat) {

  apply(
    mat,
    1L,
    function(z) {

      z <- z[
        is.finite(z)
      ]

      if (length(z) == 0L) {
        return(
          NA_real_
        )
      }

      median(
        z
      )
    }
  )
}


rankwise_quantile <- function(
    mat,
    p) {

  apply(
    mat,
    1L,
    function(z) {

      z <- z[
        is.finite(z)
      ]

      if (length(z) == 0L) {
        return(
          NA_real_
        )
      }

      as.numeric(
        stats::quantile(
          z,
          probs = p,
          names = FALSE,
          type = 7
        )
      )
    }
  )
}


rankwise_n_finite <- function(
    mat) {

  rowSums(
    is.finite(
      mat
    )
  )
}


# =============================================================================
# ONE GROUP
# =============================================================================

compute_group_geometry <- function(
    group_name,
    rank_matrix_arm,
    normalized_counts_arm,
    pooled_variance) {

  pc1 <- compute_pc1(
    rank_matrix_arm
  )

  # ---------------------------------------------------------------------------
  # EXACT RANKING RULE:
  #
  #     ascending absolute PC1 loading
  #
  # Sign is retained in the output, but magnitude determines rank.
  # ---------------------------------------------------------------------------
  rank_order <- order(
    pc1$abs_loading,
    decreasing = FALSE
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

  E <- pmax(
    pooled_variance -
    mu,
    0
  )

  P_ranked <- pc1$contribution[
    rank_order
  ]

  E_ranked <- E[
    rank_order
  ]

  mu_ranked <- mu[
    rank_order
  ]

  V_ranked <- pooled_variance[
    rank_order
  ]

  N <- length(
    rank_order
  )

  rank <- seq_len(
    N
  )

  x <- (
    rank -
    1
  ) /
    (
      N -
      1
    )

  # Exact logarithms are fit only on positive values.
  fit_P <- fit_positive_log_spline(
    x = x,
    quantity = P_ranked,
    label = paste0(
      group_name,
      " P"
    )
  )

  fit_E <- fit_positive_log_spline(
    x = x,
    quantity = E_ranked,
    label = paste0(
      group_name,
      " E"
    )
  )

  fit_mu <- fit_positive_log_spline(
    x = x,
    quantity = mu_ranked,
    label = paste0(
      group_name,
      " mu"
    )
  )

  # ---------------------------------------------------------------------------
  # Differential field 1:
  #
  #                 d log(E)/dx
  #     eta_EP =   ---------------
  #                 d log(P)/dx
  # ---------------------------------------------------------------------------
  eta_obj <- safe_derivative_ratio(
    numerator = fit_E$derivative,
    denominator = fit_P$derivative
  )

  eta_EP <- eta_obj$ratio

  eta_valid_raw <- (
    P_ranked > 0 &
    E_ranked > 0
  )

  eta_EP[
    !eta_valid_raw
  ] <- NA_real_

  # ---------------------------------------------------------------------------
  # Differential field 2:
  #
  #                d log(E)/dx
  #     p_NB =    ---------------
  #                d log(mu)/dx
  #
  # This is the local effective exponent in E = alpha * mu^p.
  # ---------------------------------------------------------------------------
  p_obj <- safe_derivative_ratio(
    numerator = fit_E$derivative,
    denominator = fit_mu$derivative
  )

  p_NB <- p_obj$ratio

  p_valid_raw <- (
    E_ranked > 0 &
    mu_ranked > 0
  )

  p_NB[
    !p_valid_raw
  ] <- NA_real_

  feature_df <- data.frame(
    group = group_name,
    rank = rank,
    normalized_rank = x,

    feature_id = rownames(
      rank_matrix_arm
    )[
      rank_order
    ],

    pc1_loading = pc1$loading[
      rank_order
    ],

    abs_pc1_loading = pc1$abs_loading[
      rank_order
    ],

    pc1_variance_contribution = P_ranked,

    group_mean_normalized = mu_ranked,

    pooled_within_group_variance = V_ranked,

    nb_excess_variance = E_ranked,

    smooth_log_P = fit_P$log_quantity,

    smooth_log_E = fit_E$log_quantity,

    smooth_log_mu = fit_mu$log_quantity,

    d_logP_dx = fit_P$derivative,

    d_logE_dx = fit_E$derivative,

    d_logmu_dx = fit_mu$derivative,

    eta_EP = eta_EP,

    p_NB = p_NB,

    stringsAsFactors = FALSE
  )

  summary_df <- data.frame(
    group = group_name,

    n_samples = ncol(
      rank_matrix_arm
    ),

    singular_value_1 = pc1$singular_value,

    eigenvalue_1 = pc1$eigenvalue,

    pc1_fraction_total_variance = pc1$pc1_fraction_total,

    spline_df_logP = fit_P$fit$df,

    spline_df_logE = fit_E$fit$df,

    spline_df_logmu = fit_mu$fit$df,

    eta_denominator_tolerance =
      eta_obj$tolerance,

    p_denominator_tolerance =
      p_obj$tolerance,

    eta_valid_rank_fraction =
      mean(
        is.finite(
          eta_EP
        )
      ),

    p_valid_rank_fraction =
      mean(
        is.finite(
          p_NB
        )
      ),

    stringsAsFactors = FALSE
  )

  list(
    features = feature_df,
    summary = summary_df
  )
}


# =============================================================================
# CONSENSUS ACROSS ALL EIGHT GROUPS
# =============================================================================

compute_consensus <- function(
    group_curves) {

  N_values <- vapply(
    group_curves,
    nrow,
    integer(1)
  )

  if (
    length(
      unique(
        N_values
      )
    ) != 1L
  ) {
    stop("All group curves must contain the same number of ranked features.")
  }

  N <- N_values[
    1L
  ]

  groups <- names(
    group_curves
  )

  make_matrix <- function(
      column_name) {

    m <- do.call(
      cbind,
      lapply(
        group_curves,
        function(df) {
          df[[
            column_name
          ]]
        }
      )
    )

    colnames(
      m
    ) <- groups

    m
  }

  eta_mat <- make_matrix(
    "eta_EP"
  )

  p_mat <- make_matrix(
    "p_NB"
  )

  dlogP_mat <- make_matrix(
    "d_logP_dx"
  )

  dlogE_mat <- make_matrix(
    "d_logE_dx"
  )

  dlogmu_mat <- make_matrix(
    "d_logmu_dx"
  )

  data.frame(
    rank = seq_len(
      N
    ),

    normalized_rank = (
      seq_len(
        N
      ) -
      1
    ) /
      (
        N -
        1
      ),

    eta_EP_median = rankwise_median(
      eta_mat
    ),

    eta_EP_Q25 = rankwise_quantile(
      eta_mat,
      0.25
    ),

    eta_EP_Q75 = rankwise_quantile(
      eta_mat,
      0.75
    ),

    eta_EP_n_valid = rankwise_n_finite(
      eta_mat
    ),

    p_NB_median = rankwise_median(
      p_mat
    ),

    p_NB_Q25 = rankwise_quantile(
      p_mat,
      0.25
    ),

    p_NB_Q75 = rankwise_quantile(
      p_mat,
      0.75
    ),

    p_NB_n_valid = rankwise_n_finite(
      p_mat
    ),

    d_logP_dx_median = rankwise_median(
      dlogP_mat
    ),

    d_logE_dx_median = rankwise_median(
      dlogE_mat
    ),

    d_logmu_dx_median = rankwise_median(
      dlogmu_mat
    ),

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# FIGURE HELPERS
# =============================================================================

save_four_panel <- function(
    plots,
    filename) {

  png(
    filename,
    width = 15,
    height = 14,
    units = "in",
    res = 260,
    bg = "white"
  )

  grid.newpage()

  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 4L,
        ncol = 1L
      )
    )
  )

  for (
    i in seq_along(
      plots
    )
  ) {

    print(
      plots[[i]],
      vp = viewport(
        layout.pos.row = i,
        layout.pos.col = 1L
      )
    )
  }

  dev.off()
}


finite_plot_range <- function(
    x) {

  z <- x[
    is.finite(
      x
    )
  ]

  if (length(z) == 0L) {
    return(
      c(
        -1,
        1
      )
    )
  }

  r <- range(
    z
  )

  if (
    !all(
      is.finite(
        r
      )
    ) ||
    diff(
      r
    ) <= 0
  ) {
    r <- c(
      min(
        z
      ) -
        1,
      max(
        z
      ) +
        1
    )
  }

  r
}


make_consensus_figure <- function(
    consensus_df,
    group_curves,
    out_file) {

  N <- nrow(
    consensus_df
  )

  all_eta <- bind_rows(
    lapply(
      names(
        group_curves
      ),
      function(g) {

        data.frame(
          group = g,
          rank = group_curves[[
            g
          ]]$rank,
          eta_EP = group_curves[[
            g
          ]]$eta_EP,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  all_p <- bind_rows(
    lapply(
      names(
        group_curves
      ),
      function(g) {

        data.frame(
          group = g,
          rank = group_curves[[
            g
          ]]$rank,
          p_NB = group_curves[[
            g
          ]]$p_NB,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  # ---------------------------------------------------------------------------
  # Panel A: derivatives that form the PC1–NB elasticity.
  # ---------------------------------------------------------------------------
  derivative_long <- consensus_df %>%
    select(
      rank,
      d_logP_dx_median,
      d_logE_dx_median
    ) %>%
    pivot_longer(
      cols = c(
        d_logP_dx_median,
        d_logE_dx_median
      ),
      names_to = "signal",
      values_to = "value"
    ) %>%
    mutate(
      signal = factor(
        signal,
        levels = c(
          "d_logP_dx_median",
          "d_logE_dx_median"
        ),
        labels = c(
          "d log(P) / dx",
          "d log(E) / dx"
        )
      )
    )

  p1 <- ggplot(
    derivative_long,
    aes(
      rank,
      value,
      color = signal
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.35
    ) +
    geom_line(
      linewidth = 0.85
    ) +
    labs(
      title = "A. Differential change in PC1 variance contribution and NB excess variance",
      subtitle = "Rank is explicitly ascending absolute PC1 loading",
      x = "PC1 rank: low |loading|  ->  high |loading|",
      y = "Median log-derivative across groups",
      color = NULL
    ) +
    annotate(
      "label",
      x = max(
        1,
        round(
          N * 0.03
        )
      ),
      y = Inf,
      label = "P_i = d1^2 * v_i1^2 / (n - 1)",
      hjust = 0,
      vjust = 1.2,
      size = 3.0,
      fill = "white"
    ) +
    annotate(
      "label",
      x = max(
        1,
        round(
          N * 0.03
        )
      ),
      y = Inf,
      label = "E_ig = max(V_i - mu_ig, 0)",
      hjust = 0,
      vjust = 3.0,
      size = 3.0,
      fill = "white"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  # ---------------------------------------------------------------------------
  # Panel B: PC1–NB elasticity.
  # ---------------------------------------------------------------------------
  p2 <- ggplot() +
    geom_line(
      data = all_eta,
      aes(
        rank,
        eta_EP,
        group = group
      ),
      linewidth = 0.30,
      alpha = 0.22
    ) +
    geom_ribbon(
      data = consensus_df,
      aes(
        x = rank,
        ymin = eta_EP_Q25,
        ymax = eta_EP_Q75
      ),
      alpha = 0.18
    ) +
    geom_line(
      data = consensus_df,
      aes(
        rank,
        eta_EP_median
      ),
      linewidth = 0.95
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.55
    ) +
    labs(
      title = "B. PC1–NB elasticity",
      subtitle = "eta_EP = 1 indicates proportional local log-change; >1 means E changes faster than P",
      x = "PC1 rank",
      y = "eta_EP(r)"
    ) +
    annotate(
      "label",
      x = max(
        1,
        round(
          N * 0.03
        )
      ),
      y = Inf,
      label = "eta_EP(r) = [d log(E)/dx] / [d log(P)/dx]",
      hjust = 0,
      vjust = 1.2,
      size = 3.1,
      fill = "white"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank()
    )

  # ---------------------------------------------------------------------------
  # Panel C: mean-variance exponent.
  # ---------------------------------------------------------------------------
  p3 <- ggplot() +
    geom_line(
      data = all_p,
      aes(
        rank,
        p_NB,
        group = group
      ),
      linewidth = 0.30,
      alpha = 0.22
    ) +
    geom_ribbon(
      data = consensus_df,
      aes(
        x = rank,
        ymin = p_NB_Q25,
        ymax = p_NB_Q75
      ),
      alpha = 0.18
    ) +
    geom_line(
      data = consensus_df,
      aes(
        rank,
        p_NB_median
      ),
      linewidth = 0.95
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.55
    ) +
    geom_hline(
      yintercept = 2,
      linetype = "dotted",
      linewidth = 0.65
    ) +
    labs(
      title = "C. Local effective negative-binomial mean–variance exponent",
      subtitle = "Under locally stable alpha: p approximately 1 is NB1-like; p approximately 2 is NB2-like",
      x = "PC1 rank",
      y = "p_NB(r)"
    ) +
    annotate(
      "label",
      x = max(
        1,
        round(
          N * 0.03
        )
      ),
      y = Inf,
      label = "E = alpha * mu^p     =>     p_NB(r) = [d log(E)/dx] / [d log(mu)/dx]",
      hjust = 0,
      vjust = 1.2,
      size = 3.0,
      fill = "white"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank()
    )

  # ---------------------------------------------------------------------------
  # Panel D: number of groups contributing finite differential estimates.
  # This makes denominator-zero regions visible rather than silently hiding them.
  # ---------------------------------------------------------------------------
  validity_long <- consensus_df %>%
    select(
      rank,
      eta_EP_n_valid,
      p_NB_n_valid
    ) %>%
    pivot_longer(
      cols = c(
        eta_EP_n_valid,
        p_NB_n_valid
      ),
      names_to = "field",
      values_to = "n_valid"
    ) %>%
    mutate(
      field = factor(
        field,
        levels = c(
          "eta_EP_n_valid",
          "p_NB_n_valid"
        ),
        labels = c(
          "eta_EP",
          "p_NB"
        )
      )
    )

  p4 <- ggplot(
    validity_long,
    aes(
      rank,
      n_valid,
      color = field
    )
  ) +
    geom_line(
      linewidth = 0.85
    ) +
    scale_y_continuous(
      breaks = 0:8,
      limits = c(
        0,
        8
      )
    ) +
    labs(
      title = "D. Support for the differential estimates",
      subtitle = "Number of experimental groups with a mathematically defined derivative ratio at each rank",
      x = "PC1 rank",
      y = "Number of groups",
      color = NULL
    ) +
    annotate(
      "label",
      x = max(
        1,
        round(
          N * 0.03
        )
      ),
      y = 8,
      label = "x = (rank - 1)/(N - 1); ratios are NA where the denominator derivative is numerically zero",
      hjust = 0,
      vjust = 1.2,
      size = 2.8,
      fill = "white"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  save_four_panel(
    list(
      p1,
      p2,
      p3,
      p4
    ),
    out_file
  )
}


make_group_figure <- function(
    feature_df,
    group_name,
    out_file) {

  N <- nrow(
    feature_df
  )

  deriv_long <- feature_df %>%
    select(
      rank,
      d_logP_dx,
      d_logE_dx,
      d_logmu_dx
    ) %>%
    pivot_longer(
      cols = c(
        d_logP_dx,
        d_logE_dx,
        d_logmu_dx
      ),
      names_to = "quantity",
      values_to = "value"
    )

  p1 <- ggplot(
    deriv_long,
    aes(
      rank,
      value,
      color = quantity
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.35
    ) +
    geom_line(
      linewidth = 0.72
    ) +
    labs(
      title = paste0(
        group_name,
        ": log-derivatives along PC1 rank"
      ),
      subtitle = "Rank = ascending absolute PC1 loading",
      x = "PC1 rank",
      y = "Derivative with respect to normalized rank x",
      color = NULL
    ) +
    annotate(
      "label",
      x = max(
        1,
        round(
          N * 0.03
        )
      ),
      y = Inf,
      label = "P_i = d1^2 * v_i1^2/(n-1);   E_ig = max(V_i - mu_ig, 0)",
      hjust = 0,
      vjust = 1.2,
      size = 2.9,
      fill = "white"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  p2 <- ggplot(
    feature_df,
    aes(
      rank,
      eta_EP
    )
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.55
    ) +
    geom_line(
      linewidth = 0.8
    ) +
    labs(
      title = paste0(
        group_name,
        ": PC1–NB elasticity"
      ),
      subtitle = "eta_EP > 1: excess variance changes faster than PC1 contribution on the local log scale",
      x = "PC1 rank",
      y = "eta_EP(r)"
    ) +
    annotate(
      "label",
      x = max(
        1,
        round(
          N * 0.03
        )
      ),
      y = Inf,
      label = "eta_EP(r) = [d log(E)/dx] / [d log(P)/dx]",
      hjust = 0,
      vjust = 1.2,
      size = 3.0,
      fill = "white"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank()
    )

  p3 <- ggplot(
    feature_df,
    aes(
      rank,
      p_NB
    )
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.55
    ) +
    geom_hline(
      yintercept = 2,
      linetype = "dotted",
      linewidth = 0.65
    ) +
    geom_line(
      linewidth = 0.8
    ) +
    labs(
      title = paste0(
        group_name,
        ": local effective NB exponent"
      ),
      subtitle = "Under locally stable alpha: p about 1 is NB1-like; p about 2 is NB2-like",
      x = "PC1 rank",
      y = "p_NB(r)"
    ) +
    annotate(
      "label",
      x = max(
        1,
        round(
          N * 0.03
        )
      ),
      y = Inf,
      label = "E = alpha * mu^p;   p_NB(r) = [d log(E)/dx] / [d log(mu)/dx]",
      hjust = 0,
      vjust = 1.2,
      size = 2.9,
      fill = "white"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank()
    )

  valid_df <- data.frame(
    rank = feature_df$rank,
    eta_defined = as.integer(
      is.finite(
        feature_df$eta_EP
      )
    ),
    p_defined = as.integer(
      is.finite(
        feature_df$p_NB
      )
    )
  ) %>%
    pivot_longer(
      cols = c(
        eta_defined,
        p_defined
      ),
      names_to = "field",
      values_to = "defined"
    )

  p4 <- ggplot(
    valid_df,
    aes(
      rank,
      defined,
      color = field
    )
  ) +
    geom_line(
      linewidth = 0.65
    ) +
    scale_y_continuous(
      breaks = c(
        0,
        1
      ),
      labels = c(
        "undefined",
        "defined"
      )
    ) +
    labs(
      title = paste0(
        group_name,
        ": mathematical support"
      ),
      subtitle = "Undefined points are retained as NA rather than replaced by extreme ratios",
      x = "PC1 rank",
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

  save_four_panel(
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
# RUN
# =============================================================================

count_mat <- read_count_matrix(
  COUNT_FILE,
  GROUP_PATTERNS
)

group_labels <- assign_sample_groups(
  sample_names = colnames(
    count_mat
  ),
  group_patterns = GROUP_PATTERNS
)

message(
  "Using count file: ",
  COUNT_FILE
)

message(
  "Count matrix dimensions: ",
  nrow(
    count_mat
  ),
  " features x ",
  ncol(
    count_mat
  ),
  " samples"
)

message(
  "Experimental groups: ",
  paste(
    levels(
      group_labels
    ),
    collapse = ", "
  )
)


# -----------------------------------------------------------------------------
# PC1 ranking matrix:
# CPM -> log1p -> arm-specific PC1 -> ascending |loading|.
# -----------------------------------------------------------------------------

rank_matrix_all <- normalize_cpm_log1p(
  count_mat
)


# -----------------------------------------------------------------------------
# Sequencing mean/variance matrix:
# DESeq2-normalized counts across the complete experiment.
# -----------------------------------------------------------------------------

deseq2_norm <- normalize_deseq2(
  count_mat = count_mat,
  group_labels = group_labels
)

normalized_counts <- deseq2_norm$normalized_counts


# -----------------------------------------------------------------------------
# Pooled within-group sequencing variance.
# -----------------------------------------------------------------------------

pooled <- compute_pooled_within_group_variance(
  normalized_counts = normalized_counts,
  group_labels = group_labels
)

pooled_variance <- pooled$variance

message(
  "Pooled within-group residual degrees of freedom: ",
  pooled$residual_df
)


# -----------------------------------------------------------------------------
# One differential-geometry trajectory per experimental group.
# -----------------------------------------------------------------------------

groups <- levels(
  group_labels
)

group_curves <- vector(
  "list",
  length(
    groups
  )
)

names(
  group_curves
) <- groups

summary_rows <- vector(
  "list",
  length(
    groups
  )
)

names(
  summary_rows
) <- groups


for (
  g in groups
) {

  idx <- which(
    group_labels == g
  )

  one <- compute_group_geometry(
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

  group_curves[[
    g
  ]] <- one$features

  summary_rows[[
    g
  ]] <- one$summary

  write.csv(
    one$features,
    file.path(
      OUT_ROOT,
      paste0(
        "Table_DifferentialGeometry_",
        g,
        ".csv"
      )
    ),
    row.names = FALSE
  )

  make_group_figure(
    feature_df = one$features,
    group_name = g,
    out_file = file.path(
      OUT_ROOT,
      paste0(
        "Figure_DifferentialGeometry_",
        g,
        ".png"
      )
    )
  )

  message(
    "Completed group ",
    g,
    ": eta defined at ",
    round(
      100 *
        one$summary$eta_valid_rank_fraction[
          1L
        ],
      1
    ),
    "% of ranks; p_NB defined at ",
    round(
      100 *
        one$summary$p_valid_rank_fraction[
          1L
        ],
      1
    ),
    "% of ranks."
  )
}


pc1_spline_summary <- bind_rows(
  summary_rows
)

write.csv(
  pc1_spline_summary,
  file.path(
    OUT_ROOT,
    "Table_PC1_and_Spline_Summary.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# Rank-wise consensus across all eight groups.
# -----------------------------------------------------------------------------

consensus_df <- compute_consensus(
  group_curves
)

write.csv(
  consensus_df,
  file.path(
    OUT_ROOT,
    "Table_Consensus_DifferentialGeometry.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# Sample normalization metadata.
# -----------------------------------------------------------------------------

write.csv(
  data.frame(
    sample = colnames(
      count_mat
    ),
    group = as.character(
      group_labels
    ),
    DESeq2_size_factor = as.numeric(
      deseq2_norm$size_factors
    ),
    stringsAsFactors = FALSE
  ),
  file.path(
    OUT_ROOT,
    "Table_Samples_and_SizeFactors.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# Main consensus mathematical figure.
# -----------------------------------------------------------------------------

make_consensus_figure(
  consensus_df = consensus_df,
  group_curves = group_curves,
  out_file = file.path(
    OUT_ROOT,
    "Figure_Consensus_PC1_NB_DifferentialGeometry.png"
  )
)


# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

message(
  "============================================================"
)

message(
  "PC1-NB DIFFERENTIAL GEOMETRY ANALYSIS COMPLETE"
)

message(
  "Ranking used: ascending absolute PC1 loading."
)

message(
  "Primary differential field:"
)

message(
  "  eta_EP(r) = [d log(E)/dx] / [d log(P)/dx]"
)

message(
  "Local effective NB exponent:"
)

message(
  "  p_NB(r) = [d log(E)/dx] / [d log(mu)/dx]"
)

message(
  "No cutoff was estimated."
)

message(
  "No bootstrap was performed."
)

message(
  "Inspect first:"
)

message(
  "  Figure_Consensus_PC1_NB_DifferentialGeometry.png"
)

message(
  "Then inspect:"
)

message(
  "  Table_Consensus_DifferentialGeometry.csv"
)

message(
  "Output directory: ",
  OUT_ROOT
)

message(
  "============================================================"
)
