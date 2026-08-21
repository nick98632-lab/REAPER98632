#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# FINAL PC1–NB CUMULATIVE-DIVERGENCE ANALYSIS
# =============================================================================
#
# PURPOSE
# -------
# 1. Rank features within each experimental group by ascending |PC1 loading|.
# 2. Quantify how PC1 variance mass and NB excess-variance mass are distributed
#    over that ranked axis.
# 3. Define cumulative divergence:
#
#       p_g(r) = P_g(r) / sum_j P_g(j)
#       q_g(r) = E_g(r) / sum_j E_g(j)
#
#       F_P,g(r) = sum_{j<=r} p_g(j)
#       F_E,g(r) = sum_{j<=r} q_g(j)
#
#       D_g(r) = F_E,g(r) - F_P,g(r)
#
#    Therefore, in the discrete ranked system:
#
#       Delta D_g(r) = q_g(r) - p_g(r)
#
# 4. Fit a continuous three-regime linear-spline model to all 8 group-specific
#    D_g(r) curves simultaneously, with TWO SHARED knots c1 and c2:
#
#       D_g(x) =
#         beta_0g + beta_1g*x
#         + gamma_1g*(x-c1)_+
#         + gamma_2g*(x-c2)_+
#
#    where x = (rank-1)/(N-1) and (z)_+ = max(z,0).
#
#    c1 and c2 are estimated jointly from the complete experiment by minimizing
#    the summed residual squared error over all eight groups.
#
# 5. Interpret:
#
#       rank < c1        = REMAINDER
#       c1 <= rank <= c2 = DIVERGENCE INTERVAL
#       rank > c2        = LEADING EDGE
#
# 6. AFTER the boundaries are defined, characterize the regions using:
#
#       E = alpha * mu^p
#
#       p ~ 1  -> NB1-like scaling
#       p ~ 2  -> NB2-like scaling
#
#    plus the original descriptive NB2-related diagnostics.
#
# 7. Produce ONE clean three-panel manuscript figure.
#
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <-
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <-
  "/root/REAPER98632/exports/PC1_NB_CUMULATIVE_DIVERGENCE_FINAL"

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

# Preserve the visual smoothing convention of the original empirical-variance
# figure. This smoothing is for display only; it does NOT determine c1 or c2.
VAR_SPLINE_SPAR <- 0.60

PNG_WIDTH_IN  <- 15
PNG_HEIGHT_IN <- 11.5
PNG_DPI       <- 300

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# COLORS
# =============================================================================

COL <- list(
  empirical   = "#117A65",
  pc1         = "#386CB0",
  nb          = "#159D91",
  divergence  = "#222222",
  fit         = "#111111",
  remainder   = "#DCEFF2",
  interval    = "#F5E8C8",
  leading     = "#DDF2EA",
  c1          = "#2166AC",
  c2          = "#1B7837",
  ribbon      = "#BDBDBD"
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


assign_groups <- function(
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

    assigned[
      idx
    ] <- group_name
  }

  if (any(is.na(assigned))) {
    stop(
      "Unassigned sample columns: ",
      paste(
        sample_names[
          is.na(assigned)
        ],
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
# DESIGN-AWARE POOLED WITHIN-GROUP VARIANCE
# =============================================================================

compute_pooled_variance <- function(
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

    resid_g <- sweep(
      xg,
      1L,
      mu_g,
      "-"
    )

    sse <- sse +
      rowSums(
        resid_g^2
      )

    residual_df <- residual_df +
      length(idx) -
      1L
  }

  if (residual_df < 2L) {
    stop("Pooled residual degrees of freedom < 2.")
  }

  V <- sse /
    residual_df

  V[
    !is.finite(V)
  ] <- 0

  V <- pmax(
    V,
    0
  )

  names(V) <- rownames(
    normalized_counts
  )

  list(
    variance = V,
    residual_df = residual_df
  )
}


# =============================================================================
# PC1 + RANK-WISE VARIANCE MASSES
# =============================================================================

compute_group_curves <- function(
    group_name,
    raw_counts_arm,
    rank_matrix_arm,
    normalized_counts_arm,
    pooled_variance) {

  # ---------------------------------------------------------------------------
  # PC1 is computed exactly in the original orientation:
  #
  #   samples x features
  #
  # and features are ranked explicitly by ascending absolute loading.
  # ---------------------------------------------------------------------------
  pca <- stats::prcomp(
    t(
      rank_matrix_arm
    ),
    center = TRUE,
    scale. = FALSE,
    rank. = 1
  )

  loading <- pca$rotation[
    ,
    1L
  ]

  loading[
    !is.finite(
      loading
    )
  ] <- 0

  abs_loading <- abs(
    loading
  )

  rank_order <- order(
    abs_loading,
    decreasing = FALSE
  )

  # prcomp$sdev[1]^2 is the PC1 eigenvalue (variance represented by PC1).
  lambda1 <- pca$sdev[
    1L
  ]^2

  # Feature-specific PC1 variance contribution:
  #
  #   P_i = lambda1 * v_i1^2
  #
  #       = d1^2 * v_i1^2 / (n-1)
  #
  P <- lambda1 *
    loading^2

  # Within-group normalized mean for the NB excess-variance term.
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

  # Empirical excess-over-Poisson variance:
  #
  #   E_ig = max(V_i - mu_ig, 0)
  #
  E <- pmax(
    pooled_variance -
    mu,
    0
  )

  # Original raw empirical variance geometry retained for Panel A.
  raw_var <- apply(
    raw_counts_arm,
    1L,
    stats::var,
    na.rm = TRUE
  )

  raw_var[
    !is.finite(raw_var)
  ] <- 0

  raw_var <- pmax(
    raw_var,
    0
  )

  P_r <- P[
    rank_order
  ]

  E_r <- E[
    rank_order
  ]

  mu_r <- mu[
    rank_order
  ]

  V_r <- pooled_variance[
    rank_order
  ]

  raw_var_r <- raw_var[
    rank_order
  ]

  loading_r <- loading[
    rank_order
  ]

  abs_loading_r <- abs_loading[
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

  # ---------------------------------------------------------------------------
  # Convert P and E into probability masses over PC1 rank.
  # ---------------------------------------------------------------------------
  P_total <- sum(
    P_r,
    na.rm = TRUE
  )

  E_total <- sum(
    E_r,
    na.rm = TRUE
  )

  if (
    !is.finite(P_total) ||
    P_total <= 0
  ) {
    stop(
      "PC1 variance mass is undefined for group ",
      group_name
    )
  }

  if (
    !is.finite(E_total) ||
    E_total <= 0
  ) {
    stop(
      "NB excess-variance mass is undefined for group ",
      group_name
    )
  }

  p_mass <- P_r /
    P_total

  q_mass <- E_r /
    E_total

  F_P <- cumsum(
    p_mass
  )

  F_E <- cumsum(
    q_mass
  )

  D <- F_E -
    F_P

  delta_D <- q_mass -
    p_mass

  # ---------------------------------------------------------------------------
  # Smooth the ORIGINAL log1p(raw empirical variance) curve for display only.
  # ---------------------------------------------------------------------------
  log_raw_var <- log1p(
    raw_var_r
  )

  var_spline <- stats::smooth.spline(
    x = rank,
    y = log_raw_var,
    spar = VAR_SPLINE_SPAR
  )

  smooth_log_raw_var <- as.numeric(
    stats::predict(
      var_spline,
      x = rank,
      deriv = 0
    )$y
  )

  # ---------------------------------------------------------------------------
  # Original moment-based NB diagnostics retained for corroboration.
  # ---------------------------------------------------------------------------
  alpha <- rep(
    0,
    N
  )

  positive_mu <- mu_r > 0

  alpha[
    positive_mu
  ] <- E_r[
    positive_mu
  ] /
    (
      mu_r[
        positive_mu
      ]^2
    )

  data.frame(
    group = group_name,
    rank = rank,
    normalized_rank = x,
    feature_id = rownames(
      rank_matrix_arm
    )[
      rank_order
    ],
    pc1_loading = loading_r,
    abs_pc1_loading = abs_loading_r,
    pc1_eigenvalue = lambda1,
    pc1_variance_contribution = P_r,
    pc1_variance_mass = p_mass,
    group_mean_normalized = mu_r,
    pooled_within_group_variance = V_r,
    nb_excess_variance = E_r,
    nb_excess_variance_mass = q_mass,
    cumulative_pc1_mass = F_P,
    cumulative_nb_mass = F_E,
    cumulative_divergence = D,
    local_mass_difference = delta_D,
    raw_empirical_variance = raw_var_r,
    log1p_raw_empirical_variance = log_raw_var,
    smooth_log1p_raw_empirical_variance = smooth_log_raw_var,
    NB2 = log1p(
      E_r
    ),
    NB2_NB1 = log1p(
      E_r
    ) -
      log1p(
        mu_r
      ),
    alpha_mu = log1p(
      alpha *
      mu_r
    ),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# SHARED TWO-KNOT CONTINUOUS LINEAR-SPLINE FIT
# =============================================================================

piecewise_basis <- function(
    x,
    c1,
    c2) {

  cbind(
    intercept = 1,
    x = x,
    hinge1 = pmax(
      x -
        c1,
      0
    ),
    hinge2 = pmax(
      x -
        c2,
      0
    )
  )
}


piecewise_sse <- function(
    par,
    x,
    D_mat,
    min_gap) {

  c1 <- par[
    1L
  ]

  c2 <- par[
    2L
  ]

  if (
    !is.finite(c1) ||
    !is.finite(c2) ||
    c1 <= 0 ||
    c2 >= 1 ||
    c2 -
      c1 <= min_gap
  ) {
    return(
      1e100
    )
  }

  X <- piecewise_basis(
    x,
    c1,
    c2
  )

  qrX <- qr(
    X
  )

  coef <- tryCatch(
    qr.coef(
      qrX,
      D_mat
    ),
    error = function(e) {
      NULL
    }
  )

  if (
    is.null(coef) ||
    any(
      !is.finite(
        coef
      )
    )
  ) {
    return(
      1e100
    )
  }

  fitted <- X %*%
    coef

  resid <- D_mat -
    fitted

  sum(
    resid^2
  )
}


fit_shared_knots <- function(
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
      unique(
        N_values
      )
    ) != 1L
  ) {
    stop("All groups must contain the same number of ranked features.")
  }

  N <- N_values[
    1L
  ]

  x_full <- (
    seq_len(
      N
    ) -
      1
  ) /
    (
      N -
      1
    )

  D_full <- do.call(
    cbind,
    lapply(
      group_curves,
      function(df) {
        df$cumulative_divergence
      }
    )
  )

  colnames(
    D_full
  ) <- groups

  # Multi-start optimization uses an evenly spaced representation of the full
  # cumulative curves only to identify the best optimization basin.
  opt_n <- min(
    N,
    5000L
  )

  opt_idx <- unique(
    as.integer(
      round(
        seq(
          1,
          N,
          length.out = opt_n
        )
      )
    )
  )

  x_opt <- x_full[
    opt_idx
  ]

  D_opt <- D_full[
    opt_idx,
    ,
    drop = FALSE
  ]

  # Only enough separation to make the two-knot model identifiable.
  min_gap <- max(
    4 /
      (
        N -
          1
      ),
    .Machine$double.eps^0.25
  )

  starts <- list(
    c(
      0.03,
      0.97
    ),
    c(
      0.08,
      0.92
    ),
    c(
      0.15,
      0.85
    ),
    c(
      0.25,
      0.75
    ),
    c(
      0.35,
      0.65
    )
  )

  coarse_results <- lapply(
    starts,
    function(start) {

      stats::optim(
        par = start,
        fn = piecewise_sse,
        x = x_opt,
        D_mat = D_opt,
        min_gap = min_gap,
        method = "Nelder-Mead",
        control = list(
          maxit = 700,
          reltol = 1e-11
        )
      )
    }
  )

  coarse_values <- vapply(
    coarse_results,
    function(z) {
      z$value
    },
    numeric(1)
  )

  if (
    !any(
      is.finite(
        coarse_values
      )
    )
  ) {
    stop("Shared-knot optimization failed at the multi-start stage.")
  }

  best_coarse <- coarse_results[[
    which.min(
      coarse_values
    )
  ]]

  # Refine the best solution on ALL ranked features and ALL eight groups.
  refined <- stats::optim(
    par = best_coarse$par,
    fn = piecewise_sse,
    x = x_full,
    D_mat = D_full,
    min_gap = min_gap,
    method = "Nelder-Mead",
    control = list(
      maxit = 1000,
      reltol = 1e-12
    )
  )

  if (
    !is.finite(
      refined$value
    )
  ) {
    stop("Full-data shared-knot refinement failed.")
  }

  c1_x <- refined$par[
    1L
  ]

  c2_x <- refined$par[
    2L
  ]

  c1_rank <- as.integer(
    round(
      1 +
        c1_x *
        (
          N -
            1
        )
    )
  )

  c2_rank <- as.integer(
    round(
      1 +
        c2_x *
        (
          N -
            1
        )
    )
  )

  c1_rank <- max(
    2L,
    min(
      N -
        2L,
      c1_rank
    )
  )

  c2_rank <- max(
    c1_rank +
      1L,
    min(
      N -
        1L,
      c2_rank
    )
  )

  # Refit coefficients at the exact integer ranks reported in the manuscript.
  c1_x_final <- (
    c1_rank -
      1
  ) /
    (
      N -
        1
    )

  c2_x_final <- (
    c2_rank -
      1
  ) /
    (
      N -
        1
    )

  X_final <- piecewise_basis(
    x_full,
    c1_x_final,
    c2_x_final
  )

  qr_final <- qr(
    X_final
  )

  coef_final <- qr.coef(
    qr_final,
    D_full
  )

  fitted_final <- X_final %*%
    coef_final

  resid_final <- D_full -
    fitted_final

  sse_final <- sum(
    resid_final^2
  )

  list(
    c1_rank = c1_rank,
    c2_rank = c2_rank,
    c1_x = c1_x_final,
    c2_x = c2_x_final,
    coefficients = coef_final,
    fitted = fitted_final,
    D_matrix = D_full,
    x = x_full,
    SSE = sse_final,
    groups = groups,
    optimizer = refined
  )
}


# =============================================================================
# REGION ASSIGNMENT + NB CHARACTERIZATION
# =============================================================================

assign_region <- function(
    rank,
    c1,
    c2) {

  ifelse(
    rank < c1,
    "REMAINDER",
    ifelse(
      rank <= c2,
      "DIVERGENCE_INTERVAL",
      "LEADING_EDGE"
    )
  )
}


estimate_nb_exponent <- function(
    mu,
    E) {

  keep <- (
    is.finite(mu) &
    is.finite(E) &
    mu > 0 &
    E > 0
  )

  if (
    sum(
      keep
    ) < 10L
  ) {
    return(
      data.frame(
        n_positive = sum(
          keep
        ),
        p = NA_real_,
        p_se = NA_real_,
        r_squared = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  fit <- stats::lm(
    log(
      E[
        keep
      ]
    ) ~
      log(
        mu[
          keep
        ]
      )
  )

  fit_summary <- summary(
    fit
  )

  data.frame(
    n_positive = sum(
      keep
    ),
    p = unname(
      coef(
        fit
      )[
        2L
      ]
    ),
    p_se = unname(
      fit_summary$coefficients[
        2L,
        2L
      ]
    ),
    r_squared = fit_summary$r.squared,
    stringsAsFactors = FALSE
  )
}


summarize_regions <- function(
    group_curves,
    c1,
    c2) {

  exponent_rows <- list()
  metric_rows <- list()

  for (
    g in names(
      group_curves
    )
  ) {

    df <- group_curves[[
      g
    ]] %>%
      mutate(
        region = assign_region(
          rank,
          c1,
          c2
        )
      )

    for (
      reg in c(
        "REMAINDER",
        "DIVERGENCE_INTERVAL",
        "LEADING_EDGE"
      )
    ) {

      sub <- df[
        df$region == reg,
        ,
        drop = FALSE
      ]

      p_fit <- estimate_nb_exponent(
        mu = sub$group_mean_normalized,
        E = sub$nb_excess_variance
      )

      exponent_rows[[
        length(
          exponent_rows
        ) +
          1L
      ]] <- data.frame(
        group = g,
        region = reg,
        p_fit,
        stringsAsFactors = FALSE
      )

      metric_rows[[
        length(
          metric_rows
        ) +
          1L
      ]] <- data.frame(
        group = g,
        region = reg,
        n_features = nrow(
          sub
        ),
        median_NB2 = median(
          sub$NB2,
          na.rm = TRUE
        ),
        median_NB2_NB1 = median(
          sub$NB2_NB1,
          na.rm = TRUE
        ),
        median_alpha_mu = median(
          sub$alpha_mu,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  exponent_df <- bind_rows(
    exponent_rows
  )

  metric_df <- bind_rows(
    metric_rows
  )

  exponent_summary <- exponent_df %>%
    group_by(
      region
    ) %>%
    summarise(
      groups_with_estimate = sum(
        is.finite(
          p
        )
      ),
      p_median = median(
        p,
        na.rm = TRUE
      ),
      p_Q25 = as.numeric(
        quantile(
          p,
          0.25,
          na.rm = TRUE,
          names = FALSE
        )
      ),
      p_Q75 = as.numeric(
        quantile(
          p,
          0.75,
          na.rm = TRUE,
          names = FALSE
        )
      ),
      .groups = "drop"
    )

  list(
    exponent_by_group = exponent_df,
    exponent_summary = exponent_summary,
    metrics_by_group = metric_df
  )
}


# =============================================================================
# CONSENSUS TABLES FOR THE ONE MANUSCRIPT FIGURE
# =============================================================================

rankwise_stat <- function(
    group_curves,
    column,
    fun) {

  mat <- do.call(
    cbind,
    lapply(
      group_curves,
      function(df) {
        df[[
          column
        ]]
      }
    )
  )

  apply(
    mat,
    1L,
    fun
  )
}


build_consensus_table <- function(
    group_curves,
    knot_fit) {

  N <- nrow(
    group_curves[[
      1L
    ]]
  )

  fitted_median <- apply(
    knot_fit$fitted,
    1L,
    median,
    na.rm = TRUE
  )

  data.frame(
    rank = seq_len(
      N
    ),

    empirical_median = rankwise_stat(
      group_curves,
      "smooth_log1p_raw_empirical_variance",
      median
    ),

    empirical_Q25 = rankwise_stat(
      group_curves,
      "smooth_log1p_raw_empirical_variance",
      function(z) {
        as.numeric(
          quantile(
            z,
            0.25,
            names = FALSE,
            na.rm = TRUE
          )
        )
      }
    ),

    empirical_Q75 = rankwise_stat(
      group_curves,
      "smooth_log1p_raw_empirical_variance",
      function(z) {
        as.numeric(
          quantile(
            z,
            0.75,
            names = FALSE,
            na.rm = TRUE
          )
        )
      }
    ),

    F_P_median = rankwise_stat(
      group_curves,
      "cumulative_pc1_mass",
      median
    ),

    F_E_median = rankwise_stat(
      group_curves,
      "cumulative_nb_mass",
      median
    ),

    D_median = rankwise_stat(
      group_curves,
      "cumulative_divergence",
      median
    ),

    D_Q25 = rankwise_stat(
      group_curves,
      "cumulative_divergence",
      function(z) {
        as.numeric(
          quantile(
            z,
            0.25,
            names = FALSE,
            na.rm = TRUE
          )
        )
      }
    ),

    D_Q75 = rankwise_stat(
      group_curves,
      "cumulative_divergence",
      function(z) {
        as.numeric(
          quantile(
            z,
            0.75,
            names = FALSE,
            na.rm = TRUE
          )
        )
      }
    ),

    local_mass_difference_median = rankwise_stat(
      group_curves,
      "local_mass_difference",
      median
    ),

    D_piecewise_fit_median = fitted_median,

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# FIGURE HELPERS
# =============================================================================

add_region_background <- function(
    p,
    c1,
    c2,
    N) {

  p +
    annotate(
      "rect",
      xmin = 1,
      xmax = c1,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$remainder,
      alpha = 0.55
    ) +
    annotate(
      "rect",
      xmin = c1,
      xmax = c2,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$interval,
      alpha = 0.42
    ) +
    annotate(
      "rect",
      xmin = c2,
      xmax = N,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$leading,
      alpha = 0.55
    )
}


get_p_summary <- function(
    exponent_summary,
    region) {

  row <- exponent_summary[
    exponent_summary$region == region,
    ,
    drop = FALSE
  ]

  if (
    nrow(
      row
    ) == 0L ||
    !is.finite(
      row$p_median[
        1L
      ]
    )
  ) {
    return(
      "NA"
    )
  }

  paste0(
    sprintf(
      "%.2f",
      row$p_median[
        1L
      ]
    ),
    " [",
    sprintf(
      "%.2f",
      row$p_Q25[
        1L
      ]
    ),
    ", ",
    sprintf(
      "%.2f",
      row$p_Q75[
        1L
      ]
    ),
    "]"
  )
}


make_main_figure <- function(
    consensus_df,
    c1,
    c2,
    exponent_summary,
    out_file) {

  N <- nrow(
    consensus_df
  )

  # ---------------------------------------------------------------------------
  # PANEL A: original empirical-variance geometry.
  # ---------------------------------------------------------------------------
  pA <- ggplot(
    consensus_df,
    aes(
      rank,
      empirical_median
    )
  )

  pA <- add_region_background(
    pA,
    c1,
    c2,
    N
  )

  pA <- pA +
    geom_ribbon(
      aes(
        ymin = empirical_Q25,
        ymax = empirical_Q75
      ),
      fill = COL$empirical,
      alpha = 0.12
    ) +
    geom_line(
      color = COL$empirical,
      linewidth = 1.0
    ) +
    geom_vline(
      xintercept = c1,
      color = COL$c1,
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = c2,
      color = COL$c2,
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    labs(
      title = "A. PC1-ranked empirical sequencing-variance geometry",
      subtitle = paste0(
        "Features ranked ascending by |PC1 loading|; shared c1 = ",
        c1,
        ", shared c2 = ",
        c2
      ),
      x = "PC1 rank: low |loading|  ->  high |loading|",
      y = "Median smoothed log(1 + empirical variance)"
    ) +
    annotate(
      "label",
      x = max(
        1,
        round(
          0.025 *
            N
        )
      ),
      y = Inf,
      label = "P[i] == lambda[1]*v[i*1]^2 == d[1]^2*v[i*1]^2/(n-1)",
      parse = TRUE,
      hjust = 0,
      vjust = 1.2,
      size = 3.0,
      fill = "white"
    ) +
    annotate(
      "text",
      x = max(
        1,
        round(
          c1 /
            2
        )
      ),
      y = -Inf,
      label = "REMAINDER",
      color = COL$c1,
      vjust = -0.7,
      fontface = "bold",
      size = 3.2
    ) +
    annotate(
      "text",
      x = round(
        (
          c1 +
            c2
        ) /
          2
      ),
      y = -Inf,
      label = "DIVERGENCE INTERVAL",
      vjust = -0.7,
      fontface = "bold",
      size = 3.2
    ) +
    annotate(
      "text",
      x = round(
        (
          c2 +
            N
        ) /
          2
      ),
      y = -Inf,
      label = "LEADING EDGE",
      color = COL$c2,
      vjust = -0.7,
      fontface = "bold",
      size = 3.2
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        face = "bold"
      )
    )

  # ---------------------------------------------------------------------------
  # PANEL B: cumulative variance mass and divergence.
  # ---------------------------------------------------------------------------
  cumulative_long <- consensus_df %>%
    select(
      rank,
      F_P_median,
      F_E_median
    ) %>%
    pivot_longer(
      cols = c(
        F_P_median,
        F_E_median
      ),
      names_to = "curve",
      values_to = "value"
    ) %>%
    mutate(
      curve = factor(
        curve,
        levels = c(
          "F_P_median",
          "F_E_median"
        ),
        labels = c(
          "Cumulative PC1 variance mass",
          "Cumulative NB excess-variance mass"
        )
      )
    )

  pB <- ggplot()

  pB <- add_region_background(
    pB,
    c1,
    c2,
    N
  )

  pB <- pB +
    geom_ribbon(
      data = consensus_df,
      aes(
        x = rank,
        ymin = D_Q25,
        ymax = D_Q75
      ),
      fill = COL$ribbon,
      alpha = 0.25
    ) +
    geom_line(
      data = cumulative_long,
      aes(
        rank,
        value,
        color = curve
      ),
      linewidth = 0.9
    ) +
    geom_line(
      data = consensus_df,
      aes(
        rank,
        D_median
      ),
      color = COL$divergence,
      linewidth = 1.05
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.45
    ) +
    geom_vline(
      xintercept = c1,
      color = COL$c1,
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = c2,
      color = COL$c2,
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    scale_color_manual(
      values = c(
        "Cumulative PC1 variance mass" = COL$pc1,
        "Cumulative NB excess-variance mass" = COL$nb
      )
    ) +
    labs(
      title = "B. Cumulative PC1–NB variance-mass divergence",
      subtitle = "Black = D(r); grey ribbon = armwise IQR of D(r)",
      x = "PC1 rank",
      y = "Dimensionless cumulative mass / divergence",
      color = NULL
    ) +
    annotate(
      "label",
      x = max(
        1,
        round(
          0.025 *
            N
        )
      ),
      y = 0.96,
      label = "p(r) == P(r)/sum[j](P(j))~~','~~q(r) == E(r)/sum[j](E(j))",
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 2.9,
      fill = "white"
    ) +
    annotate(
      "label",
      x = max(
        1,
        round(
          0.025 *
            N
        )
      ),
      y = 0.74,
      label = "D(r) == F[E](r)-F[P](r)~~','~~Delta*D(r) == q(r)-p(r)",
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 2.9,
      fill = "white"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(
        face = "bold"
      )
    )

  # ---------------------------------------------------------------------------
  # PANEL C: actual divergence + jointly fitted shared two-knot model.
  # ---------------------------------------------------------------------------
  rem_p <- get_p_summary(
    exponent_summary,
    "REMAINDER"
  )

  mid_p <- get_p_summary(
    exponent_summary,
    "DIVERGENCE_INTERVAL"
  )

  lead_p <- get_p_summary(
    exponent_summary,
    "LEADING_EDGE"
  )

  nb_box <- paste0(
    "AFTER boundaries are defined:\n",
    "E = alpha * mu^p\n",
    "p ~ 1  ->  NB1-like\n",
    "p ~ 2  ->  NB2-like\n\n",
    "Median p [IQR] across 8 groups\n",
    "Remainder:  ", rem_p, "\n",
    "Interval:   ", mid_p, "\n",
    "Leading:    ", lead_p
  )

  equation_box <-
    "D[g](x) == beta[0*g] + beta[1*g]*x + gamma[1*g]*(x-c[1])['+'] + gamma[2*g]*(x-c[2])['+']"

  pC <- ggplot(
    consensus_df,
    aes(
      rank,
      D_median
    )
  )

  pC <- add_region_background(
    pC,
    c1,
    c2,
    N
  )

  pC <- pC +
    geom_ribbon(
      aes(
        ymin = D_Q25,
        ymax = D_Q75
      ),
      fill = COL$ribbon,
      alpha = 0.25
    ) +
    geom_line(
      color = COL$nb,
      linewidth = 0.75
    ) +
    geom_line(
      aes(
        y = D_piecewise_fit_median
      ),
      color = COL$fit,
      linewidth = 1.25
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.45
    ) +
    geom_vline(
      xintercept = c1,
      color = COL$c1,
      linewidth = 0.9,
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = c2,
      color = COL$c2,
      linewidth = 0.9,
      linetype = "dashed"
    ) +
    labs(
      title = "C. Shared two-transition model and NB1/NB2 corroboration",
      subtitle = paste0(
        "Shared knots are estimated jointly across all 8 experimental groups: c1 = ",
        c1,
        ", c2 = ",
        c2
      ),
      x = "PC1 rank",
      y = "Cumulative divergence D(r)"
    ) +
    annotate(
      "label",
      x = max(
        1,
        round(
          0.025 *
            N
        )
      ),
      y = Inf,
      label = equation_box,
      parse = TRUE,
      hjust = 0,
      vjust = 1.2,
      size = 2.8,
      fill = "white"
    ) +
    annotate(
      "label",
      x = round(
        0.70 *
          N
      ),
      y = -Inf,
      label = nb_box,
      hjust = 0,
      vjust = -0.15,
      size = 2.8,
      fill = grDevices::adjustcolor(
        "white",
        alpha.f = 0.94
      )
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        face = "bold"
      )
    )

  # ---------------------------------------------------------------------------
  # One image, three panels.
  # ---------------------------------------------------------------------------
  grDevices::png(
    out_file,
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
            1.00,
            1.05,
            1.10
          ),
          "null"
        )
      )
    )
  )

  print(
    pA,
    vp = viewport(
      layout.pos.row = 1L,
      layout.pos.col = 1L
    )
  )

  print(
    pB,
    vp = viewport(
      layout.pos.row = 2L,
      layout.pos.col = 1L
    )
  )

  print(
    pC,
    vp = viewport(
      layout.pos.row = 3L,
      layout.pos.col = 1L
    )
  )

  dev.off()
}


# =============================================================================
# RUN ANALYSIS
# =============================================================================

count_mat <- read_count_matrix(
  COUNT_FILE,
  GROUP_PATTERNS
)

group_labels <- assign_groups(
  sample_names = colnames(
    count_mat
  ),
  group_patterns = GROUP_PATTERNS
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

# PC1 ranking matrix: same CPM-log1p convention as the original main analysis.
rank_matrix_all <- normalize_cpm_log1p(
  count_mat
)

# Global DESeq2 normalization for the design-aware sequencing variance.
deseq_norm <- normalize_deseq2(
  count_mat = count_mat,
  group_labels = group_labels
)

normalized_counts <- deseq_norm$normalized_counts

pooled <- compute_pooled_variance(
  normalized_counts = normalized_counts,
  group_labels = group_labels
)

pooled_variance <- pooled$variance

message(
  "Pooled within-group residual degrees of freedom: ",
  pooled$residual_df
)

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

for (
  g in groups
) {

  idx <- which(
    group_labels == g
  )

  message(
    "Computing PC1/NB mass geometry for ",
    g,
    "..."
  )

  group_curves[[
    g
  ]] <- compute_group_curves(
    group_name = g,
    raw_counts_arm = count_mat[
      ,
      idx,
      drop = FALSE
    ],
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

  write.csv(
    group_curves[[
      g
    ]],
    file.path(
      OUT_ROOT,
      paste0(
        "Table_RankedGeometry_",
        g,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

message(
  "Estimating the two shared divergence transitions jointly across all 8 groups..."
)

knot_fit <- fit_shared_knots(
  group_curves
)

C1 <- knot_fit$c1_rank

C2 <- knot_fit$c2_rank

message(
  "Shared REMAINDER boundary c1 = ",
  C1
)

message(
  "Shared LEADING-EDGE boundary c2 = ",
  C2
)

message(
  "Remainder size = ",
  C1 -
    1L
)

message(
  "Divergence interval size = ",
  C2 -
    C1 +
    1L
)

message(
  "Leading-edge size = ",
  nrow(
    group_curves[[
      1L
    ]]
  ) -
    C2
)

# NB1/NB2 characterization occurs AFTER boundary estimation.
region_results <- summarize_regions(
  group_curves = group_curves,
  c1 = C1,
  c2 = C2
)

consensus_df <- build_consensus_table(
  group_curves = group_curves,
  knot_fit = knot_fit
)

consensus_df$region <- assign_region(
  consensus_df$rank,
  C1,
  C2
)

boundary_table <- data.frame(
  remainder_boundary_c1 = C1,
  leading_edge_boundary_c2 = C2,
  remainder_size = C1 -
    1L,
  divergence_interval_size = C2 -
    C1 +
    1L,
  leading_edge_size = nrow(
    consensus_df
  ) -
    C2,
  joint_piecewise_SSE = knot_fit$SSE,
  stringsAsFactors = FALSE
)

write.csv(
  boundary_table,
  file.path(
    OUT_ROOT,
    "Table_Shared_Boundaries.csv"
  ),
  row.names = FALSE
)

write.csv(
  consensus_df,
  file.path(
    OUT_ROOT,
    "Table_Consensus_PC1_NB_Divergence.csv"
  ),
  row.names = FALSE
)

write.csv(
  region_results$exponent_by_group,
  file.path(
    OUT_ROOT,
    "Table_NB_Exponent_ByGroup_ByRegion.csv"
  ),
  row.names = FALSE
)

write.csv(
  region_results$exponent_summary,
  file.path(
    OUT_ROOT,
    "Table_NB_Exponent_Region_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  region_results$metrics_by_group,
  file.path(
    OUT_ROOT,
    "Table_NB_Diagnostics_ByGroup_ByRegion.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    sample = colnames(
      count_mat
    ),
    group = as.character(
      group_labels
    ),
    DESeq2_size_factor = as.numeric(
      deseq_norm$size_factors
    ),
    stringsAsFactors = FALSE
  ),
  file.path(
    OUT_ROOT,
    "Table_Samples_and_SizeFactors.csv"
  ),
  row.names = FALSE
)

make_main_figure(
  consensus_df = consensus_df,
  c1 = C1,
  c2 = C2,
  exponent_summary = region_results$exponent_summary,
  out_file = file.path(
    OUT_ROOT,
    "Figure_Main_PC1_NB_Cumulative_Divergence.png"
  )
)

message(
  "============================================================"
)

message(
  "ANALYSIS COMPLETE"
)

message(
  "PC1 rank = ascending absolute PC1 loading."
)

message(
  "c1 = ",
  C1,
  " | c2 = ",
  C2
)

message(
  "Main figure:"
)

message(
  file.path(
    OUT_ROOT,
    "Figure_Main_PC1_NB_Cumulative_Divergence.png"
  )
)

message(
  "Boundary table:"
)

message(
  file.path(
    OUT_ROOT,
    "Table_Shared_Boundaries.csv"
  )
)

message(
  "============================================================"
)
