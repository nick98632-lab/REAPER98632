#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# DIAGNOSTIC ANALYSIS ONLY:
# PC1 VARIANCE CONTRIBUTION vs NEGATIVE-BINOMIAL EXCESS VARIANCE
#
# IMPORTANT:
#   This script DOES NOT estimate a cutoff.
#   This script DOES NOT bootstrap.
#
# Its sole purpose is to calculate and display the mathematical object that
# should be inspected BEFORE any boundary rule is defined:
#
#   P_i   = d1^2 * v_i1^2 / (n - 1)
#   E_ig  = max(V_i - mu_ig, 0)
#
#   P*_g(r) = robustZ[ log(1 + P_g(r)) ]
#   E*_g(r) = robustZ[ log(1 + E_g(r)) ]
#
#   D_g(r) = E*_g(r) - P*_g(r)
#
#   D_cons(r) = median_g D_g(r)
#
# The consensus curve is smoothed by a GCV-selected smoothing spline and the
# first and second derivatives are calculated:
#
#   D'_cons(x)
#   D''_cons(x)
#
# where x = (rank - 1) / (N - 1) is normalized PC1 rank.
#
# PC1 RANKING IS EXPLICITLY:
#
#   ascending abs(PC1 loading)
#
# Thus:
#   low |PC1 loading|  -->  high |PC1 loading|
#   REMAINDER          -->  LEADING EDGE
#
# No boundary is selected in this script.
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <-
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <-
  "/root/REAPER98632/exports/PC1_NB_DIVERGENCE_DIAGNOSTIC"

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
# HELPERS
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
  ) / scale_value

  out
}


# =============================================================================
# READ COUNT MATRIX
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

  colnames(count_mat) <- colnames(count_df)

  storage.mode(count_mat) <- "numeric"

  count_mat[
    !is.finite(count_mat)
  ] <- 0

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
    feature_ids[blank] <- paste0(
      "__feature_row_",
      which(blank)
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
    stop("Fewer than two nonzero features remain.")
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
      "Unassigned samples: ",
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

# This matrix is used for PC1 so the original ranking convention is preserved:
# CPM -> log1p -> PCA/SVD -> ascending absolute PC1 loading.
normalize_cpm_log1p <- function(count_mat) {

  lib_size <- colSums(
    count_mat
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

  log1p(cpm)
}


# DESeq2 normalized counts are used to estimate the sequencing mean/variance
# relationship across the complete experimental design.
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
      "DESeq2 is required for this script. Install DESeq2 and rerun."
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
    DESeq2::estimateSizeFactors(
      dds
    ),
    error = function(e) {
      message(
        "Default DESeq2 size factors failed; using type='poscounts'."
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
# POOLED WITHIN-GROUP SEQUENCING VARIANCE
# =============================================================================

compute_pooled_within_group_variance <- function(
    normalized_counts,
    group_labels) {

  groups <- levels(
    factor(group_labels)
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
    stop("Insufficient pooled residual degrees of freedom.")
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
# PC1 FROM SVD
# =============================================================================

compute_pc1 <- function(
    rank_matrix_arm) {

  # rank_matrix_arm:
  #   features x samples
  #
  # X:
  #   samples x features
  X <- t(
    rank_matrix_arm
  )

  Xc <- sweep(
    X,
    2L,
    colMeans(X),
    "-"
  )

  # Efficient SVD geometry through the sample-space Gram matrix.
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

  # v1 = X^T u1 / d1
  v1 <- as.numeric(
    crossprod(
      Xc,
      u1
    )
  ) / d1

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

  # Feature-specific variance represented through PC1.
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
# ONE ARM: RANK BY abs(PC1 LOADING), THEN CALCULATE P*, E*, D
# =============================================================================

compute_group_divergence <- function(
    group_name,
    rank_matrix_arm,
    normalized_counts_arm,
    pooled_variance) {

  pc1 <- compute_pc1(
    rank_matrix_arm
  )

  # ---------------------------------------------------------------------------
  # THIS IS THE RANKING RULE:
  #
  #     ascending absolute PC1 loading
  #
  # It is deliberately written explicitly rather than ranking on P_i, even
  # though the two orders are mathematically identical within an arm.
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

  logP <- log1p(
    P_ranked
  )

  logE <- log1p(
    E_ranked
  )

  P_star <- robust_z(
    logP
  )

  E_star <- robust_z(
    logE
  )

  D <- E_star -
    P_star

  alpha_nb2 <- rep(
    0,
    length(mu_ranked)
  )

  positive_mu <- mu_ranked > 0

  alpha_nb2[
    positive_mu
  ] <- E_ranked[
    positive_mu
  ] /
    (
      mu_ranked[
        positive_mu
      ]^2
    )

  feature_df <- data.frame(
    group = group_name,
    rank = seq_along(
      rank_order
    ),
    feature_id = rownames(
      rank_matrix_arm
    )[rank_order],

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

    log1p_pc1_variance = logP,
    log1p_nb_excess_variance = logE,

    P_star = P_star,
    E_star = E_star,
    divergence = D,

    NB2_alpha = alpha_nb2,
    NB2 = log1p(
      E_ranked
    ),
    NB2_NB1 = log1p(
      E_ranked
    ) -
      log1p(
        mu_ranked
      ),
    alpha_mu = log1p(
      alpha_nb2 *
      mu_ranked
    ),

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
    stringsAsFactors = FALSE
  )

  list(
    features = feature_df,
    summary = summary_df
  )
}


# =============================================================================
# CONSENSUS ACROSS ALL 8 ARMS
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
      unique(N_values)
    ) != 1L
  ) {
    stop("All arms must contain the same number of ranked features.")
  }

  N <- N_values[
    1L
  ]

  P_mat <- do.call(
    cbind,
    lapply(
      group_curves,
      function(x) {
        x$P_star
      }
    )
  )

  E_mat <- do.call(
    cbind,
    lapply(
      group_curves,
      function(x) {
        x$E_star
      }
    )
  )

  D_mat <- do.call(
    cbind,
    lapply(
      group_curves,
      function(x) {
        x$divergence
      }
    )
  )

  colnames(P_mat) <- names(
    group_curves
  )

  colnames(E_mat) <- names(
    group_curves
  )

  colnames(D_mat) <- names(
    group_curves
  )

  out <- data.frame(
    rank = seq_len(
      N
    ),

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

    divergence_Q25 = apply(
      D_mat,
      1L,
      stats::quantile,
      probs = 0.25,
      na.rm = TRUE,
      names = FALSE
    ),

    divergence_Q75 = apply(
      D_mat,
      1L,
      stats::quantile,
      probs = 0.75,
      na.rm = TRUE,
      names = FALSE
    ),

    stringsAsFactors = FALSE
  )

  out
}


# =============================================================================
# GCV SPLINE + TRUE CALCULUS OF THE CONSENSUS DIVERGENCE
# =============================================================================

differentiate_consensus <- function(
    consensus_df) {

  N <- nrow(
    consensus_df
  )

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

  fit <- stats::smooth.spline(
    x = x,
    y = y,
    cv = FALSE
  )

  smooth_D <- as.numeric(
    stats::predict(
      fit,
      x = x,
      deriv = 0
    )$y
  )

  dD_dx <- as.numeric(
    stats::predict(
      fit,
      x = x,
      deriv = 1
    )$y
  )

  d2D_dx2 <- as.numeric(
    stats::predict(
      fit,
      x = x,
      deriv = 2
    )$y
  )

  # Also save derivatives with respect to integer rank r.
  dD_dr <- dD_dx /
    (
      N -
      1
    )

  d2D_dr2 <- d2D_dx2 /
    (
      N -
      1
    )^2

  data.frame(
    rank = consensus_df$rank,
    normalized_rank = x,
    smooth_divergence = smooth_D,
    dD_dx = dD_dx,
    d2D_dx2 = d2D_dx2,
    dD_dr = dD_dr,
    d2D_dr2 = d2D_dr2,
    stringsAsFactors = FALSE
  ) -> calculus_df

  list(
    calculus = calculus_df,
    spline = fit
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
        nrow = 4,
        ncol = 1
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
        layout.pos.col = 1
      )
    )
  }

  dev.off()
}


make_consensus_figure <- function(
    consensus_df,
    calculus_df,
    spline_fit,
    group_curves,
    out_file) {

  N <- nrow(
    consensus_df
  )

  all_D <- bind_rows(
    lapply(
      names(group_curves),
      function(g) {
        data.frame(
          group = g,
          rank = group_curves[[g]]$rank,
          divergence = group_curves[[g]]$divergence
        )
      }
    )
  )

  components_long <- consensus_df %>%
    select(
      rank,
      consensus_P_star,
      consensus_E_star
    ) %>%
    pivot_longer(
      cols = c(
        consensus_P_star,
        consensus_E_star
      ),
      names_to = "signal",
      values_to = "value"
    ) %>%
    mutate(
      signal = factor(
        signal,
        levels = c(
          "consensus_P_star",
          "consensus_E_star"
        ),
        labels = c(
          "P*: PC1 variance contribution",
          "E*: NB excess variance"
        )
      )
    )

  p1 <- ggplot(
    components_long,
    aes(
      rank,
      value,
      color = signal
    )
  ) +
    geom_line(
      linewidth = 0.9
    ) +
    labs(
      title = "A. PC1 variance contribution and NB excess variance on the same dimensionless scale",
      subtitle = "Features are ranked explicitly by ascending absolute PC1 loading",
      x = "PC1 rank: low |loading|  ->  high |loading|",
      y = "Robust standardized log signal",
      color = NULL
    ) +
    annotate(
      "label",
      x = round(
        N * 0.03
      ),
      y = Inf,
      label = "P[i] == d[1]^2*v[i*1]^2/(n-1)",
      parse = TRUE,
      hjust = 0,
      vjust = 1.2,
      size = 3.0,
      fill = "white"
    ) +
    annotate(
      "label",
      x = round(
        N * 0.03
      ),
      y = Inf,
      label = "E[i*g] == max(V[i]-mu[i*g],0)",
      parse = TRUE,
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

  p2 <- ggplot() +
    geom_line(
      data = all_D,
      aes(
        rank,
        divergence,
        group = group
      ),
      linewidth = 0.35,
      alpha = 0.28
    ) +
    geom_ribbon(
      data = consensus_df,
      aes(
        x = rank,
        ymin = divergence_Q25,
        ymax = divergence_Q75
      ),
      alpha = 0.18
    ) +
    geom_line(
      data = consensus_df,
      aes(
        rank,
        consensus_divergence
      ),
      linewidth = 0.8
    ) +
    geom_line(
      data = calculus_df,
      aes(
        rank,
        smooth_divergence
      ),
      linewidth = 1.05
    ) +
    labs(
      title = "B. PC1–NB divergence field across the ranked feature axis",
      subtitle = paste0(
        "Thin curves = 8 arms; ribbon = armwise IQR; thick smooth = GCV spline (effective df = ",
        round(
          spline_fit$df,
          1
        ),
        ")"
      ),
      x = "PC1 rank",
      y = "D(r)"
    ) +
    annotate(
      "label",
      x = round(
        N * 0.03
      ),
      y = Inf,
      label = "D[g](r) == E[g]^'*'(r)-P[g]^'*'(r)",
      parse = TRUE,
      hjust = 0,
      vjust = 1.2,
      size = 3.1,
      fill = "white"
    ) +
    annotate(
      "label",
      x = round(
        N * 0.03
      ),
      y = Inf,
      label = "D[cons](r) == median[g](D[g](r))",
      parse = TRUE,
      hjust = 0,
      vjust = 3.0,
      size = 3.1,
      fill = "white"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank()
    )

  p3 <- ggplot(
    calculus_df,
    aes(
      rank,
      dD_dx
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_line(
      linewidth = 0.9
    ) +
    labs(
      title = "C. First derivative of the smoothed divergence field",
      subtitle = "Positive values mean divergence is increasing as features move toward the high-PC1 leading edge",
      x = "PC1 rank",
      y = "dD/dx"
    ) +
    annotate(
      "label",
      x = round(
        N * 0.03
      ),
      y = Inf,
      label = "x == (r-1)/(N-1)",
      parse = TRUE,
      hjust = 0,
      vjust = 1.2,
      size = 3.1,
      fill = "white"
    ) +
    annotate(
      "label",
      x = round(
        N * 0.03
      ),
      y = Inf,
      label = "D^minute(x) == d*D/d*x",
      parse = TRUE,
      hjust = 0,
      vjust = 3.0,
      size = 3.1,
      fill = "white"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank()
    )

  p4 <- ggplot(
    calculus_df,
    aes(
      rank,
      d2D_dx2
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_line(
      linewidth = 0.9
    ) +
    labs(
      title = "D. Second derivative: curvature of PC1–NB divergence",
      subtitle = "Zero crossings and sustained curvature changes are visible here; no cutoff is selected",
      x = "PC1 rank",
      y = "d²D/dx²"
    ) +
    annotate(
      "label",
      x = round(
        N * 0.03
      ),
      y = Inf,
      label = "D^second(x) == d^2*D/d*x^2",
      parse = TRUE,
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

  x <- (
    feature_df$rank -
    1
  ) /
    (
      N -
      1
    )

  fit <- stats::smooth.spline(
    x = x,
    y = feature_df$divergence,
    cv = FALSE
  )

  smooth_D <- as.numeric(
    stats::predict(
      fit,
      x = x,
      deriv = 0
    )$y
  )

  d1 <- as.numeric(
    stats::predict(
      fit,
      x = x,
      deriv = 1
    )$y
  )

  d2 <- as.numeric(
    stats::predict(
      fit,
      x = x,
      deriv = 2
    )$y
  )

  components_long <- feature_df %>%
    select(
      rank,
      P_star,
      E_star
    ) %>%
    pivot_longer(
      cols = c(
        P_star,
        E_star
      ),
      names_to = "signal",
      values_to = "value"
    )

  p1 <- ggplot(
    components_long,
    aes(
      rank,
      value,
      color = signal
    )
  ) +
    geom_line(
      linewidth = 0.75
    ) +
    labs(
      title = paste0(
        group_name,
        ": P* versus E*"
      ),
      subtitle = "Rank = ascending absolute PC1 loading",
      x = "PC1 rank",
      y = "Robust standardized log signal",
      color = NULL
    ) +
    annotate(
      "label",
      x = round(
        N * 0.03
      ),
      y = Inf,
      label = "P[i] == d[1]^2*v[i*1]^2/(n-1)",
      parse = TRUE,
      hjust = 0,
      vjust = 1.2,
      size = 3.0,
      fill = "white"
    ) +
    annotate(
      "label",
      x = round(
        N * 0.03
      ),
      y = Inf,
      label = "E[i*g] == max(V[i]-mu[i*g],0)",
      parse = TRUE,
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

  divergence_df <- data.frame(
    rank = feature_df$rank,
    divergence = feature_df$divergence,
    smooth_divergence = smooth_D
  )

  p2 <- ggplot(
    divergence_df,
    aes(
      rank,
      divergence
    )
  ) +
    geom_line(
      linewidth = 0.35,
      alpha = 0.45
    ) +
    geom_line(
      aes(
        y = smooth_divergence
      ),
      linewidth = 1.0
    ) +
    labs(
      title = paste0(
        group_name,
        ": divergence"
      ),
      subtitle = paste0(
        "GCV spline effective df = ",
        round(
          fit$df,
          1
        )
      ),
      x = "PC1 rank",
      y = "D(r)"
    ) +
    annotate(
      "label",
      x = round(
        N * 0.03
      ),
      y = Inf,
      label = "D[g](r) == E[g]^'*'(r)-P[g]^'*'(r)",
      parse = TRUE,
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

  derivative_df <- bind_rows(
    data.frame(
      rank = feature_df$rank,
      derivative = "D'(x)",
      value = d1
    ),
    data.frame(
      rank = feature_df$rank,
      derivative = "D''(x)",
      value = d2
    )
  )

  p3 <- ggplot(
    derivative_df,
    aes(
      rank,
      value
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_line(
      linewidth = 0.75
    ) +
    facet_wrap(
      ~ derivative,
      ncol = 1,
      scales = "free_y"
    ) +
    labs(
      title = paste0(
        group_name,
        ": calculus of divergence"
      ),
      subtitle = "Derivatives are with respect to normalized rank x = (r-1)/(N-1)",
      x = "PC1 rank",
      y = NULL
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor = element_blank()
    )

  png(
    out_file,
    width = 15,
    height = 11.5,
    units = "in",
    res = 260,
    bg = "white"
  )

  grid.newpage()

  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 3,
        ncol = 1
      )
    )
  )

  print(
    p1,
    vp = viewport(
      layout.pos.row = 1,
      layout.pos.col = 1
    )
  )

  print(
    p2,
    vp = viewport(
      layout.pos.row = 2,
      layout.pos.col = 1
    )
  )

  print(
    p3,
    vp = viewport(
      layout.pos.row = 3,
      layout.pos.col = 1
    )
  )

  dev.off()

  data.frame(
    rank = feature_df$rank,
    normalized_rank = x,
    smooth_divergence = smooth_D,
    dD_dx = d1,
    d2D_dx2 = d2,
    stringsAsFactors = FALSE
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
  colnames(count_mat),
  GROUP_PATTERNS
)

message(
  "Count matrix: ",
  nrow(count_mat),
  " features x ",
  ncol(count_mat),
  " samples"
)

message(
  "Groups: ",
  paste(
    levels(group_labels),
    collapse = ", "
  )
)


# -----------------------------------------------------------------------------
# PC1 matrix: original CPM-log1p convention.
# -----------------------------------------------------------------------------

rank_matrix_all <- normalize_cpm_log1p(
  count_mat
)


# -----------------------------------------------------------------------------
# Sequencing variance matrix: DESeq2-normalized counts.
# -----------------------------------------------------------------------------

deseq_norm <- normalize_deseq2(
  count_mat,
  group_labels
)

normalized_counts <- deseq_norm$normalized_counts


# -----------------------------------------------------------------------------
# Stable variance estimated from all groups after subtracting group means.
# -----------------------------------------------------------------------------

pooled <- compute_pooled_within_group_variance(
  normalized_counts,
  group_labels
)

pooled_variance <- pooled$variance

message(
  "Pooled within-group residual df = ",
  pooled$residual_df
)


# -----------------------------------------------------------------------------
# Calculate each of the 8 PC1-ranked divergence trajectories.
# -----------------------------------------------------------------------------

groups <- levels(
  group_labels
)

group_curves <- vector(
  "list",
  length(groups)
)

names(group_curves) <- groups

pc1_summary_list <- vector(
  "list",
  length(groups)
)

names(pc1_summary_list) <- groups


for (
  g in groups
) {

  idx <- which(
    group_labels == g
  )

  one <- compute_group_divergence(
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

  group_curves[[g]] <- one$features

  pc1_summary_list[[g]] <- one$summary

  write.csv(
    one$features,
    file.path(
      OUT_ROOT,
      paste0(
        "Table_Divergence_",
        g,
        ".csv"
      )
    ),
    row.names = FALSE
  )

  group_calculus <- make_group_figure(
    feature_df = one$features,
    group_name = g,
    out_file = file.path(
      OUT_ROOT,
      paste0(
        "Figure_Diagnostic_",
        g,
        ".png"
      )
    )
  )

  write.csv(
    group_calculus,
    file.path(
      OUT_ROOT,
      paste0(
        "Table_Calculus_",
        g,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}


pc1_summary <- bind_rows(
  pc1_summary_list
)

write.csv(
  pc1_summary,
  file.path(
    OUT_ROOT,
    "Table_PC1_Summary.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# Build one consensus divergence trajectory across all 8 arms.
# -----------------------------------------------------------------------------

consensus_df <- compute_consensus(
  group_curves
)

calculus_obj <- differentiate_consensus(
  consensus_df
)

consensus_output <- consensus_df %>%
  left_join(
    calculus_obj$calculus,
    by = "rank"
  )

write.csv(
  consensus_output,
  file.path(
    OUT_ROOT,
    "Table_Consensus_Divergence_and_Calculus.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# Save sample normalization information.
# -----------------------------------------------------------------------------

write.csv(
  data.frame(
    sample = colnames(count_mat),
    group = as.character(group_labels),
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


# -----------------------------------------------------------------------------
# Main diagnostic mathematical figure.
# -----------------------------------------------------------------------------

make_consensus_figure(
  consensus_df = consensus_df,
  calculus_df = calculus_obj$calculus,
  spline_fit = calculus_obj$spline,
  group_curves = group_curves,
  out_file = file.path(
    OUT_ROOT,
    "Figure_Consensus_PC1_NB_Divergence_DIAGNOSTIC.png"
  )
)


# =============================================================================
# FINAL CONSOLE OUTPUT
# =============================================================================

message(
  "============================================================"
)

message(
  "DIAGNOSTIC PC1-NB DIVERGENCE ANALYSIS COMPLETE"
)

message(
  "PC1 ordering used: ascending abs(PC1 loading)"
)

message(
  "No cutoff was estimated."
)

message(
  "No bootstrap was performed."
)

message(
  "Inspect:"
)

message(
  "  Figure_Consensus_PC1_NB_Divergence_DIAGNOSTIC.png"
)

message(
  "  Table_Consensus_Divergence_and_Calculus.csv"
)

message(
  "Output directory: ",
  OUT_ROOT
)

message(
  "============================================================"
)
