#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# FINAL MANUSCRIPT ANALYSIS
# PC1-RANKED VARIANCE GEOMETRY, CUMULATIVE PC1-NB DIVERGENCE,
# AND TERMINAL LEADING-EDGE CUTOFF
# =============================================================================
#
# MATHEMATICAL METHODS
# --------------------
#
# 1. PC1 ranking
#
# For each experimental group, raw sequencing counts are library-size
# normalized to counts per million (CPM), transformed as log(1 + CPM), and PCA
# is performed across samples.  Feature i is ranked by the absolute value of
# its PC1 loading:
#
#       rank_i = rank(|v_i1|), ascending.
#
# Thus the x-axis runs from low PC1 contribution (remainder) to high PC1
# contribution (leading edge).  The feature-specific variance represented by
# PC1 is
#
#       P_i = lambda_1 * v_i1^2
#
# where lambda_1 is the PC1 eigenvalue.  Within one PCA, ranking by |v_i1|,
# v_i1^2, or P_i is identical.
#
#
# 2. Two explicitly different empirical variance measurements
#
# A. RAW-COUNT EMPIRICAL VARIANCE -- used for the familiar variance-geometry
#    curve and for the terminal cutoff:
#
#       s_raw,ig^2 = Var_j(Y_ij | group g)
#
#    where Y_ij is the raw read count.  This within-group sample variance is
#    reordered by the PC1 rank and displayed as log(1 + s_raw^2).
#
# B. POOLED WITHIN-GROUP EMPIRICAL VARIANCE OF DESeq2-NORMALIZED COUNTS --
#    used for the PC1-vs-NB cumulative-divergence analysis:
#
#                         sum_g sum_{j in g} (y_ij - ybar_ig)^2
#       V_i,pool =         -------------------------------------
#                                  sum_g (n_g - 1)
#
#    where y_ij is the DESeq2 size-factor-normalized count.  This is an
#    empirical pooled within-group sample variance; it is NOT a fitted DESeq2
#    dispersion parameter.  For group g,
#
#       mu_ig = mean_j(y_ij | g)
#       E_ig  = max(V_i,pool - mu_ig, 0)
#
#    E_ig is the empirical excess-over-Poisson variance signal.
#
#
# 3. Cumulative PC1-NB variance-mass divergence
#
# Within each group, PC1 variance contribution and excess count variance are
# converted into probability masses over the same PC1-ranked feature axis:
#
#       p_g(r) = P_g(r) / sum_j P_g(j)
#       q_g(r) = E_g(r) / sum_j E_g(j)
#
# Their cumulative masses are
#
#       F_P,g(r) = sum_{j <= r} p_g(j)
#       F_E,g(r) = sum_{j <= r} q_g(j)
#
# and the cumulative divergence is
#
#       D_g(r) = F_E,g(r) - F_P,g(r).
#
# Because both masses sum to one, D_g(N) = 0.  Its discrete first difference is
#
#       Delta D_g(r) = q_g(r) - p_g(r),
#
# so changes in D describe where excess sequencing variance and PC1-associated
# variance accumulate at different rates along the PC1-ranked axis.
#
#
# 4. Shared divergence-regime boundaries
#
# All eight group-specific D_g curves are fit jointly with a continuous
# two-knot linear-spline model:
#
#       D_g(x) = beta_0g + beta_1g*x
#                + gamma_1g*(x-c1)_+
#                + gamma_2g*(x-c2)_+
#
# where x = (rank-1)/(N-1) and (z)_+ = max(z,0).  Each group has its own
# coefficients, while c1 and c2 are shared across all groups.  The shared knots
# minimize the summed squared residual error over all eight curves.
#
#       rank < c1          : remainder
#       c1 <= rank <= c2   : divergence interval
#       rank > c2          : data-derived leading-edge regime
#
#
# 5. Anchor-Terminal range inside the data-derived leading-edge regime
#
# After the shared divergence model estimates c2, the older raw-count variance
# geometry is applied ONLY inside the data-derived leading edge (rank > c2).
# For each group:
#
#   a. fit the original smoothing spline (spar = 0.60) to
#          log(1 + within-group raw-count empirical variance)
#      across the complete PC1-ranked axis;
#   b. calculate the spline first and second derivatives;
#   c. identify second-derivative sign-change zero crossings;
#   d. within rank > c2, form consecutive zero-crossing intervals;
#   e. retain intervals with a positive net rise in the smoothed raw-count
#      variance curve and a positive median first derivative;
#   f. select the RIGHTMOST retained interval.
#
# Its left curvature crossing is Anchor and its right curvature crossing is
# Terminal:
#
#       Anchor_g   = left zero crossing of the rightmost sustained rising
#                    curvature interval inside rank > c2
#
#       Terminal_g = right zero crossing of that same interval.
#
# Therefore c2 defines where the leading-edge search begins, while the original
# derivative geometry identifies the local Anchor-Terminal range inside it.
#
# The manuscript's historical 5,000-feature cutoff is retained only as a
# graphical/table reference:
#
#       Ref = N - 5000 + 1.
#
# Ref has ZERO influence on c1, c2, Anchor, Terminal, regional NB scaling, or
# any other fitted quantity.
#
#
# 6. NB1/NB2 corroboration after boundaries are defined
#
# Using the pooled normalized-count excess variance above, the regional
# mean-variance relationship is summarized by
#
#       E = alpha * mu^p
#       log(E) = log(alpha) + p*log(mu).
#
# The fitted regional slope p is an empirical scaling exponent:
#
#       p near 1 : more NB1-like
#       p near 2 : more NB2-like.
#
# This regional scaling is calculated after the geometric boundaries are
# defined; it does not determine c1, c2, Ref, or Terminal.
#
#
# OUTPUTS
# -------
# Figures are written individually to OUT_ROOT/Figures and are also collected
# into OUT_ROOT/Figures_All.zip.
#
#   Figure_Overall.png
#   Figure_RT0_ZT6.png
#   Figure_RT2_ZT8.png
#   Figure_RT4_ZT10.png
#   Figure_RT8_ZT14.png
#
# Only two compact result tables are written:
#
#   Table_Key_Results.csv
#   Table_Timepoints.csv
#
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <-
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <-
  "/root/REAPER98632/exports/pc1_nb_final"

FIG_DIR <- file.path(
  OUT_ROOT,
  "Figures"
)

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

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
  RT0_ZT6  = c(control = "RT0", treatment = "ZT6"),
  RT2_ZT8  = c(control = "RT2", treatment = "ZT8"),
  RT4_ZT10 = c(control = "RT4", treatment = "ZT10"),
  RT8_ZT14 = c(control = "RT8", treatment = "ZT14")
)

PAPER_LEADING_EDGE_SIZE <- 5000L

# Exact smoothing level used by the older terminal-cutoff geometry.
TERMINAL_SPLINE_SPAR <- 0.60

# Display-only smoothing.  These do not determine any cutoff.
DISPLAY_VAR_SPAR <- 0.72
DISPLAY_MASS_SPAR <- 0.72
DISPLAY_D_SPAR <- 0.72

PNG_WIDTH_IN <- 15
PNG_DPI <- 360


# =============================================================================
# COLORS
# =============================================================================

COL <- list(
  empirical = "#117A65",
  control = "#386CB0",
  treatment = "#159D91",
  pc1 = "#386CB0",
  nb = "#159D91",
  divergence = "#222222",
  fit = "#111111",
  remainder = "#DCEFF2",
  interval = "#F5E8C8",
  leading = "#DDF2EA",
  c1 = "#2166AC",
  c2 = "#1B7837",
  ref = "#E69F00",
  terminal = "#D95F02"
)


# =============================================================================
# DATA INPUT
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

  n_bad <- sum(
    !is.finite(count_mat)
  )

  if (n_bad > 0L) {
    message(
      "Replacing ",
      n_bad,
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
      "Unassigned samples: ",
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
# NORMALIZATION AND EMPIRICAL VARIANCE
# =============================================================================

normalize_cpm_log1p <- function(
    count_mat_arm) {

  lib_size <- colSums(
    count_mat_arm,
    na.rm = TRUE
  )

  lib_size[
    !is.finite(lib_size) |
    lib_size <= 0
  ] <- 1

  cpm <- sweep(
    count_mat_arm,
    2L,
    lib_size / 1e6,
    "/"
  )

  log1p(
    cpm
  )
}


normalize_deseq2_global <- function(
    count_mat,
    group_labels) {

  if (
    !requireNamespace(
      "DESeq2",
      quietly = TRUE
    )
  ) {
    stop("DESeq2 is required for the pooled normalized-count variance.")
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
# PC1 RANKING
# =============================================================================

compute_pc1_rank <- function(
    rank_matrix_arm) {

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
    !is.finite(loading)
  ] <- 0

  abs_loading <- abs(
    loading
  )

  rank_order <- order(
    abs_loading,
    decreasing = FALSE
  )

  lambda1 <- pca$sdev[
    1L
  ]^2

  P <- lambda1 *
    loading^2

  list(
    loading = loading,
    abs_loading = abs_loading,
    rank_order = rank_order,
    lambda1 = lambda1,
    pc1_variance_contribution = P
  )
}


# =============================================================================
# ORIGINAL TERMINAL CUTOFF: RAW-COUNT VARIANCE SPLINE SECOND DERIVATIVE
# =============================================================================

find_d2_zero_crossings <- function(
    x,
    d2) {

  ok <- (
    is.finite(x) &
    is.finite(d2)
  )

  x <- x[
    ok
  ]

  d2 <- d2[
    ok
  ]

  if (length(x) < 2L) {
    return(
      numeric(0)
    )
  }

  out <- numeric(0)

  for (
    i in seq_len(
      length(x) -
        1L
    )
  ) {

    a <- d2[
      i
    ]

    b <- d2[
      i +
        1L
    ]

    if (
      a == 0 ||
      b == 0
    ) {
      next
    }

    if (
      (
        a < 0 &&
        b > 0
      ) ||
      (
        a > 0 &&
        b < 0
      )
    ) {

      frac <- abs(
        a
      ) /
        (
          abs(a) +
            abs(b)
        )

      out <- c(
        out,
        x[
          i
        ] +
          frac *
          (
            x[
              i +
                1L
            ] -
              x[
                i
              ]
          )
      )
    }
  }

  sort(
    unique(
      out
    )
  )
}


compute_raw_variance_geometry <- function(
    raw_counts_arm,
    rank_order) {

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

  raw_var_ranked <- raw_var[
    rank_order
  ]

  rank <- seq_along(
    rank_order
  )

  log_var <- log1p(
    raw_var_ranked
  )

  # Analysis spline: exact original smoothing level used for derivative
  # geometry.  This spline determines Anchor and Terminal after c2 is known.
  analysis_spline <- stats::smooth.spline(
    x = rank,
    y = log_var,
    spar = TERMINAL_SPLINE_SPAR
  )

  dense_x <- seq(
    min(rank),
    max(rank),
    length.out = max(
      5000L,
      length(rank) *
        4L
    )
  )

  dense_y <- as.numeric(
    stats::predict(
      analysis_spline,
      x = dense_x,
      deriv = 0
    )$y
  )

  dense_d1 <- as.numeric(
    stats::predict(
      analysis_spline,
      x = dense_x,
      deriv = 1
    )$y
  )

  dense_d2 <- as.numeric(
    stats::predict(
      analysis_spline,
      x = dense_x,
      deriv = 2
    )$y
  )

  dense_df <- data.frame(
    dense_rank = dense_x,
    dense_smooth_log1p_raw_variance = dense_y,
    dense_d1 = dense_d1,
    dense_d2 = dense_d2,
    stringsAsFactors = FALSE
  )

  crossings <- find_d2_zero_crossings(
    dense_df$dense_rank,
    dense_df$dense_d2
  )

  # Cleaner display spline only; this does not determine Anchor or Terminal.
  display_spline <- stats::smooth.spline(
    x = rank,
    y = log_var,
    spar = DISPLAY_VAR_SPAR
  )

  display_y <- as.numeric(
    stats::predict(
      display_spline,
      x = rank,
      deriv = 0
    )$y
  )

  curve <- data.frame(
    rank = rank,
    raw_empirical_variance = raw_var_ranked,
    log1p_raw_empirical_variance = log_var,
    display_log1p_raw_empirical_variance = display_y,
    stringsAsFactors = FALSE
  )

  list(
    curve = curve,
    dense = dense_df,
    crossings = crossings
  )
}


select_anchor_terminal_in_leading_edge <- function(
    dense_df,
    crossings,
    c2,
    total_n) {

  # Ref / the historical 5,000-feature cutoff is intentionally absent from
  # this function.  The search domain is determined only by the independently
  # estimated shared leading-edge boundary c2.

  z <- sort(
    unique(
      crossings[
        is.finite(crossings) &
        crossings >
          c2 &
        crossings <
          total_n
      ]
    )
  )

  if (length(z) < 2L) {
    stop(
      "Fewer than two raw-variance d2 zero crossings occur inside the ",
      "data-derived leading edge (rank > c2 = ",
      c2,
      ")."
    )
  }

  candidates <- vector(
    "list",
    length(z) -
      1L
  )

  for (
    k in seq_len(
      length(z) -
        1L
    )
  ) {

    a <- z[
      k
    ]

    b <- z[
      k +
        1L
    ]

    idx <- which(
      dense_df$dense_rank >=
        a &
      dense_df$dense_rank <=
        b
    )

    if (length(idx) < 2L) {
      next
    }

    y_a <- approx(
      x = dense_df$dense_rank,
      y = dense_df$dense_smooth_log1p_raw_variance,
      xout = a,
      rule = 2
    )$y

    y_b <- approx(
      x = dense_df$dense_rank,
      y = dense_df$dense_smooth_log1p_raw_variance,
      xout = b,
      rule = 2
    )$y

    delta_y <- as.numeric(
      y_b -
        y_a
    )

    median_d1 <- median(
      dense_df$dense_d1[
        idx
      ],
      na.rm = TRUE
    )

    positive_slope_fraction <- mean(
      dense_df$dense_d1[
        idx
      ] >
        0,
      na.rm = TRUE
    )

    candidates[[
      k
    ]] <- data.frame(
      anchor_crossing = a,
      terminal_crossing = b,
      delta_log_variance = delta_y,
      median_d1 = median_d1,
      positive_slope_fraction = positive_slope_fraction,
      stringsAsFactors = FALSE
    )
  }

  candidate_df <- bind_rows(
    candidates
  )

  valid <- candidate_df %>%
    filter(
      is.finite(
        delta_log_variance
      ),
      is.finite(
        median_d1
      ),
      delta_log_variance >
        0,
      median_d1 >
        0
    ) %>%
    arrange(
      desc(
        anchor_crossing
      )
    )

  if (nrow(valid) == 0L) {
    stop(
      "No sustained rising d2-bounded raw-variance interval was found ",
      "inside the data-derived leading edge (rank > c2 = ",
      c2,
      ")."
    )
  }

  selected <- valid[
    1L,
    ,
    drop = FALSE
  ]

  anchor <- as.integer(
    round(
      selected$anchor_crossing[
        1L
      ]
    )
  )

  terminal <- as.integer(
    round(
      selected$terminal_crossing[
        1L
      ]
    )
  )

  anchor <- max(
    c2 +
      1L,
    min(
      total_n -
        1L,
      anchor
    )
  )

  terminal <- max(
    anchor +
      1L,
    min(
      total_n,
      terminal
    )
  )

  if (
    anchor <=
      c2 ||
    terminal <=
      anchor
  ) {
    stop("Invalid Anchor-Terminal interval after integer conversion.")
  }

  list(
    anchor = anchor,
    terminal = terminal,
    delta_log_variance = selected$delta_log_variance[
      1L
    ],
    median_d1 = selected$median_d1[
      1L
    ],
    positive_slope_fraction = selected$positive_slope_fraction[
      1L
    ],
    candidates = candidate_df
  )
}


# =============================================================================
# GROUP-SPECIFIC CUMULATIVE DIVERGENCE
# =============================================================================

smooth_nonnegative_mass <- function(
    rank,
    mass,
    spar = DISPLAY_MASS_SPAR) {

  fit <- stats::smooth.spline(
    x = rank,
    y = mass,
    spar = spar
  )

  y <- as.numeric(
    stats::predict(
      fit,
      x = rank,
      deriv = 0
    )$y
  )

  y[
    !is.finite(y)
  ] <- 0

  y <- pmax(
    y,
    0
  )

  total <- sum(
    y
  )

  if (
    !is.finite(total) ||
    total <= 0
  ) {
    stop("Display mass could not be normalized.")
  }

  y /
    total
}


smooth_divergence_for_display <- function(
    rank,
    D,
    spar = DISPLAY_D_SPAR) {

  fit <- stats::smooth.spline(
    x = rank,
    y = D,
    spar = spar
  )

  y <- as.numeric(
    stats::predict(
      fit,
      x = rank,
      deriv = 0
    )$y
  )

  # Remove the fitted endpoint line so the display curve respects D(1)≈D(N)=0.
  endpoint_line <- seq(
    y[
      1L
    ],
    y[
      length(y)
    ],
    length.out = length(y)
  )

  y -
    endpoint_line
}


compute_group_analysis <- function(
    group_name,
    raw_counts_arm,
    normalized_counts_arm,
    pooled_variance) {

  rank_matrix <- normalize_cpm_log1p(
    raw_counts_arm
  )

  pc1 <- compute_pc1_rank(
    rank_matrix
  )

  rank_order <- pc1$rank_order

  geometry <- compute_raw_variance_geometry(
    raw_counts_arm = raw_counts_arm,
    rank_order = rank_order
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

  P_ranked <- pc1$pc1_variance_contribution[
    rank_order
  ]

  mu_ranked <- mu[
    rank_order
  ]

  V_ranked <- pooled_variance[
    rank_order
  ]

  E_ranked <- E[
    rank_order
  ]

  P_total <- sum(
    P_ranked
  )

  E_total <- sum(
    E_ranked
  )

  if (
    !is.finite(P_total) ||
    P_total <= 0
  ) {
    stop("PC1 variance mass is undefined for group ", group_name)
  }

  if (
    !is.finite(E_total) ||
    E_total <= 0
  ) {
    stop("NB excess-variance mass is undefined for group ", group_name)
  }

  p_mass <- P_ranked /
    P_total

  q_mass <- E_ranked /
    E_total

  F_P <- cumsum(
    p_mass
  )

  F_E <- cumsum(
    q_mass
  )

  D <- F_E -
    F_P

  rank <- seq_along(
    rank_order
  )

  display_p <- smooth_nonnegative_mass(
    rank,
    p_mass
  )

  display_q <- smooth_nonnegative_mass(
    rank,
    q_mass
  )

  display_F_P <- cumsum(
    display_p
  )

  display_F_E <- cumsum(
    display_q
  )

  display_D <- smooth_divergence_for_display(
    rank,
    display_F_E -
      display_F_P
  )

  df <- geometry$curve %>%
    mutate(
      group = group_name,
      feature_id = rownames(
        raw_counts_arm
      )[
        rank_order
      ],
      pc1_loading = pc1$loading[
        rank_order
      ],
      abs_pc1_loading = pc1$abs_loading[
        rank_order
      ],
      pc1_eigenvalue = pc1$lambda1,
      pc1_variance_contribution = P_ranked,
      pc1_variance_mass = p_mass,
      pooled_normalized_variance = V_ranked,
      normalized_group_mean = mu_ranked,
      nb_excess_variance = E_ranked,
      nb_excess_variance_mass = q_mass,
      cumulative_pc1_mass = F_P,
      cumulative_nb_mass = F_E,
      cumulative_divergence = D,
      local_mass_difference = q_mass -
        p_mass,
      display_F_P = display_F_P,
      display_F_E = display_F_E,
      display_D = display_D
    )

  list(
    data = df,
    dense_raw_variance_geometry = geometry$dense,
    raw_variance_crossings = geometry$crossings,
    anchor = NA_integer_,
    terminal = NA_integer_,
    anchor_terminal_delta_log_variance = NA_real_,
    anchor_terminal_median_d1 = NA_real_,
    anchor_terminal_positive_slope_fraction = NA_real_
  )
}


# =============================================================================
# SHARED TWO-KNOT MODEL
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
      c1 <=
      min_gap
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

  coef <- tryCatch(
    qr.coef(
      qr(
        X
      ),
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

  resid <- D_mat -
    X %*%
    coef

  sum(
    resid^2
  )
}


fit_shared_knots <- function(
    group_results) {

  groups <- names(
    group_results
  )

  N_values <- vapply(
    group_results,
    function(z) {
      nrow(
        z$data
      )
    },
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
      group_results,
      function(z) {
        z$data$cumulative_divergence
      }
    )
  )

  colnames(
    D_full
  ) <- groups

  opt_n <- min(
    5000L,
    N
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

  coarse <- lapply(
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

  values <- vapply(
    coarse,
    function(z) {
      z$value
    },
    numeric(1)
  )

  best <- coarse[[
    which.min(
      values
    )
  ]]

  refined <- stats::optim(
    par = best$par,
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
    stop("Shared-knot optimization failed.")
  }

  c1_rank <- as.integer(
    round(
      1 +
        refined$par[
          1L
        ] *
        (
          N -
            1
        )
    )
  )

  c2_rank <- as.integer(
    round(
      1 +
        refined$par[
          2L
        ] *
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

  c1_x <- (
    c1_rank -
      1
  ) /
    (
      N -
        1
    )

  c2_x <- (
    c2_rank -
      1
  ) /
    (
      N -
        1
    )

  X <- piecewise_basis(
    x_full,
    c1_x,
    c2_x
  )

  coef <- qr.coef(
    qr(
      X
    ),
    D_full
  )

  fitted <- X %*%
    coef

  list(
    c1 = c1_rank,
    c2 = c2_rank,
    x = x_full,
    fitted = fitted,
    SSE = sum(
      (
        D_full -
          fitted
      )^2
    ),
    groups = groups
  )
}


# =============================================================================
# REGIONAL NB SCALING
# =============================================================================

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
      NA_real_
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

  unname(
    stats::coef(
      fit
    )[
      2L
    ]
  )
}


get_region_p <- function(
    df,
    keep) {

  estimate_nb_exponent(
    mu = df$normalized_group_mean[
      keep
    ],
    E = df$nb_excess_variance[
      keep
    ]
  )
}


# =============================================================================
# RESULT TABLES
# =============================================================================

build_timepoint_table <- function(
    group_results,
    comparisons,
    c1,
    c2,
    reference_rank) {

  rows <- list()

  for (
    comparison_name in names(
      comparisons
    )
  ) {

    mapping <- comparisons[[
      comparison_name
    ]]

    for (
      arm in names(
        mapping
      )
    ) {

      g <- unname(
        mapping[[
          arm
        ]]
      )

      z <- group_results[[
        g
      ]]

      df <- z$data

      anchor <- z$anchor
      terminal <- z$terminal

      rows[[
        length(rows) +
          1L
      ]] <- data.frame(
        comparison = comparison_name,
        arm = arm,
        group = g,
        anchor_raw_variance_d2 = anchor,
        terminal_raw_variance_d2 = terminal,
        anchor_terminal_width = terminal -
          anchor +
          1L,
        terminal_tail_n = nrow(df) -
          terminal +
          1L,
        anchor_terminal_delta_log_variance =
          z$anchor_terminal_delta_log_variance,
        anchor_terminal_median_d1 =
          z$anchor_terminal_median_d1,
        p_remainder = get_region_p(
          df,
          df$rank <
            c1
        ),
        p_model_leading_edge = get_region_p(
          df,
          df$rank >
            c2
        ),
        p_paper_5000 = get_region_p(
          df,
          df$rank >=
            reference_rank
        ),
        p_anchor_terminal = get_region_p(
          df,
          df$rank >=
            anchor &
          df$rank <=
            terminal
        ),
        p_terminal_tail = get_region_p(
          df,
          df$rank >=
            terminal
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  bind_rows(
    rows
  )
}


safe_summary <- function(
    x,
    fun = median) {

  x <- x[
    is.finite(
      x
    )
  ]

  if (length(x) == 0L) {
    return(
      NA_real_
    )
  }

  fun(
    x
  )
}


build_key_table <- function(
    timepoint_table,
    N,
    c1,
    c2,
    reference_rank,
    shared_sse) {

  data.frame(
    n_features = N,
    shared_c1 = c1,
    shared_c2 = c2,
    remainder_n = c1 -
      1L,
    divergence_interval_n = c2 -
      c1 +
      1L,
    model_leading_edge_n = N -
      c2,

    anchor_median = safe_summary(
      timepoint_table$anchor_raw_variance_d2
    ),
    anchor_min = safe_summary(
      timepoint_table$anchor_raw_variance_d2,
      min
    ),
    anchor_max = safe_summary(
      timepoint_table$anchor_raw_variance_d2,
      max
    ),

    terminal_median = safe_summary(
      timepoint_table$terminal_raw_variance_d2
    ),
    terminal_min = safe_summary(
      timepoint_table$terminal_raw_variance_d2,
      min
    ),
    terminal_max = safe_summary(
      timepoint_table$terminal_raw_variance_d2,
      max
    ),

    paper_ref = reference_rank,
    paper_leading_edge_n = PAPER_LEADING_EDGE_SIZE,
    paper_ref_used_in_analysis = FALSE,

    p_remainder_median = safe_summary(
      timepoint_table$p_remainder
    ),
    p_model_leading_edge_median = safe_summary(
      timepoint_table$p_model_leading_edge
    ),
    p_paper_5000_median = safe_summary(
      timepoint_table$p_paper_5000
    ),
    p_anchor_terminal_median = safe_summary(
      timepoint_table$p_anchor_terminal
    ),
    p_terminal_tail_median = safe_summary(
      timepoint_table$p_terminal_tail
    ),

    shared_model_SSE = shared_sse,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# FIGURE DATA
# =============================================================================

rankwise_median <- function(
    group_results,
    column) {

  mat <- do.call(
    cbind,
    lapply(
      group_results,
      function(z) {
        z$data[[
          column
        ]]
      }
    )
  )

  apply(
    mat,
    1L,
    median,
    na.rm = TRUE
  )
}


build_overall_figure_data <- function(
    group_results,
    knot_fit) {

  N <- nrow(
    group_results[[
      1L
    ]]$data
  )

  data.frame(
    rank = seq_len(
      N
    ),
    raw_variance = rankwise_median(
      group_results,
      "display_log1p_raw_empirical_variance"
    ),
    F_P = rankwise_median(
      group_results,
      "display_F_P"
    ),
    F_E = rankwise_median(
      group_results,
      "display_F_E"
    ),
    D = rankwise_median(
      group_results,
      "display_D"
    ),
    D_fit = apply(
      knot_fit$fitted,
      1L,
      median,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# FIGURE STYLE
# =============================================================================

theme_manuscript <- function(
    base_size = 11.5) {

  theme_classic(
    base_size = base_size
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size +
          1.3
      ),
      plot.subtitle = element_text(
        size = base_size -
          0.4,
        margin = margin(
          b = 5
        )
      ),
      axis.title = element_text(
        face = "bold"
      ),
      axis.text = element_text(
        color = "#333333"
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.border = element_rect(
        color = "#B7B7B7",
        fill = NA,
        linewidth = 0.45
      ),
      panel.grid = element_blank(),
      plot.margin = margin(
        7,
        9,
        7,
        9
      )
    )
}


add_regions <- function(
    p,
    c1,
    c2,
    reference_rank,
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
      alpha = 0.46
    ) +
    annotate(
      "rect",
      xmin = c2,
      xmax = N,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$leading,
      alpha = 0.55
    ) +
    annotate(
      "rect",
      xmin = reference_rank,
      xmax = N,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$ref,
      alpha = 0.055
    ) +
    geom_vline(
      xintercept = c1,
      color = COL$c1,
      linetype = "dashed",
      linewidth = 0.70
    ) +
    geom_vline(
      xintercept = c2,
      color = COL$c2,
      linetype = "dashed",
      linewidth = 0.78
    ) +
    geom_vline(
      xintercept = reference_rank,
      color = COL$ref,
      linetype = "dotdash",
      linewidth = 0.78
    )
}


save_panels <- function(
    plots,
    path,
    height_in) {

  grDevices::png(
    path,
    width = PNG_WIDTH_IN,
    height = height_in,
    units = "in",
    res = PNG_DPI,
    bg = "white"
  )

  grid.newpage()

  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = length(
          plots
        ),
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
      plots[[
        i
      ]],
      vp = viewport(
        layout.pos.row = i,
        layout.pos.col = 1L
      )
    )
  }

  dev.off()
}


# =============================================================================
# OVERALL FIGURE
# =============================================================================

make_overall_figure <- function(
    overall_df,
    key_table,
    c1,
    c2,
    reference_rank,
    anchor_median,
    terminal_median,
    out_file) {

  N <- nrow(
    overall_df
  )

  # Panel A: raw-count empirical variance geometry.
  pA <- ggplot(
    overall_df,
    aes(
      rank,
      raw_variance
    )
  )

  pA <- add_regions(
    pA,
    c1,
    c2,
    reference_rank,
    N
  )

  pA <- pA +
    geom_vline(
      xintercept = anchor_median,
      color = COL$divergence,
      linetype = "solid",
      linewidth = 0.78
    ) +
    geom_vline(
      xintercept = terminal_median,
      color = COL$terminal,
      linetype = "dotted",
      linewidth = 0.85
    ) +
    geom_line(
      color = COL$empirical,
      linewidth = 1.25,
      lineend = "round"
    ) +
    labs(
      title = "A. PC1-ranked raw-count variance geometry",
      subtitle = paste0(
        "Shared c1 = ",
        c1,
        "; c2 = ",
        c2,
        "; median Anchor = ",
        round(
          anchor_median
        ),
        "; median Terminal = ",
        round(
          terminal_median
        ),
        "; Ref = ",
        reference_rank,
        " (reference only)"
      ),
      x = "PC1 rank: low |loading|  ->  high |loading|",
      y = "Smoothed log(1 + within-group raw-count variance)"
    ) +
    annotate(
      "label",
      x = round(
        0.025 *
          N
      ),
      y = Inf,
      label = "PC1 rank: log(1+CPM);  Pᵢ=λ₁vᵢ₁²\nVariance shown: within-group raw-count sample variance",
      hjust = 0,
      vjust = 1.15,
      size = 2.9,
      fill = "white"
    ) +
    theme_manuscript()

  # Panel B: cumulative normalized variance masses.
  cumulative_long <- overall_df %>%
    select(
      rank,
      F_P,
      F_E
    ) %>%
    pivot_longer(
      cols = c(
        F_P,
        F_E
      ),
      names_to = "curve",
      values_to = "value"
    ) %>%
    mutate(
      curve = factor(
        curve,
        levels = c(
          "F_P",
          "F_E"
        ),
        labels = c(
          "PC1 variance mass F_P(r)",
          "NB excess-variance mass F_E(r)"
        )
      )
    )

  pB <- ggplot()

  pB <- add_regions(
    pB,
    c1,
    c2,
    reference_rank,
    N
  )

  pB <- pB +
    geom_vline(
      xintercept = anchor_median,
      color = COL$divergence,
      linetype = "solid",
      linewidth = 0.78
    ) +
    geom_vline(
      xintercept = terminal_median,
      color = COL$terminal,
      linetype = "dotted",
      linewidth = 0.85
    ) +
    geom_line(
      data = cumulative_long,
      aes(
        rank,
        value,
        color = curve
      ),
      linewidth = 1.05,
      lineend = "round"
    ) +
    geom_line(
      data = overall_df,
      aes(
        rank,
        D
      ),
      color = COL$divergence,
      linewidth = 1.15,
      lineend = "round"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.4
    ) +
    scale_color_manual(
      values = c(
        "PC1 variance mass F_P(r)" = COL$pc1,
        "NB excess-variance mass F_E(r)" = COL$nb
      )
    ) +
    labs(
      title = "B. Cumulative PC1-NB variance-mass divergence",
      subtitle = "NB excess variance uses pooled within-group empirical variance of DESeq2-normalized counts",
      x = "PC1 rank",
      y = "Cumulative mass / D(r)",
      color = NULL
    ) +
    annotate(
      "label",
      x = round(
        0.025 *
          N
      ),
      y = 0.97,
      label = "p(r)=P(r)/ΣP;  q(r)=E(r)/ΣE;  D(r)=F_E(r)−F_P(r)\nE=max(V_pool−μ_g,0);  V_pool = pooled within-group normalized-count variance",
      hjust = 0,
      vjust = 1,
      size = 2.8,
      fill = "white"
    ) +
    theme_manuscript()

  p_summary <- paste0(
    "Empirical NB scaling:  E = α μ^p\n",
    "Remainder p = ",
    sprintf(
      "%.2f",
      key_table$p_remainder_median
    ),
    "\nModel leading edge p = ",
    sprintf(
      "%.2f",
      key_table$p_model_leading_edge_median
    ),
    "\nPaper 5,000 p = ",
    sprintf(
      "%.2f",
      key_table$p_paper_5000_median
    ),
    "\nTerminal tail p = ",
    sprintf(
      "%.2f",
      key_table$p_terminal_tail_median
    )
  )

  # Panel C: fitted shared divergence transitions.
  pC <- ggplot(
    overall_df,
    aes(
      rank,
      D
    )
  )

  pC <- add_regions(
    pC,
    c1,
    c2,
    reference_rank,
    N
  )

  pC <- pC +
    geom_vline(
      xintercept = anchor_median,
      color = COL$divergence,
      linetype = "solid",
      linewidth = 0.78
    ) +
    geom_vline(
      xintercept = terminal_median,
      color = COL$terminal,
      linetype = "dotted",
      linewidth = 0.85
    ) +
    geom_line(
      color = COL$nb,
      linewidth = 0.80,
      alpha = 0.55,
      lineend = "round"
    ) +
    geom_line(
      aes(
        y = D_fit
      ),
      color = COL$fit,
      linewidth = 1.35,
      lineend = "round"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.4
    ) +
    labs(
      title = "C. Shared divergence transitions and local Anchor-Terminal range",
      subtitle = "c2 defines the leading-edge search domain; raw-count variance curvature then identifies Anchor and Terminal; Ref is display-only",
      x = "PC1 rank",
      y = "Cumulative divergence D(r)"
    ) +
    annotate(
      "label",
      x = round(
        0.025 *
          N
      ),
      y = Inf,
      label = "D_g(x)=β₀g+β₁g x+γ₁g(x−c₁)₊+γ₂g(x−c₂)₊\nWithin rank>c2: Anchor-Terminal = rightmost sustained rising interval bounded by raw-variance d²/dr² zero crossings",
      hjust = 0,
      vjust = 1.15,
      size = 2.75,
      fill = "white"
    ) +
    annotate(
      "label",
      x = round(
        0.69 *
          N
      ),
      y = -Inf,
      label = p_summary,
      hjust = 0,
      vjust = -0.08,
      size = 2.75,
      fill = "white"
    ) +
    theme_manuscript()

  save_panels(
    plots = list(
      pA,
      pB,
      pC
    ),
    path = out_file,
    height_in = 11.5
  )
}


# =============================================================================
# TIME-POINT FIGURE
# =============================================================================

make_timepoint_figure <- function(
    comparison_name,
    mapping,
    group_results,
    timepoint_table,
    c1,
    c2,
    reference_rank,
    out_file) {

  control_group <- unname(
    mapping[[
      "control"
    ]]
  )

  treatment_group <- unname(
    mapping[[
      "treatment"
    ]]
  )

  control <- group_results[[
    control_group
  ]]$data %>%
    mutate(
      arm = paste0(
        "Control (",
        control_group,
        ")"
      )
    )

  treatment <- group_results[[
    treatment_group
  ]]$data %>%
    mutate(
      arm = paste0(
        "Treatment (",
        treatment_group,
        ")"
      )
    )

  plot_df <- bind_rows(
    control,
    treatment
  )

  N <- nrow(
    control
  )

  arm_colors <- stats::setNames(
    c(
      COL$control,
      COL$treatment
    ),
    c(
      paste0(
        "Control (",
        control_group,
        ")"
      ),
      paste0(
        "Treatment (",
        treatment_group,
        ")"
      )
    )
  )

  term_rows <- timepoint_table %>%
    filter(
      comparison ==
        comparison_name
    ) %>%
    mutate(
      arm_label = ifelse(
        arm ==
          "control",
        paste0(
          "Control (",
          group,
          ")"
        ),
        paste0(
          "Treatment (",
          group,
          ")"
        )
      )
    )

  cutoff_df <- bind_rows(
    data.frame(
      arm = term_rows$arm_label,
      event = "Anchor",
      rank = term_rows$anchor_raw_variance_d2,
      stringsAsFactors = FALSE
    ),
    data.frame(
      arm = term_rows$arm_label,
      event = "Terminal",
      rank = term_rows$terminal_raw_variance_d2,
      stringsAsFactors = FALSE
    )
  )

  # Panel A: raw-count empirical variance, close to the original geometry figure.
  pA <- ggplot(
    plot_df,
    aes(
      rank,
      display_log1p_raw_empirical_variance,
      color = arm
    )
  )

  pA <- add_regions(
    pA,
    c1,
    c2,
    reference_rank,
    N
  )

  pA <- pA +
    geom_vline(
      data = cutoff_df %>% filter(event == "Anchor"),
      aes(
        xintercept = rank,
        color = arm
      ),
      linetype = "solid",
      linewidth = 0.70,
      alpha = 0.85,
      show.legend = FALSE
    ) +
    geom_vline(
      data = cutoff_df %>% filter(event == "Terminal"),
      aes(
        xintercept = rank,
        color = arm
      ),
      linetype = "dotted",
      linewidth = 0.92,
      show.legend = FALSE
    ) +
    geom_line(
      linewidth = 1.18,
      lineend = "round"
    ) +
    scale_color_manual(
      values = arm_colors
    ) +
    labs(
      title = paste0(
        "A. ",
        comparison_name,
        ": raw-count variance geometry"
      ),
      subtitle = paste0(
        "c2 = ",
        c2,
        "; ",
        control_group,
        " A/T = ",
        term_rows$anchor_raw_variance_d2[
          term_rows$group == control_group
        ],
        "/",
        term_rows$terminal_raw_variance_d2[
          term_rows$group == control_group
        ],
        "; ",
        treatment_group,
        " A/T = ",
        term_rows$anchor_raw_variance_d2[
          term_rows$group == treatment_group
        ],
        "/",
        term_rows$terminal_raw_variance_d2[
          term_rows$group == treatment_group
        ],
        "; Ref = ",
        reference_rank,
        " (reference only)"
      ),
      x = "PC1 rank: low |loading|  ->  high |loading|",
      y = "Smoothed log(1 + within-group raw-count variance)",
      color = NULL
    ) +
    annotate(
      "label",
      x = round(
        0.025 *
          N
      ),
      y = Inf,
      label = "Rank: |PC1 loading| from log(1+CPM)\nVariance shown: within-arm raw-count sample variance",
      hjust = 0,
      vjust = 1.15,
      size = 2.85,
      fill = "white"
    ) +
    theme_manuscript()

  # Panel B: cumulative divergence.
  pB <- ggplot(
    plot_df,
    aes(
      rank,
      display_D,
      color = arm
    )
  )

  pB <- add_regions(
    pB,
    c1,
    c2,
    reference_rank,
    N
  )

  pB <- pB +
    geom_vline(
      data = cutoff_df %>% filter(event == "Anchor"),
      aes(
        xintercept = rank,
        color = arm
      ),
      linetype = "solid",
      linewidth = 0.70,
      alpha = 0.85,
      show.legend = FALSE
    ) +
    geom_vline(
      data = cutoff_df %>% filter(event == "Terminal"),
      aes(
        xintercept = rank,
        color = arm
      ),
      linetype = "dotted",
      linewidth = 0.92,
      show.legend = FALSE
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.4
    ) +
    geom_line(
      linewidth = 1.18,
      lineend = "round"
    ) +
    scale_color_manual(
      values = arm_colors
    ) +
    labs(
      title = paste0(
        "B. ",
        comparison_name,
        ": cumulative PC1-NB divergence"
      ),
      subtitle = "D(r) compares PC1 variance mass with pooled normalized-count excess-variance mass",
      x = "PC1 rank",
      y = "D(r) = F_E(r) − F_P(r)",
      color = NULL
    ) +
    annotate(
      "label",
      x = round(
        0.025 *
          N
      ),
      y = Inf,
      label = "PC1 variance: Pᵢ=λ₁vᵢ₁²\nNB excess variance: E=max(V_pool−μ_g,0), V_pool from DESeq2-normalized counts",
      hjust = 0,
      vjust = 1.15,
      size = 2.8,
      fill = "white"
    ) +
    theme_manuscript()

  # Panel C: NB scaling exponent by key region.
  p_rows <- term_rows %>%
    select(
      arm_label,
      p_remainder,
      p_model_leading_edge,
      p_anchor_terminal,
      p_paper_5000,
      p_terminal_tail
    ) %>%
    pivot_longer(
      cols = c(
        p_remainder,
        p_model_leading_edge,
        p_anchor_terminal,
        p_paper_5000,
        p_terminal_tail
      ),
      names_to = "region",
      values_to = "p"
    ) %>%
    mutate(
      region = factor(
        region,
        levels = c(
          "p_remainder",
          "p_model_leading_edge",
          "p_anchor_terminal",
          "p_paper_5000",
          "p_terminal_tail"
        ),
        labels = c(
          "Remainder",
          "Model leading edge",
          "Anchor-Terminal",
          "Paper 5,000 (reference)",
          "Terminal tail"
        )
      )
    )

  pC <- ggplot(
    p_rows,
    aes(
      x = p,
      y = region,
      color = arm_label
    )
  ) +
    geom_vline(
      xintercept = 1,
      color = "#777777",
      linetype = "dashed",
      linewidth = 0.65
    ) +
    geom_vline(
      xintercept = 2,
      color = "#777777",
      linetype = "dotted",
      linewidth = 0.75
    ) +
    geom_point(
      size = 3.3
    ) +
    scale_color_manual(
      values = arm_colors
    ) +
    labs(
      title = paste0(
        "C. ",
        comparison_name,
        ": NB mean-variance scaling"
      ),
      subtitle = "E = α μ^p using pooled within-group empirical variance of DESeq2-normalized counts",
      x = "Empirical exponent p   (1 = NB1-like; 2 = NB2-like)",
      y = NULL,
      color = NULL
    ) +
    theme_manuscript()

  save_panels(
    plots = list(
      pA,
      pB,
      pC
    ),
    path = out_file,
    height_in = 10.8
  )
}


# =============================================================================
# RUN
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
  "Count matrix: ",
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
  "Groups: ",
  paste(
    levels(
      group_labels
    ),
    collapse = ", "
  )
)

N <- nrow(
  count_mat
)

if (
  PAPER_LEADING_EDGE_SIZE >=
    N
) {
  stop("PAPER_LEADING_EDGE_SIZE must be smaller than the feature count.")
}

REFERENCE_RANK <- N -
  PAPER_LEADING_EDGE_SIZE +
  1L

deseq <- normalize_deseq2_global(
  count_mat = count_mat,
  group_labels = group_labels
)

normalized_counts <- deseq$normalized_counts

pooled <- compute_pooled_within_group_variance(
  normalized_counts = normalized_counts,
  group_labels = group_labels
)

message(
  "Pooled within-group residual degrees of freedom: ",
  pooled$residual_df
)

group_results <- vector(
  "list",
  length(
    levels(
      group_labels
    )
  )
)

names(
  group_results
) <- levels(
  group_labels
)

for (
  g in levels(
    group_labels
  )
) {

  idx <- which(
    group_labels ==
      g
  )

  if (length(idx) < 2L) {
    stop("Not enough samples in group ", g)
  }

  message(
    "Analyzing ",
    g,
    "..."
  )

  group_results[[
    g
  ]] <- compute_group_analysis(
    group_name = g,
    raw_counts_arm = count_mat[
      ,
      idx,
      drop = FALSE
    ],
    normalized_counts_arm = normalized_counts[
      ,
      idx,
      drop = FALSE
    ],
    pooled_variance = pooled$variance
  )
}

message(
  "Fitting shared two-knot cumulative-divergence model..."
)

knot_fit <- fit_shared_knots(
  group_results
)

C1 <- knot_fit$c1
C2 <- knot_fit$c2

message(
  "Selecting arm-specific Anchor-Terminal ranges inside the data-derived leading edge (rank > c2)..."
)

for (
  g in names(
    group_results
  )
) {

  at <- select_anchor_terminal_in_leading_edge(
    dense_df = group_results[[g]]$dense_raw_variance_geometry,
    crossings = group_results[[g]]$raw_variance_crossings,
    c2 = C2,
    total_n = N
  )

  group_results[[g]]$anchor <- at$anchor
  group_results[[g]]$terminal <- at$terminal
  group_results[[g]]$anchor_terminal_delta_log_variance <-
    at$delta_log_variance
  group_results[[g]]$anchor_terminal_median_d1 <-
    at$median_d1
  group_results[[g]]$anchor_terminal_positive_slope_fraction <-
    at$positive_slope_fraction
}

# Historical 5,000-feature Ref is calculated only after all fitted quantities
# above have been obtained.  It is used only in tables/figures as a reference.
if (
  REFERENCE_RANK <=
    C2
) {
  warning(
    "The historical 5,000-feature reference begins before or at shared c2; ",
    "this does not affect the fitted analysis."
  )
}

timepoint_table <- build_timepoint_table(
  group_results = group_results,
  comparisons = COMPARISONS,
  c1 = C1,
  c2 = C2,
  reference_rank = REFERENCE_RANK
)

key_table <- build_key_table(
  timepoint_table = timepoint_table,
  N = N,
  c1 = C1,
  c2 = C2,
  reference_rank = REFERENCE_RANK,
  shared_sse = knot_fit$SSE
)

write.csv(
  key_table,
  file.path(
    OUT_ROOT,
    "Table_Key_Results.csv"
  ),
  row.names = FALSE
)

write.csv(
  timepoint_table,
  file.path(
    OUT_ROOT,
    "Table_Timepoints.csv"
  ),
  row.names = FALSE
)

overall_df <- build_overall_figure_data(
  group_results = group_results,
  knot_fit = knot_fit
)

figure_paths <- character(0)

overall_path <- file.path(
  FIG_DIR,
  "Figure_Overall.png"
)

make_overall_figure(
  overall_df = overall_df,
  key_table = key_table,
  c1 = C1,
  c2 = C2,
  reference_rank = REFERENCE_RANK,
  anchor_median = key_table$anchor_median[
    1L
  ],
  terminal_median = key_table$terminal_median[
    1L
  ],
  out_file = overall_path
)

figure_paths <- c(
  figure_paths,
  overall_path
)

for (
  comparison_name in names(
    COMPARISONS
  )
) {

  fig_path <- file.path(
    FIG_DIR,
    paste0(
      "Figure_",
      comparison_name,
      ".png"
    )
  )

  make_timepoint_figure(
    comparison_name = comparison_name,
    mapping = COMPARISONS[[
      comparison_name
    ]],
    group_results = group_results,
    timepoint_table = timepoint_table,
    c1 = C1,
    c2 = C2,
    reference_rank = REFERENCE_RANK,
    out_file = fig_path
  )

  figure_paths <- c(
    figure_paths,
    fig_path
  )
}


# =============================================================================
# ZIP ALL FIGURES, WHILE KEEPING EACH PNG AVAILABLE INDIVIDUALLY
# =============================================================================

ZIP_PATH <- file.path(
  OUT_ROOT,
  "Figures_All.zip"
)

if (file.exists(ZIP_PATH)) {
  unlink(
    ZIP_PATH
  )
}

old_wd <- getwd()

zip_ok <- FALSE

tryCatch(
  {
    setwd(
      FIG_DIR
    )

    utils::zip(
      zipfile = ZIP_PATH,
      files = basename(
        figure_paths
      )
    )

    zip_ok <- file.exists(
      ZIP_PATH
    )
  },
  finally = {
    setwd(
      old_wd
    )
  }
)

if (!zip_ok) {
  warning(
    "Figure PNGs were created, but Figures_All.zip was not created."
  )
}


# =============================================================================
# FINAL CONSOLE SUMMARY
# =============================================================================

message(
  "============================================================"
)

message(
  "FINAL ANALYSIS COMPLETE"
)

message(
  "PC1 ranking: ascending absolute PC1 loading."
)

message(
  "Shared c1 = ",
  C1
)

message(
  "Shared c2 = ",
  C2
)

message(
  "Historical 5,000-feature Ref (REFERENCE ONLY; not used in fitting) = ",
  REFERENCE_RANK
)

message(
  "Arm-specific Anchor-Terminal ranges: ",
  paste0(
    timepoint_table$group,
    "=",
    timepoint_table$anchor_raw_variance_d2,
    "-",
    timepoint_table$terminal_raw_variance_d2,
    collapse = "; "
  )
)

message(
  "Figures available individually in: ",
  FIG_DIR
)

message(
  "Figure zip: ",
  ZIP_PATH
)

message(
  "Tables: Table_Key_Results.csv; Table_Timepoints.csv"
)

message(
  "============================================================"
)
