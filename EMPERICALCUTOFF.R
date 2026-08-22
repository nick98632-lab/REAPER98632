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
# PC1-NB REGIME GEOMETRY + c1-CONSTRAINED EIGENVECTOR SPLITTING
# =============================================================================
#
# PURPOSE
# -------
# This script performs three primary analyses:
#
#   A. Raw-count variance geometry along independently ranked PC1 axes.
#   B. Cumulative PC1-versus-NB excess-variance mass divergence to estimate
#      shared transition knots c1 and c2.
#   C. Data-driven top-k optimization for the control/treatment eigenvector
#      split, using c2 as the own-arm selection boundary and c1 as the
#      opposite-arm hard contamination boundary.
#
# The historical top-5,000 value is NEVER used to estimate c1, c2, Anchor,
# Terminal, the cutoff change-point, or any other fitted quantity. It is shown
# only as a prespecified reference.
#
#
# RANKING
# -------
# Each experimental arm is ranked independently by ASCENDING absolute PC1
# loading:
#
#     rank_order <- order(abs_loading, decreasing = FALSE)
#
# Therefore larger rank = stronger absolute PC1 loading.
#
#
# RAW-COUNT VARIANCE GEOMETRY
# ---------------------------
# For arm g:
#
#     s_raw,g^2(r) = Var(raw counts at rank r within arm g)
#     y_g(r)       = log(1 + s_raw,g^2(r))
#
# A smoothing spline (spar=0.60) is fit to y_g(r). ONE second derivative
# y_g''(r) is calculated.
#
# After c2:
#
#     Anchor   = first sign-change zero crossing of y_g''(r)
#     Terminal = second successive sign-change zero crossing of y_g''(r)
#
#
# CUMULATIVE PC1-NB VARIANCE-MASS DIVERGENCE
# ------------------------------------------
# Feature-level PC1 variance contribution:
#
#     P_i = lambda_1 * v_i1^2
#
# Pooled within-group empirical variance of DESeq2-normalized counts:
#
#                  sum_g sum_{j in g}(y_ij - ybar_ig)^2
#     V_pool,i =   ------------------------------------
#                           sum_g(n_g - 1)
#
# Arm-specific mean and excess-over-Poisson variance:
#
#     mu_ig = mean normalized count in arm g
#     E_ig  = max(V_pool,i - mu_ig, 0)
#
# Rank-wise masses:
#
#     p_g(r) = P_g(r) / sum P_g
#     q_g(r) = E_g(r) / sum E_g
#
# Cumulative divergence:
#
#     F_P,g(r) = cumsum[p_g(r)]
#     F_E,g(r) = cumsum[q_g(r)]
#     D_g(r)   = F_E,g(r) - F_P,g(r)
#
# A shared two-knot continuous linear-spline model across all eight arms
# estimates c1 and c2:
#
#     r < c1          : Remainder
#     c1 <= r <= c2   : Divergence interval
#     r > c2          : Leading-edge regime
#
#
# TOP-k EIGENVECTOR SPLITTING
# ---------------------------
# For each control/treatment pair and candidate top-k depth:
#
#     S_C(k) = control top-k
#     S_T(k) = treatment top-k
#
#     Joint              = S_C(k) intersection S_T(k)
#     Disjoint control   = S_C(k) \ S_T(k)
#     Disjoint treatment = S_T(k) \ S_C(k)
#
# Candidate k is restricted to:
#
#     k <= N - c2
#
# so every selected site originates from the c2-defined leading-edge regime
# of the arm that selected it.
#
# For a DISJOINT site, the opposite-arm rank is then classified as:
#
#     r_opposite > c2
#         opposite-arm Leading Edge        -> retained
#
#     c1 <= r_opposite <= c2
#         opposite-arm Divergence interval -> retained / permissible
#
#     r_opposite < c1
#         opposite-arm Remainder           -> hard cross-regime contamination
#
# Thus crossing c2 in the opposite arm is NOT automatically penalized.
# Only crossing c1 into the opposite-arm remainder is treated as contamination.
#
# For each candidate k:
#
#     G(k) =
#         Joint
#         + Disjoint with opposite-arm Leading Edge
#         + Disjoint with opposite-arm Divergence
#
#     R(k) =
#         Disjoint sites with opposite-arm rank < c1
#
# The cutoff problem is:
#
#     maximize G(k)
#     minimize R(k)
#
# subject to k <= N-c2.
#
# The number of disjoint sites landing in the opposite-arm Divergence interval
# is exported as a separate diagnostic and is NOT treated as contamination.
#
#
# DATA-DRIVEN CUTOFF CHANGE-POINT
# -------------------------------
# Non-dominated points in [R(k), G(k)] define the Pareto frontier.
#
# On that frontier the script fits:
#
#     G(R) = beta0 + beta1*R + gamma*(R - tau)+
#
# with:
#
#     slope_before = beta1
#     slope_after  = beta1 + gamma
#
# The preferred breakpoint has:
#
#     slope_after < slope_before
#
# and is selected by minimum BIC among shape-valid segmented models.
#
# A single linear G-versus-R model is also fit. Positive:
#
#     Delta_BIC = BIC_linear - BIC_segmented
#
# supports the segmented change-point description.
#
# If the segmented breakpoint is not supportable, the code uses a conservative
# fallback rather than inventing a gain/contamination weight:
#
#     largest k with zero opposite-arm remainder crossings, if available.
#
# The four pair-specific scans are also pooled:
#
#     G_total(k) = sum_m G_m(k)
#     R_total(k) = sum_m R_m(k)
#
# to obtain one common manuscript cutoff.
#
#
# POST-BOUNDARY NB MEAN-VARIANCE SCALING
# --------------------------------------
# This is NOT part of boundary or cutoff selection.
#
# After all boundaries are fixed, the script optionally estimates:
#
#     E = alpha * mu^p
#
# using lm(log(E) ~ log(mu)) for positive mu,E.
#
# Its only purpose is orthogonal corroboration: it asks whether the empirical
# count mean-variance relationship changes across the Remainder, Divergence,
# and Leading-edge regimes. Because it does not determine c1, c2, or k*, it is
# removed from the main figure family and written only as a supplemental
# diagnostic figure/table by default.
#
#
# PRIMARY FIGURES
# ---------------
#   Figure_Overall.png
#   Figure_RT0_ZT6.png
#   Figure_RT2_ZT8.png
#   Figure_RT4_ZT10.png
#   Figure_RT8_ZT14.png
#
# Each primary figure contains ONLY:
#   A. Raw-count variance geometry
#   B. Cumulative PC1-NB variance-mass divergence
#   C. c1-based top-k cutoff optimization
#
# OPTIONAL SUPPLEMENT
# -------------------
#   Figure_Supplement_NB_Scaling.png
#   Table_NB_Scaling.csv
#
# TABLES / SITE SETS
# ------------------
#   Table_Key_Results.csv
#   Table_Timepoints.csv
#   Table_Cutoff_Optimization.csv
#   Selected_LeadingEdge_Sites.csv
#   Excluded_Remainder_Crossing_Sites.csv
#   Figures_All.zip
#
# =============================================================================

# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <- "/root/REAPER98632/exports/pc1_nb_cutoff_c1_final"
FIG_DIR  <- file.path(OUT_ROOT, "Figures")

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

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

PAPER_REFERENCE_K <- 5000L

# Original raw-variance derivative smoothing.
ANALYSIS_SPLINE_SPAR <- 0.60

# Display-only smoothing. These values never determine boundaries/cutoffs.
DISPLAY_VAR_SPAR <- 0.72
DISPLAY_D_SPAR   <- 0.72

PNG_WIDTH_IN <- 15
PNG_DPI      <- 360

# Keep NB scaling as supplemental corroboration, not a primary panel.
INCLUDE_NB_SCALING_SUPPLEMENT <- TRUE


# =============================================================================
# COLORS
# =============================================================================

COL <- list(
  control = "#386CB0",
  treatment = "#159D91",
  raw = "#117A65",
  divergence = "#262626",
  fit = "#000000",

  remainder = "#DCEFF2",
  interval = "#F5E8C8",
  leading = "#DDF2EA",

  c1 = "#2166AC",
  c2 = "#1B7837",
  selected = "#7B3294",
  paper = "#E69F00",

  candidate = "#A6A6A6",
  pareto = "#6A3D9A",
  segmented = "#222222",
  divergence_disjoint = "#D9A441",
  remainder_cross = "#B2182B",

  anchor_control = "#386CB0",
  terminal_control = "#386CB0",
  anchor_treatment = "#159D91",
  terminal_treatment = "#159D91"
)


# =============================================================================
# BASIC HELPERS
# =============================================================================

theme_manuscript <- function(base_size = 11.5) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.1),
      plot.subtitle = element_text(size = base_size - 0.3, margin = margin(b = 4)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#333333"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1.2),
      legend.box = "horizontal",
      panel.border = element_rect(color = "#B7B7B7", fill = NA, linewidth = 0.45),
      panel.grid = element_blank(),
      plot.margin = margin(6, 8, 6, 8)
    )
}

save_panels <- function(plots, path, height_in) {
  grDevices::png(
    filename = path,
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
        nrow = length(plots),
        ncol = 1L
      )
    )
  )

  for (i in seq_along(plots)) {
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

safe_summary <- function(x, fun = median) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  fun(x)
}

rank_cutoff_from_k <- function(N, k) {
  if (!is.finite(k) || k < 1L || k > N) return(NA_integer_)
  as.integer(N - k + 1L)
}


# =============================================================================
# DATA INPUT
# =============================================================================

read_count_matrix <- function(path, group_patterns) {
  if (!file.exists(path)) {
    stop("Count file does not exist: ", path)
  }

  raw_df <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(raw_df) < 1L || ncol(raw_df) < 2L) {
    stop("Count file is empty or malformed.")
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
    stop("Column 1 matched a sample pattern; column 1 must contain feature IDs.")
  }

  count_df <- raw_df[, sample_idx, drop = FALSE]

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

  n_bad <- sum(!is.finite(count_mat))

  if (n_bad > 0L) {
    message("Replacing ", n_bad, " non-finite count entries with 0.")
    count_mat[!is.finite(count_mat)] <- 0
  }

  count_mat <- pmax(count_mat, 0)

  feature_ids <- trimws(as.character(raw_df[[1L]]))

  blank <- is.na(feature_ids) | feature_ids == ""

  if (any(blank)) {
    feature_ids[blank] <- paste0("__feature_row_", which(blank))
  }

  feature_ids <- make.unique(feature_ids, sep = "__dup_")
  rownames(count_mat) <- feature_ids

  keep <- rowSums(count_mat) > 0
  count_mat <- count_mat[keep, , drop = FALSE]

  if (nrow(count_mat) < 2L) {
    stop("Fewer than two nonzero features remain after filtering.")
  }

  count_mat
}

assign_groups <- function(sample_names, group_patterns) {
  assigned <- rep(NA_character_, length(sample_names))

  for (group_name in names(group_patterns)) {
    idx <- grep(group_patterns[[group_name]], sample_names)

    if (length(idx) > 0L && any(!is.na(assigned[idx]))) {
      stop("At least one sample matched more than one group pattern.")
    }

    assigned[idx] <- group_name
  }

  if (any(is.na(assigned))) {
    stop(
      "Unassigned samples: ",
      paste(sample_names[is.na(assigned)], collapse = ", ")
    )
  }

  factor(
    assigned,
    levels = names(group_patterns)
  )
}


# =============================================================================
# NORMALIZATION AND VARIANCE
# =============================================================================

normalize_cpm_log1p <- function(count_mat_arm) {
  lib_size <- colSums(count_mat_arm, na.rm = TRUE)
  lib_size[!is.finite(lib_size) | lib_size <= 0] <- 1

  cpm <- sweep(
    count_mat_arm,
    2L,
    lib_size / 1e6,
    "/"
  )

  log1p(cpm)
}

normalize_deseq2_global <- function(count_mat, group_labels) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required.")
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
      message("Default DESeq2 size factors failed; using type='poscounts'.")
      DESeq2::estimateSizeFactors(dds, type = "poscounts")
    }
  )

  list(
    normalized_counts = DESeq2::counts(dds, normalized = TRUE),
    size_factors = DESeq2::sizeFactors(dds)
  )
}

compute_pooled_within_group_variance <- function(
    normalized_counts,
    group_labels) {

  groups <- levels(factor(group_labels))

  sse <- rep(0, nrow(normalized_counts))
  residual_df <- 0L

  for (g in groups) {
    idx <- which(group_labels == g)

    if (length(idx) < 2L) next

    xg <- normalized_counts[, idx, drop = FALSE]
    mu_g <- rowMeans(xg)

    resid_g <- sweep(
      xg,
      1L,
      mu_g,
      "-"
    )

    sse <- sse + rowSums(resid_g^2)
    residual_df <- residual_df + length(idx) - 1L
  }

  if (residual_df < 2L) {
    stop("Pooled residual degrees of freedom < 2.")
  }

  V <- sse / residual_df
  V[!is.finite(V)] <- 0
  V <- pmax(V, 0)

  names(V) <- rownames(normalized_counts)

  list(
    variance = V,
    residual_df = residual_df
  )
}


# =============================================================================
# PC1 RANKING
# =============================================================================

compute_pc1_rank <- function(rank_matrix_arm) {
  pca <- stats::prcomp(
    t(rank_matrix_arm),
    center = TRUE,
    scale. = FALSE,
    rank. = 1
  )

  loading <- pca$rotation[, 1L]
  loading[!is.finite(loading)] <- 0

  abs_loading <- abs(loading)

  rank_order <- order(
    abs_loading,
    decreasing = FALSE
  )

  lambda1 <- pca$sdev[1L]^2

  P <- lambda1 * loading^2

  list(
    loading = loading,
    abs_loading = abs_loading,
    rank_order = rank_order,
    lambda1 = lambda1,
    pc1_variance_contribution = P
  )
}


# =============================================================================
# RAW-COUNT VARIANCE GEOMETRY
# =============================================================================

find_d2_zero_crossings <- function(x, d2) {
  ok <- is.finite(x) & is.finite(d2)

  x <- x[ok]
  d2 <- d2[ok]

  if (length(x) < 2L) {
    return(numeric(0))
  }

  out <- numeric(0)

  for (i in seq_len(length(x) - 1L)) {
    a <- d2[i]
    b <- d2[i + 1L]

    if (a == 0 || b == 0) next

    if (
      (a < 0 && b > 0) ||
      (a > 0 && b < 0)
    ) {
      frac <- abs(a) / (abs(a) + abs(b))

      out <- c(
        out,
        x[i] + frac * (x[i + 1L] - x[i])
      )
    }
  }

  sort(unique(out))
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

  raw_var[!is.finite(raw_var)] <- 0
  raw_var <- pmax(raw_var, 0)

  raw_var_ranked <- raw_var[rank_order]

  rank <- seq_along(rank_order)
  log_var <- log1p(raw_var_ranked)

  analysis_spline <- stats::smooth.spline(
    x = rank,
    y = log_var,
    spar = ANALYSIS_SPLINE_SPAR
  )

  dense_x <- seq(
    min(rank),
    max(rank),
    length.out = max(
      5000L,
      length(rank) * 4L
    )
  )

  dense_d2 <- as.numeric(
    stats::predict(
      analysis_spline,
      x = dense_x,
      deriv = 2
    )$y
  )

  crossings <- find_d2_zero_crossings(
    dense_x,
    dense_d2
  )

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

  list(
    curve = data.frame(
      rank = rank,
      raw_empirical_variance = raw_var_ranked,
      log1p_raw_empirical_variance = log_var,
      display_log1p_raw_empirical_variance = display_y,
      stringsAsFactors = FALSE
    ),
    crossings = crossings
  )
}

select_anchor_terminal_after_c2 <- function(
    crossings,
    c2,
    total_n) {

  z <- sort(
    unique(
      crossings[
        is.finite(crossings) &
        crossings > c2 &
        crossings < total_n
      ]
    )
  )

  if (length(z) < 2L) {
    stop(
      "Fewer than two raw-variance second-derivative zero crossings ",
      "were found after c2 = ",
      c2,
      "."
    )
  }

  anchor <- as.integer(round(z[1L]))
  terminal <- as.integer(round(z[2L]))

  anchor <- max(c2 + 1L, min(total_n - 1L, anchor))
  terminal <- max(anchor + 1L, min(total_n, terminal))

  list(
    anchor = anchor,
    terminal = terminal
  )
}


# =============================================================================
# GROUP-SPECIFIC PC1-NB DIVERGENCE
# =============================================================================

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

  endpoint_line <- seq(
    y[1L],
    y[length(y)],
    length.out = length(y)
  )

  y - endpoint_line
}

compute_group_analysis <- function(
    group_name,
    raw_counts_arm,
    normalized_counts_arm,
    pooled_variance) {

  rank_matrix <- normalize_cpm_log1p(raw_counts_arm)

  pc1 <- compute_pc1_rank(rank_matrix)
  rank_order <- pc1$rank_order

  geometry <- compute_raw_variance_geometry(
    raw_counts_arm = raw_counts_arm,
    rank_order = rank_order
  )

  mu_norm <- rowMeans(
    normalized_counts_arm,
    na.rm = TRUE
  )

  mu_norm[!is.finite(mu_norm)] <- 0
  mu_norm <- pmax(mu_norm, 0)

  P_ranked <- pc1$pc1_variance_contribution[rank_order]
  mu_ranked <- mu_norm[rank_order]
  V_pool_ranked <- pooled_variance[rank_order]

  E_ranked <- pmax(
    V_pool_ranked - mu_ranked,
    0
  )

  P_total <- sum(P_ranked)
  E_total <- sum(E_ranked)

  if (!is.finite(P_total) || P_total <= 0) {
    stop("PC1 variance mass undefined for group ", group_name)
  }

  if (!is.finite(E_total) || E_total <= 0) {
    stop("NB excess-variance mass undefined for group ", group_name)
  }

  p_mass <- P_ranked / P_total
  q_mass <- E_ranked / E_total

  F_P <- cumsum(p_mass)
  F_E <- cumsum(q_mass)

  D <- F_E - F_P

  rank <- seq_along(rank_order)

  display_D <- smooth_divergence_for_display(
    rank = rank,
    D = D
  )

  df <- geometry$curve %>%
    mutate(
      group = group_name,
      feature_id = rownames(raw_counts_arm)[rank_order],
      pc1_loading = pc1$loading[rank_order],
      abs_pc1_loading = pc1$abs_loading[rank_order],
      pc1_eigenvalue = pc1$lambda1,
      pc1_variance_contribution = P_ranked,
      pc1_variance_mass = p_mass,

      normalized_group_mean = mu_ranked,
      pooled_normalized_variance = V_pool_ranked,
      nb_excess_variance = E_ranked,
      nb_excess_variance_mass = q_mass,

      cumulative_pc1_mass = F_P,
      cumulative_nb_mass = F_E,
      cumulative_divergence = D,
      display_D = display_D
    )

  list(
    data = df,
    raw_variance_crossings = geometry$crossings,
    anchor = NA_integer_,
    terminal = NA_integer_
  )
}


# =============================================================================
# SHARED TWO-KNOT DIVERGENCE MODEL
# =============================================================================

piecewise_basis <- function(x, c1, c2) {
  cbind(
    intercept = 1,
    x = x,
    hinge1 = pmax(x - c1, 0),
    hinge2 = pmax(x - c2, 0)
  )
}

piecewise_sse <- function(
    par,
    x,
    D_mat,
    min_gap) {

  c1 <- par[1L]
  c2 <- par[2L]

  if (
    !is.finite(c1) ||
    !is.finite(c2) ||
    c1 <= 0 ||
    c2 >= 1 ||
    c2 - c1 <= min_gap
  ) {
    return(1e100)
  }

  X <- piecewise_basis(
    x = x,
    c1 = c1,
    c2 = c2
  )

  coef <- tryCatch(
    qr.coef(
      qr(X),
      D_mat
    ),
    error = function(e) NULL
  )

  if (is.null(coef) || any(!is.finite(coef))) {
    return(1e100)
  }

  resid <- D_mat - X %*% coef

  sum(resid^2)
}

fit_shared_knots <- function(group_results) {
  groups <- names(group_results)

  N_values <- vapply(
    group_results,
    function(z) nrow(z$data),
    integer(1)
  )

  if (length(unique(N_values)) != 1L) {
    stop("All groups must contain the same number of ranked features.")
  }

  N <- N_values[1L]

  x_full <- (seq_len(N) - 1) / (N - 1)

  D_full <- do.call(
    cbind,
    lapply(
      group_results,
      function(z) z$data$cumulative_divergence
    )
  )

  colnames(D_full) <- groups

  opt_n <- min(5000L, N)

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

  x_opt <- x_full[opt_idx]
  D_opt <- D_full[opt_idx, , drop = FALSE]

  min_gap <- max(
    4 / (N - 1),
    .Machine$double.eps^0.25
  )

  starts <- list(
    c(0.03, 0.97),
    c(0.08, 0.92),
    c(0.15, 0.85),
    c(0.25, 0.75),
    c(0.35, 0.65)
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
    function(z) z$value,
    numeric(1)
  )

  best <- coarse[[which.min(values)]]

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

  if (!is.finite(refined$value)) {
    stop("Shared-knot optimization failed.")
  }

  c1_rank <- as.integer(
    round(
      1 + refined$par[1L] * (N - 1)
    )
  )

  c2_rank <- as.integer(
    round(
      1 + refined$par[2L] * (N - 1)
    )
  )

  c1_rank <- max(
    2L,
    min(N - 2L, c1_rank)
  )

  c2_rank <- max(
    c1_rank + 1L,
    min(N - 1L, c2_rank)
  )

  c1_x <- (c1_rank - 1) / (N - 1)
  c2_x <- (c2_rank - 1) / (N - 1)

  X <- piecewise_basis(
    x_full,
    c1_x,
    c2_x
  )

  coef <- qr.coef(
    qr(X),
    D_full
  )

  fitted <- X %*% coef

  list(
    c1 = c1_rank,
    c2 = c2_rank,
    x = x_full,
    fitted = fitted,
    SSE = sum((D_full - fitted)^2),
    groups = groups
  )
}


# =============================================================================
# NB MEAN-VARIANCE SCALING
# =============================================================================

estimate_nb_exponent <- function(mu, E) {
  keep <- (
    is.finite(mu) &
    is.finite(E) &
    mu > 0 &
    E > 0
  )

  if (sum(keep) < 10L) {
    return(NA_real_)
  }

  fit <- stats::lm(
    log(E[keep]) ~ log(mu[keep])
  )

  unname(
    stats::coef(fit)[2L]
  )
}

get_region_p <- function(df, keep) {
  estimate_nb_exponent(
    mu = df$normalized_group_mean[keep],
    E = df$nb_excess_variance[keep]
  )
}


# =============================================================================
# TOP-k OPTIMIZATION AND JOINT / DISJOINT CLASSIFICATION
# =============================================================================

make_rank_map <- function(df) {
  stats::setNames(
    df$rank,
    df$feature_id
  )
}

rank_to_region <- function(rank, c1, c2) {
  ifelse(
    rank < c1,
    "Remainder",
    ifelse(
      rank <= c2,
      "Divergence",
      "LeadingEdge"
    )
  )
}

cumulative_activation <- function(depth, K) {
  depth <- as.integer(depth)

  keep <- (
    is.finite(depth) &
    depth >= 1L &
    depth <= K
  )

  if (!any(keep)) {
    return(rep(0L, K))
  }

  cumsum(
    tabulate(
      depth[keep],
      nbins = K
    )
  )
}

active_interval_count <- function(starts, ends, K) {
  # Count intervals active for start <= k < end.
  # end may equal K+1.
  if (length(starts) == 0L) {
    return(rep(0L, K))
  }

  starts <- as.integer(starts)
  ends <- as.integer(ends)

  valid <- (
    is.finite(starts) &
    is.finite(ends) &
    starts >= 1L &
    starts <= K &
    ends > starts
  )

  starts <- starts[valid]
  ends <- pmin(
    ends[valid],
    K + 1L
  )

  if (length(starts) == 0L) {
    return(rep(0L, K))
  }

  diff_vec <- integer(K + 1L)

  start_tab <- tabulate(
    starts,
    nbins = K + 1L
  )

  end_tab <- tabulate(
    ends,
    nbins = K + 1L
  )

  diff_vec <- (
    diff_vec +
    start_tab -
    end_tab
  )

  cumsum(diff_vec)[seq_len(K)]
}


scan_pair_cutoffs <- function(
    control_df,
    treatment_df,
    c1,
    c2,
    comparison_name,
    control_group,
    treatment_group) {

  if (nrow(control_df) != nrow(treatment_df)) {
    stop(
      "Control and treatment rankings have different feature counts."
    )
  }

  if (!setequal(
    control_df$feature_id,
    treatment_df$feature_id
  )) {
    stop(
      "Control and treatment rankings do not contain the same feature IDs."
    )
  }

  N <- nrow(control_df)

  # Own-arm selection must originate strictly beyond c2.
  K <- as.integer(N - c2)

  if (K < 1L) {
    stop(
      "No candidate top-k depth exists beyond c2."
    )
  }

  rank_control <- make_rank_map(
    control_df
  )

  rank_treatment <- make_rank_map(
    treatment_df
  )

  ids <- control_df$feature_id

  rC <- as.integer(
    unname(
      rank_control[ids]
    )
  )

  rT <- as.integer(
    unname(
      rank_treatment[ids]
    )
  )

  if (
    any(!is.finite(rC)) ||
    any(!is.finite(rT))
  ) {
    stop(
      "Non-finite rank encountered in ",
      comparison_name
    )
  }

  # Entry depth:
  # rank N enters at k=1
  # rank 1 enters at k=N
  dC <- N - rC + 1L
  dT <- N - rT + 1L

  k <- seq_len(K)

  # -----------------------------------------------------------------------
  # JOINT
  # -----------------------------------------------------------------------
  # Joint membership activates once BOTH top-k selections contain the site.
  joint_depth <- pmax(
    dC,
    dT
  )

  joint_n <- cumulative_activation(
    joint_depth,
    K
  )

  # -----------------------------------------------------------------------
  # DISJOINT, OPPOSITE ARM ALSO IN LEADING EDGE
  # -----------------------------------------------------------------------
  # Control-only while treatment has not yet admitted the site.
  idx_le_C <- which(
    dC < dT &
    dC <= K &
    dT <= K
  )

  disjoint_control_opposite_le_n <- active_interval_count(
    starts = dC[idx_le_C],
    ends = dT[idx_le_C],
    K = K
  )

  # Treatment-only while control has not yet admitted the site.
  idx_le_T <- which(
    dT < dC &
    dT <= K &
    dC <= K
  )

  disjoint_treatment_opposite_le_n <- active_interval_count(
    starts = dT[idx_le_T],
    ends = dC[idx_le_T],
    K = K
  )

  # -----------------------------------------------------------------------
  # DISJOINT, OPPOSITE ARM IN DIVERGENCE INTERVAL
  # -----------------------------------------------------------------------
  # These sites are intentionally PERMISSIBLE.
  #
  # The selecting arm is in its own top-k subset and therefore >c2.
  # The opposite arm lies between c1 and c2 and never enters a top-k list
  # because candidate k is capped at N-c2.
  idx_div_C <- which(
    dC <= K &
    rT >= c1 &
    rT <= c2
  )

  disjoint_control_opposite_divergence_n <-
    cumulative_activation(
      dC[idx_div_C],
      K
    )

  idx_div_T <- which(
    dT <= K &
    rC >= c1 &
    rC <= c2
  )

  disjoint_treatment_opposite_divergence_n <-
    cumulative_activation(
      dT[idx_div_T],
      K
    )

  # -----------------------------------------------------------------------
  # HARD CROSS-REGIME CONTAMINATION: OPPOSITE ARM IN REMAINDER
  # -----------------------------------------------------------------------
  # Only r_opposite < c1 is penalized.
  idx_rem_C <- which(
    dC <= K &
    rT < c1
  )

  remainder_cross_control_n <- cumulative_activation(
    dC[idx_rem_C],
    K
  )

  idx_rem_T <- which(
    dT <= K &
    rC < c1
  )

  remainder_cross_treatment_n <- cumulative_activation(
    dT[idx_rem_T],
    K
  )

  # -----------------------------------------------------------------------
  # AGGREGATES
  # -----------------------------------------------------------------------

  disjoint_opposite_le_n <- (
    disjoint_control_opposite_le_n +
    disjoint_treatment_opposite_le_n
  )

  disjoint_opposite_divergence_n <- (
    disjoint_control_opposite_divergence_n +
    disjoint_treatment_opposite_divergence_n
  )

  permissible_disjoint_control_n <- (
    disjoint_control_opposite_le_n +
    disjoint_control_opposite_divergence_n
  )

  permissible_disjoint_treatment_n <- (
    disjoint_treatment_opposite_le_n +
    disjoint_treatment_opposite_divergence_n
  )

  permissible_disjoint_n <- (
    permissible_disjoint_control_n +
    permissible_disjoint_treatment_n
  )

  # Benefit:
  # Joint + all permissible disjoint sites.
  good_n <- (
    joint_n +
    permissible_disjoint_n
  )

  # Cost:
  # only true crossing into opposite-arm remainder.
  remainder_cross_n <- (
    remainder_cross_control_n +
    remainder_cross_treatment_n
  )

  union_n <- (
    good_n +
    remainder_cross_n
  )

  # Independent union-size check.
  union_depth <- pmin(
    dC,
    dT
  )

  union_check <- cumulative_activation(
    union_depth,
    K
  )

  if (!all(
    union_n == union_check
  )) {
    stop(
      "Internal union-count mismatch in ",
      comparison_name
    )
  }

  data.frame(
    comparison = comparison_name,
    control_group = control_group,
    treatment_group = treatment_group,

    k = k,
    cutoff_rank = N - k + 1L,

    joint_n = joint_n,

    disjoint_control_opposite_le_n =
      disjoint_control_opposite_le_n,

    disjoint_treatment_opposite_le_n =
      disjoint_treatment_opposite_le_n,

    disjoint_opposite_le_n =
      disjoint_opposite_le_n,

    disjoint_control_opposite_divergence_n =
      disjoint_control_opposite_divergence_n,

    disjoint_treatment_opposite_divergence_n =
      disjoint_treatment_opposite_divergence_n,

    disjoint_opposite_divergence_n =
      disjoint_opposite_divergence_n,

    permissible_disjoint_control_n =
      permissible_disjoint_control_n,

    permissible_disjoint_treatment_n =
      permissible_disjoint_treatment_n,

    permissible_disjoint_n =
      permissible_disjoint_n,

    good_n = good_n,

    remainder_cross_control_n =
      remainder_cross_control_n,

    remainder_cross_treatment_n =
      remainder_cross_treatment_n,

    remainder_cross_n =
      remainder_cross_n,

    union_n = union_n,

    retained_fraction = ifelse(
      union_n > 0,
      good_n / union_n,
      NA_real_
    ),

    remainder_cross_fraction = ifelse(
      union_n > 0,
      remainder_cross_n / union_n,
      NA_real_
    ),

    divergence_disjoint_fraction = ifelse(
      union_n > 0,
      disjoint_opposite_divergence_n / union_n,
      NA_real_
    ),

    jaccard_top_k = ifelse(
      union_n > 0,
      joint_n / union_n,
      NA_real_
    ),

    stringsAsFactors = FALSE
  )
}


mark_pareto_frontier <- function(
    scan_df,
    good_col = "good_n",
    cost_col = "remainder_cross_n") {

  if (nrow(scan_df) < 1L) {
    stop(
      "Empty cutoff scan."
    )
  }

  tmp <- scan_df %>%
    transmute(
      row_id = row_number(),
      k = k,
      good = .data[[good_col]],
      cost = .data[[cost_col]]
    ) %>%
    arrange(
      cost,
      desc(good),
      desc(k)
    ) %>%
    group_by(cost) %>%
    slice(1L) %>%
    ungroup() %>%
    arrange(
      cost,
      desc(good)
    )

  running_best_before <- c(
    -Inf,
    head(
      cummax(tmp$good),
      -1L
    )
  )

  tmp$is_frontier_coord <- (
    tmp$good >
    running_best_before
  )

  frontier <- tmp %>%
    filter(
      is_frontier_coord
    ) %>%
    arrange(
      cost,
      good,
      k
    )

  key_all <- paste(
    scan_df[[cost_col]],
    scan_df[[good_col]],
    sep = "::"
  )

  key_frontier <- paste(
    frontier$cost,
    frontier$good,
    sep = "::"
  )

  out <- scan_df
  out$is_pareto <- (
    key_all %in%
    key_frontier
  )

  list(
    scan = out,
    frontier = frontier
  )
}


select_pareto_changepoint <- function(
    scan_df,
    good_col = "good_n",
    cost_col = "remainder_cross_n",
    min_segment_points = 3L) {

  marked <- mark_pareto_frontier(
    scan_df,
    good_col = good_col,
    cost_col = cost_col
  )

  frontier <- marked$frontier %>%
    arrange(
      cost,
      good,
      k
    )

  out <- marked$scan

  zero_idx <- which(
    out[[cost_col]] == 0
  )

  zero_max_k <- if (
    length(zero_idx) > 0L
  ) {
    max(
      out$k[zero_idx]
    )
  } else {
    NA_integer_
  }

  n <- nrow(frontier)

  candidates <- data.frame()
  valid <- data.frame()

  # Defaults used by any fallback.
  selected_k <- NA_integer_
  selected_good <- NA_real_
  selected_cost <- NA_real_
  breakpoint_cost <- NA_real_
  slope_before <- NA_real_
  slope_after <- NA_real_
  slope_ratio <- NA_real_
  bic_linear <- NA_real_
  bic_segmented <- NA_real_
  delta_bic <- NA_real_
  shape_supported <- FALSE
  bic_supports_segmented <- FALSE
  selection_method <- NA_character_

  frontier$segmented_fit_good <- NA_real_
  frontier$marginal_good_per_cost <- NA_real_
  frontier$is_changepoint <- FALSE

  # -----------------------------------------------------------------------
  # Need enough distinct Pareto coordinates to estimate a two-segment model.
  # -----------------------------------------------------------------------
  enough_frontier <- (
    n >=
    (2L * min_segment_points + 1L)
  )

  if (enough_frontier) {
    x <- as.numeric(
      frontier$cost
    )

    y <- as.numeric(
      frontier$good
    )

    if (
      any(!is.finite(x)) ||
      any(!is.finite(y))
    ) {
      stop(
        "Non-finite Pareto-frontier values encountered."
      )
    }

    if (
      any(diff(x) <= 0)
    ) {
      stop(
        "Pareto cost must be strictly increasing after duplicate-cost ",
        "reduction."
      )
    }

    # Single-line reference model.
    X_linear <- cbind(
      intercept = 1,
      x = x
    )

    linear_fit <- stats::lm.fit(
      x = X_linear,
      y = y
    )

    rss_linear <- sum(
      linear_fit$residuals^2
    )

    rss_linear <- max(
      rss_linear,
      .Machine$double.eps
    )

    # Two fitted regression coefficients.
    bic_linear <- (
      n * log(
        rss_linear / n
      ) +
      2 * log(n)
    )

    # Efficient exhaustive broken-stick scan.
    total_x <- sum(x)
    total_x2 <- sum(x^2)
    total_y <- sum(y)
    total_xy <- sum(x * y)
    total_y2 <- sum(y^2)

    tail_x <- rev(
      cumsum(
        rev(x)
      )
    )

    tail_x2 <- rev(
      cumsum(
        rev(x^2)
      )
    )

    tail_y <- rev(
      cumsum(
        rev(y)
      )
    )

    tail_xy <- rev(
      cumsum(
        rev(x * y)
      )
    )

    candidate_idx <- seq.int(
      min_segment_points,
      n - min_segment_points
    )

    candidate_rows <- vector(
      "list",
      length(candidate_idx)
    )

    for (jj in seq_along(
      candidate_idx
    )) {
      i <- candidate_idx[jj]
      tau <- x[i]

      n_right <- n - i

      sx_right <- tail_x[i + 1L]
      sx2_right <- tail_x2[i + 1L]
      sy_right <- tail_y[i + 1L]
      sxy_right <- tail_xy[i + 1L]

      sum_h <- (
        sx_right -
        n_right * tau
      )

      sum_xh <- (
        sx2_right -
        tau * sx_right
      )

      sum_h2 <- (
        sx2_right -
        2 * tau * sx_right +
        n_right * tau^2
      )

      sum_hy <- (
        sxy_right -
        tau * sy_right
      )

      XtX <- matrix(
        c(
          n, total_x, sum_h,
          total_x, total_x2, sum_xh,
          sum_h, sum_xh, sum_h2
        ),
        nrow = 3L,
        byrow = TRUE
      )

      Xty <- c(
        total_y,
        total_xy,
        sum_hy
      )

      beta <- tryCatch(
        solve(
          XtX,
          Xty
        ),
        error = function(e) {
          rep(
            NA_real_,
            3L
          )
        }
      )

      if (
        any(!is.finite(beta))
      ) {
        candidate_rows[[jj]] <- data.frame(
          frontier_index = i,
          tau = tau,
          k = frontier$k[i],
          rss = Inf,
          bic = Inf,
          slope_before = NA_real_,
          slope_after = NA_real_,
          slope_change = NA_real_,
          shape_valid = FALSE,
          beta0 = NA_real_,
          beta1 = NA_real_,
          gamma = NA_real_,
          stringsAsFactors = FALSE
        )

        next
      }

      rss <- (
        total_y2 -
        sum(
          beta * Xty
        )
      )

      rss <- max(
        rss,
        .Machine$double.eps
      )

      # Four effective parameters:
      # intercept, pre-slope, slope change, breakpoint.
      bic <- (
        n * log(
          rss / n
        ) +
        4 * log(n)
      )

      slope_pre <- beta[2L]
      slope_post <- (
        beta[2L] +
        beta[3L]
      )

      shape_valid <- (
        is.finite(slope_pre) &&
        is.finite(slope_post) &&
        slope_pre >= 0 &&
        slope_post >= 0 &&
        slope_post <
        slope_pre
      )

      candidate_rows[[jj]] <- data.frame(
        frontier_index = i,
        tau = tau,
        k = frontier$k[i],
        rss = rss,
        bic = bic,
        slope_before = slope_pre,
        slope_after = slope_post,
        slope_change = (
          slope_post -
          slope_pre
        ),
        shape_valid = shape_valid,
        beta0 = beta[1L],
        beta1 = beta[2L],
        gamma = beta[3L],
        stringsAsFactors = FALSE
      )
    }

    candidates <- bind_rows(
      candidate_rows
    )

    valid <- candidates %>%
      filter(
        shape_valid,
        is.finite(bic)
      )

    shape_supported <- (
      nrow(valid) > 0L
    )

    if (shape_supported) {
      chosen <- valid %>%
        arrange(
          bic,
          desc(k)
        ) %>%
        slice(1L)

      bic_segmented <- chosen$bic[1L]
      delta_bic <- (
        bic_linear -
        bic_segmented
      )

      bic_supports_segmented <- (
        is.finite(delta_bic) &&
        delta_bic > 0
      )

      if (bic_supports_segmented) {
        i_star <- as.integer(
          chosen$frontier_index[1L]
        )

        selected_k <- as.integer(
          chosen$k[1L]
        )

        selected_good <- frontier$good[i_star]
        selected_cost <- frontier$cost[i_star]
        breakpoint_cost <- chosen$tau[1L]

        slope_before <- chosen$slope_before[1L]
        slope_after <- chosen$slope_after[1L]

        slope_ratio <- ifelse(
          slope_before > 0,
          slope_after / slope_before,
          NA_real_
        )

        beta0 <- chosen$beta0[1L]
        beta1 <- chosen$beta1[1L]
        gamma <- chosen$gamma[1L]

        frontier$segmented_fit_good <- (
          beta0 +
          beta1 * x +
          gamma *
          pmax(
            x - breakpoint_cost,
            0
          )
        )

        frontier$marginal_good_per_cost <- c(
          NA_real_,
          diff(y) / diff(x)
        )

        frontier$is_changepoint[
          i_star
        ] <- TRUE

        selection_method <-
          "segmented_pareto_changepoint"
      }
    }
  }

  # -----------------------------------------------------------------------
  # Conservative fallback if there is no supported segmented change-point.
  # -----------------------------------------------------------------------
  if (!is.finite(selected_k)) {
    if (is.finite(zero_max_k)) {
      selected_k <- as.integer(
        zero_max_k
      )

      selected_row <- out %>%
        filter(
          k == selected_k
        ) %>%
        slice(1L)

      selected_good <- selected_row[[good_col]][1L]
      selected_cost <- selected_row[[cost_col]][1L]

      selection_method <-
        "largest_zero_remainder_crossing_k"
    } else {
      # No zero-cost solution.  Use the best shape-valid segmented
      # approximation if available, but explicitly flag it as descriptive.
      if (
        nrow(valid) > 0L
      ) {
        chosen <- valid %>%
          arrange(
            bic,
            desc(k)
          ) %>%
          slice(1L)

        i_star <- as.integer(
          chosen$frontier_index[1L]
        )

        selected_k <- as.integer(
          chosen$k[1L]
        )

        selected_good <- frontier$good[i_star]
        selected_cost <- frontier$cost[i_star]
        breakpoint_cost <- chosen$tau[1L]

        slope_before <- chosen$slope_before[1L]
        slope_after <- chosen$slope_after[1L]

        slope_ratio <- ifelse(
          slope_before > 0,
          slope_after / slope_before,
          NA_real_
        )

        bic_segmented <- chosen$bic[1L]
        delta_bic <- (
          bic_linear -
          bic_segmented
        )

        bic_supports_segmented <- (
          is.finite(delta_bic) &&
          delta_bic > 0
        )

        beta0 <- chosen$beta0[1L]
        beta1 <- chosen$beta1[1L]
        gamma <- chosen$gamma[1L]

        x <- as.numeric(
          frontier$cost
        )

        y <- as.numeric(
          frontier$good
        )

        frontier$segmented_fit_good <- (
          beta0 +
          beta1 * x +
          gamma *
          pmax(
            x - breakpoint_cost,
            0
          )
        )

        frontier$marginal_good_per_cost <- c(
          NA_real_,
          diff(y) / diff(x)
        )

        frontier$is_changepoint[
          i_star
        ] <- TRUE

        selection_method <-
          "segmented_descriptive_no_BIC_support"
      } else {
        # Last-resort deterministic Pareto point:
        # smallest remainder crossing count; among ties greatest benefit,
        # among ties largest k.
        fallback <- out %>%
          arrange(
            .data[[cost_col]],
            desc(.data[[good_col]]),
            desc(k)
          ) %>%
          slice(1L)

        selected_k <- as.integer(
          fallback$k[1L]
        )

        selected_good <- fallback[[good_col]][1L]
        selected_cost <- fallback[[cost_col]][1L]

        selection_method <-
          "minimum_remainder_crossing_fallback"
      }
    }
  }

  # Map frontier quantities back onto the full scan.
  frontier_key <- paste(
    frontier$cost,
    frontier$good,
    sep = "::"
  )

  out_key <- paste(
    out[[cost_col]],
    out[[good_col]],
    sep = "::"
  )

  m <- match(
    out_key,
    frontier_key
  )

  out$segmented_fit_good <-
    frontier$segmented_fit_good[m]

  out$marginal_good_per_cost <-
    frontier$marginal_good_per_cost[m]

  out$is_selected_cutoff <- (
    out$k ==
    selected_k
  )

  out$selected_method <-
    selection_method

  out$breakpoint_remainder_cross_n <-
    breakpoint_cost

  out$slope_before_breakpoint <-
    slope_before

  out$slope_after_breakpoint <-
    slope_after

  out$slope_ratio_after_before <-
    slope_ratio

  out$bic_linear <- bic_linear
  out$bic_segmented <- bic_segmented

  out$delta_BIC_linear_minus_segmented <-
    delta_bic

  out$shape_supported <-
    shape_supported

  out$bic_supports_segmented <-
    bic_supports_segmented

  list(
    scan = out,
    frontier = frontier,
    candidates = candidates,

    selected_k = as.integer(
      selected_k
    ),

    selected_good = selected_good,
    selected_cost = selected_cost,

    selection_method = selection_method,

    breakpoint_remainder_cross_n =
      breakpoint_cost,

    slope_before = slope_before,
    slope_after = slope_after,

    slope_ratio_after_before =
      slope_ratio,

    bic_linear = bic_linear,
    bic_segmented = bic_segmented,
    delta_BIC = delta_bic,

    shape_supported = shape_supported,

    bic_supports_segmented =
      bic_supports_segmented,

    max_zero_remainder_crossing_k =
      zero_max_k
  )
}


aggregate_global_cutoff_scan <- function(
    pair_scans) {

  all_pair_scan <- bind_rows(
    pair_scans
  )

  global <- all_pair_scan %>%
    group_by(k) %>%
    summarise(
      cutoff_rank = first(
        cutoff_rank
      ),

      joint_n = sum(
        joint_n
      ),

      disjoint_control_opposite_le_n = sum(
        disjoint_control_opposite_le_n
      ),

      disjoint_treatment_opposite_le_n = sum(
        disjoint_treatment_opposite_le_n
      ),

      disjoint_opposite_le_n = sum(
        disjoint_opposite_le_n
      ),

      disjoint_control_opposite_divergence_n = sum(
        disjoint_control_opposite_divergence_n
      ),

      disjoint_treatment_opposite_divergence_n = sum(
        disjoint_treatment_opposite_divergence_n
      ),

      disjoint_opposite_divergence_n = sum(
        disjoint_opposite_divergence_n
      ),

      permissible_disjoint_control_n = sum(
        permissible_disjoint_control_n
      ),

      permissible_disjoint_treatment_n = sum(
        permissible_disjoint_treatment_n
      ),

      permissible_disjoint_n = sum(
        permissible_disjoint_n
      ),

      good_n = sum(
        good_n
      ),

      remainder_cross_n = sum(
        remainder_cross_n
      ),

      union_n = sum(
        union_n
      ),

      .groups = "drop"
    ) %>%
    mutate(
      comparison = "GLOBAL",
      control_group = NA_character_,
      treatment_group = NA_character_,

      retained_fraction = ifelse(
        union_n > 0,
        good_n / union_n,
        NA_real_
      ),

      remainder_cross_fraction = ifelse(
        union_n > 0,
        remainder_cross_n / union_n,
        NA_real_
      ),

      divergence_disjoint_fraction = ifelse(
        union_n > 0,
        disjoint_opposite_divergence_n / union_n,
        NA_real_
      ),

      jaccard_top_k = ifelse(
        union_n > 0,
        joint_n / union_n,
        NA_real_
      )
    )

  # Pair-direction fields do not have a unique global direction.
  global$remainder_cross_control_n <-
    NA_integer_

  global$remainder_cross_treatment_n <-
    NA_integer_

  global
}


classify_pair_at_k <- function(
    control_df,
    treatment_df,
    k,
    c1,
    c2,
    comparison_name,
    control_group,
    treatment_group) {

  if (
    nrow(control_df) !=
    nrow(treatment_df)
  ) {
    stop(
      "Control and treatment rankings have different feature counts."
    )
  }

  if (!setequal(
    control_df$feature_id,
    treatment_df$feature_id
  )) {
    stop(
      "Control and treatment rankings do not contain the same feature IDs."
    )
  }

  N <- nrow(control_df)
  Kmax <- N - c2

  if (
    k < 1L ||
    k > Kmax
  ) {
    stop(
      "k must satisfy 1 <= k <= N-c2. Received k=",
      k,
      "; N-c2=",
      Kmax,
      "."
    )
  }

  control_top <- tail(
    control_df$feature_id,
    k
  )

  treatment_top <- tail(
    treatment_df$feature_id,
    k
  )

  union_ids <- union(
    control_top,
    treatment_top
  )

  rank_control <- make_rank_map(
    control_df
  )

  rank_treatment <- make_rank_map(
    treatment_df
  )

  rC <- as.integer(
    unname(
      rank_control[union_ids]
    )
  )

  rT <- as.integer(
    unname(
      rank_treatment[union_ids]
    )
  )

  in_control <- (
    union_ids %in%
    control_top
  )

  in_treatment <- (
    union_ids %in%
    treatment_top
  )

  joint <- (
    in_control &
    in_treatment
  )

  control_only <- (
    in_control &
    !in_treatment
  )

  treatment_only <- (
    in_treatment &
    !in_control
  )

  control_region <- rank_to_region(
    rC,
    c1,
    c2
  )

  treatment_region <- rank_to_region(
    rT,
    c1,
    c2
  )

  base_class <- ifelse(
    joint,
    "Joint",
    ifelse(
      control_only,
      paste0(
        "Disjoint_",
        control_group
      ),
      paste0(
        "Disjoint_",
        treatment_group
      )
    )
  )

  opposite_region <- rep(
    NA_character_,
    length(union_ids)
  )

  opposite_region[
    control_only
  ] <- treatment_region[
    control_only
  ]

  opposite_region[
    treatment_only
  ] <- control_region[
    treatment_only
  ]

  disjoint_opposite_leading_edge <- (
    (control_only & rT > c2) |
    (treatment_only & rC > c2)
  )

  disjoint_opposite_divergence <- (
    (
      control_only &
      rT >= c1 &
      rT <= c2
    ) |
    (
      treatment_only &
      rC >= c1 &
      rC <= c2
    )
  )

  cross_into_remainder <- (
    (
      control_only &
      rT < c1
    ) |
    (
      treatment_only &
      rC < c1
    )
  )

  retained_for_analysis <- (
    joint |
    disjoint_opposite_leading_edge |
    disjoint_opposite_divergence
  )

  analysis_class <- rep(
    NA_character_,
    length(union_ids)
  )

  analysis_class[
    joint
  ] <- "Joint"

  analysis_class[
    control_only &
    rT > c2
  ] <- paste0(
    "Disjoint_",
    control_group,
    "_OppositeLeadingEdge"
  )

  analysis_class[
    treatment_only &
    rC > c2
  ] <- paste0(
    "Disjoint_",
    treatment_group,
    "_OppositeLeadingEdge"
  )

  analysis_class[
    control_only &
    rT >= c1 &
    rT <= c2
  ] <- paste0(
    "Disjoint_",
    control_group,
    "_OppositeDivergence"
  )

  analysis_class[
    treatment_only &
    rC >= c1 &
    rC <= c2
  ] <- paste0(
    "Disjoint_",
    treatment_group,
    "_OppositeDivergence"
  )

  analysis_class[
    control_only &
    rT < c1
  ] <- paste0(
    "Excluded_",
    control_group,
    "_OppositeRemainder"
  )

  analysis_class[
    treatment_only &
    rC < c1
  ] <- paste0(
    "Excluded_",
    treatment_group,
    "_OppositeRemainder"
  )

  out <- data.frame(
    comparison = comparison_name,
    feature_id = union_ids,

    selected_k = as.integer(k),

    cutoff_rank = rank_cutoff_from_k(
      N,
      k
    ),

    control_group = control_group,
    treatment_group = treatment_group,

    base_class = base_class,
    analysis_class = analysis_class,

    control_rank = rC,
    treatment_rank = rT,

    control_region = control_region,
    treatment_region = treatment_region,

    opposite_region = opposite_region,

    control_top_k = in_control,
    treatment_top_k = in_treatment,

    disjoint_opposite_leading_edge =
      disjoint_opposite_leading_edge,

    disjoint_opposite_divergence =
      disjoint_opposite_divergence,

    cross_into_remainder =
      cross_into_remainder,

    retained_for_analysis =
      retained_for_analysis,

    stringsAsFactors = FALSE
  )

  if (
    any(
      is.na(
        out$analysis_class
      )
    )
  ) {
    stop(
      "Unclassified union site encountered in ",
      comparison_name
    )
  }

  if (
    any(
      out$retained_for_analysis &
      out$cross_into_remainder
    )
  ) {
    stop(
      "A retained site was also classified as a remainder crossing in ",
      comparison_name
    )
  }

  out
}


summarize_classification <- function(
    class_df) {

  if (
    nrow(class_df) < 1L
  ) {
    stop(
      "Cannot summarize empty classification."
    )
  }

  joint_n <- sum(
    class_df$base_class ==
    "Joint"
  )

  control_only <- (
    class_df$control_top_k &
    !class_df$treatment_top_k
  )

  treatment_only <- (
    class_df$treatment_top_k &
    !class_df$control_top_k
  )

  disjoint_control_opposite_le_n <- sum(
    control_only &
    class_df$disjoint_opposite_leading_edge
  )

  disjoint_treatment_opposite_le_n <- sum(
    treatment_only &
    class_df$disjoint_opposite_leading_edge
  )

  disjoint_control_opposite_divergence_n <- sum(
    control_only &
    class_df$disjoint_opposite_divergence
  )

  disjoint_treatment_opposite_divergence_n <- sum(
    treatment_only &
    class_df$disjoint_opposite_divergence
  )

  disjoint_opposite_le_n <- (
    disjoint_control_opposite_le_n +
    disjoint_treatment_opposite_le_n
  )

  disjoint_opposite_divergence_n <- (
    disjoint_control_opposite_divergence_n +
    disjoint_treatment_opposite_divergence_n
  )

  permissible_disjoint_control_n <- (
    disjoint_control_opposite_le_n +
    disjoint_control_opposite_divergence_n
  )

  permissible_disjoint_treatment_n <- (
    disjoint_treatment_opposite_le_n +
    disjoint_treatment_opposite_divergence_n
  )

  permissible_disjoint_n <- (
    permissible_disjoint_control_n +
    permissible_disjoint_treatment_n
  )

  good_n <- (
    joint_n +
    permissible_disjoint_n
  )

  remainder_cross_n <- sum(
    class_df$cross_into_remainder
  )

  union_n <- nrow(
    class_df
  )

  if (
    good_n +
    remainder_cross_n !=
    union_n
  ) {
    stop(
      "Classification counts do not sum to union size."
    )
  }

  data.frame(
    k = unique(
      class_df$selected_k
    )[1L],

    cutoff_rank = unique(
      class_df$cutoff_rank
    )[1L],

    joint_n = joint_n,

    disjoint_control_opposite_le_n =
      disjoint_control_opposite_le_n,

    disjoint_treatment_opposite_le_n =
      disjoint_treatment_opposite_le_n,

    disjoint_opposite_le_n =
      disjoint_opposite_le_n,

    disjoint_control_opposite_divergence_n =
      disjoint_control_opposite_divergence_n,

    disjoint_treatment_opposite_divergence_n =
      disjoint_treatment_opposite_divergence_n,

    disjoint_opposite_divergence_n =
      disjoint_opposite_divergence_n,

    permissible_disjoint_control_n =
      permissible_disjoint_control_n,

    permissible_disjoint_treatment_n =
      permissible_disjoint_treatment_n,

    permissible_disjoint_n =
      permissible_disjoint_n,

    good_n = good_n,

    remainder_cross_n =
      remainder_cross_n,

    union_n = union_n,

    retained_fraction = (
      good_n /
      union_n
    ),

    remainder_cross_fraction = (
      remainder_cross_n /
      union_n
    ),

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# RESULTS TABLES
# =============================================================================

build_arm_scaling_table <- function(
    group_results,
    comparisons,
    c1,
    c2,
    selected_k,
    paper_k) {

  rows <- list()

  for (comparison_name in names(
    comparisons
  )) {
    mapping <- comparisons[[
      comparison_name
    ]]

    for (arm in names(
      mapping
    )) {
      g <- unname(
        mapping[[arm]]
      )

      df <- group_results[[g]]$data
      N <- nrow(df)

      selected_rank <- rank_cutoff_from_k(
        N,
        selected_k
      )

      paper_rank <- rank_cutoff_from_k(
        N,
        paper_k
      )

      rows[[
        length(rows) + 1L
      ]] <- data.frame(
        comparison = comparison_name,
        arm = arm,
        group = g,

        p_remainder = get_region_p(
          df,
          df$rank < c1
        ),

        p_divergence = get_region_p(
          df,
          df$rank >= c1 &
          df$rank <= c2
        ),

        p_leading_edge = get_region_p(
          df,
          df$rank > c2
        ),

        p_selected_k = get_region_p(
          df,
          df$rank >=
          selected_rank
        ),

        p_paper_5000 = get_region_p(
          df,
          df$rank >=
          paper_rank
        ),

        stringsAsFactors = FALSE
      )
    }
  }

  bind_rows(
    rows
  )
}


build_timepoint_table <- function(
    group_results,
    comparisons,
    pair_optima,
    selected_summaries,
    paper_summaries,
    global_k,
    global_opt,
    c1,
    c2,
    paper_k,
    N) {

  rows <- list()

  for (comparison_name in names(
    comparisons
  )) {
    mapping <- comparisons[[
      comparison_name
    ]]

    control_group <- unname(
      mapping[["control"]]
    )

    treatment_group <- unname(
      mapping[["treatment"]]
    )

    pair_opt <- pair_optima[[
      comparison_name
    ]]

    selected_summary <- selected_summaries[[
      comparison_name
    ]]

    paper_summary <- paper_summaries[[
      comparison_name
    ]]

    rows[[
      length(rows) + 1L
    ]] <- data.frame(
      comparison = comparison_name,
      control_group = control_group,
      treatment_group = treatment_group,

      c1 = c1,
      c2 = c2,
      max_candidate_k = N - c2,

      control_anchor =
        group_results[[control_group]]$anchor,

      control_terminal =
        group_results[[control_group]]$terminal,

      treatment_anchor =
        group_results[[treatment_group]]$anchor,

      treatment_terminal =
        group_results[[treatment_group]]$terminal,

      pairwise_selected_k =
        pair_opt$selected_k,

      pairwise_selection_method =
        pair_opt$selection_method,

      pairwise_breakpoint_remainder_cross_n =
        pair_opt$breakpoint_remainder_cross_n,

      pairwise_slope_before =
        pair_opt$slope_before,

      pairwise_slope_after =
        pair_opt$slope_after,

      pairwise_slope_ratio_after_before =
        pair_opt$slope_ratio_after_before,

      pairwise_delta_BIC =
        pair_opt$delta_BIC,

      pairwise_shape_supported =
        pair_opt$shape_supported,

      pairwise_BIC_supports_segmented =
        pair_opt$bic_supports_segmented,

      pairwise_max_zero_remainder_crossing_k =
        pair_opt$max_zero_remainder_crossing_k,

      global_selected_k =
        global_k,

      global_selected_cutoff_rank =
        rank_cutoff_from_k(
          N,
          global_k
        ),

      global_selection_method =
        global_opt$selection_method,

      selected_joint_n =
        selected_summary$joint_n,

      selected_disjoint_opposite_le_n =
        selected_summary$disjoint_opposite_le_n,

      selected_disjoint_opposite_divergence_n =
        selected_summary$disjoint_opposite_divergence_n,

      selected_permissible_disjoint_n =
        selected_summary$permissible_disjoint_n,

      selected_good_n =
        selected_summary$good_n,

      selected_remainder_cross_n =
        selected_summary$remainder_cross_n,

      selected_union_n =
        selected_summary$union_n,

      selected_retained_fraction =
        selected_summary$retained_fraction,

      selected_remainder_cross_fraction =
        selected_summary$remainder_cross_fraction,

      paper_reference_k =
        paper_k,

      paper_reference_cutoff_rank =
        rank_cutoff_from_k(
          N,
          paper_k
        ),

      paper_joint_n =
        paper_summary$joint_n,

      paper_disjoint_opposite_le_n =
        paper_summary$disjoint_opposite_le_n,

      paper_disjoint_opposite_divergence_n =
        paper_summary$disjoint_opposite_divergence_n,

      paper_permissible_disjoint_n =
        paper_summary$permissible_disjoint_n,

      paper_good_n =
        paper_summary$good_n,

      paper_remainder_cross_n =
        paper_summary$remainder_cross_n,

      paper_union_n =
        paper_summary$union_n,

      paper_retained_fraction =
        paper_summary$retained_fraction,

      paper_remainder_cross_fraction =
        paper_summary$remainder_cross_fraction,

      stringsAsFactors = FALSE
    )
  }

  bind_rows(
    rows
  )
}


build_key_table <- function(
    timepoint_table,
    arm_scaling,
    global_opt,
    global_scan,
    N,
    c1,
    c2,
    paper_k,
    shared_sse) {

  global_selected <- global_scan %>%
    filter(
      k ==
      global_opt$selected_k
    ) %>%
    slice(1L)

  paper_row <- global_scan %>%
    filter(
      k ==
      paper_k
    ) %>%
    slice(1L)

  data.frame(
    n_features = N,

    shared_c1 = c1,
    shared_c2 = c2,

    remainder_size =
      c1 - 1L,

    divergence_interval_size =
      c2 - c1 + 1L,

    leading_edge_size =
      N - c2,

    max_candidate_k =
      N - c2,

    global_selected_k =
      global_opt$selected_k,

    global_selected_cutoff_rank =
      rank_cutoff_from_k(
        N,
        global_opt$selected_k
      ),

    global_selection_method =
      global_opt$selection_method,

    global_breakpoint_remainder_cross_n =
      global_opt$breakpoint_remainder_cross_n,

    global_slope_before =
      global_opt$slope_before,

    global_slope_after =
      global_opt$slope_after,

    global_slope_ratio_after_before =
      global_opt$slope_ratio_after_before,

    global_delta_BIC =
      global_opt$delta_BIC,

    global_shape_supported =
      global_opt$shape_supported,

    global_BIC_supports_segmented =
      global_opt$bic_supports_segmented,

    global_max_zero_remainder_crossing_k =
      global_opt$max_zero_remainder_crossing_k,

    selected_joint_n =
      global_selected$joint_n[1L],

    selected_disjoint_opposite_le_n =
      global_selected$disjoint_opposite_le_n[1L],

    selected_disjoint_opposite_divergence_n =
      global_selected$disjoint_opposite_divergence_n[1L],

    selected_permissible_disjoint_n =
      global_selected$permissible_disjoint_n[1L],

    selected_good_n =
      global_selected$good_n[1L],

    selected_remainder_cross_n =
      global_selected$remainder_cross_n[1L],

    selected_union_n =
      global_selected$union_n[1L],

    selected_retained_fraction =
      global_selected$retained_fraction[1L],

    paper_reference_k =
      paper_k,

    paper_reference_cutoff_rank =
      rank_cutoff_from_k(
        N,
        paper_k
      ),

    paper_reference_used_in_analysis =
      FALSE,

    paper_good_n = if (
      nrow(paper_row) == 1L
    ) {
      paper_row$good_n[1L]
    } else {
      NA_real_
    },

    paper_remainder_cross_n = if (
      nrow(paper_row) == 1L
    ) {
      paper_row$remainder_cross_n[1L]
    } else {
      NA_real_
    },

    pairwise_k_min = min(
      timepoint_table$pairwise_selected_k
    ),

    pairwise_k_median = median(
      timepoint_table$pairwise_selected_k
    ),

    pairwise_k_max = max(
      timepoint_table$pairwise_selected_k
    ),

    anchor_median = median(
      c(
        timepoint_table$control_anchor,
        timepoint_table$treatment_anchor
      )
    ),

    terminal_median = median(
      c(
        timepoint_table$control_terminal,
        timepoint_table$treatment_terminal
      )
    ),

    p_remainder_median =
      safe_summary(
        arm_scaling$p_remainder
      ),

    p_divergence_median =
      safe_summary(
        arm_scaling$p_divergence
      ),

    p_leading_edge_median =
      safe_summary(
        arm_scaling$p_leading_edge
      ),

    p_selected_k_median =
      safe_summary(
        arm_scaling$p_selected_k
      ),

    p_paper_5000_median =
      safe_summary(
        arm_scaling$p_paper_5000
      ),

    shared_model_SSE =
      shared_sse,

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
        z$data[[column]]
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
    group_results[[1L]]$data
  )

  data.frame(
    rank = seq_len(N),

    raw_variance = rankwise_median(
      group_results,
      "display_log1p_raw_empirical_variance"
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


add_rank_regions <- function(
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
      alpha = 0.50
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
      alpha = 0.52
    ) +
    annotate(
      "text",
      x = (1 + c1) / 2,
      y = Inf,
      label = "Remainder",
      vjust = 1.25,
      fontface = "bold",
      size = 3.0
    ) +
    annotate(
      "text",
      x = (c1 + c2) / 2,
      y = Inf,
      label = "Divergence",
      vjust = 1.25,
      fontface = "bold",
      size = 3.0
    ) +
    annotate(
      "text",
      x = (c2 + N) / 2,
      y = Inf,
      label = "Leading edge",
      vjust = 1.25,
      fontface = "bold",
      size = 3.0
    )
}


make_cutoff_panel <- function(
    scan_df,
    selected_k,
    paper_k,
    title_text,
    extra_k = NULL,
    extra_label = "Global k*") {

  selected <- scan_df %>%
    filter(
      k == selected_k
    ) %>%
    slice(1L)

  if (
    nrow(selected) != 1L
  ) {
    stop(
      "Selected k not present in cutoff scan."
    )
  }

  paper <- scan_df %>%
    filter(
      k == paper_k
    ) %>%
    slice(1L)

  frontier <- scan_df %>%
    filter(
      is_pareto
    ) %>%
    arrange(
      remainder_cross_n,
      good_n,
      k
    ) %>%
    distinct(
      remainder_cross_n,
      good_n,
      .keep_all = TRUE
    )

  p <- ggplot() +
    geom_path(
      data = scan_df,
      aes(
        x = remainder_cross_n,
        y = good_n,
        group = 1
      ),
      color = COL$candidate,
      linewidth = 0.48,
      alpha = 0.45
    ) +
    geom_path(
      data = frontier,
      aes(
        x = remainder_cross_n,
        y = good_n,
        color = "Pareto frontier",
        group = 1
      ),
      linewidth = 1.15
    )

  fit_df <- frontier %>%
    filter(
      is.finite(
        segmented_fit_good
      )
    )

  if (
    nrow(fit_df) > 1L
  ) {
    p <- p +
      geom_line(
        data = fit_df,
        aes(
          x = remainder_cross_n,
          y = segmented_fit_good,
          color = "Segmented fit",
          group = 1
        ),
        linewidth = 0.95,
        linetype = "dashed"
      )
  }

  p <- p +
    geom_point(
      data = selected,
      aes(
        x = remainder_cross_n,
        y = good_n,
        color = "Selected k*"
      ),
      shape = 18,
      size = 4.1
    ) +
    annotate(
      "text",
      x = selected$remainder_cross_n,
      y = selected$good_n,
      label = paste0(
        "k*=",
        selected_k
      ),
      vjust = 1.8,
      size = 2.9,
      fontface = "bold"
    )

  if (
    nrow(paper) == 1L
  ) {
    p <- p +
      geom_point(
        data = paper,
        aes(
          x = remainder_cross_n,
          y = good_n,
          color = "5k reference"
        ),
        shape = 1,
        size = 3.6,
        stroke = 1.15
      ) +
      annotate(
        "text",
        x = paper$remainder_cross_n,
        y = paper$good_n,
        label = "5k",
        vjust = -0.85,
        size = 2.8
      )
  }

  if (
    !is.null(extra_k) &&
    is.finite(extra_k) &&
    extra_k != selected_k
  ) {
    extra <- scan_df %>%
      filter(
        k == extra_k
      ) %>%
      slice(1L)

    if (
      nrow(extra) == 1L
    ) {
      p <- p +
        geom_point(
          data = extra,
          aes(
            x = remainder_cross_n,
            y = good_n,
            color = extra_label
          ),
          shape = 4,
          size = 3.5,
          stroke = 1.1
        )
    }
  }

  selected_div <- selected$disjoint_opposite_divergence_n

  selected_rem <- selected$remainder_cross_n

  color_values <- c(
    "Pareto frontier" = COL$pareto,
    "Segmented fit" = COL$segmented,
    "Selected k*" = COL$selected,
    "5k reference" = COL$paper,
    stats::setNames(
      COL$treatment,
      extra_label
    )
  )

  p +
    scale_color_manual(
      values = color_values,
      breaks = intersect(
        names(color_values),
        c(
          "Pareto frontier",
          "Segmented fit",
          "Selected k*",
          "5k reference",
          extra_label
        )
      )
    ) +
    labs(
      title = title_text,
      subtitle = paste0(
        "Divergence-disjoint=",
        selected_div,
        " | remainder crossings=",
        selected_rem
      ),
      x = "Sites crossing into opposite-arm remainder (rank < c1)",
      y = "Joint + permissible disjoint sites"
    ) +
    theme_manuscript()
}


# =============================================================================
# OVERALL MAIN FIGURE
# =============================================================================

make_overall_figure <- function(
    overall_df,
    global_scan,
    global_opt,
    key_table,
    c1,
    c2,
    selected_k,
    paper_k,
    out_file) {

  N <- nrow(
    overall_df
  )

  selected_rank <- rank_cutoff_from_k(
    N,
    selected_k
  )

  paper_rank <- rank_cutoff_from_k(
    N,
    paper_k
  )

  anchor_median <- key_table$anchor_median[1L]

  terminal_median <- key_table$terminal_median[1L]

  # -----------------------------------------------------------------------
  # A. Raw-count variance geometry.
  # -----------------------------------------------------------------------

  a_lines <- data.frame(
    rank = c(
      c1,
      c2,
      selected_rank,
      paper_rank
    ),
    key = c(
      "c1",
      "c2",
      paste0(
        "k*=",
        selected_k
      ),
      "5k ref"
    ),
    stringsAsFactors = FALSE
  )

  selected_key <- paste0(
    "k*=",
    selected_k
  )

  a_colors <- c(
    "Raw variance" = COL$raw,
    "c1" = COL$c1,
    "c2" = COL$c2,
    stats::setNames(
      COL$selected,
      selected_key
    ),
    "5k ref" = COL$paper
  )

  a_types <- c(
    "Raw variance" = "solid",
    "c1" = "dashed",
    "c2" = "longdash",
    stats::setNames(
      "dotdash",
      selected_key
    ),
    "5k ref" = "dotted"
  )

  pA <- ggplot()

  pA <- add_rank_regions(
    pA,
    c1,
    c2,
    N
  )

  pA <- pA +
    geom_vline(
      data = a_lines,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      linewidth = 0.85
    ) +
    geom_vline(
      xintercept = anchor_median,
      color = COL$divergence,
      linewidth = 0.68,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = terminal_median,
      color = COL$divergence,
      linetype = "dotted",
      linewidth = 0.75,
      show.legend = FALSE
    ) +
    geom_line(
      data = overall_df,
      aes(
        x = rank,
        y = raw_variance,
        color = "Raw variance",
        linetype = "Raw variance"
      ),
      linewidth = 1.20
    ) +
    annotate(
      "text",
      x = anchor_median,
      y = -Inf,
      label = "A",
      angle = 90,
      vjust = -0.30,
      size = 2.8
    ) +
    annotate(
      "text",
      x = terminal_median,
      y = -Inf,
      label = "T",
      angle = 90,
      vjust = -0.30,
      size = 2.8
    ) +
    scale_color_manual(
      values = a_colors,
      breaks = names(
        a_colors
      )
    ) +
    scale_linetype_manual(
      values = a_types,
      breaks = names(
        a_types
      )
    ) +
    labs(
      title = "A. Raw-count variance geometry",
      subtitle = "A/T = first two variance-curvature crossings after c2",
      x = "PC1 rank: low |loading| -> high |loading|",
      y = "Smoothed log(1 + raw-count variance)"
    ) +
    theme_manuscript()

  # -----------------------------------------------------------------------
  # B. Cumulative divergence.
  # -----------------------------------------------------------------------

  b_lines <- data.frame(
    rank = c(
      c1,
      c2,
      selected_rank,
      paper_rank
    ),
    key = c(
      "c1",
      "c2",
      selected_key,
      "5k ref"
    ),
    stringsAsFactors = FALSE
  )

  b_colors <- c(
    "D(r)" = COL$divergence,
    "Shared fit" = COL$fit,
    "c1" = COL$c1,
    "c2" = COL$c2,
    stats::setNames(
      COL$selected,
      selected_key
    ),
    "5k ref" = COL$paper
  )

  b_types <- c(
    "D(r)" = "solid",
    "Shared fit" = "solid",
    "c1" = "dashed",
    "c2" = "longdash",
    stats::setNames(
      "dotdash",
      selected_key
    ),
    "5k ref" = "dotted"
  )

  pB <- ggplot()

  pB <- add_rank_regions(
    pB,
    c1,
    c2,
    N
  )

  pB <- pB +
    geom_vline(
      data = b_lines,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      linewidth = 0.85
    ) +
    geom_hline(
      yintercept = 0,
      color = "#777777",
      linetype = "dotted",
      linewidth = 0.35
    ) +
    geom_line(
      data = overall_df,
      aes(
        x = rank,
        y = D,
        color = "D(r)",
        linetype = "D(r)"
      ),
      linewidth = 0.95,
      alpha = 0.72
    ) +
    geom_line(
      data = overall_df,
      aes(
        x = rank,
        y = D_fit,
        color = "Shared fit",
        linetype = "Shared fit"
      ),
      linewidth = 1.30
    ) +
    scale_color_manual(
      values = b_colors,
      breaks = names(
        b_colors
      )
    ) +
    scale_linetype_manual(
      values = b_types,
      breaks = names(
        b_types
      )
    ) +
    labs(
      title = "B. Cumulative PC1-NB variance-mass divergence",
      subtitle = "Shared c1/c2 fitted across all 8 arms",
      x = "PC1 rank",
      y = "D(r) = F_E(r) - F_P(r)"
    ) +
    theme_manuscript()

  # -----------------------------------------------------------------------
  # C. c1-based cutoff optimization.
  # -----------------------------------------------------------------------

  pC <- make_cutoff_panel(
    scan_df = global_scan,
    selected_k = selected_k,
    paper_k = paper_k,
    title_text = "C. Global top-k optimization"
  ) +
    labs(
      subtitle = paste0(
        "c1 is the hard opposite-arm boundary | method: ",
        global_opt$selection_method
      )
    )

  save_panels(
    list(
      pA,
      pB,
      pC
    ),
    out_file,
    height_in = 11.8
  )
}


# =============================================================================
# TIME-POINT MAIN FIGURES
# =============================================================================

make_timepoint_figure <- function(
    comparison_name,
    mapping,
    group_results,
    pair_scan,
    pair_optimum,
    global_k,
    c1,
    c2,
    paper_k,
    out_file) {

  control_group <- unname(
    mapping[["control"]]
  )

  treatment_group <- unname(
    mapping[["treatment"]]
  )

  control <- group_results[[
    control_group
  ]]$data

  treatment <- group_results[[
    treatment_group
  ]]$data

  N <- nrow(
    control
  )

  selected_rank <- rank_cutoff_from_k(
    N,
    global_k
  )

  paper_rank <- rank_cutoff_from_k(
    N,
    paper_k
  )

  c_anchor <- group_results[[
    control_group
  ]]$anchor

  c_terminal <- group_results[[
    control_group
  ]]$terminal

  t_anchor <- group_results[[
    treatment_group
  ]]$anchor

  t_terminal <- group_results[[
    treatment_group
  ]]$terminal

  selected_key <- paste0(
    "k*=",
    global_k
  )

  # -----------------------------------------------------------------------
  # A. Raw-count variance geometry.
  # -----------------------------------------------------------------------

  a_lines <- data.frame(
    rank = c(
      c1,
      c2,
      selected_rank,
      paper_rank
    ),
    key = c(
      "c1",
      "c2",
      selected_key,
      "5k ref"
    ),
    stringsAsFactors = FALSE
  )

  a_colors <- c(
    "Control" = COL$control,
    "Treatment" = COL$treatment,
    "c1" = COL$c1,
    "c2" = COL$c2,
    stats::setNames(
      COL$selected,
      selected_key
    ),
    "5k ref" = COL$paper
  )

  a_types <- c(
    "Control" = "solid",
    "Treatment" = "solid",
    "c1" = "dashed",
    "c2" = "longdash",
    stats::setNames(
      "dotdash",
      selected_key
    ),
    "5k ref" = "dotted"
  )

  pA <- ggplot()

  pA <- add_rank_regions(
    pA,
    c1,
    c2,
    N
  )

  pA <- pA +
    geom_vline(
      data = a_lines,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      linewidth = 0.84
    ) +
    geom_vline(
      xintercept = c_anchor,
      color = COL$anchor_control,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = c_terminal,
      color = COL$terminal_control,
      linetype = "dotted",
      linewidth = 0.75,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = t_anchor,
      color = COL$anchor_treatment,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = t_terminal,
      color = COL$terminal_treatment,
      linetype = "dotted",
      linewidth = 0.75,
      show.legend = FALSE
    ) +
    geom_line(
      data = control,
      aes(
        x = rank,
        y = display_log1p_raw_empirical_variance,
        color = "Control",
        linetype = "Control"
      ),
      linewidth = 1.10
    ) +
    geom_line(
      data = treatment,
      aes(
        x = rank,
        y = display_log1p_raw_empirical_variance,
        color = "Treatment",
        linetype = "Treatment"
      ),
      linewidth = 1.10
    ) +
    annotate(
      "text",
      x = c_anchor,
      y = -Inf,
      label = "A_C",
      angle = 90,
      vjust = -0.25,
      size = 2.5,
      color = COL$control
    ) +
    annotate(
      "text",
      x = c_terminal,
      y = -Inf,
      label = "T_C",
      angle = 90,
      vjust = -0.25,
      size = 2.5,
      color = COL$control
    ) +
    annotate(
      "text",
      x = t_anchor,
      y = -Inf,
      label = "A_T",
      angle = 90,
      vjust = -0.25,
      size = 2.5,
      color = COL$treatment
    ) +
    annotate(
      "text",
      x = t_terminal,
      y = -Inf,
      label = "T_T",
      angle = 90,
      vjust = -0.25,
      size = 2.5,
      color = COL$treatment
    ) +
    scale_color_manual(
      values = a_colors,
      breaks = names(
        a_colors
      )
    ) +
    scale_linetype_manual(
      values = a_types,
      breaks = names(
        a_types
      )
    ) +
    labs(
      title = paste0(
        "A. ",
        comparison_name,
        " raw-count variance geometry"
      ),
      subtitle = "A/T = first two variance-curvature crossings after c2",
      x = "PC1 rank: low |loading| -> high |loading|",
      y = "Smoothed log(1 + raw-count variance)"
    ) +
    theme_manuscript()

  # -----------------------------------------------------------------------
  # B. Cumulative divergence.
  # -----------------------------------------------------------------------

  b_lines <- data.frame(
    rank = c(
      c1,
      c2,
      selected_rank,
      paper_rank
    ),
    key = c(
      "c1",
      "c2",
      selected_key,
      "5k ref"
    ),
    stringsAsFactors = FALSE
  )

  b_colors <- c(
    "Control D(r)" = COL$control,
    "Treatment D(r)" = COL$treatment,
    "c1" = COL$c1,
    "c2" = COL$c2,
    stats::setNames(
      COL$selected,
      selected_key
    ),
    "5k ref" = COL$paper
  )

  b_types <- c(
    "Control D(r)" = "solid",
    "Treatment D(r)" = "solid",
    "c1" = "dashed",
    "c2" = "longdash",
    stats::setNames(
      "dotdash",
      selected_key
    ),
    "5k ref" = "dotted"
  )

  pB <- ggplot()

  pB <- add_rank_regions(
    pB,
    c1,
    c2,
    N
  )

  pB <- pB +
    geom_vline(
      data = b_lines,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      linewidth = 0.84
    ) +
    geom_hline(
      yintercept = 0,
      color = "#777777",
      linetype = "dotted",
      linewidth = 0.35
    ) +
    geom_line(
      data = control,
      aes(
        x = rank,
        y = display_D,
        color = "Control D(r)",
        linetype = "Control D(r)"
      ),
      linewidth = 1.08
    ) +
    geom_line(
      data = treatment,
      aes(
        x = rank,
        y = display_D,
        color = "Treatment D(r)",
        linetype = "Treatment D(r)"
      ),
      linewidth = 1.08
    ) +
    scale_color_manual(
      values = b_colors,
      breaks = names(
        b_colors
      )
    ) +
    scale_linetype_manual(
      values = b_types,
      breaks = names(
        b_types
      )
    ) +
    labs(
      title = paste0(
        "B. ",
        comparison_name,
        " cumulative PC1-NB divergence"
      ),
      subtitle = "Shared c1/c2",
      x = "PC1 rank",
      y = "D(r) = F_E(r) - F_P(r)"
    ) +
    theme_manuscript()

  # -----------------------------------------------------------------------
  # C. Pair-specific c1-based cutoff optimization.
  # -----------------------------------------------------------------------

  pC <- make_cutoff_panel(
    scan_df = pair_scan,
    selected_k =
      pair_optimum$selected_k,
    paper_k = paper_k,
    title_text = paste0(
      "C. ",
      comparison_name,
      " top-k optimization"
    ),
    extra_k = global_k,
    extra_label = "Global k*"
  ) +
    labs(
      subtitle = paste0(
        "Pair k*=",
        pair_optimum$selected_k,
        " | global k*=",
        global_k
      )
    )

  save_panels(
    list(
      pA,
      pB,
      pC
    ),
    out_file,
    height_in = 11.8
  )
}


# =============================================================================
# OPTIONAL SUPPLEMENT: NB MEAN-VARIANCE SCALING
# =============================================================================

make_nb_scaling_supplement <- function(
    arm_scaling,
    selected_k,
    out_file) {

  scaling <- arm_scaling %>%
    mutate(
      arm_label = ifelse(
        arm == "control",
        "Control",
        "Treatment"
      )
    ) %>%
    select(
      comparison,
      arm_label,
      p_remainder,
      p_divergence,
      p_leading_edge,
      p_selected_k
    ) %>%
    pivot_longer(
      cols = c(
        p_remainder,
        p_divergence,
        p_leading_edge,
        p_selected_k
      ),
      names_to = "region",
      values_to = "p"
    ) %>%
    mutate(
      region = factor(
        region,
        levels = c(
          "p_remainder",
          "p_divergence",
          "p_leading_edge",
          "p_selected_k"
        ),
        labels = c(
          "Remainder",
          "Divergence",
          "Leading edge",
          paste0(
            "Selected k*=",
            selected_k
          )
        )
      )
    )

  p <- ggplot(
    scaling,
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
      linewidth = 0.55
    ) +
    geom_vline(
      xintercept = 2,
      color = "#777777",
      linetype = "dotted",
      linewidth = 0.65
    ) +
    geom_point(
      size = 3.1,
      alpha = 0.88,
      position = position_dodge(
        width = 0.20
      )
    ) +
    facet_wrap(
      ~ comparison,
      ncol = 2
    ) +
    scale_color_manual(
      values = c(
        "Control" = COL$control,
        "Treatment" = COL$treatment
      )
    ) +
    labs(
      title = "Supplement. NB mean-variance scaling",
      subtitle = "Post-boundary corroboration only",
      x = "Exponent p in E = alpha * mu^p",
      y = NULL
    ) +
    theme_manuscript(
      base_size = 11
    ) +
    theme(
      strip.text = element_text(
        face = "bold"
      )
    )

  grDevices::png(
    filename = out_file,
    width = 12,
    height = 8.5,
    units = "in",
    res = PNG_DPI,
    bg = "white"
  )

  print(p)

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

N <- nrow(
  count_mat
)

if (
  PAPER_REFERENCE_K >= N
) {
  stop(
    "PAPER_REFERENCE_K must be smaller than N."
  )
}

PAPER_REFERENCE_RANK <- rank_cutoff_from_k(
  N,
  PAPER_REFERENCE_K
)

# -------------------------------------------------------------------------
# Global DESeq2 normalization and pooled within-group empirical variance.
# -------------------------------------------------------------------------

deseq <- normalize_deseq2_global(
  count_mat = count_mat,
  group_labels = group_labels
)

normalized_counts <- deseq$normalized_counts

pooled <- compute_pooled_within_group_variance(
  normalized_counts =
    normalized_counts,
  group_labels =
    group_labels
)

message(
  "Pooled within-group residual df: ",
  pooled$residual_df
)

# -------------------------------------------------------------------------
# Arm-specific PC1 ranking, raw variance geometry, and cumulative divergence.
# -------------------------------------------------------------------------

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

for (g in levels(
  group_labels
)) {
  idx <- which(
    group_labels ==
    g
  )

  if (
    length(idx) < 2L
  ) {
    stop(
      "Not enough samples in group ",
      g
    )
  }

  message(
    "Analyzing ",
    g,
    "..."
  )

  group_results[[g]] <-
    compute_group_analysis(
      group_name = g,

      raw_counts_arm =
        count_mat[
          ,
          idx,
          drop = FALSE
        ],

      normalized_counts_arm =
        normalized_counts[
          ,
          idx,
          drop = FALSE
        ],

      pooled_variance =
        pooled$variance
    )
}

# -------------------------------------------------------------------------
# Shared c1/c2 cumulative-divergence model.
# -------------------------------------------------------------------------

message(
  "Fitting shared two-knot cumulative-divergence model..."
)

knot_fit <- fit_shared_knots(
  group_results
)

C1 <- knot_fit$c1
C2 <- knot_fit$c2

MAX_CANDIDATE_K <- as.integer(
  N - C2
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
  "Own-arm leading-edge candidate limit: k <= ",
  MAX_CANDIDATE_K
)

if (
  PAPER_REFERENCE_K >
  MAX_CANDIDATE_K
) {
  warning(
    "The historical top-5,000 reference extends left of c2 and therefore ",
    "falls outside the own-arm leading-edge candidate domain. It remains ",
    "reference-only."
  )
}

# -------------------------------------------------------------------------
# Anchor / Terminal: first two raw-count variance-curvature crossings after c2.
# -------------------------------------------------------------------------

message(
  "Selecting Anchor/Terminal from first two raw-variance y''=0 crossings after c2..."
)

for (g in names(
  group_results
)) {
  at <- select_anchor_terminal_after_c2(
    crossings =
      group_results[[
        g
      ]]$raw_variance_crossings,

    c2 = C2,

    total_n = N
  )

  group_results[[g]]$anchor <-
    at$anchor

  group_results[[g]]$terminal <-
    at$terminal
}

# -------------------------------------------------------------------------
# Pairwise top-k scans.
#
# IMPORTANT:
# - selection originates only from r > c2;
# - opposite-arm divergence c1..c2 is permissible;
# - only opposite-arm remainder r < c1 is penalized.
# -------------------------------------------------------------------------

message(
  "Scanning candidate k values using c1 as the opposite-arm contamination boundary..."
)

pair_scans <- list()
pair_optima <- list()

for (comparison_name in names(
  COMPARISONS
)) {
  mapping <- COMPARISONS[[
    comparison_name
  ]]

  control_group <- unname(
    mapping[["control"]]
  )

  treatment_group <- unname(
    mapping[["treatment"]]
  )

  raw_scan <- scan_pair_cutoffs(
    control_df =
      group_results[[
        control_group
      ]]$data,

    treatment_df =
      group_results[[
        treatment_group
      ]]$data,

    c1 = C1,
    c2 = C2,

    comparison_name =
      comparison_name,

    control_group =
      control_group,

    treatment_group =
      treatment_group
  )

  opt <- select_pareto_changepoint(
    raw_scan,
    good_col = "good_n",
    cost_col = "remainder_cross_n"
  )

  pair_scans[[
    comparison_name
  ]] <- opt$scan

  pair_optima[[
    comparison_name
  ]] <- opt

  selected_row <- opt$scan %>%
    filter(
      k ==
      opt$selected_k
    ) %>%
    slice(1L)

  message(
    comparison_name,
    ": k*=",
    opt$selected_k,
    "; method=",
    opt$selection_method,
    "; retained=",
    selected_row$good_n[1L],
    "; divergence-disjoint=",
    selected_row$disjoint_opposite_divergence_n[1L],
    "; remainder crossings=",
    selected_row$remainder_cross_n[1L],
    "; Delta_BIC=",
    ifelse(
      is.finite(
        opt$delta_BIC
      ),
      signif(
        opt$delta_BIC,
        5
      ),
      "NA"
    )
  )
}

# -------------------------------------------------------------------------
# Global common cutoff from pooled pairwise counts.
# -------------------------------------------------------------------------

global_scan_raw <- aggregate_global_cutoff_scan(
  pair_scans
)

global_opt <- select_pareto_changepoint(
  global_scan_raw,
  good_col = "good_n",
  cost_col = "remainder_cross_n"
)

global_scan <- global_opt$scan

GLOBAL_K <- global_opt$selected_k

GLOBAL_CUTOFF_RANK <- rank_cutoff_from_k(
  N,
  GLOBAL_K
)

global_selected_row <- global_scan %>%
  filter(
    k ==
    GLOBAL_K
  ) %>%
  slice(1L)

message(
  "GLOBAL k*=",
  GLOBAL_K,
  " (cutoff rank ",
  GLOBAL_CUTOFF_RANK,
  ")"
)

message(
  "Global method=",
  global_opt$selection_method,
  "; retained=",
  global_selected_row$good_n[1L],
  "; divergence-disjoint=",
  global_selected_row$disjoint_opposite_divergence_n[1L],
  "; remainder crossings=",
  global_selected_row$remainder_cross_n[1L]
)

message(
  "Historical paper k=",
  PAPER_REFERENCE_K,
  " (cutoff rank ",
  PAPER_REFERENCE_RANK,
  "; reference only)"
)

# -------------------------------------------------------------------------
# Classify actual sites at GLOBAL_K.
#
# Retained:
#   Joint
#   + Disjoint / opposite Leading Edge
#   + Disjoint / opposite Divergence
#
# Excluded:
#   Disjoint / opposite Remainder
# -------------------------------------------------------------------------

selected_rows <- list()
excluded_rows <- list()

selected_summaries <- list()
paper_summaries <- list()

make_na_summary <- function(
    k,
    cutoff_rank) {

  data.frame(
    k = k,
    cutoff_rank = cutoff_rank,

    joint_n = NA_integer_,

    disjoint_control_opposite_le_n =
      NA_integer_,

    disjoint_treatment_opposite_le_n =
      NA_integer_,

    disjoint_opposite_le_n =
      NA_integer_,

    disjoint_control_opposite_divergence_n =
      NA_integer_,

    disjoint_treatment_opposite_divergence_n =
      NA_integer_,

    disjoint_opposite_divergence_n =
      NA_integer_,

    permissible_disjoint_control_n =
      NA_integer_,

    permissible_disjoint_treatment_n =
      NA_integer_,

    permissible_disjoint_n =
      NA_integer_,

    good_n = NA_integer_,

    remainder_cross_n =
      NA_integer_,

    union_n = NA_integer_,

    retained_fraction =
      NA_real_,

    remainder_cross_fraction =
      NA_real_,

    stringsAsFactors = FALSE
  )
}

for (comparison_name in names(
  COMPARISONS
)) {
  mapping <- COMPARISONS[[
    comparison_name
  ]]

  control_group <- unname(
    mapping[["control"]]
  )

  treatment_group <- unname(
    mapping[["treatment"]]
  )

  selected_all <- classify_pair_at_k(
    control_df =
      group_results[[
        control_group
      ]]$data,

    treatment_df =
      group_results[[
        treatment_group
      ]]$data,

    k = GLOBAL_K,

    c1 = C1,
    c2 = C2,

    comparison_name =
      comparison_name,

    control_group =
      control_group,

    treatment_group =
      treatment_group
  )

  selected_summaries[[
    comparison_name
  ]] <- summarize_classification(
    selected_all
  )

  selected_rows[[
    comparison_name
  ]] <- selected_all %>%
    filter(
      retained_for_analysis
    )

  excluded_rows[[
    comparison_name
  ]] <- selected_all %>%
    filter(
      cross_into_remainder
    )

  if (
    PAPER_REFERENCE_K <=
    MAX_CANDIDATE_K
  ) {
    paper_all <- classify_pair_at_k(
      control_df =
        group_results[[
          control_group
        ]]$data,

      treatment_df =
        group_results[[
          treatment_group
        ]]$data,

      k = PAPER_REFERENCE_K,

      c1 = C1,
      c2 = C2,

      comparison_name =
        comparison_name,

      control_group =
        control_group,

      treatment_group =
        treatment_group
    )

    paper_summaries[[
      comparison_name
    ]] <- summarize_classification(
      paper_all
    )
  } else {
    paper_summaries[[
      comparison_name
    ]] <- make_na_summary(
      PAPER_REFERENCE_K,
      PAPER_REFERENCE_RANK
    )
  }
}

selected_sites <- bind_rows(
  selected_rows
)

excluded_sites <- bind_rows(
  excluded_rows
)

write.csv(
  selected_sites,
  file.path(
    OUT_ROOT,
    "Selected_LeadingEdge_Sites.csv"
  ),
  row.names = FALSE
)

write.csv(
  excluded_sites,
  file.path(
    OUT_ROOT,
    "Excluded_Remainder_Crossing_Sites.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Supplemental NB scaling.
# This remains entirely post-boundary and cannot alter c1, c2, or k*.
# -------------------------------------------------------------------------

arm_scaling <- build_arm_scaling_table(
  group_results =
    group_results,

  comparisons =
    COMPARISONS,

  c1 = C1,
  c2 = C2,

  selected_k =
    GLOBAL_K,

  paper_k =
    PAPER_REFERENCE_K
)

write.csv(
  arm_scaling,
  file.path(
    OUT_ROOT,
    "Table_NB_Scaling.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Tables.
# -------------------------------------------------------------------------

timepoint_table <- build_timepoint_table(
  group_results =
    group_results,

  comparisons =
    COMPARISONS,

  pair_optima =
    pair_optima,

  selected_summaries =
    selected_summaries,

  paper_summaries =
    paper_summaries,

  global_k =
    GLOBAL_K,

  global_opt =
    global_opt,

  c1 = C1,
  c2 = C2,

  paper_k =
    PAPER_REFERENCE_K,

  N = N
)

key_table <- build_key_table(
  timepoint_table =
    timepoint_table,

  arm_scaling =
    arm_scaling,

  global_opt =
    global_opt,

  global_scan =
    global_scan,

  N = N,

  c1 = C1,
  c2 = C2,

  paper_k =
    PAPER_REFERENCE_K,

  shared_sse =
    knot_fit$SSE
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

cutoff_table <- bind_rows(
  bind_rows(
    pair_scans
  ),
  global_scan
) %>%
  arrange(
    comparison,
    k
  )

write.csv(
  cutoff_table,
  file.path(
    OUT_ROOT,
    "Table_Cutoff_Optimization.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Figures.
# -------------------------------------------------------------------------

overall_df <- build_overall_figure_data(
  group_results =
    group_results,

  knot_fit =
    knot_fit
)

figure_paths <- character(0)

overall_path <- file.path(
  FIG_DIR,
  "Figure_Overall.png"
)

make_overall_figure(
  overall_df =
    overall_df,

  global_scan =
    global_scan,

  global_opt =
    global_opt,

  key_table =
    key_table,

  c1 = C1,
  c2 = C2,

  selected_k =
    GLOBAL_K,

  paper_k =
    PAPER_REFERENCE_K,

  out_file =
    overall_path
)

figure_paths <- c(
  figure_paths,
  overall_path
)

for (comparison_name in names(
  COMPARISONS
)) {
  fig_path <- file.path(
    FIG_DIR,
    paste0(
      "Figure_",
      comparison_name,
      ".png"
    )
  )

  make_timepoint_figure(
    comparison_name =
      comparison_name,

    mapping =
      COMPARISONS[[
        comparison_name
      ]],

    group_results =
      group_results,

    pair_scan =
      pair_scans[[
        comparison_name
      ]],

    pair_optimum =
      pair_optima[[
        comparison_name
      ]],

    global_k =
      GLOBAL_K,

    c1 = C1,
    c2 = C2,

    paper_k =
      PAPER_REFERENCE_K,

    out_file =
      fig_path
  )

  figure_paths <- c(
    figure_paths,
    fig_path
  )
}

if (
  isTRUE(
    INCLUDE_NB_SCALING_SUPPLEMENT
  )
) {
  nb_supp_path <- file.path(
    FIG_DIR,
    "Figure_Supplement_NB_Scaling.png"
  )

  make_nb_scaling_supplement(
    arm_scaling =
      arm_scaling,

    selected_k =
      GLOBAL_K,

    out_file =
      nb_supp_path
  )

  figure_paths <- c(
    figure_paths,
    nb_supp_path
  )
}

# -------------------------------------------------------------------------
# Zip figures while retaining individual PNGs.
# -------------------------------------------------------------------------

ZIP_PATH <- file.path(
  OUT_ROOT,
  "Figures_All.zip"
)

if (
  file.exists(
    ZIP_PATH
  )
) {
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
      zipfile =
        ZIP_PATH,

      files =
        basename(
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

if (
  !zip_ok
) {
  warning(
    "Figure PNGs were created, but Figures_All.zip was not created."
  )
}

# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

message(
  "============================================================"
)

message(
  "FINAL c1-CONSTRAINED PC1-NB CUTOFF ANALYSIS COMPLETE"
)

message(
  "Shared c1 = ",
  C1,
  " | shared c2 = ",
  C2
)

message(
  "Remainder: rank < ",
  C1,
  " | Divergence: ",
  C1,
  "-",
  C2,
  " | Leading edge: rank > ",
  C2
)

message(
  "Candidate k domain: 1-",
  MAX_CANDIDATE_K
)

message(
  "Global k* = ",
  GLOBAL_K,
  " | method = ",
  global_opt$selection_method
)

message(
  "At global k*: retained=",
  global_selected_row$good_n[1L],
  " | divergence-disjoint=",
  global_selected_row$disjoint_opposite_divergence_n[1L],
  " | remainder crossings=",
  global_selected_row$remainder_cross_n[1L]
)

message(
  "Pair-specific k*: ",
  paste(
    names(
      pair_optima
    ),
    vapply(
      pair_optima,
      function(z) {
        z$selected_k
      },
      integer(1)
    ),
    sep = "=",
    collapse = "; "
  )
)

message(
  "Historical 5k reference used in fitting: FALSE"
)

message(
  "Primary figures: ",
  FIG_DIR
)

message(
  "Selected sites: Selected_LeadingEdge_Sites.csv"
)

message(
  "Excluded sites: Excluded_Remainder_Crossing_Sites.csv"
)

message(
  "Tables: Table_Key_Results.csv; Table_Timepoints.csv; ",
  "Table_Cutoff_Optimization.csv; Table_NB_Scaling.csv"
)

message(
  "============================================================"
)
