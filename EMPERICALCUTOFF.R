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
# PC1-NB DIVERGENCE + ORIGINAL-STYLE LEADING-EDGE CORROBORATION
# =============================================================================
#
# METHODS SUMMARY
# ---------------
#
# 1) PC1 ranking
#    For each arm, counts are normalized to CPM, transformed as log(1 + CPM),
#    and PCA is performed across samples. Features are ranked by ascending
#    absolute PC1 loading:
#
#        rank_i = rank(|v_i1|), ascending.
#
#    The x-axis therefore runs from low PC1 contribution to high PC1
#    contribution. The feature-level variance contribution for PC1 is
#
#        P_i = lambda_1 * v_i1^2
#
#    where lambda_1 is the PC1 eigenvalue.
#
#
# 2) Two variance quantities are used for two different purposes
#
#    A. Within-arm raw-count empirical variance:
#
#        s_raw,ig^2 = Var_j(Y_ij | arm g)
#
#       This is used for the original raw-variance geometry and for the
#       Anchor/Terminal markers.
#
#    B. Pooled within-group empirical variance of DESeq2-normalized counts:
#
#                          sum_g sum_{j in g} (y_ij - ybar_ig)^2
#        V_pool,i =         -------------------------------------
#                                   sum_g (n_g - 1)
#
#       where y_ij is the DESeq2 size-factor-normalized count. This is an
#       empirical pooled variance and is NOT a DESeq2 fitted dispersion.
#
#       For arm g,
#
#        mu_ig = mean_j(y_ij | arm g)
#        E_ig  = max(V_pool,i - mu_ig, 0)
#
#       This E_ig is the empirical excess-over-Poisson variance signal used in
#       the PC1-NB divergence analysis.
#
#
# 3) Shared PC1-NB divergence boundaries
#
#    Within each arm, convert PC1 variance contribution and excess variance into
#    rank-wise masses:
#
#        p_g(r) = P_g(r) / sum_j P_g(j)
#        q_g(r) = E_g(r) / sum_j E_g(j)
#
#    Their cumulative masses are
#
#        F_P,g(r) = sum_{j <= r} p_g(j)
#        F_E,g(r) = sum_{j <= r} q_g(j)
#
#    and cumulative divergence is
#
#        D_g(r) = F_E,g(r) - F_P,g(r).
#
#    A continuous two-knot linear-spline model is fit jointly across all
#    eight arms:
#
#        D_g(x) = beta_0g + beta_1g*x
#                 + gamma_1g*(x-c1)_+
#                 + gamma_2g*(x-c2)_+
#
#    with arm-specific coefficients but shared c1 and c2. The shared c2 marks
#    the onset of the broad leading-edge regime.
#
#
# 4) Anchor and Terminal inside the leading-edge regime
#
#    The original raw-count variance geometry is then applied within the
#    independently defined leading-edge regime. For each arm, features are
#    ordered by the PC1 rank above, and the within-arm raw-count variance curve
#    is transformed as
#
#        y_g(r) = log(1 + s_raw,g^2(r)).
#
#    A smoothing spline (spar = 0.60) is fit, and its second derivative is
#    computed:
#
#        y_g''(r).
#
#    Starting at c2 and moving rightward:
#
#        Anchor_g   = first sign-change zero-crossing of y_g''(r) after c2
#        Terminal_g = second successive sign-change zero-crossing of y_g''(r)
#
#    The historical 5,000-feature cutoff is retained only as a reference:
#
#        Ref = N - 5000 + 1
#
#    and does NOT influence c1, c2, Anchor, or Terminal.
#
#
# 5) Original-style corroboration
#
#    Once Anchor is identified, the original corroboration structure is
#    preserved:
#
#        RIGHT = [Anchor, N]
#        LEFT  = matched equal-sized block immediately left of Anchor
#
#    using the following raw-count feature-level summaries:
#
#        mu_raw = mean raw count
#        var_raw = raw-count empirical variance
#        NB2 = log(1 + max(var_raw - mu_raw, 0))
#        NB2-NB1 = log(1 + max(var_raw - mu_raw, 0)) - log(1 + mu_raw)
#        alpha*mu = log(1 + alpha_hat * mu_raw)
#
#    where
#
#        alpha_hat = max((var_raw - mu_raw) / mu_raw^2, 0).
#
#    These NB2 / NB2-NB1 / alpha*mu quantities are descriptive moment-based
#    corroboration metrics. They are not formal NB1-versus-NB2 likelihood
#    tests and they are not used to select any boundary.
#
#
# OUTPUTS
# -------
# Figures are written individually to:
#   OUT_ROOT/Figures/
#
# and also collected into:
#   OUT_ROOT/Figures_All.zip
#
# Figure files:
#   Figure_Overall_Divergence.png
#       Overall mathematical explanation: raw-count variance geometry,
#       cumulative PC1-NB divergence, shared c1/c2, Anchor/Terminal, and
#       the 5,000-feature manuscript reference.
#
#   Figure_New_<TIMEPOINT>.png
#       Time-point-specific raw-count variance, cumulative divergence, and
#       clearly labeled post-boundary NB scaling corroboration.
#
#   Figure_Original_<TIMEPOINT>.png
#       Rebuilt original-style geometry, matched LEFT/RIGHT raw-count
#       corroboration, and compact LEFT-versus-RIGHT summary.
#
# Every panel explicitly labels curve meaning, boundary meaning, variance
# source, region shading, and the 5,000-feature cutoff as REFERENCE ONLY.
#
# Tables:
#   Table_Key_Results.csv
#   Table_Timepoints.csv
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT   <- "/root/REAPER98632/exports/final_manuscript_dualfig"
FIG_DIR    <- file.path(OUT_ROOT, "Figures")

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

PAPER_LEADING_EDGE_SIZE <- 5000L
ANALYSIS_SPLINE_SPAR    <- 0.60

# Display-only smoothing for cleaner figures.
DISPLAY_VAR_SPAR   <- 0.72
DISPLAY_MASS_SPAR  <- 0.72
DISPLAY_DIV_SPAR   <- 0.72

PNG_WIDTH_IN <- 15
PNG_DPI      <- 360


# =============================================================================
# COLORS
# =============================================================================

COL <- list(
  control = "#386CB0",
  treatment = "#159D91",
  empirical = "#117A65",
  divergence = "#222222",
  fit = "#000000",
  pc1 = "#386CB0",
  nb = "#159D91",

  remainder = "#DCEFF2",
  interval = "#F5E8C8",
  leading = "#DDF2EA",

  left_fill = "#CBE3F8",
  right_fill = "#DDF2D5",
  interval_fill = "#A6A6A6",

  c1 = "#2166AC",
  c2 = "#1B7837",
  ref = "#E69F00",
  anchor = "#000000",
  terminal = "#D95F02",

  nb2 = "#1B9E77",
  nbgap = "#CC1E8C",
  alphamu = "#386CB0"
)


# =============================================================================
# HELPERS
# =============================================================================

theme_manuscript <- function(base_size = 11.3) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.2),
      plot.subtitle = element_text(size = base_size - 0.4, margin = margin(b = 5)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#333333"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.border = element_rect(color = "#B7B7B7", fill = NA, linewidth = 0.45),
      panel.grid = element_blank(),
      plot.margin = margin(7, 9, 7, 9)
    )
}

save_panels <- function(plots, path, height_in) {
  grDevices::png(path, width = PNG_WIDTH_IN, height = height_in, units = "in", res = PNG_DPI, bg = "white")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = length(plots), ncol = 1L)))
  for (i in seq_along(plots)) {
    print(plots[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1L))
  }
  dev.off()
}

read_count_matrix <- function(path, group_patterns) {
  raw_df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  sample_idx <- sort(unique(unlist(lapply(group_patterns, function(p) grep(p, colnames(raw_df))))))
  if (length(sample_idx) == 0L) stop("No sample columns matched GROUP_PATTERNS.")

  count_df <- raw_df[, sample_idx, drop = FALSE]
  count_mat <- do.call(cbind, lapply(count_df, function(x) suppressWarnings(as.numeric(trimws(as.character(x))))))
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

  if (nrow(count_mat) < 2L) stop("Fewer than two nonzero features remain.")

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
    stop("Unassigned samples: ", paste(sample_names[is.na(assigned)], collapse = ", "))
  }

  factor(assigned, levels = names(group_patterns))
}

normalize_cpm_log1p <- function(count_mat_arm) {
  lib_size <- colSums(count_mat_arm, na.rm = TRUE)
  lib_size[!is.finite(lib_size) | lib_size <= 0] <- 1
  cpm <- sweep(count_mat_arm, 2L, lib_size / 1e6, "/")
  log1p(cpm)
}

normalize_deseq2_global <- function(count_mat, group_labels) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required for the pooled normalized-count variance.")
  }

  col_data <- data.frame(group = factor(group_labels), row.names = colnames(count_mat))

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = col_data,
    design = ~ group
  )

  dds <- tryCatch(
    DESeq2::estimateSizeFactors(dds),
    error = function(e) DESeq2::estimateSizeFactors(dds, type = "poscounts")
  )

  list(
    normalized_counts = DESeq2::counts(dds, normalized = TRUE),
    size_factors = DESeq2::sizeFactors(dds)
  )
}

compute_pooled_within_group_variance <- function(normalized_counts, group_labels) {
  groups <- levels(factor(group_labels))
  sse <- rep(0, nrow(normalized_counts))
  residual_df <- 0L

  for (g in groups) {
    idx <- which(group_labels == g)
    if (length(idx) < 2L) next
    xg <- normalized_counts[, idx, drop = FALSE]
    mu_g <- rowMeans(xg)
    resid_g <- sweep(xg, 1L, mu_g, "-")
    sse <- sse + rowSums(resid_g^2)
    residual_df <- residual_df + length(idx) - 1L
  }

  if (residual_df < 2L) stop("Pooled residual degrees of freedom < 2.")

  V <- sse / residual_df
  V[!is.finite(V)] <- 0
  V <- pmax(V, 0)
  names(V) <- rownames(normalized_counts)

  list(variance = V, residual_df = residual_df)
}

compute_pc1_rank <- function(rank_matrix_arm) {
  pca <- stats::prcomp(t(rank_matrix_arm), center = TRUE, scale. = FALSE, rank. = 1)
  loading <- pca$rotation[, 1L]
  loading[!is.finite(loading)] <- 0
  abs_loading <- abs(loading)
  rank_order <- order(abs_loading, decreasing = FALSE)
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

find_d2_zero_crossings <- function(x, d2) {
  ok <- is.finite(x) & is.finite(d2)
  x <- x[ok]
  d2 <- d2[ok]

  if (length(x) < 2L) return(numeric(0))

  out <- numeric(0)

  for (i in seq_len(length(x) - 1L)) {
    a <- d2[i]
    b <- d2[i + 1L]
    if (a == 0 || b == 0) next

    if ((a < 0 && b > 0) || (a > 0 && b < 0)) {
      frac <- abs(a) / (abs(a) + abs(b))
      out <- c(out, x[i] + frac * (x[i + 1L] - x[i]))
    }
  }

  sort(unique(out))
}

compute_raw_variance_geometry <- function(raw_counts_arm, rank_order) {
  raw_var <- apply(raw_counts_arm, 1L, stats::var, na.rm = TRUE)
  raw_var[!is.finite(raw_var)] <- 0
  raw_var <- pmax(raw_var, 0)

  raw_var_ranked <- raw_var[rank_order]
  rank <- seq_along(rank_order)
  log_var <- log1p(raw_var_ranked)

  analysis_spline <- stats::smooth.spline(x = rank, y = log_var, spar = ANALYSIS_SPLINE_SPAR)

  dense_x <- seq(min(rank), max(rank), length.out = max(5000L, length(rank) * 4L))
  dense_y <- as.numeric(stats::predict(analysis_spline, x = dense_x, deriv = 0)$y)
  dense_d1 <- as.numeric(stats::predict(analysis_spline, x = dense_x, deriv = 1)$y)
  dense_d2 <- as.numeric(stats::predict(analysis_spline, x = dense_x, deriv = 2)$y)

  dense_df <- data.frame(
    dense_rank = dense_x,
    dense_smooth_log1p_raw_variance = dense_y,
    dense_d1 = dense_d1,
    dense_d2 = dense_d2,
    stringsAsFactors = FALSE
  )

  crossings <- find_d2_zero_crossings(dense_df$dense_rank, dense_df$dense_d2)

  display_spline <- stats::smooth.spline(x = rank, y = log_var, spar = DISPLAY_VAR_SPAR)
  display_y <- as.numeric(stats::predict(display_spline, x = rank, deriv = 0)$y)

  curve <- data.frame(
    rank = rank,
    raw_empirical_variance = raw_var_ranked,
    log1p_raw_empirical_variance = log_var,
    display_log1p_raw_empirical_variance = display_y,
    stringsAsFactors = FALSE
  )

  list(curve = curve, dense = dense_df, crossings = crossings)
}

select_anchor_terminal_after_c2 <- function(crossings, c2, total_n) {
  z <- sort(unique(crossings[is.finite(crossings) & crossings > c2 & crossings < total_n]))

  if (length(z) < 2L) {
    stop(
      "Fewer than two second-derivative sign-change zero crossings were found after c2 = ",
      c2, "."
    )
  }

  anchor <- as.integer(round(z[1L]))
  terminal <- as.integer(round(z[2L]))

  anchor <- max(c2 + 1L, min(total_n - 1L, anchor))
  terminal <- max(anchor + 1L, min(total_n, terminal))

  if (anchor <= c2 || terminal <= anchor) {
    stop("Invalid Anchor-Terminal interval after c2.")
  }

  list(anchor = anchor, terminal = terminal)
}

smooth_nonnegative_mass <- function(rank, mass, spar = DISPLAY_MASS_SPAR) {
  fit <- stats::smooth.spline(x = rank, y = mass, spar = spar)
  y <- as.numeric(stats::predict(fit, x = rank, deriv = 0)$y)
  y[!is.finite(y)] <- 0
  y <- pmax(y, 0)
  total <- sum(y)
  if (!is.finite(total) || total <= 0) stop("Display mass could not be normalized.")
  y / total
}

smooth_divergence_for_display <- function(rank, D, spar = DISPLAY_DIV_SPAR) {
  fit <- stats::smooth.spline(x = rank, y = D, spar = spar)
  y <- as.numeric(stats::predict(fit, x = rank, deriv = 0)$y)
  endpoint_line <- seq(y[1L], y[length(y)], length.out = length(y))
  y - endpoint_line
}

compute_group_analysis <- function(group_name, raw_counts_arm, normalized_counts_arm, pooled_variance) {
  rank_matrix <- normalize_cpm_log1p(raw_counts_arm)
  pc1 <- compute_pc1_rank(rank_matrix)
  rank_order <- pc1$rank_order

  geometry <- compute_raw_variance_geometry(raw_counts_arm = raw_counts_arm, rank_order = rank_order)

  # Raw-count metrics for original-style corroboration.
  raw_ranked_mat <- raw_counts_arm[rank_order, , drop = FALSE]
  mu_raw <- rowMeans(raw_ranked_mat, na.rm = TRUE)
  var_raw <- apply(raw_ranked_mat, 1L, stats::var, na.rm = TRUE)
  mu_raw[!is.finite(mu_raw)] <- 0
  var_raw[!is.finite(var_raw)] <- 0
  mu_raw <- pmax(mu_raw, 0)
  var_raw <- pmax(var_raw, 0)

  nb2_var_minus_mu <- pmax(var_raw - mu_raw, 0)
  alpha_hat <- rep(0, length(mu_raw))
  pos <- mu_raw > 0
  alpha_hat[pos] <- pmax((var_raw[pos] - mu_raw[pos]) / (mu_raw[pos]^2), 0)

  # Normalized-count mean and pooled excess variance for divergence.
  mu_norm <- rowMeans(normalized_counts_arm, na.rm = TRUE)
  mu_norm[!is.finite(mu_norm)] <- 0
  mu_norm <- pmax(mu_norm, 0)

  P_ranked <- pc1$pc1_variance_contribution[rank_order]
  mu_norm_ranked <- mu_norm[rank_order]
  V_pool_ranked <- pooled_variance[rank_order]
  E_ranked <- pmax(V_pool_ranked - mu_norm_ranked, 0)

  P_total <- sum(P_ranked)
  E_total <- sum(E_ranked)

  if (!is.finite(P_total) || P_total <= 0) stop("PC1 variance mass undefined for group ", group_name)
  if (!is.finite(E_total) || E_total <= 0) stop("NB excess-variance mass undefined for group ", group_name)

  p_mass <- P_ranked / P_total
  q_mass <- E_ranked / E_total
  F_P <- cumsum(p_mass)
  F_E <- cumsum(q_mass)
  D <- F_E - F_P

  rank <- seq_along(rank_order)

  display_p <- smooth_nonnegative_mass(rank, p_mass)
  display_q <- smooth_nonnegative_mass(rank, q_mass)
  display_F_P <- cumsum(display_p)
  display_F_E <- cumsum(display_q)
  display_D <- smooth_divergence_for_display(rank, display_F_E - display_F_P)

  df <- geometry$curve %>%
    mutate(
      group = group_name,
      feature_id = rownames(raw_counts_arm)[rank_order],
      pc1_loading = pc1$loading[rank_order],
      abs_pc1_loading = pc1$abs_loading[rank_order],
      pc1_eigenvalue = pc1$lambda1,
      pc1_variance_contribution = P_ranked,
      pc1_variance_mass = p_mass,

      pooled_normalized_variance = V_pool_ranked,
      normalized_group_mean = mu_norm_ranked,
      nb_excess_variance = E_ranked,
      nb_excess_variance_mass = q_mass,

      cumulative_pc1_mass = F_P,
      cumulative_nb_mass = F_E,
      cumulative_divergence = D,
      display_F_P = display_F_P,
      display_F_E = display_F_E,
      display_D = display_D,

      raw_mean = mu_raw,
      raw_empirical_variance = var_raw,
      NB2 = log1p(nb2_var_minus_mu),
      NB2_NB1 = log1p(nb2_var_minus_mu) - log1p(mu_raw),
      alpha_mu = log1p(alpha_hat * mu_raw)
    )

  list(
    data = df,
    dense_raw_variance_geometry = geometry$dense,
    raw_variance_crossings = geometry$crossings,
    anchor = NA_integer_,
    terminal = NA_integer_
  )
}

piecewise_basis <- function(x, c1, c2) {
  cbind(
    intercept = 1,
    x = x,
    hinge1 = pmax(x - c1, 0),
    hinge2 = pmax(x - c2, 0)
  )
}

piecewise_sse <- function(par, x, D_mat, min_gap) {
  c1 <- par[1L]
  c2 <- par[2L]

  if (!is.finite(c1) || !is.finite(c2) || c1 <= 0 || c2 >= 1 || c2 - c1 <= min_gap) {
    return(1e100)
  }

  X <- piecewise_basis(x, c1, c2)
  coef <- tryCatch(qr.coef(qr(X), D_mat), error = function(e) NULL)

  if (is.null(coef) || any(!is.finite(coef))) return(1e100)

  resid <- D_mat - X %*% coef
  sum(resid^2)
}

fit_shared_knots <- function(group_results) {
  groups <- names(group_results)
  N_values <- vapply(group_results, function(z) nrow(z$data), integer(1))
  if (length(unique(N_values)) != 1L) stop("All groups must contain the same number of ranked features.")
  N <- N_values[1L]

  x_full <- (seq_len(N) - 1) / (N - 1)
  D_full <- do.call(cbind, lapply(group_results, function(z) z$data$cumulative_divergence))
  colnames(D_full) <- groups

  opt_n <- min(5000L, N)
  opt_idx <- unique(as.integer(round(seq(1, N, length.out = opt_n))))
  x_opt <- x_full[opt_idx]
  D_opt <- D_full[opt_idx, , drop = FALSE]

  min_gap <- max(4 / (N - 1), .Machine$double.eps^0.25)

  starts <- list(c(0.03, 0.97), c(0.08, 0.92), c(0.15, 0.85), c(0.25, 0.75), c(0.35, 0.65))

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
        control = list(maxit = 700, reltol = 1e-11)
      )
    }
  )

  values <- vapply(coarse, function(z) z$value, numeric(1))
  best <- coarse[[which.min(values)]]

  refined <- stats::optim(
    par = best$par,
    fn = piecewise_sse,
    x = x_full,
    D_mat = D_full,
    min_gap = min_gap,
    method = "Nelder-Mead",
    control = list(maxit = 1000, reltol = 1e-12)
  )

  if (!is.finite(refined$value)) stop("Shared-knot optimization failed.")

  c1_rank <- as.integer(round(1 + refined$par[1L] * (N - 1)))
  c2_rank <- as.integer(round(1 + refined$par[2L] * (N - 1)))

  c1_rank <- max(2L, min(N - 2L, c1_rank))
  c2_rank <- max(c1_rank + 1L, min(N - 1L, c2_rank))

  c1_x <- (c1_rank - 1) / (N - 1)
  c2_x <- (c2_rank - 1) / (N - 1)

  X <- piecewise_basis(x_full, c1_x, c2_x)
  coef <- qr.coef(qr(X), D_full)
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

estimate_nb_exponent <- function(mu, E) {
  keep <- is.finite(mu) & is.finite(E) & mu > 0 & E > 0
  if (sum(keep) < 10L) return(NA_real_)
  fit <- stats::lm(log(E[keep]) ~ log(mu[keep]))
  unname(stats::coef(fit)[2L])
}

get_region_p <- function(df, keep) {
  estimate_nb_exponent(mu = df$normalized_group_mean[keep], E = df$nb_excess_variance[keep])
}

safe_summary <- function(x, fun = median) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  fun(x)
}

get_left_right_indices <- function(total_n, anchor) {
  right_idx <- seq.int(anchor, total_n)
  right_n <- length(right_idx)
  left_end <- anchor - 1L
  left_start <- left_end - right_n + 1L
  if (left_start < 1L) stop("LEFT block extends below rank 1 for anchor = ", anchor)
  left_idx <- seq.int(left_start, left_end)
  list(left_idx = left_idx, right_idx = right_idx)
}

summarize_original_regions <- function(df, anchor) {
  idx <- get_left_right_indices(nrow(df), anchor)
  left_df <- df[idx$left_idx, , drop = FALSE]
  right_df <- df[idx$right_idx, , drop = FALSE]

  data.frame(
    left_n = nrow(left_df),
    right_n = nrow(right_df),
    left_NB2 = median(left_df$NB2, na.rm = TRUE),
    right_NB2 = median(right_df$NB2, na.rm = TRUE),
    diff_NB2 = median(right_df$NB2, na.rm = TRUE) - median(left_df$NB2, na.rm = TRUE),
    left_gap = median(left_df$NB2_NB1, na.rm = TRUE),
    right_gap = median(right_df$NB2_NB1, na.rm = TRUE),
    diff_gap = median(right_df$NB2_NB1, na.rm = TRUE) - median(left_df$NB2_NB1, na.rm = TRUE),
    left_alpha = median(left_df$alpha_mu, na.rm = TRUE),
    right_alpha = median(right_df$alpha_mu, na.rm = TRUE),
    diff_alpha = median(right_df$alpha_mu, na.rm = TRUE) - median(left_df$alpha_mu, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

build_timepoint_table <- function(group_results, comparisons, c1, c2, reference_rank) {
  rows <- list()

  for (comparison_name in names(comparisons)) {
    mapping <- comparisons[[comparison_name]]

    for (arm in names(mapping)) {
      g <- unname(mapping[[arm]])
      z <- group_results[[g]]
      df <- z$data
      anchor <- z$anchor
      terminal <- z$terminal
      region_summary <- summarize_original_regions(df, anchor)

      rows[[length(rows) + 1L]] <- data.frame(
        comparison = comparison_name,
        arm = arm,
        group = g,
        c1 = c1,
        c2 = c2,
        anchor = anchor,
        terminal = terminal,
        reference_rank = reference_rank,
        right_n = region_summary$right_n,
        p_remainder = get_region_p(df, df$rank < c1),
        p_model_leading_edge = get_region_p(df, df$rank > c2),
        p_paper_5000 = get_region_p(df, df$rank >= reference_rank),
        left_NB2 = region_summary$left_NB2,
        right_NB2 = region_summary$right_NB2,
        left_gap = region_summary$left_gap,
        right_gap = region_summary$right_gap,
        left_alpha = region_summary$left_alpha,
        right_alpha = region_summary$right_alpha,
        stringsAsFactors = FALSE
      )
    }
  }

  bind_rows(rows)
}

build_key_table <- function(timepoint_table, N, c1, c2, reference_rank, shared_sse) {
  data.frame(
    n_features = N,
    shared_c1 = c1,
    shared_c2 = c2,
    remainder_n = c1 - 1L,
    divergence_interval_n = c2 - c1 + 1L,
    model_leading_edge_n = N - c2,

    anchor_median = safe_summary(timepoint_table$anchor),
    anchor_min = safe_summary(timepoint_table$anchor, min),
    anchor_max = safe_summary(timepoint_table$anchor, max),

    terminal_median = safe_summary(timepoint_table$terminal),
    terminal_min = safe_summary(timepoint_table$terminal, min),
    terminal_max = safe_summary(timepoint_table$terminal, max),

    paper_ref = reference_rank,
    paper_leading_edge_n = PAPER_LEADING_EDGE_SIZE,
    paper_ref_used_in_analysis = FALSE,

    p_remainder_median = safe_summary(timepoint_table$p_remainder),
    p_model_leading_edge_median = safe_summary(timepoint_table$p_model_leading_edge),
    p_paper_5000_median = safe_summary(timepoint_table$p_paper_5000),

    left_NB2_median = safe_summary(timepoint_table$left_NB2),
    right_NB2_median = safe_summary(timepoint_table$right_NB2),
    left_gap_median = safe_summary(timepoint_table$left_gap),
    right_gap_median = safe_summary(timepoint_table$right_gap),
    left_alpha_median = safe_summary(timepoint_table$left_alpha),
    right_alpha_median = safe_summary(timepoint_table$right_alpha),

    shared_model_SSE = shared_sse,
    stringsAsFactors = FALSE
  )
}

rankwise_median <- function(group_results, column) {
  mat <- do.call(cbind, lapply(group_results, function(z) z$data[[column]]))
  apply(mat, 1L, median, na.rm = TRUE)
}

build_overall_figure_data <- function(group_results, knot_fit) {
  N <- nrow(group_results[[1L]]$data)

  data.frame(
    rank = seq_len(N),
    raw_variance = rankwise_median(group_results, "display_log1p_raw_empirical_variance"),
    F_P = rankwise_median(group_results, "display_F_P"),
    F_E = rankwise_median(group_results, "display_F_E"),
    D = rankwise_median(group_results, "display_D"),
    D_fit = apply(knot_fit$fitted, 1L, median, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# FULLY LABELED MANUSCRIPT FIGURES
# =============================================================================

# Override the earlier theme with a version that makes figure keys and
# method captions readable in exported manuscript figures.
theme_manuscript <- function(base_size = 11.3) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.2),
      plot.subtitle = element_text(size = base_size - 0.4, margin = margin(b = 5)),
      plot.caption = element_text(
        size = base_size - 1.45,
        hjust = 0,
        color = "#333333",
        margin = margin(t = 5)
      ),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#333333"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size - 1.5),
      legend.box = "vertical",
      panel.border = element_rect(color = "#B7B7B7", fill = NA, linewidth = 0.45),
      panel.grid = element_blank(),
      plot.margin = margin(7, 9, 7, 9)
    )
}

divergence_region_df <- function(c1, c2, N) {
  data.frame(
    xmin = c(1, c1, c2),
    xmax = c(c1, c2, N),
    region = factor(
      c("Remainder", "Divergence interval", "Leading-edge regime"),
      levels = c("Remainder", "Divergence interval", "Leading-edge regime")
    ),
    stringsAsFactors = FALSE
  )
}

add_divergence_region_shading <- function(p, c1, c2, N) {
  region_df <- divergence_region_df(c1, c2, N)

  p +
    geom_rect(
      data = region_df,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = region
      ),
      inherit.aes = FALSE,
      alpha = 0.48
    ) +
    scale_fill_manual(
      name = "PC1-rank region",
      values = c(
        "Remainder" = COL$remainder,
        "Divergence interval" = COL$interval,
        "Leading-edge regime" = COL$leading
      ),
      drop = FALSE
    ) +
    annotate(
      "text",
      x = (1 + c1) / 2,
      y = Inf,
      label = "REMAINDER",
      vjust = 1.35,
      fontface = "bold",
      size = 2.9,
      color = "#2F5D62"
    ) +
    annotate(
      "text",
      x = (c1 + c2) / 2,
      y = Inf,
      label = "DIVERGENCE INTERVAL",
      vjust = 1.35,
      fontface = "bold",
      size = 2.9,
      color = "#8A6500"
    ) +
    annotate(
      "text",
      x = (c2 + N) / 2,
      y = Inf,
      label = "LEADING-EDGE REGIME",
      vjust = 1.35,
      fontface = "bold",
      size = 2.9,
      color = "#166B5D"
    )
}

add_combined_line_scales <- function(
    p,
    color_values,
    linetype_values,
    breaks,
    title = "Curves and boundaries") {

  p +
    scale_color_manual(
      name = title,
      values = color_values,
      breaks = breaks,
      drop = FALSE
    ) +
    scale_linetype_manual(
      name = title,
      values = linetype_values,
      breaks = breaks,
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        nrow = 1
      ),
      color = guide_legend(
        order = 2,
        nrow = 2,
        byrow = TRUE,
        override.aes = list(linewidth = 1.1)
      ),
      linetype = guide_legend(
        order = 2,
        nrow = 2,
        byrow = TRUE,
        override.aes = list(linewidth = 1.1)
      )
    )
}


# =============================================================================
# OVERALL / CONSENSUS FIGURE
# =============================================================================

make_overall_divergence_figure <- function(
    overall_df,
    key_table,
    c1,
    c2,
    reference_rank,
    anchor_median,
    terminal_median,
    out_file) {

  N <- nrow(overall_df)

  boundary_df <- data.frame(
    rank = c(
      c1,
      c2,
      anchor_median,
      terminal_median,
      reference_rank
    ),
    key = c(
      "c1: shared remainder/divergence boundary",
      "c2: shared leading-edge start",
      "Median Anchor across 8 arms",
      "Median Terminal across 8 arms",
      "5,000-feature manuscript cutoff (REFERENCE ONLY)"
    ),
    stringsAsFactors = FALSE
  )

  boundary_colors <- c(
    "c1: shared remainder/divergence boundary" = COL$c1,
    "c2: shared leading-edge start" = COL$c2,
    "Median Anchor across 8 arms" = COL$anchor,
    "Median Terminal across 8 arms" = COL$terminal,
    "5,000-feature manuscript cutoff (REFERENCE ONLY)" = COL$ref
  )

  boundary_types <- c(
    "c1: shared remainder/divergence boundary" = "dashed",
    "c2: shared leading-edge start" = "longdash",
    "Median Anchor across 8 arms" = "solid",
    "Median Terminal across 8 arms" = "dotted",
    "5,000-feature manuscript cutoff (REFERENCE ONLY)" = "dotdash"
  )

  # -----------------------------------------------------------------------
  # Panel A: raw-count empirical variance
  # -----------------------------------------------------------------------

  a_curve <- "Median smoothed within-arm RAW-COUNT empirical variance"

  a_colors <- c(
    stats::setNames(COL$empirical, a_curve),
    boundary_colors
  )

  a_types <- c(
    stats::setNames("solid", a_curve),
    boundary_types
  )

  a_breaks <- c(
    a_curve,
    "c1: shared remainder/divergence boundary",
    "c2: shared leading-edge start",
    "Median Anchor across 8 arms",
    "Median Terminal across 8 arms",
    "5,000-feature manuscript cutoff (REFERENCE ONLY)"
  )

  pA <- ggplot()
  pA <- add_divergence_region_shading(pA, c1, c2, N)

  pA <- pA +
    geom_vline(
      data = boundary_df,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      inherit.aes = FALSE,
      linewidth = 0.86
    ) +
    geom_line(
      data = overall_df,
      aes(
        x = rank,
        y = raw_variance,
        color = a_curve,
        linetype = a_curve
      ),
      linewidth = 1.28,
      lineend = "round"
    ) +
    labs(
      title = "A. PC1-ranked raw-count empirical variance geometry",
      subtitle = paste0(
        "Shared c2 = ", c2,
        " starts the data-derived leading-edge regime; Anchor/Terminal are curvature markers within that regime."
      ),
      x = "PC1 rank: low |loading|  ->  high |loading|",
      y = "Smoothed log(1 + within-arm RAW-COUNT sample variance)",
      caption = paste0(
        "Variance source: within-arm empirical sample variance of RAW read counts. ",
        "The second derivative is taken only on this smoothed raw-count variance curve. ",
        "Starting immediately after c2, Anchor is the first y''(r)=0 sign-change crossing and Terminal is the second successive crossing. ",
        "The manuscript 5,000-feature cutoff begins at rank ", reference_rank,
        " and is shown as REFERENCE ONLY; it is not used to estimate c1, c2, Anchor, or Terminal."
      )
    ) +
    theme_manuscript()

  pA <- add_combined_line_scales(
    pA,
    color_values = a_colors,
    linetype_values = a_types,
    breaks = a_breaks,
    title = "Variance curve and boundaries"
  )

  # -----------------------------------------------------------------------
  # Panel B: cumulative masses and divergence
  # -----------------------------------------------------------------------

  cumulative_long <- overall_df %>%
    select(rank, F_P, F_E, D) %>%
    pivot_longer(
      cols = c(F_P, F_E, D),
      names_to = "curve",
      values_to = "value"
    ) %>%
    mutate(
      curve = factor(
        curve,
        levels = c("F_P", "F_E", "D"),
        labels = c(
          "F_P(r): cumulative PC1 variance mass",
          "F_E(r): cumulative NB excess-variance mass",
          "D(r)=F_E(r)-F_P(r): cumulative divergence"
        )
      )
    )

  b_curve_colors <- c(
    "F_P(r): cumulative PC1 variance mass" = COL$pc1,
    "F_E(r): cumulative NB excess-variance mass" = COL$nb,
    "D(r)=F_E(r)-F_P(r): cumulative divergence" = COL$divergence
  )

  b_curve_types <- c(
    "F_P(r): cumulative PC1 variance mass" = "solid",
    "F_E(r): cumulative NB excess-variance mass" = "solid",
    "D(r)=F_E(r)-F_P(r): cumulative divergence" = "solid"
  )

  b_colors <- c(
    b_curve_colors,
    boundary_colors
  )

  b_types <- c(
    b_curve_types,
    boundary_types
  )

  b_breaks <- c(
    names(b_curve_colors),
    "c1: shared remainder/divergence boundary",
    "c2: shared leading-edge start",
    "Median Anchor across 8 arms",
    "Median Terminal across 8 arms",
    "5,000-feature manuscript cutoff (REFERENCE ONLY)"
  )

  pB <- ggplot()
  pB <- add_divergence_region_shading(pB, c1, c2, N)

  pB <- pB +
    geom_vline(
      data = boundary_df,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      inherit.aes = FALSE,
      linewidth = 0.86
    ) +
    geom_line(
      data = cumulative_long,
      aes(
        x = rank,
        y = value,
        color = curve,
        linetype = curve
      ),
      linewidth = 1.07,
      lineend = "round"
    ) +
    geom_hline(
      yintercept = 0,
      color = "#555555",
      linetype = "dotted",
      linewidth = 0.40
    ) +
    labs(
      title = "B. Cumulative PC1-NB variance-mass divergence",
      subtitle = "Purpose: estimate shared c1/c2 by comparing where PC1 variance and normalized-count excess variance accumulate across rank.",
      x = "PC1 rank",
      y = "Cumulative variance mass / divergence",
      caption = paste0(
        "PC1 variance: P_i=lambda_1*v_i1^2, normalized as p(r)=P(r)/sum(P), with F_P(r)=sum_{j<=r}p(j). ",
        "NB excess variance: E=max(V_pool-mu_g,0), where V_pool is pooled within-group empirical variance of DESeq2 size-factor-normalized counts; ",
        "q(r)=E(r)/sum(E), F_E(r)=sum_{j<=r}q(j), and D(r)=F_E(r)-F_P(r). ",
        "The shared two-knot model fitted to D(r) across all 8 arms estimates c1 and c2. ",
        "The 5,000-feature reference does not enter this fit."
      )
    ) +
    theme_manuscript()

  pB <- add_combined_line_scales(
    pB,
    color_values = b_colors,
    linetype_values = b_types,
    breaks = b_breaks,
    title = "Cumulative curves and boundaries"
  )

  # -----------------------------------------------------------------------
  # Panel C: shared fit + explicit role of NB scaling
  # -----------------------------------------------------------------------

  c_obs <- "Observed median D(r) across 8 arms"
  c_fit <- "Shared two-knot fitted D(r)"

  c_colors <- c(
    stats::setNames(COL$nb, c_obs),
    stats::setNames(COL$fit, c_fit),
    boundary_colors
  )

  c_types <- c(
    stats::setNames("solid", c_obs),
    stats::setNames("solid", c_fit),
    boundary_types
  )

  c_breaks <- c(
    c_obs,
    c_fit,
    "c1: shared remainder/divergence boundary",
    "c2: shared leading-edge start",
    "Median Anchor across 8 arms",
    "Median Terminal across 8 arms",
    "5,000-feature manuscript cutoff (REFERENCE ONLY)"
  )

  p_summary <- paste0(
    "POST-BOUNDARY CORROBORATION ONLY\n",
    "E = alpha * mu^p (pooled normalized-count variance)\n",
    "Remainder median p = ", sprintf("%.2f", key_table$p_remainder_median), "\n",
    "Data-derived leading-edge median p = ", sprintf("%.2f", key_table$p_model_leading_edge_median), "\n",
    "Paper 5,000 median p = ", sprintf("%.2f", key_table$p_paper_5000_median), " (reference subset)\n",
    "p closer to 1 = more NB1-like; p closer to 2 = more NB2-like"
  )

  pC <- ggplot()
  pC <- add_divergence_region_shading(pC, c1, c2, N)

  pC <- pC +
    geom_vline(
      data = boundary_df,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      inherit.aes = FALSE,
      linewidth = 0.86
    ) +
    geom_line(
      data = overall_df,
      aes(
        x = rank,
        y = D,
        color = c_obs,
        linetype = c_obs
      ),
      linewidth = 0.88,
      alpha = 0.60,
      lineend = "round"
    ) +
    geom_line(
      data = overall_df,
      aes(
        x = rank,
        y = D_fit,
        color = c_fit,
        linetype = c_fit
      ),
      linewidth = 1.38,
      lineend = "round"
    ) +
    geom_hline(
      yintercept = 0,
      color = "#555555",
      linetype = "dotted",
      linewidth = 0.40
    ) +
    annotate(
      "label",
      x = round(0.59 * N),
      y = -Inf,
      label = p_summary,
      hjust = 0,
      vjust = -0.10,
      size = 2.55,
      fill = "white"
    ) +
    labs(
      title = "C. Shared two-knot divergence model and post-boundary NB scaling corroboration",
      subtitle = "c1/c2 come from cumulative divergence; Anchor/Terminal come from raw-count curvature after c2; NB scaling is descriptive corroboration afterward.",
      x = "PC1 rank",
      y = "Cumulative divergence D(r)",
      caption = paste0(
        "Shared model: D_g(x)=beta_0g+beta_1g*x+gamma_1g*(x-c1)_+ + gamma_2g*(x-c2)_+. ",
        "The NB exponent p is NOT used to choose any boundary. It characterizes mean-variance scaling only after regions are fixed. ",
        "Likewise, the paper 5,000-feature cutoff is a reference subset only."
      )
    ) +
    theme_manuscript()

  pC <- add_combined_line_scales(
    pC,
    color_values = c_colors,
    linetype_values = c_types,
    breaks = c_breaks,
    title = "Model curves and boundaries"
  )

  save_panels(
    list(pA, pB, pC),
    out_file,
    height_in = 13.3
  )
}


# =============================================================================
# NEW TIME-POINT FIGURES
# =============================================================================

make_new_timepoint_figure <- function(
    comparison_name,
    mapping,
    group_results,
    timepoint_table,
    c1,
    c2,
    reference_rank,
    out_file) {

  control_group <- unname(mapping[["control"]])
  treatment_group <- unname(mapping[["treatment"]])

  control_label <- paste0("Control (", control_group, ")")
  treatment_label <- paste0("Treatment (", treatment_group, ")")

  control <- group_results[[control_group]]$data
  treatment <- group_results[[treatment_group]]$data
  N <- nrow(control)

  term_rows <- timepoint_table %>%
    filter(comparison == comparison_name)

  c_anchor <- term_rows$anchor[term_rows$group == control_group]
  c_terminal <- term_rows$terminal[term_rows$group == control_group]
  t_anchor <- term_rows$anchor[term_rows$group == treatment_group]
  t_terminal <- term_rows$terminal[term_rows$group == treatment_group]

  c_anchor_key <- paste0(control_group, " Anchor: first y''=0 after c2")
  c_terminal_key <- paste0(control_group, " Terminal: second successive y''=0")
  t_anchor_key <- paste0(treatment_group, " Anchor: first y''=0 after c2")
  t_terminal_key <- paste0(treatment_group, " Terminal: second successive y''=0")

  boundary_df <- data.frame(
    rank = c(
      c1,
      c2,
      c_anchor,
      c_terminal,
      t_anchor,
      t_terminal,
      reference_rank
    ),
    key = c(
      "c1: shared remainder/divergence boundary",
      "c2: shared leading-edge start",
      c_anchor_key,
      c_terminal_key,
      t_anchor_key,
      t_terminal_key,
      "5,000-feature manuscript cutoff (REFERENCE ONLY)"
    ),
    stringsAsFactors = FALSE
  )

  boundary_colors <- c(
    "c1: shared remainder/divergence boundary" = COL$c1,
    "c2: shared leading-edge start" = COL$c2,
    stats::setNames(COL$control, c_anchor_key),
    stats::setNames(COL$control, c_terminal_key),
    stats::setNames(COL$treatment, t_anchor_key),
    stats::setNames(COL$treatment, t_terminal_key),
    "5,000-feature manuscript cutoff (REFERENCE ONLY)" = COL$ref
  )

  boundary_types <- c(
    "c1: shared remainder/divergence boundary" = "dashed",
    "c2: shared leading-edge start" = "longdash",
    stats::setNames("solid", c_anchor_key),
    stats::setNames("dotted", c_terminal_key),
    stats::setNames("solid", t_anchor_key),
    stats::setNames("dotted", t_terminal_key),
    "5,000-feature manuscript cutoff (REFERENCE ONLY)" = "dotdash"
  )

  # -----------------------------------------------------------------------
  # Panel A
  # -----------------------------------------------------------------------

  control_curve <- paste0(control_label, ": smoothed RAW-count variance")
  treatment_curve <- paste0(treatment_label, ": smoothed RAW-count variance")

  a_colors <- c(
    stats::setNames(COL$control, control_curve),
    stats::setNames(COL$treatment, treatment_curve),
    boundary_colors
  )

  a_types <- c(
    stats::setNames("solid", control_curve),
    stats::setNames("solid", treatment_curve),
    boundary_types
  )

  a_breaks <- c(
    control_curve,
    treatment_curve,
    "c1: shared remainder/divergence boundary",
    "c2: shared leading-edge start",
    c_anchor_key,
    c_terminal_key,
    t_anchor_key,
    t_terminal_key,
    "5,000-feature manuscript cutoff (REFERENCE ONLY)"
  )

  pA <- ggplot()
  pA <- add_divergence_region_shading(pA, c1, c2, N)

  pA <- pA +
    geom_vline(
      data = boundary_df,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      inherit.aes = FALSE,
      linewidth = 0.82
    ) +
    geom_line(
      data = control,
      aes(
        x = rank,
        y = display_log1p_raw_empirical_variance,
        color = control_curve,
        linetype = control_curve
      ),
      linewidth = 1.18,
      lineend = "round"
    ) +
    geom_line(
      data = treatment,
      aes(
        x = rank,
        y = display_log1p_raw_empirical_variance,
        color = treatment_curve,
        linetype = treatment_curve
      ),
      linewidth = 1.18,
      lineend = "round"
    ) +
    labs(
      title = paste0("A. ", comparison_name, ": raw-count variance geometry"),
      subtitle = paste0(
        "Shared c2 = ", c2,
        "; ", control_group, " Anchor/Terminal = ", c_anchor, "/", c_terminal,
        "; ", treatment_group, " Anchor/Terminal = ", t_anchor, "/", t_terminal, "."
      ),
      x = "PC1 rank: low |loading|  ->  high |loading|",
      y = "Smoothed log(1 + within-arm RAW-COUNT sample variance)",
      caption = paste0(
        "Method: y(r)=smooth{log[1+s_raw^2(r)]}. The second derivative y''(r) is calculated from this RAW-count variance spline only. ",
        "After shared c2, Anchor is the first sign-change zero crossing of y''(r); Terminal is the second successive crossing. ",
        "The 5,000-feature manuscript cutoff at rank ", reference_rank, " is REFERENCE ONLY."
      )
    ) +
    theme_manuscript()

  pA <- add_combined_line_scales(
    pA,
    color_values = a_colors,
    linetype_values = a_types,
    breaks = a_breaks,
    title = "Arm curves and boundaries"
  )

  # -----------------------------------------------------------------------
  # Panel B
  # -----------------------------------------------------------------------

  control_div <- paste0(control_label, ": D(r)=F_E-F_P")
  treatment_div <- paste0(treatment_label, ": D(r)=F_E-F_P")

  b_colors <- c(
    stats::setNames(COL$control, control_div),
    stats::setNames(COL$treatment, treatment_div),
    boundary_colors
  )

  b_types <- c(
    stats::setNames("solid", control_div),
    stats::setNames("solid", treatment_div),
    boundary_types
  )

  b_breaks <- c(
    control_div,
    treatment_div,
    "c1: shared remainder/divergence boundary",
    "c2: shared leading-edge start",
    c_anchor_key,
    c_terminal_key,
    t_anchor_key,
    t_terminal_key,
    "5,000-feature manuscript cutoff (REFERENCE ONLY)"
  )

  pB <- ggplot()
  pB <- add_divergence_region_shading(pB, c1, c2, N)

  pB <- pB +
    geom_vline(
      data = boundary_df,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      inherit.aes = FALSE,
      linewidth = 0.82
    ) +
    geom_hline(
      yintercept = 0,
      color = "#555555",
      linetype = "dotted",
      linewidth = 0.40
    ) +
    geom_line(
      data = control,
      aes(
        x = rank,
        y = display_D,
        color = control_div,
        linetype = control_div
      ),
      linewidth = 1.18,
      lineend = "round"
    ) +
    geom_line(
      data = treatment,
      aes(
        x = rank,
        y = display_D,
        color = treatment_div,
        linetype = treatment_div
      ),
      linewidth = 1.18,
      lineend = "round"
    ) +
    labs(
      title = paste0("B. ", comparison_name, ": cumulative PC1-NB variance-mass divergence"),
      subtitle = "Purpose: visualize the variance-structure divergence used to estimate the shared c1/c2 regime boundaries.",
      x = "PC1 rank",
      y = "D(r) = F_E(r) - F_P(r)",
      caption = paste0(
        "F_P(r) is cumulative PC1 variance mass. F_E(r) is cumulative excess-variance mass derived from pooled within-group empirical variance ",
        "of DESeq2 size-factor-normalized counts. Shared c1/c2 are estimated jointly across all 8 arms from D(r). ",
        "Arm-specific Anchor/Terminal lines are overlaid afterward from the separate raw-count second-derivative analysis. ",
        "The 5,000-feature cutoff remains reference only."
      )
    ) +
    theme_manuscript()

  pB <- add_combined_line_scales(
    pB,
    color_values = b_colors,
    linetype_values = b_types,
    breaks = b_breaks,
    title = "Divergence curves and boundaries"
  )

  # -----------------------------------------------------------------------
  # Panel C: explicitly described as corroboration, not boundary selection
  # -----------------------------------------------------------------------

  p_rows <- term_rows %>%
    mutate(
      arm_label = ifelse(
        arm == "control",
        control_label,
        treatment_label
      )
    ) %>%
    select(
      arm_label,
      p_remainder,
      p_model_leading_edge,
      p_paper_5000
    ) %>%
    pivot_longer(
      cols = c(
        p_remainder,
        p_model_leading_edge,
        p_paper_5000
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
          "p_paper_5000"
        ),
        labels = c(
          "Remainder (rank < c1)",
          "Data-derived leading edge (rank > c2)",
          "Paper terminal 5,000 (REFERENCE ONLY)"
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
      color = "#555555",
      linetype = "dashed",
      linewidth = 0.68
    ) +
    geom_vline(
      xintercept = 2,
      color = "#555555",
      linetype = "dotted",
      linewidth = 0.78
    ) +
    geom_point(size = 3.5) +
    scale_color_manual(
      name = "Experimental arm",
      values = c(
        stats::setNames(COL$control, control_label),
        stats::setNames(COL$treatment, treatment_label)
      )
    ) +
    labs(
      title = paste0("C. ", comparison_name, ": post-boundary NB mean-variance scaling"),
      subtitle = "CORROBORATION ONLY — this panel does NOT determine c1, c2, Anchor, or Terminal.",
      x = "Empirical exponent p in E = alpha * mu^p",
      y = "Region being characterized",
      caption = paste0(
        "Variance source for E: pooled within-group empirical variance of DESeq2-normalized counts. ",
        "Interpretation: p closer to 1 is more NB1-like; p closer to 2 is more NB2-like. ",
        "This panel asks whether the independently defined regions differ in mean-variance scaling after all boundaries are fixed."
      )
    ) +
    theme_manuscript()

  save_panels(
    list(pA, pB, pC),
    out_file,
    height_in = 13.3
  )
}


# =============================================================================
# ORIGINAL-STYLE FIGURES
# =============================================================================

add_original_region_shading <- function(
    p,
    left_min,
    left_max,
    anchor,
    terminal,
    N) {

  region_df <- data.frame(
    xmin = c(left_min, anchor, anchor),
    xmax = c(left_max, N, terminal),
    region = factor(
      c(
        "Matched LEFT comparator",
        "RIGHT: Anchor to end",
        "Anchor-Terminal curvature interval"
      ),
      levels = c(
        "Matched LEFT comparator",
        "RIGHT: Anchor to end",
        "Anchor-Terminal curvature interval"
      )
    ),
    stringsAsFactors = FALSE
  )

  p +
    geom_rect(
      data = region_df,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = region
      ),
      inherit.aes = FALSE,
      alpha = 0.50
    ) +
    scale_fill_manual(
      name = "Original LEFT/RIGHT geometry",
      values = c(
        "Matched LEFT comparator" = COL$left_fill,
        "RIGHT: Anchor to end" = COL$right_fill,
        "Anchor-Terminal curvature interval" = COL$interval_fill
      ),
      drop = FALSE
    ) +
    annotate(
      "text",
      x = (left_min + left_max) / 2,
      y = Inf,
      label = "MATCHED LEFT",
      vjust = 1.35,
      fontface = "bold",
      size = 2.9,
      color = "#315A7D"
    ) +
    annotate(
      "text",
      x = (anchor + N) / 2,
      y = Inf,
      label = "RIGHT = ANCHOR TO END",
      vjust = 1.35,
      fontface = "bold",
      size = 2.9,
      color = "#3B6E36"
    )
}

make_original_arm_panels <- function(
    df,
    arm_label,
    c2,
    reference_rank,
    anchor,
    terminal) {

  N <- nrow(df)
  idx <- get_left_right_indices(N, anchor)

  left_idx <- idx$left_idx
  right_idx <- idx$right_idx
  region_summary <- summarize_original_regions(df, anchor)

  left_min <- min(left_idx)
  left_max <- max(left_idx)
  right_n <- length(right_idx)

  boundary_df <- data.frame(
    rank = c(
      c2,
      anchor,
      terminal,
      reference_rank
    ),
    key = c(
      "c2: data-derived leading-edge start",
      "Anchor: first raw-variance y''=0 crossing after c2",
      "Terminal: second successive raw-variance y''=0 crossing",
      "5,000-feature manuscript cutoff (REFERENCE ONLY)"
    ),
    stringsAsFactors = FALSE
  )

  boundary_colors <- c(
    "c2: data-derived leading-edge start" = COL$c2,
    "Anchor: first raw-variance y''=0 crossing after c2" = COL$anchor,
    "Terminal: second successive raw-variance y''=0 crossing" = COL$terminal,
    "5,000-feature manuscript cutoff (REFERENCE ONLY)" = COL$ref
  )

  boundary_types <- c(
    "c2: data-derived leading-edge start" = "longdash",
    "Anchor: first raw-variance y''=0 crossing after c2" = "solid",
    "Terminal: second successive raw-variance y''=0 crossing" = "dotted",
    "5,000-feature manuscript cutoff (REFERENCE ONLY)" = "dotdash"
  )

  # -----------------------------------------------------------------------
  # Original panel 1: raw-count geometry
  # -----------------------------------------------------------------------

  raw_curve <- "Smoothed within-arm RAW-COUNT empirical variance"

  p1_colors <- c(
    stats::setNames(COL$empirical, raw_curve),
    boundary_colors
  )

  p1_types <- c(
    stats::setNames("solid", raw_curve),
    boundary_types
  )

  p1_breaks <- c(
    raw_curve,
    "c2: data-derived leading-edge start",
    "Anchor: first raw-variance y''=0 crossing after c2",
    "Terminal: second successive raw-variance y''=0 crossing",
    "5,000-feature manuscript cutoff (REFERENCE ONLY)"
  )

  p1 <- ggplot()
  p1 <- add_original_region_shading(
    p1,
    left_min = left_min,
    left_max = left_max,
    anchor = anchor,
    terminal = terminal,
    N = N
  )

  p1 <- p1 +
    geom_vline(
      data = boundary_df,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      inherit.aes = FALSE,
      linewidth = 0.86
    ) +
    geom_line(
      data = df,
      aes(
        x = rank,
        y = display_log1p_raw_empirical_variance,
        color = raw_curve,
        linetype = raw_curve
      ),
      linewidth = 1.06,
      lineend = "round"
    ) +
    labs(
      title = paste0(arm_label, ": original-style raw-count variance geometry"),
      subtitle = paste0(
        "c2 = ", c2,
        "; Anchor = ", anchor,
        "; Terminal = ", terminal,
        "; RIGHT n = ", right_n,
        "; manuscript 5,000 reference = ", reference_rank, "."
      ),
      x = "PC1 rank: low |loading|  ->  high |loading|",
      y = "Smoothed log(1 + within-arm RAW-COUNT sample variance)",
      caption = paste0(
        "The same second-derivative curve provides both markers: after c2, Anchor is the first y''(r)=0 sign-change crossing; ",
        "Terminal is the second successive crossing. RIGHT=[Anchor,N]. LEFT is the immediately preceding equal-sized rank block. ",
        "The 5,000-feature cutoff is REFERENCE ONLY and does not select Anchor, Terminal, LEFT, or RIGHT."
      )
    ) +
    theme_manuscript()

  p1 <- add_combined_line_scales(
    p1,
    color_values = p1_colors,
    linetype_values = p1_types,
    breaks = p1_breaks,
    title = "Variance curve and boundaries"
  )

  # -----------------------------------------------------------------------
  # Original panel 2: legacy NB-related corroboration signals
  # -----------------------------------------------------------------------

  nb_long <- df %>%
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
          "NB2 excess signal",
          "NB2-NB1 contrast",
          "alpha*mu signal"
        )
      )
    )

  metric_colors <- c(
    "NB2 excess signal" = COL$nb2,
    "NB2-NB1 contrast" = COL$nbgap,
    "alpha*mu signal" = COL$alphamu
  )

  metric_types <- c(
    "NB2 excess signal" = "solid",
    "NB2-NB1 contrast" = "solid",
    "alpha*mu signal" = "solid"
  )

  p2_colors <- c(
    metric_colors,
    boundary_colors
  )

  p2_types <- c(
    metric_types,
    boundary_types
  )

  p2_breaks <- c(
    "NB2 excess signal",
    "NB2-NB1 contrast",
    "alpha*mu signal",
    "c2: data-derived leading-edge start",
    "Anchor: first raw-variance y''=0 crossing after c2",
    "Terminal: second successive raw-variance y''=0 crossing",
    "5,000-feature manuscript cutoff (REFERENCE ONLY)"
  )

  p2 <- ggplot()
  p2 <- add_original_region_shading(
    p2,
    left_min = left_min,
    left_max = left_max,
    anchor = anchor,
    terminal = terminal,
    N = N
  )

  p2 <- p2 +
    geom_vline(
      data = boundary_df,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      inherit.aes = FALSE,
      linewidth = 0.86
    ) +
    geom_line(
      data = nb_long,
      aes(
        x = rank,
        y = value,
        color = metric,
        linetype = metric
      ),
      linewidth = 0.92
    ) +
    labs(
      title = paste0(arm_label, ": original LEFT/RIGHT NB-related corroboration"),
      subtitle = "Purpose: describe how the independently selected RIGHT region differs from its matched LEFT comparator after boundaries are fixed.",
      x = "PC1 rank",
      y = "Descriptive RAW-COUNT NB-related signal",
      caption = paste0(
        "NB2 excess signal = log[1+max(s_raw^2-mu_raw,0)]. ",
        "NB2-NB1 contrast = NB2 excess signal - log(1+mu_raw). ",
        "alpha_hat=max[(s_raw^2-mu_raw)/mu_raw^2,0], and alpha*mu signal=log(1+alpha_hat*mu_raw). ",
        "These are descriptive moment-based corroboration metrics, NOT formal NB1-versus-NB2 likelihood tests and NOT boundary-selection statistics."
      )
    ) +
    theme_manuscript()

  p2 <- add_combined_line_scales(
    p2,
    color_values = p2_colors,
    linetype_values = p2_types,
    breaks = p2_breaks,
    title = "NB-related signals and boundaries"
  )

  # -----------------------------------------------------------------------
  # Original panel 3: matched LEFT versus RIGHT summary
  # -----------------------------------------------------------------------

  summary_df <- data.frame(
    metric = factor(
      c(
        "NB2 excess signal",
        "NB2-NB1 contrast",
        "alpha*mu signal"
      ),
      levels = rev(
        c(
          "NB2 excess signal",
          "NB2-NB1 contrast",
          "alpha*mu signal"
        )
      )
    ),
    LEFT = c(
      region_summary$left_NB2,
      region_summary$left_gap,
      region_summary$left_alpha
    ),
    RIGHT = c(
      region_summary$right_NB2,
      region_summary$right_gap,
      region_summary$right_alpha
    ),
    stringsAsFactors = FALSE
  )

  p3 <- ggplot(
    summary_df,
    aes(y = metric)
  ) +
    geom_segment(
      aes(
        x = LEFT,
        xend = RIGHT,
        yend = metric
      ),
      color = "#7A7A7A",
      linewidth = 0.82
    ) +
    geom_point(
      aes(
        x = LEFT,
        color = "Matched LEFT"
      ),
      size = 3.5
    ) +
    geom_point(
      aes(
        x = RIGHT,
        color = "RIGHT: Anchor to end"
      ),
      size = 3.5
    ) +
    scale_color_manual(
      name = "Comparison block",
      values = c(
        "Matched LEFT" = "#5B8FD1",
        "RIGHT: Anchor to end" = "#43A047"
      )
    ) +
    labs(
      title = paste0(arm_label, ": matched LEFT versus RIGHT summary"),
      subtitle = paste0(
        "RIGHT begins at Anchor = ", anchor,
        "; LEFT contains the same number of immediately preceding ranks (n = ", right_n, ")."
      ),
      x = "Median descriptive signal",
      y = "Metric",
      caption = "This compact panel summarizes the original LEFT-versus-RIGHT corroboration after Anchor is selected independently of the 5,000-feature reference."
    ) +
    theme_manuscript()

  list(
    p1,
    p2,
    p3
  )
}

make_original_comparison_figure <- function(
    comparison_name,
    mapping,
    group_results,
    c2,
    reference_rank,
    out_file) {

  control_group <- unname(mapping[["control"]])
  treatment_group <- unname(mapping[["treatment"]])

  control <- group_results[[control_group]]
  treatment <- group_results[[treatment_group]]

  plots <- c(
    make_original_arm_panels(
      df = control$data,
      arm_label = paste0(
        comparison_name,
        " Control (",
        control_group,
        ")"
      ),
      c2 = c2,
      reference_rank = reference_rank,
      anchor = control$anchor,
      terminal = control$terminal
    ),
    make_original_arm_panels(
      df = treatment$data,
      arm_label = paste0(
        comparison_name,
        " Treatment (",
        treatment_group,
        ")"
      ),
      c2 = c2,
      reference_rank = reference_rank,
      anchor = treatment$anchor,
      terminal = treatment$terminal
    )
  )

  save_panels(
    plots,
    out_file,
    height_in = 24.0
  )
}


# =============================================================================
# RUN
# =============================================================================

count_mat <- read_count_matrix(COUNT_FILE, GROUP_PATTERNS)
group_labels <- assign_groups(colnames(count_mat), GROUP_PATTERNS)

message("Count matrix: ", nrow(count_mat), " features x ", ncol(count_mat), " samples")
message("Groups: ", paste(levels(group_labels), collapse = ", "))

N <- nrow(count_mat)
if (PAPER_LEADING_EDGE_SIZE >= N) stop("PAPER_LEADING_EDGE_SIZE must be smaller than the feature count.")
REFERENCE_RANK <- N - PAPER_LEADING_EDGE_SIZE + 1L

deseq <- normalize_deseq2_global(count_mat = count_mat, group_labels = group_labels)
normalized_counts <- deseq$normalized_counts

pooled <- compute_pooled_within_group_variance(normalized_counts = normalized_counts, group_labels = group_labels)
message("Pooled within-group residual degrees of freedom: ", pooled$residual_df)

group_results <- vector("list", length(levels(group_labels)))
names(group_results) <- levels(group_labels)

for (g in levels(group_labels)) {
  idx <- which(group_labels == g)
  if (length(idx) < 2L) stop("Not enough samples in group ", g)

  message("Analyzing ", g, "...")
  group_results[[g]] <- compute_group_analysis(
    group_name = g,
    raw_counts_arm = count_mat[, idx, drop = FALSE],
    normalized_counts_arm = normalized_counts[, idx, drop = FALSE],
    pooled_variance = pooled$variance
  )
}

message("Fitting shared two-knot cumulative-divergence model...")
knot_fit <- fit_shared_knots(group_results)
C1 <- knot_fit$c1
C2 <- knot_fit$c2

message("Selecting Anchor and Terminal from the first two successive y''(r)=0 crossings after c2...")
for (g in names(group_results)) {
  at <- select_anchor_terminal_after_c2(
    crossings = group_results[[g]]$raw_variance_crossings,
    c2 = C2,
    total_n = N
  )
  group_results[[g]]$anchor <- at$anchor
  group_results[[g]]$terminal <- at$terminal
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

write.csv(key_table, file.path(OUT_ROOT, "Table_Key_Results.csv"), row.names = FALSE)
write.csv(timepoint_table, file.path(OUT_ROOT, "Table_Timepoints.csv"), row.names = FALSE)

overall_df <- build_overall_figure_data(group_results = group_results, knot_fit = knot_fit)

figure_paths <- character(0)

overall_path <- file.path(FIG_DIR, "Figure_Overall_Divergence.png")
make_overall_divergence_figure(
  overall_df = overall_df,
  key_table = key_table,
  c1 = C1,
  c2 = C2,
  reference_rank = REFERENCE_RANK,
  anchor_median = key_table$anchor_median[1L],
  terminal_median = key_table$terminal_median[1L],
  out_file = overall_path
)
figure_paths <- c(figure_paths, overall_path)

for (comparison_name in names(COMPARISONS)) {
  new_path <- file.path(FIG_DIR, paste0("Figure_New_", comparison_name, ".png"))
  make_new_timepoint_figure(
    comparison_name = comparison_name,
    mapping = COMPARISONS[[comparison_name]],
    group_results = group_results,
    timepoint_table = timepoint_table,
    c1 = C1,
    c2 = C2,
    reference_rank = REFERENCE_RANK,
    out_file = new_path
  )
  figure_paths <- c(figure_paths, new_path)

  old_path <- file.path(FIG_DIR, paste0("Figure_Original_", comparison_name, ".png"))
  make_original_comparison_figure(
    comparison_name = comparison_name,
    mapping = COMPARISONS[[comparison_name]],
    group_results = group_results,
    c2 = C2,
    reference_rank = REFERENCE_RANK,
    out_file = old_path
  )
  figure_paths <- c(figure_paths, old_path)
}

ZIP_PATH <- file.path(OUT_ROOT, "Figures_All.zip")
if (file.exists(ZIP_PATH)) unlink(ZIP_PATH)

old_wd <- getwd()
zip_ok <- FALSE

tryCatch(
  {
    setwd(FIG_DIR)
    utils::zip(zipfile = ZIP_PATH, files = basename(figure_paths))
    zip_ok <- file.exists(ZIP_PATH)
  },
  finally = {
    setwd(old_wd)
  }
)

if (!zip_ok) {
  warning("Figure PNGs were created, but Figures_All.zip was not created.")
}

message("============================================================")
message("FINAL ANALYSIS COMPLETE")
message("Shared c1 = ", C1)
message("Shared c2 = ", C2)
message("Historical 5,000-feature Ref (REFERENCE ONLY) = ", REFERENCE_RANK)
message(
  "Arm-specific Anchor-Terminal ranges: ",
  paste0(timepoint_table$group, "=", timepoint_table$anchor, "-", timepoint_table$terminal, collapse = "; ")
)
message("Figures directory: ", FIG_DIR)
message("Figure zip: ", ZIP_PATH)
message("Tables: Table_Key_Results.csv; Table_Timepoints.csv")
message("============================================================")
