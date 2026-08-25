#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# EMPIRICAL-k* LEFT/RIGHT NB1/NB2 CORROBORATION
# =============================================================================
#
# PURPOSE
# -------
# Use the already-established comparison-specific empirical EVS cutoffs:
#
#   RT0_ZT6  k* = 3532
#   RT2_ZT8  k* = 4617
#   RT4_ZT10 k* = 3983
#   RT8_ZT14 k* = 5664
#
# RANKING AND REFERENCE
# ---------------------
# For each arm, PCA is performed on arm-specific log1p(CPM) expression with
# centering and without feature scaling. PASs are ordered from LOW to HIGH
# absolute PC1 loading, so the comparison-specific top-k* begins at:
#
#   reference_rank = N - k* + 1.
#
# TRANSITION RANGE
# ----------------
# The original REMNB1/LEADNB2 geometry is preserved. Feature-wise empirical
# variance is calculated from RAW COUNTS across the biological replicates in
# the arm, transformed as log1p[Var(raw count)], and smoothed along the PC1
# rank axis. The nearest second-derivative sign changes immediately to the
# left and right of the k* reference define:
#
#   Anchor < k* reference < Terminal.
#
# RIGHT is Anchor through the end of the rank axis. LEFT is the immediately
# preceding block with the same number of PASs as RIGHT. Thus the displayed
# LEFT/RIGHT range is built from the comparison-specific empirical k*.
#
# DESCRIPTIVE RAW-COUNT NB2 QUANTITIES
# ------------------------------------
# For each PAS, with raw-count sample mean mu and raw-count sample variance s^2:
#
#   NB2 excess = log1p[max(s^2 - mu, 0)]
#
#   NB2-NB1    = log1p[max(s^2 - mu, 0)] - log1p(mu)
#
#   alpha_hat  = max[(s^2 - mu)/mu^2, 0]
#   alpha*mu   = log1p(alpha_hat * mu)
#
# RIGHT-LEFT values shown in the figures are DESCRIPTIVE differences between
# the corresponding region medians.
#
# FORMAL LEFT-vs-RIGHT NB2 DISPERSION LRT
# ---------------------------------------
# The formal p-value is calculated from RAW INTEGER COUNTS only.
#
#   H0: LEFT and RIGHT share one NB2 dispersion parameter alpha.
#   H1: LEFT and RIGHT have separate NB2 dispersion parameters.
#
# Feature-specific empirical means are held fixed in both models. Under the
# NB2 parameterization Var(Y)=mu+alpha*mu^2, dnbinom size=1/alpha.
#
#   LRT = 2 * (logLik_H1 - logLik_H0)
#
# H1 has one additional dispersion parameter, so the reference distribution is
# chi-square with 1 df. The figure labels this specifically as the
# "LEFT-vs-RIGHT NB2-dispersion LRT p"; it is not presented as a p-value for
# the descriptive median NB2-NB1 difference.
#
# FIGURES
# -------
# One three-panel figure is produced for each of the eight arms:
#
#   A. Smoothed log1p raw-count variance and the k*-anchored transition range.
#   B. Raw-count NB2 corroboration signals along the PC1 rank.
#   C. LEFT vs RIGHT median NB2-NB1 contrast, RIGHT-LEFT difference, and the
#      formal LEFT-vs-RIGHT NB2-dispersion LRT p-value.
#
# =============================================================================

# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT   <- "/root/REAPER98632/exports/remnb1_leadnb2_empirical_kstar_simple"

# Locked inputs from the empirical weighted-Pareto cutoff analysis.
EMPIRICAL_K <- c(
  RT0_ZT6  = 3532L,
  RT2_ZT8  = 4617L,
  RT4_ZT10 = 3983L,
  RT8_ZT14 = 5664L
)

VAR_SPLINE_SPAR <- 0.60
NB_DISPLAY_SPAR  <- 0.65

PNG_WIDTH_IN  <- 14
PNG_HEIGHT_IN <- 10.8
PNG_DPI       <- 300

COMPARISONS <- list(
  RT0_ZT6 = list(
    RT0 = "^R0_",
    ZT6 = "^ZT6_"
  ),
  RT2_ZT8 = list(
    RT2 = "^R2_",
    ZT8 = "^ZT8_"
  ),
  RT4_ZT10 = list(
    RT4 = "^R4_",
    ZT10 = "^ZT10_"
  ),
  RT8_ZT14 = list(
    RT8 = "^R8_",
    ZT14 = "^ZT14_"
  )
)

if (!identical(names(EMPIRICAL_K), names(COMPARISONS))) {
  stop("EMPIRICAL_K and COMPARISONS must contain the same comparisons in the same order.")
}

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# COLORS
# =============================================================================

COL <- list(
  var_curve = "#117A65",
  nb2       = "#1B9E77",
  nb_gap    = "#CC1E8C",
  alpha_mu  = "#386CB0",

  left_fill     = "#CBE3F8",
  right_fill    = "#DDF2D5",
  interval_fill = "#BDBDBD",

  anchor   = "#000000",
  kstar    = "#E69F00",
  terminal = "#D95F02",

  left_pt  = "#5B8FD1",
  right_pt = "#43A047"
)

TRACE_LEVELS <- c("NB2 excess", "NB2-NB1", "alpha*mu")
TRACE_COLORS <- c(
  "NB2 excess" = COL$nb2,
  "NB2-NB1" = COL$nb_gap,
  "alpha*mu" = COL$alpha_mu
)

# =============================================================================
# HELPERS
# =============================================================================

first_numeric_col_index <- function(df) {
  idx <- which(vapply(df, is.numeric, logical(1)))
  if (length(idx) == 0L) return(NA_integer_)
  idx[1]
}

read_count_matrix <- function(path) {
  raw_df <- read.csv(path, check.names = FALSE)

  if (nrow(raw_df) == 0L || ncol(raw_df) < 2L) {
    stop("Count file is empty or malformed: ", path)
  }

  first_num <- first_numeric_col_index(raw_df)
  if (is.na(first_num)) stop("No numeric count columns detected.")

  feature_ids <- make.unique(as.character(raw_df[[1L]]))
  count_df <- raw_df[, first_num:ncol(raw_df), drop = FALSE]

  count_mat <- as.matrix(count_df)
  storage.mode(count_mat) <- "numeric"
  rownames(count_mat) <- feature_ids

  count_mat[!is.finite(count_mat)] <- 0
  count_mat <- pmax(count_mat, 0)

  keep <- rowSums(count_mat) > 0
  count_mat[keep, , drop = FALSE]
}

normalize_cpm_log1p <- function(count_mat_arm) {
  lib_sizes <- colSums(count_mat_arm, na.rm = TRUE)
  lib_sizes[!is.finite(lib_sizes) | lib_sizes <= 0] <- 1
  cpm <- sweep(count_mat_arm, 2, lib_sizes / 1e6, "/")
  log1p(cpm)
}

compute_abs_pc1_loadings <- function(norm_mat_arm) {
  pca <- prcomp(t(norm_mat_arm), center = TRUE, scale. = FALSE, rank. = 1)
  out <- abs(pca$rotation[, 1L])
  out[!is.finite(out)] <- 0
  out
}

compute_ranked_variance_curve <- function(metric_mat_arm, rank_order, spar = 0.60) {
  empirical_var <- apply(metric_mat_arm, 1L, stats::var, na.rm = TRUE)
  empirical_var[!is.finite(empirical_var)] <- 0
  empirical_var <- pmax(empirical_var, 0)

  ranked_var <- empirical_var[rank_order]
  ranked_log_var <- log1p(ranked_var)
  ranks <- seq_along(rank_order)

  spline_fit <- stats::smooth.spline(x = ranks, y = ranked_log_var, spar = spar)

  smooth_y <- as.numeric(stats::predict(spline_fit, x = ranks, deriv = 0)$y)
  smooth_d2 <- as.numeric(stats::predict(spline_fit, x = ranks, deriv = 2)$y)

  dense_x <- seq(min(ranks), max(ranks), length.out = length(ranks) * 4L)
  dense_y  <- as.numeric(stats::predict(spline_fit, x = dense_x, deriv = 0)$y)
  dense_d2 <- as.numeric(stats::predict(spline_fit, x = dense_x, deriv = 2)$y)

  out <- data.frame(
    rank = ranks,
    empirical_variance = ranked_var,
    log1p_empirical_variance = ranked_log_var,
    smooth_log1p_empirical_variance = smooth_y,
    d2_spline = smooth_d2,
    stringsAsFactors = FALSE
  )

  attr(out, "dense_curve_df") <- data.frame(
    dense_rank = dense_x,
    dense_smooth_log1p_empirical_variance = dense_y,
    dense_d2 = dense_d2,
    stringsAsFactors = FALSE
  )

  out
}

find_d2_zero_crossings <- function(dense_df) {
  x <- dense_df$dense_rank
  y <- dense_df$dense_d2

  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 2L) {
    return(data.frame(crossing_rank = numeric(0), crossing_type = character(0), stringsAsFactors = FALSE))
  }

  out <- vector("list", 0L)

  for (i in seq_len(length(x) - 1L)) {
    a <- y[i]
    b <- y[i + 1L]
    xa <- x[i]
    xb <- x[i + 1L]

    if (!is.finite(a) || !is.finite(b)) next
    if (a == 0 || b == 0) next

    if ((a < 0 && b > 0) || (a > 0 && b < 0)) {
      frac <- abs(a) / (abs(a) + abs(b))
      xr <- xa + frac * (xb - xa)
      out[[length(out) + 1L]] <- data.frame(
        crossing_rank = xr,
        crossing_type = "sign_change",
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(out) == 0L) {
    return(data.frame(crossing_rank = numeric(0), crossing_type = character(0), stringsAsFactors = FALSE))
  }

  bind_rows(out) %>% distinct() %>% arrange(crossing_rank)
}

select_custom_interval <- function(zero_df, reference_rank, total_n) {
  if (nrow(zero_df) == 0L) stop("No d2 sign-change crossings found.")

  left_candidates  <- zero_df$crossing_rank[zero_df$crossing_rank < reference_rank]
  right_candidates <- zero_df$crossing_rank[zero_df$crossing_rank > reference_rank]

  if (length(left_candidates) == 0L) stop("No left d2 crossing found.")
  if (length(right_candidates) == 0L) stop("No right d2 crossing found.")

  anchor_rank   <- as.integer(round(max(left_candidates)))
  terminal_rank <- as.integer(round(min(right_candidates)))
  reference_rank <- as.integer(round(reference_rank))

  anchor_rank   <- max(1L, anchor_rank)
  terminal_rank <- min(total_n, terminal_rank)

  if (anchor_rank >= reference_rank) stop("Invalid interval: Anchor must lie left of Ref.")
  if (reference_rank >= terminal_rank) stop("Invalid interval: Terminal must lie right of Ref.")

  list(
    anchor = anchor_rank,
    ref = reference_rank,
    terminal = terminal_rank,
    interval_min = anchor_rank,
    interval_max = terminal_rank
  )
}

compute_ranked_feature_metrics <- function(metric_mat_arm, rank_order) {
  ranked_mat <- metric_mat_arm[rank_order, , drop = FALSE]

  mu <- rowMeans(ranked_mat, na.rm = TRUE)
  empirical_var <- apply(ranked_mat, 1L, stats::var, na.rm = TRUE)

  mu[!is.finite(mu)] <- 0
  empirical_var[!is.finite(empirical_var)] <- 0
  mu <- pmax(mu, 0)
  empirical_var <- pmax(empirical_var, 0)

  nb2_var_minus_mu <- pmax(empirical_var - mu, 0)

  alpha_hat <- rep(0, length(mu))
  pos <- mu > 0
  alpha_hat[pos] <- pmax((empirical_var[pos] - mu[pos]) / (mu[pos]^2), 0)
  alpha_mu_val <- alpha_hat * mu

  data.frame(
    rank = seq_along(rank_order),
    feature_id = rownames(metric_mat_arm)[rank_order],
    mu = mu,
    empirical_variance = empirical_var,
    NB2 = log1p(nb2_var_minus_mu),
    NB2_NB1 = log1p(nb2_var_minus_mu) - log1p(mu),
    alpha_mu = log1p(alpha_mu_val),
    stringsAsFactors = FALSE
  )
}

summarize_regions <- function(feature_df, anchor_rank, total_n) {
  right_idx <- seq.int(anchor_rank, total_n)
  right_n <- length(right_idx)

  left_end <- anchor_rank - 1L
  left_start <- left_end - right_n + 1L
  if (left_start < 1L) stop("LEFT block extends below rank 1.")

  left_idx <- seq.int(left_start, left_end)

  left_df <- feature_df[left_idx, , drop = FALSE]
  right_df <- feature_df[right_idx, , drop = FALSE]

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

# =============================================================================
# FORMAL LIKELIHOOD-RATIO TEST: LEFT vs RIGHT NB2 dispersion
# =============================================================================
#
# The NB2/NB2-NB1/alpha*mu quantities above are descriptive medians with no
# associated test. This section adds a formal likelihood-ratio test of
# H0: LEFT and RIGHT share one NB2 dispersion parameter (alpha)
# H1: LEFT and RIGHT have their own separate alpha
#
# Each feature's mean (mu) is held fixed at its own empirical mean, exactly
# as already computed in compute_ranked_feature_metrics; only the dispersion
# parameter alpha is estimated by maximum likelihood, under the standard NB2
# parameterization Var = mu + alpha*mu^2 (equivalently dnbinom size = 1/alpha).
# H1 has exactly one more free parameter than H0 (separate alpha_L, alpha_R
# vs one shared alpha), so:
#   LRT = 2 * (loglik_H1 - loglik_H0)  ~  chi-square(df = 1) under H0
#
# This script applies the formal likelihood-ratio test to raw integer counts only.

region_indices_from_anchor <- function(anchor_rank, total_n) {
  right_idx <- seq.int(anchor_rank, total_n)
  right_n <- length(right_idx)

  left_end <- anchor_rank - 1L
  left_start <- left_end - right_n + 1L
  if (left_start < 1L) stop("LEFT block extends below rank 1.")

  list(left_idx = seq.int(left_start, left_end), right_idx = right_idx)
}

nb2_region_loglik <- function(alpha, counts_block, mu_vec) {
  if (!is.finite(alpha) || alpha <= 0) return(-Inf)
  size <- 1 / alpha
  mu_rep <- rep(mu_vec, times = ncol(counts_block))
  # dnbinom's likelihood is defined for integer counts. Counts are rounded
  # to the nearest non-negative integer so every track uses a genuine
  # discrete NB likelihood.
  counts_vec <- round(pmax(as.vector(counts_block), 0))
  keep <- is.finite(counts_vec) & is.finite(mu_rep) & mu_rep > 0 & counts_vec >= 0
  if (!any(keep)) return(-Inf)
  sum(stats::dnbinom(counts_vec[keep], size = size, mu = mu_rep[keep], log = TRUE))
}

fit_region_alpha_mle <- function(counts_block, mu_vec, log_alpha_lower = -15, log_alpha_upper = 15) {
  # Searching on log(alpha) rather than alpha directly is robust across the
  # wide range of scales the dispersion can plausibly take (roughly 3e-7 to
  # 3e6 over this lower/upper bound).
  obj <- function(log_a) {
    a <- exp(log_a)
    -nb2_region_loglik(a, counts_block, mu_vec)
  }
  opt <- stats::optimize(obj, lower = log_alpha_lower, upper = log_alpha_upper, tol = 1e-8)
  at_bound <- (opt$minimum <= log_alpha_lower + 1e-6) || (opt$minimum >= log_alpha_upper - 1e-6)
  if (at_bound) {
    warning(
      "fit_region_alpha_mle: alpha MLE saturated at the search boundary ",
      "(log_alpha = ", signif(opt$minimum, 4), "). This fit is unreliable; ",
      "see the at_bound columns in the output table."
    )
  }
  list(alpha_mle = exp(opt$minimum), loglik = -opt$objective, at_bound = at_bound)
}

compute_region_lrt <- function(metric_matrix, feature_df, left_idx, right_idx) {
  left_ids <- feature_df$feature_id[left_idx]
  right_ids <- feature_df$feature_id[right_idx]

  left_block <- metric_matrix[left_ids, , drop = FALSE]
  right_block <- metric_matrix[right_ids, , drop = FALSE]

  left_mu <- feature_df$mu[left_idx]
  right_mu <- feature_df$mu[right_idx]

  fit_left <- fit_region_alpha_mle(left_block, left_mu)
  fit_right <- fit_region_alpha_mle(right_block, right_mu)

  pooled_block <- rbind(left_block, right_block)
  pooled_mu <- c(left_mu, right_mu)
  fit_pooled <- fit_region_alpha_mle(pooled_block, pooled_mu)

  loglik_h1 <- fit_left$loglik + fit_right$loglik
  loglik_h0 <- fit_pooled$loglik

  lrt_stat <- 2 * (loglik_h1 - loglik_h0)
  # Numerical optimization can occasionally yield a tiny negative value
  # (H1 should never fit worse than H0 at the true optimum); floor at 0.
  lrt_stat <- max(lrt_stat, 0)

  # Compute the chi-square tail in log space so extremely small p-values are
  # retained rather than silently underflowing to zero.
  lrt_log_p <- stats::pchisq(
    lrt_stat,
    df = 1,
    lower.tail = FALSE,
    log.p = TRUE
  )
  lrt_log10_p <- lrt_log_p / log(10)
  lrt_p <- if (is.finite(lrt_log_p)) {
    if (lrt_log_p > log(.Machine$double.xmin)) {
      exp(lrt_log_p)
    } else {
      0
    }
  } else if (is.infinite(lrt_log_p) && lrt_log_p < 0) {
    0
  } else {
    NA_real_
  }

  out <- data.frame(
    alpha_left_mle = fit_left$alpha_mle,
    alpha_right_mle = fit_right$alpha_mle,
    alpha_pooled_mle = fit_pooled$alpha_mle,
    alpha_left_at_bound = fit_left$at_bound,
    alpha_right_at_bound = fit_right$at_bound,
    alpha_pooled_at_bound = fit_pooled$at_bound,
    loglik_left = fit_left$loglik,
    loglik_right = fit_right$loglik,
    loglik_pooled = fit_pooled$loglik,
    lrt_stat = lrt_stat,
    lrt_df = 1L,
    lrt_p = lrt_p,
    lrt_log_p = lrt_log_p,
    lrt_log10_p = lrt_log10_p,
    # If any of the three fits saturated at the search boundary, the
    # direction call is not meaningful -- report NA rather than a
    # misleading "RIGHT_more_NB2" / "LEFT_more_NB2" label built on a
    # broken estimate.
    lrt_direction = ifelse(
      fit_left$at_bound | fit_right$at_bound | fit_pooled$at_bound,
      NA_character_,
      ifelse(fit_right$alpha_mle > fit_left$alpha_mle, "RIGHT_more_NB2", "LEFT_more_NB2")
    ),
    stringsAsFactors = FALSE
  )

  out
}

validate_method_level <- function(feature_df, interval_info, region_summary, total_n) {
  anchor <- interval_info$anchor
  ref <- interval_info$ref
  terminal <- interval_info$terminal

  if (!(anchor < ref && ref < terminal)) {
    stop("Validation failed: expected Anchor < Ref < Terminal.")
  }

  right_idx <- seq.int(anchor, total_n)
  right_n <- length(right_idx)
  left_end <- anchor - 1L
  left_start <- left_end - right_n + 1L
  if (left_start < 1L) stop("Validation failed: LEFT block extends below rank 1.")

  left_idx <- seq.int(left_start, left_end)
  if (length(left_idx) != length(right_idx)) {
    stop("Validation failed: LEFT and RIGHT blocks do not have equal size.")
  }

  left_df <- feature_df[left_idx, , drop = FALSE]
  right_df <- feature_df[right_idx, , drop = FALSE]

  expected <- list(
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
    diff_alpha = median(right_df$alpha_mu, na.rm = TRUE) - median(left_df$alpha_mu, na.rm = TRUE)
  )

  for (nm in names(expected)) {
    chk <- all.equal(as.numeric(region_summary[[nm]][1]), as.numeric(expected[[nm]]), tolerance = 1e-10)
    if (!isTRUE(chk)) stop("Validation failed for summary field: ", nm, " | ", chk)
  }

  invisible(TRUE)
}


save_three_panel_plot <- function(plot_list, filename) {
  png(filename, width = PNG_WIDTH_IN, height = PNG_HEIGHT_IN, units = "in", res = PNG_DPI, bg = "white")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 3,
    ncol = 1,
    heights = unit(c(1.18, 1.18, 0.88), "null")
  )))
  for (i in seq_along(plot_list)) {
    print(plot_list[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1))
  }
  dev.off()
}



smooth_signal_for_display <- function(rank, value, spar = NB_DISPLAY_SPAR) {
  ok <- is.finite(rank) & is.finite(value)
  out <- rep(NA_real_, length(value))

  if (sum(ok) < 8L) {
    return(out)
  }

  fit <- tryCatch(
    stats::smooth.spline(
      x = rank[ok],
      y = value[ok],
      spar = spar
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(out)
  }

  out[ok] <- as.numeric(
    stats::predict(
      fit,
      x = rank[ok],
      deriv = 0
    )$y
  )

  out
}


format_lrt_p <- function(lrt_p, lrt_log10_p) {
  if (is.finite(lrt_log10_p)) {
    exponent <- floor(lrt_log10_p)
    mantissa <- 10^(lrt_log10_p - exponent)

    if (lrt_log10_p < -4) {
      return(
        paste0(
          formatC(mantissa, format = "f", digits = 2),
          "e",
          exponent
        )
      )
    }
  }

  if (is.finite(lrt_p)) {
    return(formatC(lrt_p, format = "g", digits = 4))
  }

  "NA"
}

# =============================================================================
# FIGURE BUILDER
# =============================================================================

build_main_figure <- function(
    comparison_name,
    arm_label,
    k_star,
    variance_df,
    feature_df,
    interval_info,
    region_summary,
    lrt_result,
    out_file) {

  total_n <- nrow(feature_df)
  anchor <- interval_info$anchor
  kstar_rank <- interval_info$ref
  terminal <- interval_info$terminal

  left_n <- region_summary$left_n
  left_min <- anchor - left_n
  left_max <- anchor - 1L
  right_min <- anchor
  right_max <- total_n

  nb_display <- feature_df %>%
    transmute(
      rank = rank,
      NB2 = smooth_signal_for_display(rank, NB2),
      NB2_NB1 = smooth_signal_for_display(rank, NB2_NB1),
      alpha_mu = smooth_signal_for_display(rank, alpha_mu)
    )

  nb_long <- nb_display %>%
    pivot_longer(
      cols = c(NB2, NB2_NB1, alpha_mu),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c("NB2", "NB2_NB1", "alpha_mu"),
        labels = TRACE_LEVELS
      )
    )

  # -------------------------------------------------------------------------
  # A. k*-anchored raw-count variance geometry.
  # -------------------------------------------------------------------------

  y_top <- max(
    variance_df$smooth_log1p_empirical_variance,
    na.rm = TRUE
  )

  line_df <- data.frame(
    x = c(anchor, kstar_rank, terminal),
    event = c("Anchor", paste0("k*=", k_star), "Terminal"),
    color = c(COL$anchor, COL$kstar, COL$terminal),
    lty = c("solid", "dashed", "dotted"),
    stringsAsFactors = FALSE
  )

  p1 <- ggplot(
    variance_df,
    aes(rank, smooth_log1p_empirical_variance)
  ) +
    annotate(
      "rect",
      xmin = left_min,
      xmax = left_max,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$left_fill,
      alpha = 0.65
    ) +
    annotate(
      "rect",
      xmin = right_min,
      xmax = right_max,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$right_fill,
      alpha = 0.65
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
    )

  for (i in seq_len(nrow(line_df))) {
    p1 <- p1 +
      geom_vline(
        xintercept = line_df$x[i],
        color = line_df$color[i],
        linetype = line_df$lty[i],
        linewidth = 0.9
      )
  }

  # Direct labels replace the repeated multi-item event legend.
  label_y <- y_top * c(0.98, 0.88, 0.78)
  for (i in seq_len(nrow(line_df))) {
    p1 <- p1 +
      annotate(
        "text",
        x = line_df$x[i],
        y = label_y[i],
        label = line_df$event[i],
        color = line_df$color[i],
        angle = 90,
        hjust = 1,
        vjust = -0.35,
        size = 3.0,
        fontface = if (i == 2L) "bold" else "plain"
      )
  }

  p1 <- p1 +
    annotate(
      "text",
      x = (left_min + left_max) / 2,
      y = y_top * 0.98,
      label = paste0("LEFT (n=", region_summary$left_n, ")"),
      vjust = 1,
      size = 3.0,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = (right_min + right_max) / 2,
      y = y_top * 0.98,
      label = paste0("RIGHT (n=", region_summary$right_n, ")"),
      vjust = 1,
      size = 3.0,
      fontface = "bold"
    ) +
    labs(
      title = paste0(
        comparison_name,
        " — ",
        arm_label,
        " | empirical k* = ",
        k_star
      ),
      subtitle = paste0(
        "Raw-count variance across ",
        arm_label,
        " replicates, smoothed along ascending |PC1 loading| rank; ",
        "transition range = Anchor ",
        anchor,
        " to Terminal ",
        terminal
      ),
      x = "PAS rank by ascending absolute PC1 loading",
      y = "Smoothed log1p[Var(raw counts across arm replicates)]"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  # -------------------------------------------------------------------------
  # B. Raw-count NB2-related descriptive signals.
  # -------------------------------------------------------------------------

  p2 <- ggplot() +
    annotate(
      "rect",
      xmin = left_min,
      xmax = left_max,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$left_fill,
      alpha = 0.65
    ) +
    annotate(
      "rect",
      xmin = right_min,
      xmax = right_max,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$right_fill,
      alpha = 0.65
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
    geom_vline(
      xintercept = kstar_rank,
      color = COL$kstar,
      linetype = "dashed",
      linewidth = 0.8
    ) +
    geom_line(
      data = nb_long,
      aes(rank, value, color = metric),
      linewidth = 0.85
    ) +
    scale_color_manual(
      values = TRACE_COLORS,
      breaks = TRACE_LEVELS,
      labels = c(
        "NB2 excess" = "NB2 excess: log1p[max(s²−μ,0)]",
        "NB2-NB1" = "NB2−NB1: excess−log1p(μ)",
        "alpha*mu" = "α̂μ: log1p(α̂μ)"
      )
    ) +
    labs(
      title = "Raw-count NB2 corroboration along the same PC1 rank",
      subtitle = "Display curves are smoothed only for readability; μ and s² are PAS-level raw-count mean and sample variance within this arm",
      x = "PAS rank by ascending absolute PC1 loading",
      y = "Smoothed raw-count NB2 corroboration metric",
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  # -------------------------------------------------------------------------
  # C. Minimal LEFT-vs-RIGHT summary.
  # -------------------------------------------------------------------------

  summary_df <- data.frame(
    region = c("LEFT", "RIGHT"),
    value = c(
      region_summary$left_gap,
      region_summary$right_gap
    ),
    stringsAsFactors = FALSE
  )

  delta_gap <- region_summary$diff_gap
  p_label <- format_lrt_p(
    lrt_result$lrt_p,
    lrt_result$lrt_log10_p
  )

  lrt_direction_text <- if (
    isTRUE(
      lrt_result$alpha_right_mle >
        lrt_result$alpha_left_mle
    )
  ) {
    "αRIGHT > αLEFT"
  } else {
    "αRIGHT ≤ αLEFT"
  }

  summary_subtitle <- paste0(
    "RIGHT−LEFT median (NB2−NB1) = ",
    formatC(delta_gap, format = "f", digits = 3),
    "   |   LEFT-vs-RIGHT NB2-dispersion LRT: χ²(1)=",
    formatC(lrt_result$lrt_stat, format = "f", digits = 2),
    ", p=",
    p_label,
    "   |   ",
    lrt_direction_text
  )

  p3 <- ggplot(summary_df, aes(x = value, y = 1)) +
    geom_segment(
      aes(
        x = summary_df$value[1L],
        xend = summary_df$value[2L],
        y = 1,
        yend = 1
      ),
      color = "grey45",
      linewidth = 0.9
    ) +
    geom_point(
      data = summary_df[summary_df$region == "LEFT", , drop = FALSE],
      aes(x = value, y = 1),
      color = COL$left_pt,
      size = 4
    ) +
    geom_point(
      data = summary_df[summary_df$region == "RIGHT", , drop = FALSE],
      aes(x = value, y = 1),
      color = COL$right_pt,
      size = 4
    ) +
    geom_text(
      data = summary_df[summary_df$region == "LEFT", , drop = FALSE],
      aes(x = value, y = 1, label = paste0("LEFT  ", formatC(value, format = "f", digits = 3))),
      color = COL$left_pt,
      vjust = -1.0,
      hjust = 0.5,
      size = 3.2,
      fontface = "bold"
    ) +
    geom_text(
      data = summary_df[summary_df$region == "RIGHT", , drop = FALSE],
      aes(x = value, y = 1, label = paste0("RIGHT  ", formatC(value, format = "f", digits = 3))),
      color = COL$right_pt,
      vjust = -1.0,
      hjust = 0.5,
      size = 3.2,
      fontface = "bold"
    ) +
    scale_y_continuous(
      breaks = NULL,
      limits = c(0.82, 1.18)
    ) +
    labs(
      title = "LEFT vs RIGHT median NB2−NB1 contrast",
      subtitle = summary_subtitle,
      x = "Median raw-count NB2−NB1 contrast",
      y = NULL,
      caption = paste0(
        "Descriptive Δ uses region medians. LRT p tests H0: one shared NB2 dispersion α for LEFT and RIGHT. ",
        "Reference rank = N−k*+1 = ",
        kstar_rank,
        "."
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.caption = element_text(hjust = 0, size = 8.5)
    )

  save_three_panel_plot(
    list(p1, p2, p3),
    out_file
  )
}

# =============================================================================
# ONE ANALYSIS TRACK
# =============================================================================

run_one_arm <- function(
    comparison_name,
    arm_label,
    k_star,
    count_mat_arm,
    output_dir) {

  # PCA ranking is kept identical to the original Main track:
  # log1p(CPM) expression -> centered PCA -> ascending absolute PC1 loading.
  rank_matrix <- normalize_cpm_log1p(count_mat_arm)
  abs_loadings <- compute_abs_pc1_loadings(rank_matrix)
  rank_order <- order(abs_loadings, decreasing = FALSE)
  total_n <- length(rank_order)

  if (!is.finite(k_star) || k_star < 1L || k_star >= total_n) {
    stop(
      "Invalid empirical k* for ",
      comparison_name,
      ": ",
      k_star,
      " with N=",
      total_n
    )
  }

  # The empirical top-k* begins at this rank on the ascending loading axis.
  reference_rank <- total_n - as.integer(k_star) + 1L

  variance_df <- compute_ranked_variance_curve(
    metric_mat_arm = count_mat_arm,
    rank_order = rank_order,
    spar = VAR_SPLINE_SPAR
  )

  dense_curve_df <- attr(
    variance_df,
    "dense_curve_df"
  )

  zero_df <- find_d2_zero_crossings(
    dense_curve_df
  )

  interval_info <- select_custom_interval(
    zero_df = zero_df,
    reference_rank = reference_rank,
    total_n = total_n
  )

  feature_df <- compute_ranked_feature_metrics(
    metric_mat_arm = count_mat_arm,
    rank_order = rank_order
  )

  feature_df$abs_pc1_loading <- abs_loadings[
    rank_order
  ]

  region_summary <- summarize_regions(
    feature_df = feature_df,
    anchor_rank = interval_info$anchor,
    total_n = total_n
  )

  region_idx <- region_indices_from_anchor(
    anchor_rank = interval_info$anchor,
    total_n = total_n
  )

  # Formal likelihood test is intentionally raw-count only.
  lrt_result <- compute_region_lrt(
    metric_matrix = count_mat_arm,
    feature_df = feature_df,
    left_idx = region_idx$left_idx,
    right_idx = region_idx$right_idx
  )

  validate_method_level(
    feature_df = feature_df,
    interval_info = interval_info,
    region_summary = region_summary,
    total_n = total_n
  )

  selected_df <- data.frame(
    comparison = comparison_name,
    arm = arm_label,
    empirical_k = as.integer(k_star),
    total_ranked_PAS = total_n,
    kstar_reference_rank = reference_rank,
    Anchor = interval_info$anchor,
    Terminal = interval_info$terminal,
    transition_width = interval_info$terminal - interval_info$anchor + 1L,
    left_start = min(region_idx$left_idx),
    left_end = max(region_idx$left_idx),
    right_start = min(region_idx$right_idx),
    right_end = max(region_idx$right_idx),
    stringsAsFactors = FALSE
  )

  summary_row <- bind_cols(
    selected_df,
    region_summary,
    lrt_result
  ) %>%
    mutate(
      right_minus_left_NB2_NB1 = diff_gap,
      lrt_p_label = format_lrt_p(
        lrt_p,
        lrt_log10_p
      )
    )

  rank_path <- file.path(
    output_dir,
    paste0(
      "Table_Rank_",
      comparison_name,
      "_",
      arm_label,
      ".csv"
    )
  )

  summary_path <- file.path(
    output_dir,
    paste0(
      "Table_Summary_",
      comparison_name,
      "_",
      arm_label,
      ".csv"
    )
  )

  fig_path <- file.path(
    output_dir,
    paste0(
      "Figure_",
      comparison_name,
      "_",
      arm_label,
      ".png"
    )
  )

  rank_export <- feature_df %>%
    left_join(
      variance_df,
      by = "rank"
    ) %>%
    mutate(
      empirical_k = as.integer(k_star),
      kstar_reference_rank = reference_rank,
      region = case_when(
        rank %in% region_idx$left_idx ~ "LEFT",
        rank %in% region_idx$right_idx ~ "RIGHT",
        TRUE ~ "OTHER"
      ),
      in_transition_range = (
        rank >= interval_info$anchor &
        rank <= interval_info$terminal
      )
    )

  write.csv(
    rank_export,
    rank_path,
    row.names = FALSE
  )

  write.csv(
    summary_row,
    summary_path,
    row.names = FALSE
  )

  build_main_figure(
    comparison_name = comparison_name,
    arm_label = arm_label,
    k_star = k_star,
    variance_df = variance_df,
    feature_df = feature_df,
    interval_info = interval_info,
    region_summary = region_summary,
    lrt_result = lrt_result,
    out_file = fig_path
  )

  message(
    comparison_name,
    " ",
    arm_label,
    " | k*=",
    k_star,
    " | reference rank=",
    reference_rank,
    " | Anchor=",
    interval_info$anchor,
    " | Terminal=",
    interval_info$terminal,
    " | RIGHT-LEFT NB2-NB1=",
    signif(region_summary$diff_gap, 5),
    " | LRT p=",
    format_lrt_p(
      lrt_result$lrt_p,
      lrt_result$lrt_log10_p
    )
  )

  list(
    summary = summary_row,
    figure = fig_path,
    rank_table = rank_path,
    summary_table = summary_path
  )
}

# =============================================================================
# RUN
# =============================================================================

count_mat <- read_count_matrix(
  COUNT_FILE
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

overall_rows <- list()
figure_paths <- character(0)
table_paths <- character(0)

for (comparison_name in names(COMPARISONS)) {

  k_star <- as.integer(
    EMPIRICAL_K[[comparison_name]]
  )

  comp_dir <- file.path(
    OUT_ROOT,
    comparison_name
  )

  dir.create(
    comp_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  arm_patterns <- COMPARISONS[[
    comparison_name
  ]]

  for (arm_label in names(arm_patterns)) {

    sample_idx <- grep(
      arm_patterns[[arm_label]],
      colnames(count_mat)
    )

    if (length(sample_idx) < 2L) {
      stop(
        "Not enough samples for ",
        comparison_name,
        " ",
        arm_label
      )
    }

    count_mat_arm <- count_mat[
      ,
      sample_idx,
      drop = FALSE
    ]

    result <- run_one_arm(
      comparison_name = comparison_name,
      arm_label = arm_label,
      k_star = k_star,
      count_mat_arm = count_mat_arm,
      output_dir = comp_dir
    )

    overall_rows[[
      length(overall_rows) + 1L
    ]] <- result$summary

    figure_paths <- c(
      figure_paths,
      result$figure
    )

    table_paths <- c(
      table_paths,
      result$rank_table,
      result$summary_table
    )
  }
}

overall_summary <- bind_rows(
  overall_rows
)

# Hard validation: every comparison must carry exactly its prescribed k* in
# both arms; no global/fixed top-k can silently enter the output.
observed_k <- overall_summary %>%
  distinct(
    comparison,
    empirical_k
  )

observed_k <- observed_k[
  match(names(EMPIRICAL_K), observed_k$comparison),
  ,
  drop = FALSE
]

if (
  nrow(observed_k) != length(EMPIRICAL_K) ||
  any(is.na(observed_k$comparison)) ||
  !identical(
    as.integer(observed_k$empirical_k),
    as.integer(EMPIRICAL_K)
  )
) {
  stop(
    "Comparison-specific empirical k* validation failed."
  )
}

overall_path <- file.path(
  OUT_ROOT,
  "Table_Overall_NB1_NB2_Corroboration.csv"
)

lrt_path <- file.path(
  OUT_ROOT,
  "Table_LRT_LEFT_vs_RIGHT_NB2_Dispersion.csv"
)

cutoff_path <- file.path(
  OUT_ROOT,
  "Table_Empirical_kstar_Inputs.csv"
)

write.csv(
  overall_summary,
  overall_path,
  row.names = FALSE
)

write.csv(
  overall_summary %>%
    transmute(
      comparison = comparison,
      arm = arm,
      empirical_k = empirical_k,
      kstar_reference_rank = kstar_reference_rank,
      left_n = left_n,
      right_n = right_n,
      alpha_left_mle = alpha_left_mle,
      alpha_right_mle = alpha_right_mle,
      lrt_stat = lrt_stat,
      lrt_df = lrt_df,
      lrt_p = lrt_p,
      lrt_log10_p = lrt_log10_p,
      lrt_p_label = lrt_p_label,
      direction = lrt_direction,
      median_NB2_NB1_LEFT = left_gap,
      median_NB2_NB1_RIGHT = right_gap,
      right_minus_left_NB2_NB1 = right_minus_left_NB2_NB1
    ),
  lrt_path,
  row.names = FALSE
)

write.csv(
  data.frame(
    comparison = names(EMPIRICAL_K),
    empirical_k = as.integer(EMPIRICAL_K),
    cutoff_scope = "comparison-specific weighted-Pareto k*",
    stringsAsFactors = FALSE
  ),
  cutoff_path,
  row.names = FALSE
)

table_paths <- c(
  table_paths,
  overall_path,
  lrt_path,
  cutoff_path
)

# =============================================================================
# ZIP ARCHIVES
# =============================================================================

message(
  "Creating zip archives..."
)

tryCatch(
  {
    figures_zip_path <- file.path(
      OUT_ROOT,
      "REMNB1_LEADNB2_EmpiricalK_Figures.zip"
    )

    tables_zip_path <- file.path(
      OUT_ROOT,
      "REMNB1_LEADNB2_EmpiricalK_Tables.zip"
    )

    if (file.exists(figures_zip_path)) {
      file.remove(figures_zip_path)
    }

    if (file.exists(tables_zip_path)) {
      file.remove(tables_zip_path)
    }

    if (length(figure_paths) > 0L) {
      utils::zip(
        figures_zip_path,
        files = unique(figure_paths),
        flags = "-j"
      )

      message(
        "  Figures ZIP: ",
        figures_zip_path,
        " (",
        length(unique(figure_paths)),
        " files)"
      )
    }

    if (length(table_paths) > 0L) {
      utils::zip(
        tables_zip_path,
        files = unique(table_paths),
        flags = "-j"
      )

      message(
        "  Tables ZIP: ",
        tables_zip_path,
        " (",
        length(unique(table_paths)),
        " files)"
      )
    }
  },
  error = function(e) {
    warning(
      "ZIP creation failed; individual outputs remain intact: ",
      conditionMessage(e)
    )
  }
)

message(
  "============================================================"
)

message(
  "EMPIRICAL-k* NB1/NB2 CORROBORATION COMPLETE"
)

message(
  "Comparison-specific k*: ",
  paste(
    names(EMPIRICAL_K),
    EMPIRICAL_K,
    sep = "=",
    collapse = "; "
  )
)

message(
  "Formal LRT track: raw integer counts only."
)

message(
  "Outputs: ",
  OUT_ROOT
)

message(
  "============================================================"
)
