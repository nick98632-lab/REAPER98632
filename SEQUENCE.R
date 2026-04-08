# =============================================================================
# SEQUENCE MANUSCRIPT PIPELINE
# RIGHT-ANCHORED TERMINAL-RISE / SINGLE-NB2-BOOLEAN VERSION
#
# Features are ranked within each comparison and group using the absolute
# magnitude of the first principal-component loading derived from
# log2(x + 1)-transformed normalized counts. Late pre-terminal variance
# geometry is then used to define a compact cutoff band. The terminal
# transition is identified from the right edge inward by tracing the final
# endpoint-directed rise in the smoothed variance curve. The terminal boundary
# is the rightmost eligible pre-terminal zero-crossing immediately preceding
# that final rise. The cutoff center is the nearest eligible zero-crossing to
# its left, yielding the tightest valid late pre-terminal transition band.
#
# Negative-binomial corroboration is summarized numerically on both sides of
# the cutoff. The principal binary corroboration criterion is whether the
# post-cutoff tail region, defined from the cutoff center to the end of the
# ranked series, shows greater median NB2 support than the pre-cutoff region.
# =============================================================================


# =============================================================================
# LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(fdrtool)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(gridExtra)
  library(grid)
  library(scales)
  library(grDevices)
  library(S4Vectors)
  library(BiocParallel)
})


# =============================================================================
# PRIMARY SETTINGS
# =============================================================================

alpha_level <- 0.20
lfc_boundary <- 1.0

fdrtool_pct0 <- 0.75
fdr_clip_floor <- 1e-300
fdr_clip_ceiling <- 0.99
hc_threshold_upper_cap <- 0.95

rolling_window_fraction_d1 <- 0.015
rolling_window_fraction_d2 <- 0.025

terminal_rise_scan_back_fraction <- 0.22
terminal_rise_min_run_fraction   <- 0.035
terminal_rise_pos_frac_min       <- 0.70
terminal_rise_allow_small_neg    <- TRUE
terminal_rise_small_neg_quantile <- 0.25
terminal_rise_drop_tolerance_sd  <- 0.60

center_max_gap_fraction          <- 0.045
center_fallback_max_gap_fraction <- 0.10

nb_support_threshold_quantile    <- 0.35

comparison_pairs <- list(
  RT0_ZT6   = c("R0", "ZT6"),
  RT2_ZT8   = c("R2", "ZT8"),
  RT4_ZT10  = c("R4", "ZT10"),
  RT8_ZT14  = c("R8", "ZT14")
)

input_csv <- "/root/REAPER98632/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
output_dir <- "/root/REAPER98632/exports/variance_derivative_nb_range_final_v3"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# HELPERS
# =============================================================================

clip_probabilities <- function(x,
                               floor_value = fdr_clip_floor,
                               ceiling_value = fdr_clip_ceiling) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  pmin(pmax(x, floor_value), ceiling_value)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_quantile <- function(x, probs, na.rm = TRUE) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = probs, na.rm = na.rm, names = FALSE, type = 8))
}

rolling_mean_centered <- function(x, k) {
  n <- length(x)
  if (n == 0L) return(numeric(0))
  k <- max(1L, as.integer(k))
  if (k <= 1L) return(as.numeric(x))

  half <- floor(k / 2)
  out <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    lo <- max(1L, i - half)
    hi <- min(n, i + half)
    out[i] <- mean(x[lo:hi], na.rm = TRUE)
  }

  out
}

compute_centered_derivatives <- function(x, y,
                                         d1_window_fraction = rolling_window_fraction_d1,
                                         d2_window_fraction = rolling_window_fraction_d2) {
  n <- length(x)

  d1_raw <- c(NA_real_, diff(y) / diff(x))
  if (length(d1_raw) > 1L) d1_raw[1] <- d1_raw[2]

  d1_window <- max(5L, floor(n * d1_window_fraction))
  if (d1_window %% 2L == 0L) d1_window <- d1_window + 1L
  d1 <- rolling_mean_centered(d1_raw, d1_window)

  d2_raw <- c(NA_real_, diff(d1) / diff(x))
  if (length(d2_raw) > 1L) d2_raw[1] <- d2_raw[2]

  d2_window <- max(7L, floor(n * d2_window_fraction))
  if (d2_window %% 2L == 0L) d2_window <- d2_window + 1L
  d2 <- rolling_mean_centered(d2_raw, d2_window)

  list(d1 = d1, d2 = d2)
}

find_zero_crossings_with_sign <- function(y, x = seq_along(y)) {
  keep <- is.finite(y) & is.finite(x)
  y <- y[keep]
  x <- x[keep]

  if (length(y) < 2L) {
    return(data.frame(
      rank = numeric(0),
      idx_left = integer(0),
      idx_right = integer(0),
      y_left = numeric(0),
      y_right = numeric(0),
      crossing_type = character(0),
      stringsAsFactors = FALSE
    ))
  }

  s <- sign(y)

  for (i in seq_along(s)) {
    if (s[i] == 0) {
      left_nonzero <- if (i > 1L) tail(s[seq_len(i - 1L)][s[seq_len(i - 1L)] != 0], 1L) else numeric(0)
      right_nonzero <- if (i < length(s)) head(s[(i + 1L):length(s)][s[(i + 1L):length(s)] != 0], 1L) else numeric(0)

      if (length(left_nonzero)) {
        s[i] <- left_nonzero
      } else if (length(right_nonzero)) {
        s[i] <- right_nonzero
      }
    }
  }

  out <- list()
  j <- 1L

  for (i in seq_len(length(y) - 1L)) {
    if (!is.finite(y[i]) || !is.finite(y[i + 1L])) next

    sign_change <- sign(y[i]) != sign(y[i + 1L]) || y[i] == 0 || y[i + 1L] == 0
    if (!sign_change) next
    if ((y[i + 1L] - y[i]) == 0) next

    x0 <- x[i] - y[i] * (x[i + 1L] - x[i]) / (y[i + 1L] - y[i])

    crossing_type <- if (y[i] < 0 && y[i + 1L] > 0) {
      "neg_to_pos"
    } else if (y[i] > 0 && y[i + 1L] < 0) {
      "pos_to_neg"
    } else {
      "touch_or_flat"
    }

    out[[j]] <- data.frame(
      rank = x0,
      idx_left = i,
      idx_right = i + 1L,
      y_left = y[i],
      y_right = y[i + 1L],
      crossing_type = crossing_type,
      stringsAsFactors = FALSE
    )
    j <- j + 1L
  }

  if (!length(out)) {
    return(data.frame(
      rank = numeric(0),
      idx_left = integer(0),
      idx_right = integer(0),
      y_left = numeric(0),
      y_right = numeric(0),
      crossing_type = character(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, out)
}

compute_empirical_p <- function(wald_z,
                                pct0 = fdrtool_pct0) {
  wald_z <- as.numeric(wald_z)
  keep <- is.finite(wald_z)

  out <- rep(NA_real_, length(wald_z))
  if (!any(keep)) return(out)

  fit <- tryCatch(
    fdrtool::fdrtool(wald_z[keep],
                     statistic = "normal",
                     plot = FALSE,
                     verbose = FALSE,
                     pct0 = pct0),
    error = function(e) NULL
  )

  if (is.null(fit)) return(out)

  if (!is.null(fit$param) && !is.null(fit$param$pval)) {
    out[keep] <- fit$param$pval
  } else if (!is.null(fit$pval)) {
    out[keep] <- fit$pval
  }

  clip_probabilities(out)
}

compute_hc_threshold <- function(empirical_p,
                                 upper_cap = hc_threshold_upper_cap) {
  empirical_p <- empirical_p[is.finite(empirical_p)]
  empirical_p <- empirical_p[empirical_p > 0 & empirical_p < 1]

  if (!length(empirical_p)) return(NA_real_)

  empirical_p <- sort(empirical_p)

  hc <- tryCatch(
    fdrtool::hc.thresh(empirical_p),
    error = function(e) NA_real_
  )

  if (!is.finite(hc)) return(NA_real_)
  min(hc, upper_cap)
}

make_combined_nb_support <- function(log_nb2, log_nb1, log_alpha_mu) {
  nb_gap <- log_nb2 - log_nb1

  z1 <- as.numeric(scale(log_nb2))
  z2 <- as.numeric(scale(nb_gap))
  z3 <- as.numeric(scale(log_alpha_mu))

  combined <- rowMeans(cbind(z1, z2, z3), na.rm = TRUE)
  combined[!is.finite(combined)] <- NA_real_
  combined
}


# =============================================================================
# RIGHT-EDGE TERMINAL-RISE DETECTION
#
# The terminal transition is defined from the right edge inward by identifying
# the last endpoint-directed sustained increase in the smoothed variance curve.
# Earlier local increases are not treated as terminal if they are followed by
# a later decline and renewed ascent. This keeps the terminal geometry tied to
# the final right-edge transition rather than to an earlier interior rise.
# =============================================================================

detect_terminal_rise_start <- function(rank_axis,
                                       fitted_curve,
                                       slope_curve,
                                       support_curve = NULL,
                                       scan_back_fraction = terminal_rise_scan_back_fraction,
                                       min_run_fraction = terminal_rise_min_run_fraction,
                                       pos_frac_min = terminal_rise_pos_frac_min,
                                       allow_small_neg = terminal_rise_allow_small_neg,
                                       small_neg_quantile = terminal_rise_small_neg_quantile,
                                       drop_tolerance_sd = terminal_rise_drop_tolerance_sd) {
  n <- length(rank_axis)
  if (n < 20L) return(n)

  stopifnot(length(fitted_curve) == n, length(slope_curve) == n)

  right_n <- max(25L, floor(n * scan_back_fraction))
  min_run <- max(8L, floor(n * min_run_fraction))

  idx_window <- seq.int(max(1L, n - right_n + 1L), n)
  slope_win  <- slope_curve[idx_window]

  neg_tol <- safe_quantile(abs(slope_win[is.finite(slope_win) & slope_win < 0]), small_neg_quantile)
  if (!is.finite(neg_tol)) neg_tol <- 0

  best_start <- n - min_run + 1L

  for (start in seq.int(n - min_run + 1L, max(1L, n - right_n + 1L), by = -1L)) {
    seg_idx <- seq.int(start, n)
    seg_slope <- slope_curve[seg_idx]
    seg_fit   <- fitted_curve[seg_idx]

    good_pos <- mean(seg_slope > 0, na.rm = TRUE)
    good_net <- isTRUE(tail(seg_fit, 1L) > seg_fit[1L])

    if (allow_small_neg) {
      too_negative <- mean(seg_slope < -neg_tol, na.rm = TRUE)
    } else {
      too_negative <- mean(seg_slope < 0, na.rm = TRUE)
    }

    later_drop <- FALSE
    if (length(seg_fit) >= 5L) {
      running_max <- cummax(seg_fit)
      drop_amount <- running_max - seg_fit
      later_drop <- max(drop_amount, na.rm = TRUE) >
        stats::sd(fitted_curve, na.rm = TRUE) * drop_tolerance_sd
    }

    if (isTRUE(good_net) &&
        is.finite(good_pos) && good_pos >= pos_frac_min &&
        is.finite(too_negative) && too_negative <= (1 - pos_frac_min) &&
        !later_drop) {
      best_start <- start
      break
    }
  }

  best_start
}


# =============================================================================
# RIGHT-ANCHORED ZERO-PAIR SELECTION
#
# The terminal boundary is the rightmost eligible pre-terminal zero-crossing
# immediately preceding the final endpoint-directed rise. The cutoff center is
# the nearest eligible zero-crossing to its left. This yields the tightest valid
# late pre-terminal transition band rather than a broader band anchored by older
# zero-crossings farther to the left.
# =============================================================================

select_right_anchored_zero_pair <- function(rank_axis,
                                            d2_curve,
                                            fitted_curve,
                                            d1_curve,
                                            combined_nb_support,
                                            scan_back_fraction = terminal_rise_scan_back_fraction,
                                            center_max_gap_fraction = center_max_gap_fraction,
                                            center_fallback_max_gap_fraction = center_fallback_max_gap_fraction,
                                            nb_support_threshold_quantile = nb_support_threshold_quantile) {
  n <- length(rank_axis)

  stopifnot(length(d2_curve) == n,
            length(fitted_curve) == n,
            length(d1_curve) == n,
            length(combined_nb_support) == n)

  term_start_idx <- detect_terminal_rise_start(
    rank_axis = rank_axis,
    fitted_curve = fitted_curve,
    slope_curve = d1_curve,
    support_curve = combined_nb_support
  )
  term_start_rank <- rank_axis[term_start_idx]

  zc <- find_zero_crossings_with_sign(d2_curve, rank_axis)

  if (!nrow(zc)) {
    fallback_rank <- rank_axis[max(1L, term_start_idx - 1L)]
    return(list(
      cutoff_center_rank = fallback_rank,
      terminal_start_rank = term_start_rank,
      cutoff_band_min_rank = fallback_rank,
      cutoff_band_max_rank = term_start_rank,
      zero_table = zc
    ))
  }

  support_threshold <- safe_quantile(combined_nb_support, nb_support_threshold_quantile)
  if (!is.finite(support_threshold)) support_threshold <- -Inf

  zc$local_support <- vapply(zc$rank, function(rk) {
    idx <- which.min(abs(rank_axis - rk))
    combined_nb_support[idx]
  }, numeric(1))

  zc$is_preterminal <- zc$rank < term_start_rank
  zc$is_late <- zc$rank >= safe_quantile(rank_axis, 1 - scan_back_fraction)
  zc$is_neg_to_pos <- zc$crossing_type == "neg_to_pos"
  zc$is_supported <- is.finite(zc$local_support) & zc$local_support >= support_threshold

  terminal_candidates <- zc[zc$is_preterminal & zc$is_late & zc$is_neg_to_pos & zc$is_supported, , drop = FALSE]
  if (!nrow(terminal_candidates)) {
    terminal_candidates <- zc[zc$is_preterminal & zc$is_late & zc$is_neg_to_pos, , drop = FALSE]
  }
  if (!nrow(terminal_candidates)) {
    terminal_candidates <- zc[zc$is_preterminal & zc$is_neg_to_pos & zc$is_supported, , drop = FALSE]
  }
  if (!nrow(terminal_candidates)) {
    terminal_candidates <- zc[zc$is_preterminal & zc$is_neg_to_pos, , drop = FALSE]
  }
  if (!nrow(terminal_candidates)) {
    terminal_candidates <- zc[zc$is_preterminal, , drop = FALSE]
  }
  if (!nrow(terminal_candidates)) {
    terminal_candidates <- zc[order(zc$rank, decreasing = TRUE), , drop = FALSE]
  }

  terminal_zero <- terminal_candidates[which.max(terminal_candidates$rank), , drop = FALSE]
  terminal_rank <- terminal_zero$rank[1]

  max_gap <- max(25L, floor(n * center_max_gap_fraction))
  fallback_gap <- max(max_gap + 1L, floor(n * center_fallback_max_gap_fraction))

  earlier <- zc[zc$rank < terminal_rank & zc$is_neg_to_pos, , drop = FALSE]
  earlier <- earlier[order(earlier$rank, decreasing = TRUE), , drop = FALSE]

  close_supported <- earlier[(terminal_rank - earlier$rank) <= max_gap &
                               earlier$local_support >= support_threshold, , drop = FALSE]
  if (!nrow(close_supported)) {
    close_supported <- earlier[(terminal_rank - earlier$rank) <= max_gap, , drop = FALSE]
  }
  if (!nrow(close_supported)) {
    close_supported <- earlier[(terminal_rank - earlier$rank) <= fallback_gap &
                                 earlier$local_support >= support_threshold, , drop = FALSE]
  }
  if (!nrow(close_supported)) {
    close_supported <- earlier[(terminal_rank - earlier$rank) <= fallback_gap, , drop = FALSE]
  }
  if (!nrow(close_supported)) {
    close_supported <- earlier
  }

  if (!nrow(close_supported)) {
    center_rank <- terminal_rank
  } else {
    center_rank <- close_supported$rank[1]
  }

  list(
    cutoff_center_rank = center_rank,
    terminal_start_rank = terminal_rank,
    cutoff_band_min_rank = min(center_rank, terminal_rank),
    cutoff_band_max_rank = max(center_rank, terminal_rank),
    zero_table = zc
  )
}


# =============================================================================
# NB2 TAIL SUMMARY
#
# Negative-binomial corroboration is summarized numerically. The principal
# binary criterion is whether the post-cutoff tail region, defined from the
# cutoff center to the end of the ranked series, exhibits greater median NB2
# support than the pre-cutoff region.
# =============================================================================

build_nb2_tail_summary <- function(rank_axis,
                                   cutoff_center_rank,
                                   terminal_start_rank,
                                   cutoff_band_min_rank,
                                   cutoff_band_max_rank,
                                   pre_evs_leading_edge_size,
                                   pre_evs_remainder_size,
                                   combined_nb_support,
                                   log_nb2,
                                   nb2_minus_nb1,
                                   log_alpha_mu) {
  cutoff_idx <- which(rank_axis >= cutoff_center_rank)
  precut_idx <- which(rank_axis < cutoff_center_rank)

  if (!length(cutoff_idx)) cutoff_idx <- length(rank_axis)
  if (!length(precut_idx)) {
    precut_idx <- seq_len(max(1L, which.min(abs(rank_axis - cutoff_center_rank)) - 1L))
  }

  pre_median_log_nb2 <- safe_median(log_nb2[precut_idx])
  post_median_log_nb2 <- safe_median(log_nb2[cutoff_idx])

  data.frame(
    cutoff_center_rank = cutoff_center_rank,
    terminal_start_rank = terminal_start_rank,
    cutoff_range_rank_min = cutoff_band_min_rank,
    cutoff_range_rank_max = cutoff_band_max_rank,
    pre_evs_leading_edge_size = pre_evs_leading_edge_size,
    pre_evs_remainder_size = pre_evs_remainder_size,

    pre_cutoff_median_log_nb2 = pre_median_log_nb2,
    post_cutoff_median_log_nb2 = post_median_log_nb2,
    post_minus_pre_log_nb2 = post_median_log_nb2 - pre_median_log_nb2,

    pre_cutoff_median_nb2_nb1_contrast = safe_median(nb2_minus_nb1[precut_idx]),
    post_cutoff_median_nb2_nb1_contrast = safe_median(nb2_minus_nb1[cutoff_idx]),
    post_minus_pre_nb2_nb1_contrast =
      safe_median(nb2_minus_nb1[cutoff_idx]) - safe_median(nb2_minus_nb1[precut_idx]),

    pre_cutoff_median_log_alpha_mu = safe_median(log_alpha_mu[precut_idx]),
    post_cutoff_median_log_alpha_mu = safe_median(log_alpha_mu[cutoff_idx]),
    post_minus_pre_log_alpha_mu =
      safe_median(log_alpha_mu[cutoff_idx]) - safe_median(log_alpha_mu[precut_idx]),

    pre_cutoff_mean_combined_nb_support = safe_mean(combined_nb_support[precut_idx]),
    post_cutoff_mean_combined_nb_support = safe_mean(combined_nb_support[cutoff_idx]),
    post_minus_pre_mean_combined_nb_support =
      safe_mean(combined_nb_support[cutoff_idx]) - safe_mean(combined_nb_support[precut_idx]),

    post_cutoff_more_nb2 = post_median_log_nb2 > pre_median_log_nb2,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# PANEL TEXT
# =============================================================================

make_rank_panel_annotation <- function(label_prefix,
                                       cutoff_center_rank,
                                       terminal_start_rank,
                                       cutoff_band_min_rank,
                                       cutoff_band_max_rank,
                                       pre_evs_remainder_size,
                                       pre_evs_leading_edge_size,
                                       tail_summary_row) {
  paste(
    sprintf("%s panel", label_prefix),
    sprintf("Cutoff center rank = %s", cutoff_center_rank),
    sprintf("Cutoff range = [%s, %s]", cutoff_band_min_rank, cutoff_band_max_rank),
    sprintf("Terminal start rank = %s", terminal_start_rank),
    sprintf("Pre-EVS remainder = %s", pre_evs_remainder_size),
    sprintf("Pre-EVS leading edge = %s", pre_evs_leading_edge_size),
    sprintf("Post-cutoff med log(NB2) = %.3f", tail_summary_row$post_cutoff_median_log_nb2),
    sprintf("Pre-cutoff med log(NB2) = %.3f", tail_summary_row$pre_cutoff_median_log_nb2),
    sprintf("Post-cutoff more NB2 = %s",
            ifelse(isTRUE(tail_summary_row$post_cutoff_more_nb2), "TRUE", "FALSE")),
    sep = "\n"
  )
}

make_nb_support_annotation <- function(pre_log_nb2,
                                       post_log_nb2,
                                       pre_nb_gap,
                                       post_nb_gap,
                                       pre_log_alpha_mu,
                                       post_log_alpha_mu) {
  paste(
    "NB corroboration",
    "Post-cutoff tail is evaluated against the pre-cutoff region",
    sprintf("Pre med log(NB2) = %.3f", pre_log_nb2),
    sprintf("Post med log(NB2) = %.3f", post_log_nb2),
    sprintf("Pre med NB2-NB1 contrast = %.3f", pre_nb_gap),
    sprintf("Post med NB2-NB1 contrast = %.3f", post_nb_gap),
    sprintf("Pre med log(alpha*mu) = %.3f", pre_log_alpha_mu),
    sprintf("Post med log(alpha*mu) = %.3f", post_log_alpha_mu),
    sep = "\n"
  )
}


# =============================================================================
# PLOTTING HELPERS
# =============================================================================

plot_rank_panel <- function(df_plot,
                            zero_table,
                            comparison_label,
                            group_label,
                            tail_summary_row,
                            rank_panel_note,
                            nb_support_note,
                            out_file) {
  p1 <- ggplot(df_plot, aes(x = rank, y = abs_pc1_loading)) +
    annotate("rect",
             xmin = tail_summary_row$cutoff_range_rank_min,
             xmax = tail_summary_row$cutoff_range_rank_max,
             ymin = -Inf,
             ymax = Inf,
             alpha = 0.12) +
    geom_line(linewidth = 0.7) +
    geom_vline(xintercept = tail_summary_row$cutoff_center_rank,
               linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = tail_summary_row$terminal_start_rank,
               linetype = "dotted", linewidth = 0.5) +
    annotate("text", x = Inf, y = Inf, label = rank_panel_note,
             hjust = 1.02, vjust = 1.02, size = 3) +
    labs(
      title = paste(comparison_label, group_label, "absolute PC1 loading series"),
      subtitle = "Leading edge is on the RIGHT",
      x = "EVS rank",
      y = "|PC1 loading|"
    ) +
    theme_bw(base_size = 10)

  p2 <- ggplot(df_plot, aes(x = rank, y = fitted_log_variance)) +
    annotate("rect",
             xmin = tail_summary_row$cutoff_range_rank_min,
             xmax = tail_summary_row$cutoff_range_rank_max,
             ymin = -Inf,
             ymax = Inf,
             alpha = 0.12) +
    geom_line(linewidth = 0.7) +
    geom_point(data = df_plot[which.min(abs(df_plot$rank - tail_summary_row$cutoff_center_rank)), , drop = FALSE],
               size = 2) +
    geom_vline(xintercept = tail_summary_row$cutoff_center_rank,
               linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = tail_summary_row$terminal_start_rank,
               linetype = "dotted", linewidth = 0.5) +
    labs(
      title = paste(comparison_label, group_label, "smoothed empirical variance curve"),
      subtitle = "Selected geometric points are marked on the curve",
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    theme_bw(base_size = 10)

  p3 <- ggplot(df_plot, aes(x = rank)) +
    annotate("rect",
             xmin = tail_summary_row$cutoff_range_rank_min,
             xmax = tail_summary_row$cutoff_range_rank_max,
             ymin = -Inf,
             ymax = Inf,
             alpha = 0.12) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_line(aes(y = d2, colour = "Smoothed curvature d2"), linewidth = 0.7) +
    geom_line(aes(y = d1, colour = "Smoothed slope d1"), linewidth = 0.7) +
    geom_vline(xintercept = tail_summary_row$cutoff_center_rank,
               linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = tail_summary_row$terminal_start_rank,
               linetype = "dotted", linewidth = 0.5) +
    labs(
      title = paste(comparison_label, group_label, "derivative support"),
      subtitle = "Center and terminal start come from selected pre-terminal zero-crossings",
      x = "EVS rank",
      y = "Derivative value",
      colour = "colour"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p4_df <- data.frame(
    rank = df_plot$rank,
    combined_nb_support = df_plot$combined_nb_support,
    nb1_mu = df_plot$log_nb1,
    nb2_variance_mu = df_plot$log_nb2,
    alpha_mu = df_plot$log_alpha_mu,
    nb2_minus_nb1 = df_plot$nb2_minus_nb1
  )

  p4 <- ggplot(p4_df, aes(x = rank)) +
    annotate("rect",
             xmin = tail_summary_row$cutoff_range_rank_min,
             xmax = tail_summary_row$cutoff_range_rank_max,
             ymin = -Inf,
             ymax = Inf,
             alpha = 0.12) +
    geom_line(aes(y = combined_nb_support, colour = "Combined NB support"), linewidth = 0.6) +
    geom_line(aes(y = nb1_mu, colour = "NB1 = mu"), linewidth = 0.6) +
    geom_line(aes(y = nb2_variance_mu, colour = "NB2 = variance · mu"), linewidth = 0.6) +
    geom_line(aes(y = alpha_mu, colour = "Smoothed log(alpha*mu)"), linewidth = 0.6) +
    geom_line(aes(y = nb2_minus_nb1, colour = "Smoothed log(NB2+1) - log(NB1+1)"), linewidth = 0.6) +
    geom_vline(xintercept = tail_summary_row$cutoff_center_rank,
               linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = tail_summary_row$terminal_start_rank,
               linetype = "dotted", linewidth = 0.5) +
    annotate("text", x = min(df_plot$rank) + 0.05 * diff(range(df_plot$rank)),
             y = max(p4_df$nb2_variance_mu, na.rm = TRUE) * 0.95,
             label = nb_support_note,
             hjust = 0, vjust = 1, size = 3) +
    labs(
      title = paste(comparison_label, group_label, "NB1 / NB2 / alpha*mu support"),
      subtitle = "Right-of-cutoff elevation in NB2 and alpha*mu supports the leading-edge interpretation",
      x = "EVS rank",
      y = "Support value",
      colour = "colour"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p5 <- ggplot(df_plot, aes(x = rank, y = combined_nb_support)) +
    annotate("rect",
             xmin = tail_summary_row$cutoff_range_rank_min,
             xmax = tail_summary_row$cutoff_range_rank_max,
             ymin = -Inf,
             ymax = Inf,
             alpha = 0.12) +
    geom_hline(yintercept = safe_mean(df_plot$combined_nb_support), linetype = "dotted") +
    geom_line(linewidth = 0.7) +
    geom_vline(xintercept = tail_summary_row$cutoff_center_rank,
               linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = tail_summary_row$terminal_start_rank,
               linetype = "dotted", linewidth = 0.5) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste(
        "Final NB-supported band",
        sprintf("Center = %s", tail_summary_row$cutoff_center_rank),
        sprintf("Band = [%s, %s]",
                tail_summary_row$cutoff_range_rank_min,
                tail_summary_row$cutoff_range_rank_max),
        sprintf("Terminal start = %s", tail_summary_row$terminal_start_rank),
        sep = "\n"
      ),
      hjust = 1.02, vjust = 1.02, size = 3
    ) +
    labs(
      title = paste(comparison_label, group_label, "NB-supported cutoff band"),
      subtitle = "Band expands from the center while combined NB support stays elevated",
      x = "EVS rank",
      y = "Combined NB support"
    ) +
    theme_bw(base_size = 10)

  grob <- gridExtra::arrangeGrob(
    p1, p2, p3, p4, p5,
    ncol = 1
  )

  ggsave(out_file, grob, width = 11, height = 16, dpi = 300)
}


# =============================================================================
# INPUT DATA
# =============================================================================

if (!file.exists(input_csv)) {
  stop(sprintf("Input file not found: %s", input_csv))
}

count_df <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)

feature_col <- colnames(count_df)[1]
count_matrix <- as.matrix(count_df[, -1, drop = FALSE])
rownames(count_matrix) <- count_df[[feature_col]]

sample_names <- colnames(count_matrix)

sample_info <- data.frame(
  sample = sample_names,
  condition = ifelse(grepl("^ZT", sample_names), "trt", "untrt"),
  time_group = sub("_.*$", "", sample_names),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample


# =============================================================================
# ANALYSIS
# =============================================================================

overall_cutoff_summary_list <- list()

for (comparison_label in names(comparison_pairs)) {
  message("Processing comparison: ", comparison_label)

  pair_tokens <- comparison_pairs[[comparison_label]]
  left_prefix <- pair_tokens[1]
  right_prefix <- pair_tokens[2]

  keep_samples <- grepl(paste0("^", left_prefix, "_"), sample_names) |
    grepl(paste0("^", right_prefix, "_"), sample_names)

  counts_sub <- count_matrix[, keep_samples, drop = FALSE]
  coldata_sub <- sample_info[colnames(counts_sub), , drop = FALSE]
  coldata_sub$condition <- factor(coldata_sub$condition, levels = c("untrt", "trt"))

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_sub),
    colData = coldata_sub,
    design = ~ condition
  )

  keep_gene <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_gene, ]

  dds <- DESeq(dds, parallel = FALSE)

  res <- results(dds, alpha = alpha_level)
  coef_name <- grep("condition_trt_vs_untrt", resultsNames(dds), value = TRUE)[1]
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")

  res_df <- as.data.frame(res)
  res_shrunk_df <- as.data.frame(res_shrunk)

  res_df$feature_id <- rownames(res_df)
  res_shrunk_df$feature_id <- rownames(res_shrunk_df)

  merged_df <- res_df %>%
    select(feature_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    rename(
      log2FoldChange_raw = log2FoldChange,
      stat_wald = stat
    ) %>%
    left_join(
      res_shrunk_df %>%
        select(feature_id, log2FoldChange) %>%
        rename(log2FoldChange_shrunk = log2FoldChange),
      by = "feature_id"
    )

  merged_df$empirical_p <- compute_empirical_p(merged_df$stat_wald)
  merged_df$empirical_p <- clip_probabilities(merged_df$empirical_p)

  hc_threshold_dataset <- compute_hc_threshold(merged_df$empirical_p)
  if (!is.finite(hc_threshold_dataset)) hc_threshold_dataset <- hc_threshold_upper_cap

  hbfss_threshold <- abs(log10(hc_threshold_dataset)) * lfc_boundary
  merged_df$HBFSS <- abs(merged_df$log2FoldChange_shrunk * log10(merged_df$empirical_p))

  merged_df$deseq2_significant <- is.finite(merged_df$padj) & merged_df$padj < alpha_level
  merged_df$strong_effect <- merged_df$deseq2_significant &
    is.finite(merged_df$log2FoldChange_shrunk) &
    abs(merged_df$log2FoldChange_shrunk) >= lfc_boundary
  merged_df$weak_effect <- merged_df$deseq2_significant &
    is.finite(merged_df$log2FoldChange_shrunk) &
    abs(merged_df$log2FoldChange_shrunk) < lfc_boundary
  merged_df$hbfss_significant <- is.finite(merged_df$empirical_p) &
    merged_df$empirical_p < hc_threshold_dataset &
    is.finite(merged_df$HBFSS) &
    merged_df$HBFSS >= hbfss_threshold

  norm_counts <- counts(dds, normalized = TRUE)

  utils::write.csv(
    merged_df,
    file = file.path(output_dir, paste0(comparison_label, "_DE_results.csv")),
    row.names = FALSE
  )

  for (group_label in c("control", "treatment")) {
    group_samples <- if (group_label == "control") {
      rownames(coldata_sub)[coldata_sub$condition == "untrt"]
    } else {
      rownames(coldata_sub)[coldata_sub$condition == "trt"]
    }

    group_counts <- norm_counts[, group_samples, drop = FALSE]
    if (ncol(group_counts) < 2L) next

    group_log <- log2(group_counts + 1)

    pca <- prcomp(group_log, center = TRUE, scale. = TRUE)
    abs_pc1_loading <- abs(pca$rotation[, 1])

    ranked_idx <- order(abs_pc1_loading, decreasing = TRUE)
    ranked_feature_id <- rownames(group_log)[ranked_idx]
    rank_axis <- seq_along(ranked_feature_id)

    ranked_var <- apply(group_log[ranked_idx, , drop = FALSE], 1, var, na.rm = TRUE)
    fitted_log_variance <- log1p(ranked_var)

    derivs <- compute_centered_derivatives(rank_axis, fitted_log_variance)
    d1 <- derivs$d1
    d2 <- derivs$d2

    feature_match <- match(ranked_feature_id, merged_df$feature_id)

    ranked_base_mean <- merged_df$baseMean[feature_match]
    ranked_log_nb2 <- log1p(ranked_base_mean + pmax(ranked_var, 0))
    ranked_log_nb1 <- log1p(ranked_base_mean + sqrt(pmax(ranked_base_mean, 0)))
    ranked_log_alpha_mu <- log1p(pmax(ranked_base_mean, 0) * pmax(ranked_var, 0))
    ranked_nb2_minus_nb1 <- ranked_log_nb2 - ranked_log_nb1

    combined_nb_support <- make_combined_nb_support(
      log_nb2 = ranked_log_nb2,
      log_nb1 = ranked_log_nb1,
      log_alpha_mu = ranked_log_alpha_mu
    )

    zero_pair <- select_right_anchored_zero_pair(
      rank_axis = rank_axis,
      d2_curve = d2,
      fitted_curve = fitted_log_variance,
      d1_curve = d1,
      combined_nb_support = combined_nb_support
    )

    cutoff_center_rank <- zero_pair$cutoff_center_rank
    terminal_start_rank <- zero_pair$terminal_start_rank
    cutoff_band_min_rank <- zero_pair$cutoff_band_min_rank
    cutoff_band_max_rank <- zero_pair$cutoff_band_max_rank

    pre_evs_leading_edge_size <- cutoff_center_rank
    pre_evs_remainder_size <- length(rank_axis) - cutoff_center_rank + 1L

    message(
      sprintf(
        "%s %s cutoff center rank: %s | terminal start: %s | cutoff range: [%s, %s] | pre-EVS remainder: %s | pre-EVS leading edge: %s",
        comparison_label,
        group_label,
        cutoff_center_rank,
        terminal_start_rank,
        cutoff_band_min_rank,
        cutoff_band_max_rank,
        pre_evs_remainder_size,
        pre_evs_leading_edge_size
      )
    )

    tail_summary_row <- build_nb2_tail_summary(
      rank_axis = rank_axis,
      cutoff_center_rank = cutoff_center_rank,
      terminal_start_rank = terminal_start_rank,
      cutoff_band_min_rank = cutoff_band_min_rank,
      cutoff_band_max_rank = cutoff_band_max_rank,
      pre_evs_leading_edge_size = pre_evs_leading_edge_size,
      pre_evs_remainder_size = pre_evs_remainder_size,
      combined_nb_support = combined_nb_support,
      log_nb2 = ranked_log_nb2,
      nb2_minus_nb1 = ranked_nb2_minus_nb1,
      log_alpha_mu = ranked_log_alpha_mu
    )

    tail_summary_row$comparison <- comparison_label
    tail_summary_row$group <- group_label
    tail_summary_row$hc_threshold_dataset <- hc_threshold_dataset
    tail_summary_row$hbfss_threshold <- hbfss_threshold

    overall_cutoff_summary_list[[paste(comparison_label, group_label, sep = "__")]] <- tail_summary_row

    rank_panel_note <- make_rank_panel_annotation(
      label_prefix = paste(comparison_label, group_label),
      cutoff_center_rank = cutoff_center_rank,
      terminal_start_rank = terminal_start_rank,
      cutoff_band_min_rank = cutoff_band_min_rank,
      cutoff_band_max_rank = cutoff_band_max_rank,
      pre_evs_remainder_size = pre_evs_remainder_size,
      pre_evs_leading_edge_size = pre_evs_leading_edge_size,
      tail_summary_row = tail_summary_row
    )

    nb_support_note <- make_nb_support_annotation(
      pre_log_nb2 = tail_summary_row$pre_cutoff_median_log_nb2,
      post_log_nb2 = tail_summary_row$post_cutoff_median_log_nb2,
      pre_nb_gap = tail_summary_row$pre_cutoff_median_nb2_nb1_contrast,
      post_nb_gap = tail_summary_row$post_cutoff_median_nb2_nb1_contrast,
      pre_log_alpha_mu = tail_summary_row$pre_cutoff_median_log_alpha_mu,
      post_log_alpha_mu = tail_summary_row$post_cutoff_median_log_alpha_mu
    )

    rank_plot_df <- data.frame(
      feature_id = ranked_feature_id,
      rank = rank_axis,
      abs_pc1_loading = abs_pc1_loading[ranked_idx],
      fitted_log_variance = fitted_log_variance,
      d1 = d1,
      d2 = d2,
      log_nb2 = ranked_log_nb2,
      log_nb1 = ranked_log_nb1,
      nb2_minus_nb1 = ranked_nb2_minus_nb1,
      log_alpha_mu = ranked_log_alpha_mu,
      combined_nb_support = combined_nb_support,
      stringsAsFactors = FALSE
    )

    group_dir <- file.path(output_dir, paste0(comparison_label, "_cutoff_folder"))
    dir.create(group_dir, showWarnings = FALSE, recursive = TRUE)

    utils::write.csv(
      rank_plot_df,
      file = file.path(group_dir, paste0(comparison_label, "_", group_label, "_rank_table.csv")),
      row.names = FALSE
    )

    utils::write.csv(
      zero_pair$zero_table,
      file = file.path(group_dir, paste0(comparison_label, "_", group_label, "_zero_table.csv")),
      row.names = FALSE
    )

    plot_rank_panel(
      df_plot = rank_plot_df,
      zero_table = zero_pair$zero_table,
      comparison_label = comparison_label,
      group_label = group_label,
      tail_summary_row = tail_summary_row,
      rank_panel_note = rank_panel_note,
      nb_support_note = nb_support_note,
      out_file = file.path(group_dir, paste0(comparison_label, "_", group_label, "_rank_panel.png"))
    )
  }
}


# =============================================================================
# FINAL SUMMARY EXPORT
# =============================================================================

overall_cutoff_summary <- dplyr::bind_rows(overall_cutoff_summary_list)

overall_cutoff_summary <- overall_cutoff_summary %>%
  select(
    comparison,
    group,
    cutoff_center_rank,
    terminal_start_rank,
    cutoff_range_rank_min,
    cutoff_range_rank_max,
    pre_evs_leading_edge_size,
    pre_evs_remainder_size,
    pre_cutoff_median_log_nb2,
    post_cutoff_median_log_nb2,
    post_minus_pre_log_nb2,
    pre_cutoff_median_nb2_nb1_contrast,
    post_cutoff_median_nb2_nb1_contrast,
    post_minus_pre_nb2_nb1_contrast,
    pre_cutoff_median_log_alpha_mu,
    post_cutoff_median_log_alpha_mu,
    post_minus_pre_log_alpha_mu,
    pre_cutoff_mean_combined_nb_support,
    post_cutoff_mean_combined_nb_support,
    post_minus_pre_mean_combined_nb_support,
    post_cutoff_more_nb2,
    hc_threshold_dataset,
    hbfss_threshold
  )

utils::write.csv(
  overall_cutoff_summary,
  file = file.path(output_dir, "overall_cutoff_summary.csv"),
  row.names = FALSE
)


# =============================================================================
# MANUSCRIPT METHODS TEXT MIRROR
# =============================================================================
# Features were ranked within each comparison and group using the absolute
# magnitude of the first principal-component loading derived from
# log2(x + 1)-transformed normalized counts. The late pre-terminal cutoff
# geometry was defined on the smoothed empirical variance curve by tracing the
# final endpoint-directed rise from the right edge inward. The terminal boundary
# was taken as the rightmost eligible pre-terminal zero-crossing of the smoothed
# second derivative immediately preceding this final rise. The cutoff center was
# taken as the nearest eligible zero-crossing to the left of that terminal
# boundary, yielding the tightest valid late pre-terminal transition band.
# Negative-binomial corroboration was summarized numerically on both sides of
# the cutoff using median log(NB2), median NB2 minus NB1 contrast, median
# log(alpha·mu), and mean combined NB support. The principal binary criterion
# was whether the post-cutoff tail region, defined from the cutoff center to the
# end of the ranked series, showed greater median NB2 support than the
# pre-cutoff region.
# =============================================================================
