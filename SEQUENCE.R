#!/usr/bin/env Rscript

# =============================================================================
# VARIANCE / DERIVATIVE / NB-RANGE MANUSCRIPT SCRIPT
# -----------------------------------------------------------------------------
# Purpose
#   This script identifies a late-transition geometric interval in EVS-ranked
#   feature space for each arm of each comparison, using the smoothed empirical
#   variance curve and the second derivative around a fixed leading-edge mark.
#
# Core manuscript method implemented here
#   1. Features are EVS-ranked by absolute PC1 loading within each arm.
#   2. The leading edge is on the RIGHT (highest absolute loading ranks).
#   3. A fixed leading-edge-5000 reference rank is defined as:
#         reference_rank = total_features - 5000 + 1
#   4. The geometric interval is defined only by second-derivative zeros:
#         cutoff anchor  = nearest d2 zero immediately LEFT  of reference_rank
#         terminal start = nearest d2 zero immediately RIGHT of reference_rank
#   5. Variance is empirical per-feature variance across samples in the arm,
#      plotted as smoothed log(1 + variance) across EVS rank.
#   6. Directional NB corroboration is assessed by comparing:
#         RIGHT region = full segment from cutoff anchor to ranked-series end
#         LEFT region  = matched equal-sized block immediately to the left
#      Positive right-minus-left contrasts indicate increasing NB2-like behavior
#      from left to right. Negative contrasts indicate increasing NB1-like
#      behavior from left to right.
#   7. The geometric interval remains primary; NB quantities are corroborative.
#
# Outputs
#   - One folder per comparison with:
#       *_rank_panel.png
#       *_rank_series.csv
#       *_zero_crossings_all.csv
#       *_selected_two_zero_crossings.csv
#       *_feature_level_metrics.csv
#       *_cutoffs_summary.csv
#   - A single overall_cutoff_summary.csv in the output root.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# USER SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT   <- "/root/REAPER98632/exports/variance_derivative_nb_range_manuscript_final"

FIXED_LEADING_EDGE_SIZE <- 5000L
VARIANCE_SMOOTH_K       <- 101L
D2_SMOOTH_K             <- 51L
PNG_WIDTH               <- 2200
PNG_HEIGHT              <- 3200
PNG_RES                 <- 220

# Comparisons are defined by sample-name prefixes.
# Edit these patterns only if your sample names differ.
COMPARISONS <- list(
  RT0_ZT6  = list(control = "^R0_", treatment = "^ZT6_"),
  RT2_ZT8  = list(control = "^R2_", treatment = "^ZT8_"),
  RT4_ZT10 = list(control = "^R4_", treatment = "^ZT10_"),
  RT8_ZT14 = list(control = "^R8_", treatment = "^ZT14_")
)

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# MANUSCRIPT COLOR SYSTEM
# -----------------------------------------------------------------------------
# All figures use the same exact colors and event encodings.
# =============================================================================

COLORS <- list(
  loading_line          = "#1F78B4",
  variance_line         = "#0B6E4F",
  derivative_line       = "#2C7FB8",
  zero_baseline         = "#D95F02",
  nb1_line              = "#D4A017",
  nb2_line              = "#1B9E77",
  alpha_mu_line         = "#386CB0",
  nb_gap_line           = "#C51B7D",
  combined_support_line = "#252525",
  interval_fill         = "#BDBDBD",
  left_region_fill      = "#D9ECFF",
  right_region_fill     = "#D8F5D1",
  cutoff_anchor         = "#000000",
  fixed_le5000          = "#E69F00",
  terminal_start        = "#D95F02",
  d2_zero_left          = "#56B4E9",
  d2_zero_right         = "#CC79A7"
)

EVENT_LEVELS <- c(
  "Cutoff anchor",
  "Fixed leading-edge 5000",
  "Terminal start",
  "Left d2 zero",
  "Right d2 zero"
)

EVENT_COLORS <- c(
  "Cutoff anchor"            = COLORS$cutoff_anchor,
  "Fixed leading-edge 5000"  = COLORS$fixed_le5000,
  "Terminal start"           = COLORS$terminal_start,
  "Left d2 zero"             = COLORS$d2_zero_left,
  "Right d2 zero"            = COLORS$d2_zero_right
)

EVENT_SHAPES <- c(
  "Cutoff anchor"            = 16,
  "Fixed leading-edge 5000"  = 18,
  "Terminal start"           = 1,
  "Left d2 zero"             = 15,
  "Right d2 zero"            = 17
)

NB_METRIC_COLORS <- c(
  "NB1 = log(1 + mu)"                         = COLORS$nb1_line,
  "NB2 = log(1 + variance - mu)"             = COLORS$nb2_line,
  "alpha*mu = log(1 + alpha*mu)"             = COLORS$alpha_mu_line,
  "NB2 - NB1 contrast"                       = COLORS$nb_gap_line,
  "Combined NB support (0 to 1)"             = COLORS$combined_support_line
)

# =============================================================================
# HELPERS
# =============================================================================

safe_runmed <- function(x, k) {
  x <- as.numeric(x)
  if (length(x) < 5L) return(x)

  k <- as.integer(k)
  k <- max(5L, k)
  if ((k %% 2L) == 0L) k <- k + 1L
  if (k >= length(x)) {
    k <- max(5L, 2L * floor((length(x) - 1L) / 2L) + 1L)
  }
  stats::runmed(x, k = k, endrule = "median")
}

rescale01 <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  if (!any(ok)) return(rep(0, length(x)))
  rng <- range(x[ok], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    out <- rep(0.5, length(x))
    out[!ok] <- NA_real_
    return(out)
  }
  out <- (x - rng[1]) / (rng[2] - rng[1])
  out[!ok] <- NA_real_
  out
}

first_numeric_col_index <- function(df) {
  which(vapply(df, is.numeric, logical(1L)))[1L]
}

read_count_matrix <- function(path) {
  raw_df <- read.csv(path, check.names = FALSE)
  stopifnot(nrow(raw_df) > 0, ncol(raw_df) > 1)

  first_num <- first_numeric_col_index(raw_df)
  if (is.na(first_num)) stop("No numeric count columns detected in count file.")

  feature_ids <- raw_df[[1L]]
  count_df <- raw_df[, first_num:ncol(raw_df), drop = FALSE]

  count_mat <- as.matrix(count_df)
  storage.mode(count_mat) <- "numeric"

  rownames(count_mat) <- make.unique(as.character(feature_ids))
  count_mat <- count_mat[rowSums(is.finite(count_mat)) > 0, , drop = FALSE]
  count_mat[!is.finite(count_mat)] <- 0
  count_mat <- pmax(count_mat, 0)

  # Keep all nonzero features; features with all zero counts are uninformative.
  keep <- rowSums(count_mat) > 0
  count_mat <- count_mat[keep, , drop = FALSE]

  count_mat
}

normalize_for_ranking <- function(count_mat_arm) {
  lib_sizes <- colSums(count_mat_arm, na.rm = TRUE)
  lib_sizes[lib_sizes <= 0] <- 1
  cpm <- sweep(count_mat_arm, 2, lib_sizes / 1e6, "/")
  log1p(cpm)
}

compute_abs_pc1_loadings <- function(norm_mat_arm) {
  # samples x features PCA; feature loadings are in rotation[,1]
  pca <- prcomp(t(norm_mat_arm), center = TRUE, scale. = FALSE, rank. = 1)
  loadings <- abs(pca$rotation[, 1L])
  loadings[!is.finite(loadings)] <- 0
  loadings
}

compute_empirical_variance_curve <- function(count_mat_arm, rank_order, smooth_k) {
  empirical_var <- apply(count_mat_arm, 1L, var, na.rm = TRUE)
  empirical_var[!is.finite(empirical_var)] <- 0
  empirical_var <- pmax(empirical_var, 0)

  ranked_var <- empirical_var[rank_order]
  ranked_log_var <- log1p(ranked_var)
  ranked_log_var_smooth <- safe_runmed(ranked_log_var, smooth_k)

  data.frame(
    rank = seq_along(rank_order),
    empirical_variance = ranked_var,
    log1p_empirical_variance = ranked_log_var,
    smooth_log1p_empirical_variance = ranked_log_var_smooth
  )
}

compute_derivatives <- function(y, smooth_k) {
  y <- safe_runmed(y, smooth_k)
  d1 <- c(NA_real_, diff(y))
  d1 <- safe_runmed(ifelse(is.na(d1), 0, d1), smooth_k)
  d2 <- c(NA_real_, diff(d1))
  d2 <- safe_runmed(ifelse(is.na(d2), 0, d2), smooth_k)

  data.frame(
    d1 = d1,
    d2 = d2
  )
}

find_zero_crossings <- function(x, y) {
  stopifnot(length(x) == length(y))
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 2L) {
    return(data.frame(
      crossing_rank = numeric(0),
      type = character(0)
    ))
  }

  crossings <- list()

  # Exact zeros
  exact_idx <- which(y == 0)
  if (length(exact_idx) > 0L) {
    crossings[[length(crossings) + 1L]] <- data.frame(
      crossing_rank = x[exact_idx],
      type = "exact_zero"
    )
  }

  # Sign changes between adjacent points
  s <- sign(y)
  for (i in seq_len(length(y) - 1L)) {
    yi <- y[i]
    yj <- y[i + 1L]
    xi <- x[i]
    xj <- x[i + 1L]

    if (!is.finite(yi) || !is.finite(yj)) next
    if (yi == 0 || yj == 0) next

    if ((yi < 0 && yj > 0) || (yi > 0 && yj < 0)) {
      frac <- abs(yi) / (abs(yi) + abs(yj))
      xr <- xi + frac * (xj - xi)
      crossings[[length(crossings) + 1L]] <- data.frame(
        crossing_rank = xr,
        type = "sign_change"
      )
    }
  }

  if (length(crossings) == 0L) {
    return(data.frame(
      crossing_rank = numeric(0),
      type = character(0)
    ))
  }

  out <- bind_rows(crossings) %>%
    distinct() %>%
    arrange(crossing_rank)

  out
}

select_geometric_interval <- function(zero_df, reference_rank, total_n) {
  if (nrow(zero_df) == 0L) {
    stop("No second-derivative zero crossings found.")
  }

  left_candidates  <- zero_df$crossing_rank[zero_df$crossing_rank < reference_rank]
  right_candidates <- zero_df$crossing_rank[zero_df$crossing_rank > reference_rank]

  if (length(left_candidates) == 0L) {
    stop("No second-derivative zero immediately LEFT of fixed leading-edge-5000 reference.")
  }
  if (length(right_candidates) == 0L) {
    stop("No second-derivative zero immediately RIGHT of fixed leading-edge-5000 reference.")
  }

  cutoff_anchor_rank <- max(left_candidates)
  terminal_start_rank <- min(right_candidates)

  cutoff_anchor_rank <- as.integer(round(cutoff_anchor_rank))
  terminal_start_rank <- as.integer(round(terminal_start_rank))
  reference_rank <- as.integer(round(reference_rank))

  cutoff_anchor_rank <- min(max(cutoff_anchor_rank, 1L), total_n)
  terminal_start_rank <- min(max(terminal_start_rank, 1L), total_n)

  if (cutoff_anchor_rank >= terminal_start_rank) {
    stop("Selected cutoff anchor is not left of selected terminal start.")
  }

  list(
    cutoff_anchor_rank = cutoff_anchor_rank,
    reference_rank = reference_rank,
    terminal_start_rank = terminal_start_rank,
    interval_min_rank = cutoff_anchor_rank,
    interval_max_rank = terminal_start_rank
  )
}

compute_feature_level_metrics <- function(count_mat_arm, rank_order) {
  ranked_counts <- count_mat_arm[rank_order, , drop = FALSE]

  mu <- rowMeans(ranked_counts, na.rm = TRUE)
  empirical_var <- apply(ranked_counts, 1L, var, na.rm = TRUE)

  mu[!is.finite(mu)] <- 0
  empirical_var[!is.finite(empirical_var)] <- 0

  mu <- pmax(mu, 0)
  empirical_var <- pmax(empirical_var, 0)

  nb1_mu <- mu
  nb2_variance_minus_mu <- pmax(empirical_var - mu, 0)

  alpha_hat <- rep(0, length(mu))
  positive_mu <- mu > 0
  alpha_hat[positive_mu] <- pmax((empirical_var[positive_mu] - mu[positive_mu]) / (mu[positive_mu]^2), 0)

  alpha_mu <- alpha_hat * mu

  feature_df <- data.frame(
    rank = seq_along(rank_order),
    feature_id = rownames(count_mat_arm)[rank_order],
    mu = mu,
    empirical_variance = empirical_var,
    nb1_mu = nb1_mu,
    nb2_variance_minus_mu = nb2_variance_minus_mu,
    alpha_hat = alpha_hat,
    alpha_mu = alpha_mu,
    log_nb1 = log1p(nb1_mu),
    log_nb2 = log1p(nb2_variance_minus_mu),
    log_alpha_mu = log1p(alpha_mu),
    nb_gap = log1p(nb2_variance_minus_mu) - log1p(nb1_mu),
    stringsAsFactors = FALSE
  )

  # Combined corroborative support at the feature level
  feature_df$combined_nb_support <- rowMeans(
    cbind(
      rescale01(feature_df$log_nb2),
      rescale01(feature_df$nb_gap),
      rescale01(feature_df$log_alpha_mu)
    ),
    na.rm = TRUE
  )

  feature_df
}

summarize_directional_regions <- function(feature_df, cutoff_anchor_rank, total_n) {
  # RIGHT region is the full leading-edge side from cutoff anchor to the end.
  right_idx <- seq.int(cutoff_anchor_rank, total_n)
  right_n <- length(right_idx)

  # LEFT region is the equal-sized matched block immediately to the left.
  left_end <- cutoff_anchor_rank - 1L
  left_start <- left_end - right_n + 1L
  if (left_start < 1L) {
    stop("Matched left region would extend below rank 1. Cutoff anchor is too far left.")
  }
  left_idx <- seq.int(left_start, left_end)

  left_df  <- feature_df[left_idx,  , drop = FALSE]
  right_df <- feature_df[right_idx, , drop = FALSE]

  out <- data.frame(
    left_n = length(left_idx),
    right_n = length(right_idx),

    left_median_log_nb1 = median(left_df$log_nb1, na.rm = TRUE),
    right_median_log_nb1 = median(right_df$log_nb1, na.rm = TRUE),
    right_left_log_nb1_diff = median(right_df$log_nb1, na.rm = TRUE) - median(left_df$log_nb1, na.rm = TRUE),

    left_median_log_nb2 = median(left_df$log_nb2, na.rm = TRUE),
    right_median_log_nb2 = median(right_df$log_nb2, na.rm = TRUE),
    right_left_log_nb2_diff = median(right_df$log_nb2, na.rm = TRUE) - median(left_df$log_nb2, na.rm = TRUE),

    left_median_nb_gap = median(left_df$nb_gap, na.rm = TRUE),
    right_median_nb_gap = median(right_df$nb_gap, na.rm = TRUE),
    right_left_nb_gap_diff = median(right_df$nb_gap, na.rm = TRUE) - median(left_df$nb_gap, na.rm = TRUE),

    left_median_log_alpha_mu = median(left_df$log_alpha_mu, na.rm = TRUE),
    right_median_log_alpha_mu = median(right_df$log_alpha_mu, na.rm = TRUE),
    right_left_log_alpha_mu_diff = median(right_df$log_alpha_mu, na.rm = TRUE) - median(left_df$log_alpha_mu, na.rm = TRUE),

    center_nb_support = feature_df$combined_nb_support[cutoff_anchor_rank],
    nb_support_threshold = median(feature_df$combined_nb_support, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  out
}

make_event_df <- function(feature_df, variance_df, deriv_df, interval_info) {
  anchor_rank <- interval_info$cutoff_anchor_rank
  ref_rank    <- interval_info$reference_rank
  term_rank   <- interval_info$terminal_start_rank

  data.frame(
    event = factor(
      c("Cutoff anchor", "Fixed leading-edge 5000", "Terminal start", "Left d2 zero", "Right d2 zero"),
      levels = EVENT_LEVELS
    ),
    rank = c(anchor_rank, ref_rank, term_rank, anchor_rank, term_rank),

    loading_y = c(
      feature_df$abs_pc1_loading[anchor_rank],
      feature_df$abs_pc1_loading[ref_rank],
      feature_df$abs_pc1_loading[term_rank],
      feature_df$abs_pc1_loading[anchor_rank],
      feature_df$abs_pc1_loading[term_rank]
    ),

    variance_y = c(
      variance_df$smooth_log1p_empirical_variance[anchor_rank],
      variance_df$smooth_log1p_empirical_variance[ref_rank],
      variance_df$smooth_log1p_empirical_variance[term_rank],
      variance_df$smooth_log1p_empirical_variance[anchor_rank],
      variance_df$smooth_log1p_empirical_variance[term_rank]
    ),

    derivative_y = c(
      deriv_df$d2[anchor_rank],
      deriv_df$d2[ref_rank],
      deriv_df$d2[term_rank],
      deriv_df$d2[anchor_rank],
      deriv_df$d2[term_rank]
    ),

    support_y = c(
      feature_df$nb2_variance_minus_mu[anchor_rank],
      feature_df$nb2_variance_minus_mu[ref_rank],
      feature_df$nb2_variance_minus_mu[term_rank],
      feature_df$nb2_variance_minus_mu[anchor_rank],
      feature_df$nb2_variance_minus_mu[term_rank]
    ),

    summary_y = c(
      feature_df$combined_nb_support[anchor_rank],
      feature_df$combined_nb_support[ref_rank],
      feature_df$combined_nb_support[term_rank],
      feature_df$combined_nb_support[anchor_rank],
      feature_df$combined_nb_support[term_rank]
    )
  )
}

make_method_box_df <- function(total_n, y_top, label_text, x_frac = 0.05) {
  data.frame(
    x = total_n * x_frac,
    y = y_top,
    label = label_text,
    stringsAsFactors = FALSE
  )
}

make_summary_box_df <- function(total_n, y_top, label_text, x_frac = 0.74) {
  data.frame(
    x = total_n * x_frac,
    y = y_top,
    label = label_text,
    stringsAsFactors = FALSE
  )
}

plot_rank_panel <- function(comparison_name,
                            arm_name,
                            feature_df,
                            variance_df,
                            deriv_df,
                            interval_info,
                            region_summary,
                            out_png) {

  total_n <- nrow(feature_df)
  anchor_rank <- interval_info$cutoff_anchor_rank
  ref_rank    <- interval_info$reference_rank
  term_rank   <- interval_info$terminal_start_rank

  left_n <- region_summary$left_n
  right_n <- region_summary$right_n
  left_region_min <- anchor_rank - left_n
  left_region_max <- anchor_rank - 1L
  right_region_min <- anchor_rank
  right_region_max <- total_n

  event_df <- make_event_df(feature_df, variance_df, deriv_df, interval_info)

  loading_method_text <- paste(
    "Absolute loading panel",
    "Leading edge is on the RIGHT",
    "The grey band is the final geometric interval",
    "The legend gives the exact event colors and line types",
    sep = "\n"
  )

  variance_method_text <- paste(
    "Variance curve method",
    "A smoothing spline is fit to empirical log(1 + variance)",
    "Cutoff anchor = nearest d2 zero immediately LEFT of the fixed leading-edge 5000 mark",
    "Terminal start = nearest d2 zero immediately RIGHT of the fixed leading-edge 5000 mark",
    sep = "\n"
  )

  derivative_method_text <- paste(
    "Derivative method",
    "All d2 zeros come only from sign changes or exact zeros",
    "No amplitude threshold is used",
    "No run-length filter is used",
    "The fixed leading-edge 5000 mark lies inside the final interval",
    sep = "\n"
  )

  support_method_text <- paste(
    "NB corroboration",
    "The full RIGHT leading-edge region is the region from cutoff anchor to rank end",
    "The matched LEFT region contains the same number of genes immediately left of the cutoff anchor",
    "These values corroborate, but do not define, the geometric split",
    sep = "\n"
  )

  loading_summary_text <- paste0(
    "Cutoff anchor rank = ", anchor_rank, "\n",
    "Reference rank (5000 from right) = ", ref_rank, "\n",
    "Terminal start rank = ", term_rank, "\n",
    "Final interval = [", interval_info$interval_min_rank, ", ", interval_info$interval_max_rank, "]\n",
    "Pre-EVS remainder = ", anchor_rank - 1L, "\n",
    "Pre-EVS leading edge = ", total_n - anchor_rank + 1L
  )

  support_summary_text <- paste0(
    "LEFT median log(NB2) = ", round(region_summary$left_median_log_nb2, 3), "\n",
    "RIGHT median log(NB2) = ", round(region_summary$right_median_log_nb2, 3), "\n",
    "RIGHT-LEFT log(NB2) diff = ", round(region_summary$right_left_log_nb2_diff, 3), "\n",
    "LEFT median NB2-NB1 contrast = ", round(region_summary$left_median_nb_gap, 3), "\n",
    "RIGHT median NB2-NB1 contrast = ", round(region_summary$right_median_nb_gap, 3), "\n",
    "RIGHT-LEFT contrast diff = ", round(region_summary$right_left_nb_gap_diff, 3), "\n",
    "LEFT median log(alpha*mu) = ", round(region_summary$left_median_log_alpha_mu, 3), "\n",
    "RIGHT median log(alpha*mu) = ", round(region_summary$right_median_log_alpha_mu, 3), "\n",
    "RIGHT-LEFT log(alpha*mu) diff = ", round(region_summary$right_left_log_alpha_mu_diff, 3)
  )

  nb_summary_text <- paste0(
    "NB-supported summary\n",
    "Anchor = ", anchor_rank, "\n",
    "Reference = ", ref_rank, "\n",
    "Terminal = ", term_rank, "\n",
    "LEFT_n = ", left_n, "\n",
    "RIGHT_n = ", right_n
  )

  common_vlines <- list(
    geom_vline(xintercept = anchor_rank, color = COLORS$cutoff_anchor, linewidth = 0.65, linetype = "solid"),
    geom_vline(xintercept = ref_rank,    color = COLORS$fixed_le5000, linewidth = 0.65, linetype = "dashed"),
    geom_vline(xintercept = term_rank,   color = COLORS$terminal_start, linewidth = 0.65, linetype = "dotted")
  )

  common_event_scale <- list(
    scale_color_manual(values = EVENT_COLORS, drop = FALSE),
    scale_shape_manual(values = EVENT_SHAPES, drop = FALSE)
  )

  p1 <- ggplot(feature_df, aes(rank, abs_pc1_loading)) +
    annotate("rect",
             xmin = interval_info$interval_min_rank,
             xmax = interval_info$interval_max_rank,
             ymin = -Inf, ymax = Inf,
             alpha = 0.25, fill = COLORS$interval_fill) +
    geom_line(color = COLORS$loading_line, linewidth = 0.9) +
    common_vlines +
    geom_point(
      data = event_df,
      aes(x = rank, y = loading_y, color = event, shape = event),
      size = 2.6,
      stroke = 0.8,
      inherit.aes = FALSE
    ) +
    geom_label(
      data = make_method_box_df(total_n, max(feature_df$abs_pc1_loading, na.rm = TRUE) * 0.97, loading_method_text),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1, size = 3.0, label.size = 0.25,
      fill = alpha("white", 0.92)
    ) +
    geom_label(
      data = make_summary_box_df(total_n, max(feature_df$abs_pc1_loading, na.rm = TRUE) * 0.97, loading_summary_text),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1, size = 3.0, label.size = 0.25,
      fill = alpha("white", 0.92)
    ) +
    common_event_scale +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": absolute PC1 loading series"),
      subtitle = "Leading edge is on the RIGHT",
      x = "EVS rank",
      y = "|PC1 loading|",
      color = "Event",
      shape = "Event"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  p2 <- ggplot(variance_df, aes(rank, smooth_log1p_empirical_variance)) +
    annotate("rect",
             xmin = interval_info$interval_min_rank,
             xmax = interval_info$interval_max_rank,
             ymin = -Inf, ymax = Inf,
             alpha = 0.25, fill = COLORS$interval_fill) +
    geom_line(color = COLORS$variance_line, linewidth = 0.95) +
    common_vlines +
    geom_point(
      data = event_df,
      aes(x = rank, y = variance_y, color = event, shape = event),
      size = 2.6,
      stroke = 0.8,
      inherit.aes = FALSE
    ) +
    geom_label(
      data = make_method_box_df(total_n, max(variance_df$smooth_log1p_empirical_variance, na.rm = TRUE) * 0.97, variance_method_text),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1, size = 3.0, label.size = 0.25,
      fill = alpha("white", 0.92)
    ) +
    common_event_scale +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": smoothed empirical variance curve"),
      subtitle = "The geometric points are marked directly on the fitted curve",
      x = "EVS rank",
      y = "Fitted log(1 + variance)",
      color = "Event",
      shape = "Event"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  p3 <- ggplot(deriv_df, aes(rank, d2)) +
    annotate("rect",
             xmin = interval_info$interval_min_rank,
             xmax = interval_info$interval_max_rank,
             ymin = -Inf, ymax = Inf,
             alpha = 0.25, fill = COLORS$interval_fill) +
    geom_hline(yintercept = 0, color = COLORS$zero_baseline, linewidth = 0.65) +
    geom_line(color = COLORS$derivative_line, linewidth = 0.9) +
    common_vlines +
    geom_point(
      data = event_df,
      aes(x = rank, y = derivative_y, color = event, shape = event),
      size = 2.8,
      stroke = 0.9,
      inherit.aes = FALSE
    ) +
    geom_label(
      data = make_method_box_df(total_n, max(deriv_df$d2, na.rm = TRUE) * 0.94, derivative_method_text),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1, size = 3.0, label.size = 0.25,
      fill = alpha("white", 0.92)
    ) +
    common_event_scale +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": derivative support"),
      subtitle = "The selected d2 zero-crossings around the fixed leading-edge 5000 mark define the geometric interval",
      x = "EVS rank",
      y = "Second derivative value",
      color = "Event",
      shape = "Event"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  nb_long <- feature_df %>%
    select(rank, log_nb1, log_nb2, log_alpha_mu, nb_gap) %>%
    pivot_longer(
      cols = c(log_nb1, log_nb2, log_alpha_mu, nb_gap),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c("log_nb1", "log_nb2", "log_alpha_mu", "nb_gap"),
        labels = c(
          "NB1 = log(1 + mu)",
          "NB2 = log(1 + variance - mu)",
          "alpha*mu = log(1 + alpha*mu)",
          "NB2 - NB1 contrast"
        )
      )
    )

  p4 <- ggplot() +
    annotate("rect",
             xmin = left_region_min,
             xmax = left_region_max,
             ymin = -Inf, ymax = Inf,
             alpha = 0.22, fill = COLORS$left_region_fill) +
    annotate("rect",
             xmin = right_region_min,
             xmax = right_region_max,
             ymin = -Inf, ymax = Inf,
             alpha = 0.22, fill = COLORS$right_region_fill) +
    annotate("rect",
             xmin = interval_info$interval_min_rank,
             xmax = interval_info$interval_max_rank,
             ymin = -Inf, ymax = Inf,
             alpha = 0.18, fill = COLORS$interval_fill) +
    geom_line(data = nb_long, aes(rank, value, color = metric), linewidth = 0.75) +
    common_vlines +
    geom_point(
      data = event_df,
      aes(x = rank, y = support_y, color = event, shape = event),
      size = 2.8,
      stroke = 0.9,
      inherit.aes = FALSE
    ) +
    geom_label(
      data = make_method_box_df(total_n, max(c(feature_df$log_nb2, feature_df$log_nb1, feature_df$log_alpha_mu, feature_df$nb_gap), na.rm = TRUE) * 0.95, support_method_text),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1, size = 3.0, label.size = 0.25,
      fill = alpha("white", 0.92)
    ) +
    geom_label(
      data = make_summary_box_df(total_n, max(c(feature_df$log_nb2, feature_df$log_nb1, feature_df$log_alpha_mu, feature_df$nb_gap), na.rm = TRUE) * 0.95, support_summary_text),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1, size = 2.9, label.size = 0.25,
      fill = alpha("white", 0.92)
    ) +
    scale_color_manual(
      values = c(NB_METRIC_COLORS[names(NB_METRIC_COLORS) != "Combined NB support (0 to 1)"], EVENT_COLORS),
      breaks = c(names(NB_METRIC_COLORS)[1:4], EVENT_LEVELS)
    ) +
    scale_shape_manual(values = EVENT_SHAPES, breaks = EVENT_LEVELS) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": NB1 / NB2 / alpha*mu support"),
      subtitle = "The full RIGHT leading-edge region is compared against a matched LEFT region of equal size",
      x = "EVS rank",
      y = "Support value",
      color = NULL,
      shape = "Event"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  p5 <- ggplot(feature_df, aes(rank, combined_nb_support)) +
    annotate("rect",
             xmin = left_region_min,
             xmax = left_region_max,
             ymin = -Inf, ymax = Inf,
             alpha = 0.22, fill = COLORS$left_region_fill) +
    annotate("rect",
             xmin = right_region_min,
             xmax = right_region_max,
             ymin = -Inf, ymax = Inf,
             alpha = 0.22, fill = COLORS$right_region_fill) +
    annotate("rect",
             xmin = interval_info$interval_min_rank,
             xmax = interval_info$interval_max_rank,
             ymin = -Inf, ymax = Inf,
             alpha = 0.18, fill = COLORS$interval_fill) +
    geom_hline(yintercept = region_summary$nb_support_threshold, color = "#7F7F7F", linetype = "dashed", linewidth = 0.7) +
    geom_line(color = COLORS$combined_support_line, linewidth = 0.85) +
    common_vlines +
    geom_point(
      data = event_df,
      aes(x = rank, y = summary_y, color = event, shape = event),
      size = 2.8,
      stroke = 0.9,
      inherit.aes = FALSE
    ) +
    geom_label(
      data = make_summary_box_df(total_n, max(feature_df$combined_nb_support, na.rm = TRUE) * 0.97, nb_summary_text),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1, size = 3.0, label.size = 0.25,
      fill = alpha("white", 0.92)
    ) +
    common_event_scale +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": NB-supported summary"),
      subtitle = "The geometric interval remains primary; NB support remains corroborative",
      x = "EVS rank",
      y = "Combined NB support",
      color = "Event",
      shape = "Event"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  g <- (p1 / p2 / p3 / p4 / p5) +
    plot_layout(heights = c(1, 1, 1, 1.2, 1.1), guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(
    filename = out_png,
    plot = g,
    width = PNG_WIDTH,
    height = PNG_HEIGHT,
    units = "px",
    dpi = PNG_RES,
    bg = "white"
  )
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

count_mat <- read_count_matrix(COUNT_FILE)
message("Using count file: ", COUNT_FILE)
message("Count matrix dimensions: ", nrow(count_mat), " features x ", ncol(count_mat), " samples")

overall_rows <- list()

for (comparison_name in names(COMPARISONS)) {
  message("Processing comparison: ", comparison_name)

  comp_dir <- file.path(OUT_ROOT, paste0(comparison_name, "_cutoff_folder"))
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

  pats <- COMPARISONS[[comparison_name]]

  for (arm_name in c("control", "treatment")) {
    message("Selecting ", arm_name, " cutoff for ", comparison_name)

    sample_idx <- grep(pats[[arm_name]], colnames(count_mat))
    if (length(sample_idx) < 2L) {
      stop("Not enough samples detected for ", comparison_name, " ", arm_name,
           ". Pattern used: ", pats[[arm_name]])
    }

    count_mat_arm <- count_mat[, sample_idx, drop = FALSE]
    norm_mat_arm  <- normalize_for_ranking(count_mat_arm)

    abs_loadings <- compute_abs_pc1_loadings(norm_mat_arm)
    rank_order <- order(abs_loadings, decreasing = FALSE)  # leading edge on RIGHT
    total_n <- length(rank_order)

    if (FIXED_LEADING_EDGE_SIZE >= total_n) {
      stop("FIXED_LEADING_EDGE_SIZE is >= total feature count for ", comparison_name, " ", arm_name)
    }

    reference_rank <- total_n - FIXED_LEADING_EDGE_SIZE + 1L

    variance_df <- compute_empirical_variance_curve(
      count_mat_arm = count_mat_arm,
      rank_order = rank_order,
      smooth_k = VARIANCE_SMOOTH_K
    )

    deriv_df <- compute_derivatives(
      y = variance_df$smooth_log1p_empirical_variance,
      smooth_k = D2_SMOOTH_K
    ) %>%
      mutate(rank = seq_len(n()))

    zero_df <- find_zero_crossings(deriv_df$rank, deriv_df$d2)
    write.csv(
      zero_df,
      file = file.path(comp_dir, paste0(comparison_name, "_", arm_name, "_zero_crossings_all.csv")),
      row.names = FALSE
    )

    interval_info <- select_geometric_interval(
      zero_df = zero_df,
      reference_rank = reference_rank,
      total_n = total_n
    )

    selected_df <- data.frame(
      comparison_name = comparison_name,
      arm = arm_name,
      cutoff_anchor_rank = interval_info$cutoff_anchor_rank,
      reference_rank = interval_info$reference_rank,
      terminal_start_rank = interval_info$terminal_start_rank,
      interval_min_rank = interval_info$interval_min_rank,
      interval_max_rank = interval_info$interval_max_rank,
      stringsAsFactors = FALSE
    )

    write.csv(
      selected_df,
      file = file.path(comp_dir, paste0(comparison_name, "_", arm_name, "_selected_two_zero_crossings.csv")),
      row.names = FALSE
    )

    feature_df <- compute_feature_level_metrics(count_mat_arm, rank_order)
    feature_df$abs_pc1_loading <- abs_loadings[rank_order]

    rank_series_df <- feature_df %>%
      select(
        rank, feature_id, abs_pc1_loading, mu, empirical_variance,
        nb1_mu, nb2_variance_minus_mu, alpha_hat, alpha_mu,
        log_nb1, log_nb2, nb_gap, log_alpha_mu, combined_nb_support
      ) %>%
      left_join(variance_df, by = "rank") %>%
      left_join(deriv_df, by = "rank")

    write.csv(
      rank_series_df,
      file = file.path(comp_dir, paste0(comparison_name, "_", arm_name, "_rank_series.csv")),
      row.names = FALSE
    )

    write.csv(
      feature_df,
      file = file.path(comp_dir, paste0(comparison_name, "_", arm_name, "_feature_level_metrics.csv")),
      row.names = FALSE
    )

    region_summary <- summarize_directional_regions(
      feature_df = feature_df,
      cutoff_anchor_rank = interval_info$cutoff_anchor_rank,
      total_n = total_n
    )

    cutoff_summary <- bind_cols(selected_df, region_summary) %>%
      mutate(
        pre_evs_remainder_size = cutoff_anchor_rank - 1L,
        pre_evs_leading_edge_size = total_n - cutoff_anchor_rank + 1L
      ) %>%
      select(
        comparison_name, arm,
        cutoff_anchor_rank, reference_rank, terminal_start_rank,
        interval_min_rank, interval_max_rank,
        pre_evs_remainder_size, pre_evs_leading_edge_size,
        left_n, right_n,
        left_median_log_nb1, right_median_log_nb1, right_left_log_nb1_diff,
        left_median_log_nb2, right_median_log_nb2, right_left_log_nb2_diff,
        left_median_nb_gap, right_median_nb_gap, right_left_nb_gap_diff,
        left_median_log_alpha_mu, right_median_log_alpha_mu, right_left_log_alpha_mu_diff,
        center_nb_support, nb_support_threshold
      )

    write.csv(
      cutoff_summary,
      file = file.path(comp_dir, paste0(comparison_name, "_", arm_name, "_cutoffs_summary.csv")),
      row.names = FALSE
    )

    out_png <- file.path(comp_dir, paste0(comparison_name, "_", arm_name, "_rank_panel.png"))
    plot_rank_panel(
      comparison_name = comparison_name,
      arm_name = arm_name,
      feature_df = feature_df,
      variance_df = variance_df,
      deriv_df = deriv_df,
      interval_info = interval_info,
      region_summary = region_summary,
      out_png = out_png
    )

    message(
      tools::toTitleCase(arm_name), " cutoff anchor rank: ", interval_info$cutoff_anchor_rank,
      " | reference rank: ", interval_info$reference_rank,
      " | terminal start: ", interval_info$terminal_start_rank,
      " | interval: [", interval_info$interval_min_rank, ", ", interval_info$interval_max_rank, "]",
      " | pre-EVS remainder: ", interval_info$cutoff_anchor_rank - 1L,
      " | pre-EVS leading edge: ", total_n - interval_info$cutoff_anchor_rank + 1L
    )

    overall_rows[[length(overall_rows) + 1L]] <- cutoff_summary
  }
}

overall_summary <- bind_rows(overall_rows)
write.csv(
  overall_summary,
  file = file.path(OUT_ROOT, "overall_cutoff_summary.csv"),
  row.names = FALSE
)

message("Done. Outputs written to: ", OUT_ROOT)
