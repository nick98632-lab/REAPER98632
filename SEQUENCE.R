#!/usr/bin/env Rscript

# =============================================================================
# MANUSCRIPT-READY EVS / VARIANCE / NB2 SUPPORT FIGURE SCRIPT
# -----------------------------------------------------------------------------
# This version is intentionally simplified for manuscript use.
#
# Core argument implemented by this script:
#
# 1. A geometric transition is defined near the fixed leading-edge-5000 rank.
# 2. The cutoff anchor is the nearest d2 zero immediately LEFT of that fixed
#    leading-edge-5000 rank.
# 3. The terminal start is the nearest d2 zero immediately RIGHT of that fixed
#    leading-edge-5000 rank.
# 4. The full RIGHT leading-edge region is the region from cutoff anchor to the
#    rank end.
# 5. A matched LEFT region of equal size is taken immediately left of the
#    cutoff anchor.
# 6. NB corroboration is then evaluated by comparing RIGHT vs LEFT medians for:
#       NB1 = log(1 + mu)
#       NB2 = log(1 + variance - mu)
#       alpha*mu = log(1 + alpha*mu)
#       NB2 - NB1 contrast
#
# Manuscript design principles used here:
# - The figure is stripped of non-essential panels.
# - Only panels that support the hypothesis are retained.
# - Shading is stronger and cleaner:
#       matched LEFT region = blue
#       full RIGHT leading-edge region = green
#       final geometric interval = darker gray
# - Event markers are offset visually so they do not overlap.
# - True event lines remain at their exact ranks.
# - The derivative panel is compact and may be removed later if not useful.
# - The NB corroboration panel is the main support panel.
#
# Required packages:
#   ggplot2, dplyr, tidyr, grid, patchwork (optional; not required here)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT   <- "/root/REAPER98632/exports/variance_derivative_nb_range_manuscript_final"

FIXED_LEADING_EDGE_SIZE <- 5000L
VARIANCE_SMOOTH_K       <- 101L
DERIVATIVE_SMOOTH_K     <- 51L

PNG_WIDTH_IN  <- 14
PNG_HEIGHT_IN <- 16
PNG_DPI       <- 260

COMPARISONS <- list(
  RT0_ZT6  = list(control = "^R0_", treatment = "^ZT6_"),
  RT2_ZT8  = list(control = "^R2_", treatment = "^ZT8_"),
  RT4_ZT10 = list(control = "^R4_", treatment = "^ZT10_"),
  RT8_ZT14 = list(control = "^R8_", treatment = "^ZT14_")
)

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# COLORS
# =============================================================================

COL <- list(
  loading_line    = "#2C7FB8",
  variance_line   = "#117A65",
  derivative_line = "#386CB0",
  zero_line       = "#D95F02",

  nb1             = "#C99700",
  nb2             = "#1B9E77",
  alpha_mu        = "#386CB0",
  nb_gap          = "#CC1E8C",
  summary         = "#202020",

  left_fill       = "#BFDDF6",
  right_fill      = "#D9F0D3",
  interval_fill   = "#9E9E9E",

  cutoff_anchor   = "#000000",
  fixed_5000      = "#E69F00",
  terminal_start  = "#D95F02",
  left_zero       = "#56B4E9",
  right_zero      = "#CC79A7"
)

EVENT_LEVELS <- c(
  "Cutoff anchor",
  "Fixed leading-edge 5000",
  "Terminal start",
  "Left d2 zero",
  "Right d2 zero"
)

EVENT_COLORS <- c(
  "Cutoff anchor"           = COL$cutoff_anchor,
  "Fixed leading-edge 5000" = COL$fixed_5000,
  "Terminal start"          = COL$terminal_start,
  "Left d2 zero"            = COL$left_zero,
  "Right d2 zero"           = COL$right_zero
)

EVENT_SHAPES <- c(
  "Cutoff anchor"           = 16,
  "Fixed leading-edge 5000" = 18,
  "Terminal start"          = 1,
  "Left d2 zero"            = 15,
  "Right d2 zero"           = 17
)

EVENT_LINE_TYPES <- c(
  "Cutoff anchor"           = "solid",
  "Fixed leading-edge 5000" = "dashed",
  "Terminal start"          = "dotted",
  "Left d2 zero"            = "dashed",
  "Right d2 zero"           = "dashed"
)

NB_LEVELS <- c(
  "NB1 = log(1 + mu)",
  "NB2 = log(1 + variance - mu)",
  "alpha*mu = log(1 + alpha*mu)",
  "NB2 - NB1 contrast",
  "Combined NB support (0 to 1)"
)

NB_COLORS <- c(
  "NB1 = log(1 + mu)"               = COL$nb1,
  "NB2 = log(1 + variance - mu)"    = COL$nb2,
  "alpha*mu = log(1 + alpha*mu)"    = COL$alpha_mu,
  "NB2 - NB1 contrast"              = COL$nb_gap,
  "Combined NB support (0 to 1)"    = COL$summary
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
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(replace(out, is.na(out), 0))
  r <- range(x[ok], na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) {
    out[ok] <- 0.5
    out[!ok] <- 0
    return(out)
  }
  out[ok] <- (x[ok] - r[1]) / (r[2] - r[1])
  out[!ok] <- 0
  out
}

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

normalize_for_ranking <- function(count_mat_arm) {
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

compute_variance_curve <- function(count_mat_arm, rank_order, smooth_k) {
  empirical_var <- apply(count_mat_arm, 1L, stats::var, na.rm = TRUE)
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
  y_sm <- safe_runmed(y, smooth_k)
  d1 <- c(NA_real_, diff(y_sm))
  d1 <- safe_runmed(replace(d1, !is.finite(d1), 0), smooth_k)
  d2 <- c(NA_real_, diff(d1))
  d2 <- safe_runmed(replace(d2, !is.finite(d2), 0), smooth_k)

  data.frame(d1 = d1, d2 = d2)
}

find_zero_crossings <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  out <- list()

  exact_idx <- which(y == 0)
  if (length(exact_idx) > 0L) {
    out[[length(out) + 1L]] <- data.frame(
      crossing_rank = x[exact_idx],
      crossing_type = "exact_zero"
    )
  }

  if (length(x) >= 2L) {
    for (i in seq_len(length(x) - 1L)) {
      yi <- y[i]
      yj <- y[i + 1L]
      xi <- x[i]
      xj <- x[i + 1L]

      if (!is.finite(yi) || !is.finite(yj)) next
      if (yi == 0 || yj == 0) next

      if ((yi < 0 && yj > 0) || (yi > 0 && yj < 0)) {
        frac <- abs(yi) / (abs(yi) + abs(yj))
        xr <- xi + frac * (xj - xi)
        out[[length(out) + 1L]] <- data.frame(
          crossing_rank = xr,
          crossing_type = "sign_change"
        )
      }
    }
  }

  if (length(out) == 0L) {
    return(data.frame(crossing_rank = numeric(0), crossing_type = character(0)))
  }

  bind_rows(out) %>%
    distinct() %>%
    arrange(crossing_rank)
}

select_interval_from_reference <- function(zero_df, reference_rank, total_n) {
  if (nrow(zero_df) == 0L) stop("No d2 zero-crossings found.")

  left_candidates  <- zero_df$crossing_rank[zero_df$crossing_rank < reference_rank]
  right_candidates <- zero_df$crossing_rank[zero_df$crossing_rank > reference_rank]

  if (length(left_candidates) == 0L) stop("No left d2 zero found.")
  if (length(right_candidates) == 0L) stop("No right d2 zero found.")

  cutoff_anchor_rank <- as.integer(round(max(left_candidates)))
  terminal_start_rank <- as.integer(round(min(right_candidates)))
  reference_rank <- as.integer(round(reference_rank))

  cutoff_anchor_rank <- max(1L, cutoff_anchor_rank)
  terminal_start_rank <- min(total_n, terminal_start_rank)

  if (cutoff_anchor_rank >= terminal_start_rank) {
    stop("Invalid interval: cutoff anchor must be left of terminal start.")
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
  empirical_variance <- apply(ranked_counts, 1L, stats::var, na.rm = TRUE)

  mu[!is.finite(mu)] <- 0
  empirical_variance[!is.finite(empirical_variance)] <- 0
  mu <- pmax(mu, 0)
  empirical_variance <- pmax(empirical_variance, 0)

  nb1_mu <- mu
  nb2_var_minus_mu <- pmax(empirical_variance - mu, 0)

  alpha_hat <- rep(0, length(mu))
  pos <- mu > 0
  alpha_hat[pos] <- pmax((empirical_variance[pos] - mu[pos]) / (mu[pos]^2), 0)
  alpha_mu <- alpha_hat * mu

  df <- data.frame(
    rank = seq_along(rank_order),
    feature_id = rownames(count_mat_arm)[rank_order],
    mu = mu,
    empirical_variance = empirical_variance,
    nb1_mu = nb1_mu,
    nb2_variance_minus_mu = nb2_var_minus_mu,
    alpha_hat = alpha_hat,
    alpha_mu = alpha_mu,
    log_nb1 = log1p(nb1_mu),
    log_nb2 = log1p(nb2_var_minus_mu),
    log_alpha_mu = log1p(alpha_mu),
    nb_gap = log1p(nb2_var_minus_mu) - log1p(nb1_mu),
    stringsAsFactors = FALSE
  )

  # Manuscript directional definition:
  # Moving LEFT -> RIGHT:
  #   NB2-like behavior increases when log_nb2, nb_gap, and log_alpha_mu increase.
  #   NB1-like behavior decreases when log_nb1 decreases.
  # Combined support is therefore high when the series is more NB2-like.
  df$combined_nb_support <- rowMeans(
    cbind(
      rescale01(df$log_nb2),
      rescale01(df$nb_gap),
      rescale01(df$log_alpha_mu),
      1 - rescale01(df$log_nb1)
    ),
    na.rm = TRUE
  )

  df
}

summarize_directional_regions <- function(feature_df, cutoff_anchor_rank, total_n) {
  right_idx <- seq.int(cutoff_anchor_rank, total_n)
  right_n <- length(right_idx)

  left_end <- cutoff_anchor_rank - 1L
  left_start <- left_end - right_n + 1L
  if (left_start < 1L) stop("Matched LEFT block extends below rank 1.")

  left_idx <- seq.int(left_start, left_end)

  left_df <- feature_df[left_idx, , drop = FALSE]
  right_df <- feature_df[right_idx, , drop = FALSE]

  data.frame(
    left_n = nrow(left_df),
    right_n = nrow(right_df),

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
}

# =============================================================================
# PLOT HELPERS
# =============================================================================

event_legend_plot <- function() {
  df <- data.frame(
    x = seq_along(EVENT_LEVELS),
    y = 1,
    event = factor(EVENT_LEVELS, levels = EVENT_LEVELS)
  )

  ggplot(df, aes(x, y, color = event, shape = event)) +
    geom_point(size = 3) +
    geom_text(aes(label = event), nudge_y = -0.16, size = 3, show.legend = FALSE) +
    scale_color_manual(values = EVENT_COLORS, drop = FALSE) +
    scale_shape_manual(values = EVENT_SHAPES, drop = FALSE) +
    xlim(0.5, length(EVENT_LEVELS) + 0.5) +
    ylim(0.7, 1.15) +
    theme_void() +
    theme(legend.position = "none")
}

nb_legend_plot <- function() {
  df <- data.frame(
    x = seq_along(NB_LEVELS),
    y = 1,
    metric = factor(NB_LEVELS, levels = NB_LEVELS)
  )

  ggplot(df, aes(x, y, color = metric)) +
    geom_point(size = 3) +
    geom_text(aes(label = metric), nudge_y = -0.16, size = 3, show.legend = FALSE) +
    scale_color_manual(values = NB_COLORS, drop = FALSE) +
    xlim(0.5, length(NB_LEVELS) + 0.5) +
    ylim(0.7, 1.15) +
    theme_void() +
    theme(legend.position = "none")
}

make_event_df <- function(feature_df, variance_df, deriv_df, interval_info) {
  data.frame(
    event = factor(
      c("Cutoff anchor", "Fixed leading-edge 5000", "Terminal start", "Left d2 zero", "Right d2 zero"),
      levels = EVENT_LEVELS
    ),
    rank = c(
      interval_info$cutoff_anchor_rank,
      interval_info$reference_rank,
      interval_info$terminal_start_rank,
      interval_info$cutoff_anchor_rank,
      interval_info$terminal_start_rank
    ),
    loading_y = c(
      feature_df$abs_pc1_loading[interval_info$cutoff_anchor_rank],
      feature_df$abs_pc1_loading[interval_info$reference_rank],
      feature_df$abs_pc1_loading[interval_info$terminal_start_rank],
      feature_df$abs_pc1_loading[interval_info$cutoff_anchor_rank],
      feature_df$abs_pc1_loading[interval_info$terminal_start_rank]
    ),
    variance_y = c(
      variance_df$smooth_log1p_empirical_variance[interval_info$cutoff_anchor_rank],
      variance_df$smooth_log1p_empirical_variance[interval_info$reference_rank],
      variance_df$smooth_log1p_empirical_variance[interval_info$terminal_start_rank],
      variance_df$smooth_log1p_empirical_variance[interval_info$cutoff_anchor_rank],
      variance_df$smooth_log1p_empirical_variance[interval_info$terminal_start_rank]
    ),
    derivative_y = c(
      deriv_df$d2[interval_info$cutoff_anchor_rank],
      deriv_df$d2[interval_info$reference_rank],
      deriv_df$d2[interval_info$terminal_start_rank],
      deriv_df$d2[interval_info$cutoff_anchor_rank],
      deriv_df$d2[interval_info$terminal_start_rank]
    ),
    support_y = c(
      feature_df$log_nb2[interval_info$cutoff_anchor_rank],
      feature_df$log_nb2[interval_info$reference_rank],
      feature_df$log_nb2[interval_info$terminal_start_rank],
      feature_df$log_nb2[interval_info$cutoff_anchor_rank],
      feature_df$log_nb2[interval_info$terminal_start_rank]
    ),
    summary_y = c(
      feature_df$combined_nb_support[interval_info$cutoff_anchor_rank],
      feature_df$combined_nb_support[interval_info$reference_rank],
      feature_df$combined_nb_support[interval_info$terminal_start_rank],
      feature_df$combined_nb_support[interval_info$cutoff_anchor_rank],
      feature_df$combined_nb_support[interval_info$terminal_start_rank]
    ),
    stringsAsFactors = FALSE
  )
}

make_offset_event_df <- function(event_df) {
  event_df %>%
    mutate(
      rank_plot = case_when(
        event == "Cutoff anchor" ~ rank - 18,
        event == "Fixed leading-edge 5000" ~ rank,
        event == "Terminal start" ~ rank + 18,
        event == "Left d2 zero" ~ rank - 9,
        event == "Right d2 zero" ~ rank + 9,
        TRUE ~ rank
      )
    )
}

label_box <- function(x, y, txt) {
  data.frame(x = x, y = y, label = txt, stringsAsFactors = FALSE)
}

save_stacked_plot <- function(plot_list, filename) {
  png(filename, width = PNG_WIDTH_IN, height = PNG_HEIGHT_IN, units = "in", res = PNG_DPI, bg = "white")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = length(plot_list),
    ncol = 1,
    heights = unit(c(1, 0.22, 1, 0.22, 0.85, 0.22, 1.2, 0.22, 0.95), "null")
  )))
  for (i in seq_along(plot_list)) {
    print(plot_list[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1))
  }
  dev.off()
}

# =============================================================================
# MAIN PANEL BUILDER
# =============================================================================

plot_rank_panel <- function(comparison_name, arm_name, feature_df, variance_df, deriv_df,
                            interval_info, region_summary, out_png) {

  total_n <- nrow(feature_df)
  anchor <- interval_info$cutoff_anchor_rank
  ref    <- interval_info$reference_rank
  term   <- interval_info$terminal_start_rank

  left_n <- region_summary$left_n
  left_min <- anchor - left_n
  left_max <- anchor - 1L
  right_min <- anchor
  right_max <- total_n

  event_df <- make_event_df(feature_df, variance_df, deriv_df, interval_info)
  event_plot_df <- make_offset_event_df(event_df)

  vline_df <- data.frame(
    event = factor(
      c("Cutoff anchor", "Fixed leading-edge 5000", "Terminal start"),
      levels = EVENT_LEVELS
    ),
    xint = c(anchor, ref, term),
    stringsAsFactors = FALSE
  )

  left_text <- paste(
    "Absolute loading panel",
    "Leading edge is on the RIGHT",
    "Blue shaded region = matched LEFT comparator",
    "Green shaded region = full RIGHT leading-edge region",
    "Grey band = final geometric interval",
    sep = "\n"
  )

  geom_text <- paste(
    "Variance / d2 method",
    "Variance is empirical log(1 + variance)",
    "d2 zeros come from sign changes or exact zeros only",
    "Cutoff anchor = nearest d2 zero LEFT of fixed leading-edge 5000",
    "Terminal start = nearest d2 zero RIGHT of fixed leading-edge 5000",
    sep = "\n"
  )

  deriv_text <- paste(
    "Derivative support",
    "This panel is kept only to show where the geometric d2 points occur",
    "If this signal is visually uninformative in a final manuscript figure, remove it",
    sep = "\n"
  )

  nb_text <- paste(
    "NB corroboration",
    "RIGHT region = cutoff anchor to rank end",
    "LEFT region = equal-sized matched block immediately left of cutoff anchor",
    "RIGHT > LEFT in log(NB2), NB2-NB1 contrast, and log(alpha*mu)",
    "supports a more NB2-like right leading-edge region",
    sep = "\n"
  )

  right_summary_box <- paste0(
    "Cutoff anchor rank = ", anchor, "\n",
    "Reference rank (5000 from right) = ", ref, "\n",
    "Terminal start rank = ", term, "\n",
    "Final interval = [", anchor, ", ", term, "]\n",
    "Pre-EVS remainder = ", anchor - 1L, "\n",
    "Pre-EVS leading edge = ", total_n - anchor + 1L
  )

  nb_summary_box <- paste0(
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

  support_summary_text <- paste0(
    "NB-supported summary\n",
    "Anchor = ", anchor, "\n",
    "Reference = ", ref, "\n",
    "Terminal = ", term, "\n",
    "LEFT_n = ", left_n, "\n",
    "RIGHT_n = ", region_summary$right_n
  )

  base_bands <- list(
    annotate("rect", xmin = left_min, xmax = left_max, ymin = -Inf, ymax = Inf, fill = COL$left_fill, alpha = 0.72),
    annotate("rect", xmin = right_min, xmax = right_max, ymin = -Inf, ymax = Inf, fill = COL$right_fill, alpha = 0.72),
    annotate("rect", xmin = anchor, xmax = term, ymin = -Inf, ymax = Inf, fill = COL$interval_fill, alpha = 0.42)
  )

  base_lines <- list(
    geom_vline(data = vline_df, aes(xintercept = xint, color = event, linetype = event), linewidth = 0.9, show.legend = FALSE)
  )

  p1 <- ggplot(feature_df, aes(rank, abs_pc1_loading)) +
    base_bands +
    geom_line(color = COL$loading_line, linewidth = 1.0) +
    base_lines +
    geom_point(data = event_plot_df, aes(rank_plot, loading_y, color = event, shape = event), size = 3.0, stroke = 0.9, show.legend = FALSE) +
    geom_label(
      data = label_box(total_n * 0.05, max(feature_df$abs_pc1_loading) * 0.96, left_text),
      aes(x, y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3, linewidth = 0.25, fill = alpha("white", 0.95)
    ) +
    geom_label(
      data = label_box(total_n * 0.72, max(feature_df$abs_pc1_loading) * 0.96, right_summary_box),
      aes(x, y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3, linewidth = 0.25, fill = alpha("white", 0.95)
    ) +
    scale_color_manual(values = EVENT_COLORS, drop = FALSE) +
    scale_shape_manual(values = EVENT_SHAPES, drop = FALSE) +
    scale_linetype_manual(values = EVENT_LINE_TYPES, drop = FALSE) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": absolute PC1 loading series"),
      subtitle = "Leading edge is on the RIGHT",
      x = "EVS rank",
      y = "|PC1 loading|"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  p2 <- ggplot(variance_df, aes(rank, smooth_log1p_empirical_variance)) +
    base_bands +
    geom_line(color = COL$variance_line, linewidth = 1.0) +
    base_lines +
    geom_point(data = event_plot_df, aes(rank_plot, variance_y, color = event, shape = event), size = 3.0, stroke = 0.9, show.legend = FALSE) +
    geom_label(
      data = label_box(total_n * 0.05, max(variance_df$smooth_log1p_empirical_variance) * 0.96, geom_text),
      aes(x, y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3, linewidth = 0.25, fill = alpha("white", 0.95)
    ) +
    scale_color_manual(values = EVENT_COLORS, drop = FALSE) +
    scale_shape_manual(values = EVENT_SHAPES, drop = FALSE) +
    scale_linetype_manual(values = EVENT_LINE_TYPES, drop = FALSE) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": smoothed empirical variance curve"),
      subtitle = "The geometric points are marked directly on the fitted curve",
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  # Compact derivative support only
  p3 <- ggplot(deriv_df, aes(rank, d2)) +
    base_bands +
    geom_hline(yintercept = 0, color = COL$zero_line, linewidth = 0.8) +
    geom_line(color = COL$derivative_line, linewidth = 0.95) +
    base_lines +
    geom_point(data = event_plot_df, aes(rank_plot, derivative_y, color = event, shape = event), size = 3.0, stroke = 0.9, show.legend = FALSE) +
    geom_label(
      data = label_box(total_n * 0.06, max(deriv_df$d2, na.rm = TRUE) * 0.88, deriv_text),
      aes(x, y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3, linewidth = 0.25, fill = alpha("white", 0.95)
    ) +
    scale_color_manual(values = EVENT_COLORS, drop = FALSE) +
    scale_shape_manual(values = EVENT_SHAPES, drop = FALSE) +
    scale_linetype_manual(values = EVENT_LINE_TYPES, drop = FALSE) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": derivative support"),
      subtitle = "The selected d2 zero-crossings around the fixed leading-edge 5000 mark define the geometric interval",
      x = "EVS rank",
      y = "Second derivative value"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank())

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
    annotate("rect", xmin = left_min, xmax = left_max, ymin = -Inf, ymax = Inf, fill = COL$left_fill, alpha = 0.72) +
    annotate("rect", xmin = right_min, xmax = right_max, ymin = -Inf, ymax = Inf, fill = COL$right_fill, alpha = 0.72) +
    annotate("rect", xmin = anchor, xmax = term, ymin = -Inf, ymax = Inf, fill = COL$interval_fill, alpha = 0.42) +
    geom_line(data = nb_long, aes(rank, value, color = metric), linewidth = 0.72) +
    base_lines +
    geom_point(data = event_plot_df, aes(rank_plot, support_y, color = event, shape = event), size = 3.0, stroke = 0.9, show.legend = FALSE) +
    geom_label(
      data = label_box(total_n * 0.05, max(nb_long$value, na.rm = TRUE) * 0.95, nb_text),
      aes(x, y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3, linewidth = 0.25, fill = alpha("white", 0.95)
    ) +
    geom_label(
      data = label_box(total_n * 0.72, max(nb_long$value, na.rm = TRUE) * 0.95, nb_summary_box),
      aes(x, y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.9, linewidth = 0.25, fill = alpha("white", 0.95)
    ) +
    scale_color_manual(
      values = c(
        "NB1 = log(1 + mu)"            = COL$nb1,
        "NB2 = log(1 + variance - mu)" = COL$nb2,
        "alpha*mu = log(1 + alpha*mu)" = COL$alpha_mu,
        "NB2 - NB1 contrast"           = COL$nb_gap,
        EVENT_COLORS
      ),
      drop = FALSE
    ) +
    scale_shape_manual(values = EVENT_SHAPES, drop = FALSE) +
    scale_linetype_manual(values = EVENT_LINE_TYPES, drop = FALSE) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": NB1 / NB2 / alpha*mu support"),
      subtitle = "The full RIGHT leading-edge region is compared against a matched LEFT region of equal size",
      x = "EVS rank",
      y = "Support value"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  p5 <- ggplot(feature_df, aes(rank, combined_nb_support)) +
    annotate("rect", xmin = left_min, xmax = left_max, ymin = -Inf, ymax = Inf, fill = COL$left_fill, alpha = 0.72) +
    annotate("rect", xmin = right_min, xmax = right_max, ymin = -Inf, ymax = Inf, fill = COL$right_fill, alpha = 0.72) +
    annotate("rect", xmin = anchor, xmax = term, ymin = -Inf, ymax = Inf, fill = COL$interval_fill, alpha = 0.42) +
    geom_hline(yintercept = region_summary$nb_support_threshold, color = "#8F8F8F", linetype = "dashed", linewidth = 0.8) +
    geom_line(color = COL$summary, linewidth = 0.9) +
    base_lines +
    geom_point(data = event_plot_df, aes(rank_plot, summary_y, color = event, shape = event), size = 3.0, stroke = 0.9, show.legend = FALSE) +
    geom_label(
      data = label_box(total_n * 0.72, max(feature_df$combined_nb_support, na.rm = TRUE) * 0.95, support_summary_text),
      aes(x, y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3, linewidth = 0.25, fill = alpha("white", 0.95)
    ) +
    scale_color_manual(values = EVENT_COLORS, drop = FALSE) +
    scale_shape_manual(values = EVENT_SHAPES, drop = FALSE) +
    scale_linetype_manual(values = EVENT_LINE_TYPES, drop = FALSE) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": NB-supported summary"),
      subtitle = "The geometric interval remains primary; NB support remains corroborative",
      x = "EVS rank",
      y = "Combined NB support"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  save_stacked_plot(
    list(
      p1,
      event_legend_plot(),
      p2,
      event_legend_plot(),
      p3,
      event_legend_plot(),
      p4,
      nb_legend_plot(),
      p5
    ),
    out_png
  )
}

# =============================================================================
# RUN
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
      stop("Not enough samples for ", comparison_name, " ", arm_name)
    }

    count_mat_arm <- count_mat[, sample_idx, drop = FALSE]
    norm_mat_arm <- normalize_for_ranking(count_mat_arm)

    abs_loadings <- compute_abs_pc1_loadings(norm_mat_arm)
    rank_order <- order(abs_loadings, decreasing = FALSE)
    total_n <- length(rank_order)

    if (FIXED_LEADING_EDGE_SIZE >= total_n) {
      stop("FIXED_LEADING_EDGE_SIZE must be < total_n")
    }

    reference_rank <- total_n - FIXED_LEADING_EDGE_SIZE + 1L

    variance_df <- compute_variance_curve(count_mat_arm, rank_order, VARIANCE_SMOOTH_K)
    deriv_df <- compute_derivatives(variance_df$smooth_log1p_empirical_variance, DERIVATIVE_SMOOTH_K) %>%
      mutate(rank = seq_len(n()))

    zero_df <- find_zero_crossings(deriv_df$rank, deriv_df$d2)

    write.csv(
      zero_df,
      file = file.path(comp_dir, paste0(comparison_name, "_", arm_name, "_zero_crossings_all.csv")),
      row.names = FALSE
    )

    interval_info <- select_interval_from_reference(zero_df, reference_rank, total_n)

    selected_df <- data.frame(
      comparison_name = comparison_name,
      arm = arm_name,
      cutoff_anchor_rank = interval_info$cutoff_anchor_rank,
      fixed_leading_edge_5000_rank = interval_info$reference_rank,
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
        rank, feature_id, abs_pc1_loading,
        mu, empirical_variance, nb1_mu, nb2_variance_minus_mu,
        alpha_hat, alpha_mu, log_nb1, log_nb2, nb_gap,
        log_alpha_mu, combined_nb_support
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

    region_summary <- summarize_directional_regions(feature_df, interval_info$cutoff_anchor_rank, total_n)

    cutoff_summary <- bind_cols(selected_df, region_summary) %>%
      mutate(
        pre_evs_remainder_size = cutoff_anchor_rank - 1L,
        pre_evs_leading_edge_size = total_n - cutoff_anchor_rank + 1L
      ) %>%
      select(
        comparison_name, arm,
        cutoff_anchor_rank, fixed_leading_edge_5000_rank, terminal_start_rank,
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

    plot_rank_panel(
      comparison_name = comparison_name,
      arm_name = arm_name,
      feature_df = feature_df,
      variance_df = variance_df,
      deriv_df = deriv_df,
      interval_info = interval_info,
      region_summary = region_summary,
      out_png = file.path(comp_dir, paste0(comparison_name, "_", arm_name, "_rank_panel.png"))
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
