#!/usr/bin/env Rscript

# =============================================================================
# MANUSCRIPT SCRIPT
# CUSTOM GEOMETRY + LITERATURE-SAFE NB CORROBORATION
# -----------------------------------------------------------------------------
# LITERATURE-SAFE FRAMING
#
# This script does NOT claim that the geometric cutoff rule is a standard
# published negative-binomial method.
#
# Instead, it implements:
#
# 1. A CUSTOM geometry-based transition rule on the EVS-ranked empirical
#    variance curve.
#
# 2. A LITERATURE-MOTIVATED corroboration step using NB-related quantities
#    derived from the count mean-variance relationship:
#       NB1-like lower-order mean-linked term      ~ mu
#       NB2-like higher-order variance-linked term ~ variance - mu
#       alpha*mu term                              ~ extra-Poisson dispersion
#
# 3. A matched left-vs-right comparison:
#       RIGHT region = full leading-edge block from cutoff anchor to rank end
#       LEFT region  = equal-sized matched block immediately left of cutoff
#
# 4. Positive RIGHT-minus-LEFT differences in:
#       log(1 + variance - mu)
#       NB2 - NB1 contrast
#       log(1 + alpha*mu)
#    are interpreted as stronger NB2-like behavior toward the leading edge.
#
# MAIN FIGURE DESIGN
#
# Panel A: Custom geometric transition on the smoothed empirical variance curve
# Panel B: NB corroboration panel with only interpretable NB2-related traces
# Panel C: Compact left-vs-right median summary
#
# The derivative itself is used internally for the custom geometric rule but is
# not shown as a main manuscript panel unless explicitly requested later.
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
OUT_ROOT   <- "/root/REAPER98632/exports/manuscript_custom_geometry_nb2_support"

FIXED_LEADING_EDGE_SIZE <- 5000L
VARIANCE_SMOOTH_K       <- 101L
DERIVATIVE_SMOOTH_K     <- 51L

PNG_WIDTH_IN  <- 14
PNG_HEIGHT_IN <- 12
PNG_DPI       <- 260

# Edit patterns only if your sample names differ.
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
  variance_curve  = "#117A65",
  nb2_line        = "#1B9E77",
  nb_gap_line     = "#CC1E8C",
  alpha_mu_line   = "#386CB0",
  nb1_line        = "#C99700",

  left_fill       = "#CBE3F8",
  right_fill      = "#DDF2D5",
  interval_fill   = "#B0B0B0",

  cutoff_anchor   = "#000000",
  fixed_5000      = "#E69F00",
  terminal_start  = "#D95F02",

  summary_left    = "#5B8FD1",
  summary_right   = "#4DAF4A"
)

EVENT_LEVELS <- c(
  "Cutoff anchor",
  "Fixed leading-edge 5000",
  "Terminal start"
)

EVENT_COLORS <- c(
  "Cutoff anchor"           = COL$cutoff_anchor,
  "Fixed leading-edge 5000" = COL$fixed_5000,
  "Terminal start"          = COL$terminal_start
)

EVENT_SHAPES <- c(
  "Cutoff anchor"           = 16,
  "Fixed leading-edge 5000" = 18,
  "Terminal start"          = 1
)

EVENT_LINE_TYPES <- c(
  "Cutoff anchor"           = "solid",
  "Fixed leading-edge 5000" = "dashed",
  "Terminal start"          = "dotted"
)

NB_TRACE_COLORS <- c(
  "NB2 = log(1 + variance - mu)" = COL$nb2_line,
  "NB2 - NB1 contrast"           = COL$nb_gap_line,
  "alpha*mu = log(1 + alpha*mu)" = COL$alpha_mu_line
)

SUMMARY_METRIC_LEVELS <- c(
  "NB1 = log(1 + mu)",
  "NB2 = log(1 + variance - mu)",
  "NB2 - NB1 contrast",
  "alpha*mu = log(1 + alpha*mu)"
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
  out <- rep(0, length(x))
  if (!any(ok)) return(out)

  r <- range(x[ok], na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) {
    out[ok] <- 0.5
    return(out)
  }

  out[ok] <- (x[ok] - r[1]) / (r[2] - r[1])
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
      crossing_type = "exact_zero",
      stringsAsFactors = FALSE
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
          crossing_type = "sign_change",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(out) == 0L) {
    return(data.frame(
      crossing_rank = numeric(0),
      crossing_type = character(0),
      stringsAsFactors = FALSE
    ))
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

  # Higher combined support = more NB2-like to the right
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

make_event_df <- function(feature_df, variance_df, interval_info) {
  data.frame(
    event = factor(EVENT_LEVELS, levels = EVENT_LEVELS),
    rank = c(
      interval_info$cutoff_anchor_rank,
      interval_info$reference_rank,
      interval_info$terminal_start_rank
    ),
    variance_y = c(
      variance_df$smooth_log1p_empirical_variance[interval_info$cutoff_anchor_rank],
      variance_df$smooth_log1p_empirical_variance[interval_info$reference_rank],
      variance_df$smooth_log1p_empirical_variance[interval_info$terminal_start_rank]
    ),
    stringsAsFactors = FALSE
  )
}

make_offset_event_df <- function(event_df, total_n) {
  offset_big <- max(8L, round(total_n * 0.0025))
  event_df %>%
    mutate(
      rank_plot = case_when(
        event == "Cutoff anchor" ~ rank - offset_big,
        event == "Fixed leading-edge 5000" ~ rank,
        event == "Terminal start" ~ rank + offset_big,
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
    heights = unit(c(1.15, 0.16, 1.15, 0.16, 0.95), "null")
  )))
  for (i in seq_along(plot_list)) {
    print(plot_list[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1))
  }
  dev.off()
}

# =============================================================================
# FIGURE BUILDER
# =============================================================================

plot_rank_panel <- function(comparison_name, arm_name, feature_df, variance_df,
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

  event_df <- data.frame(
    event = factor(
      c("Cutoff anchor", "Fixed leading-edge 5000", "Terminal start"),
      levels = EVENT_LEVELS
    ),
    rank = c(anchor, ref, term),
    variance_y = c(
      variance_df$smooth_log1p_empirical_variance[anchor],
      variance_df$smooth_log1p_empirical_variance[ref],
      variance_df$smooth_log1p_empirical_variance[term]
    ),
    stringsAsFactors = FALSE
  )

  event_plot_df <- make_offset_event_df(event_df, total_n)

  vline_df <- data.frame(
    event = factor(c("Cutoff anchor", "Fixed leading-edge 5000", "Terminal start"), levels = EVENT_LEVELS),
    xint = c(anchor, ref, term),
    stringsAsFactors = FALSE
  )

  geom_text <- paste(
    "Custom geometry panel",
    "Blue region = matched LEFT comparator",
    "Green region = full RIGHT leading-edge block",
    "Grey band = final geometric interval",
    "Cutoff anchor and terminal start are custom d2-zero landmarks around the fixed leading-edge 5000 mark",
    sep = "\n"
  )

  geom_summary <- paste0(
    "Cutoff anchor rank = ", anchor, "\n",
    "Reference rank (5000 from right) = ", ref, "\n",
    "Terminal start rank = ", term, "\n",
    "Final interval = [", anchor, ", ", term, "]\n",
    "Pre-EVS remainder = ", anchor - 1L, "\n",
    "Pre-EVS leading edge = ", total_n - anchor + 1L
  )

  nb_text <- paste(
    "NB corroboration panel",
    "RIGHT region = cutoff anchor to rank end",
    "LEFT region = equal-sized matched block immediately left of cutoff anchor",
    "Higher right-side NB2, NB2-NB1 contrast, and alpha*mu support a more NB2-like leading edge",
    sep = "\n"
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

  base_bands <- list(
    annotate("rect", xmin = left_min, xmax = left_max, ymin = -Inf, ymax = Inf, fill = COL$left_fill, alpha = 0.72),
    annotate("rect", xmin = right_min, xmax = right_max, ymin = -Inf, ymax = Inf, fill = COL$right_fill, alpha = 0.72),
    annotate("rect", xmin = anchor, xmax = term, ymin = -Inf, ymax = Inf, fill = COL$interval_fill, alpha = 0.42)
  )

  base_lines <- list(
    geom_vline(data = vline_df, aes(xintercept = xint, color = event, linetype = event), linewidth = 0.9, show.legend = FALSE)
  )

  # Panel A: geometry only
  p1 <- ggplot(variance_df, aes(rank, smooth_log1p_empirical_variance)) +
    base_bands +
    geom_line(color = COL$variance_curve, linewidth = 1.0) +
    base_lines +
    geom_point(
      data = event_plot_df,
      aes(rank_plot, variance_y, color = event, shape = event),
      size = 3.2,
      stroke = 1.0,
      show.legend = FALSE
    ) +
    geom_label(
      data = label_box(total_n * 0.05, max(variance_df$smooth_log1p_empirical_variance) * 0.96, geom_text),
      aes(x, y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1,
      size = 3.0, linewidth = 0.25, fill = alpha("white", 0.95)
    ) +
    geom_label(
      data = label_box(total_n * 0.72, max(variance_df$smooth_log1p_empirical_variance) * 0.96, geom_summary),
      aes(x, y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1,
      size = 3.0, linewidth = 0.25, fill = alpha("white", 0.95)
    ) +
    scale_color_manual(values = EVENT_COLORS, drop = FALSE) +
    scale_shape_manual(values = EVENT_SHAPES, drop = FALSE) +
    scale_linetype_manual(values = EVENT_LINE_TYPES, drop = FALSE) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": custom geometric transition"),
      subtitle = "The final interval is defined around the fixed leading-edge 5000 reference on the smoothed empirical variance curve",
      x = "EVS rank",
      y = "Smoothed log(1 + empirical variance)"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  # Legend strip
  p_legend_events <- {
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

  # Panel B: NB corroboration only
  nb_long <- feature_df %>%
    select(rank, log_nb2, nb_gap, log_alpha_mu) %>%
    pivot_longer(
      cols = c(log_nb2, nb_gap, log_alpha_mu),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c("log_nb2", "nb_gap", "log_alpha_mu"),
        labels = c(
          "NB2 = log(1 + variance - mu)",
          "NB2 - NB1 contrast",
          "alpha*mu = log(1 + alpha*mu)"
        )
      )
    )

  p2 <- ggplot() +
    annotate("rect", xmin = left_min, xmax = left_max, ymin = -Inf, ymax = Inf, fill = COL$left_fill, alpha = 0.72) +
    annotate("rect", xmin = right_min, xmax = right_max, ymin = -Inf, ymax = Inf, fill = COL$right_fill, alpha = 0.72) +
    annotate("rect", xmin = anchor, xmax = term, ymin = -Inf, ymax = Inf, fill = COL$interval_fill, alpha = 0.42) +
    geom_line(data = nb_long, aes(rank, value, color = metric), linewidth = 0.9) +
    base_lines +
    geom_label(
      data = label_box(total_n * 0.05, max(nb_long$value, na.rm = TRUE) * 0.95, nb_text),
      aes(x, y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1,
      size = 3.0, linewidth = 0.25, fill = alpha("white", 0.95)
    ) +
    geom_label(
      data = label_box(total_n * 0.72, max(nb_long$value, na.rm = TRUE) * 0.95, nb_summary_box),
      aes(x, y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1,
      size = 2.9, linewidth = 0.25, fill = alpha("white", 0.95)
    ) +
    scale_color_manual(
      values = c(
        "NB2 = log(1 + variance - mu)" = COL$nb2,
        "NB2 - NB1 contrast"           = COL$nb_gap,
        "alpha*mu = log(1 + alpha*mu)" = COL$alpha_mu,
        EVENT_COLORS
      ),
      drop = FALSE
    ) +
    scale_linetype_manual(values = EVENT_LINE_TYPES, drop = FALSE) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": NB2-related corroboration"),
      subtitle = "The right leading-edge block is compared directly against the matched left block",
      x = "EVS rank",
      y = "Corroborative NB2-related signal"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  # Legend strip
  p_legend_nb <- {
    df <- data.frame(
      x = seq_along(names(NB_TRACE_COLORS)),
      y = 1,
      metric = factor(names(NB_TRACE_COLORS), levels = names(NB_TRACE_COLORS))
    )

    ggplot(df, aes(x, y, color = metric)) +
      geom_point(size = 3) +
      geom_text(aes(label = metric), nudge_y = -0.16, size = 3, show.legend = FALSE) +
      scale_color_manual(values = NB_TRACE_COLORS, drop = FALSE) +
      xlim(0.5, length(names(NB_TRACE_COLORS)) + 0.5) +
      ylim(0.7, 1.15) +
      theme_void() +
      theme(legend.position = "none")
  }

  # Panel C: compact median summary
  summary_df <- data.frame(
    metric = factor(SUMMARY_METRIC_LEVELS, levels = rev(SUMMARY_METRIC_LEVELS)),
    left = c(
      region_summary$left_median_log_nb1,
      region_summary$left_median_log_nb2,
      region_summary$left_median_nb_gap,
      region_summary$left_median_log_alpha_mu
    ),
    right = c(
      region_summary$right_median_log_nb1,
      region_summary$right_median_log_nb2,
      region_summary$right_median_nb_gap,
      region_summary$right_median_log_alpha_mu
    ),
    stringsAsFactors = FALSE
  )

  summary_long <- summary_df %>%
    pivot_longer(cols = c(left, right), names_to = "region", values_to = "value") %>%
    mutate(region = factor(region, levels = c("left", "right"), labels = c("Matched LEFT", "RIGHT leading edge")))

  p3 <- ggplot(summary_long, aes(value, metric, color = region)) +
    geom_segment(
      data = summary_df,
      aes(x = left, xend = right, y = metric, yend = metric),
      inherit.aes = FALSE,
      color = "#7F7F7F",
      linewidth = 0.7
    ) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Matched LEFT" = COL$summary_left, "RIGHT leading edge" = COL$summary_right)) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": left-versus-right median summary"),
      subtitle = "Points farther to the right indicate stronger signal in that region",
      x = "Median value",
      y = NULL,
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "plain")
    )

  save_stacked_plot(
    list(
      p1,
      p_legend_events,
      p2,
      p_legend_nb,
      p3
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
        pre_evs_leading_edge_size = total_n - cutoff_anchor_rank + 1L,
        directional_call_nb2 = ifelse(right_left_log_nb2_diff > 0, "more_NB2_like_on_right", "not_more_NB2_like_on_right"),
        directional_call_gap = ifelse(right_left_nb_gap_diff > 0, "more_NB2_like_on_right", "not_more_NB2_like_on_right"),
        directional_call_alpha_mu = ifelse(right_left_log_alpha_mu_diff > 0, "more_NB2_like_on_right", "not_more_NB2_like_on_right")
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
      interval_info = interval_info,
      region_summary = region_summary,
      out_png = file.path(comp_dir, paste0(comparison_name, "_", arm_name, "_main_figure.png"))
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
