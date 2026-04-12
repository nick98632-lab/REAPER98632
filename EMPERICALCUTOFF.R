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
# =============================================================================
#
# PURPOSE
#
# This script produces:
#
# 1. Main manuscript figures
#    - Rank: log-transformed CPM-style library-size normalization
#    - Metrics: raw counts
#
# 2. DESeq2 supplementary figures
#    - Rank: DESeq2-normalized counts, optionally VST
#    - Metrics: DESeq2-normalized counts
#
# SCIENTIFIC LOGIC
#
# A. Geometry
#    Features are ranked within each arm by absolute PC1 loading.
#    Along that rank axis, empirical variance is computed feature-wise and
#    transformed as log(1 + variance). A smoothing spline is fit to the ranked
#    variance trajectory. This smooth curve is used because a stable continuous
#    second derivative is more defensible than differentiating a noisy jagged
#    empirical series directly.
#
#    A fixed leading-edge reference is defined as the rank leaving exactly
#    FIXED_LEADING_EDGE_SIZE features on the right side, including that rank.
#
#    The final custom interval is defined by the two nearest spline-based
#    second-derivative zero-crossings flanking that fixed reference:
#      - Anchor   = nearest d2 zero immediately LEFT of the fixed reference
#      - Terminal = nearest d2 zero immediately RIGHT of the fixed reference
#
#    This is a custom geometric rule. It is not presented as a standard
#    published cutoff procedure.
#
# B. Corroboration
#    RIGHT = all ranks from Anchor through the right edge of the ranked series
#    LEFT  = equal-sized matched block immediately left of Anchor
#
#    The hypothesis is that RIGHT shows stronger NB2-like,
#    overdispersion-consistent behavior than LEFT.
#
# NB2-RELATED QUANTITIES
#
# Let mu denote empirical mean and variance denote empirical variance.
#
# 1. NB2 = log(1 + variance - mu)
#    Extra-Poisson variance signal.
#
# 2. NB2-NB1 = log(1 + variance - mu) - log(1 + mu)
#    Higher-order excess-variance signal relative to lower-order mean signal.
#
# 3. alpha*mu = log(1 + alpha*mu), where
#       alpha = max((variance - mu) / mu^2, 0)
#    Under variance = mu + alpha*mu^2, this is a normalized NB2-linked signal.
#
# These are descriptive corroborative diagnostics. They are not formal
# likelihood-ratio statistics.
#
# METHODS-LEVEL VALIDATION
#
# For every comparison arm and every analysis track, the script verifies that:
# - Anchor < Ref < Terminal
# - LEFT and RIGHT have equal size
# - reported summary medians exactly match the sliced plotted regions
#
# FIGURE RULES
#
# - 3 panels only
# - no separate legend-strip panels
# - no annotation boxes in the cutoff zone
# - short labels only inside panels
# - file names start with Figure_ or Table_
# =============================================================================

# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT   <- "/root/REAPER98632/exports/manuscript_final_clean"

FIXED_LEADING_EDGE_SIZE <- 5000L
VAR_SPLINE_SPAR         <- 0.60

PNG_WIDTH_IN  <- 14
PNG_HEIGHT_IN <- 10.8
PNG_DPI       <- 260

RUN_DESEQ2_SUPPLEMENT <- TRUE
DESEQ2_RANK_METHOD <- "normalized_log1p"   # "normalized_log1p" or "vst"

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
  var_curve = "#117A65",
  nb2       = "#1B9E77",
  nb_gap    = "#CC1E8C",
  alpha_mu  = "#386CB0",

  left_fill    = "#CBE3F8",
  right_fill   = "#DDF2D5",
  interval_fill = "#9E9E9E",

  anchor   = "#000000",
  ref      = "#E69F00",
  terminal = "#D95F02",

  left_pt  = "#5B8FD1",
  right_pt = "#43A047"
)

EVENT_LEVELS <- c("Anchor", "Ref", "Terminal")
EVENT_COLORS <- c("Anchor" = COL$anchor, "Ref" = COL$ref, "Terminal" = COL$terminal)
EVENT_SHAPES <- c("Anchor" = 16, "Ref" = 18, "Terminal" = 1)
EVENT_LTY    <- c("Anchor" = "solid", "Ref" = "dashed", "Terminal" = "dotted")

TRACE_LEVELS <- c("NB2", "NB2-NB1", "alpha*mu")
TRACE_COLORS <- c("NB2" = COL$nb2, "NB2-NB1" = COL$nb_gap, "alpha*mu" = COL$alpha_mu)

REGION_LEVELS <- c("LEFT", "RIGHT")
REGION_COLORS <- c("LEFT" = COL$left_pt, "RIGHT" = COL$right_pt)

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

compute_deseq2_matrices <- function(count_mat_arm, rank_method = "normalized_log1p") {
  if (!requireNamespace("DESeq2", quietly = TRUE)) return(NULL)
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) return(NULL)

  col_data <- data.frame(
    row.names = colnames(count_mat_arm),
    intercept = factor(rep("one", ncol(count_mat_arm)))
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_mat_arm),
    colData = col_data,
    design = ~ 1
  )

  dds <- DESeq2::estimateSizeFactors(dds)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  vst_mat <- NULL
  if (rank_method == "vst") {
    vst_obj <- tryCatch(DESeq2::vst(dds, blind = TRUE), error = function(e) NULL)
    if (!is.null(vst_obj)) {
      vst_mat <- SummarizedExperiment::assay(vst_obj)
    }
  }

  ranking_matrix <- switch(
    rank_method,
    normalized_log1p = log1p(norm_counts),
    vst = if (!is.null(vst_mat)) vst_mat else log1p(norm_counts),
    log1p(norm_counts)
  )

  list(
    normalized_counts = norm_counts,
    ranking_matrix = ranking_matrix,
    size_factors = DESeq2::sizeFactors(dds)
  )
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

  dense_x <- seq(min(ranks), max(ranks), length.out = max(5000L, length(ranks) * 4L))
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

make_display_event_df <- function(event_df, total_n) {
  offset_big <- max(12L, round(total_n * 0.006))
  event_df %>%
    mutate(
      rank_display = case_when(
        event == "Anchor" ~ rank - offset_big,
        event == "Ref" ~ rank,
        event == "Terminal" ~ rank + offset_big,
        TRUE ~ rank
      )
    )
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

compute_text_positions <- function(total_n, anchor, ref, terminal) {
  pad <- max(50L, round(total_n * 0.04))
  safe_right_limit <- max(200L, anchor - pad)

  left_x <- max(5L, round(total_n * 0.035))
  stat_x <- min(max(150L, round(total_n * 0.26)), safe_right_limit)

  if (stat_x <= left_x + 50L) {
    stat_x <- left_x + 60L
  }

  list(left_x = left_x, stat_x = stat_x)
}

# =============================================================================
# FIGURE BUILDER
# =============================================================================

build_main_figure <- function(comparison_name,
                              arm_name,
                              variance_df,
                              feature_df,
                              interval_info,
                              region_summary,
                              out_file,
                              fig_tag,
                              rank_tag,
                              metric_tag) {

  total_n <- nrow(feature_df)
  anchor <- interval_info$anchor
  ref <- interval_info$ref
  terminal <- interval_info$terminal

  left_n <- region_summary$left_n
  left_min <- anchor - left_n
  left_max <- anchor - 1L
  right_min <- anchor
  right_max <- total_n

  event_df <- data.frame(
    event = factor(EVENT_LEVELS, levels = EVENT_LEVELS),
    rank = c(anchor, ref, terminal),
    y = c(
      variance_df$smooth_log1p_empirical_variance[anchor],
      variance_df$smooth_log1p_empirical_variance[ref],
      variance_df$smooth_log1p_empirical_variance[terminal]
    ),
    stringsAsFactors = FALSE
  )
  event_display_df <- make_display_event_df(event_df, total_n)

  vline_df <- data.frame(
    event = factor(EVENT_LEVELS, levels = EVENT_LEVELS),
    xint = c(anchor, ref, terminal),
    stringsAsFactors = FALSE
  )

  nb_long <- feature_df %>%
    select(rank, NB2, NB2_NB1, alpha_mu) %>%
    pivot_longer(cols = c(NB2, NB2_NB1, alpha_mu), names_to = "metric", values_to = "value") %>%
    mutate(
      metric = factor(metric, levels = c("NB2", "NB2_NB1", "alpha_mu"),
                      labels = c("NB2", "NB2-NB1", "alpha*mu"))
    )

  pos <- compute_text_positions(total_n, anchor, ref, terminal)

  top_y <- max(variance_df$smooth_log1p_empirical_variance, na.rm = TRUE)
  mid_y <- max(nb_long$value, na.rm = TRUE)

  box1 <- paste(fig_tag, rank_tag, metric_tag, "Blue = LEFT", "Green = RIGHT", "Grey = interval", sep = "\n")
  box2 <- paste0(
    "Anchor = ", anchor, "\n",
    "Ref = ", ref, "\n",
    "Terminal = ", terminal, "\n",
    "Interval = [", anchor, ", ", terminal, "]\n",
    "LEFT n = ", region_summary$left_n, "\n",
    "RIGHT n = ", region_summary$right_n
  )

  box3 <- paste("RIGHT = Anchor to end", "LEFT = matched block", "Higher RIGHT = more NB2-like", sep = "\n")
  box4 <- paste0(
    "LEFT NB2 = ", round(region_summary$left_NB2, 3), "\n",
    "RIGHT NB2 = ", round(region_summary$right_NB2, 3), "\n",
    "RIGHT-LEFT NB2 = ", round(region_summary$diff_NB2, 3), "\n",
    "LEFT NB2-NB1 = ", round(region_summary$left_gap, 3), "\n",
    "RIGHT NB2-NB1 = ", round(region_summary$right_gap, 3), "\n",
    "RIGHT-LEFT NB2-NB1 = ", round(region_summary$diff_gap, 3), "\n",
    "LEFT alpha*mu = ", round(region_summary$left_alpha, 3), "\n",
    "RIGHT alpha*mu = ", round(region_summary$right_alpha, 3), "\n",
    "RIGHT-LEFT alpha*mu = ", round(region_summary$diff_alpha, 3)
  )

  p1 <- ggplot(variance_df, aes(rank, smooth_log1p_empirical_variance)) +
    annotate("rect", xmin = left_min, xmax = left_max, ymin = -Inf, ymax = Inf, fill = COL$left_fill, alpha = 0.70) +
    annotate("rect", xmin = right_min, xmax = right_max, ymin = -Inf, ymax = Inf, fill = COL$right_fill, alpha = 0.70) +
    annotate("rect", xmin = anchor, xmax = terminal, ymin = -Inf, ymax = Inf, fill = COL$interval_fill, alpha = 0.18) +
    geom_line(color = COL$var_curve, linewidth = 1.0) +
    geom_vline(
      data = vline_df,
      aes(xintercept = xint, color = event, linetype = event),
      linewidth = 0.9
    ) +
    geom_point(
      data = event_display_df,
      aes(rank_display, y, color = event, shape = event),
      size = 3.4,
      stroke = 1.0
    ) +
    annotate(
      "label",
      x = pos$left_x,
      y = top_y * 0.96,
      label = box1,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.25,
      fill = grDevices::adjustcolor("white", alpha.f = 0.96)
    ) +
    annotate(
      "label",
      x = pos$stat_x,
      y = top_y * 0.70,
      label = box2,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.25,
      fill = grDevices::adjustcolor("white", alpha.f = 0.96)
    ) +
    scale_color_manual(
      values = EVENT_COLORS,
      breaks = EVENT_LEVELS,
      guide = guide_legend(
        override.aes = list(
          shape = unname(EVENT_SHAPES),
          linetype = unname(EVENT_LTY),
          linewidth = 1.0,
          size = 3.4
        )
      )
    ) +
    scale_shape_manual(values = EVENT_SHAPES, guide = "none") +
    scale_linetype_manual(values = EVENT_LTY, guide = "none") +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": geometry"),
      subtitle = "Anchor, Ref, and Terminal define the custom interval",
      x = "Rank",
      y = "Smoothed log(1 + variance)",
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  p2 <- ggplot() +
    annotate("rect", xmin = left_min, xmax = left_max, ymin = -Inf, ymax = Inf, fill = COL$left_fill, alpha = 0.70) +
    annotate("rect", xmin = right_min, xmax = right_max, ymin = -Inf, ymax = Inf, fill = COL$right_fill, alpha = 0.70) +
    annotate("rect", xmin = anchor, xmax = terminal, ymin = -Inf, ymax = Inf, fill = COL$interval_fill, alpha = 0.18) +
    geom_vline(
      data = vline_df,
      aes(xintercept = xint),
      color = "grey35",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_line(
      data = nb_long,
      aes(rank, value, color = metric),
      linewidth = 0.95
    ) +
    annotate(
      "label",
      x = pos$left_x,
      y = mid_y * 0.96,
      label = box3,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.25,
      fill = grDevices::adjustcolor("white", alpha.f = 0.96)
    ) +
    annotate(
      "label",
      x = pos$stat_x,
      y = mid_y * 0.70,
      label = box4,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.25,
      fill = grDevices::adjustcolor("white", alpha.f = 0.96)
    ) +
    scale_color_manual(values = TRACE_COLORS, breaks = TRACE_LEVELS) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": corroboration"),
      subtitle = "RIGHT is compared directly against the matched LEFT block",
      x = "Rank",
      y = "NB2-related signal",
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  summary_df <- data.frame(
    metric = factor(c("NB2", "NB2-NB1", "alpha*mu"), levels = rev(c("NB2", "NB2-NB1", "alpha*mu"))),
    LEFT = c(region_summary$left_NB2, region_summary$left_gap, region_summary$left_alpha),
    RIGHT = c(region_summary$right_NB2, region_summary$right_gap, region_summary$right_alpha),
    stringsAsFactors = FALSE
  )

  p3 <- ggplot(summary_df, aes(y = metric)) +
    geom_segment(aes(x = LEFT, xend = RIGHT, yend = metric), color = "#7A7A7A", linewidth = 0.8) +
    geom_point(aes(x = LEFT, color = "LEFT"), size = 3.4) +
    geom_point(aes(x = RIGHT, color = "RIGHT"), size = 3.4) +
    scale_color_manual(values = REGION_COLORS, breaks = REGION_LEVELS) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": summary"),
      subtitle = "Points farther right indicate stronger NB2-related corroboration",
      x = "Median",
      y = NULL,
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  save_three_panel_plot(list(p1, p2, p3), out_file)
}

# =============================================================================
# ONE ANALYSIS TRACK
# =============================================================================

run_one_track <- function(comparison_name,
                          arm_name,
                          rank_matrix,
                          metric_matrix,
                          output_dir,
                          track,
                          fig_tag,
                          rank_tag,
                          metric_tag,
                          metric_name) {

  abs_loadings <- compute_abs_pc1_loadings(rank_matrix)
  rank_order <- order(abs_loadings, decreasing = FALSE)
  total_n <- length(rank_order)

  if (FIXED_LEADING_EDGE_SIZE >= total_n) {
    stop("FIXED_LEADING_EDGE_SIZE must be < total_n")
  }

  reference_rank <- total_n - FIXED_LEADING_EDGE_SIZE + 1L

  variance_df <- compute_ranked_variance_curve(
    metric_mat_arm = metric_matrix,
    rank_order = rank_order,
    spar = VAR_SPLINE_SPAR
  )

  dense_curve_df <- attr(variance_df, "dense_curve_df")
  zero_df <- find_d2_zero_crossings(dense_curve_df)
  interval_info <- select_custom_interval(zero_df, reference_rank, total_n)

  feature_df <- compute_ranked_feature_metrics(metric_matrix, rank_order)
  feature_df$abs_pc1_loading <- abs_loadings[rank_order]

  region_summary <- summarize_regions(feature_df, interval_info$anchor, total_n)

  validate_method_level(
    feature_df = feature_df,
    interval_info = interval_info,
    region_summary = region_summary,
    total_n = total_n
  )

  selected_df <- data.frame(
    comp = comparison_name,
    arm = arm_name,
    track = track,
    rank_method = rank_tag,
    metric_matrix = metric_name,
    Anchor = interval_info$anchor,
    Ref = interval_info$ref,
    Terminal = interval_info$terminal,
    IntMin = interval_info$interval_min,
    IntMax = interval_info$interval_max,
    stringsAsFactors = FALSE
  )

  cutoff_summary <- bind_cols(selected_df, region_summary) %>%
    mutate(
      LeftSize = region_summary$left_n,
      RightSize = region_summary$right_n,
      Call_NB2 = ifelse(diff_NB2 > 0, "RIGHT", "NOT_RIGHT"),
      Call_Gap = ifelse(diff_gap > 0, "RIGHT", "NOT_RIGHT"),
      Call_Alpha = ifelse(diff_alpha > 0, "RIGHT", "NOT_RIGHT")
    )

  zero_path  <- file.path(output_dir, paste0("Table_Zero_", comparison_name, "_", arm_name, "_", track, ".csv"))
  valid_path <- file.path(output_dir, paste0("Table_Valid_", comparison_name, "_", arm_name, "_", track, ".csv"))
  rank_path  <- file.path(output_dir, paste0("Table_Rank_", comparison_name, "_", arm_name, "_", track, ".csv"))
  cut_path   <- file.path(output_dir, paste0("Table_Cutoff_", comparison_name, "_", arm_name, "_", track, ".csv"))
  fig_path   <- file.path(output_dir, paste0("Figure_", track, "_", comparison_name, "_", arm_name, ".png"))

  write.csv(zero_df, zero_path, row.names = FALSE)

  write.csv(
    data.frame(
      comp = comparison_name,
      arm = arm_name,
      track = track,
      Status = "PASS",
      Anchor_lt_Ref = interval_info$anchor < interval_info$ref,
      Ref_lt_Terminal = interval_info$ref < interval_info$terminal,
      LeftN = region_summary$left_n,
      RightN = region_summary$right_n,
      stringsAsFactors = FALSE
    ),
    valid_path,
    row.names = FALSE
  )

  write.csv(
    feature_df %>% left_join(variance_df, by = "rank"),
    rank_path,
    row.names = FALSE
  )

  write.csv(cutoff_summary, cut_path, row.names = FALSE)

  build_main_figure(
    comparison_name = comparison_name,
    arm_name = arm_name,
    variance_df = variance_df,
    feature_df = feature_df,
    interval_info = interval_info,
    region_summary = region_summary,
    out_file = fig_path,
    fig_tag = fig_tag,
    rank_tag = rank_tag,
    metric_tag = metric_tag
  )

  list(
    summary = cutoff_summary,
    selected = selected_df
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
  comp_dir <- file.path(OUT_ROOT, comparison_name)
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

  pats <- COMPARISONS[[comparison_name]]

  for (arm_name in c("control", "treatment")) {
    sample_idx <- grep(pats[[arm_name]], colnames(count_mat))
    if (length(sample_idx) < 2L) {
      stop("Not enough samples for ", comparison_name, " ", arm_name)
    }

    count_mat_arm <- count_mat[, sample_idx, drop = FALSE]

    main_rank_matrix <- normalize_cpm_log1p(count_mat_arm)

    main_res <- run_one_track(
      comparison_name = comparison_name,
      arm_name = arm_name,
      rank_matrix = main_rank_matrix,
      metric_matrix = count_mat_arm,
      output_dir = comp_dir,
      track = "Main",
      fig_tag = "Main",
      rank_tag = "Rank: CPM log1p",
      metric_tag = "Metrics: raw counts",
      metric_name = "raw_counts"
    )

    overall_rows[[length(overall_rows) + 1L]] <- main_res$summary

    message(
      "[Main] ", comparison_name, " ", arm_name,
      " | Anchor=", main_res$selected$Anchor,
      " Ref=", main_res$selected$Ref,
      " Terminal=", main_res$selected$Terminal
    )

    if (RUN_DESEQ2_SUPPLEMENT) {
      deseq2_obj <- compute_deseq2_matrices(count_mat_arm, rank_method = DESEQ2_RANK_METHOD)

      if (is.null(deseq2_obj)) {
        message("[DESeq2] skipped: package not available for ", comparison_name, " ", arm_name)
      } else {
        deseq2_rank_tag <- if (DESEQ2_RANK_METHOD == "vst") {
          "Rank: DESeq2 VST"
        } else {
          "Rank: DESeq2 log1p"
        }

        deseq2_res <- run_one_track(
          comparison_name = comparison_name,
          arm_name = arm_name,
          rank_matrix = deseq2_obj$ranking_matrix,
          metric_matrix = deseq2_obj$normalized_counts,
          output_dir = comp_dir,
          track = "DESeq2",
          fig_tag = "DESeq2 supplement",
          rank_tag = deseq2_rank_tag,
          metric_tag = "Metrics: DESeq2 normalized",
          metric_name = "deseq2_normalized_counts"
        )

        write.csv(
          data.frame(
            sample = names(deseq2_obj$size_factors),
            size_factor = as.numeric(deseq2_obj$size_factors),
            stringsAsFactors = FALSE
          ),
          file.path(comp_dir, paste0("Table_SizeFactor_", comparison_name, "_", arm_name, ".csv")),
          row.names = FALSE
        )

        overall_rows[[length(overall_rows) + 1L]] <- deseq2_res$summary

        message(
          "[DESeq2] ", comparison_name, " ", arm_name,
          " | Anchor=", deseq2_res$selected$Anchor,
          " Ref=", deseq2_res$selected$Ref,
          " Terminal=", deseq2_res$selected$Terminal
        )
      }
    }
  }
}

overall_summary <- bind_rows(overall_rows)

write.csv(
  overall_summary,
  file.path(OUT_ROOT, "Table_Overall_Cutoff.csv"),
  row.names = FALSE
)

message("Done. Outputs written to: ", OUT_ROOT)
