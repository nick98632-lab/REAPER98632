#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# MANUSCRIPT-LEVEL SCRIPT
# -----------------------------------------------------------------------------
# MAIN CLAIM
#
# A custom geometry-based cutoff near a fixed leading-edge reference defines a
# right-hand leading-edge block whose matched comparison against the immediate
# left block shows stronger NB2-like, overdispersion-consistent behavior.
#
# TWO ANALYTIC TRACKS
#
# 1. Main manuscript figure
#    - Ranking normalization: log-transformed CPM-style library-size normalization
#    - Variance / corroboration metrics: computed on raw counts
#
# 2. DESeq2 sensitivity figure
#    - Ranking normalization: DESeq2 normalized counts (or VST, optional)
#    - Variance / corroboration metrics: computed on DESeq2 normalized counts
#
# LITERATURE-SAFE FRAMING
#
# - The geometric cutoff rule is custom.
# - The NB interpretation is used as corroborative diagnostics:
#     NB2-like behavior is represented by higher-order variance-linked terms.
# - We do not present the custom cutoff rule as a standard published NB method.
#
# METHODS SUMMARY
#
# Features are ranked by absolute PC1 loading within each comparison arm.
#
# Main analysis:
#   Ranking uses log-transformed CPM-style library-size normalization.
#
# DESeq2 sensitivity analysis:
#   Ranking uses DESeq2 normalized counts, or optionally VST if requested and
#   available.
#
# Along the resulting rank axis, empirical variance is computed feature-wise.
# A smoothing spline is fit to the ranked log(1 + variance) trajectory.
# The second derivative of that spline is evaluated on a dense rank grid.
#
# Let R_fixed denote the fixed leading-edge reference rank corresponding to
# FIXED_LEADING_EDGE_SIZE features from the right edge of the ranked series.
#
# The final custom geometric interval is defined by the two nearest spline-
# based second-derivative zero-crossings that flank R_fixed:
#   - cutoff anchor = nearest zero-crossing immediately LEFT of R_fixed
#   - terminal start = nearest zero-crossing immediately RIGHT of R_fixed
#
# The full right leading-edge block is then defined as all ranks from the
# cutoff anchor through the right edge of the ranked series.
# The matched left comparator block is the equal-sized block immediately to the
# left of the cutoff anchor.
#
# NB2-RELATED CORROBORATION METRICS
#
# Let mu denote the empirical mean and variance denote the empirical variance
# across samples for each ranked feature.
#
# We use three corroborative quantities:
#
# 1. log_nb2 = log(1 + variance - mu)
#    This measures extra-Poisson variance, because a Poisson baseline satisfies
#    variance = mu.
#
# 2. nb_gap = log(1 + variance - mu) - log(1 + mu)
#    This compares higher-order excess-variance signal against the lower-order
#    mean-linked signal.
#
# 3. log_alpha_mu = log(1 + alpha*mu), where
#       alpha = max((variance - mu) / mu^2, 0)
#    Under the NB2 variance identity variance = mu + alpha*mu^2, alpha*mu is a
#    normalized overdispersion-linked signal.
#
# These are used as descriptive corroborative summaries, not as formal
# likelihood-ratio test statistics.
#
# METHODS-LEVEL VALIDATION
#
# For each comparison and arm, we explicitly verify that:
# - cutoff anchor < fixed leading-edge reference < terminal start
# - the matched left block and full right block have equal size
# - the reported left and right summary medians exactly match the underlying
#   sliced feature sets used in the figure and output tables
# =============================================================================

# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT   <- "/root/REAPER98632/exports/manuscript_custom_geometry_nb2_support_final"

FIXED_LEADING_EDGE_SIZE <- 5000L
VAR_SPLINE_SPAR         <- 0.60

PNG_WIDTH_IN  <- 14
PNG_HEIGHT_IN <- 11.5
PNG_DPI       <- 260

# DESeq2 sensitivity analysis
RUN_DESEQ2_SENSITIVITY <- TRUE
DESEQ2_RANK_METHOD <- "normalized_log1p"   # options: "normalized_log1p", "vst"

# Sample naming patterns
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
  variance_curve = "#117A65",

  nb2      = "#1B9E77",
  nb_gap   = "#CC1E8C",
  alpha_mu = "#386CB0",

  left_fill     = "#CBE3F8",
  right_fill    = "#DDF2D5",
  interval_fill = "#AFAFAF",

  cutoff_anchor  = "#000000",
  fixed_ref      = "#E69F00",
  terminal_start = "#D95F02",

  left_point  = "#5B8FD1",
  right_point = "#43A047"
)

EVENT_LEVELS <- c(
  "Cutoff anchor",
  "Fixed leading-edge reference",
  "Terminal start"
)

EVENT_COLORS <- c(
  "Cutoff anchor"             = COL$cutoff_anchor,
  "Fixed leading-edge reference" = COL$fixed_ref,
  "Terminal start"            = COL$terminal_start
)

EVENT_SHAPES <- c(
  "Cutoff anchor"             = 16,
  "Fixed leading-edge reference" = 18,
  "Terminal start"            = 1
)

EVENT_LTY <- c(
  "Cutoff anchor"             = "solid",
  "Fixed leading-edge reference" = "dashed",
  "Terminal start"            = "dotted"
)

TRACE_COLORS <- c(
  "NB2 = log(1 + variance - mu)" = COL$nb2,
  "NB2 - NB1 contrast"           = COL$nb_gap,
  "alpha*mu = log(1 + alpha*mu)" = COL$alpha_mu
)

SUMMARY_LEVELS <- c(
  "NB2 = log(1 + variance - mu)",
  "NB2 - NB1 contrast",
  "alpha*mu = log(1 + alpha*mu)"
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

compute_deseq2_matrices <- function(count_mat_arm, rank_method = "normalized_log1p") {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    return(NULL)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    return(NULL)
  }

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
    vst_obj <- tryCatch(
      DESeq2::vst(dds, blind = TRUE),
      error = function(e) NULL
    )
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
    dds = dds,
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

  spline_fit <- stats::smooth.spline(
    x = ranks,
    y = ranked_log_var,
    spar = spar
  )

  smooth_y <- as.numeric(stats::predict(spline_fit, x = ranks, deriv = 0)$y)
  smooth_d1 <- as.numeric(stats::predict(spline_fit, x = ranks, deriv = 1)$y)
  smooth_d2 <- as.numeric(stats::predict(spline_fit, x = ranks, deriv = 2)$y)

  dense_x <- seq(min(ranks), max(ranks), length.out = max(5000L, length(ranks) * 4L))
  dense_y  <- as.numeric(stats::predict(spline_fit, x = dense_x, deriv = 0)$y)
  dense_d2 <- as.numeric(stats::predict(spline_fit, x = dense_x, deriv = 2)$y)

  out <- data.frame(
    rank = ranks,
    empirical_variance = ranked_var,
    log1p_empirical_variance = ranked_log_var,
    smooth_log1p_empirical_variance = smooth_y,
    d1_spline = smooth_d1,
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
    return(data.frame(
      crossing_rank = numeric(0),
      crossing_type = character(0),
      stringsAsFactors = FALSE
    ))
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

select_custom_interval <- function(zero_df, reference_rank, total_n) {
  if (nrow(zero_df) == 0L) stop("No d2 sign-change crossings found.")

  left_candidates  <- zero_df$crossing_rank[zero_df$crossing_rank < reference_rank]
  right_candidates <- zero_df$crossing_rank[zero_df$crossing_rank > reference_rank]

  if (length(left_candidates) == 0L) stop("No left d2 crossing found.")
  if (length(right_candidates) == 0L) stop("No right d2 crossing found.")

  cutoff_anchor_rank <- as.integer(round(max(left_candidates)))
  terminal_start_rank <- as.integer(round(min(right_candidates)))
  reference_rank <- as.integer(round(reference_rank))

  cutoff_anchor_rank <- max(1L, cutoff_anchor_rank)
  terminal_start_rank <- min(total_n, terminal_start_rank)

  if (cutoff_anchor_rank >= reference_rank) {
    stop("Invalid interval: cutoff anchor must lie left of the fixed leading-edge reference.")
  }
  if (reference_rank >= terminal_start_rank) {
    stop("Invalid interval: terminal start must lie right of the fixed leading-edge reference.")
  }

  list(
    cutoff_anchor_rank = cutoff_anchor_rank,
    reference_rank = reference_rank,
    terminal_start_rank = terminal_start_rank,
    interval_min_rank = cutoff_anchor_rank,
    interval_max_rank = terminal_start_rank
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
  alpha_mu <- alpha_hat * mu

  data.frame(
    rank = seq_along(rank_order),
    feature_id = rownames(metric_mat_arm)[rank_order],
    mu = mu,
    empirical_variance = empirical_var,
    log_nb1 = log1p(mu),
    log_nb2 = log1p(nb2_var_minus_mu),
    log_alpha_mu = log1p(alpha_mu),
    nb_gap = log1p(nb2_var_minus_mu) - log1p(mu),
    stringsAsFactors = FALSE
  )
}

summarize_regions <- function(feature_df, cutoff_anchor_rank, total_n) {
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

    stringsAsFactors = FALSE
  )
}

validate_method_level <- function(feature_df, interval_info, region_summary, total_n) {
  anchor <- interval_info$cutoff_anchor_rank
  ref    <- interval_info$reference_rank
  term   <- interval_info$terminal_start_rank

  if (!(anchor < ref && ref < term)) {
    stop("Validation failed: expected cutoff anchor < fixed reference < terminal start.")
  }

  right_idx <- seq.int(anchor, total_n)
  right_n <- length(right_idx)

  left_end <- anchor - 1L
  left_start <- left_end - right_n + 1L
  if (left_start < 1L) {
    stop("Validation failed: matched LEFT block extends below rank 1.")
  }

  left_idx <- seq.int(left_start, left_end)

  if (length(left_idx) != length(right_idx)) {
    stop("Validation failed: matched LEFT and RIGHT blocks do not have equal size.")
  }

  left_df <- feature_df[left_idx, , drop = FALSE]
  right_df <- feature_df[right_idx, , drop = FALSE]

  expected <- list(
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
    right_left_log_alpha_mu_diff = median(right_df$log_alpha_mu, na.rm = TRUE) - median(left_df$log_alpha_mu, na.rm = TRUE)
  )

  for (nm in names(expected)) {
    chk <- all.equal(
      as.numeric(region_summary[[nm]][1]),
      as.numeric(expected[[nm]]),
      tolerance = 1e-10
    )
    if (!isTRUE(chk)) {
      stop("Validation failed for summary field: ", nm, " | ", chk)
    }
  }

  invisible(TRUE)
}

make_offset_event_df <- function(event_df, total_n) {
  offset_big <- max(8L, round(total_n * 0.0025))
  event_df %>%
    mutate(
      rank_plot = case_when(
        event == "Cutoff anchor" ~ rank - offset_big,
        event == "Fixed leading-edge reference" ~ rank,
        event == "Terminal start" ~ rank + offset_big,
        TRUE ~ rank
      )
    )
}

label_box <- function(x, y, txt) {
  data.frame(x = x, y = y, label = txt, stringsAsFactors = FALSE)
}

save_three_panel_plot <- function(plot_list, filename) {
  png(
    filename,
    width = PNG_WIDTH_IN,
    height = PNG_HEIGHT_IN,
    units = "in",
    res = PNG_DPI,
    bg = "white"
  )
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = length(plot_list),
    ncol = 1,
    heights = unit(c(1.2, 0.18, 1.2, 0.18, 0.95), "null")
  )))
  for (i in seq_along(plot_list)) {
    print(plot_list[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1))
  }
  dev.off()
}

# =============================================================================
# FIGURE BUILDER
# =============================================================================

build_main_figure <- function(comparison_name, arm_name, variance_df, feature_df,
                              interval_info, region_summary, out_file,
                              rank_label = "Main analysis: log-transformed CPM-style library-size normalization",
                              metric_label = "Metrics computed on raw counts") {

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
    event = factor(EVENT_LEVELS, levels = EVENT_LEVELS),
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
    event = factor(EVENT_LEVELS, levels = EVENT_LEVELS),
    xint = c(anchor, ref, term),
    stringsAsFactors = FALSE
  )

  geom_text <- paste(
    "Custom geometry panel",
    rank_label,
    metric_label,
    "Blue region = matched LEFT comparator",
    "Green region = full RIGHT leading-edge block",
    "Grey band = final geometric interval",
    "Cutoff anchor and terminal start are custom d2-zero landmarks around the fixed leading-edge reference",
    sep = "\n"
  )

  geom_summary <- paste0(
    "Cutoff anchor rank = ", anchor, "\n",
    "Reference rank (", FIXED_LEADING_EDGE_SIZE, " from right) = ", ref, "\n",
    "Terminal start rank = ", term, "\n",
    "Final interval = [", anchor, ", ", term, "]\n",
    "Pre-EVS remainder = ", anchor - 1L, "\n",
    "Pre-EVS leading edge = ", total_n - anchor + 1L
  )

  p1 <- ggplot(variance_df, aes(rank, smooth_log1p_empirical_variance)) +
    annotate("rect", xmin = left_min, xmax = left_max, ymin = -Inf, ymax = Inf,
             fill = COL$left_fill, alpha = 0.72) +
    annotate("rect", xmin = right_min, xmax = right_max, ymin = -Inf, ymax = Inf,
             fill = COL$right_fill, alpha = 0.72) +
    annotate("rect", xmin = anchor, xmax = term, ymin = -Inf, ymax = Inf,
             fill = COL$interval_fill, alpha = 0.42) +
    geom_line(color = COL$variance_curve, linewidth = 1.0) +
    geom_vline(
      data = vline_df,
      aes(xintercept = xint, color = event, linetype = event),
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    geom_point(
      data = event_plot_df,
      aes(rank_plot, variance_y, color = event, shape = event),
      size = 3.2,
      stroke = 1.0,
      show.legend = FALSE
    ) +
    geom_label(
      data = label_box(total_n * 0.03, max(variance_df$smooth_log1p_empirical_variance) * 0.97, geom_text),
      aes(x, y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      linewidth = 0.25,
      fill = grDevices::adjustcolor("white", alpha.f = 0.95)
    ) +
    geom_label(
      data = label_box(total_n * 0.76, max(variance_df$smooth_log1p_empirical_variance) * 0.82, geom_summary),
      aes(x, y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      linewidth = 0.25,
      fill = grDevices::adjustcolor("white", alpha.f = 0.95)
    ) +
    scale_color_manual(values = EVENT_COLORS, drop = FALSE) +
    scale_shape_manual(values = EVENT_SHAPES, drop = FALSE) +
    scale_linetype_manual(values = EVENT_LTY, drop = FALSE) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": custom geometric transition"),
      subtitle = "The final interval is defined around the fixed leading-edge reference on the smoothed empirical variance curve",
      x = "EVS rank",
      y = "Smoothed log(1 + empirical variance)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    )

  p_event_legend <- {
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

  p2 <- ggplot() +
    annotate("rect", xmin = left_min, xmax = left_max, ymin = -Inf, ymax = Inf,
             fill = COL$left_fill, alpha = 0.72) +
    annotate("rect", xmin = right_min, xmax = right_max, ymin = -Inf, ymax = Inf,
             fill = COL$right_fill, alpha = 0.72) +
    annotate("rect", xmin = anchor, xmax = term, ymin = -Inf, ymax = Inf,
             fill = COL$interval_fill, alpha = 0.42) +
    geom_line(
      data = nb_long,
      aes(rank, value, color = metric),
      linewidth = 0.9
    ) +
    geom_vline(
      data = vline_df,
      aes(xintercept = xint, color = event, linetype = event),
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    geom_label(
      data = label_box(total_n * 0.03, max(nb_long$value, na.rm = TRUE) * 0.97, nb_text),
      aes(x, y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      linewidth = 0.25,
      fill = grDevices::adjustcolor("white", alpha.f = 0.95)
    ) +
    geom_label(
      data = label_box(total_n * 0.76, max(nb_long$value, na.rm = TRUE) * 0.82, nb_summary_box),
      aes(x, y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      linewidth = 0.25,
      fill = grDevices::adjustcolor("white", alpha.f = 0.95)
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
    scale_linetype_manual(values = EVENT_LTY, drop = FALSE) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": NB2-related corroboration"),
      subtitle = "The right leading-edge block is compared directly against the matched left block",
      x = "EVS rank",
      y = "Corroborative NB2-related signal"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    )

  p_nb_legend <- {
    df <- data.frame(
      x = seq_along(names(TRACE_COLORS)),
      y = 1,
      metric = factor(names(TRACE_COLORS), levels = names(TRACE_COLORS))
    )

    ggplot(df, aes(x, y, color = metric)) +
      geom_point(size = 3) +
      geom_text(aes(label = metric), nudge_y = -0.16, size = 3, show.legend = FALSE) +
      scale_color_manual(values = TRACE_COLORS, drop = FALSE) +
      xlim(0.5, length(names(TRACE_COLORS)) + 0.5) +
      ylim(0.7, 1.15) +
      theme_void() +
      theme(legend.position = "none")
  }

  summary_df <- data.frame(
    metric = factor(SUMMARY_LEVELS, levels = rev(SUMMARY_LEVELS)),
    left = c(
      region_summary$left_median_log_nb2,
      region_summary$left_median_nb_gap,
      region_summary$left_median_log_alpha_mu
    ),
    right = c(
      region_summary$right_median_log_nb2,
      region_summary$right_median_nb_gap,
      region_summary$right_median_log_alpha_mu
    ),
    stringsAsFactors = FALSE
  )

  p3 <- ggplot(summary_df, aes(y = metric)) +
    geom_segment(
      aes(x = left, xend = right, yend = metric),
      color = "#7A7A7A",
      linewidth = 0.8
    ) +
    geom_point(aes(x = left, color = "Matched LEFT"), size = 3.2) +
    geom_point(aes(x = right, color = "RIGHT leading edge"), size = 3.2) +
    scale_color_manual(
      values = c(
        "Matched LEFT" = COL$left_point,
        "RIGHT leading edge" = COL$right_point
      )
    ) +
    labs(
      title = paste0(comparison_name, " ", arm_name, ": left-versus-right median summary"),
      subtitle = "Points farther right indicate stronger corroborative NB2-related signal in that region",
      x = "Median value",
      y = NULL,
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  save_three_panel_plot(
    list(
      p1,
      p_event_legend,
      p2,
      p_nb_legend,
      p3
    ),
    out_file
  )
}

# =============================================================================
# ONE ANALYSIS TRACK
# =============================================================================

run_one_track <- function(comparison_name,
                          arm_name,
                          count_mat_arm,
                          rank_matrix,
                          metric_matrix,
                          output_dir,
                          figure_suffix,
                          rank_label,
                          metric_label,
                          reference_label,
                          use_metric_name) {

  abs_loadings <- compute_abs_pc1_loadings(rank_matrix)
  rank_order <- order(abs_loadings, decreasing = FALSE)
  total_n <- length(rank_order)

  if (FIXED_LEADING_EDGE_SIZE >= total_n) {
    stop("FIXED_LEADING_EDGE_SIZE must be < total_n")
  }

  # The fixed leading-edge reference is defined as the rank that leaves exactly
  # FIXED_LEADING_EDGE_SIZE features on the right side, including that reference
  # rank itself.
  #
  # The final custom interval is then defined by the nearest spline-based
  # second-derivative zero-crossing immediately left of this reference and the
  # nearest spline-based second-derivative zero-crossing immediately right of it.
  reference_rank <- total_n - FIXED_LEADING_EDGE_SIZE + 1L

  variance_df <- compute_ranked_variance_curve(
    metric_mat_arm = metric_matrix,
    rank_order = rank_order,
    spar = VAR_SPLINE_SPAR
  )

  dense_curve_df <- attr(variance_df, "dense_curve_df")
  zero_df <- find_d2_zero_crossings(dense_curve_df)

  write.csv(
    zero_df,
    file = file.path(output_dir, paste0(comparison_name, "_", arm_name, "_", figure_suffix, "_zero_crossings_all.csv")),
    row.names = FALSE
  )

  write.csv(
    dense_curve_df,
    file = file.path(output_dir, paste0(comparison_name, "_", arm_name, "_", figure_suffix, "_dense_spline_curve.csv")),
    row.names = FALSE
  )

  interval_info <- select_custom_interval(zero_df, reference_rank, total_n)

  selected_df <- data.frame(
    comparison_name = comparison_name,
    arm = arm_name,
    analysis_track = figure_suffix,
    ranking_method = rank_label,
    metric_matrix_used = use_metric_name,
    cutoff_anchor_rank = interval_info$cutoff_anchor_rank,
    fixed_leading_edge_reference_rank = interval_info$reference_rank,
    terminal_start_rank = interval_info$terminal_start_rank,
    interval_min_rank = interval_info$interval_min_rank,
    interval_max_rank = interval_info$interval_max_rank,
    stringsAsFactors = FALSE
  )

  write.csv(
    selected_df,
    file = file.path(output_dir, paste0(comparison_name, "_", arm_name, "_", figure_suffix, "_selected_two_zero_crossings.csv")),
    row.names = FALSE
  )

  feature_df <- compute_ranked_feature_metrics(metric_matrix, rank_order)
  feature_df$abs_pc1_loading <- abs_loadings[rank_order]

  rank_series_df <- feature_df %>%
    left_join(variance_df, by = "rank")

  write.csv(
    rank_series_df,
    file = file.path(output_dir, paste0(comparison_name, "_", arm_name, "_", figure_suffix, "_rank_series.csv")),
    row.names = FALSE
  )

  region_summary <- summarize_regions(feature_df, interval_info$cutoff_anchor_rank, total_n)

  validate_method_level(
    feature_df = feature_df,
    interval_info = interval_info,
    region_summary = region_summary,
    total_n = total_n
  )

  validation_df <- data.frame(
    comparison_name = comparison_name,
    arm = arm_name,
    analysis_track = figure_suffix,
    validation_status = "PASS",
    anchor_lt_reference = interval_info$cutoff_anchor_rank < interval_info$reference_rank,
    reference_lt_terminal = interval_info$reference_rank < interval_info$terminal_start_rank,
    matched_left_n = region_summary$left_n,
    matched_right_n = region_summary$right_n,
    stringsAsFactors = FALSE
  )

  write.csv(
    validation_df,
    file = file.path(output_dir, paste0(comparison_name, "_", arm_name, "_", figure_suffix, "_validation_report.csv")),
    row.names = FALSE
  )

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
    file = file.path(output_dir, paste0(comparison_name, "_", arm_name, "_", figure_suffix, "_cutoffs_summary.csv")),
    row.names = FALSE
  )

  build_main_figure(
    comparison_name = comparison_name,
    arm_name = arm_name,
    variance_df = variance_df,
    feature_df = feature_df,
    interval_info = interval_info,
    region_summary = region_summary,
    out_file = file.path(output_dir, paste0(comparison_name, "_", arm_name, "_", figure_suffix, "_main_figure.png")),
    rank_label = rank_label,
    metric_label = metric_label
  )

  list(
    selected_df = selected_df,
    cutoff_summary = cutoff_summary
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

    # -------------------------------------------------------------------------
    # MAIN ANALYSIS TRACK
    # Ranking: log-transformed CPM-style library-size normalization
    # Metrics: raw counts
    # -------------------------------------------------------------------------
    main_rank_matrix <- normalize_cpm_log1p(count_mat_arm)

    main_res <- run_one_track(
      comparison_name = comparison_name,
      arm_name = arm_name,
      count_mat_arm = count_mat_arm,
      rank_matrix = main_rank_matrix,
      metric_matrix = count_mat_arm,
      output_dir = comp_dir,
      figure_suffix = "main_cpm",
      rank_label = "Main analysis: log-transformed CPM-style library-size normalization",
      metric_label = "Metrics computed on raw counts",
      reference_label = "fixed leading-edge reference",
      use_metric_name = "raw_counts"
    )

    overall_rows[[length(overall_rows) + 1L]] <- main_res$cutoff_summary

    message(
      "[MAIN] ", tools::toTitleCase(arm_name),
      " cutoff anchor rank: ", main_res$selected_df$cutoff_anchor_rank,
      " | reference rank: ", main_res$selected_df$fixed_leading_edge_reference_rank,
      " | terminal start: ", main_res$selected_df$terminal_start_rank
    )

    # -------------------------------------------------------------------------
    # DESeq2 SENSITIVITY TRACK
    # Ranking: DESeq2 normalized counts or VST
    # Metrics: DESeq2 normalized counts
    # -------------------------------------------------------------------------
    if (RUN_DESEQ2_SENSITIVITY) {
      deseq2_obj <- compute_deseq2_matrices(
        count_mat_arm = count_mat_arm,
        rank_method = DESEQ2_RANK_METHOD
      )

      if (is.null(deseq2_obj)) {
        message("[DESEQ2] Skipped for ", comparison_name, " ", arm_name, " because DESeq2 is not available.")
      } else {
        deseq2_res <- run_one_track(
          comparison_name = comparison_name,
          arm_name = arm_name,
          count_mat_arm = count_mat_arm,
          rank_matrix = deseq2_obj$ranking_matrix,
          metric_matrix = deseq2_obj$normalized_counts,
          output_dir = comp_dir,
          figure_suffix = "deseq2_sensitivity",
          rank_label = paste0(
            "DESeq2 sensitivity: ",
            if (DESEQ2_RANK_METHOD == "vst") "VST-transformed DESeq2 normalized counts"
            else "log-transformed DESeq2 normalized counts"
          ),
          metric_label = "Metrics computed on DESeq2 normalized counts",
          reference_label = "fixed leading-edge reference",
          use_metric_name = "deseq2_normalized_counts"
        )

        write.csv(
          data.frame(
            sample = names(deseq2_obj$size_factors),
            size_factor = as.numeric(deseq2_obj$size_factors),
            stringsAsFactors = FALSE
          ),
          file = file.path(comp_dir, paste0(comparison_name, "_", arm_name, "_deseq2_size_factors.csv")),
          row.names = FALSE
        )

        overall_rows[[length(overall_rows) + 1L]] <- deseq2_res$cutoff_summary

        message(
          "[DESEQ2] ", tools::toTitleCase(arm_name),
          " cutoff anchor rank: ", deseq2_res$selected_df$cutoff_anchor_rank,
          " | reference rank: ", deseq2_res$selected_df$fixed_leading_edge_reference_rank,
          " | terminal start: ", deseq2_res$selected_df$terminal_start_rank
        )
      }
    }
  }
}

overall_summary <- bind_rows(overall_rows)

write.csv(
  overall_summary,
  file = file.path(OUT_ROOT, "overall_cutoff_summary.csv"),
  row.names = FALSE
)

message("Done. Outputs written to: ", OUT_ROOT)
