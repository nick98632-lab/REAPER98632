# =============================================================================
# FINAL MANUSCRIPT SCRIPT
# PRE-EVS LEADING-EDGE CUT-OFF SELECTION
# FULL REWRITE - RUNNABLE VERSION
#
# WHAT THIS SCRIPT DOES
# ---------------------
# 1. Ranks features by absolute PC1 loading, with the leading edge on the RIGHT.
# 2. Fits a smooth empirical variance curve on the EVS rank axis.
# 3. Uses second-derivative zero-crossings to define a geometric interval:
#      - cutoff anchor  = nearest d2 zero immediately LEFT of the fixed
#                         leading-edge-5000 reference mark
#      - terminal start = nearest d2 zero immediately RIGHT of the fixed
#                         leading-edge-5000 reference mark
# 4. Uses NB-style quantities only as corroborative evidence:
#      - RIGHT region = all genes from cutoff anchor to the end of the ranked series
#      - LEFT region  = matched left block of equal size immediately left of cutoff anchor
#      - NB1 = mu
#      - NB2 = variance - mu
#      - alpha*mu = (variance - mu) / mu when positive
# 5. Produces manuscript-style figures and CSV outputs for each comparison arm.
#
# IMPORTANT DEFINITIONS
# ---------------------
# Leading edge:
#   The high-loading region on the RIGHT side of the EVS-ranked series.
#
# Fixed leading-edge-5000 reference:
#   The rank position located 5000 genes from the RIGHT end of the ranked series.
#
# Geometric interval:
#   [cutoff anchor, terminal start]
#
# NB corroboration:
#   The geometric split is defined by the variance curve and second derivative.
#   NB quantities are used only to support why the RIGHT side is more NB2-like.
#
# FIGURE SYMBOLS
# --------------
# Filled circle = cutoff anchor
# Diamond       = fixed leading-edge-5000 mark
# Open circle   = terminal start
# Dashed line   = cutoff anchor
# Dotdash line  = fixed leading-edge-5000 mark
# Dotted line   = terminal start
# Grey band     = final geometric interval
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(gridExtra)
  library(grid)
})

options(stringsAsFactors = FALSE)
options(warn = 1)

# =============================================================================
# USER INPUTS
# =============================================================================

repo_dir <- getwd()

count_file_hint <- file.path(
  repo_dir,
  "data",
  "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
)

output_root <- file.path(
  repo_dir,
  "exports",
  "variance_derivative_nb_range_manuscript_final"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

fixed_leading_edge_k <- 5000L
variance_spline_spar <- 0.60
zero_tol <- 1e-12

comparison_table <- data.frame(
  comparison_name  = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  treatment_prefix = c("R0", "R2", "R4", "R8"),
  control_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14")
)

meta_all <- data.frame(
  id = c(
    "R0_1","R0_2","R0_3","R0_4","R0_5","ZT6_1","ZT6_2","ZT6_3","ZT6_4","ZT6_5",
    "R2_1","R2_2","R2_3","R2_4","R2_5","ZT8_1","ZT8_2","ZT8_3","ZT8_4","ZT8_5",
    "R4_1","R4_2","R4_3","R4_4","R4_5","ZT10_1","ZT10_2","ZT10_3","ZT10_4","ZT10_5",
    "R8_1","R8_2","R8_3","R8_4","R8_5","ZT14_1","ZT14_2","ZT14_3","ZT14_4","ZT14_5"
  ),
  condition = c(
    rep("trt", 5), rep("untrt", 5),
    rep("trt", 5), rep("untrt", 5),
    rep("trt", 5), rep("untrt", 5),
    rep("trt", 5), rep("untrt", 5)
  )
)
rownames(meta_all) <- meta_all$id

# =============================================================================
# COLOR SYSTEM
# =============================================================================

COL <- list(
  abs_loading   = "#111111",
  variance_fit  = "#111111",
  d1            = "#00A6D6",
  d2            = "#D1495B",
  cutoff_anchor = "#111111",
  fixed_5000    = "#8C510A",
  terminal      = "#E66101",
  interval_fill = "#BDBDBD",
  nb_support    = "#111111",
  nb1           = "#C8A100",
  nb2           = "#33A02C",
  alpha_mu      = "#1F78B4",
  nb_gap        = "#C77CFF",
  threshold     = "#666666",
  left_region   = "#4D4D4D",
  right_region  = "#66C2A5"
)

event_levels <- c("Cutoff anchor", "Fixed leading-edge 5000", "Terminal start")

event_line_types <- c(
  "Cutoff anchor" = "dashed",
  "Fixed leading-edge 5000" = "dotdash",
  "Terminal start" = "dotted"
)

event_colors <- c(
  "Cutoff anchor" = COL$cutoff_anchor,
  "Fixed leading-edge 5000" = COL$fixed_5000,
  "Terminal start" = COL$terminal
)

event_shapes <- c(
  "Cutoff anchor" = 16,
  "Fixed leading-edge 5000" = 18,
  "Terminal start" = 1
)

deriv_colors <- c(
  "First derivative d1" = COL$d1,
  "Second derivative d2" = COL$d2
)

nb_colors <- c(
  "Combined NB support" = COL$nb_support,
  "NB1 = log(1+mu)" = COL$nb1,
  "NB2 = log(1+variance-mu)" = COL$nb2,
  "log(alpha*mu)" = COL$alpha_mu,
  "log(NB2+1)-log(NB1+1)" = COL$nb_gap
)

# =============================================================================
# SAFE HELPERS
# =============================================================================

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_var <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::var(x)
}

safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_log_pos <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x) & x > 0
  out[keep] <- log(x[keep])
  out
}

safe_log1p_nonneg <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x) & x >= 0
  out[keep] <- log1p(x[keep])
  out
}

scale01 <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)

  rng <- range(x[ok], na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) <= .Machine$double.eps) {
    out[ok] <- 0
    return(out)
  }

  out[ok] <- (x[ok] - rng[1L]) / diff(rng)
  out
}

fmt_int <- function(x) {
  if (!is.finite(x)) return("NA")
  format(as.integer(round(x)), trim = TRUE, scientific = FALSE)
}

fmt_num <- function(x, digits = 3L) {
  if (!is.finite(x)) return("NA")
  format(round(x, digits), nsmall = digits, trim = TRUE, scientific = FALSE)
}

panel_left_x <- function(df) {
  xr <- range(df$rank, na.rm = TRUE)
  xr[1L] + 0.05 * diff(xr)
}

panel_right_x <- function(df) {
  xr <- range(df$rank, na.rm = TRUE)
  xr[1L] + 0.72 * diff(xr)
}

resolve_counts_file <- function(path_hint) {
  candidates <- unique(c(
    path_hint,
    file.path(getwd(), path_hint),
    file.path(getwd(), "data", basename(path_hint)),
    "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
    "/root/REAPER98632/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
  ))

  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) return(candidates[[1L]])

  found <- list.files(
    "/root/REAPER98632",
    pattern = "WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
    recursive = TRUE,
    full.names = TRUE
  )
  found <- found[file.exists(found)]
  if (length(found)) return(found[[1L]])

  stop("Count file not found.", call. = FALSE)
}

detect_feature_id_column <- function(df) {
  candidates <- c("OrigID", "feature_id", "FeatureID", "PAS", "pas_id", "GeneID", "gene_id", "id")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) return(hit[[1L]])
  names(df)[1L]
}

detect_gene_symbol_column <- function(df) {
  candidates <- c("Symbol", "symbol", "GeneSymbol", "gene_symbol", "Gene", "gene")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) return(hit[[1L]])
  NULL
}

# =============================================================================
# DATA LOADING
# =============================================================================

read_count_matrix <- function(path, sample_ids) {
  raw_df <- utils::read.csv(path, check.names = FALSE)

  feature_col <- detect_feature_id_column(raw_df)
  gene_col <- detect_gene_symbol_column(raw_df)
  sample_cols <- intersect(sample_ids, names(raw_df))

  if (!length(sample_cols)) {
    stop("No sample columns matched the metadata sample IDs.", call. = FALSE)
  }

  annot_df <- data.frame(
    feature_id = as.character(raw_df[[feature_col]])
  )
  annot_df$gene_symbol <- if (!is.null(gene_col)) as.character(raw_df[[gene_col]]) else annot_df$feature_id

  keep <- !is.na(annot_df$feature_id) & nzchar(annot_df$feature_id)
  annot_df <- annot_df[keep, , drop = FALSE]
  count_df <- raw_df[keep, sample_cols, drop = FALSE]

  count_mat <- as.matrix(count_df)
  mode(count_mat) <- "numeric"
  rownames(count_mat) <- annot_df$feature_id

  finite_rows <- rowSums(!is.finite(count_mat)) == 0L
  count_mat <- count_mat[finite_rows, , drop = FALSE]
  annot_df <- annot_df[finite_rows, , drop = FALSE]

  if (anyDuplicated(rownames(count_mat))) {
    count_mat <- rowsum(count_mat, rownames(count_mat), reorder = FALSE)
    annot_df <- annot_df[match(rownames(count_mat), annot_df$feature_id), , drop = FALSE]
  }

  list(
    count_matrix = count_mat,
    annotation = annot_df
  )
}

subset_comparison_counts <- function(count_matrix, comparison_row) {
  trt_ids <- meta_all$id[grepl(paste0("^", comparison_row$treatment_prefix, "_"), meta_all$id)]
  ctrl_ids <- meta_all$id[grepl(paste0("^", comparison_row$control_prefix, "_"), meta_all$id)]

  keep_ids <- c(ctrl_ids, trt_ids)
  missing_ids <- setdiff(keep_ids, colnames(count_matrix))
  if (length(missing_ids)) {
    stop(
      paste0(
        "Missing samples for comparison ",
        comparison_row$comparison_name,
        ": ",
        paste(missing_ids, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  list(
    counts = count_matrix[, keep_ids, drop = FALSE],
    ctrl_ids = ctrl_ids,
    trt_ids = trt_ids
  )
}

# =============================================================================
# EVS RANKING
# =============================================================================

compute_abs_pc1_loadings <- function(count_mat, sample_ids) {
  mat <- count_mat[, sample_ids, drop = FALSE]
  mat <- log2(mat + 1)
  mat <- t(mat)

  if (nrow(mat) < 2L || ncol(mat) < 2L) {
    return(rep(NA_real_, ncol(mat)))
  }

  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  abs(pca$rotation[, 1L])
}

build_feature_table <- function(count_mat, ctrl_ids, trt_ids, annot_df) {
  ctrl_loading <- compute_abs_pc1_loadings(count_mat, ctrl_ids)
  trt_loading  <- compute_abs_pc1_loadings(count_mat, trt_ids)

  ctrl_mu  <- apply(count_mat[, ctrl_ids, drop = FALSE], 1L, safe_mean)
  trt_mu   <- apply(count_mat[, trt_ids, drop = FALSE], 1L, safe_mean)

  ctrl_var <- apply(count_mat[, ctrl_ids, drop = FALSE], 1L, safe_var)
  trt_var  <- apply(count_mat[, trt_ids, drop = FALSE], 1L, safe_var)

  tibble(
    feature_id = rownames(count_mat),
    gene_symbol = annot_df$gene_symbol[match(rownames(count_mat), annot_df$feature_id)],

    ctrl_abs_loading = ctrl_loading,
    trt_abs_loading = trt_loading,

    ctrl_rank = rank(ctrl_loading, ties.method = "first"),
    trt_rank = rank(trt_loading, ties.method = "first"),

    mu_ctrl = ctrl_mu,
    mu_trt = trt_mu,

    var_ctrl = ctrl_var,
    var_trt = trt_var
  ) %>%
    mutate(
      nb1_ctrl = mu_ctrl,
      nb1_trt = mu_trt,
      nb2_ctrl = pmax(var_ctrl - mu_ctrl, 0),
      nb2_trt = pmax(var_trt - mu_trt, 0),
      alpha_mu_ctrl = ifelse(mu_ctrl > 0 & (var_ctrl - mu_ctrl) > 0, (var_ctrl - mu_ctrl) / mu_ctrl, NA_real_),
      alpha_mu_trt  = ifelse(mu_trt  > 0 & (var_trt  - mu_trt ) > 0, (var_trt  - mu_trt ) / mu_trt , NA_real_)
    )
}

build_rank_df <- function(feature_tbl, arm = c("control", "treatment")) {
  arm <- match.arg(arm)

  if (arm == "control") {
    df <- feature_tbl %>%
      transmute(
        feature_id,
        gene_symbol,
        rank = ctrl_rank,
        abs_loading = ctrl_abs_loading,
        mu = mu_ctrl,
        variance = var_ctrl,
        nb1 = nb1_ctrl,
        nb2 = nb2_ctrl,
        alpha_mu = alpha_mu_ctrl
      )
  } else {
    df <- feature_tbl %>%
      transmute(
        feature_id,
        gene_symbol,
        rank = trt_rank,
        abs_loading = trt_abs_loading,
        mu = mu_trt,
        variance = var_trt,
        nb1 = nb1_trt,
        nb2 = nb2_trt,
        alpha_mu = alpha_mu_trt
      )
  }

  df <- df %>% arrange(rank)

  df$log_variance <- safe_log1p_nonneg(df$variance)
  df$log_nb1 <- safe_log1p_nonneg(df$nb1)
  df$log_nb2 <- safe_log1p_nonneg(df$nb2)
  df$log_alpha_mu <- safe_log_pos(df$alpha_mu)
  df$nb_gap <- df$log_nb2 - df$log_nb1

  df
}

# =============================================================================
# SECOND DERIVATIVE ZERO-CROSSINGS
# =============================================================================

find_d2_zero_crossings <- function(d2, tol = 1e-12) {
  x <- as.numeric(d2)
  n <- length(x)
  if (n < 2L) return(integer(0))

  x[!is.finite(x)] <- NA_real_
  x[is.finite(x) & abs(x) <= tol] <- 0

  out <- integer(0)

  for (i in 2:n) {
    a <- x[i - 1L]
    b <- x[i]

    if (!is.finite(a) || !is.finite(b)) next

    if (a == 0) out <- c(out, i - 1L)
    if (b == 0) out <- c(out, i)

    if ((a < 0 && b > 0) || (a > 0 && b < 0)) {
      out <- c(out, i)
    }
  }

  sort(unique(out))
}

summarize_zero_tbl <- function(df, idx) {
  if (!length(idx)) {
    return(data.frame(
      index = integer(0),
      rank = numeric(0),
      d1_here = numeric(0),
      d2_here = numeric(0),
      var_fit_here = numeric(0)
    ))
  }

  data.frame(
    index = idx,
    rank = df$rank[idx],
    d1_here = df$d1[idx],
    d2_here = df$d2[idx],
    var_fit_here = df$var_fit[idx]
  )
}

# =============================================================================
# GEOMETRIC CUT-OFF SELECTION
# =============================================================================

select_cutoff_interval <- function(rank_df,
                                   fixed_leading_edge_k = 5000L,
                                   variance_spline_spar = 0.60,
                                   zero_tol = 1e-12) {
  df <- rank_df

  ok <- is.finite(df$rank) & is.finite(df$log_variance)
  if (sum(ok) < 10L) {
    stop("Not enough finite points to fit the variance spline.", call. = FALSE)
  }

  ss <- smooth.spline(
    x = df$rank[ok],
    y = df$log_variance[ok],
    spar = variance_spline_spar
  )

  df$var_fit <- NA_real_
  df$d1 <- NA_real_
  df$d2 <- NA_real_

  df$var_fit[ok] <- predict(ss, x = df$rank[ok], deriv = 0)$y
  df$d1[ok] <- predict(ss, x = df$rank[ok], deriv = 1)$y
  df$d2[ok] <- predict(ss, x = df$rank[ok], deriv = 2)$y

  total_features <- nrow(df)

  fixed_reference_rank_nominal <- max(1L, total_features - fixed_leading_edge_k + 1L)
  fixed_reference_index <- which.min(abs(df$rank - fixed_reference_rank_nominal))
  fixed_reference_rank <- df$rank[fixed_reference_index]

  zero_idx <- find_d2_zero_crossings(df$d2, tol = zero_tol)
  zero_tbl <- summarize_zero_tbl(df, zero_idx)

  left_zeros <- zero_tbl[zero_tbl$rank < fixed_reference_rank, , drop = FALSE]
  right_zeros <- zero_tbl[zero_tbl$rank > fixed_reference_rank, , drop = FALSE]

  if (nrow(left_zeros) < 1L) {
    stop("No d2 zero immediately LEFT of fixed leading-edge-5000 mark.", call. = FALSE)
  }
  if (nrow(right_zeros) < 1L) {
    stop("No d2 zero immediately RIGHT of fixed leading-edge-5000 mark.", call. = FALSE)
  }

  cutoff_anchor_row <- left_zeros[which.max(left_zeros$rank), , drop = FALSE]
  terminal_start_row <- right_zeros[which.min(right_zeros$rank), , drop = FALSE]

  cutoff_anchor_index <- cutoff_anchor_row$index
  terminal_start_index <- terminal_start_row$index

  interval_min_rank <- df$rank[cutoff_anchor_index]
  interval_max_rank <- df$rank[terminal_start_index]

  pre_evs_remainder_size <- cutoff_anchor_index - 1L
  pre_evs_leading_edge_size <- total_features - cutoff_anchor_index + 1L

  # RIGHT region = full leading edge from cutoff anchor to the end
  right_idx <- seq.int(cutoff_anchor_index, total_features)
  right_n <- length(right_idx)

  # LEFT region = same number of genes immediately left of cutoff anchor
  left_end <- cutoff_anchor_index - 1L
  left_start <- max(1L, left_end - right_n + 1L)
  left_idx <- if (left_end >= left_start) seq.int(left_start, left_end) else integer(0)

  df$nb_support <- rowMeans(
    cbind(
      scale01(df$log_nb2),
      scale01(df$nb_gap),
      scale01(df$log_alpha_mu)
    ),
    na.rm = TRUE
  )

  corrob <- data.frame(
    left_n = length(left_idx),
    right_n = length(right_idx),

    left_median_log_nb2 = if (length(left_idx)) safe_median(df$log_nb2[left_idx]) else NA_real_,
    right_median_log_nb2 = if (length(right_idx)) safe_median(df$log_nb2[right_idx]) else NA_real_,

    left_median_nb_gap = if (length(left_idx)) safe_median(df$nb_gap[left_idx]) else NA_real_,
    right_median_nb_gap = if (length(right_idx)) safe_median(df$nb_gap[right_idx]) else NA_real_,

    left_median_log_alpha_mu = if (length(left_idx)) safe_median(df$log_alpha_mu[left_idx]) else NA_real_,
    right_median_log_alpha_mu = if (length(right_idx)) safe_median(df$log_alpha_mu[right_idx]) else NA_real_,

    stringsAsFactors = FALSE
  )

  corrob$right_left_log_nb2_diff <- corrob$right_median_log_nb2 - corrob$left_median_log_nb2
  corrob$right_left_nb_gap_diff <- corrob$right_median_nb_gap - corrob$left_median_nb_gap
  corrob$right_left_log_alpha_mu_diff <- corrob$right_median_log_alpha_mu - corrob$left_median_log_alpha_mu

  nb_support_threshold <- safe_median(df$nb_support[right_idx])

  selected_points <- data.frame(
    event = factor(event_levels, levels = event_levels),
    index = c(cutoff_anchor_index, fixed_reference_index, terminal_start_index),
    rank = c(df$rank[cutoff_anchor_index], fixed_reference_rank, df$rank[terminal_start_index])
  )

  list(
    df = df,
    zero_tbl = zero_tbl,
    selected_points = selected_points,

    cutoff_anchor_index = cutoff_anchor_index,
    cutoff_anchor_rank = df$rank[cutoff_anchor_index],

    fixed_reference_index = fixed_reference_index,
    fixed_reference_rank = fixed_reference_rank,

    terminal_start_index = terminal_start_index,
    terminal_start_rank = df$rank[terminal_start_index],

    interval_min_rank = interval_min_rank,
    interval_max_rank = interval_max_rank,

    pre_evs_remainder_size = pre_evs_remainder_size,
    pre_evs_leading_edge_size = pre_evs_leading_edge_size,

    total_features = total_features,
    left_idx = left_idx,
    right_idx = right_idx,

    corrob = corrob,
    nb_support_threshold = nb_support_threshold
  )
}

# =============================================================================
# PLOT HELPERS
# =============================================================================

event_line_data <- function(sel) {
  data.frame(
    event = factor(event_levels, levels = event_levels),
    x = c(sel$cutoff_anchor_rank, sel$fixed_reference_rank, sel$terminal_start_rank)
  )
}

event_point_data <- function(sel, y_cutoff, y_reference, y_terminal) {
  data.frame(
    event = factor(event_levels, levels = event_levels),
    x = c(sel$cutoff_anchor_rank, sel$fixed_reference_rank, sel$terminal_start_rank),
    y = c(y_cutoff, y_reference, y_terminal)
  )
}

interval_band_layer <- function(sel) {
  annotate(
    "rect",
    xmin = sel$interval_min_rank,
    xmax = sel$interval_max_rank,
    ymin = -Inf,
    ymax = Inf,
    fill = COL$interval_fill,
    alpha = 0.20
  )
}

event_line_layers <- function(sel, show.legend = TRUE) {
  list(
    geom_vline(
      data = event_line_data(sel),
      aes(xintercept = x, color = event, linetype = event),
      linewidth = 0.85,
      show.legend = show.legend,
      inherit.aes = FALSE
    )
  )
}

event_point_layers <- function(point_df, show.legend = TRUE) {
  list(
    geom_point(
      data = point_df,
      aes(x = x, y = y, shape = event, color = event),
      size = 3.0,
      stroke = 1.0,
      show.legend = show.legend,
      inherit.aes = FALSE,
      fill = "white"
    )
  )
}

event_scales <- function() {
  list(
    scale_color_manual(values = event_colors, breaks = event_levels, drop = FALSE),
    scale_linetype_manual(values = event_line_types, breaks = event_levels, drop = FALSE),
    scale_shape_manual(values = event_shapes, breaks = event_levels, drop = FALSE)
  )
}

event_guides <- function() {
  guides(
    color = guide_legend(order = 1, title = "Event"),
    linetype = guide_legend(order = 1, title = "Event"),
    shape = guide_legend(order = 1, title = "Event")
  )
}

label_box_layer <- function(x, y, text, size = 2.15) {
  geom_label(
    data = data.frame(x = x, y = y, label = text),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = size,
    label.padding = unit(0.15, "lines"),
    label.r = unit(0.10, "lines"),
    linewidth = 0.25,
    fill = alpha("white", 0.96)
  )
}

base_theme <- function() {
  theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 12),
      plot.subtitle = element_text(size = 10),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin = margin(8, 55, 8, 18)
    )
}

# =============================================================================
# FIGURE CREATION
# =============================================================================

make_rank_panel <- function(sel, title_prefix, out_file) {
  df <- sel$df

  left_x <- panel_left_x(df)
  right_x <- panel_right_x(df)

  abs_ymax <- max(df$abs_loading, na.rm = TRUE)
  var_ymax <- max(df$var_fit, na.rm = TRUE)
  deriv_ymax <- max(abs(c(df$d1, df$d2)), na.rm = TRUE)
  nb_ymax <- max(c(df$nb_support, df$log_nb1, df$log_nb2, df$log_alpha_mu, df$nb_gap), na.rm = TRUE)

  abs_pts <- event_point_data(
    sel,
    y_cutoff = df$abs_loading[sel$cutoff_anchor_index],
    y_reference = df$abs_loading[sel$fixed_reference_index],
    y_terminal = df$abs_loading[sel$terminal_start_index]
  )

  var_pts <- event_point_data(
    sel,
    y_cutoff = df$var_fit[sel$cutoff_anchor_index],
    y_reference = df$var_fit[sel$fixed_reference_index],
    y_terminal = df$var_fit[sel$terminal_start_index]
  )

  deriv_pts <- event_point_data(
    sel,
    y_cutoff = df$d2[sel$cutoff_anchor_index],
    y_reference = df$d2[sel$fixed_reference_index],
    y_terminal = df$d2[sel$terminal_start_index]
  )

  nb_pts <- event_point_data(
    sel,
    y_cutoff = df$nb_support[sel$cutoff_anchor_index],
    y_reference = df$nb_support[sel$fixed_reference_index],
    y_terminal = df$nb_support[sel$terminal_start_index]
  )

  # PANEL 1
  p1 <- ggplot(df, aes(x = rank, y = abs_loading)) +
    interval_band_layer(sel) +
    geom_line(color = COL$abs_loading, linewidth = 1.0) +
    event_line_layers(sel) +
    event_point_layers(abs_pts) +
    event_scales() +
    event_guides() +
    labs(
      title = paste0(title_prefix, ": absolute PC1 loading series"),
      subtitle = "Leading edge is on the RIGHT",
      x = "EVS rank",
      y = "|PC1 loading|"
    ) +
    base_theme() +
    label_box_layer(
      left_x,
      abs_ymax,
      paste(
        "Absolute loading panel",
        "Leading edge is on the RIGHT",
        "The grey band is the final geometric interval",
        "The legend gives the exact event colors and line types",
        sep = "\n"
      )
    ) +
    label_box_layer(
      right_x,
      abs_ymax,
      paste0(
        "Cutoff anchor rank = ", fmt_int(sel$cutoff_anchor_rank),
        "\nReference rank (5000 from right) = ", fmt_int(sel$fixed_reference_rank),
        "\nTerminal start rank = ", fmt_int(sel$terminal_start_rank),
        "\nFinal interval = [", fmt_int(sel$interval_min_rank), ", ", fmt_int(sel$interval_max_rank), "]",
        "\nPre-EVS remainder = ", fmt_int(sel$pre_evs_remainder_size),
        "\nPre-EVS leading edge = ", fmt_int(sel$pre_evs_leading_edge_size)
      )
    )

  # PANEL 2
  p2 <- ggplot(df, aes(x = rank, y = var_fit)) +
    interval_band_layer(sel) +
    geom_line(color = COL$variance_fit, linewidth = 1.0) +
    event_line_layers(sel) +
    event_point_layers(var_pts) +
    event_scales() +
    event_guides() +
    labs(
      title = paste0(title_prefix, ": smoothed empirical variance curve"),
      subtitle = "The geometric points are marked directly on the fitted curve",
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    base_theme() +
    label_box_layer(
      left_x,
      var_ymax,
      paste(
        "Variance curve method",
        "A smoothing spline is fit to empirical log(1+variance).",
        "Cutoff anchor = nearest d2 zero immediately LEFT of the fixed leading-edge 5000 mark.",
        "Terminal start = nearest d2 zero immediately RIGHT of the fixed leading-edge 5000 mark.",
        sep = "\n"
      ),
      size = 2.00
    )

  # PANEL 3
  deriv_long <- bind_rows(
    data.frame(rank = df$rank, value = df$d1, quantity = "First derivative d1"),
    data.frame(rank = df$rank, value = df$d2, quantity = "Second derivative d2")
  )
  deriv_long$quantity <- factor(
    deriv_long$quantity,
    levels = c("First derivative d1", "Second derivative d2")
  )

  p3 <- ggplot(deriv_long, aes(x = rank, y = value)) +
    interval_band_layer(sel) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.55) +
    geom_line(
      data = deriv_long[deriv_long$quantity == "First derivative d1", , drop = FALSE],
      aes(color = quantity),
      linewidth = 0.95
    ) +
    geom_line(
      data = deriv_long[deriv_long$quantity == "Second derivative d2", , drop = FALSE],
      aes(color = quantity),
      linewidth = 0.95
    ) +
    event_line_layers(sel) +
    event_point_layers(deriv_pts) +
    scale_color_manual(
      values = c(deriv_colors, event_colors),
      breaks = c(names(deriv_colors), event_levels),
      drop = FALSE
    ) +
    scale_linetype_manual(values = event_line_types, breaks = event_levels, drop = FALSE) +
    scale_shape_manual(values = event_shapes, breaks = event_levels, drop = FALSE) +
    labs(
      title = paste0(title_prefix, ": derivative support"),
      subtitle = "The selected d2 zero-crossings around the fixed leading-edge 5000 mark define the geometric interval",
      x = "EVS rank",
      y = "Derivative value",
      color = NULL,
      linetype = "Event",
      shape = "Event"
    ) +
    base_theme() +
    label_box_layer(
      left_x,
      deriv_ymax,
      paste(
        "Derivative method",
        "All d2 zeros come only from sign changes or exact zeros.",
        "No amplitude threshold is used.",
        "No run-length filter is used.",
        "The fixed leading-edge 5000 mark is inside the final interval.",
        sep = "\n"
      ),
      size = 2.00
    )

  # PANEL 4
  nb_long <- bind_rows(
    data.frame(rank = df$rank, value = df$nb_support, quantity = "Combined NB support"),
    data.frame(rank = df$rank, value = df$log_nb1, quantity = "NB1 = log(1+mu)"),
    data.frame(rank = df$rank, value = df$log_nb2, quantity = "NB2 = log(1+variance-mu)"),
    data.frame(rank = df$rank, value = df$log_alpha_mu, quantity = "log(alpha*mu)"),
    data.frame(rank = df$rank, value = df$nb_gap, quantity = "log(NB2+1)-log(NB1+1)")
  )
  nb_long$quantity <- factor(nb_long$quantity, levels = names(nb_colors))

  left_region_min <- if (length(sel$left_idx)) min(df$rank[sel$left_idx]) else NA_real_
  left_region_max <- if (length(sel$left_idx)) max(df$rank[sel$left_idx]) else NA_real_
  right_region_min <- if (length(sel$right_idx)) min(df$rank[sel$right_idx]) else NA_real_
  right_region_max <- if (length(sel$right_idx)) max(df$rank[sel$right_idx]) else NA_real_

  p4 <- ggplot(nb_long, aes(x = rank, y = value)) +
    annotate(
      "rect",
      xmin = left_region_min, xmax = left_region_max,
      ymin = -Inf, ymax = Inf,
      fill = COL$left_region, alpha = 0.08
    ) +
    annotate(
      "rect",
      xmin = right_region_min, xmax = right_region_max,
      ymin = -Inf, ymax = Inf,
      fill = COL$right_region, alpha = 0.08
    ) +
    interval_band_layer(sel) +
    geom_hline(yintercept = sel$nb_support_threshold, color = COL$threshold, linetype = 3, linewidth = 0.65) +
    geom_line(aes(color = quantity), linewidth = 0.90) +
    event_line_layers(sel) +
    event_point_layers(nb_pts) +
    scale_color_manual(
      values = c(nb_colors, event_colors),
      breaks = c(names(nb_colors), event_levels),
      drop = FALSE
    ) +
    scale_linetype_manual(values = event_line_types, breaks = event_levels, drop = FALSE) +
    scale_shape_manual(values = event_shapes, breaks = event_levels, drop = FALSE) +
    labs(
      title = paste0(title_prefix, ": NB1 / NB2 / alpha*mu support"),
      subtitle = "The full RIGHT leading edge from cutoff anchor to the end is compared against a matched LEFT region of equal size",
      x = "EVS rank",
      y = "Support value",
      color = NULL,
      linetype = "Event",
      shape = "Event"
    ) +
    base_theme() +
    label_box_layer(
      left_x,
      nb_ymax,
      paste(
        "NB corroboration",
        "The full RIGHT leading edge from cutoff anchor to the end is the NB2 side.",
        "The matched LEFT comparison region has the same number of genes immediately left of the cutoff anchor.",
        "These NB values corroborate but do not define the geometric split.",
        sep = "\n"
      ),
      size = 1.95
    ) +
    label_box_layer(
      right_x,
      nb_ymax,
      paste0(
        "Left median log(NB2) = ", fmt_num(sel$corrob$left_median_log_nb2),
        "\nRight median log(NB2) = ", fmt_num(sel$corrob$right_median_log_nb2),
        "\nRight-left log(NB2) diff = ", fmt_num(sel$corrob$right_left_log_nb2_diff),
        "\nLeft median NB2-NB1 contrast = ", fmt_num(sel$corrob$left_median_nb_gap),
        "\nRight median NB2-NB1 contrast = ", fmt_num(sel$corrob$right_median_nb_gap),
        "\nRight-left NB2-NB1 diff = ", fmt_num(sel$corrob$right_left_nb_gap_diff),
        "\nLeft median log(alpha*mu) = ", fmt_num(sel$corrob$left_median_log_alpha_mu),
        "\nRight median log(alpha*mu) = ", fmt_num(sel$corrob$right_median_log_alpha_mu),
        "\nRight-left log(alpha*mu) diff = ", fmt_num(sel$corrob$right_left_log_alpha_mu_diff)
      ),
      size = 1.90
    )

  # PANEL 5
  nb_sum_df <- data.frame(rank = df$rank, value = df$nb_support)

  nb_sum_pts <- event_point_data(
    sel,
    y_cutoff = df$nb_support[sel$cutoff_anchor_index],
    y_reference = df$nb_support[sel$fixed_reference_index],
    y_terminal = df$nb_support[sel$terminal_start_index]
  )

  p5 <- ggplot(nb_sum_df, aes(x = rank, y = value)) +
    annotate(
      "rect",
      xmin = left_region_min, xmax = left_region_max,
      ymin = -Inf, ymax = Inf,
      fill = COL$left_region, alpha = 0.08
    ) +
    annotate(
      "rect",
      xmin = right_region_min, xmax = right_region_max,
      ymin = -Inf, ymax = Inf,
      fill = COL$right_region, alpha = 0.08
    ) +
    interval_band_layer(sel) +
    geom_hline(yintercept = sel$nb_support_threshold, color = COL$threshold, linetype = 3, linewidth = 0.65) +
    geom_line(color = COL$nb_support, linewidth = 1.0) +
    event_line_layers(sel) +
    event_point_layers(nb_sum_pts) +
    event_scales() +
    event_guides() +
    labs(
      title = paste0(title_prefix, ": NB-supported summary"),
      subtitle = "The geometric interval remains primary; NB support remains corroborative",
      x = "EVS rank",
      y = "Combined NB support"
    ) +
    base_theme() +
    label_box_layer(
      right_x,
      max(nb_sum_df$value, na.rm = TRUE),
      paste0(
        "NB-supported summary",
        "\nAnchor = ", fmt_int(sel$cutoff_anchor_rank),
        "\nReference = ", fmt_int(sel$fixed_reference_rank),
        "\nTerminal = ", fmt_int(sel$terminal_start_rank),
        "\nThreshold = ", fmt_num(sel$nb_support_threshold)
      ),
      size = 2.00
    )

  g <- arrangeGrob(
    p1, p2, p3, p4, p5,
    ncol = 1,
    heights = c(1, 1, 1, 1.15, 1)
  )

  ggplot2::ggsave(
    filename = out_file,
    plot = g,
    width = 16,
    height = 22,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
}

# =============================================================================
# CSV WRITERS
# =============================================================================

write_rank_series_csv <- function(sel, out_file) {
  out <- sel$df %>%
    transmute(
      rank,
      feature_id,
      gene_symbol,
      abs_loading,
      mu,
      variance,
      nb1,
      nb2,
      alpha_mu,
      log_variance,
      log_nb1,
      log_nb2,
      log_alpha_mu,
      nb_gap,
      var_fit,
      d1,
      d2,
      nb_support
    )

  utils::write.csv(out, out_file, row.names = FALSE)
}

write_zero_crossings_csv <- function(sel, out_file) {
  utils::write.csv(sel$zero_tbl, out_file, row.names = FALSE)
}

write_selected_two_zero_crossings_csv <- function(sel, out_file) {
  utils::write.csv(sel$selected_points, out_file, row.names = FALSE)
}

write_feature_level_metrics_csv <- function(sel, out_file) {
  df <- sel$df

  region <- rep("middle", nrow(df))
  if (length(sel$left_idx)) region[sel$left_idx] <- "matched_left"
  if (length(sel$right_idx)) region[sel$right_idx] <- "right_leading_edge"

  out <- df %>%
    mutate(region = region) %>%
    transmute(
      rank,
      feature_id,
      gene_symbol,
      region,
      abs_loading,
      mu,
      variance,
      nb1,
      nb2,
      alpha_mu,
      log_nb1,
      log_nb2,
      log_alpha_mu,
      nb_gap,
      nb_support,
      var_fit,
      d1,
      d2
    )

  utils::write.csv(out, out_file, row.names = FALSE)
}

write_cutoff_summary_csv <- function(summary_tbl, out_file) {
  utils::write.csv(summary_tbl, out_file, row.names = FALSE)
}

# =============================================================================
# MAIN
# =============================================================================

count_file <- resolve_counts_file(count_file_hint)
message("Using count file: ", count_file)

count_obj <- read_count_matrix(count_file, meta_all$id)
count_matrix <- count_obj$count_matrix
annot_df <- count_obj$annotation

message(
  "Count matrix dimensions: ",
  nrow(count_matrix), " features x ", ncol(count_matrix), " samples"
)

all_summary_rows <- list()

for (i in seq_len(nrow(comparison_table))) {
  cmp <- comparison_table[i, , drop = FALSE]
  cmp_name <- cmp$comparison_name

  message("Processing comparison: ", cmp_name)

  cmp_counts <- subset_comparison_counts(count_matrix, cmp)
  feature_tbl <- build_feature_table(
    count_mat = cmp_counts$counts,
    ctrl_ids = cmp_counts$ctrl_ids,
    trt_ids = cmp_counts$trt_ids,
    annot_df = annot_df
  )

  out_dir <- file.path(output_root, paste0(cmp_name, "_cutoff_folder"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (arm in c("control", "treatment")) {
    message("Selecting ", arm, " cutoff for ", cmp_name)

    rank_df <- build_rank_df(feature_tbl, arm = arm)

    sel <- select_cutoff_interval(
      rank_df = rank_df,
      fixed_leading_edge_k = fixed_leading_edge_k,
      variance_spline_spar = variance_spline_spar,
      zero_tol = zero_tol
    )

    message(
      tools::toTitleCase(arm), " cutoff anchor rank: ", fmt_int(sel$cutoff_anchor_rank),
      " | reference rank: ", fmt_int(sel$fixed_reference_rank),
      " | terminal anchor: ", fmt_int(sel$terminal_start_rank),
      " | interval: [", fmt_int(sel$interval_min_rank), ", ", fmt_int(sel$interval_max_rank), "]",
      " | pre-EVS remainder: ", fmt_int(sel$pre_evs_remainder_size),
      " | pre-EVS leading edge: ", fmt_int(sel$pre_evs_leading_edge_size)
    )

    panel_png <- file.path(out_dir, paste0(cmp_name, "_", arm, "_rank_panel.png"))
    make_rank_panel(
      sel = sel,
      title_prefix = paste0(cmp_name, " ", arm),
      out_file = panel_png
    )

    write_rank_series_csv(
      sel = sel,
      out_file = file.path(out_dir, paste0(cmp_name, "_", arm, "_rank_series.csv"))
    )

    write_zero_crossings_csv(
      sel = sel,
      out_file = file.path(out_dir, paste0(cmp_name, "_", arm, "_zero_crossings_all.csv"))
    )

    write_selected_two_zero_crossings_csv(
      sel = sel,
      out_file = file.path(out_dir, paste0(cmp_name, "_", arm, "_selected_two_zero_crossings.csv"))
    )

    write_feature_level_metrics_csv(
      sel = sel,
      out_file = file.path(out_dir, paste0(cmp_name, "_", arm, "_feature_level_metrics.csv"))
    )

    arm_summary <- data.frame(
      comparison_name = cmp_name,
      arm = arm,

      cutoff_anchor_rank = sel$cutoff_anchor_rank,
      fixed_reference_rank = sel$fixed_reference_rank,
      terminal_start_rank = sel$terminal_start_rank,

      interval_min_rank = sel$interval_min_rank,
      interval_max_rank = sel$interval_max_rank,

      pre_evs_remainder_size = sel$pre_evs_remainder_size,
      pre_evs_leading_edge_size = sel$pre_evs_leading_edge_size,

      left_n = sel$corrob$left_n,
      right_n = sel$corrob$right_n,

      left_median_log_nb2 = sel$corrob$left_median_log_nb2,
      right_median_log_nb2 = sel$corrob$right_median_log_nb2,
      right_left_log_nb2_diff = sel$corrob$right_left_log_nb2_diff,

      left_median_nb_gap = sel$corrob$left_median_nb_gap,
      right_median_nb_gap = sel$corrob$right_median_nb_gap,
      right_left_nb_gap_diff = sel$corrob$right_left_nb_gap_diff,

      left_median_log_alpha_mu = sel$corrob$left_median_log_alpha_mu,
      right_median_log_alpha_mu = sel$corrob$right_median_log_alpha_mu,
      right_left_log_alpha_mu_diff = sel$corrob$right_left_log_alpha_mu_diff,

      center_nb_support = sel$df$nb_support[sel$cutoff_anchor_index],
      nb_support_threshold = sel$nb_support_threshold
    )

    all_summary_rows[[length(all_summary_rows) + 1L]] <- arm_summary
  }
}

overall_summary <- bind_rows(all_summary_rows)

write_cutoff_summary_csv(
  summary_tbl = overall_summary,
  out_file = file.path(output_root, "overall_cutoff_summary.csv")
)

message("Done. Outputs written to: ", output_root)
