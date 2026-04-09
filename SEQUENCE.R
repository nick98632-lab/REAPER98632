# =============================================================================
# FINAL MANUSCRIPT SCRIPT
# PRE-EVS LEADING-EDGE CUTOFF SELECTION USING VARIANCE GEOMETRY + NB CORROBORATION
# =============================================================================
#
# METHOD OVERVIEW
#
# The EVS-ranked axis is ordered from LEFT to RIGHT by increasing absolute PC1
# loading within each arm separately. Therefore the leading edge is always on
# the RIGHT side of the ranked axis.
#
# For each arm independently:
#
# 1. A smooth empirical variance curve is fit over EVS rank using log(1+variance).
# 2. The first derivative d1 and second derivative d2 of that same smooth are
#    evaluated on the full EVS rank grid.
# 3. A fixed leading-edge reference mark is placed 5000 genes from the RIGHT:
#       reference_rank = total_features - 5000 + 1
# 4. All d2 zero-crossings are identified from sign changes or exact zeros.
# 5. The cutoff anchor is defined as the nearest d2 zero immediately LEFT of
#    the fixed leading-edge 5000 reference mark.
# 6. The terminal start is defined as the nearest d2 zero immediately RIGHT of
#    the same fixed leading-edge 5000 reference mark.
# 7. The final geometric interval is therefore:
#       [cutoff_anchor, terminal_start]
#    so that the fixed leading-edge 5000 reference lies inside the interval.
#
# NB corroboration is then computed relative to the cutoff anchor:
#
# - The RIGHT region is the entire leading edge from cutoff_anchor to the end
#   of the ranked axis.
# - The LEFT matched region contains the same number of genes immediately to
#   the left of the cutoff anchor.
# - The following quantities are summarized on both sides:
#       NB1 = mu
#       NB2 = variance - mu
#       NB2 - NB1 contrast on the log scale
#       log(alpha*mu)
# - These NB quantities corroborate the split but do not define the geometric
#   cutoff itself.
#
# FIGURE RULES
#
# - Filled circle = cutoff anchor
# - Diamond       = fixed leading-edge 5000 reference
# - Open circle   = terminal start
# - Dashed line    = cutoff anchor
# - Dot-dash line  = fixed leading-edge 5000 reference
# - Dotted line    = terminal start
# - Grey band      = final interval [cutoff_anchor, terminal_start]
#
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(gridExtra)
  library(grid)
})

options(warn = 1)

# =============================================================================
# USER INPUTS
# =============================================================================

repo_dir <- getwd()
count_file_hint <- file.path(repo_dir, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")

output_root <- file.path(
  repo_dir,
  "exports",
  "variance_derivative_nb_range_leadingedge5000_rewritten_final"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

leading_edge_fixed_k <- 5000L
spline_spar <- 0.60
zero_tol <- 1e-12

comparison_table <- data.frame(
  comparison_name  = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  treatment_prefix = c("R0", "R2", "R4", "R8"),
  control_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
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
  ),
  stringsAsFactors = FALSE
)
rownames(meta_all) <- meta_all$id
meta_all$condition <- factor(meta_all$condition, levels = c("untrt", "trt"))

# =============================================================================
# COLOR SYSTEM
# =============================================================================

COL <- c(
  abs_loading = "#111111",
  variance_fit = "#111111",
  d1 = "#00B4D8",
  d2 = "#D1495B",
  nb_support = "#111111",
  nb1 = "#C49A00",
  nb2 = "#2FB47C",
  amu = "#1F78FF",
  nb_gap = "#C77CFF",
  cutoff = "#111111",
  ref5000 = "#8C510A",
  terminal = "#B35806",
  interval = "#BDBDBD"
)

# =============================================================================
# BASIC HELPERS
# =============================================================================

resolve_counts_file <- function(path_hint) {
  candidates <- unique(c(
    path_hint,
    file.path(getwd(), path_hint),
    file.path(getwd(), "data", basename(path_hint)),
    "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
    "/root/REAPER98632/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
  ))

  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) > 0L) return(candidates[[1]])

  found <- list.files(
    path = "/root/REAPER98632",
    pattern = "WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
    recursive = TRUE,
    full.names = TRUE
  )
  found <- found[file.exists(found)]
  if (length(found) > 0L) return(found[[1]])

  stop("Count file not found.", call. = FALSE)
}

detect_feature_id_column <- function(df) {
  candidates <- c("OrigID", "feature_id", "FeatureID", "PAS", "pas_id", "GeneID", "gene_id", "id")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0L) return(hit[1L])
  names(df)[1L]
}

detect_gene_symbol_column <- function(df) {
  candidates <- c("Symbol", "symbol", "GeneSymbol", "gene_symbol", "Gene", "gene")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0L) return(hit[1L])
  NULL
}

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
  keep <- is.finite(x) & (x > 0)
  out[keep] <- log(x[keep])
  out
}

safe_log1p_nonneg <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x) & (x >= 0)
  out[keep] <- log1p(x[keep])
  out
}

scale01 <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)

  rng <- range(x[ok], na.rm = TRUE)
  if (!is.finite(rng[1L]) || !is.finite(rng[2L]) || abs(rng[2L] - rng[1L]) < .Machine$double.eps) {
    out[ok] <- 0
    return(out)
  }

  out[ok] <- (x[ok] - rng[1L]) / (rng[2L] - rng[1L])
  out
}

fmt_num <- function(x, digits = 3L) {
  if (!is.finite(x)) return("NA")
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

# =============================================================================
# DATA LOADING
# =============================================================================

read_count_matrix <- function(path, sample_ids) {
  raw_df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  feature_col <- detect_feature_id_column(raw_df)
  gene_col <- detect_gene_symbol_column(raw_df)
  sample_cols <- intersect(sample_ids, names(raw_df))

  if (length(sample_cols) == 0L) {
    stop("No count columns matched metadata sample IDs.", call. = FALSE)
  }

  annot_df <- data.frame(
    feature_id = as.character(raw_df[[feature_col]]),
    stringsAsFactors = FALSE
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

subset_comparison_counts <- function(count_matrix, comparison_row, meta_all) {
  trt_ids <- meta_all$id[grepl(paste0("^", comparison_row$treatment_prefix, "_"), meta_all$id)]
  ctrl_ids <- meta_all$id[grepl(paste0("^", comparison_row$control_prefix, "_"), meta_all$id)]

  keep_ids <- c(trt_ids, ctrl_ids)
  missing_ids <- setdiff(keep_ids, colnames(count_matrix))
  if (length(missing_ids) > 0L) {
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

  if (nrow(mat) < 2L) {
    return(rep(NA_real_, ncol(mat)))
  }

  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  abs(pca$rotation[, 1L])
}

build_feature_table <- function(count_mat, ctrl_ids, trt_ids, annot_df) {
  ctrl_load <- compute_abs_pc1_loadings(count_mat, ctrl_ids)
  trt_load  <- compute_abs_pc1_loadings(count_mat, trt_ids)

  ctrl_mu  <- apply(count_mat[, ctrl_ids, drop = FALSE], 1L, safe_mean)
  trt_mu   <- apply(count_mat[, trt_ids, drop = FALSE], 1L, safe_mean)
  ctrl_var <- apply(count_mat[, ctrl_ids, drop = FALSE], 1L, safe_var)
  trt_var  <- apply(count_mat[, trt_ids, drop = FALSE], 1L, safe_var)

  tibble(
    feature_id = rownames(count_mat),
    gene_symbol = annot_df$gene_symbol[match(rownames(count_mat), annot_df$feature_id)],

    ctrl_abs_loading = ctrl_load,
    trt_abs_loading  = trt_load,

    ctrl_rank = rank(ctrl_load, ties.method = "first"),
    trt_rank  = rank(trt_load, ties.method = "first"),

    mu_ctrl = ctrl_mu,
    mu_trt  = trt_mu,

    var_ctrl = ctrl_var,
    var_trt  = trt_var
  ) %>%
    mutate(
      nb1_ctrl = mu_ctrl,
      nb1_trt  = mu_trt,
      nb2_ctrl = pmax(var_ctrl - mu_ctrl, 0),
      nb2_trt  = pmax(var_trt - mu_trt, 0),
      alpha_mu_ctrl = ifelse(mu_ctrl > 0, (var_ctrl - mu_ctrl) / mu_ctrl, NA_real_),
      alpha_mu_trt  = ifelse(mu_trt > 0, (var_trt - mu_trt) / mu_trt, NA_real_)
    )
}

build_rank_df <- function(feature_tbl, arm = c("control", "treatment")) {
  arm <- match.arg(arm)

  if (arm == "control") {
    out <- feature_tbl %>%
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
      ) %>%
      arrange(rank)
  } else {
    out <- feature_tbl %>%
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
      ) %>%
      arrange(rank)
  }

  out$log_variance <- safe_log1p_nonneg(out$variance)
  out$log_nb1 <- safe_log1p_nonneg(out$nb1)
  out$log_nb2 <- safe_log1p_nonneg(out$nb2)
  out$log_alpha_mu <- safe_log_pos(out$alpha_mu)
  out$nb_gap <- out$log_nb2 - out$log_nb1

  out
}

# =============================================================================
# GEOMETRIC CUT-OFF SELECTION
# =============================================================================

find_d2_zero_crossings <- function(d2_vec, zero_tol = 1e-12) {
  x <- as.numeric(d2_vec)
  n <- length(x)
  if (n < 2L) return(integer(0))

  x[!is.finite(x)] <- NA_real_
  x[is.finite(x) & abs(x) <= zero_tol] <- 0
  s <- sign(x)
  s[!is.finite(s)] <- NA_integer_

  out <- integer(0)

  for (i in 2:n) {
    x0 <- x[i - 1L]
    x1 <- x[i]
    s0 <- s[i - 1L]
    s1 <- s[i]

    if (!is.finite(x0) || !is.finite(x1) || is.na(s0) || is.na(s1)) next

    if (s0 == 0 && s1 != 0) {
      out <- c(out, i - 1L)
    } else if (s0 != 0 && s1 == 0) {
      out <- c(out, i)
    } else if (s0 != s1) {
      out <- c(out, i)
    }
  }

  sort(unique(out))
}

summarize_zero_crossings <- function(df, idx) {
  if (!length(idx)) {
    return(data.frame(
      index = integer(0),
      rank = numeric(0),
      d1_here = numeric(0),
      d2_here = numeric(0),
      var_fit_here = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    index = idx,
    rank = df$rank[idx],
    d1_here = df$d1[idx],
    d2_here = df$d2[idx],
    var_fit_here = df$var_fit[idx],
    stringsAsFactors = FALSE
  )
}

select_interval_and_nb <- function(rank_df,
                                   spline_spar = 0.60,
                                   zero_tol = 1e-12,
                                   leading_edge_fixed_k = 5000L) {
  df <- rank_df

  ok <- is.finite(df$rank) & is.finite(df$log_variance)
  if (sum(ok) < 10L) {
    stop("Not enough finite variance points for smoothing.", call. = FALSE)
  }

  ss <- smooth.spline(x = df$rank[ok], y = df$log_variance[ok], spar = spline_spar)

  df$var_fit <- NA_real_
  df$d1 <- NA_real_
  df$d2 <- NA_real_

  df$var_fit[ok] <- predict(ss, x = df$rank[ok], deriv = 0)$y
  df$d1[ok] <- predict(ss, x = df$rank[ok], deriv = 1)$y
  df$d2[ok] <- predict(ss, x = df$rank[ok], deriv = 2)$y

  zero_idx <- find_d2_zero_crossings(df$d2, zero_tol = zero_tol)
  zero_tbl <- summarize_zero_crossings(df, zero_idx)

  total_features <- nrow(df)
  ref_rank <- max(1L, total_features - leading_edge_fixed_k + 1L)
  ref_index <- which.min(abs(df$rank - ref_rank))

  left_zeros <- zero_tbl[zero_tbl$rank < ref_rank, , drop = FALSE]
  right_zeros <- zero_tbl[zero_tbl$rank > ref_rank, , drop = FALSE]

  if (nrow(left_zeros) == 0L) {
    stop("No d2 zero found immediately LEFT of the fixed leading-edge 5000 mark.", call. = FALSE)
  }
  if (nrow(right_zeros) == 0L) {
    stop("No d2 zero found immediately RIGHT of the fixed leading-edge 5000 mark.", call. = FALSE)
  }

  cutoff_row <- left_zeros[which.max(left_zeros$rank), , drop = FALSE]
  terminal_row <- right_zeros[which.min(right_zeros$rank), , drop = FALSE]

  cutoff_index <- cutoff_row$index
  terminal_index <- terminal_row$index

  right_idx <- seq.int(cutoff_index, total_features)
  m_right <- length(right_idx)

  left_end <- cutoff_index - 1L
  left_start <- max(1L, left_end - m_right + 1L)
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

  nb_threshold <- safe_median(df$nb_support[right_idx])

  selected_points <- data.frame(
    role = c("cutoff_anchor", "fixed_leadingedge_5000", "terminal_start"),
    index = c(cutoff_index, ref_index, terminal_index),
    rank = c(df$rank[cutoff_index], df$rank[ref_index], df$rank[terminal_index]),
    stringsAsFactors = FALSE
  )

  list(
    df = df,
    zero_tbl = zero_tbl,
    selected_points = selected_points,

    cutoff_anchor_index = cutoff_index,
    cutoff_anchor_rank = df$rank[cutoff_index],

    ref5000_index = ref_index,
    ref5000_rank = df$rank[ref_index],

    terminal_start_index = terminal_index,
    terminal_start_rank = df$rank[terminal_index],

    interval_min_rank = df$rank[cutoff_index],
    interval_max_rank = df$rank[terminal_index],

    pre_evs_remainder_size = cutoff_index - 1L,
    pre_evs_leading_edge_size = total_features - cutoff_index + 1L,

    total_features = total_features,
    left_idx = left_idx,
    right_idx = right_idx,
    corrob = corrob,
    nb_threshold = nb_threshold
  )
}

# =============================================================================
# FIGURE HELPERS
# =============================================================================

make_event_df <- function(sel) {
  data.frame(
    event = factor(
      c("Cutoff anchor", "Fixed leading-edge 5000", "Terminal start"),
      levels = c("Cutoff anchor", "Fixed leading-edge 5000", "Terminal start")
    ),
    xint = c(sel$cutoff_anchor_rank, sel$ref5000_rank, sel$terminal_start_rank),
    lty = c("dashed", "dotdash", "dotted"),
    col = c(COL["cutoff"], COL["ref5000"], COL["terminal"]),
    stringsAsFactors = FALSE
  )
}

make_marker_df <- function(sel, y_cutoff, y_ref, y_term) {
  data.frame(
    event = factor(
      c("Cutoff anchor", "Fixed leading-edge 5000", "Terminal start"),
      levels = c("Cutoff anchor", "Fixed leading-edge 5000", "Terminal start")
    ),
    x = c(sel$cutoff_anchor_rank, sel$ref5000_rank, sel$terminal_start_rank),
    y = c(y_cutoff, y_ref, y_term),
    shape_label = c("filled", "diamond", "open"),
    stringsAsFactors = FALSE
  )
}

make_legend_scales <- function() {
  list(
    scale_color_manual(
      values = c(
        "Cutoff anchor" = COL["cutoff"],
        "Fixed leading-edge 5000" = COL["ref5000"],
        "Terminal start" = COL["terminal"]
      )
    ),
    scale_linetype_manual(
      values = c(
        "Cutoff anchor" = "dashed",
        "Fixed leading-edge 5000" = "dotdash",
        "Terminal start" = "dotted"
      )
    ),
    scale_shape_manual(
      values = c(
        "Cutoff anchor" = 16,
        "Fixed leading-edge 5000" = 18,
        "Terminal start" = 1
      )
    )
  )
}

label_df <- function(x, y, text) {
  data.frame(x = x, y = y, label = text, stringsAsFactors = FALSE)
}

build_rank_panel <- function(sel, title_prefix, out_file) {
  df <- sel$df
  events <- make_event_df(sel)

  abs_markers <- make_marker_df(
    sel,
    y_cutoff = df$abs_loading[df$rank == sel$cutoff_anchor_rank],
    y_ref = df$abs_loading[df$rank == sel$ref5000_rank],
    y_term = df$abs_loading[df$rank == sel$terminal_start_rank]
  )

  var_markers <- make_marker_df(
    sel,
    y_cutoff = df$var_fit[df$rank == sel$cutoff_anchor_rank],
    y_ref = df$var_fit[df$rank == sel$ref5000_rank],
    y_term = df$var_fit[df$rank == sel$terminal_start_rank]
  )

  d2_markers <- make_marker_df(
    sel,
    y_cutoff = df$d2[df$rank == sel$cutoff_anchor_rank],
    y_ref = df$d2[df$rank == sel$ref5000_rank],
    y_term = df$d2[df$rank == sel$terminal_start_rank]
  )

  nb_markers <- make_marker_df(
    sel,
    y_cutoff = df$nb_support[df$rank == sel$cutoff_anchor_rank],
    y_ref = df$nb_support[df$rank == sel$ref5000_rank],
    y_term = df$nb_support[df$rank == sel$terminal_start_rank]
  )

  x_rng <- range(df$rank, na.rm = TRUE)
  x_left <- x_rng[1L] + 0.05 * diff(x_rng)
  x_right <- x_rng[1L] + 0.72 * diff(x_rng)

  abs_ymax <- max(df$abs_loading, na.rm = TRUE)
  var_ymax <- max(df$var_fit, na.rm = TRUE)
  d_ymax <- max(c(df$d1, df$d2), na.rm = TRUE)
  nb_ymax <- max(c(df$nb_support, df$log_nb1, df$log_nb2, df$log_alpha_mu, df$nb_gap), na.rm = TRUE)

  left_abs_label <- label_df(
    x_left, abs_ymax,
    paste(
      "Absolute loading panel",
      "Leading edge is on the RIGHT",
      "The grey band is the final geometric interval",
      "The legend gives the exact event colors and line types",
      sep = "\n"
    )
  )

  right_abs_label <- label_df(
    x_right, abs_ymax,
    paste0(
      "Cutoff anchor rank = ", sel$cutoff_anchor_rank,
      "\nReference rank (5000 from right) = ", sel$ref5000_rank,
      "\nTerminal start rank = ", sel$terminal_start_rank,
      "\nFinal interval = [", sel$interval_min_rank, ", ", sel$interval_max_rank, "]",
      "\nPre-EVS remainder = ", sel$pre_evs_remainder_size,
      "\nPre-EVS leading edge = ", sel$pre_evs_leading_edge_size
    )
  )

  left_var_label <- label_df(
    x_left, var_ymax,
    paste(
      "Variance curve method",
      "A smoothing spline is fit to empirical log(1+variance).",
      "Cutoff anchor = nearest d2 zero immediately LEFT of the fixed leading-edge 5000 mark.",
      "Terminal start = nearest d2 zero immediately RIGHT of the fixed leading-edge 5000 mark.",
      sep = "\n"
    )
  )

  left_deriv_label <- label_df(
    x_left, d_ymax,
    paste(
      "Derivative method",
      "All d2 zeros come only from sign changes or exact zeros.",
      "No amplitude threshold is used.",
      "No run-length filter is used.",
      sep = "\n"
    )
  )

  left_nb_label <- label_df(
    x_left, nb_ymax,
    paste(
      "NB corroboration",
      "The full RIGHT leading edge from cutoff anchor to the end is the NB2 side.",
      "The matched LEFT region contains the same number of genes immediately left of the cutoff anchor.",
      "These NB values corroborate but do not define the geometric split.",
      sep = "\n"
    )
  )

  right_nb_label <- label_df(
    x_right, nb_ymax,
    paste0(
      "Left median log(NB2) = ", fmt_num(sel$corrob$left_median_log_nb2),
      "\nRight median log(NB2) = ", fmt_num(sel$corrob$right_median_log_nb2),
      "\nRight-left log(NB2) diff = ", fmt_num(sel$corrob$right_left_log_nb2_diff),
      "\nLeft median NB2-NB1 = ", fmt_num(sel$corrob$left_median_nb_gap),
      "\nRight median NB2-NB1 = ", fmt_num(sel$corrob$right_median_nb_gap),
      "\nRight-left NB2-NB1 diff = ", fmt_num(sel$corrob$right_left_nb_gap_diff),
      "\nLeft median log(alpha*mu) = ", fmt_num(sel$corrob$left_median_log_alpha_mu),
      "\nRight median log(alpha*mu) = ", fmt_num(sel$corrob$right_median_log_alpha_mu),
      "\nRight-left log(alpha*mu) diff = ", fmt_num(sel$corrob$right_left_log_alpha_mu_diff)
    )
  )

  right_summary_label <- label_df(
    x_right,
    max(df$nb_support, na.rm = TRUE),
    paste0(
      "NB-supported summary",
      "\nAnchor = ", sel$cutoff_anchor_rank,
      "\nReference = ", sel$ref5000_rank,
      "\nTerminal = ", sel$terminal_start_rank,
      "\nInterval = [", sel$interval_min_rank, ", ", sel$interval_max_rank, "]"
    )
  )

  event_scales <- make_legend_scales()

  p1 <- ggplot(df, aes(x = rank, y = abs_loading)) +
    annotate(
      "rect",
      xmin = sel$interval_min_rank,
      xmax = sel$interval_max_rank,
      ymin = -Inf, ymax = Inf,
      fill = COL["interval"],
      alpha = 0.18
    ) +
    geom_line(color = COL["abs_loading"], linewidth = 1.0, na.rm = TRUE) +
    geom_vline(
      data = events,
      aes(xintercept = xint, color = event, linetype = event),
      linewidth = 0.9,
      show.legend = TRUE
    ) +
    geom_point(
      data = abs_markers,
      aes(x = x, y = y, color = event, shape = event),
      size = 2.8,
      stroke = 1.0,
      show.legend = TRUE
    ) +
    event_scales[[1L]] + event_scales[[2L]] + event_scales[[3L]] +
    geom_label(
      data = left_abs_label,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 2.55,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0.10, "lines"),
      linewidth = 0.25,
      fill = alpha("white", 0.96)
    ) +
    geom_label(
      data = right_abs_label,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 2.55,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0.10, "lines"),
      linewidth = 0.25,
      fill = alpha("white", 0.96)
    ) +
    labs(
      title = paste0(title_prefix, ": absolute PC1 loading series"),
      subtitle = "Leading edge is on the RIGHT",
      x = "EVS rank",
      y = "|PC1 loading|",
      color = "Event",
      linetype = "Event",
      shape = "Event"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(8, 40, 8, 12)
    )

  p2 <- ggplot(df, aes(x = rank, y = var_fit)) +
    annotate(
      "rect",
      xmin = sel$interval_min_rank,
      xmax = sel$interval_max_rank,
      ymin = -Inf, ymax = Inf,
      fill = COL["interval"],
      alpha = 0.18
    ) +
    geom_line(color = COL["variance_fit"], linewidth = 1.0, na.rm = TRUE) +
    geom_vline(
      data = events,
      aes(xintercept = xint, color = event, linetype = event),
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    geom_point(
      data = var_markers,
      aes(x = x, y = y, color = event, shape = event),
      size = 2.8,
      stroke = 1.0,
      show.legend = FALSE
    ) +
    event_scales[[1L]] + event_scales[[2L]] + event_scales[[3L]] +
    geom_label(
      data = left_var_label,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 2.25,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0.10, "lines"),
      linewidth = 0.25,
      fill = alpha("white", 0.96)
    ) +
    labs(
      title = paste0(title_prefix, ": smoothed empirical variance curve"),
      subtitle = "The geometric points are marked directly on the fitted curve",
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      plot.margin = margin(8, 40, 8, 12)
    )

  deriv_df <- bind_rows(
    data.frame(rank = df$rank, value = df$d2, series = "Second derivative d2", stringsAsFactors = FALSE),
    data.frame(rank = df$rank, value = df$d1, series = "First derivative d1", stringsAsFactors = FALSE)
  )
  deriv_df$series <- factor(
    deriv_df$series,
    levels = c("Second derivative d2", "First derivative d1")
  )

  p3 <- ggplot(deriv_df, aes(x = rank, y = value, color = series)) +
    annotate(
      "rect",
      xmin = sel$interval_min_rank,
      xmax = sel$interval_max_rank,
      ymin = -Inf, ymax = Inf,
      fill = COL["interval"],
      alpha = 0.18
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.55) +
    geom_line(linewidth = 0.95, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "Second derivative d2" = COL["d2"],
        "First derivative d1" = COL["d1"]
      )
    ) +
    geom_vline(
      data = events,
      aes(xintercept = xint, linetype = event),
      color = c(COL["cutoff"], COL["ref5000"], COL["terminal"]),
      linewidth = 0.9,
      show.legend = TRUE
    ) +
    geom_point(
      data = d2_markers,
      aes(x = x, y = y, shape = event),
      color = c(COL["cutoff"], COL["ref5000"], COL["terminal"]),
      size = 2.8,
      stroke = 1.0,
      show.legend = TRUE
    ) +
    scale_linetype_manual(
      values = c(
        "Cutoff anchor" = "dashed",
        "Fixed leading-edge 5000" = "dotdash",
        "Terminal start" = "dotted"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Cutoff anchor" = 16,
        "Fixed leading-edge 5000" = 18,
        "Terminal start" = 1
      )
    ) +
    geom_label(
      data = left_deriv_label,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 2.25,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0.10, "lines"),
      linewidth = 0.25,
      fill = alpha("white", 0.96)
    ) +
    labs(
      title = paste0(title_prefix, ": derivative support"),
      subtitle = "The selected d2 zero-crossings around the fixed leading-edge 5000 mark define the geometric interval",
      x = "EVS rank",
      y = "Derivative value",
      color = "Derivative",
      linetype = "Event",
      shape = "Event"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(8, 40, 8, 12)
    )

  nb_long <- bind_rows(
    data.frame(rank = df$rank, value = df$nb_support, series = "Combined NB support", stringsAsFactors = FALSE),
    data.frame(rank = df$rank, value = df$log_nb1, series = "NB1 = log(1+mu)", stringsAsFactors = FALSE),
    data.frame(rank = df$rank, value = df$log_nb2, series = "NB2 = log(1+variance-mu)", stringsAsFactors = FALSE),
    data.frame(rank = df$rank, value = df$log_alpha_mu, series = "log(alpha*mu)", stringsAsFactors = FALSE),
    data.frame(rank = df$rank, value = df$nb_gap, series = "log(NB2+1)-log(NB1+1)", stringsAsFactors = FALSE)
  )
  nb_long$series <- factor(
    nb_long$series,
    levels = c(
      "Combined NB support",
      "NB1 = log(1+mu)",
      "NB2 = log(1+variance-mu)",
      "log(alpha*mu)",
      "log(NB2+1)-log(NB1+1)"
    )
  )

  p4 <- ggplot(nb_long, aes(x = rank, y = value, color = series)) +
    annotate(
      "rect",
      xmin = sel$interval_min_rank,
      xmax = sel$interval_max_rank,
      ymin = -Inf, ymax = Inf,
      fill = COL["interval"],
      alpha = 0.18
    ) +
    geom_hline(yintercept = sel$nb_threshold, color = "grey40", linetype = 3, linewidth = 0.6) +
    geom_line(linewidth = 0.90, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "Combined NB support" = COL["nb_support"],
        "NB1 = log(1+mu)" = COL["nb1"],
        "NB2 = log(1+variance-mu)" = COL["nb2"],
        "log(alpha*mu)" = COL["amu"],
        "log(NB2+1)-log(NB1+1)" = COL["nb_gap"]
      )
    ) +
    geom_vline(
      data = events,
      aes(xintercept = xint, linetype = event),
      color = c(COL["cutoff"], COL["ref5000"], COL["terminal"]),
      linewidth = 0.9,
      show.legend = TRUE
    ) +
    geom_point(
      data = nb_markers,
      aes(x = x, y = y, shape = event),
      color = c(COL["cutoff"], COL["ref5000"], COL["terminal"]),
      size = 2.8,
      stroke = 1.0,
      show.legend = TRUE
    ) +
    scale_linetype_manual(
      values = c(
        "Cutoff anchor" = "dashed",
        "Fixed leading-edge 5000" = "dotdash",
        "Terminal start" = "dotted"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Cutoff anchor" = 16,
        "Fixed leading-edge 5000" = 18,
        "Terminal start" = 1
      )
    ) +
    geom_label(
      data = left_nb_label,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 2.05,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0.10, "lines"),
      linewidth = 0.25,
      fill = alpha("white", 0.96)
    ) +
    geom_label(
      data = right_nb_label,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 1.90,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0.10, "lines"),
      linewidth = 0.25,
      fill = alpha("white", 0.96)
    ) +
    labs(
      title = paste0(title_prefix, ": NB1 / NB2 / alpha*mu support"),
      subtitle = "The full right-of-cutoff leading edge is compared against a matched left region of equal size",
      x = "EVS rank",
      y = "Support value",
      color = "NB quantity",
      linetype = "Event",
      shape = "Event"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(8, 40, 8, 12)
    )

  p5 <- ggplot(df, aes(x = rank, y = nb_support)) +
    annotate(
      "rect",
      xmin = sel$interval_min_rank,
      xmax = sel$interval_max_rank,
      ymin = -Inf, ymax = Inf,
      fill = COL["interval"],
      alpha = 0.18
    ) +
    geom_hline(yintercept = sel$nb_threshold, color = "grey40", linetype = 3, linewidth = 0.6) +
    geom_line(color = COL["nb_support"], linewidth = 1.0, na.rm = TRUE) +
    geom_vline(
      data = events,
      aes(xintercept = xint, color = event, linetype = event),
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    geom_point(
      data = nb_markers,
      aes(x = x, y = y, color = event, shape = event),
      size = 2.8,
      stroke = 1.0,
      show.legend = FALSE
    ) +
    event_scales[[1L]] + event_scales[[2L]] + event_scales[[3L]] +
    geom_label(
      data = right_summary_label,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 2.10,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0.10, "lines"),
      linewidth = 0.25,
      fill = alpha("white", 0.96)
    ) +
    labs(
      title = paste0(title_prefix, ": NB-supported summary"),
      subtitle = "The geometric interval remains primary; NB support remains corroborative",
      x = "EVS rank",
      y = "Combined NB support"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      plot.margin = margin(8, 40, 8, 12)
    )

  png(out_file, width = 3200, height = 4300, res = 260)
  grid.arrange(p1, p2, p3, p4, p5, ncol = 1)
  dev.off()
}

build_interval_overlay <- function(ctrl_sel, trt_sel, cmp_name, out_file) {
  ctrl_df <- ctrl_sel$df
  trt_df <- trt_sel$df

  overlay_df <- bind_rows(
    data.frame(rank = ctrl_df$rank, value = ctrl_df$var_fit, arm = "Control variance fit", stringsAsFactors = FALSE),
    data.frame(rank = trt_df$rank, value = trt_df$var_fit, arm = "Treatment variance fit", stringsAsFactors = FALSE)
  )
  overlay_df$arm <- factor(
    overlay_df$arm,
    levels = c("Control variance fit", "Treatment variance fit")
  )

  x_rng <- range(overlay_df$rank, na.rm = TRUE)
  y_top <- max(overlay_df$value, na.rm = TRUE)

  left_text <- label_df(
    x_rng[1L] + 0.05 * diff(x_rng),
    y_top,
    paste(
      "Overlay panel",
      "Control and treatment intervals are shown on the same fitted-variance scale.",
      "Each arm uses the same geometric rule around the fixed leading-edge 5000 reference.",
      sep = "\n"
    )
  )

  right_text <- label_df(
    x_rng[1L] + 0.68 * diff(x_rng),
    y_top,
    paste0(
      "Control interval = [", ctrl_sel$interval_min_rank, ", ", ctrl_sel$interval_max_rank, "]",
      "\nTreatment interval = [", trt_sel$interval_min_rank, ", ", trt_sel$interval_max_rank, "]",
      "\nControl reference = ", ctrl_sel$ref5000_rank,
      "\nTreatment reference = ", trt_sel$ref5000_rank
    )
  )

  p <- ggplot(overlay_df, aes(x = rank, y = value, color = arm)) +
    annotate(
      "rect",
      xmin = min(ctrl_sel$interval_min_rank, trt_sel$interval_min_rank),
      xmax = max(ctrl_sel$interval_max_rank, trt_sel$interval_max_rank),
      ymin = -Inf, ymax = Inf,
      fill = COL["interval"],
      alpha = 0.12
    ) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "Control variance fit" = "#1B9E77",
        "Treatment variance fit" = "#7570B3"
      )
    ) +
    geom_vline(xintercept = ctrl_sel$cutoff_anchor_rank, linetype = 2, linewidth = 0.85, color = "#252525") +
    geom_vline(xintercept = trt_sel$cutoff_anchor_rank, linetype = 3, linewidth = 0.85, color = "#636363") +
    geom_vline(xintercept = ctrl_sel$ref5000_rank, linetype = 4, linewidth = 0.75, color = COL["ref5000"]) +
    geom_vline(xintercept = trt_sel$ref5000_rank, linetype = 4, linewidth = 0.75, color = COL["terminal"]) +
    geom_label(
      data = left_text,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 2.55,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0.10, "lines"),
      linewidth = 0.25,
      fill = alpha("white", 0.96)
    ) +
    geom_label(
      data = right_text,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 2.55,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0.10, "lines"),
      linewidth = 0.25,
      fill = alpha("white", 0.96)
    ) +
    labs(
      title = paste0(cmp_name, ": treatment/control variance interval overlay"),
      subtitle = "Control and treatment fitted variance curves are displayed on the same axis",
      x = "EVS rank",
      y = "Fitted log(1 + variance)",
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(8, 40, 8, 12)
    )

  png(out_file, width = 3200, height = 1500, res = 260)
  print(p)
  dev.off()
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

message("SEQUENCE.R started"); flush.console()
message("Resolving count file..."); flush.console()
count_file <- resolve_counts_file(count_file_hint)
message("Using count file: ", count_file); flush.console()

message("Reading count matrix..."); flush.console()
loaded <- read_count_matrix(count_file, meta_all$id)
count_matrix <- loaded$count_matrix
annot_df <- loaded$annotation
message("Count matrix dimensions: ", nrow(count_matrix), " features x ", ncol(count_matrix), " samples"); flush.console()

overall_rows <- list()

for (i in seq_len(nrow(comparison_table))) {
  row_i <- comparison_table[i, , drop = FALSE]
  cmp_name <- row_i$comparison_name[[1]]

  message("Processing comparison: ", cmp_name); flush.console()

  cmp_dir <- file.path(output_root, paste0(cmp_name, "_cutoff_folder"))
  dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

  cmp_counts <- subset_comparison_counts(count_matrix, row_i, meta_all)
  feature_tbl <- build_feature_table(cmp_counts$counts, cmp_counts$ctrl_ids, cmp_counts$trt_ids, annot_df)

  utils::write.csv(
    feature_tbl,
    file.path(cmp_dir, paste0(cmp_name, "_feature_level_metrics.csv")),
    row.names = FALSE
  )

  ctrl_df <- build_rank_df(feature_tbl, "control")
  trt_df  <- build_rank_df(feature_tbl, "treatment")

  utils::write.csv(ctrl_df, file.path(cmp_dir, paste0(cmp_name, "_control_rank_series.csv")), row.names = FALSE)
  utils::write.csv(trt_df, file.path(cmp_dir, paste0(cmp_name, "_treatment_rank_series.csv")), row.names = FALSE)

  message("Selecting control cutoff for ", cmp_name); flush.console()
  ctrl_sel <- select_interval_and_nb(
    ctrl_df,
    spline_spar = spline_spar,
    zero_tol = zero_tol,
    leading_edge_fixed_k = leading_edge_fixed_k
  )
  message(
    "Control cutoff anchor rank: ", ctrl_sel$cutoff_anchor_rank,
    " | reference rank: ", ctrl_sel$ref5000_rank,
    " | terminal anchor: ", ctrl_sel$terminal_start_rank,
    " | interval: [", ctrl_sel$interval_min_rank, ", ", ctrl_sel$interval_max_rank, "]",
    " | pre-EVS remainder: ", ctrl_sel$pre_evs_remainder_size,
    " | pre-EVS leading edge: ", ctrl_sel$pre_evs_leading_edge_size
  ); flush.console()

  message("Selecting treatment cutoff for ", cmp_name); flush.console()
  trt_sel <- select_interval_and_nb(
    trt_df,
    spline_spar = spline_spar,
    zero_tol = zero_tol,
    leading_edge_fixed_k = leading_edge_fixed_k
  )
  message(
    "Treatment cutoff anchor rank: ", trt_sel$cutoff_anchor_rank,
    " | reference rank: ", trt_sel$ref5000_rank,
    " | terminal anchor: ", trt_sel$terminal_start_rank,
    " | interval: [", trt_sel$interval_min_rank, ", ", trt_sel$interval_max_rank, "]",
    " | pre-EVS remainder: ", trt_sel$pre_evs_remainder_size,
    " | pre-EVS leading edge: ", trt_sel$pre_evs_leading_edge_size
  ); flush.console()

  build_rank_panel(
    ctrl_sel,
    paste0(cmp_name, " control"),
    file.path(cmp_dir, paste0(cmp_name, "_control_rank_panel.png"))
  )

  build_rank_panel(
    trt_sel,
    paste0(cmp_name, " treatment"),
    file.path(cmp_dir, paste0(cmp_name, "_treatment_rank_panel.png"))
  )

  build_interval_overlay(
    ctrl_sel,
    trt_sel,
    cmp_name,
    file.path(cmp_dir, paste0(cmp_name, "_cutoff_range_panel.png"))
  )

  utils::write.csv(
    ctrl_sel$zero_tbl,
    file.path(cmp_dir, paste0(cmp_name, "_control_zero_crossings_all.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    trt_sel$zero_tbl,
    file.path(cmp_dir, paste0(cmp_name, "_treatment_zero_crossings_all.csv")),
    row.names = FALSE
  )

  utils::write.csv(
    ctrl_sel$selected_points,
    file.path(cmp_dir, paste0(cmp_name, "_control_selected_points.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    trt_sel$selected_points,
    file.path(cmp_dir, paste0(cmp_name, "_treatment_selected_points.csv")),
    row.names = FALSE
  )

  overall_rows[[cmp_name]] <- tibble(
    comparison = cmp_name,
    total_features = ctrl_sel$total_features,
    leading_edge_fixed_k = leading_edge_fixed_k,

    control_cutoff_anchor_rank = ctrl_sel$cutoff_anchor_rank,
    control_reference_rank_5000_from_right = ctrl_sel$ref5000_rank,
    control_terminal_start_rank = ctrl_sel$terminal_start_rank,
    control_interval_rank_min = ctrl_sel$interval_min_rank,
    control_interval_rank_max = ctrl_sel$interval_max_rank,
    control_pre_evs_remainder_size = ctrl_sel$pre_evs_remainder_size,
    control_pre_evs_leading_edge_size = ctrl_sel$pre_evs_leading_edge_size,
    control_left_n = ctrl_sel$corrob$left_n,
    control_right_n = ctrl_sel$corrob$right_n,
    control_left_median_log_nb2 = ctrl_sel$corrob$left_median_log_nb2,
    control_right_median_log_nb2 = ctrl_sel$corrob$right_median_log_nb2,
    control_right_left_log_nb2_diff = ctrl_sel$corrob$right_left_log_nb2_diff,
    control_left_median_nb_gap = ctrl_sel$corrob$left_median_nb_gap,
    control_right_median_nb_gap = ctrl_sel$corrob$right_median_nb_gap,
    control_right_left_nb_gap_diff = ctrl_sel$corrob$right_left_nb_gap_diff,
    control_left_median_log_alpha_mu = ctrl_sel$corrob$left_median_log_alpha_mu,
    control_right_median_log_alpha_mu = ctrl_sel$corrob$right_median_log_alpha_mu,
    control_right_left_log_alpha_mu_diff = ctrl_sel$corrob$right_left_log_alpha_mu_diff,
    control_nb_support_threshold = ctrl_sel$nb_threshold,

    treatment_cutoff_anchor_rank = trt_sel$cutoff_anchor_rank,
    treatment_reference_rank_5000_from_right = trt_sel$ref5000_rank,
    treatment_terminal_start_rank = trt_sel$terminal_start_rank,
    treatment_interval_rank_min = trt_sel$interval_min_rank,
    treatment_interval_rank_max = trt_sel$interval_max_rank,
    treatment_pre_evs_remainder_size = trt_sel$pre_evs_remainder_size,
    treatment_pre_evs_leading_edge_size = trt_sel$pre_evs_leading_edge_size,
    treatment_left_n = trt_sel$corrob$left_n,
    treatment_right_n = trt_sel$corrob$right_n,
    treatment_left_median_log_nb2 = trt_sel$corrob$left_median_log_nb2,
    treatment_right_median_log_nb2 = trt_sel$corrob$right_median_log_nb2,
    treatment_right_left_log_nb2_diff = trt_sel$corrob$right_left_log_nb2_diff,
    treatment_left_median_nb_gap = trt_sel$corrob$left_median_nb_gap,
    treatment_right_median_nb_gap = trt_sel$corrob$right_median_nb_gap,
    treatment_right_left_nb_gap_diff = trt_sel$corrob$right_left_nb_gap_diff,
    treatment_left_median_log_alpha_mu = trt_sel$corrob$left_median_log_alpha_mu,
    treatment_right_median_log_alpha_mu = trt_sel$corrob$right_median_log_alpha_mu,
    treatment_right_left_log_alpha_mu_diff = trt_sel$corrob$right_left_log_alpha_mu_diff,
    treatment_nb_support_threshold = trt_sel$nb_threshold
  )

  utils::write.csv(
    overall_rows[[cmp_name]],
    file.path(cmp_dir, paste0(cmp_name, "_cutoff_summary.csv")),
    row.names = FALSE
  )
}

overall_summary <- bind_rows(overall_rows)

utils::write.csv(
  overall_summary,
  file.path(output_root, "overall_cutoff_summary.csv"),
  row.names = FALSE
)

message("Done. Outputs written to: ", output_root); flush.console()
