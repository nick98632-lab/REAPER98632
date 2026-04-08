# =============================================================================
# SEQUENCE PRE-EVS LEADING-EDGE CUTOFF SELECTOR
# FINAL MANUSCRIPT VERSION
# =============================================================================
#
# METHODS SUMMARY
# The ranked axis is defined from LEFT to RIGHT by increasing absolute PC1
# loading. Therefore the leading edge is located on the RIGHT side of the axis.
#
# For each arm separately (control and treatment), the cutoff selection is
# defined by the following fixed manuscript rule:
#
# 1. Compute the smoothed empirical variance curve over EVS rank.
# 2. Compute the first and second derivatives of that same smooth.
# 3. Define the fixed leading-edge reference mark as 5000 genes from the
#    RIGHT edge, not rank 5000 from the left:
#         leading_edge_5000_rank = total_features - 5000 + 1
# 4. Identify all second-derivative zero-crossings.
# 5. Define the cutoff anchor as the nearest zero-crossing immediately to the
#    LEFT of the fixed leading-edge 5000 mark.
# 6. Define the terminal start as the nearest zero-crossing immediately to the
#    RIGHT of the fixed leading-edge 5000 mark.
# 7. The study interval is therefore:
#         [cutoff anchor, terminal start]
#    so that the fixed leading-edge 5000 mark lies inside the interval.
#
# NB corroboration is then computed relative to the cutoff anchor:
#
# 1. The NB2-like right-side region is ALL genes from cutoff anchor to the end
#    of the ranked series.
# 2. Let that right-side region contain m genes.
# 3. The matched left-side region is the same number m of genes immediately to
#    the left of the cutoff anchor.
# 4. Median NB2, NB2-NB1 contrast, and log(alpha*mu) are compared between the
#    right and matched-left regions and reported explicitly on the figure.
#
# IMPORTANT IMPLEMENTATION NOTES
# - No threshold is used to decide whether a zero-crossing "counts".
# - No minimum run length is used.
# - No sign-direction constraint is imposed on d2 before or after the zero.
# - The only geometric anchors are the nearest d2 zero immediately left of the
#   fixed leading-edge 5000 mark and the nearest d2 zero immediately right of
#   that same mark.
# - All figure annotations use geom_label(), not annotate("label", ...), so the
#   previous label.size warnings are removed.
# - All color legends use exact factor labels matching the plotted data.
# - Log transforms are guarded so non-positive values produce NA rather than
#   warnings.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(gridExtra)
  library(scales)
})

options(warn = 1)

# =============================================================================
# USER SETTINGS
# =============================================================================

repo_dir <- getwd()
count_file_hint <- file.path(repo_dir, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")

out_root <- file.path(
  repo_dir,
  "exports",
  "variance_derivative_nb_range_leadingedge5000_final"
)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

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

spline_spar <- 0.60
nb_smooth_window <- 151L
leading_edge_fixed_k <- 5000L
zero_tol <- 1e-10

# =============================================================================
# COLOR PALETTE
# =============================================================================
#
# METHODS SUMMARY
# Colors are fixed and reused identically in all panels so that the same
# quantity is always represented by the same color.
# =============================================================================

COLORS <- c(
  abs_loading     = "#111111",
  variance_fit    = "#111111",
  d1              = "#00B4D8",
  d2              = "#D1495B",
  nb_support      = "#111111",
  nb1             = "#D8B365",
  nb2             = "#66C2A5",
  amu             = "#1F78FF",
  nb_gap          = "#C77CFF",
  cutoff_anchor   = "#111111",
  terminal_start  = "#111111",
  ref5000         = "#8C510A",
  interval_fill   = "#BDBDBD"
)

# =============================================================================
# UTILITY HELPERS
# =============================================================================

resolve_counts_file <- function(path_hint) {
  candidates <- unique(c(
    path_hint,
    file.path(getwd(), path_hint),
    file.path(getwd(), "data", basename(path_hint)),
    "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
    "/root/REAPER98632/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
  ))

  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0L) return(existing[[1]])

  found <- list.files(
    path = "/root/REAPER98632",
    pattern = "WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
    recursive = TRUE,
    full.names = TRUE
  )
  found <- found[file.exists(found)]
  if (length(found) > 0L) return(found[[1]])

  stop(
    paste0("Count file not found. Tried: ", paste(candidates, collapse = ", ")),
    call. = FALSE
  )
}

detect_feature_id_column <- function(df) {
  candidates <- c("OrigID", "feature_id", "FeatureID", "PAS", "pas_id", "GeneID", "gene_id", "id")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0L) return(hit[1])
  names(df)[1]
}

detect_gene_symbol_column <- function(df) {
  candidates <- c("Symbol", "symbol", "GeneSymbol", "gene_symbol", "Gene", "gene")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0L) return(hit[1])
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

safe_ratio <- function(num, den) {
  if (!is.finite(num) || !is.finite(den) || den == 0) return(NA_real_)
  num / den
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

fmt_num <- function(x, digits = 3) {
  if (!is.finite(x)) return("NA")
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

scale01 <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)

  rng <- range(x[ok], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || abs(rng[2] - rng[1]) < .Machine$double.eps) {
    out[ok] <- 0
    return(out)
  }

  out[ok] <- (x[ok] - rng[1]) / (rng[2] - rng[1])
  out
}

roll_median <- function(x, k = 101L) {
  x <- as.numeric(x)
  n <- length(x)
  if (k < 1L) return(x)
  if (k %% 2L == 0L) k <- k + 1L
  h <- (k - 1L) / 2L
  out <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    lo <- max(1L, i - h)
    hi <- min(n, i + h)
    vals <- x[lo:hi]
    vals <- vals[is.finite(vals)]
    out[i] <- if (length(vals)) median(vals) else NA_real_
  }

  out
}

make_label_df <- function(x, y, label) {
  data.frame(x = x, y = y, label = label, stringsAsFactors = FALSE)
}

# =============================================================================
# DATA LOADING
# =============================================================================
#
# METHODS SUMMARY
# Raw counts are read once. Feature identifiers are retained. If duplicate
# feature identifiers exist, they are collapsed by summation.
# =============================================================================

read_count_matrix <- function(path, meta_ids) {
  raw_df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  feature_col <- detect_feature_id_column(raw_df)
  gene_col <- detect_gene_symbol_column(raw_df)
  sample_cols <- intersect(meta_ids, names(raw_df))

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

  list(count_matrix = count_mat, annot_df = annot_df)
}

subset_comparison <- function(count_matrix, comparison_row, meta_all) {
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
    count_matrix = count_matrix[, keep_ids, drop = FALSE],
    trt_ids = trt_ids,
    ctrl_ids = ctrl_ids
  )
}

# =============================================================================
# EVS RANKING
# =============================================================================
#
# METHODS SUMMARY
# Absolute PC1 loadings are computed within each arm separately after log2(count
# + 1) transformation. EVS ranks then increase from left to right.
# =============================================================================

compute_group_abs_pc1_loadings <- function(count_mat, group_cols) {
  mat <- count_mat[, group_cols, drop = FALSE]
  mat <- log2(mat + 1)
  mat <- t(mat)

  if (nrow(mat) < 2L) return(rep(NA_real_, ncol(mat)))

  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  abs(pca$rotation[, 1L])
}

build_evs_table <- function(count_mat, ctrl_cols, trt_cols, feature_ids) {
  ctrl_load <- compute_group_abs_pc1_loadings(count_mat, ctrl_cols)
  trt_load  <- compute_group_abs_pc1_loadings(count_mat, trt_cols)

  tibble(
    feature_id = feature_ids,
    ctrl_abs_loading = ctrl_load,
    trt_abs_loading = trt_load
  ) %>%
    mutate(
      rank_ctrl = rank(ctrl_abs_loading, ties.method = "first"),
      rank_trt  = rank(trt_abs_loading, ties.method = "first")
    )
}

# =============================================================================
# EMPIRICAL NB QUANTITIES
# =============================================================================
#
# METHODS SUMMARY
# Empirical mean and variance are computed directly from raw counts. Empirical
# NB-derived quantities are defined as:
#   NB1 = mu
#   NB2 = variance - mu
#   alpha*mu = (variance - mu) / mu
# Only positive values are logged.
# =============================================================================

compute_feature_metrics_empirical <- function(count_mat, ctrl_cols, trt_cols, feature_ids) {
  ctrl_mu  <- apply(count_mat[, ctrl_cols, drop = FALSE], 1L, safe_mean)
  trt_mu   <- apply(count_mat[, trt_cols, drop = FALSE], 1L, safe_mean)
  ctrl_var <- apply(count_mat[, ctrl_cols, drop = FALSE], 1L, safe_var)
  trt_var  <- apply(count_mat[, trt_cols, drop = FALSE], 1L, safe_var)

  tibble(
    feature_id = feature_ids,
    mu_ctrl = ctrl_mu,
    mu_trt = trt_mu,
    var_ctrl = ctrl_var,
    var_trt = trt_var
  ) %>%
    mutate(
      nb1_ctrl = mu_ctrl,
      nb1_trt  = mu_trt,
      nb2_ctrl = pmax(var_ctrl - mu_ctrl, 0),
      nb2_trt  = pmax(var_trt - mu_trt, 0),
      alpha_mu_ctrl = ifelse(mu_ctrl > 0, (var_ctrl - mu_ctrl) / mu_ctrl, NA_real_),
      alpha_mu_trt  = ifelse(mu_trt  > 0, (var_trt  - mu_trt)  / mu_trt,  NA_real_)
    )
}

build_rank_series <- function(full_tbl, arm = c("control", "treatment")) {
  arm <- match.arg(arm)

  if (arm == "control") {
    out <- full_tbl %>%
      transmute(
        feature_id, gene_symbol,
        rank = rank_ctrl,
        abs_loading = ctrl_abs_loading,
        mu = mu_ctrl,
        variance = var_ctrl,
        nb1 = nb1_ctrl,
        nb2 = nb2_ctrl,
        alpha_mu = alpha_mu_ctrl
      ) %>%
      arrange(rank)
  } else {
    out <- full_tbl %>%
      transmute(
        feature_id, gene_symbol,
        rank = rank_trt,
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
  out$nb_gap <- out$log_nb2 - out$log_nb1
  out$log_alpha_mu <- safe_log_pos(out$alpha_mu)

  out
}

# =============================================================================
# DERIVATIVE GEOMETRY
# =============================================================================
#
# METHODS SUMMARY
# The variance curve is smoothed by smoothing spline. The first and second
# derivatives are evaluated on the full rank grid. Zero-crossings are defined
# only by second-derivative sign changes or exact zeros; no extra thresholding,
# run-length filter, or sign-direction filter is used.
# =============================================================================

find_d2_zero_crossings <- function(d2_vec, zero_tol = 1e-10) {
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

summarize_zero_crossings <- function(df, zero_idx) {
  if (!length(zero_idx)) {
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
    index = zero_idx,
    rank = df$rank[zero_idx],
    d1_here = df$d1_sm[zero_idx],
    d2_here = df$d2_sm[zero_idx],
    var_fit_here = df$var_fit[zero_idx],
    stringsAsFactors = FALSE
  )
}

select_cutoff_interval <- function(rank_df,
                                   spline_spar = 0.60,
                                   zero_tol = 1e-10,
                                   nb_smooth_window = 151L,
                                   leading_edge_fixed_k = 5000L) {
  df <- rank_df

  ok <- is.finite(df$rank) & is.finite(df$log_variance)
  if (sum(ok) < 10L) {
    stop("Not enough finite variance points for smoothing.", call. = FALSE)
  }

  sp_fit <- smooth.spline(x = df$rank[ok], y = df$log_variance[ok], spar = spline_spar)
  pred0 <- predict(sp_fit, x = df$rank[ok], deriv = 0)
  pred1 <- predict(sp_fit, x = df$rank[ok], deriv = 1)
  pred2 <- predict(sp_fit, x = df$rank[ok], deriv = 2)

  df$var_fit <- NA_real_
  df$d1_sm <- NA_real_
  df$d2_sm <- NA_real_

  df$var_fit[ok] <- pred0$y
  df$d1_sm[ok] <- pred1$y
  df$d2_sm[ok] <- pred2$y

  df$nb_gap_sm <- roll_median(df$nb_gap, nb_smooth_window)
  df$amu_sm <- roll_median(df$log_alpha_mu, nb_smooth_window)
  df$log_nb2_sm <- roll_median(df$log_nb2, nb_smooth_window)

  gap_s <- scale01(df$nb_gap_sm)
  amu_s <- scale01(df$amu_sm)
  nb2_s <- scale01(df$log_nb2_sm)

  gap_s[!is.finite(gap_s)] <- 0
  amu_s[!is.finite(amu_s)] <- 0
  nb2_s[!is.finite(nb2_s)] <- 0

  df$nb_support <- (gap_s + amu_s + nb2_s) / 3

  zero_idx <- find_d2_zero_crossings(df$d2_sm, zero_tol = zero_tol)
  zero_tbl <- summarize_zero_crossings(df, zero_idx)

  n_features <- nrow(df)
  ref_rank <- max(1L, n_features - leading_edge_fixed_k + 1L)
  ref_index <- which.min(abs(df$rank - ref_rank))

  left_candidates <- zero_tbl[zero_tbl$rank < ref_rank, , drop = FALSE]
  right_candidates <- zero_tbl[zero_tbl$rank > ref_rank, , drop = FALSE]

  if (nrow(left_candidates) == 0L) {
    stop("No second-derivative zero-crossing was found to the LEFT of the fixed leading-edge 5000 mark.", call. = FALSE)
  }
  if (nrow(right_candidates) == 0L) {
    stop("No second-derivative zero-crossing was found to the RIGHT of the fixed leading-edge 5000 mark.", call. = FALSE)
  }

  cutoff_row <- left_candidates[which.max(left_candidates$rank), , drop = FALSE]
  terminal_row <- right_candidates[which.min(right_candidates$rank), , drop = FALSE]

  cutoff_index <- cutoff_row$index
  terminal_index <- terminal_row$index

  right_idx <- seq.int(cutoff_index, n_features)
  m_right <- length(right_idx)

  left_end <- cutoff_index - 1L
  left_start <- max(1L, left_end - m_right + 1L)
  left_idx <- if (left_end >= left_start) seq.int(left_start, left_end) else integer(0)

  corrob <- data.frame(
    left_n = length(left_idx),
    right_n = length(right_idx),

    left_median_log_nb2 = if (length(left_idx)) safe_median(df$log_nb2[left_idx]) else NA_real_,
    right_median_log_nb2 = if (length(right_idx)) safe_median(df$log_nb2[right_idx]) else NA_real_,

    left_median_nb_gap = if (length(left_idx)) safe_median(df$nb_gap_sm[left_idx]) else NA_real_,
    right_median_nb_gap = if (length(right_idx)) safe_median(df$nb_gap_sm[right_idx]) else NA_real_,

    left_median_log_alpha_mu = if (length(left_idx)) safe_median(df$amu_sm[left_idx]) else NA_real_,
    right_median_log_alpha_mu = if (length(right_idx)) safe_median(df$amu_sm[right_idx]) else NA_real_,

    left_median_nb_support = if (length(left_idx)) safe_median(df$nb_support[left_idx]) else NA_real_,
    right_median_nb_support = if (length(right_idx)) safe_median(df$nb_support[right_idx]) else NA_real_,

    stringsAsFactors = FALSE
  )

  corrob$right_left_log_nb2_diff <- corrob$right_median_log_nb2 - corrob$left_median_log_nb2
  corrob$right_left_nb_gap_diff <- corrob$right_median_nb_gap - corrob$left_median_nb_gap
  corrob$right_left_log_alpha_mu_diff <- corrob$right_median_log_alpha_mu - corrob$left_median_log_alpha_mu
  corrob$right_left_nb_support_diff <- corrob$right_median_nb_support - corrob$left_median_nb_support

  nb_threshold <- safe_median(df$nb_support[right_idx])

  selected_tbl <- data.frame(
    role = c("cutoff_anchor", "leading_edge_5000_reference", "terminal_start"),
    index = c(cutoff_index, ref_index, terminal_index),
    rank = c(df$rank[cutoff_index], df$rank[ref_index], df$rank[terminal_index]),
    d1_here = c(df$d1_sm[cutoff_index], df$d1_sm[ref_index], df$d1_sm[terminal_index]),
    d2_here = c(df$d2_sm[cutoff_index], df$d2_sm[ref_index], df$d2_sm[terminal_index]),
    var_fit_here = c(df$var_fit[cutoff_index], df$var_fit[ref_index], df$var_fit[terminal_index]),
    stringsAsFactors = FALSE
  )

  list(
    curve_df = df,
    zero_crossings = zero_tbl,
    selected_points = selected_tbl,

    cutoff_center_index = cutoff_index,
    cutoff_center_rank = df$rank[cutoff_index],

    leading_edge_5000_index = ref_index,
    leading_edge_5000_rank = df$rank[ref_index],

    terminal_start_index = terminal_index,
    terminal_start_rank = df$rank[terminal_index],

    cutoff_range_index_min = cutoff_index,
    cutoff_range_index_max = terminal_index,
    cutoff_range_rank_min = df$rank[cutoff_index],
    cutoff_range_rank_max = df$rank[terminal_index],

    center_nb_support = df$nb_support[cutoff_index],
    nb_support_threshold = nb_threshold,

    total_features = n_features,
    pre_evs_remainder_size = cutoff_index - 1L,
    pre_evs_leading_edge_size = n_features - cutoff_index + 1L,

    corrob = corrob,

    method_text = paste(
      "Cutoff anchor = nearest d2 zero immediately LEFT of the fixed leading-edge 5000 mark.",
      "Terminal start = nearest d2 zero immediately RIGHT of the fixed leading-edge 5000 mark.",
      "The fixed leading-edge 5000 mark therefore lies inside the final study interval.",
      sep = " "
    )
  )
}

# =============================================================================
# FIGURE BUILDING
# =============================================================================
#
# METHODS SUMMARY
# Every panel marks:
# - cutoff anchor with a filled circle and dashed line
# - leading-edge 5000 reference with a diamond and dot-dash line
# - terminal start with an open circle and dotted line
#
# Method descriptions are placed on the LEFT.
# Numeric summaries are placed on the RIGHT.
# =============================================================================

build_dataset_panel <- function(onset_info, title_prefix, out_file) {
  df <- onset_info$curve_df

  cutoff_x <- onset_info$cutoff_center_rank
  ref5000_x <- onset_info$leading_edge_5000_rank
  terminal_x <- onset_info$terminal_start_rank
  interval_min_x <- onset_info$cutoff_range_rank_min
  interval_max_x <- onset_info$cutoff_range_rank_max

  x_rng <- range(df$rank, na.rm = TRUE)
  x_left <- x_rng[1] + 0.05 * diff(x_rng)
  x_right <- x_rng[1] + 0.73 * diff(x_rng)

  abs_ymax <- max(df$abs_loading, na.rm = TRUE)
  var_ymax <- max(df$var_fit, na.rm = TRUE)
  deriv_ymax <- max(c(df$d1_sm, df$d2_sm), na.rm = TRUE)
  nb_ymax <- max(c(df$nb_support, df$log_nb1, df$log_nb2, df$amu_sm, df$nb_gap_sm), na.rm = TRUE)

  cutoff_pt_abs <- df[df$rank == cutoff_x, , drop = FALSE]
  ref_pt_abs <- df[df$rank == ref5000_x, , drop = FALSE]
  terminal_pt_abs <- df[df$rank == terminal_x, , drop = FALSE]

  cutoff_pt_var <- df[df$rank == cutoff_x, , drop = FALSE]
  ref_pt_var <- df[df$rank == ref5000_x, , drop = FALSE]
  terminal_pt_var <- df[df$rank == terminal_x, , drop = FALSE]

  cutoff_pt_deriv <- data.frame(rank = cutoff_x, d2 = df$d2_sm[df$rank == cutoff_x], stringsAsFactors = FALSE)
  ref_pt_deriv <- data.frame(rank = ref5000_x, d2 = df$d2_sm[df$rank == ref5000_x], stringsAsFactors = FALSE)
  terminal_pt_deriv <- data.frame(rank = terminal_x, d2 = df$d2_sm[df$rank == terminal_x], stringsAsFactors = FALSE)

  cutoff_pt_nb <- data.frame(rank = cutoff_x, value = df$nb_support[df$rank == cutoff_x], stringsAsFactors = FALSE)
  ref_pt_nb <- data.frame(rank = ref5000_x, value = df$nb_support[df$rank == ref5000_x], stringsAsFactors = FALSE)
  terminal_pt_nb <- data.frame(rank = terminal_x, value = df$nb_support[df$rank == terminal_x], stringsAsFactors = FALSE)

  corr <- onset_info$corrob

  left_abs_label <- make_label_df(
    x_left, abs_ymax,
    paste(
      "Absolute loading panel",
      "Leading edge is on the RIGHT",
      "Filled circle = cutoff anchor",
      "Diamond = fixed leading-edge 5000 mark",
      "Open circle = terminal start",
      sep = "\n"
    )
  )

  right_abs_label <- make_label_df(
    x_right, abs_ymax,
    paste0(
      "Cutoff anchor rank = ", cutoff_x,
      "\nReference rank (5000 from right) = ", ref5000_x,
      "\nTerminal start rank = ", terminal_x,
      "\nFinal interval = [", interval_min_x, ", ", interval_max_x, "]",
      "\nPre-EVS remainder = ", onset_info$pre_evs_remainder_size,
      "\nPre-EVS leading edge = ", onset_info$pre_evs_leading_edge_size
    )
  )

  left_var_label <- make_label_df(
    x_left, var_ymax,
    paste(
      "Variance curve method",
      "The smooth is fitted to empirical log(1+variance).",
      "Cutoff anchor = nearest d2 zero immediately LEFT of the fixed leading-edge 5000 mark.",
      "Terminal start = nearest d2 zero immediately RIGHT of the fixed leading-edge 5000 mark.",
      sep = "\n"
    )
  )

  left_deriv_label <- make_label_df(
    x_left, deriv_ymax,
    paste(
      "Derivative method",
      "d2 zero-crossings are selected only from sign changes or exact zeros.",
      "No amplitude threshold is used.",
      "No run-length filter is used.",
      "The fixed leading-edge 5000 mark is inside the final interval.",
      sep = "\n"
    )
  )

  left_nb_label <- make_label_df(
    x_left, nb_ymax,
    paste(
      "NB corroboration",
      "Right-side NB2 region = all genes from cutoff anchor to the end.",
      "Matched left-side NB1 region = same number of genes immediately to the left.",
      "These values corroborate, but do not define, the geometric split.",
      sep = "\n"
    )
  )

  right_nb_label <- make_label_df(
    x_right, nb_ymax,
    paste0(
      "Left median log(NB2) = ", fmt_num(corr$left_median_log_nb2),
      "\nRight median log(NB2) = ", fmt_num(corr$right_median_log_nb2),
      "\nRight-left log(NB2) diff = ", fmt_num(corr$right_left_log_nb2_diff),
      "\nLeft median NB2-NB1 = ", fmt_num(corr$left_median_nb_gap),
      "\nRight median NB2-NB1 = ", fmt_num(corr$right_median_nb_gap),
      "\nRight-left NB2-NB1 diff = ", fmt_num(corr$right_left_nb_gap_diff),
      "\nLeft median log(alpha*mu) = ", fmt_num(corr$left_median_log_alpha_mu),
      "\nRight median log(alpha*mu) = ", fmt_num(corr$right_median_log_alpha_mu),
      "\nRight-left log(alpha*mu) diff = ", fmt_num(corr$right_left_log_alpha_mu_diff)
    )
  )

  right_band_label <- make_label_df(
    x_right, max(df$nb_support, na.rm = TRUE),
    paste0(
      "Final NB-supported summary",
      "\nCenter = ", cutoff_x,
      "\nInterval = [", interval_min_x, ", ", interval_max_x, "]",
      "\nReference rank = ", ref5000_x,
      "\nTerminal start = ", terminal_x
    )
  )

  p1 <- ggplot(df, aes(x = rank, y = abs_loading)) +
    annotate(
      "rect",
      xmin = interval_min_x, xmax = interval_max_x,
      ymin = -Inf, ymax = Inf,
      fill = COLORS["interval_fill"], alpha = 0.18
    ) +
    geom_line(color = COLORS["abs_loading"], linewidth = 0.95, na.rm = TRUE) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.85, color = COLORS["cutoff_anchor"]) +
    geom_vline(xintercept = ref5000_x, linetype = 4, linewidth = 0.85, color = COLORS["ref5000"]) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.85, color = COLORS["terminal_start"]) +
    geom_point(data = cutoff_pt_abs, aes(x = rank, y = abs_loading),
               shape = 16, size = 2.5, inherit.aes = FALSE, color = COLORS["cutoff_anchor"]) +
    geom_point(data = ref_pt_abs, aes(x = rank, y = abs_loading),
               shape = 18, size = 2.8, inherit.aes = FALSE, color = COLORS["ref5000"]) +
    geom_point(data = terminal_pt_abs, aes(x = rank, y = abs_loading),
               shape = 1, size = 2.9, stroke = 1.0, inherit.aes = FALSE, color = COLORS["terminal_start"]) +
    geom_label(data = left_abs_label, aes(x = x, y = y, label = label),
               inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.55,
               label.padding = unit(0.15, "lines"), label.r = unit(0.10, "lines"),
               linewidth = 0.25, fill = alpha("white", 0.96)) +
    geom_label(data = right_abs_label, aes(x = x, y = y, label = label),
               inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.55,
               label.padding = unit(0.15, "lines"), label.r = unit(0.10, "lines"),
               linewidth = 0.25, fill = alpha("white", 0.96)) +
    labs(
      title = paste0(title_prefix, ": absolute PC1 loading series"),
      subtitle = "Leading edge is on the RIGHT",
      x = "EVS rank",
      y = "|PC1 loading|"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.margin = margin(8, 38, 8, 10))

  p2 <- ggplot(df, aes(x = rank, y = var_fit)) +
    annotate(
      "rect",
      xmin = interval_min_x, xmax = interval_max_x,
      ymin = -Inf, ymax = Inf,
      fill = COLORS["interval_fill"], alpha = 0.18
    ) +
    geom_line(color = COLORS["variance_fit"], linewidth = 1.0, na.rm = TRUE) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.85, color = COLORS["cutoff_anchor"]) +
    geom_vline(xintercept = ref5000_x, linetype = 4, linewidth = 0.85, color = COLORS["ref5000"]) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.85, color = COLORS["terminal_start"]) +
    geom_point(data = cutoff_pt_var, aes(x = rank, y = var_fit),
               shape = 16, size = 2.5, inherit.aes = FALSE, color = COLORS["cutoff_anchor"]) +
    geom_point(data = ref_pt_var, aes(x = rank, y = var_fit),
               shape = 18, size = 2.8, inherit.aes = FALSE, color = COLORS["ref5000"]) +
    geom_point(data = terminal_pt_var, aes(x = rank, y = var_fit),
               shape = 1, size = 2.9, stroke = 1.0, inherit.aes = FALSE, color = COLORS["terminal_start"]) +
    geom_label(data = left_var_label, aes(x = x, y = y, label = label),
               inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.35,
               label.padding = unit(0.15, "lines"), label.r = unit(0.10, "lines"),
               linewidth = 0.25, fill = alpha("white", 0.96)) +
    labs(
      title = paste0(title_prefix, ": smoothed empirical variance curve"),
      subtitle = "Selected geometric points are marked on the curve",
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.margin = margin(8, 38, 8, 10))

  deriv_df <- bind_rows(
    data.frame(rank = df$rank, value = df$d2_sm, series = "Smoothed curvature d2", stringsAsFactors = FALSE),
    data.frame(rank = df$rank, value = df$d1_sm, series = "Smoothed slope d1", stringsAsFactors = FALSE)
  )
  deriv_df$series <- factor(
    deriv_df$series,
    levels = c("Smoothed curvature d2", "Smoothed slope d1")
  )

  p3 <- ggplot(deriv_df, aes(x = rank, y = value, color = series)) +
    annotate(
      "rect",
      xmin = interval_min_x, xmax = interval_max_x,
      ymin = -Inf, ymax = Inf,
      fill = COLORS["interval_fill"], alpha = 0.18
    ) +
    geom_hline(yintercept = 0, linewidth = 0.55, color = "black") +
    geom_line(linewidth = 0.95, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "Smoothed curvature d2" = COLORS["d2"],
        "Smoothed slope d1" = COLORS["d1"]
      )
    ) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.85, color = COLORS["cutoff_anchor"]) +
    geom_vline(xintercept = ref5000_x, linetype = 4, linewidth = 0.85, color = COLORS["ref5000"]) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.85, color = COLORS["terminal_start"]) +
    geom_point(data = cutoff_pt_deriv, aes(x = rank, y = d2),
               shape = 16, size = 2.5, inherit.aes = FALSE, color = COLORS["cutoff_anchor"]) +
    geom_point(data = ref_pt_deriv, aes(x = rank, y = d2),
               shape = 18, size = 2.8, inherit.aes = FALSE, color = COLORS["ref5000"]) +
    geom_point(data = terminal_pt_deriv, aes(x = rank, y = d2),
               shape = 1, size = 2.9, stroke = 1.0, inherit.aes = FALSE, color = COLORS["terminal_start"]) +
    geom_label(data = left_deriv_label, aes(x = x, y = y, label = label),
               inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.30,
               label.padding = unit(0.15, "lines"), label.r = unit(0.10, "lines"),
               linewidth = 0.25, fill = alpha("white", 0.96)) +
    labs(
      title = paste0(title_prefix, ": derivative support"),
      subtitle = "Cutoff anchor and terminal start come only from the selected d2 zero-crossings around the fixed leading-edge 5000 mark",
      x = "EVS rank",
      y = "Derivative value",
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(8, 38, 8, 10)
    )

  nb_df <- bind_rows(
    data.frame(rank = df$rank, value = df$nb_support, series = "Combined NB support", stringsAsFactors = FALSE),
    data.frame(rank = df$rank, value = df$log_nb1, series = "NB1 = mu", stringsAsFactors = FALSE),
    data.frame(rank = df$rank, value = df$log_nb2, series = "NB2 = variance - mu", stringsAsFactors = FALSE),
    data.frame(rank = df$rank, value = df$amu_sm, series = "Smoothed log(alpha*mu)", stringsAsFactors = FALSE),
    data.frame(rank = df$rank, value = df$nb_gap_sm, series = "Smoothed log(NB2+1) - log(NB1+1)", stringsAsFactors = FALSE)
  )
  nb_df$series <- factor(
    nb_df$series,
    levels = c(
      "Combined NB support",
      "NB1 = mu",
      "NB2 = variance - mu",
      "Smoothed log(alpha*mu)",
      "Smoothed log(NB2+1) - log(NB1+1)"
    )
  )

  p4 <- ggplot(nb_df, aes(x = rank, y = value, color = series)) +
    annotate(
      "rect",
      xmin = interval_min_x, xmax = interval_max_x,
      ymin = -Inf, ymax = Inf,
      fill = COLORS["interval_fill"], alpha = 0.18
    ) +
    geom_hline(yintercept = onset_info$nb_support_threshold, linetype = 3, linewidth = 0.55, color = "grey40") +
    geom_line(linewidth = 0.90, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "Combined NB support" = COLORS["nb_support"],
        "NB1 = mu" = COLORS["nb1"],
        "NB2 = variance - mu" = COLORS["nb2"],
        "Smoothed log(alpha*mu)" = COLORS["amu"],
        "Smoothed log(NB2+1) - log(NB1+1)" = COLORS["nb_gap"]
      )
    ) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.85, color = COLORS["cutoff_anchor"]) +
    geom_vline(xintercept = ref5000_x, linetype = 4, linewidth = 0.85, color = COLORS["ref5000"]) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.85, color = COLORS["terminal_start"]) +
    geom_point(data = cutoff_pt_nb, aes(x = rank, y = value),
               shape = 16, size = 2.5, inherit.aes = FALSE, color = COLORS["cutoff_anchor"]) +
    geom_point(data = ref_pt_nb, aes(x = rank, y = value),
               shape = 18, size = 2.8, inherit.aes = FALSE, color = COLORS["ref5000"]) +
    geom_point(data = terminal_pt_nb, aes(x = rank, y = value),
               shape = 1, size = 2.9, stroke = 1.0, inherit.aes = FALSE, color = COLORS["terminal_start"]) +
    geom_label(data = left_nb_label, aes(x = x, y = y, label = label),
               inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.15,
               label.padding = unit(0.15, "lines"), label.r = unit(0.10, "lines"),
               linewidth = 0.25, fill = alpha("white", 0.96)) +
    geom_label(data = right_nb_label, aes(x = x, y = y, label = label),
               inherit.aes = FALSE, hjust = 0, vjust = 1, size = 1.95,
               label.padding = unit(0.15, "lines"), label.r = unit(0.10, "lines"),
               linewidth = 0.25, fill = alpha("white", 0.96)) +
    labs(
      title = paste0(title_prefix, ": NB1 / NB2 / alpha*mu support"),
      subtitle = "Right-of-cutoff elevation in NB2 and alpha*mu corroborates the leading-edge interpretation",
      x = "EVS rank",
      y = "Support value",
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(8, 38, 8, 10)
    )

  p5 <- ggplot(df, aes(x = rank, y = nb_support)) +
    annotate(
      "rect",
      xmin = interval_min_x, xmax = interval_max_x,
      ymin = -Inf, ymax = Inf,
      fill = COLORS["interval_fill"], alpha = 0.18
    ) +
    geom_hline(yintercept = onset_info$nb_support_threshold, linetype = 3, linewidth = 0.55, color = "grey40") +
    geom_line(color = COLORS["nb_support"], linewidth = 1.0, na.rm = TRUE) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.85, color = COLORS["cutoff_anchor"]) +
    geom_vline(xintercept = ref5000_x, linetype = 4, linewidth = 0.85, color = COLORS["ref5000"]) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.85, color = COLORS["terminal_start"]) +
    geom_point(data = cutoff_pt_nb, aes(x = rank, y = value),
               shape = 16, size = 2.5, inherit.aes = FALSE, color = COLORS["cutoff_anchor"]) +
    geom_point(data = ref_pt_nb, aes(x = rank, y = value),
               shape = 18, size = 2.8, inherit.aes = FALSE, color = COLORS["ref5000"]) +
    geom_point(data = terminal_pt_nb, aes(x = rank, y = value),
               shape = 1, size = 2.9, stroke = 1.0, inherit.aes = FALSE, color = COLORS["terminal_start"]) +
    geom_label(data = right_band_label, aes(x = x, y = y, label = label),
               inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.15,
               label.padding = unit(0.15, "lines"), label.r = unit(0.10, "lines"),
               linewidth = 0.25, fill = alpha("white", 0.96)) +
    labs(
      title = paste0(title_prefix, ": NB-supported summary"),
      subtitle = "The geometric interval remains primary; NB support is corroborative",
      x = "EVS rank",
      y = "Combined NB support"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.margin = margin(8, 38, 8, 10))

  png(out_file, width = 3000, height = 3900, res = 240)
  grid.arrange(p1, p2, p3, p4, p5, ncol = 1)
  dev.off()
}

build_range_panel <- function(ctrl_onset, trt_onset, title_prefix, out_file) {
  ctrl_df <- ctrl_onset$curve_df
  trt_df <- trt_onset$curve_df

  x_all <- c(ctrl_df$rank, trt_df$rank)
  y_all <- c(ctrl_df$var_fit, trt_df$var_fit)

  x_left <- min(x_all, na.rm = TRUE) + 0.05 * diff(range(x_all, na.rm = TRUE))
  x_right <- min(x_all, na.rm = TRUE) + 0.68 * diff(range(x_all, na.rm = TRUE))
  y_top <- max(y_all, na.rm = TRUE)

  label_left <- make_label_df(
    x_left, y_top,
    paste(
      "Treatment/control comparison",
      "Both arms use the same manuscript rule.",
      "Cutoff anchor = nearest d2 zero left of the fixed leading-edge 5000 mark.",
      "Terminal start = nearest d2 zero right of the fixed leading-edge 5000 mark.",
      sep = "\n"
    )
  )

  label_right <- make_label_df(
    x_right, y_top,
    paste0(
      "Control interval = [", ctrl_onset$cutoff_range_rank_min, ", ", ctrl_onset$cutoff_range_rank_max, "]",
      "\nTreatment interval = [", trt_onset$cutoff_range_rank_min, ", ", trt_onset$cutoff_range_rank_max, "]",
      "\nControl reference = ", ctrl_onset$leading_edge_5000_rank,
      "\nTreatment reference = ", trt_onset$leading_edge_5000_rank
    )
  )

  p <- ggplot() +
    annotate(
      "rect",
      xmin = min(ctrl_onset$cutoff_range_rank_min, trt_onset$cutoff_range_rank_min),
      xmax = max(ctrl_onset$cutoff_range_rank_max, trt_onset$cutoff_range_rank_max),
      ymin = -Inf, ymax = Inf,
      fill = COLORS["interval_fill"], alpha = 0.12
    ) +
    geom_line(data = ctrl_df, aes(x = rank, y = var_fit, color = "Control variance fit"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(data = trt_df, aes(x = rank, y = var_fit, color = "Treatment variance fit"), linewidth = 1.0, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "Control variance fit" = "#1B9E77",
        "Treatment variance fit" = "#7570B3"
      )
    ) +
    geom_vline(xintercept = ctrl_onset$cutoff_center_rank, linetype = 2, linewidth = 0.85, color = "#252525") +
    geom_vline(xintercept = trt_onset$cutoff_center_rank, linetype = 3, linewidth = 0.85, color = "#636363") +
    geom_vline(xintercept = ctrl_onset$leading_edge_5000_rank, linetype = 4, linewidth = 0.70, color = "#8C510A") +
    geom_vline(xintercept = trt_onset$leading_edge_5000_rank, linetype = 4, linewidth = 0.70, color = "#BF812D") +
    geom_label(data = label_left, aes(x = x, y = y, label = label),
               inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.7,
               label.padding = unit(0.15, "lines"), label.r = unit(0.10, "lines"),
               linewidth = 0.25, fill = alpha("white", 0.96)) +
    geom_label(data = label_right, aes(x = x, y = y, label = label),
               inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.7,
               label.padding = unit(0.15, "lines"), label.r = unit(0.10, "lines"),
               linewidth = 0.25, fill = alpha("white", 0.96)) +
    labs(
      title = paste0(title_prefix, ": treatment/control variance intervals"),
      subtitle = "Both arms are shown on the same variance-fit axis",
      x = "EVS rank",
      y = "Fitted log(1 + variance)",
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(8, 38, 8, 10)
    )

  png(out_file, width = 3000, height = 1400, res = 240)
  print(p)
  dev.off()
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

message("SEQUENCE.R started"); flush.console()
message("Resolving count file..."); flush.console()
count_file <- resolve_counts_file(count_file_hint)
message("Using count file: ", count_file); flush.console()

message("Reading count matrix..."); flush.console()
loaded <- read_count_matrix(count_file, meta_all$id)
count_mat <- loaded$count_matrix
annot_df <- loaded$annot_df
message("Count matrix dimensions: ", nrow(count_mat), " features x ", ncol(count_mat), " samples"); flush.console()

overall_rows <- list()

for (i in seq_len(nrow(comparison_table))) {
  comparison_row <- comparison_table[i, , drop = FALSE]
  cmp_name <- comparison_row$comparison_name[[1]]

  message("Processing comparison: ", cmp_name); flush.console()

  cmp_dir <- file.path(out_root, paste0(cmp_name, "_cutoff_folder"))
  dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

  comp <- subset_comparison(count_mat, comparison_row, meta_all)
  cmp_counts <- comp$count_matrix
  feature_ids <- rownames(cmp_counts)

  evs_tbl <- build_evs_table(cmp_counts, comp$ctrl_ids, comp$trt_ids, feature_ids)
  metrics_tbl <- compute_feature_metrics_empirical(cmp_counts, comp$ctrl_ids, comp$trt_ids, feature_ids)

  full_tbl <- evs_tbl %>%
    left_join(metrics_tbl, by = "feature_id") %>%
    left_join(annot_df, by = "feature_id")

  utils::write.csv(
    full_tbl,
    file.path(cmp_dir, paste0(cmp_name, "_feature_level_metrics.csv")),
    row.names = FALSE
  )

  ctrl_rank_df <- build_rank_series(full_tbl, "control")
  trt_rank_df  <- build_rank_series(full_tbl, "treatment")

  utils::write.csv(
    ctrl_rank_df,
    file.path(cmp_dir, paste0(cmp_name, "_control_rank_series.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    trt_rank_df,
    file.path(cmp_dir, paste0(cmp_name, "_treatment_rank_series.csv")),
    row.names = FALSE
  )

  message("Selecting control cutoff for ", cmp_name); flush.console()
  ctrl_onset <- select_cutoff_interval(
    ctrl_rank_df,
    spline_spar = spline_spar,
    zero_tol = zero_tol,
    nb_smooth_window = nb_smooth_window,
    leading_edge_fixed_k = leading_edge_fixed_k
  )
  message(
    "Control cutoff anchor rank: ", ctrl_onset$cutoff_center_rank,
    " | reference rank: ", ctrl_onset$leading_edge_5000_rank,
    " | terminal anchor: ", ctrl_onset$terminal_start_rank,
    " | interval: [", ctrl_onset$cutoff_range_rank_min, ", ", ctrl_onset$cutoff_range_rank_max, "]",
    " | pre-EVS remainder: ", ctrl_onset$pre_evs_remainder_size,
    " | pre-EVS leading edge: ", ctrl_onset$pre_evs_leading_edge_size
  ); flush.console()

  message("Selecting treatment cutoff for ", cmp_name); flush.console()
  trt_onset <- select_cutoff_interval(
    trt_rank_df,
    spline_spar = spline_spar,
    zero_tol = zero_tol,
    nb_smooth_window = nb_smooth_window,
    leading_edge_fixed_k = leading_edge_fixed_k
  )
  message(
    "Treatment cutoff anchor rank: ", trt_onset$cutoff_center_rank,
    " | reference rank: ", trt_onset$leading_edge_5000_rank,
    " | terminal anchor: ", trt_onset$terminal_start_rank,
    " | interval: [", trt_onset$cutoff_range_rank_min, ", ", trt_onset$cutoff_range_rank_max, "]",
    " | pre-EVS remainder: ", trt_onset$pre_evs_remainder_size,
    " | pre-EVS leading edge: ", trt_onset$pre_evs_leading_edge_size
  ); flush.console()

  build_dataset_panel(
    ctrl_onset,
    paste0(cmp_name, " control"),
    file.path(cmp_dir, paste0(cmp_name, "_control_rank_panel.png"))
  )

  build_dataset_panel(
    trt_onset,
    paste0(cmp_name, " treatment"),
    file.path(cmp_dir, paste0(cmp_name, "_treatment_rank_panel.png"))
  )

  build_range_panel(
    ctrl_onset,
    trt_onset,
    cmp_name,
    file.path(cmp_dir, paste0(cmp_name, "_cutoff_range_panel.png"))
  )

  utils::write.csv(
    ctrl_onset$zero_crossings,
    file.path(cmp_dir, paste0(cmp_name, "_control_zero_crossings_all.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    trt_onset$zero_crossings,
    file.path(cmp_dir, paste0(cmp_name, "_treatment_zero_crossings_all.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    ctrl_onset$selected_points,
    file.path(cmp_dir, paste0(cmp_name, "_control_selected_points.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    trt_onset$selected_points,
    file.path(cmp_dir, paste0(cmp_name, "_treatment_selected_points.csv")),
    row.names = FALSE
  )

  ctrl_corr <- ctrl_onset$corrob
  trt_corr <- trt_onset$corrob

  cutoff_summary <- tibble(
    comparison = cmp_name,
    total_features = ctrl_onset$total_features,
    leading_edge_fixed_k = leading_edge_fixed_k,

    control_method_text = ctrl_onset$method_text,
    control_cutoff_anchor_rank = ctrl_onset$cutoff_center_rank,
    control_reference_rank_5000_from_right = ctrl_onset$leading_edge_5000_rank,
    control_terminal_start_rank = ctrl_onset$terminal_start_rank,
    control_interval_rank_min = ctrl_onset$cutoff_range_rank_min,
    control_interval_rank_max = ctrl_onset$cutoff_range_rank_max,
    control_center_nb_support = ctrl_onset$center_nb_support,
    control_nb_support_threshold = ctrl_onset$nb_support_threshold,
    control_pre_evs_remainder_size = ctrl_onset$pre_evs_remainder_size,
    control_pre_evs_leading_edge_size = ctrl_onset$pre_evs_leading_edge_size,
    control_left_n = ctrl_corr$left_n,
    control_right_n = ctrl_corr$right_n,
    control_left_median_log_nb2 = ctrl_corr$left_median_log_nb2,
    control_right_median_log_nb2 = ctrl_corr$right_median_log_nb2,
    control_right_left_log_nb2_diff = ctrl_corr$right_left_log_nb2_diff,
    control_left_median_nb_gap = ctrl_corr$left_median_nb_gap,
    control_right_median_nb_gap = ctrl_corr$right_median_nb_gap,
    control_right_left_nb_gap_diff = ctrl_corr$right_left_nb_gap_diff,
    control_left_median_log_alpha_mu = ctrl_corr$left_median_log_alpha_mu,
    control_right_median_log_alpha_mu = ctrl_corr$right_median_log_alpha_mu,
    control_right_left_log_alpha_mu_diff = ctrl_corr$right_left_log_alpha_mu_diff,

    treatment_method_text = trt_onset$method_text,
    treatment_cutoff_anchor_rank = trt_onset$cutoff_center_rank,
    treatment_reference_rank_5000_from_right = trt_onset$leading_edge_5000_rank,
    treatment_terminal_start_rank = trt_onset$terminal_start_rank,
    treatment_interval_rank_min = trt_onset$cutoff_range_rank_min,
    treatment_interval_rank_max = trt_onset$cutoff_range_rank_max,
    treatment_center_nb_support = trt_onset$center_nb_support,
    treatment_nb_support_threshold = trt_onset$nb_support_threshold,
    treatment_pre_evs_remainder_size = trt_onset$pre_evs_remainder_size,
    treatment_pre_evs_leading_edge_size = trt_onset$pre_evs_leading_edge_size,
    treatment_left_n = trt_corr$left_n,
    treatment_right_n = trt_corr$right_n,
    treatment_left_median_log_nb2 = trt_corr$left_median_log_nb2,
    treatment_right_median_log_nb2 = trt_corr$right_median_log_nb2,
    treatment_right_left_log_nb2_diff = trt_corr$right_left_log_nb2_diff,
    treatment_left_median_nb_gap = trt_corr$left_median_nb_gap,
    treatment_right_median_nb_gap = trt_corr$right_median_nb_gap,
    treatment_right_left_nb_gap_diff = trt_corr$right_left_nb_gap_diff,
    treatment_left_median_log_alpha_mu = trt_corr$left_median_log_alpha_mu,
    treatment_right_median_log_alpha_mu = trt_corr$right_median_log_alpha_mu,
    treatment_right_left_log_alpha_mu_diff = trt_corr$right_left_log_alpha_mu_diff
  )

  utils::write.csv(
    cutoff_summary,
    file.path(cmp_dir, paste0(cmp_name, "_cutoff_summary.csv")),
    row.names = FALSE
  )

  overall_rows[[cmp_name]] <- cutoff_summary
}

overall_summary <- bind_rows(overall_rows)

utils::write.csv(
  overall_summary,
  file.path(out_root, "overall_cutoff_summary.csv"),
  row.names = FALSE
)

message("Done. Outputs written to: ", out_root); flush.console()
