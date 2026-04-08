# =============================================================================
# SEQUENCE STAGE 1 FINAL PRE-EVS CUTOFF SELECTOR
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# Features are ranked by absolute PC1 loading, with lower loading on the left
# and higher loading on the right. This ordered axis is treated as a geometric
# transition axis before any downstream eigenvector splitting. A smoothed
# empirical variance curve is fitted across the full ranked series, and its
# first and second derivatives are computed from the same smooth.
#
# Final geometric rule used here
# 1. Rank 5000 is retained as a fixed manuscript reference point.
# 2. The study interval is defined by the first two terminal zero-crossings of
#    the smoothed second derivative on the right side of the ranked series.
# 3. The earlier of those two terminal zero-crossings is the cutoff anchor.
# 4. The later of those two terminal zero-crossings is the terminal start.
#
# NB corroboration rule used here
# 1. The NB2-like right-side reference region is defined as all genes from the
#    cutoff anchor to the end of the ranked series.
# 2. Let that right-side region contain m genes.
# 3. The matched NB1-side reference region is defined as the same number m of
#    genes immediately to the left of the cutoff anchor, truncated at rank 1
#    if necessary.
# 4. Right-versus-left summaries of NB2, NB2-NB1 contrast, and log(alpha*mu)
#    are reported explicitly so the figure itself shows why the right side is
#    being interpreted as more NB2-like.
#
# Figure design rule used here
# 1. Method text is placed on the left side of panels.
# 2. Numeric cutoff summaries are placed on the right side of panels.
# 3. The cutoff anchor is shown with a filled circle and dashed vertical line.
# 4. The terminal start is shown with an open circle and dotted vertical line.
# 5. Rank 5000 is shown with a diamond and dot-dash vertical line.
# 6. Colors in the legends are identical to the plotted series.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(scales)
  library(tibble)
})

options(warn = 1)

# =============================================================================
# USER SETTINGS
# =============================================================================

repo_dir <- getwd()
count_file <- file.path(repo_dir, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
out_root <- file.path(repo_dir, "exports", "variance_derivative_nb_range_final_manuscript")
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
terminal_fraction_start <- 0.70
zero_tol <- 1e-8
fixed_rank_reference <- 5000L

# =============================================================================
# COLOR MAP
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# These colors are fixed so that each quantity is represented identically in
# every panel and in every legend.
# =============================================================================

COLORS <- c(
  variance_fit   = "#1b1b1b",
  d1             = "#00a6d6",
  d2             = "#d1495b",
  nb_support     = "#111111",
  nb1            = "#d9b44a",
  nb2            = "#7fcdbb",
  amu            = "#1f78ff",
  nb_gap         = "#d65cff",
  interval_fill  = "#bdbdbd",
  cutoff_anchor  = "#111111",
  terminal_anchor= "#111111",
  ref5000        = "#8c510a"
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
  if (length(hit)) return(hit[1])
  names(df)[1]
}

detect_gene_symbol_column <- function(df) {
  candidates <- c("Symbol", "symbol", "GeneSymbol", "gene_symbol", "gene", "Gene")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) return(hit[1])
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
  stats::median(x)
}

safe_ratio <- function(num, den) {
  if (!is.finite(num) || !is.finite(den) || den == 0) return(NA_real_)
  num / den
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

scale01 <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(out)

  rng <- range(x[ok], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || abs(rng[2] - rng[1]) < .Machine$double.eps) {
    out[ok] <- 0
    return(out)
  }

  out[ok] <- (x[ok] - rng[1]) / (rng[2] - rng[1])
  out
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

# =============================================================================
# DATA INGESTION
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# The raw count file is read once. Feature identifiers and, if present, gene
# symbols are retained. Only samples explicitly represented in the comparison
# metadata are used. Duplicate feature identifiers are collapsed by summation.
# =============================================================================

read_count_matrix <- function(path, meta_ids) {
  raw_df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  feature_col <- detect_feature_id_column(raw_df)
  symbol_col <- detect_gene_symbol_column(raw_df)
  sample_cols <- intersect(meta_ids, names(raw_df))

  if (length(sample_cols) == 0L) {
    stop("No count columns matched metadata sample IDs.", call. = FALSE)
  }

  annot_df <- data.frame(
    feature_id = as.character(raw_df[[feature_col]]),
    stringsAsFactors = FALSE
  )
  annot_df$gene_symbol <- if (!is.null(symbol_col)) {
    as.character(raw_df[[symbol_col]])
  } else {
    annot_df$feature_id
  }

  keep <- !is.na(annot_df$feature_id) & nzchar(annot_df$feature_id)
  annot_df <- annot_df[keep, , drop = FALSE]

  count_df <- raw_df[keep, sample_cols, drop = FALSE]
  count_mat <- as.matrix(count_df)
  mode(count_mat) <- "numeric"
  rownames(count_mat) <- annot_df$feature_id
  colnames(count_mat) <- sample_cols

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
# EVS RANKING AND EMPIRICAL NB QUANTITIES
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# EVS ranking is derived from the absolute value of the first PC loading within
# each arm after log2(count+1) transformation. Empirical mean, variance,
# empirical alpha, NB1, NB2, and alpha*mu quantities are computed directly from
# the observed counts for each arm.
# =============================================================================

compute_group_pc1_loadings <- function(count_mat, group_cols) {
  mat <- count_mat[, group_cols, drop = FALSE]
  mat <- log2(mat + 1)
  mat <- t(mat)

  if (nrow(mat) < 2L) return(rep(NA_real_, ncol(mat)))

  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  abs(pca$rotation[, 1L])
}

build_evs_table <- function(count_mat, ctrl_cols, trt_cols, feature_ids) {
  ctrl_load <- compute_group_pc1_loadings(count_mat, ctrl_cols)
  trt_load  <- compute_group_pc1_loadings(count_mat, trt_cols)

  tibble(
    feature_id = feature_ids,
    ctrl_abs_loading = ctrl_load,
    trt_abs_loading  = trt_load
  ) %>%
    mutate(
      rank_ctrl = rank(ctrl_abs_loading, ties.method = "first"),
      rank_trt  = rank(trt_abs_loading, ties.method = "first")
    )
}

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
      alpha_emp_ctrl = ifelse(mu_ctrl > 0, (var_ctrl - mu_ctrl) / (mu_ctrl ^ 2), NA_real_),
      alpha_emp_trt  = ifelse(mu_trt  > 0, (var_trt  - mu_trt)  / (mu_trt  ^ 2), NA_real_),
      nb1_ctrl = mu_ctrl,
      nb1_trt  = mu_trt,
      nb2_ctrl = pmax(var_ctrl - mu_ctrl, 0),
      nb2_trt  = pmax(var_trt  - mu_trt, 0),
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
        alpha = alpha_emp_ctrl,
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
        alpha = alpha_emp_trt,
        nb1 = nb1_trt,
        nb2 = nb2_trt,
        alpha_mu = alpha_mu_trt
      ) %>%
      arrange(rank)
  }

  out$log_variance <- ifelse(is.finite(out$variance) & out$variance >= 0, log1p(out$variance), NA_real_)
  out$log_nb1 <- ifelse(is.finite(out$nb1) & out$nb1 >= 0, log1p(out$nb1), NA_real_)
  out$log_nb2 <- ifelse(is.finite(out$nb2) & out$nb2 >= 0, log1p(out$nb2), NA_real_)
  out$nb_gap <- out$log_nb2 - out$log_nb1
  out$log_alpha_mu <- ifelse(is.finite(out$alpha_mu) & out$alpha_mu > 0, log(out$alpha_mu), NA_real_)

  out
}

# =============================================================================
# DERIVATIVE GEOMETRY
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# The variance curve is smoothed across EVS rank using a smoothing spline. The
# first and second derivatives of that smooth are then evaluated on the same
# EVS rank grid. Geometric anchors are defined only from second-derivative
# zero-crossings. No directional sign constraint, no run-length threshold, and
# no amplitude threshold are imposed. The manuscript anchor pair is the first
# two terminal second-derivative zero-crossings on the right-hand side.
# =============================================================================

find_d2_zero_crossings <- function(d2_vec, zero_tol = 1e-8) {
  x <- as.numeric(d2_vec)
  n <- length(x)
  if (n < 2L) return(integer(0))

  x[!is.finite(x)] <- NA_real_
  x[is.finite(x) & abs(x) <= zero_tol] <- 0
  s <- sign(x)
  s[!is.finite(s)] <- 0

  out <- integer(0)

  for (i in 2:n) {
    x0 <- x[i - 1L]
    x1 <- x[i]
    s0 <- s[i - 1L]
    s1 <- s[i]

    if (!is.finite(x0) || !is.finite(x1)) next

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
  if (!length(zero_idx)) return(data.frame())

  bind_rows(lapply(zero_idx, function(i) {
    data.frame(
      index = i,
      rank = df$rank[i],
      d1_here = df$d1_sm[i],
      d2_here = df$d2_sm[i],
      var_fit_here = df$var_fit[i],
      stringsAsFactors = FALSE
    )
  }))
}

choose_terminal_two_zero_crossings <- function(df,
                                               terminal_fraction_start = 0.70,
                                               zero_tol = 1e-8) {
  zero_idx <- find_d2_zero_crossings(df$d2_sm, zero_tol = zero_tol)
  zero_tbl <- summarize_zero_crossings(df, zero_idx)

  if (!nrow(zero_tbl)) {
    stop("No second-derivative zero-crossings were detected.", call. = FALSE)
  }

  terminal_start_index_min <- max(1L, floor(terminal_fraction_start * nrow(df)))
  terminal_tbl <- zero_tbl[zero_tbl$index >= terminal_start_index_min, , drop = FALSE]

  if (nrow(terminal_tbl) < 2L) {
    terminal_tbl <- zero_tbl
  }
  if (nrow(terminal_tbl) < 2L) {
    stop("Fewer than two second-derivative zero-crossings were detected.", call. = FALSE)
  }

  selected_tbl <- terminal_tbl[(nrow(terminal_tbl) - 1L):nrow(terminal_tbl), , drop = FALSE]
  rownames(selected_tbl) <- NULL

  cutoff_anchor <- selected_tbl[1, , drop = FALSE]
  terminal_anchor <- selected_tbl[2, , drop = FALSE]

  selected_roles <- data.frame(
    role = c("cutoff_anchor", "terminal_anchor"),
    index = c(cutoff_anchor$index, terminal_anchor$index),
    rank = c(cutoff_anchor$rank, terminal_anchor$rank),
    d1_here = c(cutoff_anchor$d1_here, terminal_anchor$d1_here),
    d2_here = c(cutoff_anchor$d2_here, terminal_anchor$d2_here),
    var_fit_here = c(cutoff_anchor$var_fit_here, terminal_anchor$var_fit_here),
    stringsAsFactors = FALSE
  )

  list(
    zero_crossings = zero_tbl,
    selected_zero_crossings = selected_roles,
    cutoff_anchor_index = cutoff_anchor$index,
    cutoff_anchor_rank = cutoff_anchor$rank,
    terminal_anchor_index = terminal_anchor$index,
    terminal_anchor_rank = terminal_anchor$rank,
    mode = paste(
      "Study interval is defined by the first two terminal second-derivative",
      "zero-crossings on the right; the earlier zero is the cutoff anchor and",
      "the later zero is the terminal start."
    )
  )
}

# =============================================================================
# NB SUPPORT AND RIGHT VERSUS MATCHED LEFT COMPARISON
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# Smoothed NB corroboration is built from three quantities on the ranked axis:
# smoothed NB2-NB1 contrast, smoothed log(alpha*mu), and smoothed log(NB2+1).
# These are scaled to a common 0 to 1 range and averaged. This support score
# does not define the cutoff anchor. Instead, it corroborates the geometric
# split by showing whether the right side of the split is more NB2-like than a
# matched left-side region of equal size.
# =============================================================================

compute_nb_support_score <- function(df, nb_smooth_window = 151L) {
  nb_gap_sm <- roll_median(df$nb_gap, nb_smooth_window)
  amu_sm <- roll_median(df$log_alpha_mu, nb_smooth_window)
  log_nb2_sm <- roll_median(df$log_nb2, nb_smooth_window)

  nb_gap_s <- scale01(nb_gap_sm)
  amu_s <- scale01(amu_sm)
  nb2_s <- scale01(log_nb2_sm)

  nb_gap_s[!is.finite(nb_gap_s)] <- 0
  amu_s[!is.finite(amu_s)] <- 0
  nb2_s[!is.finite(nb2_s)] <- 0

  nb_support <- (nb_gap_s + amu_s + nb2_s) / 3

  list(
    nb_gap_sm = nb_gap_sm,
    amu_sm = amu_sm,
    log_nb2_sm = log_nb2_sm,
    nb_support = nb_support
  )
}

compute_left_right_corrob <- function(df, cutoff_anchor_index) {
  n <- nrow(df)

  right_idx <- seq.int(cutoff_anchor_index, n)
  m <- length(right_idx)

  left_end <- cutoff_anchor_index - 1L
  left_start <- max(1L, left_end - m + 1L)
  left_idx <- if (left_end >= left_start) seq.int(left_start, left_end) else integer(0)

  med_left_nb2 <- if (length(left_idx)) safe_median(df$log_nb2[left_idx]) else NA_real_
  med_right_nb2 <- if (length(right_idx)) safe_median(df$log_nb2[right_idx]) else NA_real_

  med_left_gap <- if (length(left_idx)) safe_median(df$nb_gap_sm[left_idx]) else NA_real_
  med_right_gap <- if (length(right_idx)) safe_median(df$nb_gap_sm[right_idx]) else NA_real_

  med_left_amu <- if (length(left_idx)) safe_median(df$amu_sm[left_idx]) else NA_real_
  med_right_amu <- if (length(right_idx)) safe_median(df$amu_sm[right_idx]) else NA_real_

  med_left_support <- if (length(left_idx)) safe_median(df$nb_support[left_idx]) else NA_real_
  med_right_support <- if (length(right_idx)) safe_median(df$nb_support[right_idx]) else NA_real_

  data.frame(
    left_n = length(left_idx),
    right_n = length(right_idx),

    left_median_log_nb2 = med_left_nb2,
    right_median_log_nb2 = med_right_nb2,
    right_left_log_nb2_diff = med_right_nb2 - med_left_nb2,
    right_left_log_nb2_ratio = safe_ratio(med_right_nb2, med_left_nb2),

    left_median_nb_gap = med_left_gap,
    right_median_nb_gap = med_right_gap,
    right_left_nb_gap_diff = med_right_gap - med_left_gap,
    right_left_nb_gap_ratio = safe_ratio(med_right_gap, med_left_gap),

    left_median_log_alpha_mu = med_left_amu,
    right_median_log_alpha_mu = med_right_amu,
    right_left_log_alpha_mu_diff = med_right_amu - med_left_amu,
    right_left_log_alpha_mu_ratio = safe_ratio(med_right_amu, med_left_amu),

    left_median_nb_support = med_left_support,
    right_median_nb_support = med_right_support,
    right_left_nb_support_diff = med_right_support - med_left_support,
    right_left_nb_support_ratio = safe_ratio(med_right_support, med_left_support),

    right_side_more_nb2 = ifelse(is.finite(med_right_nb2 - med_left_nb2), (med_right_nb2 - med_left_nb2) > 0, NA),
    right_side_more_nb_gap = ifelse(is.finite(med_right_gap - med_left_gap), (med_right_gap - med_left_gap) > 0, NA),
    right_side_more_alpha_mu = ifelse(is.finite(med_right_amu - med_left_amu), (med_right_amu - med_left_amu) > 0, NA),

    stringsAsFactors = FALSE
  )
}

# =============================================================================
# MASTER SELECTOR
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# This function combines the smoothed variance fit, derivative geometry, and
# NB corroboration into one object for figure generation and summary export.
# The geometric cutoff anchor and terminal start are defined only from the
# selected terminal second-derivative zero-crossings.
# =============================================================================

select_cutoff_derivative_nb_range <- function(rank_df,
                                              spline_spar = 0.60,
                                              terminal_fraction_start = 0.70,
                                              zero_tol = 1e-8,
                                              nb_smooth_window = 151L,
                                              fixed_rank_reference = 5000L) {
  df <- rank_df
  x <- df$rank
  y <- df$log_variance

  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10L) {
    stop("Not enough finite variance points for smoothing.", call. = FALSE)
  }

  sp_fit <- smooth.spline(x = x[ok], y = y[ok], spar = spline_spar)
  pred0 <- predict(sp_fit, x = x[ok], deriv = 0)
  pred1 <- predict(sp_fit, x = x[ok], deriv = 1)
  pred2 <- predict(sp_fit, x = x[ok], deriv = 2)

  df$var_fit <- NA_real_
  df$d1_sm <- NA_real_
  df$d2_sm <- NA_real_
  df$var_fit[ok] <- pred0$y
  df$d1_sm[ok] <- pred1$y
  df$d2_sm[ok] <- pred2$y

  nb_obj <- compute_nb_support_score(df, nb_smooth_window = nb_smooth_window)
  df$nb_gap_sm <- nb_obj$nb_gap_sm
  df$amu_sm <- nb_obj$amu_sm
  df$log_nb2_sm <- nb_obj$log_nb2_sm
  df$nb_support <- nb_obj$nb_support

  geom_obj <- choose_terminal_two_zero_crossings(
    df = df,
    terminal_fraction_start = terminal_fraction_start,
    zero_tol = zero_tol
  )

  corrob <- compute_left_right_corrob(df, geom_obj$cutoff_anchor_index)

  fixed_ref_rank <- min(fixed_rank_reference, max(df$rank, na.rm = TRUE))
  ref_index <- which.min(abs(df$rank - fixed_ref_rank))

  total_features <- nrow(df)
  pre_evs_remainder_size <- geom_obj$cutoff_anchor_index - 1L
  pre_evs_leading_edge_size <- total_features - geom_obj$cutoff_anchor_index + 1L

  list(
    curve_df = df,
    zero_crossings = geom_obj$zero_crossings,
    selected_zero_crossings = geom_obj$selected_zero_crossings,
    mode = geom_obj$mode,

    cutoff_center_index = geom_obj$cutoff_anchor_index,
    cutoff_center_rank = geom_obj$cutoff_anchor_rank,

    terminal_start_index = geom_obj$terminal_anchor_index,
    terminal_start_rank = geom_obj$terminal_anchor_rank,

    cutoff_range_index_min = geom_obj$cutoff_anchor_index,
    cutoff_range_rank_min = geom_obj$cutoff_anchor_rank,
    cutoff_range_index_max = geom_obj$terminal_anchor_index,
    cutoff_range_rank_max = geom_obj$terminal_anchor_rank,

    fixed_rank_reference = fixed_ref_rank,
    fixed_rank_reference_index = ref_index,

    center_nb_support = df$nb_support[geom_obj$cutoff_anchor_index],
    nb_support_threshold = safe_median(df$nb_support[df$rank >= geom_obj$cutoff_anchor_rank]),

    total_features = total_features,
    pre_evs_remainder_size = pre_evs_remainder_size,
    pre_evs_leading_edge_size = pre_evs_leading_edge_size,

    corrob = corrob
  )
}

# =============================================================================
# FIGURE HELPERS
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# Left-side annotation boxes describe the method used in each panel. Right-side
# annotation boxes report only the final numerical values relevant to the
# selected anchors and interval.
# =============================================================================

make_left_method_box <- function(label_x, label_y, label_text, size = 2.7) {
  annotate(
    "label",
    x = label_x,
    y = label_y,
    hjust = 0,
    vjust = 1,
    size = size,
    label.size = 0.25,
    fill = alpha("white", 0.96),
    label = label_text
  )
}

make_right_numeric_box <- function(label_x, label_y, label_text, size = 2.7) {
  annotate(
    "label",
    x = label_x,
    y = label_y,
    hjust = 0,
    vjust = 1,
    size = size,
    label.size = 0.25,
    fill = alpha("white", 0.96),
    label = label_text
  )
}

# =============================================================================
# DATASET PANEL GENERATION
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# Each dataset panel contains five stacked subplots:
# 1. absolute loading
# 2. smoothed empirical variance
# 3. derivative support
# 4. NB corroboration quantities
# 5. combined NB support
#
# All anchors are displayed in every panel:
# filled circle and dashed line for cutoff anchor
# open circle and dotted line for terminal start
# diamond and dot-dash line for rank 5000
# =============================================================================

build_dataset_panel <- function(rank_df, onset_info, title_prefix, out_file) {
  df <- onset_info$curve_df
  cutoff_x <- onset_info$cutoff_center_rank
  terminal_x <- onset_info$terminal_start_rank
  interval_min_x <- onset_info$cutoff_range_rank_min
  interval_max_x <- onset_info$cutoff_range_rank_max
  ref_x <- onset_info$fixed_rank_reference
  corr <- onset_info$corrob

  x_rng <- range(df$rank, na.rm = TRUE)
  x_left <- x_rng[1] + 0.04 * diff(x_rng)
  x_right <- x_rng[1] + 0.72 * diff(x_rng)

  cutoff_pt_abs <- df[df$rank == cutoff_x, , drop = FALSE]
  terminal_pt_abs <- df[df$rank == terminal_x, , drop = FALSE]
  ref_pt_abs <- df[df$rank == ref_x, , drop = FALSE]

  cutoff_pt_var <- df[df$rank == cutoff_x & is.finite(df$var_fit), , drop = FALSE]
  terminal_pt_var <- df[df$rank == terminal_x & is.finite(df$var_fit), , drop = FALSE]
  ref_pt_var <- df[df$rank == ref_x & is.finite(df$var_fit), , drop = FALSE]

  deriv_pt_cut <- data.frame(rank = cutoff_x, value = df$d2_sm[df$rank == cutoff_x])
  deriv_pt_term <- data.frame(rank = terminal_x, value = df$d2_sm[df$rank == terminal_x])
  deriv_pt_ref <- data.frame(rank = ref_x, value = df$d2_sm[df$rank == ref_x])

  nb_pt_cut <- data.frame(rank = cutoff_x, value = df$nb_support[df$rank == cutoff_x])
  nb_pt_term <- data.frame(rank = terminal_x, value = df$nb_support[df$rank == terminal_x])
  nb_pt_ref <- data.frame(rank = ref_x, value = df$nb_support[df$rank == ref_x])

  left_text_abs <- paste(
    "Absolute loading panel",
    "Filled circle = cutoff anchor",
    "Open circle = terminal start",
    "Diamond = fixed rank 5000",
    sep = "\n"
  )

  right_text_abs <- paste0(
    "Cutoff anchor rank = ", cutoff_x,
    "\nTerminal start rank = ", terminal_x,
    "\nStudy interval = [", interval_min_x, ", ", interval_max_x, "]",
    "\nPre-EVS remainder = ", onset_info$pre_evs_remainder_size,
    "\nPre-EVS leading edge = ", onset_info$pre_evs_leading_edge_size
  )

  left_text_var <- paste(
    "Variance curve method",
    "Smooth empirical log(1+variance) along EVS rank",
    "The first two terminal d2 zero-crossings define the interval",
    "Earlier terminal zero = cutoff anchor",
    "Later terminal zero = terminal start",
    sep = "\n"
  )

  left_text_deriv <- paste(
    "Derivative method",
    "Zero-crossings are sign changes in smoothed d2",
    "No slope threshold is imposed",
    "No run-length rule is imposed",
    "Only the first two terminal d2 zeros are used",
    sep = "\n"
  )

  left_text_nb <- paste(
    "NB corroboration",
    "NB2-like side = all genes from cutoff anchor to the end",
    "Matched NB1 side = same number of genes immediately to the left",
    "These summaries corroborate the derivative-defined split",
    sep = "\n"
  )

  right_text_nb <- paste0(
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

  right_text_band <- paste0(
    "Final NB-supported summary",
    "\nCenter = ", cutoff_x,
    "\nInterval = [", interval_min_x, ", ", interval_max_x, "]",
    "\nTerminal start = ", terminal_x
  )

  p1 <- ggplot(df, aes(rank, abs_loading)) +
    annotate("rect", xmin = interval_min_x, xmax = interval_max_x, ymin = -Inf, ymax = Inf,
             fill = COLORS["interval_fill"], alpha = 0.18) +
    geom_line(color = COLORS["variance_fit"], linewidth = 0.95, na.rm = TRUE) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.85, color = COLORS["cutoff_anchor"]) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.85, color = COLORS["terminal_anchor"]) +
    geom_vline(xintercept = ref_x, linetype = 4, linewidth = 0.85, color = COLORS["ref5000"]) +
    geom_point(data = cutoff_pt_abs, aes(x = rank, y = abs_loading), inherit.aes = FALSE,
               shape = 16, size = 2.5, color = COLORS["cutoff_anchor"]) +
    geom_point(data = terminal_pt_abs, aes(x = rank, y = abs_loading), inherit.aes = FALSE,
               shape = 1, size = 2.8, stroke = 1.0, color = COLORS["terminal_anchor"]) +
    geom_point(data = ref_pt_abs, aes(x = rank, y = abs_loading), inherit.aes = FALSE,
               shape = 18, size = 2.8, color = COLORS["ref5000"]) +
    make_left_method_box(x_left, max(df$abs_loading, na.rm = TRUE), left_text_abs, size = 2.6) +
    make_right_numeric_box(x_right, max(df$abs_loading, na.rm = TRUE), right_text_abs, size = 2.55) +
    labs(
      title = paste0(title_prefix, ": absolute PC1 loading series"),
      subtitle = "Leading edge is on the RIGHT",
      x = "EVS rank",
      y = "|PC1 loading|"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.margin = margin(8, 34, 8, 10))

  p2 <- ggplot(df, aes(rank, var_fit)) +
    annotate("rect", xmin = interval_min_x, xmax = interval_max_x, ymin = -Inf, ymax = Inf,
             fill = COLORS["interval_fill"], alpha = 0.18) +
    geom_line(color = COLORS["variance_fit"], linewidth = 1.0, na.rm = TRUE) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.85, color = COLORS["cutoff_anchor"]) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.85, color = COLORS["terminal_anchor"]) +
    geom_vline(xintercept = ref_x, linetype = 4, linewidth = 0.85, color = COLORS["ref5000"]) +
    geom_point(data = cutoff_pt_var, aes(x = rank, y = var_fit), shape = 16, size = 2.6,
               inherit.aes = FALSE, color = COLORS["cutoff_anchor"]) +
    geom_point(data = terminal_pt_var, aes(x = rank, y = var_fit), shape = 1, size = 2.9, stroke = 1.0,
               inherit.aes = FALSE, color = COLORS["terminal_anchor"]) +
    geom_point(data = ref_pt_var, aes(x = rank, y = var_fit), shape = 18, size = 2.9,
               inherit.aes = FALSE, color = COLORS["ref5000"]) +
    make_left_method_box(x_left, max(df$var_fit, na.rm = TRUE), left_text_var, size = 2.45) +
    labs(
      title = paste0(title_prefix, ": smoothed empirical variance curve"),
      subtitle = "Selected geometric points are marked on the curve",
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.margin = margin(8, 34, 8, 10))

  deriv_df <- bind_rows(
    data.frame(rank = df$rank, value = df$d2_sm, series = "Smoothed curvature d2"),
    data.frame(rank = df$rank, value = df$d1_sm, series = "Smoothed slope d1")
  )

  p3 <- ggplot(deriv_df, aes(rank, value, color = series)) +
    annotate("rect", xmin = interval_min_x, xmax = interval_max_x, ymin = -Inf, ymax = Inf,
             fill = COLORS["interval_fill"], alpha = 0.18, inherit.aes = FALSE) +
    geom_line(linewidth = 0.95, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "Smoothed curvature d2" = COLORS["d2"],
        "Smoothed slope d1" = COLORS["d1"]
      )
    ) +
    geom_hline(yintercept = 0, linewidth = 0.55, color = "black") +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.85, color = COLORS["cutoff_anchor"]) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.85, color = COLORS["terminal_anchor"]) +
    geom_vline(xintercept = ref_x, linetype = 4, linewidth = 0.85, color = COLORS["ref5000"]) +
    geom_point(data = deriv_pt_cut, aes(x = rank, y = value), shape = 16, size = 2.6,
               inherit.aes = FALSE, color = COLORS["cutoff_anchor"]) +
    geom_point(data = deriv_pt_term, aes(x = rank, y = value), shape = 1, size = 2.9, stroke = 1.0,
               inherit.aes = FALSE, color = COLORS["terminal_anchor"]) +
    geom_point(data = deriv_pt_ref, aes(x = rank, y = value), shape = 18, size = 2.9,
               inherit.aes = FALSE, color = COLORS["ref5000"]) +
    make_left_method_box(x_left, max(deriv_df$value, na.rm = TRUE), left_text_deriv, size = 2.35) +
    labs(
      title = paste0(title_prefix, ": derivative support"),
      subtitle = "Cutoff anchor and terminal start come only from the selected terminal d2 zeros",
      x = "EVS rank",
      y = "Derivative value",
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          plot.margin = margin(8, 34, 8, 10))

  nb_df <- bind_rows(
    data.frame(rank = df$rank, value = df$nb_support, series = "Combined NB support"),
    data.frame(rank = df$rank, value = df$log_nb1, series = "NB1 = mu"),
    data.frame(rank = df$rank, value = df$log_nb2, series = "NB2 = variance - mu"),
    data.frame(rank = df$rank, value = df$amu_sm, series = "Smoothed log(alpha*mu)"),
    data.frame(rank = df$rank, value = df$nb_gap_sm, series = "Smoothed log(NB2+1) - log(NB1+1)")
  )

  p4 <- ggplot(nb_df, aes(rank, value, color = series)) +
    annotate("rect", xmin = interval_min_x, xmax = interval_max_x, ymin = -Inf, ymax = Inf,
             fill = COLORS["interval_fill"], alpha = 0.18, inherit.aes = FALSE) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_hline(yintercept = onset_info$nb_support_threshold, linetype = 3, linewidth = 0.55,
               color = "grey40") +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.85, color = COLORS["cutoff_anchor"]) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.85, color = COLORS["terminal_anchor"]) +
    geom_vline(xintercept = ref_x, linetype = 4, linewidth = 0.85, color = COLORS["ref5000"]) +
    geom_point(data = nb_pt_cut, aes(x = rank, y = value), shape = 16, size = 2.6,
               inherit.aes = FALSE, color = COLORS["cutoff_anchor"]) +
    geom_point(data = nb_pt_term, aes(x = rank, y = value), shape = 1, size = 2.9, stroke = 1.0,
               inherit.aes = FALSE, color = COLORS["terminal_anchor"]) +
    geom_point(data = nb_pt_ref, aes(x = rank, y = value), shape = 18, size = 2.9,
               inherit.aes = FALSE, color = COLORS["ref5000"]) +
    scale_color_manual(
      values = c(
        "Combined NB support" = COLORS["nb_support"],
        "NB1 = mu" = COLORS["nb1"],
        "NB2 = variance - mu" = COLORS["nb2"],
        "Smoothed log(alpha*mu)" = COLORS["amu"],
        "Smoothed log(NB2+1) - log(NB1+1)" = COLORS["nb_gap"]
      )
    ) +
    make_left_method_box(x_left, max(nb_df$value, na.rm = TRUE), left_text_nb, size = 2.2) +
    make_right_numeric_box(x_right, max(nb_df$value, na.rm = TRUE), right_text_nb, size = 2.0) +
    labs(
      title = paste0(title_prefix, ": NB1 / NB2 / alpha*mu support"),
      subtitle = "Right-of-cutoff elevation in NB2 and alpha*mu supports the leading-edge interpretation",
      x = "EVS rank",
      y = "Support value",
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          plot.margin = margin(8, 34, 8, 10))

  p5 <- ggplot(df, aes(rank, nb_support)) +
    annotate("rect", xmin = interval_min_x, xmax = interval_max_x, ymin = -Inf, ymax = Inf,
             fill = COLORS["interval_fill"], alpha = 0.18) +
    geom_line(color = COLORS["nb_support"], linewidth = 1.0, na.rm = TRUE) +
    geom_hline(yintercept = onset_info$nb_support_threshold, linetype = 3, linewidth = 0.55,
               color = "grey40") +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.85, color = COLORS["cutoff_anchor"]) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.85, color = COLORS["terminal_anchor"]) +
    geom_vline(xintercept = ref_x, linetype = 4, linewidth = 0.85, color = COLORS["ref5000"]) +
    geom_point(data = nb_pt_cut, aes(x = rank, y = value), shape = 16, size = 2.6,
               inherit.aes = FALSE, color = COLORS["cutoff_anchor"]) +
    geom_point(data = nb_pt_term, aes(x = rank, y = value), shape = 1, size = 2.9, stroke = 1.0,
               inherit.aes = FALSE, color = COLORS["terminal_anchor"]) +
    geom_point(data = nb_pt_ref, aes(x = rank, y = value), shape = 18, size = 2.9,
               inherit.aes = FALSE, color = COLORS["ref5000"]) +
    make_left_method_box(
      x_left,
      max(df$nb_support, na.rm = TRUE),
      paste(
        "Final NB support panel",
        "The shaded interval is the final manuscript study interval",
        "Combined NB support is shown for corroboration",
        sep = "\n"
      ),
      size = 2.35
    ) +
    make_right_numeric_box(x_right, max(df$nb_support, na.rm = TRUE), right_text_band, size = 2.35) +
    labs(
      title = paste0(title_prefix, ": combined NB support"),
      subtitle = "The geometric interval remains primary and NB support is corroborative",
      x = "EVS rank",
      y = "Combined NB support"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.margin = margin(8, 34, 8, 10))

  png(out_file, width = 3000, height = 3900, res = 240)
  gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol = 1)
  dev.off()
}

build_range_panel <- function(ctrl_onset, trt_onset, title_prefix, out_file) {
  ctrl_df <- ctrl_onset$curve_df
  trt_df <- trt_onset$curve_df

  x_all <- c(ctrl_df$rank, trt_df$rank)
  y_all <- c(ctrl_df$var_fit, trt_df$var_fit)

  x_left <- min(x_all, na.rm = TRUE) + 0.05 * diff(range(x_all, na.rm = TRUE))
  x_right <- min(x_all, na.rm = TRUE) + 0.68 * diff(range(x_all, na.rm = TRUE))

  p <- ggplot() +
    annotate(
      "rect",
      xmin = min(ctrl_onset$cutoff_range_rank_min, trt_onset$cutoff_range_rank_min),
      xmax = max(ctrl_onset$cutoff_range_rank_max, trt_onset$cutoff_range_rank_max),
      ymin = -Inf, ymax = Inf,
      fill = COLORS["interval_fill"], alpha = 0.12
    ) +
    geom_line(data = ctrl_df, aes(rank, var_fit, color = "Control variance fit"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(data = trt_df, aes(rank, var_fit, color = "Treatment variance fit"), linewidth = 1.0, na.rm = TRUE) +
    geom_vline(xintercept = ctrl_onset$cutoff_center_rank, linetype = 2, linewidth = 0.85, color = "#252525") +
    geom_vline(xintercept = trt_onset$cutoff_center_rank, linetype = 3, linewidth = 0.85, color = "#636363") +
    geom_vline(xintercept = ctrl_onset$terminal_start_rank, linetype = 2, linewidth = 0.55, color = "#252525") +
    geom_vline(xintercept = trt_onset$terminal_start_rank, linetype = 3, linewidth = 0.55, color = "#636363") +
    scale_color_manual(
      values = c(
        "Control variance fit" = "#1b9e77",
        "Treatment variance fit" = "#7570b3"
      )
    ) +
    make_left_method_box(
      x_left,
      max(y_all, na.rm = TRUE),
      paste(
        "Treatment/control comparison",
        "Both arms use the same manuscript rule",
        "Earlier terminal d2 zero = cutoff anchor",
        "Later terminal d2 zero = terminal start",
        sep = "\n"
      ),
      size = 2.8
    ) +
    make_right_numeric_box(
      x_right,
      max(y_all, na.rm = TRUE),
      paste0(
        "Control interval = [", ctrl_onset$cutoff_range_rank_min, ", ", ctrl_onset$cutoff_range_rank_max, "]",
        "\nTreatment interval = [", trt_onset$cutoff_range_rank_min, ", ", trt_onset$cutoff_range_rank_max, "]",
        "\nControl cutoff = ", ctrl_onset$cutoff_center_rank,
        "\nTreatment cutoff = ", trt_onset$cutoff_center_rank
      ),
      size = 2.7
    ) +
    labs(
      title = paste0(title_prefix, ": treatment/control derivative-defined intervals"),
      subtitle = "Both arms are shown on the same variance-fit axis",
      x = "EVS rank",
      y = "Fitted log(1 + variance)",
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          plot.margin = margin(8, 34, 8, 10))

  png(out_file, width = 3000, height = 1400, res = 240)
  print(p)
  dev.off()
}

# =============================================================================
# MAIN ANALYSIS LOOP
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# Each comparison is processed independently. Control and treatment arms are
# ranked separately, their own variance geometry is fitted separately, and
# their own geometric anchors are selected separately. Per-comparison feature
# tables, full zero-crossing tables, selected anchor tables, figure panels,
# and manuscript summary tables are written to disk.
# =============================================================================

message("SEQUENCE.R started"); flush.console()
message("Resolving count file..."); flush.console()
count_file <- resolve_counts_file(count_file)
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
  kept_ids <- rownames(cmp_counts)

  evs_tbl <- build_evs_table(cmp_counts, comp$ctrl_ids, comp$trt_ids, kept_ids)
  metrics_tbl <- compute_feature_metrics_empirical(cmp_counts, comp$ctrl_ids, comp$trt_ids, kept_ids)

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
  ctrl_onset <- select_cutoff_derivative_nb_range(
    ctrl_rank_df,
    spline_spar = spline_spar,
    terminal_fraction_start = terminal_fraction_start,
    zero_tol = zero_tol,
    nb_smooth_window = nb_smooth_window,
    fixed_rank_reference = fixed_rank_reference
  )
  message(
    "Control cutoff anchor rank: ", ctrl_onset$cutoff_center_rank,
    " | terminal anchor: ", ctrl_onset$terminal_start_rank,
    " | interval: [", ctrl_onset$cutoff_range_rank_min, ", ", ctrl_onset$cutoff_range_rank_max, "]",
    " | pre-EVS remainder: ", ctrl_onset$pre_evs_remainder_size,
    " | pre-EVS leading edge: ", ctrl_onset$pre_evs_leading_edge_size
  ); flush.console()

  message("Selecting treatment cutoff for ", cmp_name); flush.console()
  trt_onset <- select_cutoff_derivative_nb_range(
    trt_rank_df,
    spline_spar = spline_spar,
    terminal_fraction_start = terminal_fraction_start,
    zero_tol = zero_tol,
    nb_smooth_window = nb_smooth_window,
    fixed_rank_reference = fixed_rank_reference
  )
  message(
    "Treatment cutoff anchor rank: ", trt_onset$cutoff_center_rank,
    " | terminal anchor: ", trt_onset$terminal_start_rank,
    " | interval: [", trt_onset$cutoff_range_rank_min, ", ", trt_onset$cutoff_range_rank_max, "]",
    " | pre-EVS remainder: ", trt_onset$pre_evs_remainder_size,
    " | pre-EVS leading edge: ", trt_onset$pre_evs_leading_edge_size
  ); flush.console()

  build_dataset_panel(
    ctrl_rank_df,
    ctrl_onset,
    paste0(cmp_name, " control"),
    file.path(cmp_dir, paste0(cmp_name, "_control_rank_panel.png"))
  )
  build_dataset_panel(
    trt_rank_df,
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

  if (nrow(ctrl_onset$zero_crossings) > 0L) {
    utils::write.csv(
      ctrl_onset$zero_crossings,
      file.path(cmp_dir, paste0(cmp_name, "_control_zero_crossings_all.csv")),
      row.names = FALSE
    )
  }

  if (nrow(trt_onset$zero_crossings) > 0L) {
    utils::write.csv(
      trt_onset$zero_crossings,
      file.path(cmp_dir, paste0(cmp_name, "_treatment_zero_crossings_all.csv")),
      row.names = FALSE
    )
  }

  if (nrow(ctrl_onset$selected_zero_crossings) > 0L) {
    utils::write.csv(
      ctrl_onset$selected_zero_crossings,
      file.path(cmp_dir, paste0(cmp_name, "_control_selected_two_zero_crossings.csv")),
      row.names = FALSE
    )
  }

  if (nrow(trt_onset$selected_zero_crossings) > 0L) {
    utils::write.csv(
      trt_onset$selected_zero_crossings,
      file.path(cmp_dir, paste0(cmp_name, "_treatment_selected_two_zero_crossings.csv")),
      row.names = FALSE
    )
  }

  ctrl_corr <- ctrl_onset$corrob
  trt_corr <- trt_onset$corrob

  cutoff_summary <- tibble(
    comparison = cmp_name,
    total_features = nrow(ctrl_rank_df),
    fixed_rank_reference = fixed_rank_reference,

    control_mode = ctrl_onset$mode,
    control_cutoff_anchor_rank = ctrl_onset$cutoff_center_rank,
    control_terminal_anchor_rank = ctrl_onset$terminal_start_rank,
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
    control_right_side_more_nb2 = ctrl_corr$right_side_more_nb2,
    control_right_side_more_nb_gap = ctrl_corr$right_side_more_nb_gap,
    control_right_side_more_alpha_mu = ctrl_corr$right_side_more_alpha_mu,

    treatment_mode = trt_onset$mode,
    treatment_cutoff_anchor_rank = trt_onset$cutoff_center_rank,
    treatment_terminal_anchor_rank = trt_onset$terminal_start_rank,
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
    treatment_right_left_log_alpha_mu_diff = trt_corr$right_left_log_alpha_mu_diff,
    treatment_right_side_more_nb2 = trt_corr$right_side_more_nb2,
    treatment_right_side_more_nb_gap = trt_corr$right_side_more_nb_gap,
    treatment_right_side_more_alpha_mu = trt_corr$right_side_more_alpha_mu
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
