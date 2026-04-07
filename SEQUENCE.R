# =============================================================================
# SEQUENCE STAGE 1: DERIVATIVE-DEFINED PRE-TERMINAL CUTOFF WITH NB1/NB2 RANGE
# -----------------------------------------------------------------------------
# MANUSCRIPT-READY FINAL VERSION
#
# Purpose
#   This script identifies a pre-EVS split point along the EVS rank axis using
#   the smoothed empirical variance curve and its derivatives, then uses NB1,
#   NB2, and alpha*mu support to define a cutoff range around that derivative-
#   defined center.
#
# Orientation
#   - Left  = lowest absolute loading
#   - Right = highest absolute loading
#   - Leading edge is on the RIGHT
#
# Empirical quantities
#   mu        = mean(raw counts)
#   variance  = var(raw counts)
#   NB1       = mu
#   NB2       = variance - mu
#   alpha     = (variance - mu) / mu^2
#   alpha*mu  = (variance - mu) / mu
#
# Final selector
#   1. Smooth the empirical variance curve along EVS rank.
#   2. Compute first derivative d1 and second derivative d2.
#   3. Define terminal_start as the d2 zero-crossing immediately preceding the
#      terminal sustained positive-slope rise.
#   4. Define cutoff_center as the nearest earlier d2 zero-crossing immediately
#      preceding terminal_start and associated with decreasing-slope behavior.
#   5. Define an NB-supported cutoff range around cutoff_center using:
#        a. smoothed log(NB2 + 1) - log(NB1 + 1)
#        b. smoothed log(alpha*mu)
#      The range is the contiguous region around cutoff_center where the local
#      NB support remains elevated relative to the center.
#
# Reporting
#   cutoff_center_rank = derivative-defined center
#   cutoff_range_rank_min / cutoff_range_rank_max = NB-supported transition band
#   pre_evs_remainder_size = cutoff_center_rank - 1
#   pre_evs_leading_edge_size = N - cutoff_center_rank + 1
#
# Outputs
#   exports/variance_derivative_nb_range/
#     <comparison>_cutoff_folder/
#       <comparison>_control_rank_panel.png
#       <comparison>_treatment_rank_panel.png
#       <comparison>_cutoff_range_panel.png
#       <comparison>_feature_level_metrics.csv
#       <comparison>_control_rank_series.csv
#       <comparison>_treatment_rank_series.csv
#       <comparison>_control_zero_crossings.csv
#       <comparison>_treatment_zero_crossings.csv
#       <comparison>_cutoff_summary.csv
#     overall_cutoff_summary.csv
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

options(warn = 1)

# =============================================================================
# USER SETTINGS
# =============================================================================

repo_dir <- getwd()
count_file <- file.path(repo_dir, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
out_root <- file.path(repo_dir, "exports", "variance_derivative_nb_range")
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

# smoothing / search
spline_spar <- 0.60
left_edge_buffer <- 50L
right_edge_buffer <- 20L
search_fraction_min <- 0.20
search_fraction_max <- 0.985

# terminal-rise detection
terminal_run_fraction <- 0.20
terminal_run_min_length <- 250L
terminal_positive_slope_quantile <- 0.70

# NB range around cutoff center
nb_smooth_window <- 151L
nb_range_drop_fraction <- 0.70
nb_range_max_span_fraction <- 0.18

# =============================================================================
# HELPERS
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

find_runs <- function(cond) {
  cond[is.na(cond)] <- FALSE
  r <- rle(cond)
  ends <- cumsum(r$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  data.frame(start = starts, end = ends, value = r$values, length = r$lengths)
}

compute_group_pc1_loadings <- function(count_mat, group_cols) {
  mat <- count_mat[, group_cols, drop = FALSE]
  mat <- log2(mat + 1)
  mat <- t(mat)

  if (nrow(mat) < 2L) {
    return(rep(NA_real_, ncol(mat)))
  }

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

  out$log_variance <- NA_real_
  idx <- is.finite(out$variance) & out$variance >= 0
  out$log_variance[idx] <- log1p(out$variance[idx])

  out$log_nb1 <- NA_real_
  idx <- is.finite(out$nb1) & out$nb1 >= 0
  out$log_nb1[idx] <- log1p(out$nb1[idx])

  out$log_nb2 <- NA_real_
  idx <- is.finite(out$nb2) & out$nb2 >= 0
  out$log_nb2[idx] <- log1p(out$nb2[idx])

  out$nb_gap <- out$log_nb2 - out$log_nb1

  out$log_alpha_mu <- NA_real_
  idx <- is.finite(out$alpha_mu) & out$alpha_mu > 0
  out$log_alpha_mu[idx] <- log(out$alpha_mu[idx])

  out
}

find_d2_zero_crossings <- function(d2_vec) {
  d2_vec <- as.numeric(d2_vec)
  n <- length(d2_vec)
  out <- integer(0)

  if (n < 2L) return(out)

  s <- sign(d2_vec)
  s[!is.finite(s)] <- 0

  for (i in 2:n) {
    s0 <- s[i - 1L]
    s1 <- s[i]

    if (!is.finite(d2_vec[i - 1L]) || !is.finite(d2_vec[i])) next

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

summarize_zero_crossings <- function(df, zero_idx, terminal_start_idx, terminal_end_idx) {
  if (!length(zero_idx)) {
    return(data.frame())
  }

  tbl <- lapply(zero_idx, function(i) {
    left_d1 <- if (i > 1L) df$d1_sm[i - 1L] else NA_real_
    right_d1 <- if (i < nrow(df)) df$d1_sm[i + 1L] else NA_real_
    left_d2 <- if (i > 1L) df$d2_sm[i - 1L] else NA_real_
    right_d2 <- if (i < nrow(df)) df$d2_sm[i + 1L] else NA_real_

    data.frame(
      index = i,
      rank = df$rank[i],
      d1 = df$d1_sm[i],
      d2 = df$d2_sm[i],
      left_d1 = left_d1,
      right_d1 = right_d1,
      left_d2 = left_d2,
      right_d2 = right_d2,
      before_terminal_start = i < terminal_start_idx,
      inside_terminal_block = i >= terminal_start_idx && i <= terminal_end_idx,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(tbl)
}

find_terminal_block <- function(d1_sm, search_start, search_end,
                                terminal_run_fraction = 0.20,
                                terminal_run_min_length = 250L,
                                terminal_positive_slope_quantile = 0.70) {
  n <- length(d1_sm)
  tail_start <- max(search_start, floor((1 - terminal_run_fraction) * n))
  tail_idx <- seq.int(tail_start, search_end)
  finite_tail <- is.finite(d1_sm[tail_idx])

  if (!any(finite_tail)) {
    return(list(start = tail_start, end = search_end))
  }

  pos_cut <- as.numeric(stats::quantile(
    d1_sm[tail_idx][finite_tail],
    probs = terminal_positive_slope_quantile,
    na.rm = TRUE,
    names = FALSE
  ))

  if (!is.finite(pos_cut)) {
    pos_cut <- stats::median(d1_sm[tail_idx][finite_tail], na.rm = TRUE)
  }

  terminal_cond <- rep(FALSE, n)
  terminal_cond[tail_idx] <- is.finite(d1_sm[tail_idx]) & (d1_sm[tail_idx] >= pos_cut)

  run_tbl <- find_runs(terminal_cond[tail_idx])
  good_runs <- run_tbl[run_tbl$value & run_tbl$length >= terminal_run_min_length, , drop = FALSE]

  if (nrow(good_runs)) {
    terminal_local_start <- good_runs$start[1]
    terminal_local_end <- good_runs$end[1]
    terminal_start <- tail_start + terminal_local_start - 1L
    terminal_end <- tail_start + terminal_local_end - 1L
  } else {
    best_local <- which.max(replace(d1_sm[tail_idx], !is.finite(d1_sm[tail_idx]), -Inf))
    if (!length(best_local) || !is.finite(best_local)) best_local <- 1L
    terminal_start <- tail_idx[max(1L, best_local - terminal_run_min_length + 1L)]
    terminal_end <- tail_idx[min(length(tail_idx), best_local + terminal_run_min_length - 1L)]
  }

  list(start = terminal_start, end = terminal_end)
}

choose_terminal_start_from_geometry <- function(df, search_start, search_end,
                                                terminal_run_fraction = 0.20,
                                                terminal_run_min_length = 250L,
                                                terminal_positive_slope_quantile = 0.70) {
  block <- find_terminal_block(
    d1_sm = df$d1_sm,
    search_start = search_start,
    search_end = search_end,
    terminal_run_fraction = terminal_run_fraction,
    terminal_run_min_length = terminal_run_min_length,
    terminal_positive_slope_quantile = terminal_positive_slope_quantile
  )

  zero_idx <- find_d2_zero_crossings(df$d2_sm)
  ztbl <- summarize_zero_crossings(df, zero_idx, block$start, block$end)

  if (nrow(ztbl)) {
    candidates <- ztbl %>%
      filter(index <= block$start) %>%
      filter(is.finite(right_d1)) %>%
      filter(right_d1 > 0)

    if (nrow(candidates)) {
      terminal_idx <- max(candidates$index)
    } else {
      terminal_idx <- block$start
    }
  } else {
    terminal_idx <- block$start
  }

  list(
    terminal_start_index = terminal_idx,
    terminal_start_rank = df$rank[terminal_idx],
    terminal_end_index = block$end,
    terminal_end_rank = df$rank[block$end],
    zero_crossings = ztbl
  )
}

choose_cutoff_center_from_geometry <- function(df, terminal_start_idx, zero_tbl) {
  if (nrow(zero_tbl) == 0L) {
    idx <- max(1L, terminal_start_idx - 1L)
    return(list(
      cutoff_center_index = idx,
      cutoff_center_rank = df$rank[idx],
      mode = "fallback: terminal_start minus one"
    ))
  }

  prior_tbl <- zero_tbl %>%
    filter(index < terminal_start_idx)

  if (!nrow(prior_tbl)) {
    idx <- max(1L, terminal_start_idx - 1L)
    return(list(
      cutoff_center_index = idx,
      cutoff_center_rank = df$rank[idx],
      mode = "fallback: no earlier d2 zero-crossing"
    ))
  }

  # Prior to terminal start, select the last d2 zero-crossing whose local slope
  # behavior is consistent with a decreasing-slope turning point immediately
  # before the terminal-start transition.
  #
  # Prefer crossings with right-side slope not strongly positive and with
  # negative or weak local slope neighborhood, then choose the nearest one.
  cand1 <- prior_tbl %>%
    filter(is.finite(right_d1)) %>%
    filter(right_d1 <= 0)

  if (nrow(cand1)) {
    idx <- max(cand1$index)
    return(list(
      cutoff_center_index = idx,
      cutoff_center_rank = df$rank[idx],
      mode = "last earlier d2 zero-crossing before terminal start with nonpositive right-side slope"
    ))
  }

  cand2 <- prior_tbl %>%
    filter(is.finite(d1)) %>%
    filter(d1 <= 0)

  if (nrow(cand2)) {
    idx <- max(cand2$index)
    return(list(
      cutoff_center_index = idx,
      cutoff_center_rank = df$rank[idx],
      mode = "last earlier d2 zero-crossing before terminal start with nonpositive local slope"
    ))
  }

  idx <- max(prior_tbl$index)
  list(
    cutoff_center_index = idx,
    cutoff_center_rank = df$rank[idx],
    mode = "last earlier d2 zero-crossing before terminal start"
  )
}

compute_nb_support_score <- function(df) {
  nb_gap_sm <- roll_median(df$nb_gap, nb_smooth_window)
  amu_sm <- roll_median(df$log_alpha_mu, nb_smooth_window)

  nb_gap_s <- scale01(nb_gap_sm)
  amu_s <- scale01(amu_sm)

  nb_gap_s[!is.finite(nb_gap_s)] <- 0
  amu_s[!is.finite(amu_s)] <- 0

  nb_support <- 0.5 * nb_gap_s + 0.5 * amu_s

  list(
    nb_gap_sm = nb_gap_sm,
    amu_sm = amu_sm,
    nb_support = nb_support
  )
}

expand_nb_range_around_center <- function(df, center_idx, terminal_start_idx,
                                          nb_support,
                                          nb_range_drop_fraction = 0.70,
                                          nb_range_max_span_fraction = 0.18) {
  n <- nrow(df)
  max_span <- max(100L, floor(nb_range_max_span_fraction * n))

  center_score <- nb_support[center_idx]
  if (!is.finite(center_score)) center_score <- 0

  threshold <- nb_range_drop_fraction * center_score

  left_limit <- max(1L, center_idx - max_span)
  right_limit <- min(terminal_start_idx, center_idx + max_span)

  left_idx <- center_idx
  while (left_idx > left_limit) {
    test_idx <- left_idx - 1L
    if (!is.finite(nb_support[test_idx])) break
    if (nb_support[test_idx] < threshold) break
    left_idx <- test_idx
  }

  right_idx <- center_idx
  while (right_idx < right_limit) {
    test_idx <- right_idx + 1L
    if (!is.finite(nb_support[test_idx])) break
    if (nb_support[test_idx] < threshold) break
    right_idx <- test_idx
  }

  list(
    range_left_index = left_idx,
    range_left_rank = df$rank[left_idx],
    range_right_index = right_idx,
    range_right_rank = df$rank[right_idx],
    center_nb_support = center_score,
    nb_support_threshold = threshold
  )
}

select_cutoff_derivative_nb_range <- function(rank_df,
                                              spline_spar = 0.60,
                                              search_fraction_min = 0.20,
                                              search_fraction_max = 0.985,
                                              left_edge_buffer = 50L,
                                              right_edge_buffer = 20L,
                                              terminal_run_fraction = 0.20,
                                              terminal_run_min_length = 250L,
                                              terminal_positive_slope_quantile = 0.70,
                                              nb_range_drop_fraction = 0.70,
                                              nb_range_max_span_fraction = 0.18) {
  df <- rank_df
  n <- nrow(df)
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

  var_fit <- rep(NA_real_, n)
  d1_sm <- rep(NA_real_, n)
  d2_sm <- rep(NA_real_, n)

  var_fit[ok] <- pred0$y
  d1_sm[ok] <- pred1$y
  d2_sm[ok] <- pred2$y

  df$var_fit <- var_fit
  df$d1_sm <- d1_sm
  df$d2_sm <- d2_sm

  search_start <- max(left_edge_buffer + 1L, floor(search_fraction_min * n))
  search_end <- min(n - right_edge_buffer, floor(search_fraction_max * n))
  if (search_end <= search_start + 10L) {
    search_start <- max(1L, left_edge_buffer + 1L)
    search_end <- min(n - right_edge_buffer, n - 1L)
  }

  term <- choose_terminal_start_from_geometry(
    df = df,
    search_start = search_start,
    search_end = search_end,
    terminal_run_fraction = terminal_run_fraction,
    terminal_run_min_length = terminal_run_min_length,
    terminal_positive_slope_quantile = terminal_positive_slope_quantile
  )

  cut <- choose_cutoff_center_from_geometry(
    df = df,
    terminal_start_idx = term$terminal_start_index,
    zero_tbl = term$zero_crossings
  )

  nb_obj <- compute_nb_support_score(df)
  df$nb_gap_sm <- nb_obj$nb_gap_sm
  df$amu_sm <- nb_obj$amu_sm
  df$nb_support <- nb_obj$nb_support

  nb_rng <- expand_nb_range_around_center(
    df = df,
    center_idx = cut$cutoff_center_index,
    terminal_start_idx = term$terminal_start_index,
    nb_support = df$nb_support,
    nb_range_drop_fraction = nb_range_drop_fraction,
    nb_range_max_span_fraction = nb_range_max_span_fraction
  )

  total_features <- n
  pre_evs_remainder_size <- cut$cutoff_center_index - 1L
  pre_evs_leading_edge_size <- total_features - cut$cutoff_center_index + 1L

  list(
    curve_df = df,
    zero_crossings = term$zero_crossings,
    mode = cut$mode,

    terminal_start_index = term$terminal_start_index,
    terminal_start_rank = term$terminal_start_rank,
    terminal_end_index = term$terminal_end_index,
    terminal_end_rank = term$terminal_end_rank,

    cutoff_center_index = cut$cutoff_center_index,
    cutoff_center_rank = cut$cutoff_center_rank,

    cutoff_range_index_min = nb_rng$range_left_index,
    cutoff_range_rank_min = nb_rng$range_left_rank,
    cutoff_range_index_max = nb_rng$range_right_index,
    cutoff_range_rank_max = nb_rng$range_right_rank,

    center_nb_support = nb_rng$center_nb_support,
    nb_support_threshold = nb_rng$nb_support_threshold,

    total_features = total_features,
    pre_evs_remainder_size = pre_evs_remainder_size,
    pre_evs_leading_edge_size = pre_evs_leading_edge_size
  )
}

build_dataset_panel <- function(rank_df, onset_info, title_prefix, out_file) {
  df <- onset_info$curve_df
  center_x <- onset_info$cutoff_center_rank
  range_min_x <- onset_info$cutoff_range_rank_min
  range_max_x <- onset_info$cutoff_range_rank_max
  terminal_x <- onset_info$terminal_start_rank
  terminal_end_x <- onset_info$terminal_end_rank

  ztbl <- onset_info$zero_crossings

  label_text <- paste0(
    "Derivative-defined pre-terminal selector",
    "\nTerminal start = d2 zero-crossing before sustained increasing slope",
    "\nCutoff center = earlier d2 zero-crossing before terminal start",
    "\nNB range = region around cutoff center with elevated NB2 / alpha*mu support",
    "\nCutoff center rank = ", center_x,
    "\nCutoff range = [", range_min_x, ", ", range_max_x, "]",
    "\nTerminal start rank = ", terminal_x,
    "\nPre-EVS remainder = ", onset_info$pre_evs_remainder_size,
    "\nPre-EVS leading edge = ", onset_info$pre_evs_leading_edge_size,
    "\nCenter NB support = ", round(onset_info$center_nb_support, 3),
    "\nNB threshold = ", round(onset_info$nb_support_threshold, 3)
  )

  p1 <- ggplot(df, aes(rank, abs_loading)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    annotate(
      "rect",
      xmin = center_x,
      xmax = max(df$rank, na.rm = TRUE),
      ymin = -Inf, ymax = Inf,
      alpha = 0.04
    ) +
    annotate(
      "rect",
      xmin = range_min_x,
      xmax = range_max_x,
      ymin = -Inf, ymax = Inf,
      alpha = 0.08
    ) +
    annotate(
      "rect",
      xmin = terminal_x,
      xmax = terminal_end_x,
      ymin = -Inf, ymax = Inf,
      alpha = 0.08
    ) +
    geom_vline(xintercept = center_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    annotate(
      "label",
      x = center_x,
      y = max(df$abs_loading, na.rm = TRUE),
      label = label_text,
      hjust = 0, vjust = 1, size = 2.9
    ) +
    labs(
      title = paste0(title_prefix, ": absolute PC1 loading series"),
      subtitle = "Left = lowest loading | Right = highest loading (leading edge)",
      x = "EVS rank",
      y = "|PC1 loading|"
    ) +
    theme_bw(base_size = 10)

  p2 <- ggplot(df, aes(rank, var_fit)) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    annotate(
      "rect",
      xmin = range_min_x,
      xmax = range_max_x,
      ymin = -Inf, ymax = Inf,
      alpha = 0.08
    ) +
    annotate(
      "rect",
      xmin = terminal_x,
      xmax = terminal_end_x,
      ymin = -Inf, ymax = Inf,
      alpha = 0.08
    ) +
    geom_vline(xintercept = center_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    geom_point(
      data = df[df$rank == center_x & is.finite(df$var_fit), , drop = FALSE],
      aes(x = rank, y = var_fit),
      size = 2.4,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = df[df$rank == terminal_x & is.finite(df$var_fit), , drop = FALSE],
      aes(x = rank, y = var_fit),
      size = 2.0,
      shape = 1,
      inherit.aes = FALSE
    ) +
    labs(
      title = paste0(title_prefix, ": smoothed empirical variance curve"),
      subtitle = paste0(
        "Center = earlier d2 zero-crossing before terminal start | ",
        "Terminal start = d2 zero-crossing before sustained positive slope"
      ),
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    theme_bw(base_size = 10)

  if (nrow(ztbl) > 0L) {
    p2 <- p2 +
      geom_vline(
        data = ztbl,
        aes(xintercept = rank),
        linetype = 3,
        linewidth = 0.35,
        alpha = 0.45,
        inherit.aes = FALSE
      )
  }

  p3 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = d1_sm, color = "Smoothed slope d1"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(aes(y = d2_sm, color = "Smoothed curvature d2"), linewidth = 1.0, na.rm = TRUE) +
    annotate(
      "rect",
      xmin = range_min_x,
      xmax = range_max_x,
      ymin = -Inf, ymax = Inf,
      alpha = 0.08
    ) +
    annotate(
      "rect",
      xmin = terminal_x,
      xmax = terminal_end_x,
      ymin = -Inf, ymax = Inf,
      alpha = 0.08
    ) +
    geom_vline(xintercept = center_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(
      title = paste0(title_prefix, ": derivative support"),
      subtitle = "Zero-crossings of d2 define terminal start and cutoff center",
      x = "EVS rank",
      y = "Derivative value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p4 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = log_nb1, color = "NB1 = mu"), linewidth = 0.5, alpha = 0.35, na.rm = TRUE) +
    geom_line(aes(y = log_nb2, color = "NB2 = variance - mu"), linewidth = 0.5, alpha = 0.35, na.rm = TRUE) +
    geom_line(aes(y = nb_gap_sm, color = "Smoothed log(NB2+1) - log(NB1+1)"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(aes(y = amu_sm, color = "Smoothed log(alpha*mu)"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(aes(y = nb_support, color = "Combined NB support"), linewidth = 1.0, na.rm = TRUE) +
    annotate(
      "rect",
      xmin = range_min_x,
      xmax = range_max_x,
      ymin = -Inf, ymax = Inf,
      alpha = 0.08
    ) +
    geom_vline(xintercept = center_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    geom_hline(yintercept = onset_info$nb_support_threshold, linetype = 3, linewidth = 0.5) +
    labs(
      title = paste0(title_prefix, ": NB1 / NB2 / alpha*mu support"),
      subtitle = "NB-supported cutoff range is centered on derivative-defined cutoff center",
      x = "EVS rank",
      y = "Support value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p5 <- ggplot(df, aes(rank, nb_support)) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    annotate(
      "rect",
      xmin = range_min_x,
      xmax = range_max_x,
      ymin = -Inf, ymax = Inf,
      alpha = 0.08
    ) +
    geom_vline(xintercept = center_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    geom_hline(yintercept = onset_info$nb_support_threshold, linetype = 3, linewidth = 0.5) +
    labs(
      title = paste0(title_prefix, ": NB-supported cutoff band"),
      subtitle = "Range expands left and right from cutoff center while combined NB support stays elevated",
      x = "EVS rank",
      y = "Combined NB support"
    ) +
    theme_bw(base_size = 10)

  png(out_file, width = 2200, height = 3000, res = 200)
  gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol = 1)
  dev.off()
}

build_range_panel <- function(ctrl_onset, trt_onset, title_prefix, out_file) {
  ctrl_df <- ctrl_onset$curve_df
  trt_df <- trt_onset$curve_df

  p <- ggplot() +
    geom_line(data = ctrl_df, aes(rank, var_fit, color = "Control variance fit"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(data = trt_df, aes(rank, var_fit, color = "Treatment variance fit"), linewidth = 1.0, na.rm = TRUE) +
    annotate(
      "rect",
      xmin = min(ctrl_onset$cutoff_range_rank_min, trt_onset$cutoff_range_rank_min),
      xmax = max(ctrl_onset$cutoff_range_rank_max, trt_onset$cutoff_range_rank_max),
      ymin = -Inf, ymax = Inf,
      alpha = 0.08
    ) +
    geom_vline(xintercept = ctrl_onset$cutoff_center_rank, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = trt_onset$cutoff_center_rank, linetype = 3, linewidth = 0.8) +
    labs(
      title = paste0(title_prefix, ": treatment/control derivative-defined centers and NB-supported ranges"),
      subtitle = paste0(
        "Control center = ", ctrl_onset$cutoff_center_rank,
        " | Treatment center = ", trt_onset$cutoff_center_rank
      ),
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  png(out_file, width = 2200, height = 1200, res = 200)
  print(p)
  dev.off()
}

# =============================================================================
# MAIN
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
    search_fraction_min = search_fraction_min,
    search_fraction_max = search_fraction_max,
    left_edge_buffer = left_edge_buffer,
    right_edge_buffer = right_edge_buffer,
    terminal_run_fraction = terminal_run_fraction,
    terminal_run_min_length = terminal_run_min_length,
    terminal_positive_slope_quantile = terminal_positive_slope_quantile,
    nb_range_drop_fraction = nb_range_drop_fraction,
    nb_range_max_span_fraction = nb_range_max_span_fraction
  )
  message(
    "Control cutoff center rank: ", ctrl_onset$cutoff_center_rank,
    " | terminal start: ", ctrl_onset$terminal_start_rank,
    " | cutoff range: [", ctrl_onset$cutoff_range_rank_min, ", ", ctrl_onset$cutoff_range_rank_max, "]",
    " | pre-EVS remainder: ", ctrl_onset$pre_evs_remainder_size,
    " | pre-EVS leading edge: ", ctrl_onset$pre_evs_leading_edge_size
  ); flush.console()

  message("Selecting treatment cutoff for ", cmp_name); flush.console()
  trt_onset <- select_cutoff_derivative_nb_range(
    trt_rank_df,
    spline_spar = spline_spar,
    search_fraction_min = search_fraction_min,
    search_fraction_max = search_fraction_max,
    left_edge_buffer = left_edge_buffer,
    right_edge_buffer = right_edge_buffer,
    terminal_run_fraction = terminal_run_fraction,
    terminal_run_min_length = terminal_run_min_length,
    terminal_positive_slope_quantile = terminal_positive_slope_quantile,
    nb_range_drop_fraction = nb_range_drop_fraction,
    nb_range_max_span_fraction = nb_range_max_span_fraction
  )
  message(
    "Treatment cutoff center rank: ", trt_onset$cutoff_center_rank,
    " | terminal start: ", trt_onset$terminal_start_rank,
    " | cutoff range: [", trt_onset$cutoff_range_rank_min, ", ", trt_onset$cutoff_range_rank_max, "]",
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
      file.path(cmp_dir, paste0(cmp_name, "_control_zero_crossings.csv")),
      row.names = FALSE
    )
  }

  if (nrow(trt_onset$zero_crossings) > 0L) {
    utils::write.csv(
      trt_onset$zero_crossings,
      file.path(cmp_dir, paste0(cmp_name, "_treatment_zero_crossings.csv")),
      row.names = FALSE
    )
  }

  cutoff_summary <- tibble(
    comparison = cmp_name,
    total_features = nrow(ctrl_rank_df),

    control_cutoff_mode = ctrl_onset$mode,
    control_cutoff_center_rank = ctrl_onset$cutoff_center_rank,
    control_cutoff_fraction = ctrl_onset$cutoff_center_rank / nrow(ctrl_rank_df),
    control_cutoff_range_rank_min = ctrl_onset$cutoff_range_rank_min,
    control_cutoff_range_rank_max = ctrl_onset$cutoff_range_rank_max,
    control_terminal_start_rank = ctrl_onset$terminal_start_rank,
    control_terminal_end_rank = ctrl_onset$terminal_end_rank,
    control_center_nb_support = ctrl_onset$center_nb_support,
    control_nb_support_threshold = ctrl_onset$nb_support_threshold,
    control_pre_evs_remainder_size = ctrl_onset$pre_evs_remainder_size,
    control_pre_evs_leading_edge_size = ctrl_onset$pre_evs_leading_edge_size,

    treatment_cutoff_mode = trt_onset$mode,
    treatment_cutoff_center_rank = trt_onset$cutoff_center_rank,
    treatment_cutoff_fraction = trt_onset$cutoff_center_rank / nrow(trt_rank_df),
    treatment_cutoff_range_rank_min = trt_onset$cutoff_range_rank_min,
    treatment_cutoff_range_rank_max = trt_onset$cutoff_range_rank_max,
    treatment_terminal_start_rank = trt_onset$terminal_start_rank,
    treatment_terminal_end_rank = trt_onset$terminal_end_rank,
    treatment_center_nb_support = trt_onset$center_nb_support,
    treatment_nb_support_threshold = trt_onset$nb_support_threshold,
    treatment_pre_evs_remainder_size = trt_onset$pre_evs_remainder_size,
    treatment_pre_evs_leading_edge_size = trt_onset$pre_evs_leading_edge_size
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
