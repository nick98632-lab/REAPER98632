# =============================================================================
# SEQUENCE STAGE 1: FINAL PRE-EVS CUTOFF SELECTOR
# -----------------------------------------------------------------------------
# FINAL DRAFT
#
# What this script now does
#   1. Ranks features by absolute PC1 loading
#      Left  = lowest loading
#      Right = highest loading = leading edge
#
#   2. Builds an empirical variance curve along EVS rank
#
#   3. Uses the smoothed second derivative to identify the FINAL TWO
#      pre-terminal d2 zero-crossings before the sustained terminal rise:
#         - first  of those two zeros = cutoff center
#         - second of those two zeros = terminal start
#
#   4. Uses NB-derived quantities only to define a cutoff range around the
#      derivative-defined cutoff center, not to define the center itself
#
#   5. Makes the NB2-leading-edge evidence explicit by calculating and plotting
#      left-versus-right corroboration statistics:
#         - median smoothed NB2 on left and right
#         - median smoothed NB2-NB1 contrast on left and right
#         - median smoothed log(alpha*mu) on left and right
#         - right/left ratios and right-left differences
#
#   6. Places method explanation labels on the LEFT side of the panels so the
#      selected cutoff and terminal-start lines remain visually clear
#
# Outputs
#   exports/variance_derivative_nb_range_final/
#     <comparison>_cutoff_folder/
#       <comparison>_control_rank_panel.png
#       <comparison>_treatment_rank_panel.png
#       <comparison>_cutoff_range_panel.png
#       <comparison>_feature_level_metrics.csv
#       <comparison>_control_rank_series.csv
#       <comparison>_treatment_rank_series.csv
#       <comparison>_control_zero_crossings_all.csv
#       <comparison>_treatment_zero_crossings_all.csv
#       <comparison>_control_selected_two_zero_crossings.csv
#       <comparison>_treatment_selected_two_zero_crossings.csv
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
out_root <- file.path(repo_dir, "exports", "variance_derivative_nb_range_final")
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

# NB range
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

safe_ratio <- function(num, den) {
  if (!is.finite(num) || !is.finite(den) || den == 0) return(NA_real_)
  num / den
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
    if (!is.finite(d2_vec[i - 1L]) || !is.finite(d2_vec[i])) next

    s0 <- s[i - 1L]
    s1 <- s[i]

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

find_terminal_positive_slope_run <- function(d1_sm,
                                             search_start,
                                             search_end,
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
    s_local <- good_runs$start[1]
    e_local <- good_runs$end[1]
    return(list(
      start = tail_start + s_local - 1L,
      end   = tail_start + e_local - 1L
    ))
  }

  best_local <- which.max(replace(d1_sm[tail_idx], !is.finite(d1_sm[tail_idx]), -Inf))
  if (!length(best_local) || !is.finite(best_local)) best_local <- 1L

  list(
    start = tail_idx[max(1L, best_local - terminal_run_min_length + 1L)],
    end   = tail_idx[min(length(tail_idx), best_local + terminal_run_min_length - 1L)]
  )
}

summarize_zero_crossings <- function(df, zero_idx) {
  if (!length(zero_idx)) return(data.frame())

  out <- lapply(zero_idx, function(i) {
    data.frame(
      index = i,
      rank = df$rank[i],
      d1_here = df$d1_sm[i],
      d2_here = df$d2_sm[i],
      d1_left = if (i > 1L) df$d1_sm[i - 1L] else NA_real_,
      d1_right = if (i < nrow(df)) df$d1_sm[i + 1L] else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out)
}

choose_terminal_and_cutoff_from_first_two_zeros <- function(df,
                                                            search_start,
                                                            search_end,
                                                            terminal_run_fraction = 0.20,
                                                            terminal_run_min_length = 250L,
                                                            terminal_positive_slope_quantile = 0.70) {
  run_obj <- find_terminal_positive_slope_run(
    d1_sm = df$d1_sm,
    search_start = search_start,
    search_end = search_end,
    terminal_run_fraction = terminal_run_fraction,
    terminal_run_min_length = terminal_run_min_length,
    terminal_positive_slope_quantile = terminal_positive_slope_quantile
  )

  zero_idx <- find_d2_zero_crossings(df$d2_sm)
  zero_tbl <- summarize_zero_crossings(df, zero_idx)

  if (!nrow(zero_tbl)) {
    term_idx <- run_obj$start
    cut_idx <- max(search_start, term_idx - 1L)

    return(list(
      terminal_start_index = term_idx,
      terminal_start_rank = df$rank[term_idx],
      terminal_end_index = run_obj$end,
      terminal_end_rank = df$rank[run_obj$end],
      cutoff_center_index = cut_idx,
      cutoff_center_rank = df$rank[cut_idx],
      zero_crossings = zero_tbl,
      selected_zero_crossings = data.frame(),
      mode = "fallback: no d2 zero-crossings found"
    ))
  }

  pre_terminal_tbl <- zero_tbl %>%
    filter(index <= run_obj$start)

  if (nrow(pre_terminal_tbl) < 2L) {
    term_idx <- run_obj$start
    cut_idx <- if (nrow(pre_terminal_tbl) == 1L) pre_terminal_tbl$index[1] else max(search_start, term_idx - 1L)

    sel_tbl <- pre_terminal_tbl
    if (nrow(sel_tbl)) {
      sel_tbl$role <- "cutoff_center"
    }

    return(list(
      terminal_start_index = term_idx,
      terminal_start_rank = df$rank[term_idx],
      terminal_end_index = run_obj$end,
      terminal_end_rank = df$rank[run_obj$end],
      cutoff_center_index = cut_idx,
      cutoff_center_rank = df$rank[cut_idx],
      zero_crossings = zero_tbl,
      selected_zero_crossings = sel_tbl,
      mode = "fallback: fewer than two pre-terminal d2 zero-crossings"
    ))
  }

  pair_tbl <- tail(pre_terminal_tbl, 2)
  cut_idx <- pair_tbl$index[1]
  term_idx <- pair_tbl$index[2]
  pair_tbl$role <- c("cutoff_center", "terminal_start")

  list(
    terminal_start_index = term_idx,
    terminal_start_rank = df$rank[term_idx],
    terminal_end_index = run_obj$end,
    terminal_end_rank = df$rank[run_obj$end],
    cutoff_center_index = cut_idx,
    cutoff_center_rank = df$rank[cut_idx],
    zero_crossings = zero_tbl,
    selected_zero_crossings = pair_tbl,
    mode = "terminal_start = second of final two d2 zero-crossings before sustained rise; cutoff_center = first of those two zeros"
  )
}

compute_nb_support_score <- function(df, nb_smooth_window = 151L) {
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

compute_left_right_corrob <- function(df, cutoff_center_index) {
  n <- nrow(df)
  left_idx <- seq_len(max(1L, cutoff_center_index - 1L))
  right_idx <- seq.int(cutoff_center_index, n)

  med_left_nb2 <- stats::median(df$log_nb2[left_idx], na.rm = TRUE)
  med_right_nb2 <- stats::median(df$log_nb2[right_idx], na.rm = TRUE)

  med_left_gap <- stats::median(df$nb_gap_sm[left_idx], na.rm = TRUE)
  med_right_gap <- stats::median(df$nb_gap_sm[right_idx], na.rm = TRUE)

  med_left_amu <- stats::median(df$amu_sm[left_idx], na.rm = TRUE)
  med_right_amu <- stats::median(df$amu_sm[right_idx], na.rm = TRUE)

  med_left_support <- stats::median(df$nb_support[left_idx], na.rm = TRUE)
  med_right_support <- stats::median(df$nb_support[right_idx], na.rm = TRUE)

  data.frame(
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

    stringsAsFactors = FALSE
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
                                              nb_range_max_span_fraction = 0.18,
                                              nb_smooth_window = 151L) {
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

  df$var_fit <- NA_real_
  df$d1_sm <- NA_real_
  df$d2_sm <- NA_real_

  df$var_fit[ok] <- pred0$y
  df$d1_sm[ok] <- pred1$y
  df$d2_sm[ok] <- pred2$y

  search_start <- max(left_edge_buffer + 1L, floor(search_fraction_min * n))
  search_end <- min(n - right_edge_buffer, floor(search_fraction_max * n))
  if (search_end <= search_start + 10L) {
    search_start <- max(1L, left_edge_buffer + 1L)
    search_end <- min(n - right_edge_buffer, n - 1L)
  }

  geom_obj <- choose_terminal_and_cutoff_from_first_two_zeros(
    df = df,
    search_start = search_start,
    search_end = search_end,
    terminal_run_fraction = terminal_run_fraction,
    terminal_run_min_length = terminal_run_min_length,
    terminal_positive_slope_quantile = terminal_positive_slope_quantile
  )

  nb_obj <- compute_nb_support_score(df, nb_smooth_window = nb_smooth_window)
  df$nb_gap_sm <- nb_obj$nb_gap_sm
  df$amu_sm <- nb_obj$amu_sm
  df$nb_support <- nb_obj$nb_support

  nb_rng <- expand_nb_range_around_center(
    df = df,
    center_idx = geom_obj$cutoff_center_index,
    terminal_start_idx = geom_obj$terminal_start_index,
    nb_support = df$nb_support,
    nb_range_drop_fraction = nb_range_drop_fraction,
    nb_range_max_span_fraction = nb_range_max_span_fraction
  )

  corrob <- compute_left_right_corrob(df, geom_obj$cutoff_center_index)

  total_features <- n
  pre_evs_remainder_size <- geom_obj$cutoff_center_index - 1L
  pre_evs_leading_edge_size <- total_features - geom_obj$cutoff_center_index + 1L

  list(
    curve_df = df,
    zero_crossings = geom_obj$zero_crossings,
    selected_zero_crossings = geom_obj$selected_zero_crossings,
    mode = geom_obj$mode,

    terminal_start_index = geom_obj$terminal_start_index,
    terminal_start_rank = geom_obj$terminal_start_rank,
    terminal_end_index = geom_obj$terminal_end_index,
    terminal_end_rank = geom_obj$terminal_end_rank,

    cutoff_center_index = geom_obj$cutoff_center_index,
    cutoff_center_rank = geom_obj$cutoff_center_rank,

    cutoff_range_index_min = nb_rng$range_left_index,
    cutoff_range_rank_min = nb_rng$range_left_rank,
    cutoff_range_index_max = nb_rng$range_right_index,
    cutoff_range_rank_max = nb_rng$range_right_rank,

    center_nb_support = nb_rng$center_nb_support,
    nb_support_threshold = nb_rng$nb_support_threshold,

    total_features = total_features,
    pre_evs_remainder_size = pre_evs_remainder_size,
    pre_evs_leading_edge_size = pre_evs_leading_edge_size,

    corrob = corrob
  )
}

# =============================================================================
# FIGURES
# =============================================================================

build_dataset_panel <- function(rank_df, onset_info, title_prefix, out_file) {
  df <- onset_info$curve_df
  center_x <- onset_info$cutoff_center_rank
  range_min_x <- onset_info$cutoff_range_rank_min
  range_max_x <- onset_info$cutoff_range_rank_max
  terminal_x <- onset_info$terminal_start_rank
  terminal_end_x <- onset_info$terminal_end_rank
  corr <- onset_info$corrob

  x_rng <- range(df$rank, na.rm = TRUE)
  x_left <- x_rng[1] + 0.06 * diff(x_rng)
  x_right <- x_rng[1] + 0.72 * diff(x_rng)

  p1 <- ggplot(df, aes(rank, abs_loading)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    annotate("rect", xmin = range_min_x, xmax = range_max_x, ymin = -Inf, ymax = Inf, alpha = 0.08) +
    annotate("rect", xmin = terminal_x, xmax = terminal_end_x, ymin = -Inf, ymax = Inf, alpha = 0.08) +
    geom_vline(xintercept = center_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    annotate(
      "label",
      x = x_right,
      y = max(df$abs_loading, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label = paste0(
        "Absolute loading panel",
        "\nCutoff center rank = ", center_x,
        "\nCutoff range = [", range_min_x, ", ", range_max_x, "]",
        "\nTerminal start rank = ", terminal_x,
        "\nPre-EVS remainder = ", onset_info$pre_evs_remainder_size,
        "\nPre-EVS leading edge = ", onset_info$pre_evs_leading_edge_size
      )
    ) +
    labs(
      title = paste0(title_prefix, ": absolute PC1 loading series"),
      subtitle = "Leading edge is on the RIGHT",
      x = "EVS rank",
      y = "|PC1 loading|"
    ) +
    theme_bw(base_size = 10)

  p2 <- ggplot(df, aes(rank, var_fit)) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    annotate("rect", xmin = range_min_x, xmax = range_max_x, ymin = -Inf, ymax = Inf, alpha = 0.08) +
    annotate("rect", xmin = terminal_x, xmax = terminal_end_x, ymin = -Inf, ymax = Inf, alpha = 0.08) +
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
    annotate(
      "label",
      x = x_left,
      y = max(df$var_fit, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.7,
      label = paste0(
        "Variance curve method",
        "\nUse the FINAL TWO pre-terminal d2 zero-crossings",
        "\nZero 1 = cutoff center",
        "\nZero 2 = terminal start",
        "\nThese mark the late pre-terminal transition"
      )
    ) +
    labs(
      title = paste0(title_prefix, ": smoothed empirical variance curve"),
      subtitle = "Selected geometric points are marked on the curve",
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    theme_bw(base_size = 10)

  p3 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = d1_sm, color = "Smoothed slope d1"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(aes(y = d2_sm, color = "Smoothed curvature d2"), linewidth = 1.0, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    annotate("rect", xmin = range_min_x, xmax = range_max_x, ymin = -Inf, ymax = Inf, alpha = 0.08) +
    annotate("rect", xmin = terminal_x, xmax = terminal_end_x, ymin = -Inf, ymax = Inf, alpha = 0.08) +
    geom_vline(xintercept = center_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    annotate(
      "label",
      x = x_left,
      y = max(c(df$d1_sm, df$d2_sm), na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.7,
      label = paste0(
        "Derivative method",
        "\nCutoff center = first of final two d2 zeros",
        "\nTerminal start = second of final two d2 zeros",
        "\nBoth occur before the sustained positive-slope rise"
      )
    ) +
    labs(
      title = paste0(title_prefix, ": derivative support"),
      subtitle = "Center and terminal start come only from the selected two zero-crossings",
      x = "EVS rank",
      y = "Derivative value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p4 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = nb_support, color = "Combined NB support"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(aes(y = log_nb1, color = "NB1 = mu"), linewidth = 0.5, alpha = 0.35, na.rm = TRUE) +
    geom_line(aes(y = log_nb2, color = "NB2 = variance - mu"), linewidth = 0.5, alpha = 0.35, na.rm = TRUE) +
    geom_line(aes(y = amu_sm, color = "Smoothed log(alpha*mu)"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(aes(y = nb_gap_sm, color = "Smoothed log(NB2+1) - log(NB1+1)"), linewidth = 1.0, na.rm = TRUE) +
    annotate("rect", xmin = range_min_x, xmax = range_max_x, ymin = -Inf, ymax = Inf, alpha = 0.08) +
    geom_vline(xintercept = center_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    geom_hline(yintercept = onset_info$nb_support_threshold, linetype = 3, linewidth = 0.5) +
    annotate(
      "label",
      x = x_left,
      y = max(c(df$nb_support, df$log_nb1, df$log_nb2, df$amu_sm, df$nb_gap_sm), na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.55,
      label = paste0(
        "NB corroboration",
        "\nNB quantities DO NOT define the center",
        "\nThey define the RANGE around the derivative-selected center",
        "\nLeft med log(NB2)  = ", round(corr$left_median_log_nb2, 3),
        "\nRight med log(NB2) = ", round(corr$right_median_log_nb2, 3),
        "\nRight-left log(NB2) diff = ", round(corr$right_left_log_nb2_diff, 3),
        "\nLeft med NB2-NB1 contrast  = ", round(corr$left_median_nb_gap, 3),
        "\nRight med NB2-NB1 contrast = ", round(corr$right_median_nb_gap, 3),
        "\nRight-left contrast diff = ", round(corr$right_left_nb_gap_diff, 3),
        "\nLeft med log(alpha*mu)  = ", round(corr$left_median_log_alpha_mu, 3),
        "\nRight med log(alpha*mu) = ", round(corr$right_median_log_alpha_mu, 3),
        "\nRight-left log(alpha*mu) diff = ", round(corr$right_left_log_alpha_mu_diff, 3),
        "\nCenter NB support = ", round(onset_info$center_nb_support, 3),
        "\nNB threshold = ", round(onset_info$nb_support_threshold, 3)
      )
    ) +
    labs(
      title = paste0(title_prefix, ": NB1 / NB2 / alpha*mu support"),
      subtitle = "Right-of-cutoff elevation in NB2 and alpha*mu supports the leading-edge interpretation",
      x = "EVS rank",
      y = "Support value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p5 <- ggplot(df, aes(rank, nb_support)) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    annotate("rect", xmin = range_min_x, xmax = range_max_x, ymin = -Inf, ymax = Inf, alpha = 0.08) +
    geom_vline(xintercept = center_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    geom_hline(yintercept = onset_info$nb_support_threshold, linetype = 3, linewidth = 0.5) +
    annotate(
      "label",
      x = x_right,
      y = max(df$nb_support, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.7,
      label = paste0(
        "Final NB-supported band",
        "\nCenter = ", center_x,
        "\nBand = [", range_min_x, ", ", range_max_x, "]",
        "\nTerminal start = ", terminal_x
      )
    ) +
    labs(
      title = paste0(title_prefix, ": NB-supported cutoff band"),
      subtitle = "Band expands from the center while combined NB support stays elevated",
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

  x_all <- c(ctrl_df$rank, trt_df$rank)
  y_all <- c(ctrl_df$var_fit, trt_df$var_fit)
  x_left <- min(x_all, na.rm = TRUE) + 0.06 * diff(range(x_all, na.rm = TRUE))

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
    annotate(
      "label",
      x = x_left,
      y = max(y_all, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label = paste0(
        "Treatment / control comparison",
        "\nControl center = ", ctrl_onset$cutoff_center_rank,
        "\nTreatment center = ", trt_onset$cutoff_center_rank,
        "\nControl range = [", ctrl_onset$cutoff_range_rank_min, ", ", ctrl_onset$cutoff_range_rank_max, "]",
        "\nTreatment range = [", trt_onset$cutoff_range_rank_min, ", ", trt_onset$cutoff_range_rank_max, "]"
      )
    ) +
    labs(
      title = paste0(title_prefix, ": treatment/control derivative-defined centers and NB-supported ranges"),
      subtitle = "Both arms use the same integrated method",
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
    nb_range_max_span_fraction = nb_range_max_span_fraction,
    nb_smooth_window = nb_smooth_window
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
    nb_range_max_span_fraction = nb_range_max_span_fraction,
    nb_smooth_window = nb_smooth_window
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

    control_mode = ctrl_onset$mode,
    control_cutoff_center_rank = ctrl_onset$cutoff_center_rank,
    control_terminal_start_rank = ctrl_onset$terminal_start_rank,
    control_cutoff_range_rank_min = ctrl_onset$cutoff_range_rank_min,
    control_cutoff_range_rank_max = ctrl_onset$cutoff_range_rank_max,
    control_center_nb_support = ctrl_onset$center_nb_support,
    control_nb_support_threshold = ctrl_onset$nb_support_threshold,
    control_pre_evs_remainder_size = ctrl_onset$pre_evs_remainder_size,
    control_pre_evs_leading_edge_size = ctrl_onset$pre_evs_leading_edge_size,
    control_left_median_log_nb2 = ctrl_corr$left_median_log_nb2,
    control_right_median_log_nb2 = ctrl_corr$right_median_log_nb2,
    control_right_left_log_nb2_diff = ctrl_corr$right_left_log_nb2_diff,
    control_left_median_nb_gap = ctrl_corr$left_median_nb_gap,
    control_right_median_nb_gap = ctrl_corr$right_median_nb_gap,
    control_right_left_nb_gap_diff = ctrl_corr$right_left_nb_gap_diff,
    control_left_median_log_alpha_mu = ctrl_corr$left_median_log_alpha_mu,
    control_right_median_log_alpha_mu = ctrl_corr$right_median_log_alpha_mu,
    control_right_left_log_alpha_mu_diff = ctrl_corr$right_left_log_alpha_mu_diff,

    treatment_mode = trt_onset$mode,
    treatment_cutoff_center_rank = trt_onset$cutoff_center_rank,
    treatment_terminal_start_rank = trt_onset$terminal_start_rank,
    treatment_cutoff_range_rank_min = trt_onset$cutoff_range_rank_min,
    treatment_cutoff_range_rank_max = trt_onset$cutoff_range_rank_max,
    treatment_center_nb_support = trt_onset$center_nb_support,
    treatment_nb_support_threshold = trt_onset$nb_support_threshold,
    treatment_pre_evs_remainder_size = trt_onset$pre_evs_remainder_size,
    treatment_pre_evs_leading_edge_size = trt_onset$pre_evs_leading_edge_size,
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
