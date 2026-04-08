# =============================================================================
# SEQUENCE STAGE 1: MANUAL-5000 CUTOFF WITH LOCAL RIGHT-ANCHORED CORROBORATION
# -----------------------------------------------------------------------------
# FINAL DRAFT
#
# What this script does
#   1. Ranks features by absolute PC1 loading
#      Left  = lowest loading
#      Right = highest loading = leading edge
#
#   2. Defines the leading edge a priori as the top 5000 ranked features
#
#   3. Builds an empirical variance curve along EVS rank
#
#   4. Uses the smoothed second derivative only to corroborate the manual
#      5000-feature cutoff by selecting a compact local zero-crossing band in
#      the neighborhood of the manual cutoff
#
#   5. Summarizes NB-derived support numerically on both sides of the manual
#      cutoff and retains a single Boolean criterion:
#         - post_cutoff_more_nb2
#
#   6. Writes per-comparison folders and an overall summary table using the
#      same server/repo assumptions as the existing script
#
# Outputs
#   exports/variance_derivative_nb_range_manual5000_option2/
#     <comparison>_cutoff_folder/
#       <comparison>_control_rank_panel.png
#       <comparison>_treatment_rank_panel.png
#       <comparison>_feature_level_metrics.csv
#       <comparison>_control_rank_series.csv
#       <comparison>_treatment_rank_series.csv
#       <comparison>_control_zero_crossings_all.csv
#       <comparison>_treatment_zero_crossings_all.csv
#       <comparison>_control_selected_zero_pair.csv
#       <comparison>_treatment_selected_zero_pair.csv
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
out_root <- file.path(repo_dir, "exports", "variance_derivative_nb_range_manual5000_option2")
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

# Manual cutoff and local corroboration settings
manual_leading_edge_n <- 5000L
spline_spar <- 0.60
left_edge_buffer <- 50L
right_edge_buffer <- 20L
search_fraction_min <- 0.20
search_fraction_max <- 0.985
manual_cutoff_search_left_fraction <- 0.12
manual_cutoff_search_right_fraction <- 0.12
manual_cutoff_pair_max_width_fraction <- 0.06
local_nb_support_quantile <- 0.30
nb_smooth_window <- 151L

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

safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_var <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::var(x)
}

safe_quantile <- function(x, probs) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8))
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

smooth_rank_series <- function(df, spar = spline_spar, nb_window = nb_smooth_window) {
  df <- df %>% arrange(rank)
  fit_ok <- is.finite(df$rank) & is.finite(df$log_variance)

  if (sum(fit_ok) < 10L) {
    df$var_smooth <- df$log_variance
    df$d1_sm <- NA_real_
    df$d2_sm <- NA_real_
  } else {
    sp <- smooth.spline(x = df$rank[fit_ok], y = df$log_variance[fit_ok], spar = spar)
    pred0 <- predict(sp, x = df$rank, deriv = 0)
    pred1 <- predict(sp, x = df$rank, deriv = 1)
    pred2 <- predict(sp, x = df$rank, deriv = 2)
    df$var_smooth <- pred0$y
    df$d1_sm <- pred1$y
    df$d2_sm <- pred2$y
  }

  df$nb2_smooth <- roll_median(df$log_nb2, nb_window)
  df$nb_gap_smooth <- roll_median(df$nb_gap, nb_window)
  df$alpha_mu_smooth <- roll_median(df$log_alpha_mu, nb_window)

  df$combined_nb_support <- rowMeans(cbind(
    scale01(df$nb2_smooth),
    scale01(df$nb_gap_smooth),
    scale01(df$alpha_mu_smooth)
  ), na.rm = TRUE)

  df
}

find_zero_crossings <- function(d2_vec, rank_vec) {
  d2_vec <- as.numeric(d2_vec)
  rank_vec <- as.numeric(rank_vec)
  n <- length(d2_vec)
  if (n < 2L) return(data.frame())

  out <- list()
  j <- 1L
  for (i in 2:n) {
    y0 <- d2_vec[i - 1L]
    y1 <- d2_vec[i]
    x0 <- rank_vec[i - 1L]
    x1 <- rank_vec[i]
    if (!is.finite(y0) || !is.finite(y1) || !is.finite(x0) || !is.finite(x1)) next

    sign_change <- (sign(y0) != sign(y1)) || y0 == 0 || y1 == 0
    if (!sign_change) next

    rank_cross <- if ((y1 - y0) == 0) x1 else x0 - y0 * (x1 - x0) / (y1 - y0)
    crossing_type <- if (y0 < 0 && y1 > 0) {
      "neg_to_pos"
    } else if (y0 > 0 && y1 < 0) {
      "pos_to_neg"
    } else {
      "touch_or_flat"
    }

    out[[j]] <- data.frame(
      index_left = i - 1L,
      index_right = i,
      rank = rank_cross,
      d2_left = y0,
      d2_right = y1,
      crossing_type = crossing_type,
      stringsAsFactors = FALSE
    )
    j <- j + 1L
  }

  if (!length(out)) return(data.frame())
  bind_rows(out)
}

select_zero_pair_supporting_manual_cutoff <- function(df,
                                                      manual_cutoff_rank = manual_leading_edge_n,
                                                      search_left_fraction = manual_cutoff_search_left_fraction,
                                                      search_right_fraction = manual_cutoff_search_right_fraction,
                                                      max_width_fraction = manual_cutoff_pair_max_width_fraction,
                                                      support_quantile = local_nb_support_quantile) {
  n <- nrow(df)
  manual_cutoff_rank <- min(manual_cutoff_rank, n)

  zero_tbl <- find_zero_crossings(df$d2_sm, df$rank)
  if (!nrow(zero_tbl)) {
    return(list(
      manual_cutoff_rank = manual_cutoff_rank,
      cutoff_center_rank = manual_cutoff_rank,
      terminal_start_rank = manual_cutoff_rank,
      cutoff_range_rank_min = manual_cutoff_rank,
      cutoff_range_rank_max = manual_cutoff_rank,
      selected_zero_pair = data.frame(),
      zero_crossings = zero_tbl,
      mode = "fallback: no zero-crossings found"
    ))
  }

  zero_tbl$local_support <- vapply(zero_tbl$rank, function(rk) {
    idx <- which.min(abs(df$rank - rk))
    df$combined_nb_support[idx]
  }, numeric(1))

  support_threshold <- safe_quantile(df$combined_nb_support, support_quantile)
  if (!is.finite(support_threshold)) support_threshold <- -Inf

  left_bound <- max(1, floor(manual_cutoff_rank - n * search_left_fraction))
  right_bound <- min(n, ceiling(manual_cutoff_rank + n * search_right_fraction))
  max_width <- max(25L, floor(n * max_width_fraction))

  local_tbl <- zero_tbl %>%
    filter(rank >= left_bound, rank <= right_bound)

  if (!nrow(local_tbl)) local_tbl <- zero_tbl

  local_tbl <- local_tbl %>%
    mutate(
      is_supported = is.finite(local_support) & local_support >= support_threshold,
      is_neg_to_pos = crossing_type == "neg_to_pos"
    )

  left_tbl <- local_tbl %>%
    filter(rank <= manual_cutoff_rank, is_neg_to_pos) %>%
    arrange(desc(rank))
  right_tbl <- local_tbl %>%
    filter(rank >= manual_cutoff_rank, is_neg_to_pos) %>%
    arrange(rank)

  if (!nrow(left_tbl)) {
    left_tbl <- local_tbl %>% filter(rank <= manual_cutoff_rank) %>% arrange(desc(rank))
  }
  if (!nrow(right_tbl)) {
    right_tbl <- local_tbl %>% filter(rank >= manual_cutoff_rank) %>% arrange(rank)
  }
  if (!nrow(left_tbl)) {
    left_tbl <- zero_tbl %>% filter(rank < manual_cutoff_rank) %>% arrange(desc(rank))
  }
  if (!nrow(right_tbl)) {
    right_tbl <- zero_tbl %>% filter(rank > manual_cutoff_rank) %>% arrange(rank)
  }

  if (!nrow(left_tbl) || !nrow(right_tbl)) {
    return(list(
      manual_cutoff_rank = manual_cutoff_rank,
      cutoff_center_rank = manual_cutoff_rank,
      terminal_start_rank = manual_cutoff_rank,
      cutoff_range_rank_min = manual_cutoff_rank,
      cutoff_range_rank_max = manual_cutoff_rank,
      selected_zero_pair = data.frame(),
      zero_crossings = zero_tbl,
      mode = "fallback: no flanking pair found"
    ))
  }

  pair_tbl <- expand.grid(left_i = seq_len(nrow(left_tbl)), right_i = seq_len(nrow(right_tbl)))
  pair_tbl$left_rank <- left_tbl$rank[pair_tbl$left_i]
  pair_tbl$right_rank <- right_tbl$rank[pair_tbl$right_i]
  pair_tbl$width <- pair_tbl$right_rank - pair_tbl$left_rank
  pair_tbl$midpoint <- (pair_tbl$left_rank + pair_tbl$right_rank) / 2
  pair_tbl$midpoint_distance_to_manual <- abs(pair_tbl$midpoint - manual_cutoff_rank)
  pair_tbl$support_score <- left_tbl$local_support[pair_tbl$left_i] + right_tbl$local_support[pair_tbl$right_i]
  pair_tbl$all_supported <- left_tbl$is_supported[pair_tbl$left_i] & right_tbl$is_supported[pair_tbl$right_i]

  pair_tbl <- pair_tbl[pair_tbl$width > 0, , drop = FALSE]
  narrow_tbl <- pair_tbl[pair_tbl$width <= max_width, , drop = FALSE]
  if (nrow(narrow_tbl)) pair_tbl <- narrow_tbl

  pair_tbl <- pair_tbl[order(
    !pair_tbl$all_supported,
    pair_tbl$midpoint_distance_to_manual,
    pair_tbl$width,
    -pair_tbl$support_score
  ), , drop = FALSE]

  best <- pair_tbl[1, , drop = FALSE]

  selected_tbl <- data.frame(
    manual_cutoff_rank = manual_cutoff_rank,
    cutoff_center_rank = manual_cutoff_rank,
    terminal_start_rank = best$right_rank,
    cutoff_range_rank_min = best$left_rank,
    cutoff_range_rank_max = best$right_rank,
    midpoint_rank = best$midpoint,
    midpoint_distance_to_manual = best$midpoint_distance_to_manual,
    pair_width = best$width,
    support_score = best$support_score,
    all_supported = best$all_supported,
    stringsAsFactors = FALSE
  )

  list(
    manual_cutoff_rank = manual_cutoff_rank,
    cutoff_center_rank = manual_cutoff_rank,
    terminal_start_rank = best$right_rank,
    cutoff_range_rank_min = best$left_rank,
    cutoff_range_rank_max = best$right_rank,
    selected_zero_pair = selected_tbl,
    zero_crossings = zero_tbl,
    mode = "manual_5000_local_pair"
  )
}

summarize_manual_cutoff_support <- function(df, selected_obj, comparison_name, arm_label) {
  cutoff_rank <- selected_obj$manual_cutoff_rank
  pre_idx <- which(df$rank < cutoff_rank)
  post_idx <- which(df$rank >= cutoff_rank)

  if (!length(pre_idx)) pre_idx <- 1L
  if (!length(post_idx)) post_idx <- nrow(df)

  pre_nb2 <- safe_median(df$nb2_smooth[pre_idx])
  post_nb2 <- safe_median(df$nb2_smooth[post_idx])
  pre_gap <- safe_median(df$nb_gap_smooth[pre_idx])
  post_gap <- safe_median(df$nb_gap_smooth[post_idx])
  pre_alpha_mu <- safe_median(df$alpha_mu_smooth[pre_idx])
  post_alpha_mu <- safe_median(df$alpha_mu_smooth[post_idx])

  data.frame(
    comparison = comparison_name,
    group = arm_label,
    manual_leading_edge_n = cutoff_rank,
    cutoff_center_rank = selected_obj$cutoff_center_rank,
    terminal_start_rank = selected_obj$terminal_start_rank,
    cutoff_range_rank_min = selected_obj$cutoff_range_rank_min,
    cutoff_range_rank_max = selected_obj$cutoff_range_rank_max,
    pre_cutoff_size = length(pre_idx),
    post_cutoff_size = length(post_idx),
    pre_cutoff_median_log_nb2 = pre_nb2,
    post_cutoff_median_log_nb2 = post_nb2,
    post_minus_pre_log_nb2 = post_nb2 - pre_nb2,
    pre_cutoff_median_nb2_nb1_contrast = pre_gap,
    post_cutoff_median_nb2_nb1_contrast = post_gap,
    post_minus_pre_nb2_nb1_contrast = post_gap - pre_gap,
    pre_cutoff_median_log_alpha_mu = pre_alpha_mu,
    post_cutoff_median_log_alpha_mu = post_alpha_mu,
    post_minus_pre_log_alpha_mu = post_alpha_mu - pre_alpha_mu,
    pre_cutoff_mean_combined_nb_support = safe_mean(df$combined_nb_support[pre_idx]),
    post_cutoff_mean_combined_nb_support = safe_mean(df$combined_nb_support[post_idx]),
    post_minus_pre_mean_combined_nb_support =
      safe_mean(df$combined_nb_support[post_idx]) - safe_mean(df$combined_nb_support[pre_idx]),
    post_cutoff_more_nb2 = isTRUE(post_nb2 > pre_nb2),
    selection_mode = selected_obj$mode,
    stringsAsFactors = FALSE
  )
}

make_annotation_text <- function(summary_row) {
  paste(
    sprintf("Manual leading edge = %s", summary_row$manual_leading_edge_n),
    sprintf("Cutoff center rank = %s", summary_row$cutoff_center_rank),
    sprintf("Corroboration band = [%s, %s]",
            summary_row$cutoff_range_rank_min,
            summary_row$cutoff_range_rank_max),
    sprintf("Terminal start rank = %s", summary_row$terminal_start_rank),
    sprintf("Pre med log(NB2) = %.3f", summary_row$pre_cutoff_median_log_nb2),
    sprintf("Post med log(NB2) = %.3f", summary_row$post_cutoff_median_log_nb2),
    sprintf("Pre med NB2-NB1 = %.3f", summary_row$pre_cutoff_median_nb2_nb1_contrast),
    sprintf("Post med NB2-NB1 = %.3f", summary_row$post_cutoff_median_nb2_nb1_contrast),
    sprintf("Pre med log(alpha*mu) = %.3f", summary_row$pre_cutoff_median_log_alpha_mu),
    sprintf("Post med log(alpha*mu) = %.3f", summary_row$post_cutoff_median_log_alpha_mu),
    sprintf("Post-cutoff more NB2 = %s", ifelse(summary_row$post_cutoff_more_nb2, "TRUE", "FALSE")),
    sep = "\n"
  )
}

plot_rank_panel <- function(df, summary_row, zero_tbl, comparison_name, arm_label, out_file) {
  note_text <- make_annotation_text(summary_row)

  p1 <- ggplot(df, aes(x = rank, y = abs_loading)) +
    annotate("rect",
             xmin = summary_row$cutoff_range_rank_min,
             xmax = summary_row$cutoff_range_rank_max,
             ymin = -Inf, ymax = Inf, alpha = 0.12) +
    geom_vline(xintercept = summary_row$manual_leading_edge_n, linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = summary_row$terminal_start_rank, linetype = "dotted", linewidth = 0.5) +
    geom_line(linewidth = 0.6) +
    labs(
      title = paste(comparison_name, arm_label, "absolute PC1 loading series"),
      subtitle = "Leading edge is defined a priori as the top 5000 ranked features",
      x = "EVS rank",
      y = "|PC1 loading|"
    ) +
    theme_bw(base_size = 10)

  p2 <- ggplot(df, aes(x = rank, y = var_smooth)) +
    annotate("rect",
             xmin = summary_row$cutoff_range_rank_min,
             xmax = summary_row$cutoff_range_rank_max,
             ymin = -Inf, ymax = Inf, alpha = 0.12) +
    geom_vline(xintercept = summary_row$manual_leading_edge_n, linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = summary_row$terminal_start_rank, linetype = "dotted", linewidth = 0.5) +
    geom_line(linewidth = 0.6) +
    annotate("text", x = min(df$rank, na.rm = TRUE) + 0.03 * diff(range(df$rank, na.rm = TRUE)),
             y = max(df$var_smooth, na.rm = TRUE), label = note_text,
             hjust = 0, vjust = 1, size = 3) +
    labs(
      title = paste(comparison_name, arm_label, "smoothed empirical variance curve"),
      subtitle = "Local second-derivative band corroborates the manual 5000-feature boundary",
      x = "EVS rank",
      y = "Smoothed log(1 + variance)"
    ) +
    theme_bw(base_size = 10)

  p3 <- ggplot(df, aes(x = rank, y = d2_sm)) +
    annotate("rect",
             xmin = summary_row$cutoff_range_rank_min,
             xmax = summary_row$cutoff_range_rank_max,
             ymin = -Inf, ymax = Inf, alpha = 0.12) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_vline(xintercept = summary_row$manual_leading_edge_n, linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = summary_row$terminal_start_rank, linetype = "dotted", linewidth = 0.5) +
    geom_line(linewidth = 0.6) +
    geom_point(data = zero_tbl, aes(x = rank, y = 0), inherit.aes = FALSE, size = 1.2) +
    labs(
      title = paste(comparison_name, arm_label, "second derivative support"),
      subtitle = "Eligible local zero-crossings around the manual cutoff",
      x = "EVS rank",
      y = "Smoothed second derivative"
    ) +
    theme_bw(base_size = 10)

  p4 <- ggplot(df, aes(x = rank)) +
    geom_vline(xintercept = summary_row$manual_leading_edge_n, linetype = "dashed", linewidth = 0.5) +
    geom_line(aes(y = nb2_smooth, colour = "log(NB2)"), linewidth = 0.6) +
    geom_line(aes(y = nb_gap_smooth, colour = "NB2-NB1"), linewidth = 0.6) +
    geom_line(aes(y = alpha_mu_smooth, colour = "log(alpha*mu)"), linewidth = 0.6) +
    labs(
      title = paste(comparison_name, arm_label, "NB corroboration"),
      subtitle = "Post-cutoff support is evaluated from rank 5000 to the end of the series",
      x = "EVS rank",
      y = "Support value",
      colour = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  grob <- gridExtra::arrangeGrob(p1, p2, p3, p4, ncol = 1)
  ggsave(out_file, grob, width = 11, height = 15, dpi = 300)
}

# =============================================================================
# INPUT DATA
# =============================================================================

count_file <- resolve_counts_file(count_file)
message("Using count file: ", count_file)

count_obj <- read_count_matrix(count_file, meta_all$id)
count_matrix <- count_obj$count_matrix
annot_df <- count_obj$annot_df
feature_ids <- rownames(count_matrix)

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

overall_summary_list <- list()

for (ii in seq_len(nrow(comparison_table))) {
  comp_row <- comparison_table[ii, , drop = FALSE]
  comparison_name <- comp_row$comparison_name
  message("Processing ", comparison_name)

  sub_obj <- subset_comparison(count_matrix, comp_row, meta_all)
  count_sub <- sub_obj$count_matrix

  evs_tbl <- build_evs_table(
    count_mat = count_sub,
    ctrl_cols = sub_obj$ctrl_ids,
    trt_cols = sub_obj$trt_ids,
    feature_ids = rownames(count_sub)
  )

  feature_tbl <- compute_feature_metrics_empirical(
    count_mat = count_sub,
    ctrl_cols = sub_obj$ctrl_ids,
    trt_cols = sub_obj$trt_ids,
    feature_ids = rownames(count_sub)
  )

  full_tbl <- evs_tbl %>%
    left_join(feature_tbl, by = "feature_id") %>%
    left_join(annot_df, by = "feature_id")

  comp_dir <- file.path(out_root, paste0(comparison_name, "_cutoff_folder"))
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

  utils::write.csv(
    full_tbl,
    file = file.path(comp_dir, paste0(comparison_name, "_feature_level_metrics.csv")),
    row.names = FALSE
  )

  comp_summary_rows <- list()

  for (arm_label in c("control", "treatment")) {
    rank_df <- build_rank_series(full_tbl, arm = arm_label)
    rank_df <- smooth_rank_series(rank_df)

    selected_obj <- select_zero_pair_supporting_manual_cutoff(
      df = rank_df,
      manual_cutoff_rank = manual_leading_edge_n
    )

    summary_row <- summarize_manual_cutoff_support(
      df = rank_df,
      selected_obj = selected_obj,
      comparison_name = comparison_name,
      arm_label = arm_label
    )

    comp_summary_rows[[arm_label]] <- summary_row
    overall_summary_list[[paste(comparison_name, arm_label, sep = "__")]] <- summary_row

    utils::write.csv(
      rank_df,
      file = file.path(comp_dir, paste0(comparison_name, "_", arm_label, "_rank_series.csv")),
      row.names = FALSE
    )

    utils::write.csv(
      selected_obj$zero_crossings,
      file = file.path(comp_dir, paste0(comparison_name, "_", arm_label, "_zero_crossings_all.csv")),
      row.names = FALSE
    )

    utils::write.csv(
      selected_obj$selected_zero_pair,
      file = file.path(comp_dir, paste0(comparison_name, "_", arm_label, "_selected_zero_pair.csv")),
      row.names = FALSE
    )

    plot_rank_panel(
      df = rank_df,
      summary_row = summary_row,
      zero_tbl = selected_obj$zero_crossings,
      comparison_name = comparison_name,
      arm_label = arm_label,
      out_file = file.path(comp_dir, paste0(comparison_name, "_", arm_label, "_rank_panel.png"))
    )
  }

  comp_summary_df <- bind_rows(comp_summary_rows)
  utils::write.csv(
    comp_summary_df,
    file = file.path(comp_dir, paste0(comparison_name, "_cutoff_summary.csv")),
    row.names = FALSE
  )
}

overall_cutoff_summary <- bind_rows(overall_summary_list)
utils::write.csv(
  overall_cutoff_summary,
  file = file.path(out_root, "overall_cutoff_summary.csv"),
  row.names = FALSE
)

# =============================================================================
# MANUSCRIPT METHODS TEXT MIRROR
# =============================================================================
# Features were ranked separately within each comparison and group using the
# absolute magnitude of the first principal-component loading derived from
# log2(x + 1)-transformed counts. The leading edge was defined a priori as the
# top 5000 ranked features. Smoothed variance geometry was then used only to
# corroborate this manually specified boundary rather than to determine it.
# Specifically, the smoothed second derivative of the ranked variance curve was
# evaluated within a local neighborhood around the manual cutoff, and a compact
# corroboration band was defined by an eligible pair of nearby zero-crossings
# selected to flank or closely bracket that boundary while preserving local NB
# support. Negative-binomial corroboration was summarized numerically on both
# sides of the manual cutoff using median log(NB2), median NB2 minus NB1
# contrast, median log(alpha·mu), and mean combined NB support. The principal
# binary corroboration criterion was whether the post-cutoff region, defined
# from rank 5000 to the end of the ranked series, showed greater median NB2
# support than the pre-cutoff region.
# =============================================================================
