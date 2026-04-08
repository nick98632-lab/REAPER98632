# =============================================================================
# SEQUENCE STAGE 1: MANUAL TOP-5000 PRE-EVS CUTOFF CORROBORATION
# -----------------------------------------------------------------------------
# FINAL MANUSCRIPT VERSION
#
# What this script does
#   1. Ranks features by absolute PC1 loading
#      Left  = lowest loading
#      Right = highest loading = leading edge
#
#   2. Builds an empirical variance curve along EVS rank
#
#   3. Defines the leading edge a priori as the TOP 5000 ranked features
#      and fixes that as the manuscript cutoff center
#
#   4. Uses the smoothed second derivative only to identify a compact local
#      corroboration band around that manual cutoff
#
#   5. Uses NB-derived quantities only to define a cutoff range around the
#      manual cutoff, not to define the cutoff itself
#
#   6. Makes the leading-edge NB2 evidence explicit by calculating and plotting
#      leading-edge-versus-remainder corroboration statistics:
#         - median smoothed NB2 in leading edge and remainder
#         - median smoothed NB2-NB1 contrast in leading edge and remainder
#         - median smoothed log(alpha*mu) in leading edge and remainder
#         - median combined NB support in leading edge and remainder
#
# Outputs
#   exports/variance_derivative_nb_range_manual5000/
#     <comparison>_cutoff_folder/
#       <comparison>_control_rank_panel.png
#       <comparison>_treatment_rank_panel.png
#       <comparison>_cutoff_range_panel.png
#       <comparison>_feature_level_metrics.csv
#       <comparison>_control_rank_series.csv
#       <comparison>_treatment_rank_series.csv
#       <comparison>_control_zero_crossings_all.csv
#       <comparison>_treatment_zero_crossings_all.csv
#       <comparison>_control_selected_zero_crossings.csv
#       <comparison>_treatment_selected_zero_crossings.csv
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
out_root <- file.path(repo_dir, "exports", "variance_derivative_nb_range_manual5000")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

manual_leading_edge_n <- 26757L

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

manual_cutoff_search_left_fraction <- 0.12
manual_cutoff_search_right_fraction <- 0.12
manual_cutoff_pair_max_width_fraction <- 0.06

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

select_manual_cutoff_zero_pair <- function(df,
                                           manual_cutoff_rank,
                                           search_left_fraction = 0.12,
                                           search_right_fraction = 0.12,
                                           max_width_fraction = 0.06) {
  n <- nrow(df)

  zero_idx <- find_d2_zero_crossings(df$d2_sm)
  zero_tbl <- summarize_zero_crossings(df, zero_idx)

  if (!nrow(zero_tbl)) {
    return(list(
      cutoff_center_index = manual_cutoff_rank,
      cutoff_center_rank = df$rank[manual_cutoff_rank],
      terminal_start_index = manual_cutoff_rank,
      terminal_start_rank = df$rank[manual_cutoff_rank],
      zero_crossings = zero_tbl,
      selected_zero_crossings = data.frame(),
      mode = "fallback: no d2 zero-crossings found near manual cutoff"
    ))
  }

  left_bound <- max(1L, floor(manual_cutoff_rank - n * search_left_fraction))
  right_bound <- min(n, ceiling(manual_cutoff_rank + n * search_right_fraction))
  max_width <- max(25L, floor(n * max_width_fraction))

  local_tbl <- zero_tbl %>%
    filter(index >= left_bound, index <= right_bound)

  if (!nrow(local_tbl)) {
    local_tbl <- zero_tbl
  }

  left_tbl <- local_tbl %>%
    filter(index <= manual_cutoff_rank) %>%
    arrange(desc(index))

  right_tbl <- local_tbl %>%
    filter(index >= manual_cutoff_rank) %>%
    arrange(index)

  if (!nrow(left_tbl) && nrow(zero_tbl)) {
    left_tbl <- zero_tbl %>% filter(index < manual_cutoff_rank) %>% arrange(desc(index))
  }
  if (!nrow(right_tbl) && nrow(zero_tbl)) {
    right_tbl <- zero_tbl %>% filter(index > manual_cutoff_rank) %>% arrange(index)
  }

  if (!nrow(left_tbl) || !nrow(right_tbl)) {
    return(list(
      cutoff_center_index = manual_cutoff_rank,
      cutoff_center_rank = df$rank[manual_cutoff_rank],
      terminal_start_index = manual_cutoff_rank,
      terminal_start_rank = df$rank[manual_cutoff_rank],
      zero_crossings = zero_tbl,
      selected_zero_crossings = data.frame(),
      mode = "fallback: unable to flank manual cutoff with local zero-crossings"
    ))
  }

  pair_tbl <- do.call(
    rbind,
    lapply(seq_len(nrow(left_tbl)), function(i) {
      do.call(
        rbind,
        lapply(seq_len(nrow(right_tbl)), function(j) {
          left_i <- left_tbl$index[i]
          right_i <- right_tbl$index[j]
          width <- right_i - left_i
          if (width <= 0) return(NULL)

          data.frame(
            left_index = left_i,
            right_index = right_i,
            left_rank = df$rank[left_i],
            right_rank = df$rank[right_i],
            width = width,
            midpoint = (df$rank[left_i] + df$rank[right_i]) / 2,
            midpoint_distance_to_manual = abs(((df$rank[left_i] + df$rank[right_i]) / 2) - df$rank[manual_cutoff_rank]),
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )

  if (is.null(pair_tbl) || !nrow(pair_tbl)) {
    return(list(
      cutoff_center_index = manual_cutoff_rank,
      cutoff_center_rank = df$rank[manual_cutoff_rank],
      terminal_start_index = manual_cutoff_rank,
      terminal_start_rank = df$rank[manual_cutoff_rank],
      zero_crossings = zero_tbl,
      selected_zero_crossings = data.frame(),
      mode = "fallback: zero-crossing pair table empty"
    ))
  }

  pair_tbl_narrow <- pair_tbl %>% filter(width <= max_width)
  if (!nrow(pair_tbl_narrow)) pair_tbl_narrow <- pair_tbl

  best <- pair_tbl_narrow %>%
    arrange(midpoint_distance_to_manual, width) %>%
    slice(1)

  selected_tbl <- data.frame(
    index = c(best$left_index, best$right_index),
    rank = c(best$left_rank, best$right_rank),
    role = c("left_local_zero", "right_local_zero"),
    stringsAsFactors = FALSE
  )

  list(
    cutoff_center_index = manual_cutoff_rank,
    cutoff_center_rank = df$rank[manual_cutoff_rank],
    terminal_start_index = best$right_index,
    terminal_start_rank = best$right_rank,
    zero_crossings = zero_tbl,
    selected_zero_crossings = selected_tbl,
    mode = "manual cutoff fixed at top 5000 ranks; local d2 zero-crossings used only to corroborate that boundary"
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

compute_leading_edge_corrob <- function(df, cutoff_center_index) {
  leading_idx <- seq_len(cutoff_center_index)
  remainder_idx <- seq.int(cutoff_center_index + 1L, nrow(df))
  if (!length(remainder_idx)) remainder_idx <- nrow(df)

  med_le_nb2 <- safe_median(df$log_nb2[leading_idx])
  med_re_nb2 <- safe_median(df$log_nb2[remainder_idx])

  med_le_gap <- safe_median(df$nb_gap_sm[leading_idx])
  med_re_gap <- safe_median(df$nb_gap_sm[remainder_idx])

  med_le_amu <- safe_median(df$amu_sm[leading_idx])
  med_re_amu <- safe_median(df$amu_sm[remainder_idx])

  med_le_support <- safe_median(df$nb_support[leading_idx])
  med_re_support <- safe_median(df$nb_support[remainder_idx])

  data.frame(
    leading_edge_median_log_nb2 = med_le_nb2,
    remainder_median_log_nb2 = med_re_nb2,
    leading_minus_remainder_log_nb2 = med_le_nb2 - med_re_nb2,
    leading_remainder_log_nb2_ratio = safe_ratio(med_le_nb2, med_re_nb2),

    leading_edge_median_nb_gap = med_le_gap,
    remainder_median_nb_gap = med_re_gap,
    leading_minus_remainder_nb_gap = med_le_gap - med_re_gap,
    leading_remainder_nb_gap_ratio = safe_ratio(med_le_gap, med_re_gap),

    leading_edge_median_log_alpha_mu = med_le_amu,
    remainder_median_log_alpha_mu = med_re_amu,
    leading_minus_remainder_log_alpha_mu = med_le_amu - med_re_amu,
    leading_remainder_log_alpha_mu_ratio = safe_ratio(med_le_amu, med_re_amu),

    leading_edge_median_nb_support = med_le_support,
    remainder_median_nb_support = med_re_support,
    leading_minus_remainder_nb_support = med_le_support - med_re_support,
    leading_remainder_nb_support_ratio = safe_ratio(med_le_support, med_re_support),

    leading_edge_more_nb2 = med_le_nb2 > med_re_nb2,
    stringsAsFactors = FALSE
  )
}

select_cutoff_derivative_nb_range <- function(rank_df,
                                              manual_leading_edge_n = 5000L,
                                              spline_spar = 0.60,
                                              nb_range_drop_fraction = 0.70,
                                              nb_range_max_span_fraction = 0.18,
                                              nb_smooth_window = 151L,
                                              manual_cutoff_search_left_fraction = 0.12,
                                              manual_cutoff_search_right_fraction = 0.12,
                                              manual_cutoff_pair_max_width_fraction = 0.06) {
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

  manual_cutoff_index <- min(manual_leading_edge_n, n)

  nb_obj <- compute_nb_support_score(df, nb_smooth_window = nb_smooth_window)
  df$nb_gap_sm <- nb_obj$nb_gap_sm
  df$amu_sm <- nb_obj$amu_sm
  df$nb_support <- nb_obj$nb_support

  geom_obj <- select_manual_cutoff_zero_pair(
    df = df,
    manual_cutoff_rank = manual_cutoff_index,
    search_left_fraction = manual_cutoff_search_left_fraction,
    search_right_fraction = manual_cutoff_search_right_fraction,
    max_width_fraction = manual_cutoff_pair_max_width_fraction
  )

  nb_rng <- expand_nb_range_around_center(
    df = df,
    center_idx = geom_obj$cutoff_center_index,
    terminal_start_idx = geom_obj$terminal_start_index,
    nb_support = df$nb_support,
    nb_range_drop_fraction = nb_range_drop_fraction,
    nb_range_max_span_fraction = nb_range_max_span_fraction
  )

  corrob <- compute_leading_edge_corrob(df, geom_obj$cutoff_center_index)

  total_features <- n
  pre_evs_leading_edge_size <- geom_obj$cutoff_center_index
  pre_evs_remainder_size <- total_features - geom_obj$cutoff_center_index

  list(
    curve_df = df,
    zero_crossings = geom_obj$zero_crossings,
    selected_zero_crossings = geom_obj$selected_zero_crossings,
    mode = geom_obj$mode,

    terminal_start_index = geom_obj$terminal_start_index,
    terminal_start_rank = geom_obj$terminal_start_rank,
    terminal_end_index = geom_obj$terminal_start_index,
    terminal_end_rank = geom_obj$terminal_start_rank,

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
  corr <- onset_info$corrob

  x_rng <- range(df$rank, na.rm = TRUE)
  x_left <- x_rng[1] + 0.06 * diff(x_rng)
  x_right <- x_rng[1] + 0.72 * diff(x_rng)

  p1 <- ggplot(df, aes(rank, abs_loading)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    annotate("rect", xmin = range_min_x, xmax = range_max_x, ymin = -Inf, ymax = Inf, alpha = 0.08) +
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
        "\nManual leading edge = top ", onset_info$pre_evs_leading_edge_size,
        " ranks",
        "\nCutoff center rank = ", center_x,
        "\nLocal corroboration band = [", range_min_x, ", ", range_max_x, "]",
        "\nRight local zero = ", terminal_x,
        "\nRemainder size = ", onset_info$pre_evs_remainder_size
      )
    ) +
    labs(
      title = paste0(title_prefix, ": absolute PC1 loading series"),
      subtitle = "Leading edge is on the RIGHT and is fixed a priori at the top 5000 ranks",
      x = "EVS rank",
      y = "|PC1 loading|"
    ) +
    theme_bw(base_size = 10)

  p2 <- ggplot(df, aes(rank, var_fit)) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    annotate("rect", xmin = range_min_x, xmax = range_max_x, ymin = -Inf, ymax = Inf, alpha = 0.08) +
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
        "\nThe cutoff itself is fixed manually at the top 5000 ranks",
        "\nThe smoothed second derivative is used only to",
        "\nfind a compact local corroboration band around that cutoff"
      )
    ) +
    labs(
      title = paste0(title_prefix, ": smoothed empirical variance curve"),
      subtitle = "Selected local zero-crossings are shown only as corroboration around the manual cutoff",
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    theme_bw(base_size = 10)

  p3 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = d1_sm, color = "Smoothed slope d1"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(aes(y = d2_sm, color = "Smoothed curvature d2"), linewidth = 1.0, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    annotate("rect", xmin = range_min_x, xmax = range_max_x, ymin = -Inf, ymax = Inf, alpha = 0.08) +
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
        "\nCenter = fixed manual cutoff rank",
        "\nRight local zero = nearest corroborating d2 zero to the right",
        "\nThe zero pair supports the cutoff but does not redefine it"
      )
    ) +
    labs(
      title = paste0(title_prefix, ": derivative support"),
      subtitle = "The local d2 zero-crossings are descriptive support for the manual cutoff",
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
        "\nThe top 5000 ranks define the manual leading edge",
        "\nLeading-edge med log(NB2) = ", round(corr$leading_edge_median_log_nb2, 3),
        "\nRemainder med log(NB2) = ", round(corr$remainder_median_log_nb2, 3),
        "\nLeading-minus-remainder log(NB2) = ", round(corr$leading_minus_remainder_log_nb2, 3),
        "\nLeading-edge med NB2-NB1 = ", round(corr$leading_edge_median_nb_gap, 3),
        "\nRemainder med NB2-NB1 = ", round(corr$remainder_median_nb_gap, 3),
        "\nLeading-minus-remainder NB2-NB1 = ", round(corr$leading_minus_remainder_nb_gap, 3),
        "\nLeading-edge med log(alpha*mu) = ", round(corr$leading_edge_median_log_alpha_mu, 3),
        "\nRemainder med log(alpha*mu) = ", round(corr$remainder_median_log_alpha_mu, 3),
        "\nLeading-minus-remainder log(alpha*mu) = ", round(corr$leading_minus_remainder_log_alpha_mu, 3),
        "\nLeading-edge more NB2 = ", corr$leading_edge_more_nb2,
        "\nCenter NB support = ", round(onset_info$center_nb_support, 3),
        "\nNB threshold = ", round(onset_info$nb_support_threshold, 3)
      )
    ) +
    labs(
      title = paste0(title_prefix, ": NB1 / NB2 / alpha*mu support"),
      subtitle = "Leading-edge elevation in NB2 and alpha*mu supports the top-5000 interpretation",
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
        "\nManual cutoff = ", center_x,
        "\nBand = [", range_min_x, ", ", range_max_x, "]",
        "\nRight local zero = ", terminal_x
      )
    ) +
    labs(
      title = paste0(title_prefix, ": NB-supported cutoff band"),
      subtitle = "The band expands around the manual cutoff while combined NB support stays elevated",
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
  x_left <- min(x_all, na.rm = TRUE) + 0.06 * diff(range(x_all, na.rm = TRUE))

  p <- ggplot() +
    geom_line(data = ctrl_df, aes(rank, var_fit, color = "Control variance fit"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(data = trt_df, aes(rank, var_fit, color = "Treatment variance fit"), linewidth = 1.0, na.rm = TRUE) +
    geom_vline(xintercept = ctrl_onset$cutoff_center_rank, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = trt_onset$cutoff_center_rank, linetype = 3, linewidth = 0.8) +
    annotate(
      "rect",
      xmin = min(ctrl_onset$cutoff_range_rank_min, trt_onset$cutoff_range_rank_min),
      xmax = max(ctrl_onset$cutoff_range_rank_max, trt_onset$cutoff_range_rank_max),
      ymin = -Inf, ymax = Inf,
      alpha = 0.08
    ) +
    annotate(
      "label",
      x = x_left,
      y = max(c(ctrl_df$var_fit, trt_df$var_fit), na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label = paste0(
        "Comparison-level view",
        "\nBoth panels use the same fixed top-5000 manuscript cutoff",
        "\nControl center = ", ctrl_onset$cutoff_center_rank,
        "\nTreatment center = ", trt_onset$cutoff_center_rank,
        "\nControl local band = [", ctrl_onset$cutoff_range_rank_min, ", ", ctrl_onset$cutoff_range_rank_max, "]",
        "\nTreatment local band = [", trt_onset$cutoff_range_rank_min, ", ", trt_onset$cutoff_range_rank_max, "]"
      )
    ) +
    labs(
      title = paste0(title_prefix, ": control and treatment local corroboration bands"),
      subtitle = "The manuscript cutoff is fixed; the local derivative geometry is shown only as support",
      x = "EVS rank",
      y = "Fitted log(1 + variance)",
      color = NULL
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

count_file_resolved <- resolve_counts_file(count_file)
message("Using count file: ", count_file_resolved)

read_obj <- read_count_matrix(count_file_resolved, meta_all$id)
count_matrix <- read_obj$count_matrix
annot_df <- read_obj$annot_df

overall_summary_list <- list()

for (ii in seq_len(nrow(comparison_table))) {
  cmp <- comparison_table[ii, , drop = FALSE]
  cmp_name <- cmp$comparison_name

  message("Processing comparison: ", cmp_name)

  cmp_dir <- file.path(out_root, paste0(cmp_name, "_cutoff_folder"))
  dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

  sub_obj <- subset_comparison(count_matrix, cmp, meta_all)
  cmp_counts <- sub_obj$count_matrix

  evs_tbl <- build_evs_table(
    count_mat = cmp_counts,
    ctrl_cols = sub_obj$ctrl_ids,
    trt_cols = sub_obj$trt_ids,
    feature_ids = rownames(cmp_counts)
  )

  feat_tbl <- compute_feature_metrics_empirical(
    count_mat = cmp_counts,
    ctrl_cols = sub_obj$ctrl_ids,
    trt_cols = sub_obj$trt_ids,
    feature_ids = rownames(cmp_counts)
  )

  full_tbl <- evs_tbl %>%
    left_join(feat_tbl, by = "feature_id") %>%
    left_join(annot_df, by = "feature_id") %>%
    mutate(gene_symbol = ifelse(is.na(gene_symbol) | !nzchar(gene_symbol), feature_id, gene_symbol))

  utils::write.csv(
    full_tbl,
    file = file.path(cmp_dir, paste0(cmp_name, "_feature_level_metrics.csv")),
    row.names = FALSE
  )

  ctrl_rank_df <- build_rank_series(full_tbl, "control")
  trt_rank_df  <- build_rank_series(full_tbl, "treatment")

  utils::write.csv(
    ctrl_rank_df,
    file = file.path(cmp_dir, paste0(cmp_name, "_control_rank_series.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    trt_rank_df,
    file = file.path(cmp_dir, paste0(cmp_name, "_treatment_rank_series.csv")),
    row.names = FALSE
  )

  ctrl_onset <- select_cutoff_derivative_nb_range(
    ctrl_rank_df,
    manual_leading_edge_n = manual_leading_edge_n,
    spline_spar = spline_spar,
    nb_range_drop_fraction = nb_range_drop_fraction,
    nb_range_max_span_fraction = nb_range_max_span_fraction,
    nb_smooth_window = nb_smooth_window,
    manual_cutoff_search_left_fraction = manual_cutoff_search_left_fraction,
    manual_cutoff_search_right_fraction = manual_cutoff_search_right_fraction,
    manual_cutoff_pair_max_width_fraction = manual_cutoff_pair_max_width_fraction
  )

  trt_onset <- select_cutoff_derivative_nb_range(
    trt_rank_df,
    manual_leading_edge_n = manual_leading_edge_n,
    spline_spar = spline_spar,
    nb_range_drop_fraction = nb_range_drop_fraction,
    nb_range_max_span_fraction = nb_range_max_span_fraction,
    nb_smooth_window = nb_smooth_window,
    manual_cutoff_search_left_fraction = manual_cutoff_search_left_fraction,
    manual_cutoff_search_right_fraction = manual_cutoff_search_right_fraction,
    manual_cutoff_pair_max_width_fraction = manual_cutoff_pair_max_width_fraction
  )

  utils::write.csv(
    ctrl_onset$zero_crossings,
    file = file.path(cmp_dir, paste0(cmp_name, "_control_zero_crossings_all.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    trt_onset$zero_crossings,
    file = file.path(cmp_dir, paste0(cmp_name, "_treatment_zero_crossings_all.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    ctrl_onset$selected_zero_crossings,
    file = file.path(cmp_dir, paste0(cmp_name, "_control_selected_zero_crossings.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    trt_onset$selected_zero_crossings,
    file = file.path(cmp_dir, paste0(cmp_name, "_treatment_selected_zero_crossings.csv")),
    row.names = FALSE
  )

  build_dataset_panel(
    rank_df = ctrl_rank_df,
    onset_info = ctrl_onset,
    title_prefix = paste0(cmp_name, " control"),
    out_file = file.path(cmp_dir, paste0(cmp_name, "_control_rank_panel.png"))
  )

  build_dataset_panel(
    rank_df = trt_rank_df,
    onset_info = trt_onset,
    title_prefix = paste0(cmp_name, " treatment"),
    out_file = file.path(cmp_dir, paste0(cmp_name, "_treatment_rank_panel.png"))
  )

  build_range_panel(
    ctrl_onset = ctrl_onset,
    trt_onset = trt_onset,
    title_prefix = cmp_name,
    out_file = file.path(cmp_dir, paste0(cmp_name, "_cutoff_range_panel.png"))
  )

  ctrl_corr <- ctrl_onset$corrob
  trt_corr <- trt_onset$corrob

  summary_df <- data.frame(
    comparison = cmp_name,

    control_mode = ctrl_onset$mode,
    control_cutoff_center_rank = ctrl_onset$cutoff_center_rank,
    control_terminal_start_rank = ctrl_onset$terminal_start_rank,
    control_cutoff_range_rank_min = ctrl_onset$cutoff_range_rank_min,
    control_cutoff_range_rank_max = ctrl_onset$cutoff_range_rank_max,
    control_pre_evs_leading_edge_size = ctrl_onset$pre_evs_leading_edge_size,
    control_pre_evs_remainder_size = ctrl_onset$pre_evs_remainder_size,
    control_center_nb_support = ctrl_onset$center_nb_support,
    control_nb_support_threshold = ctrl_onset$nb_support_threshold,

    control_leading_edge_median_log_nb2 = ctrl_corr$leading_edge_median_log_nb2,
    control_remainder_median_log_nb2 = ctrl_corr$remainder_median_log_nb2,
    control_leading_minus_remainder_log_nb2 = ctrl_corr$leading_minus_remainder_log_nb2,
    control_leading_remainder_log_nb2_ratio = ctrl_corr$leading_remainder_log_nb2_ratio,

    control_leading_edge_median_nb_gap = ctrl_corr$leading_edge_median_nb_gap,
    control_remainder_median_nb_gap = ctrl_corr$remainder_median_nb_gap,
    control_leading_minus_remainder_nb_gap = ctrl_corr$leading_minus_remainder_nb_gap,
    control_leading_remainder_nb_gap_ratio = ctrl_corr$leading_remainder_nb_gap_ratio,

    control_leading_edge_median_log_alpha_mu = ctrl_corr$leading_edge_median_log_alpha_mu,
    control_remainder_median_log_alpha_mu = ctrl_corr$remainder_median_log_alpha_mu,
    control_leading_minus_remainder_log_alpha_mu = ctrl_corr$leading_minus_remainder_log_alpha_mu,
    control_leading_remainder_log_alpha_mu_ratio = ctrl_corr$leading_remainder_log_alpha_mu_ratio,

    control_leading_edge_median_nb_support = ctrl_corr$leading_edge_median_nb_support,
    control_remainder_median_nb_support = ctrl_corr$remainder_median_nb_support,
    control_leading_minus_remainder_nb_support = ctrl_corr$leading_minus_remainder_nb_support,
    control_leading_remainder_nb_support_ratio = ctrl_corr$leading_remainder_nb_support_ratio,

    control_leading_edge_more_nb2 = ctrl_corr$leading_edge_more_nb2,

    treatment_mode = trt_onset$mode,
    treatment_cutoff_center_rank = trt_onset$cutoff_center_rank,
    treatment_terminal_start_rank = trt_onset$terminal_start_rank,
    treatment_cutoff_range_rank_min = trt_onset$cutoff_range_rank_min,
    treatment_cutoff_range_rank_max = trt_onset$cutoff_range_rank_max,
    treatment_pre_evs_leading_edge_size = trt_onset$pre_evs_leading_edge_size,
    treatment_pre_evs_remainder_size = trt_onset$pre_evs_remainder_size,
    treatment_center_nb_support = trt_onset$center_nb_support,
    treatment_nb_support_threshold = trt_onset$nb_support_threshold,

    treatment_leading_edge_median_log_nb2 = trt_corr$leading_edge_median_log_nb2,
    treatment_remainder_median_log_nb2 = trt_corr$remainder_median_log_nb2,
    treatment_leading_minus_remainder_log_nb2 = trt_corr$leading_minus_remainder_log_nb2,
    treatment_leading_remainder_log_nb2_ratio = trt_corr$leading_remainder_log_nb2_ratio,

    treatment_leading_edge_median_nb_gap = trt_corr$leading_edge_median_nb_gap,
    treatment_remainder_median_nb_gap = trt_corr$remainder_median_nb_gap,
    treatment_leading_minus_remainder_nb_gap = trt_corr$leading_minus_remainder_nb_gap,
    treatment_leading_remainder_nb_gap_ratio = trt_corr$leading_remainder_nb_gap_ratio,

    treatment_leading_edge_median_log_alpha_mu = trt_corr$leading_edge_median_log_alpha_mu,
    treatment_remainder_median_log_alpha_mu = trt_corr$remainder_median_log_alpha_mu,
    treatment_leading_minus_remainder_log_alpha_mu = trt_corr$leading_minus_remainder_log_alpha_mu,
    treatment_leading_remainder_log_alpha_mu_ratio = trt_corr$leading_remainder_log_alpha_mu_ratio,

    treatment_leading_edge_median_nb_support = trt_corr$leading_edge_median_nb_support,
    treatment_remainder_median_nb_support = trt_corr$remainder_median_nb_support,
    treatment_leading_minus_remainder_nb_support = trt_corr$leading_minus_remainder_nb_support,
    treatment_leading_remainder_nb_support_ratio = trt_corr$leading_remainder_nb_support_ratio,

    treatment_leading_edge_more_nb2 = trt_corr$leading_edge_more_nb2,

    stringsAsFactors = FALSE
  )

  utils::write.csv(
    summary_df,
    file = file.path(cmp_dir, paste0(cmp_name, "_cutoff_summary.csv")),
    row.names = FALSE
  )

  overall_summary_list[[cmp_name]] <- summary_df
}

overall_cutoff_summary <- bind_rows(overall_summary_list)
utils::write.csv(
  overall_cutoff_summary,
  file = file.path(out_root, "overall_cutoff_summary.csv"),
  row.names = FALSE
)

message("Done.")
