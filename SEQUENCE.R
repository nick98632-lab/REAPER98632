# =============================================================================
# SEQUENCE STAGE 1 FINAL PRE EVS CUTOFF SELECTOR
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# Features are ranked by absolute PC1 loading from lowest loading on the left to
# highest loading on the right. The empirical variance curve is smoothed along
# this rank axis, and its second derivative is used to identify terminal
# zero-crossings on the right-hand side of the ranked series. The cutoff
# interval is defined by the last two terminal zero-crossings of the smoothed
# second derivative. The earlier zero defines the cutoff anchor and the later
# zero defines the terminal-start anchor. No minimum run-length rule, slope
# quantile threshold, or directional sign constraint on the first derivative is
# imposed. Negative-binomial support quantities are then used only to corroborate
# the derivative-defined transition and to summarize whether the region to the
# right of the cutoff behaves more like a higher-variance leading-edge regime.
#
# NB comparison rule
# The right-hand region is defined as all genes from the cutoff anchor to the end
# of the ranked series. Let that region contain m genes. The comparison left-hand
# region is defined as the same number m of genes immediately to the left of the
# cutoff anchor, truncated at the left boundary if needed. Smoothed NB2,
# smoothed NB2 minus NB1 contrast, and smoothed log(alpha*mu) are compared
# between these matched windows.
#
# Figure design
# Method descriptions are placed on the left side of panels. Cutoff markers and
# numeric labels are placed on the right side. The earlier terminal zero is
# marked by a filled circle, the later terminal zero by an open circle, and the
# fixed rank 5000 reference is marked by a diamond. The interval between the two
# anchors is shaded. All cutoffs are marked by vertical lines.
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
#
# MANUSCRIPT METHODS DESCRIPTION
# These settings control smoothing of the empirical variance curve, the right-side
# terminal region in which the final two second-derivative zero-crossings are
# sought, the smoothing of NB support quantities, and the formatting of exported
# outputs. No threshold is used to determine whether a zero-crossing is valid
# beyond a very small numeric tolerance for floating-point stability.
# =============================================================================

repo_dir <- getwd()
count_file <- file.path(repo_dir, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
out_root <- file.path(repo_dir, "exports", "variance_derivative_nb_range_final_clean")
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

# smoothing settings
spline_spar <- 0.60
nb_smooth_window <- 151L

# terminal region definition
terminal_fraction_start <- 0.70

# numeric zero tolerance
zero_tol <- 1e-6

# NB band construction
nb_range_drop_fraction <- 0.70
nb_range_max_span_fraction <- 0.18

# fixed reference marker
fixed_rank_reference <- 5000L

# =============================================================================
# BASIC HELPERS
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# These helpers perform robust numeric summaries, file resolution, rescaling, and
# rolling median smoothing. Rolling medians are used for support summaries to
# reduce the influence of isolated spikes. Zero-crossings are detected from sign
# changes after applying a minimal numeric tolerance around zero.
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
  stats::median(x, na.rm = TRUE)
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

# =============================================================================
# DATA INGESTION
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# The raw feature by sample count matrix is read from disk. Features with missing
# identifiers or non-finite values are removed. Duplicate feature IDs are
# collapsed by summation. For each comparison, treatment and control samples are
# subset directly from the count matrix while preserving feature annotations.
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
# Features are ranked separately within each arm by absolute PC1 loading computed
# from log2(count+1) transformed counts. For each feature and arm, empirical mean
# and variance are calculated, together with NB1 as the mean-like component, NB2
# as the excess variance above the mean, and alpha times mu. These empirical
# quantities provide interpretable support for whether the right-hand region of
# the ranked series behaves more like a higher variance leading edge.
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

# =============================================================================
# DERIVATIVE GEOMETRY
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# The empirical variance curve is smoothed by smoothing spline and differentiated
# analytically to obtain first and second derivatives along the full rank axis.
# Second-derivative zero-crossings are defined by sign changes after applying a
# minimal numerical tolerance around zero. The cutoff interval is then defined by
# the last two zero-crossings within the right-hand terminal portion of the rank
# axis. No threshold on first-derivative magnitude, no minimum run length, and no
# directional rule for whether the curve must cross upward or downward is used.
# =============================================================================

find_d2_zero_crossings <- function(d2_vec, zero_tol = 1e-6) {
  d2_vec <- as.numeric(d2_vec)
  n <- length(d2_vec)
  out <- integer(0)
  if (n < 2L) return(out)

  x <- d2_vec
  x[!is.finite(x)] <- NA_real_
  x[is.finite(x) & abs(x) <= zero_tol] <- 0

  s <- sign(x)
  s[!is.finite(s)] <- 0

  for (i in 2:n) {
    if (!is.finite(x[i - 1L]) || !is.finite(x[i])) next

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
      var_fit_here = df$var_fit[i],
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out)
}

choose_last_two_terminal_zeros <- function(df,
                                           terminal_fraction_start = 0.70,
                                           zero_tol = 1e-6) {
  n <- nrow(df)

  zero_idx <- find_d2_zero_crossings(df$d2_sm, zero_tol = zero_tol)
  zero_tbl <- summarize_zero_crossings(df, zero_idx)

  if (!nrow(zero_tbl)) {
    fallback2 <- n
    fallback1 <- max(1L, n - 1L)
    sel_tbl <- data.frame(
      index = c(fallback1, fallback2),
      rank = df$rank[c(fallback1, fallback2)],
      role = c("cutoff_anchor", "terminal_anchor"),
      stringsAsFactors = FALSE
    )
    return(list(
      zero_crossings = zero_tbl,
      selected_zero_crossings = sel_tbl,
      cutoff_anchor_index = fallback1,
      cutoff_anchor_rank = df$rank[fallback1],
      terminal_anchor_index = fallback2,
      terminal_anchor_rank = df$rank[fallback2],
      mode = "fallback: no detected terminal zero-crossings"
    ))
  }

  terminal_start_idx <- max(1L, floor(terminal_fraction_start * n))
  terminal_tbl <- zero_tbl[zero_tbl$index >= terminal_start_idx, , drop = FALSE]

  if (nrow(terminal_tbl) < 2L) {
    terminal_tbl <- tail(zero_tbl, min(2L, nrow(zero_tbl)))
  } else {
    terminal_tbl <- tail(terminal_tbl, 2L)
  }

  if (nrow(terminal_tbl) == 1L) {
    idx2 <- terminal_tbl$index[1]
    idx1 <- max(1L, idx2 - 1L)
    terminal_tbl <- data.frame(
      index = c(idx1, idx2),
      rank = df$rank[c(idx1, idx2)],
      stringsAsFactors = FALSE
    )
  }

  terminal_tbl$role <- c("cutoff_anchor", "terminal_anchor")

  list(
    zero_crossings = zero_tbl,
    selected_zero_crossings = terminal_tbl,
    cutoff_anchor_index = terminal_tbl$index[1],
    cutoff_anchor_rank = terminal_tbl$rank[1],
    terminal_anchor_index = terminal_tbl$index[2],
    terminal_anchor_rank = terminal_tbl$rank[2],
    mode = "last two terminal d2 zero-crossings on the right-hand side"
  )
}

# =============================================================================
# NB SUPPORT AND MATCHED LEFT RIGHT COMPARISON
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# NB support is constructed from smoothed log(NB2+1), smoothed NB2 minus NB1
# contrast, and smoothed log(alpha*mu). The cutoff interval itself is defined by
# the derivatives, not by these NB quantities. For corroboration, all genes from
# the cutoff anchor to the end of the ranked series define the right-hand region.
# If that right-hand region contains m genes, the left-hand comparison region is
# defined as the same number m of genes immediately to the left of the cutoff.
# Median smoothed NB quantities are then compared between these matched windows.
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
# This function applies the full method to a single arm. It smooths the
# empirical variance curve, computes derivatives, identifies the last two
# terminal second-derivative zero-crossings on the right-hand side, constructs
# an NB-supported local band around the cutoff anchor, and computes the matched
# left-right corroboration summaries used in the figures and output tables.
# =============================================================================

select_cutoff_derivative_nb_range <- function(rank_df,
                                              spline_spar = 0.60,
                                              terminal_fraction_start = 0.70,
                                              zero_tol = 1e-6,
                                              nb_range_drop_fraction = 0.70,
                                              nb_range_max_span_fraction = 0.18,
                                              nb_smooth_window = 151L) {
  df <- rank_df
  n <- nrow(df)
  x <- df$rank
  y <- df$log_variance

  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10L) stop("Not enough finite variance points for smoothing.", call. = FALSE)

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

  geom_obj <- choose_last_two_terminal_zeros(
    df = df,
    terminal_fraction_start = terminal_fraction_start,
    zero_tol = zero_tol
  )

  nb_rng <- expand_nb_range_around_center(
    df = df,
    center_idx = geom_obj$cutoff_anchor_index,
    terminal_start_idx = geom_obj$terminal_anchor_index,
    nb_support = df$nb_support,
    nb_range_drop_fraction = nb_range_drop_fraction,
    nb_range_max_span_fraction = nb_range_max_span_fraction
  )

  corrob <- compute_left_right_corrob(df, geom_obj$cutoff_anchor_index)

  total_features <- n
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

    nb_band_index_min = nb_rng$range_left_index,
    nb_band_rank_min = nb_rng$range_left_rank,
    nb_band_index_max = nb_rng$range_right_index,
    nb_band_rank_max = nb_rng$range_right_rank,

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
#
# MANUSCRIPT METHODS DESCRIPTION
# The figure builder places concise method notes on the left side of each panel
# so that the cutoff markers remain visually clear on the right. The earlier
# anchor is shown with a filled circle and the later terminal-side anchor with
# an open circle. A fixed 5000-rank reference is displayed with a diamond and a
# vertical line. The derivative-defined cutoff interval is shaded, and the
# NB-supported local band around the cutoff anchor is also shown.
# =============================================================================

build_dataset_panel <- function(rank_df, onset_info, title_prefix, out_file, fixed_rank_reference = 5000L) {
  df <- onset_info$curve_df
  cutoff_x <- onset_info$cutoff_center_rank
  terminal_x <- onset_info$terminal_start_rank
  interval_min_x <- onset_info$cutoff_range_rank_min
  interval_max_x <- onset_info$cutoff_range_rank_max
  nb_band_min_x <- onset_info$nb_band_rank_min
  nb_band_max_x <- onset_info$nb_band_rank_max
  corr <- onset_info$corrob

  ref_x <- min(fixed_rank_reference, max(df$rank, na.rm = TRUE))
  x_rng <- range(df$rank, na.rm = TRUE)
  x_left <- x_rng[1] + 0.06 * diff(x_rng)
  x_right <- x_rng[1] + 0.72 * diff(x_rng)

  cutoff_pt <- df[df$rank == cutoff_x & is.finite(df$var_fit), , drop = FALSE]
  terminal_pt <- df[df$rank == terminal_x & is.finite(df$var_fit), , drop = FALSE]
  ref_pt <- df[df$rank == ref_x & is.finite(df$var_fit), , drop = FALSE]

  p1 <- ggplot(df, aes(rank, abs_loading)) +
    annotate("rect", xmin = interval_min_x, xmax = interval_max_x, ymin = -Inf, ymax = Inf, alpha = 0.10) +
    annotate("rect", xmin = nb_band_min_x, xmax = nb_band_max_x, ymin = -Inf, ymax = Inf, alpha = 0.05) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    geom_vline(xintercept = ref_x, linetype = 4, linewidth = 0.8) +
    annotate(
      "label",
      x = x_left,
      y = max(df$abs_loading, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.7,
      label = paste0(
        "Absolute loading panel",
        "\nInterval = last two terminal d2 zeros on the right side",
        "\nFilled marker = earlier cutoff anchor",
        "\nOpen marker = later terminal anchor",
        "\nDiamond = fixed rank 5000 reference"
      )
    ) +
    annotate(
      "label",
      x = x_right,
      y = max(df$abs_loading, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.7,
      label = paste0(
        "Cutoff anchor rank = ", cutoff_x,
        "\nTerminal anchor rank = ", terminal_x,
        "\nInterval = [", interval_min_x, ", ", interval_max_x, "]",
        "\nNB band = [", nb_band_min_x, ", ", nb_band_max_x, "]"
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
    annotate("rect", xmin = interval_min_x, xmax = interval_max_x, ymin = -Inf, ymax = Inf, alpha = 0.10) +
    annotate("rect", xmin = nb_band_min_x, xmax = nb_band_max_x, ymin = -Inf, ymax = Inf, alpha = 0.05) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    geom_vline(xintercept = ref_x, linetype = 4, linewidth = 0.8) +
    geom_point(data = cutoff_pt, aes(x = rank, y = var_fit), shape = 16, size = 2.6, inherit.aes = FALSE) +
    geom_point(data = terminal_pt, aes(x = rank, y = var_fit), shape = 1, size = 2.8, stroke = 1.0, inherit.aes = FALSE) +
    geom_point(data = ref_pt, aes(x = rank, y = var_fit), shape = 18, size = 2.8, inherit.aes = FALSE) +
    annotate(
      "label",
      x = x_left,
      y = max(df$var_fit, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.6,
      label = paste0(
        "Variance curve method",
        "\nSmooth empirical log(1+variance)",
        "\nTake the last two terminal d2 zero-crossings",
        "\nEarlier zero = cutoff anchor",
        "\nLater zero = terminal anchor"
      )
    ) +
    labs(
      title = paste0(title_prefix, ": smoothed empirical variance curve"),
      subtitle = "Derivative-defined interval and NB-supported band are shown",
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    theme_bw(base_size = 10)

  d2_zero_tbl <- onset_info$selected_zero_crossings
  d2_zero_pts <- df[df$rank %in% d2_zero_tbl$rank, , drop = FALSE]

  p3 <- ggplot(df, aes(rank)) +
    annotate("rect", xmin = interval_min_x, xmax = interval_max_x, ymin = -Inf, ymax = Inf, alpha = 0.10) +
    geom_line(aes(y = d2_sm, color = "Smoothed curvature d2"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(aes(y = d1_sm, color = "Smoothed slope d1"), linewidth = 1.0, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    geom_vline(xintercept = ref_x, linetype = 4, linewidth = 0.8) +
    geom_point(data = data.frame(rank = cutoff_x, d2_sm = df$d2_sm[df$rank == cutoff_x]), aes(x = rank, y = d2_sm), shape = 16, size = 2.6, inherit.aes = FALSE) +
    geom_point(data = data.frame(rank = terminal_x, d2_sm = df$d2_sm[df$rank == terminal_x]), aes(x = rank, y = d2_sm), shape = 1, size = 2.8, stroke = 1.0, inherit.aes = FALSE) +
    annotate(
      "label",
      x = x_left,
      y = max(c(df$d1_sm, df$d2_sm), na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.55,
      label = paste0(
        "Derivative method",
        "\nZero-crossings are sign changes in smoothed d2",
        "\nNo first-derivative threshold is imposed",
        "\nNo minimum run-length rule is imposed",
        "\nOnly the last two terminal zeros are used"
      )
    ) +
    labs(
      title = paste0(title_prefix, ": derivative support"),
      subtitle = "Filled circle = cutoff anchor, open circle = terminal anchor",
      x = "EVS rank",
      y = "Derivative value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p4 <- ggplot(df, aes(rank)) +
    annotate("rect", xmin = interval_min_x, xmax = interval_max_x, ymin = -Inf, ymax = Inf, alpha = 0.10) +
    annotate("rect", xmin = nb_band_min_x, xmax = nb_band_max_x, ymin = -Inf, ymax = Inf, alpha = 0.05) +
    geom_line(aes(y = nb_support, color = "Combined NB support"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(aes(y = log_nb1, color = "NB1 = mu"), linewidth = 0.5, alpha = 0.35, na.rm = TRUE) +
    geom_line(aes(y = log_nb2, color = "NB2 = variance - mu"), linewidth = 0.5, alpha = 0.35, na.rm = TRUE) +
    geom_line(aes(y = amu_sm, color = "Smoothed log(alpha*mu)"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(aes(y = nb_gap_sm, color = "Smoothed log(NB2+1) - log(NB1+1)"), linewidth = 1.0, na.rm = TRUE) +
    geom_hline(yintercept = onset_info$nb_support_threshold, linetype = 3, linewidth = 0.5) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    geom_vline(xintercept = ref_x, linetype = 4, linewidth = 0.8) +
    annotate(
      "label",
      x = x_left,
      y = max(c(df$nb_support, df$log_nb1, df$log_nb2, df$amu_sm, df$nb_gap_sm), na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.35,
      label = paste0(
        "NB corroboration",
        "\nRight side = all genes from cutoff anchor to end",
        "\nLeft side = same number of genes immediately left",
        "\nLeft n = ", corr$left_n,
        "\nRight n = ", corr$right_n,
        "\nRight-left log(NB2) diff = ", round(corr$right_left_log_nb2_diff, 3),
        "\nRight-left NB2-NB1 diff = ", round(corr$right_left_nb_gap_diff, 3),
        "\nRight-left log(alpha*mu) diff = ", round(corr$right_left_log_alpha_mu_diff, 3)
      )
    ) +
    labs(
      title = paste0(title_prefix, ": NB1 / NB2 / alpha*mu support"),
      subtitle = "Matched left-right comparison around the derivative-defined cutoff",
      x = "EVS rank",
      y = "Support value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p5 <- ggplot(df, aes(rank, nb_support)) +
    annotate("rect", xmin = interval_min_x, xmax = interval_max_x, ymin = -Inf, ymax = Inf, alpha = 0.10) +
    annotate("rect", xmin = nb_band_min_x, xmax = nb_band_max_x, ymin = -Inf, ymax = Inf, alpha = 0.05) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    geom_hline(yintercept = onset_info$nb_support_threshold, linetype = 3, linewidth = 0.5) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = terminal_x, linetype = 3, linewidth = 0.8) +
    geom_vline(xintercept = ref_x, linetype = 4, linewidth = 0.8) +
    annotate(
      "label",
      x = x_right,
      y = max(df$nb_support, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.65,
      label = paste0(
        "Final interval and local NB band",
        "\nDerivative interval = [", interval_min_x, ", ", interval_max_x, "]",
        "\nNB band = [", nb_band_min_x, ", ", nb_band_max_x, "]",
        "\nCenter NB support = ", round(onset_info$center_nb_support, 3),
        "\nNB threshold = ", round(onset_info$nb_support_threshold, 3)
      )
    ) +
    labs(
      title = paste0(title_prefix, ": NB-supported cutoff band"),
      subtitle = "The NB band is local support around the derivative-defined cutoff anchor",
      x = "EVS rank",
      y = "Combined NB support"
    ) +
    theme_bw(base_size = 10)

  png(out_file, width = 2400, height = 3200, res = 220)
  gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol = 1)
  dev.off()
}

build_range_panel <- function(ctrl_onset, trt_onset, title_prefix, out_file, fixed_rank_reference = 5000L) {
  ctrl_df <- ctrl_onset$curve_df
  trt_df <- trt_onset$curve_df

  x_all <- c(ctrl_df$rank, trt_df$rank)
  y_all <- c(ctrl_df$var_fit, trt_df$var_fit)
  x_left <- min(x_all, na.rm = TRUE) + 0.06 * diff(range(x_all, na.rm = TRUE))
  ref_x <- min(fixed_rank_reference, max(x_all, na.rm = TRUE))

  p <- ggplot() +
    annotate(
      "rect",
      xmin = min(ctrl_onset$cutoff_range_rank_min, trt_onset$cutoff_range_rank_min),
      xmax = max(ctrl_onset$cutoff_range_rank_max, trt_onset$cutoff_range_rank_max),
      ymin = -Inf, ymax = Inf,
      alpha = 0.08
    ) +
    geom_line(data = ctrl_df, aes(rank, var_fit, color = "Control variance fit"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(data = trt_df, aes(rank, var_fit, color = "Treatment variance fit"), linewidth = 1.0, na.rm = TRUE) +
    geom_vline(xintercept = ctrl_onset$cutoff_center_rank, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = trt_onset$cutoff_center_rank, linetype = 3, linewidth = 0.8) +
    geom_vline(xintercept = ref_x, linetype = 4, linewidth = 0.8) +
    annotate(
      "label",
      x = x_left,
      y = max(y_all, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label = paste0(
        "Treatment / control comparison",
        "\nControl interval = [", ctrl_onset$cutoff_range_rank_min, ", ", ctrl_onset$cutoff_range_rank_max, "]",
        "\nTreatment interval = [", trt_onset$cutoff_range_rank_min, ", ", trt_onset$cutoff_range_rank_max, "]",
        "\nControl anchor = ", ctrl_onset$cutoff_center_rank,
        "\nTreatment anchor = ", trt_onset$cutoff_center_rank
      )
    ) +
    labs(
      title = paste0(title_prefix, ": treatment/control derivative-defined intervals"),
      subtitle = "Both arms use the same final method",
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  png(out_file, width = 2400, height = 1300, res = 220)
  print(p)
  dev.off()
}

# =============================================================================
# MAIN ANALYSIS LOOP
# =============================================================================
#
# MANUSCRIPT METHODS DESCRIPTION
# Each comparison arm is processed independently. Ranked empirical series are
# constructed, derivative-defined terminal anchors are identified, NB-supported
# local bands are generated, and corroboration summaries are exported. All
# zero-crossings, the selected terminal pair, feature-level tables, arm-level
# figures, and global summary tables are written to disk for manuscript use.
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
    nb_range_drop_fraction = nb_range_drop_fraction,
    nb_range_max_span_fraction = nb_range_max_span_fraction,
    nb_smooth_window = nb_smooth_window
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
    nb_range_drop_fraction = nb_range_drop_fraction,
    nb_range_max_span_fraction = nb_range_max_span_fraction,
    nb_smooth_window = nb_smooth_window
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
    file.path(cmp_dir, paste0(cmp_name, "_control_rank_panel.png")),
    fixed_rank_reference = fixed_rank_reference
  )
  build_dataset_panel(
    trt_rank_df,
    trt_onset,
    paste0(cmp_name, " treatment"),
    file.path(cmp_dir, paste0(cmp_name, "_treatment_rank_panel.png")),
    fixed_rank_reference = fixed_rank_reference
  )
  build_range_panel(
    ctrl_onset,
    trt_onset,
    cmp_name,
    file.path(cmp_dir, paste0(cmp_name, "_cutoff_range_panel.png")),
    fixed_rank_reference = fixed_rank_reference
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
    control_nb_band_rank_min = ctrl_onset$nb_band_rank_min,
    control_nb_band_rank_max = ctrl_onset$nb_band_rank_max,
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
    treatment_nb_band_rank_min = trt_onset$nb_band_rank_min,
    treatment_nb_band_rank_max = trt_onset$nb_band_rank_max,
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
  file.path(out_root, "overall_cutoff_summary.csv")),
  row.names = FALSE
)

message("Done. Outputs written to: ", out_root); flush.console()
