# =============================================================================
# SEQUENCE STAGE 1: EMPIRICAL VARIANCE + NB1/NB2 REGIME SPLIT
# RIGHT-SIDE LEADING EDGE, TERMINAL INFLECTION SELECTOR
# -----------------------------------------------------------------------------
# Orientation
#   - lowest absolute loading on the LEFT
#   - highest absolute loading on the RIGHT
#   - leading edge is on the RIGHT
#
# Empirical quantities from raw counts
#   mu        = mean(raw counts)
#   Var       = variance(raw counts)
#   alpha_emp = (Var - mu) / mu^2
#   NB1       = mu
#   NB2       = Var - mu = alpha_emp * mu^2
#   alpha_mu  = alpha_emp * mu = (Var - mu) / mu
#
# Main selector
#   - fit a smooth curve to empirical log(1 + variance)
#   - compute first and second derivatives from the smoothed curve
#   - identify the TERMINAL right-side inflection / steepening region
#   - require NB2 > NB1 support in the same neighborhood
#   - define cutoff as the LEFT EDGE of the strongest sustained terminal block
#
# Output
#   exports/variance_nb1_nb2_terminal_inflection/<comparison>_cutoff_folder/
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

# =============================================================================
# USER SETTINGS
# =============================================================================

repo_dir <- getwd()
input_dir <- file.path(repo_dir, "data")
count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
out_root <- file.path(repo_dir, "exports", "variance_nb1_nb2_terminal_inflection")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

comparison_table <- data.frame(
  comparison_name   = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  treatment_prefix  = c("R0", "R2", "R4", "R8"),
  control_prefix    = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors  = FALSE
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

loess_span <- 0.08
deriv_roll_k <- 151L
left_baseline_fraction <- 0.15
search_fraction_min <- 0.55
search_fraction_max <- 0.995
left_edge_buffer <- 50L
right_edge_buffer <- 20L
min_block_length <- 200L
score_quantile_within_search <- 0.80

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

  stop(paste0("Count file not found. Tried: ", paste(candidates, collapse = ", ")), call. = FALSE)
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
  if (length(sample_cols) == 0L) stop("No count columns matched metadata sample IDs.", call. = FALSE)

  annot_df <- data.frame(
    feature_id = as.character(raw_df[[feature_col]]),
    stringsAsFactors = FALSE
  )
  annot_df$gene_symbol <- if (!is.null(symbol_col)) as.character(raw_df[[symbol_col]]) else annot_df$feature_id

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
      paste0("Missing samples for comparison ", comparison_row$comparison_name, ": ", paste(missing_ids, collapse = ", ")),
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

robust_center_scale <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(list(z = out, center = NA_real_, scale = NA_real_))
  med <- stats::median(x[ok], na.rm = TRUE)
  madv <- stats::mad(x[ok], center = med, constant = 1, na.rm = TRUE)
  if (!is.finite(madv) || madv <= 0) {
    sdv <- stats::sd(x[ok], na.rm = TRUE)
    if (!is.finite(sdv) || sdv <= 0) sdv <- 1
    madv <- sdv
  }
  out[ok] <- (x[ok] - med) / madv
  list(z = out, center = med, scale = madv)
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

  out %>%
    mutate(
      log_variance = ifelse(is.finite(variance) & variance >= 0, log1p(variance), NA_real_),
      log_mu = ifelse(is.finite(mu) & mu > 0, log(mu), NA_real_),
      log_alpha = ifelse(is.finite(alpha) & alpha > 0, log(alpha), NA_real_),
      log_nb1 = ifelse(is.finite(nb1) & nb1 >= 0, log1p(nb1), NA_real_),
      log_nb2 = ifelse(is.finite(nb2) & nb2 >= 0, log1p(nb2), NA_real_),
      nb_gap = log_nb2 - log_nb1,
      log_alpha_mu = ifelse(is.finite(alpha_mu) & alpha_mu > 0, log(alpha_mu), NA_real_)
    )
}

find_best_block <- function(cond, score, search_start, search_end, min_block_length = 250L) {
  cond[is.na(cond)] <- FALSE
  use <- cond[search_start:search_end]
  r <- rle(use)
  ends <- cumsum(r$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)

  block_ids <- which(r$values & r$lengths >= min_block_length)
  if (!length(block_ids)) return(NULL)

  best_mean <- -Inf
  best <- NULL

  for (j in block_ids) {
    lo_local <- starts[j]
    hi_local <- ends[j]
    lo_global <- search_start + lo_local - 1L
    hi_global <- search_start + hi_local - 1L
    m <- mean(score[lo_global:hi_global], na.rm = TRUE)
    if (is.finite(m) && m > best_mean) {
      best_mean <- m
      best <- list(start = lo_global, end = hi_global, mean_score = m)
    }
  }

  best
}

select_terminal_inflection_block <- function(rank_df,
                                             loess_span = 0.08,
                                             deriv_roll_k = 151L,
                                             left_baseline_fraction = 0.15,
                                             search_fraction_min = 0.55,
                                             search_fraction_max = 0.995,
                                             left_edge_buffer = 50L,
                                             right_edge_buffer = 20L,
                                             min_block_length = 200L,
                                             score_quantile_within_search = 0.80) {
  df <- rank_df
  n <- nrow(df)
  x <- df$rank
  y <- df$log_variance

  ok <- is.finite(x) & is.finite(y)
  lo_obj <- loess(y[ok] ~ x[ok], span = loess_span, degree = 2, family = "gaussian", surface = "direct")
  var_fit <- rep(NA_real_, n)
  var_fit[ok] <- stats::predict(lo_obj, newdata = data.frame(`x[ok]` = x[ok]))

  # first derivative
  d1 <- c(NA_real_, diff(var_fit))
  d1_sm <- roll_median(d1, deriv_roll_k)

  # second derivative
  d2 <- c(NA_real_, diff(d1_sm))
  d2_sm <- roll_median(d2, deriv_roll_k)

  # smooth support quantities
  nb_gap_sm <- roll_median(df$nb_gap, deriv_roll_k)
  amu_sm <- roll_median(df$log_alpha_mu, deriv_roll_k)

  # left baseline for support z scores
  base_end <- min(max(floor(left_baseline_fraction * n), left_edge_buffer + 100L), n - right_edge_buffer)
  base_idx <- seq.int(1L + left_edge_buffer, base_end)

  d1_obj <- robust_center_scale(d1_sm[base_idx])
  d2_obj <- robust_center_scale(d2_sm[base_idx])
  gap_obj <- robust_center_scale(nb_gap_sm[base_idx])
  amu_obj <- robust_center_scale(amu_sm[base_idx])

  d1_z <- rep(NA_real_, n)
  d2_z <- rep(NA_real_, n)
  gap_z <- rep(NA_real_, n)
  amu_z <- rep(NA_real_, n)

  d1_ok <- is.finite(d1_sm)
  d2_ok <- is.finite(d2_sm)
  gap_ok <- is.finite(nb_gap_sm)
  amu_ok <- is.finite(amu_sm)

  d1_z[d1_ok] <- (d1_sm[d1_ok] - d1_obj$center) / d1_obj$scale
  d2_z[d2_ok] <- (d2_sm[d2_ok] - d2_obj$center) / d2_obj$scale
  gap_z[gap_ok] <- (nb_gap_sm[gap_ok] - gap_obj$center) / gap_obj$scale
  amu_z[amu_ok] <- (amu_sm[amu_ok] - amu_obj$center) / amu_obj$scale

  # combined terminal-inflection score:
  # high slope + high positive curvature + NB2>NB1 support + alpha*mu support
  score_raw <- rowMeans(cbind(
    scale01(d1_z),
    scale01(d2_z),
    scale01(gap_z),
    scale01(amu_z)
  ), na.rm = TRUE)
  score_sm <- roll_median(score_raw, deriv_roll_k)

  search_start <- max(base_end + 1L, floor(search_fraction_min * n))
  search_end <- min(n - right_edge_buffer, floor(search_fraction_max * n))
  search_idx <- seq.int(search_start, search_end)

  score_cut <- as.numeric(stats::quantile(
    score_sm[search_idx][is.finite(score_sm[search_idx])],
    probs = score_quantile_within_search,
    na.rm = TRUE,
    names = FALSE
  ))

  cond <- rep(FALSE, n)
  cond[search_idx] <- is.finite(score_sm[search_idx]) & (score_sm[search_idx] >= score_cut)

  # choose strongest block, then prefer the one closest to right side if tied
  use <- cond[search_start:search_end]
  r <- rle(use)
  ends <- cumsum(r$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  block_ids <- which(r$values & r$lengths >= min_block_length)

  if (!length(block_ids)) {
    block_ids <- which(r$values & r$lengths >= max(50L, floor(min_block_length / 2)))
  }

  if (!length(block_ids)) {
    best_idx <- search_idx[which.max(score_sm[search_idx])]
    block_start <- best_idx
    block_end <- best_idx
    block_mean <- score_sm[best_idx]
  } else {
    block_tbl <- lapply(block_ids, function(j) {
      lo_local <- starts[j]
      hi_local <- ends[j]
      lo_global <- search_start + lo_local - 1L
      hi_global <- search_start + hi_local - 1L
      data.frame(
        start = lo_global,
        end = hi_global,
        mean_score = mean(score_sm[lo_global:hi_global], na.rm = TRUE),
        rightness = hi_global
      )
    }) %>% bind_rows()

    # strongest block, break ties by more rightward end
    block_tbl <- block_tbl %>%
      arrange(desc(mean_score), desc(rightness))

    block_start <- block_tbl$start[1]
    block_end <- block_tbl$end[1]
    block_mean <- block_tbl$mean_score[1]
  }

  df$var_fit <- var_fit
  df$d1_sm <- d1_sm
  df$d2_sm <- d2_sm
  df$nb_gap_sm <- nb_gap_sm
  df$amu_sm <- amu_sm
  df$d1_z <- d1_z
  df$d2_z <- d2_z
  df$gap_z <- gap_z
  df$amu_z <- amu_z
  df$score_raw <- score_raw
  df$score_sm <- score_sm
  df$score_cut <- score_cut
  df$in_block <- FALSE
  df$in_block[block_start:block_end] <- TRUE

  list(
    onset_index = block_start,
    onset_rank = df$rank[block_start],
    block_end_index = block_end,
    block_end_rank = df$rank[block_end],
    mode = "left edge of strongest sustained terminal-inflection block",
    baseline_end = base_end,
    search_start = search_start,
    search_end = search_end,
    score_cut = score_cut,
    block_mean_score = block_mean,
    curve_df = df
  )
}

build_dataset_panel <- function(rank_df, onset_info, title_prefix, out_file) {
  df <- onset_info$curve_df
  cutoff_x <- onset_info$onset_rank
  block_end_x <- onset_info$block_end_rank
  label_text <- paste0(
    onset_info$mode,
    "\nStart rank = ", cutoff_x,
    "\nEnd rank = ", block_end_x
  )

  p1 <- ggplot(df, aes(rank, abs_loading)) +
    geom_line(linewidth = 0.8) +
    annotate("rect",
      xmin = cutoff_x,
      xmax = block_end_x,
      ymin = -Inf, ymax = Inf, alpha = 0.10
    ) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    annotate("label",
      x = cutoff_x,
      y = max(df$abs_loading, na.rm = TRUE),
      label = label_text,
      hjust = 1, vjust = 1, size = 3
    ) +
    labs(
      title = paste0(title_prefix, ": absolute PC1 loading series"),
      subtitle = "Left = lowest loading | Right = highest loading (leading edge)",
      x = "EVS rank",
      y = "|PC1 loading|"
    ) +
    theme_bw(base_size = 10)

  p2 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = var_fit, color = "Loess fit: log(1 + variance)"), linewidth = 1.0) +
    annotate("rect",
      xmin = cutoff_x,
      xmax = block_end_x,
      ymin = -Inf, ymax = Inf, alpha = 0.10
    ) +
    geom_point(
      data = df[df$rank == cutoff_x, , drop = FALSE],
      aes(x = rank, y = var_fit),
      size = 2.5
    ) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    labs(
      title = paste0(title_prefix, ": smoothed empirical variance curve"),
      x = "EVS rank",
      y = "Fitted log(1 + variance)"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p3 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = log_nb1, color = "NB1 = mu"), linewidth = 0.6, alpha = 0.35) +
    geom_line(aes(y = log_nb2, color = "NB2 = variance - mu"), linewidth = 0.6, alpha = 0.35) +
    geom_line(aes(y = nb_gap_sm, color = "Smoothed NB2 - NB1 gap"), linewidth = 1.0) +
    annotate("rect",
      xmin = cutoff_x,
      xmax = block_end_x,
      ymin = -Inf, ymax = Inf, alpha = 0.10
    ) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    labs(
      title = paste0(title_prefix, ": NB1 / NB2 support"),
      x = "EVS rank",
      y = "Log value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p4 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = d1_z, color = "Slope z"), linewidth = 1.0) +
    geom_line(aes(y = d2_z, color = "Curvature z"), linewidth = 1.0) +
    geom_line(aes(y = gap_z, color = "NB2-NB1 gap z"), linewidth = 1.0) +
    geom_line(aes(y = amu_z, color = "alpha*mu z"), linewidth = 1.0) +
    annotate("rect",
      xmin = cutoff_x,
      xmax = block_end_x,
      ymin = -Inf, ymax = Inf, alpha = 0.10
    ) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    labs(
      title = paste0(title_prefix, ": standardized terminal-inflection support"),
      x = "EVS rank",
      y = "Z score from left baseline"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p5 <- ggplot(df, aes(rank, score_sm)) +
    geom_line(linewidth = 1.0) +
    geom_hline(yintercept = onset_info$score_cut, linetype = 3) +
    annotate("rect",
      xmin = cutoff_x,
      xmax = block_end_x,
      ymin = -Inf, ymax = Inf, alpha = 0.10
    ) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    labs(
      title = paste0(title_prefix, ": terminal-inflection regime score"),
      subtitle = "Cutoff = left edge of strongest sustained terminal right-side block",
      x = "EVS rank",
      y = "Smoothed regime score"
    ) +
    theme_bw(base_size = 10)

  png(out_file, width = 2200, height = 3000, res = 200)
  grid.arrange(p1, p2, p3, p4, p5, ncol = 1)
  dev.off()
}

build_range_panel <- function(ctrl_onset, trt_onset, ctrl_df, trt_df, title_prefix, out_file) {
  ctrl_plot_df <- ctrl_onset$curve_df
  trt_plot_df  <- trt_onset$curve_df

  p <- ggplot() +
    geom_line(data = ctrl_plot_df, aes(rank, score_sm, color = "Control regime score"), linewidth = 1.0) +
    geom_line(data = trt_plot_df, aes(rank, score_sm, color = "Treatment regime score"), linewidth = 1.0) +
    annotate(
      "rect",
      xmin = min(ctrl_onset$onset_rank, trt_onset$onset_rank),
      xmax = max(ctrl_onset$onset_rank, trt_onset$onset_rank),
      ymin = -Inf, ymax = Inf, alpha = 0.08
    ) +
    geom_vline(xintercept = ctrl_onset$onset_rank, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = trt_onset$onset_rank, linetype = 3, linewidth = 0.8) +
    labs(
      title = paste0(title_prefix, ": treatment/control regime-change range"),
      subtitle = paste0(
        "Control rank = ", ctrl_onset$onset_rank,
        " | Treatment rank = ", trt_onset$onset_rank
      ),
      x = "EVS rank",
      y = "Terminal-inflection regime score"
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

  utils::write.csv(full_tbl, file.path(cmp_dir, paste0(cmp_name, "_feature_level_metrics.csv")), row.names = FALSE)

  ctrl_rank_df <- build_rank_series(full_tbl, "control")
  trt_rank_df  <- build_rank_series(full_tbl, "treatment")

  utils::write.csv(ctrl_rank_df, file.path(cmp_dir, paste0(cmp_name, "_control_rank_series.csv")), row.names = FALSE)
  utils::write.csv(trt_rank_df, file.path(cmp_dir, paste0(cmp_name, "_treatment_rank_series.csv")), row.names = FALSE)

  message("Selecting control terminal-inflection block for ", cmp_name); flush.console()
  ctrl_onset <- select_terminal_inflection_block(
    ctrl_rank_df,
    loess_span = loess_span,
    deriv_roll_k = deriv_roll_k,
    left_baseline_fraction = left_baseline_fraction,
    search_fraction_min = search_fraction_min,
    search_fraction_max = search_fraction_max,
    left_edge_buffer = left_edge_buffer,
    right_edge_buffer = right_edge_buffer,
    min_block_length = min_block_length,
    score_quantile_within_search = score_quantile_within_search
  )
  message("Control onset rank: ", ctrl_onset$onset_rank); flush.console()

  message("Selecting treatment terminal-inflection block for ", cmp_name); flush.console()
  trt_onset <- select_terminal_inflection_block(
    trt_rank_df,
    loess_span = loess_span,
    deriv_roll_k = deriv_roll_k,
    left_baseline_fraction = left_baseline_fraction,
    search_fraction_min = search_fraction_min,
    search_fraction_max = search_fraction_max,
    left_edge_buffer = left_edge_buffer,
    right_edge_buffer = right_edge_buffer,
    min_block_length = min_block_length,
    score_quantile_within_search = score_quantile_within_search
  )
  message("Treatment onset rank: ", trt_onset$onset_rank); flush.console()

  build_dataset_panel(
    ctrl_rank_df, ctrl_onset, paste0(cmp_name, " control"),
    file.path(cmp_dir, paste0(cmp_name, "_control_rank_panel.png"))
  )
  build_dataset_panel(
    trt_rank_df, trt_onset, paste0(cmp_name, " treatment"),
    file.path(cmp_dir, paste0(cmp_name, "_treatment_rank_panel.png"))
  )
  build_range_panel(
    ctrl_onset, trt_onset, ctrl_rank_df, trt_rank_df, cmp_name,
    file.path(cmp_dir, paste0(cmp_name, "_cutoff_range_panel.png"))
  )

  cutoff_summary <- tibble(
    comparison = cmp_name,
    control_cutoff_mode = ctrl_onset$mode,
    control_cutoff_rank = ctrl_onset$onset_rank,
    control_cutoff_fraction = ctrl_onset$onset_rank / nrow(ctrl_rank_df),
    control_block_end_rank = ctrl_onset$block_end_rank,
    treatment_cutoff_mode = trt_onset$mode,
    treatment_cutoff_rank = trt_onset$onset_rank,
    treatment_cutoff_fraction = trt_onset$onset_rank / nrow(trt_rank_df),
    treatment_block_end_rank = trt_onset$block_end_rank,
    cutoff_range_rank_min = min(ctrl_onset$onset_rank, trt_onset$onset_rank),
    cutoff_range_rank_max = max(ctrl_onset$onset_rank, trt_onset$onset_rank),
    cutoff_range_fraction_min = min(ctrl_onset$onset_rank / nrow(ctrl_rank_df), trt_onset$onset_rank / nrow(trt_rank_df)),
    cutoff_range_fraction_max = max(ctrl_onset$onset_rank / nrow(ctrl_rank_df), trt_onset$onset_rank / nrow(trt_rank_df))
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
