# =============================================================================
# SEQUENCE STAGE 1: VARIANCE + NB1/NB2 REGIME SPLIT (FAST)
# -----------------------------------------------------------------------------
# Primary driver:
#   empirical variance across EVS-ranked absolute PC1 loading series
#
# Support panels:
#   NB1 = mu
#   NB2 = max(Var - mu, 0)
#   alpha*mu = (Var - mu) / mu = NB2 / NB1   (when mu > 0)
#
# Regime interpretation:
#   - remainder is more NB1-like (linear term)
#   - leading edge is more NB2-like (quadratic/overdispersion term)
#   - the breakpoint is chosen on the variance curve, then shown across all panels
#
# Key design choices:
#   - full ranked series is plotted and demonstrated
#   - only tiny edge zones are excluded from breakpoint search to avoid artifacts
#   - breakpoint fitting is FAST via coarse-to-fine search on a smoothed curve
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

repo_dir <- getwd()
input_dir <- file.path(repo_dir, "data")
count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
out_root <- file.path(repo_dir, "exports", "variance_NB1_NB2_regime_split_fast")
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

smoother_k <- 101L
left_edge_exclusion_fraction <- 0.02
right_edge_exclusion_fraction <- 0.02
min_segment_size <- 200L
coarse_grid_step <- 100L
fine_window <- 400L

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

  annot_df <- data.frame(feature_id = as.character(raw_df[[feature_col]]), stringsAsFactors = FALSE)
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
    stop(paste0("Missing samples for comparison ", comparison_row$comparison_name, ": ", paste(missing_ids, collapse = ", ")), call. = FALSE)
  }
  list(count_matrix = count_matrix[, keep_ids, drop = FALSE], trt_ids = trt_ids, ctrl_ids = ctrl_ids)
}

safe_mean <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_var <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::var(x)
}

roll_median <- function(x, k = 101L) {
  x <- as.numeric(x)
  n <- length(x)
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
    trt_abs_loading = trt_load
  ) %>% mutate(
    rank_ctrl = rank(-ctrl_abs_loading, ties.method = "first"),
    rank_trt  = rank(-trt_abs_loading, ties.method = "first")
  )
}

compute_feature_metrics <- function(count_mat, ctrl_cols, trt_cols, feature_ids) {
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
  ) %>% mutate(
    nb1_ctrl = mu_ctrl,
    nb1_trt = mu_trt,
    nb2_ctrl = pmax(var_ctrl - mu_ctrl, 0),
    nb2_trt = pmax(var_trt - mu_trt, 0),
    alpha_mu_ctrl = ifelse(mu_ctrl > 0, (var_ctrl - mu_ctrl) / mu_ctrl, NA_real_),
    alpha_mu_trt  = ifelse(mu_trt > 0, (var_trt - mu_trt) / mu_trt, NA_real_),
    log_var_ctrl = ifelse(is.finite(var_ctrl) & var_ctrl >= 0, log1p(var_ctrl), NA_real_),
    log_var_trt  = ifelse(is.finite(var_trt) & var_trt >= 0, log1p(var_trt), NA_real_),
    log_nb1_ctrl = ifelse(is.finite(nb1_ctrl) & nb1_ctrl >= 0, log1p(nb1_ctrl), NA_real_),
    log_nb1_trt  = ifelse(is.finite(nb1_trt) & nb1_trt >= 0, log1p(nb1_trt), NA_real_),
    log_nb2_ctrl = ifelse(is.finite(nb2_ctrl) & nb2_ctrl >= 0, log1p(nb2_ctrl), NA_real_),
    log_nb2_trt  = ifelse(is.finite(nb2_trt) & nb2_trt >= 0, log1p(nb2_trt), NA_real_),
    log_alpha_mu_ctrl = ifelse(is.finite(alpha_mu_ctrl) & alpha_mu_ctrl > 0, log(alpha_mu_ctrl), NA_real_),
    log_alpha_mu_trt  = ifelse(is.finite(alpha_mu_trt) & alpha_mu_trt > 0, log(alpha_mu_trt), NA_real_)
  )
}

build_rank_series <- function(full_tbl, arm = c("control", "treatment")) {
  arm <- match.arg(arm)
  if (arm == "control") {
    full_tbl %>% transmute(
      feature_id, gene_symbol, rank = rank_ctrl, abs_loading = ctrl_abs_loading,
      variance = var_ctrl, nb1 = nb1_ctrl, nb2 = nb2_ctrl,
      log_variance = log_var_ctrl, log_nb1 = log_nb1_ctrl, log_nb2 = log_nb2_ctrl,
      log_alpha_mu = log_alpha_mu_ctrl
    ) %>% arrange(rank)
  } else {
    full_tbl %>% transmute(
      feature_id, gene_symbol, rank = rank_trt, abs_loading = trt_abs_loading,
      variance = var_trt, nb1 = nb1_trt, nb2 = nb2_trt,
      log_variance = log_var_trt, log_nb1 = log_nb1_trt, log_nb2 = log_nb2_trt,
      log_alpha_mu = log_alpha_mu_trt
    ) %>% arrange(rank)
  }
}

fit_piecewise_breakpoint_fast <- function(x, y,
                                          left_exclusion_fraction = 0.02,
                                          right_exclusion_fraction = 0.02,
                                          min_segment_size = 200L,
                                          coarse_grid_step = 100L,
                                          fine_window = 400L) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  n <- length(y)

  fitted_full <- rep(NA_real_, n)
  score <- rep(NA_real_, n)

  first_ok <- which(ok)[1]
  last_ok <- tail(which(ok), 1)
  search_start <- max(first_ok, floor(n * left_exclusion_fraction), min_segment_size)
  search_end <- min(last_ok, floor(n * (1 - right_exclusion_fraction)), n - min_segment_size)

  if (!is.finite(search_start) || !is.finite(search_end) || search_end <= search_start) {
    idx <- max(1L, min(n, floor(n * 0.10)))
    return(list(breakpoint_index = idx, breakpoint_rank = x[idx], mode = "fallback invalid search", fitted = fitted_full, rss = NA_real_, search_start = search_start, search_end = search_end, score = score))
  }

  coarse_candidates <- seq(search_start, search_end, by = coarse_grid_step)
  coarse_candidates <- coarse_candidates[coarse_candidates > min_segment_size & coarse_candidates < (n - min_segment_size)]
  if (length(coarse_candidates) < 3L) coarse_candidates <- seq(search_start, search_end)

  eval_candidate <- function(c) {
    z <- pmax(0, x - x[c])
    fit_df <- data.frame(y = y, x = x, z = z)
    fit_df <- fit_df[is.finite(fit_df$y) & is.finite(fit_df$x) & is.finite(fit_df$z), , drop = FALSE]
    if (nrow(fit_df) < (2L * min_segment_size + 5L)) return(list(rss = Inf, fit = NULL))
    fit <- lm(y ~ x + z, data = fit_df)
    list(rss = sum(resid(fit)^2), fit = fit)
  }

  best_rss <- Inf
  best_c <- coarse_candidates[1]
  best_fit <- NULL

  for (c in coarse_candidates) {
    res <- eval_candidate(c)
    score[c] <- res$rss
    if (is.finite(res$rss) && res$rss < best_rss) {
      best_rss <- res$rss
      best_c <- c
      best_fit <- res$fit
    }
  }

  fine_start <- max(search_start, best_c - fine_window)
  fine_end <- min(search_end, best_c + fine_window)
  fine_candidates <- seq(fine_start, fine_end)
  fine_candidates <- fine_candidates[fine_candidates > min_segment_size & fine_candidates < (n - min_segment_size)]

  for (c in fine_candidates) {
    res <- eval_candidate(c)
    score[c] <- res$rss
    if (is.finite(res$rss) && res$rss < best_rss) {
      best_rss <- res$rss
      best_c <- c
      best_fit <- res$fit
    }
  }

  if (!is.null(best_fit)) {
    fitted_full[ok] <- predict(best_fit, newdata = data.frame(x = x[ok], z = pmax(0, x[ok] - x[best_c])))
  }

  list(
    breakpoint_index = best_c,
    breakpoint_rank = x[best_c],
    mode = "fast global piecewise breakpoint",
    fitted = fitted_full,
    rss = best_rss,
    search_start = search_start,
    search_end = search_end,
    score = score
  )
}

build_dataset_panel <- function(rank_df, bp_info, title_prefix, out_file, smoother_k = 101L) {
  df <- rank_df %>% mutate(
    log_variance_sm = roll_median(log_variance, smoother_k),
    log_nb1_sm = roll_median(log_nb1, smoother_k),
    log_nb2_sm = roll_median(log_nb2, smoother_k),
    log_alpha_mu_sm = roll_median(log_alpha_mu, smoother_k),
    fitted_piecewise = bp_info$fitted
  )

  cutoff_x <- bp_info$breakpoint_rank
  label_text <- paste0(bp_info$mode, "\nRank = ", cutoff_x)

  p1 <- ggplot(df, aes(rank, abs_loading)) +
    geom_line(linewidth = 0.8) +
    annotate("rect", xmin = bp_info$search_start, xmax = bp_info$search_end, ymin = -Inf, ymax = Inf, alpha = 0.06) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    annotate("label", x = cutoff_x, y = max(df$abs_loading, na.rm = TRUE), label = label_text, hjust = 0, vjust = 1, size = 3) +
    labs(title = paste0(title_prefix, ": absolute PC1 loading series"), x = "EVS rank", y = "|PC1 loading|") +
    theme_bw(base_size = 10)

  p2 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = log_variance_sm, color = "Smoothed log empirical variance"), linewidth = 1.0) +
    geom_line(aes(y = fitted_piecewise, color = "Piecewise fit"), linewidth = 1.0, linetype = 2) +
    annotate("rect", xmin = bp_info$search_start, xmax = bp_info$search_end, ymin = -Inf, ymax = Inf, alpha = 0.06) +
    geom_point(data = df[df$rank == cutoff_x, , drop = FALSE], aes(x = rank, y = log_variance_sm), size = 2.5) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    labs(title = paste0(title_prefix, ": primary variance regime curve"), x = "EVS rank", y = "log(1 + variance)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p3 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = log_nb1_sm, color = "Smoothed log(NB1 = mu)"), linewidth = 1.0) +
    geom_line(aes(y = log_nb2_sm, color = "Smoothed log(NB2 = variance - mu)"), linewidth = 1.0) +
    annotate("rect", xmin = bp_info$search_start, xmax = bp_info$search_end, ymin = -Inf, ymax = Inf, alpha = 0.06) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    labs(title = paste0(title_prefix, ": NB1 vs NB2 support"), x = "EVS rank", y = "Log value") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p4 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = log_alpha_mu_sm, color = "Smoothed log(alpha*mu)"), linewidth = 1.0) +
    annotate("rect", xmin = bp_info$search_start, xmax = bp_info$search_end, ymin = -Inf, ymax = Inf, alpha = 0.06) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    labs(title = paste0(title_prefix, ": empirical NB balance support"), x = "EVS rank", y = "log(alpha*mu)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  png(out_file, width = 2200, height = 2400, res = 200)
  grid.arrange(p1, p2, p3, p4, ncol = 1)
  dev.off()
}

build_range_panel <- function(ctrl_bp, trt_bp, ctrl_df, trt_df, title_prefix, out_file, smoother_k = 101L) {
  ctrl_plot_df <- ctrl_df %>% mutate(log_variance_sm = roll_median(log_variance, smoother_k))
  trt_plot_df  <- trt_df  %>% mutate(log_variance_sm = roll_median(log_variance, smoother_k))

  p <- ggplot() +
    geom_line(data = ctrl_plot_df, aes(rank, log_variance_sm, color = "Control smoothed log variance"), linewidth = 1.0) +
    geom_line(data = trt_plot_df, aes(rank, log_variance_sm, color = "Treatment smoothed log variance"), linewidth = 1.0) +
    annotate("rect", xmin = min(ctrl_bp$breakpoint_rank, trt_bp$breakpoint_rank), xmax = max(ctrl_bp$breakpoint_rank, trt_bp$breakpoint_rank), ymin = -Inf, ymax = Inf, alpha = 0.08) +
    geom_vline(xintercept = ctrl_bp$breakpoint_rank, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = trt_bp$breakpoint_rank, linetype = 3, linewidth = 0.8) +
    labs(title = paste0(title_prefix, ": treatment/control regime-change range"), subtitle = paste0("Control rank = ", ctrl_bp$breakpoint_rank, " | Treatment rank = ", trt_bp$breakpoint_rank), x = "EVS rank", y = "Smoothed log variance") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  png(out_file, width = 2200, height = 1200, res = 200)
  print(p)
  dev.off()
}

message("SEQUENCE.R started")
message("Resolving count file...")
count_file <- resolve_counts_file(count_file)
message("Using count file: ", count_file)
message("Reading count matrix...")
loaded <- read_count_matrix(count_file, meta_all$id)
count_mat <- loaded$count_matrix
annot_df <- loaded$annot_df
message("Count matrix dimensions: ", nrow(count_mat), " features x ", ncol(count_mat), " samples")

overall_rows <- list()
for (i in seq_len(nrow(comparison_table))) {
  comparison_row <- comparison_table[i, , drop = FALSE]
  cmp_name <- comparison_row$comparison_name[[1]]
  message("Processing comparison: ", cmp_name)

  cmp_dir <- file.path(out_root, paste0(cmp_name, "_cutoff_folder"))
  dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

  comp <- subset_comparison(count_mat, comparison_row, meta_all)
  cmp_counts <- comp$count_matrix
  kept_ids <- rownames(cmp_counts)

  evs_tbl <- build_evs_table(cmp_counts, comp$ctrl_ids, comp$trt_ids, kept_ids)
  metrics_tbl <- compute_feature_metrics(cmp_counts, comp$ctrl_ids, comp$trt_ids, kept_ids)

  full_tbl <- evs_tbl %>% left_join(metrics_tbl, by = "feature_id") %>% left_join(annot_df, by = "feature_id")
  utils::write.csv(full_tbl, file.path(cmp_dir, paste0(cmp_name, "_feature_level_metrics.csv")), row.names = FALSE)

  ctrl_rank_df <- build_rank_series(full_tbl, "control")
  trt_rank_df <- build_rank_series(full_tbl, "treatment")
  utils::write.csv(ctrl_rank_df, file.path(cmp_dir, paste0(cmp_name, "_control_rank_series.csv")), row.names = FALSE)
  utils::write.csv(trt_rank_df, file.path(cmp_dir, paste0(cmp_name, "_treatment_rank_series.csv")), row.names = FALSE)

  message("Fitting control breakpoint for ", cmp_name)
  ctrl_bp <- fit_piecewise_breakpoint_fast(ctrl_rank_df$rank, roll_median(ctrl_rank_df$log_variance, smoother_k), left_edge_exclusion_fraction, right_edge_exclusion_fraction, min_segment_size, coarse_grid_step, fine_window)
  message("Control breakpoint rank: ", ctrl_bp$breakpoint_rank)

  message("Fitting treatment breakpoint for ", cmp_name)
  trt_bp <- fit_piecewise_breakpoint_fast(trt_rank_df$rank, roll_median(trt_rank_df$log_variance, smoother_k), left_edge_exclusion_fraction, right_edge_exclusion_fraction, min_segment_size, coarse_grid_step, fine_window)
  message("Treatment breakpoint rank: ", trt_bp$breakpoint_rank)

  build_dataset_panel(ctrl_rank_df, ctrl_bp, paste0(cmp_name, " control"), file.path(cmp_dir, paste0(cmp_name, "_control_rank_panel.png")), smoother_k)
  build_dataset_panel(trt_rank_df, trt_bp, paste0(cmp_name, " treatment"), file.path(cmp_dir, paste0(cmp_name, "_treatment_rank_panel.png")), smoother_k)
  build_range_panel(ctrl_bp, trt_bp, ctrl_rank_df, trt_rank_df, cmp_name, file.path(cmp_dir, paste0(cmp_name, "_cutoff_range_panel.png")), smoother_k)

  cutoff_summary <- tibble(
    comparison = cmp_name,
    control_cutoff_mode = ctrl_bp$mode,
    control_cutoff_rank = ctrl_bp$breakpoint_rank,
    control_cutoff_fraction = ctrl_bp$breakpoint_rank / nrow(ctrl_rank_df),
    treatment_cutoff_mode = trt_bp$mode,
    treatment_cutoff_rank = trt_bp$breakpoint_rank,
    treatment_cutoff_fraction = trt_bp$breakpoint_rank / nrow(trt_rank_df),
    cutoff_range_rank_min = min(ctrl_bp$breakpoint_rank, trt_bp$breakpoint_rank),
    cutoff_range_rank_max = max(ctrl_bp$breakpoint_rank, trt_bp$breakpoint_rank),
    cutoff_range_fraction_min = min(ctrl_bp$breakpoint_rank / nrow(ctrl_rank_df), trt_bp$breakpoint_rank / nrow(trt_rank_df)),
    cutoff_range_fraction_max = max(ctrl_bp$breakpoint_rank / nrow(ctrl_rank_df), trt_bp$breakpoint_rank / nrow(trt_rank_df))
  )
  utils::write.csv(cutoff_summary, file.path(cmp_dir, paste0(cmp_name, "_cutoff_summary.csv")), row.names = FALSE)
  overall_rows[[cmp_name]] <- cutoff_summary
}

overall_summary <- bind_rows(overall_rows)
utils::write.csv(overall_summary, file.path(out_root, "overall_cutoff_summary.csv"), row.names = FALSE)
message("Done. Outputs written to: ", out_root)
