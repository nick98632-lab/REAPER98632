# =============================================================================
# SEQUENCE STAGE 1: EMPIRICAL NB REGIME SHIFT
# -----------------------------------------------------------------------------
# Keep the same figure structure as the alpha*mu version, but replace DESeq2
# gene-wise dispersion with empirical dispersion from raw counts:
#
#   Var = mu + alpha * mu^2
#   alpha_emp = (Var - mu) / mu^2
#   alpha_emp * mu = (Var - mu) / mu = IOD - 1
#
# This version uses only:
#   - raw-count mean
#   - raw-count variance
#   - empirical dispersion derived from them
#
# Output:
#   exports/empirical_nb_regime_shift/<comparison>_cutoff_folder/
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
out_root <- file.path(repo_dir, "exports", "empirical_nb_regime_shift")
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
divergence_k <- 75L
leading_edge_fraction <- 0.60
left_edge_exclusion_fraction <- 0.02

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
      rank_ctrl = rank(-ctrl_abs_loading, ties.method = "first"),
      rank_trt  = rank(-trt_abs_loading, ties.method = "first")
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
      alpha_mu_ctrl  = ifelse(mu_ctrl > 0, (var_ctrl - mu_ctrl) / mu_ctrl, NA_real_),
      alpha_mu_trt   = ifelse(mu_trt  > 0, (var_trt  - mu_trt)  / mu_trt,  NA_real_),
      log_alpha_mu_ctrl = ifelse(is.finite(alpha_mu_ctrl) & alpha_mu_ctrl > 0, log(alpha_mu_ctrl), NA_real_),
      log_alpha_mu_trt  = ifelse(is.finite(alpha_mu_trt)  & alpha_mu_trt  > 0, log(alpha_mu_trt),  NA_real_)
    )
}

build_rank_series <- function(full_tbl, arm = c("control", "treatment")) {
  arm <- match.arg(arm)

  if (arm == "control") {
    out <- full_tbl %>%
      transmute(
        feature_id,
        gene_symbol,
        rank = rank_ctrl,
        abs_loading = ctrl_abs_loading,
        mu = mu_ctrl,
        variance = var_ctrl,
        alpha = alpha_emp_ctrl,
        alpha_mu = alpha_mu_ctrl,
        log_alpha_mu = log_alpha_mu_ctrl
      ) %>%
      arrange(rank)
  } else {
    out <- full_tbl %>%
      transmute(
        feature_id,
        gene_symbol,
        rank = rank_trt,
        abs_loading = trt_abs_loading,
        mu = mu_trt,
        variance = var_trt,
        alpha = alpha_emp_trt,
        alpha_mu = alpha_mu_trt,
        log_alpha_mu = log_alpha_mu_trt
      ) %>%
      arrange(rank)
  }

  out %>%
    mutate(
      log_mu = ifelse(is.finite(mu) & mu > 0, log(mu), NA_real_),
      log_alpha = ifelse(is.finite(alpha) & alpha > 0, log(alpha), NA_real_)
    )
}

select_leading_edge_cutoff <- function(rank_df, k = 101L, divergence_k = 75L,
                                       leading_edge_fraction = 0.60,
                                       left_edge_exclusion_fraction = 0.02) {
  df <- rank_df %>%
    mutate(
      log_alpha_mu_sm = roll_median(log_alpha_mu, k),
      dist_to_zero = abs(log_alpha_mu_sm)
    )

  n <- nrow(df)
  search_start <- max(2L, floor(n * left_edge_exclusion_fraction))
  search_end <- max(search_start + 2L, min(n - divergence_k, floor(n * leading_edge_fraction)))
  candidates <- integer(0)

  if (search_end >= (search_start + 1L)) {
    for (i in search_start:search_end) {
      local_lo <- max(search_start, i - 1L)
      local_hi <- min(search_end, i + 1L)
      local_vals <- df$dist_to_zero[local_lo:local_hi]
      if (!all(is.finite(local_vals)) || !is.finite(df$dist_to_zero[i])) next

      is_local_min <- df$dist_to_zero[i] <= min(local_vals, na.rm = TRUE)
      if (!is_local_min) next

      f_hi <- min(n, i + divergence_k)
      if (f_hi <= i + 1L) next

      reopen_trend <- mean(diff(df$dist_to_zero[i:f_hi]), na.rm = TRUE)
      if (is.finite(reopen_trend) && reopen_trend > 0) {
        candidates <- c(candidates, i)
      }
    }
  }

  if (length(candidates) > 0L) {
    idx <- candidates[1L]
    mode <- "minimum distance to empirical alpha*mu = 1, then reopening"
  } else {
    valid <- which(is.finite(df$dist_to_zero))
    valid <- valid[valid >= search_start & valid <= search_end]
    if (!length(valid)) {
      idx <- search_start
    } else {
      idx <- valid[which.min(df$dist_to_zero[valid])]
    }
    mode <- "minimum distance fallback"
  }

  list(
    mode = mode,
    cutoff_rank = df$rank[idx],
    cutoff_index = idx,
    cutoff_fraction = idx / n,
    curve_df = df
  )
}

build_dataset_panel <- function(rank_df, cutoff_info, title_prefix, out_file) {
  df <- cutoff_info$curve_df
  cutoff_x <- cutoff_info$cutoff_rank
  label_text <- paste0(
    cutoff_info$mode,
    "\nRank = ", cutoff_info$cutoff_rank,
    "\nFraction = ", sprintf("%.3f", cutoff_info$cutoff_fraction)
  )

  p1 <- ggplot(df, aes(rank, abs_loading)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    annotate("label", x = cutoff_x, y = max(df$abs_loading, na.rm = TRUE), label = label_text, hjust = 0, vjust = 1, size = 3) +
    labs(
      title = paste0(title_prefix, ": absolute PC1 loading series"),
      x = "EVS rank",
      y = "|PC1 loading|"
    ) +
    theme_bw(base_size = 10)

  p2 <- ggplot(df, aes(rank)) +
    geom_line(aes(y = roll_median(log_mu, 101L), color = "Smoothed log(raw-count mean)"), linewidth = 1.0) +
    geom_line(aes(y = roll_median(log_alpha, 101L), color = "Smoothed log(empirical dispersion)"), linewidth = 1.0) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    labs(
      title = paste0(title_prefix, ": raw-count mean and empirical dispersion"),
      x = "EVS rank",
      y = "Log value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p3 <- ggplot(df, aes(rank)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(aes(y = log_alpha_mu, color = "Raw log(empirical alpha*mu)"), linewidth = 0.8, alpha = 0.30) +
    geom_line(aes(y = log_alpha_mu_sm, color = "Smoothed log(empirical alpha*mu)"), linewidth = 1.0) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    labs(
      title = paste0(title_prefix, ": empirical NB balance curve log(alpha*mu)"),
      x = "EVS rank",
      y = "log(empirical alpha*mu)"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p4 <- ggplot(df, aes(rank)) +
    geom_ribbon(aes(ymin = 0, ymax = dist_to_zero), alpha = 0.20) +
    geom_line(aes(y = dist_to_zero, color = "Distance to log(empirical alpha*mu)=0"), linewidth = 1.0) +
    geom_vline(xintercept = cutoff_x, linetype = 2, linewidth = 0.8) +
    labs(
      title = paste0(title_prefix, ": shaded distance-to-balance and selected cutoff"),
      x = "EVS rank",
      y = "|log(empirical alpha*mu)|"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  png(out_file, width = 2200, height = 2400, res = 200)
  grid.arrange(p1, p2, p3, p4, ncol = 1)
  dev.off()
}

build_range_panel <- function(ctrl_cutoff, trt_cutoff, ctrl_df, trt_df, title_prefix, out_file) {
  p <- ggplot() +
    geom_line(data = ctrl_df, aes(rank, dist_to_zero, color = "Control distance to balance"), linewidth = 1.0) +
    geom_line(data = trt_df, aes(rank, dist_to_zero, color = "Treatment distance to balance"), linewidth = 1.0) +
    annotate(
      "rect",
      xmin = min(ctrl_cutoff$cutoff_rank, trt_cutoff$cutoff_rank),
      xmax = max(ctrl_cutoff$cutoff_rank, trt_cutoff$cutoff_rank),
      ymin = -Inf, ymax = Inf, alpha = 0.08
    ) +
    geom_vline(xintercept = ctrl_cutoff$cutoff_rank, linetype = 2, linewidth = 0.8) +
    geom_vline(xintercept = trt_cutoff$cutoff_rank, linetype = 3, linewidth = 0.8) +
    labs(
      title = paste0(title_prefix, ": treatment/control cutoff range"),
      subtitle = paste0(
        "Control rank = ", ctrl_cutoff$cutoff_rank,
        " | Treatment rank = ", trt_cutoff$cutoff_rank
      ),
      x = "EVS rank",
      y = "|log(empirical alpha*mu)|"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  png(out_file, width = 2200, height = 1200, res = 200)
  print(p)
  dev.off()
}

# =============================================================================
# DATA IMPORT
# =============================================================================

count_file <- resolve_counts_file(count_file)
message("Reading raw count matrix from: ", count_file)

loaded <- read_count_matrix(count_file, meta_all$id)
count_mat <- loaded$count_matrix
annot_df <- loaded$annot_df

# =============================================================================
# MAIN LOOP
# =============================================================================

overall_rows <- list()

for (i in seq_len(nrow(comparison_table))) {
  comparison_row <- comparison_table[i, , drop = FALSE]
  cmp_name <- comparison_row$comparison_name[[1]]
  message("Processing ", cmp_name, "...")

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

  ctrl_cutoff <- select_leading_edge_cutoff(
    ctrl_rank_df, smoother_k, divergence_k,
    leading_edge_fraction, left_edge_exclusion_fraction
  )
  trt_cutoff <- select_leading_edge_cutoff(
    trt_rank_df, smoother_k, divergence_k,
    leading_edge_fraction, left_edge_exclusion_fraction
  )

  build_dataset_panel(
    ctrl_rank_df, ctrl_cutoff,
    paste0(cmp_name, " control"),
    file.path(cmp_dir, paste0(cmp_name, "_control_rank_panel.png"))
  )

  build_dataset_panel(
    trt_rank_df, trt_cutoff,
    paste0(cmp_name, " treatment"),
    file.path(cmp_dir, paste0(cmp_name, "_treatment_rank_panel.png"))
  )

  build_range_panel(
    ctrl_cutoff, trt_cutoff,
    ctrl_cutoff$curve_df, trt_cutoff$curve_df,
    cmp_name,
    file.path(cmp_dir, paste0(cmp_name, "_cutoff_range_panel.png"))
  )

  cutoff_summary <- tibble(
    comparison = cmp_name,
    control_cutoff_mode = ctrl_cutoff$mode,
    control_cutoff_rank = ctrl_cutoff$cutoff_rank,
    control_cutoff_fraction = ctrl_cutoff$cutoff_fraction,
    treatment_cutoff_mode = trt_cutoff$mode,
    treatment_cutoff_rank = trt_cutoff$cutoff_rank,
    treatment_cutoff_fraction = trt_cutoff$cutoff_fraction,
    cutoff_range_rank_min = min(ctrl_cutoff$cutoff_rank, trt_cutoff$cutoff_rank),
    cutoff_range_rank_max = max(ctrl_cutoff$cutoff_rank, trt_cutoff$cutoff_rank),
    cutoff_range_fraction_min = min(ctrl_cutoff$cutoff_fraction, trt_cutoff$cutoff_fraction),
    cutoff_range_fraction_max = max(ctrl_cutoff$cutoff_fraction, trt_cutoff$cutoff_fraction)
  )

  utils::write.csv(
    cutoff_summary,
    file.path(cmp_dir, paste0(cmp_name, "_cutoff_summary.csv")),
    row.names = FALSE
  )
  overall_rows[[cmp_name]] <- cutoff_summary
}

overall_summary <- bind_rows(overall_rows)
utils::write.csv(overall_summary, file.path(out_root, "overall_cutoff_summary.csv"), row.names = FALSE)

message("Done. Outputs written to: ", out_root)
