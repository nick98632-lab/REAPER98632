
# =============================================================================
# SEQUENCE STAGE 1: EMPIRICAL + GENE-WISE NB LEADING-EDGE CUTOFF
# -----------------------------------------------------------------------------
# This script:
#   1. reads the WTTS raw count matrix
#   2. builds EVS ranks separately for treatment and control using PC1 absolute
#      loadings from log2(raw counts + 1)
#   3. computes empirical feature-level mean and variance from raw counts
#   4. computes empirical IOD and CV^2 from those raw-count moments
#   5. extracts DESeq2 gene-wise dispersion estimates only
#   6. computes NB-support IOD and CV^2 using gene-wise dispersion
#   7. summarizes treatment and control separately in percentile bins using the
#      median
#   8. defines leading-edge cutoffs from the leading-edge side, using the first
#      sustained empirical divergence region (IOD above CV^2), with a fallback
#      to the nearest-balance point if sustained divergence is absent
#   9. writes comparison-specific folders with tables and plots
#
# Main empirical quantities:
#   mu_i      = mean(raw counts)
#   var_i     = variance(raw counts)
#   IOD_i     = var_i / mu_i
#   CV2_i     = var_i / mu_i^2
#
# Gene-wise NB support quantities:
#   IOD_NB_i  = 1 + alpha_i * mu_i
#   CV2_NB_i  = 1 / mu_i + alpha_i
#
# Final reporting:
#   - treatment cutoff
#   - control cutoff
#   - cutoff range between treatment and control
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(SummarizedExperiment)
})

# =============================================================================
# USER SETTINGS
# =============================================================================

repo_dir <- getwd()
input_dir <- file.path(repo_dir, "data")
count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
out_root <- file.path(repo_dir, "exports", "leading_edge_cutoff_folders")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  treatment_prefix = c("R0", "R2", "R4", "R8"),
  control_prefix = c("ZT6", "ZT8", "ZT10", "ZT14"),
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

percentile_step <- 0.01
min_bin_n <- 5L
sustained_bins <- 3L

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
    coldata = meta_all[keep_ids, , drop = FALSE],
    trt_ids = trt_ids,
    ctrl_ids = ctrl_ids
  )
}

safe_var <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::var(x)
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
  stats::median(x)
}

rescale01 <- function(x) {
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

make_percentile_bins <- function(n, step = 0.01) {
  probs <- seq(step, 1, by = step)
  if (tail(probs, 1) < 1) probs <- c(probs, 1)
  centers <- unique(pmax(1L, pmin(n, round(probs * n))))
  dplyr::tibble(percentile = probs[seq_along(centers)], rank_index = centers)
}

compute_group_pc1_loadings <- function(count_mat, group_cols) {
  mat <- count_mat[, group_cols, drop = FALSE]
  mat <- log2(mat + 1)
  mat <- t(mat)
  if (nrow(mat) < 2L) return(rep(NA_real_, ncol(mat)))
  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  abs(pca$rotation[, 1L])
}

build_evs_table <- function(raw_counts_mat, ctrl_cols, trt_cols, feature_ids) {
  ctrl_load <- compute_group_pc1_loadings(raw_counts_mat, ctrl_cols)
  trt_load  <- compute_group_pc1_loadings(raw_counts_mat, trt_cols)

  dplyr::tibble(
    feature_id = feature_ids,
    ctrl_abs_loading = ctrl_load,
    trt_abs_loading  = trt_load
  ) %>%
    dplyr::mutate(
      rank_ctrl = rank(-ctrl_abs_loading, ties.method = "first"),
      rank_trt  = rank(-trt_abs_loading, ties.method = "first")
    )
}

extract_genewise_dispersion <- function(dds) {
  disp_df <- as.data.frame(SummarizedExperiment::mcols(dds))
  alpha <- if ("dispGeneEst" %in% names(disp_df)) disp_df$dispGeneEst else if ("dispersion" %in% names(disp_df)) disp_df$dispersion else rep(NA_real_, nrow(disp_df))
  dplyr::tibble(
    feature_id = rownames(disp_df),
    alpha_gene_wise = as.numeric(alpha)
  )
}

compute_feature_metrics <- function(raw_counts_mat, ctrl_cols, trt_cols, alpha_tbl) {
  ctrl_mu  <- apply(raw_counts_mat[, ctrl_cols, drop = FALSE], 1L, safe_mean)
  trt_mu   <- apply(raw_counts_mat[, trt_cols, drop = FALSE], 1L, safe_mean)
  ctrl_var <- apply(raw_counts_mat[, ctrl_cols, drop = FALSE], 1L, safe_var)
  trt_var  <- apply(raw_counts_mat[, trt_cols, drop = FALSE], 1L, safe_var)

  dplyr::tibble(
    feature_id = rownames(raw_counts_mat),
    mu_ctrl = ctrl_mu,
    mu_trt = trt_mu,
    variance_ctrl = ctrl_var,
    variance_trt = trt_var
  ) %>%
    dplyr::left_join(alpha_tbl, by = "feature_id") %>%
    dplyr::mutate(
      iod_emp_ctrl = ifelse(mu_ctrl > 0, variance_ctrl / mu_ctrl, NA_real_),
      iod_emp_trt  = ifelse(mu_trt  > 0, variance_trt / mu_trt, NA_real_),
      cv2_emp_ctrl = ifelse(mu_ctrl > 0, variance_ctrl / (mu_ctrl ^ 2), NA_real_),
      cv2_emp_trt  = ifelse(mu_trt  > 0, variance_trt / (mu_trt ^ 2), NA_real_),

      iod_nb_ctrl = ifelse(mu_ctrl > 0 & is.finite(alpha_gene_wise), 1 + alpha_gene_wise * mu_ctrl, NA_real_),
      iod_nb_trt  = ifelse(mu_trt  > 0 & is.finite(alpha_gene_wise), 1 + alpha_gene_wise * mu_trt,  NA_real_),
      cv2_nb_ctrl = ifelse(mu_ctrl > 0 & is.finite(alpha_gene_wise), (1 / mu_ctrl) + alpha_gene_wise, NA_real_),
      cv2_nb_trt  = ifelse(mu_trt  > 0 & is.finite(alpha_gene_wise), (1 / mu_trt) + alpha_gene_wise, NA_real_),

      diff_emp_ctrl = iod_emp_ctrl - cv2_emp_ctrl,
      diff_emp_trt  = iod_emp_trt - cv2_emp_trt,
      diff_nb_ctrl  = iod_nb_ctrl - cv2_nb_ctrl,
      diff_nb_trt   = iod_nb_trt - cv2_nb_trt,

      log_ratio_emp_ctrl = ifelse(iod_emp_ctrl > 0 & cv2_emp_ctrl > 0, log(iod_emp_ctrl / cv2_emp_ctrl), NA_real_),
      log_ratio_emp_trt  = ifelse(iod_emp_trt > 0 & cv2_emp_trt > 0, log(iod_emp_trt / cv2_emp_trt), NA_real_),
      log_ratio_nb_ctrl  = ifelse(iod_nb_ctrl > 0 & cv2_nb_ctrl > 0, log(iod_nb_ctrl / cv2_nb_ctrl), NA_real_),
      log_ratio_nb_trt   = ifelse(iod_nb_trt > 0 & cv2_nb_trt > 0, log(iod_nb_trt / cv2_nb_trt), NA_real_)
    )
}

summarize_by_percentile_median <- function(df, value_cols, percentile_step = 0.01, min_bin_n = 5L) {
  n <- nrow(df)
  bins <- make_percentile_bins(n, percentile_step)
  out <- vector("list", nrow(bins))

  for (i in seq_len(nrow(bins))) {
    lo <- if (i == 1L) 1L else bins$rank_index[i - 1L] + 1L
    hi <- bins$rank_index[i]
    chunk <- df[lo:hi, , drop = FALSE]

    row <- dplyr::tibble(
      percentile = bins$percentile[i],
      rank_index = bins$rank_index[i],
      lo_rank = lo,
      hi_rank = hi,
      n_bin = nrow(chunk)
    )

    for (nm in value_cols) {
      vals <- chunk[[nm]]
      vals <- vals[is.finite(vals)]
      row[[nm]] <- if (length(vals) >= min_bin_n) stats::median(vals) else NA_real_
    }

    out[[i]] <- row
  }

  dplyr::bind_rows(out)
}

select_leading_edge_cutoff <- function(sum_df, diff_col = "diff_emp", iod_col = "iod_emp", cv2_col = "cv2_emp", sustained_bins = 3L) {
  x <- sum_df$percentile
  d <- sum_df[[diff_col]]
  iod <- sum_df[[iod_col]]
  cv2 <- sum_df[[cv2_col]]
  n <- nrow(sum_df)

  if (n == 0L) {
    return(list(
      mode = "none",
      cutoff_percentile = NA_real_,
      cutoff_rank_index = NA_real_,
      note = "No bins available"
    ))
  }

  # first sustained positive region from the leading-edge side
  if (n >= sustained_bins) {
    for (i in seq_len(n - sustained_bins + 1L)) {
      idx <- i:(i + sustained_bins - 1L)
      ok <- is.finite(d[idx])
      if (all(ok) && all(d[idx] > 0)) {
        return(list(
          mode = "first_sustained_positive_gap",
          cutoff_percentile = x[i],
          cutoff_rank_index = sum_df$rank_index[i],
          note = paste0("First sustained empirical IOD>CV² gap across ", sustained_bins, " bins")
        ))
      }
    }
  }

  # fallback: nearest empirical balance from the leading-edge side
  valid <- which(is.finite(d))
  if (length(valid) > 0L) {
    i <- valid[which.min(abs(d[valid]))]
    return(list(
      mode = "nearest_balance_fallback",
      cutoff_percentile = x[i],
      cutoff_rank_index = sum_df$rank_index[i],
      note = "Nearest empirical balance point"
    ))
  }

  list(
    mode = "none",
    cutoff_percentile = NA_real_,
    cutoff_rank_index = NA_real_,
    note = "No valid empirical difference bins"
  )
}

build_dataset_summary <- function(full_tbl, dataset = c("control", "treatment"), percentile_step = 0.01, min_bin_n = 5L) {
  dataset <- match.arg(dataset)

  if (dataset == "control") {
    rank_tbl <- full_tbl %>%
      dplyr::arrange(rank_ctrl) %>%
      dplyr::transmute(
        feature_id,
        rank = rank_ctrl,
        mu = mu_ctrl,
        variance = variance_ctrl,
        alpha_gene_wise = alpha_gene_wise,
        iod_emp = iod_emp_ctrl,
        cv2_emp = cv2_emp_ctrl,
        iod_nb = iod_nb_ctrl,
        cv2_nb = cv2_nb_ctrl,
        diff_emp = diff_emp_ctrl,
        diff_nb = diff_nb_ctrl,
        log_ratio_emp = log_ratio_emp_ctrl,
        log_ratio_nb = log_ratio_nb_ctrl
      )
  } else {
    rank_tbl <- full_tbl %>%
      dplyr::arrange(rank_trt) %>%
      dplyr::transmute(
        feature_id,
        rank = rank_trt,
        mu = mu_trt,
        variance = variance_trt,
        alpha_gene_wise = alpha_gene_wise,
        iod_emp = iod_emp_trt,
        cv2_emp = cv2_emp_trt,
        iod_nb = iod_nb_trt,
        cv2_nb = cv2_nb_trt,
        diff_emp = diff_emp_trt,
        diff_nb = diff_nb_trt,
        log_ratio_emp = log_ratio_emp_trt,
        log_ratio_nb = log_ratio_nb_trt
      )
  }

  value_cols <- c(
    "mu", "variance", "alpha_gene_wise",
    "iod_emp", "cv2_emp", "iod_nb", "cv2_nb",
    "diff_emp", "diff_nb",
    "log_ratio_emp", "log_ratio_nb"
  )

  sum_df <- summarize_by_percentile_median(rank_tbl, value_cols, percentile_step, min_bin_n) %>%
    dplyr::mutate(
      iod_emp_scaled = rescale01(iod_emp),
      cv2_emp_scaled = rescale01(cv2_emp),
      iod_nb_scaled = rescale01(iod_nb),
      cv2_nb_scaled = rescale01(cv2_nb)
    )

  list(rank_tbl = rank_tbl, sum_df = sum_df)
}

build_dataset_panel <- function(sum_df, title_prefix, cutoff_info, out_file) {
  cutoff_label <- if (is.finite(cutoff_info$cutoff_percentile)) {
    paste0(
      "Cutoff: ", cutoff_info$mode,
      "\nPercentile = ", sprintf("%.3f", cutoff_info$cutoff_percentile),
      "\nRank = ", cutoff_info$cutoff_rank_index
    )
  } else {
    paste0("Cutoff unavailable\n", cutoff_info$note)
  }

  p1_df <- dplyr::bind_rows(
    dplyr::transmute(sum_df, percentile, value = mu, metric = "Median raw-count mean"),
    dplyr::transmute(sum_df, percentile, value = variance, metric = "Median raw-count variance"),
    dplyr::transmute(sum_df, percentile, value = alpha_gene_wise, metric = "Median gene-wise dispersion")
  )

  p2_df <- dplyr::bind_rows(
    dplyr::transmute(sum_df, percentile, value = iod_emp, metric = "Empirical IOD"),
    dplyr::transmute(sum_df, percentile, value = cv2_emp, metric = "Empirical CV²"),
    dplyr::transmute(sum_df, percentile, value = iod_nb, metric = "Gene-wise NB IOD"),
    dplyr::transmute(sum_df, percentile, value = cv2_nb, metric = "Gene-wise NB CV²")
  )

  p3_df <- dplyr::bind_rows(
    dplyr::transmute(sum_df, percentile, value = diff_emp, metric = "Empirical difference"),
    dplyr::transmute(sum_df, percentile, value = diff_nb, metric = "Gene-wise NB difference"),
    dplyr::transmute(sum_df, percentile, value = log_ratio_emp, metric = "Empirical log-ratio"),
    dplyr::transmute(sum_df, percentile, value = log_ratio_nb, metric = "Gene-wise NB log-ratio")
  )

  p1 <- ggplot(p1_df, aes(percentile, value, color = metric)) +
    geom_line(linewidth = 0.9) +
    labs(
      title = paste0(title_prefix, ": median raw-count mean, variance, and gene-wise dispersion"),
      x = "EVS percentile",
      y = "Median bin value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")

  p2 <- ggplot(p2_df, aes(percentile, value, color = metric)) +
    geom_line(linewidth = 0.9) +
    labs(
      title = paste0(title_prefix, ": empirical and gene-wise-NB IOD/CV² trajectories"),
      x = "EVS percentile",
      y = "Median bin value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
  if (is.finite(cutoff_info$cutoff_percentile)) {
    ymax <- suppressWarnings(max(p2_df$value[is.finite(p2_df$value)], na.rm = TRUE))
    if (is.finite(ymax)) {
      p2 <- p2 +
        geom_vline(xintercept = cutoff_info$cutoff_percentile, linetype = 2, linewidth = 0.8) +
        annotate("text", x = cutoff_info$cutoff_percentile, y = ymax, label = cutoff_label, hjust = 0, vjust = 1, size = 3)
    }
  }

  p3 <- ggplot(p3_df, aes(percentile, value, color = metric)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.9) +
    labs(
      title = paste0(title_prefix, ": empirical/gene-wise-NB difference and log-ratio curves"),
      x = "EVS percentile",
      y = "Median bin value"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
  if (is.finite(cutoff_info$cutoff_percentile)) {
    ymax <- suppressWarnings(max(p3_df$value[is.finite(p3_df$value)], na.rm = TRUE))
    if (is.finite(ymax)) {
      p3 <- p3 +
        geom_vline(xintercept = cutoff_info$cutoff_percentile, linetype = 2, linewidth = 0.8)
    }
  }

  png(out_file, width = 2200, height = 2000, res = 200)
  gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
  dev.off()
}

build_range_panel <- function(ctrl_sum, trt_sum, ctrl_cutoff, trt_cutoff, cmp_name, out_file) {
  diff_df <- dplyr::bind_rows(
    dplyr::transmute(ctrl_sum, percentile, value = diff_emp, dataset = "Control empirical difference"),
    dplyr::transmute(trt_sum, percentile, value = diff_emp, dataset = "Treatment empirical difference")
  )

  cut_vals <- c(ctrl_cutoff$cutoff_percentile, trt_cutoff$cutoff_percentile)
  cut_vals <- cut_vals[is.finite(cut_vals)]
  range_min <- if (length(cut_vals)) min(cut_vals) else NA_real_
  range_max <- if (length(cut_vals)) max(cut_vals) else NA_real_

  p <- ggplot(diff_df, aes(percentile, value, color = dataset)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 1.0) +
    labs(
      title = paste0(cmp_name, ": treatment/control empirical difference curves and cutoff range"),
      subtitle = "Cutoff is defined from the leading-edge side using the first sustained empirical IOD > CV² region; range spans treatment and control cutoffs",
      x = "EVS percentile",
      y = "Empirical difference (IOD - CV²)"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")

  if (is.finite(range_min) && is.finite(range_max)) {
    p <- p + annotate("rect", xmin = range_min, xmax = range_max, ymin = -Inf, ymax = Inf, alpha = 0.08)
  }
  if (is.finite(ctrl_cutoff$cutoff_percentile)) {
    p <- p + geom_vline(xintercept = ctrl_cutoff$cutoff_percentile, linetype = 2, linewidth = 0.8)
  }
  if (is.finite(trt_cutoff$cutoff_percentile)) {
    p <- p + geom_vline(xintercept = trt_cutoff$cutoff_percentile, linetype = 3, linewidth = 0.8)
  }

  ymax <- suppressWarnings(max(diff_df$value[is.finite(diff_df$value)], na.rm = TRUE))
  if (is.finite(ymax)) {
    if (is.finite(ctrl_cutoff$cutoff_percentile)) {
      p <- p + annotate(
        "text",
        x = ctrl_cutoff$cutoff_percentile,
        y = ymax,
        hjust = 0,
        vjust = 1,
        size = 3,
        label = paste0("Control cutoff\np=", sprintf("%.3f", ctrl_cutoff$cutoff_percentile), "\nr=", ctrl_cutoff$cutoff_rank_index)
      )
    }
    if (is.finite(trt_cutoff$cutoff_percentile)) {
      p <- p + annotate(
        "text",
        x = trt_cutoff$cutoff_percentile,
        y = ymax * 0.85,
        hjust = 0,
        vjust = 1,
        size = 3,
        label = paste0("Treatment cutoff\np=", sprintf("%.3f", trt_cutoff$cutoff_percentile), "\nr=", trt_cutoff$cutoff_rank_index)
      )
    }
  }

  png(out_file, width = 2200, height = 1200, res = 200)
  print(p)
  dev.off()

  list(range_min = range_min, range_max = range_max)
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

overall_summary_rows <- list()

for (i in seq_len(nrow(comparison_table))) {
  comparison_row <- comparison_table[i, , drop = FALSE]
  cmp_name <- comparison_row$comparison_name[[1]]
  message("Processing ", cmp_name, "...")

  cmp_dir <- file.path(out_root, paste0(cmp_name, "_cutoff_folder"))
  dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

  comp <- subset_comparison(count_mat, comparison_row, meta_all)
  raw_cmp_counts <- comp$count_matrix
  col_data <- comp$coldata

  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(raw_cmp_counts)),
    colData = col_data,
    design = ~ condition
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, quiet = TRUE)

  kept_ids <- rownames(dds)
  raw_cmp_counts <- raw_cmp_counts[kept_ids, , drop = FALSE]

  evs_tbl <- build_evs_table(raw_cmp_counts, comp$ctrl_ids, comp$trt_ids, kept_ids)
  alpha_tbl <- extract_genewise_dispersion(dds)
  metrics_tbl <- compute_feature_metrics(raw_cmp_counts, comp$ctrl_ids, comp$trt_ids, alpha_tbl)
  full_tbl <- dplyr::left_join(evs_tbl, metrics_tbl, by = "feature_id") %>%
    dplyr::left_join(annot_df, by = c("feature_id" = "feature_id")) %>%
    dplyr::arrange(rank_ctrl, rank_trt)

  utils::write.csv(full_tbl, file.path(cmp_dir, paste0(cmp_name, "_feature_level_metrics.csv")), row.names = FALSE)

  ctrl_obj <- build_dataset_summary(full_tbl, "control", percentile_step, min_bin_n)
  trt_obj <- build_dataset_summary(full_tbl, "treatment", percentile_step, min_bin_n)

  ctrl_cutoff <- select_leading_edge_cutoff(ctrl_obj$sum_df, "diff_emp", "iod_emp", "cv2_emp", sustained_bins)
  trt_cutoff <- select_leading_edge_cutoff(trt_obj$sum_df, "diff_emp", "iod_emp", "cv2_emp", sustained_bins)

  utils::write.csv(ctrl_obj$rank_tbl, file.path(cmp_dir, paste0(cmp_name, "_control_feature_level_metrics.csv")), row.names = FALSE)
  utils::write.csv(trt_obj$rank_tbl, file.path(cmp_dir, paste0(cmp_name, "_treatment_feature_level_metrics.csv")), row.names = FALSE)
  utils::write.csv(ctrl_obj$sum_df, file.path(cmp_dir, paste0(cmp_name, "_control_percentile_median_summary.csv")), row.names = FALSE)
  utils::write.csv(trt_obj$sum_df, file.path(cmp_dir, paste0(cmp_name, "_treatment_percentile_median_summary.csv")), row.names = FALSE)

  build_dataset_panel(
    ctrl_obj$sum_df,
    paste0(cmp_name, " control"),
    ctrl_cutoff,
    file.path(cmp_dir, paste0(cmp_name, "_control_cutoff_panel.png"))
  )

  build_dataset_panel(
    trt_obj$sum_df,
    paste0(cmp_name, " treatment"),
    trt_cutoff,
    file.path(cmp_dir, paste0(cmp_name, "_treatment_cutoff_panel.png"))
  )

  range_info <- build_range_panel(
    ctrl_obj$sum_df,
    trt_obj$sum_df,
    ctrl_cutoff,
    trt_cutoff,
    cmp_name,
    file.path(cmp_dir, paste0(cmp_name, "_cutoff_range_panel.png"))
  )

  cutoff_tbl <- dplyr::tibble(
    comparison = cmp_name,
    dataset = c("control", "treatment"),
    cutoff_mode = c(ctrl_cutoff$mode, trt_cutoff$mode),
    cutoff_percentile = c(ctrl_cutoff$cutoff_percentile, trt_cutoff$cutoff_percentile),
    cutoff_rank = c(ctrl_cutoff$cutoff_rank_index, trt_cutoff$cutoff_rank_index),
    note = c(ctrl_cutoff$note, trt_cutoff$note)
  )
  utils::write.csv(cutoff_tbl, file.path(cmp_dir, paste0(cmp_name, "_cutoff_values.csv")), row.names = FALSE)

  summary_row <- dplyr::tibble(
    comparison = cmp_name,
    n_features = nrow(full_tbl),
    control_cutoff_percentile = ctrl_cutoff$cutoff_percentile,
    control_cutoff_rank = ctrl_cutoff$cutoff_rank_index,
    treatment_cutoff_percentile = trt_cutoff$cutoff_percentile,
    treatment_cutoff_rank = trt_cutoff$cutoff_rank_index,
    cutoff_range_percentile_min = range_info$range_min,
    cutoff_range_percentile_max = range_info$range_max,
    cutoff_range_rank_min = min(c(ctrl_cutoff$cutoff_rank_index, trt_cutoff$cutoff_rank_index), na.rm = TRUE),
    cutoff_range_rank_max = max(c(ctrl_cutoff$cutoff_rank_index, trt_cutoff$cutoff_rank_index), na.rm = TRUE)
  )
  if (!is.finite(summary_row$cutoff_range_rank_min)) summary_row$cutoff_range_rank_min <- NA_real_
  if (!is.finite(summary_row$cutoff_range_rank_max)) summary_row$cutoff_range_rank_max <- NA_real_

  utils::write.csv(summary_row, file.path(cmp_dir, paste0(cmp_name, "_summary_row.csv")), row.names = FALSE)
  overall_summary_rows[[cmp_name]] <- summary_row
}

overall_summary_tbl <- dplyr::bind_rows(overall_summary_rows)
utils::write.csv(overall_summary_tbl, file.path(out_root, "overall_cutoff_summary.csv"), row.names = FALSE)
message("Done. Outputs written to: ", out_root)
