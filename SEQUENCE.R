# =============================================================================
# SEQUENCE MANUSCRIPT PIPELINE
# =============================================================================
# Purpose
# This script implements the manuscript analysis workflow for the SEQUENCE WTTS-
# Seq study. It reads raw count data from the repository data/ directory, builds
# treatment-specific and control-specific eigenvector-ranked feature tables,
# derives DESeq2-informed dispersion summaries, estimates ranked local Fourier
# wave structure for the DESeq2-derived index of dispersion (IOD) and squared
# coefficient of variation (CV2), combines the treatment and control local
# Fourier maps into one comparison-level transition score, applies one shared
# EVS cutoff per comparison, performs differential-expression testing on the
# original, leading-edge, and remainder subsets, computes Hybrid Bayesian-
# Frequentist Significance Scores (HBFSS), and exports manuscript figures and
# tables into a timestamped exports/ directory inside the same repository.
#
# Method overview
# 1. For each comparison, treatment and control datasets are prepared from the
#    raw count matrix and normalized independently for EVS ranking.
# 2. Features are ranked by absolute PC1 loading within treatment and control.
# 3. DESeq2-derived IOD and CV2 are mapped onto those ranked treatment and
#    control axes.
# 4. At percentile positions spanning the ranked dataset, local windows are
#    evaluated with low-order Fourier series for IOD and CV2.
# 5. Treatment and control local Fourier summaries are combined into one
#    comparison-level score that identifies the dominant transition interval.
# 6. One shared comparison-level cutoff rank is selected from that combined
#    Fourier score and projected back onto the normalized and raw EVS panels.
# 7. The selected cutoff defines the leading-edge and remainder subsets, which
#    are then analyzed alongside the original dataset.
# 8. Manuscript figures, wave diagnostics, and feature-level result tables are
#    written to the repository exports/ directory.
#
# Design principle
# Treatment and control are allowed to exhibit different local wave behavior,
# but downstream EVS splitting uses one shared cutoff per comparison because the
# original, leading-edge, and remainder analyses require a single comparison-
# level split. That shared cutoff is derived from the combined treatment-plus-
# control local Fourier score rather than from a treatment-only or control-only
# decision rule.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(fdrtool)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(gridExtra)
  library(grid)
  library(scales)
  library(grDevices)
  library(S4Vectors)
})

# -----------------------------------------------------------------------------
# USER INPUT
# -----------------------------------------------------------------------------

repo_dir <- getwd()
input_dir <- file.path(repo_dir, "data")
output_root <- file.path(repo_dir, "exports")
analysis_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(output_root, paste0("sequence_run_", analysis_stamp))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")

run_twas_overlap <- FALSE
twas_file <- file.path(input_dir, "3aTWAS_genes_of_11_brain_disorders.csv")

alpha_level <- 0.20
max_usable_hc_p_threshold <- 0.99
lfc_boundary <- 1.0
lfc_shrink_type <- "apeglm"   # "normal" for speed/stability; "apeglm" optional
lfc_shrink_apeglm_method <- "nbinomC"

evs_cutoff_mode_main <- "matched_curvature_wave"
evs_fixed_top_n <- 5000L
evs_fixed_rank_override <- NA_integer_
evs_use_manual_rank_first <- FALSE

evs_curvature_spar <- 0.55
evs_curvature_min_rank <- 50L
evs_curvature_max_rank_frac <- 0.50
evs_curvature_spar_iod <- 0.55
evs_curvature_spar_cv2 <- 0.55
evs_curvature_match_window <- 150L
evs_candidate_interval_halfwidth <- 75L
evs_candidate_merge_distance <- 100L
wave_local_window <- 51L
wave_w_iod_slope <- 1.00
wave_w_iod_amp   <- 0.50
wave_w_cv2_amp   <- 1.00
wave_w_cv2_slope <- 0.50

evs_candidate_min_peak_distance   <- 40L
evs_peak_quantile_grid            <- seq(0.90, 0.75, by = -0.05)
evs_peak_quantile_screen_min      <- 0.75
evs_peak_quantile_screen_max      <- 0.90
evs_candidate_max_leading_frac     <- 0.50
evs_candidate_runtime_seconds_low  <- 12
evs_candidate_runtime_seconds_high <- 25

fourier_percentile_step <- 0.01
fourier_window_fraction <- 0.12
fourier_harmonics <- 2L
fourier_top_candidate_n_for_deseq2 <- 1L
fourier_min_rank <- 1L
fourier_max_rank_frac <- 1.00
fourier_score_weight_iod_amp <- 1.0
fourier_score_weight_cv2_amp <- 1.0
fourier_score_weight_center_agreement <- 1.0
fourier_interval_fraction_of_max <- 0.90

crossing_min_percentile <- 0.05
crossing_max_percentile <- 0.95
crossing_stability_window_n <- 3L
crossing_plot_line_width <- 0.95
crossing_plot_vline_width <- 0.95
crossing_plot_hline_width <- 0.65
crossing_plot_point_size <- 1.8

figure_dpi <- 320
base_theme_size <- 10

# -----------------------------------------------------------------------------
# RAW-WALD BRACKETED CUTOFF SEARCH
# -----------------------------------------------------------------------------

raw_wald_cutoff_mode <- TRUE
raw_wald_search_coarse_step <- 25L
raw_wald_search_fine_radius <- 50L
raw_wald_search_balance_lambda <- 0.50
raw_wald_search_min_subset_n <- 200L

raw_wald_score_weight_mean <- 1.00
raw_wald_score_weight_sd   <- 1.00
raw_wald_score_weight_skew <- 0.50

# =============================================================================
# UNIFIED AESTHETIC CONSTANTS
# =============================================================================

POINT_SIZE_PRIMARY  <- 1.6
POINT_SIZE_DISP     <- 1.3
POINT_ALPHA_PRIMARY <- 0.82
POINT_ALPHA_DISP    <- 0.55
POINT_STROKE        <- 0.40

LINE_WIDTH_BOUNDARY <- 0.55
LINE_WIDTH_ZERO     <- 0.40
LINE_WIDTH_THRESH   <- 0.90

HIST_BINS  <- 60
HIST_COLOR <- "white"

plot_palette <- list(
  background   = "#BDBDBD",
  threshold    = "#8C2D04",
  hbfss        = "#E67E22",
  deseq2       = "#C0392B",
  strong       = "#C0392B",
  overlap      = "#7D3C98",
  weak         = "#4A90E2",
  intermediate = "#7F7F7F",
  histogram    = "#969696",
  control      = "#4D4D4D",
  treatment    = "#1F78B4"
)

# -----------------------------------------------------------------------------
# FIGURE EXPORT CONTROLS
# -----------------------------------------------------------------------------

export_optional_mean_histograms        <- FALSE
export_optional_empirical_p_histograms <- FALSE
export_optional_hbfss_distributions    <- FALSE
export_optional_evs_raw_histograms     <- FALSE
export_optional_evs_variance_profiles  <- TRUE

# -----------------------------------------------------------------------------
# CODE NOTATION / COLUMN GLOSSARY
# -----------------------------------------------------------------------------
# feature_id     : OrigID / PAS identifier from the WTTS count matrix.
# gene_symbol    : Symbol column from the WTTS count matrix.
# stat           : DESeq2 Wald test statistic.
# pvalue         : raw DESeq2 Wald-test p-value.
# padj           : native DESeq2 Benjamini-Hochberg adjusted p-value.
# empirical_p    : fdrtool / Strimmer empirical-null p-value.
# empirical_q    : fdrtool q-value.
# lfdr           : fdrtool local false discovery rate.
# lfc_shrunk     : apeglm-shrunken log2 fold change.
# resGA_padj     : DESeq2 adjusted p-value from the greaterAbs composite-null test.
# resLA_padj     : DESeq2 adjusted p-value from the lessAbs composite-null test.
# hc_p_threshold_dataset  : Strimmer HC threshold from hc.thresh(sort(empirical_p)).
# hbfss_threshold_dataset : abs(log10(hc_p_threshold_dataset)) * lfc_boundary.
# HBFSS / pi_valueE       : abs(apeglm-shrunken log2FC * log10(empirical_p)).
# standard_significant    : final DC2-positive call used for display.
# comparison-level cutoff : one shared combined-rank split per comparison.

# -----------------------------------------------------------------------------
# MINIMAL EMBEDDED METADATA
# -----------------------------------------------------------------------------

meta_all <- data.frame(
  id = c(
    "R0_1","R0_2","R0_3","R0_4","R0_5","ZT6_1","ZT6_2","ZT6_3","ZT6_4","ZT6_5",
    "R2_1","R2_2","R2_3","R2_4","R2_5","ZT8_1","ZT8_2","ZT8_3","ZT8_4","ZT8_5",
    "R4_1","R4_2","R4_3","R4_4","R4_5","ZT10_1","ZT10_2","ZT10_3","ZT10_4","ZT10_5",
    "R8_1","R8_2","R8_3","R8_4","R8_5","ZT14_1","ZT14_2","ZT14_3","ZT14_4","ZT14_5"
  ),
  condition = c(
    "treatment","treatment","treatment","treatment","treatment",
    "control","control","control","control","control",
    "treatment","treatment","treatment","treatment","treatment",
    "control","control","control","control","control",
    "treatment","treatment","treatment","treatment","treatment",
    "control","control","control","control","control",
    "treatment","treatment","treatment","treatment","treatment",
    "control","control","control","control","control"
  ),
  stringsAsFactors = FALSE
)

rownames(meta_all) <- meta_all$id
meta_all$condition <- factor(meta_all$condition, levels = c("control", "treatment"))
levels(meta_all$condition) <- c("untrt", "trt")

# -----------------------------------------------------------------------------
# FOUR PAIRWISE COMPARISONS
# -----------------------------------------------------------------------------

comparison_table <- data.frame(
  comparison_name  = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix    = c("R0", "R2", "R4", "R8"),
  group2_prefix    = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

dataset_key_order <- c("raw_dataset", "leading_edge_dataset", "remainder_dataset")

dataset_key_labels <- c(
  raw_dataset          = "Original dataset",
  leading_edge_dataset = "Leading-edge dataset",
  remainder_dataset    = "Remainder dataset"
)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

assert_required_columns <- function(df, required_cols, object_name = "data frame") {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns in ", object_name, ": ",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

safe_log10 <- function(x, pseudocount = 1e-12) {
  log10(pmax(x, pseudocount))
}

safe_neglog10 <- function(x, pseudocount = 1e-12) {
  -log10(pmax(x, pseudocount))
}

preprocessing_panel_label <- function(x) {
  x <- as.character(x)[1]
  if (is.na(x) || !nzchar(x)) return("Preprocessing not specified")
  dplyr::case_when(
    x %in% c("DESeq2-normalized counts", "normalized", "Normalized prior to eigenvector splitting", "Normalized before EVS") ~
      "Normalized prior to eigenvector splitting",
    x %in% c("raw counts without DESeq2 normalization", "raw_counts", "Eigenvector splitting without prior normalization", "Raw counts before EVS") ~
      "Eigenvector splitting without prior normalization",
    TRUE ~ x
  )
}

evs_preproc_short <- function(x) {
  x <- preprocessing_panel_label(x)
  if (identical(x, "Normalized prior to eigenvector splitting")) {
    "Normalized before EVS"
  } else if (identical(x, "Eigenvector splitting without prior normalization")) {
    "Raw counts before EVS"
  } else {
    x
  }
}

compute_mean_expression_table <- function(raw_counts, coldata) {
  sample_ids <- colnames(raw_counts)
  trt_ids    <- sample_ids[coldata$condition == "trt"]
  untrt_ids  <- sample_ids[coldata$condition == "untrt"]

  data.frame(
    feature_id       = as.character(rownames(raw_counts)),
    mean_trt         = rowMeans(raw_counts[, trt_ids,   drop = FALSE]),
    mean_untrt       = rowMeans(raw_counts[, untrt_ids, drop = FALSE]),
    mean_all         = rowMeans(raw_counts),
    stringsAsFactors = FALSE
  )
}

plot_mean_histogram_panel <- function(df_means, base_mean_vec, dataset_name) {
  make_hist <- function(vals, panel_title) {
    ggplot(data.frame(x = safe_log10(vals + 1)), aes(x)) +
      geom_histogram(bins = HIST_BINS, fill = plot_palette$histogram, color = HIST_COLOR) +
      labs(title = paste(dataset_name, panel_title), x = "log10(mean + 1)", y = "Count") +
      manuscript_theme()
  }

  arrangeGrob(
    make_hist(df_means$mean_trt,   "Treatment mean"),
    make_hist(df_means$mean_untrt, "Control mean"),
    make_hist(df_means$mean_all,   "Pooled mean"),
    make_hist(base_mean_vec,       "DESeq2 baseMean"),
    ncol = 2,
    top  = textGrob(
      paste(dataset_name, "Mean-Expression Histograms"),
      gp = gpar(fontface = "bold", cex = 1.2)
    )
  )
}

compact_title <- function(x, width = 58) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

compact_caption <- function(x, width = 120) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

clip_probabilities <- function(x, eps = 1e-300) {
  x <- unname(as.numeric(x))
  if (!length(x)) return(numeric(0))

  missing_idx <- is.na(x)
  pos_inf_idx <- is.infinite(x) & x > 0
  neg_inf_idx <- is.infinite(x) & x < 0

  x[pos_inf_idx] <- 1 - 1e-12
  x[neg_inf_idx] <- eps

  finite_idx <- is.finite(x) & !missing_idx
  x[finite_idx] <- pmin(pmax(x[finite_idx], eps), 1 - 1e-12)
  x[missing_idx] <- NA_real_
  x
}

run_empirical_null_fdrtool <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  stat_vec <- unname(stat_vec)

  n_stats <- length(stat_vec)
  if (n_stats == 0L) {
    return(list(
      pval = numeric(0),
      qval = numeric(0),
      lfdr = numeric(0),
      method = "empty_input"
    ))
  }

  if (n_stats < 5L) {
    fallback_p <- clip_probabilities(2 * stats::pnorm(-abs(stat_vec)))
    return(list(
      pval = fallback_p,
      qval = clip_probabilities(stats::p.adjust(fallback_p, method = "BH")),
      lfdr = rep(NA_real_, length(fallback_p)),
      method = "normal_fallback_small_n"
    ))
  }

  fit_attempts <- list(
    list(cutoff.method = "fndr", pct0 = 0.75),
    list(cutoff.method = "pct0"),
    list(cutoff.method = "fndr"),
    list(cutoff.method = "locfdr")
  )

  for (ii in seq_along(fit_attempts)) {
    attempt <- fit_attempts[[ii]]
    fit <- tryCatch(
      {
        do.call(
          fdrtool,
          c(
            list(
              x = stat_vec,
              statistic = "normal",
              plot = FALSE,
              verbose = FALSE
            ),
            attempt
          )
        )
      },
      error = function(e) {
        message(sprintf("[%s] fdrtool attempt %d failed: %s", dataset_name, ii, conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(fit) &&
        !is.null(fit$pval) && length(fit$pval) == n_stats &&
        !is.null(fit$qval) && length(fit$qval) == n_stats) {
      fit$pval <- clip_probabilities(fit$pval)
      fit$qval <- clip_probabilities(fit$qval)
      fit$lfdr <- if (!is.null(fit$lfdr)) as.numeric(fit$lfdr) else rep(NA_real_, n_stats)
      fit$method <- paste0("fdrtool_attempt_", ii)
      return(fit)
    }
  }

  message(sprintf("[%s] All fdrtool attempts failed. Falling back to two-sided normal p-values.", dataset_name))
  fallback_p <- clip_probabilities(2 * stats::pnorm(-abs(stat_vec)))

  list(
    pval = fallback_p,
    qval = clip_probabilities(stats::p.adjust(fallback_p, method = "BH")),
    lfdr = rep(NA_real_, length(fallback_p)),
    method = "normal_fallback_full"
  )
}

safe_hc_thresh <- function(empirical_p, dataset_name) {
  sorted_empirical_p <- sort(
    clip_probabilities(empirical_p),
    na.last    = NA,
    decreasing = FALSE
  )
  if (length(sorted_empirical_p) < 5) {
    return(min(max_usable_hc_p_threshold * 0.99, max(1e-6, alpha_level)))
  }

  out <- suppressWarnings(
    tryCatch(
      fdrtool::hc.thresh(as.vector(sorted_empirical_p)),
      error = function(e) {
        message(sprintf("[%s] hc.thresh failed: %s", dataset_name, conditionMessage(e)))
        NA_real_
      }
    )
  )

  out <- as.numeric(out[1])

  if (!is.finite(out) || out <= 0 || out >= max_usable_hc_p_threshold) {
    fallback_out <- suppressWarnings(as.numeric(stats::quantile(
      sorted_empirical_p,
      probs = min(0.10, max(0.05, alpha_level)),
      na.rm = TRUE,
      type = 7
    )[1]))

    fallback_out <- min(max_usable_hc_p_threshold * 0.99, max(1e-6, fallback_out))
    message(sprintf(
      "[%s] HC threshold fallback used for HBFSS calibration (primary=%s; fallback=%s)",
      dataset_name,
      ifelse(is.finite(out), signif(out, 6), "NA"),
      signif(fallback_out, 6)
    ))
    return(fallback_out)
  }

  out
}

plot_expand_xy <- function() {

  list(
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.24))),
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.30)))
  )
}

resolve_top_n_cutoff <- function(sorted_values_desc, top_n = evs_fixed_top_n) {
  n_total <- length(sorted_values_desc)
  if (n_total == 0) stop("resolve_top_n_cutoff() received an empty vector.", call. = FALSE)

  top_n_actual    <- min(max(1L, as.integer(top_n)), n_total)
  cutoff_value    <- sorted_values_desc[top_n_actual]
  cutoff_quantile <- 1 - (top_n_actual / n_total)

  list(
    top_n_actual    = top_n_actual,
    cutoff_value    = cutoff_value,
    cutoff_quantile = cutoff_quantile,
    n_total         = n_total,
    method          = "fixed_top_n",
    curve_df        = NULL,
    curvature_strength = NA_real_,
    variance_measure = "NB2_variance",
    candidate_table = data.frame(
      candidate_id = "fixed_top_n",
      rank_index = top_n_actual,
      curvature_strength = NA_real_,
      prominence = NA_real_,
      smoothed_log_nb_variance = NA_real_,
      first_derivative = NA_real_,
      second_derivative = NA_real_,
      max_rank_allowed = NA_integer_,
      cutoff_value = cutoff_value,
      cutoff_quantile = cutoff_quantile,
      selected = TRUE,
      selected_reason = "fixed_top_n",
      stringsAsFactors = FALSE
    ),
    selected_reason = "fixed_top_n"
  )
}

build_nb_variance_curve <- function(loading_tbl,
                                    mean_col = "baseMean",
                                    dispersion_col = "dispGeneEst") {
  needed <- c("feature_id", "rank", "pc1_loading_abs", mean_col, dispersion_col)
  if (!all(needed %in% colnames(loading_tbl))) return(NULL)

  df <- loading_tbl[, needed, drop = FALSE]
  colnames(df)[colnames(df) == mean_col] <- "mean_value"
  colnames(df)[colnames(df) == dispersion_col] <- "dispersion_value"

  df <- df[
    is.finite(df$mean_value) & !is.na(df$mean_value) & df$mean_value >= 0 &
      is.finite(df$dispersion_value) & !is.na(df$dispersion_value) & df$dispersion_value > 0,
    , drop = FALSE
  ]
  if (nrow(df) < 25) return(NULL)

  df <- df[order(df$rank), , drop = FALSE]
  df$nb_variance <- with(df, mean_value + dispersion_value * (mean_value ^ 2))
  df$variance_measure <- "NB2_variance"
  df$curve_metric <- "nb2_variance"
  df$curve_label <- "log10 NB2 variance"

  df <- df[is.finite(df$nb_variance) & !is.na(df$nb_variance) & df$nb_variance > 0, , drop = FALSE]
  if (nrow(df) < 25) return(NULL)

  df$log_nb_variance <- log10(df$nb_variance)
  df
}

smooth_ranked_nb_variance_curve <- function(loading_tbl,
                                            mean_col = "baseMean",
                                            dispersion_col = "dispGeneEst",
                                            spar = evs_curvature_spar) {
  df <- build_nb_variance_curve(
    loading_tbl = loading_tbl,
    mean_col = mean_col,
    dispersion_col = dispersion_col
  )
  if (is.null(df) || nrow(df) < 25) return(NULL)

  fit <- tryCatch(
    smooth.spline(x = df$rank, y = df$log_nb_variance, spar = spar),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  pred0 <- tryCatch(predict(fit, x = df$rank, deriv = 0), error = function(e) NULL)
  pred1 <- tryCatch(predict(fit, x = df$rank, deriv = 1), error = function(e) NULL)
  pred2 <- tryCatch(predict(fit, x = df$rank, deriv = 2), error = function(e) NULL)
  if (is.null(pred0) || is.null(pred1) || is.null(pred2)) return(NULL)

  df$smoothed_log_nb_variance <- as.numeric(pred0$y)
  df$first_derivative <- as.numeric(pred1$y)
  df$second_derivative <- as.numeric(pred2$y)
  df$curvature <- abs(df$second_derivative) / ((1 + df$first_derivative^2)^(3/2))
  df$curvature[!is.finite(df$curvature)] <- NA_real_
  df
}

find_local_curvature_peaks <- function(curve_df,
                                       min_rank = evs_curvature_min_rank,
                                       max_rank_frac = evs_curvature_max_rank_frac,
                                       edge_buffer = 5L) {
  if (is.null(curve_df) || nrow(curve_df) < 25) return(NULL)

  x <- as.numeric(curve_df$rank)
  y <- as.numeric(curve_df$smoothed_log_nb_variance)
  d1 <- as.numeric(curve_df$first_derivative)
  d2 <- as.numeric(curve_df$second_derivative)
  curvature <- as.numeric(curve_df$curvature)

  keep <- is.finite(x) & is.finite(y) & is.finite(d1) & is.finite(d2) & is.finite(curvature)
  x <- x[keep]
  y <- y[keep]
  d1 <- d1[keep]
  d2 <- d2[keep]
  curvature <- curvature[keep]

  if (length(x) < 25) return(NULL)

  min_rank <- max(edge_buffer + 1L, as.integer(min_rank))
  max_rank <- min(length(x) - edge_buffer, max(min_rank, floor(length(x) * max_rank_frac)))

  idx_valid <- which(x >= min_rank & x <= max_rank & is.finite(curvature))
  if (!length(idx_valid)) return(NULL)

  curv_valid <- curvature[idx_valid]
  prev_vals <- c(-Inf, head(curv_valid, -1))
  next_vals <- c(tail(curv_valid, -1), -Inf)
  local_peak_pos <- which(curv_valid >= prev_vals & curv_valid > next_vals)

  if (length(local_peak_pos) > 1L) {
    local_peak_pos <- local_peak_pos[local_peak_pos > 1L & local_peak_pos < length(curv_valid)]
  } else if (length(local_peak_pos) == 1L) {
    if (local_peak_pos <= 1L || local_peak_pos >= length(curv_valid)) {
      local_peak_pos <- integer(0)
    }
  }

  if (!length(local_peak_pos)) {
    if (length(curv_valid) > 2L) {
      interior_idx <- 2L:(length(curv_valid) - 1L)
      local_peak_pos <- interior_idx[which.max(curv_valid[interior_idx])]
    } else {
      local_peak_pos <- which.max(curv_valid)
    }
  }

  local_peak_idx <- idx_valid[local_peak_pos]

  peak_tbl <- data.frame(
    idx = local_peak_idx,
    rank_index = as.integer(round(x[local_peak_idx])),
    curvature_strength = as.numeric(curvature[local_peak_idx]),
    smoothed_log_nb_variance = as.numeric(y[local_peak_idx]),
    first_derivative = as.numeric(d1[local_peak_idx]),
    second_derivative = as.numeric(d2[local_peak_idx]),
    stringsAsFactors = FALSE
  )

  peak_tbl$prominence <- peak_tbl$curvature_strength - pmax(
    curvature[pmax(1L, peak_tbl$idx - 1L)],
    curvature[pmin(length(curvature), peak_tbl$idx + 1L)]
  )

  peak_tbl <- peak_tbl[order(peak_tbl$rank_index), , drop = FALSE]
  rownames(peak_tbl) <- NULL

  list(
    peak_table = peak_tbl,
    max_rank_allowed = as.integer(round(max_rank)),
    x = x,
    y = y,
    d1 = d1,
    d2 = d2,
    curvature = curvature
  )
}

summarize_curve_peak_grid <- function(curve_df,
                                      min_rank = evs_curvature_min_rank,
                                      max_rank_frac = evs_curvature_max_rank_frac,
                                      quantile_grid = evs_peak_quantile_grid,
                                      min_peak_distance = 1L) {
  peak_info <- find_local_curvature_peaks(
    curve_df = curve_df,
    min_rank = min_rank,
    max_rank_frac = max_rank_frac
  )
  if (is.null(peak_info) || is.null(peak_info$peak_table) || !nrow(peak_info$peak_table)) {
    return(data.frame())
  }

  peak_tbl <- peak_info$peak_table
  peak_tbl <- peak_tbl[order(peak_tbl$rank_index, -peak_tbl$curvature_strength, -peak_tbl$prominence), , drop = FALSE]
  if (!nrow(peak_tbl)) return(data.frame())

  keep_by_distance <- function(tbl, min_distance) {
    if (!nrow(tbl) || min_distance <= 1L) return(tbl)
    tbl <- tbl[order(tbl$rank_index, -tbl$curvature_strength, -tbl$prominence), , drop = FALSE]
    keep <- rep(TRUE, nrow(tbl))
    last_kept <- -Inf
    for (ii in seq_len(nrow(tbl))) {
      if ((tbl$rank_index[ii] - last_kept) < min_distance) {
        keep[ii] <- FALSE
      } else {
        last_kept <- tbl$rank_index[ii]
      }
    }
    tbl[keep, , drop = FALSE]
  }

  out <- lapply(quantile_grid, function(qv) {
    threshold <- as.numeric(stats::quantile(
      peak_tbl$curvature_strength,
      probs = qv,
      na.rm = TRUE,
      type = 7
    ))
    survivors <- peak_tbl[
      is.finite(peak_tbl$curvature_strength) &
        peak_tbl$curvature_strength >= threshold,
      ,
      drop = FALSE
    ]
    survivors <- keep_by_distance(survivors, as.integer(min_peak_distance))
    if (!nrow(survivors)) return(NULL)
    survivors <- survivors[order(survivors$rank_index, -survivors$curvature_strength, -survivors$prominence), , drop = FALSE]
    survivors$peak_quantile_cutoff <- qv
    survivors$peak_strength_threshold <- threshold
    survivors$peak_order_within_quantile <- seq_len(nrow(survivors))
    survivors$percent_rank <- 100 * survivors$rank_index / nrow(curve_df)
    survivors[, c(
      "peak_quantile_cutoff",
      "peak_strength_threshold",
      "peak_order_within_quantile",
      "rank_index",
      "percent_rank",
      "curvature_strength",
      "prominence",
      "smoothed_log_nb_variance",
      "first_derivative",
      "second_derivative"
    ), drop = FALSE]
  })

  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(data.frame())
  dplyr::bind_rows(out)
}

select_quantile_screen_candidates <- function(peak_grid_df,
                                              q_min = evs_peak_quantile_screen_min,
                                              q_max = evs_peak_quantile_screen_max) {
  if (is.null(peak_grid_df) || !nrow(peak_grid_df)) return(integer(0))

  q_lo <- min(q_min, q_max, na.rm = TRUE)
  q_hi <- max(q_min, q_max, na.rm = TRUE)

  screen_tbl <- peak_grid_df[
    is.finite(peak_grid_df$peak_quantile_cutoff) &
      peak_grid_df$peak_quantile_cutoff >= q_lo &
      peak_grid_df$peak_quantile_cutoff <= q_hi,
    , drop = FALSE
  ]
  if (!nrow(screen_tbl)) return(integer(0))

  screen_tbl <- screen_tbl[order(screen_tbl$peak_quantile_cutoff, screen_tbl$peak_order_within_quantile, screen_tbl$rank_index), , drop = FALSE]
  unique(as.integer(screen_tbl$rank_index))
}

# -----------------------------------------------------------------------------
# MANUSCRIPT METHOD: RANKED DISPERSION-METRIC TABLE
# -----------------------------------------------------------------------------

build_ranked_regime_metric_table <- function(loading_tbl,
                                             mean_col = "baseMean",
                                             dispersion_col = "dispGeneEst") {
  needed <- c("feature_id", "rank", "pc1_loading_abs", mean_col, dispersion_col)
  if (!all(needed %in% colnames(loading_tbl))) return(NULL)

  df <- loading_tbl[, needed, drop = FALSE]
  colnames(df)[colnames(df) == mean_col] <- "mean_value"
  colnames(df)[colnames(df) == dispersion_col] <- "dispersion_value"

  df <- df[
    is.finite(df$mean_value) & !is.na(df$mean_value) & df$mean_value > 0 &
      is.finite(df$dispersion_value) & !is.na(df$dispersion_value) & df$dispersion_value > 0,
    , drop = FALSE
  ]
  if (nrow(df) < 25) return(NULL)

  df <- df[order(df$rank), , drop = FALSE]
  mu <- pmax(df$mean_value, 1e-12)
  alpha <- pmax(df$dispersion_value, 1e-12)

  df$iod_nb <- 1 + alpha * mu
  df$cv2_nb <- (1 / mu) + alpha
  df$log_iod_nb <- log10(df$iod_nb)
  df$log_cv2_nb <- log10(df$cv2_nb)
  df
}

smooth_ranked_regime_curves <- function(loading_tbl,
                                        mean_col = "baseMean",
                                        dispersion_col = "dispGeneEst",
                                        spar_iod = evs_curvature_spar_iod,
                                        spar_cv2 = evs_curvature_spar_cv2) {
  df <- build_ranked_regime_metric_table(loading_tbl, mean_col, dispersion_col)
  if (is.null(df) || nrow(df) < 25) return(NULL)

  fit_iod <- tryCatch(smooth.spline(x = df$rank, y = df$log_iod_nb, spar = spar_iod), error = function(e) NULL)
  fit_cv2 <- tryCatch(smooth.spline(x = df$rank, y = df$log_cv2_nb, spar = spar_cv2), error = function(e) NULL)
  if (is.null(fit_iod) || is.null(fit_cv2)) return(NULL)

  p0_iod <- tryCatch(predict(fit_iod, x = df$rank, deriv = 0), error = function(e) NULL)
  p1_iod <- tryCatch(predict(fit_iod, x = df$rank, deriv = 1), error = function(e) NULL)
  p2_iod <- tryCatch(predict(fit_iod, x = df$rank, deriv = 2), error = function(e) NULL)
  p0_cv2 <- tryCatch(predict(fit_cv2, x = df$rank, deriv = 0), error = function(e) NULL)
  p1_cv2 <- tryCatch(predict(fit_cv2, x = df$rank, deriv = 1), error = function(e) NULL)
  p2_cv2 <- tryCatch(predict(fit_cv2, x = df$rank, deriv = 2), error = function(e) NULL)
  if (is.null(p0_iod) || is.null(p1_iod) || is.null(p2_iod) || is.null(p0_cv2) || is.null(p1_cv2) || is.null(p2_cv2)) return(NULL)

  df$smoothed_log_iod_nb <- as.numeric(p0_iod$y)
  df$iod_d1 <- as.numeric(p1_iod$y)
  df$iod_d2 <- as.numeric(p2_iod$y)
  df$iod_curvature <- abs(df$iod_d2) / ((1 + df$iod_d1^2)^(3/2))

  df$smoothed_log_cv2_nb <- as.numeric(p0_cv2$y)
  df$cv2_d1 <- as.numeric(p1_cv2$y)
  df$cv2_d2 <- as.numeric(p2_cv2$y)
  df$cv2_curvature <- abs(df$cv2_d2) / ((1 + df$cv2_d1^2)^(3/2))

  df$iod_curvature[!is.finite(df$iod_curvature)] <- NA_real_
  df$cv2_curvature[!is.finite(df$cv2_curvature)] <- NA_real_
  df
}

find_local_metric_curvature_peaks <- function(curve_df,
                                              metric = c("iod", "cv2"),
                                              min_rank = evs_curvature_min_rank,
                                              max_rank_frac = evs_curvature_max_rank_frac) {
  metric <- match.arg(metric)
  curv_col <- if (metric == "iod") "iod_curvature" else "cv2_curvature"
  d1_col <- if (metric == "iod") "iod_d1" else "cv2_d1"
  d2_col <- if (metric == "iod") "iod_d2" else "cv2_d2"
  y_col <- if (metric == "iod") "smoothed_log_iod_nb" else "smoothed_log_cv2_nb"

  x <- curve_df$rank
  y <- curve_df[[curv_col]]
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 25) return(data.frame())

  max_rank <- floor(length(x) * max_rank_frac)
  idx_valid <- which(x >= min_rank & x <= max_rank)
  if (!length(idx_valid)) return(data.frame())

  yy <- y[idx_valid]
  prev_vals <- c(-Inf, head(yy, -1))
  next_vals <- c(tail(yy, -1), -Inf)
  peak_pos <- which(yy >= prev_vals & yy > next_vals)
  if (!length(peak_pos)) peak_pos <- which.max(yy)
  out_idx <- idx_valid[peak_pos]

  out <- data.frame(
    metric = metric,
    rank_index = as.integer(round(curve_df$rank[out_idx])),
    curvature_strength = as.numeric(curve_df[[curv_col]][out_idx]),
    smoothed_value = as.numeric(curve_df[[y_col]][out_idx]),
    first_derivative = as.numeric(curve_df[[d1_col]][out_idx]),
    second_derivative = as.numeric(curve_df[[d2_col]][out_idx]),
    stringsAsFactors = FALSE
  )

  qthr <- as.numeric(stats::quantile(out$curvature_strength, probs = evs_peak_quantile_screen_min, na.rm = TRUE, type = 7))
  out <- out[out$curvature_strength >= qthr, , drop = FALSE]
  out[order(out$rank_index), , drop = FALSE]
}

match_curvature_peak_intervals <- function(iod_peaks,
                                           cv2_peaks,
                                           match_window = evs_curvature_match_window,
                                           merge_distance = evs_candidate_merge_distance) {
  if (is.null(iod_peaks) || !nrow(iod_peaks) || is.null(cv2_peaks) || !nrow(cv2_peaks)) return(data.frame())

  hits <- list()
  kk <- 1L
  for (i in seq_len(nrow(iod_peaks))) {
    for (j in seq_len(nrow(cv2_peaks))) {
      d <- abs(iod_peaks$rank_index[i] - cv2_peaks$rank_index[j])
      if (d <= match_window) {
        lo <- min(iod_peaks$rank_index[i], cv2_peaks$rank_index[j])
        hi <- max(iod_peaks$rank_index[i], cv2_peaks$rank_index[j])
        hits[[kk]] <- data.frame(
          iod_rank = iod_peaks$rank_index[i],
          cv2_rank = cv2_peaks$rank_index[j],
          interval_lo = lo,
          interval_hi = hi,
          interval_center = as.integer(round(mean(c(lo, hi)))),
          joint_strength = iod_peaks$curvature_strength[i] + cv2_peaks$curvature_strength[j],
          stringsAsFactors = FALSE
        )
        kk <- kk + 1L
      }
    }
  }
  if (!length(hits)) return(data.frame())
  out <- dplyr::bind_rows(hits)
  out <- out[order(out$interval_center, -out$joint_strength), , drop = FALSE]

  merged <- list()
  cur <- out[1, , drop = FALSE]
  kk <- 1L
  if (nrow(out) > 1L) {
    for (ii in 2:nrow(out)) {
      if ((out$interval_center[ii] - cur$interval_center[1]) <= merge_distance) {
        cur$interval_lo[1] <- min(cur$interval_lo[1], out$interval_lo[ii])
        cur$interval_hi[1] <- max(cur$interval_hi[1], out$interval_hi[ii])
        cur$interval_center[1] <- as.integer(round(mean(c(cur$interval_lo[1], cur$interval_hi[1]))))
        cur$joint_strength[1] <- max(cur$joint_strength[1], out$joint_strength[ii])
      } else {
        merged[[kk]] <- cur
        kk <- kk + 1L
        cur <- out[ii, , drop = FALSE]
      }
    }
  }
  merged[[kk]] <- cur
  dplyr::bind_rows(merged)
}

score_wave_candidate_rank <- function(curve_df, rank_index) {
  n_total <- nrow(curve_df)
  if (!is.finite(rank_index) || is.na(rank_index) || rank_index <= 2L || rank_index >= (n_total - 2L)) return(NULL)

  half_window <- max(5L, as.integer((wave_local_window - 1L) / 2L))
  left_lo  <- max(2L, rank_index - half_window)
  left_hi  <- max(left_lo, rank_index - 1L)
  right_lo <- min(n_total - 1L, rank_index + 1L)
  right_hi <- min(n_total - 1L, rank_index + half_window)
  if (left_hi <= left_lo || right_hi <= right_lo) return(NULL)

  iod_left  <- curve_df$smoothed_log_iod_nb[left_lo:left_hi]
  iod_right <- curve_df$smoothed_log_iod_nb[right_lo:right_hi]
  cv2_left  <- curve_df$smoothed_log_cv2_nb[left_lo:left_hi]
  cv2_right <- curve_df$smoothed_log_cv2_nb[right_lo:right_hi]

  iod_slope_term <- mean(diff(iod_left), na.rm = TRUE) - mean(diff(iod_right), na.rm = TRUE)
  iod_amp_term   <- sd(iod_left, na.rm = TRUE) - sd(iod_right, na.rm = TRUE)
  cv2_amp_term   <- sd(cv2_right, na.rm = TRUE) - sd(cv2_left, na.rm = TRUE)
  cv2_slope_term <- mean(diff(cv2_right), na.rm = TRUE) - mean(diff(cv2_left), na.rm = TRUE)

  wave_score <- wave_w_iod_slope * iod_slope_term +
    wave_w_iod_amp * iod_amp_term +
    wave_w_cv2_amp * cv2_amp_term +
    wave_w_cv2_slope * cv2_slope_term

  data.frame(
    rank_index = as.integer(rank_index),
    wave_score = as.numeric(wave_score),
    iod_slope_term = as.numeric(iod_slope_term),
    iod_amp_term = as.numeric(iod_amp_term),
    cv2_amp_term = as.numeric(cv2_amp_term),
    cv2_slope_term = as.numeric(cv2_slope_term),
    stringsAsFactors = FALSE
  )
}

resolve_matched_curvature_wave_cutoff <- function(loading_tbl,
                                                  mean_col = "baseMean",
                                                  dispersion_col = "dispGeneEst",
                                                  fixed_top_n = evs_fixed_top_n) {
  sorted_values_desc <- loading_tbl$pc1_loading_abs
  fallback <- resolve_top_n_cutoff(sorted_values_desc, top_n = fixed_top_n)
  n_total <- nrow(loading_tbl)

  if (isTRUE(evs_use_manual_rank_first) && is.finite(evs_fixed_rank_override) && !is.na(evs_fixed_rank_override)) {
    manual_rank <- max(1L, min(as.integer(evs_fixed_rank_override), n_total))
    return(list(
      top_n_actual = manual_rank,
      cutoff_value = loading_tbl$pc1_loading_abs[manual_rank],
      cutoff_quantile = 1 - (manual_rank / n_total),
      n_total = n_total,
      method = "manual_fixed_rank",
      curvature_strength = NA_real_,
      curve_df = NULL,
      max_rank_allowed = as.integer(n_total),
      variance_measure = "IOD_CV2_matched_curvature",
      candidate_table = data.frame(candidate_id = "manual_fixed_rank", rank_index = manual_rank, selected = TRUE, selected_reason = "manual_fixed_rank", candidate_source = "manual_fixed_rank", stringsAsFactors = FALSE),
      quantile_screen_table = data.frame(),
      matched_interval_table = data.frame(),
      selected_reason = "manual_fixed_rank"
    ))
  }

  curve_df <- smooth_ranked_regime_curves(loading_tbl, mean_col, dispersion_col)
  if (is.null(curve_df) || !nrow(curve_df)) {
    fallback$curve_df <- curve_df
    fallback$curvature_strength <- NA_real_
    fallback$method <- "fixed_top_n_fallback"
    fallback$candidate_table <- data.frame(candidate_id = "fallback_top_n", rank_index = fallback$top_n_actual, selected = TRUE, selected_reason = "fallback_top_n", candidate_source = "fallback", stringsAsFactors = FALSE)
    fallback$quantile_screen_table <- data.frame()
    fallback$matched_interval_table <- data.frame()
    return(fallback)
  }

  iod_peaks <- find_local_metric_curvature_peaks(curve_df, metric = "iod")
  cv2_peaks <- find_local_metric_curvature_peaks(curve_df, metric = "cv2")
  matched_intervals <- match_curvature_peak_intervals(iod_peaks, cv2_peaks)

  if (is.null(matched_intervals) || !nrow(matched_intervals)) {
    fallback$curve_df <- curve_df
    fallback$curvature_strength <- NA_real_
    fallback$method <- "fixed_top_n_fallback"
    fallback$candidate_table <- data.frame(candidate_id = "fallback_top_n", rank_index = fallback$top_n_actual, selected = TRUE, selected_reason = "fallback_top_n", candidate_source = "fallback", stringsAsFactors = FALSE)
    fallback$quantile_screen_table <- dplyr::bind_rows(iod_peaks, cv2_peaks)
    fallback$matched_interval_table <- matched_intervals
    return(fallback)
  }

  candidate_rows <- list()
  kk <- 1L
  for (ii in seq_len(nrow(matched_intervals))) {
    lo <- max(2L, matched_intervals$interval_center[ii] - evs_candidate_interval_halfwidth)
    hi <- min(n_total - 2L, matched_intervals$interval_center[ii] + evs_candidate_interval_halfwidth)
    for (r in seq.int(lo, hi, by = 1L)) {
      sc <- score_wave_candidate_rank(curve_df, r)
      if (!is.null(sc)) {
        sc$interval_id <- ii
        sc$joint_strength <- matched_intervals$joint_strength[ii]
        sc$candidate_source <- "matched_curvature_interval"
        candidate_rows[[kk]] <- sc
        kk <- kk + 1L
      }
    }
  }

  cand_tbl <- if (length(candidate_rows)) dplyr::bind_rows(candidate_rows) else data.frame()
  if (!nrow(cand_tbl)) {
    fallback$curve_df <- curve_df
    fallback$curvature_strength <- NA_real_
    fallback$method <- "fixed_top_n_fallback"
    fallback$candidate_table <- data.frame(candidate_id = "fallback_top_n", rank_index = fallback$top_n_actual, selected = TRUE, selected_reason = "fallback_top_n", candidate_source = "fallback", stringsAsFactors = FALSE)
    fallback$quantile_screen_table <- dplyr::bind_rows(iod_peaks, cv2_peaks)
    fallback$matched_interval_table <- matched_intervals
    return(fallback)
  }

  cand_tbl <- cand_tbl[order(-cand_tbl$wave_score, -cand_tbl$joint_strength, cand_tbl$rank_index), , drop = FALSE]
  cand_tbl <- cand_tbl[!duplicated(cand_tbl$rank_index), , drop = FALSE]
  cand_tbl$candidate_id <- paste0("candidate_", seq_len(nrow(cand_tbl)))
  cand_tbl$cutoff_value <- loading_tbl$pc1_loading_abs[cand_tbl$rank_index]
  cand_tbl$cutoff_quantile <- 1 - (cand_tbl$rank_index / n_total)
  cand_tbl$selected <- FALSE
  cand_tbl$selected_reason <- "candidate_only"
  best_idx <- 1L
  cand_tbl$selected[best_idx] <- TRUE
  cand_tbl$selected_reason[best_idx] <- "matched_curvature_wave_max_score"
  cand_tbl$curvature_strength <- cand_tbl$joint_strength
  cand_tbl$prominence <- cand_tbl$joint_strength
  cand_tbl$max_rank_allowed <- as.integer(n_total)

  list(
    top_n_actual = as.integer(cand_tbl$rank_index[best_idx]),
    cutoff_value = as.numeric(cand_tbl$cutoff_value[best_idx]),
    cutoff_quantile = as.numeric(cand_tbl$cutoff_quantile[best_idx]),
    n_total = n_total,
    method = "matched_curvature_wave",
    curvature_strength = as.numeric(cand_tbl$curvature_strength[best_idx]),
    curve_df = curve_df,
    max_rank_allowed = as.integer(n_total),
    variance_measure = "IOD_CV2_matched_curvature",
    candidate_table = cand_tbl,
    quantile_screen_table = dplyr::bind_rows(iod_peaks, cv2_peaks),
    matched_interval_table = matched_intervals,
    selected_reason = "matched_curvature_wave_max_score"
  )
}

# -----------------------------------------------------------------------------
# LOCAL FOURIER SHARED-CUTOFF METHOD
# -----------------------------------------------------------------------------

build_ranked_fourier_metric_table <- function(loading_tbl,
                                              mean_col = "baseMean",
                                              dispersion_col = "dispGeneEst") {
  needed <- c("feature_id", "rank", "pc1_loading_abs", mean_col, dispersion_col)
  if (!all(needed %in% colnames(loading_tbl))) return(NULL)

  df <- loading_tbl[, needed, drop = FALSE]
  colnames(df)[colnames(df) == mean_col] <- "mean_value"
  colnames(df)[colnames(df) == dispersion_col] <- "dispersion_value"

  df <- df[
    is.finite(df$mean_value) & !is.na(df$mean_value) & df$mean_value > 0 &
      is.finite(df$dispersion_value) & !is.na(df$dispersion_value) & df$dispersion_value > 0,
    , drop = FALSE
  ]
  if (nrow(df) < 25) return(NULL)

  df <- df[order(df$rank), , drop = FALSE]
  mu <- pmax(df$mean_value, 1e-12)
  alpha <- pmax(df$dispersion_value, 1e-12)

  df$iod_nb <- 1 + alpha * mu
  df$cv2_nb <- (1 / mu) + alpha
  df$log_iod_nb <- log10(df$iod_nb)
  df$log_cv2_nb <- log10(df$cv2_nb)
  df
}

build_percentile_windows_fourier <- function(n_total,
                                             step = fourier_percentile_step,
                                             window_fraction = fourier_window_fraction,
                                             min_rank = fourier_min_rank,
                                             max_rank_frac = fourier_max_rank_frac) {
  pct_grid <- seq(step, 1, by = step)
  center_ranks <- pmax(1L, pmin(n_total, round(pct_grid * n_total)))

  max_rank <- floor(n_total * max_rank_frac)
  keep <- center_ranks >= min_rank & center_ranks <= max_rank

  pct_grid <- pct_grid[keep]
  center_ranks <- center_ranks[keep]

  half_window <- max(5L, round((window_fraction * n_total) / 2))

  data.frame(
    percentile = pct_grid,
    center_rank = center_ranks,
    lo_rank = pmax(1L, center_ranks - half_window),
    hi_rank = pmin(n_total, center_ranks + half_window),
    stringsAsFactors = FALSE
  )
}

build_local_fourier_design <- function(x, n_harmonics = fourier_harmonics) {
  x <- as.numeric(x)
  x01 <- (x - min(x)) / max(1e-12, (max(x) - min(x)))
  out <- data.frame(x = x01)
  for (k in seq_len(n_harmonics)) {
    out[[paste0("sin_", k)]] <- sin(2 * pi * k * x01)
    out[[paste0("cos_", k)]] <- cos(2 * pi * k * x01)
  }
  out
}

fit_local_fourier <- function(rank_vec, y_vec, n_harmonics = fourier_harmonics) {
  rank_vec <- as.numeric(rank_vec)
  y_vec <- as.numeric(y_vec)

  keep <- is.finite(rank_vec) & !is.na(rank_vec) & is.finite(y_vec) & !is.na(y_vec)
  rank_vec <- rank_vec[keep]
  y_vec <- y_vec[keep]

  if (length(rank_vec) < (2 * n_harmonics + 5L)) return(NULL)

  y_sd <- suppressWarnings(stats::sd(y_vec, na.rm = TRUE))
  if (!is.finite(y_sd) || is.na(y_sd) || y_sd == 0) return(NULL)

  dd <- build_local_fourier_design(rank_vec, n_harmonics = n_harmonics)
  dd$y <- y_vec

  rhs <- paste(colnames(dd)[colnames(dd) != "y"], collapse = " + ")
  fm <- stats::as.formula(paste("y ~", rhs))

  fit <- tryCatch(stats::lm(fm, data = dd), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  fitted_y <- as.numeric(stats::predict(fit, newdata = dd))
  residual_y <- dd$y - fitted_y

  local_amplitude <- 0.5 * (max(fitted_y, na.rm = TRUE) - min(fitted_y, na.rm = TRUE))
  peak_idx <- which.max(fitted_y)
  trough_idx <- which.min(fitted_y)
  center_rank <- mean(c(rank_vec[peak_idx], rank_vec[trough_idx]))

  list(
    fit = fit,
    fitted_y = fitted_y,
    residual_y = residual_y,
    local_amplitude = local_amplitude,
    peak_rank = rank_vec[peak_idx],
    trough_rank = rank_vec[trough_idx],
    center_rank = center_rank,
    residual_sd = stats::sd(residual_y, na.rm = TRUE)
  )
}

summarize_local_fourier_window <- function(metric_df, lo_rank, hi_rank) {
  sub <- metric_df[metric_df$rank >= lo_rank & metric_df$rank <= hi_rank, , drop = FALSE]
  if (nrow(sub) < 9L) {
    return(data.frame(
      iod_amplitude = NA_real_,
      cv2_amplitude = NA_real_,
      iod_center = NA_real_,
      cv2_center = NA_real_,
      iod_residual_sd = NA_real_,
      cv2_residual_sd = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  iod_fit <- fit_local_fourier(sub$rank, sub$log_iod_nb)
  cv2_fit <- fit_local_fourier(sub$rank, sub$log_cv2_nb)
  if (is.null(iod_fit) || is.null(cv2_fit)) {
    return(data.frame(
      iod_amplitude = NA_real_,
      cv2_amplitude = NA_real_,
      iod_center = NA_real_,
      cv2_center = NA_real_,
      iod_residual_sd = NA_real_,
      cv2_residual_sd = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    iod_amplitude = iod_fit$local_amplitude,
    cv2_amplitude = cv2_fit$local_amplitude,
    iod_center = iod_fit$center_rank,
    cv2_center = cv2_fit$center_rank,
    iod_residual_sd = iod_fit$residual_sd,
    cv2_residual_sd = cv2_fit$residual_sd,
    stringsAsFactors = FALSE
  )
}

build_local_fourier_wave_map <- function(loading_tbl,
                                         mean_col = "baseMean",
                                         dispersion_col = "dispGeneEst") {
  metric_df <- build_ranked_fourier_metric_table(
    loading_tbl = loading_tbl,
    mean_col = mean_col,
    dispersion_col = dispersion_col
  )
  if (is.null(metric_df) || !nrow(metric_df)) return(NULL)

  windows <- build_percentile_windows_fourier(n_total = nrow(metric_df))
  if (!nrow(windows)) return(NULL)

  rows <- lapply(seq_len(nrow(windows)), function(i) {
    ww <- windows[i, , drop = FALSE]
    ss <- summarize_local_fourier_window(metric_df, ww$lo_rank, ww$hi_rank)
    cbind(ww, ss, stringsAsFactors = FALSE)
  })

  wave_map <- dplyr::bind_rows(rows)

  iod_amp_scaled <- if (all(is.na(wave_map$iod_amplitude))) {
    rep(NA_real_, nrow(wave_map))
  } else {
    scales::rescale(wave_map$iod_amplitude, to = c(0, 1), from = range(wave_map$iod_amplitude, na.rm = TRUE))
  }

  cv2_amp_scaled <- if (all(is.na(wave_map$cv2_amplitude))) {
    rep(NA_real_, nrow(wave_map))
  } else {
    scales::rescale(wave_map$cv2_amplitude, to = c(0, 1), from = range(wave_map$cv2_amplitude, na.rm = TRUE))
  }

  wave_map$center_distance <- abs(wave_map$iod_center - wave_map$cv2_center)
  wave_map$center_agreement <- 1 / (1 + wave_map$center_distance)
  wave_map$local_fourier_score <- (
    fourier_score_weight_iod_amp * iod_amp_scaled +
      fourier_score_weight_cv2_amp * cv2_amp_scaled +
      fourier_score_weight_center_agreement * wave_map$center_agreement
  )

  list(
    metric_df = metric_df,
    wave_map = wave_map
  )
}

combine_treatment_control_fourier_maps <- function(trt_wave_obj, ctrl_wave_obj) {
  if (is.null(trt_wave_obj) || is.null(ctrl_wave_obj)) return(NULL)
  trt_map <- trt_wave_obj$wave_map
  ctrl_map <- ctrl_wave_obj$wave_map
  if (is.null(trt_map) || is.null(ctrl_map) || !nrow(trt_map) || !nrow(ctrl_map)) return(NULL)

  keep_cols <- c("percentile", "center_rank", "iod_amplitude", "cv2_amplitude", "iod_center", "cv2_center", "local_fourier_score")
  trt_map2 <- trt_map[, keep_cols, drop = FALSE]
  ctrl_map2 <- ctrl_map[, keep_cols, drop = FALSE]

  names(trt_map2)[names(trt_map2) != "percentile"] <- paste0(names(trt_map2)[names(trt_map2) != "percentile"], "_trt")
  names(ctrl_map2)[names(ctrl_map2) != "percentile"] <- paste0(names(ctrl_map2)[names(ctrl_map2) != "percentile"], "_ctrl")

  out <- dplyr::inner_join(trt_map2, ctrl_map2, by = "percentile")
  if (!nrow(out)) return(NULL)

  out$combined_center_rank <- round((out$center_rank_trt + out$center_rank_ctrl) / 2)
  out$combined_iod_amplitude <- out$iod_amplitude_trt + out$iod_amplitude_ctrl
  out$combined_cv2_amplitude <- out$cv2_amplitude_trt + out$cv2_amplitude_ctrl

  out$combined_center_distance <- abs(
    ((out$iod_center_trt + out$iod_center_ctrl) / 2) -
      ((out$cv2_center_trt + out$cv2_center_ctrl) / 2)
  )
  out$combined_center_agreement <- 1 / (1 + out$combined_center_distance)

  iod_amp_scaled <- if (all(is.na(out$combined_iod_amplitude))) {
    rep(NA_real_, nrow(out))
  } else {
    scales::rescale(out$combined_iod_amplitude, to = c(0, 1), from = range(out$combined_iod_amplitude, na.rm = TRUE))
  }

  cv2_amp_scaled <- if (all(is.na(out$combined_cv2_amplitude))) {
    rep(NA_real_, nrow(out))
  } else {
    scales::rescale(out$combined_cv2_amplitude, to = c(0, 1), from = range(out$combined_cv2_amplitude, na.rm = TRUE))
  }

  out$combined_fourier_score <- (
    fourier_score_weight_iod_amp * iod_amp_scaled +
      fourier_score_weight_cv2_amp * cv2_amp_scaled +
      fourier_score_weight_center_agreement * out$combined_center_agreement
  )

  out
}

compute_regime_difference_curve <- function(combined_wave_df) {
  df <- as.data.frame(combined_wave_df, stringsAsFactors = FALSE)
  required_cols <- c(
    "percentile", "combined_center_rank",
    "combined_iod_amplitude", "combined_cv2_amplitude"
  )
  assert_required_columns(df, required_cols, object_name = "combined_wave_df")

  df$regime_difference <- df$combined_iod_amplitude - df$combined_cv2_amplitude
  df$regime_direction <- ifelse(
    is.na(df$regime_difference),
    NA_character_,
    ifelse(df$regime_difference > 0, "IOD_dominant",
           ifelse(df$regime_difference < 0, "CV2_dominant", "balanced"))
  )
  df
}

find_crossing_intervals <- function(diff_df,
                                    min_percentile = crossing_min_percentile,
                                    max_percentile = crossing_max_percentile) {
  df <- as.data.frame(diff_df, stringsAsFactors = FALSE)
  df <- df[order(df$percentile), , drop = FALSE]
  keep <- is.finite(df$percentile) & !is.na(df$percentile) &
    is.finite(df$regime_difference) & !is.na(df$regime_difference) &
    df$percentile >= min_percentile & df$percentile <= max_percentile
  df <- df[keep, , drop = FALSE]
  if (nrow(df) < 2L) return(data.frame())

  rows <- list()
  kk <- 1L
  for (i in seq_len(nrow(df) - 1L)) {
    y1 <- df$regime_difference[i]
    y2 <- df$regime_difference[i + 1L]
    if (!is.finite(y1) || !is.finite(y2)) next
    crossed <- (y1 == 0) || (y2 == 0) || ((y1 > 0) && (y2 < 0)) || ((y1 < 0) && (y2 > 0))
    if (!crossed) next

    x1 <- df$percentile[i]
    x2 <- df$percentile[i + 1L]
    r1 <- df$combined_center_rank[i]
    r2 <- df$combined_center_rank[i + 1L]

    if (isTRUE(all.equal(y1, y2))) {
      crossing_percentile <- mean(c(x1, x2))
      crossing_rank <- round(mean(c(r1, r2)))
    } else {
      crossing_percentile <- x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
      crossing_rank <- round(r1 + (0 - y1) * (r2 - r1) / (y2 - y1))
    }

    rows[[kk]] <- data.frame(
      crossing_id = paste0("crossing_", kk),
      idx_left = i,
      idx_right = i + 1L,
      percentile_left = x1,
      percentile_right = x2,
      rank_left = r1,
      rank_right = r2,
      regime_difference_left = y1,
      regime_difference_right = y2,
      crossing_percentile = crossing_percentile,
      crossing_rank = as.integer(crossing_rank),
      stringsAsFactors = FALSE
    )
    kk <- kk + 1L
  }
  if (!length(rows)) return(data.frame())
  dplyr::bind_rows(rows)
}

label_stable_crossings <- function(diff_df,
                                   crossing_tbl,
                                   stability_window_n = crossing_stability_window_n) {
  if (is.null(crossing_tbl) || !nrow(crossing_tbl)) return(data.frame())
  df <- as.data.frame(diff_df, stringsAsFactors = FALSE)
  out <- crossing_tbl
  out$left_window_positive_frac <- NA_real_
  out$right_window_negative_frac <- NA_real_
  out$stable_crossing <- FALSE

  for (i in seq_len(nrow(out))) {
    il <- out$idx_left[i]
    ir <- out$idx_right[i]

    left_idx <- seq.int(max(1L, il - stability_window_n + 1L), il, by = 1L)
    right_idx <- seq.int(ir, min(nrow(df), ir + stability_window_n - 1L), by = 1L)

    left_vals <- df$regime_difference[left_idx]
    right_vals <- df$regime_difference[right_idx]

    left_pos_frac <- mean(left_vals > 0, na.rm = TRUE)
    right_neg_frac <- mean(right_vals < 0, na.rm = TRUE)

    out$left_window_positive_frac[i] <- left_pos_frac
    out$right_window_negative_frac[i] <- right_neg_frac
    out$stable_crossing[i] <- isTRUE(
      is.finite(left_pos_frac) && is.finite(right_neg_frac) &&
        left_pos_frac >= 0.67 && right_neg_frac >= 0.67
    )
  }

  out
}

select_regime_shift_crossing <- function(combined_wave_df) {
  diff_df <- compute_regime_difference_curve(combined_wave_df)
  crossing_tbl <- find_crossing_intervals(diff_df)
  crossing_tbl <- label_stable_crossings(diff_df, crossing_tbl)

  if (!nrow(crossing_tbl)) {
    return(list(
      diff_df = diff_df,
      crossing_table = crossing_tbl,
      selected_crossing = NULL,
      selected_reason = "no_crossings_found"
    ))
  }

  stable_tbl <- crossing_tbl[crossing_tbl$stable_crossing, , drop = FALSE]
  if (nrow(stable_tbl)) {
    selected <- stable_tbl[order(stable_tbl$crossing_percentile), , drop = FALSE][1, , drop = FALSE]
    reason <- "first_stable_crossing"
  } else {
    selected <- crossing_tbl[order(crossing_tbl$crossing_percentile), , drop = FALSE][1, , drop = FALSE]
    reason <- "first_crossing_fallback"
  }

  list(
    diff_df = diff_df,
    crossing_table = crossing_tbl,
    selected_crossing = selected,
    selected_reason = reason
  )
}

# -----------------------------------------------------------------------------
# PATCHED SHARED COMBINED RANK TABLE + SHARED CUTOFF
# -----------------------------------------------------------------------------

build_shared_combined_loading_table <- function(fit_trt_loading_tbl,
                                                fit_ctrl_loading_tbl) {
  trt_tbl <- as.data.frame(fit_trt_loading_tbl, stringsAsFactors = FALSE)
  ctrl_tbl <- as.data.frame(fit_ctrl_loading_tbl, stringsAsFactors = FALSE)

  assert_required_columns(
    trt_tbl,
    c("feature_id", "pc1_loading", "pc1_loading_abs", "rank"),
    object_name = "fit_trt_loading_tbl"
  )
  assert_required_columns(
    ctrl_tbl,
    c("feature_id", "pc1_loading", "pc1_loading_abs", "rank"),
    object_name = "fit_ctrl_loading_tbl"
  )

  trt_use <- trt_tbl[, c("feature_id", "pc1_loading", "pc1_loading_abs", "rank"), drop = FALSE]
  ctrl_use <- ctrl_tbl[, c("feature_id", "pc1_loading", "pc1_loading_abs", "rank"), drop = FALSE]

  names(trt_use) <- c("feature_id", "pc1_loading_trt", "pc1_loading_abs_trt", "rank_trt")
  names(ctrl_use) <- c("feature_id", "pc1_loading_ctrl", "pc1_loading_abs_ctrl", "rank_ctrl")

  merged <- dplyr::full_join(trt_use, ctrl_use, by = "feature_id")

  merged$pc1_loading_abs_trt[is.na(merged$pc1_loading_abs_trt)] <- 0
  merged$pc1_loading_abs_ctrl[is.na(merged$pc1_loading_abs_ctrl)] <- 0

  merged$combined_loading <- pmax(
    merged$pc1_loading_abs_trt,
    merged$pc1_loading_abs_ctrl,
    na.rm = TRUE
  )

  merged <- merged[order(-merged$combined_loading, merged$feature_id), , drop = FALSE]
  merged$combined_rank <- seq_len(nrow(merged))
  rownames(merged) <- NULL

  merged
}

resolve_combined_fourier_cutoff <- function(fit_trt_loading_tbl,
                                            fit_ctrl_loading_tbl,
                                            fixed_top_n = evs_fixed_top_n) {
  shared_combined_tbl <- build_shared_combined_loading_table(
    fit_trt_loading_tbl  = fit_trt_loading_tbl,
    fit_ctrl_loading_tbl = fit_ctrl_loading_tbl
  )

  n_total <- nrow(shared_combined_tbl)

  fallback <- resolve_top_n_cutoff(
    sorted_values_desc = shared_combined_tbl$combined_loading,
    top_n = fixed_top_n
  )

  trt_wave_obj <- build_local_fourier_wave_map(fit_trt_loading_tbl)
  ctrl_wave_obj <- build_local_fourier_wave_map(fit_ctrl_loading_tbl)

  if (is.null(trt_wave_obj) || is.null(ctrl_wave_obj)) {
    fallback$method <- "fixed_top_n_fallback"
    fallback$trt_wave_obj <- trt_wave_obj
    fallback$ctrl_wave_obj <- ctrl_wave_obj
    fallback$shared_combined_tbl <- shared_combined_tbl
    fallback$combined_wave_map <- data.frame()
    fallback$combined_interval_table <- data.frame()
    fallback$crossing_table <- data.frame()
    fallback$selected_crossing <- NULL
    fallback$selected_reason <- "combined_wave_map_missing"
    return(fallback)
  }

  combined_wave_df <- combine_treatment_control_fourier_maps(trt_wave_obj, ctrl_wave_obj)

  if (is.null(combined_wave_df) || !nrow(combined_wave_df)) {
    fallback$method <- "fixed_top_n_fallback"
    fallback$trt_wave_obj <- trt_wave_obj
    fallback$ctrl_wave_obj <- ctrl_wave_obj
    fallback$shared_combined_tbl <- shared_combined_tbl
    fallback$combined_wave_map <- data.frame()
    fallback$combined_interval_table <- data.frame()
    fallback$crossing_table <- data.frame()
    fallback$selected_crossing <- NULL
    fallback$selected_reason <- "combined_wave_map_missing"
    return(fallback)
  }

  crossing_info <- select_regime_shift_crossing(combined_wave_df)
  selected_crossing <- crossing_info$selected_crossing

  if (is.null(selected_crossing) || !nrow(selected_crossing)) {
    fallback$method <- "fixed_top_n_fallback"
    fallback$trt_wave_obj <- trt_wave_obj
    fallback$ctrl_wave_obj <- ctrl_wave_obj
    fallback$shared_combined_tbl <- shared_combined_tbl
    fallback$combined_wave_map <- crossing_info$diff_df
    fallback$combined_interval_table <- data.frame()
    fallback$crossing_table <- crossing_info$crossing_table
    fallback$selected_crossing <- NULL
    fallback$selected_reason <- crossing_info$selected_reason
    return(fallback)
  }

  selected_rank <- as.integer(selected_crossing$crossing_rank[1])
  selected_rank <- min(max(1L, selected_rank), n_total)

  cutoff_value <- shared_combined_tbl$combined_loading[selected_rank]
  cutoff_quantile <- 1 - (selected_rank / n_total)

  candidate_table <- data.frame(
    candidate_id = selected_crossing$crossing_id[1],
    rank_index = selected_rank,
    percentile = selected_crossing$crossing_percentile[1],
    cutoff_value = cutoff_value,
    cutoff_quantile = cutoff_quantile,
    selected = TRUE,
    selected_reason = crossing_info$selected_reason,
    stringsAsFactors = FALSE
  )

  list(
    top_n_actual = selected_rank,
    cutoff_value = cutoff_value,
    cutoff_quantile = cutoff_quantile,
    n_total = n_total,
    method = "first_stable_crossing",
    trt_wave_obj = trt_wave_obj,
    ctrl_wave_obj = ctrl_wave_obj,
    shared_combined_tbl = shared_combined_tbl,
    combined_wave_map = crossing_info$diff_df,
    combined_interval_table = data.frame(),
    crossing_table = crossing_info$crossing_table,
    selected_crossing = selected_crossing,
    candidate_table = candidate_table,
    selected_reason = crossing_info$selected_reason
  )
}


select_wave_backup_cutoff <- function(wave_obj, group_label = "group") {
  if (is.null(wave_obj) || is.null(wave_obj$wave_map) || !nrow(wave_obj$wave_map)) {
    return(list(
      diff_df = data.frame(),
      crossing_table = data.frame(),
      selected_crossing = NULL,
      selected_reason = "wave_map_missing"
    ))
  }

  df <- as.data.frame(wave_obj$wave_map, stringsAsFactors = FALSE)
  assert_required_columns(
    df,
    c("percentile", "center_rank", "iod_amplitude", "cv2_amplitude"),
    object_name = paste0(group_label, " wave_map")
  )

  diff_df <- data.frame(
    percentile = as.numeric(df$percentile),
    combined_center_rank = as.integer(round(df$center_rank)),
    regime_difference = as.numeric(df$iod_amplitude - df$cv2_amplitude),
    stringsAsFactors = FALSE
  )

  crossing_tbl <- find_crossing_intervals(diff_df)
  crossing_tbl <- label_stable_crossings(diff_df, crossing_tbl)

  if (nrow(crossing_tbl)) {
    stable_tbl <- crossing_tbl[crossing_tbl$stable_crossing, , drop = FALSE]

    if (nrow(stable_tbl)) {
      selected <- stable_tbl[order(stable_tbl$crossing_percentile), , drop = FALSE][1, , drop = FALSE]
      reason <- paste0(group_label, "_first_stable_crossing")
    } else {
      selected <- crossing_tbl[order(crossing_tbl$crossing_percentile), , drop = FALSE][1, , drop = FALSE]
      reason <- paste0(group_label, "_first_crossing_fallback")
    }

    return(list(
      diff_df = diff_df,
      crossing_table = crossing_tbl,
      selected_crossing = selected,
      selected_reason = reason
    ))
  }

  finite_idx <- which(is.finite(diff_df$regime_difference) & !is.na(diff_df$regime_difference))
  if (!length(finite_idx)) {
    return(list(
      diff_df = diff_df,
      crossing_table = data.frame(),
      selected_crossing = NULL,
      selected_reason = paste0(group_label, "_no_finite_difference")
    ))
  }

  best_idx <- finite_idx[which.min(abs(diff_df$regime_difference[finite_idx]))]

  selected <- data.frame(
    crossing_id = paste0(group_label, "_nearest_zero"),
    idx_left = best_idx,
    idx_right = best_idx,
    percentile_left = diff_df$percentile[best_idx],
    percentile_right = diff_df$percentile[best_idx],
    rank_left = diff_df$combined_center_rank[best_idx],
    rank_right = diff_df$combined_center_rank[best_idx],
    regime_difference_left = diff_df$regime_difference[best_idx],
    regime_difference_right = diff_df$regime_difference[best_idx],
    crossing_percentile = diff_df$percentile[best_idx],
    crossing_rank = diff_df$combined_center_rank[best_idx],
    left_window_positive_frac = NA_real_,
    right_window_negative_frac = NA_real_,
    stable_crossing = FALSE,
    stringsAsFactors = FALSE
  )

  list(
    diff_df = diff_df,
    crossing_table = data.frame(),
    selected_crossing = selected,
    selected_reason = paste0(group_label, "_nearest_zero_fallback")
  )
}

build_group_backup_crossings <- function(combined_cutoff_info) {
  trt_backup  <- select_wave_backup_cutoff(combined_cutoff_info$trt_wave_obj,  group_label = "treatment")
  ctrl_backup <- select_wave_backup_cutoff(combined_cutoff_info$ctrl_wave_obj, group_label = "control")

  combined_cutoff_info$trt_backup_crossing  <- trt_backup$selected_crossing
  combined_cutoff_info$trt_backup_reason    <- trt_backup$selected_reason
  combined_cutoff_info$trt_backup_table     <- trt_backup$crossing_table

  combined_cutoff_info$ctrl_backup_crossing <- ctrl_backup$selected_crossing
  combined_cutoff_info$ctrl_backup_reason   <- ctrl_backup$selected_reason
  combined_cutoff_info$ctrl_backup_table    <- ctrl_backup$crossing_table

  combined_cutoff_info
}

build_backup_crossing_summary_table <- function(combined_cutoff_info, comparison_name) {
  add_row <- function(sc, label, reason) {
    if (is.null(sc) || !nrow(sc)) {
      data.frame(
        comparison_name = comparison_name,
        cutoff_type = label,
        selected_reason = reason,
        crossing_percentile = NA_real_,
        crossing_rank = NA_integer_,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        comparison_name = comparison_name,
        cutoff_type = label,
        selected_reason = reason,
        crossing_percentile = sc$crossing_percentile[1],
        crossing_rank = sc$crossing_rank[1],
        stringsAsFactors = FALSE
      )
    }
  }

  dplyr::bind_rows(
    add_row(combined_cutoff_info$selected_crossing, "shared_combined", combined_cutoff_info$selected_reason),
    add_row(combined_cutoff_info$trt_backup_crossing, "treatment_backup", combined_cutoff_info$trt_backup_reason),
    add_row(combined_cutoff_info$ctrl_backup_crossing, "control_backup", combined_cutoff_info$ctrl_backup_reason)
  )
}


safe_skewness <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) < 5) return(NA_real_)
  sx <- stats::sd(x)
  if (!is.finite(sx) || is.na(sx) || sx == 0) return(NA_real_)
  mean(((x - mean(x)) / sx)^3)
}

compute_raw_wald_normality_metrics <- function(z) {
  z <- as.numeric(z)
  z <- z[is.finite(z) & !is.na(z)]

  if (length(z) < 5) {
    return(data.frame(
      n = length(z),
      mean_abs = Inf,
      sd_dev = Inf,
      skew_abs = Inf,
      distortion = Inf,
      stringsAsFactors = FALSE
    ))
  }

  mu <- mean(z)
  sdv <- stats::sd(z)
  skw <- safe_skewness(z)

  mean_abs <- abs(mu)
  sd_dev <- abs(sdv - 1)
  skew_abs <- abs(skw)

  distortion <- raw_wald_score_weight_mean * mean_abs +
    raw_wald_score_weight_sd * sd_dev +
    raw_wald_score_weight_skew * skew_abs

  data.frame(
    n = length(z),
    mean_abs = mean_abs,
    sd_dev = sd_dev,
    skew_abs = skew_abs,
    distortion = distortion,
    stringsAsFactors = FALSE
  )
}

score_raw_wald_split <- function(rank_index, ordered_wald_df) {
  n_total <- nrow(ordered_wald_df)
  if (!is.finite(rank_index) || is.na(rank_index)) return(NULL)

  rank_index <- as.integer(rank_index)
  if (rank_index < raw_wald_search_min_subset_n) return(NULL)
  if ((n_total - rank_index) < raw_wald_search_min_subset_n) return(NULL)

  z_lead <- ordered_wald_df$wald_stat[seq_len(rank_index)]
  z_rem  <- ordered_wald_df$wald_stat[(rank_index + 1L):n_total]

  lead_metrics <- compute_raw_wald_normality_metrics(z_lead)
  rem_metrics  <- compute_raw_wald_normality_metrics(z_rem)

  if (!is.finite(lead_metrics$distortion[1]) || !is.finite(rem_metrics$distortion[1])) {
    return(NULL)
  }

  total_distortion <- lead_metrics$distortion[1] + rem_metrics$distortion[1]
  balance_penalty <- abs(lead_metrics$distortion[1] - rem_metrics$distortion[1])

  split_score <- -(total_distortion + raw_wald_search_balance_lambda * balance_penalty)

  data.frame(
    rank_index = rank_index,
    n_leading_edge = rank_index,
    n_remainder = n_total - rank_index,
    leading_mean_abs = lead_metrics$mean_abs[1],
    leading_sd_dev = lead_metrics$sd_dev[1],
    leading_skew_abs = lead_metrics$skew_abs[1],
    leading_distortion = lead_metrics$distortion[1],
    remainder_mean_abs = rem_metrics$mean_abs[1],
    remainder_sd_dev = rem_metrics$sd_dev[1],
    remainder_skew_abs = rem_metrics$skew_abs[1],
    remainder_distortion = rem_metrics$distortion[1],
    total_distortion = total_distortion,
    balance_penalty = balance_penalty,
    split_score = split_score,
    stringsAsFactors = FALSE
  )
}

build_raw_wald_rank_table <- function(full_results_df, shared_combined_tbl) {
  res_df <- as.data.frame(full_results_df, stringsAsFactors = FALSE)
  assert_required_columns(res_df, c("feature_id", "stat"), object_name = "full_results_df")
  assert_required_columns(shared_combined_tbl, c("feature_id", "combined_rank"), object_name = "shared_combined_tbl")

  res_df$feature_id <- as.character(res_df$feature_id)
  shared_combined_tbl$feature_id <- as.character(shared_combined_tbl$feature_id)

  out <- dplyr::left_join(
    shared_combined_tbl[, c("feature_id", "combined_rank"), drop = FALSE],
    res_df[, c("feature_id", "stat"), drop = FALSE],
    by = "feature_id"
  )

  colnames(out)[colnames(out) == "stat"] <- "wald_stat"
  out <- out[order(out$combined_rank), , drop = FALSE]
  out <- out[is.finite(out$wald_stat) & !is.na(out$wald_stat), , drop = FALSE]
  rownames(out) <- NULL
  out
}

run_raw_wald_bracket_search <- function(full_results_df,
                                        shared_combined_tbl,
                                        lower_rank,
                                        upper_rank,
                                        coarse_step = raw_wald_search_coarse_step,
                                        fine_radius = raw_wald_search_fine_radius) {
  ordered_wald_df <- build_raw_wald_rank_table(full_results_df, shared_combined_tbl)
  n_total <- nrow(ordered_wald_df)

  lower_rank <- max(raw_wald_search_min_subset_n, as.integer(lower_rank))
  upper_rank <- min(n_total - raw_wald_search_min_subset_n, as.integer(upper_rank))

  if (!is.finite(lower_rank) || !is.finite(upper_rank) || lower_rank >= upper_rank) {
    stop("Invalid raw-Wald bracket search interval.", call. = FALSE)
  }

  coarse_ranks <- seq.int(lower_rank, upper_rank, by = max(1L, as.integer(coarse_step)))
  if (tail(coarse_ranks, 1) != upper_rank) coarse_ranks <- c(coarse_ranks, upper_rank)

  coarse_rows <- lapply(coarse_ranks, function(rk) score_raw_wald_split(rk, ordered_wald_df))
  coarse_tbl <- dplyr::bind_rows(Filter(Negate(is.null), coarse_rows))
  if (!nrow(coarse_tbl)) stop("No valid coarse raw-Wald split candidates were produced.", call. = FALSE)

  best_coarse_rank <- coarse_tbl$rank_index[which.max(coarse_tbl$split_score)]

  fine_lo <- max(lower_rank, best_coarse_rank - as.integer(fine_radius))
  fine_hi <- min(upper_rank, best_coarse_rank + as.integer(fine_radius))
  fine_ranks <- seq.int(fine_lo, fine_hi, by = 1L)

  fine_rows <- lapply(fine_ranks, function(rk) score_raw_wald_split(rk, ordered_wald_df))
  fine_tbl <- dplyr::bind_rows(Filter(Negate(is.null), fine_rows))
  if (!nrow(fine_tbl)) stop("No valid fine raw-Wald split candidates were produced.", call. = FALSE)

  best_fine <- fine_tbl[which.max(fine_tbl$split_score), , drop = FALSE]

  list(
    ordered_wald_df = ordered_wald_df,
    coarse_table = coarse_tbl,
    fine_table = fine_tbl,
    selected_rank = as.integer(best_fine$rank_index[1]),
    selected_row = best_fine
  )
}

plot_fourier_wave_map_single <- function(wave_obj, comparison_name, group_label) {
  if (is.null(wave_obj) || is.null(wave_obj$wave_map) || !nrow(wave_obj$wave_map)) return(NULL)
  df <- wave_obj$wave_map

  backup_info <- select_wave_backup_cutoff(wave_obj, group_label = tolower(pretty_group_label(group_label)))
  sc <- backup_info$selected_crossing

  cutoff_pct  <- if (!is.null(sc) && nrow(sc)) sc$crossing_percentile[1] else NA_real_
  cutoff_rank <- if (!is.null(sc) && nrow(sc)) sc$crossing_rank[1] else NA_integer_

  line_col <- if (tolower(group_label) %in% c("treatment", "trt")) plot_palette$treatment else plot_palette$control

  p <- ggplot(df, aes(percentile)) +
    geom_line(aes(y = iod_amplitude, color = "IOD"), linewidth = 0.9) +
    geom_line(aes(y = cv2_amplitude, color = "CV²"), linewidth = 0.9) +
    scale_color_manual(values = c("IOD" = plot_palette$treatment, "CV²" = plot_palette$control)) +
    labs(
      title = paste0(pretty_group_label(group_label), " | local regime lines"),
      subtitle = compact_caption(
        "This panel includes the backup within-group cutoff. It is the first stable crossing of local IOD and local CV² for that group alone.",
        width = 88
      ),
      x = "Percentile center",
      y = "Local amplitude",
      color = NULL
    ) +
    manuscript_theme()

  if (is.finite(cutoff_pct)) {
    ymax <- max(c(df$iod_amplitude, df$cv2_amplitude), na.rm = TRUE)
    p <- p +
      geom_vline(
        xintercept = cutoff_pct,
        linetype = "dotted",
        linewidth = 0.9,
        colour = line_col
      ) +
      annotate(
        "label",
        x = cutoff_pct,
        y = ymax,
        label = paste0(
          pretty_group_label(group_label), " backup\n",
          "p = ", signif(cutoff_pct, 4), "\n",
          "rank = ", cutoff_rank
        ),
        fill = "white",
        colour = line_col,
        size = 2.8,
        label.size = 0.15,
        vjust = -0.5
      )
  }

  p
}

plot_combined_fourier_wave_map <- function(combined_cutoff_info, comparison_name) {
  df <- combined_cutoff_info$combined_wave_map
  if (is.null(df) || !nrow(df)) return(NULL)

  sc <- combined_cutoff_info$selected_crossing
  crossing_pct <- if (!is.null(sc) && nrow(sc)) sc$crossing_percentile[1] else NA_real_
  crossing_rank <- if (!is.null(sc) && nrow(sc)) sc$crossing_rank[1] else NA_integer_

  trt_sc <- combined_cutoff_info$trt_backup_crossing
  trt_pct <- if (!is.null(trt_sc) && nrow(trt_sc)) trt_sc$crossing_percentile[1] else NA_real_
  trt_rank <- if (!is.null(trt_sc) && nrow(trt_sc)) trt_sc$crossing_rank[1] else NA_integer_

  ctrl_sc <- combined_cutoff_info$ctrl_backup_crossing
  ctrl_pct <- if (!is.null(ctrl_sc) && nrow(ctrl_sc)) ctrl_sc$crossing_percentile[1] else NA_real_
  ctrl_rank <- if (!is.null(ctrl_sc) && nrow(ctrl_sc)) ctrl_sc$crossing_rank[1] else NA_integer_

  ymax <- max(c(df$combined_iod_amplitude, df$combined_cv2_amplitude), na.rm = TRUE)

  p <- ggplot(df, aes(percentile)) +
    geom_line(aes(y = combined_iod_amplitude, color = "Composite IOD"), linewidth = crossing_plot_line_width) +
    geom_line(aes(y = combined_cv2_amplitude, color = "Composite CV²"), linewidth = crossing_plot_line_width) +
    scale_color_manual(values = c("Composite IOD" = plot_palette$treatment, "Composite CV²" = plot_palette$control)) +
    labs(
      title = paste0(comparison_name, " | two-line regime crossing"),
      subtitle = compact_caption(
        "The dashed line is the selected shared comparison cutoff. The dotted blue and dotted black lines are the backup treatment-only and control-only crossings.",
        width = 90
      ),
      x = "Percentile center",
      y = "Composite local amplitude",
      color = NULL
    ) +
    manuscript_theme()

  if (is.finite(crossing_pct)) {
    p <- p +
      geom_vline(
        xintercept = crossing_pct,
        linetype = "dashed",
        linewidth = crossing_plot_vline_width,
        colour = plot_palette$threshold
      ) +
      annotate(
        "label",
        x = crossing_pct,
        y = ymax,
        label = paste0(
          "Shared cutoff\n",
          "p = ", signif(crossing_pct, 4), "\n",
          "rank = ", crossing_rank
        ),
        fill = "white",
        colour = plot_palette$threshold,
        size = 3.0,
        label.size = 0.15,
        vjust = -0.55
      )
  }

  if (is.finite(trt_pct)) {
    p <- p +
      geom_vline(
        xintercept = trt_pct,
        linetype = "dotted",
        linewidth = 0.85,
        colour = plot_palette$treatment
      ) +
      annotate(
        "label",
        x = trt_pct,
        y = ymax * 0.78,
        label = paste0(
          "Treatment backup\n",
          "p = ", signif(trt_pct, 4), "\n",
          "rank = ", trt_rank
        ),
        fill = "white",
        colour = plot_palette$treatment,
        size = 2.7,
        label.size = 0.15
      )
  }

  if (is.finite(ctrl_pct)) {
    p <- p +
      geom_vline(
        xintercept = ctrl_pct,
        linetype = "dotted",
        linewidth = 0.85,
        colour = plot_palette$control
      ) +
      annotate(
        "label",
        x = ctrl_pct,
        y = ymax * 0.58,
        label = paste0(
          "Control backup\n",
          "p = ", signif(ctrl_pct, 4), "\n",
          "rank = ", ctrl_rank
        ),
        fill = "white",
        colour = plot_palette$control,
        size = 2.7,
        label.size = 0.15
      )
  }

  p
}

plot_combined_fourier_score <- function(combined_cutoff_info, comparison_name) {
  df <- combined_cutoff_info$combined_wave_map
  cand_df <- combined_cutoff_info$candidate_table
  if (is.null(df) || !nrow(df) || !"combined_fourier_score" %in% names(df)) return(NULL)

  p <- ggplot(df, aes(percentile, combined_fourier_score)) +
    geom_col(width = 0.008, fill = plot_palette$threshold)

  if (!is.null(cand_df) && nrow(cand_df)) {
    chosen_pct <- cand_df$percentile[1]
    p <- p + geom_vline(xintercept = chosen_pct, color = plot_palette$threshold, linewidth = 0.9)
  }

  p +
    labs(
      title = paste0(comparison_name, " | descriptive score profile"),
      subtitle = compact_caption("The score profile is descriptive only. The EVS cutoff is determined by the first stable crossing of the two regime lines, not by the global score maximum.", width = 90),
      x = "Percentile center",
      y = "Score"
    ) +
    manuscript_theme()
}

plot_fourier_summary_panel <- function(combined_cutoff_info, comparison_name) {
  p_trt <- plot_fourier_wave_map_single(combined_cutoff_info$trt_wave_obj, comparison_name, "treatment")
  p_ctrl <- plot_fourier_wave_map_single(combined_cutoff_info$ctrl_wave_obj, comparison_name, "control")
  p_comb <- plot_combined_fourier_wave_map(combined_cutoff_info, comparison_name)
  p_score <- plot_combined_fourier_score(combined_cutoff_info, comparison_name)

  grobs <- list(p_trt, p_ctrl, p_comb, p_score)
  if (any(vapply(grobs, is.null, logical(1)))) return(NULL)

  arrangeGrob(
    grobs = grobs,
    ncol = 2,
    top = textGrob(
      paste0(comparison_name, " | Fourier summary panel"),
      gp = gpar(fontface = "bold", cex = 1.06)
    ),
    bottom = textGrob(
      "Top left: treatment local IOD/CV² with treatment backup cutoff. Top right: control local IOD/CV² with control backup cutoff. Bottom left: combined local IOD/CV² with shared cutoff and backup ranks. Bottom right: descriptive combined score profile.",
      gp = gpar(cex = 0.86)
    )
  )
}

plot_regime_difference_curve <- function(combined_cutoff_info, comparison_name) {
  df <- combined_cutoff_info$combined_wave_map
  if (is.null(df) || !nrow(df) || !"regime_difference" %in% names(df)) return(NULL)

  sc <- combined_cutoff_info$selected_crossing
  crossing_pct <- if (!is.null(sc) && nrow(sc)) sc$crossing_percentile[1] else NA_real_
  crossing_rank <- if (!is.null(sc) && nrow(sc)) sc$crossing_rank[1] else NA_integer_

  ggplot(df, aes(percentile, regime_difference)) +
    geom_hline(yintercept = 0, linewidth = crossing_plot_hline_width, colour = "grey50") +
    geom_line(linewidth = crossing_plot_line_width, colour = plot_palette$threshold) +
    geom_vline(xintercept = crossing_pct, linetype = "dashed", linewidth = crossing_plot_vline_width, colour = plot_palette$threshold) +
    geom_point(data = data.frame(percentile = crossing_pct, regime_difference = 0), aes(x = percentile, y = regime_difference), inherit.aes = FALSE, size = crossing_plot_point_size, colour = plot_palette$threshold) +
    annotate(
      "label",
      x = crossing_pct,
      y = 0,
      label = paste0("rank = ", crossing_rank),
      fill = "white",
      colour = plot_palette$threshold,
      size = 2.8,
      label.size = 0.15,
      vjust = -0.8
    ) +
    labs(
      title = paste0(comparison_name, " | regime-difference curve"),
      subtitle = compact_caption("Positive values indicate local IOD dominance and negative values indicate local CV² dominance. The first stable zero-crossing is the EVS cutoff.", width = 90),
      x = "Percentile center",
      y = "IOD − CV²"
    ) +
    manuscript_theme()
}

plot_crossing_summary_panel <- function(combined_cutoff_info, comparison_name) {
  p1 <- plot_combined_fourier_wave_map(combined_cutoff_info, comparison_name)
  p2 <- plot_regime_difference_curve(combined_cutoff_info, comparison_name)
  if (is.null(p1) || is.null(p2)) return(NULL)

  arrangeGrob(
    p1, p2,
    ncol = 1,
    top = textGrob(
      paste0(comparison_name, " | regime-shift crossing summary"),
      gp = gpar(fontface = "bold", cex = 1.04)
    )
  )
}

build_crossing_summary_table <- function(combined_cutoff_info, comparison_name) {
  sc <- combined_cutoff_info$selected_crossing
  if (is.null(sc) || !nrow(sc)) {
    return(data.frame(
      comparison_name = comparison_name,
      selected_reason = combined_cutoff_info$selected_reason,
      crossing_percentile = NA_real_,
      crossing_rank = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    comparison_name = comparison_name,
    selected_reason = combined_cutoff_info$selected_reason,
    crossing_id = sc$crossing_id[1],
    crossing_percentile = sc$crossing_percentile[1],
    crossing_rank = sc$crossing_rank[1],
    percentile_left = sc$percentile_left[1],
    percentile_right = sc$percentile_right[1],
    regime_difference_left = sc$regime_difference_left[1],
    regime_difference_right = sc$regime_difference_right[1],
    stringsAsFactors = FALSE
  )
}

resolve_evs_cutoff <- function(loading_tbl,
                               mean_col = "baseMean",
                               dispersion_col = "dispGeneEst",
                               fixed_top_n = evs_fixed_top_n) {
  if (isTRUE(evs_use_manual_rank_first) && is.finite(evs_fixed_rank_override) && !is.na(evs_fixed_rank_override)) {
    manual_rank <- max(1L, min(as.integer(evs_fixed_rank_override), nrow(loading_tbl)))
    return(list(
      cutoff_value = loading_tbl$pc1_loading_abs[manual_rank],
      top_n_actual = manual_rank,
      cutoff_quantile = 1 - (manual_rank / nrow(loading_tbl)),
      method = "manual_fixed_rank",
      curve_df = NULL,
      curvature_strength = NA_real_,
      candidate_table = data.frame(),
      quantile_screen_table = data.frame(),
      matched_interval_table = data.frame(),
      selected_reason = "manual_fixed_rank"
    ))
  }
  fallback <- resolve_top_n_cutoff(loading_tbl$pc1_loading_abs, top_n = fixed_top_n)
  list(
    cutoff_value = fallback$cutoff_value,
    top_n_actual = fallback$top_n_actual,
    cutoff_quantile = fallback$cutoff_quantile,
    method = "fixed_top_n_initial_placeholder",
    curve_df = NULL,
    curvature_strength = NA_real_,
    candidate_table = data.frame(),
    quantile_screen_table = data.frame(),
    matched_interval_table = data.frame(),
    selected_reason = "fixed_top_n_initial_placeholder"
  )
}

resolve_adaptive_evs_cutoff <- function(loading_tbl,
                                        mean_col = "baseMean",
                                        dispersion_col = "dispGeneEst",
                                        smoothing_spar = evs_curvature_spar,
                                        min_rank = evs_curvature_min_rank,
                                        max_rank_frac = evs_curvature_max_rank_frac,
                                        fixed_top_n = evs_fixed_top_n) {
  sorted_values_desc <- loading_tbl$pc1_loading_abs
  fallback <- resolve_top_n_cutoff(sorted_values_desc, top_n = fixed_top_n)

  curve_df <- smooth_ranked_nb_variance_curve(
    loading_tbl = loading_tbl,
    mean_col = mean_col,
    dispersion_col = dispersion_col,
    spar = smoothing_spar
  )

  if (is.null(curve_df) || !nrow(curve_df)) {
    fallback$curve_df <- curve_df
    fallback$curvature_strength <- NA_real_
    fallback$method <- "fixed_top_n_fallback"
    fallback$candidate_table <- data.frame(
      candidate_id = "fallback_top_n",
      rank_index = fallback$top_n_actual,
      curvature_strength = NA_real_,
      prominence = NA_real_,
      smoothed_log_nb_variance = NA_real_,
      first_derivative = NA_real_,
      second_derivative = NA_real_,
      max_rank_allowed = NA_integer_,
      cutoff_value = fallback$cutoff_value,
      cutoff_quantile = fallback$cutoff_quantile,
      selected = TRUE,
      selected_reason = "fallback_top_n",
      candidate_source = "fallback",
      stringsAsFactors = FALSE
    )
    fallback$selected_reason <- "fallback_top_n"
    return(fallback)
  }

  n_total <- nrow(loading_tbl)
  max_rank_allowed <- min(n_total, floor(n_total * max_rank_frac))

  quantile_screen_df <- summarize_curve_peak_grid(
    curve_df = curve_df,
    min_rank = min_rank,
    max_rank_frac = max_rank_frac,
    quantile_grid = evs_peak_quantile_grid,
    min_peak_distance = evs_candidate_min_peak_distance
  )

  quantile_screen_ranks <- select_quantile_screen_candidates(
    quantile_screen_df,
    q_min = evs_peak_quantile_screen_min,
    q_max = evs_peak_quantile_screen_max
  )

  quantile_screen_ranks <- sort(unique(as.integer(quantile_screen_ranks)))
  quantile_screen_ranks <- quantile_screen_ranks[
    is.finite(quantile_screen_ranks) &
      quantile_screen_ranks >= 1L &
      quantile_screen_ranks <= max_rank_allowed
  ]

  if (!length(quantile_screen_ranks)) {
    fallback$curve_df <- curve_df
    fallback$curvature_strength <- NA_real_
    fallback$method <- "fixed_top_n_fallback"
    fallback$candidate_table <- data.frame(
      candidate_id = "fallback_top_n",
      rank_index = fallback$top_n_actual,
      curvature_strength = NA_real_,
      prominence = NA_real_,
      smoothed_log_nb_variance = NA_real_,
      first_derivative = NA_real_,
      second_derivative = NA_real_,
      max_rank_allowed = max_rank_allowed,
      cutoff_value = fallback$cutoff_value,
      cutoff_quantile = fallback$cutoff_quantile,
      selected = TRUE,
      selected_reason = "fallback_top_n",
      candidate_source = "fallback",
      stringsAsFactors = FALSE
    )
    fallback$selected_reason <- "fallback_top_n"
    return(fallback)
  }

  peak_info <- find_local_curvature_peaks(
    curve_df = curve_df,
    min_rank = min_rank,
    max_rank_frac = max_rank_frac
  )
  peak_tbl <- if (!is.null(peak_info) && !is.null(peak_info$peak_table)) peak_info$peak_table else NULL

  candidate_tbl <- data.frame(
    rank_index = quantile_screen_ranks,
    stringsAsFactors = FALSE
  )

  if (!is.null(peak_tbl) && nrow(peak_tbl)) {
    peak_map <- peak_tbl[, c("rank_index", "curvature_strength", "prominence",
                             "smoothed_log_nb_variance", "first_derivative", "second_derivative"),
                         drop = FALSE]
    candidate_tbl <- dplyr::left_join(candidate_tbl, peak_map, by = "rank_index")
  } else {
    candidate_tbl$curvature_strength <- NA_real_
    candidate_tbl$prominence <- NA_real_
    candidate_tbl$smoothed_log_nb_variance <- NA_real_
    candidate_tbl$first_derivative <- NA_real_
    candidate_tbl$second_derivative <- NA_real_
  }

  candidate_tbl$candidate_source <- "quantile_screen_peak"
  candidate_tbl$candidate_id <- paste0("candidate_", seq_len(nrow(candidate_tbl)))
  candidate_tbl$max_rank_allowed <- as.integer(max_rank_allowed)
  candidate_tbl$cutoff_value <- loading_tbl$pc1_loading_abs[candidate_tbl$rank_index]
  candidate_tbl$cutoff_quantile <- 1 - (candidate_tbl$rank_index / n_total)
  candidate_tbl$selected <- FALSE
  candidate_tbl$selected_reason <- "candidate_only"

  candidate_tbl <- candidate_tbl[order(candidate_tbl$rank_index), , drop = FALSE]
  rownames(candidate_tbl) <- NULL

  ranking_tbl <- candidate_tbl
  ranking_tbl$curvature_rank <- rank(-ranking_tbl$curvature_strength, ties.method = "min", na.last = "keep")
  ranking_tbl$prominence_rank <- rank(-ranking_tbl$prominence, ties.method = "min", na.last = "keep")

  if (all(is.na(ranking_tbl$curvature_strength))) {
    best_idx <- 1L
    selection_reason <- "initial_quantile_screen_peak_earliest_rank_fallback"
  } else {
    best_idx <- with(
      ranking_tbl,
      order(
        dplyr::coalesce(curvature_rank, Inf),
        dplyr::coalesce(prominence_rank, Inf),
        rank_index
      )[1]
    )
    selection_reason <- "initial_quantile_screen_peak_max_curvature"
  }

  candidate_tbl$selected[best_idx] <- TRUE
  candidate_tbl$selected_reason[best_idx] <- selection_reason

  list(
    top_n_actual       = candidate_tbl$rank_index[best_idx],
    cutoff_value       = candidate_tbl$cutoff_value[best_idx],
    cutoff_quantile    = candidate_tbl$cutoff_quantile[best_idx],
    n_total            = n_total,
    method             = "variance_second_derivative_quantile_screen",
    curvature_strength = candidate_tbl$curvature_strength[best_idx],
    curve_df           = curve_df,
    max_rank_allowed   = as.integer(max_rank_allowed),
    variance_measure   = "NB2_variance",
    candidate_table    = candidate_tbl,
    quantile_screen_table = quantile_screen_df,
    selected_reason    = selection_reason
  )
}

cap_evs_candidate_table <- function(candidate_tbl,
                                    max_candidates = Inf,
                                    n_total = NA_integer_) {
  if (is.null(candidate_tbl) || !nrow(candidate_tbl)) return(candidate_tbl)

  candidate_tbl <- as.data.frame(candidate_tbl, stringsAsFactors = FALSE)

  if (!"rank_index" %in% colnames(candidate_tbl)) {
    stop("candidate_tbl must contain rank_index.", call. = FALSE)
  }

  candidate_tbl$rank_index <- suppressWarnings(as.integer(candidate_tbl$rank_index))
  candidate_tbl <- candidate_tbl[!is.na(candidate_tbl$rank_index) &
                                   is.finite(candidate_tbl$rank_index) &
                                   candidate_tbl$rank_index >= 1L, , drop = FALSE]

  if (is.finite(n_total) && !is.na(n_total)) {
    candidate_tbl <- candidate_tbl[candidate_tbl$rank_index <= as.integer(n_total), , drop = FALSE]
  }

  candidate_tbl <- candidate_tbl[order(candidate_tbl$rank_index), , drop = FALSE]
  candidate_tbl <- candidate_tbl[!duplicated(candidate_tbl$rank_index), , drop = FALSE]

  if (!"candidate_id" %in% colnames(candidate_tbl)) {
    candidate_tbl$candidate_id <- paste0("candidate_", seq_len(nrow(candidate_tbl)))
  } else {
    missing_ids <- is.na(candidate_tbl$candidate_id) | candidate_tbl$candidate_id == ""
    candidate_tbl$candidate_id <- as.character(candidate_tbl$candidate_id)
    candidate_tbl$candidate_id[missing_ids] <- paste0("candidate_", which(missing_ids))
  }

  if (!"selected" %in% colnames(candidate_tbl)) {
    candidate_tbl$selected <- FALSE
  }

  if (!"selected_reason" %in% colnames(candidate_tbl)) {
    candidate_tbl$selected_reason <- "candidate_only"
  }

  if (is.finite(max_candidates) && nrow(candidate_tbl) > max_candidates) {
    candidate_tbl <- candidate_tbl[seq_len(max_candidates), , drop = FALSE]
  }

  rownames(candidate_tbl) <- NULL
  candidate_tbl
}

run_core_analysis_candidate_score <- function(count_mat,
                                              coldata,
                                              dataset_name,
                                              annot_df = NULL) {
  if (is.null(count_mat) || !nrow(count_mat) || !ncol(count_mat)) {
    stop(sprintf("[%s] count_mat must contain at least one feature and one sample.", dataset_name), call. = FALSE)
  }

  if (is.null(rownames(count_mat)) || anyNA(rownames(count_mat)) || any(rownames(count_mat) == "")) {
    stop(sprintf("[%s] count_mat must have non-empty rownames for all features.", dataset_name), call. = FALSE)
  }

  if (is.null(annot_df)) {
    annot_df <- data.frame(
      feature_id = as.character(rownames(count_mat)),
      gene_symbol = as.character(rownames(count_mat)),
      stringsAsFactors = FALSE
    )
  }

  fit <- run_core_analysis(
    count_mat = count_mat,
    coldata = coldata,
    dataset_name = dataset_name,
    annot_df = annot_df
  )

  if (is.null(fit$results) || is.null(fit$hc_p_threshold)) {
    stop(sprintf("[%s] Candidate-scoring analysis did not return the expected results/hc_p_threshold fields.", dataset_name), call. = FALSE)
  }

  fit
}

# -----------------------------------------------------------------------------
# PATCHED APPLY LOADING CUTOFF
# -----------------------------------------------------------------------------

apply_loading_cutoff <- function(fit_obj,
                                 rank_index,
                                 selected_reason = "manual_override") {
  fit_obj <- as.list(fit_obj)
  loading_tbl <- as.data.frame(fit_obj$loading_table, stringsAsFactors = FALSE)

  if (is.null(loading_tbl) || !nrow(loading_tbl)) {
    stop("apply_loading_cutoff() requires a non-empty loading_table.", call. = FALSE)
  }

  assert_required_columns(
    loading_tbl,
    c("feature_id", "pc1_loading_abs", "rank"),
    object_name = "fit_obj$loading_table"
  )

  loading_tbl <- loading_tbl[order(loading_tbl$rank), , drop = FALSE]
  rownames(loading_tbl) <- NULL

  n_total <- nrow(loading_tbl)
  rank_index <- as.integer(rank_index)[1]

  if (!is.finite(rank_index) || is.na(rank_index)) {
    stop("apply_loading_cutoff() requires rank_index to be a finite integer.", call. = FALSE)
  }

  rank_index <- min(max(1L, rank_index), n_total)
  row_idx <- match(rank_index, loading_tbl$rank)

  if (is.na(row_idx) || length(row_idx) != 1L) {
    stop(sprintf("apply_loading_cutoff() could not locate rank %d.", rank_index), call. = FALSE)
  }

  cutoff_value <- as.numeric(loading_tbl$pc1_loading_abs[row_idx])
  cutoff_quantile <- 1 - (rank_index / n_total)

  fit_obj$loading_table <- loading_tbl
  fit_obj$cutoff <- cutoff_value
  fit_obj$top_n_used <- rank_index
  fit_obj$cutoff_quantile <- cutoff_quantile
  fit_obj$curvature_strength <- suppressWarnings({
    tbl <- fit_obj$candidate_table
    if (!is.null(tbl) && nrow(tbl) && any(tbl$rank_index == rank_index)) {
      tbl$curvature_strength[match(rank_index, tbl$rank_index)]
    } else {
      NA_real_
    }
  })
  fit_obj$cutoff_method <- paste0(fit_obj$cutoff_method, "_selected")
  fit_obj$selected_reason <- selected_reason

  fit_obj$loading_table$split_class <- ifelse(
    fit_obj$loading_table$rank <= rank_index,
    "high_loading",
    "background_loading"
  )

  if (!is.null(fit_obj$candidate_table) && nrow(fit_obj$candidate_table)) {
    fit_obj$candidate_table$selected <- fit_obj$candidate_table$rank_index == rank_index
    fit_obj$candidate_table$selected_reason <- ifelse(
      fit_obj$candidate_table$selected,
      selected_reason,
      "candidate_only"
    )
  }

  fit_obj
}

score_evs_candidate_result <- function(df) {
  if (is.null(df) || !nrow(df)) return(list(score = NA_real_, metrics = NULL))

  n_standard <- sum(df$standard_significant, na.rm = TRUE)

  evidence <- log1p(n_standard)

  list(
    score = evidence,
    metrics = c(
      n_standard = n_standard
    )
  )
}

compute_nb1_subset_moment_summary <- function(count_submatrix) {
  mat <- as.matrix(count_submatrix)
  if (is.null(mat) || !length(mat)) {
    return(c(mean_value = NA_real_, variance_value = NA_real_, iod = NA_real_, cv2 = NA_real_))
  }

  mat <- suppressWarnings(matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)))
  if (!nrow(mat) || !ncol(mat)) {
    return(c(mean_value = NA_real_, variance_value = NA_real_, iod = NA_real_, cv2 = NA_real_))
  }

  feature_mean <- rowMeans(mat, na.rm = TRUE)
  feature_var <- apply(mat, 1L, stats::var, na.rm = TRUE)
  keep <- is.finite(feature_mean) & !is.na(feature_mean) & feature_mean > 0 &
    is.finite(feature_var) & !is.na(feature_var) & feature_var >= 0

  if (!any(keep)) {
    return(c(mean_value = NA_real_, variance_value = NA_real_, iod = NA_real_, cv2 = NA_real_))
  }

  mean_value <- mean(feature_mean[keep], na.rm = TRUE)
  variance_value <- mean(feature_var[keep], na.rm = TRUE)
  iod <- if (is.finite(mean_value) && !is.na(mean_value) && mean_value > 0) variance_value / mean_value else NA_real_
  cv2 <- if (is.finite(mean_value) && !is.na(mean_value) && mean_value > 0) variance_value / (mean_value ^ 2) else NA_real_

  c(mean_value = mean_value, variance_value = variance_value, iod = iod, cv2 = cv2)
}

compute_nb2_subset_moment_summary <- function(results_df) {
  if (is.null(results_df) || !nrow(results_df)) {
    return(c(mean_value = NA_real_, dispersion_value = NA_real_, iod = NA_real_, cv2 = NA_real_))
  }

  df <- as.data.frame(results_df, stringsAsFactors = FALSE)
  if (!all(c("baseMean", "dispersion") %in% colnames(df))) {
    return(c(mean_value = NA_real_, dispersion_value = NA_real_, iod = NA_real_, cv2 = NA_real_))
  }

  df$baseMean <- suppressWarnings(as.numeric(df$baseMean))
  df$dispersion <- suppressWarnings(as.numeric(df$dispersion))
  keep <- is.finite(df$baseMean) & !is.na(df$baseMean) & df$baseMean > 0 &
    is.finite(df$dispersion) & !is.na(df$dispersion) & df$dispersion >= 0

  if (!any(keep)) {
    return(c(mean_value = NA_real_, dispersion_value = NA_real_, iod = NA_real_, cv2 = NA_real_))
  }

  mean_value <- mean(df$baseMean[keep], na.rm = TRUE)
  dispersion_value <- mean(df$dispersion[keep], na.rm = TRUE)
  iod <- mean(1 + (df$dispersion[keep] * df$baseMean[keep]), na.rm = TRUE)
  cv2 <- mean((1 / df$baseMean[keep]) + df$dispersion[keep], na.rm = TRUE)

  c(mean_value = mean_value, dispersion_value = dispersion_value, iod = iod, cv2 = cv2)
}

format_runtime_minutes <- function(seconds_value) {
  if (!is.finite(seconds_value) || is.na(seconds_value) || seconds_value < 0) return("unknown")
  if (seconds_value < 60) return(sprintf("~%d sec", as.integer(round(seconds_value))))
  sprintf("~%.1f min", seconds_value / 60)
}

subset_normalized_counts_by_ids <- function(normalized_counts,
                                            feature_ids,
                                            dataset_label,
                                            comparison_name,
                                            strict = TRUE) {
  feature_ids <- as.character(feature_ids)
  idx <- match(feature_ids, rownames(normalized_counts))
  missing_ids <- unique(feature_ids[is.na(idx)])

  if (length(missing_ids) > 0L) {
    msg <- sprintf(
      "[%s][%s] %d feature_id(s) were not found in normalized_counts after all-zero filtering. Example IDs: %s",
      comparison_name,
      dataset_label,
      length(missing_ids),
      paste(utils::head(missing_ids, 10L), collapse = ", ")
    )
    if (isTRUE(strict)) stop(msg, call. = FALSE) else warning(msg)
  }

  if (!isTRUE(strict)) idx <- idx[!is.na(idx)]

  as.matrix(normalized_counts[idx, , drop = FALSE])
}

print_evs_single_dataset_runtime_preview <- function(candidate_tbl,
                                                     comparison_name,
                                                     dataset_label,
                                                     seconds_per_run_low = evs_candidate_runtime_seconds_low,
                                                     seconds_per_run_high = evs_candidate_runtime_seconds_high) {
  candidate_tbl <- cap_evs_candidate_table(candidate_tbl)
  n_candidates <- if (is.null(candidate_tbl) || !nrow(candidate_tbl)) 0L else nrow(candidate_tbl)

  runtime_low_seconds  <- n_candidates * 2L * seconds_per_run_low
  runtime_high_seconds <- n_candidates * 2L * seconds_per_run_high

  message("\n======================================================")
  message("EVS SINGLE-DATASET CUTOFF PREVIEW: ", comparison_name, " | ", dataset_label)
  message("------------------------------------------------------")
  message("Candidate cutoffs to evaluate: ", n_candidates)
  message("Total candidate-scoring DESeq2 runs: ", n_candidates * 2L)
  message("Approximate runtime: ",
          format_runtime_minutes(runtime_low_seconds), " to ",
          format_runtime_minutes(runtime_high_seconds))
  if (n_candidates > 0L) {
    message("Candidate ranks: ", paste(candidate_tbl$rank_index, collapse = ", "))
  }
  message("======================================================\n")
}

strip_dataset_key <- function(dataset_name) {
  sub("^.*_(raw_dataset|leading_edge_dataset|remainder_dataset)$", "\\1", dataset_name)
}

comparison_label_from_dataset_name <- function(dataset_name) {
  sub("_(raw_dataset|leading_edge_dataset|remainder_dataset)$", "", dataset_name)
}

clean_gene_set <- function(x) {
  unique(tolower(trimws(x[!is.na(x) & x != ""])))
}

save_csv <- function(df, path) {
  write.csv(df, file = path, row.names = FALSE)
}

save_grob <- function(g, path, width = 16.4, height = 9.9, dpi = figure_dpi, bg = "white") {
  tryCatch(
    {
      ggsave(
        filename  = path,
        plot      = g,
        width     = width,
        height    = height,
        dpi       = dpi,
        units     = "in",
        bg        = bg,
        limitsize = FALSE
      )
      invisible(TRUE)
    },
    error = function(e) {
      warning(paste0("Panel export failed for ", basename(path), ": ", conditionMessage(e)))
      invisible(FALSE)
    }
  )
}

safe_plot_build <- function(expr, label = "plot") {
  tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) {
      warning(paste0(label, " failed: ", conditionMessage(e)))
      NULL
    }
  )
}

pretty_dataset_type <- function(dataset_key) {
  switch(
    dataset_key,
    raw_dataset          = "Original dataset",
    leading_edge_dataset = "Leading-edge dataset",
    remainder_dataset    = "Remainder dataset",
    dataset_key
  )
}

pretty_dataset_label <- function(dataset_name) {
  dataset_key     <- strip_dataset_key(dataset_name)
  comparison_name <- comparison_label_from_dataset_name(dataset_name)
  paste(comparison_name, pretty_dataset_type(dataset_key), sep = " | ")
}

pretty_group_label <- function(group_label) {
  switch(
    group_label,
    trt       = "Treatment",
    untrt     = "Control",
    treatment = "Treatment",
    control   = "Control",
    group_label
  )
}

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title        = element_text(
        face       = "bold",
        size       = base_theme_size + 1,
        hjust      = 0.5,
        lineheight = 1.00,
        margin     = margin(b = 5)
      ),
      plot.subtitle     = element_text(
        size       = base_theme_size - 1,
        hjust      = 0.5,
        lineheight = 1.00,
        margin     = margin(b = 7)
      ),
      plot.caption      = element_text(
        size       = base_theme_size - 3,
        hjust      = 0.5,
        colour     = "grey30",
        lineheight = 0.98,
        margin     = margin(t = 8)
      ),
      axis.title        = element_text(face = "bold"),
      axis.text         = element_text(colour = "black"),
      legend.title      = element_text(face = "bold"),
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.box.margin = margin(t = 3, r = 3, b = 3, l = 3),
      legend.margin     = margin(t = 2, r = 2, b = 2, l = 2),
      legend.spacing.x  = unit(5, "pt"),
      legend.spacing.y  = unit(2, "pt"),
      legend.text       = element_text(size = base_theme_size - 1),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(linewidth = 0.25, colour = "grey88"),
      plot.margin       = margin(t = 18, r = 30, b = 22, l = 24)
    )
}

make_design_formula <- function() {
  ~ condition
}

get_condition_coef <- function(dds) {
  rn  <- resultsNames(dds)
  idx <- grep("^condition_", rn)
  if (length(idx) == 0) stop("Could not identify condition coefficient in resultsNames(dds).", call. = FALSE)
  if (length(idx) > 1) message("Multiple condition coefficients found; using: ", rn[idx[1]])
  rn[idx[1]]
}

classify_effect_strength <- function(res_strong_padj, res_weak_padj, alpha = alpha_level) {
  out <- rep("intermediate", length(res_strong_padj))
  out[!is.na(res_weak_padj)   & res_weak_padj   < alpha] <- "weak_effect"
  out[!is.na(res_strong_padj) & res_strong_padj < alpha] <- "strong_effect"
  out
}

# -----------------------------------------------------------------------------
# IMPORT WTTS MASTER COUNT FILE
# -----------------------------------------------------------------------------

if (!file.exists(count_file)) stop(sprintf("Count file not found: %s", count_file), call. = FALSE)

WTTS_Seq <- read.csv(
  count_file,
  header        = TRUE,
  stringsAsFactors = FALSE,
  check.names   = FALSE
)

WTTS_Seq <- as.data.frame(WTTS_Seq, stringsAsFactors = FALSE)
WTTS_Seq$OrigID <- as.character(WTTS_Seq$OrigID)
WTTS_Seq$Symbol <- as.character(WTTS_Seq$Symbol)

assert_required_columns(WTTS_Seq, c("OrigID", "Symbol"), object_name = "WTTS count file")
assert_required_columns(WTTS_Seq, meta_all$id, object_name = "WTTS count file sample columns")

WTTS_Seq <- WTTS_Seq[!is.na(WTTS_Seq$OrigID) & !is.na(WTTS_Seq$Symbol), , drop = FALSE]

sample_na <- rowSums(is.na(WTTS_Seq[, meta_all$id, drop = FALSE])) > 0
WTTS_Seq  <- WTTS_Seq[!sample_na, , drop = FALSE]

rownames(WTTS_Seq) <- WTTS_Seq$OrigID

OrigID_Symbol <- unique(WTTS_Seq[, c("OrigID", "Symbol"), drop = FALSE])
colnames(OrigID_Symbol) <- c("feature_id", "gene_symbol")
OrigID_Symbol$feature_id  <- as.character(OrigID_Symbol$feature_id)
OrigID_Symbol$gene_symbol <- as.character(OrigID_Symbol$gene_symbol)

OrigID_Symbol <- OrigID_Symbol %>%
  dplyr::mutate(gene_symbol = dplyr::if_else(is.na(gene_symbol), "", trimws(gene_symbol))) %>%
  dplyr::arrange(feature_id, dplyr::desc(gene_symbol != ""), gene_symbol) %>%
  dplyr::distinct(feature_id, .keep_all = TRUE) %>%
  dplyr::mutate(gene_symbol = dplyr::na_if(gene_symbol, ""))

if (run_twas_overlap) {
  if (!file.exists(twas_file)) stop(sprintf("TWAS file not found: %s", twas_file), call. = FALSE)
  TWAS_Seq <- read.csv(twas_file, header = TRUE, stringsAsFactors = FALSE)
  TWAS_Seq <- as.data.frame(TWAS_Seq)
  if (ncol(TWAS_Seq) < 4) stop("TWAS file must contain at least 4 columns.", call. = FALSE)
  TWAS_data <- TWAS_Seq[, c(1, 4), drop = FALSE]
  colnames(TWAS_data) <- c("source_id", "gene_symbol")
}

# =============================================================================
# SECTION 2 OF 4
# COMPARISON PREPARATION, EVS, EVS FIGURES, PCA SUMMARIES
# =============================================================================

prepare_comparison_data <- function(comparison_name, group1_prefix, group2_prefix,
                                    WTTS_Seq, meta_all) {
  keep_ids <- grepl(paste0("^", group1_prefix, "_"), meta_all$id) |
    grepl(paste0("^", group2_prefix, "_"), meta_all$id)

  meta_sub   <- meta_all[keep_ids, , drop = FALSE]
  coldata    <- meta_sub[, c("condition"), drop = FALSE]
  sample_ids <- rownames(meta_sub)

  missing_samples <- setdiff(sample_ids, colnames(WTTS_Seq))
  if (length(missing_samples) > 0) {
    stop(
      paste(
        "Missing samples in WTTS file for", comparison_name, ":",
        paste(missing_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  count_sub <- WTTS_Seq[, sample_ids, drop = FALSE]
  stopifnot(all(colnames(count_sub) == rownames(coldata)))

  list(
    comparison_name = comparison_name,
    count_matrix    = as.matrix(count_sub),
    coldata         = coldata
  )
}

compute_condition_feature_metrics <- function(count_submatrix) {
  cd <- S4Vectors::DataFrame(row.names = colnames(count_submatrix))
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_submatrix)),
    colData   = cd,
    design    = ~ 1
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds, quiet = TRUE)
  dds <- estimateDispersionsFit(dds, quiet = TRUE)

  md <- as.data.frame(SummarizedExperiment::mcols(dds), stringsAsFactors = FALSE)
  md$feature_id <- rownames(md)
  keep_cols <- intersect(c("feature_id", "baseMean", "dispGeneEst", "dispFit", "dispersion"), colnames(md))
  md <- md[, keep_cols, drop = FALSE]
  md
}

compute_pc1_loading_table <- function(value_df, sample_names, top_n = evs_fixed_top_n,
                                      preprocessing_label = "Normalized prior to eigenvector splitting",
                                      feature_metrics = NULL,
                                      determine_cutoff = TRUE,
                                      cutoff_method_label_if_skipped = "inherited_from_normalized") {
  x <- as.matrix(value_df[, sample_names, drop = FALSE])

  max_pcs <- min(5L, ncol(x))
  pca_fit <- prcomp(t(x), scale. = FALSE, rank. = max_pcs)

  loading_abs <- abs(pca_fit$rotation[, 1])
  loading_tbl <- data.frame(
    feature_id       = names(loading_abs),
    pc1_loading      = unname(pca_fit$rotation[, 1]),
    pc1_loading_abs  = unname(loading_abs),
    stringsAsFactors = FALSE
  )

  if (!is.null(feature_metrics) && nrow(feature_metrics) > 0) {
    fm <- as.data.frame(feature_metrics, stringsAsFactors = FALSE)
    if (!"feature_id" %in% colnames(fm)) stop("feature_metrics must contain feature_id.", call. = FALSE)
    keep_cols <- intersect(c("feature_id", "baseMean", "dispGeneEst", "dispFit", "dispersion"), colnames(fm))
    fm <- fm[, keep_cols, drop = FALSE]
    fm <- fm[!duplicated(fm$feature_id), , drop = FALSE]
    loading_tbl <- dplyr::left_join(loading_tbl, fm, by = "feature_id")
  }

  loading_tbl <- loading_tbl[order(loading_tbl$pc1_loading_abs, decreasing = TRUE), ]
  loading_tbl$rank <- seq_len(nrow(loading_tbl))

  cutoff_info <- if (!isTRUE(determine_cutoff)) {
    list(
      cutoff_value          = NA_real_,
      top_n_actual          = NA_integer_,
      cutoff_quantile       = NA_real_,
      method                = cutoff_method_label_if_skipped,
      curve_df              = NULL,
      curvature_strength    = NA_real_,
      candidate_table       = data.frame(),
      quantile_screen_table = data.frame(),
      matched_interval_table = data.frame(),
      selected_reason       = cutoff_method_label_if_skipped
    )
  } else {
    resolve_evs_cutoff(
      loading_tbl,
      fixed_top_n = top_n,
      mean_col = "baseMean",
      dispersion_col = "dispGeneEst"
    )
  }
  cutoff <- cutoff_info$cutoff_value

  loading_tbl$split_class <- if (is.finite(cutoff)) {
    ifelse(
      loading_tbl$rank <= cutoff_info$top_n_actual,
      "high_loading",
      "background_loading"
    )
  } else {
    rep("background_loading", nrow(loading_tbl))
  }

  list(
    pca_fit               = pca_fit,
    loading_table         = loading_tbl,
    cutoff                = cutoff,
    top_n_used            = cutoff_info$top_n_actual,
    cutoff_quantile       = cutoff_info$cutoff_quantile,
    cutoff_method         = cutoff_info$method,
    cutoff_curve_df       = cutoff_info$curve_df,
    curvature_strength    = cutoff_info$curvature_strength,
    candidate_table       = cutoff_info$candidate_table,
    quantile_screen_table = cutoff_info$quantile_screen_table,
    matched_interval_table = cutoff_info$matched_interval_table,
    selected_reason       = cutoff_info$selected_reason,
    preprocessing_label   = preprocessing_label
  )
}

# -----------------------------------------------------------------------------
# PATCHED BUILD EIGENVECTOR SPLIT
# -----------------------------------------------------------------------------

build_eigenvector_split <- function(count_matrix, coldata, comparison_name) {
  design_formula <- make_design_formula()

  dds_init <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData   = coldata,
    design    = design_formula
  )

  dds_init <- dds_init[rowSums(counts(dds_init)) > 0, ]
  dds_init <- estimateSizeFactors(dds_init)

  retained_feature_ids <- rownames(dds_init)
  norm_counts_init <- as.data.frame(counts(dds_init, normalized = TRUE))
  raw_counts_init  <- as.data.frame(count_matrix[retained_feature_ids, , drop = FALSE])

  sample_ids <- colnames(count_matrix)
  trt_ids    <- sample_ids[coldata$condition == "trt"]
  untrt_ids  <- sample_ids[coldata$condition == "untrt"]

  feature_metrics_trt <- compute_condition_feature_metrics(count_matrix[, trt_ids, drop = FALSE])
  feature_metrics_untrt <- compute_condition_feature_metrics(count_matrix[, untrt_ids, drop = FALSE])

  fit_trt <- compute_pc1_loading_table(
    norm_counts_init,
    trt_ids,
    top_n = evs_fixed_top_n,
    preprocessing_label = "Normalized before EVS",
    feature_metrics = feature_metrics_trt,
    determine_cutoff = FALSE,
    cutoff_method_label_if_skipped = "normalized_shared_cutoff_pending"
  )

  fit_untrt <- compute_pc1_loading_table(
    norm_counts_init,
    untrt_ids,
    top_n = evs_fixed_top_n,
    preprocessing_label = "Normalized before EVS",
    feature_metrics = feature_metrics_untrt,
    determine_cutoff = FALSE,
    cutoff_method_label_if_skipped = "normalized_shared_cutoff_pending"
  )

  fit_trt_raw <- compute_pc1_loading_table(
    raw_counts_init,
    trt_ids,
    top_n = evs_fixed_top_n,
    preprocessing_label = "Raw counts before EVS",
    feature_metrics = feature_metrics_trt,
    determine_cutoff = FALSE,
    cutoff_method_label_if_skipped = "raw_cutoff_inherited_from_normalized"
  )

  fit_untrt_raw <- compute_pc1_loading_table(
    raw_counts_init,
    untrt_ids,
    top_n = evs_fixed_top_n,
    preprocessing_label = "Raw counts before EVS",
    feature_metrics = feature_metrics_untrt,
    determine_cutoff = FALSE,
    cutoff_method_label_if_skipped = "raw_cutoff_inherited_from_normalized"
  )

  combined_cutoff_info <- resolve_combined_fourier_cutoff(
    fit_trt_loading_tbl  = fit_trt$loading_table,
    fit_ctrl_loading_tbl = fit_untrt$loading_table,
    fixed_top_n          = evs_fixed_top_n
  )
  combined_cutoff_info <- build_group_backup_crossings(combined_cutoff_info)

  final_shared_rank <- as.integer(combined_cutoff_info$top_n_actual)
  final_shared_reason <- combined_cutoff_info$selected_reason
  raw_wald_search_info <- NULL

  if (isTRUE(raw_wald_cutoff_mode)) {
    full_raw_fit <- run_core_analysis(
      count_mat = count_matrix,
      coldata = coldata,
      dataset_name = paste0(comparison_name, "_raw_wald_search_reference"),
      annot_df = OrigID_Symbol
    )

    trt_backup_rank <- if (!is.null(combined_cutoff_info$trt_backup_crossing) && nrow(combined_cutoff_info$trt_backup_crossing)) {
      as.integer(combined_cutoff_info$trt_backup_crossing$crossing_rank[1])
    } else {
      as.integer(combined_cutoff_info$top_n_actual)
    }

    ctrl_backup_rank <- if (!is.null(combined_cutoff_info$ctrl_backup_crossing) && nrow(combined_cutoff_info$ctrl_backup_crossing)) {
      as.integer(combined_cutoff_info$ctrl_backup_crossing$crossing_rank[1])
    } else {
      as.integer(combined_cutoff_info$top_n_actual)
    }

    bracket_lo <- min(trt_backup_rank, ctrl_backup_rank, na.rm = TRUE)
    bracket_hi <- max(trt_backup_rank, ctrl_backup_rank, na.rm = TRUE)

    bracket_pad <- max(10L, round(0.05 * max(1L, abs(bracket_hi - bracket_lo))))
    bracket_lo <- max(1L, bracket_lo - bracket_pad)
    bracket_hi <- min(nrow(combined_cutoff_info$shared_combined_tbl), bracket_hi + bracket_pad)

    raw_wald_search_info <- run_raw_wald_bracket_search(
      full_results_df = full_raw_fit$results,
      shared_combined_tbl = combined_cutoff_info$shared_combined_tbl,
      lower_rank = bracket_lo,
      upper_rank = bracket_hi
    )

    final_shared_rank <- as.integer(raw_wald_search_info$selected_rank)
    final_shared_reason <- "raw_wald_bracket_balanced_normality_selected"
  }

  fit_trt <- apply_loading_cutoff(
    fit_trt,
    rank_index = final_shared_rank,
    selected_reason = paste0(final_shared_reason, "_applied_to_treatment")
  )
  fit_untrt <- apply_loading_cutoff(
    fit_untrt,
    rank_index = final_shared_rank,
    selected_reason = paste0(final_shared_reason, "_applied_to_control")
  )
  fit_trt_raw <- apply_loading_cutoff(
    fit_trt_raw,
    rank_index = final_shared_rank,
    selected_reason = paste0(final_shared_reason, "_projected_to_raw_treatment")
  )
  fit_trt_raw$cutoff_method <- "crossing_rank_projected_to_raw_selected"
  fit_untrt_raw <- apply_loading_cutoff(
    fit_untrt_raw,
    rank_index = final_shared_rank,
    selected_reason = paste0(final_shared_reason, "_projected_to_raw_control")
  )
  fit_untrt_raw$cutoff_method <- "crossing_rank_projected_to_raw_selected"

  primary_trt_fit <- fit_trt
  primary_untrt_fit <- fit_untrt
  independent_candidate_evals <- list()

  shared_combined_tbl <- combined_cutoff_info$shared_combined_tbl
  if (is.null(shared_combined_tbl) || !nrow(shared_combined_tbl)) {
    stop("Shared combined loading table is missing after shared cutoff resolution.", call. = FALSE)
  }

  leading_edge_ids <- as.character(
    shared_combined_tbl$feature_id[shared_combined_tbl$combined_rank <= final_shared_rank]
  )
  analyzed_feature_ids <- as.character(shared_combined_tbl$feature_id)
  remainder_ids <- setdiff(analyzed_feature_ids, leading_edge_ids)

  if (length(leading_edge_ids) == 0) {
    stop("Leading-edge dataset is empty. Check sample mapping or EVS cutoff settings.", call. = FALSE)
  }
  if (length(remainder_ids) == 0) {
    stop("Remainder dataset is empty. Check shared EVS cutoff settings.", call. = FALSE)
  }

  evs_cutoff_summary <- dplyr::bind_rows(
    data.frame(
      preprocessing = "normalized",
      group = "regime_shift_crossing",
      cutoff_mode = combined_cutoff_info$method,
      fixed_top_n_requested = evs_fixed_top_n,
      empiric_rank_selected = final_shared_rank,
      cutoff_quantile = fit_trt$cutoff_quantile,
      curvature_strength = NA_real_,
      selected_reason = combined_cutoff_info$selected_reason,
      stringsAsFactors = FALSE
    ),
    data.frame(
      preprocessing = "raw",
      group = "regime_shift_crossing",
      cutoff_mode = "crossing_rank_projected_to_raw_selected",
      fixed_top_n_requested = evs_fixed_top_n,
      empiric_rank_selected = final_shared_rank,
      cutoff_quantile = fit_trt_raw$cutoff_quantile,
      curvature_strength = NA_real_,
      selected_reason = paste0(combined_cutoff_info$selected_reason, "_projected_to_raw"),
      stringsAsFactors = FALSE
    )
  )

  if (!is.null(raw_wald_search_info) && !is.null(raw_wald_search_info$selected_row)) {
    evs_cutoff_summary <- dplyr::bind_rows(
      evs_cutoff_summary,
      data.frame(
        preprocessing = "raw_wald_rank_search",
        group = "balanced_normality_search",
        cutoff_mode = "raw_wald_bracket_balanced_normality",
        fixed_top_n_requested = evs_fixed_top_n,
        empiric_rank_selected = final_shared_rank,
        cutoff_quantile = 1 - (final_shared_rank / nrow(shared_combined_tbl)),
        curvature_strength = NA_real_,
        selected_reason = final_shared_reason,
        stringsAsFactors = FALSE
      )
    )
  }

  list(
    fit_trt               = fit_trt,
    fit_untrt             = fit_untrt,
    fit_trt_raw           = fit_trt_raw,
    fit_untrt_raw         = fit_untrt_raw,
    primary_trt_fit       = primary_trt_fit,
    primary_untrt_fit     = primary_untrt_fit,
    independent_candidate_evals = independent_candidate_evals,
    combined_cutoff_info  = combined_cutoff_info,
    shared_combined_tbl   = shared_combined_tbl,
    evs_cutoff_summary    = evs_cutoff_summary,
    raw_wald_search_info  = raw_wald_search_info,
    normalized_counts     = norm_counts_init,
    raw_dataset           = count_matrix,
    leading_edge_dataset  = count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset     = count_matrix[remainder_ids, , drop = FALSE]
  )
}

condition_shapes <- c("untrt" = 21, "trt" = 25)
condition_fills  <- c("untrt" = plot_palette$control, "trt" = plot_palette$treatment)
condition_labels <- c("untrt" = "Control", "trt" = "Treatment")

make_pca_summary_table <- function(pca_fit, dataset_name, preprocessing_label, group_label = NA_character_) {
  sdev     <- pca_fit$sdev
  variance <- sdev^2
  prop_var <- variance / sum(variance)
  cum_var  <- cumsum(prop_var)

  data.frame(
    dataset_name           = dataset_name,
    preprocessing_label    = preprocessing_label,
    group_label            = group_label,
    principal_component    = paste0("PC", seq_along(sdev)),
    standard_deviation     = as.numeric(sdev),
    variance               = as.numeric(variance),
    proportion_variance    = as.numeric(prop_var),
    cumulative_proportion  = as.numeric(cum_var),
    stringsAsFactors       = FALSE
  )
}

plot_pca_variance_profile <- function(pca_fit, dataset_name, preprocessing_label, group_label = NULL) {
  pca_tbl <- make_pca_summary_table(
    pca_fit = pca_fit,
    dataset_name = dataset_name,
    preprocessing_label = preprocessing_label,
    group_label = ifelse(is.null(group_label), NA_character_, group_label)
  )

  pca_tbl$principal_component <- factor(
    pca_tbl$principal_component,
    levels = pca_tbl$principal_component
  )

  group_label_chr <- if (length(group_label) && !is.null(group_label[1])) tolower(as.character(group_label[1])) else NA_character_
  bar_fill <- if (!is.na(group_label_chr) && group_label_chr %in% c("control", "untrt")) plot_palette$control else plot_palette$treatment

  ggplot(pca_tbl, aes(principal_component, proportion_variance)) +
    geom_col(fill = bar_fill, color = "white") +
    geom_line(
      aes(x = seq_along(principal_component), y = cumulative_proportion, group = 1),
      inherit.aes = FALSE,
      linewidth = LINE_WIDTH_BOUNDARY,
      colour = plot_palette$threshold
    ) +
    geom_point(
      aes(x = seq_along(principal_component), y = cumulative_proportion),
      inherit.aes = FALSE,
      size = 1.6,
      colour = plot_palette$threshold
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    labs(
      title = compact_title(paste(dataset_name, "|", pretty_group_label(group_label), "PCA"), width = 42),
      subtitle = NULL,
      x = "Principal component",
      y = "Variance explained"
    ) +
    manuscript_theme() +
    theme(plot.margin = margin(t = 12, r = 14, b = 14, l = 14), legend.position = "none")
}

plot_pca_scatter <- function(pca_fit, dataset_label, group_label) {
  pca_var     <- pca_fit$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)

  cond_key <- ifelse(group_label %in% c("control", "untrt"), "untrt", "trt")

  pca_df <- data.frame(
    Sample           = rownames(pca_fit$x),
    PC1              = pca_fit$x[, 1],
    PC2              = pca_fit$x[, 2],
    Condition        = cond_key,
    stringsAsFactors = FALSE
  )

  ggplot(pca_df, aes(PC1, PC2, label = Sample, shape = Condition, fill = Condition)) +
    geom_hline(yintercept = 0, linewidth = LINE_WIDTH_ZERO, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = LINE_WIDTH_ZERO, linetype = "dashed", colour = "grey70") +
    geom_point(size = 3.0, colour = "white", stroke = POINT_STROKE + 0.15) +
    geom_text_repel(
      size               = 2.0,
      max.overlaps       = 8,
      force              = 1.0,
      box.padding        = 0.22,
      point.padding      = 0.10,
      min.segment.length = 0
    ) +
    scale_shape_manual(
      values = condition_shapes,
      labels = condition_labels,
      name   = "Condition",
      guide  = guide_legend(
        override.aes = list(
          size   = 3.0,
          fill   = unname(condition_fills),
          colour = "white"
        )
      )
    ) +
    scale_fill_manual(
      values = condition_fills,
      labels = condition_labels,
      name   = "Condition",
      guide  = "none"
    ) +
    labs(
      title = pretty_group_label(group_label),
      subtitle = NULL,
      x = paste0("PC1 (", pca_var_per[1], "%)"),
      y = paste0("PC2 (", pca_var_per[2], "%)")
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(
      plot.title  = element_text(margin = margin(b = 4)),
      plot.margin = margin(t = 12, r = 18, b = 16, l = 16)
    ) +
    guides(shape = guide_legend(order = 1))
}

plot_pc1_loading_rank <- function(loading_tbl, cutoff, dataset_label, group_label,
                                  top_n_used = evs_fixed_top_n, cutoff_quantile = NA_real_,
                                  preprocessing_label = "Normalized prior to eigenvector splitting") {
  quantile_label <- if (is.finite(cutoff_quantile)) {
    paste0("Upper-tail quantile = ", signif(cutoff_quantile, 4))
  } else {
    "Upper-tail quantile = NA"
  }

  ggplot(loading_tbl, aes(rank, pc1_loading_abs)) +
    geom_line(linewidth = LINE_WIDTH_BOUNDARY, color = "grey35") +
    geom_hline(yintercept = cutoff, color = plot_palette$threshold, linewidth = LINE_WIDTH_THRESH) +
    annotate(
      "label",
      x     = max(loading_tbl$rank) * 0.70,
      y     = cutoff,
      label = paste0(
        "Top ", top_n_used, " cutoff = ", signif(cutoff, 4), "\n",
        quantile_label
      ),
      fill  = "white",
      color = plot_palette$threshold,
      vjust = -0.7,
      size  = 3.2,
      label.size = 0.15
    ) +
    labs(
      title = compact_title(
        paste(dataset_label, "|", pretty_group_label(group_label), "PC1 loading rank"),
        width = 36
      ),
      subtitle = compact_caption(
        paste0(
          evs_preproc_short(preprocessing_label),
          ". Ranked absolute PC1 loadings. ",
          "The horizontal line marks the EVS cutoff used to define the leading-edge set."
        ),
        width = 72
      ),
      x = "Ranked PAS feature",
      y = "Absolute PC1 loading"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(
      plot.title    = element_text(margin = margin(b = 4)),
      plot.subtitle = element_text(lineheight = 0.98, margin = margin(b = 6)),
      plot.margin   = margin(t = 16, r = 24, b = 18, l = 18)
    )
}

plot_eigenvector_histograms <- function(loading_tbl, cutoff, dataset_label, group_label,
                                        preprocessing_label = "Normalized prior to eigenvector splitting") {
  full_vals <- loading_tbl$pc1_loading_abs
  lead_vals <- loading_tbl$pc1_loading_abs[loading_tbl$pc1_loading_abs >= cutoff]
  rem_vals  <- loading_tbl$pc1_loading_abs[loading_tbl$pc1_loading_abs < cutoff]
  cutoff_tx <- -log(pmax(cutoff * 100, 1e-6))

  make_hist <- function(vals, panel_title) {
    ggplot(data.frame(x = -log(pmax(vals * 100, 1e-6))), aes(x)) +
      geom_histogram(bins = HIST_BINS, fill = plot_palette$histogram, color = HIST_COLOR) +
      geom_vline(xintercept = cutoff_tx, color = plot_palette$threshold, linewidth = LINE_WIDTH_THRESH) +
      labs(
        title = panel_title,
        x = "-log(|PC1 loading| × 100)",
        y = "Count"
      ) +
      manuscript_theme() +
      theme(
        plot.title  = element_text(size = base_theme_size, hjust = 0.5),
        plot.margin = margin(t = 8, r = 10, b = 10, l = 10)
      )
  }

  arrangeGrob(
    make_hist(full_vals, "All sites"),
    make_hist(rem_vals,  "Remainder"),
    make_hist(lead_vals, "Leading-edge"),
    ncol = 3,
    top = textGrob(
      paste(pretty_group_label(group_label), "|", evs_preproc_short(preprocessing_label)),
      gp = gpar(fontface = "bold", cex = 1.00)
    )
  )
}

plot_evs_candidate_curve <- function(fit_obj, comparison_name, group_label) {
  curve_df <- fit_obj$cutoff_curve_df
  if (is.null(curve_df) || !nrow(curve_df)) return(NULL)

  candidate_tbl <- fit_obj$candidate_table
  if (is.null(candidate_tbl)) candidate_tbl <- data.frame(rank_index = numeric(0), selected = logical(0), stringsAsFactors = FALSE)

  if (all(c("smoothed_log_iod_nb", "smoothed_log_cv2_nb") %in% colnames(curve_df))) {
    plot_df <- as.data.frame(curve_df)
    rng_iod <- range(plot_df$smoothed_log_iod_nb, na.rm = TRUE)
    rng_cv2 <- range(plot_df$smoothed_log_cv2_nb, na.rm = TRUE)
    scale01 <- function(x, rng) {
      den <- diff(rng)
      if (!is.finite(den) || den == 0) return(rep(0.5, length(x)))
      (x - rng[1]) / den
    }
    long_df <- dplyr::bind_rows(
      data.frame(rank = plot_df$rank, metric = "IOD", value = scale01(plot_df$smoothed_log_iod_nb, rng_iod)),
      data.frame(rank = plot_df$rank, metric = "CV2", value = scale01(plot_df$smoothed_log_cv2_nb, rng_cv2))
    )

    cand_df <- as.data.frame(candidate_tbl, stringsAsFactors = FALSE)
    if (!"selected" %in% colnames(cand_df)) cand_df$selected <- FALSE
    cand_unselected <- cand_df[!is.na(cand_df$rank_index) & !cand_df$selected, , drop = FALSE]
    cand_selected <- cand_df[!is.na(cand_df$rank_index) & cand_df$selected, , drop = FALSE]

    interval_df <- fit_obj$matched_interval_table
    if (is.null(interval_df)) interval_df <- data.frame()

    p <- ggplot(long_df, aes(rank, value, color = metric))
    if (nrow(interval_df)) {
      p <- p + geom_rect(data = interval_df, aes(xmin = interval_lo, xmax = interval_hi, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey80", alpha = 0.20, color = NA)
    }
    p <- p + geom_line(linewidth = 0.7) +
      geom_vline(data = cand_unselected, aes(xintercept = rank_index), linetype = "dashed", linewidth = 0.40, alpha = 0.45, color = plot_palette$threshold, inherit.aes = FALSE, na.rm = TRUE) +
      geom_vline(data = cand_selected, aes(xintercept = rank_index), linetype = "solid", linewidth = 0.90, alpha = 0.95, color = plot_palette$threshold, inherit.aes = FALSE, na.rm = TRUE) +
      scale_color_manual(values = c(IOD = plot_palette$treatment, CV2 = plot_palette$control)) +
      labs(
        title = pretty_group_label(group_label),
        subtitle = compact_caption(
          paste0(comparison_name, ". Smoothed ranked IOD and CV2 curves. Grey bands = matched IOD/CV2 second-derivative candidate intervals. Dashed lines = tested wave-scored ranks. Solid line = selected cutoff."),
          width = 76
        ),
        x = "Ranked PAS feature",
        y = "Scaled smoothed curve",
        color = NULL
      ) +
      manuscript_theme() +
      theme(plot.margin = margin(t = 12, r = 14, b = 14, l = 14), legend.position = "top")
    return(p)
  }

  plot_df <- as.data.frame(curve_df)
  cand_df <- dplyr::left_join(candidate_tbl, plot_df[, c("rank", "smoothed_log_nb_variance")], by = c("rank_index" = "rank"))
  cand_df$selected <- as.logical(cand_df$selected)
  cand_unselected <- cand_df[!is.na(cand_df$rank_index) & !cand_df$selected, , drop = FALSE]
  cand_selected   <- cand_df[!is.na(cand_df$rank_index) &  cand_df$selected, , drop = FALSE]

  ggplot(plot_df, aes(rank, smoothed_log_nb_variance)) +
    geom_line(linewidth = LINE_WIDTH_BOUNDARY, color = "grey35") +
    geom_vline(data = cand_unselected, aes(xintercept = rank_index), linetype = "dashed", linewidth = 0.45, alpha = 0.55, color = plot_palette$threshold, inherit.aes = FALSE, na.rm = TRUE) +
    geom_vline(data = cand_selected, aes(xintercept = rank_index), linetype = "solid", linewidth = 0.9, alpha = 0.95, color = plot_palette$threshold, inherit.aes = FALSE, na.rm = TRUE) +
    geom_point(data = cand_unselected, aes(rank_index, smoothed_log_nb_variance), shape = 21, size = 2.2, stroke = 0.35, fill = "white", color = plot_palette$threshold, inherit.aes = FALSE, na.rm = TRUE) +
    geom_point(data = cand_selected, aes(rank_index, smoothed_log_nb_variance), shape = 24, size = 2.2, stroke = 0.35, fill = "white", color = plot_palette$threshold, inherit.aes = FALSE, na.rm = TRUE) +
    labs(
      title = pretty_group_label(group_label),
      subtitle = compact_caption(paste0(comparison_name, ". Smoothed log10 EVS curve profile over ranked absolute PC1 loadings. Dashed lines = tested candidate cutoffs. Solid line = selected cutoff."), width = 74),
      x = "Ranked PAS feature",
      y = if (!is.null(plot_df$curve_label) && length(unique(plot_df$curve_label)) == 1) unique(plot_df$curve_label) else "Smoothed log10 EVS curve"
    ) +
    manuscript_theme() +
    theme(plot.margin = margin(t = 12, r = 14, b = 14, l = 14))
}

plot_primary_evs_candidate_panel <- function(evs, comparison_name) {
  p1 <- plot_evs_candidate_curve(evs$primary_trt_fit, comparison_name, "treatment")
  p2 <- plot_evs_candidate_curve(evs$primary_untrt_fit, comparison_name, "control")
  if (is.null(p1) || is.null(p2)) return(NULL)

  arrangeGrob(
    p1, p2,
    ncol = 2,
    top = textGrob(
      paste0(comparison_name, " | Primary EVS candidate cutoffs"),
      gp = gpar(fontface = "bold", cex = 1.02)
    ),
    bottom = textGrob(
      "The current manuscript workflow uses local Fourier summaries of ranked IOD and CV2 across percentile windows in treatment and control separately, combines those local treatment and control scores into one comparison-level transition score, and selects one shared cutoff rank for the comparison. This legacy curve panel is retained only as a supplemental descriptive diagnostic.",
      gp = gpar(cex = 0.86)
    )
  )
}

plot_combined_pca_scatter_panel <- function(evs, comparison_name) {
  p1 <- plot_pca_scatter(evs$fit_trt$pca_fit, comparison_name, "treatment")
  p2 <- plot_pca_scatter(evs$fit_untrt$pca_fit, comparison_name, "control")
  p3 <- plot_pca_scatter(evs$fit_trt_raw$pca_fit, comparison_name, "treatment")
  p4 <- plot_pca_scatter(evs$fit_untrt_raw$pca_fit, comparison_name, "control")

  arrangeGrob(
    p1, p2, p3, p4,
    ncol = 2,
    top = textGrob(
      paste0(comparison_name, " | EVS PCA"),
      gp = gpar(fontface = "bold", cex = 1.02)
    ),
    bottom = textGrob(
      "Top: normalized before EVS. Bottom: raw before EVS. Left: treatment. Right: control.",
      gp = gpar(cex = 0.86)
    )
  )
}

plot_combined_pca_variance_panel <- function(evs, comparison_name) {
  p1 <- plot_pca_variance_profile(evs$fit_trt$pca_fit, comparison_name, evs$fit_trt$preprocessing_label, "treatment")
  p2 <- plot_pca_variance_profile(evs$fit_untrt$pca_fit, comparison_name, evs$fit_untrt$preprocessing_label, "control")
  p3 <- plot_pca_variance_profile(evs$fit_trt_raw$pca_fit, comparison_name, evs$fit_trt_raw$preprocessing_label, "treatment")
  p4 <- plot_pca_variance_profile(evs$fit_untrt_raw$pca_fit, comparison_name, evs$fit_untrt_raw$preprocessing_label, "control")

  arrangeGrob(
    p1, p2, p3, p4,
    ncol = 2,
    top = textGrob(
      paste0(comparison_name, " | EVS PCA variance profiles"),
      gp = gpar(fontface = "bold", cex = 1.02)
    ),
    bottom = textGrob(
      "Top: normalized before EVS. Bottom: raw before EVS. Left: treatment. Right: control.",
      gp = gpar(cex = 0.86)
    )
  )
}

# =============================================================================
# SECTION 3 OF 4
# DESEQ2, HBFSS, VOLCANOES, DISPERSION, AND SUMMARY PLOTS
# =============================================================================

# -----------------------------------------------------------------------------
# PATCHED CORE DATASET ANALYSIS
# -----------------------------------------------------------------------------

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
  design_formula <- make_design_formula()

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = coldata,
    design    = design_formula
  )

  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, betaPrior = FALSE)

  res <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    alpha = alpha_level,
    independentFiltering = FALSE
  )

  res_strong <- results(
    dds,
    contrast      = c("condition", "trt", "untrt"),
    lfcThreshold  = lfc_boundary,
    altHypothesis = "greaterAbs",
    independentFiltering = FALSE
  )

  res_weak <- results(
    dds,
    contrast      = c("condition", "trt", "untrt"),
    lfcThreshold  = lfc_boundary,
    altHypothesis = "lessAbs",
    independentFiltering = FALSE
  )

  res_all_df <- as.data.frame(res, stringsAsFactors = FALSE)
  res_all_df$feature_id <- as.character(rownames(res_all_df))

  valid_stat <- is.finite(res_all_df$stat) & !is.na(res_all_df$stat)
  stat_vec   <- as.numeric(res_all_df$stat[valid_stat])
  stat_vec   <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]

  stat_mean <- mean(stat_vec, na.rm = TRUE)
  message(
    sprintf(
      "[%s] Mean Wald stat sent to fdrtool: %.4f  (computed only on finite DESeq2 Wald statistics)",
      dataset_name, stat_mean
    )
  )

  fdr_fit <- run_empirical_null_fdrtool(stat_vec, dataset_name = dataset_name)

  res_df  <- res_all_df
  n_valid <- sum(valid_stat)

  if (length(fdr_fit$pval) != n_valid || length(fdr_fit$qval) != n_valid || length(fdr_fit$lfdr) != n_valid) {
    stop(
      sprintf(
        "[%s] fdrtool output length mismatch: n_valid=%d, length(pval)=%d, length(qval)=%d, length(lfdr)=%d.",
        dataset_name, n_valid, length(fdr_fit$pval), length(fdr_fit$qval), length(fdr_fit$lfdr)
      ),
      call. = FALSE
    )
  }

  res_df$empirical_p <- NA_real_
  res_df$empirical_q <- NA_real_
  res_df$lfdr        <- NA_real_

  res_df$empirical_p[valid_stat] <- as.numeric(fdr_fit$pval)
  res_df$empirical_q[valid_stat] <- as.numeric(fdr_fit$qval)
  res_df$lfdr[valid_stat]        <- as.numeric(fdr_fit$lfdr)

  coef_name <- get_condition_coef(dds)

  message(sprintf("[%s] Starting lfcShrink(type = \"%s\")", dataset_name, lfc_shrink_type))
  shr <- tryCatch(
    {
      if (identical(lfc_shrink_type, "apeglm")) {
        suppressMessages(suppressWarnings(
          lfcShrink(
            dds,
            coef = coef_name,
            type = "apeglm",
            res = res,
            apeMethod = lfc_shrink_apeglm_method
          )
        ))
      } else {
        suppressMessages(suppressWarnings(
          lfcShrink(
            dds,
            coef = coef_name,
            type = "normal"
          )
        ))
      }
    },
    error = function(e) {
      message(sprintf("[%s] lfcShrink(type = \"%s\") failed: %s", dataset_name, lfc_shrink_type, conditionMessage(e)))
      message(sprintf("[%s] Falling back to lfcShrink(type = \"normal\")", dataset_name))
      suppressMessages(suppressWarnings(
        lfcShrink(
          dds,
          coef = coef_name,
          type = "normal"
        )
      ))
    }
  )
  message(sprintf("[%s] Finished lfcShrink()", dataset_name))
  shr_df    <- as.data.frame(shr, stringsAsFactors = FALSE)
  shr_df$feature_id <- as.character(rownames(shr_df))

  res_df <- dplyr::left_join(
    res_df,
    shr_df[, c("feature_id", "log2FoldChange")],
    by     = "feature_id",
    suffix = c("", "_shrunk")
  )

  colnames(res_df)[colnames(res_df) == "log2FoldChange_shrunk"] <- "lfc_shrunk"

  hc_p_threshold_dataset <- safe_hc_thresh(res_df$empirical_p, dataset_name = dataset_name)

  empirical_p_floored <- ifelse(
    is.na(res_df$empirical_p),
    NA_real_,
    pmax(res_df$empirical_p, 1e-300)
  )

  res_df$HBFSS <- abs(res_df$lfc_shrunk * log10(empirical_p_floored))

  if (is.na(hc_p_threshold_dataset)) {
    hbfss_threshold_dataset  <- NA_real_
    res_df$HBFSS_significant <- FALSE
    res_df$hc_empirical_pass <- FALSE
  } else {
    hbfss_threshold_dataset  <- abs(log10(hc_p_threshold_dataset)) * lfc_boundary
    res_df$hc_empirical_pass <- !is.na(res_df$empirical_p) & (res_df$empirical_p <= hc_p_threshold_dataset)
    res_df$HBFSS_significant <- ifelse(
      is.na(res_df$HBFSS),
      FALSE,
      res_df$HBFSS >= hbfss_threshold_dataset
    )
  }

  res_df$regulation_direction <- ifelse(
    is.na(res_df$lfc_shrunk),
    NA_character_,
    ifelse(
      res_df$lfc_shrunk > 0, "upregulated",
      ifelse(res_df$lfc_shrunk < 0, "downregulated", "no_change")
    )
  )

  res_df$raw_lfc_pass    <- !is.na(res_df$log2FoldChange) & (abs(res_df$log2FoldChange) >= lfc_boundary)
  res_df$shrunk_lfc_pass <- !is.na(res_df$lfc_shrunk)     & (abs(res_df$lfc_shrunk)     >= lfc_boundary)

  res_df$wald_pvalue <- res_df$pvalue
  res_df$neglog10_wald_pvalue <- safe_neglog10(res_df$wald_pvalue)

  res_strong_df            <- as.data.frame(res_strong)
  res_strong_df$feature_id <- as.character(rownames(res_strong_df))

  res_weak_df            <- as.data.frame(res_weak)
  res_weak_df$feature_id <- as.character(rownames(res_weak_df))

  res_df <- dplyr::left_join(
    res_df,
    res_strong_df[, c("feature_id", "padj")],
    by     = "feature_id",
    suffix = c("", "_strong")
  )

  res_df <- dplyr::left_join(
    res_df,
    res_weak_df[, c("feature_id", "padj")],
    by     = "feature_id",
    suffix = c("", "_weak")
  )

  colnames(res_df)[colnames(res_df) == "padj_strong"] <- "padj_strong_effect"
  colnames(res_df)[colnames(res_df) == "padj_weak"]   <- "padj_weak_effect"

  res_df$resGA_padj <- res_df$padj_strong_effect
  res_df$resLA_padj <- res_df$padj_weak_effect

  res_df$deseq2_strong_call <- !is.na(res_df$resGA_padj) &
    (res_df$resGA_padj < alpha_level) &
    res_df$shrunk_lfc_pass

  res_df$deseq2_weak_call_raw <- !is.na(res_df$resLA_padj) &
    (res_df$resLA_padj < alpha_level) &
    !res_df$shrunk_lfc_pass

  res_df$deseq2_weak_call <- res_df$deseq2_weak_call_raw & res_df$hc_empirical_pass

  res_df$deseq2_standard_call <- !is.na(res_df$padj) &
    (res_df$padj < alpha_level) &
    res_df$shrunk_lfc_pass

  res_df$standard_significant <- res_df$deseq2_strong_call | res_df$deseq2_weak_call | res_df$deseq2_standard_call
  res_df$HBFSS_only_call <- res_df$HBFSS_significant & !res_df$standard_significant
  res_df$overlap_call    <- res_df$HBFSS_significant & (res_df$deseq2_strong_call | res_df$deseq2_standard_call)

  res_df$effect_class <- dplyr::case_when(
    res_df$deseq2_strong_call ~ "strong_effect",
    res_df$deseq2_weak_call   ~ "weak_effect",
    TRUE                      ~ "intermediate"
  )

  base_mean_vec  <- res_df$baseMean[!is.na(res_df$baseMean)]
  norm_counts    <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts$feature_id <- as.character(rownames(norm_counts))

  mm            <- as.data.frame(mcols(dds))
  mm$feature_id <- as.character(rownames(mm))

  disp_cols_available <- intersect(
    c("feature_id", "dispGeneEst", "dispFit", "dispersion", "dispIter", "baseMean", "dispOutlier"),
    colnames(mm)
  )

  disp_df <- mm[, disp_cols_available, drop = FALSE]

  if ("baseMean" %in% colnames(disp_df)) {
    disp_df <- disp_df[, setdiff(colnames(disp_df), "baseMean"), drop = FALSE]
  }

  annot_df$feature_id  <- as.character(annot_df$feature_id)
  annot_df$gene_symbol <- as.character(annot_df$gene_symbol)

  annot_df <- annot_df %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(is.na(gene_symbol), "", trimws(gene_symbol))) %>%
    dplyr::arrange(feature_id, dplyr::desc(gene_symbol != ""), gene_symbol) %>%
    dplyr::distinct(feature_id, .keep_all = TRUE) %>%
    dplyr::mutate(gene_symbol = dplyr::na_if(gene_symbol, ""))

  res_df$feature_id      <- as.character(res_df$feature_id)
  norm_counts$feature_id <- as.character(norm_counts$feature_id)
  disp_df$feature_id     <- as.character(disp_df$feature_id)

  norm_counts <- norm_counts[!duplicated(norm_counts$feature_id), , drop = FALSE]
  disp_df     <- disp_df[!duplicated(disp_df$feature_id), , drop = FALSE]

  final_df <- res_df %>%
    dplyr::left_join(annot_df,    by = "feature_id") %>%
    dplyr::left_join(norm_counts, by = "feature_id") %>%
    dplyr::left_join(disp_df,     by = "feature_id")

  final_df$neglog10_padj        <- safe_neglog10(final_df$padj)
  final_df$neglog10_empirical_p <- safe_neglog10(final_df$empirical_p)

  final_df$dataset_name            <- dataset_name
  final_df$hc_p_threshold_dataset  <- hc_p_threshold_dataset
  final_df$hbfss_threshold_dataset <- hbfss_threshold_dataset

  preferred_cols <- c(
    "dataset_name", "feature_id", "gene_symbol", "baseMean",
    "lfc_shrunk", "regulation_direction",
    "lfdr",
    "pvalue", "wald_pvalue", "padj", "empirical_p", "empirical_q",
    "HBFSS", "hc_p_threshold_dataset", "hbfss_threshold_dataset",
    "resLA_padj", "resGA_padj", "hc_empirical_pass",
    "deseq2_standard_call", "deseq2_strong_call", "deseq2_weak_call_raw", "deseq2_weak_call",
    "standard_significant", "HBFSS_significant", "effect_class"
  )

  final_df <- final_df[, c(intersect(preferred_cols, names(final_df)),
                           setdiff(names(final_df), preferred_cols)), drop = FALSE]

  list(
    dds             = dds,
    results         = final_df,
    base_mean_vec   = base_mean_vec,
    hc_p_threshold  = hc_p_threshold_dataset,
    hbfss_threshold = hbfss_threshold_dataset
  )
}

# -----------------------------------------------------------------------------
# COLOUR / SHAPE SCALES
# -----------------------------------------------------------------------------

method_call_colors <- c(
  "Neither"                   = "grey70",
  "DC2 weak effect"           = plot_palette$weak,
  "DC2 strong effect"         = plot_palette$deseq2,
  "DC2 standard significance" = "#8C2D04",
  "HBFSS only"                = plot_palette$hbfss,
  "Overlap"                   = plot_palette$overlap
)

method_call_shapes <- c(
  "Neither"                   = 21,
  "DC2 weak effect"           = 22,
  "DC2 strong effect"         = 24,
  "DC2 standard significance" = 23,
  "HBFSS only"                = 23,
  "Overlap"                   = 25
)

build_plot_specific_volcano_classes <- function(df, plot_type = c("standard", "hbfss")) {
  plot_type <- match.arg(plot_type)
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  required_cols <- c(
    "gene_symbol", "lfc_shrunk",
    "wald_pvalue", "neglog10_wald_pvalue",
    "empirical_p", "neglog10_empirical_p",
    "deseq2_standard_call", "deseq2_strong_call", "deseq2_weak_call_raw",
    "HBFSS_significant", "HBFSS_only_call", "overlap_call",
    "hc_p_threshold_dataset", "hbfss_threshold_dataset",
    "HBFSS"
  )
  assert_required_columns(df, required_cols, object_name = "volcano input df")

  df$has_valid_gene_symbol <- !is.na(df$gene_symbol) &
    grepl("^[A-Za-z0-9._-]+$", trimws(df$gene_symbol))
  df$gene_symbol_plot <- ifelse(df$has_valid_gene_symbol, trimws(df$gene_symbol), NA_character_)

  strong_visible   <- !is.na(df$deseq2_strong_call)   & df$deseq2_strong_call
  weak_visible     <- !is.na(df$deseq2_weak_call_raw) & df$deseq2_weak_call_raw
  standard_visible <- !is.na(df$deseq2_standard_call) & df$deseq2_standard_call
  hbfss_visible    <- !is.na(df$HBFSS_only_call)      & df$HBFSS_only_call
  overlap_visible  <- !is.na(df$overlap_call)         & df$overlap_call

  df$effect_color_class <- dplyr::case_when(
    overlap_visible ~ "Overlap",
    strong_visible ~ "Strong effect",
    weak_visible ~ "Weak effect",
    standard_visible ~ "Standard significance",
    hbfss_visible ~ "HBFSS only",
    TRUE ~ "Intermediate effect"
  )

  df$effect_color_class <- factor(
    df$effect_color_class,
    levels = c("Strong effect", "Standard significance", "Intermediate effect", "Weak effect", "HBFSS only", "Overlap")
  )

  df$method_call_class <- dplyr::case_when(
    overlap_visible  ~ "Overlap",
    strong_visible   ~ "DC2 strong effect",
    weak_visible     ~ "DC2 weak effect",
    standard_visible ~ "DC2 standard significance",
    hbfss_visible    ~ "HBFSS only",
    TRUE             ~ "Neither"
  )

  df$method_call_class <- factor(
    df$method_call_class,
    levels = c("Neither", "DC2 weak effect", "DC2 strong effect", "DC2 standard significance", "HBFSS only", "Overlap")
  )

  df
}

select_volcano_labels <- function(df, y_col = "neglog10_empirical_p", n_labels = 10) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!nrow(df)) return(df[0, , drop = FALSE])

  df <- df[df$has_valid_gene_symbol, , drop = FALSE]
  if (!nrow(df)) return(df[0, , drop = FALSE])

  df <- df[df$method_call_class %in% c("Overlap", "DC2 strong effect", "DC2 weak effect", "DC2 standard significance", "HBFSS only"), , drop = FALSE]
  if (!nrow(df)) return(df[0, , drop = FALSE])

  df$label_priority <- dplyr::case_when(
    df$method_call_class == "Overlap"                   ~ 1,
    df$method_call_class == "DC2 strong effect"         ~ 2,
    df$method_call_class == "DC2 weak effect"           ~ 3,
    df$method_call_class == "DC2 standard significance" ~ 4,
    df$method_call_class == "HBFSS only"                ~ 5,
    TRUE ~ 9
  )

  metric_y <- suppressWarnings(as.numeric(df[[y_col]]))
  metric_y[!is.finite(metric_y)] <- -Inf

  ord <- order(df$label_priority, -metric_y, -abs(df$lfc_shrunk), na.last = TRUE)
  df <- df[ord, , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol_plot), , drop = FALSE]

  df[seq_len(min(n_labels, nrow(df))), , drop = FALSE]
}

compute_volcano_axis_limits <- function(df, x_col = "lfc_shrunk", y_col = "neglog10_empirical_p") {

  x <- suppressWarnings(as.numeric(df[[x_col]]))
  y <- suppressWarnings(as.numeric(df[[y_col]]))

  x <- x[is.finite(x) & !is.na(x)]
  y <- y[is.finite(y) & !is.na(y)]

  if (!length(x)) x <- c(-2, 2)
  if (!length(y)) y <- c(0, 5)

  x_q <- stats::quantile(abs(x), probs = 0.995, na.rm = TRUE, type = 7)
  x_lim <- max(lfc_boundary + 0.5, as.numeric(x_q))
  x_lim <- min(x_lim, max(abs(x), na.rm = TRUE))
  x_lim <- max(2.5, x_lim)

  y_q <- stats::quantile(y, probs = 0.995, na.rm = TRUE, type = 7)
  y_lim <- max(3, as.numeric(y_q) * 1.08)

  list(
    x = c(-x_lim, x_lim),
    y = c(0, y_lim)
  )
}

volcano_count_caption <- function(df) {
  method_counts <- table(factor(df$method_call_class, levels = levels(df$method_call_class)))
  paste0(
    "Neither=", method_counts["Neither"],
    " | DC2 weak=", method_counts["DC2 weak effect"],
    " | DC2 strong=", method_counts["DC2 strong effect"],
    " | DC2 standard=", method_counts["DC2 standard significance"],
    " | HBFSS=", method_counts["HBFSS only"],
    " | Overlap=", method_counts["Overlap"]
  )
}

.volcano_base_layers <- function(show_legend = TRUE) {
  list(
    geom_vline(
      xintercept = c(-lfc_boundary, lfc_boundary),
      linetype   = "dashed",
      linewidth  = LINE_WIDTH_BOUNDARY,
      colour     = plot_palette$threshold
    ),
    geom_vline(
      xintercept = 0,
      linetype   = "solid",
      linewidth  = LINE_WIDTH_ZERO,
      colour     = "grey45"
    ),
    scale_fill_manual(values = method_call_colors, drop = FALSE, name = "Interpretive tier"),
    scale_color_manual(values = method_call_colors, drop = FALSE, name = "Interpretive tier"),
    scale_shape_manual(values = method_call_shapes, drop = FALSE, name = "Interpretive tier"),
    guides(
      fill = if (show_legend) guide_legend(
        order = 1,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          size = 3.8,
          stroke = 0.8,
          alpha = 1,
          shape = unname(method_call_shapes),
          fill = unname(method_call_colors),
          colour = unname(method_call_colors)
        )
      ) else "none",
      color = "none",
      shape = "none"
    )
  )
}

volcano_label_layer <- function(lab_df) {
  if (!nrow(lab_df)) return(NULL)
  ggrepel::geom_text_repel(
    data               = lab_df,
    aes(label = gene_symbol_plot),
    size               = 1.95,
    seed               = 1,
    max.overlaps       = 35,
    force              = 1.25,
    force_pull         = 0.5,
    box.padding        = 0.38,
    point.padding      = 0.22,
    min.segment.length = 0,
    segment.alpha      = 0.6,
    segment.size       = 0.22,
    bg.color           = "white",
    bg.r               = 0.09
  )
}

volcano_point_layer <- function(df) {
  geom_point(
    data = df,
    aes(
      x = lfc_shrunk,
      y = plot_y,
      fill = method_call_class,
      color = method_call_class,
      shape = method_call_class
    ),
    alpha = POINT_ALPHA_PRIMARY,
    size = POINT_SIZE_PRIMARY + 0.45,
    stroke = POINT_STROKE + 0.14,
    na.rm = TRUE
  )
}

extract_plot_legend <- function(p) {
  g <- ggplotGrob(p)
  guide_idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
  if (!length(guide_idx)) return(NULL)
  g$grobs[[guide_idx[1]]]
}

add_hbfss_boundary_layer <- function(p, df, annotate_thresholds = TRUE) {
  threshold <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))
  hc_p_threshold <- suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1]))

  if (!is.finite(threshold) || is.na(threshold) || threshold <= 0) {
    return(p)
  }

  axis_limits <- compute_volcano_axis_limits(df, x_col = "lfc_shrunk", y_col = "neglog10_empirical_p")
  max_abs_lfc <- max(abs(axis_limits$x))

  x_abs <- seq(
    from = max(0.08, min(lfc_boundary, max_abs_lfc)),
    to = max_abs_lfc,
    length.out = 400
  )

  y_hyper <- threshold / x_abs

  boundary_df <- data.frame(
    lfc_shrunk = c(-rev(x_abs), x_abs),
    neglog10_empirical_p = c(rev(y_hyper), y_hyper),
    stringsAsFactors = FALSE
  )

  boundary_df <- boundary_df[
    is.finite(boundary_df$neglog10_empirical_p) &
      !is.na(boundary_df$neglog10_empirical_p),
    ,
    drop = FALSE
  ]

  p <- p +
    geom_path(
      data = boundary_df,
      aes(x = lfc_shrunk, y = neglog10_empirical_p),
      inherit.aes = FALSE,
      linetype = "22",
      linewidth = 0.55,
      colour = plot_palette$hbfss,
      alpha = 0.80
    )

  if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < 1) {
    y_floor <- -log10(hc_p_threshold)
    p <- p +
      geom_hline(
        yintercept = y_floor,
        linetype = "22",
        linewidth = 0.55,
        colour = plot_palette$threshold,
        alpha = 0.80
      )
  }

  if (isTRUE(annotate_thresholds)) {
    label_txt <- paste0(
      "HC p = ", ifelse(is.finite(hc_p_threshold), signif(hc_p_threshold, 4), "NA"),
      "
HBFSS = ", signif(threshold, 4), " / |LFC|"
    )

    p <- p +
      annotate(
        "label",
        x = axis_limits$x[1] * 0.96,
        y = axis_limits$y[2] * 0.96,
        label = label_txt,
        hjust = 0,
        vjust = 1,
        fill = "white",
        colour = "grey25",
        size = 2.3,
        label.size = 0.12
      )
  }

  p
}

# -----------------------------------------------------------------------------
# PATCHED VOLCANO PLOTS
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# PATCHED VOLCANO PLOTS
# -----------------------------------------------------------------------------

plot_standard_volcano <- function(df, dataset_name, show_legend = TRUE) {
  df <- build_plot_specific_volcano_classes(df, plot_type = "standard")
  df$plot_y <- df$neglog10_wald_pvalue
  lab_df <- select_volcano_labels(df, y_col = "neglog10_wald_pvalue", n_labels = 20)

  p <- ggplot(df, aes(lfc_shrunk, plot_y)) +
    volcano_point_layer(df) +
    .volcano_base_layers(show_legend = show_legend) +
    geom_hline(
      yintercept = -log10(alpha_level),
      linetype   = "dashed",
      linewidth  = LINE_WIDTH_BOUNDARY,
      colour     = plot_palette$threshold
    ) +
    labs(
      title   = pretty_dataset_label(dataset_name),
      x       = "Shrunken log2 fold change (β̂shrunk)",
      y       = expression(-log[10]("Wald p value")),
      caption = compact_caption(paste0(
        "DC2 uses Wald-test p values directly. ",
        "Strong and weak DESeq2 calls are colored only when they clear the displayed Wald p threshold on this panel. ",
        volcano_count_caption(df),
        "  |  LFC lines = ±", lfc_boundary
      ))
    ) +
    coord_cartesian(clip = "off") +
    plot_expand_xy() +
    manuscript_theme()

  if (nrow(lab_df) > 0) p <- p + volcano_label_layer(lab_df)
  p
}

plot_hbfss_volcano <- function(df, dataset_name, show_legend = TRUE) {
  df <- build_plot_specific_volcano_classes(df, plot_type = "hbfss")
  df$plot_y <- df$neglog10_empirical_p
  lab_df <- select_volcano_labels(df, y_col = "neglog10_empirical_p", n_labels = 20)

  p <- ggplot(df, aes(lfc_shrunk, plot_y)) +
    volcano_point_layer(df) +
    .volcano_base_layers(show_legend = show_legend)

  p <- add_hbfss_boundary_layer(p, df) +
    labs(
      title   = pretty_dataset_label(dataset_name),
      x       = "Shrunken log2 fold change (β̂shrunk)",
      y       = expression(-log[10](p[empirical])),
      caption = compact_caption(paste0(
        "HBFSS uses empirical p values from fdrtool. ",
        "Strong, weak, HBFSS-only, and overlap markers are colored only when they clear the displayed empirical thresholds on this panel. ",
        volcano_count_caption(df)
      ))
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    plot_expand_xy()

  if (nrow(lab_df) > 0) p <- p + volcano_label_layer(lab_df)
  p
}

plot_publication_volcano_panel <- function(df, dataset_name, show_legend = FALSE) {
  build <- build_plot_specific_volcano_classes(df, plot_type = "hbfss")
  build$plot_y <- build$neglog10_empirical_p
  top_df <- select_volcano_labels(build, y_col = "neglog10_empirical_p", n_labels = 8)

  lims <- compute_volcano_axis_limits(build, x_col = "lfc_shrunk", y_col = "plot_y")

  p <- ggplot(build, aes(lfc_shrunk, plot_y)) +
    volcano_point_layer(build) +
    .volcano_base_layers(show_legend = show_legend) +
    labs(
      title = pretty_dataset_label(dataset_name),
      x = "Shrunken log2 fold change (β̂shrunk)",
      y = expression(-log[10](p[empirical])),
      caption = compact_caption(
        paste0(
          "Representative HBFSS volcano with DC2 weak, DC2 strong, DC2 standard significance, HBFSS-only, and overlap classes. ",
          volcano_count_caption(build)
        )
      )
    ) +
    coord_cartesian(xlim = lims$x, ylim = lims$y, clip = "off") +
    manuscript_theme() +
    theme(
      legend.position = if (isTRUE(show_legend)) "bottom" else "none",
      legend.direction = "horizontal",
      plot.caption = element_text(size = base_theme_size - 3, lineheight = 0.95),
      plot.margin = margin(t = 10, r = 12, b = 12, l = 12)
    )

  p <- add_hbfss_boundary_layer(p, build, annotate_thresholds = TRUE)

  if (nrow(top_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top_df,
      aes(label = gene_symbol_plot),
      size = 1.85,
      seed = 1,
      max.overlaps = 18,
      force = 1.0,
      force_pull = 0.4,
      box.padding = 0.25,
      point.padding = 0.14,
      min.segment.length = 0,
      segment.alpha = 0.5,
      segment.size = 0.18
    )
  }

  p
}

dispersion_residual_section <- function(df, dataset_name, fig_subdir) {
  if (!all(c("dispGeneEst", "dispFit", "dispersion") %in% colnames(df))) return(NULL)

  out <- data.frame(
    feature_id                = df$feature_id,
    dispGeneEst               = df$dispGeneEst,
    dispFit                   = df$dispFit,
    dispersion                = df$dispersion,
    residual_fit_minus_final  = df$dispFit - df$dispersion,
    residual_fit_minus_gene   = df$dispFit - df$dispGeneEst,
    residual_final_minus_gene = df$dispersion - df$dispGeneEst,
    stringsAsFactors          = FALSE
  )

  make_hist <- function(vals, panel_title) {
    ggplot(data.frame(x = vals), aes(x)) +
      geom_histogram(bins = HIST_BINS, fill = plot_palette$histogram, color = HIST_COLOR) +
      labs(title = panel_title, x = "Residual", y = "Count") +
      manuscript_theme()
  }

  panel <- arrangeGrob(
    make_hist(out$residual_fit_minus_final,  "dispFit - final dispersion"),
    make_hist(out$residual_fit_minus_gene,   "dispFit - gene-wise dispersion"),
    make_hist(out$residual_final_minus_gene, "final dispersion - gene-wise dispersion"),
    ncol = 3,
    top  = textGrob(
      paste(compact_title(pretty_dataset_label(dataset_name)), "Dispersion residuals"),
      gp = gpar(fontface = "bold", cex = 1.2)
    )
  )

  save_grob(
    panel,
    file.path(fig_subdir, paste0(dataset_name, "_dispersion_residuals.png")),
    width  = 16.8,
    height = 5.8,
    dpi    = figure_dpi
  )

  out
}

plot_hbfss_distribution <- function(df, dataset_name) {
  if (is.null(df) || !nrow(df) || !"HBFSS" %in% colnames(df)) return(NULL)

  plot_df <- df[is.finite(df$HBFSS) & !is.na(df$HBFSS), , drop = FALSE]
  if (!nrow(plot_df)) return(NULL)

  ggplot(plot_df, aes(HBFSS)) +
    geom_histogram(bins = HIST_BINS, fill = plot_palette$hbfss, color = HIST_COLOR, alpha = 0.90) +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "| HBFSS distribution")),
      subtitle = NULL,
      x = "HBFSS",
      y = "PAS features"
    ) +
    manuscript_theme()
}

plot_empirical_histogram_for_panel <- function(df, dataset_name, hc_p_threshold) {
  p <- ggplot(df, aes(empirical_p)) +
    geom_histogram(bins = HIST_BINS, fill = plot_palette$histogram, color = HIST_COLOR) +
    labs(
      title    = compact_title(pretty_dataset_label(dataset_name)),
      subtitle = if (is.finite(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < max_usable_hc_p_threshold) {
        paste0("HC = ", signif(hc_p_threshold, 4))
      } else {
        "HC unavailable"
      },
      x = "Empirical-null p-value",
      y = "Count"
    ) +
    manuscript_theme()

  if (is.finite(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < max_usable_hc_p_threshold) {
    p <- p + geom_vline(
      xintercept = hc_p_threshold,
      color      = plot_palette$threshold,
      linewidth  = LINE_WIDTH_THRESH
    )
  }

  p
}

compute_dataset_pca_plot <- function(count_df, coldata, dataset_name,
                                     preprocessing = c("normalized", "raw_counts"),
                                     precomputed_matrix = NULL) {
  preprocessing <- match.arg(preprocessing)
  design_formula <- make_design_formula()

  if (!is.null(precomputed_matrix)) {
    x <- as.matrix(precomputed_matrix)
  } else if (preprocessing == "normalized") {
    dds <- DESeqDataSetFromMatrix(countData = count_df, colData = coldata, design = design_formula)
    dds <- dds[rowSums(counts(dds)) > 0, ]
    dds <- estimateSizeFactors(dds)
    x   <- counts(dds, normalized = TRUE)
  } else {
    x   <- as.matrix(count_df)
  }

  pca_fit     <- prcomp(t(x), scale. = FALSE, rank. = 2)
  pca_var     <- pca_fit$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)

  pca_df <- data.frame(
    Sample           = rownames(pca_fit$x),
    PC1              = pca_fit$x[, 1],
    PC2              = pca_fit$x[, 2],
    Condition        = as.character(coldata[rownames(pca_fit$x), "condition"]),
    stringsAsFactors = FALSE
  )

  ggplot(pca_df, aes(PC1, PC2, label = Sample, shape = Condition, fill = Condition)) +
    geom_hline(yintercept = 0, linewidth = LINE_WIDTH_ZERO, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = LINE_WIDTH_ZERO, linetype = "dashed", colour = "grey70") +
    geom_point(size = 3, colour = "white", stroke = POINT_STROKE + 0.15) +
    geom_text_repel(
      size               = 2.0,
      max.overlaps       = 8,
      force              = 1.0,
      box.padding        = 0.22,
      point.padding      = 0.10,
      min.segment.length = 0
    ) +
    scale_shape_manual(
      values = condition_shapes,
      labels = condition_labels,
      name   = "Condition",
      guide  = guide_legend(
        override.aes = list(
          size   = 3.0,
          fill   = unname(condition_fills),
          colour = "white"
        )
      )
    ) +
    scale_fill_manual(
      values = condition_fills,
      labels = condition_labels,
      name   = "Condition",
      guide  = "none"
    ) +
    labs(
      title    = compact_title(pretty_dataset_label(dataset_name)),
      subtitle = NULL,
      x = paste0("PC1 (", pca_var_per[1], "%)"),
      y = paste0("PC2 (", pca_var_per[2], "%)")
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 12, r = 18, b = 16, l = 16)
    ) +
    guides(shape = guide_legend(order = 1))
}

plot_dispersion_panel_for_dataset <- function(df, dataset_name) {
  df <- build_reviewer_volcano_classes(df)

  required_cols <- c("baseMean", "dispersion", "method_call_class")
  if (!all(required_cols %in% colnames(df))) return(NULL)

  plot_df <- data.frame(
    baseMean = suppressWarnings(as.numeric(df$baseMean)),
    dispersion = suppressWarnings(as.numeric(df$dispersion)),
    method_call_class = df$method_call_class,
    stringsAsFactors = FALSE
  )

  plot_df <- plot_df[
    is.finite(plot_df$baseMean) & !is.na(plot_df$baseMean) & plot_df$baseMean > 0 &
      is.finite(plot_df$dispersion) & !is.na(plot_df$dispersion) & plot_df$dispersion > 0,
    ,
    drop = FALSE
  ]

  if (!nrow(plot_df)) return(NULL)

  plot_df$method_call_class <- factor(plot_df$method_call_class, levels = levels(df$method_call_class))

  ggplot(plot_df, aes(baseMean, dispersion, color = method_call_class, shape = method_call_class)) +
    geom_point(alpha = POINT_ALPHA_DISP, size = POINT_SIZE_DISP, stroke = POINT_STROKE) +
    scale_x_log10(labels = label_number(accuracy = 0.1)) +
    scale_y_log10(labels = label_number(accuracy = 0.1)) +
    scale_color_manual(values = method_call_colors, drop = FALSE, name = "Interpretive tier") +
    scale_shape_manual(values = method_call_shapes, drop = FALSE, name = "Interpretive tier") +
    labs(
      title    = compact_title(pretty_dataset_label(dataset_name), width = 42),
      subtitle = NULL,
      x        = "baseMean (log10 scale)",
      y        = "Final dispersion (log10 scale)"
    ) +
    manuscript_theme() +
    volcano_guides()
}

save_cross_dataset_comparison_panels <- function(comparison_name, analysis_results, cmp_dir,
                                                 dataset_list, coldata, normalized_dataset_list = NULL) {
  keys_present <- base::intersect(as.character(dataset_key_order), as.character(names(analysis_results)))
  if (length(keys_present) == 0) return(invisible(NULL))

  safe_panel_plot <- function(expr, label) {
    safe_plot_build(expr, paste0(comparison_name, ": ", label))
  }

  make_panel <- function(grob_list, title_text, ncols = length(grob_list)) {
    grob_list <- Filter(Negate(is.null), grob_list)
    if (!length(grob_list)) return(NULL)
    do.call(
      arrangeGrob,
      c(
        grob_list,
        list(
          ncol = ncols,
          top  = textGrob(title_text, gp = gpar(fontface = "bold", cex = 1.10))
        )
      )
    )
  }

  if (isTRUE(export_optional_empirical_p_histograms)) {
    hist_grobs <- lapply(keys_present, function(k) {
      safe_panel_plot(
        plot_empirical_histogram_for_panel(
          analysis_results[[k]]$results,
          analysis_results[[k]]$summary$dataset_name[1],
          analysis_results[[k]]$summary$hc_p_threshold[1]
        ),
        paste0("empirical p histogram: ", k)
      )
    })
    hist_panel <- make_panel(hist_grobs, paste(comparison_name, "| Empirical p"))
    if (!is.null(hist_panel)) {
      save_grob(
        hist_panel,
        file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_empirical_p_histograms.png")),
        width  = 6.4 * length(Filter(Negate(is.null), hist_grobs)),
        height = 5.9
      )
    }
  }

  if (isTRUE(export_optional_hbfss_distributions)) {
    hbfss_dist_grobs <- lapply(keys_present, function(k) {
      safe_panel_plot(
        plot_hbfss_distribution(
          analysis_results[[k]]$results,
          analysis_results[[k]]$summary$dataset_name[1]
        ),
        paste0("HBFSS distribution: ", k)
      )
    })
    hbfss_panel <- make_panel(hbfss_dist_grobs, paste(comparison_name, "| HBFSS distributions"))
    if (!is.null(hbfss_panel)) {
      save_grob(
        hbfss_panel,
        file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_HBFSS_distributions.png")),
        width  = 6.4 * length(Filter(Negate(is.null), hbfss_dist_grobs)),
        height = 5.9
      )
    }
  }

  disp_grobs <- lapply(keys_present, function(k) {
    safe_panel_plot(
      plot_dispersion_panel_for_dataset(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1]
      ),
      paste0("dispersion panel: ", k)
    )
  })
  disp_panel <- make_panel(disp_grobs, paste(comparison_name, "| Dispersion"))
  if (!is.null(disp_panel)) {
    save_grob(
      disp_panel,
      file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_dispersion_panel.png")),
      width  = 6.4 * length(Filter(Negate(is.null), disp_grobs)),
      height = 5.9
    )
  }

  std_plots <- lapply(seq_along(keys_present), function(ii) {
    k <- keys_present[[ii]]
    safe_panel_plot(
      plot_standard_volcano(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1],
        show_legend = FALSE
      ),
      paste0("standard volcano: ", k)
    )
  })
  std_legend_source <- safe_panel_plot(
    plot_standard_volcano(
      analysis_results[[keys_present[1]]]$results,
      analysis_results[[keys_present[1]]]$summary$dataset_name[1],
      show_legend = TRUE
    ),
    "standard volcano legend source"
  )
  std_legend <- if (!is.null(std_legend_source)) extract_plot_legend(std_legend_source) else NULL
  std_grobs <- Filter(Negate(is.null), std_plots)
  if (length(std_grobs)) {
    std_panel <- arrangeGrob(
      grobs = c(std_grobs, list(std_legend)),
      layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
      heights = c(12, 1.8),
      top = textGrob(paste0(comparison_name, " | DESeq2 volcanoes"), gp = gpar(fontface = "bold", cex = 1.05))
    )
    save_grob(
      std_panel,
      file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_standard_volcano_panel.png")),
      width = 18,
      height = 7.8
    )
  }

  hbfss_plots <- lapply(seq_along(keys_present), function(ii) {
    k <- keys_present[[ii]]
    safe_panel_plot(
      plot_hbfss_volcano(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1],
        show_legend = FALSE
      ),
      paste0("HBFSS volcano: ", k)
    )
  })
  hbfss_legend_source <- safe_panel_plot(
    plot_hbfss_volcano(
      analysis_results[[keys_present[1]]]$results,
      analysis_results[[keys_present[1]]]$summary$dataset_name[1],
      show_legend = TRUE
    ),
    "HBFSS volcano legend source"
  )
  hbfss_legend <- if (!is.null(hbfss_legend_source)) extract_plot_legend(hbfss_legend_source) else NULL
  hbfss_grobs <- Filter(Negate(is.null), hbfss_plots)
  if (length(hbfss_grobs)) {
    hbfss_panel <- arrangeGrob(
      grobs = c(hbfss_grobs, list(hbfss_legend)),
      layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
      heights = c(12, 1.8),
      top = textGrob(paste0(comparison_name, " | HBFSS volcanoes"), gp = gpar(fontface = "bold", cex = 1.05))
    )
    save_grob(
      hbfss_panel,
      file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_HBFSS_volcano_panel.png")),
      width = 18,
      height = 7.8
    )
  }

  pub_plots <- lapply(keys_present, function(k) {
    safe_panel_plot(
      plot_publication_volcano_panel(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1],
        show_legend = FALSE
      ),
      paste0("publication volcano: ", k)
    )
  })
  pub_plots <- Filter(Negate(is.null), pub_plots)
  pub_legend_source <- safe_panel_plot(
    plot_publication_volcano_panel(
      analysis_results[[keys_present[1]]]$results,
      analysis_results[[keys_present[1]]]$summary$dataset_name[1],
      show_legend = TRUE
    ),
    "publication volcano legend source"
  )
  pub_legend <- if (!is.null(pub_legend_source)) extract_plot_legend(pub_legend_source) else NULL
  if (length(pub_plots)) {
    pub_panel <- arrangeGrob(
      grobs = c(pub_plots, list(pub_legend)),
      layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
      heights = c(12, 1.8),
      top = textGrob(paste0(comparison_name, " | HBFSS volcanoes"), gp = gpar(fontface = "bold", cex = 1.05))
    )
    save_grob(
      pub_panel,
      file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_publication_volcano_panel.png")),
      width = 18,
      height = 7.8
    )
  }

  pca_norm_grobs <- lapply(keys_present, function(k) {
    safe_panel_plot(
      compute_dataset_pca_plot(
        dataset_list[[k]],
        coldata,
        analysis_results[[k]]$summary$dataset_name[1],
        preprocessing = "normalized",
        precomputed_matrix = if (!is.null(normalized_dataset_list)) normalized_dataset_list[[k]] else NULL
      ),
      paste0("normalized PCA: ", k)
    )
  })

  pca_raw_grobs <- lapply(keys_present, function(k) {
    safe_panel_plot(
      compute_dataset_pca_plot(
        dataset_list[[k]],
        coldata,
        analysis_results[[k]]$summary$dataset_name[1],
        preprocessing = "raw_counts"
      ),
      paste0("raw-count PCA: ", k)
    )
  })

  pca_norm_grobs <- Filter(Negate(is.null), pca_norm_grobs)
  pca_raw_grobs  <- Filter(Negate(is.null), pca_raw_grobs)
  if (length(pca_norm_grobs) > 0 || length(pca_raw_grobs) > 0) {
    pca_panel <- arrangeGrob(
      grobs = c(pca_norm_grobs, pca_raw_grobs),
      ncol  = max(1L, max(length(pca_norm_grobs), length(pca_raw_grobs))),
      top   = textGrob(
        paste(comparison_name, "| Top: normalized before EVS. Bottom: raw before EVS."),
        gp = gpar(fontface = "bold", cex = 1.08)
      )
    )

    save_grob(
      pca_panel,
      file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_PCA_normalized_vs_unnormalized.png")),
      width  = 6.2 * max(1L, max(length(pca_norm_grobs), length(pca_raw_grobs))),
      height = 10.5
    )
  }

  summary_table <- dplyr::bind_rows(lapply(keys_present, function(k) {
    sm <- analysis_results[[k]]$summary
    data.frame(
      comparison_name        = comparison_name,
      dataset_key            = k,
      dataset_label          = unname(dataset_key_labels[k]),
      dataset_name           = sm$dataset_name[1],
      n_features             = sm$n_features[1],
      n_standard_significant = sm$n_standard_significant[1],
      n_HBFSS_significant    = sm$n_HBFSS_significant[1],
      n_overlap              = sum(
        analysis_results[[k]]$results$deseq2_strong_call &
          analysis_results[[k]]$results$HBFSS_significant,
        na.rm = TRUE
      ),
      n_strong_effect        = sm$n_strong_effect[1],
      n_weak_effect          = sm$n_weak_effect[1],
      hc_p_threshold         = sm$hc_p_threshold[1],
      hbfss_threshold        = sm$hbfss_threshold[1],
      evs_fixed_top_n        = evs_fixed_top_n,
      evs_cutoff_mode_main   = evs_cutoff_mode_main,
      evs_primary_preprocessing = "normalized",
      stringsAsFactors       = FALSE
    )
  }))

  save_csv(
    summary_table,
    file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_summary_for_figure_panels.csv"))
  )

  invisible(NULL)
}

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  cmp_dir     <- file.path(output_dir, comparison_name)
  tab_dir     <- file.path(cmp_dir, "tables")
  fig_raw_dir <- file.path(cmp_dir, "figures_raw")
  fig_le_dir  <- file.path(cmp_dir, "figures_leading_edge")
  fig_rem_dir <- file.path(cmp_dir, "figures_remainder")

  for (d in c(cmp_dir, tab_dir, fig_raw_dir, fig_le_dir, fig_rem_dir)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  evs <- build_eigenvector_split(count_matrix, coldata, comparison_name)

  evs_cutoff_summary_out <- evs$evs_cutoff_summary
  evs_cutoff_summary_out$comparison_name <- comparison_name
  save_csv(evs_cutoff_summary_out, file.path(tab_dir, paste0(comparison_name, "_EVS_cutoff_summary.csv")))

  if (!is.null(evs$combined_cutoff_info$trt_wave_obj$wave_map) && nrow(evs$combined_cutoff_info$trt_wave_obj$wave_map)) {
    save_csv(evs$combined_cutoff_info$trt_wave_obj$wave_map, file.path(tab_dir, paste0(comparison_name, "_treatment_local_fourier_map.csv")))
  }
  if (!is.null(evs$combined_cutoff_info$ctrl_wave_obj$wave_map) && nrow(evs$combined_cutoff_info$ctrl_wave_obj$wave_map)) {
    save_csv(evs$combined_cutoff_info$ctrl_wave_obj$wave_map, file.path(tab_dir, paste0(comparison_name, "_control_local_fourier_map.csv")))
  }
  if (!is.null(evs$combined_cutoff_info$combined_wave_map) && nrow(evs$combined_cutoff_info$combined_wave_map)) {
    save_csv(evs$combined_cutoff_info$combined_wave_map, file.path(tab_dir, paste0(comparison_name, "_first_stable_crossing_map.csv")))
  }
  if (!is.null(evs$combined_cutoff_info$crossing_table) && nrow(evs$combined_cutoff_info$crossing_table)) {
    save_csv(evs$combined_cutoff_info$crossing_table, file.path(tab_dir, paste0(comparison_name, "_all_regime_crossings.csv")))
  }
  if (!is.null(evs$combined_cutoff_info$candidate_table) && nrow(evs$combined_cutoff_info$candidate_table)) {
    save_csv(evs$combined_cutoff_info$candidate_table, file.path(tab_dir, paste0(comparison_name, "_selected_regime_crossing.csv")))
  }
  crossing_summary_df <- build_crossing_summary_table(evs$combined_cutoff_info, comparison_name)
  save_csv(crossing_summary_df, file.path(tab_dir, paste0(comparison_name, "_regime_shift_crossing_summary.csv")))
  backup_crossing_summary_df <- build_backup_crossing_summary_table(evs$combined_cutoff_info, comparison_name)
  save_csv(backup_crossing_summary_df, file.path(tab_dir, paste0(comparison_name, "_backup_crossing_summary.csv")))

  if (!is.null(evs$raw_wald_search_info)) {
    save_csv(
      evs$raw_wald_search_info$coarse_table,
      file.path(tab_dir, paste0(comparison_name, "_raw_wald_bracket_search_coarse.csv"))
    )
    save_csv(
      evs$raw_wald_search_info$fine_table,
      file.path(tab_dir, paste0(comparison_name, "_raw_wald_bracket_search_fine.csv"))
    )
    save_csv(
      evs$raw_wald_search_info$selected_row,
      file.path(tab_dir, paste0(comparison_name, "_raw_wald_bracket_search_selected.csv"))
    )
  }

  trt_fourier_plot <- safe_plot_build(
    plot_fourier_wave_map_single(evs$combined_cutoff_info$trt_wave_obj, comparison_name, "treatment"),
    paste0(comparison_name, ": treatment local fourier wave map")
  )
  if (!is.null(trt_fourier_plot)) {
    save_grob(trt_fourier_plot, file.path(cmp_dir, paste0(comparison_name, "_treatment_local_fourier_wave_map.png")), width = 12, height = 10)
  }
  ctrl_fourier_plot <- safe_plot_build(
    plot_fourier_wave_map_single(evs$combined_cutoff_info$ctrl_wave_obj, comparison_name, "control"),
    paste0(comparison_name, ": control local fourier wave map")
  )
  if (!is.null(ctrl_fourier_plot)) {
    save_grob(ctrl_fourier_plot, file.path(cmp_dir, paste0(comparison_name, "_control_local_fourier_wave_map.png")), width = 12, height = 10)
  }
  combined_fourier_plot <- safe_plot_build(
    plot_combined_fourier_wave_map(evs$combined_cutoff_info, comparison_name),
    paste0(comparison_name, ": combined local fourier wave map")
  )
  if (!is.null(combined_fourier_plot)) {
    save_grob(combined_fourier_plot, file.path(cmp_dir, paste0(comparison_name, "_combined_local_fourier_wave_map.png")), width = 12, height = 6)
  }

  combined_score_plot <- safe_plot_build(
    plot_combined_fourier_score(evs$combined_cutoff_info, comparison_name),
    paste0(comparison_name, ": combined local fourier score")
  )
  if (!is.null(combined_score_plot)) {
    save_grob(combined_score_plot, file.path(cmp_dir, paste0(comparison_name, "_combined_local_fourier_score.png")), width = 12, height = 6)
  }

  fourier_summary_panel <- safe_plot_build(
    plot_fourier_summary_panel(evs$combined_cutoff_info, comparison_name),
    paste0(comparison_name, ": fourier summary panel")
  )
  if (!is.null(fourier_summary_panel)) {
    save_grob(
      fourier_summary_panel,
      file.path(cmp_dir, paste0(comparison_name, "_fourier_summary_panel.png")),
      width = 18,
      height = 12
    )
  }

  evs_pca_scatter <- safe_plot_build(
    plot_combined_pca_scatter_panel(evs, comparison_name),
    paste0(comparison_name, ": EVS PCA scatter panel")
  )
  if (!is.null(evs_pca_scatter)) {
    save_grob(
      evs_pca_scatter,
      file.path(cmp_dir, "EVS_PCA_scatter_combined.png"),
      width = 18, height = 13
    )
  }

  if (isTRUE(export_optional_evs_raw_histograms)) {
    trt_hist_grob <- safe_plot_build(
      plot_eigenvector_histograms(
        evs$fit_trt$loading_table,
        evs$fit_trt$cutoff,
        comparison_name,
        "treatment",
        evs$fit_trt$preprocessing_label
      ),
      paste0(comparison_name, ": treatment EVS histogram")
    )

    ctrl_hist_grob <- safe_plot_build(
      plot_eigenvector_histograms(
        evs$fit_untrt$loading_table,
        evs$fit_untrt$cutoff,
        comparison_name,
        "control",
        evs$fit_untrt$preprocessing_label
      ),
      paste0(comparison_name, ": control EVS histogram")
    )

    trt_hist_grob_raw <- safe_plot_build(
      plot_eigenvector_histograms(
        evs$fit_trt_raw$loading_table,
        evs$fit_trt_raw$cutoff,
        comparison_name,
        "treatment",
        evs$fit_trt_raw$preprocessing_label
      ),
      paste0(comparison_name, ": treatment EVS histogram raw")
    )

    ctrl_hist_grob_raw <- safe_plot_build(
      plot_eigenvector_histograms(
        evs$fit_untrt_raw$loading_table,
        evs$fit_untrt_raw$cutoff,
        comparison_name,
        "control",
        evs$fit_untrt_raw$preprocessing_label
      ),
      paste0(comparison_name, ": control EVS histogram raw")
    )

    evs_hist_panel <- safe_plot_build(
      arrangeGrob(
        trt_hist_grob,
        ctrl_hist_grob,
        trt_hist_grob_raw,
        ctrl_hist_grob_raw,
        ncol = 2,
        top = textGrob(
          paste(comparison_name, "| Top: normalized before EVS. Bottom: raw before EVS. Left: treatment. Right: control."),
          gp = gpar(fontface = "bold", cex = 1.05)
        )
      ),
      paste0(comparison_name, ": EVS histogram panel")
    )

    if (!is.null(evs_hist_panel)) {
      save_grob(
        evs_hist_panel,
        file.path(cmp_dir, "EVS_histograms_combined_panel.png"),
        width = 18, height = 11
      )
    }
  }

  save_csv(
    make_pca_summary_table(evs$fit_trt$pca_fit, comparison_name, evs$fit_trt$preprocessing_label, "treatment"),
    file.path(tab_dir, paste0(comparison_name, "_EVS_treatment_PCA_summary_normalized.csv"))
  )
  save_csv(
    make_pca_summary_table(evs$fit_untrt$pca_fit, comparison_name, evs$fit_untrt$preprocessing_label, "control"),
    file.path(tab_dir, paste0(comparison_name, "_EVS_control_PCA_summary_normalized.csv"))
  )
  save_csv(
    make_pca_summary_table(evs$fit_trt_raw$pca_fit, comparison_name, evs$fit_trt_raw$preprocessing_label, "treatment"),
    file.path(tab_dir, paste0(comparison_name, "_EVS_treatment_PCA_summary_raw.csv"))
  )
  save_csv(
    make_pca_summary_table(evs$fit_untrt_raw$pca_fit, comparison_name, evs$fit_untrt_raw$preprocessing_label, "control"),
    file.path(tab_dir, paste0(comparison_name, "_EVS_control_PCA_summary_raw.csv"))
  )

  if (isTRUE(export_optional_evs_variance_profiles)) {
    variance_panel <- safe_plot_build(
      plot_combined_pca_variance_panel(evs, comparison_name),
      paste0(comparison_name, ": EVS variance panel")
    )
    if (!is.null(variance_panel)) {
      save_grob(
        variance_panel,
        file.path(cmp_dir, "EVS_PCA_variance_profiles_combined.png"),
        width = 18, height = 13
      )
    }

    pc1_rank_panel <- safe_plot_build(
      arrangeGrob(
        plot_pc1_loading_rank(evs$fit_trt$loading_table, evs$fit_trt$cutoff, comparison_name, "treatment",
                              top_n_used = evs$fit_trt$top_n_used, cutoff_quantile = evs$fit_trt$cutoff_quantile,
                              preprocessing_label = evs$fit_trt$preprocessing_label),
        plot_pc1_loading_rank(evs$fit_untrt$loading_table, evs$fit_untrt$cutoff, comparison_name, "control",
                              top_n_used = evs$fit_untrt$top_n_used, cutoff_quantile = evs$fit_untrt$cutoff_quantile,
                              preprocessing_label = evs$fit_untrt$preprocessing_label),
        plot_pc1_loading_rank(evs$fit_trt_raw$loading_table, evs$fit_trt_raw$cutoff, comparison_name, "treatment",
                              top_n_used = evs$fit_trt_raw$top_n_used, cutoff_quantile = evs$fit_trt_raw$cutoff_quantile,
                              preprocessing_label = evs$fit_trt_raw$preprocessing_label),
        plot_pc1_loading_rank(evs$fit_untrt_raw$loading_table, evs$fit_untrt_raw$cutoff, comparison_name, "control",
                              top_n_used = evs$fit_untrt_raw$top_n_used, cutoff_quantile = evs$fit_untrt_raw$cutoff_quantile,
                              preprocessing_label = evs$fit_untrt_raw$preprocessing_label),
        ncol = 2,
        top = textGrob(
          paste0(comparison_name, " | EVS PC1 loading ranks"),
          gp = gpar(fontface = "bold", cex = 1.02)
        )
      ),
      paste0(comparison_name, ": EVS PC1 rank panel")
    )

    if (!is.null(pc1_rank_panel)) {
      save_grob(
        pc1_rank_panel,
        file.path(cmp_dir, "EVS_PC1_loading_rank_combined.png"),
        width = 18, height = 13
      )
    }
  }

  if (!is.null(evs$fit_trt$candidate_table) && nrow(evs$fit_trt$candidate_table)) {
    save_csv(evs$fit_trt$candidate_table, file.path(tab_dir, paste0(comparison_name, "_EVS_treatment_candidates_normalized.csv")))
  }
  if (!is.null(evs$fit_trt$matched_interval_table) && nrow(evs$fit_trt$matched_interval_table)) {
    save_csv(evs$fit_trt$matched_interval_table, file.path(tab_dir, paste0(comparison_name, "_EVS_treatment_matched_curvature_intervals.csv")))
  }
  if (!is.null(evs$fit_untrt$candidate_table) && nrow(evs$fit_untrt$candidate_table)) {
    save_csv(evs$fit_untrt$candidate_table, file.path(tab_dir, paste0(comparison_name, "_EVS_control_candidates_normalized.csv")))
  }
  if (!is.null(evs$fit_untrt$matched_interval_table) && nrow(evs$fit_untrt$matched_interval_table)) {
    save_csv(evs$fit_untrt$matched_interval_table, file.path(tab_dir, paste0(comparison_name, "_EVS_control_matched_curvature_intervals.csv")))
  }
  if (!is.null(evs$fit_trt_raw$candidate_table) && nrow(evs$fit_trt_raw$candidate_table)) {
    save_csv(evs$fit_trt_raw$candidate_table, file.path(tab_dir, paste0(comparison_name, "_EVS_treatment_candidates_raw.csv")))
  }
  if (!is.null(evs$fit_untrt_raw$candidate_table) && nrow(evs$fit_untrt_raw$candidate_table)) {
    save_csv(evs$fit_untrt_raw$candidate_table, file.path(tab_dir, paste0(comparison_name, "_EVS_control_candidates_raw.csv")))
  }

  treatment_peak_grid_normalized <- evs$primary_trt_fit$quantile_screen_table
  if (!is.null(treatment_peak_grid_normalized) && nrow(treatment_peak_grid_normalized)) {
    save_csv(treatment_peak_grid_normalized, file.path(tab_dir, paste0(comparison_name, "_EVS_treatment_peak_grid_normalized.csv")))
  }

  control_peak_grid_normalized <- evs$primary_untrt_fit$quantile_screen_table
  if (!is.null(control_peak_grid_normalized) && nrow(control_peak_grid_normalized)) {
    save_csv(control_peak_grid_normalized, file.path(tab_dir, paste0(comparison_name, "_EVS_control_peak_grid_normalized.csv")))
  }

  if (!is.null(evs$fit_trt$quantile_screen_table) && nrow(evs$fit_trt$quantile_screen_table)) {
    save_csv(evs$fit_trt$quantile_screen_table, file.path(tab_dir, paste0(comparison_name, "_EVS_treatment_quantile_screen_candidates.csv")))
  }
  if (!is.null(evs$fit_untrt$quantile_screen_table) && nrow(evs$fit_untrt$quantile_screen_table)) {
    save_csv(evs$fit_untrt$quantile_screen_table, file.path(tab_dir, paste0(comparison_name, "_EVS_control_quantile_screen_candidates.csv")))
  }

  if (!is.null(evs$independent_candidate_evals)) {
    for (nm in names(evs$independent_candidate_evals)) {
      cand_res <- evs$independent_candidate_evals[[nm]]$candidate_results
      if (!is.null(cand_res) && nrow(cand_res)) {
        save_csv(cand_res, file.path(tab_dir, paste0(comparison_name, "_EVS_", nm, "_candidate_grid.csv")))
      }
    }
  }
  primary_candidate_panel <- safe_plot_build(
    plot_primary_evs_candidate_panel(evs, comparison_name),
    paste0(comparison_name, ": primary EVS candidate panel")
  )
  if (!is.null(primary_candidate_panel)) {
    save_grob(
      primary_candidate_panel,
      file.path(cmp_dir, "EVS_primary_candidate_cutoffs.png"),
      width = 18, height = 8
    )
  }

  dataset_list <- list(
    raw_dataset          = evs$raw_dataset,
    leading_edge_dataset = evs$leading_edge_dataset,
    remainder_dataset    = evs$remainder_dataset
  )
  normalized_dataset_list <- list(
    raw_dataset = subset_normalized_counts_by_ids(
      normalized_counts = evs$normalized_counts,
      feature_ids = rownames(evs$normalized_counts),
      dataset_label = "raw_dataset",
      comparison_name = comparison_name,
      strict = TRUE
    ),
    leading_edge_dataset = subset_normalized_counts_by_ids(
      normalized_counts = evs$normalized_counts,
      feature_ids = rownames(evs$leading_edge_dataset),
      dataset_label = "leading_edge_dataset",
      comparison_name = comparison_name,
      strict = TRUE
    ),
    remainder_dataset = subset_normalized_counts_by_ids(
      normalized_counts = evs$normalized_counts,
      feature_ids = rownames(evs$remainder_dataset),
      dataset_label = "remainder_dataset",
      comparison_name = comparison_name,
      strict = TRUE
    )
  )

  dataset_fig_dirs <- list(
    raw_dataset          = fig_raw_dir,
    leading_edge_dataset = fig_le_dir,
    remainder_dataset    = fig_rem_dir
  )

  analysis_results <- list()

  for (nm in names(dataset_list)) {
    full_dataset_name <- paste(comparison_name, nm, sep = "_")

    fit <- tryCatch(
      run_core_analysis(
        count_mat    = dataset_list[[nm]],
        coldata      = coldata,
        dataset_name = full_dataset_name,
        annot_df     = annot_df
      ),
      error = function(e) {
        warning(sprintf("[%s][%s] run_core_analysis failed: %s", comparison_name, nm, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(fit)) next

    df         <- fit$results
    fig_subdir <- dataset_fig_dirs[[nm]]

    save_csv(df, file.path(tab_dir, paste0(full_dataset_name, "_results_full.csv")))
    save_csv(subset(df, standard_significant), file.path(tab_dir, paste0(full_dataset_name, "_standard_significant.csv")))
    save_csv(subset(df, HBFSS_significant),    file.path(tab_dir, paste0(full_dataset_name, "_HBFSS_significant.csv")))
    save_csv(subset(df, effect_class == "strong_effect"), file.path(tab_dir, paste0(full_dataset_name, "_strong_effect.csv")))
    save_csv(subset(df, effect_class == "weak_effect"),   file.path(tab_dir, paste0(full_dataset_name, "_weak_effect.csv")))

    summary_row <- data.frame(
      comparison_name        = comparison_name,
      dataset_name           = full_dataset_name,
      n_features             = nrow(df),
      hc_p_threshold         = fit$hc_p_threshold,
      hbfss_threshold        = fit$hbfss_threshold,
      n_standard_significant = sum(df$standard_significant, na.rm = TRUE),
      n_HBFSS_significant    = sum(df$HBFSS_significant, na.rm = TRUE),
      n_overlap_significant  = sum(df$deseq2_strong_call & df$HBFSS_significant, na.rm = TRUE),
      n_strong_effect        = sum(df$effect_class == "strong_effect", na.rm = TRUE),
      n_weak_effect          = sum(df$effect_class == "weak_effect", na.rm = TRUE),
      evs_fixed_top_n        = evs_fixed_top_n,
      evs_cutoff_mode_main   = evs_cutoff_mode_main,
      evs_primary_preprocessing = "normalized",
      shared_cutoff_rank     = evs$combined_cutoff_info$top_n_actual,
      shared_cutoff_method   = evs$combined_cutoff_info$method,
      shared_cutoff_quantile = evs$fit_trt$cutoff_quantile,
      shared_selected_reason = evs$combined_cutoff_info$selected_reason,
      stringsAsFactors       = FALSE
    )

    save_csv(summary_row, file.path(tab_dir, paste0(full_dataset_name, "_summary.csv")))

    if (isTRUE(export_optional_mean_histograms)) {
      mean_df    <- compute_mean_expression_table(dataset_list[[nm]], coldata)
      mean_panel <- safe_plot_build(
        plot_mean_histogram_panel(mean_df, fit$base_mean_vec, full_dataset_name),
        paste0(full_dataset_name, ": mean histogram panel")
      )
      if (!is.null(mean_panel)) {
        save_grob(
          mean_panel,
          file.path(fig_subdir, paste0(full_dataset_name, "_mean_histograms.png")),
          width = 14, height = 10
        )
      }
    }

    disp_resid <- tryCatch(
      dispersion_residual_section(df, full_dataset_name, fig_subdir),
      error = function(e) {
        warning(paste0(full_dataset_name, ": dispersion residual export failed: ", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(disp_resid)) {
      save_csv(disp_resid, file.path(tab_dir, paste0(full_dataset_name, "_dispersion_residuals.csv")))
    }

    analysis_results[[nm]] <- list(
      dds         = fit$dds,
      results     = df,
      summary     = summary_row,
      dataset_mat = dataset_list[[nm]],
      fig_subdir  = fig_subdir
    )
  }

  tryCatch(
    build_pc1_feature_export(
      analysis_results = analysis_results,
      annot_df         = annot_df,
      comparison_name  = comparison_name,
      evs              = evs,
      tab_dir          = tab_dir
    ),
    error = function(e) {
      warning(paste0(comparison_name, ": PC1 feature export failed: ", conditionMessage(e)))
      NULL
    }
  )

  tryCatch(
    {
      save_cross_dataset_comparison_panels(
        comparison_name,
        analysis_results,
        cmp_dir,
        dataset_list,
        coldata,
        normalized_dataset_list = normalized_dataset_list
      )
    },
    error = function(e) {
      warning(paste0(comparison_name, ": cross-dataset panel generation failed: ", conditionMessage(e)))
    }
  )

  if (run_twas_overlap) {
    twas_genes <- clean_gene_set(TWAS_data$gene_symbol)

    get_twas_overlap <- function(result_df, dataset_nm, cmp_name, out_dir,
                                 mode = c("union", "standard_only", "HBFSS_only")) {
      mode <- match.arg(mode)

      sig_df <- switch(
        mode,
        union         = subset(result_df, standard_significant | HBFSS_significant),
        standard_only = subset(result_df, standard_significant),
        HBFSS_only    = subset(result_df, HBFSS_significant)
      )

      sig_df$gene_symbol_clean <- tolower(trimws(sig_df$gene_symbol))
      overlap_df <- subset(sig_df, gene_symbol_clean %in% twas_genes)

      summary_df <- data.frame(
        comparison_name     = cmp_name,
        dataset_name        = dataset_nm,
        selection_mode      = mode,
        n_selected_features = nrow(sig_df),
        n_overlap_features  = nrow(overlap_df),
        n_overlap_genes     = length(unique(overlap_df$gene_symbol_clean)),
        stringsAsFactors    = FALSE
      )

      save_csv(overlap_df, file.path(out_dir, paste0(dataset_nm, "_TWAS_overlap_", mode, ".csv")))
      save_csv(summary_df, file.path(out_dir, paste0(dataset_nm, "_TWAS_overlap_summary_", mode, ".csv")))

      summary_df
    }

    twas_summaries <- dplyr::bind_rows(lapply(names(analysis_results), function(nm) {
      full_nm <- paste(comparison_name, nm, sep = "_")
      dplyr::bind_rows(
        get_twas_overlap(analysis_results[[nm]]$results, full_nm, comparison_name, tab_dir, "union"),
        get_twas_overlap(analysis_results[[nm]]$results, full_nm, comparison_name, tab_dir, "standard_only"),
        get_twas_overlap(analysis_results[[nm]]$results, full_nm, comparison_name, tab_dir, "HBFSS_only")
      )
    }))

    save_csv(twas_summaries, file.path(tab_dir, "TWAS_overlap_overall_summary.csv"))
  }

  sm_list <- lapply(analysis_results, `[[`, "summary")
  sm_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, sm_list)

  if (length(sm_list) == 0) return(data.frame())

  dplyr::bind_rows(sm_list)
}

comparison_inputs <- lapply(seq_len(nrow(comparison_table)), function(i) {
  prepare_comparison_data(
    comparison_name = comparison_table$comparison_name[i],
    group1_prefix   = comparison_table$group1_prefix[i],
    group2_prefix   = comparison_table$group2_prefix[i],
    WTTS_Seq        = WTTS_Seq,
    meta_all        = meta_all
  )
})

names(comparison_inputs) <- comparison_table$comparison_name

all_summaries_list <- list()
failed_comparisons <- list()

for (cmp in names(comparison_inputs)) {
  message("\n=====================================================")
  message("Running comparison: ", cmp)
  message("=====================================================")

  input_obj <- comparison_inputs[[cmp]]

  out <- tryCatch(
    run_full_comparison_pipeline(
      comparison_name = input_obj$comparison_name,
      count_matrix    = input_obj$count_matrix,
      coldata         = input_obj$coldata,
      annot_df        = OrigID_Symbol
    ),
    error = function(e) {
      failed_comparisons[[cmp]] <<- data.frame(
        comparison_name  = cmp,
        error_message    = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      NULL
    }
  )

  if (!is.null(out) && nrow(out) > 0) {
    all_summaries_list[[cmp]] <- out
  }
}

all_summaries <- if (length(all_summaries_list) > 0) {
  dplyr::bind_rows(all_summaries_list)
} else {
  data.frame()
}

if (nrow(all_summaries) > 0) {
  save_csv(all_summaries, file.path(output_dir, "all_comparisons_overall_summary.csv"))
}

if (length(failed_comparisons) > 0) {
  failed_df <- dplyr::bind_rows(failed_comparisons)
  save_csv(failed_df, file.path(output_dir, "failed_comparisons.csv"))
}

cat("\n=====================================================\n")
cat("Pipeline complete.\n")
cat("Output directory:\n")
cat(normalizePath(output_dir), "\n")
cat("=====================================================\n\n")

if (nrow(all_summaries) > 0) {
  print(all_summaries)
} else {
  message("No comparison summaries were written. Check failed_comparisons.csv for the exact error message.")
}

session_info_txt <- capture.output(sessionInfo())
writeLines(session_info_txt, file.path(output_dir, "sessionInfo.txt"))
saveRDS(sessionInfo(), file.path(output_dir, "sessionInfo.rds"))
