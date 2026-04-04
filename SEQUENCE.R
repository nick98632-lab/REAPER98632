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

alpha_level <- 0.10
max_usable_hc_p_threshold <- 0.95
lfc_boundary <- 1.0
lfc_shrink_type <- "normal"   # "normal" for speed/stability; "apeglm" optional
lfc_shrink_apeglm_method <- "nbinomC"

evs_cutoff_mode_main <- "combined_global_fourier_last_crossing_before_divergence"
evs_fixed_top_n <- 5000L
evs_fixed_rank_override <- NA_integer_

fourier_percentile_step <- 0.01
fourier_window_fraction <- 0.12
fourier_harmonics <- 2L
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


build_global_fourier_wave_map <- function(loading_tbl,
                                          mean_col = "baseMean",
                                          dispersion_col = "dispGeneEst") {
  metric_df <- build_ranked_fourier_metric_table(
    loading_tbl = loading_tbl,
    mean_col = mean_col,
    dispersion_col = dispersion_col
  )
  if (is.null(metric_df) || !nrow(metric_df)) return(NULL)

  iod_fit <- fit_local_fourier(metric_df$rank, metric_df$log_iod_nb)
  cv2_fit <- fit_local_fourier(metric_df$rank, metric_df$log_cv2_nb)
  if (is.null(iod_fit) || is.null(cv2_fit)) return(NULL)

  n_total <- nrow(metric_df)

  wave_map <- data.frame(
    feature_id = metric_df$feature_id,
    rank = metric_df$rank,
    percentile = metric_df$rank / n_total,
    iod_fitted = as.numeric(iod_fit$fitted_y),
    cv2_fitted = as.numeric(cv2_fit$fitted_y),
    iod_residual = as.numeric(iod_fit$residual_y),
    cv2_residual = as.numeric(cv2_fit$residual_y),
    stringsAsFactors = FALSE
  )

  wave_map <- wave_map[order(wave_map$rank), , drop = FALSE]
  rownames(wave_map) <- NULL

  list(
    metric_df = metric_df,
    wave_map = wave_map,
    iod_fit = iod_fit,
    cv2_fit = cv2_fit
  )
}

combine_treatment_control_fourier_maps <- function(trt_wave_obj, ctrl_wave_obj) {
  if (is.null(trt_wave_obj) || is.null(ctrl_wave_obj)) return(NULL)
  trt_map <- trt_wave_obj$wave_map
  ctrl_map <- ctrl_wave_obj$wave_map
  if (is.null(trt_map) || is.null(ctrl_map) || !nrow(trt_map) || !nrow(ctrl_map)) return(NULL)

  trt_map2 <- trt_map[, c("feature_id", "rank", "percentile", "iod_fitted", "cv2_fitted"), drop = FALSE]
  ctrl_map2 <- ctrl_map[, c("feature_id", "rank", "percentile", "iod_fitted", "cv2_fitted"), drop = FALSE]

  names(trt_map2) <- c("feature_id", "rank_trt", "percentile_trt", "iod_fitted_trt", "cv2_fitted_trt")
  names(ctrl_map2) <- c("feature_id", "rank_ctrl", "percentile_ctrl", "iod_fitted_ctrl", "cv2_fitted_ctrl")

  out <- dplyr::inner_join(trt_map2, ctrl_map2, by = "feature_id")
  if (!nrow(out)) return(NULL)

  out$combined_center_rank <- round((out$rank_trt + out$rank_ctrl) / 2)
  out$percentile <- rowMeans(cbind(out$percentile_trt, out$percentile_ctrl), na.rm = TRUE)
  out$combined_iod_amplitude <- out$iod_fitted_trt + out$iod_fitted_ctrl
  out$combined_cv2_amplitude <- out$cv2_fitted_trt + out$cv2_fitted_ctrl

  out <- out[order(out$percentile, out$combined_center_rank, out$feature_id), , drop = FALSE]
  rownames(out) <- NULL
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

  selected <- crossing_tbl[order(crossing_tbl$crossing_percentile, decreasing = TRUE), , drop = FALSE][1, , drop = FALSE]

  list(
    diff_df = diff_df,
    crossing_table = crossing_tbl,
    selected_crossing = selected,
    selected_reason = "last_crossing_before_divergence"
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

  trt_wave_obj <- build_global_fourier_wave_map(fit_trt_loading_tbl)
  ctrl_wave_obj <- build_global_fourier_wave_map(fit_ctrl_loading_tbl)

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

  leading_edge_union_flag <- (
    (!is.na(shared_combined_tbl$rank_trt)  & shared_combined_tbl$rank_trt  <= selected_rank) |
    (!is.na(shared_combined_tbl$rank_ctrl) & shared_combined_tbl$rank_ctrl <= selected_rank)
  )
  leading_edge_union_n <- sum(leading_edge_union_flag, na.rm = TRUE)
  remainder_union_n <- n_total - leading_edge_union_n

  list(
    top_n_actual = selected_rank,
    cutoff_value = cutoff_value,
    cutoff_quantile = cutoff_quantile,
    n_total = n_total,
    rank_cutoff = selected_rank,
    leading_edge_union_n = leading_edge_union_n,
    remainder_union_n = remainder_union_n,
    method = "global_fourier_last_crossing_before_divergence",
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



plot_fourier_wave_map_single <- function(wave_obj, comparison_name, group_label) {
  if (is.null(wave_obj) || is.null(wave_obj$wave_map) || !nrow(wave_obj$wave_map)) return(NULL)
  df <- wave_obj$wave_map

  ggplot(df, aes(percentile)) +
    geom_line(aes(y = iod_fitted, color = "IOD"), linewidth = 0.9) +
    geom_line(aes(y = cv2_fitted, color = "CV²"), linewidth = 0.9) +
    scale_color_manual(values = c("IOD" = plot_palette$treatment, "CV²" = plot_palette$control)) +
    labs(
      title = paste0(pretty_group_label(group_label), " | full-axis fitted regime lines"),
      subtitle = compact_caption(
        "This panel shows the within-group full-axis Fourier fitted IOD and CV² regime lines across the entire ranked gene axis.",
        width = 88
      ),
      x = "Rank percentile",
      y = "Fourier fitted value",
      color = NULL
    ) +
    manuscript_theme()
}

plot_combined_fourier_wave_map <- function(combined_cutoff_info, comparison_name) {
  df <- combined_cutoff_info$combined_wave_map
  if (is.null(df) || !nrow(df)) return(NULL)

  sc <- combined_cutoff_info$selected_crossing
  crossing_pct <- if (!is.null(sc) && nrow(sc)) sc$crossing_percentile[1] else NA_real_
  crossing_rank <- if (!is.null(sc) && nrow(sc)) sc$crossing_rank[1] else NA_integer_
  leading_edge_union_n <- if (!is.null(combined_cutoff_info$leading_edge_union_n)) combined_cutoff_info$leading_edge_union_n else NA_integer_
  remainder_union_n <- if (!is.null(combined_cutoff_info$remainder_union_n)) combined_cutoff_info$remainder_union_n else NA_integer_

  ymax <- max(c(df$combined_iod_amplitude, df$combined_cv2_amplitude), na.rm = TRUE)

  p <- ggplot(df, aes(percentile)) +
    geom_line(aes(y = combined_iod_amplitude, color = "Composite IOD"), linewidth = crossing_plot_line_width) +
    geom_line(aes(y = combined_cv2_amplitude, color = "Composite CV²"), linewidth = crossing_plot_line_width) +
    scale_color_manual(values = c("Composite IOD" = plot_palette$treatment, "Composite CV²" = plot_palette$control)) +
    labs(
      title = paste0(comparison_name, " | treatment-control composite overlap"),
      subtitle = compact_caption(
        "The dashed line is the selected shared comparison cutoff based on the last crossing before divergence of the full-axis composite IOD and composite CV² fitted lines. The rank cutoff is projected back into the treatment and control ranked loading tables, and the leading edge is defined as the union of sites that fall within that cutoff in either group. Labels report the crossing rank and the post-union leading-edge and remainder counts.",
        width = 100
      ),
      x = "Rank percentile",
      y = "Composite Fourier fitted value",
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
          "Last crossing before divergence\n",
          "percentile = ", signif(crossing_pct, 4), "\n",
          "rank cutoff = ", crossing_rank, "\n",
          "leading-edge union n = ", leading_edge_union_n, "\n",
          "remainder n = ", remainder_union_n
        ),
        fill = "white",
        colour = plot_palette$threshold,
        size = 3.0,
        label.size = 0.15,
        vjust = -0.55
      )
  }

  p
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
      label = paste0("rank cutoff = ", crossing_rank),
      fill = "white",
      colour = plot_palette$threshold,
      size = 2.8,
      label.size = 0.15,
      vjust = -0.8
    ) +
    labs(
      title = paste0(comparison_name, " | IOD minus CV² crossing curve"),
      subtitle = compact_caption("Positive values indicate full-axis IOD dominance and negative values indicate full-axis CV² dominance. The last zero-crossing before divergence is the EVS cutoff.", width = 90),
      x = "Percentile center",
      y = "IOD − CV²"
    ) +
    manuscript_theme()
}

plot_crossing_summary_panel <- function(combined_cutoff_info, comparison_name) {
  p1 <- plot_fourier_wave_map_single(combined_cutoff_info$trt_wave_obj, comparison_name, "treatment")
  p2 <- plot_fourier_wave_map_single(combined_cutoff_info$ctrl_wave_obj, comparison_name, "control")
  p3 <- plot_combined_fourier_wave_map(combined_cutoff_info, comparison_name)
  p4 <- plot_regime_difference_curve(combined_cutoff_info, comparison_name)

  grobs <- Filter(Negate(is.null), list(p1, p2, p3, p4))
  if (!length(grobs)) return(NULL)

  arrangeGrob(
    grobs = grobs,
    ncol = 2,
    top = textGrob(
      paste0(comparison_name, " | treatment, control, overlap, and crossing panels"),
      gp = gpar(fontface = "bold", cex = 1.04)
    )
  )
}

plot_fourier_summary_panel <- function(combined_cutoff_info, comparison_name) {
  plot_crossing_summary_panel(combined_cutoff_info, comparison_name)
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
    rank_cutoff = combined_cutoff_info$rank_cutoff,
    n_total_union = combined_cutoff_info$n_total,
    leading_edge_union_n = combined_cutoff_info$leading_edge_union_n,
    remainder_union_n = combined_cutoff_info$remainder_union_n,
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
  if (is.finite(evs_fixed_rank_override) && !is.na(evs_fixed_rank_override)) {
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

  final_shared_rank <- as.integer(combined_cutoff_info$top_n_actual)
  final_shared_reason <- combined_cutoff_info$selected_reason

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


  shared_combined_tbl <- combined_cutoff_info$shared_combined_tbl
  if (is.null(shared_combined_tbl) || !nrow(shared_combined_tbl)) {
    stop("Shared combined loading table is missing after shared cutoff resolution.", call. = FALSE)
  }

  leading_edge_union_flag <- (
    (!is.na(shared_combined_tbl$rank_trt)  & shared_combined_tbl$rank_trt  <= final_shared_rank) |
    (!is.na(shared_combined_tbl$rank_ctrl) & shared_combined_tbl$rank_ctrl <= final_shared_rank)
  )
  leading_edge_ids <- as.character(
    shared_combined_tbl$feature_id[leading_edge_union_flag]
  )
  analyzed_feature_ids <- as.character(shared_combined_tbl$feature_id)
  remainder_ids <- as.character(
    shared_combined_tbl$feature_id[!leading_edge_union_flag]
  )

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

  list(
    fit_trt               = fit_trt,
    fit_untrt             = fit_untrt,
    fit_trt_raw           = fit_trt_raw,
    fit_untrt_raw         = fit_untrt_raw,
    combined_cutoff_info  = combined_cutoff_info,
    shared_combined_tbl   = shared_combined_tbl,
    evs_cutoff_summary    = evs_cutoff_summary,
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
    save_csv(evs$combined_cutoff_info$trt_wave_obj$wave_map, file.path(tab_dir, paste0(comparison_name, "_treatment_full_axis_fourier_map.csv")))
  }
  if (!is.null(evs$combined_cutoff_info$ctrl_wave_obj$wave_map) && nrow(evs$combined_cutoff_info$ctrl_wave_obj$wave_map)) {
    save_csv(evs$combined_cutoff_info$ctrl_wave_obj$wave_map, file.path(tab_dir, paste0(comparison_name, "_control_full_axis_fourier_map.csv")))
  }
  if (!is.null(evs$combined_cutoff_info$combined_wave_map) && nrow(evs$combined_cutoff_info$combined_wave_map)) {
    save_csv(evs$combined_cutoff_info$combined_wave_map, file.path(tab_dir, paste0(comparison_name, "_last_crossing_map.csv")))
  }
  if (!is.null(evs$combined_cutoff_info$crossing_table) && nrow(evs$combined_cutoff_info$crossing_table)) {
    save_csv(evs$combined_cutoff_info$crossing_table, file.path(tab_dir, paste0(comparison_name, "_all_regime_crossings.csv")))
  }
  if (!is.null(evs$combined_cutoff_info$candidate_table) && nrow(evs$combined_cutoff_info$candidate_table)) {
    save_csv(evs$combined_cutoff_info$candidate_table, file.path(tab_dir, paste0(comparison_name, "_selected_regime_crossing.csv")))
  }
  crossing_summary_df <- build_crossing_summary_table(evs$combined_cutoff_info, comparison_name)
  save_csv(crossing_summary_df, file.path(tab_dir, paste0(comparison_name, "_regime_shift_crossing_summary.csv")))

  trt_fourier_plot <- safe_plot_build(
    plot_fourier_wave_map_single(evs$combined_cutoff_info$trt_wave_obj, comparison_name, "treatment"),
    paste0(comparison_name, ": treatment full-axis fourier wave map")
  )
  if (!is.null(trt_fourier_plot)) {
    save_grob(trt_fourier_plot, file.path(cmp_dir, paste0(comparison_name, "_treatment_full_axis_fourier_wave_map.png")), width = 12, height = 10)
  }
  ctrl_fourier_plot <- safe_plot_build(
    plot_fourier_wave_map_single(evs$combined_cutoff_info$ctrl_wave_obj, comparison_name, "control"),
    paste0(comparison_name, ": control full-axis fourier wave map")
  )
  if (!is.null(ctrl_fourier_plot)) {
    save_grob(ctrl_fourier_plot, file.path(cmp_dir, paste0(comparison_name, "_control_full_axis_fourier_wave_map.png")), width = 12, height = 10)
  }
  combined_fourier_plot <- safe_plot_build(
    plot_combined_fourier_wave_map(evs$combined_cutoff_info, comparison_name),
    paste0(comparison_name, ": combined full-axis fourier wave map")
  )
  if (!is.null(combined_fourier_plot)) {
    save_grob(combined_fourier_plot, file.path(cmp_dir, paste0(comparison_name, "_combined_full_axis_fourier_wave_map.png")), width = 12, height = 6)
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

  evs_hist_panel <- safe_plot_build(
    arrangeGrob(
      plot_eigenvector_histograms(
        evs$fit_trt$loading_table,
        evs$fit_trt$cutoff,
        comparison_name,
        "treatment",
        evs$fit_trt$preprocessing_label
      ),
      plot_eigenvector_histograms(
        evs$fit_untrt$loading_table,
        evs$fit_untrt$cutoff,
        comparison_name,
        "control",
        evs$fit_untrt$preprocessing_label
      ),
      ncol = 1,
      top = textGrob(
        paste0(comparison_name, " | treatment and control EVS histograms"),
        gp = gpar(fontface = "bold", cex = 1.02)
      )
    ),
    paste0(comparison_name, ": EVS histogram panel")
  )
  if (!is.null(evs_hist_panel)) {
    save_grob(
      evs_hist_panel,
      file.path(cmp_dir, paste0(comparison_name, "_EVS_histogram_panel.png")),
      width = 18, height = 12
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
