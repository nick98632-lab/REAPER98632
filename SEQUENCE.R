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
# All inputs are read from the repository data/ directory. All outputs
# produced in a given run are written to exports/sequence_run_<timestamp>/ so
# the entire analysis can be versioned from within the same repository.
# -----------------------------------------------------------------------------

repo_dir <- getwd()
input_dir <- file.path(repo_dir, "data")
output_root <- file.path(repo_dir, "exports")
analysis_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(output_root, paste0("sequence_run_", analysis_stamp))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")

# Set to TRUE and supply twas_file to enable downstream TWAS gene overlap.
run_twas_overlap <- FALSE
twas_file        <- file.path(input_dir, "3aTWAS_genes_of_11_brain_disorders.csv")

# Significance threshold for all DESeq2 calls and effect-class classification.
alpha_level <- 0.20
max_usable_hc_p_threshold <- 0.99

# HC configuration used by HBFSS:
# The manuscript pipeline derives empirical-null p-values from DESeq2 Wald
# statistics using fdrtool, then applies a Higher Criticism cutoff to those
# empirical-null p-values. That HC p-threshold is then transformed into the
# HBFSS decision boundary as abs(log10(hc_p_threshold_dataset)) * lfc_boundary.
# In this analysis a valid resolved HC p-threshold is always written, and high
# but finite thresholds remain valid up to max_usable_hc_p_threshold.

# LFC boundary (log2 scale) used for:
# (1) strong/weak effect hypothesis tests via lfcThreshold,
# (2) the HBFSS threshold anchor: abs(log10(hc_p_threshold_dataset)) * lfc_boundary.
lfc_boundary <- 1.0

# Default number of top-ranked PAS features (by absolute PC1 loading) to retain
# when a fixed EVS cutoff is requested.

# Primary EVS pathway used for cutoff selection in the manuscript.
# The current manuscript method builds separate treatment and control loading
# tables, estimates local Fourier structure for ranked IOD and CV2 within each,
# combines those local Fourier summaries into one comparison-level transition
# score, and selects one shared cutoff rank per comparison.
#
# EVS cutoff mode.
# "matched_curvature_wave" now routes to the combined local Fourier shared-
# cutoff method retained under the legacy option name for backward
# compatibility with earlier scripts and exports.
# "variance_second_derivative" retains the earlier NB2-only curvature fallback.
# "fixed_top_n" uses a user-specified manual rank.
evs_cutoff_mode_main <- "matched_curvature_wave"

# User-adjustable fixed EVS cutoff. This is used whenever
# evs_cutoff_mode_main == "fixed_top_n" and is also retained in summaries so a
# reader can compare the manual rank choice with the empirically derived rank.
evs_fixed_top_n <- 5000
evs_fixed_rank_override <- NA_integer_
evs_use_manual_rank_first <- FALSE

# EVS curve metrics used for adaptive cutoff detection.
# The current manuscript pathway uses combined local Fourier summaries of ranked
# IOD and CV2 as the primary shared-cutoff method. The older curvature-based NB2
# functions are retained below only as explicit legacy fallback code paths and
# for method comparison, not as the primary manuscript selector.

# Legacy curvature fallback settings.
# These parameters are used only if the script is intentionally routed through
# the older NB2 second-derivative fallback. They are retained so prior methods
# can still be reproduced, but they do not define the primary local Fourier
# manuscript workflow.
evs_curvature_spar          <- 0.55
evs_curvature_min_rank      <- 50L
evs_curvature_max_rank_frac <- 0.50
evs_curvature_spar_iod      <- 0.55
evs_curvature_spar_cv2      <- 0.55
evs_curvature_match_window  <- 150L
evs_candidate_interval_halfwidth <- 75L
evs_candidate_merge_distance <- 100L
wave_local_window <- 51L
wave_w_iod_slope <- 1.00
wave_w_iod_amp   <- 0.50
wave_w_cv2_amp   <- 1.00
wave_w_cv2_slope <- 0.50

# Legacy quantile-screened second-derivative EVS cutoff.
# These settings govern the older NB2 curvature fallback only. They are kept for
# direct comparison against the manuscript's combined local Fourier shared-cutoff
# workflow and should not be interpreted as the primary manuscript method.
evs_candidate_min_peak_distance   <- 40L
evs_peak_quantile_grid            <- seq(0.90, 0.75, by = -0.05)
evs_peak_quantile_screen_min      <- 0.75
evs_peak_quantile_screen_max      <- 0.90
# Discrete candidate-scoring settings.
# The primary workflow first identifies one shared comparison-level cutoff from
# the combined local Fourier map. The settings below are then used only for the
# downstream constrained candidate-scoring step applied to the retained shared-
# cutoff candidate set.
evs_candidate_max_leading_frac     <- 0.50
evs_candidate_runtime_seconds_low  <- 12
evs_candidate_runtime_seconds_high <- 25

# Local Fourier regime-shift settings.
# The Fourier grid is evaluated over the full ranked dataset using 100
# percentile-centered overlapping local fits. Treatment and control are modeled
# separately within each metric, then combined at the waveform level to form a
# composite IOD curve and a composite CV² curve. The primary EVS cutoff is the
# first stable crossing between these two lines. The percentile-indexed score
# is retained as a descriptive summary, but the crossing itself is the regime
# shift and therefore the splitting cutoff.
fourier_percentile_step <- 0.01
fourier_window_fraction <- 0.12
fourier_harmonics <- 2L
fourier_top_candidate_n_for_deseq2 <- 1L

# Evaluate the full ranked domain with 100 overlapping percentile-centered
# local fits. The resulting waveform panels display the local Fourier-derived
# amplitudes across the entire ranked dataset rather than only the upper tail.
fourier_min_rank <- 1L
fourier_max_rank_frac <- 1.00

# The score panel is a percentile-indexed decision statistic derived from the
# composite waveform objects. It is not itself treated as the biological wave.
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

# Plot settings
figure_dpi            <- 320
base_theme_size       <- 10

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
# Keep manuscript-essential figures by default. Disable redundant diagnostics.
export_optional_mean_histograms            <- FALSE
export_optional_empirical_p_histograms     <- FALSE
export_optional_hbfss_distributions        <- FALSE
export_optional_evs_raw_histograms         <- FALSE
export_optional_evs_variance_profiles      <- TRUE

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
# max_usable_hc_p_threshold : Maximum reportable HC p-threshold used when the
#   empirical-null p-value distribution is broad. The manuscript workflow still
#   writes a numeric HC threshold up to this cap rather than treating high
#   thresholds as unavailable or non-applicable.
# HBFSS / pi_valueE       : abs(apeglm-shrunken log2FC * log10(empirical_p)).
#                           Because empirical_p lies in (0, 1], log10(empirical_p) <= 0.
#                           The absolute value is therefore essential: it converts the
#                           negative log10-scaled evidence term into a positive magnitude
#                           so larger HBFSS values reflect stronger combined effect-size
#                           and empirical-null evidence.
# standard_significant    : DESeq2-positive call: native DESeq2 padj < alpha,
#                           |rawLFC| >= lfc_boundary, and |lfc_shrunk| >=
#                           lfc_boundary.
# evs_cutoff_mode_main    : main EVS split rule. In the current manuscript
#                           workflow, "matched_curvature_wave" routes to the
#                           combined local Fourier shared-cutoff method for
#                           backward compatibility. "variance_second_derivative"
#                           retains the earlier NB2-only curvature fallback.
#                           "fixed_top_n" uses a user-specified manual rank.
# evs_fixed_top_n         : user-adjustable manual EVS rank cutoff when a fixed
#                           cutoff is requested.
# comparison-level cutoff methodology :
#                           treatment and control loading tables are built
#                           separately, local Fourier wave summaries of ranked
#                           IOD and CV2 are estimated across percentile windows,
#                           treatment and control local scores are combined into
#                           one comparison-level transition score, and one shared
#                           cutoff rank is selected for the comparison. That rank
#                           defines the leading-edge and remainder subsets used
#                           downstream.

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
      )
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
  
  if (length(stat_vec) < 5) {
    stop(sprintf("[%s] Fewer than 5 finite Wald statistics were available for fdrtool.", dataset_name))
  }
  
  fit <- tryCatch(
    fdrtool(
      stat_vec,
      statistic     = "normal",
      plot          = FALSE,
      verbose       = FALSE,
      cutoff.method = "fndr",
      pct0          = 0.75
    ),
    error = function(e1) {
      message(sprintf("[%s] Primary fdrtool call failed: %s", dataset_name, conditionMessage(e1)))
      tryCatch(
        fdrtool(
          as.vector(stat_vec),
          statistic     = "normal",
          plot          = FALSE,
          verbose       = FALSE,
          cutoff.method = "pct0",
          pct0          = 0.75
        ),
        error = function(e2) {
          stop(sprintf("[%s] fdrtool failed after retry: %s", dataset_name, conditionMessage(e2)))
        }
      )
    }
  )
  
  fit$pval <- clip_probabilities(fit$pval)
  fit$qval <- clip_probabilities(fit$qval)
  fit$lfdr <- as.numeric(fit$lfdr)
  fit
}

safe_hc_thresh <- function(empirical_p, dataset_name) {
  # Resolve one numeric Higher Criticism p-threshold from empirical-null
  # p-values derived from the DESeq2 Wald statistics. The primary estimator is
  # fdrtool::hc.thresh(). If that does not return a usable numeric value, a
  # deterministic ordered-p fallback based on the classical HC objective is
  # used. This function never returns NA for a non-empty finite p-value vector.
  sorted_empirical_p <- sort(
    clip_probabilities(empirical_p),
    na.last    = NA,
    decreasing = FALSE
  )

  if (!length(sorted_empirical_p)) {
    stop(sprintf("[%s] No finite empirical-null p-values were available for HC threshold resolution.", dataset_name))
  }

  resolve_hc_fallback <- function(p_sorted) {
    n <- length(p_sorted)
    i <- seq_len(n)
    denom <- sqrt(pmax(p_sorted * (1 - p_sorted), .Machine$double.eps))
    hc_stat <- sqrt(n) * ((i / n) - p_sorted) / denom
    hc_stat[!is.finite(hc_stat)] <- -Inf
    best_idx <- which.max(hc_stat)
    if (!length(best_idx) || !is.finite(hc_stat[best_idx])) {
      return(as.numeric(p_sorted[min(n, max(1L, floor(0.5 * n)))]))
    }
    as.numeric(p_sorted[best_idx])
  }

  primary_out <- suppressWarnings(
    tryCatch(
      fdrtool::hc.thresh(as.vector(sorted_empirical_p)),
      error = function(e) {
        message(sprintf("[%s] hc.thresh failed; switching to deterministic HC fallback: %s", dataset_name, conditionMessage(e)))
        NA_real_
      }
    )
  )

  out <- as.numeric(primary_out[1])
  if (!is.finite(out) || out <= 0 || out > 1) {
    out <- resolve_hc_fallback(sorted_empirical_p)
  }

  out <- as.numeric(clip_probabilities(out))
  out <- min(out, max_usable_hc_p_threshold)

  if (!is.finite(out) || out <= 0 || out > 1) {
    stop(sprintf("[%s] HC threshold resolution failed after primary and fallback paths.", dataset_name))
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
  if (n_total == 0) stop("resolve_top_n_cutoff() received an empty vector.")

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

  # DESeq2 / NB2 variance proxy: Var(Y) = mu + alpha * mu^2.
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

# Detect local peaks in curvature along the smoothed ranked NB2 variance
# curve. If no interior local maximum survives the peak logic, the function
# falls back to the strongest interior curvature value so the adaptive pathway
# still returns at least one interpretable candidate rather than NULL.
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
      'peak_quantile_cutoff',
      'peak_strength_threshold',
      'peak_order_within_quantile',
      'rank_index',
      'percent_rank',
      'curvature_strength',
      'prominence',
      'smoothed_log_nb_variance',
      'first_derivative',
      'second_derivative'
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
# These functions translate DESeq2 mean/dispersion summaries onto the EVS rank
# axis so the cutoff can be identified from ranked IOD and CV2 behavior rather
# than from loading magnitude alone.
#
# build_ranked_regime_metric_table():
#   Creates a ranked table containing mean, dispersion, IOD, and CV2 for every
#   retained feature in a dataset.
#
# smooth_ranked_regime_curves():
#   Smooths the ranked IOD and CV2 curves and returns their first- and second-
#   derivative summaries, which are used to identify local transition structure.
#
# find_local_metric_curvature_peaks():
#   Finds locally prominent curvature peaks for either IOD or CV2 along the
#   ranked axis.
#
# match_curvature_peak_intervals():
#   Creates plausible cutoff intervals by matching IOD and CV2 curvature peaks
#   that occur within a local tolerance window.
#
# score_wave_candidate_rank():
#   Scores an exact rank inside a matched interval using local slope and
#   amplitude contrasts from the smoothed IOD and CV2 curves.
#
# resolve_matched_curvature_wave_cutoff():
#   Executes the full manuscript cutoff workflow for one dataset and returns the
#   selected rank together with the intermediate diagnostic tables needed for
#   export and plotting.
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
  if (length(rank_vec) < (2 * n_harmonics + 5L)) return(NULL)
  if (all(!is.finite(y_vec)) || stats::sd(y_vec, na.rm = TRUE) == 0) return(NULL)

  dd <- build_local_fourier_design(rank_vec, n_harmonics = n_harmonics)
  dd$y <- as.numeric(y_vec)

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

# Select the final stable crossing in the leading-edge search region. The
# manuscript method uses the last stable IOD-versus-CV2 crossing rather than
# the first crossing so the cutoff remains at the terminal boundary of the
# leading-edge regime.
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
    selected <- stable_tbl[order(stable_tbl$crossing_percentile, decreasing = TRUE), , drop = FALSE][1, , drop = FALSE]
    reason <- "last_stable_crossing"
  } else {
    selected <- crossing_tbl[order(crossing_tbl$crossing_percentile, decreasing = TRUE), , drop = FALSE][1, , drop = FALSE]
    reason <- "last_crossing_fallback"
  }

  list(
    diff_df = diff_df,
    crossing_table = crossing_tbl,
    selected_crossing = selected,
    selected_reason = reason
  )
}

resolve_combined_fourier_cutoff <- function(fit_trt_loading_tbl,
                                            fit_ctrl_loading_tbl,
                                            fixed_top_n = evs_fixed_top_n) {
  n_total <- nrow(fit_trt_loading_tbl)
  fallback <- resolve_top_n_cutoff(fit_trt_loading_tbl$pc1_loading_abs, top_n = fixed_top_n)

  trt_wave_obj <- build_local_fourier_wave_map(fit_trt_loading_tbl)
  ctrl_wave_obj <- build_local_fourier_wave_map(fit_ctrl_loading_tbl)
  if (is.null(trt_wave_obj) || is.null(ctrl_wave_obj)) {
    fallback$method <- "fixed_top_n_fallback"
    fallback$trt_wave_obj <- trt_wave_obj
    fallback$ctrl_wave_obj <- ctrl_wave_obj
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
    fallback$combined_wave_map <- crossing_info$diff_df
    fallback$combined_interval_table <- data.frame()
    fallback$crossing_table <- crossing_info$crossing_table
    fallback$selected_crossing <- NULL
    fallback$selected_reason <- crossing_info$selected_reason
    return(fallback)
  }

  selected_rank <- as.integer(selected_crossing$crossing_rank[1])
  selected_rank <- min(max(1L, selected_rank), n_total)
  cutoff_value <- fit_trt_loading_tbl$pc1_loading_abs[selected_rank]
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
    method = "last_stable_crossing",
    trt_wave_obj = trt_wave_obj,
    ctrl_wave_obj = ctrl_wave_obj,
    combined_wave_map = crossing_info$diff_df,
    combined_interval_table = data.frame(),
    crossing_table = crossing_info$crossing_table,
    selected_crossing = selected_crossing,
    candidate_table = candidate_table,
    selected_reason = crossing_info$selected_reason
  )
}

plot_fourier_wave_map_single <- function(wave_obj, comparison_name, group_label, selected_percentile = NA_real_) {
  if (is.null(wave_obj) || is.null(wave_obj$wave_map) || !nrow(wave_obj$wave_map)) return(NULL)
  df <- wave_obj$wave_map

  p <- ggplot(df, aes(percentile)) +
    geom_line(aes(y = iod_amplitude, color = "IOD"), linewidth = 0.9) +
    geom_line(aes(y = cv2_amplitude, color = "CV²"), linewidth = 0.9)
  if (is.finite(selected_percentile)) {
    p <- p +
      geom_vline(xintercept = selected_percentile, linetype = "dashed", linewidth = crossing_plot_vline_width, colour = plot_palette$threshold) +
      annotate("label", x = selected_percentile, y = max(c(df$iod_amplitude, df$cv2_amplitude), na.rm = TRUE), label = paste0("Selected cutoff
p = ", signif(selected_percentile, 4)), fill = "white", colour = plot_palette$threshold, size = 2.8, label.size = 0.15, vjust = -0.5)
  }
  p +
    scale_color_manual(values = c("IOD" = plot_palette$treatment, "CV²" = plot_palette$control)) +
    labs(
      title = paste0(pretty_group_label(group_label), " | local regime lines"),
      subtitle = compact_caption("Two lines are shown: local IOD and local CV². The dashed vertical line marks the selected leading-edge crossing used for the final shared EVS cutoff.", width = 88),
      x = "Percentile center",
      y = "Local amplitude",
      color = NULL
    ) +
    manuscript_theme()
}

plot_combined_fourier_wave_map <- function(combined_cutoff_info, comparison_name) {
  df <- combined_cutoff_info$combined_wave_map
  if (is.null(df) || !nrow(df)) return(NULL)

  sc <- combined_cutoff_info$selected_crossing
  crossing_pct <- if (!is.null(sc) && nrow(sc)) sc$crossing_percentile[1] else NA_real_

  ggplot(df, aes(percentile)) +
    geom_line(aes(y = combined_iod_amplitude, color = "Composite IOD"), linewidth = crossing_plot_line_width) +
    geom_line(aes(y = combined_cv2_amplitude, color = "Composite CV²"), linewidth = crossing_plot_line_width) +
    geom_vline(xintercept = crossing_pct, linetype = "dashed", linewidth = crossing_plot_vline_width, colour = plot_palette$threshold) +
    annotate(
      "label",
      x = crossing_pct,
      y = max(c(df$combined_iod_amplitude, df$combined_cv2_amplitude), na.rm = TRUE),
      label = paste0("Regime-shift crossing
p = ", signif(crossing_pct, 4)),
      fill = "white",
      colour = plot_palette$threshold,
      size = 3.0,
      label.size = 0.15,
      vjust = -0.5
    ) +
    scale_color_manual(values = c("Composite IOD" = plot_palette$treatment, "Composite CV²" = plot_palette$control)) +
    labs(
      title = paste0(comparison_name, " | two-line regime crossing"),
      subtitle = compact_caption("The dashed vertical line marks the first selected crossing between the composite local IOD and composite local CV² lines within the allowed leading-edge search range. This crossing defines the EVS cutoff.", width = 90),
      x = "Percentile center",
      y = "Composite local amplitude",
      color = NULL
    ) +
    manuscript_theme()
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
      subtitle = compact_caption("The score profile is descriptive only. The EVS cutoff is determined by the selected leading-edge crossing of the two regime lines, not by the global score maximum.", width = 90),
      x = "Percentile center",
      y = "Score"
    ) +
    manuscript_theme()
}

plot_regime_difference_curve <- function(combined_cutoff_info, comparison_name) {
  df <- combined_cutoff_info$combined_wave_map
  if (is.null(df) || !nrow(df) || !"regime_difference" %in% names(df)) return(NULL)

  sc <- combined_cutoff_info$selected_crossing
  crossing_pct <- if (!is.null(sc) && nrow(sc)) sc$crossing_percentile[1] else NA_real_

  ggplot(df, aes(percentile, regime_difference)) +
    geom_hline(yintercept = 0, linewidth = crossing_plot_hline_width, colour = "grey50") +
    geom_line(linewidth = crossing_plot_line_width, colour = plot_palette$threshold) +
    geom_vline(xintercept = crossing_pct, linetype = "dashed", linewidth = crossing_plot_vline_width, colour = plot_palette$threshold) +
    geom_point(data = data.frame(percentile = crossing_pct, regime_difference = 0), aes(x = percentile, y = regime_difference), inherit.aes = FALSE, size = crossing_plot_point_size, colour = plot_palette$threshold) +
    labs(
      title = paste0(comparison_name, " | regime-difference curve"),
      subtitle = compact_caption("Positive values indicate local IOD dominance and negative values indicate local CV² dominance. The selected leading-edge zero-crossing is the EVS cutoff.", width = 90),
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

# Adaptive EVS cutoff resolution.
# Input: a PC1 loading table already joined to DESeq2 mean/dispersion metrics.
# Logic: smooth the ranked NB2 variance curve, compute curvature, identify
# local curvature peaks, screen those peaks across a quantile grid, and choose
# the strongest surviving peak (breaking ties by prominence, then earlier rank)
# as the empirical EVS cutoff candidate.
# Output: the selected cutoff plus the full curve and candidate tables needed
# for audit plots, CSV exports, and downstream discrete candidate scoring.
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
    stop("candidate_tbl must contain rank_index.")
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
    stop(sprintf("[%s] count_mat must contain at least one feature and one sample.", dataset_name))
  }

  if (is.null(rownames(count_mat)) || anyNA(rownames(count_mat)) || any(rownames(count_mat) == "")) {
    stop(sprintf("[%s] count_mat must have non-empty rownames for all features.", dataset_name))
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
    stop(sprintf("[%s] Candidate-scoring analysis did not return the expected results/hc_p_threshold fields.", dataset_name))
  }

  fit
}

apply_loading_cutoff <- function(fit_obj, rank_index, selected_reason = "manual_override") {
  fit_obj <- as.list(fit_obj)
  loading_tbl <- as.data.frame(fit_obj$loading_table, stringsAsFactors = FALSE)
  if (is.null(loading_tbl) || !nrow(loading_tbl)) {
    stop("apply_loading_cutoff() requires a non-empty loading_table.")
  }
  required_cols <- c("feature_id", "pc1_loading_abs")
  assert_required_columns(loading_tbl, required_cols, object_name = "fit_obj$loading_table")

  if (!("rank" %in% colnames(loading_tbl)) ||
      anyNA(loading_tbl$rank) ||
      !identical(as.integer(loading_tbl$rank), seq_len(nrow(loading_tbl)))) {
    loading_tbl <- loading_tbl[order(loading_tbl$pc1_loading_abs, decreasing = TRUE, na.last = NA), , drop = FALSE]
    loading_tbl$rank <- seq_len(nrow(loading_tbl))
  } else {
    loading_tbl <- loading_tbl[order(loading_tbl$rank), , drop = FALSE]
  }

  fit_obj$loading_table <- loading_tbl
  n_total <- nrow(loading_tbl)
  rank_index <- as.integer(rank_index)[1]
  if (!is.finite(rank_index) || is.na(rank_index)) {
    stop("apply_loading_cutoff() requires rank_index to be a finite integer.")
  }
  rank_index <- min(max(1L, rank_index), n_total)
  row_idx <- match(rank_index, loading_tbl$rank)
  if (is.na(row_idx) || length(row_idx) != 1L) {
    stop(sprintf(
      "apply_loading_cutoff() could not locate rank %d in the loading_table after sorting.",
      rank_index
    ))
  }
  cutoff_value <- as.numeric(loading_tbl$pc1_loading_abs[row_idx])
  if (!is.finite(cutoff_value) || is.na(cutoff_value)) {
    stop(sprintf(
      "apply_loading_cutoff() produced an invalid cutoff_value for rank %d.",
      rank_index
    ))
  }
  cutoff_quantile <- 1 - (rank_index / n_total)

  fit_obj$cutoff <- cutoff_value
  fit_obj$top_n_used <- rank_index
  fit_obj$cutoff_quantile <- cutoff_quantile
  fit_obj$curvature_strength <- suppressWarnings({
    tbl <- fit_obj$candidate_table
    if (!is.null(tbl) && nrow(tbl) && any(tbl$rank_index == rank_index)) tbl$curvature_strength[match(rank_index, tbl$rank_index)] else NA_real_
  })
  fit_obj$cutoff_method <- paste0(fit_obj$cutoff_method, "_selected")
  fit_obj$selected_reason <- selected_reason
  fit_obj$loading_table$split_class <- ifelse(
    fit_obj$loading_table$pc1_loading_abs >= cutoff_value,
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

# Candidate-cutoff dispersion gate:
# retain the existing significance objective exactly as written, but restrict the
# candidate pool to splits whose leading-edge and remainder subsets exhibit the
# expected variance ordering before ranking by significance. We compute a simple
# NB1-style moment summary directly from the subset counts and an NB2-style
# summary from DESeq2 baseMean/dispersion outputs. A cutoff is considered
# dispersion-feasible if either model supports both directional inequalities:
#   IOD_leading  > IOD_remainder
#   CV2_remainder > CV2_leading
# where IOD = variance / mean and CV2 = variance / mean^2.
# This preserves the current optimization target while preventing the selected
# cutoff from violating the intended leading-edge-versus-remainder dispersion
# pattern.
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
    if (isTRUE(strict)) stop(msg) else warning(msg)
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

# Discrete constraint methodology:
# for each dataset, exhaustively evaluate all candidate cutoff ranks, retain only
# those whose leading-edge subset contains <= 50% of all PAS features and whose
# subset dispersion ordering matches the intended leading-edge/remainder pattern,
# maximize the summed downstream evidence across leading-edge and remainder
# subsets, and then aggregate the two independently optimized normalized-dataset
# ranks independently by dataset. This is a constrained exhaustive search, not a Lagrange
# relaxation.
evaluate_single_evs_dataset_candidates <- function(count_matrix,
                                                   coldata,
                                                   comparison_name,
                                                   fit_obj,
                                                   dataset_label,
                                                   max_leading_frac = evs_candidate_max_leading_frac) {
  candidate_tbl <- cap_evs_candidate_table(fit_obj$candidate_table)

  if (is.null(candidate_tbl) || !nrow(candidate_tbl)) {
    return(list(best_candidate = NULL, candidate_results = NULL))
  }

  all_rows <- vector("list", nrow(candidate_tbl))
  feature_universe <- as.character(fit_obj$loading_table$feature_id)
  n_total <- length(feature_universe)

  for (i in seq_len(nrow(candidate_tbl))) {
    rank_i <- as.integer(candidate_tbl$rank_index[i])
    cut_i  <- as.numeric(candidate_tbl$cutoff_value[i])

    leading_edge_ids <- as.character(
      fit_obj$loading_table$feature_id[fit_obj$loading_table$pc1_loading_abs >= cut_i]
    )
    remainder_ids <- setdiff(feature_universe, leading_edge_ids)

    n_leading <- length(leading_edge_ids)
    n_remainder <- length(remainder_ids)
    frac_leading <- n_leading / max(1L, n_total)

    valid_split <- (n_leading > 0L) &&
      (n_remainder > 0L) &&
      is.finite(frac_leading) &&
      (frac_leading <= max_leading_frac)

    lead_score <- rem_score <- combined_score <- objective_score <- NA_real_
    n_standard_leading <- NA_integer_
    n_standard_remainder <- NA_integer_
    nb1_mean_leading <- nb1_variance_leading <- nb1_iod_leading <- nb1_cv2_leading <- NA_real_
    nb1_mean_remainder <- nb1_variance_remainder <- nb1_iod_remainder <- nb1_cv2_remainder <- NA_real_
    nb2_mean_leading <- nb2_dispersion_leading <- nb2_iod_leading <- nb2_cv2_leading <- NA_real_
    nb2_mean_remainder <- nb2_dispersion_remainder <- nb2_iod_remainder <- nb2_cv2_remainder <- NA_real_
    nb1_pattern_pass <- nb2_pattern_pass <- dispersion_pattern_pass <- FALSE
    error_message <- NA_character_

    if (isTRUE(valid_split)) {
      lead_fit <- tryCatch(
        run_core_analysis_candidate_score(
          count_mat = count_matrix[leading_edge_ids, , drop = FALSE],
          coldata = coldata,
          dataset_name = paste0(comparison_name, "_", dataset_label, "_leading_edge_candidate")
        ),
        error = function(e) e
      )

      rem_fit <- tryCatch(
        run_core_analysis_candidate_score(
          count_mat = count_matrix[remainder_ids, , drop = FALSE],
          coldata = coldata,
          dataset_name = paste0(comparison_name, "_", dataset_label, "_remainder_candidate")
        ),
        error = function(e) e
      )

      if (inherits(lead_fit, "error") || inherits(rem_fit, "error")) {
        valid_split <- FALSE
        error_message <- paste(
          c(
            if (inherits(lead_fit, "error")) conditionMessage(lead_fit) else NULL,
            if (inherits(rem_fit, "error")) conditionMessage(rem_fit) else NULL
          ),
          collapse = " | "
        )
      } else {
        lead_eval <- score_evs_candidate_result(lead_fit$results)
        rem_eval  <- score_evs_candidate_result(rem_fit$results)

        lead_score <- lead_eval$score
        rem_score  <- rem_eval$score
        valid_component_scores <- c(lead_score, rem_score)
        valid_component_scores <- valid_component_scores[is.finite(valid_component_scores) & !is.na(valid_component_scores)]
        combined_score <- if (length(valid_component_scores)) sum(valid_component_scores) else NA_real_
        objective_score <- combined_score

        n_standard_leading   <- as.integer(lead_eval$metrics[["n_standard"]])
        n_standard_remainder <- as.integer(rem_eval$metrics[["n_standard"]])

        nb1_leading_summary <- compute_nb1_subset_moment_summary(count_matrix[leading_edge_ids, , drop = FALSE])
        nb1_remainder_summary <- compute_nb1_subset_moment_summary(count_matrix[remainder_ids, , drop = FALSE])
        nb2_leading_summary <- compute_nb2_subset_moment_summary(lead_fit$results)
        nb2_remainder_summary <- compute_nb2_subset_moment_summary(rem_fit$results)

        nb1_mean_leading <- unname(nb1_leading_summary[["mean_value"]])
        nb1_variance_leading <- unname(nb1_leading_summary[["variance_value"]])
        nb1_iod_leading <- unname(nb1_leading_summary[["iod"]])
        nb1_cv2_leading <- unname(nb1_leading_summary[["cv2"]])
        nb1_mean_remainder <- unname(nb1_remainder_summary[["mean_value"]])
        nb1_variance_remainder <- unname(nb1_remainder_summary[["variance_value"]])
        nb1_iod_remainder <- unname(nb1_remainder_summary[["iod"]])
        nb1_cv2_remainder <- unname(nb1_remainder_summary[["cv2"]])

        nb2_mean_leading <- unname(nb2_leading_summary[["mean_value"]])
        nb2_dispersion_leading <- unname(nb2_leading_summary[["dispersion_value"]])
        nb2_iod_leading <- unname(nb2_leading_summary[["iod"]])
        nb2_cv2_leading <- unname(nb2_leading_summary[["cv2"]])
        nb2_mean_remainder <- unname(nb2_remainder_summary[["mean_value"]])
        nb2_dispersion_remainder <- unname(nb2_remainder_summary[["dispersion_value"]])
        nb2_iod_remainder <- unname(nb2_remainder_summary[["iod"]])
        nb2_cv2_remainder <- unname(nb2_remainder_summary[["cv2"]])

        nb1_pattern_pass <- isTRUE(
          is.finite(nb1_iod_leading) && !is.na(nb1_iod_leading) &&
            is.finite(nb1_iod_remainder) && !is.na(nb1_iod_remainder) &&
            is.finite(nb1_cv2_leading) && !is.na(nb1_cv2_leading) &&
            is.finite(nb1_cv2_remainder) && !is.na(nb1_cv2_remainder) &&
            (nb1_iod_leading > nb1_iod_remainder) &&
            (nb1_cv2_remainder > nb1_cv2_leading)
        )
        nb2_pattern_pass <- isTRUE(
          is.finite(nb2_iod_leading) && !is.na(nb2_iod_leading) &&
            is.finite(nb2_iod_remainder) && !is.na(nb2_iod_remainder) &&
            is.finite(nb2_cv2_leading) && !is.na(nb2_cv2_leading) &&
            is.finite(nb2_cv2_remainder) && !is.na(nb2_cv2_remainder) &&
            (nb2_iod_leading > nb2_iod_remainder) &&
            (nb2_cv2_remainder > nb2_cv2_leading)
        )
        dispersion_pattern_pass <- nb1_pattern_pass || nb2_pattern_pass
        valid_split <- isTRUE(valid_split) && isTRUE(dispersion_pattern_pass)
      }
    }

    all_rows[[i]] <- data.frame(
      comparison_name = comparison_name,
      dataset_label = dataset_label,
      candidate_id = candidate_tbl$candidate_id[i],
      rank_index = rank_i,
      cutoff_value = cut_i,
      n_leading = n_leading,
      n_remainder = n_remainder,
      frac_leading = frac_leading,
      valid_split = isTRUE(valid_split),
      lead_score = lead_score,
      remainder_score = rem_score,
      combined_score = combined_score,
      objective_score = objective_score,
      n_standard_leading = n_standard_leading,
      n_standard_remainder = n_standard_remainder,
      nb1_mean_leading = nb1_mean_leading,
      nb1_variance_leading = nb1_variance_leading,
      nb1_iod_leading = nb1_iod_leading,
      nb1_cv2_leading = nb1_cv2_leading,
      nb1_mean_remainder = nb1_mean_remainder,
      nb1_variance_remainder = nb1_variance_remainder,
      nb1_iod_remainder = nb1_iod_remainder,
      nb1_cv2_remainder = nb1_cv2_remainder,
      nb1_pattern_pass = nb1_pattern_pass,
      nb2_mean_leading = nb2_mean_leading,
      nb2_dispersion_leading = nb2_dispersion_leading,
      nb2_iod_leading = nb2_iod_leading,
      nb2_cv2_leading = nb2_cv2_leading,
      nb2_mean_remainder = nb2_mean_remainder,
      nb2_dispersion_remainder = nb2_dispersion_remainder,
      nb2_iod_remainder = nb2_iod_remainder,
      nb2_cv2_remainder = nb2_cv2_remainder,
      nb2_pattern_pass = nb2_pattern_pass,
      dispersion_pattern_pass = dispersion_pattern_pass,
      error_message = error_message,
      stringsAsFactors = FALSE
    )
  }

  candidate_results <- dplyr::bind_rows(all_rows)
  if (!nrow(candidate_results)) {
    return(list(best_candidate = NULL, candidate_results = NULL))
  }

  valid_rows <- candidate_results[
    which(candidate_results$valid_split & is.finite(candidate_results$objective_score)),
    ,
    drop = FALSE
  ]

  if (!nrow(valid_rows)) {
    return(list(best_candidate = NULL, candidate_results = candidate_results))
  }

  valid_rows <- valid_rows[
    order(-valid_rows$objective_score, -valid_rows$combined_score, valid_rows$rank_index),
    ,
    drop = FALSE
  ]

  best_candidate <- valid_rows[1, , drop = FALSE]
  best_candidate$selected <- TRUE

  candidate_results$selected <- FALSE
  candidate_results$selected[candidate_results$candidate_id == best_candidate$candidate_id] <- TRUE

  list(best_candidate = best_candidate, candidate_results = candidate_results)
}

dataset_key_order <- c("raw_dataset", "leading_edge_dataset", "remainder_dataset")

dataset_key_labels <- c(
  raw_dataset          = "Original dataset",
  leading_edge_dataset = "Leading-edge dataset",
  remainder_dataset    = "Remainder dataset"
)

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
  # Build the plot object immediately and fail loudly if the plot code errors.
  # Returning NULL is reserved for deliberate no-data conditions inside the plot
  # functions themselves, not for swallowed build failures.
  tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) {
      stop(paste0(label, " failed: ", conditionMessage(e)))
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
  if (length(idx) == 0) stop("Could not identify condition coefficient in resultsNames(dds).")
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

if (!file.exists(count_file)) stop(sprintf("Count file not found: %s", count_file))

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
  if (!file.exists(twas_file)) stop(sprintf("TWAS file not found: %s", twas_file))
  TWAS_Seq <- read.csv(twas_file, header = TRUE, stringsAsFactors = FALSE)
  TWAS_Seq <- as.data.frame(TWAS_Seq)
  if (ncol(TWAS_Seq) < 4) stop("TWAS file must contain at least 4 columns.")
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
      )
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

# Estimate DESeq2 mean-dispersion metrics within each condition subset using
# an intercept-only model. These within-condition estimates provide the DESeq2-
# derived mean and dispersion terms used downstream to construct ranked IOD and
# CV2 for the local Fourier shared-cutoff method, while avoiding contamination
# from between-condition signal during cutoff estimation.
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
    if (!"feature_id" %in% colnames(fm)) stop("feature_metrics must contain feature_id.")
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
      loading_tbl$pc1_loading_abs >= cutoff,
      "high_loading",
      "background_loading"
    )
  } else {
    rep("background_loading", nrow(loading_tbl))
  }
  
  list(
    pca_fit              = pca_fit,
    loading_table        = loading_tbl,
    cutoff               = cutoff,
    top_n_used           = cutoff_info$top_n_actual,
    cutoff_quantile      = cutoff_info$cutoff_quantile,
    cutoff_method        = cutoff_info$method,
    cutoff_curve_df      = cutoff_info$curve_df,
    curvature_strength   = cutoff_info$curvature_strength,
    candidate_table      = cutoff_info$candidate_table,
    quantile_screen_table = cutoff_info$quantile_screen_table,
    matched_interval_table = cutoff_info$matched_interval_table,
    selected_reason      = cutoff_info$selected_reason,
    preprocessing_label  = preprocessing_label
  )
}

build_eigenvector_split <- function(count_matrix, coldata, comparison_name) {
  design_formula <- make_design_formula()

  dds_init <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData   = coldata,
    design    = design_formula
  )

  dds_init <- dds_init[rowSums(counts(dds_init)) > 0, ]
  dds_init <- estimateSizeFactors(dds_init)

  # Keep raw and normalized EVS preprocessing on the same feature universe.
  # DESeq2 drops all-zero rows before normalization, so the raw panel must be
  # restricted to those same retained features before PCA/loading ranking.
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

  primary_trt_fit <- fit_trt
  primary_untrt_fit <- fit_untrt

  combined_cutoff_info <- resolve_combined_fourier_cutoff(
    fit_trt_loading_tbl = fit_trt$loading_table,
    fit_ctrl_loading_tbl = fit_untrt$loading_table,
    fixed_top_n = evs_fixed_top_n
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
  primary_trt_fit <- fit_trt
  primary_untrt_fit <- fit_untrt
  independent_candidate_evals <- list()

  trt_high   <- as.character(subset(primary_trt_fit$loading_table,   split_class == "high_loading")$feature_id)
  untrt_high <- as.character(subset(primary_untrt_fit$loading_table, split_class == "high_loading")$feature_id)

  leading_edge_ids <- union(trt_high, untrt_high)
  analyzed_feature_ids <- union(as.character(fit_trt$loading_table$feature_id), as.character(fit_untrt$loading_table$feature_id))
  remainder_ids    <- setdiff(analyzed_feature_ids, leading_edge_ids)

  if (length(leading_edge_ids) == 0) {
    stop("Leading-edge dataset is empty. Check sample mapping or EVS cutoff settings.")
  }

  evs_cutoff_summary <- dplyr::bind_rows(
    data.frame(
      preprocessing = "normalized",
      group = "regime_shift_crossing",
      cutoff_mode = combined_cutoff_info$method,
      fixed_top_n_requested = evs_fixed_top_n,
      empiric_rank_selected = final_shared_rank,
      cutoff_quantile = fit_trt$cutoff_quantile,
      treatment_cutoff_value = fit_trt$cutoff,
      control_cutoff_value = fit_untrt$cutoff,
      composite_cutoff_value = mean(c(fit_trt$cutoff, fit_untrt$cutoff), na.rm = TRUE),
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
      treatment_cutoff_value = fit_trt_raw$cutoff,
      control_cutoff_value = fit_untrt_raw$cutoff,
      composite_cutoff_value = mean(c(fit_trt_raw$cutoff, fit_untrt_raw$cutoff), na.rm = TRUE),
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
    primary_trt_fit       = primary_trt_fit,
    primary_untrt_fit     = primary_untrt_fit,
    independent_candidate_evals = independent_candidate_evals,
    combined_cutoff_info  = combined_cutoff_info,
    evs_cutoff_summary    = evs_cutoff_summary,
    normalized_counts     = norm_counts_init,
    raw_dataset           = count_matrix,
    leading_edge_dataset  = count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset     = count_matrix[remainder_ids,    , drop = FALSE]
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

# Plot the variance explained profile for one PCA fit so the manuscript can
# show how much of the ranked loading structure is captured by the leading PCs.
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
      subtitle = compact_caption(
        paste0(
          "Group = ", pretty_group_label(group_label),
          ". Preprocessing = ", evs_preproc_short(preprocessing_label),
          ". Bars show per-component variance and the overlaid line shows cumulative variance."
        ),
        width = 72
      ),
      x = "Principal component",
      y = "Variance explained"
    ) +
    manuscript_theme() +
    theme(plot.margin = margin(t = 12, r = 14, b = 14, l = 14), legend.position = "none")
}

# PCA score-distribution histograms were intentionally removed from the manuscript
# workflow because they were redundant and are not exported.

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


# -----------------------------------------------------------------------------
# MANUSCRIPT METHODS FIGURE: RANKED IOD / CV2 CUTOFF PANEL
# -----------------------------------------------------------------------------
# plot_evs_candidate_curve() visualizes the ranked dispersion structure that led
# to the empirical cutoff for a dataset. The panel overlays smoothed IOD and
# CV2 curves, shades matched curvature intervals, marks the ranks evaluated
# within those intervals, and highlights the final selected cutoff. This plot is
# intended to show both where the candidate intervals came from and where the
# exact rank was chosen within those intervals.
# -----------------------------------------------------------------------------

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
# CORE DATASET ANALYSIS
# -----------------------------------------------------------------------------
# run_core_analysis() performs the manuscript-level inference for one dataset
# after the EVS split has been defined. It runs DESeq2 hypothesis tests, applies
# apeglm shrinkage, constructs empirical-null statistics with fdrtool, computes
# HBFSS, classifies effect tiers and method provenance, and returns the figure-
# ready result table used by the volcano, dispersion, and histogram panels.
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
  
  res <- results(dds, contrast = c("condition", "trt", "untrt"), alpha = alpha_level)
  
  res_strong <- results(
    dds,
    contrast      = c("condition", "trt", "untrt"),
    lfcThreshold  = lfc_boundary,
    altHypothesis = "greaterAbs"
  )
  
  res_weak <- results(
    dds,
    contrast      = c("condition", "trt", "untrt"),
    lfcThreshold  = lfc_boundary,
    altHypothesis = "lessAbs"
  )
  
  res_all_df            <- as.data.frame(res)
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
      )
    )
  }
  
  res_df$empirical_p <- NA_real_
  res_df$empirical_q <- NA_real_
  res_df$lfdr        <- NA_real_
  
  res_df$empirical_p[valid_stat] <- as.numeric(fdr_fit$pval)
  res_df$empirical_q[valid_stat] <- as.numeric(fdr_fit$qval)
  res_df$lfdr[valid_stat]        <- as.numeric(fdr_fit$lfdr)
  
  
  coef_name <- get_condition_coef(dds)
  shr       <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)
  shr_df    <- as.data.frame(shr)
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
  
  # HBFSS significance uses empirical-null p-values plus the resolved HC
  # threshold only. No BH correction is applied on top of HBFSS.
  hbfss_threshold_dataset  <- abs(log10(hc_p_threshold_dataset)) * lfc_boundary
  res_df$HBFSS_significant <- ifelse(
    is.na(res_df$HBFSS),
    FALSE,
    res_df$HBFSS >= hbfss_threshold_dataset
  )
  
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
  
  # Classical DESeq2 significance retained for manuscript and export outputs.
  # The manuscript uses a dual-LFC gate: padj must pass, and both the raw and
  # shrunk log2 fold changes must exceed the effect-size boundary. This is more
  # stringent than padj alone and prevents calls driven by unstable large raw
  # effects that collapse after shrinkage, or by shrinkage-only effects with no
  # matching unshrunk support.
  res_df$standard_significant <- !is.na(res_df$padj) &
    (res_df$padj < alpha_level) &
    res_df$raw_lfc_pass &
    res_df$shrunk_lfc_pass
  
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
  
  res_df$effect_class <- classify_effect_strength(
    res_df$padj_strong_effect,
    res_df$padj_weak_effect,
    alpha = alpha_level
  )
  
  res_df$resGA_padj <- res_df$padj_strong_effect
  res_df$resLA_padj <- res_df$padj_weak_effect
  
  # PDF-consistent composite-null calls for volcano classification.
  res_df$deseq2_strong_call <- !is.na(res_df$resGA_padj) &
    (res_df$resGA_padj < alpha_level) &
    res_df$shrunk_lfc_pass
  
  res_df$deseq2_weak_call <- !is.na(res_df$resLA_padj) &
    (res_df$resLA_padj < alpha_level) &
    !res_df$shrunk_lfc_pass
  
  res_df$HBFSS_only_call <- res_df$HBFSS_significant & !res_df$standard_significant
  res_df$overlap_call    <- res_df$HBFSS_significant &  res_df$standard_significant
  
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
    "pvalue", "padj", "empirical_p", "empirical_q",
    "HBFSS", "hc_p_threshold_dataset", "hbfss_threshold_dataset",
    "resLA_padj", "resGA_padj",
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
  "Neither"                 = "grey70",
  "HBFSS only"              = plot_palette$hbfss,
  "DESeq2 weak only"        = plot_palette$weak,
  "DESeq2 weak + HBFSS"     = plot_palette$overlap,
  "DESeq2 strong only"      = plot_palette$deseq2,
  "DESeq2 strong + HBFSS"   = plot_palette$overlap,
  "DESeq2 standard only"    = plot_palette$deseq2,
  "DESeq2 standard + HBFSS" = plot_palette$overlap
)

method_call_shapes <- c(
  "Neither"                 = 21,
  "HBFSS only"              = 23,
  "DESeq2 weak only"        = 22,
  "DESeq2 weak + HBFSS"     = 25,
  "DESeq2 strong only"      = 24,
  "DESeq2 strong + HBFSS"   = 24,
  "DESeq2 standard only"    = 8,
  "DESeq2 standard + HBFSS" = 23
)

method_fill_colors <- method_call_colors
method_border_colors <- c(
  "Neither"                 = "grey40",
  "HBFSS only"              = "grey15",
  "DESeq2 weak only"        = "grey15",
  "DESeq2 weak + HBFSS"     = "grey15",
  "DESeq2 strong only"      = "grey15",
  "DESeq2 strong + HBFSS"   = "grey15",
  "DESeq2 standard only"    = "grey15",
  "DESeq2 standard + HBFSS" = "grey15"
)

# -----------------------------------------------------------------------------
# PLOT HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# Build volcano-plot classes with explicit precedence for manuscript figures.
# Overlap calls take priority over strong DESeq2-only calls, which in turn take
# priority over weak DESeq2-only calls and then HBFSS-only calls. This ordering
# ensures that points called by multiple methods are colored by the combined
# evidence class rather than being overwritten by a lower-priority single-method
# label.

# -----------------------------------------------------------------------------
# VOLCANO CLASSIFICATION LOGIC
# -----------------------------------------------------------------------------
# The manuscript volcanoes separate two ideas:
# 1. effect tier, which controls broad visual emphasis (strong, weak,
#    intermediate), and
# 2. method provenance, which indicates whether support came from DESeq2,
#    HBFSS, or both.
# Weak-effect coloring is intentionally broader than a narrow lower strip: any
# feature meeting the weak-effect criterion and the manuscript significance rule
# remains in the weak tier even if it also satisfies stronger HBFSS evidence.
# -----------------------------------------------------------------------------

build_reviewer_volcano_classes <- function(df) {
  df <- as.data.frame(df)

  required_cols <- c(
    "overlap_call", "deseq2_strong_call", "deseq2_weak_call", "HBFSS_only_call",
    "standard_significant", "resLA_padj", "resGA_padj", "gene_symbol"
  )
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols)) {
    stop(
      sprintf(
        "build_reviewer_volcano_classes missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  df$method_call_class <- dplyr::case_when(
    !is.na(df$overlap_call) & df$overlap_call & !is.na(df$deseq2_weak_call) & df$deseq2_weak_call ~ "DESeq2 weak + HBFSS",
    !is.na(df$overlap_call) & df$overlap_call & !is.na(df$deseq2_strong_call) & df$deseq2_strong_call ~ "DESeq2 strong + HBFSS",
    !is.na(df$overlap_call) & df$overlap_call ~ "DESeq2 standard + HBFSS",
    !is.na(df$HBFSS_only_call) & df$HBFSS_only_call ~ "HBFSS only",
    !is.na(df$deseq2_weak_call) & df$deseq2_weak_call ~ "DESeq2 weak only",
    !is.na(df$deseq2_strong_call) & df$deseq2_strong_call ~ "DESeq2 strong only",
    !is.na(df$standard_significant) & df$standard_significant ~ "DESeq2 standard only",
    TRUE ~ "Neither"
  )

  df$method_call_class <- factor(
    df$method_call_class,
    levels = c(
      "Neither",
      "HBFSS only",
      "DESeq2 weak only",
      "DESeq2 weak + HBFSS",
      "DESeq2 strong only",
      "DESeq2 strong + HBFSS",
      "DESeq2 standard only",
      "DESeq2 standard + HBFSS"
    )
  )

  df$has_valid_gene_symbol <- !is.na(df$gene_symbol) &
    grepl('^[A-Za-z0-9._-]+$', trimws(df$gene_symbol))
  df$gene_symbol_plot <- ifelse(df$has_valid_gene_symbol, trimws(df$gene_symbol), NA_character_)
  df
}

select_volcano_labels <- function(df, y_col = "neglog10_empirical_p", n_labels = 20) {
  df <- as.data.frame(df)
  if (!nrow(df)) return(df[0, , drop = FALSE])

  df <- df[df$has_valid_gene_symbol, , drop = FALSE]
  if (!nrow(df)) return(df[0, , drop = FALSE])

  keep_levels <- c("DESeq2 weak + HBFSS", "DESeq2 strong + HBFSS", "DESeq2 standard + HBFSS", "DESeq2 strong only", "HBFSS only")
  df <- df[df$method_call_class %in% keep_levels, , drop = FALSE]
  if (!nrow(df)) return(df[0, , drop = FALSE])

  df$label_priority <- dplyr::case_when(
    df$method_call_class == "DESeq2 strong + HBFSS" ~ 1,
    df$method_call_class == "DESeq2 weak + HBFSS" ~ 2,
    df$method_call_class == "DESeq2 standard + HBFSS" ~ 3,
    df$method_call_class == "DESeq2 strong only" ~ 4,
    df$method_call_class == "HBFSS only" ~ 5,
    TRUE ~ 9
  )

  metric_y <- suppressWarnings(as.numeric(df[[y_col]]))
  metric_y[!is.finite(metric_y)] <- -Inf
  metric_h <- suppressWarnings(as.numeric(df$HBFSS))
  metric_h[!is.finite(metric_h)] <- -Inf

  ord <- order(df$label_priority, -metric_h, -metric_y, -abs(df$lfc_shrunk), na.last = TRUE)
  df  <- df[ord, , drop = FALSE]
  df  <- df[!duplicated(df$gene_symbol_plot), , drop = FALSE]
  df[seq_len(min(n_labels, nrow(df))), , drop = FALSE]
}

volcano_top_caption <- function(df) {
  hc_txt <- NA_character_
  hbfss_txt <- NA_character_
  if (nrow(df)) {
    hc_val <- suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1]))
    hbfss_val <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))
    if (is.finite(hc_val)) hc_txt <- paste0("HC p threshold=", signif(hc_val, 4))
    if (is.finite(hbfss_val)) hbfss_txt <- paste0("HBFSS cutoff=", signif(hbfss_val, 4))
  }
  parts <- c(
    hc_txt,
    hbfss_txt,
    paste0("DESeq2 standard: padj<", percent(alpha_level, accuracy = 1), " and shrunk |LFC|>=", lfc_boundary),
    paste0("DESeq2 strong: greaterAbs padj<", percent(alpha_level, accuracy = 1), " with shrunk |LFC|>=", lfc_boundary),
    paste0("DESeq2 weak: lessAbs padj<", percent(alpha_level, accuracy = 1), " with shrunk |LFC|<", lfc_boundary)
  )
  paste(parts[!is.na(parts) & nzchar(parts)], collapse = " | ")
}

volcano_count_caption <- function(df) {
  method_counts <- table(factor(df$method_call_class, levels = levels(df$method_call_class)))
  paste0(
    "Neither=", method_counts["Neither"],
    " | HBFSS=", method_counts["HBFSS only"],
    " | Weak=", method_counts["DESeq2 weak only"],
    " | Weak+HBFSS=", method_counts["DESeq2 weak + HBFSS"],
    " | Strong=", method_counts["DESeq2 strong only"],
    " | Strong+HBFSS=", method_counts["DESeq2 strong + HBFSS"],
    " | Standard=", method_counts["DESeq2 standard only"],
    " | Standard+HBFSS=", method_counts["DESeq2 standard + HBFSS"]
  )
}

.volcano_base_layers <- function() {
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
    scale_fill_manual(values = method_fill_colors, drop = FALSE, name = "Interpretive tier"),
    scale_color_manual(values = method_border_colors, drop = FALSE, name = "Interpretive tier"),
    scale_shape_manual(values = method_call_shapes, drop = FALSE, name = "Interpretive tier")
  )
}

volcano_guides <- function() {
  guides(
    color = "none",
    fill  = "none",
    shape = guide_legend(
      order = 1,
      nrow  = 2,
      byrow = TRUE,
      override.aes = list(
        size   = 3.5,
        stroke = 0.72,
        alpha  = 1,
        fill   = unname(method_fill_colors[names(method_call_shapes)]),
        colour = unname(method_border_colors[names(method_call_shapes)])
      )
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

.volcano_overlap_layer <- function(df, size_add = 0.9, stroke_add = 0.3) {
  ov <- df[!is.na(df$method_call_class) & grepl("\\+ HBFSS$", as.character(df$method_call_class)), , drop = FALSE]
  if (!nrow(ov)) return(NULL)
  geom_point(
    data        = ov,
    aes(fill = method_call_class, color = method_call_class, shape = method_call_class),
    alpha       = 1,
    size        = POINT_SIZE_PRIMARY + size_add,
    stroke      = POINT_STROKE + stroke_add,
    show.legend = FALSE
  )
}

plot_standard_volcano <- function(df, dataset_name, show_legend = TRUE) {
  df     <- build_reviewer_volcano_classes(df)
  lab_df <- select_volcano_labels(df, y_col = "neglog10_padj", n_labels = 20)

  ggplot(df, aes(lfc_shrunk, neglog10_padj)) +
    geom_point(
      aes(fill = method_call_class, color = method_call_class, shape = method_call_class),
      alpha  = POINT_ALPHA_PRIMARY,
      size   = POINT_SIZE_PRIMARY + 0.45,
      stroke = POINT_STROKE + 0.12
    ) +
    .volcano_overlap_layer(df, size_add = 1.05, stroke_add = 0.42) +
    .volcano_base_layers() +
    geom_hline(
      yintercept = -log10(alpha_level),
      linetype   = "dashed",
      linewidth  = LINE_WIDTH_BOUNDARY,
      colour     = plot_palette$threshold
    ) +
    labs(
      title   = pretty_dataset_label(dataset_name),
      x       = "Shrunken log2 fold change (β̂shrunk)",
      y       = expression(-log[10](padj)),
      caption = compact_caption(paste0(
        volcano_top_caption(), "\n",
        volcano_count_caption(df),
        "  |  LFC lines = ±", lfc_boundary
      ))
    ) +
    coord_cartesian(clip = "off") +
    plot_expand_xy() +
    manuscript_theme() +
    volcano_guides() +
    volcano_label_layer(lab_df)
}


# -----------------------------------------------------------------------------
# MANUSCRIPT HBFSS DECISION BOUNDARY
# -----------------------------------------------------------------------------
# The HBFSS volcano uses a composite boundary rather than a single hyperbola.
# A feature must satisfy both the HBFSS product threshold and the empirical-p
# floor implied by the dataset-specific HC threshold. The plotted boundary is
# therefore the upper envelope of the hyperbolic HBFSS curve and the horizontal
# empirical-p cutoff.
# -----------------------------------------------------------------------------

add_hbfss_boundary_layer <- function(p, df) {
  threshold <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))
  hc_p_threshold <- suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1]))
  if (!is.finite(threshold) || is.na(threshold) || threshold <= 0) return(p)

  finite_lfc <- suppressWarnings(as.numeric(df$lfc_shrunk))
  finite_lfc <- finite_lfc[is.finite(finite_lfc) & !is.na(finite_lfc)]
  max_abs_lfc <- max(abs(finite_lfc), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || is.na(max_abs_lfc) || max_abs_lfc <= 0) max_abs_lfc <- lfc_boundary * 3

  x_abs <- seq(from = max(0.05, min(lfc_boundary, max_abs_lfc)), to = max_abs_lfc, length.out = 400)
  y_hyper <- threshold / x_abs

  if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < 1) {
    y_floor <- -log10(hc_p_threshold)
    y_boundary <- pmax(y_hyper, y_floor)
  } else {
    y_floor <- NA_real_
    y_boundary <- y_hyper
  }

  boundary_df <- data.frame(
    lfc_shrunk = c(-rev(x_abs), x_abs),
    neglog10_empirical_p = c(rev(y_boundary), y_boundary),
    stringsAsFactors = FALSE
  )
  boundary_df <- boundary_df[is.finite(boundary_df$neglog10_empirical_p) & !is.na(boundary_df$neglog10_empirical_p), , drop = FALSE]
  if (!nrow(boundary_df)) return(p)

  p <- p + geom_path(
    data = boundary_df,
    aes(x = lfc_shrunk, y = neglog10_empirical_p),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = LINE_WIDTH_BOUNDARY,
    colour = plot_palette$threshold
  )

  if (is.finite(y_floor) && !is.na(y_floor)) {
    p <- p + geom_hline(
      yintercept = y_floor,
      linetype = "dotted",
      linewidth = LINE_WIDTH_BOUNDARY,
      colour = plot_palette$threshold
    )
  }

  p
}

plot_hbfss_volcano <- function(df, dataset_name, show_legend = TRUE) {
  df     <- build_reviewer_volcano_classes(df)
  lab_df <- select_volcano_labels(df, y_col = "neglog10_empirical_p", n_labels = 20)

  p <- ggplot(df, aes(lfc_shrunk, neglog10_empirical_p)) +
    geom_point(
      aes(fill = method_call_class, color = method_call_class, shape = method_call_class),
      alpha  = POINT_ALPHA_PRIMARY,
      size   = POINT_SIZE_PRIMARY + 0.45,
      stroke = POINT_STROKE + 0.12
    ) +
    .volcano_overlap_layer(df, size_add = 1.05, stroke_add = 0.42) +
    .volcano_base_layers()

  p <- add_hbfss_boundary_layer(p, df) +
    labs(
      title   = pretty_dataset_label(dataset_name),
      x       = "Shrunken log2 fold change (β̂shrunk)",
      y       = expression(-log[10](p[empirical])),
      caption = compact_caption(paste0(
        volcano_top_caption(),
        "\n",
        volcano_count_caption(df)
      ))
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    plot_expand_xy() +
    volcano_guides()

  if (nrow(lab_df) > 0) p <- p + volcano_label_layer(lab_df)
  p
}

plot_publication_volcano_panel <- function(df, dataset_name, show_legend = TRUE) {
  build  <- build_reviewer_volcano_classes(df)
  top_df <- select_volcano_labels(build, y_col = "neglog10_empirical_p", n_labels = 20)

  p <- ggplot(build, aes(lfc_shrunk, neglog10_empirical_p)) +
    geom_point(
      aes(fill = method_call_class, color = method_call_class, shape = method_call_class),
      alpha  = POINT_ALPHA_PRIMARY,
      size   = POINT_SIZE_PRIMARY + 0.55,
      stroke = POINT_STROKE + 0.14
    ) +
    .volcano_overlap_layer(build, size_add = 1.15, stroke_add = 0.46) +
    .volcano_base_layers()

  p <- add_hbfss_boundary_layer(p, build) +
    labs(
      title   = pretty_dataset_label(dataset_name),
      x       = "Shrunken log2 fold change (β̂shrunk)",
      y       = expression(-log[10](p[empirical])),
      caption = compact_caption(paste0(
        volcano_top_caption(),
        "\n",
        volcano_count_caption(build)
      ))
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    plot_expand_xy() +
    volcano_guides()

  if (nrow(top_df) > 0) p <- p + volcano_label_layer(top_df)
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

# -----------------------------------------------------------------------------
# MANUSCRIPT COMPARISON PANELS
# -----------------------------------------------------------------------------

# Plot the empirical-null p-value distribution and annotate the resolved HC
# threshold used downstream to define the HBFSS boundary.
plot_empirical_histogram_for_panel <- function(df, dataset_name, hc_p_threshold) {
  p <- ggplot(df, aes(empirical_p)) +
    geom_histogram(bins = HIST_BINS, fill = plot_palette$histogram, color = HIST_COLOR) +
    labs(
      title    = compact_title(pretty_dataset_label(dataset_name)),
      subtitle = paste0("HC = ", signif(hc_p_threshold, 4)),
      x = "Empirical-null p-value",
      y = "Count"
    ) +
    manuscript_theme()
  
  p <- p + geom_vline(
    xintercept = hc_p_threshold,
    color      = plot_palette$threshold,
    linewidth  = LINE_WIDTH_THRESH
  )
  
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


# -----------------------------------------------------------------------------
# DISPERSION METHODS PANEL
# -----------------------------------------------------------------------------
# plot_dispersion_panel_for_dataset() documents the DESeq2 dispersion framework
# used to derive the ranked IOD and CV2 curves. By showing gene-wise estimates,
# the fitted trend, and final shrunken dispersion values, the panel links the
# count-modeling stage to the later cutoff-selection stage.
# -----------------------------------------------------------------------------

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

  extract_shared_legend <- function(p) {
    pg <- ggplotGrob(p + theme(legend.position = "bottom"))
    gtable::gtable_filter(pg, "guide-box")
  }

  make_panel_with_shared_legend <- function(plot_list, title_text, legend_plot, ncols = length(plot_list)) {
    plot_list <- Filter(Negate(is.null), plot_list)
    if (!length(plot_list)) return(NULL)
    legend_grob <- extract_shared_legend(legend_plot)
    no_legend_list <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
    top_grob <- do.call(arrangeGrob, c(no_legend_list, list(ncol = ncols)))
    arrangeGrob(
      top_grob,
      legend_grob,
      ncol = 1,
      heights = unit.c(unit(1, "null"), unit(1.1, "in")),
      top = textGrob(title_text, gp = gpar(fontface = "bold", cex = 1.10))
    )
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

  std_grobs <- lapply(keys_present, function(k) {
    safe_panel_plot(
      plot_standard_volcano(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1],
        show_legend = FALSE
      ),
      paste0("standard volcano: ", k)
    )
  })
  std_legend_plot <- plot_standard_volcano(analysis_results[[keys_present[1]]]$results, analysis_results[[keys_present[1]]]$summary$dataset_name[1], show_legend = TRUE)
  std_panel <- make_panel_with_shared_legend(std_grobs, paste(comparison_name, "| DESeq2 volcanoes"), std_legend_plot)
  if (!is.null(std_panel)) {
    save_grob(
      std_panel,
      file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_standard_volcano_panel.png")),
      width  = 6.8 * length(Filter(Negate(is.null), std_grobs)),
      height = 6.1
    )
  }

  hbfss_grobs <- lapply(keys_present, function(k) {
    safe_panel_plot(
      plot_hbfss_volcano(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1],
        show_legend = FALSE
      ),
      paste0("HBFSS volcano: ", k)
    )
  })
  hbfss_legend_plot <- plot_hbfss_volcano(analysis_results[[keys_present[1]]]$results, analysis_results[[keys_present[1]]]$summary$dataset_name[1], show_legend = TRUE)
  hbfss_panel <- make_panel_with_shared_legend(hbfss_grobs, paste(comparison_name, "| HBFSS volcanoes"), hbfss_legend_plot)
  if (!is.null(hbfss_panel)) {
    save_grob(
      hbfss_panel,
      file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_HBFSS_volcano_panel.png")),
      width  = 6.8 * length(Filter(Negate(is.null), hbfss_grobs)),
      height = 6.1
    )
  }

  pub_grobs <- lapply(keys_present, function(k) {
    safe_panel_plot(
      plot_publication_volcano_panel(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1],
        show_legend = FALSE
      ),
      paste0("publication volcano: ", k)
    )
  })
  pub_legend_plot <- plot_publication_volcano_panel(analysis_results[[keys_present[1]]]$results, analysis_results[[keys_present[1]]]$summary$dataset_name[1], show_legend = TRUE)
  pub_panel <- make_panel_with_shared_legend(pub_grobs, paste(comparison_name, "| Publication volcanoes"), pub_legend_plot)
  if (!is.null(pub_panel)) {
    save_grob(
      pub_panel,
      file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_publication_volcano_panel.png")),
      width  = 6.8 * length(Filter(Negate(is.null), pub_grobs)),
      height = 6.1
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
        analysis_results[[k]]$results$standard_significant &
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
# =============================================================================
# SECTION 4 OF 4
# FULL PIPELINE, EXPORTS, PCA TABLES, PC1 + BASEMEAN + DISPERSION OUTPUT
# =============================================================================

# Collate feature-level PC1 loadings, EVS split membership, and DESeq2
# dispersion summaries from the full dataset plus the leading-edge and remainder
# subsets into one audit table for manuscript review and reproducibility export.
build_pc1_feature_export <- function(analysis_results, annot_df, comparison_name,
                                     evs, tab_dir) {
  annot_df$feature_id  <- as.character(annot_df$feature_id)
  annot_df$gene_symbol <- as.character(annot_df$gene_symbol)
  
  base_export <- annot_df %>%
    dplyr::distinct(feature_id, .keep_all = TRUE)
  
  add_dataset_columns <- function(df, key, prefix) {
    if (is.null(analysis_results[[key]]) || is.null(analysis_results[[key]]$results)) return(df)
    src <- analysis_results[[key]]$results
    keep_cols <- intersect(
      c("feature_id", "baseMean", "dispGeneEst", "dispFit", "dispersion"),
      names(src)
    )
    src <- src[, keep_cols, drop = FALSE]
    names(src)[names(src) == "baseMean"]    <- paste0("baseMean_", prefix)
    names(src)[names(src) == "dispGeneEst"] <- paste0("dispGeneEst_", prefix)
    names(src)[names(src) == "dispFit"]     <- paste0("dispFit_", prefix)
    names(src)[names(src) == "dispersion"]  <- paste0("dispersion_", prefix)
    dplyr::left_join(df, src, by = "feature_id")
  }
  
  out <- base_export
  out <- add_dataset_columns(out, "raw_dataset",          "original")
  out <- add_dataset_columns(out, "leading_edge_dataset", "leading_edge")
  out <- add_dataset_columns(out, "remainder_dataset",    "remainder")
  
  load_tbl <- function(tbl, col_name) {
    x <- tbl[, c("feature_id", "pc1_loading"), drop = FALSE]
    names(x)[names(x) == "pc1_loading"] <- col_name
    x
  }
  
  out <- dplyr::left_join(out, load_tbl(evs$fit_trt$loading_table,      "pc1_loading_trt_normalized"), by = "feature_id")
  out <- dplyr::left_join(out, load_tbl(evs$fit_untrt$loading_table,    "pc1_loading_ctrl_normalized"), by = "feature_id")
  out <- dplyr::left_join(out, load_tbl(evs$fit_trt_raw$loading_table,  "pc1_loading_trt_raw"),        by = "feature_id")
  out <- dplyr::left_join(out, load_tbl(evs$fit_untrt_raw$loading_table,"pc1_loading_ctrl_raw"),       by = "feature_id")
  
  out$comparison_name <- comparison_name
  
  ordered_cols <- c(
    "comparison_name",
    "feature_id",
    "gene_symbol",
    "pc1_loading_trt_normalized",
    "pc1_loading_ctrl_normalized",
    "pc1_loading_trt_raw",
    "pc1_loading_ctrl_raw",
    "baseMean_original",
    "baseMean_leading_edge",
    "baseMean_remainder",
    "dispGeneEst_original",
    "dispFit_original",
    "dispersion_original",
    "dispGeneEst_leading_edge",
    "dispFit_leading_edge",
    "dispersion_leading_edge",
    "dispGeneEst_remainder",
    "dispFit_remainder",
    "dispersion_remainder"
  )
  
  out <- out[, c(intersect(ordered_cols, names(out)), setdiff(names(out), ordered_cols)), drop = FALSE]
  
  save_csv(out, file.path(tab_dir, paste0(comparison_name, "_PC1_loadings_baseMean_dispersion_export.csv")))
  out
}


# -----------------------------------------------------------------------------
# FULL COMPARISON PIPELINE
# -----------------------------------------------------------------------------
# run_full_comparison_pipeline() is the top-level comparison orchestrator. For a
# treatment-control pair, it prepares the EVS datasets, derives dataset-specific
# cutoffs, runs downstream analysis on the original/leading-edge/remainder data
# products, assembles manuscript panels, and writes all exports for that
# comparison into the repository exports/ run directory.
# -----------------------------------------------------------------------------

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
    save_csv(evs$combined_cutoff_info$combined_wave_map, file.path(tab_dir, paste0(comparison_name, "_selected_regime_crossing_map.csv")))
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
    plot_fourier_wave_map_single(evs$combined_cutoff_info$trt_wave_obj, comparison_name, "treatment", selected_percentile = if (!is.null(evs$combined_cutoff_info$selected_crossing) && nrow(evs$combined_cutoff_info$selected_crossing)) evs$combined_cutoff_info$selected_crossing$crossing_percentile[1] else NA_real_),
    paste0(comparison_name, ": treatment local fourier wave map")
  )
  if (!is.null(trt_fourier_plot)) {
    save_grob(trt_fourier_plot, file.path(cmp_dir, paste0(comparison_name, "_treatment_local_fourier_wave_map.png")), width = 12, height = 10)
  }
  ctrl_fourier_plot <- safe_plot_build(
    plot_fourier_wave_map_single(evs$combined_cutoff_info$ctrl_wave_obj, comparison_name, "control", selected_percentile = if (!is.null(evs$combined_cutoff_info$selected_crossing) && nrow(evs$combined_cutoff_info$selected_crossing)) evs$combined_cutoff_info$selected_crossing$crossing_percentile[1] else NA_real_),
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
  
  # ---------------------------------------------------------------------------
  # NEW PCA SUMMARY TABLES AND FIGURES
  # ---------------------------------------------------------------------------
  
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
      n_overlap_significant  = sum(df$standard_significant & df$HBFSS_significant, na.rm = TRUE),
      n_strong_effect        = sum(df$effect_class == "strong_effect", na.rm = TRUE),
      n_weak_effect          = sum(df$effect_class == "weak_effect", na.rm = TRUE),
      evs_fixed_top_n        = evs_fixed_top_n,
      evs_cutoff_mode_main   = evs_cutoff_mode_main,
      evs_primary_preprocessing = "normalized",
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
    
    # Optional summary-panel pathway intentionally removed.
    # Empirical p-value histograms are retained in the cross-dataset panel workflow.
    
    analysis_results[[nm]] <- list(
      dds         = fit$dds,
      results     = df,
      summary     = summary_row,
      dataset_mat = dataset_list[[nm]],
      fig_subdir  = fig_subdir
    )
  }
  
  # ---------------------------------------------------------------------------
  # NEW COMBINED PC1 + BASEMEAN + DISPERSION EXPORT
  # ---------------------------------------------------------------------------
  
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

# Reproducibility record
session_info_txt <- capture.output(sessionInfo())
writeLines(session_info_txt, file.path(output_dir, "sessionInfo.txt"))
saveRDS(sessionInfo(), file.path(output_dir, "sessionInfo.rds"))
