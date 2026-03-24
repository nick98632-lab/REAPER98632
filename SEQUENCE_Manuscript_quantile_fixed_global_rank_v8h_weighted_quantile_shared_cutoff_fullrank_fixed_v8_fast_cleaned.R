# =============================================================================
# SECTION 1 OF 4
# SETUP, MINIMAL EMBEDDED METADATA, AND HELPERS
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

# Workflow note:
# This pipeline first builds normalized and raw EVS loading tables from the same
# nonzero feature universe, then performs a discrete constrained search over
# candidate loading cutoffs in normalized space. Each candidate cutoff is scored
# using the existing significance objective exactly as before, but a candidate is
# only eligible if it also passes the dispersion-pattern gate used by the
# manuscript workflow. That gate checks whether at least one of two variance
# parameterizations (NB1 from observed subset mean/variance, or NB2 from DESeq2
# baseMean/dispersion summaries) satisfies both directional inequalities:
#   IOD_leading > IOD_remainder
#   CV2_remainder > CV2_leading
# where IOD = variance / mean and CV2 = variance / mean^2. After filtering by
# the existing split-size constraints and this dispersion gate, each
# comparison is optimized independently. The top-ranked valid candidates from
# each comparison are pooled and the final normalized cutoff rank is chosen as
# the median pooled rank, then projected onto the raw comparison panels using
# the same retained-feature rank scale.

# -----------------------------------------------------------------------------
# USER INPUT
# -----------------------------------------------------------------------------

input_dir  <- "data"
output_root <- "exports"
analysis_name <- "EVS_HBFSS_AllComparisons_Output"

count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")

DESIGN_FORMULA <- ~ condition

# Set to TRUE and supply twas_file to enable downstream TWAS gene overlap.
run_twas_overlap <- FALSE
twas_file        <- file.path(input_dir, "3aTWAS_genes_of_11_brain_disorders.csv")

# Significance threshold for all DESeq2 calls and effect-class classification.
alpha_level <- 0.10
max_usable_hc_p_threshold <- 0.95

# HC safeguard used by HBFSS:
# When hc.thresh returns a p-threshold too close to 1, that usually indicates
# little/no meaningful departure from the empirical null. Because HBFSS scales
# its cutoff as abs(log10(hc_p_threshold_dataset)) * lfc_boundary, a near-1 HC
# threshold would make the HBFSS cutoff nearly 0 and falsely inflate calls.
# Therefore, treat HC thresholds >= max_usable_hc_p_threshold as invalid and
# fall back safely instead of using them directly.

# LFC boundary (log2 scale) used for:
# (1) strong/weak effect hypothesis tests via lfcThreshold,
# (2) the HBFSS threshold anchor: abs(log10(hc_p_threshold_dataset)) * lfc_boundary.
lfc_boundary <- 1.0

# Default number of top-ranked PAS features (by absolute PC1 loading) to retain
# when a fixed EVS cutoff is requested.

# Primary EVS pathway used for cutoff selection in the manuscript.
# The cutoff is always selected in normalized space and then projected onto the
# raw comparison panels using the same final rank.

# EVS cutoff mode.
# "shared_rank_exhaustive_search" = score a coarse grid of admissible shared ranks on the
# normalized loading tables, retain the top local ranks, and aggregate them in
# quantile space to obtain the final shared cutoff.
evs_cutoff_mode_main <- "shared_rank_exhaustive_search"

# Legacy adaptive/second-derivative cutoff logic has been removed.
# The pipeline now selects EVS cutoffs only by exhaustive shared-rank search on
# the actual downstream objective.
# Discrete exhaustive-search candidate scoring settings.
# Each comparison is optimized independently by evaluating the full ranked axis
# of admissible split points rather than only a curvature-derived candidate
# subset. A candidate rank is valid only when the leading edge has higher IOD
# than the remainder and the remainder has higher CV^2 than the leading edge.
# Valid ranks are scored by:
#   log(IOD_leading / IOD_remainder) + log(CV2_remainder / CV2_leading)
# The top valid ranks from each comparison are pooled in quantile space, and
# the final manuscript cutoff is the weighted median of that pooled anchor set.
# Anchor weights are the local objective scores, so stronger local splits carry
# more influence while all comparisons still receive the same final shared
# quantile-derived cutoff location. The leading-edge fraction is intentionally
# allowed to range across the full dataset as long as both subsets remain
# non-empty, because the cutoff is chosen by the IOD/CV^2 objective rather than
# by a hard <=50% cap.
evs_candidate_max_leading_frac      <- 1.00
evss_candidate_search_mode_comment  <- "coarse_to_fine_prefix_cache"
evss_candidate_grid_coarse_n        <- 200L
evss_candidate_refinement_window    <- 250L
evss_candidate_refinement_top_k     <- 8L
evss_export_full_candidate_grids    <- FALSE
evs_top_k_ranks_per_comparison      <- 4L

# Fast runtime controls
fast_mode                           <- TRUE
export_full_candidate_grids         <- if (isTRUE(fast_mode)) FALSE else evss_export_full_candidate_grids
export_all_plots                    <- !isTRUE(fast_mode)

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
  threshold = "#8C2D04",
  hbfss     = "#E67E22",
  deseq2    = "#C0392B",
  overlap   = "#7D3C98",
  weak      = "#4A90E2",
  histogram = "#969696",
  control   = "#4D4D4D",
  treatment = "#1F78B4"
)

output_dir <- file.path(output_root, analysis_name)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# FIGURE EXPORT CONTROLS
# -----------------------------------------------------------------------------
# Keep manuscript-essential figures by default. Redundant diagnostics removed.
export_optional_evs_variance_profiles      <- FALSE

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
# max_usable_hc_p_threshold : HC safeguard for HBFSS. Near-1 HC thresholds are treated
#   as invalid because abs(log10(hc_p_threshold_dataset)) approaches 0 as the HC threshold
#   approaches 1. In the HBFSS framework, that would artificially collapse the HBFSS cutoff
#   and inflate discoveries in weak/no-signal datasets. This is not a separate fdrtool rule;
#   it is a downstream safeguard for the HBFSS transform that uses hc.thresh output.
# HBFSS / pi_valueE       : abs(apeglm-shrunken log2FC * log10(empirical_p)).
#                           Because empirical_p lies in (0, 1], log10(empirical_p) <= 0.
#                           The absolute value is therefore essential: it converts the
#                           negative log10-scaled evidence term into a positive magnitude
#                           so larger HBFSS values reflect stronger combined effect-size
#                           and empirical-null evidence.
# standard_significant    : intentionally stricter composite rule:
#                           greaterAbs composite-null adjusted p-value
#                           (resGA_padj) < alpha plus |lfc_shrunk| >=
#                           lfc_boundary. The DESeq2 greaterAbs test is run
#                           against the MLE effect with lfcThreshold =
#                           lfc_boundary, then apeglm shrinkage is applied as
#                           an additional post-hoc effect-size gate for the
#                           exported standard_significant call.
# evs_cutoff_mode_main    : main EVS split rule. The live workflow is
#                           "shared_rank_exhaustive_search". Each normalized
#                           panel is scored over the full ranked loading table,
#                           top local ranks are retained, and the final shared
#                           cutoff is obtained by score-weighted median
#                           aggregation in quantile space.
# discrete constraint methodology : each of the two normalized EVS datasets is optimized
#                           independently by exhaustive search over the full set
#                           of admissible rank cutoffs across the ranked loading
#                           table. For any candidate rank r, the leading-edge
#                           subset is defined by the union of the treatment and
#                           control high-loading PAS features, with the only hard
#                           size constraint being that both the leading-edge and
#                           remainder subsets must remain non-empty. The
#                           objective is the sum of the panel-wise
#                           log(IOD_leading / IOD_remainder) +
#                           log(CV2_remainder / CV2_leading) scores, so the best
#                           cutoff is determined directly by the IOD/CV^2
#                           contrast rather than by a <=50% leading-edge cap.
#                           The top independently optimized normalized-dataset
#                           ranks are then combined by weighted median
#                           aggregation in quantile space to obtain the final
#                           cutoff rank, which is then inherited by the raw EVS
#                           panels.

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


safe_neglog10 <- function(x, pseudocount = 1e-12) {
  -log10(pmax(x, pseudocount))
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
  sorted_empirical_p <- sort(
    clip_probabilities(empirical_p),
    na.last    = NA,
    decreasing = FALSE
  )
  if (length(sorted_empirical_p) < 5) return(NA_real_)
  
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
  
  # HBFSS-specific HC safeguard:
  # hc.thresh may legitimately return values close to 1 when the empirical-p
  # distribution is nearly uniform and shows little/no useful departure from the
  # empirical null. That is compatible with fdrtool/hc.thresh behavior. However,
  # in this pipeline HBFSS uses abs(log10(hc_p_threshold_dataset)) * lfc_boundary
  # as the score anchor. If hc_p_threshold_dataset is too close to 1, that anchor
  # collapses toward 0 and can massively inflate HBFSS calls. Therefore, near-1
  # HC thresholds are treated as invalid for HBFSS calibration and are not used.
  if (!is.finite(out) || out <= 0 || out >= max_usable_hc_p_threshold) {
    message(sprintf(
      "[%s] HC threshold rejected for HBFSS calibration (value=%s; cutoff=%s)",
      dataset_name,
      ifelse(is.finite(out), signif(out, 6), "NA"),
      max_usable_hc_p_threshold
    ))
    return(NA_real_)
  }
  
  out
}

plot_expand_xy <- function() {
  list(
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.24))),
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.30)))
  )
}

apply_loading_cutoff <- function(fit_obj, rank_index, selected_reason = "shared_rank_selected") {
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

  rank_index <- as.integer(rank_index)[1]
  if (!is.finite(rank_index) || is.na(rank_index)) {
    stop("apply_loading_cutoff() requires rank_index to be a finite integer.")
  }
  rank_index <- min(max(1L, rank_index), nrow(loading_tbl))
  row_idx <- match(rank_index, loading_tbl$rank)
  if (is.na(row_idx) || length(row_idx) != 1L) {
    stop(sprintf("apply_loading_cutoff() could not locate rank %d.", rank_index))
  }

  cutoff_value <- as.numeric(loading_tbl$pc1_loading_abs[row_idx])
  cutoff_quantile <- 1 - (rank_index / nrow(loading_tbl))

  fit_obj$loading_table <- loading_tbl
  fit_obj$cutoff <- cutoff_value
  fit_obj$top_n_used <- rank_index
  fit_obj$cutoff_quantile <- cutoff_quantile
  fit_obj$applied_rank_index <- rank_index
  fit_obj$applied_cutoff_value <- cutoff_value
  fit_obj$cutoff_method <- "shared_rank_exhaustive_search"
  fit_obj$selected_reason <- selected_reason
  fit_obj$loading_table$split_class <- ifelse(
    fit_obj$loading_table$pc1_loading_abs >= cutoff_value,
    "high_loading",
    "background_loading"
  )
  fit_obj
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
      "[%s][%s] %d feature_id(s) not found in normalized_counts.",
      comparison_name, dataset_label, length(missing_ids)
    )
    if (isTRUE(strict)) stop(msg) else warning(msg)
  }

  if (!isTRUE(strict)) idx <- idx[!is.na(idx)]
  as.matrix(normalized_counts[idx, , drop = FALSE])
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
  empty_summary <- c(
    mean_value = NA_real_,
    variance_value = NA_real_,
    iod = NA_real_,
    cv2 = NA_real_,
    iod_mean = NA_real_,
    cv2_mean = NA_real_,
    iod_sum = NA_real_,
    cv2_sum = NA_real_,
    n_features_used = NA_real_
  )

  mat <- as.matrix(count_submatrix)
  if (is.null(mat) || !length(mat)) {
    return(empty_summary)
  }

  mat <- suppressWarnings(matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)))
  if (!nrow(mat) || !ncol(mat)) {
    return(empty_summary)
  }

  feature_mean <- rowMeans(mat, na.rm = TRUE)
  feature_var <- apply(mat, 1L, stats::var, na.rm = TRUE)
  keep <- is.finite(feature_mean) & !is.na(feature_mean) & feature_mean > 0 &
    is.finite(feature_var) & !is.na(feature_var) & feature_var >= 0

  if (!any(keep)) {
    return(empty_summary)
  }

  feature_mean <- feature_mean[keep]
  feature_var <- feature_var[keep]
  feature_iod <- feature_var / feature_mean
  feature_cv2 <- feature_var / (feature_mean ^ 2)

  finite_iod <- feature_iod[is.finite(feature_iod) & !is.na(feature_iod)]
  finite_cv2 <- feature_cv2[is.finite(feature_cv2) & !is.na(feature_cv2)]

  mean_value <- mean(feature_mean, na.rm = TRUE)
  variance_value <- mean(feature_var, na.rm = TRUE)
  iod_mean <- if (length(finite_iod)) mean(finite_iod, na.rm = TRUE) else NA_real_
  cv2_mean <- if (length(finite_cv2)) mean(finite_cv2, na.rm = TRUE) else NA_real_
  iod_sum <- if (length(finite_iod)) sum(finite_iod, na.rm = TRUE) else NA_real_
  cv2_sum <- if (length(finite_cv2)) sum(finite_cv2, na.rm = TRUE) else NA_real_

  c(
    mean_value = mean_value,
    variance_value = variance_value,
    iod = iod_sum,
    cv2 = cv2_sum,
    iod_mean = iod_mean,
    cv2_mean = cv2_mean,
    iod_sum = iod_sum,
    cv2_sum = cv2_sum,
    n_features_used = length(feature_mean)
  )
}

compute_iod_cv2_logratio_score <- function(iod_leading_mean,
                                           iod_remainder_mean,
                                           cv2_leading_mean,
                                           cv2_remainder_mean,
                                           iod_leading_sum,
                                           iod_remainder_sum,
                                           cv2_leading_sum,
                                           cv2_remainder_sum,
                                           size_weight = 1) {
  mean_values <- c(iod_leading_mean, iod_remainder_mean, cv2_leading_mean, cv2_remainder_mean)
  sum_values <- c(iod_leading_sum, iod_remainder_sum, cv2_leading_sum, cv2_remainder_sum)

  if (any(!is.finite(mean_values)) || any(is.na(mean_values)) || any(mean_values <= 0)) {
    return(list(
      objective_score = NA_real_,
      mean_logratio_score = NA_real_,
      size_weighted_logratio_score = NA_real_,
      valid_direction = FALSE
    ))
  }
  if (any(!is.finite(sum_values)) || any(is.na(sum_values)) || any(sum_values <= 0)) {
    return(list(
      objective_score = NA_real_,
      mean_logratio_score = NA_real_,
      size_weighted_logratio_score = NA_real_,
      valid_direction = FALSE
    ))
  }

  valid_direction <- (iod_leading_mean > iod_remainder_mean) && (cv2_remainder_mean > cv2_leading_mean)
  if (!valid_direction) {
    return(list(
      objective_score = NA_real_,
      mean_logratio_score = NA_real_,
      size_weighted_logratio_score = NA_real_,
      valid_direction = FALSE
    ))
  }

  mean_logratio_score <- log(iod_leading_mean / iod_remainder_mean) + log(cv2_remainder_mean / cv2_leading_mean)
  size_weighted_logratio_score <- log(iod_leading_sum / iod_remainder_sum) + log(cv2_remainder_sum / cv2_leading_sum)
  objective_score <- mean_logratio_score + (as.numeric(size_weight)[1] * size_weighted_logratio_score)

  list(
    objective_score = objective_score,
    mean_logratio_score = mean_logratio_score,
    size_weighted_logratio_score = size_weighted_logratio_score,
    valid_direction = TRUE
  )
}

build_rank_prefix_cache <- function(count_matrix,
                                    fit_obj,
                                    max_leading_frac = evs_candidate_max_leading_frac) {
  loading_tbl <- as.data.frame(fit_obj$loading_table, stringsAsFactors = FALSE)
  loading_tbl <- loading_tbl[order(loading_tbl$rank), , drop = FALSE]
  n_total <- nrow(loading_tbl)

  if (!n_total) {
    return(list(loading_table = loading_tbl, n_total = 0L, tie_group_end = integer(0),
                iod_prefix = numeric(0), cv2_prefix = numeric(0), valid_prefix = integer(0),
                iod_total = 0, cv2_total = 0, valid_total = 0L, max_leading_frac = max_leading_frac))
  }

  feature_ids <- as.character(loading_tbl$feature_id)
  idx <- match(feature_ids, rownames(count_matrix))
  if (anyNA(idx)) {
    stop(sprintf("Rank prefix cache failed: %d feature_id(s) from loading table were not found in count_matrix.", sum(is.na(idx))))
  }

  mat <- as.matrix(count_matrix[idx, , drop = FALSE])
  storage.mode(mat) <- "double"
  feature_mean <- rowMeans(mat, na.rm = TRUE)
  feature_var <- apply(mat, 1L, stats::var, na.rm = TRUE)
  valid <- is.finite(feature_mean) & !is.na(feature_mean) & (feature_mean > 0) &
    is.finite(feature_var) & !is.na(feature_var) & (feature_var >= 0)

  iod <- rep(0, n_total)
  cv2 <- rep(0, n_total)
  iod[valid] <- feature_var[valid] / feature_mean[valid]
  cv2[valid] <- feature_var[valid] / (feature_mean[valid] ^ 2)

  iod_prefix <- cumsum(iod)
  cv2_prefix <- cumsum(cv2)
  valid_prefix <- cumsum(as.integer(valid))

  load_vals <- as.numeric(loading_tbl$pc1_loading_abs)
  r <- rle(load_vals)
  tie_group_end <- rep(cumsum(r$lengths), r$lengths)

  list(
    loading_table = loading_tbl,
    n_total = as.integer(n_total),
    tie_group_end = as.integer(tie_group_end),
    iod_prefix = iod_prefix,
    cv2_prefix = cv2_prefix,
    valid_prefix = valid_prefix,
    iod_total = unname(tail(iod_prefix, 1L)),
    cv2_total = unname(tail(cv2_prefix, 1L)),
    valid_total = as.integer(tail(valid_prefix, 1L)),
    max_leading_frac = max_leading_frac
  )
}

score_single_loading_table_candidate_from_cache <- function(rank_cache,
                                                            rank_index) {
  # rank_index is the requested position on the ranked loading table.
  # When that position lands inside a tie block, the actual leading-edge size is
  # expanded to tie_group_end_rank so equally ranked features stay together.
  loading_tbl <- rank_cache$loading_table
  n_total <- as.integer(rank_cache$n_total)
  rank_index <- as.integer(rank_index)[1]

  empty_row <- function(n_leading = NA_integer_, n_remainder = NA_integer_, frac_leading = NA_real_,
                        candidate_source = "prefix_cache") {
    data.frame(
      rank_index = rank_index,
      tie_group_end_rank = if (is.na(n_leading)) NA_integer_ else as.integer(n_leading),
      n_total = n_total,
      n_leading = n_leading,
      n_remainder = n_remainder,
      frac_leading = frac_leading,
      iod_leading = NA_real_,
      iod_remainder = NA_real_,
      cv2_leading = NA_real_,
      cv2_remainder = NA_real_,
      iod_mean_leading = NA_real_,
      iod_mean_remainder = NA_real_,
      cv2_mean_leading = NA_real_,
      cv2_mean_remainder = NA_real_,
      iod_sum_leading = NA_real_,
      iod_sum_remainder = NA_real_,
      cv2_sum_leading = NA_real_,
      cv2_sum_remainder = NA_real_,
      mean_logratio_score = NA_real_,
      size_weighted_logratio_score = NA_real_,
      objective_score = NA_real_,
      valid_direction = FALSE,
      valid_split = FALSE,
      candidate_source = candidate_source,
      stringsAsFactors = FALSE
    )
  }

  if (!is.finite(rank_index) || is.na(rank_index) || rank_index < 1L || rank_index > n_total) {
    return(empty_row())
  }

  tie_group_end_rank <- as.integer(rank_cache$tie_group_end[rank_index])
  n_leading <- tie_group_end_rank
  n_remainder <- as.integer(n_total - n_leading)
  frac_leading <- n_leading / max(1L, n_total)

  if (n_leading <= 0L || n_remainder <= 0L || !is.finite(frac_leading) || frac_leading > rank_cache$max_leading_frac) {
    return(empty_row(n_leading = n_leading, n_remainder = n_remainder, frac_leading = frac_leading))
  }

  iod_sum_leading <- rank_cache$iod_prefix[n_leading]
  cv2_sum_leading <- rank_cache$cv2_prefix[n_leading]
  valid_leading <- rank_cache$valid_prefix[n_leading]

  iod_sum_remainder <- rank_cache$iod_total - iod_sum_leading
  cv2_sum_remainder <- rank_cache$cv2_total - cv2_sum_leading
  valid_remainder <- rank_cache$valid_total - valid_leading

  if (valid_leading <= 0L || valid_remainder <= 0L) {
    return(empty_row(n_leading = n_leading, n_remainder = n_remainder, frac_leading = frac_leading))
  }

  iod_mean_leading <- iod_sum_leading / valid_leading
  iod_mean_remainder <- iod_sum_remainder / valid_remainder
  cv2_mean_leading <- cv2_sum_leading / valid_leading
  cv2_mean_remainder <- cv2_sum_remainder / valid_remainder

  score_obj <- compute_iod_cv2_logratio_score(
    iod_leading_mean = iod_mean_leading,
    iod_remainder_mean = iod_mean_remainder,
    cv2_leading_mean = cv2_mean_leading,
    cv2_remainder_mean = cv2_mean_remainder,
    iod_leading_sum = iod_sum_leading,
    iod_remainder_sum = iod_sum_remainder,
    cv2_leading_sum = cv2_sum_leading,
    cv2_remainder_sum = cv2_sum_remainder
  )

  data.frame(
    rank_index = rank_index,
    tie_group_end_rank = tie_group_end_rank,
    n_total = n_total,
    n_leading = n_leading,
    n_remainder = n_remainder,
    frac_leading = frac_leading,
    iod_leading = iod_sum_leading,
    iod_remainder = iod_sum_remainder,
    cv2_leading = cv2_sum_leading,
    cv2_remainder = cv2_sum_remainder,
    iod_mean_leading = iod_mean_leading,
    iod_mean_remainder = iod_mean_remainder,
    cv2_mean_leading = cv2_mean_leading,
    cv2_mean_remainder = cv2_mean_remainder,
    iod_sum_leading = iod_sum_leading,
    iod_sum_remainder = iod_sum_remainder,
    cv2_sum_leading = cv2_sum_leading,
    cv2_sum_remainder = cv2_sum_remainder,
    mean_logratio_score = score_obj$mean_logratio_score,
    size_weighted_logratio_score = score_obj$size_weighted_logratio_score,
    objective_score = score_obj$objective_score,
    valid_direction = isTRUE(score_obj$valid_direction),
    valid_split = is.finite(score_obj$objective_score) && !is.na(score_obj$objective_score),
    candidate_source = "prefix_cache",
    stringsAsFactors = FALSE
  )
}

generate_coarse_to_fine_rank_grid <- function(shared_n_total,
                                              coarse_n = evss_candidate_grid_coarse_n,
                                              refinement_window = evss_candidate_refinement_window,
                                              refinement_top_k = evss_candidate_refinement_top_k) {
  upper_rank <- as.integer(shared_n_total - 1L)
  if (!is.finite(upper_rank) || is.na(upper_rank) || upper_rank < 1L) return(integer(0))
  if (upper_rank <= max(100L, as.integer(coarse_n))) return(seq_len(upper_rank))

  coarse_n <- max(25L, min(as.integer(coarse_n), upper_rank))
  coarse_ranks <- unique(as.integer(round(seq(1, upper_rank, length.out = coarse_n))))
  coarse_ranks <- coarse_ranks[coarse_ranks >= 1L & coarse_ranks <= upper_rank]
  coarse_ranks
}

evaluate_comparison_shared_rank_candidates <- function(count_matrix,
                                                      comparison_name,
                                                      normalized_fit_map,
                                                      max_leading_frac = evs_candidate_max_leading_frac,
                                                      top_k_ranks = evs_top_k_ranks_per_comparison) {
  n_totals <- vapply(normalized_fit_map, function(fit_obj) {
    if (is.null(fit_obj$loading_table)) return(NA_integer_)
    as.integer(nrow(fit_obj$loading_table))
  }, integer(1))
  shared_n_total <- suppressWarnings(as.integer(min(n_totals[is.finite(n_totals) & !is.na(n_totals)])))
  if (!is.finite(shared_n_total) || is.na(shared_n_total) || shared_n_total < 2L) {
    empty_summary <- data.frame(
      comparison_name = comparison_name,
      selected_rank_order = integer(0),
      independently_best_rank = integer(0),
      independently_best_local_lagrange_objective = numeric(0),
      stringsAsFactors = FALSE
    )
    return(list(
      comparison_summary = empty_summary,
      candidate_grid = NULL,
      best_rank = NA_integer_,
      best_score = NA_real_,
      top_ranks = integer(0),
      top_scores = numeric(0)
    ))
  }

  rank_caches <- lapply(normalized_fit_map, function(fit_obj) {
    build_rank_prefix_cache(
      count_matrix = count_matrix,
      fit_obj = fit_obj,
      max_leading_frac = max_leading_frac
    )
  })

  score_rank_set <- function(rank_set, source_label) {
    rank_set <- unique(as.integer(rank_set))
    rank_set <- rank_set[is.finite(rank_set) & !is.na(rank_set) & rank_set >= 1L & rank_set < shared_n_total]
    if (!length(rank_set)) return(data.frame())

    per_rank_rows <- vector("list", length(rank_set))
    for (i in seq_along(rank_set)) {
      rank_i <- rank_set[i]
      panel_rows <- lapply(names(rank_caches), function(panel_name) {
        one <- score_single_loading_table_candidate_from_cache(
          rank_cache = rank_caches[[panel_name]],
          rank_index = rank_i
        )
        one$comparison_name <- comparison_name
        one$panel_name <- panel_name
        one
      })
      panel_df <- dplyr::bind_rows(panel_rows)
      valid_panel_scores <- panel_df$objective_score[panel_df$valid_split]
      shared_score <- if (nrow(panel_df) && all(panel_df$valid_split) && length(valid_panel_scores) == nrow(panel_df)) {
        sum(valid_panel_scores)
      } else {
        NA_real_
      }
      per_rank_rows[[i]] <- data.frame(
        comparison_name = comparison_name,
        rank_index = rank_i,
        local_lagrange_objective = shared_score,
        valid_in_all_panels = is.finite(shared_score) && !is.na(shared_score),
        n_valid_panels = sum(panel_df$valid_split, na.rm = TRUE),
        panel_scores = paste(sprintf("%s=%.6f", panel_df$panel_name, panel_df$objective_score), collapse = "; "),
        candidate_source = source_label,
        stringsAsFactors = FALSE
      )
    }
    dplyr::bind_rows(per_rank_rows)
  }

  coarse_rank_union <- generate_coarse_to_fine_rank_grid(shared_n_total)
  coarse_grid <- score_rank_set(coarse_rank_union, "coarse_grid")
  coarse_valid <- coarse_grid[coarse_grid$valid_in_all_panels & is.finite(coarse_grid$local_lagrange_objective), , drop = FALSE]

  refinement_rank_union <- integer(0)
  if (nrow(coarse_valid)) {
    coarse_valid <- coarse_valid[order(-coarse_valid$local_lagrange_objective, coarse_valid$rank_index), , drop = FALSE]
    n_refine_seeds <- min(max(1L, as.integer(evss_candidate_refinement_top_k)), nrow(coarse_valid))
    refine_seeds <- as.integer(coarse_valid$rank_index[seq_len(n_refine_seeds)])
    refinement_rank_union <- unique(unlist(lapply(refine_seeds, function(seed_rank) {
      seq.int(max(1L, seed_rank - as.integer(evss_candidate_refinement_window)),
              min(shared_n_total - 1L, seed_rank + as.integer(evss_candidate_refinement_window)))
    })))
  }
  refinement_grid <- score_rank_set(setdiff(refinement_rank_union, coarse_rank_union), "coarse_to_fine_refinement")

  candidate_grid <- dplyr::bind_rows(coarse_grid, refinement_grid)
  if (nrow(candidate_grid)) {
    candidate_grid <- candidate_grid[order(candidate_grid$rank_index), , drop = FALSE]
  }
  valid_grid <- candidate_grid[candidate_grid$valid_in_all_panels & is.finite(candidate_grid$local_lagrange_objective), , drop = FALSE]

  if (!nrow(valid_grid)) {
    comparison_summary <- data.frame(
      comparison_name = comparison_name,
      selected_rank_order = integer(0),
      independently_best_rank = integer(0),
      independently_best_local_lagrange_objective = numeric(0),
      stringsAsFactors = FALSE
    )
    return(list(
      comparison_summary = comparison_summary,
      candidate_grid = candidate_grid,
      best_rank = NA_integer_,
      best_score = NA_real_,
      top_ranks = integer(0),
      top_scores = numeric(0)
    ))
  }

  valid_grid <- valid_grid[order(-valid_grid$local_lagrange_objective, valid_grid$rank_index), , drop = FALSE]
  top_k <- min(as.integer(top_k_ranks), nrow(valid_grid))
  top_rows <- valid_grid[seq_len(top_k), , drop = FALSE]

  comparison_summary <- data.frame(
    comparison_name = comparison_name,
    selected_rank_order = seq_len(nrow(top_rows)),
    independently_best_rank = as.integer(top_rows$rank_index),
    independently_best_local_lagrange_objective = as.numeric(top_rows$local_lagrange_objective),
    stringsAsFactors = FALSE
  )

  list(
    comparison_summary = comparison_summary,
    candidate_grid = candidate_grid,
    best_rank = as.integer(top_rows$rank_index[1]),
    best_score = as.numeric(top_rows$local_lagrange_objective[1]),
    top_ranks = as.integer(top_rows$rank_index),
    top_scores = as.numeric(top_rows$local_lagrange_objective)
  )
}

weighted_median_numeric <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  ok <- is.finite(x) & !is.na(x) & is.finite(w) & !is.na(w) & (w > 0)
  x <- x[ok]
  w <- w[ok]
  if (!length(x)) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cum_w <- cumsum(w) / sum(w)
  idx <- which(cum_w >= 0.5)[1]
  x[idx]
}

precompute_global_evs_rank_selection <- function(comparison_inputs,
                                                 max_leading_frac = evs_candidate_max_leading_frac,
                                                 top_k_ranks_per_comparison = evs_top_k_ranks_per_comparison) {
  per_comparison <- list()
  summary_rows <- list()

  for (cmp in names(comparison_inputs)) {
    input_obj <- comparison_inputs[[cmp]]
    design_formula <- DESIGN_FORMULA

    dds_init <- DESeqDataSetFromMatrix(
      countData = round(as.matrix(input_obj$count_matrix)),
      colData   = input_obj$coldata,
      design    = design_formula
    )
    dds_init <- dds_init[rowSums(counts(dds_init)) > 0, ]
    dds_init <- estimateSizeFactors(dds_init)

    retained_feature_ids <- rownames(dds_init)
    norm_counts_init <- as.data.frame(counts(dds_init, normalized = TRUE))
    raw_counts_init  <- as.data.frame(input_obj$count_matrix[retained_feature_ids, , drop = FALSE])

    sample_ids <- colnames(input_obj$count_matrix)
    trt_ids    <- sample_ids[input_obj$coldata$condition == "trt"]
    untrt_ids  <- sample_ids[input_obj$coldata$condition == "untrt"]

    feature_metrics_trt <- compute_condition_feature_metrics(input_obj$count_matrix[, trt_ids, drop = FALSE])
    feature_metrics_untrt <- compute_condition_feature_metrics(input_obj$count_matrix[, untrt_ids, drop = FALSE])

    fit_trt <- compute_pc1_loading_table(
      norm_counts_init,
      trt_ids,
      preprocessing_label = "Normalized before EVS",
      feature_metrics = feature_metrics_trt,
      selected_reason = "shared_rank_scoring_only"
    )

    fit_untrt <- compute_pc1_loading_table(
      norm_counts_init,
      untrt_ids,
      preprocessing_label = "Normalized before EVS",
      feature_metrics = feature_metrics_untrt,
      selected_reason = "shared_rank_scoring_only"
    )

    # Shared-rank exhaustive search is the only pathway that sets the final
    # normalized cutoff during
    # precomputation. No adaptive or second-derivative cutoff selection is used anywhere in the
    # current workflow; the final split is set only by apply_loading_cutoff()
    # after shared-rank coarse-to-fine search.

    fit_trt_raw <- compute_pc1_loading_table(
      raw_counts_init,
      trt_ids,
      preprocessing_label = "Raw counts before EVS",
      feature_metrics = feature_metrics_trt,
      selected_reason = "raw_cutoff_inherited_from_normalized"
    )

    fit_untrt_raw <- compute_pc1_loading_table(
      raw_counts_init,
      untrt_ids,
      preprocessing_label = "Raw counts before EVS",
      feature_metrics = feature_metrics_untrt,
      selected_reason = "raw_cutoff_inherited_from_normalized"
    )

    normalized_fit_map <- list(
      normalized_treatment = fit_trt,
      normalized_control = fit_untrt
    )

    comparison_eval <- evaluate_comparison_shared_rank_candidates(
      count_matrix = norm_counts_init,
      comparison_name = cmp,
      normalized_fit_map = normalized_fit_map,
      max_leading_frac = max_leading_frac,
      top_k_ranks = top_k_ranks_per_comparison
    )

    per_comparison[[cmp]] <- list(
      comparison_name = cmp,
      fit_trt = fit_trt,
      fit_untrt = fit_untrt,
      fit_trt_raw = fit_trt_raw,
      fit_untrt_raw = fit_untrt_raw,
      normalized_fit_map = normalized_fit_map,
      comparison_eval = comparison_eval,
      retained_feature_ids = retained_feature_ids,
      norm_counts_init = norm_counts_init,
      raw_counts_init = raw_counts_init,
      feature_metrics_trt = feature_metrics_trt,
      feature_metrics_untrt = feature_metrics_untrt
    )
    summary_rows[[cmp]] <- comparison_eval$comparison_summary
  }

  comparison_summary_table <- dplyr::bind_rows(summary_rows)

  anchor_rank_table <- do.call(rbind, lapply(names(per_comparison), function(cmp) {
    eval_obj <- per_comparison[[cmp]]$comparison_eval
    fit_obj <- per_comparison[[cmp]]$fit_trt
    n_total <- if (!is.null(fit_obj$loading_table)) nrow(fit_obj$loading_table) else NA_integer_
    if (length(eval_obj$top_ranks) == 0L) return(NULL)
    data.frame(
      comparison_name = cmp,
      selected_rank_order = seq_along(eval_obj$top_ranks),
      anchor_rank = as.integer(eval_obj$top_ranks),
      anchor_score = as.numeric(eval_obj$top_scores),
      n_total = as.integer(n_total),
      anchor_quantile = 1 - (as.numeric(eval_obj$top_ranks) / as.numeric(n_total)),
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(anchor_rank_table)) {
    anchor_rank_table <- data.frame(
      comparison_name = character(0),
      selected_rank_order = integer(0),
      anchor_rank = integer(0),
      anchor_score = numeric(0),
      n_total = integer(0),
      anchor_quantile = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  anchor_rank_table <- anchor_rank_table[
    is.finite(anchor_rank_table$anchor_rank) & !is.na(anchor_rank_table$anchor_rank) &
      is.finite(anchor_rank_table$anchor_score) & !is.na(anchor_rank_table$anchor_score) &
      is.finite(anchor_rank_table$anchor_quantile) & !is.na(anchor_rank_table$anchor_quantile),
    , drop = FALSE
  ]

  anchor_ranks <- anchor_rank_table$anchor_rank
  anchor_quantiles <- anchor_rank_table$anchor_quantile
  anchor_scores <- anchor_rank_table$anchor_score

  weighted_median_anchor_quantile <- if (length(anchor_quantiles)) {
    weighted_median_numeric(anchor_quantiles, anchor_scores)
  } else {
    NA_real_
  }

  comparison_n_totals <- vapply(per_comparison, function(x) {
    fit_obj <- x$fit_trt
    if (is.null(fit_obj$loading_table)) return(NA_integer_)
    as.integer(nrow(fit_obj$loading_table))
  }, integer(1))

  final_shared_rank_by_comparison <- if (length(comparison_n_totals) && is.finite(weighted_median_anchor_quantile) && !is.na(weighted_median_anchor_quantile)) {
    stats::setNames(
      vapply(comparison_n_totals, function(n_total) {
        if (!is.finite(n_total) || is.na(n_total) || n_total < 1) return(NA_integer_)
        rank_i <- as.integer(round((1 - weighted_median_anchor_quantile) * n_total))
        rank_i <- max(1L, min(as.integer(n_total), rank_i))
        rank_i
      }, integer(1)),
      names(comparison_n_totals)
    )
  } else {
    stats::setNames(rep(NA_integer_, length(comparison_n_totals)), names(comparison_n_totals))
  }

  weighted_quantile_is_interpolated_estimate <- if (length(anchor_quantiles) && is.finite(weighted_median_anchor_quantile) && !is.na(weighted_median_anchor_quantile)) {
    !(weighted_median_anchor_quantile %in% anchor_quantiles)
  } else {
    NA
  }

  global_search_table <- anchor_rank_table
  if (nrow(global_search_table)) {
    global_search_table$anchor_pool_size <- nrow(global_search_table)
    global_search_table$weighted_median_anchor_quantile <- weighted_median_anchor_quantile
    global_search_table$weighted_quantile_is_interpolated_estimate <- weighted_quantile_is_interpolated_estimate
    global_search_table$aggregation_method <- "score_weighted_median_of_top_k_local_lagrange_anchor_quantiles"
    global_search_table$final_shared_rank_for_comparison <- unname(final_shared_rank_by_comparison[global_search_table$comparison_name])
  }

  final_summary <- comparison_summary_table
  if (nrow(final_summary)) {
    final_summary$top_k_ranks_per_comparison <- as.integer(top_k_ranks_per_comparison)
    final_summary$weighted_median_anchor_quantile <- weighted_median_anchor_quantile
    final_summary$weighted_quantile_is_interpolated_estimate <- weighted_quantile_is_interpolated_estimate
    final_summary$final_shared_rank <- unname(final_shared_rank_by_comparison[final_summary$comparison_name])
    final_summary$aggregation_method <- "score_weighted_median_of_top_k_local_lagrange_anchor_quantiles"
    final_summary$candidate_source <- "full_rank_exhaustive_search"
    final_summary$n_anchor_ranks_used <- nrow(anchor_rank_table)
  }

  list(
    final_shared_rank = if (length(final_shared_rank_by_comparison)) as.integer(stats::median(final_shared_rank_by_comparison[is.finite(final_shared_rank_by_comparison) & !is.na(final_shared_rank_by_comparison)])) else NA_integer_,
    final_shared_quantile = weighted_median_anchor_quantile,
    weighted_median_anchor_quantile = weighted_median_anchor_quantile,
    weighted_quantile_is_interpolated_estimate = weighted_quantile_is_interpolated_estimate,
    final_shared_rank_by_comparison = final_shared_rank_by_comparison,
    comparison_n_totals = comparison_n_totals,
    top_k_ranks_per_comparison = as.integer(top_k_ranks_per_comparison),
    comparison_summary_table = final_summary,
    global_search_table = global_search_table,
    per_comparison = per_comparison,
    n_anchor_ranks_used = nrow(anchor_rank_table)
  )
}

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
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.csv(df, file = path, row.names = FALSE)
}

save_grob <- function(g, path, width = 16.4, height = 9.9, dpi = figure_dpi, bg = "white") {
  if (!isTRUE(export_all_plots)) return(invisible(FALSE))
  tryCatch(
    {
      dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
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
    expr,
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
# an intercept-only model. This keeps the EVS NB2 variance curve anchored to
# within-condition count structure instead of mixing in between-condition signal,
# which would blur the variance profile that the leading-edge split is meant to
# detect.
compute_condition_feature_metrics <- function(count_submatrix) {
  cd <- S4Vectors::DataFrame(row.names = colnames(count_submatrix))
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_submatrix)),
    colData   = cd,
    design    = ~ 1
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, quiet = TRUE)

  md <- as.data.frame(SummarizedExperiment::mcols(dds), stringsAsFactors = FALSE)
  md$feature_id <- rownames(md)
  keep_cols <- intersect(c("feature_id", "baseMean", "dispGeneEst", "dispFit", "dispersion"), colnames(md))
  md <- md[, keep_cols, drop = FALSE]
  md
}

compute_pc1_loading_table <- function(value_df, sample_names,
                                      preprocessing_label = "Normalized before EVS",
                                      feature_metrics = NULL,
                                      selected_reason = "shared_rank_external") {
  # Both EVS arms define PCA/loading geometry on a comparable log1p-transformed
  # scale. The downstream shared-rank IOD/CV² scoring still operates on
  # normalized count-scale data.
  x <- log1p(as.matrix(value_df[, sample_names, drop = FALSE]))

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

  loading_tbl <- loading_tbl[order(loading_tbl$pc1_loading_abs, decreasing = TRUE), , drop = FALSE]
  loading_tbl$rank <- seq_len(nrow(loading_tbl))
  loading_tbl$split_class <- "background_loading"

  list(
    pca_fit               = pca_fit,
    loading_table         = loading_tbl,
    cutoff                = NA_real_,
    top_n_used            = NA_integer_,
    cutoff_quantile       = NA_real_,
    cutoff_method         = "shared_rank_exhaustive_search",
    selected_reason       = selected_reason,
    preprocessing_label   = preprocessing_label
  )
}

build_eigenvector_split <- function(count_matrix,
                                    coldata,
                                    comparison_name,
                                    final_shared_rank = NA_integer_,
                                    precomputed_selection = NULL) {
  precomputed_cmp <- NULL
  if (!is.null(precomputed_selection$per_comparison) &&
      !is.null(precomputed_selection$per_comparison[[comparison_name]])) {
    precomputed_cmp <- precomputed_selection$per_comparison[[comparison_name]]
  }

  if (!is.null(precomputed_cmp)) {
    retained_feature_ids <- as.character(precomputed_cmp$retained_feature_ids)
    norm_counts_init <- as.data.frame(precomputed_cmp$norm_counts_init)
    raw_counts_init <- as.data.frame(precomputed_cmp$raw_counts_init)
    feature_metrics_trt <- precomputed_cmp$feature_metrics_trt
    feature_metrics_untrt <- precomputed_cmp$feature_metrics_untrt

    sample_ids <- colnames(count_matrix)
    trt_ids    <- sample_ids[coldata$condition == "trt"]
    untrt_ids  <- sample_ids[coldata$condition == "untrt"]

    fit_trt <- compute_pc1_loading_table(
      norm_counts_init,
      trt_ids,
      preprocessing_label = "Normalized before EVS",
      feature_metrics = feature_metrics_trt
    )

    fit_untrt <- compute_pc1_loading_table(
      norm_counts_init,
      untrt_ids,
      preprocessing_label = "Normalized before EVS",
      feature_metrics = feature_metrics_untrt
    )

    fit_trt_raw <- compute_pc1_loading_table(
      raw_counts_init,
      trt_ids,
      preprocessing_label = "Raw counts before EVS",
      feature_metrics = feature_metrics_trt,
      selected_reason = "raw_cutoff_inherited_from_normalized"
    )

    fit_untrt_raw <- compute_pc1_loading_table(
      raw_counts_init,
      untrt_ids,
      preprocessing_label = "Raw counts before EVS",
      feature_metrics = feature_metrics_untrt,
      selected_reason = "raw_cutoff_inherited_from_normalized"
    )

    normalized_fit_map <- list(
      normalized_treatment = fit_trt,
      normalized_control = fit_untrt
    )

    comparison_eval <- precomputed_cmp$comparison_eval
  } else {
    warning("build_eigenvector_split() entered fallback recomputation path because precomputed_selection was NULL or missing for this comparison. This path is retained as a defensive fallback and is not used in the normal pipeline.")

    design_formula <- DESIGN_FORMULA

    dds_init <- DESeqDataSetFromMatrix(
      countData = round(as.matrix(count_matrix)),
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
      preprocessing_label = "Normalized before EVS",
      feature_metrics = feature_metrics_trt
    )

    fit_untrt <- compute_pc1_loading_table(
      norm_counts_init,
      untrt_ids,
      preprocessing_label = "Normalized before EVS",
      feature_metrics = feature_metrics_untrt
    )

    fit_trt_raw <- compute_pc1_loading_table(
      raw_counts_init,
      trt_ids,
      preprocessing_label = "Raw counts before EVS",
      feature_metrics = feature_metrics_trt,
      selected_reason = "raw_cutoff_inherited_from_normalized"
    )

    fit_untrt_raw <- compute_pc1_loading_table(
      raw_counts_init,
      untrt_ids,
      preprocessing_label = "Raw counts before EVS",
      feature_metrics = feature_metrics_untrt,
      selected_reason = "raw_cutoff_inherited_from_normalized"
    )

    normalized_fit_map <- list(
      normalized_treatment = fit_trt,
      normalized_control   = fit_untrt
    )

    comparison_eval <- evaluate_comparison_shared_rank_candidates(
      count_matrix = norm_counts_init,
      comparison_name = comparison_name,
      normalized_fit_map = normalized_fit_map,
      max_leading_frac = evs_candidate_max_leading_frac
    )
  }

  independently_best_rank <- comparison_eval$best_rank
  if (!is.finite(independently_best_rank) || is.na(independently_best_rank)) {
    independently_best_rank <- NA_integer_
  }

  if (!is.null(precomputed_selection) && !is.null(precomputed_selection$final_shared_rank_by_comparison)) {
    shared_rank_lookup <- precomputed_selection$final_shared_rank_by_comparison
    if (!is.null(names(shared_rank_lookup)) && comparison_name %in% names(shared_rank_lookup)) {
      final_shared_rank <- as.integer(shared_rank_lookup[[comparison_name]])
    }
  }

  if (!is.finite(final_shared_rank) || is.na(final_shared_rank)) {
    final_shared_rank <- independently_best_rank
  }
  # The final shared cutoff is now aggregated in quantile space and then
  # projected back to each comparison's ranked feature universe.

  final_rank_aggregation <- data.frame(
    comparison_name = comparison_name,
    independently_best_rank = independently_best_rank,
    independently_best_local_lagrange_objective = comparison_eval$best_score,
    weighted_median_anchor_quantile = if (!is.null(precomputed_selection$weighted_median_anchor_quantile)) precomputed_selection$weighted_median_anchor_quantile else NA_real_,
    weighted_quantile_is_interpolated_estimate = if (!is.null(precomputed_selection$weighted_quantile_is_interpolated_estimate)) precomputed_selection$weighted_quantile_is_interpolated_estimate else NA,
    final_shared_rank = final_shared_rank,
    top_k_ranks_per_comparison = if (!is.null(precomputed_selection$top_k_ranks_per_comparison)) precomputed_selection$top_k_ranks_per_comparison else evs_top_k_ranks_per_comparison,
    aggregation_method = "score_weighted_median_of_top_k_local_lagrange_anchor_quantiles",
    n_anchor_ranks_used = if (!is.null(precomputed_selection$n_anchor_ranks_used)) precomputed_selection$n_anchor_ranks_used else NA_integer_,
    stringsAsFactors = FALSE
  )

  if (!is.na(final_shared_rank)) {
    final_reason <- "score_weighted_median_of_top_k_local_lagrange_anchor_quantiles_shared_rank_selected"

    fit_trt <- apply_loading_cutoff(
      fit_trt,
      rank_index = final_shared_rank,
      selected_reason = final_reason
    )
    fit_untrt <- apply_loading_cutoff(
      fit_untrt,
      rank_index = final_shared_rank,
      selected_reason = final_reason
    )
    fit_trt_raw <- apply_loading_cutoff(
      fit_trt_raw,
      rank_index = final_shared_rank,
      selected_reason = paste0(final_reason, "_projected_to_raw")
    )
    fit_untrt_raw <- apply_loading_cutoff(
      fit_untrt_raw,
      rank_index = final_shared_rank,
      selected_reason = paste0(final_reason, "_projected_to_raw")
    )

    fit_trt_raw$cutoff_method <- "normalized_global_rank_projected_to_raw_selected"
    fit_untrt_raw$cutoff_method <- "normalized_global_rank_projected_to_raw_selected"
  }

  primary_trt_fit <- fit_trt
  primary_untrt_fit <- fit_untrt

  independent_candidate_evals <- list(
    comparison_shared_search = list(
      best_candidate = if (!is.na(independently_best_rank)) subset(comparison_eval$candidate_grid, rank_index == independently_best_rank) else NULL,
      candidate_results = comparison_eval$candidate_grid
    )
  )

  trt_high   <- as.character(subset(primary_trt_fit$loading_table,   split_class == "high_loading")$feature_id)
  untrt_high <- as.character(subset(primary_untrt_fit$loading_table, split_class == "high_loading")$feature_id)

  leading_edge_ids <- union(trt_high, untrt_high)
  analyzed_feature_ids <- union(as.character(fit_trt$loading_table$feature_id), as.character(fit_untrt$loading_table$feature_id))
  remainder_ids    <- setdiff(analyzed_feature_ids, leading_edge_ids)

  if (length(leading_edge_ids) == 0) {
    stop("Leading-edge dataset is empty. Check sample mapping or EVS cutoff settings.")
  }

  # Export a single cutoff summary row per preprocessing × group combination.
  # After the final cutoff is projected back onto the normalized loading tables,
  # the "primary" and "comparison_panel" normalized fits converge by design.
  # Collapsing to four rows avoids exporting duplicate normalized summaries.
  evs_cutoff_summary <- dplyr::bind_rows(
    data.frame(
      preprocessing = "normalized",
      group = "treatment",
      cutoff_mode = fit_trt$cutoff_method,
      empiric_rank_selected = fit_trt$top_n_used,
      cutoff_quantile = fit_trt$cutoff_quantile,
      selected_reason = fit_trt$selected_reason,
      stringsAsFactors = FALSE
    ),
    data.frame(
      preprocessing = "normalized",
      group = "control",
      cutoff_mode = fit_untrt$cutoff_method,
      empiric_rank_selected = fit_untrt$top_n_used,
      cutoff_quantile = fit_untrt$cutoff_quantile,
      selected_reason = fit_untrt$selected_reason,
      stringsAsFactors = FALSE
    ),
    data.frame(
      preprocessing = "raw",
      group = "treatment",
      cutoff_mode = fit_trt_raw$cutoff_method,
      empiric_rank_selected = fit_trt_raw$top_n_used,
      cutoff_quantile = fit_trt_raw$cutoff_quantile,
      selected_reason = fit_trt_raw$selected_reason,
      stringsAsFactors = FALSE
    ),
    data.frame(
      preprocessing = "raw",
      group = "control",
      cutoff_mode = fit_untrt_raw$cutoff_method,
      empiric_rank_selected = fit_untrt_raw$top_n_used,
      cutoff_quantile = fit_untrt_raw$cutoff_quantile,
      selected_reason = fit_untrt_raw$selected_reason,
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
    final_rank_aggregation = final_rank_aggregation,
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
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition") +
    scale_fill_manual(values  = condition_fills,  labels = condition_labels, name = "Condition") +
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
    guides(
      fill  = guide_legend(override.aes = list(size = 3.0, shape = 21, colour = "white")),
      shape = guide_legend(override.aes = list(size = 3.0, fill  = "grey70", colour = "white"))
    )
}

plot_pc1_loading_rank <- function(loading_tbl, cutoff, dataset_label, group_label,
                                  top_n_used = NA_integer_, cutoff_quantile = NA_real_,
                                  preprocessing_label = "Normalized before EVS") {
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
        "Rank ", ifelse(is.finite(top_n_used), top_n_used, NA_integer_), " cutoff = ", signif(cutoff, 4), "\n",
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
          preprocessing_label,
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

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
  design_formula <- DESIGN_FORMULA
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_mat)),
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
  
  message(
    sprintf(
      "[%s] Mean Wald stat sent to fdrtool: %.4f  (computed only on finite DESeq2 Wald statistics)",
      dataset_name, mean(stat_vec, na.rm = TRUE)
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
  
  # HBFSS significance uses empirical-null p-values plus the HC threshold only.
  # No BH correction is applied on top of HBFSS. Near-1 HC thresholds are treated
  # as invalid for HBFSS calibration because abs(log10(hc_p_threshold_dataset))
  # would become nearly 0 and would make the HBFSS cutoff spuriously permissive.
  if (is.na(hc_p_threshold_dataset)) {
    hbfss_threshold_dataset  <- NA_real_
    res_df$passes_hc_p_gate  <- FALSE
    res_df$HBFSS_significant <- FALSE
  } else {
    hbfss_threshold_dataset  <- abs(log10(hc_p_threshold_dataset)) * lfc_boundary
    res_df$passes_hc_p_gate  <- !is.na(res_df$empirical_p) &
      (res_df$empirical_p < hc_p_threshold_dataset)
    res_df$HBFSS_significant <- res_df$passes_hc_p_gate &
      !is.na(res_df$HBFSS) &
      (res_df$HBFSS >= hbfss_threshold_dataset)
  }
  
  res_df$regulation_direction <- ifelse(
    is.na(res_df$lfc_shrunk),
    NA_character_,
    ifelse(
      res_df$lfc_shrunk > 0, "upregulated",
      ifelse(res_df$lfc_shrunk < 0, "downregulated", "no_change")
    )
  )
  
  # Audit-only column retained for export; this does not gate standard_significant.
  res_df$raw_lfc_pass    <- !is.na(res_df$log2FoldChange) & (abs(res_df$log2FoldChange) >= lfc_boundary)
  res_df$shrunk_lfc_pass <- !is.na(res_df$lfc_shrunk)     & (abs(res_df$lfc_shrunk)     >= lfc_boundary)

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
  
  # Standard significance intentionally uses a stricter two-stage rule:
  # the greaterAbs composite-null adjusted p-value is the primary inferential
  # call, and apeglm-shrunken LFC is then applied as an additional post-hoc
  # effect-size direction/magnitude filter for the exported call.
  res_df$standard_significant <- !is.na(res_df$resGA_padj) &
    (res_df$resGA_padj < alpha_level) &
    res_df$shrunk_lfc_pass
  
  # PDF-consistent composite-null calls for volcano classification.
  # Retain deseq2_strong_call as an explicit export alias of the manuscript's
  # standard_significant field rather than recomputing the same expression.
  res_df$deseq2_strong_call <- res_df$standard_significant
  
  res_df$deseq2_weak_call <- !is.na(res_df$resLA_padj) &
    (res_df$resLA_padj < alpha_level) &
    !res_df$shrunk_lfc_pass
  
  res_df$HBFSS_only_call <- res_df$HBFSS_significant & !res_df$standard_significant
  res_df$overlap_call    <- res_df$HBFSS_significant &  res_df$standard_significant
  
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
    "standard_significant", "passes_hc_p_gate", "HBFSS_significant", "effect_class"
  )
  
  final_df <- final_df[, c(intersect(preferred_cols, names(final_df)),
                           setdiff(names(final_df), preferred_cols)), drop = FALSE]
  
  list(
    dds             = dds,
    results         = final_df,
    hc_p_threshold  = hc_p_threshold_dataset,
    hbfss_threshold = hbfss_threshold_dataset
  )
}

# -----------------------------------------------------------------------------
# COLOUR / SHAPE SCALES
# -----------------------------------------------------------------------------

method_call_colors <- c(
  "Neither"            = "grey70",
  "DESeq2 weak only"   = plot_palette$weak,
  "DESeq2 strong only" = plot_palette$deseq2,
  "HBFSS only"         = plot_palette$hbfss,
  "Overlap"            = plot_palette$overlap
)

method_call_shapes <- c(
  "Neither"            = 21,
  "DESeq2 weak only"   = 22,
  "DESeq2 strong only" = 24,
  "HBFSS only"         = 23,
  "Overlap"            = 25
)

method_fill_colors <- method_call_colors
method_border_colors <- c(
  "Neither"            = "grey40",
  "DESeq2 weak only"   = "grey15",
  "DESeq2 strong only" = "grey15",
  "HBFSS only"         = "grey15",
  "Overlap"            = "grey15"
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
build_reviewer_volcano_classes <- function(df) {
  df <- as.data.frame(df)

  required_cols <- c(
    "overlap_call", "deseq2_strong_call", "deseq2_weak_call", "HBFSS_only_call",
    "gene_symbol"
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
    !is.na(df$overlap_call)      & df$overlap_call      ~ "Overlap",
    !is.na(df$deseq2_strong_call)& df$deseq2_strong_call~ "DESeq2 strong only",
    !is.na(df$deseq2_weak_call)  & df$deseq2_weak_call  ~ "DESeq2 weak only",
    !is.na(df$HBFSS_only_call)   & df$HBFSS_only_call   ~ "HBFSS only",
    TRUE ~ "Neither"
  )

  df$method_call_class <- factor(
    df$method_call_class,
    levels = c("Neither", "DESeq2 weak only", "DESeq2 strong only", "HBFSS only", "Overlap")
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

  df <- df[df$method_call_class %in% c("Overlap", "DESeq2 strong only", "HBFSS only"), , drop = FALSE]
  if (!nrow(df)) return(df[0, , drop = FALSE])

  df$label_priority <- dplyr::case_when(
    df$method_call_class == "Overlap"            ~ 1,
    df$method_call_class == "DESeq2 strong only" ~ 2,
    df$method_call_class == "HBFSS only"         ~ 3,
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

volcano_top_caption <- function() {
  paste0("DESeq2: padj < ", percent(alpha_level, accuracy = 1), ". HBFSS: empirical p + HC threshold.")
}

volcano_count_caption <- function(df) {
  method_counts <- table(factor(df$method_call_class, levels = levels(df$method_call_class)))
  paste0(
    "Neither=", method_counts["Neither"],
    " | Weak=", method_counts["DESeq2 weak only"],
    " | Strong=", method_counts["DESeq2 strong only"],
    " | HBFSS=", method_counts["HBFSS only"],
    " | Overlap=", method_counts["Overlap"]
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
        fill   = unname(method_fill_colors),
        colour = unname(method_border_colors)
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
  ov <- df[!is.na(df$method_call_class) & df$method_call_class == "Overlap", , drop = FALSE]
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

plot_standard_volcano <- function(df, dataset_name) {
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

add_hbfss_boundary_layer <- function(p, df) {
  threshold <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))
  if (!is.finite(threshold) || is.na(threshold) || threshold <= 0) return(p)

  finite_lfc <- suppressWarnings(as.numeric(df$lfc_shrunk))
  finite_lfc <- finite_lfc[is.finite(finite_lfc) & !is.na(finite_lfc)]
  max_abs_lfc <- max(abs(finite_lfc), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || is.na(max_abs_lfc) || max_abs_lfc <= 0) max_abs_lfc <- lfc_boundary * 3

  # HBFSS = |lfc_shrunk * log10(empirical_p)|. Holding the threshold fixed
  # gives |log10(empirical_p)| = threshold / |lfc_shrunk|, so the decision
  # boundary drawn in volcano space is a rectangular hyperbola.
  x_abs <- seq(from = max(0.05, min(lfc_boundary, max_abs_lfc)), to = max_abs_lfc, length.out = 400)
  boundary_df <- data.frame(
    lfc_shrunk = c(-rev(x_abs), x_abs),
    neglog10_empirical_p = c(rev(threshold / x_abs), threshold / x_abs),
    stringsAsFactors = FALSE
  )
  boundary_df <- boundary_df[is.finite(boundary_df$neglog10_empirical_p) & !is.na(boundary_df$neglog10_empirical_p), , drop = FALSE]
  if (!nrow(boundary_df)) return(p)

  p + geom_path(
    data = boundary_df,
    aes(x = lfc_shrunk, y = neglog10_empirical_p),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = LINE_WIDTH_BOUNDARY,
    colour = plot_palette$threshold
  )
}

plot_hbfss_volcano_panel <- function(df, dataset_name,
                                    point_size_add = 0.45,
                                    point_stroke_add = 0.12,
                                    overlap_size_add = 1.05,
                                    overlap_stroke_add = 0.42) {
  df     <- build_reviewer_volcano_classes(df)
  lab_df <- select_volcano_labels(df, y_col = "neglog10_empirical_p", n_labels = 20)

  p <- ggplot(df, aes(lfc_shrunk, neglog10_empirical_p)) +
    geom_point(
      aes(fill = method_call_class, color = method_call_class, shape = method_call_class),
      alpha  = POINT_ALPHA_PRIMARY,
      size   = POINT_SIZE_PRIMARY + point_size_add,
      stroke = POINT_STROKE + point_stroke_add
    ) +
    .volcano_overlap_layer(df, size_add = overlap_size_add, stroke_add = overlap_stroke_add) +
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

# -----------------------------------------------------------------------------
# MANUSCRIPT COMPARISON PANELS
# -----------------------------------------------------------------------------

compute_dataset_pca_plot <- function(count_df, coldata, dataset_name,
                                     preprocessing = c("normalized", "raw_counts"),
                                     precomputed_matrix = NULL) {
  preprocessing <- match.arg(preprocessing)

  if (identical(preprocessing, "normalized")) {
    if (is.null(precomputed_matrix)) {
      stop("compute_dataset_pca_plot() requires precomputed_matrix when preprocessing='normalized'.")
    }
    x <- as.matrix(precomputed_matrix)
  } else {
    x <- log1p(as.matrix(count_df))
  }
  
  pca_fit     <- prcomp(t(x), scale. = FALSE, rank. = 2)
  pca_var     <- pca_fit$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)
  
  pca_df <- data.frame(
    Sample           = rownames(pca_fit$x),
    PC1              = pca_fit$x[, 1],
    PC2              = pca_fit$x[, 2],
    Condition        = as.character(coldata[match(rownames(pca_fit$x), rownames(coldata)), "condition"]),
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
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition") +
    scale_fill_manual(values  = condition_fills,  labels = condition_labels, name = "Condition") +
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
    guides(
      fill  = guide_legend(override.aes = list(size = 3.0, shape = 21, colour = "white")),
      shape = guide_legend(override.aes = list(size = 3.0, fill  = "grey70", colour = "white"))
    )
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
  keys_present <- names(analysis_results)
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
        analysis_results[[k]]$summary$dataset_name[1]
      ),
      paste0("standard volcano: ", k)
    )
  })
  std_panel <- make_panel(std_grobs, paste(comparison_name, "| DESeq2 volcanoes"))
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
      plot_hbfss_volcano_panel(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1],
        point_size_add = 0.45,
        point_stroke_add = 0.12,
        overlap_size_add = 1.05,
        overlap_stroke_add = 0.42
      ),
      paste0("HBFSS volcano: ", k)
    )
  })
  hbfss_panel <- make_panel(hbfss_grobs, paste(comparison_name, "| HBFSS volcanoes"))
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
      plot_hbfss_volcano_panel(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1],
        point_size_add = 0.55,
        point_stroke_add = 0.14,
        overlap_size_add = 1.15,
        overlap_stroke_add = 0.44
      ),
      paste0("publication volcano: ", k)
    )
  })
  pub_panel <- make_panel(pub_grobs, paste(comparison_name, "| Publication volcanoes"))
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

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df, final_shared_rank = NA_integer_, precomputed_selection = NULL) {
  cmp_dir     <- file.path(output_dir, comparison_name)
  tab_dir     <- file.path(cmp_dir, "tables")
  fig_raw_dir <- file.path(cmp_dir, "figures_raw")
  fig_le_dir  <- file.path(cmp_dir, "figures_leading_edge")
  fig_rem_dir <- file.path(cmp_dir, "figures_remainder")
  
  for (d in c(cmp_dir, tab_dir, fig_raw_dir, fig_le_dir, fig_rem_dir)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
  
  evs <- build_eigenvector_split(count_matrix, coldata, comparison_name, final_shared_rank = final_shared_rank, precomputed_selection = precomputed_selection)

  evs_cutoff_summary_out <- evs$evs_cutoff_summary
  evs_cutoff_summary_out$comparison_name <- comparison_name
  save_csv(evs_cutoff_summary_out, file.path(tab_dir, paste0(comparison_name, "_EVS_cutoff_summary.csv")))
  

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
  
  # ---------------------------------------------------------------------------
  # PCA SUMMARY TABLES AND FIGURES
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

  if (!is.null(evs$independent_candidate_evals)) {
    for (nm in names(evs$independent_candidate_evals)) {
      cand_res <- evs$independent_candidate_evals[[nm]]$candidate_results
      if (!is.null(cand_res) && nrow(cand_res) && isTRUE(export_full_candidate_grids)) {
        save_csv(cand_res, file.path(tab_dir, paste0(comparison_name, "_EVS_", nm, "_candidate_grid.csv")))
      }
    }
  }
  if (!is.null(evs$final_rank_aggregation) && nrow(evs$final_rank_aggregation)) {
    save_csv(evs$final_rank_aggregation, file.path(tab_dir, paste0(comparison_name, "_EVS_independent_rank_aggregation_summary.csv")))
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
      evs_cutoff_mode_main   = evs_cutoff_mode_main,
      evs_primary_preprocessing = "normalized",
      trt_top_n_used         = evs$fit_trt$top_n_used,
      ctrl_top_n_used        = evs$fit_untrt$top_n_used,
      trt_cutoff_method      = evs$fit_trt$cutoff_method,
      ctrl_cutoff_method     = evs$fit_untrt$cutoff_method,
      trt_cutoff_quantile    = evs$fit_trt$cutoff_quantile,
      ctrl_cutoff_quantile   = evs$fit_untrt$cutoff_quantile,
      trt_selected_reason    = evs$primary_trt_fit$selected_reason,
      ctrl_selected_reason   = evs$primary_untrt_fit$selected_reason,
      stringsAsFactors       = FALSE
    )
    
    save_csv(summary_row, file.path(tab_dir, paste0(full_dataset_name, "_summary.csv")))
    
    
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
      dds        = fit$dds,
      results    = df,
      summary    = summary_row,
      fig_subdir = fig_subdir
    )
  }
  
  # ---------------------------------------------------------------------------
  # COMBINED PC1 + BASEMEAN + DISPERSION EXPORT
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

global_evs_selection <- precompute_global_evs_rank_selection(comparison_inputs)

if (!is.null(global_evs_selection$comparison_summary_table) &&
    nrow(global_evs_selection$comparison_summary_table) > 0) {
  save_csv(
    global_evs_selection$comparison_summary_table,
    file.path(output_dir, "EVS_global_shared_rank_summary.csv")
  )
}

if (!is.null(global_evs_selection$global_search_table) &&
    nrow(global_evs_selection$global_search_table) > 0) {
  save_csv(
    global_evs_selection$global_search_table,
    file.path(output_dir, "EVS_global_shared_rank_search.csv")
  )
}

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
      annot_df        = OrigID_Symbol,
      final_shared_rank = global_evs_selection$final_shared_rank,
      precomputed_selection = global_evs_selection
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
