# =============================================================================
# EVS + HBFSS FULL PIPELINE
# DESEQ2-INFORMED NB REGIME CUT-OFF VERSION
# REVISED: COMPARISON-SPECIFIC COMBINED RANK AXIS
# FIXED: FEATURE-ID-BASED CUTOFF PROPAGATION
# FIXED: COARSE/FINE CANDIDATE GRID EXPORT
# FIXED: HBFSS VOLCANO DECISION-BOUNDARY DISPLAY
# FIXED: VOLCANO LABEL FILTER TO EXCLUDE PLACEHOLDER SYMBOLS
# FIXED: APEGLM SHRINKAGE USES COEF, NOT CONTRAST
# FIXED: PCA NOW USES LOG2(X + 1) INPUT
# FIXED: pretty_dataset_label() MISSING FUNCTION
# FIXED: HBFSS BOUNDARY X-DOMAIN BUG
# =============================================================================
#
# WHAT THIS SCRIPT DOES
# 1. Reads a raw WTTS-Seq count matrix.
# 2. Builds comparison-specific EVS rankings from treatment and control PC1 loadings.
# 3. Uses a DESeq2-informed changepoint procedure on a combined ranking axis to choose
#    a leading-edge cutoff.
# 4. Splits each comparison into:
#      - original dataset
#      - leading-edge dataset
#      - remainder dataset
# 5. Runs DESeq2 and HBFSS on raw counts for each dataset.
# 6. Exports summary tables, result tables, diagnostics, and publication-style plots.
#
# IMPORTANT INTERPRETIVE NOTES
# - EVS ranking may be computed on normalized counts when normalize_before_evs = TRUE,
#   but final DESeq2/HBFSS inference is always performed on raw counts.
# - PCA for EVS is performed on log2(x + 1) transformed counts to reduce domination
#   by high-abundance features.
# - The HBFSS decision rule requires BOTH:
#      empirical_p < hc_p_threshold_dataset
#      HBFSS >= hbfss_threshold_dataset
# - Therefore, the HBFSS volcano boundary must display the upper envelope of:
#      y = -log10(hc_p_threshold_dataset)
#      y = hbfss_threshold_dataset / |lfc_shrunk|
#   rather than the hyperbola alone.
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
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    library(BiocParallel)
  }
})

# -----------------------------------------------------------------------------
# GLOBAL PATHS
# -----------------------------------------------------------------------------

repo_dir <- "/root/REAPER98632"
input_dir <- file.path(repo_dir, "data")
count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")

output_root <- file.path(repo_dir, "exports")
analysis_name <- "EVS_HBFSS_AllComparisons_Output"
output_dir <- file.path(output_root, analysis_name)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# GLOBAL CONTROLS
# -----------------------------------------------------------------------------

use_fast_apeglm <- TRUE
export_all_plots <- TRUE

detected_cores <- parallel::detectCores(logical = FALSE)
if (!is.finite(detected_cores) || is.na(detected_cores)) detected_cores <- 1L
n_workers <- max(1L, min(4L, as.integer(detected_cores) - 1L))

DESIGN_FORMULA <- ~ condition

alpha_level <- 0.10
lfc_boundary <- 1.0
max_usable_hc_p_threshold <- 0.95

figure_dpi <- 320
base_theme_size <- 10

# -----------------------------------------------------------------------------
# EVS / CUT-OFF CONTROLS
# -----------------------------------------------------------------------------

normalize_before_evs <- TRUE
evs_top_k_ranks_per_comparison <- 4L
evs_dispersion_source <- "dispGeneEst"

changepoint_min_segment_size <- 250L
changepoint_smoothing_window <- 101L

changepoint_model_selection_mode <- "bic"
changepoint_fallback_modes <- c("bic", "aic", "none")

require_regime_separation <- TRUE
min_regime_separation_delta <- 0.00

combined_loading_rule <- "max"

allow_rank_fallback_if_no_valid_changepoint <- FALSE
fallback_rank_top_n <- 5000L

use_coarse_to_fine_rank_search <- FALSE
coarse_rank_step <- 25L
fine_search_half_window <- 100L

auto_push_exports <- FALSE

# -----------------------------------------------------------------------------
# PLOT CONSTANTS
# -----------------------------------------------------------------------------

POINT_SIZE_PRIMARY  <- 1.6
POINT_SIZE_DISP     <- 1.3
POINT_ALPHA_PRIMARY <- 0.82
POINT_ALPHA_DISP    <- 0.55
POINT_STROKE        <- 0.40

LINE_WIDTH_BOUNDARY <- 0.55
LINE_WIDTH_ZERO     <- 0.40
LINE_WIDTH_THRESH   <- 0.90

plot_palette <- list(
  threshold = "#8C2D04",
  hbfss     = "#E67E22",
  deseq2    = "#C0392B",
  overlap   = "#7D3C98",
  weak      = "#4A90E2",
  control   = "#4D4D4D",
  treatment = "#1F78B4",
  iod       = "#2E86C1",
  cv2       = "#AF601A"
)

# -----------------------------------------------------------------------------
# PARALLEL HELPER
# -----------------------------------------------------------------------------

get_bpparam <- function(workers = n_workers) {
  workers <- as.integer(workers)[1]
  if (!is.finite(workers) || is.na(workers) || workers < 1L) workers <- 1L

  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    return(NULL)
  }

  if (.Platform$OS.type == "unix" && workers > 1L) {
    return(BiocParallel::MulticoreParam(workers = workers, progressbar = FALSE))
  }

  if (workers > 1L) {
    return(BiocParallel::SnowParam(workers = workers, progressbar = FALSE, type = "SOCK"))
  }

  BiocParallel::SerialParam(progressbar = FALSE)
}

# -----------------------------------------------------------------------------
# BASIC HELPERS
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

warn_if_noninteger_counts <- function(mat, label = "count matrix", tolerance = 1e-8) {
  x <- as.matrix(mat)
  bad <- is.finite(x) & !is.na(x) & abs(x - round(x)) > tolerance
  if (any(bad)) {
    warning(sprintf(
      "[%s] Non-integer values detected. Values will be rounded before DESeq2.",
      label
    ))
  }
}

safe_neglog10 <- function(x, pseudocount = 1e-12) {
  -log10(pmax(x, pseudocount))
}

compact_title <- function(x, width = 58) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

pretty_dataset_label <- function(x) {
  gsub("_", " ", as.character(x), fixed = TRUE)
}

clip_probabilities <- function(x, eps = 1e-300) {
  x <- unname(as.numeric(x))
  if (!length(x)) return(numeric(0))

  missing_idx <- is.na(x)
  neg_inf_idx <- is.infinite(x) & x < 0
  pos_inf_idx <- is.infinite(x) & x > 0
  finite_idx <- is.finite(x) & !missing_idx

  x[neg_inf_idx] <- eps
  x[pos_inf_idx] <- 1 - 1e-12
  x[finite_idx] <- pmin(pmax(x[finite_idx], eps), 1 - 1e-12)
  x[missing_idx] <- NA_real_
  x
}

fast_centered_rolling_mean <- function(x, k = changepoint_smoothing_window) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0L) return(numeric(0))
  if (!is.finite(k) || is.na(k) || k < 2L) return(x)

  k <- as.integer(k)
  if ((k %% 2L) == 0L) k <- k + 1L
  if (k <= 1L) return(x)

  if (k >= n) {
    k <- if ((n %% 2L) == 0L) n - 1L else n
    if ((k %% 2L) == 0L) k <- k - 1L
    if (k < 3L) return(x)
  }

  kernel <- rep(1 / k, k)
  left_pad  <- rep(x[1], k %/% 2L)
  right_pad <- rep(x[n], k %/% 2L)
  x_pad <- c(left_pad, x, right_pad)

  out <- stats::filter(x_pad, filter = kernel, sides = 2, method = "convolution")
  out <- as.numeric(out[(length(left_pad) + 1L):(length(left_pad) + n)])

  if (anyNA(out)) {
    bad <- is.na(out)
    out[bad] <- x[bad]
  }
  out
}

resolve_rank_or_fail <- function(rank_index, n_total, label = "cutoff") {
  n_total <- as.integer(n_total)[1]

  if (!is.finite(n_total) || is.na(n_total) || n_total < 2L) {
    stop(sprintf("[%s] Cannot resolve cutoff rank because n_total < 2 (n_total=%s).", label, n_total))
  }

  if (is.finite(rank_index) && !is.na(rank_index)) {
    return(max(1L, min(n_total - 1L, as.integer(rank_index))))
  }

  if (!isTRUE(allow_rank_fallback_if_no_valid_changepoint)) {
    stop(sprintf("[%s] No valid changepoint rank could be determined, and fallback is disabled.", label))
  }

  fallback_rank <- as.integer(fallback_rank_top_n)
  fallback_rank <- max(1L, min(n_total - 1L, fallback_rank))
  warning(sprintf("[%s] No valid changepoint rank found. Falling back to fixed top-N rank %d.", label, fallback_rank))
  fallback_rank
}

coalesce_numeric <- function(x, value = 0) {
  x <- as.numeric(x)
  x[is.na(x)] <- value
  x
}

resolve_model_selection_ladder <- function(primary_mode = changepoint_model_selection_mode,
                                           fallback_modes = changepoint_fallback_modes) {
  modes <- unique(c(primary_mode, fallback_modes))
  modes[modes %in% c("bic", "aic", "none")]
}

compute_combined_loading_score <- function(abs_trt, abs_ctrl, rule = combined_loading_rule) {
  abs_trt <- coalesce_numeric(abs_trt, 0)
  abs_ctrl <- coalesce_numeric(abs_ctrl, 0)

  if (identical(rule, "max")) return(pmax(abs_trt, abs_ctrl))
  if (identical(rule, "sum")) return(abs_trt + abs_ctrl)
  if (identical(rule, "l2"))  return(sqrt(abs_trt^2 + abs_ctrl^2))
  stop(sprintf("Unknown combined_loading_rule: %s", rule))
}

build_rank_grid <- function(shared_n_total,
                            min_segment_size = changepoint_min_segment_size,
                            use_coarse_to_fine = use_coarse_to_fine_rank_search,
                            coarse_step = coarse_rank_step) {
  start <- as.integer(min_segment_size)
  stop_i <- as.integer(shared_n_total - min_segment_size)
  if (stop_i < start) return(integer(0))

  if (!isTRUE(use_coarse_to_fine)) {
    return(seq.int(from = start, to = stop_i, by = 1L))
  }

  coarse_step <- max(1L, as.integer(coarse_step))
  grid <- unique(c(seq.int(from = start, to = stop_i, by = coarse_step), stop_i))
  as.integer(grid)
}

match_feature_to_rank <- function(feature_id, feature_vector) {
  if (!length(feature_id) || is.na(feature_id) || !nzchar(feature_id)) return(NA_integer_)
  idx <- match(as.character(feature_id)[1], as.character(feature_vector))
  if (!length(idx) || is.na(idx)) return(NA_integer_)
  as.integer(idx)
}

derive_group_cutoff_from_feature <- function(loading_table, selected_feature_id, value_col = "pc1_loading_abs") {
  idx <- match_feature_to_rank(selected_feature_id, loading_table$feature_id)
  if (!is.finite(idx) || is.na(idx)) {
    return(list(
      cutoff = NA_real_,
      top_n_used = NA_integer_,
      cutoff_quantile = NA_real_,
      split_class = rep("background_loading", nrow(loading_table))
    ))
  }

  split_class <- ifelse(seq_len(nrow(loading_table)) <= idx, "high_loading", "background_loading")
  cutoff_val <- as.numeric(loading_table[[value_col]][idx])
  cutoff_quantile <- 1 - (idx / nrow(loading_table))

  list(
    cutoff = cutoff_val,
    top_n_used = as.integer(idx),
    cutoff_quantile = cutoff_quantile,
    split_class = split_class
  )
}

resolve_selected_cutoff <- function(comparison_name, precomputed_selection) {
  cmp_obj <- precomputed_selection$per_comparison[[comparison_name]]
  if (is.null(cmp_obj)) {
    stop(sprintf("No precomputed comparison object found for %s", comparison_name))
  }

  cmp_eval <- cmp_obj$comparison_eval
  combined_tbl <- cmp_obj$combined_fit$loading_table
  n_total <- nrow(combined_tbl)

  selected_feature_id <- cmp_eval$best_feature_id
  selected_full_rank <- cmp_eval$best_rank

  feature_ok <- is.character(selected_feature_id) &&
    length(selected_feature_id) == 1L &&
    !is.na(selected_feature_id) &&
    nzchar(selected_feature_id) &&
    selected_feature_id %in% combined_tbl$feature_id

  rank_ok <- is.finite(selected_full_rank) &&
    !is.na(selected_full_rank) &&
    as.integer(selected_full_rank) >= 1L &&
    as.integer(selected_full_rank) <= n_total &&
    identical(
      as.character(combined_tbl$feature_id[as.integer(selected_full_rank)]),
      as.character(selected_feature_id)
    )

  if (!feature_ok || !rank_ok) {
    stop(sprintf(
      "[%s] Precomputed cutoff resolution failed because best_feature_id and best_rank were not internally consistent.",
      comparison_name
    ))
  }

  filtered_rank <- match_feature_to_rank(
    selected_feature_id,
    cmp_obj$combined_rank_nb_table$feature_id
  )

  list(
    cutoff_feature_id = as.character(selected_feature_id),
    full_rank = as.integer(selected_full_rank),
    filtered_rank = as.integer(filtered_rank),
    fallback_used = FALSE
  )
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

ranking_space_label <- function(normalize_before_evs_flag = normalize_before_evs) {
  if (isTRUE(normalize_before_evs_flag)) "Normalized before EVS" else "Raw counts before EVS"
}

final_analysis_label <- function() {
  "Final DESeq2 and HBFSS performed on raw counts"
}

reorder_result_columns <- function(final_df) {
  preferred_front_cols <- c(
    "site_num", "gene_name", "description", "Ensembl.ID", "feature_id",
    "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
    "log2FoldChange_shrunken", "empirical_p", "pi_valueE", "HBFSS",
    "effect_class", "standard_significant", "HBFSS_significant"
  )

  keep_front <- intersect(preferred_front_cols, names(final_df))
  keep_rest  <- setdiff(names(final_df), keep_front)
  final_df[, c(keep_front, keep_rest), drop = FALSE]
}

# -----------------------------------------------------------------------------
# THEME HELPERS
# -----------------------------------------------------------------------------

base_plot_theme <- function(base_size = base_theme_size) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.subtitle = element_text(hjust = 0, size = base_size),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 0.2, colour = "black"),
      legend.title = element_text(face = "bold", size = base_size - 0.1),
      legend.text = element_text(size = base_size - 0.2),
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(8, 16, 12, 10),
      legend.spacing.y = unit(0.12, "cm"),
      legend.box.spacing = unit(0.18, "cm")
    )
}

pca_theme <- function(base_size = base_theme_size) {
  base_plot_theme(base_size) +
    theme(
      legend.position = "right",
      legend.box = "vertical"
    )
}

volcano_theme <- function(base_size = base_theme_size) {
  base_plot_theme(base_size) +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.34, "cm"),
      legend.spacing.y = unit(0.10, "cm")
    )
}

dispersion_theme <- function(base_size = base_theme_size) {
  base_plot_theme(base_size) +
    theme(
      legend.position = "right",
      legend.box = "vertical"
    )
}

nb_theme <- function(base_size = base_theme_size) {
  base_plot_theme(base_size) +
    theme(
      legend.position = "right"
    )
}

# -----------------------------------------------------------------------------
# DESEQ2 / HBFSS HELPERS
# -----------------------------------------------------------------------------

safe_shrink_lfc <- function(dds, contrast_vector, coef_name = NULL, use_fast = use_fast_apeglm) {
  if (is.null(coef_name) || !length(coef_name) || is.na(coef_name) || !nzchar(coef_name)) {
    stop("safe_shrink_lfc() requires a valid DESeq2 coefficient name when type='apeglm'.")
  }

  if (isTRUE(use_fast)) {
    out <- tryCatch(
      lfcShrink(dds, coef = coef_name, type = "apeglm", apeMethod = "nbinomC"),
      error = function(e) e
    )
    if (!inherits(out, "error")) return(out)
    message("Fast apeglm shrinkage failed; retrying with default apeglm.")
  }

  lfcShrink(dds, coef = coef_name, type = "apeglm")
}

compute_hc_p_threshold_dataset <- function(stat_vec, max_threshold = max_usable_hc_p_threshold) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec)]
  if (!length(stat_vec)) return(NA_real_)

  fd <- tryCatch(
    fdrtool(stat_vec, statistic = "normal", plot = FALSE, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(fd)) return(NA_real_)

  p_local <- fd$param[["cutoff"]]
  if (!is.finite(p_local) || is.na(p_local)) return(NA_real_)

  p_local <- 2 * pnorm(-abs(p_local))
  if (!is.finite(p_local) || is.na(p_local) || p_local <= 0 || p_local >= max_threshold) {
    return(NA_real_)
  }
  p_local
}

compute_pi_valueE <- function(lfc_shrunk, empirical_p) {
  lfc_shrunk <- as.numeric(lfc_shrunk)
  empirical_p <- clip_probabilities(empirical_p)

  signed_component <- sign(lfc_shrunk) * pmax(abs(lfc_shrunk) - lfc_boundary, 0)
  signed_component * safe_neglog10(empirical_p)
}

classify_effect_strength <- function(abs_lfc) {
  ifelse(abs_lfc >= lfc_boundary, "strong_effect", "weak_effect")
}

# -----------------------------------------------------------------------------
# INPUT / ANNOTATION
# -----------------------------------------------------------------------------

load_count_matrix <- function(path) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"ID" %in% names(df)) {
    stop("Count matrix must contain an 'ID' column.")
  }

  feature_ids <- make.unique(as.character(df$ID))
  count_cols <- setdiff(names(df), "ID")
  mat <- as.matrix(df[, count_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- feature_ids
  mat
}

extract_annotation <- function(path) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  assert_required_columns(df, "ID", "raw annotation file")

  df$feature_id <- make.unique(as.character(df$ID))

  candidate_gene_cols <- c("gene_name", "Gene", "gene", "gene.symbol", "SYMBOL")
  gene_col <- candidate_gene_cols[candidate_gene_cols %in% names(df)][1]
  if (is.na(gene_col)) {
    df$gene_name <- df$feature_id
  } else {
    df$gene_name <- as.character(df[[gene_col]])
  }

  keep_cols <- intersect(
    c("feature_id", "gene_name", "site_num", "description", "Ensembl.ID"),
    names(df)
  )
  unique(df[, keep_cols, drop = FALSE])
}

# -----------------------------------------------------------------------------
# COMPARISON SETUP
# -----------------------------------------------------------------------------

build_comparison_inputs <- function(count_mat) {
  sample_names <- colnames(count_mat)

  comparison_specs <- list(
    RT0_ZT6  = list(trt = c("RT0_1", "RT0_2", "RT0_3"), untrt = c("ZT6_1", "ZT6_2", "ZT6_3")),
    RT2_ZT8  = list(trt = c("RT2_1", "RT2_2", "RT2_3"), untrt = c("ZT8_1", "ZT8_2", "ZT8_3")),
    RT4_ZT10 = list(trt = c("RT4_1", "RT4_2", "RT4_3"), untrt = c("ZT10_1", "ZT10_2", "ZT10_3")),
    RT8_ZT14 = list(trt = c("RT8_1", "RT8_2", "RT8_3"), untrt = c("ZT14_1", "ZT14_2", "ZT14_3"))
  )

  out <- list()

  for (nm in names(comparison_specs)) {
    spec <- comparison_specs[[nm]]
    missing_samples <- setdiff(c(spec$trt, spec$untrt), sample_names)
    if (length(missing_samples)) {
      stop(sprintf(
        "Comparison %s is missing required columns: %s",
        nm, paste(missing_samples, collapse = ", ")
      ))
    }

    cm <- count_mat[, c(spec$trt, spec$untrt), drop = FALSE]
    coldata <- data.frame(
      row.names = colnames(cm),
      condition = factor(
        c(rep("trt", length(spec$trt)), rep("untrt", length(spec$untrt))),
        levels = c("untrt", "trt")
      )
    )

    out[[nm]] <- list(
      count_matrix = cm,
      coldata = coldata
    )
  }

  out
}

# -----------------------------------------------------------------------------
# PCA HELPERS
# -----------------------------------------------------------------------------

compute_pca_from_counts <- function(count_matrix_subset, rank_mode_label = "counts") {
  mat <- as.matrix(count_matrix_subset)
  storage.mode(mat) <- "numeric"

  if (nrow(mat) < 2L) stop("PCA requires at least 2 features.")
  if (ncol(mat) < 2L) stop("PCA requires at least 2 samples.")

  mat <- log2(mat + 1)

  row_sds <- apply(mat, 1, sd, na.rm = TRUE)
  keep <- is.finite(row_sds) & !is.na(row_sds) & row_sds > 0
  if (sum(keep) < 2L) {
    stop("Not enough variable features remain after removing zero-variance rows for PCA.")
  }

  mat <- mat[keep, , drop = FALSE]

  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  variance_explained <- (pca$sdev^2) / sum(pca$sdev^2)

  list(
    pca = pca,
    feature_ids = rownames(mat),
    variance_explained = variance_explained,
    rank_mode_label = rank_mode_label
  )
}

compute_pc1_loading_table <- function(count_matrix_all, sample_ids, preprocessing_label = "Input used for PCA") {
  cm <- count_matrix_all[, sample_ids, drop = FALSE]
  pca_obj <- compute_pca_from_counts(cm, rank_mode_label = preprocessing_label)

  rot <- pca_obj$pca$rotation[, 1, drop = FALSE]
  loading_tbl <- data.frame(
    feature_id = rownames(rot),
    pc1_loading = as.numeric(rot[, 1]),
    pc1_loading_abs = abs(as.numeric(rot[, 1])),
    stringsAsFactors = FALSE
  )
  loading_tbl <- loading_tbl[order(-loading_tbl$pc1_loading_abs, loading_tbl$feature_id), , drop = FALSE]
  loading_tbl$rank <- seq_len(nrow(loading_tbl))

  scores_df <- as.data.frame(pca_obj$pca$x[, 1:2, drop = FALSE])
  scores_df$sample <- rownames(scores_df)
  scores_df$group <- sample_ids[match(scores_df$sample, sample_ids)]

  list(
    loading_table = loading_tbl,
    sample_scores = scores_df,
    variance_explained = pca_obj$variance_explained,
    preprocessing_label = preprocessing_label
  )
}

build_combined_pc1_loading_table <- function(fit_trt,
                                             fit_untrt,
                                             preprocessing_label,
                                             combined_rule = combined_loading_rule) {
  trt_tbl <- fit_trt$loading_table[, c("feature_id", "pc1_loading", "pc1_loading_abs", "rank")]
  names(trt_tbl) <- c("feature_id", "pc1_loading_trt", "pc1_loading_abs_trt", "rank_trt")

  untrt_tbl <- fit_untrt$loading_table[, c("feature_id", "pc1_loading", "pc1_loading_abs", "rank")]
  names(untrt_tbl) <- c("feature_id", "pc1_loading_untrt", "pc1_loading_abs_untrt", "rank_untrt")

  combined_tbl <- merge(trt_tbl, untrt_tbl, by = "feature_id", all = TRUE, sort = FALSE)

  combined_tbl$pc1_loading_abs_trt <- coalesce_numeric(combined_tbl$pc1_loading_abs_trt, 0)
  combined_tbl$pc1_loading_abs_untrt <- coalesce_numeric(combined_tbl$pc1_loading_abs_untrt, 0)

  combined_tbl$combined_pc1_loading_abs <- compute_combined_loading_score(
    abs_trt = combined_tbl$pc1_loading_abs_trt,
    abs_ctrl = combined_tbl$pc1_loading_abs_untrt,
    rule = combined_rule
  )

  combined_tbl <- combined_tbl[
    order(-combined_tbl$combined_pc1_loading_abs, combined_tbl$feature_id),
    ,
    drop = FALSE
  ]
  combined_tbl$rank <- seq_len(nrow(combined_tbl))
  combined_tbl$split_class <- "background_loading"

  list(
    loading_table = combined_tbl,
    preprocessing_label = preprocessing_label,
    combined_rule = combined_rule,
    cutoff = NA_real_,
    top_n_used = NA_integer_,
    cutoff_quantile = NA_real_,
    cutoff_method = "combined_comparison_specific_rank",
    selected_reason = NA_character_
  )
}

apply_combined_loading_cutoff <- function(combined_fit, rank_index, selected_reason) {
  loading_tbl <- combined_fit$loading_table
  n_total <- nrow(loading_tbl)
  rank_index <- resolve_rank_or_fail(rank_index, n_total, label = "combined_rank_cutoff")

  cutoff_value <- as.numeric(loading_tbl$combined_pc1_loading_abs[rank_index])
  loading_tbl$split_class <- ifelse(
    loading_tbl$rank <= rank_index,
    "high_loading",
    "background_loading"
  )

  combined_fit$loading_table <- loading_tbl
  combined_fit$cutoff <- cutoff_value
  combined_fit$top_n_used <- as.integer(rank_index)
  combined_fit$cutoff_quantile <- 1 - (rank_index / nrow(loading_tbl))
  combined_fit$selected_reason <- selected_reason
  combined_fit
}

# -----------------------------------------------------------------------------
# DISPERSION / NB REGIME HELPERS
# -----------------------------------------------------------------------------

extract_dispersion_column <- function(res_df, preferred = evs_dispersion_source) {
  candidates <- c(preferred, setdiff(c("dispGeneEst", "dispersion"), preferred))
  candidates <- unique(candidates)
  hit <- candidates[candidates %in% names(res_df)][1]
  if (is.na(hit)) {
    stop("Could not find a usable dispersion column in DESeq2 results.")
  }
  hit
}

fit_free_slope_segment <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 3L) {
    return(list(rss = Inf, slope = NA_real_, intercept = NA_real_))
  }

  fit <- lm(y ~ x)
  rss <- sum(resid(fit)^2)
  slope <- as.numeric(coef(fit)[["x"]])
  intercept <- as.numeric(coef(fit)[["(Intercept)"]])

  list(rss = rss, slope = slope, intercept = intercept)
}

compute_segmented_model_score <- function(split_rss,
                                          null_rss,
                                          n_obs,
                                          split_n_params = 5L,
                                          null_n_params = 4L,
                                          mode = changepoint_model_selection_mode) {
  split_rss <- as.numeric(split_rss)[1]
  null_rss  <- as.numeric(null_rss)[1]
  n_obs     <- as.integer(n_obs)[1]

  if (!is.finite(split_rss) || !is.finite(null_rss) || is.na(split_rss) || is.na(null_rss)) {
    return(NA_real_)
  }
  if (n_obs < 3L || split_rss <= 0 || null_rss <= 0) {
    return(NA_real_)
  }

  mode <- match.arg(mode, choices = c("bic", "aic", "none"))

  if (identical(mode, "none")) {
    return(log(null_rss / split_rss))
  }

  penalty_const <- if (identical(mode, "bic")) log(n_obs) else 2
  split_score <- n_obs * log(split_rss / n_obs) + penalty_const * split_n_params
  null_score  <- n_obs * log(null_rss  / n_obs) + penalty_const * null_n_params
  null_score - split_score
}

build_full_dataset_deseq2_dispersion_table <- function(count_matrix,
                                                       coldata,
                                                       dataset_label = "full_dataset_for_cutoff") {
  warn_if_noninteger_counts(count_matrix, label = dataset_label)

  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_matrix)),
    colData = coldata,
    design = DESIGN_FORMULA
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, fitType = "parametric")

  norm_counts <- counts(dds, normalized = TRUE)
  cond_levels <- as.character(coldata$condition)
  idx_trt <- which(cond_levels == "trt")
  idx_untrt <- which(cond_levels == "untrt")

  if (!length(idx_trt) || !length(idx_untrt)) {
    stop(sprintf("[%s] Both trt and untrt groups are required for full-dataset cutoff scoring.", dataset_label))
  }

  mean_trt <- rowMeans(norm_counts[, idx_trt, drop = FALSE], na.rm = TRUE)
  mean_untrt <- rowMeans(norm_counts[, idx_untrt, drop = FALSE], na.rm = TRUE)
  mean_avg <- rowMeans(cbind(mean_trt, mean_untrt), na.rm = TRUE)

  disp_gene <- mcols(dds)$dispGeneEst
  disp_shrunk <- dispersions(dds)

  disp_tbl <- data.frame(
    feature_id = rownames(dds),
    baseMean = as.numeric(mean_avg),
    mean_trt = as.numeric(mean_trt),
    mean_untrt = as.numeric(mean_untrt),
    dispGeneEst = as.numeric(disp_gene),
    dispersion = as.numeric(disp_shrunk),
    stringsAsFactors = FALSE
  )

  disp_col <- extract_dispersion_column(disp_tbl, preferred = evs_dispersion_source)
  alpha <- pmax(as.numeric(disp_tbl[[disp_col]]), 1e-12)
  mu <- pmax(as.numeric(disp_tbl$baseMean), 1e-12)

  disp_tbl$dispersion_used <- alpha
  disp_tbl$iod_nb <- 1 + alpha * mu
  disp_tbl$cv2_nb <- (1 / mu) + alpha
  disp_tbl$log_baseMean <- log10(mu)
  disp_tbl$log_iod_nb <- log10(pmax(disp_tbl$iod_nb, 1e-12))
  disp_tbl$log_cv2_nb <- log10(pmax(disp_tbl$cv2_nb, 1e-12))

  disp_tbl[order(disp_tbl$feature_id), , drop = FALSE]
}

build_ranked_nb_regime_table <- function(rank_order_ids, deseq2_disp_tbl) {
  rank_order_ids <- as.character(rank_order_ids)
  rank_map <- data.frame(
    feature_id = rank_order_ids,
    filtered_rank_index = seq_along(rank_order_ids),
    stringsAsFactors = FALSE
  )

  merged <- merge(rank_map, deseq2_disp_tbl, by = "feature_id", all.x = TRUE, sort = FALSE)
  merged <- merged[order(merged$filtered_rank_index), , drop = FALSE]

  merged$original_rank_index <- merged$filtered_rank_index
  merged$log_baseMean_smooth <- fast_centered_rolling_mean(merged$log_baseMean)
  merged$log_iod_nb_smooth   <- fast_centered_rolling_mean(merged$log_iod_nb)
  merged$log_cv2_nb_smooth   <- fast_centered_rolling_mean(merged$log_cv2_nb)

  merged
}

# -----------------------------------------------------------------------------
# EVS CUT-OFF SELECTION
# -----------------------------------------------------------------------------

fit_nb_regime_split_score <- function(rank_tbl,
                                      rank_index,
                                      min_segment_size = changepoint_min_segment_size,
                                      mode = changepoint_model_selection_mode,
                                      min_delta = min_regime_separation_delta,
                                      require_separation = require_regime_separation,
                                      wave_window = changepoint_smoothing_window,
                                      w_nb2 = 1.0,
                                      w_nb1 = 1.0,
                                      w_wave = 0.50,
                                      w_sep = 0.50) {
  n_total <- nrow(rank_tbl)

  empty_row <- function(rank_index_value = rank_index) {
    data.frame(
      filtered_rank_index = as.integer(rank_index_value)[1],
      cutoff_feature_id = NA_character_,
      original_rank_index = NA_integer_,
      objective_score = NA_real_,
      valid_split = FALSE,
      lead_rss = NA_real_,
      rem_rss = NA_real_,
      lead_slope = NA_real_,
      rem_slope = NA_real_,
      null_rss = NA_real_,
      improvement_vs_null = NA_real_,
      regime_plausible = FALSE,
      lead_iod_mean = NA_real_,
      rem_iod_mean = NA_real_,
      lead_cv2_mean = NA_real_,
      rem_cv2_mean = NA_real_,
      nb2_left_score = NA_real_,
      nb1_right_score = NA_real_,
      wave_score = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  if (!is.finite(rank_index) || is.na(rank_index) || n_total < (2L * min_segment_size + 1L)) {
    return(empty_row(rank_index))
  }

  rank_index <- as.integer(rank_index)[1]
  if (rank_index <= 1L || rank_index >= n_total) {
    return(empty_row(rank_index))
  }

  n_lead <- rank_index
  n_rem  <- n_total - rank_index
  if (n_lead < min_segment_size || n_rem < min_segment_size) {
    return(empty_row(rank_index))
  }

  lead_idx <- seq_len(n_lead)
  rem_idx  <- (n_lead + 1L):n_total
  lead_df <- rank_tbl[lead_idx, , drop = FALSE]
  rem_df  <- rank_tbl[rem_idx,  , drop = FALSE]

  fit_lead <- fit_free_slope_segment(lead_df$log_baseMean_smooth, lead_df$log_iod_nb_smooth)
  fit_rem  <- fit_free_slope_segment(rem_df$log_baseMean_smooth, rem_df$log_cv2_nb_smooth)
  split_rss <- fit_lead$rss + fit_rem$rss

  fit_null_iod <- fit_free_slope_segment(rank_tbl$log_baseMean_smooth, rank_tbl$log_iod_nb_smooth)
  fit_null_cv2 <- fit_free_slope_segment(rank_tbl$log_baseMean_smooth, rank_tbl$log_cv2_nb_smooth)
  null_rss <- fit_null_iod$rss + fit_null_cv2$rss

  base_score <- compute_segmented_model_score(
    split_rss = split_rss,
    null_rss  = null_rss,
    n_obs     = n_total,
    split_n_params = 5L,
    null_n_params  = 4L,
    mode = mode
  )

  lead_iod_mean <- mean(lead_df$log_iod_nb_smooth, na.rm = TRUE)
  rem_iod_mean  <- mean(rem_df$log_iod_nb_smooth,  na.rm = TRUE)
  lead_cv2_mean <- mean(lead_df$log_cv2_nb_smooth, na.rm = TRUE)
  rem_cv2_mean  <- mean(rem_df$log_cv2_nb_smooth,  na.rm = TRUE)

  nb2_left_score <- compute_segmented_model_score(
    split_rss = fit_lead$rss,
    null_rss  = fit_null_iod$rss,
    n_obs     = nrow(lead_df),
    split_n_params = 2L,
    null_n_params  = 2L,
    mode = "none"
  )

  nb1_right_score <- compute_segmented_model_score(
    split_rss = fit_rem$rss,
    null_rss  = fit_null_cv2$rss,
    n_obs     = nrow(rem_df),
    split_n_params = 2L,
    null_n_params  = 2L,
    mode = "none"
  )

  B <- rank_tbl$log_iod_nb_smooth
  C <- rank_tbl$log_cv2_nb_smooth
  dB <- c(NA_real_, diff(B))
  half_window <- max(5L, as.integer((wave_window - 1L) / 2L))
  left_lo  <- max(2L, rank_index - half_window)
  left_hi  <- max(left_lo, rank_index - 1L)
  right_lo <- min(n_total - 1L, rank_index + 1L)
  right_hi <- min(n_total - 1L, rank_index + half_window)

  if (left_hi <= left_lo || right_hi <= right_lo) {
    return(empty_row(rank_index))
  }

  iod_slope_jump <- mean(dB[right_lo:right_hi], na.rm = TRUE) - mean(dB[left_lo:left_hi], na.rm = TRUE)
  cv2_amp_left   <- sd(C[left_lo:left_hi],  na.rm = TRUE)
  cv2_amp_right  <- sd(C[right_lo:right_hi], na.rm = TRUE)
  cv2_amp_drop   <- cv2_amp_left - cv2_amp_right
  wave_score     <- iod_slope_jump + cv2_amp_drop

  regime_plausible <- is.finite(lead_iod_mean) &&
    is.finite(rem_iod_mean) &&
    is.finite(lead_cv2_mean) &&
    is.finite(rem_cv2_mean) &&
    is.finite(nb2_left_score) &&
    is.finite(nb1_right_score) &&
    ((lead_iod_mean - rem_iod_mean) > min_delta) &&
    ((rem_cv2_mean - lead_cv2_mean) > min_delta) &&
    (fit_lead$slope > 0) &&
    (fit_rem$slope < 0) &&
    (cv2_amp_drop > 0)

  if (!isTRUE(require_separation)) {
    regime_plausible <- TRUE
  }

  joint_score <-
    w_nb2 * nb2_left_score +
    w_nb1 * nb1_right_score +
    w_wave * wave_score +
    w_sep  * ((lead_iod_mean - rem_iod_mean) + (rem_cv2_mean - lead_cv2_mean))

  valid_split <- is.finite(base_score) &&
    is.finite(joint_score) &&
    (base_score > 0) &&
    (nb2_left_score > 0) &&
    (nb1_right_score > 0) &&
    regime_plausible

  data.frame(
    filtered_rank_index = rank_index,
    cutoff_feature_id = as.character(rank_tbl$feature_id[rank_index]),
    original_rank_index = as.integer(rank_tbl$original_rank_index[rank_index]),
    objective_score = as.numeric(joint_score)[1],
    valid_split = valid_split,
    lead_rss = fit_lead$rss,
    rem_rss = fit_rem$rss,
    lead_slope = fit_lead$slope,
    rem_slope = fit_rem$slope,
    null_rss = null_rss,
    improvement_vs_null = as.numeric(base_score)[1],
    regime_plausible = regime_plausible,
    lead_iod_mean = lead_iod_mean,
    rem_iod_mean = rem_iod_mean,
    lead_cv2_mean = lead_cv2_mean,
    rem_cv2_mean = rem_cv2_mean,
    nb2_left_score = as.numeric(nb2_left_score)[1],
    nb1_right_score = as.numeric(nb1_right_score)[1],
    wave_score = as.numeric(wave_score)[1],
    stringsAsFactors = FALSE
  )
}

evaluate_comparison_shared_rank_candidates <- function(raw_count_matrix,
                                                       comparison_name,
                                                       combined_fit,
                                                       coldata,
                                                       top_k = evs_top_k_ranks_per_comparison,
                                                       min_segment_size = changepoint_min_segment_size,
                                                       model_modes = resolve_model_selection_ladder()) {
  combined_rank_ids <- as.character(combined_fit$loading_table$feature_id[order(combined_fit$loading_table$rank)])

  full_disp_tbl <- build_full_dataset_deseq2_dispersion_table(
    count_matrix = raw_count_matrix,
    coldata = coldata,
    dataset_label = paste0(comparison_name, "_full_dataset_for_cutoff")
  )

  combined_nb_tbl <- build_ranked_nb_regime_table(
    rank_order_ids = combined_rank_ids,
    deseq2_disp_tbl = full_disp_tbl
  )

  shared_n_total <- nrow(combined_nb_tbl)

  if (!is.finite(shared_n_total) || is.na(shared_n_total) || shared_n_total <= (2L * min_segment_size)) {
    return(list(
      comparison_summary = data.frame(),
      candidate_grid = data.frame(),
      mode_diagnostics = data.frame(),
      best_rank = NA_integer_,
      best_score = NA_real_,
      best_mode = NA_character_,
      best_feature_id = NA_character_,
      best_filtered_rank = NA_integer_,
      top_ranks = integer(0),
      top_scores = numeric(0),
      combined_rank_nb_table = combined_nb_tbl,
      full_dispersion_table = full_disp_tbl
    ))
  }

  rank_grid <- build_rank_grid(
    shared_n_total = shared_n_total,
    min_segment_size = min_segment_size,
    use_coarse_to_fine = use_coarse_to_fine_rank_search,
    coarse_step = coarse_rank_step
  )

  all_mode_tables <- list()
  mode_diag_tables <- list()
  winning_mode <- NA_character_
  winning_valid_grid <- data.frame()

  for (mode_i in model_modes) {
    score_rows <- lapply(rank_grid, function(r) {
      s <- fit_nb_regime_split_score(
        rank_tbl = combined_nb_tbl,
        rank_index = r,
        min_segment_size = min_segment_size,
        mode = mode_i,
        min_delta = min_regime_separation_delta,
        require_separation = require_regime_separation
      )

      data.frame(
        comparison_name = comparison_name,
        model_selection_mode = mode_i,
        search_stage = "coarse",
        filtered_rank_index = s$filtered_rank_index[1],
        cutoff_feature_id = s$cutoff_feature_id[1],
        original_rank_index = s$original_rank_index[1],
        objective_score = s$objective_score[1],
        valid_split = s$valid_split[1],
        improvement_vs_null = s$improvement_vs_null[1],
        regime_plausible = s$regime_plausible[1],
        lead_iod_mean = s$lead_iod_mean[1],
        rem_iod_mean = s$rem_iod_mean[1],
        lead_cv2_mean = s$lead_cv2_mean[1],
        rem_cv2_mean = s$rem_cv2_mean[1],
        lead_slope = s$lead_slope[1],
        rem_slope = s$rem_slope[1],
        stringsAsFactors = FALSE
      )
    })

    coarse_grid <- dplyr::bind_rows(score_rows)
    coarse_grid <- coarse_grid[order(-coarse_grid$objective_score, coarse_grid$filtered_rank_index), , drop = FALSE]

    export_grid <- coarse_grid

    mode_diag_tables[[paste0(mode_i, "_coarse")]] <- data.frame(
      comparison_name = comparison_name,
      model_selection_mode = paste0(mode_i, "_coarse"),
      n_candidates = nrow(coarse_grid),
      n_positive_score = sum(is.finite(coarse_grid$objective_score) & coarse_grid$objective_score > 0, na.rm = TRUE),
      n_regime_plausible = sum(coarse_grid$regime_plausible, na.rm = TRUE),
      n_valid_split = sum(coarse_grid$valid_split, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    valid_grid <- coarse_grid[
      coarse_grid$valid_split &
        is.finite(coarse_grid$objective_score) &
        !is.na(coarse_grid$cutoff_feature_id),
      ,
      drop = FALSE
    ]

    if (nrow(valid_grid) > 0 && isTRUE(use_coarse_to_fine_rank_search) && fine_search_half_window > 0L) {
      coarse_best <- as.integer(valid_grid$filtered_rank_index[1])
      lo <- max(min_segment_size, coarse_best - as.integer(fine_search_half_window))
      hi <- min(shared_n_total - min_segment_size, coarse_best + as.integer(fine_search_half_window))
      fine_grid_ranks <- seq.int(lo, hi, by = 1L)

      fine_rows <- lapply(fine_grid_ranks, function(r) {
        s <- fit_nb_regime_split_score(
          rank_tbl = combined_nb_tbl,
          rank_index = r,
          min_segment_size = min_segment_size,
          mode = mode_i,
          min_delta = min_regime_separation_delta,
          require_separation = require_regime_separation
        )

        data.frame(
          comparison_name = comparison_name,
          model_selection_mode = mode_i,
          search_stage = "fine",
          filtered_rank_index = s$filtered_rank_index[1],
          cutoff_feature_id = s$cutoff_feature_id[1],
          original_rank_index = s$original_rank_index[1],
          objective_score = s$objective_score[1],
          valid_split = s$valid_split[1],
          improvement_vs_null = s$improvement_vs_null[1],
          regime_plausible = s$regime_plausible[1],
          lead_iod_mean = s$lead_iod_mean[1],
          rem_iod_mean = s$rem_iod_mean[1],
          lead_cv2_mean = s$lead_cv2_mean[1],
          rem_cv2_mean = s$rem_cv2_mean[1],
          lead_slope = s$lead_slope[1],
          rem_slope = s$rem_slope[1],
          stringsAsFactors = FALSE
        )
      })

      fine_grid_df <- dplyr::bind_rows(fine_rows)
      fine_grid_df <- fine_grid_df[order(-fine_grid_df$objective_score, fine_grid_df$filtered_rank_index), , drop = FALSE]

      export_grid <- dplyr::bind_rows(coarse_grid, fine_grid_df)

      mode_diag_tables[[paste0(mode_i, "_fine")]] <- data.frame(
        comparison_name = comparison_name,
        model_selection_mode = paste0(mode_i, "_fine"),
        n_candidates = nrow(fine_grid_df),
        n_positive_score = sum(is.finite(fine_grid_df$objective_score) & fine_grid_df$objective_score > 0, na.rm = TRUE),
        n_regime_plausible = sum(fine_grid_df$regime_plausible, na.rm = TRUE),
        n_valid_split = sum(fine_grid_df$valid_split, na.rm = TRUE),
        stringsAsFactors = FALSE
    )

      valid_grid <- fine_grid_df[
        fine_grid_df$valid_split &
          is.finite(fine_grid_df$objective_score) &
          !is.na(fine_grid_df$cutoff_feature_id),
        ,
        drop = FALSE
      ]
    }

    all_mode_tables[[mode_i]] <- export_grid

    if (nrow(valid_grid) > 0) {
      valid_grid <- valid_grid[order(-valid_grid$objective_score, valid_grid$filtered_rank_index), , drop = FALSE]
      if (!nrow(winning_valid_grid) || valid_grid$objective_score[1] > winning_valid_grid$objective_score[1]) {
        winning_mode <- mode_i
        winning_valid_grid <- valid_grid
      }
    }
  }

  candidate_grid <- dplyr::bind_rows(all_mode_tables)
  mode_diagnostics <- dplyr::bind_rows(mode_diag_tables)

  if (!nrow(winning_valid_grid)) {
    return(list(
      comparison_summary = data.frame(
        comparison_name = comparison_name,
        selected_rank_order = 1L,
        independently_best_rank = NA_integer_,
        independently_best_local_lagrange_objective = NA_real_,
        model_selection_mode_used = ifelse(length(model_modes), tail(model_modes, 1), NA_character_),
        filtered_rank_index = NA_integer_,
        regime_plausible = NA,
        cutoff_feature_id = NA_character_,
        stringsAsFactors = FALSE
      ),
      candidate_grid = candidate_grid,
      mode_diagnostics = mode_diagnostics,
      best_rank = NA_integer_,
      best_score = NA_real_,
      best_mode = "fallback_only",
      best_feature_id = NA_character_,
      best_filtered_rank = NA_integer_,
      top_ranks = integer(0),
      top_scores = numeric(0),
      combined_rank_nb_table = combined_nb_tbl,
      full_dispersion_table = full_disp_tbl
    ))
  }

  top_k <- min(as.integer(top_k), nrow(winning_valid_grid))
  top_rows <- winning_valid_grid[seq_len(top_k), , drop = FALSE]
  top_rows$mapped_full_rank <- match(top_rows$cutoff_feature_id, combined_fit$loading_table$feature_id)

  comparison_summary <- data.frame(
    comparison_name = comparison_name,
    selected_rank_order = seq_len(nrow(top_rows)),
    independently_best_rank = as.integer(top_rows$mapped_full_rank),
    independently_best_local_lagrange_objective = as.numeric(top_rows$objective_score),
    model_selection_mode_used = winning_mode,
    filtered_rank_index = as.integer(top_rows$filtered_rank_index),
    regime_plausible = as.logical(top_rows$regime_plausible),
    cutoff_feature_id = as.character(top_rows$cutoff_feature_id),
    stringsAsFactors = FALSE
  )

  list(
    comparison_summary = comparison_summary,
    candidate_grid = candidate_grid,
    mode_diagnostics = mode_diagnostics,
    best_rank = as.integer(top_rows$mapped_full_rank[1]),
    best_score = as.numeric(top_rows$objective_score[1]),
    best_mode = winning_mode,
    best_feature_id = as.character(top_rows$cutoff_feature_id[1]),
    best_filtered_rank = as.integer(top_rows$filtered_rank_index[1]),
    top_ranks = as.integer(top_rows$mapped_full_rank),
    top_scores = as.numeric(top_rows$objective_score),
    combined_rank_nb_table = combined_nb_tbl,
    full_dispersion_table = full_disp_tbl
  )
}

# =============================================================================
# GLOBAL PRECOMPUTE
# =============================================================================

precompute_global_evs_rank_selection <- function(comparison_inputs,
                                                 top_k_ranks_per_comparison = evs_top_k_ranks_per_comparison) {
  per_comparison <- list()
  summary_rows <- list()
  mode_diag_rows <- list()

  for (cmp in names(comparison_inputs)) {
    input_obj <- comparison_inputs[[cmp]]

    warn_if_noninteger_counts(input_obj$count_matrix, label = paste0(cmp, " precompute count matrix"))

    dds_init <- DESeqDataSetFromMatrix(
      countData = round(as.matrix(input_obj$count_matrix)),
      colData = input_obj$coldata,
      design = DESIGN_FORMULA
    )
    dds_init <- dds_init[rowSums(counts(dds_init)) > 0, ]
    dds_init <- estimateSizeFactors(dds_init)

    retained_feature_ids <- rownames(dds_init)
    norm_counts_init <- as.data.frame(counts(dds_init, normalized = TRUE))
    raw_counts_init  <- as.data.frame(input_obj$count_matrix[retained_feature_ids, , drop = FALSE])

    sample_ids <- colnames(input_obj$count_matrix)
    trt_ids   <- sample_ids[input_obj$coldata$condition == "trt"]
    untrt_ids <- sample_ids[input_obj$coldata$condition == "untrt"]

    if (normalize_before_evs) {
      counts_for_evs <- norm_counts_init
      evs_preproc_label <- "Normalized before EVS"
    } else {
      counts_for_evs <- raw_counts_init
      evs_preproc_label <- "Raw counts before EVS"
    }

    fit_trt <- compute_pc1_loading_table(counts_for_evs, trt_ids, preprocessing_label = evs_preproc_label)
    fit_untrt <- compute_pc1_loading_table(counts_for_evs, untrt_ids, preprocessing_label = evs_preproc_label)

    fit_trt_raw_view <- compute_pc1_loading_table(raw_counts_init, trt_ids, preprocessing_label = "Raw counts before EVS")
    fit_untrt_raw_view <- compute_pc1_loading_table(raw_counts_init, untrt_ids, preprocessing_label = "Raw counts before EVS")

    combined_fit <- build_combined_pc1_loading_table(
      fit_trt = fit_trt,
      fit_untrt = fit_untrt,
      preprocessing_label = evs_preproc_label,
      combined_rule = combined_loading_rule
    )

    comparison_eval <- evaluate_comparison_shared_rank_candidates(
      raw_count_matrix = raw_counts_init,
      comparison_name = cmp,
      combined_fit = combined_fit,
      coldata = input_obj$coldata,
      top_k = top_k_ranks_per_comparison,
      model_modes = resolve_model_selection_ladder()
    )

    per_comparison[[cmp]] <- list(
      comparison_name = cmp,
      retained_feature_ids = retained_feature_ids,
      norm_counts_init = norm_counts_init,
      raw_counts_init = raw_counts_init,
      counts_for_evs = counts_for_evs,
      fit_trt = fit_trt,
      fit_untrt = fit_untrt,
      fit_trt_raw_view = fit_trt_raw_view,
      fit_untrt_raw_view = fit_untrt_raw_view,
      combined_fit = combined_fit,
      comparison_eval = comparison_eval,
      combined_rank_nb_table = comparison_eval$combined_rank_nb_table,
      full_dispersion_table = comparison_eval$full_dispersion_table
    )

    summary_rows[[cmp]] <- comparison_eval$comparison_summary
    mode_diag_rows[[cmp]] <- comparison_eval$mode_diagnostics
  }

  comparison_summary_table <- dplyr::bind_rows(summary_rows)
  mode_diagnostics_table <- dplyr::bind_rows(mode_diag_rows)

  list(
    comparison_summary_table = comparison_summary_table,
    mode_diagnostics_table = mode_diagnostics_table,
    per_comparison = per_comparison
  )
}

# =============================================================================
# FINALIZED EVS SPLIT
# =============================================================================

build_eigenvector_split <- function(count_matrix,
                                    coldata,
                                    comparison_name,
                                    selection_obj = NULL,
                                    precomputed_selection = NULL) {
  precomputed_cmp <- NULL
  if (!is.null(precomputed_selection$per_comparison) &&
      !is.null(precomputed_selection$per_comparison[[comparison_name]])) {
    precomputed_cmp <- precomputed_selection$per_comparison[[comparison_name]]
  }

  if (is.null(precomputed_cmp)) {
    stop("build_eigenvector_split() requires precomputed_selection for the finalized workflow.")
  }

  fit_trt <- precomputed_cmp$fit_trt
  fit_untrt <- precomputed_cmp$fit_untrt
  fit_trt_raw_view <- precomputed_cmp$fit_trt_raw_view
  fit_untrt_raw_view <- precomputed_cmp$fit_untrt_raw_view
  combined_fit <- precomputed_cmp$combined_fit
  comparison_eval <- precomputed_cmp$comparison_eval

  if (nrow(combined_fit$loading_table) < 2L) {
    stop(sprintf("[%s] Combined loading table has fewer than 2 features; cannot define a cutoff.", comparison_name))
  }

  if (is.null(selection_obj)) {
    selection_obj <- resolve_selected_cutoff(comparison_name, precomputed_selection)
  }

  selected_feature_id <- as.character(selection_obj$cutoff_feature_id)[1]
  final_full_rank <- as.integer(selection_obj$full_rank)[1]
  final_filtered_rank <- as.integer(selection_obj$filtered_rank)[1]

  if (!is.finite(final_full_rank) || is.na(final_full_rank)) {
    stop(sprintf("[%s] Failed to resolve a valid full-rank cutoff.", comparison_name))
  }

  final_reason <- if (isTRUE(selection_obj$fallback_used)) {
    "fallback_rank_used_after_failed_cutoff_resolution"
  } else if (normalize_before_evs) {
    "joint_nb2_nb1_wave_cutoff_selected_from_full_raw_dataset_after_normalized_EVS_combined_ranking"
  } else {
    "joint_nb2_nb1_wave_cutoff_selected_from_full_raw_dataset_after_raw_EVS_combined_ranking"
  }

  combined_fit <- apply_combined_loading_cutoff(combined_fit, final_full_rank, final_reason)

  high_ids <- as.character(subset(combined_fit$loading_table, split_class == "high_loading")$feature_id)

  analyzed_feature_ids <- intersect(
    as.character(combined_fit$loading_table$feature_id),
    rownames(count_matrix)
  )

  dropped_from_combined <- setdiff(as.character(combined_fit$loading_table$feature_id), rownames(count_matrix))
  if (length(dropped_from_combined) > 0L) {
    warning(sprintf(
      "[%s] %d features were present in the combined loading table but absent from count_matrix and were dropped before final dataset splitting.",
      comparison_name, length(dropped_from_combined)
    ))
  }

  leading_edge_ids <- intersect(high_ids, rownames(count_matrix))
  remainder_ids <- setdiff(analyzed_feature_ids, leading_edge_ids)

  if (!length(leading_edge_ids)) {
    stop("Leading-edge dataset is empty after final cutoff.")
  }

  trt_proj <- derive_group_cutoff_from_feature(
    loading_table = fit_trt$loading_table,
    selected_feature_id = selected_feature_id,
    value_col = "pc1_loading_abs"
  )
  fit_trt$cutoff <- trt_proj$cutoff
  fit_trt$top_n_used <- trt_proj$top_n_used
  fit_trt$cutoff_quantile <- trt_proj$cutoff_quantile
  fit_trt$selected_reason <- "feature_membership_projected_from_combined_cutoff"
  fit_trt$loading_table$split_class <- trt_proj$split_class

  untrt_proj <- derive_group_cutoff_from_feature(
    loading_table = fit_untrt$loading_table,
    selected_feature_id = selected_feature_id,
    value_col = "pc1_loading_abs"
  )
  fit_untrt$cutoff <- untrt_proj$cutoff
  fit_untrt$top_n_used <- untrt_proj$top_n_used
  fit_untrt$cutoff_quantile <- untrt_proj$cutoff_quantile
  fit_untrt$selected_reason <- "feature_membership_projected_from_combined_cutoff"
  fit_untrt$loading_table$split_class <- untrt_proj$split_class

  trt_raw_proj <- derive_group_cutoff_from_feature(
    loading_table = fit_trt_raw_view$loading_table,
    selected_feature_id = selected_feature_id,
    value_col = "pc1_loading_abs"
  )
  fit_trt_raw_view$cutoff <- trt_raw_proj$cutoff
  fit_trt_raw_view$top_n_used <- trt_raw_proj$top_n_used
  fit_trt_raw_view$cutoff_quantile <- trt_raw_proj$cutoff_quantile
  fit_trt_raw_view$selected_reason <- "feature_membership_projected_from_combined_cutoff"
  fit_trt_raw_view$loading_table$split_class <- trt_raw_proj$split_class

  untrt_raw_proj <- derive_group_cutoff_from_feature(
    loading_table = fit_untrt_raw_view$loading_table,
    selected_feature_id = selected_feature_id,
    value_col = "pc1_loading_abs"
  )
  fit_untrt_raw_view$cutoff <- untrt_raw_proj$cutoff
  fit_untrt_raw_view$top_n_used <- untrt_raw_proj$top_n_used
  fit_untrt_raw_view$cutoff_quantile <- untrt_raw_proj$cutoff_quantile
  fit_untrt_raw_view$selected_reason <- "feature_membership_projected_from_combined_cutoff"
  fit_untrt_raw_view$loading_table$split_class <- untrt_raw_proj$split_class

  final_rank_aggregation <- data.frame(
    comparison_name = comparison_name,
    independently_best_rank = comparison_eval$best_rank,
    independently_best_local_lagrange_objective = comparison_eval$best_score,
    independently_best_feature_id = comparison_eval$best_feature_id,
    model_selection_mode_used = comparison_eval$best_mode,
    final_selected_full_rank = final_full_rank,
    final_selected_filtered_rank = final_filtered_rank,
    final_selected_feature_id = selected_feature_id,
    fallback_used = isTRUE(selection_obj$fallback_used),
    stringsAsFactors = FALSE
  )

  list(
    fit_trt = fit_trt,
    fit_untrt = fit_untrt,
    fit_trt_raw_view = fit_trt_raw_view,
    fit_untrt_raw_view = fit_untrt_raw_view,
    combined_fit = combined_fit,
    leading_edge_ids = leading_edge_ids,
    remainder_ids = remainder_ids,
    selected_feature_id = selected_feature_id,
    final_full_rank = final_full_rank,
    final_filtered_rank = final_filtered_rank,
    final_rank_aggregation = final_rank_aggregation,
    comparison_eval = comparison_eval
  )
}

# -----------------------------------------------------------------------------
# DIAGNOSTIC EXPORTS
# -----------------------------------------------------------------------------

build_evs_candidate_export <- function(precomputed_selection,
                                       out_dir = file.path(output_dir, "tables")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  comparison_summary <- precomputed_selection$comparison_summary_table
  if (nrow(comparison_summary)) {
    save_csv(
      comparison_summary,
      file.path(out_dir, "EVS_candidate_cutoffs_by_comparison.csv")
    )
  }

  mode_diagnostics <- precomputed_selection$mode_diagnostics_table
  if (nrow(mode_diagnostics)) {
    save_csv(
      mode_diagnostics,
      file.path(out_dir, "EVS_model_mode_diagnostics.csv")
    )
  }

  for (cmp in names(precomputed_selection$per_comparison)) {
    cmp_obj <- precomputed_selection$per_comparison[[cmp]]
    cmp_dir <- file.path(out_dir, cmp)
    dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

    if (nrow(cmp_obj$comparison_eval$candidate_grid)) {
      save_csv(
        cmp_obj$comparison_eval$candidate_grid,
        file.path(cmp_dir, paste0(cmp, "_EVS_candidate_grid.csv"))
      )
    }

    if (nrow(cmp_obj$comparison_eval$comparison_summary)) {
      save_csv(
        cmp_obj$comparison_eval$comparison_summary,
        file.path(cmp_dir, paste0(cmp, "_EVS_independent_rank_aggregation_summary.csv"))
      )
    }

    if (nrow(cmp_obj$comparison_eval$mode_diagnostics)) {
      save_csv(
        cmp_obj$comparison_eval$mode_diagnostics,
        file.path(cmp_dir, paste0(cmp, "_EVS_mode_diagnostics.csv"))
      )
    }
  }

  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# PLOT HELPERS FOR PCA / EVS
# -----------------------------------------------------------------------------

variance_df_from_fit <- function(fit_obj, group_label, comparison_name) {
  ve <- fit_obj$variance_explained
  if (length(ve) < 2L) ve <- c(ve, rep(NA_real_, 2L - length(ve)))

  data.frame(
    comparison_name = comparison_name,
    group = group_label,
    PC = factor(c("PC1", "PC2"), levels = c("PC1", "PC2")),
    variance_explained = ve[1:2],
    preprocessing_label = fit_obj$preprocessing_label,
    stringsAsFactors = FALSE
  )
}

build_pca_variance_plot <- function(fit_trt, fit_untrt, comparison_name) {
  df <- dplyr::bind_rows(
    variance_df_from_fit(fit_trt, "Treatment", comparison_name),
    variance_df_from_fit(fit_untrt, "Control", comparison_name)
  )

  ggplot(df, aes(x = PC, y = variance_explained, fill = group)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.64, alpha = 0.92) +
    scale_fill_manual(values = c(
      Treatment = plot_palette$treatment,
      Control   = plot_palette$control
    )) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(
      title = compact_title(paste(comparison_name, "PCA variance explained")),
      subtitle = paste(
        ranking_space_label(),
        "| Separate PCA performed in treatment and control groups"
      ),
      x = NULL,
      y = "Explained variance",
      fill = NULL
    ) +
    pca_theme()
}

build_pca_scatter_plot <- function(fit_trt, fit_untrt, comparison_name) {
  score_df <- dplyr::bind_rows(
    transform(fit_trt$sample_scores, biological_group = "Treatment"),
    transform(fit_untrt$sample_scores, biological_group = "Control")
  )

  names(score_df)[names(score_df) == "PC1"] <- "PC1_score"
  names(score_df)[names(score_df) == "PC2"] <- "PC2_score"

  ggplot(score_df, aes(x = PC1_score, y = PC2_score, colour = biological_group, label = sample)) +
    geom_hline(yintercept = 0, linewidth = LINE_WIDTH_ZERO, linetype = "dotted", colour = "grey45") +
    geom_vline(xintercept = 0, linewidth = LINE_WIDTH_ZERO, linetype = "dotted", colour = "grey45") +
    geom_point(size = 3.0, alpha = 0.88) +
    geom_text(vjust = -0.72, size = 3.2, show.legend = FALSE) +
    scale_colour_manual(values = c(
      Treatment = plot_palette$treatment,
      Control   = plot_palette$control
    )) +
    labs(
      title = compact_title(paste(comparison_name, "PCA sample scores")),
      subtitle = paste(
        ranking_space_label(),
        "| PCA scores shown separately for treatment and control"
      ),
      x = "PC1 score",
      y = "PC2 score",
      colour = NULL
    ) +
    pca_theme()
}

build_pc1_loading_rank_plot <- function(loading_tbl,
                                        comparison_name,
                                        group_label,
                                        value_col = "pc1_loading_abs",
                                        cutoff_rank = NA_integer_) {
  df <- loading_tbl
  yval <- df[[value_col]]

  ggplot(df, aes(x = rank, y = yval)) +
    geom_line(linewidth = 0.45, colour = if (group_label == "Treatment") plot_palette$treatment else plot_palette$control) +
    {
      if (is.finite(cutoff_rank) && !is.na(cutoff_rank)) {
        geom_vline(
          xintercept = cutoff_rank,
          linetype = "dashed",
          linewidth = LINE_WIDTH_THRESH,
          colour = plot_palette$threshold
        )
      }
    } +
    labs(
      title = compact_title(paste(comparison_name, group_label, "PC1 loading rank profile")),
      subtitle = paste(ranking_space_label(), "| Ranked by absolute PC1 loading"),
      x = "Rank",
      y = "Absolute PC1 loading"
    ) +
    pca_theme()
}

build_pc1_loading_hist_plot <- function(loading_tbl, comparison_name, group_label, value_col = "pc1_loading_abs") {
  df <- loading_tbl

  ggplot(df, aes(x = .data[[value_col]])) +
    geom_histogram(
      bins = 60,
      fill = if (group_label == "Treatment") plot_palette$treatment else plot_palette$control,
      colour = "white",
      linewidth = 0.15,
      alpha = 0.92
    ) +
    labs(
      title = compact_title(paste(comparison_name, group_label, "PC1 loading histogram")),
      subtitle = paste(ranking_space_label(), "| Distribution of absolute PC1 loadings"),
      x = "Absolute PC1 loading",
      y = "Number of features"
    ) +
    pca_theme()
}

build_combined_loading_rank_plot <- function(combined_fit, comparison_name) {
  df <- combined_fit$loading_table

  ggplot(df, aes(x = rank, y = combined_pc1_loading_abs)) +
    geom_line(linewidth = 0.52, colour = plot_palette$threshold) +
    {
      if (is.finite(combined_fit$top_n_used) && !is.na(combined_fit$top_n_used)) {
        geom_vline(
          xintercept = combined_fit$top_n_used,
          linetype = "dashed",
          linewidth = LINE_WIDTH_THRESH,
          colour = plot_palette$cv2
        )
      }
    } +
    labs(
      title = compact_title(paste(comparison_name, "Combined EVS loading rank profile")),
      subtitle = paste(
        ranking_space_label(),
        "| Combined treatment/control ranking rule:",
        combined_fit$combined_rule
      ),
      x = "Combined rank",
      y = "Combined absolute PC1 loading"
    ) +
    pca_theme()
}

build_nb_rank_profile_plot <- function(combined_nb_tbl, comparison_name, cutoff_rank = NA_integer_) {
  long_df <- dplyr::bind_rows(
    data.frame(
      rank = combined_nb_tbl$filtered_rank_index,
      value = combined_nb_tbl$log_iod_nb_smooth,
      metric = "Leading-edge IOD",
      stringsAsFactors = FALSE
    ),
    data.frame(
      rank = combined_nb_tbl$filtered_rank_index,
      value = combined_nb_tbl$log_cv2_nb_smooth,
      metric = "Remainder CV²",
      stringsAsFactors = FALSE
    )
  )

  ggplot(long_df, aes(x = rank, y = value, colour = metric)) +
    geom_line(linewidth = 0.55, alpha = 0.94) +
    {
      if (is.finite(cutoff_rank) && !is.na(cutoff_rank)) {
        geom_vline(
          xintercept = cutoff_rank,
          linetype = "dashed",
          linewidth = LINE_WIDTH_THRESH,
          colour = plot_palette$threshold
        )
      }
    } +
    scale_colour_manual(values = c(
      "Leading-edge IOD" = plot_palette$iod,
      "Remainder CV²" = plot_palette$cv2
    )) +
    labs(
      title = compact_title(paste(comparison_name, "NB regime rank profiles on final combined rank")),
      subtitle = paste(
        "Full-dataset DESeq2 dispersion structure mapped onto the filtered combined EVS rank order"
      ),
      x = "Filtered combined rank",
      y = "Smoothed log-regime metric",
      colour = NULL
    ) +
    nb_theme()
}

build_nb_segment_fit_plot <- function(combined_nb_tbl, comparison_name, cutoff_rank = NA_integer_) {
  n_total <- nrow(combined_nb_tbl)
  if (!is.finite(cutoff_rank) || is.na(cutoff_rank)) {
    stop(sprintf("[%s] cutoff_rank is NA in build_nb_segment_fit_plot().", comparison_name))
  }
  cutoff_rank <- max(1L, min(n_total - 1L, as.integer(cutoff_rank)))

  plot_df <- dplyr::bind_rows(
    data.frame(
      x = combined_nb_tbl$log_baseMean_smooth[seq_len(cutoff_rank)],
      y = combined_nb_tbl$log_iod_nb_smooth[seq_len(cutoff_rank)],
      regime = "Leading-edge IOD",
      stringsAsFactors = FALSE
    ),
    data.frame(
      x = combined_nb_tbl$log_baseMean_smooth[(cutoff_rank + 1L):n_total],
      y = combined_nb_tbl$log_cv2_nb_smooth[(cutoff_rank + 1L):n_total],
      regime = "Remainder CV²",
      stringsAsFactors = FALSE
    )
  )

  ggplot(plot_df, aes(x = x, y = y, colour = regime)) +
    geom_point(size = POINT_SIZE_DISP, alpha = POINT_ALPHA_DISP) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.80) +
    scale_colour_manual(values = c(
      "Leading-edge IOD" = plot_palette$iod,
      "Remainder CV²" = plot_palette$cv2
    )) +
    labs(
      title = compact_title(paste(comparison_name, "NB regime segment-fit diagnostic")),
      subtitle = paste(
        "Leading-edge through filtered rank", cutoff_rank,
        "scored on IOD; remainder scored on CV²"
      ),
      x = "Smoothed log(baseMean)",
      y = "Smoothed log-regime metric",
      colour = NULL
    ) +
    nb_theme()
}

# -----------------------------------------------------------------------------
# DESEQ2 ANALYSIS
# -----------------------------------------------------------------------------

run_deseq2_analysis <- function(count_matrix_subset,
                                coldata_subset,
                                dataset_name,
                                annotation_df) {
  warn_if_noninteger_counts(count_matrix_subset, label = dataset_name)

  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_matrix_subset)),
    colData = coldata_subset,
    design = DESIGN_FORMULA
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]

  if (nrow(dds) < 2L) {
    stop(sprintf("[%s] Fewer than 2 features remain after zero-row filtering.", dataset_name))
  }

  bp <- get_bpparam()
  dds <- DESeq(dds, parallel = !is.null(bp), BPPARAM = bp)

  res <- results(dds, contrast = c("condition", "trt", "untrt"), alpha = alpha_level)

  coef_name <- resultsNames(dds)[
    grepl("^condition_", resultsNames(dds)) &
      grepl("_trt_vs_untrt$", resultsNames(dds))
  ][1]

  if (!length(coef_name) || is.na(coef_name) || !nzchar(coef_name)) {
    stop(sprintf("[%s] Could not find DESeq2 coefficient name for trt_vs_untrt.", dataset_name))
  }

  res_shrunk <- safe_shrink_lfc(
    dds = dds,
    contrast_vector = c("condition", "trt", "untrt"),
    coef_name = coef_name
  )

  res_df <- as.data.frame(res)
  res_df$feature_id <- rownames(res_df)

  shrunk_df <- as.data.frame(res_shrunk)
  shrunk_df$feature_id <- rownames(shrunk_df)
  names(shrunk_df)[names(shrunk_df) == "log2FoldChange"] <- "log2FoldChange_shrunken"

  merged <- dplyr::left_join(res_df, shrunk_df[, c("feature_id", "log2FoldChange_shrunken")], by = "feature_id")
  merged <- dplyr::left_join(annotation_df, merged, by = "feature_id")

  fd <- tryCatch(
    fdrtool(as.numeric(merged$stat), statistic = "normal", plot = FALSE, verbose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(fd) && length(fd$pval) == nrow(merged)) {
    merged$empirical_p <- clip_probabilities(fd$pval)

    hc_cutoff_z <- fd$param[["cutoff"]]
    if (is.finite(hc_cutoff_z) && !is.na(hc_cutoff_z)) {
      hc_p_threshold_dataset <- 2 * pnorm(-abs(hc_cutoff_z))
    } else {
      hc_p_threshold_dataset <- NA_real_
    }
  } else {
    merged$empirical_p <- clip_probabilities(merged$pvalue)
    hc_p_threshold_dataset <- compute_hc_p_threshold_dataset(merged$stat)
  }

  if (!is.finite(hc_p_threshold_dataset) ||
      is.na(hc_p_threshold_dataset) ||
      hc_p_threshold_dataset <= 0 ||
      hc_p_threshold_dataset >= max_usable_hc_p_threshold) {
    hc_p_threshold_dataset <- alpha_level
  }

  hbfss_threshold_dataset <- abs(log10(hc_p_threshold_dataset)) * lfc_boundary
  merged$pi_valueE <- compute_pi_valueE(merged$log2FoldChange_shrunken, merged$empirical_p)
  merged$HBFSS <- abs(merged$log2FoldChange_shrunken) * safe_neglog10(merged$empirical_p)

  merged$standard_significant <- !is.na(merged$padj) &
    is.finite(merged$padj) &
    (merged$padj < alpha_level) &
    (abs(merged$log2FoldChange_shrunken) >= lfc_boundary)

  merged$HBFSS_significant <- !is.na(merged$empirical_p) &
    is.finite(merged$empirical_p) &
    (merged$empirical_p < hc_p_threshold_dataset) &
    (merged$HBFSS >= hbfss_threshold_dataset)

  merged$effect_class <- classify_effect_strength(abs(merged$log2FoldChange_shrunken))
  merged <- reorder_result_columns(merged)

  summary_df <- data.frame(
    dataset_name = dataset_name,
    n_features = nrow(merged),
    n_standard_significant = sum(merged$standard_significant, na.rm = TRUE),
    n_hbfss_significant = sum(merged$HBFSS_significant, na.rm = TRUE),
    alpha_level = alpha_level,
    lfc_boundary = lfc_boundary,
    hc_p_threshold_dataset = hc_p_threshold_dataset,
    hbfss_threshold_dataset = hbfss_threshold_dataset,
    stringsAsFactors = FALSE
  )

  list(
    dds = dds,
    results_df = merged,
    summary_df = summary_df,
    hc_p_threshold_dataset = hc_p_threshold_dataset,
    hbfss_threshold_dataset = hbfss_threshold_dataset
  )
}

# -----------------------------------------------------------------------------
# VOLCANO HELPERS
# -----------------------------------------------------------------------------

is_valid_gene_label <- function(x) {
  x <- trimws(as.character(x))
  !(is.na(x) |
      x == "" |
      x %in% c("-", ".", "---", "NA", "N/A", "na", "n/a", "null", "NULL"))
}

build_volcano_classification <- function(df,
                                         hc_p_threshold_dataset,
                                         hbfss_threshold_dataset) {
  out <- df
  out$neglog10_empirical_p <- safe_neglog10(out$empirical_p)
  out$neglog10_padj <- safe_neglog10(out$padj)

  out$volcano_class <- "Not significant"
  out$volcano_class[out$standard_significant] <- "DESeq2 only"
  out$volcano_class[out$HBFSS_significant] <- "HBFSS only"
  out$volcano_class[out$standard_significant & out$HBFSS_significant] <- "Overlap"

  good_labels <- is_valid_gene_label(out$gene_name)
  out$label_candidate <- ifelse(good_labels & out$HBFSS_significant, trimws(out$gene_name), NA_character_)

  out
}

add_hbfss_boundary_layer <- function(p,
                                     x_range,
                                     hc_p_threshold_dataset,
                                     hbfss_threshold_dataset) {
  x_abs <- seq(
    from = 1e-4,
    to = max(abs(x_range), na.rm = TRUE),
    length.out = 800
  )

  hyperbola_y <- hbfss_threshold_dataset / x_abs
  floor_y <- rep(safe_neglog10(hc_p_threshold_dataset), length(x_abs))
  boundary_y <- pmax(hyperbola_y, floor_y)

  boundary_df <- dplyr::bind_rows(
    data.frame(x = -rev(x_abs), y = rev(boundary_y)),
    data.frame(x = x_abs, y = boundary_y)
  )

  p +
    geom_path(
      data = boundary_df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      colour = plot_palette$hbfss,
      linewidth = LINE_WIDTH_BOUNDARY,
      linetype = "solid"
    ) +
    geom_hline(
      yintercept = safe_neglog10(hc_p_threshold_dataset),
      colour = plot_palette$hbfss,
      linewidth = LINE_WIDTH_ZERO,
      linetype = "dotted"
    )
}

build_hbfss_volcano_plot <- function(df,
                                     dataset_name,
                                     hc_p_threshold_dataset,
                                     hbfss_threshold_dataset) {
  plot_df <- build_volcano_classification(df, hc_p_threshold_dataset, hbfss_threshold_dataset)

  x_vals <- plot_df$log2FoldChange_shrunken
  x_max <- max(abs(x_vals[is.finite(x_vals)]), na.rm = TRUE)
  if (!is.finite(x_max) || is.na(x_max) || x_max <= 0) x_max <- 2.5
  x_lim <- c(-1, 1) * max(2.5, x_max * 1.08)

  y_vals <- plot_df$neglog10_empirical_p
  y_max_data <- max(y_vals[is.finite(y_vals)], na.rm = TRUE)
  y_floor <- safe_neglog10(hc_p_threshold_dataset)
  y_max <- max(y_floor * 1.15, y_max_data * 1.05, 5)

  p <- ggplot(
    plot_df,
    aes(x = log2FoldChange_shrunken, y = neglog10_empirical_p, colour = volcano_class)
  ) +
    geom_point(size = POINT_SIZE_PRIMARY, alpha = POINT_ALPHA_PRIMARY) +
    scale_colour_manual(values = c(
      "Not significant" = plot_palette$weak,
      "DESeq2 only" = plot_palette$deseq2,
      "HBFSS only" = plot_palette$hbfss,
      "Overlap" = plot_palette$overlap
    )) +
    geom_vline(
      xintercept = c(-lfc_boundary, lfc_boundary),
      linewidth = LINE_WIDTH_ZERO,
      linetype = "dashed",
      colour = plot_palette$threshold
    ) +
    coord_cartesian(xlim = x_lim, ylim = c(0, y_max)) +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "HBFSS volcano")),
      subtitle = paste(
        final_analysis_label(),
        "| Boundary = max( empirical-p floor, HBFSS hyperbola )"
      ),
      x = "Shrunken log2 fold change",
      y = expression(-log[10]("empirical p")),
      colour = NULL
    ) +
    volcano_theme()

  p <- add_hbfss_boundary_layer(
    p = p,
    x_range = x_lim,
    hc_p_threshold_dataset = hc_p_threshold_dataset,
    hbfss_threshold_dataset = hbfss_threshold_dataset
  )

  top_labels <- plot_df %>%
    dplyr::filter(!is.na(label_candidate)) %>%
    dplyr::arrange(desc(HBFSS)) %>%
    dplyr::slice_head(n = 20)

  if (nrow(top_labels) > 0) {
    p <- p + geom_text_repel(
      data = top_labels,
      aes(label = label_candidate),
      size = 3.0,
      max.overlaps = Inf,
      box.padding = 0.22,
      point.padding = 0.12,
      segment.alpha = 0.55,
      show.legend = FALSE
    )
  }

  p
}

build_deseq2_volcano_plot <- function(df, dataset_name) {
  plot_df <- df
  plot_df$neglog10_padj <- safe_neglog10(plot_df$padj)

  x_vals <- plot_df$log2FoldChange_shrunken
  x_max <- max(abs(x_vals[is.finite(x_vals)]), na.rm = TRUE)
  if (!is.finite(x_max) || is.na(x_max) || x_max <= 0) x_max <- 2.5
  x_lim <- c(-1, 1) * max(2.5, x_max * 1.08)

  p <- ggplot(
    plot_df,
    aes(
      x = log2FoldChange_shrunken,
      y = neglog10_padj,
      colour = standard_significant
    )
  ) +
    geom_point(size = POINT_SIZE_PRIMARY, alpha = POINT_ALPHA_PRIMARY) +
    scale_colour_manual(values = c(
      `FALSE` = plot_palette$weak,
      `TRUE`  = plot_palette$deseq2
    )) +
    geom_vline(
      xintercept = c(-lfc_boundary, lfc_boundary),
      linewidth = LINE_WIDTH_ZERO,
      linetype = "dashed",
      colour = plot_palette$threshold
    ) +
    geom_hline(
      yintercept = safe_neglog10(alpha_level),
      linewidth = LINE_WIDTH_ZERO,
      linetype = "dotted",
      colour = plot_palette$deseq2
    ) +
    coord_cartesian(xlim = x_lim) +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "DESeq2 volcano")),
      subtitle = paste(
        final_analysis_label(),
        "| Standard threshold: padj < 0.10 and |LFC| >= 1"
      ),
      x = "Shrunken log2 fold change",
      y = expression(-log[10]("padj")),
      colour = NULL
    ) +
    volcano_theme()

  top_labels <- plot_df %>%
    dplyr::filter(
      standard_significant,
      is_valid_gene_label(gene_name)
    ) %>%
    dplyr::arrange(desc(abs(log2FoldChange_shrunken) * neglog10_padj)) %>%
    dplyr::slice_head(n = 20)

  if (nrow(top_labels) > 0) {
    p <- p + geom_text_repel(
      data = top_labels,
      aes(label = trimws(gene_name)),
      size = 3.0,
      max.overlaps = Inf,
      box.padding = 0.22,
      point.padding = 0.12,
      segment.alpha = 0.55,
      show.legend = FALSE
    )
  }

  p
}

# -----------------------------------------------------------------------------
# DISPERSION PLOT
# -----------------------------------------------------------------------------

build_dispersion_plot <- function(dds, dataset_name) {
  disp_df <- dplyr::bind_rows(
    data.frame(
      baseMean = mcols(dds)$baseMean,
      dispersion_value = mcols(dds)$dispGeneEst,
      dispersion_type = "Gene-wise estimate",
      stringsAsFactors = FALSE
    ),
    data.frame(
      baseMean = mcols(dds)$baseMean,
      dispersion_value = dispersions(dds),
      dispersion_type = "Shrunk estimate",
      stringsAsFactors = FALSE
    )
  )

  ggplot(disp_df, aes(x = baseMean, y = dispersion_value, colour = dispersion_type)) +
    geom_point(size = POINT_SIZE_DISP, alpha = POINT_ALPHA_DISP) +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_manual(values = c(
      "Gene-wise estimate" = plot_palette$weak,
      "Shrunk estimate" = plot_palette$deseq2
    )) +
    labs(
      title = compact_title(paste(pretty_dataset_label(dataset_name), "DESeq2 dispersion fit")),
      subtitle = final_analysis_label(),
      x = "baseMean",
      y = "Dispersion",
      colour = NULL
    ) +
    dispersion_theme()
}

# -----------------------------------------------------------------------------
# EXPORT HELPERS
# -----------------------------------------------------------------------------

export_dataset_outputs <- function(analysis_obj,
                                   out_dir,
                                   dataset_name) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  save_csv(
    analysis_obj$results_df,
    file.path(out_dir, paste0(dataset_name, "_results.csv"))
  )

  save_csv(
    analysis_obj$summary_df,
    file.path(out_dir, paste0(dataset_name, "_summary.csv"))
  )

  p_hbfss <- safe_plot_build(
    build_hbfss_volcano_plot(
      df = analysis_obj$results_df,
      dataset_name = dataset_name,
      hc_p_threshold_dataset = analysis_obj$hc_p_threshold_dataset,
      hbfss_threshold_dataset = analysis_obj$hbfss_threshold_dataset
    ),
    label = paste(dataset_name, "HBFSS volcano")
  )

  if (!is.null(p_hbfss)) {
    save_grob(
      p_hbfss,
      file.path(out_dir, paste0(dataset_name, "_HBFSS_volcano.png"))
    )
  }

  p_deseq2 <- safe_plot_build(
    build_deseq2_volcano_plot(
      df = analysis_obj$results_df,
      dataset_name = dataset_name
    ),
    label = paste(dataset_name, "DESeq2 volcano")
  )

  if (!is.null(p_deseq2)) {
    save_grob(
      p_deseq2,
      file.path(out_dir, paste0(dataset_name, "_DESeq2_volcano.png"))
    )
  }

  p_disp <- safe_plot_build(
    build_dispersion_plot(
      dds = analysis_obj$dds,
      dataset_name = dataset_name
    ),
    label = paste(dataset_name, "dispersion plot")
  )

  if (!is.null(p_disp)) {
    save_grob(
      p_disp,
      file.path(out_dir, paste0(dataset_name, "_dispersion.png"))
    )
  }

  invisible(TRUE)
}

export_comparison_evs_outputs <- function(split_obj,
                                          comparison_name,
                                          out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  save_csv(
    split_obj$fit_trt$loading_table,
    file.path(out_dir, paste0(comparison_name, "_Treatment_EVS_loadings.csv"))
  )
  save_csv(
    split_obj$fit_untrt$loading_table,
    file.path(out_dir, paste0(comparison_name, "_Control_EVS_loadings.csv"))
  )
  save_csv(
    split_obj$combined_fit$loading_table,
    file.path(out_dir, paste0(comparison_name, "_Combined_EVS_loadings.csv"))
  )
  save_csv(
    split_obj$final_rank_aggregation,
    file.path(out_dir, paste0(comparison_name, "_EVS_cutoff_summary.csv"))
  )

  save_csv(
    data.frame(feature_id = split_obj$leading_edge_ids, stringsAsFactors = FALSE),
    file.path(out_dir, paste0(comparison_name, "_leading_edge_ids.csv"))
  )
  save_csv(
    data.frame(feature_id = split_obj$remainder_ids, stringsAsFactors = FALSE),
    file.path(out_dir, paste0(comparison_name, "_remainder_ids.csv"))
  )

  p_var <- safe_plot_build(
    build_pca_variance_plot(split_obj$fit_trt, split_obj$fit_untrt, comparison_name),
    label = paste(comparison_name, "PCA variance plot")
  )
  if (!is.null(p_var)) {
    save_grob(p_var, file.path(out_dir, paste0(comparison_name, "_PCA_variance.png")))
  }

  p_scat <- safe_plot_build(
    build_pca_scatter_plot(split_obj$fit_trt, split_obj$fit_untrt, comparison_name),
    label = paste(comparison_name, "PCA scatter plot")
  )
  if (!is.null(p_scat)) {
    save_grob(p_scat, file.path(out_dir, paste0(comparison_name, "_PCA_scatter.png")))
  }

  p_trt_rank <- safe_plot_build(
    build_pc1_loading_rank_plot(
      split_obj$fit_trt$loading_table,
      comparison_name,
      "Treatment",
      cutoff_rank = split_obj$fit_trt$top_n_used
    ),
    label = paste(comparison_name, "Treatment rank plot")
  )
  if (!is.null(p_trt_rank)) {
    save_grob(p_trt_rank, file.path(out_dir, paste0(comparison_name, "_Treatment_loading_rank.png")))
  }

  p_untrt_rank <- safe_plot_build(
    build_pc1_loading_rank_plot(
      split_obj$fit_untrt$loading_table,
      comparison_name,
      "Control",
      cutoff_rank = split_obj$fit_untrt$top_n_used
    ),
    label = paste(comparison_name, "Control rank plot")
  )
  if (!is.null(p_untrt_rank)) {
    save_grob(p_untrt_rank, file.path(out_dir, paste0(comparison_name, "_Control_loading_rank.png")))
  }

  p_trt_hist <- safe_plot_build(
    build_pc1_loading_hist_plot(split_obj$fit_trt$loading_table, comparison_name, "Treatment"),
    label = paste(comparison_name, "Treatment histogram")
  )
  if (!is.null(p_trt_hist)) {
    save_grob(p_trt_hist, file.path(out_dir, paste0(comparison_name, "_Treatment_loading_hist.png")))
  }

  p_untrt_hist <- safe_plot_build(
    build_pc1_loading_hist_plot(split_obj$fit_untrt$loading_table, comparison_name, "Control"),
    label = paste(comparison_name, "Control histogram")
  )
  if (!is.null(p_untrt_hist)) {
    save_grob(p_untrt_hist, file.path(out_dir, paste0(comparison_name, "_Control_loading_hist.png")))
  }

  p_combined <- safe_plot_build(
    build_combined_loading_rank_plot(split_obj$combined_fit, comparison_name),
    label = paste(comparison_name, "Combined rank plot")
  )
  if (!is.null(p_combined)) {
    save_grob(p_combined, file.path(out_dir, paste0(comparison_name, "_Combined_loading_rank.png")))
  }

  p_nb_rank <- safe_plot_build(
    build_nb_rank_profile_plot(
      split_obj$comparison_eval$combined_rank_nb_table,
      comparison_name,
      cutoff_rank = split_obj$final_filtered_rank
    ),
    label = paste(comparison_name, "NB rank profile")
  )
  if (!is.null(p_nb_rank)) {
    save_grob(p_nb_rank, file.path(out_dir, paste0(comparison_name, "_EVS_NB_regime_rank_profiles_combined.png")))
  }

  p_nb_seg <- safe_plot_build(
    build_nb_segment_fit_plot(
      split_obj$comparison_eval$combined_rank_nb_table,
      comparison_name,
      cutoff_rank = split_obj$final_filtered_rank
    ),
    label = paste(comparison_name, "NB segment fit")
  )
  if (!is.null(p_nb_seg)) {
    save_grob(p_nb_seg, file.path(out_dir, paste0(comparison_name, "_EVS_NB_regime_segment_fit_combined.png")))
  }

  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# GIT PUSH HELPER
# -----------------------------------------------------------------------------

push_exports_to_git <- function(repo_path = repo_dir) {
  if (!isTRUE(auto_push_exports)) return(invisible(FALSE))

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(repo_path)

  git_add_status <- system("git add exports")
  if (!identical(git_add_status, 0L)) stop("git add exports failed.")

  git_commit_status <- system(sprintf("git commit -m %s", shQuote("Update EVS/HBFSS exports")))
  if (!(identical(git_commit_status, 0L) || identical(git_commit_status, 1L))) {
    stop("git commit failed.")
  }

  git_push_status <- system("git push")
  if (!identical(git_push_status, 0L)) stop("git push failed.")

  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

main <- function() {
  count_mat <- load_count_matrix(count_file)
  annotation_df <- extract_annotation(count_file)
  comparison_inputs <- build_comparison_inputs(count_mat)

  precomputed_selection <- precompute_global_evs_rank_selection(comparison_inputs)
  build_evs_candidate_export(
    precomputed_selection = precomputed_selection,
    out_dir = file.path(output_dir, "tables")
  )

  all_dataset_summaries <- list()

  for (comparison_name in names(comparison_inputs)) {
    message(sprintf("Processing %s ...", comparison_name))

    cmp_input <- comparison_inputs[[comparison_name]]
    cmp_dir <- file.path(output_dir, comparison_name)
    dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

    split_obj <- build_eigenvector_split(
      count_matrix = cmp_input$count_matrix,
      coldata = cmp_input$coldata,
      comparison_name = comparison_name,
      precomputed_selection = precomputed_selection
    )

    export_comparison_evs_outputs(
      split_obj = split_obj,
      comparison_name = comparison_name,
      out_dir = cmp_dir
    )

    retained_ids <- precomputed_selection$per_comparison[[comparison_name]]$retained_feature_ids

    dataset_map <- list(
      raw_dataset = retained_ids,
      leading_edge_dataset = split_obj$leading_edge_ids,
      remainder_dataset = split_obj$remainder_ids
    )

    for (dataset_key in names(dataset_map)) {
      dataset_ids <- intersect(dataset_map[[dataset_key]], rownames(cmp_input$count_matrix))
      dataset_name <- paste(comparison_name, dataset_key, sep = "_")
      dataset_dir <- file.path(cmp_dir, dataset_key)
      dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

      cm_subset <- cmp_input$count_matrix[dataset_ids, , drop = FALSE]

      analysis_obj <- run_deseq2_analysis(
        count_matrix_subset = cm_subset,
        coldata_subset = cmp_input$coldata,
        dataset_name = dataset_name,
        annotation_df = annotation_df
      )

      export_dataset_outputs(
        analysis_obj = analysis_obj,
        out_dir = dataset_dir,
        dataset_name = dataset_name
      )

      all_dataset_summaries[[dataset_name]] <- analysis_obj$summary_df
    }
  }

  if (length(all_dataset_summaries)) {
    combined_summary <- dplyr::bind_rows(all_dataset_summaries)
    combined_summary_path <- file.path(output_dir, "EVS_HBFSS_dataset_summary_all_comparisons.csv")
    save_csv(combined_summary, combined_summary_path)
    message(sprintf("Combined dataset summary written to: %s", combined_summary_path))
  }

  push_exports_to_git(repo_dir)
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# RUN
# -----------------------------------------------------------------------------

main()
