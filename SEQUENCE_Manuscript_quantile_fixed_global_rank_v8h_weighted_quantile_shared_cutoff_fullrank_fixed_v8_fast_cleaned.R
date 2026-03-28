# =============================================================================
# EVS + HBFSS FULL PIPELINE
# DESEQ2-INFORMED NB REGIME CUT-OFF VERSION
# REVISED: COMPARISON-SPECIFIC COMBINED RANK AXIS
#
# METHOD SUMMARY
# 1. Read raw counts.
# 2. For each comparison, normalize BEFORE EVS for PCA ranking if requested.
# 3. Compute PC1 loadings separately in treatment and control.
# 4. Build one comparison-specific combined rank using:
#      combined_loading = max(|loading_treatment|, |loading_control|)
#    unless another combined_loading_rule is selected.
# 5. Run DESeq2 ON THE FULL RAW DATASET for that comparison to estimate:
#      - baseMean
#      - dispersion (or dispGeneEst)
#    Then derive:
#      - IOD_NB  = 1 + alpha * mu
#      - CV2_NB  = 1/mu + alpha
#    using the average of condition-specific normalized means for scoring.
# 6. For each candidate cutoff rank on the combined rank axis:
#      - leading edge is scored on log(IOD_NB)
#      - remainder   is scored on log(CV2_NB)
#    using free-slope linear fits in log(baseMean_scoring) space.
# 7. Candidate split selection uses automatic fallback:
#      BIC -> AIC -> pure RSS improvement -> deterministic fallback
# 8. Regime plausibility is based on level separation:
#      - leading edge mean log(IOD_NB) > remainder mean log(IOD_NB)
#      - remainder mean log(CV2_NB)  > leading edge mean log(CV2_NB)
# 9. Freeze the cutoff and rerun DESeq2 separately on:
#      - original dataset
#      - leading-edge dataset
#      - remainder dataset
# 10. Export full tables and figures.
#
# NOTES
# - EVS ranking and DESeq2 cutoff scoring are conceptually separated:
#     * EVS ranking uses normalized counts if normalize_before_evs = TRUE
#     * cutoff scoring uses DESeq2 model-based dispersion from raw counts
# - Final DESeq2 analyses always use raw counts.
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
evs_dispersion_source <- "dispGeneEst"   # options: "dispGeneEst", "dispersion"

# fixed workflow label
evs_final_rank_mode <- "comparison_specific_combined_rank"

changepoint_min_segment_size <- 250L
changepoint_smoothing_window <- 101L

# Primary preference + automatic fallback ladder
changepoint_model_selection_mode <- "bic"   # options: "bic", "aic", "none"
changepoint_fallback_modes <- c("bic", "aic", "none")

# Regime plausibility checks
require_regime_separation <- TRUE
min_regime_separation_delta <- 0.00

# Combined comparison-level ranking rule
combined_loading_rule <- "max"  # options: "max", "sum", "l2"

# Deterministic fallback
allow_rank_fallback_if_no_valid_changepoint <- TRUE
fallback_rank_mode <- "top_n"   # options: "top_n", "fraction"
fallback_rank_top_n <- 5000L
fallback_rank_fraction <- 0.10

# Optional candidate-grid thinning for speed
use_coarse_to_fine_rank_search <- FALSE
coarse_rank_step <- 25L
fine_search_half_window <- 100L

# Git push of exports
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
    if (k < 3L) return(x)
  }

  kernel <- rep(1 / k, k)
  left_pad  <- rep(x[1], k %/% 2L)
  right_pad <- rep(x[n], k %/% 2L)
  x_pad <- c(left_pad, x, right_pad)

  out <- stats::filter(x_pad, filter = kernel, sides = 2, method = "convolution")
  out <- as.numeric(out[(length(left_pad) + 1L):(length(left_pad) + n)])

  if (anyNA(out)) out[is.na(out)] <- x[is.na(out)]
  out
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
  if (!is.finite(idx) || is.na(idx) || idx < 1L) return(NA_real_)
  x[idx]
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
    stop(sprintf(
      "[%s] No valid changepoint rank could be determined, and fallback is disabled.",
      label
    ))
  }

  if (identical(fallback_rank_mode, "top_n")) {
    fallback_rank <- as.integer(fallback_rank_top_n)
    fallback_rank <- max(1L, min(n_total - 1L, fallback_rank))
    warning(sprintf(
      "[%s] No valid changepoint rank found. Falling back to fixed top-N rank %d.",
      label, fallback_rank
    ))
    return(fallback_rank)
  }

  if (identical(fallback_rank_mode, "fraction")) {
    fallback_rank <- as.integer(round(fallback_rank_fraction * n_total))
    fallback_rank <- max(1L, min(n_total - 1L, fallback_rank))
    warning(sprintf(
      "[%s] No valid changepoint rank found. Falling back to rank %d (fraction %.3f of %d features).",
      label, fallback_rank, fallback_rank_fraction, n_total
    ))
    return(fallback_rank)
  }

  stop(sprintf("[%s] Unknown fallback_rank_mode: %s", label, fallback_rank_mode))
}

make_safe_rank <- function(rank_index, n_total) {
  if (!is.finite(rank_index) || is.na(rank_index) || n_total < 2L) return(NA_integer_)
  rank_index <- as.integer(rank_index)[1]
  max(1L, min(as.integer(n_total) - 1L, rank_index))
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

  if (identical(rule, "max")) {
    return(pmax(abs_trt, abs_ctrl))
  }
  if (identical(rule, "sum")) {
    return(abs_trt + abs_ctrl)
  }
  if (identical(rule, "l2")) {
    return(sqrt(abs_trt^2 + abs_ctrl^2))
  }
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

plot_expand_xy <- function() {
  list(
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.24))),
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.30)))
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

strip_dataset_key <- function(dataset_name) {
  sub("^.*_(raw_dataset|leading_edge_dataset|remainder_dataset)$", "\\1", dataset_name)
}

comparison_label_from_dataset_name <- function(dataset_name) {
  sub("_(raw_dataset|leading_edge_dataset|remainder_dataset)$", "", dataset_name)
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
    combined  = "Combined",
    group_label
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
    "feature_id",
    "gene_symbol",
    "dataset_name",
    "regulation_direction",
    "effect_class",
    "standard_significant",
    "HBFSS_significant",
    "HBFSS_only_call",
    "deseq2_strong_call",
    "deseq2_weak_call",
    "overlap_call",
    "lfc_shrunk",
    "log2FoldChange",
    "HBFSS",
    "padj",
    "empirical_p",
    "empirical_q",
    "lfdr",
    "padj_strong_effect",
    "padj_weak_effect",
    "baseMean",
    "dispGeneEst",
    "dispFit",
    "dispersion"
  )

  front <- intersect(preferred_front_cols, names(final_df))
  back  <- setdiff(names(final_df), front)
  final_df[, c(front, back), drop = FALSE]
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

# -----------------------------------------------------------------------------
# EMPIRICAL NULL / HBFSS HELPERS
# -----------------------------------------------------------------------------

run_empirical_null_fdrtool <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  stat_vec <- unname(stat_vec)

  if (length(stat_vec) < 5L) {
    stop(sprintf("[%s] Fewer than 5 finite Wald statistics were available for fdrtool.", dataset_name))
  }

  fit <- tryCatch(
    fdrtool::fdrtool(
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
        fdrtool::fdrtool(
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
    na.last = NA,
    decreasing = FALSE
  )
  if (length(sorted_empirical_p) < 5L) return(NA_real_)

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

# -----------------------------------------------------------------------------
# METADATA
# -----------------------------------------------------------------------------

meta_all <- data.frame(
  id = c(
    "R0_1","R0_2","R0_3","R0_4","R0_5","ZT6_1","ZT6_2","ZT6_3","ZT6_4","ZT6_5",
    "R2_1","R2_2","R2_3","R2_4","R2_5","ZT8_1","ZT8_2","ZT8_3","ZT8_4","ZT8_5",
    "R4_1","R4_2","R4_3","R4_4","R4_5","ZT10_1","ZT10_2","ZT10_3","ZT10_4","ZT10_5",
    "R8_1","R8_2","R8_3","R8_4","R8_5","ZT14_1","ZT14_2","ZT14_3","ZT14_4","ZT14_5"
  ),
  condition = c(
    rep("treatment", 5), rep("control", 5),
    rep("treatment", 5), rep("control", 5),
    rep("treatment", 5), rep("control", 5),
    rep("treatment", 5), rep("control", 5)
  ),
  stringsAsFactors = FALSE
)

rownames(meta_all) <- meta_all$id
meta_all$condition <- factor(meta_all$condition, levels = c("control", "treatment"))
levels(meta_all$condition) <- c("untrt", "trt")

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix   = c("R0", "R2", "R4", "R8"),
  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# DATA IMPORT
# -----------------------------------------------------------------------------

WTTS_Seq <- read.csv(
  count_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

WTTS_Seq <- as.data.frame(WTTS_Seq, stringsAsFactors = FALSE)
WTTS_Seq$OrigID <- as.character(WTTS_Seq$OrigID)
WTTS_Seq$Symbol <- as.character(WTTS_Seq$Symbol)

assert_required_columns(WTTS_Seq, c("OrigID", "Symbol"), object_name = "WTTS count file")
assert_required_columns(WTTS_Seq, meta_all$id, object_name = "WTTS count file sample columns")

WTTS_Seq <- WTTS_Seq[!is.na(WTTS_Seq$OrigID), , drop = FALSE]
WTTS_Seq <- WTTS_Seq[rowSums(is.na(WTTS_Seq[, meta_all$id, drop = FALSE])) == 0, , drop = FALSE]
rownames(WTTS_Seq) <- WTTS_Seq$OrigID

OrigID_Symbol <- unique(WTTS_Seq[, c("OrigID", "Symbol"), drop = FALSE])
colnames(OrigID_Symbol) <- c("feature_id", "gene_symbol")
OrigID_Symbol$feature_id <- as.character(OrigID_Symbol$feature_id)
OrigID_Symbol$gene_symbol <- as.character(OrigID_Symbol$gene_symbol)

OrigID_Symbol <- OrigID_Symbol %>%
  dplyr::mutate(gene_symbol = dplyr::if_else(is.na(gene_symbol), "", trimws(gene_symbol))) %>%
  dplyr::arrange(feature_id, dplyr::desc(gene_symbol != ""), gene_symbol) %>%
  dplyr::distinct(feature_id, .keep_all = TRUE) %>%
  dplyr::mutate(gene_symbol = dplyr::na_if(gene_symbol, ""))

# -----------------------------------------------------------------------------
# COMPARISON PREP
# -----------------------------------------------------------------------------

prepare_comparison_data <- function(comparison_name, group1_prefix, group2_prefix, WTTS_Seq, meta_all) {
  keep_ids <- grepl(paste0("^", group1_prefix, "_"), meta_all$id) |
    grepl(paste0("^", group2_prefix, "_"), meta_all$id)

  meta_sub <- meta_all[keep_ids, , drop = FALSE]
  coldata <- meta_sub[, "condition", drop = FALSE]
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
  count_sub <- count_sub[, rownames(coldata), drop = FALSE]
  stopifnot(all(colnames(count_sub) == rownames(coldata)))

  list(
    comparison_name = comparison_name,
    count_matrix = as.matrix(count_sub),
    coldata = coldata
  )
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

# =============================================================================
# EVS PCA HELPERS
# =============================================================================

compute_pc1_loading_table <- function(value_df, sample_names,
                                      preprocessing_label,
                                      selected_reason = "external_rank") {
  x <- log1p(as.matrix(value_df[, sample_names, drop = FALSE]))
  storage.mode(x) <- "double"

  if (nrow(x) < 2L || ncol(x) < 2L) {
    stop(sprintf(
      "PCA failed for '%s': need at least 2 features and 2 samples, found %d features and %d samples.",
      preprocessing_label, nrow(x), ncol(x)
    ))
  }

  max_pcs <- min(5L, ncol(x))
  if (!is.finite(max_pcs) || is.na(max_pcs) || max_pcs < 1L) {
    stop(sprintf("PCA failed for '%s': invalid rank requested.", preprocessing_label))
  }

  pca_fit <- prcomp(t(x), scale. = FALSE, rank. = max_pcs)

  if (is.null(pca_fit$rotation) || nrow(pca_fit$rotation) < 1L || ncol(pca_fit$rotation) < 1L) {
    stop(sprintf("PCA failed for '%s': rotation matrix was empty.", preprocessing_label))
  }

  loading_abs <- abs(pca_fit$rotation[, 1])
  loading_tbl <- data.frame(
    feature_id = names(loading_abs),
    pc1_loading = unname(pca_fit$rotation[, 1]),
    pc1_loading_abs = unname(loading_abs),
    stringsAsFactors = FALSE
  )

  loading_tbl <- loading_tbl[order(loading_tbl$pc1_loading_abs, decreasing = TRUE), , drop = FALSE]
  loading_tbl$rank <- seq_len(nrow(loading_tbl))
  loading_tbl$split_class <- "background_loading"

  list(
    pca_fit = pca_fit,
    loading_table = loading_tbl,
    cutoff = NA_real_,
    top_n_used = NA_integer_,
    cutoff_quantile = NA_real_,
    cutoff_method = "deseq2_nb_regime_group_rank",
    selected_reason = selected_reason,
    preprocessing_label = preprocessing_label
  )
}

build_combined_pc1_loading_table <- function(fit_trt, fit_untrt,
                                             preprocessing_label,
                                             combined_rule = combined_loading_rule) {
  trt_tbl <- fit_trt$loading_table[, c("feature_id", "pc1_loading", "pc1_loading_abs"), drop = FALSE]
  ctrl_tbl <- fit_untrt$loading_table[, c("feature_id", "pc1_loading", "pc1_loading_abs"), drop = FALSE]

  names(trt_tbl)  <- c("feature_id", "pc1_loading_trt",  "pc1_loading_abs_trt")
  names(ctrl_tbl) <- c("feature_id", "pc1_loading_ctrl", "pc1_loading_abs_ctrl")

  combined_tbl <- dplyr::full_join(trt_tbl, ctrl_tbl, by = "feature_id")
  combined_tbl$pc1_loading_abs_trt  <- coalesce_numeric(combined_tbl$pc1_loading_abs_trt,  0)
  combined_tbl$pc1_loading_abs_ctrl <- coalesce_numeric(combined_tbl$pc1_loading_abs_ctrl, 0)
  combined_tbl$pc1_loading_trt      <- coalesce_numeric(combined_tbl$pc1_loading_trt,      0)
  combined_tbl$pc1_loading_ctrl     <- coalesce_numeric(combined_tbl$pc1_loading_ctrl,     0)

  combined_tbl$combined_pc1_loading_abs <- compute_combined_loading_score(
    abs_trt = combined_tbl$pc1_loading_abs_trt,
    abs_ctrl = combined_tbl$pc1_loading_abs_ctrl,
    rule = combined_rule
  )

  combined_tbl <- combined_tbl[order(combined_tbl$combined_pc1_loading_abs, decreasing = TRUE), , drop = FALSE]
  combined_tbl$rank <- seq_len(nrow(combined_tbl))
  combined_tbl$split_class <- "background_loading"

  list(
    loading_table = combined_tbl,
    cutoff = NA_real_,
    top_n_used = NA_integer_,
    cutoff_quantile = NA_real_,
    cutoff_method = "deseq2_nb_regime_combined_rank",
    selected_reason = "comparison_specific_combined_rank",
    preprocessing_label = preprocessing_label,
    combined_rule = combined_rule
  )
}

apply_combined_loading_cutoff <- function(fit_obj, rank_index, selected_reason = "combined_rank_selected") {
  loading_tbl <- fit_obj$loading_table
  loading_tbl <- loading_tbl[order(loading_tbl$rank), , drop = FALSE]

  rank_index <- resolve_rank_or_fail(rank_index, nrow(loading_tbl), label = selected_reason)
  cutoff_value <- as.numeric(loading_tbl$combined_pc1_loading_abs[rank_index])

  loading_tbl$split_class <- ifelse(
    loading_tbl$rank <= rank_index,
    "high_loading",
    "background_loading"
  )

  fit_obj$loading_table <- loading_tbl
  fit_obj$cutoff <- cutoff_value
  fit_obj$top_n_used <- rank_index
  fit_obj$cutoff_quantile <- 1 - (rank_index / nrow(loading_tbl))
  fit_obj$cutoff_method <- "deseq2_nb_regime_combined_rank"
  fit_obj$selected_reason <- selected_reason
  fit_obj
}

# =============================================================================
# DESEQ2-INFORMED NB REGIME CUT-OFF HELPERS
# =============================================================================

build_full_dataset_deseq2_dispersion_table <- function(count_matrix, coldata,
                                                       dataset_label = "full_dataset_for_cutoff") {
  count_matrix <- as.matrix(count_matrix)
  warn_if_noninteger_counts(count_matrix, label = dataset_label)

  dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData   = coldata,
    design    = DESIGN_FORMULA
  )

  dds <- dds[rowSums(counts(dds)) > 0, ]
  if (nrow(dds) < 5L) {
    stop(sprintf("[%s] Too few nonzero features for DESeq2-based cutoff scoring.", dataset_label))
  }

  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, quiet = TRUE)

  mm <- as.data.frame(mcols(dds))
  mm$feature_id <- rownames(mm)

  keep_cols <- intersect(
    c("feature_id", "baseMean", "dispGeneEst", "dispFit", "dispersion"),
    colnames(mm)
  )
  mm <- mm[, keep_cols, drop = FALSE]

  required_cols <- c("feature_id", "baseMean", evs_dispersion_source)
  assert_required_columns(mm, required_cols, object_name = paste0(dataset_label, " dispersion table"))

  names(mm)[names(mm) == evs_dispersion_source] <- "alpha_used"

  mm$baseMean <- as.numeric(mm$baseMean)
  mm$alpha_used <- as.numeric(mm$alpha_used)

  norm_counts <- counts(dds, normalized = TRUE)
  trt_cols    <- colnames(norm_counts)[coldata$condition == "trt"]
  untrt_cols  <- colnames(norm_counts)[coldata$condition == "untrt"]

  mm$mean_trt <- rowMeans(norm_counts[mm$feature_id, trt_cols, drop = FALSE])
  mm$mean_untrt <- rowMeans(norm_counts[mm$feature_id, untrt_cols, drop = FALSE])

  # comparison-level scoring mean
  mm$baseMean_scoring <- (mm$mean_trt + mm$mean_untrt) / 2

  mm <- mm[
    is.finite(mm$baseMean_scoring) & !is.na(mm$baseMean_scoring) & (mm$baseMean_scoring > 0) &
      is.finite(mm$alpha_used) & !is.na(mm$alpha_used) & (mm$alpha_used > 0),
    ,
    drop = FALSE
  ]

  if (nrow(mm) < 5L) {
    stop(sprintf("[%s] Too few valid features after filtering scoring mean/dispersion.", dataset_label))
  }

  mm$iod_nb <- 1 + (mm$alpha_used * mm$baseMean_scoring)
  mm$cv2_nb <- (1 / mm$baseMean_scoring) + mm$alpha_used

  mm$log_baseMean <- log(mm$baseMean_scoring)
  mm$log_iod_nb   <- log(mm$iod_nb)
  mm$log_cv2_nb   <- log(mm$cv2_nb)

  mm
}

build_ranked_nb_regime_table <- function(rank_order_ids,
                                         deseq2_disp_tbl,
                                         smoothing_window = changepoint_smoothing_window) {
  rank_order_ids <- as.character(rank_order_ids)

  original_rank_tbl <- data.frame(
    feature_id = rank_order_ids,
    original_rank_index = seq_along(rank_order_ids),
    stringsAsFactors = FALSE
  )

  tbl <- dplyr::left_join(original_rank_tbl, deseq2_disp_tbl, by = "feature_id")
  tbl <- tbl[
    is.finite(tbl$baseMean_scoring) & !is.na(tbl$baseMean_scoring) &
      is.finite(tbl$alpha_used) & !is.na(tbl$alpha_used),
    ,
    drop = FALSE
  ]

  if (nrow(tbl) < 5L) return(data.frame())

  tbl$filtered_rank_index <- seq_len(nrow(tbl))

  if (is.finite(smoothing_window) && !is.na(smoothing_window) && smoothing_window > 1L) {
    tbl$log_baseMean_smooth <- fast_centered_rolling_mean(tbl$log_baseMean, smoothing_window)
    tbl$log_iod_nb_smooth   <- fast_centered_rolling_mean(tbl$log_iod_nb, smoothing_window)
    tbl$log_cv2_nb_smooth   <- fast_centered_rolling_mean(tbl$log_cv2_nb, smoothing_window)
  } else {
    tbl$log_baseMean_smooth <- tbl$log_baseMean
    tbl$log_iod_nb_smooth   <- tbl$log_iod_nb
    tbl$log_cv2_nb_smooth   <- tbl$log_cv2_nb
  }

  tbl
}

fit_free_slope_segment <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  ok <- is.finite(x) & !is.na(x) & is.finite(y) & !is.na(y)
  x <- x[ok]
  y <- y[ok]

  n <- length(x)
  if (n < 3L) {
    return(list(
      slope = NA_real_,
      intercept = NA_real_,
      rss = Inf
    ))
  }

  x_mean <- mean(x)
  y_mean <- mean(y)
  sxx <- sum((x - x_mean)^2)

  if (!is.finite(sxx) || sxx <= 0) {
    slope_hat <- 0
  } else {
    slope_hat <- sum((x - x_mean) * (y - y_mean)) / sxx
  }

  intercept_hat <- y_mean - slope_hat * x_mean
  resid <- y - (intercept_hat + slope_hat * x)
  rss <- sum(resid^2)

  list(
    slope = slope_hat,
    intercept = intercept_hat,
    rss = rss
  )
}

compute_segmented_model_score <- function(split_rss, null_rss, n_obs,
                                          split_n_params = 5L,
                                          null_n_params = 4L,
                                          mode = changepoint_model_selection_mode) {
  split_rss <- as.numeric(split_rss)[1]
  null_rss  <- as.numeric(null_rss)[1]
  n_obs     <- as.integer(n_obs)[1]

  if (!is.finite(split_rss) || !is.finite(null_rss) || n_obs < 5L) {
    return(NA_real_)
  }

  split_rss <- max(split_rss, 1e-12)
  null_rss  <- max(null_rss,  1e-12)

  if (identical(mode, "none")) {
    return(null_rss - split_rss)
  }

  if (identical(mode, "bic")) {
    split_ic <- n_obs * log(split_rss / n_obs) + split_n_params * log(n_obs)
    null_ic  <- n_obs * log(null_rss  / n_obs) + null_n_params  * log(n_obs)
    return(null_ic - split_ic)
  }

  if (identical(mode, "aic")) {
    split_ic <- n_obs * log(split_rss / n_obs) + 2 * split_n_params
    null_ic  <- n_obs * log(null_rss  / n_obs) + 2 * null_n_params
    return(null_ic - split_ic)
  }

  stop(sprintf("Unknown changepoint model selection mode: %s", mode))
}

fit_nb_regime_split_score <- function(rank_tbl,
                                      rank_index,
                                      min_segment_size = changepoint_min_segment_size,
                                      mode = changepoint_model_selection_mode,
                                      min_delta = min_regime_separation_delta,
                                      require_separation = require_regime_separation) {
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
      stringsAsFactors = FALSE
    )
  }

  if (!is.finite(rank_index) || is.na(rank_index) || n_total < (2L * min_segment_size + 1L)) {
    return(empty_row(rank_index))
  }

  rank_index <- as.integer(rank_index)[1]
  if (rank_index <= 0L || rank_index >= n_total) {
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
  fit_rem  <- fit_free_slope_segment(rem_df$log_baseMean_smooth,  rem_df$log_cv2_nb_smooth)
  split_rss <- fit_lead$rss + fit_rem$rss

  fit_null_iod <- fit_free_slope_segment(rank_tbl$log_baseMean_smooth, rank_tbl$log_iod_nb_smooth)
  fit_null_cv2 <- fit_free_slope_segment(rank_tbl$log_baseMean_smooth, rank_tbl$log_cv2_nb_smooth)
  null_rss <- fit_null_iod$rss + fit_null_cv2$rss

  improvement_vs_null <- compute_segmented_model_score(
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

  regime_plausible <- is.finite(lead_iod_mean) &&
    is.finite(rem_iod_mean) &&
    is.finite(lead_cv2_mean) &&
    is.finite(rem_cv2_mean) &&
    ((lead_iod_mean - rem_iod_mean) > min_delta) &&
    ((rem_cv2_mean - lead_cv2_mean) > min_delta)

  if (!isTRUE(require_separation)) regime_plausible <- TRUE

  valid_split <- is.finite(improvement_vs_null) &&
    (improvement_vs_null > 0) &&
    regime_plausible

  data.frame(
    filtered_rank_index = rank_index,
    cutoff_feature_id = as.character(rank_tbl$feature_id[rank_index]),
    original_rank_index = as.integer(rank_tbl$original_rank_index[rank_index]),
    objective_score = as.numeric(improvement_vs_null)[1],
    valid_split = valid_split,
    lead_rss = fit_lead$rss,
    rem_rss = fit_rem$rss,
    lead_slope = fit_lead$slope,
    rem_slope = fit_rem$slope,
    null_rss = null_rss,
    improvement_vs_null = as.numeric(improvement_vs_null)[1],
    regime_plausible = regime_plausible,
    lead_iod_mean = lead_iod_mean,
    rem_iod_mean = rem_iod_mean,
    lead_cv2_mean = lead_cv2_mean,
    rem_cv2_mean = rem_cv2_mean,
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
  combined_rank_ids <- as.character(
    combined_fit$loading_table$feature_id[order(combined_fit$loading_table$rank)]
  )

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

    mode_grid <- dplyr::bind_rows(score_rows)
    mode_grid <- mode_grid[order(-mode_grid$objective_score, mode_grid$filtered_rank_index), , drop = FALSE]
    all_mode_tables[[mode_i]] <- mode_grid

    mode_diag_tables[[mode_i]] <- data.frame(
      comparison_name = comparison_name,
      model_selection_mode = mode_i,
      n_candidates = nrow(mode_grid),
      n_positive_score = sum(is.finite(mode_grid$objective_score) & mode_grid$objective_score > 0, na.rm = TRUE),
      n_regime_plausible = sum(mode_grid$regime_plausible, na.rm = TRUE),
      n_valid_split = sum(mode_grid$valid_split, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    valid_grid <- mode_grid[
      mode_grid$valid_split &
        is.finite(mode_grid$objective_score) &
        !is.na(mode_grid$cutoff_feature_id),
      ,
      drop = FALSE
    ]

    if (nrow(valid_grid) > 0) {
      if (isTRUE(use_coarse_to_fine_rank_search) && fine_search_half_window > 0L) {
        coarse_best <- as.integer(valid_grid$filtered_rank_index[1])
        lo <- max(min_segment_size, coarse_best - as.integer(fine_search_half_window))
        hi <- min(shared_n_total - min_segment_size, coarse_best + as.integer(fine_search_half_window))
        fine_grid <- seq.int(lo, hi, by = 1L)

        fine_rows <- lapply(fine_grid, function(r) {
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
        mode_grid <- fine_grid_df
        all_mode_tables[[mode_i]] <- fine_grid_df
        mode_diag_tables[[mode_i]] <- data.frame(
          comparison_name = comparison_name,
          model_selection_mode = paste0(mode_i, "_refined"),
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

      winning_mode <- mode_i
      winning_valid_grid <- valid_grid
      break
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
        stringsAsFactors = FALSE
      ),
      candidate_grid = candidate_grid,
      mode_diagnostics = mode_diagnostics,
      best_rank = NA_integer_,
      best_score = NA_real_,
      best_mode = "fallback_only",
      best_feature_id = NA_character_,
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

choose_final_rank_for_comparison <- function(comparison_name, precomputed_selection) {
  cmp_obj <- precomputed_selection$per_comparison[[comparison_name]]
  if (is.null(cmp_obj)) {
    stop(sprintf("No precomputed comparison object found for %s", comparison_name))
  }

  cmp_eval <- cmp_obj$comparison_eval
  n_total <- nrow(cmp_obj$combined_fit$loading_table)

  rank_to_use <- cmp_eval$best_rank

  if (!is.finite(rank_to_use) || is.na(rank_to_use)) {
    rank_to_use <- resolve_rank_or_fail(
      rank_index = NA_integer_,
      n_total = n_total,
      label = paste0(comparison_name, "_fallback")
    )
  }

  make_safe_rank(rank_to_use, n_total)
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

  independently_best_rank <- comparison_eval$best_rank

  if (!is.finite(final_shared_rank) || is.na(final_shared_rank)) {
    final_shared_rank <- independently_best_rank
  }
  final_shared_rank <- make_safe_rank(final_shared_rank, nrow(combined_fit$loading_table))

  if (!is.finite(final_shared_rank) || is.na(final_shared_rank)) {
    final_shared_rank <- resolve_rank_or_fail(
      rank_index = NA_integer_,
      n_total = nrow(combined_fit$loading_table),
      label = paste0(comparison_name, "_final_rank")
    )
  }

  final_reason <- if (normalize_before_evs) {
    "deseq2_nb_regime_cutoff_selected_from_full_raw_dataset_after_normalized_EVS_combined_ranking"
  } else {
    "deseq2_nb_regime_cutoff_selected_from_full_raw_dataset_after_raw_EVS_combined_ranking"
  }

  combined_fit <- apply_combined_loading_cutoff(combined_fit, final_shared_rank, final_reason)

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

  # project the combined cutoff value into the per-group plots
  trt_cutoff_val <- combined_fit$cutoff
  ctrl_cutoff_val <- combined_fit$cutoff

  fit_trt$cutoff <- trt_cutoff_val
  fit_trt$top_n_used <- sum(fit_trt$loading_table$pc1_loading_abs >= trt_cutoff_val, na.rm = TRUE)
  fit_trt$cutoff_quantile <- 1 - (fit_trt$top_n_used / nrow(fit_trt$loading_table))
  fit_trt$selected_reason <- "projected_from_combined_cutoff"

  fit_untrt$cutoff <- ctrl_cutoff_val
  fit_untrt$top_n_used <- sum(fit_untrt$loading_table$pc1_loading_abs >= ctrl_cutoff_val, na.rm = TRUE)
  fit_untrt$cutoff_quantile <- 1 - (fit_untrt$top_n_used / nrow(fit_untrt$loading_table))
  fit_untrt$selected_reason <- "projected_from_combined_cutoff"

  fit_trt_raw_view$cutoff <- trt_cutoff_val
  fit_trt_raw_view$top_n_used <- sum(fit_trt_raw_view$loading_table$pc1_loading_abs >= trt_cutoff_val, na.rm = TRUE)
  fit_trt_raw_view$cutoff_quantile <- 1 - (fit_trt_raw_view$top_n_used / nrow(fit_trt_raw_view$loading_table))
  fit_trt_raw_view$selected_reason <- "projected_from_combined_cutoff"

  fit_untrt_raw_view$cutoff <- ctrl_cutoff_val
  fit_untrt_raw_view$top_n_used <- sum(fit_untrt_raw_view$loading_table$pc1_loading_abs >= ctrl_cutoff_val, na.rm = TRUE)
  fit_untrt_raw_view$cutoff_quantile <- 1 - (fit_untrt_raw_view$top_n_used / nrow(fit_untrt_raw_view$loading_table))
  fit_untrt_raw_view$selected_reason <- "projected_from_combined_cutoff"

  final_rank_aggregation <- data.frame(
    comparison_name = comparison_name,
    independently_best_rank = independently_best_rank,
    independently_best_local_lagrange_objective = comparison_eval$best_score,
    final_rank = final_shared_rank,
    final_rank_mode = "comparison_specific_combined_rank",
    aggregation_method = if (normalize_before_evs) {
      "deseq2_nb_regime_scoring_with_normalized_combined_EVS_ranking"
    } else {
      "deseq2_nb_regime_scoring_with_raw_combined_EVS_ranking"
    },
    model_selection_mode_used = comparison_eval$best_mode,
    combined_loading_rule = combined_loading_rule,
    stringsAsFactors = FALSE
  )

  evs_cutoff_summary <- data.frame(
    preprocessing = combined_fit$preprocessing_label,
    group = "comparison_specific_combined_rank",
    cutoff_mode = combined_fit$cutoff_method,
    empiric_rank_selected = combined_fit$top_n_used,
    cutoff_quantile = combined_fit$cutoff_quantile,
    selected_reason = combined_fit$selected_reason,
    combined_loading_rule = combined_loading_rule,
    model_selection_mode_used = comparison_eval$best_mode,
    stringsAsFactors = FALSE
  )

  list(
    fit_trt = fit_trt,
    fit_untrt = fit_untrt,
    fit_trt_raw_view = fit_trt_raw_view,
    fit_untrt_raw_view = fit_untrt_raw_view,
    combined_fit = combined_fit,
    final_rank_aggregation = final_rank_aggregation,
    evs_cutoff_summary = evs_cutoff_summary,
    raw_dataset = count_matrix,
    leading_edge_dataset = count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset = if (length(remainder_ids)) count_matrix[remainder_ids, , drop = FALSE] else count_matrix[0, , drop = FALSE]
  )
}

# =============================================================================
# CORE DESEQ2 / HBFSS ANALYSIS
# =============================================================================

classify_effect_strength <- function(res_strong_padj, res_weak_padj, alpha = alpha_level) {
  out <- rep("intermediate", length(res_strong_padj))
  out[!is.na(res_weak_padj)   & res_weak_padj   < alpha] <- "weak_effect"
  out[!is.na(res_strong_padj) & res_strong_padj < alpha] <- "strong_effect"
  out
}

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
  bp <- get_bpparam()
  use_parallel <- !is.null(bp) && n_workers > 1L

  count_mat <- as.matrix(count_mat)
  if (nrow(count_mat) < 1L) {
    stop(sprintf("[%s] Dataset is empty after filtering/splitting.", dataset_name))
  }

  warn_if_noninteger_counts(count_mat, label = dataset_name)

  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData   = coldata,
    design    = DESIGN_FORMULA
  )

  dds <- dds[rowSums(counts(dds)) > 0, ]
  if (nrow(dds) < 1L) {
    stop(sprintf("[%s] No nonzero features remained after row-sum filtering.", dataset_name))
  }

  dds <- DESeq(
    dds,
    betaPrior = FALSE,
    quiet = TRUE,
    parallel = use_parallel,
    BPPARAM = bp
  )

  res <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    alpha = alpha_level,
    parallel = use_parallel,
    BPPARAM = bp
  )

  res_strong <- results(
    dds,
    contrast      = c("condition", "trt", "untrt"),
    lfcThreshold  = lfc_boundary,
    altHypothesis = "greaterAbs",
    parallel      = use_parallel,
    BPPARAM       = bp
  )

  res_weak <- results(
    dds,
    contrast      = c("condition", "trt", "untrt"),
    lfcThreshold  = lfc_boundary,
    altHypothesis = "lessAbs",
    parallel      = use_parallel,
    BPPARAM       = bp
  )

  res_all_df <- as.data.frame(res)
  res_all_df$feature_id <- as.character(rownames(res_all_df))

  valid_stat <- is.finite(res_all_df$stat) & !is.na(res_all_df$stat)
  stat_vec <- unname(as.numeric(res_all_df$stat[valid_stat]))

  fdr_fit <- tryCatch(
    run_empirical_null_fdrtool(stat_vec, dataset_name = dataset_name),
    error = function(e) {
      stop(sprintf("[%s] Empirical-null estimation failed: %s", dataset_name, conditionMessage(e)))
    }
  )

  res_df <- res_all_df
  n_valid <- sum(valid_stat)

  len_p <- length(fdr_fit$pval)
  len_q <- length(fdr_fit$qval)
  len_l <- length(fdr_fit$lfdr)

  if (len_p != n_valid || len_q != n_valid || len_l != n_valid) {
    stop(sprintf(
      "[%s] fdrtool output length mismatch: expected %d finite statistics, got p=%d, q=%d, lfdr=%d.",
      dataset_name, n_valid, len_p, len_q, len_l
    ))
  }

  res_df$empirical_p <- NA_real_
  res_df$empirical_q <- NA_real_
  res_df$lfdr <- NA_real_

  res_df$empirical_p[valid_stat] <- as.numeric(fdr_fit$pval)
  res_df$empirical_q[valid_stat] <- as.numeric(fdr_fit$qval)
  res_df$lfdr[valid_stat] <- as.numeric(fdr_fit$lfdr)

  rn <- resultsNames(dds)
  coef_idx <- grep("^condition_", rn)
  if (length(coef_idx) < 1L) {
    stop(sprintf("[%s] No condition coefficient found in resultsNames(dds): %s", dataset_name, paste(rn, collapse = ", ")))
  }
  coef_name <- rn[coef_idx[1]]

  shrink_args <- list(
    dds      = dds,
    coef     = coef_name,
    type     = "apeglm",
    res      = res,
    parallel = use_parallel,
    BPPARAM  = bp
  )

  if (isTRUE(use_fast_apeglm)) {
    shrink_args$apeMethod <- "nbinomC"
  }

  shr <- tryCatch(
    do.call(lfcShrink, shrink_args),
    error = function(e) {
      message(sprintf("[%s] Fast apeglm shrinkage failed or is unsupported; retrying without apeMethod. Error: %s",
                      dataset_name, conditionMessage(e)))
      shrink_args$apeMethod <- NULL
      do.call(lfcShrink, shrink_args)
    }
  )

  shr_df <- as.data.frame(shr)
  shr_df$feature_id <- as.character(rownames(shr_df))

  res_df <- dplyr::left_join(
    res_df,
    shr_df[, c("feature_id", "log2FoldChange")],
    by = "feature_id",
    suffix = c("", "_shrunk")
  )
  colnames(res_df)[colnames(res_df) == "log2FoldChange_shrunk"] <- "lfc_shrunk"

  hc_p_threshold_dataset <- safe_hc_thresh(res_df$empirical_p, dataset_name = dataset_name)
  hbfss_threshold_dataset <- NA_real_

  empirical_p_floored <- ifelse(
    is.na(res_df$empirical_p),
    NA_real_,
    pmax(res_df$empirical_p, 1e-300)
  )

  res_df$HBFSS <- NA_real_
  hbfss_ok <- !is.na(res_df$lfc_shrunk) & !is.na(empirical_p_floored) & is.finite(empirical_p_floored)
  res_df$HBFSS[hbfss_ok] <- abs(res_df$lfc_shrunk[hbfss_ok] * log10(empirical_p_floored[hbfss_ok]))

  if (is.na(hc_p_threshold_dataset)) {
    res_df$passes_hc_p_gate  <- FALSE
    res_df$HBFSS_significant <- FALSE
  } else {
    hbfss_threshold_dataset  <- abs(log10(hc_p_threshold_dataset)) * lfc_boundary
    res_df$passes_hc_p_gate  <- !is.na(res_df$empirical_p) &
      !is.na(res_df$lfc_shrunk) &
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

  res_df$raw_lfc_pass <- !is.na(res_df$log2FoldChange) & (abs(res_df$log2FoldChange) >= lfc_boundary)
  res_df$shrunk_lfc_pass <- !is.na(res_df$lfc_shrunk) & (abs(res_df$lfc_shrunk) >= lfc_boundary)

  res_strong_df <- as.data.frame(res_strong)
  res_strong_df$feature_id <- as.character(rownames(res_strong_df))

  res_weak_df <- as.data.frame(res_weak)
  res_weak_df$feature_id <- as.character(rownames(res_weak_df))

  res_df <- dplyr::left_join(
    res_df,
    res_strong_df[, c("feature_id", "padj")],
    by = "feature_id",
    suffix = c("", "_strong")
  )

  res_df <- dplyr::left_join(
    res_df,
    res_weak_df[, c("feature_id", "padj")],
    by = "feature_id",
    suffix = c("", "_weak")
  )

  colnames(res_df)[colnames(res_df) == "padj_strong"] <- "padj_strong_effect"
  colnames(res_df)[colnames(res_df) == "padj_weak"] <- "padj_weak_effect"

  res_df$effect_class <- classify_effect_strength(
    res_df$padj_strong_effect,
    res_df$padj_weak_effect,
    alpha = alpha_level
  )

  res_df$standard_significant <- !is.na(res_df$padj_strong_effect) &
    (res_df$padj_strong_effect < alpha_level) &
    res_df$shrunk_lfc_pass

  res_df$deseq2_strong_call <- res_df$standard_significant
  res_df$deseq2_weak_call <- !is.na(res_df$padj_weak_effect) &
    (res_df$padj_weak_effect < alpha_level) &
    !res_df$shrunk_lfc_pass

  res_df$HBFSS_only_call <- res_df$HBFSS_significant & !res_df$standard_significant
  res_df$overlap_call <- res_df$HBFSS_significant & res_df$standard_significant

  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts$feature_id <- as.character(rownames(norm_counts))

  mm <- as.data.frame(mcols(dds))
  mm$feature_id <- as.character(rownames(mm))

  disp_cols_available <- intersect(
    c("feature_id", "dispGeneEst", "dispFit", "dispersion", "dispIter", "baseMean", "dispOutlier"),
    colnames(mm)
  )
  disp_df <- mm[, disp_cols_available, drop = FALSE]

  if ("baseMean" %in% colnames(disp_df)) {
    disp_df <- disp_df[, setdiff(colnames(disp_df), "baseMean"), drop = FALSE]
  }

  annot_df$feature_id <- as.character(annot_df$feature_id)
  annot_df$gene_symbol <- as.character(annot_df$gene_symbol)

  annot_df <- annot_df %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(is.na(gene_symbol), "", trimws(gene_symbol))) %>%
    dplyr::arrange(feature_id, dplyr::desc(gene_symbol != ""), gene_symbol) %>%
    dplyr::distinct(feature_id, .keep_all = TRUE) %>%
    dplyr::mutate(gene_symbol = dplyr::na_if(gene_symbol, ""))

  res_df$feature_id <- as.character(res_df$feature_id)
  norm_counts$feature_id <- as.character(norm_counts$feature_id)
  disp_df$feature_id <- as.character(disp_df$feature_id)

  norm_counts <- norm_counts[!duplicated(norm_counts$feature_id), , drop = FALSE]
  disp_df <- disp_df[!duplicated(disp_df$feature_id), , drop = FALSE]

  final_df <- res_df %>%
    dplyr::left_join(annot_df, by = "feature_id") %>%
    dplyr::left_join(norm_counts, by = "feature_id") %>%
    dplyr::left_join(disp_df, by = "feature_id")

  final_df$neglog10_padj <- safe_neglog10(final_df$padj)
  final_df$neglog10_empirical_p <- safe_neglog10(final_df$empirical_p)

  final_df$dataset_name <- dataset_name
  final_df$hc_p_threshold_dataset <- hc_p_threshold_dataset
  final_df$hbfss_threshold_dataset <- hbfss_threshold_dataset
  final_df <- reorder_result_columns(final_df)

  list(
    dds = dds,
    results = final_df,
    hc_p_threshold = hc_p_threshold_dataset,
    hbfss_threshold = hbfss_threshold_dataset
  )
}

# =============================================================================
# PLOTTING HELPERS
# =============================================================================

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

  pca_tbl$pc_index <- seq_len(nrow(pca_tbl))

  group_label_chr <- if (length(group_label) && !is.null(group_label[1])) {
    tolower(as.character(group_label[1]))
  } else {
    NA_character_
  }

  bar_fill <- if (!is.na(group_label_chr) && group_label_chr %in% c("control", "untrt")) {
    plot_palette$control
  } else {
    plot_palette$treatment
  }

  ggplot(pca_tbl, aes(pc_index, proportion_variance)) +
    geom_col(fill = bar_fill, color = "white") +
    geom_line(
      aes(y = cumulative_proportion, group = 1),
      linewidth = LINE_WIDTH_BOUNDARY,
      colour = plot_palette$threshold
    ) +
    geom_point(
      aes(y = cumulative_proportion),
      size = 1.6,
      colour = plot_palette$threshold
    ) +
    scale_x_continuous(
      breaks = pca_tbl$pc_index,
      labels = pca_tbl$principal_component
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      title = compact_title(
        paste(dataset_name, "|", pretty_group_label(group_label), "| PCA variance profile"),
        width = 42
      ),
      subtitle = compact_caption(
        paste(preprocessing_label, "|", final_analysis_label()),
        width = 72
      ),
      x = "Principal component",
      y = "Variance explained"
    ) +
    manuscript_theme() +
    theme(
      plot.margin = margin(t = 12, r = 14, b = 14, l = 14),
      legend.position = "none"
    )
}

plot_pca_scatter <- function(pca_fit, dataset_label, group_label, preprocessing_label = NULL) {
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

  subtitle_txt <- if (!is.null(preprocessing_label)) {
    compact_caption(paste(preprocessing_label, "|", final_analysis_label()), width = 72)
  } else {
    NULL
  }

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
      title = paste(dataset_label, "|", pretty_group_label(group_label), "| PCA scatter"),
      subtitle = subtitle_txt,
      x = paste0("PC1 (", pca_var_per[1], "%)"),
      y = paste0("PC2 (", pca_var_per[2], "%)")
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(
      plot.title  = element_text(margin = margin(b = 4)),
      plot.margin = margin(t = 12, r = 18, b = 16, l = 16)
    )
}

plot_pc1_loading_rank <- function(loading_tbl, cutoff, dataset_label, group_label,
                                  top_n_used = NA_integer_, cutoff_quantile = NA_real_,
                                  preprocessing_label = "Normalized before EVS",
                                  value_col = "pc1_loading_abs") {
  p <- ggplot(loading_tbl, aes(rank, .data[[value_col]])) +
    geom_line(linewidth = LINE_WIDTH_BOUNDARY, color = "grey35") +
    labs(
      title = compact_title(
        paste(dataset_label, "|", pretty_group_label(group_label), "| PC1 loading rank"),
        width = 36
      ),
      subtitle = compact_caption(
        paste0(
          preprocessing_label,
          ". Ranked loadings. ", final_analysis_label(), "."
        ),
        width = 72
      ),
      x = "Ranked PAS feature",
      y = "Loading magnitude"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()

  if (is.finite(cutoff) && !is.na(cutoff)) {
    quantile_label <- if (is.finite(cutoff_quantile)) {
      paste0("Upper-tail quantile = ", signif(cutoff_quantile, 4))
    } else {
      "Upper-tail quantile = NA"
    }

    p <- p +
      geom_hline(yintercept = cutoff, color = plot_palette$threshold, linewidth = LINE_WIDTH_THRESH) +
      annotate(
        "label",
        x     = max(loading_tbl$rank) * 0.70,
        y     = cutoff,
        label = paste0(
          "Projected cutoff = ", signif(cutoff, 4), "\n",
          "Approx N = ", ifelse(is.finite(top_n_used), top_n_used, NA_integer_), "\n",
          quantile_label
        ),
        fill  = "white",
        color = plot_palette$threshold,
        vjust = -0.7,
        size  = 3.2,
        label.size = 0.15
      )
  }

  p
}

plot_pc1_loading_histogram <- function(loading_tbl, cutoff, dataset_label, group_label,
                                       preprocessing_label = "Normalized before EVS",
                                       value_col = "pc1_loading_abs") {
  p <- ggplot(loading_tbl, aes(.data[[value_col]])) +
    geom_histogram(bins = 60, fill = "grey70", color = "white") +
    labs(
      title = compact_title(
        paste(dataset_label, "|", pretty_group_label(group_label), "| PC1 loading distribution"),
        width = 42
      ),
      subtitle = compact_caption(
        paste(
          preprocessing_label,
          "|", final_analysis_label()
        ),
        width = 78
      ),
      x = "Loading magnitude",
      y = "Feature count"
    ) +
    manuscript_theme()

  if (is.finite(cutoff) && !is.na(cutoff)) {
    p <- p + geom_vline(
      xintercept = cutoff,
      colour = plot_palette$threshold,
      linewidth = 0.9,
      linetype = "dashed"
    )
  }

  p
}

plot_nb_regime_rank_profile <- function(rank_tbl, cutoff_rank, dataset_label, group_label) {
  if (is.null(rank_tbl) || !nrow(rank_tbl)) return(NULL)
  cutoff_rank <- make_safe_rank(cutoff_rank, nrow(rank_tbl))
  if (!is.finite(cutoff_rank) || is.na(cutoff_rank)) return(NULL)

  plot_df <- data.frame(
    filtered_rank_index = rank_tbl$filtered_rank_index,
    log_IOD = rank_tbl$log_iod_nb_smooth,
    log_CV2 = rank_tbl$log_cv2_nb_smooth,
    stringsAsFactors = FALSE
  )

  iod_df <- data.frame(
    filtered_rank_index = plot_df$filtered_rank_index,
    value = plot_df$log_IOD,
    regime = "Leading-edge IOD",
    stringsAsFactors = FALSE
  )

  cv2_df <- data.frame(
    filtered_rank_index = plot_df$filtered_rank_index,
    value = plot_df$log_CV2,
    regime = "Remainder CV²",
    stringsAsFactors = FALSE
  )

  long_df <- rbind(iod_df, cv2_df)

  ggplot(long_df, aes(filtered_rank_index, value, color = regime)) +
    geom_line(linewidth = 0.65) +
    geom_vline(xintercept = cutoff_rank, linetype = "dashed", linewidth = 0.8, colour = plot_palette$threshold) +
    annotate(
      "label",
      x = cutoff_rank,
      y = max(long_df$value, na.rm = TRUE),
      label = paste0("Selected rank = ", cutoff_rank),
      vjust = -0.4,
      size = 3.0,
      fill = "white",
      colour = plot_palette$threshold,
      label.size = 0.15
    ) +
    scale_color_manual(
      values = c(
        "Leading-edge IOD" = plot_palette$iod,
        "Remainder CV²" = plot_palette$cv2
      ),
      name = "NB regime"
    ) +
    labs(
      title = compact_title(paste(dataset_label, "|", pretty_group_label(group_label), "| NB-regime ranked profile"), width = 42),
      subtitle = compact_caption(
        paste("Full-dataset DESeq2 dispersion structure mapped onto the combined EVS rank order.", final_analysis_label()),
        width = 72
      ),
      x = "Filtered ranked PAS feature",
      y = "Smoothed log-scale regime metric"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()
}

plot_nb_regime_scatter_with_segment_fits <- function(rank_tbl, cutoff_rank, dataset_label, group_label) {
  if (is.null(rank_tbl) || !nrow(rank_tbl)) return(NULL)
  cutoff_rank <- make_safe_rank(cutoff_rank, nrow(rank_tbl))
  if (!is.finite(cutoff_rank) || is.na(cutoff_rank) || cutoff_rank <= 0L || cutoff_rank >= nrow(rank_tbl)) return(NULL)

  lead_df <- rank_tbl[seq_len(cutoff_rank), , drop = FALSE]
  rem_df  <- rank_tbl[(cutoff_rank + 1L):nrow(rank_tbl), , drop = FALSE]

  fit_lead <- fit_free_slope_segment(lead_df$log_baseMean_smooth, lead_df$log_iod_nb_smooth)
  fit_rem  <- fit_free_slope_segment(rem_df$log_baseMean_smooth, rem_df$log_cv2_nb_smooth)

  lead_line <- data.frame(
    x = range(lead_df$log_baseMean_smooth, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  lead_line$y <- fit_lead$intercept + fit_lead$slope * lead_line$x
  lead_line$regime <- "Leading-edge IOD"

  rem_line <- data.frame(
    x = range(rem_df$log_baseMean_smooth, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  rem_line$y <- fit_rem$intercept + fit_rem$slope * rem_line$x
  rem_line$regime <- "Remainder CV²"

  point_df <- rbind(
    data.frame(
      x = lead_df$log_baseMean_smooth,
      y = lead_df$log_iod_nb_smooth,
      regime = "Leading-edge IOD",
      stringsAsFactors = FALSE
    ),
    data.frame(
      x = rem_df$log_baseMean_smooth,
      y = rem_df$log_cv2_nb_smooth,
      regime = "Remainder CV²",
      stringsAsFactors = FALSE
    )
  )

  line_df <- rbind(lead_line, rem_line)

  ggplot(point_df, aes(x, y, color = regime)) +
    geom_point(alpha = 0.55, size = 1.2) +
    geom_line(data = line_df, aes(x, y, color = regime), linewidth = 0.9, inherit.aes = FALSE) +
    scale_color_manual(
      values = c(
        "Leading-edge IOD" = plot_palette$iod,
        "Remainder CV²" = plot_palette$cv2
      ),
      name = "NB regime"
    ) +
    labs(
      title = compact_title(paste(dataset_label, "|", pretty_group_label(group_label), "| segment-fit view"), width = 42),
      subtitle = paste0("Leading-edge ranked top ", cutoff_rank, " features scored on IOD; remainder scored on CV²."),
      x = "Smoothed log(baseMean)",
      y = "Smoothed log-regime metric"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()
}

plot_combined_pca_scatter_panel <- function(evs, comparison_name) {
  p1 <- plot_pca_scatter(evs$fit_trt$pca_fit, comparison_name, "treatment", evs$fit_trt$preprocessing_label)
  p2 <- plot_pca_scatter(evs$fit_untrt$pca_fit, comparison_name, "control", evs$fit_untrt$preprocessing_label)
  p3 <- plot_pca_scatter(evs$fit_trt_raw_view$pca_fit, comparison_name, "treatment", evs$fit_trt_raw_view$preprocessing_label)
  p4 <- plot_pca_scatter(evs$fit_untrt_raw_view$pca_fit, comparison_name, "control", evs$fit_untrt_raw_view$preprocessing_label)

  arrangeGrob(
    p1, p2, p3, p4,
    ncol = 2,
    top = textGrob(
      paste0(comparison_name, " | EVS PCA"),
      gp = gpar(fontface = "bold", cex = 1.02)
    ),
    bottom = textGrob(
      paste0(
        "Top: ", if (normalize_before_evs) "Normalized before EVS." else "Raw counts before EVS.",
        " Bottom: Raw-count EVS view. Left: treatment. Right: control."
      ),
      gp = gpar(cex = 0.86)
    )
  )
}

plot_combined_pca_variance_panel <- function(evs, comparison_name) {
  p1 <- plot_pca_variance_profile(evs$fit_trt$pca_fit, comparison_name, evs$fit_trt$preprocessing_label, "treatment")
  p2 <- plot_pca_variance_profile(evs$fit_untrt$pca_fit, comparison_name, evs$fit_untrt$preprocessing_label, "control")
  p3 <- plot_pca_variance_profile(evs$fit_trt_raw_view$pca_fit, comparison_name, evs$fit_trt_raw_view$preprocessing_label, "treatment")
  p4 <- plot_pca_variance_profile(evs$fit_untrt_raw_view$pca_fit, comparison_name, evs$fit_untrt_raw_view$preprocessing_label, "control")

  arrangeGrob(
    p1, p2, p3, p4,
    ncol = 2,
    top = textGrob(
      paste0(comparison_name, " | EVS PCA variance profiles"),
      gp = gpar(fontface = "bold", cex = 1.02)
    ),
    bottom = textGrob(
      paste0(
        "Top: ", if (normalize_before_evs) "Normalized before EVS." else "Raw counts before EVS.",
        " Bottom: Raw-count EVS view. Left: treatment. Right: control."
      ),
      gp = gpar(cex = 0.86)
    )
  )
}

plot_combined_pc1_histogram_panel <- function(evs, comparison_name) {
  p1 <- plot_pc1_loading_histogram(
    evs$fit_trt$loading_table,
    evs$fit_trt$cutoff,
    comparison_name,
    "treatment",
    evs$fit_trt$preprocessing_label,
    value_col = "pc1_loading_abs"
  )

  p2 <- plot_pc1_loading_histogram(
    evs$fit_untrt$loading_table,
    evs$fit_untrt$cutoff,
    comparison_name,
    "control",
    evs$fit_untrt$preprocessing_label,
    value_col = "pc1_loading_abs"
  )

  p3 <- plot_pc1_loading_histogram(
    evs$fit_trt_raw_view$loading_table,
    evs$fit_trt_raw_view$cutoff,
    comparison_name,
    "treatment",
    evs$fit_trt_raw_view$preprocessing_label,
    value_col = "pc1_loading_abs"
  )

  p4 <- plot_pc1_loading_histogram(
    evs$fit_untrt_raw_view$loading_table,
    evs$fit_untrt_raw_view$cutoff,
    comparison_name,
    "control",
    evs$fit_untrt_raw_view$preprocessing_label,
    value_col = "pc1_loading_abs"
  )

  arrangeGrob(
    p1, p2, p3, p4,
    ncol = 2,
    top = textGrob(
      paste0(comparison_name, " | EVS PC1 loading histograms"),
      gp = gpar(fontface = "bold", cex = 1.02)
    ),
    bottom = textGrob(
      paste0(
        "Top: ", ranking_space_label(), ". Bottom: raw-count EVS view. ",
        "Left: treatment. Right: control."
      ),
      gp = gpar(cex = 0.86)
    )
  )
}

plot_combined_rank_panel <- function(evs, comparison_name) {
  p <- plot_pc1_loading_rank(
    evs$combined_fit$loading_table,
    evs$combined_fit$cutoff,
    comparison_name,
    "combined",
    top_n_used = evs$combined_fit$top_n_used,
    cutoff_quantile = evs$combined_fit$cutoff_quantile,
    preprocessing_label = evs$combined_fit$preprocessing_label,
    value_col = "combined_pc1_loading_abs"
  )

  arrangeGrob(
    p,
    ncol = 1,
    top = textGrob(
      paste0(comparison_name, " | Comparison-specific combined rank"),
      gp = gpar(fontface = "bold", cex = 1.02)
    ),
    bottom = textGrob(
      paste0("Combined rule: ", evs$combined_fit$combined_rule, " over |PC1 loading| treatment and control."),
      gp = gpar(cex = 0.86)
    )
  )
}

method_call_colors <- c(
  "Neither"            = "grey80",
  "DESeq2 weak only"   = "#4A90E2",
  "DESeq2 strong only" = "#C0392B",
  "HBFSS only"         = "#E67E22",
  "Overlap"            = "#7D3C98"
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
  "Neither"            = "grey45",
  "DESeq2 weak only"   = "black",
  "DESeq2 strong only" = "black",
  "HBFSS only"         = "black",
  "Overlap"            = "black"
)

build_reviewer_volcano_classes <- function(df) {
  df <- as.data.frame(df)

  df$method_call_class <- dplyr::case_when(
    !is.na(df$overlap_call)       & df$overlap_call       ~ "Overlap",
    !is.na(df$deseq2_strong_call) & df$deseq2_strong_call ~ "DESeq2 strong only",
    !is.na(df$deseq2_weak_call)   & df$deseq2_weak_call   ~ "DESeq2 weak only",
    !is.na(df$HBFSS_only_call)    & df$HBFSS_only_call    ~ "HBFSS only",
    TRUE ~ "Neither"
  )

  df$method_call_class <- factor(
    df$method_call_class,
    levels = c("Neither", "DESeq2 weak only", "DESeq2 strong only", "HBFSS only", "Overlap")
  )

  df$has_valid_gene_symbol <- !is.na(df$gene_symbol) &
    grepl("^[A-Za-z0-9._-]+$", trimws(df$gene_symbol))
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
      title   = paste(pretty_dataset_label(dataset_name), "| Standard volcano"),
      subtitle = compact_caption(paste("Call categories reflect DESeq2 and HBFSS agreement patterns.", final_analysis_label()), width = 76),
      x       = "Shrunken log2 fold change (β̂shrunk)",
      y       = expression(-log[10](padj))
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

plot_hbfss_volcano_panel <- function(df, dataset_name) {
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
      title   = paste(pretty_dataset_label(dataset_name), "| HBFSS volcano"),
      subtitle = compact_caption("Dashed curve marks the HBFSS decision boundary implied by the empirical HC threshold.", width = 76),
      x       = "Shrunken log2 fold change (β̂shrunk)",
      y       = expression(-log[10](p[empirical]))
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    plot_expand_xy() +
    volcano_guides()

  if (nrow(lab_df) > 0) p <- p + volcano_label_layer(lab_df)
  p
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
      title    = compact_title(paste(pretty_dataset_label(dataset_name), "| Dispersion profile"), width = 42),
      subtitle = compact_caption(final_analysis_label(), width = 72),
      x        = "baseMean (log10 scale)",
      y        = "Final dispersion (log10 scale)"
    ) +
    manuscript_theme() +
    volcano_guides()
}

# =============================================================================
# EXPORT HELPERS
# =============================================================================

save_cross_dataset_comparison_panels <- function(comparison_name, analysis_results, cmp_dir) {
  if (!isTRUE(export_all_plots)) return(invisible(NULL))

  keys_present <- names(analysis_results)[vapply(
    analysis_results,
    function(x) {
      is.list(x) &&
        !is.null(x$results) &&
        is.data.frame(x$results) &&
        nrow(x$results) > 0 &&
        !is.null(x$summary) &&
        is.data.frame(x$summary) &&
        nrow(x$summary) > 0 &&
        "dataset_name" %in% names(x$summary)
    },
    logical(1)
  )]

  if (length(keys_present) == 0) return(invisible(NULL))

  make_panel <- function(grob_list, title_text, ncols = NULL) {
    grob_list <- Filter(Negate(is.null), grob_list)
    if (!length(grob_list)) return(NULL)
    if (is.null(ncols)) ncols <- length(grob_list)

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

  std_grobs <- lapply(keys_present, function(k) {
    safe_plot_build(
      plot_standard_volcano(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1]
      ),
      paste0(comparison_name, ": standard volcano: ", k)
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
    safe_plot_build(
      plot_hbfss_volcano_panel(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1]
      ),
      paste0(comparison_name, ": HBFSS volcano: ", k)
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

  disp_grobs <- lapply(keys_present, function(k) {
    safe_plot_build(
      plot_dispersion_panel_for_dataset(
        analysis_results[[k]]$results,
        analysis_results[[k]]$summary$dataset_name[1]
      ),
      paste0(comparison_name, ": dispersion panel: ", k)
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

  invisible(NULL)
}

build_pc1_feature_export <- function(analysis_results, annot_df, comparison_name, evs, tab_dir) {
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

  load_tbl_trt <- evs$fit_trt$loading_table[, c("feature_id", "pc1_loading"), drop = FALSE]
  names(load_tbl_trt)[2] <- "pc1_loading_treatment_ranking_space"

  load_tbl_ctrl <- evs$fit_untrt$loading_table[, c("feature_id", "pc1_loading"), drop = FALSE]
  names(load_tbl_ctrl)[2] <- "pc1_loading_control_ranking_space"

  load_tbl_trt_raw <- evs$fit_trt_raw_view$loading_table[, c("feature_id", "pc1_loading"), drop = FALSE]
  names(load_tbl_trt_raw)[2] <- "pc1_loading_treatment_raw_count_view"

  load_tbl_ctrl_raw <- evs$fit_untrt_raw_view$loading_table[, c("feature_id", "pc1_loading"), drop = FALSE]
  names(load_tbl_ctrl_raw)[2] <- "pc1_loading_control_raw_count_view"

  load_tbl_combined <- evs$combined_fit$loading_table[, c(
    "feature_id", "pc1_loading_abs_trt", "pc1_loading_abs_ctrl", "combined_pc1_loading_abs"
  ), drop = FALSE]

  out <- dplyr::left_join(out, load_tbl_trt, by = "feature_id")
  out <- dplyr::left_join(out, load_tbl_ctrl, by = "feature_id")
  out <- dplyr::left_join(out, load_tbl_trt_raw, by = "feature_id")
  out <- dplyr::left_join(out, load_tbl_ctrl_raw, by = "feature_id")
  out <- dplyr::left_join(out, load_tbl_combined, by = "feature_id")

  out$comparison_name <- comparison_name
  out$ranking_space <- paste0(ranking_space_label(), " | comparison-specific combined rank")
  out$final_analysis_space <- "Raw counts"
  out$combined_loading_rule <- evs$combined_fit$combined_rule
  save_csv(out, file.path(tab_dir, paste0(comparison_name, "_PC1_loadings_baseMean_dispersion_export.csv")))
  out
}

plot_global_cutoff_audit <- function(global_evs_selection) {
  if (is.null(global_evs_selection$comparison_summary_table) ||
      !nrow(global_evs_selection$comparison_summary_table)) {
    return(NULL)
  }

  df <- global_evs_selection$comparison_summary_table
  df <- df[df$selected_rank_order == 1, , drop = FALSE]
  df$comparison_name <- factor(df$comparison_name, levels = unique(df$comparison_name))

  ggplot(df, aes(comparison_name, independently_best_rank)) +
    geom_point(size = 3.0, shape = 21, fill = plot_palette$overlap, color = "black") +
    geom_line(aes(group = 1), linewidth = 0.5, color = "grey40") +
    labs(
      title = "EVS selected cutoff ranks by comparison",
      subtitle = "Comparison-specific combined rank axis with fallback ladder (BIC -> AIC -> RSS -> deterministic fallback).",
      x = "Comparison",
      y = "Selected cutoff rank"
    ) +
    manuscript_theme()
}

# =============================================================================
# FULL COMPARISON PIPELINE
# =============================================================================

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df,
                                         precomputed_selection = NULL) {
  cmp_dir <- file.path(output_dir, comparison_name)
  tab_dir <- file.path(cmp_dir, "tables")
  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

  final_rank_to_use <- choose_final_rank_for_comparison(
    comparison_name = comparison_name,
    precomputed_selection = precomputed_selection
  )

  evs <- build_eigenvector_split(
    count_matrix = count_matrix,
    coldata = coldata,
    comparison_name = comparison_name,
    final_shared_rank = final_rank_to_use,
    precomputed_selection = precomputed_selection
  )

  message(sprintf("[%s] final rank mode: %s", comparison_name, evs_final_rank_mode))
  message(sprintf("[%s] final rank used: %d", comparison_name, final_rank_to_use))
  message(sprintf("[%s] raw rows: %d", comparison_name, nrow(evs$raw_dataset)))
  message(sprintf("[%s] leading-edge rows: %d", comparison_name, nrow(evs$leading_edge_dataset)))
  message(sprintf("[%s] remainder rows: %d", comparison_name, nrow(evs$remainder_dataset)))

  evs_cutoff_summary_out <- evs$evs_cutoff_summary
  evs_cutoff_summary_out$comparison_name <- comparison_name
  evs_cutoff_summary_out$final_rank_mode <- evs_final_rank_mode
  evs_cutoff_summary_out$final_rank_used <- final_rank_to_use
  save_csv(evs_cutoff_summary_out, file.path(tab_dir, paste0(comparison_name, "_EVS_cutoff_summary.csv")))

  if (!is.null(evs$final_rank_aggregation) && nrow(evs$final_rank_aggregation)) {
    save_csv(
      evs$final_rank_aggregation,
      file.path(tab_dir, paste0(comparison_name, "_EVS_final_rank_aggregation_summary.csv"))
    )
  }

  if (is.null(precomputed_selection$per_comparison[[comparison_name]]) ||
      is.null(precomputed_selection$per_comparison[[comparison_name]]$comparison_eval)) {
    stop(sprintf("[%s] comparison_eval was missing from precomputed_selection.", comparison_name))
  }
  cmp_eval <- precomputed_selection$per_comparison[[comparison_name]]$comparison_eval

  if (!is.null(cmp_eval$candidate_grid) && nrow(cmp_eval$candidate_grid) > 0) {
    save_csv(
      cmp_eval$candidate_grid,
      file.path(tab_dir, paste0(comparison_name, "_EVS_changepoint_candidate_grid.csv"))
    )
  }

  if (!is.null(cmp_eval$mode_diagnostics) && nrow(cmp_eval$mode_diagnostics) > 0) {
    save_csv(
      cmp_eval$mode_diagnostics,
      file.path(tab_dir, paste0(comparison_name, "_EVS_mode_diagnostics.csv"))
    )
  }

  combined_nb_tbl <- precomputed_selection$per_comparison[[comparison_name]]$combined_rank_nb_table
  full_disp_tbl <- precomputed_selection$per_comparison[[comparison_name]]$full_dispersion_table

  if (!is.null(combined_nb_tbl) && nrow(combined_nb_tbl) > 0) {
    save_csv(combined_nb_tbl, file.path(tab_dir, paste0(comparison_name, "_combined_ranked_nb_regime_table.csv")))
  }
  if (!is.null(full_disp_tbl) && nrow(full_disp_tbl) > 0) {
    save_csv(full_disp_tbl, file.path(tab_dir, paste0(comparison_name, "_full_dataset_deseq2_dispersion_table.csv")))
  }

  if (isTRUE(export_all_plots)) {
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
                              preprocessing_label = evs$fit_trt$preprocessing_label,
                              value_col = "pc1_loading_abs"),
        plot_pc1_loading_rank(evs$fit_untrt$loading_table, evs$fit_untrt$cutoff, comparison_name, "control",
                              top_n_used = evs$fit_untrt$top_n_used, cutoff_quantile = evs$fit_untrt$cutoff_quantile,
                              preprocessing_label = evs$fit_untrt$preprocessing_label,
                              value_col = "pc1_loading_abs"),
        plot_pc1_loading_rank(evs$fit_trt_raw_view$loading_table, evs$fit_trt_raw_view$cutoff, comparison_name, "treatment",
                              top_n_used = evs$fit_trt_raw_view$top_n_used, cutoff_quantile = evs$fit_trt_raw_view$cutoff_quantile,
                              preprocessing_label = evs$fit_trt_raw_view$preprocessing_label,
                              value_col = "pc1_loading_abs"),
        plot_pc1_loading_rank(evs$fit_untrt_raw_view$loading_table, evs$fit_untrt_raw_view$cutoff, comparison_name, "control",
                              top_n_used = evs$fit_untrt_raw_view$top_n_used, cutoff_quantile = evs$fit_untrt_raw_view$cutoff_quantile,
                              preprocessing_label = evs$fit_untrt_raw_view$preprocessing_label,
                              value_col = "pc1_loading_abs"),
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

    hist_panel <- safe_plot_build(
      plot_combined_pc1_histogram_panel(evs, comparison_name),
      paste0(comparison_name, ": EVS PC1 histogram panel")
    )

    if (!is.null(hist_panel)) {
      save_grob(
        hist_panel,
        file.path(cmp_dir, "EVS_PC1_loading_histograms_combined.png"),
        width = 18,
        height = 13
      )
    }

    combined_rank_panel <- safe_plot_build(
      plot_combined_rank_panel(evs, comparison_name),
      paste0(comparison_name, ": combined rank panel")
    )

    if (!is.null(combined_rank_panel)) {
      save_grob(
        combined_rank_panel,
        file.path(cmp_dir, "EVS_comparison_specific_combined_rank.png"),
        width = 12,
        height = 8
      )
    }

    nb_rank_profile_panel <- safe_plot_build(
      arrangeGrob(
        plot_nb_regime_rank_profile(combined_nb_tbl, final_rank_to_use, comparison_name, "combined"),
        ncol = 1,
        top = textGrob(
          paste0(comparison_name, " | DESeq2-informed NB regime ranked profile (comparison-specific combined rank)"),
          gp = gpar(fontface = "bold", cex = 1.02)
        )
      ),
      paste0(comparison_name, ": NB regime ranked profile panel")
    )

    if (!is.null(nb_rank_profile_panel)) {
      save_grob(
        nb_rank_profile_panel,
        file.path(cmp_dir, "EVS_NB_regime_rank_profiles_combined.png"),
        width = 12, height = 7
      )
    }

    nb_segment_fit_panel <- safe_plot_build(
      arrangeGrob(
        plot_nb_regime_scatter_with_segment_fits(combined_nb_tbl, final_rank_to_use, comparison_name, "combined"),
        ncol = 1,
        top = textGrob(
          paste0(comparison_name, " | DESeq2-informed IOD/CV² segment-fit diagnostics (comparison-specific combined rank)"),
          gp = gpar(fontface = "bold", cex = 1.02)
        )
      ),
      paste0(comparison_name, ": NB regime segment-fit panel")
    )

    if (!is.null(nb_segment_fit_panel)) {
      save_grob(
        nb_segment_fit_panel,
        file.path(cmp_dir, "EVS_NB_regime_segment_fit_combined.png"),
        width = 12, height = 7
      )
    }
  }

  dataset_list <- list(
    raw_dataset          = evs$raw_dataset,
    leading_edge_dataset = evs$leading_edge_dataset,
    remainder_dataset    = evs$remainder_dataset
  )

  analysis_results <- list()

  for (nm in names(dataset_list)) {
    dataset_mat <- as.matrix(dataset_list[[nm]])
    full_dataset_name <- paste(comparison_name, nm, sep = "_")

    message(sprintf("[%s][%s] rows=%d cols=%d", comparison_name, nm, nrow(dataset_mat), ncol(dataset_mat)))

    if (nrow(dataset_mat) < 1L) {
      warning(sprintf("[%s][%s] Skipping empty dataset.", comparison_name, nm))
      next
    }

    fit <- run_core_analysis(
      count_mat = dataset_mat,
      coldata = coldata,
      dataset_name = full_dataset_name,
      annot_df = annot_df
    )

    df <- fit$results

    save_csv(df, file.path(tab_dir, paste0(full_dataset_name, "_results_full.csv")))
    save_csv(subset(df, standard_significant), file.path(tab_dir, paste0(full_dataset_name, "_standard_significant.csv")))
    save_csv(subset(df, HBFSS_significant),    file.path(tab_dir, paste0(full_dataset_name, "_HBFSS_significant.csv")))
    save_csv(subset(df, effect_class == "strong_effect"), file.path(tab_dir, paste0(full_dataset_name, "_strong_effect.csv")))
    save_csv(subset(df, effect_class == "weak_effect"),   file.path(tab_dir, paste0(full_dataset_name, "_weak_effect.csv")))

    summary_row <- data.frame(
      comparison_name        = comparison_name,
      dataset_name           = full_dataset_name,
      dataset_type           = nm,
      ranking_space          = paste0(ranking_space_label(), " | comparison-specific combined rank"),
      final_analysis_space   = "Raw counts",
      n_features             = nrow(df),
      hc_p_threshold         = fit$hc_p_threshold,
      hbfss_threshold        = fit$hbfss_threshold,
      n_standard_significant = sum(df$standard_significant, na.rm = TRUE),
      n_HBFSS_significant    = sum(df$HBFSS_significant, na.rm = TRUE),
      n_overlap_significant  = sum(df$standard_significant & df$HBFSS_significant, na.rm = TRUE),
      n_strong_effect        = sum(df$effect_class == "strong_effect", na.rm = TRUE),
      n_weak_effect          = sum(df$effect_class == "weak_effect", na.rm = TRUE),
      evs_selection_space    = paste0(ranking_space_label(), " | comparison-specific combined rank"),
      evs_dispersion_source  = evs_dispersion_source,
      evs_final_rank_mode    = "comparison_specific_combined_rank",
      final_rank_used        = final_rank_to_use,
      combined_loading_rule  = combined_loading_rule,
      model_selection_mode_used = cmp_eval$best_mode,
      cutoff_method          = evs$combined_fit$cutoff_method,
      cutoff_quantile        = evs$combined_fit$cutoff_quantile,
      selected_reason        = evs$combined_fit$selected_reason,
      fallback_rank_mode     = fallback_rank_mode,
      stringsAsFactors       = FALSE
    )

    save_csv(summary_row, file.path(tab_dir, paste0(full_dataset_name, "_summary.csv")))

    analysis_results[[nm]] <- list(
      dds = fit$dds,
      results = df,
      summary = summary_row
    )
  }

  tryCatch(
    build_pc1_feature_export(
      analysis_results = analysis_results,
      annot_df = annot_df,
      comparison_name = comparison_name,
      evs = evs,
      tab_dir = tab_dir
    ),
    error = function(e) {
      warning(paste0(comparison_name, ": PC1 feature export failed: ", conditionMessage(e)))
      NULL
    }
  )

  if (isTRUE(export_all_plots)) {
    tryCatch(
      save_cross_dataset_comparison_panels(
        comparison_name = comparison_name,
        analysis_results = analysis_results,
        cmp_dir = cmp_dir
      ),
      error = function(e) {
        warning(paste0(comparison_name, ": cross-dataset panel generation failed: ", conditionMessage(e)))
      }
    )
  }

  sm_list <- lapply(analysis_results, `[[`, "summary")
  sm_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, sm_list)

  if (!length(sm_list)) return(data.frame())
  dplyr::bind_rows(sm_list)
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

global_evs_selection <- precompute_global_evs_rank_selection(comparison_inputs)

if (!is.null(global_evs_selection$comparison_summary_table) &&
    nrow(global_evs_selection$comparison_summary_table) > 0) {
  save_csv(
    global_evs_selection$comparison_summary_table,
    file.path(output_dir, "EVS_candidate_cutoffs_by_comparison.csv")
  )
}

if (!is.null(global_evs_selection$mode_diagnostics_table) &&
    nrow(global_evs_selection$mode_diagnostics_table) > 0) {
  save_csv(
    global_evs_selection$mode_diagnostics_table,
    file.path(output_dir, "EVS_mode_diagnostics_all_comparisons.csv")
  )
}

if (isTRUE(export_all_plots)) {
  cutoff_audit_plot <- safe_plot_build(
    plot_global_cutoff_audit(global_evs_selection),
    "Global cutoff audit plot"
  )

  if (!is.null(cutoff_audit_plot)) {
    save_grob(
      cutoff_audit_plot,
      file.path(output_dir, "EVS_cutoff_audit_plot.png"),
      width = 10,
      height = 6
    )
  }
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
      count_matrix = input_obj$count_matrix,
      coldata = input_obj$coldata,
      annot_df = OrigID_Symbol,
      precomputed_selection = global_evs_selection
    ),
    error = function(e) e
  )

  if (inherits(out, "error")) {
    failed_comparisons[[cmp]] <- data.frame(
      comparison_name = cmp,
      error_message = conditionMessage(out),
      stringsAsFactors = FALSE
    )
  } else if (!is.null(out) && nrow(out) > 0) {
    all_summaries_list[[cmp]] <- out
  }
}

all_summaries <- if (length(all_summaries_list) > 0) {
  dplyr::bind_rows(all_summaries_list)
} else {
  data.frame()
}

if (nrow(all_summaries) > 0) {
  save_csv(all_summaries, file.path(output_dir, "EVS_overall_pipeline_summary.csv"))
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

# -----------------------------------------------------------------------------
# AUTO-PUSH EXPORTS TO GITHUB
# -----------------------------------------------------------------------------

cmd <- paste(
  "cd", shQuote(repo_dir), "&&",
  "git add -A exports &&",
  "branch=$(git rev-parse --abbrev-ref HEAD) &&",
  "if ! git diff --cached --quiet; then",
  "git commit -m", shQuote("Auto-update pipeline outputs"), "&&",
  "git push origin \"$branch\";",
  "fi"
)

if (isTRUE(auto_push_exports)) {
  status <- system(cmd)

  if (status == 0) {
    message("Git push successful.")
  } else {
    warning("Git push failed. Check authentication.")
  }
} else {
  message("auto_push_exports = FALSE; skipping git commit/push.")
}
