# =============================================================================
# SEQUENCE MANUSCRIPT PIPELINE (REPAIRED)
# =============================================================================
# This script is a repaired, runnable manuscript pipeline for the SEQUENCE WTTS-
# Seq dataset. It preserves the core workflow you were using:
#   1. Build treatment and control EVS loading tables.
#   2. Estimate ranked IOD and CV2 structure within each condition.
#   3. Select one shared EVS cutoff per comparison from the first stable
#      crossing of the combined treatment plus control local IOD and local CV2
#      curves.
#   4. Split each comparison into original, leading-edge, and remainder sets.
#   5. Run DESeq2, empirical-null fitting, HC thresholding, and HBFSS.
#   6. Export figures and tables into one stable repository folder.
#
# Important repaired behavior:
#   - The EVS shared cutoff is based on the first stable IOD minus CV2 crossing.
#   - Weak DC2 calls are gated by the HC threshold computed on Wald p-values.
#   - HBFSS uses empirical p-values from fdrtool and requires both:
#         HBFSS >= HBFSS threshold
#         empirical_p <= empirical HC threshold
#   - Volcano markers and the legend use consistent colored markers.
#   - Cross-comparison outputs are written into one stable exports folder.
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
analysis_name <- "EVS_HBFSS_AllComparisons_Output"
output_dir <- file.path(output_root, analysis_name)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
run_twas_overlap <- FALSE
twas_file <- file.path(input_dir, "3aTWAS_genes_of_11_brain_disorders.csv")

alpha_level <- 0.10
lfc_boundary <- 1.0
max_usable_hc_p_threshold <- 0.95
fdrtool_pct0 <- 0.75
fdr_clip_floor <- 1e-300
fdr_clip_ceiling <- 0.99

# Shared EVS cutoff settings
crossing_min_percentile <- 0.05
crossing_max_percentile <- 0.95
crossing_stability_window_n <- 3L
crossing_rule <- "first_left_crossing"

fourier_percentile_step <- 0.01
fourier_window_fraction <- 0.12
fourier_harmonics <- 2L
fourier_min_rank <- 1L
fourier_max_rank_frac <- 1.00
fourier_score_weight_iod_amp <- 1.0
fourier_score_weight_cv2_amp <- 1.0
fourier_score_weight_center_agreement <- 1.0

evs_fixed_top_n <- 5000L

figure_dpi <- 320
base_theme_size <- 10

POINT_SIZE_PRIMARY <- 1.6
POINT_SIZE_DISP <- 1.3
POINT_ALPHA_PRIMARY <- 0.82
POINT_ALPHA_DISP <- 0.55
POINT_STROKE <- 0.40
LINE_WIDTH_BOUNDARY <- 0.55
LINE_WIDTH_ZERO <- 0.40
LINE_WIDTH_THRESH <- 0.90
HIST_BINS <- 60
HIST_COLOR <- "white"

plot_palette <- list(
  background = "#BDBDBD",
  threshold = "#8C2D04",
  hbfss = "#E67E22",
  deseq2 = "#C0392B",
  strong = "#C0392B",
  overlap = "#7D3C98",
  weak = "#4A90E2",
  intermediate = "#7F7F7F",
  histogram = "#969696",
  control = "#4D4D4D",
  treatment = "#1F78B4"
)

export_optional_mean_histograms <- FALSE
export_optional_empirical_p_histograms <- FALSE
export_optional_hbfss_distributions <- FALSE
export_optional_evs_raw_histograms <- FALSE
export_optional_evs_variance_profiles <- TRUE

figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# EMBEDDED METADATA
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

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix = c("R0", "R2", "R4", "R8"),
  group2_prefix = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

dataset_key_order <- c("raw_dataset", "leading_edge_dataset", "remainder_dataset")
dataset_key_labels <- c(
  raw_dataset = "Original dataset",
  leading_edge_dataset = "Leading-edge dataset",
  remainder_dataset = "Remainder dataset"
)

# -----------------------------------------------------------------------------
# HELPERS
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

clip_probabilities <- function(x, eps = fdr_clip_floor, upper = fdr_clip_ceiling) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  x <- pmax(x, eps, na.rm = FALSE)
  x <- pmin(x, upper, na.rm = FALSE)
  x
}

compact_title <- function(x, width = 58) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

compact_caption <- function(x, width = 120) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

strip_dataset_key <- function(dataset_name) {
  sub("^.*_(raw_dataset|leading_edge_dataset|remainder_dataset)$", "\\1", dataset_name)
}

comparison_label_from_dataset_name <- function(dataset_name) {
  sub("_(raw_dataset|leading_edge_dataset|remainder_dataset)$", "", dataset_name)
}

pretty_dataset_type <- function(dataset_key) {
  switch(
    dataset_key,
    raw_dataset = "Original dataset",
    leading_edge_dataset = "Leading-edge dataset",
    remainder_dataset = "Remainder dataset",
    dataset_key
  )
}

pretty_dataset_label <- function(dataset_name) {
  if (!grepl("_(raw_dataset|leading_edge_dataset|remainder_dataset)$", dataset_name)) {
    return(dataset_name)
  }
  dataset_key <- strip_dataset_key(dataset_name)
  comparison_name <- comparison_label_from_dataset_name(dataset_name)
  paste(comparison_name, pretty_dataset_type(dataset_key), sep = " | ")
}

pretty_group_label <- function(group_label) {
  switch(
    as.character(group_label),
    trt = "Treatment",
    untrt = "Control",
    treatment = "Treatment",
    control = "Control",
    as.character(group_label)
  )
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

save_csv <- function(df, path) {
  utils::write.csv(df, file = path, row.names = FALSE)
}

save_grob <- function(g, path, width = 8.8, height = 6.5, dpi = figure_dpi, bg = "white") {
  if (is.null(g)) return(invisible(FALSE))
  ggplot2::ggsave(
    filename = path,
    plot = g,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = bg,
    limitsize = FALSE
  )
  invisible(TRUE)
}

safe_plot_build <- function(expr, label = "plot") {
  tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) {
      warning(sprintf("%s failed: %s", label, conditionMessage(e)), call. = FALSE)
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

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_theme_size + 1, hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = base_theme_size - 1, hjust = 0.5, margin = margin(b = 7)),
      plot.caption = element_text(size = base_theme_size - 3, hjust = 0.5, colour = "grey30", margin = margin(t = 8)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = base_theme_size - 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      plot.margin = margin(t = 18, r = 30, b = 22, l = 24)
    )
}

format_runtime_minutes <- function(seconds_value) {
  if (!is.finite(seconds_value) || is.na(seconds_value) || seconds_value < 0) return("unknown")
  if (seconds_value < 60) return(sprintf("~%d sec", as.integer(round(seconds_value))))
  sprintf("~%.1f min", seconds_value / 60)
}

make_design_formula <- function() {
  ~ condition
}

get_condition_coef <- function(dds) {
  rn <- resultsNames(dds)
  idx <- grep("^condition_", rn)
  if (length(idx) == 0) stop("Could not identify condition coefficient in resultsNames(dds).", call. = FALSE)
  rn[idx[1]]
}

clean_gene_set <- function(x) {
  unique(tolower(trimws(x[!is.na(x) & x != ""])))
}

# -----------------------------------------------------------------------------
# EMPIRICAL NULL HELPERS
# -----------------------------------------------------------------------------
run_empirical_null_fdrtool <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  stat_vec <- unname(stat_vec)

  if (length(stat_vec) < 50) {
    stop(sprintf("[%s] Too few finite Wald statistics for stable empirical-null fitting.", dataset_name), call. = FALSE)
  }

  stat_sd <- suppressWarnings(stats::sd(stat_vec, na.rm = TRUE))
  if (!is.finite(stat_sd) || is.na(stat_sd) || stat_sd == 0) {
    stop(sprintf("[%s] Wald statistics have zero or invalid spread.", dataset_name), call. = FALSE)
  }

  fit <- tryCatch(
    fdrtool(stat_vec, statistic = "normal", plot = FALSE, verbose = FALSE, cutoff.method = "fndr", pct0 = fdrtool_pct0),
    error = function(e1) {
      message(sprintf("[%s] Primary fdrtool call failed: %s", dataset_name, conditionMessage(e1)))
      tryCatch(
        fdrtool(stat_vec, statistic = "normal", plot = FALSE, verbose = FALSE, cutoff.method = "pct0", pct0 = fdrtool_pct0),
        error = function(e2) {
          stop(sprintf("[%s] fdrtool failed after retry: %s", dataset_name, conditionMessage(e2)), call. = FALSE)
        }
      )
    }
  )

  fit$pval <- clip_probabilities(fit$pval)
  fit$qval <- clip_probabilities(fit$qval)
  fit$lfdr <- as.numeric(fit$lfdr)
  fit
}

safe_hc_thresh <- function(p_vec, dataset_name, upper_cap = max_usable_hc_p_threshold) {
  sorted_p <- sort(clip_probabilities(p_vec), na.last = NA, decreasing = FALSE)
  if (length(sorted_p) < 5) return(NA_real_)

  out <- suppressWarnings(
    tryCatch(
      fdrtool::hc.thresh(as.vector(sorted_p)),
      error = function(e) {
        message(sprintf("[%s] hc.thresh failed: %s", dataset_name, conditionMessage(e)))
        NA_real_
      }
    )
  )

  out <- as.numeric(out[1])
  if (!is.finite(out) || out <= 0 || out >= upper_cap) return(NA_real_)
  out
}

resolve_hc_threshold <- function(p_vec, dataset_name, source_label) {
  primary_hc <- safe_hc_thresh(p_vec, dataset_name = paste(dataset_name, source_label, sep = " | "))
  if (is.finite(primary_hc) && !is.na(primary_hc)) {
    return(list(threshold = primary_hc, source = paste0(source_label, "_hc_primary")))
  }

  p_valid <- clip_probabilities(p_vec)
  p_valid <- p_valid[is.finite(p_valid) & !is.na(p_valid)]
  if (length(p_valid) >= 5) {
    fallback_hc <- suppressWarnings(stats::quantile(p_valid, probs = alpha_level, na.rm = TRUE, type = 7))
    fallback_hc <- as.numeric(fallback_hc[1])
    if (is.finite(fallback_hc) && !is.na(fallback_hc) && fallback_hc > 0 && fallback_hc < max_usable_hc_p_threshold) {
      return(list(threshold = fallback_hc, source = paste0(source_label, "_quantile_fallback")))
    }
  }

  list(threshold = NA_real_, source = paste0(source_label, "_missing"))
}

# -----------------------------------------------------------------------------
# IMPORT WTTS FILE
# -----------------------------------------------------------------------------
if (!file.exists(count_file)) {
  stop(sprintf("Count file not found: %s", count_file), call. = FALSE)
}

WTTS_Seq <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
WTTS_Seq <- as.data.frame(WTTS_Seq, stringsAsFactors = FALSE)
WTTS_Seq$OrigID <- as.character(WTTS_Seq$OrigID)
WTTS_Seq$Symbol <- as.character(WTTS_Seq$Symbol)

assert_required_columns(WTTS_Seq, c("OrigID", "Symbol"), object_name = "WTTS count file")
assert_required_columns(WTTS_Seq, meta_all$id, object_name = "WTTS count file sample columns")

WTTS_Seq <- WTTS_Seq[!is.na(WTTS_Seq$OrigID), , drop = FALSE]
sample_na <- rowSums(is.na(WTTS_Seq[, meta_all$id, drop = FALSE])) > 0
WTTS_Seq <- WTTS_Seq[!sample_na, , drop = FALSE]
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
# COMPARISON PREPARATION
# -----------------------------------------------------------------------------
prepare_comparison_data <- function(comparison_name, group1_prefix, group2_prefix, WTTS_Seq, meta_all) {
  keep_ids <- grepl(paste0("^", group1_prefix, "_"), meta_all$id) |
    grepl(paste0("^", group2_prefix, "_"), meta_all$id)

  meta_sub <- meta_all[keep_ids, , drop = FALSE]
  coldata <- meta_sub[, c("condition"), drop = FALSE]
  sample_ids <- rownames(meta_sub)

  missing_samples <- setdiff(sample_ids, colnames(WTTS_Seq))
  if (length(missing_samples) > 0) {
    stop(paste("Missing samples in WTTS file for", comparison_name, ":", paste(missing_samples, collapse = ", ")), call. = FALSE)
  }

  count_sub <- WTTS_Seq[, sample_ids, drop = FALSE]
  stopifnot(all(colnames(count_sub) == rownames(coldata)))

  count_mat <- as.matrix(count_sub)
  storage.mode(count_mat) <- "numeric"

  list(
    comparison_name = comparison_name,
    count_matrix = count_mat,
    coldata = coldata
  )
}

compute_condition_feature_metrics <- function(count_submatrix) {
  cd <- S4Vectors::DataFrame(row.names = colnames(count_submatrix))
  dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(count_submatrix)), colData = cd, design = ~ 1)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds, quiet = TRUE)
  dds <- estimateDispersionsFit(dds, quiet = TRUE)
  md <- as.data.frame(SummarizedExperiment::mcols(dds), stringsAsFactors = FALSE)
  md$feature_id <- rownames(md)
  keep_cols <- intersect(c("feature_id", "baseMean", "dispGeneEst", "dispFit", "dispersion"), colnames(md))
  md[, keep_cols, drop = FALSE]
}

compute_pc1_loading_table <- function(value_df, sample_names, preprocessing_label, feature_metrics = NULL) {
  x <- as.matrix(value_df[, sample_names, drop = FALSE])
  x_log <- log2(x + 1)
  pca_fit <- prcomp(t(x_log), scale. = FALSE, rank. = min(5L, ncol(x_log)))
  loading_abs <- abs(pca_fit$rotation[, 1])
  loading_tbl <- data.frame(
    feature_id = names(loading_abs),
    pc1_loading = unname(pca_fit$rotation[, 1]),
    pc1_loading_abs = unname(loading_abs),
    stringsAsFactors = FALSE
  )

  if (!is.null(feature_metrics) && nrow(feature_metrics) > 0) {
    fm <- as.data.frame(feature_metrics, stringsAsFactors = FALSE)
    keep_cols <- intersect(c("feature_id", "baseMean", "dispGeneEst", "dispFit", "dispersion"), colnames(fm))
    fm <- fm[, keep_cols, drop = FALSE]
    fm <- fm[!duplicated(fm$feature_id), , drop = FALSE]
    loading_tbl <- dplyr::left_join(loading_tbl, fm, by = "feature_id")
  }

  loading_tbl <- loading_tbl[order(loading_tbl$pc1_loading_abs, decreasing = TRUE), , drop = FALSE]
  loading_tbl$rank <- seq_len(nrow(loading_tbl))
  loading_tbl$split_class <- "background_loading"

  list(
    pca_fit = pca_fit,
    loading_table = loading_tbl,
    cutoff = NA_real_,
    top_n_used = NA_integer_,
    cutoff_quantile = NA_real_,
    cutoff_method = "shared_cutoff_pending",
    selected_reason = "shared_cutoff_pending",
    preprocessing_label = preprocessing_label
  )
}

# -----------------------------------------------------------------------------
# FOURIER SHARED-CUTOFF METHOD
# -----------------------------------------------------------------------------
build_ranked_fourier_metric_table <- function(loading_tbl, mean_col = "baseMean", dispersion_col = "dispGeneEst") {
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
  if (nrow(df) < 25L) return(NULL)

  df <- df[order(df$rank), , drop = FALSE]
  mu <- pmax(df$mean_value, 1e-12)
  alpha <- pmax(df$dispersion_value, 1e-12)
  df$iod_nb <- 1 + alpha * mu
  df$cv2_nb <- (1 / mu) + alpha
  df$log_iod_nb <- log10(df$iod_nb)
  df$log_cv2_nb <- log10(df$cv2_nb)
  df
}

build_percentile_windows_fourier <- function(n_total, step = fourier_percentile_step, window_fraction = fourier_window_fraction, min_rank = fourier_min_rank, max_rank_frac = fourier_max_rank_frac) {
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
  if (!is.finite(stats::sd(y_vec, na.rm = TRUE)) || stats::sd(y_vec, na.rm = TRUE) == 0) return(NULL)

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

build_local_fourier_wave_map <- function(loading_tbl, mean_col = "baseMean", dispersion_col = "dispGeneEst") {
  metric_df <- build_ranked_fourier_metric_table(loading_tbl, mean_col, dispersion_col)
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

  list(metric_df = metric_df, wave_map = wave_map)
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
  required_cols <- c("percentile", "combined_center_rank", "combined_iod_amplitude", "combined_cv2_amplitude")
  assert_required_columns(df, required_cols, object_name = "combined_wave_df")
  df$regime_difference <- df$combined_iod_amplitude - df$combined_cv2_amplitude
  df$regime_direction <- ifelse(
    is.na(df$regime_difference),
    NA_character_,
    ifelse(df$regime_difference > 0, "IOD_dominant", ifelse(df$regime_difference < 0, "CV2_dominant", "balanced"))
  )
  df
}

find_crossing_intervals <- function(diff_df, min_percentile = crossing_min_percentile, max_percentile = crossing_max_percentile) {
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

    if (identical(crossing_rule, "first_left_crossing")) {
      crossing_percentile <- x1
      crossing_rank <- r1
    } else if (isTRUE(all.equal(y1, y2))) {
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

label_stable_crossings <- function(diff_df, crossing_tbl, stability_window_n = crossing_stability_window_n) {
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
    left_pos_frac <- mean(df$regime_difference[left_idx] > 0, na.rm = TRUE)
    right_neg_frac <- mean(df$regime_difference[right_idx] < 0, na.rm = TRUE)
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
    return(list(diff_df = diff_df, crossing_table = crossing_tbl, selected_crossing = NULL, selected_reason = "no_crossings_found"))
  }

  stable_tbl <- crossing_tbl[crossing_tbl$stable_crossing, , drop = FALSE]
  if (nrow(stable_tbl)) {
    selected <- stable_tbl[order(stable_tbl$crossing_percentile), , drop = FALSE][1, , drop = FALSE]
    reason <- "first_stable_crossing"
  } else {
    selected <- crossing_tbl[order(crossing_tbl$crossing_percentile), , drop = FALSE][1, , drop = FALSE]
    reason <- "first_crossing_fallback"
  }

  list(diff_df = diff_df, crossing_table = crossing_tbl, selected_crossing = selected, selected_reason = reason)
}

build_shared_combined_loading_table <- function(fit_trt_loading_tbl, fit_ctrl_loading_tbl) {
  trt_tbl <- as.data.frame(fit_trt_loading_tbl, stringsAsFactors = FALSE)
  ctrl_tbl <- as.data.frame(fit_ctrl_loading_tbl, stringsAsFactors = FALSE)
  assert_required_columns(trt_tbl, c("feature_id", "pc1_loading", "pc1_loading_abs", "rank"), object_name = "fit_trt_loading_tbl")
  assert_required_columns(ctrl_tbl, c("feature_id", "pc1_loading", "pc1_loading_abs", "rank"), object_name = "fit_ctrl_loading_tbl")

  trt_use <- trt_tbl[, c("feature_id", "pc1_loading", "pc1_loading_abs", "rank"), drop = FALSE]
  ctrl_use <- ctrl_tbl[, c("feature_id", "pc1_loading", "pc1_loading_abs", "rank"), drop = FALSE]
  names(trt_use) <- c("feature_id", "pc1_loading_trt", "pc1_loading_abs_trt", "rank_trt")
  names(ctrl_use) <- c("feature_id", "pc1_loading_ctrl", "pc1_loading_abs_ctrl", "rank_ctrl")
  merged <- dplyr::full_join(trt_use, ctrl_use, by = "feature_id")
  merged$pc1_loading_abs_trt[is.na(merged$pc1_loading_abs_trt)] <- 0
  merged$pc1_loading_abs_ctrl[is.na(merged$pc1_loading_abs_ctrl)] <- 0
  merged$combined_loading <- pmax(merged$pc1_loading_abs_trt, merged$pc1_loading_abs_ctrl, na.rm = TRUE)
  merged <- merged[order(-merged$combined_loading, merged$feature_id), , drop = FALSE]
  merged$combined_rank <- seq_len(nrow(merged))
  rownames(merged) <- NULL
  merged
}

resolve_top_n_cutoff <- function(sorted_values_desc, top_n = evs_fixed_top_n) {
  n_total <- length(sorted_values_desc)
  if (n_total == 0) stop("resolve_top_n_cutoff() received an empty vector.", call. = FALSE)
  top_n_actual <- min(max(1L, as.integer(top_n)), n_total)
  cutoff_value <- sorted_values_desc[top_n_actual]
  cutoff_quantile <- 1 - (top_n_actual / n_total)
  list(
    top_n_actual = top_n_actual,
    cutoff_value = cutoff_value,
    cutoff_quantile = cutoff_quantile,
    n_total = n_total,
    method = "fixed_top_n",
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

resolve_combined_fourier_cutoff <- function(fit_trt_loading_tbl, fit_ctrl_loading_tbl, fixed_top_n = evs_fixed_top_n) {
  shared_combined_tbl <- build_shared_combined_loading_table(fit_trt_loading_tbl, fit_ctrl_loading_tbl)
  n_total <- nrow(shared_combined_tbl)
  fallback <- resolve_top_n_cutoff(shared_combined_tbl$combined_loading, top_n = fixed_top_n)

  trt_wave_obj <- build_local_fourier_wave_map(fit_trt_loading_tbl)
  ctrl_wave_obj <- build_local_fourier_wave_map(fit_ctrl_loading_tbl)

  if (is.null(trt_wave_obj) || is.null(ctrl_wave_obj)) {
    fallback$method <- "fixed_top_n_fallback"
    fallback$trt_wave_obj <- trt_wave_obj
    fallback$ctrl_wave_obj <- ctrl_wave_obj
    fallback$shared_combined_tbl <- shared_combined_tbl
    fallback$combined_wave_map <- data.frame()
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
    crossing_table = crossing_info$crossing_table,
    selected_crossing = selected_crossing,
    candidate_table = candidate_table,
    selected_reason = crossing_info$selected_reason
  )
}

select_wave_backup_cutoff <- function(wave_obj, group_label = "group") {
  if (is.null(wave_obj) || is.null(wave_obj$wave_map) || !nrow(wave_obj$wave_map)) {
    return(list(diff_df = data.frame(), crossing_table = data.frame(), selected_crossing = NULL, selected_reason = "wave_map_missing"))
  }
  df <- as.data.frame(wave_obj$wave_map, stringsAsFactors = FALSE)
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
    return(list(diff_df = diff_df, crossing_table = crossing_tbl, selected_crossing = selected, selected_reason = reason))
  }

  list(diff_df = diff_df, crossing_table = data.frame(), selected_crossing = NULL, selected_reason = paste0(group_label, "_no_crossing"))
}

build_group_backup_crossings <- function(combined_cutoff_info) {
  trt_backup <- select_wave_backup_cutoff(combined_cutoff_info$trt_wave_obj, group_label = "treatment")
  ctrl_backup <- select_wave_backup_cutoff(combined_cutoff_info$ctrl_wave_obj, group_label = "control")
  combined_cutoff_info$trt_backup_crossing <- trt_backup$selected_crossing
  combined_cutoff_info$trt_backup_reason <- trt_backup$selected_reason
  combined_cutoff_info$trt_backup_table <- trt_backup$crossing_table
  combined_cutoff_info$ctrl_backup_crossing <- ctrl_backup$selected_crossing
  combined_cutoff_info$ctrl_backup_reason <- ctrl_backup$selected_reason
  combined_cutoff_info$ctrl_backup_table <- ctrl_backup$crossing_table
  combined_cutoff_info
}

apply_loading_cutoff <- function(fit_obj, rank_index, selected_reason = "manual_override") {
  fit_obj <- as.list(fit_obj)
  loading_tbl <- as.data.frame(fit_obj$loading_table, stringsAsFactors = FALSE)
  assert_required_columns(loading_tbl, c("feature_id", "pc1_loading_abs", "rank"), object_name = "fit_obj$loading_table")
  loading_tbl <- loading_tbl[order(loading_tbl$rank), , drop = FALSE]
  n_total <- nrow(loading_tbl)
  rank_index <- as.integer(rank_index)[1]
  rank_index <- min(max(1L, rank_index), n_total)
  row_idx <- match(rank_index, loading_tbl$rank)
  cutoff_value <- as.numeric(loading_tbl$pc1_loading_abs[row_idx])
  cutoff_quantile <- 1 - (rank_index / n_total)
  fit_obj$loading_table <- loading_tbl
  fit_obj$cutoff <- cutoff_value
  fit_obj$top_n_used <- rank_index
  fit_obj$cutoff_quantile <- cutoff_quantile
  fit_obj$cutoff_method <- paste0(fit_obj$cutoff_method, "_selected")
  fit_obj$selected_reason <- selected_reason
  fit_obj$loading_table$split_class <- ifelse(fit_obj$loading_table$rank <= rank_index, "high_loading", "background_loading")
  fit_obj
}

# -----------------------------------------------------------------------------
# EVS SPLIT
# -----------------------------------------------------------------------------
build_eigenvector_split <- function(count_matrix, coldata, comparison_name) {
  dds_init <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = make_design_formula())
  dds_init <- dds_init[rowSums(counts(dds_init)) > 0, ]
  dds_init <- estimateSizeFactors(dds_init)

  retained_feature_ids <- rownames(dds_init)
  norm_counts_init <- as.data.frame(counts(dds_init, normalized = TRUE))
  raw_counts_init <- as.data.frame(count_matrix[retained_feature_ids, , drop = FALSE])

  sample_ids <- colnames(count_matrix)
  trt_ids <- sample_ids[coldata$condition == "trt"]
  untrt_ids <- sample_ids[coldata$condition == "untrt"]

  feature_metrics_trt <- compute_condition_feature_metrics(count_matrix[, trt_ids, drop = FALSE])
  feature_metrics_untrt <- compute_condition_feature_metrics(count_matrix[, untrt_ids, drop = FALSE])

  fit_trt <- compute_pc1_loading_table(norm_counts_init, trt_ids, preprocessing_label = "Normalized before EVS", feature_metrics = feature_metrics_trt)
  fit_untrt <- compute_pc1_loading_table(norm_counts_init, untrt_ids, preprocessing_label = "Normalized before EVS", feature_metrics = feature_metrics_untrt)
  fit_trt_raw <- compute_pc1_loading_table(raw_counts_init, trt_ids, preprocessing_label = "Raw counts before EVS", feature_metrics = feature_metrics_trt)
  fit_untrt_raw <- compute_pc1_loading_table(raw_counts_init, untrt_ids, preprocessing_label = "Raw counts before EVS", feature_metrics = feature_metrics_untrt)

  combined_cutoff_info <- resolve_combined_fourier_cutoff(
    fit_trt_loading_tbl = fit_trt$loading_table,
    fit_ctrl_loading_tbl = fit_untrt$loading_table,
    fixed_top_n = evs_fixed_top_n
  )
  combined_cutoff_info <- build_group_backup_crossings(combined_cutoff_info)

  final_shared_rank <- as.integer(combined_cutoff_info$top_n_actual)
  final_shared_reason <- combined_cutoff_info$selected_reason

  fit_trt <- apply_loading_cutoff(fit_trt, rank_index = final_shared_rank, selected_reason = paste0(final_shared_reason, "_applied_to_treatment"))
  fit_untrt <- apply_loading_cutoff(fit_untrt, rank_index = final_shared_rank, selected_reason = paste0(final_shared_reason, "_applied_to_control"))
  fit_trt_raw <- apply_loading_cutoff(fit_trt_raw, rank_index = final_shared_rank, selected_reason = paste0(final_shared_reason, "_projected_to_raw_treatment"))
  fit_untrt_raw <- apply_loading_cutoff(fit_untrt_raw, rank_index = final_shared_rank, selected_reason = paste0(final_shared_reason, "_projected_to_raw_control"))
  fit_trt_raw$cutoff_method <- "crossing_rank_projected_to_raw_selected"
  fit_untrt_raw$cutoff_method <- "crossing_rank_projected_to_raw_selected"

  shared_combined_tbl <- combined_cutoff_info$shared_combined_tbl
  if (is.null(shared_combined_tbl) || !nrow(shared_combined_tbl)) {
    stop("Shared combined loading table is missing after shared cutoff resolution.", call. = FALSE)
  }

  leading_edge_ids <- as.character(shared_combined_tbl$feature_id[shared_combined_tbl$combined_rank <= final_shared_rank])
  analyzed_feature_ids <- as.character(shared_combined_tbl$feature_id)
  remainder_ids <- setdiff(analyzed_feature_ids, leading_edge_ids)

  if (length(leading_edge_ids) == 0) stop("Leading-edge dataset is empty.", call. = FALSE)
  if (length(remainder_ids) == 0) stop("Remainder dataset is empty.", call. = FALSE)

  evs_cutoff_summary <- dplyr::bind_rows(
    data.frame(
      preprocessing = "normalized",
      group = "regime_shift_crossing",
      cutoff_mode = combined_cutoff_info$method,
      fixed_top_n_requested = evs_fixed_top_n,
      empiric_rank_selected = final_shared_rank,
      cutoff_quantile = fit_trt$cutoff_quantile,
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
      selected_reason = paste0(combined_cutoff_info$selected_reason, "_projected_to_raw"),
      stringsAsFactors = FALSE
    )
  )

  list(
    fit_trt = fit_trt,
    fit_untrt = fit_untrt,
    fit_trt_raw = fit_trt_raw,
    fit_untrt_raw = fit_untrt_raw,
    combined_cutoff_info = combined_cutoff_info,
    shared_combined_tbl = shared_combined_tbl,
    evs_cutoff_summary = evs_cutoff_summary,
    normalized_counts = norm_counts_init,
    raw_dataset = count_matrix,
    leading_edge_dataset = count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset = count_matrix[remainder_ids, , drop = FALSE]
  )
}

# -----------------------------------------------------------------------------
# PCA AND EVS FIGURES
# -----------------------------------------------------------------------------
condition_shapes <- c("untrt" = 21, "trt" = 25)
condition_fills <- c("untrt" = plot_palette$control, "trt" = plot_palette$treatment)
condition_labels <- c("untrt" = "Control", "trt" = "Treatment")

make_pca_summary_table <- function(pca_fit, dataset_name, preprocessing_label, group_label = NA_character_) {
  sdev <- pca_fit$sdev
  variance <- sdev^2
  prop_var <- variance / sum(variance)
  cum_var <- cumsum(prop_var)
  data.frame(
    dataset_name = dataset_name,
    preprocessing_label = preprocessing_label,
    group_label = group_label,
    principal_component = paste0("PC", seq_along(sdev)),
    standard_deviation = as.numeric(sdev),
    variance = as.numeric(variance),
    proportion_variance = as.numeric(prop_var),
    cumulative_proportion = as.numeric(cum_var),
    stringsAsFactors = FALSE
  )
}

plot_pca_variance_profile <- function(pca_fit, dataset_name, preprocessing_label, group_label = NULL) {
  pca_tbl <- make_pca_summary_table(pca_fit, dataset_name, preprocessing_label, ifelse(is.null(group_label), NA_character_, group_label))
  pca_tbl$principal_component <- factor(pca_tbl$principal_component, levels = pca_tbl$principal_component)
  group_label_chr <- if (length(group_label) && !is.null(group_label[1])) tolower(as.character(group_label[1])) else NA_character_
  bar_fill <- if (!is.na(group_label_chr) && group_label_chr %in% c("control", "untrt")) plot_palette$control else plot_palette$treatment

  ggplot(pca_tbl, aes(principal_component, proportion_variance)) +
    geom_col(fill = bar_fill, color = "white") +
    geom_line(aes(x = seq_along(principal_component), y = cumulative_proportion, group = 1), inherit.aes = FALSE, linewidth = LINE_WIDTH_BOUNDARY, colour = plot_palette$threshold) +
    geom_point(aes(x = seq_along(principal_component), y = cumulative_proportion), inherit.aes = FALSE, size = 1.6, colour = plot_palette$threshold) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(title = compact_title(paste(dataset_name, "|", pretty_group_label(group_label), "PCA"), width = 42), x = "Principal component", y = "Variance explained") +
    manuscript_theme() +
    theme(plot.margin = margin(t = 12, r = 14, b = 14, l = 14), legend.position = "none")
}

plot_pca_scatter <- function(pca_fit, dataset_label, group_label) {
  pca_var <- pca_fit$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)
  cond_key <- ifelse(group_label %in% c("control", "untrt"), "untrt", "trt")
  pca_df <- data.frame(Sample = rownames(pca_fit$x), PC1 = pca_fit$x[, 1], PC2 = pca_fit$x[, 2], Condition = cond_key, stringsAsFactors = FALSE)

  ggplot(pca_df, aes(PC1, PC2, label = Sample, shape = Condition, fill = Condition)) +
    geom_hline(yintercept = 0, linewidth = LINE_WIDTH_ZERO, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = LINE_WIDTH_ZERO, linetype = "dashed", colour = "grey70") +
    geom_point(size = 3.0, colour = "white", stroke = POINT_STROKE + 0.15) +
    geom_text_repel(size = 2.0, max.overlaps = 8, force = 1.0, box.padding = 0.22, point.padding = 0.10, min.segment.length = 0) +
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition", guide = guide_legend(override.aes = list(size = 3.0, fill = unname(condition_fills), colour = "white"))) +
    scale_fill_manual(values = condition_fills, labels = condition_labels, name = "Condition", guide = "none") +
    labs(title = pretty_group_label(group_label), x = paste0("PC1 (", pca_var_per[1], "%)"), y = paste0("PC2 (", pca_var_per[2], "%)")) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(plot.title = element_text(margin = margin(b = 4)), plot.margin = margin(t = 12, r = 18, b = 16, l = 16)) +
    guides(shape = guide_legend(order = 1))
}

plot_fourier_wave_map_single <- function(wave_obj, comparison_name, group_label) {
  if (is.null(wave_obj) || is.null(wave_obj$wave_map) || !nrow(wave_obj$wave_map)) return(NULL)
  df <- wave_obj$wave_map
  backup_info <- select_wave_backup_cutoff(wave_obj, group_label = tolower(pretty_group_label(group_label)))
  sc <- backup_info$selected_crossing
  cutoff_pct <- if (!is.null(sc) && nrow(sc)) sc$crossing_percentile[1] else NA_real_
  cutoff_rank <- if (!is.null(sc) && nrow(sc)) sc$crossing_rank[1] else NA_integer_
  line_col <- if (tolower(group_label) %in% c("treatment", "trt")) plot_palette$treatment else plot_palette$control

  p <- ggplot(df, aes(percentile)) +
    geom_line(aes(y = iod_amplitude, color = "IOD"), linewidth = 0.9) +
    geom_line(aes(y = cv2_amplitude, color = "CV²"), linewidth = 0.9) +
    scale_color_manual(values = c("IOD" = plot_palette$treatment, "CV²" = plot_palette$control)) +
    labs(
      title = paste0(pretty_group_label(group_label), " | local regime lines"),
      subtitle = compact_caption("This panel includes the backup within-group cutoff.", width = 88),
      x = "Percentile center",
      y = "Local amplitude",
      color = NULL
    ) +
    manuscript_theme()

  if (is.finite(cutoff_pct)) {
    ymax <- max(c(df$iod_amplitude, df$cv2_amplitude), na.rm = TRUE)
    p <- p +
      geom_vline(xintercept = cutoff_pct, linetype = "dotted", linewidth = 0.9, colour = line_col) +
      annotate("label", x = cutoff_pct, y = ymax, label = paste0(pretty_group_label(group_label), " backup\np = ", signif(cutoff_pct, 4), "\nrank = ", cutoff_rank), fill = "white", colour = line_col, size = 2.8, label.size = 0.15, vjust = -0.5)
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
    geom_line(aes(y = combined_iod_amplitude, color = "Composite IOD"), linewidth = 0.95) +
    geom_line(aes(y = combined_cv2_amplitude, color = "Composite CV²"), linewidth = 0.95) +
    scale_color_manual(values = c("Composite IOD" = plot_palette$treatment, "Composite CV²" = plot_palette$control)) +
    labs(
      title = paste0(comparison_name, " | two-line regime crossing"),
      subtitle = compact_caption("Dashed = shared cutoff. Dotted = backup treatment and control crossings.", width = 90),
      x = "Percentile center",
      y = "Composite local amplitude",
      color = NULL
    ) +
    manuscript_theme()

  if (is.finite(crossing_pct)) {
    p <- p + geom_vline(xintercept = crossing_pct, linetype = "dashed", linewidth = 0.95, colour = plot_palette$threshold) +
      annotate("label", x = crossing_pct, y = ymax, label = paste0("Shared cutoff\np = ", signif(crossing_pct, 4), "\nrank = ", crossing_rank), fill = "white", colour = plot_palette$threshold, size = 3.0, label.size = 0.15, vjust = -0.55)
  }
  if (is.finite(trt_pct)) {
    p <- p + geom_vline(xintercept = trt_pct, linetype = "dotted", linewidth = 0.85, colour = plot_palette$treatment) +
      annotate("label", x = trt_pct, y = ymax * 0.78, label = paste0("Treatment backup\np = ", signif(trt_pct, 4), "\nrank = ", trt_rank), fill = "white", colour = plot_palette$treatment, size = 2.7, label.size = 0.15)
  }
  if (is.finite(ctrl_pct)) {
    p <- p + geom_vline(xintercept = ctrl_pct, linetype = "dotted", linewidth = 0.85, colour = plot_palette$control) +
      annotate("label", x = ctrl_pct, y = ymax * 0.58, label = paste0("Control backup\np = ", signif(ctrl_pct, 4), "\nrank = ", ctrl_rank), fill = "white", colour = plot_palette$control, size = 2.7, label.size = 0.15)
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
    geom_hline(yintercept = 0, linewidth = 0.65, colour = "grey50") +
    geom_line(linewidth = 0.95, colour = plot_palette$threshold) +
    geom_vline(xintercept = crossing_pct, linetype = "dashed", linewidth = 0.95, colour = plot_palette$threshold) +
    geom_point(data = data.frame(percentile = crossing_pct, regime_difference = 0), aes(x = percentile, y = regime_difference), inherit.aes = FALSE, size = 1.8, colour = plot_palette$threshold) +
    annotate("label", x = crossing_pct, y = 0, label = paste0("rank = ", crossing_rank), fill = "white", colour = plot_palette$threshold, size = 2.8, label.size = 0.15, vjust = -0.8) +
    labs(title = paste0(comparison_name, " | regime-difference curve"), subtitle = compact_caption("Positive = IOD dominant, negative = CV2 dominant. First stable zero-crossing = EVS cutoff.", width = 90), x = "Percentile center", y = "IOD - CV2") +
    manuscript_theme()
}

plot_crossing_summary_panel <- function(combined_cutoff_info, comparison_name) {
  p1 <- plot_combined_fourier_wave_map(combined_cutoff_info, comparison_name)
  p2 <- plot_regime_difference_curve(combined_cutoff_info, comparison_name)
  if (is.null(p1) || is.null(p2)) return(NULL)
  arrangeGrob(p1, p2, ncol = 1, top = textGrob(paste0(comparison_name, " | regime-shift crossing summary"), gp = gpar(fontface = "bold", cex = 1.04)))
}

build_crossing_summary_table <- function(combined_cutoff_info, comparison_name) {
  sc <- combined_cutoff_info$selected_crossing
  if (is.null(sc) || !nrow(sc)) {
    return(data.frame(comparison_name = comparison_name, selected_reason = combined_cutoff_info$selected_reason, crossing_percentile = NA_real_, crossing_rank = NA_integer_, stringsAsFactors = FALSE))
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

build_backup_crossing_summary_table <- function(combined_cutoff_info, comparison_name) {
  add_row <- function(sc, label, reason) {
    if (is.null(sc) || !nrow(sc)) {
      data.frame(comparison_name = comparison_name, cutoff_type = label, selected_reason = reason, crossing_percentile = NA_real_, crossing_rank = NA_integer_, stringsAsFactors = FALSE)
    } else {
      data.frame(comparison_name = comparison_name, cutoff_type = label, selected_reason = reason, crossing_percentile = sc$crossing_percentile[1], crossing_rank = sc$crossing_rank[1], stringsAsFactors = FALSE)
    }
  }
  dplyr::bind_rows(
    add_row(combined_cutoff_info$selected_crossing, "shared_combined", combined_cutoff_info$selected_reason),
    add_row(combined_cutoff_info$trt_backup_crossing, "treatment_backup", combined_cutoff_info$trt_backup_reason),
    add_row(combined_cutoff_info$ctrl_backup_crossing, "control_backup", combined_cutoff_info$ctrl_backup_reason)
  )
}

# -----------------------------------------------------------------------------
# CORE DESEQ2 ANALYSIS
# -----------------------------------------------------------------------------
run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
  dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = make_design_formula())
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, betaPrior = FALSE, quiet = TRUE)

  res <- results(dds, contrast = c("condition", "trt", "untrt"), alpha = alpha_level, independentFiltering = FALSE)
  res_strong <- results(dds, contrast = c("condition", "trt", "untrt"), lfcThreshold = lfc_boundary, altHypothesis = "greaterAbs", independentFiltering = FALSE)
  res_weak <- results(dds, contrast = c("condition", "trt", "untrt"), lfcThreshold = lfc_boundary, altHypothesis = "lessAbs", independentFiltering = FALSE)

  res_df <- as.data.frame(res, stringsAsFactors = FALSE)
  res_df$feature_id <- as.character(rownames(res_df))

  valid_stat <- is.finite(res_df$stat) & !is.na(res_df$stat)
  stat_vec <- as.numeric(res_df$stat[valid_stat])
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]

  message(sprintf("[%s] Mean Wald stat sent to fdrtool: %.4f", dataset_name, mean(stat_vec, na.rm = TRUE)))
  fdr_fit <- run_empirical_null_fdrtool(stat_vec, dataset_name = dataset_name)

  res_df$empirical_p <- NA_real_
  res_df$empirical_q <- NA_real_
  res_df$lfdr <- NA_real_
  res_df$empirical_p[valid_stat] <- as.numeric(fdr_fit$pval)
  res_df$empirical_q[valid_stat] <- as.numeric(fdr_fit$qval)
  res_df$lfdr[valid_stat] <- as.numeric(fdr_fit$lfdr)

  coef_name <- get_condition_coef(dds)
  shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)
  shr_df <- as.data.frame(shr, stringsAsFactors = FALSE)
  shr_df$feature_id <- as.character(rownames(shr_df))
  res_df <- dplyr::left_join(res_df, shr_df[, c("feature_id", "log2FoldChange", "lfcSE"), drop = FALSE], by = "feature_id", suffix = c("", "_shrunk"))
  colnames(res_df)[colnames(res_df) == "log2FoldChange_shrunk"] <- "lfc_shrunk"
  colnames(res_df)[colnames(res_df) == "lfcSE_shrunk"] <- "lfcSE_shrunk"

  res_strong_df <- as.data.frame(res_strong, stringsAsFactors = FALSE)
  res_strong_df$feature_id <- as.character(rownames(res_strong_df))
  res_weak_df <- as.data.frame(res_weak, stringsAsFactors = FALSE)
  res_weak_df$feature_id <- as.character(rownames(res_weak_df))

  res_df <- dplyr::left_join(res_df, res_strong_df[, c("feature_id", "padj"), drop = FALSE], by = "feature_id", suffix = c("", "_strong"))
  res_df <- dplyr::left_join(res_df, res_weak_df[, c("feature_id", "padj"), drop = FALSE], by = "feature_id", suffix = c("", "_weak"))
  colnames(res_df)[colnames(res_df) == "padj_strong"] <- "padj_strong_effect"
  colnames(res_df)[colnames(res_df) == "padj_weak"] <- "padj_weak_effect"
  res_df$resGA_padj <- res_df$padj_strong_effect
  res_df$resLA_padj <- res_df$padj_weak_effect

  res_df$wald_pvalue <- res_df$pvalue
  res_df$neglog10_wald_pvalue <- safe_neglog10(res_df$wald_pvalue)
  res_df$neglog10_empirical_p <- safe_neglog10(res_df$empirical_p)

  hc_empirical_info <- resolve_hc_threshold(res_df$empirical_p, dataset_name, source_label = "empirical")
  hc_wald_info <- resolve_hc_threshold(res_df$wald_pvalue, dataset_name, source_label = "wald")

  hc_empirical_threshold <- hc_empirical_info$threshold
  hc_wald_threshold <- hc_wald_info$threshold

  empirical_p_floored <- ifelse(is.na(res_df$empirical_p), NA_real_, pmax(res_df$empirical_p, 1e-300))
  res_df$HBFSS <- abs(res_df$lfc_shrunk * log10(empirical_p_floored))

  if (is.na(hc_empirical_threshold)) {
    hbfss_threshold_dataset <- NA_real_
    res_df$hc_empirical_pass <- FALSE
    res_df$HBFSS_significant <- FALSE
  } else {
    hbfss_threshold_dataset <- abs(log10(hc_empirical_threshold)) * lfc_boundary
    res_df$hc_empirical_pass <- !is.na(res_df$empirical_p) & (res_df$empirical_p <= hc_empirical_threshold)
    res_df$HBFSS_significant <- !is.na(res_df$HBFSS) & (res_df$HBFSS >= hbfss_threshold_dataset) & res_df$hc_empirical_pass
  }

  if (is.na(hc_wald_threshold)) {
    res_df$hc_wald_pass <- FALSE
  } else {
    res_df$hc_wald_pass <- !is.na(res_df$wald_pvalue) & (res_df$wald_pvalue <= hc_wald_threshold)
  }

  res_df$regulation_direction <- ifelse(is.na(res_df$lfc_shrunk), NA_character_, ifelse(res_df$lfc_shrunk > 0, "upregulated", ifelse(res_df$lfc_shrunk < 0, "downregulated", "no_change")))
  res_df$raw_lfc_pass <- !is.na(res_df$log2FoldChange) & (abs(res_df$log2FoldChange) >= lfc_boundary)
  res_df$shrunk_lfc_pass <- !is.na(res_df$lfc_shrunk) & (abs(res_df$lfc_shrunk) >= lfc_boundary)

  res_df$deseq2_strong_call <- !is.na(res_df$resGA_padj) & (res_df$resGA_padj < alpha_level) & res_df$shrunk_lfc_pass
  res_df$deseq2_weak_call_raw <- !is.na(res_df$resLA_padj) & (res_df$resLA_padj < alpha_level) & !res_df$shrunk_lfc_pass
  res_df$deseq2_weak_call <- res_df$deseq2_weak_call_raw & res_df$hc_wald_pass
  res_df$standard_significant <- res_df$deseq2_strong_call | res_df$deseq2_weak_call
  res_df$HBFSS_only_call <- res_df$HBFSS_significant & !res_df$standard_significant
  res_df$overlap_call <- res_df$HBFSS_significant & res_df$deseq2_strong_call

  res_df$effect_class <- dplyr::case_when(
    res_df$deseq2_strong_call ~ "strong_effect",
    res_df$deseq2_weak_call ~ "weak_effect",
    TRUE ~ "intermediate"
  )

  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts$feature_id <- as.character(rownames(norm_counts))
  mm <- as.data.frame(mcols(dds), stringsAsFactors = FALSE)
  mm$feature_id <- as.character(rownames(mm))
  disp_cols_available <- intersect(c("feature_id", "baseMean", "dispGeneEst", "dispFit", "dispersion", "dispIter", "dispOutlier"), colnames(mm))
  disp_df <- mm[, disp_cols_available, drop = FALSE]

  annot_df$feature_id <- as.character(annot_df$feature_id)
  annot_df$gene_symbol <- as.character(annot_df$gene_symbol)
  annot_df <- annot_df %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(is.na(gene_symbol), "", trimws(gene_symbol))) %>%
    dplyr::arrange(feature_id, dplyr::desc(gene_symbol != ""), gene_symbol) %>%
    dplyr::distinct(feature_id, .keep_all = TRUE) %>%
    dplyr::mutate(gene_symbol = dplyr::na_if(gene_symbol, ""))

  final_df <- res_df %>%
    dplyr::left_join(annot_df, by = "feature_id") %>%
    dplyr::left_join(norm_counts, by = "feature_id") %>%
    dplyr::left_join(disp_df, by = "feature_id")

  final_df$dataset_name <- dataset_name
  final_df$hc_empirical_threshold_dataset <- hc_empirical_threshold
  final_df$hc_wald_threshold_dataset <- hc_wald_threshold
  final_df$hc_empirical_threshold_source <- hc_empirical_info$source
  final_df$hc_wald_threshold_source <- hc_wald_info$source
  final_df$hbfss_threshold_dataset <- hbfss_threshold_dataset

  preferred_cols <- c(
    "dataset_name", "feature_id", "gene_symbol", "baseMean",
    "lfc_shrunk", "regulation_direction", "lfdr",
    "pvalue", "wald_pvalue", "padj", "empirical_p", "empirical_q",
    "HBFSS", "hc_wald_threshold_dataset", "hc_empirical_threshold_dataset", "hbfss_threshold_dataset",
    "resLA_padj", "resGA_padj", "hc_wald_pass", "hc_empirical_pass",
    "deseq2_strong_call", "deseq2_weak_call_raw", "deseq2_weak_call",
    "standard_significant", "HBFSS_significant", "effect_class"
  )
  final_df <- final_df[, c(intersect(preferred_cols, names(final_df)), setdiff(names(final_df), preferred_cols)), drop = FALSE]

  list(
    dds = dds,
    results = final_df,
    base_mean_vec = final_df$baseMean[!is.na(final_df$baseMean)],
    hc_empirical_threshold = hc_empirical_threshold,
    hc_wald_threshold = hc_wald_threshold,
    hbfss_threshold = hbfss_threshold_dataset
  )
}

# -----------------------------------------------------------------------------
# VOLCANO HELPERS
# -----------------------------------------------------------------------------
volcano_class_colors <- c(
  background = plot_palette$background,
  weak = plot_palette$weak,
  strong = plot_palette$strong,
  hbfss = plot_palette$hbfss,
  overlap = plot_palette$overlap
)

volcano_class_labels <- c(
  background = "Background",
  weak = "DC2 weak + HC",
  strong = "DC2 strong",
  hbfss = "HBFSS only",
  overlap = "DC2 strong + HBFSS"
)

build_plot_specific_volcano_classes <- function(df, plot_type = c("standard", "hbfss")) {
  plot_type <- match.arg(plot_type)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  required_cols <- c("gene_symbol", "lfc_shrunk", "wald_pvalue", "neglog10_wald_pvalue", "empirical_p", "neglog10_empirical_p", "deseq2_strong_call", "deseq2_weak_call", "HBFSS_significant", "HBFSS_only_call", "overlap_call")
  assert_required_columns(df, required_cols, object_name = "volcano input df")

  df$has_valid_gene_symbol <- !is.na(df$gene_symbol) & grepl("^[A-Za-z0-9._-]+$", trimws(df$gene_symbol))
  df$gene_symbol_plot <- ifelse(df$has_valid_gene_symbol, trimws(df$gene_symbol), NA_character_)

  if (plot_type == "standard") {
    df$plot_y <- df$neglog10_wald_pvalue
  } else {
    df$plot_y <- df$neglog10_empirical_p
  }

  df$volcano_class <- dplyr::case_when(
    df$overlap_call ~ "overlap",
    df$deseq2_strong_call ~ "strong",
    df$deseq2_weak_call ~ "weak",
    df$HBFSS_only_call ~ "hbfss",
    TRUE ~ "background"
  )
  df$volcano_class <- factor(df$volcano_class, levels = c("background", "weak", "strong", "hbfss", "overlap"))
  df
}

select_volcano_labels <- function(df, y_col = "plot_y", n_labels = 20) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!nrow(df)) return(df[0, , drop = FALSE])
  df <- df[df$has_valid_gene_symbol, , drop = FALSE]
  if (!nrow(df)) return(df[0, , drop = FALSE])
  df <- df[df$volcano_class %in% c("overlap", "strong", "weak", "hbfss"), , drop = FALSE]
  if (!nrow(df)) return(df[0, , drop = FALSE])

  pri <- dplyr::case_when(
    df$volcano_class == "overlap" ~ 1L,
    df$volcano_class == "strong" ~ 2L,
    df$volcano_class == "weak" ~ 3L,
    df$volcano_class == "hbfss" ~ 4L,
    TRUE ~ 9L
  )

  yy <- suppressWarnings(as.numeric(df[[y_col]]))
  yy[!is.finite(yy)] <- -Inf
  hh <- suppressWarnings(as.numeric(df$HBFSS))
  hh[!is.finite(hh)] <- -Inf

  ord <- order(pri, -yy, -hh, -abs(df$lfc_shrunk), na.last = TRUE)
  df <- df[ord, , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol_plot), , drop = FALSE]
  df[seq_len(min(n_labels, nrow(df))), , drop = FALSE]
}

volcano_label_layer <- function(lab_df) {
  if (!nrow(lab_df)) return(NULL)
  ggrepel::geom_text_repel(
    data = lab_df,
    aes(label = gene_symbol_plot),
    size = 1.95,
    seed = 1,
    max.overlaps = 35,
    force = 1.25,
    force_pull = 0.5,
    box.padding = 0.38,
    point.padding = 0.22,
    min.segment.length = 0,
    segment.alpha = 0.6,
    segment.size = 0.22
  )
}

add_hbfss_boundary_layer <- function(p, df) {
  threshold <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))
  hc_p_threshold <- suppressWarnings(as.numeric(df$hc_empirical_threshold_dataset[1]))
  if (!is.finite(threshold) || is.na(threshold) || threshold <= 0) return(p)

  finite_lfc <- suppressWarnings(as.numeric(df$lfc_shrunk))
  finite_lfc <- finite_lfc[is.finite(finite_lfc) & !is.na(finite_lfc)]
  max_abs_lfc <- max(abs(finite_lfc), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || is.na(max_abs_lfc) || max_abs_lfc <= 0) max_abs_lfc <- lfc_boundary * 3

  x_abs <- seq(from = max(0.05, min(lfc_boundary, max_abs_lfc)), to = max_abs_lfc, length.out = 400)
  y_hyper <- threshold / x_abs

  boundary_df <- data.frame(
    lfc_shrunk = c(-rev(x_abs), x_abs),
    neglog10_empirical_p = c(rev(y_hyper), y_hyper),
    stringsAsFactors = FALSE
  )
  boundary_df <- boundary_df[is.finite(boundary_df$neglog10_empirical_p), , drop = FALSE]

  if (nrow(boundary_df)) {
    obs_y_max <- suppressWarnings(max(df$neglog10_empirical_p[is.finite(df$neglog10_empirical_p)], na.rm = TRUE))
    if (!is.finite(obs_y_max)) obs_y_max <- max(y_hyper[is.finite(y_hyper)], na.rm = TRUE)
    label_y <- min(max(y_hyper, na.rm = TRUE), obs_y_max * 0.92)
    if (!is.finite(label_y)) label_y <- max(y_hyper, na.rm = TRUE)

    p <- p +
      geom_path(data = boundary_df, aes(x = lfc_shrunk, y = neglog10_empirical_p), inherit.aes = FALSE, linetype = "dashed", linewidth = LINE_WIDTH_BOUNDARY, colour = plot_palette$hbfss) +
      annotate("label", x = max_abs_lfc * 0.82, y = label_y, label = paste0("HBFSS = ", signif(threshold, 4), " / |LFC|"), fill = "white", colour = plot_palette$hbfss, size = 2.8, label.size = 0.15)
  }

  if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold) && hc_p_threshold > 0 && hc_p_threshold < 1) {
    p <- p +
      geom_hline(yintercept = -log10(hc_p_threshold), linetype = "dotted", linewidth = LINE_WIDTH_BOUNDARY, colour = plot_palette$threshold) +
      annotate("label", x = 0, y = -log10(hc_p_threshold), label = paste0("HC empirical p = ", signif(hc_p_threshold, 4)), fill = "white", colour = plot_palette$threshold, size = 2.8, label.size = 0.15, vjust = -0.7)
  }
  p
}

plot_standard_volcano <- function(df, dataset_name) {
  df <- build_plot_specific_volcano_classes(df, plot_type = "standard")
  lab_df <- select_volcano_labels(df, y_col = "plot_y", n_labels = 20)

  p <- ggplot(df, aes(lfc_shrunk, plot_y)) +
    geom_point(aes(color = volcano_class), alpha = POINT_ALPHA_PRIMARY, size = POINT_SIZE_PRIMARY + 0.45, shape = 16) +
    geom_vline(xintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", linewidth = LINE_WIDTH_BOUNDARY, colour = plot_palette$threshold) +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = LINE_WIDTH_ZERO, colour = "grey45") +
    geom_hline(yintercept = -log10(alpha_level), linetype = "dashed", linewidth = LINE_WIDTH_BOUNDARY, colour = plot_palette$threshold) +
    scale_color_manual(values = volcano_class_colors, labels = volcano_class_labels, drop = FALSE, name = "Call class") +
    labs(title = pretty_dataset_label(dataset_name), x = "Shrunken log2 fold change", y = expression(-log[10](Wald~p)), caption = compact_caption("Colored markers match the final call classes. Weak DC2 calls require the Wald HC threshold.")) +
    coord_cartesian(clip = "off") +
    plot_expand_xy() +
    manuscript_theme() +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 3.8, alpha = 1)))

  hc_wald <- suppressWarnings(as.numeric(df$hc_wald_threshold_dataset[1]))
  if (is.finite(hc_wald) && !is.na(hc_wald) && hc_wald > 0 && hc_wald < 1) {
    p <- p +
      geom_hline(yintercept = -log10(hc_wald), linetype = "dotted", linewidth = LINE_WIDTH_BOUNDARY, colour = plot_palette$weak) +
      annotate("label", x = 0, y = -log10(hc_wald), label = paste0("HC Wald p = ", signif(hc_wald, 4)), fill = "white", colour = plot_palette$weak, size = 2.8, label.size = 0.15, vjust = -0.7)
  }

  if (nrow(lab_df) > 0) p <- p + volcano_label_layer(lab_df)
  p
}

plot_hbfss_volcano <- function(df, dataset_name) {
  df <- build_plot_specific_volcano_classes(df, plot_type = "hbfss")
  lab_df <- select_volcano_labels(df, y_col = "plot_y", n_labels = 20)

  p <- ggplot(df, aes(lfc_shrunk, plot_y)) +
    geom_point(aes(color = volcano_class), alpha = POINT_ALPHA_PRIMARY, size = POINT_SIZE_PRIMARY + 0.45, shape = 16) +
    geom_vline(xintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", linewidth = LINE_WIDTH_BOUNDARY, colour = plot_palette$threshold) +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = LINE_WIDTH_ZERO, colour = "grey45") +
    scale_color_manual(values = volcano_class_colors, labels = volcano_class_labels, drop = FALSE, name = "Call class") +
    labs(title = pretty_dataset_label(dataset_name), x = "Shrunken log2 fold change", y = expression(-log[10](p[empirical])), caption = compact_caption("HBFSS uses empirical p values. The HBFSS boundary and empirical HC threshold are both shown.")) +
    coord_cartesian(clip = "off") +
    plot_expand_xy() +
    manuscript_theme() +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 3.8, alpha = 1)))

  p <- add_hbfss_boundary_layer(p, df)
  if (nrow(lab_df) > 0) p <- p + volcano_label_layer(lab_df)
  p
}

plot_publication_volcano_panel <- function(df, dataset_name) {
  plot_hbfss_volcano(df, dataset_name)
}

plot_dispersion_panel_for_dataset <- function(df, dataset_name) {
  required_cols <- c("baseMean", "dispersion", "volcano_color")
  plot_df <- data.frame(
    baseMean = suppressWarnings(as.numeric(df$baseMean)),
    dispersion = suppressWarnings(as.numeric(df$dispersion)),
    class = build_plot_specific_volcano_classes(df, plot_type = "hbfss")$volcano_class,
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df[is.finite(plot_df$baseMean) & plot_df$baseMean > 0 & is.finite(plot_df$dispersion) & plot_df$dispersion > 0, , drop = FALSE]
  if (!nrow(plot_df)) return(NULL)
  ggplot(plot_df, aes(baseMean, dispersion, color = class)) +
    geom_point(alpha = POINT_ALPHA_DISP, size = POINT_SIZE_DISP) +
    scale_x_log10(labels = label_number(accuracy = 0.1)) +
    scale_y_log10(labels = label_number(accuracy = 0.1)) +
    scale_color_manual(values = volcano_class_colors, labels = volcano_class_labels, drop = FALSE, name = "Call class") +
    labs(title = compact_title(pretty_dataset_label(dataset_name), width = 42), x = "baseMean (log10 scale)", y = "Final dispersion (log10 scale)") +
    manuscript_theme() +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 3.8, alpha = 1)))
}

compute_dataset_pca_plot <- function(count_df, coldata, dataset_name, preprocessing = c("normalized", "raw_counts"), precomputed_matrix = NULL) {
  preprocessing <- match.arg(preprocessing)
  if (!is.null(precomputed_matrix)) {
    x <- as.matrix(precomputed_matrix)
  } else if (preprocessing == "normalized") {
    dds <- DESeqDataSetFromMatrix(countData = count_df, colData = coldata, design = make_design_formula())
    dds <- dds[rowSums(counts(dds)) > 0, ]
    dds <- estimateSizeFactors(dds)
    x <- counts(dds, normalized = TRUE)
  } else {
    x <- as.matrix(count_df)
  }

  pca_fit <- prcomp(t(x), scale. = FALSE, rank. = 2)
  pca_var <- pca_fit$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)

  pca_df <- data.frame(Sample = rownames(pca_fit$x), PC1 = pca_fit$x[, 1], PC2 = pca_fit$x[, 2], Condition = as.character(coldata[rownames(pca_fit$x), "condition"]), stringsAsFactors = FALSE)
  ggplot(pca_df, aes(PC1, PC2, label = Sample, shape = Condition, fill = Condition)) +
    geom_hline(yintercept = 0, linewidth = LINE_WIDTH_ZERO, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = LINE_WIDTH_ZERO, linetype = "dashed", colour = "grey70") +
    geom_point(size = 3, colour = "white", stroke = POINT_STROKE + 0.15) +
    geom_text_repel(size = 2.0, max.overlaps = 8, force = 1.0, box.padding = 0.22, point.padding = 0.10, min.segment.length = 0) +
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition", guide = guide_legend(override.aes = list(size = 3.0, fill = unname(condition_fills), colour = "white"))) +
    scale_fill_manual(values = condition_fills, labels = condition_labels, name = "Condition", guide = "none") +
    labs(title = compact_title(pretty_dataset_label(dataset_name)), x = paste0("PC1 (", pca_var_per[1], "%)"), y = paste0("PC2 (", pca_var_per[2], "%)")) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(legend.position = "none", plot.margin = margin(t = 12, r = 18, b = 16, l = 16))
}

subset_normalized_counts_by_ids <- function(normalized_counts, feature_ids, dataset_label, comparison_name, strict = TRUE) {
  feature_ids <- as.character(feature_ids)
  idx <- match(feature_ids, rownames(normalized_counts))
  missing_ids <- unique(feature_ids[is.na(idx)])
  if (length(missing_ids) > 0L) {
    msg <- sprintf("[%s][%s] %d feature_id(s) were not found in normalized_counts.", comparison_name, dataset_label, length(missing_ids))
    if (isTRUE(strict)) stop(msg, call. = FALSE) else warning(msg)
  }
  if (!isTRUE(strict)) idx <- idx[!is.na(idx)]
  as.matrix(normalized_counts[idx, , drop = FALSE])
}

# -----------------------------------------------------------------------------
# EXPORT BUILDERS
# -----------------------------------------------------------------------------
build_pc1_feature_export <- function(analysis_results, annot_df, comparison_name, evs, tab_dir) {
  annot_df$feature_id <- as.character(annot_df$feature_id)
  annot_df$gene_symbol <- as.character(annot_df$gene_symbol)
  base_export <- annot_df %>% dplyr::distinct(feature_id, .keep_all = TRUE)

  add_dataset_columns <- function(df, key, prefix) {
    if (is.null(analysis_results[[key]]) || is.null(analysis_results[[key]]$results)) return(df)
    src <- analysis_results[[key]]$results
    keep_cols <- intersect(c("feature_id", "baseMean", "dispGeneEst", "dispFit", "dispersion"), names(src))
    src <- src[, keep_cols, drop = FALSE]
    names(src)[names(src) == "baseMean"] <- paste0("baseMean_", prefix)
    names(src)[names(src) == "dispGeneEst"] <- paste0("dispGeneEst_", prefix)
    names(src)[names(src) == "dispFit"] <- paste0("dispFit_", prefix)
    names(src)[names(src) == "dispersion"] <- paste0("dispersion_", prefix)
    dplyr::left_join(df, src, by = "feature_id")
  }

  load_tbl <- function(tbl, col_name) {
    x <- tbl[, c("feature_id", "pc1_loading"), drop = FALSE]
    names(x)[names(x) == "pc1_loading"] <- col_name
    x
  }

  out <- base_export
  out <- add_dataset_columns(out, "raw_dataset", "original")
  out <- add_dataset_columns(out, "leading_edge_dataset", "leading_edge")
  out <- add_dataset_columns(out, "remainder_dataset", "remainder")
  out <- dplyr::left_join(out, load_tbl(evs$fit_trt$loading_table, "pc1_loading_trt_normalized"), by = "feature_id")
  out <- dplyr::left_join(out, load_tbl(evs$fit_untrt$loading_table, "pc1_loading_ctrl_normalized"), by = "feature_id")
  out <- dplyr::left_join(out, load_tbl(evs$fit_trt_raw$loading_table, "pc1_loading_trt_raw"), by = "feature_id")
  out <- dplyr::left_join(out, load_tbl(evs$fit_untrt_raw$loading_table, "pc1_loading_ctrl_raw"), by = "feature_id")
  out$comparison_name <- comparison_name
  save_csv(out, file.path(tab_dir, paste0(comparison_name, "_PC1_loadings_baseMean_dispersion_export.csv")))
  out
}

save_cross_dataset_comparison_panels <- function(comparison_name, analysis_results, cmp_dir, dataset_list, coldata, normalized_dataset_list = NULL) {
  keys_present <- base::intersect(as.character(dataset_key_order), as.character(names(analysis_results)))
  if (!length(keys_present)) return(invisible(NULL))

  make_panel <- function(grobs, title_text, file_name, width_unit = 6.8, height = 6.1) {
    grobs <- Filter(Negate(is.null), grobs)
    if (!length(grobs)) return(invisible(NULL))
    panel <- do.call(arrangeGrob, c(grobs, list(ncol = length(grobs), top = textGrob(title_text, gp = gpar(fontface = "bold", cex = 1.10)))))
    save_grob(panel, file.path(cmp_dir, file_name), width = width_unit * length(grobs), height = height)
  }

  std_grobs <- lapply(keys_present, function(k) plot_standard_volcano(analysis_results[[k]]$results, analysis_results[[k]]$summary$dataset_name[1]))
  make_panel(std_grobs, paste(comparison_name, "| DESeq2 volcanoes"), paste0(comparison_name, "_cross_dataset_standard_volcano_panel.png"))

  hbfss_grobs <- lapply(keys_present, function(k) plot_hbfss_volcano(analysis_results[[k]]$results, analysis_results[[k]]$summary$dataset_name[1]))
  make_panel(hbfss_grobs, paste(comparison_name, "| HBFSS volcanoes"), paste0(comparison_name, "_cross_dataset_HBFSS_volcano_panel.png"))

  pub_grobs <- lapply(keys_present, function(k) plot_publication_volcano_panel(analysis_results[[k]]$results, analysis_results[[k]]$summary$dataset_name[1]))
  make_panel(pub_grobs, paste(comparison_name, "| Publication volcanoes"), paste0(comparison_name, "_cross_dataset_publication_volcano_panel.png"))

  disp_grobs <- lapply(keys_present, function(k) plot_dispersion_panel_for_dataset(analysis_results[[k]]$results, analysis_results[[k]]$summary$dataset_name[1]))
  make_panel(disp_grobs, paste(comparison_name, "| Dispersion"), paste0(comparison_name, "_cross_dataset_dispersion_panel.png"), width_unit = 6.4, height = 5.9)

  pca_norm_grobs <- lapply(keys_present, function(k) compute_dataset_pca_plot(dataset_list[[k]], coldata, analysis_results[[k]]$summary$dataset_name[1], preprocessing = "normalized", precomputed_matrix = if (!is.null(normalized_dataset_list)) normalized_dataset_list[[k]] else NULL))
  pca_raw_grobs <- lapply(keys_present, function(k) compute_dataset_pca_plot(dataset_list[[k]], coldata, analysis_results[[k]]$summary$dataset_name[1], preprocessing = "raw_counts"))
  pca_norm_grobs <- Filter(Negate(is.null), pca_norm_grobs)
  pca_raw_grobs <- Filter(Negate(is.null), pca_raw_grobs)
  if (length(pca_norm_grobs) > 0 || length(pca_raw_grobs) > 0) {
    pca_panel <- arrangeGrob(grobs = c(pca_norm_grobs, pca_raw_grobs), ncol = max(1L, max(length(pca_norm_grobs), length(pca_raw_grobs))), top = textGrob(paste(comparison_name, "| Top: normalized before EVS. Bottom: raw before EVS."), gp = gpar(fontface = "bold", cex = 1.08)))
    save_grob(pca_panel, file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_PCA_normalized_vs_unnormalized.png")), width = 6.2 * max(1L, max(length(pca_norm_grobs), length(pca_raw_grobs))), height = 10.5)
  }

  summary_table <- dplyr::bind_rows(lapply(keys_present, function(k) {
    sm <- analysis_results[[k]]$summary
    data.frame(
      comparison_name = comparison_name,
      dataset_key = k,
      dataset_label = unname(dataset_key_labels[k]),
      dataset_name = sm$dataset_name[1],
      n_features = sm$n_features[1],
      n_standard_significant = sm$n_standard_significant[1],
      n_HBFSS_significant = sm$n_HBFSS_significant[1],
      n_overlap = sm$n_overlap_significant[1],
      n_strong_effect = sm$n_strong_effect[1],
      n_weak_effect = sm$n_weak_effect[1],
      hc_empirical_threshold = sm$hc_empirical_threshold[1],
      hc_wald_threshold = sm$hc_wald_threshold[1],
      hbfss_threshold = sm$hbfss_threshold[1],
      evs_fixed_top_n = evs_fixed_top_n,
      stringsAsFactors = FALSE
    )
  }))
  save_csv(summary_table, file.path(cmp_dir, paste0(comparison_name, "_cross_dataset_summary_for_figure_panels.csv")))
  invisible(NULL)
}

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  cmp_dir <- file.path(output_dir, comparison_name)
  tab_dir <- file.path(cmp_dir, "tables")
  fig_raw_dir <- file.path(cmp_dir, "figures_raw")
  fig_le_dir <- file.path(cmp_dir, "figures_leading_edge")
  fig_rem_dir <- file.path(cmp_dir, "figures_remainder")
  for (d in c(cmp_dir, tab_dir, fig_raw_dir, fig_le_dir, fig_rem_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

  evs <- build_eigenvector_split(count_matrix, coldata, comparison_name)
  save_csv(transform(evs$evs_cutoff_summary, comparison_name = comparison_name), file.path(tab_dir, paste0(comparison_name, "_EVS_cutoff_summary.csv")))
  if (!is.null(evs$combined_cutoff_info$trt_wave_obj$wave_map) && nrow(evs$combined_cutoff_info$trt_wave_obj$wave_map)) save_csv(evs$combined_cutoff_info$trt_wave_obj$wave_map, file.path(tab_dir, paste0(comparison_name, "_treatment_local_fourier_map.csv")))
  if (!is.null(evs$combined_cutoff_info$ctrl_wave_obj$wave_map) && nrow(evs$combined_cutoff_info$ctrl_wave_obj$wave_map)) save_csv(evs$combined_cutoff_info$ctrl_wave_obj$wave_map, file.path(tab_dir, paste0(comparison_name, "_control_local_fourier_map.csv")))
  if (!is.null(evs$combined_cutoff_info$combined_wave_map) && nrow(evs$combined_cutoff_info$combined_wave_map)) save_csv(evs$combined_cutoff_info$combined_wave_map, file.path(tab_dir, paste0(comparison_name, "_first_stable_crossing_map.csv")))
  if (!is.null(evs$combined_cutoff_info$crossing_table) && nrow(evs$combined_cutoff_info$crossing_table)) save_csv(evs$combined_cutoff_info$crossing_table, file.path(tab_dir, paste0(comparison_name, "_all_regime_crossings.csv")))
  save_csv(build_crossing_summary_table(evs$combined_cutoff_info, comparison_name), file.path(tab_dir, paste0(comparison_name, "_regime_shift_crossing_summary.csv")))
  save_csv(build_backup_crossing_summary_table(evs$combined_cutoff_info, comparison_name), file.path(tab_dir, paste0(comparison_name, "_backup_crossing_summary.csv")))

  save_grob(safe_plot_build(plot_fourier_wave_map_single(evs$combined_cutoff_info$trt_wave_obj, comparison_name, "treatment"), paste0(comparison_name, ": treatment local fourier wave map")), file.path(cmp_dir, paste0(comparison_name, "_treatment_local_fourier_wave_map.png")), width = 12, height = 10)
  save_grob(safe_plot_build(plot_fourier_wave_map_single(evs$combined_cutoff_info$ctrl_wave_obj, comparison_name, "control"), paste0(comparison_name, ": control local fourier wave map")), file.path(cmp_dir, paste0(comparison_name, "_control_local_fourier_wave_map.png")), width = 12, height = 10)
  save_grob(safe_plot_build(plot_combined_fourier_wave_map(evs$combined_cutoff_info, comparison_name), paste0(comparison_name, ": combined local fourier wave map")), file.path(cmp_dir, paste0(comparison_name, "_combined_local_fourier_wave_map.png")), width = 12, height = 6)
  save_grob(safe_plot_build(plot_regime_difference_curve(evs$combined_cutoff_info, comparison_name), paste0(comparison_name, ": regime difference curve")), file.path(cmp_dir, paste0(comparison_name, "_combined_regime_difference_curve.png")), width = 12, height = 6)
  save_grob(safe_plot_build(plot_crossing_summary_panel(evs$combined_cutoff_info, comparison_name), paste0(comparison_name, ": crossing summary panel")), file.path(cmp_dir, paste0(comparison_name, "_crossing_summary_panel.png")), width = 10, height = 10)

  evs_pca_scatter <- safe_plot_build(arrangeGrob(plot_pca_scatter(evs$fit_trt$pca_fit, comparison_name, "treatment"), plot_pca_scatter(evs$fit_untrt$pca_fit, comparison_name, "control"), plot_pca_scatter(evs$fit_trt_raw$pca_fit, comparison_name, "treatment"), plot_pca_scatter(evs$fit_untrt_raw$pca_fit, comparison_name, "control"), ncol = 2, top = textGrob(paste0(comparison_name, " | EVS PCA"), gp = gpar(fontface = "bold", cex = 1.02)), bottom = textGrob("Top: normalized before EVS. Bottom: raw before EVS. Left: treatment. Right: control.", gp = gpar(cex = 0.86))), paste0(comparison_name, ": EVS PCA scatter panel"))
  save_grob(evs_pca_scatter, file.path(cmp_dir, "EVS_PCA_scatter_combined.png"), width = 18, height = 13)

  if (isTRUE(export_optional_evs_variance_profiles)) {
    variance_panel <- safe_plot_build(arrangeGrob(plot_pca_variance_profile(evs$fit_trt$pca_fit, comparison_name, evs$fit_trt$preprocessing_label, "treatment"), plot_pca_variance_profile(evs$fit_untrt$pca_fit, comparison_name, evs$fit_untrt$preprocessing_label, "control"), plot_pca_variance_profile(evs$fit_trt_raw$pca_fit, comparison_name, evs$fit_trt_raw$preprocessing_label, "treatment"), plot_pca_variance_profile(evs$fit_untrt_raw$pca_fit, comparison_name, evs$fit_untrt_raw$preprocessing_label, "control"), ncol = 2, top = textGrob(paste0(comparison_name, " | EVS PCA variance profiles"), gp = gpar(fontface = "bold", cex = 1.02)), bottom = textGrob("Top: normalized before EVS. Bottom: raw before EVS. Left: treatment. Right: control.", gp = gpar(cex = 0.86))), paste0(comparison_name, ": EVS variance panel"))
    save_grob(variance_panel, file.path(cmp_dir, "EVS_PCA_variance_profiles_combined.png"), width = 18, height = 13)
  }

  dataset_list <- list(raw_dataset = evs$raw_dataset, leading_edge_dataset = evs$leading_edge_dataset, remainder_dataset = evs$remainder_dataset)
  normalized_dataset_list <- list(
    raw_dataset = subset_normalized_counts_by_ids(evs$normalized_counts, rownames(evs$normalized_counts), "raw_dataset", comparison_name, strict = TRUE),
    leading_edge_dataset = subset_normalized_counts_by_ids(evs$normalized_counts, rownames(evs$leading_edge_dataset), "leading_edge_dataset", comparison_name, strict = TRUE),
    remainder_dataset = subset_normalized_counts_by_ids(evs$normalized_counts, rownames(evs$remainder_dataset), "remainder_dataset", comparison_name, strict = TRUE)
  )
  dataset_fig_dirs <- list(raw_dataset = fig_raw_dir, leading_edge_dataset = fig_le_dir, remainder_dataset = fig_rem_dir)

  analysis_results <- list()
  for (nm in names(dataset_list)) {
    full_dataset_name <- paste(comparison_name, nm, sep = "_")
    fit <- run_core_analysis(count_mat = dataset_list[[nm]], coldata = coldata, dataset_name = full_dataset_name, annot_df = annot_df)
    df <- fit$results
    fig_subdir <- dataset_fig_dirs[[nm]]

    save_csv(df, file.path(tab_dir, paste0(full_dataset_name, "_results_full.csv")))
    save_csv(subset(df, standard_significant), file.path(tab_dir, paste0(full_dataset_name, "_standard_significant.csv")))
    save_csv(subset(df, HBFSS_significant), file.path(tab_dir, paste0(full_dataset_name, "_HBFSS_significant.csv")))
    save_csv(subset(df, effect_class == "strong_effect"), file.path(tab_dir, paste0(full_dataset_name, "_strong_effect.csv")))
    save_csv(subset(df, effect_class == "weak_effect"), file.path(tab_dir, paste0(full_dataset_name, "_weak_effect.csv")))

    summary_row <- data.frame(
      comparison_name = comparison_name,
      dataset_name = full_dataset_name,
      n_features = nrow(df),
      hc_empirical_threshold = fit$hc_empirical_threshold,
      hc_wald_threshold = fit$hc_wald_threshold,
      hbfss_threshold = fit$hbfss_threshold,
      n_standard_significant = sum(df$standard_significant, na.rm = TRUE),
      n_HBFSS_significant = sum(df$HBFSS_significant, na.rm = TRUE),
      n_overlap_significant = sum(df$standard_significant & df$HBFSS_significant, na.rm = TRUE),
      n_strong_effect = sum(df$effect_class == "strong_effect", na.rm = TRUE),
      n_weak_effect = sum(df$effect_class == "weak_effect", na.rm = TRUE),
      evs_fixed_top_n = evs_fixed_top_n,
      shared_cutoff_rank = evs$combined_cutoff_info$top_n_actual,
      shared_cutoff_method = evs$combined_cutoff_info$method,
      shared_cutoff_quantile = evs$fit_trt$cutoff_quantile,
      shared_selected_reason = evs$combined_cutoff_info$selected_reason,
      stringsAsFactors = FALSE
    )
    save_csv(summary_row, file.path(tab_dir, paste0(full_dataset_name, "_summary.csv")))

    save_grob(plot_standard_volcano(df, full_dataset_name), file.path(fig_subdir, paste0(full_dataset_name, "_standard_volcano.png")), width = 8.8, height = 6.5)
    save_grob(plot_hbfss_volcano(df, full_dataset_name), file.path(fig_subdir, paste0(full_dataset_name, "_HBFSS_volcano.png")), width = 8.8, height = 6.5)
    save_grob(plot_publication_volcano_panel(df, full_dataset_name), file.path(fig_subdir, paste0(full_dataset_name, "_publication_volcano.png")), width = 8.8, height = 6.5)
    save_grob(plot_dispersion_panel_for_dataset(df, full_dataset_name), file.path(fig_subdir, paste0(full_dataset_name, "_dispersion.png")), width = 8.6, height = 6.2)

    analysis_results[[nm]] <- list(dds = fit$dds, results = df, summary = summary_row, dataset_mat = dataset_list[[nm]], fig_subdir = fig_subdir)
  }

  build_pc1_feature_export(analysis_results, annot_df, comparison_name, evs, tab_dir)
  save_cross_dataset_comparison_panels(comparison_name, analysis_results, cmp_dir, dataset_list, coldata, normalized_dataset_list = normalized_dataset_list)

  if (run_twas_overlap) {
    if (!file.exists(twas_file)) stop(sprintf("TWAS file not found: %s", twas_file), call. = FALSE)
    TWAS_Seq <- read.csv(twas_file, header = TRUE, stringsAsFactors = FALSE)
    if (ncol(TWAS_Seq) < 4) stop("TWAS file must contain at least 4 columns.", call. = FALSE)
    TWAS_data <- TWAS_Seq[, c(1, 4), drop = FALSE]
    colnames(TWAS_data) <- c("source_id", "gene_symbol")
    twas_genes <- clean_gene_set(TWAS_data$gene_symbol)

    get_twas_overlap <- function(result_df, dataset_nm, cmp_name, out_dir, mode = c("union", "standard_only", "HBFSS_only")) {
      mode <- match.arg(mode)
      sig_df <- switch(mode,
        union = subset(result_df, standard_significant | HBFSS_significant),
        standard_only = subset(result_df, standard_significant),
        HBFSS_only = subset(result_df, HBFSS_significant)
      )
      sig_df$gene_symbol_clean <- tolower(trimws(sig_df$gene_symbol))
      overlap_df <- subset(sig_df, gene_symbol_clean %in% twas_genes)
      summary_df <- data.frame(comparison_name = cmp_name, dataset_name = dataset_nm, selection_mode = mode, n_selected_features = nrow(sig_df), n_overlap_features = nrow(overlap_df), n_overlap_genes = length(unique(overlap_df$gene_symbol_clean)), stringsAsFactors = FALSE)
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
  if (!length(sm_list)) return(data.frame())
  dplyr::bind_rows(sm_list)
}

# -----------------------------------------------------------------------------
# RUN ALL COMPARISONS
# -----------------------------------------------------------------------------
comparison_inputs <- lapply(seq_len(nrow(comparison_table)), function(i) {
  prepare_comparison_data(
    comparison_name = comparison_table$comparison_name[i],
    group1_prefix = comparison_table$group1_prefix[i],
    group2_prefix = comparison_table$group2_prefix[i],
    WTTS_Seq = WTTS_Seq,
    meta_all = meta_all
  )
})
names(comparison_inputs) <- comparison_table$comparison_name

all_summaries_list <- list()
failed_comparisons <- list()
run_start <- Sys.time()

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
      annot_df = OrigID_Symbol
    ),
    error = function(e) {
      failed_comparisons[[cmp]] <<- data.frame(comparison_name = cmp, error_message = conditionMessage(e), stringsAsFactors = FALSE)
      NULL
    }
  )

  if (!is.null(out) && nrow(out) > 0) all_summaries_list[[cmp]] <- out
}

all_summaries <- if (length(all_summaries_list) > 0) dplyr::bind_rows(all_summaries_list) else data.frame()
if (nrow(all_summaries) > 0) save_csv(all_summaries, file.path(output_dir, "all_comparisons_overall_summary.csv"))
if (length(failed_comparisons) > 0) save_csv(dplyr::bind_rows(failed_comparisons), file.path(output_dir, "failed_comparisons.csv"))

runtime_seconds <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
run_manifest <- data.frame(
  analysis_name = analysis_name,
  runtime_readable = format_runtime_minutes(runtime_seconds),
  runtime_seconds = runtime_seconds,
  alpha_level = alpha_level,
  lfc_boundary = lfc_boundary,
  evs_fixed_top_n = evs_fixed_top_n,
  output_dir = output_dir,
  stringsAsFactors = FALSE
)
save_csv(run_manifest, file.path(output_dir, "run_manifest.csv"))

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
