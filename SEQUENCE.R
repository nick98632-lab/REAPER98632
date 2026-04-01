
# =============================================================================
# SEQUENCE MANUSCRIPT PIPELINE
# FULL REPAIRED SCRIPT
# =============================================================================
# PURPOSE
# This script implements the final manuscript analysis workflow for the
# SEQUENCE WTTS Seq dataset.
#
# PRIMARY ANALYSIS RULES
# 1. The EVS cutoff is chosen by the shared Fourier regime crossing method.
# 2. The only backup cutoff is a user supplied manual fixed rank.
# 3. DESeq2 significance uses Benjamini Hochberg adjusted p value < 0.20.
# 4. Volcano plots must reflect the actual analysis classes exactly.
# 5. Cross comparison exports are written to one stable repository folder.
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

# =============================================================================
# USER INPUT AND REPOSITORY PATHS
# =============================================================================

repo_dir <- getwd()
input_dir <- file.path(repo_dir, "data")
output_root <- file.path(repo_dir, "exports")

analysis_name <- "EVS_HBFSS_AllComparisons_Output"
output_dir <- file.path(output_root, analysis_name)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

figure_dir <- file.path(output_dir, "figures")
table_dir  <- file.path(output_dir, "tables")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir,  showWarnings = FALSE, recursive = TRUE)

count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")

run_twas_overlap <- FALSE
twas_file <- file.path(input_dir, "3aTWAS_genes_of_11_brain_disorders.csv")

# =============================================================================
# PRIMARY SIGNIFICANCE AND CUTOFF SETTINGS
# =============================================================================

alpha_level <- 0.20
lfc_boundary <- 1.0

evs_cutoff_mode_main   <- "fourier_shared_cutoff"
evs_backup_cutoff_mode <- "fixed_top_n"

evs_fixed_top_n <- 5000L
evs_fixed_rank_override <- NA_integer_
evs_use_manual_rank_first <- FALSE

crossing_min_percentile <- 0.00
crossing_max_percentile <- 0.99
crossing_rule <- "first_left_crossing"

use_group_crossing_bracket <- TRUE
group_bracket_fallback <- "midpoint"

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

crossing_stability_window_n <- 3L
crossing_plot_line_width <- 0.95
crossing_plot_vline_width <- 0.95
crossing_plot_hline_width <- 0.65
crossing_plot_point_size <- 1.8

fdr_clip_floor <- 1e-300
fdr_clip_ceiling <- 0.99
hc_threshold_upper_cap <- 0.99
fdrtool_pct0 <- 0.75
max_usable_hc_p_threshold <- 0.95

# =============================================================================
# FIGURE EXPORT CONTROLS
# =============================================================================

export_optional_mean_histograms        <- FALSE
export_optional_empirical_p_histograms <- FALSE
export_optional_hbfss_distributions    <- FALSE
export_optional_evs_raw_histograms     <- FALSE
export_optional_evs_variance_profiles  <- FALSE

# =============================================================================
# PLOT CONSTANTS
# =============================================================================

figure_dpi <- 320
base_theme_size <- 10

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

# =============================================================================
# EMBEDDED METADATA
# =============================================================================

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
  comparison_name  = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix    = c("R0", "R2", "R4", "R8"),
  group2_prefix    = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

dataset_key_order <- c("original", "leading_edge", "remainder")

# =============================================================================
# LABELS AND DISPLAY DICTIONARIES
# =============================================================================

condition_labels <- c(
  "untrt" = "Control",
  "trt"   = "Treatment"
)

condition_shapes <- c(
  "untrt" = 21,
  "trt"   = 24
)

condition_fills <- c(
  "untrt" = plot_palette$control,
  "trt"   = plot_palette$treatment
)

dataset_pretty_labels <- c(
  "original"     = "Original dataset",
  "leading_edge" = "Leading edge dataset",
  "remainder"    = "Remainder dataset"
)

pretty_group_label <- function(x) {
  x <- as.character(x)[1]
  dplyr::case_when(
    identical(x, "treatment") ~ "Treatment",
    identical(x, "control")   ~ "Control",
    identical(x, "trt")       ~ "Treatment",
    identical(x, "untrt")     ~ "Control",
    TRUE ~ x
  )
}

pretty_dataset_label <- function(x) {
  x <- as.character(x)[1]
  if (x %in% names(dataset_pretty_labels)) {
    return(dataset_pretty_labels[[x]])
  }
  x
}

preprocessing_panel_label <- function(x) {
  x <- as.character(x)[1]
  if (is.na(x) || !nzchar(x)) return("Preprocessing not specified")
  dplyr::case_when(
    x %in% c(
      "DESeq2-normalized counts",
      "normalized",
      "Normalized prior to eigenvector splitting",
      "Normalized before EVS"
    ) ~ "Normalized prior to eigenvector splitting",
    x %in% c(
      "raw counts without DESeq2 normalization",
      "raw_counts",
      "Eigenvector splitting without prior normalization",
      "Raw counts before EVS"
    ) ~ "Eigenvector splitting without prior normalization",
    TRUE ~ x
  )
}

# =============================================================================
# GENERAL HELPER FUNCTIONS
# =============================================================================

assert_required_columns <- function(df, required_cols, object_name = "data frame") {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns in ",
        object_name,
        ": ",
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

compact_title <- function(x, width = 58) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

compact_caption <- function(x, width = 120) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

clip_probabilities <- function(x,
                               eps = fdr_clip_floor,
                               upper = fdr_clip_ceiling) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  x <- pmax(x, eps, na.rm = FALSE)
  x <- pmin(x, upper, na.rm = FALSE)
  x
}

make_design_formula <- function() {
  as.formula("~ condition")
}

save_csv <- function(df, file_path, row.names = FALSE) {
  utils::write.csv(df, file = file_path, row.names = row.names)
}

save_grob <- function(grob_obj,
                      file_path,
                      width,
                      height,
                      dpi = figure_dpi,
                      bg = "white") {
  if (is.null(grob_obj)) return(invisible(NULL))
  ggplot2::ggsave(
    filename = file_path,
    plot = grob_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )
  invisible(file_path)
}

safe_plot_build <- function(expr, label = "plot") {
  tryCatch(
    expr,
    error = function(e) {
      warning(sprintf("%s failed: %s", label, conditionMessage(e)), call. = FALSE)
      NULL
    }
  )
}

manuscript_theme <- function(base_size = base_theme_size) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", linewidth = 0.5),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.title = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = rel(0.95)),
      plot.caption = element_text(size = rel(0.85), hjust = 0),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      legend.key = element_rect(fill = "white", colour = NA),
      strip.background = element_rect(fill = "grey95", colour = "grey50"),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(t = 10, r = 14, b = 12, l = 12)
    )
}

format_runtime_minutes <- function(seconds_value) {
  if (!is.finite(seconds_value) || is.na(seconds_value) || seconds_value < 0) {
    return("unknown")
  }
  if (seconds_value < 60) {
    return(sprintf("~%d sec", as.integer(round(seconds_value))))
  }
  sprintf("~%.1f min", seconds_value / 60)
}

# =============================================================================
# MEAN EXPRESSION AND BASIC QC HELPERS
# =============================================================================

compute_mean_expression_table <- function(raw_counts, coldata) {
  sample_ids <- colnames(raw_counts)
  trt_ids    <- sample_ids[coldata$condition == "trt"]
  untrt_ids  <- sample_ids[coldata$condition == "untrt"]

  data.frame(
    feature_id       = as.character(rownames(raw_counts)),
    mean_trt         = rowMeans(raw_counts[, trt_ids, drop = FALSE]),
    mean_untrt       = rowMeans(raw_counts[, untrt_ids, drop = FALSE]),
    mean_all         = rowMeans(raw_counts),
    stringsAsFactors = FALSE
  )
}

plot_mean_histogram_panel <- function(df_means, base_mean_vec, dataset_name) {
  make_hist <- function(vals, panel_title) {
    ggplot(data.frame(x = safe_log10(vals + 1)), aes(x)) +
      geom_histogram(
        bins = HIST_BINS,
        fill = plot_palette$histogram,
        color = HIST_COLOR
      ) +
      labs(
        title = paste(dataset_name, panel_title),
        x = "log10(mean + 1)",
        y = "Count"
      ) +
      manuscript_theme()
  }

  arrangeGrob(
    make_hist(df_means$mean_trt,   "Treatment mean"),
    make_hist(df_means$mean_untrt, "Control mean"),
    make_hist(df_means$mean_all,   "Pooled mean"),
    make_hist(base_mean_vec,       "DESeq2 baseMean"),
    ncol = 2,
    top = textGrob(
      paste(dataset_name, "Mean expression histograms"),
      gp = gpar(fontface = "bold", cex = 1.2)
    )
  )
}

# =============================================================================
# EMPIRICAL NULL AND HC HELPERS
# =============================================================================

run_empirical_null_fdrtool <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  stat_vec <- unname(stat_vec)

  if (length(stat_vec) < 5) {
    stop(
      sprintf("[%s] Fewer than 5 finite Wald statistics were available for fdrtool.", dataset_name),
      call. = FALSE
    )
  }

  fit <- tryCatch(
    fdrtool(
      stat_vec,
      statistic = "normal",
      plot = FALSE,
      verbose = FALSE,
      cutoff.method = "fndr",
      pct0 = fdrtool_pct0
    ),
    error = function(e1) {
      message(sprintf("[%s] Primary fdrtool call failed: %s", dataset_name, conditionMessage(e1)))
      tryCatch(
        fdrtool(
          as.vector(stat_vec),
          statistic = "normal",
          plot = FALSE,
          verbose = FALSE,
          cutoff.method = "pct0",
          pct0 = fdrtool_pct0
        ),
        error = function(e2) {
          stop(
            sprintf("[%s] fdrtool failed after retry: %s", dataset_name, conditionMessage(e2)),
            call. = FALSE
          )
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

  if (length(sorted_empirical_p) < 5) {
    return(NA_real_)
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

# =============================================================================
# MANUAL BACKUP CUTOFF HELPER
# =============================================================================

resolve_top_n_cutoff <- function(sorted_values_desc, top_n = evs_fixed_top_n) {
  n_total <- length(sorted_values_desc)

  if (n_total == 0) {
    stop("resolve_top_n_cutoff() received an empty vector.", call. = FALSE)
  }

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

# =============================================================================
# PCA HELPER
# =============================================================================

compute_dataset_pca_plot <- function(count_df,
                                     coldata,
                                     dataset_name,
                                     preprocessing = c("normalized", "raw_counts"),
                                     precomputed_matrix = NULL) {
  preprocessing <- match.arg(preprocessing)
  design_formula <- make_design_formula()

  if (!is.null(precomputed_matrix)) {
    x <- as.matrix(precomputed_matrix)
  } else if (preprocessing == "normalized") {
    dds <- DESeqDataSetFromMatrix(
      countData = count_df,
      colData = coldata,
      design = design_formula
    )
    dds <- dds[rowSums(counts(dds)) > 0, ]
    dds <- estimateSizeFactors(dds)
    x <- counts(dds, normalized = TRUE)
  } else {
    x <- as.matrix(count_df)
  }

  pca_fit <- prcomp(t(x), scale. = FALSE, rank. = 2)
  pca_var <- pca_fit$sdev ^ 2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)

  pca_df <- data.frame(
    Sample = rownames(pca_fit$x),
    PC1 = pca_fit$x[, 1],
    PC2 = pca_fit$x[, 2],
    Condition = as.character(coldata[rownames(pca_fit$x), "condition"]),
    stringsAsFactors = FALSE
  )

  ggplot(pca_df, aes(PC1, PC2, label = Sample, shape = Condition, fill = Condition)) +
    geom_hline(
      yintercept = 0,
      linewidth = LINE_WIDTH_ZERO,
      linetype = "dashed",
      colour = "grey70"
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = LINE_WIDTH_ZERO,
      linetype = "dashed",
      colour = "grey70"
    ) +
    geom_point(size = 3, colour = "white", stroke = POINT_STROKE + 0.15) +
    geom_text_repel(
      size = 2.0,
      max.overlaps = 8,
      force = 1.0,
      box.padding = 0.22,
      point.padding = 0.10,
      min.segment.length = 0
    ) +
    scale_shape_manual(
      values = condition_shapes,
      labels = condition_labels,
      name = "Condition",
      guide = guide_legend(
        override.aes = list(
          size = 3.0,
          fill = unname(condition_fills),
          colour = "white"
        )
      )
    ) +
    scale_fill_manual(
      values = condition_fills,
      labels = condition_labels,
      name = "Condition",
      guide = "none"
    ) +
    labs(
      title = compact_title(pretty_dataset_label(dataset_name)),
      subtitle = preprocessing_panel_label(preprocessing),
      x = paste0("PC1 (", pca_var_per[1], "%)"),
      y = paste0("PC2 (", pca_var_per[2], "%)")
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 12, r = 18, b = 16, l = 16)
    )
}

# =============================================================================
# PART 2
# =============================================================================

estimate_feature_metrics_for_loading <- function(count_df, coldata) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_df)),
    colData   = coldata,
    design    = ~ 1
  )

  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds, quiet = TRUE)
  dds <- estimateDispersionsFit(dds, quiet = TRUE)

  md <- as.data.frame(SummarizedExperiment::mcols(dds), stringsAsFactors = FALSE)
  md$feature_id <- rownames(md)

  keep_cols <- intersect(
    c("feature_id", "baseMean", "dispGeneEst", "dispFit", "dispersion"),
    colnames(md)
  )

  md <- md[, keep_cols, drop = FALSE]
  md <- md[!duplicated(md$feature_id), , drop = FALSE]
  md
}

build_ranked_fourier_metric_table <- function(loading_tbl,
                                              mean_col = "baseMean",
                                              dispersion_col = "dispGeneEst") {
  needed <- c("feature_id", "rank", "pc1_loading_abs", mean_col, dispersion_col)
  assert_required_columns(loading_tbl, needed, "loading_tbl")

  df <- loading_tbl[, needed, drop = FALSE]
  colnames(df)[colnames(df) == mean_col] <- "mean_value"
  colnames(df)[colnames(df) == dispersion_col] <- "dispersion_value"

  df <- df[
    is.finite(df$mean_value) &
      !is.na(df$mean_value) &
      df$mean_value > 0 &
      is.finite(df$dispersion_value) &
      !is.na(df$dispersion_value) &
      df$dispersion_value > 0,
    ,
    drop = FALSE
  ]

  if (nrow(df) < 25) {
    return(NULL)
  }

  df <- df[order(df$rank), , drop = FALSE]

  mu    <- pmax(df$mean_value, 1e-12)
  alpha <- pmax(df$dispersion_value, 1e-12)

  df$iod_nb     <- 1 + alpha * mu
  df$cv2_nb     <- (1 / mu) + alpha
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

  max_rank <- max(1L, floor(n_total * max_rank_frac))
  keep <- center_ranks >= min_rank & center_ranks <= max_rank

  pct_grid <- pct_grid[keep]
  center_ranks <- center_ranks[keep]

  half_window <- max(5L, round((window_fraction * n_total) / 2))

  data.frame(
    percentile  = pct_grid,
    center_rank = center_ranks,
    lo_rank     = pmax(1L, center_ranks - half_window),
    hi_rank     = pmin(n_total, center_ranks + half_window),
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
  fm  <- stats::as.formula(paste("y ~", rhs))

  fit <- tryCatch(
    stats::lm(fm, data = dd),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  fitted_y   <- as.numeric(stats::predict(fit, newdata = dd))
  residual_y <- dd$y - fitted_y

  local_amplitude <- 0.5 * (max(fitted_y, na.rm = TRUE) - min(fitted_y, na.rm = TRUE))
  peak_idx   <- which.max(fitted_y)
  trough_idx <- which.min(fitted_y)
  center_rank <- mean(c(rank_vec[peak_idx], rank_vec[trough_idx]))

  list(
    fit             = fit,
    fitted_y        = fitted_y,
    residual_y      = residual_y,
    local_amplitude = local_amplitude,
    peak_rank       = rank_vec[peak_idx],
    trough_rank     = rank_vec[trough_idx],
    center_rank     = center_rank,
    residual_sd     = stats::sd(residual_y, na.rm = TRUE)
  )
}

summarize_local_fourier_window <- function(metric_df, lo_rank, hi_rank) {
  sub <- metric_df[metric_df$rank >= lo_rank & metric_df$rank <= hi_rank, , drop = FALSE]

  if (nrow(sub) < 9L) {
    return(data.frame(
      iod_amplitude   = NA_real_,
      cv2_amplitude   = NA_real_,
      iod_center      = NA_real_,
      cv2_center      = NA_real_,
      iod_residual_sd = NA_real_,
      cv2_residual_sd = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  iod_fit <- fit_local_fourier(sub$rank, sub$log_iod_nb)
  cv2_fit <- fit_local_fourier(sub$rank, sub$log_cv2_nb)

  if (is.null(iod_fit) || is.null(cv2_fit)) {
    return(data.frame(
      iod_amplitude   = NA_real_,
      cv2_amplitude   = NA_real_,
      iod_center      = NA_real_,
      cv2_center      = NA_real_,
      iod_residual_sd = NA_real_,
      cv2_residual_sd = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    iod_amplitude   = iod_fit$local_amplitude,
    cv2_amplitude   = cv2_fit$local_amplitude,
    iod_center      = iod_fit$center_rank,
    cv2_center      = cv2_fit$center_rank,
    iod_residual_sd = iod_fit$residual_sd,
    cv2_residual_sd = cv2_fit$residual_sd,
    stringsAsFactors = FALSE
  )
}

build_local_fourier_wave_map <- function(loading_tbl,
                                         mean_col = "baseMean",
                                         dispersion_col = "dispGeneEst") {
  metric_df <- build_ranked_fourier_metric_table(
    loading_tbl     = loading_tbl,
    mean_col        = mean_col,
    dispersion_col  = dispersion_col
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
    scales::rescale(
      wave_map$iod_amplitude,
      to = c(0, 1),
      from = range(wave_map$iod_amplitude, na.rm = TRUE)
    )
  }

  cv2_amp_scaled <- if (all(is.na(wave_map$cv2_amplitude))) {
    rep(NA_real_, nrow(wave_map))
  } else {
    scales::rescale(
      wave_map$cv2_amplitude,
      to = c(0, 1),
      from = range(wave_map$cv2_amplitude, na.rm = TRUE)
    )
  }

  wave_map$center_distance  <- abs(wave_map$iod_center - wave_map$cv2_center)
  wave_map$center_agreement <- 1 / (1 + wave_map$center_distance)
  wave_map$local_fourier_score <- (
    fourier_score_weight_iod_amp * iod_amp_scaled +
      fourier_score_weight_cv2_amp * cv2_amp_scaled +
      fourier_score_weight_center_agreement * wave_map$center_agreement
  )

  list(
    metric_df = metric_df,
    wave_map  = wave_map
  )
}

combine_treatment_control_fourier_maps <- function(trt_wave_obj, ctrl_wave_obj) {
  if (is.null(trt_wave_obj) || is.null(ctrl_wave_obj)) return(NULL)

  trt_map  <- trt_wave_obj$wave_map
  ctrl_map <- ctrl_wave_obj$wave_map

  if (is.null(trt_map) || is.null(ctrl_map) || !nrow(trt_map) || !nrow(ctrl_map)) {
    return(NULL)
  }

  keep_cols <- c(
    "percentile", "center_rank", "iod_amplitude", "cv2_amplitude",
    "iod_center", "cv2_center", "iod_residual_sd", "cv2_residual_sd",
    "center_agreement", "local_fourier_score"
  )

  trt_map  <- trt_map[, keep_cols, drop = FALSE]
  ctrl_map <- ctrl_map[, keep_cols, drop = FALSE]

  names(trt_map)  <- paste0("trt_",  names(trt_map))
  names(ctrl_map) <- paste0("ctrl_", names(ctrl_map))

  names(trt_map)[names(trt_map) == "trt_percentile"] <- "percentile"
  names(ctrl_map)[names(ctrl_map) == "ctrl_percentile"] <- "percentile"

  df <- dplyr::inner_join(trt_map, ctrl_map, by = "percentile")
  if (!nrow(df)) return(NULL)

  df$center_rank <- round((df$trt_center_rank + df$ctrl_center_rank) / 2)

  df$combined_iod_amplitude <- rowMeans(
    cbind(df$trt_iod_amplitude, df$ctrl_iod_amplitude),
    na.rm = TRUE
  )

  df$combined_cv2_amplitude <- rowMeans(
    cbind(df$trt_cv2_amplitude, df$ctrl_cv2_amplitude),
    na.rm = TRUE
  )

  df$combined_center_agreement <- rowMeans(
    cbind(df$trt_center_agreement, df$ctrl_center_agreement),
    na.rm = TRUE
  )

  df$combined_fourier_score <- (
    fourier_score_weight_iod_amp * scales::rescale(
      df$combined_iod_amplitude,
      to = c(0, 1),
      from = range(df$combined_iod_amplitude, na.rm = TRUE)
    ) +
      fourier_score_weight_cv2_amp * scales::rescale(
        df$combined_cv2_amplitude,
        to = c(0, 1),
        from = range(df$combined_cv2_amplitude, na.rm = TRUE)
      ) +
      fourier_score_weight_center_agreement * df$combined_center_agreement
  )

  df$regime_difference <- df$combined_iod_amplitude - df$combined_cv2_amplitude
  df
}

find_regime_difference_crossings <- function(diff_df) {
  stopifnot("regime_difference" %in% names(diff_df))
  stopifnot("percentile" %in% names(diff_df))
  stopifnot("center_rank" %in% names(diff_df))

  rd <- diff_df$regime_difference
  pp <- diff_df$percentile
  rr <- diff_df$center_rank

  crossings <- list()
  kk <- 1L

  for (i in seq_len(nrow(diff_df) - 1L)) {
    left_val  <- rd[i]
    right_val <- rd[i + 1L]

    if (!is.finite(left_val) || !is.finite(right_val)) {
      next
    }

    is_crossing <- (left_val == 0) || (right_val == 0) || (left_val * right_val < 0)
    if (!is_crossing) next

    if (identical(crossing_rule, "first_left_crossing")) {
      crossing_pct  <- pp[i]
      crossing_rank <- rr[i]
    } else {
      crossing_pct  <- mean(c(pp[i], pp[i + 1L]))
      crossing_rank <- as.integer(round(mean(c(rr[i], rr[i + 1L]))))
    }

    crossings[[kk]] <- data.frame(
      crossing_id             = paste0("crossing_", kk),
      index_left              = i,
      index_right             = i + 1L,
      percentile_left         = pp[i],
      percentile_right        = pp[i + 1L],
      crossing_percentile     = crossing_pct,
      crossing_rank           = crossing_rank,
      regime_difference_left  = left_val,
      regime_difference_right = right_val,
      stringsAsFactors = FALSE
    )

    kk <- kk + 1L
  }

  if (!length(crossings)) {
    return(data.frame())
  }

  dplyr::bind_rows(crossings)
}

flag_stable_crossings <- function(diff_df,
                                  crossing_tbl,
                                  stability_window_n = crossing_stability_window_n) {
  if (is.null(crossing_tbl) || !nrow(crossing_tbl)) {
    return(crossing_tbl)
  }

  rd <- diff_df$regime_difference
  crossing_tbl$stable_crossing <- FALSE

  for (i in seq_len(nrow(crossing_tbl))) {
    left_idx  <- crossing_tbl$index_left[i]
    right_idx <- crossing_tbl$index_right[i]

    left_window_lo  <- max(1L, left_idx  - stability_window_n + 1L)
    left_window_hi  <- left_idx
    right_window_lo <- right_idx
    right_window_hi <- min(length(rd), right_idx + stability_window_n - 1L)

    left_window  <- rd[left_window_lo:left_window_hi]
    right_window <- rd[right_window_lo:right_window_hi]

    left_ok  <- all(is.finite(left_window))  && all(left_window  >= 0)
    right_ok <- all(is.finite(right_window)) && all(right_window <= 0)

    crossing_tbl$stable_crossing[i] <- isTRUE(left_ok && right_ok)
  }

  crossing_tbl
}

select_regime_shift_crossing <- function(combined_wave_df) {
  diff_df <- combined_wave_df[order(combined_wave_df$percentile), , drop = FALSE]
  crossing_tbl <- find_regime_difference_crossings(diff_df)

  if (is.null(crossing_tbl) || !nrow(crossing_tbl)) {
    return(list(
      diff_df            = diff_df,
      crossing_table     = data.frame(),
      selected_crossing  = NULL,
      selected_reason    = "no_crossings_found"
    ))
  }

  crossing_tbl <- flag_stable_crossings(diff_df, crossing_tbl)

  stable_tbl <- crossing_tbl[crossing_tbl$stable_crossing, , drop = FALSE]
  if (nrow(stable_tbl)) {
    selected <- stable_tbl[order(stable_tbl$crossing_percentile), , drop = FALSE][1, , drop = FALSE]
    reason <- "first_stable_crossing"
  } else {
    selected <- crossing_tbl[order(crossing_tbl$crossing_percentile), , drop = FALSE][1, , drop = FALSE]
    reason <- "first_crossing_fallback"
  }

  list(
    diff_df           = diff_df,
    crossing_table    = crossing_tbl,
    selected_crossing = selected,
    selected_reason   = reason
  )
}

resolve_combined_fourier_cutoff <- function(fit_trt_loading_tbl,
                                            fit_ctrl_loading_tbl,
                                            fixed_top_n = evs_fixed_top_n) {
  n_total  <- nrow(fit_trt_loading_tbl)
  fallback <- resolve_top_n_cutoff(
    fit_trt_loading_tbl$pc1_loading_abs,
    top_n = fixed_top_n
  )

  trt_wave_obj  <- build_local_fourier_wave_map(fit_trt_loading_tbl)
  ctrl_wave_obj <- build_local_fourier_wave_map(fit_ctrl_loading_tbl)

  if (is.null(trt_wave_obj) || is.null(ctrl_wave_obj)) {
    fallback$method                <- "fixed_top_n_fallback"
    fallback$trt_wave_obj          <- trt_wave_obj
    fallback$ctrl_wave_obj         <- ctrl_wave_obj
    fallback$combined_wave_map     <- data.frame()
    fallback$combined_interval_table <- data.frame()
    fallback$crossing_table        <- data.frame()
    fallback$selected_crossing     <- NULL
    fallback$selected_reason       <- "combined_wave_map_missing"
    return(fallback)
  }

  combined_wave_df <- combine_treatment_control_fourier_maps(trt_wave_obj, ctrl_wave_obj)
  if (is.null(combined_wave_df) || !nrow(combined_wave_df)) {
    fallback$method                <- "fixed_top_n_fallback"
    fallback$trt_wave_obj          <- trt_wave_obj
    fallback$ctrl_wave_obj         <- ctrl_wave_obj
    fallback$combined_wave_map     <- data.frame()
    fallback$combined_interval_table <- data.frame()
    fallback$crossing_table        <- data.frame()
    fallback$selected_crossing     <- NULL
    fallback$selected_reason       <- "combined_wave_map_missing"
    return(fallback)
  }

  crossing_info <- select_regime_shift_crossing(combined_wave_df)
  selected_crossing <- crossing_info$selected_crossing

  if (is.null(selected_crossing) || !nrow(selected_crossing)) {
    fallback$method                <- "fixed_top_n_fallback"
    fallback$trt_wave_obj          <- trt_wave_obj
    fallback$ctrl_wave_obj         <- ctrl_wave_obj
    fallback$combined_wave_map     <- crossing_info$diff_df
    fallback$combined_interval_table <- data.frame()
    fallback$crossing_table        <- crossing_info$crossing_table
    fallback$selected_crossing     <- NULL
    fallback$selected_reason       <- crossing_info$selected_reason
    return(fallback)
  }

  selected_rank <- as.integer(selected_crossing$crossing_rank[1])
  selected_rank <- min(max(1L, selected_rank), n_total)

  cutoff_value    <- fit_trt_loading_tbl$pc1_loading_abs[selected_rank]
  cutoff_quantile <- 1 - (selected_rank / n_total)

  candidate_table <- data.frame(
    candidate_id     = selected_crossing$crossing_id[1],
    rank_index       = selected_rank,
    percentile       = selected_crossing$crossing_percentile[1],
    cutoff_value     = cutoff_value,
    cutoff_quantile  = cutoff_quantile,
    selected         = TRUE,
    selected_reason  = crossing_info$selected_reason,
    stringsAsFactors = FALSE
  )

  list(
    top_n_actual          = selected_rank,
    cutoff_value          = cutoff_value,
    cutoff_quantile       = cutoff_quantile,
    n_total               = n_total,
    method                = "first_stable_crossing",
    trt_wave_obj          = trt_wave_obj,
    ctrl_wave_obj         = ctrl_wave_obj,
    combined_wave_map     = crossing_info$diff_df,
    combined_interval_table = data.frame(),
    crossing_table        = crossing_info$crossing_table,
    selected_crossing     = selected_crossing,
    candidate_table       = candidate_table,
    selected_reason       = crossing_info$selected_reason
  )
}

plot_combined_fourier_wave_map <- function(combined_cutoff_info, comparison_name) {
  df <- combined_cutoff_info$combined_wave_map
  if (is.null(df) || !nrow(df)) return(NULL)

  sc <- combined_cutoff_info$selected_crossing
  crossing_pct  <- if (!is.null(sc) && nrow(sc)) sc$crossing_percentile[1] else NA_real_
  crossing_rank <- if (!is.null(sc) && nrow(sc)) sc$crossing_rank[1] else NA_integer_

  ymax <- max(c(df$combined_iod_amplitude, df$combined_cv2_amplitude), na.rm = TRUE)

  ggplot(df, aes(percentile)) +
    geom_line(aes(y = combined_iod_amplitude, color = "Composite IOD"), linewidth = crossing_plot_line_width) +
    geom_line(aes(y = combined_cv2_amplitude, color = "Composite CV²"), linewidth = crossing_plot_line_width) +
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
        "Cutoff percentile = ", signif(crossing_pct, 4),
        "\nCutoff rank = ", crossing_rank
      ),
      fill = "white",
      colour = plot_palette$threshold,
      size = 3.0,
      label.size = 0.15,
      vjust = -0.5
    ) +
    scale_color_manual(
      values = c(
        "Composite IOD" = plot_palette$treatment,
        "Composite CV²" = plot_palette$control
      )
    ) +
    labs(
      title = paste0(comparison_name, " | two-line regime crossing"),
      subtitle = compact_caption(
        "The dashed vertical line marks the first stable crossing between the composite local IOD and composite local CV² lines. This crossing defines the shared EVS cutoff.",
        width = 92
      ),
      x = "Percentile center",
      y = "Composite local amplitude",
      color = NULL
    ) +
    manuscript_theme()
}

plot_combined_fourier_score <- function(combined_cutoff_info, comparison_name) {
  df <- combined_cutoff_info$combined_wave_map
  cand_df <- combined_cutoff_info$candidate_table

  if (is.null(df) || !nrow(df) || !"combined_fourier_score" %in% names(df)) {
    return(NULL)
  }

  p <- ggplot(df, aes(percentile, combined_fourier_score)) +
    geom_col(width = 0.008, fill = plot_palette$threshold)

  if (!is.null(cand_df) && nrow(cand_df)) {
    chosen_pct <- cand_df$percentile[1]
    p <- p + geom_vline(
      xintercept = chosen_pct,
      color = plot_palette$threshold,
      linewidth = 0.9
    )
  }

  p +
    labs(
      title = paste0(comparison_name, " | descriptive score profile"),
      subtitle = compact_caption(
        "This score curve is descriptive only. The selected EVS cutoff comes from the first stable crossing of the two composite regime lines.",
        width = 92
      ),
      x = "Percentile center",
      y = "Composite score"
    ) +
    manuscript_theme()
}

plot_regime_difference_curve <- function(combined_cutoff_info, comparison_name) {
  df <- combined_cutoff_info$combined_wave_map
  if (is.null(df) || !nrow(df) || !"regime_difference" %in% names(df)) {
    return(NULL)
  }

  sc <- combined_cutoff_info$selected_crossing
  crossing_pct  <- if (!is.null(sc) && nrow(sc)) sc$crossing_percentile[1] else NA_real_
  crossing_rank <- if (!is.null(sc) && nrow(sc)) sc$crossing_rank[1] else NA_integer_

  ggplot(df, aes(percentile, regime_difference)) +
    geom_hline(yintercept = 0, linewidth = crossing_plot_hline_width, colour = "grey50") +
    geom_line(linewidth = crossing_plot_line_width, colour = plot_palette$threshold) +
    geom_vline(
      xintercept = crossing_pct,
      linetype = "dashed",
      linewidth = crossing_plot_vline_width,
      colour = plot_palette$threshold
    ) +
    geom_point(
      data = data.frame(percentile = crossing_pct, regime_difference = 0),
      aes(x = percentile, y = regime_difference),
      inherit.aes = FALSE,
      size = crossing_plot_point_size,
      colour = plot_palette$threshold
    ) +
    annotate(
      "label",
      x = crossing_pct,
      y = 0,
      label = paste0("Rank = ", crossing_rank),
      fill = "white",
      colour = plot_palette$threshold,
      size = 2.9,
      label.size = 0.15,
      vjust = -1.0
    ) +
    labs(
      title = paste0(comparison_name, " | regime-difference curve"),
      subtitle = compact_caption(
        "Positive values indicate local IOD dominance and negative values indicate local CV² dominance. The first stable zero crossing defines the shared EVS cutoff.",
        width = 92
      ),
      x = "Percentile center",
      y = "IOD minus CV²"
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
      comparison_name    = comparison_name,
      selected_reason    = combined_cutoff_info$selected_reason,
      crossing_percentile = NA_real_,
      crossing_rank      = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    comparison_name         = comparison_name,
    selected_reason         = combined_cutoff_info$selected_reason,
    crossing_id             = sc$crossing_id[1],
    crossing_percentile     = sc$crossing_percentile[1],
    crossing_rank           = sc$crossing_rank[1],
    percentile_left         = sc$percentile_left[1],
    percentile_right        = sc$percentile_right[1],
    regime_difference_left  = sc$regime_difference_left[1],
    regime_difference_right = sc$regime_difference_right[1],
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# PART 3
# =============================================================================

get_condition_coef <- function(dds) {
  rn <- resultsNames(dds)
  hit <- grep("condition.*trt.*untrt", rn, value = TRUE)
  if (!length(hit)) {
    stop("Could not locate the DESeq2 treatment-versus-control coefficient.", call. = FALSE)
  }
  hit[1]
}

classify_effect_strength <- function(res_strong_padj, res_weak_padj, alpha = alpha_level) {
  dplyr::case_when(
    !is.na(res_strong_padj) & (res_strong_padj < alpha) ~ "strong_effect",
    !is.na(res_weak_padj)   & (res_weak_padj   < alpha) ~ "weak_effect",
    TRUE ~ "intermediate_or_nonsignificant"
  )
}

read_sequence_count_matrix <- function(file_path = count_file) {
  raw_df <- utils::read.csv(
    file_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!nrow(raw_df)) {
    stop("The raw count matrix is empty.", call. = FALSE)
  }

  first_col <- names(raw_df)[1]
  names(raw_df)[1] <- "feature_id"

  annotation_candidates <- c(
    "gene_symbol", "symbol", "gene", "Gene", "GENE",
    "feature_label", "annotation", "Annotation"
  )

  present_annotation_cols <- intersect(annotation_candidates, names(raw_df))
  sample_cols <- intersect(colnames(raw_df), rownames(meta_all))

  if (!length(sample_cols)) {
    stop(
      "No sample columns from the embedded metadata were found in the raw count file.",
      call. = FALSE
    )
  }

  annotation_df <- raw_df[, c("feature_id", present_annotation_cols), drop = FALSE]
  annotation_df$feature_id <- as.character(annotation_df$feature_id)
  annotation_df <- annotation_df[!duplicated(annotation_df$feature_id), , drop = FALSE]

  count_df <- raw_df[, c("feature_id", sample_cols), drop = FALSE]
  count_df$feature_id <- as.character(count_df$feature_id)

  for (cc in sample_cols) {
    count_df[[cc]] <- suppressWarnings(as.numeric(count_df[[cc]]))
    count_df[[cc]][!is.finite(count_df[[cc]]) | is.na(count_df[[cc]])] <- 0
  }

  count_mat <- as.matrix(count_df[, sample_cols, drop = FALSE])
  rownames(count_mat) <- count_df$feature_id
  storage.mode(count_mat) <- "numeric"

  keep_rows <- !is.na(rownames(count_mat)) & nzchar(rownames(count_mat))
  count_mat <- count_mat[keep_rows, , drop = FALSE]
  annotation_df <- annotation_df[match(rownames(count_mat), annotation_df$feature_id), , drop = FALSE]

  list(
    count_mat = count_mat,
    annotation_df = annotation_df,
    source_first_column = first_col,
    sample_cols = sample_cols
  )
}

get_comparison_sample_ids <- function(group1_prefix, group2_prefix, sample_ids) {
  group1_ids <- grep(paste0("^", group1_prefix, "_"), sample_ids, value = TRUE)
  group2_ids <- grep(paste0("^", group2_prefix, "_"), sample_ids, value = TRUE)

  if (!length(group1_ids) || !length(group2_ids)) {
    stop(
      sprintf(
        "Comparison sample extraction failed for %s versus %s.",
        group1_prefix, group2_prefix
      ),
      call. = FALSE
    )
  }

  list(
    treatment_ids = group1_ids,
    control_ids   = group2_ids,
    all_ids       = c(group1_ids, group2_ids)
  )
}

build_comparison_count_set <- function(full_count_mat,
                                       full_annotation_df,
                                       comparison_name,
                                       group1_prefix,
                                       group2_prefix) {
  sample_ids <- colnames(full_count_mat)
  cmp_ids <- get_comparison_sample_ids(group1_prefix, group2_prefix, sample_ids)

  count_df <- full_count_mat[, cmp_ids$all_ids, drop = FALSE]
  coldata  <- meta_all[cmp_ids$all_ids, , drop = FALSE]

  if (!all(colnames(count_df) == rownames(coldata))) {
    stop(
      sprintf("[%s] Count columns and coldata rows are not aligned.", comparison_name),
      call. = FALSE
    )
  }

  keep <- rowSums(count_df) > 0
  count_df <- count_df[keep, , drop = FALSE]

  annotation_df <- full_annotation_df[match(rownames(count_df), full_annotation_df$feature_id), , drop = FALSE]

  list(
    count_df = count_df,
    coldata = coldata,
    annotation_df = annotation_df,
    comparison_name = comparison_name,
    treatment_ids = cmp_ids$treatment_ids,
    control_ids = cmp_ids$control_ids
  )
}

compute_pc1_loading_table <- function(count_df,
                                      condition_vector,
                                      preprocess_mode = c("normalized", "raw_counts"),
                                      feature_annotation = NULL,
                                      comparison_name = NULL,
                                      group_label = NULL) {
  preprocess_mode <- match.arg(preprocess_mode)

  x <- as.matrix(count_df)
  storage.mode(x) <- "numeric"

  if (preprocess_mode == "normalized") {
    tmp_coldata <- data.frame(
      condition = factor(as.character(condition_vector), levels = c("untrt", "trt")),
      row.names = colnames(x),
      stringsAsFactors = FALSE
    )

    dds_tmp <- DESeqDataSetFromMatrix(
      countData = round(x),
      colData = tmp_coldata,
      design = ~ 1
    )
    dds_tmp <- dds_tmp[rowSums(counts(dds_tmp)) > 0, ]
    dds_tmp <- estimateSizeFactors(dds_tmp)
    x <- counts(dds_tmp, normalized = TRUE)
  }

  x_log <- log2(x + 1)

  pca_fit <- prcomp(x_log, center = TRUE, scale. = FALSE, rank. = 1)
  pc1_load <- as.numeric(pca_fit$rotation[, 1])
  names(pc1_load) <- rownames(pca_fit$rotation)

  loading_tbl <- data.frame(
    feature_id       = names(pc1_load),
    pc1_loading      = pc1_load,
    pc1_loading_abs  = abs(pc1_load),
    stringsAsFactors = FALSE
  )

  loading_tbl <- loading_tbl[order(-loading_tbl$pc1_loading_abs, loading_tbl$feature_id), , drop = FALSE]
  loading_tbl$rank <- seq_len(nrow(loading_tbl))

  if (!is.null(feature_annotation) && "feature_id" %in% names(feature_annotation)) {
    loading_tbl <- dplyr::left_join(
      loading_tbl,
      feature_annotation,
      by = "feature_id"
    )
  }

  if (!is.null(comparison_name)) {
    loading_tbl$comparison_name <- comparison_name
  }
  if (!is.null(group_label)) {
    loading_tbl$group_label <- group_label
  }

  loading_tbl
}

build_comparison_loading_pair <- function(comparison_obj,
                                          preprocess_mode = c("normalized", "raw_counts")) {
  preprocess_mode <- match.arg(preprocess_mode)

  count_df  <- comparison_obj$count_df
  coldata   <- comparison_obj$coldata
  annot_df  <- comparison_obj$annotation_df
  cmp_name  <- comparison_obj$comparison_name

  trt_ids   <- rownames(coldata)[coldata$condition == "trt"]
  ctrl_ids  <- rownames(coldata)[coldata$condition == "untrt"]

  trt_tbl <- compute_pc1_loading_table(
    count_df = count_df[, trt_ids, drop = FALSE],
    condition_vector = coldata[trt_ids, "condition"],
    preprocess_mode = preprocess_mode,
    feature_annotation = annot_df,
    comparison_name = cmp_name,
    group_label = "treatment"
  )

  ctrl_tbl <- compute_pc1_loading_table(
    count_df = count_df[, ctrl_ids, drop = FALSE],
    condition_vector = coldata[ctrl_ids, "condition"],
    preprocess_mode = preprocess_mode,
    feature_annotation = annot_df,
    comparison_name = cmp_name,
    group_label = "control"
  )

  list(
    treatment_loading = trt_tbl,
    control_loading   = ctrl_tbl
  )
}

merge_loading_pair_to_combined_rank <- function(trt_tbl, ctrl_tbl) {
  keep_cols <- intersect(
    c("feature_id", "pc1_loading", "pc1_loading_abs", "rank"),
    names(trt_tbl)
  )

  trt_use <- trt_tbl[, keep_cols, drop = FALSE]
  ctrl_use <- ctrl_tbl[, keep_cols, drop = FALSE]

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

  merged
}

compute_thresholded_effect_tests <- function(res_df, lfc_threshold = lfc_boundary) {
  assert_required_columns(
    res_df,
    c("feature_id", "lfc_shrunk", "lfcSE_shrunk"),
    "res_df for thresholded effect tests"
  )

  beta <- as.numeric(res_df$lfc_shrunk)
  se   <- as.numeric(res_df$lfcSE_shrunk)

  valid <- is.finite(beta) & !is.na(beta) & is.finite(se) & !is.na(se) & se > 0

  z_greater <- rep(NA_real_, length(beta))
  p_greater <- rep(NA_real_, length(beta))

  z_less <- rep(NA_real_, length(beta))
  p_less <- rep(NA_real_, length(beta))

  z_greater[valid] <- (abs(beta[valid]) - lfc_threshold) / se[valid]
  p_greater[valid] <- stats::pnorm(z_greater[valid], lower.tail = FALSE)

  z_less[valid] <- (lfc_threshold - abs(beta[valid])) / se[valid]
  p_less[valid] <- stats::pnorm(z_less[valid], lower.tail = FALSE)

  data.frame(
    feature_id = res_df$feature_id,
    resGA_z = z_greater,
    resGA_pvalue = clip_probabilities(p_greater),
    resGA_padj = clip_probabilities(stats::p.adjust(p_greater, method = "BH")),
    resLA_z = z_less,
    resLA_pvalue = clip_probabilities(p_less),
    resLA_padj = clip_probabilities(stats::p.adjust(p_less, method = "BH")),
    stringsAsFactors = FALSE
  )
}

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
  design_formula <- make_design_formula()

  count_mat <- round(as.matrix(count_mat))
  storage.mode(count_mat) <- "numeric"

  if (nrow(count_mat) == 0 || ncol(count_mat) == 0) {
    stop(sprintf("[%s] count_mat is empty before DESeq2.", dataset_name), call. = FALSE)
  }

  keep_nonzero <- rowSums(count_mat, na.rm = TRUE) > 0
  count_mat <- count_mat[keep_nonzero, , drop = FALSE]

  if (nrow(count_mat) == 0) {
    stop(sprintf("[%s] all rows were zero after split and zero-row filtering.", dataset_name), call. = FALSE)
  }

  if (!is.null(annot_df) && "feature_id" %in% names(annot_df)) {
    annot_df <- annot_df[match(rownames(count_mat), annot_df$feature_id), , drop = FALSE]
  }

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
  shr       <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)
  shr_df    <- as.data.frame(shr)
  shr_df$feature_id <- as.character(rownames(shr_df))

  res_df <- dplyr::left_join(
    res_df,
    shr_df[, c("feature_id", "log2FoldChange", "lfcSE"), drop = FALSE],
    by     = "feature_id",
    suffix = c("", "_shrunk")
  )

  if ("log2FoldChange_shrunk" %in% colnames(res_df)) {
    colnames(res_df)[colnames(res_df) == "log2FoldChange_shrunk"] <- "lfc_shrunk"
  } else {
    res_df$lfc_shrunk <- res_df$log2FoldChange
  }

  if ("lfcSE_shrunk" %in% colnames(res_df)) {
    colnames(res_df)[colnames(res_df) == "lfcSE_shrunk"] <- "lfcSE_shrunk"
  } else if ("lfcSE" %in% colnames(res_df)) {
    res_df$lfcSE_shrunk <- res_df$lfcSE
  }

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
  } else {
    hbfss_threshold_dataset  <- abs(log10(hc_p_threshold_dataset)) * lfc_boundary
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
  if ("gene_symbol" %in% names(annot_df)) {
    annot_df$gene_symbol <- as.character(annot_df$gene_symbol)
  } else {
    annot_df$gene_symbol <- NA_character_
  }

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

split_by_evs_rank <- function(count_df,
                              annotation_df,
                              merged_loading_tbl,
                              cutoff_rank,
                              comparison_name) {
  assert_required_columns(
    merged_loading_tbl,
    c("feature_id", "combined_rank"),
    "merged_loading_tbl for split_by_evs_rank"
  )

  cutoff_rank <- max(1L, min(as.integer(cutoff_rank), nrow(merged_loading_tbl)))

  leading_ids <- merged_loading_tbl$feature_id[merged_loading_tbl$combined_rank <= cutoff_rank]
  remainder_ids <- merged_loading_tbl$feature_id[merged_loading_tbl$combined_rank > cutoff_rank]

  leading_counts <- count_df[rownames(count_df) %in% leading_ids, , drop = FALSE]
  remainder_counts <- count_df[rownames(count_df) %in% remainder_ids, , drop = FALSE]

  leading_annot <- annotation_df[match(rownames(leading_counts), annotation_df$feature_id), , drop = FALSE]
  remainder_annot <- annotation_df[match(rownames(remainder_counts), annotation_df$feature_id), , drop = FALSE]

  list(
    original = list(
      count_df = count_df,
      annotation_df = annotation_df,
      dataset_name = "original",
      comparison_name = comparison_name
    ),
    leading_edge = list(
      count_df = leading_counts,
      annotation_df = leading_annot,
      dataset_name = "leading_edge",
      comparison_name = comparison_name
    ),
    remainder = list(
      count_df = remainder_counts,
      annotation_df = remainder_annot,
      dataset_name = "remainder",
      comparison_name = comparison_name
    )
  )
}

append_loading_metadata_to_results <- function(res_df, merged_loading_tbl) {
  keep_cols <- intersect(
    c(
      "feature_id",
      "pc1_loading_trt", "pc1_loading_abs_trt", "rank_trt",
      "pc1_loading_ctrl", "pc1_loading_abs_ctrl", "rank_ctrl",
      "combined_loading", "combined_rank"
    ),
    names(merged_loading_tbl)
  )

  dplyr::left_join(
    res_df,
    merged_loading_tbl[, keep_cols, drop = FALSE],
    by = "feature_id"
  )
}

run_single_dataset_analysis <- function(dataset_obj,
                                        coldata,
                                        merged_loading_tbl,
                                        comparison_name,
                                        dataset_name) {
  core <- run_core_analysis(
    count_mat   = round(as.matrix(dataset_obj$count_df)),
    coldata     = coldata,
    dataset_name = paste(comparison_name, dataset_name, sep = " | "),
    annot_df    = dataset_obj$annotation_df
  )

  res_df <- core$results
  res_df <- append_loading_metadata_to_results(res_df, merged_loading_tbl)

  res_df
}

build_dataset_summary_table <- function(res_df, comparison_name, dataset_name) {
  hc_vals <- stats::na.omit(res_df$hc_p_threshold_dataset)
  hbfss_vals <- stats::na.omit(res_df$hbfss_threshold_dataset)

  data.frame(
    comparison_name        = comparison_name,
    dataset_name           = dataset_name,
    n_features             = nrow(res_df),
    n_standard_significant = sum(res_df$standard_significant, na.rm = TRUE),
    n_deseq2_strong        = sum(res_df$deseq2_strong_call, na.rm = TRUE),
    n_deseq2_weak          = sum(res_df$deseq2_weak_call, na.rm = TRUE),
    n_hbfss                = sum(res_df$HBFSS_significant, na.rm = TRUE),
    n_overlap_strong_hbfss = sum(res_df$overlap_call, na.rm = TRUE),
    hc_p_threshold_dataset = if (length(hc_vals)) hc_vals[1] else NA_real_,
    hbfss_threshold_dataset = if (length(hbfss_vals)) hbfss_vals[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}

build_minimal_results_export <- function(res_df) {
  keep_cols <- intersect(
    c(
      "feature_id",
      "gene_symbol",
      "dataset_name",
      "baseMean",
      "lfc_shrunk",
      "regulation_direction",
      "lfdr",
      "pvalue",
      "padj",
      "empirical_p",
      "empirical_q",
      "HBFSS",
      "hc_p_threshold_dataset",
      "hbfss_threshold_dataset",
      "resLA_padj",
      "resGA_padj",
      "standard_significant",
      "deseq2_strong_call",
      "deseq2_weak_call",
      "HBFSS_significant",
      "effect_class",
      "combined_rank",
      "combined_loading"
    ),
    names(res_df)
  )

  out <- res_df[, keep_cols, drop = FALSE]
  out <- out[order(out$combined_rank, out$feature_id), , drop = FALSE]
  out
}

# =============================================================================
# PART 4
# =============================================================================

volcano_color_values <- c(
  "background" = plot_palette$background,
  "strong"     = plot_palette$strong,
  "weak"       = plot_palette$weak,
  "hbfss"      = plot_palette$hbfss,
  "overlap"    = plot_palette$overlap
)

volcano_color_labels <- c(
  "background" = "Background or nonsignificant",
  "strong"     = "DESeq2 strong effect",
  "weak"       = "DESeq2 weak effect",
  "hbfss"      = "HBFSS only",
  "overlap"    = "DESeq2 standard and HBFSS"
)

volcano_shape_values <- c(
  "background"  = 16,
  "deseq2_only" = 17,
  "weak_effect" = 15,
  "hbfss_only"  = 18,
  "overlap"     = 23
)

volcano_shape_labels <- c(
  "background"  = "Background",
  "deseq2_only" = "DESeq2 standard only",
  "weak_effect" = "DESeq2 weak effect",
  "hbfss_only"  = "HBFSS only",
  "overlap"     = "DESeq2 standard and HBFSS"
)

derive_plot_classes <- function(res_df) {
  df <- res_df

  if (!"standard_significant" %in% names(df)) {
    df$standard_significant <- FALSE
  }
  if (!"HBFSS_significant" %in% names(df)) {
    df$HBFSS_significant <- FALSE
  }
  if (!"deseq2_strong_call" %in% names(df)) {
    df$deseq2_strong_call <- FALSE
  }
  if (!"deseq2_weak_call" %in% names(df)) {
    df$deseq2_weak_call <- FALSE
  }
  if (!"overlap_call" %in% names(df)) {
    df$overlap_call <- df$standard_significant & df$HBFSS_significant
  }
  if (!"HBFSS_only_call" %in% names(df)) {
    df$HBFSS_only_call <- df$HBFSS_significant & !df$standard_significant
  }

  df$volcano_display_class <- dplyr::case_when(
    df$overlap_call ~ "strong_and_hbfss",
    df$deseq2_strong_call & !df$HBFSS_significant ~ "strong_deseq2_only",
    df$deseq2_weak_call ~ "weak_deseq2_effect",
    df$HBFSS_only_call ~ "hbfss_only",
    TRUE ~ "background"
  )

  df$volcano_color <- dplyr::case_when(
    df$volcano_display_class == "strong_and_hbfss"   ~ "overlap",
    df$volcano_display_class == "strong_deseq2_only" ~ "strong",
    df$volcano_display_class == "weak_deseq2_effect" ~ "weak",
    df$volcano_display_class == "hbfss_only"         ~ "hbfss",
    TRUE                                             ~ "background"
  )

  df$volcano_shape <- dplyr::case_when(
    df$volcano_display_class == "strong_and_hbfss"   ~ "overlap",
    df$volcano_display_class == "strong_deseq2_only" ~ "deseq2_only",
    df$volcano_display_class == "weak_deseq2_effect" ~ "weak_effect",
    df$volcano_display_class == "hbfss_only"         ~ "hbfss_only",
    TRUE                                             ~ "background"
  )

  df$label_candidate <- df$volcano_display_class != "background"

  df
}

pick_volcano_labels <- function(res_df,
                                max_labels = 18L,
                                ranking_mode = c("HBFSS_then_lfc", "lfc_then_p")) {
  ranking_mode <- match.arg(ranking_mode)

  df <- derive_plot_classes(res_df)
  lab_df <- df[df$label_candidate, , drop = FALSE]
  if (!nrow(lab_df)) return(lab_df[0, , drop = FALSE])

  if (!"gene_symbol" %in% names(lab_df)) {
    lab_df$gene_symbol <- NA_character_
  }

  usable_label <- !is.na(lab_df$gene_symbol) & nzchar(lab_df$gene_symbol) &
    !grepl("^LOC|^NA$|^\\.$|^-$", lab_df$gene_symbol)

  lab_df$plot_label <- ifelse(usable_label, lab_df$gene_symbol, lab_df$feature_id)

  signal_vec <- if ("HBFSS" %in% names(lab_df)) {
    ifelse(is.finite(lab_df$HBFSS), lab_df$HBFSS, -Inf)
  } else {
    ifelse(is.finite(lab_df$neglog10_empirical_p), lab_df$neglog10_empirical_p, -Inf)
  }

  if (ranking_mode == "HBFSS_then_lfc") {
    ord <- order(
      -signal_vec,
      -abs(ifelse(is.finite(lab_df$lfc_shrunk), lab_df$lfc_shrunk, 0)),
      ifelse(is.finite(lab_df$padj), lab_df$padj, Inf)
    )
  } else {
    ord <- order(
      -abs(ifelse(is.finite(lab_df$lfc_shrunk), lab_df$lfc_shrunk, 0)),
      ifelse(is.finite(lab_df$padj), lab_df$padj, Inf)
    )
  }

  lab_df <- lab_df[ord, , drop = FALSE]
  lab_df <- lab_df[!duplicated(lab_df$plot_label), , drop = FALSE]
  utils::head(lab_df, max_labels)
}

prepare_volcano_df <- function(res_df) {
  df <- derive_plot_classes(res_df)

  assert_required_columns(
    df,
    c("feature_id", "lfc_shrunk", "empirical_p", "padj", "volcano_color", "volcano_shape"),
    "res_df for volcano plotting"
  )

  df$x_value <- as.numeric(df$lfc_shrunk)
  df$y_value <- safe_neglog10(ifelse(is.finite(df$empirical_p), df$empirical_p, df$padj))
  df$y_value[!is.finite(df$y_value) | is.na(df$y_value)] <- 0

  df$volcano_color <- factor(
    df$volcano_color,
    levels = c("background", "weak", "strong", "hbfss", "overlap")
  )

  df$volcano_shape <- factor(
    df$volcano_shape,
    levels = c("background", "weak_effect", "deseq2_only", "hbfss_only", "overlap")
  )

  df
}

plot_volcano_manuscript <- function(res_df,
                                    comparison_name,
                                    dataset_name,
                                    max_labels = 18L) {
  df <- prepare_volcano_df(res_df)
  lab_df <- pick_volcano_labels(df, max_labels = max_labels)

  thr_vals <- stats::na.omit(df$hbfss_threshold_dataset)
  hbfss_thr <- if (length(thr_vals)) thr_vals[1] else NA_real_

  x_max <- max(abs(df$x_value), na.rm = TRUE)
  x_max <- max(1.25, x_max * 1.08)

  p <- ggplot(df, aes(x = x_value, y = y_value)) +
    geom_vline(
      xintercept = c(-lfc_boundary, lfc_boundary),
      linetype = "dashed",
      linewidth = LINE_WIDTH_THRESH,
      colour = plot_palette$threshold
    ) +
    geom_hline(
      yintercept = 0,
      linewidth = LINE_WIDTH_ZERO,
      colour = "grey80"
    ) +
    geom_point(
      aes(color = volcano_color, shape = volcano_shape),
      size = POINT_SIZE_PRIMARY,
      alpha = POINT_ALPHA_PRIMARY,
      stroke = POINT_STROKE
    ) +
    scale_color_manual(
      values = volcano_color_values,
      labels = volcano_color_labels,
      breaks = c("weak", "strong", "hbfss", "overlap", "background"),
      name = "Color class"
    ) +
    scale_shape_manual(
      values = volcano_shape_values,
      labels = volcano_shape_labels,
      breaks = c("weak_effect", "deseq2_only", "hbfss_only", "overlap", "background"),
      name = "Method support"
    ) +
    labs(
      title = paste0(comparison_name, " | ", pretty_dataset_label(dataset_name)),
      subtitle = compact_caption(
        paste0(
          "Blue marks DESeq2 weak effects below the absolute log2 fold-change boundary of ",
          lfc_boundary,
          ". Red marks DESeq2 strong effects. Orange marks HBFSS only. Purple marks DESeq2 standard significance with HBFSS overlap."
        ),
        width = 105
      ),
      x = "Shrunken log2 fold change",
      y = expression(-log[10]("empirical p value"))
    ) +
    coord_cartesian(xlim = c(-x_max, x_max), clip = "off") +
    manuscript_theme() +
    theme(
      plot.margin = margin(t = 12, r = 20, b = 14, l = 14)
    )

  if (is.finite(hbfss_thr) && !is.na(hbfss_thr) && hbfss_thr > 0) {
    p <- p +
      geom_hline(
        yintercept = abs(log10(hbfss_thr)),
        linetype = "dotted",
        linewidth = LINE_WIDTH_BOUNDARY,
        colour = plot_palette$hbfss
      ) +
      annotate(
        "label",
        x = x_max * 0.78,
        y = abs(log10(hbfss_thr)),
        label = paste0("HC p = ", signif(hbfss_thr, 4)),
        colour = plot_palette$hbfss,
        fill = "white",
        size = 2.8,
        label.size = 0.15,
        vjust = -0.6
      )
  }

  if (nrow(lab_df)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = lab_df,
        aes(label = plot_label),
        size = 2.15,
        max.overlaps = Inf,
        box.padding = 0.22,
        point.padding = 0.12,
        segment.size = 0.22,
        min.segment.length = 0,
        show.legend = FALSE
      )
  }

  p
}

plot_dispersion_manuscript <- function(res_df,
                                       comparison_name,
                                       dataset_name) {
  df <- derive_plot_classes(res_df)

  assert_required_columns(
    df,
    c("baseMean", "dispGeneEst", "dispFit", "volcano_color", "volcano_shape"),
    "res_df for dispersion plotting"
  )

  df$log_baseMean   <- safe_log10(df$baseMean + 1)
  df$log_dispGeneEst <- safe_log10(df$dispGeneEst)
  df$log_dispFit     <- safe_log10(df$dispFit)

  fit_line_df <- df[
    is.finite(df$log_baseMean) & is.finite(df$log_dispFit),
    c("log_baseMean", "log_dispFit"),
    drop = FALSE
  ]
  fit_line_df <- fit_line_df[order(fit_line_df$log_baseMean), , drop = FALSE]

  ggplot(df, aes(log_baseMean, log_dispGeneEst)) +
    geom_point(
      aes(color = volcano_color, shape = volcano_shape),
      size = POINT_SIZE_DISP,
      alpha = POINT_ALPHA_DISP,
      stroke = POINT_STROKE
    ) +
    geom_line(
      data = fit_line_df,
      aes(x = log_baseMean, y = log_dispFit),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.65
    ) +
    scale_color_manual(
      values = volcano_color_values,
      labels = volcano_color_labels,
      breaks = c("weak", "strong", "hbfss", "overlap", "background"),
      name = "Color class"
    ) +
    scale_shape_manual(
      values = volcano_shape_values,
      labels = volcano_shape_labels,
      breaks = c("weak_effect", "deseq2_only", "hbfss_only", "overlap", "background"),
      name = "Method support"
    ) +
    labs(
      title = paste0(comparison_name, " | ", pretty_dataset_label(dataset_name), " dispersion"),
      subtitle = compact_caption(
        "Points are gene-wise dispersion estimates colored by final volcano display class. The black line is the fitted DESeq2 dispersion trend.",
        width = 92
      ),
      x = "log10(baseMean + 1)",
      y = "log10(gene-wise dispersion)"
    ) +
    manuscript_theme()
}

plot_pca_panel_pair <- function(dataset_obj,
                                coldata,
                                comparison_name,
                                dataset_name) {
  p_norm <- safe_plot_build(
    compute_dataset_pca_plot(
      count_df = dataset_obj$count_df,
      coldata = coldata,
      dataset_name = dataset_name,
      preprocessing = "normalized"
    ),
    label = paste(comparison_name, dataset_name, "PCA normalized")
  )

  p_raw <- safe_plot_build(
    compute_dataset_pca_plot(
      count_df = dataset_obj$count_df,
      coldata = coldata,
      dataset_name = dataset_name,
      preprocessing = "raw_counts"
    ),
    label = paste(comparison_name, dataset_name, "PCA raw")
  )

  if (is.null(p_norm) || is.null(p_raw)) return(NULL)

  arrangeGrob(
    p_norm, p_raw,
    ncol = 2,
    top = textGrob(
      paste0(comparison_name, " | ", pretty_dataset_label(dataset_name), " PCA"),
      gp = gpar(fontface = "bold", cex = 1.05)
    )
  )
}

plot_empirical_p_histogram <- function(res_df,
                                       comparison_name,
                                       dataset_name) {
  if (!export_optional_empirical_p_histograms) return(NULL)
  if (!"empirical_p" %in% names(res_df)) return(NULL)

  df <- data.frame(empirical_p = res_df$empirical_p)
  df <- df[is.finite(df$empirical_p) & !is.na(df$empirical_p), , drop = FALSE]
  if (!nrow(df)) return(NULL)

  hc_vals <- stats::na.omit(res_df$hc_p_threshold_dataset)
  hc_p_threshold <- if (length(hc_vals)) hc_vals[1] else NA_real_

  p <- ggplot(df, aes(empirical_p)) +
    geom_histogram(
      bins = HIST_BINS,
      fill = plot_palette$histogram,
      color = HIST_COLOR
    ) +
    labs(
      title = paste0(comparison_name, " | ", pretty_dataset_label(dataset_name), " empirical p"),
      subtitle = compact_caption(
        "Histogram of empirical-null p values derived from the DESeq2 Wald statistic distribution.",
        width = 92
      ),
      x = "Empirical p value",
      y = "Count"
    ) +
    manuscript_theme()

  if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold)) {
    p <- p +
      geom_vline(
        xintercept = hc_p_threshold,
        linetype = "dashed",
        linewidth = 0.9,
        colour = plot_palette$threshold
      ) +
      annotate(
        "label",
        x = hc_p_threshold,
        y = Inf,
        label = paste0("HC p = ", signif(hc_p_threshold, 4)),
        vjust = 1.4,
        fill = "white",
        colour = plot_palette$threshold,
        size = 2.8,
        label.size = 0.15
      )
  }

  p
}

plot_hbfss_histogram <- function(res_df,
                                 comparison_name,
                                 dataset_name) {
  if (!export_optional_hbfss_distributions) return(NULL)
  if (!"HBFSS" %in% names(res_df)) return(NULL)

  df <- data.frame(HBFSS = res_df$HBFSS)
  df <- df[is.finite(df$HBFSS) & !is.na(df$HBFSS), , drop = FALSE]
  if (!nrow(df)) return(NULL)

  thr_vals <- stats::na.omit(res_df$hbfss_threshold_dataset)
  thr <- if (length(thr_vals)) thr_vals[1] else NA_real_

  p <- ggplot(df, aes(HBFSS)) +
    geom_histogram(
      bins = HIST_BINS,
      fill = plot_palette$histogram,
      color = HIST_COLOR
    ) +
    labs(
      title = paste0(comparison_name, " | ", pretty_dataset_label(dataset_name), " HBFSS"),
      subtitle = compact_caption(
        "HBFSS equals the absolute value of shrunken log2 fold change multiplied by log10 empirical p signal.",
        width = 92
      ),
      x = "HBFSS",
      y = "Count"
    ) +
    manuscript_theme()

  if (is.finite(thr) && !is.na(thr)) {
    p <- p +
      geom_vline(
        xintercept = thr,
        linetype = "dashed",
        linewidth = 0.9,
        colour = plot_palette$threshold
      ) +
      annotate(
        "label",
        x = thr,
        y = Inf,
        label = paste0("HBFSS threshold = ", signif(thr, 4)),
        vjust = 1.4,
        fill = "white",
        colour = plot_palette$threshold,
        size = 2.8,
        label.size = 0.15
      )
  }

  p
}

build_manuscript_gene_table <- function(res_df,
                                        comparison_name,
                                        dataset_name) {
  df <- derive_plot_classes(res_df)

  keep_cols <- intersect(
    c(
      "feature_id",
      "gene_symbol",
      "dataset_name",
      "baseMean",
      "lfc_shrunk",
      "regulation_direction",
      "lfdr",
      "pvalue",
      "padj",
      "empirical_p",
      "empirical_q",
      "HBFSS",
      "hc_p_threshold_dataset",
      "hbfss_threshold_dataset",
      "resLA_padj",
      "resGA_padj",
      "standard_significant",
      "deseq2_strong_call",
      "deseq2_weak_call",
      "HBFSS_significant",
      "effect_class",
      "volcano_display_class",
      "combined_rank",
      "combined_loading"
    ),
    names(df)
  )

  out <- df[, keep_cols, drop = FALSE]

  out <- out[order(
    factor(out$volcano_display_class,
           levels = c(
             "strong_and_hbfss",
             "strong_deseq2_only",
             "weak_deseq2_effect",
             "hbfss_only",
             "background"
           )),
    out$padj,
    -abs(out$lfc_shrunk)
  ), , drop = FALSE]

  out
}

build_manuscript_class_counts <- function(res_df,
                                          comparison_name,
                                          dataset_name) {
  df <- derive_plot_classes(res_df)

  tab <- table(
    factor(
      df$volcano_display_class,
      levels = c(
        "strong_and_hbfss",
        "strong_deseq2_only",
        "weak_deseq2_effect",
        "hbfss_only",
        "background"
      )
    )
  )

  data.frame(
    comparison_name = comparison_name,
    dataset_name = dataset_name,
    display_class = names(tab),
    n_features = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

build_cutoff_manifest_row <- function(comparison_name,
                                      cutoff_info) {
  data.frame(
    comparison_name = comparison_name,
    cutoff_method = cutoff_info$method,
    selected_reason = cutoff_info$selected_reason,
    top_n_actual = cutoff_info$top_n_actual,
    cutoff_quantile = cutoff_info$cutoff_quantile,
    cutoff_value = cutoff_info$cutoff_value,
    stringsAsFactors = FALSE
  )
}

export_single_dataset_tables <- function(res_df,
                                         comparison_name,
                                         dataset_name,
                                         table_dir_local = table_dir) {
  gene_tbl  <- build_manuscript_gene_table(res_df, comparison_name, dataset_name)
  class_tbl <- build_manuscript_class_counts(res_df, comparison_name, dataset_name)
  sum_tbl   <- build_dataset_summary_table(res_df, comparison_name, dataset_name)

  base_stub <- paste(comparison_name, dataset_name, sep = "__")

  save_csv(
    gene_tbl,
    file.path(table_dir_local, paste0(base_stub, "__manuscript_gene_table.csv"))
  )
  save_csv(
    class_tbl,
    file.path(table_dir_local, paste0(base_stub, "__display_class_counts.csv"))
  )
  save_csv(
    sum_tbl,
    file.path(table_dir_local, paste0(base_stub, "__dataset_summary.csv"))
  )

  invisible(
    list(
      gene_table = gene_tbl,
      class_table = class_tbl,
      summary_table = sum_tbl
    )
  )
}

export_single_dataset_figures <- function(res_df,
                                          dataset_obj,
                                          coldata,
                                          comparison_name,
                                          dataset_name,
                                          figure_dir_local = figure_dir) {
  base_stub <- paste(comparison_name, dataset_name, sep = "__")

  volc <- safe_plot_build(
    plot_volcano_manuscript(res_df, comparison_name, dataset_name),
    label = paste(base_stub, "volcano")
  )

  disp <- safe_plot_build(
    plot_dispersion_manuscript(res_df, comparison_name, dataset_name),
    label = paste(base_stub, "dispersion")
  )

  pca_panel <- safe_plot_build(
    plot_pca_panel_pair(dataset_obj, coldata, comparison_name, dataset_name),
    label = paste(base_stub, "pca")
  )

  emp_hist <- safe_plot_build(
    plot_empirical_p_histogram(res_df, comparison_name, dataset_name),
    label = paste(base_stub, "empirical p histogram")
  )

  hbfss_hist <- safe_plot_build(
    plot_hbfss_histogram(res_df, comparison_name, dataset_name),
    label = paste(base_stub, "HBFSS histogram")
  )

  if (!is.null(volc)) {
    save_grob(
      volc,
      file.path(figure_dir_local, paste0(base_stub, "__volcano.png")),
      width = 8.8,
      height = 6.5
    )
  }

  if (!is.null(disp)) {
    save_grob(
      disp,
      file.path(figure_dir_local, paste0(base_stub, "__dispersion.png")),
      width = 8.6,
      height = 6.2
    )
  }

  if (!is.null(pca_panel)) {
    save_grob(
      pca_panel,
      file.path(figure_dir_local, paste0(base_stub, "__pca.png")),
      width = 10.8,
      height = 5.2
    )
  }

  if (!is.null(emp_hist)) {
    save_grob(
      emp_hist,
      file.path(figure_dir_local, paste0(base_stub, "__empirical_p_histogram.png")),
      width = 7.2,
      height = 5.2
    )
  }

  if (!is.null(hbfss_hist)) {
    save_grob(
      hbfss_hist,
      file.path(figure_dir_local, paste0(base_stub, "__HBFSS_histogram.png")),
      width = 7.2,
      height = 5.2
    )
  }

  invisible(
    list(
      volcano = volc,
      dispersion = disp,
      pca = pca_panel,
      empirical_p_hist = emp_hist,
      HBFSS_hist = hbfss_hist
    )
  )
}

export_comparison_cutoff_figures <- function(combined_cutoff_info,
                                             comparison_name,
                                             figure_dir_local = figure_dir,
                                             table_dir_local = table_dir) {
  crossing_panel <- safe_plot_build(
    plot_crossing_summary_panel(combined_cutoff_info, comparison_name),
    label = paste(comparison_name, "crossing panel")
  )

  score_plot <- safe_plot_build(
    plot_combined_fourier_score(combined_cutoff_info, comparison_name),
    label = paste(comparison_name, "score plot")
  )

  crossing_tbl <- build_crossing_summary_table(combined_cutoff_info, comparison_name)
  cutoff_tbl   <- build_cutoff_manifest_row(comparison_name, combined_cutoff_info)

  if (!is.null(crossing_panel)) {
    save_grob(
      crossing_panel,
      file.path(figure_dir_local, paste0(comparison_name, "__regime_crossing_summary.png")),
      width = 9.0,
      height = 9.0
    )
  }

  if (!is.null(score_plot)) {
    save_grob(
      score_plot,
      file.path(figure_dir_local, paste0(comparison_name, "__descriptive_score_profile.png")),
      width = 8.4,
      height = 5.5
    )
  }

  save_csv(
    crossing_tbl,
    file.path(table_dir_local, paste0(comparison_name, "__crossing_summary.csv"))
  )

  save_csv(
    cutoff_tbl,
    file.path(table_dir_local, paste0(comparison_name, "__cutoff_manifest.csv"))
  )

  invisible(
    list(
      crossing_table = crossing_tbl,
      cutoff_table = cutoff_tbl
    )
  )
}

bind_dataset_summary_rows <- function(summary_rows) {
  if (!length(summary_rows)) return(data.frame())
  dplyr::bind_rows(summary_rows)
}

bind_class_count_rows <- function(class_rows) {
  if (!length(class_rows)) return(data.frame())
  dplyr::bind_rows(class_rows)
}

bind_cutoff_manifest_rows <- function(cutoff_rows) {
  if (!length(cutoff_rows)) return(data.frame())
  dplyr::bind_rows(cutoff_rows)
}

# =============================================================================
# PART 5
# =============================================================================

run_single_comparison_pipeline <- function(comparison_row,
                                           full_count_mat,
                                           full_annotation_df) {
  cmp_name <- comparison_row$comparison_name[[1]]
  grp1     <- comparison_row$group1_prefix[[1]]
  grp2     <- comparison_row$group2_prefix[[1]]

  message("")
  message("============================================================")
  message("Running comparison: ", cmp_name, " (", grp1, " vs ", grp2, ")")
  message("============================================================")

  comparison_obj <- build_comparison_count_set(
    full_count_mat       = full_count_mat,
    full_annotation_df   = full_annotation_df,
    comparison_name      = cmp_name,
    group1_prefix        = grp1,
    group2_prefix        = grp2
  )

  loading_pair <- build_comparison_loading_pair(
    comparison_obj = comparison_obj,
    preprocess_mode = "normalized"
  )

  feature_metric_tbl <- estimate_feature_metrics_for_loading(
    count_df = comparison_obj$count_df,
    coldata  = comparison_obj$coldata
  )

  trt_loading_tbl <- dplyr::left_join(
    loading_pair$treatment_loading,
    feature_metric_tbl,
    by = "feature_id"
  )

  ctrl_loading_tbl <- dplyr::left_join(
    loading_pair$control_loading,
    feature_metric_tbl,
    by = "feature_id"
  )

  merged_loading_tbl <- merge_loading_pair_to_combined_rank(
    trt_tbl  = trt_loading_tbl,
    ctrl_tbl = ctrl_loading_tbl
  )

  cutoff_info <- resolve_combined_fourier_cutoff(
    fit_trt_loading_tbl  = trt_loading_tbl,
    fit_ctrl_loading_tbl = ctrl_loading_tbl,
    fixed_top_n          = evs_fixed_top_n
  )

  if (isTRUE(evs_use_manual_rank_first) &&
      is.finite(evs_fixed_rank_override) &&
      !is.na(evs_fixed_rank_override)) {
    manual_rank <- max(1L, min(as.integer(evs_fixed_rank_override), nrow(merged_loading_tbl)))
    cutoff_info$top_n_actual <- manual_rank
    cutoff_info$method <- "manual_fixed_rank"
    cutoff_info$selected_reason <- "manual_fixed_rank"
  }

  cutoff_info$top_n_actual <- max(1L, min(as.integer(cutoff_info$top_n_actual), nrow(merged_loading_tbl)))
  cutoff_info$cutoff_value <- merged_loading_tbl$combined_loading[cutoff_info$top_n_actual]
  cutoff_info$cutoff_quantile <- 1 - (cutoff_info$top_n_actual / nrow(merged_loading_tbl))

  cutoff_rank <- cutoff_info$top_n_actual

  dataset_split <- split_by_evs_rank(
    count_df            = comparison_obj$count_df,
    annotation_df       = comparison_obj$annotation_df,
    merged_loading_tbl  = merged_loading_tbl,
    cutoff_rank         = cutoff_rank,
    comparison_name     = cmp_name
  )

  dataset_results <- list()
  dataset_summary_rows <- list()
  dataset_class_rows <- list()

  for (dataset_key in dataset_key_order) {
    dataset_obj <- dataset_split[[dataset_key]]

    res_df <- run_single_dataset_analysis(
      dataset_obj         = dataset_obj,
      coldata             = comparison_obj$coldata,
      merged_loading_tbl  = merged_loading_tbl,
      comparison_name     = cmp_name,
      dataset_name        = dataset_key
    )

    dataset_results[[dataset_key]] <- res_df

    export_single_dataset_tables(
      res_df          = res_df,
      comparison_name = cmp_name,
      dataset_name    = dataset_key,
      table_dir_local = table_dir
    )

    export_single_dataset_figures(
      res_df           = res_df,
      dataset_obj      = dataset_obj,
      coldata          = comparison_obj$coldata,
      comparison_name  = cmp_name,
      dataset_name     = dataset_key,
      figure_dir_local = figure_dir
    )

    dataset_summary_rows[[dataset_key]] <- build_dataset_summary_table(
      res_df          = res_df,
      comparison_name = cmp_name,
      dataset_name    = dataset_key
    )

    dataset_class_rows[[dataset_key]] <- build_manuscript_class_counts(
      res_df          = res_df,
      comparison_name = cmp_name,
      dataset_name    = dataset_key
    )
  }

  export_comparison_cutoff_figures(
    combined_cutoff_info = cutoff_info,
    comparison_name      = cmp_name,
    figure_dir_local     = figure_dir,
    table_dir_local      = table_dir
  )

  save_csv(
    merged_loading_tbl,
    file.path(table_dir, paste0(cmp_name, "__combined_loading_rank_table.csv"))
  )

  save_csv(
    trt_loading_tbl,
    file.path(table_dir, paste0(cmp_name, "__treatment_loading_table.csv"))
  )

  save_csv(
    ctrl_loading_tbl,
    file.path(table_dir, paste0(cmp_name, "__control_loading_table.csv"))
  )

  list(
    comparison_name       = cmp_name,
    comparison_obj        = comparison_obj,
    merged_loading_tbl    = merged_loading_tbl,
    treatment_loading_tbl = trt_loading_tbl,
    control_loading_tbl   = ctrl_loading_tbl,
    cutoff_info           = cutoff_info,
    dataset_split         = dataset_split,
    dataset_results       = dataset_results,
    dataset_summary_rows  = dataset_summary_rows,
    dataset_class_rows    = dataset_class_rows
  )
}

build_run_manifest <- function(runtime_seconds,
                               comparison_results) {
  cutoff_rows <- lapply(comparison_results, function(x) {
    build_cutoff_manifest_row(
      comparison_name = x$comparison_name,
      cutoff_info     = x$cutoff_info
    )
  })

  cutoff_manifest <- bind_cutoff_manifest_rows(cutoff_rows)

  data.frame(
    cutoff_manifest_present          = nrow(cutoff_manifest) > 0,
    analysis_name                    = analysis_name,
    runtime_readable                 = format_runtime_minutes(runtime_seconds),
    runtime_seconds                  = runtime_seconds,
    deseq2_alpha_level               = alpha_level,
    lfc_boundary                     = lfc_boundary,
    main_cutoff_mode                 = evs_cutoff_mode_main,
    backup_cutoff_mode               = evs_backup_cutoff_mode,
    manual_fixed_rank_override       = if (is.na(evs_fixed_rank_override)) NA_integer_ else as.integer(evs_fixed_rank_override),
    fixed_top_n_backup               = evs_fixed_top_n,
    fourier_percentile_step          = fourier_percentile_step,
    fourier_window_fraction          = fourier_window_fraction,
    fourier_harmonics                = fourier_harmonics,
    crossing_stability_window_n      = crossing_stability_window_n,
    crossing_rule                    = crossing_rule,
    use_group_crossing_bracket       = use_group_crossing_bracket,
    group_bracket_fallback           = group_bracket_fallback,
    n_comparisons                    = length(comparison_results),
    output_dir                       = output_dir,
    stringsAsFactors = FALSE
  )
}

write_readme_summary <- function(comparison_results,
                                 runtime_seconds,
                                 output_file = file.path(output_dir, "README_run_summary.txt")) {
  lines <- c(
    "SEQUENCE manuscript analysis run summary",
    "",
    paste0("Analysis name: ", analysis_name),
    paste0("Runtime: ", format_runtime_minutes(runtime_seconds)),
    paste0("Output directory: ", output_dir),
    "",
    "Core rules",
    paste0("  DESeq2 adjusted p value threshold: ", alpha_level),
    paste0("  Strong effect boundary |LFC|: ", lfc_boundary),
    paste0("  Primary EVS cutoff method: ", evs_cutoff_mode_main),
    paste0("  Backup EVS cutoff method: ", evs_backup_cutoff_mode),
    paste0("  Crossing stability window: ", crossing_stability_window_n),
    "",
    "Comparison cutoffs"
  )

  for (res in comparison_results) {
    ci <- res$cutoff_info
    lines <- c(
      lines,
      paste0(
        "  ", res$comparison_name,
        " | method = ", ci$method,
        " | selected_reason = ", ci$selected_reason,
        " | rank = ", ci$top_n_actual,
        " | quantile = ", signif(ci$cutoff_quantile, 4),
        " | cutoff_value = ", signif(ci$cutoff_value, 6)
      )
    )
  }

  lines <- c(
    lines,
    "",
    "Export families",
    "  figures/",
    "    comparison-level regime crossing panels",
    "    volcano plots",
    "    dispersion plots",
    "    PCA panels",
    "  tables/",
    "    manuscript gene tables",
    "    display class count tables",
    "    dataset summary tables",
    "    cutoff manifest and crossing summaries",
    "    EVS loading rank tables"
  )

  writeLines(lines, con = output_file)
  invisible(output_file)
}

run_sequence_manuscript_pipeline <- function() {
  run_start <- Sys.time()

  imported <- read_sequence_count_matrix(count_file)
  full_count_mat     <- imported$count_mat
  full_annotation_df <- imported$annotation_df

  all_comparison_results <- list()

  for (i in seq_len(nrow(comparison_table))) {
    cmp_row <- comparison_table[i, , drop = FALSE]

    cmp_result <- run_single_comparison_pipeline(
      comparison_row     = cmp_row,
      full_count_mat     = full_count_mat,
      full_annotation_df = full_annotation_df
    )

    all_comparison_results[[cmp_result$comparison_name]] <- cmp_result
  }

  runtime_seconds <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))

  all_dataset_summary_rows <- unlist(
    lapply(all_comparison_results, function(x) x$dataset_summary_rows),
    recursive = FALSE
  )

  all_class_rows <- unlist(
    lapply(all_comparison_results, function(x) x$dataset_class_rows),
    recursive = FALSE
  )

  cutoff_rows <- lapply(all_comparison_results, function(x) {
    build_cutoff_manifest_row(
      comparison_name = x$comparison_name,
      cutoff_info     = x$cutoff_info
    )
  })

  dataset_summary_tbl <- bind_dataset_summary_rows(all_dataset_summary_rows)
  class_count_tbl     <- bind_class_count_rows(all_class_rows)
  cutoff_manifest_tbl <- bind_cutoff_manifest_rows(cutoff_rows)
  run_manifest_tbl    <- build_run_manifest(runtime_seconds, all_comparison_results)

  save_csv(
    dataset_summary_tbl,
    file.path(table_dir, "ALL_COMPARISONS__dataset_summary.csv")
  )

  save_csv(
    class_count_tbl,
    file.path(table_dir, "ALL_COMPARISONS__display_class_counts.csv")
  )

  save_csv(
    cutoff_manifest_tbl,
    file.path(table_dir, "ALL_COMPARISONS__cutoff_manifest.csv")
  )

  save_csv(
    run_manifest_tbl,
    file.path(output_dir, "run_manifest.csv")
  )

  write_readme_summary(
    comparison_results = all_comparison_results,
    runtime_seconds    = runtime_seconds,
    output_file        = file.path(output_dir, "README_run_summary.txt")
  )

  message("")
  message("Analysis complete.")
  message("Output directory: ", output_dir)
  message("Runtime: ", format_runtime_minutes(runtime_seconds))

  invisible(
    list(
      output_dir = output_dir,
      figure_dir = figure_dir,
      table_dir = table_dir,
      runtime_seconds = runtime_seconds,
      comparison_results = all_comparison_results,
      dataset_summary_table = dataset_summary_tbl,
      display_class_count_table = class_count_tbl,
      cutoff_manifest_table = cutoff_manifest_tbl,
      run_manifest_table = run_manifest_tbl
    )
  )
}

run_twas_overlap_summary <- function(all_results_object,
                                     twas_csv = twas_file) {
  if (!isTRUE(run_twas_overlap)) {
    message("TWAS overlap disabled.")
    return(invisible(NULL))
  }

  if (!file.exists(twas_csv)) {
    warning("TWAS file was requested but not found: ", twas_csv, call. = FALSE)
    return(invisible(NULL))
  }

  twas_df <- utils::read.csv(twas_csv, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(twas_df)) {
    warning("TWAS file is empty.", call. = FALSE)
    return(invisible(NULL))
  }

  gene_col <- names(twas_df)[1]
  twas_genes <- unique(as.character(twas_df[[gene_col]]))
  twas_genes <- twas_genes[!is.na(twas_genes) & nzchar(twas_genes)]

  overlap_rows <- list()

  for (cmp_name in names(all_results_object$comparison_results)) {
    cmp_res <- all_results_object$comparison_results[[cmp_name]]

    for (dataset_key in names(cmp_res$dataset_results)) {
      res_df <- cmp_res$dataset_results[[dataset_key]]

      if (!"gene_symbol" %in% names(res_df)) next

      plot_df <- derive_plot_classes(res_df)

      hit_df <- plot_df[
        !is.na(plot_df$gene_symbol) &
          plot_df$gene_symbol %in% twas_genes &
          plot_df$volcano_display_class != "background",
        ,
        drop = FALSE
      ]

      if (!nrow(hit_df)) next

      keep_cols <- intersect(
        c(
          "feature_id",
          "gene_symbol",
          "dataset_name",
          "lfc_shrunk",
          "padj",
          "empirical_p",
          "HBFSS",
          "volcano_display_class",
          "combined_rank"
        ),
        names(hit_df)
      )

      overlap_rows[[paste(cmp_name, dataset_key, sep = "__")]] <- hit_df[, keep_cols, drop = FALSE]
    }
  }

  if (!length(overlap_rows)) {
    message("No TWAS overlap rows were found among non-background calls.")
    return(invisible(NULL))
  }

  overlap_tbl <- dplyr::bind_rows(overlap_rows)
  save_csv(
    overlap_tbl,
    file.path(table_dir, "ALL_COMPARISONS__TWAS_overlap_hits.csv")
  )

  invisible(overlap_tbl)
}

sequence_pipeline_results <- run_sequence_manuscript_pipeline()

if (isTRUE(run_twas_overlap)) {
  run_twas_overlap_summary(sequence_pipeline_results, twas_csv = twas_file)
}
