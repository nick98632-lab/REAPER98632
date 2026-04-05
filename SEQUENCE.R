# =============================================================================
# SEQUENCE STAGE 1
# EVS + LOCAL FOURIER WAVE MAP AT EVERY GENE + EIGENVECTOR SPLITTING
# =============================================================================
# This script does only the Stage 1 workflow:
# 1. Read the raw count matrix.
# 2. Normalize counts with DESeq2.
# 3. Build treatment and control EVS loading tables from PC1 absolute loadings.
# 4. Build local Fourier half-range amplitude maps at every gene rank for NB-derived IOD and CV2.
# 5. Rescale each within-group amplitude trajectory to 0 to 1, then build the treatment-control composite amplitude overlap curve on the shared EVS axis.
# 6. Select the last local half-range amplitude crossing before divergence.
# 7. Perform union-based eigenvector splitting at that cutoff rank.
# 8. Export leading-edge and remainder datasets, tables, and panels.
#
# Interpretation: the local-amplitude crossing marks a transition in the dominant
# oscillatory regime of the EVS-ranked data, such as a switch from IOD-dominant
# to CV2-dominant local structure or the reverse.
#
# This script does not include:
# - percentile gating
# - fallback cutoff logic
# - stability screening
# - higher criticism
# - HBFSS
# - LFC shrinkage
# - volcano plots
# - completed differential expression significance calling
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(S4Vectors)
})

# =============================================================================
# USER SETTINGS
# =============================================================================
# These settings define the analysis scope and the local Fourier model used for
# Stage 1. The goal of this stage is not endpoint significance calling. Rather,
# it is to identify a data-defined structural transition along the EVS-ranked
# axis, then use that transition to split the feature space into a leading edge
# and a remainder for downstream analysis.

repo_dir <- getwd()
input_dir <- file.path(repo_dir, "data")
output_root <- file.path(repo_dir, "exports")
analysis_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(output_root, paste0("sequence_stage1_evs_fourier_clean_", analysis_stamp))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
figure_dpi <- 320
base_theme_size <- 10
fourier_harmonics <- 2L
local_window_fraction <- 0.12
local_percentile_step <- 0.01
min_features_required <- 25L

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

plot_colors <- list(
  iod = "#1F78B4",
  cv2 = "#4D4D4D",
  diff = "#8C2D04",
  hist = "#8C8C8C"
)

# =============================================================================
# BASIC HELPERS
# =============================================================================
# These helper functions keep I/O, plotting, and validation behavior simple and
# explicit so that the methodological code below remains readable.


rescale_to_unit_interval <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  xmin <- min(x[ok])
  xmax <- max(x[ok])
  if (!is.finite(xmin) || !is.finite(xmax)) return(out)
  if (identical(xmax, xmin) || abs(xmax - xmin) < .Machine$double.eps) {
    out[ok] <- 0
    return(out)
  }
  out[ok] <- (x[ok] - xmin) / (xmax - xmin)
  out
}
assert_columns <- function(df, cols, object_name) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      paste0(object_name, " is missing required columns: ", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
}

plain_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88")
    )
}

save_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE)
}

save_plot <- function(p, path, width = 12, height = 8) {
  ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = figure_dpi, bg = "white")
}

save_grob <- function(g, path, width = 14, height = 10) {
  ggplot2::ggsave(path, plot = g, width = width, height = height, dpi = figure_dpi, bg = "white")
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

# =============================================================================
# INPUT AND PREP
# =============================================================================
# The input stage reads the raw count matrix, identifies the feature identifier
# columns, subsets the comparison-specific samples, and generates normalized
# counts for EVS ranking. Group-specific negative-binomial summary quantities are
# then estimated for later construction of the IOD and CV2 local-amplitude maps.

read_count_matrix <- function(path, meta_ids) {
  if (!file.exists(path)) {
    stop(paste0("Count file not found: ", path), call. = FALSE)
  }

  raw_df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  feature_col <- detect_feature_id_column(raw_df)
  symbol_col <- detect_gene_symbol_column(raw_df)
  sample_cols <- intersect(meta_ids, names(raw_df))

  if (length(sample_cols) == 0L) {
    stop("No count columns matched the metadata sample IDs.", call. = FALSE)
  }

  annot_df <- data.frame(
    feature_id = as.character(raw_df[[feature_col]]),
    stringsAsFactors = FALSE
  )
  annot_df$gene_symbol <- if (!is.null(symbol_col)) as.character(raw_df[[symbol_col]]) else annot_df$feature_id
  keep <- !is.na(annot_df$feature_id) & nzchar(annot_df$feature_id)
  annot_df <- annot_df[keep, , drop = FALSE]

  count_mat <- as.matrix(raw_df[keep, sample_cols, drop = FALSE])
  mode(count_mat) <- "numeric"
  rownames(count_mat) <- annot_df$feature_id
  colnames(count_mat) <- sample_cols

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

compute_normalized_counts <- function(count_matrix, coldata) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_matrix)),
    colData = coldata,
    design = ~ condition
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  counts(dds, normalized = TRUE)
}

# Estimate group-specific negative-binomial mean and dispersion terms.
# These quantities define IOD = 1 + alpha*mu and CV2 = 1/mu + alpha at the
# feature level before local oscillatory modeling.
compute_group_feature_metrics <- function(count_submatrix) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_submatrix)),
    colData = S4Vectors::DataFrame(row.names = colnames(count_submatrix)),
    design = ~ 1
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds, quiet = TRUE)

  md <- as.data.frame(SummarizedExperiment::mcols(dds), stringsAsFactors = FALSE)
  md$feature_id <- rownames(md)
  md[, intersect(c("feature_id", "baseMean", "dispGeneEst"), names(md)), drop = FALSE]
}

# =============================================================================
# EVS LOADING TABLES
# =============================================================================
# EVS is represented here through PC1 loading geometry within each group. The
# absolute PC1 loading rank defines the within-group ordering. The treatment and
# control rank tables are then merged onto a shared EVS axis, which is used for
# local-amplitude comparison and final eigenvector splitting.

# Build the within-group EVS ranking table. Features are ordered by the absolute
# magnitude of the first principal-component loading computed from log2(normalized
# count + 1) data.
compute_pc1_loading_table <- function(norm_counts, sample_ids, feature_metrics, annot_df) {
  x <- log2(as.matrix(norm_counts[, sample_ids, drop = FALSE]) + 1)
  pca_fit <- prcomp(t(x), scale. = FALSE, rank. = 2)
  loading_vec <- pca_fit$rotation[, 1]

  loading_tbl <- data.frame(
    feature_id = names(loading_vec),
    pc1_loading = as.numeric(loading_vec),
    pc1_loading_abs = abs(as.numeric(loading_vec)),
    stringsAsFactors = FALSE
  )

  loading_tbl <- loading_tbl %>%
    left_join(feature_metrics, by = "feature_id") %>%
    left_join(annot_df, by = "feature_id") %>%
    arrange(desc(pc1_loading_abs), feature_id)

  loading_tbl$rank <- seq_len(nrow(loading_tbl))
  loading_tbl
}

build_shared_evs_table <- function(trt_loading_tbl, ctrl_loading_tbl) {
  trt_use <- trt_loading_tbl[, c("feature_id", "gene_symbol", "pc1_loading", "pc1_loading_abs", "rank"), drop = FALSE]
  ctrl_use <- ctrl_loading_tbl[, c("feature_id", "gene_symbol", "pc1_loading", "pc1_loading_abs", "rank"), drop = FALSE]

  names(trt_use) <- c("feature_id", "gene_symbol_trt", "pc1_loading_trt", "pc1_loading_abs_trt", "rank_trt")
  names(ctrl_use) <- c("feature_id", "gene_symbol_ctrl", "pc1_loading_ctrl", "pc1_loading_abs_ctrl", "rank_ctrl")

  merged <- full_join(trt_use, ctrl_use, by = "feature_id")
  merged$gene_symbol <- ifelse(
    !is.na(merged$gene_symbol_trt) & nzchar(merged$gene_symbol_trt),
    merged$gene_symbol_trt,
    merged$gene_symbol_ctrl
  )

  merged <- merged %>%
    mutate(
      rank_trt = as.integer(rank_trt),
      rank_ctrl = as.integer(rank_ctrl),
      combined_rank = round(rowMeans(cbind(rank_trt, rank_ctrl), na.rm = TRUE))
    ) %>%
    arrange(combined_rank)

  merged
}

# =============================================================================
# LOCAL FOURIER AMPLITUDE MAPS AT EVERY GENE
# =============================================================================
# This is the core structural measurement stage. For each EVS rank, the script
# fits a local truncated Fourier model within a neighborhood centered at that
# rank. The resulting local amplitude is interpreted as the strength of the
# local oscillatory regime for the metric of interest. This is done separately
# for IOD and CV2 within treatment and control.

build_fourier_design <- function(x, n_harmonics = fourier_harmonics) {
  x <- as.numeric(x)
  x01 <- (x - min(x)) / max(1e-12, (max(x) - min(x)))
  out <- data.frame(x = x01)
  for (k in seq_len(n_harmonics)) {
    out[[paste0("sin_", k)]] <- sin(2 * pi * k * x01)
    out[[paste0("cos_", k)]] <- cos(2 * pi * k * x01)
  }
  out
}

# Fit a one-harmonic local Fourier model inside one window and summarize the
# fitted oscillation by its midpoint and half-range amplitude.
fit_local_fourier_window_summary <- function(rank_vec, y_vec, n_harmonics = fourier_harmonics) {
  rank_vec <- as.numeric(rank_vec)
  y_vec <- as.numeric(y_vec)

  keep <- is.finite(rank_vec) & !is.na(rank_vec) & is.finite(y_vec) & !is.na(y_vec)
  rank_vec <- rank_vec[keep]
  y_vec <- y_vec[keep]

  if (length(rank_vec) < (2 * n_harmonics + 5L)) return(NULL)

  y_sd <- suppressWarnings(stats::sd(y_vec, na.rm = TRUE))
  if (!is.finite(y_sd) || is.na(y_sd) || y_sd == 0) return(NULL)

  dd <- build_fourier_design(rank_vec, n_harmonics = n_harmonics)
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

build_percentile_windows <- function(n_total,
                                     step = local_percentile_step,
                                     window_fraction = local_window_fraction) {
  pct_grid <- seq(step, 1, by = step)
  center_ranks <- pmax(1L, pmin(n_total, round(pct_grid * n_total)))

  half_window <- max(5L, round((window_fraction * n_total) / 2))

  data.frame(
    percentile = pct_grid,
    center_rank = center_ranks,
    lo_rank = pmax(1L, center_ranks - half_window),
    hi_rank = pmin(n_total, center_ranks + half_window),
    stringsAsFactors = FALSE
  )
}

summarize_percentile_window <- function(metric_df, lo_rank, hi_rank) {
  sub <- metric_df[metric_df$rank >= lo_rank & metric_df$rank <= hi_rank, , drop = FALSE]
  if (nrow(sub) < 9L) {
    return(data.frame(
      iod_local_amplitude = NA_real_,
      cv2_local_amplitude = NA_real_,
      iod_local_midpoint = NA_real_,
      cv2_local_midpoint = NA_real_,
      iod_peak_rank = NA_integer_,
      iod_trough_rank = NA_integer_,
      cv2_peak_rank = NA_integer_,
      cv2_trough_rank = NA_integer_,
      iod_center = NA_real_,
      cv2_center = NA_real_,
      iod_residual_sd = NA_real_,
      cv2_residual_sd = NA_real_,
      local_window_n = nrow(sub),
      stringsAsFactors = FALSE
    ))
  }

  iod_fit <- fit_local_fourier_window_summary(sub$rank, sub$log_iod_nb)
  cv2_fit <- fit_local_fourier_window_summary(sub$rank, sub$log_cv2_nb)
  if (is.null(iod_fit) || is.null(cv2_fit)) {
    return(data.frame(
      iod_local_amplitude = NA_real_,
      cv2_local_amplitude = NA_real_,
      iod_local_midpoint = NA_real_,
      cv2_local_midpoint = NA_real_,
      iod_peak_rank = NA_integer_,
      iod_trough_rank = NA_integer_,
      cv2_peak_rank = NA_integer_,
      cv2_trough_rank = NA_integer_,
      iod_center = NA_real_,
      cv2_center = NA_real_,
      iod_residual_sd = NA_real_,
      cv2_residual_sd = NA_real_,
      local_window_n = nrow(sub),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    iod_local_amplitude = iod_fit$local_amplitude,
    cv2_local_amplitude = cv2_fit$local_amplitude,
    iod_local_midpoint = 0.5 * (max(iod_fit$fitted_y, na.rm = TRUE) + min(iod_fit$fitted_y, na.rm = TRUE)),
    cv2_local_midpoint = 0.5 * (max(cv2_fit$fitted_y, na.rm = TRUE) + min(cv2_fit$fitted_y, na.rm = TRUE)),
    iod_peak_rank = iod_fit$peak_rank,
    iod_trough_rank = iod_fit$trough_rank,
    cv2_peak_rank = cv2_fit$peak_rank,
    cv2_trough_rank = cv2_fit$trough_rank,
    iod_center = iod_fit$center_rank,
    cv2_center = cv2_fit$center_rank,
    iod_residual_sd = iod_fit$residual_sd,
    cv2_residual_sd = cv2_fit$residual_sd,
    local_window_n = nrow(sub),
    stringsAsFactors = FALSE
  )
}

build_group_wave_map <- function(loading_tbl) {
  assert_columns(loading_tbl, c("feature_id", "rank", "baseMean", "dispGeneEst"), "loading_tbl")

  df <- loading_tbl %>%
    filter(is.finite(baseMean), baseMean > 0, is.finite(dispGeneEst), dispGeneEst > 0) %>%
    arrange(rank)

  if (nrow(df) < min_features_required) {
    stop("Too few analyzable features after filtering for baseMean and dispGeneEst.", call. = FALSE)
  }

  mu <- pmax(df$baseMean, 1e-12)
  alpha <- pmax(df$dispGeneEst, 1e-12)
  df$iod_nb <- 1 + alpha * mu
  df$cv2_nb <- (1 / mu) + alpha
  df$log_iod_nb <- log10(df$iod_nb)
  df$log_cv2_nb <- log10(df$cv2_nb)

  windows <- build_percentile_windows(nrow(df), step = local_percentile_step, window_fraction = local_window_fraction)
  rows <- lapply(seq_len(nrow(windows)), function(i) {
    ww <- windows[i, , drop = FALSE]
    ss <- summarize_percentile_window(df, ww$lo_rank, ww$hi_rank)
    cbind(ww, ss, stringsAsFactors = FALSE)
  })
  wave_map <- dplyr::bind_rows(rows)

  wave_map$iod_local_amplitude_scaled <- if (all(is.na(wave_map$iod_local_amplitude))) {
    rep(NA_real_, nrow(wave_map))
  } else {
    scales::rescale(wave_map$iod_local_amplitude, to = c(0, 1), from = range(wave_map$iod_local_amplitude, na.rm = TRUE))
  }

  wave_map$cv2_local_amplitude_scaled <- if (all(is.na(wave_map$cv2_local_amplitude))) {
    rep(NA_real_, nrow(wave_map))
  } else {
    scales::rescale(wave_map$cv2_local_amplitude, to = c(0, 1), from = range(wave_map$cv2_local_amplitude, na.rm = TRUE))
  }

  wave_map$center_distance <- abs(wave_map$iod_center - wave_map$cv2_center)
  wave_map$center_agreement <- 1 / (1 + wave_map$center_distance)

  wave_map
}

# =============================================================================
# CROSSING AND EVS SPLITTING
# =============================================================================
# The treatment and control local-amplitude maps are first rescaled to the unit
# interval within each group. The treatment-control comparison is then performed
# on these 0 to 1 amplitude trajectories so that the crossing identifies a
# change in relative local oscillatory dominance rather than raw-magnitude scale.
# The selected cutoff is the last local-amplitude crossing before persistent
# divergence. That rank is then projected back onto the treatment and control
# EVS tables, and the leading edge is defined as the union of all features
# retained by either side at or above that EVS threshold.

build_composite_wave_map <- function(trt_wave_map, ctrl_wave_map) {
  trt_use <- trt_wave_map[, c("percentile", "center_rank", "iod_local_amplitude", "cv2_local_amplitude", "iod_local_amplitude_scaled", "cv2_local_amplitude_scaled", "iod_center", "cv2_center"), drop = FALSE]
  ctrl_use <- ctrl_wave_map[, c("percentile", "center_rank", "iod_local_amplitude", "cv2_local_amplitude", "iod_local_amplitude_scaled", "cv2_local_amplitude_scaled", "iod_center", "cv2_center"), drop = FALSE]
  names(trt_use) <- c("percentile", "center_rank_trt", "iod_local_amplitude_trt", "cv2_local_amplitude_trt", "iod_local_amplitude_scaled_trt", "cv2_local_amplitude_scaled_trt", "iod_center_trt", "cv2_center_trt")
  names(ctrl_use) <- c("percentile", "center_rank_ctrl", "iod_local_amplitude_ctrl", "cv2_local_amplitude_ctrl", "iod_local_amplitude_scaled_ctrl", "cv2_local_amplitude_scaled_ctrl", "iod_center_ctrl", "cv2_center_ctrl")

  out <- inner_join(trt_use, ctrl_use, by = "percentile") %>%
    mutate(
      combined_rank = round((center_rank_trt + center_rank_ctrl) / 2),
      composite_iod_amplitude = iod_local_amplitude_trt + iod_local_amplitude_ctrl,
      composite_cv2_amplitude = cv2_local_amplitude_trt + cv2_local_amplitude_ctrl,
      composite_iod_amplitude_scaled = iod_local_amplitude_scaled_trt + iod_local_amplitude_scaled_ctrl,
      composite_cv2_amplitude_scaled = cv2_local_amplitude_scaled_trt + cv2_local_amplitude_scaled_ctrl,
      combined_center_distance = abs(((iod_center_trt + iod_center_ctrl) / 2) - ((cv2_center_trt + cv2_center_ctrl) / 2)),
      combined_center_agreement = 1 / (1 + combined_center_distance),
      regime_difference = composite_iod_amplitude - composite_cv2_amplitude,
      regime_difference_scaled = composite_iod_amplitude_scaled - composite_cv2_amplitude_scaled
    ) %>%
    arrange(percentile)

  out
}

# Find all zero crossings of the composite local-amplitude difference function
# D(r) = A_IOD(r) - A_CV2(r) on the shared EVS axis.
find_all_crossings <- function(composite_wave_map) {
  df <- composite_wave_map %>%
    filter(is.finite(percentile), is.finite(combined_rank), is.finite(regime_difference)) %>%
    arrange(percentile)

  if (nrow(df) < 2L) {
    stop("Not enough points to evaluate crossings on the local wave map.", call. = FALSE)
  }

  out <- vector("list", nrow(df) - 1L)
  idx <- 1L
  for (i in seq_len(nrow(df) - 1L)) {
    y1 <- df$regime_difference[i]
    y2 <- df$regime_difference[i + 1L]
    x1 <- df$percentile[i]
    x2 <- df$percentile[i + 1L]
    r1 <- df$combined_rank[i]
    r2 <- df$combined_rank[i + 1L]

    crossed <- (y1 == 0) || (y2 == 0) || ((y1 > 0) && (y2 < 0)) || ((y1 < 0) && (y2 > 0))
    if (!crossed) next

    crossing_percentile <- if (isTRUE(all.equal(y1, y2))) mean(c(x1, x2)) else x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
    crossing_rank <- if (isTRUE(all.equal(y1, y2))) round(mean(c(r1, r2))) else round(r1 + (0 - y1) * (r2 - r1) / (y2 - y1))

    out[[idx]] <- data.frame(
      crossing_id = paste0("crossing_", idx),
      percentile_left = x1,
      percentile_right = x2,
      rank_left = r1,
      rank_right = r2,
      crossing_percentile = crossing_percentile,
      crossing_rank = as.integer(crossing_rank),
      regime_difference_left = y1,
      regime_difference_right = y2,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }

  out <- out[seq_len(idx - 1L)]
  if (!length(out)) {
    return(data.frame(
      crossing_id = character(), percentile_left = numeric(), percentile_right = numeric(), rank_left = integer(), rank_right = integer(), crossing_percentile = numeric(), crossing_rank = integer(), regime_difference_left = numeric(), regime_difference_right = numeric(), stringsAsFactors = FALSE
    ))
  }

  bind_rows(out) %>% arrange(crossing_percentile)
}

# Select the rightmost crossing. This is the operational definition of the last
# local-amplitude crossing before divergence used for Stage 1 splitting.
select_last_crossing_before_divergence <- function(crossing_tbl) {
  if (!nrow(crossing_tbl)) {
    return(data.frame(
      crossing_id = NA_character_,
      percentile_left = NA_real_,
      percentile_right = NA_real_,
      rank_left = NA_integer_,
      rank_right = NA_integer_,
      crossing_percentile = NA_real_,
      crossing_rank = NA_integer_,
      regime_difference_left = NA_real_,
      regime_difference_right = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  crossing_tbl %>% arrange(desc(crossing_percentile)) %>% slice(1)
}

# Project the selected shared cutoff back onto the treatment and control EVS
# rank tables, then define the leading edge as the union of retained treatment
# and control features. All remaining features are assigned to the remainder.
build_eigenvector_split <- function(shared_evs_tbl, rank_cutoff) {
  membership_tbl <- shared_evs_tbl %>%
    mutate(
      in_treatment_leading_edge = !is.na(rank_trt) & rank_trt <= rank_cutoff,
      in_control_leading_edge = !is.na(rank_ctrl) & rank_ctrl <= rank_cutoff,
      in_leading_edge_union = in_treatment_leading_edge | in_control_leading_edge,
      in_remainder = !in_leading_edge_union
    )

  leading_edge_ids <- as.character(membership_tbl$feature_id[membership_tbl$in_leading_edge_union])
  remainder_ids <- as.character(membership_tbl$feature_id[membership_tbl$in_remainder])

  if (!length(leading_edge_ids)) {
    stop("Leading edge is empty after eigenvector splitting.", call. = FALSE)
  }

  list(
    membership_tbl = membership_tbl,
    leading_edge_ids = leading_edge_ids,
    remainder_ids = remainder_ids,
    leading_edge_union_n = length(leading_edge_ids),
    remainder_union_n = length(remainder_ids),
    n_total = nrow(membership_tbl)
  )
}

run_stage1_method <- function(trt_loading_tbl, ctrl_loading_tbl, comparison_name = "") {
  cat("Building treatment local Fourier wave map", if (nzchar(comparison_name)) paste0(" for ", comparison_name) else "", "...\n", sep = "")
  trt_wave_map <- build_group_wave_map(trt_loading_tbl)
  cat("Building control local Fourier wave map", if (nzchar(comparison_name)) paste0(" for ", comparison_name) else "", "...\n", sep = "")
  ctrl_wave_map <- build_group_wave_map(ctrl_loading_tbl)
  cat("Building shared EVS table", if (nzchar(comparison_name)) paste0(" for ", comparison_name) else "", "...\n", sep = "")
  shared_evs_tbl <- build_shared_evs_table(trt_loading_tbl, ctrl_loading_tbl)
  cat("Building composite local Fourier overlap", if (nzchar(comparison_name)) paste0(" for ", comparison_name) else "", "...\n", sep = "")
  composite_wave_map <- build_composite_wave_map(trt_wave_map, ctrl_wave_map)
  cat("Finding crossings", if (nzchar(comparison_name)) paste0(" for ", comparison_name) else "", "...\n", sep = "")
  crossing_tbl <- find_all_crossings(composite_wave_map)
  selected_crossing <- select_last_crossing_before_divergence(crossing_tbl)

  if (!nrow(crossing_tbl)) {
    return(list(
      trt_wave_map = trt_wave_map,
      ctrl_wave_map = ctrl_wave_map,
      shared_evs_tbl = shared_evs_tbl,
      composite_wave_map = composite_wave_map,
      crossing_tbl = crossing_tbl,
      selected_crossing = selected_crossing,
      rank_cutoff = NA_integer_,
      split_obj = NULL,
      leading_edge_ids = character(),
      remainder_ids = character(),
      leading_edge_union_n = NA_integer_,
      remainder_union_n = NA_integer_,
      n_total = nrow(shared_evs_tbl),
      selected_reason = "no_crossing_detected"
    ))
  }

  rank_cutoff <- max(1L, as.integer(selected_crossing$crossing_rank[1]))
  split_obj <- build_eigenvector_split(shared_evs_tbl, rank_cutoff)

  list(
    trt_wave_map = trt_wave_map,
    ctrl_wave_map = ctrl_wave_map,
    shared_evs_tbl = shared_evs_tbl,
    composite_wave_map = composite_wave_map,
    crossing_tbl = crossing_tbl,
    selected_crossing = selected_crossing,
    rank_cutoff = rank_cutoff,
    split_obj = split_obj,
    leading_edge_ids = split_obj$leading_edge_ids,
    remainder_ids = split_obj$remainder_ids,
    leading_edge_union_n = split_obj$leading_edge_union_n,
    remainder_union_n = split_obj$remainder_union_n,
    n_total = split_obj$n_total,
    selected_reason = "last_local_amplitude_crossing_before_divergence"
  )
}

# =============================================================================
# EXPORT HELPERS
# =============================================================================
# Export helpers write manuscript-readable tables and preserve the split
# datasets needed for later downstream analysis.

subset_matrix_by_ids <- function(mat, ids) {
  keep_ids <- intersect(ids, rownames(mat))
  mat[keep_ids, , drop = FALSE]
}

matrix_to_export_table <- function(mat, annot_df) {
  out <- as.data.frame(mat, stringsAsFactors = FALSE)
  out$feature_id <- rownames(out)
  out <- out %>% left_join(annot_df, by = "feature_id")
  out[, c("feature_id", "gene_symbol", setdiff(names(out), c("feature_id", "gene_symbol"))), drop = FALSE]
}

build_stage1_summary_table <- function(stage1_obj, comparison_name) {
  sc <- stage1_obj$selected_crossing
  data.frame(
    comparison_name = comparison_name,
    selected_reason = stage1_obj$selected_reason,
    crossing_id = sc$crossing_id[1],
    crossing_rank = sc$crossing_rank[1],
    rank_cutoff = stage1_obj$rank_cutoff,
    n_total = stage1_obj$n_total,
    leading_edge_union_n = stage1_obj$leading_edge_union_n,
    remainder_union_n = stage1_obj$remainder_union_n,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# PLOTS
# =============================================================================
# The figures are intended to make the method visually interpretable. They show
# within-group local amplitudes, the shared treatment-control overlap, the
# difference curve whose zero identifies the crossing, and the EVS geometry used
# for the final union-based split.

plot_group_wave_map <- function(group_wave_map, group_label, comparison_name) {
  ggplot(group_wave_map, aes(percentile)) +
    geom_line(aes(y = iod_local_amplitude_scaled, color = "IOD"), linewidth = 0.9) +
    geom_line(aes(y = cv2_local_amplitude_scaled, color = "CV²"), linewidth = 0.9) +
    scale_color_manual(values = c("IOD" = plot_colors$iod, "CV²" = plot_colors$cv2)) +
    labs(
      title = paste0(comparison_name, " | ", group_label, " local Fourier amplitude map"),
      subtitle = "Evaluated on a 0.5% percentile-center grid with 12% window width; amplitudes were rescaled to 0 to 1 within group",
      x = "Percentile center",
      y = "Scaled local amplitude (0-1)"
    ) +
    plain_theme()
}

plot_composite_overlap <- function(stage1_obj, comparison_name) {
  df <- stage1_obj$composite_wave_map
  sc <- stage1_obj$selected_crossing
  ymax <- max(c(df$composite_iod_amplitude, df$composite_cv2_amplitude), na.rm = TRUE)

  p <- ggplot(df, aes(percentile)) +
    geom_line(aes(y = composite_iod_amplitude, color = "Composite IOD"), linewidth = 0.9) +
    geom_line(aes(y = composite_cv2_amplitude, color = "Composite CV²"), linewidth = 0.9) +
    scale_color_manual(values = c("Composite IOD" = plot_colors$iod, "Composite CV²" = plot_colors$cv2)) +
    labs(
      title = paste0(comparison_name, " | treatment-control local amplitude overlap"),
      subtitle = "Composite amplitudes on the shared 0.5% percentile-center grid after within-group 0 to 1 rescaling",
      x = "Percentile center",
      y = "Composite local amplitude"
    ) +
    plain_theme()

  if (isTRUE(stage1_obj$selected_reason == "last_local_amplitude_crossing_before_divergence") && is.finite(sc$crossing_rank[1])) {
    p <- p +
      geom_vline(xintercept = approx(x = df$combined_rank, y = df$percentile, xout = sc$crossing_rank[1], ties = "ordered")$y, linetype = "dashed", linewidth = 0.9, colour = plot_colors$diff) +
      annotate(
        "label",
        x = approx(x = df$combined_rank, y = df$percentile, xout = sc$crossing_rank[1], ties = "ordered")$y,
        y = ymax,
        label = paste0(
          "Last local-amplitude crossing before divergence
",
          "rank cutoff = ", stage1_obj$rank_cutoff, "
",
          "leading edge n = ", stage1_obj$leading_edge_union_n, "
",
          "remainder n = ", stage1_obj$remainder_union_n
        ),
        fill = "white",
        colour = plot_colors$diff,
        size = 3,
        vjust = -0.5
      )
  } else {
    p <- p +
      annotate(
        "label",
        x = median(df$percentile, na.rm = TRUE),
        y = ymax,
        label = "No crossing detected
Diagnostic plot only",
        fill = "white",
        colour = plot_colors$diff,
        size = 3,
        vjust = -0.5
      )
  }

  p
}

plot_difference_curve <- function(stage1_obj, comparison_name) {
  df <- stage1_obj$composite_wave_map
  sc <- stage1_obj$selected_crossing

  p <- ggplot(df, aes(percentile, regime_difference)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.5) +
    geom_line(colour = plot_colors$diff, linewidth = 0.9) +
    labs(
      title = paste0(comparison_name, " | local half-range amplitude difference (IOD minus CV²)"),
      subtitle = "Computed on the 0.5% percentile-center grid; the selected cutoff is the last zero crossing before persistent divergence",
      x = "Percentile center",
      y = "IOD amplitude − CV² amplitude"
    ) +
    plain_theme()

  if (isTRUE(stage1_obj$selected_reason == "last_local_amplitude_crossing_before_divergence") && is.finite(sc$crossing_rank[1])) {
    p <- p +
      geom_vline(xintercept = approx(x = df$combined_rank, y = df$percentile, xout = sc$crossing_rank[1], ties = "ordered")$y, linetype = "dashed", linewidth = 0.9, colour = plot_colors$diff) +
      geom_point(data = data.frame(percentile = approx(x = df$combined_rank, y = df$percentile, xout = sc$crossing_rank[1], ties = "ordered")$y, regime_difference = 0), aes(percentile, regime_difference), inherit.aes = FALSE, size = 2) +
      annotate(
        "label",
        x = approx(x = df$combined_rank, y = df$percentile, xout = sc$crossing_rank[1], ties = "ordered")$y,
        y = 0,
        label = paste0("rank cutoff = ", stage1_obj$rank_cutoff),
        fill = "white",
        colour = plot_colors$diff,
        size = 3,
        vjust = -0.8
      )
  } else {
    p <- p +
      annotate(
        "label",
        x = median(df$percentile, na.rm = TRUE),
        y = 0,
        label = "No crossing detected",
        fill = "white",
        colour = plot_colors$diff,
        size = 3,
        vjust = -0.8
      )
  }

  p
}

plot_loading_rank_curve <- function(loading_tbl, rank_cutoff, group_label, comparison_name) {
  ggplot(loading_tbl, aes(rank, pc1_loading_abs)) +
    geom_line(colour = plot_colors$hist, linewidth = 0.5) +
    geom_vline(xintercept = rank_cutoff, linetype = "dashed", linewidth = 0.9, colour = plot_colors$diff) +
    labs(
      title = paste0(comparison_name, " | ", group_label, " EVS loading rank curve"),
      subtitle = "Absolute PC1 loading defines the within-group EVS rank",
      x = "EVS rank",
      y = "Absolute PC1 loading"
    ) +
    plain_theme()
}

plot_loading_histogram <- function(loading_tbl, group_label, comparison_name) {
  ggplot(loading_tbl, aes(pc1_loading_abs)) +
    geom_histogram(bins = 60, fill = plot_colors$hist, colour = "white") +
    labs(
      title = paste0(comparison_name, " | ", group_label, " EVS loading histogram"),
      subtitle = "Distribution of absolute PC1 loading magnitudes",
      x = "Absolute PC1 loading",
      y = "Count"
    ) +
    plain_theme()
}

build_fourier_panel <- function(stage1_obj, comparison_name) {
  p1 <- plot_group_wave_map(stage1_obj$trt_wave_map, "treatment", comparison_name)
  p2 <- plot_group_wave_map(stage1_obj$ctrl_wave_map, "control", comparison_name)
  p3 <- plot_composite_overlap(stage1_obj, comparison_name)
  p4 <- plot_difference_curve(stage1_obj, comparison_name)

  gridExtra::arrangeGrob(
    grobs = list(p1, p2, p3, p4),
    ncol = 2,
    top = grid::textGrob(
      paste0(comparison_name, " | local Fourier amplitude panels and crossing summary"),
      gp = grid::gpar(fontface = "bold", cex = 1.05)
    )
  )
}

build_evs_panel <- function(trt_loading_tbl, ctrl_loading_tbl, rank_cutoff, comparison_name) {
  p1 <- plot_loading_rank_curve(trt_loading_tbl, rank_cutoff, "treatment", comparison_name)
  p2 <- plot_loading_rank_curve(ctrl_loading_tbl, rank_cutoff, "control", comparison_name)
  p3 <- plot_loading_histogram(trt_loading_tbl, "treatment", comparison_name)
  p4 <- plot_loading_histogram(ctrl_loading_tbl, "control", comparison_name)

  gridExtra::arrangeGrob(
    grobs = list(p1, p2, p3, p4),
    ncol = 2,
    top = grid::textGrob(
      paste0(comparison_name, " | EVS geometry panels"),
      gp = grid::gpar(fontface = "bold", cex = 1.05)
    )
  )
}

# =============================================================================
# COMPARISON RUNNER
# =============================================================================

run_one_comparison <- function(comparison_row, count_matrix, annot_df) {
  comparison_name <- comparison_row$comparison_name
  cat("\n--- Running ", comparison_name, " ---\n", sep = "")

  comp <- subset_comparison(count_matrix, comparison_row, meta_all)
  cat("Estimating comparison-normalized counts...\n")
  norm_counts <- compute_normalized_counts(comp$count_matrix, comp$coldata)
  cat("Estimating treatment NB metrics...\n")
  trt_metrics <- compute_group_feature_metrics(comp$count_matrix[, comp$trt_ids, drop = FALSE])
  cat("Estimating control NB metrics...\n")
  ctrl_metrics <- compute_group_feature_metrics(comp$count_matrix[, comp$ctrl_ids, drop = FALSE])
  trt_loading_tbl <- compute_pc1_loading_table(norm_counts, comp$trt_ids, trt_metrics, annot_df)
  ctrl_loading_tbl <- compute_pc1_loading_table(norm_counts, comp$ctrl_ids, ctrl_metrics, annot_df)

  stage1_obj <- run_stage1_method(trt_loading_tbl, ctrl_loading_tbl, comparison_name = comparison_name)

  comparison_dir <- file.path(output_dir, comparison_name)
  table_dir <- file.path(comparison_dir, "tables")
  figure_dir <- file.path(comparison_dir, "figures")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  save_csv(trt_loading_tbl, file.path(table_dir, paste0(comparison_name, "_treatment_loading_table.csv")))
  save_csv(ctrl_loading_tbl, file.path(table_dir, paste0(comparison_name, "_control_loading_table.csv")))
  save_csv(stage1_obj$trt_wave_map, file.path(table_dir, paste0(comparison_name, "_treatment_local_fourier_amplitude_map.csv")))
  save_csv(stage1_obj$ctrl_wave_map, file.path(table_dir, paste0(comparison_name, "_control_local_fourier_amplitude_map.csv")))
  save_csv(stage1_obj$shared_evs_tbl, file.path(table_dir, paste0(comparison_name, "_shared_evs_table.csv")))
  save_csv(stage1_obj$composite_wave_map, file.path(table_dir, paste0(comparison_name, "_composite_local_fourier_amplitude_map.csv")))
  save_csv(stage1_obj$crossing_tbl, file.path(table_dir, paste0(comparison_name, "_crossing_table.csv")))
  save_csv(build_stage1_summary_table(stage1_obj, comparison_name), file.path(table_dir, paste0(comparison_name, "_stage1_summary.csv")))

  save_plot(plot_group_wave_map(stage1_obj$trt_wave_map, "treatment", comparison_name), file.path(figure_dir, paste0(comparison_name, "_treatment_local_fourier_amplitude_map.png")))
  save_plot(plot_group_wave_map(stage1_obj$ctrl_wave_map, "control", comparison_name), file.path(figure_dir, paste0(comparison_name, "_control_local_fourier_amplitude_map.png")))
  save_plot(plot_composite_overlap(stage1_obj, comparison_name), file.path(figure_dir, paste0(comparison_name, "_treatment_control_local_amplitude_overlap.png")))
  save_plot(plot_difference_curve(stage1_obj, comparison_name), file.path(figure_dir, paste0(comparison_name, "_iod_minus_cv2_amplitude_curve.png")))
  save_grob(build_fourier_panel(stage1_obj, comparison_name), file.path(figure_dir, paste0(comparison_name, "_fourier_panel.png")))
  save_grob(build_evs_panel(trt_loading_tbl, ctrl_loading_tbl, stage1_obj$rank_cutoff, comparison_name), file.path(figure_dir, paste0(comparison_name, "_evs_panel.png")))

  if (!is.null(stage1_obj$split_obj)) {
    leading_raw <- subset_matrix_by_ids(comp$count_matrix, stage1_obj$leading_edge_ids)
    remainder_raw <- subset_matrix_by_ids(comp$count_matrix, stage1_obj$remainder_ids)
    leading_norm <- subset_matrix_by_ids(norm_counts, stage1_obj$leading_edge_ids)
    remainder_norm <- subset_matrix_by_ids(norm_counts, stage1_obj$remainder_ids)

    save_csv(stage1_obj$split_obj$membership_tbl, file.path(table_dir, paste0(comparison_name, "_split_membership_table.csv")))
    save_csv(matrix_to_export_table(leading_raw, annot_df), file.path(table_dir, paste0(comparison_name, "_leading_edge_raw_counts.csv")))
    save_csv(matrix_to_export_table(remainder_raw, annot_df), file.path(table_dir, paste0(comparison_name, "_remainder_raw_counts.csv")))
    save_csv(matrix_to_export_table(leading_norm, annot_df), file.path(table_dir, paste0(comparison_name, "_leading_edge_normalized_counts.csv")))
    save_csv(matrix_to_export_table(remainder_norm, annot_df), file.path(table_dir, paste0(comparison_name, "_remainder_normalized_counts.csv")))
    if (!length(stage1_obj$remainder_ids)) {
      save_csv(data.frame(note = "Selected cutoff places all features in the leading-edge union; remainder export is intentionally empty.", stringsAsFactors = FALSE), file.path(table_dir, paste0(comparison_name, "_remainder_status.csv")))
    }
  } else {
    save_csv(data.frame(note = "No local-amplitude crossing detected; split datasets were not created.", stringsAsFactors = FALSE), file.path(table_dir, paste0(comparison_name, "_split_status.csv")))
  }

  build_stage1_summary_table(stage1_obj, comparison_name)
}

# =============================================================================
# DRIVER
# =============================================================================

main <- function() {
  loaded <- read_count_matrix(count_file, meta_all$id)
  count_matrix <- loaded$count_matrix
  annot_df <- loaded$annot_df

  all_summary <- bind_rows(lapply(seq_len(nrow(comparison_table)), function(i) {
    run_one_comparison(comparison_table[i, , drop = FALSE], count_matrix, annot_df)
  }))

  save_csv(all_summary, file.path(output_dir, "all_comparisons_stage1_summary.csv"))
  cat("\nCompleted Stage 1 EVS + local Fourier amplitude rewrite. Output directory:\n", output_dir, "\n", sep = "")
}

main()
