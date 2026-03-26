# =============================================================================
# PART 1 OF 4
# SETUP, GLOBAL CONTROLS, METADATA, IMPORT, AND CORE HELPERS
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
# GLOBAL CONTROLS
# -----------------------------------------------------------------------------

fast_mode <- TRUE
compute_optional_diagnostics <- FALSE
use_fast_apeglm <- TRUE
export_all_plots <- TRUE

n_workers <- max(
  1L,
  min(4L, parallel::detectCores(logical = FALSE) - 1L)
)

output_root <- "exports"
analysis_name <- "EVS_HBFSS_AllComparisons_Output"

input_dir <- "/root/REAPER98632/data"
count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")

output_dir <- file.path(output_root, analysis_name)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

DESIGN_FORMULA <- ~ condition

alpha_level <- 0.10
lfc_boundary <- 1.0
max_usable_hc_p_threshold <- 0.95

figure_dpi <- 320
base_theme_size <- 10

# -----------------------------------------------------------------------------
# EVS SEARCH CONTROLS
# -----------------------------------------------------------------------------

evs_cutoff_mode_main <- "shared_rank_exhaustive_search"
evs_candidate_max_leading_frac <- 1.00
evs_top_k_ranks_per_comparison <- 4L

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
# FAST ROW VARIANCE
# -----------------------------------------------------------------------------

fast_row_vars <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::rowVars(mat, na.rm = TRUE))
  }

  apply(mat, 1L, stats::var, na.rm = TRUE)
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
    na.last = NA,
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

dataset_key_labels <- c(
  raw_dataset          = "Original dataset",
  leading_edge_dataset = "Leading-edge dataset",
  remainder_dataset    = "Remainder dataset"
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
# PART 2 OF 4
# EVS HELPERS, SHARED-RANK SEARCH, AND GLOBAL PRECOMPUTE
# =============================================================================

compute_pc1_loading_table <- function(value_df, sample_names,
                                      preprocessing_label = "Normalized before EVS",
                                      selected_reason = "shared_rank_external") {
  x <- log1p(as.matrix(value_df[, sample_names, drop = FALSE]))
  max_pcs <- min(5L, ncol(x))
  pca_fit <- prcomp(t(x), scale. = FALSE, rank. = max_pcs)

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
    cutoff_method = "shared_rank_exhaustive_search",
    selected_reason = selected_reason,
    preprocessing_label = preprocessing_label
  )
}

build_rank_prefix_cache <- function(count_matrix, fit_obj,
                                    max_leading_frac = evs_candidate_max_leading_frac) {
  loading_tbl <- as.data.frame(fit_obj$loading_table, stringsAsFactors = FALSE)
  loading_tbl <- loading_tbl[order(loading_tbl$rank), , drop = FALSE]
  n_total <- nrow(loading_tbl)

  if (!n_total) {
    return(list(
      loading_table = loading_tbl,
      n_total = 0L,
      iod_prefix = numeric(0),
      cv2_prefix = numeric(0),
      valid_prefix = integer(0),
      iod_total = 0,
      cv2_total = 0,
      valid_total = 0L,
      max_leading_frac = max_leading_frac
    ))
  }

  feature_ids <- as.character(loading_tbl$feature_id)
  idx <- match(feature_ids, rownames(count_matrix))
  if (anyNA(idx)) {
    stop(sprintf(
      "Rank prefix cache failed: %d feature_id(s) from loading table were not found in count_matrix.",
      sum(is.na(idx))
    ))
  }

  mat <- as.matrix(count_matrix[idx, , drop = FALSE])
  storage.mode(mat) <- "double"

  feature_mean <- rowMeans(mat, na.rm = TRUE)
  feature_var  <- fast_row_vars(mat)

  valid <- is.finite(feature_mean) & !is.na(feature_mean) & (feature_mean > 0) &
    is.finite(feature_var) & !is.na(feature_var) & (feature_var >= 0)

  iod <- rep(0, n_total)
  cv2 <- rep(0, n_total)
  iod[valid] <- feature_var[valid] / feature_mean[valid]
  cv2[valid] <- feature_var[valid] / (feature_mean[valid]^2)

  list(
    loading_table = loading_tbl,
    n_total = as.integer(n_total),
    iod_prefix = cumsum(iod),
    cv2_prefix = cumsum(cv2),
    valid_prefix = cumsum(as.integer(valid)),
    iod_total = sum(iod),
    cv2_total = sum(cv2),
    valid_total = sum(valid),
    max_leading_frac = max_leading_frac
  )
}

score_single_rank_from_cache <- function(rank_cache, rank_index) {
  n_total <- rank_cache$n_total
  n_leading <- as.integer(rank_index)
  n_remainder <- n_total - n_leading

  if (n_leading < 1L || n_remainder < 1L) {
    return(data.frame(
      rank_index = rank_index,
      objective_score = NA_real_,
      valid_split = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  frac_leading <- n_leading / max(1L, n_total)
  if (!is.finite(frac_leading) || frac_leading > rank_cache$max_leading_frac) {
    return(data.frame(
      rank_index = rank_index,
      objective_score = NA_real_,
      valid_split = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  iod_sum_leading <- rank_cache$iod_prefix[n_leading]
  cv2_sum_leading <- rank_cache$cv2_prefix[n_leading]
  valid_leading <- rank_cache$valid_prefix[n_leading]

  iod_sum_remainder <- rank_cache$iod_total - iod_sum_leading
  cv2_sum_remainder <- rank_cache$cv2_total - cv2_sum_leading
  valid_remainder <- rank_cache$valid_total - valid_leading

  if (valid_leading <= 0L || valid_remainder <= 0L) {
    return(data.frame(
      rank_index = rank_index,
      objective_score = NA_real_,
      valid_split = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  iod_mean_leading <- iod_sum_leading / valid_leading
  iod_mean_remainder <- iod_sum_remainder / valid_remainder
  cv2_mean_leading <- cv2_sum_leading / valid_leading
  cv2_mean_remainder <- cv2_sum_remainder / valid_remainder

  valid_direction <- (iod_mean_leading > iod_mean_remainder) &&
    (cv2_mean_remainder > cv2_mean_leading)

  if (!valid_direction) {
    return(data.frame(
      rank_index = rank_index,
      objective_score = NA_real_,
      valid_split = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  mean_logratio_score <- log(iod_mean_leading / iod_mean_remainder) +
    log(cv2_mean_remainder / cv2_mean_leading)

  size_logratio_score <- log(iod_sum_leading / iod_sum_remainder) +
    log(cv2_sum_remainder / cv2_sum_leading)

  objective_score <- mean_logratio_score + size_logratio_score

  data.frame(
    rank_index = rank_index,
    objective_score = objective_score,
    valid_split = TRUE,
    stringsAsFactors = FALSE
  )
}

evaluate_comparison_shared_rank_candidates <- function(normalized_counts,
                                                       comparison_name,
                                                       fit_trt,
                                                       fit_untrt,
                                                       top_k = evs_top_k_ranks_per_comparison) {
  cache_trt <- build_rank_prefix_cache(normalized_counts, fit_trt)
  cache_untrt <- build_rank_prefix_cache(normalized_counts, fit_untrt)

  shared_n_total <- min(cache_trt$n_total, cache_untrt$n_total)

  if (!is.finite(shared_n_total) || is.na(shared_n_total) || shared_n_total <= 1L) {
    candidate_grid <- data.frame(
      comparison_name = comparison_name,
      rank_index = integer(0),
      local_lagrange_objective = numeric(0),
      valid_in_all_panels = logical(0),
      stringsAsFactors = FALSE
    )

    return(list(
      comparison_summary = data.frame(),
      candidate_grid = candidate_grid,
      best_rank = NA_integer_,
      best_score = NA_real_,
      top_ranks = integer(0),
      top_scores = numeric(0)
    ))
  }

  rank_grid <- seq_len(shared_n_total - 1L)

  score_rows <- lapply(rank_grid, function(r) {
    s1 <- score_single_rank_from_cache(cache_trt, r)
    s2 <- score_single_rank_from_cache(cache_untrt, r)

    both_valid <- isTRUE(s1$valid_split[1]) && isTRUE(s2$valid_split[1])
    score <- if (both_valid) s1$objective_score[1] + s2$objective_score[1] else NA_real_

    data.frame(
      comparison_name = comparison_name,
      rank_index = r,
      local_lagrange_objective = score,
      valid_in_all_panels = both_valid,
      stringsAsFactors = FALSE
    )
  })

  candidate_grid <- bind_rows(score_rows)
  valid_grid <- candidate_grid[candidate_grid$valid_in_all_panels & is.finite(candidate_grid$local_lagrange_objective), , drop = FALSE]
  valid_grid <- valid_grid[order(-valid_grid$local_lagrange_objective, valid_grid$rank_index), , drop = FALSE]

  if (!nrow(valid_grid)) {
    return(list(
      comparison_summary = data.frame(),
      candidate_grid = candidate_grid,
      best_rank = NA_integer_,
      best_score = NA_real_,
      top_ranks = integer(0),
      top_scores = numeric(0)
    ))
  }

  top_k <- min(as.integer(top_k), nrow(valid_grid))
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

apply_loading_cutoff <- function(fit_obj, rank_index, selected_reason = "shared_rank_selected") {
  loading_tbl <- fit_obj$loading_table
  loading_tbl <- loading_tbl[order(loading_tbl$rank), , drop = FALSE]

  rank_index <- max(1L, min(as.integer(rank_index), nrow(loading_tbl)))
  cutoff_value <- as.numeric(loading_tbl$pc1_loading_abs[rank_index])

  loading_tbl$split_class <- ifelse(
    loading_tbl$rank <= rank_index,
    "high_loading",
    "background_loading"
  )

  fit_obj$loading_table <- loading_tbl
  fit_obj$cutoff <- cutoff_value
  fit_obj$top_n_used <- rank_index
  fit_obj$cutoff_quantile <- 1 - (rank_index / nrow(loading_tbl))
  fit_obj$cutoff_method <- "shared_rank_exhaustive_search"
  fit_obj$selected_reason <- selected_reason
  fit_obj
}

precompute_global_evs_rank_selection <- function(comparison_inputs,
                                                 top_k_ranks_per_comparison = evs_top_k_ranks_per_comparison) {
  per_comparison <- list()
  summary_rows <- list()

  for (cmp in names(comparison_inputs)) {
    input_obj <- comparison_inputs[[cmp]]

    dds_init <- DESeqDataSetFromMatrix(
      countData = round(as.matrix(input_obj$count_matrix)),
      colData = input_obj$coldata,
      design = DESIGN_FORMULA
    )
    dds_init <- dds_init[rowSums(counts(dds_init)) > 0, ]
    dds_init <- estimateSizeFactors(dds_init)

    retained_feature_ids <- rownames(dds_init)
    norm_counts_init <- as.data.frame(counts(dds_init, normalized = TRUE))
    raw_counts_init <- as.data.frame(input_obj$count_matrix[retained_feature_ids, , drop = FALSE])

    sample_ids <- colnames(input_obj$count_matrix)
    trt_ids <- sample_ids[input_obj$coldata$condition == "trt"]
    untrt_ids <- sample_ids[input_obj$coldata$condition == "untrt"]

    fit_trt <- compute_pc1_loading_table(norm_counts_init, trt_ids)
    fit_untrt <- compute_pc1_loading_table(norm_counts_init, untrt_ids)
    fit_trt_raw <- compute_pc1_loading_table(raw_counts_init, trt_ids, preprocessing_label = "Raw counts before EVS")
    fit_untrt_raw <- compute_pc1_loading_table(raw_counts_init, untrt_ids, preprocessing_label = "Raw counts before EVS")

    comparison_eval <- evaluate_comparison_shared_rank_candidates(
      normalized_counts = norm_counts_init,
      comparison_name = cmp,
      fit_trt = fit_trt,
      fit_untrt = fit_untrt,
      top_k = top_k_ranks_per_comparison
    )

    per_comparison[[cmp]] <- list(
      comparison_name = cmp,
      retained_feature_ids = retained_feature_ids,
      norm_counts_init = norm_counts_init,
      raw_counts_init = raw_counts_init,
      fit_trt = fit_trt,
      fit_untrt = fit_untrt,
      fit_trt_raw = fit_trt_raw,
      fit_untrt_raw = fit_untrt_raw,
      comparison_eval = comparison_eval
    )

    summary_rows[[cmp]] <- comparison_eval$comparison_summary
  }

  comparison_summary_table <- bind_rows(summary_rows)

  anchor_rank_table <- do.call(rbind, lapply(names(per_comparison), function(cmp) {
    eval_obj <- per_comparison[[cmp]]$comparison_eval
    fit_obj <- per_comparison[[cmp]]$fit_trt
    n_total <- nrow(fit_obj$loading_table)

    if (!length(eval_obj$top_ranks)) return(NULL)

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

  if (is.null(anchor_rank_table)) anchor_rank_table <- data.frame()

  weighted_median_anchor_quantile <- if (nrow(anchor_rank_table)) {
    weighted_median_numeric(anchor_rank_table$anchor_quantile, anchor_rank_table$anchor_score)
  } else {
    NA_real_
  }

  comparison_n_totals <- vapply(per_comparison, function(x) nrow(x$fit_trt$loading_table), integer(1))

  final_shared_rank_by_comparison <- if (is.finite(weighted_median_anchor_quantile) && !is.na(weighted_median_anchor_quantile)) {
    stats::setNames(
      vapply(comparison_n_totals, function(n_total) {
        rank_i <- as.integer(round((1 - weighted_median_anchor_quantile) * n_total))
        max(1L, min(as.integer(n_total) - 1L, rank_i))
      }, integer(1)),
      names(comparison_n_totals)
    )
  } else {
    stats::setNames(rep(NA_integer_, length(comparison_n_totals)), names(comparison_n_totals))
  }

  overall_median_rank <- if (length(final_shared_rank_by_comparison) > 0 &&
                             any(is.finite(final_shared_rank_by_comparison))) {
    as.integer(stats::median(final_shared_rank_by_comparison, na.rm = TRUE))
  } else {
    NA_integer_
  }

  list(
    final_shared_rank = overall_median_rank,
    final_shared_quantile = weighted_median_anchor_quantile,
    final_shared_rank_by_comparison = final_shared_rank_by_comparison,
    comparison_summary_table = comparison_summary_table,
    global_search_table = anchor_rank_table,
    per_comparison = per_comparison
  )
}

# =============================================================================
# PART 3 OF 4
# FINALIZED EVS SPLIT + POST-CUTOFF DESEQ2/HBFSS + FIGURE HELPERS
# =============================================================================

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

  norm_counts_init <- as.data.frame(precomputed_cmp$norm_counts_init)
  fit_trt <- precomputed_cmp$fit_trt
  fit_untrt <- precomputed_cmp$fit_untrt
  fit_trt_raw <- precomputed_cmp$fit_trt_raw
  fit_untrt_raw <- precomputed_cmp$fit_untrt_raw
  comparison_eval <- precomputed_cmp$comparison_eval

  independently_best_rank <- comparison_eval$best_rank
  if (!is.finite(final_shared_rank) || is.na(final_shared_rank)) {
    final_shared_rank <- independently_best_rank
  }

  final_reason <- "comparison_specific_rank_from_score_weighted_median_of_top_k_local_lagrange_anchor_quantiles"

  fit_trt <- apply_loading_cutoff(fit_trt, final_shared_rank, final_reason)
  fit_untrt <- apply_loading_cutoff(fit_untrt, final_shared_rank, final_reason)
  fit_trt_raw <- apply_loading_cutoff(fit_trt_raw, final_shared_rank, paste0(final_reason, "_same_rank_raw_view"))
  fit_untrt_raw <- apply_loading_cutoff(fit_untrt_raw, final_shared_rank, paste0(final_reason, "_same_rank_raw_view"))

  trt_high <- as.character(subset(fit_trt$loading_table, split_class == "high_loading")$feature_id)
  untrt_high <- as.character(subset(fit_untrt$loading_table, split_class == "high_loading")$feature_id)

  leading_edge_ids <- union(trt_high, untrt_high)
  analyzed_feature_ids <- union(as.character(fit_trt$loading_table$feature_id), as.character(fit_untrt$loading_table$feature_id))
  remainder_ids <- setdiff(analyzed_feature_ids, leading_edge_ids)

  if (!length(leading_edge_ids)) {
    stop("Leading-edge dataset is empty after final cutoff.")
  }

  final_rank_aggregation <- data.frame(
    comparison_name = comparison_name,
    independently_best_rank = independently_best_rank,
    independently_best_local_lagrange_objective = comparison_eval$best_score,
    final_shared_rank = final_shared_rank,
    aggregation_method = "comparison_specific_rank_from_score_weighted_median_of_top_k_local_lagrange_anchor_quantiles",
    stringsAsFactors = FALSE
  )

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
    fit_trt = fit_trt,
    fit_untrt = fit_untrt,
    fit_trt_raw = fit_trt_raw,
    fit_untrt_raw = fit_untrt_raw,
    primary_trt_fit = fit_trt,
    primary_untrt_fit = fit_untrt,
    final_rank_aggregation = final_rank_aggregation,
    evs_cutoff_summary = evs_cutoff_summary,
    normalized_counts = norm_counts_init,
    raw_dataset = count_matrix,
    leading_edge_dataset = count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset = count_matrix[remainder_ids, , drop = FALSE]
  )
}

classify_effect_strength <- function(res_strong_padj, res_weak_padj, alpha = alpha_level) {
  out <- rep("intermediate", length(res_strong_padj))
  out[!is.na(res_weak_padj)   & res_weak_padj   < alpha] <- "weak_effect"
  out[!is.na(res_strong_padj) & res_strong_padj < alpha] <- "strong_effect"
  out
}

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
  bp <- get_bpparam()
  use_parallel <- !is.null(bp) && n_workers > 1L

  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_mat)),
    colData   = coldata,
    design    = DESIGN_FORMULA
  )

  dds <- dds[rowSums(counts(dds)) > 0, ]
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
  stat_vec <- as.numeric(res_all_df$stat[valid_stat])
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]

  fdr_fit <- run_empirical_null_fdrtool(stat_vec, dataset_name = dataset_name)

  res_df <- res_all_df
  n_valid <- sum(valid_stat)

  if (length(fdr_fit$pval) != n_valid || length(fdr_fit$qval) != n_valid || length(fdr_fit$lfdr) != n_valid) {
    stop(sprintf("[%s] fdrtool output length mismatch.", dataset_name))
  }

  res_df$empirical_p <- NA_real_
  res_df$empirical_q <- NA_real_
  res_df$lfdr <- NA_real_

  res_df$empirical_p[valid_stat] <- as.numeric(fdr_fit$pval)
  res_df$empirical_q[valid_stat] <- as.numeric(fdr_fit$qval)
  res_df$lfdr[valid_stat] <- as.numeric(fdr_fit$lfdr)

  coef_name <- resultsNames(dds)[grep("^condition_", resultsNames(dds))[1]]

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

  shr <- do.call(lfcShrink, shrink_args)
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

  empirical_p_floored <- ifelse(
    is.na(res_df$empirical_p),
    NA_real_,
    pmax(res_df$empirical_p, 1e-300)
  )

  res_df$HBFSS <- abs(res_df$lfc_shrunk * log10(empirical_p_floored))

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

  res_df$resGA_padj <- res_df$padj_strong_effect
  res_df$resLA_padj <- res_df$padj_weak_effect

  res_df$standard_significant <- !is.na(res_df$resGA_padj) &
    (res_df$resGA_padj < alpha_level) &
    res_df$shrunk_lfc_pass

  res_df$deseq2_strong_call <- res_df$standard_significant
  res_df$deseq2_weak_call <- !is.na(res_df$resLA_padj) &
    (res_df$resLA_padj < alpha_level) &
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

  list(
    dds = dds,
    results = final_df,
    hc_p_threshold = hc_p_threshold_dataset,
    hbfss_threshold = hbfss_threshold_dataset
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
          "The horizontal line marks the loading value at the selected rank."
        ),
        width = 72
      ),
      x = "Ranked PAS feature",
      y = "Absolute PC1 loading"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()
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
# PART 4 OF 4
# VOLCANO/DISPERSION PANELS, EXPORTS, FULL PIPELINE, AND MAIN EXECUTION
# =============================================================================

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
      title   = pretty_dataset_label(dataset_name),
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
      title   = pretty_dataset_label(dataset_name),
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
      title    = compact_title(pretty_dataset_label(dataset_name), width = 42),
      subtitle = NULL,
      x        = "baseMean (log10 scale)",
      y        = "Final dispersion (log10 scale)"
    ) +
    manuscript_theme() +
    volcano_guides()
}

save_cross_dataset_comparison_panels <- function(comparison_name, analysis_results, cmp_dir) {
  if (!isTRUE(export_all_plots)) return(invisible(NULL))

  keys_present <- names(analysis_results)
  if (length(keys_present) == 0) return(invisible(NULL))

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
  save_csv(out, file.path(tab_dir, paste0(comparison_name, "_PC1_loadings_baseMean_dispersion_export.csv")))
  out
}

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df,
                                         final_shared_rank = NA_integer_,
                                         precomputed_selection = NULL) {
  cmp_dir <- file.path(output_dir, comparison_name)
  tab_dir <- file.path(cmp_dir, "tables")
  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

  evs <- build_eigenvector_split(
    count_matrix = count_matrix,
    coldata = coldata,
    comparison_name = comparison_name,
    final_shared_rank = final_shared_rank,
    precomputed_selection = precomputed_selection
  )

  evs_cutoff_summary_out <- evs$evs_cutoff_summary
  evs_cutoff_summary_out$comparison_name <- comparison_name
  save_csv(evs_cutoff_summary_out, file.path(tab_dir, paste0(comparison_name, "_EVS_cutoff_summary.csv")))

  if (!is.null(evs$final_rank_aggregation) && nrow(evs$final_rank_aggregation)) {
    save_csv(
      evs$final_rank_aggregation,
      file.path(tab_dir, paste0(comparison_name, "_EVS_independent_rank_aggregation_summary.csv"))
    )
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

  analysis_results <- list()

  for (nm in names(dataset_list)) {
    full_dataset_name <- paste(comparison_name, nm, sep = "_")

    fit <- tryCatch(
      run_core_analysis(
        count_mat = dataset_list[[nm]],
        coldata = coldata,
        dataset_name = full_dataset_name,
        annot_df = annot_df
      ),
      error = function(e) {
        warning(sprintf("[%s][%s] run_core_analysis failed: %s", comparison_name, nm, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(fit)) next

    df <- fit$results

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

# -----------------------------------------------------------------------------
# MAIN EXECUTION
# -----------------------------------------------------------------------------

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

if (!is.null(global_evs_selection$final_shared_rank_by_comparison) &&
    length(global_evs_selection$final_shared_rank_by_comparison) > 0) {
  save_csv(
    data.frame(
      comparison_name = names(global_evs_selection$final_shared_rank_by_comparison),
      final_shared_rank_by_comparison = as.integer(global_evs_selection$final_shared_rank_by_comparison),
      final_shared_quantile = global_evs_selection$final_shared_quantile,
      overall_median_rank_for_audit = global_evs_selection$final_shared_rank,
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "EVS_final_rank_projection_by_comparison.csv")
  )
}

all_summaries_list <- list()
failed_comparisons <- list()

for (cmp in names(comparison_inputs)) {
  message("\n=====================================================")
  message("Running comparison: ", cmp)
  message("=====================================================")

  input_obj <- comparison_inputs[[cmp]]

  cmp_rank <- NA_integer_
  if (!is.null(global_evs_selection$final_shared_rank_by_comparison) &&
      cmp %in% names(global_evs_selection$final_shared_rank_by_comparison)) {
    cmp_rank <- as.integer(global_evs_selection$final_shared_rank_by_comparison[[cmp]])
  }

  out <- tryCatch(
    run_full_comparison_pipeline(
      comparison_name = input_obj$comparison_name,
      count_matrix = input_obj$count_matrix,
      coldata = input_obj$coldata,
      annot_df = OrigID_Symbol,
      final_shared_rank = cmp_rank,
      precomputed_selection = global_evs_selection
    ),
    error = function(e) {
      failed_comparisons[[cmp]] <<- data.frame(
        comparison_name = cmp,
        error_message = conditionMessage(e),
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
