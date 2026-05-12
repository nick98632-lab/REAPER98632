#!/usr/bin/env Rscript

# =============================================================================
# SEQUENCE FINAL SUBMISSION PIPELINE
# Clean rewrite: DESeq2 + fixed top-N EVS + empirical-null HC/HBFSS
#
# This script intentionally keeps only the essential manuscript workflow:
#   1. import raw WTTS-Seq counts
#   2. run original/raw DESeq2 analysis
#   3. run NormEVS and RawEVS fixed top-N eigenvector splitting
#   4. run DESeq2 + empirical-null HC/HBFSS on Original, Leading Edge, Remainder
#   5. export tables, volcanoes, PCA plots, PC1 variance figures, counts, manifest
#   6. optionally regenerate simulation validation outputs
#
# No legacy changepoint logic.
# No hidden cutoff fallback.
# EVS cutoff is explicitly fixed top-N per treatment and control PCA rank.
# =============================================================================

options(stringsAsFactors = FALSE)

# =============================================================================
# PACKAGES
# =============================================================================

required_packages <- c(
  "DESeq2", "apeglm", "fdrtool", "ggplot2", "ggrepel",
  "dplyr", "gridExtra", "grid", "scales", "S4Vectors", "SummarizedExperiment"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop("Install missing package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

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
  library(S4Vectors)
  library(SummarizedExperiment)
})

# =============================================================================
# SETTINGS
# =============================================================================

count_file_candidates <- c(
  file.path("data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"),
  "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
)

alpha_level <- 0.10
strong_alpha_level <- 0.10
weak_alpha_level <- 0.10
lfc_boundary <- 1.0

evs_cutoff_mode <- "fixed_top_n"
top_n_target <- 5000L
analysis_tracks <- c("NormEVS", "RawEVS")

hc_invalid_at_or_above <- 0.95
calculation_probability_floor <- .Machine$double.xmin
plot_probability_floor <- 1e-16

figure_dpi <- 600
base_theme_size <- 9
n_top_labels <- 6L
n_top_labels_per_class <- 2L

reset_output_dir <- FALSE
run_simulation_validation <- FALSE

simulation_seed <- 42L
simulation_n_features <- 10000L
simulation_n_samples_per_group <- 8L
simulation_base_mean <- 200
simulation_dispersion_null <- 0.10
simulation_dispersion_de <- 0.15
simulation_de_fractions <- c(0.05, 0.10, 0.20)
simulation_lfc_magnitudes <- c(0.25, 0.5, 1.0, 2.0)
simulation_weak_lfc_min <- 0.20
simulation_weak_lfc_max <- 0.80
simulation_null_inflation <- c(0.0, 0.10)
simulation_n_reps <- 20L

class_levels <- c("BG", "Weak", "Strong", "Std", "HBFSS")

class_colors <- c(
  BG = "#BDBDBD",
  Weak = "#0072B2",
  Strong = "#E31A1C",
  Std = "#33A02C",
  HBFSS = "#6A3D9A"
)

class_labels <- c(
  BG = "Background",
  Weak = "Weak effect",
  Strong = "Strong effect",
  Std = "Standard DESeq2",
  HBFSS = "HBFSS-only"
)

class_shapes <- c(BG = 21, Weak = 24, Strong = 22, Std = 23, HBFSS = 25)
class_sizes <- c(BG = 0.55, Weak = 1.10, Strong = 1.15, Std = 1.10, HBFSS = 1.15)
class_alphas <- c(BG = 0.24, Weak = 0.92, Strong = 0.95, Std = 0.90, HBFSS = 0.95)

threshold_color <- "#A65628"
treatment_color <- "#1F78B4"
control_color <- "#4D4D4D"

sample_metadata <- data.frame(
  id = c(
    "R0_1", "R0_2", "R0_3", "R0_4", "R0_5",
    "ZT6_1", "ZT6_2", "ZT6_3", "ZT6_4", "ZT6_5",
    "R2_1", "R2_2", "R2_3", "R2_4", "R2_5",
    "ZT8_1", "ZT8_2", "ZT8_3", "ZT8_4", "ZT8_5",
    "R4_1", "R4_2", "R4_3", "R4_4", "R4_5",
    "ZT10_1", "ZT10_2", "ZT10_3", "ZT10_4", "ZT10_5",
    "R8_1", "R8_2", "R8_3", "R8_4", "R8_5",
    "ZT14_1", "ZT14_2", "ZT14_3", "ZT14_4", "ZT14_5"
  ),
  condition = c(
    rep("trt", 5), rep("untrt", 5),
    rep("trt", 5), rep("untrt", 5),
    rep("trt", 5), rep("untrt", 5),
    rep("trt", 5), rep("untrt", 5)
  ),
  stringsAsFactors = FALSE
)
rownames(sample_metadata) <- sample_metadata$id
sample_metadata$condition <- factor(sample_metadata$condition, levels = c("untrt", "trt"))

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  treatment_prefix = c("R0", "R2", "R4", "R8"),
  control_prefix = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

# =============================================================================
# PATHS AND BASIC HELPERS
# =============================================================================

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0L) return(NA_character_)
  normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = FALSE)
}

find_repo_root <- function() {
  starts <- unique(c(dirname(script_path()), getwd()))
  starts <- starts[!is.na(starts) & dir.exists(starts)]

  walk_up <- function(start) {
    current <- normalizePath(start, winslash = "/", mustWork = TRUE)
    repeat {
      has_git <- dir.exists(file.path(current, ".git"))
      has_data <- file.exists(file.path(current, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"))
      if (has_git || has_data) return(current)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
    NA_character_
  }

  hits <- vapply(starts, walk_up, character(1))
  hits <- hits[!is.na(hits)]
  if (length(hits) > 0L) return(hits[1])
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

repo_root <- find_repo_root()
output_dir <- file.path(repo_root, "exports", "manuscript_final_clean")
figure_dir <- file.path(output_dir, "manuscript_figures")
simulation_dir <- file.path(output_dir, "simulation_validation")

if (isTRUE(reset_output_dir) && dir.exists(output_dir)) {
  unlink(output_dir, recursive = TRUE, force = TRUE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(simulation_dir, recursive = TRUE, showWarnings = FALSE)

cleanup_rplots_pdf <- function() {
  stray <- file.path(repo_root, "Rplots.pdf")
  if (file.exists(stray)) unlink(stray, force = TRUE)
}

resolve_file <- function(candidates, label) {
  candidates <- unique(c(file.path(repo_root, candidates), candidates))
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L) {
    stop("Could not find ", label, ". Tried: ", paste(candidates, collapse = " | "), call. = FALSE)
  }
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

stop_missing_columns <- function(df, required, label) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0L) {
    stop(label, " missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

clean_count_column <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)), fixed = TRUE)))
}

as_integer_count_matrix <- function(x, label) {
  row_ids <- rownames(x)
  col_ids <- colnames(x)
  x <- as.matrix(x)
  suppressWarnings(storage.mode(x) <- "numeric")

  if (any(is.na(x) | !is.finite(x))) stop(label, " has NA or non-finite count values.", call. = FALSE)
  if (any(x < 0)) stop(label, " has negative count values.", call. = FALSE)

  rounded <- round(x)
  if (any(abs(x - rounded) > 1e-6)) {
    warning(label, " had non-integer values; rounded before DESeq2.", call. = FALSE)
  }

  storage.mode(rounded) <- "integer"
  rownames(rounded) <- row_ids
  colnames(rounded) <- col_ids
  rounded
}

clip_probability <- function(x, floor = calculation_probability_floor) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  ok <- !is.na(x)
  x[ok] <- pmin(pmax(x[ok], floor), 1 - 1e-12)
  x
}

format_compact_number <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep("NA", length(x))
  finite <- is.finite(x) & !is.na(x)
  ax <- abs(x[finite])
  val <- x[finite]
  out[finite] <- ifelse(
    ax >= 1e9, paste0(signif(val / 1e9, digits), "B"),
    ifelse(
      ax >= 1e6, paste0(signif(val / 1e6, digits), "M"),
      ifelse(
        ax >= 1e3, paste0(signif(val / 1e3, digits), "K"),
        ifelse(ax > 0 & ax < 0.001, formatC(val, format = "e", digits = digits - 1), as.character(signif(val, digits)))
      )
    )
  )
  out
}

safe_csv_name <- function(...) {
  paste0(gsub("[^A-Za-z0-9_.-]+", "_", paste(..., sep = "_")), ".csv")
}

save_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE)
  invisible(path)
}

analysis_table_path <- function(comparison_name, analysis_label, table_label) {
  path <- file.path(output_dir, comparison_name, analysis_label, "tables")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  file.path(path, safe_csv_name("Table", comparison_name, analysis_label, table_label))
}

rename_columns_existing <- function(df, map) {
  old_names <- names(map)
  hit <- old_names[old_names %in% names(df)]
  if (length(hit) > 0L) {
    names(df)[match(hit, names(df))] <- unname(map[hit])
  }
  df
}

save_plot <- function(plot_obj, path, width, height, export_pdf = TRUE) {
  if (is.null(plot_obj)) return(invisible(NULL))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = figure_dpi,
    bg = "white",
    limitsize = FALSE
  )

  if (isTRUE(export_pdf)) {
    pdf_path <- sub("\\.[^.]+$", ".pdf", path)
    if (!identical(pdf_path, path)) {
      ggplot2::ggsave(
        filename = pdf_path,
        plot = plot_obj,
        width = width,
        height = height,
        units = "in",
        device = "pdf",
        bg = "white",
        limitsize = FALSE
      )
    }
  }

  invisible(path)
}

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_theme_size + 0.5, hjust = 0.5),
      plot.subtitle = element_text(size = base_theme_size - 1.2, hjust = 0.5),
      plot.caption = element_text(size = base_theme_size - 4, color = "grey30", hjust = 0.5),
      axis.title = element_text(face = "bold", size = base_theme_size - 0.1),
      axis.text = element_text(color = "black", size = base_theme_size - 0.8),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_theme_size - 0.2),
      legend.text = element_text(size = base_theme_size - 0.4),
      legend.spacing.x = unit(4, "pt"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey88"),
      strip.text = element_text(face = "bold", size = base_theme_size - 0.2),
      plot.margin = margin(6, 8, 6, 6)
    )
}

get_shared_legend <- function(p) {
  grob <- ggplotGrob(p + theme(legend.position = "bottom"))
  idx <- which(vapply(grob$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(idx) == 0L) return(NULL)
  grob$grobs[[idx[1]]]
}

arrange_with_one_legend <- function(plot_list, title, ncol) {
  plot_list <- Filter(Negate(is.null), plot_list)
  if (length(plot_list) == 0L) return(NULL)

  legend <- get_shared_legend(plot_list[[1]])
  body <- do.call(
    gridExtra::arrangeGrob,
    c(lapply(plot_list, function(p) p + theme(legend.position = "none")), list(ncol = ncol))
  )
  title_grob <- grid::textGrob(title, gp = grid::gpar(fontface = "bold", cex = 1.00))

  if (is.null(legend)) {
    return(gridExtra::arrangeGrob(body, ncol = 1, top = title_grob))
  }

  gridExtra::arrangeGrob(body, legend, ncol = 1, heights = c(12, 0.95), top = title_grob)
}

write_run_session_info <- function() {
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "SessionInfo.txt"))
  invisible(TRUE)
}

fail_pipeline <- function(message_text, status = 1L) {
  try(write_run_session_info(), silent = TRUE)
  try(cleanup_rplots_pdf(), silent = TRUE)
  message(message_text)
  quit(save = "no", status = status, runLast = FALSE)
}

options(error = function() {
  try(write_run_session_info(), silent = TRUE)
  try(cleanup_rplots_pdf(), silent = TRUE)
  traceback(2)
  quit(save = "no", status = 1L, runLast = FALSE)
})

# =============================================================================
# DATA IMPORT
# =============================================================================

read_count_data <- function() {
  count_file <- resolve_file(count_file_candidates, "WTTS count matrix")
  message("Using count file: ", count_file)
  message("Output directory: ", output_dir)

  df <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  stop_missing_columns(df, c("OrigID", "Symbol"), "Count matrix")
  stop_missing_columns(df, sample_metadata$id, "Count matrix")

  df$OrigID <- trimws(as.character(df$OrigID))
  df$Symbol <- trimws(as.character(df$Symbol))

  for (sample_id in sample_metadata$id) {
    df[[sample_id]] <- clean_count_column(df[[sample_id]])
  }

  valid <- !is.na(df$OrigID) & nzchar(df$OrigID)
  valid <- valid & rowSums(is.na(df[, sample_metadata$id, drop = FALSE])) == 0L
  df <- df[valid, , drop = FALSE]

  df$feature_id <- make.unique(df$OrigID, sep = "_dup")
  rownames(df) <- df$feature_id

  annotation <- data.frame(
    feature_id = df$feature_id,
    orig_id = df$OrigID,
    gene_symbol = ifelse(nzchar(df$Symbol), df$Symbol, NA_character_),
    stringsAsFactors = FALSE
  ) %>%
    distinct(feature_id, .keep_all = TRUE)

  list(count_df = df, annotation_df = annotation)
}

input_data <- read_count_data()
count_df <- input_data$count_df
annotation_df <- input_data$annotation_df

prepare_comparison <- function(i) {
  row <- comparison_table[i, , drop = FALSE]

  treatment_ids <- sample_metadata$id[grepl(paste0("^", row$treatment_prefix, "_"), sample_metadata$id)]
  control_ids <- sample_metadata$id[grepl(paste0("^", row$control_prefix, "_"), sample_metadata$id)]
  sample_ids <- c(treatment_ids, control_ids)

  if (length(treatment_ids) < 2L || length(control_ids) < 2L) {
    stop("Comparison ", row$comparison_name, " does not have enough samples.", call. = FALSE)
  }

  coldata <- sample_metadata[sample_ids, "condition", drop = FALSE]
  counts <- count_df[, sample_ids, drop = FALSE]
  rownames(counts) <- count_df$feature_id
  counts <- as_integer_count_matrix(counts, paste0(row$comparison_name, " count matrix"))

  if (!identical(colnames(counts), rownames(coldata))) {
    stop("Sample order mismatch for ", row$comparison_name, call. = FALSE)
  }

  list(
    comparison_name = row$comparison_name,
    count_matrix = counts,
    coldata = coldata
  )
}

# =============================================================================
# EVS
# =============================================================================

make_norm_evs_matrix <- function(count_matrix, coldata) {
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)

  transformed <- tryCatch(
    as.matrix(SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))),
    error = function(e) NULL
  )

  if (is.null(transformed)) {
    transformed <- log2(as.matrix(counts(dds, normalized = TRUE)) + 1)
  }

  storage.mode(transformed) <- "numeric"
  transformed
}

get_evs_matrix <- function(count_matrix, coldata, track) {
  if (identical(track, "RawEVS")) {
    x <- matrix(as.numeric(count_matrix), nrow = nrow(count_matrix), dimnames = dimnames(count_matrix))
    return(log2(x + 1))
  }

  if (identical(track, "NormEVS")) {
    return(make_norm_evs_matrix(count_matrix, coldata))
  }

  stop("Unknown EVS track: ", track, call. = FALSE)
}

condition_pc1_rank <- function(evs_matrix, sample_ids, condition_label) {
  x <- evs_matrix[, sample_ids, drop = FALSE]
  storage.mode(x) <- "numeric"

  keep <- rowSums(is.finite(x) & !is.na(x)) == ncol(x)
  keep <- keep & apply(x, 1, stats::var, na.rm = TRUE) > 0

  if (!any(keep)) stop("No variable features for EVS PCA in ", condition_label, call. = FALSE)

  x <- x[keep, , drop = FALSE]
  pca <- stats::prcomp(t(x), center = TRUE, scale. = FALSE)

  loading <- pca$rotation[, 1]
  pc1_score_variance <- as.numeric(pca$sdev[1]^2)

  rank_df <- data.frame(
    feature_id = names(loading),
    pc1_loading = as.numeric(loading),
    pc1_loading_abs = abs(as.numeric(loading)),
    pc1_score_variance = pc1_score_variance,
    pc1_variance_contribution = (as.numeric(loading)^2) * pc1_score_variance,
    stringsAsFactors = FALSE
  )

  rank_df <- rank_df[order(rank_df$pc1_loading_abs, decreasing = TRUE), , drop = FALSE]
  rank_df$evs_rank <- seq_len(nrow(rank_df))

  top_n_used <- min(top_n_target, nrow(rank_df))
  rank_df$evs_selected <- rank_df$evs_rank <= top_n_used
  rank_df$top_n_used <- top_n_used
  rank_df$loading_cutoff_at_top_n <- rank_df$pc1_loading_abs[top_n_used]
  rank_df$pc1_contribution_cutoff_at_top_n <- rank_df$pc1_variance_contribution[top_n_used]

  list(pca = pca, rank_df = rank_df, top_n_used = top_n_used)
}

build_evs_split <- function(count_matrix, coldata, comparison_name, track) {
  if (!identical(evs_cutoff_mode, "fixed_top_n")) {
    stop("Unsupported evs_cutoff_mode. This clean submission script only implements fixed_top_n.", call. = FALSE)
  }

  evs_matrix <- get_evs_matrix(count_matrix, coldata, track)

  treatment_samples <- rownames(coldata)[coldata$condition == "trt"]
  control_samples <- rownames(coldata)[coldata$condition == "untrt"]
  all_features <- rownames(count_matrix)

  treatment_rank <- condition_pc1_rank(evs_matrix, treatment_samples, "Treatment")
  control_rank <- condition_pc1_rank(evs_matrix, control_samples, "Control")

  treatment_top <- treatment_rank$rank_df$feature_id[treatment_rank$rank_df$evs_selected]
  control_top <- control_rank$rank_df$feature_id[control_rank$rank_df$evs_selected]

  leading_ids <- intersect(union(treatment_top, control_top), all_features)
  remainder_ids <- setdiff(all_features, leading_ids)

  if (length(leading_ids) == 0L) {
    stop("EVS leading edge is empty for ", comparison_name, " ", track, call. = FALSE)
  }

  if (length(remainder_ids) == 0L) {
    stop("EVS remainder is empty for ", comparison_name, " ", track, ". Reduce top_n_target.", call. = FALSE)
  }

  trt_rank <- treatment_rank$rank_df
  ctl_rank <- control_rank$rank_df

  names(trt_rank)[names(trt_rank) != "feature_id"] <- paste0("trt_", names(trt_rank)[names(trt_rank) != "feature_id"])
  names(ctl_rank)[names(ctl_rank) != "feature_id"] <- paste0("ctl_", names(ctl_rank)[names(ctl_rank) != "feature_id"])

  joint_rank <- data.frame(feature_id = all_features, stringsAsFactors = FALSE) %>%
    left_join(trt_rank, by = "feature_id") %>%
    left_join(ctl_rank, by = "feature_id") %>%
    mutate(
      in_leading_edge = feature_id %in% leading_ids,
      in_remainder = feature_id %in% remainder_ids
    ) %>%
    left_join(annotation_df, by = "feature_id")

  summary <- data.frame(
    comparison_name = comparison_name,
    track = track,
    evs_cutoff_mode = evs_cutoff_mode,
    top_n_target = top_n_target,
    treatment_top_n_used = treatment_rank$top_n_used,
    control_top_n_used = control_rank$top_n_used,
    treatment_loading_cutoff = treatment_rank$rank_df$loading_cutoff_at_top_n[1],
    control_loading_cutoff = control_rank$rank_df$loading_cutoff_at_top_n[1],
    leading_edge_n = length(leading_ids),
    remainder_n = length(remainder_ids),
    evs_input = ifelse(
      track == "NormEVS",
      "DESeq2 VST counts for PCA/EVS; log2 normalized-count fallback",
      "log2 raw counts + 1 for PCA/EVS"
    ),
    stringsAsFactors = FALSE
  )

  list(
    comparison_name = comparison_name,
    track = track,
    treatment_rank = treatment_rank,
    control_rank = control_rank,
    joint_rank = joint_rank,
    summary = summary,
    leading_ids = leading_ids,
    remainder_ids = remainder_ids,
    leading_matrix = count_matrix[leading_ids, , drop = FALSE],
    remainder_matrix = count_matrix[remainder_ids, , drop = FALSE]
  )
}

# =============================================================================
# DESEQ2 + EMPIRICAL NULL + HC/HBFSS
# =============================================================================

condition_coef_name <- function(dds) {
  nm <- resultsNames(dds)
  if ("condition_trt_vs_untrt" %in% nm) return("condition_trt_vs_untrt")
  hit <- grep("condition.*trt.*vs.*untrt", nm, value = TRUE)
  if (length(hit) > 0L) return(hit[1])
  stop("Could not find trt-vs-untrt coefficient. Available: ", paste(nm, collapse = ", "), call. = FALSE)
}

extract_eta0 <- function(fit) {
  if (is.null(fit) || is.null(fit$param)) return(NA_real_)
  param <- fit$param

  if (is.matrix(param) || is.data.frame(param)) {
    if ("eta0" %in% colnames(param)) return(suppressWarnings(as.numeric(param[1, "eta0"])))
    if ("eta0" %in% rownames(param)) return(suppressWarnings(as.numeric(param["eta0", 1])))
  }

  if ("eta0" %in% names(param)) return(suppressWarnings(as.numeric(param[["eta0"]][1])))
  NA_real_
}

fit_empirical_null <- function(statistic, label) {
  valid <- is.finite(statistic) & !is.na(statistic)
  if (sum(valid) < 5L) stop(label, " has fewer than five finite Wald statistics.", call. = FALSE)

  z <- statistic[valid]

  run_one <- function(method) {
    fit <- tryCatch(
      fdrtool::fdrtool(
        z,
        statistic = "normal",
        plot = FALSE,
        verbose = FALSE,
        cutoff.method = method,
        pct0 = 0.75
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) return(NULL)
    if (is.null(fit$pval) || length(fit$pval) != length(z)) return(NULL)

    eta0 <- extract_eta0(fit)
    if (!is.finite(eta0) || eta0 <= 0) return(NULL)

    fit
  }

  fit <- run_one("fndr")
  if (is.null(fit)) {
    message(label, ": empirical-null fndr fit failed; retrying pct0.")
    fit <- run_one("pct0")
  }
  if (is.null(fit)) stop(label, ": fdrtool failed with fndr and pct0.", call. = FALSE)

  n <- length(statistic)
  empirical_p <- rep(NA_real_, n)
  empirical_q <- rep(NA_real_, n)
  empirical_lfdr <- rep(NA_real_, n)

  qval <- if (!is.null(fit$qval) && length(fit$qval) == length(z)) fit$qval else rep(NA_real_, length(z))
  lfdr <- if (!is.null(fit$lfdr) && length(fit$lfdr) == length(z)) fit$lfdr else rep(NA_real_, length(z))

  empirical_p[valid] <- clip_probability(fit$pval)
  empirical_q[valid] <- clip_probability(qval)
  empirical_lfdr[valid] <- suppressWarnings(as.numeric(lfdr))

  empirical_bh <- rep(NA_real_, n)
  ok <- !is.na(empirical_p) & is.finite(empirical_p)
  empirical_bh[ok] <- p.adjust(empirical_p[ok], method = "BH")

  list(
    empirical_p = empirical_p,
    empirical_q = empirical_q,
    empirical_lfdr = empirical_lfdr,
    empirical_bh = empirical_bh
  )
}

hc_threshold <- function(empirical_p) {
  p <- sort(clip_probability(empirical_p), decreasing = FALSE, na.last = NA)
  if (length(p) < 5L) return(NA_real_)

  threshold <- tryCatch(
    suppressWarnings(as.numeric(fdrtool::hc.thresh(p)[1])),
    error = function(e) NA_real_
  )

  if (!is.finite(threshold) || is.na(threshold) || threshold <= 0 || threshold >= 1) return(NA_real_)
  if (threshold >= hc_invalid_at_or_above) return(NA_real_)

  threshold
}

classify_results <- function(df, hc_p, hbfss_cutoff) {
  df$standard_flag <- !is.na(df$padj) &
    df$padj < alpha_level &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= lfc_boundary

  df$strong_alt_flag <- !is.na(df$greaterAbs_padj) &
    df$greaterAbs_padj < strong_alpha_level &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= lfc_boundary

  df$strong_flag <- df$strong_alt_flag

  df$lessAbs_alt_flag <- !is.na(df$lessAbs_padj) &
    df$lessAbs_padj < weak_alpha_level &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) < lfc_boundary

  df$hc_pass <- !is.na(hc_p) &
    !is.na(df$empirical_p) &
    is.finite(df$empirical_p) &
    df$empirical_p <= hc_p

  df$HBFSS <- abs(df$lfc_shrunk) * df$neglog10_empirical_p_calc

  df$hbfss_raw_flag <- !is.na(hbfss_cutoff) &
    !is.na(df$HBFSS) &
    is.finite(df$HBFSS) &
    df$HBFSS >= hbfss_cutoff &
    df$hc_pass

  df$hbfss_flag <- df$standard_flag | df$hbfss_raw_flag

  df$weak_region_hbfss_flag <- df$lessAbs_alt_flag & df$hbfss_raw_flag
  df$weak_flag <- df$weak_region_hbfss_flag

  df$display_strong <- df$strong_flag
  df$display_standard <- df$standard_flag & !df$display_strong
  df$display_weak <- df$weak_flag & !df$display_strong & !df$display_standard
  df$display_hbfss <- df$hbfss_flag & !df$standard_flag & !df$display_weak

  df$final_class <- "BG"
  df$final_class[df$display_hbfss] <- "HBFSS"
  df$final_class[df$display_weak] <- "Weak"
  df$final_class[df$display_standard] <- "Std"
  df$final_class[df$display_strong] <- "Strong"
  df$final_class <- factor(df$final_class, levels = class_levels)

  df
}

run_deseq2_hbfss <- function(count_matrix, coldata, comparison_name, analysis_label, dataset_key) {
  label <- paste(comparison_name, analysis_label, sep = "_")

  dds <- DESeqDataSetFromMatrix(
    countData = as_integer_count_matrix(count_matrix, paste0(label, " DESeq2 input")),
    colData = coldata,
    design = ~ condition
  )

  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, betaPrior = FALSE)

  coef_name <- condition_coef_name(dds)

  standard <- results(dds, contrast = c("condition", "trt", "untrt"), alpha = alpha_level)

  strong <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "greaterAbs",
    alpha = strong_alpha_level
  )

  weak <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "lessAbs",
    alpha = weak_alpha_level
  )

  shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")

  df <- as.data.frame(standard)
  df$feature_id <- rownames(df)

  empirical <- fit_empirical_null(df$stat, label)
  df$empirical_p <- empirical$empirical_p
  df$empirical_q <- empirical$empirical_q
  df$empirical_lfdr <- empirical$empirical_lfdr
  df$empirical_bh <- empirical$empirical_bh

  df$empirical_p_calc <- clip_probability(df$empirical_p, floor = calculation_probability_floor)
  df$empirical_p_plot <- clip_probability(df$empirical_p, floor = plot_probability_floor)
  df$neglog10_empirical_p_calc <- -log10(df$empirical_p_calc)
  df$neglog10_empirical_p_plot <- -log10(df$empirical_p_plot)
  df$neglog10_empirical_p <- df$neglog10_empirical_p_plot

  strong_df <- data.frame(
    feature_id = rownames(strong),
    greaterAbs_pvalue = strong$pvalue,
    greaterAbs_padj = strong$padj,
    stringsAsFactors = FALSE
  )

  weak_df <- data.frame(
    feature_id = rownames(weak),
    lessAbs_pvalue = weak$pvalue,
    lessAbs_padj = weak$padj,
    stringsAsFactors = FALSE
  )

  shrink_df <- data.frame(
    feature_id = rownames(shrunk),
    lfc_shrunk = as.data.frame(shrunk)$log2FoldChange,
    stringsAsFactors = FALSE
  )

  df <- df %>%
    left_join(strong_df, by = "feature_id") %>%
    left_join(weak_df, by = "feature_id") %>%
    left_join(shrink_df, by = "feature_id") %>%
    left_join(annotation_df, by = "feature_id")

  hc_p <- hc_threshold(df$empirical_p)
  hbfss_cutoff <- if (is.na(hc_p)) NA_real_ else -log10(hc_p) * lfc_boundary

  df <- classify_results(df, hc_p, hbfss_cutoff)

  normalized_counts <- as.data.frame(counts(dds, normalized = TRUE))
  normalized_counts$feature_id <- rownames(normalized_counts)

  dispersion_df <- as.data.frame(S4Vectors::mcols(dds))
  dispersion_df$feature_id <- rownames(dispersion_df)
  dispersion_keep <- intersect(
    c("feature_id", "dispGeneEst", "dispFit", "dispersion", "dispIter", "dispOutlier"),
    colnames(dispersion_df)
  )
  dispersion_df <- dispersion_df[, dispersion_keep, drop = FALSE]

  df <- df %>%
    left_join(normalized_counts, by = "feature_id") %>%
    left_join(dispersion_df, by = "feature_id")

  df$comparison_name <- comparison_name
  df$analysis_label <- analysis_label
  df$dataset_key <- dataset_key
  df$hc_p_threshold_dataset <- hc_p
  df$hbfss_threshold_dataset <- hbfss_cutoff
  df$regulation_direction <- ifelse(
    is.na(df$lfc_shrunk),
    NA_character_,
    ifelse(df$lfc_shrunk > 0, "upregulated", ifelse(df$lfc_shrunk < 0, "downregulated", "no_change"))
  )

  preferred <- c(
    "comparison_name", "analysis_label", "dataset_key",
    "feature_id", "orig_id", "gene_symbol",
    "baseMean", "log2FoldChange", "lfc_shrunk", "regulation_direction",
    "stat", "pvalue", "padj",
    "greaterAbs_pvalue", "greaterAbs_padj",
    "lessAbs_pvalue", "lessAbs_padj",
    "empirical_p", "empirical_p_calc", "empirical_p_plot",
    "empirical_bh", "empirical_q", "empirical_lfdr",
    "neglog10_empirical_p_calc", "neglog10_empirical_p_plot", "neglog10_empirical_p",
    "HBFSS", "hc_p_threshold_dataset", "hbfss_threshold_dataset",
    "hc_pass", "standard_flag", "strong_alt_flag", "strong_flag",
    "lessAbs_alt_flag", "weak_region_hbfss_flag", "weak_flag",
    "hbfss_raw_flag", "hbfss_flag",
    "display_standard", "display_strong", "display_weak", "display_hbfss",
    "final_class"
  )
  df <- df[, c(intersect(preferred, colnames(df)), setdiff(colnames(df), preferred)), drop = FALSE]

  summary <- data.frame(
    comparison_name = comparison_name,
    analysis_label = analysis_label,
    dataset_key = dataset_key,
    n_features = nrow(df),
    n_lessAbs_alt = sum(df$lessAbs_alt_flag, na.rm = TRUE),
    n_weak = sum(df$weak_region_hbfss_flag, na.rm = TRUE),
    n_strong_alt = sum(df$strong_alt_flag, na.rm = TRUE),
    n_strong = sum(df$strong_flag, na.rm = TRUE),
    n_standard = sum(df$standard_flag, na.rm = TRUE),
    n_hbfss_raw = sum(df$hbfss_raw_flag, na.rm = TRUE),
    n_hbfss_total = sum(df$hbfss_flag, na.rm = TRUE),
    n_hbfss_overlap_standard = sum(df$hbfss_flag & df$standard_flag, na.rm = TRUE),
    n_hbfss_only = sum(df$hbfss_flag & !df$standard_flag, na.rm = TRUE),
    n_hbfss_display_only = sum(df$display_hbfss, na.rm = TRUE),
    hc_p_threshold = hc_p,
    hbfss_threshold = hbfss_cutoff,
    alpha_level = alpha_level,
    strong_alpha_level = strong_alpha_level,
    weak_alpha_level = weak_alpha_level,
    lfc_boundary = lfc_boundary,
    stringsAsFactors = FALSE
  )

  list(dds = dds, results = df, summary = summary)
}

concise_analysis_summary <- function(df) {
  map <- c(
    comparison_name = "comparison",
    analysis_label = "analysis",
    dataset_key = "dataset",
    n_features = "n_features",
    n_lessAbs_alt = "lessAbs_significant",
    n_weak = "weak_effect_hbfss_region",
    n_strong_alt = "greaterAbs_significant",
    n_strong = "strong_effect",
    n_standard = "standard_effect",
    n_hbfss_raw = "hbfss_raw_geometric",
    n_hbfss_total = "hbfss_total",
    n_hbfss_overlap_standard = "hbfss_overlap_standard",
    n_hbfss_only = "hbfss_only",
    n_hbfss_display_only = "hbfss_display_only",
    hc_p_threshold = "hc_p_threshold",
    hbfss_threshold = "hbfss_cutoff",
    alpha_level = "wald_bh_alpha",
    strong_alpha_level = "greaterAbs_bh_alpha",
    weak_alpha_level = "lessAbs_bh_alpha",
    lfc_boundary = "lfc_boundary"
  )

  out <- rename_columns_existing(df, map)

  keep <- c(
    "comparison", "analysis", "dataset", "n_features",
    "lessAbs_significant", "weak_effect_hbfss_region",
    "greaterAbs_significant", "strong_effect",
    "standard_effect", "hbfss_raw_geometric", "hbfss_total",
    "hbfss_overlap_standard", "hbfss_only", "hbfss_display_only",
    "hc_p_threshold", "hbfss_cutoff",
    "wald_bh_alpha", "greaterAbs_bh_alpha", "lessAbs_bh_alpha",
    "lfc_boundary"
  )
  out[, intersect(keep, names(out)), drop = FALSE]
}

concise_evs_summary <- function(df) {
  map <- c(
    comparison_name = "comparison",
    track = "evs_mode",
    evs_cutoff_mode = "evs_cutoff_mode",
    top_n_target = "top_n_target",
    treatment_top_n_used = "treatment_top_n",
    control_top_n_used = "control_top_n",
    treatment_loading_cutoff = "treatment_loading_cutoff",
    control_loading_cutoff = "control_loading_cutoff",
    leading_edge_n = "leading_edge_features",
    remainder_n = "remainder_features",
    evs_input = "input_matrix"
  )

  out <- rename_columns_existing(df, map)

  keep <- c(
    "comparison", "evs_mode", "evs_cutoff_mode", "top_n_target",
    "treatment_top_n", "control_top_n",
    "treatment_loading_cutoff", "control_loading_cutoff",
    "leading_edge_features", "remainder_features", "input_matrix"
  )
  out[, intersect(keep, names(out)), drop = FALSE]
}

# =============================================================================
# FIGURES
# =============================================================================

volcano_data <- function(df) {
  plot_df <- df[
    is.finite(df$lfc_shrunk) & !is.na(df$lfc_shrunk) &
      is.finite(df$neglog10_empirical_p) & !is.na(df$neglog10_empirical_p),
    ,
    drop = FALSE
  ]

  plot_df$gene_label <- ifelse(
    !is.na(plot_df$gene_symbol) & nzchar(trimws(plot_df$gene_symbol)),
    trimws(plot_df$gene_symbol),
    as.character(plot_df$feature_id)
  )

  plot_df$final_class <- factor(as.character(plot_df$final_class), levels = class_levels)

  draw_order <- c(BG = 1, HBFSS = 2, Weak = 3, Std = 4, Strong = 5)
  plot_df$draw_order <- unname(draw_order[as.character(plot_df$final_class)])
  plot_df$draw_order[is.na(plot_df$draw_order)] <- 1
  plot_df[order(plot_df$draw_order, plot_df$neglog10_empirical_p), , drop = FALSE]
}

volcano_labels <- function(plot_df) {
  labels <- plot_df[as.character(plot_df$final_class) != "BG", , drop = FALSE]
  if (nrow(labels) == 0L) return(labels)

  labels <- labels[!duplicated(labels$gene_label), , drop = FALSE]
  class_priority <- c("Strong", "Weak", "Std", "HBFSS")
  picked <- list()

  for (cls in class_priority) {
    sub <- labels[as.character(labels$final_class) == cls, , drop = FALSE]
    if (nrow(sub) == 0L) next
    sub <- sub[order(-sub$HBFSS, sub$empirical_p, -abs(sub$lfc_shrunk), na.last = TRUE), , drop = FALSE]
    picked[[cls]] <- sub[seq_len(min(n_top_labels_per_class, nrow(sub))), , drop = FALSE]
  }

  labels <- if (length(picked) > 0L) do.call(rbind, picked) else labels[0, , drop = FALSE]
  if (nrow(labels) == 0L) return(labels)

  labels <- labels[
    order(match(as.character(labels$final_class), class_priority), -labels$HBFSS, labels$empirical_p, -abs(labels$lfc_shrunk), na.last = TRUE),
    ,
    drop = FALSE
  ]
  labels[seq_len(min(n_top_labels, nrow(labels))), , drop = FALSE]
}

hbfss_boundary <- function(plot_df, hc_p, hbfss_cutoff, y_limit) {
  if (!is.finite(hbfss_cutoff) || is.na(hbfss_cutoff) || hbfss_cutoff <= 0) return(NULL)

  hc_y <- if (!is.na(hc_p) && is.finite(hc_p) && hc_p > 0 && hc_p < 1) -log10(hc_p) else NA_real_
  x_limit <- max(abs(plot_df$lfc_shrunk), lfc_boundary * 1.1, na.rm = TRUE)
  x_start <- max(0.05, hbfss_cutoff / max(y_limit, 1e-6))

  if (!is.finite(x_limit) || x_start >= x_limit) return(NULL)

  x_abs <- seq(x_start, x_limit, length.out = 600)
  y <- hbfss_cutoff / x_abs

  if (!is.na(hc_y) && is.finite(hc_y)) y <- pmax(y, hc_y)
  keep <- is.finite(y) & y >= 0 & y <= y_limit
  if (!any(keep)) return(NULL)

  x_abs <- x_abs[keep]
  y <- y[keep]

  rbind(
    data.frame(x = -rev(x_abs), y = rev(y), stringsAsFactors = FALSE),
    data.frame(x = x_abs, y = y, stringsAsFactors = FALSE)
  )
}

plot_volcano <- function(df, title) {
  plot_df <- volcano_data(df)
  if (nrow(plot_df) == 0L) stop("No finite rows for volcano plot: ", title, call. = FALSE)

  label_df <- volcano_labels(plot_df)
  hc_p <- suppressWarnings(as.numeric(plot_df$hc_p_threshold_dataset[1]))
  hbfss_cutoff <- suppressWarnings(as.numeric(plot_df$hbfss_threshold_dataset[1]))
  hc_y <- if (!is.na(hc_p) && is.finite(hc_p) && hc_p > 0 && hc_p < 1) -log10(hc_p) else NA_real_

  y_limit <- max(plot_df$neglog10_empirical_p, hc_y, na.rm = TRUE) * 1.06
  if (!is.finite(y_limit) || y_limit <= 0) y_limit <- 1

  boundary <- hbfss_boundary(plot_df, hc_p, hbfss_cutoff, y_limit)

  x_min <- min(plot_df$lfc_shrunk, na.rm = TRUE)
  x_max <- max(plot_df$lfc_shrunk, na.rm = TRUE)
  x_span <- max(x_max - x_min, 1e-6)

  threshold_label_y <- if (!is.na(hc_y) && is.finite(hc_y)) {
    max(0.06 * y_limit, hc_y - 0.06 * y_limit)
  } else {
    0.08 * y_limit
  }

  caption <- paste0(
    "Weak=", sum(df$weak_flag, na.rm = TRUE),
    "  Strong=", sum(df$strong_flag, na.rm = TRUE),
    "  Std=", sum(df$standard_flag, na.rm = TRUE),
    "  HBFSS total=", sum(df$hbfss_flag, na.rm = TRUE),
    "  HBFSS-only=", sum(df$display_hbfss, na.rm = TRUE)
  )

  p <- ggplot(
    plot_df,
    aes(
      x = lfc_shrunk,
      y = neglog10_empirical_p,
      color = final_class,
      fill = final_class,
      shape = final_class,
      size = final_class,
      alpha = final_class
    )
  ) +
    geom_point(stroke = 0.32) +
    geom_vline(xintercept = c(-lfc_boundary, lfc_boundary), linetype = "dashed", linewidth = 0.50, color = threshold_color) +
    geom_vline(xintercept = 0, linewidth = 0.30, color = "grey55") +
    scale_color_manual(values = class_colors, breaks = class_levels, labels = class_labels[class_levels], drop = FALSE, name = NULL) +
    scale_fill_manual(values = class_colors, breaks = class_levels, labels = class_labels[class_levels], drop = FALSE, name = NULL) +
    scale_shape_manual(values = class_shapes, breaks = class_levels, labels = class_labels[class_levels], drop = FALSE, name = NULL) +
    scale_size_manual(values = class_sizes, breaks = class_levels, guide = "none") +
    scale_alpha_manual(values = class_alphas, breaks = class_levels, guide = "none") +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = unname(class_shapes[class_levels]),
          color = unname(class_colors[class_levels]),
          fill = unname(class_colors[class_levels]),
          size = rep(3.4, length(class_levels)),
          alpha = rep(1, length(class_levels)),
          stroke = rep(0.80, length(class_levels))
        ),
        nrow = 1,
        byrow = TRUE
      ),
      fill = "none",
      shape = "none",
      size = "none",
      alpha = "none"
    ) +
    labs(
      title = title,
      x = "Shrunken log2 fold change",
      y = expression(-log[10]("empirical p")),
      caption = caption
    ) +
    coord_cartesian(ylim = c(0, y_limit), clip = "off") +
    manuscript_theme()

  if (!is.na(hc_y) && is.finite(hc_y)) {
    p <- p +
      geom_hline(yintercept = hc_y, linetype = "dotted", linewidth = 0.60, color = threshold_color) +
      annotate(
        "text",
        x = x_min + 0.03 * x_span,
        y = threshold_label_y,
        label = paste0("HC p=", signif(hc_p, 3)),
        hjust = 0,
        vjust = 1,
        size = 1.90,
        color = threshold_color
      )
  }

  if (!is.null(boundary)) {
    p <- p + geom_line(data = boundary, aes(x = x, y = y), inherit.aes = FALSE, color = class_colors[["HBFSS"]], linewidth = 0.75)
  }

  if (!is.na(hbfss_cutoff) && is.finite(hbfss_cutoff)) {
    p <- p + annotate(
      "text",
      x = x_min + 0.34 * x_span,
      y = threshold_label_y,
      label = paste0("HBFSS cutoff=", signif(hbfss_cutoff, 3)),
      hjust = 0,
      vjust = 1,
      size = 1.90,
      color = class_colors[["HBFSS"]]
    )
  }

  if (nrow(label_df) > 0L) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      aes(x = lfc_shrunk, y = neglog10_empirical_p, label = gene_label, color = final_class),
      inherit.aes = FALSE,
      show.legend = FALSE,
      size = 1.45,
      seed = 1,
      max.overlaps = Inf,
      force = 0.85,
      force_pull = 0.22,
      box.padding = 0.20,
      point.padding = 0.10,
      min.segment.length = 0,
      segment.alpha = 0.50,
      segment.size = 0.22
    )
  }

  p
}

pca_support_components <- function(fit_obj) {
  dds <- fit_obj$dds
  if (is.null(dds)) return(NULL)

  x <- tryCatch(
    as.matrix(SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))),
    error = function(e) NULL
  )

  if (is.null(x)) {
    x <- log2(as.matrix(DESeq2::counts(dds, normalized = TRUE)) + 1)
  }

  storage.mode(x) <- "numeric"

  keep <- rowSums(is.finite(x) & !is.na(x)) == ncol(x)
  keep <- keep & apply(x, 1, stats::var, na.rm = TRUE) > 0
  x <- x[keep, , drop = FALSE]

  if (nrow(x) < 2L || ncol(x) < 3L) return(NULL)

  pca <- stats::prcomp(t(x), center = TRUE, scale. = FALSE, rank. = 2)
  var_pct <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100, 1)

  pca_df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    condition = as.character(SummarizedExperiment::colData(dds)[rownames(pca$x), "condition"]),
    stringsAsFactors = FALSE
  )

  pca_df$condition <- factor(
    ifelse(pca_df$condition == "trt", "Treatment", "Control"),
    levels = c("Control", "Treatment")
  )

  list(
    pca = pca,
    var_pct = var_pct,
    pca_df = pca_df,
    score_variance = c(PC1 = stats::var(pca_df$PC1, na.rm = TRUE), PC2 = stats::var(pca_df$PC2, na.rm = TRUE))
  )
}

plot_pca_support <- function(fit_obj, title) {
  comp <- pca_support_components(fit_obj)
  if (is.null(comp)) return(NULL)

  pca_df <- comp$pca_df
  var_pct <- comp$var_pct

  ggplot(pca_df, aes(PC1, PC2, label = sample_id, shape = condition, fill = condition)) +
    geom_hline(yintercept = 0, linewidth = 0.22, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.22, linetype = "dashed", color = "grey70") +
    geom_point(size = 2.2, color = "white", stroke = 0.45) +
    ggrepel::geom_text_repel(
      size = 1.9,
      max.overlaps = 12,
      force = 0.7,
      box.padding = 0.18,
      point.padding = 0.10,
      min.segment.length = 0,
      segment.alpha = 0.42,
      segment.size = 0.14
    ) +
    scale_shape_manual(values = c(Control = 21, Treatment = 24), drop = FALSE, name = NULL) +
    scale_fill_manual(values = c(Control = control_color, Treatment = treatment_color), drop = FALSE, name = NULL) +
    labs(
      title = title,
      x = paste0("PC1 ", var_pct[1], "%"),
      y = paste0("PC2 ", var_pct[2], "%"),
      caption = NULL
    ) +
    manuscript_theme()
}

analysis_label_pretty <- function(analysis_label) {
  if (analysis_label == "Raw") return("Original")
  if (grepl("^Lead_", analysis_label)) return(paste("Lead", sub("^Lead_", "", analysis_label)))
  if (grepl("^Rem_", analysis_label)) return(paste("Remainder", sub("^Rem_", "", analysis_label)))
  analysis_label
}

extract_pca_variance_rows <- function(fit_obj, comparison_name, analysis_label) {
  comp <- pca_support_components(fit_obj)
  if (is.null(comp)) return(NULL)

  dataset_group <- if (analysis_label == "Raw") {
    "Original"
  } else if (grepl("^Lead_", analysis_label)) {
    "Lead"
  } else if (grepl("^Rem_", analysis_label)) {
    "Remainder"
  } else {
    analysis_label
  }

  evs_mode <- if (analysis_label == "Raw") "Original" else sub("^[^_]+_", "", analysis_label)

  data.frame(
    comparison_name = comparison_name,
    analysis_label = analysis_label,
    analysis_pretty = analysis_label_pretty(analysis_label),
    dataset_group = factor(dataset_group, levels = c("Original", "Lead", "Remainder")),
    evs_mode = factor(evs_mode, levels = c("Original", "NormEVS", "RawEVS")),
    component = factor(c("PC1", "PC2"), levels = c("PC1", "PC2")),
    variance_explained_pct = as.numeric(comp$var_pct[c(1, 2)]),
    score_variance = as.numeric(comp$score_variance[c("PC1", "PC2")]),
    stringsAsFactors = FALSE
  )
}

pca_pc1_evs_plot_df <- function(pca_var_df) {
  if (nrow(pca_var_df) == 0L) return(data.frame())

  pc1 <- pca_var_df[pca_var_df$component == "PC1", , drop = FALSE]
  original_rows <- pc1[pc1$dataset_group == "Original" & pc1$evs_mode == "Original", , drop = FALSE]
  split_rows <- pc1[pc1$dataset_group %in% c("Lead", "Remainder") & pc1$evs_mode %in% c("NormEVS", "RawEVS"), , drop = FALSE]

  if (nrow(original_rows) == 0L || nrow(split_rows) == 0L) return(data.frame())

  out_rows <- list()

  for (mode_name in c("NormEVS", "RawEVS")) {
    for (comparison_name in unique(split_rows$comparison_name)) {
      original_one <- original_rows[original_rows$comparison_name == comparison_name, , drop = FALSE]
      split_one <- split_rows[split_rows$comparison_name == comparison_name & split_rows$evs_mode == mode_name, , drop = FALSE]
      if (nrow(original_one) == 0L || nrow(split_one) == 0L) next

      original_one <- original_one[1, , drop = FALSE]
      original_one$evs_mode <- mode_name
      original_one$analysis_label <- paste0("Original_for_", mode_name)
      original_one$analysis_pretty <- "Original"

      block <- bind_rows(original_one, split_one)
      original_variance <- suppressWarnings(as.numeric(original_one$score_variance[1]))
      block$original_pc1_score_variance <- original_variance
      block$pc1_score_variance <- suppressWarnings(as.numeric(block$score_variance))
      block$pc1_score_variance_ratio_to_original <- if (is.finite(original_variance) && original_variance > 0) {
        block$pc1_score_variance / original_variance
      } else {
        NA_real_
      }
      block$pc1_score_variance_percent_of_original <- 100 * block$pc1_score_variance_ratio_to_original
      block$dataset_group <- factor(as.character(block$dataset_group), levels = c("Original", "Lead", "Remainder"))
      block$evs_mode <- factor(as.character(block$evs_mode), levels = c("NormEVS", "RawEVS"))

      out_rows[[length(out_rows) + 1L]] <- block
    }
  }

  if (length(out_rows) == 0L) return(data.frame())
  bind_rows(out_rows)
}

plot_pc1_evs_variance_absolute <- function(pca_var_df) {
  plot_df <- pca_pc1_evs_plot_df(pca_var_df)
  if (nrow(plot_df) == 0L) return(NULL)

  plot_df$pc1_score_variance <- pmax(suppressWarnings(as.numeric(plot_df$pc1_score_variance)), 0)

  ggplot(plot_df, aes(dataset_group, pc1_score_variance, fill = dataset_group)) +
    geom_col(width = 0.70, color = "grey25", linewidth = 0.16) +
    geom_text(
      aes(label = format_compact_number(pc1_score_variance, digits = 3)),
      vjust = -0.25,
      size = 1.65,
      color = "grey20"
    ) +
    facet_grid(evs_mode ~ comparison_name, scales = "free_y") +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      labels = function(x) format_compact_number(x, digits = 3),
      expand = expansion(mult = c(0.02, 0.18))
    ) +
    scale_fill_manual(values = c(Original = "#9E9E9E", Lead = "#E31A1C", Remainder = "#1F78B4"), drop = FALSE, name = NULL) +
    labs(
      title = "PC1 score variance by dataset after EVS",
      x = NULL,
      y = "PC1 score variance, absolute scale",
      caption = "Original, EVS leading edge, and EVS remainder are compared directly. PC2 percent variance is not used."
    ) +
    manuscript_theme() +
    theme(
      strip.text = element_text(size = base_theme_size - 0.6),
      axis.text.x = element_text(size = base_theme_size - 0.9),
      legend.position = "bottom"
    )
}

pc1_variance_contribution_rows <- function(expr_matrix, sample_ids, scope_label) {
  x <- expr_matrix[, sample_ids, drop = FALSE]
  storage.mode(x) <- "numeric"

  keep <- rowSums(is.finite(x) & !is.na(x)) == ncol(x)
  keep <- keep & apply(x, 1, stats::var, na.rm = TRUE) > 0
  x <- x[keep, , drop = FALSE]

  if (nrow(x) < 2L || ncol(x) < 3L) return(data.frame())

  pca <- stats::prcomp(t(x), center = TRUE, scale. = FALSE, rank. = 1)
  loading <- as.numeric(pca$rotation[, 1])
  pc1_score_variance <- as.numeric(pca$sdev[1]^2)

  data.frame(
    feature_id = rownames(pca$rotation),
    condition = scope_label,
    pc1_loading = loading,
    pc1_loading_abs = abs(loading),
    pc1_score_variance = pc1_score_variance,
    pc1_variance_contribution = (loading^2) * pc1_score_variance,
    stringsAsFactors = FALSE
  )
}

evs_loading_distribution_rows <- function(full_count_matrix, evs_obj, coldata, comparison_name, track) {
  extract_one <- function(mat, dataset_group) {
    source_matrix <- get_evs_matrix(mat, coldata, track)
    trt_ids <- rownames(coldata)[coldata$condition == "trt"]
    ctl_ids <- rownames(coldata)[coldata$condition == "untrt"]
    all_ids <- rownames(coldata)

    bind_rows(
      pc1_variance_contribution_rows(source_matrix, all_ids, "All"),
      pc1_variance_contribution_rows(source_matrix, ctl_ids, "Control"),
      pc1_variance_contribution_rows(source_matrix, trt_ids, "Treatment")
    ) %>%
      mutate(
        comparison_name = comparison_name,
        track = track,
        dataset_group = dataset_group
      )
  }

  out <- bind_rows(
    extract_one(full_count_matrix, "Original"),
    extract_one(evs_obj$leading_matrix, "Lead"),
    extract_one(evs_obj$remainder_matrix, "Remainder")
  )

  cutoff_rows <- data.frame(
    comparison_name = comparison_name,
    track = track,
    condition = c("Treatment", "Control"),
    pc1_loading_cutoff = c(evs_obj$summary$treatment_loading_cutoff[1], evs_obj$summary$control_loading_cutoff[1]),
    pc1_score_variance = c(
      as.numeric(evs_obj$treatment_rank$pca$sdev[1]^2),
      as.numeric(evs_obj$control_rank$pca$sdev[1]^2)
    ),
    stringsAsFactors = FALSE
  )

  cutoff_rows$pc1_variance_contribution_cutoff <- (cutoff_rows$pc1_loading_cutoff^2) * cutoff_rows$pc1_score_variance
  cutoff_rows <- cutoff_rows[rep(seq_len(nrow(cutoff_rows)), each = 3), , drop = FALSE]
  cutoff_rows$dataset_group <- rep(c("Original", "Lead", "Remainder"), times = 2L)

  attr(out, "cutoffs") <- cutoff_rows
  out
}

plot_loading_histograms <- function(load_df, title_text) {
  if (nrow(load_df) == 0L) return(NULL)

  cutoff_df <- attr(load_df, "cutoffs")
  if (is.null(cutoff_df)) cutoff_df <- data.frame()

  load_df$dataset_group <- factor(load_df$dataset_group, levels = c("Original", "Lead", "Remainder"))
  load_df$condition <- factor(load_df$condition, levels = c("All", "Control", "Treatment"))
  load_df$pc1_variance_contribution <- pmax(suppressWarnings(as.numeric(load_df$pc1_variance_contribution)), 0)

  if (nrow(cutoff_df) > 0L) {
    cutoff_df$dataset_group <- factor(cutoff_df$dataset_group, levels = c("Original", "Lead", "Remainder"))
    cutoff_df$condition <- factor(cutoff_df$condition, levels = c("All", "Control", "Treatment"))
  }

  p <- ggplot(load_df, aes(pc1_variance_contribution, color = condition, fill = condition)) +
    geom_histogram(
      bins = 60,
      position = "identity",
      alpha = 0.32,
      linewidth = 0.22,
      boundary = 0
    ) +
    facet_grid(dataset_group ~ comparison_name, scales = "free_y") +
    scale_x_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      labels = function(x) format_compact_number(x, digits = 3)
    ) +
    scale_color_manual(values = c(All = "#7570B3", Control = control_color, Treatment = treatment_color), breaks = c("All", "Control", "Treatment"), drop = FALSE, name = NULL) +
    scale_fill_manual(values = c(All = "#7570B3", Control = control_color, Treatment = treatment_color), breaks = c("All", "Control", "Treatment"), drop = FALSE, name = NULL) +
    guides(color = guide_legend(override.aes = list(fill = c("#7570B3", control_color, treatment_color), alpha = 0.85, linewidth = 0.8)), fill = "none") +
    labs(
      title = title_text,
      x = "Per-feature contribution to PC1 score variance",
      y = "Feature frequency",
      caption = "Dashed vertical lines mark fixed top-N treatment/control EVS loading cutoffs projected onto PC1 variance contribution."
    ) +
    manuscript_theme() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = base_theme_size - 0.6),
      axis.text = element_text(size = base_theme_size - 1.0)
    )

  if (nrow(cutoff_df) > 0L) {
    p <- p + geom_vline(
      data = cutoff_df,
      aes(xintercept = pc1_variance_contribution_cutoff, color = condition),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.38,
      alpha = 0.85
    )
  }

  p
}

plot_discovery_counts <- function(summary_df) {
  long_df <- bind_rows(lapply(seq_len(nrow(summary_df)), function(i) {
    row <- summary_df[i, , drop = FALSE]
    data.frame(
      comparison_name = row$comparison_name,
      analysis_label = row$analysis_label,
      dataset_key = row$dataset_key,
      Class = factor(c("Weak", "Strong", "Std", "HBFSS"), levels = c("Weak", "Strong", "Std", "HBFSS")),
      Count = as.numeric(c(row$n_weak, row$n_strong, row$n_standard, row$n_hbfss_display_only)),
      stringsAsFactors = FALSE
    )
  }))

  save_csv(
    rename_columns_existing(
      long_df,
      c(comparison_name = "comparison", analysis_label = "analysis", dataset_key = "dataset", Class = "class", Count = "count")
    ),
    file.path(output_dir, "Counts_Long.csv")
  )

  ggplot(long_df, aes(comparison_name, Count, fill = Class)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.74, color = "grey25", linewidth = 0.15) +
    facet_wrap(~ analysis_label, scales = "free_y", ncol = 3) +
    scale_fill_manual(
      values = class_colors[c("Weak", "Strong", "Std", "HBFSS")],
      breaks = c("Weak", "Strong", "Std", "HBFSS"),
      labels = class_labels[c("Weak", "Strong", "Std", "HBFSS")],
      drop = FALSE,
      name = NULL
    ) +
    labs(title = "Significant feature counts", x = NULL, y = "Count", caption = NULL) +
    manuscript_theme() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

# =============================================================================
# SIMULATION VALIDATION
# =============================================================================

# Simulation design:
# - Truth is known at the feature level.
# - Counts are generated from a negative-binomial model with heterogeneous means,
#   heterogeneous dispersions, sample-specific library factors, bidirectional DE,
#   and optional null-statistic inflation.
# - Every successful replicate exports method-level confusion matrices.
# - Every failed replicate exports the parameter combination and error message.
# - HBFSS is never silently rescued with a hidden cutoff. If HC is invalid,
#   HBFSS raw calls are zero for that replicate and the invalid threshold is recorded.

export_simulation_feature_results <- FALSE

simulation_lfc_label <- function(lfc_magnitude, simulation_profile) {
  profile <- as.character(simulation_profile)
  lfc <- suppressWarnings(as.numeric(lfc_magnitude))
  ifelse(
    profile == "weak_mixture",
    paste0("Weak mix\nLFC ", simulation_weak_lfc_min, "-", simulation_weak_lfc_max),
    paste0("Fixed LFC\n", lfc)
  )
}

simulation_lfc_levels <- function() {
  c(
    paste0("Fixed LFC\n", simulation_lfc_magnitudes),
    paste0("Weak mix\nLFC ", simulation_weak_lfc_min, "-", simulation_weak_lfc_max)
  )
}

simulate_sequence_counts <- function(n_features, n_samples, base_mean, disp_null,
                                     de_fraction, lfc_magnitude, disp_de,
                                     lfc_profile = c("fixed", "weak_mixture")) {
  lfc_profile <- match.arg(lfc_profile)

  n_de <- round(n_features * de_fraction)
  n_de <- max(2L, min(n_de, n_features - 2L))
  n_null <- n_features - n_de

  feature_id <- paste0("sim_feature_", seq_len(n_features))
  de_id <- seq_len(n_de)
  null_id <- seq.int(n_de + 1L, n_features)

  # Heterogeneous baseline expression prevents the simulation from being an
  # unrealistically single-mean test. The median remains near simulation_base_mean.
  base_mu <- stats::rgamma(
    n_features,
    shape = 2.5,
    scale = base_mean / 2.5
  )
  base_mu <- pmax(base_mu, 2)

  # Feature-level dispersion heterogeneity around null/DE dispersion targets.
  dispersion <- rep(disp_null, n_features)
  dispersion[de_id] <- disp_de
  dispersion <- dispersion * exp(stats::rnorm(n_features, mean = 0, sd = 0.25))
  dispersion <- pmin(pmax(dispersion, 0.01), 1.50)

  de_direction <- sample(c(-1, 1), n_de, replace = TRUE)

  if (identical(lfc_profile, "weak_mixture")) {
    lfc_abs <- stats::runif(n_de, min = simulation_weak_lfc_min, max = simulation_weak_lfc_max)
  } else {
    if (!is.finite(lfc_magnitude)) {
      stop("Fixed-LFC simulation requires finite lfc_magnitude.", call. = FALSE)
    }
    lfc_abs <- rep(lfc_magnitude, n_de)
  }

  true_lfc <- rep(0, n_features)
  true_lfc[de_id] <- lfc_abs * de_direction

  control_mu <- base_mu
  treatment_mu <- base_mu * 2^true_lfc

  # Mild sample-level library factors create realistic normalization work.
  control_lib <- exp(stats::rnorm(n_samples, mean = 0, sd = 0.12))
  treatment_lib <- exp(stats::rnorm(n_samples, mean = 0, sd = 0.12))
  control_lib <- control_lib / exp(mean(log(control_lib)))
  treatment_lib <- treatment_lib / exp(mean(log(treatment_lib)))

  generate_group <- function(mu, lib_factor, dispersion_vector) {
    out <- matrix(0L, nrow = length(mu), ncol = length(lib_factor))
    for (j in seq_along(lib_factor)) {
      sample_mu <- pmax(mu * lib_factor[j], 1e-3)
      out[, j] <- stats::rnbinom(
        n = length(sample_mu),
        mu = sample_mu,
        size = 1 / dispersion_vector
      )
    }
    out
  }

  counts <- cbind(
    generate_group(control_mu, control_lib, dispersion),
    generate_group(treatment_mu, treatment_lib, dispersion)
  )
  storage.mode(counts) <- "integer"
  rownames(counts) <- feature_id
  colnames(counts) <- c(paste0("ctrl_", seq_len(n_samples)), paste0("trt_", seq_len(n_samples)))

  truth <- data.frame(
    feature_id = feature_id,
    is_de = seq_len(n_features) %in% de_id,
    is_null = seq_len(n_features) %in% null_id,
    true_lfc = true_lfc,
    true_abs_lfc = abs(true_lfc),
    true_direction = ifelse(true_lfc > 0, "up", ifelse(true_lfc < 0, "down", "null")),
    true_weak = seq_len(n_features) %in% de_id & abs(true_lfc) < lfc_boundary,
    true_strong = seq_len(n_features) %in% de_id & abs(true_lfc) >= lfc_boundary,
    base_mean = base_mu,
    dispersion = dispersion,
    simulation_profile = lfc_profile,
    stringsAsFactors = FALSE
  )

  list(counts = counts, truth = truth)
}

inflate_null_wald_statistics <- function(wald, is_de, inflation_fraction, inflation_sd = 1.75) {
  if (!is.finite(inflation_fraction) || inflation_fraction <= 0) return(wald)

  out <- wald
  null_idx <- which(!is_de & is.finite(out) & !is.na(out))
  n_inflate <- round(length(null_idx) * inflation_fraction)
  if (n_inflate <= 0L) return(out)

  target <- sample(null_idx, n_inflate, replace = FALSE)
  out[target] <- stats::rnorm(n_inflate, mean = 0, sd = inflation_sd)
  out
}

simulation_metrics <- function(predicted, actual) {
  predicted <- !is.na(predicted) & predicted
  actual <- !is.na(actual) & actual

  tp <- sum(predicted & actual)
  fp <- sum(predicted & !actual)
  fn <- sum(!predicted & actual)
  tn <- sum(!predicted & !actual)

  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  specificity <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  fdr <- if ((tp + fp) == 0) NA_real_ else fp / (tp + fp)
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA_real_ else {
    2 * precision * recall / (precision + recall)
  }

  data.frame(
    tp = tp,
    fp = fp,
    fn = fn,
    tn = tn,
    precision = precision,
    recall = recall,
    specificity = specificity,
    f1 = f1,
    fdr = fdr,
    discovery_count = sum(predicted),
    truth_count = sum(actual),
    stringsAsFactors = FALSE
  )
}

simulation_metric_row <- function(template, method, truth_target, predicted, actual,
                                  hc_p, hbfss_cutoff, diagnostics) {
  met <- simulation_metrics(predicted, actual)

  cbind(
    template,
    data.frame(
      method = method,
      truth_target = truth_target,
      hc_p_threshold = hc_p,
      hbfss_cutoff = hbfss_cutoff,
      hc_valid = is.finite(hc_p) & !is.na(hc_p),
      n_hc_pass = diagnostics$n_hc_pass,
      n_standard = diagnostics$n_standard,
      n_empirical_bh = diagnostics$n_empirical_bh,
      n_lessAbs = diagnostics$n_lessAbs,
      n_greaterAbs = diagnostics$n_greaterAbs,
      n_hbfss_raw = diagnostics$n_hbfss_raw,
      n_hbfss_total = diagnostics$n_hbfss_total,
      n_weak_hbfss = diagnostics$n_weak_hbfss,
      stringsAsFactors = FALSE
    ),
    met
  )
}

sequence_simulation_analysis <- function(sim_obj, null_inflation) {
  counts <- sim_obj$counts
  n_samp <- ncol(counts) / 2L

  coldata <- data.frame(
    condition = factor(c(rep("untrt", n_samp), rep("trt", n_samp)), levels = c("untrt", "trt")),
    row.names = colnames(counts)
  )

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, betaPrior = FALSE, quiet = TRUE)

  coef_name <- condition_coef_name(dds)

  standard <- results(dds, contrast = c("condition", "trt", "untrt"), alpha = alpha_level)

  greater_abs <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "greaterAbs",
    alpha = strong_alpha_level
  )

  less_abs <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "lessAbs",
    alpha = weak_alpha_level
  )

  shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)

  df <- as.data.frame(standard)
  df$feature_id <- rownames(df)

  df <- df %>%
    left_join(
      data.frame(
        feature_id = rownames(shrunk),
        lfc_shrunk = as.data.frame(shrunk)$log2FoldChange,
        stringsAsFactors = FALSE
      ),
      by = "feature_id"
    ) %>%
    left_join(
      data.frame(
        feature_id = rownames(greater_abs),
        greaterAbs_pvalue = greater_abs$pvalue,
        greaterAbs_padj = greater_abs$padj,
        stringsAsFactors = FALSE
      ),
      by = "feature_id"
    ) %>%
    left_join(
      data.frame(
        feature_id = rownames(less_abs),
        lessAbs_pvalue = less_abs$pvalue,
        lessAbs_padj = less_abs$padj,
        stringsAsFactors = FALSE
      ),
      by = "feature_id"
    ) %>%
    left_join(sim_obj$truth, by = "feature_id")

  wald_for_empirical_null <- inflate_null_wald_statistics(
    wald = df$stat,
    is_de = df$is_de,
    inflation_fraction = null_inflation
  )

  empirical <- fit_empirical_null(wald_for_empirical_null, "simulation")
  df$wald_for_empirical_null <- wald_for_empirical_null
  df$empirical_p <- empirical$empirical_p
  df$empirical_bh <- empirical$empirical_bh
  df$empirical_q <- empirical$empirical_q
  df$empirical_lfdr <- empirical$empirical_lfdr
  df$empirical_p_calc <- clip_probability(df$empirical_p, floor = calculation_probability_floor)
  df$neglog10_empirical_p_calc <- -log10(df$empirical_p_calc)

  hc_p <- hc_threshold(df$empirical_p)
  hbfss_cutoff <- if (is.na(hc_p)) NA_real_ else -log10(hc_p) * lfc_boundary

  df$hc_pass <- !is.na(hc_p) &
    !is.na(df$empirical_p) &
    is.finite(df$empirical_p) &
    df$empirical_p <= hc_p

  df$HBFSS <- abs(df$lfc_shrunk) * df$neglog10_empirical_p_calc

  df$standard_sig <- !is.na(df$padj) &
    df$padj < alpha_level &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= lfc_boundary

  df$greaterAbs_sig <- !is.na(df$greaterAbs_padj) &
    df$greaterAbs_padj < strong_alpha_level &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= lfc_boundary

  df$lessAbs_sig <- !is.na(df$lessAbs_padj) &
    df$lessAbs_padj < weak_alpha_level &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) < lfc_boundary

  df$empirical_bh_sig <- !is.na(df$empirical_bh) &
    df$empirical_bh < alpha_level &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= lfc_boundary

  df$hbfss_raw_sig <- !is.na(hbfss_cutoff) &
    !is.na(df$HBFSS) &
    is.finite(df$HBFSS) &
    df$HBFSS >= hbfss_cutoff &
    df$hc_pass

  df$hbfss_total_sig <- df$standard_sig | df$hbfss_raw_sig
  df$weak_region_hbfss_sig <- df$lessAbs_sig & df$hbfss_raw_sig

  diagnostics <- list(
    n_hc_pass = sum(df$hc_pass, na.rm = TRUE),
    n_standard = sum(df$standard_sig, na.rm = TRUE),
    n_empirical_bh = sum(df$empirical_bh_sig, na.rm = TRUE),
    n_lessAbs = sum(df$lessAbs_sig, na.rm = TRUE),
    n_greaterAbs = sum(df$greaterAbs_sig, na.rm = TRUE),
    n_hbfss_raw = sum(df$hbfss_raw_sig, na.rm = TRUE),
    n_hbfss_total = sum(df$hbfss_total_sig, na.rm = TRUE),
    n_weak_hbfss = sum(df$weak_region_hbfss_sig, na.rm = TRUE)
  )

  list(
    results = df,
    hc_p = hc_p,
    hbfss_cutoff = hbfss_cutoff,
    diagnostics = diagnostics
  )
}

build_simulation_metric_rows <- function(out, template) {
  df <- out$results

  all_de <- df$is_de
  weak_de <- df$true_weak
  strong_de <- df$true_strong
  null_truth <- df$is_null

  bind_rows(
    simulation_metric_row(template, "DESeq2_BH", "all_de", df$standard_sig, all_de, out$hc_p, out$hbfss_cutoff, out$diagnostics),
    simulation_metric_row(template, "Empirical_BH", "all_de", df$empirical_bh_sig, all_de, out$hc_p, out$hbfss_cutoff, out$diagnostics),
    simulation_metric_row(template, "GreaterAbs", "strong_de", df$greaterAbs_sig, strong_de, out$hc_p, out$hbfss_cutoff, out$diagnostics),
    simulation_metric_row(template, "HBFSS_total", "all_de", df$hbfss_total_sig, all_de, out$hc_p, out$hbfss_cutoff, out$diagnostics),
    simulation_metric_row(template, "HBFSS_raw", "all_de", df$hbfss_raw_sig, all_de, out$hc_p, out$hbfss_cutoff, out$diagnostics),
    simulation_metric_row(template, "DESeq2_BH", "weak_de", df$standard_sig, weak_de, out$hc_p, out$hbfss_cutoff, out$diagnostics),
    simulation_metric_row(template, "LessAbs", "weak_de", df$lessAbs_sig, weak_de, out$hc_p, out$hbfss_cutoff, out$diagnostics),
    simulation_metric_row(template, "HBFSS_weak_region", "weak_de", df$weak_region_hbfss_sig, weak_de, out$hc_p, out$hbfss_cutoff, out$diagnostics),
    simulation_metric_row(template, "HBFSS_raw", "weak_de", df$hbfss_raw_sig, weak_de, out$hc_p, out$hbfss_cutoff, out$diagnostics),
    simulation_metric_row(template, "DESeq2_BH", "null_false_positive", df$standard_sig, null_truth, out$hc_p, out$hbfss_cutoff, out$diagnostics),
    simulation_metric_row(template, "HBFSS_raw", "null_false_positive", df$hbfss_raw_sig, null_truth, out$hc_p, out$hbfss_cutoff, out$diagnostics)
  )
}

simulation_wide_from_long <- function(metric_long) {
  if (nrow(metric_long) == 0L) return(data.frame())

  key_cols <- c(
    "de_fraction", "lfc_magnitude", "simulation_profile", "null_inflation",
    "replicate", "lfc_label", "de_label", "inflation_label",
    "hc_p_threshold", "hbfss_cutoff", "hc_valid",
    "n_hc_pass", "n_standard", "n_empirical_bh", "n_lessAbs",
    "n_greaterAbs", "n_hbfss_raw", "n_hbfss_total", "n_weak_hbfss"
  )

  base <- metric_long[!duplicated(metric_long[, key_cols]), key_cols, drop = FALSE]

  make_metric <- function(method_name, target_name, metric_name, output_name) {
    sub <- metric_long[
      metric_long$method == method_name & metric_long$truth_target == target_name,
      c(key_cols, metric_name),
      drop = FALSE
    ]
    names(sub)[names(sub) == metric_name] <- output_name
    sub
  }

  pieces <- list(
    make_metric("DESeq2_BH", "all_de", "precision", "standard_precision"),
    make_metric("DESeq2_BH", "all_de", "recall", "standard_recall"),
    make_metric("DESeq2_BH", "all_de", "f1", "standard_f1"),
    make_metric("DESeq2_BH", "all_de", "fdr", "standard_fdr"),
    make_metric("Empirical_BH", "all_de", "precision", "empirical_bh_precision"),
    make_metric("Empirical_BH", "all_de", "recall", "empirical_bh_recall"),
    make_metric("Empirical_BH", "all_de", "f1", "empirical_bh_f1"),
    make_metric("Empirical_BH", "all_de", "fdr", "empirical_bh_fdr"),
    make_metric("HBFSS_total", "all_de", "precision", "hbfss_precision"),
    make_metric("HBFSS_total", "all_de", "recall", "hbfss_recall"),
    make_metric("HBFSS_total", "all_de", "f1", "hbfss_f1"),
    make_metric("HBFSS_total", "all_de", "fdr", "hbfss_fdr"),
    make_metric("DESeq2_BH", "weak_de", "recall", "standard_weak_recall"),
    make_metric("LessAbs", "weak_de", "recall", "lessAbs_weak_recall"),
    make_metric("HBFSS_weak_region", "weak_de", "recall", "hbfss_weak_recall")
  )

  out <- base
  for (piece in pieces) {
    out <- left_join(out, piece, by = key_cols)
  }
  out
}

simulation_method_summary <- function(metric_long) {
  metric_long %>%
    group_by(
      de_fraction, lfc_magnitude, simulation_profile, null_inflation,
      lfc_label, de_label, inflation_label, method, truth_target
    ) %>%
    summarise(
      n_replicates = n(),
      precision_mean = mean(precision, na.rm = TRUE),
      precision_median = median(precision, na.rm = TRUE),
      recall_mean = mean(recall, na.rm = TRUE),
      recall_median = median(recall, na.rm = TRUE),
      f1_mean = mean(f1, na.rm = TRUE),
      f1_median = median(f1, na.rm = TRUE),
      fdr_mean = mean(fdr, na.rm = TRUE),
      fdr_median = median(fdr, na.rm = TRUE),
      discovery_count_mean = mean(discovery_count, na.rm = TRUE),
      truth_count_mean = mean(truth_count, na.rm = TRUE),
      hc_valid_fraction = mean(hc_valid, na.rm = TRUE),
      hbfss_cutoff_median = median(hbfss_cutoff, na.rm = TRUE),
      .groups = "drop"
    )
}

finite_summary <- function(x, fun) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  fun(x)
}

simulation_threshold_summary <- function(metric_wide) {
  if (nrow(metric_wide) == 0L) return(data.frame())

  metric_wide %>%
    group_by(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label) %>%
    summarise(
      n_replicates = n(),
      hc_valid_fraction = mean(hc_valid, na.rm = TRUE),
      hc_p_threshold_median = median(hc_p_threshold, na.rm = TRUE),
      hc_p_threshold_min = finite_summary(hc_p_threshold, min),
      hc_p_threshold_max = finite_summary(hc_p_threshold, max),
      hbfss_cutoff_median = median(hbfss_cutoff, na.rm = TRUE),
      hbfss_cutoff_min = finite_summary(hbfss_cutoff, min),
      hbfss_cutoff_max = finite_summary(hbfss_cutoff, max),
      n_hc_pass_mean = mean(n_hc_pass, na.rm = TRUE),
      n_standard_mean = mean(n_standard, na.rm = TRUE),
      n_empirical_bh_mean = mean(n_empirical_bh, na.rm = TRUE),
      n_lessAbs_mean = mean(n_lessAbs, na.rm = TRUE),
      n_greaterAbs_mean = mean(n_greaterAbs, na.rm = TRUE),
      n_hbfss_raw_mean = mean(n_hbfss_raw, na.rm = TRUE),
      n_hbfss_total_mean = mean(n_hbfss_total, na.rm = TRUE),
      n_weak_hbfss_mean = mean(n_weak_hbfss, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_simulation_metric_boxplot <- function(metric_long, metric, title, y_label,
                                           methods, target, alpha_line = FALSE) {
  plot_df <- metric_long[
    metric_long$method %in% methods & metric_long$truth_target == target,
    ,
    drop = FALSE
  ]
  if (nrow(plot_df) == 0L) return(NULL)

  plot_df$method <- factor(plot_df$method, levels = methods)

  method_colors <- c(
    DESeq2_BH = "#999999",
    Empirical_BH = treatment_color,
    GreaterAbs = class_colors[["Strong"]],
    LessAbs = class_colors[["Weak"]],
    HBFSS_total = class_colors[["HBFSS"]],
    HBFSS_raw = "#54278F",
    HBFSS_weak_region = class_colors[["Weak"]]
  )

  p <- ggplot(plot_df, aes(method, .data[[metric]], fill = method)) +
    geom_boxplot(outlier.size = 0.35, width = 0.62, linewidth = 0.22, na.rm = TRUE) +
    facet_grid(inflation_label + de_label ~ lfc_label) +
    scale_fill_manual(values = method_colors[methods], breaks = methods, drop = FALSE, name = NULL) +
    labs(title = title, x = NULL, y = y_label, caption = NULL) +
    manuscript_theme() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(size = base_theme_size - 1.0),
      legend.position = "bottom"
    )

  if (alpha_line) {
    p <- p + geom_hline(yintercept = alpha_level, linetype = "dashed", linewidth = 0.30, color = "grey35")
  }

  p
}

plot_simulation_delta_heatmap <- function(metric_long, metric, title, fill_label,
                                          hbfss_method = "HBFSS_total",
                                          baseline_method = "DESeq2_BH",
                                          target = "all_de") {
  plot_df <- metric_long[metric_long$truth_target == target & metric_long$method %in% c(hbfss_method, baseline_method), , drop = FALSE]
  if (nrow(plot_df) == 0L) return(NULL)

  scenario_cols <- c(
    "de_fraction", "lfc_magnitude", "simulation_profile", "null_inflation",
    "lfc_label", "de_label", "inflation_label"
  )

  summary_long <- plot_df %>%
    group_by(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label, method) %>%
    summarise(mean_value = mean(.data[[metric]], na.rm = TRUE), .groups = "drop")

  hbfss_df <- summary_long[summary_long$method == hbfss_method, c(scenario_cols, "mean_value"), drop = FALSE]
  base_df <- summary_long[summary_long$method == baseline_method, c(scenario_cols, "mean_value"), drop = FALSE]

  if (nrow(hbfss_df) == 0L || nrow(base_df) == 0L) return(NULL)

  names(hbfss_df)[names(hbfss_df) == "mean_value"] <- "hbfss_mean"
  names(base_df)[names(base_df) == "mean_value"] <- "baseline_mean"

  summary <- left_join(hbfss_df, base_df, by = scenario_cols)
  summary <- summary[is.finite(summary$hbfss_mean) & is.finite(summary$baseline_mean), , drop = FALSE]
  if (nrow(summary) == 0L) return(NULL)

  summary$delta <- summary$hbfss_mean - summary$baseline_mean
  summary$label <- sprintf(
    "H %.2f\nD %.2f\nDelta %.2f",
    summary$hbfss_mean,
    summary$baseline_mean,
    summary$delta
  )

  ggplot(summary, aes(lfc_label, de_label, fill = delta)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = label), size = 2.05, lineheight = 0.88) +
    facet_wrap(~ inflation_label) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, name = fill_label) +
    labs(
      title = title,
      x = "Simulated effect profile",
      y = "True DE fraction",
      caption = "Tile text: H = HBFSS mean, D = DESeq2 BH mean, Delta = HBFSS - DESeq2."
    ) +
    manuscript_theme() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

plot_simulation_fdr_heatmap <- function(metric_long, target = "all_de") {
  methods <- c("DESeq2_BH", "Empirical_BH", "HBFSS_total", "HBFSS_raw")
  plot_df <- metric_long[metric_long$truth_target == target & metric_long$method %in% methods, , drop = FALSE]
  if (nrow(plot_df) == 0L) return(NULL)

  summary <- plot_df %>%
    group_by(method, de_label, lfc_label, inflation_label) %>%
    summarise(mean_fdr = mean(fdr, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      method = factor(method, levels = methods),
      label = sprintf("%.2f", mean_fdr)
    )

  ggplot(summary, aes(lfc_label, de_label, fill = mean_fdr)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = label), size = 2.15) +
    facet_grid(method ~ inflation_label) +
    scale_fill_gradient(low = "white", high = "#B2182B", name = "Mean FDR") +
    labs(
      title = "Simulation observed FDR by method",
      x = "Simulated effect profile",
      y = "True DE fraction",
      caption = paste0("Nominal alpha = ", alpha_level, ". Observed FDR is descriptive and should not be called formal FDR control unless supported by the setting.")
    ) +
    manuscript_theme() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

plot_simulation_weak_recall <- function(metric_long) {
  methods <- c("DESeq2_BH", "LessAbs", "HBFSS_weak_region", "HBFSS_raw")
  plot_df <- metric_long[metric_long$truth_target == "weak_de" & metric_long$method %in% methods, , drop = FALSE]
  if (nrow(plot_df) == 0L) return(NULL)

  summary <- plot_df %>%
    group_by(de_label, lfc_label, inflation_label, method) %>%
    summarise(
      mean_recall = mean(recall, na.rm = TRUE),
      sd_recall = sd(recall, na.rm = TRUE),
      .groups = "drop"
    )

  summary$method <- factor(summary$method, levels = methods)

  method_colors <- c(
    DESeq2_BH = "#999999",
    LessAbs = "#56B4E9",
    HBFSS_weak_region = class_colors[["Weak"]],
    HBFSS_raw = class_colors[["HBFSS"]]
  )

  ggplot(summary, aes(lfc_label, mean_recall, color = method, group = method)) +
    geom_line(linewidth = 0.60, position = position_dodge(width = 0.28)) +
    geom_point(size = 1.85, position = position_dodge(width = 0.28)) +
    geom_errorbar(
      aes(ymin = pmax(mean_recall - sd_recall, 0), ymax = pmin(mean_recall + sd_recall, 1)),
      width = 0.14,
      linewidth = 0.28,
      position = position_dodge(width = 0.28)
    ) +
    facet_grid(inflation_label ~ de_label) +
    scale_color_manual(values = method_colors[methods], breaks = methods, drop = FALSE, name = NULL) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.12))) +
    labs(
      title = "Weak-effect recall",
      x = "Sub-boundary simulated effect profile",
      y = "Recall",
      caption = paste0("Weak-mixture simulations draw true absolute LFC from ", simulation_weak_lfc_min, " to ", simulation_weak_lfc_max, " with n = ", simulation_n_samples_per_group, " per group.")
    ) +
    manuscript_theme() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

plot_simulation_threshold_stability <- function(metric_wide) {
  if (nrow(metric_wide) == 0L) return(NULL)

  plot_df <- metric_wide[is.finite(metric_wide$hbfss_cutoff) & !is.na(metric_wide$hbfss_cutoff), , drop = FALSE]
  if (nrow(plot_df) == 0L) return(NULL)

  ggplot(plot_df, aes(de_label, hbfss_cutoff, fill = de_label)) +
    geom_boxplot(outlier.size = 0.35, linewidth = 0.25, na.rm = TRUE) +
    facet_grid(inflation_label ~ lfc_label, scales = "free_y") +
    scale_fill_manual(values = scales::grey_pal(start = 0.35, end = 0.75)(length(unique(plot_df$de_label))), name = NULL) +
    labs(
      title = "HBFSS cutoff stability across simulation replicates",
      x = NULL,
      y = "HBFSS cutoff",
      caption = "Only replicates with valid HC thresholds are shown."
    ) +
    manuscript_theme() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.text = element_text(size = base_theme_size - 1.0)
    )
}

plot_simulation_discovery_counts <- function(metric_wide) {
  if (nrow(metric_wide) == 0L) return(NULL)

  long <- bind_rows(
    data.frame(metric_wide[, c("de_fraction", "lfc_magnitude", "simulation_profile", "null_inflation", "lfc_label", "de_label", "inflation_label")], method = "DESeq2_BH", count = metric_wide$n_standard),
    data.frame(metric_wide[, c("de_fraction", "lfc_magnitude", "simulation_profile", "null_inflation", "lfc_label", "de_label", "inflation_label")], method = "Empirical_BH", count = metric_wide$n_empirical_bh),
    data.frame(metric_wide[, c("de_fraction", "lfc_magnitude", "simulation_profile", "null_inflation", "lfc_label", "de_label", "inflation_label")], method = "HBFSS_raw", count = metric_wide$n_hbfss_raw),
    data.frame(metric_wide[, c("de_fraction", "lfc_magnitude", "simulation_profile", "null_inflation", "lfc_label", "de_label", "inflation_label")], method = "HBFSS_total", count = metric_wide$n_hbfss_total),
    data.frame(metric_wide[, c("de_fraction", "lfc_magnitude", "simulation_profile", "null_inflation", "lfc_label", "de_label", "inflation_label")], method = "HBFSS_weak_region", count = metric_wide$n_weak_hbfss)
  )

  long$method <- factor(long$method, levels = c("DESeq2_BH", "Empirical_BH", "HBFSS_raw", "HBFSS_total", "HBFSS_weak_region"))

  ggplot(long, aes(method, count, fill = method)) +
    geom_boxplot(outlier.size = 0.35, linewidth = 0.25, na.rm = TRUE) +
    facet_grid(inflation_label + de_label ~ lfc_label, scales = "free_y") +
    scale_fill_manual(
      values = c(
        DESeq2_BH = "#999999",
        Empirical_BH = treatment_color,
        HBFSS_raw = "#54278F",
        HBFSS_total = class_colors[["HBFSS"]],
        HBFSS_weak_region = class_colors[["Weak"]]
      ),
      name = NULL
    ) +
    labs(
      title = "Simulation discovery counts by method",
      x = NULL,
      y = "Discovered features",
      caption = NULL
    ) +
    manuscript_theme() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(size = base_theme_size - 1.0)
    )
}

run_sequence_simulation_validation <- function() {
  set.seed(simulation_seed)
  dir.create(simulation_dir, recursive = TRUE, showWarnings = FALSE)

  old_files <- list.files(simulation_dir, full.names = TRUE)
  if (length(old_files) > 0L) unlink(old_files, recursive = TRUE, force = TRUE)

  fixed_grid <- expand.grid(
    de_fraction = simulation_de_fractions,
    lfc_magnitude = simulation_lfc_magnitudes,
    null_inflation = simulation_null_inflation,
    simulation_profile = "fixed",
    stringsAsFactors = FALSE
  )

  weak_grid <- expand.grid(
    de_fraction = simulation_de_fractions,
    lfc_magnitude = NA_real_,
    null_inflation = simulation_null_inflation,
    simulation_profile = "weak_mixture",
    stringsAsFactors = FALSE
  )

  sim_grid <- bind_rows(fixed_grid, weak_grid)
  total_runs <- nrow(sim_grid) * simulation_n_reps

  message("Running simulation validation: ", total_runs, " planned replicates")

  metric_rows <- list()
  failure_rows <- list()
  feature_rows <- list()
  run_index <- 0L

  for (grid_i in seq_len(nrow(sim_grid))) {
    grid_row <- sim_grid[grid_i, , drop = FALSE]

    for (rep_i in seq_len(simulation_n_reps)) {
      run_index <- run_index + 1L

      if (run_index %% 10L == 0L || run_index == 1L || run_index == total_runs) {
        message("Simulation progress: ", run_index, "/", total_runs)
      }

      template <- data.frame(
        de_fraction = grid_row$de_fraction,
        lfc_magnitude = grid_row$lfc_magnitude,
        simulation_profile = grid_row$simulation_profile,
        null_inflation = grid_row$null_inflation,
        replicate = rep_i,
        stringsAsFactors = FALSE
      )
      template$lfc_label <- simulation_lfc_label(template$lfc_magnitude, template$simulation_profile)
      template$de_label <- paste0(template$de_fraction * 100, "% DE")
      template$inflation_label <- paste0(template$null_inflation * 100, "% null inflation")

      out <- tryCatch({
        sim_obj <- simulate_sequence_counts(
          n_features = simulation_n_features,
          n_samples = simulation_n_samples_per_group,
          base_mean = simulation_base_mean,
          disp_null = simulation_dispersion_null,
          de_fraction = grid_row$de_fraction,
          lfc_magnitude = grid_row$lfc_magnitude,
          disp_de = simulation_dispersion_de,
          lfc_profile = grid_row$simulation_profile
        )

        sequence_simulation_analysis(sim_obj, null_inflation = grid_row$null_inflation)
      }, error = function(e) {
        failure_rows[[length(failure_rows) + 1L]] <<- cbind(
          template,
          data.frame(error_message = conditionMessage(e), stringsAsFactors = FALSE)
        )
        NULL
      })

      if (is.null(out)) next

      metric_rows[[length(metric_rows) + 1L]] <- build_simulation_metric_rows(out, template)

      if (isTRUE(export_simulation_feature_results)) {
        feature_export <- out$results
        feature_export$de_fraction <- grid_row$de_fraction
        feature_export$lfc_magnitude <- grid_row$lfc_magnitude
        feature_export$simulation_profile <- grid_row$simulation_profile
        feature_export$null_inflation <- grid_row$null_inflation
        feature_export$replicate <- rep_i
        feature_rows[[length(feature_rows) + 1L]] <- feature_export
      }
    }
  }

  failures <- if (length(failure_rows) > 0L) bind_rows(failure_rows) else data.frame()
  save_csv(failures, file.path(simulation_dir, "Simulation_Failures.csv"))

  metric_long <- if (length(metric_rows) > 0L) bind_rows(metric_rows) else data.frame()
  if (nrow(metric_long) == 0L) {
    stop("Simulation validation produced no successful replicates. See Simulation_Failures.csv.", call. = FALSE)
  }

  metric_long$lfc_label <- factor(metric_long$lfc_label, levels = simulation_lfc_levels())
  metric_long$method <- factor(
    metric_long$method,
    levels = c("DESeq2_BH", "Empirical_BH", "GreaterAbs", "LessAbs", "HBFSS_raw", "HBFSS_total", "HBFSS_weak_region")
  )

  metric_wide <- simulation_wide_from_long(metric_long)
  metric_summary <- simulation_method_summary(metric_long)
  threshold_summary <- simulation_threshold_summary(metric_wide)

  save_csv(metric_long, file.path(simulation_dir, "Simulation_RunMetrics_Long.csv"))
  save_csv(metric_wide, file.path(simulation_dir, "Simulation_Raw_Wide.csv"))
  save_csv(metric_summary, file.path(simulation_dir, "Simulation_MethodSummary.csv"))
  save_csv(threshold_summary, file.path(simulation_dir, "Simulation_ThresholdSummary.csv"))

  # Backward-compatible filenames from the older draft.
  save_csv(metric_wide, file.path(simulation_dir, "Simulation_Raw.csv"))
  save_csv(
    metric_wide %>%
      group_by(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label) %>%
      summarise(
        n_replicates = n(),
        standard_f1_mean = mean(standard_f1, na.rm = TRUE),
        empirical_bh_f1_mean = mean(empirical_bh_f1, na.rm = TRUE),
        hbfss_f1_mean = mean(hbfss_f1, na.rm = TRUE),
        standard_fdr_mean = mean(standard_fdr, na.rm = TRUE),
        empirical_bh_fdr_mean = mean(empirical_bh_fdr, na.rm = TRUE),
        hbfss_fdr_mean = mean(hbfss_fdr, na.rm = TRUE),
        standard_weak_recall_mean = mean(standard_weak_recall, na.rm = TRUE),
        lessAbs_weak_recall_mean = mean(lessAbs_weak_recall, na.rm = TRUE),
        hbfss_weak_recall_mean = mean(hbfss_weak_recall, na.rm = TRUE),
        hc_valid_fraction = mean(hc_valid, na.rm = TRUE),
        hc_p_threshold_median = median(hc_p_threshold, na.rm = TRUE),
        hbfss_cutoff_median = median(hbfss_cutoff, na.rm = TRUE),
        .groups = "drop"
      ),
    file.path(simulation_dir, "Simulation_Summary.csv")
  )

  if (isTRUE(export_simulation_feature_results) && length(feature_rows) > 0L) {
    save_csv(bind_rows(feature_rows), file.path(simulation_dir, "Simulation_FeatureResults.csv"))
  }

  save_plot(
    plot_simulation_metric_boxplot(
      metric_long,
      metric = "f1",
      title = "Simulation F1 score by method",
      y_label = "F1",
      methods = c("DESeq2_BH", "Empirical_BH", "HBFSS_total", "HBFSS_raw"),
      target = "all_de"
    ),
    file.path(simulation_dir, "Simulation_F1_Boxplot.png"),
    width = 11.8,
    height = 9.2
  )

  save_plot(
    plot_simulation_metric_boxplot(
      metric_long,
      metric = "fdr",
      title = "Simulation observed FDR by method",
      y_label = "FDR",
      methods = c("DESeq2_BH", "Empirical_BH", "HBFSS_total", "HBFSS_raw"),
      target = "all_de",
      alpha_line = TRUE
    ),
    file.path(simulation_dir, "Simulation_FDR_Boxplot.png"),
    width = 11.8,
    height = 9.2
  )

  save_plot(
    plot_simulation_delta_heatmap(
      metric_long,
      metric = "f1",
      title = "Simulation F1 gain: HBFSS vs DESeq2 BH",
      fill_label = "Delta F1",
      hbfss_method = "HBFSS_total",
      baseline_method = "DESeq2_BH",
      target = "all_de"
    ),
    file.path(simulation_dir, "Simulation_F1_DeltaHeatmap.png"),
    width = 11.8,
    height = 6.4
  )

  save_plot(
    plot_simulation_fdr_heatmap(metric_long, target = "all_de"),
    file.path(simulation_dir, "Simulation_FDR_Heatmap.png"),
    width = 12.8,
    height = 7.2
  )

  save_plot(
    plot_simulation_weak_recall(metric_long),
    file.path(simulation_dir, "Simulation_WeakRecall.png"),
    width = 10.8,
    height = 6.2
  )

  save_plot(
    plot_simulation_threshold_stability(metric_wide),
    file.path(simulation_dir, "Simulation_CutoffStability.png"),
    width = 12.8,
    height = 7.2
  )

  save_plot(
    plot_simulation_discovery_counts(metric_wide),
    file.path(simulation_dir, "Simulation_DiscoveryCounts.png"),
    width = 12.8,
    height = 8.6
  )

  writeLines(capture.output(sessionInfo()), file.path(simulation_dir, "SessionInfo_Simulation.txt"))

  message("Simulation successful replicates: ", length(unique(paste(metric_long$de_fraction, metric_long$lfc_magnitude, metric_long$simulation_profile, metric_long$null_inflation, metric_long$replicate))))
  message("Simulation failures: ", nrow(failures))

  invisible(metric_summary)
}



# =============================================================================
# METHODS AND MANIFEST
# =============================================================================

write_sequence_methods <- function() {
  methods_lines <- c(
    "# SEQUENCE methods summary",
    "",
    "## Parameters used in this run",
    paste0("- alpha_level: ", alpha_level),
    paste0("- strong_alpha_level: ", strong_alpha_level),
    paste0("- weak_alpha_level: ", weak_alpha_level),
    paste0("- LFC boundary: ", lfc_boundary),
    paste0("- EVS cutoff mode: ", evs_cutoff_mode),
    paste0("- EVS top-N target: ", top_n_target),
    paste0("- HC invalid threshold: ", hc_invalid_at_or_above),
    paste0("- Calculation probability floor: ", format(calculation_probability_floor, scientific = TRUE)),
    paste0("- Plot probability floor: ", plot_probability_floor),
    paste0("- Comparisons: ", paste(comparison_table$comparison_name, collapse = ", ")),
    paste0("- Analysis tracks: ", paste(analysis_tracks, collapse = ", ")),
    paste0("- Figure DPI: ", figure_dpi, " with parallel PDF export"),
    "",
    "## Input and preprocessing",
    "Raw read-count matrices are imported for each RT/ZT comparison. Raw integer counts are retained for DESeq2 differential testing. NormEVS uses DESeq2 variance-stabilized expression with log2(normalized counts + 1) fallback. RawEVS uses log2(raw counts + 1).",
    "",
    "## Eigenvector splitting",
    paste0("PC1 is calculated separately in treatment and control samples. Features are ranked by absolute PC1 loading. The leading edge is the union of the top-", top_n_target, " treatment-ranked and top-", top_n_target, " control-ranked features. The remainder is every feature not in that union. This is fixed top-N EVS; no changepoint, AIC/BIC, or hidden fallback is used."),
    "",
    "## Differential expression",
    paste0("DESeq2 is run with design ~ condition, reference level untrt. Standard effects: BH-adjusted Wald padj < ", alpha_level, " and |apeglm-shrunken LFC| >= ", lfc_boundary, ". Strong effects: greaterAbs alternative-hypothesis calls at the same LFC boundary. Weak effects: lessAbs support, |shrunken LFC| < ", lfc_boundary, ", and HBFSS raw significance."),
    "",
    "## Empirical null and HBFSS",
    "Wald statistics are passed to fdrtool using statistic = normal. fndr is attempted first; pct0 with pct0 = 0.75 is used only if fndr fails or returns invalid eta0. Higher criticism supplies the dataset-specific empirical-p threshold unless invalid. HBFSS cutoff = -log10(HC threshold) x LFC boundary. HBFSS total is reported as standard DESeq2 plus HBFSS-only calls.",
    "",
    "## PCA and EVS figures",
    "PC1_Variance_EVS compares absolute PC1 score variance across Original, Lead, and Remainder datasets. PC1_Variance_Distributions shows per-feature contribution to PC1 score variance for All samples, Control samples, and Treatment samples, with treatment/control EVS cutoff lines.",
    "",
    "## Simulation validation",
    paste0("Simulation regeneration is disabled by default. If run_simulation_validation is TRUE, heterogeneous negative-binomial counts are simulated with feature-level mean and dispersion variation, sample-level library factors, known weak/strong truth labels, explicit replicate failure logging, and method-level confusion matrices across DE fractions (", paste(simulation_de_fractions, collapse = ", "), "), fixed LFC magnitudes (", paste(simulation_lfc_magnitudes, collapse = ", "), "), a weak-mixture profile with |LFC| drawn from ", simulation_weak_lfc_min, " to ", simulation_weak_lfc_max, ", null-inflation settings (", paste(simulation_null_inflation, collapse = ", "), "), n = ", simulation_n_samples_per_group, " samples per group, and ", simulation_n_reps, " replicates per condition. Seed: ", simulation_seed, ".")
  )

  writeLines(methods_lines, file.path(output_dir, "METHODS_SEQUENCE_PIPELINE.md"))
}

write_manifest <- function() {
  exported_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  exported_files <- exported_files[file.info(exported_files)$isdir %in% FALSE]

  manifest <- data.frame(
    file = sub(
      paste0("^", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "/?"),
      "",
      normalizePath(exported_files, winslash = "/", mustWork = FALSE)
    ),
    size_bytes = file.info(exported_files)$size,
    stringsAsFactors = FALSE
  )

  save_csv(manifest, file.path(output_dir, "Manifest.csv"))
  manifest
}

# =============================================================================
# PIPELINE EXECUTION
# =============================================================================

cleanup_rplots_pdf()

analysis_store <- list()
evs_store <- list()
comparison_store <- list()
summary_rows <- list()
failure_rows <- list()

run_store_save <- function(count_matrix, coldata, comparison_name, analysis_label, dataset_key) {
  fit <- run_deseq2_hbfss(count_matrix, coldata, comparison_name, analysis_label, dataset_key)

  save_csv(fit$results, analysis_table_path(comparison_name, analysis_label, "Results"))
  save_csv(concise_analysis_summary(fit$summary), analysis_table_path(comparison_name, analysis_label, "Summary"))

  key <- paste(comparison_name, analysis_label, sep = "__")
  analysis_store[[key]] <<- fit
  summary_rows[[key]] <<- fit$summary

  invisible(fit)
}

for (i in seq_len(nrow(comparison_table))) {
  comparison <- prepare_comparison(i)
  comparison_name <- comparison$comparison_name
  comparison_store[[comparison_name]] <- comparison

  message("\n=====================================================")
  message("Running comparison: ", comparison_name)
  message("=====================================================")

  tryCatch({
    run_store_save(comparison$count_matrix, comparison$coldata, comparison_name, "Raw", "Raw")

    for (track in analysis_tracks) {
      evs <- build_evs_split(comparison$count_matrix, comparison$coldata, comparison_name, track)
      evs_key <- paste(comparison_name, track, sep = "__")
      evs_store[[evs_key]] <- evs

      save_csv(evs$joint_rank, analysis_table_path(comparison_name, track, "EVS_Rank_Table"))
      save_csv(concise_evs_summary(evs$summary), analysis_table_path(comparison_name, track, "EVS_Summary"))

      run_store_save(evs$leading_matrix, comparison$coldata, comparison_name, paste0("Lead_", track), "Lead")
      run_store_save(evs$remainder_matrix, comparison$coldata, comparison_name, paste0("Rem_", track), "Rem")
    }
  }, error = function(e) {
    failure_rows[[comparison_name]] <<- data.frame(
      comparison_name = comparison_name,
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
    message("FAILED: ", comparison_name, ": ", conditionMessage(e))
  })
}

summary_df <- if (length(summary_rows) > 0L) bind_rows(summary_rows) else data.frame()

if (nrow(summary_df) > 0L) {
  save_csv(concise_analysis_summary(summary_df), file.path(output_dir, "Summary_Overall.csv"))
}

if (length(failure_rows) > 0L) {
  save_csv(bind_rows(failure_rows), file.path(output_dir, "Failures.csv"))
  fail_pipeline("One or more comparisons failed. See Failures.csv. Manifest was not written because the run is incomplete.")
}

if (nrow(summary_df) == 0L) {
  fail_pipeline("No successful analyses were completed. Manifest was not written because the run is incomplete.")
}

# =============================================================================
# CURATED FIGURE EXPORT
# =============================================================================

comparison_order <- comparison_table$comparison_name

raw_plots <- list()
for (comparison_name in comparison_order) {
  key <- paste(comparison_name, "Raw", sep = "__")
  if (!is.null(analysis_store[[key]])) {
    raw_plots[[comparison_name]] <- plot_volcano(analysis_store[[key]]$results, comparison_name)
  }
}
raw_panel <- arrange_with_one_legend(raw_plots, "Volcano: Raw", ncol = length(comparison_order))
save_plot(raw_panel, file.path(figure_dir, "Volcano_Raw.png"), width = 18.0, height = 5.8)

for (dataset_prefix in c("Lead", "Rem")) {
  plots <- list()

  for (track in analysis_tracks) {
    for (comparison_name in comparison_order) {
      analysis_label <- paste0(dataset_prefix, "_", track)
      key <- paste(comparison_name, analysis_label, sep = "__")

      if (!is.null(analysis_store[[key]])) {
        plots[[paste(comparison_name, track, sep = "_")]] <- plot_volcano(
          analysis_store[[key]]$results,
          paste0(comparison_name, "\n", track)
        )
      }
    }
  }

  panel_title <- if (dataset_prefix == "Lead") "Volcano: Lead" else "Volcano: Remainder"
  panel_file <- if (dataset_prefix == "Lead") "Volcano_Lead.png" else "Volcano_Rem.png"
  panel <- arrange_with_one_legend(plots, panel_title, ncol = length(comparison_order))
  save_plot(panel, file.path(figure_dir, panel_file), width = 18.0, height = 9.6)
}

count_plot <- plot_discovery_counts(summary_df)
save_plot(count_plot, file.path(figure_dir, "Counts.png"), width = 14.0, height = 8.0)

pca_raw_plots <- list()
for (comparison_name in comparison_order) {
  key <- paste(comparison_name, "Raw", sep = "__")
  if (!is.null(analysis_store[[key]])) {
    pca_raw_plots[[comparison_name]] <- plot_pca_support(analysis_store[[key]], comparison_name)
  }
}
pca_raw_panel <- arrange_with_one_legend(pca_raw_plots, "PCA: Raw", ncol = length(comparison_order))
save_plot(pca_raw_panel, file.path(figure_dir, "PCA_Raw.png"), width = 18.0, height = 5.8)

for (dataset_prefix in c("Lead", "Rem")) {
  pca_plots <- list()

  for (track in analysis_tracks) {
    for (comparison_name in comparison_order) {
      analysis_label <- paste0(dataset_prefix, "_", track)
      key <- paste(comparison_name, analysis_label, sep = "__")

      if (!is.null(analysis_store[[key]])) {
        pca_plots[[paste(comparison_name, track, sep = "_")]] <- plot_pca_support(
          analysis_store[[key]],
          paste0(comparison_name, "\n", track)
        )
      }
    }
  }

  panel_title <- if (dataset_prefix == "Lead") "PCA: Lead" else "PCA: Remainder"
  panel_file <- if (dataset_prefix == "Lead") "PCA_Lead.png" else "PCA_Rem.png"
  pca_panel <- arrange_with_one_legend(pca_plots, panel_title, ncol = length(comparison_order))
  save_plot(pca_panel, file.path(figure_dir, panel_file), width = 18.0, height = 9.4)
}

pca_var_rows <- list()
for (comparison_name in comparison_order) {
  for (analysis_label in c("Raw", "Lead_NormEVS", "Lead_RawEVS", "Rem_NormEVS", "Rem_RawEVS")) {
    key <- paste(comparison_name, analysis_label, sep = "__")
    if (!is.null(analysis_store[[key]])) {
      pca_var_rows[[key]] <- extract_pca_variance_rows(analysis_store[[key]], comparison_name, analysis_label)
    }
  }
}

pca_var_df <- if (length(pca_var_rows) > 0L) bind_rows(pca_var_rows) else data.frame()

if (nrow(pca_var_df) > 0L) {
  save_csv(
    rename_columns_existing(
      pca_var_df,
      c(
        comparison_name = "comparison",
        analysis_label = "analysis",
        analysis_pretty = "panel",
        dataset_group = "dataset",
        evs_mode = "evs_mode",
        component = "component",
        variance_explained_pct = "variance_explained_pct",
        score_variance = "score_variance"
      )
    ),
    file.path(output_dir, "Summary_PCA.csv")
  )

  pca_pc1_compare_df <- pca_pc1_evs_plot_df(pca_var_df)

  if (nrow(pca_pc1_compare_df) > 0L) {
    save_csv(
      rename_columns_existing(
        pca_pc1_compare_df,
        c(
          comparison_name = "comparison",
          analysis_label = "analysis",
          analysis_pretty = "panel",
          dataset_group = "dataset",
          evs_mode = "evs_mode",
          pc1_score_variance = "pc1_score_variance",
          original_pc1_score_variance = "original_pc1_score_variance",
          pc1_score_variance_ratio_to_original = "pc1_score_variance_ratio_to_original",
          pc1_score_variance_percent_of_original = "pc1_score_variance_percent_of_original"
        )
      ),
      file.path(output_dir, "Summary_PCA_PC1_ScoreVariance_EVS.csv")
    )
  }

  save_plot(
    plot_pc1_evs_variance_absolute(pca_var_df),
    file.path(figure_dir, "PC1_Variance_EVS.png"),
    width = 15.8,
    height = 8.8
  )
}

load_rows <- list()
cutoff_rows <- list()

for (track in analysis_tracks) {
  for (comparison_name in comparison_order) {
    key <- paste(comparison_name, track, sep = "__")
    comp_obj <- comparison_store[[comparison_name]]

    if (!is.null(evs_store[[key]]) && !is.null(comp_obj)) {
      tmp <- evs_loading_distribution_rows(
        full_count_matrix = comp_obj$count_matrix,
        evs_obj = evs_store[[key]],
        coldata = comp_obj$coldata,
        comparison_name = comparison_name,
        track = track
      )

      load_rows[[key]] <- tmp
      cut <- attr(tmp, "cutoffs")
      if (!is.null(cut) && nrow(cut) > 0L) cutoff_rows[[key]] <- cut
    }
  }
}

load_df <- if (length(load_rows) > 0L) bind_rows(load_rows) else data.frame()
cutoff_df <- if (length(cutoff_rows) > 0L) bind_rows(cutoff_rows) else data.frame()

if (nrow(load_df) > 0L) {
  attr(load_df, "cutoffs") <- cutoff_df

  save_csv(
    rename_columns_existing(
      load_df,
      c(
        comparison_name = "comparison",
        track = "evs_mode",
        dataset_group = "dataset",
        condition = "scope",
        pc1_loading_abs = "abs_pc1_loading",
        pc1_score_variance = "pc1_score_variance",
        pc1_variance_contribution = "pc1_variance_contribution"
      )
    ),
    file.path(output_dir, "Summary_PC1VarianceDistributions.csv")
  )

  if (nrow(cutoff_df) > 0L) {
    save_csv(
      rename_columns_existing(
        cutoff_df,
        c(
          comparison_name = "comparison",
          track = "evs_mode",
          dataset_group = "dataset",
          condition = "scope",
          pc1_loading_cutoff = "pc1_loading_cutoff",
          pc1_score_variance = "pc1_score_variance",
          pc1_variance_contribution_cutoff = "pc1_variance_contribution_cutoff"
        )
      ),
      file.path(output_dir, "Summary_PC1VarianceCutoffs.csv")
    )
  }

  save_plot(
    plot_loading_histograms(load_df, "PC1 variance-contribution frequency after EVS"),
    file.path(figure_dir, "PC1_Variance_Distributions.png"),
    width = 15.8,
    height = 10.4
  )
}

if (isTRUE(run_simulation_validation)) {
  simulation_summary <- run_sequence_simulation_validation()
  print(simulation_summary)
}

# =============================================================================
# COMPLETION
# =============================================================================

write_sequence_methods()
write_run_session_info()
cleanup_rplots_pdf()

manifest <- write_manifest()

cat("\n=====================================================\n")
cat("SEQUENCE pipeline complete.\n")
cat("Output directory: ", output_dir, "\n", sep = "")
cat("Exported files: ", nrow(manifest), "\n", sep = "")
cat("=====================================================\n\n")

print(summary_df)
