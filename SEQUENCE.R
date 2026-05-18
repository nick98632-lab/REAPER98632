#!/usr/bin/env Rscript

# =============================================================================
# SEQUENCE FINAL MANUSCRIPT PIPELINE
# Complete clean rewrite: DESeq2 + fixed top-N EVS + empirical-null HC/HBFSS
# =============================================================================
# This script is intentionally organized as a single executable manuscript
# pipeline. It does four things only:
#   1. Imports the WTTS-Seq raw count matrix.
#   2. Runs Original, NormEVS-Leading, NormEVS-Remainder, RawEVS-Leading,
#      and RawEVS-Remainder DESeq2/HBFSS analyses for each comparison.
#   3. Exports result tables, summaries, manuscript figures, session info,
#      methods notes, and a manifest.
#   4. Optionally runs simulation validation and optional Git commit/push.
# =============================================================================

options(stringsAsFactors = FALSE)
options(width = 140)

# =============================================================================
# USER SETTINGS
# =============================================================================

count_file_candidates <- c(
  file.path("data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"),
  "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
)

output_folder_name <- "manuscript_final_clean"
reset_output_dir <- FALSE

# The three alpha thresholds are intentionally matched at 0.10 so standard Wald,
# greaterAbs, and lessAbs calls are compared under the same false-discovery target.
alpha_standard <- 0.10
alpha_strong <- 0.10
alpha_weak <- 0.10
lfc_boundary <- 1.0

# Fixed top-N EVS. The union of treatment top-N and control top-N is the Leading Edge.
evs_top_n <- 5000L
evs_modes <- c("NormEVS", "RawEVS")

# Empirical-null / HC behavior. If strict_empirical_null = TRUE, fdrtool failure stops the run.
strict_empirical_null <- TRUE
hc_invalid_below <- 0.001
hc_invalid_at_or_above <- 0.95
calculation_p_floor <- .Machine$double.xmin
plot_p_floor <- 1e-16

# Figure behavior.
figure_dpi <- 600
base_theme_size <- 9
export_pdf_also <- TRUE
save_individual_figures <- TRUE
label_top_n_total <- 8L
label_top_n_per_class <- 2L

# Large raw per-feature PC1 variance export. Keep FALSE for GitHub-safe runs.
# The PC1 variance figure and cutoff summary still generate from in-memory data.
export_full_pc1_variance_distributions <- FALSE
remove_skipped_large_exports_from_disk <- TRUE

# Optional simulation and Git behavior.
# Simulation is ON, but configured as a fast smoke test that proves execution without stalling.
# For manuscript validation, set simulation_n_reps >= 25 and usually set simulation_use_apeglm_shrinkage <- TRUE.
run_simulation_validation <- TRUE
export_simulation_feature_results <- FALSE
simulation_use_apeglm_shrinkage <- FALSE
simulation_seed <- 42L
simulation_manuscript_min_reps <- 25L

# For a deeper validation run, increase n_features, n_reps, and the grids explicitly.
simulation_n_features <- 500L
simulation_n_samples_per_group <- 6L
simulation_base_mean <- 200
simulation_dispersion_null <- 0.10
simulation_dispersion_de <- 0.15
simulation_de_fractions <- c(0.10)
simulation_lfc_magnitudes <- c(1.00)
simulation_weak_lfc_min <- 0.20
simulation_weak_lfc_max <- 0.80
simulation_null_inflation <- c(0.00)
simulation_n_reps <- 1L
simulation_progress_every <- max(1L, as.integer(1L))
simulation_checkpoint_every <- max(1L, as.integer(1L))

# Optional Git push. This is GitHub-safe by default: no force-adding ignored exports.
git_push_after_success <- TRUE
git_remote_name <- "origin"
git_branch_name <- NA_character_
git_commit_message <- paste0("Update SEQUENCE manuscript outputs ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
git_pull_rebase_before_push <- TRUE
git_stage_pipeline_script <- TRUE
git_allow_no_change_success <- TRUE
git_fail_if_preexisting_staged_changes <- FALSE
git_max_file_size_bytes <- 95 * 1024^2
# TRUE repairs the common failed-push state where an unpushed local commit already contains
# an oversized export blob. Set FALSE in shared branches or CI if automatic history rewrite is not desired.
git_rewrite_unpushed_large_export_commits <- TRUE
git_exclude_large_export_patterns <- c(
  "exports/manuscript_final_clean/Summary_PC1VarianceDistributions.csv",
  "exports/manuscript_final_clean/Summary_PC1VarianceDistributions*.csv",
  "exports/manuscript_final_clean/simulation_validation/Simulation_FeatureResults.csv",
  "exports/manuscript_final_clean/simulation_validation/Simulation_FeatureResults*.csv"
)

# =============================================================================
# STUDY DESIGN
# =============================================================================

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
# PACKAGES
# =============================================================================

required_packages <- c(
  "DESeq2", "apeglm", "fdrtool", "ggplot2", "ggrepel", "dplyr",
  "gridExtra", "grid", "scales", "S4Vectors", "SummarizedExperiment"
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
# CONSTANTS FOR DISPLAY
# =============================================================================

class_levels <- c("BG", "Weak", "Strong", "Std", "HBFSS")
class_labels <- c(
  BG = "Background",
  Weak = "Weak effect",
  Strong = "Strong effect",
  Std = "Standard DESeq2",
  HBFSS = "HBFSS-only"
)
class_colors <- c(
  BG = "#BDBDBD",
  Weak = "#0072B2",
  Strong = "#E31A1C",
  Std = "#33A02C",
  HBFSS = "#6A3D9A"
)
class_shapes <- c(BG = 21, Weak = 24, Strong = 22, Std = 23, HBFSS = 25)
class_sizes <- c(BG = 0.55, Weak = 1.10, Strong = 1.15, Std = 1.10, HBFSS = 1.15)
class_alphas <- c(BG = 0.24, Weak = 0.92, Strong = 0.95, Std = 0.90, HBFSS = 0.95)

threshold_color <- "#A65628"
treatment_color <- "#1F78B4"
control_color <- "#4D4D4D"

# =============================================================================
# PATHS, LOGGING, AND GENERAL HELPERS
# =============================================================================

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0L) return(NA_character_)
  normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = FALSE)
}

find_repo_root <- function() {
  candidates <- unique(c(dirname(script_path()), getwd()))
  candidates <- candidates[!is.na(candidates) & dir.exists(candidates)]

  walk_up <- function(start) {
    current <- normalizePath(start, winslash = "/", mustWork = TRUE)
    repeat {
      if (dir.exists(file.path(current, ".git"))) return(current)
      if (file.exists(file.path(current, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"))) return(current)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
    NA_character_
  }

  hits <- vapply(candidates, walk_up, character(1))
  hits <- hits[!is.na(hits)]
  if (length(hits) > 0L) return(hits[1])
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

repo_root <- find_repo_root()
output_dir <- file.path(repo_root, "exports", output_folder_name)
figure_dir <- file.path(output_dir, "manuscript_figures")
simulation_dir <- file.path(output_dir, "simulation_validation")
log_file <- file.path(output_dir, "Pipeline_Log.txt")

if (isTRUE(reset_output_dir) && dir.exists(output_dir)) {
  unlink(output_dir, recursive = TRUE, force = TRUE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(simulation_dir, recursive = TRUE, showWarnings = FALSE)

log_message <- function(...) {
  txt <- paste0(...)
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", stamp, "] ", txt)
  message(line)
  cat(line, "\n", file = log_file, append = TRUE)
}

cleanup_rplots_pdf <- function() {
  stray <- file.path(repo_root, "Rplots.pdf")
  if (file.exists(stray)) unlink(stray, force = TRUE)
}

close_open_graphics_devices <- function() {
  while (grDevices::dev.cur() > 1L) {
    try(grDevices::dev.off(), silent = TRUE)
  }
  invisible(TRUE)
}

cleanup_run_artifacts <- function() {
  close_open_graphics_devices()
  cleanup_rplots_pdf()
  invisible(TRUE)
}

resolve_file <- function(candidates, label) {
  full_candidates <- unique(c(file.path(repo_root, candidates), candidates))
  hits <- full_candidates[file.exists(full_candidates)]
  if (length(hits) == 0L) {
    stop("Could not find ", label, ". Tried: ", paste(full_candidates, collapse = " | "), call. = FALSE)
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
  mat <- as.matrix(x)
  suppressWarnings(storage.mode(mat) <- "numeric")

  if (any(is.na(mat) | !is.finite(mat))) stop(label, " has NA or non-finite count values.", call. = FALSE)
  if (any(mat < 0)) stop(label, " has negative count values.", call. = FALSE)

  rounded <- round(mat)
  if (any(abs(mat - rounded) > 1e-6)) {
    warning(label, " had non-integer count values; values were rounded for DESeq2.", call. = FALSE)
  }

  storage.mode(rounded) <- "integer"
  rownames(rounded) <- row_ids
  colnames(rounded) <- col_ids
  rounded
}

clip_probability <- function(x, floor_value = calculation_p_floor) {
  y <- suppressWarnings(as.numeric(x))
  y[!is.finite(y)] <- NA_real_
  ok <- !is.na(y)
  y[ok] <- pmin(pmax(y[ok], floor_value), 1 - 1e-12)
  y
}

safe_neglog10 <- function(p, floor_value = plot_p_floor) {
  -log10(clip_probability(p, floor_value = floor_value))
}

format_compact_number <- function(x, digits = 3) {
  y <- suppressWarnings(as.numeric(x))
  out <- rep("NA", length(y))
  ok <- is.finite(y) & !is.na(y)
  ay <- abs(y[ok])
  val <- y[ok]
  out[ok] <- ifelse(
    ay >= 1e9, paste0(signif(val / 1e9, digits), "B"),
    ifelse(
      ay >= 1e6, paste0(signif(val / 1e6, digits), "M"),
      ifelse(
        ay >= 1e3, paste0(signif(val / 1e3, digits), "K"),
        ifelse(ay > 0 & ay < 0.001, formatC(val, format = "e", digits = digits - 1), as.character(signif(val, digits)))
      )
    )
  )
  out
}

mean_finite <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y <- y[is.finite(y) & !is.na(y)]
  if (length(y) == 0L) return(NA_real_)
  mean(y)
}

median_finite <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y <- y[is.finite(y) & !is.na(y)]
  if (length(y) == 0L) return(NA_real_)
  stats::median(y)
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
  if (length(hit) > 0L) names(df)[match(hit, names(df))] <- unname(map[hit])
  df
}

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_theme_size + 0.5, hjust = 0.5),
      plot.subtitle = element_text(size = base_theme_size - 1.2, hjust = 0.5),
      plot.caption = element_text(size = base_theme_size - 3.5, color = "grey30", hjust = 0.5),
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

is_grid_object <- function(x) {
  inherits(x, c("grob", "gtable", "gTree"))
}

draw_to_device <- function(plot_obj) {
  grid::grid.newpage()
  if (inherits(plot_obj, "ggplot")) {
    print(plot_obj)
  } else if (is_grid_object(plot_obj)) {
    grid::grid.draw(plot_obj)
  } else {
    stop("Unsupported figure object class: ", paste(class(plot_obj), collapse = ", "), call. = FALSE)
  }
}

save_figure <- function(plot_obj, path, width, height, export_pdf = export_pdf_also) {
  if (is.null(plot_obj)) return(invisible(NULL))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  ok_png <- FALSE
  ok_pdf <- FALSE

  tryCatch({
    if (inherits(plot_obj, "ggplot")) {
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
    } else {
      grDevices::png(filename = path, width = width, height = height, units = "in", res = figure_dpi, bg = "white")
      on.exit(grDevices::dev.off(), add = TRUE)
      draw_to_device(plot_obj)
      grDevices::dev.off()
      on.exit(NULL, add = FALSE)
    }
    ok_png <- TRUE
  }, error = function(e) {
    stop("PNG figure export failed for ", basename(path), ": ", conditionMessage(e), call. = FALSE)
  })

  if (isTRUE(export_pdf)) {
    pdf_path <- sub("\\.[^.]+$", ".pdf", path)
    tryCatch({
      grDevices::pdf(file = pdf_path, width = width, height = height, onefile = TRUE, useDingbats = FALSE)
      on.exit(grDevices::dev.off(), add = TRUE)
      draw_to_device(plot_obj)
      grDevices::dev.off()
      on.exit(NULL, add = FALSE)
      ok_pdf <- TRUE
    }, error = function(e) {
      stop("PDF figure export failed for ", basename(pdf_path), ": ", conditionMessage(e), call. = FALSE)
    })
  }

  log_message("Figure saved: ", path, if (ok_pdf) " + PDF" else "")
  invisible(ok_png)
}

safe_save_figure <- function(plot_obj, path, width, height) {
  tryCatch({
    save_figure(plot_obj, path, width, height)
    data.frame(file = basename(path), status = "saved", error_message = NA_character_, stringsAsFactors = FALSE)
  }, error = function(e) {
    log_message("FIGURE FAILED: ", basename(path), " | ", conditionMessage(e))
    data.frame(file = basename(path), status = "failed", error_message = conditionMessage(e), stringsAsFactors = FALSE)
  })
}

get_shared_legend <- function(p) {
  if (is.null(p)) return(NULL)
  grob <- ggplotGrob(p + theme(legend.position = "bottom"))
  idx <- which(vapply(grob$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(idx) == 0L) return(NULL)
  grob$grobs[[idx[1]]]
}

with_panel_mode <- function(expr) {
  old <- getOption("sequence.panel.mode", FALSE)
  options(sequence.panel.mode = TRUE)
  on.exit(options(sequence.panel.mode = old), add = TRUE)
  force(expr)
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

  if (is.null(legend)) return(gridExtra::arrangeGrob(body, ncol = 1, top = title_grob))
  gridExtra::arrangeGrob(body, legend, ncol = 1, heights = c(12, 0.95), top = title_grob)
}

write_run_session_info <- function() {
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "SessionInfo.txt"))
  invisible(TRUE)
}

fail_pipeline <- function(message_text, status = 1L) {
  try(write_run_session_info(), silent = TRUE)
  try(cleanup_run_artifacts(), silent = TRUE)
  log_message(message_text)
  quit(save = "no", status = status, runLast = FALSE)
}

options(error = function() {
  try(write_run_session_info(), silent = TRUE)
  try(cleanup_run_artifacts(), silent = TRUE)
  traceback(2)
  quit(save = "no", status = 1L, runLast = FALSE)
})

# =============================================================================
# DATA IMPORT
# =============================================================================

read_count_data <- function() {
  count_file <- resolve_file(count_file_candidates, "WTTS count matrix")
  log_message("Using count file: ", count_file)
  log_message("Output directory: ", output_dir)

  raw <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)

  stop_missing_columns(raw, c("OrigID", "Symbol"), "Count matrix")
  stop_missing_columns(raw, sample_metadata$id, "Count matrix")

  raw$OrigID <- trimws(as.character(raw$OrigID))
  raw$Symbol <- trimws(as.character(raw$Symbol))

  for (sample_id in sample_metadata$id) {
    raw[[sample_id]] <- clean_count_column(raw[[sample_id]])
  }

  count_block <- raw[, sample_metadata$id, drop = FALSE]
  keep <- !is.na(raw$OrigID) & nzchar(raw$OrigID) & rowSums(is.na(count_block)) == 0
  raw <- raw[keep, , drop = FALSE]

  raw$feature_id <- make.unique(raw$OrigID, sep = "_dup")
  rownames(raw) <- raw$feature_id

  annotation <- data.frame(
    feature_id = raw$feature_id,
    orig_id = raw$OrigID,
    gene_symbol = ifelse(nzchar(raw$Symbol), raw$Symbol, NA_character_),
    stringsAsFactors = FALSE
  )
  annotation <- dplyr::distinct(annotation, feature_id, .keep_all = TRUE)

  list(count_df = raw, annotation_df = annotation)
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
    stop("Comparison ", row$comparison_name, " does not have enough treatment/control samples.", call. = FALSE)
  }

  coldata <- sample_metadata[sample_ids, "condition", drop = FALSE]
  counts_mat <- count_df[, sample_ids, drop = FALSE]
  rownames(counts_mat) <- count_df$feature_id
  counts_mat <- as_integer_count_matrix(counts_mat, paste0(row$comparison_name, " count matrix"))

  if (!identical(colnames(counts_mat), rownames(coldata))) {
    stop("Sample order mismatch for ", row$comparison_name, call. = FALSE)
  }

  list(comparison_name = row$comparison_name, count_matrix = counts_mat, coldata = coldata)
}

# =============================================================================
# EVS
# =============================================================================

make_vst_matrix <- function(count_matrix, coldata) {
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)

  transformed <- tryCatch(
    as.matrix(SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))),
    error = function(e) NULL
  )

  if (is.null(transformed)) {
    transformed <- log2(as.matrix(DESeq2::counts(dds, normalized = TRUE)) + 1)
  }

  storage.mode(transformed) <- "numeric"
  transformed
}

get_evs_matrix <- function(count_matrix, coldata, mode) {
  if (identical(mode, "RawEVS")) {
    mat <- matrix(as.numeric(count_matrix), nrow = nrow(count_matrix), dimnames = dimnames(count_matrix))
    return(log2(mat + 1))
  }
  if (identical(mode, "NormEVS")) return(make_vst_matrix(count_matrix, coldata))
  stop("Unknown EVS mode: ", mode, call. = FALSE)
}

condition_pc1_rank <- function(evs_matrix, sample_ids, condition_label) {
  x <- evs_matrix[, sample_ids, drop = FALSE]
  storage.mode(x) <- "numeric"

  keep <- rowSums(is.finite(x) & !is.na(x)) == ncol(x)
  keep <- keep & apply(x, 1, stats::var, na.rm = TRUE) > 0
  if (!any(keep)) stop("No variable EVS features for ", condition_label, ".", call. = FALSE)

  x <- x[keep, , drop = FALSE]
  pca <- stats::prcomp(t(x), center = TRUE, scale. = FALSE)
  loading <- pca$rotation[, 1]
  pc1_score_variance <- as.numeric(pca$sdev[1]^2)

  rank_df <- data.frame(
    feature_id = names(loading),
    pc1_loading = as.numeric(loading),
    pc1_loading_abs = abs(as.numeric(loading)),
    pc1_score_variance = pc1_score_variance,
    pc1_variance_contribution = as.numeric(loading)^2 * pc1_score_variance,
    stringsAsFactors = FALSE
  )

  rank_df <- rank_df[order(rank_df$pc1_loading_abs, decreasing = TRUE), , drop = FALSE]
  rank_df$evs_rank <- seq_len(nrow(rank_df))
  top_n_used <- min(evs_top_n, nrow(rank_df))
  rank_df$evs_selected <- rank_df$evs_rank <= top_n_used
  rank_df$top_n_used <- top_n_used
  rank_df$loading_cutoff_at_top_n <- rank_df$pc1_loading_abs[top_n_used]
  rank_df$pc1_contribution_cutoff_at_top_n <- rank_df$pc1_variance_contribution[top_n_used]

  list(pca = pca, rank_df = rank_df, top_n_used = top_n_used)
}

build_evs_split <- function(count_matrix, coldata, comparison_name, mode) {
  evs_matrix <- get_evs_matrix(count_matrix, coldata, mode)
  treatment_samples <- rownames(coldata)[coldata$condition == "trt"]
  control_samples <- rownames(coldata)[coldata$condition == "untrt"]
  all_features <- rownames(count_matrix)

  trt <- condition_pc1_rank(evs_matrix, treatment_samples, "Treatment")
  ctl <- condition_pc1_rank(evs_matrix, control_samples, "Control")

  trt_top <- trt$rank_df$feature_id[trt$rank_df$evs_selected]
  ctl_top <- ctl$rank_df$feature_id[ctl$rank_df$evs_selected]

  leading_ids <- intersect(union(trt_top, ctl_top), all_features)
  remainder_ids <- setdiff(all_features, leading_ids)

  if (length(leading_ids) == 0L) stop("EVS leading edge is empty for ", comparison_name, " ", mode, ".", call. = FALSE)
  if (length(remainder_ids) == 0L) {
    stop(
      "EVS remainder is empty for ", comparison_name, " ", mode,
      ". Fixed evs_top_n=", evs_top_n, " selected every feature. Reduce evs_top_n explicitly.",
      call. = FALSE
    )
  }

  trt_rank <- trt$rank_df
  ctl_rank <- ctl$rank_df
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
    evs_mode = mode,
    evs_cutoff_mode = "fixed_top_n_union",
    top_n_target = evs_top_n,
    treatment_top_n_used = trt$top_n_used,
    control_top_n_used = ctl$top_n_used,
    treatment_loading_cutoff = trt$rank_df$loading_cutoff_at_top_n[1],
    control_loading_cutoff = ctl$rank_df$loading_cutoff_at_top_n[1],
    leading_edge_n = length(leading_ids),
    remainder_n = length(remainder_ids),
    evs_input = ifelse(mode == "NormEVS", "DESeq2 VST matrix", "log2(raw counts + 1) matrix"),
    stringsAsFactors = FALSE
  )

  list(
    comparison_name = comparison_name,
    evs_mode = mode,
    treatment_rank = trt,
    control_rank = ctl,
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

  run_fit <- function(method) {
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

  method_used <- "fndr"
  fit <- run_fit("fndr")
  if (is.null(fit)) {
    method_used <- "pct0"
    fit <- run_fit("pct0")
  }

  if (is.null(fit)) {
    if (isTRUE(strict_empirical_null)) {
      stop(label, ": fdrtool empirical-null fit failed with fndr and pct0.", call. = FALSE)
    }
    method_used <- "theoretical_normal_fallback"
    p_raw <- 2 * stats::pnorm(-abs(z))
    fit <- list(pval = p_raw, qval = p.adjust(p_raw, method = "BH"), lfdr = rep(NA_real_, length(z)))
  }

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
    empirical_bh = empirical_bh,
    empirical_null_method = method_used
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
  if (threshold < hc_invalid_below) return(NA_real_)
  if (threshold >= hc_invalid_at_or_above) return(NA_real_)
  threshold
}

classify_results <- function(df, hc_p, hbfss_cutoff) {
  df$standard_flag <- !is.na(df$padj) &
    df$padj < alpha_standard &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= lfc_boundary

  df$strong_alt_flag <- !is.na(df$greaterAbs_padj) &
    df$greaterAbs_padj < alpha_strong &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= lfc_boundary

  df$lessAbs_alt_flag <- !is.na(df$lessAbs_padj) &
    df$lessAbs_padj < alpha_weak &
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

  df$strong_flag <- df$strong_alt_flag
  df$hbfss_flag <- df$standard_flag | df$hbfss_raw_flag
  df$weak_region_hbfss_flag <- df$lessAbs_alt_flag & df$hbfss_raw_flag
  df$weak_flag <- df$weak_region_hbfss_flag

  df$display_strong <- df$strong_flag
  df$display_standard <- df$standard_flag & !df$display_strong
  df$display_weak <- df$weak_flag & !df$display_strong & !df$display_standard
  df$display_hbfss <- df$hbfss_flag & !df$standard_flag & !df$display_weak & !df$display_strong

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
  log_message("DESeq2/HBFSS: ", label, " | features=", nrow(count_matrix))

  dds <- DESeqDataSetFromMatrix(
    countData = as_integer_count_matrix(count_matrix, paste0(label, " DESeq2 input")),
    colData = coldata,
    design = ~ condition
  )

  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, betaPrior = FALSE, quiet = TRUE)
  coef_name <- condition_coef_name(dds)

  standard <- results(dds, contrast = c("condition", "trt", "untrt"), alpha = alpha_standard)
  strong <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "greaterAbs",
    alpha = alpha_strong
  )
  weak <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "lessAbs",
    alpha = alpha_weak
  )
  shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)

  df <- as.data.frame(standard)
  df$feature_id <- rownames(df)

  empirical <- fit_empirical_null(df$stat, label)
  df$empirical_p <- empirical$empirical_p
  df$empirical_q <- empirical$empirical_q
  df$empirical_lfdr <- empirical$empirical_lfdr
  df$empirical_bh <- empirical$empirical_bh
  df$empirical_null_method <- empirical$empirical_null_method

  df$empirical_p_calc <- clip_probability(df$empirical_p, floor_value = calculation_p_floor)
  df$empirical_p_plot <- clip_probability(df$empirical_p, floor_value = plot_p_floor)
  df$neglog10_empirical_p_calc <- safe_neglog10(df$empirical_p, floor_value = calculation_p_floor)
  df$neglog10_empirical_p_plot <- safe_neglog10(df$empirical_p, floor_value = plot_p_floor)
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
  # The HBFSS cutoff scales with the LFC boundary so changing the effect-size gate
  # changes both the vertical and geometric HBFSS decision surface deliberately.
  hbfss_cutoff <- if (is.na(hc_p)) NA_real_ else -log10(hc_p) * lfc_boundary
  df <- classify_results(df, hc_p, hbfss_cutoff)

  normalized_counts <- as.data.frame(counts(dds, normalized = TRUE))
  normalized_counts$feature_id <- rownames(normalized_counts)

  dispersion_df <- as.data.frame(S4Vectors::mcols(dds))
  dispersion_df$feature_id <- rownames(dispersion_df)
  keep_disp <- intersect(c("feature_id", "dispGeneEst", "dispFit", "dispersion", "dispIter", "dispOutlier"), colnames(dispersion_df))
  dispersion_df <- dispersion_df[, keep_disp, drop = FALSE]

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
    "greaterAbs_pvalue", "greaterAbs_padj", "lessAbs_pvalue", "lessAbs_padj",
    "empirical_null_method", "empirical_p", "empirical_p_calc", "empirical_p_plot",
    "empirical_bh", "empirical_q", "empirical_lfdr",
    "neglog10_empirical_p_calc", "neglog10_empirical_p_plot", "neglog10_empirical_p",
    "HBFSS", "hc_p_threshold_dataset", "hbfss_threshold_dataset",
    "hc_pass", "standard_flag", "strong_alt_flag", "strong_flag",
    "lessAbs_alt_flag", "weak_region_hbfss_flag", "weak_flag",
    "hbfss_raw_flag", "hbfss_flag",
    "display_standard", "display_strong", "display_weak", "display_hbfss", "final_class"
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
    n_display_weak = sum(df$display_weak, na.rm = TRUE),
    n_display_strong = sum(df$display_strong, na.rm = TRUE),
    n_display_standard = sum(df$display_standard, na.rm = TRUE),
    n_display_hbfss = sum(df$display_hbfss, na.rm = TRUE),
    hc_p_threshold = hc_p,
    hbfss_threshold = hbfss_cutoff,
    empirical_null_method = unique(df$empirical_null_method)[1],
    alpha_standard = alpha_standard,
    alpha_strong = alpha_strong,
    alpha_weak = alpha_weak,
    lfc_boundary = lfc_boundary,
    stringsAsFactors = FALSE
  )

  list(dds = dds, results = df, summary = summary)
}

concise_analysis_summary <- function(df) {
  map <- c(
    comparison_name = "comparison", analysis_label = "analysis", dataset_key = "dataset",
    n_features = "n_features", n_lessAbs_alt = "lessAbs_significant",
    n_weak = "weak_effect_hbfss_region", n_strong_alt = "greaterAbs_significant",
    n_strong = "strong_effect", n_standard = "standard_effect",
    n_hbfss_raw = "hbfss_raw_geometric", n_hbfss_total = "hbfss_total",
    n_hbfss_overlap_standard = "hbfss_overlap_standard", n_hbfss_only = "hbfss_only",
    n_display_weak = "display_weak_effect", n_display_strong = "display_strong_effect",
    n_display_standard = "display_standard_only", n_display_hbfss = "display_hbfss_only",
    hc_p_threshold = "hc_p_threshold", hbfss_threshold = "hbfss_cutoff",
    empirical_null_method = "empirical_null_method", alpha_standard = "wald_bh_alpha",
    alpha_strong = "greaterAbs_bh_alpha", alpha_weak = "lessAbs_bh_alpha", lfc_boundary = "lfc_boundary"
  )
  out <- rename_columns_existing(df, map)
  keep <- c(
    "comparison", "analysis", "dataset", "n_features",
    "lessAbs_significant", "weak_effect_hbfss_region", "greaterAbs_significant", "strong_effect",
    "standard_effect", "hbfss_raw_geometric", "hbfss_total", "hbfss_overlap_standard", "hbfss_only",
    "display_weak_effect", "display_strong_effect", "display_standard_only", "display_hbfss_only",
    "hc_p_threshold", "hbfss_cutoff", "empirical_null_method",
    "wald_bh_alpha", "greaterAbs_bh_alpha", "lessAbs_bh_alpha", "lfc_boundary"
  )
  out[, intersect(keep, names(out)), drop = FALSE]
}

concise_evs_summary <- function(df) {
  map <- c(
    comparison_name = "comparison", evs_mode = "evs_mode", evs_cutoff_mode = "evs_cutoff_mode",
    top_n_target = "top_n_target", treatment_top_n_used = "treatment_top_n",
    control_top_n_used = "control_top_n", treatment_loading_cutoff = "treatment_loading_cutoff",
    control_loading_cutoff = "control_loading_cutoff", leading_edge_n = "leading_edge_features",
    remainder_n = "remainder_features", evs_input = "input_matrix"
  )
  out <- rename_columns_existing(df, map)
  keep <- c(
    "comparison", "evs_mode", "evs_cutoff_mode", "top_n_target",
    "treatment_top_n", "control_top_n", "treatment_loading_cutoff", "control_loading_cutoff",
    "leading_edge_features", "remainder_features", "input_matrix"
  )
  out[, intersect(keep, names(out)), drop = FALSE]
}

# =============================================================================
# VOLCANO FIGURES
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
    picked[[cls]] <- sub[seq_len(min(label_top_n_per_class, nrow(sub))), , drop = FALSE]
  }

  labels <- if (length(picked) > 0L) do.call(rbind, picked) else labels[0, , drop = FALSE]
  if (nrow(labels) == 0L) return(labels)
  labels <- labels[order(match(as.character(labels$final_class), class_priority), -labels$HBFSS, labels$empirical_p, na.last = TRUE), , drop = FALSE]
  labels[seq_len(min(label_top_n_total, nrow(labels))), , drop = FALSE]
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
    data.frame(x = -rev(x_abs), y = rev(y), side = "negative", stringsAsFactors = FALSE),
    data.frame(x = x_abs, y = y, side = "positive", stringsAsFactors = FALSE)
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
  threshold_label_y <- if (!is.na(hc_y) && is.finite(hc_y)) max(0.06 * y_limit, hc_y - 0.06 * y_limit) else 0.08 * y_limit

  caption <- paste0(
    "Weak=", sum(df$weak_flag, na.rm = TRUE),
    "  Strong=", sum(df$strong_flag, na.rm = TRUE),
    "  Std=", sum(df$standard_flag, na.rm = TRUE),
    "  HBFSS total=", sum(df$hbfss_flag, na.rm = TRUE),
    "  HBFSS-only=", sum(df$display_hbfss, na.rm = TRUE)
  )
  in_panel_mode <- isTRUE(getOption("sequence.panel.mode", FALSE))
  show_caption <- !in_panel_mode
  show_gene_labels <- !in_panel_mode
  show_threshold_labels <- !in_panel_mode

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
      fill = "none", shape = "none", size = "none", alpha = "none"
    ) +
    labs(title = title, x = "Shrunken log2 fold change", y = expression(-log[10]("empirical p")), caption = if (show_caption) caption else NULL) +
    coord_cartesian(ylim = c(0, y_limit), clip = "off") +
    manuscript_theme()

  if (!is.na(hc_y) && is.finite(hc_y)) {
    p <- p + geom_hline(yintercept = hc_y, linetype = "dotted", linewidth = 0.60, color = threshold_color)
    if (show_threshold_labels) {
      p <- p + annotate(
        "text", x = x_min + 0.03 * x_span, y = threshold_label_y,
        label = paste0("HC p=", signif(hc_p, 3)), hjust = 0, vjust = 1,
        size = 1.90, color = threshold_color
      )
    }
  }

  if (!is.null(boundary)) {
    p <- p + geom_line(
      data = boundary,
      aes(x = x, y = y, group = side),
      inherit.aes = FALSE,
      color = class_colors[["HBFSS"]],
      linewidth = 0.75
    )
  }

  if (show_threshold_labels && !is.na(hbfss_cutoff) && is.finite(hbfss_cutoff)) {
    p <- p + annotate(
      "text", x = x_min + 0.34 * x_span, y = threshold_label_y,
      label = paste0("HBFSS cutoff=", signif(hbfss_cutoff, 3)), hjust = 0, vjust = 1,
      size = 1.90, color = class_colors[["HBFSS"]]
    )
  }

  if (show_gene_labels && nrow(label_df) > 0L) {
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

# =============================================================================
# PCA AND EVS FIGURES
# =============================================================================

pca_support_components <- function(fit_obj) {
  dds <- fit_obj$dds
  if (is.null(dds)) return(NULL)

  x <- tryCatch(
    as.matrix(SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))),
    error = function(e) NULL
  )
  if (is.null(x)) x <- log2(as.matrix(DESeq2::counts(dds, normalized = TRUE)) + 1)
  storage.mode(x) <- "numeric"

  keep <- rowSums(is.finite(x) & !is.na(x)) == ncol(x)
  keep <- keep & apply(x, 1, stats::var, na.rm = TRUE) > 0
  x <- x[keep, , drop = FALSE]
  if (nrow(x) < 2L || ncol(x) < 3L) return(NULL)

  pca <- stats::prcomp(t(x), center = TRUE, scale. = FALSE)
  total_variance <- sum(pca$sdev^2)
  var_pct <- round((pca$sdev^2 / total_variance) * 100, 1)

  pca_df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    condition = as.character(SummarizedExperiment::colData(dds)[rownames(pca$x), "condition"]),
    stringsAsFactors = FALSE
  )
  pca_df$condition <- factor(ifelse(pca_df$condition == "trt", "Treatment", "Control"), levels = c("Control", "Treatment"))

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
    geom_point(size = 2.35, color = "grey15", stroke = 0.35) +
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
    labs(title = title, x = paste0("PC1 ", var_pct[1], "%"), y = paste0("PC2 ", var_pct[2], "%")) +
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

  dataset_group <- if (analysis_label == "Raw") "Original" else if (grepl("^Lead_", analysis_label)) "Lead" else if (grepl("^Rem_", analysis_label)) "Remainder" else analysis_label
  evs_mode <- if (analysis_label == "Raw") "Original" else sub("^[^_]+_", "", analysis_label)

  data.frame(
    comparison_name = comparison_name,
    analysis_label = analysis_label,
    analysis_pretty = analysis_label_pretty(analysis_label),
    dataset_group = factor(dataset_group, levels = c("Original", "Lead", "Remainder")),
    evs_mode = factor(evs_mode, levels = c("Original", evs_modes)),
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
  split_rows <- pc1[pc1$dataset_group %in% c("Lead", "Remainder") & pc1$evs_mode %in% evs_modes, , drop = FALSE]
  if (nrow(original_rows) == 0L || nrow(split_rows) == 0L) return(data.frame())

  out_rows <- list()
  for (mode_name in evs_modes) {
    for (comp_name in unique(as.character(split_rows$comparison_name))) {
      original_one <- original_rows[original_rows$comparison_name == comp_name, , drop = FALSE]
      split_one <- split_rows[split_rows$comparison_name == comp_name & as.character(split_rows$evs_mode) == mode_name, , drop = FALSE]
      if (nrow(original_one) == 0L || nrow(split_one) == 0L) next

      original_one <- original_one[1, , drop = FALSE]
      original_one$evs_mode <- mode_name
      original_one$analysis_label <- paste0("Original_for_", mode_name)
      original_one$analysis_pretty <- "Original"

      block <- bind_rows(original_one, split_one)
      original_variance <- suppressWarnings(as.numeric(original_one$score_variance[1]))
      block$original_pc1_score_variance <- original_variance
      block$pc1_score_variance <- suppressWarnings(as.numeric(block$score_variance))
      block$pc1_score_variance_ratio_to_original <- if (is.finite(original_variance) && original_variance > 0) block$pc1_score_variance / original_variance else NA_real_
      block$pc1_score_variance_percent_of_original <- 100 * block$pc1_score_variance_ratio_to_original
      block$dataset_group <- factor(as.character(block$dataset_group), levels = c("Original", "Lead", "Remainder"))
      block$evs_mode <- factor(as.character(block$evs_mode), levels = evs_modes)
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
    geom_text(aes(label = format_compact_number(pc1_score_variance, digits = 3)), vjust = -0.25, size = 1.65, color = "grey20") +
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
      caption = "Original, EVS leading edge, and EVS remainder are compared directly."
    ) +
    manuscript_theme() +
    theme(strip.text = element_text(size = base_theme_size - 0.6), axis.text.x = element_text(size = base_theme_size - 0.9))
}

pc1_variance_contribution_rows <- function(expr_matrix, sample_ids, scope_label) {
  x <- expr_matrix[, sample_ids, drop = FALSE]
  storage.mode(x) <- "numeric"
  keep <- rowSums(is.finite(x) & !is.na(x)) == ncol(x)
  keep <- keep & apply(x, 1, stats::var, na.rm = TRUE) > 0
  x <- x[keep, , drop = FALSE]
  if (nrow(x) < 2L || ncol(x) < 3L) return(data.frame())

  pca <- stats::prcomp(t(x), center = TRUE, scale. = FALSE)
  loading <- as.numeric(pca$rotation[, 1])
  pc1_score_variance <- as.numeric(pca$sdev[1]^2)

  data.frame(
    feature_id = rownames(pca$rotation),
    scope = scope_label,
    pc1_loading = loading,
    pc1_loading_abs = abs(loading),
    pc1_score_variance = pc1_score_variance,
    pc1_variance_contribution = loading^2 * pc1_score_variance,
    stringsAsFactors = FALSE
  )
}

evs_loading_distribution_rows <- function(full_count_matrix, evs_obj, coldata, comparison_name, mode) {
  extract_one <- function(mat, dataset_group) {
    source_matrix <- get_evs_matrix(mat, coldata, mode)
    trt_ids <- rownames(coldata)[coldata$condition == "trt"]
    ctl_ids <- rownames(coldata)[coldata$condition == "untrt"]
    all_ids <- rownames(coldata)

    bind_rows(
      pc1_variance_contribution_rows(source_matrix, all_ids, "All"),
      pc1_variance_contribution_rows(source_matrix, ctl_ids, "Control"),
      pc1_variance_contribution_rows(source_matrix, trt_ids, "Treatment")
    ) %>%
      mutate(comparison_name = comparison_name, evs_mode = mode, dataset_group = dataset_group)
  }

  out <- bind_rows(
    extract_one(full_count_matrix, "Original"),
    extract_one(evs_obj$leading_matrix, "Lead"),
    extract_one(evs_obj$remainder_matrix, "Remainder")
  )

  cutoff_rows <- data.frame(
    comparison_name = comparison_name,
    evs_mode = mode,
    scope = c("Treatment", "Control"),
    pc1_loading_cutoff = c(evs_obj$summary$treatment_loading_cutoff[1], evs_obj$summary$control_loading_cutoff[1]),
    pc1_score_variance = c(as.numeric(evs_obj$treatment_rank$pca$sdev[1]^2), as.numeric(evs_obj$control_rank$pca$sdev[1]^2)),
    stringsAsFactors = FALSE
  )
  cutoff_rows$pc1_variance_contribution_cutoff <- cutoff_rows$pc1_loading_cutoff^2 * cutoff_rows$pc1_score_variance
  cutoff_rows$dataset_group <- "Original"

  attr(out, "cutoffs") <- cutoff_rows
  out
}

plot_loading_histograms <- function(load_df, title_text) {
  if (nrow(load_df) == 0L) return(NULL)
  cutoff_df <- attr(load_df, "cutoffs")
  if (is.null(cutoff_df)) cutoff_df <- data.frame()

  load_df$evs_mode <- factor(load_df$evs_mode, levels = evs_modes)
  load_df$dataset_group <- factor(load_df$dataset_group, levels = c("Original", "Lead", "Remainder"))
  load_df$scope <- factor(load_df$scope, levels = c("All", "Control", "Treatment"))
  load_df$pc1_variance_contribution <- pmax(suppressWarnings(as.numeric(load_df$pc1_variance_contribution)), 0)

  if (nrow(cutoff_df) > 0L) {
    cutoff_df$evs_mode <- factor(cutoff_df$evs_mode, levels = evs_modes)
    cutoff_df$dataset_group <- factor(cutoff_df$dataset_group, levels = c("Original", "Lead", "Remainder"))
    cutoff_df$scope <- factor(cutoff_df$scope, levels = c("All", "Control", "Treatment"))
  }

  p <- ggplot(load_df, aes(pc1_variance_contribution, color = scope, fill = scope)) +
    geom_histogram(bins = 60, position = "identity", alpha = 0.32, linewidth = 0.22, boundary = 0) +
    facet_grid(evs_mode + dataset_group ~ comparison_name, scales = "free_y") +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10), labels = function(x) format_compact_number(x, digits = 3)) +
    scale_color_manual(values = c(All = "#7570B3", Control = control_color, Treatment = treatment_color), breaks = c("All", "Control", "Treatment"), drop = FALSE, name = NULL) +
    scale_fill_manual(values = c(All = "#7570B3", Control = control_color, Treatment = treatment_color), breaks = c("All", "Control", "Treatment"), drop = FALSE, name = NULL) +
    guides(color = guide_legend(override.aes = list(fill = c("#7570B3", control_color, treatment_color), alpha = 0.85, linewidth = 0.8)), fill = "none") +
    labs(
      title = title_text,
      x = "Per-feature contribution to PC1 score variance",
      y = "Feature frequency",
      caption = "Dashed vertical lines mark fixed top-N treatment/control EVS loading cutoffs in the original EVS input space only."
    ) +
    manuscript_theme() +
    theme(legend.position = "bottom", strip.text = element_text(size = base_theme_size - 0.8), axis.text = element_text(size = base_theme_size - 1.0))

  if (nrow(cutoff_df) > 0L) {
    p <- p + geom_vline(
      data = cutoff_df,
      aes(xintercept = pc1_variance_contribution_cutoff, color = scope),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.38,
      alpha = 0.85
    )
  }

  p
}

plot_discovery_counts <- function(summary_df) {
  analysis_levels <- c("Raw", "Lead_NormEVS", "Lead_RawEVS", "Rem_NormEVS", "Rem_RawEVS")
  analysis_labels <- c(
    Raw = "Original",
    Lead_NormEVS = "Lead NormEVS",
    Lead_RawEVS = "Lead RawEVS",
    Rem_NormEVS = "Remainder NormEVS",
    Rem_RawEVS = "Remainder RawEVS"
  )

  required <- c("comparison_name", "analysis_label", "dataset_key", "n_display_weak", "n_display_strong", "n_display_standard", "n_display_hbfss")
  missing <- setdiff(required, names(summary_df))
  if (length(missing) > 0L) stop("Discovery-count plot missing summary column(s): ", paste(missing, collapse = ", "), call. = FALSE)

  long_df <- bind_rows(lapply(seq_len(nrow(summary_df)), function(i) {
    row <- summary_df[i, , drop = FALSE]
    data.frame(
      comparison_name = row$comparison_name,
      analysis_label = row$analysis_label,
      dataset_key = row$dataset_key,
      Class = factor(c("Weak", "Strong", "Std", "HBFSS"), levels = c("Weak", "Strong", "Std", "HBFSS")),
      Count = as.numeric(c(row$n_display_weak, row$n_display_strong, row$n_display_standard, row$n_display_hbfss)),
      stringsAsFactors = FALSE
    )
  }))

  long_df$comparison_name <- factor(long_df$comparison_name, levels = comparison_table$comparison_name)
  long_df$analysis_label <- factor(long_df$analysis_label, levels = analysis_levels, labels = analysis_labels[analysis_levels])

  save_csv(rename_columns_existing(long_df, c(comparison_name = "comparison", analysis_label = "analysis", dataset_key = "dataset", Class = "class", Count = "count")), file.path(output_dir, "Counts_Long.csv"))

  max_count <- max(long_df$Count, na.rm = TRUE)
  count_breaks <- c(0, 1, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000)
  count_breaks <- count_breaks[count_breaks <= max_count * 1.15 | count_breaks <= 10]

  ggplot(long_df, aes(analysis_label, Count, fill = Class)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.72, color = "grey25", linewidth = 0.14) +
    facet_wrap(~ comparison_name, ncol = length(comparison_table$comparison_name)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = count_breaks, labels = function(x) format_compact_number(x, digits = 3), expand = expansion(mult = c(0.02, 0.15))) +
    scale_fill_manual(values = class_colors[c("Weak", "Strong", "Std", "HBFSS")], breaks = c("Weak", "Strong", "Std", "HBFSS"), labels = c(Weak = "Weak effect", Strong = "Strong effect", Std = "Standard-only", HBFSS = "HBFSS-only"), drop = FALSE, name = NULL) +
    labs(
      title = "Significant feature counts",
      x = NULL,
      y = "Final display-class count, pseudo-log scale",
      caption = "Counts are mutually exclusive final display classes matching the volcano legend."
    ) +
    manuscript_theme() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = base_theme_size - 1.5), strip.text = element_text(size = base_theme_size - 0.7))
}

# =============================================================================
# SIMULATION VALIDATION
# =============================================================================

simulation_lfc_label <- function(lfc_magnitude, simulation_profile) {
  profile <- as.character(simulation_profile)
  lfc <- suppressWarnings(as.numeric(lfc_magnitude))
  ifelse(profile == "weak_mixture", paste0("Weak mix\nLFC ", simulation_weak_lfc_min, "-", simulation_weak_lfc_max), paste0("Fixed LFC\n", lfc))
}

simulation_lfc_levels <- function() {
  c(paste0("Fixed LFC\n", simulation_lfc_magnitudes), paste0("Weak mix\nLFC ", simulation_weak_lfc_min, "-", simulation_weak_lfc_max))
}

simulate_sequence_counts <- function(n_features, n_samples, base_mean, disp_null, de_fraction, lfc_magnitude, disp_de, lfc_profile = c("fixed", "weak_mixture")) {
  lfc_profile <- match.arg(lfc_profile)
  n_de <- round(n_features * de_fraction)
  n_de <- max(2L, min(n_de, n_features - 2L))
  feature_id <- paste0("sim_feature_", seq_len(n_features))
  de_id <- seq_len(n_de)

  base_mu <- stats::rgamma(n_features, shape = 2.5, scale = base_mean / 2.5)
  base_mu <- pmax(base_mu, 2)
  dispersion <- rep(disp_null, n_features)
  dispersion[de_id] <- disp_de
  dispersion <- dispersion * exp(stats::rnorm(n_features, mean = 0, sd = 0.25))
  dispersion <- pmin(pmax(dispersion, 0.01), 1.50)

  de_direction <- sample(c(-1, 1), n_de, replace = TRUE)
  lfc_abs <- if (identical(lfc_profile, "weak_mixture")) {
    stats::runif(n_de, min = simulation_weak_lfc_min, max = simulation_weak_lfc_max)
  } else {
    rep(lfc_magnitude, n_de)
  }

  true_lfc <- rep(0, n_features)
  true_lfc[de_id] <- lfc_abs * de_direction
  control_mu <- base_mu
  treatment_mu <- base_mu * 2^true_lfc

  control_lib <- exp(stats::rnorm(n_samples, mean = 0, sd = 0.12))
  treatment_lib <- exp(stats::rnorm(n_samples, mean = 0, sd = 0.12))
  control_lib <- control_lib / exp(mean(log(control_lib)))
  treatment_lib <- treatment_lib / exp(mean(log(treatment_lib)))

  generate_group <- function(mu, lib_factor, dispersion_vector) {
    out <- matrix(0L, nrow = length(mu), ncol = length(lib_factor))
    for (j in seq_along(lib_factor)) {
      sample_mu <- pmax(mu * lib_factor[j], 1e-3)
      out[, j] <- stats::rnbinom(n = length(sample_mu), mu = sample_mu, size = 1 / dispersion_vector)
    }
    out
  }

  counts_mat <- cbind(generate_group(control_mu, control_lib, dispersion), generate_group(treatment_mu, treatment_lib, dispersion))
  storage.mode(counts_mat) <- "integer"
  rownames(counts_mat) <- feature_id
  colnames(counts_mat) <- c(paste0("ctrl_", seq_len(n_samples)), paste0("trt_", seq_len(n_samples)))

  truth <- data.frame(
    feature_id = feature_id,
    is_de = seq_len(n_features) %in% de_id,
    true_lfc = true_lfc,
    true_abs_lfc = abs(true_lfc),
    true_weak = seq_len(n_features) %in% de_id & abs(true_lfc) < lfc_boundary,
    true_strong = seq_len(n_features) %in% de_id & abs(true_lfc) >= lfc_boundary,
    base_mean = base_mu,
    dispersion = dispersion,
    simulation_profile = lfc_profile,
    stringsAsFactors = FALSE
  )

  list(counts = counts_mat, truth = truth)
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
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA_real_ else 2 * precision * recall / (precision + recall)
  data.frame(tp = tp, fp = fp, fn = fn, tn = tn, precision = precision, recall = recall, specificity = specificity, f1 = f1, fdr = fdr, discovery_count = sum(predicted), truth_count = sum(actual), stringsAsFactors = FALSE)
}

simulation_metric_row <- function(template, method, truth_target, predicted, actual, hc_p, hbfss_cutoff) {
  cbind(template, data.frame(method = method, truth_target = truth_target, hc_p_threshold = hc_p, hbfss_cutoff = hbfss_cutoff, hc_valid = is.finite(hc_p) & !is.na(hc_p), stringsAsFactors = FALSE), simulation_metrics(predicted, actual))
}

sequence_simulation_analysis <- function(sim_obj, null_inflation) {
  counts_mat <- sim_obj$counts
  n_samp <- as.integer(ncol(counts_mat) / 2L)

  if (!is.matrix(counts_mat) && !is.data.frame(counts_mat)) {
    stop("Simulation counts must be a matrix-like object.", call. = FALSE)
  }
  if (ncol(counts_mat) %% 2L != 0L) {
    stop("Simulation count matrix must contain paired control/treatment sample columns.", call. = FALSE)
  }
  if (n_samp < 2L) {
    stop("Simulation requires at least two samples per group.", call. = FALSE)
  }

  counts_mat <- as_integer_count_matrix(counts_mat, "simulation DESeq2 input")

  coldata <- data.frame(
    condition = factor(
      c(rep("untrt", n_samp), rep("trt", n_samp)),
      levels = c("untrt", "trt")
    ),
    row.names = colnames(counts_mat)
  )

  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = coldata,
    design = ~ condition
  )

  dds <- dds[rowSums(counts(dds)) > 0, ]
  if (nrow(dds) < 5L) {
    stop("Simulation DESeq2 object has fewer than five nonzero features after filtering.", call. = FALSE)
  }

  dds <- DESeq(dds, betaPrior = FALSE, quiet = TRUE)
  coef_name <- condition_coef_name(dds)

  standard <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    alpha = alpha_standard
  )

  greater_abs <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "greaterAbs",
    alpha = alpha_strong
  )

  less_abs <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "lessAbs",
    alpha = alpha_weak
  )

  df <- as.data.frame(standard)
  df$feature_id <- rownames(df)

  if (isTRUE(simulation_use_apeglm_shrinkage)) {
    shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
    shrink_df <- data.frame(
      feature_id = rownames(shrunk),
      lfc_shrunk = as.data.frame(shrunk)$log2FoldChange,
      stringsAsFactors = FALSE
    )
  } else {
    shrink_df <- data.frame(
      feature_id = rownames(df),
      lfc_shrunk = df$log2FoldChange,
      stringsAsFactors = FALSE
    )
  }

  df <- df %>%
    left_join(shrink_df, by = "feature_id") %>%
    left_join(
      data.frame(
        feature_id = rownames(greater_abs),
        greaterAbs_padj = greater_abs$padj,
        stringsAsFactors = FALSE
      ),
      by = "feature_id"
    ) %>%
    left_join(
      data.frame(
        feature_id = rownames(less_abs),
        lessAbs_padj = less_abs$padj,
        stringsAsFactors = FALSE
      ),
      by = "feature_id"
    ) %>%
    left_join(sim_obj$truth, by = "feature_id")

  df$is_de[is.na(df$is_de)] <- FALSE
  df$true_weak[is.na(df$true_weak)] <- FALSE
  df$true_strong[is.na(df$true_strong)] <- FALSE

  wald_for_empirical_null <- inflate_null_wald_statistics(
    wald = df$stat,
    is_de = df$is_de,
    inflation_fraction = null_inflation
  )

  empirical <- tryCatch(
    fit_empirical_null(wald_for_empirical_null, "simulation"),
    error = function(e) {
      log_message("Simulation empirical-null fallback used: ", conditionMessage(e))

      p_raw <- rep(NA_real_, length(wald_for_empirical_null))
      ok <- is.finite(wald_for_empirical_null) & !is.na(wald_for_empirical_null)
      p_raw[ok] <- 2 * stats::pnorm(-abs(wald_for_empirical_null[ok]))
      p_raw <- clip_probability(p_raw)

      empirical_bh <- rep(NA_real_, length(p_raw))
      ok_p <- is.finite(p_raw) & !is.na(p_raw)
      empirical_bh[ok_p] <- p.adjust(p_raw[ok_p], method = "BH")

      list(
        empirical_p = p_raw,
        empirical_q = empirical_bh,
        empirical_lfdr = rep(NA_real_, length(p_raw)),
        empirical_bh = empirical_bh,
        empirical_null_method = "theoretical_normal_simulation_fallback"
      )
    }
  )

  df$empirical_p <- empirical$empirical_p
  df$empirical_bh <- empirical$empirical_bh
  df$empirical_null_method <- empirical$empirical_null_method
  df$neglog10_empirical_p_calc <- safe_neglog10(
    df$empirical_p,
    floor_value = calculation_p_floor
  )

  hc_p <- hc_threshold(df$empirical_p)
  # Mirrors the manuscript HBFSS cutoff coupling to the LFC boundary.
  hbfss_cutoff <- if (is.na(hc_p)) NA_real_ else -log10(hc_p) * lfc_boundary

  df$hc_pass <- !is.na(hc_p) &
    !is.na(df$empirical_p) &
    is.finite(df$empirical_p) &
    df$empirical_p <= hc_p

  df$HBFSS <- abs(df$lfc_shrunk) * df$neglog10_empirical_p_calc

  df$standard_sig <- !is.na(df$padj) &
    df$padj < alpha_standard &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= lfc_boundary

  df$greaterAbs_sig <- !is.na(df$greaterAbs_padj) &
    df$greaterAbs_padj < alpha_strong &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= lfc_boundary

  df$lessAbs_sig <- !is.na(df$lessAbs_padj) &
    df$lessAbs_padj < alpha_weak &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) < lfc_boundary

  df$empirical_bh_sig <- !is.na(df$empirical_bh) &
    df$empirical_bh < alpha_standard &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= lfc_boundary

  df$hbfss_raw_sig <- !is.na(hbfss_cutoff) &
    !is.na(df$HBFSS) &
    is.finite(df$HBFSS) &
    df$HBFSS >= hbfss_cutoff &
    df$hc_pass

  df$hbfss_total_sig <- df$standard_sig | df$hbfss_raw_sig
  df$weak_region_hbfss_sig <- df$lessAbs_sig & df$hbfss_raw_sig

  list(
    results = df,
    hc_p = hc_p,
    hbfss_cutoff = hbfss_cutoff
  )
}

build_simulation_metric_rows <- function(out, template) {
  df <- out$results
  bind_rows(
    simulation_metric_row(template, "DESeq2_BH", "all_de", df$standard_sig, df$is_de, out$hc_p, out$hbfss_cutoff),
    simulation_metric_row(template, "Empirical_BH", "all_de", df$empirical_bh_sig, df$is_de, out$hc_p, out$hbfss_cutoff),
    simulation_metric_row(template, "GreaterAbs", "strong_de", df$greaterAbs_sig, df$true_strong, out$hc_p, out$hbfss_cutoff),
    simulation_metric_row(template, "LessAbs", "weak_de", df$lessAbs_sig, df$true_weak, out$hc_p, out$hbfss_cutoff),
    simulation_metric_row(template, "HBFSS_total", "all_de", df$hbfss_total_sig, df$is_de, out$hc_p, out$hbfss_cutoff),
    simulation_metric_row(template, "HBFSS_raw", "all_de", df$hbfss_raw_sig, df$is_de, out$hc_p, out$hbfss_cutoff),
    simulation_metric_row(template, "HBFSS_weak_region", "weak_de", df$weak_region_hbfss_sig, df$true_weak, out$hc_p, out$hbfss_cutoff)
  )
}

plot_simulation_metric_boxplot <- function(metric_long, metric, title, y_label, methods, target, alpha_line = FALSE) {
  plot_df <- metric_long[metric_long$method %in% methods & metric_long$truth_target == target, , drop = FALSE]
  if (nrow(plot_df) == 0L) return(NULL)
  plot_df$method <- factor(plot_df$method, levels = methods)
  method_colors <- c(DESeq2_BH = "#999999", Empirical_BH = treatment_color, GreaterAbs = class_colors[["Strong"]], LessAbs = class_colors[["Weak"]], HBFSS_total = class_colors[["HBFSS"]], HBFSS_raw = "#54278F", HBFSS_weak_region = class_colors[["Weak"]])

  p <- ggplot(plot_df, aes(method, .data[[metric]], fill = method)) +
    geom_boxplot(outlier.size = 0.35, width = 0.62, linewidth = 0.22, na.rm = TRUE) +
    facet_grid(inflation_label + de_label ~ lfc_label) +
    scale_fill_manual(values = method_colors[methods], breaks = methods, drop = FALSE, name = NULL) +
    labs(title = title, x = NULL, y = y_label) +
    manuscript_theme() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(size = base_theme_size - 1.0), legend.position = "bottom")
  if (alpha_line) p <- p + geom_hline(yintercept = alpha_standard, linetype = "dashed", linewidth = 0.30, color = "grey35")
  p
}

run_sequence_simulation_validation <- function() {
  set.seed(simulation_seed)
  dir.create(simulation_dir, recursive = TRUE, showWarnings = FALSE)
  old_files <- list.files(simulation_dir, full.names = TRUE)
  if (length(old_files) > 0L) unlink(old_files, recursive = TRUE, force = TRUE)

  # Default smoke test: one fixed-effect run plus one weak-mixture run.
  fixed_grid <- expand.grid(de_fraction = simulation_de_fractions, lfc_magnitude = simulation_lfc_magnitudes, null_inflation = simulation_null_inflation, simulation_profile = "fixed", stringsAsFactors = FALSE)
  weak_grid <- expand.grid(de_fraction = simulation_de_fractions, lfc_magnitude = NA_real_, null_inflation = simulation_null_inflation, simulation_profile = "weak_mixture", stringsAsFactors = FALSE)
  sim_grid <- bind_rows(fixed_grid, weak_grid)
  total_runs <- nrow(sim_grid) * simulation_n_reps
  log_message("Running simulation validation: ", total_runs, " planned replicates")
  progress_every <- max(1L, suppressWarnings(as.integer(simulation_progress_every)), na.rm = TRUE)
  checkpoint_every <- max(1L, suppressWarnings(as.integer(simulation_checkpoint_every)), na.rm = TRUE)
  log_message("Simulation config: features=", simulation_n_features, ", samples/group=", simulation_n_samples_per_group, ", reps/grid=", simulation_n_reps)
  if (simulation_n_reps < simulation_manuscript_min_reps) {
    log_message("Simulation mode: smoke test only; increase simulation_n_reps to at least ", simulation_manuscript_min_reps, " for manuscript validation figures.")
  }
  if (!isTRUE(simulation_use_apeglm_shrinkage)) {
    log_message("Simulation uses unshrunk DESeq2 log2 fold changes for speed; set simulation_use_apeglm_shrinkage <- TRUE for exact classification parity with the main manuscript analysis.")
  }
  simulation_start <- Sys.time()

  metric_rows <- list()
  failure_rows <- list()
  feature_rows <- list()
  run_index <- 0L

  for (grid_i in seq_len(nrow(sim_grid))) {
    grid_row <- sim_grid[grid_i, , drop = FALSE]
    for (rep_i in seq_len(simulation_n_reps)) {
      run_index <- run_index + 1L
      if (run_index %% progress_every == 0L || run_index == 1L || run_index == total_runs) {
        elapsed_min <- as.numeric(difftime(Sys.time(), simulation_start, units = "mins"))
        eta_min <- if (run_index > 0L) elapsed_min * (total_runs - run_index) / run_index else NA_real_
        log_message(
          "Simulation progress: ", run_index, "/", total_runs,
          " | elapsed=", signif(elapsed_min, 3), " min",
          " | ETA=", signif(eta_min, 3), " min"
        )
      }

      template <- data.frame(de_fraction = grid_row$de_fraction, lfc_magnitude = grid_row$lfc_magnitude, simulation_profile = grid_row$simulation_profile, null_inflation = grid_row$null_inflation, replicate = rep_i, stringsAsFactors = FALSE)
      template$lfc_label <- simulation_lfc_label(template$lfc_magnitude, template$simulation_profile)
      template$de_label <- paste0(template$de_fraction * 100, "% DE")
      template$inflation_label <- paste0(template$null_inflation * 100, "% null inflation")

      sim_result <- tryCatch({
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
        list(ok = TRUE, out = sequence_simulation_analysis(sim_obj, null_inflation = grid_row$null_inflation), error_message = NA_character_)
      }, error = function(e) {
        list(ok = FALSE, out = NULL, error_message = conditionMessage(e))
      })

      if (!isTRUE(sim_result$ok)) {
        failure_rows[[length(failure_rows) + 1L]] <- cbind(template, data.frame(error_message = sim_result$error_message, stringsAsFactors = FALSE))
        next
      }

      out <- sim_result$out
      metric_rows[[length(metric_rows) + 1L]] <- build_simulation_metric_rows(out, template)
      if (length(metric_rows) > 0L && (length(metric_rows) %% checkpoint_every == 0L || run_index == total_runs)) {
        save_csv(bind_rows(metric_rows), file.path(simulation_dir, "Simulation_RunMetrics_Checkpoint.csv"))
      }

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
  if (nrow(metric_long) == 0L) stop("Simulation validation produced no successful replicates. See Simulation_Failures.csv.", call. = FALSE)
  metric_long$lfc_label <- factor(metric_long$lfc_label, levels = simulation_lfc_levels())
  metric_long$method <- factor(metric_long$method, levels = c("DESeq2_BH", "Empirical_BH", "GreaterAbs", "LessAbs", "HBFSS_raw", "HBFSS_total", "HBFSS_weak_region"))

  metric_summary <- metric_long %>%
    group_by(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label, method, truth_target) %>%
    summarise(
      n_replicates = n(),
      precision_mean = mean_finite(precision),
      recall_mean = mean_finite(recall),
      f1_mean = mean_finite(f1),
      fdr_mean = mean_finite(fdr),
      discovery_count_mean = mean_finite(discovery_count),
      truth_count_mean = mean_finite(truth_count),
      hc_valid_fraction = mean(as.numeric(hc_valid), na.rm = TRUE),
      hbfss_cutoff_median = median_finite(hbfss_cutoff),
      .groups = "drop"
    )

  save_csv(metric_long, file.path(simulation_dir, "Simulation_RunMetrics_Long.csv"))
  save_csv(metric_summary, file.path(simulation_dir, "Simulation_MethodSummary.csv"))

  if (isTRUE(export_simulation_feature_results) && length(feature_rows) > 0L) save_csv(bind_rows(feature_rows), file.path(simulation_dir, "Simulation_FeatureResults.csv"))

  safe_save_figure(plot_simulation_metric_boxplot(metric_long, "f1", "Simulation F1 score by method", "F1", c("DESeq2_BH", "Empirical_BH", "HBFSS_total", "HBFSS_raw"), "all_de"), file.path(simulation_dir, "Simulation_F1_Boxplot.png"), 11.8, 9.2)
  safe_save_figure(plot_simulation_metric_boxplot(metric_long, "fdr", "Simulation observed FDR by method", "FDR", c("DESeq2_BH", "Empirical_BH", "HBFSS_total", "HBFSS_raw"), "all_de", alpha_line = TRUE), file.path(simulation_dir, "Simulation_FDR_Boxplot.png"), 11.8, 9.2)

  writeLines(capture.output(sessionInfo()), file.path(simulation_dir, "SessionInfo_Simulation.txt"))
  log_message("Simulation successful metric blocks: ", length(metric_rows))
  log_message("Simulation failures: ", nrow(failures))
  invisible(metric_summary)
}

# =============================================================================
# METHODS, MANIFEST, AND GIT
# =============================================================================

write_package_versions <- function() {
  versions <- data.frame(
    package = required_packages,
    version = vapply(required_packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1)),
    stringsAsFactors = FALSE
  )
  save_csv(versions, file.path(output_dir, "PackageVersions.csv"))
  invisible(versions)
}

write_sequence_methods <- function() {
  methods_lines <- c(
    "# SEQUENCE methods summary",
    "",
    "## Run parameters",
    paste0("- alpha_standard: ", alpha_standard),
    paste0("- alpha_strong: ", alpha_strong),
    paste0("- alpha_weak: ", alpha_weak),
    paste0("- LFC boundary: ", lfc_boundary),
    paste0("- EVS cutoff: fixed top-N union; N = ", evs_top_n),
    paste0("- HC invalid lower threshold: ", hc_invalid_below),
    paste0("- HC invalid upper threshold: ", hc_invalid_at_or_above),
    paste0("- Strict empirical null: ", strict_empirical_null),
    paste0("- Comparisons: ", paste(comparison_table$comparison_name, collapse = ", ")),
    paste0("- EVS modes: ", paste(evs_modes, collapse = ", ")),
    "",
    "## Input and preprocessing",
    "Raw read-count matrices are imported for each comparison. Raw integer counts are retained for DESeq2. NormEVS uses DESeq2 VST expression; RawEVS uses log2(raw counts + 1).",
    "",
    "## Eigenvector splitting",
    paste0("PC1 is calculated separately in treatment and control samples. Features are ranked by absolute PC1 loading. The leading edge is the union of the top-", evs_top_n, " treatment-ranked and top-", evs_top_n, " control-ranked features. The remainder contains all other features."),
    "",
    "## Differential expression and HBFSS",
    paste0("DESeq2 design is ~ condition with untrt as the reference. Standard effects use BH padj < ", alpha_standard, " and |apeglm-shrunken LFC| >= ", lfc_boundary, ". Strong effects use greaterAbs. Weak effects use lessAbs, |shrunken LFC| < boundary, and HBFSS raw significance. All three alpha thresholds are intentionally matched at 0.10 for direct comparison. HBFSS = |shrunken LFC| x -log10(empirical p). The HBFSS cutoff is -log10(HC p threshold) x LFC boundary."),
    "",
    "## Figure export",
    "Figures are exported as PNG and PDF. Panel grobs are rendered through explicit grid devices so arranged multi-panel figures are saved reliably. Combined panel figures suppress repeated per-plot captions, gene labels, and threshold text while retaining threshold lines and one shared legend.",
    "",
    "## Simulation validation",
    paste0("Simulation validation is enabled as a smoke-test run by default: n_features=", simulation_n_features, ", n_samples_per_group=", simulation_n_samples_per_group, ", n_reps=", simulation_n_reps, ", simulation_use_apeglm_shrinkage=", simulation_use_apeglm_shrinkage, ". This smoke-test proves execution. For manuscript validation, set n_reps >= ", simulation_manuscript_min_reps, " and consider simulation_use_apeglm_shrinkage <- TRUE to match the main analysis exactly.")
  )
  writeLines(methods_lines, file.path(output_dir, "METHODS_SEQUENCE_PIPELINE.md"))
}

write_manifest <- function() {
  exported_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  exported_files <- exported_files[file.info(exported_files)$isdir %in% FALSE]
  root <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  manifest <- data.frame(
    file = sub(paste0("^", root, "/?"), "", normalizePath(exported_files, winslash = "/", mustWork = FALSE)),
    size_bytes = file.info(exported_files)$size,
    stringsAsFactors = FALSE
  )
  save_csv(manifest, file.path(output_dir, "Manifest.csv"))
  manifest
}

path_relative_to_repo <- function(path) {
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  abs_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  if (startsWith(abs_path, prefix)) return(substr(abs_path, nchar(prefix) + 1L, nchar(abs_path)))
  abs_path
}

git_command <- function(args, allow_failure = FALSE, echo_output = TRUE) {
  if (Sys.which("git") == "") stop("Git executable was not found on PATH.", call. = FALSE)
  cmd_display <- paste("git", "-C", repo_root, paste(args, collapse = " "))
  log_message("$ ", cmd_display)
  out <- suppressWarnings(system2("git", args = c("-C", repo_root, args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  status <- as.integer(status)
  if (isTRUE(echo_output) && length(out) > 0L) log_message(paste(out, collapse = "\n"))
  if (!isTRUE(allow_failure) && status != 0L) {
    stop("Git command failed with status ", status, ": ", cmd_display, "\n", paste(out, collapse = "\n"), call. = FALSE)
  }
  list(status = status, output = out)
}

git_output_first_line <- function(args, allow_failure = FALSE) {
  res <- git_command(args, allow_failure = allow_failure, echo_output = FALSE)
  if (res$status != 0L || length(res$output) == 0L) return(NA_character_)
  trimws(res$output[1])
}

git_has_staged_changes <- function() {
  res <- git_command(c("diff", "--cached", "--quiet"), allow_failure = TRUE, echo_output = FALSE)
  if (res$status == 0L) return(FALSE)
  if (res$status == 1L) return(TRUE)
  stop("Unable to inspect staged Git changes.", call. = FALSE)
}

ensure_gitignore_patterns <- function(patterns) {
  patterns <- unique(patterns[nzchar(patterns)])
  if (length(patterns) == 0L) return(invisible(FALSE))
  ignore_path <- file.path(repo_root, ".gitignore")
  existing <- if (file.exists(ignore_path)) readLines(ignore_path, warn = FALSE) else character(0)
  missing <- setdiff(patterns, existing)
  if (length(missing) > 0L) {
    cat(paste0(missing, collapse = "\n"), "\n", file = ignore_path, append = TRUE, sep = "")
    log_message("Updated .gitignore with GitHub-safe export exclusions: ", paste(missing, collapse = ", "))
  }
  invisible(length(missing) > 0L)
}

remove_skipped_large_export_files <- function() {
  if (!isTRUE(remove_skipped_large_exports_from_disk)) return(invisible(FALSE))
  targets <- c(
    file.path(output_dir, "Summary_PC1VarianceDistributions.csv"),
    file.path(simulation_dir, "Simulation_FeatureResults.csv")
  )
  removed <- character(0)
  for (target in targets) {
    if (file.exists(target)) {
      unlink(target, force = TRUE)
      removed <- c(removed, path_relative_to_repo(target))
    }
  }
  if (length(removed) > 0L) log_message("Removed large skipped export file(s) from disk: ", paste(removed, collapse = ", "))
  invisible(length(removed) > 0L)
}

git_remove_cached_excluded_exports <- function() {
  patterns <- unique(git_exclude_large_export_patterns)
  patterns <- patterns[nzchar(patterns)]
  if (length(patterns) == 0L) return(invisible(FALSE))
  git_command(c("rm", "-r", "--cached", "--ignore-unmatch", "--", patterns), allow_failure = TRUE)
  invisible(TRUE)
}

git_staged_paths <- function() {
  res <- git_command(c("diff", "--cached", "--name-only"), allow_failure = TRUE, echo_output = FALSE)
  if (res$status != 0L || length(res$output) == 0L) return(character(0))
  unique(trimws(res$output[nzchar(res$output)]))
}

git_unstage_oversized_files <- function(max_bytes = git_max_file_size_bytes) {
  staged <- git_staged_paths()
  if (length(staged) == 0L) return(invisible(character(0)))
  oversized <- character(0)
  for (rel in staged) {
    abs_path <- file.path(repo_root, rel)
    if (file.exists(abs_path)) {
      sz <- suppressWarnings(file.info(abs_path)$size)
      if (is.finite(sz) && !is.na(sz) && sz > max_bytes) oversized <- c(oversized, rel)
    }
  }
  oversized <- unique(oversized)
  if (length(oversized) > 0L) {
    git_command(c("reset", "-q", "HEAD", "--", oversized), allow_failure = FALSE)
    log_message("Unstaged oversized file(s) to protect GitHub push: ", paste(oversized, collapse = ", "))
  }
  invisible(oversized)
}

git_unpushed_large_blobs <- function(upstream_ref, max_bytes = git_max_file_size_bytes) {
  if (is.na(upstream_ref) || !nzchar(upstream_ref)) return(data.frame())
  rev_range <- paste0(upstream_ref, "..HEAD")
  revs <- git_command(c("rev-list", "--objects", rev_range), allow_failure = TRUE, echo_output = FALSE)
  if (revs$status != 0L || length(revs$output) == 0L) return(data.frame())

  rows <- list()
  for (line in revs$output) {
    line <- trimws(line)
    if (!nzchar(line)) next
    parts <- strsplit(line, " ", fixed = TRUE)[[1]]
    sha <- parts[1]
    rel <- if (length(parts) > 1L) paste(parts[-1], collapse = " ") else NA_character_
    if (is.na(rel) || !nzchar(rel)) next
    sz <- git_command(c("cat-file", "-s", sha), allow_failure = TRUE, echo_output = FALSE)
    if (sz$status != 0L || length(sz$output) == 0L) next
    bytes <- suppressWarnings(as.numeric(trimws(sz$output[1])))
    if (is.finite(bytes) && !is.na(bytes) && bytes > max_bytes) {
      rows[[length(rows) + 1L]] <- data.frame(object = sha, file = rel, size_bytes = bytes, stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0L) return(data.frame())
  bind_rows(rows)
}

git_rewrite_unpushed_large_export_history <- function(upstream_ref) {
  if (!isTRUE(git_rewrite_unpushed_large_export_commits)) return(invisible(FALSE))
  large <- git_unpushed_large_blobs(upstream_ref)
  if (nrow(large) == 0L) return(invisible(FALSE))

  export_large <- large[grepl("^exports/", large$file), , drop = FALSE]
  if (nrow(export_large) == 0L) {
    stop(
      "Unpushed commit history contains oversized non-export file(s), so automatic rewrite was not attempted: ",
      paste(large$file, collapse = ", "),
      call. = FALSE
    )
  }

  log_message(
    "Oversized export blob(s) found in unpushed commit history: ",
    paste(paste0(export_large$file, " (", round(export_large$size_bytes / 1024^2, 1), " MB)"), collapse = ", ")
  )
  log_message("Rewriting unpushed export commits with git reset --mixed ", upstream_ref, " so the large blobs are removed before push. Worktree files are preserved and will be restaged safely.")
  git_command(c("reset", "--mixed", upstream_ref), allow_failure = FALSE)
  git_remove_cached_excluded_exports()
  invisible(TRUE)
}

git_commit_and_push <- function() {
  if (!isTRUE(git_push_after_success)) {
    log_message("Git push disabled: git_push_after_success is FALSE.")
    return(invisible(FALSE))
  }

  git_root <- git_output_first_line(c("rev-parse", "--show-toplevel"), allow_failure = TRUE)
  if (is.na(git_root) || !nzchar(git_root)) stop("Git push requested, but this run is not inside a Git repository.", call. = FALSE)
  git_root <- normalizePath(git_root, winslash = "/", mustWork = TRUE)
  expected_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  if (!identical(git_root, expected_root)) stop("Git root mismatch. repo_root is ", expected_root, " but git root is ", git_root, call. = FALSE)

  ensure_gitignore_patterns(git_exclude_large_export_patterns)
  remove_skipped_large_export_files()

  if (isTRUE(git_fail_if_preexisting_staged_changes) && git_has_staged_changes()) {
    stop("Git index already contains staged changes before the SEQUENCE push step. Commit or unstage them first.", call. = FALSE)
  }

  remote_url <- git_output_first_line(c("remote", "get-url", git_remote_name), allow_failure = TRUE)
  if (is.na(remote_url) || !nzchar(remote_url)) stop("Git remote '", git_remote_name, "' is not configured.", call. = FALSE)

  current_branch <- git_output_first_line(c("rev-parse", "--abbrev-ref", "HEAD"), allow_failure = FALSE)
  if (identical(current_branch, "HEAD") || is.na(current_branch) || !nzchar(current_branch)) stop("Git is in detached HEAD state.", call. = FALSE)
  target_branch <- if (!is.na(git_branch_name) && nzchar(git_branch_name)) git_branch_name else current_branch

  git_command(c("fetch", git_remote_name), allow_failure = FALSE)
  upstream <- git_output_first_line(c("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"), allow_failure = TRUE)
  if (is.na(upstream) || !nzchar(upstream)) {
    remote_branch_ref <- paste0(git_remote_name, "/", target_branch)
    remote_branch_sha <- git_output_first_line(c("rev-parse", "--verify", remote_branch_ref), allow_failure = TRUE)
    if (!is.na(remote_branch_sha) && nzchar(remote_branch_sha)) upstream <- remote_branch_ref
  }
  if (isTRUE(git_pull_rebase_before_push) && !is.na(upstream) && nzchar(upstream)) git_command(c("pull", "--rebase", "--autostash", git_remote_name, target_branch), allow_failure = FALSE)

  upstream_after_pull <- git_output_first_line(c("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"), allow_failure = TRUE)
  if (is.na(upstream_after_pull) || !nzchar(upstream_after_pull)) upstream_after_pull <- upstream
  if (!is.na(upstream_after_pull) && nzchar(upstream_after_pull)) git_rewrite_unpushed_large_export_history(upstream_after_pull)

  git_remove_cached_excluded_exports()

  stage_paths <- c(".gitignore", output_dir)
  if (isTRUE(git_stage_pipeline_script)) {
    sp <- script_path()
    if (!is.na(sp) && file.exists(sp)) stage_paths <- c(sp, stage_paths)
  }
  stage_paths <- unique(vapply(stage_paths, path_relative_to_repo, character(1)))
  git_command(c("add", "--", stage_paths), allow_failure = FALSE)

  git_remove_cached_excluded_exports()
  oversized <- git_unstage_oversized_files()
  if (length(oversized) > 0L) {
    ensure_gitignore_patterns(oversized)
    git_command(c("add", "--", ".gitignore"), allow_failure = FALSE)
  }

  if (!git_has_staged_changes()) {
    log_message("No changed SEQUENCE files to commit.")
    if (isTRUE(git_allow_no_change_success)) return(invisible(FALSE))
    stop("No changed files to commit and git_allow_no_change_success is FALSE.", call. = FALSE)
  }

  git_command(c("status", "--short"), allow_failure = FALSE)
  git_command(c("commit", "-m", git_commit_message), allow_failure = FALSE)

  final_upstream <- git_output_first_line(c("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"), allow_failure = TRUE)
  if (is.na(final_upstream) || !nzchar(final_upstream)) {
    remote_branch_ref <- paste0(git_remote_name, "/", target_branch)
    remote_branch_sha <- git_output_first_line(c("rev-parse", "--verify", remote_branch_ref), allow_failure = TRUE)
    if (!is.na(remote_branch_sha) && nzchar(remote_branch_sha)) final_upstream <- remote_branch_ref
  }
  if (!is.na(final_upstream) && nzchar(final_upstream)) {
    large_after_commit <- git_unpushed_large_blobs(final_upstream)
    if (nrow(large_after_commit) > 0L) {
      stop(
        "Refusing to push because unpushed history still contains oversized file(s): ",
        paste(paste0(large_after_commit$file, " (", round(large_after_commit$size_bytes / 1024^2, 1), " MB)"), collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (!is.na(final_upstream) && nzchar(final_upstream) && identical(target_branch, current_branch)) {
    git_command(c("push"), allow_failure = FALSE)
  } else {
    git_command(c("push", "-u", git_remote_name, paste0("HEAD:", target_branch)), allow_failure = FALSE)
  }
  log_message("GitHub push complete.")
  invisible(TRUE)
}

# =============================================================================
# PIPELINE EXECUTION
# =============================================================================

cleanup_run_artifacts()
cat("", file = log_file)

analysis_store <- list()
evs_store <- list()
comparison_store <- list()
summary_rows <- list()
failure_rows <- list()
figure_status_rows <- list()

run_store_save <- function(count_matrix, coldata, comparison_name, analysis_label, dataset_key) {
  fit <- run_deseq2_hbfss(count_matrix, coldata, comparison_name, analysis_label, dataset_key)
  save_csv(fit$results, analysis_table_path(comparison_name, analysis_label, "Results"))
  save_csv(concise_analysis_summary(fit$summary), analysis_table_path(comparison_name, analysis_label, "Summary"))
  key <- paste(comparison_name, analysis_label, sep = "__")
  list(key = key, fit = fit, summary = fit$summary)
}


for (i in seq_len(nrow(comparison_table))) {
  comparison <- prepare_comparison(i)
  comparison_name <- comparison$comparison_name
  comparison_store[[comparison_name]] <- comparison

  log_message("=====================================================")
  log_message("Running comparison: ", comparison_name)
  log_message("=====================================================")

  comparison_failure <- tryCatch({
    stored <- run_store_save(comparison$count_matrix, comparison$coldata, comparison_name, "Raw", "Raw")
    analysis_store[[stored$key]] <- stored$fit
    summary_rows[[stored$key]] <- stored$summary

    for (mode in evs_modes) {
      evs <- build_evs_split(comparison$count_matrix, comparison$coldata, comparison_name, mode)
      evs_key <- paste(comparison_name, mode, sep = "__")
      evs_store[[evs_key]] <- evs

      save_csv(evs$joint_rank, analysis_table_path(comparison_name, mode, "EVS_Rank_Table"))
      save_csv(concise_evs_summary(evs$summary), analysis_table_path(comparison_name, mode, "EVS_Summary"))

      stored <- run_store_save(evs$leading_matrix, comparison$coldata, comparison_name, paste0("Lead_", mode), "Lead")
      analysis_store[[stored$key]] <- stored$fit
      summary_rows[[stored$key]] <- stored$summary

      stored <- run_store_save(evs$remainder_matrix, comparison$coldata, comparison_name, paste0("Rem_", mode), "Rem")
      analysis_store[[stored$key]] <- stored$fit
      summary_rows[[stored$key]] <- stored$summary
    }
    NULL
  }, error = function(e) {
    data.frame(comparison_name = comparison_name, error_message = conditionMessage(e), stringsAsFactors = FALSE)
  })

  if (!is.null(comparison_failure)) {
    failure_rows[[comparison_name]] <- comparison_failure
    log_message("FAILED: ", comparison_name, ": ", comparison_failure$error_message[1])
  }
}

summary_df <- if (length(summary_rows) > 0L) bind_rows(summary_rows) else data.frame()
if (nrow(summary_df) > 0L) save_csv(concise_analysis_summary(summary_df), file.path(output_dir, "Summary_Overall.csv"))

if (length(failure_rows) > 0L) {
  save_csv(bind_rows(failure_rows), file.path(output_dir, "Failures.csv"))
  fail_pipeline("One or more comparisons failed. See Failures.csv. Figure export was not attempted because the analysis run is incomplete.")
}
if (nrow(summary_df) == 0L) fail_pipeline("No successful analyses were completed.")

# =============================================================================
# FIGURE EXPORT
# =============================================================================

comparison_order <- comparison_table$comparison_name

if (isTRUE(save_individual_figures)) {
  for (comparison_name in comparison_order) {
    for (analysis_label in c("Raw", "Lead_NormEVS", "Lead_RawEVS", "Rem_NormEVS", "Rem_RawEVS")) {
      key <- paste(comparison_name, analysis_label, sep = "__")
      if (!is.null(analysis_store[[key]])) {
        plot_title <- paste0(comparison_name, "\n", analysis_label_pretty(analysis_label))
        out_name <- paste0("Volcano_", comparison_name, "_", analysis_label, ".png")
        figure_status_rows[[length(figure_status_rows) + 1L]] <- safe_save_figure(plot_volcano(analysis_store[[key]]$results, plot_title), file.path(figure_dir, out_name), 6.4, 5.3)
      }
    }
  }
}

raw_plots <- with_panel_mode({
  plots <- list()
  for (comparison_name in comparison_order) {
    key <- paste(comparison_name, "Raw", sep = "__")
    if (!is.null(analysis_store[[key]])) plots[[comparison_name]] <- plot_volcano(analysis_store[[key]]$results, comparison_name)
  }
  plots
})
raw_panel <- arrange_with_one_legend(raw_plots, "Volcano: Original", ncol = length(comparison_order))
figure_status_rows[[length(figure_status_rows) + 1L]] <- safe_save_figure(raw_panel, file.path(figure_dir, "Volcano_Original.png"), 18.0, 5.8)

for (dataset_prefix in c("Lead", "Rem")) {
  plots <- with_panel_mode({
    tmp <- list()
    for (mode in evs_modes) {
      for (comparison_name in comparison_order) {
        analysis_label <- paste0(dataset_prefix, "_", mode)
        key <- paste(comparison_name, analysis_label, sep = "__")
        if (!is.null(analysis_store[[key]])) tmp[[paste(comparison_name, mode, sep = "_")]] <- plot_volcano(analysis_store[[key]]$results, paste0(comparison_name, "\n", mode))
      }
    }
    tmp
  })
  panel_title <- if (dataset_prefix == "Lead") "Volcano: Leading Edge" else "Volcano: Remainder"
  panel_file <- if (dataset_prefix == "Lead") "Volcano_Lead.png" else "Volcano_Remainder.png"
  panel <- arrange_with_one_legend(plots, panel_title, ncol = length(comparison_order))
  figure_status_rows[[length(figure_status_rows) + 1L]] <- safe_save_figure(panel, file.path(figure_dir, panel_file), 18.0, 9.6)
}

figure_status_rows[[length(figure_status_rows) + 1L]] <- safe_save_figure(plot_discovery_counts(summary_df), file.path(figure_dir, "Counts.png"), 18.0, 6.2)

if (isTRUE(save_individual_figures)) {
  for (comparison_name in comparison_order) {
    for (analysis_label in c("Raw", "Lead_NormEVS", "Lead_RawEVS", "Rem_NormEVS", "Rem_RawEVS")) {
      key <- paste(comparison_name, analysis_label, sep = "__")
      if (!is.null(analysis_store[[key]])) {
        p <- plot_pca_support(analysis_store[[key]], paste0(comparison_name, "\n", analysis_label_pretty(analysis_label)))
        out_name <- paste0("PCA_", comparison_name, "_", analysis_label, ".png")
        figure_status_rows[[length(figure_status_rows) + 1L]] <- safe_save_figure(p, file.path(figure_dir, out_name), 5.9, 4.8)
      }
    }
  }
}

pca_raw_plots <- list()
for (comparison_name in comparison_order) {
  key <- paste(comparison_name, "Raw", sep = "__")
  if (!is.null(analysis_store[[key]])) pca_raw_plots[[comparison_name]] <- plot_pca_support(analysis_store[[key]], comparison_name)
}
pca_raw_panel <- arrange_with_one_legend(pca_raw_plots, "PCA: Original", ncol = length(comparison_order))
figure_status_rows[[length(figure_status_rows) + 1L]] <- safe_save_figure(pca_raw_panel, file.path(figure_dir, "PCA_Original.png"), 18.0, 5.8)

for (dataset_prefix in c("Lead", "Rem")) {
  pca_plots <- list()
  for (mode in evs_modes) {
    for (comparison_name in comparison_order) {
      analysis_label <- paste0(dataset_prefix, "_", mode)
      key <- paste(comparison_name, analysis_label, sep = "__")
      if (!is.null(analysis_store[[key]])) pca_plots[[paste(comparison_name, mode, sep = "_")]] <- plot_pca_support(analysis_store[[key]], paste0(comparison_name, "\n", mode))
    }
  }
  panel_title <- if (dataset_prefix == "Lead") "PCA: Leading Edge" else "PCA: Remainder"
  panel_file <- if (dataset_prefix == "Lead") "PCA_Lead.png" else "PCA_Remainder.png"
  figure_status_rows[[length(figure_status_rows) + 1L]] <- safe_save_figure(arrange_with_one_legend(pca_plots, panel_title, ncol = length(comparison_order)), file.path(figure_dir, panel_file), 18.0, 9.4)
}

pca_var_rows <- list()
for (comparison_name in comparison_order) {
  for (analysis_label in c("Raw", "Lead_NormEVS", "Lead_RawEVS", "Rem_NormEVS", "Rem_RawEVS")) {
    key <- paste(comparison_name, analysis_label, sep = "__")
    if (!is.null(analysis_store[[key]])) pca_var_rows[[key]] <- extract_pca_variance_rows(analysis_store[[key]], comparison_name, analysis_label)
  }
}

pca_var_df <- if (length(pca_var_rows) > 0L) bind_rows(pca_var_rows) else data.frame()
if (nrow(pca_var_df) > 0L) {
  save_csv(rename_columns_existing(pca_var_df, c(comparison_name = "comparison", analysis_label = "analysis", analysis_pretty = "panel", dataset_group = "dataset", evs_mode = "evs_mode", component = "component", variance_explained_pct = "variance_explained_pct", score_variance = "score_variance")), file.path(output_dir, "Summary_PCA.csv"))
  pca_pc1_compare_df <- pca_pc1_evs_plot_df(pca_var_df)
  if (nrow(pca_pc1_compare_df) > 0L) save_csv(rename_columns_existing(pca_pc1_compare_df, c(comparison_name = "comparison", analysis_label = "analysis", analysis_pretty = "panel", dataset_group = "dataset", evs_mode = "evs_mode", pc1_score_variance = "pc1_score_variance", original_pc1_score_variance = "original_pc1_score_variance", pc1_score_variance_ratio_to_original = "pc1_score_variance_ratio_to_original", pc1_score_variance_percent_of_original = "pc1_score_variance_percent_of_original")), file.path(output_dir, "Summary_PCA_PC1_ScoreVariance_EVS.csv"))
  figure_status_rows[[length(figure_status_rows) + 1L]] <- safe_save_figure(plot_pc1_evs_variance_absolute(pca_var_df), file.path(figure_dir, "PC1_Variance_EVS.png"), 15.8, 8.8)
}

load_rows <- list()
cutoff_rows <- list()
for (mode in evs_modes) {
  for (comparison_name in comparison_order) {
    key <- paste(comparison_name, mode, sep = "__")
    comp_obj <- comparison_store[[comparison_name]]
    if (!is.null(evs_store[[key]]) && !is.null(comp_obj)) {
      tmp <- evs_loading_distribution_rows(comp_obj$count_matrix, evs_store[[key]], comp_obj$coldata, comparison_name, mode)
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
  pc1_distribution_path <- file.path(output_dir, "Summary_PC1VarianceDistributions.csv")
  if (isTRUE(export_full_pc1_variance_distributions)) {
    save_csv(rename_columns_existing(load_df, c(comparison_name = "comparison", evs_mode = "evs_mode", dataset_group = "dataset", scope = "scope", pc1_loading_abs = "abs_pc1_loading", pc1_score_variance = "pc1_score_variance", pc1_variance_contribution = "pc1_variance_contribution")), pc1_distribution_path)
  } else {
    if (isTRUE(remove_skipped_large_exports_from_disk) && file.exists(pc1_distribution_path)) unlink(pc1_distribution_path, force = TRUE)
    log_message("Skipped Summary_PC1VarianceDistributions.csv because export_full_pc1_variance_distributions is FALSE.")
  }
  if (nrow(cutoff_df) > 0L) save_csv(rename_columns_existing(cutoff_df, c(comparison_name = "comparison", evs_mode = "evs_mode", dataset_group = "dataset", scope = "scope", pc1_loading_cutoff = "pc1_loading_cutoff", pc1_score_variance = "pc1_score_variance", pc1_variance_contribution_cutoff = "pc1_variance_contribution_cutoff")), file.path(output_dir, "Summary_PC1VarianceCutoffs.csv"))
  figure_status_rows[[length(figure_status_rows) + 1L]] <- safe_save_figure(plot_loading_histograms(load_df, "PC1 variance-contribution frequency after EVS"), file.path(figure_dir, "PC1_Variance_Distributions.png"), 15.8, 12.0)
}

figure_status_df <- if (length(figure_status_rows) > 0L) bind_rows(figure_status_rows) else data.frame()
save_csv(figure_status_df, file.path(output_dir, "Figure_Export_Status.csv"))
if (nrow(figure_status_df) > 0L && any(figure_status_df$status == "failed")) {
  save_csv(figure_status_df[figure_status_df$status == "failed", , drop = FALSE], file.path(output_dir, "Figure_Export_Failures.csv"))
  fail_pipeline("One or more figure exports failed. See Figure_Export_Failures.csv.")
}

if (isTRUE(run_simulation_validation)) {
  simulation_summary <- run_sequence_simulation_validation()
  print(simulation_summary)
}

# =============================================================================
# COMPLETION
# =============================================================================

write_sequence_methods()
write_package_versions()
write_run_session_info()
cleanup_run_artifacts()
remove_skipped_large_export_files()
manifest <- write_manifest()
git_commit_and_push()

cat("\n=====================================================\n")
cat("SEQUENCE pipeline complete.\n")
cat("Output directory: ", output_dir, "\n", sep = "")
cat("Exported files: ", nrow(manifest), "\n", sep = "")
cat("=====================================================\n\n")
print(summary_df)
