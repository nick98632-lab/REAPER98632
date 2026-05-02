#!/usr/bin/env Rscript

# =============================================================================
# SEQUENCE FINAL PIPELINE
# DESeq2 + Eigenvector Splitting + empirical-null HC/HBFSS
# This file is a complete rewrite, not a patch.
# No legacy changepoint logic, no hidden cutoff fallback, no log transform for EVS.
# Standard effects use DESeq2 Wald BH < 10% with |shrunken LFC| >= 1.
# Strong effects are DESeq2 greaterAbs alternative-hypothesis calls: BH < 10% with |shrunken LFC| >= 1.
# Weak effects are DESeq2 lessAbs alternative-hypothesis calls: BH < 10%, |shrunken LFC| < 1, plus the HBFSS parabolic cutoff.
# =============================================================================

required_packages <- c(
  "DESeq2", "apeglm", "fdrtool", "ggplot2", "ggrepel",
  "dplyr", "gridExtra", "grid", "scales", "S4Vectors", "SummarizedExperiment"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop("Install missing package(s): ", paste(missing_packages, collapse = ", "))
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

options(stringsAsFactors = FALSE)

# =============================================================================
# SETTINGS
# =============================================================================

count_file_candidates <- c(
  file.path("data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"),
  "WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
)

alpha_level <- 0.10
strong_alpha_level <- 0.10
weak_alpha_level <- 0.10
lfc_boundary <- 1.0
top_n_target <- 5000L
hc_invalid_at_or_above <- 0.95
calculation_probability_floor <- .Machine$double.xmin
plot_probability_floor <- 1e-16
figure_dpi <- 320
base_theme_size <- 10
n_top_labels <- 8L
n_top_labels_per_class <- 2L

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
  Weak = "Weak",
  Strong = "Strong",
  Std = "Standard",
  HBFSS = "HBFSS-only"
)
class_shapes <- c(BG = 21, Weak = 24, Strong = 22, Std = 23, HBFSS = 25)
class_sizes <- c(BG = 0.55, Weak = 1.10, Strong = 1.15, Std = 1.10, HBFSS = 1.15)
class_alphas <- c(BG = 0.24, Weak = 0.92, Strong = 0.95, Std = 0.90, HBFSS = 0.95)

threshold_color <- "#A65628"
treatment_color <- "#1F78B4"
control_color <- "#4D4D4D"
pca_component_colors <- c(PC1 = "#1F78B4", PC2 = "#E31A1C")

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
  )
)
rownames(sample_metadata) <- sample_metadata$id
sample_metadata$condition <- factor(sample_metadata$condition, levels = c("untrt", "trt"))

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  treatment_prefix = c("R0", "R2", "R4", "R8"),
  control_prefix = c("ZT6", "ZT8", "ZT10", "ZT14")
)

analysis_tracks <- c("NormEVS", "RawEVS")

# =============================================================================
# BASIC HELPERS
# =============================================================================

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0L) return(NA_character_)
  normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = FALSE)
}

find_repo_root <- function() {
  candidates <- unique(c(dirname(script_path()), getwd(), dirname(getwd()), "/root/REAPER98632"))
  candidates <- candidates[!is.na(candidates) & dir.exists(candidates)]

  for (candidate in candidates) {
    if (dir.exists(file.path(candidate, ".git"))) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"))) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

repo_root <- find_repo_root()
output_dir <- file.path(repo_root, "exports", "manuscript_final_clean")
figure_dir <- file.path(output_dir, "manuscript_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

resolve_file <- function(candidates, label) {
  candidates <- c(file.path(repo_root, candidates), candidates)
  hits <- unique(candidates[file.exists(candidates)])
  if (length(hits) == 0L) {
    stop("Could not find ", label, ". Tried: ", paste(unique(candidates), collapse = " | "))
  }
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

stop_missing_columns <- function(df, required, label) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0L) {
    stop(label, " missing column(s): ", paste(missing, collapse = ", "))
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

  if (any(!is.finite(x) | is.na(x))) stop(label, " has NA or non-finite count values.")
  if (any(x < 0)) stop(label, " has negative count values.")

  rounded <- round(x)
  if (any(abs(x - rounded) > 1e-6)) {
    warning(label, " had non-integer values; rounded for DESeq2.")
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

neglog10_probability <- function(x, floor = plot_probability_floor) {
  -log10(clip_probability(x, floor = floor))
}

save_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE)
  invisible(path)
}

concise_analysis_summary <- function(df) {
  out <- df
  map <- c(
    comparison_name = "comp",
    analysis_label = "label",
    dataset_key = "set",
    n_features = "n",
    n_lessAbs_alt = "lessAbs",
    n_lessAbs_alt_HBFSS = "weakH",
    n_weak = "weak",
    n_strong_alt = "strongAlt",
    n_strong = "strong",
    n_standard = "std",
    n_hbfss_total = "H",
    n_hbfss_added_display = "Honly",
    hc_p_threshold = "HCp",
    hbfss_threshold = "Htau",
    alpha_level = "alpha",
    strong_alpha_level = "alphaS",
    weak_alpha_level = "alphaW",
    lfc_boundary = "lfc"
  )
  hit <- intersect(names(map), names(out))
  names(out)[match(hit, names(out))] <- unname(map[hit])
  out
}

concise_evs_summary <- function(df) {
  out <- df
  map <- c(
    comparison_name = "comp",
    track = "mode",
    top_n_target = "topN",
    treatment_top_n_used = "trtN",
    control_top_n_used = "ctlN",
    treatment_loading_cutoff = "trtCut",
    control_loading_cutoff = "ctlCut",
    leading_edge_n = "leadN",
    remainder_n = "remN",
    evs_input = "input"
  )
  hit <- intersect(names(map), names(out))
  names(out)[match(hit, names(out))] <- unname(map[hit])
  out
}

safe_csv_name <- function(...) {
  paste0(gsub("[^A-Za-z0-9_.-]+", "_", paste(..., sep = "_")), ".csv")
}

analysis_table_path <- function(comparison_name, analysis_label, table_label) {
  path <- file.path(output_dir, comparison_name, analysis_label, "tables")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  file.path(path, safe_csv_name("Table", comparison_name, analysis_label, table_label))
}

save_plot <- function(plot_obj, path, width, height) {
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
  invisible(path)
}

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_theme_size + 1, hjust = 0.5),
      plot.subtitle = element_text(size = base_theme_size - 1, hjust = 0.5),
      plot.caption = element_text(size = base_theme_size - 4, color = "grey30", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_theme_size - 1),
      legend.spacing.x = unit(4, "pt"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey88"),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(8, 10, 8, 8)
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
  title_grob <- grid::textGrob(title, gp = grid::gpar(fontface = "bold", cex = 1.08))

  if (is.null(legend)) return(gridExtra::arrangeGrob(body, ncol = 1, top = title_grob))
  gridExtra::arrangeGrob(body, legend, ncol = 1, heights = c(12, 1.1), top = title_grob)
}

# =============================================================================
# DATA IMPORT
# =============================================================================

count_file <- resolve_file(count_file_candidates, "WTTS count matrix")
message("Using count file: ", count_file)
message("Output directory: ", output_dir)

count_df <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
count_df <- as.data.frame(count_df, stringsAsFactors = FALSE)

stop_missing_columns(count_df, c("OrigID", "Symbol"), "Count matrix")
stop_missing_columns(count_df, sample_metadata$id, "Count matrix")

count_df$OrigID <- trimws(as.character(count_df$OrigID))
count_df$Symbol <- trimws(as.character(count_df$Symbol))

for (sample_id in sample_metadata$id) {
  count_df[[sample_id]] <- clean_count_column(count_df[[sample_id]])
}

valid_rows <- !is.na(count_df$OrigID) & nzchar(count_df$OrigID)
valid_rows <- valid_rows & rowSums(is.na(count_df[, sample_metadata$id, drop = FALSE])) == 0L
count_df <- count_df[valid_rows, , drop = FALSE]

count_df$feature_id <- make.unique(count_df$OrigID, sep = "_dup")
rownames(count_df) <- count_df$feature_id

annotation_df <- data.frame(
  feature_id = count_df$feature_id,
  orig_id = count_df$OrigID,
  gene_symbol = ifelse(nzchar(count_df$Symbol), count_df$Symbol, NA_character_)
) %>%
  distinct(feature_id, .keep_all = TRUE)

# =============================================================================
# COMPARISON DATA
# =============================================================================

prepare_comparison <- function(i) {
  row <- comparison_table[i, , drop = FALSE]
  treatment_ids <- sample_metadata$id[grepl(paste0("^", row$treatment_prefix, "_"), sample_metadata$id)]
  control_ids <- sample_metadata$id[grepl(paste0("^", row$control_prefix, "_"), sample_metadata$id)]
  sample_ids <- c(treatment_ids, control_ids)

  if (length(treatment_ids) < 2L || length(control_ids) < 2L) {
    stop("Comparison ", row$comparison_name, " does not have enough samples.")
  }

  coldata <- sample_metadata[sample_ids, "condition", drop = FALSE]
  counts <- count_df[, sample_ids, drop = FALSE]
  rownames(counts) <- count_df$feature_id
  counts <- as_integer_count_matrix(counts, paste0(row$comparison_name, " count matrix"))

  if (!identical(colnames(counts), rownames(coldata))) {
    stop("Sample order mismatch for ", row$comparison_name)
  }

  list(
    comparison_name = row$comparison_name,
    count_matrix = counts,
    coldata = coldata
  )
}

# =============================================================================
# EIGENVECTOR SPLITTING
# =============================================================================

make_norm_evs_matrix <- function(count_matrix, coldata) {
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds)
  as.matrix(counts(dds, normalized = TRUE))
}

get_evs_matrix <- function(count_matrix, coldata, track) {
  if (track == "RawEVS") {
    x <- matrix(as.numeric(count_matrix), nrow = nrow(count_matrix), dimnames = dimnames(count_matrix))
    return(x)
  }

  if (track == "NormEVS") {
    return(make_norm_evs_matrix(count_matrix, coldata))
  }

  stop("Unknown EVS track: ", track)
}

condition_pc1_rank <- function(evs_matrix, sample_ids, condition_name) {
  x <- evs_matrix[, sample_ids, drop = FALSE]
  storage.mode(x) <- "numeric"

  variable_feature <- apply(x, 1, stats::var, na.rm = TRUE) > 0
  if (!any(variable_feature)) stop("No variable features for EVS PCA in ", condition_name)

  x <- x[variable_feature, , drop = FALSE]
  pca <- stats::prcomp(t(x), center = TRUE, scale. = FALSE)
  loading <- pca$rotation[, 1]

  rank_df <- data.frame(
    feature_id = names(loading),
    pc1_loading = as.numeric(loading),
    pc1_loading_abs = abs(as.numeric(loading))
  )
  rank_df <- rank_df[order(rank_df$pc1_loading_abs, decreasing = TRUE), , drop = FALSE]
  rank_df$evs_rank <- seq_len(nrow(rank_df))

  top_n_used <- min(top_n_target, nrow(rank_df))
  rank_df$evs_selected <- rank_df$evs_rank <= top_n_used
  rank_df$top_n_used <- top_n_used
  rank_df$loading_cutoff_at_top_n <- rank_df$pc1_loading_abs[top_n_used]

  list(pca = pca, rank_df = rank_df, top_n_used = top_n_used)
}

build_evs_split <- function(count_matrix, coldata, comparison_name, track) {
  evs_matrix <- get_evs_matrix(count_matrix, coldata, track)
  all_features <- rownames(count_matrix)
  treatment_samples <- rownames(coldata)[coldata$condition == "trt"]
  control_samples <- rownames(coldata)[coldata$condition == "untrt"]

  treatment_rank <- condition_pc1_rank(evs_matrix, treatment_samples, "treatment")
  control_rank <- condition_pc1_rank(evs_matrix, control_samples, "control")

  treatment_top <- treatment_rank$rank_df$feature_id[treatment_rank$rank_df$evs_selected]
  control_top <- control_rank$rank_df$feature_id[control_rank$rank_df$evs_selected]
  leading_ids <- intersect(union(treatment_top, control_top), all_features)
  remainder_ids <- setdiff(all_features, leading_ids)

  if (length(leading_ids) == 0L) stop("EVS leading edge empty for ", comparison_name, " ", track)
  if (length(remainder_ids) == 0L) {
    stop("EVS remainder empty for ", comparison_name, " ", track, ". Reduce top_n_target.")
  }

  trt_rank <- treatment_rank$rank_df
  ctl_rank <- control_rank$rank_df
  names(trt_rank)[names(trt_rank) != "feature_id"] <- paste0("trt_", names(trt_rank)[names(trt_rank) != "feature_id"])
  names(ctl_rank)[names(ctl_rank) != "feature_id"] <- paste0("ctl_", names(ctl_rank)[names(ctl_rank) != "feature_id"])

  joint_rank <- data.frame(feature_id = all_features) %>%
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
    top_n_target = top_n_target,
    treatment_top_n_used = treatment_rank$top_n_used,
    control_top_n_used = control_rank$top_n_used,
    treatment_loading_cutoff = treatment_rank$rank_df$loading_cutoff_at_top_n[1],
    control_loading_cutoff = control_rank$rank_df$loading_cutoff_at_top_n[1],
    leading_edge_n = length(leading_ids),
    remainder_n = length(remainder_ids),
    evs_input = ifelse(track == "NormEVS", "DESeq2 normalized counts; no log transform", "Raw counts; no log transform")
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
  names <- resultsNames(dds)
  if ("condition_trt_vs_untrt" %in% names) return("condition_trt_vs_untrt")
  hit <- grep("condition.*trt.*vs.*untrt", names, value = TRUE)
  if (length(hit) > 0L) return(hit[1])
  stop("Could not find trt-vs-untrt coefficient. Available: ", paste(names, collapse = ", "))
}

empirical_null <- function(statistic, label) {
  valid <- is.finite(statistic) & !is.na(statistic)
  if (sum(valid) < 5L) stop(label, " has fewer than five finite Wald statistics.")

  fit <- fdrtool::fdrtool(
    statistic[valid],
    statistic = "normal",
    plot = FALSE,
    verbose = FALSE,
    cutoff.method = "fndr",
    pct0 = 0.75
  )

  p <- rep(NA_real_, length(statistic))
  q <- rep(NA_real_, length(statistic))
  lfdr <- rep(NA_real_, length(statistic))
  p[valid] <- clip_probability(fit$pval)
  q[valid] <- clip_probability(fit$qval)
  lfdr[valid] <- suppressWarnings(as.numeric(fit$lfdr))

  bh <- rep(NA_real_, length(statistic))
  ok <- !is.na(p) & is.finite(p)
  bh[ok] <- p.adjust(p[ok], method = "BH")

  list(empirical_p = p, empirical_q = q, empirical_lfdr = lfdr, empirical_bh = bh)
}

hc_threshold <- function(empirical_p) {
  p <- sort(clip_probability(empirical_p), decreasing = FALSE, na.last = NA)
  if (length(p) < 5L) return(NA_real_)

  threshold <- suppressWarnings(tryCatch(fdrtool::hc.thresh(p), error = function(e) NA_real_))
  threshold <- as.numeric(threshold[1])

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

  df$hc_pass <- !is.na(hc_p) &
    !is.na(df$empirical_p) &
    is.finite(df$empirical_p) &
    df$empirical_p <= hc_p

  df$HBFSS <- abs(df$lfc_shrunk) * df$neglog10_empirical_p_calc

  df$hbfss_flag <- !is.na(hbfss_cutoff) &
    !is.na(df$HBFSS) &
    is.finite(df$HBFSS) &
    df$HBFSS >= hbfss_cutoff &
    df$hc_pass

  df$lessAbs_alt_flag <- !is.na(df$lessAbs_padj) &
    df$lessAbs_padj < weak_alpha_level &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) < lfc_boundary

  # Weak effects are DESeq2 lessAbs alternative-hypothesis calls under the LFC boundary.
  # For volcano display they must also pass the same HBFSS parabolic cutoff used in the manuscript boundary.
  df$weak_flag <- df$lessAbs_alt_flag & df$hbfss_flag
  df$weak_hbfss_flag <- df$weak_flag

  df$display_strong <- df$strong_flag
  df$display_standard <- df$standard_flag & !df$display_strong
  df$display_weak <- df$weak_flag & !df$display_strong & !df$display_standard
  # Volcano/display-only HBFSS points are restricted to HBFSS-significant features that are not already BG/Weak/Strong/Std.
  df$display_hbfss <- df$hbfss_flag & !df$display_strong & !df$display_standard & !df$display_weak

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

  strong_df <- data.frame(
    feature_id = rownames(strong),
    greaterAbs_pvalue = strong$pvalue,
    greaterAbs_padj = strong$padj
  )

  weak_df <- data.frame(
    feature_id = rownames(weak),
    lessAbs_pvalue = weak$pvalue,
    lessAbs_padj = weak$padj
  )

  shrink_df <- data.frame(
    feature_id = rownames(shrunk),
    lfc_shrunk = as.data.frame(shrunk)$log2FoldChange
  )

  empirical <- empirical_null(df$stat, label)
  df$empirical_p <- empirical$empirical_p
  df$empirical_q <- empirical$empirical_q
  df$empirical_lfdr <- empirical$empirical_lfdr
  df$empirical_bh <- empirical$empirical_bh

  # Separate numerical floors:
  # calculation floor protects HBFSS from zero p-values without changing thresholds.
  # plot floor prevents extreme p-values from stretching figures.
  # The dataset-specific HC threshold is used for hc_pass, the HC line, and the HBFSS cutoff;
  # it is not used as a p-value floor because that would flatten all more-significant points.
  df$empirical_p_calc <- clip_probability(df$empirical_p, floor = calculation_probability_floor)
  df$empirical_p_plot <- clip_probability(df$empirical_p, floor = plot_probability_floor)
  df$neglog10_empirical_p_calc <- -log10(df$empirical_p_calc)
  df$neglog10_empirical_p_plot <- -log10(df$empirical_p_plot)
  df$neglog10_empirical_p <- df$neglog10_empirical_p_plot

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
  keep_disp_cols <- intersect(
    c("feature_id", "dispGeneEst", "dispFit", "dispersion", "dispIter", "dispOutlier"),
    colnames(dispersion_df)
  )
  dispersion_df <- dispersion_df[, keep_disp_cols, drop = FALSE]

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
    "empirical_p", "empirical_p_calc", "empirical_p_plot", "empirical_bh", "empirical_q", "empirical_lfdr",
    "neglog10_empirical_p_calc", "neglog10_empirical_p_plot", "neglog10_empirical_p", "HBFSS",
    "hc_p_threshold_dataset", "hbfss_threshold_dataset",
    "hc_pass", "standard_flag", "strong_alt_flag", "strong_flag",
    "lessAbs_alt_flag", "weak_hbfss_flag", "weak_flag", "hbfss_flag",
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
    n_weak = sum(df$weak_flag, na.rm = TRUE),
    n_strong_alt = sum(df$strong_alt_flag, na.rm = TRUE),
    n_strong = sum(df$strong_flag, na.rm = TRUE),
    n_standard = sum(df$standard_flag, na.rm = TRUE),
    n_hbfss_total = sum(df$hbfss_flag, na.rm = TRUE),
    n_hbfss_added_display = sum(df$display_hbfss, na.rm = TRUE),
    hc_p_threshold = hc_p,
    hbfss_threshold = hbfss_cutoff,
    alpha_level = alpha_level,
    strong_alpha_level = strong_alpha_level,
    weak_alpha_level = weak_alpha_level,
    lfc_boundary = lfc_boundary
  )

  list(dds = dds, results = df, summary = summary)
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
  plot_df <- plot_df[order(plot_df$draw_order, plot_df$neglog10_empirical_p), , drop = FALSE]
  plot_df
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
  labels <- labels[order(match(as.character(labels$final_class), class_priority), -labels$HBFSS, labels$empirical_p, -abs(labels$lfc_shrunk), na.last = TRUE), , drop = FALSE]
  labels[seq_len(min(n_top_labels, nrow(labels))), , drop = FALSE]
}

hbfss_boundary <- function(plot_df, hc_p, hbfss_cutoff, y_limit) {
  if (!is.finite(hbfss_cutoff) || is.na(hbfss_cutoff) || hbfss_cutoff <= 0) return(NULL)

  hc_y <- if (!is.na(hc_p) && is.finite(hc_p) && hc_p > 0 && hc_p < 1) -log10(hc_p) else NA_real_
  x_limit <- max(abs(plot_df$lfc_shrunk), lfc_boundary * 1.1, na.rm = TRUE)
  x_start <- max(0.05, hbfss_cutoff / max(y_limit, 1e-6))
  x_abs <- seq(x_start, x_limit, length.out = 600)
  y <- hbfss_cutoff / x_abs

  if (!is.na(hc_y) && is.finite(hc_y)) y <- pmax(y, hc_y)
  keep <- is.finite(y) & y >= 0 & y <= y_limit
  if (!any(keep)) return(NULL)

  x_abs <- x_abs[keep]
  y <- y[keep]

  rbind(
    data.frame(x = -rev(x_abs), y = rev(y)),
    data.frame(x = x_abs, y = y)
  )
}

plot_volcano <- function(df, title) {
  plot_df <- volcano_data(df)
  if (nrow(plot_df) == 0L) stop("No finite rows for volcano plot: ", title)

  label_df <- volcano_labels(plot_df)
  hc_p <- suppressWarnings(as.numeric(plot_df$hc_p_threshold_dataset[1]))
  hbfss_cutoff <- suppressWarnings(as.numeric(plot_df$hbfss_threshold_dataset[1]))
  hc_y <- if (!is.na(hc_p) && is.finite(hc_p) && hc_p > 0 && hc_p < 1) -log10(hc_p) else NA_real_
  y_limit <- max(plot_df$neglog10_empirical_p, hc_y, na.rm = TRUE) * 1.06
  boundary <- hbfss_boundary(plot_df, hc_p, hbfss_cutoff, y_limit)
  x_min <- min(plot_df$lfc_shrunk, na.rm = TRUE)
  x_max <- max(plot_df$lfc_shrunk, na.rm = TRUE)
  x_span <- max(x_max - x_min, 1e-6)
  threshold_label_y <- if (!is.na(hc_y) && is.finite(hc_y)) max(0.06 * y_limit, hc_y - 0.06 * y_limit) else 0.08 * y_limit
  hc_label_x <- x_min + 0.03 * x_span
  hbfss_label_x <- x_min + 0.34 * x_span

  caption <- paste0(
    "W=", sum(df$weak_flag, na.rm = TRUE),
    "  S=", sum(df$strong_flag, na.rm = TRUE),
    "  Std=", sum(df$standard_flag, na.rm = TRUE),
    "  H=", sum(df$hbfss_flag, na.rm = TRUE),
    "  H+=", sum(df$display_hbfss, na.rm = TRUE)
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
      color = "none",
      shape = "none",
      fill = guide_legend(
        override.aes = list(
          shape = unname(class_shapes[class_levels]),
          color = unname(class_colors[class_levels]),
          fill = unname(class_colors[class_levels]),
          size = rep(3.1, length(class_levels)),
          alpha = rep(1, length(class_levels)),
          stroke = rep(0.55, length(class_levels))
        ),
        nrow = 1
      )
    ) +
    labs(
      title = title,
      x = "Shrunken log2 fold change",
      y = expression(-log[10]("empirical p")),
      caption = caption
    ) +
    coord_cartesian(ylim = c(0, y_limit), clip = "off") +
    manuscript_theme() +
    theme(plot.caption = element_text(size = base_theme_size - 4))

  if (!is.na(hc_y) && is.finite(hc_y)) {
    p <- p + geom_hline(yintercept = hc_y, linetype = "dotted", linewidth = 0.60, color = threshold_color)
  }

  if (!is.null(boundary)) {
    p <- p + geom_line(data = boundary, aes(x = x, y = y), inherit.aes = FALSE, color = class_colors[["HBFSS"]], linewidth = 0.75)
  }

  if (!is.na(hc_p) && is.finite(hc_p) && !is.na(hc_y) && is.finite(hc_y)) {
    p <- p + annotate(
      "text",
      x = hc_label_x,
      y = threshold_label_y,
      label = paste0("HC p=", signif(hc_p, 3)),
      hjust = 0,
      vjust = 1,
      size = 2.15,
      color = threshold_color
    )
  }

  if (!is.na(hbfss_cutoff) && is.finite(hbfss_cutoff)) {
    p <- p + annotate(
      "text",
      x = hbfss_label_x,
      y = threshold_label_y,
      label = paste0("HBFSS \u03c4=", signif(hbfss_cutoff, 3)),
      hjust = 0,
      vjust = 1,
      size = 2.15,
      color = class_colors[["HBFSS"]]
    )
  }

  if (nrow(label_df) > 0L) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      aes(x = lfc_shrunk, y = neglog10_empirical_p, label = gene_label, color = final_class),
      inherit.aes = FALSE,
      show.legend = FALSE,
      size = 1.60,
      seed = 1,
      max.overlaps = Inf,
      force = 1.0,
      force_pull = 0.30,
      box.padding = 0.28,
      point.padding = 0.14,
      min.segment.length = 0,
      segment.alpha = 0.60,
      segment.size = 0.22
    )
  }

  p
}

plot_evs_rank <- function(evs_obj) {
  trt <- evs_obj$treatment_rank$rank_df %>%
    transmute(condition = "Treatment", evs_rank = evs_rank, pc1_loading_abs = pc1_loading_abs)
  ctl <- evs_obj$control_rank$rank_df %>%
    transmute(condition = "Control", evs_rank = evs_rank, pc1_loading_abs = pc1_loading_abs)
  plot_df <- bind_rows(trt, ctl)

  ggplot(plot_df, aes(evs_rank, pc1_loading_abs, color = condition)) +
    geom_line(linewidth = 0.45, alpha = 0.90) +
    geom_vline(xintercept = top_n_target, linetype = "dashed", linewidth = 0.55, color = threshold_color) +
    scale_color_manual(values = c(Treatment = treatment_color, Control = control_color), name = NULL) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title = paste0(evs_obj$comparison_name, "\n", evs_obj$track),
      x = "Rank",
      y = "Abs PC1",
      caption = paste0("topN=", top_n_target, "  Lead=", length(evs_obj$leading_ids), "  Rem=", length(evs_obj$remainder_ids))
    ) +
    manuscript_theme()
}

plot_evs_loading_histogram <- function(evs_obj) {
  trt <- evs_obj$treatment_rank$rank_df %>%
    transmute(condition = "Treatment", pc1_loading_abs = pc1_loading_abs, cutoff = loading_cutoff_at_top_n)
  ctl <- evs_obj$control_rank$rank_df %>%
    transmute(condition = "Control", pc1_loading_abs = pc1_loading_abs, cutoff = loading_cutoff_at_top_n)
  plot_df <- bind_rows(trt, ctl)
  cutoff_df <- plot_df %>% distinct(condition, cutoff)

  ggplot(plot_df, aes(pc1_loading_abs, fill = condition)) +
    geom_histogram(bins = 80, color = "grey30", linewidth = 0.10, alpha = 0.72, position = "identity") +
    geom_vline(data = cutoff_df, aes(xintercept = cutoff, color = condition), linetype = "dashed", linewidth = 0.60, show.legend = FALSE) +
    scale_fill_manual(values = c(Treatment = treatment_color, Control = control_color), name = NULL) +
    scale_color_manual(values = c(Treatment = treatment_color, Control = control_color)) +
    labs(
      title = paste0(evs_obj$comparison_name, "\n", evs_obj$track),
      x = "Abs PC1",
      y = "n",
      caption = NULL
    ) +
    manuscript_theme()
}

plot_discovery_counts <- function(summary_df) {
  long_df <- bind_rows(lapply(seq_len(nrow(summary_df)), function(i) {
    row <- summary_df[i, , drop = FALSE]
    data.frame(
      comparison_name = row$comparison_name,
      analysis_label = row$analysis_label,
      dataset_key = row$dataset_key,
      Class = factor(c("Weak", "Strong", "Std", "HBFSS"), levels = c("Weak", "Strong", "Std", "HBFSS")),
      # Histogram/discovery counts use all HBFSS-significant genes, not only the HBFSS-only display subset.
      Count = as.numeric(c(row$n_weak, row$n_strong, row$n_standard, row$n_hbfss_total))
    )
  }))

  save_csv(dplyr::rename(long_df, comp = comparison_name, label = analysis_label, set = dataset_key, cls = Class, n = Count), file.path(output_dir, "Table_Final_Discovery_Counts_Long.csv"))

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
    labs(
      title = "Counts",
      x = NULL,
      y = "n",
      caption = "BH<10%; H=all HBFSS; H+=H-only."
    ) +
    manuscript_theme() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}


pca_support_components <- function(fit_obj) {
  dds <- fit_obj$dds
  if (is.null(dds)) return(NULL)

  x <- as.matrix(DESeq2::counts(dds, normalized = TRUE))
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
    evs_mode = factor(evs_mode, levels = c("Original", "NormEVS", "RawEVS")),
    component = factor(c("PC1", "PC2"), levels = c("PC1", "PC2")),
    variance_explained_pct = as.numeric(comp$var_pct[c(1, 2)]),
    score_variance = as.numeric(comp$score_variance[c("PC1", "PC2")]),
    stringsAsFactors = FALSE
  )
}

plot_pca_variance_summary <- function(pca_var_df, value_col, title, y_label, caption) {
  if (nrow(pca_var_df) == 0L) return(NULL)

  ggplot(pca_var_df, aes(comparison_name, .data[[value_col]], fill = component)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.74, color = "grey25", linewidth = 0.15) +
    facet_wrap(~ analysis_pretty, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = pca_component_colors, drop = FALSE, name = NULL) +
    labs(
      title = title,
      x = NULL,
      y = y_label,
      caption = caption
    ) +
    manuscript_theme() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

loading_profile_from_matrix <- function(count_matrix, coldata, comparison_name, analysis_label, track) {
  evs_matrix <- get_evs_matrix(count_matrix, coldata, track)
  treatment_samples <- rownames(coldata)[coldata$condition == "trt"]
  control_samples <- rownames(coldata)[coldata$condition == "untrt"]

  treatment_rank <- condition_pc1_rank(evs_matrix, treatment_samples, "treatment")$rank_df
  control_rank <- condition_pc1_rank(evs_matrix, control_samples, "control")$rank_df

  bind_rows(
    data.frame(
      comparison_name = comparison_name,
      analysis_label = analysis_label,
      analysis_pretty = analysis_label_pretty(analysis_label),
      track = track,
      condition = factor("Treatment", levels = c("Control", "Treatment")),
      pc1_loading_abs = treatment_rank$pc1_loading_abs,
      stringsAsFactors = FALSE
    ),
    data.frame(
      comparison_name = comparison_name,
      analysis_label = analysis_label,
      analysis_pretty = analysis_label_pretty(analysis_label),
      track = track,
      condition = factor("Control", levels = c("Control", "Treatment")),
      pc1_loading_abs = control_rank$pc1_loading_abs,
      stringsAsFactors = FALSE
    )
  )
}

plot_loading_distribution_panel <- function(profile_df, title) {
  if (nrow(profile_df) == 0L) return(NULL)

  ggplot(profile_df, aes(pc1_loading_abs, fill = condition)) +
    geom_histogram(bins = 80, color = "grey30", linewidth = 0.10, alpha = 0.65, position = "identity") +
    scale_fill_manual(values = c(Control = control_color, Treatment = treatment_color), drop = FALSE, name = NULL) +
    labs(
      title = title,
      x = "Abs PC1",
      y = "n",
      caption = NULL
    ) +
    manuscript_theme()
}

pca_score_distribution_rows <- function(fit_obj, comparison_name, analysis_label) {
  comp <- pca_support_components(fit_obj)
  if (is.null(comp)) return(NULL)

  bind_rows(
    data.frame(
      comparison_name = comparison_name,
      analysis_label = analysis_label,
      component = factor("PC1", levels = c("PC1", "PC2")),
      condition = comp$pca_df$condition,
      score = comp$pca_df$PC1,
      stringsAsFactors = FALSE
    ),
    data.frame(
      comparison_name = comparison_name,
      analysis_label = analysis_label,
      component = factor("PC2", levels = c("PC1", "PC2")),
      condition = comp$pca_df$condition,
      score = comp$pca_df$PC2,
      stringsAsFactors = FALSE
    )
  )
}

plot_pca_score_distribution_panel <- function(score_df, title) {
  if (nrow(score_df) == 0L) return(NULL)

  ggplot(score_df, aes(score, fill = condition)) +
    geom_histogram(bins = 16, color = "grey30", linewidth = 0.10, alpha = 0.66, position = "identity") +
    facet_wrap(~ component, scales = "free", ncol = 1) +
    scale_fill_manual(values = c(Control = control_color, Treatment = treatment_color), drop = FALSE, name = NULL) +
    labs(
      title = title,
      x = "Score",
      y = "n",
      caption = NULL
    ) +
    manuscript_theme()
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

plot_empirical_hbfss_support <- function(df, title) {
  req <- c("empirical_p", "HBFSS", "final_class", "hc_p_threshold_dataset", "hbfss_threshold_dataset")
  if (length(setdiff(req, colnames(df))) > 0L) return(NULL)

  plot_df <- df[
    is.finite(df$empirical_p) & !is.na(df$empirical_p) &
      is.finite(df$HBFSS) & !is.na(df$HBFSS),
    ,
    drop = FALSE
  ]
  if (nrow(plot_df) == 0L) return(NULL)

  if (!"empirical_p_plot" %in% colnames(plot_df)) {
    plot_df$empirical_p_plot <- clip_probability(plot_df$empirical_p, floor = plot_probability_floor)
  }
  plot_df$final_class <- factor(as.character(plot_df$final_class), levels = class_levels)

  hc_p <- suppressWarnings(as.numeric(plot_df$hc_p_threshold_dataset[1]))
  hbfss_cutoff <- suppressWarnings(as.numeric(plot_df$hbfss_threshold_dataset[1]))

  p <- ggplot(
    plot_df,
    aes(empirical_p_plot, HBFSS, color = final_class, fill = final_class, shape = final_class)
  ) +
    geom_point(alpha = 0.62, size = 1.00, stroke = 0.25) +
    scale_color_manual(values = class_colors, breaks = class_levels, labels = class_labels[class_levels], drop = FALSE, name = NULL) +
    scale_fill_manual(values = class_colors, breaks = class_levels, labels = class_labels[class_levels], drop = FALSE, name = NULL) +
    scale_shape_manual(values = class_shapes, breaks = class_levels, labels = class_labels[class_levels], drop = FALSE, name = NULL) +
    guides(
      color = "none",
      shape = "none",
      fill = guide_legend(
        override.aes = list(
          shape = unname(class_shapes[class_levels]),
          color = unname(class_colors[class_levels]),
          fill = unname(class_colors[class_levels]),
          size = rep(3.1, length(class_levels)),
          alpha = rep(1, length(class_levels)),
          stroke = rep(0.55, length(class_levels))
        ),
        nrow = 1
      )
    ) +
    scale_x_log10(labels = scales::label_scientific()) +
    labs(
      title = title,
      x = "Empirical p",
      y = "HBFSS",
      caption = paste0(
        if (!is.na(hc_p) && is.finite(hc_p)) paste0("HCp=", signif(hc_p, 3)) else "HCp=NA",
        "  ",
        if (!is.na(hbfss_cutoff) && is.finite(hbfss_cutoff)) paste0("Hτ=", signif(hbfss_cutoff, 3)) else "Hτ=NA"
      )
    ) +
    manuscript_theme()

  if (!is.na(hc_p) && is.finite(hc_p) && hc_p > 0 && hc_p < 1) {
    p <- p + geom_vline(xintercept = hc_p, color = threshold_color, linewidth = 0.55, linetype = "dotted")
  }

  if (!is.na(hbfss_cutoff) && is.finite(hbfss_cutoff)) {
    p <- p + geom_hline(yintercept = hbfss_cutoff, color = threshold_color, linewidth = 0.55, linetype = "dashed")
  }

  p
}


topn_loading_mass_rows <- function(evs_obj) {
  one_side <- function(rank_df, condition_label) {
    total <- sum(rank_df$pc1_loading_abs, na.rm = TRUE)
    selected <- sum(rank_df$pc1_loading_abs[rank_df$evs_selected], na.rm = TRUE)
    data.frame(
      comparison_name = evs_obj$comparison_name,
      mode = evs_obj$track,
      condition = factor(condition_label, levels = c("Control", "Treatment")),
      topN_fraction = ifelse(total > 0, selected / total, NA_real_),
      stringsAsFactors = FALSE
    )
  }

  bind_rows(
    one_side(evs_obj$control_rank$rank_df, "Control"),
    one_side(evs_obj$treatment_rank$rank_df, "Treatment")
  )
}

plot_topn_loading_mass <- function(df) {
  if (nrow(df) == 0L) return(NULL)

  ggplot(df, aes(comparison_name, topN_fraction, fill = condition)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.70, color = "grey25", linewidth = 0.12) +
    facet_wrap(~ mode, ncol = 2) +
    scale_fill_manual(values = c(Control = control_color, Treatment = treatment_color), drop = FALSE, name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
    labs(
      title = "TopN loading mass",
      x = NULL,
      y = "TopN / total",
      caption = NULL
    ) +
    manuscript_theme() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

plot_evs_benefit_summary <- function(pca_var_df) {
  if (nrow(pca_var_df) == 0L) return(NULL)

  pc1 <- pca_var_df[pca_var_df$component == "PC1" & pca_var_df$dataset_group %in% c("Lead", "Remainder"), , drop = FALSE]
  if (nrow(pc1) == 0L) return(NULL)

  wide <- reshape(
    pc1[, c("comparison_name", "evs_mode", "dataset_group", "score_variance")],
    idvar = c("comparison_name", "evs_mode"),
    timevar = "dataset_group",
    direction = "wide"
  )

  lead_col <- "score_variance.Lead"
  rem_col <- "score_variance.Remainder"
  if (!all(c(lead_col, rem_col) %in% colnames(wide))) return(NULL)

  wide$rem_lead_ratio <- wide[[rem_col]] / wide[[lead_col]]
  wide <- wide[is.finite(wide$rem_lead_ratio), , drop = FALSE]
  if (nrow(wide) == 0L) return(NULL)

  wide$evs_mode <- factor(as.character(wide$evs_mode), levels = c("NormEVS", "RawEVS"))

  ggplot(wide, aes(comparison_name, rem_lead_ratio, fill = evs_mode)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey45") +
    geom_col(position = position_dodge(width = 0.78), width = 0.70, color = "grey25", linewidth = 0.12) +
    scale_fill_manual(values = c(NormEVS = treatment_color, RawEVS = control_color), drop = FALSE, name = NULL) +
    labs(
      title = "EVS benefit",
      x = NULL,
      y = "Rem / Lead PC1 var",
      caption = "Lower = cleaner residual."
    ) +
    manuscript_theme() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}


# =============================================================================
# PIPELINE EXECUTION
# =============================================================================

old_figures <- list.files(figure_dir, pattern = "^Figure_Manuscript_.*\\.png$", full.names = TRUE)
if (length(old_figures) > 0L) unlink(old_figures)

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
    run_store_save(
      comparison$count_matrix,
      comparison$coldata,
      comparison_name,
      "Raw",
      "Raw"
    )

    for (track in analysis_tracks) {
      evs <- build_evs_split(comparison$count_matrix, comparison$coldata, comparison_name, track)
      evs_key <- paste(comparison_name, track, sep = "__")
      evs_store[[evs_key]] <- evs

      save_csv(evs$joint_rank, analysis_table_path(comparison_name, track, "EVS_Rank_Table"))
      save_csv(concise_evs_summary(evs$summary), analysis_table_path(comparison_name, track, "EVS_Summary"))

      run_store_save(
        evs$leading_matrix,
        comparison$coldata,
        comparison_name,
        paste0("Lead_", track),
        "Lead"
      )

      run_store_save(
        evs$remainder_matrix,
        comparison$coldata,
        comparison_name,
        paste0("Rem_", track),
        "Rem"
      )
    }
  }, error = function(e) {
    failure_rows[[comparison_name]] <<- data.frame(
      comparison_name = comparison_name,
      error_message = conditionMessage(e)
    )
    message("FAILED: ", comparison_name, ": ", conditionMessage(e))
  })
}

summary_df <- if (length(summary_rows) > 0L) bind_rows(summary_rows) else data.frame()
if (nrow(summary_df) > 0L) save_csv(concise_analysis_summary(summary_df), file.path(output_dir, "Table_Overall_Summary.csv"))

if (length(failure_rows) > 0L) {
  failure_df <- bind_rows(failure_rows)
  save_csv(failure_df, file.path(output_dir, "Table_Failed.csv"))
}

# =============================================================================
# MANUSCRIPT FIGURE EXPORT
# =============================================================================

comparison_order <- comparison_table$comparison_name

if (nrow(summary_df) > 0L) {
  raw_plots <- list()
  for (comparison_name in comparison_order) {
    key <- paste(comparison_name, "Raw", sep = "__")
    if (!is.null(analysis_store[[key]])) raw_plots[[comparison_name]] <- plot_volcano(analysis_store[[key]]$results, comparison_name)
  }

  raw_panel <- arrange_with_one_legend(raw_plots, "Volcano: Raw", ncol = length(comparison_order))
  if (!is.null(raw_panel)) {
    save_plot(raw_panel, file.path(figure_dir, "Figure_Manuscript_Volcano_Raw_AllComparisons.png"), width = 18.0, height = 5.8)
  }

  for (dataset_prefix in c("Lead", "Rem")) {
    plots <- list()
    for (track in analysis_tracks) {
      for (comparison_name in comparison_order) {
        analysis_label <- paste0(dataset_prefix, "_", track)
        key <- paste(comparison_name, analysis_label, sep = "__")
        if (!is.null(analysis_store[[key]])) {
          plots[[paste(comparison_name, track, sep = "_")]] <- plot_volcano(analysis_store[[key]]$results, paste0(comparison_name, "\n", track))
        }
      }
    }

    title <- if (dataset_prefix == "Lead") "Volcano: Lead" else "Volcano: Rem"

    panel <- arrange_with_one_legend(plots, title, ncol = length(comparison_order))
    if (!is.null(panel)) {
      save_plot(
        panel,
        file.path(figure_dir, paste0("Figure_Manuscript_Volcano_", dataset_prefix, "_AllComparisons_NormEVS_vs_RawEVS.png")),
        width = 18.0,
        height = 9.8
      )
    }
  }

  count_plot <- plot_discovery_counts(summary_df)
  save_plot(count_plot, file.path(figure_dir, "Figure_Manuscript_Discovery_Counts.png"), width = 14.5, height = 8.2)

  rank_plots <- list()
  hist_plots <- list()
  for (track in analysis_tracks) {
    for (comparison_name in comparison_order) {
      key <- paste(comparison_name, track, sep = "__")
      if (!is.null(evs_store[[key]])) {
        rank_plots[[key]] <- plot_evs_rank(evs_store[[key]])
        hist_plots[[key]] <- plot_evs_loading_histogram(evs_store[[key]])
      }
    }
  }

  rank_panel <- arrange_with_one_legend(rank_plots, "EVS loading rank", ncol = length(comparison_order))
  if (!is.null(rank_panel)) {
    save_plot(rank_panel, file.path(figure_dir, "Figure_Manuscript_EVS_PC1_Loading_Rank.png"), width = 18.0, height = 9.4)
  }

  hist_panel <- arrange_with_one_legend(hist_plots, "EVS loading hist", ncol = length(comparison_order))
  if (!is.null(hist_panel)) {
    save_plot(hist_panel, file.path(figure_dir, "Figure_Manuscript_EVS_PC1_Loading_Histogram.png"), width = 18.0, height = 9.4)
  }

  pca_raw_plots <- list()
  empirical_raw_plots <- list()
  for (comparison_name in comparison_order) {
    key <- paste(comparison_name, "Raw", sep = "__")
    if (!is.null(analysis_store[[key]])) {
      pca_raw_plots[[comparison_name]] <- plot_pca_support(analysis_store[[key]], comparison_name)
      empirical_raw_plots[[comparison_name]] <- plot_empirical_hbfss_support(analysis_store[[key]]$results, comparison_name)
    }
  }

  pca_raw_panel <- arrange_with_one_legend(pca_raw_plots, "PCA: Raw", ncol = length(comparison_order))
  if (!is.null(pca_raw_panel)) {
    save_plot(pca_raw_panel, file.path(figure_dir, "Figure_Manuscript_PCA_Raw_AllComparisons.png"), width = 18.0, height = 5.8)
  }

  empirical_raw_panel <- arrange_with_one_legend(empirical_raw_plots, "Empirical p vs HBFSS: Raw", ncol = length(comparison_order))
  if (!is.null(empirical_raw_panel)) {
    save_plot(empirical_raw_panel, file.path(figure_dir, "Figure_Manuscript_Empirical_HBFSS_Raw_AllComparisons.png"), width = 18.0, height = 5.8)
  }

  for (dataset_prefix in c("Lead", "Rem")) {
    pca_plots <- list()
    empirical_plots <- list()

    for (track in analysis_tracks) {
      for (comparison_name in comparison_order) {
        analysis_label <- paste0(dataset_prefix, "_", track)
        key <- paste(comparison_name, analysis_label, sep = "__")
        if (!is.null(analysis_store[[key]])) {
          short_title <- paste0(comparison_name, "\n", track)
          pca_plots[[paste(comparison_name, track, sep = "_")]] <- plot_pca_support(analysis_store[[key]], short_title)
          empirical_plots[[paste(comparison_name, track, sep = "_")]] <- plot_empirical_hbfss_support(analysis_store[[key]]$results, short_title)
        }
      }
    }

    dataset_file <- if (dataset_prefix == "Lead") "Lead" else "Rem"

    pca_panel <- arrange_with_one_legend(
      pca_plots,
      paste0("PCA: ", ifelse(dataset_prefix == "Lead", "Lead", "Rem")),
      ncol = length(comparison_order)
    )
    if (!is.null(pca_panel)) {
      save_plot(
        pca_panel,
        file.path(figure_dir, paste0("Figure_Manuscript_PCA_", dataset_file, "_AllComparisons_NormEVS_vs_RawEVS.png")),
        width = 18.0,
        height = 9.6
      )
    }

    empirical_panel <- arrange_with_one_legend(
      empirical_plots,
      paste0("Empirical p vs HBFSS: ", ifelse(dataset_prefix == "Lead", "Lead", "Rem")),
      ncol = length(comparison_order)
    )
    if (!is.null(empirical_panel)) {
      save_plot(
        empirical_panel,
        file.path(figure_dir, paste0("Figure_Manuscript_Empirical_HBFSS_", dataset_file, "_AllComparisons_NormEVS_vs_RawEVS.png")),
        width = 18.0,
        height = 9.6
      )
    }
  }

  for (dataset_prefix in c("Raw", "Lead", "Rem")) {
    score_plots <- list()

    if (dataset_prefix == "Raw") {
      for (comparison_name in comparison_order) {
        key <- paste(comparison_name, "Raw", sep = "__")
        if (!is.null(analysis_store[[key]])) {
          score_df <- pca_score_distribution_rows(analysis_store[[key]], comparison_name, "Raw")
          score_plots[[comparison_name]] <- plot_pca_score_distribution_panel(score_df, comparison_name)
        }
      }
      score_title <- "PCA score hist: Raw"
      score_file <- "Raw"
      score_height <- 6.8
    } else {
      for (track in analysis_tracks) {
        for (comparison_name in comparison_order) {
          analysis_label <- paste0(dataset_prefix, "_", track)
          key <- paste(comparison_name, analysis_label, sep = "__")
          if (!is.null(analysis_store[[key]])) {
            score_df <- pca_score_distribution_rows(analysis_store[[key]], comparison_name, analysis_label)
            score_plots[[paste(comparison_name, track, sep = "_")]] <- plot_pca_score_distribution_panel(score_df, paste0(comparison_name, "
", track))
          }
        }
      }
      score_title <- if (dataset_prefix == "Lead") "PCA score hist: Lead" else "PCA score hist: Rem"
      score_file <- dataset_prefix
      score_height <- 10.0
    }

    score_panel <- arrange_with_one_legend(score_plots, score_title, ncol = length(comparison_order))
    if (!is.null(score_panel)) {
      save_plot(
        score_panel,
        file.path(figure_dir, paste0("Figure_Manuscript_PCA_Score_Hist_", score_file, "_AllComparisons", ifelse(score_file == "Raw", "", "_NormEVS_vs_RawEVS"), ".png")),
        width = 18.0,
        height = score_height
      )
    }
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
    save_csv(dplyr::rename(pca_var_df, comp = comparison_name, label = analysis_label, panel = analysis_pretty, set = dataset_group, mode = evs_mode, PC = component, varPct = variance_explained_pct, scoreVar = score_variance), file.path(output_dir, "Table_PCA_Variance_Summary.csv"))

    pca_var_plot <- plot_pca_variance_summary(
      pca_var_df,
      value_col = "variance_explained_pct",
      title = "PCA variance explained",
      y_label = "% var",
      caption = NULL
    )
    if (!is.null(pca_var_plot)) {
      save_plot(pca_var_plot, file.path(figure_dir, "Figure_Manuscript_PCA_Variance_Explained.png"), width = 15.0, height = 9.2)
    }

    pca_score_var_plot <- plot_pca_variance_summary(
      pca_var_df,
      value_col = "score_variance",
      title = "PCA score variance",
      y_label = "Score var",
      caption = NULL
    )
    if (!is.null(pca_score_var_plot)) {
      save_plot(pca_score_var_plot, file.path(figure_dir, "Figure_Manuscript_PCA_Score_Variance.png"), width = 15.0, height = 9.2)
    }

    evs_benefit_plot <- plot_evs_benefit_summary(pca_var_df)
    if (!is.null(evs_benefit_plot)) {
      save_plot(evs_benefit_plot, file.path(figure_dir, "Figure_Manuscript_EVS_Benefit_Summary.png"), width = 10.5, height = 5.8)
    }
  }

  loading_mass_rows <- lapply(evs_store, topn_loading_mass_rows)
  loading_mass_df <- if (length(loading_mass_rows) > 0L) bind_rows(loading_mass_rows) else data.frame()
  if (nrow(loading_mass_df) > 0L) {
    save_csv(dplyr::rename(loading_mass_df, comp = comparison_name, cond = condition, frac = topN_fraction), file.path(output_dir, "Table_TopN_Loading_Mass.csv"))
    loading_mass_plot <- plot_topn_loading_mass(loading_mass_df)
    if (!is.null(loading_mass_plot)) {
      save_plot(loading_mass_plot, file.path(figure_dir, "Figure_Manuscript_TopN_Loading_Mass.png"), width = 10.5, height = 5.8)
    }
  }

  for (dataset_prefix in c("Raw", "Lead", "Rem")) {
    loading_plots <- list()

    if (dataset_prefix == "Raw") {
      for (track in analysis_tracks) {
        for (comparison_name in comparison_order) {
          comparison <- comparison_store[[comparison_name]]
          if (!is.null(comparison)) {
            prof <- loading_profile_from_matrix(comparison$count_matrix, comparison$coldata, comparison_name, "Raw", track)
            loading_plots[[paste(comparison_name, track, sep = "_")]] <- plot_loading_distribution_panel(prof, paste0(comparison_name, "\n", track))
          }
        }
      }
      panel_title <- "Loading: Raw"
      file_stub <- "Raw"
    } else {
      for (track in analysis_tracks) {
        for (comparison_name in comparison_order) {
          comparison <- comparison_store[[comparison_name]]
          evs <- evs_store[[paste(comparison_name, track, sep = "__")]]
          if (!is.null(comparison) && !is.null(evs)) {
            mat <- if (dataset_prefix == "Lead") evs$leading_matrix else evs$remainder_matrix
            prof <- loading_profile_from_matrix(mat, comparison$coldata, comparison_name, paste0(dataset_prefix, "_", track), track)
            loading_plots[[paste(comparison_name, track, sep = "_")]] <- plot_loading_distribution_panel(prof, paste0(comparison_name, "\n", track))
          }
        }
      }
      panel_title <- if (dataset_prefix == "Lead") {
        "Loading: Lead"
      } else {
        "Loading: Rem"
      }
      file_stub <- dataset_prefix
    }

    loading_panel <- arrange_with_one_legend(loading_plots, panel_title, ncol = length(comparison_order))
    if (!is.null(loading_panel)) {
      save_plot(
        loading_panel,
        file.path(figure_dir, paste0("Figure_Manuscript_Loading_Distribution_", file_stub, "_AllComparisons_NormEVS_vs_RawEVS.png")),
        width = 18.0,
        height = ifelse(dataset_prefix == "Raw", 9.6, 9.6)
      )
    }
  }
}

# =============================================================================
# MANIFEST AND COMPLETION
# =============================================================================

exported_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
exported_files <- exported_files[file.info(exported_files)$isdir %in% FALSE]
manifest <- data.frame(
  f = sub(paste0("^", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(exported_files, winslash = "/", mustWork = FALSE)),
  bytes = file.info(exported_files)$size
)
save_csv(manifest, file.path(output_dir, "Table_Export_Manifest.csv"))

cat("\n=====================================================\n")
cat("SEQUENCE pipeline complete.\n")
cat("Output directory: ", output_dir, "\n", sep = "")
cat("Exported files: ", nrow(manifest), "\n", sep = "")
cat("=====================================================\n\n")

if (length(failure_rows) > 0L) {
  stop("One or more comparisons failed. See Table_Failed.csv.")
}

if (nrow(summary_df) == 0L) {
  stop("No successful analyses were completed.")
}

print(summary_df)
