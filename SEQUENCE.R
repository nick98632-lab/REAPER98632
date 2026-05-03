#!/usr/bin/env Rscript

# =============================================================================
# SEQUENCE FINAL SUBMISSION PIPELINE
# DESeq2 + Eigenvector Splitting + empirical-null HC/HBFSS
# Final manuscript version. Complete rewrite, not a patch.
# No legacy changepoint logic and no hidden cutoff fallback.
# PCA score-variance figures compare Original vs Leading Edge vs Remainder using PC1 score variance only; EVS histograms show per-feature PC1 variance contributions.
# Standard effects use DESeq2 Wald BH < 10% with |shrunken LFC| >= 1.
# Strong effects are DESeq2 greaterAbs alternative-hypothesis calls: BH < 10% with |shrunken LFC| >= 1.
# Weak class requires DESeq2 lessAbs support plus |shrunken LFC| < 1 plus HBFSS passage.
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
base_theme_size <- 9
n_top_labels <- 6L
n_top_labels_per_class <- 2L

run_simulation_validation <- TRUE
reset_output_dir <- TRUE
simulation_seed <- 42L
simulation_n_features <- 10000L
simulation_n_samples_per_group <- 5L
simulation_base_mean <- 200
simulation_dispersion_null <- 0.10
simulation_dispersion_de <- 0.15
simulation_de_fractions <- c(0.05, 0.10, 0.20)
simulation_lfc_magnitudes <- c(0.5, 1.0, 2.0)
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
simulation_dir <- file.path(output_dir, "simulation_validation")

if (isTRUE(reset_output_dir) && dir.exists(output_dir)) {
  unlink(output_dir, recursive = TRUE, force = TRUE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(simulation_dir, recursive = TRUE, showWarnings = FALSE)

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

format_compact_number <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_character_, length(x))
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
  out[!finite] <- "NA"
  out
}


save_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE)
  invisible(path)
}

concise_analysis_summary <- function(df) {
  out <- df
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
  hit <- intersect(names(map), names(out))
  names(out)[match(hit, names(out))] <- unname(map[hit])

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
  keep <- intersect(keep, names(out))
  out[, keep, drop = FALSE]
}

concise_evs_summary <- function(df) {
  out <- df
  map <- c(
    comparison_name = "comparison",
    track = "evs_mode",
    top_n_target = "top_n_target",
    treatment_top_n_used = "treatment_top_n",
    control_top_n_used = "control_top_n",
    treatment_loading_cutoff = "treatment_loading_cutoff",
    control_loading_cutoff = "control_loading_cutoff",
    leading_edge_n = "leading_edge_features",
    remainder_n = "remainder_features",
    evs_input = "input_matrix"
  )
  hit <- intersect(names(map), names(out))
  names(out)[match(hit, names(out))] <- unname(map[hit])

  keep <- c(
    "comparison", "evs_mode", "top_n_target",
    "treatment_top_n", "control_top_n",
    "treatment_loading_cutoff", "control_loading_cutoff",
    "leading_edge_features", "remainder_features", "input_matrix"
  )
  keep <- intersect(keep, names(out))
  out[, keep, drop = FALSE]
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

  if (is.null(legend)) return(gridExtra::arrangeGrob(body, ncol = 1, top = title_grob))
  gridExtra::arrangeGrob(body, legend, ncol = 1, heights = c(12, 0.95), top = title_grob)
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
  if (track == "RawEVS") {
    x <- matrix(as.numeric(count_matrix), nrow = nrow(count_matrix), dimnames = dimnames(count_matrix))
    return(log2(x + 1))
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
  pc1_score_variance <- as.numeric(pca$sdev[1]^2)
  pc1_variance_contribution <- (as.numeric(loading)^2) * pc1_score_variance

  rank_df <- data.frame(
    feature_id = names(loading),
    pc1_loading = as.numeric(loading),
    pc1_loading_abs = abs(as.numeric(loading)),
    pc1_score_variance = pc1_score_variance,
    pc1_variance_contribution = pc1_variance_contribution
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
    evs_input = ifelse(track == "NormEVS", "DESeq2 VST counts for PCA/EVS; log2 normalized-count fallback", "log2 raw counts + 1 for RawEVS PCA/EVS")
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

  # HBFSS is reported as a superset of the ordinary DESeq2 standard calls.
  # This makes the summary arithmetic explicit:
  # hbfss_total = standard_effect + hbfss_only.
  df$hbfss_flag <- df$standard_flag | df$hbfss_raw_flag

  df$weak_region_hbfss_flag <- df$lessAbs_alt_flag & df$hbfss_flag
  df$weak_flag <- df$weak_region_hbfss_flag
  df$weak_hbfss_flag <- df$weak_flag

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
    "lessAbs_alt_flag", "weak_region_hbfss_flag", "weak_hbfss_flag", "weak_flag", "hbfss_raw_flag", "hbfss_flag",
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
      fill = "none",
      shape = guide_legend(
        override.aes = list(
          shape = unname(class_shapes[class_levels]),
          color = unname(class_colors[class_levels]),
          fill = unname(class_colors[class_levels]),
          size = rep(3.2, length(class_levels)),
          alpha = rep(1, length(class_levels)),
          stroke = rep(0.70, length(class_levels))
        ),
        nrow = 1,
        byrow = TRUE
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
      size = 1.90,
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

  cutoff_rows <- list()
  if (!is.null(evs_obj$treatment_rank$pca)) {
    cutoff_rows[[length(cutoff_rows) + 1L]] <- data.frame(
      comparison_name = comparison_name,
      track = track,
      condition = "Treatment",
      pc1_loading_cutoff = evs_obj$summary$treatment_loading_cutoff[1],
      pc1_score_variance = as.numeric(evs_obj$treatment_rank$pca$sdev[1]^2),
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(evs_obj$control_rank$pca)) {
    cutoff_rows[[length(cutoff_rows) + 1L]] <- data.frame(
      comparison_name = comparison_name,
      track = track,
      condition = "Control",
      pc1_loading_cutoff = evs_obj$summary$control_loading_cutoff[1],
      pc1_score_variance = as.numeric(evs_obj$control_rank$pca$sdev[1]^2),
      stringsAsFactors = FALSE
    )
  }
  cutoffs <- if (length(cutoff_rows) > 0L) bind_rows(cutoff_rows) else data.frame()
  if (nrow(cutoffs) > 0L) {
    cutoffs$pc1_variance_contribution_cutoff <- (cutoffs$pc1_loading_cutoff^2) * cutoffs$pc1_score_variance
    cutoffs <- cutoffs[rep(seq_len(nrow(cutoffs)), each = 3), , drop = FALSE]
    cutoffs$dataset_group <- rep(c("Original", "Lead", "Remainder"), times = nrow(cutoffs) / 3)
    attr(out, "cutoffs") <- cutoffs
  }

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
      alpha = 0.18,
      linewidth = 0.16,
      boundary = 0
    ) +
    facet_grid(dataset_group ~ comparison_name, scales = "free_y") +
    scale_x_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      labels = function(x) format_compact_number(x, digits = 3)
    ) +
    scale_color_manual(values = c(All = "#7570B3", Control = control_color, Treatment = treatment_color), drop = FALSE, name = NULL) +
    scale_fill_manual(values = c(All = "#7570B3", Control = control_color, Treatment = treatment_color), drop = FALSE, name = NULL) +
    labs(
      title = title_text,
      x = "Per-feature contribution to PC1 score variance",
      y = "Feature frequency",
      caption = "Dashed vertical lines mark the treatment/control top-N EVS cutoff projected onto PC1 variance contribution. Original, leading-edge, and remainder distributions are shown directly; no percent-of-original scaling is used."
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

pca_pc1_evs_plot_df <- function(pca_var_df) {
  if (nrow(pca_var_df) == 0L) return(data.frame())

  pc1 <- pca_var_df[pca_var_df$component == "PC1", , drop = FALSE]
  if (nrow(pc1) == 0L) return(data.frame())

  pc1$dataset_group <- as.character(pc1$dataset_group)
  pc1$evs_mode <- as.character(pc1$evs_mode)

  original_rows <- pc1[pc1$dataset_group == "Original" & pc1$evs_mode == "Original", , drop = FALSE]
  split_rows <- pc1[pc1$dataset_group %in% c("Lead", "Remainder") & pc1$evs_mode %in% c("NormEVS", "RawEVS"), , drop = FALSE]
  if (nrow(split_rows) == 0L || nrow(original_rows) == 0L) return(data.frame())

  out_rows <- list()
  idx <- 0L
  for (mode_name in c("NormEVS", "RawEVS")) {
    mode_split <- split_rows[split_rows$evs_mode == mode_name, , drop = FALSE]
    if (nrow(mode_split) == 0L) next

    for (comparison_name in unique(mode_split$comparison_name)) {
      original_one <- original_rows[original_rows$comparison_name == comparison_name, , drop = FALSE]
      split_one <- mode_split[mode_split$comparison_name == comparison_name, , drop = FALSE]
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
      block$pc1_score_variance_log10 <- log10(pmax(block$pc1_score_variance, 0) + 1)
      block$dataset_group <- factor(as.character(block$dataset_group), levels = c("Original", "Lead", "Remainder"))
      block$evs_mode <- factor(as.character(block$evs_mode), levels = c("NormEVS", "RawEVS"))
      idx <- idx + 1L
      out_rows[[idx]] <- block
    }
  }

  if (length(out_rows) == 0L) return(data.frame())
  bind_rows(out_rows)
}

plot_pc1_evs_variance_absolute <- function(pca_var_df) {
  plot_df <- pca_pc1_evs_plot_df(pca_var_df)
  if (nrow(plot_df) == 0L) return(NULL)

  plot_df$dataset_group <- factor(as.character(plot_df$dataset_group), levels = c("Original", "Lead", "Remainder"))
  plot_df$evs_mode <- factor(as.character(plot_df$evs_mode), levels = c("NormEVS", "RawEVS"))
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
    scale_fill_manual(
      values = c(Original = "#9E9E9E", Lead = "#E31A1C", Remainder = "#1F78B4"),
      drop = FALSE,
      name = NULL
    ) +
    labs(
      title = "PC1 score variance by dataset after EVS",
      x = NULL,
      y = "PC1 score variance, absolute scale",
      caption = "Bars compare the original dataset, EVS leading edge, and EVS remainder directly. PC2 percent-variance is not used in this figure."
    ) +
    manuscript_theme() +
    theme(
      strip.text = element_text(size = base_theme_size - 0.6),
      axis.text.x = element_text(size = base_theme_size - 0.9),
      legend.position = "bottom"
    )
}

plot_pc1_evs_variance_relative <- function(pca_var_df) {
  # Kept only for backward compatibility with older script calls.
  # It now returns the absolute PC1 score-variance plot because the relative 100% plot obscured the remainder.
  plot_pc1_evs_variance_absolute(pca_var_df)
}

plot_pc1_evs_variance <- function(pca_var_df) {
  plot_pc1_evs_variance_absolute(pca_var_df)
}

plot_discovery_counts <- function(summary_df) {
  long_df <- bind_rows(lapply(seq_len(nrow(summary_df)), function(i) {
    row <- summary_df[i, , drop = FALSE]
    data.frame(
      comparison_name = row$comparison_name,
      analysis_label = row$analysis_label,
      dataset_key = row$dataset_key,
      Class = factor(c("Weak", "Strong", "Std", "HBFSS"), levels = c("Weak", "Strong", "Std", "HBFSS")),
      Count = as.numeric(c(row$n_weak, row$n_strong, row$n_standard, row$n_hbfss_display_only))
    )
  }))

  save_csv(dplyr::rename(long_df, comparison = comparison_name, analysis = analysis_label, dataset = dataset_key, class = Class, count = Count), file.path(output_dir, "Counts_Long.csv"))

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
      title = "Significant feature counts",
      x = NULL,
      y = "Count",
      caption = NULL
    ) +
    manuscript_theme() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
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








# =============================================================================
# SIMULATION VALIDATION
# =============================================================================

simulate_sequence_counts <- function(n_features, n_samples, base_mean, disp_null, de_fraction, lfc_magnitude, disp_de) {
  n_de <- round(n_features * de_fraction)
  n_null <- n_features - n_de

  de_direction <- ifelse(seq_len(n_de) <= n_de / 2, 1, -1)
  de_control_mean <- rep(base_mean, n_de)
  de_treatment_mean <- base_mean * 2^(lfc_magnitude * de_direction)

  null_control_mean <- rep(base_mean, n_null)
  null_treatment_mean <- rep(base_mean, n_null)

  generate_nb <- function(mu, dispersion) {
    out <- matrix(0L, nrow = length(mu), ncol = n_samples)
    for (j in seq_len(n_samples)) {
      out[, j] <- rnbinom(length(mu), mu = mu, size = 1 / dispersion)
    }
    out
  }

  counts <- rbind(
    cbind(generate_nb(de_control_mean, disp_de), generate_nb(de_treatment_mean, disp_de)),
    cbind(generate_nb(null_control_mean, disp_null), generate_nb(null_treatment_mean, disp_null))
  )
  storage.mode(counts) <- "integer"

  colnames(counts) <- c(paste0("ctrl_", seq_len(n_samples)), paste0("trt_", seq_len(n_samples)))
  rownames(counts) <- paste0("sim_feature_", seq_len(n_features))

  data.frame(
    feature_id = rownames(counts),
    is_de = c(rep(TRUE, n_de), rep(FALSE, n_null)),
    true_lfc = c(lfc_magnitude * de_direction, rep(0, n_null)),
    stringsAsFactors = FALSE
  ) -> truth

  list(counts = counts, truth = truth)
}

inflate_null_wald_statistics <- function(wald, is_de, inflation_fraction, inflation_sd = 1.5) {
  if (inflation_fraction <= 0) return(wald)

  null_idx <- which(!is_de & is.finite(wald) & !is.na(wald))
  n_inflate <- round(length(null_idx) * inflation_fraction)
  if (n_inflate <= 0L) return(wald)

  target <- sample(null_idx, n_inflate, replace = FALSE)
  wald[target] <- rnorm(n_inflate, mean = 0, sd = inflation_sd)
  wald
}

sequence_simulation_analysis <- function(sim_obj, null_inflation) {
  counts <- sim_obj$counts
  n_samp <- ncol(counts) / 2

  coldata <- data.frame(
    condition = factor(c(rep("untrt", n_samp), rep("trt", n_samp)), levels = c("untrt", "trt")),
    row.names = colnames(counts)
  )

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, betaPrior = FALSE, quiet = TRUE)

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
  shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)

  df <- as.data.frame(standard)
  df$feature_id <- rownames(df)
  df$lfc_shrunk <- as.data.frame(shrunk)$log2FoldChange
  df$greaterAbs_padj <- as.data.frame(strong)$padj
  df$lessAbs_padj <- as.data.frame(weak)$padj
  df <- left_join(df, sim_obj$truth, by = "feature_id")

  wald <- inflate_null_wald_statistics(df$stat, df$is_de, null_inflation)
  empirical <- empirical_null(wald, "simulation")
  df$empirical_p <- empirical$empirical_p
  df$empirical_bh <- empirical$empirical_bh

  df$empirical_p_calc <- clip_probability(df$empirical_p, floor = calculation_probability_floor)
  df$neglog10_empirical_p_calc <- -log10(df$empirical_p_calc)

  hc_p <- hc_threshold(df$empirical_p)
  hbfss_cutoff <- if (is.na(hc_p)) NA_real_ else -log10(hc_p) * lfc_boundary

  df$hc_pass <- !is.na(hc_p) &
    !is.na(df$empirical_p) &
    is.finite(df$empirical_p) &
    df$empirical_p <= hc_p

  df$HBFSS <- abs(df$lfc_shrunk) * df$neglog10_empirical_p_calc
  df$hbfss_raw_sig <- !is.na(hbfss_cutoff) &
    !is.na(df$HBFSS) &
    is.finite(df$HBFSS) &
    df$HBFSS >= hbfss_cutoff &
    df$hc_pass

  df$standard_sig <- !is.na(df$padj) &
    df$padj < alpha_level &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= lfc_boundary

  df$empirical_bh_sig <- !is.na(df$empirical_bh) &
    df$empirical_bh < alpha_level &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= lfc_boundary

  df$hbfss_sig <- df$standard_sig | df$hbfss_raw_sig

  df$lessAbs_sig <- !is.na(df$lessAbs_padj) &
    df$lessAbs_padj < weak_alpha_level &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) < lfc_boundary

  df$weak_region_hbfss_sig <- df$lessAbs_sig & df$hbfss_sig

  list(results = df, hc_p = hc_p, hbfss_cutoff = hbfss_cutoff)
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
  fdr <- if ((tp + fp) == 0) NA_real_ else fp / (tp + fp)
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA_real_ else {
    2 * precision * recall / (precision + recall)
  }

  data.frame(tp = tp, fp = fp, fn = fn, tn = tn, precision = precision, recall = recall, f1 = f1, fdr = fdr)
}

simulation_method_summary <- function(results_df, metric) {
  rows <- list(
    data.frame(
      de_fraction = results_df$de_fraction,
      lfc_magnitude = results_df$lfc_magnitude,
      null_inflation = results_df$null_inflation,
      method = "DESeq2 BH",
      value = results_df[[paste0("standard_", metric)]]
    ),
    data.frame(
      de_fraction = results_df$de_fraction,
      lfc_magnitude = results_df$lfc_magnitude,
      null_inflation = results_df$null_inflation,
      method = "Empirical BH",
      value = results_df[[paste0("empirical_bh_", metric)]]
    ),
    data.frame(
      de_fraction = results_df$de_fraction,
      lfc_magnitude = results_df$lfc_magnitude,
      null_inflation = results_df$null_inflation,
      method = "HBFSS",
      value = results_df[[paste0("hbfss_", metric)]]
    )
  )
  bind_rows(rows) %>%
    group_by(de_fraction, lfc_magnitude, null_inflation, method) %>%
    summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      lfc_label = factor(paste0("LFC ", lfc_magnitude), levels = paste0("LFC ", simulation_lfc_magnitudes)),
      de_label = paste0(de_fraction * 100, "% DE"),
      inflation_label = paste0(null_inflation * 100, "% inflation"),
      method = factor(method, levels = c("DESeq2 BH", "Empirical BH", "HBFSS"))
    )
}

plot_simulation_metric <- function(summary_df, title, y_label, alpha_line = FALSE) {
  method_colors <- c("DESeq2 BH" = "#999999", "Empirical BH" = treatment_color, "HBFSS" = class_colors[["HBFSS"]])

  p <- ggplot(summary_df, aes(lfc_label, mean, color = method, group = method)) +
    geom_line(linewidth = 0.55, position = position_dodge(width = 0.30)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.30)) +
    geom_errorbar(
      aes(ymin = pmax(mean - sd, 0), ymax = mean + sd),
      width = 0.12,
      linewidth = 0.25,
      position = position_dodge(width = 0.30)
    ) +
    facet_grid(inflation_label ~ de_label) +
    scale_color_manual(values = method_colors, drop = FALSE, name = NULL) +
    labs(title = title, x = NULL, y = y_label, caption = NULL) +
    manuscript_theme() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  if (alpha_line) {
    p <- p + geom_hline(yintercept = alpha_level, linetype = "dashed", linewidth = 0.30, color = "grey35")
  }

  p
}

plot_simulation_weak_recall <- function(results_df) {
  weak_df <- results_df[results_df$lfc_magnitude < lfc_boundary, , drop = FALSE]
  if (nrow(weak_df) == 0L) return(NULL)

  long <- bind_rows(
    data.frame(
      de_fraction = weak_df$de_fraction,
      null_inflation = weak_df$null_inflation,
      method = "DESeq2 BH",
      recall = weak_df$standard_weak_recall
    ),
    data.frame(
      de_fraction = weak_df$de_fraction,
      null_inflation = weak_df$null_inflation,
      method = "HBFSS weak-region",
      recall = weak_df$hbfss_weak_recall
    )
  ) %>%
    mutate(
      de_label = paste0(de_fraction * 100, "% DE"),
      inflation_label = paste0(null_inflation * 100, "% inflation"),
      method = factor(method, levels = c("DESeq2 BH", "HBFSS weak-region"))
    )

  ggplot(long, aes(de_label, recall, color = method, group = method)) +
    geom_point(size = 1.7, position = position_dodge(width = 0.25)) +
    geom_line(linewidth = 0.55, position = position_dodge(width = 0.25)) +
    facet_wrap(~ inflation_label) +
    scale_color_manual(values = c("DESeq2 BH" = "#999999", "HBFSS weak-region" = class_colors[["Weak"]]), drop = FALSE, name = NULL) +
    labs(title = "Weak-effect recall", x = NULL, y = "Recall", caption = NULL) +
    manuscript_theme()
}

plot_simulation_cutoffs <- function(results_df) {
  plot_df <- results_df[is.finite(results_df$hbfss_cutoff) & !is.na(results_df$hbfss_cutoff), , drop = FALSE]
  if (nrow(plot_df) == 0L) return(NULL)

  plot_df <- plot_df %>%
    mutate(
      lfc_label = factor(paste0("LFC ", lfc_magnitude), levels = paste0("LFC ", simulation_lfc_magnitudes)),
      de_label = paste0(de_fraction * 100, "% DE")
    )

  ggplot(plot_df, aes(de_label, hbfss_cutoff, fill = de_label)) +
    geom_boxplot(outlier.size = 0.35, linewidth = 0.25) +
    facet_wrap(~ lfc_label) +
    scale_fill_manual(values = scales::grey_pal(start = 0.35, end = 0.75)(length(unique(plot_df$de_label))), name = NULL) +
    labs(title = "HBFSS cutoff stability", x = NULL, y = "HBFSS cutoff", caption = NULL) +
    manuscript_theme() +
    theme(legend.position = "none")
}


simulation_long_metric_rows <- function(results_df, metric) {
  bind_rows(
    data.frame(
      de_fraction = results_df$de_fraction,
      lfc_magnitude = results_df$lfc_magnitude,
      null_inflation = results_df$null_inflation,
      replicate = results_df$replicate,
      method = "DESeq2 BH",
      value = results_df[[paste0("standard_", metric)]],
      stringsAsFactors = FALSE
    ),
    data.frame(
      de_fraction = results_df$de_fraction,
      lfc_magnitude = results_df$lfc_magnitude,
      null_inflation = results_df$null_inflation,
      replicate = results_df$replicate,
      method = "Empirical BH",
      value = results_df[[paste0("empirical_bh_", metric)]],
      stringsAsFactors = FALSE
    ),
    data.frame(
      de_fraction = results_df$de_fraction,
      lfc_magnitude = results_df$lfc_magnitude,
      null_inflation = results_df$null_inflation,
      replicate = results_df$replicate,
      method = "HBFSS",
      value = results_df[[paste0("hbfss_", metric)]],
      stringsAsFactors = FALSE
    )
  ) %>%
    mutate(
      method = factor(method, levels = c("DESeq2 BH", "Empirical BH", "HBFSS")),
      lfc_label = factor(
        paste0("LFC = ", lfc_magnitude),
        levels = paste0("LFC = ", simulation_lfc_magnitudes)
      ),
      de_label = paste0(de_fraction * 100, "% DE"),
      inflation_label = paste0(null_inflation * 100, "% null inflation")
    )
}

plot_simulation_metric_boxplot <- function(results_df, metric, title, y_label, alpha_line = FALSE) {
  long <- simulation_long_metric_rows(results_df, metric)
  method_colors <- c("DESeq2 BH" = "#999999", "Empirical BH" = treatment_color, "HBFSS" = class_colors[["HBFSS"]])

  p <- ggplot(long, aes(method, value, fill = method)) +
    geom_boxplot(outlier.size = 0.45, width = 0.62, linewidth = 0.22, na.rm = TRUE) +
    facet_grid(inflation_label + de_label ~ lfc_label) +
    scale_fill_manual(values = method_colors, drop = FALSE, name = NULL) +
    labs(title = title, x = NULL, y = y_label, caption = NULL) +
    manuscript_theme() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(size = base_theme_size - 0.8),
      legend.position = "bottom"
    )

  if (alpha_line) {
    p <- p + geom_hline(yintercept = alpha_level, linetype = "dashed", linewidth = 0.30, color = "grey35")
  }

  p
}

run_sequence_simulation_validation <- function() {
  set.seed(simulation_seed)
  dir.create(simulation_dir, recursive = TRUE, showWarnings = FALSE)

  old_sim_files <- list.files(simulation_dir, full.names = TRUE)
  if (length(old_sim_files) > 0L) unlink(old_sim_files, recursive = TRUE, force = TRUE)

  sim_grid <- expand.grid(
    de_fraction = simulation_de_fractions,
    lfc_magnitude = simulation_lfc_magnitudes,
    null_inflation = simulation_null_inflation,
    stringsAsFactors = FALSE
  )

  total_runs <- nrow(sim_grid) * simulation_n_reps
  message("Running simulation validation: ", total_runs, " runs")

  result_rows <- list()
  run_index <- 0L

  for (grid_i in seq_len(nrow(sim_grid))) {
    grid_row <- sim_grid[grid_i, , drop = FALSE]

    for (rep_i in seq_len(simulation_n_reps)) {
      run_index <- run_index + 1L
      if (run_index %% 10L == 0L) message("Simulation progress: ", run_index, "/", total_runs)

      sim_obj <- simulate_sequence_counts(
        n_features = simulation_n_features,
        n_samples = simulation_n_samples_per_group,
        base_mean = simulation_base_mean,
        disp_null = simulation_dispersion_null,
        de_fraction = grid_row$de_fraction,
        lfc_magnitude = grid_row$lfc_magnitude,
        disp_de = simulation_dispersion_de
      )

      out <- tryCatch(
        sequence_simulation_analysis(sim_obj, null_inflation = grid_row$null_inflation),
        error = function(e) {
          message("Simulation failed: ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(out)) next

      df <- out$results
      true_de <- df$is_de
      true_weak <- df$is_de & abs(df$true_lfc) < lfc_boundary

      standard_all <- simulation_metrics(df$standard_sig, true_de)
      empirical_bh_all <- simulation_metrics(df$empirical_bh_sig, true_de)
      hbfss_all <- simulation_metrics(df$hbfss_sig, true_de)
      standard_weak <- simulation_metrics(df$standard_sig, true_weak)
      hbfss_weak <- simulation_metrics(df$weak_region_hbfss_sig, true_weak)

      result_rows[[length(result_rows) + 1L]] <- data.frame(
        de_fraction = grid_row$de_fraction,
        lfc_magnitude = grid_row$lfc_magnitude,
        null_inflation = grid_row$null_inflation,
        replicate = rep_i,
        hc_p_threshold = out$hc_p,
        hbfss_cutoff = out$hbfss_cutoff,
        standard_precision = standard_all$precision,
        standard_recall = standard_all$recall,
        standard_f1 = standard_all$f1,
        standard_fdr = standard_all$fdr,
        empirical_bh_precision = empirical_bh_all$precision,
        empirical_bh_recall = empirical_bh_all$recall,
        empirical_bh_f1 = empirical_bh_all$f1,
        empirical_bh_fdr = empirical_bh_all$fdr,
        hbfss_precision = hbfss_all$precision,
        hbfss_recall = hbfss_all$recall,
        hbfss_f1 = hbfss_all$f1,
        hbfss_fdr = hbfss_all$fdr,
        standard_weak_recall = standard_weak$recall,
        hbfss_weak_recall = hbfss_weak$recall
      )
    }
  }

  simulation_results <- if (length(result_rows) > 0L) bind_rows(result_rows) else data.frame()
  if (nrow(simulation_results) == 0L) stop("Simulation validation produced no successful runs.")

  save_csv(simulation_results, file.path(simulation_dir, "Simulation_Raw.csv"))

  simulation_summary <- simulation_results %>%
    group_by(de_fraction, lfc_magnitude, null_inflation) %>%
    summarise(
      n_replicates = n(),
      standard_f1_mean = mean(standard_f1, na.rm = TRUE),
      empirical_bh_f1_mean = mean(empirical_bh_f1, na.rm = TRUE),
      hbfss_f1_mean = mean(hbfss_f1, na.rm = TRUE),
      standard_fdr_mean = mean(standard_fdr, na.rm = TRUE),
      empirical_bh_fdr_mean = mean(empirical_bh_fdr, na.rm = TRUE),
      hbfss_fdr_mean = mean(hbfss_fdr, na.rm = TRUE),
      standard_weak_recall_mean = mean(standard_weak_recall, na.rm = TRUE),
      hbfss_weak_recall_mean = mean(hbfss_weak_recall, na.rm = TRUE),
      hc_p_threshold_median = median(hc_p_threshold, na.rm = TRUE),
      hbfss_cutoff_median = median(hbfss_cutoff, na.rm = TRUE),
      .groups = "drop"
    )

  save_csv(simulation_summary, file.path(simulation_dir, "Simulation_Summary.csv"))

  save_csv(simulation_long_metric_rows(simulation_results, "f1"), file.path(simulation_dir, "Simulation_F1_Long.csv"))
  save_csv(simulation_long_metric_rows(simulation_results, "fdr"), file.path(simulation_dir, "Simulation_FDR_Long.csv"))

  save_plot(
    plot_simulation_metric_boxplot(simulation_results, "f1", "Simulation F1 score by method", "F1", alpha_line = FALSE),
    file.path(simulation_dir, "Simulation_F1_Boxplot.png"),
    width = 11.8,
    height = 9.2
  )

  save_plot(
    plot_simulation_metric_boxplot(simulation_results, "fdr", "Simulation observed FDR by method", "FDR", alpha_line = TRUE),
    file.path(simulation_dir, "Simulation_FDR_Boxplot.png"),
    width = 11.8,
    height = 9.2
  )

  weak_plot <- plot_simulation_weak_recall(simulation_results)
  if (!is.null(weak_plot)) {
    save_plot(weak_plot, file.path(simulation_dir, "Simulation_WeakRecall.png"), width = 8.6, height = 5.4)
  }

  cutoff_plot <- plot_simulation_cutoffs(simulation_results)
  if (!is.null(cutoff_plot)) {
    save_plot(cutoff_plot, file.path(simulation_dir, "Simulation_Cutoffs.png"), width = 9.0, height = 5.4)
  }

  writeLines(capture.output(sessionInfo()), file.path(simulation_dir, "SessionInfo_Simulation.txt"))
  invisible(simulation_summary)
}


# =============================================================================
# PIPELINE EXECUTION
# =============================================================================

old_figures <- list.files(figure_dir, pattern = "\\.png$", full.names = TRUE)
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
if (nrow(summary_df) > 0L) save_csv(concise_analysis_summary(summary_df), file.path(output_dir, "Summary_Overall.csv"))

if (length(failure_rows) > 0L) {
  failure_df <- bind_rows(failure_rows)
  save_csv(failure_df, file.path(output_dir, "Failures.csv"))
}

# =============================================================================
# CURATED FIGURE EXPORT
# =============================================================================

comparison_order <- comparison_table$comparison_name

if (nrow(summary_df) > 0L) {
  raw_plots <- list()
  for (comparison_name in comparison_order) {
    key <- paste(comparison_name, "Raw", sep = "__")
    if (!is.null(analysis_store[[key]])) {
      raw_plots[[comparison_name]] <- plot_volcano(analysis_store[[key]]$results, comparison_name)
    }
  }
  raw_panel <- arrange_with_one_legend(raw_plots, "Volcano: Raw", ncol = length(comparison_order))
  if (!is.null(raw_panel)) save_plot(raw_panel, file.path(figure_dir, "Volcano_Raw.png"), width = 18.0, height = 5.8)

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

    panel_title <- if (dataset_prefix == "Lead") "Volcano: Lead" else "Volcano: Rem"
    panel_file <- if (dataset_prefix == "Lead") "Volcano_Lead.png" else "Volcano_Rem.png"
    panel <- arrange_with_one_legend(plots, panel_title, ncol = length(comparison_order))
    if (!is.null(panel)) save_plot(panel, file.path(figure_dir, panel_file), width = 18.0, height = 9.6)
  }

  count_plot <- plot_discovery_counts(summary_df)
  if (!is.null(count_plot)) save_plot(count_plot, file.path(figure_dir, "Counts.png"), width = 14.0, height = 8.0)

  pca_raw_plots <- list()
  for (comparison_name in comparison_order) {
    key <- paste(comparison_name, "Raw", sep = "__")
    if (!is.null(analysis_store[[key]])) {
      pca_raw_plots[[comparison_name]] <- plot_pca_support(analysis_store[[key]], comparison_name)
    }
  }
  pca_raw_panel <- arrange_with_one_legend(pca_raw_plots, "PCA: Raw", ncol = length(comparison_order))
  if (!is.null(pca_raw_panel)) save_plot(pca_raw_panel, file.path(figure_dir, "PCA_Raw.png"), width = 18.0, height = 5.8)

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

    panel_title <- if (dataset_prefix == "Lead") "PCA: Lead" else "PCA: Rem"
    panel_file <- if (dataset_prefix == "Lead") "PCA_Lead.png" else "PCA_Rem.png"
    pca_panel <- arrange_with_one_legend(pca_plots, panel_title, ncol = length(comparison_order))
    if (!is.null(pca_panel)) save_plot(pca_panel, file.path(figure_dir, panel_file), width = 18.0, height = 9.4)
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
      dplyr::rename(
        pca_var_df,
        comparison = comparison_name,
        analysis = analysis_label,
        panel = analysis_pretty,
        dataset = dataset_group,
        evs_mode = evs_mode,
        component = component,
        variance_explained_pct = variance_explained_pct,
        score_variance = score_variance
      ),
      file.path(output_dir, "Summary_PCA.csv")
    )

    pca_pc1_compare_df <- pca_pc1_evs_plot_df(pca_var_df)
    if (nrow(pca_pc1_compare_df) > 0L) {
      save_csv(
        dplyr::rename(
          pca_pc1_compare_df,
          comparison = comparison_name,
          analysis = analysis_label,
          panel = analysis_pretty,
          dataset = dataset_group,
          evs_mode = evs_mode,
          pc1_score_variance = pc1_score_variance,
          original_pc1_score_variance = original_pc1_score_variance,
          pc1_score_variance_ratio_to_original = pc1_score_variance_ratio_to_original,
          pc1_score_variance_percent_of_original = pc1_score_variance_percent_of_original
        ),
        file.path(output_dir, "Summary_PCA_PC1_ScoreVariance_EVS.csv")
      )
    }

    pc1_var_plot <- plot_pc1_evs_variance_absolute(pca_var_df)
    if (!is.null(pc1_var_plot)) {
      save_plot(pc1_var_plot, file.path(figure_dir, "PCA_Variance.png"), width = 15.8, height = 8.8)
      save_plot(pc1_var_plot, file.path(figure_dir, "PC1_Variance_EVS.png"), width = 15.8, height = 8.8)
    }
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
        load_rows[[paste(comparison_name, track, sep = "__")]] <- tmp
        cut <- attr(tmp, "cutoffs")
        if (!is.null(cut) && nrow(cut) > 0L) cutoff_rows[[paste(comparison_name, track, sep = "__")]] <- cut
      }
    }
  }
  load_df <- if (length(load_rows) > 0L) bind_rows(load_rows) else data.frame()
  cutoff_df <- if (length(cutoff_rows) > 0L) bind_rows(cutoff_rows) else data.frame()
  if (nrow(load_df) > 0L) {
    attr(load_df, "cutoffs") <- cutoff_df
    save_csv(
      dplyr::rename(
        load_df,
        comparison = comparison_name,
        evs_mode = track,
        dataset = dataset_group,
        scope = condition,
        abs_pc1_loading = pc1_loading_abs,
        pc1_score_variance = pc1_score_variance,
        pc1_variance_contribution = pc1_variance_contribution
      ),
      file.path(output_dir, "Summary_PC1VarianceDistributions.csv")
    )
    if (nrow(cutoff_df) > 0L) {
      save_csv(
        dplyr::rename(
          cutoff_df,
          comparison = comparison_name,
          evs_mode = track,
          dataset = dataset_group,
          scope = condition,
          pc1_loading_cutoff = pc1_loading_cutoff,
          pc1_score_variance = pc1_score_variance,
          pc1_variance_contribution_cutoff = pc1_variance_contribution_cutoff
        ),
        file.path(output_dir, "Summary_PC1VarianceCutoffs.csv")
      )
    }
    load_plot <- plot_loading_histograms(load_df, "PC1 variance-contribution frequency after EVS")
    if (!is.null(load_plot)) {
      save_plot(load_plot, file.path(figure_dir, "EVS_LoadingRank.png"), width = 15.8, height = 10.4)
      save_plot(load_plot, file.path(figure_dir, "PC1_Variance_Distributions.png"), width = 15.8, height = 10.4)
    }
  }
}

if (isTRUE(run_simulation_validation)) {
  simulation_summary <- run_sequence_simulation_validation()
  print(simulation_summary)
}


write_sequence_methods <- function() {
  methods_lines <- c(
    "# SEQUENCE methods summary",
    "",
    "## Input and preprocessing",
    "Raw read-count matrices are imported for each RT/ZT comparison. Raw integer counts are retained for DESeq2 differential testing. PCA/EVS calculations are performed on variance-stabilized expression for NormEVS, with log2(normalized counts + 1) as a fallback if VST fails. RawEVS uses log2(raw counts + 1).",
    "",
    "## Eigenvector splitting",
    "For each comparison and EVS mode, PC1 is calculated separately in treatment and control samples. Features are ranked by absolute PC1 loading. The leading edge is the union of the top-N treatment-ranked and top-N control-ranked features. The remainder is every feature not in that union. No AIC/BIC or legacy changepoint fallback is used.",
    "",
    "## Differential expression",
    "DESeq2 is run with design ~ condition and reference level untrt. Standard effects are BH-adjusted Wald results with padj < alpha and |apeglm-shrunken log2 fold-change| >= 1. Strong effects are DESeq2 greaterAbs alternative-hypothesis calls at the same LFC boundary. Weak HBFSS-region effects require DESeq2 lessAbs support, |shrunken LFC| < 1, and HBFSS significance.",
    "",
    "## Empirical null and HBFSS",
    "The Wald statistic distribution is passed to fdrtool using statistic = normal. The fndr cutoff method is attempted first, with pct0 fallback. Empirical p-values are clipped only for numerical stability. Higher criticism supplies the dataset-specific empirical-p threshold unless the returned threshold is invalid or near 1. The HBFSS cutoff is -log10(HC p-threshold) multiplied by the LFC boundary. HBFSS is reported as a superset of standard DESeq2 calls so hbfss_total = standard_effect + hbfss_only by construction.",
    "",
    "## PCA variance figures",
    "The PCA variance figure uses PC1 score variance only. It directly compares Original, Lead, and Remainder datasets on an absolute pseudo-log scale, avoiding the prior percent-of-original plot that compressed the remainder to near zero. The PC1 variance-distribution figure shows per-feature contribution to PC1 score variance for All samples, Control samples, and Treatment samples, with dashed treatment/control EVS top-N cutoff lines.",
    "",
    "## Simulation validation",
    "The simulation validation generates negative-binomial count data across DE fractions, LFC magnitudes, and null-inflation settings. It runs the same DESeq2/apeglm/fdrtool/HC/HBFSS logic and exports raw results, summary metrics, F1, FDR, weak-effect recall, and cutoff-stability figures."
  )
  writeLines(methods_lines, file.path(output_dir, "METHODS_SEQUENCE_PIPELINE.md"))
}

write_sequence_methods()

# =============================================================================
# MANIFEST AND COMPLETION
# =============================================================================

exported_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
exported_files <- exported_files[file.info(exported_files)$isdir %in% FALSE]
manifest <- data.frame(
  file = sub(paste0("^", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(exported_files, winslash = "/", mustWork = FALSE)),
  size_bytes = file.info(exported_files)$size
)
save_csv(manifest, file.path(output_dir, "Manifest.csv"))

cat("\n=====================================================\n")
cat("SEQUENCE pipeline complete.\n")
cat("Output directory: ", output_dir, "\n", sep = "")
cat("Exported files: ", nrow(manifest), "\n", sep = "")
cat("=====================================================\n\n")

if (length(failure_rows) > 0L) {
  stop("One or more comparisons failed. See Failures.csv.")
}

if (nrow(summary_df) == 0L) {
  stop("No successful analyses were completed.")
}

print(summary_df)
