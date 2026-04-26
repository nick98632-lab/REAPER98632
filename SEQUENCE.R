#!/usr/bin/env Rscript

# =============================================================================
# FINAL MANUSCRIPT PIPELINE
# EVS + DESeq2 + empirical-null calibration + higher criticism + HBFSS
# =============================================================================
#
# This script is written as a manuscript supplement and a runnable analysis
# pipeline. It is intentionally verbose. The comments are part of the scientific
# record: they define the analytic choices, the mathematical thresholds, the
# figure logic, and the export structure.
#
# -----------------------------------------------------------------------------
# STUDY UNIT
# -----------------------------------------------------------------------------
# The unit of analysis is the PAS feature represented by OrigID in the WTTS-Seq
# count matrix. Gene Symbol is retained only as annotation. The pipeline does
# not collapse PAS features to genes before modeling.
#
# -----------------------------------------------------------------------------
# PAIRWISE COMPARISONS
# -----------------------------------------------------------------------------
# Four matched circadian comparisons are evaluated:
#
#   RT0_ZT6
#   RT2_ZT8
#   RT4_ZT10
#   RT8_ZT14
#
# Each comparison uses the same two-level condition design:
#
#   design = ~ condition
#
# where:
#
#   untrt = control arm
#   trt   = treatment arm
#
# -----------------------------------------------------------------------------
# EVS TRACKS
# -----------------------------------------------------------------------------
# The script deliberately evaluates two eigenvector-splitting tracks:
#
#   1. NormEVS
#      PCA / PC1-loading ranking is performed after DESeq2 size-factor
#      normalization of the full comparison matrix.
#
#   2. RawEVS
#      PCA / PC1-loading ranking is performed directly on raw counts, without
#      normalization before eigenvector splitting.
#
# For each track, the downstream DESeq2 analysis is still run independently on
# each resulting dataset. Thus the two tracks compare the effect of the EVS
# preprocessing/ranking input, not the downstream DESeq2 modeling framework.
#
# -----------------------------------------------------------------------------
# DATASETS PRODUCED PER COMPARISON AND PER EVS TRACK
# -----------------------------------------------------------------------------
#
#   Raw dataset:
#     the full comparison count matrix.
#
#   Leading-edge dataset:
#     the union of the top-N PAS features by absolute PC1 loading from each
#     condition-specific PCA.
#
#   Remainder dataset:
#     all PAS features not included in the leading-edge union.
#
# -----------------------------------------------------------------------------
# DESEQ2 AND EMPIRICAL-NULL CALIBRATION
# -----------------------------------------------------------------------------
# DESeq2 is fit independently to Raw, Lead, and Rem datasets for each comparison
# and EVS track. The finite DESeq2 Wald statistics are passed to fdrtool under a
# normal empirical-null model. fdrtool returns empirical-null p-values, q-values,
# and local FDR values. Higher criticism is applied directly to sorted
# empirical-null p-values:
#
#   hc_p_threshold_dataset = fdrtool::hc.thresh(sort(empirical_p))
#
# This threshold is used as a dataset-specific evidence threshold on volcano
# plots. Nothing below this empirical-null HC threshold is colored as a
# manuscript-positive category.
#
# -----------------------------------------------------------------------------
# HBFSS DEFINITION
# -----------------------------------------------------------------------------
# HBFSS combines empirical-null evidence and apeglm-shrunken effect size:
#
#   HBFSS_i = | lfc_shrunk_i * log10(empirical_p_i) |
#
# Since -log10(empirical_p_i) is the y-axis of the manuscript volcano plots,
# this can also be interpreted geometrically as:
#
#   HBFSS_i = | lfc_shrunk_i | * [-log10(empirical_p_i)]
#
# The dataset-specific HBFSS threshold is:
#
#   hbfss_threshold_dataset = abs(log10(hc_p_threshold_dataset)) * lfc_boundary
#
# On a volcano plot, the HBFSS decision boundary is:
#
#   y = hbfss_threshold_dataset / |x|
#
# where x is the shrunken log2 fold change and y is -log10(empirical_p).
# The curve is drawn only inside the plotted coordinate range to avoid a
# nonprofessional vertical asymptote running off the graph near x = 0.
#
# -----------------------------------------------------------------------------
# MATHEMATICAL EFFECT DEFINITIONS
# -----------------------------------------------------------------------------
# Let beta denote the DESeq2 log2 fold change and let c = lfc_boundary.
# In this script c is set to 1.0, corresponding to a two-fold change.
#
# Strong composite-null hypothesis, abbreviated Strong CNH:
#
#   H0,strong : |beta| <= c
#   HA,strong : |beta| >  c
#
# In DESeq2 this is evaluated with:
#
#   altHypothesis = "greaterAbs"
#   lfcThreshold  = c
#
# A PAS is plotted as Strong CNH only if all of the following hold:
#
#   resGA_padj < alpha_level
#   |lfc_shrunk| >= c
#   empirical_p <= hc_p_threshold_dataset
#
# Weak composite-null hypothesis, abbreviated Weak CNH:
#
#   H0,weak : |beta| >= c
#   HA,weak : |beta| <  c
#
# In DESeq2 this is evaluated with:
#
#   altHypothesis = "lessAbs"
#   lfcThreshold  = c
#
# A PAS is plotted as Weak CNH only if all of the following hold:
#
#   resLA_padj < alpha_level
#   |lfc_shrunk| < c
#   empirical_p <= hc_p_threshold_dataset
#
# Standard DESeq2:
#
# A PAS is plotted as Standard only if all of the following hold:
#
#   padj < alpha_level
#   |raw log2FoldChange| >= c
#   |lfc_shrunk| >= c
#   empirical_p <= hc_p_threshold_dataset
#
# HBFSS:
#
# A PAS is plotted as HBFSS only if all of the following hold:
#
#   HBFSS >= hbfss_threshold_dataset
#   empirical_p <= hc_p_threshold_dataset
#
# -----------------------------------------------------------------------------
# FINAL VOLCANO DISPLAY RULES
# -----------------------------------------------------------------------------
# The plotted classes are exactly:
#
#   Background
#   Weak CNH
#   Strong CNH
#   Standard
#   HBFSS
#
# There is one legend for each volcano panel. Overlap counts are written in the
# plot caption and on-panel summary text. Overlap is not a separate color class.
#
# Priority for coloring is:
#
#   HBFSS > Standard > Strong CNH > Weak CNH > Background
#
# This preserves a single clean legend while still reporting overlap counts.
#
# -----------------------------------------------------------------------------
# EXPORT RULES
# -----------------------------------------------------------------------------
# All tables begin with Table_.
# All figures begin with Figure_.
# Files are written into:
#
#   exports/manuscript_final_clean
#
# This directory is inside the GitHub repository when the script is run from the
# repository or from /root/REAPER98632. At the end, the script writes a manifest
# and stages the exports with git add when a git repository is detected.
# Committing and pushing are available but disabled by default.
#
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
})

options(stringsAsFactors = FALSE)

# =============================================================================
# SECTION 1 OF 5
# USER SETTINGS, METADATA, PATHS, AND GENERAL HELPERS
# =============================================================================

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------

count_file_candidates <- c(
  "WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  file.path("data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"),
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
)

run_twas_overlap <- FALSE

twas_file_candidates <- c(
  "3aTWAS_genes_of_11_brain_disorders.csv",
  file.path("data", "3aTWAS_genes_of_11_brain_disorders.csv"),
  "/root/REAPER98632/data/3aTWAS_genes_of_11_brain_disorders.csv"
)

alpha_level <- 0.10

lfc_boundary <- 1.0

top_n_target <- 5000L

figure_dpi <- 320

base_theme_size <- 10

n_top_labels_volcano <- 18

# Git export behavior.
# The script always writes exports into the repository exports folder when it
# can detect the repository. Staging is enabled by default. Commit and push are
# disabled by default so that the user can review output before publishing.
GIT_STAGE_EXPORTS <- TRUE
GIT_COMMIT_EXPORTS <- FALSE
GIT_PUSH_EXPORTS <- FALSE
GIT_COMMIT_MESSAGE <- "Refresh final manuscript EVS HBFSS exports"

# -----------------------------------------------------------------------------
# Palette
# -----------------------------------------------------------------------------

plot_palette <- list(
  background = "#BDBDBD",
  threshold = "#A65628",
  hc = "#A65628",
  hbfss_line = "#6A3D9A",
  weak = "#4EA3F1",
  strong = "#1F78B4",
  standard = "#33A02C",
  hbfss = "#6A3D9A",
  control = "#4D4D4D",
  treatment = "#1F78B4",
  histogram = "#969696"
)

# -----------------------------------------------------------------------------
# Short names used in file exports
# -----------------------------------------------------------------------------

dataset_short <- c(
  raw_dataset = "Raw",
  leading_edge_dataset = "Lead",
  remainder_dataset = "Rem"
)

track_short <- c(
  normalized_evs = "NormEVS",
  raw_evs = "RawEVS"
)

dataset_key_order <- c(
  "raw_dataset",
  "leading_edge_dataset",
  "remainder_dataset"
)

dataset_key_labels <- c(
  raw_dataset = "Original dataset",
  leading_edge_dataset = "Leading-edge dataset",
  remainder_dataset = "Remainder dataset"
)

# -----------------------------------------------------------------------------
# Embedded sample metadata
# -----------------------------------------------------------------------------

meta_all <- data.frame(
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
    "treatment", "treatment", "treatment", "treatment", "treatment",
    "control", "control", "control", "control", "control",
    "treatment", "treatment", "treatment", "treatment", "treatment",
    "control", "control", "control", "control", "control",
    "treatment", "treatment", "treatment", "treatment", "treatment",
    "control", "control", "control", "control", "control",
    "treatment", "treatment", "treatment", "treatment", "treatment",
    "control", "control", "control", "control", "control"
  ),
  stringsAsFactors = FALSE
)

rownames(meta_all) <- meta_all$id

meta_all$condition <- factor(
  meta_all$condition,
  levels = c("control", "treatment")
)

levels(meta_all$condition) <- c("untrt", "trt")

# -----------------------------------------------------------------------------
# Four pairwise comparisons
# -----------------------------------------------------------------------------

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix   = c("R0", "R2", "R4", "R8"),
  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Repository and file helpers
# -----------------------------------------------------------------------------

resolve_existing_file <- function(candidates, label) {
  hits <- candidates[file.exists(candidates)]

  if (length(hits) == 0L) {
    stop(
      "Could not find ", label, ". Tried: ",
      paste(candidates, collapse = " | ")
    )
  }

  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0L) {
    candidate <- sub("^--file=", "", file_arg[1])
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  NA_character_
}

find_repo_root <- function() {
  candidates <- c(
    getwd(),
    dirname(getwd()),
    "/root/REAPER98632"
  )

  script_path <- get_script_path()

  if (!is.na(script_path)) {
    candidates <- c(dirname(script_path), candidates)
  }

  candidates <- unique(candidates[file.exists(candidates) | dir.exists(candidates)])

  for (cand in candidates) {
    if (dir.exists(file.path(cand, ".git"))) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }

  for (cand in candidates) {
    if (file.exists(file.path(cand, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")) ||
        file.exists(file.path(cand, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"))) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

repo_root <- find_repo_root()

output_dir <- file.path(
  repo_root,
  "exports",
  "manuscript_final_clean"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fig_file <- function(dir, cmp, track_key, ds_key, tag) {
  file.path(
    dir,
    paste0(
      "Figure_",
      cmp,
      "_",
      unname(track_short[track_key]),
      "_",
      unname(dataset_short[ds_key]),
      "_",
      tag,
      ".png"
    )
  )
}

tab_file <- function(dir, cmp, track_key, ds_key, tag) {
  file.path(
    dir,
    paste0(
      "Table_",
      cmp,
      "_",
      unname(track_short[track_key]),
      "_",
      unname(dataset_short[ds_key]),
      "_",
      tag,
      ".csv"
    )
  )
}

# -----------------------------------------------------------------------------
# Generic numeric and string helpers
# -----------------------------------------------------------------------------

assert_required_columns <- function(df, required_cols, object_name = "data frame") {
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns in ",
      object_name,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

safe_log10 <- function(x, pseudocount = 1e-12) {
  log10(pmax(x, pseudocount))
}

safe_neglog10 <- function(x, pseudocount = 1e-12) {
  -log10(pmax(x, pseudocount))
}

clip_probabilities <- function(x, eps = 1e-300) {
  x <- unname(as.numeric(x))

  if (!length(x)) {
    return(numeric(0))
  }

  bad <- !is.finite(x) | is.na(x)
  x[bad] <- NA_real_

  good <- !is.na(x)
  x[good] <- pmin(pmax(x[good], eps), 1 - 1e-12)

  x
}

finite_plot_df <- function(df, x_col, y_col) {
  keep <- is.finite(df[[x_col]]) &
    !is.na(df[[x_col]]) &
    is.finite(df[[y_col]]) &
    !is.na(df[[y_col]])

  df[keep, , drop = FALSE]
}

compact_title <- function(x, width = 54) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

compact_caption <- function(x, width = 118) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

save_csv <- function(df, path) {
  write.csv(df, file = path, row.names = FALSE)
}

save_plot <- function(p, path, width = 12.0, height = 8.5, dpi = figure_dpi, bg = "white") {
  ggplot2::ggsave(
    filename = path,
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = bg,
    limitsize = FALSE
  )
}

save_grob <- function(g, path, width = 14.0, height = 8.5, dpi = figure_dpi, bg = "white") {
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
  parts <- strsplit(dataset_name, "_", fixed = TRUE)[[1]]

  if (length(parts) < 4L) {
    return(dataset_name)
  }

  comparison_name <- paste(parts[1], parts[2], sep = "_")
  dataset_key <- paste(parts[3:length(parts)], collapse = "_")

  paste(comparison_name, pretty_dataset_type(dataset_key), sep = " | ")
}

pretty_group_label <- function(group_label) {
  switch(
    group_label,
    trt = "Treatment",
    untrt = "Control",
    treatment = "Treatment",
    control = "Control",
    group_label
  )
}

clean_gene_set <- function(x) {
  unique(tolower(trimws(x[!is.na(x) & x != ""])))
}

# -----------------------------------------------------------------------------
# Statistical helper functions
# -----------------------------------------------------------------------------

make_design_formula <- function(coldata) {
  ~ condition
}

get_condition_coef <- function(dds) {
  rn <- DESeq2::resultsNames(dds)
  idx <- grep("^condition_", rn)

  if (length(idx) == 0L) {
    stop("Could not identify condition coefficient in resultsNames(dds).")
  }

  rn[idx[1]]
}

classify_effect_strength <- function(res_strong_padj, res_weak_padj, alpha = alpha_level) {
  out <- rep("intermediate", length(res_strong_padj))

  out[!is.na(res_weak_padj) & res_weak_padj < alpha] <- "weak_effect"
  out[!is.na(res_strong_padj) & res_strong_padj < alpha] <- "strong_effect"

  out
}

resolve_top_n_cutoff <- function(sorted_values_desc, top_n = top_n_target) {
  n_total <- length(sorted_values_desc)

  if (n_total == 0L) {
    stop("resolve_top_n_cutoff() received an empty vector.")
  }

  top_n_actual <- min(max(1L, as.integer(top_n)), n_total)
  cutoff_value <- sorted_values_desc[top_n_actual]
  cutoff_quantile <- 1 - (top_n_actual / n_total)

  list(
    top_n_actual = top_n_actual,
    cutoff_value = cutoff_value,
    cutoff_quantile = cutoff_quantile,
    n_total = n_total
  )
}

run_empirical_null_fdrtool <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  stat_vec <- unname(stat_vec)

  if (length(stat_vec) < 5L) {
    stop(
      sprintf(
        "[%s] Fewer than 5 finite Wald statistics were available for fdrtool.",
        dataset_name
      )
    )
  }

  fit <- tryCatch(
    fdrtool::fdrtool(
      stat_vec,
      statistic = "normal",
      plot = FALSE,
      verbose = FALSE,
      cutoff.method = "fndr",
      pct0 = 0.75
    ),
    error = function(e1) {
      message(
        sprintf(
          "[%s] Primary fdrtool call failed: %s",
          dataset_name,
          conditionMessage(e1)
        )
      )

      tryCatch(
        fdrtool::fdrtool(
          as.vector(stat_vec),
          statistic = "normal",
          plot = FALSE,
          verbose = FALSE,
          cutoff.method = "pct0",
          pct0 = 0.75
        ),
        error = function(e2) {
          stop(
            sprintf(
              "[%s] fdrtool failed after retry: %s",
              dataset_name,
              conditionMessage(e2)
            )
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

  if (length(sorted_empirical_p) < 5L) {
    return(NA_real_)
  }

  out <- suppressWarnings(
    tryCatch(
      fdrtool::hc.thresh(as.vector(sorted_empirical_p)),
      error = function(e) {
        message(
          sprintf(
            "[%s] hc.thresh failed: %s",
            dataset_name,
            conditionMessage(e)
          )
        )
        NA_real_
      }
    )
  )

  out <- as.numeric(out[1])

  if (!is.finite(out) || is.na(out) || out <= 0 || out >= 1) {
    return(NA_real_)
  }

  out
}

# -----------------------------------------------------------------------------
# Plotting conventions
# -----------------------------------------------------------------------------

condition_shapes <- c(
  untrt = 21,
  trt = 24
)

condition_fills <- c(
  untrt = plot_palette$control,
  trt = plot_palette$treatment
)

condition_labels <- c(
  untrt = "Control",
  trt = "Treatment"
)

final_class_levels <- c(
  "Background",
  "Weak CNH",
  "Strong CNH",
  "Standard",
  "HBFSS"
)

final_class_colors <- c(
  "Background" = plot_palette$background,
  "Weak CNH" = plot_palette$weak,
  "Strong CNH" = plot_palette$strong,
  "Standard" = plot_palette$standard,
  "HBFSS" = plot_palette$hbfss
)

final_class_shapes <- c(
  "Background" = 16,
  "Weak CNH" = 16,
  "Strong CNH" = 17,
  "Standard" = 15,
  "HBFSS" = 18
)

plot_expand_xy <- function() {
  list(
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.10))),
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))
  )
}

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_theme_size + 1,
        hjust = 0.5,
        lineheight = 1.00,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = base_theme_size - 1,
        hjust = 0.5,
        lineheight = 1.00,
        margin = margin(b = 5)
      ),
      plot.caption = element_text(
        size = base_theme_size - 3,
        hjust = 0.5,
        colour = "grey30",
        lineheight = 0.98,
        margin = margin(t = 6)
      ),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(1, 1, 1, 1),
      legend.spacing.x = unit(4, "pt"),
      legend.spacing.y = unit(1, "pt"),
      legend.text = element_text(size = base_theme_size - 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      plot.margin = margin(10, 12, 10, 10)
    )
}

shared_panel_legend <- function(plot_obj) {
  g <- ggplotGrob(plot_obj + theme(legend.position = "bottom"))
  guide_idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")

  if (length(guide_idx) == 0L) {
    return(NULL)
  }

  g$grobs[[guide_idx[1]]]
}

strip_legend <- function(p) {
  p + theme(legend.position = "none")
}

assemble_one_legend_panel <- function(plot_list, panel_title, ncol = length(plot_list), width_legend = TRUE) {
  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0L) {
    return(NULL)
  }

  legend <- shared_panel_legend(plot_list[[1]])
  no_legend <- lapply(plot_list, strip_legend)

  row <- do.call(
    gridExtra::arrangeGrob,
    c(no_legend, list(ncol = ncol))
  )

  if (is.null(legend)) {
    return(
      gridExtra::arrangeGrob(
        row,
        ncol = 1,
        top = grid::textGrob(
          panel_title,
          gp = grid::gpar(fontface = "bold", cex = 1.15)
        )
      )
    )
  }

  gridExtra::arrangeGrob(
    row,
    legend,
    ncol = 1,
    heights = c(12, 1.4),
    top = grid::textGrob(
      panel_title,
      gp = grid::gpar(fontface = "bold", cex = 1.15)
    )
  )
}

write_export_manifest <- function(root_dir) {
  files <- list.files(
    root_dir,
    recursive = TRUE,
    full.names = TRUE
  )

  files <- files[file.info(files)$isdir %in% FALSE]

  root_norm <- normalizePath(
    root_dir,
    winslash = "/",
    mustWork = FALSE
  )

  file_norm <- normalizePath(
    files,
    winslash = "/",
    mustWork = FALSE
  )

  manifest <- data.frame(
    file = sub(paste0("^", root_norm, "/?"), "", file_norm),
    size_bytes = file.info(files)$size,
    stringsAsFactors = FALSE
  )

  write.csv(
    manifest,
    file.path(root_dir, "Table_Export_Manifest.csv"),
    row.names = FALSE
  )

  manifest
}

git_run <- function(args) {
  tryCatch(
    system2("git", args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      message("Git command failed: ", paste(args, collapse = " "))
      message(conditionMessage(e))
      character(0)
    }
  )
}

finalize_git_exports <- function(root_dir, repo_root) {
  if (!dir.exists(file.path(repo_root, ".git"))) {
    message("No .git directory detected. Exports were written but not staged.")
    return(invisible(FALSE))
  }

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(repo_root)

  rel_export_dir <- sub(
    paste0("^", normalizePath(repo_root, winslash = "/", mustWork = TRUE), "/"),
    "",
    normalizePath(root_dir, winslash = "/", mustWork = TRUE)
  )

  if (isTRUE(GIT_STAGE_EXPORTS)) {
    git_run(c("add", rel_export_dir))

    script_path <- get_script_path()
    if (!is.na(script_path) && file.exists(script_path)) {
      git_run(c("add", script_path))
    }

    message("Git status after staging exports:")
    print(git_run(c("status", "--short")))
  }

  if (isTRUE(GIT_COMMIT_EXPORTS)) {
    git_run(c("commit", "-m", GIT_COMMIT_MESSAGE))
  }

  if (isTRUE(GIT_PUSH_EXPORTS)) {
    branch <- trimws(system2("git", c("branch", "--show-current"), stdout = TRUE))
    if (nzchar(branch)) {
      git_run(c("push", "origin", branch))
    } else {
      message("Could not determine current branch; skipping git push.")
    }
  }

  invisible(TRUE)
}

# =============================================================================
# SECTION 2 OF 5
# IMPORT COUNT MATRIX AND ANNOTATION
# =============================================================================

count_file <- resolve_existing_file(count_file_candidates, "WTTS count file")

message("Using count file: ", count_file)
message("Repository root: ", repo_root)
message("Output directory: ", output_dir)

WTTS_Seq <- read.csv(
  count_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

WTTS_Seq <- as.data.frame(
  WTTS_Seq,
  stringsAsFactors = FALSE
)

WTTS_Seq$OrigID <- as.character(WTTS_Seq$OrigID)
WTTS_Seq$Symbol <- as.character(WTTS_Seq$Symbol)

assert_required_columns(
  WTTS_Seq,
  c("OrigID", "Symbol"),
  object_name = "WTTS count file"
)

assert_required_columns(
  WTTS_Seq,
  meta_all$id,
  object_name = "WTTS count file sample columns"
)

WTTS_Seq <- WTTS_Seq[
  !is.na(WTTS_Seq$OrigID) & !is.na(WTTS_Seq$Symbol),
  ,
  drop = FALSE
]

sample_na <- rowSums(is.na(WTTS_Seq[, meta_all$id, drop = FALSE])) > 0

WTTS_Seq <- WTTS_Seq[!sample_na, , drop = FALSE]

rownames(WTTS_Seq) <- make.unique(WTTS_Seq$OrigID)

OrigID_Symbol <- unique(WTTS_Seq[, c("OrigID", "Symbol"), drop = FALSE])

colnames(OrigID_Symbol) <- c("feature_id", "gene_symbol")

OrigID_Symbol$feature_id <- as.character(OrigID_Symbol$feature_id)
OrigID_Symbol$gene_symbol <- as.character(OrigID_Symbol$gene_symbol)

OrigID_Symbol <- OrigID_Symbol %>%
  dplyr::mutate(
    gene_symbol = dplyr::if_else(
      is.na(gene_symbol),
      "",
      trimws(gene_symbol)
    )
  ) %>%
  dplyr::arrange(
    feature_id,
    dplyr::desc(gene_symbol != ""),
    gene_symbol
  ) %>%
  dplyr::distinct(
    feature_id,
    .keep_all = TRUE
  ) %>%
  dplyr::mutate(
    gene_symbol = dplyr::na_if(gene_symbol, "")
  )

if (run_twas_overlap) {
  twas_file <- resolve_existing_file(twas_file_candidates, "TWAS file")

  TWAS_Seq <- read.csv(
    twas_file,
    header = TRUE,
    stringsAsFactors = FALSE
  )

  TWAS_Seq <- as.data.frame(TWAS_Seq)

  if (ncol(TWAS_Seq) < 4L) {
    stop("TWAS file must contain at least 4 columns.")
  }

  TWAS_data <- TWAS_Seq[, c(1, 4), drop = FALSE]
  colnames(TWAS_data) <- c("source_id", "gene_symbol")
}

# =============================================================================
# SECTION 3 OF 5
# EVS CONSTRUCTION AND EVS SUPPORT FIGURES
# =============================================================================

prepare_comparison_data <- function(comparison_name, group1_prefix, group2_prefix, WTTS_Seq, meta_all) {
  keep_ids <- grepl(paste0("^", group1_prefix, "_"), meta_all$id) |
    grepl(paste0("^", group2_prefix, "_"), meta_all$id)

  meta_sub <- meta_all[keep_ids, , drop = FALSE]

  coldata <- meta_sub[, c("condition"), drop = FALSE]

  sample_ids <- rownames(meta_sub)

  missing_samples <- setdiff(sample_ids, colnames(WTTS_Seq))

  if (length(missing_samples) > 0L) {
    stop(
      "Missing samples in WTTS file for ",
      comparison_name,
      ": ",
      paste(missing_samples, collapse = ", ")
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

make_rank_matrix_for_track <- function(count_matrix, coldata, track_key) {
  track_key <- match.arg(track_key, c("normalized_evs", "raw_evs"))

  if (track_key == "raw_evs") {
    return(
      list(
        rank_matrix = as.data.frame(count_matrix),
        preprocessing_label = "Raw counts prior to EVS"
      )
    )
  }

  dds_init <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData = coldata,
    design = make_design_formula(coldata)
  )

  dds_init <- dds_init[rowSums(DESeq2::counts(dds_init)) > 0, ]

  dds_init <- DESeq2::estimateSizeFactors(dds_init)

  norm_counts <- as.data.frame(
    DESeq2::counts(dds_init, normalized = TRUE)
  )

  list(
    rank_matrix = norm_counts,
    preprocessing_label = "DESeq2-normalized counts prior to EVS"
  )
}

compute_pc1_loading_table <- function(value_df, sample_names, top_n = top_n_target, preprocessing_label = "Normalized prior to EVS") {
  x <- as.matrix(value_df[, sample_names, drop = FALSE])

  pca_fit <- stats::prcomp(
    t(x),
    center = TRUE,
    scale. = FALSE,
    rank. = 2
  )

  loading_abs <- abs(pca_fit$rotation[, 1])

  loading_tbl <- data.frame(
    feature_id = names(loading_abs),
    pc1_loading_abs = unname(loading_abs),
    stringsAsFactors = FALSE
  )

  loading_tbl <- loading_tbl[
    order(loading_tbl$pc1_loading_abs, decreasing = TRUE),
    ,
    drop = FALSE
  ]

  loading_tbl$rank <- seq_len(nrow(loading_tbl))

  cutoff_info <- resolve_top_n_cutoff(
    loading_tbl$pc1_loading_abs,
    top_n = top_n
  )

  cutoff <- cutoff_info$cutoff_value

  loading_tbl$split_class <- ifelse(
    loading_tbl$pc1_loading_abs >= cutoff,
    "high_loading",
    "background_loading"
  )

  list(
    pca_fit = pca_fit,
    loading_table = loading_tbl,
    cutoff = cutoff,
    top_n_used = cutoff_info$top_n_actual,
    cutoff_quantile = cutoff_info$cutoff_quantile,
    preprocessing_label = preprocessing_label
  )
}

build_eigenvector_split <- function(count_matrix, coldata, track_key) {
  track_key <- match.arg(track_key, c("normalized_evs", "raw_evs"))

  rank_obj <- make_rank_matrix_for_track(
    count_matrix = count_matrix,
    coldata = coldata,
    track_key = track_key
  )

  rank_matrix <- rank_obj$rank_matrix
  preprocessing_label <- rank_obj$preprocessing_label

  sample_ids <- colnames(count_matrix)

  trt_ids <- sample_ids[coldata$condition == "trt"]
  untrt_ids <- sample_ids[coldata$condition == "untrt"]

  fit_trt <- compute_pc1_loading_table(
    rank_matrix,
    trt_ids,
    top_n = top_n_target,
    preprocessing_label = preprocessing_label
  )

  fit_untrt <- compute_pc1_loading_table(
    rank_matrix,
    untrt_ids,
    top_n = top_n_target,
    preprocessing_label = preprocessing_label
  )

  trt_high <- as.character(
    subset(
      fit_trt$loading_table,
      split_class == "high_loading"
    )$feature_id
  )

  untrt_high <- as.character(
    subset(
      fit_untrt$loading_table,
      split_class == "high_loading"
    )$feature_id
  )

  leading_edge_ids <- union(trt_high, untrt_high)

  remainder_ids <- setdiff(
    rownames(count_matrix),
    leading_edge_ids
  )

  if (length(leading_edge_ids) == 0L) {
    stop("Leading-edge dataset is empty. Check sample mapping or top_n_target.")
  }

  if (length(remainder_ids) == 0L) {
    stop("Remainder dataset is empty. top_n_target may be too large.")
  }

  list(
    track_key = track_key,
    preprocessing_label = preprocessing_label,
    rank_matrix = rank_matrix,
    fit_trt = fit_trt,
    fit_untrt = fit_untrt,
    leading_edge_ids = leading_edge_ids,
    remainder_ids = remainder_ids,
    raw_dataset = count_matrix,
    leading_edge_dataset = count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset = count_matrix[remainder_ids, , drop = FALSE]
  )
}

plot_pca_scatter <- function(pca_fit, dataset_label, group_label, preprocessing_label) {
  pca_var <- pca_fit$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)

  cond_key <- ifelse(
    group_label %in% c("control", "untrt"),
    "untrt",
    "trt"
  )

  pca_df <- data.frame(
    Sample = rownames(pca_fit$x),
    PC1 = pca_fit$x[, 1],
    PC2 = pca_fit$x[, 2],
    Condition = cond_key,
    stringsAsFactors = FALSE
  )

  ggplot(
    pca_df,
    aes(
      PC1,
      PC2,
      label = Sample,
      shape = Condition,
      fill = Condition
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      linetype = "dashed",
      colour = "grey70"
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.3,
      linetype = "dashed",
      colour = "grey70"
    ) +
    geom_point(
      size = 2.8,
      colour = "white",
      stroke = 0.55
    ) +
    ggrepel::geom_text_repel(
      size = 1.9,
      max.overlaps = 8,
      force = 1.0,
      box.padding = 0.22,
      point.padding = 0.10,
      min.segment.length = 0
    ) +
    scale_shape_manual(
      values = condition_shapes,
      labels = condition_labels,
      name = "Condition"
    ) +
    scale_fill_manual(
      values = condition_fills,
      labels = condition_labels,
      name = "Condition"
    ) +
    labs(
      title = compact_title(
        paste(dataset_label, "|", pretty_group_label(group_label), "PCA"),
        width = 44
      ),
      subtitle = paste0(
        preprocessing_label,
        " · PC1 = ",
        pca_var_per[1],
        "%; PC2 = ",
        pca_var_per[2],
        "%."
      ),
      x = paste0("PC1 (", pca_var_per[1], "%)"),
      y = paste0("PC2 (", pca_var_per[2], "%)")
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()
}

plot_pc1_loading_rank <- function(loading_tbl, cutoff, dataset_label, group_label, top_n_used, cutoff_quantile, preprocessing_label) {
  quantile_label <- if (is.finite(cutoff_quantile)) {
    paste0("Upper-tail quantile = ", signif(cutoff_quantile, 4))
  } else {
    "Upper-tail quantile = NA"
  }

  ggplot(
    loading_tbl,
    aes(rank, pc1_loading_abs)
  ) +
    geom_line(
      linewidth = 0.4,
      color = "grey35"
    ) +
    geom_hline(
      yintercept = cutoff,
      color = plot_palette$threshold,
      linewidth = 0.9
    ) +
    annotate(
      "text",
      x = max(loading_tbl$rank) * 0.72,
      y = cutoff,
      label = paste0(
        "Top ",
        top_n_used,
        " cutoff = ",
        signif(cutoff, 4),
        "\n",
        quantile_label
      ),
      color = plot_palette$threshold,
      vjust = -0.8,
      size = 3.7
    ) +
    labs(
      title = compact_title(
        paste(dataset_label, "|", pretty_group_label(group_label), "PC1 loading rank"),
        width = 52
      ),
      subtitle = paste(
        "Absolute PC1 loading ranked within condition.",
        preprocessing_label
      ),
      x = "Ranked PAS feature",
      y = "Absolute PC1 loading"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()
}

plot_eigenvector_histograms <- function(loading_tbl, cutoff, dataset_label, group_label) {
  full_vals <- loading_tbl$pc1_loading_abs
  lead_vals <- loading_tbl$pc1_loading_abs[loading_tbl$pc1_loading_abs >= cutoff]
  rem_vals <- loading_tbl$pc1_loading_abs[loading_tbl$pc1_loading_abs < cutoff]

  cutoff_tx <- -log(pmax(cutoff * 100, 1e-6))

  p1 <- ggplot(
    data.frame(x = -log(pmax(full_vals * 100, 1e-6))),
    aes(x)
  ) +
    geom_histogram(
      bins = 50,
      fill = plot_palette$histogram,
      color = "white"
    ) +
    geom_vline(
      xintercept = cutoff_tx,
      color = plot_palette$threshold,
      linewidth = 0.9
    ) +
    labs(
      title = "All sites",
      x = "-log(|PC1 loading| x 100)",
      y = "Count"
    ) +
    manuscript_theme()

  p2 <- ggplot(
    data.frame(x = -log(pmax(rem_vals * 100, 1e-6))),
    aes(x)
  ) +
    geom_histogram(
      bins = 50,
      fill = plot_palette$histogram,
      color = "white"
    ) +
    geom_vline(
      xintercept = cutoff_tx,
      color = plot_palette$threshold,
      linewidth = 0.9
    ) +
    labs(
      title = "Remainder",
      x = "-log(|PC1 loading| x 100)",
      y = "Count"
    ) +
    manuscript_theme()

  p3 <- ggplot(
    data.frame(x = -log(pmax(lead_vals * 100, 1e-6))),
    aes(x)
  ) +
    geom_histogram(
      bins = 50,
      fill = plot_palette$histogram,
      color = "white"
    ) +
    geom_vline(
      xintercept = cutoff_tx,
      color = plot_palette$threshold,
      linewidth = 0.9
    ) +
    labs(
      title = "Leading edge",
      x = "-log(|PC1 loading| x 100)",
      y = "Count"
    ) +
    manuscript_theme()

  gridExtra::arrangeGrob(
    p1,
    p2,
    p3,
    ncol = 3,
    top = grid::textGrob(
      paste(
        dataset_label,
        "|",
        pretty_group_label(group_label),
        "EVS loading distributions"
      ),
      gp = grid::gpar(fontface = "bold", cex = 1.15)
    )
  )
}

compute_mean_expression_table <- function(raw_counts, coldata) {
  sample_ids <- colnames(raw_counts)
  trt_ids <- sample_ids[coldata$condition == "trt"]
  untrt_ids <- sample_ids[coldata$condition == "untrt"]

  data.frame(
    feature_id = as.character(rownames(raw_counts)),
    mean_trt = rowMeans(raw_counts[, trt_ids, drop = FALSE]),
    mean_untrt = rowMeans(raw_counts[, untrt_ids, drop = FALSE]),
    mean_all = rowMeans(raw_counts),
    stringsAsFactors = FALSE
  )
}

plot_mean_histogram_panel <- function(df_means, base_mean_vec, dataset_name) {
  p1 <- ggplot(
    data.frame(x = safe_log10(df_means$mean_trt + 1)),
    aes(x)
  ) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    labs(title = "Treatment mean", x = "log10(mean + 1)", y = "Count") +
    manuscript_theme()

  p2 <- ggplot(
    data.frame(x = safe_log10(df_means$mean_untrt + 1)),
    aes(x)
  ) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    labs(title = "Control mean", x = "log10(mean + 1)", y = "Count") +
    manuscript_theme()

  p3 <- ggplot(
    data.frame(x = safe_log10(df_means$mean_all + 1)),
    aes(x)
  ) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    labs(title = "Pooled mean", x = "log10(mean + 1)", y = "Count") +
    manuscript_theme()

  p4 <- ggplot(
    data.frame(x = safe_log10(base_mean_vec + 1)),
    aes(x)
  ) +
    geom_histogram(bins = 50, fill = plot_palette$histogram, color = "white") +
    labs(title = "DESeq2 baseMean", x = "log10(baseMean + 1)", y = "Count") +
    manuscript_theme()

  gridExtra::arrangeGrob(
    p1,
    p2,
    p3,
    p4,
    ncol = 2,
    top = grid::textGrob(
      paste(dataset_name, "mean-expression histograms"),
      gp = grid::gpar(fontface = "bold", cex = 1.15)
    )
  )
}

# =============================================================================
# SECTION 4 OF 5
# DESEQ2, EMPIRICAL-NULL CALIBRATION, HBFSS, AND FIGURES
# =============================================================================

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
  design_formula <- make_design_formula(coldata)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = coldata,
    design = design_formula
  )

  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]

  dds <- DESeq2::DESeq(
    dds,
    betaPrior = FALSE
  )

  res <- DESeq2::results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    alpha = alpha_level
  )

  res_strong <- DESeq2::results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "greaterAbs"
  )

  res_weak <- DESeq2::results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "lessAbs"
  )

  res_all_df <- as.data.frame(res)
  res_all_df$feature_id <- as.character(rownames(res_all_df))

  valid_stat <- is.finite(res_all_df$stat) & !is.na(res_all_df$stat)

  stat_vec <- as.numeric(res_all_df$stat[valid_stat])
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]

  if (length(stat_vec) < 5L) {
    stop(
      sprintf(
        "[%s] Fewer than 5 finite Wald statistics were available for fdrtool.",
        dataset_name
      )
    )
  }

  fdr_fit <- run_empirical_null_fdrtool(
    stat_vec,
    dataset_name = dataset_name
  )

  res_df <- res_all_df

  n_valid <- sum(valid_stat)

  if (length(fdr_fit$pval) != n_valid ||
      length(fdr_fit$qval) != n_valid ||
      length(fdr_fit$lfdr) != n_valid) {
    stop(sprintf("[%s] fdrtool output length mismatch.", dataset_name))
  }

  res_df$empirical_p <- NA_real_
  res_df$empirical_q <- NA_real_
  res_df$lfdr <- NA_real_

  res_df$empirical_p[valid_stat] <- as.numeric(fdr_fit$pval)
  res_df$empirical_q[valid_stat] <- as.numeric(fdr_fit$qval)
  res_df$lfdr[valid_stat] <- as.numeric(fdr_fit$lfdr)

  res_df$empirical_bh <- NA_real_

  valid_empirical <- is.finite(res_df$empirical_p) & !is.na(res_df$empirical_p)

  if (any(valid_empirical)) {
    res_df$empirical_bh[valid_empirical] <- p.adjust(
      res_df$empirical_p[valid_empirical],
      method = "BH"
    )
  }

  res_df$pval <- res_df$empirical_p
  res_df$padjc <- res_df$empirical_bh
  res_df$qval <- res_df$empirical_q

  coef_name <- get_condition_coef(dds)

  shr <- DESeq2::lfcShrink(
    dds,
    coef = coef_name,
    type = "apeglm",
    res = res
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

  hc_p_threshold_dataset <- safe_hc_thresh(
    res_df$empirical_p,
    dataset_name = dataset_name
  )

  res_df$gene_empirical_pvalue <- res_df$empirical_p

  res_df$HBFSS <- abs(
    res_df$lfc_shrunk *
      log10(pmax(res_df$gene_empirical_pvalue, 1e-300))
  )

  if (is.na(hc_p_threshold_dataset) ||
      !is.finite(hc_p_threshold_dataset) ||
      hc_p_threshold_dataset <= 0 ||
      hc_p_threshold_dataset >= 1) {
    hbfss_threshold_dataset <- NA_real_
    res_df$HBFSS_core_pass <- FALSE
  } else {
    hbfss_threshold_dataset <- abs(log10(hc_p_threshold_dataset)) * lfc_boundary
    res_df$HBFSS_core_pass <- res_df$HBFSS >= hbfss_threshold_dataset
  }

  res_df$regulation_direction <- ifelse(
    is.na(res_df$lfc_shrunk),
    NA_character_,
    ifelse(
      res_df$lfc_shrunk > 0,
      "upregulated",
      ifelse(res_df$lfc_shrunk < 0, "downregulated", "no_change")
    )
  )

  res_df$raw_lfc_pass <- !is.na(res_df$log2FoldChange) &
    abs(res_df$log2FoldChange) >= lfc_boundary

  res_df$shrunk_lfc_pass <- !is.na(res_df$lfc_shrunk) &
    abs(res_df$lfc_shrunk) >= lfc_boundary

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

  res_df$resGA_padj <- res_df$padj_strong_effect
  res_df$resLA_padj <- res_df$padj_weak_effect

  res_df$effect_class <- classify_effect_strength(
    res_df$padj_strong_effect,
    res_df$padj_weak_effect,
    alpha = alpha_level
  )

  res_df$hc_pass <- !is.na(hc_p_threshold_dataset) &
    is.finite(hc_p_threshold_dataset) &
    !is.na(res_df$empirical_p) &
    is.finite(res_df$empirical_p) &
    res_df$empirical_p <= hc_p_threshold_dataset

  res_df$weak_cnh_flag <- !is.na(res_df$resLA_padj) &
    res_df$resLA_padj < alpha_level &
    !is.na(res_df$lfc_shrunk) &
    abs(res_df$lfc_shrunk) < lfc_boundary &
    res_df$hc_pass

  res_df$strong_cnh_flag <- !is.na(res_df$resGA_padj) &
    res_df$resGA_padj < alpha_level &
    !is.na(res_df$lfc_shrunk) &
    abs(res_df$lfc_shrunk) >= lfc_boundary &
    res_df$hc_pass

  res_df$standard_significant <- !is.na(res_df$padj) &
    res_df$padj < alpha_level &
    res_df$raw_lfc_pass &
    res_df$shrunk_lfc_pass

  res_df$standard_flag <- res_df$standard_significant &
    res_df$hc_pass

  res_df$hbfss_flag <- !is.na(res_df$HBFSS_core_pass) &
    res_df$HBFSS_core_pass &
    res_df$hc_pass

  res_df$HBFSS_significant <- res_df$hbfss_flag

  res_df$weak_hbfss_overlap <- res_df$weak_cnh_flag & res_df$hbfss_flag
  res_df$strong_hbfss_overlap <- res_df$strong_cnh_flag & res_df$hbfss_flag
  res_df$standard_hbfss_overlap <- res_df$standard_flag & res_df$hbfss_flag

  res_df$any_overlap <- res_df$weak_hbfss_overlap |
    res_df$strong_hbfss_overlap |
    res_df$standard_hbfss_overlap

  res_df$final_class <- "Background"
  res_df$final_class[res_df$weak_cnh_flag] <- "Weak CNH"
  res_df$final_class[res_df$strong_cnh_flag] <- "Strong CNH"
  res_df$final_class[res_df$standard_flag] <- "Standard"
  res_df$final_class[res_df$hbfss_flag] <- "HBFSS"

  res_df$final_class <- factor(
    res_df$final_class,
    levels = final_class_levels
  )

  base_mean_vec <- res_df$baseMean[!is.na(res_df$baseMean)]

  norm_counts <- as.data.frame(
    DESeq2::counts(dds, normalized = TRUE)
  )

  norm_counts$feature_id <- as.character(rownames(norm_counts))

  mm <- as.data.frame(S4Vectors::mcols(dds))
  mm$feature_id <- as.character(rownames(mm))

  disp_cols_available <- intersect(
    c(
      "feature_id",
      "dispGeneEst",
      "dispFit",
      "dispersion",
      "dispIter",
      "baseMean",
      "dispOutlier"
    ),
    colnames(mm)
  )

  disp_df <- mm[, disp_cols_available, drop = FALSE]

  annot_df$feature_id <- as.character(annot_df$feature_id)
  annot_df$gene_symbol <- as.character(annot_df$gene_symbol)

  annot_df <- annot_df %>%
    dplyr::mutate(
      gene_symbol = dplyr::if_else(
        is.na(gene_symbol),
        "",
        trimws(gene_symbol)
      )
    ) %>%
    dplyr::arrange(
      feature_id,
      dplyr::desc(gene_symbol != ""),
      gene_symbol
    ) %>%
    dplyr::distinct(
      feature_id,
      .keep_all = TRUE
    ) %>%
    dplyr::mutate(
      gene_symbol = dplyr::na_if(gene_symbol, "")
    )

  norm_counts <- norm_counts[!duplicated(norm_counts$feature_id), , drop = FALSE]
  disp_df <- disp_df[!duplicated(disp_df$feature_id), , drop = FALSE]

  final_df <- res_df %>%
    dplyr::left_join(annot_df, by = "feature_id") %>%
    dplyr::left_join(norm_counts, by = "feature_id") %>%
    dplyr::left_join(disp_df, by = "feature_id")

  final_df$neglog10_padj <- safe_neglog10(final_df$padj)
  final_df$neglog10_empirical_p <- safe_neglog10(final_df$empirical_p)
  final_df$neglog10_padjc <- safe_neglog10(final_df$empirical_bh)

  final_df$dataset_name <- dataset_name
  final_df$hc_p_threshold_dataset <- hc_p_threshold_dataset
  final_df$hbfss_threshold_dataset <- hbfss_threshold_dataset
  final_df$Apeglm_L2FC <- final_df$lfc_shrunk
  final_df$pi_valueE <- final_df$HBFSS

  preferred_cols <- c(
    "dataset_name",
    "feature_id",
    "gene_symbol",
    "baseMean",
    "log2FoldChange",
    "lfc_shrunk",
    "Apeglm_L2FC",
    "regulation_direction",
    "stat",
    "pvalue",
    "padj",
    "empirical_p",
    "empirical_q",
    "empirical_bh",
    "lfdr",
    "HBFSS",
    "pi_valueE",
    "hc_p_threshold_dataset",
    "hbfss_threshold_dataset",
    "resLA_padj",
    "resGA_padj",
    "weak_cnh_flag",
    "strong_cnh_flag",
    "standard_flag",
    "hbfss_flag",
    "final_class",
    "weak_hbfss_overlap",
    "strong_hbfss_overlap",
    "standard_hbfss_overlap",
    "any_overlap"
  )

  final_df <- final_df[
    ,
    c(
      intersect(preferred_cols, names(final_df)),
      setdiff(names(final_df), preferred_cols)
    ),
    drop = FALSE
  ]

  list(
    dds = dds,
    results = final_df,
    base_mean_vec = base_mean_vec,
    hc_p_threshold = hc_p_threshold_dataset,
    hbfss_threshold = hbfss_threshold_dataset
  )
}

build_final_volcano_df <- function(df, y_col = "neglog10_empirical_p") {
  df <- finite_plot_df(
    df,
    "lfc_shrunk",
    y_col
  )

  df$has_valid_gene_symbol <- !is.na(df$gene_symbol) &
    grepl("[A-Za-z0-9]", trimws(df$gene_symbol))

  df$gene_symbol_plot <- ifelse(
    df$has_valid_gene_symbol,
    trimws(df$gene_symbol),
    NA_character_
  )

  df
}

select_final_volcano_labels <- function(df, y_col, n_labels = n_top_labels_volcano) {
  if (!nrow(df)) {
    return(df[0, , drop = FALSE])
  }

  df <- df[df$has_valid_gene_symbol, , drop = FALSE]

  if (!nrow(df)) {
    return(df[0, , drop = FALSE])
  }

  class_priority <- c(
    "HBFSS" = 1,
    "Standard" = 2,
    "Strong CNH" = 3,
    "Weak CNH" = 4,
    "Background" = 5
  )

  df$label_priority <- class_priority[as.character(df$final_class)]

  metric <- suppressWarnings(as.numeric(df[[y_col]]))
  metric[!is.finite(metric)] <- -Inf

  ord <- order(
    df$label_priority,
    -metric,
    -abs(df$lfc_shrunk),
    na.last = TRUE
  )

  df <- df[ord, , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol_plot), , drop = FALSE]

  df[seq_len(min(n_labels, nrow(df))), , drop = FALSE]
}

make_hbfss_boundary_df <- function(plot_df, hbfss_threshold, y_limit) {
  if (!is.finite(hbfss_threshold) || is.na(hbfss_threshold) || hbfss_threshold <= 0) {
    return(NULL)
  }

  x_max <- max(
    max(abs(plot_df$lfc_shrunk), na.rm = TRUE),
    lfc_boundary * 1.1
  )

  x_min <- max(0.05, hbfss_threshold / max(y_limit, 1e-6))

  x_abs <- seq(
    x_min,
    x_max,
    length.out = 600
  )

  y_curve <- hbfss_threshold / x_abs

  keep <- is.finite(y_curve) &
    y_curve >= 0 &
    y_curve <= y_limit

  if (!any(keep)) {
    return(NULL)
  }

  x_abs <- x_abs[keep]
  y_curve <- y_curve[keep]

  rbind(
    data.frame(x = -rev(x_abs), y = rev(y_curve)),
    data.frame(x = x_abs, y = y_curve)
  )
}

plot_final_volcano <- function(df, dataset_name, short_title = NULL) {
  plot_df <- build_final_volcano_df(
    df,
    y_col = "neglog10_empirical_p"
  )

  if (!nrow(plot_df)) {
    stop("No finite volcano plotting rows for ", dataset_name)
  }

  lab_df <- select_final_volcano_labels(
    plot_df,
    y_col = "neglog10_empirical_p"
  )

  hc_raw <- suppressWarnings(
    as.numeric(df$hc_p_threshold_dataset[1])
  )

  hbfss_raw <- suppressWarnings(
    as.numeric(df$hbfss_threshold_dataset[1])
  )

  hc_y <- if (is.finite(hc_raw) && !is.na(hc_raw) && hc_raw > 0 && hc_raw < 1) {
    safe_neglog10(hc_raw)
  } else {
    NA_real_
  }

  y_limit <- max(
    plot_df$neglog10_empirical_p,
    na.rm = TRUE
  ) * 1.05

  boundary_df <- make_hbfss_boundary_df(
    plot_df = plot_df,
    hbfss_threshold = hbfss_raw,
    y_limit = y_limit
  )

  weak_n <- sum(df$weak_cnh_flag, na.rm = TRUE)
  strong_n <- sum(df$strong_cnh_flag, na.rm = TRUE)
  std_n <- sum(df$standard_flag, na.rm = TRUE)
  hbfss_n <- sum(df$hbfss_flag, na.rm = TRUE)
  overlap_n <- sum(df$any_overlap, na.rm = TRUE)

  count_text <- paste0(
    "Weak=", weak_n,
    "  Strong=", strong_n,
    "  Std=", std_n,
    "  HBFSS=", hbfss_n,
    "  Overlap=", overlap_n
  )

  plot_title <- if (is.null(short_title)) {
    compact_title(pretty_dataset_label(dataset_name), width = 42)
  } else {
    short_title
  }

  p <- ggplot(
    plot_df,
    aes(
      lfc_shrunk,
      neglog10_empirical_p,
      color = final_class,
      shape = final_class
    )
  ) +
    geom_point(
      alpha = 0.80,
      size = 1.35,
      stroke = 0.30
    ) +
    scale_color_manual(
      values = final_class_colors,
      drop = FALSE,
      name = "Class"
    ) +
    scale_shape_manual(
      values = final_class_shapes,
      drop = FALSE,
      name = "Class"
    ) +
    geom_vline(
      xintercept = c(-lfc_boundary, lfc_boundary),
      linetype = "dashed",
      linewidth = 0.55,
      colour = plot_palette$threshold
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      linewidth = 0.35,
      colour = "grey55"
    ) +
    labs(
      title = plot_title,
      subtitle = "x = shrunken log2FC; y = -log10(empirical p)",
      x = "Shrunken log2FC",
      y = expression(-log[10]("Empirical p")),
      caption = count_text
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    plot_expand_xy() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin = margin(10, 12, 10, 10)
    )

  if (!is.na(hc_y) && is.finite(hc_y)) {
    x_rng <- range(plot_df$lfc_shrunk, na.rm = TRUE)

    p <- p +
      geom_hline(
        yintercept = hc_y,
        linetype = "dotted",
        linewidth = 0.65,
        color = plot_palette$hc
      ) +
      annotate(
        "text",
        x = x_rng[1] + diff(x_rng) * 0.04,
        y = hc_y,
        label = paste0("HC=", signif(hc_raw, 3)),
        color = plot_palette$hc,
        hjust = 0,
        vjust = -0.35,
        size = 2.9
      )
  }

  if (!is.null(boundary_df)) {
    p <- p +
      geom_line(
        data = boundary_df,
        aes(x, y),
        inherit.aes = FALSE,
        color = plot_palette$hbfss_line,
        linewidth = 0.80
      )

    x_rng <- range(plot_df$lfc_shrunk, na.rm = TRUE)
    y_rng <- range(plot_df$neglog10_empirical_p, na.rm = TRUE)

    p <- p +
      annotate(
        "text",
        x = x_rng[2] - diff(x_rng) * 0.24,
        y = y_rng[1] + diff(y_rng) * 0.14,
        label = paste0("HBFSS=", signif(hbfss_raw, 3)),
        color = plot_palette$hbfss_line,
        hjust = 0,
        vjust = 0,
        size = 2.9
      )
  }

  x_rng <- range(plot_df$lfc_shrunk, na.rm = TRUE)
  y_rng <- range(plot_df$neglog10_empirical_p, na.rm = TRUE)

  p <- p +
    annotate(
      "text",
      x = 0,
      y = y_rng[2] - diff(y_rng) * 0.08,
      label = count_text,
      size = 3.0,
      fontface = "plain"
    ) +
    annotate(
      "text",
      x = -lfc_boundary,
      y = y_rng[1] + diff(y_rng) * 0.06,
      label = "LFC=-1",
      color = plot_palette$threshold,
      hjust = 1.05,
      size = 2.8
    ) +
    annotate(
      "text",
      x = lfc_boundary,
      y = y_rng[1] + diff(y_rng) * 0.06,
      label = "LFC=1",
      color = plot_palette$threshold,
      hjust = -0.05,
      size = 2.8
    )

  if (nrow(lab_df) > 0L) {
    p <- p +
      ggrepel::geom_text_repel(
        data = lab_df,
        aes(label = gene_symbol_plot),
        size = 1.8,
        seed = 1,
        max.overlaps = 20,
        force = 1.15,
        force_pull = 0.4,
        box.padding = 0.30,
        point.padding = 0.14,
        min.segment.length = 0,
        segment.alpha = 0.55,
        segment.size = 0.20
      )
  }

  p
}

plot_empirical_p_histogram <- function(df, dataset_name, hc_p_threshold) {
  subtitle_text <- if (is.na(hc_p_threshold) || hc_p_threshold <= 0 || hc_p_threshold >= 1) {
    "Empirical-null p-values; HC threshold unavailable"
  } else {
    paste0("Empirical-null p-values; HC threshold = ", signif(hc_p_threshold, 4))
  }

  p <- ggplot(
    df,
    aes(empirical_p)
  ) +
    geom_histogram(
      bins = 60,
      fill = plot_palette$histogram,
      color = "white"
    ) +
    labs(
      title = compact_title(
        paste(pretty_dataset_label(dataset_name), "| empirical-null p-value distribution")
      ),
      subtitle = subtitle_text,
      x = "Empirical-null p-value",
      y = "Count"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()

  if (!is.na(hc_p_threshold) &&
      is.finite(hc_p_threshold) &&
      hc_p_threshold > 0 &&
      hc_p_threshold < 1) {
    p <- p +
      geom_vline(
        xintercept = hc_p_threshold,
        color = plot_palette$hc,
        linewidth = 0.9
      )
  }

  p
}

plot_hc_profile <- function(df, dataset_name, hc_p_threshold) {
  pvals <- sort(
    df$empirical_p[
      is.finite(df$empirical_p) &
        !is.na(df$empirical_p) &
        df$empirical_p > 0 &
        df$empirical_p < 1
    ]
  )

  if (length(pvals) < 5L) {
    return(NULL)
  }

  n <- length(pvals)
  i <- seq_len(n)

  v <- (i / n) * (1 - (i / n)) / n
  v[v == 0] <- min(v[v > 0])

  hc_score <- abs((i / n) - pvals) / sqrt(v)

  hc_df <- data.frame(
    rank = i,
    empirical_p = pvals,
    hc_score = hc_score
  )

  peak_df <- hc_df[which.max(hc_df$hc_score), , drop = FALSE]

  p <- ggplot(
    hc_df,
    aes(empirical_p, hc_score)
  ) +
    geom_line(
      linewidth = 0.5,
      colour = "grey30"
    ) +
    geom_point(
      data = peak_df,
      aes(empirical_p, hc_score),
      inherit.aes = FALSE,
      size = 2.2,
      colour = plot_palette$hc
    ) +
    labs(
      title = compact_title(
        paste(pretty_dataset_label(dataset_name), "| higher-criticism profile")
      ),
      subtitle = paste0(
        "HC threshold = ",
        ifelse(
          is.finite(hc_p_threshold) && !is.na(hc_p_threshold),
          signif(hc_p_threshold, 4),
          "NA"
        )
      ),
      x = "Sorted empirical p-value",
      y = "HC score"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()

  if (!is.na(hc_p_threshold) &&
      is.finite(hc_p_threshold) &&
      hc_p_threshold > 0 &&
      hc_p_threshold < 1) {
    p <- p +
      geom_vline(
        xintercept = hc_p_threshold,
        color = plot_palette$hc,
        linewidth = 0.9
      )
  }

  p
}

plot_hbfss_distribution <- function(df, dataset_name, hbfss_threshold) {
  p <- ggplot(
    df,
    aes(HBFSS)
  ) +
    geom_histogram(
      bins = 70,
      fill = plot_palette$histogram,
      color = "white"
    ) +
    labs(
      title = compact_title(
        paste(pretty_dataset_label(dataset_name), "| HBFSS distribution")
      ),
      subtitle = if (is.finite(hbfss_threshold) && !is.na(hbfss_threshold)) {
        paste0("Threshold = ", signif(hbfss_threshold, 4))
      } else {
        "Threshold unavailable"
      },
      x = "HBFSS",
      y = "Count"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()

  if (!is.na(hbfss_threshold) && is.finite(hbfss_threshold)) {
    p <- p +
      geom_vline(
        xintercept = hbfss_threshold,
        color = plot_palette$hbfss_line,
        linewidth = 0.9
      )
  }

  p
}

plot_shrunken_ma <- function(df, dataset_name) {
  if (!("baseMean" %in% colnames(df))) {
    return(NULL)
  }

  ma_df <- df[
    is.finite(df$baseMean) &
      !is.na(df$baseMean) &
      is.finite(df$lfc_shrunk) &
      !is.na(df$lfc_shrunk),
    ,
    drop = FALSE
  ]

  ggplot(
    ma_df,
    aes(safe_log10(baseMean + 1), lfc_shrunk)
  ) +
    geom_point(
      alpha = 0.55,
      size = 1.0,
      color = "grey40"
    ) +
    geom_hline(
      yintercept = c(-lfc_boundary, lfc_boundary),
      linetype = "dashed",
      color = plot_palette$threshold,
      linewidth = 0.7
    ) +
    labs(
      title = compact_title(
        paste(pretty_dataset_label(dataset_name), "| MA plot")
      ),
      subtitle = "Shrunken effect size versus baseMean",
      x = expression(log[10]("baseMean + 1")),
      y = "Shrunken log2 fold change"
    ) +
    manuscript_theme()
}

plot_dispersion_cloud <- function(dds, dataset_name, fig_subdir, cmp_short, track_key, ds_key) {
  outfile <- fig_file(
    fig_subdir,
    cmp_short,
    track_key,
    ds_key,
    "DispEst"
  )

  png(
    outfile,
    width = 2200,
    height = 1700,
    res = figure_dpi
  )

  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  par(
    mar = c(5.0, 5.0, 4.2, 2.0),
    mgp = c(2.6, 0.9, 0),
    cex.main = 1.1,
    cex.lab = 1.0
  )

  DESeq2::plotDispEsts(
    dds,
    main = paste(
      compact_title(pretty_dataset_label(dataset_name)),
      "dispersion estimates"
    )
  )
}

plot_dispersion_panel_for_dataset <- function(df, dataset_name) {
  if (!all(c("baseMean", "dispersion", "final_class") %in% names(df))) {
    return(NULL)
  }

  disp_df <- df[
    is.finite(df$baseMean) &
      !is.na(df$baseMean) &
      is.finite(df$dispersion) &
      !is.na(df$dispersion),
    ,
    drop = FALSE
  ]

  ggplot(
    disp_df,
    aes(
      baseMean,
      dispersion,
      color = final_class,
      shape = final_class
    )
  ) +
    geom_point(
      alpha = 0.55,
      size = 1.15,
      stroke = 0.28
    ) +
    scale_x_log10(labels = label_number(accuracy = 0.1)) +
    scale_y_log10(labels = label_number(accuracy = 0.1)) +
    scale_color_manual(
      values = final_class_colors,
      drop = FALSE,
      name = "Class"
    ) +
    scale_shape_manual(
      values = final_class_shapes,
      drop = FALSE,
      name = "Class"
    ) +
    labs(
      title = compact_title(pretty_dataset_label(dataset_name), width = 42),
      subtitle = "Final dispersion vs mean",
      x = "baseMean",
      y = "Dispersion"
    ) +
    manuscript_theme()
}

plot_dispersion_relationships <- function(df, dataset_name) {
  plots <- list()

  if (all(c("baseMean", "dispersion", "final_class") %in% names(df))) {
    tmp <- df[
      is.finite(df$baseMean) &
        !is.na(df$baseMean) &
        is.finite(df$dispersion) &
        !is.na(df$dispersion),
      ,
      drop = FALSE
    ]

    plots[[length(plots) + 1L]] <- ggplot(
      tmp,
      aes(
        baseMean,
        dispersion,
        color = final_class,
        shape = final_class
      )
    ) +
      geom_point(
        alpha = 0.55,
        size = 1.05,
        stroke = 0.28
      ) +
      scale_x_log10(labels = comma_format()) +
      scale_y_log10() +
      scale_color_manual(values = final_class_colors, drop = FALSE, name = "Class") +
      scale_shape_manual(values = final_class_shapes, drop = FALSE, name = "Class") +
      labs(title = "Final dispersion vs mean", x = "baseMean", y = "Final dispersion") +
      manuscript_theme()
  }

  if (all(c("baseMean", "dispFit", "final_class") %in% names(df))) {
    tmp <- df[
      is.finite(df$baseMean) &
        !is.na(df$baseMean) &
        is.finite(df$dispFit) &
        !is.na(df$dispFit),
      ,
      drop = FALSE
    ]

    plots[[length(plots) + 1L]] <- ggplot(
      tmp,
      aes(
        baseMean,
        dispFit,
        color = final_class,
        shape = final_class
      )
    ) +
      geom_point(
        alpha = 0.55,
        size = 1.05,
        stroke = 0.28
      ) +
      scale_x_log10(labels = comma_format()) +
      scale_y_log10() +
      scale_color_manual(values = final_class_colors, drop = FALSE, name = "Class") +
      scale_shape_manual(values = final_class_shapes, drop = FALSE, name = "Class") +
      labs(title = "Trend dispersion vs mean", x = "baseMean", y = "Trend dispersion") +
      manuscript_theme()
  }

  if (all(c("baseMean", "dispGeneEst", "final_class") %in% names(df))) {
    tmp <- df[
      is.finite(df$baseMean) &
        !is.na(df$baseMean) &
        is.finite(df$dispGeneEst) &
        !is.na(df$dispGeneEst),
      ,
      drop = FALSE
    ]

    plots[[length(plots) + 1L]] <- ggplot(
      tmp,
      aes(
        baseMean,
        dispGeneEst,
        color = final_class,
        shape = final_class
      )
    ) +
      geom_point(
        alpha = 0.55,
        size = 1.05,
        stroke = 0.28
      ) +
      scale_x_log10(labels = comma_format()) +
      scale_y_log10() +
      scale_color_manual(values = final_class_colors, drop = FALSE, name = "Class") +
      scale_shape_manual(values = final_class_shapes, drop = FALSE, name = "Class") +
      labs(title = "Gene-wise dispersion vs mean", x = "baseMean", y = "Gene-wise dispersion") +
      manuscript_theme()
  }

  if (length(plots) == 0L) {
    return(NULL)
  }

  assemble_one_legend_panel(
    plots,
    panel_title = paste(
      compact_title(pretty_dataset_label(dataset_name)),
      "dispersion relationships"
    ),
    ncol = min(2, length(plots))
  )
}

dispersion_residual_section <- function(df, dataset_name, fig_subdir, tab_dir, cmp_short, track_key, ds_key) {
  if (!all(c("dispGeneEst", "dispFit", "dispersion") %in% colnames(df))) {
    return(NULL)
  }

  out <- data.frame(
    feature_id = df$feature_id,
    dispGeneEst = df$dispGeneEst,
    dispFit = df$dispFit,
    dispersion = df$dispersion,
    residual_fit_minus_final = df$dispFit - df$dispersion,
    residual_fit_minus_gene = df$dispFit - df$dispGeneEst,
    residual_final_minus_gene = df$dispersion - df$dispGeneEst,
    stringsAsFactors = FALSE
  )

  p1 <- ggplot(out, aes(residual_fit_minus_final)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, color = "white") +
    labs(title = "dispFit - final", x = "Residual", y = "Count") +
    manuscript_theme()

  p2 <- ggplot(out, aes(residual_fit_minus_gene)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, color = "white") +
    labs(title = "dispFit - gene-wise", x = "Residual", y = "Count") +
    manuscript_theme()

  p3 <- ggplot(out, aes(residual_final_minus_gene)) +
    geom_histogram(bins = 60, fill = plot_palette$histogram, color = "white") +
    labs(title = "final - gene-wise", x = "Residual", y = "Count") +
    manuscript_theme()

  panel <- gridExtra::arrangeGrob(
    p1,
    p2,
    p3,
    ncol = 3,
    top = grid::textGrob(
      paste(
        compact_title(pretty_dataset_label(dataset_name)),
        "dispersion residuals"
      ),
      gp = grid::gpar(fontface = "bold", cex = 1.15)
    )
  )

  save_grob(
    panel,
    fig_file(fig_subdir, cmp_short, track_key, ds_key, "DispResid"),
    width = 15.8,
    height = 5.2
  )

  save_csv(
    out,
    tab_file(tab_dir, cmp_short, track_key, ds_key, "DispResid")
  )

  out
}

# =============================================================================
# SECTION 5 OF 5
# COMPARISON PANELS, BATCH EXECUTION, EXPORTS, AND GIT STAGING
# =============================================================================

compute_dataset_pca_plot <- function(count_df, coldata, dataset_name, preprocessing = c("normalized", "raw_counts")) {
  preprocessing <- match.arg(preprocessing)

  if (preprocessing == "normalized") {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(count_df),
      colData = coldata,
      design = make_design_formula(coldata)
    )

    dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
    dds <- DESeq2::estimateSizeFactors(dds)
    x <- DESeq2::counts(dds, normalized = TRUE)
    preprocessing_label <- "Normalized prior to EVS"
  } else {
    x <- as.matrix(count_df)
    preprocessing_label <- "Raw prior to EVS"
  }

  pca_fit <- stats::prcomp(
    t(x),
    center = TRUE,
    scale. = FALSE,
    rank. = 2
  )

  pca_var <- pca_fit$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)

  pca_df <- data.frame(
    Sample = rownames(pca_fit$x),
    PC1 = pca_fit$x[, 1],
    PC2 = pca_fit$x[, 2],
    Condition = as.character(coldata[rownames(pca_fit$x), "condition"]),
    stringsAsFactors = FALSE
  )

  ggplot(
    pca_df,
    aes(
      PC1,
      PC2,
      label = Sample,
      shape = Condition,
      fill = Condition
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey70") +
    geom_point(size = 2.8, colour = "white", stroke = 0.55) +
    ggrepel::geom_text_repel(
      size = 1.9,
      max.overlaps = 8,
      force = 1.0,
      box.padding = 0.22,
      point.padding = 0.10,
      min.segment.length = 0
    ) +
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition") +
    scale_fill_manual(values = condition_fills, labels = condition_labels, name = "Condition") +
    labs(
      title = compact_title(pretty_dataset_label(dataset_name)),
      subtitle = paste0(
        preprocessing_label,
        " · PC1 = ",
        pca_var_per[1],
        "%; PC2 = ",
        pca_var_per[2],
        "%."
      ),
      x = paste0("PC1 (", pca_var_per[1], "% variance)"),
      y = paste0("PC2 (", pca_var_per[2], "% variance)")
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme()
}

save_cross_dataset_comparison_panels <- function(comparison_name, track_key, analysis_results, cmp_dir, dataset_list, coldata) {
  keys_present <- base::intersect(
    as.character(dataset_key_order),
    as.character(names(analysis_results))
  )

  if (length(keys_present) == 0L) {
    return(invisible(NULL))
  }

  volcano_plots <- lapply(keys_present, function(k) {
    ds_name <- analysis_results[[k]]$summary$dataset_name[1]
    short_title <- paste0(
      comparison_name,
      " ",
      unname(track_short[track_key]),
      " ",
      unname(dataset_short[k])
    )

    plot_final_volcano(
      analysis_results[[k]]$results,
      ds_name,
      short_title = short_title
    )
  })

  volcano_panel <- assemble_one_legend_panel(
    volcano_plots,
    panel_title = paste(
      comparison_name,
      unname(track_short[track_key]),
      "volcano panel"
    ),
    ncol = length(volcano_plots)
  )

  save_grob(
    volcano_panel,
    file.path(
      cmp_dir,
      paste0(
        "Figure_",
        comparison_name,
        "_",
        unname(track_short[track_key]),
        "_Panel_Volcano.png"
      )
    ),
    width = 5.1 * length(volcano_plots),
    height = 6.4
  )

  disp_plots <- lapply(keys_present, function(k) {
    plot_dispersion_panel_for_dataset(
      analysis_results[[k]]$results,
      analysis_results[[k]]$summary$dataset_name[1]
    )
  })

  disp_plots <- Filter(Negate(is.null), disp_plots)

  if (length(disp_plots) > 0L) {
    disp_panel <- assemble_one_legend_panel(
      disp_plots,
      panel_title = paste(
        comparison_name,
        unname(track_short[track_key]),
        "dispersion panel"
      ),
      ncol = length(disp_plots)
    )

    save_grob(
      disp_panel,
      file.path(
        cmp_dir,
        paste0(
          "Figure_",
          comparison_name,
          "_",
          unname(track_short[track_key]),
          "_Panel_Disp.png"
        )
      ),
      width = 5.1 * length(disp_plots),
      height = 6.0
    )
  }

  pca_norm_grobs <- lapply(keys_present, function(k) {
    compute_dataset_pca_plot(
      dataset_list[[k]],
      coldata,
      analysis_results[[k]]$summary$dataset_name[1],
      preprocessing = "normalized"
    )
  })

  pca_raw_grobs <- lapply(keys_present, function(k) {
    compute_dataset_pca_plot(
      dataset_list[[k]],
      coldata,
      analysis_results[[k]]$summary$dataset_name[1],
      preprocessing = "raw_counts"
    )
  })

  pca_panel <- gridExtra::arrangeGrob(
    grobs = c(pca_norm_grobs, pca_raw_grobs),
    ncol = length(keys_present),
    top = grid::textGrob(
      paste(
        comparison_name,
        unname(track_short[track_key]),
        "PCA comparison. Top row normalized prior to EVS; bottom row raw prior to EVS."
      ),
      gp = grid::gpar(fontface = "bold", cex = 1.05)
    )
  )

  save_grob(
    pca_panel,
    file.path(
      cmp_dir,
      paste0(
        "Figure_",
        comparison_name,
        "_",
        unname(track_short[track_key]),
        "_Panel_PCA.png"
      )
    ),
    width = 5.2 * length(keys_present),
    height = 9.2
  )

  hist_grobs <- lapply(keys_present, function(k) {
    plot_empirical_p_histogram(
      analysis_results[[k]]$results,
      analysis_results[[k]]$summary$dataset_name[1],
      analysis_results[[k]]$summary$hc_p_threshold[1]
    )
  })

  hist_panel <- gridExtra::arrangeGrob(
    grobs = hist_grobs,
    ncol = length(hist_grobs),
    top = grid::textGrob(
      paste(
        comparison_name,
        unname(track_short[track_key]),
        "empirical-null p-value histograms"
      ),
      gp = grid::gpar(fontface = "bold", cex = 1.12)
    )
  )

  save_grob(
    hist_panel,
    file.path(
      cmp_dir,
      paste0(
        "Figure_",
        comparison_name,
        "_",
        unname(track_short[track_key]),
        "_Panel_EmpHist.png"
      )
    ),
    width = 5.0 * length(hist_grobs),
    height = 4.8
  )

  summary_table <- dplyr::bind_rows(
    lapply(keys_present, function(k) {
      sm <- analysis_results[[k]]$summary

      data.frame(
        comparison_name = comparison_name,
        track_key = track_key,
        track_label = unname(track_short[track_key]),
        dataset_key = k,
        dataset_label = unname(dataset_key_labels[k]),
        dataset_name = sm$dataset_name[1],
        n_features = sm$n_features[1],
        n_weak_cnh = sm$n_weak_cnh[1],
        n_strong_cnh = sm$n_strong_cnh[1],
        n_standard = sm$n_standard[1],
        n_hbfss = sm$n_hbfss[1],
        n_overlap = sm$n_overlap[1],
        hc_p_threshold = sm$hc_p_threshold[1],
        hbfss_threshold = sm$hbfss_threshold[1],
        stringsAsFactors = FALSE
      )
    })
  )

  save_csv(
    summary_table,
    file.path(
      cmp_dir,
      paste0(
        "Table_",
        comparison_name,
        "_",
        unname(track_short[track_key]),
        "_Panel_Summary.csv"
      )
    )
  )

  invisible(NULL)
}

run_one_evs_track <- function(comparison_name, track_key, count_matrix, coldata, annot_df, cmp_dir) {
  track_dir <- file.path(
    cmp_dir,
    unname(track_short[track_key])
  )

  tab_dir <- file.path(track_dir, "tables")
  fig_raw_dir <- file.path(track_dir, "fig_raw")
  fig_lead_dir <- file.path(track_dir, "fig_lead")
  fig_rem_dir <- file.path(track_dir, "fig_rem")

  for (d in c(track_dir, tab_dir, fig_raw_dir, fig_lead_dir, fig_rem_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  evs <- build_eigenvector_split(
    count_matrix = count_matrix,
    coldata = coldata,
    track_key = track_key
  )

  save_plot(
    plot_pca_scatter(
      evs$fit_trt$pca_fit,
      paste(comparison_name, unname(track_short[track_key])),
      "treatment",
      evs$fit_trt$preprocessing_label
    ),
    file.path(
      track_dir,
      paste0(
        "Figure_",
        comparison_name,
        "_",
        unname(track_short[track_key]),
        "_EVS_Trt_PCA.png"
      )
    ),
    width = 8.3,
    height = 6.2
  )

  save_plot(
    plot_pca_scatter(
      evs$fit_untrt$pca_fit,
      paste(comparison_name, unname(track_short[track_key])),
      "control",
      evs$fit_untrt$preprocessing_label
    ),
    file.path(
      track_dir,
      paste0(
        "Figure_",
        comparison_name,
        "_",
        unname(track_short[track_key]),
        "_EVS_Ctrl_PCA.png"
      )
    ),
    width = 8.3,
    height = 6.2
  )

  save_plot(
    plot_pc1_loading_rank(
      evs$fit_trt$loading_table,
      evs$fit_trt$cutoff,
      paste(comparison_name, unname(track_short[track_key])),
      "treatment",
      top_n_used = evs$fit_trt$top_n_used,
      cutoff_quantile = evs$fit_trt$cutoff_quantile,
      preprocessing_label = evs$fit_trt$preprocessing_label
    ),
    file.path(
      track_dir,
      paste0(
        "Figure_",
        comparison_name,
        "_",
        unname(track_short[track_key]),
        "_EVS_Trt_Rank.png"
      )
    ),
    width = 9.0,
    height = 6.2
  )

  save_plot(
    plot_pc1_loading_rank(
      evs$fit_untrt$loading_table,
      evs$fit_untrt$cutoff,
      paste(comparison_name, unname(track_short[track_key])),
      "control",
      top_n_used = evs$fit_untrt$top_n_used,
      cutoff_quantile = evs$fit_untrt$cutoff_quantile,
      preprocessing_label = evs$fit_untrt$preprocessing_label
    ),
    file.path(
      track_dir,
      paste0(
        "Figure_",
        comparison_name,
        "_",
        unname(track_short[track_key]),
        "_EVS_Ctrl_Rank.png"
      )
    ),
    width = 9.0,
    height = 6.2
  )

  save_grob(
    plot_eigenvector_histograms(
      evs$fit_trt$loading_table,
      evs$fit_trt$cutoff,
      paste(comparison_name, unname(track_short[track_key])),
      "treatment"
    ),
    file.path(
      track_dir,
      paste0(
        "Figure_",
        comparison_name,
        "_",
        unname(track_short[track_key]),
        "_EVS_Trt_Hist.png"
      )
    ),
    width = 14.0,
    height = 4.8
  )

  save_grob(
    plot_eigenvector_histograms(
      evs$fit_untrt$loading_table,
      evs$fit_untrt$cutoff,
      paste(comparison_name, unname(track_short[track_key])),
      "control"
    ),
    file.path(
      track_dir,
      paste0(
        "Figure_",
        comparison_name,
        "_",
        unname(track_short[track_key]),
        "_EVS_Ctrl_Hist.png"
      )
    ),
    width = 14.0,
    height = 4.8
  )

  dataset_list <- list(
    raw_dataset = evs$raw_dataset,
    leading_edge_dataset = evs$leading_edge_dataset,
    remainder_dataset = evs$remainder_dataset
  )

  dataset_fig_dirs <- list(
    raw_dataset = fig_raw_dir,
    leading_edge_dataset = fig_lead_dir,
    remainder_dataset = fig_rem_dir
  )

  analysis_results <- list()

  for (nm in names(dataset_list)) {
    full_dataset_name <- paste(
      comparison_name,
      unname(track_short[track_key]),
      nm,
      sep = "_"
    )

    fig_subdir <- dataset_fig_dirs[[nm]]

    fit <- tryCatch(
      run_core_analysis(
        count_mat = dataset_list[[nm]],
        coldata = coldata,
        dataset_name = full_dataset_name,
        annot_df = annot_df
      ),
      error = function(e) {
        save_csv(
          data.frame(
            comparison_name = comparison_name,
            track_key = track_key,
            dataset_name = full_dataset_name,
            stage = "run_core_analysis",
            error_message = conditionMessage(e),
            stringsAsFactors = FALSE
          ),
          tab_file(tab_dir, comparison_name, track_key, nm, "Error")
        )

        stop(e)
      }
    )

    df <- fit$results

    save_csv(
      df,
      tab_file(tab_dir, comparison_name, track_key, nm, "Results")
    )

    save_csv(
      subset(df, weak_cnh_flag),
      tab_file(tab_dir, comparison_name, track_key, nm, "WeakCNH")
    )

    save_csv(
      subset(df, strong_cnh_flag),
      tab_file(tab_dir, comparison_name, track_key, nm, "StrongCNH")
    )

    save_csv(
      subset(df, standard_flag),
      tab_file(tab_dir, comparison_name, track_key, nm, "Standard")
    )

    save_csv(
      subset(df, hbfss_flag),
      tab_file(tab_dir, comparison_name, track_key, nm, "HBFSS")
    )

    save_csv(
      subset(df, any_overlap),
      tab_file(tab_dir, comparison_name, track_key, nm, "Overlap")
    )

    summary_row <- data.frame(
      comparison_name = comparison_name,
      track_key = track_key,
      track_label = unname(track_short[track_key]),
      dataset_key = nm,
      dataset_name = full_dataset_name,
      n_features = nrow(df),
      hc_p_threshold = fit$hc_p_threshold,
      hbfss_threshold = fit$hbfss_threshold,
      n_weak_cnh = sum(df$weak_cnh_flag, na.rm = TRUE),
      n_strong_cnh = sum(df$strong_cnh_flag, na.rm = TRUE),
      n_standard = sum(df$standard_flag, na.rm = TRUE),
      n_hbfss = sum(df$hbfss_flag, na.rm = TRUE),
      n_overlap = sum(df$any_overlap, na.rm = TRUE),
      n_weak_hbfss_overlap = sum(df$weak_hbfss_overlap, na.rm = TRUE),
      n_strong_hbfss_overlap = sum(df$strong_hbfss_overlap, na.rm = TRUE),
      n_standard_hbfss_overlap = sum(df$standard_hbfss_overlap, na.rm = TRUE),
      top_n_target = top_n_target,
      trt_cutoff_quantile = evs$fit_trt$cutoff_quantile,
      ctrl_cutoff_quantile = evs$fit_untrt$cutoff_quantile,
      leading_edge_n = length(evs$leading_edge_ids),
      remainder_n = length(evs$remainder_ids),
      stringsAsFactors = FALSE
    )

    save_csv(
      summary_row,
      tab_file(tab_dir, comparison_name, track_key, nm, "Summary")
    )

    mean_df <- compute_mean_expression_table(
      dataset_list[[nm]],
      coldata
    )

    mean_panel <- plot_mean_histogram_panel(
      mean_df,
      fit$base_mean_vec,
      full_dataset_name
    )

    save_grob(
      mean_panel,
      fig_file(fig_subdir, comparison_name, track_key, nm, "MeanHist"),
      width = 12.0,
      height = 8.5
    )

    p_vol <- plot_final_volcano(
      df,
      full_dataset_name,
      short_title = paste(
        unname(track_short[track_key]),
        unname(dataset_short[nm])
      )
    )

    p_emp <- plot_empirical_p_histogram(
      df,
      full_dataset_name,
      fit$hc_p_threshold
    )

    p_hc <- plot_hc_profile(
      df,
      full_dataset_name,
      fit$hc_p_threshold
    )

    p_hbfss_dist <- plot_hbfss_distribution(
      df,
      full_dataset_name,
      fit$hbfss_threshold
    )

    p_ma <- plot_shrunken_ma(
      df,
      full_dataset_name
    )

    save_plot(
      p_vol,
      fig_file(fig_subdir, comparison_name, track_key, nm, "Volcano"),
      width = 8.4,
      height = 6.4
    )

    save_plot(
      p_emp,
      fig_file(fig_subdir, comparison_name, track_key, nm, "EmpHist"),
      width = 8.4,
      height = 6.0
    )

    if (!is.null(p_hc)) {
      save_plot(
        p_hc,
        fig_file(fig_subdir, comparison_name, track_key, nm, "HCProfile"),
        width = 8.4,
        height = 6.0
      )
    }

    save_plot(
      p_hbfss_dist,
      fig_file(fig_subdir, comparison_name, track_key, nm, "HBFSSDist"),
      width = 8.4,
      height = 6.0
    )

    if (!is.null(p_ma)) {
      save_plot(
        p_ma,
        fig_file(fig_subdir, comparison_name, track_key, nm, "MA"),
        width = 8.4,
        height = 6.0
      )
    }

    plot_dispersion_cloud(
      fit$dds,
      full_dataset_name,
      fig_subdir,
      comparison_name,
      track_key,
      nm
    )

    disp_panel <- plot_dispersion_relationships(
      df,
      full_dataset_name
    )

    if (!is.null(disp_panel)) {
      save_grob(
        disp_panel,
        fig_file(fig_subdir, comparison_name, track_key, nm, "DispRel"),
        width = 12.0,
        height = 7.6
      )
    }

    dispersion_residual_section(
      df,
      full_dataset_name,
      fig_subdir,
      tab_dir,
      comparison_name,
      track_key,
      nm
    )

    summary_components <- Filter(
      Negate(is.null),
      list(
        p_vol,
        p_emp,
        p_hc,
        p_hbfss_dist,
        p_ma
      )
    )

    summary_panel <- do.call(
      gridExtra::arrangeGrob,
      c(
        summary_components,
        list(
          ncol = 2,
          top = grid::textGrob(
            paste(
              pretty_dataset_label(full_dataset_name),
              "summary panel"
            ),
            gp = grid::gpar(fontface = "bold", cex = 1.15)
          )
        )
      )
    )

    save_grob(
      summary_panel,
      fig_file(fig_subdir, comparison_name, track_key, nm, "Summary"),
      width = 14.0,
      height = 12.0
    )

    analysis_results[[nm]] <- list(
      dds = fit$dds,
      results = df,
      summary = summary_row,
      dataset_mat = dataset_list[[nm]],
      fig_subdir = fig_subdir
    )
  }

  tryCatch(
    {
      save_cross_dataset_comparison_panels(
        comparison_name = comparison_name,
        track_key = track_key,
        analysis_results = analysis_results,
        cmp_dir = track_dir,
        dataset_list = dataset_list,
        coldata = coldata
      )
    },
    error = function(e) {
      warning(
        paste0(
          comparison_name,
          " ",
          track_key,
          ": cross-dataset panel generation failed: ",
          conditionMessage(e)
        )
      )
    }
  )

  if (run_twas_overlap) {
    twas_genes <- clean_gene_set(TWAS_data$gene_symbol)

    get_twas_overlap <- function(result_df, dataset_nm, cmp_name, out_dir, mode = c("union", "standard_only", "hbfss_only")) {
      mode <- match.arg(mode)

      sig_df <- switch(
        mode,
        union = subset(result_df, standard_flag | hbfss_flag),
        standard_only = subset(result_df, standard_flag),
        hbfss_only = subset(result_df, hbfss_flag)
      )

      sig_df$gene_symbol_clean <- tolower(trimws(sig_df$gene_symbol))

      overlap_df <- subset(
        sig_df,
        gene_symbol_clean %in% twas_genes
      )

      summary_df <- data.frame(
        comparison_name = cmp_name,
        dataset_name = dataset_nm,
        selection_mode = mode,
        n_selected_features = nrow(sig_df),
        n_overlap_features = nrow(overlap_df),
        n_overlap_genes = length(unique(overlap_df$gene_symbol_clean)),
        stringsAsFactors = FALSE
      )

      mode_tag <- switch(
        mode,
        union = "Union",
        standard_only = "StdOnly",
        hbfss_only = "HBFSSOnly"
      )

      write.csv(
        overlap_df,
        file.path(
          out_dir,
          paste0("Table_", dataset_nm, "_TWAS_", mode_tag, ".csv")
        ),
        row.names = FALSE
      )

      write.csv(
        summary_df,
        file.path(
          out_dir,
          paste0("Table_", dataset_nm, "_TWAS_", mode_tag, "_Summary.csv")
        ),
        row.names = FALSE
      )

      summary_df
    }

    twas_summaries <- dplyr::bind_rows(
      lapply(names(analysis_results), function(nm) {
        full_nm <- paste(
          comparison_name,
          unname(track_short[track_key]),
          nm,
          sep = "_"
        )

        dplyr::bind_rows(
          get_twas_overlap(analysis_results[[nm]]$results, full_nm, comparison_name, tab_dir, "union"),
          get_twas_overlap(analysis_results[[nm]]$results, full_nm, comparison_name, tab_dir, "standard_only"),
          get_twas_overlap(analysis_results[[nm]]$results, full_nm, comparison_name, tab_dir, "hbfss_only")
        )
      })
    )

    save_csv(
      twas_summaries,
      file.path(
        tab_dir,
        paste0(
          "Table_",
          comparison_name,
          "_",
          unname(track_short[track_key]),
          "_TWAS_Summary.csv"
        )
      )
    )
  }

  sm_list <- lapply(analysis_results, `[[`, "summary")
  sm_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, sm_list)

  if (length(sm_list) == 0L) {
    return(data.frame())
  }

  dplyr::bind_rows(sm_list)
}

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  cmp_dir <- file.path(
    output_dir,
    comparison_name
  )

  dir.create(
    cmp_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  track_keys <- c(
    "normalized_evs",
    "raw_evs"
  )

  track_summaries <- list()

  for (track_key in track_keys) {
    message(
      "Running ",
      comparison_name,
      " ",
      unname(track_short[track_key])
    )

    out <- run_one_evs_track(
      comparison_name = comparison_name,
      track_key = track_key,
      count_matrix = count_matrix,
      coldata = coldata,
      annot_df = annot_df,
      cmp_dir = cmp_dir
    )

    if (!is.null(out) && nrow(out) > 0L) {
      track_summaries[[track_key]] <- out
    }
  }

  if (length(track_summaries) == 0L) {
    return(data.frame())
  }

  dplyr::bind_rows(track_summaries)
}

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
      failed_comparisons[[cmp]] <<- data.frame(
        comparison_name = cmp,
        error_message = conditionMessage(e),
        stringsAsFactors = FALSE
      )

      NULL
    }
  )

  if (!is.null(out) && nrow(out) > 0L) {
    all_summaries_list[[cmp]] <- out
  }
}

all_summaries <- if (length(all_summaries_list) > 0L) {
  dplyr::bind_rows(all_summaries_list)
} else {
  data.frame()
}

if (nrow(all_summaries) > 0L) {
  save_csv(
    all_summaries,
    file.path(output_dir, "Table_Overall_Summary.csv")
  )
}

if (length(failed_comparisons) > 0L) {
  failed_df <- dplyr::bind_rows(failed_comparisons)

  save_csv(
    failed_df,
    file.path(output_dir, "Table_Failed.csv")
  )
}

manifest_df <- write_export_manifest(output_dir)

finalize_git_exports(
  root_dir = output_dir,
  repo_root = repo_root
)

cat("\n=====================================================\n")
cat("Pipeline complete.\n")
cat("Repository root:\n")
cat(repo_root, "\n")
cat("Output directory:\n")
cat(output_dir, "\n")
cat("Exported files:\n")
cat(nrow(manifest_df), "\n")
cat("=====================================================\n\n")

if (nrow(manifest_df) == 0L) {
  stop("Pipeline finished but no exported files were found in output_dir.")
}

if (length(failed_comparisons) > 0L) {
  stop("One or more comparisons failed. See Table_Failed.csv.")
}

if (nrow(all_summaries) > 0L) {
  print(all_summaries)
}
