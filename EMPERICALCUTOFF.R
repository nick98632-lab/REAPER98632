#!/usr/bin/env Rscript

# =============================================================================
# FINAL MANUSCRIPT PIPELINE — JOURNAL REVIEW VERSION
# EVS + DESeq2 + empirical-null calibration + higher criticism + HBFSS
# =============================================================================
#
# This script is written as a manuscript supplement and a runnable analysis
# pipeline. The scientific explanations are intentionally extensive, but the
# executable export path is intentionally lean: the final run writes result
# tables, compact summary tables, and journal-ready manuscript figures only.
# The figure set is restricted to discovery volcanoes, PCA support, EVS loading
# support, empirical/HBFSS support, and discovery-count summaries.
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
#      The full comparison matrix is normalized with DESeq2 median-of-ratios
#      size factors before condition-specific PCA / PC1-loading ranking.
#
#   2. RawEVS
#      Condition-specific PCA / PC1-loading ranking is performed directly on
#      raw counts, without normalization before eigenvector splitting.
#
# For each track, the EVS step ONLY selects feature IDs for the leading-edge
# and remainder split. The downstream DESeq2 input is always the raw count
# subset corresponding to those feature IDs. DESeq2 is then fit independently
# to Raw, Lead, and Rem datasets, including its own size-factor normalization
# after the split. Thus NormEVS versus RawEVS compares only the data used for
# PC1-loading ranking before the split, not the downstream DESeq2 modeling or
# normalization framework.
#
# -----------------------------------------------------------------------------
# DATASETS PRODUCED PER COMPARISON AND PER EVS TRACK
# -----------------------------------------------------------------------------
#
#   Raw dataset:
#     the full comparison count matrix.
#
#   Leading-edge dataset:
#     exactly the top 5,000 PAS features by absolute PC1 loading are selected
#     independently from the treatment and control condition eigenvectors. Joint
#     passes (in both top-5,000 sets) and disjoint passes (unique to either set)
#     are combined to form the leading edge.
#
#   Remainder dataset:
#     all PAS features not present in either condition-specific top-5,000 set.
#
# EVS follows the supplied Figure 1 / Figure 1B workflow: median-of-ratios
# normalization precedes PCA in the normalized track; PCA is performed with
# prcomp, PC1 loadings are converted to absolute values, and the highest 5,000
# loadings are selected separately in each condition. No adaptive top-N
# reduction and no pre-PCA log2(x + 1) transformation are used. RawEVS applies
# the same PCA/ranking/splitting rule directly to the raw-count matrix.
#
# -----------------------------------------------------------------------------
# DESEQ2, EMPIRICAL-NULL CALIBRATION, AND FINAL DECISION RULES
# -----------------------------------------------------------------------------
# DESeq2 is fit independently to the Original, Leading-edge, and Remainder
# datasets. For the EVS-derived datasets, the split is defined either from
# DESeq2-normalized counts (NormEVS) or raw counts (RawEVS), but downstream
# DESeq2 always receives the corresponding raw-count subset and performs its own
# size-factor normalization.
#
# The finite DESeq2 Wald statistics are calibrated with the Strimmer fdrtool
# empirical-null model. The fdrtool empirical-null p-values are the p-values used
# in the HBFSS calculation. Higher criticism is applied to the sorted empirical
# p-values:
#
#   HCp = fdrtool::hc.thresh(sort(empirical_p))
#
# The dataset-specific HBFSS cutoff is derived from HCp and the manuscript
# effect boundary c = 1:
#
#   Htau = |-log10(HCp)| * c
#
# HBFSS is calculated from the apeglm-shrunken log2 fold change:
#
#   HBFSS_i = |apeglm_LFC_i| * [-log10(empirical_p_i)]
#
# A PAS is HBFSS-significant when:
#
#   HBFSS_i > Htau
#
# There is NO additional padj/q-value gate on HBFSS and HCp is NOT applied as a
# second independent significance requirement. HCp is the empirical threshold
# used to derive Htau. The volcano therefore shows both the HCp reference line
# and the HBFSS hyperbolic decision boundary:
#
#   y = Htau / |x|
#
# where x is the apeglm-shrunken log2 fold change and
# y = -log10(empirical_p).
#
# -----------------------------------------------------------------------------
# DESEQ2 EFFECT DEFINITIONS
# -----------------------------------------------------------------------------
# Let beta denote the apeglm-shrunken log2 fold change used for manuscript
# effect interpretation and let c = lfc_boundary = 1.0.
#
# Standard DESeq2:
#   DESeq2 Wald test with explicit Benjamini-Hochberg adjustment.
#   Significant when padj < 0.10 and |apeglm_LFC| >= 1.
#
# Strong composite-null hypothesis (Strong CNH):
#   H0,strong : |beta| <= 1
#   HA,strong : |beta| >  1
#   DESeq2: altHypothesis = "greaterAbs", lfcThreshold = 1
#   Pass when resGA_padj < 0.10 and |apeglm_LFC| >= 1.
#
# Weak composite-null hypothesis (Weak CNH):
#   H0,weak : |beta| >= 1
#   HA,weak : |beta| <  1
#   DESeq2: altHypothesis = "lessAbs", lfcThreshold = 1
#   Pass when resLA_padj < 0.10 and |apeglm_LFC| < 1.
#
# A Weak CNH pass is NOT treated as a final significant weak-effect discovery by
# itself. Final weak-effect significance requires independent HBFSS support:
#
#   weak_significant = weak_cnh_flag & hbfss_flag
#
# Thus weak-effect PASs can be inspected as Weak CNH passes, while only the
# HBFSS-overlapping subset enters the significant-site tables and Weak discovery
# counts.
#
# -----------------------------------------------------------------------------
# HBFSS / DESEQ2 OVERLAP AND VOLCANO DISPLAY
# -----------------------------------------------------------------------------
# Method flags are deliberately kept non-mutually-exclusive for counting:
#
#   standard_flag
#   strong_cnh_flag
#   weak_cnh_flag
#   hbfss_flag
#
# Overlap is:
#
#   hbfss_flag & (standard_flag | strong_cnh_flag | weak_cnh_flag)
#
# Final significant sites are those supported by HBFSS, Standard DESeq2, or
# Strong CNH. Weak CNH contributes to final significance only through
# weak_cnh_flag & hbfss_flag.
#
# Volcano display preserves this logic. All HBFSS-positive points, including
# HBFSS overlaps, use the HBFSS purple color. Marker shape identifies the
# overlapping DESeq2 class:
#
#   HBFSS-only        = purple star
#   HBFSS ∩ Weak      = purple triangle
#   HBFSS ∩ Strong    = purple square
#   HBFSS ∩ Standard  = purple diamond
#
# Non-HBFSS Strong and Standard discoveries retain red and green markers.
# Weak-CNH-only passes are shown in blue for transparency but are explicitly
# labeled as unsupported weak passes and are excluded from final significant
# tables/counts. Background is gray.
#
# Each volcano labels at most the top 20 FINAL significant PASs (gene symbol
# when available, otherwise PAS/feature ID). The concise plot annotation reports:
#
#   H = total HBFSS
#   Std = Standard DESeq2
#   Str = Strong CNH
#   Wk = HBFSS-supported Weak CNH
#   Ovlp = HBFSS ∩ any DESeq2 class
#   HCp = higher-criticism empirical-p threshold
#   Htau = dataset-specific HBFSS cutoff
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

required_packages <- c(
  "DESeq2",
  "apeglm",
  "fdrtool",
  "ggplot2",
  "ggrepel",
  "dplyr",
  "gridExtra",
  "grid",
  "scales",
  "grDevices",
  "S4Vectors"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Required R package(s) are not installed: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this manuscript pipeline."
  )
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

# DESeq2 standard and composite-null tests use explicit Benjamini-Hochberg
# adjustment at 10%. HBFSS does not use DESeq2 padj/q-values as a decision gate.
alpha_level <- 0.10
standard_alpha_level <- 0.10

lfc_boundary <- 1.0

# EVS method specification: select exactly the 5,000 highest absolute PC1
# loadings independently from each condition eigenvector.
top_n_target <- 5000L

figure_dpi <- 320

base_theme_size <- 10

n_top_labels_volcano <- 20L

# Manuscript export behavior. When TRUE, the script writes only the focused
# paper-ready figure panels into exports/manuscript_final_clean/manuscript_figures.
# Tables are deliberately concise: significant sites plus compact method/count
# summaries, with the unsplit Original dataset represented once.
EXPORT_ONLY_PAPER_FIGURES <- TRUE
EXPORT_SUPPORT_FIGURES <- TRUE

# Git export behavior.
# The script always writes exports into the repository exports folder when it
# can detect the repository. Staging is enabled by default. Commit and push are
# disabled by default so that the user can review output before publishing.
GIT_STAGE_EXPORTS <- TRUE
GIT_COMMIT_EXPORTS <- FALSE
GIT_PUSH_EXPORTS <- FALSE
GIT_COMMIT_MESSAGE <- "Refresh final manuscript EVS HBFSS exports"

# -----------------------------------------------------------------------------
# Decision-rule audit notes for manuscript review
# -----------------------------------------------------------------------------
# DESeq2 padj is explicitly Benjamini-Hochberg adjusted at 10% for Standard,
# Strong CNH, and Weak CNH tests. This is BH FDR control (not Bonferroni/Holm).
#
# HBFSS is mathematically separate from those DESeq2 adjusted-p decisions:
#
#   standard_flag = padj < 0.10 and |apeglm_LFC| >= 1
#   strong_cnh_flag = greaterAbs padj < 0.10 and |apeglm_LFC| >= 1
#   weak_cnh_flag = lessAbs padj < 0.10 and |apeglm_LFC| < 1
#   hbfss_flag = HBFSS > Htau
#   weak_significant_flag = weak_cnh_flag & hbfss_flag
#
# HCp is used to derive Htau. It is displayed as a reference threshold but is
# not imposed as a second HBFSS gate. HBFSS/DESeq2 overlaps remain explicit.
#
# -----------------------------------------------------------------------------
# Palette
# -----------------------------------------------------------------------------

plot_palette <- list(
  background = "#BDBDBD",
  threshold = "#A65628",
  hc = "#A65628",
  hbfss_line = "#6A3D9A",
  weak = "#4EA3F1",
  strong = "#E31A1C",
  standard = "#33A02C",
  hbfss = "#6A3D9A",
  overlap = "#54278F",
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

paper_fig_dir <- file.path(output_dir, "manuscript_figures")
dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)

paper_registry <- list()

registry_key <- function(comparison_name, track_key) {
  paste(comparison_name, track_key, sep = "__")
}

clean_manuscript_figure_exports <- function() {
  old_files <- list.files(
    paper_fig_dir,
    pattern = "^Figure_Manuscript_.*\\.png$",
    full.names = TRUE
  )

  if (length(old_files) > 0L) {
    unlink(old_files)
  }

  invisible(TRUE)
}

clean_manuscript_table_exports <- function() {
  old_files <- list.files(
    output_dir,
    pattern = "^Table_.*\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(old_files) > 0L) {
    unlink(old_files)
  }

  old_methods <- file.path(output_dir, "Methods_Statistical_Decision_Rules.txt")
  if (file.exists(old_methods)) {
    unlink(old_methods)
  }

  invisible(TRUE)
}


should_write_figure <- function(path) {
  if (!isTRUE(EXPORT_ONLY_PAPER_FIGURES)) {
    return(TRUE)
  }

  target_dir <- normalizePath(
    paper_fig_dir,
    winslash = "/",
    mustWork = FALSE
  )

  path_dir <- normalizePath(
    dirname(path),
    winslash = "/",
    mustWork = FALSE
  )

  startsWith(path_dir, target_dir)
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

save_grob <- function(g, path, width = 14.0, height = 8.5, dpi = figure_dpi, bg = "white") {
  if (!should_write_figure(path)) {
    return(invisible(NULL))
  }

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

  if (length(parts) >= 5L && parts[3] %in% unname(track_short)) {
    track_label <- parts[3]
    dataset_key <- paste(parts[4:length(parts)], collapse = "_")
    return(paste(comparison_name, track_label, pretty_dataset_type(dataset_key), sep = " | "))
  }

  dataset_key <- paste(parts[3:length(parts)], collapse = "_")

  paste(comparison_name, pretty_dataset_type(dataset_key), sep = " | ")
}

clean_gene_set <- function(x) {
  unique(tolower(trimws(x[!is.na(x) & x != ""])))
}

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
  top_n <- as.integer(top_n)

  if (n_total == 0L) {
    stop("resolve_top_n_cutoff() received an empty vector.")
  }

  if (!is.finite(top_n) || is.na(top_n) || top_n < 1L) {
    stop("top_n must be a positive integer.")
  }

  # The EVS method is explicitly fixed at the requested top-N. Do not shrink N
  # to preserve a remainder or replace the rank rule with a loading-value cutoff.
  if (n_total < top_n) {
    stop(
      "EVS requires exactly ", top_n,
      " PAS features per condition, but only ", n_total,
      " features are available for ranking."
    )
  }

  list(
    top_n_actual = top_n,
    cutoff_value = sorted_values_desc[top_n],
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

  # Retain the hc.thresh result whenever it is a valid probability. No
  # additional near-1 exclusion is imposed. HCp = 1 is allowed and yields
  # Htau = 0 under the stated HBFSS rule.
  if (!is.finite(out) ||
      is.na(out) ||
      out <= 0 ||
      out > 1) {
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
  "BG",
  "Weak pass",
  "Strong",
  "Std",
  "HBFSS",
  "HBFSS∩Weak",
  "HBFSS∩Strong",
  "HBFSS∩Std"
)

final_class_colors <- c(
  "BG" = plot_palette$background,
  "Weak pass" = plot_palette$weak,
  "Strong" = plot_palette$strong,
  "Std" = plot_palette$standard,
  "HBFSS" = plot_palette$hbfss,
  "HBFSS∩Weak" = plot_palette$hbfss,
  "HBFSS∩Strong" = plot_palette$hbfss,
  "HBFSS∩Std" = plot_palette$hbfss
)

final_class_shapes <- c(
  "BG" = 16,
  "Weak pass" = 17,
  "Strong" = 15,
  "Std" = 18,
  "HBFSS" = 8,
  "HBFSS∩Weak" = 17,
  "HBFSS∩Strong" = 15,
  "HBFSS∩Std" = 18
)

final_class_sizes <- c(
  "BG" = 0.52,
  "Weak pass" = 0.85,
  "Strong" = 1.20,
  "Std" = 1.10,
  "HBFSS" = 1.30,
  "HBFSS∩Weak" = 1.35,
  "HBFSS∩Strong" = 1.35,
  "HBFSS∩Std" = 1.35
)

final_class_alphas <- c(
  "BG" = 0.24,
  "Weak pass" = 0.58,
  "Strong" = 0.94,
  "Std" = 0.90,
  "HBFSS" = 0.96,
  "HBFSS∩Weak" = 1.00,
  "HBFSS∩Strong" = 1.00,
  "HBFSS∩Std" = 1.00
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

coerce_count_column <- function(x) {
  suppressWarnings(
    as.numeric(
      gsub(
        ",",
        "",
        trimws(as.character(x)),
        fixed = TRUE
      )
    )
  )
}

for (sid in meta_all$id) {
  WTTS_Seq[[sid]] <- coerce_count_column(WTTS_Seq[[sid]])
}

WTTS_Seq <- WTTS_Seq[
  !is.na(WTTS_Seq$OrigID) & nzchar(trimws(WTTS_Seq$OrigID)),
  ,
  drop = FALSE
]

sample_na <- rowSums(is.na(WTTS_Seq[, meta_all$id, drop = FALSE])) > 0

if (any(sample_na)) {
  message("Removing ", sum(sample_na), " rows with missing or nonnumeric sample counts.")
}

WTTS_Seq <- WTTS_Seq[!sample_na, , drop = FALSE]

WTTS_Seq$feature_id <- make.unique(as.character(WTTS_Seq$OrigID), sep = "_dup")
rownames(WTTS_Seq) <- WTTS_Seq$feature_id

OrigID_Symbol <- data.frame(
  feature_id = WTTS_Seq$feature_id,
  orig_id = WTTS_Seq$OrigID,
  gene_symbol = WTTS_Seq$Symbol,
  stringsAsFactors = FALSE
)

OrigID_Symbol$feature_id <- as.character(OrigID_Symbol$feature_id)
OrigID_Symbol$orig_id <- as.character(OrigID_Symbol$orig_id)
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
  rownames(count_sub) <- rownames(WTTS_Seq)

  stopifnot(all(colnames(count_sub) == rownames(coldata)))

  list(
    comparison_name = comparison_name,
    count_matrix = coerce_raw_count_matrix_for_deseq2(
      count_sub,
      context = paste0(comparison_name, " raw comparison matrix")
    ),
    coldata = coldata
  )
}


coerce_raw_count_matrix_for_deseq2 <- function(count_mat, context = "count matrix") {
  original_rownames <- rownames(count_mat)
  original_colnames <- colnames(count_mat)

  count_mat <- as.matrix(count_mat)

  if (!is.numeric(count_mat)) {
    suppressWarnings(storage.mode(count_mat) <- "numeric")
  }

  if (!is.numeric(count_mat)) {
    stop(context, " must be numeric raw counts before DESeq2.")
  }

  if (any(!is.finite(count_mat) | is.na(count_mat))) {
    stop(context, " contains NA or non-finite values after numeric coercion.")
  }

  if (any(count_mat < 0)) {
    stop(context, " contains negative values; DESeq2 requires non-negative counts.")
  }

  rounded <- round(count_mat)

  if (any(abs(count_mat - rounded) > 1e-6)) {
    warning(context, " contained non-integer values; values were rounded for DESeq2.")
  }

  if (any(rounded > .Machine$integer.max)) {
    stop(context, " contains counts larger than R integer storage can represent.")
  }

  storage.mode(rounded) <- "integer"
  rownames(rounded) <- original_rownames
  colnames(rounded) <- original_colnames

  rounded
}

make_rank_matrix_for_track <- function(count_matrix, coldata, track_key) {
  track_key <- match.arg(track_key, c("normalized_evs", "raw_evs"))

  if (track_key == "raw_evs") {
    return(
      list(
        rank_matrix = as.data.frame(count_matrix),
        preprocessing_label = "Raw counts prior to EVS; PCA is performed directly on raw counts"
      )
    )
  }

  dds_init <- DESeq2::DESeqDataSetFromMatrix(
    countData = coerce_raw_count_matrix_for_deseq2(
      count_matrix,
      context = "NormEVS pre-split full comparison matrix"
    ),
    colData = coldata,
    design = make_design_formula(coldata)
  )

  dds_init <- DESeq2::estimateSizeFactors(dds_init)

  norm_counts <- as.data.frame(
    DESeq2::counts(dds_init, normalized = TRUE)
  )

  list(
    rank_matrix = norm_counts,
    preprocessing_label = "Median-of-ratios normalized counts prior to EVS; PCA is performed directly on normalized counts"
  )
}

compute_pc1_loading_table <- function(value_df, sample_names, top_n = top_n_target, preprocessing_label = "Normalized prior to EVS") {
  x <- as.matrix(value_df[, sample_names, drop = FALSE])
  storage.mode(x) <- "numeric"

  if (ncol(x) < 2L) {
    stop("EVS PCA requires at least two samples in each condition.")
  }

  # Supplied EVS method: run PCA directly on the pre-split matrix for the
  # condition, extract the PC1 feature eigenvector/loading, take absolute values,
  # rank descending, and select the fixed top-5,000. No log transformation is
  # applied before PCA.
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

  loading_tbl$split_class <- ifelse(
    loading_tbl$rank <= cutoff_info$top_n_actual,
    "high_loading",
    "background_loading"
  )

  list(
    pca_fit = pca_fit,
    loading_table = loading_tbl,
    cutoff = cutoff_info$cutoff_value,
    top_n_used = cutoff_info$top_n_actual,
    preprocessing_label = preprocessing_label
  )
}

build_eigenvector_split <- function(count_matrix, coldata, track_key) {
  track_key <- match.arg(track_key, c("normalized_evs", "raw_evs"))

  # The EVS matrix is used only to choose feature IDs. The normalized track
  # follows the supplied workflow: median-of-ratios normalization of the full
  # RT/ZT comparison matrix, followed by condition-specific PCA/PC1 eigenvectors.
  # RawEVS repeats the same PCA/ranking rule without pre-split normalization.
  # Downstream DESeq2 always receives the corresponding raw-count subset.
  raw_count_matrix <- coerce_raw_count_matrix_for_deseq2(
    count_matrix,
    context = paste0(track_key, " full comparison matrix before EVS")
  )

  rank_obj <- make_rank_matrix_for_track(
    count_matrix = raw_count_matrix,
    coldata = coldata,
    track_key = track_key
  )

  rank_matrix <- rank_obj$rank_matrix
  preprocessing_label <- rank_obj$preprocessing_label

  sample_ids <- colnames(raw_count_matrix)
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
    fit_trt$loading_table$feature_id[
      fit_trt$loading_table$rank <= top_n_target
    ]
  )

  untrt_high <- as.character(
    fit_untrt$loading_table$feature_id[
      fit_untrt$loading_table$rank <= top_n_target
    ]
  )

  if (length(trt_high) != top_n_target || length(untrt_high) != top_n_target) {
    stop("EVS failed to select exactly the top 5,000 PAS features per condition.")
  }

  # Supplied Figure 1 / Figure 1B classification.
  joint_ids <- intersect(trt_high, untrt_high)
  disjoint_trt_ids <- setdiff(trt_high, untrt_high)
  disjoint_untrt_ids <- setdiff(untrt_high, trt_high)

  leading_edge_ids <- union(trt_high, untrt_high)
  remainder_ids <- setdiff(rownames(raw_count_matrix), leading_edge_ids)

  if (length(leading_edge_ids) == 0L) {
    stop("Leading-edge dataset is empty. Check sample mapping or EVS inputs.")
  }

  if (length(remainder_ids) == 0L) {
    stop("Remainder dataset is empty after the fixed top-5,000 EVS split.")
  }

  list(
    track_key = track_key,
    preprocessing_label = preprocessing_label,
    rank_matrix = rank_matrix,
    fit_trt = fit_trt,
    fit_untrt = fit_untrt,
    trt_top_ids = trt_high,
    untrt_top_ids = untrt_high,
    joint_ids = joint_ids,
    disjoint_trt_ids = disjoint_trt_ids,
    disjoint_untrt_ids = disjoint_untrt_ids,
    leading_edge_ids = leading_edge_ids,
    remainder_ids = remainder_ids,
    downstream_deseq2_input = "raw count subsets; DESeq2 size-factor normalization is run after EVS split",
    raw_dataset = raw_count_matrix,
    leading_edge_dataset = raw_count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset = raw_count_matrix[remainder_ids, , drop = FALSE]
  )
}

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
  design_formula <- make_design_formula(coldata)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = coerce_raw_count_matrix_for_deseq2(
      count_mat,
      context = paste0(dataset_name, " DESeq2 input after EVS split")
    ),
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
    alpha = standard_alpha_level,
    pAdjustMethod = "BH"
  )

  res_strong <- DESeq2::results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "greaterAbs",
    alpha = alpha_level,
    pAdjustMethod = "BH"
  )

  res_weak <- DESeq2::results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "lessAbs",
    alpha = alpha_level,
    pAdjustMethod = "BH"
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
    type = "apeglm"
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
      hc_p_threshold_dataset > 1) {
    hbfss_threshold_dataset <- NA_real_
    res_df$HBFSS_core_pass <- FALSE
  } else {
    hbfss_threshold_dataset <- abs(log10(hc_p_threshold_dataset)) * lfc_boundary
    res_df$HBFSS_core_pass <- !is.na(res_df$HBFSS) &
      is.finite(res_df$HBFSS) &
      res_df$HBFSS > hbfss_threshold_dataset
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

  # Final decision flags.
  # DESeq2 Standard/Strong/Weak use BH-adjusted p-values at 10%. HBFSS is
  # significant from the score cutoff alone (HBFSS > Htau); HCp is already
  # incorporated into Htau and is not an additional gate.
  res_df$weak_cnh_flag <- !is.na(res_df$resLA_padj) &
    res_df$resLA_padj < alpha_level &
    !is.na(res_df$lfc_shrunk) &
    abs(res_df$lfc_shrunk) < lfc_boundary

  res_df$strong_cnh_flag <- !is.na(res_df$resGA_padj) &
    res_df$resGA_padj < alpha_level &
    !is.na(res_df$lfc_shrunk) &
    abs(res_df$lfc_shrunk) >= lfc_boundary

  res_df$standard_significant <- !is.na(res_df$padj) &
    res_df$padj < standard_alpha_level &
    !is.na(res_df$lfc_shrunk) &
    abs(res_df$lfc_shrunk) >= lfc_boundary

  res_df$standard_flag <- res_df$standard_significant

  res_df$hbfss_flag <- !is.na(res_df$HBFSS_core_pass) &
    res_df$HBFSS_core_pass

  res_df$HBFSS_significant <- res_df$hbfss_flag

  res_df$weak_hbfss_overlap <- res_df$weak_cnh_flag & res_df$hbfss_flag
  res_df$strong_hbfss_overlap <- res_df$strong_cnh_flag & res_df$hbfss_flag
  res_df$standard_hbfss_overlap <- res_df$standard_flag & res_df$hbfss_flag

  res_df$weak_significant_flag <- res_df$weak_hbfss_overlap
  res_df$any_overlap <- res_df$hbfss_flag &
    (res_df$weak_cnh_flag | res_df$strong_cnh_flag | res_df$standard_flag)

  # A final significant site is supported by HBFSS, Standard DESeq2, or Strong
  # CNH. Weak CNH alone is shown for context but is not a final discovery.
  res_df$final_significant_flag <- res_df$hbfss_flag |
    res_df$standard_flag |
    res_df$strong_cnh_flag

  # Mutually exclusive display classes. Every HBFSS-positive site is purple;
  # marker shape identifies its DESeq2 overlap. Weak-only passes remain blue and
  # are explicitly excluded from final significance.
  res_df$final_class <- "BG"
  res_df$final_class[res_df$weak_cnh_flag & !res_df$hbfss_flag] <- "Weak pass"
  res_df$final_class[res_df$standard_flag & !res_df$hbfss_flag] <- "Std"
  res_df$final_class[res_df$strong_cnh_flag & !res_df$hbfss_flag] <- "Strong"
  res_df$final_class[res_df$hbfss_flag] <- "HBFSS"
  res_df$final_class[res_df$standard_hbfss_overlap] <- "HBFSS∩Std"
  res_df$final_class[res_df$strong_hbfss_overlap] <- "HBFSS∩Strong"
  res_df$final_class[res_df$weak_hbfss_overlap] <- "HBFSS∩Weak"

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
      "dispOutlier"
    ),
    colnames(mm)
  )

  disp_df <- mm[, disp_cols_available, drop = FALSE]

  annot_df$feature_id <- as.character(annot_df$feature_id)
  if (!"gene_symbol" %in% names(annot_df)) {
    annot_df$gene_symbol <- NA_character_
  }
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
    "orig_id",
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
    "weak_significant_flag",
    "final_significant_flag",
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
    as.character(df$feature_id)
  )

  df$has_valid_gene_symbol <- !is.na(df$gene_symbol_plot) &
    grepl("[A-Za-z0-9]", trimws(df$gene_symbol_plot))

  df$final_class <- factor(
    as.character(df$final_class),
    levels = final_class_levels
  )

  draw_order <- c(
    "BG" = 1,
    "Weak pass" = 2,
    "Std" = 3,
    "Strong" = 4,
    "HBFSS" = 5,
    "HBFSS∩Std" = 6,
    "HBFSS∩Strong" = 7,
    "HBFSS∩Weak" = 8
  )

  df$final_draw_order <- unname(draw_order[as.character(df$final_class)])
  df$final_draw_order[is.na(df$final_draw_order)] <- 1

  df <- df[
    order(df$final_draw_order, df$neglog10_empirical_p),
    ,
    drop = FALSE
  ]

  df
}

select_final_volcano_labels <- function(df, y_col, n_labels = n_top_labels_volcano) {
  if (!nrow(df)) {
    return(df[0, , drop = FALSE])
  }

  # Label only final significant sites. Weak-CNH-only passes are deliberately
  # excluded because they require HBFSS overlap to count as final weak evidence.
  lab_df <- df[
    df$has_valid_gene_symbol &
      !is.na(df$final_significant_flag) &
      df$final_significant_flag,
    ,
    drop = FALSE
  ]

  if (!nrow(lab_df)) {
    return(lab_df[0, , drop = FALSE])
  }

  emp <- suppressWarnings(as.numeric(lab_df$empirical_p))
  emp[!is.finite(emp)] <- Inf

  hscore <- suppressWarnings(as.numeric(lab_df$HBFSS))
  hscore[!is.finite(hscore)] <- -Inf

  # Rank the plotted discoveries by empirical evidence first, then HBFSS and
  # absolute apeglm effect size. Labels are site-based: up to 20 PASs are labeled
  # even when multiple PASs map to the same gene.
  ord <- order(
    emp,
    -hscore,
    -abs(lab_df$lfc_shrunk),
    na.last = TRUE
  )

  lab_df <- lab_df[ord, , drop = FALSE]

  lab_df[seq_len(min(as.integer(n_labels), nrow(lab_df))), , drop = FALSE]
}


make_hbfss_boundary_df <- function(plot_df, hbfss_threshold, y_limit) {
  if (!is.finite(hbfss_threshold) || is.na(hbfss_threshold) || hbfss_threshold < 0) {
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

  # HBFSS significance is score-threshold driven. HCp is shown separately as a
  # reference line because it is used to calculate Htau, not as a second gate.
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

plot_final_volcano <- function(df, dataset_name, short_title = NULL, label_genes = TRUE) {
  plot_df <- build_final_volcano_df(
    df,
    y_col = "neglog10_empirical_p"
  )

  if (!nrow(plot_df)) {
    stop("No finite volcano plotting rows for ", dataset_name)
  }

  lab_df <- if (isTRUE(label_genes)) {
    select_final_volcano_labels(
      plot_df,
      y_col = "neglog10_empirical_p"
    )
  } else {
    plot_df[0, , drop = FALSE]
  }

  hc_raw <- suppressWarnings(
    as.numeric(df$hc_p_threshold_dataset[1])
  )

  hbfss_raw <- suppressWarnings(
    as.numeric(df$hbfss_threshold_dataset[1])
  )

  hc_y <- if (is.finite(hc_raw) && !is.na(hc_raw) && hc_raw > 0 && hc_raw <= 1) {
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

  weak_n <- sum(df$weak_significant_flag, na.rm = TRUE)
  strong_n <- sum(df$strong_cnh_flag, na.rm = TRUE)
  std_n <- sum(df$standard_flag, na.rm = TRUE)
  hbfss_total_n <- sum(df$hbfss_flag, na.rm = TRUE)
  overlap_n <- sum(df$any_overlap, na.rm = TRUE)

  hc_label <- if (is.finite(hc_raw) && !is.na(hc_raw)) {
    paste0("  HCp=", signif(hc_raw, 3))
  } else {
    ""
  }

  hbfss_label <- if (is.finite(hbfss_raw) && !is.na(hbfss_raw)) {
    paste0("  Hτ=", signif(hbfss_raw, 3))
  } else {
    ""
  }

  count_text <- paste0(
    "HBFSS=", hbfss_total_n,
    "  Std=", std_n,
    "  Str=", strong_n,
    "  Wk=", weak_n,
    "  Ovlp=", overlap_n,
    hc_label,
    hbfss_label
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
      shape = final_class,
      size = final_class,
      alpha = final_class
    )
  ) +
    geom_point(
      stroke = 0.30
    ) +
    scale_color_manual(
      values = final_class_colors,
      breaks = final_class_levels,
      drop = FALSE,
      name = "Class",
      guide = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          shape = unname(final_class_shapes[final_class_levels]),
          size = rep(3.2, length(final_class_levels)),
          alpha = rep(1.0, length(final_class_levels)),
          stroke = rep(0.55, length(final_class_levels))
        )
      )
    ) +
    scale_shape_manual(
      values = final_class_shapes,
      breaks = final_class_levels,
      drop = FALSE,
      name = "Class",
      guide = "none"
    ) +
    scale_size_manual(
      values = final_class_sizes,
      breaks = final_class_levels,
      guide = "none"
    ) +
    scale_alpha_manual(
      values = final_class_alphas,
      breaks = final_class_levels,
      guide = "none"
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
      x = "Shrunken log2FC",
      y = expression(-log[10]("Empirical p")),
      caption = compact_caption(count_text, width = 84)
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    plot_expand_xy() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      plot.caption = element_text(
        size = base_theme_size - 2,
        hjust = 0.5,
        margin = margin(t = 4)
      ),
      plot.caption.position = "plot",
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
  }

  x_rng <- range(plot_df$lfc_shrunk, na.rm = TRUE)
  y_rng <- range(plot_df$neglog10_empirical_p, na.rm = TRUE)

  if (nrow(lab_df) > 0L) {
    p <- p +
      ggrepel::geom_text_repel(
        data = lab_df,
        aes(
          x = lfc_shrunk,
          y = neglog10_empirical_p,
          label = gene_symbol_plot,
          color = final_class
        ),
        inherit.aes = FALSE,
        show.legend = FALSE,
        size = 2.05,
        seed = 1,
        max.overlaps = Inf,
        force = 1.15,
        force_pull = 0.35,
        box.padding = 0.30,
        point.padding = 0.14,
        min.segment.length = 0,
        segment.alpha = 0.60,
        segment.size = 0.22
      )
  }

  p
}

run_one_evs_track <- function(comparison_name, track_key, count_matrix, coldata, annot_df, cmp_dir) {
  track_dir <- file.path(
    cmp_dir,
    unname(track_short[track_key])
  )

  tab_dir <- file.path(cmp_dir, "tables")

  for (d in c(track_dir, tab_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  evs <- build_eigenvector_split(
    count_matrix = count_matrix,
    coldata = coldata,
    track_key = track_key
  )

  message(
    "EVS track ", unname(track_short[track_key]),
    ": split IDs selected using ", evs$preprocessing_label,
    "; fixed top-5,000 is selected independently from each condition PC1 eigenvector; downstream EVS Lead/Rem DESeq2 inputs are raw count subsets and are normalized by DESeq2 after the split; Original is analyzed once without EVS."
  )

  # EVS diagnostics are registered here and exported later as manuscript support
  # panels: PCA structure, PC1-loading rank curves, and PC1-loading histograms
  # showing the fixed top-5,000 boundary and the joint/disjoint union split.

  dataset_list <- list(
    raw_dataset = evs$raw_dataset,
    leading_edge_dataset = evs$leading_edge_dataset,
    remainder_dataset = evs$remainder_dataset
  )

  paper_registry[[registry_key(comparison_name, track_key)]] <<- list(
    comparison_name = comparison_name,
    track_key = track_key,
    evs = evs,
    dataset_list = dataset_list,
    coldata = coldata
  )

  analysis_results <- list()

  # The unsplit Original dataset is track-independent and is analyzed once under
  # NormEVS registry bookkeeping. RawEVS therefore contributes only Lead and Rem.
  analysis_dataset_list <- dataset_list
  if (identical(track_key, "raw_evs")) {
    analysis_dataset_list$raw_dataset <- NULL
  }

  for (nm in names(analysis_dataset_list)) {
    full_dataset_name <- paste(
      comparison_name,
      unname(track_short[track_key]),
      nm,
      sep = "_"
    )

    fit <- tryCatch(
      run_core_analysis(
        count_mat = analysis_dataset_list[[nm]],
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

    summary_row <- data.frame(
      comparison_name = comparison_name,
      track_key = track_key,
      track_label = unname(track_short[track_key]),
      dataset_key = nm,
      dataset_name = full_dataset_name,
      evs_split_input = evs$preprocessing_label,
      deseq2_input = evs$downstream_deseq2_input,
      n_features = nrow(df),
      HCp = fit$hc_p_threshold,
      Htau = fit$hbfss_threshold,
      Std_BH_FDR = standard_alpha_level,
      CNH_BH_FDR = alpha_level,
      n_HBFSS = sum(df$hbfss_flag, na.rm = TRUE),
      n_Std = sum(df$standard_flag, na.rm = TRUE),
      n_Strong = sum(df$strong_cnh_flag, na.rm = TRUE),
      n_Weak = sum(df$weak_significant_flag, na.rm = TRUE),
      n_Overlap = sum(df$any_overlap, na.rm = TRUE),
      n_HBFSS_Weak = sum(df$weak_hbfss_overlap, na.rm = TRUE),
      n_HBFSS_Strong = sum(df$strong_hbfss_overlap, na.rm = TRUE),
      n_HBFSS_Std = sum(df$standard_hbfss_overlap, na.rm = TRUE),
      n_Weak_CNH_pass = sum(df$weak_cnh_flag, na.rm = TRUE),
      top_n_per_condition = top_n_target,
      joint_top5000_n = length(evs$joint_ids),
      disjoint_trt_top5000_n = length(evs$disjoint_trt_ids),
      disjoint_untrt_top5000_n = length(evs$disjoint_untrt_ids),
      leading_edge_n = length(evs$leading_edge_ids),
      remainder_n = length(evs$remainder_ids),
      stringsAsFactors = FALSE
    )

    # Final manuscript figures are generated once at the paper-panel stage
    # after all comparisons have completed. No single-dataset diagnostic plots
    # are generated here.

    analysis_results[[nm]] <- list(
      dds = fit$dds,
      results = df,
      summary = summary_row,
      dataset_mat = analysis_dataset_list[[nm]]
    )
  }

  registry_obj <- paper_registry[[registry_key(comparison_name, track_key)]]
  registry_obj$analysis_results <- analysis_results
  paper_registry[[registry_key(comparison_name, track_key)]] <<- registry_obj

  # Cross-dataset diagnostic panels are intentionally omitted from final
  # manuscript export; final paper panels are built at the end of the run.

  if (run_twas_overlap) {
    twas_genes <- clean_gene_set(TWAS_data$gene_symbol)

    get_twas_overlap <- function(result_df, dataset_nm, cmp_name, out_dir, mode = c("union", "standard_only", "hbfss_only")) {
      mode <- match.arg(mode)

      sig_df <- switch(
        mode,
        union = subset(result_df, final_significant_flag),
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

analysis_view_label <- function(track_key, dataset_key) {
  if (identical(dataset_key, "raw_dataset")) {
    return("Original (No EVS)")
  }

  paste0(
    unname(track_short[track_key]),
    " ",
    ifelse(dataset_key == "leading_edge_dataset", "Lead", "Rem")
  )
}

get_registered_result <- function(comparison_name, track_key, dataset_key) {
  obj <- paper_registry[[registry_key(comparison_name, track_key)]]

  if (is.null(obj) || is.null(obj$analysis_results) ||
      is.null(obj$analysis_results[[dataset_key]])) {
    return(NULL)
  }

  obj$analysis_results[[dataset_key]]$results
}

build_support_label <- function(df) {
  if (!nrow(df)) {
    return(character(0))
  }

  vapply(seq_len(nrow(df)), function(i) {
    tags <- character(0)

    if (isTRUE(df$hbfss_flag[i])) tags <- c(tags, "HBFSS")
    if (isTRUE(df$standard_flag[i])) tags <- c(tags, "Std")
    if (isTRUE(df$strong_cnh_flag[i])) tags <- c(tags, "Strong")
    if (isTRUE(df$weak_significant_flag[i])) tags <- c(tags, "Weak")

    if (!length(tags)) "" else paste(tags, collapse = "+")
  }, character(1))
}

comparison_analysis_views <- function() {
  data.frame(
    track_key = c(
      "normalized_evs",
      "normalized_evs", "normalized_evs",
      "raw_evs", "raw_evs"
    ),
    dataset_key = c(
      "raw_dataset",
      "leading_edge_dataset", "remainder_dataset",
      "leading_edge_dataset", "remainder_dataset"
    ),
    stringsAsFactors = FALSE
  )
}

export_comparison_manuscript_tables <- function(comparison_name) {
  cmp_dir <- file.path(output_dir, comparison_name)
  tab_dir <- file.path(cmp_dir, "tables")
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

  views <- comparison_analysis_views()
  sig_rows <- list()
  summary_rows <- list()

  for (i in seq_len(nrow(views))) {
    track_key <- views$track_key[i]
    dataset_key <- views$dataset_key[i]
    df <- get_registered_result(comparison_name, track_key, dataset_key)

    if (is.null(df) || !nrow(df)) {
      next
    }

    analysis_label <- analysis_view_label(track_key, dataset_key)
    hc <- suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1]))
    htau <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      Analysis = analysis_label,
      N = nrow(df),
      HCp = hc,
      Htau = htau,
      HBFSS = sum(df$hbfss_flag, na.rm = TRUE),
      Std = sum(df$standard_flag, na.rm = TRUE),
      Strong = sum(df$strong_cnh_flag, na.rm = TRUE),
      Weak = sum(df$weak_significant_flag, na.rm = TRUE),
      Ovlp = sum(df$any_overlap, na.rm = TRUE),
      BH_FDR = standard_alpha_level,
      stringsAsFactors = FALSE
    )

    sig <- df[
      !is.na(df$final_significant_flag) & df$final_significant_flag,
      ,
      drop = FALSE
    ]

    if (!nrow(sig)) {
      next
    }

    pas <- if ("orig_id" %in% names(sig)) as.character(sig$orig_id) else as.character(sig$feature_id)
    bad_pas <- is.na(pas) | !nzchar(trimws(pas))
    pas[bad_pas] <- as.character(sig$feature_id[bad_pas])

    gene <- if ("gene_symbol" %in% names(sig)) as.character(sig$gene_symbol) else rep(NA_character_, nrow(sig))
    gene[is.na(gene) | !nzchar(trimws(gene))] <- NA_character_

    sig_out <- data.frame(
      Analysis = analysis_label,
      PAS = pas,
      Gene = gene,
      Direction = as.character(sig$regulation_direction),
      Apeglm_LFC = as.numeric(sig$lfc_shrunk),
      EmpP = as.numeric(sig$empirical_p),
      HBFSS = as.numeric(sig$HBFSS),
      HCp = hc,
      Htau = htau,
      Std_BH = as.numeric(sig$padj),
      Strong_BH = as.numeric(sig$resGA_padj),
      Weak_BH = as.numeric(sig$resLA_padj),
      Support = build_support_label(sig),
      stringsAsFactors = FALSE
    )

    sig_out <- sig_out[
      order(
        sig_out$Analysis,
        sig_out$EmpP,
        -sig_out$HBFSS,
        -abs(sig_out$Apeglm_LFC),
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]

    sig_rows[[length(sig_rows) + 1L]] <- sig_out
  }

  sig_table <- if (length(sig_rows)) dplyr::bind_rows(sig_rows) else data.frame()
  summary_table <- if (length(summary_rows)) dplyr::bind_rows(summary_rows) else data.frame()

  save_csv(
    sig_table,
    file.path(tab_dir, paste0("Table_", comparison_name, "_Significant_Sites.csv"))
  )

  save_csv(
    summary_table,
    file.path(tab_dir, paste0("Table_", comparison_name, "_Summary.csv"))
  )

  invisible(list(significant = sig_table, summary = summary_table))
}

build_overall_manuscript_summary <- function() {
  rows <- list()
  views <- comparison_analysis_views()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (i in seq_len(nrow(views))) {
      track_key <- views$track_key[i]
      dataset_key <- views$dataset_key[i]
      df <- get_registered_result(comparison_name, track_key, dataset_key)

      if (is.null(df) || !nrow(df)) {
        next
      }

      rows[[length(rows) + 1L]] <- data.frame(
        Comparison = comparison_name,
        Analysis = analysis_view_label(track_key, dataset_key),
        N = nrow(df),
        HCp = suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1])),
        Htau = suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1])),
        HBFSS = sum(df$hbfss_flag, na.rm = TRUE),
        Std = sum(df$standard_flag, na.rm = TRUE),
        Strong = sum(df$strong_cnh_flag, na.rm = TRUE),
        Weak = sum(df$weak_significant_flag, na.rm = TRUE),
        Ovlp = sum(df$any_overlap, na.rm = TRUE),
        BH_FDR = standard_alpha_level,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) data.frame() else dplyr::bind_rows(rows)
}

write_methods_note <- function() {
  methods_text <- c(
    "Statistical decision rules",
    "",
    "Eigenvector splitting followed the supplied Figure 1 / Figure 1B workflow. For NormEVS, the full comparison raw-count matrix was normalized by DESeq2 median-of-ratios size factors before EVS; RawEVS omitted this pre-split normalization.",
    sprintf("Within each EVS track, treatment and control samples were decomposed separately by prcomp; absolute PC1 feature loadings were ranked and exactly the top %d PASs from each condition were selected. No pre-PCA log2 transformation and no adaptive reduction of top-N were used.", top_n_target),
    "Top-set intersection defined Joint passes; condition-specific set differences defined Disjoint passes; Joint plus both Disjoint sets formed the Leading Edge, and every other PAS formed the Remainder.",
    "PAS features were modeled with DESeq2 using design ~ condition. Log2 fold changes used for effect interpretation and HBFSS were shrunken with apeglm.",
    sprintf("Standard DESeq2 significance required BH-adjusted padj < %.2f and |apeglm LFC| >= %.1f.", standard_alpha_level, lfc_boundary),
    sprintf("Strong CNH used DESeq2 greaterAbs with lfcThreshold = %.1f and required BH-adjusted padj < %.2f plus |apeglm LFC| >= %.1f.", lfc_boundary, alpha_level, lfc_boundary),
    sprintf("Weak CNH used DESeq2 lessAbs with lfcThreshold = %.1f and BH-adjusted padj < %.2f plus |apeglm LFC| < %.1f; a weak-effect PAS was counted as a final significant weak discovery only when it also passed HBFSS.", lfc_boundary, alpha_level, lfc_boundary),
    "Finite DESeq2 Wald statistics were calibrated with the Strimmer fdrtool empirical-null model; its empirical-null p-values were used in HBFSS.",
    sprintf("For each dataset, HCp was obtained with hc.thresh(sorted empirical p-values), Htau = |-log10(HCp)| x %.1f, and HBFSS = |apeglm LFC| x [-log10(empirical p)].", lfc_boundary),
    "HBFSS significance required HBFSS > Htau. No DESeq2 padj/q-value, no separate empirical-p <= HCp gate, and no additional upper validity cutoff on HCp were imposed; HCp was used directly to derive Htau whenever hc.thresh returned a valid probability in (0, 1].",
    "HBFSS/DESeq2 overlap was defined as HBFSS plus any Standard, Strong CNH, or Weak CNH support. On volcano plots all HBFSS-positive sites are purple; marker shape identifies Standard, Strong, or Weak overlap. Weak-CNH-only passes are blue but are excluded from final significant-site tables.",
    sprintf("Each volcano labels at most the top %d final significant PASs, using gene symbol when available and PAS ID otherwise.", n_top_labels_volcano),
    "The five manuscript analysis views per comparison are Original (No EVS), NormEVS Lead, NormEVS Rem, RawEVS Lead, and RawEVS Rem. EVS preprocessing affects only PC1-based feature membership; downstream DESeq2 always uses the corresponding raw-count subset and performs its own size-factor normalization."
  )

  writeLines(
    methods_text,
    con = file.path(output_dir, "Methods_Statistical_Decision_Rules.txt"),
    useBytes = TRUE
  )

  invisible(TRUE)
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

  export_comparison_manuscript_tables(comparison_name)

  dplyr::bind_rows(track_summaries)
}


# -----------------------------------------------------------------------------
# Paper-ready multi-comparison volcano panels
# -----------------------------------------------------------------------------

read_result_table_for_panel <- function(comparison_name, track_key, dataset_key) {
  df <- get_registered_result(
    comparison_name = comparison_name,
    track_key = track_key,
    dataset_key = dataset_key
  )

  if (is.null(df)) {
    warning(
      "Missing registered result for manuscript panel: ",
      comparison_name, " / ", track_key, " / ", dataset_key
    )
    return(NULL)
  }

  df$final_class <- factor(
    as.character(df$final_class),
    levels = final_class_levels
  )

  df
}


save_paper_volcano_panels <- function() {
  # Original data are not affected by the EVS preprocessing path, so only one
  # original-data volcano is exported per comparison. Leading-edge and remainder
  # panels retain both NormEVS and RawEVS rows because the split membership can
  # differ between those two EVS inputs.
  dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)

  comparison_order <- as.character(comparison_table$comparison_name)

  for (dataset_key in dataset_key_order) {
    plots <- list()

    if (identical(dataset_key, "raw_dataset")) {
      track_order_use <- "normalized_evs"
      panel_title <- "Original dataset volcanoes across all comparisons"
      output_suffix <- "AllComparisons"
      panel_ncol <- length(comparison_order)
      panel_width <- 18.0
      panel_height <- 6.6
    } else {
      track_order_use <- c("normalized_evs", "raw_evs")
      panel_title <- paste(
        pretty_dataset_type(dataset_key),
        "volcano comparison: NormEVS vs RawEVS across all comparisons"
      )
      output_suffix <- "AllComparisons_NormEVS_vs_RawEVS"
      panel_ncol <- length(comparison_order)
      panel_width <- 18.0
      panel_height <- 11.0
    }

    for (track_key in track_order_use) {
      for (comparison_name in comparison_order) {
        df <- read_result_table_for_panel(
          comparison_name = comparison_name,
          track_key = track_key,
          dataset_key = dataset_key
        )

        if (is.null(df)) {
          next
        }

        dataset_name <- paste(
          comparison_name,
          unname(track_short[track_key]),
          dataset_key,
          sep = "_"
        )

        short_title <- if (identical(dataset_key, "raw_dataset")) {
          comparison_name
        } else {
          paste0(comparison_name, "\n", unname(track_short[track_key]))
        }

        plots[[length(plots) + 1L]] <- plot_final_volcano(
          df = df,
          dataset_name = dataset_name,
          short_title = short_title,
          label_genes = TRUE
        ) +
          theme(
            plot.title = element_text(size = base_theme_size, face = "bold"),
            axis.title = element_text(size = base_theme_size - 1),
            axis.text = element_text(size = base_theme_size - 2)
          )
      }
    }

    if (length(plots) == 0L) {
      next
    }

    panel <- assemble_one_legend_panel(
      plots,
      panel_title = panel_title,
      ncol = panel_ncol
    )

    save_grob(
      panel,
      file.path(
        paper_fig_dir,
        paste0(
          "Figure_Manuscript_Volcano_",
          unname(dataset_short[dataset_key]),
          "_",
          output_suffix,
          ".png"
        )
      ),
      width = panel_width,
      height = panel_height
    )
  }

  invisible(TRUE)
}


make_single_legend_scale <- function() {
  list(
    scale_color_manual(
      values = final_class_colors,
      breaks = final_class_levels,
      drop = FALSE,
      name = "Class",
      guide = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          shape = unname(final_class_shapes[final_class_levels]),
          size = rep(3.0, length(final_class_levels)),
          alpha = rep(1.0, length(final_class_levels)),
          stroke = rep(0.55, length(final_class_levels))
        )
      )
    ),
    scale_shape_manual(
      values = final_class_shapes,
      breaks = final_class_levels,
      drop = FALSE,
      name = "Class",
      guide = "none"
    )
  )
}

compute_pca_support_plot <- function(count_df, coldata, short_title) {
  count_mat <- coerce_raw_count_matrix_for_deseq2(
    count_df,
    context = paste0(short_title, " PCA count matrix")
  )

  common_samples <- intersect(colnames(count_mat), rownames(coldata))
  count_mat <- count_mat[, common_samples, drop = FALSE]
  coldata <- coldata[common_samples, , drop = FALSE]

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = make_design_formula(coldata)
  )

  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]

  if (nrow(dds) < 2L || ncol(dds) < 3L) {
    return(NULL)
  }

  dds <- DESeq2::estimateSizeFactors(dds)
  x <- log2(DESeq2::counts(dds, normalized = TRUE) + 1)

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
    Condition = factor(
      as.character(coldata[rownames(pca_fit$x), "condition"]),
      levels = c("untrt", "trt")
    ),
    stringsAsFactors = FALSE
  )

  ggplot(
    pca_df,
    aes(PC1, PC2, label = Sample, shape = Condition, fill = Condition)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey70") +
    geom_point(size = 2.6, colour = "white", stroke = 0.5) +
    ggrepel::geom_text_repel(
      size = 1.7,
      max.overlaps = 10,
      force = 0.9,
      box.padding = 0.18,
      point.padding = 0.10,
      min.segment.length = 0,
      segment.alpha = 0.45,
      segment.size = 0.16
    ) +
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition") +
    scale_fill_manual(values = condition_fills, labels = condition_labels, name = "Condition") +
    labs(
      title = short_title,
      x = paste0("PC1 ", pca_var_per[1], "%"),
      y = paste0("PC2 ", pca_var_per[2], "%"),
      caption = "DESeq2-normalized after dataset definition"
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(
      legend.position = "bottom",
      plot.caption = element_text(size = base_theme_size - 3, hjust = 0.5),
      plot.margin = margin(8, 8, 8, 8)
    )
}

plot_evs_rank_support <- function(evs, comparison_name, track_key) {
  trt_df <- evs$fit_trt$loading_table
  ctrl_df <- evs$fit_untrt$loading_table

  trt_df$Group <- "Treatment"
  ctrl_df$Group <- "Control"

  rank_df <- rbind(
    trt_df[, c("rank", "pc1_loading_abs", "Group")],
    ctrl_df[, c("rank", "pc1_loading_abs", "Group")]
  )

  rank_df <- rank_df[
    is.finite(rank_df$pc1_loading_abs) &
      !is.na(rank_df$pc1_loading_abs) &
      rank_df$pc1_loading_abs > 0,
    ,
    drop = FALSE
  ]

  cutoff_df <- data.frame(
    Group = c("Treatment", "Control"),
    cutoff = c(evs$fit_trt$cutoff, evs$fit_untrt$cutoff),
    stringsAsFactors = FALSE
  )

  ggplot(
    rank_df,
    aes(rank, pc1_loading_abs, color = Group)
  ) +
    geom_line(linewidth = 0.45, alpha = 0.95) +
    geom_hline(
      data = cutoff_df,
      aes(yintercept = cutoff, color = Group),
      linetype = "dashed",
      linewidth = 0.55,
      inherit.aes = FALSE
    ) +
    scale_x_log10(labels = scales::label_number()) +
    scale_y_log10(labels = scales::label_number()) +
    scale_color_manual(
      values = c("Treatment" = plot_palette$treatment, "Control" = plot_palette$control),
      name = "Group"
    ) +
    labs(
      title = paste0(comparison_name, "\n", unname(track_short[track_key])),
      x = "EVS rank",
      y = "|PC1 loading|",
      caption = paste0(
        "Top 5,000/condition; Joint=", length(evs$joint_ids),
        "  Disj-T=", length(evs$disjoint_trt_ids),
        "  Disj-C=", length(evs$disjoint_untrt_ids),
        "  Lead=", length(evs$leading_edge_ids),
        "  Rem=", length(evs$remainder_ids)
      )
    ) +
    manuscript_theme() +
    theme(
      legend.position = "bottom",
      plot.caption = element_text(size = base_theme_size - 3, hjust = 0.5),
      plot.margin = margin(8, 8, 8, 8)
    )
}

plot_evs_loading_histogram_support <- function(evs, comparison_name, track_key) {
  trt_df <- evs$fit_trt$loading_table
  ctrl_df <- evs$fit_untrt$loading_table

  trt_df$Group <- "Treatment"
  ctrl_df$Group <- "Control"

  hist_df <- rbind(
    trt_df[, c("pc1_loading_abs", "split_class", "Group")],
    ctrl_df[, c("pc1_loading_abs", "split_class", "Group")]
  )

  hist_df <- hist_df[
    is.finite(hist_df$pc1_loading_abs) &
      !is.na(hist_df$pc1_loading_abs) &
      hist_df$pc1_loading_abs > 0,
    ,
    drop = FALSE
  ]

  cutoff_df <- data.frame(
    Group = c("Treatment", "Control"),
    cutoff = c(evs$fit_trt$cutoff, evs$fit_untrt$cutoff),
    stringsAsFactors = FALSE
  )

  ggplot(
    hist_df,
    aes(pc1_loading_abs, fill = Group, color = Group)
  ) +
    geom_histogram(
      bins = 60,
      alpha = 0.28,
      position = "identity",
      linewidth = 0.15
    ) +
    geom_vline(
      data = cutoff_df,
      aes(xintercept = cutoff, color = Group),
      linetype = "dashed",
      linewidth = 0.60,
      inherit.aes = FALSE
    ) +
    scale_x_log10(labels = scales::label_number()) +
    scale_fill_manual(
      values = c("Treatment" = plot_palette$treatment, "Control" = plot_palette$control),
      name = "Group"
    ) +
    scale_color_manual(
      values = c("Treatment" = plot_palette$treatment, "Control" = plot_palette$control),
      name = "Group"
    ) +
    labs(
      title = paste0(comparison_name, "\n", unname(track_short[track_key])),
      x = "|PC1 loading|",
      y = "Feature count",
      caption = paste0(
        "Fixed top 5,000/condition; Lead = Joint + Disjoint treatment + Disjoint control (n=",
        length(evs$leading_edge_ids), ")"
      )
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(
      legend.position = "bottom",
      plot.caption = element_text(size = base_theme_size - 3, hjust = 0.5),
      plot.margin = margin(8, 8, 8, 8)
    )
}

plot_empirical_hbfss_support <- function(df, short_title) {
  req <- c("empirical_p", "HBFSS", "final_class", "hc_p_threshold_dataset", "hbfss_threshold_dataset")
  if (length(setdiff(req, names(df))) > 0L) {
    return(NULL)
  }

  plot_df <- df[
    is.finite(df$empirical_p) &
      !is.na(df$empirical_p) &
      is.finite(df$HBFSS) &
      !is.na(df$HBFSS),
    ,
    drop = FALSE
  ]

  if (nrow(plot_df) == 0L) {
    return(NULL)
  }

  plot_df$empirical_p <- pmax(plot_df$empirical_p, 1e-300)
  plot_df$final_class <- factor(as.character(plot_df$final_class), levels = final_class_levels)

  hc_raw <- suppressWarnings(as.numeric(plot_df$hc_p_threshold_dataset[1]))
  hbfss_raw <- suppressWarnings(as.numeric(plot_df$hbfss_threshold_dataset[1]))

  caption_text <- paste0(
    if (is.finite(hc_raw) && !is.na(hc_raw)) paste0("HCp=", signif(hc_raw, 3)) else "HCp=NA",
    "  ",
    if (is.finite(hbfss_raw) && !is.na(hbfss_raw)) paste0("Hτ=", signif(hbfss_raw, 3)) else "Hτ=NA"
  )

  p <- ggplot(
    plot_df,
    aes(empirical_p, HBFSS, color = final_class, shape = final_class)
  ) +
    geom_point(alpha = 0.65, size = 1.05, stroke = 0.25) +
    make_single_legend_scale() +
    scale_x_log10(labels = scales::label_scientific()) +
    labs(
      title = short_title,
      x = "Empirical p",
      y = "HBFSS",
      caption = caption_text
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(
      legend.position = "bottom",
      plot.caption = element_text(size = base_theme_size - 3, hjust = 0.5),
      plot.margin = margin(8, 8, 8, 8)
    )

  if (is.finite(hc_raw) && !is.na(hc_raw) && hc_raw > 0 && hc_raw <= 1) {
    p <- p +
      geom_vline(
        xintercept = hc_raw,
        color = plot_palette$hc,
        linewidth = 0.55,
        linetype = "dotted"
      )
  }

  if (is.finite(hbfss_raw) && !is.na(hbfss_raw)) {
    p <- p +
      geom_hline(
        yintercept = hbfss_raw,
        color = plot_palette$hbfss_line,
        linewidth = 0.55,
        linetype = "dashed"
      )
  }

  p
}

save_paper_pca_panels <- function() {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) {
    return(invisible(FALSE))
  }

  comparison_order <- as.character(comparison_table$comparison_name)

  for (dataset_key in dataset_key_order) {
    plots <- list()

    if (identical(dataset_key, "raw_dataset")) {
      track_order_use <- "normalized_evs"
      panel_title <- "PCA structure: original dataset"
      output_suffix <- "Raw_AllComparisons"
      panel_height <- 5.8
    } else {
      track_order_use <- c("normalized_evs", "raw_evs")
      panel_title <- paste0("PCA structure: ", pretty_dataset_type(dataset_key), " NormEVS vs RawEVS")
      output_suffix <- paste0(unname(dataset_short[dataset_key]), "_AllComparisons_NormEVS_vs_RawEVS")
      panel_height <- 11.0
    }

    for (track_key in track_order_use) {
      for (comparison_name in comparison_order) {
        obj <- paper_registry[[registry_key(comparison_name, track_key)]]

        if (is.null(obj) || is.null(obj$dataset_list[[dataset_key]])) {
          next
        }

        short_title <- if (identical(dataset_key, "raw_dataset")) {
          comparison_name
        } else {
          paste0(comparison_name, "\n", unname(track_short[track_key]))
        }

        plots[[length(plots) + 1L]] <- compute_pca_support_plot(
          count_df = obj$dataset_list[[dataset_key]],
          coldata = obj$coldata,
          short_title = short_title
        )
      }
    }

    plots <- Filter(Negate(is.null), plots)

    if (length(plots) == 0L) {
      next
    }

    panel <- assemble_one_legend_panel(
      plots,
      panel_title = panel_title,
      ncol = length(comparison_order)
    )

    save_grob(
      panel,
      file.path(
        paper_fig_dir,
        paste0("Figure_Manuscript_PCA_", output_suffix, ".png")
      ),
      width = 18.0,
      height = panel_height
    )
  }

  invisible(TRUE)
}

save_paper_evs_rank_panel <- function() {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) {
    return(invisible(FALSE))
  }

  plots <- list()

  for (track_key in c("normalized_evs", "raw_evs")) {
    for (comparison_name in as.character(comparison_table$comparison_name)) {
      obj <- paper_registry[[registry_key(comparison_name, track_key)]]

      if (is.null(obj) || is.null(obj$evs)) {
        next
      }

      plots[[length(plots) + 1L]] <- plot_evs_rank_support(
        evs = obj$evs,
        comparison_name = comparison_name,
        track_key = track_key
      )
    }
  }

  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 0L) {
    return(invisible(FALSE))
  }

  panel <- assemble_one_legend_panel(
    plots,
    panel_title = "EVS PC1-loading rank support: NormEVS vs RawEVS",
    ncol = length(as.character(comparison_table$comparison_name))
  )

  save_grob(
    panel,
    file.path(paper_fig_dir, "Figure_Manuscript_EVS_PC1_Loading_Rank.png"),
    width = 18.0,
    height = 9.4
  )

  invisible(TRUE)
}

save_paper_evs_loading_histogram_panel <- function() {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) {
    return(invisible(FALSE))
  }

  plots <- list()

  for (track_key in c("normalized_evs", "raw_evs")) {
    for (comparison_name in as.character(comparison_table$comparison_name)) {
      obj <- paper_registry[[registry_key(comparison_name, track_key)]]

      if (is.null(obj) || is.null(obj$evs)) {
        next
      }

      plots[[length(plots) + 1L]] <- plot_evs_loading_histogram_support(
        evs = obj$evs,
        comparison_name = comparison_name,
        track_key = track_key
      )
    }
  }

  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 0L) {
    return(invisible(FALSE))
  }

  panel <- assemble_one_legend_panel(
    plots,
    panel_title = "EVS PC1-loading distribution support: NormEVS vs RawEVS",
    ncol = length(as.character(comparison_table$comparison_name))
  )

  save_grob(
    panel,
    file.path(paper_fig_dir, "Figure_Manuscript_EVS_PC1_Loading_Histogram.png"),
    width = 18.0,
    height = 9.4
  )

  invisible(TRUE)
}

save_paper_empirical_hbfss_panels <- function() {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) {
    return(invisible(FALSE))
  }

  comparison_order <- as.character(comparison_table$comparison_name)

  for (dataset_key in dataset_key_order) {
    plots <- list()

    if (identical(dataset_key, "raw_dataset")) {
      track_order_use <- "normalized_evs"
      panel_title <- "Empirical p and HBFSS support: original dataset"
      output_suffix <- "Raw_AllComparisons"
      panel_height <- 5.8
    } else {
      track_order_use <- c("normalized_evs", "raw_evs")
      panel_title <- paste0("Empirical p and HBFSS support: ", pretty_dataset_type(dataset_key), " NormEVS vs RawEVS")
      output_suffix <- paste0(unname(dataset_short[dataset_key]), "_AllComparisons_NormEVS_vs_RawEVS")
      panel_height <- 11.0
    }

    for (track_key in track_order_use) {
      for (comparison_name in comparison_order) {
        df <- read_result_table_for_panel(
          comparison_name = comparison_name,
          track_key = track_key,
          dataset_key = dataset_key
        )

        if (is.null(df)) {
          next
        }

        short_title <- if (identical(dataset_key, "raw_dataset")) {
          comparison_name
        } else {
          paste0(comparison_name, "\n", unname(track_short[track_key]))
        }

        plots[[length(plots) + 1L]] <- plot_empirical_hbfss_support(
          df = df,
          short_title = short_title
        )
      }
    }

    plots <- Filter(Negate(is.null), plots)

    if (length(plots) == 0L) {
      next
    }

    panel <- assemble_one_legend_panel(
      plots,
      panel_title = panel_title,
      ncol = length(comparison_order)
    )

    save_grob(
      panel,
      file.path(
        paper_fig_dir,
        paste0("Figure_Manuscript_Empirical_HBFSS_", output_suffix, ".png")
      ),
      width = 18.0,
      height = panel_height
    )
  }

  invisible(TRUE)
}

build_discovery_long_table <- function(summary_df) {
  if (!is.data.frame(summary_df) || nrow(summary_df) == 0L) {
    return(data.frame())
  }

  rows <- lapply(seq_len(nrow(summary_df)), function(i) {
    sm <- summary_df[i, , drop = FALSE]

    data.frame(
      Comparison = sm$Comparison,
      Analysis = sm$Analysis,
      Class = factor(
        c("HBFSS", "Std", "Strong", "Weak(H∩)", "Overlap"),
        levels = c("HBFSS", "Std", "Strong", "Weak(H∩)", "Overlap")
      ),
      Count = as.numeric(c(
        sm$HBFSS,
        sm$Std,
        sm$Strong,
        sm$Weak,
        sm$Ovlp
      )),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}


save_discovery_count_panel <- function(summary_df) {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) {
    return(invisible(FALSE))
  }

  long_df <- build_discovery_long_table(summary_df)

  if (nrow(long_df) == 0L) {
    return(invisible(FALSE))
  }

  p <- ggplot(
    long_df,
    aes(Comparison, Count, fill = Class)
  ) +
    geom_col(
      position = position_dodge(width = 0.82),
      width = 0.74,
      color = "grey25",
      linewidth = 0.15
    ) +
    facet_wrap(~ Analysis, scales = "free_y", ncol = 3) +
    scale_fill_manual(
      values = c(
        "HBFSS" = plot_palette$hbfss,
        "Std" = plot_palette$standard,
        "Strong" = plot_palette$strong,
        "Weak(H∩)" = plot_palette$weak,
        "Overlap" = plot_palette$overlap
      ),
      drop = FALSE,
      name = "Method"
    ) +
    labs(
      title = "Discovery counts by dataset, EVS mode, and method",
      x = NULL,
      y = "Significant sites",
      caption = "Method totals are not mutually exclusive. Weak(H∩) counts only Weak CNH passes with HBFSS support; Overlap is HBFSS with any DESeq2 class."
    ) +
    manuscript_theme() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(face = "bold")
    )

  save_grob(
    p,
    file.path(paper_fig_dir, "Figure_Manuscript_Discovery_Counts.png"),
    width = 14.5,
    height = 8.2
  )

  invisible(TRUE)
}

save_paper_support_figures <- function(summary_df) {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) {
    return(invisible(FALSE))
  }

  support_steps <- list(
    PCA = function() save_paper_pca_panels(),
    EVS_Rank = function() save_paper_evs_rank_panel(),
    EVS_Histogram = function() save_paper_evs_loading_histogram_panel(),
    Empirical_HBFSS = function() save_paper_empirical_hbfss_panels(),
    Counts = function() save_discovery_count_panel(summary_df)
  )

  for (nm in names(support_steps)) {
    tryCatch(
      support_steps[[nm]](),
      error = function(e) {
        warning("Support figure generation failed at ", nm, ": ", conditionMessage(e))
      }
    )
  }

  invisible(TRUE)
}

clean_manuscript_table_exports()

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

overall_summary <- build_overall_manuscript_summary()

if (nrow(overall_summary) > 0L) {
  save_csv(
    overall_summary,
    file.path(output_dir, "Table_Overall_Summary.csv")
  )
}

write_methods_note()
clean_manuscript_figure_exports()

tryCatch(
  {
    save_paper_volcano_panels()
    save_paper_support_figures(overall_summary)
  },
  error = function(e) {
    warning("Paper figure generation failed: ", conditionMessage(e))
  }
)

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

if (nrow(overall_summary) > 0L) {
  print(overall_summary)
}
