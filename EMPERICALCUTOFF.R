#!/usr/bin/env Rscript

PIPELINE_BUILD <- "SEQUENCE_JOURNAL_READY_GITHUB_2026-08-22"

# =============================================================================
# SEQUENCE MANUSCRIPT ANALYSIS
# WTTS-Seq PAS analysis with eigenvector splitting, DESeq2, apeglm,
# empirical-null calibration, higher criticism, HBFSS, and 3'aTWAS overlap
# =============================================================================
#
# ANALYSIS UNIT
# Each OrigID is analyzed as an individual polyadenylation-site (PAS) feature.
# Gene symbols are retained as annotation and are not used to collapse PASs
# before differential-expression testing.
#
# COMPARISONS
#   RT0 vs ZT6
#   RT2 vs ZT8
#   RT4 vs ZT10
#   RT8 vs ZT14
#
# ANALYSIS VIEWS PER COMPARISON
#   Original (No EVS)
#   NormEVS Lead
#   NormEVS Rem
#   RawEVS Lead
#   RawEVS Rem
#
# EIGENVECTOR SPLITTING
# NormEVS uses DESeq2 median-of-ratios normalized counts before PCA. RawEVS
# uses raw counts before PCA. Within each condition, prcomp is applied directly
# to the corresponding feature-by-sample matrix, absolute PC1 feature loadings
# are ranked, and the 5,000 highest-loading PASs are selected. PASs present in
# both condition-specific top-5,000 sets are Joint; PASs present in only one set
# are Disjoint. Joint plus both Disjoint sets form the Leading Edge. All other
# PASs form the Remainder. Downstream DESeq2 always receives raw counts for the
# selected PAS subset and estimates its own size factors and dispersions.
#
# DIFFERENTIAL EXPRESSION AND EFFECT TESTS
# DESeq2 uses design ~ condition with trt relative to untrt. The ordinary Wald
# p-value is adjusted by Benjamini-Hochberg at FDR 10%. The manuscript Standard
# effect is the ordinary DESeq2 BH-significant result restricted to PASs with
# |apeglm-shrunken LFC| >= 1, so Standard markers cannot occur inside the stated
# effect boundary. Strong (decoupled) calls use DESeq2 greaterAbs with
# lfcThreshold = 1, BH padj < 0.10, and the same |apeglm LFC| >= 1 reporting
# boundary. Weak-CNH support uses DESeq2 lessAbs with lfcThreshold = 1 and BH
# padj < 0.20; the final Weak category additionally requires |apeglm LFC| < 1
# and HBFSS significance. apeglm-shrunken LFC is the reported effect estimate,
# HBFSS effect term, and x-coordinate in all significance figures.
#
# EMPIRICAL NULL, HIGHER CRITICISM, AND HBFSS
# Finite DESeq2 Wald statistics are calibrated with fdrtool using a normal
# empirical-null model. HBFSS uses the resulting empirical p-values. Higher
# criticism is applied to the sorted empirical p-values to obtain HCp. With the
# manuscript LFC boundary c = 1:
#
#   Htau = -log10(HCp) * c
#   HBFSS = |apeglm LFC| * [-log10(empirical p)]
#
# A PAS is HBFSS-significant when HBFSS > Htau. HCp is used to derive Htau and
# is not imposed as an additional significance gate.
#
# VOLCANO FIGURES
# Volcano x-axis: apeglm-shrunken log2 fold change.
# Volcano y-axis: -log10(empirical p) from the fdrtool empirical-null model.
# Standard, Strong, Weak, and HBFSS are plotted as separate method layers using
# one fixed color/marker key across every significance figure. Weak markers are
# shown only for final Weak discoveries (lessAbs + HBFSS), never for lessAbs-only
# PASs. Method overlap is reported numerically and by superimposed method markers;
# it is not treated as a fifth significance method. Each volcano labels at most
# the top 20 final significant PASs.
#
# 3'aTWAS COMPARISON
# After all WTTS analyses and manuscript figures are complete, human 3'aTWAS
# gene symbols are mapped to rat orthologs. TWAS ortholog overlap is evaluated
# against every final WTTS significance method (Standard, Strong, Weak, HBFSS)
# in every analysis view. The TWAS exports identify the human TWAS symbol, rat
# ortholog, significant PASs, method support, and whether the rat gene contains
# multiple WTTS PAS features consistent with alternative polyadenylation.
# =============================================================================

required_packages <- c(
  "DESeq2",
  "apeglm",
  "fdrtool",
  "ggplot2",
  "ggrepel",
  "dplyr",
  "tidyr",
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
  library(tidyr)
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
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  "WTTS-Seq_2022.2_DE_raw_read_numbers(20260822-183312).csv",
  "/mnt/data/WTTS-Seq_2022.2_DE_raw_read_numbers(20260822-183312).csv"
)

twas_file_candidates <- c(
  "3aTWAS_genes_of_11_brain_disorders.csv",
  file.path("data", "3aTWAS_genes_of_11_brain_disorders.csv"),
  "/root/REAPER98632/data/3aTWAS_genes_of_11_brain_disorders.csv",
  "/mnt/data/3aTWAS_genes_of_11_brain_disorders.csv"
)

TWAS_TARGET_SPECIES <- "rat"
TWAS_ORTHOLOG_MIN_SUPPORT <- 1L

# DESeq2 Standard and Strong tests use Benjamini-Hochberg FDR 10%.
# Weak-CNH uses a more permissive BH screen and becomes a final Weak discovery
# only after independent HBFSS support. HBFSS uses no DESeq2 adjusted-p-value gate.
BH_FDR_STANDARD <- 0.10
BH_FDR_STRONG <- 0.10
BH_FDR_WEAK <- 0.20

lfc_boundary <- 1.0


# Published WTTS-seq sleep-pressure/TWAS reference from Flores et al. (2024),
# Table 3. These values are used only as a runtime reproducibility assertion;
# they do not alter any statistical call or create an additional result table.
PUBLISHED_2024_TWAS_DE_REFERENCE <- data.frame(
  PAS = c(
    "467701", "467700", "191750", "218787", "549537", "35839",
    "205547", "205543", "261774", "519450", "78129"
  ),
  Gene = c(
    "Mark3", "Mark3", "Ndfip2", "Wac", "Sf3b1", "Mtrf1l",
    "Gatad2a", "Gatad2a", "Spg7", "Snx19", "Plekhm1"
  ),
  Paper_LFC = c(
    3.18, 3.00, -1.15, -1.51, -1.30, -1.60,
    1.03, 1.42, 2.17, -1.18, 1.96
  ),
  Paper_p = c(
    5.81e-09, 7.13e-08, 2.55e-08, 1.63e-04, 6.47e-04, 9.64e-04,
    9.24e-05, 2.29e-03, 2.75e-04, 3.65e-04, 2.59e-03
  ),
  Paper_padj = c(
    6.94e-06, 4.95e-05, 2.19e-05, 1.30e-02, 3.15e-02, 4.08e-02,
    9.16e-03, 7.24e-02, 1.89e-02, 2.20e-02, 7.74e-02
  ),
  stringsAsFactors = FALSE
)

# EVS method specification: select exactly the 5,000 highest absolute PC1
# loadings independently from each condition eigenvector.
top_n_target <- 5000L

figure_dpi <- 320

base_theme_size <- 10

n_top_labels_volcano <- 20L

# Manuscript export behavior. When TRUE, the script writes only the focused
# paper-ready figure panels into the manuscript output tree.
# Tables are deliberately concise: significant sites plus compact method/count
# summaries, with the unsplit Original dataset represented once.
EXPORT_ONLY_PAPER_FIGURES <- FALSE
EXPORT_SUPPORT_FIGURES <- TRUE
EXPORT_INDIVIDUAL_VIEW_FIGURES <- FALSE

# Repository publication. After the complete analysis, export audit, and ZIP
# creation succeed, the manuscript output tree and the running pipeline script
# are committed and pushed to the current Git branch so the results are directly
# accessible from the GitHub repository.
PUBLISH_TO_GITHUB <- TRUE
GIT_REMOTE <- "origin"
GIT_COMMIT_PREFIX <- "Update SEQUENCE manuscript outputs"
GITHUB_MAX_FILE_MB <- 95

# -----------------------------------------------------------------------------
# Statistical decision rules
# -----------------------------------------------------------------------------
#   Std       = ordinary DESeq2 Wald BH padj < 0.10 AND |apeglm LFC| >= 1
#   Strong    = DESeq2 greaterAbs(lfcThreshold = 1) BH padj < 0.10
#               AND |apeglm LFC| >= 1
#   Weak-CNH  = DESeq2 lessAbs(lfcThreshold = 1) BH padj < 0.20
#   HBFSS     = |apeglm LFC| * [-log10(empirical p)] > Htau
#   Weak      = Weak-CNH AND |apeglm LFC| < 1 AND HBFSS
#   Overlap   = HBFSS AND (Std OR Strong OR Weak-CNH)
#
# HCp is the higher-criticism empirical-p threshold used to calculate Htau.
# It is not applied again as a second HBFSS significance gate.
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
  "sequence_journal_ready"
)

if (dir.exists(output_dir)) {
  unlink(output_dir, recursive = TRUE, force = TRUE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paper_fig_dir <- file.path(output_dir, "Combined_Figures")
summary_table_dir <- file.path(output_dir, "Summary_Tables")
dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_table_dir, recursive = TRUE, showWarnings = FALSE)

paper_registry <- list()

registry_key <- function(comparison_name, track_key) {
  paste(comparison_name, track_key, sep = "__")
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
    stop(sprintf("[%s] Fewer than 5 finite Wald statistics were available for fdrtool.", dataset_name))
  }

  fit <- fdrtool::fdrtool(
    stat_vec,
    statistic = "normal",
    plot = FALSE,
    verbose = FALSE,
    cutoff.method = "fndr"
  )
  fit$pval <- clip_probabilities(fit$pval)
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

significance_method_levels <- c(
  "Standard",
  "Strong",
  "Weak",
  "HBFSS"
)

significance_method_labels <- c(
  "Standard" = "Std",
  "Strong" = "Str",
  "Weak" = "Weak",
  "HBFSS" = "HBFSS"
)

significance_method_sizes <- c(
  "Standard" = 2.85,
  "Strong" = 3.15,
  "Weak" = 3.05,
  "HBFSS" = 2.55
)

significance_method_colors <- c(
  "Standard" = plot_palette$standard,
  "Strong" = plot_palette$strong,
  "Weak" = plot_palette$weak,
  "HBFSS" = plot_palette$hbfss
)

# Filled DESeq2 symbols plus a star for HBFSS keep the shared key visually
# obvious in every panel, including reduced mobile views and exported legends.
significance_method_shapes <- c(
  "Standard" = 18,
  "Strong" = 15,
  "Weak" = 17,
  "HBFSS" = 8
)

build_significance_plot_long <- function(df) {
  rows <- list()

  add_method <- function(flag_col, method_name) {
    keep <- !is.na(df[[flag_col]]) & df[[flag_col]]
    if (!any(keep)) return(NULL)
    out <- df[keep, , drop = FALSE]
    out$Method <- method_name
    out
  }

  rows[["Standard"]] <- add_method("standard_flag", "Standard")
  rows[["Strong"]] <- add_method("strong_cnh_flag", "Strong")
  rows[["Weak"]] <- add_method("weak_significant_flag", "Weak")
  rows[["HBFSS"]] <- add_method("hbfss_flag", "HBFSS")
  rows <- Filter(Negate(is.null), rows)

  if (!length(rows)) {
    out <- df[0, , drop = FALSE]
    out$Method <- factor(character(0), levels = significance_method_levels)
    return(out)
  }

  out <- dplyr::bind_rows(rows)
  out$Method <- factor(out$Method, levels = significance_method_levels)

  draw_rank <- c(HBFSS = 1L, Standard = 2L, Strong = 3L, Weak = 4L)
  out$.draw_rank <- unname(draw_rank[as.character(out$Method)])
  out <- out[order(out$.draw_rank), , drop = FALSE]
  out$.draw_rank <- NULL
  out
}

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


create_all_figures_zip <- function() {
  # The manuscript figure archive contains only nonredundant, comparison-level
  # panels plus the single global discovery-count summary.
  panel_files <- list.files(
    output_dir,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "Figure_.*_(PCA_EVS|Volcano_5Views|TWAS_5Views)\\.png$"
  )

  figure_files <- c(
    panel_files[
      grepl("/Panels/", gsub("\\\\", "/", panel_files))
    ],
    file.path(paper_fig_dir, "Figure_Manuscript_Discovery_Counts.png")
  )

  figure_files <- figure_files[
    file.exists(figure_files) &
      !(file.info(figure_files)$isdir %in% TRUE)
  ]

  if (!length(figure_files)) {
    stop("No curated manuscript PNG panels were found for SEQUENCE_ALL_FIGURES.zip.")
  }

  zip_path <- file.path(output_dir, "SEQUENCE_ALL_FIGURES.zip")
  if (file.exists(zip_path)) unlink(zip_path, force = TRUE)

  root_norm <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  file_norm <- normalizePath(unique(figure_files), winslash = "/", mustWork = TRUE)
  relative_files <- sort(substring(file_norm, nchar(root_norm) + 2L))

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)
  utils::zip(zipfile = basename(zip_path), files = relative_files, flags = "-q")

  if (!file.exists(zip_path)) stop("Figure ZIP creation failed: ", zip_path)
  normalizePath(zip_path, winslash = "/", mustWork = TRUE)
}


create_all_tables_zip <- function() {
  # The table archive contains concise manuscript-level results while preserving
  # a time-group-specific organization.
  wanted <- character(0)

  add_if_exists <- function(path) {
    if (file.exists(path)) wanted <<- c(wanted, path)
  }

  add_if_exists(file.path(summary_table_dir, "Table_DE_Method_Counts.csv"))
  add_if_exists(file.path(summary_table_dir, "Table_EVS_PCA_Evidence.csv"))

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    add_if_exists(file.path(
      output_dir,
      comparison_name,
      paste0("Table_", comparison_name, "_Significant_Sites_All_Views.csv")
    ))
    add_if_exists(file.path(
      output_dir,
      comparison_name,
      "Table_Method_Counts_All_Views.csv"
    ))
    add_if_exists(file.path(
      output_dir,
      comparison_name,
      paste0("Table_", comparison_name, "_EVS_PCA_Evidence.csv")
    ))
    add_if_exists(file.path(
      output_dir,
      comparison_name,
      paste0("Table_", comparison_name, "_TWAS_PAS_All_Views.csv")
    ))
    add_if_exists(file.path(
      output_dir,
      comparison_name,
      paste0("Table_", comparison_name, "_TWAS_Genes_All_Views.csv")
    ))
    add_if_exists(file.path(
      output_dir,
      comparison_name,
      paste0("Table_", comparison_name, "_TWAS_Method_Counts_All_Views.csv")
    ))
  }

  add_if_exists(file.path(output_dir, "TWAS", "tables", "Table_TWAS_Overlap_PAS.csv"))
  add_if_exists(file.path(output_dir, "TWAS", "tables", "Table_TWAS_Overlap_Genes.csv"))
  add_if_exists(file.path(output_dir, "TWAS", "tables", "Table_TWAS_Method_Counts.csv"))

  wanted <- sort(unique(wanted))
  if (!length(wanted)) stop("No manuscript tables were found for SEQUENCE_ALL_TABLES.zip.")

  zip_path <- file.path(output_dir, "SEQUENCE_ALL_TABLES.zip")
  if (file.exists(zip_path)) unlink(zip_path, force = TRUE)

  root_norm <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  file_norm <- normalizePath(wanted, winslash = "/", mustWork = TRUE)
  relative_files <- substring(file_norm, nchar(root_norm) + 2L)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)
  utils::zip(zipfile = basename(zip_path), files = relative_files, flags = "-q")

  if (!file.exists(zip_path)) stop("Table ZIP creation failed: ", zip_path)
  normalizePath(zip_path, winslash = "/", mustWork = TRUE)
}


# -----------------------------------------------------------------------------
# Repository publication helpers
# -----------------------------------------------------------------------------

git_exec <- function(args, fail = TRUE) {
  if (!nzchar(Sys.which("git"))) {
    if (isTRUE(fail)) stop("git is not installed or is not on PATH.")
    return(list(status = 127L, output = "git not found"))
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(repo_root)

  out <- suppressWarnings(
    system2("git", args = args, stdout = TRUE, stderr = TRUE)
  )
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L

  if (isTRUE(fail) && status != 0L) {
    stop(
      "Git command failed: git ", paste(args, collapse = " "), "\n",
      paste(out, collapse = "\n")
    )
  }

  list(status = as.integer(status), output = as.character(out))
}

repo_relative_path <- function(path) {
  root_norm <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)

  prefix <- paste0(root_norm, "/")
  if (!startsWith(path_norm, prefix) && path_norm != root_norm) {
    return(NA_character_)
  }

  if (path_norm == root_norm) return(".")
  substring(path_norm, nchar(prefix) + 1L)
}

write_repository_index <- function(figure_zip, table_zip) {
  index_path <- file.path(output_dir, "README.md")

  comparison_lines <- unlist(lapply(
    as.character(comparison_table$comparison_name),
    function(cmp) {
      c(
        paste0("### ", cmp),
        paste0("- [Panels](./", cmp, "/Panels/)"),
        paste0("- [All-view significant PAS table](./", cmp, "/Table_", cmp, "_Significant_Sites_All_Views.csv)"),
        paste0("- [All-view method counts](./", cmp, "/Table_Method_Counts_All_Views.csv)"),
        paste0("- [EVS/PCA evidence](./", cmp, "/Table_", cmp, "_EVS_PCA_Evidence.csv)"),
        paste0("- [TWAS PAS overlap](./", cmp, "/Table_", cmp, "_TWAS_PAS_All_Views.csv)"),
        paste0("- [TWAS gene overlap](./", cmp, "/Table_", cmp, "_TWAS_Genes_All_Views.csv)"),
        ""
      )
    }
  ))

  lines <- c(
    "# SEQUENCE manuscript outputs",
    "",
    paste0("Build: `", PIPELINE_BUILD, "`"),
    "",
    "This directory is generated only after the complete analysis and final export audit pass.",
    "",
    "## Downloadable archives",
    paste0("- [All manuscript figures](./", basename(figure_zip), ")"),
    paste0("- [All manuscript tables](./", basename(table_zip), ")"),
    "",
    "## Manuscript methods",
    "- [Methods_Manuscript.md](./Methods_Manuscript.md)",
    "",
    "## Time-group results",
    comparison_lines,
    "## Combined summaries",
    "- [Discovery counts](./Summary_Tables/Table_DE_Method_Counts.csv)",
    "- [EVS/PCA evidence](./Summary_Tables/Table_EVS_PCA_Evidence.csv)",
    "- [Combined TWAS PAS overlap](./TWAS/tables/Table_TWAS_Overlap_PAS.csv)",
    "- [Combined TWAS gene overlap](./TWAS/tables/Table_TWAS_Overlap_Genes.csv)",
    "- [Combined TWAS method counts](./TWAS/tables/Table_TWAS_Method_Counts.csv)",
    ""
  )

  writeLines(lines, index_path, useBytes = TRUE)
  normalizePath(index_path, winslash = "/", mustWork = TRUE)
}

copy_pipeline_snapshot <- function() {
  script_path <- get_script_path()
  if (is.na(script_path) || !file.exists(script_path)) return(NA_character_)

  snapshot_path <- file.path(
    output_dir,
    paste0("Pipeline_", basename(script_path))
  )

  ok <- file.copy(script_path, snapshot_path, overwrite = TRUE)
  if (!isTRUE(ok)) stop("Could not copy pipeline snapshot to output directory.")
  normalizePath(snapshot_path, winslash = "/", mustWork = TRUE)
}

assert_github_file_sizes <- function(paths) {
  paths <- unique(paths[file.exists(paths)])
  if (!length(paths)) return(invisible(TRUE))

  expanded <- unlist(lapply(paths, function(p) {
    if (dir.exists(p)) {
      list.files(p, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    } else {
      p
    }
  }))

  expanded <- unique(expanded[file.exists(expanded) & !dir.exists(expanded)])
  if (!length(expanded)) return(invisible(TRUE))

  size_mb <- file.info(expanded)$size / (1024^2)
  too_large <- expanded[is.finite(size_mb) & size_mb > GITHUB_MAX_FILE_MB]

  if (length(too_large)) {
    stop(
      "Repository publication stopped because GitHub-sized file limit would be exceeded: ",
      paste0(basename(too_large), collapse = ", "),
      ". Reduce these files or use Git LFS before publishing."
    )
  }

  invisible(TRUE)
}

github_repo_base_from_remote <- function(remote_url) {
  x <- trimws(remote_url)
  x <- sub("\\.git$", "", x)

  if (grepl("^https://github\\.com/", x)) return(x)
  if (grepl("^http://github\\.com/", x)) return(sub("^http://", "https://", x))

  if (grepl("^git@github\\.com:", x)) {
    return(sub("^git@github\\.com:", "https://github.com/", x))
  }

  if (grepl("^ssh://git@github\\.com/", x)) {
    return(sub("^ssh://git@github\\.com/", "https://github.com/", x))
  }

  NA_character_
}

publish_results_to_repository <- function(figure_zip, table_zip) {
  if (!isTRUE(PUBLISH_TO_GITHUB)) {
    return(list(published = FALSE, commit = NA_character_, url = NA_character_))
  }

  if (!dir.exists(file.path(repo_root, ".git"))) {
    stop("PUBLISH_TO_GITHUB is TRUE, but repository root has no .git directory: ", repo_root)
  }

  # Create a repository landing page and an exact snapshot of the pipeline that
  # generated these outputs before staging anything.
  write_repository_index(figure_zip, table_zip)
  copy_pipeline_snapshot()

  output_rel <- repo_relative_path(output_dir)
  if (is.na(output_rel)) stop("Output directory is not inside the Git repository.")

  script_path <- get_script_path()
  script_rel <- if (!is.na(script_path)) repo_relative_path(script_path) else NA_character_
  publish_paths <- unique(c(output_rel, script_rel[!is.na(script_rel)]))

  assert_github_file_sizes(c(output_dir, script_path[!is.na(script_path)]))

  branch_res <- git_exec(c("symbolic-ref", "--quiet", "--short", "HEAD"), fail = FALSE)
  if (branch_res$status != 0L || !length(branch_res$output) || !nzchar(branch_res$output[1])) {
    stop("Repository publication requires a named Git branch; detached HEAD is not supported.")
  }
  branch <- trimws(branch_res$output[1])

  remote_check <- git_exec(c("remote", "get-url", GIT_REMOTE), fail = FALSE)
  if (remote_check$status != 0L || !length(remote_check$output)) {
    stop("Git remote '", GIT_REMOTE, "' is not configured.")
  }
  remote_url <- trimws(remote_check$output[1])

  # Force-add generated exports even when an exports rule exists in .gitignore.
  git_exec(c("add", "-f", "--", publish_paths))

  diff_res <- git_exec(c("diff", "--cached", "--quiet", "--", publish_paths), fail = FALSE)

  commit_created <- FALSE
  if (diff_res$status == 1L) {
    commit_message <- paste0(
      GIT_COMMIT_PREFIX,
      " | ",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )

    # --only ensures unrelated files that happened to be staged before the run
    # are not included in the automated manuscript-output commit.
    git_exec(c(
      "commit", "--only", "-m", shQuote(commit_message), "--", publish_paths
    ))
    commit_created <- TRUE
  } else if (diff_res$status != 0L) {
    stop("Could not determine whether manuscript outputs changed in the Git index.")
  }

  # Push the current HEAD explicitly to the branch that ran the analysis. This
  # works even when the branch does not yet have an upstream tracking branch.
  git_exec(c("push", GIT_REMOTE, paste0("HEAD:", branch)))

  commit_res <- git_exec(c("rev-parse", "HEAD"))
  commit_hash <- trimws(commit_res$output[1])

  repo_base <- github_repo_base_from_remote(remote_url)
  github_url <- if (!is.na(repo_base)) {
    paste0(
      repo_base,
      "/tree/",
      utils::URLencode(branch, reserved = FALSE),
      "/",
      utils::URLencode(output_rel, reserved = FALSE)
    )
  } else {
    NA_character_
  }

  message(
    "Repository publication complete. Branch: ", branch,
    " | commit: ", substr(commit_hash, 1, 12),
    if (isTRUE(commit_created)) " | new commit created" else " | no new commit required"
  )

  list(
    published = TRUE,
    branch = branch,
    commit = commit_hash,
    url = github_url,
    output_rel = output_rel
  )
}

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




# -----------------------------------------------------------------------------
# Gene/PAS annotation used by the final 3'aTWAS analysis
# -----------------------------------------------------------------------------

valid_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & nzchar(x) & x != "-" & grepl("[A-Za-z]", x)
}

gene_key <- function(x) {
  x <- trimws(as.character(x))
  x[!valid_gene_symbol(x)] <- NA_character_
  toupper(x)
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  paste(sort(x), collapse = sep)
}

collapse_counts_by_label <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)

  tab <- table(x)
  tab <- tab[order(names(tab))]
  paste0(names(tab), "(", as.integer(tab), ")", collapse = ", ")
}

safe_min_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else min(x)
}

safe_max_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else max(x)
}

WTTS_gene_pas_summary <- OrigID_Symbol %>%
  dplyr::mutate(gene_key = gene_key(gene_symbol)) %>%
  dplyr::filter(!is.na(gene_key)) %>%
  dplyr::group_by(gene_key) %>%
  dplyr::summarise(
    WTTS_Gene = dplyr::first(gene_symbol[valid_gene_symbol(gene_symbol)]),
    WTTS_PAS_n = dplyr::n_distinct(orig_id),
    .groups = "drop"
  ) %>%
  dplyr::mutate(APA_multi_PAS = WTTS_PAS_n >= 2L)

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
  count_sub <- count_sub[rowSums(as.matrix(count_sub)) > 0, , drop = FALSE]

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
    downstream_deseq2_input = "raw-count feature subsets for DESeq2",
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
    alpha = BH_FDR_STANDARD,
    pAdjustMethod = "BH"
  )

  res_strong <- DESeq2::results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "greaterAbs",
    alpha = BH_FDR_STRONG,
    pAdjustMethod = "BH"
  )

  res_weak <- DESeq2::results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "lessAbs",
    alpha = BH_FDR_WEAK,
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

  if (length(fdr_fit$pval) != n_valid) {
    stop(sprintf("[%s] fdrtool empirical-p output length mismatch.", dataset_name))
  }

  # Preserve DESeq2's ordinary Wald pvalue and BH-adjusted padj exactly as
  # returned by DESeq2. Empirical p-values are stored in a separate column and
  # are used only for HBFSS/higher-criticism calculations and plotting.
  res_df$empirical_p <- NA_real_
  res_df$empirical_p[valid_stat] <- as.numeric(fdr_fit$pval)

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

  res_df$HBFSS <- abs(res_df$lfc_shrunk) *
    (-log10(pmax(res_df$empirical_p, 1e-300)))

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
    res_strong_df[, c("feature_id", "pvalue", "padj")],
    by = "feature_id",
    suffix = c("", "_strong")
  )

  res_df <- dplyr::left_join(
    res_df,
    res_weak_df[, c("feature_id", "pvalue", "padj")],
    by = "feature_id",
    suffix = c("", "_weak")
  )

  colnames(res_df)[colnames(res_df) == "pvalue_strong"] <- "pvalue_strong_effect"
  colnames(res_df)[colnames(res_df) == "padj_strong"] <- "padj_strong_effect"
  colnames(res_df)[colnames(res_df) == "pvalue_weak"] <- "pvalue_weak_effect"
  colnames(res_df)[colnames(res_df) == "padj_weak"] <- "padj_weak_effect"

  res_df$resGA_pvalue <- res_df$pvalue_strong_effect
  res_df$resGA_padj <- res_df$padj_strong_effect
  res_df$resLA_pvalue <- res_df$pvalue_weak_effect
  res_df$resLA_padj <- res_df$padj_weak_effect

  # Final decision flags. DESeq2 ordinary/greaterAbs/lessAbs p-values are
  # kept separate from empirical-null p-values. The manuscript Standard effect
  # combines the ordinary DESeq2 Wald/BH result with the prespecified effect
  # reporting boundary applied to the apeglm-shrunken LFC.
  finite_lfc <- !is.na(res_df$lfc_shrunk) & is.finite(res_df$lfc_shrunk)
  abs_shrunk_lfc <- abs(res_df$lfc_shrunk)

  res_df$standard_flag <- finite_lfc &
    abs_shrunk_lfc >= lfc_boundary &
    !is.na(res_df$padj) &
    res_df$padj < BH_FDR_STANDARD

  res_df$strong_cnh_flag <- finite_lfc &
    abs_shrunk_lfc >= lfc_boundary &
    !is.na(res_df$resGA_padj) &
    res_df$resGA_padj < BH_FDR_STRONG

  res_df$weak_cnh_flag <- finite_lfc &
    abs_shrunk_lfc < lfc_boundary &
    !is.na(res_df$resLA_padj) &
    res_df$resLA_padj < BH_FDR_WEAK

  res_df$hbfss_flag <- finite_lfc &
    !is.na(res_df$HBFSS_core_pass) &
    res_df$HBFSS_core_pass

  # lessAbs alone is not reported as differential expression. Final Weak is the
  # intersection of sub-boundary lessAbs evidence and HBFSS signal evidence.
  res_df$weak_significant_flag <- res_df$weak_cnh_flag & res_df$hbfss_flag

  res_df$standard_hbfss_overlap <- res_df$standard_flag & res_df$hbfss_flag
  res_df$strong_hbfss_overlap <- res_df$strong_cnh_flag & res_df$hbfss_flag
  res_df$weak_hbfss_overlap <- res_df$weak_cnh_flag & res_df$hbfss_flag
  res_df$any_overlap <- res_df$hbfss_flag &
    (res_df$standard_flag | res_df$strong_cnh_flag | res_df$weak_cnh_flag)

  res_df$final_significant_flag <- res_df$standard_flag |
    res_df$strong_cnh_flag |
    res_df$weak_significant_flag |
    res_df$hbfss_flag

  # Standard/Strong are outside the |LFC|=1 reporting boundary, whereas final
  # Weak is inside it. These sets must therefore be disjoint by construction.
  if (any(res_df$standard_flag & res_df$weak_significant_flag, na.rm = TRUE)) {
    stop("Standard and Weak classifications overlapped in ", dataset_name)
  }
  if (any(res_df$strong_cnh_flag & res_df$weak_significant_flag, na.rm = TRUE)) {
    stop("Strong and Weak classifications overlapped in ", dataset_name)
  }
  if (any(res_df$weak_significant_flag & !res_df$hbfss_flag, na.rm = TRUE)) {
    stop("A final Weak PAS lacked HBFSS support in ", dataset_name)
  }

  # Exact decision-rule checks. These are runtime assertions only; they do not
  # create validation tables or alter the reported results.
  expected_standard <- finite_lfc &
    abs_shrunk_lfc >= lfc_boundary &
    !is.na(res_df$padj) &
    res_df$padj < BH_FDR_STANDARD

  expected_strong <- finite_lfc &
    abs_shrunk_lfc >= lfc_boundary &
    !is.na(res_df$resGA_padj) &
    res_df$resGA_padj < BH_FDR_STRONG

  expected_weak_cnh <- finite_lfc &
    abs_shrunk_lfc < lfc_boundary &
    !is.na(res_df$resLA_padj) &
    res_df$resLA_padj < BH_FDR_WEAK

  expected_weak <- expected_weak_cnh & res_df$hbfss_flag
  expected_overlap <- res_df$hbfss_flag &
    (expected_standard | expected_strong | expected_weak_cnh)

  stopifnot(identical(res_df$standard_flag, expected_standard))
  stopifnot(identical(res_df$strong_cnh_flag, expected_strong))
  stopifnot(identical(res_df$weak_cnh_flag, expected_weak_cnh))
  stopifnot(identical(res_df$weak_significant_flag, expected_weak))
  stopifnot(identical(res_df$any_overlap, expected_overlap))

  expected_hbfss <- abs(res_df$lfc_shrunk) *
    (-log10(pmax(res_df$empirical_p, 1e-300)))
  both_na <- is.na(expected_hbfss) & is.na(res_df$HBFSS)
  both_finite <- is.finite(expected_hbfss) & is.finite(res_df$HBFSS)
  close_enough <- rep(FALSE, length(expected_hbfss))
  close_enough[both_finite] <- abs(expected_hbfss[both_finite] - res_df$HBFSS[both_finite]) <=
    1e-12 * pmax(1, abs(expected_hbfss[both_finite]))
  if (!all(both_na | close_enough)) stop("HBFSS arithmetic validation failed for ", dataset_name)

  if (is.finite(hc_p_threshold_dataset) && !is.na(hc_p_threshold_dataset)) {
    expected_htau <- -log10(hc_p_threshold_dataset) * lfc_boundary
    if (!isTRUE(all.equal(hbfss_threshold_dataset, expected_htau, tolerance = 1e-12))) {
      stop("HBFSS threshold arithmetic validation failed for ", dataset_name)
    }
  }

  if (is.finite(hbfss_threshold_dataset) && !is.na(hbfss_threshold_dataset)) {
    expected_hbfss_flag <- finite_lfc & !is.na(res_df$HBFSS) &
      is.finite(res_df$HBFSS) & res_df$HBFSS > hbfss_threshold_dataset
    stopifnot(identical(res_df$hbfss_flag, expected_hbfss_flag))
  }

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

  final_df$dataset_name <- dataset_name
  final_df$hc_p_threshold_dataset <- hc_p_threshold_dataset
  final_df$hbfss_threshold_dataset <- hbfss_threshold_dataset

  preferred_cols <- c(
    "dataset_name",
    "feature_id",
    "orig_id",
    "gene_symbol",
    "baseMean",
    "log2FoldChange",
    "lfc_shrunk",
    "regulation_direction",
    "stat",
    "pvalue",
    "padj",
    "empirical_p",
    "HBFSS",
    "hc_p_threshold_dataset",
    "hbfss_threshold_dataset",
    "resLA_pvalue",
    "resLA_padj",
    "resGA_pvalue",
    "resGA_padj",
    "weak_cnh_flag",
    "strong_cnh_flag",
    "standard_flag",
    "hbfss_flag",
    "weak_significant_flag",
    "final_significant_flag",
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
  df <- finite_plot_df(df, "lfc_shrunk", y_col)

  # Volcano labels are gene symbols only. Numeric PAS/feature identifiers are
  # never substituted into the figure when a gene symbol is absent.
  df$gene_symbol_plot <- if ("gene_symbol" %in% names(df)) {
    trimws(as.character(df$gene_symbol))
  } else {
    rep(NA_character_, nrow(df))
  }

  df$has_valid_gene_symbol <- !is.na(df$gene_symbol_plot) &
    nzchar(df$gene_symbol_plot) &
    grepl("[A-Za-z]", df$gene_symbol_plot) &
    !df$gene_symbol_plot %in% c("-", ".", "NA", "N/A")

  df[order(df$neglog10_empirical_p, na.last = TRUE), , drop = FALSE]
}

select_final_volcano_labels <- function(df, y_col, n_labels = n_top_labels_volcano) {
  if (!nrow(df)) return(df[0, , drop = FALSE])

  lab_df <- df[
    df$has_valid_gene_symbol &
      !is.na(df$final_significant_flag) &
      df$final_significant_flag,
    ,
    drop = FALSE
  ]

  if (!nrow(lab_df)) return(lab_df[0, , drop = FALSE])

  # Rank final discoveries by the evidence actually used on the HBFSS plotting
  # coordinates, then by effect magnitude. This does not change significance.
  emp <- suppressWarnings(as.numeric(lab_df$empirical_p))
  emp[!is.finite(emp)] <- Inf
  hscore <- suppressWarnings(as.numeric(lab_df$HBFSS))
  hscore[!is.finite(hscore)] <- -Inf

  ord <- order(emp, -hscore, -abs(lab_df$lfc_shrunk), na.last = TRUE)
  lab_df <- lab_df[ord, , drop = FALSE]
  lab_df <- lab_df[seq_len(min(as.integer(n_labels), nrow(lab_df))), , drop = FALSE]

  pas_id <- if ("orig_id" %in% names(lab_df)) {
    as.character(lab_df$orig_id)
  } else {
    as.character(lab_df$feature_id)
  }
  bad_pas <- is.na(pas_id) | !nzchar(trimws(pas_id))
  pas_id[bad_pas] <- as.character(lab_df$feature_id[bad_pas])

  gene_label <- as.character(lab_df$gene_symbol_plot)
  duplicate_gene <- duplicated(gene_label) | duplicated(gene_label, fromLast = TRUE)
  lab_df$plot_label <- ifelse(
    duplicate_gene,
    paste0(gene_label, " [", pas_id, "]"),
    gene_label
  )

  lab_df
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
  plot_df <- build_final_volcano_df(df, y_col = "neglog10_empirical_p")
  if (!nrow(plot_df)) stop("No finite volcano plotting rows for ", dataset_name)

  method_df <- build_significance_plot_long(plot_df)
  lab_df <- if (isTRUE(label_genes)) {
    select_final_volcano_labels(plot_df, y_col = "neglog10_empirical_p")
  } else {
    plot_df[0, , drop = FALSE]
  }

  hc_raw <- suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1]))
  htau <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))
  hc_y <- if (is.finite(hc_raw) && !is.na(hc_raw) && hc_raw > 0 && hc_raw <= 1) {
    safe_neglog10(hc_raw)
  } else NA_real_

  y_limit <- max(plot_df$neglog10_empirical_p, na.rm = TRUE) * 1.05
  boundary_df <- make_hbfss_boundary_df(plot_df, htau, y_limit)

  std_n <- sum(df$standard_flag, na.rm = TRUE)
  str_n <- sum(df$strong_cnh_flag, na.rm = TRUE)
  wk_n <- sum(df$weak_significant_flag, na.rm = TRUE)
  hbfss_n <- sum(df$hbfss_flag, na.rm = TRUE)
  ovlp_n <- sum(df$any_overlap, na.rm = TRUE)

  count_text <- paste0(
    "Std=", std_n,
    "  Str=", str_n,
    "  Wk=", wk_n,
    "  HBFSS=", hbfss_n,
    "  Ovlp=", ovlp_n,
    if (is.finite(hc_raw) && !is.na(hc_raw)) paste0("  HCp=", signif(hc_raw, 3)) else "",
    if (is.finite(htau) && !is.na(htau)) paste0("  Hτ=", signif(htau, 3)) else ""
  )

  plot_title <- if (is.null(short_title)) {
    compact_title(pretty_dataset_label(dataset_name), width = 42)
  } else short_title

  p <- ggplot() +
    geom_point(
      data = plot_df,
      aes(x = lfc_shrunk, y = neglog10_empirical_p),
      color = plot_palette$background,
      shape = 16,
      size = 0.48,
      alpha = 0.22
    ) +
    geom_point(
      data = method_df,
      aes(
        x = lfc_shrunk,
        y = neglog10_empirical_p,
        color = Method,
        shape = Method,
        size = Method
      ),
      alpha = 0.98,
      stroke = 0.90
    ) +
    scale_color_manual(
      values = significance_method_colors,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method",
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          shape = unname(significance_method_shapes[significance_method_levels]),
          color = unname(significance_method_colors[significance_method_levels]),
          size = rep(3.4, length(significance_method_levels)),
          alpha = rep(1, length(significance_method_levels)),
          stroke = rep(0.85, length(significance_method_levels))
        )
      )
    ) +
    scale_shape_manual(
      values = significance_method_shapes,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method",
      guide = "none"
    ) +
    scale_size_manual(
      values = significance_method_sizes,
      breaks = significance_method_levels,
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
      x = "apeglm shrunken log2FC",
      y = expression(-log[10]("Empirical p")),
      caption = compact_caption(count_text, width = 96)
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    plot_expand_xy() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.caption = element_text(
        size = base_theme_size - 2,
        hjust = 0.5,
        margin = margin(t = 4)
      ),
      plot.caption.position = "plot",
      plot.margin = margin(10, 12, 10, 10)
    )

  if (!is.na(hc_y) && is.finite(hc_y)) {
    p <- p + geom_hline(
      yintercept = hc_y,
      linetype = "dotted",
      linewidth = 0.65,
      color = plot_palette$hc
    )
  }

  if (!is.null(boundary_df)) {
    p <- p + geom_line(
      data = boundary_df,
      aes(x, y),
      inherit.aes = FALSE,
      color = plot_palette$hbfss_line,
      linewidth = 0.80
    )
  }

  if (nrow(lab_df) > 0L) {
    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      aes(x = lfc_shrunk, y = neglog10_empirical_p, label = plot_label),
      inherit.aes = FALSE,
      show.legend = FALSE,
      size = 2.05,
      color = "black",
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

run_one_evs_track <- function(comparison_name, track_key, count_matrix, coldata, annot_df) {
  evs <- build_eigenvector_split(
    count_matrix = count_matrix,
    coldata = coldata,
    track_key = track_key
  )

  message(
    "EVS track ", unname(track_short[track_key]),
    ": fixed top-5,000 selected independently from each condition PC1 eigenvector; downstream Lead/Rem DESeq2 analyses use the corresponding raw-count subsets."
  )

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

  # Original is independent of the EVS preprocessing path and is analyzed once.
  analysis_dataset_list <- dataset_list
  if (identical(track_key, "raw_evs")) analysis_dataset_list$raw_dataset <- NULL

  analysis_results <- list()
  for (nm in names(analysis_dataset_list)) {
    dataset_name <- paste(comparison_name, unname(track_short[track_key]), nm, sep = "_")
    fit <- run_core_analysis(
      count_mat = analysis_dataset_list[[nm]],
      coldata = coldata,
      dataset_name = dataset_name,
      annot_df = annot_df
    )

    analysis_results[[nm]] <- list(
      dds = fit$dds,
      results = fit$results,
      dataset_mat = analysis_dataset_list[[nm]]
    )
  }

  reg <- paper_registry[[registry_key(comparison_name, track_key)]]
  reg$analysis_results <- analysis_results
  paper_registry[[registry_key(comparison_name, track_key)]] <<- reg
  TRUE
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

    if (isTRUE(df$standard_flag[i])) tags <- c(tags, "Std")
    if (isTRUE(df$strong_cnh_flag[i])) tags <- c(tags, "Strong")
    if (isTRUE(df$weak_significant_flag[i])) tags <- c(tags, "Weak")
    if (isTRUE(df$hbfss_flag[i])) tags <- c(tags, "HBFSS")

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

analysis_view_folder_name <- function(track_key, dataset_key) {
  if (identical(dataset_key, "raw_dataset")) return("Original_No_EVS")
  paste0(
    unname(track_short[track_key]),
    "_",
    ifelse(dataset_key == "leading_edge_dataset", "Lead", "Rem")
  )
}

analysis_view_paths <- function(comparison_name, track_key, dataset_key) {
  base <- file.path(
    output_dir,
    comparison_name,
    analysis_view_folder_name(track_key, dataset_key)
  )
  fig <- file.path(base, "figures")
  tab <- file.path(base, "tables")
  dir.create(fig, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab, recursive = TRUE, showWarnings = FALSE)
  list(base = base, figures = fig, tables = tab)
}

export_comparison_manuscript_tables <- function(comparison_name) {
  views <- comparison_analysis_views()
  comparison_counts <- list()
  comparison_sig_rows <- list()

  for (i in seq_len(nrow(views))) {
    track_key <- views$track_key[i]
    dataset_key <- views$dataset_key[i]
    df <- get_registered_result(comparison_name, track_key, dataset_key)
    if (is.null(df) || !nrow(df)) next

    paths <- analysis_view_paths(comparison_name, track_key, dataset_key)
    analysis_label <- analysis_view_label(track_key, dataset_key)
    hc <- suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1]))
    htau <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))

    method_counts <- data.frame(
      Comparison = comparison_name,
      Analysis = analysis_label,
      PAS_tested = nrow(df),
      Std = sum(df$standard_flag, na.rm = TRUE),
      Strong = sum(df$strong_cnh_flag, na.rm = TRUE),
      Weak = sum(df$weak_significant_flag, na.rm = TRUE),
      HBFSS = sum(df$hbfss_flag, na.rm = TRUE),
      Ovlp = sum(df$any_overlap, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    expected_overlap <- sum(
      df$hbfss_flag & (df$standard_flag | df$strong_cnh_flag | df$weak_cnh_flag),
      na.rm = TRUE
    )
    if (method_counts$Ovlp[1] != expected_overlap) {
      stop("Overlap-count export mismatch: ", comparison_name, " / ", analysis_label)
    }

    save_csv(method_counts, file.path(paths$tables, "Table_Method_Counts.csv"))
    comparison_counts[[length(comparison_counts) + 1L]] <- method_counts

    sig <- df[
      !is.na(df$final_significant_flag) & df$final_significant_flag,
      ,
      drop = FALSE
    ]

    expected_sig_n <- sum(df$final_significant_flag, na.rm = TRUE)
    if (nrow(sig) != expected_sig_n) {
      stop("Significant-PAS export mismatch: ", comparison_name, " / ", analysis_label)
    }

    if (nrow(sig)) {
      pas <- if ("orig_id" %in% names(sig)) as.character(sig$orig_id) else as.character(sig$feature_id)
      bad_pas <- is.na(pas) | !nzchar(trimws(pas))
      pas[bad_pas] <- as.character(sig$feature_id[bad_pas])

      gene <- if ("gene_symbol" %in% names(sig)) as.character(sig$gene_symbol) else rep(NA_character_, nrow(sig))
      gene[is.na(gene) | !nzchar(trimws(gene))] <- NA_character_

      sig_out <- data.frame(
        Comparison = comparison_name,
        Analysis = analysis_label,
        PAS = pas,
        Gene = gene,
        Direction = as.character(sig$regulation_direction),
        Apeglm_LFC = as.numeric(sig$lfc_shrunk),
        Std_p = as.numeric(sig$pvalue),
        Std_BH = as.numeric(sig$padj),
        Strong_p = as.numeric(sig$resGA_pvalue),
        Strong_BH = as.numeric(sig$resGA_padj),
        Weak_p = as.numeric(sig$resLA_pvalue),
        Weak_BH = as.numeric(sig$resLA_padj),
        EmpP = as.numeric(sig$empirical_p),
        HBFSS = as.numeric(sig$HBFSS),
        HCp = hc,
        Htau = htau,
        Std = as.logical(sig$standard_flag),
        Strong = as.logical(sig$strong_cnh_flag),
        Weak = as.logical(sig$weak_significant_flag),
        HBFSS_sig = as.logical(sig$hbfss_flag),
        Ovlp = as.logical(sig$any_overlap),
        Support = build_support_label(sig),
        stringsAsFactors = FALSE
      )

      sig_out <- sig_out[
        order(sig_out$EmpP, -sig_out$HBFSS, -abs(sig_out$Apeglm_LFC), na.last = TRUE),
        ,
        drop = FALSE
      ]
    } else {
      sig_out <- data.frame()
    }

    save_csv(sig_out, file.path(paths$tables, "Table_Significant_Sites.csv"))
    if (nrow(sig_out)) comparison_sig_rows[[length(comparison_sig_rows) + 1L]] <- sig_out

    # Individual view figures are optional. Journal-ready review uses the
    # comparison-level panels generated after all five analysis views are fit.
    if (isTRUE(EXPORT_INDIVIDUAL_VIEW_FIGURES)) {
      dataset_name <- paste(comparison_name, unname(track_short[track_key]), dataset_key, sep = "_")
      volcano <- plot_final_volcano(
        df,
        dataset_name = dataset_name,
        short_title = paste0(comparison_name, " | ", analysis_label),
        label_genes = TRUE
      )
      save_grob(volcano, file.path(paths$figures, "Volcano.png"), width = 8.2, height = 6.4)
      save_grob(volcano, file.path(paths$figures, "Volcano.pdf"), width = 8.2, height = 6.4)
    }
  }

  comparison_counts <- if (length(comparison_counts)) dplyr::bind_rows(comparison_counts) else data.frame()
  if (nrow(comparison_counts)) {
    save_csv(
      comparison_counts,
      file.path(output_dir, comparison_name, "Table_Method_Counts_All_Views.csv")
    )
  }

  comparison_sig <- if (length(comparison_sig_rows)) {
    dplyr::bind_rows(comparison_sig_rows)
  } else {
    data.frame()
  }
  save_csv(
    comparison_sig,
    file.path(
      output_dir,
      comparison_name,
      paste0("Table_", comparison_name, "_Significant_Sites_All_Views.csv")
    )
  )

  invisible(comparison_counts)
}

build_overall_manuscript_summary <- function() {
  rows <- list()
  views <- comparison_analysis_views()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (i in seq_len(nrow(views))) {
      track_key <- views$track_key[i]
      dataset_key <- views$dataset_key[i]
      df <- get_registered_result(comparison_name, track_key, dataset_key)
      if (is.null(df) || !nrow(df)) next

      rows[[length(rows) + 1L]] <- data.frame(
        Comparison = comparison_name,
        Analysis = analysis_view_label(track_key, dataset_key),
        PAS_tested = nrow(df),
        Std = sum(df$standard_flag, na.rm = TRUE),
        Strong = sum(df$strong_cnh_flag, na.rm = TRUE),
        Weak = sum(df$weak_significant_flag, na.rm = TRUE),
        HBFSS = sum(df$hbfss_flag, na.rm = TRUE),
        Ovlp = sum(df$any_overlap, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) data.frame() else dplyr::bind_rows(rows)
}

write_methods_note <- function() {
  methods_text <- c(
    "# Manuscript Methods",
    "",
    "## Study unit and comparisons",
    "The analytical unit was the polyadenylation-site (PAS) feature defined by OrigID in the WTTS-Seq raw-count matrix. PASs were retained as separate observations throughout differential testing and gene symbols were retained as annotation. Four comparisons were analyzed: RT0 versus ZT6, RT2 versus ZT8, RT4 versus ZT10, and RT8 versus ZT14. Within each pairwise comparison, PASs with zero counts across every treatment and control sample were removed before EVS and differential-expression analysis; all remaining nonzero PASs were retained.",
    "",
    "## Eigenvector splitting",
    sprintf("NormEVS used DESeq2 median-of-ratios normalized counts before PCA, whereas RawEVS used raw counts before PCA. Within each condition, prcomp was applied to the transposed feature-by-sample matrix with centering and without scaling. PASs were ranked by absolute PC1 loading, and the %d highest-loading PASs were selected independently from the RT and ZT condition eigenvectors.", top_n_target),
    "PASs present in both condition-specific top-loading sets were classified as Joint; PASs present in only one set were classified as Disjoint for that condition. Joint plus both Disjoint sets formed the Leading Edge, and all other PASs formed the Remainder. EVS defined PAS membership only. Downstream DESeq2 analyses received the corresponding raw-count subset and estimated normalization and dispersion parameters within that analysis view. The five reported views were Original (No EVS), NormEVS Lead, NormEVS Rem, RawEVS Lead, and RawEVS Rem.",
    "",
    "## DESeq2 model and log2 fold-change shrinkage",
    "Each analysis view was modeled independently with DESeq2 using a negative-binomial generalized linear model with design ~ condition and ZT/untrt as the reference. DESeq2 estimated median-of-ratios size factors, gene-wise dispersions, the mean-dispersion relationship, final dispersions, and Wald statistics for the RT/trt coefficient. Log2 fold changes were then shrunken with apeglm. The apeglm-shrunken log2 fold change was the reported effect estimate, the effect term in HBFSS, the direction indicator, and the x-coordinate for all significance figures.",
    "",
    "## Standard effect",
    sprintf("The Standard effect used the ordinary two-sided DESeq2 Wald p-value with Benjamini-Hochberg adjustment at FDR %.2f and the prespecified manuscript effect boundary on the apeglm-shrunken estimate. A PAS was reported as Standard when padj < %.2f and |apeglm LFC| >= %.1f. The ordinary DESeq2 Wald p-value/padj were retained separately from the empirical-null p-value used for HBFSS.", BH_FDR_STANDARD, BH_FDR_STANDARD, lfc_boundary),
    "",
    "## Strong (decoupled) composite-null effect",
    sprintf("Strong (decoupled) effects were tested with DESeq2 results(..., lfcThreshold=%.1f, altHypothesis='greaterAbs'), followed by Benjamini-Hochberg adjustment at FDR %.2f. Reported Strong PASs additionally had |apeglm LFC| >= %.1f so the plotted/reported shrunken effect remained outside the stated effect boundary.", lfc_boundary, BH_FDR_STRONG, lfc_boundary),
    "",
    "## Weak composite-null effect",
    sprintf("Weak-effect support was tested with DESeq2 results(..., lfcThreshold=%.1f, altHypothesis='lessAbs'). Weak-CNH support required Benjamini-Hochberg adjusted p-value < %.2f and |apeglm LFC| < %.1f. A lessAbs rejection alone was not reported as differential expression. A PAS was reported as final Weak only when the same sub-boundary PAS also passed HBFSS.", lfc_boundary, BH_FDR_WEAK, lfc_boundary),
    "",
    "## Empirical-null calibration, higher criticism, and HBFSS",
    "The ordinary DESeq2 Wald statistics were supplied to fdrtool with statistic='normal' to estimate an empirical null distribution and empirical p-values. These empirical p-values were distinct from the ordinary DESeq2 Wald p-values and from the greaterAbs/lessAbs p-values. Higher criticism was applied to the empirical p-values with fdrtool::hc.thresh to obtain the dataset-specific HCp threshold.",
    sprintf("For each PAS, HBFSS = |apeglm-shrunken LFC| x [-log10(empirical p)]. With c=%.1f, Htau = -log10(HCp) x c. A PAS was HBFSS-significant when HBFSS > Htau. HBFSS significance did not use a DESeq2 adjusted-p-value gate.", lfc_boundary),
    "",
    "## Overlap",
    "Overlap was the number of unique HBFSS-significant PASs that also satisfied at least one DESeq2 criterion (Standard, Strong, or Weak-CNH). Overlap was a comparison quantity, not a separate significance test. Because final Weak required HBFSS, every final Weak PAS contributed to the overlap set.",
    "",
    "## PCA comparison after eigenvector splitting",
    "PCA comparison figures used the same matrices that defined EVS: full-comparison median-of-ratios normalized counts for NormEVS and raw counts for RawEVS. No additional log transformation or post-split re-normalization was applied. Original, Leading Edge, and Remainder were compared within the same preprocessing scale. For each view, the full PCA eigenspectrum was used to calculate the percentage of total variance explained by PC1 and PC2. PC1 eigenvalue retention was calculated as the first eigenvalue of the view divided by the first eigenvalue of the corresponding Original matrix. Original-PC1 loading energy captured by a feature subset was calculated as the sum of squared Original PC1 loadings for that subset divided by the sum of squared Original PC1 loadings for all PASs; Leading Edge and Remainder fractions were required to sum to 100%. Treatment-control separation was summarized in the PC1-PC2 plane as centroid distance divided by pooled within-group root-mean-square distance. The PCA/EVS manuscript panel combined the six Original/Leading Edge/Remainder PCA views with a 100% stacked summary of how Original PC1 squared-loading energy partitioned between Leading Edge and Remainder; exact numerical PCA/EVS metrics were retained in the corresponding evidence table.",
    "",
    "## Volcano figures and tables",
    sprintf("All significance volcanoes used the HBFSS plotting coordinates: apeglm-shrunken log2 fold change on the x-axis and -log10(empirical-null p) on the y-axis. Standard DESeq2, Strong greaterAbs, and Weak lessAbs significance were determined from their own DESeq2 p-values and BH-adjusted p-values, then their markers were projected onto these common HBFSS coordinates only for visual comparison. The HBFSS boundary y=Htau/|LFC| and HCp reference were drawn on the same axes. The same key was used throughout: Std = green diamond, Strong = red square, Weak = blue triangle, HBFSS = purple star. Multiple method markers were superimposed at the same PAS coordinate. LessAbs-only PASs remained background. Each volcano reported the number of significant PASs for Std, Strong, Weak, HBFSS, and their HBFSS/DESeq2 overlap and labeled at most %d final significant PASs.", n_top_labels_volcano),
    "Each comparison/view folder contained one significant-PAS table and one method-count table. Significant-PAS tables contained only the union of final Standard, Strong, Weak, and HBFSS discoveries and reported the DESeq2 p-values/BH-adjusted p-values used by the applicable DESeq2 tests together with empirical p, HBFSS, HCp, Htau, apeglm LFC, and explicit method indicators.",
    "",
    "## 3'aTWAS ortholog overlap",
    "Human 3'aTWAS gene symbols were mapped to rat gene symbols by combining database-supported babelgene human-to-rat ortholog mappings (top=FALSE) with direct case-insensitive symbol-equivalent matches present in the WTTS annotation. Each source TWAS row was assigned a unique record identifier before ortholog expansion. The expected many-to-many gene/ortholog relationship was retained explicitly, exact human-rat mapping duplicates were removed, and TWAS record/transcript counts were calculated from distinct source identifiers so ortholog expansion could not inflate counts. For every comparison and analysis view, mapped TWAS orthologs were intersected with Standard, Strong, Weak, and HBFSS WTTS discoveries. TWAS tables reported human TWAS symbols, rat orthologs, mapping source, TWAS disorder annotations, distinct TWAS records and PAS/transcript IDs, TWAS multi-PAS/APA status, total WTTS PAS count, WTTS multi-PAS/APA status, significant WTTS PAS identifiers, method-specific PAS identifiers/counts, significant multi-PAS status, and the method(s) and analysis view in which significance was observed.",
    "",
    "## Output organization",
    "Original (No EVS), NormEVS Lead, NormEVS Rem, RawEVS Lead, and RawEVS Rem were written to separate folders within each comparison for view-specific tables. Each RT/ZT comparison also contained combined all-view differential-expression and 3'aTWAS tables. Manuscript figures were restricted to nonredundant comparison-level PCA/EVS, volcano, and 3'aTWAS panels plus one global discovery-count panel. The curated figure archive SEQUENCE_ALL_FIGURES.zip contained those PNG panels only. The separate archive SEQUENCE_ALL_TABLES.zip contained the overall differential-expression method-count table, quantitative EVS/PCA evidence table, comparison-specific all-view significant-PAS/method-count tables, comparison-specific all-view TWAS PAS/gene/count tables, and the combined 3'aTWAS PAS-, gene-, and method-summary tables. After the final export audit passed, the complete manuscript output directory, downloadable figure and table archives, methods file, repository index, and the exact pipeline script used for the run were committed and pushed to the active repository branch."
  )

  writeLines(methods_text, file.path(output_dir, "Methods_Manuscript.md"), useBytes = TRUE)
  invisible(TRUE)
}

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  dir.create(file.path(output_dir, comparison_name), recursive = TRUE, showWarnings = FALSE)

  for (track_key in c("normalized_evs", "raw_evs")) {
    message("Running ", comparison_name, " ", unname(track_short[track_key]))
    ok <- run_one_evs_track(
      comparison_name = comparison_name,
      track_key = track_key,
      count_matrix = count_matrix,
      coldata = coldata,
      annot_df = annot_df
    )
    if (!isTRUE(ok)) stop("Analysis track failed: ", comparison_name, " / ", track_key)
  }

  export_comparison_manuscript_tables(comparison_name)
  TRUE
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

  df
}


save_paper_volcano_panels <- function() {
  dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)
  views <- comparison_analysis_views()
  short_view <- c(
    "Original (No EVS)" = "Orig",
    "NormEVS Lead" = "N-Lead",
    "NormEVS Rem" = "N-Rem",
    "RawEVS Lead" = "R-Lead",
    "RawEVS Rem" = "R-Rem"
  )

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    plots <- list()
    for (i in seq_len(nrow(views))) {
      track_key <- views$track_key[i]
      dataset_key <- views$dataset_key[i]
      df <- read_result_table_for_panel(comparison_name, track_key, dataset_key)
      if (is.null(df)) next
      analysis_label <- analysis_view_label(track_key, dataset_key)
      dataset_name <- paste(comparison_name, unname(track_short[track_key]), dataset_key, sep = "_")

      plots[[length(plots) + 1L]] <- plot_final_volcano(
        df = df,
        dataset_name = dataset_name,
        short_title = unname(short_view[analysis_label]),
        label_genes = TRUE
      ) +
        theme(
          plot.title = element_text(size = base_theme_size, face = "bold"),
          axis.title = element_text(size = base_theme_size - 1),
          axis.text = element_text(size = base_theme_size - 2),
          plot.caption = element_text(size = base_theme_size - 3)
        )
    }

    plots <- Filter(Negate(is.null), plots)
    if (!length(plots)) next

    panel <- assemble_one_legend_panel(
      plots,
      panel_title = paste0(comparison_name, " | significance across five analysis views"),
      ncol = 5
    )

    panel_dir <- file.path(output_dir, comparison_name, "Panels")
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
    base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_Volcano_5Views"))
    save_grob(panel, paste0(base, ".png"), width = 25.0, height = 7.6)
    save_grob(panel, paste0(base, ".pdf"), width = 25.0, height = 7.6)
  }

  invisible(TRUE)
}

pca_variance_stats <- function(value_df, coldata) {
  x <- as.matrix(value_df)
  storage.mode(x) <- "numeric"

  common_samples <- intersect(colnames(x), rownames(coldata))
  x <- x[, common_samples, drop = FALSE]
  coldata_use <- coldata[common_samples, , drop = FALSE]

  keep <- rowSums(is.finite(x)) == ncol(x) & apply(x, 1, stats::var) > 0
  x_use <- x[keep, , drop = FALSE]

  if (nrow(x_use) < 2L || ncol(x_use) < 3L) return(NULL)

  # Full PCA is used here so PC1/PC2 explained-variance percentages use the
  # complete non-zero eigenspectrum rather than only the first two components.
  fit <- stats::prcomp(
    t(x_use),
    center = TRUE,
    scale. = FALSE
  )

  eig <- fit$sdev^2
  total_var <- sum(eig)
  if (!is.finite(total_var) || total_var <= 0) return(NULL)

  score_df <- data.frame(
    Sample = rownames(fit$x),
    PC1 = fit$x[, 1],
    PC2 = fit$x[, 2],
    Condition = factor(
      as.character(coldata_use[rownames(fit$x), "condition"]),
      levels = c("untrt", "trt")
    ),
    stringsAsFactors = FALSE
  )

  centroids <- stats::aggregate(cbind(PC1, PC2) ~ Condition, data = score_df, FUN = mean)
  centroid_distance <- NA_real_
  if (nrow(centroids) == 2L) {
    centroid_distance <- sqrt(
      (centroids$PC1[1] - centroids$PC1[2])^2 +
        (centroids$PC2[1] - centroids$PC2[2])^2
    )
  }

  centroid_lookup <- merge(
    score_df,
    centroids,
    by = "Condition",
    suffixes = c("", "_centroid"),
    sort = FALSE
  )
  within_distance <- sqrt(
    (centroid_lookup$PC1 - centroid_lookup$PC1_centroid)^2 +
      (centroid_lookup$PC2 - centroid_lookup$PC2_centroid)^2
  )
  within_group_rms <- sqrt(mean(within_distance^2, na.rm = TRUE))
  separation_ratio <- if (
    is.finite(centroid_distance) && is.finite(within_group_rms) && within_group_rms > 0
  ) {
    centroid_distance / within_group_rms
  } else {
    NA_real_
  }

  list(
    fit = fit,
    coldata = coldata_use,
    score_df = score_df,
    centroids = centroids,
    n_features_input = nrow(x),
    n_features_pca = nrow(x_use),
    pc1_var = eig[1],
    pc2_var = if (length(eig) >= 2L) eig[2] else NA_real_,
    total_var = total_var,
    pc1_fraction = eig[1] / total_var,
    pc2_fraction = if (length(eig) >= 2L) eig[2] / total_var else NA_real_,
    centroid_distance = centroid_distance,
    within_group_rms = within_group_rms,
    separation_ratio = separation_ratio
  )
}

pc1_loading_energy_pct <- function(original_fit, feature_ids) {
  if (is.null(original_fit) || is.null(original_fit$rotation) || !ncol(original_fit$rotation)) {
    return(NA_real_)
  }

  load <- original_fit$rotation[, 1]
  denom <- sum(load^2, na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)

  ids <- intersect(as.character(feature_ids), names(load))
  100 * sum(load[ids]^2, na.rm = TRUE) / denom
}

pca_track_matrix <- function(registry_obj, dataset_key) {
  if (is.null(registry_obj) || is.null(registry_obj$evs)) return(NULL)
  evs <- registry_obj$evs
  mat <- as.data.frame(evs$rank_matrix)

  if (identical(dataset_key, "raw_dataset")) return(mat)
  if (identical(dataset_key, "leading_edge_dataset")) {
    return(mat[evs$leading_edge_ids, , drop = FALSE])
  }
  if (identical(dataset_key, "remainder_dataset")) {
    return(mat[evs$remainder_ids, , drop = FALSE])
  }
  stop("Unknown PCA dataset key: ", dataset_key)
}

compute_pca_support_plot <- function(value_df, coldata, short_title,
                                     original_stats = NULL) {
  stats_obj <- pca_variance_stats(value_df, coldata)
  if (is.null(stats_obj)) return(NULL)

  pca_df <- stats_obj$score_df
  centroid_df <- stats_obj$centroids

  retained <- if (
    !is.null(original_stats) &&
      is.finite(original_stats$pc1_var) &&
      original_stats$pc1_var > 0
  ) {
    100 * stats_obj$pc1_var / original_stats$pc1_var
  } else {
    100
  }

  energy <- if (!is.null(original_stats)) {
    pc1_loading_energy_pct(original_stats$fit, rownames(value_df))
  } else {
    100
  }

  # Segments from each sample to its condition centroid visualize within-group
  # dispersion; centroid crosses summarize treatment/control separation.
  segment_df <- merge(
    pca_df,
    centroid_df,
    by = "Condition",
    suffixes = c("", "_centroid"),
    sort = FALSE
  )

  cap <- paste0(
    "n=", stats_obj$n_features_input,
    " | PC1=", round(100 * stats_obj$pc1_fraction, 1), "%",
    " | λ1=", round(retained, 1), "% Orig",
    " | E1=", round(energy, 1), "%",
    " | Sep=", ifelse(is.finite(stats_obj$separation_ratio),
                       format(round(stats_obj$separation_ratio, 2), trim = TRUE), "NA")
  )

  ggplot(
    pca_df,
    aes(PC1, PC2, label = Sample, shape = Condition, fill = Condition)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey78") +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey78") +
    geom_segment(
      data = segment_df,
      aes(
        x = PC1,
        y = PC2,
        xend = PC1_centroid,
        yend = PC2_centroid,
        color = Condition
      ),
      inherit.aes = FALSE,
      linewidth = 0.32,
      alpha = 0.42,
      show.legend = FALSE
    ) +
    geom_point(size = 2.9, colour = "white", stroke = 0.55) +
    geom_point(
      data = centroid_df,
      aes(PC1, PC2, color = Condition),
      inherit.aes = FALSE,
      shape = 4,
      stroke = 1.15,
      size = 3.5,
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      size = 1.65,
      max.overlaps = 10,
      force = 0.9,
      box.padding = 0.16,
      point.padding = 0.08,
      min.segment.length = 0,
      segment.alpha = 0.45,
      segment.size = 0.16
    ) +
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition") +
    scale_fill_manual(values = condition_fills, labels = condition_labels, name = "Condition") +
    scale_color_manual(values = condition_fills, guide = "none") +
    labs(
      title = short_title,
      x = "PC1",
      y = "PC2",
      caption = cap
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(
      legend.position = "bottom",
      plot.caption = element_text(size = base_theme_size - 3, hjust = 0.5),
      plot.margin = margin(8, 10, 8, 10)
    )
}

build_evs_pca_evidence_table <- function(comparison_name) {
  rows <- list()

  for (track_key in c("normalized_evs", "raw_evs")) {
    obj <- paper_registry[[registry_key(comparison_name, track_key)]]
    if (is.null(obj) || is.null(obj$evs)) next

    original_mat <- pca_track_matrix(obj, "raw_dataset")
    original_stats <- pca_variance_stats(original_mat, obj$coldata)
    if (is.null(original_stats)) next

    track_rows <- list()

    for (dataset_key in c("raw_dataset", "leading_edge_dataset", "remainder_dataset")) {
      mat <- pca_track_matrix(obj, dataset_key)
      st <- pca_variance_stats(mat, obj$coldata)
      if (is.null(st)) next

      energy <- pc1_loading_energy_pct(original_stats$fit, rownames(mat))
      dataset_label <- c(
        raw_dataset = "Orig",
        leading_edge_dataset = "Lead",
        remainder_dataset = "Rem"
      )[[dataset_key]]

      track_rows[[dataset_key]] <- data.frame(
        Comparison = comparison_name,
        Track = unname(track_short[track_key]),
        Dataset = dataset_label,
        PAS_n = nrow(mat),
        PCA_PAS_n = st$n_features_pca,
        PC1_explained_pct = 100 * st$pc1_fraction,
        PC2_explained_pct = 100 * st$pc2_fraction,
        PC1_eigenvalue = st$pc1_var,
        PC1_eigenvalue_retained_pct = 100 * st$pc1_var / original_stats$pc1_var,
        Original_PC1_loading_energy_pct = energy,
        Centroid_distance_PC1_PC2 = st$centroid_distance,
        Within_group_RMS_PC1_PC2 = st$within_group_rms,
        Separation_ratio = st$separation_ratio,
        stringsAsFactors = FALSE
      )
    }

    track_tbl <- dplyr::bind_rows(track_rows)

    # Lead and Rem form a partition of the Original feature set. Their shares of
    # Original PC1 squared-loading energy must therefore sum to 100% apart from
    # floating-point tolerance. This assertion catches membership or ID errors.
    lr <- track_tbl[track_tbl$Dataset %in% c("Lead", "Rem"), , drop = FALSE]
    if (nrow(lr) == 2L && all(is.finite(lr$Original_PC1_loading_energy_pct))) {
      energy_sum <- sum(lr$Original_PC1_loading_energy_pct)
      if (abs(energy_sum - 100) > 1e-6) {
        stop(
          comparison_name, " ", unname(track_short[track_key]),
          ": Lead + Rem Original-PC1 loading energy did not sum to 100%."
        )
      }
    }

    rows[[length(rows) + 1L]] <- track_tbl
  }

  if (!length(rows)) data.frame() else dplyr::bind_rows(rows)
}


plot_evs_pca_evidence <- function(comparison_name) {
  tbl <- build_evs_pca_evidence_table(comparison_name)
  if (!nrow(tbl)) return(NULL)

  split_tbl <- tbl[tbl$Dataset %in% c("Lead", "Rem"), , drop = FALSE]
  if (!nrow(split_tbl)) return(NULL)

  split_tbl$Dataset <- factor(
    split_tbl$Dataset,
    levels = c("Lead", "Rem"),
    labels = c("Leading Edge", "Remainder")
  )

  ggplot(
    split_tbl,
    aes(
      x = Track,
      y = Original_PC1_loading_energy_pct,
      fill = Dataset
    )
  ) +
    geom_col(
      width = 0.62,
      position = "stack",
      color = "white",
      linewidth = 0.35
    ) +
    geom_text(
      aes(
        label = ifelse(
          Original_PC1_loading_energy_pct >= 1,
          paste0(round(Original_PC1_loading_energy_pct, 1), "%"),
          paste0(format(round(Original_PC1_loading_energy_pct, 2), nsmall = 2), "%")
        )
      ),
      position = position_stack(vjust = 0.5),
      size = 3.2,
      fontface = "bold"
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      title = paste0(comparison_name, " | original PC1 signal partition after EVS"),
      subtitle = "Squared-loading energy of the Original PC1 eigenvector",
      x = NULL,
      y = "Original PC1 loading energy (%)",
      fill = NULL
    ) +
    manuscript_theme() +
    theme(
      legend.position = "bottom",
      plot.margin = margin(8, 12, 8, 12)
    )
}


save_paper_pca_panels <- function() {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) return(invisible(FALSE))

  all_evidence <- list()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    plots <- list()

    for (track_key in c("normalized_evs", "raw_evs")) {
      obj <- paper_registry[[registry_key(comparison_name, track_key)]]
      if (is.null(obj) || is.null(obj$evs)) next

      original_mat <- pca_track_matrix(obj, "raw_dataset")
      original_stats <- pca_variance_stats(original_mat, obj$coldata)
      if (is.null(original_stats)) next

      track_abbr <- if (track_key == "normalized_evs") "Norm" else "Raw"

      for (dataset_key in c("raw_dataset", "leading_edge_dataset", "remainder_dataset")) {
        mat <- pca_track_matrix(obj, dataset_key)
        ds_abbr <- c(
          raw_dataset = "Orig",
          leading_edge_dataset = "Lead",
          remainder_dataset = "Rem"
        )[[dataset_key]]

        plots[[length(plots) + 1L]] <- compute_pca_support_plot(
          value_df = mat,
          coldata = obj$coldata,
          short_title = paste0(track_abbr, "-", ds_abbr),
          original_stats = original_stats
        )
      }
    }

    plots <- Filter(Negate(is.null), plots)
    panel_dir <- file.path(output_dir, comparison_name, "Panels")
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

    evidence_tbl <- build_evs_pca_evidence_table(comparison_name)
    if (nrow(evidence_tbl)) {
      all_evidence[[comparison_name]] <- evidence_tbl
      save_csv(
        evidence_tbl,
        file.path(output_dir, comparison_name, paste0("Table_", comparison_name, "_EVS_PCA_Evidence.csv"))
      )
    }

    pevidence <- plot_evs_pca_evidence(comparison_name)

    if (length(plots)) {
      pca_grid <- assemble_one_legend_panel(
        plots,
        panel_title = paste0(comparison_name, " | PCA structure before and after EVS"),
        ncol = 3
      )

      combined <- if (!is.null(pevidence)) {
        gridExtra::arrangeGrob(
          pca_grid,
          pevidence,
          ncol = 1,
          heights = c(2.8, 1.0)
        )
      } else {
        pca_grid
      }

      base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_PCA_EVS"))
      save_grob(combined, paste0(base, ".png"), width = 17.0, height = 13.6)
      save_grob(combined, paste0(base, ".pdf"), width = 17.0, height = 13.6)
    }
  }

  evidence_all <- if (length(all_evidence)) dplyr::bind_rows(all_evidence) else data.frame()
  if (nrow(evidence_all)) {
    save_csv(evidence_all, file.path(summary_table_dir, "Table_EVS_PCA_Evidence.csv"))
  }

  invisible(evidence_all)
}

build_discovery_long_table <- function(summary_df) {
  if (!is.data.frame(summary_df) || nrow(summary_df) == 0L) return(data.frame())

  rows <- lapply(seq_len(nrow(summary_df)), function(i) {
    sm <- summary_df[i, , drop = FALSE]
    data.frame(
      Comparison = sm$Comparison,
      Analysis = sm$Analysis,
      Method = factor(
        c("Standard", "Strong", "Weak", "HBFSS"),
        levels = significance_method_levels
      ),
      Count = as.numeric(c(sm$Std, sm$Strong, sm$Weak, sm$HBFSS)),
      Ovlp = as.numeric(sm$Ovlp),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}


save_discovery_count_panel <- function(summary_df) {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) return(invisible(FALSE))
  long_df <- build_discovery_long_table(summary_df)
  if (!nrow(long_df)) return(invisible(FALSE))

  long_df$Method <- factor(long_df$Method, levels = significance_method_levels)

  p <- ggplot(
    long_df,
    aes(Comparison, Count, color = Method, shape = Method)
  ) +
    geom_point(
      position = position_dodge(width = 0.58),
      size = 3.2,
      alpha = 0.98,
      stroke = 0.90
    ) +
    geom_text(
      aes(label = Count),
      position = position_dodge(width = 0.58),
      vjust = -0.65,
      size = 2.7,
      color = "black",
      show.legend = FALSE
    ) +
    facet_wrap(~ Analysis, scales = "free_y", ncol = 5) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    scale_color_manual(
      values = significance_method_colors,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method",
      guide = guide_legend(override.aes = list(size = 3.0, stroke = 0.95))
    ) +
    scale_shape_manual(
      values = significance_method_shapes,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method"
    ) +
    labs(
      title = "Significant PAS counts by method",
      x = NULL,
      y = "Significant PASs",
      caption = "Ovlp is reported separately; method totals are not mutually exclusive."
    ) +
    manuscript_theme() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(face = "bold")
    )

  save_grob(p, file.path(paper_fig_dir, "Figure_Manuscript_Discovery_Counts.png"), width = 20.0, height = 5.8)
  save_grob(p, file.path(paper_fig_dir, "Figure_Manuscript_Discovery_Counts.pdf"), width = 20.0, height = 5.8)
  invisible(TRUE)
}


save_paper_support_figures <- function(summary_df) {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) {
    return(invisible(FALSE))
  }

  support_steps <- list(
    PCA = function() save_paper_pca_panels(),
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

ensure_babelgene <- function() {
  if (requireNamespace("babelgene", quietly = TRUE)) {
    return(TRUE)
  }

  message("Installing CRAN package 'babelgene' for database-supported human-to-rat ortholog mapping...")
  suppressWarnings(
    try(
      utils::install.packages(
        "babelgene",
        repos = "https://cloud.r-project.org",
        quiet = TRUE
      ),
      silent = TRUE
    )
  )

  if (!requireNamespace("babelgene", quietly = TRUE)) {
    warning(
      "babelgene could not be installed. TWAS overlap will still run using direct case-insensitive human/rat symbol matches, but database-supported non-identical ortholog symbols cannot be added in this run."
    )
    return(FALSE)
  }

  TRUE
}


read_twas_study <- function(path) {
  twas_raw <- utils::read.csv(
    path,
    skip = 1,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  assert_required_columns(
    twas_raw,
    c("Disease", "PANEL", "Transcript ID", "Gene symbol", "3'aTWAS.Z", "3'aTWAS.P"),
    object_name = "3'aTWAS file"
  )

  twas_raw$TWAS_Source_Row <- seq_len(nrow(twas_raw))

  twas_raw %>%
    dplyr::transmute(
      TWAS_Record_ID = paste0("TWAS_", TWAS_Source_Row),
      TWAS_Source_Row = as.integer(TWAS_Source_Row),
      Disease = trimws(as.character(Disease)),
      Panel = trimws(as.character(PANEL)),
      TWAS_Transcript = trimws(as.character(`Transcript ID`)),
      Human_TWAS_Gene = trimws(as.character(`Gene symbol`)),
      TWAS_Z = suppressWarnings(as.numeric(`3'aTWAS.Z`)),
      TWAS_P = suppressWarnings(as.numeric(`3'aTWAS.P`)),
      COLOC_PP4 = if ("COLOC.PP4" %in% names(twas_raw)) suppressWarnings(as.numeric(COLOC.PP4)) else NA_real_
    ) %>%
    dplyr::filter(valid_gene_symbol(Human_TWAS_Gene)) %>%
    dplyr::distinct(TWAS_Record_ID, .keep_all = TRUE)
}

build_twas_ortholog_map <- function(twas) {
  human_genes <- sort(unique(twas$Human_TWAS_Gene))

  # Direct symbol equivalence is retained for genes whose human symbol differs
  # from the rat WTTS symbol only by capitalization. This preserves obvious
  # one-to-one symbol matches and is combined with database-supported orthology.
  rat_symbols <- sort(unique(OrigID_Symbol$gene_symbol[valid_gene_symbol(OrigID_Symbol$gene_symbol)]))
  rat_lookup <- data.frame(
    gene_key = gene_key(rat_symbols),
    Rat_Ortholog = rat_symbols,
    stringsAsFactors = FALSE
  ) %>% dplyr::filter(!is.na(gene_key)) %>% dplyr::distinct(gene_key, .keep_all = TRUE)

  direct <- data.frame(
    Human_TWAS_Gene = human_genes,
    gene_key = gene_key(human_genes),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::inner_join(rat_lookup, by = "gene_key") %>%
    dplyr::transmute(
      Human_TWAS_Gene,
      Rat_Ortholog,
      Ortholog_support_n = NA_integer_,
      Ortholog_support = "case-insensitive symbol match",
      Mapping_source = "symbol"
    )

  babel <- data.frame()
  if (ensure_babelgene()) {
    orth <- tryCatch(
      babelgene::orthologs(
        genes = human_genes,
        species = TWAS_TARGET_SPECIES,
        human = TRUE,
        min_support = TWAS_ORTHOLOG_MIN_SUPPORT,
        top = FALSE
      ),
      error = function(e) {
        warning("babelgene ortholog mapping failed; continuing with direct symbol matches: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(orth) && nrow(orth)) {
      orth <- as.data.frame(orth, stringsAsFactors = FALSE)
      assert_required_columns(
        orth,
        c("human_symbol", "symbol", "support_n"),
        object_name = "babelgene ortholog output"
      )

      babel <- data.frame(
        Human_TWAS_Gene = as.character(orth$human_symbol),
        Rat_Ortholog = as.character(orth$symbol),
        Ortholog_support_n = suppressWarnings(as.integer(orth$support_n)),
        Ortholog_support = if ("support" %in% names(orth)) as.character(orth$support) else NA_character_,
        Mapping_source = "babelgene",
        stringsAsFactors = FALSE
      ) %>%
        dplyr::filter(
          Human_TWAS_Gene %in% human_genes,
          valid_gene_symbol(Rat_Ortholog)
        )
    }
  }

  dplyr::bind_rows(babel, direct) %>%
    dplyr::mutate(gene_key = gene_key(Rat_Ortholog)) %>%
    dplyr::filter(!is.na(gene_key)) %>%
    dplyr::arrange(Human_TWAS_Gene, dplyr::desc(Ortholog_support_n), Mapping_source) %>%
    dplyr::distinct(Human_TWAS_Gene, Rat_Ortholog, .keep_all = TRUE)
}


build_twas_rat_metadata <- function(twas, mapping) {
  # A human TWAS gene can legitimately have multiple TWAS records and can also
  # map to more than one supported rat ortholog. The relationship is therefore
  # explicitly many-to-many. Counts are subsequently computed from distinct
  # source-record IDs/transcripts so ortholog expansion cannot inflate them.
  mapped <- twas %>%
    dplyr::left_join(
      mapping,
      by = "Human_TWAS_Gene",
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(!is.na(gene_key), valid_gene_symbol(Rat_Ortholog)) %>%
    dplyr::distinct(
      TWAS_Record_ID,
      Human_TWAS_Gene,
      Rat_Ortholog,
      .keep_all = TRUE
    )

  mapped %>%
    dplyr::group_by(gene_key, Rat_Ortholog) %>%
    dplyr::summarise(
      Human_TWAS_Genes = collapse_unique(Human_TWAS_Gene),
      TWAS_Disorders = collapse_counts_by_label(Disease),
      TWAS_n = dplyr::n_distinct(TWAS_Record_ID),
      TWAS_Panel_n = dplyr::n_distinct(Panel, na.rm = TRUE),
      TWAS_PAS_n = dplyr::n_distinct(TWAS_Transcript, na.rm = TRUE),
      TWAS_APA = dplyr::n_distinct(TWAS_Transcript, na.rm = TRUE) >= 2L,
      TWAS_min_P = safe_min_numeric(TWAS_P),
      TWAS_max_abs_Z = safe_max_numeric(abs(TWAS_Z)),
      TWAS_max_COLOC_PP4 = safe_max_numeric(COLOC_PP4),
      Ortholog_support_n = safe_max_numeric(Ortholog_support_n),
      Ortholog_mapping_source = collapse_unique(Mapping_source),
      .groups = "drop"
    ) %>%
    dplyr::left_join(WTTS_gene_pas_summary, by = "gene_key")
}

collect_twas_analysis_results <- function() {
  views <- comparison_analysis_views()
  rows <- list()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (i in seq_len(nrow(views))) {
      track_key <- views$track_key[i]
      dataset_key <- views$dataset_key[i]
      df <- get_registered_result(comparison_name, track_key, dataset_key)
      if (is.null(df) || !nrow(df)) next

      pas <- if ("orig_id" %in% names(df)) as.character(df$orig_id) else as.character(df$feature_id)
      bad_pas <- is.na(pas) | !nzchar(trimws(pas))
      pas[bad_pas] <- as.character(df$feature_id[bad_pas])

      rows[[length(rows) + 1L]] <- data.frame(
        Comparison = comparison_name,
        Analysis = analysis_view_label(track_key, dataset_key),
        feature_id = as.character(df$feature_id),
        PAS = pas,
        WTTS_Gene = as.character(df$gene_symbol),
        gene_key = gene_key(df$gene_symbol),
        Apeglm_LFC = suppressWarnings(as.numeric(df$lfc_shrunk)),
        Std_p = suppressWarnings(as.numeric(df$pvalue)),
        Std_BH = suppressWarnings(as.numeric(df$padj)),
        Strong_p = suppressWarnings(as.numeric(df$resGA_pvalue)),
        Strong_BH = suppressWarnings(as.numeric(df$resGA_padj)),
        Weak_p = suppressWarnings(as.numeric(df$resLA_pvalue)),
        Weak_BH = suppressWarnings(as.numeric(df$resLA_padj)),
        EmpP = suppressWarnings(as.numeric(df$empirical_p)),
        HBFSS_score = suppressWarnings(as.numeric(df$HBFSS)),
        HCp = suppressWarnings(as.numeric(df$hc_p_threshold_dataset)),
        Htau = suppressWarnings(as.numeric(df$hbfss_threshold_dataset)),
        Std = !is.na(df$standard_flag) & df$standard_flag,
        Strong = !is.na(df$strong_cnh_flag) & df$strong_cnh_flag,
        Weak = !is.na(df$weak_significant_flag) & df$weak_significant_flag,
        HBFSS = !is.na(df$hbfss_flag) & df$hbfss_flag,
        Ovlp = !is.na(df$any_overlap) & df$any_overlap,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) return(data.frame())
  dplyr::bind_rows(rows) %>% dplyr::filter(!is.na(gene_key))
}


build_twas_overlap_pas_table <- function(all_results, twas_rat_meta) {
  if (!nrow(all_results) || !nrow(twas_rat_meta)) return(data.frame())

  all_results %>%
    dplyr::filter(Std | Strong | Weak | HBFSS) %>%
    dplyr::inner_join(
      twas_rat_meta,
      by = "gene_key",
      relationship = "many-to-many"
    ) %>%
    dplyr::mutate(
      Support = vapply(seq_len(dplyr::n()), function(i) {
        tags <- c(
          if (Std[i]) "Std" else NULL,
          if (Strong[i]) "Str" else NULL,
          if (Weak[i]) "Weak" else NULL,
          if (HBFSS[i]) "HBFSS" else NULL
        )
        paste(tags, collapse = "+")
      }, character(1)),
      Analysis = factor(
        Analysis,
        levels = c(
          "Original (No EVS)",
          "NormEVS Lead", "NormEVS Rem",
          "RawEVS Lead", "RawEVS Rem"
        )
      )
    ) %>%
    dplyr::transmute(
      Comparison,
      Analysis = as.character(Analysis),
      Human_TWAS = Human_TWAS_Genes,
      Rat_Ortholog,
      Ortholog_mapping_source,
      TWAS_Disorders,
      TWAS_n,
      TWAS_PAS_n,
      TWAS_APA,
      WTTS_PAS_n,
      WTTS_APA = APA_multi_PAS,
      WTTS_DE_PAS = PAS,
      Apeglm_LFC,
      Std,
      Strong,
      Weak,
      HBFSS,
      Ovlp,
      Support,
      Std_BH,
      Strong_BH,
      Weak_BH,
      EmpP,
      HBFSS_score
    ) %>%
    dplyr::distinct(
      Comparison,
      Analysis,
      Rat_Ortholog,
      WTTS_DE_PAS,
      .keep_all = TRUE
    ) %>%
    dplyr::arrange(
      Comparison,
      factor(
        Analysis,
        levels = c(
          "Original (No EVS)",
          "NormEVS Lead", "NormEVS Rem",
          "RawEVS Lead", "RawEVS Rem"
        )
      ),
      Rat_Ortholog,
      WTTS_DE_PAS
    )
}


build_twas_overlap_gene_table <- function(pas_table) {
  if (!nrow(pas_table)) return(data.frame())

  pas_table %>%
    dplyr::mutate(
      Std_flag = Std,
      Strong_flag = Strong,
      Weak_flag = Weak,
      HBFSS_flag = HBFSS,
      Ovlp_flag = Ovlp
    ) %>%
    dplyr::group_by(
      Comparison,
      Analysis,
      Human_TWAS,
      Rat_Ortholog,
      Ortholog_mapping_source,
      TWAS_Disorders,
      TWAS_n,
      TWAS_PAS_n,
      TWAS_APA,
      WTTS_PAS_n,
      WTTS_APA
    ) %>%
    dplyr::summarise(
      Std = any(Std_flag, na.rm = TRUE),
      Strong = any(Strong_flag, na.rm = TRUE),
      Weak = any(Weak_flag, na.rm = TRUE),
      HBFSS = any(HBFSS_flag, na.rm = TRUE),
      Ovlp = any(Ovlp_flag, na.rm = TRUE),
      Std_PAS_n = dplyr::n_distinct(WTTS_DE_PAS[Std_flag]),
      Strong_PAS_n = dplyr::n_distinct(WTTS_DE_PAS[Strong_flag]),
      Weak_PAS_n = dplyr::n_distinct(WTTS_DE_PAS[Weak_flag]),
      HBFSS_PAS_n = dplyr::n_distinct(WTTS_DE_PAS[HBFSS_flag]),
      Std_PAS_IDs = collapse_unique(WTTS_DE_PAS[Std_flag]),
      Strong_PAS_IDs = collapse_unique(WTTS_DE_PAS[Strong_flag]),
      Weak_PAS_IDs = collapse_unique(WTTS_DE_PAS[Weak_flag]),
      HBFSS_PAS_IDs = collapse_unique(WTTS_DE_PAS[HBFSS_flag]),
      WTTS_DE_PAS_n = dplyr::n_distinct(WTTS_DE_PAS),
      WTTS_DE_APA = dplyr::n_distinct(WTTS_DE_PAS) >= 2L,
      WTTS_DE_PAS_IDs = collapse_unique(WTTS_DE_PAS),
      Methods = paste(
        c(
          if (any(Std_flag, na.rm = TRUE)) "Std" else NULL,
          if (any(Strong_flag, na.rm = TRUE)) "Str" else NULL,
          if (any(Weak_flag, na.rm = TRUE)) "Weak" else NULL,
          if (any(HBFSS_flag, na.rm = TRUE)) "HBFSS" else NULL
        ),
        collapse = "+"
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      Comparison,
      factor(
        Analysis,
        levels = c(
          "Original (No EVS)",
          "NormEVS Lead", "NormEVS Rem",
          "RawEVS Lead", "RawEVS Rem"
        )
      ),
      Rat_Ortholog
    )
}


build_twas_overlap_summary <- function(gene_table) {
  view_levels <- c(
    "Original (No EVS)",
    "NormEVS Lead",
    "NormEVS Rem",
    "RawEVS Lead",
    "RawEVS Rem"
  )

  grid <- expand.grid(
    Comparison = as.character(comparison_table$comparison_name),
    Analysis = view_levels,
    stringsAsFactors = FALSE
  )

  if (!nrow(gene_table)) {
    out <- grid
    for (nm in c("Std", "Strong", "Weak", "HBFSS", "Ovlp", "TWAS_APA", "WTTS_APA", "WTTS_DE_APA")) {
      out[[nm]] <- 0L
    }
    return(out)
  }

  observed <- gene_table %>%
    dplyr::group_by(Comparison, Analysis) %>%
    dplyr::summarise(
      Std = sum(Std, na.rm = TRUE),
      Strong = sum(Strong, na.rm = TRUE),
      Weak = sum(Weak, na.rm = TRUE),
      HBFSS = sum(HBFSS, na.rm = TRUE),
      Ovlp = sum(Ovlp, na.rm = TRUE),
      TWAS_APA = sum(TWAS_APA, na.rm = TRUE),
      WTTS_APA = sum(WTTS_APA, na.rm = TRUE),
      WTTS_DE_APA = sum(WTTS_DE_APA, na.rm = TRUE),
      .groups = "drop"
    )

  grid %>%
    dplyr::left_join(observed, by = c("Comparison", "Analysis")) %>%
    dplyr::mutate(
      dplyr::across(
        c(Std, Strong, Weak, HBFSS, Ovlp, TWAS_APA, WTTS_APA, WTTS_DE_APA),
        ~ as.integer(dplyr::coalesce(.x, 0))
      ),
      Analysis = factor(Analysis, levels = view_levels)
    ) %>%
    dplyr::arrange(
      match(Comparison, as.character(comparison_table$comparison_name)),
      Analysis
    ) %>%
    dplyr::mutate(Analysis = as.character(Analysis))
}


plot_twas_comparison_gene_panel <- function(gene_table, comparison_name) {
  analysis_map <- c(
    "Original (No EVS)" = "Orig",
    "NormEVS Lead" = "N-Lead",
    "NormEVS Rem" = "N-Rem",
    "RawEVS Lead" = "R-Lead",
    "RawEVS Rem" = "R-Rem"
  )
  analysis_levels <- unname(analysis_map)

  sub <- gene_table[gene_table$Comparison == comparison_name, , drop = FALSE]

  if (!nrow(sub)) {
    empty <- data.frame(
      Analysis = factor(analysis_levels, levels = analysis_levels),
      y = 0
    )
    return(
      ggplot(empty, aes(Analysis, y)) +
        geom_blank() +
        annotate(
          "text",
          x = 3,
          y = 0,
          label = "No 3'aTWAS/WTTS significant-gene overlap",
          size = 4
        ) +
        scale_x_discrete(drop = FALSE) +
        scale_y_continuous(NULL, breaks = NULL) +
        labs(
          title = paste0(comparison_name, " | 3'aTWAS ortholog overlap"),
          x = NULL,
          y = NULL
        ) +
        manuscript_theme() +
        theme(legend.position = "none")
    )
  }

  method_rows <- list()
  for (method in c("Std", "Strong", "Weak", "HBFSS")) {
    keep <- !is.na(sub[[method]]) & sub[[method]]
    if (!any(keep)) next
    tmp <- sub[keep, , drop = FALSE]
    tmp$Method <- c(
      Std = "Standard",
      Strong = "Strong",
      Weak = "Weak",
      HBFSS = "HBFSS"
    )[[method]]
    method_rows[[method]] <- tmp
  }

  if (!length(method_rows)) return(NULL)

  long <- dplyr::bind_rows(method_rows) %>%
    dplyr::mutate(
      Analysis = factor(unname(analysis_map[Analysis]), levels = analysis_levels),
      Method = factor(Method, levels = significance_method_levels),
      Gene_label = ifelse(
        toupper(Rat_Ortholog) == toupper(Human_TWAS),
        Rat_Ortholog,
        paste0(Rat_Ortholog, " [", Human_TWAS, "]")
      )
    )

  gene_order <- sort(unique(long$Gene_label))
  long$Gene_label <- factor(long$Gene_label, levels = rev(gene_order))
  n_gene <- length(gene_order)

  observed_counts <- sub %>%
    dplyr::group_by(Analysis) %>%
    dplyr::summarise(
      Std = sum(Std, na.rm = TRUE),
      Strong = sum(Strong, na.rm = TRUE),
      Weak = sum(Weak, na.rm = TRUE),
      HBFSS = sum(HBFSS, na.rm = TRUE),
      Ovlp = sum(Ovlp, na.rm = TRUE),
      .groups = "drop"
    )

  count_grid <- data.frame(
    Analysis = names(analysis_map),
    Short = unname(analysis_map),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(observed_counts, by = "Analysis") %>%
    dplyr::mutate(
      dplyr::across(
        c(Std, Strong, Weak, HBFSS, Ovlp),
        ~ as.integer(dplyr::coalesce(.x, 0))
      ),
      Axis_label = paste0(
        Short,
        "\nStd=", Std,
        " Str=", Strong,
        " W=", Weak,
        " H=", HBFSS,
        " O=", Ovlp
      )
    )

  x_labels <- stats::setNames(count_grid$Axis_label, count_grid$Short)

  ggplot(long, aes(Analysis, Gene_label, color = Method, shape = Method)) +
    geom_point(
      position = position_dodge(width = 0.48),
      size = 3.0,
      alpha = 0.98,
      stroke = 0.90
    ) +
    scale_x_discrete(drop = FALSE, labels = x_labels) +
    scale_color_manual(
      values = significance_method_colors,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method",
      guide = guide_legend(override.aes = list(size = 3.3, stroke = 1.0))
    ) +
    scale_shape_manual(
      values = significance_method_shapes,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method"
    ) +
    labs(
      title = paste0(comparison_name, " | 3'aTWAS ortholog overlap"),
      subtitle = paste0(
        "Significant WTTS ortholog genes across five analysis views",
        " | genes=", dplyr::n_distinct(sub$Rat_Ortholog)
      ),
      x = NULL,
      y = "Rat ortholog"
    ) +
    manuscript_theme() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      axis.text.y = element_text(size = if (n_gene > 30) 6.2 else 7.5),
      legend.position = "bottom",
      plot.margin = margin(10, 18, 10, 18)
    )
}


export_twas_by_view <- function(pas_table, gene_table, summary_table) {
  views <- comparison_analysis_views()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    pas_cmp <- pas_table[pas_table$Comparison == comparison_name, , drop = FALSE]
    gene_cmp <- gene_table[gene_table$Comparison == comparison_name, , drop = FALSE]
    summary_cmp <- summary_table[summary_table$Comparison == comparison_name, , drop = FALSE]

    save_csv(
      pas_cmp,
      file.path(
        output_dir,
        comparison_name,
        paste0("Table_", comparison_name, "_TWAS_PAS_All_Views.csv")
      )
    )
    save_csv(
      gene_cmp,
      file.path(
        output_dir,
        comparison_name,
        paste0("Table_", comparison_name, "_TWAS_Genes_All_Views.csv")
      )
    )
    save_csv(
      summary_cmp,
      file.path(
        output_dir,
        comparison_name,
        paste0("Table_", comparison_name, "_TWAS_Method_Counts_All_Views.csv")
      )
    )

    for (i in seq_len(nrow(views))) {
      track_key <- views$track_key[i]
      dataset_key <- views$dataset_key[i]
      analysis_label <- analysis_view_label(track_key, dataset_key)
      paths <- analysis_view_paths(comparison_name, track_key, dataset_key)

      pas_sub <- pas_cmp[pas_cmp$Analysis == analysis_label, , drop = FALSE]
      gene_sub <- gene_cmp[gene_cmp$Analysis == analysis_label, , drop = FALSE]

      save_csv(pas_sub, file.path(paths$tables, "Table_TWAS_PAS.csv"))
      save_csv(gene_sub, file.path(paths$tables, "Table_TWAS_Genes.csv"))
    }
  }

  invisible(TRUE)
}


save_twas_comparison_panels <- function(gene_table) {
  for (comparison_name in as.character(comparison_table$comparison_name)) {
    p <- plot_twas_comparison_gene_panel(gene_table, comparison_name)
    if (is.null(p)) next

    n_gene <- dplyr::n_distinct(
      gene_table$Rat_Ortholog[gene_table$Comparison == comparison_name]
    )
    height <- min(14.0, max(5.6, 4.5 + 0.23 * n_gene))

    panel_dir <- file.path(output_dir, comparison_name, "Panels")
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

    base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_TWAS_5Views"))
    save_grob(p, paste0(base, ".png"), width = 12.5, height = height)
    save_grob(p, paste0(base, ".pdf"), width = 12.5, height = height)
  }
  invisible(TRUE)
}



audit_published_2024_twas_reference <- function(all_results) {
  ref <- PUBLISHED_2024_TWAS_DE_REFERENCE

  obs <- all_results %>%
    dplyr::filter(
      Comparison == "RT4_ZT10",
      Analysis == "RawEVS Rem",
      as.character(PAS) %in% ref$PAS
    ) %>%
    dplyr::select(
      PAS,
      WTTS_Gene,
      Apeglm_LFC,
      Std_p,
      Std_BH,
      Std
    )

  merged <- ref %>%
    dplyr::left_join(obs, by = "PAS")

  if (nrow(merged) != nrow(ref) ||
      any(is.na(merged$WTTS_Gene)) ||
      dplyr::n_distinct(merged$Gene) != 9L ||
      !all(!is.na(merged$Std) & merged$Std)) {
    stop(
      "Published 2024 WTTS/TWAS benchmark was not reproduced: ",
      "expected 11 Standard PASs representing 9 genes in RT4_ZT10 RawEVS Rem."
    )
  }

  gene_ok <- toupper(merged$Gene) == toupper(merged$WTTS_Gene)
  lfc_ok <- abs(merged$Apeglm_LFC - merged$Paper_LFC) <= 0.02

  p_ok <- abs(merged$Std_p - merged$Paper_p) <= pmax(
    5e-8,
    abs(merged$Paper_p) * 0.05
  )

  padj_ok <- abs(merged$Std_BH - merged$Paper_padj) <= pmax(
    5e-6,
    abs(merged$Paper_padj) * 0.05
  )

  if (!all(gene_ok & lfc_ok & p_ok & padj_ok)) {
    bad <- merged$PAS[!(gene_ok & lfc_ok & p_ok & padj_ok)]
    stop(
      "Published 2024 WTTS/TWAS numeric benchmark mismatch for PAS(s): ",
      paste(bad, collapse = ", ")
    )
  }

  message(
    "Published 2024 WTTS/TWAS benchmark reproduced: ",
    "11 PASs / 9 genes in RT4_ZT10 RawEVS Rem."
  )

  invisible(TRUE)
}


audit_final_exports <- function(overall_summary, twas_results) {
  expected_views <- c(
    "Original (No EVS)",
    "NormEVS Lead",
    "NormEVS Rem",
    "RawEVS Lead",
    "RawEVS Rem"
  )

  if (!is.data.frame(overall_summary) ||
      nrow(overall_summary) != nrow(comparison_table) * length(expected_views)) {
    stop("Final export audit failed: differential-expression summary does not contain all comparison/view rows.")
  }

  twas_summary <- twas_results$summary
  if (!is.data.frame(twas_summary) ||
      nrow(twas_summary) != nrow(comparison_table) * length(expected_views)) {
    stop("Final export audit failed: TWAS summary does not contain all five views for every comparison.")
  }

  twas_pas <- twas_results$pas
  if (nrow(twas_pas)) {
    twas_pas_key <- paste(
      twas_pas$Comparison,
      twas_pas$Analysis,
      twas_pas$Rat_Ortholog,
      twas_pas$WTTS_DE_PAS,
      sep = "||"
    )
    if (anyDuplicated(twas_pas_key)) {
      stop("Final export audit failed: duplicate TWAS PAS rows were detected.")
    }
  }

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    panel_dir <- file.path(output_dir, comparison_name, "Panels")

    required_panels <- file.path(
      panel_dir,
      c(
        paste0("Figure_", comparison_name, "_PCA_EVS.png"),
        paste0("Figure_", comparison_name, "_Volcano_5Views.png"),
        paste0("Figure_", comparison_name, "_TWAS_5Views.png")
      )
    )

    missing_panels <- required_panels[
      !file.exists(required_panels) |
        is.na(file.info(required_panels)$size) |
        file.info(required_panels)$size <= 0
    ]

    if (length(missing_panels)) {
      stop(
        "Final export audit failed: missing/empty panel(s) for ",
        comparison_name, ": ",
        paste(basename(missing_panels), collapse = ", ")
      )
    }

    required_tables <- c(
      file.path(
        output_dir,
        comparison_name,
        paste0("Table_", comparison_name, "_Significant_Sites_All_Views.csv")
      ),
      file.path(
        output_dir,
        comparison_name,
        "Table_Method_Counts_All_Views.csv"
      ),
      file.path(
        output_dir,
        comparison_name,
        paste0("Table_", comparison_name, "_TWAS_PAS_All_Views.csv")
      ),
      file.path(
        output_dir,
        comparison_name,
        paste0("Table_", comparison_name, "_TWAS_Genes_All_Views.csv")
      ),
      file.path(
        output_dir,
        comparison_name,
        paste0("Table_", comparison_name, "_TWAS_Method_Counts_All_Views.csv")
      )
    )

    missing_tables <- required_tables[!file.exists(required_tables)]
    if (length(missing_tables)) {
      stop(
        "Final export audit failed: missing table(s) for ",
        comparison_name, ": ",
        paste(basename(missing_tables), collapse = ", ")
      )
    }

    sig_tbl <- utils::read.csv(
      file.path(
        output_dir,
        comparison_name,
        paste0("Table_", comparison_name, "_Significant_Sites_All_Views.csv")
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    count_tbl <- utils::read.csv(
      file.path(
        output_dir,
        comparison_name,
        "Table_Method_Counts_All_Views.csv"
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (!setequal(count_tbl$Analysis, expected_views)) {
      stop("Final export audit failed: incomplete DE method-count views for ", comparison_name, ".")
    }

    for (analysis_name in expected_views) {
      one_count <- count_tbl[count_tbl$Analysis == analysis_name, , drop = FALSE]
      one_sig <- sig_tbl[sig_tbl$Analysis == analysis_name, , drop = FALSE]

      if (nrow(one_count) != 1L) {
        stop("Final export audit failed: duplicate/missing DE count row for ", comparison_name, " / ", analysis_name, ".")
      }

      observed <- c(
        Std = sum(one_sig$Std, na.rm = TRUE),
        Strong = sum(one_sig$Strong, na.rm = TRUE),
        Weak = sum(one_sig$Weak, na.rm = TRUE),
        HBFSS = sum(one_sig$HBFSS_sig, na.rm = TRUE),
        Ovlp = sum(one_sig$Ovlp, na.rm = TRUE)
      )

      reported <- as.numeric(one_count[1, c("Std", "Strong", "Weak", "HBFSS", "Ovlp")])
      names(reported) <- c("Std", "Strong", "Weak", "HBFSS", "Ovlp")

      if (!identical(as.integer(observed), as.integer(reported))) {
        stop(
          "Final export audit failed: significant-site counts disagree for ",
          comparison_name, " / ", analysis_name, "."
        )
      }
    }

    twas_count_path <- file.path(
      output_dir,
      comparison_name,
      paste0("Table_", comparison_name, "_TWAS_Method_Counts_All_Views.csv")
    )
    twas_count_tbl <- utils::read.csv(
      twas_count_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (!setequal(twas_count_tbl$Analysis, expected_views) ||
        nrow(twas_count_tbl) != length(expected_views)) {
      stop("Final export audit failed: incomplete TWAS views for ", comparison_name, ".")
    }
  }

  discovery_panel <- file.path(
    paper_fig_dir,
    "Figure_Manuscript_Discovery_Counts.png"
  )
  if (!file.exists(discovery_panel) ||
      is.na(file.info(discovery_panel)$size) ||
      file.info(discovery_panel)$size <= 0) {
    stop("Final export audit failed: global discovery-count panel is missing or empty.")
  }

  message(
    "Final export audit passed: all comparison/view tables, curated panels, ",
    "method counts, TWAS rows, and published benchmark are internally consistent."
  )

  invisible(TRUE)
}

run_twas_overlap_analysis <- function() {
  twas_file <- resolve_existing_file(twas_file_candidates, "3'aTWAS file")
  message("Running final 3'aTWAS ortholog overlap: ", twas_file)

  twas <- read_twas_study(twas_file)
  mapping <- build_twas_ortholog_map(twas)
  rat_meta <- build_twas_rat_metadata(twas, mapping)
  all_results <- collect_twas_analysis_results()

  audit_published_2024_twas_reference(all_results)

  pas_table <- build_twas_overlap_pas_table(all_results, rat_meta)
  gene_table <- build_twas_overlap_gene_table(pas_table)
  summary_table <- build_twas_overlap_summary(gene_table)

  twas_dir <- file.path(output_dir, "TWAS")
  twas_tab_dir <- file.path(twas_dir, "tables")
  dir.create(twas_tab_dir, recursive = TRUE, showWarnings = FALSE)

  save_csv(pas_table, file.path(twas_tab_dir, "Table_TWAS_Overlap_PAS.csv"))
  save_csv(gene_table, file.path(twas_tab_dir, "Table_TWAS_Overlap_Genes.csv"))
  save_csv(summary_table, file.path(twas_tab_dir, "Table_TWAS_Method_Counts.csv"))

  export_twas_by_view(pas_table, gene_table, summary_table)
  save_twas_comparison_panels(gene_table)

  invisible(
    list(
      mapping = mapping,
      pas = pas_table,
      genes = gene_table,
      summary = summary_table
    )
  )
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

  if (is.null(out) && is.null(failed_comparisons[[cmp]])) {
    failed_comparisons[[cmp]] <- data.frame(
      comparison_name = cmp,
      error_message = "Comparison did not complete.",
      stringsAsFactors = FALSE
    )
  }
}

if (length(failed_comparisons) > 0L) {
  failed_df <- dplyr::bind_rows(failed_comparisons)
  writeLines(
    apply(failed_df, 1, function(x) paste(x, collapse = " | ")),
    file.path(output_dir, "Failed_Comparisons.txt")
  )
  stop("One or more comparisons failed. See Failed_Comparisons.txt.")
}

overall_summary <- build_overall_manuscript_summary()
if (nrow(overall_summary) > 0L) {
  save_csv(overall_summary, file.path(summary_table_dir, "Table_DE_Method_Counts.csv"))
}

write_methods_note()
save_paper_volcano_panels()
save_paper_support_figures(overall_summary)

# Final analytical stage: 3'aTWAS ortholog overlap with every reported WTTS method.
twas_results <- run_twas_overlap_analysis()

audit_final_exports(overall_summary, twas_results)

figure_zip <- create_all_figures_zip()
table_zip <- create_all_tables_zip()

# Publish only after every analysis, export, audit, and archive step succeeds.
repo_publish <- publish_results_to_repository(figure_zip, table_zip)

cat("\n=====================================================\n")
cat("Pipeline complete.\n")
cat("Build: ", PIPELINE_BUILD, "\n", sep = "")
cat("Repository root:\n", repo_root, "\n", sep = "")
cat("Output directory:\n", output_dir, "\n", sep = "")
cat("3'aTWAS overlap: complete\n")
cat("Figure ZIP:\n", figure_zip, "\n", sep = "")
cat("Table ZIP:\n", table_zip, "\n", sep = "")
if (isTRUE(repo_publish$published)) {
  cat("Git branch:\n", repo_publish$branch, "\n", sep = "")
  cat("Git commit:\n", repo_publish$commit, "\n", sep = "")
  if (!is.na(repo_publish$url) && nzchar(repo_publish$url)) {
    cat("GitHub output folder:\n", repo_publish$url, "\n", sep = "")
  } else {
    cat("Repository output path:\n", repo_publish$output_rel, "\n", sep = "")
  }
}
cat("=====================================================\n\n")

if (nrow(overall_summary) > 0L) print(overall_summary)
