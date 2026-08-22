#!/usr/bin/env Rscript

# =============================================================================
# SEQUENCE MANUSCRIPT ANALYSIS
# WTTS-Seq EVS + DESeq2 + apeglm + empirical-null HBFSS + 3'aTWAS overlap
# =============================================================================

options(stringsAsFactors = FALSE)

# =============================================================================
# 1. PACKAGES
# =============================================================================

required_packages <- c(
  "DESeq2",
  "apeglm",
  "fdrtool",
  "ggplot2",
  "ggrepel",
  "dplyr",
  "tidyr",
  "babelgene"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install required R package(s) before running: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
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
  library(babelgene)
})

# =============================================================================
# 2. USER SETTINGS
# =============================================================================

COUNT_FILE_CANDIDATES <- c(
  "WTTS-Seq_2022.2_DE_raw_read_numbers(20260822-183312).csv",
  "WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  file.path("data", "WTTS-Seq_2022.2_DE_raw_read_numbers(20260822-183312).csv"),
  file.path("data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"),
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  "/mnt/data/WTTS-Seq_2022.2_DE_raw_read_numbers(20260822-183312).csv"
)

TWAS_FILE_CANDIDATES <- c(
  "3aTWAS_genes_of_11_brain_disorders.csv",
  file.path("data", "3aTWAS_genes_of_11_brain_disorders.csv"),
  "/root/REAPER98632/data/3aTWAS_genes_of_11_brain_disorders.csv",
  "/mnt/data/3aTWAS_genes_of_11_brain_disorders.csv"
)

OUTPUT_DIR <- file.path("exports", "SEQUENCE_manuscript_complete")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
TABLE_DIR  <- file.path(OUTPUT_DIR, "tables")
VALIDATION_DIR <- file.path(OUTPUT_DIR, "validation")

DESEQ2_FDR <- 0.10
LFC_BOUNDARY <- 1.0
EVS_TOP_N <- 5000L
VOLCANO_LABEL_N <- 20L
FIGURE_DPI <- 320L
EXPORT_PDF <- TRUE

# Validation-module controls.
RUN_EVS_CUTOFF_SENSITIVITY <- FALSE
EMPIRICAL_EVS_TOP_N <- 4077L

RUN_SYNTHETIC_VALIDATION <- FALSE
SIMULATION_SEED <- 42L
SIMULATION_REPLICATES <- 25L
SIMULATION_FEATURES <- 5000L
SIMULATION_SAMPLES_PER_GROUP <- 5L
SIMULATION_DE_FRACTION <- 0.10
SIMULATION_STRONG_FRACTION_WITHIN_DE <- 0.50
SIMULATION_WEAK_LFC_RANGE <- c(0.20, 0.80)
SIMULATION_STRONG_LFC_RANGE <- c(1.20, 2.00)

# Human 3'aTWAS symbols are mapped to rat orthologs because the WTTS symbols are
# Rattus norvegicus gene symbols.
TWAS_TARGET_SPECIES <- "rat"
TWAS_ORTHOLOG_MIN_SUPPORT <- 1L

# =============================================================================
# 3. STUDY METADATA
# =============================================================================

sample_metadata <- data.frame(
  sample_id = c(
    paste0("R0_", 1:5), paste0("ZT6_", 1:5),
    paste0("R2_", 1:5), paste0("ZT8_", 1:5),
    paste0("R4_", 1:5), paste0("ZT10_", 1:5),
    paste0("R8_", 1:5), paste0("ZT14_", 1:5)
  ),
  condition = c(
    rep("treatment", 5), rep("control", 5),
    rep("treatment", 5), rep("control", 5),
    rep("treatment", 5), rep("control", 5),
    rep("treatment", 5), rep("control", 5)
  ),
  stringsAsFactors = FALSE
)
rownames(sample_metadata) <- sample_metadata$sample_id
sample_metadata$condition <- factor(
  sample_metadata$condition,
  levels = c("control", "treatment")
)

comparison_table <- data.frame(
  comparison = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  treatment_prefix = c("R0", "R2", "R4", "R8"),
  control_prefix = c("ZT6", "ZT8", "ZT10", "ZT14"),
  treatment_label = c("RT0", "RT2", "RT4", "RT8"),
  control_label = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

analysis_view_levels <- c(
  "Original",
  "NormEVS Lead",
  "NormEVS Rem",
  "RawEVS Lead",
  "RawEVS Rem"
)

# =============================================================================
# 4. GENERAL HELPERS
# =============================================================================

resolve_existing_file <- function(candidates, label) {
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L) {
    stop(
      label, " was not found. Checked:\n  ",
      paste(candidates, collapse = "\n  "),
      call. = FALSE
    )
  }
  normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
}

ensure_output_dirs <- function() {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(VALIDATION_DIR, recursive = TRUE, showWarnings = FALSE)
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

save_figure <- function(plot, stem, width, height, directory = FIGURE_DIR) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(directory, paste0(stem, ".png"))
  ggplot2::ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = FIGURE_DPI,
    bg = "white",
    limitsize = FALSE
  )
  if (isTRUE(EXPORT_PDF)) {
    pdf_path <- file.path(directory, paste0(stem, ".pdf"))
    ggplot2::ggsave(
      filename = pdf_path,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      bg = "white",
      limitsize = FALSE
    )
  }
  invisible(png_path)
}


# Manuscript-ready description of every statistical method executed by the pipeline.
write_manuscript_methods <- function() {
  methods <- c(
    "# Manuscript statistical methods",
    "",
    "## Study unit and comparisons",
    "WTTS-Seq was analyzed at the individual polyadenylation-site (PAS) level using OrigID as the feature identifier and Gene Symbol as annotation. Four pairwise comparisons were evaluated: RT0 versus ZT6, RT2 versus ZT8, RT4 versus ZT10, and RT8 versus ZT14. Positive log2 fold change indicates higher abundance in the RT condition.",
    "",
    "## Eigenvector splitting",
    paste0("For each comparison, NormEVS normalized the complete raw-count matrix with DESeq2 median-of-ratios size factors before eigenvector splitting, whereas RawEVS used raw counts directly. Within each track, RT and ZT samples were decomposed separately with prcomp(center=TRUE, scale.=FALSE). PASs were ranked by absolute PC1 feature loading. Exactly the top ", EVS_TOP_N, " PASs from the RT eigenvector and the top ", EVS_TOP_N, " PASs from the ZT eigenvector were selected. Their intersection defined Joint PASs; PASs unique to either condition defined Disjoint PASs; Joint plus both Disjoint sets formed the Leading Edge; all other PASs formed the Remainder. EVS selected feature IDs only, and downstream DESeq2 models were fit to the corresponding raw-count subsets."),
    "",
    "## Differential testing and effect-size shrinkage",
    paste0("Each Original, NormEVS Lead, NormEVS Rem, RawEVS Lead, and RawEVS Rem dataset was analyzed independently with DESeq2 using design ~ condition and ZT/control as the reference. Standard DESeq2 significance required Benjamini-Hochberg adjusted p < ", DESEQ2_FDR, " and |apeglm-shrunken log2 fold change| >= ", LFC_BOUNDARY, ". Strong composite-null testing used lfcThreshold=", LFC_BOUNDARY, " and altHypothesis='greaterAbs' with BH adjusted p < ", DESEQ2_FDR, ". Weak composite-null testing used the same LFC threshold with altHypothesis='lessAbs' and BH adjusted p < ", DESEQ2_FDR, ". Apeglm-shrunken log2 fold changes were used for reported LFC values and HBFSS calculation."),
    "",
    "## Empirical null, higher criticism, and HBFSS",
    paste0("Finite DESeq2 Wald statistics were calibrated with Strimmer's fdrtool using statistic='normal'. The resulting empirical-null p-values were supplied to higher criticism as HCp = hc.thresh(sort(empirical_p)). The dataset-specific HBFSS cutoff was Htau = -log10(HCp) x ", LFC_BOUNDARY, ". For PAS i, HBFSS_i = |apeglm_LFC_i| x [-log10(empirical_p_i)]. HBFSS significance required HBFSS_i > Htau. HCp was used to derive Htau; HBFSS significance was determined by the score threshold HBFSS_i > Htau."),
    "",
    "## Weak composite-null support and final discoveries",
    "Weak composite-null PASs were considered final weak-effect discoveries only when they also exceeded the dataset-specific HBFSS threshold. Final PAS tables retained HBFSS-significant PASs, Standard-DESeq2-significant PASs, and Strong-CNH-significant PASs; weak-CNH-only PASs were not treated as final discoveries.",
    "",
    "## Volcano visualization",
    paste0("Volcano plots used apeglm-shrunken LFC on the x-axis and -log10(empirical p) on the y-axis. The HBFSS boundary y = Htau/|LFC|, the HCp reference, and LFC boundaries at +/-", LFC_BOUNDARY, " were displayed. HBFSS-positive points shared one color, with marker shape encoding HBFSS-only, HBFSS+Standard, HBFSS+Strong, or HBFSS+Weak support. A maximum of ", VOLCANO_LABEL_N, " final significant PASs were labeled per analysis panel."),
    "",
    "## 3'aTWAS ortholog overlap",
    paste0("Human gene symbols in the supplied 3'aTWAS study were mapped to Rattus norvegicus ortholog symbols using babelgene::orthologs with target species='", TWAS_TARGET_SPECIES, "', min_support=", TWAS_ORTHOLOG_MIN_SUPPORT, ", and top=TRUE so the highest-supported ortholog mapping was carried into the WTTS comparison. For every comparison and each of the five analysis views, TWAS orthologs were intersected separately with Standard DESeq2 and HBFSS significant WTTS genes. The TWAS overlap tables report the human TWAS gene(s), rat ortholog, TWAS disease annotations, method-specific significant PAS counts, and the WTTS PAS multiplicity for that gene. APA_multi_PAS denotes WTTS genes represented by at least two distinct PAS features; DE_multi_PAS denotes TWAS-overlap genes with at least two significant PASs in the same comparison/view."),
    "",
    "## Optional validation modules",
    paste0("EVS cutoff sensitivity enabled: ", RUN_EVS_CUTOFF_SENSITIVITY, ". When enabled, the normalized EVS leading edge is recomputed and analyzed independently at K=", EVS_TOP_N, " and K=", EMPIRICAL_EVS_TOP_N, "."),
    paste0("Fully artificial negative-binomial validation enabled: ", RUN_SYNTHETIC_VALIDATION, ". When enabled, it uses ", SIMULATION_FEATURES, " PASs, ", SIMULATION_SAMPLES_PER_GROUP, " samples per group, ", SIMULATION_REPLICATES, " replicates, known null/weak/strong truth labels, and the same DESeq2, apeglm, empirical-null, higher-criticism, HBFSS, Standard, Strong, and Weak decision functions as the manuscript analysis.")
  )

  writeLines(methods, file.path(OUTPUT_DIR, "METHODS_MANUSCRIPT.md"))
}



assert_columns <- function(df, cols, label) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0L) {
    stop(
      label, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

coerce_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)), fixed = TRUE)))
}

coerce_count_matrix <- function(x, label) {
  rn <- rownames(x)
  cn <- colnames(x)
  mat <- as.matrix(x)
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (any(!is.finite(mat) | is.na(mat))) {
    stop(label, " contains non-finite or missing counts.", call. = FALSE)
  }
  if (any(mat < 0)) {
    stop(label, " contains negative counts.", call. = FALSE)
  }
  rounded <- round(mat)
  if (any(abs(mat - rounded) > 1e-6)) {
    warning(label, " contained non-integer values; values were rounded for DESeq2.")
  }
  if (any(rounded > .Machine$integer.max)) {
    stop(label, " contains values above R integer range.", call. = FALSE)
  }
  storage.mode(rounded) <- "integer"
  rownames(rounded) <- rn
  colnames(rounded) <- cn
  rounded
}

valid_gene_symbol <- function(x) {
  y <- trimws(as.character(x))
  !is.na(y) & nzchar(y) & y != "-"
}

gene_key <- function(x) {
  y <- trimws(as.character(x))
  y[!valid_gene_symbol(y)] <- NA_character_
  toupper(y)
}

safe_neglog10 <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p) & !is.na(p) & p >= 0 & p <= 1
  out[ok] <- -log10(pmax(p[ok], .Machine$double.xmin))
  out
}

safe_min <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

collapse_unique <- function(x) {
  x <- sort(unique(trimws(as.character(x))))
  x <- x[!is.na(x) & nzchar(x)]
  paste(x, collapse = ";")
}

format_threshold <- function(x) {
  ifelse(is.finite(x) & !is.na(x), formatC(x, digits = 3, format = "e"), "NA")
}

manuscript_theme <- function() {
  ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.20, color = "grey90"),
      strip.background = ggplot2::element_rect(fill = "grey96", color = "grey75"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9),
      legend.position = "bottom"
    )
}

# =============================================================================
# 5. INPUT WTTS DATA AND PAS ANNOTATION
# =============================================================================

ensure_output_dirs()
write_manuscript_methods()

count_file <- resolve_existing_file(COUNT_FILE_CANDIDATES, "WTTS count file")
message("WTTS count file: ", count_file)

wtts <- utils::read.csv(
  count_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

assert_columns(wtts, c("OrigID", "Symbol"), "WTTS count file")
assert_columns(wtts, sample_metadata$sample_id, "WTTS count file")

for (sid in sample_metadata$sample_id) {
  wtts[[sid]] <- coerce_numeric(wtts[[sid]])
}

wtts <- wtts[
  !is.na(wtts$OrigID) & nzchar(trimws(as.character(wtts$OrigID))),
  ,
  drop = FALSE
]

bad_count_row <- rowSums(is.na(wtts[, sample_metadata$sample_id, drop = FALSE])) > 0
if (any(bad_count_row)) {
  stop(
    "WTTS input contains ", sum(bad_count_row),
    " row(s) with missing/non-numeric sample counts.",
    call. = FALSE
  )
}

wtts$feature_id <- make.unique(as.character(wtts$OrigID), sep = "_dup")
rownames(wtts) <- wtts$feature_id

pas_annotation <- data.frame(
  feature_id = wtts$feature_id,
  PAS_ID = as.character(wtts$OrigID),
  Gene = trimws(as.character(wtts$Symbol)),
  Chromosome = if ("Chromosome" %in% names(wtts)) as.character(wtts$Chromosome) else NA_character_,
  Strand = if ("Strand" %in% names(wtts)) as.character(wtts$Strand) else NA_character_,
  Peak = if ("Peak" %in% names(wtts)) as.character(wtts$Peak) else NA_character_,
  stringsAsFactors = FALSE
)
pas_annotation$Gene[!valid_gene_symbol(pas_annotation$Gene)] <- NA_character_
pas_annotation$gene_key <- gene_key(pas_annotation$Gene)

wtts_gene_pas_summary <- pas_annotation %>%
  dplyr::filter(!is.na(gene_key)) %>%
  dplyr::group_by(gene_key) %>%
  dplyr::summarise(
    WTTS_Gene = dplyr::first(Gene),
    WTTS_PAS_n = dplyr::n_distinct(PAS_ID),
    .groups = "drop"
  ) %>%
  dplyr::mutate(APA_multi_PAS = WTTS_PAS_n >= 2L)

# =============================================================================
# 6. COMPARISON PREPARATION
# =============================================================================

prepare_comparison <- function(comparison_row) {
  trt_pattern <- paste0("^", comparison_row$treatment_prefix, "_")
  ctl_pattern <- paste0("^", comparison_row$control_prefix, "_")

  sample_ids <- sample_metadata$sample_id[
    grepl(trt_pattern, sample_metadata$sample_id) |
      grepl(ctl_pattern, sample_metadata$sample_id)
  ]

  coldata <- sample_metadata[sample_ids, "condition", drop = FALSE]
  count_mat <- wtts[, sample_ids, drop = FALSE]
  rownames(count_mat) <- wtts$feature_id

  list(
    comparison = comparison_row$comparison,
    treatment_label = comparison_row$treatment_label,
    control_label = comparison_row$control_label,
    counts = coerce_count_matrix(count_mat, paste0(comparison_row$comparison, " counts")),
    coldata = coldata
  )
}

# =============================================================================
# 7. EIGENVECTOR SPLITTING
# =============================================================================

normalized_matrix_for_evs <- function(raw_counts, coldata) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = coldata,
    design = ~ condition
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  as.matrix(DESeq2::counts(dds, normalized = TRUE))
}

pc1_loading_rank <- function(matrix_for_evs, sample_ids, top_n) {
  x <- as.matrix(matrix_for_evs[, sample_ids, drop = FALSE])
  storage.mode(x) <- "numeric"

  if (nrow(x) < top_n) {
    stop(
      "EVS requires at least ", top_n,
      " PAS features but only ", nrow(x), " are available.",
      call. = FALSE
    )
  }
  if (ncol(x) < 2L) {
    stop("Condition-specific PCA requires at least two samples.", call. = FALSE)
  }

  pca <- stats::prcomp(t(x), center = TRUE, scale. = FALSE)
  load <- abs(pca$rotation[, 1])

  tbl <- data.frame(
    feature_id = names(load),
    PC1_abs_loading = unname(load),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(dplyr::desc(PC1_abs_loading), feature_id) %>%
    dplyr::mutate(
      Rank = dplyr::row_number(),
      Selected = Rank <= top_n
    )

  list(
    pca = pca,
    table = tbl,
    top_ids = tbl$feature_id[tbl$Selected],
    loading_at_rank = tbl$PC1_abs_loading[top_n]
  )
}

build_evs <- function(raw_counts, coldata, track, treatment_label, control_label, top_n = EVS_TOP_N) {
  track <- match.arg(track, c("NormEVS", "RawEVS"))

  evs_matrix <- if (track == "NormEVS") {
    normalized_matrix_for_evs(raw_counts, coldata)
  } else {
    raw_counts
  }

  trt_samples <- rownames(coldata)[coldata$condition == "treatment"]
  ctl_samples <- rownames(coldata)[coldata$condition == "control"]

  trt <- pc1_loading_rank(evs_matrix, trt_samples, top_n)
  ctl <- pc1_loading_rank(evs_matrix, ctl_samples, top_n)

  joint <- intersect(trt$top_ids, ctl$top_ids)
  disjoint_trt <- setdiff(trt$top_ids, ctl$top_ids)
  disjoint_ctl <- setdiff(ctl$top_ids, trt$top_ids)
  leading <- union(trt$top_ids, ctl$top_ids)
  remainder <- setdiff(rownames(raw_counts), leading)

  membership <- data.frame(
    feature_id = rownames(raw_counts),
    EVS_Class = "Remainder",
    stringsAsFactors = FALSE
  )
  membership$EVS_Class[membership$feature_id %in% joint] <- "Joint"
  membership$EVS_Class[membership$feature_id %in% disjoint_trt] <- paste0("Disjoint_", treatment_label)
  membership$EVS_Class[membership$feature_id %in% disjoint_ctl] <- paste0("Disjoint_", control_label)

  rank_support <- dplyr::bind_rows(
    trt$table %>%
      dplyr::mutate(Track = track, Condition = treatment_label),
    ctl$table %>%
      dplyr::mutate(Track = track, Condition = control_label)
  )

  list(
    track = track,
    top_n = top_n,
    treatment_top_ids = trt$top_ids,
    control_top_ids = ctl$top_ids,
    joint_ids = joint,
    disjoint_treatment_ids = disjoint_trt,
    disjoint_control_ids = disjoint_ctl,
    leading_ids = leading,
    remainder_ids = remainder,
    membership = membership,
    rank_support = rank_support,
    treatment_loading_at_rank = trt$loading_at_rank,
    control_loading_at_rank = ctl$loading_at_rank,
    leading_counts = raw_counts[leading, , drop = FALSE],
    remainder_counts = raw_counts[remainder, , drop = FALSE]
  )
}

# =============================================================================
# 8. DESEQ2, APEGLM, EMPIRICAL NULL, HIGHER CRITICISM, AND HBFSS
# =============================================================================

get_condition_coef <- function(dds) {
  rn <- DESeq2::resultsNames(dds)
  hit <- grep("^condition_.*_vs_.*$", rn, value = TRUE)
  if (length(hit) != 1L) {
    stop(
      "Could not uniquely identify the condition coefficient. resultsNames: ",
      paste(rn, collapse = ", "),
      call. = FALSE
    )
  }
  hit[[1]]
}

empirical_null_from_wald <- function(wald_stat, dataset_label) {
  valid <- is.finite(wald_stat) & !is.na(wald_stat)
  z <- as.numeric(wald_stat[valid])

  if (length(z) < 5L) {
    stop(dataset_label, ": fewer than five finite Wald statistics for fdrtool.", call. = FALSE)
  }

  fit <- tryCatch(
    fdrtool::fdrtool(
      z,
      statistic = "normal",
      plot = FALSE,
      verbose = FALSE,
      cutoff.method = "fndr"
    ),
    error = function(e) {
      stop(dataset_label, ": fdrtool failed: ", conditionMessage(e), call. = FALSE)
    }
  )

  empirical_p <- rep(NA_real_, length(wald_stat))
  empirical_values <- as.numeric(fit$pval)
  if (length(empirical_values) != length(z)) {
    stop(dataset_label, ": fdrtool empirical-p output length did not match the Wald input.", call. = FALSE)
  }
  empirical_p[valid] <- empirical_values

  if (any(empirical_p[valid] < 0 | empirical_p[valid] > 1, na.rm = TRUE)) {
    stop(dataset_label, ": fdrtool returned invalid empirical p-values.", call. = FALSE)
  }

  hc_input <- sort(empirical_p[is.finite(empirical_p) & !is.na(empirical_p)], na.last = NA)
  hc_input <- pmax(hc_input, .Machine$double.xmin)

  hc_p <- suppressWarnings(
    tryCatch(
      as.numeric(fdrtool::hc.thresh(hc_input))[1],
      error = function(e) {
        stop(dataset_label, ": higher-criticism threshold calculation failed: ", conditionMessage(e), call. = FALSE)
      }
    )
  )

  if (!is.finite(hc_p) || is.na(hc_p) || hc_p <= 0 || hc_p > 1) {
    stop(dataset_label, ": higher criticism did not return a valid probability in (0, 1].", call. = FALSE)
  }

  list(empirical_p = empirical_p, hc_p = hc_p)
}

build_support_label <- function(std, strong, weak, hbfss) {
  out <- character(length(std))
  for (i in seq_along(out)) {
    pieces <- character(0)
    if (isTRUE(hbfss[i])) pieces <- c(pieces, "HBFSS")
    if (isTRUE(std[i])) pieces <- c(pieces, "Std")
    if (isTRUE(strong[i])) pieces <- c(pieces, "Strong")
    if (isTRUE(weak[i] && hbfss[i])) pieces <- c(pieces, "Weak")
    out[i] <- if (length(pieces) == 0L) "None" else paste(pieces, collapse = "+")
  }
  out
}

run_differential_analysis <- function(
    raw_counts,
    coldata,
    comparison,
    analysis_view,
    evs_membership = NULL,
    annotation_df = pas_annotation,
    gene_pas_summary = wtts_gene_pas_summary) {
  raw_counts <- coerce_count_matrix(
    raw_counts,
    paste0(comparison, " | ", analysis_view, " DESeq2 input")
  )

  keep <- rowSums(raw_counts) > 0
  raw_counts <- raw_counts[keep, , drop = FALSE]

  if (nrow(raw_counts) < 5L) {
    stop(comparison, " | ", analysis_view, ": fewer than five nonzero PASs.", call. = FALSE)
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = coldata,
    design = ~ condition
  )
  dds <- DESeq2::DESeq(dds, betaPrior = FALSE, quiet = TRUE)

  std <- DESeq2::results(
    dds,
    contrast = c("condition", "treatment", "control"),
    alpha = DESEQ2_FDR,
    pAdjustMethod = "BH"
  )

  strong <- DESeq2::results(
    dds,
    contrast = c("condition", "treatment", "control"),
    lfcThreshold = LFC_BOUNDARY,
    altHypothesis = "greaterAbs",
    alpha = DESEQ2_FDR,
    pAdjustMethod = "BH"
  )

  weak <- DESeq2::results(
    dds,
    contrast = c("condition", "treatment", "control"),
    lfcThreshold = LFC_BOUNDARY,
    altHypothesis = "lessAbs",
    alpha = DESEQ2_FDR,
    pAdjustMethod = "BH"
  )

  coef_name <- get_condition_coef(dds)
  shr <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")

  std_df <- as.data.frame(std)
  std_df$feature_id <- rownames(std_df)
  strong_df <- as.data.frame(strong)
  strong_df$feature_id <- rownames(strong_df)
  weak_df <- as.data.frame(weak)
  weak_df$feature_id <- rownames(weak_df)
  shr_df <- as.data.frame(shr)
  shr_df$feature_id <- rownames(shr_df)

  result <- data.frame(
    feature_id = std_df$feature_id,
    Wald_stat = std_df$stat,
    DESeq2_padj = std_df$padj,
    Strong_padj = strong_df$padj[match(std_df$feature_id, strong_df$feature_id)],
    Weak_padj = weak_df$padj[match(std_df$feature_id, weak_df$feature_id)],
    Apeglm_LFC = shr_df$log2FoldChange[match(std_df$feature_id, shr_df$feature_id)],
    stringsAsFactors = FALSE
  )

  emp <- empirical_null_from_wald(result$Wald_stat, paste(comparison, analysis_view, sep = " | "))
  result$Empirical_p <- emp$empirical_p

  hc_p <- emp$hc_p
  h_tau <- -log10(hc_p) * LFC_BOUNDARY

  result$HBFSS <- abs(result$Apeglm_LFC) * safe_neglog10(result$Empirical_p)

  result$Std <- !is.na(result$DESeq2_padj) &
    result$DESeq2_padj < DESEQ2_FDR &
    !is.na(result$Apeglm_LFC) &
    abs(result$Apeglm_LFC) >= LFC_BOUNDARY

  result$Strong <- !is.na(result$Strong_padj) &
    result$Strong_padj < DESEQ2_FDR &
    !is.na(result$Apeglm_LFC) &
    abs(result$Apeglm_LFC) >= LFC_BOUNDARY

  result$Weak_CNH <- !is.na(result$Weak_padj) &
    result$Weak_padj < DESEQ2_FDR &
    !is.na(result$Apeglm_LFC) &
    abs(result$Apeglm_LFC) < LFC_BOUNDARY

  result$HBFSS_sig <- !is.na(result$HBFSS) &
    is.finite(result$HBFSS) &
    result$HBFSS > h_tau

  result$Weak <- result$Weak_CNH & result$HBFSS_sig
  result$Overlap <- result$HBFSS_sig & (result$Std | result$Strong | result$Weak_CNH)
  result$Final_sig <- result$HBFSS_sig | result$Std | result$Strong
  result$Direction <- ifelse(
    is.na(result$Apeglm_LFC), NA_character_,
    ifelse(result$Apeglm_LFC > 0, "RT_higher", ifelse(result$Apeglm_LFC < 0, "ZT_higher", "No_change"))
  )

  result$Support <- build_support_label(result$Std, result$Strong, result$Weak_CNH, result$HBFSS_sig)
  result$HCp <- hc_p
  result$Htau <- h_tau
  result$Comparison <- comparison
  result$Analysis <- analysis_view

  result <- result %>%
    dplyr::left_join(annotation_df, by = "feature_id") %>%
    dplyr::left_join(gene_pas_summary[, c("gene_key", "WTTS_PAS_n", "APA_multi_PAS")], by = "gene_key")

  if (is.null(evs_membership)) {
    result$EVS_Class <- "Original"
  } else {
    result <- result %>%
      dplyr::left_join(evs_membership, by = "feature_id")
    result$EVS_Class[is.na(result$EVS_Class)] <- "Remainder"
  }

  result$neglog10_Empirical_p <- safe_neglog10(result$Empirical_p)

  list(
    results = result,
    HCp = hc_p,
    Htau = h_tau,
    n_tested = nrow(result)
  )
}

# =============================================================================
# 9. RESULT SUMMARIES AND TABLE EXPORT
# =============================================================================

summarize_result <- function(result_obj, split_info = NULL) {
  df <- result_obj$results
  row <- data.frame(
    Comparison = unique(df$Comparison)[1],
    Analysis = unique(df$Analysis)[1],
    PAS_tested = nrow(df),
    H = sum(df$HBFSS_sig, na.rm = TRUE),
    Std = sum(df$Std, na.rm = TRUE),
    Str = sum(df$Strong, na.rm = TRUE),
    Wk = sum(df$Weak, na.rm = TRUE),
    Ovlp = sum(df$Overlap, na.rm = TRUE),
    Weak_CNH_pass = sum(df$Weak_CNH, na.rm = TRUE),
    H_Std = sum(df$HBFSS_sig & df$Std, na.rm = TRUE),
    H_Str = sum(df$HBFSS_sig & df$Strong, na.rm = TRUE),
    H_Wk = sum(df$HBFSS_sig & df$Weak_CNH, na.rm = TRUE),
    HCp = result_obj$HCp,
    Htau = result_obj$Htau,
    BH_FDR = DESEQ2_FDR,
    LFC_boundary = LFC_BOUNDARY,
    stringsAsFactors = FALSE
  )

  if (is.null(split_info)) {
    row$EVS_K <- NA_integer_
    row$Joint_n <- NA_integer_
    row$Disjoint_RT_n <- NA_integer_
    row$Disjoint_ZT_n <- NA_integer_
    row$Leading_n <- NA_integer_
    row$Remainder_n <- NA_integer_
  } else {
    row$EVS_K <- split_info$top_n
    row$Joint_n <- length(split_info$joint_ids)
    row$Disjoint_RT_n <- length(split_info$disjoint_treatment_ids)
    row$Disjoint_ZT_n <- length(split_info$disjoint_control_ids)
    row$Leading_n <- length(split_info$leading_ids)
    row$Remainder_n <- length(split_info$remainder_ids)
  }

  row
}

significant_table_for_comparison <- function(comparison_results) {
  out <- dplyr::bind_rows(lapply(comparison_results, function(x) x$results)) %>%
    dplyr::filter(Final_sig) %>%
    dplyr::mutate(
      Analysis = factor(Analysis, levels = analysis_view_levels)
    ) %>%
    dplyr::arrange(Analysis, dplyr::desc(HBFSS), DESeq2_padj, PAS_ID) %>%
    dplyr::select(
      Comparison,
      Analysis,
      PAS_ID,
      Gene,
      EVS_Class,
      Direction,
      Apeglm_LFC,
      Wald_stat,
      Empirical_p,
      HBFSS,
      HCp,
      Htau,
      DESeq2_padj,
      Strong_padj,
      Weak_padj,
      HBFSS_sig,
      Std,
      Strong,
      Weak,
      Overlap,
      Support,
      WTTS_PAS_n,
      APA_multi_PAS
    )

  out$Analysis <- as.character(out$Analysis)
  out
}

# =============================================================================
# 10. VOLCANO PLOTS
# =============================================================================

volcano_class <- function(df) {
  cls <- rep("Background", nrow(df))
  cls[df$Weak_CNH & !df$HBFSS_sig] <- "Weak CNH only"
  cls[df$Std & !df$HBFSS_sig] <- "Standard only"
  cls[df$Strong & !df$HBFSS_sig] <- "Strong only"
  cls[df$HBFSS_sig] <- "HBFSS only"
  cls[df$HBFSS_sig & df$Std] <- "HBFSS + Standard"
  cls[df$HBFSS_sig & df$Strong] <- "HBFSS + Strong"
  cls[df$HBFSS_sig & df$Weak_CNH] <- "HBFSS + Weak"
  cls
}

VOLCANO_CLASSES <- c(
  "Background",
  "Weak CNH only",
  "Standard only",
  "Strong only",
  "HBFSS only",
  "HBFSS + Standard",
  "HBFSS + Strong",
  "HBFSS + Weak"
)

VOLCANO_COLORS <- c(
  "Background" = "#BDBDBD",
  "Weak CNH only" = "#4EA3F1",
  "Standard only" = "#33A02C",
  "Strong only" = "#E31A1C",
  "HBFSS only" = "#6A3D9A",
  "HBFSS + Standard" = "#6A3D9A",
  "HBFSS + Strong" = "#6A3D9A",
  "HBFSS + Weak" = "#6A3D9A"
)

VOLCANO_SHAPES <- c(
  "Background" = 16,
  "Weak CNH only" = 17,
  "Standard only" = 18,
  "Strong only" = 15,
  "HBFSS only" = 8,
  "HBFSS + Standard" = 18,
  "HBFSS + Strong" = 15,
  "HBFSS + Weak" = 17
)

prepare_volcano_data <- function(comparison_results) {
  plot_df <- dplyr::bind_rows(lapply(comparison_results, function(x) x$results)) %>%
    dplyr::filter(
      is.finite(Apeglm_LFC), !is.na(Apeglm_LFC),
      is.finite(neglog10_Empirical_p), !is.na(neglog10_Empirical_p)
    ) %>%
    dplyr::mutate(Analysis = factor(Analysis, levels = analysis_view_levels))

  plot_df$PlotClass <- factor(volcano_class(plot_df), levels = VOLCANO_CLASSES)
  plot_df <- plot_df %>%
    dplyr::mutate(
      PlotOrder = dplyr::case_when(
        HBFSS_sig ~ 3L,
        Final_sig ~ 2L,
        Weak_CNH ~ 1L,
        TRUE ~ 0L
      )
    ) %>%
    dplyr::arrange(Analysis, PlotOrder)

  label_df <- plot_df %>%
    dplyr::filter(Final_sig) %>%
    dplyr::group_by(Analysis) %>%
    dplyr::arrange(
      dplyr::desc(neglog10_Empirical_p),
      dplyr::desc(abs(Apeglm_LFC)),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = VOLCANO_LABEL_N) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Label = ifelse(
        valid_gene_symbol(Gene),
        paste0(Gene, " [", PAS_ID, "]"),
        paste0("PAS ", PAS_ID)
      )
    )

  summary_df <- dplyr::bind_rows(lapply(comparison_results, summarize_result)) %>%
    dplyr::mutate(
      Analysis = factor(Analysis, levels = analysis_view_levels),
      Annotation = paste0(
        "H=", H,
        "  Std=", Std,
        "  Str=", Str,
        "  Wk=", Wk,
        "  Ovlp=", Ovlp,
        "\nHCp=", format_threshold(HCp),
        "  Htau=", format_threshold(Htau)
      )
    )

  boundary_rows <- lapply(seq_len(nrow(summary_df)), function(i) {
    a <- as.character(summary_df$Analysis[i])
    h_tau <- summary_df$Htau[i]
    panel_df <- plot_df[as.character(plot_df$Analysis) == a, , drop = FALSE]

    x_max <- safe_max(abs(panel_df$Apeglm_LFC))
    y_max <- safe_max(panel_df$neglog10_Empirical_p)
    if (!is.finite(x_max) || is.na(x_max) || x_max <= 0) x_max <- 2
    if (!is.finite(y_max) || is.na(y_max) || y_max <= 0) y_max <- 1
    x_max <- max(2, x_max)

    if (!is.finite(h_tau) || is.na(h_tau) || h_tau < 0) {
      return(data.frame())
    }

    x_min <- if (h_tau == 0) 0.02 else max(0.02, h_tau / (1.05 * y_max))
    x_min <- min(x_min, x_max * 0.98)
    pos <- seq(x_min, x_max, length.out = 400)

    dplyr::bind_rows(
      data.frame(
        Analysis = factor(a, levels = analysis_view_levels),
        Side = "negative",
        x = -rev(pos),
        y = h_tau / rev(pos),
        stringsAsFactors = FALSE
      ),
      data.frame(
        Analysis = factor(a, levels = analysis_view_levels),
        Side = "positive",
        x = pos,
        y = h_tau / pos,
        stringsAsFactors = FALSE
      )
    )
  })
  boundary_df <- dplyr::bind_rows(boundary_rows)

  hc_df <- summary_df %>%
    dplyr::filter(is.finite(HCp), !is.na(HCp)) %>%
    dplyr::mutate(HC_y = -log10(HCp))

  list(
    points = plot_df,
    labels = label_df,
    annotation = summary_df,
    boundary = boundary_df,
    hc = hc_df
  )
}

plot_volcano_panel <- function(comparison_results, comparison, treatment_label, control_label) {
  v <- prepare_volcano_data(comparison_results)

  p <- ggplot2::ggplot(
    v$points,
    ggplot2::aes(x = Apeglm_LFC, y = neglog10_Empirical_p)
  ) +
    ggplot2::geom_vline(
      xintercept = c(-LFC_BOUNDARY, LFC_BOUNDARY),
      linetype = "dotted",
      linewidth = 0.30,
      color = "grey45"
    ) +
    ggplot2::geom_hline(
      data = v$hc,
      ggplot2::aes(yintercept = HC_y),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.35,
      color = "#A65628"
    ) +
    ggplot2::geom_line(
      data = v$boundary,
      ggplot2::aes(x = x, y = y, group = interaction(Analysis, Side)),
      inherit.aes = FALSE,
      linewidth = 0.45,
      color = "#6A3D9A"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = PlotClass, shape = PlotClass),
      size = 1.45,
      alpha = 0.78,
      stroke = 0.30
    ) +
    ggrepel::geom_text_repel(
      data = v$labels,
      ggplot2::aes(label = Label),
      size = 2.15,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.15,
      seed = 42,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = v$annotation,
      ggplot2::aes(x = -Inf, y = Inf, label = Annotation),
      inherit.aes = FALSE,
      hjust = -0.03,
      vjust = 1.08,
      size = 2.45,
      lineheight = 0.95
    ) +
    ggplot2::facet_wrap(~ Analysis, ncol = 3, scales = "free") +
    ggplot2::scale_color_manual(values = VOLCANO_COLORS, drop = FALSE, name = NULL) +
    ggplot2::scale_shape_manual(values = VOLCANO_SHAPES, drop = FALSE, name = NULL) +
    ggplot2::labs(
      title = paste0(comparison, " differential PAS discovery"),
      subtitle = paste0(
        "Positive apeglm LFC = ", treatment_label, " higher than ", control_label,
        "; HBFSS uses empirical-null p-values"
      ),
      x = "apeglm-shrunken log2 fold change",
      y = expression(-log[10]("empirical p"))
    ) +
    manuscript_theme() +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 8),
      panel.spacing = grid::unit(0.75, "lines")
    )

  save_figure(
    p,
    paste0("Figure_", comparison, "_Volcanoes"),
    width = 16.5,
    height = 10.5
  )

  p
}

# =============================================================================
# 11. EVS SUPPORT FIGURE
# =============================================================================

plot_evs_rank_panel <- function(norm_evs, raw_evs, comparison) {
  df <- dplyr::bind_rows(norm_evs$rank_support, raw_evs$rank_support) %>%
    dplyr::mutate(
      Track = factor(Track, levels = c("NormEVS", "RawEVS")),
      PC1_abs_loading_plot = pmax(PC1_abs_loading, .Machine$double.xmin)
    )

  shade_df <- df %>% dplyr::distinct(Track, Condition)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = Rank, y = PC1_abs_loading_plot)
  ) +
    ggplot2::geom_rect(
      data = shade_df,
      ggplot2::aes(xmin = 1, xmax = EVS_TOP_N, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey92",
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.35, color = "#2C7FB8") +
    ggplot2::geom_vline(
      xintercept = EVS_TOP_N,
      linetype = "dashed",
      linewidth = 0.40,
      color = "#D95F02"
    ) +
    ggplot2::facet_grid(Track ~ Condition, scales = "free_y") +
    ggplot2::labs(
      title = paste0(comparison, " eigenvector-splitting rank profiles"),
      subtitle = "Absolute PC1 feature loadings ranked independently within RT and ZT; shaded region = top 5,000",
      x = "PC1 loading rank (highest to lowest)",
      y = "Absolute PC1 loading"
    ) +
    manuscript_theme()

  save_figure(
    p,
    paste0("Figure_", comparison, "_EVS_Loading_Ranks"),
    width = 12.5,
    height = 7.5
  )

  p
}

# =============================================================================
# 12. MAIN WTTS PIPELINE
# =============================================================================

all_results <- list()
all_summary_rows <- list()
comparison_outputs <- list()

for (i in seq_len(nrow(comparison_table))) {
  cmp_row <- comparison_table[i, , drop = FALSE]
  cmp <- prepare_comparison(cmp_row)

  message("Running ", cmp$comparison)

  norm_evs <- build_evs(
    cmp$counts,
    cmp$coldata,
    track = "NormEVS",
    treatment_label = cmp$treatment_label,
    control_label = cmp$control_label,
    top_n = EVS_TOP_N
  )

  raw_evs <- build_evs(
    cmp$counts,
    cmp$coldata,
    track = "RawEVS",
    treatment_label = cmp$treatment_label,
    control_label = cmp$control_label,
    top_n = EVS_TOP_N
  )

  original <- run_differential_analysis(
    cmp$counts,
    cmp$coldata,
    cmp$comparison,
    "Original",
    evs_membership = NULL
  )

  norm_lead <- run_differential_analysis(
    norm_evs$leading_counts,
    cmp$coldata,
    cmp$comparison,
    "NormEVS Lead",
    evs_membership = norm_evs$membership
  )

  norm_rem <- run_differential_analysis(
    norm_evs$remainder_counts,
    cmp$coldata,
    cmp$comparison,
    "NormEVS Rem",
    evs_membership = norm_evs$membership
  )

  raw_lead <- run_differential_analysis(
    raw_evs$leading_counts,
    cmp$coldata,
    cmp$comparison,
    "RawEVS Lead",
    evs_membership = raw_evs$membership
  )

  raw_rem <- run_differential_analysis(
    raw_evs$remainder_counts,
    cmp$coldata,
    cmp$comparison,
    "RawEVS Rem",
    evs_membership = raw_evs$membership
  )

  cmp_results <- list(
    Original = original,
    NormEVS_Lead = norm_lead,
    NormEVS_Rem = norm_rem,
    RawEVS_Lead = raw_lead,
    RawEVS_Rem = raw_rem
  )

  summary_cmp <- dplyr::bind_rows(
    summarize_result(original, NULL),
    summarize_result(norm_lead, norm_evs),
    summarize_result(norm_rem, norm_evs),
    summarize_result(raw_lead, raw_evs),
    summarize_result(raw_rem, raw_evs)
  ) %>%
    dplyr::mutate(Analysis = factor(Analysis, levels = analysis_view_levels)) %>%
    dplyr::arrange(Analysis)
  summary_cmp$Analysis <- as.character(summary_cmp$Analysis)

  sig_cmp <- significant_table_for_comparison(cmp_results)

  write_csv(
    sig_cmp,
    file.path(TABLE_DIR, paste0("Table_", cmp$comparison, "_Significant_PAS.csv"))
  )
  write_csv(
    summary_cmp,
    file.path(TABLE_DIR, paste0("Table_", cmp$comparison, "_Summary.csv"))
  )

  plot_volcano_panel(
    cmp_results,
    cmp$comparison,
    cmp$treatment_label,
    cmp$control_label
  )

  plot_evs_rank_panel(norm_evs, raw_evs, cmp$comparison)

  for (nm in names(cmp_results)) {
    key <- paste(cmp$comparison, nm, sep = "__")
    all_results[[key]] <- cmp_results[[nm]]$results
  }

  all_summary_rows[[cmp$comparison]] <- summary_cmp
  comparison_outputs[[cmp$comparison]] <- list(
    comparison = cmp,
    norm_evs = norm_evs,
    raw_evs = raw_evs,
    results = cmp_results
  )
}

all_summary <- dplyr::bind_rows(all_summary_rows)
write_csv(all_summary, file.path(TABLE_DIR, "Table_All_Comparisons_Summary.csv"))

# Overall discovery-count figure.
discovery_long <- all_summary %>%
  dplyr::select(Comparison, Analysis, H, Std, Str, Wk, Ovlp) %>%
  tidyr::pivot_longer(
    cols = c(H, Std, Str, Wk, Ovlp),
    names_to = "Method",
    values_to = "Significant_PAS"
  ) %>%
  dplyr::mutate(
    Analysis = factor(Analysis, levels = analysis_view_levels),
    Method = factor(Method, levels = c("H", "Std", "Str", "Wk", "Ovlp"))
  )

p_discovery <- ggplot2::ggplot(
  discovery_long,
  ggplot2::aes(x = Analysis, y = Significant_PAS, fill = Method)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78), width = 0.70) +
  ggplot2::facet_wrap(~ Comparison, ncol = 2, scales = "free_y") +
  ggplot2::scale_fill_manual(
    values = c(
      "H" = "#6A3D9A",
      "Std" = "#33A02C",
      "Str" = "#E31A1C",
      "Wk" = "#4EA3F1",
      "Ovlp" = "#54278F"
    ),
    labels = c(
      "H" = "HBFSS",
      "Std" = "Standard",
      "Str" = "Strong",
      "Wk" = "Weak (HBFSS-supported)",
      "Ovlp" = "HBFSS-DESeq2 overlap"
    ),
    name = NULL
  ) +
  ggplot2::labs(
    title = "Significant PAS counts by analysis view",
    x = NULL,
    y = "Significant PASs"
  ) +
  manuscript_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

save_figure(p_discovery, "Figure_All_Comparisons_Discovery_Counts", 14.5, 9.0)

# =============================================================================
# 13. OPTIONAL EVS K = 5,000 VS K = 4,077 SENSITIVITY
# =============================================================================

run_evs_cutoff_sensitivity <- function() {
  rows <- list()
  overlap_rows <- list()

  for (cmp_name in names(comparison_outputs)) {
    obj <- comparison_outputs[[cmp_name]]
    cmp <- obj$comparison

    evs_list <- list(
      K5000 = build_evs(
        cmp$counts, cmp$coldata, "NormEVS",
        cmp$treatment_label, cmp$control_label,
        top_n = EVS_TOP_N
      ),
      K4077 = build_evs(
        cmp$counts, cmp$coldata, "NormEVS",
        cmp$treatment_label, cmp$control_label,
        top_n = EMPIRICAL_EVS_TOP_N
      )
    )

    fit_list <- lapply(names(evs_list), function(k_name) {
      evs <- evs_list[[k_name]]
      fit <- run_differential_analysis(
        evs$leading_counts,
        cmp$coldata,
        cmp$comparison,
        paste0("Sensitivity ", k_name),
        evs_membership = evs$membership
      )
      list(evs = evs, fit = fit)
    })
    names(fit_list) <- names(evs_list)

    for (k_name in names(fit_list)) {
      fit <- fit_list[[k_name]]$fit
      evs <- fit_list[[k_name]]$evs
      df <- fit$results
      rows[[length(rows) + 1L]] <- data.frame(
        Comparison = cmp_name,
        K = evs$top_n,
        Leading_n = length(evs$leading_ids),
        Joint_n = length(evs$joint_ids),
        H = sum(df$HBFSS_sig, na.rm = TRUE),
        Std = sum(df$Std, na.rm = TRUE),
        Str = sum(df$Strong, na.rm = TRUE),
        Wk = sum(df$Weak, na.rm = TRUE),
        Ovlp = sum(df$Overlap, na.rm = TRUE),
        HCp = fit$HCp,
        Htau = fit$Htau,
        stringsAsFactors = FALSE
      )
    }

    a <- fit_list$K5000$fit$results
    b <- fit_list$K4077$fit$results

    for (method in c("HBFSS_sig", "Std", "Strong", "Weak")) {
      a_ids <- a$feature_id[!is.na(a[[method]]) & a[[method]]]
      b_ids <- b$feature_id[!is.na(b[[method]]) & b[[method]]]
      overlap_rows[[length(overlap_rows) + 1L]] <- data.frame(
        Comparison = cmp_name,
        Method = method,
        Shared = length(intersect(a_ids, b_ids)),
        K5000_only = length(setdiff(a_ids, b_ids)),
        K4077_only = length(setdiff(b_ids, a_ids)),
        stringsAsFactors = FALSE
      )
    }
  }

  summary_df <- dplyr::bind_rows(rows)
  overlap_df <- dplyr::bind_rows(overlap_rows)

  write_csv(summary_df, file.path(VALIDATION_DIR, "Table_EVS_K5000_vs_K4077_Summary.csv"))
  write_csv(overlap_df, file.path(VALIDATION_DIR, "Table_EVS_K5000_vs_K4077_Overlap.csv"))

  plot_df <- summary_df %>%
    dplyr::select(Comparison, K, H, Std, Str, Wk, Ovlp) %>%
    tidyr::pivot_longer(c(H, Std, Str, Wk, Ovlp), names_to = "Method", values_to = "Significant_PAS") %>%
    dplyr::mutate(K = factor(K, levels = c(EVS_TOP_N, EMPIRICAL_EVS_TOP_N)))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Comparison, y = Significant_PAS, fill = K)
  ) +
    ggplot2::geom_col(position = "dodge", width = 0.72) +
    ggplot2::facet_wrap(~ Method, scales = "free_y") +
    ggplot2::labs(
      title = "EVS cutoff sensitivity: K = 5,000 versus K = 4,077",
      x = NULL,
      y = "Significant PASs",
      fill = "EVS K"
    ) +
    manuscript_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

  save_figure(
    p,
    "Figure_EVS_K5000_vs_K4077",
    13.5,
    8.0,
    directory = VALIDATION_DIR
  )

  invisible(list(summary = summary_df, overlap = overlap_df))
}

if (isTRUE(RUN_EVS_CUTOFF_SENSITIVITY)) {
  message("Running EVS cutoff sensitivity: K=", EVS_TOP_N, " vs K=", EMPIRICAL_EVS_TOP_N)
  evs_cutoff_sensitivity_results <- run_evs_cutoff_sensitivity()
}

# =============================================================================
# 14. OPTIONAL FULLY ARTIFICIAL NEGATIVE-BINOMIAL VALIDATION
# =============================================================================

simulate_nb_dataset <- function(seed) {
  set.seed(seed)

  n_features <- SIMULATION_FEATURES
  n_group <- SIMULATION_SAMPLES_PER_GROUP
  n_total <- 2L * n_group

  base_mean <- exp(stats::rnorm(n_features, log(80), 1.0))
  dispersion <- pmax(stats::rgamma(n_features, shape = 2.5, rate = 18), 0.01)
  library_factor <- exp(stats::rnorm(n_total, 0, 0.20))

  n_de <- max(2L, round(n_features * SIMULATION_DE_FRACTION))
  de_ids <- sample(seq_len(n_features), n_de, replace = FALSE)
  n_strong <- max(1L, round(n_de * SIMULATION_STRONG_FRACTION_WITHIN_DE))
  strong_ids <- sample(de_ids, n_strong, replace = FALSE)
  weak_ids <- setdiff(de_ids, strong_ids)

  true_lfc <- rep(0, n_features)
  sign_vec <- sample(c(-1, 1), n_de, replace = TRUE)
  names(sign_vec) <- de_ids

  if (length(strong_ids) > 0L) {
    true_lfc[strong_ids] <- sign_vec[as.character(strong_ids)] *
      stats::runif(length(strong_ids), SIMULATION_STRONG_LFC_RANGE[1], SIMULATION_STRONG_LFC_RANGE[2])
  }
  if (length(weak_ids) > 0L) {
    true_lfc[weak_ids] <- sign_vec[as.character(weak_ids)] *
      stats::runif(length(weak_ids), SIMULATION_WEAK_LFC_RANGE[1], SIMULATION_WEAK_LFC_RANGE[2])
  }

  counts <- matrix(0L, nrow = n_features, ncol = n_total)
  for (j in seq_len(n_total)) {
    is_treatment <- j > n_group
    mu <- base_mean * library_factor[j]
    if (is_treatment) mu <- mu * (2 ^ true_lfc)
    size <- 1 / dispersion
    counts[, j] <- stats::rnbinom(n_features, mu = mu, size = size)
  }

  rownames(counts) <- paste0("SIM_", seq_len(n_features))
  colnames(counts) <- c(paste0("C_", seq_len(n_group)), paste0("T_", seq_len(n_group)))

  coldata <- data.frame(
    condition = factor(
      c(rep("control", n_group), rep("treatment", n_group)),
      levels = c("control", "treatment")
    ),
    row.names = colnames(counts)
  )

  annotation <- data.frame(
    feature_id = rownames(counts),
    PAS_ID = rownames(counts),
    Gene = rownames(counts),
    Chromosome = NA_character_,
    Strand = NA_character_,
    Peak = NA_character_,
    gene_key = rownames(counts),
    stringsAsFactors = FALSE
  )

  truth <- data.frame(
    feature_id = rownames(counts),
    True_LFC = true_lfc,
    True_DE = seq_len(n_features) %in% de_ids,
    True_Strong = seq_len(n_features) %in% strong_ids,
    True_Weak = seq_len(n_features) %in% weak_ids,
    stringsAsFactors = FALSE
  )

  list(counts = counts, coldata = coldata, annotation = annotation, truth = truth)
}

binary_metrics <- function(pred, truth) {
  pred <- as.logical(pred)
  truth <- as.logical(truth)
  tp <- sum(pred & truth, na.rm = TRUE)
  fp <- sum(pred & !truth, na.rm = TRUE)
  fn <- sum(!pred & truth, na.rm = TRUE)
  precision <- if ((tp + fp) == 0L) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0L) NA_real_ else tp / (tp + fn)
  f1 <- if (!is.finite(precision) || !is.finite(recall) || (precision + recall) == 0) {
    NA_real_
  } else {
    2 * precision * recall / (precision + recall)
  }
  fdr <- if ((tp + fp) == 0L) 0 else fp / (tp + fp)
  data.frame(TP = tp, FP = fp, FN = fn, Precision = precision, Recall = recall, F1 = f1, FDR = fdr)
}

run_simulation_once <- function(rep_id) {
  sim <- simulate_nb_dataset(SIMULATION_SEED + rep_id)

  sim_gene_summary <- data.frame(
    gene_key = sim$annotation$gene_key,
    WTTS_Gene = sim$annotation$Gene,
    WTTS_PAS_n = 1L,
    APA_multi_PAS = FALSE,
    stringsAsFactors = FALSE
  )

  fit <- run_differential_analysis(
    sim$counts,
    sim$coldata,
    comparison = paste0("Simulation_", rep_id),
    analysis_view = "Original",
    evs_membership = NULL,
    annotation_df = sim$annotation,
    gene_pas_summary = sim_gene_summary
  )

  df <- fit$results %>% dplyr::left_join(sim$truth, by = "feature_id")

  method_rows <- list(
    cbind(data.frame(Method = "Standard", Truth = "All_DE"), binary_metrics(df$Std, df$True_DE)),
    cbind(data.frame(Method = "HBFSS", Truth = "All_DE"), binary_metrics(df$HBFSS_sig, df$True_DE)),
    cbind(data.frame(Method = "Strong", Truth = "Strong_DE"), binary_metrics(df$Strong, df$True_Strong)),
    cbind(data.frame(Method = "Weak_HBFSS", Truth = "Weak_DE"), binary_metrics(df$Weak, df$True_Weak))
  )

  dplyr::bind_rows(method_rows) %>%
    dplyr::mutate(Replicate = rep_id, HCp = fit$HCp, Htau = fit$Htau)
}

run_synthetic_validation <- function() {
  out_dir <- file.path(VALIDATION_DIR, "synthetic")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- lapply(seq_len(SIMULATION_REPLICATES), function(r) {
    message("Synthetic validation replicate ", r, "/", SIMULATION_REPLICATES)
    tryCatch(
      run_simulation_once(r),
      error = function(e) {
        data.frame(
          Method = NA_character_, Truth = NA_character_, TP = NA, FP = NA, FN = NA,
          Precision = NA_real_, Recall = NA_real_, F1 = NA_real_, FDR = NA_real_,
          Replicate = r, HCp = NA_real_, Htau = NA_real_, Error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  })

  run_df <- dplyr::bind_rows(rows)
  write_csv(run_df, file.path(out_dir, "Table_Synthetic_Validation_Replicates.csv"))

  ok <- run_df %>% dplyr::filter(!is.na(Method))
  summary_df <- ok %>%
    dplyr::group_by(Method, Truth) %>%
    dplyr::summarise(
      Replicates = dplyr::n(),
      F1_mean = mean(F1, na.rm = TRUE),
      F1_sd = stats::sd(F1, na.rm = TRUE),
      FDR_mean = mean(FDR, na.rm = TRUE),
      Recall_mean = mean(Recall, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv(summary_df, file.path(out_dir, "Table_Synthetic_Validation_Summary.csv"))

  if (nrow(ok) > 0L) {
    metric_long <- ok %>%
      dplyr::select(Method, Replicate, F1, FDR, Recall) %>%
      tidyr::pivot_longer(c(F1, FDR, Recall), names_to = "Metric", values_to = "Value")

    p <- ggplot2::ggplot(metric_long, ggplot2::aes(x = Method, y = Value, fill = Method)) +
      ggplot2::geom_boxplot(width = 0.65, outlier.size = 0.6) +
      ggplot2::facet_wrap(~ Metric, scales = "free_y") +
      ggplot2::labs(
        title = "Fully artificial negative-binomial validation",
        x = NULL,
        y = NULL,
        fill = NULL
      ) +
      manuscript_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

    ggplot2::ggsave(
      filename = file.path(out_dir, "Figure_Synthetic_Validation.png"),
      plot = p,
      width = 11.5,
      height = 7.5,
      units = "in",
      dpi = FIGURE_DPI,
      bg = "white"
    )
    if (isTRUE(EXPORT_PDF)) {
      ggplot2::ggsave(
        filename = file.path(out_dir, "Figure_Synthetic_Validation.pdf"),
        plot = p,
        width = 11.5,
        height = 7.5,
        units = "in",
        bg = "white"
      )
    }
  }

  invisible(list(replicates = run_df, summary = summary_df))
}

if (isTRUE(RUN_SYNTHETIC_VALIDATION)) {
  synthetic_validation_results <- run_synthetic_validation()
}

# =============================================================================
# 15. FINAL ANALYTICAL STAGE: 3'aTWAS ORTHOLOG OVERLAP
# =============================================================================

read_twas_study <- function(path) {
  # The supplied CSV has one descriptive title row followed by the true header.
  twas <- utils::read.csv(
    path,
    skip = 1,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  assert_columns(
    twas,
    c("Disease", "PANEL", "Transcript ID", "Gene symbol", "3'aTWAS.Z", "3'aTWAS.P"),
    "3'aTWAS file"
  )

  twas <- twas %>%
    dplyr::transmute(
      Disease = trimws(as.character(Disease)),
      Panel = trimws(as.character(PANEL)),
      TWAS_Transcript = trimws(as.character(`Transcript ID`)),
      Human_TWAS_Gene = trimws(as.character(`Gene symbol`)),
      TWAS_Z = suppressWarnings(as.numeric(`3'aTWAS.Z`)),
      TWAS_P = suppressWarnings(as.numeric(`3'aTWAS.P`)),
      COLOC_PP4 = if ("COLOC.PP4" %in% names(twas)) suppressWarnings(as.numeric(COLOC.PP4)) else NA_real_
    ) %>%
    dplyr::filter(!is.na(Human_TWAS_Gene), nzchar(Human_TWAS_Gene))

  twas
}

build_twas_ortholog_map <- function(twas) {
  human_genes <- sort(unique(twas$Human_TWAS_Gene))

  orth <- babelgene::orthologs(
    genes = human_genes,
    species = TWAS_TARGET_SPECIES,
    human = TRUE,
    min_support = TWAS_ORTHOLOG_MIN_SUPPORT,
    top = TRUE
  )

  orth <- as.data.frame(orth, stringsAsFactors = FALSE)
  assert_columns(orth, c("human_symbol", "symbol", "support_n"), "babelgene ortholog output")

  orth_best <- data.frame(
    Human_TWAS_Gene = as.character(orth$human_symbol),
    Rat_Ortholog = as.character(orth$symbol),
    Ortholog_support_n = as.integer(orth$support_n),
    Ortholog_support = if ("support" %in% names(orth)) as.character(orth$support) else NA_character_,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(Human_TWAS_Gene, dplyr::desc(Ortholog_support_n), Rat_Ortholog) %>%
    dplyr::distinct(Human_TWAS_Gene, .keep_all = TRUE)

  mapping <- data.frame(Human_TWAS_Gene = human_genes, stringsAsFactors = FALSE) %>%
    dplyr::left_join(orth_best, by = "Human_TWAS_Gene") %>%
    dplyr::mutate(gene_key = gene_key(Rat_Ortholog))

  mapping
}

build_twas_rat_metadata <- function(twas, mapping) {
  mapped_rows <- twas %>%
    dplyr::left_join(mapping, by = "Human_TWAS_Gene") %>%
    dplyr::filter(!is.na(gene_key))

  mapped_rows %>%
    dplyr::group_by(gene_key) %>%
    dplyr::summarise(
      Rat_Ortholog = dplyr::first(Rat_Ortholog[!is.na(Rat_Ortholog)]),
      Human_TWAS_Genes = collapse_unique(Human_TWAS_Gene),
      TWAS_Diseases = collapse_unique(Disease),
      TWAS_Panels_n = dplyr::n_distinct(Panel),
      TWAS_Records_n = dplyr::n(),
      TWAS_min_P = safe_min(TWAS_P),
      TWAS_max_abs_Z = safe_max(abs(TWAS_Z)),
      TWAS_max_COLOC_PP4 = safe_max(COLOC_PP4),
      Ortholog_support_n = safe_max(Ortholog_support_n),
      .groups = "drop"
    ) %>%
    dplyr::left_join(wtts_gene_pas_summary, by = "gene_key")
}

build_twas_overlap_gene_table <- function(all_results, twas_rat_meta) {
  all_df <- dplyr::bind_rows(all_results) %>%
    dplyr::filter(!is.na(gene_key), Std | HBFSS_sig) %>%
    dplyr::inner_join(twas_rat_meta, by = "gene_key")

  if (nrow(all_df) == 0L) {
    return(data.frame())
  }

  all_df %>%
    dplyr::group_by(
      Comparison,
      Analysis,
      gene_key,
      Rat_Ortholog,
      Human_TWAS_Genes,
      TWAS_Diseases,
      TWAS_Panels_n,
      TWAS_Records_n,
      TWAS_min_P,
      TWAS_max_abs_Z,
      TWAS_max_COLOC_PP4,
      Ortholog_support_n,
      WTTS_PAS_n,
      APA_multi_PAS
    ) %>%
    dplyr::summarise(
      Standard = any(Std, na.rm = TRUE),
      HBFSS = any(HBFSS_sig, na.rm = TRUE),
      Standard_PAS_n = dplyr::n_distinct(PAS_ID[Std]),
      HBFSS_PAS_n = dplyr::n_distinct(PAS_ID[HBFSS_sig]),
      Shared_PAS_n = dplyr::n_distinct(PAS_ID[Std & HBFSS_sig]),
      Union_significant_PAS_n = dplyr::n_distinct(PAS_ID[Std | HBFSS_sig]),
      Significant_PAS_IDs = collapse_unique(PAS_ID[Std | HBFSS_sig]),
      Best_abs_Apeglm_LFC = safe_max(abs(Apeglm_LFC[Std | HBFSS_sig])),
      Min_DESeq2_padj = safe_min(DESeq2_padj[Std]),
      Max_HBFSS = safe_max(HBFSS[HBFSS_sig]),
      Min_HBFSS_empirical_p = safe_min(Empirical_p[HBFSS_sig]),
      HCp = dplyr::first(HCp),
      Htau = dplyr::first(Htau),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Both_methods = Standard & HBFSS,
      Method_support = dplyr::case_when(
        Both_methods ~ "Both",
        Standard ~ "Standard only",
        HBFSS ~ "HBFSS only",
        TRUE ~ NA_character_
      ),
      DE_multi_PAS = Union_significant_PAS_n >= 2L,
      Analysis = factor(Analysis, levels = analysis_view_levels)
    ) %>%
    dplyr::arrange(Comparison, Analysis, Rat_Ortholog) %>%
    dplyr::mutate(Analysis = as.character(Analysis))
}

build_twas_overlap_summary <- function(gene_table) {
  if (nrow(gene_table) == 0L) return(data.frame())

  gene_table %>%
    dplyr::group_by(Comparison, Analysis) %>%
    dplyr::summarise(
      TWAS_Standard_genes = sum(Standard, na.rm = TRUE),
      TWAS_HBFSS_genes = sum(HBFSS, na.rm = TRUE),
      TWAS_Both_genes = sum(Both_methods, na.rm = TRUE),
      TWAS_Standard_only_genes = sum(Standard & !HBFSS, na.rm = TRUE),
      TWAS_HBFSS_only_genes = sum(HBFSS & !Standard, na.rm = TRUE),
      TWAS_Union_genes = dplyr::n(),
      TWAS_Standard_APA_genes = sum(Standard & APA_multi_PAS, na.rm = TRUE),
      TWAS_HBFSS_APA_genes = sum(HBFSS & APA_multi_PAS, na.rm = TRUE),
      TWAS_DE_multi_PAS_genes = sum(DE_multi_PAS, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Analysis = factor(Analysis, levels = analysis_view_levels)) %>%
    dplyr::arrange(Comparison, Analysis) %>%
    dplyr::mutate(Analysis = as.character(Analysis))
}

plot_twas_overlap_counts <- function(summary_df) {
  if (nrow(summary_df) == 0L) return(NULL)

  plot_df <- summary_df %>%
    dplyr::select(
      Comparison, Analysis,
      TWAS_Standard_only_genes, TWAS_HBFSS_only_genes, TWAS_Both_genes
    ) %>%
    tidyr::pivot_longer(
      c(TWAS_Standard_only_genes, TWAS_HBFSS_only_genes, TWAS_Both_genes),
      names_to = "Category",
      values_to = "TWAS_genes"
    ) %>%
    dplyr::mutate(
      Analysis = factor(Analysis, levels = analysis_view_levels),
      Category = factor(
        Category,
        levels = c(
          "TWAS_Standard_only_genes",
          "TWAS_HBFSS_only_genes",
          "TWAS_Both_genes"
        ),
        labels = c("Standard only", "HBFSS only", "Both methods")
      )
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Analysis, y = TWAS_genes, fill = Category)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.76), width = 0.68) +
    ggplot2::geom_text(
      ggplot2::aes(label = TWAS_genes),
      position = ggplot2::position_dodge(width = 0.76),
      vjust = -0.25,
      size = 2.7
    ) +
    ggplot2::facet_wrap(~ Comparison, ncol = 2, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = c(
        "Standard only" = "#33A02C",
        "HBFSS only" = "#6A3D9A",
        "Both methods" = "#1F78B4"
      ),
      name = NULL
    ) +
    ggplot2::labs(
      title = "3'aTWAS ortholog overlap: Standard DESeq2 versus HBFSS",
      x = NULL,
      y = "Unique TWAS ortholog genes"
    ) +
    manuscript_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

  save_figure(p, "Figure_TWAS_Overlap_Method_Comparison", 14.0, 8.8)
  p
}

plot_twas_gene_matrix <- function(gene_table, comparison) {
  df <- gene_table %>% dplyr::filter(Comparison == comparison)
  if (nrow(df) == 0L) return(NULL)

  long <- df %>%
    dplyr::select(
      Comparison, Analysis, Rat_Ortholog, Human_TWAS_Genes,
      APA_multi_PAS, Standard, HBFSS, Standard_PAS_n, HBFSS_PAS_n
    ) %>%
    tidyr::pivot_longer(
      cols = c(Standard, HBFSS),
      names_to = "Method",
      values_to = "Identified"
    ) %>%
    dplyr::filter(Identified) %>%
    dplyr::mutate(
      Significant_PAS_n = ifelse(Method == "Standard", Standard_PAS_n, HBFSS_PAS_n),
      Method = factor(Method, levels = c("Standard", "HBFSS"), labels = c("Standard DESeq2", "HBFSS")),
      Analysis = factor(Analysis, levels = analysis_view_levels),
      Gene_label = paste0(
        Rat_Ortholog,
        " [", Human_TWAS_Genes, "]",
        ifelse(APA_multi_PAS, " *", "")
      )
    )

  gene_order <- long %>%
    dplyr::group_by(Gene_label) %>%
    dplyr::summarise(
      Views = dplyr::n_distinct(Analysis),
      PAS = sum(Significant_PAS_n, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Views, PAS, Gene_label) %>%
    dplyr::pull(Gene_label)

  long$Gene_label <- factor(long$Gene_label, levels = gene_order)

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = Analysis, y = Gene_label, color = Method, shape = Method, size = Significant_PAS_n)
  ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.42),
      alpha = 0.90,
      stroke = 0.45
    ) +
    ggplot2::scale_color_manual(
      values = c("Standard DESeq2" = "#33A02C", "HBFSS" = "#6A3D9A"),
      name = NULL
    ) +
    ggplot2::scale_shape_manual(
      values = c("Standard DESeq2" = 18, "HBFSS" = 8),
      name = NULL
    ) +
    ggplot2::scale_size_continuous(range = c(2.0, 5.5), name = "Significant PASs") +
    ggplot2::labs(
      title = paste0(comparison, " 3'aTWAS ortholog genes"),
      subtitle = "* WTTS gene is represented by >=2 PAS features",
      x = NULL,
      y = "Rat ortholog [human TWAS gene]"
    ) +
    manuscript_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 7.2)
    )

  n_genes <- dplyr::n_distinct(long$Gene_label)
  height <- min(24, max(6.5, 3.0 + 0.22 * n_genes))
  save_figure(p, paste0("Figure_", comparison, "_TWAS_Ortholog_Genes"), 13.5, height)
  p
}

twas_file <- resolve_existing_file(TWAS_FILE_CANDIDATES, "3'aTWAS file")
message("3'aTWAS file: ", twas_file)

twas <- read_twas_study(twas_file)
twas_mapping <- build_twas_ortholog_map(twas)
twas_rat_meta <- build_twas_rat_metadata(twas, twas_mapping)
twas_gene_overlap <- build_twas_overlap_gene_table(all_results, twas_rat_meta)
twas_overlap_summary <- build_twas_overlap_summary(twas_gene_overlap)

write_csv(twas_gene_overlap, file.path(TABLE_DIR, "Table_TWAS_Overlap_Genes.csv"))
write_csv(twas_overlap_summary, file.path(TABLE_DIR, "Table_TWAS_Overlap_Summary.csv"))

plot_twas_overlap_counts(twas_overlap_summary)
for (cmp_name in comparison_table$comparison) {
  plot_twas_gene_matrix(twas_gene_overlap, cmp_name)
}

# =============================================================================
# 16. REPRODUCIBILITY RECORD
# =============================================================================

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "SessionInfo.txt"))

cat("\n============================================================\n")
cat("SEQUENCE manuscript analysis complete\n")
cat("Output directory: ", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Main comparisons: ", paste(comparison_table$comparison, collapse = ", "), "\n", sep = "")
cat("Analysis views: ", paste(analysis_view_levels, collapse = ", "), "\n", sep = "")
cat("TWAS overlap stage: complete\n")
cat("============================================================\n")
