#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(width = 140)

# =============================================================================
# HBFSS COUNT-DATA SIMULATION ARM
# Standalone manuscript simulation for DESeq2, empirical-null p-values,
# higher-criticism thresholding, and HBFSS classification.
# =============================================================================

read_bool <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

read_int <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || !is.finite(value) || value < minimum) return(as.integer(default))
  as.integer(value)
}

read_num <- function(name, default, minimum = -Inf, maximum = Inf) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || !is.finite(value) || value < minimum || value > maximum) return(as.numeric(default))
  as.numeric(value)
}

# =============================================================================
# USER SETTINGS
# =============================================================================

OUTPUT_FOLDER_NAME <- "manuscript_final_clean"
SIMULATION_FOLDER_NAME <- "simulation_validation"

RESET_OUTPUT <- read_bool("SEQUENCE_RESET_SIM_OUTPUT", TRUE)
WRITE_PDF <- read_bool("SEQUENCE_WRITE_PDF", TRUE)
GIT_PUSH_AFTER_SUCCESS <- read_bool("SEQUENCE_GIT_PUSH", TRUE)

SEED <- read_int("SEQUENCE_SIM_SEED", 42L)
N_FEATURES <- read_int("SEQUENCE_SIM_FEATURES", 5000L)
N_SAMPLES_PER_GROUP <- read_int("SEQUENCE_SIM_SAMPLES_PER_GROUP", 6L)
N_REPS <- read_int("SEQUENCE_SIM_REPS", 12L)

BASE_MEAN <- read_num("SEQUENCE_SIM_BASE_MEAN", 200)
DISP_NULL <- read_num("SEQUENCE_SIM_DISP_NULL", 0.10)
DISP_DE <- read_num("SEQUENCE_SIM_DISP_DE", 0.15)

ALPHA_STANDARD <- 0.10
ALPHA_STRONG <- 0.10
ALPHA_WEAK <- 0.10
LFC_BOUNDARY <- 1.00

PROB_FLOOR <- .Machine$double.xmin
PLOT_FLOOR <- 1e-16
HC_INVALID_BELOW <- 0.001
HC_INVALID_AT_OR_ABOVE <- 0.95

DE_FRACTIONS <- c(0.05, 0.10)
FIXED_LFC_MAGNITUDES <- c(0.50, 1.00, 2.00)
WEAK_LFC_MIN <- 0.20
WEAK_LFC_MAX <- 0.80
NULL_INFLATION_FRACTIONS <- c(0.00, 0.10)

FIGURE_DPI <- read_int("SEQUENCE_SIM_DPI", 300L)
BASE_THEME_SIZE <- 9

MAX_FAILURE_FRACTION <- read_num("SEQUENCE_SIM_MAX_FAILURE_FRACTION", 0.05, minimum = 0, maximum = 1)

GIT_REMOTE <- Sys.getenv("SEQUENCE_GIT_REMOTE", unset = "origin")
GIT_BRANCH <- Sys.getenv("SEQUENCE_GIT_BRANCH", unset = NA_character_)
GIT_COMMIT_MESSAGE <- Sys.getenv(
  "SEQUENCE_GIT_COMMIT_MESSAGE",
  unset = paste0("Update HBFSS simulation outputs ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)

# =============================================================================
# PACKAGES
# =============================================================================

REQUIRED_PACKAGES <- c("DESeq2", "apeglm", "fdrtool", "ggplot2", "dplyr", "tidyr", "grid")

missing_packages <- REQUIRED_PACKAGES[
  !vapply(REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop("Install missing package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(fdrtool)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

# =============================================================================
# PATHS
# =============================================================================

current_script <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0L) return(NA_character_)
  normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = FALSE)
}

find_repo_root <- function() {
  starts <- unique(c(dirname(current_script()), getwd()))
  starts <- starts[!is.na(starts) & dir.exists(starts)]

  for (start in starts) {
    here <- normalizePath(start, winslash = "/", mustWork = TRUE)
    repeat {
      if (dir.exists(file.path(here, ".git"))) return(here)
      parent <- dirname(here)
      if (identical(parent, here)) break
      here <- parent
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

REPO_ROOT <- find_repo_root()
OUTPUT_DIR <- file.path(REPO_ROOT, "exports", OUTPUT_FOLDER_NAME)
SIM_DIR <- file.path(OUTPUT_DIR, SIMULATION_FOLDER_NAME)
FIG_DIR <- file.path(SIM_DIR, "figures")
PNG_DIR <- file.path(FIG_DIR, "png")
PDF_DIR <- file.path(FIG_DIR, "pdf")
LOG_FILE <- file.path(SIM_DIR, "simulation_run.log")
STATUS_FILE <- file.path(SIM_DIR, "simulation_status.txt")

if (isTRUE(RESET_OUTPUT) && dir.exists(SIM_DIR)) {
  unlink(SIM_DIR, recursive = TRUE, force = TRUE)
}

dir.create(SIM_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PNG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PDF_DIR, recursive = TRUE, showWarnings = FALSE)

repo_relative <- function(path) {
  root <- normalizePath(REPO_ROOT, winslash = "/", mustWork = TRUE)
  full <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  if (startsWith(full, prefix)) return(substr(full, nchar(prefix) + 1L, nchar(full)))
  full
}

log_msg <- function(...) {
  text <- paste0(...)
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", text)
  message(line)
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

write_status <- function(status, detail = "") {
  lines <- c(
    paste0("status=", status),
    paste0("time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("detail=", detail)
  )
  writeLines(lines, STATUS_FILE)
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE)
  invisible(path)
}

# =============================================================================
# NUMERIC HELPERS
# =============================================================================

clip_probability <- function(x, floor_value = PROB_FLOOR) {
  y <- suppressWarnings(as.numeric(x))
  y[!is.finite(y)] <- NA_real_
  ok <- !is.na(y)
  y[ok] <- pmin(pmax(y[ok], floor_value), 1 - 1e-12)
  y
}

neglog10_safe <- function(p, floor_value = PLOT_FLOOR) {
  -log10(clip_probability(p, floor_value = floor_value))
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

min_finite <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y <- y[is.finite(y) & !is.na(y)]
  if (length(y) == 0L) return(NA_real_)
  min(y)
}

max_finite <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y <- y[is.finite(y) & !is.na(y)]
  if (length(y) == 0L) return(NA_real_)
  max(y)
}

# =============================================================================
# SIMULATION DATA GENERATION
# =============================================================================

make_simulation_grid <- function() {
  fixed_grid <- expand.grid(
    de_fraction = DE_FRACTIONS,
    lfc_magnitude = FIXED_LFC_MAGNITUDES,
    null_inflation = NULL_INFLATION_FRACTIONS,
    profile = "fixed",
    stringsAsFactors = FALSE
  )

  weak_grid <- expand.grid(
    de_fraction = DE_FRACTIONS,
    lfc_magnitude = NA_real_,
    null_inflation = NULL_INFLATION_FRACTIONS,
    profile = "weak_mixture",
    stringsAsFactors = FALSE
  )

  dplyr::bind_rows(fixed_grid, weak_grid)
}

lfc_label_from_profile <- function(profile, lfc_magnitude) {
  if (identical(as.character(profile), "weak_mixture")) {
    return(paste0("Weak mix\n|LFC| ", WEAK_LFC_MIN, "-", WEAK_LFC_MAX))
  }
  paste0("Fixed\n|LFC| ", lfc_magnitude)
}

simulate_counts <- function(n_features, n_samples_per_group, de_fraction, lfc_magnitude, profile) {
  n_de <- max(2L, min(n_features - 2L, round(n_features * de_fraction)))
  feature_id <- paste0("sim_", seq_len(n_features))

  de_index <- sample(seq_len(n_features), n_de, replace = FALSE)

  base_mu <- stats::rgamma(n_features, shape = 2.5, scale = BASE_MEAN / 2.5)
  base_mu <- pmax(base_mu, 2)

  dispersion <- rep(DISP_NULL, n_features)
  dispersion[de_index] <- DISP_DE
  dispersion <- dispersion * exp(stats::rnorm(n_features, mean = 0, sd = 0.25))
  dispersion <- pmin(pmax(dispersion, 0.01), 1.50)

  direction <- sample(c(-1, 1), n_de, replace = TRUE)

  abs_lfc <- if (identical(profile, "weak_mixture")) {
    stats::runif(n_de, min = WEAK_LFC_MIN, max = WEAK_LFC_MAX)
  } else {
    rep(lfc_magnitude, n_de)
  }

  true_lfc <- rep(0, n_features)
  true_lfc[de_index] <- direction * abs_lfc

  control_mu <- base_mu
  treatment_mu <- base_mu * 2^true_lfc

  control_lib <- exp(stats::rnorm(n_samples_per_group, mean = 0, sd = 0.12))
  treatment_lib <- exp(stats::rnorm(n_samples_per_group, mean = 0, sd = 0.12))
  control_lib <- control_lib / exp(mean(log(control_lib)))
  treatment_lib <- treatment_lib / exp(mean(log(treatment_lib)))

  draw_group <- function(mu, lib_factor, disp) {
    out <- matrix(0L, nrow = length(mu), ncol = length(lib_factor))
    for (j in seq_along(lib_factor)) {
      out[, j] <- stats::rnbinom(
        n = length(mu),
        mu = pmax(mu * lib_factor[j], 1e-3),
        size = 1 / disp
      )
    }
    out
  }

  counts <- cbind(
    draw_group(control_mu, control_lib, dispersion),
    draw_group(treatment_mu, treatment_lib, dispersion)
  )

  storage.mode(counts) <- "integer"
  rownames(counts) <- feature_id
  colnames(counts) <- c(
    paste0("ctrl_", seq_len(n_samples_per_group)),
    paste0("trt_", seq_len(n_samples_per_group))
  )

  truth <- data.frame(
    feature_id = feature_id,
    is_de = seq_len(n_features) %in% de_index,
    true_lfc = true_lfc,
    true_abs_lfc = abs(true_lfc),
    true_weak = seq_len(n_features) %in% de_index & abs(true_lfc) < LFC_BOUNDARY,
    true_strong = seq_len(n_features) %in% de_index & abs(true_lfc) >= LFC_BOUNDARY,
    base_mean = base_mu,
    dispersion = dispersion,
    stringsAsFactors = FALSE
  )

  list(counts = counts, truth = truth)
}

inflate_null_wald <- function(wald, is_de, fraction, inflation_sd = 1.75) {
  if (!is.finite(fraction) || fraction <= 0) return(wald)

  out <- wald
  null_index <- which(!is_de & is.finite(out) & !is.na(out))
  n_target <- round(length(null_index) * fraction)

  if (n_target <= 0L) return(out)

  target <- sample(null_index, n_target, replace = FALSE)
  out[target] <- stats::rnorm(n_target, mean = 0, sd = inflation_sd)
  out
}

# =============================================================================
# DESEQ2, EMPIRICAL NULL, HC, HBFSS
# =============================================================================

condition_coef_name <- function(dds) {
  names_available <- resultsNames(dds)

  if ("condition_trt_vs_untrt" %in% names_available) {
    return("condition_trt_vs_untrt")
  }

  hit <- grep("condition.*trt.*vs.*untrt", names_available, value = TRUE)
  if (length(hit) > 0L) return(hit[1])

  stop("Could not identify treatment coefficient. Available names: ", paste(names_available, collapse = ", "), call. = FALSE)
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

fit_empirical_null <- function(wald) {
  valid <- is.finite(wald) & !is.na(wald)

  if (sum(valid) < 5L) {
    p <- rep(NA_real_, length(wald))
    return(list(empirical_p = p, empirical_bh = p, method = "failed_too_few_statistics"))
  }

  z <- wald[valid]

  run_fdrtool <- function(method_name) {
    fit <- tryCatch(
      fdrtool::fdrtool(
        z,
        statistic = "normal",
        plot = FALSE,
        verbose = FALSE,
        cutoff.method = method_name,
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
  fit <- run_fdrtool("fndr")

  if (is.null(fit)) {
    method_used <- "pct0"
    fit <- run_fdrtool("pct0")
  }

  if (is.null(fit)) {
    method_used <- "theoretical_normal"
    p_z <- 2 * stats::pnorm(-abs(z))
    fit <- list(pval = p_z)
  }

  empirical_p <- rep(NA_real_, length(wald))
  empirical_bh <- rep(NA_real_, length(wald))

  empirical_p[valid] <- clip_probability(fit$pval)

  ok <- !is.na(empirical_p) & is.finite(empirical_p)
  empirical_bh[ok] <- p.adjust(empirical_p[ok], method = "BH")

  list(empirical_p = empirical_p, empirical_bh = empirical_bh, method = method_used)
}

hc_threshold <- function(empirical_p) {
  p <- clip_probability(empirical_p)
  p <- sort(p[is.finite(p) & !is.na(p)], decreasing = FALSE)

  if (length(p) < 5L) return(NA_real_)

  threshold <- tryCatch(
    suppressWarnings(as.numeric(fdrtool::hc.thresh(p)[1])),
    error = function(e) NA_real_
  )

  if (!is.finite(threshold) || is.na(threshold)) return(NA_real_)
  if (threshold <= 0 || threshold >= 1) return(NA_real_)
  if (threshold < HC_INVALID_BELOW) return(NA_real_)
  if (threshold >= HC_INVALID_AT_OR_ABOVE) return(NA_real_)

  threshold
}

run_one_analysis <- function(sim, null_inflation) {
  counts_mat <- sim$counts

  coldata <- data.frame(
    condition = factor(
      c(rep("untrt", N_SAMPLES_PER_GROUP), rep("trt", N_SAMPLES_PER_GROUP)),
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
    stop("Fewer than five nonzero features after DESeq2 filtering.", call. = FALSE)
  }

  dds <- DESeq(dds, betaPrior = FALSE, quiet = TRUE)
  coef_name <- condition_coef_name(dds)

  standard <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    alpha = ALPHA_STANDARD
  )

  greater_abs <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = LFC_BOUNDARY,
    altHypothesis = "greaterAbs",
    alpha = ALPHA_STRONG
  )

  less_abs <- results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = LFC_BOUNDARY,
    altHypothesis = "lessAbs",
    alpha = ALPHA_WEAK
  )

  shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)

  df <- as.data.frame(standard)
  df$feature_id <- rownames(df)
  df$lfc_shrunk <- as.data.frame(shrunk)$log2FoldChange
  df$greaterAbs_padj <- as.data.frame(greater_abs)$padj
  df$lessAbs_padj <- as.data.frame(less_abs)$padj

  df <- dplyr::left_join(df, sim$truth, by = "feature_id")

  df$is_de[is.na(df$is_de)] <- FALSE
  df$true_weak[is.na(df$true_weak)] <- FALSE
  df$true_strong[is.na(df$true_strong)] <- FALSE

  inflated_wald <- inflate_null_wald(
    wald = df$stat,
    is_de = df$is_de,
    fraction = null_inflation
  )

  empirical <- fit_empirical_null(inflated_wald)

  df$empirical_p <- empirical$empirical_p
  df$empirical_bh <- empirical$empirical_bh

  hc_p <- hc_threshold(df$empirical_p)
  hbfss_cutoff <- if (is.na(hc_p)) NA_real_ else -log10(hc_p) * LFC_BOUNDARY

  df$neglog10_empirical_p <- neglog10_safe(df$empirical_p, floor_value = PROB_FLOOR)

  df$HBFSS <- abs(df$lfc_shrunk) * df$neglog10_empirical_p

  df$hc_pass <- !is.na(hc_p) &
    !is.na(df$empirical_p) &
    is.finite(df$empirical_p) &
    df$empirical_p <= hc_p

  df$standard_sig <- !is.na(df$padj) &
    df$padj < ALPHA_STANDARD &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= LFC_BOUNDARY

  df$empirical_bh_sig <- !is.na(df$empirical_bh) &
    df$empirical_bh < ALPHA_STANDARD &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= LFC_BOUNDARY

  df$greaterAbs_sig <- !is.na(df$greaterAbs_padj) &
    df$greaterAbs_padj < ALPHA_STRONG &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) >= LFC_BOUNDARY

  df$lessAbs_sig <- !is.na(df$lessAbs_padj) &
    df$lessAbs_padj < ALPHA_WEAK &
    !is.na(df$lfc_shrunk) &
    abs(df$lfc_shrunk) < LFC_BOUNDARY

  df$hbfss_raw_sig <- !is.na(hbfss_cutoff) &
    !is.na(df$HBFSS) &
    is.finite(df$HBFSS) &
    df$HBFSS >= hbfss_cutoff &
    df$hc_pass

  df$hbfss_total_sig <- df$standard_sig | df$hbfss_raw_sig

  df$hbfss_weak_region_sig <- df$lessAbs_sig & df$hbfss_raw_sig

  list(
    results = df,
    hc_p = hc_p,
    hbfss_cutoff = hbfss_cutoff,
    empirical_method = empirical$method
  )
}

# =============================================================================
# METRICS
# =============================================================================

metric_values <- function(predicted, truth) {
  predicted <- !is.na(predicted) & predicted
  truth <- !is.na(truth) & truth

  tp <- sum(predicted & truth)
  fp <- sum(predicted & !truth)
  fn <- sum(!predicted & truth)
  tn <- sum(!predicted & !truth)

  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  fdr <- if ((tp + fp) == 0) NA_real_ else fp / (tp + fp)
  specificity <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
    NA_real_
  } else {
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
    truth_count = sum(truth),
    stringsAsFactors = FALSE
  )
}

metric_row <- function(template, method, target, predicted, truth, out) {
  cbind(
    template,
    data.frame(
      method = method,
      truth_target = target,
      hc_p_threshold = out$hc_p,
      hbfss_cutoff = out$hbfss_cutoff,
      hc_valid = is.finite(out$hc_p) & !is.na(out$hc_p),
      empirical_method = out$empirical_method,
      stringsAsFactors = FALSE
    ),
    metric_values(predicted, truth)
  )
}

collect_metrics <- function(out, template) {
  df <- out$results

  dplyr::bind_rows(
    metric_row(template, "DESeq2_BH", "all_de", df$standard_sig, df$is_de, out),
    metric_row(template, "Empirical_BH", "all_de", df$empirical_bh_sig, df$is_de, out),
    metric_row(template, "HBFSS_total", "all_de", df$hbfss_total_sig, df$is_de, out),
    metric_row(template, "HBFSS_raw", "all_de", df$hbfss_raw_sig, df$is_de, out),
    metric_row(template, "GreaterAbs", "strong_de", df$greaterAbs_sig, df$true_strong, out),
    metric_row(template, "HBFSS_raw", "strong_de", df$hbfss_raw_sig, df$true_strong, out),
    metric_row(template, "LessAbs", "weak_de", df$lessAbs_sig, df$true_weak, out),
    metric_row(template, "HBFSS_weak_region", "weak_de", df$hbfss_weak_region_sig, df$true_weak, out)
  )
}

# =============================================================================
# FIGURES
# =============================================================================

METHOD_COLORS <- c(
  DESeq2_BH = "#999999",
  Empirical_BH = "#1F78B4",
  HBFSS_total = "#6A3D9A",
  HBFSS_raw = "#54278F",
  GreaterAbs = "#E31A1C",
  LessAbs = "#0072B2",
  HBFSS_weak_region = "#0072B2"
)

plot_theme <- function() {
  theme_bw(base_size = BASE_THEME_SIZE) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = BASE_THEME_SIZE + 1),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey88"),
      strip.text = element_text(face = "bold", size = BASE_THEME_SIZE - 0.5)
    )
}

metric_boxplot <- function(metrics, metric_name, title, y_label, methods, target, add_alpha_line = FALSE, drop_zero_truth = FALSE) {
  plot_df <- metrics[
    as.character(metrics$method) %in% methods &
      as.character(metrics$truth_target) == target,
    ,
    drop = FALSE
  ]

  if (isTRUE(drop_zero_truth)) {
    plot_df <- plot_df[!is.na(plot_df$truth_count) & plot_df$truth_count > 0, , drop = FALSE]
  }

  if (nrow(plot_df) == 0L) return(NULL)

  plot_df$method <- factor(as.character(plot_df$method), levels = methods)
  plot_df$value <- plot_df[[metric_name]]

  p <- ggplot(plot_df, aes(x = method, y = value, fill = method)) +
    geom_boxplot(outlier.size = 0.4, width = 0.62, linewidth = 0.25, na.rm = TRUE) +
    facet_grid(inflation_label + de_label ~ lfc_label) +
    scale_fill_manual(values = METHOD_COLORS[methods], breaks = methods, drop = FALSE) +
    labs(title = title, x = NULL, y = y_label) +
    plot_theme() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  if (isTRUE(add_alpha_line)) {
    p <- p + geom_hline(yintercept = ALPHA_STANDARD, linetype = "dashed", linewidth = 0.30)
  }

  p
}

hc_valid_plot <- function(qc) {
  if (nrow(qc) == 0L) return(NULL)

  ggplot(qc, aes(x = de_label, y = hc_valid_fraction, fill = de_label)) +
    geom_col(width = 0.62, linewidth = 0.25) +
    facet_grid(inflation_label ~ lfc_label) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    labs(
      title = "Higher-criticism valid-threshold fraction",
      x = "Differential-expression fraction",
      y = "Valid HC threshold fraction"
    ) +
    plot_theme() +
    theme(legend.position = "none")
}

hbfss_cutoff_plot <- function(qc) {
  plot_df <- qc[is.finite(qc$hbfss_cutoff_median) & !is.na(qc$hbfss_cutoff_median), , drop = FALSE]
  if (nrow(plot_df) == 0L) return(NULL)

  ggplot(plot_df, aes(x = de_label, y = hbfss_cutoff_median, fill = de_label)) +
    geom_col(width = 0.62, linewidth = 0.25) +
    facet_grid(inflation_label ~ lfc_label) +
    labs(
      title = "Median HBFSS cutoff across simulation conditions",
      x = "Differential-expression fraction",
      y = "Median cutoff"
    ) +
    plot_theme() +
    theme(legend.position = "none")
}

save_plot <- function(plot_obj, name, width = 13, height = 8) {
  if (is.null(plot_obj)) {
    stop("Cannot save missing plot: ", name, call. = FALSE)
  }

  root_png <- file.path(SIM_DIR, paste0(name, ".png"))
  archive_png <- file.path(PNG_DIR, paste0(name, ".png"))
  archive_pdf <- file.path(PDF_DIR, paste0(name, ".pdf"))

  ggsave(
    filename = archive_png,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = FIGURE_DPI,
    bg = "white",
    limitsize = FALSE
  )

  file.copy(archive_png, root_png, overwrite = TRUE)

  if (isTRUE(WRITE_PDF)) {
    ggsave(
      filename = archive_pdf,
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      bg = "white",
      limitsize = FALSE
    )
  }

  for (path in c(root_png, archive_png, if (isTRUE(WRITE_PDF)) archive_pdf else character(0))) {
    if (!file.exists(path)) stop("Figure was not written: ", path, call. = FALSE)
    size <- file.info(path)$size
    if (!is.finite(size) || is.na(size) || size <= 0) stop("Figure is empty: ", path, call. = FALSE)
  }

  data.frame(
    figure = name,
    png = repo_relative(root_png),
    png_archive = repo_relative(archive_png),
    pdf = if (isTRUE(WRITE_PDF)) repo_relative(archive_pdf) else NA_character_,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# OUTPUT TEXT
# =============================================================================

write_methods_note <- function() {
  lines <- c(
    "HBFSS simulation methods",
    "",
    paste0("Negative-binomial count matrices were simulated with ", N_FEATURES, " features and ", N_SAMPLES_PER_GROUP, " samples per condition."),
    paste0("Differential-expression fractions were ", paste(DE_FRACTIONS, collapse = ", "), "."),
    paste0("Fixed absolute log2 fold-change magnitudes were ", paste(FIXED_LFC_MAGNITUDES, collapse = ", "), "."),
    paste0("The weak-mixture condition sampled true absolute log2 fold changes uniformly from ", WEAK_LFC_MIN, " to ", WEAK_LFC_MAX, "."),
    paste0("Null Wald-statistic inflation fractions were ", paste(NULL_INFLATION_FRACTIONS, collapse = ", "), "."),
    paste0("Each grid condition used ", N_REPS, " replicate simulations."),
    "",
    "Differential features were randomly selected in each replicate.",
    "Null inflation was applied only to null features.",
    "DESeq2 was run with design ~ condition.",
    "Log2 fold changes were shrunken using standard apeglm through DESeq2::lfcShrink(type = 'apeglm').",
    "No apeMethod shortcut was used.",
    "",
    "Empirical-null p-values were estimated from DESeq2 Wald statistics using fdrtool.",
    "Higher criticism thresholds were estimated with fdrtool::hc.thresh.",
    paste0("HC thresholds below ", HC_INVALID_BELOW, " or at/above ", HC_INVALID_AT_OR_ABOVE, " were marked invalid."),
    "",
    "HBFSS was defined as abs(apeglm-shrunken log2 fold change) * -log10(empirical p).",
    "The HBFSS cutoff was defined as -log10(HC p-threshold) * LFC boundary.",
    "No artificial HBFSS floor or alternate cutoff was added.",
    "",
    "The simulation reports F1, precision, recall, observed FDR, discovery count, weak-effect recovery, strong-effect recovery, HC valid-threshold fraction, and HBFSS cutoff stability."
  )

  writeLines(lines, file.path(SIM_DIR, "Methods_Simulation.txt"))
}

format_number <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(x) || is.na(x)) return("NA")
  formatC(x, format = "f", digits = digits)
}

write_interpretation <- function(summary_table, qc_table, failures) {
  overall <- summary_table %>%
    group_by(method, truth_target) %>%
    summarise(
      precision_mean = mean_finite(precision_mean),
      recall_mean = mean_finite(recall_mean),
      f1_mean = mean_finite(f1_mean),
      fdr_mean = mean_finite(fdr_mean),
      discovery_count_mean = mean_finite(discovery_count_mean),
      truth_count_mean = mean_finite(truth_count_mean),
      .groups = "drop"
    ) %>%
    arrange(truth_target, desc(f1_mean), method)

  write_csv(overall, file.path(SIM_DIR, "Simulation_Overall_MethodMeans.csv"))

  best_all <- overall[overall$truth_target == "all_de" & is.finite(overall$f1_mean), , drop = FALSE]
  best_method <- if (nrow(best_all) > 0L) as.character(best_all$method[which.max(best_all$f1_mean)]) else "NA"

  hbfss_total <- overall[overall$method == "HBFSS_total" & overall$truth_target == "all_de", , drop = FALSE]
  deseq2 <- overall[overall$method == "DESeq2_BH" & overall$truth_target == "all_de", , drop = FALSE]

  delta_f1 <- if (nrow(hbfss_total) && nrow(deseq2)) hbfss_total$f1_mean[1] - deseq2$f1_mean[1] else NA_real_
  delta_recall <- if (nrow(hbfss_total) && nrow(deseq2)) hbfss_total$recall_mean[1] - deseq2$recall_mean[1] else NA_real_
  delta_fdr <- if (nrow(hbfss_total) && nrow(deseq2)) hbfss_total$fdr_mean[1] - deseq2$fdr_mean[1] else NA_real_

  hc_mean <- mean_finite(qc_table$hc_valid_fraction)

  lines <- c(
    "HBFSS simulation interpretation",
    "",
    "Metric definitions:",
    "F1 balances precision and recall.",
    "Precision is the fraction of called features that are true positives.",
    "Recall is the fraction of true positives recovered.",
    "Observed FDR is the fraction of called features that are false positives.",
    "Discovery count shows how aggressive each method is.",
    "HC valid fraction shows how often higher criticism produced a usable threshold.",
    "",
    paste0("Best mean all-DE F1 method: ", best_method),
    paste0("HBFSS_total minus DESeq2_BH mean F1: ", format_number(delta_f1)),
    paste0("HBFSS_total minus DESeq2_BH mean recall: ", format_number(delta_recall)),
    paste0("HBFSS_total minus DESeq2_BH mean observed FDR: ", format_number(delta_fdr)),
    paste0("Mean HC valid-threshold fraction: ", format_number(hc_mean)),
    paste0("Failed replicate count: ", nrow(failures)),
    "",
    "Interpretation guide:",
    "If HBFSS improves recall with stable or improved FDR, it supports the claim that HBFSS recovers additional true positives without excessive false discovery burden.",
    "If HBFSS improves recall but also increases FDR, it should be described as a sensitivity-oriented tradeoff.",
    "If HC valid fraction is low in any condition, those conditions should be interpreted cautiously because the HBFSS boundary depends on a valid HC threshold.",
    "Weak-effect plots evaluate sub-boundary true effects.",
    "Strong-effect plots evaluate above-boundary true effects.",
    "",
    "Primary files:",
    "Simulation_RunMetrics_Long.csv",
    "Simulation_MethodSummary.csv",
    "Simulation_HC_QC_Summary.csv",
    "Simulation_Overall_MethodMeans.csv",
    "Simulation_F1_Boxplot.png",
    "Simulation_Precision_Boxplot.png",
    "Simulation_Recall_Boxplot.png",
    "Simulation_FDR_Boxplot.png",
    "Simulation_DiscoveryCount_Boxplot.png",
    "Simulation_WeakEffect_F1_Boxplot.png",
    "Simulation_StrongEffect_F1_Boxplot.png",
    "Simulation_HC_ValidFraction.png",
    "Simulation_HBFSS_Cutoff.png"
  )

  writeLines(lines, file.path(SIM_DIR, "Interpretation_Simulation.txt"))
  writeLines(lines, file.path(SIM_DIR, "Interpretation_Simulation.md"))
}

write_manifest <- function() {
  files <- list.files(SIM_DIR, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files)]

  manifest <- data.frame(
    file = vapply(files, repo_relative, character(1)),
    size_bytes = file.info(files)$size,
    stringsAsFactors = FALSE
  )

  write_csv(manifest, file.path(SIM_DIR, "Manifest_Simulation.csv"))
  manifest
}

# =============================================================================
# GIT
# =============================================================================

git_available <- function() {
  nzchar(Sys.which("git"))
}

git_run <- function(args, allow_failure = FALSE) {
  if (!git_available()) stop("Git is not available on PATH.", call. = FALSE)

  out <- suppressWarnings(system2(
    "git",
    args = c("-C", REPO_ROOT, args),
    stdout = TRUE,
    stderr = TRUE,
    timeout = 180
  ))

  status <- attr(out, "status")
  if (is.null(status)) status <- 0L

  if (length(out) > 0L) log_msg(paste(out, collapse = "\n"))

  if (status != 0L && !allow_failure) {
    stop("Git command failed: git ", paste(args, collapse = " "), "\n", paste(out, collapse = "\n"), call. = FALSE)
  }

  list(status = status, output = out)
}

git_first_line <- function(args, allow_failure = FALSE) {
  result <- git_run(args, allow_failure = allow_failure)
  if (result$status != 0L || length(result$output) == 0L) return(NA_character_)
  trimws(result$output[1])
}

push_outputs <- function() {
  if (!isTRUE(GIT_PUSH_AFTER_SUCCESS)) {
    log_msg("Git push skipped by setting.")
    return(invisible(FALSE))
  }

  if (!dir.exists(file.path(REPO_ROOT, ".git"))) {
    log_msg("Git push skipped because repository root has no .git directory.")
    return(invisible(FALSE))
  }

  remote_url <- git_first_line(c("remote", "get-url", GIT_REMOTE), allow_failure = TRUE)
  if (is.na(remote_url) || !nzchar(remote_url)) {
    stop("Git remote is not configured: ", GIT_REMOTE, call. = FALSE)
  }

  current_branch <- git_first_line(c("rev-parse", "--abbrev-ref", "HEAD"))
  if (is.na(current_branch) || !nzchar(current_branch) || identical(current_branch, "HEAD")) {
    stop("Git is in detached HEAD state.", call. = FALSE)
  }

  target_branch <- if (!is.na(GIT_BRANCH) && nzchar(GIT_BRANCH)) GIT_BRANCH else current_branch

  script_file <- current_script()
  stage_paths <- c(repo_relative(SIM_DIR))

  if (!is.na(script_file) && file.exists(script_file)) {
    stage_paths <- c(repo_relative(script_file), stage_paths)
  }

  git_run(c("fetch", GIT_REMOTE), allow_failure = TRUE)
  git_run(c("add", "--", unique(stage_paths)))

  changed <- git_run(c("diff", "--cached", "--quiet"), allow_failure = TRUE)
  if (changed$status == 0L) {
    log_msg("No Git changes to commit.")
    return(invisible(FALSE))
  }

  git_run(c("status", "--short"))
  git_run(c("commit", "-m", GIT_COMMIT_MESSAGE))
  git_run(c("push", "-u", GIT_REMOTE, paste0("HEAD:", target_branch)))

  log_msg("Git push complete.")
  invisible(TRUE)
}

# =============================================================================
# MAIN RUN
# =============================================================================

run_simulation <- function() {
  set.seed(SEED)

  grid <- make_simulation_grid()
  total_runs <- nrow(grid) * N_REPS

  log_msg("Simulation started.")
  log_msg("Repository root: ", REPO_ROOT)
  log_msg("Output directory: ", SIM_DIR)
  log_msg("Features: ", N_FEATURES)
  log_msg("Samples per group: ", N_SAMPLES_PER_GROUP)
  log_msg("Replicates per condition: ", N_REPS)
  log_msg("Total planned runs: ", total_runs)

  write_status("running", paste0("0/", total_runs))

  metric_blocks <- list()
  failure_blocks <- list()
  run_counter <- 0L
  start_time <- Sys.time()

  for (grid_i in seq_len(nrow(grid))) {
    grid_row <- grid[grid_i, , drop = FALSE]

    for (rep_i in seq_len(N_REPS)) {
      run_counter <- run_counter + 1L

      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      eta <- if (run_counter > 0) elapsed * (total_runs - run_counter) / run_counter else NA_real_

      log_msg(
        "Run ", run_counter, "/", total_runs,
        " | profile=", grid_row$profile,
        " | DE=", grid_row$de_fraction,
        " | LFC=", ifelse(is.na(grid_row$lfc_magnitude), "weak_mix", grid_row$lfc_magnitude),
        " | inflation=", grid_row$null_inflation,
        " | replicate=", rep_i,
        " | elapsed_min=", signif(elapsed, 3),
        " | eta_min=", signif(eta, 3)
      )

      write_status("running", paste0(run_counter, "/", total_runs))

      template <- data.frame(
        de_fraction = grid_row$de_fraction,
        lfc_magnitude = grid_row$lfc_magnitude,
        profile = grid_row$profile,
        null_inflation = grid_row$null_inflation,
        replicate = rep_i,
        de_label = paste0(grid_row$de_fraction * 100, "% DE"),
        lfc_label = lfc_label_from_profile(grid_row$profile, grid_row$lfc_magnitude),
        inflation_label = ifelse(grid_row$null_inflation == 0, "No null inflation", paste0(grid_row$null_inflation * 100, "% null inflation")),
        stringsAsFactors = FALSE
      )

      result <- tryCatch({
        sim <- simulate_counts(
          n_features = N_FEATURES,
          n_samples_per_group = N_SAMPLES_PER_GROUP,
          de_fraction = grid_row$de_fraction,
          lfc_magnitude = grid_row$lfc_magnitude,
          profile = grid_row$profile
        )

        out <- run_one_analysis(sim, null_inflation = grid_row$null_inflation)
        collect_metrics(out, template)
      }, error = function(e) {
        failure_blocks[[length(failure_blocks) + 1L]] <<- cbind(
          template,
          data.frame(error_message = conditionMessage(e), stringsAsFactors = FALSE)
        )
        NULL
      })

      if (!is.null(result)) {
        metric_blocks[[length(metric_blocks) + 1L]] <- result
      }

      if (length(metric_blocks) > 0L && (run_counter %% 10L == 0L || run_counter == total_runs)) {
        write_csv(dplyr::bind_rows(metric_blocks), file.path(SIM_DIR, "Simulation_RunMetrics_Long.partial.csv"))
      }
    }
  }

  failures <- if (length(failure_blocks) > 0L) dplyr::bind_rows(failure_blocks) else data.frame()
  write_csv(failures, file.path(SIM_DIR, "Simulation_Failures.csv"))

  metrics <- if (length(metric_blocks) > 0L) dplyr::bind_rows(metric_blocks) else data.frame()

  if (nrow(metrics) == 0L) {
    stop("No successful simulation replicates. See Simulation_Failures.csv.", call. = FALSE)
  }

  failure_fraction <- nrow(failures) / total_runs

  metrics$method <- factor(
    as.character(metrics$method),
    levels = c("DESeq2_BH", "Empirical_BH", "HBFSS_total", "HBFSS_raw", "GreaterAbs", "LessAbs", "HBFSS_weak_region")
  )

  metrics$truth_target <- factor(
    as.character(metrics$truth_target),
    levels = c("all_de", "strong_de", "weak_de")
  )

  write_csv(metrics, file.path(SIM_DIR, "Simulation_RunMetrics_Long.csv"))

  summary_table <- metrics %>%
    group_by(
      de_fraction,
      lfc_magnitude,
      profile,
      null_inflation,
      de_label,
      lfc_label,
      inflation_label,
      method,
      truth_target,
      empirical_method
    ) %>%
    summarise(
      n_rows = n(),
      precision_mean = mean_finite(precision),
      precision_sd = sd(precision, na.rm = TRUE),
      recall_mean = mean_finite(recall),
      recall_sd = sd(recall, na.rm = TRUE),
      f1_mean = mean_finite(f1),
      f1_sd = sd(f1, na.rm = TRUE),
      fdr_mean = mean_finite(fdr),
      fdr_sd = sd(fdr, na.rm = TRUE),
      discovery_count_mean = mean_finite(discovery_count),
      truth_count_mean = mean_finite(truth_count),
      hc_valid_fraction = mean(as.numeric(hc_valid), na.rm = TRUE),
      hc_p_threshold_median = median_finite(hc_p_threshold),
      hbfss_cutoff_median = median_finite(hbfss_cutoff),
      .groups = "drop"
    )

  write_csv(summary_table, file.path(SIM_DIR, "Simulation_MethodSummary.csv"))

  qc_table <- metrics %>%
    group_by(
      de_fraction,
      lfc_magnitude,
      profile,
      null_inflation,
      de_label,
      lfc_label,
      inflation_label,
      replicate,
      empirical_method
    ) %>%
    summarise(
      hc_p_threshold = first(hc_p_threshold),
      hbfss_cutoff = first(hbfss_cutoff),
      hc_valid = first(hc_valid),
      .groups = "drop"
    ) %>%
    group_by(
      de_fraction,
      lfc_magnitude,
      profile,
      null_inflation,
      de_label,
      lfc_label,
      inflation_label,
      empirical_method
    ) %>%
    summarise(
      n_replicates = n(),
      hc_valid_fraction = mean(as.numeric(hc_valid), na.rm = TRUE),
      hc_p_threshold_min = min_finite(ifelse(hc_valid, hc_p_threshold, NA_real_)),
      hc_p_threshold_median = median_finite(ifelse(hc_valid, hc_p_threshold, NA_real_)),
      hc_p_threshold_max = max_finite(ifelse(hc_valid, hc_p_threshold, NA_real_)),
      hbfss_cutoff_median = median_finite(ifelse(hc_valid, hbfss_cutoff, NA_real_)),
      .groups = "drop"
    )

  write_csv(qc_table, file.path(SIM_DIR, "Simulation_HC_QC_Summary.csv"))

  all_methods <- c("DESeq2_BH", "Empirical_BH", "HBFSS_total", "HBFSS_raw")
  strong_methods <- c("GreaterAbs", "HBFSS_raw")
  weak_methods <- c("LessAbs", "HBFSS_weak_region")

  figure_table <- dplyr::bind_rows(
    save_plot(metric_boxplot(metrics, "f1", "Simulation F1 score", "F1", all_methods, "all_de"), "Simulation_F1_Boxplot"),
    save_plot(metric_boxplot(metrics, "precision", "Simulation precision", "Precision", all_methods, "all_de"), "Simulation_Precision_Boxplot"),
    save_plot(metric_boxplot(metrics, "recall", "Simulation recall", "Recall", all_methods, "all_de"), "Simulation_Recall_Boxplot"),
    save_plot(metric_boxplot(metrics, "fdr", "Simulation observed FDR", "Observed FDR", all_methods, "all_de", add_alpha_line = TRUE), "Simulation_FDR_Boxplot"),
    save_plot(metric_boxplot(metrics, "discovery_count", "Simulation discovery count", "Discovery count", all_methods, "all_de"), "Simulation_DiscoveryCount_Boxplot"),
    save_plot(metric_boxplot(metrics, "f1", "Weak-effect F1 score", "F1", weak_methods, "weak_de"), "Simulation_WeakEffect_F1_Boxplot"),
    save_plot(metric_boxplot(metrics, "precision", "Weak-effect precision", "Precision", weak_methods, "weak_de"), "Simulation_WeakEffect_Precision_Boxplot"),
    save_plot(metric_boxplot(metrics, "recall", "Weak-effect recall", "Recall", weak_methods, "weak_de"), "Simulation_WeakEffect_Recall_Boxplot"),
    save_plot(metric_boxplot(metrics, "f1", "Strong-effect F1 score", "F1", strong_methods, "strong_de", drop_zero_truth = TRUE), "Simulation_StrongEffect_F1_Boxplot"),
    save_plot(metric_boxplot(metrics, "precision", "Strong-effect precision", "Precision", strong_methods, "strong_de", drop_zero_truth = TRUE), "Simulation_StrongEffect_Precision_Boxplot"),
    save_plot(metric_boxplot(metrics, "recall", "Strong-effect recall", "Recall", strong_methods, "strong_de", drop_zero_truth = TRUE), "Simulation_StrongEffect_Recall_Boxplot"),
    save_plot(metric_boxplot(metrics, "fdr", "Strong-effect observed FDR", "Observed FDR", strong_methods, "strong_de", add_alpha_line = TRUE, drop_zero_truth = TRUE), "Simulation_StrongEffect_FDR_Boxplot"),
    save_plot(hc_valid_plot(qc_table), "Simulation_HC_ValidFraction", width = 12, height = 7),
    save_plot(hbfss_cutoff_plot(qc_table), "Simulation_HBFSS_Cutoff", width = 12, height = 7)
  )

  write_csv(figure_table, file.path(SIM_DIR, "Figure_Export_Status.csv"))

  write_methods_note()
  write_interpretation(summary_table, qc_table, failures)
  writeLines(capture.output(sessionInfo()), file.path(SIM_DIR, "SessionInfo_Simulation.txt"))

  if (file.exists(file.path(SIM_DIR, "Simulation_RunMetrics_Long.partial.csv"))) {
    unlink(file.path(SIM_DIR, "Simulation_RunMetrics_Long.partial.csv"), force = TRUE)
  }

  manifest <- write_manifest()

  if (failure_fraction > MAX_FAILURE_FRACTION) {
    stop(
      "Failure fraction exceeded allowed maximum. Failures: ",
      nrow(failures),
      "/",
      total_runs,
      " = ",
      round(failure_fraction, 4),
      ". Figures were written but Git push was stopped.",
      call. = FALSE
    )
  }

  list(
    metrics = metrics,
    summary = summary_table,
    qc = qc_table,
    figures = figure_table,
    failures = failures,
    manifest = manifest
  )
}

main <- function() {
  write_status("starting", "initializing")

  result <- run_simulation()

  write_status("figures_complete", "simulation outputs written")

  push_outputs()

  write_status("complete", "simulation finished")

  cat("\n")
  cat("=====================================================\n")
  cat("HBFSS simulation complete\n")
  cat("Output directory: ", SIM_DIR, "\n", sep = "")
  cat("Metric rows: ", nrow(result$metrics), "\n", sep = "")
  cat("Figures: ", nrow(result$figures), "\n", sep = "")
  cat("Failures: ", nrow(result$failures), "\n", sep = "")
  cat("=====================================================\n\n")

  print(result$summary)

  invisible(result)
}

tryCatch(
  main(),
  error = function(e) {
    write_status("failed", conditionMessage(e))
    log_msg("FAILED: ", conditionMessage(e))
    try(writeLines(capture.output(sessionInfo()), file.path(SIM_DIR, "SessionInfo_Simulation_FAILED.txt")), silent = TRUE)
    if (!interactive()) quit(save = "no", status = 1L, runLast = FALSE)
    stop(e)
  }
)
