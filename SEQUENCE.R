#!/usr/bin/env Rscript

# =============================================================================
# SEQUENCE / HBFSS SIMULATION VALIDATION - METHOD-ALIGNED FINAL SCRIPT
#
# Purpose:
#   Standalone simulation comparing DESeq2 Wald/BH and DESeq2 composite-null
#   approaches against empirical-null / higher-criticism HBFSS calls.
#
# Key definitions preserved:
#   DESeq2 uses Wald p-values and DESeq2's intrinsic BH-adjusted p-values.
#   DESeq2 strong/weak arms use DESeq2 composite null tests:
#       greaterAbs for strong effects
#       lessAbs for sub-boundary effects
#   fdrtool estimates empirical p-values from the Wald statistics.
#   higher criticism chooses one empirical-p threshold per simulated dataset.
#   HBFSS = abs(apeglm-shrunken log2FC) * -log10(empirical p).
#   HBFSS cutoff = -log10(HC p-threshold) * LFC boundary.
#
# Main outputs are direct, interpretable figures rather than large walls of
# boxplots. The script is self-contained and does not depend on the main
# manuscript pipeline.
# =============================================================================

options(stringsAsFactors = FALSE)
options(width = 140)

# =============================================================================
# SETTINGS
# =============================================================================

output_folder_name <- "manuscript_final_clean"
simulation_folder_name <- "simulation_validation"
reset_simulation_dir <- TRUE

alpha_standard <- 0.10
alpha_strong <- 0.10
alpha_weak <- 0.10
lfc_boundary <- 1.0

hc_invalid_below <- 0.001
hc_invalid_at_or_above <- 0.95
calculation_p_floor <- .Machine$double.xmin
plot_p_floor <- 1e-16
strict_empirical_null <- FALSE

simulation_seed <- as.integer(Sys.getenv("SEQUENCE_SIM_SEED", "42"))
simulation_n_features <- as.integer(Sys.getenv("SEQUENCE_SIM_FEATURES", "1000"))
simulation_n_samples_per_group <- as.integer(Sys.getenv("SEQUENCE_SIM_SAMPLES_PER_GROUP", "6"))
simulation_n_reps <- as.integer(Sys.getenv("SEQUENCE_SIM_REPS", "20"))
simulation_progress_every <- as.integer(Sys.getenv("SEQUENCE_SIM_PROGRESS_EVERY", "5"))
simulation_checkpoint_every <- as.integer(Sys.getenv("SEQUENCE_SIM_CHECKPOINT_EVERY", "10"))

simulation_base_mean <- 200
simulation_dispersion_null <- 0.10
simulation_dispersion_de <- 0.15
simulation_de_fractions <- c(0.05, 0.10)
simulation_lfc_magnitudes <- c(0.50, 1.00, 2.00)
simulation_weak_lfc_min <- 0.20
simulation_weak_lfc_max <- 0.80
simulation_null_inflation <- c(0.00, 0.10)
simulation_fit_type <- "local"

figure_dpi <- as.integer(Sys.getenv("SEQUENCE_SIM_DPI", "600"))
figure_width_twocol <- 6.7
base_theme_size <- 7
export_pdf_also <- TRUE
export_tiff_also <- TRUE
stamp_figures <- FALSE

# Git is optional. Default FALSE to avoid turning a successful simulation into a
# failed run because of credentials/network state. Set SEQUENCE_GIT_PUSH=true to push.
git_push_after_success <- tolower(Sys.getenv("SEQUENCE_GIT_PUSH", "false")) %in% c("1", "true", "yes", "y")
git_remote_name <- Sys.getenv("SEQUENCE_GIT_REMOTE", "origin")
git_branch_name <- Sys.getenv("SEQUENCE_GIT_BRANCH", NA_character_)
git_commit_message <- Sys.getenv(
  "SEQUENCE_GIT_COMMIT_MESSAGE",
  paste0("Update SEQUENCE simulation validation outputs ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)

# =============================================================================
# PACKAGES
# =============================================================================

required_packages <- c("DESeq2", "apeglm", "fdrtool", "ggplot2", "dplyr", "scales", "grid")
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
  library(dplyr)
  library(scales)
  library(grid)
})

# =============================================================================
# PATHS / LOGGING
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
  for (start in starts) {
    current <- normalizePath(start, winslash = "/", mustWork = TRUE)
    repeat {
      if (dir.exists(file.path(current, ".git"))) return(current)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

repo_root <- find_repo_root()
output_dir <- file.path(repo_root, "exports", output_folder_name)
simulation_dir <- file.path(output_dir, simulation_folder_name)
figures_dir <- file.path(simulation_dir, "figures")
figures_png_dir <- file.path(figures_dir, "png")
figures_pdf_dir <- file.path(figures_dir, "pdf")
figures_tiff_dir <- file.path(figures_dir, "tiff")
log_file <- file.path(simulation_dir, "Simulation_Log.txt")
status_file <- file.path(simulation_dir, "Simulation_Status.txt")

if (isTRUE(reset_simulation_dir) && dir.exists(simulation_dir)) {
  unlink(simulation_dir, recursive = TRUE, force = TRUE)
}
for (d in c(simulation_dir, figures_png_dir, figures_pdf_dir, figures_tiff_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_message <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(...))
  message(line)
  cat(line, "\n", file = log_file, append = TRUE)
}

write_status <- function(status, detail = "") {
  writeLines(
    c(
      paste0("status=", status),
      paste0("time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("detail=", detail)
    ),
    status_file
  )
}

repo_relative <- function(path) {
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  full <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  if (startsWith(full, prefix)) return(substr(full, nchar(prefix) + 1L, nchar(full)))
  full
}

save_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE)
  invisible(path)
}

cleanup_rplots_pdf <- function() {
  stray <- file.path(repo_root, "Rplots.pdf")
  if (file.exists(stray)) unlink(stray, force = TRUE)
}

positive_integer <- function(x, default, label) {
  y <- suppressWarnings(as.integer(x[1]))
  if (!is.finite(y) || is.na(y) || y < 1L) {
    warning(label, " was invalid; using ", default, ".", call. = FALSE)
    return(as.integer(default))
  }
  as.integer(y)
}

simulation_seed <- positive_integer(simulation_seed, 42L, "simulation_seed")
simulation_n_features <- positive_integer(simulation_n_features, 1000L, "simulation_n_features")
simulation_n_samples_per_group <- positive_integer(simulation_n_samples_per_group, 6L, "simulation_n_samples_per_group")
simulation_n_reps <- positive_integer(simulation_n_reps, 20L, "simulation_n_reps")
simulation_progress_every <- positive_integer(simulation_progress_every, 5L, "simulation_progress_every")
simulation_checkpoint_every <- positive_integer(simulation_checkpoint_every, 10L, "simulation_checkpoint_every")
figure_dpi <- positive_integer(figure_dpi, 600L, "figure_dpi")

# =============================================================================
# NUMERIC HELPERS
# =============================================================================

finite_vals <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[is.finite(y) & !is.na(y)]
}

mean_finite <- function(x) {
  v <- finite_vals(x)
  if (!length(v)) NA_real_ else mean(v)
}

median_finite <- function(x) {
  v <- finite_vals(x)
  if (!length(v)) NA_real_ else stats::median(v)
}

min_finite <- function(x) {
  v <- finite_vals(x)
  if (!length(v)) NA_real_ else min(v)
}

max_finite <- function(x) {
  v <- finite_vals(x)
  if (!length(v)) NA_real_ else max(v)
}

quantile_finite <- function(x, prob) {
  v <- finite_vals(x)
  if (!length(v)) NA_real_ else as.numeric(stats::quantile(v, probs = prob, names = FALSE, type = 7, na.rm = TRUE))
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

as_integer_count_matrix <- function(x, label) {
  row_ids <- rownames(x)
  col_ids <- colnames(x)
  mat <- as.matrix(x)
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (any(is.na(mat) | !is.finite(mat))) stop(label, " has NA or non-finite count values.", call. = FALSE)
  if (any(mat < 0)) stop(label, " has negative count values.", call. = FALSE)
  rounded <- round(mat)
  storage.mode(rounded) <- "integer"
  rownames(rounded) <- row_ids
  colnames(rounded) <- col_ids
  rounded
}

# =============================================================================
# LABELS / THEME
# =============================================================================

method_label_map <- c(
  DESeq2_BH = "DESeq2 BH",
  DESeq2_Strong = "DESeq2 greaterAbs",
  DESeq2_Weak = "DESeq2 weak composite",
  Empirical_BH = "Empirical BH",
  HBFSS_All = "HBFSS all",
  HBFSS_Strong = "HBFSS strong",
  HBFSS_Weak = "HBFSS weak"
)

# Okabe-Ito / grayscale-conscious choices.
method_color_map <- c(
  DESeq2_BH = "#000000",
  DESeq2_Strong = "#D55E00",
  DESeq2_Weak = "#E69F00",
  Empirical_BH = "#0072B2",
  HBFSS_All = "#009E73",
  HBFSS_Strong = "#009E73",
  HBFSS_Weak = "#56B4E9"
)

lfc_label <- function(lfc_magnitude, profile) {
  if (identical(as.character(profile), "weak_mixture")) {
    paste0("Weak mix\n|LFC| ", simulation_weak_lfc_min, "-", simulation_weak_lfc_max)
  } else {
    paste0("Fixed\n|LFC| ", formatC(as.numeric(lfc_magnitude), format = "fg", digits = 3))
  }
}

inflation_label <- function(x) {
  ifelse(as.numeric(x) > 0, paste0(as.numeric(x) * 100, "% null inflation"), "No null inflation")
}

condition_label <- function(lfc, infl, de) {
  paste(gsub("\n", " ", lfc), infl, paste0(de * 100, "% DE"), sep = " | ")
}

figure_stamp <- function() {
  if (!isTRUE(stamp_figures)) return(NULL)
  sha <- tryCatch(
    system2("git", args = c("-C", repo_root, "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  sha <- if (length(sha) > 0L && !is.na(sha[1])) sha[1] else "NA"
  paste0("seed: ", simulation_seed, " | git: ", sha)
}

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title = element_text(face = "bold", size = 8.5, hjust = 0.5),
      plot.subtitle = element_text(size = 7, hjust = 0.5),
      plot.caption = element_text(size = 5.5, hjust = 1),
      axis.title = element_text(face = "bold", size = 7),
      axis.text = element_text(color = "black", size = 6),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.20, color = "grey88"),
      strip.text = element_text(face = "bold", size = 6.2),
      plot.margin = margin(4, 5, 4, 4)
    )
}

# =============================================================================
# FIGURE EXPORT
# =============================================================================

check_nonempty_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " was not written: ", path, call. = FALSE)
  file_size <- suppressWarnings(file.info(path)$size)
  if (!is.finite(file_size) || is.na(file_size) || file_size <= 0) {
    stop(label, " is empty: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

save_figure <- function(plot_obj, figure_name, width = figure_width_twocol, height = 4.8) {
  if (is.null(plot_obj) || !inherits(plot_obj, "ggplot")) {
    stop("Invalid plot object for ", figure_name, call. = FALSE)
  }

  root_png <- file.path(simulation_dir, paste0(figure_name, ".png"))
  png_path <- file.path(figures_png_dir, paste0(figure_name, ".png"))
  pdf_path <- file.path(figures_pdf_dir, paste0(figure_name, ".pdf"))
  tiff_path <- file.path(figures_tiff_dir, paste0(figure_name, ".tiff"))

  ggplot2::ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = figure_dpi,
    bg = "white",
    limitsize = FALSE
  )
  check_nonempty_file(png_path, "PNG figure")
  file.copy(png_path, root_png, overwrite = TRUE)
  check_nonempty_file(root_png, "Root PNG figure")

  if (isTRUE(export_pdf_also)) {
    ggplot2::ggsave(
      filename = pdf_path,
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      bg = "white",
      limitsize = FALSE
    )
    check_nonempty_file(pdf_path, "PDF figure")
  } else {
    pdf_path <- NA_character_
  }

  tiff_status <- "not_requested"
  tiff_error <- NA_character_
  if (isTRUE(export_tiff_also)) {
    tiff_status <- tryCatch({
      grDevices::tiff(filename = tiff_path, width = width, height = height, units = "in", res = figure_dpi, compression = "lzw")
      print(plot_obj)
      grDevices::dev.off()
      check_nonempty_file(tiff_path, "TIFF figure")
      "saved"
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      tiff_error <<- conditionMessage(e)
      log_message("TIFF export failed for ", figure_name, ": ", tiff_error)
      "failed"
    })
  } else {
    tiff_path <- NA_character_
  }

  data.frame(
    figure = figure_name,
    png_file = repo_relative(root_png),
    png_archive_file = repo_relative(png_path),
    pdf_file = if (isTRUE(export_pdf_also)) repo_relative(pdf_path) else NA_character_,
    tiff_file = if (isTRUE(export_tiff_also) && identical(tiff_status, "saved")) repo_relative(tiff_path) else NA_character_,
    tiff_status = tiff_status,
    tiff_error = tiff_error,
    status = "saved",
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

safe_save_figure <- function(plot_obj, figure_name, width = figure_width_twocol, height = 4.8) {
  tryCatch(
    save_figure(plot_obj, figure_name, width = width, height = height),
    error = function(e) {
      msg <- conditionMessage(e)
      log_message("Figure export failed for ", figure_name, ": ", msg)
      data.frame(
        figure = figure_name,
        png_file = NA_character_,
        png_archive_file = NA_character_,
        pdf_file = NA_character_,
        tiff_file = NA_character_,
        tiff_status = if (isTRUE(export_tiff_also)) "not_attempted" else "not_requested",
        tiff_error = NA_character_,
        status = "failed",
        error_message = msg,
        stringsAsFactors = FALSE
      )
    }
  )
}

# =============================================================================
# EMPIRICAL NULL / HC
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
      fdrtool::fdrtool(z, statistic = "normal", plot = FALSE, verbose = FALSE, cutoff.method = method, pct0 = 0.75),
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
    if (isTRUE(strict_empirical_null)) stop(label, ": empirical-null fit failed.", call. = FALSE)
    method_used <- "theoretical_normal_fallback"
    p_raw <- 2 * stats::pnorm(-abs(z))
    fit <- list(pval = p_raw, qval = p.adjust(p_raw, method = "BH"), lfdr = rep(NA_real_, length(z)))
  }

  empirical_p <- rep(NA_real_, length(statistic))
  empirical_p[valid] <- clip_probability(fit$pval)
  empirical_bh <- rep(NA_real_, length(statistic))
  ok <- !is.na(empirical_p) & is.finite(empirical_p)
  empirical_bh[ok] <- p.adjust(empirical_p[ok], method = "BH")

  list(empirical_p = empirical_p, empirical_bh = empirical_bh, empirical_null_method = method_used)
}

hc_threshold <- function(empirical_p) {
  p <- sort(clip_probability(empirical_p), decreasing = FALSE, na.last = NA)
  if (length(p) < 5L) return(NA_real_)
  threshold <- tryCatch(suppressWarnings(as.numeric(fdrtool::hc.thresh(p)[1])), error = function(e) NA_real_)
  if (!is.finite(threshold) || is.na(threshold) || threshold <= 0 || threshold >= 1) return(NA_real_)
  if (threshold < hc_invalid_below) return(NA_real_)
  if (threshold >= hc_invalid_at_or_above) return(NA_real_)
  threshold
}

# =============================================================================
# SIMULATION ENGINE
# =============================================================================

simulate_sequence_counts <- function(n_features, n_samples, base_mean, disp_null, de_fraction, lfc_magnitude, disp_de, lfc_profile = c("fixed", "weak_mixture")) {
  lfc_profile <- match.arg(lfc_profile)
  n_de <- max(2L, min(n_features - 2L, round(n_features * de_fraction)))
  feature_id <- paste0("sim_feature_", seq_len(n_features))
  de_id <- sample(seq_len(n_features), n_de, replace = FALSE)

  base_mu <- pmax(stats::rgamma(n_features, shape = 2.5, scale = base_mean / 2.5), 2)
  dispersion <- rep(disp_null, n_features)
  dispersion[de_id] <- disp_de
  dispersion <- pmin(pmax(dispersion * exp(stats::rnorm(n_features, mean = 0, sd = 0.25)), 0.01), 1.50)

  direction <- sample(c(-1, 1), n_de, replace = TRUE)
  abs_lfc <- if (identical(lfc_profile, "weak_mixture")) {
    stats::runif(n_de, min = simulation_weak_lfc_min, max = simulation_weak_lfc_max)
  } else {
    rep(lfc_magnitude, n_de)
  }
  true_lfc <- rep(0, n_features)
  true_lfc[de_id] <- direction * abs_lfc

  control_mu <- base_mu
  treatment_mu <- base_mu * 2^true_lfc
  control_lib <- exp(stats::rnorm(n_samples, mean = 0, sd = 0.12))
  treatment_lib <- exp(stats::rnorm(n_samples, mean = 0, sd = 0.12))
  control_lib <- control_lib / exp(mean(log(control_lib)))
  treatment_lib <- treatment_lib / exp(mean(log(treatment_lib)))

  generate_group <- function(mu, lib_factor, dispersion_vector) {
    out <- matrix(0L, nrow = length(mu), ncol = length(lib_factor))
    for (j in seq_along(lib_factor)) {
      out[, j] <- stats::rnbinom(n = length(mu), mu = pmax(mu * lib_factor[j], 1e-3), size = 1 / dispersion_vector)
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
    true_cusp = seq_len(n_features) %in% de_id & abs(true_lfc) > 0 & abs(true_lfc) < lfc_boundary,
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

simulation_metric_row <- function(template, method, truth_target, predicted, actual, hc_p, hbfss_cutoff, empirical_null_method) {
  cbind(
    template,
    data.frame(method = method, truth_target = truth_target, hc_p_threshold = hc_p, hbfss_cutoff = hbfss_cutoff, hc_valid = is.finite(hc_p) & !is.na(hc_p), empirical_null_method = empirical_null_method, stringsAsFactors = FALSE),
    simulation_metrics(predicted, actual)
  )
}

sequence_simulation_analysis <- function(sim_obj, null_inflation) {
  counts_mat <- as_integer_count_matrix(sim_obj$counts, "simulation DESeq2 input")
  n_samp <- as.integer(ncol(counts_mat) / 2L)
  coldata <- data.frame(condition = factor(c(rep("untrt", n_samp), rep("trt", n_samp)), levels = c("untrt", "trt")), row.names = colnames(counts_mat))

  dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  if (nrow(dds) < 5L) stop("Simulation DESeq2 object has fewer than five nonzero features after filtering.", call. = FALSE)
  dds <- DESeq(dds, betaPrior = FALSE, quiet = TRUE, fitType = simulation_fit_type)
  coef_name <- condition_coef_name(dds)

  standard <- results(dds, contrast = c("condition", "trt", "untrt"), alpha = alpha_standard)
  greater_abs <- results(dds, contrast = c("condition", "trt", "untrt"), lfcThreshold = lfc_boundary, altHypothesis = "greaterAbs", alpha = alpha_strong)
  less_abs <- results(dds, contrast = c("condition", "trt", "untrt"), lfcThreshold = lfc_boundary, altHypothesis = "lessAbs", alpha = alpha_weak)
  shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)

  df <- as.data.frame(standard)
  df$feature_id <- rownames(df)
  df$lfc_shrunk <- as.data.frame(shrunk)$log2FoldChange
  df$greaterAbs_padj <- as.data.frame(greater_abs)$padj
  df$lessAbs_padj <- as.data.frame(less_abs)$padj
  df <- left_join(df, sim_obj$truth, by = "feature_id")
  df$is_de[is.na(df$is_de)] <- FALSE
  df$true_cusp[is.na(df$true_cusp)] <- FALSE
  df$true_strong[is.na(df$true_strong)] <- FALSE

  inflated_wald <- inflate_null_wald_statistics(df$stat, df$is_de, null_inflation)
  empirical <- fit_empirical_null(inflated_wald, "simulation")

  df$empirical_p <- empirical$empirical_p
  df$empirical_bh <- empirical$empirical_bh
  df$neglog10_empirical_p_calc <- safe_neglog10(df$empirical_p, floor_value = calculation_p_floor)

  hc_p <- hc_threshold(df$empirical_p)
  hbfss_cutoff <- if (is.na(hc_p)) NA_real_ else -log10(hc_p) * lfc_boundary

  df$hc_pass <- !is.na(hc_p) & !is.na(df$empirical_p) & is.finite(df$empirical_p) & df$empirical_p <= hc_p
  df$HBFSS <- abs(df$lfc_shrunk) * df$neglog10_empirical_p_calc

  # DESeq2 standard Wald/BH call. No LFC filter is applied here by design;
  # LFC-boundary questions are tested by the composite-null strong/cusp arms.
  df$DESeq2_BH_sig <- !is.na(df$padj) & df$padj < alpha_standard
  df$DESeq2_Strong_sig <- !is.na(df$greaterAbs_padj) & df$greaterAbs_padj < alpha_strong & !is.na(df$lfc_shrunk) & abs(df$lfc_shrunk) >= lfc_boundary
  df$DESeq2_Weak_sig <- !is.na(df$lessAbs_padj) & df$lessAbs_padj < alpha_weak & !is.na(df$lfc_shrunk) & abs(df$lfc_shrunk) < lfc_boundary
  df$Empirical_BH_sig <- !is.na(df$empirical_bh) & df$empirical_bh < alpha_standard

  # HBFSS calls from empirical p-values and higher criticism.
  df$HBFSS_All_sig <- !is.na(hbfss_cutoff) & !is.na(df$HBFSS) & is.finite(df$HBFSS) & df$HBFSS >= hbfss_cutoff & df$hc_pass
  df$HBFSS_Strong_sig <- df$HBFSS_All_sig & !is.na(df$greaterAbs_padj) & df$greaterAbs_padj < alpha_strong & !is.na(df$lfc_shrunk) & abs(df$lfc_shrunk) >= lfc_boundary
  df$HBFSS_Weak_sig <- df$hc_pass & !is.na(df$lessAbs_padj) & df$lessAbs_padj < alpha_weak & !is.na(df$lfc_shrunk) & abs(df$lfc_shrunk) < lfc_boundary

  list(results = df, hc_p = hc_p, hbfss_cutoff = hbfss_cutoff, empirical_null_method = empirical$empirical_null_method)
}

build_metric_rows <- function(out, template) {
  df <- out$results
  bind_rows(
    simulation_metric_row(template, "DESeq2_BH", "all_de", df$DESeq2_BH_sig, df$is_de, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "Empirical_BH", "all_de", df$Empirical_BH_sig, df$is_de, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "HBFSS_All", "all_de", df$HBFSS_All_sig, df$is_de, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "DESeq2_Strong", "strong_de", df$DESeq2_Strong_sig, df$true_strong, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "HBFSS_Strong", "strong_de", df$HBFSS_Strong_sig, df$true_strong, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "DESeq2_Weak", "cusp_lfc", df$DESeq2_Weak_sig, df$true_cusp, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "HBFSS_Weak", "cusp_lfc", df$HBFSS_Weak_sig, df$true_cusp, out$hc_p, out$hbfss_cutoff, out$empirical_null_method)
  )
}

# =============================================================================
# SUMMARIES
# =============================================================================

summarise_metrics <- function(metric_long) {
  metric_long %>%
    group_by(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label, condition_label, method, truth_target, empirical_null_method) %>%
    summarise(
      n_replicates = n(),
      precision_mean = mean_finite(precision), precision_median = median_finite(precision), precision_q25 = quantile_finite(precision, 0.25), precision_q75 = quantile_finite(precision, 0.75),
      recall_mean = mean_finite(recall), recall_median = median_finite(recall), recall_q25 = quantile_finite(recall, 0.25), recall_q75 = quantile_finite(recall, 0.75),
      f1_mean = mean_finite(f1), f1_median = median_finite(f1), f1_q25 = quantile_finite(f1, 0.25), f1_q75 = quantile_finite(f1, 0.75),
      fdr_mean = mean_finite(fdr), fdr_median = median_finite(fdr), fdr_q25 = quantile_finite(fdr, 0.25), fdr_q75 = quantile_finite(fdr, 0.75),
      discovery_count_mean = mean_finite(discovery_count), discovery_count_median = median_finite(discovery_count), discovery_count_q25 = quantile_finite(discovery_count, 0.25), discovery_count_q75 = quantile_finite(discovery_count, 0.75),
      truth_count_mean = mean_finite(truth_count),
      hc_valid_fraction = mean(as.numeric(hc_valid), na.rm = TRUE),
      hc_p_threshold_median = median_finite(hc_p_threshold),
      hbfss_cutoff_median = median_finite(hbfss_cutoff),
      .groups = "drop"
    )
}

summarise_qc <- function(metric_long) {
  metric_long %>%
    group_by(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label, condition_label, replicate, empirical_null_method) %>%
    summarise(hc_p_threshold = first(hc_p_threshold), hbfss_cutoff = first(hbfss_cutoff), hc_valid = first(hc_valid), .groups = "drop") %>%
    group_by(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label, condition_label, empirical_null_method) %>%
    summarise(
      n_replicates = n(),
      hc_valid_fraction = mean(as.numeric(hc_valid), na.rm = TRUE),
      hc_p_threshold_min = min_finite(ifelse(hc_valid, hc_p_threshold, NA_real_)),
      hc_p_threshold_median = median_finite(ifelse(hc_valid, hc_p_threshold, NA_real_)),
      hc_p_threshold_max = max_finite(ifelse(hc_valid, hc_p_threshold, NA_real_)),
      hbfss_cutoff_median = median_finite(ifelse(hc_valid, hbfss_cutoff, NA_real_)),
      .groups = "drop"
    )
}

comparison_summary <- function(summary_df, target, baseline_method, test_method) {
  base <- summary_df %>%
    filter(truth_target == target, method == baseline_method) %>%
    select(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label, condition_label,
           baseline_precision = precision_mean, baseline_recall = recall_mean, baseline_f1 = f1_mean, baseline_fdr = fdr_mean, baseline_discovery_count = discovery_count_mean, baseline_truth_count = truth_count_mean)
  test <- summary_df %>%
    filter(truth_target == target, method == test_method) %>%
    select(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label, condition_label,
           test_precision = precision_mean, test_recall = recall_mean, test_f1 = f1_mean, test_fdr = fdr_mean, test_discovery_count = discovery_count_mean, test_truth_count = truth_count_mean)
  out <- inner_join(base, test, by = c("de_fraction", "lfc_magnitude", "simulation_profile", "null_inflation", "lfc_label", "de_label", "inflation_label", "condition_label"))
  if (nrow(out) == 0L) return(out)
  out$delta_precision <- out$test_precision - out$baseline_precision
  out$delta_recall <- out$test_recall - out$baseline_recall
  out$delta_f1 <- out$test_f1 - out$baseline_f1
  out$delta_fdr <- out$test_fdr - out$baseline_fdr
  out$delta_discovery_count <- out$test_discovery_count - out$baseline_discovery_count
  out$baseline_method <- baseline_method
  out$test_method <- test_method
  out$truth_target <- target
  out
}

# =============================================================================
# FIGURES
# =============================================================================

metric_labels <- c(f1 = "F1", recall = "Recall", precision = "Precision", fdr = "Observed FDR", discovery_count = "Discovery count")

stack_metric_summary <- function(summary_df, methods, target, metrics) {
  df <- summary_df[summary_df$method %in% methods & summary_df$truth_target == target & summary_df$truth_count_mean > 0, , drop = FALSE]
  if (nrow(df) == 0L) return(data.frame())
  bind_rows(lapply(metrics, function(m) {
    data.frame(
      lfc_label = df$lfc_label,
      de_label = df$de_label,
      inflation_label = df$inflation_label,
      condition_label = df$condition_label,
      method = df$method,
      method_label = unname(method_label_map[df$method]),
      metric = m,
      metric_label = unname(metric_labels[m]),
      value = df[[paste0(m, "_median")]],
      q25 = df[[paste0(m, "_q25")]],
      q75 = df[[paste0(m, "_q75")]],
      stringsAsFactors = FALSE
    )
  }))
}

plot_delta_heatmap <- function(delta_df, title, subtitle) {
  if (nrow(delta_df) == 0L) return(NULL)
  plot_df <- bind_rows(
    data.frame(condition_label = delta_df$condition_label, metric = "ΔF1", value = delta_df$delta_f1),
    data.frame(condition_label = delta_df$condition_label, metric = "ΔRecall", value = delta_df$delta_recall),
    data.frame(condition_label = delta_df$condition_label, metric = "ΔPrecision", value = delta_df$delta_precision),
    data.frame(condition_label = delta_df$condition_label, metric = "ΔFDR", value = delta_df$delta_fdr),
    data.frame(condition_label = delta_df$condition_label, metric = "ΔDiscoveries", value = delta_df$delta_discovery_count)
  )
  plot_df$condition_label <- factor(plot_df$condition_label, levels = rev(unique(delta_df$condition_label)))
  plot_df$metric <- factor(plot_df$metric, levels = c("ΔF1", "ΔRecall", "ΔPrecision", "ΔFDR", "ΔDiscoveries"))
  lim <- max(abs(plot_df$value), na.rm = TRUE)
  if (!is.finite(lim) || is.na(lim) || lim <= 0) lim <- 1
  ggplot(plot_df, aes(x = metric, y = condition_label, fill = value)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = formatC(value, format = "f", digits = 2)), size = 2.0) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-lim, lim), name = "Delta") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL, caption = figure_stamp()) +
    manuscript_theme() +
    theme(axis.text.y = element_text(size = 5.2))
}

plot_head_to_head <- function(summary_df, methods, target, metrics, title, subtitle) {
  plot_df <- stack_metric_summary(summary_df, methods, target, metrics)
  if (nrow(plot_df) == 0L) return(NULL)
  plot_df$method_label <- factor(plot_df$method_label, levels = unname(method_label_map[methods]))
  plot_df$metric_label <- factor(plot_df$metric_label, levels = unname(metric_labels[metrics]))
  dodge <- position_dodge(width = 0.45)
  p <- ggplot(plot_df, aes(x = de_label, y = value, fill = method_label, color = method_label, group = method_label)) +
    geom_crossbar(aes(ymin = q25, ymax = q75), width = 0.35, fatten = 1.1, position = dodge, linewidth = 0.25, alpha = 0.85) +
    facet_grid(metric_label ~ inflation_label + lfc_label, scales = "free_y") +
    scale_fill_manual(values = unname(method_color_map[methods]), breaks = unname(method_label_map[methods])) +
    scale_color_manual(values = unname(method_color_map[methods]), breaks = unname(method_label_map[methods])) +
    labs(title = title, subtitle = subtitle, x = "DE fraction", y = "Median with IQR", caption = figure_stamp()) +
    manuscript_theme()
  if (all(metrics %in% c("f1", "recall", "precision", "fdr"))) {
    p <- p + scale_y_continuous(labels = scales::label_percent(accuracy = 1))
  }
  p
}

plot_hc_qc <- function(qc_df, value_col, title, subtitle, fill_label) {
  if (nrow(qc_df) == 0L) return(NULL)
  plot_df <- qc_df
  plot_df$value <- plot_df[[value_col]]
  ggplot(plot_df, aes(x = de_label, y = inflation_label, fill = value)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = formatC(value, format = "f", digits = 2)), size = 2.0) +
    facet_wrap(~ lfc_label, nrow = 1) +
    scale_fill_gradient(low = "white", high = "#0072B2", name = fill_label) +
    labs(title = title, subtitle = subtitle, x = "DE fraction", y = NULL, caption = figure_stamp()) +
    manuscript_theme()
}

# =============================================================================
# MAIN RUN
# =============================================================================

run_simulation_validation <- function() {
  set.seed(simulation_seed)
  fixed_grid <- expand.grid(de_fraction = simulation_de_fractions, lfc_magnitude = simulation_lfc_magnitudes, null_inflation = simulation_null_inflation, simulation_profile = "fixed", stringsAsFactors = FALSE)
  weak_grid <- expand.grid(de_fraction = simulation_de_fractions, lfc_magnitude = NA_real_, null_inflation = simulation_null_inflation, simulation_profile = "weak_mixture", stringsAsFactors = FALSE)
  sim_grid <- bind_rows(fixed_grid, weak_grid)
  total_runs <- nrow(sim_grid) * simulation_n_reps

  log_message("Simulation started: ", total_runs, " planned runs")
  write_status("running", paste0("0/", total_runs))

  metric_rows <- list()
  failure_rows <- list()
  run_index <- 0L
  start_time <- Sys.time()
  checkpoint_path <- file.path(simulation_dir, "Simulation_RunMetrics_Checkpoint.csv")

  for (grid_i in seq_len(nrow(sim_grid))) {
    grid_row <- sim_grid[grid_i, , drop = FALSE]
    for (rep_i in seq_len(simulation_n_reps)) {
      run_index <- run_index + 1L
      if (run_index %% simulation_progress_every == 0L || run_index == 1L || run_index == total_runs) {
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
        eta <- elapsed * (total_runs - run_index) / max(run_index, 1L)
        log_message("Simulation progress: ", run_index, "/", total_runs, " | elapsed=", signif(elapsed, 3), " min | ETA=", signif(eta, 3), " min")
        write_status("running", paste0(run_index, "/", total_runs))
      }

      template <- data.frame(de_fraction = grid_row$de_fraction, lfc_magnitude = grid_row$lfc_magnitude, simulation_profile = grid_row$simulation_profile, null_inflation = grid_row$null_inflation, replicate = rep_i, stringsAsFactors = FALSE)
      template$lfc_label <- lfc_label(template$lfc_magnitude, template$simulation_profile)
      template$de_label <- paste0(template$de_fraction * 100, "% DE")
      template$inflation_label <- inflation_label(template$null_inflation)
      template$condition_label <- condition_label(template$lfc_label, template$inflation_label, template$de_fraction)

      out <- tryCatch({
        sim_obj <- simulate_sequence_counts(simulation_n_features, simulation_n_samples_per_group, simulation_base_mean, simulation_dispersion_null, grid_row$de_fraction, grid_row$lfc_magnitude, simulation_dispersion_de, grid_row$simulation_profile)
        sequence_simulation_analysis(sim_obj, null_inflation = grid_row$null_inflation)
      }, error = function(e) {
        failure_rows[[length(failure_rows) + 1L]] <<- cbind(template, data.frame(error_message = conditionMessage(e), stringsAsFactors = FALSE))
        NULL
      })
      if (is.null(out)) next
      metric_rows[[length(metric_rows) + 1L]] <- build_metric_rows(out, template)
      if ((run_index %% simulation_checkpoint_every == 0L || run_index == total_runs) && length(metric_rows) > 0L) save_csv(bind_rows(metric_rows), checkpoint_path)
    }
  }

  failures <- if (length(failure_rows) > 0L) bind_rows(failure_rows) else data.frame()
  save_csv(failures, file.path(simulation_dir, "Simulation_Failures.csv"))
  metrics <- if (length(metric_rows) > 0L) bind_rows(metric_rows) else data.frame()
  if (nrow(metrics) == 0L) stop("Simulation produced no successful replicates.", call. = FALSE)

  metric_summary <- summarise_metrics(metrics)
  qc_summary <- summarise_qc(metrics)
  save_csv(metrics, file.path(simulation_dir, "Simulation_RunMetrics_Long.csv"))
  save_csv(metric_summary, file.path(simulation_dir, "Simulation_Condition_Method_Summary.csv"))
  save_csv(qc_summary, file.path(simulation_dir, "Simulation_HC_QC_Summary.csv"))
  if (file.exists(checkpoint_path)) unlink(checkpoint_path, force = TRUE)

  delta_all <- comparison_summary(metric_summary, "all_de", "DESeq2_BH", "HBFSS_All")
  delta_strong <- comparison_summary(metric_summary, "strong_de", "DESeq2_Strong", "HBFSS_Strong")
  delta_weak <- comparison_summary(metric_summary, "cusp_lfc", "DESeq2_Weak", "HBFSS_Weak")
  save_csv(delta_all, file.path(simulation_dir, "Simulation_Delta_All_HBFSSvsDESeq2.csv"))
  save_csv(delta_strong, file.path(simulation_dir, "Simulation_Delta_Strong_HBFSSvsDESeq2Composite.csv"))
  save_csv(delta_weak, file.path(simulation_dir, "Simulation_Delta_Weak_HBFSSvsDESeq2Composite.csv"))

  figures <- bind_rows(
    safe_save_figure(plot_delta_heatmap(delta_all, "All-DE performance: HBFSS versus DESeq2 BH", "Positive values favor HBFSS; positive ΔFDR indicates additional false-discovery cost"), "Figure_01_AllDE_DeltaHeatmap", 6.7, 5.2),
    safe_save_figure(plot_head_to_head(metric_summary, c("DESeq2_Strong", "HBFSS_Strong"), "strong_de", c("f1", "recall", "precision", "fdr"), "Strong-effect comparison", "DESeq2 greaterAbs composite-null calls versus HBFSS strong calls"), "Figure_02_StrongEffect_HeadToHead", 6.7, 5.4),
    safe_save_figure(plot_head_to_head(metric_summary, c("DESeq2_Weak", "HBFSS_Weak"), "cusp_lfc", c("f1", "recall", "precision", "fdr"), "Composite-null cusp comparison", "DESeq2 lessAbs composite-null calls versus HC-supported HBFSS cusp calls"), "Figure_03_CompositeCusp_HeadToHead", 6.7, 5.4),
    safe_save_figure(plot_hc_qc(qc_summary, "hc_valid_fraction", "Higher-criticism threshold validity", "Fraction of replicates with usable HC threshold", "Valid fraction"), "Figure_04_HC_Validity", 6.7, 3.0),
    safe_save_figure(plot_hc_qc(qc_summary, "hbfss_cutoff_median", "Median HBFSS cutoff", "HC-derived HBFSS cutoff by simulation condition", "Median cutoff"), "Figure_05_HBFSS_Cutoff", 6.7, 3.0)
  )
  save_csv(figures, file.path(simulation_dir, "Figure_Export_Status.csv"))
  failed_figures <- figures[figures$status != "saved", , drop = FALSE]
  if (nrow(failed_figures) > 0L) {
    save_csv(failed_figures, file.path(simulation_dir, "Figure_Export_Failures.csv"))
    stop("One or more required figures failed to export. See Figure_Export_Failures.csv.", call. = FALSE)
  }

  list(metrics = metrics, metric_summary = metric_summary, qc_summary = qc_summary, delta_all = delta_all, delta_strong = delta_strong, delta_weak = delta_weak, failures = failures, figures = figures)
}

# =============================================================================
# TEXT OUTPUTS / GIT
# =============================================================================

fmt_num <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(x) || is.na(x)) return("NA")
  formatC(x, format = "f", digits = digits)
}

fmt_delta <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(x) || is.na(x)) return("NA")
  paste0(ifelse(x >= 0, "+", ""), formatC(x, format = "f", digits = digits))
}

write_methods_and_interpretation <- function(result) {
  methods_lines <- c(
    "# Simulation methods",
    "",
    paste0("Negative-binomial count matrices were simulated with ", simulation_n_features, " features and ", simulation_n_samples_per_group, " samples per condition."),
    paste0("DE fractions: ", paste(simulation_de_fractions, collapse = ", "), ". Fixed |LFC| values: ", paste(simulation_lfc_magnitudes, collapse = ", "), ". Weak mixture |LFC| range: ", simulation_weak_lfc_min, "-", simulation_weak_lfc_max, "."),
    paste0("Each condition used ", simulation_n_reps, " replicate simulations. DESeq2 fitType was '", simulation_fit_type, "'."),
    "The simulation truth target for the lessAbs/cusp arm is true_cusp = 0 < |true LFC| < LFC boundary. This excludes exact null features and treats the weak arm as a boundary-adjacent composite-null classification problem, not conventional differential-expression recovery.",
    "",
    "DESeq2_BH uses DESeq2 Wald-test p-values with DESeq2's BH-adjusted padj only; it intentionally does not impose an LFC filter because the strong and cusp arms are handled by DESeq2's composite-null tests.",
    "DESeq2_Strong uses the DESeq2 greaterAbs composite-null test plus apeglm-shrunken |LFC| >= boundary.",
    "DESeq2_Weak uses the DESeq2 lessAbs composite-null test and apeglm-shrunken |LFC| < boundary. It is not a standard Wald/BH differential-expression call.",
    "Empirical p-values were estimated from the DESeq2 Wald statistic by fdrtool. Higher criticism selected one empirical-p threshold per simulated dataset.",
    "HBFSS = abs(apeglm-shrunken log2FC) * -log10(empirical p). HBFSS cutoff = -log10(HC p-threshold) * LFC boundary. No HBFSS floor was added.",
    "HBFSS_All uses the HC-derived HBFSS boundary. HBFSS_Strong is the HBFSS all call restricted to greaterAbs-supported strong effects. HBFSS_Weak is the HC-supported, lessAbs-supported cusp/sub-boundary call."
  )
  writeLines(methods_lines, file.path(simulation_dir, "Methods_Simulation.txt"))

  interp <- c(
    "# Simulation interpretation",
    "",
    "The figures are organized around direct comparisons. Figure 01 compares HBFSS_All against DESeq2_BH. Figure 02 compares HBFSS_Strong against the DESeq2 greaterAbs composite-null strong-effect method. Figure 03 compares HBFSS_Weak against the DESeq2 lessAbs composite-null cusp method. Figures 04 and 05 document HC threshold behavior.",
    "",
    paste0("Mean all-DE ΔF1: ", fmt_delta(mean_finite(result$delta_all$delta_f1))),
    paste0("Mean all-DE ΔRecall: ", fmt_delta(mean_finite(result$delta_all$delta_recall))),
    paste0("Mean all-DE ΔFDR: ", fmt_delta(mean_finite(result$delta_all$delta_fdr))),
    paste0("Mean strong-effect ΔF1: ", fmt_delta(mean_finite(result$delta_strong$delta_f1))),
    paste0("Mean composite-cusp ΔF1: ", fmt_delta(mean_finite(result$delta_weak$delta_f1))),
    paste0("Mean HC valid-threshold fraction: ", fmt_num(mean_finite(result$qc_summary$hc_valid_fraction))),
    "",
    "Positive ΔF1, ΔRecall, or ΔPrecision favors HBFSS. Positive ΔFDR indicates that the gain carries a false-discovery cost. The results should be described as sensitivity and discovery tradeoffs rather than as automatic universal FDR improvement."
  )
  writeLines(interp, file.path(simulation_dir, "Interpretation_Simulation.md"))
  writeLines(interp, file.path(simulation_dir, "Interpretation_Simulation.txt"))
}

write_manifest <- function() {
  files <- list.files(simulation_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files)]
  manifest <- data.frame(file = vapply(files, repo_relative, character(1)), size_bytes = file.info(files)$size, stringsAsFactors = FALSE)
  save_csv(manifest, file.path(simulation_dir, "Manifest_Simulation.csv"))
  manifest
}

git_run <- function(args, allow_failure = FALSE) {
  if (Sys.which("git") == "") stop("Git executable was not found on PATH.", call. = FALSE)
  out <- suppressWarnings(system2("git", args = c("-C", repo_root, args), stdout = TRUE, stderr = TRUE, timeout = 180))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (length(out) > 0L) log_message(paste(out, collapse = "\n"))
  if (status != 0L && !allow_failure) stop("Git command failed: git ", paste(args, collapse = " "), "\n", paste(out, collapse = "\n"), call. = FALSE)
  list(status = status, output = out)
}

push_outputs <- function() {
  if (!isTRUE(git_push_after_success)) {
    log_message("Git push skipped. Set SEQUENCE_GIT_PUSH=true to enable.")
    return(invisible(FALSE))
  }
  if (!dir.exists(file.path(repo_root, ".git"))) stop("No .git directory found for push.", call. = FALSE)
  current_branch <- trimws(system2("git", args = c("-C", repo_root, "rev-parse", "--abbrev-ref", "HEAD"), stdout = TRUE))
  target_branch <- if (!is.na(git_branch_name) && nzchar(git_branch_name)) git_branch_name else current_branch
  paths_to_stage <- c(repo_relative(simulation_dir))
  sp <- script_path()
  if (!is.na(sp) && file.exists(sp)) {
    paths_to_stage <- c(paths_to_stage, repo_relative(sp))
  }
  git_run(c("add", "--", unique(paths_to_stage)))
  changed <- git_run(c("diff", "--cached", "--quiet"), allow_failure = TRUE)
  if (changed$status == 0L) {
    log_message("No Git changes to commit.")
    return(invisible(FALSE))
  }
  git_run(c("commit", "-m", git_commit_message))
  git_run(c("push", "-u", git_remote_name, paste0("HEAD:", target_branch)))
  log_message("Git push complete.")
  invisible(TRUE)
}

main <- function() {
  write_status("starting", "initializing")
  log_message("Method-aligned SEQUENCE simulation started.")
  log_message("Repository root: ", repo_root)
  log_message("Output directory: ", simulation_dir)
  result <- run_simulation_validation()
  write_methods_and_interpretation(result)
  writeLines(capture.output(sessionInfo()), file.path(simulation_dir, "SessionInfo_Simulation.txt"))
  manifest <- write_manifest()
  cleanup_rplots_pdf()
  write_status("outputs_complete", paste0("files=", nrow(manifest)))
  push_outputs()
  write_status("complete", "simulation finished")
  cat("\nSimulation complete\n")
  cat("Output directory: ", simulation_dir, "\n", sep = "")
  cat("Exported files: ", nrow(manifest), "\n", sep = "")
  invisible(result)
}

tryCatch(
  main(),
  error = function(e) {
    write_status("failed", conditionMessage(e))
    log_message("SIMULATION FAILED - DO NOT PUSH: ", conditionMessage(e))
    try(writeLines(capture.output(sessionInfo()), file.path(simulation_dir, "SessionInfo_Simulation_Failed.txt")), silent = TRUE)
    try(cleanup_rplots_pdf(), silent = TRUE)
    if (!interactive()) quit(save = "no", status = 1L, runLast = FALSE)
    stop(e)
  }
)
