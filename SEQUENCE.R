#!/usr/bin/env Rscript

# =============================================================================
# SEQUENCE SIMULATION VALIDATION - INTERPRETABLE OVERHAUL
# Standalone manuscript simulation with direct, interpretable figures for
# DESeq2, empirical-null p-values, higher criticism, and HBFSS.
# =============================================================================

options(stringsAsFactors = FALSE)
options(width = 140)

# =============================================================================
# USER SETTINGS
# =============================================================================

output_folder_name <- "manuscript_final_clean"
reset_simulation_dir <- TRUE

alpha_standard <- 0.10
alpha_strong <- 0.10
alpha_weak <- 0.10
lfc_boundary <- 1.0

strict_empirical_null <- FALSE
hc_invalid_at_or_above <- 0.95
hc_invalid_below <- 0.001
calculation_p_floor <- .Machine$double.xmin
plot_p_floor <- 1e-16

figure_dpi <- 600
base_theme_size <- 7
figure_width_twocol <- 6.7
export_pdf_also <- TRUE
export_tiff_also <- TRUE
stamp_figures <- FALSE

simulation_seed <- as.integer(Sys.getenv("SEQUENCE_SIM_SEED", "42"))
simulation_n_features <- as.integer(Sys.getenv("SEQUENCE_SIM_FEATURES", "1000"))
simulation_n_samples_per_group <- as.integer(Sys.getenv("SEQUENCE_SIM_SAMPLES_PER_GROUP", "6"))
simulation_n_reps <- as.integer(Sys.getenv("SEQUENCE_SIM_REPS", "20"))
simulation_progress_every <- as.integer(Sys.getenv("SEQUENCE_SIM_PROGRESS_EVERY", "5"))
simulation_checkpoint_every <- as.integer(Sys.getenv("SEQUENCE_SIM_CHECKPOINT_EVERY", "5"))

simulation_base_mean <- 200
simulation_dispersion_null <- 0.10
simulation_dispersion_de <- 0.15
simulation_de_fractions <- c(0.05, 0.10)
simulation_lfc_magnitudes <- c(0.50, 1.00, 2.00)
simulation_weak_lfc_min <- 0.20
simulation_weak_lfc_max <- 0.80
simulation_null_inflation <- c(0.00, 0.10)
simulation_use_apeglm_shrinkage <- TRUE
simulation_fit_type <- "local"
simulation_fail_on_failed_replicates <- TRUE
export_simulation_feature_results <- FALSE

# Optional GitHub push after successful simulation run.
git_push_after_success <- TRUE
git_remote_name <- "origin"
git_branch_name <- NA_character_
git_commit_message <- paste0("Update SEQUENCE simulation validation outputs ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
git_pull_rebase_before_push <- TRUE
git_stage_pipeline_script <- TRUE
git_allow_no_change_success <- TRUE
git_max_file_size_bytes <- 95 * 1024^2
git_exclude_patterns <- c(
  "exports/manuscript_final_clean/simulation_validation/Simulation_FeatureResults.csv"
)

# =============================================================================
# PACKAGES
# =============================================================================

required_packages <- c("DESeq2", "fdrtool", "ggplot2", "dplyr", "grid")
if (isTRUE(simulation_use_apeglm_shrinkage)) required_packages <- c(required_packages, "apeglm")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install missing package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(fdrtool)
  library(ggplot2)
  library(dplyr)
  library(grid)
  if (isTRUE(simulation_use_apeglm_shrinkage)) library(apeglm)
})

# =============================================================================
# PATHS AND LOGGING
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
simulation_dir <- file.path(output_dir, "simulation_validation")
log_file <- file.path(simulation_dir, "Simulation_Log.txt")

if (isTRUE(reset_simulation_dir) && dir.exists(simulation_dir)) {
  unlink(simulation_dir, recursive = TRUE, force = TRUE)
}
dir.create(simulation_dir, recursive = TRUE, showWarnings = FALSE)
figures_dir <- file.path(simulation_dir, "figures")
figures_png_dir <- file.path(figures_dir, "png")
figures_pdf_dir <- file.path(figures_dir, "pdf")
figures_tiff_dir <- file.path(figures_dir, "tiff")
dir.create(figures_png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_tiff_dir, recursive = TRUE, showWarnings = FALSE)

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
simulation_n_reps <- positive_integer(simulation_n_reps, 25L, "simulation_n_reps")
simulation_progress_every <- positive_integer(simulation_progress_every, 5L, "simulation_progress_every")
simulation_checkpoint_every <- positive_integer(simulation_checkpoint_every, 5L, "simulation_checkpoint_every")

# =============================================================================
# GENERAL HELPERS
# =============================================================================

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

finite_vals <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[is.finite(y) & !is.na(y)]
}

mean_finite <- function(x) {
  v <- finite_vals(x)
  if (!length(v)) return(NA_real_)
  mean(v)
}

median_finite <- function(x) {
  v <- finite_vals(x)
  if (!length(v)) return(NA_real_)
  stats::median(v)
}

quantile_finite <- function(x, prob) {
  v <- finite_vals(x)
  if (!length(v)) return(NA_real_)
  as.numeric(stats::quantile(v, probs = prob, names = FALSE, type = 7, na.rm = TRUE))
}

min_finite <- function(x) {
  v <- finite_vals(x)
  if (!length(v)) return(NA_real_)
  min(v)
}

max_finite <- function(x) {
  v <- finite_vals(x)
  if (!length(v)) return(NA_real_)
  max(v)
}

percent_label <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  ifelse(is.na(y) | !is.finite(y), NA_character_, paste0(formatC(100 * y, format = "f", digits = 0), "%"))
}

percent_delta_label <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  ifelse(is.na(y) | !is.finite(y), "NA", paste0(ifelse(y >= 0, "+", ""), formatC(100 * y, format = "f", digits = 0), "%"))
}

save_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE)
  invisible(path)
}

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_theme_size + 0.5, hjust = 0.5),
      axis.title = element_text(face = "bold", size = base_theme_size - 0.1),
      axis.text = element_text(color = "black", size = base_theme_size - 0.8),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_theme_size - 0.2),
      legend.text = element_text(size = base_theme_size - 0.4),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey88"),
      strip.text = element_text(face = "bold", size = base_theme_size - 0.8),
      plot.margin = margin(6, 8, 6, 6)
    )
}

draw_to_device <- function(plot_obj) {
  grid::grid.newpage()
  if (inherits(plot_obj, "ggplot")) {
    print(plot_obj)
  } else {
    stop("Unsupported figure object class: ", paste(class(plot_obj), collapse = ", "), call. = FALSE)
  }
}

check_nonempty_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " was not written: ", path, call. = FALSE)
  file_size <- suppressWarnings(file.info(path)$size)
  if (!is.finite(file_size) || is.na(file_size) || file_size <= 0) {
    stop(label, " is empty: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

git_sha_short <- function() {
  if (Sys.which("git") == "" || !dir.exists(file.path(repo_root, ".git"))) return("NA")
  out <- suppressWarnings(system2("git", args = c("-C", repo_root, "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L || length(out) == 0L) return("NA")
  trimws(out[1])
}

apply_figure_stamp <- function(plot_obj) {
  if (!isTRUE(stamp_figures)) return(plot_obj)
  plot_obj + labs(caption = paste0("git: ", git_sha_short(), " | seed: ", simulation_seed, " | reps: ", simulation_n_reps))
}

save_figure <- function(plot_obj, figure_name, width, height, export_pdf = export_pdf_also) {
  if (is.null(plot_obj)) stop("Figure object is NULL for ", figure_name, call. = FALSE)
  if (!inherits(plot_obj, "ggplot")) stop("Figure object is not a ggplot for ", figure_name, call. = FALSE)

  plot_obj <- apply_figure_stamp(plot_obj)

  png_path <- file.path(figures_png_dir, paste0(figure_name, ".png"))
  pdf_path <- file.path(figures_pdf_dir, paste0(figure_name, ".pdf"))
  tiff_path <- file.path(figures_tiff_dir, paste0(figure_name, ".tiff"))

  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
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

  if (isTRUE(export_pdf)) {
    grDevices::pdf(file = pdf_path, width = width, height = height, onefile = TRUE, useDingbats = FALSE)
    tryCatch({
      draw_to_device(plot_obj)
    }, finally = {
      grDevices::dev.off()
    })
    check_nonempty_file(pdf_path, "PDF figure")
  } else {
    pdf_path <- NA_character_
  }

  if (isTRUE(export_tiff_also)) {
    ggplot2::ggsave(
      filename = tiff_path,
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      dpi = figure_dpi,
      compression = "lzw",
      bg = "white",
      limitsize = FALSE
    )
    check_nonempty_file(tiff_path, "TIFF figure")
  } else {
    tiff_path <- NA_character_
  }

  root_png_path <- file.path(simulation_dir, paste0(figure_name, ".png"))
  file.copy(png_path, root_png_path, overwrite = TRUE)
  check_nonempty_file(root_png_path, "Root PNG figure")

  log_message("Figure saved: ", root_png_path, " | ", png_path, if (isTRUE(export_pdf)) paste0(" | ", pdf_path) else "", if (isTRUE(export_tiff_also)) paste0(" | ", tiff_path) else "")
  data.frame(
    figure = figure_name,
    png_file = path_relative_to_repo(root_png_path),
    png_archive_file = path_relative_to_repo(png_path),
    pdf_file = if (isTRUE(export_pdf)) path_relative_to_repo(pdf_path) else NA_character_,
    tiff_file = if (isTRUE(export_tiff_also)) path_relative_to_repo(tiff_path) else NA_character_,
    status = "saved",
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

safe_save_figure <- function(plot_obj, figure_name, width, height) {
  tryCatch({
    save_figure(plot_obj, figure_name, width, height)
  }, error = function(e) {
    log_message("FIGURE FAILED: ", figure_name, " | ", conditionMessage(e))
    data.frame(
      figure = figure_name,
      png_file = NA_character_,
      png_archive_file = NA_character_,
      pdf_file = NA_character_,
      tiff_file = NA_character_,
      status = "failed",
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}

method_label_map <- c(
  DESeq2_BH = "DESeq2 BH",
  Empirical_BH = "Empirical BH",
  GreaterAbs = "GreaterAbs",
  LessAbs = "LessAbs",
  HBFSS_total = "HBFSS total",
  HBFSS_raw = "HBFSS raw",
  HBFSS_weak_region = "HBFSS weak"
)

method_color_map <- c(
  DESeq2_BH = "#000000",
  Empirical_BH = "#E69F00",
  GreaterAbs = "#D55E00",
  LessAbs = "#CC79A7",
  HBFSS_total = "#56B4E9",
  HBFSS_raw = "#009E73",
  HBFSS_weak_region = "#0072B2"
)

condition_short_label <- function(inflation_label, de_label) {
  paste0(de_label, " | ", inflation_label)
}

condition_long_label <- function(lfc_label, inflation_label, de_label) {
  paste0(gsub("\n", " | ", lfc_label), " | ", inflation_label, " | ", de_label)
}

# =============================================================================
# EMPIRICAL NULL AND HC/HBFSS
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

# =============================================================================
# SIMULATION CORE
# =============================================================================

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

simulate_sequence_counts <- function(n_features, n_samples, base_mean, disp_null, de_fraction, lfc_magnitude, disp_de, lfc_profile = c("fixed", "weak_mixture")) {
  lfc_profile <- match.arg(lfc_profile)
  n_de <- round(n_features * de_fraction)
  n_de <- max(2L, min(n_de, n_features - 2L))
  feature_id <- paste0("sim_feature_", seq_len(n_features))
  de_id <- sample(seq_len(n_features), n_de, replace = FALSE)

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

  counts_mat <- cbind(
    generate_group(control_mu, control_lib, dispersion),
    generate_group(treatment_mu, treatment_lib, dispersion)
  )
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

simulation_metric_row <- function(template, method, truth_target, predicted, actual, hc_p, hbfss_cutoff, empirical_null_method) {
  cbind(
    template,
    data.frame(
      method = method,
      truth_target = truth_target,
      hc_p_threshold = hc_p,
      hbfss_cutoff = hbfss_cutoff,
      hc_valid = is.finite(hc_p) & !is.na(hc_p),
      empirical_null_method = empirical_null_method,
      stringsAsFactors = FALSE
    ),
    simulation_metrics(predicted, actual)
  )
}

sequence_simulation_analysis <- function(sim_obj, null_inflation) {
  counts_mat <- sim_obj$counts

  if (!is.matrix(counts_mat) && !is.data.frame(counts_mat)) {
    stop("Simulation counts must be a matrix-like object.", call. = FALSE)
  }
  if (ncol(counts_mat) %% 2L != 0L) {
    stop("Simulation count matrix must contain paired control/treatment sample columns.", call. = FALSE)
  }

  n_samp <- as.integer(ncol(counts_mat) / 2L)
  if (n_samp < 2L) {
    stop("Simulation requires at least two samples per group.", call. = FALSE)
  }

  counts_mat <- as_integer_count_matrix(counts_mat, "simulation DESeq2 input")

  coldata <- data.frame(
    condition = factor(c(rep("untrt", n_samp), rep("trt", n_samp)), levels = c("untrt", "trt")),
    row.names = colnames(counts_mat)
  )

  dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  if (nrow(dds) < 5L) stop("Simulation DESeq2 object has fewer than five nonzero features after filtering.", call. = FALSE)

  dds <- DESeq(dds, betaPrior = FALSE, quiet = TRUE, fitType = simulation_fit_type)
  coef_name <- condition_coef_name(dds)

  standard <- results(dds, contrast = c("condition", "trt", "untrt"), alpha = alpha_standard)
  greater_abs <- results(dds, contrast = c("condition", "trt", "untrt"), lfcThreshold = lfc_boundary, altHypothesis = "greaterAbs", alpha = alpha_strong)
  less_abs <- results(dds, contrast = c("condition", "trt", "untrt"), lfcThreshold = lfc_boundary, altHypothesis = "lessAbs", alpha = alpha_weak)

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
      data.frame(feature_id = rownames(greater_abs), greaterAbs_padj = greater_abs$padj, stringsAsFactors = FALSE),
      by = "feature_id"
    ) %>%
    left_join(
      data.frame(feature_id = rownames(less_abs), lessAbs_padj = less_abs$padj, stringsAsFactors = FALSE),
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

  empirical <- fit_empirical_null(wald_for_empirical_null, "simulation")

  df$empirical_p <- empirical$empirical_p
  df$empirical_bh <- empirical$empirical_bh
  df$empirical_null_method <- empirical$empirical_null_method
  df$neglog10_empirical_p_calc <- safe_neglog10(df$empirical_p, floor_value = calculation_p_floor)

  hc_p <- hc_threshold(df$empirical_p)
  hbfss_cutoff <- if (is.na(hc_p)) NA_real_ else -log10(hc_p) * lfc_boundary

  df$hc_pass <- !is.na(hc_p) & !is.na(df$empirical_p) & is.finite(df$empirical_p) & df$empirical_p <= hc_p

  # Manuscript definition. Do not add a cutoff floor or change the boundary.
  df$HBFSS <- abs(df$lfc_shrunk) * df$neglog10_empirical_p_calc

  df$standard_sig <- !is.na(df$padj) & df$padj < alpha_standard & !is.na(df$lfc_shrunk) & abs(df$lfc_shrunk) >= lfc_boundary
  df$greaterAbs_sig <- !is.na(df$greaterAbs_padj) & df$greaterAbs_padj < alpha_strong & !is.na(df$lfc_shrunk) & abs(df$lfc_shrunk) >= lfc_boundary
  df$lessAbs_sig <- !is.na(df$lessAbs_padj) & df$lessAbs_padj < alpha_weak & !is.na(df$lfc_shrunk) & abs(df$lfc_shrunk) < lfc_boundary
  df$empirical_bh_sig <- !is.na(df$empirical_bh) & df$empirical_bh < alpha_standard & !is.na(df$lfc_shrunk) & abs(df$lfc_shrunk) >= lfc_boundary

  df$hbfss_raw_sig <- !is.na(hbfss_cutoff) & !is.na(df$HBFSS) & is.finite(df$HBFSS) & df$HBFSS >= hbfss_cutoff & df$hc_pass
  df$hbfss_total_sig <- df$standard_sig | df$hbfss_raw_sig
  df$weak_region_hbfss_sig <- df$lessAbs_sig & df$hc_pass

  list(
    results = df,
    hc_p = hc_p,
    hbfss_cutoff = hbfss_cutoff,
    empirical_null_method = empirical$empirical_null_method
  )
}

build_simulation_metric_rows <- function(out, template) {
  df <- out$results
  bind_rows(
    simulation_metric_row(template, "DESeq2_BH", "all_de", df$standard_sig, df$is_de, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "Empirical_BH", "all_de", df$empirical_bh_sig, df$is_de, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "GreaterAbs", "strong_de", df$greaterAbs_sig, df$true_strong, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "LessAbs", "weak_de", df$lessAbs_sig, df$true_weak, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "HBFSS_total", "all_de", df$hbfss_total_sig, df$is_de, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "HBFSS_raw", "all_de", df$hbfss_raw_sig, df$is_de, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "HBFSS_raw", "strong_de", df$hbfss_raw_sig, df$true_strong, out$hc_p, out$hbfss_cutoff, out$empirical_null_method),
    simulation_metric_row(template, "HBFSS_weak_region", "weak_de", df$weak_region_hbfss_sig, df$true_weak, out$hc_p, out$hbfss_cutoff, out$empirical_null_method)
  )
}

summarise_metric_long <- function(metric_long) {
  metric_long %>%
    group_by(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label, condition_label, condition_short, method, truth_target, empirical_null_method) %>%
    summarise(
      n_replicates = n(),
      precision_mean = mean_finite(precision),
      precision_median = median_finite(precision),
      precision_q25 = quantile_finite(precision, 0.25),
      precision_q75 = quantile_finite(precision, 0.75),
      recall_mean = mean_finite(recall),
      recall_median = median_finite(recall),
      recall_q25 = quantile_finite(recall, 0.25),
      recall_q75 = quantile_finite(recall, 0.75),
      f1_mean = mean_finite(f1),
      f1_median = median_finite(f1),
      f1_q25 = quantile_finite(f1, 0.25),
      f1_q75 = quantile_finite(f1, 0.75),
      fdr_mean = mean_finite(fdr),
      fdr_median = median_finite(fdr),
      fdr_q25 = quantile_finite(fdr, 0.25),
      fdr_q75 = quantile_finite(fdr, 0.75),
      discovery_count_mean = mean_finite(discovery_count),
      discovery_count_median = median_finite(discovery_count),
      discovery_count_q25 = quantile_finite(discovery_count, 0.25),
      discovery_count_q75 = quantile_finite(discovery_count, 0.75),
      truth_count_mean = mean_finite(truth_count),
      hc_valid_fraction = mean(as.numeric(hc_valid), na.rm = TRUE),
      hc_p_threshold_median = median_finite(hc_p_threshold),
      hbfss_cutoff_median = median_finite(hbfss_cutoff),
      .groups = "drop"
    )
}

summarise_qc <- function(metric_long) {
  metric_long %>%
    group_by(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label, condition_label, condition_short, replicate, empirical_null_method) %>%
    summarise(
      hc_p_threshold = first(hc_p_threshold),
      hbfss_cutoff = first(hbfss_cutoff),
      hc_valid = first(hc_valid),
      .groups = "drop"
    ) %>%
    group_by(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label, condition_label, condition_short, empirical_null_method) %>%
    summarise(
      n_replicates = n(),
      hc_valid_fraction = mean(as.numeric(hc_valid), na.rm = TRUE),
      hc_p_threshold_min = min_finite(ifelse(hc_valid, hc_p_threshold, NA_real_)),
      hc_p_threshold_median = median_finite(ifelse(hc_valid, hc_p_threshold, NA_real_)),
      hc_p_threshold_max = max_finite(ifelse(hc_valid, hc_p_threshold, NA_real_)),
      hbfss_cutoff_median = median_finite(hbfss_cutoff),
      .groups = "drop"
    )
}

comparison_summary <- function(metric_summary, target, baseline_method, test_method) {
  base <- metric_summary %>%
    filter(truth_target == target, method == baseline_method) %>%
    select(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label, condition_label, condition_short,
           baseline_precision = precision_mean,
           baseline_recall = recall_mean,
           baseline_f1 = f1_mean,
           baseline_fdr = fdr_mean,
           baseline_discovery_count = discovery_count_mean,
           baseline_truth_count = truth_count_mean)

  test <- metric_summary %>%
    filter(truth_target == target, method == test_method) %>%
    select(de_fraction, lfc_magnitude, simulation_profile, null_inflation, lfc_label, de_label, inflation_label, condition_label, condition_short,
           test_precision = precision_mean,
           test_recall = recall_mean,
           test_f1 = f1_mean,
           test_fdr = fdr_mean,
           test_discovery_count = discovery_count_mean,
           test_truth_count = truth_count_mean)

  out <- inner_join(
    base,
    test,
    by = c("de_fraction", "lfc_magnitude", "simulation_profile", "null_inflation", "lfc_label", "de_label", "inflation_label", "condition_label", "condition_short")
  )
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

stack_metrics_for_plot <- function(summary_df, metrics, method_subset = NULL) {
  if (!is.null(method_subset)) {
    summary_df <- summary_df[summary_df$method %in% method_subset, , drop = FALSE]
  }
  rows <- lapply(metrics, function(metric_name) {
    data.frame(
      de_fraction = summary_df$de_fraction,
      lfc_magnitude = summary_df$lfc_magnitude,
      simulation_profile = summary_df$simulation_profile,
      null_inflation = summary_df$null_inflation,
      lfc_label = summary_df$lfc_label,
      de_label = summary_df$de_label,
      inflation_label = summary_df$inflation_label,
      condition_label = summary_df$condition_label,
      condition_short = summary_df$condition_short,
      method = summary_df$method,
      truth_target = summary_df$truth_target,
      metric = metric_name,
      value = summary_df[[paste0(metric_name, "_median")]],
      q25 = summary_df[[paste0(metric_name, "_q25")]],
      q75 = summary_df[[paste0(metric_name, "_q75")]],
      truth_count_mean = summary_df$truth_count_mean,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

stack_deltas_for_heatmap <- function(delta_df) {
  bind_rows(
    data.frame(condition_label = delta_df$condition_label, metric = "Delta F1", value = delta_df$delta_f1, stringsAsFactors = FALSE),
    data.frame(condition_label = delta_df$condition_label, metric = "Delta Recall", value = delta_df$delta_recall, stringsAsFactors = FALSE),
    data.frame(condition_label = delta_df$condition_label, metric = "Delta Precision", value = delta_df$delta_precision, stringsAsFactors = FALSE),
    data.frame(condition_label = delta_df$condition_label, metric = "Delta FDR", value = delta_df$delta_fdr, stringsAsFactors = FALSE)
  )
}

plot_all_de_delta_heatmap <- function(delta_df) {
  if (nrow(delta_df) == 0L) return(NULL)
  plot_df <- stack_deltas_for_heatmap(delta_df)
  plot_df$condition_label <- factor(plot_df$condition_label, levels = rev(unique(delta_df$condition_label)))
  plot_df$metric <- factor(plot_df$metric, levels = c("Delta F1", "Delta Recall", "Delta Precision", "Delta FDR"))
  plot_df <- plot_df %>%
    group_by(metric) %>%
    mutate(metric_max_abs = max(abs(value), na.rm = TRUE)) %>%
    ungroup()
  plot_df$fill_value <- ifelse(is.finite(plot_df$metric_max_abs) & plot_df$metric_max_abs > 0, plot_df$value / plot_df$metric_max_abs, 0)
  plot_df$label <- percent_delta_label(plot_df$value)

  ggplot(plot_df, aes(x = metric, y = condition_label, fill = fill_value)) +
    geom_tile(color = "white", linewidth = 0.30) +
    geom_text(aes(label = label), size = 1.9) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-1, 1), name = "Direction") +
    labs(
      title = "Overall HBFSS benefit versus DESeq2 BH",
      subtitle = "Positive percentages favor HBFSS; color is scaled within each metric",
      x = NULL,
      y = NULL
    ) +
    manuscript_theme() +
    theme(axis.text.y = element_text(size = 5.3))
}

plot_head_to_head <- function(summary_df, methods, target, metrics, title, subtitle = NULL) {
  plot_df <- stack_metrics_for_plot(summary_df, metrics = metrics, method_subset = methods)
  plot_df <- plot_df[plot_df$truth_target == target, , drop = FALSE]
  plot_df <- plot_df[!is.na(plot_df$truth_count_mean) & plot_df$truth_count_mean > 0, , drop = FALSE]
  if (nrow(plot_df) == 0L) return(NULL)

  metric_label_lookup <- c(f1 = "F1", recall = "Recall", precision = "Precision", fdr = "Observed FDR", discovery_count = "Discovery count")
  plot_df$metric <- factor(plot_df$metric, levels = metrics, labels = unname(metric_label_lookup[metrics]))
  plot_df$method_label <- factor(method_label_map[plot_df$method], levels = method_label_map[methods])

  dodge <- position_dodge(width = 0.55)

  ggplot(plot_df, aes(x = de_label, y = value, color = method_label, group = method_label)) +
    geom_crossbar(aes(ymin = q25, ymax = q75), position = dodge, width = 0.32, linewidth = 0.28, fatten = 1.2, alpha = 0.95) +
    facet_grid(metric ~ inflation_label + lfc_label, scales = "free_y") +
    scale_color_manual(values = setNames(unname(method_color_map[methods]), method_label_map[methods]), breaks = method_label_map[methods], name = NULL) +
    scale_y_continuous(labels = percent_label) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Differential-expression fraction",
      y = "Median with IQR"
    ) +
    manuscript_theme() +
    theme(legend.position = "bottom")
}

plot_qc_heatmap <- function(qc_df, value_col, title, subtitle = NULL, fill_label = NULL) {
  if (nrow(qc_df) == 0L) return(NULL)
  plot_df <- qc_df
  plot_df$value <- plot_df[[value_col]]
  plot_df$de_label <- factor(plot_df$de_label, levels = c("5% DE", "10% DE"))
  plot_df$inflation_label <- factor(plot_df$inflation_label, levels = c("No null inflation", "10% null inflation"))
  plot_df$label <- if (identical(value_col, "hc_valid_fraction")) percent_label(plot_df$value) else formatC(plot_df$value, format = "f", digits = 2)

  p <- ggplot(plot_df, aes(x = de_label, y = inflation_label, fill = value)) +
    geom_tile(color = "white", linewidth = 0.30) +
    geom_text(aes(label = label), size = 2.0) +
    facet_wrap(~ lfc_label, nrow = 1) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Differential-expression fraction",
      y = NULL,
      fill = fill_label
    ) +
    manuscript_theme()

  if (identical(value_col, "hc_valid_fraction")) {
    p <- p + scale_fill_gradient(limits = c(0, 1), low = "white", high = "#2166AC")
  } else {
    p <- p + scale_fill_gradient(low = "white", high = "#54278F")
  }
  p
}

run_sequence_simulation_validation <- function() {
  set.seed(simulation_seed)

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

  log_message("Running SEQUENCE simulation validation: ", total_runs, " planned replicates")
  log_message(
    "Simulation config: features=", simulation_n_features,
    ", samples/group=", simulation_n_samples_per_group,
    ", reps/grid=", simulation_n_reps,
    ", apeglm_shrinkage=", simulation_use_apeglm_shrinkage
  )

  simulation_start <- Sys.time()
  metric_rows <- list()
  failure_rows <- list()
  feature_rows <- list()
  checkpoint_path <- file.path(simulation_dir, "Simulation_RunMetrics_Checkpoint.csv")
  run_index <- 0L

  for (grid_i in seq_len(nrow(sim_grid))) {
    grid_row <- sim_grid[grid_i, , drop = FALSE]

    for (rep_i in seq_len(simulation_n_reps)) {
      run_index <- run_index + 1L

      if (run_index %% simulation_progress_every == 0L || run_index == 1L || run_index == total_runs) {
        elapsed_min <- as.numeric(difftime(Sys.time(), simulation_start, units = "mins"))
        eta_min <- if (run_index > 0L) elapsed_min * (total_runs - run_index) / run_index else NA_real_
        log_message(
          "Simulation progress: ", run_index, "/", total_runs,
          " | elapsed=", signif(elapsed_min, 3), " min",
          " | ETA=", signif(eta_min, 3), " min"
        )
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
      template$inflation_label <- ifelse(template$null_inflation > 0, paste0(template$null_inflation * 100, "% null inflation"), "No null inflation")
      template$condition_short <- condition_short_label(template$inflation_label, template$de_label)
      template$condition_label <- condition_long_label(template$lfc_label, template$inflation_label, template$de_label)

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
        failure_rows[[length(failure_rows) + 1L]] <<- cbind(template, data.frame(error_message = conditionMessage(e), stringsAsFactors = FALSE))
        NULL
      })

      if (is.null(out)) next

      metric_rows[[length(metric_rows) + 1L]] <- build_simulation_metric_rows(out, template)

      if (run_index %% simulation_checkpoint_every == 0L || run_index == total_runs) {
        if (length(metric_rows) > 0L) {
          save_csv(bind_rows(metric_rows), checkpoint_path)
        }
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

  if (nrow(failures) > 0L) {
    log_message("Simulation replicate failure(s) recorded: ", nrow(failures), ". Figures will still be exported from successful replicates before final failure handling.")
  }

  metric_long <- if (length(metric_rows) > 0L) bind_rows(metric_rows) else data.frame()
  if (nrow(metric_long) == 0L) stop("Simulation validation produced no successful replicates. See Simulation_Failures.csv.", call. = FALSE)

  planned_lfc_levels <- simulation_lfc_levels()
  actual_lfc_levels <- unique(as.character(metric_long$lfc_label))
  metric_long$lfc_label <- factor(as.character(metric_long$lfc_label), levels = unique(c(planned_lfc_levels, actual_lfc_levels)))
  metric_long$method <- factor(
    as.character(metric_long$method),
    levels = c("DESeq2_BH", "Empirical_BH", "GreaterAbs", "LessAbs", "HBFSS_raw", "HBFSS_total", "HBFSS_weak_region")
  )

  metric_summary <- summarise_metric_long(metric_long)
  qc_summary <- summarise_qc(metric_long)

  save_csv(metric_long, file.path(simulation_dir, "Simulation_RunMetrics_Long.csv"))
  save_csv(metric_summary, file.path(simulation_dir, "Simulation_Condition_Method_Summary.csv"))
  save_csv(qc_summary, file.path(simulation_dir, "Simulation_HC_QC_Summary.csv"))
  if (file.exists(checkpoint_path)) unlink(checkpoint_path, force = TRUE)

  if (isTRUE(export_simulation_feature_results) && length(feature_rows) > 0L) {
    save_csv(bind_rows(feature_rows), file.path(simulation_dir, "Simulation_FeatureResults.csv"))
  }

  all_de_delta <- comparison_summary(metric_summary, "all_de", "DESeq2_BH", "HBFSS_total")
  strong_delta <- comparison_summary(metric_summary, "strong_de", "GreaterAbs", "HBFSS_raw")
  weak_delta <- comparison_summary(metric_summary, "weak_de", "LessAbs", "HBFSS_weak_region")

  save_csv(all_de_delta, file.path(simulation_dir, "Simulation_Delta_AllDE_HBFSSvsDESeq2.csv"))
  save_csv(strong_delta, file.path(simulation_dir, "Simulation_Delta_Strong_HBFSSvsGreaterAbs.csv"))
  save_csv(weak_delta, file.path(simulation_dir, "Simulation_Delta_Weak_HBFSSvsLessAbs.csv"))

  figure_status <- bind_rows(
    safe_save_figure(
      plot_all_de_delta_heatmap(all_de_delta),
      "Figure_01_AllDE_DeltaHeatmap",
      figure_width_twocol,
      5.8
    ),
    safe_save_figure(
      plot_head_to_head(
        metric_summary,
        methods = c("GreaterAbs", "HBFSS_raw"),
        target = "strong_de",
        metrics = c("f1", "recall", "precision", "fdr"),
        title = "Strong-effect head-to-head comparison",
        subtitle = "Median and IQR across replicate simulations"
      ),
      "Figure_02_StrongEffect_HeadToHead",
      figure_width_twocol,
      5.8
    ),
    safe_save_figure(
      plot_head_to_head(
        metric_summary,
        methods = c("LessAbs", "HBFSS_weak_region"),
        target = "weak_de",
        metrics = c("f1", "recall", "precision"),
        title = "Weak-effect head-to-head comparison",
        subtitle = "HBFSS weak = LessAbs-supported genes that also pass higher criticism"
      ),
      "Figure_03_WeakEffect_HeadToHead",
      figure_width_twocol,
      5.6
    ),
    safe_save_figure(
      plot_qc_heatmap(
        qc_summary,
        value_col = "hc_valid_fraction",
        title = "Higher-criticism threshold validity",
        subtitle = "Fraction of replicates in which HC returned a valid threshold",
        fill_label = "Valid fraction"
      ),
      "Figure_04_HC_Validity",
      figure_width_twocol,
      3.0
    ),
    safe_save_figure(
      plot_qc_heatmap(
        qc_summary,
        value_col = "hbfss_cutoff_median",
        title = "Median HBFSS cutoff across simulation conditions",
        subtitle = "Condition-level median of the HC-derived HBFSS decision threshold",
        fill_label = "Median cutoff"
      ),
      "Figure_05_HBFSS_Cutoff",
      figure_width_twocol,
      3.0
    )
  )
  save_csv(figure_status, file.path(simulation_dir, "Figure_Export_Status.csv"))

  if (nrow(figure_status) == 0L || any(figure_status$status != "saved")) {
    save_csv(figure_status[figure_status$status != "saved", , drop = FALSE], file.path(simulation_dir, "Figure_Export_Failures.csv"))
    stop("One or more simulation figures failed. See Figure_Export_Failures.csv.", call. = FALSE)
  }

  expected_png <- file.path(simulation_dir, paste0(figure_status$figure, ".png"))
  missing_png <- expected_png[!file.exists(expected_png) | file.info(expected_png)$size <= 0]
  if (length(missing_png) > 0L) {
    stop("Simulation figure verification failed. Missing or empty PNG file(s): ", paste(missing_png, collapse = ", "), call. = FALSE)
  }

  writeLines(
    c(
      "SEQUENCE simulation figures",
      paste0("Root figure folder: ", simulation_dir),
      paste0("PNG archive folder: ", figures_png_dir),
      paste0("PDF archive folder: ", figures_pdf_dir),
      paste0("TIFF archive folder: ", figures_tiff_dir),
      "",
      paste0("PNG files: ", paste(basename(expected_png), collapse = "; "))
    ),
    file.path(simulation_dir, "README_FIGURES.txt")
  )

  writeLines(capture.output(sessionInfo()), file.path(simulation_dir, "SessionInfo_Simulation.txt"))

  if (nrow(failures) > 0L && isTRUE(simulation_fail_on_failed_replicates)) {
    stop("Simulation validation had failed replicate(s). Figures were exported, but Git push was stopped. See Simulation_Failures.csv.", call. = FALSE)
  }
  log_message("Simulation successful replicates: ", length(metric_rows), "/", total_runs)
  log_message("Simulation failed replicates: ", nrow(failures))

  invisible(list(metric_long = metric_long, metric_summary = metric_summary, qc_summary = qc_summary, all_de_delta = all_de_delta, strong_delta = strong_delta, weak_delta = weak_delta, figure_status = figure_status))
}

# =============================================================================
# METHODS, MANIFEST, AND GIT
# =============================================================================

write_simulation_methods <- function() {
  lines <- c(
    "# SEQUENCE simulation validation methods",
    "",
    paste0("The standalone simulation generated negative-binomial count matrices with ", simulation_n_features, " features and ", simulation_n_samples_per_group, " control plus ", simulation_n_samples_per_group, " treatment samples per replicate."),
    paste0("Baseline feature means were sampled from a gamma distribution with mean scale centered around ", simulation_base_mean, "; null dispersion was ", simulation_dispersion_null, " and differential-feature dispersion was ", simulation_dispersion_de, "."),
    paste0("The simulation grid used differential-expression fractions of ", paste(simulation_de_fractions, collapse = ", "), ", fixed absolute log2 fold-change magnitudes of ", paste(simulation_lfc_magnitudes, collapse = ", "), ", and a weak-effect mixture sampled uniformly from absolute log2 fold-change ", simulation_weak_lfc_min, " to ", simulation_weak_lfc_max, "."),
    paste0("Null Wald statistic inflation fractions were ", paste(simulation_null_inflation, collapse = ", "), "."),
    paste0("Each grid cell was repeated ", simulation_n_reps, " times."),
    "",
    "For each simulated count matrix, DESeq2 was run with design ~ condition. Simulated differential features were randomly sampled per replicate rather than assigned by row position. Standard DESeq2 calls used Benjamini-Hochberg adjusted p < alpha_standard and absolute apeglm-shrunken log2 fold change >= lfc_boundary. GreaterAbs and LessAbs tests used the same lfc_boundary with alpha_strong and alpha_weak, respectively.",
    "",
    "Empirical-null p-values were fit from the DESeq2 Wald statistics using fdrtool. Higher criticism was applied to the empirical p-value distribution using fdrtool::hc.thresh. HC thresholds that were non-finite, outside (0,1), below hc_invalid_below, or >= hc_invalid_at_or_above were marked invalid.",
    "",
    "HBFSS was calculated exactly as abs(apeglm-shrunken log2 fold change) * -log10(empirical p). The HBFSS cutoff was calculated exactly as -log10(HC p-threshold) * lfc_boundary. No artificial minimum HBFSS floor was added.",
    "The weak-mixture grid contains true differential features below the lfc_boundary by design. The weak-effect HBFSS arm is defined as LessAbs-supported sub-boundary genes that also pass the HC empirical-p threshold, rather than requiring weak effects to cross the strong-effect HBFSS curve.",
    "",
    "Simulation outputs include per-run metric tables, method summary tables, HC quality-control summaries, F1, precision, recall, FDR, discovery-count, weak-region, strong-effect, and HC-valid-fraction figures in both root-level PNG form and archived PNG/PDF folders, figure export status, session info, and this methods note. Per-feature simulation results are optional and disabled by default to keep GitHub pushes small."
  )
  writeLines(lines, file.path(simulation_dir, "Methods_Simulation.txt"))
  invisible(TRUE)
}


fmt_number <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(x) || is.na(x)) return("NA")
  formatC(x, format = "f", digits = digits)
}

fmt_delta <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(x) || is.na(x)) return("NA")
  paste0(ifelse(x >= 0, "+", ""), formatC(x, format = "f", digits = digits))
}

fmt_pct <- function(x, digits = 1) {
  x <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(x) || is.na(x)) return("NA")
  paste0(formatC(100 * x, format = "f", digits = digits), "%")
}

lookup_overall <- function(overall, method_name, target_name) {
  row <- overall[as.character(overall$method) == method_name & as.character(overall$truth_target) == target_name, , drop = FALSE]
  if (nrow(row) == 0L) return(NULL)
  row[1, , drop = FALSE]
}

delta_metric <- function(overall, method_a, method_b, target_name, metric_col) {
  a <- lookup_overall(overall, method_a, target_name)
  b <- lookup_overall(overall, method_b, target_name)
  if (is.null(a) || is.null(b)) return(NA_real_)
  suppressWarnings(as.numeric(a[[metric_col]][1]) - as.numeric(b[[metric_col]][1]))
}

markdown_table <- function(df, columns, headers) {
  if (is.null(df) || nrow(df) == 0L) return(character(0))
  out <- c(
    paste0("|", paste(headers, collapse = "|"), "|"),
    paste0("|", paste(rep("---", length(headers)), collapse = "|"), "|")
  )
  for (i in seq_len(nrow(df))) {
    values <- vapply(columns, function(col) as.character(df[[col]][i]), character(1))
    out <- c(out, paste0("|", paste(values, collapse = "|"), "|"))
  }
  out
}

write_simulation_interpretation <- function(simulation_result) {
  metric_summary <- simulation_result$metric_summary
  qc_summary <- simulation_result$qc_summary
  metric_long <- simulation_result$metric_long

  if (is.null(metric_summary) || nrow(metric_summary) == 0L) {
    stop("Cannot write interpretation because Simulation_Condition_Method_Summary is empty.", call. = FALSE)
  }

  overall <- metric_summary %>%
    group_by(method, truth_target) %>%
    summarise(
      n_grid_cells = n(),
      precision_mean = mean_finite(precision_mean),
      recall_mean = mean_finite(recall_mean),
      f1_mean = mean_finite(f1_mean),
      fdr_mean = mean_finite(fdr_mean),
      discovery_count_mean = mean_finite(discovery_count_mean),
      truth_count_mean = mean_finite(truth_count_mean),
      hc_valid_fraction = mean_finite(hc_valid_fraction),
      .groups = "drop"
    ) %>%
    arrange(truth_target, desc(f1_mean), method)

  overall_print <- overall %>%
    mutate(
      precision = vapply(precision_mean, fmt_number, character(1), digits = 3),
      recall = vapply(recall_mean, fmt_number, character(1), digits = 3),
      F1 = vapply(f1_mean, fmt_number, character(1), digits = 3),
      FDR = vapply(fdr_mean, fmt_number, character(1), digits = 3),
      discoveries = vapply(discovery_count_mean, fmt_number, character(1), digits = 1),
      truth = vapply(truth_count_mean, fmt_number, character(1), digits = 1)
    ) %>%
    select(method, truth_target, precision, recall, F1, FDR, discoveries, truth)

  save_csv(overall, file.path(simulation_dir, "Simulation_Interpretation_OverallMetrics.csv"))

  best_all <- overall[overall$truth_target == "all_de" & is.finite(overall$f1_mean), , drop = FALSE]
  best_all <- best_all[order(-best_all$f1_mean), , drop = FALSE]
  best_all_method <- if (nrow(best_all) > 0L) as.character(best_all$method[1]) else "NA"

  hbfss_total <- lookup_overall(overall, "HBFSS_total", "all_de")
  deseq2 <- lookup_overall(overall, "DESeq2_BH", "all_de")
  empirical_bh <- lookup_overall(overall, "Empirical_BH", "all_de")
  hbfss_raw <- lookup_overall(overall, "HBFSS_raw", "all_de")

  strong_hbfss <- lookup_overall(overall, "HBFSS_raw", "strong_de")
  strong_greater <- lookup_overall(overall, "GreaterAbs", "strong_de")
  weak_hbfss <- lookup_overall(overall, "HBFSS_weak_region", "weak_de")
  weak_less <- lookup_overall(overall, "LessAbs", "weak_de")

  hc_overall_fraction <- mean_finite(qc_summary$hc_valid_fraction)
  hc_min_fraction <- min_finite(qc_summary$hc_valid_fraction)
  hc_median_threshold <- median_finite(qc_summary$hc_p_threshold_median)
  hc_median_cutoff <- median_finite(qc_summary$hbfss_cutoff_median)

  failed_path <- file.path(simulation_dir, "Simulation_Failures.csv")
  failures <- if (file.exists(failed_path)) tryCatch(read.csv(failed_path), error = function(e) data.frame()) else data.frame()
  failed_count <- if (is.data.frame(failures)) nrow(failures) else 0L
  total_success <- length(unique(paste(metric_long$de_fraction, metric_long$lfc_magnitude, metric_long$simulation_profile, metric_long$null_inflation, metric_long$replicate, sep = "|")))

  lines <- c(
    "# SEQUENCE simulation validation interpretation",
    "",
    "## What this simulation tests",
    "This standalone simulation asks whether the SEQUENCE HBFSS framework behaves as intended when the truth is known. It generates negative-binomial count data, randomly assigns true differential features, runs DESeq2, estimates empirical-null p-values from Wald statistics, applies higher criticism, calculates HBFSS exactly as `abs(apeglm-shrunken log2 fold change) * -log10(empirical p)`, and compares calls against the known simulated truth.",
    "",
    "## How to read the metrics",
    "Precision is the fraction of called features that are truly differential. Recall is the fraction of true differential features recovered. F1 is the harmonic mean of precision and recall, so it rewards methods that balance both. FDR is the observed false-discovery proportion among called features; lower is better. Discovery count shows how many features each method calls. The HC valid-fraction figure shows how often higher criticism produced a usable threshold under each simulation condition.",
    "",
    "## Main automated summary",
    paste0("Successful simulation grid/replicate runs represented in the metric table: ", total_success, "."),
    paste0("Failed replicate records: ", failed_count, "."),
    paste0("Best average all-DE F1 method across exported grid summaries: ", best_all_method, "."),
    paste0("Mean HC valid-threshold fraction across conditions: ", fmt_pct(hc_overall_fraction), "; minimum condition-level HC valid fraction: ", fmt_pct(hc_min_fraction), "."),
    paste0("Median condition-level HC p-threshold: ", fmt_number(hc_median_threshold, 4), "; median HBFSS cutoff: ", fmt_number(hc_median_cutoff, 4), "."),
    "",
    "## Direct all-DE comparison: HBFSS_total versus DESeq2_BH",
    paste0("Delta F1: ", fmt_delta(delta_metric(overall, "HBFSS_total", "DESeq2_BH", "all_de", "f1_mean")), "."),
    paste0("Delta recall: ", fmt_delta(delta_metric(overall, "HBFSS_total", "DESeq2_BH", "all_de", "recall_mean")), "."),
    paste0("Delta precision: ", fmt_delta(delta_metric(overall, "HBFSS_total", "DESeq2_BH", "all_de", "precision_mean")), "."),
    paste0("Delta FDR: ", fmt_delta(delta_metric(overall, "HBFSS_total", "DESeq2_BH", "all_de", "fdr_mean")), "."),
    paste0("Delta mean discovery count: ", fmt_delta(delta_metric(overall, "HBFSS_total", "DESeq2_BH", "all_de", "discovery_count_mean"), digits = 1), "."),
    "",
    "Interpretation: positive delta F1 or recall means HBFSS_total improved balance or sensitivity compared with the standard DESeq2 BH workflow. A positive delta FDR means that gain came with more false discoveries; a negative delta FDR means HBFSS_total was cleaner on average.",
    "",
    "## Strong-effect comparison",
    paste0("GreaterAbs average F1: ", if (is.null(strong_greater)) "NA" else fmt_number(strong_greater$f1_mean), "; HBFSS_raw strong-effect average F1: ", if (is.null(strong_hbfss)) "NA" else fmt_number(strong_hbfss$f1_mean), "."),
    paste0("GreaterAbs average recall: ", if (is.null(strong_greater)) "NA" else fmt_number(strong_greater$recall_mean), "; HBFSS_raw strong-effect average recall: ", if (is.null(strong_hbfss)) "NA" else fmt_number(strong_hbfss$recall_mean), "."),
    "",
    "Interpretation: this is the validation arm for large-effect discoveries. It directly tests whether the HBFSS boundary recovers strong effects compared with DESeq2's `greaterAbs` alternative-hypothesis test.",
    "",
    "## Weak-region comparison",
    paste0("LessAbs average F1: ", if (is.null(weak_less)) "NA" else fmt_number(weak_less$f1_mean), "; HBFSS_weak_region average F1: ", if (is.null(weak_hbfss)) "NA" else fmt_number(weak_hbfss$f1_mean), "."),
    paste0("LessAbs average recall: ", if (is.null(weak_less)) "NA" else fmt_number(weak_less$recall_mean), "; HBFSS_weak_region average recall: ", if (is.null(weak_hbfss)) "NA" else fmt_number(weak_hbfss$recall_mean), "."),
    "",
    "Interpretation: this is the validation arm for sub-boundary effects. The weak-mixture grid intentionally creates true effects below the log2 fold-change boundary, so strong-effect truth counts are zero there by design and strong-effect boxes are omitted for those cells.",
    "",
    "## Overall method table",
    markdown_table(overall_print, c("method", "truth_target", "precision", "recall", "F1", "FDR", "discoveries", "truth"), c("Method", "Truth target", "Precision", "Recall", "F1", "FDR", "Mean discoveries", "Mean truth count")),
    "",
    "## Figure guide",
    "`Figure_01_AllDE_DeltaHeatmap` is the primary overall result. It directly shows HBFSS total minus DESeq2 BH for F1, recall, precision, observed FDR, and discoveries. Positive values favor HBFSS except for FDR, where positive values indicate a false-discovery cost.",
    "`Figure_02_StrongEffect_HeadToHead` compares GreaterAbs with HBFSS raw for strong effects using median and IQR values across replicate simulations.",
    "`Figure_03_WeakEffect_HeadToHead` compares LessAbs with the repaired HBFSS weak classifier, defined as LessAbs-supported genes that also pass the HC empirical-p threshold.",
    "`Figure_04_HC_Validity` and `Figure_05_HBFSS_Cutoff` are QC figures showing threshold availability and cutoff stability.",
    "",
    "## Submission note",
    "For manuscript reporting, cite Figure 01 first, then Figure 02 and Figure 03 for strong- and weak-effect head-to-head validation. Use the HC validity and cutoff figures as QC support. Do not interpret F1 alone; report recall, precision, and FDR tradeoffs together."
  )

  writeLines(lines, file.path(simulation_dir, "Interpretation_Simulation.md"))
  writeLines(lines, file.path(simulation_dir, "Interpretation_Simulation.txt"))
  log_message("Simulation interpretation written: ", file.path(simulation_dir, "Interpretation_Simulation.md"))
  invisible(TRUE)
}


write_manifest <- function() {
  exported_files <- list.files(simulation_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  exported_files <- exported_files[file.exists(exported_files)]
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  manifest <- data.frame(
    file = sub(paste0("^", root, "/?"), "", normalizePath(exported_files, winslash = "/", mustWork = FALSE)),
    size_bytes = file.info(exported_files)$size,
    stringsAsFactors = FALSE
  )
  save_csv(manifest, file.path(simulation_dir, "Manifest_Simulation.csv"))
  manifest
}

path_relative_to_repo <- function(path) {
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  abs_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  if (startsWith(abs_path, prefix)) return(substr(abs_path, nchar(prefix) + 1L, nchar(abs_path)))
  abs_path
}

is_under_repo <- function(path) {
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  abs_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  startsWith(abs_path, paste0(root, "/")) || identical(abs_path, root)
}

git_command <- function(args, allow_failure = FALSE, echo_output = TRUE) {
  if (Sys.which("git") == "") stop("Git executable was not found on PATH.", call. = FALSE)
  cmd_display <- paste("git", "-C", repo_root, paste(args, collapse = " "))
  if (isTRUE(echo_output)) log_message("$ ", cmd_display)
  out <- suppressWarnings(system2("git", args = c("-C", repo_root, args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (isTRUE(echo_output) && length(out) > 0L) log_message(paste(out, collapse = "\n"))
  if (status != 0L && !isTRUE(allow_failure)) {
    stop("Git command failed with status ", status, ": ", cmd_display, "\n", paste(out, collapse = "\n"), call. = FALSE)
  }
  list(status = as.integer(status), output = out)
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
  if (length(missing) == 0L) return(invisible(FALSE))

  prefix <- if (file.exists(ignore_path) && file.info(ignore_path)$size > 0L) "\n" else ""
  cat(prefix, paste(missing, collapse = "\n"), "\n", file = ignore_path, append = TRUE, sep = "")
  log_message("Updated .gitignore with: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

git_remove_cached_excluded_files <- function() {
  patterns <- unique(git_exclude_patterns[nzchar(git_exclude_patterns)])
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
  if (length(staged) == 0L) return(character(0))

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
    ensure_gitignore_patterns(oversized)
    git_command(c("add", "--", ".gitignore"), allow_failure = FALSE)
    log_message("Unstaged oversized file(s): ", paste(oversized, collapse = ", "))
  }
  oversized
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
  if (!identical(git_root, expected_root)) stop("Git root mismatch. repo_root is ", expected_root, " but Git root is ", git_root, call. = FALSE)

  remote_url <- git_output_first_line(c("remote", "get-url", git_remote_name), allow_failure = TRUE)
  if (is.na(remote_url) || !nzchar(remote_url)) stop("Git remote '", git_remote_name, "' is not configured.", call. = FALSE)

  current_branch <- git_output_first_line(c("rev-parse", "--abbrev-ref", "HEAD"), allow_failure = FALSE)
  if (identical(current_branch, "HEAD") || is.na(current_branch) || !nzchar(current_branch)) stop("Git is in detached HEAD state.", call. = FALSE)
  target_branch <- if (!is.na(git_branch_name) && nzchar(git_branch_name)) git_branch_name else current_branch

  ensure_gitignore_patterns(git_exclude_patterns)

  git_command(c("fetch", git_remote_name), allow_failure = FALSE)
  upstream <- git_output_first_line(c("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"), allow_failure = TRUE)
  if (is.na(upstream) || !nzchar(upstream)) {
    remote_branch_ref <- paste0(git_remote_name, "/", target_branch)
    remote_branch_sha <- git_output_first_line(c("rev-parse", "--verify", remote_branch_ref), allow_failure = TRUE)
    if (!is.na(remote_branch_sha) && nzchar(remote_branch_sha)) upstream <- remote_branch_ref
  }

  if (isTRUE(git_pull_rebase_before_push) && !is.na(upstream) && nzchar(upstream)) {
    git_command(c("pull", "--rebase", "--autostash", git_remote_name, target_branch), allow_failure = FALSE)
  }

  git_remove_cached_excluded_files()

  stage_paths <- c(".gitignore", path_relative_to_repo(simulation_dir))
  if (isTRUE(git_stage_pipeline_script)) {
    sp <- script_path()
    if (!is.na(sp) && file.exists(sp) && is_under_repo(sp)) {
      stage_paths <- c(path_relative_to_repo(sp), stage_paths)
    } else {
      log_message("Pipeline script was not staged because it is not inside the Git repository.")
    }
  }
  stage_paths <- unique(stage_paths[nzchar(stage_paths)])

  git_command(c("add", "--", stage_paths), allow_failure = FALSE)
  git_remove_cached_excluded_files()
  git_unstage_oversized_files()

  if (!git_has_staged_changes()) {
    log_message("No changed simulation files to commit.")
    if (isTRUE(git_allow_no_change_success)) return(invisible(FALSE))
    stop("No changed files to commit and git_allow_no_change_success is FALSE.", call. = FALSE)
  }

  git_command(c("status", "--short"), allow_failure = FALSE)
  git_command(c("commit", "-m", git_commit_message), allow_failure = FALSE)

  final_upstream <- git_output_first_line(c("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"), allow_failure = TRUE)
  if (!is.na(final_upstream) && nzchar(final_upstream) && identical(target_branch, current_branch)) {
    git_command(c("push"), allow_failure = FALSE)
  } else {
    git_command(c("push", "-u", git_remote_name, paste0("HEAD:", target_branch)), allow_failure = FALSE)
  }

  log_message("GitHub push complete.")
  invisible(TRUE)
}

# =============================================================================
# EXECUTION
# =============================================================================

main <- function() {
  log_message("SEQUENCE standalone simulation started.")
  log_message("Repository root: ", repo_root)
  log_message("Simulation output directory: ", simulation_dir)

  simulation_result <- run_sequence_simulation_validation()
  write_simulation_methods()
  write_simulation_interpretation(simulation_result)
  manifest <- write_manifest()
  cleanup_rplots_pdf()
  git_commit_and_push()

  cat("\n=====================================================\n")
  cat("SEQUENCE standalone simulation complete.\n")
  cat("Output directory: ", simulation_dir, "\n", sep = "")
  cat("Exported files: ", nrow(manifest), "\n", sep = "")
  cat("=====================================================\n\n")

  print(simulation_result$metric_summary)
  invisible(simulation_result)
}

tryCatch(
  main(),
  error = function(e) {
    log_message("SIMULATION FAILED - DO NOT PUSH: ", conditionMessage(e))
    try(writeLines(capture.output(sessionInfo()), file.path(simulation_dir, "SessionInfo_Simulation_Failed.txt")), silent = TRUE)
    try(cleanup_rplots_pdf(), silent = TRUE)
    if (!interactive()) quit(save = "no", status = 1L, runLast = FALSE)
    stop(e)
  }
)
