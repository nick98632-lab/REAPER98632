#!/usr/bin/env Rscript
# =============================================================================
# SEQUENCE MANUSCRIPT PIPELINE
# Likelihood-stratified moderation of DESeq2 dispersion estimates for WTTS-Seq
# polyadenylation sites (PASs), with eigenvector splitting (EVS), apeglm,
# empirical-null calibration, higher criticism, HBFSS, and 3'aTWAS overlap.
# =============================================================================
#
# RUN
#   Rscript SEQUENCE_FINAL_MANUSCRIPT.R
#   Optional environment variables:
#     SEQUENCE_CORES      cores for the cutoff scan and stability refits (default 2)
#     SEQUENCE_REPO_ROOT  folder holding WTTS-Seq_2022.2_DE_raw_read_numbers.csv
#     SEQUENCE_MULTI_ROOT output folder (default <repo>/exports/sequence_final_manuscript)
#
# The script runs as a master process that relaunches itself once per cutoff
# method (Fixed 5,000, then Likelihood k*), then writes the cross-method
# comparison and one ZIP of everything.
#
# ANALYSIS UNIT
# Each OrigID is analyzed as an individual PAS. Gene symbols are annotation
# only and are never used to collapse PASs before testing.
#
# COMPARISONS (5 RT vs 5 ZT samples each)
#   RT0 vs ZT6, RT2 vs ZT8, RT4 vs ZT10, RT8 vs ZT14
#
# ANALYSIS VIEWS PER COMPARISON
#   Original (No EVS), NormEVS Lead, NormEVS Rem, RawEVS Lead, RawEVS Rem
#
# EIGENVECTOR SPLITTING
# Within each condition, prcomp (centred, unscaled) is applied to the PAS-by-
# sample matrix: DESeq2 median-of-ratios normalized counts for NormEVS, raw
# counts for RawEVS. PASs are ranked by |PC1 loading| and the top k are taken
# from each condition. PASs in both top-k sets are Joint, in one set Disjoint;
# Joint plus Disjoint form the Leading Edge (LE) and all other PASs the
# Remainder (REM). LE and REM are separate datasets: DESeq2 receives their raw
# counts and estimates its own size factors, dispersion trend and prior for
# each. The same k is used for NormEVS and RawEVS.
#
# CHOICE OF k
#   Fixed 5,000    prespecified comparator.
#   Likelihood k*  per comparison, the k that maximizes the negative-binomial
#                  likelihood of the counts when the dispersion prior (trend
#                  and width) may differ between LE and REM. See the section
#                  LIKELIHOOD-STRATIFIED MODERATION below. k* is fixed before
#                  any differential test; no p-value or discovery count is used.
#
# DIFFERENTIAL EXPRESSION
# DESeq2 design ~ condition (trt = RT relative to untrt = ZT).
#   Standard  ordinary Wald, BH padj < 0.10 and |apeglm LFC| >= 1
#   Strong    greaterAbs(lfcThreshold = 1), BH padj < 0.10 and |apeglm LFC| >= 1
#   Weak      lessAbs(lfcThreshold = 1), BH padj < 0.20, |apeglm LFC| < 1, and HBFSS
#   HBFSS     |apeglm LFC| x [-log10(empirical p)] > Htau
# Empirical p-values come from fdrtool (normal empirical null, cutoff 'fndr')
# applied to the Wald statistics. Higher criticism is maximized over the lowest
# 10% of ordered empirical p-values (alpha0 = 0.10) and must be positive;
# Htau = -log10(HCp) x 1.
#
# FIGURES (per cutoff method in <method>/Combined_Figures; cross-method in
# Cross_Method_Comparison)
#   Figure 1   choosing k by maximum likelihood            (Likelihood k* run)
#   Figure 2   what the split does to DESeq2's dispersions (each run)
#   Figure S1  precision and stability of k*               (Likelihood k* run)
#   Figure S3  Fixed 5,000 vs Likelihood k*                (master)
#   Figure S2  FDR / FNDR calibration: SEQUENCE_FDR_FNDR_CALIBRATION.R
# Volcano, PCA/EVS evidence, HBFSS support and 3'aTWAS figures as before.
# =============================================================================

PIPELINE_BUILD <- "SEQUENCE_FINAL_MANUSCRIPT_2026-10-31"


# =============================================================================
# CUTOFF METHODS, SHARED FIGURE STYLE, AND THE CROSS-METHOD FIGURE (S3)
# =============================================================================

COMPARISON_NAMES <- c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14")
METHOD_LEVELS <- c("Fixed 5,000", "Likelihood k*")
METHOD_SLUGS <- c("Fixed 5,000" = "Fixed_5000", "Likelihood k*" = "Likelihood_k")
METHOD_KEYS <- c("Fixed 5,000" = "fixed5000", "Likelihood k*" = "likelihood")
METHOD_COLS <- c("Fixed 5,000" = "#9A9A9A", "Likelihood k*" = "#2E7D5B")
ORIGINAL_LABEL <- "Original (No EVS)"
SET_COLS <- c("Original (No EVS)" = "#2B2B2B", "Leading Edge" = "#1F6FB2", "Remainder" = "#C0392B")
SET_TREND_COLS <- c("Original (No EVS)" = "#2B2B2B", "Leading Edge" = "#0B3C6E", "Remainder" = "#7B1F14")
COMPOSITION_COLS <- c("Joint (top-k in both)" = "#0B3C6E", "Disjoint (top-k in one)" = "#BFD3E8")

comp_lab <- function(cmp = COMPARISON_NAMES) stats::setNames(sub("_", " / ", cmp), cmp)

fig_theme <- function(base = 7.5) {
  theme_minimal(base_size = base, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
      axis.line = element_line(linewidth = 0.3, colour = "grey35"),
      axis.ticks = element_line(linewidth = 0.3, colour = "grey35"),
      axis.text = element_text(size = base - 1.3, colour = "grey20"),
      axis.title = element_text(size = base - 0.4, colour = "grey10"),
      strip.background = element_rect(fill = "grey94", colour = NA),
      strip.text = element_text(face = "bold", size = base - 0.3, colour = "grey10",
                                margin = margin(2.5, 2, 2.5, 2)),
      panel.spacing = unit(4, "mm"),
      legend.title = element_blank(), legend.text = element_text(size = base - 1),
      legend.key.size = unit(3, "mm"), legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
      plot.title = element_text(size = base + 0.6, face = "bold", margin = margin(b = 1.5)),
      plot.subtitle = element_text(size = base - 1, colour = "grey35", margin = margin(b = 4)),
      plot.title.position = "plot", plot.tag = element_text(size = base + 4, face = "bold"),
      plot.margin = margin(4, 5, 4, 4))
}

# PDF (vector) and 600-dpi PNG of a ggplot, patchwork or grob.
save_figure <- function(p, out_dir, stem, w_mm, h_mm) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  w <- w_mm / 25.4; h <- h_mm / 25.4
  draw <- function() if (inherits(p, "ggplot")) print(p) else { grid::grid.newpage(); grid::grid.draw(p) }
  cairo <- isTRUE(capabilities("cairo"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  if (cairo) grDevices::cairo_pdf(pdf_path, width = w, height = h) else
    grDevices::pdf(pdf_path, width = w, height = h, useDingbats = FALSE)
  tryCatch(draw(), finally = grDevices::dev.off())
  grDevices::png(file.path(out_dir, paste0(stem, ".png")), width = w, height = h, units = "in", res = 600,
                 bg = "white", type = if (cairo) "cairo" else getOption("bitmapType"))
  tryCatch(draw(), finally = grDevices::dev.off())
  invisible(TRUE)
}

# rows: list; each element is a ggplot or list(plots = list(...), widths = c(...)).
assemble_rows <- function(rows, heights, title) {
  title_theme <- theme(plot.title = element_text(size = 10, face = "bold", margin = margin(b = 6)))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    built <- lapply(rows, function(r) if (inherits(r, "ggplot")) r else
      patchwork::wrap_plots(r$plots, nrow = 1, widths = r$widths))
    return(patchwork::wrap_plots(built, ncol = 1, heights = heights) +
             patchwork::plot_annotation(title = title, theme = title_theme))
  }
  built <- lapply(rows, function(r) if (inherits(r, "ggplot")) ggplot2::ggplotGrob(r) else
    gridExtra::arrangeGrob(grobs = lapply(r$plots, ggplot2::ggplotGrob), nrow = 1, widths = r$widths))
  gridExtra::arrangeGrob(grobs = built, ncol = 1, heights = heights,
                         top = grid::textGrob(title, x = 0.01, hjust = 0,
                                              gp = grid::gpar(fontsize = 10, fontface = "bold")))
}

master_plotting_ready <- function() {
  ok <- requireNamespace("ggplot2", quietly = TRUE)
  if (!ok) warning("ggplot2 is not installed: cross-method figure skipped (tables still written).")
  if (ok) suppressPackageStartupMessages(library(ggplot2))
  ok
}

try_fig <- function(expr, what) {
  tryCatch(expr, error = function(e) {
    warning("Figure skipped (", what, "): ", conditionMessage(e), call. = FALSE)
    invisible(NULL)
  })
}

read_run_table <- function(multi_root, method_label, file) {
  f <- file.path(multi_root, METHOD_SLUGS[[method_label]], "Summary_Tables", file)
  if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) else NULL
}

# Supplementary Figure S3: (a) where each cutoff sits on the likelihood,
# (b) pooled prior width at each cutoff, (c) discoveries (reference only).
save_figure_S3 <- function(multi_root, out_dir, discoveries) {
  if (!master_plotting_ready()) return(invisible(FALSE))
  cmp <- COMPARISON_NAMES
  k_lab <- function(x) ifelse(is.na(x), NA_character_, ifelse(x == 0, "0", paste0(x / 1000, "k")))
  scan <- read_run_table(multi_root, "Likelihood k*", "Table_Likelihood_Cutoff_Scan.csv")
  sel <- read_run_table(multi_root, "Likelihood k*", "Table_Likelihood_Cutoff_Selected.csv")
  if (is.null(scan) || is.null(sel)) stop("Likelihood scan tables are missing.")
  scan$Comparison <- factor(scan$Comparison, levels = cmp)
  marks <- rbind(data.frame(Comparison = sel$Comparison, k = sel$k_star, Method = "Likelihood k*"),
                 data.frame(Comparison = sel$Comparison, k = 5000L, Method = "Fixed 5,000"))
  marks$y <- scan$delta_loglik_vs_Original[match(paste(marks$Comparison, marks$k), paste(scan$Comparison, scan$k))]
  marks <- marks[is.finite(marks$y), ]
  marks$Comparison <- factor(marks$Comparison, levels = cmp)
  marks$Method <- factor(marks$Method, levels = METHOD_LEVELS)
  pa <- ggplot(scan, aes(k, delta_loglik_vs_Original)) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.3, colour = "grey45") +
    geom_line(linewidth = 0.55, colour = "grey15") +
    geom_point(data = marks, aes(y = y, colour = Method, shape = Method), size = 2.4) +
    facet_wrap(~ Comparison, nrow = 1, scales = "free_y", labeller = labeller(Comparison = comp_lab())) +
    scale_colour_manual(values = METHOD_COLS) +
    scale_shape_manual(values = c("Fixed 5,000" = 15, "Likelihood k*" = 16)) +
    scale_x_continuous(labels = k_lab, breaks = scales::breaks_pretty(4)) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "k (top-|PC1| PASs per condition)", y = "Log-likelihood gain\nover Original (No EVS)", tag = "a",
         title = "Where each cutoff sits on the likelihood", subtitle = "Higher = the split explains the counts better.") +
    fig_theme()

  pw <- do.call(rbind, lapply(METHOD_LEVELS, function(m) {
    d <- read_run_table(multi_root, m, "Table_Prior_Width_ML.csv")
    if (is.null(d)) NULL else data.frame(Comparison = d$Comparison, Series = m, W = d$W_split,
                                         W0 = d$tau2_Original)
  }))
  pb <- NULL
  if (!is.null(pw) && nrow(pw)) {
    o <- unique(pw[, c("Comparison", "W0")])
    b <- rbind(data.frame(Comparison = o$Comparison, Series = ORIGINAL_LABEL, W = o$W0),
               pw[, c("Comparison", "Series", "W")])
    b$Series <- factor(b$Series, levels = c(ORIGINAL_LABEL, METHOD_LEVELS))
    b$Comparison <- factor(b$Comparison, levels = cmp)
    pb <- ggplot(b, aes(Comparison, W, fill = Series)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.74) +
      scale_fill_manual(values = c("Original (No EVS)" = unname(SET_COLS[ORIGINAL_LABEL]), METHOD_COLS)) +
      scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
      labs(x = NULL, y = "Pooled prior width W", tag = "b", title = "Prior width (lower is better)",
           subtitle = "Pooled Leading Edge + Remainder width; Original (No EVS) for reference.") +
      fig_theme() + theme(panel.grid.major.x = element_blank())
  }

  pc <- NULL
  if (!is.null(discoveries) && nrow(discoveries)) {
    d <- discoveries[discoveries$Track == "NormEVS", ]
    sums <- stats::aggregate(cbind(Std, HBFSS) ~ Cutoff_Method, data = d, FUN = sum)
    o <- discoveries[discoveries$Set == "Original" & discoveries$Cutoff_Method == METHOD_LEVELS[1], ]
    sums <- rbind(data.frame(Cutoff_Method = ORIGINAL_LABEL, Std = sum(o$Std), HBFSS = sum(o$HBFSS)), sums)
    long <- rbind(data.frame(Series = sums$Cutoff_Method, Call = "Standard", n = sums$Std),
                  data.frame(Series = sums$Cutoff_Method, Call = "HBFSS", n = sums$HBFSS))
    long$Series <- factor(long$Series, levels = c(ORIGINAL_LABEL, METHOD_LEVELS))
    long$Call <- factor(long$Call, levels = c("Standard", "HBFSS"))
    pc <- ggplot(long, aes(Call, n, fill = Series)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.74) +
      geom_text(aes(label = format(n, big.mark = ",")), position = position_dodge(width = 0.8),
                vjust = -0.4, size = 1.9) +
      scale_fill_manual(values = c("Original (No EVS)" = unname(SET_COLS[ORIGINAL_LABEL]), METHOD_COLS)) +
      scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.12))) +
      labs(x = NULL, y = "Significant PASs\n(NormEVS LE + REM, all comparisons)", tag = "c",
           title = "Discoveries (reference only)", subtitle = "Not used to choose k.") +
      fig_theme() + theme(panel.grid.major.x = element_blank())
  }
  bottom <- Filter(Negate(is.null), list(pb, pc))
  rows <- list(pa)
  if (length(bottom)) rows[[2]] <- list(plots = bottom, widths = rep(1, length(bottom)))
  fig <- assemble_rows(rows, heights = c(1, 1.1)[seq_along(rows)],
                       title = "Supplementary Figure S3. Fixed 5,000 vs Likelihood k*")
  save_figure(fig, out_dir, "Figure_S3_Fixed_vs_Likelihood", 180, 140)
  invisible(TRUE)
}

# Per cutoff method, EVS track and comparison: final significant PASs
# (Std | Strong | Weak | HBFSS) in the Leading Edge and Remainder, and the
# Original (No EVS) reference. Leading Edge and Remainder are disjoint.
collect_cross_method_discoveries <- function(multi_root) {
  views <- c("NormEVS Lead", "NormEVS Rem", "RawEVS Lead", "RawEVS Rem", ORIGINAL_LABEL)
  flag_cols <- c(Std = "Std", Strong = "Strong", Weak = "Weak", HBFSS = "HBFSS_sig")
  rows <- list()
  for (m in METHOD_LEVELS) {
    for (cn in COMPARISON_NAMES) {
      f <- file.path(multi_root, METHOD_SLUGS[[m]], cn, paste0("Table_", cn, "_Significant_Sites_All_Views.csv"))
      sig <- if (file.exists(f) && file.info(f)$size > 0) {
        tryCatch(utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
      } else data.frame()
      for (v in views) {
        s <- if (nrow(sig) && "Analysis" %in% names(sig)) sig[sig$Analysis == v, , drop = FALSE] else sig[0, , drop = FALSE]
        cnt <- vapply(flag_cols, function(col) {
          if (nrow(s) && col %in% names(s)) sum(as.logical(s[[col]]), na.rm = TRUE) else 0L
        }, integer(1))
        rows[[length(rows) + 1L]] <- data.frame(
          Cutoff_Method = m, Comparison = cn, Analysis = v, Any_significant = nrow(s),
          Std = cnt[["Std"]], Strong = cnt[["Strong"]], Weak = cnt[["Weak"]], HBFSS = cnt[["HBFSS"]])
      }
    }
  }
  out <- do.call(rbind, rows)
  out$Track <- ifelse(grepl("^NormEVS", out$Analysis), "NormEVS",
                      ifelse(grepl("^RawEVS", out$Analysis), "RawEVS", "Original"))
  out$Set <- ifelse(grepl("Lead$", out$Analysis), "Leading Edge",
                    ifelse(grepl("Rem$", out$Analysis), "Remainder", "Original"))
  out
}

# =============================================================================
# MASTER PROCESS: one child run per cutoff method, then cross-method outputs
# =============================================================================

SEQUENCE_CHILD_RUN <- identical(Sys.getenv("SEQUENCE_CHILD_RUN", unset = "0"), "1")

if (!SEQUENCE_CHILD_RUN) {
  script_args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", script_args, value = TRUE)
  if (!length(script_arg)) stop("Run this pipeline with Rscript so the master process can relaunch itself.")
  this_script <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)

  count_name <- "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
  root_candidates <- unique(c(Sys.getenv("SEQUENCE_REPO_ROOT", unset = ""), getwd(), dirname(getwd()),
                              "/root/REAPER98632", dirname(this_script)))
  root_candidates <- root_candidates[nzchar(root_candidates)]
  repo_guess <- NULL
  for (cand in root_candidates) {
    if (file.exists(file.path(cand, count_name)) || file.exists(file.path(cand, "data", count_name))) {
      repo_guess <- normalizePath(cand, winslash = "/", mustWork = TRUE)
      break
    }
  }
  if (is.null(repo_guess)) repo_guess <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

  multi_root <- Sys.getenv("SEQUENCE_MULTI_ROOT",
                           unset = file.path(repo_guess, "exports", "sequence_final_manuscript"))
  if (dir.exists(multi_root)) unlink(multi_root, recursive = TRUE, force = TRUE)
  dir.create(multi_root, recursive = TRUE, showWarnings = FALSE)

  # Likelihood k* is estimated inside its own run; its row is filled afterwards.
  cutoff_manifest <- data.frame(
    Cutoff_Method = unname(METHOD_SLUGS),
    RT0_ZT6 = c(5000L, NA_integer_), RT2_ZT8 = c(5000L, NA_integer_),
    RT4_ZT10 = c(5000L, NA_integer_), RT8_ZT14 = c(5000L, NA_integer_),
    Source = c("Prespecified fixed comparator",
               "Maximum likelihood of the varying-dispersion (Leading Edge / Remainder) model, estimated in its run"))
  utils::write.csv(cutoff_manifest, file.path(multi_root, "Cutoff_Method_Manifest.csv"), row.names = FALSE)

  for (m in METHOD_LEVELS) {
    message("=====================================================")
    message("Starting cutoff method: ", m)
    message("=====================================================")
    status <- system2(
      command = file.path(R.home("bin"), "Rscript"), args = shQuote(this_script),
      env = c("SEQUENCE_CHILD_RUN=1", paste0("SEQUENCE_CUTOFF_METHOD=", METHOD_KEYS[[m]]),
              paste0("SEQUENCE_MULTI_ROOT=", multi_root), paste0("SEQUENCE_REPO_ROOT=", repo_guess)))
    if (!identical(as.integer(status), 0L)) stop("Cutoff-method run failed: ", m, " (exit status ", status, ")")
  }

  cross_dir <- file.path(multi_root, "Cross_Method_Comparison")
  dir.create(cross_dir, recursive = TRUE, showWarnings = FALSE)
  sel <- read_run_table(multi_root, "Likelihood k*", "Table_Likelihood_Cutoff_Selected.csv")
  if (!is.null(sel)) {
    for (cn in intersect(COMPARISON_NAMES, sel$Comparison)) {
      cutoff_manifest[cutoff_manifest$Cutoff_Method == METHOD_SLUGS[["Likelihood k*"]], cn] <-
        as.integer(sel$k_star[sel$Comparison == cn])
    }
    utils::write.csv(cutoff_manifest, file.path(multi_root, "Cutoff_Method_Manifest.csv"), row.names = FALSE)
    utils::write.csv(sel, file.path(cross_dir, "Table_Likelihood_Cutoff_Selected.csv"), row.names = FALSE)
  }
  discoveries <- collect_cross_method_discoveries(multi_root)
  utils::write.csv(discoveries, file.path(cross_dir, "Table_Cross_Method_Discoveries.csv"), row.names = FALSE)
  try_fig(save_figure_S3(multi_root, cross_dir, discoveries), "Figure S3")

  final_zip <- file.path(dirname(multi_root), paste0(basename(multi_root), ".zip"))
  if (file.exists(final_zip)) unlink(final_zip, force = TRUE)
  oldwd <- getwd()
  setwd(dirname(multi_root))
  zip_rel <- list.files(multi_root, recursive = TRUE)
  zip_rel <- zip_rel[!grepl("SEQUENCE_ALL_(FIGURES|TABLES)\\.zip$", zip_rel)]
  utils::zip(zipfile = basename(final_zip), files = file.path(basename(multi_root), zip_rel), flags = "-q")
  setwd(oldwd)
  if (!file.exists(final_zip)) stop("Final ZIP creation failed: ", final_zip)

  cat("\n=====================================================\n")
  cat("All runs complete: Fixed 5,000 and Likelihood k*.\n")
  if (!is.null(sel)) for (i in seq_len(nrow(sel))) cat(sprintf(
    "  %-9s k* = %6d (support %d-%d)  log-lik gain vs Original = %.1f, vs k = 5,000 = %.1f  prior width %.0f%% narrower\n",
    sel$Comparison[i], sel$k_star[i], sel$support_lo[i], sel$support_hi[i], sel$delta_loglik_vs_Original[i],
    sel$delta_loglik_k_star_vs_5000[i], sel$prior_width_reduction_pct[i]))
  tot <- stats::aggregate(cbind(Std, HBFSS) ~ Cutoff_Method, FUN = sum,
                          data = discoveries[discoveries$Track == "NormEVS", ])
  for (i in seq_len(nrow(tot))) cat(sprintf("  NormEVS LE + REM, %-13s Std = %5d  HBFSS = %5d\n",
                                            tot$Cutoff_Method[i], tot$Std[i], tot$HBFSS[i]))
  cat("Final ZIP:\n", normalizePath(final_zip, winslash = "/", mustWork = TRUE), "\n", sep = "")
  cat("=====================================================\n\n")
  quit(save = "no", status = 0L)
}

# Active cutoff method of this child run.
CUTOFF_METHOD_KEY <- Sys.getenv("SEQUENCE_CUTOFF_METHOD", unset = "fixed5000")
CUTOFF_METHOD_LABEL <- names(METHOD_KEYS)[match(CUTOFF_METHOD_KEY, METHOD_KEYS)]
if (is.na(CUTOFF_METHOD_LABEL)) stop("Unknown SEQUENCE_CUTOFF_METHOD: ", CUTOFF_METHOD_KEY)
CUTOFF_METHOD_SLUG <- METHOD_SLUGS[[CUTOFF_METHOD_LABEL]]
CUTOFF_METHOD_BASIS <- c(
  "Fixed 5,000" = "Prespecified fixed top-5,000 PASs per condition",
  "Likelihood k*" = "Per comparison, the k maximizing the negative-binomial likelihood of the counts when the dispersion prior (trend and width) may differ between Leading Edge and Remainder (varying-dispersion model fitted by EM); no p-values used"
)[[CUTOFF_METHOD_LABEL]]


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

# Higher Criticism searches only the lower tail of the ordered empirical-null
# p-value distribution. alpha0 = 0.10 means that the maximum HC score is sought
# only among the lowest 10% of empirical p-values. A non-positive maximum HC
# score yields no HC threshold and therefore no HBFSS discoveries in that view.
HC_ALPHA0 <- 0.10

if (!is.numeric(HC_ALPHA0) || length(HC_ALPHA0) != 1L ||
    !is.finite(HC_ALPHA0) || HC_ALPHA0 <= 0 || HC_ALPHA0 > 1) {
  stop("HC_ALPHA0 must be a single finite number in (0, 1].")
}

# EVS selection size k is comparison-specific: prespecified (Fixed 5,000) or
# estimated by maximum likelihood before any test (Likelihood k*).
EMPIRICAL_EVS_CUTOFF_BASIS <- paste0(
  CUTOFF_METHOD_BASIS,
  "; locked before downstream DE significance testing"
)

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
# HC is maximized only over the lowest HC_ALPHA0 fraction of ordered empirical
# p-values and must attain a strictly positive maximum. HCp is not applied again
# as a second HBFSS significance gate.
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

cutoff_manifest_path <- file.path(
  Sys.getenv("SEQUENCE_MULTI_ROOT", unset = ""),
  "Cutoff_Method_Manifest.csv"
)
if (!nzchar(Sys.getenv("SEQUENCE_MULTI_ROOT", unset = "")) || !file.exists(cutoff_manifest_path)) {
  stop("Cutoff manifest is missing. Run this script normally through the master process; do not launch a child run directly.")
}
cutoff_manifest_runtime <- utils::read.csv(cutoff_manifest_path, stringsAsFactors=FALSE, check.names=FALSE)
required_manifest_methods <- c("Fixed_5000")
if (!all(required_manifest_methods %in% cutoff_manifest_runtime$Cutoff_Method)) {
  stop("Cutoff manifest does not contain all required methods: ", paste(required_manifest_methods, collapse=", "))
}
manifest_row <- function(method) cutoff_manifest_runtime[match(method, cutoff_manifest_runtime$Cutoff_Method),,drop=FALSE]
fixed_row <- manifest_row("Fixed_5000")

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix   = c("R0", "R2", "R4", "R8"),
  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),
  fixed_5000_k    = as.integer(unlist(fixed_row[c("RT0_ZT6","RT2_ZT8","RT4_ZT10","RT8_ZT14")], use.names=FALSE)),
  stringsAsFactors = FALSE
)

get_active_evs_cutoff <- function(comparison_name) {
  idx <- match(as.character(comparison_name), comparison_table$comparison_name)
  if (is.na(idx)) stop("No EVS cutoff is defined for comparison: ", comparison_name)
  k <- switch(
    CUTOFF_METHOD_KEY,
    fixed5000 = comparison_table$fixed_5000_k[idx],
    likelihood = if (is.null(LIKELIHOOD_K[[as.character(comparison_name)]])) NA_integer_ else LIKELIHOOD_K[[as.character(comparison_name)]]
  )
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || !is.finite(k) || k < 1L) {
    stop("Invalid EVS cutoff for ", comparison_name, " under ", CUTOFF_METHOD_LABEL)
  }
  k
}


if (anyDuplicated(comparison_table$comparison_name)) stop("comparison_table contains duplicated comparison names.")
cutoff_cols <- c("fixed_5000_k")
if (any(vapply(comparison_table[cutoff_cols], function(x) any(is.na(x) | x < 1L), logical(1)))) {
  stop("Every comparison must have a positive cutoff for every cutoff method.")
}

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

multi_output_root <- Sys.getenv(
  "SEQUENCE_MULTI_ROOT",
  unset = file.path(repo_root, "exports", "sequence_final_manuscript")
)
output_dir <- file.path(multi_output_root, CUTOFF_METHOD_SLUG)

if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE, force = TRUE)

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



resolve_top_n_cutoff <- function(sorted_values_desc, top_n) {
  n_total <- length(sorted_values_desc)
  top_n <- as.integer(top_n)

  if (n_total == 0L) {
    stop("resolve_top_n_cutoff() received an empty vector.")
  }

  if (!is.finite(top_n) || is.na(top_n) || top_n < 1L) {
    stop("top_n must be a positive integer.")
  }

  # Use the empirically predetermined comparison-specific k* exactly. Do not
  # silently shrink k* or replace the rank rule with a loading-value cutoff.
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

  n_total <- length(sorted_empirical_p)

  if (n_total < 5L) {
    message(sprintf("[%s] Fewer than 5 empirical p-values were available for Higher Criticism.", dataset_name))
    return(NA_real_)
  }

  hc_scores <- suppressWarnings(
    tryCatch(
      fdrtool::hc.score(as.vector(sorted_empirical_p)),
      error = function(e) {
        message(
          sprintf(
            "[%s] hc.score failed: %s",
            dataset_name,
            conditionMessage(e)
          )
        )
        rep(NA_real_, n_total)
      }
    )
  )

  if (length(hc_scores) != n_total) {
    message(sprintf("[%s] Higher-Criticism score length mismatch.", dataset_name))
    return(NA_real_)
  }

  # Match fdrtool::hc.thresh(alpha0 = HC_ALPHA0): maximize HC only over the
  # lowest alpha0 fraction of ordered empirical p-values. The explicit score
  # calculation additionally permits the required positive-HC safeguard.
  # fdrtool::hc.thresh maximizes over 1:ceiling(alpha0 * n); use the same range.
  n_search <- max(1L, min(n_total, ceiling(HC_ALPHA0 * n_total)))
  search_idx <- seq_len(n_search)

  valid_idx <- search_idx[
    is.finite(hc_scores[search_idx]) &
      !is.na(hc_scores[search_idx]) &
      is.finite(sorted_empirical_p[search_idx]) &
      !is.na(sorted_empirical_p[search_idx]) &
      sorted_empirical_p[search_idx] > 0 &
      sorted_empirical_p[search_idx] < 1
  ]

  if (length(valid_idx) == 0L) {
    message(sprintf("[%s] No valid Higher-Criticism search points were available.", dataset_name))
    return(NA_real_)
  }

  best_idx <- valid_idx[which.max(hc_scores[valid_idx])]
  best_hc <- as.numeric(hc_scores[best_idx])
  hc_p <- as.numeric(sorted_empirical_p[best_idx])

  # Runtime consistency audit: when fdrtool::hc.thresh succeeds, our selected
  # lower-tail threshold must match hc.thresh(alpha0 = HC_ALPHA0). We calculate
  # scores explicitly only so that a non-positive maximum can be rejected.
  package_hc_p <- suppressWarnings(
    tryCatch(
      as.numeric(
        fdrtool::hc.thresh(
          as.vector(sorted_empirical_p),
          alpha0 = HC_ALPHA0,
          plot = FALSE
        )[1]
      ),
      error = function(e) NA_real_
    )
  )

  if (is.finite(package_hc_p) && !is.na(package_hc_p) &&
      !isTRUE(all.equal(hc_p, package_hc_p, tolerance = 1e-12))) {
    # Adopt the package threshold (the reference definition) and record it.
    pkg_idx <- which(sorted_empirical_p == package_hc_p)[1]
    warning(sprintf(
      "[%s] HC audit: explicit HCp %.17g differs from fdrtool::hc.thresh(alpha0=%.3f) %.17g; using the fdrtool value.",
      dataset_name, hc_p, HC_ALPHA0, package_hc_p), call. = FALSE)
    hc_p <- package_hc_p
    if (!is.na(pkg_idx)) best_hc <- as.numeric(hc_scores[pkg_idx])
  }

  # A non-positive maximum indicates no excess of small empirical p-values in
  # the prespecified HC search region. In that case no HC/HBFSS threshold is
  # asserted for the dataset/view.
  if (!is.finite(best_hc) || is.na(best_hc) || best_hc <= 0) {
    message(
      sprintf(
        "[%s] No positive Higher-Criticism signal within the lowest %.1f%% of empirical p-values; HBFSS disabled for this view.",
        dataset_name,
        100 * HC_ALPHA0
      )
    )
    return(NA_real_)
  }

  if (!is.finite(hc_p) || is.na(hc_p) || hc_p <= 0 || hc_p >= 1) {
    return(NA_real_)
  }

  message(
    sprintf(
      "[%s] HC calibration: alpha0=%.3f, HCmax=%.6f, HCp=%.8g, search=%d/%d empirical p-values.",
      dataset_name,
      HC_ALPHA0,
      best_hc,
      hc_p,
      n_search,
      n_total
    )
  )

  hc_p
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

# A panel's real data can legitimately have zero rows for a method (e.g. no
# Weak-CNH discoveries in a given view). When that happens to be true of the
# specific panel a legend gets borrowed from, ggplot silently omits that
# method's marker glyph from the legend key even with drop = FALSE. To
# guarantee every manuscript legend always shows all four method markers,
# build the shared legend from a small synthetic dataset that always
# contains exactly one row per method, rather than reusing a real panel.
build_full_method_legend <- function() {
  dummy <- data.frame(
    x = rep(0, length(significance_method_levels)),
    y = rep(0, length(significance_method_levels)),
    Method = factor(significance_method_levels, levels = significance_method_levels)
  )

  p <- ggplot(dummy, aes(x = x, y = y, color = Method, shape = Method, size = Method)) +
    geom_point(alpha = 0.98, stroke = 0.90) +
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
      drop = FALSE,
      guide = "none"
    ) +
    scale_size_manual(
      values = significance_method_sizes,
      breaks = significance_method_levels,
      guide = "none"
    ) +
    manuscript_theme() +
    theme(legend.position = "bottom")

  shared_panel_legend(p)
}



strip_legend <- function(p) {
  p + theme(legend.position = "none")
}

assemble_one_legend_panel <- function(plot_list, panel_title, ncol = length(plot_list), width_legend = TRUE) {
  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0L) {
    return(NULL)
  }

  legend <- build_full_method_legend()
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
  # The manuscript figure archive contains curated PNG panels only. Individual
  # per-view figures and duplicate PDF renditions remain available in the output
  # tree when generated, but are excluded from the archive to avoid redundancy.
  panel_patterns <- c(
    "/Panels/.*\\.png$",
    "/Combined_Figures/.*\\.png$",
    "/TWAS/figures/Figure_TWAS_Gene_Support\\.png$"
  )

  all_png <- list.files(
    output_dir,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.png$",
    ignore.case = TRUE
  )
  all_png <- all_png[file.info(all_png)$isdir %in% FALSE]

  normalized <- gsub("\\\\", "/", all_png)
  keep <- rep(FALSE, length(all_png))
  for (pat in panel_patterns) keep <- keep | grepl(pat, normalized)
  figure_files <- all_png[keep]

  # The global TWAS method-count figure is intentionally excluded because the
  # per-comparison five-view TWAS panels already contain the same count summary.
  figure_files <- figure_files[
    !grepl("Figure_TWAS_Method_Counts\\.png$", figure_files)
  ]

  if (!length(figure_files)) {
    stop("No curated manuscript PNG panels were found for SEQUENCE_ALL_FIGURES.zip.")
  }

  zip_path <- file.path(output_dir, "SEQUENCE_ALL_FIGURES.zip")
  if (file.exists(zip_path)) unlink(zip_path, force = TRUE)

  root_norm <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  file_norm <- normalizePath(figure_files, winslash = "/", mustWork = TRUE)
  relative_files <- substring(file_norm, nchar(root_norm) + 2L)
  relative_files <- sort(unique(relative_files))

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)
  utils::zip(zipfile = basename(zip_path), files = relative_files, flags = "-q")

  if (!file.exists(zip_path)) stop("Figure ZIP creation failed: ", zip_path)
  normalizePath(zip_path, winslash = "/", mustWork = TRUE)
}

create_all_tables_zip <- function() {
  # The manuscript table archive contains the concise tables used to interpret
  # significance, EVS/PCA evidence, and 3'aTWAS overlap. Per-view working tables
  # remain in their comparison folders but are not duplicated in the archive.
  wanted <- character(0)

  add_if_exists <- function(path) {
    if (file.exists(path)) wanted <<- c(wanted, path)
  }

  add_if_exists(file.path(summary_table_dir, "Table_DE_Method_Counts.csv"))
  for (tab in c("Table_EVS_Cutoffs.csv", "Table_Likelihood_Cutoff_Selected.csv", "Table_Likelihood_Cutoff_Scan.csv",
                "Table_Varying_Dispersion_Trends.csv", "Table_k_Stability.csv", "Table_Prior_Width_ML.csv",
                "Table_Dispersion_After_EVS.csv", "Table_Shrinkage_By_View.csv", "Table_Final_Dispersion_Score.csv",
                "Table_EVS_Composition.csv")) {
    add_if_exists(file.path(summary_table_dir, tab))
  }
  add_if_exists(file.path(summary_table_dir, "Table_EVS_Split_Audit.csv"))
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




# -----------------------------------------------------------------------------
# Gene/PAS annotation used by the final 3'aTWAS analysis
# -----------------------------------------------------------------------------

valid_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & nzchar(x) & x != "-" & grepl("[A-Za-z0-9]", x)
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

compute_pc1_loading_table <- function(value_df, sample_names, top_n, preprocessing_label = "Normalized prior to EVS") {
  x <- as.matrix(value_df[, sample_names, drop = FALSE])
  storage.mode(x) <- "numeric"

  if (ncol(x) < 2L) {
    stop("EVS PCA requires at least two samples in each condition.")
  }

  # Refined EVS method: run PCA directly on the pre-split matrix for the
  # condition, extract the PC1 feature eigenvector/loading, take absolute values,
  # rank descending, and select the comparison-specific empirical k*. No log
  # transformation is applied before PCA.
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

build_eigenvector_split <- function(comparison_name, count_matrix, coldata, track_key) {
  track_key <- match.arg(track_key, c("normalized_evs", "raw_evs"))
  empirical_k <- get_active_evs_cutoff(comparison_name)

  # The EVS matrix is used only to choose feature IDs. The normalized track
  # follows the supplied workflow: median-of-ratios normalization of the full
  # RT/ZT comparison matrix, followed by condition-specific PCA/PC1 eigenvectors.
  # RawEVS repeats the same PCA/ranking rule without pre-split normalization.
  # Both tracks use the same comparison-specific empirical k* so preprocessing
  # is the only difference in the selection rule. Downstream DESeq2 always
  # receives the corresponding raw-count subset.
  raw_count_matrix <- coerce_raw_count_matrix_for_deseq2(
    count_matrix,
    context = paste0(track_key, " full comparison matrix before EVS")
  )

  if (nrow(raw_count_matrix) < empirical_k) {
    stop(
      comparison_name, " ", track_key,
      ": empirical EVS k*=", empirical_k,
      " exceeds the ", nrow(raw_count_matrix),
      " PASs available after zero-row filtering."
    )
  }

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
    top_n = empirical_k,
    preprocessing_label = preprocessing_label
  )

  fit_untrt <- compute_pc1_loading_table(
    rank_matrix,
    untrt_ids,
    top_n = empirical_k,
    preprocessing_label = preprocessing_label
  )

  trt_high <- as.character(
    fit_trt$loading_table$feature_id[
      fit_trt$loading_table$rank <= empirical_k
    ]
  )

  untrt_high <- as.character(
    fit_untrt$loading_table$feature_id[
      fit_untrt$loading_table$rank <= empirical_k
    ]
  )

  if (length(trt_high) != empirical_k || length(untrt_high) != empirical_k) {
    stop(
      comparison_name, " ", track_key,
      ": EVS failed to select exactly empirical k*=", empirical_k,
      " PAS features per condition."
    )
  }

  # Joint/disjoint classification from the two independently ranked PC1 axes.
  joint_ids <- intersect(trt_high, untrt_high)
  disjoint_trt_ids <- setdiff(trt_high, untrt_high)
  disjoint_untrt_ids <- setdiff(untrt_high, trt_high)

  leading_edge_ids <- union(trt_high, untrt_high)
  remainder_ids <- setdiff(rownames(raw_count_matrix), leading_edge_ids)

  if (length(leading_edge_ids) == 0L) {
    stop("Leading-edge dataset is empty. Check sample mapping or EVS inputs.")
  }

  if (length(remainder_ids) == 0L) {
    stop(
      comparison_name, " ", track_key,
      ": remainder dataset is empty after empirical EVS k*=", empirical_k, "."
    )
  }

  if (length(intersect(leading_edge_ids, remainder_ids)) != 0L ||
      !setequal(union(leading_edge_ids, remainder_ids), rownames(raw_count_matrix))) {
    stop(comparison_name, " ", track_key, ": EVS Lead/Rem partition audit failed.")
  }

  list(
    comparison_name = comparison_name,
    track_key = track_key,
    empirical_evs_k = empirical_k,
    cutoff_basis = EMPIRICAL_EVS_CUTOFF_BASIS,
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
    hbfss_threshold = hbfss_threshold_dataset,
    fdrtool_param = fdr_fit$param
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

plot_final_volcano <- function(df, dataset_name, short_title = NULL, label_genes = TRUE, y_limit_override = NULL, n_labels = n_top_labels_volcano) {
  plot_df <- build_final_volcano_df(df, y_col = "neglog10_empirical_p")
  if (!nrow(plot_df)) stop("No finite volcano plotting rows for ", dataset_name)

  method_df <- build_significance_plot_long(plot_df)
  lab_df <- if (isTRUE(label_genes)) {
    select_final_volcano_labels(plot_df, y_col = "neglog10_empirical_p", n_labels = n_labels)
  } else {
    plot_df[0, , drop = FALSE]
  }

  hc_raw <- suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1]))
  htau <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))
  hc_y <- if (is.finite(hc_raw) && !is.na(hc_raw) && hc_raw > 0 && hc_raw <= 1) {
    safe_neglog10(hc_raw)
  } else NA_real_

  # A shared y_limit_override (computed once across every comparison/view)
  # keeps all manuscript volcano panels on the same vertical scale, so
  # significance magnitudes are visually comparable across timepoints
  # instead of each panel silently rescaling to its own local maximum.
  y_limit <- if (!is.null(y_limit_override) && is.finite(y_limit_override)) {
    y_limit_override
  } else {
    max(plot_df$neglog10_empirical_p, na.rm = TRUE) * 1.05
  }
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
    coord_cartesian(clip = "off", ylim = c(0, y_limit)) +
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
      # A finite max.overlaps lets ggrepel silently drop the labels it
      # genuinely cannot place without collision, rather than forcing every
      # requested label onto the page and producing illegible overlapping
      # text in dense regions (this matters most in the narrow 5-panel
      # manuscript figures, where each subplot has limited width).
      max.overlaps = 15,
      force = 2.2,
      force_pull = 0.25,
      box.padding = 0.45,
      point.padding = 0.18,
      min.segment.length = 0,
      segment.alpha = 0.60,
      segment.size = 0.22
    )
  }

  p
}

run_one_evs_track <- function(comparison_name, track_key, count_matrix, coldata, annot_df) {
  evs <- build_eigenvector_split(
    comparison_name = comparison_name,
    count_matrix = count_matrix,
    coldata = coldata,
    track_key = track_key
  )

  message(
    comparison_name, " ", unname(track_short[track_key]),
    ": empirical EVS k*=", evs$empirical_evs_k,
    " selected independently from each condition PC1 eigenvector; ",
    "Lead union n=", length(evs$leading_edge_ids),
    ", Rem n=", length(evs$remainder_ids), "."
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
      EVS_k_per_condition = if (identical(dataset_key, "raw_dataset")) NA_integer_ else get_active_evs_cutoff(comparison_name),
      PAS_tested = nrow(df),
      HC_alpha0 = HC_ALPHA0,
      HCp = hc,
      Htau = htau,
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
        EVS_k_per_condition = if (identical(dataset_key, "raw_dataset")) NA_integer_ else get_active_evs_cutoff(comparison_name),
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
        HC_alpha0 = HC_ALPHA0,
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
        EVS_k_per_condition = if (identical(dataset_key, "raw_dataset")) NA_integer_ else get_active_evs_cutoff(comparison_name),
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
    paste0(
      "This run used the ", CUTOFF_METHOD_LABEL, " EVS selection rule. ",
      "The active per-comparison cutoffs were: ",
      paste(
        paste0(comparison_table$comparison_name, " k=", vapply(comparison_table$comparison_name, get_active_evs_cutoff, integer(1))),
        collapse = ", "
      ),
      ". These cutoffs were locked before downstream differential-expression testing. ",
      "The same active cutoff was applied to NormEVS and RawEVS within each comparison so the two EVS tracks differed only in the matrix used for PC1 ranking, not in the number of PASs selected per condition. ",
      "NormEVS used DESeq2 median-of-ratios normalized counts before PCA, whereas RawEVS used raw counts before PCA. ",
      "Within each condition, prcomp was applied to the transposed feature-by-sample matrix with centering and without scaling. ",
      "PASs were ranked by absolute PC1 loading, and the comparison-specific k* highest-loading PASs were selected independently from the RT and ZT condition eigenvectors."
    ),
    "PASs present in both condition-specific top-k* sets were classified as Joint; PASs present in only one set were classified as Disjoint for that condition. Joint plus both Disjoint sets formed the Leading Edge, and all other PASs formed the Remainder. EVS defined PAS membership only. Downstream DESeq2 analyses received the corresponding raw-count subset and estimated normalization and dispersion parameters within that analysis view. The five reported views were Original (No EVS), NormEVS Lead, NormEVS Rem, RawEVS Lead, and RawEVS Rem.",
    "",
    "## Choice of the EVS cutoff k (likelihood-stratified moderation)",
    "Two cutoffs were analyzed: a prespecified k = 5,000 and a likelihood cutoff k*. For each comparison, k* was the k that maximized the marginal negative-binomial likelihood of the counts under DESeq2's empirical-Bayes dispersion model with the dispersion prior allowed to differ between the Leading Edge and the Remainder, a varying-dispersion regression in the sense of Tzougas (2020) fitted by EM. For PAS i in group g (Leading Edge or Remainder at cutoff k), y_ij | alpha_i ~ NB(mu_ij, alpha_i) and log alpha_i ~ N(f_g(log mean_i), tau_g^2), where f_g is a natural cubic spline in log mean (4 df) and tau_g^2 is the prior width. Each group was treated as its own dataset with its own median-of-ratios size factors, and mu_ij = s_j q_i,c(j) with q the per-condition normalized mean. The per-PAS likelihood of alpha was the Cox-Reid adjusted profile likelihood used by DESeq2 for gene-wise dispersions; alpha_i was integrated over a grid of 121 log-dispersion values from 1e-4 to 100. The E-step computed the posterior mean and second moment of log alpha_i; the M-step refitted f_g by least squares and set tau_g^2 to the mean posterior squared deviation, iterating until the log-likelihood changed by less than 0.001. The log-likelihood l(k) = l_LE(k) + l_REM(k) was compared with the Original (No EVS) model (one trend and one prior width for all PASs). All admissible k (both groups with at least 300 PASs) were scanned every 250, then every 10 within 250 of the best value, then at every single k within 10 of the best value; k* = argmax l(k), with support interval {k : l(k) >= l(k*) - 2}. The cutoff was fixed before any differential test and no p-value or discovery count entered its choice. The pooled prior width W = (n_LE tau_LE^2 + n_REM tau_REM^2) / (n_LE + n_REM) was compared with the Original prior width. Stability was assessed by refitting DESeq2 dispersions at k* - 100, -25, -5, -1, +1, +5, +25 and +100 and comparing final dispersions with those at k*.",
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
    sprintf(
      paste0(
        "Finite ordinary DESeq2 Wald statistics were supplied directly to fdrtool with statistic='normal' and cutoff.method='fndr' to estimate a zero-centered normal empirical null and the corresponding empirical-null p-values. ",
        "These empirical p-values were distinct from the ordinary DESeq2 Wald p-values and from the greaterAbs/lessAbs p-values. ",
        "Higher-Criticism scores were then calculated from the ordered empirical-null p-values with fdrtool::hc.score. ",
        "The HC maximization region was prespecified as the lowest %.0f%% of ordered empirical p-values (alpha0=%.2f), matching the lower-tail restriction provided by fdrtool::hc.thresh(alpha0=...). ",
        "Within that region, HCp was the empirical p-value at the maximum HC score. A dataset/view was assigned an HCp threshold only when that maximum HC score was strictly positive; if the maximum HC score was non-positive or no valid search point existed, HCp and Htau were left undefined and no PAS in that dataset/view was called HBFSS-significant."
      ),
      100 * HC_ALPHA0,
      HC_ALPHA0
    ),
    sprintf("For each PAS, HBFSS = |apeglm-shrunken LFC| x [-log10(empirical p)]. With c=%.1f, Htau = -log10(HCp) x c when a valid positive-HC threshold existed. A PAS was HBFSS-significant when HBFSS > Htau. HCp served only to derive Htau and was not imposed as an additional per-PAS significance gate. HBFSS significance did not use a DESeq2 adjusted-p-value gate.", lfc_boundary),
    "",
    "## Overlap",
    "Overlap was the number of unique HBFSS-significant PASs that also satisfied at least one DESeq2 criterion (Standard, Strong, or Weak-CNH). Overlap was a comparison quantity, not a separate significance test. Because final Weak required HBFSS, every final Weak PAS contributed to the overlap set.",
    "",
    "## PCA comparison after eigenvector splitting",
    "PCA comparison figures used the same matrices that defined EVS: full-comparison median-of-ratios normalized counts for NormEVS and raw counts for RawEVS. No additional log transformation or post-split re-normalization was applied. Original, Leading Edge, and Remainder were compared within the same preprocessing scale. For each view, the full PCA eigenspectrum was used to calculate the percentage of total variance explained by PC1 and PC2. PC1 eigenvalue retention was calculated as the first eigenvalue of the view divided by the first eigenvalue of the corresponding Original matrix. Original-PC1 loading energy captured by a feature subset was calculated as the sum of squared Original PC1 loadings for that subset divided by the sum of squared Original PC1 loadings for all PASs. Leading Edge and Remainder loading-energy fractions were required to sum to 100%. Treatment-control separation was summarized in the PC1-PC2 plane as centroid distance divided by pooled within-group root-mean-square distance. PCA panels reported PAS count, PC1 explained variance, PC1 eigenvalue retention, Original-PC1 loading-energy fraction, and the separation ratio.",
    "",
    "## Volcano figures and tables",
    sprintf("All significance volcanoes used the HBFSS plotting coordinates: apeglm-shrunken log2 fold change on the x-axis and -log10(empirical-null p) on the y-axis. Standard DESeq2, Strong greaterAbs, and Weak lessAbs significance were determined from their own DESeq2 p-values and BH-adjusted p-values, then their markers were projected onto these common HBFSS coordinates only for visual comparison. The HBFSS boundary y=Htau/|LFC| and HCp reference were drawn on the same axes. The same key was used throughout: Std = green diamond, Strong = red square, Weak = blue triangle, HBFSS = purple star. Multiple method markers were superimposed at the same PAS coordinate. LessAbs-only PASs remained background. Each volcano reported the number of significant PASs for Std, Strong, Weak, HBFSS, and their HBFSS/DESeq2 overlap and labeled at most %d final significant PASs.", n_top_labels_volcano),
    "Each comparison/view folder contained one significant-PAS table and one method-count table. Significant-PAS tables contained only the union of final Standard, Strong, Weak, and HBFSS discoveries and reported the DESeq2 p-values/BH-adjusted p-values used by the applicable DESeq2 tests together with empirical p, HBFSS, HC alpha0, HCp, Htau, apeglm LFC, and explicit method indicators. Method-count tables also reported HC alpha0, HCp, and Htau for view-level calibration auditing.",
    "",
    "## 3'aTWAS ortholog overlap",
    "Human 3'aTWAS gene symbols were mapped to rat gene symbols by combining database-supported babelgene human-to-rat ortholog mappings (top=FALSE) with direct case-insensitive symbol-equivalent matches present in the WTTS annotation. For every comparison and analysis view, mapped TWAS orthologs were intersected with Standard, Strong, Weak, and HBFSS WTTS discoveries. Concise TWAS tables reported the human TWAS symbol, rat ortholog, total TWAS records, distinct TWAS transcript count, TWAS multi-transcript/APA status, total WTTS PAS count, WTTS multi-PAS/APA status, significant WTTS PAS identifiers, significant multi-PAS status, and the method(s) and analysis view in which significance was observed.",
    "",
    "## Output organization",
    "Original (No EVS), NormEVS Lead, NormEVS Rem, RawEVS Lead, and RawEVS Rem were written to separate folders within each comparison for view-specific tables. Manuscript figures were organized as comparison-level panels so related volcano, PCA/EVS, and 3'aTWAS results could be reviewed side by side without redundant individual images. The curated figure archive SEQUENCE_ALL_FIGURES.zip contained manuscript PNG panels only. The separate archive SEQUENCE_ALL_TABLES.zip contained the locked EVS cutoff table, the likelihood cutoff scan and selection tables, the maximum-likelihood prior widths, the dispersion-after-EVS tables, the EVS split audit table, overall differential-expression method-count table, quantitative EVS/PCA evidence table, one all-view significant-PAS table and method-count table for each comparison, and the combined 3'aTWAS PAS-, gene-, and method-summary tables."
  )

  writeLines(methods_text, file.path(output_dir, "Methods_Manuscript.md"), useBytes = TRUE)
  invisible(TRUE)
}


# =============================================================================
# LIKELIHOOD-STRATIFIED MODERATION: CHOOSING k BY MAXIMUM LIKELIHOOD
# =============================================================================
# DESeq2 shrinks every gene-wise dispersion toward ONE mean-dispersion trend
# with ONE prior width. Eigenvector splitting (EVS) separates PASs into a
# Leading Edge (LE) and a Remainder (REM), and each is fitted as its own
# dataset, so each gets its own trend and its own prior. The cutoff k decides
# which PASs go where. k is chosen by maximum likelihood under DESeq2's own
# empirical-Bayes dispersion model, with the dispersion prior allowed to
# differ between the two groups (a varying-dispersion regression in the sense
# of Tzougas 2020, fitted by EM, with negative-binomial counts):
#
#   y_ij | alpha_i ~ NB(mu_ij, alpha_i)
#   log alpha_i    ~ N( f_g(log mubar_i), tau_g^2 ),     g = LE or REM at cutoff k
#
#   PAS i is in the LE at cutoff k when it is among the top k by |PC1 loading|
#   in either condition (NormEVS ranking), i.e. depth_i <= k.
#   f_g      natural cubic spline in log mean (LIK_SPLINE_DF df): group trend
#   tau_g^2  group prior width
#   mu_ij    s_j * q_i,cond(j), with the group's own median-of-ratios size
#            factors s_j (each group is its own dataset, as downstream)
#
# The per-PAS likelihood of alpha is the Cox-Reid adjusted profile likelihood
# that DESeq2 uses for its gene-wise estimates. alpha_i is integrated out on a
# fixed grid of log alpha, and (f_g, tau_g^2) are estimated by EM:
#   E-step  posterior of log alpha_i on the grid: E[log alpha_i], E[(log alpha_i)^2]
#   M-step  f_g = least-squares spline fit to E[log alpha_i]
#           tau_g^2 = mean of E[(log alpha_i - f_g(log mubar_i))^2]
#
# l(k) = l_LE(k) + l_REM(k) is the marginal log-likelihood of all counts. The
# Original (No EVS) model (one trend, one width, all PASs) gives l_0. Every k
# uses the same counts and the same number of parameters, so
#   k* = argmax_k l(k).
# Search, three stages: every LIK_COARSE_STEP (250) over all admissible k (both
# groups keep at least LIK_MIN_GROUP PASs); every LIK_FINE_STEP (10) within
# +/- LIK_FINE_HALF_WIDTH (250) of the best coarse k; then every single k
# within +/- LIK_EXACT_HALF_WIDTH (10) of the best k so far, so k* is exact.
# Support interval: the k with l(k) >= l(k*) - 2. No p-value or discovery
# count enters the choice.

LIK_LOG_ALPHA_GRID <- seq(log(1e-4), log(100), length.out = 121L)
LIK_SPLINE_DF <- 4L
LIK_MIN_GROUP <- 300L
LIK_COARSE_STEP <- 250L
LIK_FINE_STEP <- 10L
LIK_FINE_HALF_WIDTH <- 250L
LIK_EXACT_HALF_WIDTH <- 10L
LIK_EM_TOL <- 1e-3
LIK_EM_MAX_ITER <- 500L
LIK_STABILITY_STEPS <- c(-100L, -25L, -5L, -1L, 1L, 5L, 25L, 100L)
DISP_MIN_GENE_EST <- 1e-6

analysis_cores <- function() {
  requested <- suppressWarnings(as.integer(Sys.getenv("SEQUENCE_CORES", unset = "2")))
  available <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(requested) || requested < 1L) requested <- 1L
  if (is.na(available) || available < 1L) available <- 1L
  if (.Platform$OS.type == "windows") return(1L)
  as.integer(min(requested, available))
}

# Run fun over xs on `cores` forked workers, in chunks, with progress messages.
# Elements that fail return NULL and are dropped.
parallel_map <- function(xs, fun, cores, label) {
  chunks <- split(xs, ceiling(seq_along(xs) / max(1L, 8L * cores)))
  out <- list(); done <- 0L
  for (ch in chunks) {
    res <- if (cores > 1L) {
      parallel::mclapply(ch, fun, mc.cores = cores, mc.preschedule = FALSE)
    } else lapply(ch, fun)
    res <- res[!vapply(res, function(x) is.null(x) || inherits(x, "try-error"), logical(1))]
    out <- c(out, res)
    done <- done + length(ch)
    message(sprintf("    %s: %d / %d", label, done, length(xs)))
  }
  out
}

# depth_i = min(rank of |PC1 loading| among trt samples, rank among untrt
# samples) on the NormEVS matrix. The Leading Edge at cutoff k is exactly
# {i : depth_i <= k}, the same set build_eigenvector_split() selects.
evs_depth <- function(raw, coldata) {
  rank_obj <- make_rank_matrix_for_track(raw, coldata, "normalized_evs")
  ids <- colnames(raw)
  trt <- ids[coldata$condition == "trt"]
  untrt <- ids[coldata$condition == "untrt"]
  lt <- compute_pc1_loading_table(rank_obj$rank_matrix, trt, top_n = 1L)$loading_table
  lu <- compute_pc1_loading_table(rank_obj$rank_matrix, untrt, top_n = 1L)$loading_table
  d <- pmin(stats::setNames(lt$rank, lt$feature_id)[rownames(raw)],
            stats::setNames(lu$rank, lu$feature_id)[rownames(raw)])
  stats::setNames(as.integer(d), rownames(raw))
}

# One group = one dataset: its own size factors, per-condition means, and log mean.
nb_group_data <- function(raw_sub, condition) {
  condition <- factor(condition)
  sf <- DESeq2::estimateSizeFactorsForMatrix(raw_sub)
  nrm <- sweep(raw_sub, 2, sf, "/")
  q <- vapply(levels(condition), function(l) rowMeans(nrm[, condition == l, drop = FALSE]),
              numeric(nrow(raw_sub)))
  q <- matrix(q, nrow = nrow(raw_sub))
  mu <- sweep(q[, as.integer(condition), drop = FALSE], 2, sf, "*")
  y <- raw_sub
  storage.mode(y) <- "double"
  list(y = y, mu = mu, cond = as.integer(condition), log_mean = log(rowMeans(nrm)))
}

# Cox-Reid adjusted NB log-likelihood of every PAS at every grid value of log alpha.
nb_cr_loglik_grid <- function(gd, grid = LIK_LOG_ALPHA_GRID) {
  y <- gd$y; mu <- gd$mu; m <- ncol(y)
  ylogmu <- ifelse(y > 0, y * log(pmax(mu, 1e-300)), 0)
  const <- rowSums(ylogmu - lgamma(y + 1))
  cond_cols <- split(seq_len(m), gd$cond)
  L <- matrix(NA_real_, nrow(y), length(grid))
  for (q in seq_along(grid)) {
    a <- exp(grid[q]); r <- 1 / a
    ll <- rowSums(lgamma(y + r) - y * log(r + mu) - r * log1p(a * mu)) - m * lgamma(r) + const
    w <- mu / (1 + a * mu)
    cr <- 0
    for (cc in cond_cols) cr <- cr + log(pmax(rowSums(w[, cc, drop = FALSE]), 1e-300))
    L[, q] <- ll - 0.5 * cr
  }
  L
}

# EM for the log-normal dispersion prior with a spline trend in log mean.
em_dispersion_prior <- function(L, log_mean, grid = LIK_LOG_ALPHA_GRID, df = LIK_SPLINE_DF,
                                tol = LIK_EM_TOL, max_iter = LIK_EM_MAX_ITER) {
  n <- nrow(L)
  spline_basis <- splines::ns(log_mean, df = df)
  qrB <- qr(cbind(1, spline_basis))
  dgrid <- grid[2] - grid[1]
  theta0 <- grid[max.col(L, ties.method = "first")]
  centre <- qr.fitted(qrB, theta0)
  tau2 <- max(stats::var(theta0 - centre) - 0.3, 0.05)
  ll_old <- -Inf
  it <- 0L
  repeat {
    it <- it + 1L
    lp <- L - outer(-centre, grid, "+")^2 / (2 * tau2)
    mx <- lp[cbind(seq_len(n), max.col(lp, ties.method = "first"))]
    w <- exp(lp - mx)
    s <- rowSums(w)
    ll <- sum(mx + log(s)) + n * (log(dgrid) - 0.5 * log(2 * pi * tau2))
    w <- w / s
    e1 <- drop(w %*% grid)
    e2 <- drop(w %*% (grid^2))
    if (ll - ll_old < tol || it >= max_iter) break
    ll_old <- ll
    centre <- qr.fitted(qrB, e1)
    tau2 <- max(mean(e2 - 2 * e1 * centre + centre^2), 1e-4)
  }
  list(loglik = ll, tau2 = tau2, coef = qr.coef(qrB, centre), spline = spline_basis,
       log_mean = log_mean, posterior_mean = e1, iterations = it, converged = it < max_iter)
}

fit_dispersion_group <- function(raw_sub, condition) {
  gd <- nb_group_data(raw_sub, condition)
  fit <- em_dispersion_prior(nb_cr_loglik_grid(gd), gd$log_mean)
  fit$n <- nrow(raw_sub)
  fit
}

# Fitted trend exp(f_g) on a grid of means, for plotting.
dispersion_trend_curve <- function(fit, n_points = 100L) {
  rng <- stats::quantile(fit$log_mean, c(0.005, 0.995), names = FALSE)
  x <- seq(rng[1], rng[2], length.out = n_points)
  B <- cbind(1, stats::predict(fit$spline, x))
  data.frame(baseMean = exp(x), trend = exp(drop(B %*% fit$coef)))
}

fit_split_model <- function(raw, condition, depth, k, keep_fits = FALSE) {
  le <- depth <= k
  if (sum(le) < LIK_MIN_GROUP || sum(!le) < LIK_MIN_GROUP) return(NULL)
  a <- fit_dispersion_group(raw[le, , drop = FALSE], condition)
  b <- fit_dispersion_group(raw[!le, , drop = FALSE], condition)
  row <- data.frame(k = as.integer(k), LE_n = a$n, REM_n = b$n, loglik = a$loglik + b$loglik,
                    tau2_LE = a$tau2, tau2_REM = b$tau2,
                    em_iterations_LE = a$iterations, em_iterations_REM = b$iterations,
                    em_converged = a$converged && b$converged)
  if (keep_fits) list(row = row, LE = a, REM = b) else list(row = row)
}

admissible_k_range <- function(depth) {
  sd <- unname(sort(depth)); N <- length(sd)
  if (N < 2L * LIK_MIN_GROUP) stop("Too few PASs for a Leading Edge / Remainder split.")
  c(lo = sd[LIK_MIN_GROUP], hi = sd[N - LIK_MIN_GROUP + 1L] - 1L)
}

# Full scan for one comparison. Returns the scan table, the selected k*,
# fitted trends at k*, and the EVS depth used.
likelihood_k_scan <- function(raw, coldata, label, cores = analysis_cores()) {
  raw <- raw[rowSums(raw) > 0, , drop = FALSE]
  cond <- as.character(coldata$condition)
  depth <- evs_depth(raw, coldata)
  orig <- fit_dispersion_group(raw, cond)
  rng <- admissible_k_range(depth)
  eval_k <- function(k) tryCatch(fit_split_model(raw, cond, depth, k)$row, error = function(e) {
    message("    ", label, " k = ", k, " skipped: ", conditionMessage(e)); NULL })

  coarse_k <- seq(LIK_COARSE_STEP, rng[["hi"]], by = LIK_COARSE_STEP)
  coarse_k <- coarse_k[coarse_k >= rng[["lo"]]]
  message(sprintf("  %s: likelihood scan over k = %s-%s (%d coarse cutoffs, %d core%s)",
                  label, format(min(coarse_k), big.mark = ","), format(max(coarse_k), big.mark = ","),
                  length(coarse_k), cores, if (cores > 1L) "s" else ""))
  rows <- do.call(rbind, parallel_map(coarse_k, eval_k, cores, paste(label, "coarse scan")))
  if (is.null(rows) || !nrow(rows)) stop(label, ": likelihood scan produced no admissible cutoffs.")
  kb <- rows$k[which.max(rows$loglik)]
  fine_k <- seq(max(rng[["lo"]], kb - LIK_FINE_HALF_WIDTH), min(rng[["hi"]], kb + LIK_FINE_HALF_WIDTH),
                by = LIK_FINE_STEP)
  fine_k <- setdiff(fine_k, rows$k)
  if (length(fine_k)) {
    rows <- rbind(rows, do.call(rbind, parallel_map(fine_k, eval_k, cores, paste(label, "fine scan"))))
  }
  kb <- rows$k[which.max(rows$loglik)]
  exact_k <- setdiff(seq(max(rng[["lo"]], kb - LIK_EXACT_HALF_WIDTH), min(rng[["hi"]], kb + LIK_EXACT_HALF_WIDTH)),
                     rows$k)
  if (length(exact_k)) {
    rows <- rbind(rows, do.call(rbind, parallel_map(exact_k, eval_k, cores, paste(label, "exact scan"))))
  }
  rows <- rows[order(rows$k), ]
  rows$loglik_Original <- orig$loglik
  rows$delta_loglik_vs_Original <- rows$loglik - orig$loglik
  rows$tau2_Original <- orig$tau2
  rows$W_split <- (rows$LE_n * rows$tau2_LE + rows$REM_n * rows$tau2_REM) / (rows$LE_n + rows$REM_n)

  best <- rows[which.max(rows$loglik), ]
  near <- rows[abs(rows$k - best$k) <= LIK_FINE_HALF_WIDTH, ]
  support <- near$k[near$loglik >= best$loglik - 2]
  at5k <- rows[rows$k == 5000L, ]
  n_obs <- nrow(raw) * ncol(raw)
  p_extra <- LIK_SPLINE_DF + 2L  # one more spline trend (df + intercept) and one more prior width
  selected <- data.frame(
    k_star = best$k, support_lo = min(support), support_hi = max(support),
    LE_n = best$LE_n, REM_n = best$REM_n,
    loglik_k_star = best$loglik, loglik_Original = orig$loglik,
    delta_loglik_vs_Original = best$delta_loglik_vs_Original,
    delta_BIC_vs_Original = -2 * best$delta_loglik_vs_Original + p_extra * log(n_obs),
    delta_loglik_k_star_vs_5000 = if (nrow(at5k)) best$loglik - at5k$loglik else NA_real_,
    tau2_Original = orig$tau2, tau2_LE = best$tau2_LE, tau2_REM = best$tau2_REM,
    W_split = best$W_split, prior_width_reduction_pct = 100 * (1 - best$W_split / orig$tau2),
    n_cutoffs_evaluated = nrow(rows))

  final <- fit_split_model(raw, cond, depth, best$k, keep_fits = TRUE)
  trends <- rbind(
    cbind(Set = ORIGINAL_LABEL, dispersion_trend_curve(orig)),
    cbind(Set = "Leading Edge", dispersion_trend_curve(final$LE)),
    cbind(Set = "Remainder", dispersion_trend_curve(final$REM)))
  message(sprintf("  %s: k* = %s (support %s-%s); log-likelihood gain over Original = %s; vs k = 5,000 = %s",
                  label, format(best$k, big.mark = ","), format(min(support), big.mark = ","),
                  format(max(support), big.mark = ","), format(round(best$delta_loglik_vs_Original), big.mark = ","),
                  if (nrow(at5k)) format(round(best$loglik - at5k$loglik, 1)) else "NA"))
  list(scan = rows, selected = selected, trends = trends, depth = depth, k_star = as.integer(best$k))
}

# Prior widths (ML, same model) at a given k, used for the Fixed 5,000 run.
ml_prior_widths_at_k <- function(raw, coldata, k) {
  raw <- raw[rowSums(raw) > 0, , drop = FALSE]
  cond <- as.character(coldata$condition)
  depth <- evs_depth(raw, coldata)
  orig <- fit_dispersion_group(raw, cond)
  sp <- fit_split_model(raw, cond, depth, k)$row
  if (is.null(sp)) stop("k = ", k, " leaves fewer than ", LIK_MIN_GROUP, " PASs in a group.")
  W <- (sp$LE_n * sp$tau2_LE + sp$REM_n * sp$tau2_REM) / (sp$LE_n + sp$REM_n)
  data.frame(k = as.integer(k), LE_n = sp$LE_n, REM_n = sp$REM_n, tau2_Original = orig$tau2,
             tau2_LE = sp$tau2_LE, tau2_REM = sp$tau2_REM, W_split = W,
             prior_width_reduction_pct = 100 * (1 - W / orig$tau2),
             delta_loglik_vs_Original = sp$loglik - orig$loglik)
}

# Stability of DESeq2's final dispersions when k moves away from k*.
deseq2_final_dispersions <- function(raw_sub, coldata) {
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_sub, colData = coldata,
                                        design = make_design_formula(coldata))
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- suppressMessages(DESeq2::estimateDispersions(dds, quiet = TRUE))
  stats::setNames(DESeq2::dispersions(dds), rownames(dds))
}

k_stability_test <- function(raw, coldata, depth, k_star, label, cores = analysis_cores()) {
  raw <- raw[names(depth), , drop = FALSE]
  rng <- admissible_k_range(depth)
  ks <- k_star + c(0L, LIK_STABILITY_STEPS)
  ks <- ks[ks >= rng[["lo"]] & ks <= rng[["hi"]]]
  message(sprintf("  %s: stability of final dispersions at %d cutoffs around k*", label, length(ks)))
  fit_k <- function(k) {
    le <- depth <= k
    tryCatch(c(deseq2_final_dispersions(raw[le, , drop = FALSE], coldata),
               deseq2_final_dispersions(raw[!le, , drop = FALSE], coldata)),
             error = function(e) NULL)
  }
  fits <- parallel_map(as.list(ks), fit_k, cores, paste(label, "stability refits"))
  if (length(fits) != length(ks)) stop(label, ": a stability refit failed.")
  ref <- fits[[1]]
  in_le_ref <- depth[names(ref)] <= k_star
  do.call(rbind, lapply(seq_along(ks)[-1], function(i) {
    d <- fits[[i]][names(ref)]
    lr <- abs(log2(d / ref))
    switched <- (depth[names(ref)] <= ks[i]) != in_le_ref
    ok <- is.finite(lr)
    data.frame(delta_k = ks[i] - k_star, k = ks[i], n_switched = sum(switched),
               median_abs_log2_change_switched = if (any(switched & ok)) stats::median(lr[switched & ok]) else NA_real_,
               median_abs_log2_change_stayed = stats::median(lr[!switched & ok]),
               p95_abs_log2_change_stayed = stats::quantile(lr[!switched & ok], 0.95, names = FALSE))
  }))
}

# ---- per-run state -----------------------------------------------------------

LIKELIHOOD_K <- list()
likelihood_scan_rows <- list()
likelihood_selected_rows <- list()
likelihood_trend_rows <- list()
stability_rows <- list()
prior_width_rows <- list()
dispersion_metrics_rows <- list()
dispersion_trend_rows <- list()
dispersion_pas_rows <- list()
composition_rows <- list()

# Likelihood run: choose k* for this comparison (before any split or test).
run_likelihood_cutoff <- function(comparison_name, count_matrix, coldata) {
  raw <- coerce_raw_count_matrix_for_deseq2(count_matrix, context = paste0(comparison_name, " likelihood scan"))
  res <- likelihood_k_scan(raw, coldata, comparison_name)
  LIKELIHOOD_K[[comparison_name]] <<- res$k_star
  likelihood_scan_rows[[comparison_name]] <<- cbind(Comparison = comparison_name, res$scan)
  likelihood_selected_rows[[comparison_name]] <<- cbind(Comparison = comparison_name, res$selected)
  likelihood_trend_rows[[comparison_name]] <<- cbind(Comparison = comparison_name, res$trends)
  prior_width_rows[[comparison_name]] <<- data.frame(
    Comparison = comparison_name, Cutoff_Method = CUTOFF_METHOD_LABEL,
    res$selected[, c("LE_n", "REM_n", "tau2_Original", "tau2_LE", "tau2_REM", "W_split",
                     "prior_width_reduction_pct", "delta_loglik_vs_Original")],
    k = res$k_star)
  stab <- tryCatch(k_stability_test(raw[rowSums(raw) > 0, , drop = FALSE], coldata, res$depth,
                                    res$k_star, comparison_name),
                   error = function(e) { warning("Stability test skipped for ", comparison_name, ": ",
                                                 conditionMessage(e), call. = FALSE); NULL })
  if (!is.null(stab)) stability_rows[[comparison_name]] <<- cbind(Comparison = comparison_name, stab)
  invisible(res$k_star)
}

# Fixed run: prior widths at the prespecified k, from the same model.
run_fixed_prior_widths <- function(comparison_name, count_matrix, coldata) {
  raw <- coerce_raw_count_matrix_for_deseq2(count_matrix, context = paste0(comparison_name, " prior widths"))
  k <- get_active_evs_cutoff(comparison_name)
  pw <- ml_prior_widths_at_k(raw, coldata, k)
  prior_width_rows[[comparison_name]] <<- data.frame(
    Comparison = comparison_name, Cutoff_Method = CUTOFF_METHOD_LABEL,
    pw[, c("LE_n", "REM_n", "tau2_Original", "tau2_LE", "tau2_REM", "W_split",
           "prior_width_reduction_pct", "delta_loglik_vs_Original")], k = k)
  invisible(TRUE)
}

# ---- DESeq2 dispersions after the split (read from the fitted objects) ------

VIEW_SET <- c("Original (No EVS)" = "Original (No EVS)", "NormEVS Lead" = "Leading Edge",
              "NormEVS Rem" = "Remainder")

`%||%` <- function(a, b) if (is.null(a)) b else a

dispersion_view_metrics <- function(dds) {
  mc <- S4Vectors::mcols(dds)
  fin <- DESeq2::dispersions(dds)
  ok <- !mc$allZero & is.finite(mc$dispGeneEst) & mc$dispGeneEst >= DISP_MIN_GENE_EST &
        is.finite(mc$dispFit) & mc$dispFit > 0 & is.finite(fin) & mc$baseMean > 0
  dfun <- DESeq2::dispersionFunction(dds)
  X <- stats::model.matrix(DESeq2::design(dds), as.data.frame(SummarizedExperiment::colData(dds)))
  v <- attr(dfun, "varLogDispEsts")
  outl <- if (is.null(mc$dispOutlier)) rep(FALSE, length(ok)) else (mc$dispOutlier %in% TRUE)
  data.frame(
    n_PAS = sum(ok),
    var_log_disp_residuals = if (length(v) == 1L) as.numeric(v) else NA_real_,
    sampling_var_log_disp = trigamma((ncol(dds) - ncol(X)) / 2),
    deseq2_disp_prior_var = as.numeric(attr(dfun, "dispPriorVar") %||% NA_real_),
    fit_type = as.character(attr(dfun, "fitType") %||% NA_character_),
    dispersion_outlier_pct = 100 * mean(outl[ok]),
    median_final_dispersion = stats::median(fin[ok]))
}

collect_dispersion_after_evs <- function(comparison_name) {
  views <- comparison_analysis_views()
  norm_reg <- paper_registry[[registry_key(comparison_name, "normalized_evs")]]
  for (i in seq_len(nrow(views))) {
    obj <- paper_registry[[registry_key(comparison_name, views$track_key[i])]]
    dds <- obj$analysis_results[[views$dataset_key[i]]]$dds
    if (is.null(dds)) next
    label <- analysis_view_label(views$track_key[i], views$dataset_key[i])
    m <- tryCatch(dispersion_view_metrics(dds), error = function(e) NULL)
    if (!is.null(m)) {
      dispersion_metrics_rows[[length(dispersion_metrics_rows) + 1L]] <<-
        cbind(data.frame(Comparison = comparison_name, View = label), m)
    }
    if (!label %in% names(VIEW_SET)) next
    mc <- S4Vectors::mcols(dds)
    keep <- !mc$allZero & mc$baseMean > 0
    outl <- if (is.null(mc$dispOutlier)) rep(FALSE, nrow(dds)) else (mc$dispOutlier %in% TRUE)
    dispersion_pas_rows[[length(dispersion_pas_rows) + 1L]] <<- data.frame(
      Comparison = comparison_name, View = label, feature_id = rownames(dds)[keep],
      baseMean = mc$baseMean[keep], dispGeneEst = mc$dispGeneEst[keep], dispFit = mc$dispFit[keep],
      dispersion = DESeq2::dispersions(dds)[keep], dispOutlier = outl[keep])
    bm <- mc$baseMean[keep]
    rng <- stats::quantile(bm, c(0.005, 0.995), names = FALSE)
    grid_bm <- 10^seq(log10(rng[1]), log10(rng[2]), length.out = 100)
    dispersion_trend_rows[[length(dispersion_trend_rows) + 1L]] <<- data.frame(
      Comparison = comparison_name, View = label, baseMean = grid_bm,
      trend = DESeq2::dispersionFunction(dds)(grid_bm))
  }
  evs <- norm_reg$evs
  if (!is.null(evs)) {
    composition_rows[[comparison_name]] <<- data.frame(
      Comparison = comparison_name, k = evs$empirical_evs_k,
      Class = c("Joint (top-k in both)", "Disjoint (top-k in one)"),
      n = c(length(evs$joint_ids), length(evs$disjoint_trt_ids) + length(evs$disjoint_untrt_ids)))
  }
  invisible(TRUE)
}

# Final (MAP) dispersion of each PAS in the NormEVS Leading Edge / Remainder
# against the same PAS in the Original fit. The Wald SE ratio is approximated
# at the PAS's mean as sqrt[(1/mu + alpha_split) / (1/mu + alpha_Original)].
final_dispersion_pairs <- function(pas) {
  o <- pas[pas$View == ORIGINAL_LABEL, c("Comparison", "feature_id", "baseMean", "dispersion")]
  names(o)[3:4] <- c("baseMean_Original", "disp_Original")
  s <- pas[pas$View %in% c("NormEVS Lead", "NormEVS Rem"), c("Comparison", "View", "feature_id", "dispersion")]
  names(s)[4] <- "disp_split"
  d <- merge(s, o, by = c("Comparison", "feature_id"))
  d <- d[is.finite(d$disp_split) & d$disp_split > 0 & is.finite(d$disp_Original) & d$disp_Original > 0 &
         is.finite(d$baseMean_Original) & d$baseMean_Original > 0, ]
  d$Set <- unname(VIEW_SET[d$View])
  d$log2_ratio <- log2(d$disp_split / d$disp_Original)
  inv_mu <- 1 / d$baseMean_Original
  d$log_SE_ratio <- 0.5 * log((inv_mu + d$disp_split) / (inv_mu + d$disp_Original))
  d[, c("Comparison", "Set", "feature_id", "baseMean_Original", "disp_Original", "disp_split",
        "log2_ratio", "log_SE_ratio")]
}

final_dispersion_score <- function(pairs) {
  pct <- function(x) 100 * (exp(mean(x)) - 1)
  do.call(rbind, lapply(split(pairs, pairs$Comparison), function(d) {
    le <- d$Set == "Leading Edge"
    data.frame(Comparison = d$Comparison[1], Cutoff_Method = CUTOFF_METHOD_LABEL,
               LE_n = sum(le), REM_n = sum(!le),
               SE_change_LE_pct = pct(d$log_SE_ratio[le]), SE_change_REM_pct = pct(d$log_SE_ratio[!le]),
               median_log2_final_dispersion_ratio_LE = stats::median(d$log2_ratio[le]),
               median_log2_final_dispersion_ratio_REM = stats::median(d$log2_ratio[!le]))
  }))
}

# How far DESeq2 moves each gene-wise estimate toward its trend:
# pull = log(gene-wise / final) / log(gene-wise / trend), 0 = not moved, 1 = on the trend.
# Dispersion outliers (> 2 prior SDs above the trend) keep their gene-wise estimate.
shrinkage_by_view <- function(pas) {
  d <- pas[pas$View %in% names(VIEW_SET) & is.finite(pas$dispGeneEst) & pas$dispGeneEst >= DISP_MIN_GENE_EST &
           is.finite(pas$dispFit) & pas$dispFit > 0 & is.finite(pas$dispersion) & pas$dispersion > 0, ]
  do.call(rbind, lapply(split(d, list(d$Comparison, d$View), drop = TRUE), function(x) {
    gap <- log(x$dispGeneEst / x$dispFit)
    use <- !x$dispOutlier & abs(gap) > 0.1
    pull <- pmin(pmax(log(x$dispGeneEst[use] / x$dispersion[use]) / gap[use], 0), 1)
    data.frame(Comparison = x$Comparison[1], Cutoff_Method = CUTOFF_METHOD_LABEL, View = x$View[1],
               Set = unname(VIEW_SET[x$View[1]]), n_PAS = nrow(x),
               mean_pull_toward_trend = mean(pull), outlier_pct = 100 * mean(x$dispOutlier))
  }))
}

# ---- figure helpers ----------------------------------------------------------

comp_levels <- function() as.character(comparison_table$comparison_name)
log_lab <- function() scales::trans_format("log10", scales::math_format(10^.x))
k_lab <- function(x) ifelse(is.na(x), NA_character_, ifelse(x == 0, "0", paste0(x / 1000, "k")))
as_comp <- function(x) factor(x, levels = comp_levels())

# Stratified subsample across the mean range so plots show every abundance.
sample_across_mean <- function(d, n, bins = 20L, seed = 20261031L) {
  if (nrow(d) <= n) return(d)
  set.seed(seed)
  b <- cut(log(d$baseMean), bins, labels = FALSE)
  per <- max(1L, ceiling(n / bins))
  idx <- unlist(lapply(split(seq_len(nrow(d)), b), function(i) if (length(i) > per) sample(i, per) else i))
  d[idx, , drop = FALSE]
}

# ---- Figure 1: choosing k ------------------------------------------------------

panel_loglik <- function(scan, sel, tag) {
  scan$Comparison <- as_comp(scan$Comparison); sel$Comparison <- as_comp(sel$Comparison)
  pk <- data.frame(Comparison = sel$Comparison, k = sel$k_star, y = sel$delta_loglik_vs_Original,
                   lab = paste0("k* = ", format(sel$k_star, big.mark = ",")))
  p5 <- scan[scan$k == 5000L, c("Comparison", "k", "delta_loglik_vs_Original")]
  ggplot(scan, aes(k, delta_loglik_vs_Original)) +
    geom_area(fill = METHOD_COLS[["Likelihood k*"]], alpha = 0.07) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.3, colour = "grey45") +
    geom_vline(xintercept = 5000, linetype = "22", linewidth = 0.3, colour = METHOD_COLS[["Fixed 5,000"]]) +
    geom_line(linewidth = 0.55, colour = "grey15") +
    geom_point(data = p5, shape = 15, size = 1.9, colour = METHOD_COLS[["Fixed 5,000"]]) +
    geom_point(data = pk, aes(k, y), shape = 21, size = 2.8, stroke = 0.7, colour = "white",
               fill = METHOD_COLS[["Likelihood k*"]]) +
    geom_label(data = pk, aes(k, y, label = lab), hjust = -0.15, vjust = 1.3, size = 2.2,
               fontface = "bold", colour = METHOD_COLS[["Likelihood k*"]], label.size = 0,
               label.padding = unit(0.6, "mm"), fill = "white") +
    facet_wrap(~ Comparison, nrow = 1, scales = "free_y", labeller = labeller(Comparison = comp_lab())) +
    scale_x_continuous(labels = k_lab, breaks = scales::breaks_pretty(4)) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "k (top-|PC1| PASs per condition)", y = "Log-likelihood gain\nover Original (No EVS)", tag = tag,
         title = "Choosing k by maximum likelihood",
         subtitle = "Negative-binomial likelihood of the counts under the varying-dispersion model. Green: k*; grey square and dashed line: Fixed 5,000.") +
    fig_theme()
}

panel_width_vs_k <- function(scan, sel, tag) {
  long <- rbind(data.frame(Comparison = scan$Comparison, k = scan$k, w = scan$tau2_LE, Line = "Leading Edge"),
                data.frame(Comparison = scan$Comparison, k = scan$k, w = scan$tau2_REM, Line = "Remainder"),
                data.frame(Comparison = scan$Comparison, k = scan$k, w = scan$W_split, Line = "Pooled W"))
  long$Line <- factor(long$Line, levels = c("Leading Edge", "Remainder", "Pooled W"))
  long$Comparison <- as_comp(long$Comparison)
  o <- unique(scan[, c("Comparison", "tau2_Original")]); o$Comparison <- as_comp(o$Comparison)
  s <- sel; s$Comparison <- as_comp(s$Comparison)
  ggplot(long, aes(k, w, colour = Line, linetype = Line)) +
    geom_hline(data = o, aes(yintercept = tau2_Original), linetype = "11", linewidth = 0.45, colour = "grey15") +
    geom_vline(data = s, aes(xintercept = k_star), linewidth = 0.45, colour = METHOD_COLS[["Likelihood k*"]]) +
    geom_vline(xintercept = 5000, linetype = "22", linewidth = 0.3, colour = METHOD_COLS[["Fixed 5,000"]]) +
    geom_line(linewidth = 0.65) +
    facet_wrap(~ Comparison, nrow = 1, labeller = labeller(Comparison = comp_lab())) +
    scale_colour_manual(values = c("Leading Edge" = unname(SET_COLS["Leading Edge"]),
                                   "Remainder" = unname(SET_COLS["Remainder"]), "Pooled W" = "grey15")) +
    scale_linetype_manual(values = c("Leading Edge" = "solid", "Remainder" = "solid", "Pooled W" = "42")) +
    scale_x_continuous(labels = k_lab, breaks = scales::breaks_pretty(4)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(x = "k (top-|PC1| PASs per condition)", y = "Prior width \u03c4\u00b2", tag = tag,
         title = "One prior width for the Leading Edge, one for the Remainder",
         subtitle = "Maximum-likelihood prior widths as k changes. Dotted: Original (No EVS). Green: k*. Grey dashed: Fixed 5,000.") +
    fig_theme()
}

panel_varying_dispersion_fit <- function(pas, trends, sel, tag, max_pts = 2500L) {
  d <- pas[pas$View %in% c("NormEVS Lead", "NormEVS Rem") & is.finite(pas$dispGeneEst) &
           pas$dispGeneEst >= DISP_MIN_GENE_EST, c("Comparison", "View", "baseMean", "dispGeneEst")]
  d <- do.call(rbind, lapply(split(d, list(d$Comparison, d$View), drop = TRUE), sample_across_mean, n = max_pts))
  d$Set <- factor(unname(VIEW_SET[d$View]), levels = c("Remainder", "Leading Edge"))
  d <- d[order(d$Set), ]
  d$Comparison <- as_comp(d$Comparison)
  t <- trends; t$Comparison <- as_comp(t$Comparison)
  t$Set <- factor(t$Set, levels = c("Leading Edge", "Remainder", ORIGINAL_LABEL))
  lab <- sel; lab$Comparison <- as_comp(lab$Comparison)
  lab$txt <- sprintf("Gain over one trend:\n%s log-lik units\nΔBIC = %s",
                     format(round(lab$delta_loglik_vs_Original), big.mark = ","),
                     format(round(lab$delta_BIC_vs_Original), big.mark = ","))
  yl <- stats::quantile(d$dispGeneEst, c(0.002, 0.998), names = FALSE)
  ggplot(d, aes(baseMean, dispGeneEst)) +
    geom_point(data = d[d$Set == "Remainder", ], colour = "#E8836F", size = 0.35, alpha = 0.35, stroke = 0) +
    geom_point(data = d[d$Set == "Leading Edge", ], colour = "#5B9BD5", size = 0.35, alpha = 0.35, stroke = 0) +
    geom_line(data = t, aes(baseMean, trend, group = Set), colour = "white", linewidth = 1.5) +
    geom_line(data = t, aes(baseMean, trend, linetype = Set, colour = Set), linewidth = 0.75) +
    geom_label(data = lab, aes(x = Inf, y = Inf, label = txt), inherit.aes = FALSE, hjust = 1.04, vjust = 1.08,
               size = 1.9, lineheight = 0.95, label.size = 0.15, label.padding = unit(0.8, "mm"), fill = "white") +
    facet_wrap(~ Comparison, nrow = 1, labeller = labeller(Comparison = comp_lab())) +
    scale_colour_manual(values = c("Leading Edge" = unname(SET_TREND_COLS["Leading Edge"]),
                                   "Remainder" = unname(SET_TREND_COLS["Remainder"]),
                                   "Original (No EVS)" = unname(SET_TREND_COLS["Original (No EVS)"])),
                        breaks = c("Leading Edge", "Remainder", ORIGINAL_LABEL),
                        labels = c("Leading Edge trend", "Remainder trend", "Single trend, Original (No EVS)")) +
    scale_linetype_manual(values = c("Leading Edge" = "solid", "Remainder" = "solid", "Original (No EVS)" = "11"),
                          breaks = c("Leading Edge", "Remainder", ORIGINAL_LABEL),
                          labels = c("Leading Edge trend", "Remainder trend", "Single trend, Original (No EVS)")) +
    scale_x_log10(labels = log_lab()) + scale_y_log10(labels = log_lab()) +
    coord_cartesian(ylim = yl) +
    labs(x = "Mean of normalized counts", y = "Gene-wise dispersion", tag = tag,
         title = "The fitted varying-dispersion model at k*",
         subtitle = "log \u03b1 = f_g(log mean) + u,  u ~ N(0, \u03c4\u00b2_g), one f_g and \u03c4\u00b2_g per group. Points: DESeq2 gene-wise estimates (blue LE, red REM). Original (No EVS): one trend for all PASs.") +
    fig_theme()
}

# ---- Figure 2: what the split does to DESeq2 ---------------------------------

panel_shrinkage_arrows <- function(pas, trd, met, tag, n_per_group = 160L, n_outliers = 30L) {
  d <- pas[pas$View %in% c("NormEVS Lead", "NormEVS Rem") & is.finite(pas$dispGeneEst) &
           pas$dispGeneEst >= DISP_MIN_GENE_EST & is.finite(pas$dispersion) & pas$dispersion > 0, ]
  pick <- function(x) {
    out <- x[x$dispOutlier, , drop = FALSE]
    if (nrow(out) > n_outliers) { set.seed(7L); out <- out[sample.int(nrow(out), n_outliers), ] }
    rbind(sample_across_mean(x[!x$dispOutlier, , drop = FALSE], n_per_group), out)
  }
  d <- do.call(rbind, lapply(split(d, list(d$Comparison, d$View), drop = TRUE), pick))
  d$Set <- unname(VIEW_SET[d$View]); d$Comparison <- as_comp(d$Comparison)
  pulled <- d[!d$dispOutlier, ]; outl <- d[d$dispOutlier, ]
  t <- trd[trd$View %in% c("NormEVS Lead", "NormEVS Rem"), ]
  t$Set <- unname(VIEW_SET[t$View]); t$Comparison <- as_comp(t$Comparison)
  # DESeq2 outlier rule: log(gene-wise) > log(trend) + 2 * sqrt(varLogDispEsts), per group
  v <- met[met$View %in% c("NormEVS Lead", "NormEVS Rem"), c("Comparison", "View", "var_log_disp_residuals")]
  t$threshold <- t$trend * exp(2 * sqrt(v$var_log_disp_residuals[match(paste(t$Comparison, t$View),
                                                                        paste(v$Comparison, v$View))]))
  ggplot() +
    geom_line(data = t, aes(baseMean, trend, colour = Set), linewidth = 0.7, show.legend = FALSE) +
    geom_line(data = t[is.finite(t$threshold), ], aes(baseMean, threshold, colour = Set), linewidth = 0.4,
              linetype = "22", show.legend = FALSE) +
    geom_segment(data = pulled, aes(x = baseMean, xend = baseMean, y = dispGeneEst, yend = dispersion, colour = Set),
                 linewidth = 0.25, alpha = 0.75, arrow = arrow(length = unit(0.6, "mm"), type = "closed")) +
    geom_point(data = outl, aes(baseMean, dispGeneEst, colour = Set,
                                shape = "Outlier: above its group's dashed 2-SD line, kept at its gene-wise estimate"),
               size = 1.3, stroke = 0.5, fill = NA) +
    facet_wrap(~ Comparison, nrow = 1, labeller = labeller(Comparison = comp_lab())) +
    scale_colour_manual(values = c("Leading Edge" = unname(SET_COLS["Leading Edge"]),
                                   "Remainder" = unname(SET_COLS["Remainder"])),
                        breaks = c("Leading Edge", "Remainder"),
                        labels = c("Leading Edge: pulled to its own trend", "Remainder: pulled to its own trend")) +
    scale_shape_manual(values = 21) +
    scale_x_log10(labels = log_lab()) + scale_y_log10(labels = log_lab()) +
    labs(x = "Mean of normalized counts", y = "Dispersion", tag = tag,
         title = "Which PASs are pulled toward their trend, and which are not",
         subtitle = "Arrows: DESeq2 gene-wise estimate to final (MAP) dispersion, pulled to the group's own trend (solid). Open circles: outliers above the group's 2-SD line (dashed), not shrunk.") +
    guides(shape = guide_legend(order = 2, override.aes = list(colour = "grey20")), colour = guide_legend(order = 1)) +
    fig_theme()
}

set_bars_scale <- function() {
  scale_fill_manual(values = SET_COLS, breaks = names(SET_COLS))
}

panel_shrinkage <- function(shr, tag) {
  d <- shr; d$Comparison <- as_comp(d$Comparison); d$Set <- factor(d$Set, levels = names(SET_COLS))
  ggplot(d, aes(Comparison, mean_pull_toward_trend, fill = Set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.74) +
    geom_text(aes(label = sprintf("%.1f%%", outlier_pct), y = mean_pull_toward_trend + 0.015),
              position = position_dodge(width = 0.8), angle = 90, hjust = 0, size = 1.8, colour = "grey30") +
    set_bars_scale() +
    scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.22))) +
    labs(x = NULL, y = "Mean pull toward trend\n(1 = fully shrunk)", tag = tag,
         title = "How much each group is shrunk",
         subtitle = "Above bars: % dispersion outliers left unshrunk.") +
    fig_theme() + theme(panel.grid.major.x = element_blank(), legend.position = "none")
}

panel_prior_width <- function(pw, tag) {
  m <- rbind(data.frame(Comparison = pw$Comparison, Set = ORIGINAL_LABEL, w = pw$tau2_Original),
             data.frame(Comparison = pw$Comparison, Set = "Leading Edge", w = pw$tau2_LE),
             data.frame(Comparison = pw$Comparison, Set = "Remainder", w = pw$tau2_REM))
  m$Set <- factor(m$Set, levels = names(SET_COLS)); m$Comparison <- as_comp(m$Comparison)
  s <- pw; s$Comparison <- as_comp(s$Comparison); s$x <- as.numeric(s$Comparison)
  s$top <- pmax(s$tau2_Original, s$tau2_LE, s$tau2_REM)
  s$lab <- sprintf("W %.2f\n%s%.0f%%", s$W_split, ifelse(s$prior_width_reduction_pct >= 0, "−", "+"),
                   abs(s$prior_width_reduction_pct))
  ggplot(m, aes(Comparison, w, fill = Set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.74) +
    geom_segment(data = s, aes(x = x - 0.42, xend = x + 0.42, y = W_split, yend = W_split),
                 inherit.aes = FALSE, colour = METHOD_COLS[["Likelihood k*"]], linewidth = 1.1) +
    geom_text(data = s, aes(x = x, y = top, label = lab), inherit.aes = FALSE, vjust = -0.3, size = 2.0,
              lineheight = 0.9, fontface = "bold", colour = METHOD_COLS[["Likelihood k*"]]) +
    set_bars_scale() +
    scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.3))) +
    labs(x = NULL, y = "Prior width \u03c4\u00b2", tag = tag,
         title = "The split narrows the prior",
         subtitle = "Green line: pooled width W; % narrower than Original.") +
    fig_theme() + theme(panel.grid.major.x = element_blank())
}

panel_composition <- function(comp, tag) {
  d <- comp; d$Class <- factor(d$Class, levels = rev(names(COMPOSITION_COLS)))
  tot <- stats::aggregate(n ~ Comparison + k, data = d, FUN = sum)
  lev <- rev(comp_levels())
  ylab <- stats::setNames(sprintf("%s\nk = %s", sub("_", " / ", tot$Comparison), format(tot$k, big.mark = ",")),
                          tot$Comparison)
  d$Comparison <- factor(d$Comparison, levels = lev); tot$Comparison <- factor(tot$Comparison, levels = lev)
  ggplot(d, aes(n, Comparison, fill = Class)) +
    geom_col(width = 0.66) +
    geom_text(data = tot, aes(n, Comparison, label = format(n, big.mark = ",")), inherit.aes = FALSE,
              hjust = -0.12, size = 2.0, fontface = "bold") +
    scale_fill_manual(values = COMPOSITION_COLS, breaks = names(COMPOSITION_COLS)) +
    guides(fill = guide_legend(nrow = 2)) +
    scale_y_discrete(labels = function(x) unname(ylab[x])) +
    scale_x_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.2))) +
    labs(x = "PASs in the Leading Edge", y = NULL, tag = tag,
         title = "What the Leading Edge is made of", subtitle = "Joint vs disjoint PASs.") +
    fig_theme() + theme(panel.grid.major.y = element_blank())
}

panel_final_change <- function(pairs, fscore, tag) {
  d <- pairs; d$Comparison <- as_comp(d$Comparison)
  d$Set <- factor(d$Set, levels = c("Leading Edge", "Remainder"))
  yl <- stats::quantile(d$log2_ratio, c(0.005, 0.995), na.rm = TRUE, names = FALSE)
  f <- rbind(data.frame(Comparison = fscore$Comparison, Set = "Leading Edge", v = fscore$SE_change_LE_pct),
             data.frame(Comparison = fscore$Comparison, Set = "Remainder", v = fscore$SE_change_REM_pct))
  f$Comparison <- as_comp(f$Comparison); f$Set <- factor(f$Set, levels = levels(d$Set))
  f$lab <- sprintf("SE %+.1f%%", f$v)
  ggplot(d, aes(Comparison, log2_ratio, fill = Set)) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.3, colour = "grey45") +
    geom_violin(position = position_dodge(width = 0.85), width = 0.8, linewidth = 0, scale = "width", alpha = 0.85) +
    geom_boxplot(aes(group = interaction(Comparison, Set)), position = position_dodge(width = 0.85),
                 width = 0.1, outlier.shape = NA, linewidth = 0.3, colour = "grey10", fill = "white") +
    geom_text(data = f, aes(Comparison, y = yl[2] + 0.22 * diff(yl), label = lab, colour = Set),
              position = position_dodge(width = 0.85), size = 2.1, fontface = "bold", show.legend = FALSE) +
    coord_cartesian(ylim = c(yl[1], yl[2] + 0.3 * diff(yl))) +
    scale_fill_manual(values = SET_COLS[c("Leading Edge", "Remainder")]) +
    scale_colour_manual(values = SET_COLS[c("Leading Edge", "Remainder")]) +
    scale_x_discrete(labels = comp_lab()) +
    labs(x = NULL, y = "log2 final dispersion\n(split / Original)", tag = tag,
         title = "Consequence: each PAS's final dispersion moves to its own population's level",
         subtitle = "Leading Edge rises (no longer shrunk toward too low a trend); Remainder falls (sharper tests). Labels: mean change in Wald standard error.") +
    fig_theme() + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(face = "bold"))
}

# ---- Supplementary S1: precision and stability of k* ------------------------

panel_zoom <- function(scan, sel, tag) {
  s <- sel; s$Comparison <- as_comp(s$Comparison)
  z <- do.call(rbind, lapply(split(scan, scan$Comparison), function(x) {
    ks <- s$k_star[s$Comparison == x$Comparison[1]]
    x <- x[abs(x$k - ks) <= LIK_FINE_HALF_WIDTH, ]
    x$rel <- x$loglik - max(x$loglik)
    x
  }))
  z$Comparison <- as_comp(z$Comparison)
  s$lab <- sprintf("k* = %s\nsupport %s–%s", format(s$k_star, big.mark = ","),
                   format(s$support_lo, big.mark = ","), format(s$support_hi, big.mark = ","))
  ggplot(z, aes(k, rel)) +
    geom_rect(data = s, aes(xmin = support_lo - LIK_FINE_STEP / 2, xmax = support_hi + LIK_FINE_STEP / 2,
                            ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
              fill = METHOD_COLS[["Likelihood k*"]], alpha = 0.15) +
    geom_hline(yintercept = -2, linetype = "22", linewidth = 0.3, colour = "grey45") +
    geom_line(linewidth = 0.4, colour = "grey20") + geom_point(size = 0.5, colour = "grey20") +
    geom_point(data = s, aes(k_star, 0), inherit.aes = FALSE, shape = 21, size = 2.6, stroke = 0.7,
               colour = "white", fill = METHOD_COLS[["Likelihood k*"]]) +
    geom_text(data = s, aes(x = Inf, y = -Inf, label = lab), inherit.aes = FALSE, hjust = 1.05, vjust = -0.3,
              size = 2.0, lineheight = 0.95) +
    facet_wrap(~ Comparison, nrow = 1, scales = "free", labeller = labeller(Comparison = comp_lab())) +
    scale_x_continuous(labels = scales::label_comma(), breaks = scales::breaks_pretty(3)) +
    labs(x = "k", y = "Log-likelihood − maximum", tag = tag,
         title = "Zoom on the peak: how sharply the data pin down k",
         subtitle = sprintf("Every %d within ±%d of k*. Shaded: support interval (within 2 log-likelihood units of the maximum).",
                            LIK_FINE_STEP, LIK_FINE_HALF_WIDTH)) +
    fig_theme()
}

panel_stability <- function(stab, tag) {
  lev <- sprintf("%+d", sort(unique(stab$delta_k)))
  long <- rbind(data.frame(Comparison = stab$Comparison, dk = sprintf("%+d", stab$delta_k),
                           v = stab$median_abs_log2_change_switched, n = stab$n_switched,
                           Group = "PASs that switched group"),
                data.frame(Comparison = stab$Comparison, dk = sprintf("%+d", stab$delta_k),
                           v = stab$median_abs_log2_change_stayed, n = NA, Group = "All other PASs"))
  long <- long[is.finite(long$v) & long$v > 0, ]
  long$dk <- factor(long$dk, levels = lev); long$Comparison <- as_comp(long$Comparison)
  long$Group <- factor(long$Group, levels = c("PASs that switched group", "All other PASs"))
  sw <- long[long$Group == "PASs that switched group", ]
  ggplot(long, aes(dk, v, colour = Group, group = Group, shape = Group)) +
    geom_line(linewidth = 0.5) + geom_point(size = 1.5) +
    geom_text(data = sw, aes(label = n), vjust = -0.9, size = 1.8, show.legend = FALSE) +
    facet_wrap(~ Comparison, nrow = 1, labeller = labeller(Comparison = comp_lab())) +
    scale_y_log10(labels = log_lab(), expand = expansion(mult = c(0.08, 0.25))) +
    scale_colour_manual(values = c("PASs that switched group" = "#6A3D9A", "All other PASs" = "#33A02C")) +
    scale_shape_manual(values = c("PASs that switched group" = 16, "All other PASs" = 15)) +
    labs(x = "Change in k from k*", y = "Median |log2 change| in\nfinal (MAP) dispersion", tag = tag,
         title = "Stability: does moving k by a few points change the final dispersions?",
         subtitle = "Each k refit end to end with DESeq2 (size factors, trends, prior, MAP). Numbers: PASs that crossed the boundary.") +
    fig_theme() + theme(axis.text.x = element_text(size = 5.5))
}

# ---- outputs ------------------------------------------------------------------

save_likelihood_and_dispersion_outputs <- function() {
  out <- summary_table_dir
  pas <- if (length(dispersion_pas_rows)) dplyr::bind_rows(dispersion_pas_rows) else NULL
  trd <- if (length(dispersion_trend_rows)) dplyr::bind_rows(dispersion_trend_rows) else NULL
  met <- if (length(dispersion_metrics_rows)) dplyr::bind_rows(dispersion_metrics_rows) else NULL
  comp <- if (length(composition_rows)) dplyr::bind_rows(composition_rows) else NULL
  pw <- if (length(prior_width_rows)) dplyr::bind_rows(prior_width_rows) else NULL
  if (!is.null(met)) save_csv(met, file.path(out, "Table_Dispersion_After_EVS.csv"))
  if (!is.null(pw)) save_csv(pw, file.path(out, "Table_Prior_Width_ML.csv"))
  if (!is.null(comp)) save_csv(comp, file.path(out, "Table_EVS_Composition.csv"))
  if (!is.null(trd)) save_csv(trd, file.path(out, "Table_Dispersion_Trends_By_View.csv"))
  if (!is.null(pas)) save_csv(pas, file.path(out, "Table_Dispersion_Estimates_By_View.csv"))
  pairs <- fscore <- shr <- NULL
  if (!is.null(pas)) {
    pairs <- final_dispersion_pairs(pas); fscore <- final_dispersion_score(pairs); shr <- shrinkage_by_view(pas)
    save_csv(pairs, file.path(out, "Table_Final_Dispersion_Change_By_PAS.csv"))
    save_csv(fscore, file.path(out, "Table_Final_Dispersion_Score.csv"))
    save_csv(shr, file.path(out, "Table_Shrinkage_By_View.csv"))
  }
  scan <- sel <- trends <- stab <- NULL
  if (length(likelihood_scan_rows)) {
    scan <- dplyr::bind_rows(likelihood_scan_rows); sel <- dplyr::bind_rows(likelihood_selected_rows)
    trends <- dplyr::bind_rows(likelihood_trend_rows)
    save_csv(scan, file.path(out, "Table_Likelihood_Cutoff_Scan.csv"))
    save_csv(sel, file.path(out, "Table_Likelihood_Cutoff_Selected.csv"))
    save_csv(trends, file.path(out, "Table_Varying_Dispersion_Trends.csv"))
  }
  if (length(stability_rows)) {
    stab <- dplyr::bind_rows(stability_rows)
    save_csv(stab, file.path(out, "Table_k_Stability.csv"))
  }

  if (!is.null(scan) && !is.null(pas)) try_fig({
    fig1 <- assemble_rows(list(panel_loglik(scan, sel, "a"), panel_width_vs_k(scan, sel, "b"),
                               panel_varying_dispersion_fit(pas, trends, sel, "c")),
                          heights = c(1, 1, 1.1),
                          title = "Figure 1. Choosing the eigenvector split by maximum likelihood")
    save_figure(fig1, paper_fig_dir, "Figure_1_Choosing_k", 180, 190)
  }, "Figure 1")

  if (!is.null(pas) && !is.null(pw) && !is.null(comp) && !is.null(met)) try_fig({
    fig2 <- assemble_rows(
      list(panel_shrinkage_arrows(pas, trd, met, "a"),
           list(plots = list(panel_shrinkage(shr, "b"), panel_prior_width(pw, "c"), panel_composition(comp, "d")),
                widths = c(1, 1, 1)),
           panel_final_change(pairs, fscore, "e")),
      heights = c(1.35, 1.1, 0.95),
      title = paste0("Figure 2. What likelihood-stratified moderation does to DESeq2's dispersion estimates (",
                     CUTOFF_METHOD_LABEL, ")"))
    save_figure(fig2, paper_fig_dir, paste0("Figure_2_Dispersion_After_EVS_", CUTOFF_METHOD_SLUG), 180, 215)
  }, "Figure 2")

  if (!is.null(scan) && !is.null(stab)) try_fig({
    figS1 <- assemble_rows(list(panel_zoom(scan, sel, "a"), panel_stability(stab, "b")), heights = c(1, 1),
                           title = "Supplementary Figure S1. Precision and stability of k*")
    save_figure(figS1, paper_fig_dir, "Figure_S1_k_Stability", 180, 120)
  }, "Figure S1")
  invisible(TRUE)
}

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  dir.create(file.path(output_dir, comparison_name), recursive = TRUE, showWarnings = FALSE)

  # The cutoff is fixed before any split or test: Likelihood k* is estimated
  # from this comparison's counts (no p-values); Fixed 5,000 is prespecified.
  if (identical(CUTOFF_METHOD_KEY, "likelihood")) {
    run_likelihood_cutoff(comparison_name, count_matrix, coldata)
  } else {
    tryCatch(run_fixed_prior_widths(comparison_name, count_matrix, coldata),
             error = function(e) warning("Prior widths skipped for ", comparison_name, ": ",
                                         conditionMessage(e), call. = FALSE))
  }

  for (track_key in c("normalized_evs", "raw_evs")) {
    message("Running ", comparison_name, " ", unname(track_short[track_key]))
    ok <- run_one_evs_track(comparison_name = comparison_name, track_key = track_key,
                            count_matrix = count_matrix, coldata = coldata, annot_df = annot_df)
    if (!isTRUE(ok)) stop("Analysis track failed: ", comparison_name, " / ", track_key)
  }

  tryCatch(collect_dispersion_after_evs(comparison_name),
           error = function(e) warning("Dispersion summaries skipped for ", comparison_name, ": ",
                                       conditionMessage(e), call. = FALSE))
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


compute_global_volcano_y_limit <- function() {
  views <- comparison_analysis_views()
  running_max <- NA_real_

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (i in seq_len(nrow(views))) {
      df <- read_result_table_for_panel(
        comparison_name = comparison_name,
        track_key = views$track_key[i],
        dataset_key = views$dataset_key[i]
      )
      if (is.null(df)) next

      plot_df <- tryCatch(
        build_final_volcano_df(df, y_col = "neglog10_empirical_p"),
        error = function(e) NULL
      )
      if (is.null(plot_df) || !nrow(plot_df)) next

      this_max <- suppressWarnings(max(plot_df$neglog10_empirical_p, na.rm = TRUE))
      if (is.finite(this_max)) {
        running_max <- if (is.na(running_max)) this_max else max(running_max, this_max)
      }
    }
  }

  if (is.na(running_max)) return(NULL)
  running_max * 1.05
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

  # Computed once so every comparison/view volcano in the manuscript shares
  # the same y-axis scale (see plot_final_volcano's y_limit_override).
  shared_y_limit <- compute_global_volcano_y_limit()

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
        label_genes = TRUE,
        y_limit_override = shared_y_limit,
        n_labels = 10L
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
      panel_title = paste0(comparison_name, " | ", CUTOFF_METHOD_LABEL, " | significance across ", length(plots), " analysis view", ifelse(length(plots) == 1L, "", "s")),
      ncol = length(plots)
    )

    panel_dir <- file.path(output_dir, comparison_name, "Panels")
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
    base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_Volcano_", length(plots), "Views"))
    panel_width <- max(10.0, 5.0 * length(plots))
    save_grob(panel, paste0(base, ".png"), width = panel_width, height = 7.6)
    save_grob(panel, paste0(base, ".pdf"), width = panel_width, height = 7.6)
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

build_evs_split_audit_table <- function() {
  rows <- list()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (track_key in c("normalized_evs", "raw_evs")) {
      obj <- paper_registry[[registry_key(comparison_name, track_key)]]
      if (is.null(obj) || is.null(obj$evs)) next

      evs <- obj$evs
      total_n <- nrow(evs$raw_dataset)
      lead_n <- length(evs$leading_edge_ids)
      rem_n <- length(evs$remainder_ids)

      if (length(evs$trt_top_ids) != evs$empirical_evs_k ||
          length(evs$untrt_top_ids) != evs$empirical_evs_k ||
          lead_n + rem_n != total_n) {
        stop(comparison_name, " ", unname(track_short[track_key]),
             ": empirical EVS split audit failed.")
      }

      rows[[length(rows) + 1L]] <- data.frame(
        Comparison = comparison_name,
        Track = unname(track_short[track_key]),
        Empirical_EVS_k_per_condition = evs$empirical_evs_k,
        TRT_selected_n = length(evs$trt_top_ids),
        UNTRT_selected_n = length(evs$untrt_top_ids),
        Joint_n = length(evs$joint_ids),
        TRT_disjoint_n = length(evs$disjoint_trt_ids),
        UNTRT_disjoint_n = length(evs$disjoint_untrt_ids),
        Leading_edge_union_n = lead_n,
        Remainder_n = rem_n,
        Total_nonzero_PAS_n = total_n,
        Leading_edge_pct = 100 * lead_n / total_n,
        Remainder_pct = 100 * rem_n / total_n,
        TRT_PC1_loading_cutoff = as.numeric(evs$fit_trt$cutoff),
        UNTRT_PC1_loading_cutoff = as.numeric(evs$fit_untrt$cutoff),
        Cutoff_basis = evs$cutoff_basis,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) data.frame() else dplyr::bind_rows(rows)
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
        Empirical_EVS_k_per_condition = obj$evs$empirical_evs_k,
        Leading_edge_union_n = length(obj$evs$leading_edge_ids),
        Remainder_n = length(obj$evs$remainder_ids),
        TRT_PC1_loading_cutoff = as.numeric(obj$evs$fit_trt$cutoff),
        UNTRT_PC1_loading_cutoff = as.numeric(obj$evs$fit_untrt$cutoff),
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

  tbl$Dataset <- factor(tbl$Dataset, levels = c("Orig", "Lead", "Rem"))

  # Three directly interpretable quantities are shown together: PC1 explained
  # variance, PC1 eigenvalue retained relative to Orig, and the fraction of the
  # Original PC1 squared-loading energy contained in each feature subset.
  long <- tbl %>%
    dplyr::select(
      Comparison, Track, Dataset,
      PC1_explained_pct,
      PC1_eigenvalue_retained_pct,
      Original_PC1_loading_energy_pct
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        PC1_explained_pct,
        PC1_eigenvalue_retained_pct,
        Original_PC1_loading_energy_pct
      ),
      names_to = "Metric",
      values_to = "Percent"
    ) %>%
    dplyr::mutate(
      Metric = factor(
        Metric,
        levels = c(
          "PC1_explained_pct",
          "PC1_eigenvalue_retained_pct",
          "Original_PC1_loading_energy_pct"
        ),
        labels = c("PC1 explained", "λ1 vs Orig", "Orig PC1 energy")
      )
    )

  ggplot(long, aes(Dataset, Percent, group = Metric, shape = Metric)) +
    geom_hline(yintercept = 100, linewidth = 0.3, linetype = "dashed", colour = "grey70") +
    geom_line(aes(linetype = Metric), linewidth = 0.55, position = position_dodge(width = 0.08)) +
    geom_point(size = 2.6, position = position_dodge(width = 0.08)) +
    geom_text(
      aes(label = paste0(round(Percent, 1), "%")),
      position = position_dodge(width = 0.08),
      vjust = -0.7,
      size = 2.5,
      check_overlap = TRUE
    ) +
    facet_wrap(~ Track, nrow = 1) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    labs(
      title = paste0(comparison_name, " | quantitative EVS/PC1 evidence"),
      x = NULL,
      y = "Percent",
      shape = NULL,
      linetype = NULL
    ) +
    manuscript_theme() +
    theme(legend.position = "bottom")
}

plot_empirical_hbfss_support <- function(df, short_title) {
  req <- c("empirical_p", "HBFSS", "hc_p_threshold_dataset", "hbfss_threshold_dataset")
  if (length(setdiff(req, names(df))) > 0L) return(NULL)

  plot_df <- df[
    is.finite(df$empirical_p) & !is.na(df$empirical_p) &
      is.finite(df$HBFSS) & !is.na(df$HBFSS),
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(NULL)

  plot_df$empirical_p <- pmax(plot_df$empirical_p, 1e-300)
  method_df <- build_significance_plot_long(plot_df)

  hc_raw <- suppressWarnings(as.numeric(plot_df$hc_p_threshold_dataset[1]))
  htau <- suppressWarnings(as.numeric(plot_df$hbfss_threshold_dataset[1]))

  p <- ggplot() +
    geom_point(
      data = plot_df,
      aes(empirical_p, HBFSS),
      color = plot_palette$background,
      size = 0.45,
      alpha = 0.20
    ) +
    geom_point(
      data = method_df,
      aes(empirical_p, HBFSS, color = Method, shape = Method),
      size = 2.25,
      alpha = 0.98,
      stroke = 0.80
    ) +
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
    scale_x_log10(labels = scales::label_scientific()) +
    labs(
      title = short_title,
      x = "Empirical p",
      y = "HBFSS",
      caption = paste0(
        if (is.finite(hc_raw) && !is.na(hc_raw)) paste0("HCp=", signif(hc_raw, 3)) else "HCp=NA",
        "  ",
        if (is.finite(htau) && !is.na(htau)) paste0("Hτ=", signif(htau, 3)) else "Hτ=NA"
      )
    ) +
    manuscript_theme() +
    theme(legend.position = "bottom")

  if (is.finite(hc_raw) && !is.na(hc_raw) && hc_raw > 0 && hc_raw <= 1) {
    p <- p + geom_vline(xintercept = hc_raw, color = plot_palette$hc, linewidth = 0.55, linetype = "dotted")
  }
  if (is.finite(htau) && !is.na(htau)) {
    p <- p + geom_hline(yintercept = htau, color = plot_palette$hbfss_line, linewidth = 0.55, linetype = "dashed")
  }
  p
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

    if (length(plots)) {
      panel <- assemble_one_legend_panel(
        plots,
        panel_title = paste0(comparison_name, " | ", CUTOFF_METHOD_LABEL, " | PCA structure before and after EVS"),
        ncol = 3
      )
      base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_PCA_EVS_2x3"))
      save_grob(panel, paste0(base, ".png"), width = 16.5, height = 10.8)
      save_grob(panel, paste0(base, ".pdf"), width = 16.5, height = 10.8)
    }

    evidence_tbl <- build_evs_pca_evidence_table(comparison_name)
    if (nrow(evidence_tbl)) {
      all_evidence[[comparison_name]] <- evidence_tbl
      save_csv(
        evidence_tbl,
        file.path(output_dir, comparison_name, paste0("Table_", comparison_name, "_EVS_PCA_Evidence.csv"))
      )
    }

    pevidence <- plot_evs_pca_evidence(comparison_name)
    if (!is.null(pevidence)) {
      base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_EVS_PC1_Evidence"))
      save_grob(pevidence, paste0(base, ".png"), width = 11.0, height = 5.2)
      save_grob(pevidence, paste0(base, ".pdf"), width = 11.0, height = 5.2)
    }
  }

  evidence_all <- if (length(all_evidence)) dplyr::bind_rows(all_evidence) else data.frame()
  if (nrow(evidence_all)) {
    save_csv(evidence_all, file.path(summary_table_dir, "Table_EVS_PCA_Evidence.csv"))
  }

  invisible(evidence_all)
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
      panel_title <- paste0("HBFSS calibration | Orig | ", CUTOFF_METHOD_LABEL)
      output_suffix <- "Raw_AllComparisons"
      panel_height <- 5.8
    } else {
      track_order_use <- c("normalized_evs", "raw_evs")
      panel_title <- paste0("HBFSS calibration | ", unname(dataset_short[dataset_key]), " (Norm vs Raw) | ", CUTOFF_METHOD_LABEL)
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
    facet_wrap(~ Analysis, scales = "free_y", nrow = 1) +
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
      title = paste0("Significant PAS counts by method | ", CUTOFF_METHOD_LABEL),
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

  save_grob(p, file.path(paper_fig_dir, "Figure_Manuscript_Discovery_Counts.png"), width = 21.0, height = 6.5)
  save_grob(p, file.path(paper_fig_dir, "Figure_Manuscript_Discovery_Counts.pdf"), width = 21.0, height = 6.5)
  invisible(TRUE)
}

save_paper_support_figures <- function(summary_df) {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) {
    return(invisible(FALSE))
  }

  support_steps <- list(
    PCA = function() save_paper_pca_panels(),
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


# =============================================================================
# FINAL 3'aTWAS ORTHOLOG COMPARISON
# =============================================================================

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
  twas <- utils::read.csv(
    path,
    skip = 1,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  assert_required_columns(
    twas,
    c("Disease", "PANEL", "Transcript ID", "Gene symbol", "3'aTWAS.Z", "3'aTWAS.P"),
    object_name = "3'aTWAS file"
  )

  twas %>%
    dplyr::transmute(
      Disease = trimws(as.character(Disease)),
      Panel = trimws(as.character(PANEL)),
      TWAS_Transcript = trimws(as.character(`Transcript ID`)),
      Human_TWAS_Gene = trimws(as.character(`Gene symbol`)),
      TWAS_Z = suppressWarnings(as.numeric(`3'aTWAS.Z`)),
      TWAS_P = suppressWarnings(as.numeric(`3'aTWAS.P`)),
      COLOC_PP4 = if ("COLOC.PP4" %in% names(twas)) suppressWarnings(as.numeric(COLOC.PP4)) else NA_real_
    ) %>%
    dplyr::filter(valid_gene_symbol(Human_TWAS_Gene))
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
  mapped <- twas %>%
    dplyr::left_join(mapping, by = "Human_TWAS_Gene") %>%
    dplyr::filter(!is.na(gene_key))

  mapped %>%
    dplyr::group_by(gene_key) %>%
    dplyr::summarise(
      Rat_Ortholog = dplyr::first(Rat_Ortholog[valid_gene_symbol(Rat_Ortholog)]),
      Human_TWAS_Genes = collapse_unique(Human_TWAS_Gene),
      TWAS_Diseases = collapse_unique(Disease),
      TWAS_n = dplyr::n(),
      TWAS_Panel_n = dplyr::n_distinct(Panel),
      TWAS_Transcript_n = dplyr::n_distinct(TWAS_Transcript),
      TWAS_APA = dplyr::n_distinct(TWAS_Transcript) >= 2L,
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
        HC_alpha0 = HC_ALPHA0,
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
    dplyr::inner_join(twas_rat_meta, by = "gene_key") %>%
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
      Analysis = factor(Analysis, levels = c(
        "Original (No EVS)",
        "NormEVS Lead", "NormEVS Rem",
        "RawEVS Lead", "RawEVS Rem"
      ))
    ) %>%
    dplyr::transmute(
      Comparison,
      Analysis = as.character(Analysis),
      Human_TWAS = Human_TWAS_Genes,
      Rat_Ortholog,
      TWAS_n,
      TWAS_Transcript_n,
      TWAS_APA,
      WTTS_PAS_n,
      WTTS_APA = APA_multi_PAS,
      WTTS_DE_PAS = PAS,
      Apeglm_LFC,
      Std, Strong, Weak, HBFSS, Ovlp,
      Support,
      Std_BH,
      Strong_BH,
      Weak_BH,
      EmpP,
      HBFSS_score
    ) %>%
    dplyr::arrange(
      Comparison,
      factor(Analysis, levels = c(
        "Original (No EVS)",
        "NormEVS Lead", "NormEVS Rem",
        "RawEVS Lead", "RawEVS Rem"
      )),
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
      TWAS_n,
      TWAS_Transcript_n,
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
      factor(Analysis, levels = c(
        "Original (No EVS)",
        "NormEVS Lead", "NormEVS Rem",
        "RawEVS Lead", "RawEVS Rem"
      )),
      Rat_Ortholog
    )
}

build_twas_overlap_summary <- function(gene_table) {
  if (!nrow(gene_table)) return(data.frame())

  gene_table %>%
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
}

plot_twas_gene_support <- function(gene_table, figure_dir) {
  if (!nrow(gene_table)) return(invisible(NULL))

  analysis_map <- c(
    "Original (No EVS)"="Orig",
    "NormEVS Lead"="N-Lead",
    "NormEVS Rem"="N-Rem",
    "RawEVS Lead"="R-Lead",
    "RawEVS Rem"="R-Rem"
  )
  analysis_levels <- unname(analysis_map)

  method_rows <- list()
  for (method in c("Std", "Strong", "Weak", "HBFSS")) {
    keep <- !is.na(gene_table[[method]]) & gene_table[[method]]
    if (!any(keep)) next
    tmp <- gene_table[keep, , drop = FALSE]
    tmp$Method <- c(Std="Standard", Strong="Strong", Weak="Weak", HBFSS="HBFSS")[[method]]
    method_rows[[method]] <- tmp
  }
  if (!length(method_rows)) return(invisible(NULL))

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

  long$Gene_label <- factor(long$Gene_label, levels = rev(sort(unique(long$Gene_label))))

  p <- ggplot(long, aes(Analysis, Gene_label, color = Method, shape = Method)) +
    geom_point(position = position_dodge(width = 0.42), size = 2.7, alpha = 0.98, stroke = 0.85) +
    facet_wrap(~ Comparison, ncol = 2, scales = "free_y") +
    scale_color_manual(
      values = significance_method_colors,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method",
      guide = guide_legend(override.aes = list(size = 3.2, stroke = 0.95))
    ) +
    scale_shape_manual(
      values = significance_method_shapes,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method"
    ) +
    labs(
      title = paste0("3'aTWAS ortholog genes identified in WTTS-Seq | ", CUTOFF_METHOD_LABEL),
      x = NULL,
      y = "Rat ortholog [human TWAS]"
    ) +
    manuscript_theme() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      axis.text.y = element_text(size = 6.8),
      legend.position = "bottom",
      plot.margin = margin(10, 18, 10, 18)
    )

  max_genes <- max(table(long$Comparison))
  height <- min(24, max(8.0, 5.0 + 0.17 * max_genes))
  save_grob(p, file.path(figure_dir, "Figure_TWAS_Gene_Support.png"), width = 15.5, height = height)
  save_grob(p, file.path(figure_dir, "Figure_TWAS_Gene_Support.pdf"), width = 15.5, height = height)
  invisible(p)
}

plot_twas_view_counts <- function(summary_row, title) {
  if (!nrow(summary_row)) return(NULL)
  df <- data.frame(
    Method = factor(c("Standard", "Strong", "Weak", "HBFSS"), levels = significance_method_levels),
    n = c(summary_row$Std[1], summary_row$Strong[1], summary_row$Weak[1], summary_row$HBFSS[1]),
    stringsAsFactors = FALSE
  )

  ggplot(df, aes(Method, n, color = Method, shape = Method)) +
    geom_point(size = 3.2, stroke = 0.90) +
    geom_text(aes(label = n), vjust = -0.8, size = 3.0, show.legend = FALSE) +
    scale_color_manual(values = significance_method_colors, breaks = significance_method_levels,
                       labels = unname(significance_method_labels[significance_method_levels]), name = "Method") +
    scale_shape_manual(values = significance_method_shapes, breaks = significance_method_levels,
                       labels = unname(significance_method_labels[significance_method_levels]), name = "Method") +
    scale_x_discrete(labels = unname(significance_method_labels[significance_method_levels])) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
    labs(
      title = compact_title(title, width = 42),
      subtitle = paste0("Ovlp=", summary_row$Ovlp[1],
                        "  TWAS-APA=", summary_row$TWAS_APA[1],
                        "  WTTS-DE-APA=", summary_row$WTTS_DE_APA[1]),
      x = NULL,
      y = "TWAS genes"
    ) +
    manuscript_theme() +
    theme(legend.position = "none", plot.margin = margin(10, 14, 10, 14))
}

export_twas_by_view <- function(pas_table, gene_table, summary_table) {
  views <- comparison_analysis_views()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (i in seq_len(nrow(views))) {
      track_key <- views$track_key[i]
      dataset_key <- views$dataset_key[i]
      analysis_label <- analysis_view_label(track_key, dataset_key)
      paths <- analysis_view_paths(comparison_name, track_key, dataset_key)

      pas_sub <- if (nrow(pas_table)) {
        pas_table[
          pas_table$Comparison == comparison_name & pas_table$Analysis == analysis_label,
          ,
          drop = FALSE
        ]
      } else pas_table

      gene_sub <- if (nrow(gene_table)) {
        gene_table[
          gene_table$Comparison == comparison_name & gene_table$Analysis == analysis_label,
          ,
          drop = FALSE
        ]
      } else gene_table

      summary_sub <- summary_table[
        summary_table$Comparison == comparison_name & summary_table$Analysis == analysis_label,
        ,
        drop = FALSE
      ]

      save_csv(pas_sub, file.path(paths$tables, "Table_TWAS_PAS.csv"))
      save_csv(gene_sub, file.path(paths$tables, "Table_TWAS_Genes.csv"))

      if (isTRUE(EXPORT_INDIVIDUAL_VIEW_FIGURES)) {
        p <- plot_twas_view_counts(
          summary_sub,
          paste0(comparison_name, " | ", analysis_label, " | 3'aTWAS")
        )
        if (!is.null(p)) {
          save_grob(p, file.path(paths$figures, "TWAS_Overlap.png"), width = 7.2, height = 5.4)
          save_grob(p, file.path(paths$figures, "TWAS_Overlap.pdf"), width = 7.2, height = 5.4)
        }
      }
    }
  }

  invisible(TRUE)
}

save_twas_comparison_panels <- function(summary_table) {
  for (comparison_name in as.character(comparison_table$comparison_name)) {
    sub <- summary_table[summary_table$Comparison == comparison_name, , drop = FALSE]
    if (!nrow(sub)) next

    preferred_views <- c("Original (No EVS)", "NormEVS Lead", "NormEVS Rem", "RawEVS Lead", "RawEVS Rem")
    views <- preferred_views[preferred_views %in% unique(as.character(sub$Analysis))]
    plots <- lapply(views, function(v) {
      one <- sub[sub$Analysis == v, , drop = FALSE]
      if (!nrow(one)) return(NULL)
      short <- c(
        "Original (No EVS)"="Orig",
        "NormEVS Lead"="N-Lead",
        "NormEVS Rem"="N-Rem",
        "RawEVS Lead"="R-Lead",
        "RawEVS Rem"="R-Rem"
      )[[v]]
      plot_twas_view_counts(one, short)
    })
    plots <- Filter(Negate(is.null), plots)
    n_views <- length(plots)
    if (!n_views) next

    panel <- assemble_one_legend_panel(
      lapply(plots, function(p) p + theme(legend.position = "bottom")),
      panel_title = paste0(
        comparison_name, " | ", CUTOFF_METHOD_LABEL, " | 3'aTWAS overlap across ",
        n_views, " view", ifelse(n_views == 1L, "", "s")
      ),
      ncol = n_views
    )
    panel_dir <- file.path(output_dir, comparison_name, "Panels")
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
    base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_TWAS_", n_views, "Views"))
    save_grob(panel, paste0(base, ".png"), width = 20.0, height = 6.0)
    save_grob(panel, paste0(base, ".pdf"), width = 20.0, height = 6.0)
  }
  invisible(TRUE)
}

run_twas_overlap_analysis <- function() {
  twas_file <- resolve_existing_file(twas_file_candidates, "3'aTWAS file")
  message("Running final 3'aTWAS ortholog overlap: ", twas_file)

  twas <- read_twas_study(twas_file)
  mapping <- build_twas_ortholog_map(twas)
  rat_meta <- build_twas_rat_metadata(twas, mapping)
  all_results <- collect_twas_analysis_results()

  pas_table <- build_twas_overlap_pas_table(all_results, rat_meta)
  gene_table <- build_twas_overlap_gene_table(pas_table)
  summary_table <- build_twas_overlap_summary(gene_table)

  twas_dir <- file.path(output_dir, "TWAS")
  twas_fig_dir <- file.path(twas_dir, "figures")
  twas_tab_dir <- file.path(twas_dir, "tables")
  dir.create(twas_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(twas_tab_dir, recursive = TRUE, showWarnings = FALSE)

  save_csv(pas_table, file.path(twas_tab_dir, "Table_TWAS_Overlap_PAS.csv"))
  save_csv(gene_table, file.path(twas_tab_dir, "Table_TWAS_Overlap_Genes.csv"))
  save_csv(summary_table, file.path(twas_tab_dir, "Table_TWAS_Method_Counts.csv"))

  plot_twas_gene_support(gene_table, twas_fig_dir)
  export_twas_by_view(pas_table, gene_table, summary_table)
  save_twas_comparison_panels(summary_table)

  invisible(list(mapping = mapping, pas = pas_table, genes = gene_table, summary = summary_table))
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

empirical_cutoff_export <- comparison_table[, c("comparison_name", "group1_prefix", "group2_prefix", "fixed_5000_k"), drop = FALSE]
names(empirical_cutoff_export) <- c("Comparison", "RT_prefix", "ZT_prefix", "Fixed_5000_k")
empirical_cutoff_export$Active_Cutoff_Method <- CUTOFF_METHOD_LABEL
empirical_cutoff_export$Active_k_per_condition <- vapply(
  empirical_cutoff_export$Comparison,
  function(cn) tryCatch(get_active_evs_cutoff(cn), error = function(e) NA_integer_), integer(1)
)
empirical_cutoff_export$Basis <- EMPIRICAL_EVS_CUTOFF_BASIS
empirical_cutoff_export$Applied_to <- "NormEVS and RawEVS"
save_csv(
  empirical_cutoff_export,
  file.path(summary_table_dir, "Table_EVS_Cutoffs.csv")
)

evs_split_audit <- build_evs_split_audit_table()
if (nrow(evs_split_audit) > 0L) {
  save_csv(
    evs_split_audit,
    file.path(summary_table_dir, "Table_EVS_Split_Audit.csv")
  )
}

write_methods_note()
save_paper_volcano_panels()
save_paper_support_figures(overall_summary)
tryCatch(save_likelihood_and_dispersion_outputs(),
         error = function(e) warning("Likelihood and dispersion outputs skipped: ", conditionMessage(e), call. = FALSE))

# Final analytical stage: 3'aTWAS ortholog overlap with every reported WTTS method.
twas_results <- run_twas_overlap_analysis()

figure_zip <- create_all_figures_zip()
table_zip <- create_all_tables_zip()

cat("\n=====================================================\n")
cat("Pipeline complete.\n")
cat("Build: ", PIPELINE_BUILD, "\n", sep = "")
cat("Cutoff method: ", CUTOFF_METHOD_LABEL, "\n", sep = "")
cat("Repository root:\n", repo_root, "\n", sep = "")
cat("Output directory:\n", output_dir, "\n", sep = "")
cat("3'aTWAS overlap: complete\n")
cat("Figure ZIP:\n", figure_zip, "\n", sep = "")
cat("Table ZIP:\n", table_zip, "\n", sep = "")
cat("=====================================================\n\n")

if (nrow(overall_summary) > 0L) print(overall_summary)
