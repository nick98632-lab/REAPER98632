#!/usr/bin/env Rscript
# =============================================================================
# SEQUENCE MANUSCRIPT PIPELINE
# Eigenvector splitting (EVS) of WTTS-Seq polyadenylation sites (PASs) into a
# Leading Edge and a Remainder, each moderated by DESeq2 as its own dataset,
# with the EVS cutoff k chosen to give DESeq2 the tightest dispersion prior; apeglm, empirical-null
# calibration, higher criticism, HBFSS, and 3'aTWAS overlap.
# =============================================================================
#
# RUN
#   Rscript SEQUENCE_EVS_MANUSCRIPT.R
#   Optional environment variables:
#     SEQUENCE_REPO_ROOT  folder holding WTTS-Seq_2022.2_DE_raw_read_numbers.csv
#     SEQUENCE_MULTI_ROOT output folder (default <repo>/exports/sequence_evs_manuscript)
#
# The script runs as a master process. The master relaunches itself once per
# cutoff method (Fixed 5,000, then Optimal k), then writes the cross-method
# comparison (Figure S3) and one ZIP of everything.
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
# CHOICE OF k (fixed before any differential test; no DE p-value is used)
#   Fixed 5,000      the prespecified k of the paper.
#   Optimal k        per comparison: the k whose split gives DESeq2 the
#                    narrowest dispersion prior, pooled over Leading Edge and
#                    Remainder (minimum of W(k)). See OPTIMAL k below.
#
# DISPERSION TREND
# Every DESeq2 fit uses a0 + a1 x mean^(-b): DESeq2's parametric form (b = 1)
# with the exponent estimated, fitted on the log scale to match DESeq2's
# log-normal dispersion prior. See the section DESeq2 DISPERSION TREND below.
#
# DIFFERENTIAL EXPRESSION
# DESeq2 design ~ condition (trt = RT relative to untrt = ZT).
#   Standard  ordinary Wald, BH padj < 0.10 and |apeglm LFC| >= 1
#   Strong    greaterAbs(lfcThreshold = 1), BH padj < 0.10 and |apeglm LFC| >= 1
#   Weak      lessAbs(lfcThreshold = 1), BH padj < 0.20, |apeglm LFC| < 1, and HBFSS
#   HBFSS     |apeglm LFC| x [-log10(empirical p)] > Htau
# BH adjustment is DESeq2's own, within each view (with its independent
# filtering). Empirical p-values come from fdrtool (normal empirical null,
# cutoff 'fndr') applied to the Wald statistics, with ONE null per analysis:
# the Original's own statistics, or the pooled Leading Edge + Remainder
# statistics of an EVS analysis. Higher criticism is maximized over the lowest
# 10% of the analysis' ordered empirical p-values (alpha0 = 0.10) and must be
# positive; Htau = -log10(HCp) x 1, shared by the Leading Edge and Remainder.
# See the section DIFFERENTIAL EXPRESSION: MODEL STAGE AND SIGNIFICANCE STAGE.
#
# FIGURES (per cutoff method in <method>/Combined_Figures; master figures in
# Cross_Method_Comparison)
#   Figure 1   choosing k: W(k) and its minimum              (Optimal k run)
#   Figure 2   DESeq2 dispersions (gene-wise, trend, final), trend exponent b
#              and prior width: Original vs Leading Edge and Remainder (each run)
#   Figure S3  Fixed 5,000 vs Optimal k                      (master)
#   Figure S2  FDR / FNDR calibration: SEQUENCE_EVS_CALIBRATION.R
# Plus volcano, PCA/EVS evidence, HBFSS support and 3'aTWAS figures.
# =============================================================================

PIPELINE_BUILD <- "SEQUENCE_EVS_MANUSCRIPT_2026-11-17"


# =============================================================================
# CUTOFF METHODS, SHARED FIGURE STYLE, AND THE CROSS-METHOD FIGURE (S3)
# =============================================================================

COMPARISON_NAMES <- c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14")
METHOD_LEVELS <- c("Fixed 5,000", "Optimal k")
METHOD_SLUGS <- c("Fixed 5,000" = "Fixed_5000", "Optimal k" = "Optimal_k")
METHOD_KEYS <- c("Fixed 5,000" = "fixed5000", "Optimal k" = "optimal")
METHOD_COLS <- c("Fixed 5,000" = "#9A9A9A", "Optimal k" = "#1F6FB2")
ORIGINAL_LABEL <- "Original (No EVS)"
SET_COLS <- c("Original (No EVS)" = "#2B2B2B", "Leading Edge" = "#1F6FB2", "Remainder" = "#C0392B")

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

# Supplementary Figure S3: (a) the cutoff k chosen by each method,
# (b) DESeq2's pooled dispersion prior variance at each cutoff, (c) discoveries
# (reference only; never used to choose k).
save_figure_S3 <- function(multi_root, out_dir, manifest, discoveries) {
  if (!master_plotting_ready()) return(invisible(FALSE))
  cmp <- COMPARISON_NAMES
  ks <- do.call(rbind, lapply(METHOD_LEVELS, function(m) {
    r <- read_run_table(multi_root, m, "Table_Active_Cutoff.csv")
    if (is.null(r)) return(NULL)
    data.frame(Comparison = r$Comparison, Method = m, k = as.integer(r$k))
  }))
  ks$Comparison <- factor(ks$Comparison, levels = rev(cmp))
  ks$Method <- factor(ks$Method, levels = METHOD_LEVELS)
  pa <- ggplot(ks, aes(k, Comparison, colour = Method, shape = Method)) +
    geom_line(aes(group = Comparison), colour = "grey70", linewidth = 0.5) +
    geom_point(size = 2.4) +
    geom_text(aes(label = format(k, big.mark = ",")), vjust = -1.1, size = 2.0, show.legend = FALSE) +
    scale_colour_manual(values = METHOD_COLS) +
    scale_shape_manual(values = c("Fixed 5,000" = 15, "Optimal k" = 16)) +
    scale_y_discrete(labels = function(x) sub("_", " / ", x), expand = expansion(add = c(0.5, 0.9))) +
    scale_x_continuous(labels = scales::label_comma(), expand = expansion(mult = 0.12)) +
    labs(x = "k (top-|PC1| PASs per condition)", y = NULL, tag = "a", title = "The cutoff chosen",
         subtitle = "k per condition: Optimal k (minimum pooled prior variance) and the paper's 5,000.") +
    fig_theme() + theme(panel.grid.major.y = element_blank())

  pw <- do.call(rbind, lapply(METHOD_LEVELS, function(m) {
    d <- read_run_table(multi_root, m, "Table_Dispersion_After_EVS.csv")
    if (is.null(d)) return(NULL)
    do.call(rbind, lapply(cmp, function(cn) {
      x <- d[d$Comparison == cn, , drop = FALSE]
      o <- x[x$View == ORIGINAL_LABEL, , drop = FALSE]
      sp <- x[x$View %in% c("NormEVS Lead", "NormEVS Rem"), , drop = FALSE]
      if (!nrow(o) || nrow(sp) != 2L) return(NULL)
      data.frame(Comparison = cn, Series = m, W = sum(sp$n_PAS * sp$prior_var_before_floor) / sum(sp$n_PAS),
                 W0 = o$prior_var_before_floor[1])
    }))
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
      labs(x = NULL, y = "DESeq2 prior variance\n(pooled LE + REM)", tag = "b",
           title = "Dispersion prior width",
           subtitle = "Before DESeq2's 0.25 floor; NormEVS, n-weighted over LE and REM; Original for reference.") +
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
  fig <- assemble_rows(rows, heights = c(0.8, 1.1)[seq_along(rows)],
                       title = "Supplementary Figure S3. Fixed 5,000 vs Optimal k")
  save_figure(fig, out_dir, "Figure_S3_Fixed_vs_Optimal", 180, 135)
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
# MASTER PROCESS: one child run per cutoff method, cross-method outputs
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
                           unset = file.path(repo_guess, "exports", "sequence_evs_manuscript"))
  if (dir.exists(multi_root)) unlink(multi_root, recursive = TRUE, force = TRUE)
  dir.create(multi_root, recursive = TRUE, showWarnings = FALSE)

  count_path <- if (file.exists(file.path(repo_guess, count_name))) file.path(repo_guess, count_name) else
    file.path(repo_guess, "data", count_name)
  if (!file.exists(count_path)) stop("Cannot find ", count_name, ". Set SEQUENCE_REPO_ROOT.")

  cross_dir <- file.path(multi_root, "Cross_Method_Comparison")
  dir.create(cross_dir, recursive = TRUE, showWarnings = FALSE)

  # Fixed_5000 holds k; Optimal k is chosen inside its own run, per comparison,
  # before any split or test (see OPTIMAL k).
  cutoff_manifest <- data.frame(
    Cutoff_Method = unname(METHOD_SLUGS),
    RT0_ZT6 = c(5000L, NA_integer_), RT2_ZT8 = c(5000L, NA_integer_),
    RT4_ZT10 = c(5000L, NA_integer_), RT8_ZT14 = c(5000L, NA_integer_),
    Value = c("k per condition", "k per condition, chosen in its run"),
    Source = c("Prespecified k of the paper",
               "Minimum of the pooled DESeq2 dispersion prior variance W(k) over Leading Edge and Remainder"))
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

  discoveries <- collect_cross_method_discoveries(multi_root)
  utils::write.csv(discoveries, file.path(cross_dir, "Table_Cross_Method_Discoveries.csv"), row.names = FALSE)
  try_fig(save_figure_S3(multi_root, cross_dir, cutoff_manifest, discoveries), "Figure S3")

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
  cat("All runs complete: Fixed 5,000 and Optimal k.\n")
  act <- read_run_table(multi_root, "Optimal k", "Table_Active_Cutoff.csv")
  if (!is.null(act)) for (i in seq_len(nrow(act)))
    cat(sprintf("  %-9s Optimal k = %5d (Leading Edge %5d, pooled prior variance %.4f)\n",
                act$Comparison[i], act$k[i], act$Leading_Edge[i], act$Pooled_prior_variance_W[i]))
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
  "Optimal k" = "Per comparison, the k whose NormEVS split (union of RT and ZT top-k sets) minimizes the pooled DESeq2 dispersion prior variance over Leading Edge and Remainder; no DE p-values used"
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
# chosen before any test as the minimum of W(k) (Optimal k).
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
required_manifest_methods <- unname(METHOD_SLUGS)
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

# Optimal k per comparison, set by choose_optimal_k() at the start of each
# comparison's run; the search itself is kept for Figure 1 and its table.
OPTIMAL_K <- list()
optimal_k_search_rows <- list()

get_active_evs_cutoff <- function(comparison_name) {
  idx <- match(as.character(comparison_name), comparison_table$comparison_name)
  if (is.na(idx)) stop("No EVS cutoff is defined for comparison: ", comparison_name)
  k <- switch(
    CUTOFF_METHOD_KEY,
    fixed5000 = comparison_table$fixed_5000_k[idx],
    optimal = OPTIMAL_K[[as.character(comparison_name)]] %||%
      stop("Optimal k has not been set for ", comparison_name)
  )
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || !is.finite(k) || k < 1L) {
    stop("Invalid EVS cutoff for ", comparison_name, " under ", CUTOFF_METHOD_LABEL)
  }
  k
}


if (anyDuplicated(comparison_table$comparison_name)) stop("comparison_table contains duplicated comparison names.")
cutoff_cols <- "fixed_5000_k"
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
  unset = file.path(repo_root, "exports", "sequence_evs_manuscript")
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

  # Use the comparison-specific k exactly. Do not
  # silently shrink k or replace the rank rule with a loading-value cutoff.
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

assemble_one_legend_panel <- function(plot_list, panel_title, ncol = length(plot_list),
                                      legend = build_full_method_legend()) {
  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0L) {
    return(NULL)
  }

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
  for (tab in c("Table_EVS_Cutoffs.csv", "Table_Empirical_Null_By_Analysis.csv",
                "Table_Dispersion_After_EVS.csv", "Table_Active_Cutoff.csv",
                "Table_EVS_Composition.csv", "Table_Optimal_k_Search.csv")) {
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
  # rank descending, and select the comparison-specific k. No log
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
  # Both tracks use the same comparison-specific k so preprocessing
  # is the only difference in the selection rule. Downstream DESeq2 always
  # receives the corresponding raw-count subset.
  raw_count_matrix <- coerce_raw_count_matrix_for_deseq2(
    count_matrix,
    context = paste0(track_key, " full comparison matrix before EVS")
  )

  if (nrow(raw_count_matrix) < empirical_k) {
    stop(
      comparison_name, " ", track_key,
      ": EVS k=", empirical_k,
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
      ": EVS failed to select exactly k=", empirical_k,
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
      ": remainder dataset is empty after EVS k=", empirical_k, "."
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

# =============================================================================
# DESeq2 DISPERSION TREND: PARAMETRIC FORM, EXPONENT ESTIMATED, FITTED ON THE
# LOG SCALE
# =============================================================================
# DESeq2 shrinks each PAS's gene-wise dispersion toward a trend in the mean
# count, using the prior
#     log alpha_i ~ N( log alpha_tr(mean_i), sigma_d^2 ).
# Its default ("parametric") trend is
#     alpha_tr(mean) = a0 + a1 / mean                       (Love et al. 2014)
# fitted by an iterative gamma-family GLM of the gene-wise estimates on 1/mean.
# Two problems with that default appear in these data:
#
# 1. Shape. Nearly all well-expressed PASs have large |PC1| loadings and move
#    to the Leading Edge. The Remainder's dispersions fall faster than 1/mean
#    with no plateau, the fit returns a0 < 0, and DESeq2 silently switches that
#    group alone to a local-regression trend.
# 2. Height. The gamma GLM puts the curve through the arithmetic mean of the
#    dispersions at each count, which a minority of very noisy PASs pulls up.
#    About 64% of PASs then lie below their trend, in every view including the
#    Original, and are shrunk upward. The prior above is a normal distribution
#    for LOG dispersions, so its centre is the mean of the log dispersions.
#
# Every DESeq2 fit in this pipeline (Original, Leading Edge, Remainder, RawEVS
# views, calibration) therefore uses
#     alpha_tr(mean) = a0 + a1 * mean^(-b)
# with the exponent b estimated (b = 1 is DESeq2's form; dispersion as a power
# of the mean is the NB-P model, Di et al. 2011), and a0, a1 fitted on the log
# scale, consistent with DESeq2's prior:
#   * For a given b, a0 and a1 minimize the sum of squared log residuals,
#     sum_i [ log(dispGeneEst_i) - log(a0 + a1 * mean_i^(-b)) ]^2,
#     with a0, a1 >= TREND_MIN_COEF (so the trend is always positive). As in
#     DESeq2's own fit, PASs whose estimate is outside 1e-4 to 15 times the
#     current curve are set aside and the fit is repeated until the
#     coefficients stop changing (at most TREND_MAX_ITER + 1 passes).
#   * b minimizes the median absolute log residual, the criterion DESeq2's
#     authors recommend for comparing trend fits: every TREND_B_COARSE_STEP
#     over TREND_B_RANGE, then every TREND_B_FINE_STEP within one coarse step
#     of the best value.
# With this fit about half of the PASs lie on each side of their trend.
#
# The trend is handed to DESeq2 with dispersionFunction(), which computes the
# fitted dispersions, the residual variance of log dispersions and the prior
# variance as it does for its built-in trends. MAP shrinkage, the
# dispersion-outlier rule and the Wald test then run unchanged.
TREND_MIN_GENE_DISP <- 1e-6          # PASs used for the fit: dispGeneEst >= 100 x DESeq2's minDisp
TREND_B_RANGE <- c(0.5, 3.5)
TREND_B_COARSE_STEP <- 0.1
TREND_B_FINE_STEP <- 0.01
TREND_MIN_COEF <- 1e-8               # lower bound for a0 and a1
TREND_MAX_COEF <- c(a0 = 1e3, a1 = 1e8)
TREND_MAX_ITER <- 10L                 # refits after setting aside extreme PASs, as in DESeq2

# The trend as a function of the mean, holding only its three coefficients.
make_trend_function <- function(a0, a1, b) {
  force(a0); force(a1); force(b)
  trend <- function(means) a0 + a1 * means^(-b)
  attr(trend, "fitType") <- "parametric power trend, log-scale fit"
  attr(trend, "coefficients") <- c(a0 = a0, a1 = a1, b = b)
  trend
}

# a0 and a1 for one exponent b: least squares on the log scale, refitted after
# setting aside PASs outside 1e-4 to 15 times the current curve until the
# coefficients settle. Coefficients are optimized as log(a0), log(a1) so they
# stay positive. Returns the coefficients and the median absolute log residual
# over all PASs, or NULL if the fit does not settle.
fit_trend_at_exponent <- function(means, disps, b) {
  x <- means^(-b)
  y <- log(disps)
  lower <- rep(log(TREND_MIN_COEF), 2)
  upper <- log(TREND_MAX_COEF)
  theta <- log(c(0.1, 1))
  for (iter in seq_len(TREND_MAX_ITER + 1L)) {
    ratio <- disps / (exp(theta[1]) + exp(theta[2]) * x)
    good <- ratio > 1e-4 & ratio < 15
    if (sum(good) < 3L) return(NULL)
    xg <- x[good]
    yg <- y[good]
    sse <- function(t) sum((yg - log(exp(t[1]) + exp(t[2]) * xg))^2)
    grad <- function(t) {
      a <- exp(t)
      f <- a[1] + a[2] * xg
      r <- yg - log(f)
      c(-2 * sum(r / f) * a[1], -2 * sum(r * xg / f) * a[2])
    }
    opt <- tryCatch(stats::optim(theta, sse, grad, method = "L-BFGS-B", lower = lower, upper = upper),
                    error = function(e) NULL)
    if (is.null(opt) || any(!is.finite(opt$par))) return(NULL)
    old <- theta
    theta <- opt$par
    if (sum((theta - old)^2) < 1e-6) {
      a <- exp(theta)
      return(c(a0 = a[1], a1 = a[2], b = b,
               mar = stats::median(abs(y - log(a[1] + a[2] * x)))))
    }
  }
  NULL
}

# Best exponent: coarse grid, then a fine grid around the best value.
fit_dispersion_trend <- function(means, disps) {
  fits_for <- function(bs) Filter(Negate(is.null), lapply(bs, function(b) fit_trend_at_exponent(means, disps, b)))
  best_of <- function(fits) fits[[which.min(vapply(fits, function(f) f[["mar"]], numeric(1)))]]
  fits <- fits_for(seq(TREND_B_RANGE[1], TREND_B_RANGE[2], by = TREND_B_COARSE_STEP))
  if (!length(fits)) stop("the dispersion trend could not be fitted at any exponent.")
  best <- best_of(fits)
  fine <- seq(best[["b"]] - TREND_B_COARSE_STEP, best[["b"]] + TREND_B_COARSE_STEP, by = TREND_B_FINE_STEP)
  best <- best_of(c(list(best), fits_for(fine[fine > 0])))
  make_trend_function(best[["a0"]], best[["a1"]], best[["b"]])
}

# Fit the trend to a DESeqDataSet that already has gene-wise dispersions and
# store it with DESeq2's dispersionFunction() setter (which also fills dispFit,
# the residual variance and the dispersion prior variance).
set_dispersion_trend <- function(dds) {
  mc <- S4Vectors::mcols(dds)
  if (any(mc$allZero)) stop("remove PASs with no counts before fitting the dispersion trend.")
  use <- !mc$allZero & is.finite(mc$dispGeneEst) & mc$dispGeneEst >= TREND_MIN_GENE_DISP
  trend <- fit_dispersion_trend(mc$baseMean[use], mc$dispGeneEst[use])
  DESeq2::`dispersionFunction<-`(dds, value = trend)
}

# DESeq2's dispersion steps with the trend above: size factors (each dataset its
# own median-of-ratios), gene-wise dispersions, trend, MAP shrinkage. These are
# the steps DESeq() runs, with the trend fit replaced.
estimate_dispersions_with_trend <- function(dds) {
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- suppressMessages(DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE))
  dds <- set_dispersion_trend(dds)
  suppressMessages(DESeq2::estimateDispersionsMAP(dds, quiet = TRUE))
}

# =============================================================================
# DIFFERENTIAL EXPRESSION: MODEL STAGE (per analysis view) AND SIGNIFICANCE
# STAGE (per analysis)
# =============================================================================
# Each analysis view (Original; Leading Edge; Remainder) is its own DESeq2
# dataset with its own size factors, dispersion trend, prior and shrinkage
# (model stage, fit_view_model). Significance is then called once per
# ANALYSIS (significance stage, call_significance):
#   Original (No EVS)      one view  -> one empirical null
#   NormEVS (LE + REM)     two views -> one empirical null for the pair
#   RawEVS  (LE + REM)     two views -> one empirical null for the pair
# The Wald statistics of all views in an analysis are pooled, fdrtool fits one
# empirical null to the pool, and higher criticism sets one HBFSS threshold
# from the pooled empirical p-values. EVS therefore changes how dispersions are
# estimated, and the Original and EVS analyses are each judged against a null
# fitted to all of their PASs. (Fitting a separate null within the Leading Edge
# gave it a wider null than the Original's and set its HBFSS threshold above the
# Standard criterion for borderline sites.)

# Model stage: DESeq2 with the power trend, the three Wald tests, apeglm.
fit_view_model <- function(count_mat, coldata, dataset_name) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = coerce_raw_count_matrix_for_deseq2(
      count_mat,
      context = paste0(dataset_name, " DESeq2 input after EVS split")
    ),
    colData = coldata,
    design = make_design_formula(coldata)
  )
  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]

  # DESeq2 with the power trend fitted above, then the Wald test (no beta prior;
  # with 5 replicates per group DESeq() would not replace outliers either).
  dds <- estimate_dispersions_with_trend(dds)
  dds <- DESeq2::nbinomWaldTest(dds, quiet = TRUE)

  contrast <- c("condition", "trt", "untrt")
  res <- DESeq2::results(dds, contrast = contrast, alpha = BH_FDR_STANDARD, pAdjustMethod = "BH")
  res_strong <- DESeq2::results(dds, contrast = contrast, lfcThreshold = lfc_boundary,
                                altHypothesis = "greaterAbs", alpha = BH_FDR_STRONG, pAdjustMethod = "BH")
  res_weak <- DESeq2::results(dds, contrast = contrast, lfcThreshold = lfc_boundary,
                              altHypothesis = "lessAbs", alpha = BH_FDR_WEAK, pAdjustMethod = "BH")

  res_df <- as.data.frame(res)
  res_df$feature_id <- as.character(rownames(res_df))

  shr <- DESeq2::lfcShrink(dds, coef = get_condition_coef(dds), type = "apeglm")
  shr_df <- data.frame(feature_id = rownames(shr), lfc_shrunk = shr$log2FoldChange)
  res_df <- dplyr::left_join(res_df, shr_df, by = "feature_id")

  res_df$regulation_direction <- ifelse(
    is.na(res_df$lfc_shrunk), NA_character_,
    ifelse(res_df$lfc_shrunk > 0, "upregulated",
           ifelse(res_df$lfc_shrunk < 0, "downregulated", "no_change")))

  res_df$resGA_pvalue <- res_strong$pvalue[match(res_df$feature_id, rownames(res_strong))]
  res_df$resGA_padj <- res_strong$padj[match(res_df$feature_id, rownames(res_strong))]
  res_df$resLA_pvalue <- res_weak$pvalue[match(res_df$feature_id, rownames(res_weak))]
  res_df$resLA_padj <- res_weak$padj[match(res_df$feature_id, rownames(res_weak))]

  list(dds = dds, res_df = res_df, dataset_name = dataset_name)
}

# Significance stage for one analysis (a named list of model-stage fits): one
# fdrtool empirical null fitted to the pooled Wald statistics, one higher-
# criticism threshold from the pooled empirical p-values, then each view's
# HBFSS scores and calls. Returns the per-view results and a one-row summary of
# the shared null.
call_significance <- function(fits, analysis_label, annot_df) {
  valid <- lapply(fits, function(f) is.finite(f$res_df$stat))
  n_view <- vapply(valid, sum, integer(1))
  pooled <- unlist(Map(function(f, v) f$res_df$stat[v], fits, valid), use.names = FALSE)
  if (length(pooled) < 5L) {
    stop(sprintf("[%s] Fewer than 5 finite Wald statistics were available for fdrtool.", analysis_label))
  }
  fdr_fit <- run_empirical_null_fdrtool(pooled, dataset_name = analysis_label)
  if (length(fdr_fit$pval) != length(pooled)) {
    stop(sprintf("[%s] fdrtool empirical-p output length mismatch.", analysis_label))
  }
  emp_by_view <- split(as.numeric(fdr_fit$pval), factor(rep(seq_along(fits), n_view), levels = seq_along(fits)))

  hc_p <- safe_hc_thresh(fdr_fit$pval, dataset_name = analysis_label)
  hc_ok <- !is.na(hc_p) && is.finite(hc_p) && hc_p > 0 && hc_p <= 1
  htau <- if (hc_ok) abs(log10(hc_p)) * lfc_boundary else NA_real_
  param <- fdr_fit$param
  null_summary <- data.frame(
    n_wald_statistics = length(pooled),
    null_sd = if ("sd" %in% colnames(param)) unname(param[1, "sd"]) else NA_real_,
    null_fraction_eta0 = if ("eta0" %in% colnames(param)) unname(param[1, "eta0"]) else NA_real_,
    hc_p_threshold = hc_p,
    hbfss_threshold = htau)

  views <- list()
  for (i in seq_along(fits)) {
    f <- fits[[i]]
    f$res_df$empirical_p <- NA_real_
    f$res_df$empirical_p[valid[[i]]] <- emp_by_view[[i]]
    views[[names(fits)[i]]] <- classify_view(f, hc_p, htau, annot_df)
  }
  list(views = views, null = null_summary)
}

# Model stage for every view of one analysis, then the shared significance stage.
run_core_analysis <- function(count_mats, coldata, dataset_names, annot_df, analysis_label) {
  fits <- Map(function(m, nm) fit_view_model(m, coldata, nm), count_mats, dataset_names)
  names(fits) <- names(count_mats)
  call_significance(fits, analysis_label, annot_df)
}

# HBFSS scores and the Standard / Strong / Weak / HBFSS calls for one view,
# given its empirical p-values and the analysis' shared HC threshold.
classify_view <- function(fit, hc_p_threshold, hbfss_threshold, annot_df) {
  dds <- fit$dds
  res_df <- fit$res_df
  dataset_name <- fit$dataset_name

  res_df$HBFSS <- abs(res_df$lfc_shrunk) * (-log10(pmax(res_df$empirical_p, 1e-300)))
  res_df$HBFSS_core_pass <- if (is.na(hbfss_threshold)) FALSE else
    !is.na(res_df$HBFSS) & is.finite(res_df$HBFSS) & res_df$HBFSS > hbfss_threshold

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

  if (is.finite(hc_p_threshold) && !is.na(hc_p_threshold)) {
    expected_htau <- -log10(hc_p_threshold) * lfc_boundary
    if (!isTRUE(all.equal(hbfss_threshold, expected_htau, tolerance = 1e-12))) {
      stop("HBFSS threshold arithmetic validation failed for ", dataset_name)
    }
  }

  if (is.finite(hbfss_threshold) && !is.na(hbfss_threshold)) {
    expected_hbfss_flag <- finite_lfc & !is.na(res_df$HBFSS) &
      is.finite(res_df$HBFSS) & res_df$HBFSS > hbfss_threshold
    stopifnot(identical(res_df$hbfss_flag, expected_hbfss_flag))
  }

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
  final_df$hc_p_threshold_dataset <- hc_p_threshold
  final_df$hbfss_threshold_dataset <- hbfss_threshold

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

  list(dds = dds, results = final_df)
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
    ": EVS k=", evs$empirical_evs_k,
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

  # One analysis = the views that share one empirical null: the Original on
  # its own, and the Leading Edge with the Remainder.
  analyses <- list(split = c("leading_edge_dataset", "remainder_dataset"))
  if (!is.null(analysis_dataset_list$raw_dataset)) analyses <- c(list(original = "raw_dataset"), analyses)

  analysis_results <- list()
  for (a in names(analyses)) {
    keys <- analyses[[a]]
    analysis_label <- if (a == "original") ORIGINAL_LABEL else paste(unname(track_short[track_key]), "LE + REM")
    out <- run_core_analysis(
      count_mats = analysis_dataset_list[keys],
      coldata = coldata,
      dataset_names = paste(comparison_name, unname(track_short[track_key]), keys, sep = "_"),
      annot_df = annot_df,
      analysis_label = paste(comparison_name, analysis_label)
    )
    for (nm in keys) {
      analysis_results[[nm]] <- list(dds = out$views[[nm]]$dds, results = out$views[[nm]]$results,
                                     dataset_mat = analysis_dataset_list[[nm]])
    }
    empirical_null_rows[[length(empirical_null_rows) + 1L]] <<- cbind(
      data.frame(Comparison = comparison_name, Analysis = analysis_label), out$null)
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
      "PASs were ranked by absolute PC1 loading, and the comparison-specific k highest-loading PASs were selected independently from the RT and ZT condition eigenvectors."
    ),
    "PASs present in both condition-specific top-k sets were classified as Joint; PASs present in only one set were classified as Disjoint for that condition. Joint plus both Disjoint sets formed the Leading Edge, and all other PASs formed the Remainder. EVS defined PAS membership only. Downstream DESeq2 analyses received the corresponding raw-count subset and estimated normalization and dispersion parameters within that analysis view. The five reported views were Original (No EVS), NormEVS Lead, NormEVS Rem, RawEVS Lead, and RawEVS Rem.",
    "",
    "## Choice of the EVS cutoff k",
    sprintf("Two cutoffs were analyzed: k = 5,000, the prespecified cutoff of the original study, and an Optimal k chosen for each comparison before any differential test. EVS increases detection power through dispersion estimation: the Leading Edge and Remainder each receive their own mean-dispersion trend, and DESeq2 shrinks each PAS's gene-wise dispersion toward its own set's trend under a log-normal prior, log alpha_i ~ N(log alpha_tr(mean_i), sigma_d^2). The narrower this prior, the more information each PAS borrows and the more precise its final dispersion and Wald statistic. DESeq2 estimates sigma_d^2 as the squared median absolute deviation of log gene-wise dispersions about the trend minus their sampling variance, trigamma((m - p)/2) for m samples and p coefficients (Love et al. 2014). For a candidate k, the split was built exactly as in the analysis (NormEVS ranking; top-k PASs by absolute PC1 loading in each condition; Leading Edge = their union; Remainder = all other PASs), and each set was processed with the analysis' own DESeq2 steps (size factors, gene-wise dispersions and the power trend described below). The pooled prior variance W(k) = (n_LE sigma_LE^2 + n_REM sigma_REM^2)/(n_LE + n_REM) was computed before DESeq2's lower bound of 0.25, which would otherwise mask differences between splits, and the Optimal k was the value minimizing W(k). Too small a k leaves structured PASs in the Remainder and widens its prior; too large a k places unstructured PASs in the Leading Edge and widens its prior; the minimum lies where the two dispersion laws separate. To locate the global minimum, W(k) was evaluated at %d values of k evenly spaced on the logarithmic scale from %d to %.0f%% of the PASs in the comparison (or less, so that at least %.0f%% of the PASs remained in the Remainder), and the minimum was then refined between the neighbours of the best value by Brent's method (stats::optimize, tolerance %d). Only counts and dispersion estimates entered W(k); no condition effect, p-value or significance call was used. The same k was used for the NormEVS and RawEVS tracks.", OPTIMAL_K_GRID_POINTS, OPTIMAL_K_MIN, 100 * OPTIMAL_K_MAX_FRACTION, 100 * OPTIMAL_K_MIN_REMAINDER_FRACTION, as.integer(OPTIMAL_K_TOL)),
    "",
    "## DESeq2 model and log2 fold-change shrinkage",
    "Each analysis view was modeled independently with DESeq2 using a negative-binomial generalized linear model with design ~ condition and ZT/untrt as the reference. DESeq2 estimated median-of-ratios size factors, gene-wise dispersions, the mean-dispersion trend, final dispersions, and Wald statistics for the RT/trt coefficient. The mean-dispersion trend in every view, including the Original, was alpha_tr(mean) = a0 + a1 mean^(-b), DESeq2's parametric form (b = 1; Love et al. 2014) with the exponent estimated, and was fitted on the log scale. DESeq2's default fit failed for the Remainder (a0 < 0), because nearly all well-expressed PASs fall in the Leading Edge and the Remainder's dispersions decline faster than 1/mean without a plateau; DESeq2 would then have substituted a local-regression trend for that group only. DESeq2's gamma-GLM fit also places the trend at the arithmetic mean of the dispersions, above about two thirds of PASs, whereas its prior is normal on the log scale (log alpha_i ~ N(log alpha_tr, sigma_d^2)), whose centre is the mean of the log dispersions. Dispersion as a power of the mean is the NB-P model (Di et al. 2011). For each candidate b, a0 and a1 (both constrained to be at least 1e-8) minimized the sum of squared differences between log gene-wise dispersions and log trend; as in DESeq2's fit, PASs outside 1e-4 to 15 times the current curve were set aside and the fit repeated until the coefficients converged. b was chosen, every 0.1 from 0.5 to 3.5 and then every 0.01 around the best value, to minimize the median absolute residual of log gene-wise dispersion about log trend. With this fit about half of the PASs lay on each side of their trend. The trend was supplied to DESeq2 through dispersionFunction(), and DESeq2's prior-variance estimate, MAP shrinkage, dispersion-outlier rule and Wald test were applied unchanged. Size factors, gene-wise dispersions, MAP dispersions and Wald tests were computed with estimateSizeFactors, estimateDispersionsGeneEst, estimateDispersionsMAP and nbinomWaldTest, the steps DESeq() runs. Log2 fold changes were then shrunken with apeglm. The apeglm-shrunken log2 fold change was the reported effect estimate, the effect term in HBFSS, the direction indicator, and the x-coordinate for all significance figures.",
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
        "The empirical null was estimated once per analysis: once for the Original (No EVS) from its Wald statistics, and once for each EVS analysis from the pooled Wald statistics of its Leading Edge and Remainder, which were fitted as separate DESeq2 datasets. EVS therefore changed how dispersions were estimated, while the Original and each EVS analysis were judged against a null fitted to all of their PASs; a null fitted within the Leading Edge alone was wider than the Original's and made HBFSS stricter than the Standard criterion for borderline sites. ",
        "The finite ordinary DESeq2 Wald statistics of the analysis were supplied to fdrtool with statistic='normal' and cutoff.method='fndr' to estimate a zero-centered normal empirical null and the corresponding empirical-null p-values. The fitted null standard deviation and null fraction (eta0) were reported for every analysis. ",
        "These empirical p-values were distinct from the ordinary DESeq2 Wald p-values and from the greaterAbs/lessAbs p-values. ",
        "Higher-Criticism scores were then calculated from the analysis' ordered empirical-null p-values (pooled over Leading Edge and Remainder for EVS) with fdrtool::hc.score, giving one HCp and one Htau per analysis. ",
        "The HC maximization region was prespecified as the lowest %.0f%% of ordered empirical p-values (alpha0=%.2f), matching the lower-tail restriction provided by fdrtool::hc.thresh(alpha0=...). ",
        "Within that region, HCp was the empirical p-value at the maximum HC score. An analysis was assigned an HCp threshold only when that maximum HC score was strictly positive; if the maximum HC score was non-positive or no valid search point existed, HCp and Htau were left undefined and no PAS in that analysis was called HBFSS-significant."
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
    "Original (No EVS), NormEVS Lead, NormEVS Rem, RawEVS Lead, and RawEVS Rem were written to separate folders within each comparison for view-specific tables. Every value of k evaluated in the Optimal k search (Leading Edge and Remainder sizes, W(k), prior variance and trend exponent b of each set) was written to Table_Optimal_k_Search.csv and plotted in Figure 1 of the Optimal k run; Figure S3 was written to Cross_Method_Comparison by the master process. Manuscript figures were organized as comparison-level panels so related volcano, PCA/EVS, and 3'aTWAS results could be reviewed side by side without redundant individual images. The curated figure archive SEQUENCE_ALL_FIGURES.zip contained manuscript PNG panels only. The separate archive SEQUENCE_ALL_TABLES.zip contained the locked EVS cutoff table, the empirical-null table for every analysis, the dispersion-after-EVS tables, the EVS split audit table, overall differential-expression method-count table, quantitative EVS/PCA evidence table, one all-view significant-PAS table and method-count table for each comparison, and the combined 3'aTWAS PAS-, gene-, and method-summary tables."
  )

  writeLines(methods_text, file.path(output_dir, "Methods_Manuscript.md"), useBytes = TRUE)
  invisible(TRUE)
}


# =============================================================================
# DESeq2 DISPERSIONS AFTER THE SPLIT: SUMMARIES AND FIGURE 2
# =============================================================================
# Everything here is read from the fitted DESeq2 objects of each view; nothing
# is refitted.

DISP_MIN_GENE_EST <- 1e-6

# ---- per-run state -----------------------------------------------------------

dispersion_metrics_rows <- list()
dispersion_trend_rows <- list()
dispersion_pas_rows <- list()
composition_rows <- list()
active_cutoff_rows <- list()
empirical_null_rows <- list()

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
    prior_var_before_floor = if (length(v) == 1L) as.numeric(v) - trigamma((ncol(dds) - ncol(X)) / 2) else NA_real_,
    trend_a0 = attr(dfun, "coefficients")[["a0"]],
    trend_a1 = attr(dfun, "coefficients")[["a1"]],
    trend_b = attr(dfun, "coefficients")[["b"]],
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
      trend = pmin(DESeq2::dispersionFunction(dds)(grid_bm), max(10, ncol(dds))))  # DESeq2's dispersion ceiling
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

# ---- figure helpers ----------------------------------------------------------

comp_levels <- function() as.character(comparison_table$comparison_name)
log_lab <- function() scales::trans_format("log10", scales::math_format(10^.x))
as_comp <- function(x) factor(x, levels = comp_levels())

# ---- Figure 2: DESeq2's dispersions and the width of its prior ---------------

set_bars_scale <- function() {
  scale_fill_manual(values = SET_COLS, breaks = names(SET_COLS))
}

# DESeq2's dispersion plot for the Original, Leading Edge and Remainder:
# gene-wise estimates (grey), final MAP estimates (colour) and the fitted
# parametric trend a0 + a1 * mean^(-b) (orange, capped at DESeq2's maximum, 10).
panel_dispersion_grid <- function(pas, trd, met, tag, n_points = 3000L) {
  d <- pas[is.finite(pas$dispGeneEst) & pas$dispGeneEst >= DISP_MIN_GENE_EST & is.finite(pas$dispersion), ]
  set.seed(20261115L)
  d <- do.call(rbind, lapply(split(d, list(d$Comparison, d$View), drop = TRUE), function(x)
    x[sort(sample.int(nrow(x), min(nrow(x), n_points))), , drop = FALSE]))
  d$Set <- factor(unname(VIEW_SET[d$View]), levels = names(SET_COLS)); d$Comparison <- as_comp(d$Comparison)
  t <- trd; t$Set <- factor(unname(VIEW_SET[t$View]), levels = names(SET_COLS)); t$Comparison <- as_comp(t$Comparison)
  m <- met[met$View %in% names(VIEW_SET), ]
  m$Set <- factor(unname(VIEW_SET[m$View]), levels = names(SET_COLS)); m$Comparison <- as_comp(m$Comparison)
  m$label <- sprintf("n = %s\ntrend %.3f + %.2f x mean^-%.2f\nprior variance %.2f", format(m$n_PAS, big.mark = ","),
                     m$trend_a0, m$trend_a1, m$trend_b, m$prior_var_before_floor)
  ggplot(d, aes(baseMean)) +
    geom_point(aes(y = dispGeneEst), colour = "grey72", size = 0.12, alpha = 0.5) +
    geom_point(aes(y = dispersion, colour = Set), size = 0.12, alpha = 0.6) +
    geom_line(data = t, aes(baseMean, trend), colour = "#E8A317", linewidth = 0.6) +
    geom_text(data = m, aes(x = Inf, y = Inf, label = label), hjust = 1.03, vjust = 1.15, size = 1.8, lineheight = 0.95) +
    facet_grid(Set ~ Comparison, labeller = labeller(Comparison = comp_lab())) +
    scale_colour_manual(values = SET_COLS, guide = "none") +
    scale_x_log10(labels = log_lab()) + scale_y_log10(labels = log_lab()) +
    coord_cartesian(ylim = c(1e-3, 20)) +
    labs(x = "Mean of normalized counts", y = "Dispersion", tag = tag,
         title = "DESeq2 dispersion estimates: Original vs Leading Edge and Remainder",
         subtitle = "Grey: gene-wise estimates. Colour: final (shrunken) estimates. Orange: parametric trend a0 + a1 x mean^-b, fitted on the log scale.") +
    fig_theme() + theme(panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.3))
}

# Exponent b of the fitted trend a0 + a1 * mean^(-b), Original vs each NormEVS
# set: a different b means the Leading Edge and Remainder follow different
# dispersion laws.
panel_trend_exponent <- function(met, tag) {
  d <- met[met$View %in% names(VIEW_SET), c("Comparison", "View", "trend_b")]
  d$Set <- factor(unname(VIEW_SET[d$View]), levels = names(SET_COLS)); d$Comparison <- as_comp(d$Comparison)
  ggplot(d, aes(Comparison, trend_b, fill = Set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.74) +
    geom_text(aes(label = sprintf("%.2f", trend_b)), position = position_dodge(width = 0.8), vjust = -0.4, size = 1.8) +
    set_bars_scale() +
    scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Trend exponent b", tag = tag,
         title = "Shape of the dispersion trend",
         subtitle = "b in a0 + a1 x mean^-b: Leading Edge and Remainder follow different laws.") +
    fig_theme() + theme(panel.grid.major.x = element_blank())
}

# DESeq2's dispersion prior variance before its 0.25 floor (residual variance
# of log gene-wise dispersions about the trend minus their sampling variance)
# for the Original and each NormEVS set; W = n-weighted pool, the quantity the
# Optimal k minimizes.
panel_prior_width <- function(met, tag) {
  d <- met[met$View %in% names(VIEW_SET), c("Comparison", "View", "n_PAS", "prior_var_before_floor")]
  d$Set <- factor(unname(VIEW_SET[d$View]), levels = names(SET_COLS)); d$Comparison <- as_comp(d$Comparison)
  s <- do.call(rbind, lapply(split(d, d$Comparison, drop = TRUE), function(x) {
    sp <- x[x$Set != ORIGINAL_LABEL, ]
    w <- sum(sp$n_PAS * sp$prior_var_before_floor) / sum(sp$n_PAS)
    w0 <- x$prior_var_before_floor[x$Set == ORIGINAL_LABEL]
    data.frame(Comparison = x$Comparison[1], W = w, top = max(x$prior_var_before_floor),
               pct = if (length(w0)) 100 * (1 - w / w0) else NA_real_)
  }))
  s$x <- as.numeric(s$Comparison)
  s$lab <- ifelse(is.finite(s$pct), sprintf("W %.2f\n%.0f%% %s", s$W, abs(s$pct), ifelse(s$pct >= 0, "narrower", "wider")),
                  sprintf("W %.2f", s$W))
  ggplot(d, aes(Comparison, prior_var_before_floor, fill = Set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.74) +
    geom_segment(data = s, aes(x = x - 0.42, xend = x + 0.42, y = W, yend = W),
                 inherit.aes = FALSE, colour = "grey10", linewidth = 0.9) +
    geom_text(data = s, aes(x = x, y = top, label = lab), inherit.aes = FALSE, vjust = -0.3, size = 2.0,
              lineheight = 0.9, fontface = "bold", colour = "grey10") +
    set_bars_scale() +
    scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    labs(x = NULL, y = "DESeq2 prior variance\n(before the 0.25 floor)", tag = tag,
         title = "Width of the dispersion prior",
         subtitle = "Black line: Leading Edge and Remainder pooled (n-weighted), vs the Original.") +
    fig_theme() + theme(panel.grid.major.x = element_blank())
}

# ---- outputs ------------------------------------------------------------------

save_dispersion_outputs <- function() {
  out <- summary_table_dir
  pas <- if (length(dispersion_pas_rows)) dplyr::bind_rows(dispersion_pas_rows) else NULL
  trd <- if (length(dispersion_trend_rows)) dplyr::bind_rows(dispersion_trend_rows) else NULL
  met <- if (length(dispersion_metrics_rows)) dplyr::bind_rows(dispersion_metrics_rows) else NULL
  comp <- if (length(composition_rows)) dplyr::bind_rows(composition_rows) else NULL
  if (!is.null(met)) save_csv(met, file.path(out, "Table_Dispersion_After_EVS.csv"))
  if (length(empirical_null_rows)) {
    save_csv(dplyr::bind_rows(empirical_null_rows), file.path(out, "Table_Empirical_Null_By_Analysis.csv"))
  }
  if (!is.null(comp)) save_csv(comp, file.path(out, "Table_EVS_Composition.csv"))
  if (!is.null(trd)) save_csv(trd, file.path(out, "Table_Dispersion_Trends_By_View.csv"))
  if (!is.null(pas)) save_csv(pas, file.path(out, "Table_Dispersion_Estimates_By_View.csv"))
  if (length(active_cutoff_rows)) save_csv(dplyr::bind_rows(active_cutoff_rows), file.path(out, "Table_Active_Cutoff.csv"))
  if (length(optimal_k_search_rows)) {
    save_csv(dplyr::bind_rows(optimal_k_search_rows), file.path(out, "Table_Optimal_k_Search.csv"))
    try_fig(save_figure_1_optimal_k(), "Figure 1")
  }

  if (!is.null(pas) && !is.null(trd) && !is.null(met)) try_fig({
    fig2 <- assemble_rows(list(panel_dispersion_grid(pas, trd, met, "a"),
                               list(plots = list(panel_trend_exponent(met, "b"), panel_prior_width(met, "c")),
                                    widths = c(1, 1))),
                          heights = c(2.2, 0.9),
                          title = paste0("Figure 2. DESeq2 dispersions, trend shape and prior width, Original vs split (",
                                         CUTOFF_METHOD_LABEL, ")"))
    save_figure(fig2, paper_fig_dir, paste0("Figure_2_Dispersion_", CUTOFF_METHOD_SLUG), 180, 215)
  }, "Figure 2")
  invisible(TRUE)
}

# =============================================================================
# OPTIMAL k: THE SPLIT THAT GIVES DESeq2 THE TIGHTEST DISPERSION PRIOR
# =============================================================================
# EVS gains power in one way: the Leading Edge and the Remainder each get their
# own dispersion trend, and DESeq2 shrinks every PAS's gene-wise dispersion
# toward its own set's trend with the prior
#     log alpha_i ~ N( log alpha_tr(mean_i), sigma_d^2 ).
# The narrower this prior, the more information each PAS borrows from the
# others and the more precise its final dispersion and Wald test. DESeq2
# estimates its width as
#     sigma_d^2 = s_r^2 - trigamma((m - p) / 2),
# where s_r^2 is the squared MAD of log gene-wise dispersions about the trend
# (varLogDispEsts) and trigamma((m - p) / 2) their sampling variance (m samples,
# p model coefficients). DESeq2 then floors sigma_d^2 at 0.25 for safety; the
# criterion below uses the estimate before that floor, which would otherwise
# hide differences between splits.
#
# For a given k the split is built exactly as in the analysis (NormEVS
# ranking: top-k PASs by |PC1 loading| in RT and in ZT; Leading Edge = their
# union; Remainder = all other PASs), and each set runs the analysis' own
# DESeq2 steps up to the prior (size factors, gene-wise dispersions, the power
# trend of this pipeline). The pooled prior variance is
#     W(k) = [ n_LE sigma_LE^2 + n_REM sigma_REM^2 ] / (n_LE + n_REM).
# The Optimal k minimizes W(k). If k is too small, structured PASs stay in the
# Remainder and widen its prior; if k is too large, flat PASs enter the Leading
# Edge and widen its prior; the minimum is where the two dispersion laws
# separate. Only counts and DESeq2's dispersion estimates enter W(k): no
# condition effect, p-value or call is used, so k is fixed before any test.
#
# Search. W(k) need not be smooth or have a single minimum, so the minimum is
# first bracketed and then refined: W is evaluated at OPTIMAL_K_GRID_POINTS
# values of k evenly spaced on the log scale from OPTIMAL_K_MIN to
# OPTIMAL_K_MAX_FRACTION x the number of PASs, or the largest k that leaves at
# least OPTIMAL_K_MIN_REMAINDER_FRACTION of the PASs in the Remainder if that is
# smaller (plus k = 5,000, the paper's cutoff, for reference); Brent's method
# (stats::optimize) then refines the minimum between the two neighbours of the
# best of those values, to within OPTIMAL_K_TOL. The Optimal k is the evaluated k with the smallest W. The same
# k is used for the RawEVS track, exactly as k = 5,000 is.
OPTIMAL_K_MIN <- 250L
OPTIMAL_K_MAX_FRACTION <- 0.5
OPTIMAL_K_MIN_REMAINDER_FRACTION <- 0.1
OPTIMAL_K_GRID_POINTS <- 12L
OPTIMAL_K_TOL <- 10
OPTIMAL_K_REFERENCE <- 5000L

# Better |PC1| rank of each PAS over the two conditions (NormEVS ranking, as
# build_eigenvector_split uses). A PAS is in the Leading Edge at k exactly
# when this depth is <= k.
evs_union_depth <- function(raw, coldata) {
  rank_obj <- make_rank_matrix_for_track(raw, coldata, "normalized_evs")
  ids <- colnames(raw)
  rank_in <- function(cond) {
    lt <- compute_pc1_loading_table(rank_obj$rank_matrix, ids[coldata$condition == cond], top_n = 1L)$loading_table
    stats::setNames(lt$rank, lt$feature_id)[rownames(raw)]
  }
  depth <- pmin(rank_in("trt"), rank_in("untrt"))
  if (anyNA(depth)) stop("PC1 ranks are missing for some PASs.")
  depth
}

# DESeq2's dispersion prior for one set of PASs: the analysis' own steps
# (fit_view_model) up to the trend. Returns the set size, the prior variance
# before DESeq2's 0.25 floor, and the trend coefficients.
dispersion_prior_of_set <- function(raw_subset, coldata) {
  dds <- DESeq2::DESeqDataSetFromMatrix(raw_subset, colData = coldata, design = make_design_formula(coldata))
  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- suppressMessages(DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE))
  dds <- set_dispersion_trend(dds)
  f <- DESeq2::dispersionFunction(dds)
  X <- stats::model.matrix(DESeq2::design(dds), as.data.frame(SummarizedExperiment::colData(dds)))
  cf <- attr(f, "coefficients")
  c(n = nrow(dds), prior_var = as.numeric(attr(f, "varLogDispEsts")) - trigamma((ncol(dds) - ncol(X)) / 2),
    a0 = cf[["a0"]], a1 = cf[["a1"]], b = cf[["b"]])
}

choose_optimal_k <- function(comparison_name, count_matrix, coldata) {
  raw <- coerce_raw_count_matrix_for_deseq2(count_matrix, context = paste0(comparison_name, " optimal k"))
  depth <- evs_union_depth(raw, coldata)
  union_at <- cumsum(tabulate(depth, nrow(raw)))                 # Leading Edge size at each k
  fits_remainder <- function(k) nrow(raw) - union_at[k] >= OPTIMAL_K_MIN_REMAINDER_FRACTION * nrow(raw)
  k_max <- min(as.integer(floor(OPTIMAL_K_MAX_FRACTION * nrow(raw))),
               max(c(0L, which(union_at <= (1 - OPTIMAL_K_MIN_REMAINDER_FRACTION) * nrow(raw)))))
  if (k_max <= OPTIMAL_K_MIN) stop(comparison_name, ": too few PASs (", nrow(raw), ") to search for k.")
  evaluated <- list()
  W <- function(k) {
    k <- as.integer(round(k))
    key <- as.character(k)
    if (is.null(evaluated[[key]])) {
      le <- depth <= k
      a <- dispersion_prior_of_set(raw[le, , drop = FALSE], coldata)
      r <- dispersion_prior_of_set(raw[!le, , drop = FALSE], coldata)
      evaluated[[key]] <<- data.frame(
        Comparison = comparison_name, k = k, Leading_Edge = as.integer(a[["n"]]), Remainder = as.integer(r[["n"]]),
        W = (a[["n"]] * a[["prior_var"]] + r[["n"]] * r[["prior_var"]]) / (a[["n"]] + r[["n"]]),
        prior_var_LE = a[["prior_var"]], prior_var_REM = r[["prior_var"]],
        b_LE = a[["b"]], b_REM = r[["b"]], Search_step = "bracket")
    }
    evaluated[[key]]$W
  }
  grid <- unique(as.integer(round(exp(seq(log(OPTIMAL_K_MIN), log(k_max), length.out = OPTIMAL_K_GRID_POINTS)))))
  w_grid <- vapply(grid, W, numeric(1))
  i <- which.min(w_grid)
  if (i %in% c(1L, length(grid))) warning(comparison_name, ": the minimum of W(k) is at the edge of the search interval (k = ",
                                          grid[i], ").", call. = FALSE)
  n_bracket <- length(evaluated)
  stats::optimize(W, interval = c(grid[max(1L, i - 1L)], grid[min(length(grid), i + 1L)]), tol = OPTIMAL_K_TOL)
  refined <- names(evaluated)[seq_along(evaluated) > n_bracket]
  for (key in refined) evaluated[[key]]$Search_step <- "refine"
  if (OPTIMAL_K_REFERENCE <= nrow(raw) && fits_remainder(OPTIMAL_K_REFERENCE)) {
    ref_new <- is.null(evaluated[[as.character(OPTIMAL_K_REFERENCE)]])
    W(OPTIMAL_K_REFERENCE)
    if (ref_new) evaluated[[as.character(OPTIMAL_K_REFERENCE)]]$Search_step <- "reference (paper)"
  }
  s <- do.call(rbind, evaluated)
  s <- s[order(s$k), , drop = FALSE]
  rownames(s) <- NULL
  k <- s$k[which.min(s$W)]
  s$Chosen <- s$k == k
  OPTIMAL_K[[comparison_name]] <<- as.integer(k)
  optimal_k_search_rows[[comparison_name]] <<- s
  ch <- s[s$Chosen, ]
  message(sprintf("  %s: optimal k = %d (Leading Edge %d, Remainder %d; pooled prior variance %.4f, b LE %.2f / REM %.2f; %d evaluations)",
                  comparison_name, k, ch$Leading_Edge, ch$Remainder, ch$W, ch$b_LE, ch$b_REM, nrow(s)))
  invisible(k)
}

# Figure 1: W(k) for each comparison, its minimum and the paper's 5,000.
save_figure_1_optimal_k <- function() {
  s <- dplyr::bind_rows(optimal_k_search_rows)
  if (!nrow(s)) return(invisible(FALSE))
  s$Comparison <- as_comp(s$Comparison)
  ch <- s[s$Chosen, , drop = FALSE]
  ref <- s[s$k == OPTIMAL_K_REFERENCE, , drop = FALSE]
  lab <- merge(ch[, c("Comparison", "k", "W", "Leading_Edge")],
               data.frame(Comparison = ref$Comparison, W_ref = ref$W), all.x = TRUE)
  lab$label <- sprintf("k = %s (Leading Edge %s)\nW = %.3f%s", format(lab$k, big.mark = ","),
                       format(lab$Leading_Edge, big.mark = ","), lab$W,
                       ifelse(is.finite(lab$W_ref), sprintf("; at 5,000: %.3f", lab$W_ref), ""))
  col <- unname(METHOD_COLS[["Optimal k"]])
  p <- ggplot(s, aes(k, W)) +
    geom_vline(xintercept = OPTIMAL_K_REFERENCE, colour = unname(METHOD_COLS[["Fixed 5,000"]]), linetype = "22", linewidth = 0.4) +
    geom_line(colour = col, linewidth = 0.5) +
    geom_point(colour = col, size = 0.9) +
    geom_point(data = ch, colour = col, size = 2.6, shape = 21, fill = "white", stroke = 0.9) +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = label), hjust = -0.03, vjust = 1.2, size = 2.0,
              colour = col, fontface = "bold", lineheight = 0.95) +
    facet_wrap(~ Comparison, nrow = 1, scales = "free_y", labeller = labeller(Comparison = comp_lab())) +
    scale_x_log10(labels = scales::label_comma(), expand = expansion(mult = c(0.25, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    labs(x = "k (top-|PC1| PASs per condition; log scale)", y = "Pooled DESeq2 prior variance W(k)",
         title = "Figure 1. Choosing k: the split that gives DESeq2 the tightest dispersion prior",
         subtitle = paste0("W(k): size-weighted prior variance of log dispersion about each set's own trend (Leading Edge, Remainder), ",
                           "before DESeq2's 0.25 floor. Open circle: minimum (Optimal k). Dashed: the paper's 5,000. No p-values used.")) +
    fig_theme() + theme(panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.3))
  save_figure(p, paper_fig_dir, "Figure_1_Choosing_k", 180, 70)
  invisible(TRUE)
}

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  dir.create(file.path(output_dir, comparison_name), recursive = TRUE, showWarnings = FALSE)

  # The cutoff is fixed before any split or test: 5,000 (prespecified), or the
  # k that minimizes the pooled DESeq2 prior variance (Optimal k).
  if (identical(CUTOFF_METHOD_KEY, "optimal")) choose_optimal_k(comparison_name, count_matrix, coldata)
  ch <- if (identical(CUTOFF_METHOD_KEY, "optimal")) {
    s <- optimal_k_search_rows[[comparison_name]]
    s[s$Chosen, , drop = FALSE]
  } else NULL
  active_cutoff_rows[[comparison_name]] <<- data.frame(
    Comparison = comparison_name, Cutoff_Method = CUTOFF_METHOD_LABEL, k = get_active_evs_cutoff(comparison_name),
    Leading_Edge = if (is.null(ch)) NA_integer_ else ch$Leading_Edge,
    Pooled_prior_variance_W = if (is.null(ch)) NA_real_ else ch$W)
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
        ncol = 3,
        legend = grid::textGrob(
          "Blue triangles: RT samples. Grey circles: ZT samples. X: condition centroid; lines join each sample to its centroid.",
          gp = grid::gpar(fontsize = 11, col = "grey20"))
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
tryCatch(save_dispersion_outputs(),
         error = function(e) warning("Dispersion outputs skipped: ", conditionMessage(e), call. = FALSE))

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
