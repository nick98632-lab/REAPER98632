#!/usr/bin/env Rscript
# =============================================================================
# SEQUENCE MANUSCRIPT PIPELINE
# Eigenvector splitting (EVS) of WTTS-Seq polyadenylation sites (PASs) into a
# Leading Edge and a Remainder, each moderated by DESeq2 as its own dataset.
# One EVS cutoff k* per comparison, chosen without any significance test: two
# dispersion trends are fitted on the samples outside the comparison, and k*
# places the PC1 split where it best matches them. Then one differential-
# expression analysis at that k* (DESeq2 Wald tests, apeglm, empirical-null
# calibration, higher criticism, HBFSS) and 3'aTWAS overlap.
# =============================================================================
#
# RUN
#   Rscript SEQUENCE_EVS_MANUSCRIPT.R
#   Optional environment variables:
#     SEQUENCE_REPO_ROOT       folder holding WTTS-Seq_2022.2_DE_raw_read_numbers.csv
#     SEQUENCE_MULTI_ROOT      output folder (default <repo>/exports/sequence_evs_manuscript)
#     SEQUENCE_NULL_METHOD     fndr (default), pct0 or theoretical (HBFSS null)
#     SEQUENCE_PCA_METHODS     EVS PCA preprocessing methods to run, comma-
#                              separated (default VST,CPM; the first is the
#                              primary analysis; NORMALIZED is also accepted)
#     EVS_PCA_METHOD           set by the master for each run (VST, CPM or
#                              NORMALIZED; default VST when run directly)
#
# The script runs as a master process. The master relaunches itself once per
# PCA preprocessing method (one BLAS/OpenMP thread each; outputs in
# EVS_PCA_<METHOD>), then writes the cross-method tables
# (EVS_PCA_Method_Comparison/) and one ZIP of everything.
#
# ANALYSIS UNIT
# Each OrigID is analyzed as an individual PAS. Gene symbols are annotation
# only and are never used to collapse PASs before testing.
#
# COMPARISONS (5 RT vs 5 ZT samples each)
#   RT0 vs ZT6, RT2 vs ZT8, RT4 vs ZT10, RT8 vs ZT14
#
# ANALYSIS VIEWS PER COMPARISON
#   Original (No EVS), NormEVS Lead, NormEVS Rem
#
# EIGENVECTOR SPLITTING
# The experiment-wide count matrix (all 40 samples, every PAS with a nonzero
# count, before any split) is transformed once by one PCA preprocessing method
# (EVS_PCA_METHOD):
#   VST         DESeq2 varianceStabilizingTransformation(blind = TRUE) with
#               experiment-wide median-of-ratios size factors (default)
#   CPM         log1p(counts per million), library sizes = column totals
#   NORMALIZED  experiment-wide DESeq2 median-of-ratios normalized counts
# These matrices rank PASs only; none is ever passed to DESeq2. The matrix is
# used as is (no further log). Within each arm (condition), prcomp (centred,
# unscaled) is applied to the PAS-by-sample matrix; PASs are ranked by
# |PC1 loading| and the top k* are taken
# from each condition. PASs in both top-k sets are Joint, in one set Disjoint;
# Joint plus Disjoint form the Leading Edge (LE) and all other PASs the
# Remainder (REM). LE and REM are separate datasets: DESeq2 receives their raw
# counts with the comparison's shared size factors and
# estimates its own dispersion trend and prior. The PCA preprocessing changes
# only which PASs form the Leading Edge; DESeq2 always receives raw counts, the
# same shared size factors and the same negative-binomial model.
#
# CHOICE OF k*
#   Per comparison, without any test: a two-component mixture of mean-
#   dispersion trends (a0 + a1 x mean^-b each) is fitted by EM to DESeq2
#   gene-wise dispersions of the 30 samples outside the comparison, giving each
#   PAS the likelihood ratio of fitting the Leading-Edge vs the Remainder
#   trend; k* maximizes the summed log likelihood ratio of the PC1 split's
#   Leading Edge, S(k), along the PC1 union depth. The regime cutoff (shared Remainder / Divergence /
#   Leading-Edge rank regimes) is reported as a reference. Significance testing
#   runs once, after k* is chosen. See TWO-TREND CUTOFF k* and REGIME CUTOFF k*.
#
# DISPERSION TREND
# Every DESeq2 fit uses a0 + a1 x mean^(-b): DESeq2's parametric form (b = 1)
# with the exponent estimated, fitted on the log scale to match DESeq2's
# log-normal dispersion prior and then corrected for the downward sampling bias
# of log gene-wise dispersions. See the section DESeq2 DISPERSION TREND below.
#
# DIFFERENTIAL EXPRESSION
# DESeq2 design ~ condition (trt = RT relative to untrt = ZT).
#   Standard  ordinary Wald, BH padj < 0.10 and |apeglm LFC| >= 1
#   Strong    greaterAbs(lfcThreshold = 1), BH padj < 0.10 and |apeglm LFC| >= 1
#   Weak      lessAbs(lfcThreshold = 1), BH padj < 0.20, |apeglm LFC| < 1, and HBFSS
#   HBFSS     |apeglm LFC| x [-log10(empirical p)] > Htau
# Multiple testing (Standard, Strong, Weak): one BH correction per analysis and
# test, over the Original's raw p-values, or over the Leading Edge and Remainder
# raw p-values combined, without independent filtering; Cook's outlier PASs are
# identified once, in the Original's fit, and their p-values are removed in
# every view. The Original and the EVS views therefore test the same PASs under
# the same rules. Empirical p-values
# come from fdrtool (normal empirical null; cutoff SEQUENCE_NULL_METHOD, default
# 'fndr') applied to the Wald statistics of each view separately: the
# Original, the Leading Edge and the Remainder each have their own null. Higher
# criticism is maximized over the lowest 10% (alpha0 = 0.10) of the ordered
# empirical p-values of the comparison's Original and must be positive;
# Htau = -log10(HCp), one threshold per comparison, applied to the Original
# and to the Leading Edge and Remainder alike.
# See the section DIFFERENTIAL EXPRESSION: MODEL STAGE AND SIGNIFICANCE STAGE.
#
# FIGURES (per run in EVS_PCA_<METHOD>/Combined_Figures; master figures in
# EVS_PCA_Method_Comparison/)
#   Figure 1   two-trend cutoff: held-out dispersions and the two trends, the
#              share of PASs on the Leading-Edge trend along the PC1 union
#              depth, and the allocation score S(k) with k*        (per run)
#   Figure 2   DESeq2 dispersions (gene-wise, trend, final), trend exponent b
#              and prior width: Original vs Leading Edge and Remainder (per run)
#   Figure 4   PCA methods compared: k*, share of PASs on their own trend,
#              Leading-Edge overlap and discoveries vs the Original  (master)
#   Figure S1  rank regimes: D(r) of the eight arms, shared knots c1, c2
#              (master, all methods; per run)
#   Figure S3  regime reference cutoff: Pareto frontier and regime k* (per run)
#   Figure S2  FDR / FNDR calibration: SEQUENCE_EVS_CALIBRATION.R
# Plus volcano, PCA/EVS evidence, HBFSS support and 3'aTWAS figures.
# =============================================================================

PIPELINE_BUILD <- "SEQUENCE_EVS_MANUSCRIPT_2026-11-22_TWOTREND"


# =============================================================================
# CUTOFF METHOD AND SHARED FIGURE STYLE
# =============================================================================

COMPARISON_NAMES <- c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14")
# One cutoff method: held-out two-trend allocation (TWO-TREND CUTOFF k*).
METHOD_SLUGS <- c("Two-trend k*" = "TwoTrend_k")
METHOD_KEYS <- c("Two-trend k*" = "twotrend")

# EVS PCA preprocessing methods (see EIGENVECTOR SPLITTING). One run per method;
# each run fits its own regimes and chooses its own k*.
PCA_METHOD_LEVELS <- c("VST", "CPM", "NORMALIZED")
PCA_METHOD_DEFAULT <- c("VST", "CPM")
PCA_METHOD_LABELS <- c(VST = "experiment-wide DESeq2 VST (blind = TRUE)",
                       CPM = "experiment-wide log1p(CPM)",
                       NORMALIZED = "experiment-wide DESeq2 median-of-ratios normalized counts")
pca_run_folder <- function(method) paste0("EVS_PCA_", method)
EVS_PCA_METHOD <- toupper(trimws(Sys.getenv("EVS_PCA_METHOD", unset = "VST")))
if (!EVS_PCA_METHOD %in% PCA_METHOD_LEVELS) {
  stop("EVS_PCA_METHOD must be one of: ", paste(PCA_METHOD_LEVELS, collapse = ", "))
}
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

try_fig <- function(expr, what) {
  tryCatch(expr, error = function(e) {
    warning("Figure skipped (", what, "): ", conditionMessage(e), call. = FALSE)
    invisible(NULL)
  })
}

REGIME_BANDS <- c(Remainder = "#E8ECF2", Divergence = "#FBF0D0", "Leading edge" = "#DCEFE4")
REGIME_EDGES <- c(Remainder = "#5B6672", Divergence = "#8A6D12", "Leading edge" = "#2E8B62")

# D(r) of every arm with the across-arm median, the shared two-knot fit and the
# three regimes; one panel per PCA method. Used by each run and by the master.
plot_regime_definition <- function(curves, knots) {
  med <- stats::aggregate(cbind(D, Fitted) ~ EVS_PCA_Method + Rank, data = curves, FUN = stats::median)
  bands <- do.call(rbind, lapply(seq_len(nrow(knots)), function(i) data.frame(
    EVS_PCA_Method = knots$EVS_PCA_Method[i], Regime = factor(names(REGIME_BANDS), levels = names(REGIME_BANDS)),
    xmin = c(1, knots$c1[i], knots$c2[i]), xmax = c(knots$c1[i], knots$c2[i], knots$N[i]),
    n = c(knots$c1[i] - 1, knots$c2[i] - knots$c1[i] + 1, knots$N[i] - knots$c2[i]))))
  bands$lab <- sprintf("%s\n%s", bands$Regime, format(bands$n, big.mark = ","))
  top <- max(curves$D, na.rm = TRUE); bot <- min(curves$D, na.rm = TRUE); pad <- top - bot
  kl <- do.call(rbind, lapply(seq_len(nrow(knots)), function(i) data.frame(
    EVS_PCA_Method = knots$EVS_PCA_Method[i], x = c(knots$c1[i], knots$c2[i]),
    lab = sprintf(c("c1 = %s", "c2 = %s"), format(c(knots$c1[i], knots$c2[i]), big.mark = ",")))))
  key <- c("Individual arms", "Median observed", "Shared two-knot fit")
  ggplot() +
    geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Regime), alpha = 1) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dotted", linewidth = 0.2) +
    geom_line(data = curves, aes(Rank, D, group = Arm, colour = key[1]), linewidth = 0.15) +
    geom_line(data = med, aes(Rank, D, colour = key[2]), linewidth = 0.55) +
    geom_line(data = med, aes(Rank, Fitted, colour = key[3]), linewidth = 0.45, linetype = "22") +
    geom_vline(data = kl, aes(xintercept = x), colour = "grey35", linetype = "22", linewidth = 0.3) +
    geom_text(data = bands, aes(x = (xmin + xmax) / 2, y = top + 0.30 * pad, label = lab),
              size = 1.9, lineheight = 0.9, fontface = "bold", colour = REGIME_EDGES[as.character(bands$Regime)]) +
    geom_text(data = kl, aes(x = x, y = top + 0.08 * pad, label = lab), size = 1.8, colour = "grey25",
              hjust = c(1.05, -0.05)[(seq_len(nrow(kl)) - 1) %% 2 + 1]) +
    facet_wrap(~ EVS_PCA_Method, nrow = 1, labeller = labeller(EVS_PCA_Method = function(x) paste0(x, "-EVS"))) +
    scale_fill_manual(values = REGIME_BANDS, guide = "none") +
    scale_colour_manual(values = stats::setNames(c("#9AA3AD", "#CC79A7", "black"), key), name = NULL) +
    scale_x_continuous(labels = scales::label_comma(), expand = c(0, 0)) +
    coord_cartesian(ylim = c(bot - 0.04 * pad, top + 0.40 * pad)) +
    labs(x = "PAS rank by |PC1 loading| within arm (ascending)", y = "D(r) = F_E(r) - F_P(r)",
         title = "Remainder, Divergence and Leading-Edge rank regimes",
         subtitle = "Cumulative excess-variance mass minus cumulative PC1 mass, all eight arms; dashed: shared two-knot fit.") +
    fig_theme() + theme(legend.position = "bottom")
}

read_run_table <- function(multi_root, pca_method, file) {
  f <- file.path(multi_root, pca_run_folder(pca_method), "Summary_Tables", file)
  if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) else NULL
}

# Cross-method tables (master): per comparison and PCA method, the selected k,
# set sizes, cutoff criteria and discoveries with their differences from the
# Original; the overlap of the PC1 rankings between methods; cross-method checks.
write_pca_method_comparison <- function(multi_root, methods) {
  out_dir <- file.path(multi_root, "EVS_PCA_Method_Comparison")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  calls <- c(Standard = "Std", Strong = "Strong", Weak = "Weak", HBFSS = "HBFSS")
  rows <- list(); ranks <- list(); val <- list(); orig <- list()
  addv <- function(check, cmp, ok, detail) val[[length(val) + 1L]] <<- data.frame(
    Check = check, Comparison = cmp, Status = if (isTRUE(ok)) "PASS" else "FAIL", Detail = detail)
  for (pm in methods) {
    act <- read_run_table(multi_root, pm, "Table_Active_Cutoff.csv")
    cnt <- read_run_table(multi_root, pm, "Table_DE_Method_Counts.csv")
    rk <- read_run_table(multi_root, pm, "Table_EVS_PCA_Ranking.csv")
    if (is.null(act) || is.null(cnt)) { warning("Missing run tables for ", pm, call. = FALSE); next }
    if (!is.null(rk)) ranks[[pm]] <- rk
    for (cmp in unique(act$Comparison)) {
      a <- act[act$Comparison == cmp, , drop = FALSE]
      addv("Exactly one k selected per comparison", cmp, nrow(a) == 1L && is.finite(a$k[1]),
           sprintf("%s: k = %s", pm, paste(a$k, collapse = ",")))
      o <- cnt[cnt$Comparison == cmp & cnt$Analysis == "Original (No EVS)", , drop = FALSE]
      e <- cnt[cnt$Comparison == cmp & cnt$Analysis %in% c("NormEVS Lead", "NormEVS Rem"), , drop = FALSE]
      orig[[paste(cmp, pm)]] <- data.frame(Comparison = cmp, Method = pm, o[, unname(calls), drop = FALSE])
      d <- data.frame(Comparison = cmp, EVS_PCA_Method = pm, PCA_Preprocessing = PCA_METHOD_LABELS[[pm]],
                      Selected_k = a$k[1], Leading_Edge = a$Leading_Edge[1], Remainder = a$Remainder[1],
                      Expected_On_Trend = a$Expected_On_Trend[1],
                      Expected_On_Trend_All_Remainder = a$Expected_On_Trend_All_Remainder[1],
                      Regime_k = a$Regime_k[1], Expected_On_Trend_Regime_k = a$Expected_On_Trend_Regime_k[1],
                      Expected_On_Trend_k5000 = a$Expected_On_Trend_k5000[1],
                      LE_trend_b = a$LE_trend_b[1], REM_trend_b = a$REM_trend_b[1])
      for (cl in names(calls)) {
        ov <- if (nrow(o)) sum(o[[calls[[cl]]]], na.rm = TRUE) else NA_real_
        ev <- if (nrow(e)) sum(e[[calls[[cl]]]], na.rm = TRUE) else NA_real_
        d[[paste0("Original_", cl)]] <- ov
        d[[paste0("EVS_", cl)]] <- ev
        d[[paste0("Diff_", cl)]] <- ev - ov
      }
      rows[[length(rows) + 1L]] <- d
    }
  }
  comparison <- do.call(rbind, rows)
  if (!is.null(comparison)) {
    comparison <- comparison[order(comparison$Comparison, match(comparison$EVS_PCA_Method, methods)), , drop = FALSE]
    utils::write.csv(comparison, file.path(out_dir, "EVS_PCA_Method_Comparison.csv"), row.names = FALSE)
  }
  if (length(orig) > 1L) {
    ot <- do.call(rbind, orig)
    for (cmp in unique(ot$Comparison)) {
      x <- ot[ot$Comparison == cmp, unname(calls), drop = FALSE]
      addv("Original (No EVS) discoveries identical across PCA methods", cmp,
           nrow(unique(x)) == 1L, paste(ot$Method[ot$Comparison == cmp], collapse = ", "))
    }
  }
  overlap <- list()
  if (length(ranks) > 1L) {
    pairs <- utils::combn(names(ranks), 2L, simplify = FALSE)
    for (cmp in unique(unlist(lapply(ranks, `[[`, "Comparison")))) for (p in pairs) {
      A <- ranks[[p[1]]][ranks[[p[1]]]$Comparison == cmp, , drop = FALSE]
      B <- ranks[[p[2]]][ranks[[p[2]]]$Comparison == cmp, , drop = FALSE]
      same_ids <- nrow(A) == nrow(B) && setequal(A$feature_id, B$feature_id)
      addv("Identical PAS identifiers across PCA methods", cmp, same_ids, paste(p, collapse = " vs "))
      if (!same_ids) next
      B <- B[match(A$feature_id, B$feature_id), , drop = FALSE]
      kA <- comparison$Selected_k[comparison$Comparison == cmp & comparison$EVS_PCA_Method == p[1]]
      kB <- comparison$Selected_k[comparison$Comparison == cmp & comparison$EVS_PCA_Method == p[2]]
      le <- function(depth, k) if (length(k) == 1L && is.finite(k)) depth <= k else rep(NA, length(depth))
      leA <- le(A$Union_Depth, kA); leB <- le(B$Union_Depth, kB)
      top <- function(K, ra, rb) if (K <= nrow(A)) sum(ra <= K & rb <= K) / K else NA_real_
      base <- data.frame(Comparison = cmp, Method_A = p[1], Method_B = p[2], n_PAS = nrow(A),
                         Spearman_Rank_RT = stats::cor(A$Rank_RT, B$Rank_RT, method = "spearman"),
                         Spearman_Rank_ZT = stats::cor(A$Rank_ZT, B$Rank_ZT, method = "spearman"),
                         Spearman_Union_Depth = stats::cor(A$Union_Depth, B$Union_Depth, method = "spearman"),
                         k_A = if (length(kA)) kA else NA_integer_, k_B = if (length(kB)) kB else NA_integer_,
                         Leading_Edge_A = sum(leA), Leading_Edge_B = sum(leB), Leading_Edge_Shared = sum(leA & leB),
                         Leading_Edge_Jaccard = sum(leA & leB) / sum(leA | leB))
      for (K in c(500L, 1000L, 2500L, 5000L)) {
        base[[paste0("TopK", K, "_Overlap_RT")]] <- top(K, A$Rank_RT, B$Rank_RT)
        base[[paste0("TopK", K, "_Overlap_ZT")]] <- top(K, A$Rank_ZT, B$Rank_ZT)
        inA <- A$Union_Depth <= K; inB <- B$Union_Depth <= K
        base[[paste0("Union_k", K, "_Jaccard")]] <- if (any(inA | inB)) sum(inA & inB) / sum(inA | inB) else NA_real_
      }
      overlap[[length(overlap) + 1L]] <- base
    }
  }
  overlap <- do.call(rbind, overlap)
  if (!is.null(overlap)) utils::write.csv(overlap, file.path(out_dir, "EVS_PCA_Ranking_Overlap.csv"), row.names = FALSE)
  validation <- do.call(rbind, val)
  if (!is.null(validation)) utils::write.csv(validation, file.path(out_dir, "EVS_PCA_Method_Validation.csv"), row.names = FALSE)
  try_fig(write_regime_master_figures(multi_root, methods, comparison, out_dir), "master regime figures")
  invisible(list(comparison = comparison, overlap = overlap, validation = validation))
}

# Master figures: Figure S1 (rank regimes of every PCA method) and Figure 4
# (methods compared: k*, share of PASs on their own trend, Leading-Edge
# overlap, discoveries relative to the Original).
write_regime_master_figures <- function(multi_root, methods, comparison, out_dir) {
  suppressPackageStartupMessages({ library(ggplot2) })
  rd <- function(f) do.call(rbind, lapply(methods, function(pm) read_run_table(multi_root, pm, f)))
  curves <- rd("Table_Regime_D_Curves.csv"); knots <- rd("Table_Regime_Knots.csv")
  tsum <- rd("Table_TwoTrend_Cutoff_Summary.csv"); sites <- rd("Table_TwoTrend_Leading_Edge_Sites.csv")
  ord <- function(x) factor(x, levels = methods)
  if (!is.null(curves) && !is.null(knots)) {
    curves$EVS_PCA_Method <- ord(curves$EVS_PCA_Method); knots$EVS_PCA_Method <- ord(knots$EVS_PCA_Method)
    save_figure(plot_regime_definition(curves, knots) +
                  labs(title = "Figure S1. Remainder, Divergence and Leading-Edge rank regimes"),
                out_dir, "Figure_S1_Regime_Definition", 180, 95)
  }
  if (is.null(tsum)) return(invisible(FALSE))
  mcols <- stats::setNames(c("#D55E00", "#0072B2", "#7A7F86")[seq_along(methods)], methods)
  tsum$EVS_PCA_Method <- ord(tsum$EVS_PCA_Method)
  tsum$Comparison <- factor(tsum$Comparison, levels = COMPARISON_NAMES, labels = comp_lab())
  pa <- ggplot(tsum, aes(Comparison, k, fill = EVS_PCA_Method)) +
    geom_col(position = position_dodge(0.75), width = 0.68) +
    geom_text(aes(label = format(k, big.mark = ",")), position = position_dodge(0.75), vjust = -0.4, size = 1.9) +
    scale_fill_manual(values = mcols, name = NULL) +
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(c(0, 0.15))) +
    labs(x = NULL, y = "k* (PASs per arm)", title = "Two-trend cutoffs", tag = "a") + fig_theme()
  oc <- do.call(rbind, lapply(c("Two-trend k*" = "Expected_Correct_k", "Regime k*" = "Expected_Correct_Regime_k",
                                 "k = 5,000" = "Expected_Correct_k5000", "All Remainder" = "Expected_Correct_All_Remainder"),
                               function(col) data.frame(tsum[, c("Comparison", "EVS_PCA_Method")], v = tsum[[col]])))
  oc$Cutoff <- factor(rep(c("Two-trend k*", "Regime k*", "k = 5,000", "All Remainder"), each = nrow(tsum)),
                      levels = c("Two-trend k*", "Regime k*", "k = 5,000", "All Remainder"))
  pb <- ggplot(oc, aes(Cutoff, v, fill = Cutoff)) + geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", v)), vjust = -0.4, size = 1.6) +
    facet_grid(EVS_PCA_Method ~ Comparison) +
    scale_fill_manual(values = c("Two-trend k*" = "#CC79A7", "Regime k*" = "#D55E00", "k = 5,000" = "grey55",
                                 "All Remainder" = "grey80"), name = NULL) +
    coord_cartesian(ylim = c(0.5, 0.85)) +
    labs(x = NULL, y = "Expected share of PASs on their own trend", title = "Allocation accuracy", tag = "b") +
    fig_theme() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  rows <- list(list(plots = list(pa, pb), widths = c(1, 1.4)))
  pc <- NULL
  if (!is.null(sites) && length(methods) > 1L) {
    jac <- do.call(rbind, lapply(COMPARISON_NAMES, function(cmp) do.call(rbind, lapply(
      utils::combn(methods, 2L, simplify = FALSE), function(p) {
        A <- sites$feature_id[sites$Comparison == cmp & sites$EVS_PCA_Method == p[1]]
        B <- sites$feature_id[sites$Comparison == cmp & sites$EVS_PCA_Method == p[2]]
        data.frame(Comparison = cmp, Pair = paste(p, collapse = " vs "),
                   Jaccard = if (length(union(A, B))) length(intersect(A, B)) / length(union(A, B)) else NA_real_)
      }))))
    jac$Comparison <- factor(jac$Comparison, levels = COMPARISON_NAMES, labels = comp_lab())
    pc <- ggplot(jac, aes(Comparison, Jaccard, fill = Pair)) + geom_col(position = position_dodge(0.75), width = 0.68) +
      geom_text(aes(label = sprintf("%.2f", Jaccard)), position = position_dodge(0.75), vjust = -0.4, size = 1.9) +
      scale_y_continuous(limits = c(0, 1)) + scale_fill_grey(start = 0.35, end = 0.7, name = NULL) +
      labs(x = NULL, y = "Jaccard index", title = "Leading-Edge overlap between methods", tag = "c") + fig_theme()
  }
  pd <- NULL
  if (!is.null(comparison)) {
    dd <- do.call(rbind, lapply(c("Standard", "HBFSS"), function(cl) data.frame(
      Comparison = comparison$Comparison, EVS_PCA_Method = comparison$EVS_PCA_Method, Call = cl,
      Diff = comparison[[paste0("Diff_", cl)]])))
    dd$Comparison <- factor(dd$Comparison, levels = COMPARISON_NAMES, labels = comp_lab())
    dd$EVS_PCA_Method <- ord(dd$EVS_PCA_Method)
    pd <- ggplot(dd, aes(Comparison, Diff, fill = EVS_PCA_Method)) + geom_hline(yintercept = 0, colour = "grey30") +
      geom_col(position = position_dodge(0.75), width = 0.68) +
      geom_text(aes(label = sprintf("%+d", as.integer(Diff)), vjust = ifelse(Diff >= 0, -0.4, 1.3)),
                position = position_dodge(0.75), size = 1.8) +
      facet_wrap(~ Call, nrow = 1, scales = "free_y") + scale_fill_manual(values = mcols, name = NULL) +
      scale_y_continuous(expand = expansion(c(0.15, 0.15))) +
      labs(x = NULL, y = "EVS (LE + REM) minus Original", title = "Discoveries relative to the Original", tag = "d") +
      fig_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }
  second <- Filter(Negate(is.null), list(pc, pd))
  if (length(second)) rows[[2]] <- list(plots = second, widths = if (length(second) == 2L) c(1, 1.3) else 1)
  save_figure(assemble_rows(rows, heights = rep(1, length(rows)),
                            title = "Figure 4. PCA preprocessing methods compared under the two-trend cutoff"),
              out_dir, "Figure_4_PCA_Methods_Summary", 180, 150)
  invisible(TRUE)
}

# =============================================================================
# MASTER PROCESS: one analysis run per PCA preprocessing method, the
# cross-method comparison, then one ZIP of everything
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

  # The cutoff k is chosen inside the run, per comparison, before any test.
  cutoff_manifest <- data.frame(
    Cutoff_Method = unname(METHOD_SLUGS),
    RT0_ZT6 = NA_integer_, RT2_ZT8 = NA_integer_, RT4_ZT10 = NA_integer_, RT8_ZT14 = NA_integer_,
    Value = "k per condition, chosen in the run",
    Source = "Two-trend rule: two mean-dispersion trends fitted by EM on the 30 samples outside each comparison; k* maximizes the agreement of the PC1 split with them; no differential-expression test")
  utils::write.csv(cutoff_manifest, file.path(multi_root, "Cutoff_Method_Manifest.csv"), row.names = FALSE)

  pca_methods <- unique(toupper(trimws(strsplit(Sys.getenv("SEQUENCE_PCA_METHODS",
                                                           unset = paste(PCA_METHOD_DEFAULT, collapse = ",")), ",")[[1]])))
  if (!length(pca_methods) || !all(pca_methods %in% PCA_METHOD_LEVELS)) {
    stop("SEQUENCE_PCA_METHODS must be a comma-separated subset of: ", paste(PCA_METHOD_LEVELS, collapse = ", "))
  }
  # Inputs shared by every method's run (raw counts, size factors, gene-wise
  # dispersions); written by the first run, reused and checked by the others.
  shared_cache <- file.path(multi_root, "Shared_Cache")
  dir.create(shared_cache, showWarnings = FALSE)

  # One BLAS / OpenMP thread per run (unless set by the user).
  thread_vars <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
  thread_env <- paste0(thread_vars[!nzchar(Sys.getenv(thread_vars))], "=1")
  for (pm in pca_methods) {
    message("=====================================================")
    message("Starting the analysis run: EVS PCA preprocessing ", pm)
    message("=====================================================")
    status <- system2(
      command = file.path(R.home("bin"), "Rscript"), args = shQuote(this_script),
      env = c("SEQUENCE_CHILD_RUN=1", paste0("SEQUENCE_CUTOFF_METHOD=", METHOD_KEYS[[1]]),
              paste0("EVS_PCA_METHOD=", pm),
              paste0("SEQUENCE_MULTI_ROOT=", multi_root), paste0("SEQUENCE_REPO_ROOT=", repo_guess), thread_env))
    if (!identical(as.integer(status), 0L)) stop("Analysis run failed for ", pm, " (exit status ", status, ")")
  }
  unlink(shared_cache, recursive = TRUE, force = TRUE)
  cross <- tryCatch(write_pca_method_comparison(multi_root, pca_methods),
                    error = function(e) { warning("PCA method comparison skipped: ", conditionMessage(e), call. = FALSE); NULL })

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
  cat("Analysis complete.\n")
  if (!is.null(cross$comparison)) {
    cc <- cross$comparison
    for (i in seq_len(nrow(cc)))
      cat(sprintf("  %-9s %-10s k* = %5d  LE %5d  EVS Std %5d (%+5d)  HBFSS %5d (%+5d)\n",
                  cc$Comparison[i], cc$EVS_PCA_Method[i], cc$Selected_k[i], cc$Leading_Edge[i],
                  cc$EVS_Standard[i], cc$Diff_Standard[i], cc$EVS_HBFSS[i], cc$Diff_HBFSS[i]))
  }
  for (pm in pca_methods) {
    chk <- read_run_table(multi_root, pm, "Table_Validation_Checks.csv")
    if (is.null(chk)) next
    cat(sprintf("  %-10s validation checks: %d PASS, %d FAIL\n", pm, sum(chk$Status == "PASS"), sum(chk$Status == "FAIL")))
    bad <- chk[chk$Status == "FAIL", , drop = FALSE]
    for (i in seq_len(nrow(bad))) cat("    FAIL:", bad$Comparison[i], "-", bad$Check[i], "-", bad$Detail[i], "\n")
  }
  if (!is.null(cross$validation)) {
    v <- cross$validation
    cat(sprintf("  Cross-method checks: %d PASS, %d FAIL\n", sum(v$Status == "PASS"), sum(v$Status == "FAIL")))
    for (i in which(v$Status == "FAIL")) cat("    FAIL:", v$Comparison[i], "-", v$Check[i], "-", v$Detail[i], "\n")
  }
  cat("Final ZIP:\n", normalizePath(final_zip, winslash = "/", mustWork = TRUE), "\n", sep = "")
  cat("=====================================================\n\n")
  quit(save = "no", status = 0L)
}

# Active cutoff method of this child run.
CUTOFF_METHOD_KEY <- Sys.getenv("SEQUENCE_CUTOFF_METHOD", unset = "twotrend")
CUTOFF_METHOD_LABEL <- names(METHOD_KEYS)[match(CUTOFF_METHOD_KEY, METHOD_KEYS)]
if (is.na(CUTOFF_METHOD_LABEL)) stop("Unknown SEQUENCE_CUTOFF_METHOD: ", CUTOFF_METHOD_KEY)
CUTOFF_METHOD_BASIS <- paste0("EVS PCA preprocessing ", EVS_PCA_METHOD, ". Per comparison, two mean-dispersion trends fitted by EM to DESeq2 gene-wise dispersions of the 30 samples outside the comparison; k* maximizes the allocation score of the PC1 split against them (regime k* reported as reference); no differential-expression test used")


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
# score yields no HC threshold and therefore no HBFSS discoveries in that
# analysis. HC is computed once per analysis, over the empirical p-values of all
# its views (call_significance).
HC_ALPHA0 <- 0.10

if (!is.numeric(HC_ALPHA0) || length(HC_ALPHA0) != 1L ||
    !is.finite(HC_ALPHA0) || HC_ALPHA0 <= 0 || HC_ALPHA0 > 1) {
  stop("HC_ALPHA0 must be a single finite number in (0, 1].")
}

# EVS selection size k* is comparison-specific, chosen before any test by the
# regime rule (REGIME CUTOFF k*).
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

track_short <- c(normalized_evs = "NormEVS")

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

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix   = c("R0", "R2", "R4", "R8"),
  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

# k* per comparison, set by compute_twotrend_cutoffs() before any comparison runs.
OPTIMAL_K <- list()
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(as.integer(Sys.getenv("OPENBLAS_NUM_THREADS", unset = "1")))
  RhpcBLASctl::omp_set_num_threads(as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1")))
}
# Run-wide validation checks, and the significance-call counter used to show
# that the cutoff rule performs no differential-expression test.
# The Original's Cook's outlier PASs per comparison: applied in every view.
original_cooks_cache <- list()
# The Original refitted with DESeq2's default parametric trend (reference only).
default_trend_rows <- list()
validation_rows <- list()
pca_ranking_rows <- list()   # per comparison: PC1 ranks and union depth of this run's PCA method
# Per comparison: the raw matrix, the regime ranks and union depth, shared size
# factors and gene-wise dispersions (prepare_comparison_cache), reused by the
# analysis.
comparison_cache <- list()
# wald: DESeq2::nbinomWaldTest fits; null: call_significance (fdrtool / HC / HBFSS)
SIGNIFICANCE_CALLS <- list2env(list(wald = 0L, null = 0L))

get_active_evs_cutoff <- function(comparison_name) {
  idx <- match(as.character(comparison_name), comparison_table$comparison_name)
  if (is.na(idx)) stop("No EVS cutoff is defined for comparison: ", comparison_name)
  k <- OPTIMAL_K[[as.character(comparison_name)]] %||%
    stop("The cutoff k has not been chosen for ", comparison_name)
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || !is.finite(k) || k < 1L) {
    stop("Invalid EVS cutoff for ", comparison_name, " under ", CUTOFF_METHOD_LABEL)
  }
  k
}


if (anyDuplicated(comparison_table$comparison_name)) stop("comparison_table contains duplicated comparison names.")

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
output_dir <- file.path(multi_output_root, pca_run_folder(EVS_PCA_METHOD))
SHARED_CACHE_DIR <- if (nzchar(Sys.getenv("SEQUENCE_MULTI_ROOT", unset = ""))) file.path(multi_output_root, "Shared_Cache") else ""

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

# Null for the HBFSS empirical p-values (SEQUENCE_NULL_METHOD):
#   fndr         fdrtool's censored empirical null, central region chosen by the
#                false-nondiscovery-rate rule (fdrtool default; the paper's choice)
#   pct0         fdrtool's censored empirical null fitted to the central
#                EMPIRICAL_NULL_PCT0 fraction of the statistics (narrower centre,
#                less affected by many small real effects)
#   theoretical  N(0, 1), no empirical correction (reference only)
EMPIRICAL_NULL_METHOD <- Sys.getenv("SEQUENCE_NULL_METHOD", unset = "fndr")
EMPIRICAL_NULL_PCT0 <- 0.5
if (!EMPIRICAL_NULL_METHOD %in% c("fndr", "pct0", "theoretical")) {
  stop("SEQUENCE_NULL_METHOD must be fndr, pct0 or theoretical.")
}

run_empirical_null_fdrtool <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  stat_vec <- unname(stat_vec)

  if (length(stat_vec) < 5L) {
    stop(sprintf("[%s] Fewer than 5 finite Wald statistics were available for fdrtool.", dataset_name))
  }

  if (identical(EMPIRICAL_NULL_METHOD, "theoretical")) {
    return(list(pval = clip_probabilities(2 * stats::pnorm(-abs(stat_vec))),
                param = matrix(c(1, NA_real_), nrow = 1, dimnames = list(NULL, c("sd", "eta0")))))
  }
  fit <- fdrtool::fdrtool(
    stat_vec,
    statistic = "normal",
    plot = FALSE,
    verbose = FALSE,
    cutoff.method = EMPIRICAL_NULL_METHOD,
    pct0 = EMPIRICAL_NULL_PCT0
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
  # asserted for the analysis.
  if (!is.finite(best_hc) || is.na(best_hc) || best_hc <= 0) {
    message(
      sprintf(
        "[%s] No positive Higher-Criticism signal within the lowest %.1f%% of empirical p-values; HBFSS disabled for this analysis.",
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
  # per-comparison three-view TWAS panels already contain the same count summary.
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
                "Table_EVS_Composition.csv", "Table_TwoTrend_Cutoff_Summary.csv", "Table_TwoTrend_Heldout_Fits.csv",
                "Table_Regime_Cutoff_Summary.csv", "Table_Regime_Knots.csv",
                "Table_Validation_Checks.csv", "Table_Original_DESeq2_Default_Trend.csv")) {
    add_if_exists(file.path(summary_table_dir, tab))
  }
  add_if_exists(file.path(summary_table_dir, "Table_EVS_Split_Audit.csv"))
  add_if_exists(file.path(summary_table_dir, "Table_EVS_PCA_Evidence.csv"))
  add_if_exists(file.path(summary_table_dir, "Table_EVS_PCA_Ranking.csv"))

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

# Per comparison and PCA method: the comparison's slice of the experiment-wide
# PCA input and the per-condition PC1 fits, cached (keyed "<comparison>|<method>").
PCA_RANK_CACHE <- list()

# The matrix the EVS PCA is run on: the experiment-wide matrix of the PCA
# method (REGIME CUTOFF k*) restricted to the comparison's PASs and samples,
# the same values the regime rule ranked. Used only to rank PASs; no
# transformed matrix is ever passed to DESeq2. (size_factors is unused: the
# ranking uses experiment-wide normalization.)
make_rank_matrix_for_track <- function(count_matrix, coldata, track_key, method = EVS_PCA_METHOD,
                                       size_factors = NULL, cache_key = NULL) {
  track_key <- match.arg(track_key, "normalized_evs")
  method <- match.arg(method, PCA_METHOD_LEVELS)
  key <- if (is.null(cache_key)) NULL else paste(cache_key, method, sep = "|")
  hit <- if (is.null(key)) NULL else PCA_RANK_CACHE[[key]]
  if (!is.null(hit) && identical(hit$ids, rownames(count_matrix)) && identical(hit$samples, colnames(count_matrix))) {
    return(hit$rank_obj)
  }
  full <- experiment_rank_matrix(method)
  miss <- c(setdiff(rownames(count_matrix), rownames(full)), setdiff(colnames(count_matrix), colnames(full)))
  if (length(miss)) stop("The experiment-wide ", method, " matrix lacks comparison PASs or samples.")
  x <- full[rownames(count_matrix), colnames(count_matrix), drop = FALSE]
  rank_obj <- list(
    rank_matrix = as.data.frame(x),
    method = method,
    preprocessing_label = paste0(PCA_METHOD_LABELS[[method]], " prior to EVS; ",
                                 "PCA is performed directly on this matrix (no further transformation)")
  )
  if (!is.null(key)) PCA_RANK_CACHE[[key]] <<- list(ids = rownames(count_matrix), samples = colnames(count_matrix),
                                                   rank_obj = rank_obj, fits = NULL)
  rank_obj
}

# Per-condition PC1 fits with the full ranking (cached with the matrix).
pc1_condition_fits <- function(rank_obj, coldata, cache_key = NULL) {
  key <- if (is.null(cache_key)) NULL else paste(cache_key, rank_obj$method, sep = "|")
  if (!is.null(key) && !is.null(PCA_RANK_CACHE[[key]]$fits)) return(PCA_RANK_CACHE[[key]]$fits)
  ids <- colnames(rank_obj$rank_matrix)
  fits <- lapply(c(trt = "trt", untrt = "untrt"), function(cond)
    compute_pc1_loading_table(rank_obj$rank_matrix, ids[coldata$condition == cond], top_n = 1L,
                              preprocessing_label = rank_obj$preprocessing_label))
  if (!is.null(key) && !is.null(PCA_RANK_CACHE[[key]])) PCA_RANK_CACHE[[key]]$fits <<- fits
  fits
}

# A cached PC1 fit at cutoff k: identical to compute_pc1_loading_table(top_n = k).
pc1_fit_at_k <- function(fit, k) {
  lt <- fit$loading_table
  ci <- resolve_top_n_cutoff(lt$pc1_loading_abs, top_n = k)
  lt$split_class <- ifelse(lt$rank <= ci$top_n_actual, "high_loading", "background_loading")
  fit$loading_table <- lt
  fit$cutoff <- ci$cutoff_value
  fit$top_n_used <- ci$top_n_actual
  fit
}

compute_pc1_loading_table <- function(value_df, sample_names, top_n, preprocessing_label = "Normalized prior to EVS") {
  x <- as.matrix(value_df[, sample_names, drop = FALSE])
  storage.mode(x) <- "numeric"

  if (ncol(x) < 2L) {
    stop("EVS PCA requires at least two samples in each condition.")
  }

  # Refined EVS method: run PCA directly on the pre-split matrix for the
  # condition, extract the PC1 feature eigenvector/loading, take absolute values,
  # rank descending, and select the comparison-specific k. The matrix is used
  # as supplied (make_rank_matrix_for_track); nothing is log-transformed again.
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
  track_key <- match.arg(track_key, "normalized_evs")
  empirical_k <- get_active_evs_cutoff(comparison_name)

  # The EVS matrix is used only to choose feature IDs. NormEVS
  # follows the supplied workflow: median-of-ratios normalization of the full
  # RT/ZT comparison matrix, followed by condition-specific PCA/PC1 eigenvectors.
  # Downstream DESeq2 always receives the corresponding raw-count subset.
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

  # The comparison slice of the experiment-wide PCA input and its PC1 rankings (cached).
  rank_obj <- make_rank_matrix_for_track(
    count_matrix = raw_count_matrix,
    coldata = coldata,
    track_key = track_key,
    size_factors = comparison_cache[[comparison_name]]$sf,
    cache_key = comparison_name
  )

  rank_matrix <- rank_obj$rank_matrix
  preprocessing_label <- rank_obj$preprocessing_label

  fits <- pc1_condition_fits(rank_obj, coldata, comparison_name)
  fit_trt <- pc1_fit_at_k(fits$trt, empirical_k)
  fit_untrt <- pc1_fit_at_k(fits$untrt, empirical_k)

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
    pca_method = rank_obj$method,
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
# Every DESeq2 fit in this pipeline (Original, Leading Edge, Remainder,
# calibration) therefore uses
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
# Sampling-bias correction. The log-scale fit estimates the average LOG
# gene-wise dispersion at each mean, and a log gene-wise estimate is biased low:
# with nu = m - p residual degrees of freedom (m samples, p coefficients),
# alpha_hat ~ alpha * chi^2_nu / nu, the same sampling model DESeq2 uses for the
# prior variance (trigamma(nu/2)), so
#     E[log alpha_hat] = log alpha + psi(nu/2) - log(nu/2)
# (psi = digamma; -0.130 for nu = 8). The prior must be centred on the true log
# dispersion, so a0 and a1 are multiplied by exp(-(psi(nu/2) - log(nu/2)))
# (1.139 for nu = 8). b and the prior variance (a squared MAD of log
# residuals, which a constant shift does not change) are unaffected.
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
  fine <- fine[fine >= TREND_B_RANGE[1] - 1e-12 & fine <= TREND_B_RANGE[2] + 1e-12]   # stay within TREND_B_RANGE
  best <- best_of(c(list(best), fits_for(fine)))
  trend <- make_trend_function(best[["a0"]], best[["a1"]], best[["b"]])
  # Fit diagnostics: b at the edge of its range, or a0 / a1 at a bound, mean
  # the trend shape is constrained by the settings rather than the data.
  attr(trend, "b_at_range_edge") <- min(abs(best[["b"]] - TREND_B_RANGE)) < TREND_B_FINE_STEP / 2
  attr(trend, "coef_at_bound") <- any(c(best[["a0"]], best[["a1"]]) <= TREND_MIN_COEF * 1.0001) ||
    best[["a0"]] >= TREND_MAX_COEF[["a0"]] * 0.9999 || best[["a1"]] >= TREND_MAX_COEF[["a1"]] * 0.9999
  attr(trend, "n_exponents_fitted") <- length(fits)
  trend
}

# Fit the trend to a DESeqDataSet that already has gene-wise dispersions and
# store it with DESeq2's dispersionFunction() setter (which also fills dispFit,
# the residual variance and the dispersion prior variance).
set_dispersion_trend <- function(dds) {
  mc <- S4Vectors::mcols(dds)
  if (any(mc$allZero)) stop("remove PASs with no counts before fitting the dispersion trend.")
  use <- !mc$allZero & is.finite(mc$dispGeneEst) & mc$dispGeneEst >= TREND_MIN_GENE_DISP
  log_trend <- fit_dispersion_trend(mc$baseMean[use], mc$dispGeneEst[use])
  log_fit <- attr(log_trend, "coefficients")
  X <- stats::model.matrix(DESeq2::design(dds), as.data.frame(SummarizedExperiment::colData(dds)))
  nu <- ncol(dds) - ncol(X)
  if (nu < 1) stop("the dispersion trend needs more samples than model coefficients.")
  log_bias <- digamma(nu / 2) - log(nu / 2)
  trend <- make_trend_function(log_fit[["a0"]] * exp(-log_bias), log_fit[["a1"]] * exp(-log_bias), log_fit[["b"]])
  attr(trend, "log_sampling_bias") <- log_bias
  for (a in c("b_at_range_edge", "coef_at_bound", "n_exponents_fitted")) attr(trend, a) <- attr(log_trend, a)
  DESeq2::`dispersionFunction<-`(dds, value = trend)
}

# DESeq2's dispersion steps with the trend above: size factors, gene-wise
# dispersions, trend, MAP shrinkage. These are the steps DESeq() runs, with the
# trend fit replaced. Size factors: the comparison's shared factors, estimated
# once on the full comparison matrix, so a PAS has identical normalized counts
# in the Original and in its EVS partition. Supplied factors are never
# re-estimated downstream (NULL, used only for a whole comparison matrix,
# estimates them from that matrix).
apply_size_factors <- function(dds, size_factors = NULL) {
  if (is.null(size_factors)) return(DESeq2::estimateSizeFactors(dds))
  if (!setequal(names(size_factors), colnames(dds)) || anyNA(size_factors) || any(size_factors <= 0)) {
    stop("Shared size factors do not match the samples of this dataset.")
  }
  sizeFactors(dds) <- unname(size_factors[colnames(dds)])
  dds
}

# Shared size factors of a comparison: DESeq2 median-of-ratios on the full
# comparison matrix, used by the Original and every EVS view.
comparison_size_factors <- function(raw, coldata) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    coerce_raw_count_matrix_for_deseq2(raw, context = "shared size factors"),
    colData = coldata, design = make_design_formula(coldata))
  sf <- DESeq2::sizeFactors(DESeq2::estimateSizeFactors(dds))
  stats::setNames(as.numeric(sf), colnames(raw))
}

# Gene-wise dispersions computed once per comparison.
# A PAS's gene-wise (Cox-Reid) dispersion depends only on its own counts, the
# design and the size factors, so with shared size factors it is the same in
# the Original and in any EVS partition: it is estimated once on the full
# comparison matrix and reused for every split. Only the trend, the prior and
# the MAP dispersions depend on the split. Every run checks the equality on its
# final analysis (validate_track_normalization).
comparison_gene_dds <- function(raw, coldata, size_factors) {
  if (is.null(size_factors)) return(NULL)
  dds <- DESeq2::DESeqDataSetFromMatrix(coerce_raw_count_matrix_for_deseq2(raw, context = "gene-wise dispersions"),
                                        colData = coldata, design = make_design_formula(coldata))
  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
  dds <- apply_size_factors(dds, size_factors)
  suppressMessages(DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE))
}

estimate_dispersions_with_trend <- function(dds, size_factors = NULL, gene_dds = NULL) {
  if (!is.null(gene_dds)) {
    if (!all(rownames(dds) %in% rownames(gene_dds))) stop("Precomputed gene-wise dispersions lack some PASs.")
    dds <- gene_dds[rownames(dds), ]
  } else {
    dds <- apply_size_factors(dds, size_factors)
    dds <- suppressMessages(DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE))
  }
  dds <- set_dispersion_trend(dds)
  suppressMessages(DESeq2::estimateDispersionsMAP(dds, quiet = TRUE))
}

# =============================================================================
# DIFFERENTIAL EXPRESSION: MODEL STAGE (per analysis view) AND SIGNIFICANCE
# STAGE (per analysis view)
# =============================================================================
# Each analysis view (Original; Leading Edge; Remainder) is its own DESeq2
# dataset with its own dispersion trend, prior and shrinkage (size factors are
# shared from the full comparison matrix)
# (model stage, fit_view_model), and its own significance calibration
# (significance stage, call_significance): fdrtool fits an empirical null to
# the view's own Wald statistics, and higher criticism sets one HBFSS threshold
# per comparison from the Original's empirical p-values, applied to every view. The Leading Edge and Remainder are
# different populations (different expression ranges, dispersion laws and
# fractions of changed PASs), and a null shared by the two would be too strict
# for one and too loose for the other (Efron 2008, "Simultaneous inference:
# when should hypothesis testing problems be combined?").

# Model stage: DESeq2 with the power trend, the three Wald tests, apeglm.
fit_view_model <- function(count_mat, coldata, dataset_name, size_factors = NULL, gene_dds = NULL,
                           cooks_outliers = NULL, trend = c("power", "deseq2_default")) {
  trend <- match.arg(trend)
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
  # trend = "deseq2_default" (reference baseline only): DESeq2's own
  # parametric gamma-GLM trend, with DESeq2's automatic fall-back to a local fit.
  dds <- if (identical(trend, "power")) estimate_dispersions_with_trend(dds, size_factors, gene_dds) else
    suppressMessages(DESeq2::estimateDispersions(apply_size_factors(dds, size_factors), fitType = "parametric", quiet = TRUE))
  dds <- DESeq2::nbinomWaldTest(dds, quiet = TRUE)
  SIGNIFICANCE_CALLS$wald <- SIGNIFICANCE_CALLS$wald + 1L

  # The three DESeq2 tests with independentFiltering = FALSE and
  # cooksCutoff = FALSE. The Cook's outlier rule is applied below with ONE
  # outlier list per comparison (the Original's), and the BH correction is
  # applied later, once per analysis (call_significance), so the Original and
  # the EVS views test the same PASs under the same rules.
  contrast <- c("condition", "trt", "untrt")
  tests <- list(std = list(alpha = BH_FDR_STANDARD),
                strong = list(lfcThreshold = lfc_boundary, altHypothesis = "greaterAbs", alpha = BH_FDR_STRONG),
                weak = list(lfcThreshold = lfc_boundary, altHypothesis = "lessAbs", alpha = BH_FDR_WEAK))
  prim <- lapply(tests, function(a) do.call(DESeq2::results, c(list(dds, contrast = contrast, pAdjustMethod = "BH",
                                                                    independentFiltering = FALSE, cooksCutoff = FALSE), a)))

  res_df <- as.data.frame(prim$std)
  res_df$feature_id <- as.character(rownames(res_df))
  # Cook's outliers by DESeq2's default rule (Cook's distance above the 0.99
  # quantile of F(p, m - p)), detected in the Original's fit (cooks_outliers =
  # NULL) and supplied as that list to every EVS view of the comparison.
  if (is.null(cooks_outliers)) {
    cooks_p <- DESeq2::results(dds, contrast = contrast, independentFiltering = FALSE)$pvalue
    res_df$cooks_outlier <- is.na(cooks_p) & !is.na(prim$std$pvalue)
  } else {
    res_df$cooks_outlier <- res_df$feature_id %in% cooks_outliers
  }
  bh_view <- function(p) { q <- rep(NA_real_, length(p)); ok <- is.finite(p); q[ok] <- stats::p.adjust(p[ok], "BH"); q }
  res_df$pvalue[res_df$cooks_outlier] <- NA_real_
  res_df$padj <- bh_view(res_df$pvalue)      # replaced by the analysis-wide BH in call_significance

  shr <- DESeq2::lfcShrink(dds, coef = get_condition_coef(dds), type = "apeglm")
  shr_df <- data.frame(feature_id = rownames(shr), lfc_shrunk = shr$log2FoldChange)
  res_df <- dplyr::left_join(res_df, shr_df, by = "feature_id")

  res_df$regulation_direction <- ifelse(
    is.na(res_df$lfc_shrunk), NA_character_,
    ifelse(res_df$lfc_shrunk > 0, "upregulated",
           ifelse(res_df$lfc_shrunk < 0, "downregulated", "no_change")))

  i_ga <- match(res_df$feature_id, rownames(prim$strong))
  i_la <- match(res_df$feature_id, rownames(prim$weak))
  res_df$resGA_pvalue <- ifelse(res_df$cooks_outlier, NA_real_, prim$strong$pvalue[i_ga])
  res_df$resLA_pvalue <- ifelse(res_df$cooks_outlier, NA_real_, prim$weak$pvalue[i_la])
  res_df$resGA_padj <- bh_view(res_df$resGA_pvalue)
  res_df$resLA_padj <- bh_view(res_df$resLA_pvalue)

  list(dds = dds, res_df = res_df, dataset_name = dataset_name)
}

# Significance stage for one analysis (a named list of model-stage fits):
#   1. one fdrtool empirical null per view, fitted to that view's own Wald
#      statistics, so each view's empirical p-values are calibrated to its own
#      population (the Leading Edge and the Remainder differ);
#   2. one higher-criticism threshold per analysis, from the empirical p-values
#      of all its views together (the Original alone; the Leading Edge and the
#      Remainder combined), as the BH correction is one per analysis. Every
#      view of an analysis therefore has the same HCp and Htau, and the Original
#      and the split each get one threshold from one family of all PASs;
#   3. each view's HBFSS scores and calls.
# Returns the per-view results and one summary row per view.
# hc_p_fixed: when given (the EVS analysis), the HC threshold of the
# comparison's Original is used instead of a new one, so the Original and the
# Leading Edge + Remainder share one HBFSS threshold per comparison.
call_significance <- function(fits, analysis_label, annot_df, lean = FALSE, hc_p_fixed = NULL) {
  SIGNIFICANCE_CALLS$null <- SIGNIFICANCE_CALLS$null + 1L
  fits <- combine_bh_across_views(fits)
  null_rows <- list()
  for (i in seq_along(fits)) {
    f <- fits[[i]]
    view_label <- paste(analysis_label, names(fits)[i])
    valid <- is.finite(f$res_df$stat)
    if (sum(valid) < 5L) {
      stop(sprintf("[%s] Fewer than 5 finite Wald statistics were available for fdrtool.", view_label))
    }
    fdr_fit <- run_empirical_null_fdrtool(f$res_df$stat[valid], dataset_name = view_label)
    if (length(fdr_fit$pval) != sum(valid)) {
      stop(sprintf("[%s] fdrtool empirical-p output length mismatch.", view_label))
    }
    param <- fdr_fit$param
    null_rows[[i]] <- data.frame(
      View = names(fits)[i],
      n_wald_statistics = sum(valid),
      null_sd = if ("sd" %in% colnames(param)) unname(param[1, "sd"]) else NA_real_,
      null_fraction_eta0 = if ("eta0" %in% colnames(param)) unname(param[1, "eta0"]) else NA_real_)
    fits[[i]]$res_df$empirical_p <- NA_real_
    fits[[i]]$res_df$empirical_p[valid] <- as.numeric(fdr_fit$pval)
  }
  all_p <- unlist(lapply(fits, function(f) f$res_df$empirical_p), use.names = FALSE)
  all_p <- all_p[is.finite(all_p)]
  hc_source <- if (is.null(hc_p_fixed)) "this analysis" else "the comparison's Original"
  hc_p <- if (is.null(hc_p_fixed)) safe_hc_thresh(all_p, dataset_name = paste(analysis_label, "(all views)")) else hc_p_fixed
  hc_ok <- !is.na(hc_p) && is.finite(hc_p) && hc_p > 0 && hc_p <= 1
  htau <- if (hc_ok) abs(log10(hc_p)) * lfc_boundary else NA_real_
  null_rows <- do.call(rbind, null_rows)
  null_rows$hc_p_threshold <- hc_p
  null_rows$hc_threshold_source <- hc_source
  null_rows$hbfss_threshold <- htau
  null_rows$n_hc_empirical_p <- length(all_p)
  views <- list()
  for (i in seq_along(fits)) views[[names(fits)[i]]] <- classify_view(fits[[i]], hc_p, htau, annot_df, lean = lean)
  list(views = views, null = null_rows)
}

# One BH correction per analysis for each DESeq2 test: over the Original's raw
# p-values, or over the Leading Edge and Remainder raw p-values combined (one
# hypothesis family; the views are disjoint). No independent filtering.
combine_bh_across_views <- function(fits) {
  n <- vapply(fits, function(f) nrow(f$res_df), integer(1))
  idx <- rep(seq_along(fits), n)
  for (col in c("pvalue", "resGA_pvalue", "resLA_pvalue")) {
    p <- unlist(lapply(fits, function(f) f$res_df[[col]]), use.names = FALSE)
    q <- rep(NA_real_, length(p))
    ok <- is.finite(p)
    q[ok] <- stats::p.adjust(p[ok], method = "BH")
    for (i in seq_along(fits)) fits[[i]]$res_df[[sub("pvalue", "padj", col)]] <- q[idx == i]
  }
  fits
}

# Reference baseline: the Original with DESeq2's default parametric trend
# (gamma GLM, b = 1, DESeq2's automatic local fall-back), all other steps and
# testing rules unchanged, next to the Original with the power trend. Shows
# how much of any EVS gain over standard DESeq2 practice comes from the trend.
record_default_trend_baseline <- function(comparison_name, raw, coldata, annot_df, sf, power_out) {
  fit <- fit_view_model(raw, coldata, paste(comparison_name, "Original (DESeq2 default trend)"),
                        size_factors = sf, trend = "deseq2_default")
  sig <- call_significance(list(raw_dataset = fit), paste(comparison_name, "Original (DESeq2 default trend)"),
                           annot_df, lean = TRUE)
  row <- function(label, res, dds, nl, fit_type) {
    cnt <- count_discoveries(discovery_flags(list(raw_dataset = res)))
    f <- DESeq2::dispersionFunction(dds)
    X <- stats::model.matrix(DESeq2::design(dds), as.data.frame(SummarizedExperiment::colData(dds)))
    data.frame(Comparison = comparison_name, Analysis = label, Trend = fit_type,
               Median_Final_Dispersion = stats::median(S4Vectors::mcols(dds)$dispersion, na.rm = TRUE),
               Share_PAS_Below_Trend = mean(S4Vectors::mcols(dds)$dispGeneEst < S4Vectors::mcols(dds)$dispFit, na.rm = TRUE),
               Prior_Variance = as.numeric(attr(f, "varLogDispEsts")) - trigamma((ncol(dds) - ncol(X)) / 2),
               Null_SD = nl$null_sd[1], Htau = nl$hbfss_threshold[1],
               Standard = cnt[["Standard"]], Strong = cnt[["Strong"]], Weak = cnt[["Weak"]], HBFSS = cnt[["HBFSS"]])
  }
  ft <- attr(DESeq2::dispersionFunction(sig$views$raw_dataset$dds), "fitType")
  default_trend_rows[[comparison_name]] <<- rbind(
    row(ORIGINAL_LABEL, power_out$views$raw_dataset$results, power_out$views$raw_dataset$dds,
        power_out$null, "power trend (this pipeline)"),
    row("Original (DESeq2 default trend)", sig$views$raw_dataset$results, sig$views$raw_dataset$dds, sig$null,
        paste("DESeq2", if (is.null(ft)) "parametric" else ft)))
  invisible(TRUE)
}

# The Original's Cook's outlier list of a comparison (fitted once; the EVS views
# of that comparison use it instead of their own).
remember_original_cooks <- function(comparison_name, results) {
  original_cooks_cache[[comparison_name]] <<- as.character(results$feature_id[results$cooks_outlier %in% TRUE])
  invisible(original_cooks_cache[[comparison_name]])
}
original_cooks_outliers <- function(comparison_name, raw, coldata, size_factors) {
  if (!is.null(original_cooks_cache[[comparison_name]])) return(original_cooks_cache[[comparison_name]])
  fit <- fit_view_model(raw, coldata, paste(comparison_name, "Original (Cook's outliers)"), size_factors = size_factors)
  remember_original_cooks(comparison_name, fit$res_df)
}

run_core_analysis <- function(count_mats, coldata, dataset_names, annot_df, analysis_label, size_factors = NULL,
                              gene_dds = NULL, cooks_outliers = NULL, hc_p_fixed = NULL) {
  fits <- Map(function(m, nm) fit_view_model(m, coldata, nm, size_factors = size_factors, gene_dds = gene_dds,
                                             cooks_outliers = cooks_outliers), count_mats, dataset_names)
  names(fits) <- names(count_mats)
  call_significance(fits, analysis_label, annot_df, hc_p_fixed = hc_p_fixed)
}

# HBFSS scores and the Standard / Strong / Weak / HBFSS calls for one view,
# given its empirical p-values and the analysis's HC threshold.
# lean = TRUE (used only by the default-trend baseline) returns the decision table
# without normalized counts, annotation and dispersion columns; every flag and
# runtime check is identical. The final analysis always uses lean = FALSE.
classify_view <- function(fit, hc_p_threshold, hbfss_threshold, annot_df, lean = FALSE) {
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
  if (lean) {
    res_df$dataset_name <- dataset_name
    res_df$hc_p_threshold_dataset <- hc_p_threshold
    res_df$hbfss_threshold_dataset <- hbfss_threshold
    return(list(dds = dds, results = res_df))
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

  analysis_dataset_list <- dataset_list

  # Two analyses: the Original, then the Leading Edge with the Remainder. Each
  # view gets its own empirical null; each analysis one BH. One HC / HBFSS
  # threshold per comparison: set by the Original and reused by the EVS
  # analysis (call_significance, hc_p_fixed).
  analyses <- list(original = "raw_dataset", split = c("leading_edge_dataset", "remainder_dataset"))

  # Shared size factors and gene-wise dispersions: estimated once on the full
  # comparison matrix (prepare_comparison_cache) and reused by every view.
  cache <- comparison_cache[[comparison_name]]
  if (is.null(cache) || !identical(rownames(cache$raw), rownames(evs$raw_dataset)))
    stop(comparison_name, ": the comparison cache does not match the analysed matrix")
  sf <- cache$sf

  analysis_results <- list()
  original_hc_p <- NULL
  for (a in names(analyses)) {
    keys <- analyses[[a]]
    analysis_label <- if (a == "original") ORIGINAL_LABEL else paste(unname(track_short[track_key]), "LE + REM")
    # Cook's outliers: the Original uses its own; the EVS views use the Original's.
    cooks <- if (a == "original") NULL else original_cooks_outliers(comparison_name, evs$raw_dataset, coldata, sf)
    out <- run_core_analysis(
      count_mats = analysis_dataset_list[keys],
      coldata = coldata,
      dataset_names = paste(comparison_name, unname(track_short[track_key]), keys, sep = "_"),
      annot_df = annot_df,
      analysis_label = paste(comparison_name, analysis_label),
      size_factors = sf,
      gene_dds = cache$gene_dds,
      cooks_outliers = cooks,
      hc_p_fixed = if (a == "original") NULL else original_hc_p
    )
    if (a == "original") {
      original_hc_p <- out$null$hc_p_threshold[1]
      if (is.null(original_hc_p) || !is.finite(original_hc_p)) original_hc_p <- NA_real_
    } else {
      add_validation("One HBFSS threshold per comparison (Original and EVS share the Original's HC threshold)",
                     comparison_name, identical(unique(out$null$hc_p_threshold), original_hc_p),
                     sprintf("HC p = %.4g; Htau = %.3f", original_hc_p, out$null$hbfss_threshold[1]))
      remember_original_cooks(comparison_name, out$views$raw_dataset$results)
      # Reference baseline: the Original with DESeq2's default trend (once).
      tryCatch(record_default_trend_baseline(comparison_name, evs$raw_dataset, coldata, annot_df, sf, out),
               error = function(e) warning(comparison_name, ": DESeq2 default-trend baseline failed: ",
                                           conditionMessage(e), call. = FALSE))
    }
    for (nm in keys) {
      analysis_results[[nm]] <- list(dds = out$views[[nm]]$dds, results = out$views[[nm]]$results,
                                     dataset_mat = analysis_dataset_list[[nm]])
    }
    null_tab <- out$null
    null_tab$View <- vapply(null_tab$View, function(k) analysis_view_label(track_key, k), character(1))
    empirical_null_rows[[length(empirical_null_rows) + 1L]] <<- cbind(
      data.frame(Comparison = comparison_name, Analysis = analysis_label), null_tab)
  }

  reg <- paper_registry[[registry_key(comparison_name, track_key)]]
  reg$analysis_results <- analysis_results
  reg$size_factors <- sf
  paper_registry[[registry_key(comparison_name, track_key)]] <<- reg
  validate_track_normalization(comparison_name, track_key, analysis_results, sf, evs, coldata)
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
    track_key = "normalized_evs",
    dataset_key = c("raw_dataset", "leading_edge_dataset", "remainder_dataset"),
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
    # comparison-level panels generated after all three analysis views are fit.
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
      "The per-comparison cutoffs were: ",
      paste(
        paste0(comparison_table$comparison_name, " k=", vapply(comparison_table$comparison_name, get_active_evs_cutoff, integer(1))),
        collapse = ", "
      ),
      ". These cutoffs were chosen as described under Choice of the EVS cutoff k*, before any differential-expression test. ",
      "PCA preprocessing (EVS_PCA_METHOD = ", EVS_PCA_METHOD, "): ",
      switch(EVS_PCA_METHOD,
             VST = "DESeq2 varianceStabilizingTransformation (blind = TRUE) of the experiment-wide matrix (all 40 samples; every PAS with a nonzero count), with experiment-wide median-of-ratios size factors. ",
             CPM = "log1p of counts per million of the experiment-wide matrix (all 40 samples), with library sizes equal to the column totals. ",
             NORMALIZED = "experiment-wide DESeq2 median-of-ratios normalized counts (all 40 samples), without transformation. "),
      "This matrix was computed once, before any split, used only to rank PASs (never passed to DESeq2) and not transformed again. Rankings, cutoff rule, size factors, raw-count DESeq2 models and significance rules were identical across the preprocessing methods (VST, the primary analysis; CPM), each fitting its own regimes and choosing its own k* (EVS_PCA_Method_Comparison.csv, EVS_PCA_Ranking_Overlap.csv). ",
      "Within each condition, prcomp was applied to the transposed feature-by-sample matrix with centering and without scaling. ",
      "PASs were ranked by absolute PC1 loading, and the comparison-specific k highest-loading PASs were selected independently from the RT and ZT condition eigenvectors."
    ),
    "PASs present in both condition-specific top-k sets were classified as Joint; PASs present in only one set were classified as Disjoint for that condition. Joint plus both Disjoint sets formed the Leading Edge, and all other PASs formed the Remainder. EVS defined PAS membership only. Downstream DESeq2 analyses received the corresponding raw-count subset and estimated dispersion parameters within that analysis view. The three reported views were Original (No EVS), NormEVS Lead and NormEVS Rem.",
    "",
    "## Normalization",
    "Size factors were estimated once per comparison by DESeq2's median-of-ratios method on the full comparison matrix and supplied unchanged to the Original, Leading Edge and Remainder models so each PAS had identical normalized counts in the Original and in its EVS partition; this was verified for every view. Because a PAS's gene-wise dispersion depends only on its own counts, the design and the size factors, it is then also identical in every partition; it was therefore estimated once per comparison and reused by every view (only the trend, prior and MAP dispersions depend on the split), and its equality across views was verified. Estimating size factors within a partition would let the median-of-ratios reference change with the partition, confounding the effect of the split with a change of normalization.",
    "",
    "## Choice of the EVS cutoff k*",
    sprintf("One cutoff k* was selected for each comparison without any Wald test, p-value or discovery count, and without using the samples being tested. For a comparison (for example RT0 vs ZT6), the other six arms (30 samples) were analysed with DESeq2 (design ~ arm, experiment-wide median-of-ratios size factors) to obtain gene-wise dispersion estimates and mean counts. The log10 gene-wise dispersions were modelled as a two-component mixture, each component with mean trend log10(a0 + a1 x mean^-b) and its own Gaussian residual variance, fitted by expectation-maximization (weighted least-squares M-step, Nelder-Mead; convergence at a relative log-likelihood change below %g, at most %d iterations). As in DESeq2's own trend fit, estimates below %g were excluded from the fit; those PASs and PASs without counts in the held-out samples received the mixing proportion as their posterior. For each PAS, the likelihood ratio LR_i of its held-out dispersion under the Leading-Edge trend relative to the Remainder trend measured which trend it fitted better (equivalently q_i = LR_i / (1 + LR_i), the posterior with equal weights for the two trends; the mixing proportions were not used, so a PAS was judged by its fit alone and not penalized because the Leading-Edge trend is the smaller population). Excluded PASs had LR_i = 1 (no information). These quantities were estimated independently of the comparison's samples.", TWOTREND_EM_TOL, TWOTREND_EM_MAX_ITER, TWOTREND_MIN_DISP),
    sprintf("The comparison's PASs were ordered by their PC1 union depth (the smaller of their positions from the top of the |PC1| ranking of the two arms; the Leading Edge at cutoff k is exactly the PASs with depth <= k). For every k, the allocation score S(k) = sum over the Leading Edge of log q_i + sum over the Remainder of log(1 - q_i), equal up to a constant to the summed log likelihood ratio of the Leading Edge, measured how well the split agreed with the trend each PAS fits better; it was evaluated exactly for all k by cumulative sums. k* is therefore the cutoff at which further top-|PC1| PASs no longer fit the Leading-Edge trend better on balance. k* maximized S(k) among cutoffs leaving at least %d PASs in each set (k up to half the comparison's PASs); the mixture component used as the Leading-Edge trend was the labelling with the larger maximum, which was the component enriched among top-|PC1| PASs. The share of Leading-Edge PASs that fit the Leading-Edge trend better was reported, with the expected share of PASs on their own trend under the fitted mixture, (sum over LE of p_i + sum over REM of (1 - p_i)) / n with p_i the posterior including the mixing proportions, at k*, at the regime reference cutoff, at k = 5,000 and with all PASs in the Remainder. The held-out fit does not depend on the PCA preprocessing and was computed once per comparison.", TWOTREND_MIN_SET),
    "As a reference cutoff, shared rank regimes were also fitted. Within each arm, PCA of the experiment-wide preprocessing matrix ranked the N PASs by ascending |PC1 loading|; with PC1 contribution P_i = lambda_1 v_i1^2 and excess variance E_i = max(V_pool,i - mu_i,g, 0) (pooled within-arm variance of experiment-wide DESeq2-normalized counts minus the arm mean), the divergence curve D_g(r) = F_E,g(r) - F_P,g(r) of each arm was fitted jointly by a continuous two-knot piecewise-linear model, whose knots c1 < c2 defined Remainder, Divergence and Leading-Edge rank regimes. The regime k* was the Pareto-optimal k (Joint plus permissible Disjoint top-k PASs against opposite-arm Remainder crossings) farthest from the chord joining the scaled frontier endpoints (Figures S1 and S3).",
    "No significance computation occurred during cutoff selection (verified by a counter of Wald fits and significance calls); the held-out samples never overlapped the tested samples, and the use of k* and of the scored Leading Edge in the analysis were verified (Table_Validation_Checks.csv). Because the cutoff rule does not guarantee calibration, calibration was assessed separately (Figure S2).",
    "",
    "## DESeq2 model and log2 fold-change shrinkage",
    "Each analysis view was modeled independently with DESeq2 using a negative-binomial generalized linear model with design ~ condition and ZT/untrt as the reference. DESeq2 used the size factors described under Normalization and estimated gene-wise dispersions, the mean-dispersion trend, final dispersions, and Wald statistics for the RT/trt coefficient. The mean-dispersion trend in every view, including the Original, was alpha_tr(mean) = a0 + a1 mean^(-b), DESeq2's parametric form (b = 1; Love et al. 2014) with the exponent estimated, and was fitted on the log scale. DESeq2's default fit failed for the Remainder (a0 < 0), because nearly all well-expressed PASs fall in the Leading Edge and the Remainder's dispersions decline faster than 1/mean without a plateau; DESeq2 would then have substituted a local-regression trend for that group only. As a reference only, the Original was also refitted with DESeq2's default parametric trend (gamma GLM with b = 1 and DESeq2's automatic local fall-back), with every other step and testing rule unchanged (Table_Original_DESeq2_Default_Trend.csv), to show how the trend alone changes the Original's calls; the comparison of EVS with the Original used the same power trend in every view. DESeq2's gamma-GLM fit also places the trend at the arithmetic mean of the dispersions, above about two thirds of PASs, whereas its prior is normal on the log scale (log alpha_i ~ N(log alpha_tr, sigma_d^2)), whose centre is the mean of the log dispersions. Dispersion as a power of the mean is the NB-P model (Di et al. 2011). For each candidate b, a0 and a1 (both constrained to be at least 1e-8) minimized the sum of squared differences between log gene-wise dispersions and log trend; as in DESeq2's fit, PASs outside 1e-4 to 15 times the current curve were set aside and the fit repeated until the coefficients converged. b was chosen, every 0.1 from 0.5 to 3.5 and then every 0.01 around the best value, to minimize the median absolute residual of log gene-wise dispersion about log trend. Because a log gene-wise dispersion estimate is biased low by psi(nu/2) - log(nu/2) under the chi-square sampling model DESeq2 uses for its prior variance (nu = m - p residual degrees of freedom; -0.130 for 10 samples and 2 coefficients), a0 and a1 were then multiplied by exp(-(psi(nu/2) - log(nu/2))) (1.139), so that the prior was centred on the true rather than the estimated log dispersion; b and the prior variance were unchanged by this correction. The trend was supplied to DESeq2 through dispersionFunction(), and DESeq2's prior-variance estimate, MAP shrinkage, dispersion-outlier rule and Wald test were applied unchanged. Size factors, gene-wise dispersions, MAP dispersions and Wald tests were computed with estimateSizeFactors, estimateDispersionsGeneEst, estimateDispersionsMAP and nbinomWaldTest, the steps DESeq() runs. Log2 fold changes were then shrunken with apeglm. The apeglm-shrunken log2 fold change was the reported effect estimate, the effect term in HBFSS, the direction indicator, and the x-coordinate for all significance figures.",
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
    "## Multiple testing",
    "For each of the three DESeq2 tests (Wald, greaterAbs, lessAbs), the Benjamini-Hochberg correction was applied once per analysis: over the raw p-values of the Original, and over the raw p-values of the Leading Edge and the Remainder combined, which are disjoint and together contain every PAS of the comparison. DESeq2's independent filtering was not used, because it chooses a different mean-count filter in each view from that view's own results. Cook's-distance outliers were identified once per comparison, by DESeq2's default rule (Cook's distance above the 0.99 quantile of the F(p, m - p) distribution) in the Original's fit, and the p-values of those PASs were removed in every view, so the Original and the EVS analyses tested the same PASs. HBFSS used the Wald statistics, a separate empirical null for each view and one higher-criticism threshold per analysis (below), so neither the BH family nor the Cook's rule applied to it.",
    "",
    "## Empirical-null calibration, higher criticism, and HBFSS",
    sprintf(
      paste0(
        "The empirical null was estimated separately for each analysis view: the Original (No EVS), the Leading Edge and the Remainder, each from its own Wald statistics. The Leading Edge and Remainder are different populations, differing in expression range, dispersion law and fraction of changed PASs, and a null shared between them would be too strict for one and too loose for the other (Efron 2008). ",
        "The finite ordinary DESeq2 Wald statistics of each view were ", if (identical(EMPIRICAL_NULL_METHOD, "theoretical")) "referred to the theoretical N(0,1) null (reference analysis; no empirical correction). " else paste0("supplied to fdrtool with statistic='normal' and cutoff.method='", EMPIRICAL_NULL_METHOD, "'", if (identical(EMPIRICAL_NULL_METHOD, "pct0")) sprintf(" (central %.0f%%%% of the statistics)", 100 * EMPIRICAL_NULL_PCT0) else "", " to estimate a zero-centered normal empirical null and the corresponding empirical-null p-values. "), "The fitted null standard deviation and null fraction (eta0) were reported for every view. ",
        "These empirical p-values were distinct from the ordinary DESeq2 Wald p-values and from the greaterAbs/lessAbs p-values. ",
        "Higher-Criticism scores were then calculated once per comparison with fdrtool::hc.score, from the ordered empirical-null p-values of the Original, giving one HCp and one Htau per comparison that was applied to the Original and to the Leading Edge and Remainder alike. The HC maximum is sensitive to small changes in the p-value distribution, so separate thresholds for the Original and the split would make their HBFSS counts differ mainly through where each threshold fell; with one threshold, HBFSS differences between the Original and EVS reflect only each PAS's shrunken fold change and its empirical p-value, which was calibrated by its own view's empirical null. ",
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
    "PCA comparison figures used the same matrix that defined EVS in the run (the full-comparison matrix of the run's PCA preprocessing method). No additional transformation or post-split re-normalization was applied. Original, Leading Edge, and Remainder were compared within the same preprocessing scale. For each view, the full PCA eigenspectrum was used to calculate the percentage of total variance explained by PC1 and PC2. PC1 eigenvalue retention was calculated as the first eigenvalue of the view divided by the first eigenvalue of the corresponding Original matrix. Original-PC1 loading energy captured by a feature subset was calculated as the sum of squared Original PC1 loadings for that subset divided by the sum of squared Original PC1 loadings for all PASs. Leading Edge and Remainder loading-energy fractions were required to sum to 100%. Treatment-control separation was summarized in the PC1-PC2 plane as centroid distance divided by pooled within-group root-mean-square distance. PCA panels reported PAS count, PC1 explained variance, PC1 eigenvalue retention, Original-PC1 loading-energy fraction, and the separation ratio.",
    "",
    "## Volcano figures and tables",
    sprintf("All significance volcanoes used the HBFSS plotting coordinates: apeglm-shrunken log2 fold change on the x-axis and -log10(empirical-null p) on the y-axis. Standard DESeq2, Strong greaterAbs, and Weak lessAbs significance were determined from their own DESeq2 p-values and BH-adjusted p-values, then their markers were projected onto these common HBFSS coordinates only for visual comparison. The HBFSS boundary y=Htau/|LFC| and HCp reference were drawn on the same axes. The same key was used throughout: Std = green diamond, Strong = red square, Weak = blue triangle, HBFSS = purple star. Multiple method markers were superimposed at the same PAS coordinate. LessAbs-only PASs remained background. Each volcano reported the number of significant PASs for Std, Strong, Weak, HBFSS, and their HBFSS/DESeq2 overlap and labeled at most %d final significant PASs.", n_top_labels_volcano),
    "Each comparison/view folder contained one significant-PAS table and one method-count table. Significant-PAS tables contained only the union of final Standard, Strong, Weak, and HBFSS discoveries and reported the DESeq2 p-values/BH-adjusted p-values used by the applicable DESeq2 tests together with empirical p, HBFSS, HC alpha0, HCp, Htau, apeglm LFC, and explicit method indicators. Method-count tables also reported HC alpha0, HCp, and Htau for view-level calibration auditing.",
    "",
    "## 3'aTWAS ortholog overlap",
    "Human 3'aTWAS gene symbols were mapped to rat gene symbols by combining database-supported babelgene human-to-rat ortholog mappings (top=FALSE) with direct case-insensitive symbol-equivalent matches present in the WTTS annotation. For every comparison and analysis view, mapped TWAS orthologs were intersected with Standard, Strong, Weak, and HBFSS WTTS discoveries. Concise TWAS tables reported the human TWAS symbol, rat ortholog, total TWAS records, distinct TWAS transcript count, TWAS multi-transcript/APA status, total WTTS PAS count, WTTS multi-PAS/APA status, significant WTTS PAS identifiers, significant multi-PAS status, and the method(s) and analysis view in which significance was observed.",
    "",
    "## Output organization",
    "Original (No EVS), NormEVS Lead and NormEVS Rem were written to separate folders within each comparison for view-specific tables. The two-trend cutoff was written to Table_TwoTrend_Heldout_Fits.csv (trend parameters, shares, EM convergence), Table_TwoTrend_PAS_Posteriors.csv (held-out mean and dispersion, P(Leading-Edge trend) and union depth of every PAS), Table_TwoTrend_Score_Scans.csv (S(k) and the expected share on-trend at every k), Table_TwoTrend_Leading_Edge_Sites.csv and Table_TwoTrend_Cutoff_Summary.csv, and plotted in Figure 1. The regime reference was written to Table_Regime_Knots.csv, Table_Regime_D_Curves.csv, Table_Regime_Cutoff_Scans.csv, Table_Regime_Pareto_Frontiers.csv, Table_Regime_TopK_Sites.csv and Table_Regime_Cutoff_Summary.csv, and plotted in Figures S1 and S3. The automated checks were written to Table_Validation_Checks.csv. Manuscript figures were organized as comparison-level panels so related volcano, PCA/EVS, and 3'aTWAS results could be reviewed side by side without redundant individual images. The curated figure archive SEQUENCE_ALL_FIGURES.zip contained manuscript PNG panels only. The separate archive SEQUENCE_ALL_TABLES.zip contained the locked EVS cutoff table, the empirical-null table for every analysis, the dispersion-after-EVS tables, the EVS split audit table, overall differential-expression method-count table, quantitative EVS/PCA evidence table, one all-view significant-PAS table and method-count table for each comparison, and the combined 3'aTWAS PAS-, gene-, and method-summary tables."
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
# EVS aims to reduce.
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
  save_regime_outputs()
  save_twotrend_outputs()

  if (!is.null(pas) && !is.null(trd) && !is.null(met)) try_fig({
    fig2 <- assemble_rows(list(panel_dispersion_grid(pas, trd, met, "a"),
                               list(plots = list(panel_trend_exponent(met, "b"), panel_prior_width(met, "c")),
                                    widths = c(1, 1))),
                          heights = c(2.2, 0.9),
                          title = paste0("Figure 2. DESeq2 dispersions, trend shape and prior width, Original vs split (",
                                         CUTOFF_METHOD_LABEL, ")"))
    save_figure(fig2, paper_fig_dir, paste0("Figure_2_Dispersion_", pca_run_folder(EVS_PCA_METHOD)), 180, 215)
  }, "Figure 2")
  invisible(TRUE)
}

# =============================================================================
# REGIME CUTOFF k*: SHARED RANK REGIMES AND COUNT-BASED PARETO SELECTION
# =============================================================================
# No DESeq2 fit, p-value or discovery count enters the choice of k.
#
# 1. Experiment-wide PCA input. All 40 samples and every PAS with a nonzero
#    count anywhere (the PAS universe, N PASs) are transformed once by the
#    run's EVS_PCA_METHOD (VST, CPM or NORMALIZED; see EIGENVECTOR SPLITTING).
#    Within each of the eight arms (RT0, RT2, RT4, RT8, ZT6, ZT8, ZT10, ZT14;
#    five samples each) prcomp (centred, unscaled) ranks PASs by |PC1 loading|
#    (rank 1 = smallest |loading|, rank N = largest).
# 2. Divergence curve per arm. PC1 mass p_i = lambda_1 v_i1^2 / sum; excess-
#    variance mass q_i = E_i / sum, with E_i = max(V_pool,i - mu_i,arm, 0):
#    the pooled within-arm variance of the experiment-wide DESeq2-normalized
#    counts minus the arm mean (variance beyond Poisson). Along the rank order,
#    D(r) = F_E(r) - F_P(r), the cumulative excess-variance mass minus the
#    cumulative PC1 mass.
# 3. Shared regimes. The eight D(r) curves are fitted jointly by one
#    continuous two-knot piecewise-linear model (least squares, Nelder-Mead
#    over the knots from four starts). Ranks below c1 form the Remainder
#    regime, ranks c1..c2 the Divergence regime, ranks above c2 the
#    Leading-Edge regime.
# 4. Per comparison, for every k = 1..N - c2: G(k) = Joint top-k PASs plus
#    Disjoint ones whose rank in the opposite arm is in the Leading-Edge or
#    Divergence regime; R(k) = Disjoint PASs whose opposite-arm rank falls in
#    the Remainder regime. On the Pareto frontier of (R, G), both are min-max
#    scaled; k* is the frontier point farthest from the chord joining the two
#    frontier endpoints (ties: larger G, smaller R, larger k).
# The split then takes the top k* |PC1| PASs of each arm of the comparison,
# exactly the sets scored here.

REGIME_ARMS <- c("R0", "ZT6", "R2", "ZT8", "R4", "ZT10", "R8", "ZT14")
REGIME_ARM_LABELS <- c(R0 = "RT0", ZT6 = "ZT6", R2 = "RT2", ZT8 = "ZT8", R4 = "RT4", ZT10 = "ZT10",
                       R8 = "RT8", ZT14 = "ZT14")
EXPERIMENT_CACHE <- new.env(parent = emptyenv())
REGIME <- NULL      # set by compute_regime_cutoffs()

# Experiment-wide raw counts: every PAS with a nonzero count in any sample.
experiment_counts <- function() {
  if (!is.null(EXPERIMENT_CACHE$raw)) return(EXPERIMENT_CACHE$raw)
  x <- as.matrix(WTTS_Seq[, meta_all$id, drop = FALSE])
  storage.mode(x) <- "numeric"
  rownames(x) <- rownames(WTTS_Seq)
  x[!is.finite(x)] <- 0
  x <- x[rowSums(x) > 0, , drop = FALSE]
  EXPERIMENT_CACHE$raw <- coerce_raw_count_matrix_for_deseq2(x, context = "experiment-wide count matrix")
  EXPERIMENT_CACHE$raw
}

experiment_arm <- function(samples) {
  g <- rep(NA_character_, length(samples))
  for (a in REGIME_ARMS) {
    hit <- grepl(paste0("^", a, "_"), samples)
    if (any(!is.na(g[hit]))) stop("A sample matched more than one arm: ", a)
    g[hit] <- a
  }
  if (anyNA(g)) stop("Samples without an arm: ", paste(samples[is.na(g)], collapse = ", "))
  factor(g, levels = REGIME_ARMS)
}

# Experiment-wide DESeq2 object with median-of-ratios size factors (used for
# the PCA input and the excess variance only; every DE test uses the
# comparison's own shared size factors).
experiment_dds <- function() {
  if (!is.null(EXPERIMENT_CACHE$dds)) return(EXPERIMENT_CACHE$dds)
  raw <- experiment_counts()
  dds <- DESeq2::DESeqDataSetFromMatrix(raw, colData = data.frame(arm = experiment_arm(colnames(raw)),
                                                                   row.names = colnames(raw)), design = ~ arm)
  dds <- tryCatch(DESeq2::estimateSizeFactors(dds),
                  error = function(e) DESeq2::estimateSizeFactors(dds, type = "poscounts"))
  EXPERIMENT_CACHE$dds <- dds
  dds
}

# The experiment-wide PCA input of a method (cached).
experiment_rank_matrix <- function(method = EVS_PCA_METHOD) {
  method <- match.arg(method, PCA_METHOD_LEVELS)
  key <- paste0("rank_", method)
  if (!is.null(EXPERIMENT_CACHE[[key]])) return(EXPERIMENT_CACHE[[key]])
  raw <- experiment_counts()
  x <- switch(method,
    VST = SummarizedExperiment::assay(suppressMessages(DESeq2::varianceStabilizingTransformation(experiment_dds(),
                                                                                                blind = TRUE))),
    NORMALIZED = DESeq2::counts(experiment_dds(), normalized = TRUE),
    CPM = {
      lib <- colSums(raw)
      lib[!is.finite(lib) | lib <= 0] <- 1
      log1p(sweep(raw, 2, lib / 1e6, "/"))
    })
  x <- as.matrix(x)
  dimnames(x) <- dimnames(raw)
  if (any(!is.finite(x))) stop("The experiment-wide ", method, " matrix has non-finite values.")
  EXPERIMENT_CACHE[[key]] <- x
  x
}

# One arm's ranking and divergence curve (rank 1 = smallest |PC1 loading|).
build_regime_arm <- function(rank_matrix, norm, pooled_var, arm_samples) {
  z <- rank_matrix[, arm_samples, drop = FALSE]
  p <- stats::prcomp(t(z), center = TRUE, scale. = FALSE, rank. = 1)
  loading <- p$rotation[, 1]
  loading[!is.finite(loading)] <- 0
  ord <- order(abs(loading), decreasing = FALSE)
  P <- (p$sdev[1]^2 * loading^2)[ord]
  E <- pmax(pooled_var - rowMeans(norm[, arm_samples, drop = FALSE]), 0)[ord]
  if (sum(P) <= 0 || sum(E) <= 0) stop("Undefined PC1 or excess-variance mass for an arm.")
  data.frame(rank = seq_along(ord), feature_id = rownames(rank_matrix)[ord],
             abs_pc1_loading = abs(loading[ord]), pc1_mass = P / sum(P), excess_mass = E / sum(E),
             F_P = cumsum(P / sum(P)), F_E = cumsum(E / sum(E)), D = cumsum(E / sum(E)) - cumsum(P / sum(P)),
             stringsAsFactors = FALSE)
}

regime_basis <- function(x, c1, c2) cbind(1, x, pmax(x - c1, 0), pmax(x - c2, 0))

# Shared two-knot piecewise-linear fit of all arms' D(r).
fit_shared_knots <- function(arms) {
  N <- nrow(arms[[1]])
  if (length(unique(vapply(arms, nrow, integer(1)))) != 1L) stop("All arms must use the same PAS universe.")
  x_full <- (seq_len(N) - 1) / (N - 1)
  D_full <- do.call(cbind, lapply(arms, function(z) z$D))
  idx <- unique(as.integer(round(seq(1, N, length.out = min(4000L, N)))))
  x <- x_full[idx]; D <- D_full[idx, , drop = FALSE]
  min_gap <- max(4 / (N - 1), 1e-4)
  objective <- function(par) {
    c1 <- par[1]; c2 <- par[2]
    if (!is.finite(c1) || !is.finite(c2) || c1 <= 0 || c2 >= 1 || c2 - c1 <= min_gap) return(1e100)
    X <- regime_basis(x, c1, c2)
    b <- tryCatch(qr.coef(qr(X), D), error = function(e) NULL)
    if (is.null(b) || any(!is.finite(b))) return(1e100)
    sum((D - X %*% b)^2)
  }
  starts <- list(c(.05, .95), c(.12, .88), c(.22, .78), c(.32, .68))
  fits <- lapply(starts, function(st) stats::optim(st, objective, method = "Nelder-Mead",
                                                   control = list(maxit = 500, reltol = 1e-10)))
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  c1 <- as.integer(round(1 + best$par[1] * (N - 1)))
  c2 <- as.integer(round(1 + best$par[2] * (N - 1)))
  c1 <- max(2L, min(N - 2L, c1)); c2 <- max(c1 + 1L, min(N - 1L, c2))
  X_full <- regime_basis(x_full, (c1 - 1) / (N - 1), (c2 - 1) / (N - 1))
  fitted <- X_full %*% qr.coef(qr(X_full), D_full)
  list(c1 = c1, c2 = c2, N = N, arms = names(arms), fitted = fitted, SSE = sum((D_full - fitted)^2))
}

regime_rank_map <- function(df) stats::setNames(df$rank, df$feature_id)

regime_add_events <- function(pos, weight, L) {
  out <- numeric(L)
  ok <- is.finite(pos) & is.finite(weight) & pos >= 1 & pos <= L
  if (!any(ok)) return(out)
  z <- rowsum(weight[ok], group = as.integer(pos[ok]), reorder = FALSE)
  out[as.integer(rownames(z))] <- z[, 1]
  out
}

# Weight active for start <= k < end, k = 1..K.
regime_active_interval <- function(start, end, weight, K) {
  ok <- is.finite(start) & is.finite(end) & is.finite(weight) & start >= 1 & start <= K & end > start
  if (!any(ok)) return(numeric(K))
  delta <- regime_add_events(as.integer(start[ok]), weight[ok], K + 1L) -
    regime_add_events(pmin(as.integer(end[ok]), K + 1L), weight[ok], K + 1L)
  cumsum(delta)[seq_len(K)]
}

# G(k), R(k) and the top-k union composition for every k = 1..N - c2.
scan_regime_pair <- function(arm_a, arm_b, c1, c2) {
  N <- nrow(arm_a); K <- N - c2
  if (K < 1L) stop("No candidate k exists beyond c2.")
  ids <- arm_a$feature_id
  rA <- as.integer(regime_rank_map(arm_a)[ids]); rB <- as.integer(regime_rank_map(arm_b)[ids])
  dA <- N - rA + 1L; dB <- N - rB + 1L          # position from the top (1 = largest |loading|)
  joint_n <- cumsum(regime_add_events(pmax(dA, dB), rep(1, N), K))
  one_side <- function(d_self, d_other, r_other) {
    i <- which(d_self < d_other & d_self <= K)
    s <- d_self[i]; e <- pmin(d_other[i], K + 1L); o <- r_other[i]
    list(n = regime_active_interval(s, e, rep(1, length(i)), K),
         rem = regime_active_interval(s, e, as.numeric(o < c1), K),
         div = regime_active_interval(s, e, as.numeric(o >= c1 & o <= c2), K),
         le = regime_active_interval(s, e, as.numeric(o > c2), K))
  }
  a <- one_side(dA, dB, rB); b <- one_side(dB, dA, rA)
  opp_le <- a$le + b$le; opp_div <- a$div + b$div; rem <- a$rem + b$rem
  good <- joint_n + opp_le + opp_div
  union_n <- 2 * seq_len(K) - joint_n
  if (!isTRUE(all.equal(good + rem, union_n))) stop("Internal union-count mismatch in the regime scan.")
  data.frame(k = seq_len(K), joint_n = joint_n, disjoint_n = a$n + b$n, opposite_le_n = opp_le,
             opposite_divergence_n = opp_div, remainder_cross_n = rem, good_n = good, union_n = union_n,
             retained_fraction = ifelse(union_n > 0, good / union_n, NA_real_),
             remainder_cross_fraction = ifelse(union_n > 0, rem / union_n, NA_real_))
}

# Pareto frontier of (R, G) and the maximum endpoint-chord deviation.
select_regime_k <- function(scan) {
  tmp <- data.frame(k = scan$k, good = scan$good_n, cost = scan$remainder_cross_n)
  tmp <- tmp[order(tmp$cost, -tmp$good, -tmp$k), , drop = FALSE]
  tmp <- tmp[!duplicated(tmp$cost), , drop = FALSE]
  best_before <- c(-Inf, utils::head(cummax(tmp$good), -1L))
  fr <- tmp[tmp$good > best_before, , drop = FALSE]
  fr <- fr[order(fr$cost, fr$good, fr$k), , drop = FALSE]
  if (nrow(fr) < 2L) stop("At least two Pareto-optimal points are needed for the endpoint chord.")
  sc <- function(x, r) if (diff(r) == 0) rep(0, length(x)) else (x - r[1]) / diff(r)
  gr <- range(fr$good); cr <- range(fr$cost)
  fr$good_norm <- sc(fr$good, gr); fr$remainder_norm <- sc(fr$cost, cr)
  x1 <- fr$remainder_norm[1]; y1 <- fr$good_norm[1]
  x2 <- fr$remainder_norm[nrow(fr)]; y2 <- fr$good_norm[nrow(fr)]
  den <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
  if (!is.finite(den) || den <= .Machine$double.eps) stop("Coincident Pareto endpoints: the chord is undefined.")
  dev <- function(xn, yn) abs((y2 - y1) * xn - (x2 - x1) * yn + x2 * y1 - y2 * x1) / den
  fr$endpoint_deviation <- dev(fr$remainder_norm, fr$good_norm)
  top <- fr[abs(fr$endpoint_deviation - max(fr$endpoint_deviation)) < 1e-12, , drop = FALSE]
  top <- top[order(-top$good, top$cost, -top$k), , drop = FALSE]
  scan$is_pareto <- paste(scan$remainder_cross_n, scan$good_n) %in% paste(fr$cost, fr$good)
  scan$good_norm <- sc(scan$good_n, gr); scan$remainder_norm <- sc(scan$remainder_cross_n, cr)
  scan$endpoint_deviation <- dev(scan$remainder_norm, scan$good_norm)
  scan$Chosen <- scan$k == top$k[1]
  list(scan = scan, frontier = fr, k = as.integer(top$k[1]))
}

# Joint / Disjoint classes and opposite-arm regimes of the top-k union.
classify_regime_sites <- function(k, comparison, arm_a, arm_b, c1, c2) {
  rA <- regime_rank_map(arm_a); rB <- regime_rank_map(arm_b)
  topA <- utils::tail(arm_a$feature_id, k); topB <- utils::tail(arm_b$feature_id, k)
  ids <- union(topA, topB)
  a <- as.integer(rA[ids]); b <- as.integer(rB[ids])
  inA <- ids %in% topA; inB <- ids %in% topB
  joint <- inA & inB
  opp <- ifelse(inA & !inB, b, ifelse(inB & !inA, a, NA_integer_))
  region <- ifelse(joint, "Joint", ifelse(opp > c2, "Opposite Leading Edge",
                                          ifelse(opp >= c1, "Opposite Divergence", "Opposite Remainder")))
  data.frame(Comparison = comparison, feature_id = ids, k = k, RT_rank = a, ZT_rank = b,
             Class = ifelse(joint, "Joint", ifelse(inA, "RT only", "ZT only")), Opposite_Region = region,
             High_Confidence = region != "Opposite Remainder", stringsAsFactors = FALSE)
}

# The cutoffs of every comparison for one PCA method (run once per run).
compute_regime_cutoffs <- function(method = EVS_PCA_METHOD) {
  raw <- experiment_counts()
  arm <- experiment_arm(colnames(raw))
  norm <- DESeq2::counts(experiment_dds(), normalized = TRUE)
  sse <- numeric(nrow(norm)); df <- 0L
  for (a in levels(arm)) {
    z <- norm[, arm == a, drop = FALSE]
    sse <- sse + rowSums((z - rowMeans(z))^2); df <- df + ncol(z) - 1L
  }
  pooled_var <- pmax(sse / df, 0)
  rm <- experiment_rank_matrix(method)
  arms <- lapply(stats::setNames(REGIME_ARMS, REGIME_ARMS), function(a)
    build_regime_arm(rm, norm, pooled_var, colnames(raw)[arm == a]))
  knot <- fit_shared_knots(arms)
  message(sprintf("Regime cutoffs (%s): N = %d PASs, c1 = %d, c2 = %d, candidate k = 1..%d",
                  method, knot$N, knot$c1, knot$c2, knot$N - knot$c2))
  scans <- list(); frontiers <- list(); sites <- list(); summary <- list()
  for (i in seq_len(nrow(comparison_table))) {
    cmp <- comparison_table$comparison_name[i]
    A <- arms[[comparison_table$group1_prefix[i]]]; B <- arms[[comparison_table$group2_prefix[i]]]
    sel <- select_regime_k(scan_regime_pair(A, B, knot$c1, knot$c2))
    k <- sel$k
    row <- sel$scan[sel$scan$k == k, , drop = FALSE]
    REGIME_K[[cmp]] <<- k
    scans[[cmp]] <- cbind(Comparison = cmp, EVS_PCA_Method = method, sel$scan)
    frontiers[[cmp]] <- cbind(Comparison = cmp, EVS_PCA_Method = method, sel$frontier)
    sites[[cmp]] <- classify_regime_sites(k, cmp, A, B, knot$c1, knot$c2)
    summary[[cmp]] <- data.frame(
      Comparison = cmp, EVS_PCA_Method = method, N = knot$N, c1 = knot$c1, c2 = knot$c2,
      Remainder_regime = knot$c1 - 1L, Divergence_regime = knot$c2 - knot$c1 + 1L,
      Leading_Edge_regime = knot$N - knot$c2, k = k, k_over_LE_regime = k / (knot$N - knot$c2),
      G = row$good_n, R = row$remainder_cross_n, Joint = row$joint_n, Opposite_LE = row$opposite_le_n,
      Opposite_Divergence = row$opposite_divergence_n, Union = row$union_n,
      Endpoint_Deviation = row$endpoint_deviation, Pareto_Points = nrow(sel$frontier), Knot_SSE = knot$SSE)
    message(sprintf("  %s: regime k* = %d (union %d; Joint %d, opposite LE %d, opposite Divergence %d, Remainder crossings %d)",
                    cmp, k, row$union_n, row$joint_n, row$opposite_le_n, row$opposite_divergence_n,
                    row$remainder_cross_n))
  }
  REGIME <<- list(method = method, arms = arms, knot = knot, scans = scans, frontiers = frontiers,
                  sites = sites, summary = dplyr::bind_rows(summary))
  invisible(REGIME)
}

# Position from the top (1 = largest |PC1 loading|) of each comparison PAS in
# its two arms, and the union depth: a PAS is in the Leading Edge at k exactly
# when its union depth is <= k.
regime_depths <- function(comparison_name, ids) {
  i <- match(comparison_name, comparison_table$comparison_name)
  N <- REGIME$knot$N
  pos <- function(arm) N - as.integer(regime_rank_map(REGIME$arms[[arm]])[ids]) + 1L
  rt <- pos(comparison_table$group1_prefix[i]); zt <- pos(comparison_table$group2_prefix[i])
  if (anyNA(rt) || anyNA(zt)) stop(comparison_name, ": PASs missing from the experiment-wide ranking")
  list(RT = stats::setNames(rt, ids), ZT = stats::setNames(zt, ids), depth = stats::setNames(pmin(rt, zt), ids))
}

# Comparison inputs computed once and reused by the analysis: shared size
# factors and gene-wise dispersions (also shared across PCA-method runs
# through Shared_Cache, with an identity check).
prepare_comparison_cache <- function(comparison_name, count_matrix, coldata) {
  raw <- coerce_raw_count_matrix_for_deseq2(count_matrix, context = paste0(comparison_name, " comparison matrix"))
  sf <- comparison_size_factors(raw, coldata)
  shared_file <- if (nzchar(SHARED_CACHE_DIR)) file.path(SHARED_CACHE_DIR, paste0(comparison_name, ".rds")) else ""
  shared_check <- NA
  shared_detail <- "first PCA method run for this comparison (shared inputs written)"
  gdds <- NULL
  if (nzchar(shared_file) && file.exists(shared_file)) {
    sh <- readRDS(shared_file)
    shared_check <- identical(rownames(sh$raw), rownames(raw)) && identical(colnames(sh$raw), colnames(raw)) &&
      identical(unname(sh$raw), unname(raw)) && identical(sh$coldata, coldata) && identical(sh$sf, sf)
    shared_detail <- sprintf("%s vs %s run: %d PASs, samples %s", EVS_PCA_METHOD, sh$method, nrow(raw),
                             paste(colnames(raw), collapse = ","))
    if (isTRUE(shared_check)) gdds <- sh$gene_dds
  }
  if (is.null(gdds)) gdds <- comparison_gene_dds(raw, coldata, sf)
  if (nzchar(shared_file) && !file.exists(shared_file)) {
    saveRDS(list(raw = raw, coldata = coldata, sf = sf, gene_dds = gdds, method = EVS_PCA_METHOD), shared_file)
  }
  dp <- regime_depths(comparison_name, rownames(raw))
  comparison_cache[[comparison_name]] <<- list(raw = raw, coldata = coldata, depth = dp$depth, rank_RT = dp$RT,
                                               rank_ZT = dp$ZT, sf = sf, gene_dds = gdds,
                                               shared_check = shared_check, shared_detail = shared_detail)
  invisible(TRUE)
}

# ---- regime outputs (tables and figures of this run) --------------------------

save_regime_outputs <- function() {
  if (is.null(REGIME)) return(invisible(FALSE))
  out <- summary_table_dir
  m <- REGIME$method; kn <- REGIME$knot; N <- kn$N
  save_csv(REGIME$summary, file.path(out, "Table_Regime_Cutoff_Summary.csv"))
  save_csv(dplyr::bind_rows(REGIME$scans), file.path(out, "Table_Regime_Cutoff_Scans.csv"))
  save_csv(dplyr::bind_rows(REGIME$frontiers), file.path(out, "Table_Regime_Pareto_Frontiers.csv"))
  save_csv(dplyr::bind_rows(REGIME$sites), file.path(out, "Table_Regime_TopK_Sites.csv"))
  thin <- unique(c(1L, as.integer(round(seq(1, N, length.out = 2000L))), N))
  curves <- dplyr::bind_rows(lapply(seq_along(REGIME$arms), function(j) data.frame(
    EVS_PCA_Method = m, Arm = REGIME_ARM_LABELS[[names(REGIME$arms)[j]]], Rank = thin,
    D = REGIME$arms[[j]]$D[thin], Fitted = kn$fitted[thin, j])))
  save_csv(curves, file.path(out, "Table_Regime_D_Curves.csv"))
  save_csv(data.frame(EVS_PCA_Method = m, N = N, c1 = kn$c1, c2 = kn$c2, SSE = kn$SSE),
           file.path(out, "Table_Regime_Knots.csv"))
  try_fig(save_figure(plot_regime_definition(curves, data.frame(EVS_PCA_Method = m, N = N, c1 = kn$c1, c2 = kn$c2)),
                      paper_fig_dir, paste0("Figure_S1_Regime_Definition_", m), 180, 95), "Figure S1 (regimes)")
  try_fig(save_figure(plot_regime_pareto(dplyr::bind_rows(REGIME$scans), dplyr::bind_rows(REGIME$frontiers), m),
                      paper_fig_dir, paste0("Figure_S3_Regime_Pareto_Cutoffs_", m), 180, 150), "Figure S3 (regime Pareto)")
  invisible(TRUE)
}

# Per comparison: all candidate k in (R, G), the Pareto frontier, the endpoint
# chord and k*.
plot_regime_pareto <- function(scans, frontiers, method) {
  scans$Comparison <- as_comp(scans$Comparison); frontiers$Comparison <- as_comp(frontiers$Comparison)
  thin <- do.call(rbind, lapply(split(scans, scans$Comparison, drop = TRUE), function(s)
    s[unique(round(seq(1, nrow(s), length.out = min(nrow(s), 2500L)))), , drop = FALSE]))
  chord <- do.call(rbind, lapply(split(frontiers, frontiers$Comparison, drop = TRUE), function(f)
    f[c(1L, nrow(f)), c("Comparison", "cost", "good")]))
  sel <- scans[scans$Chosen, , drop = FALSE]
  sel$lab <- sprintf("k* = %s", format(sel$k, big.mark = ","))
  ggplot() +
    geom_path(data = thin, aes(remainder_cross_n, good_n, colour = "All candidate k"), linewidth = 0.3) +
    geom_line(data = frontiers, aes(cost, good, colour = "Pareto frontier"), linewidth = 0.6) +
    geom_line(data = chord, aes(cost, good, colour = "Endpoint chord"), linewidth = 0.4, linetype = "22") +
    geom_point(data = sel, aes(remainder_cross_n, good_n), shape = 23, size = 2, fill = "#CC79A7", colour = "white") +
    geom_text(data = sel, aes(remainder_cross_n, good_n, label = lab), hjust = 1.1, vjust = -0.6, size = 2.3,
              fontface = "bold", colour = "#CC79A7") +
    facet_wrap(~ Comparison, nrow = 2, scales = "free", labeller = labeller(Comparison = comp_lab())) +
    scale_colour_manual(values = c("All candidate k" = "grey78", "Pareto frontier" = "#D55E00",
                                   "Endpoint chord" = "grey45"), name = NULL) +
    scale_x_continuous(labels = scales::label_comma()) + scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "Opposite-arm Remainder crossings, R(k)", y = "Joint + permissible Disjoint, G(k)",
         title = paste0("Regime reference cutoff per comparison (", method, "-EVS)"),
         subtitle = "Regime k* is the Pareto-optimal k farthest from the chord joining the scaled frontier endpoints.") +
    fig_theme() + theme(legend.position = "bottom")
}

# =============================================================================
# TWO-TREND CUTOFF k*: HELD-OUT DISPERSION TRENDS DETERMINE WHERE k SITS
# =============================================================================
# The cutoff is placed where the PC1 split best agrees with two dispersion
# trends that are estimated first, without any k and without the samples
# being tested. No Wald test, p-value or discovery count is used.
#
# 1. Held-out dispersions. For a comparison (e.g. RT0 vs ZT6) the other six
#    arms (30 samples) are used: DESeq2 gene-wise dispersion estimates with the
#    design ~ arm (each arm its own mean) and experiment-wide median-of-ratios
#    size factors; mean = DESeq2 baseMean of those samples.
# 2. Two trends. log10 gene-wise dispersion is modelled as a two-component
#    mixture: component c has mean trend log10(a0_c + a1_c * mean^-b_c) and
#    Gaussian residual variance s2_c; mixing proportion pi. EM (weighted least
#    squares M-step, Nelder-Mead) runs to convergence. Estimates below
#    100 x DESeq2's minimum dispersion are left out of the fit (as in DESeq2's
#    own trend fit) and, with all-zero PASs, get the prior pi as posterior.
#    Output: for each PAS, how much better it fits the Leading-Edge trend than
#    the Remainder trend: the likelihood ratio LR_i (equivalently q_i =
#    LR_i / (1 + LR_i), the posterior with equal weights for the two trends).
#    The mixing proportions are not used, so a PAS is judged only by its fit,
#    not penalized because the Leading-Edge trend is the smaller population.
# 3. Where k sits. Along the comparison's PC1 union depth (the regime
#    ranking of the run's PCA method, so the Leading Edge at k is
#    {depth <= k}), the allocation score is
#        S(k) = sum_{LE(k)} log q_i + sum_{REM(k)} log(1 - q_i)
#             = const + sum_{LE(k)} log LR_i,
#    evaluated exactly for every k by cumulative sums. k* = argmax S(k) over
#    k with at least TWOTREND_MIN_SET PASs in each set. Which mixture
#    component is the Leading-Edge trend is the labelling with the larger
#    maximum S (it is the component enriched among top-|PC1| PASs).
# The held-out fit does not depend on the PCA method and is computed once per
# comparison (Shared_Cache). The regime k* (REGIME CUTOFF k*) is reported
# alongside as a reference.

TWOTREND_EM_MAX_ITER <- 3000L
TWOTREND_EM_TOL <- 1e-10         # relative change in the log-likelihood
TWOTREND_MIN_DISP <- 1e-6        # 100 x DESeq2 minDisp: lower estimates are not fitted
TWOTREND_MIN_SET <- 200L         # minimum PASs in the Leading Edge and in the Remainder
TWOTREND_K_REFERENCE <- 5000L
TWOTREND <- NULL                 # set by compute_twotrend_cutoffs()
REGIME_K <- list()               # regime k* per comparison (reference)

twotrend_curve <- function(par, x) log10(exp(par[1]) + exp(par[2]) * x^(-par[3]))

twotrend_fit_weighted <- function(x, d, w, par0) {
  f <- function(p) sum(w * (d - twotrend_curve(p, x))^2)
  stats::optim(par0, f, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-12))$par
}

# Two-component mixture of log10 dispersion about two power trends (EM).
twotrend_em <- function(x, d) {
  p1 <- twotrend_fit_weighted(x, d, rep(1, length(d)), c(log(0.05), 0, 1))
  res <- d - twotrend_curve(p1, x)
  w <- as.numeric(res > stats::median(res))     # start: above vs below one trend
  pars <- list(twotrend_fit_weighted(x, d, w, p1), twotrend_fit_weighted(x, d, 1 - w, p1))
  pi1 <- mean(w); ll_old <- -Inf; converged <- FALSE; it <- 0L
  comp_ll <- function(par, s2, pr) log(pr) - 0.5 * log(2 * pi * s2) - (d - twotrend_curve(par, x))^2 / (2 * s2)
  for (it in seq_len(TWOTREND_EM_MAX_ITER)) {
    s2 <- c(sum(w * (d - twotrend_curve(pars[[1]], x))^2) / sum(w),
            sum((1 - w) * (d - twotrend_curve(pars[[2]], x))^2) / sum(1 - w))
    l1 <- comp_ll(pars[[1]], s2[1], pi1); l2 <- comp_ll(pars[[2]], s2[2], 1 - pi1)
    mx <- pmax(l1, l2)
    ll <- sum(mx + log(exp(l1 - mx) + exp(l2 - mx)))
    w <- stats::plogis(l1 - l2)
    pi1 <- mean(w)
    pars <- list(twotrend_fit_weighted(x, d, w, pars[[1]]), twotrend_fit_weighted(x, d, 1 - w, pars[[2]]))
    if (is.finite(ll_old) && abs(ll - ll_old) <= TWOTREND_EM_TOL * abs(ll)) { converged <- TRUE; break }
    ll_old <- ll
  }
  list(pars = pars, s2 = s2, pi = c(pi1, 1 - pi1), loglik = ll, iterations = it, converged = converged)
}

# Posterior of component 1 for any (x, d), from a fitted mixture.
twotrend_posterior <- function(mix, x, d) {
  l <- lapply(1:2, function(c) log(mix$pi[c]) - 0.5 * log(mix$s2[c]) -
                (d - twotrend_curve(mix$pars[[c]], x))^2 / (2 * mix$s2[c]))
  stats::plogis(l[[1]] - l[[2]])
}

# Held-out gene-wise dispersions and the two-trend mixture for one comparison
# (independent of the PCA method; cached across method runs).
twotrend_heldout_fit <- function(comparison_name, tested_arms) {
  f <- if (nzchar(SHARED_CACHE_DIR)) file.path(SHARED_CACHE_DIR, paste0("twotrend_", comparison_name, ".rds")) else ""
  if (nzchar(f) && file.exists(f)) {
    fit <- readRDS(f); fit$from_cache <- TRUE
    return(fit)
  }
  raw <- experiment_counts()
  arm <- experiment_arm(colnames(raw))
  ho <- colnames(raw)[!as.character(arm) %in% tested_arms]
  cd <- data.frame(arm = droplevels(arm[match(ho, colnames(raw))]), row.names = ho)
  dds <- DESeq2::DESeqDataSetFromMatrix(raw[, ho, drop = FALSE], colData = cd, design = ~ arm)
  DESeq2::sizeFactors(dds) <- DESeq2::sizeFactors(experiment_dds())[ho]
  dds <- suppressMessages(DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE))
  mc <- S4Vectors::mcols(dds)
  x <- mc$baseMean; disp <- mc$dispGeneEst
  zero <- mc$allZero %in% TRUE | !is.finite(x) | x <= 0 | !is.finite(disp)
  fit_set <- !zero & disp >= TWOTREND_MIN_DISP
  mix <- twotrend_em(x[fit_set], log10(disp[fit_set]))
  post1 <- rep(mix$pi[1], nrow(raw))
  post1[fit_set] <- twotrend_posterior(mix, x[fit_set], log10(disp[fit_set]))
  # Which trend each PAS fits better: the likelihood ratio alone (equal weights
  # for the two trends, i.e. without the mixing proportions); 0.5 = no information.
  fit1 <- rep(0.5, nrow(raw))
  fit1[fit_set] <- twotrend_posterior(utils::modifyList(mix, list(pi = c(0.5, 0.5))), x[fit_set], log10(disp[fit_set]))
  fit <- list(comparison = comparison_name, heldout_samples = ho, tested_arms = tested_arms,
              feature_id = rownames(raw), mean = x, disp = disp, fit_set = fit_set, mix = mix,
              post1 = stats::setNames(post1, rownames(raw)), fit1 = stats::setNames(fit1, rownames(raw)),
              from_cache = FALSE, method_first = EVS_PCA_METHOD)
  if (nzchar(f)) saveRDS(fit, f)
  fit
}

# S(k) for every k along the comparison's union depth, for one labelling:
# q = likelihood-only probability of fitting the Leading-Edge trend better, so
# S(k) - S(0) = sum over LE(k) of the log likelihood ratio (LE vs REM trend);
# p_post = posterior (with the mixing proportions) for the expected share of
# PASs on their own trend.
twotrend_score_scan <- function(q_le, depth, k_max, p_post) {
  q <- pmin(pmax(q_le, 1e-6), 1 - 1e-6)
  o <- order(depth)
  dd <- depth[o]
  cs <- cumsum(log(q[o]) - log1p(-q[o]))
  ends <- c(which(diff(dd) > 0), length(dd))
  base <- sum(log1p(-q))
  n <- length(q)
  out <- data.frame(k = dd[ends], LE_n = ends, Remainder_n = n - ends, Score = cs[ends] + base,
                    Share_Fit_LE_Better_in_LE = cumsum(q[o] > 0.5)[ends] / ends,
                    Expected_Correct = (cumsum(p_post[o])[ends] + (sum(1 - p_post) - cumsum(1 - p_post[o])[ends])) / n)
  out[out$k <= k_max, , drop = FALSE]
}

# Expected share of PASs on their own trend for a given Leading Edge.
twotrend_expected_correct <- function(p_le, in_le) (sum(p_le[in_le]) + sum(1 - p_le[!in_le])) / length(p_le)

# k* of every comparison for one PCA method (after compute_regime_cutoffs()).
compute_twotrend_cutoffs <- function(method = EVS_PCA_METHOD) {
  if (is.null(REGIME)) stop("compute_regime_cutoffs() must run first (it provides the PC1 ranking).")
  raw <- experiment_counts()
  fits <- list(); scans <- list(); post <- list(); summary <- list(); le_sites <- list()
  for (i in seq_len(nrow(comparison_table))) {
    cmp <- comparison_table$comparison_name[i]
    g <- c(comparison_table$group1_prefix[i], comparison_table$group2_prefix[i])
    fit <- twotrend_heldout_fit(cmp, g)
    tested <- colnames(raw)[as.character(experiment_arm(colnames(raw))) %in% g]
    ids <- rownames(raw)[rowSums(raw[, tested, drop = FALSE]) > 0]
    depth <- regime_depths(cmp, ids)$depth
    p1 <- fit$post1[ids]; q1 <- fit$fit1[ids]
    k_max <- max(1L, floor(length(ids) / 2))
    cand <- lapply(1:2, function(c) {
      p <- if (c == 1L) p1 else 1 - p1
      q <- if (c == 1L) q1 else 1 - q1
      s <- twotrend_score_scan(q, depth, k_max, p)
      s <- s[s$LE_n >= TWOTREND_MIN_SET & s$Remainder_n >= TWOTREND_MIN_SET, , drop = FALSE]
      if (!nrow(s)) stop(cmp, ": no k leaves both sets with at least ", TWOTREND_MIN_SET, " PASs.")
      list(c = c, p = p, q = q, scan = s, best = max(s$Score))
    })
    pick <- cand[[which.max(vapply(cand, `[[`, numeric(1), "best"))]]
    s <- pick$scan; p_le <- pick$p; q_le <- pick$q
    j <- which(s$Score == max(s$Score))[1]
    k <- as.integer(s$k[j])
    le_c <- pick$c; rem_c <- 3L - le_c
    OPTIMAL_K[[cmp]] <<- k
    s$Chosen <- seq_len(nrow(s)) == j
    kr <- REGIME_K[[cmp]]
    ec <- function(kk) twotrend_expected_correct(p_le, depth <= kk)
    par_le <- fit$mix$pars[[le_c]]; par_rem <- fit$mix$pars[[rem_c]]
    fits[[cmp]] <- data.frame(
      Comparison = cmp, Heldout_Samples = paste(fit$heldout_samples, collapse = ","),
      n_PAS_fitted = sum(fit$fit_set), EM_iterations = fit$mix$iterations, EM_converged = fit$mix$converged,
      LogLik = fit$mix$loglik,
      LE_a0 = exp(par_le[1]), LE_a1 = exp(par_le[2]), LE_b = par_le[3], LE_sd = sqrt(fit$mix$s2[le_c]),
      LE_share = fit$mix$pi[le_c],
      REM_a0 = exp(par_rem[1]), REM_a1 = exp(par_rem[2]), REM_b = par_rem[3], REM_sd = sqrt(fit$mix$s2[rem_c]),
      REM_share = fit$mix$pi[rem_c])
    scans[[cmp]] <- cbind(Comparison = cmp, EVS_PCA_Method = method, s)
    post[[cmp]] <- data.frame(Comparison = cmp, feature_id = ids, Heldout_Mean = fit$mean[match(ids, fit$feature_id)],
                              Heldout_Dispersion = fit$disp[match(ids, fit$feature_id)],
                              In_Trend_Fit = fit$fit_set[match(ids, fit$feature_id)],
                              P_LE_Trend = unname(p_le), Fit_P_LE_Trend = unname(q_le),
                              Fits_LE_Trend_Better = unname(q_le > 0.5),
                              Union_Depth = unname(depth), In_Leading_Edge = unname(depth <= k))
    le_sites[[cmp]] <- data.frame(Comparison = cmp, EVS_PCA_Method = method, feature_id = ids[depth <= k])
    summary[[cmp]] <- data.frame(
      Comparison = cmp, EVS_PCA_Method = method, k = k, Leading_Edge = sum(depth <= k), Remainder = sum(depth > k),
      Score = s$Score[j], Expected_Correct_k = ec(k),
      Regime_k = if (is.null(kr)) NA_integer_ else kr,
      Expected_Correct_Regime_k = if (is.null(kr)) NA_real_ else ec(kr),
      Expected_Correct_k5000 = ec(TWOTREND_K_REFERENCE),
      Expected_Correct_All_Remainder = mean(1 - p_le),
      LE_trend_b = par_le[3], REM_trend_b = par_rem[3], LE_trend_share = fit$mix$pi[le_c],
      Mean_P_LE_in_LE = mean(p_le[depth <= k]), Mean_P_LE_in_REM = mean(p_le[depth > k]),
      Share_Fit_LE_Better_in_LE = mean(q_le[depth <= k] > 0.5), Share_Fit_LE_Better_in_REM = mean(q_le[depth > k] > 0.5),
      Share_Fit_LE_Better_All = mean(q_le > 0.5),
      EM_converged = fit$mix$converged)
    message(sprintf("  %s: two-trend k* = %d (Leading Edge %d, %.0f%% fit the LE trend better vs %.0f%% of all PASs; LE trend b = %.2f, REM trend b = %.2f%s)",
                    cmp, k, sum(depth <= k), 100 * mean(q_le[depth <= k] > 0.5), 100 * mean(q_le > 0.5), par_le[3], par_rem[3],
                    if (fit$from_cache) "; held-out fit reused" else ""))
  }
  TWOTREND <<- list(method = method, fits = dplyr::bind_rows(fits), scans = scans, posteriors = post,
                    summary = dplyr::bind_rows(summary), le_sites = dplyr::bind_rows(le_sites))
  invisible(TWOTREND)
}

# ---- two-trend outputs (tables and Figure 1 of this run) ----------------------

save_twotrend_outputs <- function() {
  if (is.null(TWOTREND)) return(invisible(FALSE))
  out <- summary_table_dir
  save_csv(TWOTREND$summary, file.path(out, "Table_TwoTrend_Cutoff_Summary.csv"))
  save_csv(TWOTREND$fits, file.path(out, "Table_TwoTrend_Heldout_Fits.csv"))
  save_csv(dplyr::bind_rows(TWOTREND$scans), file.path(out, "Table_TwoTrend_Score_Scans.csv"))
  save_csv(dplyr::bind_rows(TWOTREND$posteriors), file.path(out, "Table_TwoTrend_PAS_Posteriors.csv"))
  save_csv(TWOTREND$le_sites, file.path(out, "Table_TwoTrend_Leading_Edge_Sites.csv"))
  try_fig(save_figure(plot_twotrend_figure(), paper_fig_dir, paste0("Figure_1_TwoTrend_Cutoff_", TWOTREND$method),
                      180, 165), "Figure 1 (two-trend cutoff)")
  invisible(TRUE)
}

TWOTREND_COLS <- c("Leading-Edge trend" = "#1F6FB2", "Remainder trend" = "#C0392B")
TWOTREND_KCOLS <- c("Two-trend k*" = "#CC79A7", "Regime k*" = "#D55E00", "k = 5,000" = "grey45")

# Figure 1: (a) held-out dispersions and the two trends, (b) share of PASs on
# the Leading-Edge trend along the PC1 union depth, (c) the allocation score
# S(k) with k*, the regime k* and k = 5,000.
plot_twotrend_figure <- function() {
  post <- dplyr::bind_rows(TWOTREND$posteriors)
  fits <- TWOTREND$fits; summ <- TWOTREND$summary
  post$Comparison <- as_comp(post$Comparison); fits$Comparison <- as_comp(fits$Comparison)
  summ$Comparison <- as_comp(summ$Comparison)
  pts <- post[post$In_Trend_Fit %in% TRUE, , drop = FALSE]
  set.seed(20261122L)
  pts <- do.call(rbind, lapply(split(pts, pts$Comparison, drop = TRUE), function(z)
    z[sort(sample.int(nrow(z), min(nrow(z), 12000L))), , drop = FALSE]))
  pts <- pts[order(pts$Fit_P_LE_Trend), , drop = FALSE]
  curves <- do.call(rbind, lapply(seq_len(nrow(fits)), function(i) {
    xr <- range(post$Heldout_Mean[post$Comparison == fits$Comparison[i] & post$In_Trend_Fit %in% TRUE])
    xs <- exp(seq(log(max(xr[1], 0.1)), log(xr[2]), length.out = 200))
    rbind(data.frame(Comparison = fits$Comparison[i], x = xs, Trend = "Leading-Edge trend",
                     y = fits$LE_a0[i] + fits$LE_a1[i] * xs^(-fits$LE_b[i])),
          data.frame(Comparison = fits$Comparison[i], x = xs, Trend = "Remainder trend",
                     y = fits$REM_a0[i] + fits$REM_a1[i] * xs^(-fits$REM_b[i])))
  }))
  blab <- data.frame(Comparison = fits$Comparison,
                     lab = sprintf("LE: b = %.2f, %.0f%%\nREM: b = %.2f, %.0f%%", fits$LE_b, 100 * fits$LE_share,
                                   fits$REM_b, 100 * fits$REM_share))
  pa <- ggplot(pts, aes(Heldout_Mean, Heldout_Dispersion)) +
    geom_point(aes(colour = Fit_P_LE_Trend), size = 0.05, alpha = 0.45, stroke = 0) +
    geom_line(data = curves[curves$Trend == "Leading-Edge trend", ], aes(x, y), colour = TWOTREND_COLS[[1]],
              linewidth = 0.7) +
    geom_line(data = curves[curves$Trend == "Remainder trend", ], aes(x, y), colour = TWOTREND_COLS[[2]],
              linewidth = 0.7, linetype = "22") +
    geom_text(data = blab, aes(label = lab), x = Inf, y = Inf, hjust = 1.04, vjust = 1.15,
              size = 1.9, lineheight = 0.9, inherit.aes = FALSE) +
    facet_wrap(~ Comparison, nrow = 1, labeller = labeller(Comparison = comp_lab())) +
    scale_x_log10(labels = scales::label_number(drop0trailing = TRUE, big.mark = ",")) +
    scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
    scale_colour_gradient2(low = "#C0392B", mid = "grey85", high = "#1F6FB2", midpoint = 0.5, limits = c(0, 1),
                           name = "Fits LE trend\n(likelihood share)") +
    labs(x = "Mean, held-out samples", y = "Gene-wise dispersion,\nheld-out samples", tag = "a",
         title = "Two dispersion trends fitted on the 30 samples not in the comparison",
         subtitle = "Solid blue: Leading-Edge trend; dashed red: Remainder trend; up to 12,000 PASs shown per panel.") +
    fig_theme() + theme(legend.position = "right", legend.key.height = unit(3, "mm"))
  bins <- do.call(rbind, lapply(split(post, post$Comparison, drop = TRUE), function(z) {
    br <- unique(round(exp(seq(0, log(max(z$Union_Depth)), length.out = 26))))
    b <- cut(z$Union_Depth, c(0, br), labels = FALSE)
    z$Fits <- as.numeric(z$Fits_LE_Trend_Better)
    a <- stats::aggregate(cbind(Fits, Union_Depth) ~ b, data = data.frame(b = b, z), FUN = mean)
    a$n <- as.integer(table(b)[as.character(a$b)])
    data.frame(Comparison = z$Comparison[1], a[a$n >= 20, , drop = FALSE])
  }))
  klines <- rbind(data.frame(Comparison = summ$Comparison, k = summ$k, Cutoff = "Two-trend k*"),
                  data.frame(Comparison = summ$Comparison, k = summ$Regime_k, Cutoff = "Regime k*"),
                  data.frame(Comparison = summ$Comparison, k = TWOTREND_K_REFERENCE, Cutoff = "k = 5,000"))
  klines <- klines[is.finite(klines$k), , drop = FALSE]
  klines$Cutoff <- factor(klines$Cutoff, levels = names(TWOTREND_KCOLS))
  share <- data.frame(Comparison = summ$Comparison, s = summ$Share_Fit_LE_Better_All)
  pb <- ggplot(bins, aes(Union_Depth, Fits)) +
    geom_hline(yintercept = 0.5, colour = "grey40", linetype = "22", linewidth = 0.3) +
    geom_hline(data = share, aes(yintercept = s), colour = "grey60", linetype = "dotted", linewidth = 0.3) +
    geom_vline(data = klines, aes(xintercept = k, colour = Cutoff), linewidth = 0.45, linetype = "22") +
    geom_line(linewidth = 0.5) + geom_point(size = 0.7) +
    facet_wrap(~ Comparison, nrow = 1, labeller = labeller(Comparison = comp_lab())) +
    scale_x_log10(labels = scales::label_comma()) + scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = TWOTREND_KCOLS, name = NULL) +
    labs(x = "PC1 union depth (log scale)", y = "Share fitting the LE trend better", tag = "b",
         title = "Top-|PC1| PASs fit the Leading-Edge trend better; the share falls with depth",
         subtitle = "Dotted: share among all PASs of the comparison; dashed: 0.5.") +
    fig_theme() + theme(legend.position = "right")
  sc <- dplyr::bind_rows(TWOTREND$scans)
  sc$Comparison <- as_comp(sc$Comparison)
  sc <- do.call(rbind, lapply(split(sc, sc$Comparison, drop = TRUE), function(z) {
    z$Score_rel <- z$Score - max(z$Score)
    z[unique(round(seq(1, nrow(z), length.out = min(nrow(z), 3000L)))), , drop = FALSE]
  }))
  klab <- data.frame(Comparison = summ$Comparison, k = summ$k,
                     lab = sprintf("k* = %s\nLE %s PASs", format(summ$k, big.mark = ","),
                                   format(summ$Leading_Edge, big.mark = ",")))
  pc <- ggplot(sc, aes(k, Score_rel)) +
    geom_vline(data = klines, aes(xintercept = k, colour = Cutoff), linewidth = 0.45, linetype = "22") +
    geom_line(linewidth = 0.5) +
    geom_text(data = klab, aes(x = k, y = -Inf, label = lab), hjust = -0.08, vjust = -0.3, size = 1.9,
              colour = TWOTREND_KCOLS[["Two-trend k*"]], fontface = "bold") +
    facet_wrap(~ Comparison, nrow = 1, scales = "free_y", labeller = labeller(Comparison = comp_lab())) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_colour_manual(values = TWOTREND_KCOLS, name = NULL) +
    labs(x = "k (top-|PC1| PASs per arm)", y = "S(k) - max S", tag = "c",
         title = "Allocation score: k* places the split where it best matches the two trends",
         subtitle = "S(k) = sum over the Leading Edge of the log likelihood ratio, Leading-Edge vs Remainder trend.") +
    fig_theme() + theme(legend.position = "right")
  assemble_rows(list(list(plots = list(pa), widths = 1), list(plots = list(pb), widths = 1),
                     list(plots = list(pc), widths = 1)), heights = c(1.15, 1, 1),
                title = paste0("Figure 1. Two-trend cutoff k* (", TWOTREND$method, "-EVS ranking)"))
}

# =============================================================================
# VALIDATION CHECKS AND DISCOVERY COUNTS
# =============================================================================

NORMALIZATION_MODE <- "shared"   # size factors from the full comparison matrix, in every view

add_validation <- function(check, comparison, passed, detail = "") {
  validation_rows[[length(validation_rows) + 1L]] <<- data.frame(
    Check = check, Comparison = comparison, Cutoff_Method = CUTOFF_METHOD_LABEL, EVS_PCA_Method = EVS_PCA_METHOD,
    Status = if (isTRUE(passed)) "PASS" else if (is.na(passed)) "NOT APPLICABLE" else "FAIL",
    Detail = detail)
  if (identical(passed, FALSE)) warning("Validation failed (", check, ", ", comparison, "): ", detail, call. = FALSE)
  invisible(passed)
}

size_factors_match <- function(dds, size_factors) {
  if (is.null(size_factors)) return(NA)
  isTRUE(all.equal(unname(sizeFactors(dds)), unname(size_factors[colnames(dds)]), tolerance = 1e-12))
}

# Per-PAS calls of one analysis (a named list of view results): Standard,
# Strong and Weak from one BH over all its views' raw p-values (as in
# call_significance), HBFSS as calibrated per view.
discovery_flags <- function(res_list) {
  cols <- c("feature_id", "pvalue", "resGA_pvalue", "resLA_pvalue", "lfc_shrunk",
            "standard_flag", "strong_cnh_flag", "weak_significant_flag", "hbfss_flag")
  d <- do.call(rbind, lapply(names(res_list), function(v) {
    r <- res_list[[v]]
    miss <- setdiff(cols, names(r))
    if (length(miss)) stop("Result table ", v, " lacks: ", paste(miss, collapse = ", "))
    data.frame(View = v, r[, cols, drop = FALSE], stringsAsFactors = FALSE)
  }))
  if (anyDuplicated(d$feature_id)) stop("A PAS appears in more than one view of the same analysis.")
  bh <- function(p) {
    q <- rep(NA_real_, length(p))
    ok <- is.finite(p)
    q[ok] <- stats::p.adjust(p[ok], method = "BH")
    q
  }
  lfc_ok <- is.finite(d$lfc_shrunk)
  big <- lfc_ok & abs(d$lfc_shrunk) >= lfc_boundary
  small <- lfc_ok & abs(d$lfc_shrunk) < lfc_boundary
  d$padj_global <- bh(d$pvalue)
  d$padjGA_global <- bh(d$resGA_pvalue)
  d$padjLA_global <- bh(d$resLA_pvalue)
  d$HBFSS <- d$hbfss_flag %in% TRUE
  d$Standard <- big & !is.na(d$padj_global) & d$padj_global < BH_FDR_STANDARD
  d$Strong <- big & !is.na(d$padjGA_global) & d$padjGA_global < BH_FDR_STRONG
  d$Weak <- small & !is.na(d$padjLA_global) & d$padjLA_global < BH_FDR_WEAK & d$HBFSS
  d
}

count_discoveries <- function(flags) {
  c(Standard = sum(flags$Standard), Strong = sum(flags$Strong), Weak = sum(flags$Weak), HBFSS = sum(flags$HBFSS),
    n_tested = nrow(flags), n_with_pvalue = sum(is.finite(flags$pvalue)))
}

# Shared normalization reached every view, the samples match, the partition is
# exact, and a PAS has identical normalized counts in the
# Original and in its EVS partition.
validate_track_normalization <- function(comparison_name, track_key, analysis_results, sf, evs, coldata) {
  tr <- unname(track_short[track_key])
  ids <- rownames(evs$raw_dataset)
  add_validation("Every PAS is assigned exactly once (LE / REM)", comparison_name,
                 !length(intersect(evs$leading_edge_ids, evs$remainder_ids)) &&
                   setequal(c(evs$leading_edge_ids, evs$remainder_ids), ids) &&
                   length(evs$leading_edge_ids) + length(evs$remainder_ids) == length(ids), tr)
  same_samples <- all(vapply(analysis_results, function(a) identical(colnames(a$dds), rownames(coldata)), logical(1)))
  add_validation("Original and EVS sample identities match", comparison_name, same_samples, tr)
  ev <- c("leading_edge_dataset", "remainder_dataset")
  oc <- original_cooks_cache[[comparison_name]]
  cooks_ok <- !is.null(oc) && all(vapply(ev, function(v) {
    r <- analysis_results[[v]]$results
    out <- r$cooks_outlier %in% TRUE
    setequal(r$feature_id[out], intersect(oc, r$feature_id)) && all(is.na(r$pvalue[out]))
  }, logical(1)))
  add_validation("Cook's outliers are the Original's in every view", comparison_name, cooks_ok,
                 sprintf("%s; %d outlier PASs in the Original", tr, length(oc)))
  ab <- do.call(rbind, lapply(ev, function(v) analysis_results[[v]]$results[, c("pvalue", "padj")]))
  q <- rep(NA_real_, nrow(ab))
  okp <- is.finite(ab$pvalue)
  q[okp] <- stats::p.adjust(ab$pvalue[okp], method = "BH")
  add_validation("Standard padj is one BH over LE + REM, without independent filtering", comparison_name,
                 isTRUE(all.equal(q, ab$padj)), tr)
  add_validation("Shared normalization is propagated", comparison_name,
                 all(vapply(analysis_results, function(a) isTRUE(size_factors_match(a$dds, sf)), logical(1))), tr)
  o0 <- analysis_results$raw_dataset %||% {
    reg0 <- paper_registry[[registry_key(comparison_name, "normalized_evs")]]
    if (is.null(reg0)) NULL else reg0$analysis_results$raw_dataset
  }
  if (!is.null(o0)) {
    og <- S4Vectors::mcols(o0$dds)$dispGeneEst
    names(og) <- rownames(o0$dds)
    same_gw <- all(vapply(c("leading_edge_dataset", "remainder_dataset"), function(v) {
      d <- analysis_results[[v]]$dds
      isTRUE(all.equal(unname(S4Vectors::mcols(d)$dispGeneEst), unname(og[rownames(d)]), tolerance = 1e-8))
    }, logical(1)))
    add_validation("Gene-wise dispersions are identical in the Original and the partition (computed once per comparison)",
                   comparison_name, same_gw, tr)
  }
  o <- analysis_results$raw_dataset %||% {
    reg <- paper_registry[[registry_key(comparison_name, "normalized_evs")]]
    if (is.null(reg)) NULL else reg$analysis_results$raw_dataset
  }
  if (!is.null(o)) {
    on <- DESeq2::counts(o$dds, normalized = TRUE)
    same <- all(vapply(c("leading_edge_dataset", "remainder_dataset"), function(v) {
      vn <- DESeq2::counts(analysis_results[[v]]$dds, normalized = TRUE)
      isTRUE(all.equal(unname(vn), unname(on[rownames(vn), colnames(vn), drop = FALSE]), tolerance = 1e-10))
    }, logical(1)))
    add_validation("Normalized counts are identical in the Original and the EVS partition", comparison_name, same, tr)
  }
  invisible(TRUE)
}

# The chosen k* is the one used: the split has k = k*, and its Leading Edge is
# exactly the set the two-trend score evaluated at k*.
validate_selected_k_used <- function(comparison_name) {
  k <- OPTIMAL_K[[comparison_name]]
  reg <- paper_registry[[registry_key(comparison_name, "normalized_evs")]]
  cache <- comparison_cache[[comparison_name]]
  srow <- if (is.null(TWOTREND)) NULL else TWOTREND$summary[TWOTREND$summary$Comparison == comparison_name, , drop = FALSE]
  sites <- if (is.null(TWOTREND)) NULL else TWOTREND$le_sites$feature_id[TWOTREND$le_sites$Comparison == comparison_name]
  if (is.null(reg) || is.null(cache) || is.null(srow) || nrow(srow) != 1L) {
    add_validation("The selected k* is used throughout the analysis", comparison_name, FALSE, "missing cutoff fit or analysis")
    return(invisible(FALSE))
  }
  evs <- reg$evs
  ks <- c(selected = k, split = as.integer(evs$empirical_evs_k), active = get_active_evs_cutoff(comparison_name),
          table = as.integer(srow$k))
  add_validation("The selected k* is used throughout the analysis", comparison_name,
                 length(unique(ks)) == 1L, paste(names(ks), ks, sep = " = ", collapse = "; "))
  add_validation("The analysed Leading Edge is the set scored at k*", comparison_name,
                 setequal(evs$leading_edge_ids, sites) &&
                   setequal(sites, rownames(cache$raw)[cache$depth[rownames(cache$raw)] <= k]),
                 sprintf("Leading Edge %d; scored set %d", length(evs$leading_edge_ids), length(sites)))
  invisible(TRUE)
}

# PCA preprocessing checks for this run's method, plus the ranking table.
validate_pca_preprocessing <- function(comparison_name) {
  cc <- comparison_cache[[comparison_name]]
  raw <- cc$raw; coldata <- cc$coldata
  add_validation("PAS IDs, sample order, raw counts and size factors identical across PCA methods",
                 comparison_name, cc$shared_check, cc$shared_detail)
  ro <- make_rank_matrix_for_track(raw, coldata, "normalized_evs", cache_key = comparison_name)
  x <- as.matrix(ro$rank_matrix)
  full <- experiment_rank_matrix(ro$method)
  add_validation("PCA input is finite, matches the raw PAS IDs and samples, and equals the experiment-wide matrix",
                 comparison_name,
                 all(is.finite(x)) && identical(rownames(x), rownames(raw)) && identical(colnames(x), colnames(raw)) &&
                   identical(unname(x), unname(full[rownames(raw), colnames(raw), drop = FALSE])),
                 sprintf("%s: %d PASs x %d samples; range %.3g to %.3g", ro$method, nrow(x), ncol(x),
                         min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
  fits <- pc1_condition_fits(ro, coldata, comparison_name)
  for (cond in c("trt", "untrt")) {
    v <- fits[[cond]]$pca_fit$sdev^2
    tot <- sum(v)
    ok <- is.finite(tot) && tot > 1e-12 && v[1] > 0 && (length(v) < 2L || v[1] > v[2] * (1 + 1e-8))
    add_validation("PCA has sufficient variance (PC1 nonzero and distinct from PC2)", comparison_name, ok,
                   sprintf("%s %s: PC1 %.1f%% of variance, PC1/PC2 variance ratio %.2f", ro$method,
                           if (cond == "trt") "RT" else "ZT", 100 * v[1] / tot,
                           if (length(v) > 1L) v[1] / v[2] else NA_real_))
  }
  rk <- lapply(fits, function(f) stats::setNames(f$loading_table$rank, f$loading_table$feature_id)[rownames(raw)])
  k <- OPTIMAL_K[[comparison_name]]
  same_top <- setequal(names(rk$trt)[rk$trt <= k], names(cc$rank_RT)[cc$rank_RT <= k]) &&
    setequal(names(rk$untrt)[rk$untrt <= k], names(cc$rank_ZT)[cc$rank_ZT <= k])
  add_validation("Comparison PC1 top-k* sets equal the regime rule's arm top-k* sets", comparison_name, same_top,
                 sprintf("%s, k* = %d", ro$method, k))
  pca_ranking_rows[[comparison_name]] <<- data.frame(
    Comparison = comparison_name, EVS_PCA_Method = ro$method, feature_id = rownames(raw),
    Rank_RT = unname(rk$trt), Rank_ZT = unname(rk$untrt), Union_Depth = unname(cc$depth[rownames(raw)]))
  invisible(TRUE)
}

# Validation table, the DESeq2 default-trend baseline and the PC1 rankings
# (end of run).
save_run_checks <- function() {
  out <- summary_table_dir
  if (length(validation_rows)) save_csv(dplyr::bind_rows(validation_rows), file.path(out, "Table_Validation_Checks.csv"))
  if (length(default_trend_rows)) save_csv(dplyr::bind_rows(default_trend_rows), file.path(out, "Table_Original_DESeq2_Default_Trend.csv"))
  if (length(pca_ranking_rows)) save_csv(dplyr::bind_rows(pca_ranking_rows), file.path(out, "Table_EVS_PCA_Ranking.csv"))
  invisible(TRUE)
}


run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  dir.create(file.path(output_dir, comparison_name), recursive = TRUE, showWarnings = FALSE)

  # Stage 1: k* was chosen by compute_regime_cutoffs() (no test of any kind);
  # comparison inputs reused by the analysis are prepared here.
  calls_before <- c(SIGNIFICANCE_CALLS$wald, SIGNIFICANCE_CALLS$null)
  prepare_comparison_cache(comparison_name, count_matrix, coldata)
  calls_after <- c(SIGNIFICANCE_CALLS$wald, SIGNIFICANCE_CALLS$null)
  add_validation("Comparison preparation performed no significance testing", comparison_name,
                 identical(calls_before, calls_after),
                 sprintf("Wald fits %d -> %d; significance calls %d -> %d",
                         calls_before[1], calls_after[1], calls_before[2], calls_after[2]))
  ts <- TWOTREND$summary[TWOTREND$summary$Comparison == comparison_name, , drop = FALSE]
  tf <- TWOTREND$fits[TWOTREND$fits$Comparison == comparison_name, , drop = FALSE]
  le_n <- sum(comparison_cache[[comparison_name]]$depth <= get_active_evs_cutoff(comparison_name))
  active_cutoff_rows[[comparison_name]] <<- data.frame(
    Comparison = comparison_name, Cutoff_Method = CUTOFF_METHOD_LABEL, k = get_active_evs_cutoff(comparison_name),
    Leading_Edge = le_n, Remainder = nrow(comparison_cache[[comparison_name]]$raw) - le_n,
    Expected_On_Trend = ts$Expected_Correct_k, Expected_On_Trend_All_Remainder = ts$Expected_Correct_All_Remainder,
    Regime_k = ts$Regime_k, Expected_On_Trend_Regime_k = ts$Expected_Correct_Regime_k,
    Expected_On_Trend_k5000 = ts$Expected_Correct_k5000, LE_trend_b = ts$LE_trend_b, REM_trend_b = ts$REM_trend_b,
    LE_trend_share = ts$LE_trend_share, EM_converged = tf$EM_converged,
    EVS_PCA_Method = EVS_PCA_METHOD, PCA_Preprocessing = PCA_METHOD_LABELS[[EVS_PCA_METHOD]],
    Normalization_Mode = NORMALIZATION_MODE)
  add_validation("Two-trend fit used only samples outside the comparison", comparison_name,
                 !any(strsplit(tf$Heldout_Samples, ",")[[1]] %in% colnames(comparison_cache[[comparison_name]]$raw)) &&
                   length(strsplit(tf$Heldout_Samples, ",")[[1]]) == 30L,
                 sprintf("held-out samples: %s", tf$Heldout_Samples))
  add_validation("Two-trend mixture (EM) converged", comparison_name, isTRUE(tf$EM_converged),
                 sprintf("%d iterations; LE b = %.2f, REM b = %.2f", tf$EM_iterations, tf$LE_b, tf$REM_b))
  validate_pca_preprocessing(comparison_name)

  # Stage 2: one significance analysis at the selected k (Original; LE + REM).
  message("Running ", comparison_name, " NormEVS at k* = ", OPTIMAL_K[[comparison_name]])
  ok <- run_one_evs_track(comparison_name = comparison_name, track_key = "normalized_evs",
                          count_matrix = count_matrix, coldata = coldata, annot_df = annot_df)
  if (!isTRUE(ok)) stop("Analysis failed: ", comparison_name)
  validate_selected_k_used(comparison_name)
  add_validation("Exactly one k* selected for this comparison and PCA method", comparison_name,
                 sum(TWOTREND$scans[[comparison_name]]$Chosen) == 1L && length(OPTIMAL_K[[comparison_name]]) == 1L,
                 sprintf("%s: k* = %d", EVS_PCA_METHOD, OPTIMAL_K[[comparison_name]]))

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
    "NormEVS Rem" = "N-Rem"
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
    for (track_key in "normalized_evs") {
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

  for (track_key in "normalized_evs") {
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

    for (track_key in "normalized_evs") {
      obj <- paper_registry[[registry_key(comparison_name, track_key)]]
      if (is.null(obj) || is.null(obj$evs)) next

      original_mat <- pca_track_matrix(obj, "raw_dataset")
      original_stats <- pca_variance_stats(original_mat, obj$coldata)
      if (is.null(original_stats)) next


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
          short_title = ds_abbr,
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
      base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_PCA_EVS"))
      save_grob(panel, paste0(base, ".png"), width = 16.5, height = 6.0)
      save_grob(panel, paste0(base, ".pdf"), width = 16.5, height = 6.0)
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
      track_order_use <- "normalized_evs"
      panel_title <- paste0("HBFSS calibration | ", unname(dataset_short[dataset_key]), " | ", CUTOFF_METHOD_LABEL)
      output_suffix <- paste0(unname(dataset_short[dataset_key]), "_AllComparisons")
      panel_height <- 5.8
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
        "NormEVS Lead", "NormEVS Rem"
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
        "NormEVS Lead", "NormEVS Rem"
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
        "NormEVS Lead", "NormEVS Rem"
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
    "NormEVS Rem"="N-Rem"
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

    preferred_views <- c("Original (No EVS)", "NormEVS Lead", "NormEVS Rem")
    views <- preferred_views[preferred_views %in% unique(as.character(sub$Analysis))]
    plots <- lapply(views, function(v) {
      one <- sub[sub$Analysis == v, , drop = FALSE]
      if (!nrow(one)) return(NULL)
      short <- c(
        "Original (No EVS)"="Orig",
        "NormEVS Lead"="N-Lead",
        "NormEVS Rem"="N-Rem"
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

# Cutoffs for every comparison, before any test: the PC1 ranking and the
# regime reference cutoff (compute_regime_cutoffs), then the two held-out
# dispersion trends and k* (compute_twotrend_cutoffs).
regime_calls_before <- c(SIGNIFICANCE_CALLS$wald, SIGNIFICANCE_CALLS$null)
compute_regime_cutoffs(EVS_PCA_METHOD)
compute_twotrend_cutoffs(EVS_PCA_METHOD)
regime_calls_after <- c(SIGNIFICANCE_CALLS$wald, SIGNIFICANCE_CALLS$null)
add_validation("Cutoff selection performed no significance testing", "All comparisons",
               identical(regime_calls_before, regime_calls_after),
               sprintf("Wald fits %d -> %d; significance calls %d -> %d", regime_calls_before[1],
                       regime_calls_after[1], regime_calls_before[2], regime_calls_after[2]))

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
  overall_summary$EVS_PCA_Method <- EVS_PCA_METHOD
  save_csv(overall_summary, file.path(summary_table_dir, "Table_DE_Method_Counts.csv"))
}

empirical_cutoff_export <- comparison_table[, c("comparison_name", "group1_prefix", "group2_prefix"), drop = FALSE]
names(empirical_cutoff_export) <- c("Comparison", "RT_prefix", "ZT_prefix")
empirical_cutoff_export$Active_Cutoff_Method <- CUTOFF_METHOD_LABEL
empirical_cutoff_export$Active_k_per_condition <- vapply(
  empirical_cutoff_export$Comparison,
  function(cn) tryCatch(get_active_evs_cutoff(cn), error = function(e) NA_integer_), integer(1)
)
empirical_cutoff_export$Basis <- EMPIRICAL_EVS_CUTOFF_BASIS
empirical_cutoff_export$Applied_to <- "NormEVS"
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
tryCatch(save_run_checks(),
         error = function(e) warning("Validation outputs skipped: ", conditionMessage(e), call. = FALSE))

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
