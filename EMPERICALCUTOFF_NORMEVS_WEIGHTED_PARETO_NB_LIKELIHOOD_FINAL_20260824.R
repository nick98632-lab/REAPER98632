#!/usr/bin/env Rscript
# =============================================================================
# EVS_Manuscript_Figures.R
# Publication-quality figure layer for the CPM-EVS vs DESeq2-EVS pipeline.
#
# USAGE
# -----
# In CPM_vs_DESeq2_EVS_MANUSCRIPT_FINAL.R, delete the block running from
#   "# FOUR COMPOSITE MANUSCRIPT FIGURES"
# down to and including the writeLines() call that creates Figure_Legends.md,
# then insert in its place:
#
#   source("EVS_Manuscript_Figures.R")
#   evs_render_all(
#     method_objects = method_objects,
#     summary_df     = summary_df,
#     overlap_df     = overlap_df,
#     comparisons    = COMPARISONS,
#     fig_dir        = FIG_DIR,
#     out_root       = OUT_ROOT
#   )
#
# The old figure helpers (theme_pub, save_composite, make_*_panel) may be
# deleted; nothing here depends on them. The analysis code is untouched.
#
# Check rendering without the count matrix:
#   Rscript EVS_Manuscript_Figures.R --selftest
#
# DEPENDENCIES: ggplot2 (>= 3.4), dplyr, tidyr, grid, scales, grDevices.
# scales and gtable ship with ggplot2, so nothing new needs installing.
# patchwork is used for panel alignment IF it happens to be installed; without
# it the module falls back to a gtable aligner and still runs.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

# =============================================================================
# 1. DESIGN SYSTEM
# =============================================================================
# Figures are authored at final print size. A figure drawn 15 in wide and then
# reduced to a 180 mm journal column shrinks all type by ~60%, which is the
# single most common reason draft figures fail production. Everything below is
# sized for 180 mm (double column) reproduced at 100%.

EVS_FIG_W_MM   <- 180    # double-column width
EVS_BASE_PT    <- 7      # body text, in points, at final size
EVS_PNG_DPI    <- 600
EVS_WRITE_TIFF <- FALSE  # set TRUE if the journal requires TIFF

mm2in <- function(mm) mm / 25.4
pt2mm <- function(pt) pt / .pt   # ggplot geom text size is in mm, not points

# Okabe-Ito: colourblind-safe and still separable in greyscale.
EVS_COL <- list(
  cpm       = "#0072B2",  # blue
  deseq     = "#D55E00",  # vermillion
  median    = "#CC79A7",  # reddish purple
  fit       = "#000000",
  arms      = "#9AA3AD",
  knee      = "#CC79A7",
  chord     = "#B9C0C7",
  shared    = "#7A7F86",
  joint     = "#009E73",
  opp_le    = "#56B4E9",
  opp_div   = "#F0E442",
  opp_rem   = "#E69F00",
  band_rem  = "#E8ECF2",  # pale regime fills ...
  band_div  = "#FBF0D0",
  band_lead = "#DCEFE4",
  edge_rem  = "#5B6672",  # ... with saturated partners for rules and labels
  edge_div  = "#8A6D12",
  edge_lead = "#2E8B62"
)

EVS_METHOD_COL <- c("CPM-EVS" = EVS_COL$cpm, "DESeq2-EVS" = EVS_COL$deseq)

evs_method_label <- function(m) ifelse(m == "CPM_EVS", "CPM-EVS", "DESeq2-EVS")

# "RT0_ZT6" -> "RT0 vs ZT6". Underscores do not belong on a published axis.
evs_comparison_label <- function(x) sub("_", " vs ", x, fixed = TRUE)

evs_num <- function(x, digits = 0) {
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

evs_comma <- function(v) {
  format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
}

# legend.position.inside arrived in ggplot2 3.5.0; older versions take a
# numeric legend.position. Keeps the module working on both.
evs_legend_inside <- function(x, y, just, direction = "vertical") {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside",
          legend.position.inside = c(x, y),
          legend.justification = just,
          legend.direction = direction)
  } else {
    theme(legend.position = c(x, y),
          legend.justification = just,
          legend.direction = direction)
  }
}

evs_theme <- function(base_pt = EVS_BASE_PT) {
  theme_classic(base_size = base_pt) +
    theme(
      plot.title        = element_text(size = base_pt + 0.5, face = "bold",
                                       hjust = 0, margin = margin(b = 2.5)),
      plot.tag          = element_text(size = base_pt + 2, face = "bold"),
      plot.tag.position = "topleft",
      axis.title        = element_text(size = base_pt),
      axis.title.x      = element_text(margin = margin(t = 2.5)),
      axis.title.y      = element_text(margin = margin(r = 2.5)),
      axis.text         = element_text(size = base_pt - 0.5, colour = "grey15"),
      axis.line         = element_line(linewidth = 0.25, colour = "grey20"),
      axis.ticks        = element_line(linewidth = 0.25, colour = "grey20"),
      axis.ticks.length = unit(1.1, "pt"),
      legend.position   = "none",
      legend.title      = element_blank(),
      legend.text       = element_text(size = base_pt - 1),
      legend.key.size   = unit(6, "pt"),
      legend.spacing.x  = unit(2, "pt"),
      legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
      legend.margin     = margin(1, 2, 1, 2),
      panel.background  = element_blank(),
      plot.background   = element_blank(),
      plot.margin       = margin(3, 4, 3, 3)
    )
}

evs_theme_inset <- function() {
  theme_classic(base_size = 5) +
    theme(
      axis.title        = element_text(size = 5, colour = "grey20"),
      axis.title.x      = element_text(margin = margin(t = 0.5)),
      axis.title.y      = element_text(margin = margin(r = 0.5)),
      axis.text         = element_text(size = 4.4, colour = "grey25"),
      axis.line         = element_line(linewidth = 0.2, colour = "grey30"),
      axis.ticks        = element_line(linewidth = 0.2, colour = "grey30"),
      axis.ticks.length = unit(0.8, "pt"),
      legend.position   = "none",
      plot.background   = element_rect(fill = "white", colour = "grey75",
                                       linewidth = 0.2),
      plot.margin       = margin(2, 3, 1, 1)
    )
}

# An 8-arm D(r) plot over ~32,000 ranks carries ~256,000 line vertices.
# Thinning is visually identical at print size and keeps the PDF small enough
# for a submission portal.
evs_thin <- function(n, out = 2000L) {
  if (n <= out) return(seq_len(n))
  unique(c(1L, as.integer(round(seq(1, n, length.out = out))), n))
}

# =============================================================================
# 2. DEVICES AND PANEL COMPOSITION
# =============================================================================

evs_open_device <- function(file, width_in, height_in, kind) {
  cairo_ok <- isTRUE(capabilities("cairo"))
  if (kind == "pdf") {
    if (cairo_ok) grDevices::cairo_pdf(file, width = width_in, height = height_in)
    else          grDevices::pdf(file, width = width_in, height = height_in)
  } else if (kind == "png") {
    if (cairo_ok) {
      grDevices::png(file, width = width_in, height = height_in, units = "in",
                     res = EVS_PNG_DPI, bg = "white", type = "cairo")
    } else {
      grDevices::png(file, width = width_in, height = height_in, units = "in",
                     res = EVS_PNG_DPI, bg = "white")
    }
  } else if (kind == "tiff") {
    if (cairo_ok) {
      grDevices::tiff(file, width = width_in, height = height_in, units = "in",
                      res = EVS_PNG_DPI, bg = "white", compression = "lzw",
                      type = "cairo")
    } else {
      grDevices::tiff(file, width = width_in, height = height_in, units = "in",
                      res = EVS_PNG_DPI, bg = "white", compression = "lzw")
    }
  } else {
    stop("Unknown device: ", kind)
  }
}

# Align panel edges across the grid so axes line up column-wise and row-wise.
# Uses gtable only (ships with ggplot2), so no patchwork/cowplot dependency.
evs_align <- function(plots, nrow, ncol) {
  gs <- lapply(plots, ggplot2::ggplotGrob)

  same_w <- length(unique(vapply(gs, function(g) length(g$widths), integer(1)))) == 1L
  same_h <- length(unique(vapply(gs, function(g) length(g$heights), integer(1)))) == 1L

  if (same_w) {
    for (j in seq_len(ncol)) {
      idx <- which(((seq_along(gs) - 1L) %% ncol) + 1L == j)
      if (length(idx) > 1L) {
        w <- do.call(grid::unit.pmax, lapply(gs[idx], function(g) g$widths))
        for (i in idx) gs[[i]]$widths <- w
      }
    }
  }
  if (same_h) {
    for (r in seq_len(nrow)) {
      idx <- which(((seq_along(gs) - 1L) %/% ncol) + 1L == r)
      if (length(idx) > 1L) {
        h <- do.call(grid::unit.pmax, lapply(gs[idx], function(g) g$heights))
        for (i in idx) gs[[i]]$heights <- h
      }
    }
  }
  if (!same_w || !same_h) {
    message("  note: heterogeneous panel layouts; drawing without full alignment")
  }
  gs
}

evs_draw_grid <- function(grobs, nrow, ncol) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow, ncol)))
  for (i in seq_along(grobs)) {
    r  <- ((i - 1L) %/% ncol) + 1L
    cc <- ((i - 1L) %% ncol) + 1L
    pushViewport(viewport(layout.pos.row = r, layout.pos.col = cc))
    grid.draw(grobs[[i]])
    popViewport()
  }
  popViewport()
}

evs_save <- function(plots, file_base, nrow, ncol,
                     width_mm = EVS_FIG_W_MM, height_mm = 150) {
  # patchwork aligns panels correctly even when some are faceted and some are
  # not (Figure 4C is faceted, A/B/D are not). The gtable path below only
  # aligns gtables of identical structure, so it is the fallback.
  use_pw <- requireNamespace("patchwork", quietly = TRUE)
  grobs  <- if (use_pw) NULL else evs_align(plots, nrow, ncol)

  w <- mm2in(width_mm)
  h <- mm2in(height_mm)
  kinds <- c("pdf", "png", if (isTRUE(EVS_WRITE_TIFF)) "tiff")

  for (k in kinds) {
    f <- paste0(file_base, ".", k)
    evs_open_device(f, w, h, k)
    ok <- tryCatch({
      if (use_pw) {
        print(patchwork::wrap_plots(plots, nrow = nrow, ncol = ncol))
      } else {
        evs_draw_grid(grobs, nrow, ncol)
      }
      TRUE
    }, error = function(e) {
      message("  ERROR drawing ", basename(f), ": ", conditionMessage(e))
      FALSE
    })
    grDevices::dev.off()
    if (ok) message("  wrote ", basename(f))
  }
  invisible(file_base)
}

# =============================================================================
# 3. FIGURE 1 PANELS — regime definition
# =============================================================================

evs_arm_matrix <- function(arms) {
  m <- do.call(cbind, lapply(arms, function(z) z$D))
  if (is.null(colnames(m))) colnames(m) <- paste0("arm", seq_len(ncol(m)))
  m
}

panel_regime <- function(method, arms, knot, tag, show_legend = FALSE) {
  N   <- knot$N
  c1  <- knot$c1
  c2  <- knot$c2
  idx <- evs_thin(N, 1800L)

  Dm <- evs_arm_matrix(arms)[idx, , drop = FALSE]
  rr <- idx

  n_arms    <- ncol(Dm)
  lab_arms  <- sprintf("Individual arms (n = %d)", n_arms)
  lab_med   <- "Median observed"
  lab_fit   <- "Shared two-knot fit"
  key_levels <- c(lab_arms, lab_med, lab_fit)

  arms_long <- data.frame(
    rank  = rep(rr, times = n_arms),
    D     = as.vector(Dm),
    arm   = rep(colnames(Dm), each = length(rr)),
    key   = lab_arms,
    stringsAsFactors = FALSE
  )
  med <- data.frame(rank = rr, D = apply(Dm, 1, median), key = lab_med)
  fit <- data.frame(rank = rr,
                    D = apply(knot$fitted[idx, , drop = FALSE], 1, median),
                    key = lab_fit)

  yr   <- range(c(arms_long$D, fit$D), na.rm = TRUE)
  ypad <- diff(yr)
  if (!is.finite(ypad) || ypad <= 0) ypad <- 1
  ylo <- yr[1] - 0.04 * ypad
  yhi <- yr[2] + 0.32 * ypad          # headroom for the regime labels

  bands <- data.frame(
    xmin = c(1, c1, c2),
    xmax = c(c1, c2, N),
    fill = c(EVS_COL$band_rem, EVS_COL$band_div, EVS_COL$band_lead),
    lab  = c("Remainder", "Divergence", "Leading edge"),
    col  = c(EVS_COL$edge_rem, EVS_COL$edge_div, EVS_COL$edge_lead),
    stringsAsFactors = FALSE
  )
  bands$xmid <- (bands$xmin + bands$xmax) / 2
  bands$y    <- yr[2] + 0.265 * ypad

  knots <- data.frame(
    x   = c(c1, c2),
    y   = yr[2] + 0.135 * ypad,
    col = c(EVS_COL$edge_rem, EVS_COL$edge_lead),
    hj  = c(1.06, -0.06),
    lab = c(sprintf("italic(c)[1] * \" = \" * \"%s\"", evs_num(c1)),
            sprintf("italic(c)[2] * \" = \" * \"%s\"", evs_num(c2))),
    stringsAsFactors = FALSE
  )

  p <- ggplot() +
    annotate("rect", xmin = bands$xmin, xmax = bands$xmax,
             ymin = -Inf, ymax = Inf, fill = bands$fill) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dotted",
               linewidth = 0.2) +
    geom_line(data = arms_long,
              aes(rank, D, group = arm, colour = key, linetype = key,
                  linewidth = key)) +
    geom_line(data = med,
              aes(rank, D, colour = key, linetype = key, linewidth = key)) +
    geom_line(data = fit,
              aes(rank, D, colour = key, linetype = key, linewidth = key)) +
    geom_vline(xintercept = c1, colour = EVS_COL$edge_rem,
               linetype = "22", linewidth = 0.3) +
    geom_vline(xintercept = c2, colour = EVS_COL$edge_lead,
               linetype = "22", linewidth = 0.3) +
    geom_text(data = bands, aes(x = xmid, y = y, label = lab),
              colour = bands$col, size = pt2mm(5.8), vjust = 1,
              inherit.aes = FALSE) +
    geom_text(data = knots, aes(x = x, y = y, label = lab),
              colour = knots$col, size = pt2mm(5.6), hjust = knots$hj,
              vjust = 1, fontface = "bold", parse = TRUE, inherit.aes = FALSE) +
    scale_colour_manual(breaks = key_levels,
                        values = setNames(c(EVS_COL$arms, EVS_COL$median,
                                            EVS_COL$fit), key_levels)) +
    scale_linetype_manual(breaks = key_levels,
                          values = setNames(c("solid", "solid", "22"), key_levels)) +
    scale_linewidth_manual(breaks = key_levels,
                           values = setNames(c(0.15, 0.55, 0.42), key_levels)) +
    scale_x_continuous(labels = evs_comma) +
    coord_cartesian(xlim = c(1, N), ylim = c(ylo, yhi), expand = FALSE) +
    labs(
      tag   = tag,
      title = evs_method_label(method),
      x     = "PAS rank by |PC1 loading| (ascending)",
      y     = expression(italic(D)(italic(r)) ==
                           italic(F)[E](italic(r)) - italic(F)[P](italic(r)))
    ) +
    evs_theme()

  if (show_legend) {
    p <- p +
      guides(colour    = guide_legend(order = 1),
             linetype  = guide_legend(order = 1),
             linewidth = guide_legend(order = 1)) +
      evs_legend_inside(0.02, 0.60, c(0, 1))
  }
  p
}

panel_residuals <- function(method_objects, tag, n_bins = 180L) {
  res <- bind_rows(lapply(names(method_objects), function(m) {
    ob <- method_objects[[m]]
    D  <- evs_arm_matrix(ob$arms)
    R  <- D - ob$knot$fitted
    N  <- nrow(R)
    bin <- cut(seq_len(N), breaks = n_bins, labels = FALSE)
    data.frame(
      rank   = rep(seq_len(N), times = ncol(R)),
      bin    = rep(bin, times = ncol(R)),
      resid  = as.vector(R),
      method = evs_method_label(m),
      stringsAsFactors = FALSE
    ) %>%
      group_by(method, bin) %>%
      summarise(rank = mean(rank),
                lo   = stats::quantile(resid, 0.25, names = FALSE),
                mid  = stats::median(resid),
                hi   = stats::quantile(resid, 0.75, names = FALSE),
                .groups = "drop")
  }))
  res$method <- factor(res$method, levels = names(EVS_METHOD_COL))

  ggplot(res, aes(rank, mid, colour = method, fill = method)) +
    geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.25) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, colour = NA) +
    geom_line(aes(linetype = method), linewidth = 0.4) +
    scale_colour_manual(values = EVS_METHOD_COL) +
    scale_fill_manual(values = EVS_METHOD_COL) +
    scale_linetype_manual(values = c("CPM-EVS" = "solid", "DESeq2-EVS" = "22")) +
    scale_x_continuous(labels = evs_comma) +
    labs(
      tag   = tag,
      title = "Shared-knot fit residuals",
      x     = "PAS rank by |PC1 loading| (ascending)",
      y     = expression(paste("Residual, ",
                               italic(D)[italic(g)](italic(r)) -
                                 hat(italic(D))[italic(g)](italic(r))))
    ) +
    evs_theme() +
    evs_legend_inside(0.98, 0.98, c(1, 1), direction = "horizontal")
}

panel_regime_widths <- function(method_objects, tag) {
  d <- bind_rows(lapply(seq_along(method_objects), function(i) {
    m <- names(method_objects)[i]
    k <- method_objects[[m]]$knot
    data.frame(
      method = evs_method_label(m),
      y      = length(method_objects) - i + 1,
      xmin   = c(1, k$c1, k$c2),
      xmax   = c(k$c1, k$c2, k$N),
      regime = c("Remainder", "Divergence", "Leading edge"),
      n      = c(k$c1 - 1L, k$c2 - k$c1 + 1L, k$N - k$c2),
      c1     = k$c1, c2 = k$c2, N = k$N,
      stringsAsFactors = FALSE
    )
  }))
  d$regime <- factor(d$regime, levels = c("Remainder", "Divergence", "Leading edge"))
  d$xmid   <- (d$xmin + d$xmax) / 2

  kn <- d %>%
    distinct(method, y, c1, c2) %>%
    tidyr::pivot_longer(c("c1", "c2"), names_to = "which", values_to = "x") %>%
    mutate(
      lab = ifelse(which == "c1",
                   sprintf("italic(c)[1] * \" = \" * \"%s\"", evs_num(x)),
                   sprintf("italic(c)[2] * \" = \" * \"%s\"", evs_num(x))),
      col = ifelse(which == "c1", EVS_COL$edge_rem, EVS_COL$edge_lead)
    )

  Nmax <- max(d$N)
  h    <- 0.24
  lab_df <- distinct(d, method, y)

  ggplot(d) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - h, ymax = y + h,
                  fill = regime), colour = "grey45", linewidth = 0.2) +
    geom_text(aes(x = xmid, y = y, label = evs_comma(n)), size = pt2mm(5.6)) +
    geom_text(data = kn, aes(x = x, y = y + h + 0.09, label = lab),
              colour = kn$col, size = pt2mm(5.2), vjust = 0,
              parse = TRUE, inherit.aes = FALSE) +
    geom_text(data = lab_df, aes(x = -0.015 * Nmax, y = y, label = method),
              hjust = 1, size = pt2mm(6.2), inherit.aes = FALSE) +
    scale_fill_manual(values = c("Remainder"    = EVS_COL$band_rem,
                                 "Divergence"   = EVS_COL$band_div,
                                 "Leading edge" = EVS_COL$band_lead)) +
    scale_x_continuous(labels = evs_comma, breaks = scales::pretty_breaks(4)) +
    coord_cartesian(xlim = c(-0.40 * Nmax, Nmax),
                    ylim = c(0.30, max(d$y) + 0.58),
                    expand = FALSE, clip = "off") +
    labs(tag = tag, title = "Regime widths (PAS counts)",
         x = "PAS rank by |PC1 loading| (ascending)", y = NULL) +
    evs_theme() +
    theme(
      axis.line.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.text.y       = element_blank(),
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.background = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 1))
}

# =============================================================================
# 4. FIGURES 2-3 PANELS — Pareto frontier and knee selection
# =============================================================================

evs_norm01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!is.finite(diff(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

# Reproduces pareto_knee() exactly, so the inset shows the criterion that was
# actually optimised rather than a redrawn approximation.
evs_chord_distance <- function(frontier) {
  x <- evs_norm01(frontier$weighted_penalty)
  y <- evs_norm01(frontier$weighted_benefit)
  n <- length(x)
  if (n < 2L) return(rep(0, n))
  x1 <- x[1]; y1 <- y[1]; x2 <- x[n]; y2 <- y[n]
  den <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
  if (!is.finite(den) || den == 0) return(y - x)
  abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / den
}

panel_pareto <- function(method, comparison, scan, frontier, kstar, tag,
                         method_col) {
  fr <- frontier[order(frontier$weighted_penalty), , drop = FALSE]
  fr$chord <- evs_chord_distance(fr)

  j <- which.min(abs(fr$k - kstar))
  if (!length(j)) stop("k* not found on the frontier for ", comparison)
  sel <- fr[j, , drop = FALSE]

  n <- nrow(fr)
  chord_line <- data.frame(
    x = c(fr$weighted_penalty[1], fr$weighted_penalty[n]),
    y = c(fr$weighted_benefit[1],  fr$weighted_benefit[n])
  )
  # Vertical drop from k* to the endpoint chord: a visual cue for the gap the
  # knee criterion maximises.
  drop_y <- if (diff(chord_line$x) == 0) {
    mean(chord_line$y)
  } else {
    stats::approx(chord_line$x, chord_line$y,
                  xout = sel$weighted_penalty, rule = 2)$y
  }

  fr_plot <- fr[evs_thin(n, 1500L), , drop = FALSE]

  xr <- range(fr$weighted_penalty)
  yr <- range(fr$weighted_benefit)
  xd <- if (diff(xr) > 0) diff(xr) else 1
  yd <- if (diff(yr) > 0) diff(yr) else 1
  xlim <- c(xr[1] - 0.03 * xd, xr[2] + 0.05 * xd)
  ylim <- c(yr[1] - 0.05 * yd, yr[2] + 0.10 * yd)

  cmax <- max(fr$chord, na.rm = TRUE)
  if (!is.finite(cmax) || cmax <= 0) cmax <- 1

  ins <- ggplot(fr_plot, aes(k, chord)) +
    geom_line(colour = "grey25", linewidth = 0.25) +
    geom_vline(xintercept = sel$k, colour = EVS_COL$knee,
               linetype = "22", linewidth = 0.3) +
    scale_x_continuous(labels = evs_comma,
                       breaks = range(pretty(fr_plot$k, 3))) +
    scale_y_continuous(breaks = c(0, signif(cmax, 2)),
                       limits = c(0, cmax * 1.2)) +
    labs(x = expression(italic(k)), y = "dist. to chord") +
    evs_theme_inset()

  ggplot() +
    geom_line(data = chord_line, aes(x, y), colour = EVS_COL$chord,
              linetype = "22", linewidth = 0.3) +
    geom_line(data = fr_plot, aes(weighted_penalty, weighted_benefit),
              colour = method_col, linewidth = 0.5) +
    annotate("segment",
             x = sel$weighted_penalty, xend = sel$weighted_penalty,
             y = sel$weighted_benefit, yend = drop_y,
             colour = "grey45", linewidth = 0.25) +
    annotation_custom(
      ggplotGrob(ins),
      xmin = xlim[1] + 0.50 * diff(xlim), xmax = xlim[1] + 1.00 * diff(xlim),
      ymin = ylim[1] + 0.04 * diff(ylim), ymax = ylim[1] + 0.46 * diff(ylim)
    ) +
    geom_point(data = sel, aes(weighted_penalty, weighted_benefit),
               shape = 23, size = 1.5, stroke = 0.3,
               fill = EVS_COL$knee, colour = "white") +
    annotate("text",
             x = sel$weighted_penalty - 0.015 * diff(xlim),
             y = sel$weighted_benefit + 0.035 * diff(ylim),
             label = sprintf("italic(k)^\"*\" * \" = \" * \"%s\"", evs_num(sel$k)),
             parse = TRUE, hjust = 1, vjust = 0,
             size = pt2mm(6.2), colour = EVS_COL$knee, fontface = "bold") +
    scale_x_continuous(labels = evs_comma) +
    scale_y_continuous(labels = evs_comma) +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(
      tag   = tag,
      title = evs_comparison_label(comparison),
      x     = expression(paste("Weighted Remainder disagreement, ", italic(R)[w])),
      y     = expression(paste("Weighted retained-site benefit, ", italic(G)[w]))
    ) +
    evs_theme()
}

# =============================================================================
# 5. FIGURE 4 PANELS — cross-method summary
# =============================================================================

evs_prep_summary <- function(summary_df, comparisons) {
  summary_df %>%
    mutate(
      method     = factor(evs_method_label(method), levels = names(EVS_METHOD_COL)),
      comparison = factor(evs_comparison_label(comparison),
                          levels = evs_comparison_label(names(comparisons)))
    )
}

panel_kstar <- function(summary_df, comparisons, tag) {
  d <- evs_prep_summary(summary_df, comparisons)
  ggplot(d, aes(comparison, selected_k, fill = method)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62) +
    geom_text(aes(label = evs_comma(selected_k)),
              position = position_dodge(width = 0.72),
              vjust = -0.45, size = pt2mm(5.4)) +
    scale_fill_manual(values = EVS_METHOD_COL) +
    scale_y_continuous(labels = evs_comma,
                       expand = expansion(mult = c(0, 0.16))) +
    labs(tag = tag, title = "Dataset-specific cutoffs", x = NULL,
         y = expression(paste("Selected cutoff, ", italic(k)^"*", " (PAS)"))) +
    evs_theme() +
    evs_legend_inside(0.99, 0.99, c(1, 1), direction = "horizontal")
}

panel_hc_membership <- function(overlap_df, comparisons, tag) {
  d <- overlap_df %>%
    mutate(
      cpm_only   = CPM_high_confidence    - overlap_high_confidence,
      shared     = overlap_high_confidence,
      deseq_only = DESeq2_high_confidence - overlap_high_confidence,
      comparison = factor(evs_comparison_label(comparison),
                          levels = evs_comparison_label(names(comparisons)))
    )

  tot <- d %>%
    transmute(comparison,
              total   = cpm_only + shared + deseq_only,
              jaccard = jaccard_high_confidence)

  long <- d %>%
    select(comparison, cpm_only, shared, deseq_only) %>%
    tidyr::pivot_longer(-comparison, names_to = "set", values_to = "n") %>%
    mutate(set = factor(recode(set,
                               cpm_only   = "CPM-EVS only",
                               shared     = "Shared",
                               deseq_only = "DESeq2-EVS only"),
                        levels = c("CPM-EVS only", "Shared", "DESeq2-EVS only")))

  ggplot(long, aes(comparison, n, fill = set)) +
    geom_col(width = 0.56, colour = "white", linewidth = 0.2) +
    geom_text(aes(label = evs_comma(n)),
              position = position_stack(vjust = 0.5),
              size = pt2mm(5.0), colour = "white") +
    geom_text(data = tot, inherit.aes = FALSE,
              aes(x = comparison, y = total,
                  label = sprintf("italic(J) * \" = \" * \"%.2f\"", jaccard)),
              parse = TRUE, vjust = -0.6, size = pt2mm(5.4), fontface = "bold") +
    scale_fill_manual(values = c("CPM-EVS only"    = EVS_COL$cpm,
                                 "Shared"          = EVS_COL$shared,
                                 "DESeq2-EVS only" = EVS_COL$deseq)) +
    scale_y_continuous(labels = evs_comma,
                       expand = expansion(mult = c(0, 0.20))) +
    labs(tag = tag, title = "High-confidence EVS membership", x = NULL,
         y = "High-confidence PAS (count)") +
    evs_theme() +
    evs_legend_inside(0.5, 1.0, c(0.5, 1), direction = "horizontal")
}

panel_composition <- function(summary_df, comparisons, tag) {
  d <- evs_prep_summary(summary_df, comparisons) %>%
    select(method, comparison, joint_n, opposite_le_n,
           opposite_divergence_n, remainder_cross_n) %>%
    tidyr::pivot_longer(-c("method", "comparison"),
                        names_to = "cls", values_to = "n") %>%
    mutate(cls = factor(recode(cls,
                               joint_n               = "Joint",
                               opposite_le_n         = "Opposite LE",
                               opposite_divergence_n = "Opposite Divergence",
                               remainder_cross_n     = "Opposite Remainder"),
                        levels = c("Joint", "Opposite LE",
                                   "Opposite Divergence", "Opposite Remainder"))) %>%
    group_by(method, comparison) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  ggplot(d, aes(method, pct, fill = cls)) +
    geom_col(width = 0.62, colour = "white", linewidth = 0.2) +
    geom_text(aes(label = ifelse(pct >= 8, sprintf("%.0f", pct), "")),
              position = position_stack(vjust = 0.5), size = pt2mm(4.8)) +
    facet_wrap(~ comparison, nrow = 1, strip.position = "bottom") +
    scale_fill_manual(values = c("Joint"               = EVS_COL$joint,
                                 "Opposite LE"         = EVS_COL$opp_le,
                                 "Opposite Divergence" = EVS_COL$opp_div,
                                 "Opposite Remainder"  = EVS_COL$opp_rem)) +
    scale_x_discrete(labels = c("CPM-EVS" = "CPM", "DESeq2-EVS" = "DESeq2")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(tag = tag,
         title = expression(paste("Site composition at ", italic(k)^"*")),
         x = NULL,
         y = expression(paste("Composition of top-", italic(k)^"*", " union (%)"))) +
    evs_theme() +
    theme(
      strip.background  = element_blank(),
      strip.placement   = "outside",
      strip.text        = element_text(size = EVS_BASE_PT - 0.5,
                                       margin = margin(t = 1.5)),
      panel.spacing.x   = unit(3, "pt"),
      axis.text.x       = element_text(size = EVS_BASE_PT - 1.6, angle = 90,
                                       hjust = 1, vjust = 0.5),
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.background = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}

panel_le_fraction <- function(summary_df, comparisons, tag) {
  d <- evs_prep_summary(summary_df, comparisons) %>%
    mutate(pct = 100 * k_over_leading_edge)
  ggplot(d, aes(comparison, pct, fill = method)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62) +
    geom_text(aes(label = sprintf("%.1f", pct)),
              position = position_dodge(width = 0.72),
              vjust = -0.45, size = pt2mm(5.4)) +
    scale_fill_manual(values = EVS_METHOD_COL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(tag = tag, title = "Fraction of candidate leading edge selected",
         x = NULL,
         y = expression(paste(italic(k)^"*", " as % of leading-edge domain"))) +
    evs_theme() +
    evs_legend_inside(0.99, 0.99, c(1, 1), direction = "horizontal")
}

# =============================================================================
# 6. AUTO-GENERATED FIGURE LEGENDS
# =============================================================================
# Legends are written with the run's real numbers, so the manuscript text
# cannot drift out of sync with the figures.

evs_write_legends <- function(method_objects, summary_df, overlap_df,
                              comparisons, out_root) {

  kn <- lapply(method_objects, function(o) o$knot)
  N  <- kn[[1]]$N

  knot_txt <- paste(vapply(names(kn), function(m) {
    k <- kn[[m]]
    sprintf("%s, c1 = %s and c2 = %s", evs_method_label(m),
            evs_num(k$c1), evs_num(k$c2))
  }, character(1)), collapse = "; ")

  # What fraction of candidate k are actually non-dominated? If the frontier
  # contains nearly all of them, the knee (not dominance filtering) is doing
  # the work, and the legend should say so rather than imply a sharp trade-off.
  fr_frac <- vapply(names(method_objects), function(m) {
    ob <- method_objects[[m]]
    mean(vapply(names(comparisons), function(nm) {
      nrow(ob$frontiers[[nm]]) / nrow(ob$scans[[nm]])
    }, numeric(1)))
  }, numeric(1))

  kstar_of <- function(m) {
    vapply(names(comparisons),
           function(nm) method_objects[[m]]$pairs[[nm]]$kstar, numeric(1))
  }
  kstar_txt <- paste(vapply(names(method_objects), function(m) {
    sprintf("%s, k* = %s", evs_method_label(m),
            paste(evs_num(kstar_of(m)), collapse = ", "))
  }, character(1)), collapse = "; ")

  jac <- sprintf("%.2f-%.2f",
                 min(overlap_df$jaccard_high_confidence, na.rm = TRUE),
                 max(overlap_df$jaccard_high_confidence, na.rm = TRUE))

  L <- c(
    "# Figure legends",
    "",
    "<!-- Auto-generated by EVS_Manuscript_Figures.R. Do not hand-edit: edit the",
    "     generator instead, so the numbers can never drift from the figures. -->",
    "",
    "## Abbreviations",
    "",
    "CPM, counts per million; c1 and c2, shared lower and upper rank knots;",
    "D(r), cumulative PC1-excess-variance divergence at rank r; DESeq2,",
    "median-of-ratios normalisation; EVS, excess-variance selection; F_E(r),",
    "cumulative excess-variance mass; F_P(r), cumulative PC1 variance-contribution",
    "mass; G_w, weighted retained-site benefit; HC, high confidence; J, Jaccard",
    "index; k*, selected per-comparison cutoff; LE, leading edge; N, size of the",
    "experiment-wide PAS universe; PAS, polyadenylation site; PC1, first principal",
    "component; R_w, distance-weighted Remainder disagreement; RT, [AUTHORS: define];",
    "ZT, zeitgeber time.",
    "",
    "> Two abbreviations could not be recovered from the analysis code and must be",
    "> set by the authors before submission: **EVS** (written above as",
    "> \"excess-variance selection\") and **RT**. Edit `evs_write_legends()` once",
    "> the intended expansions are fixed, and they will propagate to every run.",
    "",
    "---",
    "",
    "**Figure 1. Definition of the Remainder, Divergence and Leading-Edge rank",
    "regimes under two normalisation strategies.**",
    sprintf(paste0(
      "(A, B) Cumulative PC1-excess-variance divergence, D(r) = F_E(r) - F_P(r), ",
      "against PAS rank by ascending absolute PC1 loading, for (A) CPM-EVS and ",
      "(B) DESeq2-EVS. Grey lines show all eight experimental arms, the coloured ",
      "line the across-arm median, and the dashed black line the shared two-knot ",
      "piecewise-linear fit estimated jointly across arms by minimising the summed ",
      "squared residual. Vertical dashed lines mark the shared knots (%s), which ",
      "partition the universe of N = %s PAS into the Remainder (rank < c1), ",
      "Divergence (c1 <= rank <= c2) and Leading Edge (rank > c2), shaded left to ",
      "right. (C) Residuals of the shared-knot fit, D_g(r) - Dhat_g(r), binned by ",
      "rank; lines show the across-arm median and shaded envelopes the ",
      "interquartile range. (D) Widths of the three regimes for each method, ",
      "annotated with PAS counts. c1 and c2 define rank regimes common to all arms ",
      "and are not themselves the selected cutoff k*."),
      knot_txt, evs_num(N)),
    "",
    "---",
    "",
    "**Figure 2. Per-comparison cutoff selection under CPM-EVS.**",
    sprintf(paste0(
      "(A-D) Pareto frontier of weighted retained-site benefit (G_w) against ",
      "distance-weighted Remainder disagreement (R_w) over all candidate cutoffs ",
      "k, shown separately for each RT/ZT comparison; each comparison is optimised ",
      "independently over 1 <= k <= N - c2. The grey dashed line joins the frontier ",
      "endpoints, and the filled diamond marks k*, the point of maximum ",
      "perpendicular distance from that chord once both axes are rescaled to ",
      "[0, 1]; the vertical segment shows the corresponding gap. Insets plot that ",
      "distance against k, so the selected maximum is directly visible. ",
      "Non-dominated candidates make up %.0f%% of the k values scanned under ",
      "CPM-EVS, so the frontier is close to monotone and the knee criterion, ",
      "rather than dominance filtering, determines the cutoff. Selected values: ",
      "k* = %s for %s respectively; full diagnostics are in Table 1."),
      100 * fr_frac[["CPM_EVS"]],
      paste(evs_num(kstar_of("CPM_EVS")), collapse = ", "),
      paste(evs_comparison_label(names(comparisons)), collapse = ", ")),
    "",
    "---",
    "",
    "**Figure 3. Per-comparison cutoff selection under DESeq2-EVS.**",
    sprintf(paste0(
      "Panels, axes and annotations are as in Figure 2, computed from ",
      "log1p-transformed DESeq2 median-of-ratios normalised counts. Non-dominated ",
      "candidates make up %.0f%% of the k values scanned. Selected values: ",
      "k* = %s."),
      100 * fr_frac[["DESeq2_EVS"]],
      paste(evs_num(kstar_of("DESeq2_EVS")), collapse = ", ")),
    "",
    "---",
    "",
    "**Figure 4. CPM-EVS and DESeq2-EVS compared across all four RT/ZT",
    "comparisons.**",
    sprintf(paste0(
      "(A) Selected cutoff k* for each comparison and method (%s). (B) ",
      "High-confidence PAS membership, partitioned into sites recovered only by ",
      "CPM-EVS, only by DESeq2-EVS, or by both; J is the Jaccard index of the two ",
      "high-confidence sets (range across comparisons, %s). High-confidence sites ",
      "are Joint sites together with disjoint sites whose opposite-arm rank falls ",
      "in the Leading Edge or Divergence regime. (C) Composition of the selected ",
      "top-k* union at k*, as a percentage of union size; segment labels are ",
      "omitted below 8%%. Opposite-Remainder sites are retained in the exported ",
      "union with a penalty flag rather than discarded. (D) k* as a percentage of ",
      "the candidate Leading-Edge domain, N - c2. Bars in A and D are annotated ",
      "with their values; exact counts are given in Tables 1 and 2."),
      kstar_txt, jac),
    "",
    "---",
    "",
    "*Figures were authored at 180 mm final width and are supplied as vector PDF",
    "and 600 dpi raster. Colours follow the Okabe-Ito palette and remain",
    "distinguishable under common forms of colour-vision deficiency and in",
    "greyscale.*"
  )

  writeLines(L, file.path(out_root, "Figure_Legends.md"))
  message("  wrote Figure_Legends.md")
  invisible(TRUE)
}

# =============================================================================
# 7. TOP-LEVEL RENDER
# =============================================================================

evs_render_all <- function(method_objects, summary_df, overlap_df,
                           comparisons, fig_dir, out_root) {

  dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  tags <- LETTERS[1:4]

  message("Figure 1 ...")
  f1 <- list(
    panel_regime("CPM_EVS", method_objects[["CPM_EVS"]]$arms,
                 method_objects[["CPM_EVS"]]$knot, "A", show_legend = TRUE),
    panel_regime("DESeq2_EVS", method_objects[["DESeq2_EVS"]]$arms,
                 method_objects[["DESeq2_EVS"]]$knot, "B"),
    panel_residuals(method_objects, "C"),
    panel_regime_widths(method_objects, "D")
  )
  evs_save(f1, file.path(fig_dir, "Figure_1_Regime_Definition"),
           nrow = 2, ncol = 2, height_mm = 150)

  for (m in names(method_objects)) {
    fig_no <- if (m == "CPM_EVS") 2L else 3L
    message("Figure ", fig_no, " ...")
    ob  <- method_objects[[m]]
    col <- unname(EVS_METHOD_COL[evs_method_label(m)])
    pl  <- lapply(seq_along(comparisons), function(i) {
      nm <- names(comparisons)[i]
      panel_pareto(m, nm, ob$scans[[nm]], ob$frontiers[[nm]],
                   ob$pairs[[nm]]$kstar, tags[i], col)
    })
    evs_save(pl, file.path(fig_dir,
                           sprintf("Figure_%d_%s_Pareto_Cutoffs", fig_no, m)),
             nrow = 2, ncol = 2, height_mm = 145)
  }

  message("Figure 4 ...")
  f4 <- list(
    panel_kstar(summary_df, comparisons, "A"),
    panel_hc_membership(overlap_df, comparisons, "B"),
    panel_composition(summary_df, comparisons, "C"),
    panel_le_fraction(summary_df, comparisons, "D")
  )
  evs_save(f4, file.path(fig_dir, "Figure_4_CPM_vs_DESeq2_EVS_Summary"),
           nrow = 2, ncol = 2, height_mm = 150)

  evs_write_legends(method_objects, summary_df, overlap_df, comparisons, out_root)
  message("Figures complete: ", fig_dir)
  invisible(TRUE)
}

# =============================================================================
# 8. SELF-TEST (synthetic data; verifies rendering without the count matrix)
# =============================================================================

evs_selftest <- function(dir = file.path(getwd(), "evs_figure_selftest")) {
  set.seed(11)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  N <- 32104L
  comps <- list(RT0_ZT6 = 1, RT2_ZT8 = 2, RT4_ZT10 = 3, RT8_ZT14 = 4)
  knots <- list(CPM_EVS = c(10451L, 23509L), DESeq2_EVS = c(9120L, 21860L))

  mk_method <- function(m) {
    c1 <- knots[[m]][1]; c2 <- knots[[m]][2]
    x  <- (seq_len(N) - 1) / (N - 1)
    x1 <- (c1 - 1) / (N - 1); x2 <- (c2 - 1) / (N - 1)
    base <- -0.05 + 1.20 * x - 0.55 * pmax(x - x1, 0) - 2.25 * pmax(x - x2, 0)
    base <- base - base[N]
    arms <- lapply(1:8, function(g) {
      d <- base + rnorm(N, 0, 0.004) + 0.02 * sin(x * 9 + g)
      data.frame(rank = seq_len(N), D = d - d[N])
    })
    names(arms) <- paste0("arm", 1:8)
    fitted <- matrix(rep(base, 8), ncol = 8)
    knot <- list(c1 = c1, c2 = c2, N = N, fitted = fitted, SSE = NA_real_)

    K <- N - c2
    scans <- list(); fronts <- list(); pairs <- list()
    for (nm in names(comps)) {
      k   <- seq_len(K)
      pen <- cummax(1300 * (k / K)^1.35)
      ben <- cummax(5900 * (k / K)^0.92)
      s <- data.frame(k = k, weighted_penalty = pen, weighted_benefit = ben,
                      joint_n = round(k * 0.20), opposite_le_n = round(k * 0.22),
                      opposite_divergence_n = round(k * 0.80),
                      remainder_cross_n = round(k * 0.58),
                      union_n = 2 * k - round(k * 0.20))
      xx <- (pen - min(pen)) / diff(range(pen))
      yy <- (ben - min(ben)) / diff(range(ben))
      dd <- abs((yy[K] - yy[1]) * xx - (xx[K] - xx[1]) * yy +
                  xx[K] * yy[1] - yy[K] * xx[1]) /
        sqrt((yy[K] - yy[1])^2 + (xx[K] - xx[1])^2)
      scans[[nm]]  <- s
      fronts[[nm]] <- s
      pairs[[nm]]  <- list(kstar = k[which.max(dd)])
    }
    list(arms = arms, knot = knot, scans = scans, frontiers = fronts, pairs = pairs)
  }

  mo <- list(CPM_EVS = mk_method("CPM_EVS"), DESeq2_EVS = mk_method("DESeq2_EVS"))

  summary_df <- bind_rows(lapply(names(mo), function(m) {
    bind_rows(lapply(names(comps), function(nm) {
      ks <- mo[[m]]$pairs[[nm]]$kstar
      s  <- mo[[m]]$scans[[nm]][ks, ]
      data.frame(method = m, comparison = nm, N = N,
                 c1 = mo[[m]]$knot$c1, c2 = mo[[m]]$knot$c2,
                 selected_k = ks,
                 k_over_leading_edge = ks / (N - mo[[m]]$knot$c2),
                 joint_n = s$joint_n, opposite_le_n = s$opposite_le_n,
                 opposite_divergence_n = s$opposite_divergence_n,
                 remainder_cross_n = s$remainder_cross_n)
    }))
  }))

  overlap_df <- bind_rows(lapply(names(comps), function(nm) {
    a  <- summary_df$joint_n[summary_df$method == "CPM_EVS" &
                               summary_df$comparison == nm] * 4
    b  <- summary_df$joint_n[summary_df$method == "DESeq2_EVS" &
                               summary_df$comparison == nm] * 4
    ov <- round(min(a, b) * 0.68)
    data.frame(comparison = nm, CPM_high_confidence = a,
               DESeq2_high_confidence = b, overlap_high_confidence = ov,
               union_high_confidence = a + b - ov,
               jaccard_high_confidence = ov / (a + b - ov))
  }))

  evs_render_all(mo, summary_df, overlap_df, comps, dir, dir)
  message("Self-test output: ", dir)
  invisible(dir)
}

if (!interactive() &&
    any(grepl("--selftest", commandArgs(trailingOnly = TRUE), fixed = TRUE))) {
  evs_selftest()
}
