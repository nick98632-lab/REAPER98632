#!/usr/bin/env Rscript

# =============================================================================
# ANCHOR SWEEP: LEFT/RIGHT NB2 DISPERSION ACROSS THE PC1-RANK AXIS
# =============================================================================
#
# This script is standalone. It does not modify SEQUENCE.R, EMPERICALCUTOFF.R,
# REGIME_DIAGNOSTIC.r, REMNB1_LEADNB2.R, or LFC_BALANCE_DIAGNOSTIC.R.
#
# For each arm of each comparison, features are ranked by absolute PC1
# loading. At a candidate anchor rank:
#   RIGHT = all ranks from the anchor through the most extreme end of the
#           ranking
#   LEFT  = an equal-sized matched block immediately below the anchor
# A likelihood-ratio test compares a shared NB2 dispersion parameter (alpha)
# for LEFT and RIGHT against separate alphas for each (see REMNB1_LEADNB2.R
# for the full LRT specification; this script uses the same raw-count
# likelihood). Because RIGHT always extends to the most extreme rank, LEFT
# (equal-sized) is only feasible once the anchor is past roughly the
# midpoint of the ranking; the swept anchor range is restricted accordingly.
#
# This is run at a systematic range of anchor positions per arm, not only
# the established leading-edge boundary (paper-reference k = 5000), to
# characterize how the LEFT/RIGHT dispersion gap (diff_alpha) varies across
# the whole rank axis, and where it is maximized relative to the
# established boundary.
#
# OUTPUTS:
#   Table_Anchor_Sweep.csv            one row per anchor position tested
#   Table_Peak_vs_Established.csv     established boundary vs. the anchor
#                                      that maximizes diff_alpha, per arm
#   Figure_Anchor_Sweep_<comparison>.png
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Settings (kept consistent with REMNB1_LEADNB2.R / REGIME_DIAGNOSTIC.r)
# -----------------------------------------------------------------------------

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT   <- "/root/REAPER98632/exports/anchor_sweep_diagnostic"
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

GROUP_PATTERNS <- c(
  RT0  = "^R0_",
  ZT6  = "^ZT6_",
  RT2  = "^R2_",
  ZT8  = "^ZT8_",
  RT4  = "^R4_",
  ZT10 = "^ZT10_",
  RT8  = "^R8_",
  ZT14 = "^ZT14_"
)

COMPARISONS <- list(
  RT0_ZT6  = c(control = "RT0", treatment = "ZT6"),
  RT2_ZT8  = c(control = "RT2", treatment = "ZT8"),
  RT4_ZT10 = c(control = "RT4", treatment = "ZT10"),
  RT8_ZT14 = c(control = "RT8", treatment = "ZT14")
)

# Established boundary reference: paper-reference k = 5000. total_n varies
# slightly per arm (a handful of all-zero features get dropped per arm), so
# the reference anchor is computed per arm as total_n - 5000 + 1, matching
# REMNB1_LEADNB2.R's reference_rank convention exactly.
PAPER_REFERENCE_K <- 5000L

# Sweep step size across the feasible anchor range.
SWEEP_STEP <- 1000L

# -----------------------------------------------------------------------------
# Data loading and ranking (mirrors REMNB1_LEADNB2.R conventions exactly)
# -----------------------------------------------------------------------------

read_count_matrix <- function(path, group_patterns) {
  if (!file.exists(path)) stop("Count file does not exist: ", path)
  raw_df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(raw_df) < 1L || ncol(raw_df) < 2L) stop("Count file is empty or malformed.")

  sample_idx <- sort(unique(unlist(
    lapply(group_patterns, function(pattern) grep(pattern, colnames(raw_df)))
  )))
  if (length(sample_idx) == 0L) stop("No sample columns matched GROUP_PATTERNS.")

  feature_ids <- make.unique(as.character(raw_df[[1]]))
  count_mat <- as.matrix(raw_df[, sample_idx, drop = FALSE])
  storage.mode(count_mat) <- "numeric"
  rownames(count_mat) <- feature_ids

  count_mat[!is.finite(count_mat)] <- 0
  count_mat <- pmax(count_mat, 0)

  keep <- rowSums(count_mat) > 0
  count_mat[keep, , drop = FALSE]
}

normalize_cpm_log1p <- function(count_mat_arm) {
  lib_sizes <- colSums(count_mat_arm, na.rm = TRUE)
  lib_sizes[!is.finite(lib_sizes) | lib_sizes <= 0] <- 1
  cpm <- sweep(count_mat_arm, 2, lib_sizes / 1e6, "/")
  log1p(cpm)
}

compute_abs_pc1_loadings <- function(norm_mat_arm) {
  pca <- stats::prcomp(t(norm_mat_arm), center = TRUE, scale. = FALSE, rank. = 1)
  out <- abs(pca$rotation[, 1L])
  out[!is.finite(out)] <- 0
  out
}

# -----------------------------------------------------------------------------
# LRT machinery (identical to REMNB1_LEADNB2.R after its optimizer fix)
# -----------------------------------------------------------------------------

nb2_region_loglik <- function(alpha, counts_block, mu_vec) {
  if (!is.finite(alpha) || alpha <= 0) return(-Inf)
  size <- 1 / alpha
  mu_rep <- rep(mu_vec, times = ncol(counts_block))
  counts_vec <- as.vector(counts_block)
  keep <- is.finite(counts_vec) & is.finite(mu_rep) & mu_rep > 0 & counts_vec >= 0
  if (!any(keep)) return(-Inf)
  sum(stats::dnbinom(counts_vec[keep], size = size, mu = mu_rep[keep], log = TRUE))
}

fit_region_alpha_mle <- function(counts_block, mu_vec, log_alpha_lower = -15, log_alpha_upper = 15) {
  obj <- function(log_a) -nb2_region_loglik(exp(log_a), counts_block, mu_vec)
  opt <- stats::optimize(obj, lower = log_alpha_lower, upper = log_alpha_upper, tol = 1e-8)
  at_bound <- (opt$minimum <= log_alpha_lower + 1e-6) || (opt$minimum >= log_alpha_upper - 1e-6)
  list(alpha_mle = exp(opt$minimum), loglik = -opt$objective, at_bound = at_bound)
}

lrt_at_anchor <- function(count_mat_arm, mu_by_rank, rank_order, anchor_rank, total_n) {
  right_idx <- seq.int(anchor_rank, total_n)
  right_n <- length(right_idx)
  left_end <- anchor_rank - 1L
  left_start <- left_end - right_n + 1L
  if (left_start < 1L) return(NULL)
  left_idx <- seq.int(left_start, left_end)

  left_ids <- rank_order[left_idx]
  right_ids <- rank_order[right_idx]

  left_block <- count_mat_arm[left_ids, , drop = FALSE]
  right_block <- count_mat_arm[right_ids, , drop = FALSE]
  left_mu <- mu_by_rank[left_idx]
  right_mu <- mu_by_rank[right_idx]

  fit_left <- fit_region_alpha_mle(left_block, left_mu)
  fit_right <- fit_region_alpha_mle(right_block, right_mu)

  pooled_block <- rbind(left_block, right_block)
  pooled_mu <- c(left_mu, right_mu)
  fit_pooled <- fit_region_alpha_mle(pooled_block, pooled_mu)

  lrt_stat <- max(2 * ((fit_left$loglik + fit_right$loglik) - fit_pooled$loglik), 0)
  lrt_p <- stats::pchisq(lrt_stat, df = 1, lower.tail = FALSE)
  any_at_bound <- fit_left$at_bound || fit_right$at_bound || fit_pooled$at_bound

  data.frame(
    anchor_rank = anchor_rank,
    left_n = length(left_idx),
    right_n = length(right_idx),
    alpha_left = fit_left$alpha_mle,
    alpha_right = fit_right$alpha_mle,
    diff_alpha = fit_right$alpha_mle - fit_left$alpha_mle,
    lrt_stat = lrt_stat,
    lrt_p = lrt_p,
    any_at_bound = any_at_bound,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------

message("Reading count matrix...")
count_mat <- read_count_matrix(COUNT_FILE, GROUP_PATTERNS)

sweep_rows <- list()
peak_rows <- list()

for (comparison_name in names(COMPARISONS)) {
  pair <- COMPARISONS[[comparison_name]]

  for (arm_name in c("control", "treatment")) {
    arm_label <- pair[[arm_name]]
    arm_ids <- grep(GROUP_PATTERNS[[arm_label]], colnames(count_mat), value = TRUE)
    if (length(arm_ids) < 2L) next

    count_mat_arm <- count_mat[, arm_ids, drop = FALSE]
    norm_mat_arm <- normalize_cpm_log1p(count_mat_arm)

    abs_loadings <- compute_abs_pc1_loadings(norm_mat_arm)
    rank_order_idx <- order(abs_loadings, decreasing = FALSE)
    rank_order <- rownames(count_mat_arm)[rank_order_idx]
    total_n <- length(rank_order)

    mu_by_rank <- rowMeans(count_mat_arm[rank_order, , drop = FALSE])

    established_anchor <- total_n - PAPER_REFERENCE_K + 1L

    # Feasible range: anchor must be high enough that an equal-sized LEFT
    # block still fits below rank 1.
    candidate_anchors <- seq(
      from = max(SWEEP_STEP, round(total_n * 0.5)),
      to = total_n - 100L,
      by = SWEEP_STEP
    )
    candidate_anchors <- sort(unique(c(candidate_anchors, established_anchor)))

    message(sprintf(
      "[%s %s] total_n=%d, established_anchor=%d, sweeping %d candidate anchors",
      comparison_name, arm_name, total_n, established_anchor, length(candidate_anchors)
    ))

    arm_rows <- list()
    for (anchor in candidate_anchors) {
      res <- lrt_at_anchor(count_mat_arm, mu_by_rank, rank_order, anchor, total_n)
      if (is.null(res)) next
      res$comparison <- comparison_name
      res$arm <- arm_name
      res$is_established_boundary <- (anchor == established_anchor)
      arm_rows[[length(arm_rows) + 1L]] <- res
      sweep_rows[[length(sweep_rows) + 1L]] <- res
    }

    arm_table <- dplyr::bind_rows(arm_rows)
    if (nrow(arm_table) > 0L) {
      peak_row <- arm_table[which.max(arm_table$diff_alpha), ]
      established_row <- arm_table[arm_table$is_established_boundary, ][1, ]

      peak_rows[[length(peak_rows) + 1L]] <- data.frame(
        comparison = comparison_name,
        arm = arm_name,
        established_anchor = established_row$anchor_rank,
        established_diff_alpha = established_row$diff_alpha,
        peak_anchor = peak_row$anchor_rank,
        peak_diff_alpha = peak_row$diff_alpha,
        rank_distance = established_row$anchor_rank - peak_row$anchor_rank,
        pct_of_peak_achieved = 100 * established_row$diff_alpha / peak_row$diff_alpha,
        stringsAsFactors = FALSE
      )
    }
  }
}

sweep_table <- dplyr::bind_rows(sweep_rows)
write.csv(sweep_table, file.path(OUT_ROOT, "Table_Anchor_Sweep.csv"), row.names = FALSE)

peak_table <- dplyr::bind_rows(peak_rows)
write.csv(peak_table, file.path(OUT_ROOT, "Table_Peak_vs_Established.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Figures: one per comparison, faceted by arm
# -----------------------------------------------------------------------------

MARKER_LEVELS <- c("Sweep point", "At search bound (unreliable)", "Established boundary", "Empirical peak")
MARKER_COLORS <- c(
  "Sweep point" = "#3B6FA0",
  "At search bound (unreliable)" = "#B0B0B0",
  "Established boundary" = "#C0392B",
  "Empirical peak" = "#E8A33D"
)
MARKER_SHAPES <- c(
  "Sweep point" = 16,
  "At search bound (unreliable)" = 4,
  "Established boundary" = 17,
  "Empirical peak" = 18
)
MARKER_SIZES <- c(
  "Sweep point" = 1.8,
  "At search bound (unreliable)" = 2.4,
  "Established boundary" = 3.6,
  "Empirical peak" = 3.6
)

for (comparison_name in names(COMPARISONS)) {
  sub <- sweep_table[sweep_table$comparison == comparison_name, , drop = FALSE]
  if (nrow(sub) == 0L) next

  marker_df <- sub %>%
    group_by(arm) %>%
    mutate(
      marker_type = case_when(
        is_established_boundary ~ "Established boundary",
        anchor_rank == anchor_rank[which.max(diff_alpha)] ~ "Empirical peak",
        any_at_bound ~ "At search bound (unreliable)",
        TRUE ~ "Sweep point"
      )
    ) %>%
    ungroup() %>%
    mutate(marker_type = factor(marker_type, levels = MARKER_LEVELS))

  subtitle_text <- paste(
    strwrap(
      "Points colored by type below. A curve that stays flat across the sweep would indicate no boundary-specific effect.",
      width = 70
    ),
    collapse = "\n"
  )

  p <- ggplot(marker_df, aes(x = anchor_rank, y = diff_alpha)) +
    geom_line(color = "grey60", linewidth = 0.6) +
    geom_point(aes(color = marker_type, shape = marker_type, size = marker_type)) +
    facet_wrap(~arm, ncol = 1, scales = "free_y") +
    scale_color_manual(values = MARKER_COLORS, breaks = MARKER_LEVELS, drop = FALSE, name = NULL) +
    scale_shape_manual(values = MARKER_SHAPES, breaks = MARKER_LEVELS, drop = FALSE, name = NULL) +
    scale_size_manual(values = MARKER_SIZES, breaks = MARKER_LEVELS, drop = FALSE, guide = "none") +
    guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3.2))) +
    labs(
      title = paste0(comparison_name, ": LEFT/RIGHT dispersion gap across the anchor sweep"),
      subtitle = subtitle_text,
      x = "Anchor rank",
      y = "diff_alpha (alpha_right - alpha_left)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12),
      plot.subtitle = element_text(size = 9)
    )

  ggsave(
    file.path(OUT_ROOT, paste0("Figure_Anchor_Sweep_", comparison_name, ".png")),
    p, width = 10, height = 7.5, dpi = 300
  )
}

message("Diagnostic complete. Outputs written to: ", OUT_ROOT)
message("  Table_Anchor_Sweep.csv")
message("  Table_Peak_vs_Established.csv")
message("  Figure_Anchor_Sweep_<comparison>.png")

# -----------------------------------------------------------------------------
# Zip archives: everything in one download for figures and for tables.
# Wrapped in tryCatch so that if the zip utility is unavailable, the actual
# diagnostic results above are still kept -- only the packaging step is lost.
# -----------------------------------------------------------------------------

message("Creating zip archives...")

zip_result <- tryCatch(
  {
    figure_files <- list.files(OUT_ROOT, pattern = "\\.png$", full.names = TRUE)
    table_files <- list.files(OUT_ROOT, pattern = "\\.csv$", full.names = TRUE)

    if (length(figure_files) > 0L) {
      figures_zip_path <- file.path(OUT_ROOT, "AnchorSweep_Figures.zip")
      if (file.exists(figures_zip_path)) file.remove(figures_zip_path)
      utils::zip(figures_zip_path, files = figure_files, flags = "-j")
      message("  AnchorSweep_Figures.zip (", length(figure_files), " files)")
    } else {
      message("  No .png files found; skipping AnchorSweep_Figures.zip")
    }

    if (length(table_files) > 0L) {
      tables_zip_path <- file.path(OUT_ROOT, "AnchorSweep_Tables.zip")
      if (file.exists(tables_zip_path)) file.remove(tables_zip_path)
      utils::zip(tables_zip_path, files = table_files, flags = "-j")
      message("  AnchorSweep_Tables.zip (", length(table_files), " files)")
    } else {
      message("  No .csv files found; skipping AnchorSweep_Tables.zip")
    }

    TRUE
  },
  error = function(e) {
    message(
      "Zip archive creation failed (individual files above are still ",
      "intact and usable): ", conditionMessage(e)
    )
    FALSE
  }
)
