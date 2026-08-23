#!/usr/bin/env Rscript

# =============================================================================
# ANCHOR SWEEP DIAGNOSTIC: is the LEFT/RIGHT LRT signal specific to the
# established leading-edge boundary, or just a smooth global gradient?
# =============================================================================
#
# This script is standalone. It does not modify SEQUENCE.R, EMPERICALCUTOFF.R,
# REGIME_DIAGNOSTIC.r, REMNB1_LEADNB2.R, or LFC_BALANCE_DIAGNOSTIC.R.
#
# MOTIVATING QUESTION
# REMNB1_LEADNB2.R's likelihood-ratio test compares a LEFT block against a
# RIGHT block that are defined relative to one specific anchor: the
# established leading-edge boundary (near the paper-reference/empirical k*
# cutoff). It found RIGHT more NB2-like than LEFT there, significantly, in
# every arm tested. The open question: is that specific to the established
# boundary, or would picking almost any anchor along the PC1-rank axis show
# the same pattern, because NB2-ness increases smoothly and monotonically
# across the whole ranking with no distinct feature at that particular point?
#
# WHAT THIS SCRIPT DOES
# For each arm (control, treatment) of each comparison, this sweeps the
# anchor rank across a systematic range of positions and, at each one,
# reruns the same LEFT/RIGHT alpha-MLE likelihood-ratio test used in
# REMNB1_LEADNB2.R (raw counts only -- the track with a trustworthy fit).
#
# IMPORTANT DESIGN NOTE, unchanged from REMNB1_LEADNB2.R: RIGHT is defined
# as [anchor, total_n] -- i.e. anchor all the way to the most extreme end of
# the PC1-rank axis, not a small local window next to the anchor. LEFT is an
# equal-sized block immediately below anchor. This means RIGHT always
# contains the genuine extreme tail no matter where the anchor sits; moving
# the anchor down only dilutes RIGHT with more mid-rank genes. If the true
# NB2 signal is concentrated specifically in the extreme tail (supporting a
# real, localized transition near the established boundary), diff_alpha and
# the LRT statistic should visibly weaken as the anchor moves away from that
# boundary. If they stay roughly flat across the whole sweep, that supports
# the skeptical reading: a smooth global gradient, not a distinct boundary.
#
# Because RIGHT always reaches to total_n, LEFT (equal-sized) becomes
# infeasible once the anchor is much below roughly total_n / 2 (LEFT would
# need to extend below rank 1). The sweep range below is therefore
# restricted to anchors where both blocks fit.
#
# OUTPUTS:
#   Table_Anchor_Sweep.csv
#   Figure_Anchor_Sweep_<comparison>_<arm>.png
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
# the reference anchor is computed per arm below as total_n - 5000 + 1,
# matching REMNB1_LEADNB2.R's reference_rank convention exactly.
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
    # block still fits below rank 1. right_n = total_n - anchor + 1, and we
    # need anchor - right_n >= 1, i.e. anchor >= (total_n + 1) / 2
    # (approximately) -- solved numerically below per anchor instead of by
    # formula, to stay exactly consistent with lrt_at_anchor's own check.
    candidate_anchors <- seq(
      from = max(SWEEP_STEP, round(total_n * 0.5)),
      to = total_n - 100L,
      by = SWEEP_STEP
    )
    # Always include the established boundary itself as a reference point.
    candidate_anchors <- sort(unique(c(candidate_anchors, established_anchor)))

    message(sprintf(
      "[%s %s] total_n=%d, established_anchor=%d, sweeping %d candidate anchors",
      comparison_name, arm_name, total_n, established_anchor, length(candidate_anchors)
    ))

    for (anchor in candidate_anchors) {
      res <- lrt_at_anchor(count_mat_arm, mu_by_rank, rank_order, anchor, total_n)
      if (is.null(res)) next
      res$comparison <- comparison_name
      res$arm <- arm_name
      res$is_established_boundary <- (anchor == established_anchor)
      sweep_rows[[length(sweep_rows) + 1L]] <- res
    }
  }
}

sweep_table <- dplyr::bind_rows(sweep_rows)
write.csv(sweep_table, file.path(OUT_ROOT, "Table_Anchor_Sweep.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Figures: one per comparison, faceted by arm
# -----------------------------------------------------------------------------

for (comparison_name in names(COMPARISONS)) {
  sub <- sweep_table[sweep_table$comparison == comparison_name, , drop = FALSE]
  if (nrow(sub) == 0L) next

  marker_df <- sub[sub$is_established_boundary, , drop = FALSE]

  p <- ggplot(sub, aes(x = anchor_rank, y = diff_alpha)) +
    geom_line(color = "#3B6FA0", linewidth = 0.7) +
    geom_point(aes(shape = any_at_bound), size = 1.6, color = "#3B6FA0") +
    geom_point(data = marker_df, aes(x = anchor_rank, y = diff_alpha),
               color = "#C0392B", size = 3.2, shape = 17) +
    facet_wrap(~arm, ncol = 1, scales = "free_y") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
    labs(
      title = paste0(comparison_name, ": diff_alpha (RIGHT - LEFT) across the anchor sweep"),
      subtitle = "Red triangle = established leading-edge boundary. A flat curve across the sweep would support a smooth global gradient rather than a boundary-specific effect.",
      x = "Anchor rank",
      y = "alpha_right_mle - alpha_left_mle"
    ) +
    theme_bw(base_size = 11)

  ggsave(
    file.path(OUT_ROOT, paste0("Figure_Anchor_Sweep_", comparison_name, ".png")),
    p, width = 9, height = 7, dpi = 300
  )
}

message("Diagnostic complete. Outputs written to: ", OUT_ROOT)
message("  Table_Anchor_Sweep.csv")
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
