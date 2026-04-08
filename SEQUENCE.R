# =============================================================================
# SEQUENCE MANUSCRIPT PIPELINE
# MANUAL 5000-CUTOFF CORROBORATION VERSION
#
# Purpose:
# This version treats the leading-edge cutoff as a fixed manual choice
# (default = 5000 features). The variance-geometry analysis is then used to
# identify a local two-zero corroboration band that best supports that manual
# cutoff, rather than letting the zero-selection logic determine the cutoff
# itself. The geometric corroboration is therefore centered on the question:
# which nearby late-transition second-derivative zero pair best supports the
# chosen leading-edge boundary?
#
# Core interpretation:
# - The leading edge is defined manually as the top N ranked features.
# - The corroboration band is defined by two nearby eligible second-derivative
#   zero-crossings selected in the neighborhood of the manual cutoff.
# - The principal binary criterion is whether the post-cutoff region exhibits
#   greater median NB2 support than the pre-cutoff region.
# - Other NB support quantities are retained as continuous descriptive summaries.
# =============================================================================


# =============================================================================
# LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(fdrtool)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(gridExtra)
  library(grid)
  library(scales)
  library(grDevices)
  library(S4Vectors)
  library(BiocParallel)
})


# =============================================================================
# PRIMARY SETTINGS
# =============================================================================

alpha_level <- 0.20
lfc_boundary <- 1.0

manual_leading_edge_n <- 5000L

fdrtool_pct0 <- 0.75
fdr_clip_floor <- 1e-300
fdr_clip_ceiling <- 0.99
hc_threshold_upper_cap <- 0.95

rolling_window_fraction_d1 <- 0.015
rolling_window_fraction_d2 <- 0.025

manual_cutoff_search_left_fraction  <- 0.12
manual_cutoff_search_right_fraction <- 0.12

manual_cutoff_center_gap_fraction   <- 0.035
manual_cutoff_fallback_gap_fraction <- 0.07

near_cutoff_support_quantile <- 0.30

comparison_pairs <- list(
  RT0_ZT6   = c("R0", "ZT6"),
  RT2_ZT8   = c("R2", "ZT8"),
  RT4_ZT10  = c("R4", "ZT10"),
  RT8_ZT14  = c("R8", "ZT14")
)

input_csv <- "WTTS-Seq_2022_DE_raw_read_numbers.csv"
output_dir <- "exports/manual_5000_cutoff_version"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

clip_probabilities <- function(x,
                               floor_value = fdr_clip_floor,
                               ceiling_value = fdr_clip_ceiling) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  x <- pmin(pmax(x, floor_value), ceiling_value)
  x
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_quantile <- function(x, probs, na.rm = TRUE) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = probs, na.rm = na.rm, names = FALSE, type = 8))
}

rolling_mean_centered <- function(x, k) {
  n <- length(x)
  if (n == 0L) return(numeric(0))
  k <- max(1L, as.integer(k))
  if (k <= 1L) return(as.numeric(x))

  half <- floor(k / 2)
  out <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    lo <- max(1L, i - half)
    hi <- min(n, i + half)
    out[i] <- mean(x[lo:hi], na.rm = TRUE)
  }

  out
}

compute_centered_derivatives <- function(x, y,
                                         d1_window_fraction = rolling_window_fraction_d1,
                                         d2_window_fraction = rolling_window_fraction_d2) {
  n <- length(x)

  d1_raw <- c(NA_real_, diff(y) / diff(x))
  d1_raw[1] <- d1_raw[2]

  d1_window <- max(5L, floor(n * d1_window_fraction))
  if (d1_window %% 2L == 0L) d1_window <- d1_window + 1L
  d1 <- rolling_mean_centered(d1_raw, d1_window)

  d2_raw <- c(NA_real_, diff(d1) / diff(x))
  d2_raw[1] <- d2_raw[2]

  d2_window <- max(7L, floor(n * d2_window_fraction))
  if (d2_window %% 2L == 0L) d2_window <- d2_window + 1L
  d2 <- rolling_mean_centered(d2_raw, d2_window)

  list(
    d1 = d1,
    d2 = d2
  )
}

find_zero_crossings_with_sign <- function(y, x = seq_along(y)) {
  keep <- is.finite(y) & is.finite(x)
  y <- y[keep]
  x <- x[keep]

  if (length(y) < 2L) {
    return(data.frame(
      rank = numeric(0),
      idx_left = integer(0),
      idx_right = integer(0),
      y_left = numeric(0),
      y_right = numeric(0),
      crossing_type = character(0),
      stringsAsFactors = FALSE
    ))
  }

  s <- sign(y)

  for (i in seq_along(s)) {
    if (s[i] == 0) {
      left_nonzero <- if (i > 1L) tail(s[seq_len(i - 1L)][s[seq_len(i - 1L)] != 0], 1L) else numeric(0)
      right_nonzero <- if (i < length(s)) head(s[(i + 1L):length(s)][s[(i + 1L):length(s)] != 0], 1L) else numeric(0)

      if (length(left_nonzero)) {
        s[i] <- left_nonzero
      } else if (length(right_nonzero)) {
        s[i] <- right_nonzero
      }
    }
  }

  out <- list()
  j <- 1L

  for (i in seq_len(length(y) - 1L)) {
    if (!is.finite(y[i]) || !is.finite(y[i + 1L])) next

    sign_change <- sign(y[i]) != sign(y[i + 1L]) || y[i] == 0 || y[i + 1L] == 0
    if (!sign_change) next
    if ((y[i + 1L] - y[i]) == 0) next

    x0 <- x[i] - y[i] * (x[i + 1L] - x[i]) / (y[i + 1L] - y[i])

    crossing_type <- if (y[i] < 0 && y[i + 1L] > 0) {
      "neg_to_pos"
    } else if (y[i] > 0 && y[i + 1L] < 0) {
      "pos_to_neg"
    } else {
      "touch_or_flat"
    }

    out[[j]] <- data.frame(
      rank = x0,
      idx_left = i,
      idx_right = i + 1L,
      y_left = y[i],
      y_right = y[i + 1L],
      crossing_type = crossing_type,
      stringsAsFactors = FALSE
    )
    j <- j + 1L
  }

  if (!length(out)) {
    return(data.frame(
      rank = numeric(0),
      idx_left = integer(0),
      idx_right = integer(0),
      y_left = numeric(0),
      y_right = numeric(0),
      crossing_type = character(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, out)
}

compute_hc_threshold <- function(empirical_p,
                                 upper_cap = hc_threshold_upper_cap) {
  empirical_p <- empirical_p[is.finite(empirical_p)]
  empirical_p <- empirical_p[empirical_p > 0 & empirical_p < 1]

  if (!length(empirical_p)) return(NA_real_)

  empirical_p <- sort(empirical_p)
  hc <- tryCatch(
    fdrtool::hc.thresh(empirical_p),
    error = function(e) NA_real_
  )

  if (!is.finite(hc)) return(NA_real_)
  min(hc, upper_cap)
}

compute_empirical_p <- function(wald_z,
                                pct0 = fdrtool_pct0) {
  wald_z <- as.numeric(wald_z)
  keep <- is.finite(wald_z)

  out <- rep(NA_real_, length(wald_z))
  if (!any(keep)) return(out)

  fit <- tryCatch(
    fdrtool::fdrtool(wald_z[keep], statistic = "normal", plot = FALSE, verbose = FALSE, pct0 = pct0),
    error = function(e) NULL
  )

  if (is.null(fit)) return(out)

  if (!is.null(fit$param)) {
    if (!is.null(fit$param$pval)) {
      out[keep] <- fit$param$pval
    } else if (!is.null(fit$pval)) {
      out[keep] <- fit$pval
    }
  }

  if (all(!is.finite(out[keep])) && !is.null(fit$pval)) {
    out[keep] <- fit$pval
  }

  clip_probabilities(out)
}

make_combined_nb_support <- function(log_nb2, log_nb1, log_alpha_mu) {
  nb_gap <- log_nb2 - log_nb1

  z1 <- as.numeric(scale(log_nb2))
  z2 <- as.numeric(scale(nb_gap))
  z3 <- as.numeric(scale(log_alpha_mu))

  combined <- rowMeans(cbind(z1, z2, z3), na.rm = TRUE)
  combined[!is.finite(combined)] <- NA_real_
  combined
}


# =============================================================================
# MANUAL-5000 ZERO-PAIR CORROBORATION
#
# Manuscript method:
# The leading edge is defined manually as the first 5000 ranked features.
# Geometric corroboration is then obtained by selecting a local pair of eligible
# second-derivative zero-crossings in the neighborhood of the manual cutoff.
# The selected pair is the pair that most tightly and symmetrically brackets the
# manual cutoff while preserving proximity and local NB support.
# =============================================================================

select_zero_pair_for_manual_cutoff <- function(rank_axis,
                                               d2_curve,
                                               combined_nb_support,
                                               manual_cutoff_rank,
                                               search_left_fraction = manual_cutoff_search_left_fraction,
                                               search_right_fraction = manual_cutoff_search_right_fraction,
                                               center_gap_fraction = manual_cutoff_center_gap_fraction,
                                               fallback_gap_fraction = manual_cutoff_fallback_gap_fraction,
                                               support_quantile = near_cutoff_support_quantile) {
  n <- length(rank_axis)

  zc <- find_zero_crossings_with_sign(d2_curve, rank_axis)
  if (!nrow(zc)) {
    return(list(
      cutoff_center_rank = manual_cutoff_rank,
      terminal_start_rank = manual_cutoff_rank,
      cutoff_range_rank_min = manual_cutoff_rank,
      cutoff_range_rank_max = manual_cutoff_rank,
      zero_table = zc
    ))
  }

  zc$local_support <- vapply(zc$rank, function(rk) {
    idx <- which.min(abs(rank_axis - rk))
    combined_nb_support[idx]
  }, numeric(1))

  support_threshold <- safe_quantile(combined_nb_support, support_quantile)
  if (!is.finite(support_threshold)) support_threshold <- -Inf

  search_left_width  <- max(50L, floor(n * search_left_fraction))
  search_right_width <- max(50L, floor(n * search_right_fraction))

  left_bound  <- max(1L, manual_cutoff_rank - search_left_width)
  right_bound <- min(n, manual_cutoff_rank + search_right_width)

  zc$is_neg_to_pos <- zc$crossing_type == "neg_to_pos"
  zc$is_supported  <- is.finite(zc$local_support) & zc$local_support >= support_threshold
  zc$is_in_window  <- zc$rank >= left_bound & zc$rank <= right_bound

  candidates <- zc[zc$is_neg_to_pos & zc$is_in_window, , drop = FALSE]
  if (!nrow(candidates)) {
    candidates <- zc[zc$is_in_window, , drop = FALSE]
  }
  if (!nrow(candidates)) {
    candidates <- zc
  }

  left_candidates <- candidates[candidates$rank <= manual_cutoff_rank, , drop = FALSE]
  right_candidates <- candidates[candidates$rank >= manual_cutoff_rank, , drop = FALSE]

  left_candidates <- left_candidates[order(left_candidates$rank, decreasing = TRUE), , drop = FALSE]
  right_candidates <- right_candidates[order(right_candidates$rank, decreasing = FALSE), , drop = FALSE]

  max_gap <- max(50L, floor(n * center_gap_fraction))
  fallback_gap <- max(max_gap + 1L, floor(n * fallback_gap_fraction))

  left_supported <- left_candidates[left_candidates$is_supported, , drop = FALSE]
  right_supported <- right_candidates[right_candidates$is_supported, , drop = FALSE]

  if (nrow(left_supported) && nrow(right_supported)) {
    pair_grid <- expand.grid(
      left_i = seq_len(nrow(left_supported)),
      right_i = seq_len(nrow(right_supported))
    )
    pair_grid$left_rank <- left_supported$rank[pair_grid$left_i]
    pair_grid$right_rank <- right_supported$rank[pair_grid$right_i]
    pair_grid$width <- pair_grid$right_rank - pair_grid$left_rank
    pair_grid$center <- (pair_grid$left_rank + pair_grid$right_rank) / 2
    pair_grid$cutoff_offset <- abs(pair_grid$center - manual_cutoff_rank)
    pair_grid$support_score <- left_supported$local_support[pair_grid$left_i] +
      right_supported$local_support[pair_grid$right_i]

    pair_grid <- pair_grid[pair_grid$width > 0, , drop = FALSE]

    close_pairs <- pair_grid[pair_grid$width <= max_gap, , drop = FALSE]
    if (!nrow(close_pairs)) {
      close_pairs <- pair_grid[pair_grid$width <= fallback_gap, , drop = FALSE]
    }
    if (!nrow(close_pairs)) {
      close_pairs <- pair_grid
    }

    close_pairs <- close_pairs[order(close_pairs$cutoff_offset,
                                     close_pairs$width,
                                     -close_pairs$support_score), , drop = FALSE]

    best_pair <- close_pairs[1, , drop = FALSE]

    left_rank <- best_pair$left_rank[1]
    right_rank <- best_pair$right_rank[1]
  } else {
    left_rank <- if (nrow(left_candidates)) left_candidates$rank[1] else manual_cutoff_rank
    right_rank <- if (nrow(right_candidates)) right_candidates$rank[1] else manual_cutoff_rank
  }

  list(
    cutoff_center_rank = left_rank,
    terminal_start_rank = right_rank,
    cutoff_range_rank_min = min(left_rank, right_rank),
    cutoff_range_rank_max = max(left_rank, right_rank),
    zero_table = zc
  )
}


# =============================================================================
# NB2 CORROBORATION SUMMARY
#
# Manuscript method:
# The principal binary corroboration criterion is whether the post-cutoff region
# exhibits greater median NB2 support than the pre-cutoff region. Other support
# quantities are reported numerically.
# =============================================================================

build_nb2_summary <- function(rank_axis,
                              cutoff_center_rank,
                              terminal_start_rank,
                              cutoff_range_rank_min,
                              cutoff_range_rank_max,
                              manual_leading_edge_n,
                              combined_nb_support,
                              log_nb2,
                              log_nb1,
                              log_alpha_mu) {
  pre_idx <- which(rank_axis < cutoff_center_rank)
  post_idx <- which(rank_axis >= cutoff_center_rank)

  if (!length(pre_idx)) pre_idx <- 1L
  if (!length(post_idx)) post_idx <- length(rank_axis)

  nb_gap <- log_nb2 - log_nb1

  pre_median_log_nb2 <- safe_median(log_nb2[pre_idx])
  post_median_log_nb2 <- safe_median(log_nb2[post_idx])

  data.frame(
    manual_leading_edge_n = manual_leading_edge_n,
    cutoff_center_rank = cutoff_center_rank,
    terminal_start_rank = terminal_start_rank,
    cutoff_range_rank_min = cutoff_range_rank_min,
    cutoff_range_rank_max = cutoff_range_rank_max,

    pre_cutoff_size = length(pre_idx),
    post_cutoff_size = length(post_idx),

    pre_cutoff_median_log_nb2 = pre_median_log_nb2,
    post_cutoff_median_log_nb2 = post_median_log_nb2,
    post_minus_pre_log_nb2 = post_median_log_nb2 - pre_median_log_nb2,

    pre_cutoff_median_nb2_nb1_contrast = safe_median(nb_gap[pre_idx]),
    post_cutoff_median_nb2_nb1_contrast = safe_median(nb_gap[post_idx]),
    post_minus_pre_nb2_nb1_contrast =
      safe_median(nb_gap[post_idx]) - safe_median(nb_gap[pre_idx]),

    pre_cutoff_median_log_alpha_mu = safe_median(log_alpha_mu[pre_idx]),
    post_cutoff_median_log_alpha_mu = safe_median(log_alpha_mu[post_idx]),
    post_minus_pre_log_alpha_mu =
      safe_median(log_alpha_mu[post_idx]) - safe_median(log_alpha_mu[pre_idx]),

    pre_cutoff_mean_combined_nb_support = safe_mean(combined_nb_support[pre_idx]),
    post_cutoff_mean_combined_nb_support = safe_mean(combined_nb_support[post_idx]),
    post_minus_pre_mean_combined_nb_support =
      safe_mean(combined_nb_support[post_idx]) - safe_mean(combined_nb_support[pre_idx]),

    post_cutoff_more_nb2 = post_median_log_nb2 > pre_median_log_nb2,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# PLOT HELPERS
# =============================================================================

make_annotation_block <- function(summary_row, label_prefix) {
  paste(
    sprintf("%s", label_prefix),
    sprintf("Manual leading edge = %s", summary_row$manual_leading_edge_n),
    sprintf("Cutoff center rank = %s", summary_row$cutoff_center_rank),
    sprintf("Terminal start rank = %s", summary_row$terminal_start_rank),
    sprintf("Cutoff range = [%s, %s]",
            summary_row$cutoff_range_rank_min,
            summary_row$cutoff_range_rank_max),
    sprintf("Pre med log(NB2) = %.3f", summary_row$pre_cutoff_median_log_nb2),
    sprintf("Post med log(NB2) = %.3f", summary_row$post_cutoff_median_log_nb2),
    sprintf("Pre med NB2-NB1 = %.3f", summary_row$pre_cutoff_median_nb2_nb1_contrast),
    sprintf("Post med NB2-NB1 = %.3f", summary_row$post_cutoff_median_nb2_nb1_contrast),
    sprintf("Pre med log(alpha*mu) = %.3f", summary_row$pre_cutoff_median_log_alpha_mu),
    sprintf("Post med log(alpha*mu) = %.3f", summary_row$post_cutoff_median_log_alpha_mu),
    sprintf("Post-cutoff more NB2 = %s",
            ifelse(isTRUE(summary_row$post_cutoff_more_nb2), "TRUE", "FALSE")),
    sep = "\n"
  )
}

plot_variance_geometry_panel <- function(df_plot,
                                         comparison_label,
                                         group_label,
                                         summary_row,
                                         annotation_text,
                                         out_file) {
  p <- ggplot(df_plot, aes(x = rank, y = fitted_log_variance)) +
    annotate("rect",
             xmin = summary_row$cutoff_range_rank_min,
             xmax = summary_row$cutoff_range_rank_max,
             ymin = -Inf,
             ymax = Inf,
             alpha = 0.12) +
    geom_line(linewidth = 0.7) +
    geom_vline(xintercept = summary_row$cutoff_center_rank, linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = summary_row$terminal_start_rank, linetype = "dotted", linewidth = 0.5) +
    labs(
      title = paste(comparison_label, group_label, "variance geometry"),
      subtitle = "Manual leading-edge cutoff corroborated by nearby second-derivative zero pair",
      x = "Rank",
      y = "Smoothed log(variance + 1)"
    ) +
    annotate("text",
             x = Inf, y = Inf,
             label = annotation_text,
             hjust = 1.02, vjust = 1.02,
             size = 3) +
    theme_bw(base_size = 10)

  ggsave(out_file, p, width = 10, height = 7, dpi = 300)
}

plot_second_derivative_panel <- function(df_plot,
                                         comparison_label,
                                         group_label,
                                         summary_row,
                                         zero_table,
                                         out_file) {
  p <- ggplot(df_plot, aes(x = rank, y = d2)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_vline(xintercept = summary_row$cutoff_center_rank, linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = summary_row$terminal_start_rank, linetype = "dotted", linewidth = 0.5) +
    annotate("rect",
             xmin = summary_row$cutoff_range_rank_min,
             xmax = summary_row$cutoff_range_rank_max,
             ymin = -Inf,
             ymax = Inf,
             alpha = 0.12) +
    geom_point(
      data = zero_table,
      aes(x = rank, y = 0),
      inherit.aes = FALSE,
      size = 1.5
    ) +
    labs(
      title = paste(comparison_label, group_label, "second derivative"),
      subtitle = "Eligible zero-crossings in the neighborhood of the manual cutoff",
      x = "Rank",
      y = "Smoothed second derivative"
    ) +
    theme_bw(base_size = 10)

  ggsave(out_file, p, width = 10, height = 7, dpi = 300)
}

plot_nb_support_panel <- function(df_plot,
                                  comparison_label,
                                  group_label,
                                  summary_row,
                                  out_file) {
  p <- ggplot(df_plot, aes(x = rank)) +
    geom_line(aes(y = log_nb2, linetype = "log(NB2)"), linewidth = 0.6) +
    geom_line(aes(y = nb_gap, linetype = "NB2-NB1"), linewidth = 0.6) +
    geom_line(aes(y = log_alpha_mu, linetype = "log(alpha*mu)"), linewidth = 0.6) +
    geom_vline(xintercept = summary_row$cutoff_center_rank, linetype = "dashed", linewidth = 0.5) +
    labs(
      title = paste(comparison_label, group_label, "NB corroboration"),
      subtitle = "Post-cutoff corroboration is evaluated against the pre-cutoff region",
      x = "Rank",
      y = "Support metric",
      linetype = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top")

  ggsave(out_file, p, width = 10, height = 7, dpi = 300)
}


# =============================================================================
# INPUT DATA
# =============================================================================

count_df <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)

feature_col <- colnames(count_df)[1]
count_matrix <- as.matrix(count_df[, -1, drop = FALSE])
rownames(count_matrix) <- count_df[[feature_col]]

sample_names <- colnames(count_matrix)

sample_info <- data.frame(
  sample = sample_names,
  condition = ifelse(grepl("^ZT", sample_names), "trt", "untrt"),
  time_group = sub("_.*$", "", sample_names),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample


# =============================================================================
# PER-COMPARISON ANALYSIS
# =============================================================================

overall_summary_list <- list()
per_feature_export_list <- list()

for (comparison_label in names(comparison_pairs)) {

  message("Processing comparison: ", comparison_label)

  pair_tokens <- comparison_pairs[[comparison_label]]
  left_prefix <- pair_tokens[1]
  right_prefix <- pair_tokens[2]

  keep_samples <- grepl(paste0("^", left_prefix, "_"), sample_names) |
    grepl(paste0("^", right_prefix, "_"), sample_names)

  counts_sub <- count_matrix[, keep_samples, drop = FALSE]
  coldata_sub <- sample_info[colnames(counts_sub), , drop = FALSE]
  coldata_sub$condition <- factor(coldata_sub$condition, levels = c("untrt", "trt"))

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_sub),
    colData = coldata_sub,
    design = ~ condition
  )

  keep_gene <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_gene, ]

  dds <- DESeq(dds, parallel = FALSE)

  res <- results(dds, alpha = alpha_level)
  coef_name <- grep("condition_trt_vs_untrt", resultsNames(dds), value = TRUE)[1]

  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")

  res_df <- as.data.frame(res)
  res_shrunk_df <- as.data.frame(res_shrunk)
  norm_counts <- counts(dds, normalized = TRUE)

  res_df$feature_id <- rownames(res_df)
  res_shrunk_df$feature_id <- rownames(res_shrunk_df)

  merged_df <- res_df %>%
    select(feature_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    rename(
      log2FoldChange_raw = log2FoldChange,
      stat_wald = stat
    ) %>%
    left_join(
      res_shrunk_df %>%
        select(feature_id, log2FoldChange) %>%
        rename(log2FoldChange_shrunk = log2FoldChange),
      by = "feature_id"
    )

  merged_df$empirical_p <- compute_empirical_p(merged_df$stat_wald)
  merged_df$empirical_p <- clip_probabilities(merged_df$empirical_p)

  hc_threshold_dataset <- compute_hc_threshold(merged_df$empirical_p)

  if (!is.finite(hc_threshold_dataset)) {
    hc_threshold_dataset <- 0.95
  }

  hbfss_threshold <- abs(log10(hc_threshold_dataset)) * lfc_boundary
  merged_df$HBFSS <- abs(merged_df$log2FoldChange_shrunk * log10(merged_df$empirical_p))

  merged_df$deseq2_significant <- is.finite(merged_df$padj) & merged_df$padj < alpha_level
  merged_df$strong_effect <- merged_df$deseq2_significant &
    is.finite(merged_df$log2FoldChange_shrunk) &
    abs(merged_df$log2FoldChange_shrunk) >= lfc_boundary
  merged_df$weak_effect <- merged_df$deseq2_significant &
    is.finite(merged_df$log2FoldChange_shrunk) &
    abs(merged_df$log2FoldChange_shrunk) < lfc_boundary
  merged_df$hbfss_significant <- is.finite(merged_df$empirical_p) &
    merged_df$empirical_p < hc_threshold_dataset &
    is.finite(merged_df$HBFSS) &
    merged_df$HBFSS >= hbfss_threshold

  vst_mat <- assay(varianceStabilizingTransformation(dds, blind = TRUE))
  row_variances <- apply(vst_mat, 1, var, na.rm = TRUE)

  # ---------------------------------------------------------------------------
  # Build group-specific EVS-style ranked geometry for control and treatment
  # ---------------------------------------------------------------------------

  for (group_label in c("control", "treatment")) {

    group_samples <- if (group_label == "control") {
      rownames(coldata_sub)[coldata_sub$condition == "untrt"]
    } else {
      rownames(coldata_sub)[coldata_sub$condition == "trt"]
    }

    group_counts <- norm_counts[, group_samples, drop = FALSE]

    if (ncol(group_counts) < 2L) next

    group_log <- log2(group_counts + 1)

    pca <- prcomp(group_log, center = TRUE, scale. = TRUE)
    abs_pc1_loading <- abs(pca$rotation[, 1])

    ranked_idx <- order(abs_pc1_loading, decreasing = TRUE)
    ranked_feature_id <- colnames(group_log)[ranked_idx]

    if (is.null(ranked_feature_id)) {
      ranked_feature_id <- rownames(group_log)[ranked_idx]
    }

    feature_ids_ranked <- rownames(group_log)[ranked_idx]

    ranked_var <- apply(group_log[ranked_idx, , drop = FALSE], 1, var, na.rm = TRUE)
    rank_axis <- seq_along(ranked_var)

    fitted_log_variance <- log1p(ranked_var)

    derivs <- compute_centered_derivatives(rank_axis, fitted_log_variance)
    d1 <- derivs$d1
    d2 <- derivs$d2

    feature_match <- match(rownames(group_log)[ranked_idx], merged_df$feature_id)

    ranked_base_mean <- merged_df$baseMean[feature_match]
    ranked_log_nb2 <- log1p(ranked_base_mean + pmax(ranked_var, 0))
    ranked_log_nb1 <- log1p(ranked_base_mean + sqrt(pmax(ranked_base_mean, 0)))
    ranked_log_alpha_mu <- log1p(pmax(ranked_base_mean, 0) * pmax(ranked_var, 0))

    combined_nb_support <- make_combined_nb_support(
      log_nb2 = ranked_log_nb2,
      log_nb1 = ranked_log_nb1,
      log_alpha_mu = ranked_log_alpha_mu
    )

    zero_pair <- select_zero_pair_for_manual_cutoff(
      rank_axis = rank_axis,
      d2_curve = d2,
      combined_nb_support = combined_nb_support,
      manual_cutoff_rank = min(manual_leading_edge_n, length(rank_axis))
    )

    cutoff_center_rank <- zero_pair$cutoff_center_rank
    terminal_start_rank <- zero_pair$terminal_start_rank
    cutoff_range_rank_min <- zero_pair$cutoff_range_rank_min
    cutoff_range_rank_max <- zero_pair$cutoff_range_rank_max

    summary_row <- build_nb2_summary(
      rank_axis = rank_axis,
      cutoff_center_rank = cutoff_center_rank,
      terminal_start_rank = terminal_start_rank,
      cutoff_range_rank_min = cutoff_range_rank_min,
      cutoff_range_rank_max = cutoff_range_rank_max,
      manual_leading_edge_n = min(manual_leading_edge_n, length(rank_axis)),
      combined_nb_support = combined_nb_support,
      log_nb2 = ranked_log_nb2,
      log_nb1 = ranked_log_nb1,
      log_alpha_mu = ranked_log_alpha_mu
    )

    summary_row$comparison <- comparison_label
    summary_row$group <- group_label
    summary_row$hc_threshold_dataset <- hc_threshold_dataset
    summary_row$hbfss_threshold <- hbfss_threshold

    overall_summary_list[[paste(comparison_label, group_label, sep = "__")]] <- summary_row

    per_feature_df <- data.frame(
      comparison = comparison_label,
      group = group_label,
      feature_id = rownames(group_log)[ranked_idx],
      rank = rank_axis,
      fitted_log_variance = fitted_log_variance,
      d1 = d1,
      d2 = d2,
      log_nb2 = ranked_log_nb2,
      log_nb1 = ranked_log_nb1,
      nb_gap = ranked_log_nb2 - ranked_log_nb1,
      log_alpha_mu = ranked_log_alpha_mu,
      combined_nb_support = combined_nb_support,
      stringsAsFactors = FALSE
    )

    per_feature_export_list[[paste(comparison_label, group_label, sep = "__")]] <- per_feature_df

    annotation_text <- make_annotation_block(
      summary_row = summary_row,
      label_prefix = paste(comparison_label, group_label)
    )

    variance_plot_file <- file.path(
      output_dir,
      paste0(comparison_label, "_", group_label, "_variance_geometry.png")
    )

    d2_plot_file <- file.path(
      output_dir,
      paste0(comparison_label, "_", group_label, "_second_derivative.png")
    )

    nb_plot_file <- file.path(
      output_dir,
      paste0(comparison_label, "_", group_label, "_nb_support.png")
    )

    plot_variance_geometry_panel(
      df_plot = per_feature_df,
      comparison_label = comparison_label,
      group_label = group_label,
      summary_row = summary_row,
      annotation_text = annotation_text,
      out_file = variance_plot_file
    )

    plot_second_derivative_panel(
      df_plot = per_feature_df,
      comparison_label = comparison_label,
      group_label = group_label,
      summary_row = summary_row,
      zero_table = zero_pair$zero_table,
      out_file = d2_plot_file
    )

    plot_nb_support_panel(
      df_plot = per_feature_df,
      comparison_label = comparison_label,
      group_label = group_label,
      summary_row = summary_row,
      out_file = nb_plot_file
    )
  }

  # ---------------------------------------------------------------------------
  # Comparison-level DE results export
  # ---------------------------------------------------------------------------

  merged_df$comparison <- comparison_label
  merged_df$hc_threshold_dataset <- hc_threshold_dataset
  merged_df$hbfss_threshold <- hbfss_threshold

  utils::write.csv(
    merged_df,
    file = file.path(output_dir, paste0(comparison_label, "_DE_results.csv")),
    row.names = FALSE
  )
}


# =============================================================================
# FINAL EXPORTS
# =============================================================================

overall_cutoff_summary <- bind_rows(overall_summary_list)
overall_cutoff_summary <- overall_cutoff_summary %>%
  select(
    comparison,
    group,
    manual_leading_edge_n,
    cutoff_center_rank,
    terminal_start_rank,
    cutoff_range_rank_min,
    cutoff_range_rank_max,
    pre_cutoff_size,
    post_cutoff_size,
    pre_cutoff_median_log_nb2,
    post_cutoff_median_log_nb2,
    post_minus_pre_log_nb2,
    pre_cutoff_median_nb2_nb1_contrast,
    post_cutoff_median_nb2_nb1_contrast,
    post_minus_pre_nb2_nb1_contrast,
    pre_cutoff_median_log_alpha_mu,
    post_cutoff_median_log_alpha_mu,
    post_minus_pre_log_alpha_mu,
    pre_cutoff_mean_combined_nb_support,
    post_cutoff_mean_combined_nb_support,
    post_minus_pre_mean_combined_nb_support,
    post_cutoff_more_nb2,
    hc_threshold_dataset,
    hbfss_threshold
  )

utils::write.csv(
  overall_cutoff_summary,
  file = file.path(output_dir, "overall_cutoff_summary.csv"),
  row.names = FALSE
)

per_feature_all <- bind_rows(per_feature_export_list)

utils::write.csv(
  per_feature_all,
  file = file.path(output_dir, "per_feature_rank_geometry.csv"),
  row.names = FALSE
)


# =============================================================================
# METHODS TEXT FOR MANUSCRIPT MIRRORING
# =============================================================================
# Features were ranked within each comparison and group using the magnitude of
# the first principal-component loading derived from log2(x + 1)-transformed
# normalized counts. The leading edge was defined manually as the top 5000
# ranked features. Smoothed variance geometry was then used only to corroborate
# this manual boundary rather than to determine it. Specifically, the smoothed
# second derivative of the ranked variance curve was evaluated in the local
# neighborhood of the manual cutoff, and a nearby pair of eligible zero-
# crossings was selected to provide a compact geometric transition band that
# best bracketed the chosen leading-edge boundary while preserving local
# negative-binomial support. Negative-binomial corroboration was summarized on
# both sides of the cutoff using median log(NB2), median NB2-minus-NB1 contrast,
# median log(alpha·mu), and mean combined NB support. The principal binary
# corroboration criterion was whether the post-cutoff region exhibited greater
# median NB2 support than the pre-cutoff region.
# =============================================================================
