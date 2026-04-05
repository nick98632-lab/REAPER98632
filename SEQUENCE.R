#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

# =============================================================================
# STAGE 1 NB PARTS ANALYSIS
# Empirical negative-binomial decomposition along the EVS-ranked axis
#
# Core identities:
#   Var(X) = mu + alpha * mu^2
#   IOD    = Var(X) / mu   = 1 + alpha * mu
#   CV^2   = Var(X) / mu^2 = 1 / mu + alpha
#
# This script estimates mu and alpha empirically from the replicate counts,
# builds an EVS rank using absolute PC1 loadings, and then plots the two NB
# parts directly across rank and percentile summaries.
# =============================================================================

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
input_dir  <- "."
output_dir <- file.path("exports", paste0("nb_parts_run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

count_file <- file.path(input_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
if (!file.exists(count_file)) {
  stop("Count file not found: ", count_file, call. = FALSE)
}

# ----------------------------------------------------------------------------
# Embedded metadata copied from the working SEQUENCE file
# ----------------------------------------------------------------------------
meta_all <- data.frame(
  id = c(
    "R0_1","R0_2","R0_3","R0_4","R0_5","ZT6_1","ZT6_2","ZT6_3","ZT6_4","ZT6_5",
    "R2_1","R2_2","R2_3","R2_4","R2_5","ZT8_1","ZT8_2","ZT8_3","ZT8_4","ZT8_5",
    "R4_1","R4_2","R4_3","R4_4","R4_5","ZT10_1","ZT10_2","ZT10_3","ZT10_4","ZT10_5",
    "R8_1","R8_2","R8_3","R8_4","R8_5","ZT14_1","ZT14_2","ZT14_3","ZT14_4","ZT14_5"
  ),
  condition = c(
    "treatment","treatment","treatment","treatment","treatment",
    "control","control","control","control","control",
    "treatment","treatment","treatment","treatment","treatment",
    "control","control","control","control","control",
    "treatment","treatment","treatment","treatment","treatment",
    "control","control","control","control","control",
    "treatment","treatment","treatment","treatment","treatment",
    "control","control","control","control","control"
  ),
  stringsAsFactors = FALSE
)
rownames(meta_all) <- meta_all$id
meta_all$condition <- factor(meta_all$condition, levels = c("control", "treatment"))
levels(meta_all$condition) <- c("untrt", "trt")

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix   = c("R0", "R2", "R4", "R8"),
  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

safe_name_vector <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

read_count_matrix <- function(path) {
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(dat) < 2) stop("Count file has fewer than 2 columns.", call. = FALSE)

  feature_col <- 1L
  feature_ids <- safe_name_vector(dat[[feature_col]])
  if (anyDuplicated(feature_ids)) {
    dup_idx <- duplicated(feature_ids)
    feature_ids[dup_idx] <- paste0(feature_ids[dup_idx], "__dup", seq_len(sum(dup_idx)))
  }

  count_df <- dat[, -feature_col, drop = FALSE]
  sample_cols <- intersect(colnames(count_df), meta_all$id)
  if (length(sample_cols) == 0) {
    stop("No sample columns from meta_all were found in the count file.", call. = FALSE)
  }

  count_df <- count_df[, sample_cols, drop = FALSE]
  count_mat <- as.matrix(count_df)
  storage.mode(count_mat) <- "numeric"
  rownames(count_mat) <- feature_ids
  count_mat
}

subset_comparison <- function(count_mat, comparison_name) {
  row <- comparison_table[comparison_table$comparison_name == comparison_name, , drop = FALSE]
  if (nrow(row) != 1) stop("Unknown comparison: ", comparison_name, call. = FALSE)

  prefixes <- c(row$group1_prefix, row$group2_prefix)
  keep_ids <- meta_all$id[Reduce(`|`, lapply(prefixes, function(px) startsWith(meta_all$id, px)))]
  keep_ids <- intersect(keep_ids, colnames(count_mat))
  if (length(keep_ids) == 0) stop("No columns found for comparison: ", comparison_name, call. = FALSE)

  count_sub <- count_mat[, keep_ids, drop = FALSE]
  coldata   <- meta_all[keep_ids, , drop = FALSE]
  stopifnot(identical(colnames(count_sub), rownames(coldata)))
  list(count = count_sub, coldata = coldata)
}

compute_size_factors_median_ratio <- function(count_mat) {
  log_counts <- log(count_mat)
  log_counts[!is.finite(log_counts)] <- NA_real_
  geom_means <- exp(rowMeans(log_counts, na.rm = TRUE))
  valid <- is.finite(geom_means) & geom_means > 0
  if (!any(valid)) return(setNames(rep(1, ncol(count_mat)), colnames(count_mat)))
  ratios <- sweep(count_mat[valid, , drop = FALSE], 1, geom_means[valid], "/")
  sf <- apply(ratios, 2, function(x) median(x[is.finite(x) & x > 0], na.rm = TRUE))
  sf[!is.finite(sf) | sf <= 0] <- 1
  sf / exp(mean(log(sf)))
}

normalize_counts <- function(count_mat, size_factors) {
  sweep(count_mat, 2, size_factors[colnames(count_mat)], "/")
}

compute_condition_feature_metrics <- function(norm_mat) {
  mu  <- rowMeans(norm_mat, na.rm = TRUE)
  v   <- apply(norm_mat, 1, var, na.rm = TRUE)
  iod <- ifelse(mu > 0, v / mu, NA_real_)
  cv2 <- ifelse(mu > 0, v / (mu ^ 2), NA_real_)
  alpha_hat <- ifelse(mu > 0, pmax((v - mu) / (mu ^ 2), 0), NA_real_)

  data.frame(
    feature_id = rownames(norm_mat),
    mean_norm  = mu,
    var_norm   = v,
    iod_emp    = iod,
    cv2_emp    = cv2,
    alpha_hat  = alpha_hat,
    iod_nb     = ifelse(is.finite(alpha_hat) & is.finite(mu), 1 + alpha_hat * mu, NA_real_),
    cv2_nb     = ifelse(is.finite(alpha_hat) & is.finite(mu) & mu > 0, (1 / mu) + alpha_hat, NA_real_),
    poisson_part = ifelse(is.finite(mu) & mu > 0, 1 / mu, NA_real_),
    stringsAsFactors = FALSE
  )
}

compute_pc1_abs_loadings <- function(norm_mat) {
  x <- log2(norm_mat + 1)
  x <- x[rowSums(is.finite(x)) == ncol(x), , drop = FALSE]
  if (nrow(x) < 2) stop("Not enough finite rows for PCA.", call. = FALSE)
  x_centered <- sweep(x, 1, rowMeans(x), "-")
  pca <- prcomp(t(x_centered), center = TRUE, scale. = FALSE)
  # feature loadings on PC1 from rotation in feature space via SVD equivalent
  s <- svd(x_centered)
  loadings <- s$u[, 1]
  names(loadings) <- rownames(x_centered)
  abs(loadings)
}

build_evs_rank_table <- function(count_sub, coldata) {
  trt_ids   <- rownames(coldata)[coldata$condition == "trt"]
  untrt_ids <- rownames(coldata)[coldata$condition == "untrt"]

  sf        <- compute_size_factors_median_ratio(count_sub)
  norm_mat  <- normalize_counts(count_sub, sf)

  trt_metrics   <- compute_condition_feature_metrics(norm_mat[, trt_ids, drop = FALSE])
  untrt_metrics <- compute_condition_feature_metrics(norm_mat[, untrt_ids, drop = FALSE])

  trt_load   <- compute_pc1_abs_loadings(norm_mat[, trt_ids, drop = FALSE])
  untrt_load <- compute_pc1_abs_loadings(norm_mat[, untrt_ids, drop = FALSE])

  trt_metrics$abs_pc1_loading   <- trt_load[trt_metrics$feature_id]
  untrt_metrics$abs_pc1_loading <- untrt_load[untrt_metrics$feature_id]

  trt_metrics$rank_trt     <- rank(-trt_metrics$abs_pc1_loading, ties.method = "first")
  untrt_metrics$rank_untrt <- rank(-untrt_metrics$abs_pc1_loading, ties.method = "first")

  merged <- merge(trt_metrics, untrt_metrics,
                  by = "feature_id", all = TRUE,
                  suffixes = c("_trt", "_untrt"), sort = FALSE)

  merged$combined_rank <- pmin(merged$rank_trt %||% Inf, merged$rank_untrt %||% Inf, na.rm = TRUE)
  merged <- merged[order(merged$combined_rank, merged$feature_id), , drop = FALSE]
  merged$combined_rank <- seq_len(nrow(merged))
  rownames(merged) <- NULL
  list(rank_table = merged, norm_mat = norm_mat, size_factors = sf)
}

percentile_summary <- function(rank_tbl, step = 0.01, window = 0.12) {
  n <- nrow(rank_tbl)
  centers <- seq(step, 1, by = step)
  out <- vector("list", length(centers))

  for (i in seq_along(centers)) {
    p <- centers[i]
    center_rank <- max(1L, min(n, round(p * n)))
    half_width  <- max(1L, round((window * n) / 2))
    lo <- max(1L, center_rank - half_width)
    hi <- min(n, center_rank + half_width)
    idx <- lo:hi
    sub <- rank_tbl[idx, , drop = FALSE]

    mean_nonmissing <- function(x) {
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      mean(x)
    }

    out[[i]] <- data.frame(
      percentile = p,
      center_rank = center_rank,
      lo_rank = lo,
      hi_rank = hi,
      n_window = length(idx),
      mu_trt = mean_nonmissing(sub$mean_norm_trt),
      mu_untrt = mean_nonmissing(sub$mean_norm_untrt),
      alpha_trt = mean_nonmissing(sub$alpha_hat_trt),
      alpha_untrt = mean_nonmissing(sub$alpha_hat_untrt),
      iod_emp_trt = mean_nonmissing(sub$iod_emp_trt),
      iod_emp_untrt = mean_nonmissing(sub$iod_emp_untrt),
      cv2_emp_trt = mean_nonmissing(sub$cv2_emp_trt),
      cv2_emp_untrt = mean_nonmissing(sub$cv2_emp_untrt),
      iod_nb_trt = mean_nonmissing(sub$iod_nb_trt),
      iod_nb_untrt = mean_nonmissing(sub$iod_nb_untrt),
      cv2_nb_trt = mean_nonmissing(sub$cv2_nb_trt),
      cv2_nb_untrt = mean_nonmissing(sub$cv2_nb_untrt),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

find_last_crossing_before_divergence <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2) return(NULL)

  s <- sign(y)
  s[s == 0] <- NA
  for (i in seq_along(s)) {
    if (is.na(s[i])) {
      prev <- if (i > 1) na.omit(s[seq_len(i - 1)]) else numeric(0)
      nextv <- if (i < length(s)) na.omit(s[(i + 1):length(s)]) else numeric(0)
      s[i] <- if (length(prev)) tail(prev, 1) else if (length(nextv)) nextv[1] else 0
    }
  }

  cross_idx <- which(diff(s) != 0)
  if (!length(cross_idx)) return(NULL)
  i <- tail(cross_idx, 1)
  x1 <- x[i]; x2 <- x[i + 1]
  y1 <- y[i]; y2 <- y[i + 1]
  xr <- if (isTRUE(all.equal(y1, y2))) x2 else x1 - y1 * (x2 - x1) / (y2 - y1)
  list(index = i, x = xr, x_left = x1, x_right = x2)
}

plot_nb_curves <- function(df, comparison_name, suffix, out_dir) {
  mkline <- function(title, y1, y2, lab1, lab2, filename, crossing = TRUE) {
    dd <- data.frame(percentile = df$percentile, y1 = y1, y2 = y2)
    g <- ggplot(dd, aes(percentile)) +
      geom_line(aes(y = y1, color = lab1), linewidth = 0.9) +
      geom_line(aes(y = y2, color = lab2), linewidth = 0.9) +
      labs(title = paste0(comparison_name, " — ", title), x = "Percentile", y = "Value", color = NULL) +
      theme_bw(base_size = 11)
    if (crossing) {
      d <- dd$y1 - dd$y2
      cr <- find_last_crossing_before_divergence(dd$percentile, d)
      if (!is.null(cr)) {
        g <- g + geom_vline(xintercept = cr$x, linetype = "dashed") +
          annotate("text", x = cr$x, y = max(c(dd$y1, dd$y2), na.rm = TRUE),
                   label = sprintf("Last crossing = %.3f", cr$x), vjust = -0.4, size = 3)
      }
    }
    ggsave(file.path(out_dir, filename), g, width = 9, height = 5, dpi = 300)
  }

  mkline("Treatment NB parts", df$iod_nb_trt, df$cv2_nb_trt,
         "IOD = 1 + alpha*mu", "CV^2 = 1/mu + alpha",
         paste0(comparison_name, "_", suffix, "_treatment_nb_parts.png"))

  mkline("Control NB parts", df$iod_nb_untrt, df$cv2_nb_untrt,
         "IOD = 1 + alpha*mu", "CV^2 = 1/mu + alpha",
         paste0(comparison_name, "_", suffix, "_control_nb_parts.png"))

  df$iod_nb_comb <- rowMeans(cbind(df$iod_nb_trt, df$iod_nb_untrt), na.rm = TRUE)
  df$cv2_nb_comb <- rowMeans(cbind(df$cv2_nb_trt, df$cv2_nb_untrt), na.rm = TRUE)

  mkline("Combined NB parts", df$iod_nb_comb, df$cv2_nb_comb,
         "Mean IOD NB part", "Mean CV^2 NB part",
         paste0(comparison_name, "_", suffix, "_combined_nb_parts.png"))

  diff_df <- data.frame(
    percentile = df$percentile,
    treatment_diff = df$iod_nb_trt - df$cv2_nb_trt,
    control_diff = df$iod_nb_untrt - df$cv2_nb_untrt,
    combined_diff = df$iod_nb_comb - df$cv2_nb_comb
  )
  gdiff <- ggplot(diff_df, aes(percentile)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(aes(y = treatment_diff, color = "Treatment"), linewidth = 0.9) +
    geom_line(aes(y = control_diff, color = "Control"), linewidth = 0.9) +
    geom_line(aes(y = combined_diff, color = "Combined"), linewidth = 0.9) +
    labs(title = paste0(comparison_name, " — NB difference curves"),
         x = "Percentile", y = "IOD - CV^2", color = NULL) +
    theme_bw(base_size = 11)
  cr <- find_last_crossing_before_divergence(diff_df$percentile, diff_df$combined_diff)
  if (!is.null(cr)) {
    gdiff <- gdiff + geom_vline(xintercept = cr$x, linetype = "dashed") +
      annotate("text", x = cr$x, y = max(diff_df$combined_diff, na.rm = TRUE),
               label = sprintf("Last combined crossing = %.3f", cr$x), vjust = -0.4, size = 3)
  }
  ggsave(file.path(out_dir, paste0(comparison_name, "_", suffix, "_difference_curves.png")), gdiff,
         width = 9, height = 5, dpi = 300)

  invisible(df)
}

# ----------------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------------
message("Reading count matrix...")
count_mat <- read_count_matrix(count_file)

summary_rows <- list()
for (cmp in comparison_table$comparison_name) {
  message("Running comparison: ", cmp)
  cmp_dir <- file.path(output_dir, cmp)
  dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

  sub_obj <- subset_comparison(count_mat, cmp)
  evs_obj <- build_evs_rank_table(sub_obj$count, sub_obj$coldata)
  rank_tbl <- evs_obj$rank_table
  write.csv(rank_tbl, file.path(cmp_dir, paste0(cmp, "_evs_rank_table.csv")), row.names = FALSE)

  sum_df <- percentile_summary(rank_tbl, step = 0.01, window = 0.12)
  write.csv(sum_df, file.path(cmp_dir, paste0(cmp, "_nb_percentile_summary.csv")), row.names = FALSE)

  plot_nb_curves(sum_df, cmp, "nb", cmp_dir)

  combined_diff <- rowMeans(cbind(sum_df$iod_nb_trt, sum_df$iod_nb_untrt), na.rm = TRUE) -
    rowMeans(cbind(sum_df$cv2_nb_trt, sum_df$cv2_nb_untrt), na.rm = TRUE)
  cr <- find_last_crossing_before_divergence(sum_df$percentile, combined_diff)
  summary_rows[[cmp]] <- data.frame(
    comparison = cmp,
    crossing_percentile = if (is.null(cr)) NA_real_ else cr$x,
    crossing_rank = if (is.null(cr)) NA_real_ else round(cr$x * nrow(rank_tbl)),
    n_features = nrow(rank_tbl),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(output_dir, "nb_crossing_summary.csv"), row.names = FALSE)
message("Done. Outputs written to: ", output_dir)
