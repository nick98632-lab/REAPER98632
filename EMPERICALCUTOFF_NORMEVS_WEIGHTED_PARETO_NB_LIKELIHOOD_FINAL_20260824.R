#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# NORMEVS FINAL PIPELINE
# DESeq2 median-of-ratios -> log1p(normalized counts) -> arm-specific PC1
# -> size-factor-aware excess-variance geometry -> comparison-specific knots
# -> weighted Pareto top-k selection -> post-selection NB1/NB2 corroboration
# =============================================================================

# ------------------------------- SETTINGS -----------------------------------
COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT <- "/root/REAPER98632/exports/normevs_final"

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

BENEFIT_WEIGHT <- 1.0
CONTAMINATION_WEIGHT <- 1.0
KNOT_COARSE_GRID_POINTS <- 4096L
DISPLAY_VAR_SPAR <- 0.72
DISPLAY_D_SPAR <- 0.72
DISPLAY_NB_SPAR <- 0.68
PNG_DPI <- 360
EXPORT_PDF <- TRUE
NB_LOG_ALPHA_LOWER <- -14
NB_LOG_ALPHA_UPPER <- 8

COL <- list(
  control = "#386CB0",
  treatment = "#159D91",
  remainder = "#DCE6F2",
  interval = "#FFF0B3",
  leading = "#D8F3E7",
  c1 = "#D73027",
  c2 = "#1A9850",
  selected = "#B5179E",
  candidate = "#A7A7A7",
  pareto = "#5E3C99",
  nb1 = "#E69F00",
  nb2 = "#0072B2",
  left = "#5B8FD1",
  right = "#43A047"
)

# ----------------------------- OUTPUT SETUP ---------------------------------
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
SUMMARY_FIG_DIR <- file.path(OUT_ROOT, "Summary", "Figures")
SUMMARY_TAB_DIR <- file.path(OUT_ROOT, "Summary", "Tables")
dir.create(SUMMARY_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ BASIC HELPERS -------------------------------
theme_manuscript <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.3),
      plot.subtitle = element_text(size = base_size - 0.2, margin = margin(b = 5)),
      plot.caption = element_text(size = base_size - 1.5, color = "grey25", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#222222"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(color = "#B7B7B7", fill = NA, linewidth = 0.45),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(8, 10, 8, 10)
    )
}

save_figure <- function(plot_obj, png_path, width = 15, height = 11) {
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(png_path, plot_obj, width = width, height = height, units = "in",
         dpi = PNG_DPI, bg = "white", limitsize = FALSE)
  if (isTRUE(EXPORT_PDF)) {
    pdf_path <- sub("\\.png$", ".pdf", png_path, ignore.case = TRUE)
    ggsave(pdf_path, plot_obj, width = width, height = height, units = "in",
           device = "pdf", bg = "white", limitsize = FALSE)
  }
  invisible(png_path)
}

save_grid_2x2 <- function(plots, png_path, width = 16, height = 12.5) {
  stopifnot(length(plots) == 4L)
  draw_once <- function(device_fun) {
    device_fun()
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 2)))
    for (i in seq_along(plots)) {
      r <- if (i <= 2L) 1L else 2L
      c <- if (i %% 2L == 1L) 1L else 2L
      print(plots[[i]], vp = viewport(layout.pos.row = r, layout.pos.col = c))
    }
    dev.off()
  }
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  draw_once(function() png(png_path, width = width, height = height, units = "in",
                           res = PNG_DPI, bg = "white"))
  if (isTRUE(EXPORT_PDF)) {
    pdf_path <- sub("\\.png$", ".pdf", png_path, ignore.case = TRUE)
    draw_once(function() pdf(pdf_path, width = width, height = height,
                             onefile = TRUE, useDingbats = FALSE))
  }
  invisible(png_path)
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

rank_cutoff_from_k <- function(N, k) {
  if (!is.finite(k) || k < 1L || k > N) return(NA_integer_)
  as.integer(N - k + 1L)
}

verify_outputs <- function(paths) {
  bad <- paths[!file.exists(paths) | is.na(file.info(paths)$size) | file.info(paths)$size <= 0]
  if (length(bad)) stop("Missing or empty expected outputs:\n", paste(bad, collapse = "\n"))
  invisible(TRUE)
}

# ------------------------------ DATA INPUT ----------------------------------
read_count_data <- function(path, group_patterns) {
  if (!file.exists(path)) stop("Count file does not exist: ", path)
  raw_df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(raw_df) < 1L || ncol(raw_df) < 2L) stop("Count file is empty or malformed.")

  sample_idx <- sort(unique(unlist(lapply(group_patterns, function(p) grep(p, colnames(raw_df))))))
  if (!length(sample_idx)) stop("No sample columns matched GROUP_PATTERNS.")
  if (1L %in% sample_idx) stop("Column 1 must contain feature IDs, not sample counts.")

  feature_id <- trimws(as.character(raw_df[[1L]]))
  blank <- is.na(feature_id) | feature_id == ""
  feature_id[blank] <- paste0("__feature_row_", which(blank))
  feature_id <- make.unique(feature_id, sep = "__dup_")

  non_sample_idx <- setdiff(seq_len(ncol(raw_df)), sample_idx)
  symbol_candidates <- non_sample_idx[
    tolower(colnames(raw_df)[non_sample_idx]) %in% c("symbol", "gene_symbol", "genesymbol", "gene")
  ]
  if (length(symbol_candidates)) {
    gene_symbol <- as.character(raw_df[[symbol_candidates[1L]]])
  } else if (length(non_sample_idx) >= 2L) {
    gene_symbol <- as.character(raw_df[[non_sample_idx[2L]]])
  } else {
    gene_symbol <- rep(NA_character_, nrow(raw_df))
  }

  count_df <- raw_df[, sample_idx, drop = FALSE]
  count_mat <- do.call(cbind, lapply(count_df, function(x) suppressWarnings(as.numeric(trimws(as.character(x))))))
  colnames(count_mat) <- colnames(count_df)
  rownames(count_mat) <- feature_id
  storage.mode(count_mat) <- "numeric"
  count_mat[!is.finite(count_mat)] <- 0
  count_mat <- pmax(count_mat, 0)

  annotation <- data.frame(feature_id = feature_id, gene_symbol = gene_symbol, stringsAsFactors = FALSE)
  list(counts = count_mat, annotation = annotation)
}

# --------------------------- DESEQ2 NORMALIZATION ----------------------------
normalize_deseq2_comparison <- function(count_mat, group_labels) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("DESeq2 is required.")
  col_data <- data.frame(group = factor(group_labels), row.names = colnames(count_mat))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(count_mat), colData = col_data, design = ~ group)

  # Strict median-of-ratios normalization. If this fails because no feature has
  # a positive geometric mean, stop rather than silently switching methods.
  dds <- DESeq2::estimateSizeFactors(dds, type = "ratio")
  sf <- DESeq2::sizeFactors(dds)
  if (any(!is.finite(sf)) || any(sf <= 0)) stop("Invalid DESeq2 size factor.")

  list(
    normalized_counts = DESeq2::counts(dds, normalized = TRUE),
    size_factors = sf
  )
}

# -------------------------- PC1 / VARIANCE GEOMETRY --------------------------
compute_pc1_rank <- function(normalized_counts_arm) {
  # FINAL EVS INPUT: DESeq2 median-of-ratios normalized counts, then log1p.
  rank_matrix <- log1p(normalized_counts_arm)

  pca <- prcomp(t(rank_matrix), center = TRUE, scale. = FALSE, rank. = 1)
  loading <- pca$rotation[, 1L]
  loading[!is.finite(loading)] <- 0
  abs_loading <- abs(loading)
  rank_order <- order(abs_loading, decreasing = FALSE)
  lambda1 <- pca$sdev[1L]^2
  P <- lambda1 * loading^2

  list(
    loading = loading,
    abs_loading = abs_loading,
    rank_order = rank_order,
    lambda1 = lambda1,
    pc1_variance_contribution = P
  )
}

compute_raw_variance_geometry <- function(raw_counts_arm, rank_order) {
  raw_var <- apply(raw_counts_arm, 1L, var, na.rm = TRUE)
  raw_var[!is.finite(raw_var)] <- 0
  raw_var <- pmax(raw_var, 0)
  raw_var_ranked <- raw_var[rank_order]
  rank <- seq_along(rank_order)
  y <- log1p(raw_var_ranked)
  fit <- smooth.spline(rank, y, spar = DISPLAY_VAR_SPAR)
  yhat <- as.numeric(predict(fit, x = rank, deriv = 0)$y)
  data.frame(
    rank = rank,
    raw_empirical_variance = raw_var_ranked,
    log1p_raw_empirical_variance = y,
    display_log1p_raw_empirical_variance = yhat,
    stringsAsFactors = FALSE
  )
}

smooth_divergence_for_display <- function(rank, D, spar = DISPLAY_D_SPAR) {
  fit <- smooth.spline(rank, D, spar = spar)
  y <- as.numeric(predict(fit, x = rank, deriv = 0)$y)
  y - seq(y[1L], y[length(y)], length.out = length(y))
}

compute_group_analysis <- function(group_name, raw_counts_arm, normalized_counts_arm, size_factors_arm) {
  pc1 <- compute_pc1_rank(normalized_counts_arm)
  o <- pc1$rank_order
  geometry <- compute_raw_variance_geometry(raw_counts_arm, o)

  mu <- rowMeans(normalized_counts_arm, na.rm = TRUE)
  v <- apply(normalized_counts_arm, 1L, var, na.rm = TRUE)
  mu[!is.finite(mu)] <- 0
  v[!is.finite(v)] <- 0
  mu <- pmax(mu, 0)
  v <- pmax(v, 0)

  inv_sf_mean <- mean(1 / as.numeric(size_factors_arm))
  poisson_baseline <- mu * inv_sf_mean
  excess <- pmax(v - poisson_baseline, 0)

  P_ranked <- pc1$pc1_variance_contribution[o]
  E_ranked <- excess[o]
  P_total <- sum(P_ranked)
  E_total <- sum(E_ranked)
  if (!is.finite(P_total) || P_total <= 0) stop("PC1 mass undefined for group ", group_name)
  if (!is.finite(E_total) || E_total <= 0) stop("Excess-variance mass undefined for group ", group_name)

  p_mass <- P_ranked / P_total
  q_mass <- E_ranked / E_total
  F_P <- cumsum(p_mass)
  F_E <- cumsum(q_mass)
  D <- F_E - F_P
  rank <- seq_along(o)

  df <- geometry %>%
    mutate(
      group = group_name,
      feature_id = rownames(raw_counts_arm)[o],
      pc1_loading = pc1$loading[o],
      abs_pc1_loading = pc1$abs_loading[o],
      pc1_eigenvalue = pc1$lambda1,
      pc1_variance_contribution = P_ranked,
      pc1_variance_mass = p_mass,
      normalized_group_mean = mu[o],
      normalized_group_variance = v[o],
      mean_inverse_size_factor = inv_sf_mean,
      poisson_variance_baseline_normalized = poisson_baseline[o],
      excess_variance = E_ranked,
      excess_variance_mass = q_mass,
      cumulative_pc1_mass = F_P,
      cumulative_excess_variance_mass = F_E,
      cumulative_divergence = D,
      display_D = smooth_divergence_for_display(rank, D)
    )

  list(data = df, rank_order = o)
}

# ---------------------------- SHARED TWO-KNOT FIT ----------------------------
piecewise_basis <- function(x, c1, c2) {
  cbind(intercept = 1, x = x, hinge1 = pmax(x - c1, 0), hinge2 = pmax(x - c2, 0))
}

piecewise_sse <- function(par, x, D_mat, min_gap) {
  c1 <- par[1L]; c2 <- par[2L]
  if (!is.finite(c1) || !is.finite(c2) || c1 <= 0 || c2 >= 1 || c2 - c1 <= min_gap) return(1e100)
  X <- piecewise_basis(x, c1, c2)
  coef <- tryCatch(qr.coef(qr(X), D_mat), error = function(e) NULL)
  if (is.null(coef) || any(!is.finite(coef))) return(1e100)
  resid <- D_mat - X %*% coef
  sum(resid^2)
}

fit_shared_knots <- function(group_results) {
  groups <- names(group_results)
  N_values <- vapply(group_results, function(z) nrow(z$data), integer(1))
  if (length(unique(N_values)) != 1L) stop("Arms must have the same feature universe.")
  N <- N_values[1L]
  x_full <- (seq_len(N) - 1) / (N - 1)
  D_full <- do.call(cbind, lapply(group_results, function(z) z$data$cumulative_divergence))
  colnames(D_full) <- groups

  opt_n <- min(KNOT_COARSE_GRID_POINTS, N)
  opt_idx <- unique(as.integer(round(seq(1, N, length.out = opt_n))))
  x_opt <- x_full[opt_idx]
  D_opt <- D_full[opt_idx, , drop = FALSE]
  min_gap <- max(4 / (N - 1), .Machine$double.eps^0.25)
  starts <- list(c(0.03,0.97), c(0.08,0.92), c(0.15,0.85), c(0.25,0.75), c(0.35,0.65))

  coarse <- lapply(starts, function(start) {
    optim(start, piecewise_sse, x = x_opt, D_mat = D_opt, min_gap = min_gap,
          method = "Nelder-Mead", control = list(maxit = 700, reltol = 1e-11))
  })
  best <- coarse[[which.min(vapply(coarse, function(z) z$value, numeric(1)))]]
  refined <- optim(best$par, piecewise_sse, x = x_full, D_mat = D_full, min_gap = min_gap,
                   method = "Nelder-Mead", control = list(maxit = 1000, reltol = 1e-12))
  if (!is.finite(refined$value)) stop("Shared-knot optimization failed.")

  c1 <- as.integer(round(1 + refined$par[1L] * (N - 1)))
  c2 <- as.integer(round(1 + refined$par[2L] * (N - 1)))
  c1 <- max(2L, min(N - 2L, c1))
  c2 <- max(c1 + 1L, min(N - 1L, c2))

  c1x <- (c1 - 1) / (N - 1)
  c2x <- (c2 - 1) / (N - 1)
  X <- piecewise_basis(x_full, c1x, c2x)
  coef <- qr.coef(qr(X), D_full)
  fitted <- X %*% coef

  list(c1 = c1, c2 = c2, x = x_full, fitted = fitted,
       SSE = sum((D_full - fitted)^2), groups = groups)
}

# ------------------------------ PARETO SCAN ---------------------------------
make_rank_map <- function(df) setNames(df$rank, df$feature_id)

rank_to_region <- function(rank, c1, c2) {
  ifelse(rank < c1, "Remainder", ifelse(rank <= c2, "Divergence", "LeadingEdge"))
}

cumulative_activation <- function(depth, K) {
  depth <- as.integer(depth)
  keep <- is.finite(depth) & depth >= 1L & depth <= K
  if (!any(keep)) return(rep(0L, K))
  cumsum(tabulate(depth[keep], nbins = K))
}

active_interval_count <- function(starts, ends, K) {
  if (!length(starts)) return(rep(0L, K))
  starts <- as.integer(starts); ends <- as.integer(ends)
  valid <- is.finite(starts) & is.finite(ends) & starts >= 1L & starts <= K & ends > starts
  starts <- starts[valid]; ends <- pmin(ends[valid], K + 1L)
  if (!length(starts)) return(rep(0L, K))
  d <- integer(K + 1L)
  d <- d + tabulate(starts, nbins = K + 1L) - tabulate(ends, nbins = K + 1L)
  cumsum(d)[seq_len(K)]
}

scan_pair_cutoffs <- function(control_df, treatment_df, c1, c2, comparison_name, control_group, treatment_group) {
  if (nrow(control_df) != nrow(treatment_df)) stop("Arm feature counts differ.")
  if (!setequal(control_df$feature_id, treatment_df$feature_id)) stop("Arm feature IDs differ.")
  N <- nrow(control_df)
  K <- as.integer(N - c2)
  if (K < 1L) stop("No candidate top-k depth exists beyond c2.")

  ids <- control_df$feature_id
  rC <- as.integer(unname(make_rank_map(control_df)[ids]))
  rT <- as.integer(unname(make_rank_map(treatment_df)[ids]))
  dC <- N - rC + 1L
  dT <- N - rT + 1L
  k <- seq_len(K)

  joint_n <- cumulative_activation(pmax(dC, dT), K)

  idx_le_C <- which(dC < dT & dC <= K & dT <= K)
  idx_le_T <- which(dT < dC & dT <= K & dC <= K)
  disjoint_control_opposite_le_n <- active_interval_count(dC[idx_le_C], dT[idx_le_C], K)
  disjoint_treatment_opposite_le_n <- active_interval_count(dT[idx_le_T], dC[idx_le_T], K)

  idx_div_C <- which(dC <= K & rT >= c1 & rT <= c2)
  idx_div_T <- which(dT <= K & rC >= c1 & rC <= c2)
  disjoint_control_opposite_divergence_n <- cumulative_activation(dC[idx_div_C], K)
  disjoint_treatment_opposite_divergence_n <- cumulative_activation(dT[idx_div_T], K)

  idx_rem_C <- which(dC <= K & rT < c1)
  idx_rem_T <- which(dT <= K & rC < c1)
  remainder_cross_control_n <- cumulative_activation(dC[idx_rem_C], K)
  remainder_cross_treatment_n <- cumulative_activation(dT[idx_rem_T], K)

  disjoint_opposite_le_n <- disjoint_control_opposite_le_n + disjoint_treatment_opposite_le_n
  disjoint_opposite_divergence_n <- disjoint_control_opposite_divergence_n + disjoint_treatment_opposite_divergence_n
  permissible_disjoint_n <- disjoint_opposite_le_n + disjoint_opposite_divergence_n
  good_n <- joint_n + permissible_disjoint_n
  remainder_cross_n <- remainder_cross_control_n + remainder_cross_treatment_n
  union_n <- good_n + remainder_cross_n
  union_check <- cumulative_activation(pmin(dC, dT), K)
  if (!all(union_n == union_check)) stop("Internal union-count mismatch in ", comparison_name)

  data.frame(
    comparison = comparison_name,
    control_group = control_group,
    treatment_group = treatment_group,
    k = k,
    cutoff_rank = N - k + 1L,
    joint_n = joint_n,
    disjoint_control_opposite_le_n = disjoint_control_opposite_le_n,
    disjoint_treatment_opposite_le_n = disjoint_treatment_opposite_le_n,
    disjoint_opposite_le_n = disjoint_opposite_le_n,
    disjoint_control_opposite_divergence_n = disjoint_control_opposite_divergence_n,
    disjoint_treatment_opposite_divergence_n = disjoint_treatment_opposite_divergence_n,
    disjoint_opposite_divergence_n = disjoint_opposite_divergence_n,
    permissible_disjoint_n = permissible_disjoint_n,
    good_n = good_n,
    remainder_cross_control_n = remainder_cross_control_n,
    remainder_cross_treatment_n = remainder_cross_treatment_n,
    remainder_cross_n = remainder_cross_n,
    union_n = union_n,
    retained_fraction = ifelse(union_n > 0, good_n / union_n, NA_real_),
    remainder_cross_fraction = ifelse(union_n > 0, remainder_cross_n / union_n, NA_real_),
    jaccard_top_k = ifelse(union_n > 0, joint_n / union_n, NA_real_),
    stringsAsFactors = FALSE
  )
}

mark_pareto_frontier <- function(scan_df, good_col = "good_n", cost_col = "remainder_cross_n") {
  tmp <- scan_df %>%
    transmute(row_id = row_number(), k = k, good = .data[[good_col]], cost = .data[[cost_col]]) %>%
    arrange(cost, desc(good), desc(k)) %>%
    group_by(cost) %>% slice(1L) %>% ungroup() %>%
    arrange(cost, desc(good))
  running_best_before <- c(-Inf, head(cummax(tmp$good), -1L))
  frontier <- tmp %>% filter(good > running_best_before) %>% arrange(cost, good, k)
  key_all <- paste(scan_df[[cost_col]], scan_df[[good_col]], sep = "::")
  key_frontier <- paste(frontier$cost, frontier$good, sep = "::")
  out <- scan_df
  out$is_pareto <- key_all %in% key_frontier
  list(scan = out, frontier = frontier)
}

select_weighted_pareto_optimum <- function(scan_df, good_col = "good_n", cost_col = "remainder_cross_n",
                                            benefit_weight = BENEFIT_WEIGHT,
                                            contamination_weight = CONTAMINATION_WEIGHT) {
  if (benefit_weight < 0 || contamination_weight < 0 || benefit_weight + contamination_weight <= 0)
    stop("Pareto weights must be non-negative and not both zero.")

  marked <- mark_pareto_frontier(scan_df, good_col, cost_col)
  frontier <- marked$frontier
  if (!nrow(frontier)) stop("No Pareto-optimal points.")

  gr <- range(frontier$good, na.rm = TRUE)
  cr <- range(frontier$cost, na.rm = TRUE)
  ng <- function(x) if (diff(gr) == 0) rep(1, length(x)) else (x - gr[1L]) / diff(gr)
  nr <- function(x) if (diff(cr) == 0) rep(0, length(x)) else (x - cr[1L]) / diff(cr)

  frontier$good_norm <- ng(frontier$good)
  frontier$remainder_norm <- nr(frontier$cost)
  frontier$weighted_utility <- benefit_weight * frontier$good_norm - contamination_weight * frontier$remainder_norm
  best <- max(frontier$weighted_utility, na.rm = TRUE)
  chosen <- frontier %>%
    filter(abs(weighted_utility - best) < 1e-12) %>%
    arrange(desc(good), cost, desc(k)) %>% slice(1L)

  out <- marked$scan
  out$good_norm <- ng(out[[good_col]])
  out$remainder_norm <- nr(out[[cost_col]])
  out$weighted_utility <- benefit_weight * out$good_norm - contamination_weight * out$remainder_norm
  out$is_selected_weighted <- out$k == chosen$k[1L]

  list(
    scan = out,
    frontier = frontier,
    selected_k = as.integer(chosen$k[1L]),
    selected_good = chosen$good[1L],
    selected_cost = chosen$cost[1L],
    selected_good_norm = chosen$good_norm[1L],
    selected_remainder_norm = chosen$remainder_norm[1L],
    selected_utility = chosen$weighted_utility[1L]
  )
}

# -------------------------- FINAL EVS CLASSIFICATION -------------------------
classify_pair_at_k <- function(control_df, treatment_df, k, c1, c2,
                               comparison_name, control_group, treatment_group) {
  N <- nrow(control_df)
  if (k < 1L || k > N - c2) stop("Selected k is outside the admissible domain.")

  control_top <- tail(control_df$feature_id, k)
  treatment_top <- tail(treatment_df$feature_id, k)
  union_ids <- union(control_top, treatment_top)
  rC <- as.integer(unname(make_rank_map(control_df)[union_ids]))
  rT <- as.integer(unname(make_rank_map(treatment_df)[union_ids]))
  inC <- union_ids %in% control_top
  inT <- union_ids %in% treatment_top
  joint <- inC & inT
  c_only <- inC & !inT
  t_only <- inT & !inC

  control_region <- rank_to_region(rC, c1, c2)
  treatment_region <- rank_to_region(rT, c1, c2)
  opposite_region <- rep(NA_character_, length(union_ids))
  opposite_region[c_only] <- treatment_region[c_only]
  opposite_region[t_only] <- control_region[t_only]

  remainder_cost <- (c_only & rT < c1) | (t_only & rC < c1)
  opposite_div <- (c_only & rT >= c1 & rT <= c2) | (t_only & rC >= c1 & rC <= c2)
  opposite_le <- (c_only & rT > c2) | (t_only & rC > c2)

  analysis_class <- ifelse(
    joint, "Joint",
    ifelse(
      remainder_cost,
      ifelse(c_only, paste0("Disjoint_", control_group, "_OppositeRemainder_Cost"),
             paste0("Disjoint_", treatment_group, "_OppositeRemainder_Cost")),
      ifelse(
        opposite_div,
        ifelse(c_only, paste0("Disjoint_", control_group, "_OppositeDivergence"),
               paste0("Disjoint_", treatment_group, "_OppositeDivergence")),
        ifelse(c_only, paste0("Disjoint_", control_group, "_OppositeLeadingEdge"),
               paste0("Disjoint_", treatment_group, "_OppositeLeadingEdge"))
      )
    )
  )

  data.frame(
    comparison = comparison_name,
    feature_id = union_ids,
    selected_k = as.integer(k),
    cutoff_rank = rank_cutoff_from_k(N, k),
    control_group = control_group,
    treatment_group = treatment_group,
    base_class = ifelse(joint, "Joint", ifelse(c_only, paste0("Disjoint_", control_group), paste0("Disjoint_", treatment_group))),
    analysis_class = analysis_class,
    control_rank = rC,
    treatment_rank = rT,
    control_region = control_region,
    treatment_region = treatment_region,
    opposite_region = opposite_region,
    control_top_k = inC,
    treatment_top_k = inT,
    disjoint_opposite_leading_edge = opposite_le,
    disjoint_opposite_divergence = opposite_div,
    pareto_remainder_crossing_cost = remainder_cost,
    evs_membership = "LeadingEdge",
    stringsAsFactors = FALSE
  )
}

# ----------------------- POST-SELECTION NB CORROBORATION ---------------------
smooth_metric_for_display <- function(rank, value, spar = DISPLAY_NB_SPAR) {
  ok <- is.finite(rank) & is.finite(value)
  out <- rep(NA_real_, length(value))
  if (sum(ok) < 8L) return(out)
  fit <- tryCatch(smooth.spline(rank[ok], value[ok], spar = spar), error = function(e) NULL)
  if (is.null(fit)) return(out)
  out[ok] <- as.numeric(predict(fit, x = rank[ok], deriv = 0)$y)
  out
}

compute_ranked_nb_metrics <- function(normalized_counts_arm, size_factors_arm, ranked_feature_ids) {
  x <- normalized_counts_arm[ranked_feature_ids, , drop = FALSE]
  mu <- rowMeans(x, na.rm = TRUE)
  v <- apply(x, 1L, var, na.rm = TRUE)
  mu[!is.finite(mu)] <- 0; v[!is.finite(v)] <- 0
  mu <- pmax(mu, 0); v <- pmax(v, 0)
  inv_sf_mean <- mean(1 / as.numeric(size_factors_arm))
  poisson_baseline <- mu * inv_sf_mean
  excess <- pmax(v - poisson_baseline, 0)
  alpha_hat <- rep(0, length(mu))
  pos <- mu > 0
  alpha_hat[pos] <- pmax(excess[pos] / (mu[pos]^2), 0)

  data.frame(
    rank = seq_along(ranked_feature_ids),
    feature_id = ranked_feature_ids,
    normalized_mean = mu,
    normalized_variance = v,
    poisson_variance_baseline = poisson_baseline,
    excess_variance_signal = log1p(excess),
    alpha_hat_nb2 = alpha_hat,
    alpha_mu_signal = log1p(alpha_hat * mu),
    stringsAsFactors = FALSE
  )
}

assign_matched_regions <- function(rank_df, k) {
  N <- nrow(rank_df)
  right_start <- rank_cutoff_from_k(N, k)
  left_end <- right_start - 1L
  left_start <- left_end - k + 1L
  if (left_start < 1L) {
    rank_df$corroboration_region <- "Other"
    rank_df$corroboration_region[rank_df$rank >= right_start] <- "RIGHT"
    attr(rank_df, "matched_left_available") <- FALSE
    return(rank_df)
  }
  rank_df$corroboration_region <- "Other"
  rank_df$corroboration_region[rank_df$rank >= left_start & rank_df$rank <= left_end] <- "LEFT"
  rank_df$corroboration_region[rank_df$rank >= right_start] <- "RIGHT"
  attr(rank_df, "matched_left_available") <- TRUE
  rank_df
}

fit_nb_dispersion_model <- function(y, mu, model = c("NB1", "NB2")) {
  model <- match.arg(model)
  y <- as.numeric(y); mu <- as.numeric(mu)
  ok <- is.finite(y) & is.finite(mu) & y >= 0 & mu > 0
  y <- y[ok]; mu <- pmax(mu[ok], 1e-10)
  if (length(y) < 10L) return(list(model = model, alpha = NA_real_, logLik = NA_real_, n_obs = length(y)))

  nll <- function(log_alpha) {
    alpha <- exp(log_alpha)
    size <- if (model == "NB2") rep(1 / alpha, length(mu)) else mu / alpha
    size <- pmax(size, 1e-10)
    ll <- suppressWarnings(dnbinom(round(y), mu = mu, size = size, log = TRUE))
    if (any(!is.finite(ll))) return(1e100)
    -sum(ll)
  }
  opt <- optimize(nll, c(NB_LOG_ALPHA_LOWER, NB_LOG_ALPHA_UPPER))
  list(model = model, alpha = exp(opt$minimum), logLik = -opt$objective, n_obs = length(y))
}

fit_nb1_nb2_region <- function(raw_counts_arm, normalized_counts_arm, size_factors_arm,
                               feature_ids, arm_label, region_label) {
  raw_sub <- raw_counts_arm[feature_ids, , drop = FALSE]
  norm_sub <- normalized_counts_arm[feature_ids, , drop = FALSE]
  mu_norm <- rowMeans(norm_sub, na.rm = TRUE)
  mu_norm[!is.finite(mu_norm)] <- 0
  mu_norm <- pmax(mu_norm, 1e-10)
  sf <- as.numeric(size_factors_arm[colnames(raw_sub)])
  if (any(!is.finite(sf)) || any(sf <= 0)) stop("Invalid size factor in NB likelihood.")
  mu_mat <- outer(mu_norm, sf, "*")
  y <- as.vector(raw_sub)
  mu <- as.vector(mu_mat)

  nb1 <- fit_nb_dispersion_model(y, mu, "NB1")
  nb2 <- fit_nb_dispersion_model(y, mu, "NB2")
  delta <- nb2$logLik - nb1$logLik

  data.frame(
    arm = arm_label,
    region = region_label,
    n_features = length(feature_ids),
    n_observations = nb1$n_obs,
    alpha_NB1 = nb1$alpha,
    alpha_NB2 = nb2$alpha,
    logLik_NB1 = nb1$logLik,
    logLik_NB2 = nb2$logLik,
    two_delta_logLik_NB2_minus_NB1 = 2 * delta,
    log10_likelihood_ratio_NB2_to_NB1 = delta / log(10),
    preferred_model = ifelse(is.finite(delta) & delta > 0, "NB2",
                             ifelse(is.finite(delta) & delta < 0, "NB1", "Tie")),
    stringsAsFactors = FALSE
  )
}

build_corroboration_for_arm <- function(arm_label, raw_counts_arm, normalized_counts_arm,
                                        size_factors_arm, rank_df, k) {
  nb_df <- compute_ranked_nb_metrics(normalized_counts_arm, size_factors_arm, rank_df$feature_id)
  nb_df$abs_pc1_loading <- rank_df$abs_pc1_loading
  nb_df <- assign_matched_regions(nb_df, k)
  nb_df$excess_variance_signal_smooth <- smooth_metric_for_display(nb_df$rank, nb_df$excess_variance_signal)
  nb_df$alpha_mu_signal_smooth <- smooth_metric_for_display(nb_df$rank, nb_df$alpha_mu_signal)

  regions <- if (isTRUE(attr(nb_df, "matched_left_available"))) c("LEFT", "RIGHT") else "RIGHT"
  region_summary <- nb_df %>%
    filter(corroboration_region %in% regions) %>%
    group_by(corroboration_region) %>%
    summarise(
      n = n(),
      median_excess_variance_signal = median(excess_variance_signal, na.rm = TRUE),
      median_alpha_mu_signal = median(alpha_mu_signal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(arm = arm_label)

  likelihood <- bind_rows(lapply(regions, function(region_name) {
    ids <- nb_df$feature_id[nb_df$corroboration_region == region_name]
    fit_nb1_nb2_region(raw_counts_arm, normalized_counts_arm, size_factors_arm,
                       ids, arm_label, region_name)
  }))

  list(rank_data = nb_df, region_summary = region_summary, likelihood = likelihood,
       matched_left_available = isTRUE(attr(nb_df, "matched_left_available")))
}

# ------------------------------- FIGURES ------------------------------------
add_rank_regions <- function(p, c1, c2, N) {
  p +
    annotate("rect", xmin = 1, xmax = c1, ymin = -Inf, ymax = Inf, fill = COL$remainder, alpha = 0.42) +
    annotate("rect", xmin = c1, xmax = c2, ymin = -Inf, ymax = Inf, fill = COL$interval, alpha = 0.38) +
    annotate("rect", xmin = c2, xmax = N, ymin = -Inf, ymax = Inf, fill = COL$leading, alpha = 0.40)
}

make_pareto_panel <- function(scan_df, selected_k, comparison_name) {
  selected <- scan_df %>% filter(k == selected_k) %>% slice(1L)
  frontier <- scan_df %>% filter(is_pareto) %>% arrange(remainder_cross_n, good_n, k) %>%
    distinct(remainder_cross_n, good_n, .keep_all = TRUE)
  ggplot() +
    geom_path(data = scan_df, aes(remainder_cross_n, good_n, group = 1), color = COL$candidate, linewidth = 0.45, alpha = 0.55) +
    geom_point(data = scan_df, aes(remainder_cross_n, good_n), color = COL$candidate, size = 1.0, alpha = 0.45) +
    geom_path(data = frontier, aes(remainder_cross_n, good_n, group = 1), color = COL$pareto, linewidth = 1.25) +
    geom_point(data = selected, aes(remainder_cross_n, good_n), shape = 23, fill = COL$selected,
               color = COL$selected, size = 4.7, stroke = 1.1) +
    annotate("label", x = selected$remainder_cross_n, y = selected$good_n,
             label = paste0("k* = ", selected_k, "\nG = ", selected$good_n, "\nR = ", selected$remainder_cross_n,
                            "\nU = ", formatC(selected$weighted_utility, digits = 4, format = "f")),
             hjust = -0.05, vjust = 1.1, size = 3.1, label.size = 0.25, fill = "white") +
    labs(title = paste0("D. ", comparison_name, ": weighted Pareto cutoff"),
         subtitle = "Equal weights: U(k) = Gnorm(k) - Rnorm(k)",
         x = "Opposite-arm Remainder crossings, R(k)",
         y = "Joint + permissible Disjoint, G(k)") +
    theme_manuscript()
}

make_cutoff_framework_figure <- function(comparison_name, control_group, treatment_group,
                                         group_results, knot_fit, scan_df, selected_k, c1, c2, out_file) {
  control <- group_results[[control_group]]$data
  treatment <- group_results[[treatment_group]]$data
  N <- nrow(control)
  cutoff_rank <- rank_cutoff_from_k(N, selected_k)
  long <- bind_rows(control %>% mutate(arm = control_group), treatment %>% mutate(arm = treatment_group))
  arm_colors <- c(setNames(COL$control, control_group), setNames(COL$treatment, treatment_group))

  pA <- add_rank_regions(ggplot(), c1, c2, N) +
    geom_vline(xintercept = c(c1, c2, cutoff_rank), linetype = c("dashed", "longdash", "dotdash"),
               color = c(COL$c1, COL$c2, COL$selected), linewidth = 0.8) +
    geom_line(data = long, aes(rank, abs_pc1_loading, color = arm), linewidth = 0.9) +
    scale_color_manual(values = arm_colors) +
    labs(title = paste0("A. ", comparison_name, ": absolute PC1-loading rank"),
         subtitle = "PCA input = log1p(DESeq2 median-of-ratios normalized counts)",
         x = "PC1 rank: low |loading| to high |loading|", y = "|PC1 loading|") +
    theme_manuscript()

  pB <- add_rank_regions(ggplot(), c1, c2, N) +
    geom_vline(xintercept = c(c1, c2, cutoff_rank), linetype = c("dashed", "longdash", "dotdash"),
               color = c(COL$c1, COL$c2, COL$selected), linewidth = 0.8) +
    geom_line(data = long, aes(rank, display_log1p_raw_empirical_variance, color = arm), linewidth = 0.9) +
    scale_color_manual(values = arm_colors) +
    labs(title = paste0("B. ", comparison_name, ": raw-count variance geometry"),
         subtitle = "Display only; does not determine c1, c2, or k*",
         x = "PC1 rank", y = "Smoothed log(1 + raw-count variance)") +
    theme_manuscript()

  fit_long <- bind_rows(lapply(seq_along(knot_fit$groups), function(j) {
    data.frame(rank = seq_len(N), arm = knot_fit$groups[j], fitted_D = knot_fit$fitted[, j])
  }))
  pC <- add_rank_regions(ggplot(), c1, c2, N) +
    geom_vline(xintercept = c(c1, c2, cutoff_rank), linetype = c("dashed", "longdash", "dotdash"),
               color = c(COL$c1, COL$c2, COL$selected), linewidth = 0.8) +
    geom_hline(yintercept = 0, color = "grey55", linetype = "dotted", linewidth = 0.35) +
    geom_line(data = long, aes(rank, cumulative_divergence, color = arm), linewidth = 0.65, alpha = 0.5) +
    geom_line(data = fit_long, aes(rank, fitted_D, color = arm), linewidth = 1.15) +
    scale_color_manual(values = arm_colors) +
    labs(title = paste0("C. ", comparison_name, ": PC1–variance cumulative divergence"),
         subtitle = "Shared two-knot fit within this comparison",
         x = "PC1 rank", y = "D(r) = F_E(r) - F_P(r)") +
    theme_manuscript()

  pD <- make_pareto_panel(scan_df, selected_k, comparison_name)
  save_grid_2x2(list(pA, pB, pC, pD), out_file)
}

make_nb_corroboration_figure <- function(comparison_name, selected_k, arm_results, out_file) {
  rank_data <- bind_rows(lapply(names(arm_results), function(a) arm_results[[a]]$rank_data %>% mutate(arm = a)))
  likelihood <- bind_rows(lapply(names(arm_results), function(a) arm_results[[a]]$likelihood))
  region_summary <- bind_rows(lapply(names(arm_results), function(a) arm_results[[a]]$region_summary))
  N <- max(rank_data$rank)
  cutoff_rank <- rank_cutoff_from_k(N, selected_k)

  pA <- ggplot(rank_data, aes(rank, excess_variance_signal_smooth, color = arm)) +
    geom_vline(xintercept = cutoff_rank, color = COL$selected, linetype = "dotdash", linewidth = 0.9) +
    geom_line(linewidth = 0.85) +
    labs(title = paste0("A. ", comparison_name, ": size-factor-aware excess variance"),
         subtitle = "Smoothed post-selection corroboration signal",
         x = "PC1 rank", y = "log(1 + excess variance)") +
    theme_manuscript()

  pB <- ggplot(rank_data, aes(rank, alpha_mu_signal_smooth, color = arm)) +
    geom_vline(xintercept = cutoff_rank, color = COL$selected, linetype = "dotdash", linewidth = 0.9) +
    geom_line(linewidth = 0.85) +
    labs(title = paste0("B. ", comparison_name, ": NB2 alpha*mu signal"),
         subtitle = "Moment estimate on the DESeq2-normalized scale",
         x = "PC1 rank", y = "log(1 + alpha_hat * mu)") +
    theme_manuscript()

  pC <- ggplot(region_summary, aes(corroboration_region, median_alpha_mu_signal, fill = corroboration_region)) +
    geom_col(width = 0.65) + facet_wrap(~ arm, nrow = 1) +
    scale_fill_manual(values = c("LEFT" = COL$left, "RIGHT" = COL$right), na.value = COL$right) +
    labs(title = paste0("C. ", comparison_name, ": matched-block corroboration"),
         subtitle = "LEFT is shown only when a full equal-sized preceding block exists",
         x = NULL, y = "Median alpha*mu signal") +
    theme_manuscript() + theme(legend.position = "none")

  likelihood$arm_region <- factor(paste(likelihood$arm, likelihood$region, sep = " / "),
                                  levels = rev(paste(likelihood$arm, likelihood$region, sep = " / ")))
  pD <- ggplot(likelihood, aes(two_delta_logLik_NB2_minus_NB1, arm_region, fill = preferred_model)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.55) +
    geom_col(width = 0.68) +
    scale_fill_manual(values = c(NB1 = COL$nb1, NB2 = COL$nb2, Tie = "grey60")) +
    labs(title = paste0("D. ", comparison_name, ": NB2 versus NB1 likelihood evidence"),
         subtitle = "Descriptive likelihood contrast; positive values favor NB2",
         x = "2 × (logLik_NB2 - logLik_NB1)", y = NULL) +
    theme_manuscript()

  save_grid_2x2(list(pA, pB, pC, pD), out_file)
}

make_summary_figure <- function(cutoff_summary, out_file) {
  df <- cutoff_summary %>% mutate(comparison = factor(comparison, levels = comparison))
  p <- ggplot(df, aes(comparison, selected_k)) +
    geom_col(fill = COL$selected, width = 0.65) +
    geom_text(aes(label = paste0("k*=", selected_k)), vjust = -0.35, fontface = "bold") +
    labs(title = "Comparison-specific empirical EVS cutoffs",
         subtitle = "DESeq2 median-of-ratios -> log1p -> PC1 -> weighted Pareto",
         x = NULL, y = "Selected top-k* per arm") +
    theme_manuscript() + theme(legend.position = "none") +
    expand_limits(y = max(df$selected_k) * 1.15)
  save_figure(p, out_file, width = 11.5, height = 7.2)
}

# ------------------------- METHODS / FIGURE LEGENDS --------------------------
write_methods_manuscript <- function(cutoff_summary) {
  cutoff_text <- paste(paste0(cutoff_summary$comparison, " k*=", cutoff_summary$selected_k), collapse = ", ")
  lines <- c(
    "# NormEVS empirical cutoff methods",
    "",
    "## Preprocessing and PC1 ranking",
    "",
    paste0(
      "Each RT/ZT comparison was analyzed independently. PASs with zero counts across all samples in the comparison were removed. ",
      "DESeq2 median-of-ratios size factors were estimated within the comparison. Normalized counts were calculated by dividing raw counts by the DESeq2 size factors. ",
      "For EVS ranking, log1p was applied to the normalized counts, and PCA was then performed independently within each arm with centering and without feature scaling. ",
      "PASs were ordered from lowest to highest absolute PC1 loading."
    ),
    "",
    "## PC1 and size-factor-aware variance mass",
    "",
    paste0(
      "For PAS i in arm g, PC1 contribution was P_ig=lambda_1g*loading_ig^2. ",
      "For the normalized count X_ij=Y_ij/s_j, the Poisson contribution to variance is mu_ig/s_j. ",
      "Accordingly, the arm-level Poisson variance baseline was mu_ig times the mean inverse size factor across samples in that arm. ",
      "Observed normalized-count sample variance was compared with this baseline, and nonnegative excess variance was defined as E_ig=max[V_ig - mu_ig*mean_j(1/s_j),0]. ",
      "P and E were normalized separately to probability masses and accumulated along the absolute-PC1-loading rank to obtain D_g(r)=F_E,g(r)-F_P,g(r)."
    ),
    "",
    "## Variance regimes and weighted Pareto cutoff",
    "",
    paste0(
      "Within each comparison, a shared two-knot continuous linear spline was fitted jointly to the two arm-specific D_g(r) curves. ",
      "The fitted knots defined Remainder (rank<c1), Divergence (c1<=rank<=c2), and Leading-edge (rank>c2) regimes. ",
      "Candidate top-k values were restricted to 1<=k<=N-c2. For each k, benefit G(k) equaled Joint PASs plus Disjoint PASs whose opposite-arm rank was in Leading-edge or Divergence. ",
      "Cost R(k) equaled Disjoint PASs whose opposite-arm rank was in Remainder. Pareto-optimal candidates maximized G while minimizing R. ",
      "On the Pareto frontier, G and R were min-max normalized and the equal-weight utility U(k)=G_norm(k)-R_norm(k) was maximized. ",
      "Ties were resolved by greater G, lower R, then larger k."
    ),
    "",
    "## Final EVS membership",
    "",
    paste0(
      "After k* was selected, the k* highest absolute-PC1-loading PASs were selected independently in both arms. ",
      "The final Leading Edge was the union of these two top-k* sets. Opposite-arm Remainder crossings contributed to the Pareto cost but did not override union membership. ",
      "All PASs outside the union were assigned to the Remainder."
    ),
    "",
    "## Post-selection NB corroboration",
    "",
    paste0(
      "NB corroboration was calculated only after k* was fixed. Moment-based excess variance and alpha estimates used the same size-factor-aware normalized-count variance baseline. ",
      "For likelihood corroboration, raw-count expectations were reconstructed as the arm-specific normalized mean multiplied by each sample's DESeq2 size factor. ",
      "Conditional on these plug-in means, NB1 and NB2 each fitted one dispersion parameter by maximum likelihood. NB1 used Var(Y)=mu+alpha*mu and size=mu/alpha; NB2 used Var(Y)=mu+alpha*mu^2 and size=1/alpha. ",
      "Support was summarized descriptively by log-likelihoods, 2*(logLik_NB2-logLik_NB1), and log10(L_NB2/L_NB1); no chi-square likelihood-ratio test was assumed because the two variance functions are not nested in that form."
    ),
    "",
    "## Selected empirical cutoffs",
    "",
    cutoff_text
  )
  writeLines(lines, file.path(OUT_ROOT, "Methods_Manuscript.md"))
}

write_figure_legends <- function() {
  lines <- c(
    "# Figure legends",
    "",
    "## Figure 1. Comparison-specific NormEVS cutoff framework",
    "Panel A shows absolute PC1-loading rank after PCA of log1p(DESeq2 median-of-ratios normalized counts). Panel B shows raw-count variance geometry for descriptive context. Panel C shows cumulative divergence between size-factor-aware excess-variance mass and PC1 variance-contribution mass with the fitted c1/c2 boundaries. Panel D shows all candidate cutoffs, the Pareto frontier, and the selected equal-weight optimum k*.",
    "",
    "## Figure 2. Post-selection NB corroboration",
    "Panels A-B show size-factor-aware excess-variance and alpha*mu signals along the PC1 rank. Panel C compares the selected RIGHT block with the immediately preceding equal-sized LEFT block when such a full block exists. Panel D shows descriptive NB2-versus-NB1 likelihood evidence after k* has been fixed."
  )
  writeLines(lines, file.path(OUT_ROOT, "Figure_Legends.md"))
}

# --------------------------------- RUN ---------------------------------------
input <- read_count_data(COUNT_FILE, GROUP_PATTERNS)
all_counts <- input$counts
annotation <- input$annotation

cutoff_rows <- list()
likelihood_rows <- list()
leading_rows <- list()
remainder_rows <- list()
expected_figures <- character(0)
expected_tables <- character(0)

for (comparison_name in names(COMPARISONS)) {
  message("============================================================")
  message("Analyzing comparison: ", comparison_name)

  mapping <- COMPARISONS[[comparison_name]]
  control_group <- unname(mapping[["control"]])
  treatment_group <- unname(mapping[["treatment"]])
  control_idx <- grep(GROUP_PATTERNS[[control_group]], colnames(all_counts))
  treatment_idx <- grep(GROUP_PATTERNS[[treatment_group]], colnames(all_counts))
  if (length(control_idx) < 2L || length(treatment_idx) < 2L)
    stop("Each arm requires at least two samples: ", comparison_name)

  sample_idx <- c(control_idx, treatment_idx)
  count_mat <- all_counts[, sample_idx, drop = FALSE]
  keep <- rowSums(count_mat) > 0
  count_mat <- count_mat[keep, , drop = FALSE]
  annotation_cmp <- annotation[match(rownames(count_mat), annotation$feature_id), , drop = FALSE]
  N <- nrow(count_mat)
  if (N < 20L) stop("Too few PASs after filtering: ", comparison_name)

  group_labels <- factor(c(rep(control_group, length(control_idx)), rep(treatment_group, length(treatment_idx))),
                         levels = c(control_group, treatment_group))
  names(group_labels) <- colnames(count_mat)

  deseq <- normalize_deseq2_comparison(count_mat, group_labels)
  normalized_counts <- deseq$normalized_counts
  size_factors <- deseq$size_factors

  group_results <- list()
  for (g in c(control_group, treatment_group)) {
    idx <- which(group_labels == g)
    group_results[[g]] <- compute_group_analysis(
      g,
      raw_counts_arm = count_mat[, idx, drop = FALSE],
      normalized_counts_arm = normalized_counts[, idx, drop = FALSE],
      size_factors_arm = size_factors[colnames(count_mat)[idx]]
    )
  }

  knot_fit <- fit_shared_knots(group_results)
  c1 <- knot_fit$c1
  c2 <- knot_fit$c2

  raw_scan <- scan_pair_cutoffs(group_results[[control_group]]$data,
                                group_results[[treatment_group]]$data,
                                c1, c2, comparison_name, control_group, treatment_group)
  opt <- select_weighted_pareto_optimum(raw_scan)
  scan_df <- opt$scan
  selected_k <- opt$selected_k
  cutoff_rank <- rank_cutoff_from_k(N, selected_k)

  class_df <- classify_pair_at_k(group_results[[control_group]]$data,
                                 group_results[[treatment_group]]$data,
                                 selected_k, c1, c2, comparison_name, control_group, treatment_group)

  lead_ids <- class_df$feature_id
  rem_ids <- setdiff(rownames(count_mat), lead_ids)
  lead_table <- class_df %>% left_join(annotation_cmp, by = "feature_id") %>%
    select(comparison, feature_id, gene_symbol, everything())

  rCmap <- make_rank_map(group_results[[control_group]]$data)
  rTmap <- make_rank_map(group_results[[treatment_group]]$data)
  remainder_table <- data.frame(
    comparison = comparison_name,
    feature_id = rem_ids,
    gene_symbol = annotation_cmp$gene_symbol[match(rem_ids, annotation_cmp$feature_id)],
    evs_membership = "Remainder",
    selected_k = selected_k,
    cutoff_rank = cutoff_rank,
    control_rank = as.integer(rCmap[rem_ids]),
    treatment_rank = as.integer(rTmap[rem_ids]),
    stringsAsFactors = FALSE
  )

  arm_results <- list()
  for (g in c(control_group, treatment_group)) {
    idx <- which(group_labels == g)
    arm_results[[g]] <- build_corroboration_for_arm(
      g,
      raw_counts_arm = count_mat[, idx, drop = FALSE],
      normalized_counts_arm = normalized_counts[, idx, drop = FALSE],
      size_factors_arm = size_factors[colnames(count_mat)[idx]],
      rank_df = group_results[[g]]$data,
      k = selected_k
    )
  }

  likelihood_table <- bind_rows(lapply(names(arm_results), function(g) arm_results[[g]]$likelihood)) %>%
    mutate(comparison = comparison_name, selected_k = selected_k, cutoff_rank = cutoff_rank, .before = 1)
  region_summary <- bind_rows(lapply(names(arm_results), function(g) arm_results[[g]]$region_summary)) %>%
    mutate(comparison = comparison_name, selected_k = selected_k, cutoff_rank = cutoff_rank, .before = 1)

  selected_scan <- scan_df %>% filter(k == selected_k) %>% slice(1L)
  cutoff_row <- data.frame(
    comparison = comparison_name,
    N = N,
    control_group = control_group,
    treatment_group = treatment_group,
    c1 = c1,
    c2 = c2,
    max_candidate_k = N - c2,
    selected_k = selected_k,
    cutoff_rank = cutoff_rank,
    selected_k_fraction_of_all_features = selected_k / N,
    selected_k_fraction_of_leading_edge_domain = selected_k / (N - c2),
    weighted_utility = selected_scan$weighted_utility,
    good_n = selected_scan$good_n,
    remainder_cross_n = selected_scan$remainder_cross_n,
    joint_n = selected_scan$joint_n,
    permissible_disjoint_n = selected_scan$permissible_disjoint_n,
    leading_edge_union_n = length(lead_ids),
    remainder_n = length(rem_ids),
    matched_left_available_control = arm_results[[control_group]]$matched_left_available,
    matched_left_available_treatment = arm_results[[treatment_group]]$matched_left_available,
    knot_fit_SSE = knot_fit$SSE,
    stringsAsFactors = FALSE
  )

  comp_fig_dir <- file.path(OUT_ROOT, comparison_name, "Figures")
  comp_tab_dir <- file.path(OUT_ROOT, comparison_name, "Tables")
  dir.create(comp_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(comp_tab_dir, recursive = TRUE, showWarnings = FALSE)

  cutoff_fig <- file.path(comp_fig_dir, paste0("Figure_", comparison_name, "_Cutoff_Framework.png"))
  nb_fig <- file.path(comp_fig_dir, paste0("Figure_", comparison_name, "_NB_Corroboration.png"))
  make_cutoff_framework_figure(comparison_name, control_group, treatment_group,
                               group_results, knot_fit, scan_df, selected_k, c1, c2, cutoff_fig)
  make_nb_corroboration_figure(comparison_name, selected_k, arm_results, nb_fig)

  table_paths <- c(
    cutoff = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_Cutoff_Summary.csv")),
    scan = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_Weighted_Pareto_Scan.csv")),
    lead = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_Leading_Edge_Sites.csv")),
    rem = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_Remainder_Sites.csv")),
    rank = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_PC1_Variance_Rank_Data.csv")),
    nb_regions = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_NB_Corroboration_Regions.csv")),
    nb_likelihood = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_NB1_NB2_Likelihood.csv")),
    sf = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_DESeq2_Size_Factors.csv"))
  )

  rank_table <- bind_rows(group_results[[control_group]]$data %>% mutate(arm = control_group),
                          group_results[[treatment_group]]$data %>% mutate(arm = treatment_group)) %>%
    left_join(annotation_cmp, by = "feature_id")

  write_csv(cutoff_row, table_paths[["cutoff"]])
  write_csv(scan_df, table_paths[["scan"]])
  write_csv(lead_table, table_paths[["lead"]])
  write_csv(remainder_table, table_paths[["rem"]])
  write_csv(rank_table, table_paths[["rank"]])
  write_csv(region_summary, table_paths[["nb_regions"]])
  write_csv(likelihood_table, table_paths[["nb_likelihood"]])
  write_csv(data.frame(comparison = comparison_name,
                       sample = names(size_factors),
                       group = as.character(group_labels[names(size_factors)]),
                       size_factor = as.numeric(size_factors)), table_paths[["sf"]])

  figs <- c(cutoff_fig, nb_fig)
  if (EXPORT_PDF) figs <- c(figs, sub("\\.png$", ".pdf", cutoff_fig), sub("\\.png$", ".pdf", nb_fig))
  verify_outputs(c(figs, unname(table_paths)))
  expected_figures <- c(expected_figures, figs)
  expected_tables <- c(expected_tables, unname(table_paths))

  cutoff_rows[[comparison_name]] <- cutoff_row
  likelihood_rows[[comparison_name]] <- likelihood_table
  leading_rows[[comparison_name]] <- lead_table
  remainder_rows[[comparison_name]] <- remainder_table

  message(comparison_name, ": k*=", selected_k, " | c1=", c1, " | c2=", c2,
          " | k*/(N-c2)=", signif(selected_k / (N - c2), 4),
          " | Lead union=", length(lead_ids), " | Remainder=", length(rem_ids),
          " | G=", selected_scan$good_n, " | R=", selected_scan$remainder_cross_n)
}

# --------------------------- CROSS-COMPARISON OUTPUTS ------------------------
cutoff_summary <- bind_rows(cutoff_rows)
likelihood_all <- bind_rows(likelihood_rows)
leading_all <- bind_rows(leading_rows)
remainder_all <- bind_rows(remainder_rows)

summary_paths <- c(
  cutoff = file.path(SUMMARY_TAB_DIR, "Table_Empirical_Cutoffs.csv"),
  likelihood = file.path(SUMMARY_TAB_DIR, "Table_NB1_NB2_Likelihood_All.csv"),
  lead = file.path(SUMMARY_TAB_DIR, "Table_Leading_Edge_Sites_All.csv"),
  rem = file.path(SUMMARY_TAB_DIR, "Table_Remainder_Sites_All.csv")
)
write_csv(cutoff_summary, summary_paths[["cutoff"]])
write_csv(likelihood_all, summary_paths[["likelihood"]])
write_csv(leading_all, summary_paths[["lead"]])
write_csv(remainder_all, summary_paths[["rem"]])
expected_tables <- c(expected_tables, unname(summary_paths))

summary_fig <- file.path(SUMMARY_FIG_DIR, "Figure_Empirical_Cutoff_Summary.png")
make_summary_figure(cutoff_summary, summary_fig)
expected_figures <- c(expected_figures, summary_fig)
if (EXPORT_PDF) expected_figures <- c(expected_figures, sub("\\.png$", ".pdf", summary_fig))

write_methods_manuscript(cutoff_summary)
write_figure_legends()

manifest_paths <- unique(c(expected_figures, expected_tables,
                           file.path(OUT_ROOT, "Methods_Manuscript.md"),
                           file.path(OUT_ROOT, "Figure_Legends.md")))
manifest <- data.frame(
  relative_path = sub(paste0("^", normalizePath(OUT_ROOT, winslash = "/", mustWork = TRUE), "/?"), "",
                      normalizePath(manifest_paths, winslash = "/", mustWork = TRUE)),
  type = ifelse(grepl("\\.(png|pdf)$", manifest_paths, ignore.case = TRUE), "figure",
                ifelse(grepl("\\.csv$", manifest_paths, ignore.case = TRUE), "table", "text")),
  size_bytes = file.info(manifest_paths)$size,
  stringsAsFactors = FALSE
)
manifest_path <- file.path(OUT_ROOT, "Table_Export_Manifest.csv")
write_csv(manifest, manifest_path)
expected_tables <- c(expected_tables, manifest_path)
verify_outputs(c(expected_figures, expected_tables,
                 file.path(OUT_ROOT, "Methods_Manuscript.md"),
                 file.path(OUT_ROOT, "Figure_Legends.md")))

# -------------------------------- ZIP FILES ----------------------------------
create_zip_from_files <- function(zip_path, files) {
  files <- unique(files[file.exists(files)])
  if (!length(files)) stop("No files available for zip: ", zip_path)
  root <- normalizePath(OUT_ROOT, winslash = "/", mustWork = TRUE)
  abs_files <- normalizePath(files, winslash = "/", mustWork = TRUE)
  rel <- sub(paste0("^", root, "/?"), "", abs_files)
  if (file.exists(zip_path)) unlink(zip_path)
  old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(OUT_ROOT)
  if (requireNamespace("zip", quietly = TRUE)) zip::zipr(zip_path, rel) else utils::zip(zip_path, rel, flags = "-r9X")
  if (!file.exists(zip_path) || file.info(zip_path)$size <= 0) stop("Zip creation failed: ", zip_path)
  invisible(zip_path)
}

fig_zip <- file.path(OUT_ROOT, "Figures_All.zip")
tab_zip <- file.path(OUT_ROOT, "Tables_All.zip")
all_zip <- file.path(OUT_ROOT, "NormEVS_All_Outputs.zip")
create_zip_from_files(fig_zip, expected_figures)
create_zip_from_files(tab_zip, expected_tables)
create_zip_from_files(all_zip, c(expected_figures, expected_tables,
                                 file.path(OUT_ROOT, "Methods_Manuscript.md"),
                                 file.path(OUT_ROOT, "Figure_Legends.md")))
verify_outputs(c(fig_zip, tab_zip, all_zip))

# ------------------------------ FINAL SUMMARY -------------------------------
message("============================================================")
message("NORMEVS FINAL ANALYSIS COMPLETE")
message("EVS PCA input: log1p(DESeq2 median-of-ratios normalized counts).")
message("Variance mass: size-factor-aware excess variance on normalized counts.")
message("Selected k*: ", paste(cutoff_summary$comparison, cutoff_summary$selected_k, sep = "=", collapse = "; "))
message("Methods: ", file.path(OUT_ROOT, "Methods_Manuscript.md"))
message("Figures ZIP: ", fig_zip)
message("Tables ZIP: ", tab_zip)
message("Complete ZIP: ", all_zip)
message("============================================================")
