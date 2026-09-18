#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# DUAL EVS EMPIRICAL CUTOFF ANALYSIS
# RawEVS + NormEVS
# =============================================================================
#
# PRIMARY OBJECTIVE
# -----------------
# Estimate one empirical eigenvector-splitting cutoff k* independently for each
# RT/ZT comparison under TWO EVS definitions:
#
#   1) RawEVS:
#      raw counts -> log1p(raw counts) -> arm-specific PCA -> PC1 ranking
#
#   2) NormEVS:
#      raw counts -> DESeq2 median-of-ratios normalization
#                 -> log1p(normalized counts)
#                 -> arm-specific PCA -> PC1 ranking
#
# The two methods are calibrated independently. A cutoff estimated under one
# representation is never imposed on the other.
#
# Comparisons:
#   RT0_ZT6, RT2_ZT8, RT4_ZT10, RT8_ZT14
#
# =============================================================================
# METHODS OVERVIEW
# =============================================================================
#
# A. FEATURE FILTERING
# --------------------
# Filtering is comparison-specific. A PAS is retained if its total raw count
# across the samples in that RT/ZT comparison is > 0.
#
# B. DESEQ2 NORMALIZATION
# -----------------------
# DESeq2 median-of-ratios size factors are estimated independently within each
# comparison. These normalized counts are used in NormEVS and in the
# size-factor-aware variance geometry for both methods.
#
# C. EVS RANKING
# --------------
# RawEVS PCA input:
#
#     X_raw = log(1 + raw count)
#
# NormEVS PCA input:
#
#     X_norm = log(1 + DESeq2-normalized count)
#
# PCA is performed independently within each arm with:
#
#     center = TRUE
#     scale. = FALSE
#
# PASs are ranked from LOW to HIGH absolute PC1 loading. Thus rank N is the
# strongest PC1-loading PAS.
#
# For PAS i in arm g:
#
#     P_ig = lambda_1g * loading_ig^2
#
# D. SIZE-FACTOR-AWARE EXCESS VARIANCE
# ------------------------------------
# Let:
#
#     Y_ij = raw count
#     s_j  = DESeq2 size factor
#     X_ij = Y_ij / s_j
#
# If raw counts are Poisson with E[Y_ij] = s_j * mu_i, then:
#
#     Var(X_ij) = mu_i / s_j
#
# Therefore the expected Poisson variance on the normalized-count scale is
# approximated within arm g as:
#
#     V_Pois,ig = mu_ig * mean_j(1 / s_j)
#
# The pooled within-group normalized empirical variance V_pool,i is calculated
# across the two groups in the comparison. For arm g:
#
#     E_ig = max(V_pool,i - V_Pois,ig, 0)
#
# This E term is used only for variance-mass geometry. PCA rankings differ
# between RawEVS and NormEVS.
#
# E. PC1-VARIANCE MASS DIVERGENCE
# -------------------------------
# Along each method-specific arm-specific PC1 rank:
#
#     p_g(r) = P_g(r) / sum_r P_g(r)
#     q_g(r) = E_g(r) / sum_r E_g(r)
#
#     F_P,g(r) = cumulative p_g(r)
#     F_E,g(r) = cumulative q_g(r)
#
#     D_g(r) = F_E,g(r) - F_P,g(r)
#
# Within each comparison and each EVS method, a shared two-knot continuous
# linear spline is fitted jointly to the two arm-specific D_g(r) curves.
#
# Regimes:
#
#     rank < c1          : Remainder
#     c1 <= rank <= c2   : Divergence
#     rank > c2          : Leading Edge
#
# F. WEIGHTED PARETO SELECTION
# ----------------------------
# Candidate top-k depths satisfy:
#
#     1 <= k <= N - c2
#
# For each k:
#
#     G(k) = Joint
#          + Disjoint with opposite arm in Leading Edge
#          + Disjoint with opposite arm in Divergence
#
#     R(k) = Disjoint with opposite arm in Remainder
#
# The Pareto frontier maximizes G and minimizes R.
#
# Frontier coordinates are min-max normalized:
#
#     U(k) = G_norm(k) - R_norm(k)
#
# Equal weights are used by default. k* maximizes U(k), with ties resolved by:
#
#     1. greater G
#     2. lower R
#     3. larger k
#
# G. FINAL EVS MEMBERSHIP
# -----------------------
# The final Leading Edge is the union of the two arm-specific top-k* sets.
# Opposite-arm Remainder crossings remain in the union; they are recorded as a
# Pareto cost rather than post-hoc deleted.
#
# H. POST-SELECTION NB1/NB2 CORROBORATION
# ---------------------------------------
# After k* is fixed, each arm's RIGHT block is its selected top-k* region. LEFT
# is the immediately preceding equal-sized block when available. If k* > N/2,
# matched-block corroboration is marked unavailable instead of changing k*.
#
# Raw-count moment summaries:
#
#     excess              = max(var - mean, 0)
#     NB2 excess signal   = log(1 + excess)
#     NB2-NB1 contrast    = log(1 + excess) - log(1 + mean)
#     alpha_hat           = max[(var - mean)/mean^2, 0]
#     alpha*mu signal     = log(1 + alpha_hat * mean)
#
# Likelihood corroboration uses raw counts and DESeq2 size factors. Expected
# raw means are reconstructed from arm-specific normalized means:
#
#     mu_ij = normalized_mean_i * s_j
#
# Conditional on those plug-in means:
#
#     NB1: Var(Y) = mu + alpha*mu
#          size   = mu/alpha
#
#     NB2: Var(Y) = mu + alpha*mu^2
#          size   = 1/alpha
#
# NB1 and NB2 each fit one dispersion parameter alpha by maximum likelihood.
# Positive 2*(logLik_NB2-logLik_NB1) and positive log10 likelihood ratio favor
# NB2. These quantities do not affect k*.
#
# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <- "/root/REAPER98632/exports/dual_raw_norm_evs_empirical_cutoff_final"

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

EVS_METHODS <- c("RawEVS", "NormEVS")

BENEFIT_WEIGHT       <- 1.0
CONTAMINATION_WEIGHT <- 1.0

KNOT_COARSE_GRID_POINTS <- 4096L

DISPLAY_VAR_SPAR <- 0.72
DISPLAY_D_SPAR   <- 0.72
DISPLAY_NB_SPAR  <- 0.68

PNG_DPI <- 360
EXPORT_PDF <- TRUE

NB_LOG_ALPHA_LOWER <- -14
NB_LOG_ALPHA_UPPER <- 8

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

for (method in EVS_METHODS) {
  dir.create(file.path(OUT_ROOT, method), recursive = TRUE, showWarnings = FALSE)
}

SUMMARY_FIG_DIR <- file.path(OUT_ROOT, "Summary", "Figures")
SUMMARY_TAB_DIR <- file.path(OUT_ROOT, "Summary", "Tables")
dir.create(SUMMARY_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# COLORS / THEME
# =============================================================================

COL <- list(
  control = "#386CB0",
  treatment = "#159D91",
  raw_evs = "#D95F02",
  norm_evs = "#1B9E77",
  raw = "#646464",
  divergence = "#6A3D9A",
  fit = "#111111",
  remainder = "#DCE6F2",
  interval = "#FFF0B3",
  leading = "#D8F3E7",
  c1 = "#D73027",
  c2 = "#1A9850",
  selected = "#B5179E",
  candidate_line = "#A7A7A7",
  pareto = "#5E3C99",
  nb2 = "#1B9E77",
  nbgap = "#CC1E8C",
  alphamu = "#386CB0",
  left = "#5B8FD1",
  right = "#43A047",
  nb1 = "#E69F00",
  nb2fit = "#0072B2"
)

theme_manuscript <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.4),
      plot.subtitle = element_text(size = base_size - 0.2, margin = margin(b = 5)),
      plot.caption = element_text(size = base_size - 1.6, color = "grey25", hjust = 0),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "#222222", size = base_size - 0.7),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_size - 0.4),
      legend.text = element_text(size = base_size - 0.8),
      legend.box = "vertical",
      panel.border = element_rect(color = "#B7B7B7", fill = NA, linewidth = 0.45),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 0.2),
      plot.margin = margin(8, 10, 8, 10)
    )
}

save_figure <- function(plot_obj, png_path, width = 15, height = 11) {
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = PNG_DPI,
    bg = "white",
    limitsize = FALSE
  )

  if (isTRUE(EXPORT_PDF)) {
    pdf_path <- sub("\\.png$", ".pdf", png_path, ignore.case = TRUE)
    ggplot2::ggsave(
      filename = pdf_path,
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      bg = "white",
      device = "pdf",
      limitsize = FALSE
    )
  }

  invisible(png_path)
}

save_grid_2x2 <- function(plots, png_path, width = 16, height = 12.5) {
  if (length(plots) != 4L) stop("save_grid_2x2 requires exactly four plots.")

  draw_once <- function(device_fun) {
    device_fun()
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(
          nrow = 2L,
          ncol = 2L,
          widths = unit(c(1, 1), "null"),
          heights = unit(c(1, 1), "null")
        )
      )
    )
    for (i in seq_along(plots)) {
      r <- if (i <= 2L) 1L else 2L
      c <- if (i %% 2L == 1L) 1L else 2L
      print(plots[[i]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))
    }
    grDevices::dev.off()
  }

  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)

  draw_once(function() {
    grDevices::png(
      filename = png_path,
      width = width,
      height = height,
      units = "in",
      res = PNG_DPI,
      bg = "white"
    )
  })

  if (isTRUE(EXPORT_PDF)) {
    pdf_path <- sub("\\.png$", ".pdf", png_path, ignore.case = TRUE)
    draw_once(function() {
      grDevices::pdf(
        file = pdf_path,
        width = width,
        height = height,
        onefile = TRUE,
        useDingbats = FALSE
      )
    })
  }

  invisible(png_path)
}

rank_cutoff_from_k <- function(N, k) {
  if (!is.finite(k) || k < 1L || k > N) return(NA_integer_)
  as.integer(N - k + 1L)
}

# =============================================================================
# INPUT
# =============================================================================

read_count_data <- function(path, group_patterns) {
  if (!file.exists(path)) stop("Count file does not exist: ", path)

  raw_df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  if (nrow(raw_df) < 1L || ncol(raw_df) < 2L) {
    stop("Count file is empty or malformed.")
  }

  sample_idx <- sort(unique(unlist(lapply(
    group_patterns,
    function(pattern) grep(pattern, colnames(raw_df))
  ))))

  if (!length(sample_idx)) stop("No sample columns matched GROUP_PATTERNS.")
  if (1L %in% sample_idx) stop("Column 1 must contain feature IDs, not samples.")

  feature_id <- trimws(as.character(raw_df[[1L]]))
  blank <- is.na(feature_id) | feature_id == ""
  if (any(blank)) feature_id[blank] <- paste0("__feature_row_", which(blank))
  feature_id <- make.unique(feature_id, sep = "__dup_")

  non_sample_idx <- setdiff(seq_len(ncol(raw_df)), sample_idx)
  symbol_candidates <- non_sample_idx[
    tolower(colnames(raw_df)[non_sample_idx]) %in%
      c("symbol", "gene_symbol", "genesymbol", "gene")
  ]

  if (length(symbol_candidates)) {
    gene_symbol <- as.character(raw_df[[symbol_candidates[1L]]])
  } else if (length(non_sample_idx) >= 2L) {
    gene_symbol <- as.character(raw_df[[non_sample_idx[2L]]])
  } else {
    gene_symbol <- rep(NA_character_, nrow(raw_df))
  }

  count_df <- raw_df[, sample_idx, drop = FALSE]
  count_mat <- do.call(cbind, lapply(
    count_df,
    function(x) suppressWarnings(as.numeric(trimws(as.character(x))))
  ))

  colnames(count_mat) <- colnames(count_df)
  rownames(count_mat) <- feature_id
  storage.mode(count_mat) <- "numeric"

  count_mat[!is.finite(count_mat)] <- 0
  count_mat <- pmax(count_mat, 0)

  annotation <- data.frame(
    feature_id = feature_id,
    gene_symbol = gene_symbol,
    stringsAsFactors = FALSE
  )

  list(counts = count_mat, annotation = annotation)
}

normalize_deseq2_comparison <- function(count_mat, group_labels) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required.")
  }

  col_data <- data.frame(
    group = factor(group_labels),
    row.names = colnames(count_mat)
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = col_data,
    design = ~ group
  )

  dds <- tryCatch(
    DESeq2::estimateSizeFactors(dds),
    error = function(e) {
      message("Default DESeq2 size factors failed; retrying with type='poscounts'.")
      DESeq2::estimateSizeFactors(dds, type = "poscounts")
    }
  )

  list(
    normalized_counts = DESeq2::counts(dds, normalized = TRUE),
    size_factors = DESeq2::sizeFactors(dds)
  )
}

compute_pooled_within_group_variance <- function(normalized_counts, group_labels) {
  groups <- levels(factor(group_labels))
  sse <- rep(0, nrow(normalized_counts))
  residual_df <- 0L

  for (g in groups) {
    idx <- which(group_labels == g)
    if (length(idx) < 2L) next

    xg <- normalized_counts[, idx, drop = FALSE]
    mu_g <- rowMeans(xg)
    resid_g <- sweep(xg, 1L, mu_g, "-")

    sse <- sse + rowSums(resid_g^2)
    residual_df <- residual_df + length(idx) - 1L
  }

  if (residual_df < 2L) stop("Pooled residual degrees of freedom < 2.")

  V <- sse / residual_df
  V[!is.finite(V)] <- 0
  V <- pmax(V, 0)
  names(V) <- rownames(normalized_counts)

  list(variance = V, residual_df = residual_df)
}

# =============================================================================
# PCA / METHOD-SPECIFIC RANKING
# =============================================================================

build_rank_matrix <- function(method, raw_counts_arm, normalized_counts_arm) {
  if (method == "RawEVS") {
    return(log1p(raw_counts_arm))
  }

  if (method == "NormEVS") {
    return(log1p(normalized_counts_arm))
  }

  stop("Unknown EVS method: ", method)
}

compute_pc1_rank <- function(rank_matrix_arm) {
  if (ncol(rank_matrix_arm) < 2L) {
    stop("At least two samples are required for PCA.")
  }

  pca <- stats::prcomp(
    t(rank_matrix_arm),
    center = TRUE,
    scale. = FALSE,
    rank. = 1
  )

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
  raw_var <- apply(raw_counts_arm, 1L, stats::var, na.rm = TRUE)
  raw_var[!is.finite(raw_var)] <- 0
  raw_var <- pmax(raw_var, 0)

  raw_var_ranked <- raw_var[rank_order]
  rank <- seq_along(rank_order)
  log_var <- log1p(raw_var_ranked)

  display_spline <- stats::smooth.spline(
    x = rank,
    y = log_var,
    spar = DISPLAY_VAR_SPAR
  )

  display_y <- as.numeric(
    stats::predict(display_spline, x = rank, deriv = 0)$y
  )

  data.frame(
    rank = rank,
    raw_empirical_variance = raw_var_ranked,
    log1p_raw_empirical_variance = log_var,
    display_log1p_raw_empirical_variance = display_y,
    stringsAsFactors = FALSE
  )
}

smooth_divergence_for_display <- function(rank, D, spar = DISPLAY_D_SPAR) {
  fit <- stats::smooth.spline(x = rank, y = D, spar = spar)

  y <- as.numeric(stats::predict(fit, x = rank, deriv = 0)$y)

  endpoint_line <- seq(
    y[1L],
    y[length(y)],
    length.out = length(y)
  )

  y - endpoint_line
}

compute_group_analysis <- function(
    method,
    group_name,
    raw_counts_arm,
    normalized_counts_arm,
    size_factors_arm,
    pooled_variance) {

  rank_matrix <- build_rank_matrix(
    method = method,
    raw_counts_arm = raw_counts_arm,
    normalized_counts_arm = normalized_counts_arm
  )

  pc1 <- compute_pc1_rank(rank_matrix)
  rank_order <- pc1$rank_order

  geometry <- compute_raw_variance_geometry(
    raw_counts_arm = raw_counts_arm,
    rank_order = rank_order
  )

  mu_norm <- rowMeans(normalized_counts_arm, na.rm = TRUE)
  mu_norm[!is.finite(mu_norm)] <- 0
  mu_norm <- pmax(mu_norm, 0)

  sf <- as.numeric(size_factors_arm)
  if (any(!is.finite(sf)) || any(sf <= 0)) {
    stop("Invalid DESeq2 size factor in ", group_name)
  }

  poisson_scale_factor <- mean(1 / sf)

  P_ranked <- pc1$pc1_variance_contribution[rank_order]
  mu_ranked <- mu_norm[rank_order]
  V_pool_ranked <- pooled_variance[rank_order]

  poisson_variance_ranked <- mu_ranked * poisson_scale_factor
  E_ranked <- pmax(V_pool_ranked - poisson_variance_ranked, 0)

  P_total <- sum(P_ranked)
  E_total <- sum(E_ranked)

  if (!is.finite(P_total) || P_total <= 0) {
    stop("PC1 variance mass undefined for ", method, " / ", group_name)
  }

  if (!is.finite(E_total) || E_total <= 0) {
    stop("Excess-variance mass undefined for ", method, " / ", group_name)
  }

  p_mass <- P_ranked / P_total
  q_mass <- E_ranked / E_total

  F_P <- cumsum(p_mass)
  F_E <- cumsum(q_mass)
  D <- F_E - F_P
  rank <- seq_along(rank_order)

  display_D <- smooth_divergence_for_display(rank, D)

  df <- geometry %>%
    mutate(
      evs_method = method,
      group = group_name,
      feature_id = rownames(raw_counts_arm)[rank_order],
      pc1_loading = pc1$loading[rank_order],
      abs_pc1_loading = pc1$abs_loading[rank_order],
      pc1_eigenvalue = pc1$lambda1,
      pc1_variance_contribution = P_ranked,
      pc1_variance_mass = p_mass,
      normalized_group_mean = mu_ranked,
      pooled_normalized_variance = V_pool_ranked,
      poisson_reference_variance = poisson_variance_ranked,
      nb_excess_variance = E_ranked,
      nb_excess_variance_mass = q_mass,
      cumulative_pc1_mass = F_P,
      cumulative_nb_mass = F_E,
      cumulative_divergence = D,
      display_D = display_D
    )

  list(
    data = df,
    rank_order = rank_order,
    poisson_scale_factor = poisson_scale_factor
  )
}

# =============================================================================
# SHARED TWO-KNOT FIT
# =============================================================================

piecewise_basis <- function(x, c1, c2) {
  cbind(
    intercept = 1,
    x = x,
    hinge1 = pmax(x - c1, 0),
    hinge2 = pmax(x - c2, 0)
  )
}

piecewise_sse <- function(par, x, D_mat, min_gap) {
  c1 <- par[1L]
  c2 <- par[2L]

  if (
    !is.finite(c1) ||
    !is.finite(c2) ||
    c1 <= 0 ||
    c2 >= 1 ||
    c2 - c1 <= min_gap
  ) return(1e100)

  X <- piecewise_basis(x, c1, c2)

  coef <- tryCatch(
    qr.coef(qr(X), D_mat),
    error = function(e) NULL
  )

  if (is.null(coef) || any(!is.finite(coef))) return(1e100)

  resid <- D_mat - X %*% coef
  sum(resid^2)
}

fit_shared_knots <- function(group_results) {
  groups <- names(group_results)

  N_values <- vapply(
    group_results,
    function(z) nrow(z$data),
    integer(1)
  )

  if (length(unique(N_values)) != 1L) {
    stop("All groups must contain the same number of ranked features.")
  }

  N <- N_values[1L]
  if (N < 10L) stop("Too few ranked features for two-knot fitting.")

  x_full <- (seq_len(N) - 1) / (N - 1)

  D_full <- do.call(
    cbind,
    lapply(group_results, function(z) z$data$cumulative_divergence)
  )
  colnames(D_full) <- groups

  opt_n <- min(KNOT_COARSE_GRID_POINTS, N)
  opt_idx <- unique(as.integer(round(seq(1, N, length.out = opt_n))))

  x_opt <- x_full[opt_idx]
  D_opt <- D_full[opt_idx, , drop = FALSE]

  min_gap <- max(4 / (N - 1), .Machine$double.eps^0.25)

  starts <- list(
    c(0.03, 0.97),
    c(0.08, 0.92),
    c(0.15, 0.85),
    c(0.25, 0.75),
    c(0.35, 0.65)
  )

  coarse <- lapply(starts, function(start) {
    stats::optim(
      par = start,
      fn = piecewise_sse,
      x = x_opt,
      D_mat = D_opt,
      min_gap = min_gap,
      method = "Nelder-Mead",
      control = list(maxit = 700, reltol = 1e-11)
    )
  })

  values <- vapply(coarse, function(z) z$value, numeric(1))
  best <- coarse[[which.min(values)]]

  refined <- stats::optim(
    par = best$par,
    fn = piecewise_sse,
    x = x_full,
    D_mat = D_full,
    min_gap = min_gap,
    method = "Nelder-Mead",
    control = list(maxit = 1000, reltol = 1e-12)
  )

  if (!is.finite(refined$value)) stop("Shared-knot optimization failed.")

  c1_rank <- as.integer(round(1 + refined$par[1L] * (N - 1)))
  c2_rank <- as.integer(round(1 + refined$par[2L] * (N - 1)))

  c1_rank <- max(2L, min(N - 2L, c1_rank))
  c2_rank <- max(c1_rank + 1L, min(N - 1L, c2_rank))

  c1_x <- (c1_rank - 1) / (N - 1)
  c2_x <- (c2_rank - 1) / (N - 1)

  X <- piecewise_basis(x_full, c1_x, c2_x)
  coef <- qr.coef(qr(X), D_full)
  fitted <- X %*% coef

  list(
    c1 = c1_rank,
    c2 = c2_rank,
    x = x_full,
    fitted = fitted,
    SSE = sum((D_full - fitted)^2),
    groups = groups
  )
}

# =============================================================================
# PARETO SCAN
# =============================================================================

make_rank_map <- function(df) {
  stats::setNames(df$rank, df$feature_id)
}

rank_to_region <- function(rank, c1, c2) {
  ifelse(
    rank < c1,
    "Remainder",
    ifelse(rank <= c2, "Divergence", "LeadingEdge")
  )
}

cumulative_activation <- function(depth, K) {
  depth <- as.integer(depth)
  keep <- is.finite(depth) & depth >= 1L & depth <= K

  if (!any(keep)) return(rep(0L, K))

  cumsum(tabulate(depth[keep], nbins = K))
}

active_interval_count <- function(starts, ends, K) {
  if (length(starts) == 0L) return(rep(0L, K))

  starts <- as.integer(starts)
  ends <- as.integer(ends)

  valid <- (
    is.finite(starts) &
    is.finite(ends) &
    starts >= 1L &
    starts <= K &
    ends > starts
  )

  starts <- starts[valid]
  ends <- pmin(ends[valid], K + 1L)

  if (!length(starts)) return(rep(0L, K))

  start_tab <- tabulate(starts, nbins = K + 1L)
  end_tab <- tabulate(ends, nbins = K + 1L)

  cumsum(start_tab - end_tab)[seq_len(K)]
}

scan_pair_cutoffs <- function(
    control_df,
    treatment_df,
    c1,
    c2,
    comparison_name,
    control_group,
    treatment_group,
    method) {

  if (nrow(control_df) != nrow(treatment_df)) {
    stop("Control and treatment rankings have different feature counts.")
  }

  if (!setequal(control_df$feature_id, treatment_df$feature_id)) {
    stop("Control and treatment rankings do not contain the same feature IDs.")
  }

  N <- nrow(control_df)
  K <- as.integer(N - c2)

  if (K < 1L) stop("No candidate top-k depth exists beyond c2.")

  rank_control <- make_rank_map(control_df)
  rank_treatment <- make_rank_map(treatment_df)

  ids <- control_df$feature_id
  rC <- as.integer(unname(rank_control[ids]))
  rT <- as.integer(unname(rank_treatment[ids]))

  if (any(!is.finite(rC)) || any(!is.finite(rT))) {
    stop("Non-finite rank encountered in ", comparison_name, " / ", method)
  }

  dC <- N - rC + 1L
  dT <- N - rT + 1L
  k <- seq_len(K)

  joint_depth <- pmax(dC, dT)
  joint_n <- cumulative_activation(joint_depth, K)

  idx_le_C <- which(dC < dT & dC <= K & dT <= K)
  disjoint_control_opposite_le_n <- active_interval_count(
    dC[idx_le_C], dT[idx_le_C], K
  )

  idx_le_T <- which(dT < dC & dT <= K & dC <= K)
  disjoint_treatment_opposite_le_n <- active_interval_count(
    dT[idx_le_T], dC[idx_le_T], K
  )

  idx_div_C <- which(dC <= K & rT >= c1 & rT <= c2)
  disjoint_control_opposite_divergence_n <- cumulative_activation(
    dC[idx_div_C], K
  )

  idx_div_T <- which(dT <= K & rC >= c1 & rC <= c2)
  disjoint_treatment_opposite_divergence_n <- cumulative_activation(
    dT[idx_div_T], K
  )

  idx_rem_C <- which(dC <= K & rT < c1)
  remainder_cross_control_n <- cumulative_activation(dC[idx_rem_C], K)

  idx_rem_T <- which(dT <= K & rC < c1)
  remainder_cross_treatment_n <- cumulative_activation(dT[idx_rem_T], K)

  disjoint_opposite_le_n <- (
    disjoint_control_opposite_le_n +
    disjoint_treatment_opposite_le_n
  )

  disjoint_opposite_divergence_n <- (
    disjoint_control_opposite_divergence_n +
    disjoint_treatment_opposite_divergence_n
  )

  permissible_disjoint_control_n <- (
    disjoint_control_opposite_le_n +
    disjoint_control_opposite_divergence_n
  )

  permissible_disjoint_treatment_n <- (
    disjoint_treatment_opposite_le_n +
    disjoint_treatment_opposite_divergence_n
  )

  permissible_disjoint_n <- (
    permissible_disjoint_control_n +
    permissible_disjoint_treatment_n
  )

  good_n <- joint_n + permissible_disjoint_n

  remainder_cross_n <- (
    remainder_cross_control_n +
    remainder_cross_treatment_n
  )

  union_n <- good_n + remainder_cross_n

  union_depth <- pmin(dC, dT)
  union_check <- cumulative_activation(union_depth, K)

  if (!all(union_n == union_check)) {
    stop("Internal union-count mismatch in ", comparison_name, " / ", method)
  }

  data.frame(
    evs_method = method,
    comparison = comparison_name,
    control_group = control_group,
    treatment_group = treatment_group,
    k = k,
    cutoff_rank = N - k + 1L,
    joint_n = joint_n,
    disjoint_control_opposite_le_n = disjoint_control_opposite_le_n,
    disjoint_treatment_opposite_le_n = disjoint_treatment_opposite_le_n,
    disjoint_opposite_le_n = disjoint_opposite_le_n,
    disjoint_control_opposite_divergence_n =
      disjoint_control_opposite_divergence_n,
    disjoint_treatment_opposite_divergence_n =
      disjoint_treatment_opposite_divergence_n,
    disjoint_opposite_divergence_n = disjoint_opposite_divergence_n,
    permissible_disjoint_control_n = permissible_disjoint_control_n,
    permissible_disjoint_treatment_n = permissible_disjoint_treatment_n,
    permissible_disjoint_n = permissible_disjoint_n,
    good_n = good_n,
    remainder_cross_control_n = remainder_cross_control_n,
    remainder_cross_treatment_n = remainder_cross_treatment_n,
    remainder_cross_n = remainder_cross_n,
    union_n = union_n,
    retained_fraction = ifelse(union_n > 0, good_n / union_n, NA_real_),
    remainder_cross_fraction = ifelse(
      union_n > 0, remainder_cross_n / union_n, NA_real_
    ),
    divergence_disjoint_fraction = ifelse(
      union_n > 0,
      disjoint_opposite_divergence_n / union_n,
      NA_real_
    ),
    jaccard_top_k = ifelse(union_n > 0, joint_n / union_n, NA_real_),
    stringsAsFactors = FALSE
  )
}

mark_pareto_frontier <- function(
    scan_df,
    good_col = "good_n",
    cost_col = "remainder_cross_n") {

  if (nrow(scan_df) < 1L) stop("Empty cutoff scan.")

  tmp <- scan_df %>%
    transmute(
      row_id = row_number(),
      k = k,
      good = .data[[good_col]],
      cost = .data[[cost_col]]
    ) %>%
    arrange(cost, desc(good), desc(k)) %>%
    group_by(cost) %>%
    slice(1L) %>%
    ungroup() %>%
    arrange(cost, desc(good))

  running_best_before <- c(-Inf, head(cummax(tmp$good), -1L))
  tmp$is_frontier_coord <- tmp$good > running_best_before

  frontier <- tmp %>%
    filter(is_frontier_coord) %>%
    arrange(cost, good, k)

  key_all <- paste(scan_df[[cost_col]], scan_df[[good_col]], sep = "::")
  key_frontier <- paste(frontier$cost, frontier$good, sep = "::")

  out <- scan_df
  out$is_pareto <- key_all %in% key_frontier

  list(scan = out, frontier = frontier)
}

select_weighted_pareto_optimum <- function(
    scan_df,
    good_col = "good_n",
    cost_col = "remainder_cross_n",
    benefit_weight = BENEFIT_WEIGHT,
    contamination_weight = CONTAMINATION_WEIGHT) {

  if (
    !is.finite(benefit_weight) ||
    !is.finite(contamination_weight) ||
    benefit_weight < 0 ||
    contamination_weight < 0 ||
    (benefit_weight + contamination_weight) <= 0
  ) {
    stop("Pareto weights must be finite, non-negative, and not both zero.")
  }

  marked <- mark_pareto_frontier(scan_df, good_col, cost_col)
  frontier <- marked$frontier %>% arrange(cost, good, k)

  if (nrow(frontier) < 1L) stop("No Pareto-optimal cutoff points identified.")

  good_range <- range(frontier$good, na.rm = TRUE)
  cost_range <- range(frontier$cost, na.rm = TRUE)

  normalize_good <- function(x) {
    if (diff(good_range) == 0) rep(1, length(x))
    else (x - good_range[1L]) / diff(good_range)
  }

  normalize_cost <- function(x) {
    if (diff(cost_range) == 0) rep(0, length(x))
    else (x - cost_range[1L]) / diff(cost_range)
  }

  frontier$good_norm <- normalize_good(frontier$good)
  frontier$remainder_norm <- normalize_cost(frontier$cost)

  frontier$weighted_utility <- (
    benefit_weight * frontier$good_norm -
    contamination_weight * frontier$remainder_norm
  )

  best_utility <- max(frontier$weighted_utility, na.rm = TRUE)

  chosen <- frontier %>%
    filter(abs(weighted_utility - best_utility) < 1e-12) %>%
    arrange(desc(good), cost, desc(k)) %>%
    slice(1L)

  selected_k <- as.integer(chosen$k[1L])

  out <- marked$scan
  out$good_norm <- normalize_good(out[[good_col]])
  out$remainder_norm <- normalize_cost(out[[cost_col]])
  out$weighted_utility <- (
    benefit_weight * out$good_norm -
    contamination_weight * out$remainder_norm
  )
  out$is_selected_weighted <- out$k == selected_k

  zero_idx <- which(out[[cost_col]] == 0)
  zero_max_k <- if (length(zero_idx)) max(out$k[zero_idx]) else NA_integer_

  list(
    scan = out,
    frontier = frontier,
    selected_k = selected_k,
    selected_good = chosen$good[1L],
    selected_cost = chosen$cost[1L],
    selected_good_norm = chosen$good_norm[1L],
    selected_remainder_norm = chosen$remainder_norm[1L],
    selected_utility = chosen$weighted_utility[1L],
    benefit_weight = benefit_weight,
    contamination_weight = contamination_weight,
    selection_method = "weighted_pareto_utility",
    max_zero_remainder_crossing_k = zero_max_k
  )
}

# =============================================================================
# FINAL CLASSIFICATION
# =============================================================================

classify_pair_at_k <- function(
    control_df,
    treatment_df,
    k,
    c1,
    c2,
    comparison_name,
    control_group,
    treatment_group,
    method) {

  if (nrow(control_df) != nrow(treatment_df)) {
    stop("Control and treatment rankings have different feature counts.")
  }

  if (!setequal(control_df$feature_id, treatment_df$feature_id)) {
    stop("Control and treatment rankings do not contain the same feature IDs.")
  }

  N <- nrow(control_df)
  Kmax <- N - c2

  if (k < 1L || k > Kmax) {
    stop(
      "k must satisfy 1 <= k <= N-c2. Received k=",
      k, "; N-c2=", Kmax, "."
    )
  }

  control_top <- tail(control_df$feature_id, k)
  treatment_top <- tail(treatment_df$feature_id, k)
  union_ids <- union(control_top, treatment_top)

  rank_control <- make_rank_map(control_df)
  rank_treatment <- make_rank_map(treatment_df)

  rC <- as.integer(unname(rank_control[union_ids]))
  rT <- as.integer(unname(rank_treatment[union_ids]))

  in_control <- union_ids %in% control_top
  in_treatment <- union_ids %in% treatment_top

  joint <- in_control & in_treatment
  control_only <- in_control & !in_treatment
  treatment_only <- in_treatment & !in_control

  control_region <- rank_to_region(rC, c1, c2)
  treatment_region <- rank_to_region(rT, c1, c2)

  base_class <- ifelse(
    joint,
    "Joint",
    ifelse(
      control_only,
      paste0("Disjoint_", control_group),
      paste0("Disjoint_", treatment_group)
    )
  )

  opposite_region <- rep(NA_character_, length(union_ids))
  opposite_region[control_only] <- treatment_region[control_only]
  opposite_region[treatment_only] <- control_region[treatment_only]

  disjoint_opposite_leading_edge <- (
    (control_only & rT > c2) |
    (treatment_only & rC > c2)
  )

  disjoint_opposite_divergence <- (
    (control_only & rT >= c1 & rT <= c2) |
    (treatment_only & rC >= c1 & rC <= c2)
  )

  cross_into_remainder <- (
    (control_only & rT < c1) |
    (treatment_only & rC < c1)
  )

  analysis_class <- rep(NA_character_, length(union_ids))
  analysis_class[joint] <- "Joint"

  analysis_class[control_only & rT > c2] <- paste0(
    "Disjoint_", control_group, "_OppositeLeadingEdge"
  )
  analysis_class[treatment_only & rC > c2] <- paste0(
    "Disjoint_", treatment_group, "_OppositeLeadingEdge"
  )
  analysis_class[control_only & rT >= c1 & rT <= c2] <- paste0(
    "Disjoint_", control_group, "_OppositeDivergence"
  )
  analysis_class[treatment_only & rC >= c1 & rC <= c2] <- paste0(
    "Disjoint_", treatment_group, "_OppositeDivergence"
  )
  analysis_class[control_only & rT < c1] <- paste0(
    "Disjoint_", control_group, "_OppositeRemainder_Cost"
  )
  analysis_class[treatment_only & rC < c1] <- paste0(
    "Disjoint_", treatment_group, "_OppositeRemainder_Cost"
  )

  out <- data.frame(
    evs_method = method,
    comparison = comparison_name,
    feature_id = union_ids,
    selected_k = as.integer(k),
    cutoff_rank = rank_cutoff_from_k(N, k),
    control_group = control_group,
    treatment_group = treatment_group,
    base_class = base_class,
    analysis_class = analysis_class,
    control_rank = rC,
    treatment_rank = rT,
    control_region = control_region,
    treatment_region = treatment_region,
    opposite_region = opposite_region,
    control_top_k = in_control,
    treatment_top_k = in_treatment,
    disjoint_opposite_leading_edge = disjoint_opposite_leading_edge,
    disjoint_opposite_divergence = disjoint_opposite_divergence,
    cross_into_remainder = cross_into_remainder,
    retained_for_analysis = TRUE,
    stringsAsFactors = FALSE
  )

  if (any(is.na(out$analysis_class))) {
    stop("Unclassified union site encountered in ", comparison_name, " / ", method)
  }

  out
}

summarize_classification <- function(class_df) {
  if (nrow(class_df) < 1L) stop("Cannot summarize empty classification.")

  joint_n <- sum(class_df$base_class == "Joint")

  control_only <- class_df$control_top_k & !class_df$treatment_top_k
  treatment_only <- class_df$treatment_top_k & !class_df$control_top_k

  disjoint_control_opposite_le_n <- sum(
    control_only & class_df$disjoint_opposite_leading_edge
  )
  disjoint_treatment_opposite_le_n <- sum(
    treatment_only & class_df$disjoint_opposite_leading_edge
  )

  disjoint_control_opposite_divergence_n <- sum(
    control_only & class_df$disjoint_opposite_divergence
  )
  disjoint_treatment_opposite_divergence_n <- sum(
    treatment_only & class_df$disjoint_opposite_divergence
  )

  disjoint_opposite_le_n <- (
    disjoint_control_opposite_le_n +
    disjoint_treatment_opposite_le_n
  )

  disjoint_opposite_divergence_n <- (
    disjoint_control_opposite_divergence_n +
    disjoint_treatment_opposite_divergence_n
  )

  permissible_disjoint_control_n <- (
    disjoint_control_opposite_le_n +
    disjoint_control_opposite_divergence_n
  )

  permissible_disjoint_treatment_n <- (
    disjoint_treatment_opposite_le_n +
    disjoint_treatment_opposite_divergence_n
  )

  permissible_disjoint_n <- (
    permissible_disjoint_control_n +
    permissible_disjoint_treatment_n
  )

  good_n <- joint_n + permissible_disjoint_n
  remainder_cross_n <- sum(class_df$cross_into_remainder)
  union_n <- nrow(class_df)

  if (good_n + remainder_cross_n != union_n) {
    stop("Classification counts do not sum to union size.")
  }

  data.frame(
    k = unique(class_df$selected_k)[1L],
    cutoff_rank = unique(class_df$cutoff_rank)[1L],
    joint_n = joint_n,
    disjoint_opposite_le_n = disjoint_opposite_le_n,
    disjoint_opposite_divergence_n = disjoint_opposite_divergence_n,
    permissible_disjoint_n = permissible_disjoint_n,
    good_n = good_n,
    remainder_cross_n = remainder_cross_n,
    union_n = union_n,
    retained_fraction = good_n / union_n,
    remainder_cross_fraction = remainder_cross_n / union_n,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# NB1 / NB2 CORROBORATION
# =============================================================================

compute_ranked_nb_metrics <- function(
    raw_counts_arm,
    normalized_counts_arm,
    ranked_feature_ids) {

  raw_x <- raw_counts_arm[ranked_feature_ids, , drop = FALSE]
  norm_x <- normalized_counts_arm[ranked_feature_ids, , drop = FALSE]

  raw_mu <- rowMeans(raw_x, na.rm = TRUE)
  raw_var <- apply(raw_x, 1L, stats::var, na.rm = TRUE)
  norm_mu <- rowMeans(norm_x, na.rm = TRUE)
  norm_var <- apply(norm_x, 1L, stats::var, na.rm = TRUE)

  raw_mu[!is.finite(raw_mu)] <- 0
  raw_var[!is.finite(raw_var)] <- 0
  norm_mu[!is.finite(norm_mu)] <- 0
  norm_var[!is.finite(norm_var)] <- 0

  raw_mu <- pmax(raw_mu, 0)
  raw_var <- pmax(raw_var, 0)
  norm_mu <- pmax(norm_mu, 0)
  norm_var <- pmax(norm_var, 0)

  excess <- pmax(raw_var - raw_mu, 0)

  alpha_hat <- rep(0, length(raw_mu))
  pos <- raw_mu > 0
  alpha_hat[pos] <- pmax(excess[pos] / (raw_mu[pos]^2), 0)

  data.frame(
    rank = seq_along(ranked_feature_ids),
    feature_id = ranked_feature_ids,
    raw_mean = raw_mu,
    raw_variance = raw_var,
    normalized_mean = norm_mu,
    normalized_variance = norm_var,
    nb2_excess_signal = log1p(excess),
    nb2_nb1_contrast = log1p(excess) - log1p(raw_mu),
    alpha_hat = alpha_hat,
    alpha_mu_signal = log1p(alpha_hat * raw_mu),
    stringsAsFactors = FALSE
  )
}

assign_matched_regions <- function(rank_df, k) {
  N <- nrow(rank_df)
  cutoff_rank <- rank_cutoff_from_k(N, k)

  right_start <- cutoff_rank
  right_end <- N
  left_end <- right_start - 1L
  left_start <- left_end - k + 1L

  rank_df$corroboration_region <- "Other"

  if (left_start < 1L) {
    attr(rank_df, "matched_block_available") <- FALSE
    attr(rank_df, "left_start") <- NA_integer_
    attr(rank_df, "left_end") <- NA_integer_
    attr(rank_df, "right_start") <- right_start
    attr(rank_df, "right_end") <- right_end

    rank_df$corroboration_region[
      rank_df$rank >= right_start & rank_df$rank <= right_end
    ] <- "RIGHT"

    return(rank_df)
  }

  rank_df$corroboration_region[
    rank_df$rank >= left_start & rank_df$rank <= left_end
  ] <- "LEFT"

  rank_df$corroboration_region[
    rank_df$rank >= right_start & rank_df$rank <= right_end
  ] <- "RIGHT"

  attr(rank_df, "matched_block_available") <- TRUE
  attr(rank_df, "left_start") <- left_start
  attr(rank_df, "left_end") <- left_end
  attr(rank_df, "right_start") <- right_start
  attr(rank_df, "right_end") <- right_end

  rank_df
}

smooth_metric_for_display <- function(rank, value, spar = DISPLAY_NB_SPAR) {
  ok <- is.finite(rank) & is.finite(value)
  out <- rep(NA_real_, length(value))
  if (sum(ok) < 8L) return(out)

  fit <- tryCatch(
    stats::smooth.spline(rank[ok], value[ok], spar = spar),
    error = function(e) NULL
  )

  if (is.null(fit)) return(out)

  out[ok] <- as.numeric(
    stats::predict(fit, x = rank[ok], deriv = 0)$y
  )

  out
}

fit_nb_dispersion_model <- function(y, mu, model = c("NB1", "NB2")) {
  model <- match.arg(model)

  y <- as.numeric(y)
  mu <- as.numeric(mu)

  ok <- is.finite(y) & is.finite(mu) & y >= 0 & mu > 0
  y <- y[ok]
  mu <- pmax(mu[ok], 1e-10)

  if (length(y) < 10L) {
    return(list(
      model = model,
      alpha = NA_real_,
      logLik = NA_real_,
      n_obs = length(y)
    ))
  }

  neg_loglik <- function(log_alpha) {
    alpha <- exp(log_alpha)

    size <- if (model == "NB2") {
      rep(1 / alpha, length(mu))
    } else {
      mu / alpha
    }

    size <- pmax(size, 1e-10)

    ll <- suppressWarnings(
      stats::dnbinom(
        x = round(y),
        mu = mu,
        size = size,
        log = TRUE
      )
    )

    if (any(!is.finite(ll))) return(1e100)
    -sum(ll)
  }

  opt <- stats::optimize(
    f = neg_loglik,
    interval = c(NB_LOG_ALPHA_LOWER, NB_LOG_ALPHA_UPPER)
  )

  list(
    model = model,
    alpha = exp(opt$minimum),
    logLik = -opt$objective,
    n_obs = length(y)
  )
}

fit_nb1_nb2_region <- function(
    raw_counts_arm,
    normalized_counts_arm,
    size_factors,
    feature_ids,
    arm_label,
    region_label) {

  if (!length(feature_ids)) {
    return(data.frame(
      arm = arm_label,
      region = region_label,
      n_features = 0L,
      n_observations = 0L,
      alpha_NB1 = NA_real_,
      alpha_NB2 = NA_real_,
      logLik_NB1 = NA_real_,
      logLik_NB2 = NA_real_,
      two_delta_logLik_NB2_minus_NB1 = NA_real_,
      log10_likelihood_ratio_NB2_to_NB1 = NA_real_,
      preferred_model = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  raw_sub <- raw_counts_arm[feature_ids, , drop = FALSE]
  norm_sub <- normalized_counts_arm[feature_ids, , drop = FALSE]

  feature_mu_norm <- rowMeans(norm_sub, na.rm = TRUE)
  feature_mu_norm[!is.finite(feature_mu_norm)] <- 0
  feature_mu_norm <- pmax(feature_mu_norm, 1e-10)

  sf <- as.numeric(size_factors[colnames(raw_sub)])
  if (any(!is.finite(sf)) || any(sf <= 0)) {
    stop("Invalid DESeq2 size factor in likelihood calculation.")
  }

  mu_mat <- outer(feature_mu_norm, sf, "*")
  dimnames(mu_mat) <- dimnames(raw_sub)

  y <- as.vector(raw_sub)
  mu <- as.vector(mu_mat)

  nb1 <- fit_nb_dispersion_model(y, mu, "NB1")
  nb2 <- fit_nb_dispersion_model(y, mu, "NB2")

  delta_ll <- nb2$logLik - nb1$logLik

  data.frame(
    arm = arm_label,
    region = region_label,
    n_features = length(feature_ids),
    n_observations = nb1$n_obs,
    alpha_NB1 = nb1$alpha,
    alpha_NB2 = nb2$alpha,
    logLik_NB1 = nb1$logLik,
    logLik_NB2 = nb2$logLik,
    two_delta_logLik_NB2_minus_NB1 = 2 * delta_ll,
    log10_likelihood_ratio_NB2_to_NB1 = delta_ll / log(10),
    preferred_model = ifelse(
      is.finite(delta_ll) & delta_ll > 0,
      "NB2",
      ifelse(is.finite(delta_ll) & delta_ll < 0, "NB1", "Tie")
    ),
    stringsAsFactors = FALSE
  )
}

build_corroboration_for_arm <- function(
    arm_label,
    raw_counts_arm,
    normalized_counts_arm,
    size_factors,
    rank_df,
    k) {

  nb_df <- compute_ranked_nb_metrics(
    raw_counts_arm = raw_counts_arm,
    normalized_counts_arm = normalized_counts_arm,
    ranked_feature_ids = rank_df$feature_id
  )

  nb_df$abs_pc1_loading <- rank_df$abs_pc1_loading
  nb_df <- assign_matched_regions(nb_df, k)

  matched_available <- isTRUE(attr(nb_df, "matched_block_available"))

  for (metric in c(
    "nb2_excess_signal",
    "nb2_nb1_contrast",
    "alpha_mu_signal"
  )) {
    nb_df[[paste0(metric, "_smooth")]] <- smooth_metric_for_display(
      nb_df$rank,
      nb_df[[metric]]
    )
  }

  region_summary <- nb_df %>%
    filter(corroboration_region %in% c("LEFT", "RIGHT")) %>%
    group_by(corroboration_region) %>%
    summarise(
      n = n(),
      median_nb2_excess = median(nb2_excess_signal, na.rm = TRUE),
      median_nb2_nb1 = median(nb2_nb1_contrast, na.rm = TRUE),
      median_alpha_mu = median(alpha_mu_signal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      arm = arm_label,
      matched_left_available = matched_available
    )

  lrt_rows <- list()

  regions_to_fit <- if (matched_available) c("LEFT", "RIGHT") else "RIGHT"

  for (region_name in regions_to_fit) {
    ids <- nb_df$feature_id[
      nb_df$corroboration_region == region_name
    ]

    lrt_rows[[region_name]] <- fit_nb1_nb2_region(
      raw_counts_arm = raw_counts_arm,
      normalized_counts_arm = normalized_counts_arm,
      size_factors = size_factors,
      feature_ids = ids,
      arm_label = arm_label,
      region_label = region_name
    )
  }

  list(
    rank_data = nb_df,
    region_summary = region_summary,
    likelihood = bind_rows(lrt_rows),
    matched_left_available = matched_available
  )
}

# =============================================================================
# FIGURES
# =============================================================================

add_rank_regions <- function(p, c1, c2, N) {
  p +
    annotate(
      "rect",
      xmin = 1, xmax = c1,
      ymin = -Inf, ymax = Inf,
      fill = COL$remainder, alpha = 0.42
    ) +
    annotate(
      "rect",
      xmin = c1, xmax = c2,
      ymin = -Inf, ymax = Inf,
      fill = COL$interval, alpha = 0.38
    ) +
    annotate(
      "rect",
      xmin = c2, xmax = N,
      ymin = -Inf, ymax = Inf,
      fill = COL$leading, alpha = 0.40
    )
}

make_boundary_lines <- function(c1, c2, cutoff_rank, k) {
  data.frame(
    rank = c(c1, c2, cutoff_rank),
    key = c("c1", "c2", paste0("Selected k* = ", k)),
    stringsAsFactors = FALSE
  )
}

make_pareto_panel <- function(scan_df, selected_k, comparison_name, method) {
  selected <- scan_df %>%
    filter(k == selected_k) %>%
    slice(1L)

  frontier <- scan_df %>%
    filter(is_pareto) %>%
    arrange(remainder_cross_n, good_n, k) %>%
    distinct(remainder_cross_n, good_n, .keep_all = TRUE)

  ggplot() +
    geom_path(
      data = scan_df,
      aes(remainder_cross_n, good_n, group = 1),
      color = COL$candidate_line,
      linewidth = 0.45,
      alpha = 0.55
    ) +
    geom_point(
      data = scan_df,
      aes(remainder_cross_n, good_n),
      color = COL$candidate_line,
      size = 1.0,
      alpha = 0.45
    ) +
    geom_path(
      data = frontier,
      aes(
        remainder_cross_n,
        good_n,
        color = "Pareto frontier",
        group = 1
      ),
      linewidth = 1.25
    ) +
    geom_point(
      data = selected,
      aes(
        remainder_cross_n,
        good_n,
        color = "Weighted optimum"
      ),
      shape = 23,
      fill = COL$selected,
      size = 4.7,
      stroke = 1.1
    ) +
    annotate(
      "label",
      x = selected$remainder_cross_n,
      y = selected$good_n,
      label = paste0(
        "k* = ", selected_k,
        "\nG = ", selected$good_n,
        "\nR = ", selected$remainder_cross_n,
        "\nU = ", formatC(selected$weighted_utility, digits = 4, format = "f")
      ),
      hjust = -0.05,
      vjust = 1.1,
      size = 3.1,
      label.size = 0.25,
      fill = "white"
    ) +
    scale_color_manual(
      name = NULL,
      values = c(
        "Pareto frontier" = COL$pareto,
        "Weighted optimum" = COL$selected
      )
    ) +
    labs(
      title = paste0("D. ", method, " / ", comparison_name, ": weighted Pareto cutoff"),
      subtitle = "Equal normalized weights: U(k) = Gnorm(k) - Rnorm(k)",
      x = "Opposite-arm Remainder crossings, R(k)",
      y = "Joint + permissible Disjoint, G(k)",
      caption = "Selection is restricted to the Pareto frontier."
    ) +
    theme_manuscript() +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))
}

make_cutoff_framework_figure <- function(
    method,
    comparison_name,
    control_group,
    treatment_group,
    group_results,
    knot_fit,
    scan_df,
    selected_k,
    c1,
    c2,
    out_file) {

  control <- group_results[[control_group]]$data
  treatment <- group_results[[treatment_group]]$data
  N <- nrow(control)
  cutoff_rank <- rank_cutoff_from_k(N, selected_k)

  long <- bind_rows(
    control %>% mutate(arm = control_group),
    treatment %>% mutate(arm = treatment_group)
  )

  boundary_df <- make_boundary_lines(c1, c2, cutoff_rank, selected_k)

  boundary_colors <- c(
    "c1" = COL$c1,
    "c2" = COL$c2,
    stats::setNames(COL$selected, paste0("Selected k* = ", selected_k))
  )

  boundary_types <- c(
    "c1" = "dashed",
    "c2" = "longdash",
    stats::setNames("dotdash", paste0("Selected k* = ", selected_k))
  )

  arm_colors <- c(
    stats::setNames(COL$control, control_group),
    stats::setNames(COL$treatment, treatment_group)
  )

  pca_subtitle <- if (method == "RawEVS") {
    "PCA input = log1p(raw counts); no library-size normalization before EVS"
  } else {
    "PCA input = log1p(DESeq2 median-of-ratios normalized counts)"
  }

  pA <- add_rank_regions(ggplot(), c1, c2, N) +
    geom_vline(
      data = boundary_df,
      aes(xintercept = rank, color = key, linetype = key),
      linewidth = 0.8
    ) +
    geom_line(
      data = long,
      aes(rank, abs_pc1_loading, color = arm),
      linewidth = 0.9
    ) +
    scale_color_manual(name = NULL, values = c(arm_colors, boundary_colors)) +
    scale_linetype_manual(
      name = NULL,
      values = boundary_types,
      na.translate = FALSE
    ) +
    labs(
      title = paste0("A. ", method, " / ", comparison_name, ": absolute PC1 loading"),
      subtitle = pca_subtitle,
      x = "PC1 rank: low |loading| to high |loading|",
      y = "|PC1 loading|"
    ) +
    theme_manuscript() +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 1)
    )

  pB <- add_rank_regions(ggplot(), c1, c2, N) +
    geom_vline(
      data = boundary_df,
      aes(xintercept = rank, color = key, linetype = key),
      linewidth = 0.8
    ) +
    geom_line(
      data = long,
      aes(rank, display_log1p_raw_empirical_variance, color = arm),
      linewidth = 0.9
    ) +
    scale_color_manual(name = NULL, values = c(arm_colors, boundary_colors)) +
    scale_linetype_manual(
      name = NULL,
      values = boundary_types,
      na.translate = FALSE
    ) +
    labs(
      title = paste0("B. ", method, " / ", comparison_name, ": raw-count variance geometry"),
      subtitle = "Raw-count variance displayed along the method-specific PC1 rank",
      x = "PC1 rank",
      y = "Smoothed log(1 + raw-count variance)"
    ) +
    theme_manuscript() +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 1)
    )

  fit_mat <- knot_fit$fitted

  fit_long <- bind_rows(lapply(seq_along(knot_fit$groups), function(j) {
    data.frame(
      rank = seq_len(N),
      arm = knot_fit$groups[j],
      fitted_D = fit_mat[, j],
      stringsAsFactors = FALSE
    )
  }))

  div_long <- long %>%
    select(rank, arm, cumulative_divergence)

  pC <- add_rank_regions(ggplot(), c1, c2, N) +
    geom_vline(
      data = boundary_df,
      aes(xintercept = rank, color = key, linetype = key),
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = 0,
      color = "grey55",
      linetype = "dotted",
      linewidth = 0.35
    ) +
    geom_line(
      data = div_long,
      aes(rank, cumulative_divergence, color = arm),
      linewidth = 0.65,
      alpha = 0.5
    ) +
    geom_line(
      data = fit_long,
      aes(rank, fitted_D, color = arm),
      linewidth = 1.15
    ) +
    scale_color_manual(name = NULL, values = c(arm_colors, boundary_colors)) +
    scale_linetype_manual(
      name = NULL,
      values = boundary_types,
      na.translate = FALSE
    ) +
    labs(
      title = paste0("C. ", method, " / ", comparison_name, ": PC1-variance divergence"),
      subtitle = "Two-knot fit shared by the two arms within this comparison",
      x = "PC1 rank",
      y = "D(r) = F_E(r) - F_P(r)"
    ) +
    theme_manuscript() +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 1)
    )

  pD <- make_pareto_panel(
    scan_df = scan_df,
    selected_k = selected_k,
    comparison_name = comparison_name,
    method = method
  )

  save_grid_2x2(
    list(pA, pB, pC, pD),
    out_file,
    width = 16,
    height = 12.5
  )
}

make_method_comparison_figure <- function(summary_df, out_file) {
  df <- summary_df %>%
    mutate(
      evs_method = factor(evs_method, levels = c("RawEVS", "NormEVS")),
      comparison = factor(comparison, levels = names(COMPARISONS))
    )

  p <- ggplot(
    df,
    aes(
      x = comparison,
      y = selected_k,
      fill = evs_method
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.75),
      width = 0.68
    ) +
    geom_text(
      aes(
        label = paste0(
          "k*=", selected_k,
          "\n", sprintf("%.1f%% LE", 100 * k_over_leading_edge)
        )
      ),
      position = position_dodge(width = 0.75),
      vjust = -0.35,
      size = 3.0,
      fontface = "bold"
    ) +
    scale_fill_manual(
      name = "EVS method",
      values = c(
        "RawEVS" = COL$raw_evs,
        "NormEVS" = COL$norm_evs
      )
    ) +
    labs(
      title = "RawEVS versus NormEVS empirical cutoffs",
      subtitle = "Each method is independently calibrated within each RT/ZT comparison",
      x = NULL,
      y = "Selected top-k* per arm",
      caption = "Labels show selected k* and the fraction of the fitted Leading-edge candidate domain selected."
    ) +
    theme_manuscript(base_size = 12.5) +
    expand_limits(y = max(df$selected_k, na.rm = TRUE) * 1.18)

  save_figure(p, out_file, width = 12.5, height = 7.5)
}

# =============================================================================
# METHODS / EXPORTS
# =============================================================================

write_methods_manuscript <- function(summary_df) {
  raw_text <- summary_df %>%
    filter(evs_method == "RawEVS") %>%
    transmute(x = paste0(comparison, " k*=", selected_k)) %>%
    pull(x) %>%
    paste(collapse = ", ")

  norm_text <- summary_df %>%
    filter(evs_method == "NormEVS") %>%
    transmute(x = paste0(comparison, " k*=", selected_k)) %>%
    pull(x) %>%
    paste(collapse = ", ")

  lines <- c(
    "# Dual RawEVS and NormEVS empirical cutoff methods",
    "",
    "## Comparison-specific preprocessing",
    "",
    paste0(
      "Empirical EVS cutoffs were estimated independently for RT0_ZT6, RT2_ZT8, RT4_ZT10, and RT8_ZT14 under two EVS representations. ",
      "Within each comparison, PASs with zero counts across all samples were removed. ",
      "DESeq2 median-of-ratios size factors were estimated independently within the comparison."
    ),
    "",
    "## RawEVS",
    "",
    paste0(
      "For RawEVS, arm-specific PCA was performed on log1p-transformed raw counts. ",
      "No library-size normalization was applied before PCA. ",
      "This transformation reduced extreme count-scale leverage while retaining the unnormalized raw-count representation."
    ),
    "",
    "## NormEVS",
    "",
    paste0(
      "For NormEVS, raw counts were first normalized using the DESeq2 median-of-ratios size factors estimated within the comparison. ",
      "The normalized counts were then transformed as log(1+x) before arm-specific PCA."
    ),
    "",
    "## PC1 ranking",
    "",
    paste0(
      "For both methods, PCA was centered and not feature-scaled. ",
      "PASs were ordered from lowest to highest absolute PC1 loading. ",
      "For PAS i in arm g, PC1 variance contribution was P_ig = lambda_1g * loading_ig^2."
    ),
    "",
    "## Size-factor-aware excess variance",
    "",
    paste0(
      "Variance-mass geometry was calculated on the DESeq2-normalized count scale. ",
      "Let Y_ij denote the raw count, s_j the DESeq2 size factor, and X_ij=Y_ij/s_j. ",
      "Under a Poisson reference with E[Y_ij]=s_j*mu_i, Var(X_ij)=mu_i/s_j. ",
      "Accordingly, the arm-specific Poisson reference variance was approximated as V_Pois,ig=mu_ig*mean_j(1/s_j). ",
      "A pooled within-group normalized empirical variance V_pool,i was calculated across the two groups in each comparison, and excess variance was E_ig=max(V_pool,i-V_Pois,ig,0)."
    ),
    "",
    "## PC1-variance divergence and regimes",
    "",
    paste0(
      "Within each EVS method and comparison, PC1 contribution and excess variance were normalized separately to rank-wise probability masses. ",
      "Their cumulative distributions were F_P,g(r) and F_E,g(r), and cumulative divergence was D_g(r)=F_E,g(r)-F_P,g(r). ",
      "A shared two-knot continuous linear spline was fitted jointly to the two arm-specific D_g(r) curves. ",
      "The fitted knots defined c1 and c2, with rank<c1 classified as Remainder, c1<=rank<=c2 as Divergence, and rank>c2 as Leading Edge."
    ),
    "",
    "## Weighted Pareto cutoff",
    "",
    paste0(
      "Candidate top-k values were restricted to 1<=k<=N-c2. ",
      "For each k, benefit G(k) was the number of Joint PASs plus Disjoint PASs whose opposite-arm rank lay in the Leading Edge or Divergence interval. ",
      "Cost R(k) was the number of Disjoint PASs whose opposite-arm rank lay in the Remainder. ",
      "Pareto-optimal candidates maximized G while minimizing R. ",
      "Within the Pareto frontier, G and R were min-max normalized and equal weights were used: U(k)=G_norm(k)-R_norm(k). ",
      "The comparison-specific empirical cutoff k* maximized U(k), with ties resolved by greater G, lower R, then larger k."
    ),
    "",
    "## Final membership",
    "",
    paste0(
      "After k* was selected independently for each method and comparison, the k* highest absolute-PC1-loading PASs were selected independently in the two arms. ",
      "The final Leading Edge was the union of these two top-k* sets. ",
      "Opposite-arm Remainder crossings contributed to the Pareto cost but did not override union membership."
    ),
    "",
    "## NB1/NB2 corroboration",
    "",
    paste0(
      "After k* was fixed, raw-count NB1/NB2 corroboration was calculated and did not contribute to cutoff selection. ",
      "When an immediately preceding equal-sized LEFT rank block was available, it was compared with the selected RIGHT block. ",
      "Conditional on size-factor-adjusted plug-in means, NB1 and NB2 each fitted one dispersion parameter by maximum likelihood. ",
      "NB1 used Var(Y)=mu+alpha*mu and size=mu/alpha; NB2 used Var(Y)=mu+alpha*mu^2 and size=1/alpha. ",
      "Model support was summarized by log-likelihoods, 2*(logLik_NB2-logLik_NB1), and log10(L_NB2/L_NB1)."
    ),
    "",
    "## Selected cutoffs",
    "",
    paste0("RawEVS: ", raw_text),
    "",
    paste0("NormEVS: ", norm_text)
  )

  writeLines(lines, file.path(OUT_ROOT, "Methods_Manuscript.md"))
}

write_figure_legends <- function() {
  lines <- c(
    "# Figure legends",
    "",
    "## Method-specific cutoff framework",
    "For each EVS method and RT/ZT comparison, Panel A shows absolute PC1 loading along the method-specific rank. RawEVS uses log1p(raw counts) for PCA; NormEVS uses log1p(DESeq2 median-of-ratios normalized counts). Panel B shows raw-count variance along the same method-specific rank. Panel C shows cumulative divergence D(r)=F_E(r)-F_P(r) and the fitted two-knot regime boundaries. Panel D shows all candidate cutoff coordinates, the Pareto frontier, and the selected weighted-Pareto k*. The fitted c1 and c2 boundaries and selected k* are estimated independently for each method and comparison.",
    "",
    "## RawEVS versus NormEVS summary",
    "The cross-method summary compares the independently calibrated RawEVS and NormEVS top-k* values for each RT/ZT comparison. Labels report k* and k*/(N-c2), the fraction of the fitted Leading-edge candidate domain selected by the Pareto optimum."
  )

  writeLines(lines, file.path(OUT_ROOT, "Figure_Legends.md"))
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

verify_outputs <- function(expected_paths) {
  missing <- expected_paths[
    !file.exists(expected_paths) |
      is.na(file.info(expected_paths)$size) |
      file.info(expected_paths)$size <= 0
  ]

  if (length(missing)) {
    stop(
      "Expected output files are missing or empty:\n",
      paste(missing, collapse = "\n")
    )
  }

  invisible(TRUE)
}

create_zip_from_files <- function(zip_path, files) {
  files <- unique(files[file.exists(files)])
  if (!length(files)) stop("No files available for zip: ", zip_path)

  root_norm <- normalizePath(OUT_ROOT, winslash = "/", mustWork = TRUE)
  file_norm <- normalizePath(files, winslash = "/", mustWork = TRUE)

  rel <- sub(
    paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", root_norm), "/?"),
    "",
    file_norm
  )

  if (file.exists(zip_path)) unlink(zip_path)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(OUT_ROOT)

  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(zipfile = zip_path, files = rel)
  } else {
    utils::zip(zipfile = zip_path, files = rel, flags = "-r9X")
  }

  setwd(old_wd)

  if (!file.exists(zip_path) || file.info(zip_path)$size <= 0) {
    stop("Zip creation failed: ", zip_path)
  }

  invisible(zip_path)
}

# =============================================================================
# RUN
# =============================================================================

input <- read_count_data(COUNT_FILE, GROUP_PATTERNS)
all_counts <- input$counts
annotation <- input$annotation

summary_rows <- list()
likelihood_rows <- list()
leading_rows <- list()
remainder_rows <- list()
pareto_cost_rows <- list()
all_scan_rows <- list()

expected_figure_paths <- character(0)
expected_table_paths <- character(0)

for (comparison_name in names(COMPARISONS)) {

  message("============================================================")
  message("Analyzing comparison: ", comparison_name)

  mapping <- COMPARISONS[[comparison_name]]
  control_group <- unname(mapping[["control"]])
  treatment_group <- unname(mapping[["treatment"]])

  control_idx <- grep(GROUP_PATTERNS[[control_group]], colnames(all_counts))
  treatment_idx <- grep(GROUP_PATTERNS[[treatment_group]], colnames(all_counts))

  if (length(control_idx) < 2L || length(treatment_idx) < 2L) {
    stop("Each comparison arm requires at least two samples: ", comparison_name)
  }

  sample_idx <- c(control_idx, treatment_idx)
  count_mat <- all_counts[, sample_idx, drop = FALSE]

  keep <- rowSums(count_mat) > 0
  count_mat <- count_mat[keep, , drop = FALSE]

  annotation_cmp <- annotation[
    match(rownames(count_mat), annotation$feature_id),
    ,
    drop = FALSE
  ]

  N <- nrow(count_mat)
  if (N < 20L) stop("Too few PASs after filtering: ", comparison_name)

  group_labels <- factor(
    c(
      rep(control_group, length(control_idx)),
      rep(treatment_group, length(treatment_idx))
    ),
    levels = c(control_group, treatment_group)
  )

  names(group_labels) <- colnames(count_mat)

  deseq <- normalize_deseq2_comparison(
    count_mat = count_mat,
    group_labels = group_labels
  )

  normalized_counts <- deseq$normalized_counts
  size_factors <- deseq$size_factors

  pooled <- compute_pooled_within_group_variance(
    normalized_counts = normalized_counts,
    group_labels = group_labels
  )

  for (method in EVS_METHODS) {

    message("  Method: ", method)

    group_results <- list()

    for (g in c(control_group, treatment_group)) {
      idx <- which(group_labels == g)

      group_results[[g]] <- compute_group_analysis(
        method = method,
        group_name = g,
        raw_counts_arm = count_mat[, idx, drop = FALSE],
        normalized_counts_arm = normalized_counts[, idx, drop = FALSE],
        size_factors_arm = size_factors[colnames(count_mat)[idx]],
        pooled_variance = pooled$variance
      )
    }

    knot_fit <- fit_shared_knots(group_results)
    c1 <- knot_fit$c1
    c2 <- knot_fit$c2

    raw_scan <- scan_pair_cutoffs(
      control_df = group_results[[control_group]]$data,
      treatment_df = group_results[[treatment_group]]$data,
      c1 = c1,
      c2 = c2,
      comparison_name = comparison_name,
      control_group = control_group,
      treatment_group = treatment_group,
      method = method
    )

    opt <- select_weighted_pareto_optimum(
      raw_scan,
      good_col = "good_n",
      cost_col = "remainder_cross_n",
      benefit_weight = BENEFIT_WEIGHT,
      contamination_weight = CONTAMINATION_WEIGHT
    )

    scan_df <- opt$scan
    selected_k <- opt$selected_k
    cutoff_rank <- rank_cutoff_from_k(N, selected_k)

    class_df <- classify_pair_at_k(
      control_df = group_results[[control_group]]$data,
      treatment_df = group_results[[treatment_group]]$data,
      k = selected_k,
      c1 = c1,
      c2 = c2,
      comparison_name = comparison_name,
      control_group = control_group,
      treatment_group = treatment_group,
      method = method
    )

    class_summary <- summarize_classification(class_df)

    lead_ids <- unique(class_df$feature_id)
    rem_ids <- setdiff(rownames(count_mat), lead_ids)

    lead_table <- class_df %>%
      left_join(annotation_cmp, by = "feature_id") %>%
      mutate(
        evs_membership = "LeadingEdge",
        pareto_remainder_crossing_cost = cross_into_remainder
      ) %>%
      select(
        evs_method,
        comparison,
        feature_id,
        gene_symbol,
        evs_membership,
        selected_k,
        cutoff_rank,
        control_group,
        treatment_group,
        base_class,
        analysis_class,
        control_rank,
        treatment_rank,
        control_region,
        treatment_region,
        opposite_region,
        control_top_k,
        treatment_top_k,
        pareto_remainder_crossing_cost
      )

    rem_control_rank <- make_rank_map(group_results[[control_group]]$data)
    rem_treatment_rank <- make_rank_map(group_results[[treatment_group]]$data)

    remainder_table <- data.frame(
      evs_method = method,
      comparison = comparison_name,
      feature_id = rem_ids,
      gene_symbol = annotation_cmp$gene_symbol[
        match(rem_ids, annotation_cmp$feature_id)
      ],
      evs_membership = "Remainder",
      selected_k = selected_k,
      cutoff_rank = cutoff_rank,
      control_rank = as.integer(rem_control_rank[rem_ids]),
      treatment_rank = as.integer(rem_treatment_rank[rem_ids]),
      stringsAsFactors = FALSE
    )

    cost_table <- lead_table %>%
      filter(pareto_remainder_crossing_cost)

    arm_results <- list()

    for (g in c(control_group, treatment_group)) {
      idx <- which(group_labels == g)

      arm_results[[g]] <- build_corroboration_for_arm(
        arm_label = g,
        raw_counts_arm = count_mat[, idx, drop = FALSE],
        normalized_counts_arm = normalized_counts[, idx, drop = FALSE],
        size_factors = size_factors[colnames(count_mat)[idx]],
        rank_df = group_results[[g]]$data,
        k = selected_k
      )
    }

    likelihood_table <- bind_rows(lapply(
      names(arm_results),
      function(g) arm_results[[g]]$likelihood
    )) %>%
      mutate(
        evs_method = method,
        comparison = comparison_name,
        selected_k = selected_k,
        cutoff_rank = cutoff_rank
      ) %>%
      select(
        evs_method,
        comparison,
        selected_k,
        cutoff_rank,
        everything()
      )

    region_summary <- bind_rows(lapply(
      names(arm_results),
      function(g) arm_results[[g]]$region_summary
    )) %>%
      mutate(
        evs_method = method,
        comparison = comparison_name,
        selected_k = selected_k,
        cutoff_rank = cutoff_rank
      ) %>%
      select(
        evs_method,
        comparison,
        selected_k,
        cutoff_rank,
        everything()
      )

    selected_scan_row <- scan_df %>%
      filter(k == selected_k) %>%
      slice(1L)

    lead_domain <- N - c2

    cutoff_row <- data.frame(
      evs_method = method,
      comparison = comparison_name,
      N = N,
      control_group = control_group,
      treatment_group = treatment_group,
      c1 = c1,
      c2 = c2,
      remainder_size = c1 - 1L,
      divergence_size = c2 - c1 + 1L,
      leading_edge_candidate_size = lead_domain,
      selected_k = selected_k,
      cutoff_rank = cutoff_rank,
      k_over_N = selected_k / N,
      k_over_leading_edge = selected_k / lead_domain,
      benefit_weight = BENEFIT_WEIGHT,
      contamination_weight = CONTAMINATION_WEIGHT,
      weighted_utility = selected_scan_row$weighted_utility,
      good_n = selected_scan_row$good_n,
      remainder_cross_n = selected_scan_row$remainder_cross_n,
      joint_n = selected_scan_row$joint_n,
      permissible_disjoint_n = selected_scan_row$permissible_disjoint_n,
      divergence_disjoint_n = selected_scan_row$disjoint_opposite_divergence_n,
      leading_edge_union_n = length(lead_ids),
      remainder_n = length(rem_ids),
      pc1_energy_control_selected_fraction =
        sum(tail(
          group_results[[control_group]]$data$pc1_variance_contribution,
          selected_k
        )) /
        sum(group_results[[control_group]]$data$pc1_variance_contribution),
      pc1_energy_treatment_selected_fraction =
        sum(tail(
          group_results[[treatment_group]]$data$pc1_variance_contribution,
          selected_k
        )) /
        sum(group_results[[treatment_group]]$data$pc1_variance_contribution),
      pooled_within_group_residual_df = pooled$residual_df,
      knot_fit_SSE = knot_fit$SSE,
      control_poisson_scale_factor =
        group_results[[control_group]]$poisson_scale_factor,
      treatment_poisson_scale_factor =
        group_results[[treatment_group]]$poisson_scale_factor,
      matched_left_control_available =
        arm_results[[control_group]]$matched_left_available,
      matched_left_treatment_available =
        arm_results[[treatment_group]]$matched_left_available,
      stringsAsFactors = FALSE
    )

    method_root <- file.path(OUT_ROOT, method)
    comp_fig_dir <- file.path(method_root, comparison_name, "Figures")
    comp_tab_dir <- file.path(method_root, comparison_name, "Tables")

    dir.create(comp_fig_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(comp_tab_dir, recursive = TRUE, showWarnings = FALSE)

    cutoff_fig <- file.path(
      comp_fig_dir,
      paste0("Figure_", method, "_", comparison_name, "_Cutoff_Framework.png")
    )

    make_cutoff_framework_figure(
      method = method,
      comparison_name = comparison_name,
      control_group = control_group,
      treatment_group = treatment_group,
      group_results = group_results,
      knot_fit = knot_fit,
      scan_df = scan_df,
      selected_k = selected_k,
      c1 = c1,
      c2 = c2,
      out_file = cutoff_fig
    )

    table_paths <- c(
      Cutoff_Summary = file.path(
        comp_tab_dir,
        paste0("Table_", method, "_", comparison_name, "_Cutoff_Summary.csv")
      ),
      Cutoff_Scan = file.path(
        comp_tab_dir,
        paste0("Table_", method, "_", comparison_name, "_Weighted_Pareto_Scan.csv")
      ),
      Leading_Edge = file.path(
        comp_tab_dir,
        paste0("Table_", method, "_", comparison_name, "_Leading_Edge_Sites.csv")
      ),
      Remainder = file.path(
        comp_tab_dir,
        paste0("Table_", method, "_", comparison_name, "_Remainder_Sites.csv")
      ),
      Pareto_Cost = file.path(
        comp_tab_dir,
        paste0("Table_", method, "_", comparison_name, "_Pareto_Remainder_Crossing_Cost_Sites.csv")
      ),
      PC1_NB_Rank = file.path(
        comp_tab_dir,
        paste0("Table_", method, "_", comparison_name, "_PC1_Variance_Rank_Data.csv")
      ),
      NB_Region_Summary = file.path(
        comp_tab_dir,
        paste0("Table_", method, "_", comparison_name, "_NB_Corroboration_Regions.csv")
      ),
      NB_Likelihood = file.path(
        comp_tab_dir,
        paste0("Table_", method, "_", comparison_name, "_NB1_NB2_Likelihood.csv")
      ),
      Size_Factors = file.path(
        comp_tab_dir,
        paste0("Table_", method, "_", comparison_name, "_DESeq2_Size_Factors.csv")
      )
    )

    write_csv(cutoff_row, table_paths[["Cutoff_Summary"]])
    write_csv(scan_df, table_paths[["Cutoff_Scan"]])
    write_csv(lead_table, table_paths[["Leading_Edge"]])
    write_csv(remainder_table, table_paths[["Remainder"]])
    write_csv(cost_table, table_paths[["Pareto_Cost"]])

    rank_table <- bind_rows(
      group_results[[control_group]]$data %>% mutate(arm = control_group),
      group_results[[treatment_group]]$data %>% mutate(arm = treatment_group)
    ) %>%
      left_join(annotation_cmp, by = "feature_id")

    write_csv(rank_table, table_paths[["PC1_NB_Rank"]])
    write_csv(region_summary, table_paths[["NB_Region_Summary"]])
    write_csv(likelihood_table, table_paths[["NB_Likelihood"]])

    write_csv(
      data.frame(
        comparison = comparison_name,
        sample = names(size_factors),
        group = as.character(group_labels[names(size_factors)]),
        size_factor = as.numeric(size_factors),
        stringsAsFactors = FALSE
      ),
      table_paths[["Size_Factors"]]
    )

    expected_figs_this <- cutoff_fig

    if (isTRUE(EXPORT_PDF)) {
      expected_figs_this <- c(
        expected_figs_this,
        sub("\\.png$", ".pdf", cutoff_fig)
      )
    }

    verify_outputs(c(expected_figs_this, unname(table_paths)))

    expected_figure_paths <- c(expected_figure_paths, expected_figs_this)
    expected_table_paths <- c(expected_table_paths, unname(table_paths))

    key <- paste(method, comparison_name, sep = "__")

    summary_rows[[key]] <- cutoff_row
    likelihood_rows[[key]] <- likelihood_table
    leading_rows[[key]] <- lead_table
    remainder_rows[[key]] <- remainder_table
    pareto_cost_rows[[key]] <- cost_table
    all_scan_rows[[key]] <- scan_df

    message(
      "    k*=", selected_k,
      " | c1=", c1,
      " | c2=", c2,
      " | k*/N=", signif(selected_k / N, 4),
      " | k*/LE=", signif(selected_k / lead_domain, 4),
      " | union=", length(lead_ids),
      " | G=", selected_scan_row$good_n,
      " | R=", selected_scan_row$remainder_cross_n
    )
  }
}

# =============================================================================
# CROSS-METHOD / CROSS-COMPARISON SUMMARY
# =============================================================================

cutoff_summary <- bind_rows(summary_rows)
likelihood_all <- bind_rows(likelihood_rows)
leading_all <- bind_rows(leading_rows)
remainder_all <- bind_rows(remainder_rows)
pareto_cost_all <- bind_rows(pareto_cost_rows)
scan_all <- bind_rows(all_scan_rows)

summary_cutoff_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_Dual_EVS_Empirical_Cutoffs.csv"
)

summary_likelihood_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_Dual_EVS_NB1_NB2_Likelihood_All.csv"
)

summary_lead_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_Dual_EVS_Leading_Edge_Sites_All.csv"
)

summary_rem_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_Dual_EVS_Remainder_Sites_All.csv"
)

summary_cost_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_Dual_EVS_Pareto_Remainder_Crossing_Cost_Sites_All.csv"
)

summary_scan_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_Dual_EVS_Weighted_Pareto_Scans_All.csv"
)

write_csv(cutoff_summary, summary_cutoff_path)
write_csv(likelihood_all, summary_likelihood_path)
write_csv(leading_all, summary_lead_path)
write_csv(remainder_all, summary_rem_path)
write_csv(pareto_cost_all, summary_cost_path)
write_csv(scan_all, summary_scan_path)

summary_fig_path <- file.path(
  SUMMARY_FIG_DIR,
  "Figure_RawEVS_vs_NormEVS_Cutoff_Summary.png"
)

make_method_comparison_figure(
  cutoff_summary,
  summary_fig_path
)

expected_summary_figs <- summary_fig_path
if (isTRUE(EXPORT_PDF)) {
  expected_summary_figs <- c(
    expected_summary_figs,
    sub("\\.png$", ".pdf", summary_fig_path)
  )
}

expected_figure_paths <- c(expected_figure_paths, expected_summary_figs)

expected_table_paths <- c(
  expected_table_paths,
  summary_cutoff_path,
  summary_likelihood_path,
  summary_lead_path,
  summary_rem_path,
  summary_cost_path,
  summary_scan_path
)

write_methods_manuscript(cutoff_summary)
write_figure_legends()

manifest_paths <- unique(c(
  expected_figure_paths,
  expected_table_paths,
  file.path(OUT_ROOT, "Methods_Manuscript.md"),
  file.path(OUT_ROOT, "Figure_Legends.md")
))

manifest <- data.frame(
  relative_path = sub(
    paste0(
      "^",
      gsub(
        "([][{}()+*^$|\\\\?.])",
        "\\\\\\1",
        normalizePath(OUT_ROOT, winslash = "/", mustWork = TRUE)
      ),
      "/?"
    ),
    "",
    normalizePath(
      manifest_paths,
      winslash = "/",
      mustWork = TRUE
    )
  ),
  type = ifelse(
    grepl("\\.(png|pdf)$", manifest_paths, ignore.case = TRUE),
    "figure",
    ifelse(
      grepl("\\.csv$", manifest_paths, ignore.case = TRUE),
      "table",
      "methods"
    )
  ),
  size_bytes = file.info(manifest_paths)$size,
  stringsAsFactors = FALSE
)

manifest_path <- file.path(
  OUT_ROOT,
  "Table_Export_Manifest.csv"
)

write_csv(manifest, manifest_path)
expected_table_paths <- c(expected_table_paths, manifest_path)

verify_outputs(c(
  expected_figure_paths,
  expected_table_paths,
  file.path(OUT_ROOT, "Methods_Manuscript.md"),
  file.path(OUT_ROOT, "Figure_Legends.md")
))

# =============================================================================
# ZIP ARCHIVES
# =============================================================================

raw_files <- c(
  expected_figure_paths[grepl("/RawEVS/", expected_figure_paths, fixed = TRUE)],
  expected_table_paths[grepl("/RawEVS/", expected_table_paths, fixed = TRUE)]
)

norm_files <- c(
  expected_figure_paths[grepl("/NormEVS/", expected_figure_paths, fixed = TRUE)],
  expected_table_paths[grepl("/NormEVS/", expected_table_paths, fixed = TRUE)]
)

summary_files <- c(
  expected_figure_paths[grepl("/Summary/", expected_figure_paths, fixed = TRUE)],
  expected_table_paths[grepl("/Summary/", expected_table_paths, fixed = TRUE)],
  file.path(OUT_ROOT, "Methods_Manuscript.md"),
  file.path(OUT_ROOT, "Figure_Legends.md"),
  manifest_path
)

raw_zip <- file.path(OUT_ROOT, "RawEVS_All_Outputs.zip")
norm_zip <- file.path(OUT_ROOT, "NormEVS_All_Outputs.zip")
summary_zip <- file.path(OUT_ROOT, "Dual_EVS_Summary_Outputs.zip")
all_zip <- file.path(OUT_ROOT, "Dual_RawEVS_NormEVS_All_Outputs.zip")

create_zip_from_files(raw_zip, raw_files)
create_zip_from_files(norm_zip, norm_files)
create_zip_from_files(summary_zip, summary_files)

create_zip_from_files(
  all_zip,
  c(
    expected_figure_paths,
    expected_table_paths,
    file.path(OUT_ROOT, "Methods_Manuscript.md"),
    file.path(OUT_ROOT, "Figure_Legends.md")
  )
)

verify_outputs(c(raw_zip, norm_zip, summary_zip, all_zip))

# =============================================================================
# FINAL CONSOLE SUMMARY
# =============================================================================

message("============================================================")
message("DUAL RawEVS + NormEVS ANALYSIS COMPLETE")
message("")
message("RawEVS PCA input: log1p(raw counts)")
message("NormEVS PCA input: log1p(DESeq2 median-of-ratios normalized counts)")
message("Each method has independent comparison-specific c1/c2 and k*.")
message("")

for (method in EVS_METHODS) {
  sub <- cutoff_summary %>% filter(evs_method == method)

  message(
    method, " selected k*: ",
    paste(
      sub$comparison,
      sub$selected_k,
      sep = "=",
      collapse = "; "
    )
  )
}

message("")
message("Cross-method summary: ", summary_cutoff_path)
message("RawEVS ZIP: ", raw_zip)
message("NormEVS ZIP: ", norm_zip)
message("Summary ZIP: ", summary_zip)
message("Complete ZIP: ", all_zip)
message("Methods: ", file.path(OUT_ROOT, "Methods_Manuscript.md"))
message("Figure legends: ", file.path(OUT_ROOT, "Figure_Legends.md"))
message("============================================================")
