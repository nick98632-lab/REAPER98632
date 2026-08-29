#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

SCRIPT_BUILD <- "EMPIRICAL_CUTOFF_COMPARISON_SPECIFIC_KSTAR_v9"

# =============================================================================
# EMPIRICAL EVS CUTOFF: WEIGHTED PARETO SELECTION + NB LOG-RATIO EVIDENCE
# =============================================================================
#
# SCOPE
# -----
# One empirical eigenvector-splitting (EVS) cutoff k* is estimated separately
# for each comparison:
#
#   RT0_ZT6, RT2_ZT8, RT4_ZT10, RT8_ZT14.
#
# The feature universe, the DESeq2 normalization and variance reference, and
# the c1/c2 regime boundaries are estimated once across the full 40-sample
# RT/ZT matrix and define a rank geometry common to all eight arms. Within that
# geometry each comparison receives its own weighted-Pareto scan and its own
# k*.
#
# PREPROCESSING AND PC1 RANKING
# -----------------------------
# PASs with zero counts across all 40 samples are removed once. Within each
# arm, raw counts are library-size normalized to counts per million (CPM) and
# log1p-transformed, and PCA is applied to that expression matrix with
# centering and without feature scaling. The PC1 loading values are not
# transformed: PASs are ranked from lowest to highest absolute PC1 loading.
# DESeq2 median-of-ratios size factors are estimated once across all 40 samples
# and supply the pooled within-group negative-binomial variance term.
#
# CUTOFF CALCULATION
# ------------------
# For PAS i in arm g, the PC1 variance contribution is
#
#   P_ig = lambda_1g * loading_ig^2.
#
# DESeq2-normalized counts give the pooled within-group variance V_pool,i, and
# for each arm the excess-over-Poisson variance is
#
#   E_ig = max(V_pool,i - mu_ig, 0).
#
# P and E are normalized to rank-wise probability masses and accumulated along
# the absolute-PC1-loading rank:
#
#   D_g(r) = F_E,g(r) - F_P,g(r).
#
# A single two-knot continuous linear spline is fitted jointly to the D_g(r)
# curves of all eight arms. The fitted knots define:
#
#   rank < c1          : Remainder regime
#   c1 <= rank <= c2   : Divergence interval
#   rank > c2          : Leading-edge regime
#
# Candidate top-k values are restricted to 1 <= k <= N-c2. For each k,
#
#   G(k) = Joint
#          + Disjoint with opposite arm in Leading Edge
#          + Disjoint with opposite arm in Divergence
#
#   R(k) = Disjoint with opposite arm in Remainder.
#
# The Pareto frontier maximizes G while minimizing R. On the frontier, G and R
# are min-max normalized and the equal-weight utility is
#
#   U(k) = G_norm(k) - R_norm(k).
#
# The empirical cutoff is the Pareto-optimal k that maximizes U(k), with ties
# resolved by greater G, then lower R, then larger k.
#
# EVS MEMBERSHIP
# --------------
# k* is applied independently to both arm-specific absolute-PC1-loading
# rankings. The Leading Edge is the union of the two top-k* sets; Joint and all
# Disjoint PASs are members. The Remainder is the complement of that union.
# Opposite-arm Remainder crossings enter the Pareto cost R(k) and are reported
# as a flag on the membership table.
#
# NB SCALING CORROBORATION
# ------------------------
# After c1, c2 and k* are fixed, each arm's selected top-k* block (RIGHT) is
# compared with the matched equal-sized block immediately below it (LEFT).
# Under
#
#   E = alpha * mu^p,   fitted as   log E = log alpha + p log mu,
#
# the exponent p is estimated within each block with a 95% confidence interval.
# p near 1 is NB1-like and p near 2 is NB2-like.
#
# The dispersion change is tested under NB2, Var(Y) = mu + alpha*mu^2, as a
# nested hypothesis: H0 holds one shared alpha across both blocks and H1 gives
# each block its own. H1 adds exactly one free parameter, so
#
#   LR = 2*(logLik_H1 - logLik_H0) ~ chi-square(df = 1).
#
# A larger exponent in the selected block, with alpha_RIGHT > alpha_LEFT and a
# large LR, is the evidence that the regime changes at the selected cutoff.
#
# OUTPUTS
# -------
# Three manuscript figures in <OUT_ROOT>/Figures/:
#
#   Figure_1_Rank_Geometry.png            (A) excess variance along the rank,
#                                         (B) the two cumulative masses,
#                                         (C) divergence and the shared fit
#   Figure_2_Cutoff_Selection.png         G_norm, R_norm and U(k) across
#                                         candidate k with k* and its numbers,
#                                         one panel per comparison
#   Figure_3_NB_Scaling_Corroboration.png the exponent p in the matched blocks
#                                         with LR and its p-value, one row per
#                                         arm
#
# Tables:
#
#   Summary/Tables/Table_Empirical_Cutoffs.csv   k*, c1, c2, G, R, U per
#                                                comparison
#   Summary/Tables/Table_NB_Evidence.csv         exponents with CIs, alphas,
#                                                LR and p per arm
#   <comparison>/Tables/Table_..._EVS_Membership.csv   the PAS assignment
#
# Methods_Manuscript.md and Figure_Legends.md are written from the same
# constants used by the code.
# =============================================================================

# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT   <- "/root/REAPER98632/exports/empirical_cutoff_comparison_specific_kstar_v9"

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

# Equal normalized weights in the weighted Pareto utility.
BENEFIT_WEIGHT       <- 1.0
CONTAMINATION_WEIGHT <- 1.0

# Numerical optimization uses a coarse starting grid and then refits the
# selected knot solution on every rank. This value controls computation only;
# the final c1/c2 fit is evaluated on the full rank series.
KNOT_COARSE_GRID_POINTS <- 5000L

# Figure export.
PNG_DPI <- 360
EXPORT_PDF <- TRUE

# Likelihood fitting bounds on log(alpha).
NB_LOG_ALPHA_LOWER <- -14
NB_LOG_ALPHA_UPPER <- 8


dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
FIG_DIR <- file.path(OUT_ROOT, "Figures")
SUMMARY_TAB_DIR <- file.path(OUT_ROOT, "Summary", "Tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# COLORS / FIGURE THEME
# =============================================================================

COL <- list(
  candidate_line = "#A7A7A7",
  utility        = "#5E3C99",
  benefit        = "#1B7837",
  cost           = "#D95F02",
  selected       = "#B5179E",
  control        = "#386CB0",
  treatment      = "#159D91",
  remainder      = "#DCE6F2",
  interval       = "#FFF0B3",
  leading        = "#D8F3E7",
  left           = "#386CB0",
  right          = "#159D91"
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
      legend.spacing.y = unit(2, "pt"),
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
      print(
        plots[[i]],
        vp = grid::viewport(layout.pos.row = r, layout.pos.col = c)
      )
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
# INPUT / NORMALIZATION
# =============================================================================

read_count_data <- function(path, group_patterns) {
  if (!file.exists(path)) stop("Count file does not exist: ", path)

  raw_df <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

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

  list(
    counts = count_mat,
    annotation = annotation
  )
}

assign_groups <- function(sample_names, group_patterns) {
  assigned <- rep(NA_character_, length(sample_names))

  for (group_name in names(group_patterns)) {
    idx <- grep(group_patterns[[group_name]], sample_names)
    if (length(idx) > 0L && any(!is.na(assigned[idx]))) {
      stop("At least one sample matched more than one group pattern.")
    }
    assigned[idx] <- group_name
  }

  if (any(is.na(assigned))) {
    stop("Unassigned samples: ", paste(sample_names[is.na(assigned)], collapse = ", "))
  }

  out <- factor(assigned, levels = names(group_patterns))
  names(out) <- sample_names
  out
}

# Library-size normalization to counts per million followed by log1p. This is
# the expression matrix that enters PCA for the ranking step.
normalize_cpm_log1p <- function(count_mat_arm) {
  lib_size <- colSums(count_mat_arm, na.rm = TRUE)
  lib_size[!is.finite(lib_size) | lib_size <= 0] <- 1

  cpm <- sweep(
    count_mat_arm,
    2L,
    lib_size / 1e6,
    "/"
  )

  log1p(cpm)
}

normalize_deseq2_comparison <- function(count_mat, group_labels) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required for normalized-before-EVS cutoff calibration.")
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
      message("Default DESeq2 size factors failed; using type='poscounts'.")
      DESeq2::estimateSizeFactors(dds, type = "poscounts")
    }
  )

  list(
    normalized_counts = DESeq2::counts(dds, normalized = TRUE),
    size_factors = DESeq2::sizeFactors(dds)
  )
}


compute_pooled_within_group_variance <- function(
    normalized_counts,
    group_labels) {

  groups <- levels(factor(group_labels))

  sse <- rep(0, nrow(normalized_counts))
  residual_df <- 0L

  for (g in groups) {
    idx <- which(group_labels == g)

    if (length(idx) < 2L) next

    xg <- normalized_counts[, idx, drop = FALSE]
    mu_g <- rowMeans(xg)

    resid_g <- sweep(
      xg,
      1L,
      mu_g,
      "-"
    )

    sse <- sse + rowSums(resid_g^2)
    residual_df <- residual_df + length(idx) - 1L
  }

  if (residual_df < 2L) {
    stop("Pooled residual degrees of freedom < 2.")
  }

  V <- sse / residual_df
  V[!is.finite(V)] <- 0
  V <- pmax(V, 0)

  names(V) <- rownames(normalized_counts)

  list(
    variance = V,
    residual_df = residual_df
  )
}

compute_pc1_rank <- function(rank_matrix_arm) {
  pca <- stats::prcomp(
    t(rank_matrix_arm),
    center = TRUE,
    scale. = FALSE,
    rank. = 1
  )

  loading <- pca$rotation[, 1L]
  loading[!is.finite(loading)] <- 0

  abs_loading <- abs(loading)

  rank_order <- order(
    abs_loading,
    decreasing = FALSE
  )

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

compute_raw_variance_geometry <- function(
    raw_counts_arm,
    rank_order) {

  raw_mu <- rowMeans(raw_counts_arm, na.rm = TRUE)

  raw_var <- apply(
    raw_counts_arm,
    1L,
    stats::var,
    na.rm = TRUE
  )

  raw_mu[!is.finite(raw_mu)] <- 0
  raw_var[!is.finite(raw_var)] <- 0
  raw_mu <- pmax(raw_mu, 0)
  raw_var <- pmax(raw_var, 0)

  raw_mu_ranked <- raw_mu[rank_order]
  raw_var_ranked <- raw_var[rank_order]

  data.frame(
    rank = seq_along(rank_order),
    raw_empirical_mean = raw_mu_ranked,
    raw_empirical_variance = raw_var_ranked,
    stringsAsFactors = FALSE
  )
}


compute_group_analysis <- function(
    group_name,
    raw_counts_arm,
    normalized_counts_arm,
    pooled_variance) {

  # PCA is applied to log1p(CPM) expression. The PC1 loading values themselves
  # are not transformed; PASs are ranked by absolute PC1 loading.
  rank_matrix <- normalize_cpm_log1p(raw_counts_arm)

  pc1 <- compute_pc1_rank(rank_matrix)
  rank_order <- pc1$rank_order

  geometry <- compute_raw_variance_geometry(
    raw_counts_arm = raw_counts_arm,
    rank_order = rank_order
  )

  mu_norm <- rowMeans(normalized_counts_arm, na.rm = TRUE)
  mu_norm[!is.finite(mu_norm)] <- 0
  mu_norm <- pmax(mu_norm, 0)

  P_ranked <- pc1$pc1_variance_contribution[rank_order]
  mu_ranked <- mu_norm[rank_order]
  V_pool_ranked <- pooled_variance[rank_order]

  E_ranked <- pmax(V_pool_ranked - mu_ranked, 0)

  P_total <- sum(P_ranked)
  E_total <- sum(E_ranked)

  if (!is.finite(P_total) || P_total <= 0) {
    stop("PC1 variance mass undefined for group ", group_name)
  }

  if (!is.finite(E_total) || E_total <= 0) {
    stop("NB excess-variance mass undefined for group ", group_name)
  }

  p_mass <- P_ranked / P_total
  q_mass <- E_ranked / E_total

  F_P <- cumsum(p_mass)
  F_E <- cumsum(q_mass)
  D <- F_E - F_P
  rank <- seq_along(rank_order)

  df <- geometry %>%
    mutate(
      group = group_name,
      feature_id = rownames(raw_counts_arm)[rank_order],
      pc1_loading = pc1$loading[rank_order],
      abs_pc1_loading = pc1$abs_loading[rank_order],
      pc1_eigenvalue = pc1$lambda1,
      pc1_variance_contribution = P_ranked,
      pc1_variance_mass = p_mass,
      normalized_group_mean = mu_ranked,
      pooled_normalized_variance = V_pool_ranked,
      nb_excess_variance = E_ranked,
      nb_excess_variance_mass = q_mass,
      cumulative_pc1_mass = F_P,
      cumulative_nb_mass = F_E,
      cumulative_divergence = D
    )

  list(
    data = df,
    rank_order = rank_order
  )
}


piecewise_basis <- function(x, c1, c2) {
  cbind(
    intercept = 1,
    x = x,
    hinge1 = pmax(x - c1, 0),
    hinge2 = pmax(x - c2, 0)
  )
}

piecewise_sse <- function(
    par,
    x,
    D_mat,
    min_gap) {

  c1 <- par[1L]
  c2 <- par[2L]

  if (
    !is.finite(c1) ||
    !is.finite(c2) ||
    c1 <= 0 ||
    c2 >= 1 ||
    c2 - c1 <= min_gap
  ) {
    return(1e100)
  }

  X <- piecewise_basis(
    x = x,
    c1 = c1,
    c2 = c2
  )

  coef <- tryCatch(
    qr.coef(
      qr(X),
      D_mat
    ),
    error = function(e) NULL
  )

  if (is.null(coef) || any(!is.finite(coef))) {
    return(1e100)
  }

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

  x_full <- (seq_len(N) - 1) / (N - 1)

  D_full <- do.call(
    cbind,
    lapply(
      group_results,
      function(z) z$data$cumulative_divergence
    )
  )

  colnames(D_full) <- groups

  opt_n <- min(KNOT_COARSE_GRID_POINTS, N)

  opt_idx <- unique(
    as.integer(
      round(
        seq(
          1,
          N,
          length.out = opt_n
        )
      )
    )
  )

  x_opt <- x_full[opt_idx]
  D_opt <- D_full[opt_idx, , drop = FALSE]

  min_gap <- max(
    4 / (N - 1),
    .Machine$double.eps^0.25
  )

  starts <- list(
    c(0.03, 0.97),
    c(0.08, 0.92),
    c(0.15, 0.85),
    c(0.25, 0.75),
    c(0.35, 0.65)
  )

  coarse <- lapply(
    starts,
    function(start) {
      stats::optim(
        par = start,
        fn = piecewise_sse,
        x = x_opt,
        D_mat = D_opt,
        min_gap = min_gap,
        method = "Nelder-Mead",
        control = list(
          maxit = 700,
          reltol = 1e-11
        )
      )
    }
  )

  values <- vapply(
    coarse,
    function(z) z$value,
    numeric(1)
  )

  best <- coarse[[which.min(values)]]

  refined <- stats::optim(
    par = best$par,
    fn = piecewise_sse,
    x = x_full,
    D_mat = D_full,
    min_gap = min_gap,
    method = "Nelder-Mead",
    control = list(
      maxit = 1000,
      reltol = 1e-12
    )
  )

  if (!is.finite(refined$value)) {
    stop("Shared-knot optimization failed.")
  }

  c1_rank <- as.integer(
    round(
      1 + refined$par[1L] * (N - 1)
    )
  )

  c2_rank <- as.integer(
    round(
      1 + refined$par[2L] * (N - 1)
    )
  )

  c1_rank <- max(
    2L,
    min(N - 2L, c1_rank)
  )

  c2_rank <- max(
    c1_rank + 1L,
    min(N - 1L, c2_rank)
  )

  c1_x <- (c1_rank - 1) / (N - 1)
  c2_x <- (c2_rank - 1) / (N - 1)

  X <- piecewise_basis(
    x_full,
    c1_x,
    c2_x
  )

  coef <- qr.coef(
    qr(X),
    D_full
  )

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

make_rank_map <- function(df) {
  stats::setNames(
    df$rank,
    df$feature_id
  )
}

rank_to_region <- function(rank, c1, c2) {
  ifelse(
    rank < c1,
    "Remainder",
    ifelse(
      rank <= c2,
      "Divergence",
      "LeadingEdge"
    )
  )
}

cumulative_activation <- function(depth, K) {
  depth <- as.integer(depth)

  keep <- (
    is.finite(depth) &
    depth >= 1L &
    depth <= K
  )

  if (!any(keep)) {
    return(rep(0L, K))
  }

  cumsum(
    tabulate(
      depth[keep],
      nbins = K
    )
  )
}

active_interval_count <- function(starts, ends, K) {
  # Count intervals active for start <= k < end.
  # end may equal K+1.
  if (length(starts) == 0L) {
    return(rep(0L, K))
  }

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
  ends <- pmin(
    ends[valid],
    K + 1L
  )

  if (length(starts) == 0L) {
    return(rep(0L, K))
  }

  diff_vec <- integer(K + 1L)

  start_tab <- tabulate(
    starts,
    nbins = K + 1L
  )

  end_tab <- tabulate(
    ends,
    nbins = K + 1L
  )

  diff_vec <- (
    diff_vec +
    start_tab -
    end_tab
  )

  cumsum(diff_vec)[seq_len(K)]
}

scan_pair_cutoffs <- function(
    control_df,
    treatment_df,
    c1,
    c2,
    comparison_name,
    control_group,
    treatment_group) {

  if (nrow(control_df) != nrow(treatment_df)) {
    stop(
      "Control and treatment rankings have different feature counts."
    )
  }

  if (!setequal(
    control_df$feature_id,
    treatment_df$feature_id
  )) {
    stop(
      "Control and treatment rankings do not contain the same feature IDs."
    )
  }

  N <- nrow(control_df)

  # Own-arm selection must originate strictly beyond c2.
  K <- as.integer(N - c2)

  if (K < 1L) {
    stop(
      "No candidate top-k depth exists beyond c2."
    )
  }

  rank_control <- make_rank_map(
    control_df
  )

  rank_treatment <- make_rank_map(
    treatment_df
  )

  ids <- control_df$feature_id

  rC <- as.integer(
    unname(
      rank_control[ids]
    )
  )

  rT <- as.integer(
    unname(
      rank_treatment[ids]
    )
  )

  if (
    any(!is.finite(rC)) ||
    any(!is.finite(rT))
  ) {
    stop(
      "Non-finite rank encountered in ",
      comparison_name
    )
  }

  # Entry depth:
  # rank N enters at k=1
  # rank 1 enters at k=N
  dC <- N - rC + 1L
  dT <- N - rT + 1L

  k <- seq_len(K)

  # -----------------------------------------------------------------------
  # JOINT
  # -----------------------------------------------------------------------
  # Joint membership activates once BOTH top-k selections contain the site.
  joint_depth <- pmax(
    dC,
    dT
  )

  joint_n <- cumulative_activation(
    joint_depth,
    K
  )

  # -----------------------------------------------------------------------
  # DISJOINT, OPPOSITE ARM ALSO IN LEADING EDGE
  # -----------------------------------------------------------------------
  # Control-only while treatment has not yet admitted the site.
  idx_le_C <- which(
    dC < dT &
    dC <= K &
    dT <= K
  )

  disjoint_control_opposite_le_n <- active_interval_count(
    starts = dC[idx_le_C],
    ends = dT[idx_le_C],
    K = K
  )

  # Treatment-only while control has not yet admitted the site.
  idx_le_T <- which(
    dT < dC &
    dT <= K &
    dC <= K
  )

  disjoint_treatment_opposite_le_n <- active_interval_count(
    starts = dT[idx_le_T],
    ends = dC[idx_le_T],
    K = K
  )

  # -----------------------------------------------------------------------
  # DISJOINT, OPPOSITE ARM IN DIVERGENCE INTERVAL
  # -----------------------------------------------------------------------
  # These sites are permissible.
  #
  # The selecting arm is in its own top-k subset and therefore >c2.
  # The opposite arm lies between c1 and c2 and never enters a top-k list
  # because candidate k is capped at N-c2.
  idx_div_C <- which(
    dC <= K &
    rT >= c1 &
    rT <= c2
  )

  disjoint_control_opposite_divergence_n <-
    cumulative_activation(
      dC[idx_div_C],
      K
    )

  idx_div_T <- which(
    dT <= K &
    rC >= c1 &
    rC <= c2
  )

  disjoint_treatment_opposite_divergence_n <-
    cumulative_activation(
      dT[idx_div_T],
      K
    )

  # -----------------------------------------------------------------------
  # HARD CROSS-REGIME CONTAMINATION: OPPOSITE ARM IN REMAINDER
  # -----------------------------------------------------------------------
  # Only r_opposite < c1 is penalized.
  idx_rem_C <- which(
    dC <= K &
    rT < c1
  )

  remainder_cross_control_n <- cumulative_activation(
    dC[idx_rem_C],
    K
  )

  idx_rem_T <- which(
    dT <= K &
    rC < c1
  )

  remainder_cross_treatment_n <- cumulative_activation(
    dT[idx_rem_T],
    K
  )

  # -----------------------------------------------------------------------
  # AGGREGATES
  # -----------------------------------------------------------------------

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

  # Benefit:
  # Joint + all permissible disjoint sites.
  good_n <- (
    joint_n +
    permissible_disjoint_n
  )

  # Cost:
  # only true crossing into opposite-arm remainder.
  remainder_cross_n <- (
    remainder_cross_control_n +
    remainder_cross_treatment_n
  )

  union_n <- (
    good_n +
    remainder_cross_n
  )

  # Independent union-size check.
  union_depth <- pmin(
    dC,
    dT
  )

  union_check <- cumulative_activation(
    union_depth,
    K
  )

  if (!all(
    union_n == union_check
  )) {
    stop(
      "Internal union-count mismatch in ",
      comparison_name
    )
  }

  data.frame(
    comparison = comparison_name,
    control_group = control_group,
    treatment_group = treatment_group,

    k = k,
    cutoff_rank = N - k + 1L,

    joint_n = joint_n,

    disjoint_control_opposite_le_n =
      disjoint_control_opposite_le_n,

    disjoint_treatment_opposite_le_n =
      disjoint_treatment_opposite_le_n,

    disjoint_opposite_le_n =
      disjoint_opposite_le_n,

    disjoint_control_opposite_divergence_n =
      disjoint_control_opposite_divergence_n,

    disjoint_treatment_opposite_divergence_n =
      disjoint_treatment_opposite_divergence_n,

    disjoint_opposite_divergence_n =
      disjoint_opposite_divergence_n,

    permissible_disjoint_control_n =
      permissible_disjoint_control_n,

    permissible_disjoint_treatment_n =
      permissible_disjoint_treatment_n,

    permissible_disjoint_n =
      permissible_disjoint_n,

    good_n = good_n,

    remainder_cross_control_n =
      remainder_cross_control_n,

    remainder_cross_treatment_n =
      remainder_cross_treatment_n,

    remainder_cross_n =
      remainder_cross_n,

    union_n = union_n,

    retained_fraction = ifelse(
      union_n > 0,
      good_n / union_n,
      NA_real_
    ),

    remainder_cross_fraction = ifelse(
      union_n > 0,
      remainder_cross_n / union_n,
      NA_real_
    ),

    divergence_disjoint_fraction = ifelse(
      union_n > 0,
      disjoint_opposite_divergence_n / union_n,
      NA_real_
    ),

    jaccard_top_k = ifelse(
      union_n > 0,
      joint_n / union_n,
      NA_real_
    ),

    stringsAsFactors = FALSE
  )
}

mark_pareto_frontier <- function(
    scan_df,
    good_col = "good_n",
    cost_col = "remainder_cross_n") {

  if (nrow(scan_df) < 1L) {
    stop(
      "Empty cutoff scan."
    )
  }

  tmp <- scan_df %>%
    transmute(
      row_id = row_number(),
      k = k,
      good = .data[[good_col]],
      cost = .data[[cost_col]]
    ) %>%
    arrange(
      cost,
      desc(good),
      desc(k)
    ) %>%
    group_by(cost) %>%
    slice(1L) %>%
    ungroup() %>%
    arrange(
      cost,
      desc(good)
    )

  running_best_before <- c(
    -Inf,
    head(
      cummax(tmp$good),
      -1L
    )
  )

  tmp$is_frontier_coord <- (
    tmp$good >
    running_best_before
  )

  frontier <- tmp %>%
    filter(
      is_frontier_coord
    ) %>%
    arrange(
      cost,
      good,
      k
    )

  key_all <- paste(
    scan_df[[cost_col]],
    scan_df[[good_col]],
    sep = "::"
  )

  key_frontier <- paste(
    frontier$cost,
    frontier$good,
    sep = "::"
  )

  out <- scan_df
  out$is_pareto <- (
    key_all %in%
    key_frontier
  )

  list(
    scan = out,
    frontier = frontier
  )
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

  marked <- mark_pareto_frontier(
    scan_df,
    good_col = good_col,
    cost_col = cost_col
  )

  frontier <- marked$frontier %>%
    arrange(
      cost,
      good,
      k
    )

  if (nrow(frontier) < 1L) {
    stop("No Pareto-optimal cutoff points were identified.")
  }

  good_range <- range(
    frontier$good,
    na.rm = TRUE
  )

  cost_range <- range(
    frontier$cost,
    na.rm = TRUE
  )

  normalize_good <- function(x) {
    if (diff(good_range) == 0) {
      rep(1, length(x))
    } else {
      (x - good_range[1L]) / diff(good_range)
    }
  }

  normalize_cost <- function(x) {
    if (diff(cost_range) == 0) {
      rep(0, length(x))
    } else {
      (x - cost_range[1L]) / diff(cost_range)
    }
  }

  frontier$good_norm <- normalize_good(
    frontier$good
  )

  frontier$remainder_norm <- normalize_cost(
    frontier$cost
  )

  frontier$weighted_utility <- (
    benefit_weight * frontier$good_norm -
    contamination_weight * frontier$remainder_norm
  )

  best_utility <- max(
    frontier$weighted_utility,
    na.rm = TRUE
  )

  chosen <- frontier %>%
    filter(
      abs(
        weighted_utility - best_utility
      ) < 1e-12
    ) %>%
    arrange(
      desc(good),
      cost,
      desc(k)
    ) %>%
    slice(1L)

  selected_k <- as.integer(
    chosen$k[1L]
  )

  out <- marked$scan

  # Use the same frontier-derived normalization for every candidate row so the
  # exported table can show the utility landscape. Selection itself is still
  # restricted to the Pareto frontier.
  out$good_norm <- normalize_good(
    out[[good_col]]
  )

  out$remainder_norm <- normalize_cost(
    out[[cost_col]]
  )

  out$weighted_utility <- (
    benefit_weight * out$good_norm -
    contamination_weight * out$remainder_norm
  )

  out$is_selected_weighted <- (
    out$k == selected_k
  )

  zero_idx <- which(
    out[[cost_col]] == 0
  )

  zero_max_k <- if (
    length(zero_idx) > 0L
  ) {
    max(
      out$k[zero_idx]
    )
  } else {
    NA_integer_
  }

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

classify_pair_at_k <- function(
    control_df,
    treatment_df,
    k,
    c1,
    c2,
    comparison_name,
    control_group,
    treatment_group) {

  if (
    nrow(control_df) !=
    nrow(treatment_df)
  ) {
    stop(
      "Control and treatment rankings have different feature counts."
    )
  }

  if (!setequal(
    control_df$feature_id,
    treatment_df$feature_id
  )) {
    stop(
      "Control and treatment rankings do not contain the same feature IDs."
    )
  }

  N <- nrow(control_df)
  Kmax <- N - c2

  if (
    k < 1L ||
    k > Kmax
  ) {
    stop(
      "k must satisfy 1 <= k <= N-c2. Received k=",
      k,
      "; N-c2=",
      Kmax,
      "."
    )
  }

  control_top <- tail(
    control_df$feature_id,
    k
  )

  treatment_top <- tail(
    treatment_df$feature_id,
    k
  )

  union_ids <- union(
    control_top,
    treatment_top
  )

  rank_control <- make_rank_map(
    control_df
  )

  rank_treatment <- make_rank_map(
    treatment_df
  )

  rC <- as.integer(
    unname(
      rank_control[union_ids]
    )
  )

  rT <- as.integer(
    unname(
      rank_treatment[union_ids]
    )
  )

  in_control <- (
    union_ids %in%
    control_top
  )

  in_treatment <- (
    union_ids %in%
    treatment_top
  )

  joint <- (
    in_control &
    in_treatment
  )

  control_only <- (
    in_control &
    !in_treatment
  )

  treatment_only <- (
    in_treatment &
    !in_control
  )

  control_region <- rank_to_region(
    rC,
    c1,
    c2
  )

  treatment_region <- rank_to_region(
    rT,
    c1,
    c2
  )

  base_class <- ifelse(
    joint,
    "Joint",
    ifelse(
      control_only,
      paste0(
        "Disjoint_",
        control_group
      ),
      paste0(
        "Disjoint_",
        treatment_group
      )
    )
  )

  opposite_region <- rep(
    NA_character_,
    length(union_ids)
  )

  opposite_region[
    control_only
  ] <- treatment_region[
    control_only
  ]

  opposite_region[
    treatment_only
  ] <- control_region[
    treatment_only
  ]

  disjoint_opposite_leading_edge <- (
    (control_only & rT > c2) |
    (treatment_only & rC > c2)
  )

  disjoint_opposite_divergence <- (
    (
      control_only &
      rT >= c1 &
      rT <= c2
    ) |
    (
      treatment_only &
      rC >= c1 &
      rC <= c2
    )
  )

  cross_into_remainder <- (
    (
      control_only &
      rT < c1
    ) |
    (
      treatment_only &
      rC < c1
    )
  )

  # EVS membership is the union of the two arm-specific top-k sets.
  # Opposite-arm Remainder crossings are a Pareto cost term only; they remain
  # in the final Leading Edge if selected by either arm.
  retained_for_analysis <- rep(TRUE, length(union_ids))

  analysis_class <- rep(
    NA_character_,
    length(union_ids)
  )

  analysis_class[
    joint
  ] <- "Joint"

  analysis_class[
    control_only &
    rT > c2
  ] <- paste0(
    "Disjoint_",
    control_group,
    "_OppositeLeadingEdge"
  )

  analysis_class[
    treatment_only &
    rC > c2
  ] <- paste0(
    "Disjoint_",
    treatment_group,
    "_OppositeLeadingEdge"
  )

  analysis_class[
    control_only &
    rT >= c1 &
    rT <= c2
  ] <- paste0(
    "Disjoint_",
    control_group,
    "_OppositeDivergence"
  )

  analysis_class[
    treatment_only &
    rC >= c1 &
    rC <= c2
  ] <- paste0(
    "Disjoint_",
    treatment_group,
    "_OppositeDivergence"
  )

  analysis_class[
    control_only &
    rT < c1
  ] <- paste0(
    "Disjoint_",
    control_group,
    "_OppositeRemainder_Cost"
  )

  analysis_class[
    treatment_only &
    rC < c1
  ] <- paste0(
    "Disjoint_",
    treatment_group,
    "_OppositeRemainder_Cost"
  )

  out <- data.frame(
    comparison = comparison_name,
    feature_id = union_ids,

    selected_k = as.integer(k),

    cutoff_rank = rank_cutoff_from_k(
      N,
      k
    ),

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

    disjoint_opposite_leading_edge =
      disjoint_opposite_leading_edge,

    disjoint_opposite_divergence =
      disjoint_opposite_divergence,

    cross_into_remainder =
      cross_into_remainder,

    retained_for_analysis =
      retained_for_analysis,

    stringsAsFactors = FALSE
  )

  if (
    any(
      is.na(
        out$analysis_class
      )
    )
  ) {
    stop(
      "Unclassified union site encountered in ",
      comparison_name
    )
  }

  out
}


# =============================================================================
# NB SCALING CORROBORATION AT THE SELECTED CUTOFF
# =============================================================================
#
# After k* is fixed, each arm's selected top-k* block (RIGHT) is compared with
# the matched equal-sized block immediately below it (LEFT).
#
#   E = alpha * mu^p     fitted as   log E = log alpha + p log mu
#
# p near 1 is NB1-like, p near 2 is NB2-like.
#
# The change is tested under NB2 as a nested hypothesis: H0 holds one shared
# dispersion across both blocks, H1 gives each block its own.
#
#   LR = 2 * (logLik_H1 - logLik_H0) ~ chi-square(df = 1)

assign_matched_regions <- function(rank_df, k) {
  N <- nrow(rank_df)
  cutoff_rank <- rank_cutoff_from_k(N, k)

  right_start <- cutoff_rank
  left_end <- right_start - 1L
  left_start <- left_end - k + 1L

  if (left_start < 1L) {
    stop("Matched LEFT block is unavailable for k=", k, "; N=", N, ".")
  }

  rank_df$block <- "Other"
  rank_df$block[rank_df$rank >= left_start & rank_df$rank <= left_end] <- "LEFT"
  rank_df$block[rank_df$rank >= right_start & rank_df$rank <= N] <- "RIGHT"

  rank_df
}

# NB scaling exponent within one block.
fit_nb_scaling_exponent <- function(block_df) {
  mu <- block_df$raw_empirical_mean
  excess <- pmax(block_df$raw_empirical_variance - mu, 0)

  ok <- is.finite(mu) & is.finite(excess) & mu > 0 & excess > 0

  if (sum(ok) < 20L) {
    return(c(p = NA_real_, lower = NA_real_, upper = NA_real_, n = sum(ok)))
  }

  fit <- stats::lm(log(excess[ok]) ~ log(mu[ok]))
  ci <- suppressWarnings(stats::confint(fit, 2L, level = 0.95))

  c(
    p = unname(stats::coef(fit)[2L]),
    lower = unname(ci[1L]),
    upper = unname(ci[2L]),
    n = sum(ok)
  )
}

nb2_block_loglik <- function(alpha, y, mu) {
  if (!is.finite(alpha) || alpha <= 0) return(-Inf)

  ll <- suppressWarnings(
    stats::dnbinom(x = y, mu = mu, size = max(1 / alpha, 1e-10), log = TRUE)
  )

  if (any(!is.finite(ll))) return(-Inf)
  sum(ll)
}

fit_nb2_alpha <- function(y, mu) {
  opt <- stats::optimize(
    function(log_alpha) -nb2_block_loglik(exp(log_alpha), y, mu),
    interval = c(NB_LOG_ALPHA_LOWER, NB_LOG_ALPHA_UPPER),
    tol = 1e-8
  )

  list(
    alpha = exp(opt$minimum),
    logLik = -opt$objective,
    at_bound = (opt$minimum <= NB_LOG_ALPHA_LOWER + 1e-6) ||
      (opt$minimum >= NB_LOG_ALPHA_UPPER - 1e-6)
  )
}

# Expected counts: each PAS's arm-specific normalized mean scaled by the
# per-sample DESeq2 size factor. Counts are rounded because the negative
# binomial likelihood is defined on the integers.
block_observations <- function(
    raw_counts_arm,
    normalized_counts_arm,
    size_factors,
    feature_ids) {

  raw_sub <- raw_counts_arm[feature_ids, , drop = FALSE]
  norm_sub <- normalized_counts_arm[feature_ids, , drop = FALSE]

  mu_feature <- pmax(rowMeans(norm_sub, na.rm = TRUE), 1e-10)
  mu_feature[!is.finite(mu_feature)] <- 1e-10

  sf <- as.numeric(size_factors[colnames(raw_sub)])
  if (any(!is.finite(sf)) || any(sf <= 0)) {
    stop("Invalid DESeq2 size factor in the NB likelihood calculation.")
  }

  y <- round(pmax(as.vector(raw_sub), 0))
  mu <- as.vector(outer(mu_feature, sf, "*"))

  ok <- is.finite(y) & is.finite(mu) & mu > 0
  list(y = y[ok], mu = mu[ok])
}

build_arm_nb_evidence <- function(
    arm_label,
    raw_counts_arm,
    normalized_counts_arm,
    size_factors,
    rank_df,
    k) {

  blocks <- assign_matched_regions(rank_df, k)

  left_df <- blocks[blocks$block == "LEFT", , drop = FALSE]
  right_df <- blocks[blocks$block == "RIGHT", , drop = FALSE]

  p_left <- fit_nb_scaling_exponent(left_df)
  p_right <- fit_nb_scaling_exponent(right_df)

  obs_left <- block_observations(
    raw_counts_arm, normalized_counts_arm, size_factors, left_df$feature_id
  )
  obs_right <- block_observations(
    raw_counts_arm, normalized_counts_arm, size_factors, right_df$feature_id
  )

  fit_left <- fit_nb2_alpha(obs_left$y, obs_left$mu)
  fit_right <- fit_nb2_alpha(obs_right$y, obs_right$mu)
  fit_pooled <- fit_nb2_alpha(
    c(obs_left$y, obs_right$y),
    c(obs_left$mu, obs_right$mu)
  )

  lr_stat <- max(2 * (fit_left$logLik + fit_right$logLik - fit_pooled$logLik), 0)
  lr_p <- stats::pchisq(lr_stat, df = 1L, lower.tail = FALSE)

  at_bound <- fit_left$at_bound || fit_right$at_bound || fit_pooled$at_bound

  data.frame(
    arm = arm_label,
    left_n_features = nrow(left_df),
    right_n_features = nrow(right_df),
    p_exponent_left = unname(p_left["p"]),
    p_exponent_left_lower = unname(p_left["lower"]),
    p_exponent_left_upper = unname(p_left["upper"]),
    p_exponent_left_n = unname(p_left["n"]),
    p_exponent_right = unname(p_right["p"]),
    p_exponent_right_lower = unname(p_right["lower"]),
    p_exponent_right_upper = unname(p_right["upper"]),
    p_exponent_right_n = unname(p_right["n"]),
    diff_p_exponent = unname(p_right["p"]) - unname(p_left["p"]),
    alpha_NB2_left = fit_left$alpha,
    alpha_NB2_right = fit_right$alpha,
    alpha_NB2_pooled = fit_pooled$alpha,
    logLik_separate = fit_left$logLik + fit_right$logLik,
    logLik_shared = fit_pooled$logLik,
    likelihood_ratio_stat = lr_stat,
    lr_df = 1L,
    lr_p = lr_p,
    direction = ifelse(
      at_bound,
      NA_character_,
      ifelse(fit_right$alpha > fit_left$alpha, "RIGHT_more_NB2", "LEFT_more_NB2")
    ),
    any_alpha_at_bound = at_bound,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# FIGURE BUILDERS
# =============================================================================

REGION_LEVELS <- c("Remainder", "Divergence", "Leading edge")

REGION_FILLS <- c(
  "Remainder"    = COL$remainder,
  "Divergence"   = COL$interval,
  "Leading edge" = COL$leading
)

# Figure 1. Three stacked panels realizing the method on the real data:
#   A. excess variance E along the absolute PC1 loading rank
#   B. the two normalized cumulative variance masses
#   C. their divergence and the shared two-transition model fixing c1 and c2
make_divergence_figure <- function(group_results, knot_fit, c1, c2, n_bins = 400L) {
  N <- nrow(group_results[[1L]]$data)
  arm_levels <- names(group_results)

  per_arm <- bind_rows(lapply(arm_levels, function(g) {
    d <- group_results[[g]]$data
    data.frame(
      rank = d$rank,
      arm = g,
      E = d$nb_excess_variance,
      F_P = d$cumulative_pc1_mass,
      F_E = d$cumulative_nb_mass,
      D = d$cumulative_divergence,
      stringsAsFactors = FALSE
    )
  })) %>%
    mutate(arm = factor(arm, levels = arm_levels))

  fitted <- bind_rows(lapply(seq_along(knot_fit$groups), function(j) {
    data.frame(
      rank = seq_len(N),
      arm = knot_fit$groups[j],
      D = knot_fit$fitted[, j],
      stringsAsFactors = FALSE
    )
  })) %>%
    mutate(arm = factor(arm, levels = arm_levels))

  region_df <- data.frame(
    xmin = c(1, c1, c2),
    xmax = c(c1, c2, N),
    region = factor(REGION_LEVELS, levels = REGION_LEVELS),
    stringsAsFactors = FALSE
  )

  region_label_df <- data.frame(
    x = c(c1 / 2, (c1 + c2) / 2, (c2 + N) / 2),
    label = c("REMAINDER", "DIVERGENCE INTERVAL", "LEADING EDGE"),
    stringsAsFactors = FALSE
  )

  base_layers <- function(p) {
    p +
      geom_rect(
        data = region_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = region),
        inherit.aes = FALSE,
        alpha = 0.45,
        show.legend = FALSE
      ) +
      scale_fill_manual(values = REGION_FILLS, drop = FALSE) +
      geom_vline(xintercept = c1, color = "#D73027", linetype = "dashed", linewidth = 0.8) +
      geom_vline(xintercept = c2, color = "#1A9850", linetype = "dashed", linewidth = 0.8)
  }

  arm_scale <- scale_color_brewer(name = "Arm", palette = "Dark2")

  # ---- Panel A: E along the rank, averaged within equal-width rank bins.
  bin_width <- max(1L, as.integer(ceiling(N / n_bins)))

  binned <- per_arm %>%
    mutate(bin = ((rank - 1L) %/% bin_width) + 1L) %>%
    group_by(arm, bin) %>%
    summarise(
      rank = mean(rank),
      E = mean(E, na.rm = TRUE),
      .groups = "drop"
    )

  pA <- base_layers(ggplot()) +
    geom_line(
      data = binned,
      aes(rank, E, color = arm),
      linewidth = 0.7
    ) +
    geom_text(
      data = region_label_df,
      aes(x = x, y = Inf, label = label),
      inherit.aes = FALSE,
      vjust = 1.6,
      size = 3.2,
      fontface = "bold",
      color = "grey25"
    ) +
    annotate(
      "label",
      x = N * 0.985, y = Inf, hjust = 1, vjust = 1.15,
      size = 3.0, label.size = 0.3, fill = "white", lineheight = 1.2,
      label = paste0(
        "P_i = lambda_1 * v_i1^2\n",
        "E_i = max(pooled variance - mean, 0)\n",
        "PASs ordered ascending by |v_i1|"
      )
    ) +
    arm_scale +
    labs(
      title = "A. Rank PASs by absolute PC1 loading",
      subtitle = paste0(
        "Excess-over-Poisson variance E along the rank, averaged within ",
        format(n_bins, big.mark = ","), " equal-width rank bins for display."
      ),
      x = "PC1 rank:  low |loading|  ->  high |loading|",
      y = "Excess variance E"
    ) +
    theme_manuscript(base_size = 11) +
    guides(color = guide_legend(nrow = 1))

  # ---- Panel B: the two normalized cumulative masses.
  mass_long <- bind_rows(
    per_arm %>% transmute(rank, arm, value = F_P, series = "F_P(r): PC1 variance mass"),
    per_arm %>% transmute(rank, arm, value = F_E, series = "F_E(r): NB excess-variance mass")
  ) %>%
    mutate(
      series = factor(
        series,
        levels = c("F_P(r): PC1 variance mass", "F_E(r): NB excess-variance mass")
      )
    )

  pB <- base_layers(ggplot()) +
    geom_line(
      data = mass_long,
      aes(rank, value, color = series, group = interaction(arm, series)),
      linewidth = 0.55,
      alpha = 0.75
    ) +
    annotate(
      "label",
      x = N * 0.015, y = 1.0, hjust = 0, vjust = 1,
      size = 3.0, label.size = 0.3, fill = "white", lineheight = 1.2,
      label = paste0(
        "p(r) = P(r) / sum_j P(j)        q(r) = E(r) / sum_j E(j)\n",
        "F_P(r) = sum_{j<=r} p(j)        F_E(r) = sum_{j<=r} q(j)"
      )
    ) +
    scale_color_manual(
      name = NULL,
      values = c(
        "F_P(r): PC1 variance mass" = "#3B5BA5",
        "F_E(r): NB excess-variance mass" = "#159D91"
      )
    ) +
    scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(
      title = "B. Convert variance into rank-wise distributions and cumulative mass",
      subtitle = "Both masses integrate to 1 and are therefore directly comparable. All eight arms shown.",
      x = "PC1 rank",
      y = "Cumulative mass"
    ) +
    theme_manuscript(base_size = 11) +
    guides(color = guide_legend(nrow = 1))

  # ---- Panel C: divergence and the shared two-transition fit.
  pC <- base_layers(ggplot()) +
    geom_hline(yintercept = 0, color = "grey45", linetype = "dotted", linewidth = 0.4) +
    geom_line(
      data = per_arm,
      aes(rank, D, group = arm),
      color = "grey35",
      linewidth = 0.4,
      alpha = 0.40
    ) +
    geom_line(
      data = fitted,
      aes(rank, D, color = arm),
      linewidth = 1.0
    ) +
    annotate(
      "label",
      x = N * 0.015, y = max(fitted$D, na.rm = TRUE), hjust = 0, vjust = 1,
      size = 3.0, label.size = 0.3, fill = "white", lineheight = 1.2,
      label = paste0(
        "D(r) = F_E(r) - F_P(r)\n",
        "D_g(x) = b0g + b1g*x + g1g*(x-c1)+ + g2g*(x-c2)+\n",
        "c1 and c2 fitted jointly across all 8 arms"
      )
    ) +
    annotate(
      "label", x = c1, y = -Inf, vjust = -0.25, size = 3.1, fontface = "bold",
      label.size = 0.25, fill = "white",
      label = paste0("c1 = ", format(c1, big.mark = ","))
    ) +
    annotate(
      "label", x = c2, y = -Inf, vjust = -0.25, size = 3.1, fontface = "bold",
      label.size = 0.25, fill = "white",
      label = paste0("c2 = ", format(c2, big.mark = ","))
    ) +
    arm_scale +
    labs(
      title = "C. Shared two-transition model fixes the regime boundaries",
      subtitle = paste0(
        "Grey: observed D(r) per arm. Coloured: the single jointly fitted spline.  ",
        "N = ", format(N, big.mark = ","),
        ",  candidate k = 1 to ", format(N - c2, big.mark = ",")
      ),
      x = "PC1 rank",
      y = "D(r) = F_E(r) - F_P(r)"
    ) +
    theme_manuscript(base_size = 11) +
    guides(color = guide_legend(nrow = 1))

  list(pA, pB, pC)
}

save_stacked <- function(plots, png_path, width, height) {
  draw_once <- function(device_fun) {
    device_fun()
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(nrow = length(plots), ncol = 1L)
      )
    )
    for (i in seq_along(plots)) {
      print(plots[[i]], vp = grid::viewport(layout.pos.row = i, layout.pos.col = 1L))
    }
    grDevices::dev.off()
  }

  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)

  draw_once(function() {
    grDevices::png(
      filename = png_path, width = width, height = height,
      units = "in", res = PNG_DPI, bg = "white"
    )
  })

  if (isTRUE(EXPORT_PDF)) {
    draw_once(function() {
      grDevices::pdf(
        file = sub("\\.png$", ".pdf", png_path),
        width = width, height = height, onefile = TRUE, useDingbats = FALSE
      )
    })
  }

  invisible(png_path)
}

# Figure 2 panel: the weighted-Pareto objective itself, for one comparison.
# Both normalized components and their difference are drawn on one axis so the
# figure shows the function being maximized, not just its argmax:
#
#   G(k) = Joint + Disjoint with opposite arm in Leading edge or Divergence
#   R(k) = Disjoint with opposite arm in Remainder
#   U(k) = G_norm(k) - R_norm(k)
make_pareto_panel <- function(scan_df, selected_k, comparison_name, max_candidate_k) {
  selected <- scan_df %>%
    filter(k == selected_k) %>%
    slice(1L)

  curves <- bind_rows(
    scan_df %>% transmute(k, value = good_norm, term = "G_norm(k)  benefit"),
    scan_df %>% transmute(k, value = remainder_norm, term = "R_norm(k)  cost"),
    scan_df %>% transmute(k, value = weighted_utility, term = "U(k) = G_norm - R_norm")
  ) %>%
    mutate(
      term = factor(
        term,
        levels = c(
          "G_norm(k)  benefit",
          "R_norm(k)  cost",
          "U(k) = G_norm - R_norm"
        )
      )
    )

  ggplot(curves, aes(k, value, color = term, linewidth = term)) +
    geom_hline(yintercept = 0, color = "grey60", linetype = "dotted", linewidth = 0.35) +
    geom_vline(
      xintercept = selected_k,
      color = COL$selected,
      linetype = "dashed",
      linewidth = 0.8
    ) +
    geom_line() +
    geom_point(
      data = selected %>% mutate(term = "U(k) = G_norm - R_norm"),
      aes(k, weighted_utility),
      inherit.aes = FALSE,
      shape = 23,
      fill = COL$selected,
      color = COL$selected,
      size = 4.4,
      stroke = 1.0
    ) +
    annotate(
      "label",
      x = -Inf,
      y = Inf,
      hjust = -0.03,
      vjust = 1.05,
      size = 2.9,
      label.size = 0.3,
      fill = "white",
      lineheight = 1.2,
      label = paste0(
        "G(k) = Joint + Disjoint with opposite arm in Leading edge or Divergence\n",
        "R(k) = Disjoint with opposite arm in Remainder\n",
        "U(k) = G_norm(k) - R_norm(k),  equal weights,  k* = argmax U(k)"
      )
    ) +
    annotate(
      "label",
      x = Inf,
      y = -Inf,
      hjust = 1.02,
      vjust = -0.08,
      size = 3.1,
      label.size = 0.25,
      fill = "white",
      lineheight = 1.15,
      label = paste0(
        "k* = ", format(selected_k, big.mark = ","), "\n",
        "G(k*) = ", format(selected$good_n, big.mark = ","), "\n",
        "R(k*) = ", format(selected$remainder_cross_n, big.mark = ","), "\n",
        "U(k*) = ", formatC(selected$weighted_utility, format = "f", digits = 4), "\n",
        "scanned k = 1 to ", format(max_candidate_k, big.mark = ",")
      )
    ) +
    scale_color_manual(
      name = NULL,
      values = c(
        "G_norm(k)  benefit" = COL$benefit,
        "R_norm(k)  cost" = COL$cost,
        "U(k) = G_norm - R_norm" = COL$utility
      )
    ) +
    scale_linewidth_manual(
      name = NULL,
      values = c(
        "G_norm(k)  benefit" = 0.7,
        "R_norm(k)  cost" = 0.7,
        "U(k) = G_norm - R_norm" = 1.3
      ),
      guide = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.04))) +
    labs(
      title = paste0(comparison_name, ": weighted Pareto optimization"),
      x = "Candidate top-k",
      y = "Normalized benefit, cost, and utility"
    ) +
    theme_manuscript(base_size = 12) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))
}

make_grid_figure <- function(panels, out_file) {
  if (length(panels) != 4L) {
    stop("This figure expects one panel per comparison.")
  }

  save_grid_2x2(
    panels,
    out_file,
    width = 14.0,
    height = 10.0
  )
}

# Figure 3. The corroboration step of the method: after the boundaries and k*
# are fixed, the NB scaling exponent is estimated in each matched block, and
# the dispersion change is tested by likelihood ratio. One row per arm, with
# the likelihood-ratio statistic and its p-value at that comparison's cutoff.
make_evidence_figure <- function(evidence_all) {
  df <- evidence_all %>%
    mutate(
      arm_label = paste0(comparison, "  |  ", arm, "   k* = ", format(selected_k, big.mark = ",")),
      stat_label = paste0(
        "LR = ", formatC(likelihood_ratio_stat, format = "f", digits = 1),
        ",  ",
        ifelse(
          !is.finite(lr_p) | is.na(lr_p),
          "p = NA",
          ifelse(
            lr_p < 1e-300,
            "p < 1e-300",
            paste0("p = ", formatC(lr_p, format = "e", digits = 1))
          )
        )
      )
    ) %>%
    arrange(comparison, arm)

  df$arm_label <- factor(df$arm_label, levels = rev(unique(df$arm_label)))

  region_colors <- c(
    "LEFT: matched block below k*" = COL$left,
    "RIGHT: selected top-k* block" = COL$right
  )

  ggplot(df, aes(y = arm_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#B2182B", linewidth = 0.6) +
    geom_vline(xintercept = 2, linetype = "dashed", color = "#2166AC", linewidth = 0.6) +
    annotate("text", x = 1, y = Inf, label = "p = 1   NB1-like", vjust = 1.4, hjust = -0.06,
             size = 3.1, fontface = "bold", color = "#B2182B") +
    annotate("text", x = 2, y = Inf, label = "p = 2   NB2-like", vjust = 1.4, hjust = -0.06,
             size = 3.1, fontface = "bold", color = "#2166AC") +
    geom_segment(
      aes(x = p_exponent_left, xend = p_exponent_right, yend = arm_label),
      color = "grey55",
      linewidth = 0.85,
      arrow = arrow(length = unit(0.11, "in"), type = "closed")
    ) +
    geom_errorbarh(
      aes(xmin = p_exponent_left_lower, xmax = p_exponent_left_upper),
      height = 0.16, color = COL$left, linewidth = 0.6
    ) +
    geom_errorbarh(
      aes(xmin = p_exponent_right_lower, xmax = p_exponent_right_upper),
      height = 0.16, color = COL$right, linewidth = 0.6
    ) +
    geom_point(aes(x = p_exponent_left, color = names(region_colors)[1]), size = 3.3) +
    geom_point(aes(x = p_exponent_right, color = names(region_colors)[2]), size = 3.3) +
    geom_text(
      aes(
        x = pmax(p_exponent_left_upper, p_exponent_right_upper),
        label = stat_label
      ),
      hjust = -0.12,
      size = 3.0
    ) +
    annotate(
      "label",
      x = -Inf, y = -Inf, hjust = -0.03, vjust = -0.2,
      size = 3.0, label.size = 0.3, fill = "white", lineheight = 1.2,
      label = paste0(
        "E = alpha * mu^p     fitted as   log E = log alpha + p log mu\n",
        "Bars are 95% CI on p.\n",
        "LR = 2(logLik separate alpha - logLik shared alpha) ~ chi-square(1)"
      )
    ) +
    scale_color_manual(name = NULL, values = region_colors) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.30))) +
    labs(
      title = "NB scaling corroboration at each selected cutoff",
      subtitle = "At each comparison's own k*, RIGHT is the selected top-k* block and LEFT the matched equal-sized block below it.",
      x = "NB scaling exponent p",
      y = NULL
    ) +
    theme_manuscript(base_size = 11) +
    theme(legend.position = "bottom")
}


# =============================================================================
# METHODS / EXPORTS
# =============================================================================

write_methods_manuscript <- function(cutoff_summary, evidence_all) {
  cutoff_text <- paste(
    paste0(cutoff_summary$comparison, " k*=", cutoff_summary$selected_k),
    collapse = ", "
  )

  lines <- c(
    "# Empirical eigenvector-splitting cutoff",
    "",
    "## Ranking PASs by absolute PC1 loading",
    "",
    paste0(
      "An empirical eigenvector-splitting (EVS) cutoff k* was estimated separately for each of four comparisons: RT0_ZT6, RT2_ZT8, RT4_ZT10, and RT8_ZT14. ",
      "PASs with zero counts across the complete 40-sample RT/ZT matrix were removed once, giving a single feature universe of N PASs shared by all eight arms. ",
      "Within each arm, raw counts were library-size normalized to counts per million and log1p-transformed, and principal component analysis was applied to that matrix with centering and without feature scaling. ",
      "The loadings were not transformed; PASs were ordered ascending by absolute PC1 loading. ",
      "For PAS i, the PC1 variance contribution was P_i = lambda_1 * v_i1^2. ",
      "DESeq2 median-of-ratios size factors were estimated once across all 40 samples, pooled within-group variance was computed from the resulting normalized counts, and the excess-over-Poisson variance was E_i = max(pooled variance - mean, 0)."
    ),
    "",
    "## Rank-wise distributions and cumulative divergence",
    "",
    paste0(
      "P and E were each normalized to rank-wise masses p(r) = P(r)/sum_j P(j) and q(r) = E(r)/sum_j E(j), so both integrate to 1 and are directly comparable. ",
      "The masses were accumulated into F_P(r) = sum_{j<=r} p(j) and F_E(r) = sum_{j<=r} q(j), and their difference gave the cumulative divergence D(r) = F_E(r) - F_P(r)."
    ),
    "",
    "## Shared two-transition model",
    "",
    paste0(
      "A two-knot continuous linear spline D_g(x) = beta_0g + beta_1g*x + gamma_1g*(x - c1)+ + gamma_2g*(x - c2)+ was fitted to the eight arm-specific divergence curves by least squares, ",
      "with the knots c1 and c2 shared across all eight arms and estimated numerically. ",
      "Ranks below c1 defined the Remainder regime, ranks from c1 through c2 the Divergence interval, and ranks above c2 the Leading-edge regime."
    ),
    "",
    "## Weighted Pareto cutoff",
    "",
    paste0(
      "Candidate top-k values were restricted to 1 <= k <= N - c2, so a selecting arm drew only from its own Leading-edge regime. ",
      "At each k, a PAS was Joint if both arms placed it in their top k and Disjoint if only one arm did. ",
      "The benefit G(k) was the number of Joint PASs plus Disjoint PASs whose opposite-arm rank fell in the Leading-edge regime or the Divergence interval; ",
      "the cost R(k) was the number of Disjoint PASs whose opposite-arm rank fell in the Remainder regime. ",
      "Candidates for which no other candidate simultaneously increased G and decreased R formed the Pareto frontier, on which G and R were min-max normalized and combined with equal weights as U(k) = G_norm(k) - R_norm(k). ",
      "The cutoff k* maximized U(k), with ties resolved by greater G, then lower R, then larger k. ",
      "It was applied independently to both arm-specific rankings: the Leading Edge was the union of the two top-k* sets and the Remainder its complement."
    ),
    "",
    "## NB scaling corroboration",
    "",
    paste0(
      "Once the boundaries and k* were fixed, each arm's selected top-k* block (RIGHT) was compared with the matched equal-sized block immediately below it (LEFT). ",
      "Under E = alpha * mu^p, the scaling exponent was estimated within each block by regressing log E on log mu across the PASs of that block, with a 95% confidence interval on p; ",
      "p near 1 corresponds to NB1-like scaling and p near 2 to NB2-like scaling."
    ),
    "",
    paste0(
      "The dispersion change was tested as a nested hypothesis under NB2 with Var(Y) = mu + alpha*mu^2. ",
      "Expected counts were each PAS's arm-specific normalized mean scaled by the sample DESeq2 size factor, and observed counts were rounded to non-negative integers. ",
      "H0 held one shared alpha across the pooled LEFT and RIGHT observations and H1 gave each block its own alpha, adding exactly one free parameter, ",
      "so the likelihood ratio LR = 2*(logLik_H1 - logLik_H0) was referred to a chi-square distribution with one degree of freedom. ",
      "The exponent and the likelihood-ratio test were computed after k* was fixed and did not enter the Pareto optimization."
    ),
    "",
    "## Results",
    "",
    cutoff_text,
    "",
    paste0(
      "The NB scaling exponent was larger in the selected block than in the matched block below it in ",
      sum(evidence_all$diff_p_exponent > 0, na.rm = TRUE),
      " of ", nrow(evidence_all), " arms, rising from a median of ",
      formatC(stats::median(evidence_all$p_exponent_left, na.rm = TRUE), format = "f", digits = 2),
      " to ",
      formatC(stats::median(evidence_all$p_exponent_right, na.rm = TRUE), format = "f", digits = 2),
      ". The fitted NB2 dispersion was higher in the selected block in ",
      sum(evidence_all$alpha_NB2_right > evidence_all$alpha_NB2_left, na.rm = TRUE),
      " of ", nrow(evidence_all),
      " arms, with likelihood-ratio statistics from ",
      formatC(min(evidence_all$likelihood_ratio_stat, na.rm = TRUE), format = "f", digits = 1),
      " to ",
      formatC(max(evidence_all$likelihood_ratio_stat, na.rm = TRUE), format = "f", digits = 1),
      " on one degree of freedom and a maximum p-value of ",
      formatC(max(evidence_all$lr_p, na.rm = TRUE), format = "e", digits = 2),
      "."
    )
  )

  writeLines(lines, file.path(OUT_ROOT, "Methods_Manuscript.md"))
}

write_figure_legends <- function() {
  lines <- c(
    "# Figure legends",
    "",
    "## Figure 1. Rank geometry and the shared regime boundaries",
    paste0(
      "(A) PASs ordered ascending by absolute PC1 loading, showing the excess-over-Poisson variance E averaged within equal-width rank bins for display. ",
      "The PC1 variance contribution is P_i = lambda_1 * v_i1^2 and E_i = max(pooled within-group variance - mean, 0). ",
      "(B) P and E normalized to rank-wise masses p(r) and q(r) and accumulated into F_P(r) and F_E(r). Both integrate to 1 and are therefore directly comparable; all eight RT/ZT arms are shown. ",
      "(C) Cumulative divergence D(r) = F_E(r) - F_P(r) in grey per arm, with the two-transition spline D_g(x) = beta_0g + beta_1g*x + gamma_1g*(x - c1)+ + gamma_2g*(x - c2)+ fitted jointly across all eight arms in colour. ",
      "The shared knots c1 and c2, annotated with their fitted values, partition the rank axis into the Remainder regime, the Divergence interval, and the Leading-edge regime, shaded and labelled in every panel."
    ),
    "",
    "## Figure 2. Weighted Pareto optimization of the cutoff",
    paste0(
      "One panel per comparison, showing the objective rather than only its solution. ",
      "At each candidate top-k, G(k) counts Joint PASs plus Disjoint PASs whose opposite-arm rank lies in the Leading-edge regime or the Divergence interval, and R(k) counts Disjoint PASs whose opposite-arm rank lies in the Remainder regime. ",
      "G and R are min-max normalized over the Pareto frontier and combined with equal weights as U(k) = G_norm(k) - R_norm(k). ",
      "Green and orange trace the normalized components and purple the utility; the diamond and dashed line mark k* = argmax U(k). ",
      "Each panel is annotated with the definitions and with k*, G(k*), R(k*), U(k*), and the scanned candidate range."
    ),
    "",
    "## Figure 3. NB scaling corroboration at each selected cutoff",
    paste0(
      "One row per arm, labelled with that comparison's k*. RIGHT is the selected top-k* block and LEFT the matched equal-sized block immediately below it. ",
      "Points give the scaling exponent p estimated within each block under E = alpha * mu^p by regressing log E on log mu, with 95% confidence intervals. ",
      "Dashed guides mark p = 1, the NB1-like case in which excess variance is proportional to the mean, and p = 2, the NB2-like case in which it is proportional to the squared mean. ",
      "Each row is annotated with LR, the likelihood ratio 2*(logLik with separate LEFT and RIGHT NB2 dispersions - logLik with one shared dispersion), and its p-value on one degree of freedom. ",
      "Rightward arrows with large LR show that the selected block scales more NB2-like than the block immediately below the cutoff."
    )
  )

  writeLines(lines, file.path(OUT_ROOT, "Figure_Legends.md"))
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

create_zip_from_files <- function(zip_path, files) {
  files <- unique(files[file.exists(files)])
  files <- files[normalizePath(files) != normalizePath(zip_path, mustWork = FALSE)]

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
    zip::zipr(
      zipfile = zip_path,
      files = rel
    )
  } else {
    utils::zip(
      zipfile = zip_path,
      files = rel,
      flags = "-r9X"
    )
  }

  setwd(old_wd)

  if (!file.exists(zip_path) || file.info(zip_path)$size <= 0) {
    stop("Zip creation failed: ", zip_path)
  }

  invisible(zip_path)
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



# =============================================================================
# RUN ANALYSIS
# =============================================================================

input <- read_count_data(
  COUNT_FILE,
  GROUP_PATTERNS
)

all_counts <- input$counts
annotation <- input$annotation

# Feature universe: zero filtering is applied once across all 40 RT/ZT samples,
# so N and every arm-specific rank axis are identical.
experiment_keep <- rowSums(all_counts) > 0
all_counts <- all_counts[experiment_keep, , drop = FALSE]
annotation <- annotation[
  match(rownames(all_counts), annotation$feature_id),
  ,
  drop = FALSE
]

N_EXPERIMENT <- nrow(all_counts)
if (N_EXPERIMENT < 20L) stop("Too few PASs after global nonzero filtering.")

experiment_group_labels <- assign_groups(
  colnames(all_counts),
  GROUP_PATTERNS
)

# Normalization and variance reference: DESeq2 size factors are estimated once
# across all 40 samples, and pooled within-group variance is then estimated
# across all eight arms.
experiment_deseq <- normalize_deseq2_comparison(
  count_mat = all_counts,
  group_labels = experiment_group_labels
)

experiment_normalized_counts <- experiment_deseq$normalized_counts
experiment_size_factors <- experiment_deseq$size_factors

experiment_pooled <- compute_pooled_within_group_variance(
  normalized_counts = experiment_normalized_counts,
  group_labels = experiment_group_labels
)

# Arm-specific PC1 rankings: PCA is run on log1p(CPM) expression within each
# arm and PASs are ranked by absolute PC1 loading, which is left untransformed.
experiment_group_results <- vector("list", length(levels(experiment_group_labels)))
names(experiment_group_results) <- levels(experiment_group_labels)

for (g in levels(experiment_group_labels)) {
  idx <- which(experiment_group_labels == g)
  if (length(idx) < 2L) stop("Not enough samples in group ", g)

  experiment_group_results[[g]] <- compute_group_analysis(
    group_name = g,
    raw_counts_arm = all_counts[, idx, drop = FALSE],
    normalized_counts_arm = experiment_normalized_counts[, idx, drop = FALSE],
    pooled_variance = experiment_pooled$variance
  )
}

# Regime boundaries: one c1/c2 fit is estimated jointly across all eight arms
# and held fixed while each comparison is scanned for its own k*.
shared_regime_fit <- fit_shared_knots(experiment_group_results)
SHARED_C1 <- shared_regime_fit$c1
SHARED_C2 <- shared_regime_fit$c2
SHARED_MAX_CANDIDATE_K <- as.integer(N_EXPERIMENT - SHARED_C2)

message("Experiment-wide feature count N = ", N_EXPERIMENT)
message("Shared eight-arm c1 = ", SHARED_C1)
message("Shared eight-arm c2 = ", SHARED_C2)
message("Shared regime-derived candidate ceiling k <= ", SHARED_MAX_CANDIDATE_K)

cutoff_rows <- list()
comparison_optima <- list()
evidence_rows <- list()
pareto_panels <- list()
expected_figure_paths <- character(0)
expected_table_paths <- character(0)

for (comparison_name in names(COMPARISONS)) {
  message("============================================================")
  message("Analyzing comparison: ", comparison_name)

  mapping <- COMPARISONS[[comparison_name]]
  control_group <- unname(mapping[["control"]])
  treatment_group <- unname(mapping[["treatment"]])

  control_idx <- which(experiment_group_labels == control_group)
  treatment_idx <- which(experiment_group_labels == treatment_group)

  if (length(control_idx) < 2L || length(treatment_idx) < 2L) {
    stop("Each comparison arm requires at least two samples: ", comparison_name)
  }

  annotation_cmp <- annotation
  N <- N_EXPERIMENT
  pooled <- experiment_pooled
  group_results <- experiment_group_results
  knot_fit <- shared_regime_fit
  c1 <- SHARED_C1
  c2 <- SHARED_C2

  raw_scan <- scan_pair_cutoffs(
    control_df = group_results[[control_group]]$data,
    treatment_df = group_results[[treatment_group]]$data,
    c1 = c1,
    c2 = c2,
    comparison_name = comparison_name,
    control_group = control_group,
    treatment_group = treatment_group
  )

  opt <- select_weighted_pareto_optimum(
    raw_scan,
    good_col = "good_n",
    cost_col = "remainder_cross_n",
    benefit_weight = BENEFIT_WEIGHT,
    contamination_weight = CONTAMINATION_WEIGHT
  )

  scan_df <- opt$scan
  selected_k <- as.integer(opt$selected_k)

  comparison_optima[[comparison_name]] <- opt

  if (!comparison_name %in% names(COMPARISONS)) {
    stop("Unexpected comparison while storing k*: ", comparison_name)
  }
  if (length(selected_k) != 1L || !is.finite(selected_k) || selected_k < 1L) {
    stop("Invalid comparison-specific k* for ", comparison_name)
  }

  cutoff_rank <- rank_cutoff_from_k(N, selected_k)

  class_df <- classify_pair_at_k(
    control_df = group_results[[control_group]]$data,
    treatment_df = group_results[[treatment_group]]$data,
    k = selected_k,
    c1 = c1,
    c2 = c2,
    comparison_name = comparison_name,
    control_group = control_group,
    treatment_group = treatment_group
  )

  # The Leading Edge is the union of the two arm-specific top-k* sets; the
  # Remainder is its complement. Opposite-arm Remainder crossings are a Pareto
  # cost term and are carried as a flag, not as an exclusion.
  lead_ids <- unique(class_df$feature_id)
  rem_ids <- setdiff(rownames(all_counts), lead_ids)

  control_rank_map <- make_rank_map(group_results[[control_group]]$data)
  treatment_rank_map <- make_rank_map(group_results[[treatment_group]]$data)

  lead_rows <- class_df %>%
    left_join(annotation_cmp, by = "feature_id") %>%
    mutate(
      evs_membership = "LeadingEdge",
      pareto_remainder_crossing_cost = cross_into_remainder
    ) %>%
    select(
      comparison,
      feature_id,
      gene_symbol,
      evs_membership,
      selected_k,
      cutoff_rank,
      control_group,
      treatment_group,
      site_class = analysis_class,
      control_rank,
      treatment_rank,
      control_region,
      treatment_region,
      control_top_k,
      treatment_top_k,
      pareto_remainder_crossing_cost
    )

  remainder_rows <- data.frame(
    comparison = comparison_name,
    feature_id = rem_ids,
    gene_symbol = annotation_cmp$gene_symbol[
      match(rem_ids, annotation_cmp$feature_id)
    ],
    evs_membership = "Remainder",
    selected_k = selected_k,
    cutoff_rank = cutoff_rank,
    control_group = control_group,
    treatment_group = treatment_group,
    site_class = "Remainder",
    control_rank = as.integer(control_rank_map[rem_ids]),
    treatment_rank = as.integer(treatment_rank_map[rem_ids]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      control_region = rank_to_region(control_rank, c1, c2),
      treatment_region = rank_to_region(treatment_rank, c1, c2),
      control_top_k = FALSE,
      treatment_top_k = FALSE,
      pareto_remainder_crossing_cost = FALSE
    )

  membership_table <- bind_rows(
    lead_rows,
    remainder_rows[, colnames(lead_rows), drop = FALSE]
  )

  if (nrow(membership_table) != N) {
    stop("EVS membership table does not cover every PAS in ", comparison_name)
  }

  arm_results <- list()

  for (g in c(control_group, treatment_group)) {
    idx_global <- which(experiment_group_labels == g)

    arm_results[[g]] <- build_arm_nb_evidence(
      arm_label = g,
      raw_counts_arm = all_counts[, idx_global, drop = FALSE],
      normalized_counts_arm = experiment_normalized_counts[, idx_global, drop = FALSE],
      size_factors = experiment_size_factors[colnames(all_counts)[idx_global]],
      rank_df = group_results[[g]]$data,
      k = selected_k
    )
  }

  # One row per arm at this comparison's k*: the NB scaling exponent in each
  # matched block and the likelihood-ratio test of the dispersion change.
  evidence_table <- bind_rows(arm_results) %>%
    mutate(
      comparison = comparison_name,
      selected_k = selected_k,
      cutoff_rank = cutoff_rank
    ) %>%
    select(
      comparison,
      selected_k,
      cutoff_rank,
      everything()
    )


  selected_scan_row <- scan_df %>%
    filter(k == selected_k) %>%
    slice(1L)

  cutoff_row <- data.frame(
    comparison = comparison_name,
    cutoff_scope = "comparison_specific_weighted_pareto",
    regime_boundary_scope = "shared_across_eight_arms",
    N = N,
    control_group = control_group,
    treatment_group = treatment_group,
    c1 = c1,
    c2 = c2,
    max_candidate_k = N - c2,
    selected_k = selected_k,
    cutoff_rank = cutoff_rank,
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
    pc1_energy_control_lead_fraction =
      sum(tail(group_results[[control_group]]$data$pc1_variance_contribution, selected_k)) /
      sum(group_results[[control_group]]$data$pc1_variance_contribution),
    pc1_energy_treatment_lead_fraction =
      sum(tail(group_results[[treatment_group]]$data$pc1_variance_contribution, selected_k)) /
      sum(group_results[[treatment_group]]$data$pc1_variance_contribution),
    pooled_within_group_residual_df = pooled$residual_df,
    knot_fit_SSE = knot_fit$SSE,
    stringsAsFactors = FALSE
  )

  comp_tab_dir <- file.path(OUT_ROOT, comparison_name, "Tables")
  dir.create(comp_tab_dir, recursive = TRUE, showWarnings = FALSE)

  pareto_panels[[comparison_name]] <- make_pareto_panel(
    scan_df = scan_df,
    selected_k = selected_k,
    comparison_name = comparison_name,
    max_candidate_k = SHARED_MAX_CANDIDATE_K
  )

  table_paths <- c(
    EVS_Membership = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_EVS_Membership.csv"))
  )

  write_csv(membership_table, table_paths[["EVS_Membership"]])

  verify_outputs(unname(table_paths))
  expected_table_paths <- c(expected_table_paths, unname(table_paths))

  cutoff_rows[[comparison_name]] <- cutoff_row
  evidence_rows[[comparison_name]] <- evidence_table

  message(
    comparison_name,
    ": k*=", selected_k,
    " | Lead union=", length(lead_ids),
    " | Remainder=", length(rem_ids),
    " | G=", selected_scan_row$good_n,
    " | R=", selected_scan_row$remainder_cross_n,
    " | p_left->p_right=",
    paste(
      formatC(evidence_table$p_exponent_left, format = "f", digits = 2),
      "->",
      formatC(evidence_table$p_exponent_right, format = "f", digits = 2),
      collapse = ", "
    )
  )
}

# =============================================================================
# CROSS-COMPARISON SUMMARY OUTPUTS
# =============================================================================

cutoff_summary <- bind_rows(cutoff_rows)

# One optimized k* per comparison is required before any summary is written.
expected_comparisons <- names(COMPARISONS)
observed_comparisons <- as.character(cutoff_summary$comparison)

if (nrow(cutoff_summary) != length(expected_comparisons)) {
  stop(
    "Expected ", length(expected_comparisons),
    " comparison-specific cutoff rows but found ", nrow(cutoff_summary), "."
  )
}

if (anyDuplicated(observed_comparisons)) {
  stop("Duplicate comparison rows detected in cutoff_summary.")
}

if (!setequal(observed_comparisons, expected_comparisons)) {
  stop(
    "Comparison-specific cutoff set mismatch. Expected: ",
    paste(expected_comparisons, collapse = ", "),
    "; observed: ", paste(observed_comparisons, collapse = ", ")
  )
}

if (!setequal(names(comparison_optima), expected_comparisons)) {
  stop("Not every comparison has its own stored weighted-Pareto optimum.")
}

comparison_specific_k <- stats::setNames(
  as.integer(cutoff_summary$selected_k),
  cutoff_summary$comparison
)[expected_comparisons]

message(
  "Independent comparison-specific empirical k*: ",
  paste(
    names(comparison_specific_k),
    comparison_specific_k,
    sep = "=",
    collapse = "; "
  )
)

evidence_all <- bind_rows(evidence_rows)

summary_cutoff_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_Empirical_Cutoffs.csv"
)

summary_evidence_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_NB_Evidence.csv"
)

write_csv(cutoff_summary, summary_cutoff_path)
write_csv(evidence_all, summary_evidence_path)

# =============================================================================
# MANUSCRIPT FIGURES
# =============================================================================
#
# Figure 1: (A) excess variance along the PC1 rank; (B) the two normalized
#           cumulative variance masses; (C) their divergence and the shared
#           two-transition model that fixes c1 and c2.
# Figure 2: the weighted Pareto objective across candidate k, one panel per
#           comparison, with the selected k* and its numbers.
# Figure 3: the NB scaling exponent in the matched blocks at each cutoff, with
#           the likelihood-ratio statistic and p-value.

fig1_path <- file.path(FIG_DIR, "Figure_1_Rank_Geometry.png")
fig2_path <- file.path(FIG_DIR, "Figure_2_Cutoff_Selection.png")
fig3_path <- file.path(FIG_DIR, "Figure_3_NB_Scaling_Corroboration.png")

save_stacked(
  make_divergence_figure(
    group_results = experiment_group_results,
    knot_fit = shared_regime_fit,
    c1 = SHARED_C1,
    c2 = SHARED_C2
  ),
  fig1_path,
  width = 13.0,
  height = 13.5
)

make_grid_figure(
  pareto_panels[names(COMPARISONS)],
  fig2_path
)

save_figure(
  make_evidence_figure(evidence_all),
  fig3_path,
  width = 13.0,
  height = 7.5
)

expected_figure_paths <- c(fig1_path, fig2_path, fig3_path)
if (isTRUE(EXPORT_PDF)) {
  expected_figure_paths <- c(
    expected_figure_paths,
    sub("\\.png$", ".pdf", fig1_path),
    sub("\\.png$", ".pdf", fig2_path),
    sub("\\.png$", ".pdf", fig3_path)
  )
}

expected_table_paths <- c(
  expected_table_paths,
  summary_cutoff_path,
  summary_evidence_path
)

write_methods_manuscript(cutoff_summary, evidence_all)
write_figure_legends()

verify_outputs(c(
  expected_figure_paths,
  expected_table_paths,
  file.path(OUT_ROOT, "Methods_Manuscript.md"),
  file.path(OUT_ROOT, "Figure_Legends.md")
))

# =============================================================================
# ZIP ARCHIVES
# =============================================================================

fig_zip <- file.path(OUT_ROOT, "Figures_All.zip")
all_zip <- file.path(OUT_ROOT, "Empirical_Cutoff_All_Outputs.zip")

create_zip_from_files(
  fig_zip,
  expected_figure_paths
)

create_zip_from_files(
  all_zip,
  c(
    expected_figure_paths,
    expected_table_paths,
    file.path(OUT_ROOT, "Methods_Manuscript.md"),
    file.path(OUT_ROOT, "Figure_Legends.md")
  )
)

verify_outputs(c(fig_zip, all_zip))

# =============================================================================
# FINAL CONSOLE SUMMARY
# =============================================================================

message("============================================================")
message("EMPIRICAL EVS CUTOFF ANALYSIS COMPLETE")
message("PCA input for ranking: arm-specific log1p(CPM) expression.")
message("One empirical k* is estimated independently for each comparison.")
message(
  "Selected k*: ",
  paste(
    cutoff_summary$comparison,
    cutoff_summary$selected_k,
    sep = "=",
    collapse = "; "
  )
)
message("Figures: 3 manuscript figures in ", FIG_DIR)
message("Figures ZIP: ", fig_zip)
message("Complete ZIP: ", all_zip)
message("Build: ", SCRIPT_BUILD)
message("Methods: ", file.path(OUT_ROOT, "Methods_Manuscript.md"))
message("Figure legends: ", file.path(OUT_ROOT, "Figure_Legends.md"))
message("============================================================")
