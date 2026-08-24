#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

SCRIPT_BUILD <- "EMPERICALCUTOFF_FINAL_COMPARISON_SPECIFIC_KSTAR_2026-08-23_v2"

# CRITICAL SCOPE RULE
# -------------------
# There is NO single global EVS k* in this analysis. The experiment-wide
# feature universe, normalization/variance reference, and shared c1/c2 regime
# model define a common rank geometry only. Each RT/ZT comparison then receives
# its own independent weighted-Pareto scan and its own selected k*.

# =============================================================================
# EMPIRICAL EVS CUTOFF: NORMEVS PC1-NB GEOMETRY + WEIGHTED PARETO
# =============================================================================
#
# PURPOSE
# -------
# Determine one empirical eigenvector-splitting cutoff k* independently for
# each comparison:
#
#   RT0_ZT6, RT2_ZT8, RT4_ZT10, RT8_ZT14.
#
# The cutoff-estimation core preserves the validated historical procedure.
# PASs are filtered once across the complete RT/ZT matrix. For PC1 ranking, raw
# counts are converted to CPM and log1p-transformed within each arm before PCA.
# PCA is centered and not feature-scaled; the loading scores themselves are NOT
# transformed, and PASs are ranked by raw absolute PC1 loading. DESeq2
# median-of-ratios normalization is estimated globally across all 40 samples for
# the pooled within-group NB variance term.
#
# CUTOFF CALCULATION
# ------------------
# For PAS i in arm g, the PC1 variance contribution is
#
#   P_ig = lambda_1g * loading_ig^2.
#
# Within each comparison, DESeq2-normalized counts are also used to estimate
# pooled within-group variance V_pool,i. For each arm,
#
#   E_ig = max(V_pool,i - mu_ig, 0)
#
# is the excess-over-Poisson variance. PC1 contribution P and excess variance E
# are normalized to rank-wise probability masses and accumulated along the
# absolute-PC1-loading rank:
#
#   D_g(r) = F_E,g(r) - F_P,g(r).
#
# A single shared two-knot continuous linear spline is fitted jointly to the
# D_g(r) curves from all eight RT/ZT arms. The fitted knots define:
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
# The Pareto frontier is formed by maximizing G while minimizing R. On the
# frontier, G and R are min-max normalized and the equal-weight utility is
#
#   U(k) = G_norm(k) - R_norm(k).
#
# The comparison-specific empirical cutoff is the Pareto-optimal k that
# maximizes U(k), with ties resolved by greater G, lower R, then larger k.
#
# FINAL EVS MEMBERSHIP
# --------------------
# The selected k* is applied independently to both arm-specific absolute-PC1
# loading rankings. The final Leading Edge is the union of the two top-k* sets;
# Joint and all Disjoint PASs remain members of the Leading Edge. Remainder is
# the complement of that union. Opposite-arm Remainder crossings are used only
# as the Pareto cost R(k); they do not override union membership.
#
# NB1/NB2 CORROBORATION
# ---------------------
# After k* is locked, each arm's selected top-k* region (RIGHT) is compared with
# an immediately preceding equal-sized rank block (LEFT). Moment-based
# corroboration is calculated from raw count moments:
#
#   NB2 excess signal = log[1 + max(variance - mean, 0)]
#   NB2-NB1 contrast  = NB2 excess signal - log(1 + mean)
#   alpha_hat         = max[(variance - mean)/mean^2, 0]
#   alpha*mu signal   = log(1 + alpha_hat*mean).
#
# Likelihood-based corroboration is calculated from raw counts with expected
# means adjusted by the comparison-specific DESeq2 size factors. NB1 and NB2
# models each fit one dispersion parameter by maximum likelihood:
#
#   NB1: Var(Y) = mu + alpha*mu
#   NB2: Var(Y) = mu + alpha*mu^2.
#
# For the negative-binomial PMF, NB1 uses size = mu/alpha and NB2 uses
# size = 1/alpha. The exported likelihood evidence includes log-likelihoods,
# 2*(logLik_NB2-logLik_NB1) and log10(L_NB2/L_NB1).
# Positive 2-delta-log-likelihood and positive log10 likelihood ratio favor NB2.
# These corroboration quantities are calculated only after k* is selected and
# do not enter the Pareto optimization.
#
# OUTPUT ORGANIZATION
# -------------------
# Each comparison receives:
#
#   <OUT_ROOT>/<comparison>/Figures/
#   <OUT_ROOT>/<comparison>/Tables/
#
# with two manuscript figures and all supporting tables. Root-level Summary/
# folders contain cross-comparison summaries. Every generated figure is
# included in Figures_All.zip, every generated CSV table in Tables_All.zip, and
# the complete output tree in Empirical_Cutoff_All_Outputs.zip.
#
# Methods_Manuscript.md is generated from the same constants and equations used
# by the code so the methods description and implementation remain aligned.
# =============================================================================

# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT   <- "/root/REAPER98632/exports/empirical_cutoff_comparison_specific_kstar_final_20260823_v2"

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

# Validated regression targets from the pre-regression WTTS cutoff run.
# These values are NOT used to select k*. They are checked only after k* has
# been recomputed from the data, so an unintended implementation change fails
# loudly instead of silently replacing the validated cutoffs.
VALIDATED_REFERENCE_K <- c(
  RT0_ZT6  = 3532L,
  RT2_ZT8  = 4617L,
  RT4_ZT10 = 3983L,
  RT8_ZT14 = 5664L
)
ENFORCE_VALIDATED_K_REGRESSION <- TRUE


# Numerical optimization uses a coarse starting grid and then refits the
# selected knot solution on every rank. This value controls computation only;
# the final c1/c2 fit is evaluated on the full rank series.
KNOT_COARSE_GRID_POINTS <- 5000L

# Display-only smoothing. These values do not determine k*.
DISPLAY_VAR_SPAR <- 0.72
DISPLAY_D_SPAR   <- 0.72
DISPLAY_NB_SPAR  <- 0.68

# Figure export.
PNG_DPI <- 360
EXPORT_PDF <- TRUE

# Likelihood fitting bounds on log(alpha).
NB_LOG_ALPHA_LOWER <- -14
NB_LOG_ALPHA_UPPER <- 8


dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
SUMMARY_FIG_DIR <- file.path(OUT_ROOT, "Summary", "Figures")
SUMMARY_TAB_DIR <- file.path(OUT_ROOT, "Summary", "Tables")
dir.create(SUMMARY_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# COLORS / FIGURE THEME
# =============================================================================

COL <- list(
  control = "#386CB0",
  treatment = "#159D91",
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

  raw_var <- apply(
    raw_counts_arm,
    1L,
    stats::var,
    na.rm = TRUE
  )

  raw_var[!is.finite(raw_var)] <- 0
  raw_var <- pmax(raw_var, 0)

  raw_var_ranked <- raw_var[rank_order]

  rank <- seq_along(rank_order)
  log_var <- log1p(raw_var_ranked)

  # Display-only smoothing. This curve is descriptive and is not used to
  # determine c1, c2, or k*.
  display_spline <- stats::smooth.spline(
    x = rank,
    y = log_var,
    spar = DISPLAY_VAR_SPAR
  )

  display_y <- as.numeric(
    stats::predict(
      display_spline,
      x = rank,
      deriv = 0
    )$y
  )

  data.frame(
    rank = rank,
    raw_empirical_variance = raw_var_ranked,
    log1p_raw_empirical_variance = log_var,
    display_log1p_raw_empirical_variance = display_y,
    stringsAsFactors = FALSE
  )
}

smooth_divergence_for_display <- function(
    rank,
    D,
    spar = DISPLAY_D_SPAR) {

  fit <- stats::smooth.spline(
    x = rank,
    y = D,
    spar = spar
  )

  y <- as.numeric(
    stats::predict(
      fit,
      x = rank,
      deriv = 0
    )$y
  )

  endpoint_line <- seq(
    y[1L],
    y[length(y)],
    length.out = length(y)
  )

  y - endpoint_line
}


compute_group_analysis <- function(
    group_name,
    raw_counts_arm,
    normalized_counts_arm,
    pooled_variance) {

  # RESTORED CORRECT RANKING CORE:
  # Library-size normalize raw counts to CPM, apply log1p to the EXPRESSION
  # matrix, then run PCA. The PC1 loading values themselves are NOT logged or
  # otherwise transformed; PASs are ranked by the raw absolute PC1 loading.
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

  display_D <- smooth_divergence_for_display(
    rank = rank,
    D = D
  )

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
      cumulative_divergence = D,
      display_D = display_D
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
  # These sites are intentionally PERMISSIBLE.
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

summarize_classification <- function(
    class_df) {

  if (
    nrow(class_df) < 1L
  ) {
    stop(
      "Cannot summarize empty classification."
    )
  }

  joint_n <- sum(
    class_df$base_class ==
    "Joint"
  )

  control_only <- (
    class_df$control_top_k &
    !class_df$treatment_top_k
  )

  treatment_only <- (
    class_df$treatment_top_k &
    !class_df$control_top_k
  )

  disjoint_control_opposite_le_n <- sum(
    control_only &
    class_df$disjoint_opposite_leading_edge
  )

  disjoint_treatment_opposite_le_n <- sum(
    treatment_only &
    class_df$disjoint_opposite_leading_edge
  )

  disjoint_control_opposite_divergence_n <- sum(
    control_only &
    class_df$disjoint_opposite_divergence
  )

  disjoint_treatment_opposite_divergence_n <- sum(
    treatment_only &
    class_df$disjoint_opposite_divergence
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

  good_n <- (
    joint_n +
    permissible_disjoint_n
  )

  remainder_cross_n <- sum(
    class_df$cross_into_remainder
  )

  union_n <- nrow(
    class_df
  )

  if (
    good_n +
    remainder_cross_n !=
    union_n
  ) {
    stop(
      "Classification counts do not sum to union size."
    )
  }

  data.frame(
    k = unique(
      class_df$selected_k
    )[1L],

    cutoff_rank = unique(
      class_df$cutoff_rank
    )[1L],

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

    remainder_cross_n =
      remainder_cross_n,

    union_n = union_n,

    retained_fraction = (
      good_n /
      union_n
    ),

    remainder_cross_fraction = (
      remainder_cross_n /
      union_n
    ),

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

  # NB corroboration follows the raw-count moment formulation.
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

  if (left_start < 1L) {
    stop(
      "Matched LEFT block is unavailable for k=", k,
      "; N=", N, ". Reduce k or revise the comparison domain."
    )
  }

  rank_df$corroboration_region <- "Other"
  rank_df$corroboration_region[
    rank_df$rank >= left_start & rank_df$rank <= left_end
  ] <- "LEFT"
  rank_df$corroboration_region[
    rank_df$rank >= right_start & rank_df$rank <= right_end
  ] <- "RIGHT"

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
  out[ok] <- as.numeric(stats::predict(fit, x = rank[ok], deriv = 0)$y)
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
  two_delta_ll <- 2 * delta_ll
  log10_lr <- delta_ll / log(10)

  data.frame(
    arm = arm_label,
    region = region_label,
    n_features = length(feature_ids),
    n_observations = nb1$n_obs,
    alpha_NB1 = nb1$alpha,
    alpha_NB2 = nb2$alpha,
    logLik_NB1 = nb1$logLik,
    logLik_NB2 = nb2$logLik,
    two_delta_logLik_NB2_minus_NB1 = two_delta_ll,
    log10_likelihood_ratio_NB2_to_NB1 = log10_lr,
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
    mutate(arm = arm_label)

  lrt_rows <- list()

  for (region_name in c("LEFT", "RIGHT")) {
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
    likelihood = bind_rows(lrt_rows)
  )
}

# =============================================================================
# FIGURE BUILDERS
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

make_pareto_panel <- function(scan_df, selected_k, comparison_name) {
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
      title = paste0("D. ", comparison_name, ": weighted Pareto cutoff"),
      subtitle = "Equal normalized weights: U(k) = Gnorm(k) - Rnorm(k)",
      x = "Opposite-arm Remainder crossings, R(k)  [cost]",
      y = "Joint + permissible Disjoint, G(k)  [benefit]",
      caption = "The selected k* maximizes weighted utility on the Pareto frontier; opposite-arm Divergence crossings are permissible."
    ) +
    theme_manuscript() +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))
}

make_cutoff_framework_figure <- function(
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
    scale_color_manual(
      name = NULL,
      values = c(arm_colors, boundary_colors)
    ) +
    scale_linetype_manual(
      name = NULL,
      values = boundary_types,
      na.translate = FALSE
    ) +
    labs(
      title = paste0("A. ", comparison_name, ": absolute PC1 loading rank"),
      subtitle = "log1p(CPM) expression entered into PCA; untransformed absolute PC1 loadings ranked",
      x = "PC1 rank: low |loading| to high |loading|",
      y = "|PC1 loading|",
      caption = "c1 = Remainder/Divergence boundary; c2 = Divergence/Leading-edge boundary; k* = weighted-Pareto top-k cutoff."
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
    scale_color_manual(
      name = NULL,
      values = c(arm_colors, boundary_colors)
    ) +
    scale_linetype_manual(
      name = NULL,
      values = boundary_types,
      na.translate = FALSE
    ) +
    labs(
      title = paste0("B. ", comparison_name, ": raw-count variance geometry"),
      subtitle = "Raw-count sample variance shown on a log(1+x) display scale along the same PC1 rank",
      x = "PC1 rank",
      y = "Smoothed log(1 + raw-count variance)",
      caption = "Region shading follows c1 and c2; the selected k* is shown as the top-k boundary on the ranked axis."
    ) +
    theme_manuscript() +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 1)
    )

  fit_mat <- knot_fit$fitted
  pair_fit_idx <- match(c(control_group, treatment_group), knot_fit$groups)
  if (anyNA(pair_fit_idx)) {
    stop("Pair arms were not found in the shared eight-arm knot fit.")
  }
  fit_long <- bind_rows(lapply(pair_fit_idx, function(j) {
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
    scale_color_manual(
      name = NULL,
      values = c(arm_colors, boundary_colors)
    ) +
    scale_linetype_manual(
      name = NULL,
      values = boundary_types,
      na.translate = FALSE
    ) +
    labs(
      title = paste0("C. ", comparison_name, ": PC1-NB cumulative divergence"),
      subtitle = "One two-knot fit shared across all eight RT/ZT arms; pair shown here",
      x = "PC1 rank",
      y = "D(r) = F_E(r) - F_P(r)",
      caption = "F_E is cumulative excess-over-Poisson variance mass; F_P is cumulative PC1 variance-contribution mass."
    ) +
    theme_manuscript() +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 1)
    )

  pD <- make_pareto_panel(
    scan_df = scan_df,
    selected_k = selected_k,
    comparison_name = comparison_name
  )

  save_grid_2x2(
    list(pA, pB, pC, pD),
    out_file,
    width = 16,
    height = 12.5
  )
}

make_nb_corroboration_figure <- function(
    comparison_name,
    control_group,
    treatment_group,
    selected_k,
    arm_results,
    out_file) {

  rank_data <- bind_rows(lapply(names(arm_results), function(a) {
    arm_results[[a]]$rank_data %>% mutate(arm = a)
  }))

  region_summary <- bind_rows(lapply(
    names(arm_results),
    function(a) arm_results[[a]]$region_summary
  ))

  likelihood <- bind_rows(lapply(
    names(arm_results),
    function(a) arm_results[[a]]$likelihood
  ))

  N <- max(rank_data$rank)
  cutoff_rank <- rank_cutoff_from_k(N, selected_k)
  left_start <- cutoff_rank - selected_k
  left_end <- cutoff_rank - 1L

  metric_long <- bind_rows(
    rank_data %>%
      transmute(
        rank, arm,
        metric = "NB2 excess",
        value = nb2_excess_signal_smooth
      ),
    rank_data %>%
      transmute(
        rank, arm,
        metric = "NB2-NB1",
        value = nb2_nb1_contrast_smooth
      ),
    rank_data %>%
      transmute(
        rank, arm,
        metric = "alpha*mu",
        value = alpha_mu_signal_smooth
      )
  )

  metric_colors <- c(
    "NB2 excess" = COL$nb2,
    "NB2-NB1" = COL$nbgap,
    "alpha*mu" = COL$alphamu
  )

  pA <- ggplot(metric_long, aes(rank, value, color = metric)) +
    annotate(
      "rect",
      xmin = left_start, xmax = left_end,
      ymin = -Inf, ymax = Inf,
      fill = COL$left, alpha = 0.12
    ) +
    annotate(
      "rect",
      xmin = cutoff_rank, xmax = N,
      ymin = -Inf, ymax = Inf,
      fill = COL$right, alpha = 0.12
    ) +
    geom_vline(
      xintercept = cutoff_rank,
      color = COL$selected,
      linetype = "dotdash",
      linewidth = 0.9
    ) +
    geom_line(linewidth = 0.85) +
    facet_wrap(~ arm, ncol = 1, scales = "free_y") +
    scale_color_manual(
      name = "Corroboration metric",
      values = metric_colors
    ) +
    labs(
      title = paste0("A. ", comparison_name, ": NB1/NB2 corroboration along EVS rank"),
      subtitle = paste0(
        "RIGHT = selected top-k* (k*=", selected_k,
        "); LEFT = immediately preceding matched block"
      ),
      x = "Absolute-PC1-loading rank",
      y = "Smoothed corroboration signal"
    ) +
    theme_manuscript() +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))

  # Mean-variance relationship with region/model fits.
  mv <- rank_data %>%
    filter(corroboration_region %in% c("LEFT", "RIGHT")) %>%
    select(
      arm,
      corroboration_region,
      raw_mean,
      raw_variance
    ) %>%
    rename(region = corroboration_region)

  # Deterministic thinning for readability only; model fits use all features.
  mv_plot <- mv %>%
    group_by(arm, region) %>%
    arrange(raw_mean) %>%
    mutate(plot_keep = row_number() %% max(1L, ceiling(n() / 1200L)) == 0L) %>%
    filter(plot_keep) %>%
    ungroup()

  curve_rows <- list()
  idx <- 1L

  for (i in seq_len(nrow(likelihood))) {
    rr <- likelihood[i, , drop = FALSE]
    arm_i <- rr$arm
    region_i <- rr$region

    mu_vals <- mv$raw_mean[
      mv$arm == arm_i & mv$region == region_i
    ]
    mu_vals <- mu_vals[is.finite(mu_vals) & mu_vals > 0]
    if (!length(mu_vals)) next

    grid_mu <- exp(seq(
      log(max(min(mu_vals), 1e-6)),
      log(max(mu_vals)),
      length.out = 180L
    ))

    curve_rows[[idx]] <- data.frame(
      arm = arm_i,
      region = region_i,
      raw_mean = grid_mu,
      predicted_variance = grid_mu + rr$alpha_NB1 * grid_mu,
      model = "NB1",
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L

    curve_rows[[idx]] <- data.frame(
      arm = arm_i,
      region = region_i,
      raw_mean = grid_mu,
      predicted_variance = grid_mu + rr$alpha_NB2 * grid_mu^2,
      model = "NB2",
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }

  curve_df <- bind_rows(curve_rows)

  mv_plot <- mv_plot %>%
    mutate(
      display_log1p_mean = log1p(raw_mean),
      display_log1p_variance = log1p(raw_variance)
    )

  curve_df <- curve_df %>%
    mutate(
      display_log1p_mean = log1p(raw_mean),
      display_log1p_variance = log1p(predicted_variance)
    )

  pB <- ggplot(
    mv_plot,
    aes(display_log1p_mean, display_log1p_variance)
  ) +
    geom_point(
      size = 0.8,
      alpha = 0.20,
      color = "grey35"
    ) +
    geom_line(
      data = curve_df,
      aes(
        display_log1p_mean,
        display_log1p_variance,
        color = model
      ),
      linewidth = 1.0
    ) +
    facet_grid(arm ~ region, scales = "free") +
    scale_color_manual(
      name = "Fitted variance model",
      values = c("NB1" = COL$nb1, "NB2" = COL$nb2fit)
    ) +
    labs(
      title = paste0("B. ", comparison_name, ": observed mean-variance relationship"),
      subtitle = "Points are PAS-level raw-count moments; lines use maximum-likelihood NB1/NB2 dispersion estimates",
      x = "log(1 + raw-count mean)",
      y = "log(1 + raw-count variance)"
    ) +
    theme_manuscript() +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))

  summary_long <- region_summary %>%
    pivot_longer(
      cols = c(
        median_nb2_excess,
        median_nb2_nb1,
        median_alpha_mu
      ),
      names_to = "metric",
      values_to = "median_value"
    ) %>%
    mutate(
      metric = recode(
        metric,
        median_nb2_excess = "NB2 excess",
        median_nb2_nb1 = "NB2-NB1",
        median_alpha_mu = "alpha*mu"
      ),
      region = factor(
        corroboration_region,
        levels = c("LEFT", "RIGHT")
      )
    )

  pC <- ggplot(
    summary_long,
    aes(
      x = median_value,
      y = metric,
      color = region
    )
  ) +
    geom_point(size = 3.3) +
    facet_wrap(~ arm, ncol = 1, scales = "free_x") +
    scale_color_manual(
      name = "Matched rank block",
      values = c("LEFT" = COL$left, "RIGHT" = COL$right)
    ) +
    labs(
      title = paste0("C. ", comparison_name, ": matched LEFT versus RIGHT medians"),
      subtitle = "The Pareto-selected cutoff is fixed before these summaries are calculated",
      x = "Median corroboration signal",
      y = NULL
    ) +
    theme_manuscript() +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))

  likelihood$arm_region <- paste(likelihood$arm, likelihood$region, sep = " / ")
  likelihood$arm_region <- factor(
    likelihood$arm_region,
    levels = rev(likelihood$arm_region)
  )

  pD <- ggplot(
    likelihood,
    aes(
      x = two_delta_logLik_NB2_minus_NB1,
      y = arm_region,
      fill = preferred_model
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.55
    ) +
    geom_col(width = 0.68) +
    geom_text(
      aes(
        label = paste0(
          "log10 LR = ",
          formatC(
            log10_likelihood_ratio_NB2_to_NB1,
            format = "f",
            digits = 2
          )
        )
      ),
      hjust = ifelse(
        likelihood$two_delta_logLik_NB2_minus_NB1 >= 0,
        -0.05,
        1.05
      ),
      size = 3.0
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.24))) +
    scale_fill_manual(
      name = "Higher likelihood",
      values = c(
        "NB1" = COL$nb1,
        "NB2" = COL$nb2fit,
        "Tie" = "grey60"
      )
    ) +
    labs(
      title = paste0("D. ", comparison_name, ": NB2 versus NB1 likelihood evidence"),
      subtitle = "Positive 2ΔlogL and positive log10(L_NB2/L_NB1) favor NB2",
      x = "2 × (logLik_NB2 - logLik_NB1)",
      y = NULL
    ) +
    theme_manuscript() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))

  save_grid_2x2(
    list(pA, pB, pC, pD),
    out_file,
    width = 16,
    height = 12.5
  )
}

make_summary_figure <- function(cutoff_summary, out_file) {
  df <- cutoff_summary %>%
    mutate(
      comparison = factor(comparison, levels = comparison),
      label = paste0(
        "k*=", selected_k,
        "\nG=", good_n,
        ", R=", remainder_cross_n
      )
    )

  p <- ggplot(
    df,
    aes(comparison, selected_k)
  ) +
    geom_col(
      fill = COL$selected,
      width = 0.65
    ) +
    geom_text(
      aes(label = label),
      vjust = -0.35,
      size = 3.5,
      fontface = "bold"
    ) +
    labs(
      title = "Comparison-specific empirical EVS cutoffs",
      subtitle = "Normalized-before-EVS PC1-NB geometry followed by equal-weight Pareto optimization",
      x = NULL,
      y = "Selected top-k* per arm"
    ) +
    theme_manuscript(base_size = 12.5) +
    theme(legend.position = "none") +
    expand_limits(
      y = max(df$selected_k, na.rm = TRUE) * 1.15
    )

  save_figure(
    p,
    out_file,
    width = 11.5,
    height = 7.2
  )
}

# =============================================================================
# METHODS / EXPORTS
# =============================================================================

write_methods_manuscript <- function(cutoff_summary) {
  cutoff_text <- paste(
    paste0(
      cutoff_summary$comparison,
      " k*=",
      cutoff_summary$selected_k
    ),
    collapse = ", "
  )

  lines <- c(
    "# Empirical eigenvector-splitting cutoff methods",
    "",
    "## Comparison-specific preprocessing and PC1 ranking",
    "",
    paste0(
      "Empirical EVS cutoffs were estimated independently for RT0_ZT6, RT2_ZT8, RT4_ZT10, and RT8_ZT14; no single global EVS k* was estimated or applied across comparisons. ",
      "PASs with zero counts across the complete 40-sample RT/ZT matrix were removed once before cutoff estimation. ",
      "For PC1 ranking, raw counts were library-size normalized to counts per million (CPM) within each arm and log1p-transformed before PCA. ",
      "PCA was performed independently within each arm with centering and without feature scaling. The PC1 loading values themselves were not transformed; PASs were ordered from lowest to highest raw absolute PC1 loading."
    ),
    "",
    "## PC1-NB variance-mass divergence",
    "",
    paste0(
      "For PAS i in arm g, PC1 variance contribution was P_ig = lambda_1g * loading_ig^2, where lambda_1g is the PC1 eigenvalue. ",
      "DESeq2 median-of-ratios size factors were estimated once across all 40 samples, and pooled within-group variance was calculated from those normalized counts across all eight RT/ZT arms. ",
      "For each arm, excess-over-Poisson variance was E_ig = max(V_pool,i - mu_ig, 0). ",
      "P and E were normalized separately to rank-wise probability masses and accumulated along the absolute-PC1-loading rank. ",
      "Cumulative divergence was D_g(r) = F_E,g(r) - F_P,g(r)."
    ),
    "",
    "## Comparison-specific variance regimes",
    "",
    paste0(
      "A single shared two-knot continuous linear spline was fitted jointly to the cumulative-divergence curves from all eight RT/ZT arms. ",
      "The shared fitted knots defined c1 and c2 for every comparison. Ranks below c1 were classified as the Remainder regime, ranks from c1 through c2 as the Divergence interval, and ranks above c2 as the Leading-edge regime."
    ),
    "",
    "## Weighted Pareto cutoff",
    "",
    paste0(
      "For each RT/ZT pair, candidate top-k values were restricted to 1 <= k <= N-c2 using the globally shared c2 so each selecting arm contributed PASs from its own Leading-edge regime. ",
      "For each k, the benefit G(k) was the number of Joint PASs plus Disjoint PASs whose opposite-arm rank lay in either the Leading-edge regime or Divergence interval. ",
      "The cost R(k) was the number of Disjoint PASs whose opposite-arm rank lay in the Remainder regime. ",
      "Pareto-optimal candidates were those for which no other candidate simultaneously increased G and decreased R. ",
      "Within the Pareto frontier, G and R were min-max normalized and equal weights were used: U(k) = G_norm(k) - R_norm(k). ",
      "The empirical cutoff k* maximized U(k); ties were resolved by greater G, then lower R, then larger k."
    ),
    "",
    "## EVS membership",
    "",
    paste0(
      "After k* was selected, the k* highest absolute-PC1-loading PASs were selected independently in the two arms. ",
      "The Leading Edge was the union of the two top-k* sets, including Joint and all Disjoint PASs, and the Remainder was the complement. ",
      "Opposite-arm Remainder crossings contributed to the Pareto cost but did not override union membership."
    ),
    "",
    "## NB1/NB2 corroboration",
    "",
    paste0(
      "After k* was fixed, each arm's selected top-k* block was designated RIGHT and compared with the immediately preceding equal-sized rank block, designated LEFT. ",
      "From raw counts, PAS-level sample mean and variance were used to calculate NB2 excess signal = log[1+max(variance-mean,0)], ",
      "NB2-NB1 contrast = NB2 excess signal - log(1+mean), and alpha*mu signal = log(1+alpha_hat*mean), ",
      "where alpha_hat = max[(variance-mean)/mean^2,0]."
    ),
    "",
    paste0(
      "Likelihood-based NB1/NB2 corroboration used the raw counts and comparison-specific DESeq2 size factors. ",
      "Expected raw-count means were obtained by multiplying each PAS's arm-specific normalized mean by the sample size factor. ",
      "NB1 and NB2 models each fitted one dispersion parameter alpha by maximum likelihood. ",
      "NB1 used Var(Y)=mu+alpha*mu and negative-binomial size=mu/alpha; NB2 used Var(Y)=mu+alpha*mu^2 and size=1/alpha. ",
      "Model support was summarized by the fitted log-likelihoods, 2*(logLik_NB2-logLik_NB1), and log10(L_NB2/L_NB1). ",
      "The NB corroboration analysis was performed only after cutoff selection and did not contribute to k*."
    ),
    "",
    "## Selected empirical cutoffs",
    "",
    cutoff_text
  )

  writeLines(
    lines,
    file.path(OUT_ROOT, "Methods_Manuscript.md")
  )
}

write_figure_legends <- function() {
  lines <- c(
    "# Figure legends",
    "",
    "## Figure 1. Comparison-specific normalized-EVS cutoff framework",
    "For each RT/ZT comparison, Panel A shows PASs ordered from low to high raw absolute PC1 loading after PCA was applied to arm-specific log1p(CPM) expression; the loading values themselves are not transformed. Panel B shows the corresponding raw-count variance geometry on the same PC1-ranked axis; this smoothed display curve is descriptive. Panel C shows the cumulative divergence D(r)=F_E(r)-F_P(r) between normalized excess-over-Poisson variance mass and PC1 variance-contribution mass, with the single c1 and c2 regime boundaries shared across all eight arms. Panel D shows the complete weighted-Pareto candidate set, the Pareto frontier, and the selected comparison-specific k*. Benefit G(k) retains Joint PASs and permissible Disjoint PASs; cost R(k) counts Disjoint PASs whose opposite-arm rank lies in the Remainder regime. The selected k* maximizes U(k)=G_norm(k)-R_norm(k) on the Pareto frontier.",
    "",
    "## Figure 2. NB1/NB2 and likelihood corroboration",
    "For each RT/ZT comparison, the independently selected weighted-Pareto k* is fixed before corroboration. Panel A shows the NB2 excess, NB2-NB1, and alpha*mu signals along the absolute-PC1-loading rank for the matched LEFT block and selected RIGHT block. Panel B shows PAS-level raw-count mean-variance observations with maximum-likelihood NB1 and NB2 variance-model curves. Panel C compares the matched LEFT and RIGHT median corroboration signals. Panel D shows 2*(logLik_NB2-logLik_NB1) and labels log10(L_NB2/L_NB1) for each arm and rank block; positive values favor NB2. These corroboration quantities do not enter the Pareto optimization or change k*.",
    "",
    "## Summary figure. Comparison-specific empirical EVS cutoffs",
    "The summary figure displays the independently selected top-k* value for RT0_ZT6, RT2_ZT8, RT4_ZT10, and RT8_ZT14 together with the retained-benefit G(k*) and Remainder-crossing cost R(k*) at each selected weighted-Pareto optimum. No global k* is substituted for these four comparison-specific optima."
  )

  writeLines(
    lines,
    file.path(OUT_ROOT, "Figure_Legends.md")
  )
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

# EXPERIMENT-WIDE FEATURE UNIVERSE (NOT A GLOBAL k*):
# zero filtering occurs once across all 40 RT/ZT samples, not separately inside
# each comparison. This keeps N and every arm-specific rank axis identical.
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

# EXPERIMENT-WIDE NORMALIZATION/VARIANCE REFERENCE (NOT A GLOBAL k*):
# size factors are estimated once across all 40 samples; pooled within-group
# variance is then estimated across all eight arms.
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

# ARM-SPECIFIC PC1 RANKINGS:
# each arm uses log1p(CPM) expression for PCA; the absolute PC1 loadings are
# ranked directly without transforming the loading scores themselves.
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

# SHARED EIGHT-ARM REGIME MODEL (BOUNDARIES ONLY; k* REMAINS COMPARISON-SPECIFIC):
# one c1/c2 fit is estimated jointly across all eight arms and is then held
# fixed while each of the four RT/ZT pairs receives its own Pareto-optimal k*.
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
corroboration_rows <- list()
selected_site_rows <- list()
remainder_site_rows <- list()
pareto_cost_rows <- list()
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

  sample_idx <- c(control_idx, treatment_idx)
  count_mat <- all_counts[, sample_idx, drop = FALSE]
  normalized_counts <- experiment_normalized_counts[, sample_idx, drop = FALSE]
  size_factors <- experiment_size_factors[colnames(count_mat)]

  group_labels <- factor(
    as.character(experiment_group_labels[sample_idx]),
    levels = c(control_group, treatment_group)
  )
  names(group_labels) <- colnames(count_mat)

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

  # k* is owned by THIS comparison only. It is never pooled across comparisons
  # and is never replaced by an experiment-wide/global k.
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

  class_summary <- summarize_classification(class_df)

  # Preserve the current requested EVS membership rule: the final Leading Edge
  # is the full union of the two arm-specific top-k sets. Remainder crossings
  # remain a Pareto cost term and do not override union membership.
  lead_ids <- unique(class_df$feature_id)
  rem_ids <- setdiff(rownames(all_counts), lead_ids)

  lead_table <- class_df %>%
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

  # NB1/NB2 corroboration remains strictly post-selection and cannot change k*.
  arm_results <- list()

  for (g in c(control_group, treatment_group)) {
    idx_global <- which(experiment_group_labels == g)

    arm_results[[g]] <- build_corroboration_for_arm(
      arm_label = g,
      raw_counts_arm = all_counts[, idx_global, drop = FALSE],
      normalized_counts_arm = experiment_normalized_counts[, idx_global, drop = FALSE],
      size_factors = experiment_size_factors[colnames(all_counts)[idx_global]],
      rank_df = group_results[[g]]$data,
      k = selected_k
    )
  }

  likelihood_table <- bind_rows(lapply(
    names(arm_results),
    function(g) arm_results[[g]]$likelihood
  )) %>%
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

  region_summary <- bind_rows(lapply(
    names(arm_results),
    function(g) arm_results[[g]]$region_summary
  )) %>%
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

  comp_fig_dir <- file.path(OUT_ROOT, comparison_name, "Figures")
  comp_tab_dir <- file.path(OUT_ROOT, comparison_name, "Tables")
  dir.create(comp_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(comp_tab_dir, recursive = TRUE, showWarnings = FALSE)

  cutoff_fig <- file.path(
    comp_fig_dir,
    paste0("Figure_", comparison_name, "_Cutoff_Framework.png")
  )

  nb_fig <- file.path(
    comp_fig_dir,
    paste0("Figure_", comparison_name, "_NB_Corroboration.png")
  )

  make_cutoff_framework_figure(
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

  make_nb_corroboration_figure(
    comparison_name = comparison_name,
    control_group = control_group,
    treatment_group = treatment_group,
    selected_k = selected_k,
    arm_results = arm_results,
    out_file = nb_fig
  )

  table_paths <- c(
    Cutoff_Summary = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_Cutoff_Summary.csv")),
    Cutoff_Scan = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_Weighted_Pareto_Scan.csv")),
    Leading_Edge = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_Leading_Edge_Sites.csv")),
    Remainder = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_Remainder_Sites.csv")),
    Pareto_Cost = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_Pareto_Remainder_Crossing_Cost_Sites.csv")),
    PC1_NB_Rank = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_PC1_NB_Rank_Data.csv")),
    NB_Region_Summary = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_NB_Corroboration_Regions.csv")),
    NB_Likelihood = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_NB1_NB2_Likelihood.csv")),
    Size_Factors = file.path(comp_tab_dir, paste0("Table_", comparison_name, "_DESeq2_Size_Factors.csv"))
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
      group = as.character(experiment_group_labels[names(size_factors)]),
      size_factor = as.numeric(size_factors),
      stringsAsFactors = FALSE
    ),
    table_paths[["Size_Factors"]]
  )

  expected_figs_this <- c(cutoff_fig, nb_fig)
  if (isTRUE(EXPORT_PDF)) {
    expected_figs_this <- c(
      expected_figs_this,
      sub("\\.png$", ".pdf", cutoff_fig),
      sub("\\.png$", ".pdf", nb_fig)
    )
  }

  verify_outputs(c(expected_figs_this, unname(table_paths)))

  expected_figure_paths <- c(expected_figure_paths, expected_figs_this)
  expected_table_paths <- c(expected_table_paths, unname(table_paths))

  cutoff_rows[[comparison_name]] <- cutoff_row
  corroboration_rows[[comparison_name]] <- likelihood_table
  selected_site_rows[[comparison_name]] <- lead_table
  remainder_site_rows[[comparison_name]] <- remainder_table
  pareto_cost_rows[[comparison_name]] <- cost_table

  message(
    comparison_name,
    ": comparison-specific k*=", selected_k,
    " | shared c1=", c1,
    " | shared c2=", c2,
    " | Lead union=", length(lead_ids),
    " | Remainder=", length(rem_ids),
    " | G=", selected_scan_row$good_n,
    " | R=", selected_scan_row$remainder_cross_n
  )
}

# =============================================================================
# CROSS-COMPARISON SUMMARY OUTPUTS
# =============================================================================

cutoff_summary <- bind_rows(cutoff_rows)

# HARD SCOPE CHECKS: one and only one independently optimized k* per comparison.
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

# Regression guard: compare recomputed k* values with the validated historical
# result. This never enters the optimization and never overwrites selected_k.
if (isTRUE(ENFORCE_VALIDATED_K_REGRESSION)) {
  observed_k <- comparison_specific_k
  missing_cmp <- setdiff(names(VALIDATED_REFERENCE_K), names(observed_k))
  if (length(missing_cmp)) {
    stop(
      "Regression check could not find comparisons: ",
      paste(missing_cmp, collapse = ", ")
    )
  }

  observed_k <- observed_k[names(VALIDATED_REFERENCE_K)]
  if (!identical(unname(observed_k), unname(VALIDATED_REFERENCE_K))) {
    stop(
      "EMPIRICAL CUTOFF REGRESSION DETECTED. Recomputed k*: ",
      paste(names(observed_k), observed_k, sep = "=", collapse = ", "),
      ". Validated k*: ",
      paste(names(VALIDATED_REFERENCE_K), VALIDATED_REFERENCE_K, sep = "=", collapse = ", "),
      ". Do not propagate these outputs downstream until the implementation change is explained."
    )
  }

  message(
    "Validated cutoff regression check passed: ",
    paste(names(observed_k), observed_k, sep = "=", collapse = ", ")
  )
}

likelihood_all <- bind_rows(corroboration_rows)
leading_all <- bind_rows(selected_site_rows)
remainder_all <- bind_rows(remainder_site_rows)
pareto_cost_all <- bind_rows(pareto_cost_rows)

summary_cutoff_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_Empirical_Cutoffs.csv"
)

summary_likelihood_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_NB1_NB2_Likelihood_All.csv"
)

summary_lead_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_Leading_Edge_Sites_All.csv"
)

summary_rem_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_Remainder_Sites_All.csv"
)

summary_cost_path <- file.path(
  SUMMARY_TAB_DIR,
  "Table_Pareto_Remainder_Crossing_Cost_Sites_All.csv"
)

write_csv(cutoff_summary, summary_cutoff_path)
write_csv(likelihood_all, summary_likelihood_path)
write_csv(leading_all, summary_lead_path)
write_csv(remainder_all, summary_rem_path)
write_csv(pareto_cost_all, summary_cost_path)

summary_fig_path <- file.path(
  SUMMARY_FIG_DIR,
  "Figure_Empirical_Cutoff_Summary.png"
)

make_summary_figure(
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

expected_figure_paths <- c(
  expected_figure_paths,
  expected_summary_figs
)

expected_table_paths <- c(
  expected_table_paths,
  summary_cutoff_path,
  summary_likelihood_path,
  summary_lead_path,
  summary_rem_path,
  summary_cost_path
)

write_methods_manuscript(cutoff_summary)
write_figure_legends()

# Manifest before ZIP creation.
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

fig_zip <- file.path(OUT_ROOT, "Figures_All.zip")
tab_zip <- file.path(OUT_ROOT, "Tables_All.zip")
all_zip <- file.path(OUT_ROOT, "Empirical_Cutoff_All_Outputs.zip")

create_zip_from_files(
  fig_zip,
  expected_figure_paths
)

create_zip_from_files(
  tab_zip,
  expected_table_paths
)

all_files_for_zip <- c(
  expected_figure_paths,
  expected_table_paths,
  file.path(OUT_ROOT, "Methods_Manuscript.md"),
  file.path(OUT_ROOT, "Figure_Legends.md")
)

create_zip_from_files(
  all_zip,
  all_files_for_zip
)

verify_outputs(c(fig_zip, tab_zip, all_zip))

# =============================================================================
# FINAL CONSOLE SUMMARY
# =============================================================================

message("============================================================")
message("EMPIRICAL EVS CUTOFF ANALYSIS COMPLETE")
message("PCA input: DESeq2 median-of-ratios normalized counts entered directly.")
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
message("Figures per comparison: 2 manuscript figures (+ PDF copies when enabled).")
message("All per-comparison figures/tables remain in their comparison folders.")
message("Figures ZIP: ", fig_zip)
message("Tables ZIP: ", tab_zip)
message("Complete ZIP: ", all_zip)
message("Build: ", SCRIPT_BUILD)
message("Methods: ", file.path(OUT_ROOT, "Methods_Manuscript.md"))
message("Figure legends: ", file.path(OUT_ROOT, "Figure_Legends.md"))
message("============================================================")
