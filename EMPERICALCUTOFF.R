#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# BUILD: 2026-08-24_PAIRWISE_ONLY_CLEAN_KEY
# =============================================================================
# PC1-NB REGIME GEOMETRY + COMPARISON-SPECIFIC WEIGHTED PARETO EVS
# =============================================================================
#
# PURPOSE
# -------
# Identify variance regimes along independently ranked PC1 axes and determine
# a separate comparison-specific top-k eigenvector-splitting cutoff for each
# RT/ZT pair. Each k* maximizes retention of Joint and permissible Disjoint sites
# while minimizing recruitment of sites whose opposite-arm rank falls into the
# Remainder.
#
#
# 1. ARM-SPECIFIC PC1 RANKING
# ---------------------------
# For each experimental arm g, counts are converted to CPM and transformed:
#
#     X_ig = log(1 + CPM_ig)
#
# PCA/SVD is performed independently within each arm.  Let v_i1 denote the
# loading of feature i on PC1. Features are ranked by ASCENDING absolute PC1
# loading:
#
#     rank_order = order(|v_i1|, decreasing = FALSE)
#
# Thus:
#
#     low rank  = small |PC1 loading|
#     high rank = large |PC1 loading|
#
# Feature-level contribution to PC1 variance is:
#
#     P_i = lambda_1 * v_i1^2
#
# where lambda_1 is the PC1 eigenvalue.
#
#
# 2. RAW-COUNT VARIANCE GEOMETRY
# ------------------------------
# Within each arm, the empirical sample variance of the raw counts is:
#
#     s_raw,g^2(i) = Var_j(Y_ij | g)
#
# After ordering features by PC1 rank:
#
#     y_g(r) = log[1 + s_raw,g^2(r)]
#
# A smoothing spline is used to display the empirical variance geometry along
# the ranked axis.
#
#
# 3. POOLED WITHIN-GROUP NORMALIZED VARIANCE
# ------------------------------------------
# Counts are normalized globally with DESeq2 size factors.  The pooled
# within-group empirical variance for feature i is:
#
#                    sum_g sum_{j in g}(y_ij - ybar_ig)^2
#     V_pool,i =     -------------------------------------
#                              sum_g(n_g - 1)
#
# where y_ij is the DESeq2-normalized count.
#
# For each arm g:
#
#     mu_ig = mean_j(y_ij | g)
#
# and the excess-over-Poisson variance is:
#
#     E_ig = max(V_pool,i - mu_ig, 0)
#
#
# 4. CUMULATIVE PC1-NB VARIANCE-MASS DIVERGENCE
# ----------------------------------------------
# PC1 variance contribution and NB excess variance are converted to rank-wise
# probability masses:
#
#     p_g(r) = P_g(r) / sum_r P_g(r)
#
#     q_g(r) = E_g(r) / sum_r E_g(r)
#
# Their cumulative distributions are:
#
#     F_P,g(r) = sum_{j <= r} p_g(j)
#
#     F_E,g(r) = sum_{j <= r} q_g(j)
#
# and cumulative divergence is:
#
#     D_g(r) = F_E,g(r) - F_P,g(r)
#
# A shared two-knot continuous linear spline is fit jointly to D_g(r) for all
# eight arms:
#
#     D_g(x) =
#       beta_0g
#       + beta_1g*x
#       + gamma_1g*(x-c1)_+
#       + gamma_2g*(x-c2)_+
#
# with:
#
#     x = (rank - 1)/(N - 1)
#
# and shared c1 and c2 across arms.
#
# Rank regimes are defined as:
#
#     rank < c1          = Remainder
#
#     c1 <= rank <= c2   = Divergence interval
#
#     rank > c2          = Leading-edge regime
#
#
# 5. TOP-k EIGENVECTOR SPLITTING
# ------------------------------
# For each control/treatment pair and candidate depth k:
#
#     S_C(k) = top-k features in the control PC1 ranking
#
#     S_T(k) = top-k features in the treatment PC1 ranking
#
# Joint and Disjoint classes are:
#
#     Joint      = S_C(k) intersection S_T(k)
#
#     Disjoint C = S_C(k) \ S_T(k)
#
#     Disjoint T = S_T(k) \ S_C(k)
#
# Candidate k is constrained by the selecting arm's Leading-edge regime:
#
#     1 <= k <= N - c2
#
# Therefore every site entering S_C(k) or S_T(k) originates from rank > c2
# in the arm that selects it.
#
#
# 6. OPPOSITE-ARM CLASSIFICATION OF DISJOINT SITES
# ------------------------------------------------
# For a Disjoint site, its rank in the opposite arm determines whether it is
# retained:
#
#     r_opposite > c2
#         = opposite-arm Leading edge
#         = retained
#
#     c1 <= r_opposite <= c2
#         = opposite-arm Divergence interval
#         = retained
#
#     r_opposite < c1
#         = opposite-arm Remainder
#         = contamination
#
# Thus the Divergence interval is permissible.  The Remainder boundary c1 is
# the contamination boundary.
#
# For candidate k:
#
#     G(k) =
#       N_Joint(k)
#       + N_Disjoint,opposite-LE(k)
#       + N_Disjoint,opposite-Divergence(k)
#
#     R(k) =
#       N_Disjoint,opposite-Remainder(k)
#
#
# 7. WEIGHTED PARETO OPTIMIZATION
# -------------------------------
# Each candidate k is represented by the pair:
#
#     [R(k), G(k)]
#
# A candidate is Pareto-optimal if no other candidate has both:
#
#     G(k') >= G(k)
#
# and
#
#     R(k') <= R(k)
#
# with at least one strict inequality.
#
# On the Pareto frontier, benefit and contamination are normalized:
#
#     G_norm(k) =
#       [G(k) - G_min] / [G_max - G_min]
#
#     R_norm(k) =
#       [R(k) - R_min] / [R_max - R_min]
#
# Weighted utility is:
#
#     U(k) =
#       w_G * G_norm(k)
#       - w_R * R_norm(k)
#
# with:
#
#     w_G = 1
#     w_R = 1
#
# The selected cutoff is:
#
#     k* = argmax_k U(k)
#
# over Pareto-optimal candidates.
#
# Ties are resolved by:
#
#     1. greater G(k)
#     2. lower R(k)
#     3. larger k
#
#
# 8. FINAL COMPARISON-SPECIFIC CUTOFFS
# --------------------------------------
# The four RT/ZT comparisons retain their own independently selected weighted-
# Pareto k*. Each comparison-specific optimum is carried forward to its own
# site classification, tables, and figure.
#
#
# FIGURES
# -------
#   Figure_Overall.png
#   Figure_RT0_ZT6.png
#   Figure_RT2_ZT8.png
#   Figure_RT4_ZT10.png
#   Figure_RT8_ZT14.png
#
# Each comparison-specific figure contains:
#
#   A. Raw-count variance geometry
#   B. Cumulative PC1-NB variance-mass divergence
#   C. Comparison-specific weighted Pareto top-k optimization
#
# Figure_Overall contains only the shared experiment-wide variance/divergence
# geometry (Panels A-B).
#
#
# OUTPUTS
# -------
#   Table_Key_Results.csv
#   Table_Timepoints.csv
#   Table_Cutoff_Optimization.csv
#   Selected_LeadingEdge_Sites.csv
#   Excluded_Remainder_Crossing_Sites.csv
#   Figures_All.zip
#
# =============================================================================
# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <- "/root/REAPER98632/exports/pc1_nb_manuscript_pairwise_only_clean_key"
FIG_DIR  <- file.path(OUT_ROOT, "Figures")

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

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


# Display-only smoothing. These values never determine boundaries/cutoffs.
DISPLAY_VAR_SPAR <- 0.72
DISPLAY_D_SPAR   <- 0.72

PNG_WIDTH_IN <- 15
PNG_DPI      <- 360

# Equal normalized weights for retained-site benefit and remainder
# contamination.
BENEFIT_WEIGHT <- 1.0
CONTAMINATION_WEIGHT <- 1.0


# =============================================================================
# COLORS
# =============================================================================

COL <- list(
  control = "#0072B2",
  treatment = "#009E73",
  raw = "#00A6A6",
  divergence = "#6A3D9A",
  fit = "#252525",

  remainder = "#DCE6F2",
  interval = "#FFF0B3",
  leading = "#D8F3E7",

  c1 = "#D73027",
  c2 = "#1A9850",
  selected = "#B5179E",
  candidate_line = "#9AA0A6",
  pareto = "#5E3C99"
)


# =============================================================================
# BASIC HELPERS
# =============================================================================

theme_manuscript <- function(base_size = 11.5) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.1),
      plot.subtitle = element_text(size = base_size - 0.3, margin = margin(b = 4)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#333333"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1.2),
      legend.box = "horizontal",
      panel.border = element_rect(color = "#B7B7B7", fill = NA, linewidth = 0.45),
      panel.grid = element_blank(),
      plot.margin = margin(6, 8, 6, 8)
    )
}

save_panels <- function(
    plots,
    path,
    height_in,
    key_grob = NULL,
    key_height_in = 0.48) {

  grDevices::png(
    filename = path,
    width = PNG_WIDTH_IN,
    height = height_in,
    units = "in",
    res = PNG_DPI,
    bg = "white"
  )

  grid.newpage()

  if (is.null(key_grob)) {
    layout_heights <- unit(rep(1, length(plots)), "null")
    n_rows <- length(plots)
  } else {
    layout_heights <- unit.c(
      unit(rep(1, length(plots)), "null"),
      unit(key_height_in, "in")
    )
    n_rows <- length(plots) + 1L
  }

  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = n_rows,
        ncol = 1L,
        heights = layout_heights
      )
    )
  )

  for (i in seq_along(plots)) {
    print(
      plots[[i]],
      vp = viewport(
        layout.pos.row = i,
        layout.pos.col = 1L
      )
    )
  }

  if (!is.null(key_grob)) {
    pushViewport(
      viewport(
        layout.pos.row = n_rows,
        layout.pos.col = 1L
      )
    )
    grid.draw(key_grob)
    popViewport()
  }

  dev.off()
}

make_figure_key_grob <- function(items, fontsize = 9.4) {
  n <- length(items)
  grobs <- list()

  for (i in seq_along(items)) {
    item <- items[[i]]
    center <- (i - 0.5) / n
    x0 <- center - 0.070
    x1 <- center - 0.020
    xt <- center - 0.005

    grobs[[length(grobs) + 1L]] <- segmentsGrob(
      x0 = unit(x0, "npc"),
      x1 = unit(x1, "npc"),
      y0 = unit(0.5, "npc"),
      y1 = unit(0.5, "npc"),
      gp = gpar(
        col = item$color,
        lty = item$lty,
        lwd = item$lwd %||% 2.0
      )
    )

    if (!is.null(item$pch)) {
      grobs[[length(grobs) + 1L]] <- pointsGrob(
        x = unit((x0 + x1) / 2, "npc"),
        y = unit(0.5, "npc"),
        pch = item$pch,
        size = unit(2.5, "mm"),
        gp = gpar(
          col = item$color,
          fill = item$fill %||% item$color
        )
      )
    }

    grobs[[length(grobs) + 1L]] <- textGrob(
      label = item$label,
      x = unit(xt, "npc"),
      y = unit(0.5, "npc"),
      just = c("left", "center"),
      gp = gpar(
        col = "#222222",
        fontsize = fontsize
      )
    )
  }

  do.call(grobTree, grobs)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

make_comparison_key <- function(selected_k) {
  make_figure_key_grob(
    list(
      list(label = "Control", color = COL$control, lty = "solid"),
      list(label = "Treatment", color = COL$treatment, lty = "solid"),
      list(label = "c1", color = COL$c1, lty = "dashed"),
      list(label = "c2", color = COL$c2, lty = "longdash"),
      list(
        label = paste0("k*=", selected_k),
        color = COL$selected,
        lty = "dotdash",
        pch = 23L,
        fill = COL$selected
      ),
      list(label = "Pareto frontier", color = COL$pareto, lty = "solid")
    )
  )
}

make_overall_key <- function() {
  make_figure_key_grob(
    list(
      list(label = "Raw variance", color = COL$raw, lty = "solid"),
      list(label = "D(r)", color = COL$divergence, lty = "solid"),
      list(label = "Shared fit", color = COL$fit, lty = "solid"),
      list(label = "c1", color = COL$c1, lty = "dashed"),
      list(label = "c2", color = COL$c2, lty = "longdash")
    )
  )
}

rank_cutoff_from_k <- function(N, k) {
  if (!is.finite(k) || k < 1L || k > N) return(NA_integer_)
  as.integer(N - k + 1L)
}


# =============================================================================
# DATA INPUT
# =============================================================================

read_count_matrix <- function(path, group_patterns) {
  if (!file.exists(path)) {
    stop("Count file does not exist: ", path)
  }

  raw_df <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(raw_df) < 1L || ncol(raw_df) < 2L) {
    stop("Count file is empty or malformed.")
  }

  sample_idx <- sort(
    unique(
      unlist(
        lapply(
          group_patterns,
          function(pattern) {
            grep(pattern, colnames(raw_df))
          }
        )
      )
    )
  )

  if (length(sample_idx) == 0L) {
    stop("No sample columns matched GROUP_PATTERNS.")
  }

  if (1L %in% sample_idx) {
    stop("Column 1 matched a sample pattern; column 1 must contain feature IDs.")
  }

  count_df <- raw_df[, sample_idx, drop = FALSE]

  count_mat <- do.call(
    cbind,
    lapply(
      count_df,
      function(x) {
        suppressWarnings(
          as.numeric(
            trimws(
              as.character(x)
            )
          )
        )
      }
    )
  )

  colnames(count_mat) <- colnames(count_df)
  storage.mode(count_mat) <- "numeric"

  n_bad <- sum(!is.finite(count_mat))

  if (n_bad > 0L) {
    message("Replacing ", n_bad, " non-finite count entries with 0.")
    count_mat[!is.finite(count_mat)] <- 0
  }

  count_mat <- pmax(count_mat, 0)

  feature_ids <- trimws(as.character(raw_df[[1L]]))

  blank <- is.na(feature_ids) | feature_ids == ""

  if (any(blank)) {
    feature_ids[blank] <- paste0("__feature_row_", which(blank))
  }

  feature_ids <- make.unique(feature_ids, sep = "__dup_")
  rownames(count_mat) <- feature_ids

  keep <- rowSums(count_mat) > 0
  count_mat <- count_mat[keep, , drop = FALSE]

  if (nrow(count_mat) < 2L) {
    stop("Fewer than two nonzero features remain after filtering.")
  }

  count_mat
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
    stop(
      "Unassigned samples: ",
      paste(sample_names[is.na(assigned)], collapse = ", ")
    )
  }

  factor(
    assigned,
    levels = names(group_patterns)
  )
}


# =============================================================================
# NORMALIZATION AND VARIANCE
# =============================================================================

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

normalize_deseq2_global <- function(count_mat, group_labels) {
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


# =============================================================================
# PC1 RANKING
# =============================================================================

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


# =============================================================================
# RAW-COUNT VARIANCE GEOMETRY
# =============================================================================

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


# =============================================================================
# GROUP-SPECIFIC PC1-NB DIVERGENCE
# =============================================================================

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

  rank_matrix <- normalize_cpm_log1p(raw_counts_arm)

  pc1 <- compute_pc1_rank(rank_matrix)
  rank_order <- pc1$rank_order

  geometry <- compute_raw_variance_geometry(
    raw_counts_arm = raw_counts_arm,
    rank_order = rank_order
  )

  mu_norm <- rowMeans(
    normalized_counts_arm,
    na.rm = TRUE
  )

  mu_norm[!is.finite(mu_norm)] <- 0
  mu_norm <- pmax(mu_norm, 0)

  P_ranked <- pc1$pc1_variance_contribution[rank_order]
  mu_ranked <- mu_norm[rank_order]
  V_pool_ranked <- pooled_variance[rank_order]

  E_ranked <- pmax(
    V_pool_ranked - mu_ranked,
    0
  )

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
    data = df
  )
}


# =============================================================================
# SHARED TWO-KNOT DIVERGENCE MODEL
# =============================================================================

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

  # 5,000 here is ONLY the knot-search grid size; it is not an EVS cutoff.
  opt_n <- min(5000L, N)

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


# =============================================================================
# TOP-k OPTIMIZATION AND JOINT / DISJOINT CLASSIFICATION
# =============================================================================

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

  retained_for_analysis <- (
    joint |
    disjoint_opposite_leading_edge |
    disjoint_opposite_divergence
  )

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
    "Excluded_",
    control_group,
    "_OppositeRemainder"
  )

  analysis_class[
    treatment_only &
    rC < c1
  ] <- paste0(
    "Excluded_",
    treatment_group,
    "_OppositeRemainder"
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

  if (
    any(
      out$retained_for_analysis &
      out$cross_into_remainder
    )
  ) {
    stop(
      "A retained site was also classified as a remainder crossing in ",
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
# RESULTS TABLES
# =============================================================================

build_timepoint_table <- function(
    comparisons,
    pair_optima,
    selected_summaries,
    c1,
    c2,
    N) {

  rows <- list()

  for (comparison_name in names(comparisons)) {
    mapping <- comparisons[[comparison_name]]
    pair_opt <- pair_optima[[comparison_name]]
    selected_summary <- selected_summaries[[comparison_name]]

    rows[[length(rows) + 1L]] <- data.frame(
      comparison = comparison_name,
      control_group = unname(mapping[["control"]]),
      treatment_group = unname(mapping[["treatment"]]),
      c1 = c1,
      c2 = c2,
      max_candidate_k = N - c2,
      selected_k = pair_opt$selected_k,
      selected_cutoff_rank = rank_cutoff_from_k(N, pair_opt$selected_k),
      weighted_utility = pair_opt$selected_utility,
      good_norm = pair_opt$selected_good_norm,
      remainder_norm = pair_opt$selected_remainder_norm,
      max_zero_remainder_crossing_k = pair_opt$max_zero_remainder_crossing_k,
      benefit_weight = pair_opt$benefit_weight,
      contamination_weight = pair_opt$contamination_weight,
      selected_joint_n = selected_summary$joint_n,
      selected_disjoint_opposite_le_n = selected_summary$disjoint_opposite_le_n,
      selected_disjoint_opposite_divergence_n = selected_summary$disjoint_opposite_divergence_n,
      selected_permissible_disjoint_n = selected_summary$permissible_disjoint_n,
      selected_good_n = selected_summary$good_n,
      selected_remainder_cross_n = selected_summary$remainder_cross_n,
      selected_union_n = selected_summary$union_n,
      selected_retained_fraction = selected_summary$retained_fraction,
      selected_remainder_cross_fraction = selected_summary$remainder_cross_fraction,
      stringsAsFactors = FALSE
    )
  }

  bind_rows(rows)
}

build_key_table <- function(
    timepoint_table,
    N,
    c1,
    c2,
    shared_sse) {

  k_map <- stats::setNames(
    timepoint_table$selected_k,
    timepoint_table$comparison
  )

  data.frame(
    n_features = N,
    shared_c1 = c1,
    shared_c2 = c2,
    remainder_size = c1 - 1L,
    divergence_interval_size = c2 - c1 + 1L,
    leading_edge_size = N - c2,
    max_candidate_k = N - c2,
    benefit_weight = BENEFIT_WEIGHT,
    contamination_weight = CONTAMINATION_WEIGHT,
    RT0_ZT6_k = unname(k_map[["RT0_ZT6"]]),
    RT2_ZT8_k = unname(k_map[["RT2_ZT8"]]),
    RT4_ZT10_k = unname(k_map[["RT4_ZT10"]]),
    RT8_ZT14_k = unname(k_map[["RT8_ZT14"]]),
    shared_model_SSE = shared_sse,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# FIGURE DATA
# =============================================================================

rankwise_median <- function(
    group_results,
    column) {

  mat <- do.call(
    cbind,
    lapply(
      group_results,
      function(z) {
        z$data[[column]]
      }
    )
  )

  apply(
    mat,
    1L,
    median,
    na.rm = TRUE
  )
}


build_overall_figure_data <- function(
    group_results,
    knot_fit) {

  N <- nrow(
    group_results[[1L]]$data
  )

  data.frame(
    rank = seq_len(N),

    raw_variance = rankwise_median(
      group_results,
      "display_log1p_raw_empirical_variance"
    ),

    D = rankwise_median(
      group_results,
      "display_D"
    ),

    D_fit = apply(
      knot_fit$fitted,
      1L,
      median,
      na.rm = TRUE
    ),

    stringsAsFactors = FALSE
  )
}


add_rank_regions <- function(
    p,
    c1,
    c2,
    N) {

  p +
    annotate(
      "rect",
      xmin = 1,
      xmax = c1,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$remainder,
      alpha = 0.50
    ) +
    annotate(
      "rect",
      xmin = c1,
      xmax = c2,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$interval,
      alpha = 0.46
    ) +
    annotate(
      "rect",
      xmin = c2,
      xmax = N,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$leading,
      alpha = 0.52
    ) +
    annotate(
      "text",
      x = (1 + c1) / 2,
      y = Inf,
      label = "Remainder",
      vjust = 1.25,
      fontface = "bold",
      size = 3.0
    ) +
    annotate(
      "text",
      x = (c1 + c2) / 2,
      y = Inf,
      label = "Divergence",
      vjust = 1.25,
      fontface = "bold",
      size = 3.0
    ) +
    annotate(
      "text",
      x = (c2 + N) / 2,
      y = Inf,
      label = "Leading edge",
      vjust = 1.25,
      fontface = "bold",
      size = 3.0
    )
}


make_cutoff_panel <- function(
    scan_df,
    selected_k,
    title_text) {

  selected <- scan_df %>%
    filter(k == selected_k) %>%
    slice(1L)

  if (nrow(selected) != 1L) {
    stop("Selected k not present in cutoff scan.")
  }

  frontier <- scan_df %>%
    filter(is_pareto) %>%
    arrange(remainder_cross_n, good_n, k) %>%
    distinct(remainder_cross_n, good_n, .keep_all = TRUE)

  ggplot() +
    geom_path(
      data = scan_df,
      aes(
        x = remainder_cross_n,
        y = good_n,
        group = 1
      ),
      color = COL$candidate_line,
      linewidth = 0.45,
      alpha = 0.40
    ) +
    geom_point(
      data = scan_df,
      aes(
        x = remainder_cross_n,
        y = good_n,
        fill = k
      ),
      shape = 21,
      color = "white",
      stroke = 0.15,
      size = 1.55,
      alpha = 0.88
    ) +
    geom_path(
      data = frontier,
      aes(
        x = remainder_cross_n,
        y = good_n,
        group = 1
      ),
      color = COL$pareto,
      linewidth = 1.30
    ) +
    geom_point(
      data = selected,
      aes(
        x = remainder_cross_n,
        y = good_n
      ),
      shape = 23,
      color = COL$selected,
      fill = COL$selected,
      size = 4.7,
      stroke = 1.15
    ) +
    annotate(
      "text",
      x = selected$remainder_cross_n,
      y = selected$good_n,
      label = paste0("k*=", selected_k),
      vjust = 1.85,
      size = 3.0,
      fontface = "bold",
      color = COL$selected
    ) +
    scale_fill_gradientn(
      colors = c(
        "#2C7BB6",
        "#00A6CA",
        "#00CCBC",
        "#90EB9D",
        "#F9D057",
        "#F29E2E",
        "#D7191C"
      ),
      guide = "none"
    ) +
    labs(
      title = title_text,
      subtitle = "Weighted utility: U = Gnorm - Rnorm (1:1)",
      x = "Opposite-arm remainder crossings (rank < c1)",
      y = "Joint + permissible disjoint sites"
    ) +
    theme_manuscript() +
    theme(legend.position = "none")
}


# =============================================================================
# OVERALL MAIN FIGURE
# =============================================================================

make_overall_figure <- function(
    overall_df,
    c1,
    c2,
    out_file) {

  N <- nrow(overall_df)

  boundary_lines <- data.frame(
    rank = c(c1, c2),
    key = c("c1", "c2"),
    stringsAsFactors = FALSE
  )

  boundary_colors <- c(
    "c1" = COL$c1,
    "c2" = COL$c2
  )

  boundary_types <- c(
    "c1" = "dashed",
    "c2" = "longdash"
  )

  pA <- add_rank_regions(ggplot(), c1, c2, N) +
    geom_vline(
      data = boundary_lines,
      aes(xintercept = rank, color = key, linetype = key),
      linewidth = 0.85
    ) +
    geom_line(
      data = overall_df,
      aes(x = rank, y = raw_variance),
      color = COL$raw,
      linewidth = 1.20
    ) +
    scale_color_manual(values = boundary_colors) +
    scale_linetype_manual(values = boundary_types) +
    labs(
      title = "A. Raw-count variance geometry",
      subtitle = "Empirical variance along the PC1-ranked axis",
      x = "PC1 rank: low |loading| -> high |loading|",
      y = "Smoothed log(1 + raw-count variance)"
    ) +
    theme_manuscript() +
    theme(legend.position = "none")

  pB <- add_rank_regions(ggplot(), c1, c2, N) +
    geom_vline(
      data = boundary_lines,
      aes(xintercept = rank, color = key, linetype = key),
      linewidth = 0.85
    ) +
    geom_hline(
      yintercept = 0,
      color = "#777777",
      linetype = "dotted",
      linewidth = 0.35
    ) +
    geom_line(
      data = overall_df,
      aes(x = rank, y = D),
      color = COL$divergence,
      linewidth = 0.95,
      alpha = 0.72
    ) +
    geom_line(
      data = overall_df,
      aes(x = rank, y = D_fit),
      color = COL$fit,
      linewidth = 1.30
    ) +
    scale_color_manual(values = boundary_colors) +
    scale_linetype_manual(values = boundary_types) +
    labs(
      title = "B. Cumulative PC1-NB variance-mass divergence",
      subtitle = "Shared c1/c2 fitted across all 8 arms",
      x = "PC1 rank",
      y = "D(r) = F_E(r) - F_P(r)"
    ) +
    theme_manuscript() +
    theme(legend.position = "none")

  save_panels(
    list(pA, pB),
    out_file,
    height_in = 8.3,
    key_grob = make_overall_key()
  )
}


# =============================================================================
# TIME-POINT MAIN FIGURES
# =============================================================================

make_timepoint_figure <- function(
    comparison_name,
    mapping,
    group_results,
    pair_scan,
    pair_optimum,
    c1,
    c2,
    out_file) {

  control_group <- unname(mapping[["control"]])
  treatment_group <- unname(mapping[["treatment"]])
  control <- group_results[[control_group]]$data
  treatment <- group_results[[treatment_group]]$data
  N <- nrow(control)

  selected_k <- pair_optimum$selected_k
  selected_rank <- rank_cutoff_from_k(N, selected_k)
  selected_key <- paste0("k*=", selected_k)

  boundary_lines <- data.frame(
    rank = c(c1, c2, selected_rank),
    key = c("c1", "c2", selected_key),
    stringsAsFactors = FALSE
  )

  boundary_colors <- c(
    "c1" = COL$c1,
    "c2" = COL$c2,
    stats::setNames(COL$selected, selected_key)
  )

  boundary_types <- c(
    "c1" = "dashed",
    "c2" = "longdash",
    stats::setNames("dotdash", selected_key)
  )

  pA <- add_rank_regions(ggplot(), c1, c2, N) +
    geom_vline(
      data = boundary_lines,
      aes(xintercept = rank, color = key, linetype = key),
      linewidth = 0.84
    ) +
    geom_line(
      data = control,
      aes(x = rank, y = display_log1p_raw_empirical_variance),
      color = COL$control,
      linewidth = 1.10
    ) +
    geom_line(
      data = treatment,
      aes(x = rank, y = display_log1p_raw_empirical_variance),
      color = COL$treatment,
      linewidth = 1.10
    ) +
    scale_color_manual(values = boundary_colors) +
    scale_linetype_manual(values = boundary_types) +
    labs(
      title = paste0("A. ", comparison_name, " raw-count variance geometry"),
      subtitle = "Empirical variance along independently ranked PC1 axes",
      x = "PC1 rank: low |loading| -> high |loading|",
      y = "Smoothed log(1 + raw-count variance)"
    ) +
    theme_manuscript() +
    theme(legend.position = "none")

  pB <- add_rank_regions(ggplot(), c1, c2, N) +
    geom_vline(
      data = boundary_lines,
      aes(xintercept = rank, color = key, linetype = key),
      linewidth = 0.84
    ) +
    geom_hline(
      yintercept = 0,
      color = "#777777",
      linetype = "dotted",
      linewidth = 0.35
    ) +
    geom_line(
      data = control,
      aes(x = rank, y = display_D),
      color = COL$control,
      linewidth = 1.08
    ) +
    geom_line(
      data = treatment,
      aes(x = rank, y = display_D),
      color = COL$treatment,
      linewidth = 1.08
    ) +
    scale_color_manual(values = boundary_colors) +
    scale_linetype_manual(values = boundary_types) +
    labs(
      title = paste0("B. ", comparison_name, " cumulative PC1-NB divergence"),
      subtitle = "Shared c1/c2",
      x = "PC1 rank",
      y = "D(r) = F_E(r) - F_P(r)"
    ) +
    theme_manuscript() +
    theme(legend.position = "none")

  pC <- make_cutoff_panel(
    scan_df = pair_scan,
    selected_k = selected_k,
    title_text = paste0("C. ", comparison_name, " weighted Pareto optimization")
  )

  save_panels(
    list(pA, pB, pC),
    out_file,
    height_in = 11.8,
    key_grob = make_comparison_key(selected_k)
  )
}


# =============================================================================
# RUN ANALYSIS
# =============================================================================

count_mat <- read_count_matrix(
  COUNT_FILE,
  GROUP_PATTERNS
)

group_labels <- assign_groups(
  colnames(count_mat),
  GROUP_PATTERNS
)

message(
  "Count matrix: ",
  nrow(count_mat),
  " features x ",
  ncol(count_mat),
  " samples"
)

message(
  "Groups: ",
  paste(
    levels(group_labels),
    collapse = ", "
  )
)

N <- nrow(
  count_mat
)


# -------------------------------------------------------------------------
# Global DESeq2 normalization and pooled within-group empirical variance.
# -------------------------------------------------------------------------

deseq <- normalize_deseq2_global(
  count_mat = count_mat,
  group_labels = group_labels
)

normalized_counts <- deseq$normalized_counts

pooled <- compute_pooled_within_group_variance(
  normalized_counts =
    normalized_counts,
  group_labels =
    group_labels
)

message(
  "Pooled within-group residual df: ",
  pooled$residual_df
)

# -------------------------------------------------------------------------
# Arm-specific PC1 ranking, raw variance geometry, and cumulative divergence.
# -------------------------------------------------------------------------

group_results <- vector(
  "list",
  length(
    levels(
      group_labels
    )
  )
)

names(
  group_results
) <- levels(
  group_labels
)

for (g in levels(
  group_labels
)) {
  idx <- which(
    group_labels ==
    g
  )

  if (
    length(idx) < 2L
  ) {
    stop(
      "Not enough samples in group ",
      g
    )
  }

  message(
    "Analyzing ",
    g,
    "..."
  )

  group_results[[g]] <-
    compute_group_analysis(
      group_name = g,

      raw_counts_arm =
        count_mat[
          ,
          idx,
          drop = FALSE
        ],

      normalized_counts_arm =
        normalized_counts[
          ,
          idx,
          drop = FALSE
        ],

      pooled_variance =
        pooled$variance
    )
}

# -------------------------------------------------------------------------
# Shared c1/c2 cumulative-divergence model.
# -------------------------------------------------------------------------

message(
  "Fitting shared two-knot cumulative-divergence model..."
)

knot_fit <- fit_shared_knots(
  group_results
)

C1 <- knot_fit$c1
C2 <- knot_fit$c2

MAX_CANDIDATE_K <- as.integer(
  N - C2
)

message(
  "Shared c1 = ",
  C1
)

message(
  "Shared c2 = ",
  C2
)

message(
  "Own-arm leading-edge candidate limit: k <= ",
  MAX_CANDIDATE_K
)


# -------------------------------------------------------------------------
# Pairwise top-k scans.
#
# IMPORTANT:
# - selection originates only from r > c2;
# - opposite-arm divergence c1..c2 is permissible;
# - only opposite-arm remainder r < c1 is penalized.
# -------------------------------------------------------------------------

message(
  "Scanning candidate k values using c1 as the opposite-arm contamination boundary..."
)

pair_scans <- list()
pair_optima <- list()

for (comparison_name in names(
  COMPARISONS
)) {
  mapping <- COMPARISONS[[
    comparison_name
  ]]

  control_group <- unname(
    mapping[["control"]]
  )

  treatment_group <- unname(
    mapping[["treatment"]]
  )

  raw_scan <- scan_pair_cutoffs(
    control_df =
      group_results[[
        control_group
      ]]$data,

    treatment_df =
      group_results[[
        treatment_group
      ]]$data,

    c1 = C1,
    c2 = C2,

    comparison_name =
      comparison_name,

    control_group =
      control_group,

    treatment_group =
      treatment_group
  )

  opt <- select_weighted_pareto_optimum(
    raw_scan,
    good_col = "good_n",
    cost_col = "remainder_cross_n"
  )

  pair_scans[[
    comparison_name
  ]] <- opt$scan

  pair_optima[[
    comparison_name
  ]] <- opt

  selected_row <- opt$scan %>%
    filter(
      k ==
      opt$selected_k
    ) %>%
    slice(1L)

  message(
    comparison_name,
    ": weighted k*=",
    opt$selected_k,
    "; retained=",
    selected_row$good_n[1L],
    "; divergence-disjoint=",
    selected_row$disjoint_opposite_divergence_n[1L],
    "; remainder crossings=",
    selected_row$remainder_cross_n[1L],
    "; utility=",
    signif(opt$selected_utility, 5)
  )
}

# -------------------------------------------------------------------------
# Classify and export sites at each comparison-specific weighted-Pareto k*.
# Each comparison uses its own weighted-Pareto optimum downstream.
# -------------------------------------------------------------------------

selected_rows <- list()
excluded_rows <- list()

selected_summaries <- list()

for (comparison_name in names(COMPARISONS)) {
  mapping <- COMPARISONS[[comparison_name]]
  control_group <- unname(mapping[["control"]])
  treatment_group <- unname(mapping[["treatment"]])
  selected_k <- pair_optima[[comparison_name]]$selected_k

  selected_all <- classify_pair_at_k(
    control_df = group_results[[control_group]]$data,
    treatment_df = group_results[[treatment_group]]$data,
    k = selected_k,
    c1 = C1,
    c2 = C2,
    comparison_name = comparison_name,
    control_group = control_group,
    treatment_group = treatment_group
  )

  selected_summaries[[comparison_name]] <- summarize_classification(selected_all)

  selected_rows[[comparison_name]] <- selected_all %>%
    filter(retained_for_analysis)

  excluded_rows[[comparison_name]] <- selected_all %>%
    filter(cross_into_remainder)
}

selected_sites <- bind_rows(
  selected_rows
)

excluded_sites <- bind_rows(
  excluded_rows
)

write.csv(
  selected_sites,
  file.path(
    OUT_ROOT,
    "Selected_LeadingEdge_Sites.csv"
  ),
  row.names = FALSE
)

write.csv(
  excluded_sites,
  file.path(
    OUT_ROOT,
    "Excluded_Remainder_Crossing_Sites.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Tables.
# -------------------------------------------------------------------------

timepoint_table <- build_timepoint_table(
  comparisons = COMPARISONS,
  pair_optima = pair_optima,
  selected_summaries = selected_summaries,
  c1 = C1,
  c2 = C2,
  N = N
)

key_table <- build_key_table(
  timepoint_table = timepoint_table,
  N = N,
  c1 = C1,
  c2 = C2,
  shared_sse = knot_fit$SSE
)

write.csv(
  key_table,
  file.path(OUT_ROOT, "Table_Key_Results.csv"),
  row.names = FALSE
)

write.csv(
  timepoint_table,
  file.path(OUT_ROOT, "Table_Timepoints.csv"),
  row.names = FALSE
)

cutoff_table <- bind_rows(pair_scans) %>%
  arrange(comparison, k)

write.csv(
  cutoff_table,
  file.path(
    OUT_ROOT,
    "Table_Cutoff_Optimization.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Figures.
# -------------------------------------------------------------------------

overall_df <- build_overall_figure_data(
  group_results = group_results,
  knot_fit = knot_fit
)

figure_paths <- character(0)

overall_path <- file.path(FIG_DIR, "Figure_Overall.png")

make_overall_figure(
  overall_df = overall_df,
  c1 = C1,
  c2 = C2,
  out_file = overall_path
)

figure_paths <- c(figure_paths, overall_path)

for (comparison_name in names(COMPARISONS)) {
  fig_path <- file.path(
    FIG_DIR,
    paste0("Figure_", comparison_name, ".png")
  )

  make_timepoint_figure(
    comparison_name = comparison_name,
    mapping = COMPARISONS[[comparison_name]],
    group_results = group_results,
    pair_scan = pair_scans[[comparison_name]],
    pair_optimum = pair_optima[[comparison_name]],
    c1 = C1,
    c2 = C2,
    out_file = fig_path
  )

  figure_paths <- c(figure_paths, fig_path)
}

# -------------------------------------------------------------------------
# Zip figures while retaining individual PNGs.
# -------------------------------------------------------------------------

ZIP_PATH <- file.path(
  OUT_ROOT,
  "Figures_All.zip"
)

if (
  file.exists(
    ZIP_PATH
  )
) {
  unlink(
    ZIP_PATH
  )
}

old_wd <- getwd()
zip_ok <- FALSE

tryCatch(
  {
    setwd(
      FIG_DIR
    )

    utils::zip(
      zipfile =
        ZIP_PATH,

      files =
        basename(
          figure_paths
        )
    )

    zip_ok <- file.exists(
      ZIP_PATH
    )
  },
  finally = {
    setwd(
      old_wd
    )
  }
)

if (
  !zip_ok
) {
  warning(
    "Figure PNGs were created, but Figures_All.zip was not created."
  )
}

# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

message(
  "============================================================"
)

message(
  "FINAL WEIGHTED-PARETO PC1-NB CUTOFF ANALYSIS COMPLETE"
)

message(
  "Shared c1 = ",
  C1,
  " | shared c2 = ",
  C2
)

message(
  "Remainder: rank < ",
  C1,
  " | Divergence: ",
  C1,
  "-",
  C2,
  " | Leading edge: rank > ",
  C2
)

message(
  "Candidate k domain: 1-",
  MAX_CANDIDATE_K
)


message(
  "Comparison-specific weighted k*: ",
  paste(
    names(
      pair_optima
    ),
    vapply(
      pair_optima,
      function(z) {
        z$selected_k
      },
      integer(1)
    ),
    sep = "=",
    collapse = "; "
  )
)


message(
  "Figures: ",
  FIG_DIR
)

message(
  "Selected sites: Selected_LeadingEdge_Sites.csv"
)

message(
  "Excluded sites: Excluded_Remainder_Crossing_Sites.csv"
)

message(
  "Tables: Table_Key_Results.csv; Table_Timepoints.csv; ",
  "Table_Cutoff_Optimization.csv"
)

message(
  "============================================================"
)
