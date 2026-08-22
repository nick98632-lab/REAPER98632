#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# FINAL MANUSCRIPT ANALYSIS
# PC1-NB LEADING-EDGE REGIME + DATA-DRIVEN TOP-k EIGENVECTOR SPLITTING
# =============================================================================
#
# CORE LOGIC
# ----------
#
# 1. Rank features independently within each experimental arm by ascending
#    absolute PC1 loading:
#
#       rank_i = rank(|v_i1|), ascending
#
#    Thus larger rank = stronger PC1 contribution.
#
#
# 2. Define two empirical variance quantities explicitly.
#
#    RAW-count empirical variance, within one arm:
#
#       s_raw,ig^2 = Var_j(Y_ij | arm g)
#
#    This quantity is used ONLY for the familiar raw-variance geometry and
#    Anchor/Terminal curvature markers.
#
#
#    Pooled within-group empirical variance of DESeq2-normalized counts:
#
#                           sum_g sum_{j in g} (y_ij - ybar_ig)^2
#       V_pool,i =          -------------------------------------
#                                    sum_g (n_g - 1)
#
#    This is an empirical pooled within-group variance; it is NOT a DESeq2
#    fitted dispersion parameter.
#
#    For arm g:
#
#       mu_ig = mean_j(y_ij | arm g)
#       E_ig  = max(V_pool,i - mu_ig, 0)
#
#    E_ig is the excess-over-Poisson variance signal used for PC1-NB
#    cumulative divergence.
#
#
# 3. PC1 variance contribution and cumulative variance-mass divergence:
#
#       P_i = lambda_1 * v_i1^2
#
#       p_g(r) = P_g(r) / sum P_g
#       q_g(r) = E_g(r) / sum E_g
#
#       F_P,g(r) = cumulative sum of p_g(r)
#       F_E,g(r) = cumulative sum of q_g(r)
#
#       D_g(r) = F_E,g(r) - F_P,g(r)
#
#    A shared two-knot continuous linear-spline model is fit jointly to all
#    eight arm-specific D_g(r) curves.  The second shared knot c2 marks the
#    onset of the broad data-derived leading-edge regime:
#
#       LE_g = { feature i : rank_g(i) > c2 }
#
#
# 4. Raw-count Anchor / Terminal geometry:
#
#       y_g(r) = smooth{ log[1 + s_raw,g^2(r)] }
#
#    The second derivative y_g''(r) is calculated from this raw-count variance
#    spline only.
#
#       Anchor_g   = first y_g''(r)=0 sign-change crossing after c2
#       Terminal_g = second successive y_g''(r)=0 sign-change crossing
#
#
# 5. DATA-DRIVEN TOP-k CUTOFF FOR EIGENVECTOR SPLITTING
#
#    For each control/treatment pair and candidate k:
#
#       S_C(k) = top-k features in the control PC1 ranking
#       S_T(k) = top-k features in the treatment PC1 ranking
#       U(k)   = S_C(k) union S_T(k)
#
#    Joint and disjoint status is still defined by top-k membership:
#
#       Joint              = S_C(k) intersection S_T(k)
#       Disjoint control   = S_C(k) \ S_T(k)
#       Disjoint treatment = S_T(k) \ S_C(k)
#
#    A candidate k is ELIGIBLE only when EVERY feature admitted by either
#    top-k list remains inside the broader c2-defined leading-edge regime in
#    BOTH condition-specific rankings:
#
#       U(k) subset of LE_C intersection LE_T
#
#    Equivalently, for every i in U(k):
#
#       rank_C(i) > c2  AND  rank_T(i) > c2
#
#    This preserves biologically meaningful disjoint sites: a site can fail the
#    stricter top-k threshold in the opposite arm while still remaining inside
#    that arm's broader leading-edge regime.
#
#    Eligibility is monotone in k.  Once a newly admitted feature falls at or
#    below c2 in either arm, that k and all larger k are ineligible.
#
#    For each comparison:
#
#       k_pair* = largest eligible k
#
#    To retain ONE common manuscript cutoff across all four comparisons:
#
#       k_common* = min(k_pair*)
#
#    k_common* is therefore the largest single top-k depth that satisfies the
#    leading-edge eligibility rule in every control/treatment comparison.
#
#
# 6. Historical top-5,000 cutoff:
#
#       paper_ref_rank = N - 5000 + 1
#
#    The 5,000 cutoff is retained only as a reference.  It does NOT determine
#    c1, c2, Anchor, Terminal, k_pair*, or k_common*.
#
#
# 7. Post-boundary NB scaling corroboration:
#
#       E = alpha * mu^p
#
#    p is fit after regions are defined:
#
#       p near 1 -> more NB1-like
#       p near 2 -> more NB2-like
#
#    This is corroboration only and does not select any boundary or cutoff.
#
#
# OUTPUTS
# -------
# Figures:
#   Figure_Overall.png
#   Figure_RT0_ZT6.png
#   Figure_RT2_ZT8.png
#   Figure_RT4_ZT10.png
#   Figure_RT8_ZT14.png
#
# Tables:
#   Table_Key_Results.csv
#   Table_Timepoints.csv
#
# Downstream selected sites:
#   Selected_LeadingEdge_Sites.csv
#
# Figure archive:
#   Figures_All.zip
#
# =============================================================================


# =============================================================================
# SETTINGS
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"

OUT_ROOT <- "/root/REAPER98632/exports/pc1_nb_cutoff_final"
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

PAPER_REFERENCE_K <- 5000L

# Original raw-variance derivative smoothing.
ANALYSIS_SPLINE_SPAR <- 0.60

# Display-only smoothing. These values never determine boundaries/cutoffs.
DISPLAY_VAR_SPAR <- 0.72
DISPLAY_D_SPAR   <- 0.72

PNG_WIDTH_IN <- 15
PNG_DPI      <- 360


# =============================================================================
# COLORS
# =============================================================================

COL <- list(
  control = "#386CB0",
  treatment = "#159D91",
  raw = "#117A65",
  divergence = "#262626",
  fit = "#000000",

  remainder = "#DCEFF2",
  interval = "#F5E8C8",
  leading = "#DDF2EA",

  c1 = "#2166AC",
  c2 = "#1B7837",
  selected = "#7B3294",
  paper = "#E69F00",

  anchor_control = "#386CB0",
  terminal_control = "#386CB0",
  anchor_treatment = "#159D91",
  terminal_treatment = "#159D91"
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

save_panels <- function(plots, path, height_in) {
  grDevices::png(
    filename = path,
    width = PNG_WIDTH_IN,
    height = height_in,
    units = "in",
    res = PNG_DPI,
    bg = "white"
  )

  grid.newpage()

  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = length(plots),
        ncol = 1L
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

  dev.off()
}

safe_summary <- function(x, fun = median) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  fun(x)
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

find_d2_zero_crossings <- function(x, d2) {
  ok <- is.finite(x) & is.finite(d2)

  x <- x[ok]
  d2 <- d2[ok]

  if (length(x) < 2L) {
    return(numeric(0))
  }

  out <- numeric(0)

  for (i in seq_len(length(x) - 1L)) {
    a <- d2[i]
    b <- d2[i + 1L]

    if (a == 0 || b == 0) next

    if (
      (a < 0 && b > 0) ||
      (a > 0 && b < 0)
    ) {
      frac <- abs(a) / (abs(a) + abs(b))

      out <- c(
        out,
        x[i] + frac * (x[i + 1L] - x[i])
      )
    }
  }

  sort(unique(out))
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

  analysis_spline <- stats::smooth.spline(
    x = rank,
    y = log_var,
    spar = ANALYSIS_SPLINE_SPAR
  )

  dense_x <- seq(
    min(rank),
    max(rank),
    length.out = max(
      5000L,
      length(rank) * 4L
    )
  )

  dense_d2 <- as.numeric(
    stats::predict(
      analysis_spline,
      x = dense_x,
      deriv = 2
    )$y
  )

  crossings <- find_d2_zero_crossings(
    dense_x,
    dense_d2
  )

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

  list(
    curve = data.frame(
      rank = rank,
      raw_empirical_variance = raw_var_ranked,
      log1p_raw_empirical_variance = log_var,
      display_log1p_raw_empirical_variance = display_y,
      stringsAsFactors = FALSE
    ),
    crossings = crossings
  )
}

select_anchor_terminal_after_c2 <- function(
    crossings,
    c2,
    total_n) {

  z <- sort(
    unique(
      crossings[
        is.finite(crossings) &
        crossings > c2 &
        crossings < total_n
      ]
    )
  )

  if (length(z) < 2L) {
    stop(
      "Fewer than two raw-variance second-derivative zero crossings ",
      "were found after c2 = ",
      c2,
      "."
    )
  }

  anchor <- as.integer(round(z[1L]))
  terminal <- as.integer(round(z[2L]))

  anchor <- max(c2 + 1L, min(total_n - 1L, anchor))
  terminal <- max(anchor + 1L, min(total_n, terminal))

  list(
    anchor = anchor,
    terminal = terminal
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

  df <- geometry$curve %>%
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
    raw_variance_crossings = geometry$crossings,
    anchor = NA_integer_,
    terminal = NA_integer_
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
# NB MEAN-VARIANCE SCALING
# =============================================================================

estimate_nb_exponent <- function(mu, E) {
  keep <- (
    is.finite(mu) &
    is.finite(E) &
    mu > 0 &
    E > 0
  )

  if (sum(keep) < 10L) {
    return(NA_real_)
  }

  fit <- stats::lm(
    log(E[keep]) ~ log(mu[keep])
  )

  unname(
    stats::coef(fit)[2L]
  )
}

get_region_p <- function(df, keep) {
  estimate_nb_exponent(
    mu = df$normalized_group_mean[keep],
    E = df$nb_excess_variance[keep]
  )
}


# =============================================================================
# TOP-k ELIGIBILITY AND JOINT / DISJOINT CLASSIFICATION
# =============================================================================

make_rank_map <- function(df) {
  stats::setNames(
    df$rank,
    df$feature_id
  )
}

find_pairwise_max_eligible_k <- function(
    control_df,
    treatment_df,
    c2) {

  if (nrow(control_df) != nrow(treatment_df)) {
    stop("Control and treatment rankings have different feature counts.")
  }

  if (!setequal(control_df$feature_id, treatment_df$feature_id)) {
    stop("Control and treatment rankings do not contain the same feature IDs.")
  }

  N <- nrow(control_df)

  max_leading_edge_k <- N - c2

  if (max_leading_edge_k < 1L) {
    stop("No features exist beyond c2.")
  }

  control_desc <- rev(control_df$feature_id)
  treatment_desc <- rev(treatment_df$feature_id)

  rank_control <- make_rank_map(control_df)
  rank_treatment <- make_rank_map(treatment_df)

  first_invalid_k <- NA_integer_
  failure_records <- list()

  for (k in seq_len(max_leading_edge_k)) {
    new_ids <- unique(
      c(
        control_desc[k],
        treatment_desc[k]
      )
    )

    rc <- unname(rank_control[new_ids])
    rt <- unname(rank_treatment[new_ids])

    bad <- (
      !is.finite(rc) |
      !is.finite(rt) |
      rc <= c2 |
      rt <= c2
    )

    if (any(bad)) {
      first_invalid_k <- k

      bad_ids <- new_ids[bad]

      failure_records[[1L]] <- data.frame(
        feature_id = bad_ids,
        control_rank = unname(rank_control[bad_ids]),
        treatment_rank = unname(rank_treatment[bad_ids]),
        fails_control_leading_edge =
          unname(rank_control[bad_ids]) <= c2,
        fails_treatment_leading_edge =
          unname(rank_treatment[bad_ids]) <= c2,
        stringsAsFactors = FALSE
      )

      break
    }
  }

  if (is.na(first_invalid_k)) {
    k_star <- max_leading_edge_k
    first_invalid_k <- max_leading_edge_k + 1L
  } else {
    k_star <- first_invalid_k - 1L
  }

  if (k_star < 1L) {
    stop(
      "No positive top-k cutoff satisfies the cross-condition leading-edge ",
      "eligibility rule."
    )
  }

  failures <- if (length(failure_records) == 0L) {
    data.frame(
      feature_id = character(0),
      control_rank = integer(0),
      treatment_rank = integer(0),
      fails_control_leading_edge = logical(0),
      fails_treatment_leading_edge = logical(0),
      stringsAsFactors = FALSE
    )
  } else {
    bind_rows(failure_records)
  }

  list(
    k_star = as.integer(k_star),
    cutoff_rank = rank_cutoff_from_k(N, k_star),
    first_invalid_k = as.integer(first_invalid_k),
    max_leading_edge_k = as.integer(max_leading_edge_k),
    failures = failures
  )
}

classify_pair_at_k <- function(
    control_df,
    treatment_df,
    k,
    c2,
    comparison_name,
    control_group,
    treatment_group) {

  N <- nrow(control_df)

  if (k < 1L || k > N) {
    stop("Invalid k for ", comparison_name, ": ", k)
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

  rank_control <- make_rank_map(control_df)
  rank_treatment <- make_rank_map(treatment_df)

  in_control <- union_ids %in% control_top
  in_treatment <- union_ids %in% treatment_top

  class <- ifelse(
    in_control & in_treatment,
    "Joint",
    ifelse(
      in_control,
      paste0("Disjoint_", control_group),
      paste0("Disjoint_", treatment_group)
    )
  )

  out <- data.frame(
    comparison = comparison_name,
    feature_id = union_ids,
    class = class,
    control_group = control_group,
    treatment_group = treatment_group,
    selected_k = k,
    control_rank = unname(rank_control[union_ids]),
    treatment_rank = unname(rank_treatment[union_ids]),
    control_top_k = in_control,
    treatment_top_k = in_treatment,
    control_in_c2_leading_edge =
      unname(rank_control[union_ids]) > c2,
    treatment_in_c2_leading_edge =
      unname(rank_treatment[union_ids]) > c2,
    stringsAsFactors = FALSE
  )

  out$eligible_both_c2 <- (
    out$control_in_c2_leading_edge &
    out$treatment_in_c2_leading_edge
  )

  out
}

summarize_classification <- function(class_df) {
  class_counts <- table(class_df$class)

  joint_n <- if ("Joint" %in% names(class_counts)) {
    unname(class_counts[["Joint"]])
  } else {
    0L
  }

  disjoint_n <- nrow(class_df) - joint_n

  control_only_n <- sum(
    class_df$control_top_k &
    !class_df$treatment_top_k
  )

  treatment_only_n <- sum(
    !class_df$control_top_k &
    class_df$treatment_top_k
  )

  union_n <- nrow(class_df)

  k <- unique(class_df$selected_k)

  jaccard <- if (union_n > 0L) {
    joint_n / union_n
  } else {
    NA_real_
  }

  data.frame(
    k = k,
    joint_n = joint_n,
    disjoint_control_n = control_only_n,
    disjoint_treatment_n = treatment_only_n,
    disjoint_total_n = disjoint_n,
    union_n = union_n,
    jaccard = jaccard,
    all_union_sites_inside_both_c2_regimes =
      all(class_df$eligible_both_c2),
    offending_n =
      sum(!class_df$eligible_both_c2),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# RESULTS TABLES
# =============================================================================

build_arm_scaling_table <- function(
    group_results,
    comparisons,
    c1,
    c2,
    common_k,
    paper_k) {

  rows <- list()

  for (comparison_name in names(comparisons)) {
    mapping <- comparisons[[comparison_name]]

    for (arm in names(mapping)) {
      g <- unname(mapping[[arm]])
      df <- group_results[[g]]$data
      N <- nrow(df)

      selected_rank <- rank_cutoff_from_k(
        N,
        common_k
      )

      paper_rank <- rank_cutoff_from_k(
        N,
        paper_k
      )

      rows[[length(rows) + 1L]] <- data.frame(
        comparison = comparison_name,
        arm = arm,
        group = g,
        p_remainder = get_region_p(
          df,
          df$rank < c1
        ),
        p_leading_edge = get_region_p(
          df,
          df$rank > c2
        ),
        p_selected_k = get_region_p(
          df,
          df$rank >= selected_rank
        ),
        p_paper_5000 = get_region_p(
          df,
          df$rank >= paper_rank
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  bind_rows(rows)
}

build_timepoint_table <- function(
    group_results,
    comparisons,
    pair_cutoffs,
    selected_sites,
    paper_summaries,
    common_k,
    c1,
    c2,
    paper_k,
    N) {

  rows <- list()

  for (comparison_name in names(comparisons)) {
    mapping <- comparisons[[comparison_name]]

    control_group <- unname(mapping[["control"]])
    treatment_group <- unname(mapping[["treatment"]])

    pair_info <- pair_cutoffs[[comparison_name]]

    selected_df <- selected_sites %>%
      filter(comparison == comparison_name)

    selected_summary <- summarize_classification(
      selected_df
    )

    paper_summary <- paper_summaries[[comparison_name]]

    rows[[length(rows) + 1L]] <- data.frame(
      comparison = comparison_name,
      control_group = control_group,
      treatment_group = treatment_group,

      c1 = c1,
      c2 = c2,

      control_anchor =
        group_results[[control_group]]$anchor,
      control_terminal =
        group_results[[control_group]]$terminal,
      treatment_anchor =
        group_results[[treatment_group]]$anchor,
      treatment_terminal =
        group_results[[treatment_group]]$terminal,

      pairwise_max_eligible_k =
        pair_info$k_star,
      pairwise_cutoff_rank =
        pair_info$cutoff_rank,
      first_invalid_k =
        pair_info$first_invalid_k,

      common_selected_k = common_k,
      common_selected_cutoff_rank =
        rank_cutoff_from_k(N, common_k),

      paper_reference_k = paper_k,
      paper_reference_cutoff_rank =
        rank_cutoff_from_k(N, paper_k),
      paper_5000_eligible =
        paper_k <= pair_info$k_star,

      selected_joint_n =
        selected_summary$joint_n,
      selected_disjoint_control_n =
        selected_summary$disjoint_control_n,
      selected_disjoint_treatment_n =
        selected_summary$disjoint_treatment_n,
      selected_union_n =
        selected_summary$union_n,
      selected_jaccard =
        selected_summary$jaccard,

      paper_joint_n =
        paper_summary$joint_n,
      paper_disjoint_control_n =
        paper_summary$disjoint_control_n,
      paper_disjoint_treatment_n =
        paper_summary$disjoint_treatment_n,
      paper_union_n =
        paper_summary$union_n,
      paper_offending_n =
        paper_summary$offending_n,

      stringsAsFactors = FALSE
    )
  }

  bind_rows(rows)
}

build_key_table <- function(
    timepoint_table,
    arm_scaling,
    N,
    c1,
    c2,
    common_k,
    paper_k,
    shared_sse) {

  data.frame(
    n_features = N,

    shared_c1 = c1,
    shared_c2 = c2,
    c2_leading_edge_n = N - c2,

    common_selected_k = common_k,
    common_selected_cutoff_rank =
      rank_cutoff_from_k(N, common_k),

    pairwise_k_min =
      min(timepoint_table$pairwise_max_eligible_k),
    pairwise_k_median =
      median(timepoint_table$pairwise_max_eligible_k),
    pairwise_k_max =
      max(timepoint_table$pairwise_max_eligible_k),

    paper_reference_k = paper_k,
    paper_reference_cutoff_rank =
      rank_cutoff_from_k(N, paper_k),
    paper_reference_used_in_analysis = FALSE,
    paper_reference_eligible_all_comparisons =
      all(timepoint_table$paper_5000_eligible),

    anchor_median =
      median(
        c(
          timepoint_table$control_anchor,
          timepoint_table$treatment_anchor
        )
      ),

    terminal_median =
      median(
        c(
          timepoint_table$control_terminal,
          timepoint_table$treatment_terminal
        )
      ),

    p_remainder_median =
      safe_summary(arm_scaling$p_remainder),

    p_leading_edge_median =
      safe_summary(arm_scaling$p_leading_edge),

    p_selected_k_median =
      safe_summary(arm_scaling$p_selected_k),

    p_paper_5000_median =
      safe_summary(arm_scaling$p_paper_5000),

    shared_model_SSE = shared_sse,

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# FIGURE DATA
# =============================================================================

rankwise_median <- function(group_results, column) {
  mat <- do.call(
    cbind,
    lapply(
      group_results,
      function(z) z$data[[column]]
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
      alpha = 0.48
    ) +
    annotate(
      "rect",
      xmin = c1,
      xmax = c2,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$interval,
      alpha = 0.42
    ) +
    annotate(
      "rect",
      xmin = c2,
      xmax = N,
      ymin = -Inf,
      ymax = Inf,
      fill = COL$leading,
      alpha = 0.50
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


# =============================================================================
# OVERALL FIGURE
# =============================================================================

make_overall_figure <- function(
    overall_df,
    key_table,
    arm_scaling,
    c1,
    c2,
    common_k,
    paper_k,
    out_file) {

  N <- nrow(overall_df)

  selected_rank <- rank_cutoff_from_k(
    N,
    common_k
  )

  paper_rank <- rank_cutoff_from_k(
    N,
    paper_k
  )

  anchor_median <- key_table$anchor_median[1L]
  terminal_median <- key_table$terminal_median[1L]

  paper_status <- if (
    key_table$paper_reference_eligible_all_comparisons[1L]
  ) {
    "eligible"
  } else {
    "ineligible"
  }

  # -----------------------------------------------------------------------
  # A. RAW-COUNT VARIANCE GEOMETRY
  # -----------------------------------------------------------------------

  a_lines <- data.frame(
    rank = c(
      c2,
      selected_rank,
      paper_rank
    ),
    key = c(
      "c2",
      paste0("k*=", common_k),
      "5k ref"
    ),
    stringsAsFactors = FALSE
  )

  a_colors <- c(
    "Raw variance" = COL$raw,
    "c2" = COL$c2,
    stats::setNames(
      COL$selected,
      paste0("k*=", common_k)
    ),
    "5k ref" = COL$paper
  )

  a_types <- c(
    "Raw variance" = "solid",
    "c2" = "longdash",
    stats::setNames(
      "dotdash",
      paste0("k*=", common_k)
    ),
    "5k ref" = "dotted"
  )

  pA <- ggplot()
  pA <- add_rank_regions(
    pA,
    c1,
    c2,
    N
  )

  pA <- pA +
    geom_vline(
      data = a_lines,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      linewidth = 0.85
    ) +
    geom_vline(
      xintercept = anchor_median,
      color = COL$divergence,
      linewidth = 0.70,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = terminal_median,
      color = COL$divergence,
      linetype = "dotted",
      linewidth = 0.75,
      show.legend = FALSE
    ) +
    geom_line(
      data = overall_df,
      aes(
        x = rank,
        y = raw_variance,
        color = "Raw variance",
        linetype = "Raw variance"
      ),
      linewidth = 1.20,
      lineend = "round"
    ) +
    annotate(
      "text",
      x = anchor_median,
      y = -Inf,
      label = "A",
      angle = 90,
      vjust = -0.35,
      size = 2.8
    ) +
    annotate(
      "text",
      x = terminal_median,
      y = -Inf,
      label = "T",
      angle = 90,
      vjust = -0.35,
      size = 2.8
    ) +
    scale_color_manual(
      values = a_colors,
      breaks = names(a_colors)
    ) +
    scale_linetype_manual(
      values = a_types,
      breaks = names(a_types)
    ) +
    labs(
      title = "A. Raw-count variance geometry",
      subtitle = paste0(
        "A/T = first two y''=0 crossings after c2; 5k reference = ",
        paper_status
      ),
      x = "PC1 rank: low |loading| -> high |loading|",
      y = "Smoothed log(1 + raw-count variance)"
    ) +
    theme_manuscript()

  # -----------------------------------------------------------------------
  # B. CUMULATIVE DIVERGENCE
  # -----------------------------------------------------------------------

  b_lines <- data.frame(
    rank = c(
      c1,
      c2,
      selected_rank,
      paper_rank
    ),
    key = c(
      "c1",
      "c2",
      paste0("k*=", common_k),
      "5k ref"
    ),
    stringsAsFactors = FALSE
  )

  b_colors <- c(
    "D(r)" = COL$divergence,
    "Shared fit" = COL$fit,
    "c1" = COL$c1,
    "c2" = COL$c2,
    stats::setNames(
      COL$selected,
      paste0("k*=", common_k)
    ),
    "5k ref" = COL$paper
  )

  b_types <- c(
    "D(r)" = "solid",
    "Shared fit" = "solid",
    "c1" = "dashed",
    "c2" = "longdash",
    stats::setNames(
      "dotdash",
      paste0("k*=", common_k)
    ),
    "5k ref" = "dotted"
  )

  pB <- ggplot()
  pB <- add_rank_regions(
    pB,
    c1,
    c2,
    N
  )

  pB <- pB +
    geom_vline(
      data = b_lines,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
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
      aes(
        x = rank,
        y = D,
        color = "D(r)",
        linetype = "D(r)"
      ),
      linewidth = 0.95,
      alpha = 0.70
    ) +
    geom_line(
      data = overall_df,
      aes(
        x = rank,
        y = D_fit,
        color = "Shared fit",
        linetype = "Shared fit"
      ),
      linewidth = 1.35
    ) +
    scale_color_manual(
      values = b_colors,
      breaks = names(b_colors)
    ) +
    scale_linetype_manual(
      values = b_types,
      breaks = names(b_types)
    ) +
    labs(
      title = "B. Cumulative PC1-NB variance-mass divergence",
      subtitle = "D(r)=F_E(r)-F_P(r); shared c1/c2 fitted across all 8 arms",
      x = "PC1 rank",
      y = "Cumulative divergence D(r)"
    ) +
    theme_manuscript()

  # -----------------------------------------------------------------------
  # C. POST-BOUNDARY NB SCALING
  # -----------------------------------------------------------------------

  scaling_plot <- arm_scaling %>%
    select(
      group,
      p_remainder,
      p_leading_edge,
      p_selected_k,
      p_paper_5000
    ) %>%
    pivot_longer(
      cols = c(
        p_remainder,
        p_leading_edge,
        p_selected_k,
        p_paper_5000
      ),
      names_to = "region",
      values_to = "p"
    ) %>%
    mutate(
      region = factor(
        region,
        levels = c(
          "p_remainder",
          "p_leading_edge",
          "p_selected_k",
          "p_paper_5000"
        ),
        labels = c(
          "Remainder",
          "Leading edge",
          paste0("Selected k*=", common_k),
          "5k reference"
        )
      )
    )

  pC <- ggplot(
    scaling_plot,
    aes(
      x = p,
      y = region
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "#777777",
      linewidth = 0.55
    ) +
    geom_vline(
      xintercept = 2,
      linetype = "dotted",
      color = "#777777",
      linewidth = 0.65
    ) +
    geom_point(
      alpha = 0.72,
      size = 2.4,
      position = position_jitter(
        height = 0.08,
        width = 0
      )
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 4.0,
      shape = 18
    ) +
    labs(
      title = "C. Post-boundary NB mean-variance scaling",
      subtitle = "Corroboration only: p~1 NB1-like; p~2 NB2-like",
      x = "Exponent p in E = alpha * mu^p",
      y = NULL
    ) +
    theme_manuscript() +
    theme(
      legend.position = "none"
    )

  save_panels(
    list(
      pA,
      pB,
      pC
    ),
    out_file,
    height_in = 11.7
  )
}


# =============================================================================
# TIME-POINT FIGURES
# =============================================================================

make_timepoint_figure <- function(
    comparison_name,
    mapping,
    group_results,
    arm_scaling,
    timepoint_table,
    c1,
    c2,
    common_k,
    paper_k,
    out_file) {

  control_group <- unname(
    mapping[["control"]]
  )

  treatment_group <- unname(
    mapping[["treatment"]]
  )

  control <- group_results[[control_group]]$data
  treatment <- group_results[[treatment_group]]$data

  N <- nrow(control)

  selected_rank <- rank_cutoff_from_k(
    N,
    common_k
  )

  paper_rank <- rank_cutoff_from_k(
    N,
    paper_k
  )

  tp <- timepoint_table %>%
    filter(
      comparison == comparison_name
    )

  pair_k <- tp$pairwise_max_eligible_k[1L]

  paper_status <- if (
    tp$paper_5000_eligible[1L]
  ) {
    "eligible"
  } else {
    "ineligible"
  }

  c_anchor <- group_results[[control_group]]$anchor
  c_terminal <- group_results[[control_group]]$terminal

  t_anchor <- group_results[[treatment_group]]$anchor
  t_terminal <- group_results[[treatment_group]]$terminal

  # -----------------------------------------------------------------------
  # A. RAW-COUNT VARIANCE GEOMETRY
  # -----------------------------------------------------------------------

  a_lines <- data.frame(
    rank = c(
      c2,
      selected_rank,
      paper_rank
    ),
    key = c(
      "c2",
      paste0("k*=", common_k),
      "5k ref"
    ),
    stringsAsFactors = FALSE
  )

  a_colors <- c(
    "Control" = COL$control,
    "Treatment" = COL$treatment,
    "c2" = COL$c2,
    stats::setNames(
      COL$selected,
      paste0("k*=", common_k)
    ),
    "5k ref" = COL$paper
  )

  a_types <- c(
    "Control" = "solid",
    "Treatment" = "solid",
    "c2" = "longdash",
    stats::setNames(
      "dotdash",
      paste0("k*=", common_k)
    ),
    "5k ref" = "dotted"
  )

  pA <- ggplot()
  pA <- add_rank_regions(
    pA,
    c1,
    c2,
    N
  )

  pA <- pA +
    geom_vline(
      data = a_lines,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
      linewidth = 0.84
    ) +
    geom_vline(
      xintercept = c_anchor,
      color = COL$anchor_control,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = c_terminal,
      color = COL$terminal_control,
      linetype = "dotted",
      linewidth = 0.75,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = t_anchor,
      color = COL$anchor_treatment,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = t_terminal,
      color = COL$terminal_treatment,
      linetype = "dotted",
      linewidth = 0.75,
      show.legend = FALSE
    ) +
    geom_line(
      data = control,
      aes(
        x = rank,
        y = display_log1p_raw_empirical_variance,
        color = "Control",
        linetype = "Control"
      ),
      linewidth = 1.12
    ) +
    geom_line(
      data = treatment,
      aes(
        x = rank,
        y = display_log1p_raw_empirical_variance,
        color = "Treatment",
        linetype = "Treatment"
      ),
      linewidth = 1.12
    ) +
    annotate(
      "text",
      x = c_anchor,
      y = -Inf,
      label = "A_C",
      angle = 90,
      vjust = -0.25,
      size = 2.5,
      color = COL$control
    ) +
    annotate(
      "text",
      x = c_terminal,
      y = -Inf,
      label = "T_C",
      angle = 90,
      vjust = -0.25,
      size = 2.5,
      color = COL$control
    ) +
    annotate(
      "text",
      x = t_anchor,
      y = -Inf,
      label = "A_T",
      angle = 90,
      vjust = -0.25,
      size = 2.5,
      color = COL$treatment
    ) +
    annotate(
      "text",
      x = t_terminal,
      y = -Inf,
      label = "T_T",
      angle = 90,
      vjust = -0.25,
      size = 2.5,
      color = COL$treatment
    ) +
    scale_color_manual(
      values = a_colors,
      breaks = names(a_colors)
    ) +
    scale_linetype_manual(
      values = a_types,
      breaks = names(a_types)
    ) +
    labs(
      title = paste0(
        "A. ",
        comparison_name,
        " raw-count variance geometry"
      ),
      subtitle = paste0(
        "Pair k*=", pair_k,
        "; common k*=", common_k,
        "; 5k=", paper_status,
        "; A/T = first two y''=0 crossings after c2"
      ),
      x = "PC1 rank: low |loading| -> high |loading|",
      y = "Smoothed log(1 + raw-count variance)"
    ) +
    theme_manuscript()

  # -----------------------------------------------------------------------
  # B. CUMULATIVE DIVERGENCE
  # -----------------------------------------------------------------------

  b_lines <- data.frame(
    rank = c(
      c1,
      c2,
      selected_rank,
      paper_rank
    ),
    key = c(
      "c1",
      "c2",
      paste0("k*=", common_k),
      "5k ref"
    ),
    stringsAsFactors = FALSE
  )

  b_colors <- c(
    "Control D(r)" = COL$control,
    "Treatment D(r)" = COL$treatment,
    "c1" = COL$c1,
    "c2" = COL$c2,
    stats::setNames(
      COL$selected,
      paste0("k*=", common_k)
    ),
    "5k ref" = COL$paper
  )

  b_types <- c(
    "Control D(r)" = "solid",
    "Treatment D(r)" = "solid",
    "c1" = "dashed",
    "c2" = "longdash",
    stats::setNames(
      "dotdash",
      paste0("k*=", common_k)
    ),
    "5k ref" = "dotted"
  )

  pB <- ggplot()
  pB <- add_rank_regions(
    pB,
    c1,
    c2,
    N
  )

  pB <- pB +
    geom_vline(
      data = b_lines,
      aes(
        xintercept = rank,
        color = key,
        linetype = key
      ),
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
      aes(
        x = rank,
        y = display_D,
        color = "Control D(r)",
        linetype = "Control D(r)"
      ),
      linewidth = 1.08
    ) +
    geom_line(
      data = treatment,
      aes(
        x = rank,
        y = display_D,
        color = "Treatment D(r)",
        linetype = "Treatment D(r)"
      ),
      linewidth = 1.08
    ) +
    scale_color_manual(
      values = b_colors,
      breaks = names(b_colors)
    ) +
    scale_linetype_manual(
      values = b_types,
      breaks = names(b_types)
    ) +
    labs(
      title = paste0(
        "B. ",
        comparison_name,
        " cumulative PC1-NB divergence"
      ),
      subtitle = "D(r)=F_E(r)-F_P(r); shared c1/c2 from all 8 arms",
      x = "PC1 rank",
      y = "Cumulative divergence D(r)"
    ) +
    theme_manuscript()

  # -----------------------------------------------------------------------
  # C. POST-BOUNDARY NB SCALING
  # -----------------------------------------------------------------------

  scaling <- arm_scaling %>%
    filter(
      comparison == comparison_name
    ) %>%
    mutate(
      arm_label = ifelse(
        arm == "control",
        "Control",
        "Treatment"
      )
    ) %>%
    select(
      arm_label,
      p_remainder,
      p_leading_edge,
      p_selected_k,
      p_paper_5000
    ) %>%
    pivot_longer(
      cols = c(
        p_remainder,
        p_leading_edge,
        p_selected_k,
        p_paper_5000
      ),
      names_to = "region",
      values_to = "p"
    ) %>%
    mutate(
      region = factor(
        region,
        levels = c(
          "p_remainder",
          "p_leading_edge",
          "p_selected_k",
          "p_paper_5000"
        ),
        labels = c(
          "Remainder",
          "Leading edge",
          paste0("Selected k*=", common_k),
          "5k reference"
        )
      )
    )

  pC <- ggplot(
    scaling,
    aes(
      x = p,
      y = region,
      color = arm_label
    )
  ) +
    geom_vline(
      xintercept = 1,
      color = "#777777",
      linetype = "dashed",
      linewidth = 0.55
    ) +
    geom_vline(
      xintercept = 2,
      color = "#777777",
      linetype = "dotted",
      linewidth = 0.65
    ) +
    geom_point(size = 3.5) +
    scale_color_manual(
      values = c(
        "Control" = COL$control,
        "Treatment" = COL$treatment
      )
    ) +
    labs(
      title = paste0(
        "C. ",
        comparison_name,
        " post-boundary NB scaling"
      ),
      subtitle = "Corroboration only: p~1 NB1-like; p~2 NB2-like",
      x = "Exponent p in E = alpha * mu^p",
      y = NULL
    ) +
    theme_manuscript()

  save_panels(
    list(
      pA,
      pB,
      pC
    ),
    out_file,
    height_in = 11.8
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

N <- nrow(count_mat)

if (PAPER_REFERENCE_K >= N) {
  stop("PAPER_REFERENCE_K must be smaller than N.")
}

PAPER_REFERENCE_RANK <- rank_cutoff_from_k(
  N,
  PAPER_REFERENCE_K
)

deseq <- normalize_deseq2_global(
  count_mat = count_mat,
  group_labels = group_labels
)

normalized_counts <- deseq$normalized_counts

pooled <- compute_pooled_within_group_variance(
  normalized_counts = normalized_counts,
  group_labels = group_labels
)

message(
  "Pooled within-group residual df: ",
  pooled$residual_df
)

group_results <- vector(
  "list",
  length(levels(group_labels))
)

names(group_results) <- levels(group_labels)

for (g in levels(group_labels)) {
  idx <- which(group_labels == g)

  if (length(idx) < 2L) {
    stop("Not enough samples in group ", g)
  }

  message("Analyzing ", g, "...")

  group_results[[g]] <- compute_group_analysis(
    group_name = g,
    raw_counts_arm = count_mat[, idx, drop = FALSE],
    normalized_counts_arm = normalized_counts[, idx, drop = FALSE],
    pooled_variance = pooled$variance
  )
}

message("Fitting shared two-knot cumulative-divergence model...")

knot_fit <- fit_shared_knots(
  group_results
)

C1 <- knot_fit$c1
C2 <- knot_fit$c2

message("Shared c1 = ", C1)
message("Shared c2 = ", C2)

message(
  "Selecting Anchor/Terminal from first two raw-variance y''=0 crossings after c2..."
)

for (g in names(group_results)) {
  at <- select_anchor_terminal_after_c2(
    crossings =
      group_results[[g]]$raw_variance_crossings,
    c2 = C2,
    total_n = N
  )

  group_results[[g]]$anchor <- at$anchor
  group_results[[g]]$terminal <- at$terminal
}

message(
  "Finding pairwise maximum eligible top-k values..."
)

pair_cutoffs <- list()

for (comparison_name in names(COMPARISONS)) {
  mapping <- COMPARISONS[[comparison_name]]

  control_group <- unname(
    mapping[["control"]]
  )

  treatment_group <- unname(
    mapping[["treatment"]]
  )

  pair_cutoffs[[comparison_name]] <-
    find_pairwise_max_eligible_k(
      control_df =
        group_results[[control_group]]$data,
      treatment_df =
        group_results[[treatment_group]]$data,
      c2 = C2
    )

  message(
    comparison_name,
    ": pairwise k* = ",
    pair_cutoffs[[comparison_name]]$k_star,
    "; first invalid k = ",
    pair_cutoffs[[comparison_name]]$first_invalid_k
  )
}

PAIR_K_VALUES <- vapply(
  pair_cutoffs,
  function(z) z$k_star,
  integer(1)
)

COMMON_K <- min(
  PAIR_K_VALUES
)

COMMON_CUTOFF_RANK <- rank_cutoff_from_k(
  N,
  COMMON_K
)

message(
  "Largest common eligible k across all comparisons = ",
  COMMON_K,
  " (cutoff rank ",
  COMMON_CUTOFF_RANK,
  ")"
)

message(
  "Historical paper k = ",
  PAPER_REFERENCE_K,
  " (cutoff rank ",
  PAPER_REFERENCE_RANK,
  ")"
)

# -------------------------------------------------------------------------
# Classify selected Joint / Disjoint sites at the common eligible k.
# -------------------------------------------------------------------------

selected_site_rows <- list()
paper_summaries <- list()

for (comparison_name in names(COMPARISONS)) {
  mapping <- COMPARISONS[[comparison_name]]

  control_group <- unname(
    mapping[["control"]]
  )

  treatment_group <- unname(
    mapping[["treatment"]]
  )

  selected_df <- classify_pair_at_k(
    control_df =
      group_results[[control_group]]$data,
    treatment_df =
      group_results[[treatment_group]]$data,
    k = COMMON_K,
    c2 = C2,
    comparison_name = comparison_name,
    control_group = control_group,
    treatment_group = treatment_group
  )

  if (!all(selected_df$eligible_both_c2)) {
    stop(
      "Internal error: common selected k is not eligible for ",
      comparison_name,
      "."
    )
  }

  selected_site_rows[[comparison_name]] <- selected_df

  paper_df <- classify_pair_at_k(
    control_df =
      group_results[[control_group]]$data,
    treatment_df =
      group_results[[treatment_group]]$data,
    k = PAPER_REFERENCE_K,
    c2 = C2,
    comparison_name = comparison_name,
    control_group = control_group,
    treatment_group = treatment_group
  )

  paper_summaries[[comparison_name]] <-
    summarize_classification(
      paper_df
    )
}

selected_sites <- bind_rows(
  selected_site_rows
)

write.csv(
  selected_sites,
  file.path(
    OUT_ROOT,
    "Selected_LeadingEdge_Sites.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Post-boundary NB scaling.
# -------------------------------------------------------------------------

arm_scaling <- build_arm_scaling_table(
  group_results = group_results,
  comparisons = COMPARISONS,
  c1 = C1,
  c2 = C2,
  common_k = COMMON_K,
  paper_k = PAPER_REFERENCE_K
)

# -------------------------------------------------------------------------
# Compact manuscript result tables.
# -------------------------------------------------------------------------

timepoint_table <- build_timepoint_table(
  group_results = group_results,
  comparisons = COMPARISONS,
  pair_cutoffs = pair_cutoffs,
  selected_sites = selected_sites,
  paper_summaries = paper_summaries,
  common_k = COMMON_K,
  c1 = C1,
  c2 = C2,
  paper_k = PAPER_REFERENCE_K,
  N = N
)

key_table <- build_key_table(
  timepoint_table = timepoint_table,
  arm_scaling = arm_scaling,
  N = N,
  c1 = C1,
  c2 = C2,
  common_k = COMMON_K,
  paper_k = PAPER_REFERENCE_K,
  shared_sse = knot_fit$SSE
)

write.csv(
  key_table,
  file.path(
    OUT_ROOT,
    "Table_Key_Results.csv"
  ),
  row.names = FALSE
)

write.csv(
  timepoint_table,
  file.path(
    OUT_ROOT,
    "Table_Timepoints.csv"
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

overall_path <- file.path(
  FIG_DIR,
  "Figure_Overall.png"
)

make_overall_figure(
  overall_df = overall_df,
  key_table = key_table,
  arm_scaling = arm_scaling,
  c1 = C1,
  c2 = C2,
  common_k = COMMON_K,
  paper_k = PAPER_REFERENCE_K,
  out_file = overall_path
)

figure_paths <- c(
  figure_paths,
  overall_path
)

for (comparison_name in names(COMPARISONS)) {
  fig_path <- file.path(
    FIG_DIR,
    paste0(
      "Figure_",
      comparison_name,
      ".png"
    )
  )

  make_timepoint_figure(
    comparison_name = comparison_name,
    mapping = COMPARISONS[[comparison_name]],
    group_results = group_results,
    arm_scaling = arm_scaling,
    timepoint_table = timepoint_table,
    c1 = C1,
    c2 = C2,
    common_k = COMMON_K,
    paper_k = PAPER_REFERENCE_K,
    out_file = fig_path
  )

  figure_paths <- c(
    figure_paths,
    fig_path
  )
}

# -------------------------------------------------------------------------
# Zip all figures while keeping each PNG individually available.
# -------------------------------------------------------------------------

ZIP_PATH <- file.path(
  OUT_ROOT,
  "Figures_All.zip"
)

if (file.exists(ZIP_PATH)) {
  unlink(ZIP_PATH)
}

old_wd <- getwd()
zip_ok <- FALSE

tryCatch(
  {
    setwd(FIG_DIR)

    utils::zip(
      zipfile = ZIP_PATH,
      files = basename(figure_paths)
    )

    zip_ok <- file.exists(ZIP_PATH)
  },
  finally = {
    setwd(old_wd)
  }
)

if (!zip_ok) {
  warning(
    "Figure PNGs were created, but Figures_All.zip was not created."
  )
}


# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

message("============================================================")
message("FINAL PC1-NB LEADING-EDGE ANALYSIS COMPLETE")
message("Shared c1 = ", C1)
message("Shared c2 = ", C2)
message("c2-defined leading-edge size per arm = ", N - C2)

message(
  "Pairwise maximum eligible k*: ",
  paste(
    names(PAIR_K_VALUES),
    PAIR_K_VALUES,
    sep = "=",
    collapse = "; "
  )
)

message(
  "Largest common eligible k* = ",
  COMMON_K,
  " (rank >= ",
  COMMON_CUTOFF_RANK,
  ")"
)

message(
  "Historical 5,000 reference = rank >= ",
  PAPER_REFERENCE_RANK,
  "; eligible in all comparisons = ",
  all(timepoint_table$paper_5000_eligible)
)

message(
  "Anchor-Terminal ranges: ",
  paste(
    names(group_results),
    vapply(
      group_results,
      function(z) {
        paste0(
          z$anchor,
          "-",
          z$terminal
        )
      },
      character(1)
    ),
    sep = "=",
    collapse = "; "
  )
)

message("Figures: ", FIG_DIR)
message("Figure zip: ", ZIP_PATH)
message("Selected sites: Selected_LeadingEdge_Sites.csv")
message("Tables: Table_Key_Results.csv; Table_Timepoints.csv")
message("============================================================")
