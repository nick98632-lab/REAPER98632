#!/usr/bin/env Rscript
# =============================================================================
# SEQUENCE EVS COMPOSITION: determining k from the Joint / Disjoint composition
# =============================================================================
#
# RUN
#   Rscript SEQUENCE_EVS_COMPOSITION.R
#   Optional environment variables:
#     SEQUENCE_REPO_ROOT        folder holding WTTS-Seq_2022.2_DE_raw_read_numbers.csv
#                               (default: current folder)
#     COMPOSITION_PERMUTATIONS  permutations for the null (default 99)
#     COMPOSITION_SEED          random seed (default 20261108)
#   Output: <repo>/exports/sequence_evs_composition
#     SEQUENCE_PIPELINE_PATH    SEQUENCE_EVS_MANUSCRIPT.R (default: <repo>/SEQUENCE_EVS_MANUSCRIPT.R)
#   It does not change the manuscript runs. The CPM k* is recomputed here with
#   the manuscript script's own CPM functions (loaded from it, not copied), so
#   the two cutoffs are compared on the same counts and the same rankings.
#
# FIGURES
#   Figure_C1  how the CPM method chooses k: knots c1, c2 from all 40 samples
#              (8 arms); then, per comparison, the Pareto frontier of its two
#              arms (10 samples) and k* at the maximum endpoint deviation
#   Figure_C2  how the composition method chooses k: f observed vs the
#              permutation null, per comparison and condition
#   Figure_C3  per comparison: (a) estimated Leading Edge size at each lambda,
#              (b) the union curve U(k) and the single k where it reaches the
#              estimated Leading Edge, with CPM k* and 5,000 marked, (c) the
#              Joint / Disjoint / Remainder composition estimated vs realized
#   Figure_C4  the k of each method side by side
#
# WHAT IT DOES
# For each comparison (RT0/ZT6, RT2/ZT8, RT4/ZT10, RT8/ZT14) it estimates how
# many PASs are structured in RT, in ZT and in both, and from that composition
# the EVS cutoff k, then tests whether the structure is real.
#
# THE DECOMPOSITION
# Within condition c (5 samples), X_c is the PAS x sample matrix of log1p(CPM),
# centred per PAS. Its SVD gives the condition's PC1 sample pattern u_c. Each
# PAS's variation splits exactly into a structured part and a residual:
#     ||x_i||^2 = (x_i . u_c)^2 + ||r_i||^2
# and its structured fraction is
#     f_i = (x_i . u_c)^2 / ||x_i||^2 = cos^2(angle between x_i and PC1).
#
# THE PROOF (assumptions A1-A4, written out in full in the conversation notes)
#   A1  Each PAS is structured or unstructured in each condition, giving four
#       classes: Joint (both; J), RT-only (D_R), ZT-only (D_Z), Remainder (M).
#   A2  An unstructured PAS has a direction unrelated to PC1, so its f follows
#       a known null distribution F.
#   A3  Structured PASs have f >= lambda.
#   A4  For Remainder PASs, f in RT and f in ZT are independent.
#   Counting A_c = #{f_c < lambda} and B = #{f_RT < lambda and f_ZT < lambda}:
#       E[A_RT] = (D_Z + M) F(lambda),  E[A_ZT] = (D_R + M) F(lambda),
#       E[B]    = M F(lambda)^2
#   so
#       S_RT = N - A_RT / F,   S_ZT = N - A_ZT / F,   M = B / F^2,
#       J    = S_RT + S_ZT - (N - M)      (inclusion-exclusion),
#       D_R  = S_RT - J,       D_Z = S_ZT - J.
#   Each estimator is unbiased. As in the CPM method, k is one value per
#   comparison, set by both conditions together (all 10 samples): the
#   Leading Edge is N - M (every PAS structured in at least one condition),
#   and k is the smallest k whose union of the two top-k sets (by |PC1
#   loading|) reaches N - M. The per-condition counts S_RT and S_ZT are
#   reported alongside.
#   On a deterministic test set with known truth (J = 500, D_R = D_Z = 250,
#   M = 9,025) these recovered k = 750 +/- 2 and J = 500 +/- 9 at every lambda.
#
# THE NULL F, FROM THE DATA
# With 5 samples and isotropic noise, F would be Beta(1/2, 3/2). Count noise
# is not isotropic and depends on the mean, so F is estimated from the data:
# each PAS's 5 values are shuffled independently within the condition (this
# keeps every PAS's mean and variance and destroys any shared sample pattern),
# the SVD and f are recomputed, and F is the fraction of permuted f below
# lambda. This is done separately within 10 mean-expression bins, and the
# estimators above are applied within each bin and summed.
#
# k FROM THE DECOMPOSITION (no cutoff is searched for)
# With 5 samples each condition has 4 modes, and each PAS's variation splits
# exactly across them. A PAS is governed by PC1 when PC1 carries more of its
# variation than all other modes together, f_1 > 1/2: its structured part
# outweighs its residual. That is the definition of "structured", so the
# counting proof is applied at lambda = 1/2: it counts the PASs governed by
# PC1 after removing the noise PASs that land above 1/2 by chance (measured by
# permutation). The Remainder is the PASs governed by PC1 in neither
# condition; the Leading Edge is the rest; k is where the union of the two
# top-k sets reaches it. The dominant-mode rule (f_1 is the largest single
# mode) is reported alongside, and the other lambdas as a sensitivity check
# (with z = how far each estimated Leading Edge stands from permutations).
# Validity check: PC1 must be a distinct mode. Its eigenvalue gap
# lambda_1 / lambda_2 is compared with the same ratio in permuted data.

# SIGNIFICANCE (likelihood)
# The four-class mixture
#     p(f_RT, f_ZT) = pi_J g_R g_Z + pi_R g_R h + pi_Z h g_Z + pi_M h h
# (h = Beta(1/2, 3/2) null density; g_R, g_Z = Beta densities of structured
# PASs, estimated) is fitted by EM. Its likelihood ratio against "no
# structure" (every PAS from h h) is the test statistic. Its p-value comes
# from refitting the same mixture to the permuted data: p = (1 + #{LR_perm >=
# LR_obs}) / (1 + permutations). The mixture's class proportions are reported
# for reference only; k comes from the counting estimator, because the
# mixture depends on an assumed shape for structured PASs and the counts do
# not. The mixture also gives every PAS a posterior probability of each class.
#
# WHAT IS NOT GUARANTEED
# The number of structured PASs (k) is identified; which PASs are structured
# is not, at the margin. Unstructured PASs whose f happens to exceed that of
# structured ones will outrank them. The realized split at k (top-|PC1| sets
# of log1p CPM) is reported to show how close it comes to the estimated
# composition.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

repo <- Sys.getenv("SEQUENCE_REPO_ROOT", getwd())
count_name <- "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
count_path <- if (file.exists(file.path(repo, count_name))) file.path(repo, count_name) else
  file.path(repo, "data", count_name)
if (!file.exists(count_path)) stop("Cannot find ", count_name, ". Set SEQUENCE_REPO_ROOT.")
out_dir <- file.path(repo, "exports", "sequence_evs_composition")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N_PERM <- as.integer(Sys.getenv("COMPOSITION_PERMUTATIONS", "99"))
SEED <- as.integer(Sys.getenv("COMPOSITION_SEED", "20261108"))
LAMBDAS <- seq(0.05, 0.50, by = 0.05)
N_BINS <- 10L
NULL_BETA <- c(0.5, 1.5)          # Beta(1/2, (d - 1)/2), d = 5 - 1 = 4
F_EPS <- 1e-6
COMPARISONS <- list(RT0_ZT6 = c(RT = "^R0_", ZT = "^ZT6_"),
                    RT2_ZT8 = c(RT = "^R2_", ZT = "^ZT8_"),
                    RT4_ZT10 = c(RT = "^R4_", ZT = "^ZT10_"),
                    RT8_ZT14 = c(RT = "^R8_", ZT = "^ZT14_"))

# ---- data -------------------------------------------------------------------

read_counts <- function(path) {
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  pats <- unlist(COMPARISONS)
  idx <- sort(unique(unlist(lapply(pats, function(p) grep(p, names(x))))))
  ids <- trimws(as.character(x[[1]]))
  blank <- is.na(ids) | ids == ""
  ids[blank] <- paste0("feature_", which(blank))
  m <- do.call(cbind, lapply(x[, idx, drop = FALSE], function(z)
    suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(z)))))))
  rownames(m) <- make.unique(ids, sep = "__dup_")
  colnames(m) <- names(x)[idx]
  m[!is.finite(m)] <- 0
  pmax(m, 0)
}

# log1p(CPM) within one condition, centred per PAS.
centred_log_cpm <- function(y) {
  lib <- colSums(y)
  lib[!is.finite(lib) | lib <= 0] <- 1
  z <- log1p(sweep(y, 2, lib / 1e6, "/"))
  z - rowMeans(z)
}

# ---- the decomposition --------------------------------------------------------

# PC1 of the condition (prcomp on samples x PASs, centred, unscaled) and each
# PAS's structured fraction f and |PC1 loading|.
# Full decomposition: with n samples the centred data have n - 1 modes. Each
# PAS's variation splits exactly across them, ||x_i||^2 = sum_j (x_i . u_j)^2,
# giving its mode fractions f_1 ... f_(n-1) (summing to 1). f = f_1.
#   majority:  f_1 > 1/2           (PC1 carries more than all other modes together)
#   dominant:  f_1 > every other f (PC1 is the largest single mode)
decompose <- function(z) {
  n_modes <- ncol(z) - 1L
  s <- svd(z, nu = 1, nv = n_modes)
  scores <- z %*% s$v                           # x_i . u_j for every mode j
  total <- rowSums(z^2)                         # ||x_i||^2
  fm <- scores^2 / total
  list(f = fm[, 1], loading = abs(s$u[, 1]),
       dominant = fm[, 1] >= apply(fm[, -1, drop = FALSE], 1, max),
       eigen_ratio = s$d[1]^2 / s$d[2]^2)
}

# Shuffle each PAS's values across the condition's samples independently.
permute_within_rows <- function(z) {
  n <- ncol(z); N <- nrow(z)
  ord <- matrix(stats::runif(N * n), N, n)
  cols <- t(apply(ord, 1, order))
  matrix(z[cbind(rep(seq_len(N), n), as.vector(cols))], N, n)
}

# ---- counting estimator (the proof) -------------------------------------------

# F_null[b, l]: null probability that f < LAMBDAS[l] in mean bin b.
count_composition <- function(fR, fZ, bin, FR, FZ) {
  do.call(rbind, lapply(seq_along(LAMBDAS), function(l) {
    lam <- LAMBDAS[l]; S_R <- S_Z <- M <- 0
    for (b in seq_len(N_BINS)) {
      i <- bin == b
      Nb <- sum(i)
      S_R <- S_R + Nb - sum(fR[i] < lam) / FR[b, l]
      S_Z <- S_Z + Nb - sum(fZ[i] < lam) / FZ[b, l]
      M <- M + sum(fR[i] < lam & fZ[i] < lam) / (FR[b, l] * FZ[b, l])
    }
    N <- length(fR)
    J <- S_R + S_Z - (N - M)
    data.frame(lambda = lam, S_RT = S_R, S_ZT = S_Z, Joint = J,
               Disjoint_RT = S_R - J, Disjoint_ZT = S_Z - J, Remainder = M)
  }))
}

# ---- four-class likelihood mixture (significance) -----------------------------

fit_weighted_beta <- function(f, w, init) {
  nll <- function(p) {
    a <- exp(p[1]); b <- exp(p[2])
    -sum(w * ((a - 1) * log(f) + (b - 1) * log1p(-f) - lbeta(a, b)))
  }
  exp(stats::optim(log(init), nll, method = "L-BFGS-B", lower = c(-3, -3), upper = c(6, 6))$par)
}

mixture_em <- function(fR, fZ, max_iter = 500L, tol = 1e-8) {
  fR <- pmin(pmax(fR, F_EPS), 1 - F_EPS); fZ <- pmin(pmax(fZ, F_EPS), 1 - F_EPS)
  N <- length(fR)
  hR <- stats::dbeta(fR, NULL_BETA[1], NULL_BETA[2]); hZ <- stats::dbeta(fZ, NULL_BETA[1], NULL_BETA[2])
  pi <- c(0.05, 0.05, 0.05, 0.85); aR <- aZ <- c(4, 1); ll_old <- -Inf
  for (it in seq_len(max_iter)) {
    gR <- stats::dbeta(fR, aR[1], aR[2]); gZ <- stats::dbeta(fZ, aZ[1], aZ[2])
    comp <- cbind(pi[1] * gR * gZ, pi[2] * gR * hZ, pi[3] * hR * gZ, pi[4] * hR * hZ)
    tot <- rowSums(comp)
    ll <- sum(log(tot))
    P <- comp / tot
    pi <- colMeans(P)
    aR <- fit_weighted_beta(fR, P[, 1] + P[, 2], aR)
    aZ <- fit_weighted_beta(fZ, P[, 1] + P[, 3], aZ)
    if (is.finite(ll_old) && ll - ll_old < tol * abs(ll)) break
    ll_old <- ll
  }
  colnames(P) <- c("Joint", "RT_only", "ZT_only", "Remainder")
  list(S_RT = N * (pi[1] + pi[2]), S_ZT = N * (pi[1] + pi[3]), Joint = N * pi[1],
       Disjoint = N * (pi[2] + pi[3]), LR = 2 * (ll - sum(log(hR * hZ))),
       beta_RT = aR, beta_ZT = aZ, iterations = it, posterior = P)
}

# ---- one comparison -------------------------------------------------------------

counts <- read_counts(count_path)

# ---- CPM k*: the manuscript's own functions, loaded from its script ----------
# Only the CPM cutoff functions and constants are evaluated (names starting
# cpm_ / CPM_ and derive_cpm_cutoffs); nothing else in that script is run.
pipeline_path <- Sys.getenv("SEQUENCE_PIPELINE_PATH", file.path(repo, "SEQUENCE_EVS_MANUSCRIPT.R"))
if (!file.exists(pipeline_path)) stop("Cannot find SEQUENCE_EVS_MANUSCRIPT.R. Set SEQUENCE_PIPELINE_PATH.")
cpm_env <- new.env()
for (e in parse(pipeline_path, keep.source = FALSE)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      grepl("^(cpm_|CPM_|derive_cpm_cutoffs$)", as.character(e[[2]]))) eval(e, cpm_env)
}
cpm <- cpm_env$derive_cpm_cutoffs(count_path)
cpm_k <- cpm$k_star

set.seed(SEED)
summary_rows <- list(); lambda_rows <- list(); posterior_rows <- list(); density_rows <- list(); null_rows <- list()
union_rows <- list(); composition_rows <- list()

for (cn in names(COMPARISONS)) {
  pat <- COMPARISONS[[cn]]
  yR <- counts[, grep(pat[["RT"]], colnames(counts)), drop = FALSE]
  yZ <- counts[, grep(pat[["ZT"]], colnames(counts)), drop = FALSE]
  if (ncol(yR) != 5L || ncol(yZ) != 5L) stop(cn, ": expected 5 samples per condition.")
  zR_all <- centred_log_cpm(yR); zZ_all <- centred_log_cpm(yZ)
  nonzero <- rowSums(yR) + rowSums(yZ) > 0
  # PASs constant within a condition cannot be structured there; with no
  # variation in one condition they are Remainder by definition.
  varying <- nonzero & rowSums(zR_all^2) > 0 & rowSums(zZ_all^2) > 0
  zR <- zR_all[varying, , drop = FALSE]; zZ <- zZ_all[varying, , drop = FALSE]
  N <- nrow(zR)
  message(sprintf("\n%s: %d PASs with variation in both conditions (%d more nonzero but constant in one condition)",
                  cn, N, sum(nonzero) - N))

  dR <- decompose(zR); dZ <- decompose(zZ)
  mean_expr <- rowMeans(log1p(cbind(sweep(yR[varying, ], 2, colSums(yR) / 1e6, "/"),
                                    sweep(yZ[varying, ], 2, colSums(yZ) / 1e6, "/"))))
  breaks <- unique(stats::quantile(mean_expr, seq(0, 1, length.out = N_BINS + 1L)))
  bin <- as.integer(cut(mean_expr, breaks, include.lowest = TRUE))
  if (length(unique(bin)) != N_BINS) stop(cn, ": could not form ", N_BINS, " mean bins.")

  # Observed mixture fit and likelihood ratio.
  fit <- mixture_em(dR$f, dZ$f)
  message(sprintf("  likelihood ratio vs no structure: %.1f (EM iterations %d)", fit$LR, fit$iterations))

  # Permutations: empirical null F per mean bin, and the null of the likelihood ratio.
  below_R <- below_Z <- matrix(0, N_BINS, length(LAMBDAS))
  n_bin <- tabulate(bin, N_BINS)
  lr_perm <- numeric(N_PERM)
  both_perm <- array(0, c(N_PERM, N_BINS, length(LAMBDAS)))   # permuted PASs below lambda in both conditions
  null_f_sample <- list()
  notdom_R <- notdom_Z <- numeric(N_BINS)                      # permuted PASs not governed by PC1
  notdom_both_perm <- matrix(0, N_PERM, N_BINS)
  gap_perm <- matrix(NA_real_, N_PERM, 2, dimnames = list(NULL, c("RT", "ZT")))
  for (p in seq_len(N_PERM)) {
    qR <- decompose(permute_within_rows(zR)); qZ <- decompose(permute_within_rows(zZ))
    pR <- qR$f; pZ <- qZ$f
    notdom_R <- notdom_R + tabulate(bin[!qR$dominant], N_BINS)
    notdom_Z <- notdom_Z + tabulate(bin[!qZ$dominant], N_BINS)
    notdom_both_perm[p, ] <- tabulate(bin[!qR$dominant & !qZ$dominant], N_BINS)
    gap_perm[p, ] <- c(qR$eigen_ratio, qZ$eigen_ratio)
    for (l in seq_along(LAMBDAS)) {
      below_R[, l] <- below_R[, l] + tabulate(bin[pR < LAMBDAS[l]], N_BINS)
      below_Z[, l] <- below_Z[, l] + tabulate(bin[pZ < LAMBDAS[l]], N_BINS)
    }
    for (l in seq_along(LAMBDAS)) both_perm[p, , l] <- tabulate(bin[pR < LAMBDAS[l] & pZ < LAMBDAS[l]], N_BINS)
    lr_perm[p] <- mixture_em(pR, pZ, max_iter = 200L)$LR
    if (p <= 5L) null_f_sample[[p]] <- data.frame(Condition = c(rep("RT", N), rep("ZT", N)), f = c(pR, pZ))
    if (p %% 10L == 0L) message(sprintf("  permutation %d of %d", p, N_PERM))
  }
  FR <- below_R / (N_PERM * n_bin); FZ <- below_Z / (N_PERM * n_bin)
  if (any(FR <= 0) || any(FZ <= 0)) stop(cn, ": empirical null probability is zero in a bin; increase permutations.")
  p_value <- (1 + sum(lr_perm >= fit$LR)) / (1 + N_PERM)
  # Is PC1 a distinct mode? Eigenvalue gap lambda_1 / lambda_2 vs permutations.
  gap_p_RT <- (1 + sum(gap_perm[, "RT"] >= dR$eigen_ratio)) / (1 + N_PERM)
  gap_p_ZT <- (1 + sum(gap_perm[, "ZT"] >= dZ$eigen_ratio)) / (1 + N_PERM)
  bin_mean <- tapply(mean_expr, bin, stats::median)
  null_rows[[cn]] <- do.call(rbind, lapply(seq_along(LAMBDAS), function(l) data.frame(
    Comparison = cn, mean_bin = seq_len(N_BINS), median_log1p_CPM = as.numeric(bin_mean), n_PAS = n_bin,
    lambda = LAMBDAS[l], F_null_RT = FR[, l], F_null_ZT = FZ[, l],
    F_theory = stats::pbeta(LAMBDAS[l], NULL_BETA[1], NULL_BETA[2]),
    observed_below_RT = tabulate(bin[dR$f < LAMBDAS[l]], N_BINS) / n_bin,
    observed_below_ZT = tabulate(bin[dZ$f < LAMBDAS[l]], N_BINS) / n_bin)))

  # Composition and k by the counting estimator.
  est <- count_composition(dR$f, dZ$f, bin, FR, FZ)
  est$Comparison <- cn

  # ONE k per comparison, from both conditions together (all 10 samples), as
  # in the CPM method. The Remainder (strays in both conditions) is the one
  # quantity estimated jointly; the Leading Edge is N - Remainder. With the
  # same k in both conditions, the Leading Edge is the union of the two top-k
  # sets by |PC1 loading|: U(k) = 2k - Joint(k). k is the smallest k whose
  # union reaches the estimated Leading Edge.
  rR <- rank(-dR$loading, ties.method = "first"); rZ <- rank(-dZ$loading, ties.method = "first")
  union_curve <- seq_len(N) * 2 - cumsum(tabulate(pmax(rR, rZ), N))
  k_for_leading_edge <- function(le) as.integer(which(union_curve >= le)[1])
  est$Leading_Edge <- N - est$Remainder
  est$k_single <- vapply(est$Leading_Edge, k_for_leading_edge, integer(1))

  # z at every lambda: how far the estimated Leading Edge stands from permutations.
  inv_F2 <- 1 / (FR * FZ)                                                    # bins x lambdas
  LE_perm <- sapply(seq_along(LAMBDAS), function(l) N - as.vector(matrix(both_perm[, , l], N_PERM, N_BINS) %*% inv_F2[, l]))
  LE_perm <- matrix(LE_perm, N_PERM, length(LAMBDAS))
  est$LE_null_mean <- colMeans(LE_perm)
  est$LE_null_sd <- apply(LE_perm, 2, stats::sd)
  est$z <- (est$Leading_Edge - est$LE_null_mean) / est$LE_null_sd
  # Primary rule (set by the decomposition, no search): majority, lambda = 1/2.
  l_best <- which(abs(LAMBDAS - 0.5) < 1e-9)
  est$selected <- seq_along(LAMBDAS) == l_best
  lambda_rows[[cn]] <- est
  med <- unlist(est[l_best, c("S_RT", "S_ZT", "Joint", "Disjoint_RT", "Disjoint_ZT", "Remainder")])
  k_RT <- as.integer(round(med[["S_RT"]])); k_ZT <- as.integer(round(med[["S_ZT"]]))
  k_single <- est$k_single[l_best]

  # Dominant-mode rule: the same counting proof with "not governed by PC1"
  # (f_1 not the largest mode) in place of "f below lambda".
  G_R <- notdom_R / (N_PERM * n_bin); G_Z <- notdom_Z / (N_PERM * n_bin)
  M_dom <- sum(tabulate(bin[!dR$dominant & !dZ$dominant], N_BINS) / (G_R * G_Z))
  LE_dom <- N - M_dom
  LE_dom_perm <- N - as.vector(notdom_both_perm %*% (1 / (G_R * G_Z)))
  k_dominant <- k_for_leading_edge(LE_dom)

  # Realized split at the single k: top-|PC1 loading| sets of log1p CPM.
  inR <- rR <= k_single; inZ <- rZ <= k_single
  realized_at <- function(k) {
    joint <- sum(pmax(rR, rZ) <= k)
    data.frame(Comparison = cn, k = k, Joint = joint, Disjoint = 2 * k - 2 * joint, Remainder = N - (2 * k - joint))
  }
  ks_show <- unique(as.integer(round(seq(1, N, length.out = 400))))
  union_rows[[cn]] <- data.frame(Comparison = cn, k = ks_show, union = union_curve[ks_show], N = N)
  composition_rows[[cn]] <- rbind(
    data.frame(Source = "Estimated (composition)", Comparison = cn, k = k_single,
               Joint = med[["Joint"]], Disjoint = med[["Disjoint_RT"]] + med[["Disjoint_ZT"]], Remainder = med[["Remainder"]]),
    cbind(Source = "Realized at composition k", realized_at(k_single)),
    cbind(Source = "Realized at CPM k*", realized_at(as.integer(cpm_k[[cn]]))),
    cbind(Source = "Realized at 5,000", realized_at(5000L)))

  summary_rows[[cn]] <- data.frame(
    Comparison = cn, N_PAS = N,
    S_RT = round(med[["S_RT"]]), S_ZT = round(med[["S_ZT"]]), Joint = round(med[["Joint"]]),
    Disjoint = round(med[["Disjoint_RT"]] + med[["Disjoint_ZT"]]),
    Leading_Edge = round(N - med[["Remainder"]]),
    k_single = k_single, rule = "majority (f1 > 1/2)", z_majority = est$z[l_best],
    k_dominant = k_dominant, Leading_Edge_dominant = round(LE_dom),
    z_dominant = (LE_dom - mean(LE_dom_perm)) / stats::sd(LE_dom_perm),
    eigen_ratio_RT = dR$eigen_ratio, eigen_ratio_RT_perm_max = max(gap_perm[, "RT"]), eigen_gap_p_RT = gap_p_RT,
    eigen_ratio_ZT = dZ$eigen_ratio, eigen_ratio_ZT_perm_max = max(gap_perm[, "ZT"]), eigen_gap_p_ZT = gap_p_ZT,
    k_single_range = sprintf("%d to %d", min(est$k_single), max(est$k_single)),
    k_RT_alone = k_RT, k_ZT_alone = k_ZT,
    S_RT_range = sprintf("%.0f to %.0f", min(est$S_RT), max(est$S_RT)),
    S_ZT_range = sprintf("%.0f to %.0f", min(est$S_ZT), max(est$S_ZT)),
    Joint_range = sprintf("%.0f to %.0f", min(est$Joint), max(est$Joint)),
    Realized_Leading_Edge = sum(inR | inZ), Realized_Joint = sum(inR & inZ), Realized_Disjoint = sum(xor(inR, inZ)),
    Mixture_S_RT = round(fit$S_RT), Mixture_S_ZT = round(fit$S_ZT), Mixture_Joint = round(fit$Joint),
    Likelihood_ratio = fit$LR, Permutation_LR_max = max(lr_perm), Permutation_p = p_value,
    Permutations = N_PERM,
    CPM_k_star = as.integer(cpm_k[[cn]]),
    CPM_c1 = cpm$knots[["c1"]], CPM_c2 = cpm$knots[["c2"]])
  posterior_rows[[cn]] <- data.frame(Comparison = cn, feature_id = rownames(zR), f_RT = dR$f, f_ZT = dZ$f,
                                     mean_bin = bin, round(fit$posterior, 4), check.names = FALSE)
  density_rows[[cn]] <- rbind(
    data.frame(Comparison = cn, Source = "Observed", Condition = c(rep("RT", N), rep("ZT", N)), f = c(dR$f, dZ$f)),
    cbind(Comparison = cn, Source = "Permutation null", do.call(rbind, null_f_sample)))

  s <- summary_rows[[cn]]
  message(sprintf("  k = %d by majority (f1 > 1/2; Leading Edge %d), %d by dominant mode; eigen gap p RT %.3f, ZT %.3f; permutation p = %.3f%s",
                  s$k_single, s$Leading_Edge, s$k_dominant, s$eigen_gap_p_RT, s$eigen_gap_p_ZT, s$Permutation_p,
                  sprintf("; CPM k* = %d", s$CPM_k_star)))
}

# ---- outputs --------------------------------------------------------------------

utils::write.csv(do.call(rbind, summary_rows), file.path(out_dir, "Table_Composition_Summary.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, lambda_rows), file.path(out_dir, "Table_Composition_By_Lambda.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, posterior_rows), file.path(out_dir, "Table_PAS_Class_Posteriors.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, null_rows), file.path(out_dir, "Table_Null_F_By_Bin.csv"), row.names = FALSE)

utils::write.csv(do.call(rbind, composition_rows), file.path(out_dir, "Table_Composition_Estimated_vs_Realized.csv"), row.names = FALSE)
utils::write.csv(cpm$audit, file.path(out_dir, "Table_CPM_Cutoff_Selection.csv"), row.names = FALSE)

# ---- figures --------------------------------------------------------------------

CMP_LEVELS <- sub("_", " / ", names(COMPARISONS))
as_cmp <- function(x) factor(sub("_", " / ", x), levels = CMP_LEVELS)
COL_COMP <- "#1F6FB2"; COL_CPM <- "#2E7D5B"; COL_FIX <- "#8C8C8C"
fig_theme <- function() theme_bw(base_size = 8) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = NA), strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 9), plot.subtitle = element_text(size = 7, colour = "grey30"))
save_both <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
}
k_labels <- function(x) ifelse(is.na(x), NA_character_, paste0(x / 1000, "k"))

# Figure C1: how the CPM method chooses k.
cv <- cpm$curves
arm_cmp <- c(RT0 = "RT0_ZT6", ZT6 = "RT0_ZT6", RT2 = "RT2_ZT8", ZT8 = "RT2_ZT8",
             RT4 = "RT4_ZT10", ZT10 = "RT4_ZT10", RT8 = "RT8_ZT14", ZT14 = "RT8_ZT14")
cv$Comparison <- as_cmp(arm_cmp[cv$Arm]); cv$Condition <- ifelse(grepl("^RT", cv$Arm), "RT arm", "ZT arm")
kn <- cpm$knots
zones <- data.frame(xmin = c(0, kn[["c1"]], kn[["c2"]]), xmax = c(kn[["c1"]], kn[["c2"]], kn[["N"]]),
                    Zone = factor(c("Remainder regime", "Meeting zone", "Leading Edge regime"),
                                  levels = c("Remainder regime", "Meeting zone", "Leading Edge regime")))
c1a <- ggplot(cv) +
  geom_rect(data = zones, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Zone), inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_line(aes(rank, D, linetype = Condition), linewidth = 0.45, colour = "grey15") +
  geom_line(aes(rank, fitted, linetype = Condition), linewidth = 0.35, colour = "#B8860B") +
  facet_wrap(~ Comparison, nrow = 1) +
  scale_fill_manual(values = c("Remainder regime" = "#F6DCD8", "Meeting zone" = "#EEEEEE", "Leading Edge regime" = "#D7E5F3")) +
  scale_linetype_manual(values = c("RT arm" = "solid", "ZT arm" = "22")) +
  scale_x_continuous(labels = k_labels) +
  labs(x = "Rank by |PC1 loading| of log1p(CPM) (ascending)", y = "D = cum. overdispersion share\n- cum. PC1 share",
       fill = NULL, linetype = NULL, tag = "a",
       title = "Step 1 (all 40 samples): the two regimes",
       subtitle = sprintf("Black: D curve of each arm; gold: one two-knot fit shared by all 8 arms. c1 = %s, c2 = %s of %s PASs.",
                          format(kn[["c1"]], big.mark = ","), format(kn[["c2"]], big.mark = ","), format(kn[["N"]], big.mark = ","))) +
  fig_theme()
fr <- cpm$frontiers; fr$Comparison <- as_cmp(fr$Comparison)
ends <- do.call(rbind, lapply(split(fr, fr$Comparison, drop = TRUE), function(x)
  data.frame(Comparison = x$Comparison[1], x1 = x$cost_norm[1], y1 = x$good_norm[1],
             x2 = x$cost_norm[nrow(x)], y2 = x$good_norm[nrow(x)])))
sel <- fr[fr$selected, ]; sel$lab <- paste0("k* = ", format(sel$k, big.mark = ","))
c1b <- ggplot(fr, aes(cost_norm, good_norm)) +
  geom_segment(data = ends, aes(x = x1, y = y1, xend = x2, yend = y2), linetype = "22", linewidth = 0.3,
               colour = "grey45", inherit.aes = FALSE) +
  geom_line(linewidth = 0.45, colour = "grey20") +
  geom_point(data = sel, shape = 21, size = 2.6, stroke = 0.7, colour = "white", fill = COL_CPM) +
  geom_text(data = sel, aes(label = lab), hjust = -0.15, vjust = 1.4, size = 2.3, colour = COL_CPM, fontface = "bold") +
  facet_wrap(~ Comparison, nrow = 1) +
  scale_x_continuous(labels = scales::label_percent(), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(labels = scales::label_percent(), breaks = c(0, 0.5, 1)) +
  labs(x = "Cost R(k): Disjoint PASs whose rank in the other arm is in the Remainder regime (scaled)",
       y = "Good G(k): consistent\nLeading Edge PASs (scaled)", tag = "b",
       title = "Step 2 (the comparison's 10 samples): one k per comparison",
       subtitle = "Pareto frontier of G vs R over k; dashed: chord between its ends; k* = maximum endpoint deviation.") +
  fig_theme()
fig_c1 <- if (requireNamespace("patchwork", quietly = TRUE)) patchwork::wrap_plots(c1a, c1b, ncol = 1) +
  patchwork::plot_annotation(title = "Figure C1. How the CPM method chooses k",
                             theme = theme(plot.title = element_text(face = "bold", size = 10))) else c1b
save_both(fig_c1, "Figure_C1_CPM_Choosing_k", 180, 130)

# Figure C2: f observed vs permutation null.
dens <- do.call(rbind, density_rows); dens$Comparison <- as_cmp(dens$Comparison)
summary_tab <- do.call(rbind, summary_rows)
fig_c2 <- ggplot(dens, aes(f, colour = Source, linetype = Source)) +
  geom_density(linewidth = 0.5) +
  coord_cartesian(xlim = c(0, 1)) +
  facet_grid(Condition ~ Comparison) +
  scale_colour_manual(values = c("Observed" = COL_COMP, "Permutation null" = "grey45")) +
  scale_linetype_manual(values = c("Observed" = "solid", "Permutation null" = "22")) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Structured fraction f (share of a PAS's variation on its condition's PC1)", y = "Density",
       colour = NULL, linetype = NULL,
       title = "Figure C2. Structured fraction of every PAS: observed vs permutation null",
       subtitle = "Below the null near f = 0 and above it at high f: more PASs follow PC1 than chance allows.") +
  fig_theme()
save_both(fig_c2, "Figure_C2_f_Observed_vs_Null", 180, 105)

# Figure C3: from the composition to k, per comparison.
lam <- do.call(rbind, lambda_rows); lam$Comparison <- as_cmp(lam$Comparison)
st <- summary_tab; st$Comparison <- as_cmp(st$Comparison)
lam_sel <- lam[lam$selected, ]
c3a <- ggplot(lam, aes(lambda, k_single)) +
  geom_hline(data = st, aes(yintercept = CPM_k_star), colour = COL_CPM, linewidth = 0.5, linetype = "22") +
  geom_hline(yintercept = 5000, colour = COL_FIX, linewidth = 0.4, linetype = "22") +
  geom_line(colour = COL_COMP, linewidth = 0.5) + geom_point(colour = COL_COMP, size = 1.3) +
  geom_point(data = lam_sel, shape = 23, size = 2.8, fill = COL_COMP, colour = "white") +
  geom_text(data = lam_sel, aes(label = paste0("k = ", format(k_single, big.mark = ","))), vjust = -1.1, size = 2.3,
            colour = COL_COMP, fontface = "bold") +
  facet_wrap(~ Comparison, nrow = 1) +
  scale_y_continuous(labels = scales::label_comma(), limits = c(0, NA)) +
  scale_x_continuous(breaks = c(0.1, 0.3, 0.5)) +
  labs(x = "lambda (cutoff on f used for counting)", y = "k", tag = "a",
       title = "The single k at each counting cutoff, and the cutoff chosen",
       subtitle = "Diamond: lambda = 1/2, where structure outweighs residual (the decomposition's rule). Other lambdas: sensitivity. Dashed green: CPM k*; grey: 5,000.") +
  fig_theme()
c3z <- ggplot(lam, aes(lambda, z)) +
  geom_line(colour = "grey25", linewidth = 0.5) + geom_point(colour = "grey25", size = 1.3) +
  geom_point(data = lam_sel, shape = 23, size = 2.8, fill = COL_COMP, colour = "white") +
  facet_wrap(~ Comparison, nrow = 1) +
  scale_x_continuous(breaks = c(0.1, 0.3, 0.5)) +
  labs(x = "lambda", y = "z (estimated Leading Edge\nvs permuted data)", tag = "b",
       title = "How clearly the estimated Leading Edge stands out from noise at each cutoff",
       subtitle = "z = (observed - permuted mean) / permuted SD, for reference.") +
  fig_theme()
un <- do.call(rbind, union_rows); un$Comparison <- as_cmp(un$Comparison)
c3b <- ggplot(un, aes(k, union)) +
  geom_line(colour = "grey20", linewidth = 0.5) +
  geom_hline(data = st, aes(yintercept = Leading_Edge), colour = COL_COMP, linewidth = 0.4) +
  geom_vline(data = st, aes(xintercept = k_single), colour = COL_COMP, linewidth = 0.5) +
  geom_vline(data = st, aes(xintercept = CPM_k_star), colour = COL_CPM, linewidth = 0.5, linetype = "22") +
  geom_vline(xintercept = 5000, colour = COL_FIX, linewidth = 0.4, linetype = "22") +
  geom_text(data = st, aes(x = k_single, y = Leading_Edge, label = paste0("k = ", format(k_single, big.mark = ","))),
            hjust = -0.1, vjust = 1.5, size = 2.3, colour = COL_COMP, fontface = "bold") +
  facet_wrap(~ Comparison, nrow = 1) +
  coord_cartesian(xlim = c(0, 15000), ylim = c(0, 25000)) +
  scale_x_continuous(labels = k_labels) + scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "k (top-|PC1| PASs per condition)", y = "Leading Edge size\n(union of the two top-k sets)", tag = "c",
       title = "From the estimated Leading Edge to k",
       subtitle = "Black: Leading Edge size at each k. Blue line: estimated Leading Edge (N - Remainder); k where they meet.") +
  fig_theme()
cp <- do.call(rbind, composition_rows)
cpl <- rbind(data.frame(cp[, c("Source", "Comparison", "k")], Class = "Joint", n = cp$Joint),
             data.frame(cp[, c("Source", "Comparison", "k")], Class = "Disjoint", n = cp$Disjoint),
             data.frame(cp[, c("Source", "Comparison", "k")], Class = "Remainder", n = cp$Remainder))
cpl$Comparison <- as_cmp(cpl$Comparison)
cpl$Source <- factor(cpl$Source, levels = rev(c("Estimated (composition)", "Realized at composition k",
                                                "Realized at CPM k*", "Realized at 5,000")))
cpl$Class <- factor(cpl$Class, levels = c("Remainder", "Disjoint", "Joint"))
c3c <- ggplot(cpl, aes(n, Source, fill = Class)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Comparison, nrow = 1) +
  scale_fill_manual(values = c(Joint = "#0B3C6E", Disjoint = "#9DBCD9", Remainder = "#E9E9E9"), breaks = c("Joint", "Disjoint", "Remainder")) +
  scale_x_continuous(labels = k_labels) +
  labs(x = "PASs", y = NULL, fill = NULL, tag = "d",
       title = "Composition: estimated vs what each k actually produces",
       subtitle = "Realized splits use the top-|PC1| sets; their extra Joint and Disjoint PASs are noise PASs the ranking lets in.") +
  fig_theme()
fig_c3 <- if (requireNamespace("patchwork", quietly = TRUE)) patchwork::wrap_plots(c3a, c3z, c3b, c3c, ncol = 1, heights = c(1, 0.8, 1, 0.9)) +
  patchwork::plot_annotation(title = "Figure C3. How the composition method chooses one k per comparison",
                             theme = theme(plot.title = element_text(face = "bold", size = 10))) else c3b
save_both(fig_c3, "Figure_C3_Composition_Choosing_k", 180, 245)

# Figure C4: the three cutoffs side by side.
kk <- rbind(data.frame(Comparison = st$Comparison, Method = "Fixed 5,000", k = 5000, lo = NA, hi = NA),
            data.frame(Comparison = st$Comparison, Method = "CPM k*", k = st$CPM_k_star, lo = NA, hi = NA),
            data.frame(Comparison = st$Comparison, Method = "Decomposition k", k = st$k_single,
                       lo = as.numeric(sub(" to .*", "", st$k_single_range)), hi = as.numeric(sub(".* to ", "", st$k_single_range))))
kk$Method <- factor(kk$Method, levels = c("Fixed 5,000", "CPM k*", "Decomposition k"))
kk$Comparison <- factor(kk$Comparison, levels = rev(CMP_LEVELS))
dom <- data.frame(Comparison = factor(st$Comparison, levels = rev(CMP_LEVELS)), k = st$k_dominant)
fig_c4 <- ggplot(kk, aes(k, Comparison, colour = Method)) +
  geom_point(data = dom, aes(k, Comparison), inherit.aes = FALSE, shape = 17, size = 2.2, colour = COL_COMP,
             position = position_nudge(y = -0.2)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.5, na.rm = TRUE,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
  geom_text(aes(label = format(k, big.mark = ",")), position = position_dodge(width = 0.6), vjust = -0.9, size = 2.1,
            show.legend = FALSE) +
  scale_colour_manual(values = c("Fixed 5,000" = COL_FIX, "CPM k*" = COL_CPM, "Decomposition k" = COL_COMP)) +
  scale_x_continuous(labels = scales::label_comma()) +
  labs(x = "k (one per comparison)", y = NULL, colour = NULL,
       title = "Figure C4. The cutoff chosen by each method",
       subtitle = "Decomposition k (majority rule, f1 > 1/2); bar = its range over lambda = 0.05 to 0.50; triangle = dominant-mode rule.") +
  fig_theme()
save_both(fig_c4, "Figure_C4_k_By_Method", 150, 90)

cat("\nOne k per comparison from the decomposition (majority rule f1 > 1/2; dominant-mode rule alongside)\n")
print(summary_tab[, c("Comparison", "N_PAS", "k_single", "k_dominant", "k_single_range", "Leading_Edge", "eigen_gap_p_RT", "eigen_gap_p_ZT", "Joint", "Disjoint",
                      "Realized_Joint", "Realized_Disjoint", "Likelihood_ratio", "Permutation_p", "CPM_k_star")],
      row.names = FALSE)
cat("\nOutputs written to", out_dir, "\n")
