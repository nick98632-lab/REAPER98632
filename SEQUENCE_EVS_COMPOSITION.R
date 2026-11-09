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
#   Standalone: it does not change, or depend on, the manuscript runs. If the
#   manuscript run exists, its CPM k* values are read for comparison.
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
#   Each estimator is unbiased. The EVS split contains exactly the structured
#   PASs when each condition's top-k set holds its structured PASs, so
#       k_RT = S_RT = J + D_R,   k_ZT = S_ZT = J + D_Z.
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
# estimators above are applied within each bin and summed. Estimates are
# reported for lambda = 0.1 to 0.5; the summary uses their median.
#
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
LAMBDAS <- c(0.1, 0.2, 0.3, 0.4, 0.5)
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
decompose <- function(z) {
  s <- svd(z, nu = 1, nv = 1)
  score <- as.vector(z %*% s$v[, 1])            # x_i . u_c
  total <- rowSums(z^2)                         # ||x_i||^2
  list(f = score^2 / total, loading = abs(s$u[, 1]))
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
cpm_manifest <- file.path(repo, "exports", "sequence_evs_manuscript", "Cutoff_Method_Manifest.csv")
cpm_k <- if (file.exists(cpm_manifest)) {
  mf <- utils::read.csv(cpm_manifest, check.names = FALSE)
  unlist(mf[mf$Cutoff_Method == "CPM_k", names(COMPARISONS)])
} else NULL

set.seed(SEED)
summary_rows <- list(); lambda_rows <- list(); posterior_rows <- list(); density_rows <- list()

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
  null_f_sample <- list()
  for (p in seq_len(N_PERM)) {
    pR <- decompose(permute_within_rows(zR))$f
    pZ <- decompose(permute_within_rows(zZ))$f
    for (l in seq_along(LAMBDAS)) {
      below_R[, l] <- below_R[, l] + tabulate(bin[pR < LAMBDAS[l]], N_BINS)
      below_Z[, l] <- below_Z[, l] + tabulate(bin[pZ < LAMBDAS[l]], N_BINS)
    }
    lr_perm[p] <- mixture_em(pR, pZ, max_iter = 200L)$LR
    if (p <= 5L) null_f_sample[[p]] <- data.frame(Condition = c(rep("RT", N), rep("ZT", N)), f = c(pR, pZ))
    if (p %% 10L == 0L) message(sprintf("  permutation %d of %d", p, N_PERM))
  }
  FR <- below_R / (N_PERM * n_bin); FZ <- below_Z / (N_PERM * n_bin)
  if (any(FR <= 0) || any(FZ <= 0)) stop(cn, ": empirical null probability is zero in a bin; increase permutations.")
  p_value <- (1 + sum(lr_perm >= fit$LR)) / (1 + N_PERM)

  # Composition and k by the counting estimator.
  est <- count_composition(dR$f, dZ$f, bin, FR, FZ)
  est$Comparison <- cn
  lambda_rows[[cn]] <- est
  med <- vapply(est[, c("S_RT", "S_ZT", "Joint", "Disjoint_RT", "Disjoint_ZT", "Remainder")], stats::median, numeric(1))
  k_RT <- as.integer(round(med[["S_RT"]])); k_ZT <- as.integer(round(med[["S_ZT"]]))

  # Realized split at those k: top-|PC1 loading| sets of log1p CPM.
  rR <- rank(-dR$loading, ties.method = "first"); rZ <- rank(-dZ$loading, ties.method = "first")
  inR <- rR <= k_RT; inZ <- rZ <= k_ZT

  summary_rows[[cn]] <- data.frame(
    Comparison = cn, N_PAS = N,
    S_RT = round(med[["S_RT"]]), S_ZT = round(med[["S_ZT"]]), Joint = round(med[["Joint"]]),
    Disjoint = round(med[["Disjoint_RT"]] + med[["Disjoint_ZT"]]),
    Leading_Edge = round(med[["S_RT"]] + med[["S_ZT"]] - med[["Joint"]]),
    k_RT = k_RT, k_ZT = k_ZT, k_common = as.integer(round((k_RT + k_ZT) / 2)),
    S_RT_range = sprintf("%.0f to %.0f", min(est$S_RT), max(est$S_RT)),
    S_ZT_range = sprintf("%.0f to %.0f", min(est$S_ZT), max(est$S_ZT)),
    Joint_range = sprintf("%.0f to %.0f", min(est$Joint), max(est$Joint)),
    Realized_Leading_Edge = sum(inR | inZ), Realized_Joint = sum(inR & inZ), Realized_Disjoint = sum(xor(inR, inZ)),
    Mixture_S_RT = round(fit$S_RT), Mixture_S_ZT = round(fit$S_ZT), Mixture_Joint = round(fit$Joint),
    Likelihood_ratio = fit$LR, Permutation_LR_max = max(lr_perm), Permutation_p = p_value,
    Permutations = N_PERM,
    CPM_k_star = if (is.null(cpm_k)) NA_integer_ else as.integer(cpm_k[[cn]]))
  posterior_rows[[cn]] <- data.frame(Comparison = cn, feature_id = rownames(zR), f_RT = dR$f, f_ZT = dZ$f,
                                     mean_bin = bin, round(fit$posterior, 4), check.names = FALSE)
  density_rows[[cn]] <- rbind(
    data.frame(Comparison = cn, Source = "Observed", Condition = c(rep("RT", N), rep("ZT", N)), f = c(dR$f, dZ$f)),
    cbind(Comparison = cn, Source = "Permutation null", do.call(rbind, null_f_sample)))

  s <- summary_rows[[cn]]
  message(sprintf("  k_RT = %d, k_ZT = %d (Joint %d, Disjoint %d, Leading Edge %d); permutation p = %.3f%s",
                  s$k_RT, s$k_ZT, s$Joint, s$Disjoint, s$Leading_Edge, s$Permutation_p,
                  if (is.na(s$CPM_k_star)) "" else sprintf("; CPM k* = %d", s$CPM_k_star)))
}

# ---- outputs --------------------------------------------------------------------

summary_tab <- do.call(rbind, summary_rows)
utils::write.csv(summary_tab, file.path(out_dir, "Table_Composition_Summary.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, lambda_rows), file.path(out_dir, "Table_Composition_By_Lambda.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, posterior_rows), file.path(out_dir, "Table_PAS_Class_Posteriors.csv"), row.names = FALSE)

dens <- do.call(rbind, density_rows)
dens$Comparison <- factor(sub("_", " / ", dens$Comparison), levels = sub("_", " / ", names(COMPARISONS)))
lab <- summary_tab
lab$Comparison <- factor(sub("_", " / ", lab$Comparison), levels = levels(dens$Comparison))
lab <- rbind(data.frame(Comparison = lab$Comparison, Condition = "RT",
                        text = sprintf("k = %s", format(lab$k_RT, big.mark = ","))),
             data.frame(Comparison = lab$Comparison, Condition = "ZT",
                        text = sprintf("k = %s", format(lab$k_ZT, big.mark = ","))))
fig <- ggplot(dens, aes(f, colour = Source, linetype = Source)) +
  geom_density(linewidth = 0.5) +
  coord_cartesian(xlim = c(0, 1)) +
  geom_text(data = lab, aes(x = 0.97, y = Inf, label = text), inherit.aes = FALSE,
            hjust = 1, vjust = 1.6, size = 2.4, fontface = "bold") +
  facet_grid(Condition ~ Comparison) +
  scale_colour_manual(values = c("Observed" = "#1F6FB2", "Permutation null" = "grey45")) +
  scale_linetype_manual(values = c("Observed" = "solid", "Permutation null" = "22")) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Structured fraction f = cos^2(angle to PC1)", y = "Density", colour = NULL, linetype = NULL,
       title = "Structured fraction of every PAS: observed vs permutation null",
       subtitle = "The excess of observed over null at high f is the structured PASs; k = their estimated number per condition.") +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = NA), strip.text = element_text(face = "bold"))
ggsave(file.path(out_dir, "Figure_Composition_f_Observed_vs_Null.png"), fig, width = 180, height = 105, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "Figure_Composition_f_Observed_vs_Null.pdf"), fig, width = 180, height = 105, units = "mm")

cat("\nComposition estimate of k (median over lambda = 0.1 to 0.5)\n")
print(summary_tab[, c("Comparison", "N_PAS", "k_RT", "k_ZT", "k_common", "Joint", "Disjoint", "Leading_Edge",
                      "Realized_Joint", "Realized_Disjoint", "Likelihood_ratio", "Permutation_p", "CPM_k_star")],
      row.names = FALSE)
cat("\nOutputs written to", out_dir, "\n")
