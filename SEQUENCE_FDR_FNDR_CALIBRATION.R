#!/usr/bin/env Rscript
# =============================================================================
# SEQUENCE FDR / FNDR CALIBRATION  (Supplementary Figure S2)
# =============================================================================
# Simulated counts with known truth, calibrated to each comparison's real data.
# On every simulated dataset the full procedure runs end to end, using the
# manuscript pipeline's own functions (loaded from SEQUENCE_FINAL_MANUSCRIPT.R),
# so every decision rule is identical to the analysis:
#   Original (No EVS)    one DESeq2 fit on all PASs
#   EVS, Fixed 5,000     Leading Edge / Remainder at the fixed cutoff
#   EVS, Likelihood k*   Leading Edge / Remainder at k* re-chosen on the simulated data
# Standard and HBFSS calls come from run_core_analysis() (DESeq2 with the
# log-scale power trend a0 + a1 x mean^(-b), apeglm, fdrtool empirical null,
# higher criticism), exactly as in the manuscript. The truth's dispersion trends are
# fitted the same way.
#
# Truth scenarios (means, size factors and dispersion model fitted on the real
# comparison; 5 RT vs 5 ZT samples):
#   one_population   one dispersion trend and prior width for all PASs
#                    (DESeq2's one-population assumption holds)
#   two_populations  the Leading Edge / Remainder trends and prior widths
#                    fitted on the real data at its k*
#   null             one population and no true differences
#   structured       two populations, plus the real data's sample structure:
#                    each condition's PC1 sample pattern from the real counts is
#                    carried by every PAS in proportion to its real loading
#                    (log scale), so PC1 in the simulated data reflects real
#                    structure, as in the analysis. Each PAS's total dispersion
#                    is kept at its truth value: the structured part is taken
#                    out of the NB dispersion, alpha_NB = (alpha - v) / (1 + v),
#                    with v the structured squared CV (floored at 5% of alpha).
#   structured_null  structured, with no true differences
# In the non-null scenarios a fraction SEQUENCE_CAL_DE_FRACTION of PASs gets a
# true log2 fold change of random sign and size 0.5 + Gamma(2, 0.5).
#
# Metrics, per simulated dataset, method and call (Standard, HBFSS):
#   FDR    false calls / calls                       (lower is better)
#   FNDR   truly changed PASs not called / not called (lower is better)
#   power  truly changed PASs called / truly changed  (higher is better)
#   false calls per null dataset                       (lower is better)
#   fdrtool empirical-null sd and eta0 per view        (sd near 1 = well scaled)
#   RMSE of log final dispersion against the truth     (lower is better)
#
# Runtime is set by the number of simulations and PASs:
#   SEQUENCE_CAL_SIMS          datasets per scenario per comparison (default 2)
#   SEQUENCE_CAL_PAS           PASs sampled per comparison (default 10000)
#   SEQUENCE_CAL_DE_FRACTION   truly changed fraction (default 0.10)
#   SEQUENCE_CORES             cores for the likelihood scan (default 2)
#   SEQUENCE_CAL_SCENARIOS     comma-separated subset to run (default: all five).
#                              Rows already in the output tables for scenarios
#                              not rerun are kept, so a subset can be added to an
#                              earlier run.
# With the defaults (40 datasets): about 3.5 hours on 2 cores.
# The Fixed 5,000 cutoff is scaled to the subsample: k = 5000 x n_sub / n_PAS.
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(DESeq2); library(apeglm); library(fdrtool); library(ggplot2); library(dplyr)
  library(grid); library(scales)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
pipeline_path <- Sys.getenv("SEQUENCE_PIPELINE", unset = file.path(script_dir, "SEQUENCE_FINAL_MANUSCRIPT.R"))
if (!file.exists(pipeline_path)) stop("Cannot find the manuscript pipeline: ", pipeline_path)

# Load every function and every constant (ALL_CAPS, plus lfc_boundary and
# track_short) from the pipeline without running it.
load_pipeline_definitions <- function(path) {
  for (e in parse(path, keep.source = FALSE)) {
    if (!is.call(e) || !identical(e[[1]], as.name("<-")) || !is.name(e[[2]])) next
    nm <- as.character(e[[2]]); rhs <- e[[3]]
    is_fun <- is.call(rhs) && identical(rhs[[1]], as.name("function"))
    is_const <- grepl("^[A-Z][A-Z0-9_]*$", nm) || nm %in% c("lfc_boundary", "track_short")
    if (is_fun || is_const) tryCatch(eval(e, globalenv()), error = function(err) NULL)
  }
}
load_pipeline_definitions(pipeline_path)

env_int <- function(name, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(v) || v < 1L) default else v
}
N_SIMS <- env_int("SEQUENCE_CAL_SIMS", 2L)
N_PAS <- env_int("SEQUENCE_CAL_PAS", 10000L)
DE_FRACTION <- suppressWarnings(as.numeric(Sys.getenv("SEQUENCE_CAL_DE_FRACTION", unset = "0.10")))
if (!is.finite(DE_FRACTION) || DE_FRACTION <= 0 || DE_FRACTION >= 1) DE_FRACTION <- 0.10
CORES <- analysis_cores()
ALL_SCENARIOS <- c("one_population", "two_populations", "null", "structured", "structured_null")
SCENARIOS <- trimws(strsplit(Sys.getenv("SEQUENCE_CAL_SCENARIOS", unset = paste(ALL_SCENARIOS, collapse = ",")), ",")[[1]])
if (!length(SCENARIOS) || !all(SCENARIOS %in% ALL_SCENARIOS)) {
  stop("SEQUENCE_CAL_SCENARIOS must be a comma-separated subset of: ", paste(ALL_SCENARIOS, collapse = ", "))
}
NULL_SCENARIOS <- c("null", "structured_null")
STRUCTURED_SCENARIOS <- c("structured", "structured_null")
TWO_POP_SCENARIOS <- c("two_populations", "structured", "structured_null")
CAL_METHODS <- c(ORIGINAL_LABEL, "EVS, Fixed 5,000", "EVS, Likelihood k*")
CAL_COLS <- c("Original (No EVS)" = "#2B2B2B", "EVS, Fixed 5,000" = "#9A9A9A", "EVS, Likelihood k*" = "#2E7D5B")

# ---- data -----------------------------------------------------------------------

count_name <- "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
roots <- unique(c(Sys.getenv("SEQUENCE_REPO_ROOT", unset = ""), getwd(), dirname(getwd()), script_dir,
                  "/root/REAPER98632"))
roots <- roots[nzchar(roots)]
count_file <- NULL
for (r in roots) for (f in c(file.path(r, count_name), file.path(r, "data", count_name))) {
  if (is.null(count_file) && file.exists(f)) count_file <- normalizePath(f)
}
if (is.null(count_file)) stop("Cannot find ", count_name, ". Set SEQUENCE_REPO_ROOT.")
repo_root <- dirname(count_file)
if (basename(repo_root) == "data") repo_root <- dirname(repo_root)
out_dir <- Sys.getenv("SEQUENCE_CAL_OUT", unset = file.path(repo_root, "exports", "sequence_fdr_fndr_calibration"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wtts <- utils::read.csv(count_file, check.names = FALSE)
wtts <- wtts[!is.na(wtts$OrigID) & nzchar(trimws(as.character(wtts$OrigID))), , drop = FALSE]
rownames(wtts) <- make.unique(as.character(wtts$OrigID), sep = "_dup")
CAL_COMPARISONS <- list(RT0_ZT6 = c("R0", "ZT6"), RT2_ZT8 = c("R2", "ZT8"),
                        RT4_ZT10 = c("R4", "ZT10"), RT8_ZT14 = c("R8", "ZT14"))

comparison_counts <- function(prefixes) {
  ids <- c(paste0(prefixes[1], "_", 1:5), paste0(prefixes[2], "_", 1:5))
  x <- vapply(ids, function(s) suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(wtts[[s]]))))),
              numeric(nrow(wtts)))
  rownames(x) <- rownames(wtts)
  x <- x[stats::complete.cases(x), , drop = FALSE]
  x <- coerce_raw_count_matrix_for_deseq2(x, context = "calibration counts")
  x[rowSums(x) > 0, , drop = FALSE]
}

make_coldata <- function(ids) {
  cd <- data.frame(condition = factor(rep(c("trt", "untrt"), each = 5), levels = c("untrt", "trt")))
  rownames(cd) <- ids
  cd
}


# Real sample structure: within each condition, the rank-1 (PC1) component of
# centred log normalized counts, as a multiplier with mean 1 per PAS and
# condition; v = its squared coefficient of variation per PAS.
structure_multipliers <- function(raw, sf, condition) {
  nrm <- sweep(raw, 2, sf, "/")
  M <- matrix(1, nrow(raw), ncol(raw), dimnames = dimnames(raw))
  v <- numeric(nrow(raw))
  levs <- unique(as.character(condition))
  for (cl in levs) {
    cols <- which(condition == cl)
    X <- log(nrm[, cols, drop = FALSE] + 0.5)
    X <- X - rowMeans(X)
    s1 <- svd(X, nu = 1, nv = 1)
    m <- exp(s1$d[1] * (s1$u[, 1] %o% s1$v[, 1]))
    m <- m / rowMeans(m)
    M[, cols] <- m
    v <- v + rowMeans((m - 1)^2) / length(levs)
  }
  list(M = M, v = stats::setNames(v, rownames(raw)))
}

# ---- truth, fitted once per comparison on the real (subsampled) counts --------

fit_truth <- function(raw, coldata, label) {
  cond <- as.character(coldata$condition)
  sf <- DESeq2::estimateSizeFactorsForMatrix(raw)
  q <- rowMeans(sweep(raw, 2, sf, "/"))
  orig <- fit_dispersion_group(raw, cond)
  scan <- likelihood_k_scan(raw, coldata, paste(label, "(real data)"), cores = CORES)
  split <- fit_split_model(raw, cond, scan$depth[rownames(raw)], scan$k_star, keep_fits = TRUE)
  le <- scan$depth[rownames(raw)] <= scan$k_star
  centre_two <- numeric(nrow(raw)); names(centre_two) <- rownames(raw)
  centre_two[le] <- split$LE$centre
  centre_two[!le] <- split$REM$centre
  st <- structure_multipliers(raw, sf, cond)
  list(sf = sf, q = q, in_le = le, M = st$M, v = st$v,
       one = list(centre = orig$centre, tau2 = rep(orig$tau2, nrow(raw))),
       two = list(centre = centre_two, tau2 = ifelse(le, split$LE$tau2, split$REM$tau2)),
       k_star_real = scan$k_star, tau2 = c(Original = orig$tau2, LE = split$LE$tau2, REM = split$REM$tau2))
}

simulate_counts <- function(truth, scenario, ids) {
  disp_model <- if (scenario %in% TWO_POP_SCENARIOS) truth$two else truth$one
  n <- length(truth$q)
  log_alpha <- disp_model$centre + stats::rnorm(n, 0, sqrt(disp_model$tau2))
  alpha_nb <- exp(log_alpha)
  mult <- 1
  if (scenario %in% STRUCTURED_SCENARIOS) {
    alpha_nb <- pmax((alpha_nb - truth$v) / (1 + truth$v), 0.05 * alpha_nb)
    mult <- truth$M
  }
  lfc <- numeric(n)
  if (!scenario %in% NULL_SCENARIOS) {
    de <- sample.int(n, round(DE_FRACTION * n))
    lfc[de] <- sample(c(-1, 1), length(de), replace = TRUE) * (0.5 + stats::rgamma(length(de), shape = 2, scale = 0.5))
  }
  trt <- rep(c(TRUE, FALSE), each = 5)
  mu <- outer(truth$q, truth$sf) * ifelse(matrix(trt, n, 10, byrow = TRUE), 2^(lfc / 2), 2^(-lfc / 2)) * mult
  y <- matrix(stats::rnbinom(n * 10, size = rep(1 / alpha_nb, 10), mu = mu), n, 10)
  dimnames(y) <- list(names(truth$q), ids)
  storage.mode(y) <- "integer"
  keep <- rowSums(y) > 0
  list(y = y[keep, , drop = FALSE], lfc = lfc[keep], log_alpha = log_alpha[keep], in_le = truth$in_le[keep])
}

# ---- one simulated dataset: run all three methods and score them --------------

score_calls <- function(called, changed) {
  n_call <- sum(called); n_not <- sum(!called)
  data.frame(n_called = n_call,
             FDR = if (n_call) sum(called & !changed) / n_call else 0,
             FNDR = if (n_not) sum(!called & changed) / n_not else 0,
             power = if (any(changed)) sum(called & changed) / sum(changed) else NA_real_,
             false_calls = sum(called & !changed))
}

run_simulation <- function(sim, coldata, k_fixed, label) {
  raw <- sim$y; ids <- rownames(raw)
  annot <- data.frame(feature_id = ids, gene_symbol = NA_character_)
  fit_view <- function(rows, view) {
    r <- run_core_analysis(raw[rows, , drop = FALSE], coldata, paste(label, view), annot)
    list(res = r$results, param = r$fdrtool_param)
  }
  depth <- evs_depth(raw, coldata)
  k_lik <- likelihood_k_scan(raw, coldata, label, cores = CORES)$k_star
  views <- list(list(method = ORIGINAL_LABEL, set = ORIGINAL_LABEL, rows = ids))
  for (m in list(c("EVS, Fixed 5,000", k_fixed), c("EVS, Likelihood k*", k_lik))) {
    le <- depth[ids] <= as.integer(m[2])
    views[[length(views) + 1L]] <- list(method = m[1], set = "Leading Edge", rows = ids[le])
    views[[length(views) + 1L]] <- list(method = m[1], set = "Remainder", rows = ids[!le])
  }
  fits <- lapply(views, function(v) tryCatch(fit_view(v$rows, paste(v$method, v$set)), error = function(e) {
    message("  ", label, " ", v$method, " ", v$set, " failed: ", conditionMessage(e)); NULL }))
  changed <- stats::setNames(sim$lfc != 0, ids)
  true_alpha <- stats::setNames(sim$log_alpha, ids)
  true_le <- stats::setNames(sim$in_le, ids)
  calls <- list(); nulls <- list(); disp <- list()
  for (m in CAL_METHODS) {
    idx <- which(vapply(views, function(v) v$method == m, logical(1)))
    in_le <- stats::setNames(rep(NA, length(ids)), ids)
    if (m != ORIGINAL_LABEL) in_le[unlist(lapply(views[idx], function(v) if (v$set == "Leading Edge") v$rows))] <- TRUE
    if (m != ORIGINAL_LABEL) in_le[is.na(in_le)] <- FALSE
    if (any(vapply(fits[idx], is.null, logical(1)))) next
    res <- do.call(rbind, lapply(fits[idx], function(f) f$res[, c("feature_id", "standard_flag", "hbfss_flag", "dispersion")]))
    res <- res[match(ids, res$feature_id), ]
    for (call in c("Standard", "HBFSS")) {
      flag <- if (call == "Standard") res$standard_flag else res$hbfss_flag
      flag[is.na(flag)] <- FALSE
      by_group <- function(g, truth_value) if (m == ORIGINAL_LABEL) NA_integer_ else
        sum(flag & changed == truth_value & in_le == (g == "LE"))
      calls[[length(calls) + 1L]] <- data.frame(
        Method = m, Call = call, score_calls(flag, changed),
        k = if (m == "EVS, Fixed 5,000") k_fixed else if (m == "EVS, Likelihood k*") k_lik else NA_integer_,
        false_calls_LE = by_group("LE", FALSE), false_calls_REM = by_group("REM", FALSE),
        true_calls_LE = by_group("LE", TRUE), true_calls_REM = by_group("REM", TRUE))
    }
    for (j in idx) {
      p <- fits[[j]]$param
      nulls[[length(nulls) + 1L]] <- data.frame(Method = m, View = views[[j]]$set,
        null_sd = if (!is.null(p) && "sd" %in% colnames(p)) p[1, "sd"] else NA_real_,
        eta0 = if (!is.null(p) && "eta0" %in% colnames(p)) p[1, "eta0"] else NA_real_)
    }
    err <- log(res$dispersion) - true_alpha
    ok <- is.finite(err)
    for (g in c("Leading Edge", "Remainder")) {
      sel <- ok & (true_le == (g == "Leading Edge"))
      disp[[length(disp) + 1L]] <- data.frame(Method = m, True_group = g, rmse_log_dispersion = sqrt(mean(err[sel]^2)))
    }
  }
  list(calls = do.call(rbind, calls), nulls = do.call(rbind, nulls), disp = do.call(rbind, disp),
       k = c(fixed = k_fixed, likelihood = k_lik))
}

# ---- run -------------------------------------------------------------------------

all_calls <- list(); all_nulls <- list(); all_disp <- list(); truth_rows <- list()
n_total <- length(CAL_COMPARISONS) * length(SCENARIOS) * N_SIMS
done <- 0L
t0 <- Sys.time()
for (cn in names(CAL_COMPARISONS)) {
  raw_full <- comparison_counts(CAL_COMPARISONS[[cn]])
  ids <- colnames(raw_full)
  coldata <- make_coldata(ids)
  set.seed(20261031L)
  raw <- raw_full[sort(sample.int(nrow(raw_full), min(N_PAS, nrow(raw_full)))), , drop = FALSE]
  k_fixed <- as.integer(round(5000 * nrow(raw) / nrow(raw_full)))
  message("\n", cn, ": fitting the truth on ", format(nrow(raw), big.mark = ","), " real PASs")
  truth <- fit_truth(raw, coldata, cn)
  truth_rows[[cn]] <- data.frame(Comparison = cn, n_PAS = nrow(raw), k_fixed_scaled = k_fixed,
                                 k_star_real = truth$k_star_real, tau2_Original = truth$tau2[["Original"]],
                                 tau2_LE = truth$tau2[["LE"]], tau2_REM = truth$tau2[["REM"]])
  for (sc in SCENARIOS) for (s in seq_len(N_SIMS)) {
    set.seed(20261031L + 1000L * match(cn, names(CAL_COMPARISONS)) + 100L * match(sc, SCENARIOS) + s)
    label <- sprintf("%s %s sim %d", cn, sc, s)
    sim <- simulate_counts(truth, sc, ids)
    out <- tryCatch(run_simulation(sim, coldata, k_fixed, label), error = function(e) {
      message("  ", label, " failed: ", conditionMessage(e)); NULL })
    done <- done + 1L
    if (!is.null(out)) {
      tag <- data.frame(Comparison = cn, Scenario = sc, Sim = s)
      all_calls[[label]] <- cbind(tag, out$calls)
      all_nulls[[label]] <- cbind(tag, out$nulls)
      all_disp[[label]] <- cbind(tag, out$disp)
    }
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("  [%d / %d] %s done; %.0f min elapsed, about %.0f min left", done, n_total, label,
                    el, el / done * (n_total - done)))
  }
}

calls <- do.call(rbind, all_calls); nulls <- do.call(rbind, all_nulls); disp <- do.call(rbind, all_disp)
truth_tab <- do.call(rbind, truth_rows)
if (is.null(calls) || !nrow(calls)) stop("No simulation completed.")

# Keep earlier rows for scenarios not rerun this time.
keep_earlier <- function(new, file) {
  f <- file.path(out_dir, file)
  if (!file.exists(f)) return(new)
  old <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Scenario" %in% names(old)) return(new)
  dplyr::bind_rows(old[!old$Scenario %in% SCENARIOS, , drop = FALSE], new)
}
calls <- keep_earlier(calls, "Table_S2_Calls_By_Simulation.csv")
nulls <- keep_earlier(nulls, "Table_S2_fdrtool_By_View.csv")
disp <- keep_earlier(disp, "Table_S2_Dispersion_Error.csv")
utils::write.csv(truth_tab, file.path(out_dir, "Table_S2_Truth_Parameters.csv"), row.names = FALSE)
utils::write.csv(calls, file.path(out_dir, "Table_S2_Calls_By_Simulation.csv"), row.names = FALSE)
utils::write.csv(nulls, file.path(out_dir, "Table_S2_fdrtool_By_View.csv"), row.names = FALSE)
utils::write.csv(disp, file.path(out_dir, "Table_S2_Dispersion_Error.csv"), row.names = FALSE)

mean_ci <- function(d, value, by) {
  out <- do.call(rbind, lapply(split(d, d[by], drop = TRUE), function(x) {
    v <- x[[value]][is.finite(x[[value]])]
    se <- if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA_real_
    cbind(x[1, by, drop = FALSE], mean = if (length(v)) mean(v) else NA_real_, lo = mean(v) - 1.96 * se,
          hi = mean(v) + 1.96 * se, n_sims = length(v))
  }))
  rownames(out) <- NULL
  out
}
summary_tab <- rbind(
  cbind(Metric = "FDR", mean_ci(calls, "FDR", c("Scenario", "Method", "Call"))),
  cbind(Metric = "FNDR", mean_ci(calls, "FNDR", c("Scenario", "Method", "Call"))),
  cbind(Metric = "Power", mean_ci(calls, "power", c("Scenario", "Method", "Call"))),
  cbind(Metric = "False calls", mean_ci(calls, "false_calls", c("Scenario", "Method", "Call"))))
group_tab <- do.call(rbind, lapply(c("false_calls_LE", "false_calls_REM", "true_calls_LE", "true_calls_REM"), function(v)
  cbind(Metric = v, mean_ci(calls[calls$Method != ORIGINAL_LABEL, ], v, c("Scenario", "Method", "Call")))))
utils::write.csv(group_tab, file.path(out_dir, "Table_S2_Calls_By_Group.csv"), row.names = FALSE)
utils::write.csv(summary_tab, file.path(out_dir, "Table_S2_Calibration_Summary.csv"), row.names = FALSE)

# ---- Figure S2 ----------------------------------------------------------------------

scen_lab <- c(one_population = "one pop.", two_populations = "two pops.", structured = "structured")
scen_lab <- scen_lab[names(scen_lab) %in% unique(calls$Scenario)]
bar_panel <- function(metric, tag, title, subtitle, ylab, pct = TRUE) {
  d <- summary_tab[summary_tab$Metric == metric & summary_tab$Scenario %in% names(scen_lab), ]
  d$Group <- factor(paste0(d$Call, "\n(", scen_lab[d$Scenario], ")"),
                    levels = as.vector(t(outer(c("Standard", "HBFSS"), paste0("\n(", scen_lab, ")"), paste0))))
  d$Method <- factor(d$Method, levels = CAL_METHODS)
  ggplot(d, aes(Group, mean, fill = Method)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.76) +
    geom_errorbar(aes(ymin = pmax(lo, 0), ymax = hi), position = position_dodge(width = 0.82), width = 0,
                  linewidth = 0.35, na.rm = TRUE) +
    scale_fill_manual(values = CAL_COLS) +
    scale_y_continuous(labels = if (pct) scales::label_percent(accuracy = 1) else scales::label_comma(),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = ylab, tag = tag, title = title, subtitle = subtitle) +
    fig_theme() + theme(panel.grid.major.x = element_blank(), legend.position = "none")
}
pa <- bar_panel("FDR", "a", "False discoveries (lower is better)", "Of the sites called, the fraction truly unchanged.",
                "False discovery rate")
pb <- bar_panel("FNDR", "b", "Missed sites (lower is better)", "Of the sites not called, the fraction truly changed.",
                "False non-discovery rate")
pc <- bar_panel("Power", "c", "True sites found (higher is better)", "Fraction of truly changed sites called.", "Power")

dn <- summary_tab[summary_tab$Metric == "False calls" & summary_tab$Scenario %in% NULL_SCENARIOS, ]
null_lab <- c(null = "no structure", structured_null = "structured")
dn$Group <- factor(paste0(dn$Call, "\n(", null_lab[dn$Scenario], ")"),
                   levels = as.vector(t(outer(c("Standard", "HBFSS"), paste0("\n(", null_lab, ")"), paste0))))
dn$Method <- factor(dn$Method, levels = CAL_METHODS)
pd <- ggplot(dn, aes(Group, mean, fill = Method)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.76) +
  geom_errorbar(aes(ymin = pmax(lo, 0), ymax = hi), position = position_dodge(width = 0.82), width = 0,
                linewidth = 0.35, na.rm = TRUE) +
  scale_fill_manual(values = CAL_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "False calls per dataset", tag = "d", title = "Nothing truly changed (lower is better)",
       subtitle = "Every call here is a false positive.") +
  fig_theme() + theme(panel.grid.major.x = element_blank(), legend.position = "none")

ref_scen <- if ("structured" %in% nulls$Scenario) "structured" else "two_populations"
ne <- nulls[nulls$Method %in% c(ORIGINAL_LABEL, "EVS, Likelihood k*") & nulls$Scenario == ref_scen, ]
se_tab <- mean_ci(ne, "null_sd", "View")
eta_tab <- mean_ci(ne, "eta0", "View")
se_tab$View <- factor(se_tab$View, levels = names(SET_COLS))
pe <- ggplot(se_tab, aes(View, mean, fill = View)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.35, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = 0.35, colour = "grey30") +
  scale_fill_manual(values = SET_COLS, guide = "none") +
  coord_cartesian(ylim = c(min(0.8, min(se_tab$lo, na.rm = TRUE)), max(1.15, max(se_tab$hi, na.rm = TRUE)))) +
  labs(x = NULL, y = "fdrtool null spread (σ)", tag = "e", title = "Test statistics well scaled? (near 1 is better)",
       subtitle = sprintf("Likelihood k*, %s truth. fdrtool null fraction %s (truth %.2f).",
                          if (ref_scen == "structured") "structured" else "two-population",
                          paste(sprintf("%.2f", eta_tab$mean[match(names(SET_COLS), eta_tab$View)]), collapse = " / "),
                          1 - DE_FRACTION)) +
  fig_theme() + theme(panel.grid.major.x = element_blank())

df_tab <- mean_ci(disp[disp$Scenario == ref_scen & disp$Method %in% c(ORIGINAL_LABEL, "EVS, Likelihood k*"), ],
                  "rmse_log_dispersion", c("Method", "True_group"))
df_tab$Method <- factor(df_tab$Method, levels = CAL_METHODS)
pf <- ggplot(df_tab, aes(True_group, mean, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.8), width = 0, linewidth = 0.35,
                na.rm = TRUE) +
  scale_fill_manual(values = CAL_COLS) +
  scale_x_discrete(labels = function(x) paste(x, "PASs")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Dispersion error\n(RMSE of log final dispersion)", tag = "f",
       title = "Dispersion error (lower is better)",
       subtitle = sprintf("Distance from the true dispersions (%s truth).", if (ref_scen == "structured") "structured" else "two-population")) +
  fig_theme() + theme(panel.grid.major.x = element_blank(), legend.position = "none")

pa <- pa + theme(legend.position = "top", legend.justification = "left")
fig <- assemble_rows(list(list(plots = list(pa, pb, pc), widths = c(1, 1, 1)),
                          list(plots = list(pd, pe, pf), widths = c(1, 1, 1))),
                     heights = c(1, 1),
                     title = "Supplementary Figure S2. FDR and FNDR calibration (simulations with known truth)")
save_figure(fig, out_dir, "Figure_S2_FDR_FNDR_Calibration", 200, 150)

cat("\n=====================================================\n")
cat("Calibration complete:", nrow(calls) / (2 * length(CAL_METHODS)), "simulated datasets in the tables.\n")
print(summary_tab[summary_tab$Metric %in% c("FDR", "Power", "False calls"), c("Metric", "Scenario", "Method", "Call", "mean")], digits = 3)
print(group_tab[, c("Metric", "Scenario", "Method", "Call", "mean")], digits = 3)
cat("Output: ", normalizePath(out_dir), "\n", sep = "")
cat("=====================================================\n")
