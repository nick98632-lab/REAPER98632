# =============================================================================
# SEQUENCE 10 ACTIVE REBUILD
# Current manuscript method only
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(fdrtool)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(gridExtra)
  library(grid)
})

# -----------------------------------------------------------------------------
# paths and settings
# -----------------------------------------------------------------------------

repo_dir <- getwd()
in_dir <- file.path(repo_dir, "data")
out_root <- file.path(repo_dir, "exports")
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(out_root, paste0("run_", stamp))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

count_file <- file.path(in_dir, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
if (!file.exists(count_file)) stop("count file not found: ", count_file, call. = FALSE)

alpha_level <- 0.20
lfc_cut <- 1.0
hc_cap <- 0.99
shrink_type <- "apeglm"
wald_coarse_step <- 25L
wald_fine_radius <- 50L
wald_min_n <- 200L
w_mean <- 1.0
w_sd <- 1.0
w_skew <- 0.5
wave_starts <- 10L
fig_dpi <- 320
base_size <- 10
hist_bins <- 60L

pal <- list(
  thr = "#8C2D04",
  hbfss = "#E67E22",
  strong = "#C0392B",
  weak = "#4A90E2",
  overlap = "#7D3C98",
  grey = "#B3B3B3",
  trt = "#1F78B4",
  ctrl = "#4D4D4D"
)

# -----------------------------------------------------------------------------
# metadata
# -----------------------------------------------------------------------------

meta_all <- data.frame(
  id = c(
    "R0_1","R0_2","R0_3","R0_4","R0_5","ZT6_1","ZT6_2","ZT6_3","ZT6_4","ZT6_5",
    "R2_1","R2_2","R2_3","R2_4","R2_5","ZT8_1","ZT8_2","ZT8_3","ZT8_4","ZT8_5",
    "R4_1","R4_2","R4_3","R4_4","R4_5","ZT10_1","ZT10_2","ZT10_3","ZT10_4","ZT10_5",
    "R8_1","R8_2","R8_3","R8_4","R8_5","ZT14_1","ZT14_2","ZT14_3","ZT14_4","ZT14_5"
  ),
  cond = c(
    rep("trt", 5), rep("untrt", 5),
    rep("trt", 5), rep("untrt", 5),
    rep("trt", 5), rep("untrt", 5),
    rep("trt", 5), rep("untrt", 5)
  ),
  stringsAsFactors = FALSE
)
rownames(meta_all) <- meta_all$id
meta_all$cond <- factor(meta_all$cond, levels = c("untrt", "trt"))

cmp_tbl <- data.frame(
  cmp = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  trt = c("R0", "R2", "R4", "R8"),
  ctrl = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# helpers
# -----------------------------------------------------------------------------

log_step <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
  flush.console()
}

fmt <- function(x, d = 4) {
  if (!is.finite(x) || is.na(x)) return("NA")
  formatC(x, digits = d, format = "fg", flag = "#")
}

assert_cols <- function(df, cols, obj = "data frame") {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop(obj, " missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
}

clip_p <- function(x, eps = 1e-300) {
  x <- as.numeric(x)
  x[is.finite(x)] <- pmin(pmax(x[is.finite(x)], eps), 1 - 1e-12)
  x[!is.finite(x)] <- NA_real_
  x
}

safe_skew <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) < 5L) return(NA_real_)
  sx <- stats::sd(x)
  if (!is.finite(sx) || sx <= 0) return(NA_real_)
  mean(((x - mean(x)) / sx)^3)
}

theme_seq <- function() {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1),
      plot.caption = element_text(size = base_size - 3),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88")
    )
}

save_plot <- function(p, path, w = 12, h = 8) {
  ggsave(path, p, width = w, height = h, dpi = fig_dpi, units = "in", limitsize = FALSE, bg = "white")
}

save_grob <- function(g, path, w = 14, h = 10) {
  png(path, width = w, height = h, units = "in", res = fig_dpi)
  grid::grid.newpage(); grid::grid.draw(g)
  dev.off()
}

roll_mean_dense <- function(x, k = 151L) {
  x <- as.numeric(x)
  n <- length(x)
  if (!n) return(numeric(0))
  k <- max(3L, as.integer(k))
  if (k %% 2L == 0L) k <- k + 1L
  h <- (k - 1L) / 2L
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    lo <- max(1L, i - h)
    hi <- min(n, i + h)
    vals <- x[lo:hi]
    vals <- vals[is.finite(vals) & !is.na(vals)]
    if (length(vals)) out[i] <- mean(vals)
  }
  out
}

# -----------------------------------------------------------------------------
# load data
# -----------------------------------------------------------------------------

wt <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
assert_cols(wt, c("OrigID", "Symbol"), "count file")
assert_cols(wt, meta_all$id, "count file")
wt$OrigID <- as.character(wt$OrigID)
wt$Symbol <- as.character(wt$Symbol)
wt <- wt[!is.na(wt$OrigID), , drop = FALSE]
wt <- wt[rowSums(is.na(wt[, meta_all$id, drop = FALSE])) == 0, , drop = FALSE]
rownames(wt) <- wt$OrigID

annot <- unique(wt[, c("OrigID", "Symbol"), drop = FALSE])
names(annot) <- c("feature_id", "gene")
annot$feature_id <- as.character(annot$feature_id)
annot$gene <- as.character(annot$gene)

# -----------------------------------------------------------------------------
# DE / Strimmer
# -----------------------------------------------------------------------------

get_coef <- function(dds) {
  rn <- resultsNames(dds)
  idx <- grep("^cond_", rn)
  if (!length(idx)) stop("condition coefficient not found", call. = FALSE)
  rn[idx[1]]
}

run_fdr <- function(stat, ds_name) {
  z <- as.numeric(stat)
  z <- z[is.finite(z) & !is.na(z)]
  if (length(z) < 5L) stop("[", ds_name, "] too few finite stats for fdrtool", call. = FALSE)
  ft <- suppressWarnings(suppressMessages(fdrtool(z, statistic = "normal", plot = FALSE, verbose = FALSE)))
  p <- clip_p(ft$pval)
  hc <- suppressWarnings(fdrtool::hc.thresh(p))
  hc <- as.numeric(hc)
  if (!is.finite(hc) || is.na(hc)) hc <- NA_real_
  hc <- min(hc, hc_cap, na.rm = TRUE)
  list(ft = ft, p = p, hc = hc)
}

run_ds <- function(cnt, meta, ds_name) {
  dds <- DESeqDataSetFromMatrix(round(as.matrix(cnt)), colData = meta, design = ~ cond)
  dds <- suppressMessages(DESeq(dds, quiet = TRUE))
  coef_name <- get_coef(dds)
  res <- suppressMessages(results(dds, alpha = alpha_level, independentFiltering = FALSE))
  shr <- suppressMessages(suppressWarnings(lfcShrink(dds, coef = coef_name, type = shrink_type, res = res)))
  out <- as.data.frame(res)
  out$lfc <- as.numeric(shr$log2FoldChange)
  out$feature_id <- rownames(out)
  out$baseMean <- as.numeric(out$baseMean)
  out$stat <- as.numeric(out$stat)
  out$pvalue <- as.numeric(out$pvalue)
  out$padj <- as.numeric(out$padj)
  out$dispGeneEst <- mcols(dds)$dispGeneEst[match(out$feature_id, rownames(dds))]
  out$dispGeneEst <- as.numeric(out$dispGeneEst)

  fd <- run_fdr(out$stat, ds_name)
  stat_all <- as.numeric(out$stat)
  keep <- is.finite(stat_all) & !is.na(stat_all)
  out$emp_p <- NA_real_
  out$emp_q <- NA_real_
  out$lfdr <- NA_real_
  out$emp_p[keep] <- fd$ft$pval
  out$emp_q[keep] <- fd$ft$qval
  out$lfdr[keep] <- fd$ft$lfdr
  out$hc <- fd$hc
  out$hbfss_thr <- if (is.finite(fd$hc) && !is.na(fd$hc)) abs(log10(fd$hc)) * lfc_cut else NA_real_
  out$hbfss <- abs(out$lfc) * abs(log10(pmax(out$emp_p, 1e-300)))

  out$std <- !is.na(out$padj) & out$padj < alpha_level
  out$weak <- out$std & abs(out$lfc) < lfc_cut
  out$strong <- out$std & abs(out$lfc) >= lfc_cut
  out$hbfss_sig <- if (is.finite(fd$hc) && !is.na(fd$hc)) {
    !is.na(out$emp_p) & out$emp_p <= fd$hc & !is.na(out$hbfss_thr) & out$hbfss >= out$hbfss_thr
  } else {
    rep(FALSE, nrow(out))
  }
  out$hbfss_only <- out$hbfss_sig & !out$std
  out$overlap <- out$hbfss_sig & out$std
  out$dir <- ifelse(out$lfc > 0, "up", ifelse(out$lfc < 0, "down", "flat"))
  out <- dplyr::left_join(out, annot, by = "feature_id")
  out
}

# -----------------------------------------------------------------------------
# EVS ranking
# -----------------------------------------------------------------------------

get_ids <- function(prefix) grep(paste0("^", prefix, "_"), meta_all$id, value = TRUE)

get_loads <- function(cnt) {
  x <- log2(as.matrix(cnt) + 1)
  pc <- prcomp(t(x), center = TRUE, scale. = FALSE)
  load <- pc$rotation[, 1]
  data.frame(feature_id = rownames(x), load = as.numeric(load), stringsAsFactors = FALSE)
}

build_rank_tbl <- function(cnt, trt_ids, ctrl_ids) {
  trt <- get_loads(cnt[, trt_ids, drop = FALSE])
  ctrl <- get_loads(cnt[, ctrl_ids, drop = FALSE])
  names(trt)[2] <- "load_trt"
  names(ctrl)[2] <- "load_ctrl"
  x <- full_join(trt, ctrl, by = "feature_id")
  x$load_trt <- as.numeric(x$load_trt)
  x$load_ctrl <- as.numeric(x$load_ctrl)
  x$load_max <- pmax(abs(x$load_trt), abs(x$load_ctrl), na.rm = TRUE)
  x <- x[order(-x$load_max, x$feature_id), , drop = FALSE]
  x$rank <- seq_len(nrow(x))
  x$rank_trt <- rank(-abs(x$load_trt), ties.method = "first", na.last = "keep")
  x$rank_ctrl <- rank(-abs(x$load_ctrl), ties.method = "first", na.last = "keep")
  x
}

# -----------------------------------------------------------------------------
# Damped wave method
# -----------------------------------------------------------------------------

build_rank_nb <- function(rank_tbl, raw_ref) {
  x <- left_join(rank_tbl, raw_ref[, c("feature_id", "baseMean", "dispGeneEst", "stat")], by = "feature_id")
  x$mu <- as.numeric(x$baseMean)
  x$alpha <- as.numeric(x$dispGeneEst)
  x$stat <- as.numeric(x$stat)
  x <- x[is.finite(x$rank) & is.finite(x$mu) & x$mu > 0 & is.finite(x$alpha) & x$alpha >= 0, , drop = FALSE]
  x$x <- (x$rank - min(x$rank)) / (max(x$rank) - min(x$rank))
  x$iod <- 1 + x$alpha * x$mu
  x$cv2 <- (1 / x$mu) + x$alpha
  x$delta <- x$iod - x$cv2
  x
}

wave_fun <- function(par, x) {
  b0 <- par[1]; b1 <- par[2]; b2 <- par[3]
  A <- par[4]; lam <- par[5]; f <- par[6]; phi <- par[7]
  b0 + b1 * x + b2 * x^2 + A * exp(-lam * x) * cos(2 * pi * f * x + phi)
}

init_wave <- function(x, y) {
  fit <- lm(y ~ x + I(x^2))
  b <- coef(fit)
  b0 <- ifelse(is.finite(b[1]), b[1], mean(y, na.rm = TRUE))
  b1 <- ifelse(length(b) >= 2 && is.finite(b[2]), b[2], 0)
  b2 <- ifelse(length(b) >= 3 && is.finite(b[3]), b[3], 0)
  res <- y - (b0 + b1 * x + b2 * x^2)
  A0 <- 0.5 * diff(range(res, na.rm = TRUE))
  if (!is.finite(A0) || A0 <= 0) A0 <- stats::sd(res, na.rm = TRUE)
  if (!is.finite(A0) || A0 <= 0) A0 <- 0.1
  s <- sign(res)
  s[!is.finite(s)] <- 0
  zc <- sum(abs(diff(s)) > 0)
  f0 <- max(1, min(25, zc / 2))
  c(b0, b1, b2, A0, 1, f0, 0)
}

fit_wave <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 30L) stop("too few points for damped wave fit", call. = FALSE)
  base <- init_wave(x, y)
  ysd <- stats::sd(y, na.rm = TRUE); if (!is.finite(ysd) || ysd <= 0) ysd <- 1
  obj <- function(par) {
    fit <- wave_fun(par, x)
    if (any(!is.finite(fit))) return(1e30)
    pen <- 0
    if (par[4] < 0) pen <- pen + 1e6 * abs(par[4])
    if (par[5] < 0) pen <- pen + 1e6 * abs(par[5])
    if (par[6] <= 0) pen <- pen + 1e6 * abs(par[6])
    sum((y - fit)^2) + pen
  }
  best <- NULL; best_val <- Inf
  set.seed(1)
  for (i in seq_len(wave_starts)) {
    j <- c(rnorm(1, 0, 0.15 * ysd), rnorm(1, 0, 0.5 * ysd), rnorm(1, 0, 0.5 * ysd),
           abs(rnorm(1, 0, 0.5 * abs(base[4]))), abs(rnorm(1, 0, 1.2)), abs(rnorm(1, 0, 2.0)), runif(1, -pi, pi))
    p0 <- base + j
    p0[4] <- max(1e-6, abs(p0[4])); p0[5] <- max(1e-6, abs(p0[5])); p0[6] <- max(1e-6, abs(p0[6]))
    fit_try <- try(optim(p0, obj, method = "Nelder-Mead", control = list(maxit = 20000, reltol = 1e-12)), silent = TRUE)
    if (inherits(fit_try, "try-error")) next
    if (is.finite(fit_try$value) && fit_try$value < best_val) { best <- fit_try; best_val <- fit_try$value }
  }
  if (is.null(best)) stop("damped wave fit failed", call. = FALSE)
  fit <- wave_fun(best$par, x)
  rss <- sum((y - fit)^2)
  tss <- sum((y - mean(y))^2)
  r2 <- if (tss > 0) 1 - rss / tss else NA_real_
  list(par = best$par, fitted = fit, r2 = r2)
}

predict_wave <- function(fit, x) wave_fun(fit$par, x)

find_cross <- function(rank_df, iod_hat, cv2_hat) {
  d <- iod_hat - cv2_hat
  zidx <- which(abs(d) < 1e-8)
  if (length(zidx)) return(as.integer(rank_df$rank[zidx[1]]))
  s <- sign(d)
  idx <- which(s[-1] != s[-length(s)])
  if (!length(idx)) return(as.integer(rank_df$rank[which.min(abs(d))]))
  i <- idx[1]
  as.integer(round(mean(c(rank_df$rank[i], rank_df$rank[i + 1L]))))
}

build_ref_profiles <- function(rank_df) {
  # treatment-specific and control-specific profiles use same NB quantities on their own loading ranks
  make_prof <- function(rank_col) {
    d <- rank_df[order(rank_df[[rank_col]]), c("feature_id", "mu", "alpha", "stat", rank_col), drop = FALSE]
    names(d)[names(d) == rank_col] <- "rank"
    d$x <- (d$rank - min(d$rank)) / (max(d$rank) - min(d$rank))
    d$iod <- 1 + d$alpha * d$mu
    d$cv2 <- (1 / d$mu) + d$alpha
    iod_fit <- fit_wave(d$x, d$iod)
    cv2_fit <- fit_wave(d$x, d$cv2)
    iod_hat <- predict_wave(iod_fit, d$x)
    cv2_hat <- predict_wave(cv2_fit, d$x)
    cut <- find_cross(d, iod_hat, cv2_hat)
    list(df = d, iod_fit = iod_fit, cv2_fit = cv2_fit, iod_hat = iod_hat, cv2_hat = cv2_hat, cut = cut)
  }
  list(trt = make_prof("rank_trt"), ctrl = make_prof("rank_ctrl"), comb = make_prof("rank"))
}

# -----------------------------------------------------------------------------
# raw Wald scan
# -----------------------------------------------------------------------------

dist_one <- function(z) {
  z <- as.numeric(z)
  z <- z[is.finite(z) & !is.na(z)]
  if (length(z) < 5L) return(list(total = Inf, m = NA_real_, sd = NA_real_, skew = NA_real_))
  m <- abs(mean(z))
  s <- abs(stats::sd(z) - 1)
  k <- abs(safe_skew(z))
  list(total = w_mean * m + w_sd * s + w_skew * k, m = m, sd = s, skew = k)
}

scan_wald <- function(rank_df, left_rank, right_rank) {
  r <- sort(unique(c(left_rank, right_rank)))
  lo <- max(wald_min_n, r[1])
  hi <- min(nrow(rank_df) - wald_min_n, r[2])
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
    lo <- max(wald_min_n, floor(nrow(rank_df) * 0.10))
    hi <- min(nrow(rank_df) - wald_min_n, ceiling(nrow(rank_df) * 0.25))
  }
  cand <- seq(lo, hi, by = wald_coarse_step)
  eval_rank <- function(k) {
    z1 <- rank_df$stat[seq_len(k)]
    z2 <- rank_df$stat[(k + 1L):nrow(rank_df)]
    d1 <- dist_one(z1); d2 <- dist_one(z2)
    data.frame(rank = k, total = d1$total + d2$total, lead_n = length(z1), rem_n = length(z2), stringsAsFactors = FALSE)
  }
  coarse <- bind_rows(lapply(cand, eval_rank))
  best <- coarse$rank[which.min(coarse$total)]
  fine <- seq(max(lo, best - wald_fine_radius), min(hi, best + wald_fine_radius), by = 1L)
  fine_tbl <- bind_rows(lapply(fine, eval_rank))
  sel <- fine_tbl[which.min(fine_tbl$total), , drop = FALSE]
  list(coarse = coarse, fine = fine_tbl, sel = sel)
}

# -----------------------------------------------------------------------------
# split builder
# -----------------------------------------------------------------------------

build_split <- function(cnt, rank_tbl, final_rank) {
  lead_ids <- rank_tbl$feature_id[(rank_tbl$rank_trt <= final_rank) | (rank_tbl$rank_ctrl <= final_rank)]
  lead_ids <- unique(lead_ids)
  rem_ids <- setdiff(rank_tbl$feature_id, lead_ids)
  list(
    raw = cnt,
    lead = cnt[rownames(cnt) %in% lead_ids, , drop = FALSE],
    rem = cnt[rownames(cnt) %in% rem_ids, , drop = FALSE],
    lead_ids = lead_ids,
    rem_ids = rem_ids
  )
}

# -----------------------------------------------------------------------------
# figures
# -----------------------------------------------------------------------------

plot_cut_panel <- function(cmp, ref, cut_rank, path) {
  d <- ref$comb$df
  iod_hat <- ref$comb$iod_hat
  cv2_hat <- ref$comb$cv2_hat
  delta_hat <- iod_hat - cv2_hat
  old <- par(no.readonly = TRUE); on.exit(par(old), add = TRUE)
  png(path, width = 14, height = 10, units = "in", res = fig_dpi)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(2,2), mar = c(4.2,4.8,3.2,1.2), oma = c(1.3,1.2,2.8,0.3), mgp = c(2.3,0.7,0), tcl = -0.25)

  yr <- range(c(iod_hat, cv2_hat), na.rm = TRUE)
  plot(d$rank, iod_hat, type = "l", lwd = 2, col = pal$trt, xlab = "Rank", ylab = "amplitude", main = "A. Rank-linked regime profiles", ylim = yr)
  lines(d$rank, cv2_hat, lwd = 2, col = pal$hbfss)
  abline(v = cut_rank, lty = 2, lwd = 2, col = pal$trt)
  legend("topright", legend = c("IOD(r)", "CV2(r)"), col = c(pal$trt, pal$hbfss), lwd = 2, bty = "n", cex = 0.9)

  plot(d$rank, iod_hat, type = "l", lwd = 2, col = pal$trt, xlab = "Rank", ylab = "IOD(r)", main = "B. IOD as a function of rank")
  abline(v = cut_rank, lty = 2, lwd = 2, col = pal$trt)
  legend("topleft", legend = c("IOD(r) = Var(Y[r]) / mu[r] = 1 + alpha[r] * mu[r]"), bty = "n", cex = 0.85)

  plot(d$rank, cv2_hat, type = "l", lwd = 2, col = pal$trt, xlab = "Rank", ylab = "CV2(r)", main = "C. CV2 as a function of rank")
  abline(v = cut_rank, lty = 2, lwd = 2, col = pal$trt)
  legend("topleft", legend = c("CV2(r) = Var(Y[r]) / mu[r]^2 = 1 / mu[r] + alpha[r]"), bty = "n", cex = 0.85)

  plot(d$rank, delta_hat, type = "l", lwd = 2, col = pal$trt, xlab = "Rank", ylab = "IOD(r) - CV2(r)", main = "D. Crossing condition")
  abline(h = 0, lty = 3, lwd = 1.5, col = pal$trt)
  abline(v = cut_rank, lty = 2, lwd = 2, col = pal$trt)
  legend("topleft", legend = c("Delta(r) = IOD(r) - CV2(r)", "Crossing occurs where Delta(r) = 0"), bty = "n", cex = 0.85)

  mtext(paste0(cmp, " | cut panel"), outer = TRUE, side = 3, line = 0.9, cex = 1.2, font = 2)
  mtext("Dashed vertical line = projected cutoff across all panels.", outer = TRUE, side = 1, line = 0.1, cex = 0.95)
}

plot_wald_panel <- function(cmp, rank_df, scan, path) {
  sel_rank <- as.integer(scan$sel$rank[1])
  z_raw <- rank_df$stat
  z_lead <- rank_df$stat[seq_len(sel_rank)]
  z_rem <- rank_df$stat[(sel_rank + 1L):nrow(rank_df)]

  hist_df <- function(z, fill) {
    ggplot(data.frame(z = z), aes(z)) +
      geom_histogram(bins = hist_bins, fill = fill, color = "white") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = pal$thr) +
      theme_seq()
  }
  p1 <- hist_df(z_raw, pal$grey) + labs(title = paste0(cmp, " | Raw"), subtitle = paste0("n = ", length(z_raw), " | mean = ", fmt(mean(z_raw)), " | sd = ", fmt(sd(z_raw)), " | skew = ", fmt(safe_skew(z_raw))), x = "Wald statistic", y = "Count")
  p2 <- ggplot(scan$fine, aes(rank, total)) + geom_line(colour = pal$ctrl) + geom_point(size = 0.8, colour = pal$ctrl) + geom_vline(xintercept = sel_rank, colour = pal$thr) + theme_seq() + labs(title = paste0(cmp, " | raw-Wald scan"), subtitle = "Minimum total distortion selects final rank", x = "Candidate rank", y = "Total distortion")
  p3 <- hist_df(z_lead, pal$trt) + labs(title = paste0(cmp, " | Lead"), subtitle = paste0("n = ", length(z_lead), " | mean = ", fmt(mean(z_lead)), " | sd = ", fmt(sd(z_lead)), " | skew = ", fmt(safe_skew(z_lead))), x = "Wald statistic", y = "Count")
  p4 <- hist_df(z_rem, pal$ctrl) + labs(title = paste0(cmp, " | Rem"), subtitle = paste0("n = ", length(z_rem), " | mean = ", fmt(mean(z_rem)), " | sd = ", fmt(sd(z_rem)), " | skew = ", fmt(safe_skew(z_rem))), x = "Wald statistic", y = "Count")
  g <- arrangeGrob(p1, p2, p3, p4, ncol = 2, top = textGrob(paste0(cmp, " | wald panel"), gp = gpar(fontface = "bold", cex = 1.1)))
  save_grob(g, path, 14, 10)
}

make_volc_class <- function(df) {
  cls <- rep("Neither", nrow(df))
  cls[df$weak] <- "DC2 weak"
  cls[df$strong] <- "DC2 strong"
  cls[df$hbfss_only] <- "HBFSS only"
  cls[df$overlap] <- "Overlap"
  factor(cls, levels = c("Neither", "DC2 weak", "DC2 strong", "HBFSS only", "Overlap"))
}

plot_volc_one <- function(df, title_txt) {
  df$cls <- make_volc_class(df)
  use <- !is.na(df$emp_p) & is.finite(df$emp_p)
  dd <- df[use, , drop = FALSE]
  ggplot(dd, aes(lfc, -log10(pmax(emp_p, 1e-300)))) +
    geom_point(aes(color = cls, shape = cls), alpha = 0.8, size = 1.4) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = pal$thr) +
    geom_hline(yintercept = if (is.finite(dd$hc[1]) && !is.na(dd$hc[1])) -log10(dd$hc[1]) else NA_real_, linetype = "dashed", colour = pal$thr) +
    scale_color_manual(values = c("Neither" = pal$grey, "DC2 weak" = pal$weak, "DC2 strong" = pal$strong, "HBFSS only" = pal$hbfss, "Overlap" = pal$overlap)) +
    scale_shape_manual(values = c("Neither" = 16, "DC2 weak" = 15, "DC2 strong" = 17, "HBFSS only" = 18, "Overlap" = 25)) +
    theme_seq() +
    labs(title = title_txt, x = "Shrunken log2 fold change", y = expression(-log[10](p[empirical])), color = NULL, shape = NULL)
}

plot_volc_panel <- function(cmp, raw_df, lead_df, rem_df, path) {
  p1 <- plot_volc_one(raw_df, paste0(cmp, " | Raw"))
  p2 <- plot_volc_one(lead_df, paste0(cmp, " | Lead"))
  p3 <- plot_volc_one(rem_df, paste0(cmp, " | Rem"))
  g <- arrangeGrob(p1 + theme(legend.position = "none"), p2 + theme(legend.position = "none"), p3,
                   ncol = 3, top = textGrob(paste0(cmp, " | volcano panel"), gp = gpar(fontface = "bold", cex = 1.1)))
  save_grob(g, path, 16, 6)
}

# -----------------------------------------------------------------------------
# summaries
# -----------------------------------------------------------------------------

count_calls <- function(df) {
  c(std = sum(df$std, na.rm = TRUE), weak = sum(df$weak, na.rm = TRUE), strong = sum(df$strong, na.rm = TRUE), hbfss = sum(df$hbfss_sig, na.rm = TRUE), overlap = sum(df$overlap, na.rm = TRUE))
}

# -----------------------------------------------------------------------------
# comparison runner
# -----------------------------------------------------------------------------

run_cmp <- function(cmp_row) {
  cmp <- cmp_row$cmp
  cmp_dir <- file.path(out_dir, cmp)
  fig_dir <- file.path(cmp_dir, "fig")
  tab_dir <- file.path(cmp_dir, "tab")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

  log_step(cmp, " loaded")

  trt_ids <- get_ids(cmp_row$trt)
  ctrl_ids <- get_ids(cmp_row$ctrl)
  ids <- c(trt_ids, ctrl_ids)
  cnt <- wt[, ids, drop = FALSE]
  rownames(cnt) <- wt$OrigID
  cnt <- as.matrix(cnt)
  mode(cnt) <- "numeric"
  meta <- meta_all[ids, , drop = FALSE]

  rank_tbl <- build_rank_tbl(cnt, trt_ids, ctrl_ids)
  log_step(cmp, " rank built")

  raw_res <- run_ds(cnt, meta, paste0(cmp, "_raw"))
  log_step(cmp, " raw done")

  rank_df <- build_rank_nb(rank_tbl, raw_res)
  ref <- build_ref_profiles(rank_df)
  trt_ref <- ref$trt$cut
  ctrl_ref <- ref$ctrl$cut
  comb_ref <- ref$comb$cut
  log_step(cmp, " refs | trt=", trt_ref, " | ctrl=", ctrl_ref, " | comb=", comb_ref)

  scan <- scan_wald(rank_df, min(trt_ref, ctrl_ref, comb_ref, na.rm = TRUE), max(trt_ref, ctrl_ref, comb_ref, na.rm = TRUE))
  final_rank <- as.integer(scan$sel$rank[1])
  log_step(cmp, " final rank=", final_rank, " | method=raw_wald")

  split <- build_split(cnt, rank_tbl, final_rank)
  log_step(cmp, " split built | lead=", nrow(split$lead), " | rem=", nrow(split$rem))

  lead_res <- run_ds(split$lead, meta, paste0(cmp, "_lead"))
  log_step(cmp, " lead done")
  rem_res <- run_ds(split$rem, meta, paste0(cmp, "_rem"))
  log_step(cmp, " rem done")

  plot_cut_panel(cmp, ref, comb_ref, file.path(fig_dir, "cut_panel.png"))
  log_step(cmp, " cut panel saved")
  plot_wald_panel(cmp, rank_df, scan, file.path(fig_dir, "wald_panel.png"))
  log_step(cmp, " wald panel saved")
  plot_volc_panel(cmp, raw_res, lead_res, rem_res, file.path(fig_dir, "volc_panel.png"))
  log_step(cmp, " volc panel saved")

  raw_out <- raw_res[, c("feature_id", "gene", "baseMean", "lfc", "stat", "pvalue", "padj", "emp_p", "emp_q", "lfdr", "hc", "hbfss_thr", "hbfss", "std", "weak", "strong", "hbfss_sig", "hbfss_only", "overlap", "dir")]
  lead_out <- lead_res[, c("feature_id", "gene", "baseMean", "lfc", "stat", "pvalue", "padj", "emp_p", "emp_q", "lfdr", "hc", "hbfss_thr", "hbfss", "std", "weak", "strong", "hbfss_sig", "hbfss_only", "overlap", "dir")]
  rem_out <- rem_res[, c("feature_id", "gene", "baseMean", "lfc", "stat", "pvalue", "padj", "emp_p", "emp_q", "lfdr", "hc", "hbfss_thr", "hbfss", "std", "weak", "strong", "hbfss_sig", "hbfss_only", "overlap", "dir")]
  write.csv(raw_out, file.path(tab_dir, "raw.csv"), row.names = FALSE)
  write.csv(lead_out, file.path(tab_dir, "lead.csv"), row.names = FALSE)
  write.csv(rem_out, file.path(tab_dir, "rem.csv"), row.names = FALSE)

  cut_out <- data.frame(cmp = cmp, method = "wave_raw_wald", trt_ref = trt_ref, ctrl_ref = ctrl_ref, comb_ref = comb_ref, rank = final_rank, stringsAsFactors = FALSE)
  write.csv(cut_out, file.path(tab_dir, "cut.csv"), row.names = FALSE)
  write.csv(scan$fine, file.path(tab_dir, "scan.csv"), row.names = FALSE)

  wave_out <- ref$comb$df
  wave_out$iod_fit <- ref$comb$iod_hat
  wave_out$cv2_fit <- ref$comb$cv2_hat
  wave_out$delta_fit <- ref$comb$iod_hat - ref$comb$cv2_hat
  write.csv(wave_out, file.path(tab_dir, "wave.csv"), row.names = FALSE)
  log_step(cmp, " tables saved")

  rc <- count_calls(raw_res); lc <- count_calls(lead_res); mc <- count_calls(rem_res)
  sum_row <- data.frame(
    cmp = cmp,
    method = "wave_raw_wald",
    rank = final_rank,
    n_raw = nrow(raw_res),
    n_lead = nrow(lead_res),
    n_rem = nrow(rem_res),
    raw_std = rc["std"], raw_weak = rc["weak"], raw_strong = rc["strong"], raw_hbfss = rc["hbfss"], raw_overlap = rc["overlap"],
    lead_std = lc["std"], lead_weak = lc["weak"], lead_strong = lc["strong"], lead_hbfss = lc["hbfss"], lead_overlap = lc["overlap"],
    rem_std = mc["std"], rem_weak = mc["weak"], rem_strong = mc["strong"], rem_hbfss = mc["hbfss"], rem_overlap = mc["overlap"],
    stringsAsFactors = FALSE
  )

  list(sum = sum_row)
}

# -----------------------------------------------------------------------------
# run all
# -----------------------------------------------------------------------------

all_sum <- list()
failed <- list()
for (i in seq_len(nrow(cmp_tbl))) {
  row <- cmp_tbl[i, , drop = FALSE]
  cmp <- row$cmp[1]
  log_step("Running ", cmp)
  x <- try(run_cmp(row), silent = TRUE)
  if (inherits(x, "try-error")) {
    failed[[length(failed) + 1L]] <- data.frame(cmp = cmp, err = as.character(attr(x, "condition")$message), stringsAsFactors = FALSE)
    next
  }
  all_sum[[length(all_sum) + 1L]] <- x$sum
}

if (length(all_sum)) write.csv(bind_rows(all_sum), file.path(out_dir, "cmp_sum.csv"), row.names = FALSE)
if (length(failed)) write.csv(bind_rows(failed), file.path(out_dir, "failed.csv"), row.names = FALSE)

log_step("done")
