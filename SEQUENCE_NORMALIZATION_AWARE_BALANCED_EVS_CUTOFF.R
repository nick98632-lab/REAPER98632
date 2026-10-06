# =============================================================================
# HC10 REVISION — 2026-08-23
# Higher Criticism is restricted to the lowest 10% of ordered empirical-null
# p-values (HC_ALPHA0 = 0.10) and requires HCmax > 0. This revision writes to
# a separate output tree so it cannot overwrite the prior empirical-EVS run.
# =============================================================================

#!/usr/bin/env Rscript

PIPELINE_BUILD <- "SEQUENCE_LIKELIHOOD_CUTOFF_2026-10-04"

# =============================================================================
# SEQUENCE MANUSCRIPT ANALYSIS
# WTTS-Seq PAS analysis with eigenvector splitting, DESeq2, apeglm,
# empirical-null calibration, higher criticism, HBFSS, and 3'aTWAS overlap
# =============================================================================
#
# ANALYSIS UNIT
# Each OrigID is analyzed as an individual polyadenylation-site (PAS) feature.
# Gene symbols are retained as annotation and are not used to collapse PASs
# before differential-expression testing.
#
# COMPARISONS
#   RT0 vs ZT6
#   RT2 vs ZT8
#   RT4 vs ZT10
#   RT8 vs ZT14
#
# ANALYSIS VIEWS PER COMPARISON
#   Original (No EVS)
#   NormEVS Lead
#   NormEVS Rem
#   RawEVS Lead
#   RawEVS Rem
#
# EIGENVECTOR SPLITTING
# NormEVS uses DESeq2 median-of-ratios normalized counts before PCA. RawEVS
# uses raw counts before PCA. Within each comparison, the EVS selection size k*
# is either the prespecified fixed 5,000 comparator or a comparison-specific
# empirical cutoff computed internally from the raw count matrix before any
# downstream DE testing. Empirical cutoffs are generated independently for
# log1p(CPM)-EVS and DESeq2 VST-EVS using the shared-regime, count-based Pareto
# endpoint-chord-maximum algorithm defined below. The same comparison-specific
# active k* is applied to NormEVS and
# RawEVS so those tracks differ only in the matrix used for PC1 ranking, not in
# the number of PASs admitted per condition. Within each condition, prcomp is
# applied directly to the corresponding feature-by-sample matrix, absolute PC1
# feature loadings are ranked, and the k* highest-loading PASs are selected. PASs
# present in both condition-specific top-k* sets are Joint; PASs present in only
# one set are Disjoint. Joint plus both Disjoint sets form the Leading Edge. All
# other PASs form the Remainder. Downstream DESeq2 always receives raw counts for
# the selected PAS subset. Size factors are estimated once per comparison on all
# 10 RT/ZT samples and passed to every view's DESeq2 fit (Original, Lead, Rem),
# so all views share one normalization; dispersions are re-estimated per view.
#
# DIFFERENTIAL EXPRESSION AND EFFECT TESTS
# DESeq2 uses design ~ condition with trt relative to untrt. The ordinary Wald
# p-value is adjusted by Benjamini-Hochberg at FDR 10%. The manuscript Standard
# effect is the ordinary DESeq2 BH-significant result restricted to PASs with
# |apeglm-shrunken LFC| >= 1, so Standard markers cannot occur inside the stated
# effect boundary. Strong (decoupled) calls use DESeq2 greaterAbs with
# lfcThreshold = 1, BH padj < 0.10, and the same |apeglm LFC| >= 1 reporting
# boundary. Weak-CNH support uses DESeq2 lessAbs with lfcThreshold = 1 and BH
# padj < 0.20; the final Weak category additionally requires |apeglm LFC| < 1
# and HBFSS significance. apeglm-shrunken LFC is the reported effect estimate,
# HBFSS effect term, and x-coordinate in all significance figures.
#
# EMPIRICAL NULL, HIGHER CRITICISM, AND HBFSS
# Finite DESeq2 Wald statistics are calibrated with fdrtool using a normal
# empirical-null model. HBFSS uses the resulting empirical p-values. Empirical
# Higher-Criticism scores are calculated from the sorted empirical p-values.
# The HC maximization is restricted a priori to the lowest 10% of ordered
# empirical p-values (alpha0 = 0.10), preventing the threshold from being chosen
# from the uninformative p~1 boundary. A dataset/view receives an HC threshold
# only when the maximum HC score within that search region is strictly positive;
# otherwise HCp and Htau are undefined and HBFSS significance is disabled for
# that dataset/view. With the manuscript LFC boundary c = 1:
#
#   Htau = -log10(HCp) * c
#   HBFSS = |apeglm LFC| * [-log10(empirical p)]
#
# When a valid positive-HC threshold exists, a PAS is HBFSS-significant when
# HBFSS > Htau. HCp is used to derive Htau and is not imposed as an additional
# significance gate.
#
# VOLCANO FIGURES
# Volcano x-axis: apeglm-shrunken log2 fold change.
# Volcano y-axis: -log10(empirical p) from the fdrtool empirical-null model.
# Standard, Strong, Weak, and HBFSS are plotted as separate method layers using
# one fixed color/marker key across every significance figure. Weak markers are
# shown only for final Weak discoveries (lessAbs + HBFSS), never for lessAbs-only
# PASs. Method overlap is reported numerically and by superimposed method markers;
# it is not treated as a fifth significance method. Each volcano labels at most
# the top 20 final significant PASs.
#
# 3'aTWAS COMPARISON
# After all WTTS analyses and manuscript figures are complete, human 3'aTWAS
# gene symbols are mapped to rat orthologs. TWAS ortholog overlap is evaluated
# against every final WTTS significance method (Standard, Strong, Weak, HBFSS)
# in every analysis view. The TWAS exports identify the human TWAS symbol, rat
# ortholog, significant PASs, method support, and whether the rat gene contains
# multiple WTTS PAS features consistent with alternative polyadenylation.
# =============================================================================



# =============================================================================
# INTERNAL EMPIRICAL CUTOFF ENGINE
# =============================================================================
# This is deliberately run by the MASTER process before any downstream
# differential-expression analysis.  No empirical k* values are hard-coded.
#
# For each preprocessing branch (log1p(CPM)-EVS and DESeq2 VST-EVS):
#   1) rank PASs within each of the 8 arms by ascending |PC1 loading|;
#   2) construct P_i = lambda_1 * v_i1^2 and the globally normalized
#      excess-over-Poisson reference E_i = max(V_pool,i - mu_i,g, 0);
#   3) fit one shared c1/c2 two-knot cumulative-divergence model across all
#      8 arm-specific D_g(r)=F_E,g(r)-F_P,g(r) curves;
#   4) for each RT/ZT comparison and candidate top-k, count
#         G(k) = Joint + Disjoint whose opposite-arm rank is LE or Divergence
#         R(k) = Disjoint whose opposite-arm rank is Remainder (< c1);
#   5) retain the nondominated Pareto frontier, min-max normalize G and R on
#      that frontier, connect its two endpoints with a chord, and select
#         k* = argmax perpendicular distance to the endpoint chord.
#
# The resulting k* values are written to Cutoff_Method_Manifest.csv and are
# then read back by each child run.  The fixed 5,000 comparator remains fully
# independent of empirical cutoff estimation.
# =============================================================================

SEQCUT_GROUP_PATTERNS <- c(
  RT0="^R0_", ZT6="^ZT6_", RT2="^R2_", ZT8="^ZT8_",
  RT4="^R4_", ZT10="^ZT10_", RT8="^R8_", ZT14="^ZT14_"
)

SEQCUT_COMPARISONS <- list(
  RT0_ZT6=c(control="RT0", treatment="ZT6"),
  RT2_ZT8=c(control="RT2", treatment="ZT8"),
  RT4_ZT10=c(control="RT4", treatment="ZT10"),
  RT8_ZT14=c(control="RT8", treatment="ZT14")
)

seqcut_read_counts <- function(path) {
  if (!file.exists(path)) stop("Count file not found for empirical cutoff derivation: ", path)
  x <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  idx <- sort(unique(unlist(lapply(SEQCUT_GROUP_PATTERNS, function(p) grep(p, names(x))))))
  if (!length(idx)) stop("No sample columns matched SEQCUT_GROUP_PATTERNS.")

  ids <- trimws(as.character(x[[1]]))
  blank <- is.na(ids) | ids == ""
  ids[blank] <- paste0("feature_", which(blank))
  ids <- make.unique(ids, sep="__dup_")

  m <- do.call(cbind, lapply(x[,idx,drop=FALSE], function(z)
    suppressWarnings(as.numeric(as.character(z)))
  ))
  rownames(m) <- ids
  colnames(m) <- names(x)[idx]
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  m <- pmax(m, 0)
  m <- m[rowSums(m) > 0, , drop=FALSE]
  if (nrow(m) < 20L) stop("Too few nonzero PASs for empirical cutoff derivation.")
  m
}

seqcut_assign_groups <- function(samples) {
  g <- rep(NA_character_, length(samples)); names(g) <- samples
  for (nm in names(SEQCUT_GROUP_PATTERNS)) {
    idx <- grep(SEQCUT_GROUP_PATTERNS[[nm]], samples)
    if (any(!is.na(g[idx]))) stop("A sample matched more than one empirical-cutoff group.")
    g[idx] <- nm
  }
  if (anyNA(g)) stop("Unassigned cutoff samples: ", paste(names(g)[is.na(g)], collapse=", "))
  factor(g, levels=names(SEQCUT_GROUP_PATTERNS))
}

seqcut_normalize_and_vst <- function(counts, groups) {
  if (!requireNamespace("DESeq2", quietly=TRUE)) stop("DESeq2 is required to derive empirical cutoffs.")
  if (!requireNamespace("SummarizedExperiment", quietly=TRUE)) stop("SummarizedExperiment is required for VST empirical cutoffs.")

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=round(counts),
    colData=data.frame(group=groups, row.names=colnames(counts)),
    design=~group
  )
  dds <- tryCatch(
    DESeq2::estimateSizeFactors(dds),
    error=function(e) DESeq2::estimateSizeFactors(dds, type="poscounts")
  )
  vst_obj <- DESeq2::varianceStabilizingTransformation(dds, blind=TRUE)
  vst_mat <- SummarizedExperiment::assay(vst_obj)
  norm_mat <- DESeq2::counts(dds, normalized=TRUE)

  if (!identical(dim(vst_mat), dim(counts)) ||
      !identical(rownames(vst_mat), rownames(counts)) ||
      !identical(colnames(vst_mat), colnames(counts))) {
    stop("VST matrix is not aligned with the raw count matrix.")
  }
  if (any(!is.finite(vst_mat))) stop("VST matrix contains non-finite values.")
  list(norm=norm_mat, vst=vst_mat)
}

seqcut_cpm_log1p <- function(x) {
  lib <- colSums(x)
  lib[!is.finite(lib) | lib <= 0] <- 1
  log1p(sweep(x, 2, lib/1e6, "/"))
}

seqcut_pooled_var <- function(norm, groups) {
  sse <- rep(0, nrow(norm)); df <- 0L
  for (g in levels(groups)) {
    z <- norm[,groups==g,drop=FALSE]
    mu <- rowMeans(z)
    sse <- sse + rowSums((z-mu)^2)
    df <- df + ncol(z)-1L
  }
  if (df < 2L) stop("Insufficient pooled residual df for cutoff derivation.")
  pmax(sse/df, 0)
}

seqcut_pc1_rank <- function(x) {
  p <- stats::prcomp(t(x), center=TRUE, scale.=FALSE, rank.=1)
  loading <- p$rotation[,1]
  loading[!is.finite(loading)] <- 0
  list(
    order=order(abs(loading), decreasing=FALSE),
    loading=loading,
    contribution=(p$sdev[1]^2)*loading^2
  )
}

seqcut_build_arm <- function(method, raw, norm, vst, pooled_var) {
  rank_matrix <- switch(
    method,
    CPM_EVS=seqcut_cpm_log1p(raw),
    VST_EVS=vst,
    stop("Unknown empirical cutoff method: ", method)
  )
  pc <- seqcut_pc1_rank(rank_matrix)
  ord <- pc$order
  P <- pc$contribution[ord]
  mu <- rowMeans(norm)[ord]
  E <- pmax(pooled_var[ord]-mu, 0)
  if (sum(P) <= 0 || sum(E) <= 0) stop("Undefined PC1/NB mass during empirical cutoff derivation.")
  p <- P/sum(P); q <- E/sum(E)
  data.frame(
    rank=seq_along(ord),
    feature_id=rownames(raw)[ord],
    D=cumsum(q)-cumsum(p),
    stringsAsFactors=FALSE
  )
}

seqcut_piecewise_basis <- function(x,c1,c2) {
  cbind(1, x, pmax(x-c1,0), pmax(x-c2,0))
}

seqcut_fit_shared_knots <- function(arms) {
  N <- nrow(arms[[1]])
  if (length(unique(vapply(arms,nrow,integer(1)))) != 1L) stop("Cutoff arms do not share one PAS universe.")
  x_full <- (seq_len(N)-1)/(N-1)
  D_full <- do.call(cbind,lapply(arms,function(z) z$D))

  opt_n <- min(4000L,N)
  idx <- unique(as.integer(round(seq(1,N,length.out=opt_n))))
  x <- x_full[idx]
  D <- D_full[idx,,drop=FALSE]
  min_gap <- max(4/(N-1),1e-4)

  objective <- function(par) {
    c1 <- par[1]; c2 <- par[2]
    if (!is.finite(c1) || !is.finite(c2) || c1<=0 || c2>=1 || c2-c1<=min_gap) return(1e100)
    X <- seqcut_piecewise_basis(x,c1,c2)
    b <- tryCatch(qr.coef(qr(X),D), error=function(e) NULL)
    if (is.null(b) || any(!is.finite(b))) return(1e100)
    sum((D-X%*%b)^2)
  }

  starts <- list(c(.05,.95),c(.12,.88),c(.22,.78),c(.32,.68))
  fits <- lapply(starts,function(st) stats::optim(st,objective,method="Nelder-Mead",control=list(maxit=500,reltol=1e-10)))
  best <- fits[[which.min(vapply(fits,function(z) z$value,numeric(1)))]]

  c1 <- as.integer(round(1+best$par[1]*(N-1)))
  c2 <- as.integer(round(1+best$par[2]*(N-1)))
  c1 <- max(2L,min(N-2L,c1)); c2 <- max(c1+1L,min(N-1L,c2))

  X_full <- seqcut_piecewise_basis(x_full,(c1-1)/(N-1),(c2-1)/(N-1))
  coef <- qr.coef(qr(X_full),D_full)
  fitted <- X_full%*%coef
  list(c1=c1,c2=c2,N=N,SSE=sum((D_full-fitted)^2))
}

seqcut_rank_map <- function(df) setNames(df$rank,df$feature_id)

seqcut_add_events <- function(pos,weight,L) {
  out <- numeric(L)
  ok <- is.finite(pos) & is.finite(weight) & pos>=1 & pos<=L
  if (!any(ok)) return(out)
  z <- rowsum(weight[ok], group=as.integer(pos[ok]), reorder=FALSE)
  out[as.integer(rownames(z))] <- z[,1]
  out
}

seqcut_activate_from <- function(depth,weight,K) cumsum(seqcut_add_events(depth,weight,K))

seqcut_active_interval <- function(start,end,weight,K) {
  delta <- numeric(K+1L)
  ok <- is.finite(start) & is.finite(end) & is.finite(weight) & start>=1 & start<=K & end>start
  if (!any(ok)) return(numeric(K))
  ss <- as.integer(start[ok]); ee <- pmin(as.integer(end[ok]),K+1L); w <- weight[ok]
  delta <- delta + seqcut_add_events(ss,w,K+1L)
  delta <- delta - seqcut_add_events(ee,w,K+1L)
  cumsum(delta)[seq_len(K)]
}

seqcut_scan_pair <- function(control_df,treatment_df,c1,c2) {
  N <- nrow(control_df); K <- N-c2
  if (K < 1L) stop("No candidate k exists beyond c2.")

  rC_map <- seqcut_rank_map(control_df); rT_map <- seqcut_rank_map(treatment_df)
  ids <- control_df$feature_id
  rC <- as.integer(rC_map[ids]); rT <- as.integer(rT_map[ids])
  if (anyNA(rC) || anyNA(rT)) stop("Cross-arm PAS rank mapping failed.")
  dC <- N-rC+1L; dT <- N-rT+1L
  joint_n <- seqcut_activate_from(pmax(dC,dT),rep(1,N),K)

  idxC <- which(dC<dT & dC<=K)
  startC <- dC[idxC]; endC <- pmin(dT[idxC],K+1L); oppC <- rT[idxC]
  disC_rem <- seqcut_active_interval(startC,endC,as.numeric(oppC<c1),K)
  disC_div <- seqcut_active_interval(startC,endC,as.numeric(oppC>=c1 & oppC<=c2),K)
  disC_le  <- seqcut_active_interval(startC,endC,as.numeric(oppC>c2),K)

  idxT <- which(dT<dC & dT<=K)
  startT <- dT[idxT]; endT <- pmin(dC[idxT],K+1L); oppT <- rC[idxT]
  disT_rem <- seqcut_active_interval(startT,endT,as.numeric(oppT<c1),K)
  disT_div <- seqcut_active_interval(startT,endT,as.numeric(oppT>=c1 & oppT<=c2),K)
  disT_le  <- seqcut_active_interval(startT,endT,as.numeric(oppT>c2),K)

  opposite_le_n <- disC_le+disT_le
  opposite_divergence_n <- disC_div+disT_div
  remainder_cross_n <- disC_rem+disT_rem
  good_n <- joint_n+opposite_le_n+opposite_divergence_n
  union_n <- 2*seq_len(K)-joint_n
  if (!all(good_n+remainder_cross_n == union_n)) stop("Internal union-count mismatch in cutoff scan.")

  data.frame(
    k=seq_len(K), joint_n=joint_n, opposite_le_n=opposite_le_n,
    opposite_divergence_n=opposite_divergence_n,
    remainder_cross_n=remainder_cross_n, good_n=good_n, union_n=union_n,
    stringsAsFactors=FALSE
  )
}

seqcut_pareto_endpoint_max <- function(scan_df) {
  tmp <- data.frame(row_id=seq_len(nrow(scan_df)), k=scan_df$k,
                    good=scan_df$good_n, cost=scan_df$remainder_cross_n)
  tmp <- tmp[order(tmp$cost,-tmp$good,-tmp$k),,drop=FALSE]
  tmp <- tmp[!duplicated(tmp$cost),,drop=FALSE]
  prior_best <- c(-Inf, head(cummax(tmp$good),-1L))
  frontier <- tmp[tmp$good > prior_best,,drop=FALSE]
  frontier <- frontier[order(frontier$cost,frontier$good,frontier$k),,drop=FALSE]
  if (nrow(frontier) < 2L) stop("At least two Pareto points are required for the endpoint chord.")

  gr <- range(frontier$good); cr <- range(frontier$cost)
  norm_good <- function(x) if (diff(gr)==0) rep(0,length(x)) else (x-gr[1])/diff(gr)
  norm_cost <- function(x) if (diff(cr)==0) rep(0,length(x)) else (x-cr[1])/diff(cr)
  frontier$good_norm <- norm_good(frontier$good)
  frontier$remainder_norm <- norm_cost(frontier$cost)

  x1 <- frontier$remainder_norm[1]; y1 <- frontier$good_norm[1]
  x2 <- frontier$remainder_norm[nrow(frontier)]; y2 <- frontier$good_norm[nrow(frontier)]
  denom <- sqrt((y2-y1)^2+(x2-x1)^2)
  if (!is.finite(denom) || denom <= .Machine$double.eps) stop("Endpoint chord is undefined.")
  frontier$endpoint_deviation <- abs(
    (y2-y1)*frontier$remainder_norm - (x2-x1)*frontier$good_norm + x2*y1-y2*x1
  )/denom

  best <- max(frontier$endpoint_deviation,na.rm=TRUE)
  cand <- frontier[abs(frontier$endpoint_deviation-best)<1e-12,,drop=FALSE]
  cand <- cand[order(-cand$good,cand$cost,-cand$k),,drop=FALSE]
  list(selected_k=as.integer(cand$k[1]), selected=cand[1,,drop=FALSE], frontier=frontier)
}

derive_sequence_empirical_cutoffs <- function(count_path) {
  message("Deriving CPM-EVS and VST-EVS empirical cutoffs from the raw count matrix...")
  counts <- seqcut_read_counts(count_path)
  groups <- seqcut_assign_groups(colnames(counts))
  norm_obj <- seqcut_normalize_and_vst(counts,groups)
  pooled_var <- seqcut_pooled_var(norm_obj$norm,groups)

  method_results <- list(); audit_rows <- list(); details <- list()
  for (method in c("CPM_EVS","VST_EVS")) {
    arms <- list()
    for (g in levels(groups)) {
      idx <- which(groups==g)
      arms[[g]] <- seqcut_build_arm(
        method=method,
        raw=counts[,idx,drop=FALSE],
        norm=norm_obj$norm[,idx,drop=FALSE],
        vst=norm_obj$vst[,idx,drop=FALSE],
        pooled_var=pooled_var
      )
    }
    knot <- seqcut_fit_shared_knots(arms)
    message("  ",method,": shared c1=",knot$c1," | c2=",knot$c2," | SSE=",signif(knot$SSE,7))
    ks <- integer(length(SEQCUT_COMPARISONS)); names(ks) <- names(SEQCUT_COMPARISONS)
    det <- list(knot = knot, arms = arms, scans = list(), frontiers = list())

    for (nm in names(SEQCUT_COMPARISONS)) {
      mp <- SEQCUT_COMPARISONS[[nm]]
      scan <- seqcut_scan_pair(arms[[mp[["control"]]]],arms[[mp[["treatment"]]]],knot$c1,knot$c2)
      opt <- seqcut_pareto_endpoint_max(scan)
      kstar <- opt$selected_k
      det$scans[[nm]] <- scan
      det$frontiers[[nm]] <- opt$frontier
      sel <- scan[scan$k==kstar,,drop=FALSE]
      ks[nm] <- kstar
      audit_rows[[paste(method,nm,sep="__")]] <- data.frame(
        Method=method, Comparison=nm, N=knot$N, c1=knot$c1, c2=knot$c2,
        Candidate_k_max=knot$N-knot$c2, k_star=kstar,
        G=sel$good_n, R=sel$remainder_cross_n,
        Joint=sel$joint_n, Opposite_LE=sel$opposite_le_n,
        Opposite_Divergence=sel$opposite_divergence_n,
        Endpoint_deviation=opt$selected$endpoint_deviation,
        Knot_SSE=knot$SSE, stringsAsFactors=FALSE
      )
      message("    ",nm,": k*=",kstar," | G=",sel$good_n," | R=",sel$remainder_cross_n,
              " | endpoint deviation=",round(opt$selected$endpoint_deviation,6))
    }
    method_results[[method]] <- ks
    details[[method]] <- det
  }

  audit <- do.call(rbind,audit_rows); rownames(audit) <- NULL
  list(cpm=method_results[["CPM_EVS"]], vst=method_results[["VST_EVS"]], audit=audit, details=details)
}


# =============================================================================
# NORMALIZED-COUNT REGIME KNOTS (used only to label the actual NormEVS split)
# =============================================================================
# Same shared two-knot fit as the CPM/VST engine, on median-of-ratios
# normalized counts (the NormEVS representation). No cutoff is chosen here.

seqcut_build_arm_norm <- function(norm, pooled_var) {
  pc <- seqcut_pc1_rank(norm)
  ord <- pc$order
  P <- pc$contribution[ord]
  mu <- rowMeans(norm)[ord]
  E <- pmax(pooled_var[ord] - mu, 0)
  if (sum(P) <= 0 || sum(E) <= 0) stop("Undefined PC1 or excess-variance mass on normalized counts.")
  data.frame(rank = seq_along(ord), feature_id = rownames(norm)[ord],
             D = cumsum(E / sum(E)) - cumsum(P / sum(P)), stringsAsFactors = FALSE)
}

derive_normcount_knots <- function(count_path) {
  message("Fitting shared regime knots on normalized-count PC1 rankings...")
  counts <- seqcut_read_counts(count_path)
  groups <- seqcut_assign_groups(colnames(counts))
  norm_obj <- seqcut_normalize_and_vst(counts, groups)
  pooled_var <- seqcut_pooled_var(norm_obj$norm, groups)
  arms <- list()
  for (g in levels(groups)) {
    arms[[g]] <- seqcut_build_arm_norm(norm_obj$norm[, groups == g, drop = FALSE], pooled_var)
  }
  knot <- seqcut_fit_shared_knots(arms)
  message("  NormEVS knots: c1=", knot$c1, " | c2=", knot$c2, " | N=", knot$N)
  list(c1 = knot$c1, c2 = knot$c2, N = knot$N)
}

# =============================================================================
# CUTOFF METHODS, SHARED FIGURE STYLE AND THE FIXED-vs-LIKELIHOOD COMPARISON
# =============================================================================
# Two cutoff methods are run:
#   Fixed 5,000    prespecified comparator (as in the manuscript)
#   Likelihood k*  per comparison, the k that maximizes the likelihood of the
#                  Leading Edge / Remainder dispersion model (two DESeq2 trends,
#                  two prior widths) on the actual NormEVS split. Computed inside
#                  its own child run before any testing; see likelihood_cutoff_scan().
#
# Why the split helps DESeq2: DESeq2 shrinks every gene-wise dispersion toward
# one mean-dispersion trend with one prior width for all PASs. When PASs at the
# same mean differ in dispersion, that prior is centred wrongly and too wide.
# EVS separates the populations, so each subset gets its own trend and prior.

METHOD_LEVELS <- c("Fixed 5,000", "Likelihood k*")
METHOD_SLUGS  <- c("Fixed 5,000" = "Fixed_5000", "Likelihood k*" = "Likelihood_Width")
METHOD_COLS   <- c("Fixed 5,000" = "#8C8C8C", "Likelihood k*" = "#2E7D5B")
SPLIT_COLS <- c("Original (no EVS)" = "grey15", "Leading Edge" = "#1F6FB2", "Remainder" = "#C0392B")
SPLIT_POINT_COLS <- c("Leading Edge" = "#5B9BD5", "Remainder" = "#E8836F")
SPLIT_TREND_COLS <- c("Leading Edge" = "#0B3C6E", "Remainder" = "#8E1B10", "Original (no EVS)" = "black")
SPLIT_VIEW_SET <- c("Original (No EVS)" = "Original (no EVS)", "NormEVS Lead" = "Leading Edge",
                    "NormEVS Rem" = "Remainder")

# Normalized-count regime knots, used only to label the actual Leading Edge.
write_cutoff_data <- function(norm_knots, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data.frame(c1 = norm_knots$c1, c2 = norm_knots$c2, N = norm_knots$N),
                   file.path(out_dir, "NormEVS_Regime_Knots.csv"), row.names = FALSE)
  invisible(TRUE)
}

collect_method_tables <- function(multi_root, slugs, labels, file) {
  rows <- list()
  for (j in seq_along(slugs)) {
    f <- file.path(multi_root, slugs[j], "Summary_Tables", file)
    if (!file.exists(f)) next
    d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    d$Cutoff_Method <- labels[j]
    rows[[length(rows) + 1L]] <- d
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}

master_theme <- function() {
  theme_classic(base_size = 7.5, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.3, colour = "grey20"),
          axis.ticks = element_line(linewidth = 0.3, colour = "grey20"),
          axis.text = element_text(colour = "grey10", size = 6.5), axis.title = element_text(size = 7.3),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7.5),
          panel.spacing.x = unit(4, "mm"), panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"),
          legend.title = element_blank(), legend.text = element_text(size = 6.8),
          legend.key.size = unit(3, "mm"), legend.position = "none",
          plot.title = element_text(size = 7.8, face = "bold", margin = margin(b = 1)),
          plot.subtitle = element_text(size = 6.8, colour = "grey30", margin = margin(b = 3)),
          plot.title.position = "plot", plot.tag = element_text(size = 10, face = "bold"),
          plot.margin = margin(3, 4, 3, 3))
}

master_save <- function(p, out_dir, stem, w_mm, h_mm) {
  pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = w_mm, height = h_mm, units = "mm", device = pdf_dev)
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w_mm, height = h_mm, units = "mm",
         dpi = 600, device = "png", bg = "white")
}

master_plotting_ready <- function() {
  ok <- requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)
  if (!ok) warning("ggplot2/patchwork not installed: master figure skipped (tables still written).")
  if (ok) suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  ok
}

try_fig <- function(expr, what) {
  tryCatch(expr, error = function(e) {
    warning("Figure skipped (", what, "): ", conditionMessage(e), call. = FALSE)
    invisible(NULL)
  })
}

# Fixed 5,000 vs Likelihood k*: (a) where each cutoff sits on the likelihood
# curve, (b) the resulting dispersion prior widths, (c) discoveries (reference).
save_cutoff_comparison_figure <- function(multi_root, out_dir) {
  if (!master_plotting_ready()) return(invisible(FALSE))
  cmp <- names(SEQCUT_COMPARISONS)
  lab <- labeller(Comparison = setNames(sub("_", " / ", cmp), cmp))
  ml_dir <- file.path(multi_root, METHOD_SLUGS[["Likelihood k*"]], "Summary_Tables")
  scan <- utils::read.csv(file.path(ml_dir, "Table_Likelihood_Cutoff_Scan.csv"), stringsAsFactors = FALSE)
  sel <- utils::read.csv(file.path(ml_dir, "Table_Likelihood_Cutoff_Selected.csv"), stringsAsFactors = FALSE)
  scan$Comparison <- factor(scan$Comparison, levels = cmp); sel$Comparison <- factor(sel$Comparison, levels = cmp)
  marks <- rbind(
    data.frame(Comparison = sel$Comparison, k = sel$k_star, Method = "Likelihood k*"),
    data.frame(Comparison = sel$Comparison, k = 5000L, Method = "Fixed 5,000"))
  marks$y <- scan$delta_loglik_vs_Original[match(paste(marks$Comparison, marks$k), paste(scan$Comparison, scan$k))]
  marks$Method <- factor(marks$Method, levels = METHOD_LEVELS)

  pa <- ggplot(scan, aes(k, delta_loglik_vs_Original)) +
    geom_rect(data = sel, aes(xmin = support_lo, xmax = support_hi, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = METHOD_COLS[["Likelihood k*"]], alpha = 0.12) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.35, colour = "grey30") +
    geom_line(linewidth = 0.6, colour = "grey15") +
    geom_point(data = marks, aes(y = y, colour = Method), size = 2.4) +
    facet_wrap(~ Comparison, nrow = 1, scales = "free_y", labeller = lab) +
    scale_colour_manual(values = METHOD_COLS) +
    scale_x_continuous(labels = scales::label_comma(), breaks = scales::breaks_pretty(3)) +
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0.08, 0.12))) +
    labs(x = "k (top-|PC1| PASs taken from each condition)", y = "Log-likelihood gain over\nthe unsplit DESeq2 model",
         title = "Where each cutoff sits on the likelihood of the Leading Edge / Remainder dispersion model",
         subtitle = "Shaded: likelihood support interval around k*. Above 0 = the split model fits better than DESeq2's single prior.") +
    master_theme() + theme(legend.position = "top", legend.justification = "right")

  sc <- collect_method_tables(multi_root, unname(METHOD_SLUGS), METHOD_LEVELS, "Table_Dispersion_Split_Score.csv")
  pb <- NULL
  if (!is.null(sc) && nrow(sc)) {
    o <- unique(sc[, c("Comparison", "w_Original")]); names(o)[2] <- "W"; o$Series <- "Original (no EVS)"
    b <- rbind(o, data.frame(Comparison = sc$Comparison, W = sc$W_split, Series = sc$Cutoff_Method))
    b$Series <- factor(b$Series, levels = c("Original (no EVS)", METHOD_LEVELS))
    b$Comparison <- factor(b$Comparison, levels = cmp)
    pb <- ggplot(b, aes(Comparison, W, fill = Series)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.74) +
      geom_text(aes(label = sprintf("%.2f", W)), position = position_dodge(width = 0.8), vjust = -0.4, size = 2.0) +
      scale_fill_manual(values = c("Original (no EVS)" = "grey15", METHOD_COLS)) +
      scale_x_discrete(labels = setNames(sub("_", " / ", cmp), cmp)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = NULL, y = "Pooled dispersion prior\nwidth W (LE + REM)",
           title = "Dispersion prior width: unsplit DESeq2 vs EVS at each cutoff",
           subtitle = "Lower = narrower prior, stronger and better-centred shrinkage of each PAS's dispersion.") +
      master_theme() + theme(legend.position = "top", legend.justification = "right")
  }

  disc <- collect_cross_method_discoveries(multi_root, unname(METHOD_SLUGS), METHOD_LEVELS)
  pc <- NULL
  if (!is.null(disc) && nrow(disc)) {
    d <- disc[disc$Track == "NormEVS" & disc$Set %in% c("Leading Edge", "Remainder"), ]
    d <- aggregate(Any_significant ~ Cutoff_Method + Comparison + Set, data = d, FUN = sum)
    d$Cutoff_Method <- factor(d$Cutoff_Method, levels = METHOD_LEVELS)
    d$Comparison <- factor(d$Comparison, levels = cmp)
    d$Set <- factor(d$Set, levels = c("Remainder", "Leading Edge"))
    tot <- aggregate(Any_significant ~ Cutoff_Method + Comparison, data = d, FUN = sum)
    o <- aggregate(Any_significant ~ Comparison, data = disc[disc$Set == "Original", ], FUN = max)
    o$Comparison <- factor(o$Comparison, levels = cmp)
    pc <- ggplot(d, aes(Cutoff_Method, Any_significant, fill = Set)) +
      geom_col(width = 0.62) +
      geom_text(data = tot, aes(Cutoff_Method, Any_significant, label = scales::comma(Any_significant)),
                inherit.aes = FALSE, vjust = -0.4, size = 2.0) +
      geom_hline(data = o, aes(yintercept = Any_significant), linetype = "22", linewidth = 0.35) +
      facet_wrap(~ Comparison, nrow = 1, scales = "free_y", labeller = lab) +
      scale_fill_manual(values = c("Leading Edge" = unname(SPLIT_COLS["Leading Edge"]),
                                   "Remainder" = unname(SPLIT_COLS["Remainder"])),
                        breaks = c("Leading Edge", "Remainder")) +
      scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.18))) +
      labs(x = NULL, y = "Significant PASs\n(NormEVS LE + REM)",
           title = "Significant PASs (reference only; not used to choose k)",
           subtitle = "Dashed line: Original (no EVS).") +
      master_theme() + theme(legend.position = "top", legend.justification = "right")
  }
  parts <- list(pa, pb, pc); keep <- !vapply(parts, is.null, logical(1))
  fig <- wrap_plots(parts[keep], ncol = 1) + plot_layout(heights = c(1, 0.9, 1)[keep]) +
    plot_annotation(tag_levels = "a")
  master_save(fig, out_dir, "Figure_Cutoff_Comparison_Fixed_vs_Likelihood", 180, 200)
  invisible(TRUE)
}

# =============================================================================
# CROSS-CUTOFF-METHOD DISCOVERY SUMMARY (written by the master after all runs)
# =============================================================================
# Reads every child run's per-comparison significant-site table and reports,
# for each cutoff method, EVS track and comparison, the number of final
# significant PASs (Std | Strong | Weak | HBFSS) in the Leading Edge and in the
# Remainder, their combined total, and the Original (No EVS) reference.
# Leading Edge and Remainder are disjoint, so their counts add without overlap.

collect_cross_method_discoveries <- function(multi_root, method_slugs, method_labels) {
  views <- c("NormEVS Lead", "NormEVS Rem", "RawEVS Lead", "RawEVS Rem", "Original (No EVS)")
  flag_cols <- c(Std = "Std", Strong = "Strong", Weak = "Weak", HBFSS = "HBFSS_sig")
  rows <- list()
  for (j in seq_along(method_slugs)) {
    for (cn in names(SEQCUT_COMPARISONS)) {
      f <- file.path(multi_root, method_slugs[j], cn, paste0("Table_", cn, "_Significant_Sites_All_Views.csv"))
      sig <- if (file.exists(f) && file.info(f)$size > 0) {
        tryCatch(utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) data.frame())
      } else data.frame()
      for (v in views) {
        s <- if (nrow(sig) && "Analysis" %in% names(sig)) sig[sig$Analysis == v, , drop = FALSE] else sig[0, , drop = FALSE]
        cnt <- vapply(flag_cols, function(col) {
          if (nrow(s) && col %in% names(s)) sum(as.logical(s[[col]]), na.rm = TRUE) else 0L
        }, integer(1))
        rows[[length(rows) + 1L]] <- data.frame(
          Cutoff_Method = method_labels[j], Comparison = cn, Analysis = v,
          Any_significant = nrow(s), Std = cnt[["Std"]], Strong = cnt[["Strong"]],
          Weak = cnt[["Weak"]], HBFSS = cnt[["HBFSS"]], stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$Track <- ifelse(grepl("^NormEVS", out$Analysis), "NormEVS",
                      ifelse(grepl("^RawEVS", out$Analysis), "RawEVS", "Original"))
  out$Set <- ifelse(grepl("Lead$", out$Analysis), "Leading Edge",
                    ifelse(grepl("Rem$", out$Analysis), "Remainder", "Original"))
  out
}


# =============================================================================
# MULTI-CUTOFF ORCHESTRATION
# Runs this same script once per cutoff method (Fixed 5,000 and Likelihood k*),
# then creates one
# final ZIP containing all figures, tables, audit files, and TWAS results.
# =============================================================================

SEQUENCE_CHILD_RUN <- identical(Sys.getenv("SEQUENCE_CHILD_RUN", unset = "0"), "1")

if (!SEQUENCE_CHILD_RUN) {
  script_args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", script_args, value = TRUE)
  if (!length(script_arg)) stop("Run this pipeline with Rscript so the master process can relaunch itself.")
  this_script <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)

  root_candidates <- unique(c(
    Sys.getenv("SEQUENCE_REPO_ROOT", unset = ""),
    getwd(), dirname(getwd()), "/root/REAPER98632", dirname(this_script)
  ))
  root_candidates <- root_candidates[nzchar(root_candidates)]
  repo_guess <- NULL
  for (cand in root_candidates) {
    if (file.exists(file.path(cand, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")) ||
        file.exists(file.path(cand, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"))) {
      repo_guess <- normalizePath(cand, winslash = "/", mustWork = TRUE)
      break
    }
  }
  if (is.null(repo_guess) && dir.exists("/root/REAPER98632")) {
    repo_guess <- normalizePath("/root/REAPER98632", winslash = "/", mustWork = TRUE)
  }
  if (is.null(repo_guess)) repo_guess <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

  multi_root <- Sys.getenv(
    "SEQUENCE_MULTI_ROOT",
    unset = file.path(repo_guess, "exports", "sequence_final_cutoff_methods")
  )
  if (dir.exists(multi_root)) unlink(multi_root, recursive = TRUE, force = TRUE)
  dir.create(multi_root, recursive = TRUE, showWarnings = FALSE)

  count_path <- if (file.exists(file.path(repo_guess, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"))) {
    file.path(repo_guess, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
  } else {
    file.path(repo_guess, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
  }

  norm_knots <- derive_normcount_knots(count_path)
  write_cutoff_data(norm_knots, file.path(multi_root, "Cutoff_Data"))
  cmp_names <- names(SEQCUT_COMPARISONS)
  # Likelihood k* is estimated inside its own child run (it needs the actual
  # NormEVS split and DESeq2 dispersion fits); its row is NA here by design.
  cutoff_manifest <- data.frame(
    Cutoff_Method = c("Fixed_5000", "Likelihood_Width"),
    RT0_ZT6 = c(5000L, NA_integer_), RT2_ZT8 = c(5000L, NA_integer_),
    RT4_ZT10 = c(5000L, NA_integer_), RT8_ZT14 = c(5000L, NA_integer_),
    Source = c("Prespecified fixed comparator (manuscript)",
               "Estimated per comparison in its child run: maximum likelihood of the Leading Edge / Remainder dispersion model on the actual NormEVS split (two DESeq2 trends, two prior widths)"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(cutoff_manifest, file.path(multi_root, "Cutoff_Method_Manifest.csv"), row.names = FALSE)

  run_cutoff_child <- function(method) {
    message("=====================================================")
    message("Starting cutoff method: ", method)
    status <- system2(
      command = file.path(R.home("bin"), "Rscript"),
      args = shQuote(this_script),
      env = c(
        "SEQUENCE_CHILD_RUN=1",
        paste0("SEQUENCE_CUTOFF_METHOD=", method),
        paste0("SEQUENCE_MULTI_ROOT=", multi_root),
        paste0("SEQUENCE_REPO_ROOT=", repo_guess)
      )
    )
    if (!identical(as.integer(status), 0L)) {
      stop("Cutoff-method run failed: ", method, " (exit status ", status, ")")
    }
    invisible(TRUE)
  }
  for (method in c("fixed5000", "ml_width")) run_cutoff_child(method)

  cross_dir <- file.path(multi_root, "Cross_Method_Comparison")
  dir.create(cross_dir, recursive = TRUE, showWarnings = FALSE)
  ml_sel_path <- file.path(multi_root, "Likelihood_Width", "Summary_Tables", "Table_Likelihood_Cutoff_Selected.csv")
  ml_sel <- if (file.exists(ml_sel_path)) utils::read.csv(ml_sel_path, stringsAsFactors = FALSE) else NULL
  if (!is.null(ml_sel)) {
    for (cn in intersect(cmp_names, ml_sel$Comparison)) {
      cutoff_manifest[cutoff_manifest$Cutoff_Method == "Likelihood_Width", cn] <- as.integer(ml_sel$k_star[ml_sel$Comparison == cn])
    }
    utils::write.csv(cutoff_manifest, file.path(multi_root, "Cutoff_Method_Manifest.csv"), row.names = FALSE)
    utils::write.csv(ml_sel, file.path(cross_dir, "Table_Likelihood_Cutoff_Selected.csv"), row.names = FALSE)
  }
  try_fig(save_cutoff_comparison_figure(multi_root, cross_dir), "Fixed 5,000 vs Likelihood k*")

  final_zip <- file.path(dirname(multi_root), "SEQUENCE_FINAL_ALL_CUTOFF_METHODS.zip")
  if (file.exists(final_zip)) unlink(final_zip, force = TRUE)
  oldwd <- getwd(); on.exit(setwd(oldwd), add = TRUE)
  setwd(dirname(multi_root))
  zip_rel <- list.files(multi_root, recursive = TRUE, all.files = FALSE)
  zip_rel <- zip_rel[!grepl("SEQUENCE_ALL_(FIGURES|TABLES)\\.zip$", zip_rel)]
  zip_items <- file.path(basename(multi_root), zip_rel)
  zip_items <- zip_items[file.exists(zip_items)]
  utils::zip(zipfile = basename(final_zip), files = zip_items, flags = "-q")
  if (!file.exists(final_zip)) stop("Final multi-cutoff ZIP creation failed: ", final_zip)

  cat("\n=====================================================\n")
  cat("All runs complete: Fixed 5,000 and Likelihood k*.\n")
  if (!is.null(ml_sel)) for (i in seq_len(nrow(ml_sel))) cat(sprintf(
    "  %-9s Likelihood k* = %6d (support %d-%d)   log-lik gain vs unsplit = %.1f, vs k = 5,000 = %.1f   prior width %.0f%% narrower\n",
    ml_sel$Comparison[i], ml_sel$k_star[i], ml_sel$support_lo[i], ml_sel$support_hi[i],
    ml_sel$delta_loglik_vs_Original[i], ml_sel$delta_loglik_k_star_vs_5000[i], ml_sel$prior_width_reduction_pct[i]))
  cat("Final combined ZIP:\n", normalizePath(final_zip, winslash = "/", mustWork = TRUE), "\n", sep = "")
  cat("=====================================================\n\n")
  quit(save = "no", status = 0L)
}

# Active child-run cutoff method.
CUTOFF_METHOD_KEY <- Sys.getenv("SEQUENCE_CUTOFF_METHOD", unset = "fixed5000")
CUTOFF_METHOD_INFO <- switch(
  CUTOFF_METHOD_KEY,
  fixed5000 = list(slug = "Fixed_5000", label = "Fixed 5,000", basis = "Prespecified fixed top-5,000 PASs per condition"),
  ml_width = list(slug = "Likelihood_Width", label = "Likelihood k*", basis = "Per comparison, the k maximizing the likelihood of a two-group (Leading Edge / Remainder) DESeq2 dispersion model on the actual NormEVS split, after the varying-dispersion framework of Tzougas (2020); no p-values used"),
  stop("Unknown SEQUENCE_CUTOFF_METHOD: ", CUTOFF_METHOD_KEY)
)
CUTOFF_METHOD_SLUG <- CUTOFF_METHOD_INFO$slug
CUTOFF_METHOD_LABEL <- CUTOFF_METHOD_INFO$label

required_packages <- c(
  "DESeq2",
  "apeglm",
  "fdrtool",
  "ggplot2",
  "ggrepel",
  "dplyr",
  "tidyr",
  "gridExtra",
  "grid",
  "scales",
  "grDevices",
  "S4Vectors"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Required R package(s) are not installed: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this manuscript pipeline."
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(fdrtool)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(grid)
  library(scales)
  library(grDevices)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# SECTION 1 OF 5
# USER SETTINGS, METADATA, PATHS, AND GENERAL HELPERS
# =============================================================================

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------

count_file_candidates <- c(
  "WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  file.path("data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"),
  "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv",
  "WTTS-Seq_2022.2_DE_raw_read_numbers(20260822-183312).csv",
  "/mnt/data/WTTS-Seq_2022.2_DE_raw_read_numbers(20260822-183312).csv"
)

twas_file_candidates <- c(
  "3aTWAS_genes_of_11_brain_disorders.csv",
  file.path("data", "3aTWAS_genes_of_11_brain_disorders.csv"),
  "/root/REAPER98632/data/3aTWAS_genes_of_11_brain_disorders.csv",
  "/mnt/data/3aTWAS_genes_of_11_brain_disorders.csv"
)

TWAS_TARGET_SPECIES <- "rat"
TWAS_ORTHOLOG_MIN_SUPPORT <- 1L

# DESeq2 Standard and Strong tests use Benjamini-Hochberg FDR 10%.
# Weak-CNH uses a more permissive BH screen and becomes a final Weak discovery
# only after independent HBFSS support. HBFSS uses no DESeq2 adjusted-p-value gate.
BH_FDR_STANDARD <- 0.10
BH_FDR_STRONG <- 0.10
BH_FDR_WEAK <- 0.20

lfc_boundary <- 1.0

# Higher Criticism searches only the lower tail of the ordered empirical-null
# p-value distribution. alpha0 = 0.10 means that the maximum HC score is sought
# only among the lowest 10% of empirical p-values. A non-positive maximum HC
# score yields no HC threshold and therefore no HBFSS discoveries in that view.
HC_ALPHA0 <- 0.10

if (!is.numeric(HC_ALPHA0) || length(HC_ALPHA0) != 1L ||
    !is.finite(HC_ALPHA0) || HC_ALPHA0 <= 0 || HC_ALPHA0 > 1) {
  stop("HC_ALPHA0 must be a single finite number in (0, 1].")
}

# EVS selection size is comparison-specific and is defined in comparison_table
# below from the empirically estimated weighted-Pareto optimum (k*). No global
# active cutoff is selected by the multi-cutoff orchestration layer.
EMPIRICAL_EVS_CUTOFF_BASIS <- paste0(
  CUTOFF_METHOD_INFO$basis,
  "; locked before downstream DE significance testing"
)

figure_dpi <- 320

base_theme_size <- 10

n_top_labels_volcano <- 20L

# Manuscript export behavior. When TRUE, the script writes only the focused
# paper-ready figure panels into the manuscript output tree.
# Tables are deliberately concise: significant sites plus compact method/count
# summaries, with the unsplit Original dataset represented once.
EXPORT_ONLY_PAPER_FIGURES <- FALSE
EXPORT_SUPPORT_FIGURES <- TRUE
EXPORT_INDIVIDUAL_VIEW_FIGURES <- FALSE

# -----------------------------------------------------------------------------
# Statistical decision rules
# -----------------------------------------------------------------------------
#   Std       = ordinary DESeq2 Wald BH padj < 0.10 AND |apeglm LFC| >= 1
#   Strong    = DESeq2 greaterAbs(lfcThreshold = 1) BH padj < 0.10
#               AND |apeglm LFC| >= 1
#   Weak-CNH  = DESeq2 lessAbs(lfcThreshold = 1) BH padj < 0.20
#   HBFSS     = |apeglm LFC| * [-log10(empirical p)] > Htau
#   Weak      = Weak-CNH AND |apeglm LFC| < 1 AND HBFSS
#   Overlap   = HBFSS AND (Std OR Strong OR Weak-CNH)
#
# HCp is the higher-criticism empirical-p threshold used to calculate Htau.
# HC is maximized only over the lowest HC_ALPHA0 fraction of ordered empirical
# p-values and must attain a strictly positive maximum. HCp is not applied again
# as a second HBFSS significance gate.
#
# -----------------------------------------------------------------------------
# Palette
# -----------------------------------------------------------------------------

plot_palette <- list(
  background = "#BDBDBD",
  threshold = "#A65628",
  hc = "#A65628",
  hbfss_line = "#6A3D9A",
  weak = "#4EA3F1",
  strong = "#E31A1C",
  standard = "#33A02C",
  hbfss = "#6A3D9A",
  overlap = "#54278F",
  control = "#4D4D4D",
  treatment = "#1F78B4",
  histogram = "#969696"
)

# -----------------------------------------------------------------------------
# Short names used in file exports
# -----------------------------------------------------------------------------

dataset_short <- c(
  raw_dataset = "Raw",
  leading_edge_dataset = "Lead",
  remainder_dataset = "Rem"
)

track_short <- c(
  normalized_evs = "NormEVS",
  raw_evs = "RawEVS"
)

dataset_key_order <- c(
  "raw_dataset",
  "leading_edge_dataset",
  "remainder_dataset"
)

dataset_key_labels <- c(
  raw_dataset = "Original dataset",
  leading_edge_dataset = "Leading-edge dataset",
  remainder_dataset = "Remainder dataset"
)

# -----------------------------------------------------------------------------
# Embedded sample metadata
# -----------------------------------------------------------------------------

meta_all <- data.frame(
  id = c(
    "R0_1", "R0_2", "R0_3", "R0_4", "R0_5",
    "ZT6_1", "ZT6_2", "ZT6_3", "ZT6_4", "ZT6_5",
    "R2_1", "R2_2", "R2_3", "R2_4", "R2_5",
    "ZT8_1", "ZT8_2", "ZT8_3", "ZT8_4", "ZT8_5",
    "R4_1", "R4_2", "R4_3", "R4_4", "R4_5",
    "ZT10_1", "ZT10_2", "ZT10_3", "ZT10_4", "ZT10_5",
    "R8_1", "R8_2", "R8_3", "R8_4", "R8_5",
    "ZT14_1", "ZT14_2", "ZT14_3", "ZT14_4", "ZT14_5"
  ),
  condition = c(
    "treatment", "treatment", "treatment", "treatment", "treatment",
    "control", "control", "control", "control", "control",
    "treatment", "treatment", "treatment", "treatment", "treatment",
    "control", "control", "control", "control", "control",
    "treatment", "treatment", "treatment", "treatment", "treatment",
    "control", "control", "control", "control", "control",
    "treatment", "treatment", "treatment", "treatment", "treatment",
    "control", "control", "control", "control", "control"
  ),
  stringsAsFactors = FALSE
)

rownames(meta_all) <- meta_all$id

meta_all$condition <- factor(
  meta_all$condition,
  levels = c("control", "treatment")
)

levels(meta_all$condition) <- c("untrt", "trt")

# -----------------------------------------------------------------------------
# Four pairwise comparisons
# -----------------------------------------------------------------------------

cutoff_manifest_path <- file.path(
  Sys.getenv("SEQUENCE_MULTI_ROOT", unset = ""),
  "Cutoff_Method_Manifest.csv"
)
if (!nzchar(Sys.getenv("SEQUENCE_MULTI_ROOT", unset = "")) || !file.exists(cutoff_manifest_path)) {
  stop("Dynamic cutoff manifest is missing. Run this script normally through the master process; do not launch a child run directly.")
}
cutoff_manifest_runtime <- utils::read.csv(cutoff_manifest_path, stringsAsFactors=FALSE, check.names=FALSE)
required_manifest_methods <- c("Fixed_5000")
if (!all(required_manifest_methods %in% cutoff_manifest_runtime$Cutoff_Method)) {
  stop("Cutoff manifest does not contain all required methods: ", paste(required_manifest_methods, collapse=", "))
}
manifest_row <- function(method) cutoff_manifest_runtime[match(method, cutoff_manifest_runtime$Cutoff_Method),,drop=FALSE]
fixed_row <- manifest_row("Fixed_5000")

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix   = c("R0", "R2", "R4", "R8"),
  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),
  fixed_5000_k    = as.integer(unlist(fixed_row[c("RT0_ZT6","RT2_ZT8","RT4_ZT10","RT8_ZT14")], use.names=FALSE)),
  stringsAsFactors = FALSE
)

get_active_evs_cutoff <- function(comparison_name) {
  idx <- match(as.character(comparison_name), comparison_table$comparison_name)
  if (is.na(idx)) stop("No EVS cutoff is defined for comparison: ", comparison_name)
  k <- switch(
    CUTOFF_METHOD_KEY,
    fixed5000 = comparison_table$fixed_5000_k[idx],
    ml_width = if (is.null(LIKELIHOOD_K[[as.character(comparison_name)]])) NA_integer_ else LIKELIHOOD_K[[as.character(comparison_name)]]
  )
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || !is.finite(k) || k < 1L) {
    stop("Invalid EVS cutoff for ", comparison_name, " under ", CUTOFF_METHOD_LABEL)
  }
  k
}

# Backward-compatible function name used throughout the established pipeline.
get_empirical_evs_cutoff <- get_active_evs_cutoff

if (anyDuplicated(comparison_table$comparison_name)) stop("comparison_table contains duplicated comparison names.")
cutoff_cols <- c("fixed_5000_k")
if (any(vapply(comparison_table[cutoff_cols], function(x) any(is.na(x) | x < 1L), logical(1)))) {
  stop("Every comparison must have a positive cutoff for every cutoff method.")
}

# -----------------------------------------------------------------------------
# Repository and file helpers
# -----------------------------------------------------------------------------

resolve_existing_file <- function(candidates, label) {
  hits <- candidates[file.exists(candidates)]

  if (length(hits) == 0L) {
    stop(
      "Could not find ", label, ". Tried: ",
      paste(candidates, collapse = " | ")
    )
  }

  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0L) {
    candidate <- sub("^--file=", "", file_arg[1])
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  NA_character_
}

find_repo_root <- function() {
  candidates <- c(
    getwd(),
    dirname(getwd()),
    "/root/REAPER98632"
  )

  script_path <- get_script_path()

  if (!is.na(script_path)) {
    candidates <- c(dirname(script_path), candidates)
  }

  candidates <- unique(candidates[file.exists(candidates) | dir.exists(candidates)])

  for (cand in candidates) {
    if (dir.exists(file.path(cand, ".git"))) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }

  for (cand in candidates) {
    if (file.exists(file.path(cand, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")) ||
        file.exists(file.path(cand, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"))) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

repo_root <- find_repo_root()

multi_output_root <- Sys.getenv(
  "SEQUENCE_MULTI_ROOT",
  unset = file.path(repo_root, "exports", "sequence_final_cutoff_methods")
)
output_dir <- file.path(multi_output_root, CUTOFF_METHOD_SLUG)

if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE, force = TRUE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paper_fig_dir <- file.path(output_dir, "Combined_Figures")
summary_table_dir <- file.path(output_dir, "Summary_Tables")
dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_table_dir, recursive = TRUE, showWarnings = FALSE)

paper_registry <- list()

registry_key <- function(comparison_name, track_key) {
  paste(comparison_name, track_key, sep = "__")
}

should_write_figure <- function(path) {
  if (!isTRUE(EXPORT_ONLY_PAPER_FIGURES)) {
    return(TRUE)
  }

  target_dir <- normalizePath(
    paper_fig_dir,
    winslash = "/",
    mustWork = FALSE
  )

  path_dir <- normalizePath(
    dirname(path),
    winslash = "/",
    mustWork = FALSE
  )

  startsWith(path_dir, target_dir)
}



# -----------------------------------------------------------------------------
# Generic numeric and string helpers
# -----------------------------------------------------------------------------

assert_required_columns <- function(df, required_cols, object_name = "data frame") {
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns in ",
      object_name,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

safe_neglog10 <- function(x, pseudocount = 1e-12) {
  -log10(pmax(x, pseudocount))
}

clip_probabilities <- function(x, eps = 1e-300) {
  x <- unname(as.numeric(x))

  if (!length(x)) {
    return(numeric(0))
  }

  bad <- !is.finite(x) | is.na(x)
  x[bad] <- NA_real_

  good <- !is.na(x)
  x[good] <- pmin(pmax(x[good], eps), 1 - 1e-12)

  x
}

finite_plot_df <- function(df, x_col, y_col) {
  keep <- is.finite(df[[x_col]]) &
    !is.na(df[[x_col]]) &
    is.finite(df[[y_col]]) &
    !is.na(df[[y_col]])

  df[keep, , drop = FALSE]
}

compact_title <- function(x, width = 54) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

compact_caption <- function(x, width = 118) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

save_csv <- function(df, path) {
  write.csv(df, file = path, row.names = FALSE)
}

save_grob <- function(g, path, width = 14.0, height = 8.5, dpi = figure_dpi, bg = "white") {
  if (!should_write_figure(path)) {
    return(invisible(NULL))
  }

  ggplot2::ggsave(
    filename = path,
    plot = g,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = bg,
    limitsize = FALSE
  )
}

pretty_dataset_type <- function(dataset_key) {
  switch(
    dataset_key,
    raw_dataset = "Original dataset",
    leading_edge_dataset = "Leading-edge dataset",
    remainder_dataset = "Remainder dataset",
    dataset_key
  )
}

pretty_dataset_label <- function(dataset_name) {
  parts <- strsplit(dataset_name, "_", fixed = TRUE)[[1]]

  if (length(parts) < 4L) {
    return(dataset_name)
  }

  comparison_name <- paste(parts[1], parts[2], sep = "_")

  if (length(parts) >= 5L && parts[3] %in% unname(track_short)) {
    track_label <- parts[3]
    dataset_key <- paste(parts[4:length(parts)], collapse = "_")
    return(paste(comparison_name, track_label, pretty_dataset_type(dataset_key), sep = " | "))
  }

  dataset_key <- paste(parts[3:length(parts)], collapse = "_")

  paste(comparison_name, pretty_dataset_type(dataset_key), sep = " | ")
}


make_design_formula <- function(coldata) {
  ~ condition
}

get_condition_coef <- function(dds) {
  rn <- DESeq2::resultsNames(dds)
  idx <- grep("^condition_", rn)

  if (length(idx) == 0L) {
    stop("Could not identify condition coefficient in resultsNames(dds).")
  }

  rn[idx[1]]
}



resolve_top_n_cutoff <- function(sorted_values_desc, top_n) {
  n_total <- length(sorted_values_desc)
  top_n <- as.integer(top_n)

  if (n_total == 0L) {
    stop("resolve_top_n_cutoff() received an empty vector.")
  }

  if (!is.finite(top_n) || is.na(top_n) || top_n < 1L) {
    stop("top_n must be a positive integer.")
  }

  # Use the empirically predetermined comparison-specific k* exactly. Do not
  # silently shrink k* or replace the rank rule with a loading-value cutoff.
  if (n_total < top_n) {
    stop(
      "EVS requires exactly ", top_n,
      " PAS features per condition, but only ", n_total,
      " features are available for ranking."
    )
  }

  list(
    top_n_actual = top_n,
    cutoff_value = sorted_values_desc[top_n],
    n_total = n_total
  )
}

run_empirical_null_fdrtool <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  stat_vec <- unname(stat_vec)

  if (length(stat_vec) < 5L) {
    stop(sprintf("[%s] Fewer than 5 finite Wald statistics were available for fdrtool.", dataset_name))
  }

  fit <- fdrtool::fdrtool(
    stat_vec,
    statistic = "normal",
    plot = FALSE,
    verbose = FALSE,
    cutoff.method = "fndr"
  )
  fit$pval <- clip_probabilities(fit$pval)
  fit
}

safe_hc_thresh <- function(empirical_p, dataset_name) {
  sorted_empirical_p <- sort(
    clip_probabilities(empirical_p),
    na.last = NA,
    decreasing = FALSE
  )

  n_total <- length(sorted_empirical_p)

  if (n_total < 5L) {
    message(sprintf("[%s] Fewer than 5 empirical p-values were available for Higher Criticism.", dataset_name))
    return(NA_real_)
  }

  hc_scores <- suppressWarnings(
    tryCatch(
      fdrtool::hc.score(as.vector(sorted_empirical_p)),
      error = function(e) {
        message(
          sprintf(
            "[%s] hc.score failed: %s",
            dataset_name,
            conditionMessage(e)
          )
        )
        rep(NA_real_, n_total)
      }
    )
  )

  if (length(hc_scores) != n_total) {
    message(sprintf("[%s] Higher-Criticism score length mismatch.", dataset_name))
    return(NA_real_)
  }

  # Match fdrtool::hc.thresh(alpha0 = HC_ALPHA0): maximize HC only over the
  # lowest alpha0 fraction of ordered empirical p-values. The explicit score
  # calculation additionally permits the required positive-HC safeguard.
  # fdrtool::hc.thresh maximizes over 1:ceiling(alpha0 * n); use the same range.
  n_search <- max(1L, min(n_total, ceiling(HC_ALPHA0 * n_total)))
  search_idx <- seq_len(n_search)

  valid_idx <- search_idx[
    is.finite(hc_scores[search_idx]) &
      !is.na(hc_scores[search_idx]) &
      is.finite(sorted_empirical_p[search_idx]) &
      !is.na(sorted_empirical_p[search_idx]) &
      sorted_empirical_p[search_idx] > 0 &
      sorted_empirical_p[search_idx] < 1
  ]

  if (length(valid_idx) == 0L) {
    message(sprintf("[%s] No valid Higher-Criticism search points were available.", dataset_name))
    return(NA_real_)
  }

  best_idx <- valid_idx[which.max(hc_scores[valid_idx])]
  best_hc <- as.numeric(hc_scores[best_idx])
  hc_p <- as.numeric(sorted_empirical_p[best_idx])

  # Runtime consistency audit: when fdrtool::hc.thresh succeeds, our selected
  # lower-tail threshold must match hc.thresh(alpha0 = HC_ALPHA0). We calculate
  # scores explicitly only so that a non-positive maximum can be rejected.
  package_hc_p <- suppressWarnings(
    tryCatch(
      as.numeric(
        fdrtool::hc.thresh(
          as.vector(sorted_empirical_p),
          alpha0 = HC_ALPHA0,
          plot = FALSE
        )[1]
      ),
      error = function(e) NA_real_
    )
  )

  if (is.finite(package_hc_p) && !is.na(package_hc_p) &&
      !isTRUE(all.equal(hc_p, package_hc_p, tolerance = 1e-12))) {
    # Adopt the package threshold (the reference definition) and record it.
    pkg_idx <- which(sorted_empirical_p == package_hc_p)[1]
    warning(sprintf(
      "[%s] HC audit: explicit HCp %.17g differs from fdrtool::hc.thresh(alpha0=%.3f) %.17g; using the fdrtool value.",
      dataset_name, hc_p, HC_ALPHA0, package_hc_p), call. = FALSE)
    hc_p <- package_hc_p
    if (!is.na(pkg_idx)) best_hc <- as.numeric(hc_scores[pkg_idx])
  }

  # A non-positive maximum indicates no excess of small empirical p-values in
  # the prespecified HC search region. In that case no HC/HBFSS threshold is
  # asserted for the dataset/view.
  if (!is.finite(best_hc) || is.na(best_hc) || best_hc <= 0) {
    message(
      sprintf(
        "[%s] No positive Higher-Criticism signal within the lowest %.1f%% of empirical p-values; HBFSS disabled for this view.",
        dataset_name,
        100 * HC_ALPHA0
      )
    )
    return(NA_real_)
  }

  if (!is.finite(hc_p) || is.na(hc_p) || hc_p <= 0 || hc_p >= 1) {
    return(NA_real_)
  }

  message(
    sprintf(
      "[%s] HC calibration: alpha0=%.3f, HCmax=%.6f, HCp=%.8g, search=%d/%d empirical p-values.",
      dataset_name,
      HC_ALPHA0,
      best_hc,
      hc_p,
      n_search,
      n_total
    )
  )

  hc_p
}

# -----------------------------------------------------------------------------
# Plotting conventions
# -----------------------------------------------------------------------------

condition_shapes <- c(
  untrt = 21,
  trt = 24
)

condition_fills <- c(
  untrt = plot_palette$control,
  trt = plot_palette$treatment
)

condition_labels <- c(
  untrt = "Control",
  trt = "Treatment"
)

significance_method_levels <- c(
  "Standard",
  "Strong",
  "Weak",
  "HBFSS"
)

significance_method_labels <- c(
  "Standard" = "Std",
  "Strong" = "Str",
  "Weak" = "Weak",
  "HBFSS" = "HBFSS"
)

significance_method_sizes <- c(
  "Standard" = 2.85,
  "Strong" = 3.15,
  "Weak" = 3.05,
  "HBFSS" = 2.55
)

significance_method_colors <- c(
  "Standard" = plot_palette$standard,
  "Strong" = plot_palette$strong,
  "Weak" = plot_palette$weak,
  "HBFSS" = plot_palette$hbfss
)

# Filled DESeq2 symbols plus a star for HBFSS keep the shared key visually
# obvious in every panel, including reduced mobile views and exported legends.
significance_method_shapes <- c(
  "Standard" = 18,
  "Strong" = 15,
  "Weak" = 17,
  "HBFSS" = 8
)

build_significance_plot_long <- function(df) {
  rows <- list()

  add_method <- function(flag_col, method_name) {
    keep <- !is.na(df[[flag_col]]) & df[[flag_col]]
    if (!any(keep)) return(NULL)
    out <- df[keep, , drop = FALSE]
    out$Method <- method_name
    out
  }

  rows[["Standard"]] <- add_method("standard_flag", "Standard")
  rows[["Strong"]] <- add_method("strong_cnh_flag", "Strong")
  rows[["Weak"]] <- add_method("weak_significant_flag", "Weak")
  rows[["HBFSS"]] <- add_method("hbfss_flag", "HBFSS")
  rows <- Filter(Negate(is.null), rows)

  if (!length(rows)) {
    out <- df[0, , drop = FALSE]
    out$Method <- factor(character(0), levels = significance_method_levels)
    return(out)
  }

  out <- dplyr::bind_rows(rows)
  out$Method <- factor(out$Method, levels = significance_method_levels)

  draw_rank <- c(HBFSS = 1L, Standard = 2L, Strong = 3L, Weak = 4L)
  out$.draw_rank <- unname(draw_rank[as.character(out$Method)])
  out <- out[order(out$.draw_rank), , drop = FALSE]
  out$.draw_rank <- NULL
  out
}

plot_expand_xy <- function() {
  list(
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.10))),
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))
  )
}

manuscript_theme <- function() {
  theme_bw(base_size = base_theme_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_theme_size + 1,
        hjust = 0.5,
        lineheight = 1.00,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = base_theme_size - 1,
        hjust = 0.5,
        lineheight = 1.00,
        margin = margin(b = 5)
      ),
      plot.caption = element_text(
        size = base_theme_size - 3,
        hjust = 0.5,
        colour = "grey30",
        lineheight = 0.98,
        margin = margin(t = 6)
      ),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(1, 1, 1, 1),
      legend.spacing.x = unit(4, "pt"),
      legend.spacing.y = unit(1, "pt"),
      legend.text = element_text(size = base_theme_size - 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      plot.margin = margin(10, 12, 10, 10)
    )
}

shared_panel_legend <- function(plot_obj) {
  g <- ggplotGrob(plot_obj + theme(legend.position = "bottom"))
  guide_idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")

  if (length(guide_idx) == 0L) {
    return(NULL)
  }

  g$grobs[[guide_idx[1]]]
}

# A panel's real data can legitimately have zero rows for a method (e.g. no
# Weak-CNH discoveries in a given view). When that happens to be true of the
# specific panel a legend gets borrowed from, ggplot silently omits that
# method's marker glyph from the legend key even with drop = FALSE. To
# guarantee every manuscript legend always shows all four method markers,
# build the shared legend from a small synthetic dataset that always
# contains exactly one row per method, rather than reusing a real panel.
build_full_method_legend <- function() {
  dummy <- data.frame(
    x = rep(0, length(significance_method_levels)),
    y = rep(0, length(significance_method_levels)),
    Method = factor(significance_method_levels, levels = significance_method_levels)
  )

  p <- ggplot(dummy, aes(x = x, y = y, color = Method, shape = Method, size = Method)) +
    geom_point(alpha = 0.98, stroke = 0.90) +
    scale_color_manual(
      values = significance_method_colors,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method",
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          shape = unname(significance_method_shapes[significance_method_levels]),
          color = unname(significance_method_colors[significance_method_levels]),
          size = rep(3.4, length(significance_method_levels)),
          alpha = rep(1, length(significance_method_levels)),
          stroke = rep(0.85, length(significance_method_levels))
        )
      )
    ) +
    scale_shape_manual(
      values = significance_method_shapes,
      breaks = significance_method_levels,
      drop = FALSE,
      guide = "none"
    ) +
    scale_size_manual(
      values = significance_method_sizes,
      breaks = significance_method_levels,
      guide = "none"
    ) +
    manuscript_theme() +
    theme(legend.position = "bottom")

  shared_panel_legend(p)
}



strip_legend <- function(p) {
  p + theme(legend.position = "none")
}

assemble_one_legend_panel <- function(plot_list, panel_title, ncol = length(plot_list), width_legend = TRUE) {
  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0L) {
    return(NULL)
  }

  legend <- build_full_method_legend()
  no_legend <- lapply(plot_list, strip_legend)

  row <- do.call(
    gridExtra::arrangeGrob,
    c(no_legend, list(ncol = ncol))
  )

  if (is.null(legend)) {
    return(
      gridExtra::arrangeGrob(
        row,
        ncol = 1,
        top = grid::textGrob(
          panel_title,
          gp = grid::gpar(fontface = "bold", cex = 1.15)
        )
      )
    )
  }

  gridExtra::arrangeGrob(
    row,
    legend,
    ncol = 1,
    heights = c(12, 1.4),
    top = grid::textGrob(
      panel_title,
      gp = grid::gpar(fontface = "bold", cex = 1.15)
    )
  )
}

create_all_figures_zip <- function() {
  # The manuscript figure archive contains curated PNG panels only. Individual
  # per-view figures and duplicate PDF renditions remain available in the output
  # tree when generated, but are excluded from the archive to avoid redundancy.
  panel_patterns <- c(
    "/Panels/.*\\.png$",
    "/Combined_Figures/.*\\.png$",
    "/TWAS/figures/Figure_TWAS_Gene_Support\\.png$"
  )

  all_png <- list.files(
    output_dir,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.png$",
    ignore.case = TRUE
  )
  all_png <- all_png[file.info(all_png)$isdir %in% FALSE]

  normalized <- gsub("\\\\", "/", all_png)
  keep <- rep(FALSE, length(all_png))
  for (pat in panel_patterns) keep <- keep | grepl(pat, normalized)
  figure_files <- all_png[keep]

  # The global TWAS method-count figure is intentionally excluded because the
  # per-comparison five-view TWAS panels already contain the same count summary.
  figure_files <- figure_files[
    !grepl("Figure_TWAS_Method_Counts\\.png$", figure_files)
  ]

  if (!length(figure_files)) {
    stop("No curated manuscript PNG panels were found for SEQUENCE_ALL_FIGURES.zip.")
  }

  zip_path <- file.path(output_dir, "SEQUENCE_ALL_FIGURES.zip")
  if (file.exists(zip_path)) unlink(zip_path, force = TRUE)

  root_norm <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  file_norm <- normalizePath(figure_files, winslash = "/", mustWork = TRUE)
  relative_files <- substring(file_norm, nchar(root_norm) + 2L)
  relative_files <- sort(unique(relative_files))

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)
  utils::zip(zipfile = basename(zip_path), files = relative_files, flags = "-q")

  if (!file.exists(zip_path)) stop("Figure ZIP creation failed: ", zip_path)
  normalizePath(zip_path, winslash = "/", mustWork = TRUE)
}

create_all_tables_zip <- function() {
  # The manuscript table archive contains the concise tables used to interpret
  # significance, EVS/PCA evidence, and 3'aTWAS overlap. Per-view working tables
  # remain in their comparison folders but are not duplicated in the archive.
  wanted <- character(0)

  add_if_exists <- function(path) {
    if (file.exists(path)) wanted <<- c(wanted, path)
  }

  add_if_exists(file.path(summary_table_dir, "Table_DE_Method_Counts.csv"))
  add_if_exists(file.path(summary_table_dir, "Table_EVS_Empirical_Cutoffs.csv"))
  add_if_exists(file.path(summary_table_dir, "Table_EVS_Split_Audit.csv"))
  add_if_exists(file.path(summary_table_dir, "Table_EVS_PCA_Evidence.csv"))

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    add_if_exists(file.path(
      output_dir,
      comparison_name,
      paste0("Table_", comparison_name, "_Significant_Sites_All_Views.csv")
    ))
    add_if_exists(file.path(
      output_dir,
      comparison_name,
      "Table_Method_Counts_All_Views.csv"
    ))
  }

  add_if_exists(file.path(output_dir, "TWAS", "tables", "Table_TWAS_Overlap_PAS.csv"))
  add_if_exists(file.path(output_dir, "TWAS", "tables", "Table_TWAS_Overlap_Genes.csv"))
  add_if_exists(file.path(output_dir, "TWAS", "tables", "Table_TWAS_Method_Counts.csv"))

  wanted <- sort(unique(wanted))
  if (!length(wanted)) stop("No manuscript tables were found for SEQUENCE_ALL_TABLES.zip.")

  zip_path <- file.path(output_dir, "SEQUENCE_ALL_TABLES.zip")
  if (file.exists(zip_path)) unlink(zip_path, force = TRUE)

  root_norm <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  file_norm <- normalizePath(wanted, winslash = "/", mustWork = TRUE)
  relative_files <- substring(file_norm, nchar(root_norm) + 2L)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)
  utils::zip(zipfile = basename(zip_path), files = relative_files, flags = "-q")

  if (!file.exists(zip_path)) stop("Table ZIP creation failed: ", zip_path)
  normalizePath(zip_path, winslash = "/", mustWork = TRUE)
}

# =============================================================================
# SECTION 2 OF 5
# IMPORT COUNT MATRIX AND ANNOTATION
# =============================================================================

count_file <- resolve_existing_file(count_file_candidates, "WTTS count file")

message("Using count file: ", count_file)
message("Repository root: ", repo_root)
message("Output directory: ", output_dir)

WTTS_Seq <- read.csv(
  count_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

WTTS_Seq <- as.data.frame(
  WTTS_Seq,
  stringsAsFactors = FALSE
)

WTTS_Seq$OrigID <- as.character(WTTS_Seq$OrigID)
WTTS_Seq$Symbol <- as.character(WTTS_Seq$Symbol)

assert_required_columns(
  WTTS_Seq,
  c("OrigID", "Symbol"),
  object_name = "WTTS count file"
)

assert_required_columns(
  WTTS_Seq,
  meta_all$id,
  object_name = "WTTS count file sample columns"
)

coerce_count_column <- function(x) {
  suppressWarnings(
    as.numeric(
      gsub(
        ",",
        "",
        trimws(as.character(x)),
        fixed = TRUE
      )
    )
  )
}

for (sid in meta_all$id) {
  WTTS_Seq[[sid]] <- coerce_count_column(WTTS_Seq[[sid]])
}

WTTS_Seq <- WTTS_Seq[
  !is.na(WTTS_Seq$OrigID) & nzchar(trimws(WTTS_Seq$OrigID)),
  ,
  drop = FALSE
]

sample_na <- rowSums(is.na(WTTS_Seq[, meta_all$id, drop = FALSE])) > 0

if (any(sample_na)) {
  message("Removing ", sum(sample_na), " rows with missing or nonnumeric sample counts.")
}

WTTS_Seq <- WTTS_Seq[!sample_na, , drop = FALSE]

WTTS_Seq$feature_id <- make.unique(as.character(WTTS_Seq$OrigID), sep = "_dup")
rownames(WTTS_Seq) <- WTTS_Seq$feature_id

OrigID_Symbol <- data.frame(
  feature_id = WTTS_Seq$feature_id,
  orig_id = WTTS_Seq$OrigID,
  gene_symbol = WTTS_Seq$Symbol,
  stringsAsFactors = FALSE
)

OrigID_Symbol$feature_id <- as.character(OrigID_Symbol$feature_id)
OrigID_Symbol$orig_id <- as.character(OrigID_Symbol$orig_id)
OrigID_Symbol$gene_symbol <- as.character(OrigID_Symbol$gene_symbol)

OrigID_Symbol <- OrigID_Symbol %>%
  dplyr::mutate(
    gene_symbol = dplyr::if_else(
      is.na(gene_symbol),
      "",
      trimws(gene_symbol)
    )
  ) %>%
  dplyr::arrange(
    feature_id,
    dplyr::desc(gene_symbol != ""),
    gene_symbol
  ) %>%
  dplyr::distinct(
    feature_id,
    .keep_all = TRUE
  ) %>%
  dplyr::mutate(
    gene_symbol = dplyr::na_if(gene_symbol, "")
  )




# -----------------------------------------------------------------------------
# Gene/PAS annotation used by the final 3'aTWAS analysis
# -----------------------------------------------------------------------------

valid_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & nzchar(x) & x != "-" & grepl("[A-Za-z0-9]", x)
}

gene_key <- function(x) {
  x <- trimws(as.character(x))
  x[!valid_gene_symbol(x)] <- NA_character_
  toupper(x)
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  paste(sort(x), collapse = sep)
}

safe_min_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else min(x)
}

safe_max_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else max(x)
}

WTTS_gene_pas_summary <- OrigID_Symbol %>%
  dplyr::mutate(gene_key = gene_key(gene_symbol)) %>%
  dplyr::filter(!is.na(gene_key)) %>%
  dplyr::group_by(gene_key) %>%
  dplyr::summarise(
    WTTS_Gene = dplyr::first(gene_symbol[valid_gene_symbol(gene_symbol)]),
    WTTS_PAS_n = dplyr::n_distinct(orig_id),
    .groups = "drop"
  ) %>%
  dplyr::mutate(APA_multi_PAS = WTTS_PAS_n >= 2L)

# =============================================================================
# SECTION 3 OF 5
# EVS CONSTRUCTION AND EVS SUPPORT FIGURES
# =============================================================================

prepare_comparison_data <- function(comparison_name, group1_prefix, group2_prefix, WTTS_Seq, meta_all) {
  keep_ids <- grepl(paste0("^", group1_prefix, "_"), meta_all$id) |
    grepl(paste0("^", group2_prefix, "_"), meta_all$id)

  meta_sub <- meta_all[keep_ids, , drop = FALSE]

  coldata <- meta_sub[, c("condition"), drop = FALSE]

  sample_ids <- rownames(meta_sub)

  missing_samples <- setdiff(sample_ids, colnames(WTTS_Seq))

  if (length(missing_samples) > 0L) {
    stop(
      "Missing samples in WTTS file for ",
      comparison_name,
      ": ",
      paste(missing_samples, collapse = ", ")
    )
  }

  count_sub <- WTTS_Seq[, sample_ids, drop = FALSE]
  rownames(count_sub) <- rownames(WTTS_Seq)
  count_sub <- count_sub[rowSums(as.matrix(count_sub)) > 0, , drop = FALSE]

  stopifnot(all(colnames(count_sub) == rownames(coldata)))

  list(
    comparison_name = comparison_name,
    count_matrix = coerce_raw_count_matrix_for_deseq2(
      count_sub,
      context = paste0(comparison_name, " raw comparison matrix")
    ),
    coldata = coldata
  )
}


coerce_raw_count_matrix_for_deseq2 <- function(count_mat, context = "count matrix") {
  original_rownames <- rownames(count_mat)
  original_colnames <- colnames(count_mat)

  count_mat <- as.matrix(count_mat)

  if (!is.numeric(count_mat)) {
    suppressWarnings(storage.mode(count_mat) <- "numeric")
  }

  if (!is.numeric(count_mat)) {
    stop(context, " must be numeric raw counts before DESeq2.")
  }

  if (any(!is.finite(count_mat) | is.na(count_mat))) {
    stop(context, " contains NA or non-finite values after numeric coercion.")
  }

  if (any(count_mat < 0)) {
    stop(context, " contains negative values; DESeq2 requires non-negative counts.")
  }

  rounded <- round(count_mat)

  if (any(abs(count_mat - rounded) > 1e-6)) {
    warning(context, " contained non-integer values; values were rounded for DESeq2.")
  }

  if (any(rounded > .Machine$integer.max)) {
    stop(context, " contains counts larger than R integer storage can represent.")
  }

  storage.mode(rounded) <- "integer"
  rownames(rounded) <- original_rownames
  colnames(rounded) <- original_colnames

  rounded
}

# Median-of-ratios size factors for the full 10-sample RT/ZT comparison, reused
# by every downstream DESeq2 view (Original, Lead, Rem) so they share one
# normalization instead of re-estimating it from a selected PAS subset.
comparison_size_factors <- function(count_matrix, coldata) {
  dds_sf <- DESeq2::DESeqDataSetFromMatrix(
    countData = coerce_raw_count_matrix_for_deseq2(count_matrix, context = "comparison size-factor matrix"),
    colData = coldata,
    design = make_design_formula(coldata)
  )
  DESeq2::sizeFactors(DESeq2::estimateSizeFactors(dds_sf))
}

make_rank_matrix_for_track <- function(count_matrix, coldata, track_key) {
  track_key <- match.arg(track_key, c("normalized_evs", "raw_evs"))

  if (track_key == "raw_evs") {
    return(
      list(
        rank_matrix = as.data.frame(count_matrix),
        preprocessing_label = "Raw counts prior to EVS; PCA is performed directly on raw counts",
        size_factors = comparison_size_factors(count_matrix, coldata)
      )
    )
  }

  dds_init <- DESeq2::DESeqDataSetFromMatrix(
    countData = coerce_raw_count_matrix_for_deseq2(
      count_matrix,
      context = "NormEVS pre-split full comparison matrix"
    ),
    colData = coldata,
    design = make_design_formula(coldata)
  )

  dds_init <- DESeq2::estimateSizeFactors(dds_init)

  norm_counts <- as.data.frame(
    DESeq2::counts(dds_init, normalized = TRUE)
  )

  list(
    rank_matrix = norm_counts,
    preprocessing_label = "Median-of-ratios normalized counts prior to EVS; PCA is performed directly on normalized counts",
    size_factors = DESeq2::sizeFactors(dds_init)
  )
}

compute_pc1_loading_table <- function(value_df, sample_names, top_n, preprocessing_label = "Normalized prior to EVS") {
  x <- as.matrix(value_df[, sample_names, drop = FALSE])
  storage.mode(x) <- "numeric"

  if (ncol(x) < 2L) {
    stop("EVS PCA requires at least two samples in each condition.")
  }

  # Refined EVS method: run PCA directly on the pre-split matrix for the
  # condition, extract the PC1 feature eigenvector/loading, take absolute values,
  # rank descending, and select the comparison-specific empirical k*. No log
  # transformation is applied before PCA.
  pca_fit <- stats::prcomp(
    t(x),
    center = TRUE,
    scale. = FALSE,
    rank. = 2
  )

  loading_abs <- abs(pca_fit$rotation[, 1])

  loading_tbl <- data.frame(
    feature_id = names(loading_abs),
    pc1_loading_abs = unname(loading_abs),
    stringsAsFactors = FALSE
  )

  loading_tbl <- loading_tbl[
    order(loading_tbl$pc1_loading_abs, decreasing = TRUE),
    ,
    drop = FALSE
  ]

  loading_tbl$rank <- seq_len(nrow(loading_tbl))

  cutoff_info <- resolve_top_n_cutoff(
    loading_tbl$pc1_loading_abs,
    top_n = top_n
  )

  loading_tbl$split_class <- ifelse(
    loading_tbl$rank <= cutoff_info$top_n_actual,
    "high_loading",
    "background_loading"
  )

  list(
    pca_fit = pca_fit,
    loading_table = loading_tbl,
    cutoff = cutoff_info$cutoff_value,
    top_n_used = cutoff_info$top_n_actual,
    preprocessing_label = preprocessing_label
  )
}

build_eigenvector_split <- function(comparison_name, count_matrix, coldata, track_key) {
  track_key <- match.arg(track_key, c("normalized_evs", "raw_evs"))
  empirical_k <- get_empirical_evs_cutoff(comparison_name)

  # The EVS matrix is used only to choose feature IDs. The normalized track
  # follows the supplied workflow: median-of-ratios normalization of the full
  # RT/ZT comparison matrix, followed by condition-specific PCA/PC1 eigenvectors.
  # RawEVS repeats the same PCA/ranking rule without pre-split normalization.
  # Both tracks use the same comparison-specific empirical k* so preprocessing
  # is the only difference in the selection rule. Downstream DESeq2 always
  # receives the corresponding raw-count subset.
  raw_count_matrix <- coerce_raw_count_matrix_for_deseq2(
    count_matrix,
    context = paste0(track_key, " full comparison matrix before EVS")
  )

  if (nrow(raw_count_matrix) < empirical_k) {
    stop(
      comparison_name, " ", track_key,
      ": empirical EVS k*=", empirical_k,
      " exceeds the ", nrow(raw_count_matrix),
      " PASs available after zero-row filtering."
    )
  }

  rank_obj <- make_rank_matrix_for_track(
    count_matrix = raw_count_matrix,
    coldata = coldata,
    track_key = track_key
  )

  rank_matrix <- rank_obj$rank_matrix
  preprocessing_label <- rank_obj$preprocessing_label

  sample_ids <- colnames(raw_count_matrix)
  trt_ids <- sample_ids[coldata$condition == "trt"]
  untrt_ids <- sample_ids[coldata$condition == "untrt"]

  fit_trt <- compute_pc1_loading_table(
    rank_matrix,
    trt_ids,
    top_n = empirical_k,
    preprocessing_label = preprocessing_label
  )

  fit_untrt <- compute_pc1_loading_table(
    rank_matrix,
    untrt_ids,
    top_n = empirical_k,
    preprocessing_label = preprocessing_label
  )

  trt_high <- as.character(
    fit_trt$loading_table$feature_id[
      fit_trt$loading_table$rank <= empirical_k
    ]
  )

  untrt_high <- as.character(
    fit_untrt$loading_table$feature_id[
      fit_untrt$loading_table$rank <= empirical_k
    ]
  )

  if (length(trt_high) != empirical_k || length(untrt_high) != empirical_k) {
    stop(
      comparison_name, " ", track_key,
      ": EVS failed to select exactly empirical k*=", empirical_k,
      " PAS features per condition."
    )
  }

  # Joint/disjoint classification from the two independently ranked PC1 axes.
  joint_ids <- intersect(trt_high, untrt_high)
  disjoint_trt_ids <- setdiff(trt_high, untrt_high)
  disjoint_untrt_ids <- setdiff(untrt_high, trt_high)

  leading_edge_ids <- union(trt_high, untrt_high)
  remainder_ids <- setdiff(rownames(raw_count_matrix), leading_edge_ids)

  if (length(leading_edge_ids) == 0L) {
    stop("Leading-edge dataset is empty. Check sample mapping or EVS inputs.")
  }

  if (length(remainder_ids) == 0L) {
    stop(
      comparison_name, " ", track_key,
      ": remainder dataset is empty after empirical EVS k*=", empirical_k, "."
    )
  }

  if (length(intersect(leading_edge_ids, remainder_ids)) != 0L ||
      !setequal(union(leading_edge_ids, remainder_ids), rownames(raw_count_matrix))) {
    stop(comparison_name, " ", track_key, ": EVS Lead/Rem partition audit failed.")
  }

  list(
    comparison_name = comparison_name,
    track_key = track_key,
    empirical_evs_k = empirical_k,
    cutoff_basis = EMPIRICAL_EVS_CUTOFF_BASIS,
    preprocessing_label = preprocessing_label,
    rank_matrix = rank_matrix,
    fit_trt = fit_trt,
    fit_untrt = fit_untrt,
    trt_top_ids = trt_high,
    untrt_top_ids = untrt_high,
    joint_ids = joint_ids,
    disjoint_trt_ids = disjoint_trt_ids,
    disjoint_untrt_ids = disjoint_untrt_ids,
    leading_edge_ids = leading_edge_ids,
    remainder_ids = remainder_ids,
    downstream_deseq2_input = "raw-count feature subsets for DESeq2",
    size_factors = rank_obj$size_factors,
    raw_dataset = raw_count_matrix,
    leading_edge_dataset = raw_count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset = raw_count_matrix[remainder_ids, , drop = FALSE]
  )
}

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df, size_factors = NULL) {
  design_formula <- make_design_formula(coldata)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = coerce_raw_count_matrix_for_deseq2(
      count_mat,
      context = paste0(dataset_name, " DESeq2 input after EVS split")
    ),
    colData = coldata,
    design = design_formula
  )

  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]

  # Shared comparison-level normalization: DESeq() keeps pre-set size factors.
  if (!is.null(size_factors)) {
    if (!all(colnames(dds) %in% names(size_factors))) {
      stop(dataset_name, ": comparison size factors do not cover every sample.")
    }
    DESeq2::sizeFactors(dds) <- size_factors[colnames(dds)]
  }

  dds <- DESeq2::DESeq(
    dds,
    betaPrior = FALSE
  )

  res <- DESeq2::results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    alpha = BH_FDR_STANDARD,
    pAdjustMethod = "BH"
  )

  res_strong <- DESeq2::results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "greaterAbs",
    alpha = BH_FDR_STRONG,
    pAdjustMethod = "BH"
  )

  res_weak <- DESeq2::results(
    dds,
    contrast = c("condition", "trt", "untrt"),
    lfcThreshold = lfc_boundary,
    altHypothesis = "lessAbs",
    alpha = BH_FDR_WEAK,
    pAdjustMethod = "BH"
  )

  res_all_df <- as.data.frame(res)
  res_all_df$feature_id <- as.character(rownames(res_all_df))

  valid_stat <- is.finite(res_all_df$stat) & !is.na(res_all_df$stat)

  stat_vec <- as.numeric(res_all_df$stat[valid_stat])
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]

  if (length(stat_vec) < 5L) {
    stop(
      sprintf(
        "[%s] Fewer than 5 finite Wald statistics were available for fdrtool.",
        dataset_name
      )
    )
  }

  fdr_fit <- run_empirical_null_fdrtool(
    stat_vec,
    dataset_name = dataset_name
  )

  res_df <- res_all_df

  n_valid <- sum(valid_stat)

  if (length(fdr_fit$pval) != n_valid) {
    stop(sprintf("[%s] fdrtool empirical-p output length mismatch.", dataset_name))
  }

  # Preserve DESeq2's ordinary Wald pvalue and BH-adjusted padj exactly as
  # returned by DESeq2. Empirical p-values are stored in a separate column and
  # are used only for HBFSS/higher-criticism calculations and plotting.
  res_df$empirical_p <- NA_real_
  res_df$empirical_p[valid_stat] <- as.numeric(fdr_fit$pval)

  coef_name <- get_condition_coef(dds)

  shr <- DESeq2::lfcShrink(
    dds,
    coef = coef_name,
    type = "apeglm"
  )

  shr_df <- as.data.frame(shr)
  shr_df$feature_id <- as.character(rownames(shr_df))

  res_df <- dplyr::left_join(
    res_df,
    shr_df[, c("feature_id", "log2FoldChange")],
    by = "feature_id",
    suffix = c("", "_shrunk")
  )

  colnames(res_df)[colnames(res_df) == "log2FoldChange_shrunk"] <- "lfc_shrunk"

  hc_p_threshold_dataset <- safe_hc_thresh(
    res_df$empirical_p,
    dataset_name = dataset_name
  )

  res_df$HBFSS <- abs(res_df$lfc_shrunk) *
    (-log10(pmax(res_df$empirical_p, 1e-300)))

  if (is.na(hc_p_threshold_dataset) ||
      !is.finite(hc_p_threshold_dataset) ||
      hc_p_threshold_dataset <= 0 ||
      hc_p_threshold_dataset > 1) {
    hbfss_threshold_dataset <- NA_real_
    res_df$HBFSS_core_pass <- FALSE
  } else {
    hbfss_threshold_dataset <- abs(log10(hc_p_threshold_dataset)) * lfc_boundary
    res_df$HBFSS_core_pass <- !is.na(res_df$HBFSS) &
      is.finite(res_df$HBFSS) &
      res_df$HBFSS > hbfss_threshold_dataset
  }

  res_df$regulation_direction <- ifelse(
    is.na(res_df$lfc_shrunk),
    NA_character_,
    ifelse(
      res_df$lfc_shrunk > 0,
      "upregulated",
      ifelse(res_df$lfc_shrunk < 0, "downregulated", "no_change")
    )
  )

  res_strong_df <- as.data.frame(res_strong)
  res_strong_df$feature_id <- as.character(rownames(res_strong_df))

  res_weak_df <- as.data.frame(res_weak)
  res_weak_df$feature_id <- as.character(rownames(res_weak_df))

  res_df <- dplyr::left_join(
    res_df,
    res_strong_df[, c("feature_id", "pvalue", "padj")],
    by = "feature_id",
    suffix = c("", "_strong")
  )

  res_df <- dplyr::left_join(
    res_df,
    res_weak_df[, c("feature_id", "pvalue", "padj")],
    by = "feature_id",
    suffix = c("", "_weak")
  )

  colnames(res_df)[colnames(res_df) == "pvalue_strong"] <- "pvalue_strong_effect"
  colnames(res_df)[colnames(res_df) == "padj_strong"] <- "padj_strong_effect"
  colnames(res_df)[colnames(res_df) == "pvalue_weak"] <- "pvalue_weak_effect"
  colnames(res_df)[colnames(res_df) == "padj_weak"] <- "padj_weak_effect"

  res_df$resGA_pvalue <- res_df$pvalue_strong_effect
  res_df$resGA_padj <- res_df$padj_strong_effect
  res_df$resLA_pvalue <- res_df$pvalue_weak_effect
  res_df$resLA_padj <- res_df$padj_weak_effect

  # Final decision flags. DESeq2 ordinary/greaterAbs/lessAbs p-values are
  # kept separate from empirical-null p-values. The manuscript Standard effect
  # combines the ordinary DESeq2 Wald/BH result with the prespecified effect
  # reporting boundary applied to the apeglm-shrunken LFC.
  finite_lfc <- !is.na(res_df$lfc_shrunk) & is.finite(res_df$lfc_shrunk)
  abs_shrunk_lfc <- abs(res_df$lfc_shrunk)

  res_df$standard_flag <- finite_lfc &
    abs_shrunk_lfc >= lfc_boundary &
    !is.na(res_df$padj) &
    res_df$padj < BH_FDR_STANDARD

  res_df$strong_cnh_flag <- finite_lfc &
    abs_shrunk_lfc >= lfc_boundary &
    !is.na(res_df$resGA_padj) &
    res_df$resGA_padj < BH_FDR_STRONG

  res_df$weak_cnh_flag <- finite_lfc &
    abs_shrunk_lfc < lfc_boundary &
    !is.na(res_df$resLA_padj) &
    res_df$resLA_padj < BH_FDR_WEAK

  res_df$hbfss_flag <- finite_lfc &
    !is.na(res_df$HBFSS_core_pass) &
    res_df$HBFSS_core_pass

  # lessAbs alone is not reported as differential expression. Final Weak is the
  # intersection of sub-boundary lessAbs evidence and HBFSS signal evidence.
  res_df$weak_significant_flag <- res_df$weak_cnh_flag & res_df$hbfss_flag

  res_df$standard_hbfss_overlap <- res_df$standard_flag & res_df$hbfss_flag
  res_df$strong_hbfss_overlap <- res_df$strong_cnh_flag & res_df$hbfss_flag
  res_df$weak_hbfss_overlap <- res_df$weak_cnh_flag & res_df$hbfss_flag
  res_df$any_overlap <- res_df$hbfss_flag &
    (res_df$standard_flag | res_df$strong_cnh_flag | res_df$weak_cnh_flag)

  res_df$final_significant_flag <- res_df$standard_flag |
    res_df$strong_cnh_flag |
    res_df$weak_significant_flag |
    res_df$hbfss_flag

  # Standard/Strong are outside the |LFC|=1 reporting boundary, whereas final
  # Weak is inside it. These sets must therefore be disjoint by construction.
  if (any(res_df$standard_flag & res_df$weak_significant_flag, na.rm = TRUE)) {
    stop("Standard and Weak classifications overlapped in ", dataset_name)
  }
  if (any(res_df$strong_cnh_flag & res_df$weak_significant_flag, na.rm = TRUE)) {
    stop("Strong and Weak classifications overlapped in ", dataset_name)
  }
  if (any(res_df$weak_significant_flag & !res_df$hbfss_flag, na.rm = TRUE)) {
    stop("A final Weak PAS lacked HBFSS support in ", dataset_name)
  }

  # Exact decision-rule checks. These are runtime assertions only; they do not
  # create validation tables or alter the reported results.
  expected_standard <- finite_lfc &
    abs_shrunk_lfc >= lfc_boundary &
    !is.na(res_df$padj) &
    res_df$padj < BH_FDR_STANDARD

  expected_strong <- finite_lfc &
    abs_shrunk_lfc >= lfc_boundary &
    !is.na(res_df$resGA_padj) &
    res_df$resGA_padj < BH_FDR_STRONG

  expected_weak_cnh <- finite_lfc &
    abs_shrunk_lfc < lfc_boundary &
    !is.na(res_df$resLA_padj) &
    res_df$resLA_padj < BH_FDR_WEAK

  expected_weak <- expected_weak_cnh & res_df$hbfss_flag
  expected_overlap <- res_df$hbfss_flag &
    (expected_standard | expected_strong | expected_weak_cnh)

  stopifnot(identical(res_df$standard_flag, expected_standard))
  stopifnot(identical(res_df$strong_cnh_flag, expected_strong))
  stopifnot(identical(res_df$weak_cnh_flag, expected_weak_cnh))
  stopifnot(identical(res_df$weak_significant_flag, expected_weak))
  stopifnot(identical(res_df$any_overlap, expected_overlap))

  expected_hbfss <- abs(res_df$lfc_shrunk) *
    (-log10(pmax(res_df$empirical_p, 1e-300)))
  both_na <- is.na(expected_hbfss) & is.na(res_df$HBFSS)
  both_finite <- is.finite(expected_hbfss) & is.finite(res_df$HBFSS)
  close_enough <- rep(FALSE, length(expected_hbfss))
  close_enough[both_finite] <- abs(expected_hbfss[both_finite] - res_df$HBFSS[both_finite]) <=
    1e-12 * pmax(1, abs(expected_hbfss[both_finite]))
  if (!all(both_na | close_enough)) stop("HBFSS arithmetic validation failed for ", dataset_name)

  if (is.finite(hc_p_threshold_dataset) && !is.na(hc_p_threshold_dataset)) {
    expected_htau <- -log10(hc_p_threshold_dataset) * lfc_boundary
    if (!isTRUE(all.equal(hbfss_threshold_dataset, expected_htau, tolerance = 1e-12))) {
      stop("HBFSS threshold arithmetic validation failed for ", dataset_name)
    }
  }

  if (is.finite(hbfss_threshold_dataset) && !is.na(hbfss_threshold_dataset)) {
    expected_hbfss_flag <- finite_lfc & !is.na(res_df$HBFSS) &
      is.finite(res_df$HBFSS) & res_df$HBFSS > hbfss_threshold_dataset
    stopifnot(identical(res_df$hbfss_flag, expected_hbfss_flag))
  }

  base_mean_vec <- res_df$baseMean[!is.na(res_df$baseMean)]

  norm_counts <- as.data.frame(
    DESeq2::counts(dds, normalized = TRUE)
  )

  norm_counts$feature_id <- as.character(rownames(norm_counts))

  mm <- as.data.frame(S4Vectors::mcols(dds))
  mm$feature_id <- as.character(rownames(mm))

  disp_cols_available <- intersect(
    c(
      "feature_id",
      "dispGeneEst",
      "dispFit",
      "dispersion",
      "dispIter",
      "dispOutlier"
    ),
    colnames(mm)
  )

  disp_df <- mm[, disp_cols_available, drop = FALSE]

  annot_df$feature_id <- as.character(annot_df$feature_id)
  if (!"gene_symbol" %in% names(annot_df)) {
    annot_df$gene_symbol <- NA_character_
  }
  annot_df$gene_symbol <- as.character(annot_df$gene_symbol)

  annot_df <- annot_df %>%
    dplyr::mutate(
      gene_symbol = dplyr::if_else(
        is.na(gene_symbol),
        "",
        trimws(gene_symbol)
      )
    ) %>%
    dplyr::arrange(
      feature_id,
      dplyr::desc(gene_symbol != ""),
      gene_symbol
    ) %>%
    dplyr::distinct(
      feature_id,
      .keep_all = TRUE
    ) %>%
    dplyr::mutate(
      gene_symbol = dplyr::na_if(gene_symbol, "")
    )

  norm_counts <- norm_counts[!duplicated(norm_counts$feature_id), , drop = FALSE]
  disp_df <- disp_df[!duplicated(disp_df$feature_id), , drop = FALSE]

  final_df <- res_df %>%
    dplyr::left_join(annot_df, by = "feature_id") %>%
    dplyr::left_join(norm_counts, by = "feature_id") %>%
    dplyr::left_join(disp_df, by = "feature_id")

  final_df$neglog10_padj <- safe_neglog10(final_df$padj)
  final_df$neglog10_empirical_p <- safe_neglog10(final_df$empirical_p)

  final_df$dataset_name <- dataset_name
  final_df$hc_p_threshold_dataset <- hc_p_threshold_dataset
  final_df$hbfss_threshold_dataset <- hbfss_threshold_dataset

  preferred_cols <- c(
    "dataset_name",
    "feature_id",
    "orig_id",
    "gene_symbol",
    "baseMean",
    "log2FoldChange",
    "lfc_shrunk",
    "regulation_direction",
    "stat",
    "pvalue",
    "padj",
    "empirical_p",
    "HBFSS",
    "hc_p_threshold_dataset",
    "hbfss_threshold_dataset",
    "resLA_pvalue",
    "resLA_padj",
    "resGA_pvalue",
    "resGA_padj",
    "weak_cnh_flag",
    "strong_cnh_flag",
    "standard_flag",
    "hbfss_flag",
    "weak_significant_flag",
    "final_significant_flag",
    "weak_hbfss_overlap",
    "strong_hbfss_overlap",
    "standard_hbfss_overlap",
    "any_overlap"
  )

  final_df <- final_df[
    ,
    c(
      intersect(preferred_cols, names(final_df)),
      setdiff(names(final_df), preferred_cols)
    ),
    drop = FALSE
  ]

  list(
    dds = dds,
    results = final_df,
    base_mean_vec = base_mean_vec,
    hc_p_threshold = hc_p_threshold_dataset,
    hbfss_threshold = hbfss_threshold_dataset
  )
}

build_final_volcano_df <- function(df, y_col = "neglog10_empirical_p") {
  df <- finite_plot_df(df, "lfc_shrunk", y_col)

  # Volcano labels are gene symbols only. Numeric PAS/feature identifiers are
  # never substituted into the figure when a gene symbol is absent.
  df$gene_symbol_plot <- if ("gene_symbol" %in% names(df)) {
    trimws(as.character(df$gene_symbol))
  } else {
    rep(NA_character_, nrow(df))
  }

  df$has_valid_gene_symbol <- !is.na(df$gene_symbol_plot) &
    nzchar(df$gene_symbol_plot) &
    grepl("[A-Za-z]", df$gene_symbol_plot) &
    !df$gene_symbol_plot %in% c("-", ".", "NA", "N/A")

  df[order(df$neglog10_empirical_p, na.last = TRUE), , drop = FALSE]
}

select_final_volcano_labels <- function(df, y_col, n_labels = n_top_labels_volcano) {
  if (!nrow(df)) return(df[0, , drop = FALSE])

  lab_df <- df[
    df$has_valid_gene_symbol &
      !is.na(df$final_significant_flag) &
      df$final_significant_flag,
    ,
    drop = FALSE
  ]

  if (!nrow(lab_df)) return(lab_df[0, , drop = FALSE])

  # Rank final discoveries by the evidence actually used on the HBFSS plotting
  # coordinates, then by effect magnitude. This does not change significance.
  emp <- suppressWarnings(as.numeric(lab_df$empirical_p))
  emp[!is.finite(emp)] <- Inf
  hscore <- suppressWarnings(as.numeric(lab_df$HBFSS))
  hscore[!is.finite(hscore)] <- -Inf

  ord <- order(emp, -hscore, -abs(lab_df$lfc_shrunk), na.last = TRUE)
  lab_df <- lab_df[ord, , drop = FALSE]
  lab_df <- lab_df[seq_len(min(as.integer(n_labels), nrow(lab_df))), , drop = FALSE]

  pas_id <- if ("orig_id" %in% names(lab_df)) {
    as.character(lab_df$orig_id)
  } else {
    as.character(lab_df$feature_id)
  }
  bad_pas <- is.na(pas_id) | !nzchar(trimws(pas_id))
  pas_id[bad_pas] <- as.character(lab_df$feature_id[bad_pas])

  gene_label <- as.character(lab_df$gene_symbol_plot)
  duplicate_gene <- duplicated(gene_label) | duplicated(gene_label, fromLast = TRUE)
  lab_df$plot_label <- ifelse(
    duplicate_gene,
    paste0(gene_label, " [", pas_id, "]"),
    gene_label
  )

  lab_df
}


make_hbfss_boundary_df <- function(plot_df, hbfss_threshold, y_limit) {
  if (!is.finite(hbfss_threshold) || is.na(hbfss_threshold) || hbfss_threshold < 0) {
    return(NULL)
  }

  x_max <- max(
    max(abs(plot_df$lfc_shrunk), na.rm = TRUE),
    lfc_boundary * 1.1
  )

  x_min <- max(0.05, hbfss_threshold / max(y_limit, 1e-6))

  x_abs <- seq(
    x_min,
    x_max,
    length.out = 600
  )

  # HBFSS significance is score-threshold driven. HCp is shown separately as a
  # reference line because it is used to calculate Htau, not as a second gate.
  y_curve <- hbfss_threshold / x_abs

  keep <- is.finite(y_curve) &
    y_curve >= 0 &
    y_curve <= y_limit

  if (!any(keep)) {
    return(NULL)
  }

  x_abs <- x_abs[keep]
  y_curve <- y_curve[keep]

  rbind(
    data.frame(x = -rev(x_abs), y = rev(y_curve)),
    data.frame(x = x_abs, y = y_curve)
  )
}

plot_final_volcano <- function(df, dataset_name, short_title = NULL, label_genes = TRUE, y_limit_override = NULL, n_labels = n_top_labels_volcano) {
  plot_df <- build_final_volcano_df(df, y_col = "neglog10_empirical_p")
  if (!nrow(plot_df)) stop("No finite volcano plotting rows for ", dataset_name)

  method_df <- build_significance_plot_long(plot_df)
  lab_df <- if (isTRUE(label_genes)) {
    select_final_volcano_labels(plot_df, y_col = "neglog10_empirical_p", n_labels = n_labels)
  } else {
    plot_df[0, , drop = FALSE]
  }

  hc_raw <- suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1]))
  htau <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))
  hc_y <- if (is.finite(hc_raw) && !is.na(hc_raw) && hc_raw > 0 && hc_raw <= 1) {
    safe_neglog10(hc_raw)
  } else NA_real_

  # A shared y_limit_override (computed once across every comparison/view)
  # keeps all manuscript volcano panels on the same vertical scale, so
  # significance magnitudes are visually comparable across timepoints
  # instead of each panel silently rescaling to its own local maximum.
  y_limit <- if (!is.null(y_limit_override) && is.finite(y_limit_override)) {
    y_limit_override
  } else {
    max(plot_df$neglog10_empirical_p, na.rm = TRUE) * 1.05
  }
  boundary_df <- make_hbfss_boundary_df(plot_df, htau, y_limit)

  std_n <- sum(df$standard_flag, na.rm = TRUE)
  str_n <- sum(df$strong_cnh_flag, na.rm = TRUE)
  wk_n <- sum(df$weak_significant_flag, na.rm = TRUE)
  hbfss_n <- sum(df$hbfss_flag, na.rm = TRUE)
  ovlp_n <- sum(df$any_overlap, na.rm = TRUE)

  count_text <- paste0(
    "Std=", std_n,
    "  Str=", str_n,
    "  Wk=", wk_n,
    "  HBFSS=", hbfss_n,
    "  Ovlp=", ovlp_n,
    if (is.finite(hc_raw) && !is.na(hc_raw)) paste0("  HCp=", signif(hc_raw, 3)) else "",
    if (is.finite(htau) && !is.na(htau)) paste0("  Hτ=", signif(htau, 3)) else ""
  )

  plot_title <- if (is.null(short_title)) {
    compact_title(pretty_dataset_label(dataset_name), width = 42)
  } else short_title

  p <- ggplot() +
    geom_point(
      data = plot_df,
      aes(x = lfc_shrunk, y = neglog10_empirical_p),
      color = plot_palette$background,
      shape = 16,
      size = 0.48,
      alpha = 0.22
    ) +
    geom_point(
      data = method_df,
      aes(
        x = lfc_shrunk,
        y = neglog10_empirical_p,
        color = Method,
        shape = Method,
        size = Method
      ),
      alpha = 0.98,
      stroke = 0.90
    ) +
    scale_color_manual(
      values = significance_method_colors,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method",
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          shape = unname(significance_method_shapes[significance_method_levels]),
          color = unname(significance_method_colors[significance_method_levels]),
          size = rep(3.4, length(significance_method_levels)),
          alpha = rep(1, length(significance_method_levels)),
          stroke = rep(0.85, length(significance_method_levels))
        )
      )
    ) +
    scale_shape_manual(
      values = significance_method_shapes,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method",
      guide = "none"
    ) +
    scale_size_manual(
      values = significance_method_sizes,
      breaks = significance_method_levels,
      guide = "none"
    ) +
    geom_vline(
      xintercept = c(-lfc_boundary, lfc_boundary),
      linetype = "dashed",
      linewidth = 0.55,
      colour = plot_palette$threshold
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      linewidth = 0.35,
      colour = "grey55"
    ) +
    labs(
      title = plot_title,
      x = "apeglm shrunken log2FC",
      y = expression(-log[10]("Empirical p")),
      caption = compact_caption(count_text, width = 96)
    ) +
    coord_cartesian(clip = "off", ylim = c(0, y_limit)) +
    manuscript_theme() +
    plot_expand_xy() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.caption = element_text(
        size = base_theme_size - 2,
        hjust = 0.5,
        margin = margin(t = 4)
      ),
      plot.caption.position = "plot",
      plot.margin = margin(10, 12, 10, 10)
    )

  if (!is.na(hc_y) && is.finite(hc_y)) {
    p <- p + geom_hline(
      yintercept = hc_y,
      linetype = "dotted",
      linewidth = 0.65,
      color = plot_palette$hc
    )
  }

  if (!is.null(boundary_df)) {
    p <- p + geom_line(
      data = boundary_df,
      aes(x, y),
      inherit.aes = FALSE,
      color = plot_palette$hbfss_line,
      linewidth = 0.80
    )
  }

  if (nrow(lab_df) > 0L) {
    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      aes(x = lfc_shrunk, y = neglog10_empirical_p, label = plot_label),
      inherit.aes = FALSE,
      show.legend = FALSE,
      size = 2.05,
      color = "black",
      seed = 1,
      # A finite max.overlaps lets ggrepel silently drop the labels it
      # genuinely cannot place without collision, rather than forcing every
      # requested label onto the page and producing illegible overlapping
      # text in dense regions (this matters most in the narrow 5-panel
      # manuscript figures, where each subplot has limited width).
      max.overlaps = 15,
      force = 2.2,
      force_pull = 0.25,
      box.padding = 0.45,
      point.padding = 0.18,
      min.segment.length = 0,
      segment.alpha = 0.60,
      segment.size = 0.22
    )
  }

  p
}

run_one_evs_track <- function(comparison_name, track_key, count_matrix, coldata, annot_df) {
  evs <- build_eigenvector_split(
    comparison_name = comparison_name,
    count_matrix = count_matrix,
    coldata = coldata,
    track_key = track_key
  )

  message(
    comparison_name, " ", unname(track_short[track_key]),
    ": empirical EVS k*=", evs$empirical_evs_k,
    " selected independently from each condition PC1 eigenvector; ",
    "Lead union n=", length(evs$leading_edge_ids),
    ", Rem n=", length(evs$remainder_ids), "."
  )

  dataset_list <- list(
    raw_dataset = evs$raw_dataset,
    leading_edge_dataset = evs$leading_edge_dataset,
    remainder_dataset = evs$remainder_dataset
  )

  paper_registry[[registry_key(comparison_name, track_key)]] <<- list(
    comparison_name = comparison_name,
    track_key = track_key,
    evs = evs,
    dataset_list = dataset_list,
    coldata = coldata
  )

  # Original is independent of the EVS preprocessing path and is analyzed once.
  analysis_dataset_list <- dataset_list
  if (identical(track_key, "raw_evs")) analysis_dataset_list$raw_dataset <- NULL

  analysis_results <- list()
  for (nm in names(analysis_dataset_list)) {
    dataset_name <- paste(comparison_name, unname(track_short[track_key]), nm, sep = "_")
    fit <- run_core_analysis(
      count_mat = analysis_dataset_list[[nm]],
      coldata = coldata,
      dataset_name = dataset_name,
      annot_df = annot_df,
      size_factors = evs$size_factors
    )

    analysis_results[[nm]] <- list(
      dds = fit$dds,
      results = fit$results,
      dataset_mat = analysis_dataset_list[[nm]]
    )
  }

  reg <- paper_registry[[registry_key(comparison_name, track_key)]]
  reg$analysis_results <- analysis_results
  paper_registry[[registry_key(comparison_name, track_key)]] <<- reg
  TRUE
}

analysis_view_label <- function(track_key, dataset_key) {
  if (identical(dataset_key, "raw_dataset")) {
    return("Original (No EVS)")
  }

  paste0(
    unname(track_short[track_key]),
    " ",
    ifelse(dataset_key == "leading_edge_dataset", "Lead", "Rem")
  )
}

get_registered_result <- function(comparison_name, track_key, dataset_key) {
  obj <- paper_registry[[registry_key(comparison_name, track_key)]]

  if (is.null(obj) || is.null(obj$analysis_results) ||
      is.null(obj$analysis_results[[dataset_key]])) {
    return(NULL)
  }

  obj$analysis_results[[dataset_key]]$results
}

build_support_label <- function(df) {
  if (!nrow(df)) {
    return(character(0))
  }

  vapply(seq_len(nrow(df)), function(i) {
    tags <- character(0)

    if (isTRUE(df$standard_flag[i])) tags <- c(tags, "Std")
    if (isTRUE(df$strong_cnh_flag[i])) tags <- c(tags, "Strong")
    if (isTRUE(df$weak_significant_flag[i])) tags <- c(tags, "Weak")
    if (isTRUE(df$hbfss_flag[i])) tags <- c(tags, "HBFSS")

    if (!length(tags)) "" else paste(tags, collapse = "+")
  }, character(1))
}

comparison_analysis_views <- function() {
  data.frame(
    track_key = c(
      "normalized_evs",
      "normalized_evs", "normalized_evs",
      "raw_evs", "raw_evs"
    ),
    dataset_key = c(
      "raw_dataset",
      "leading_edge_dataset", "remainder_dataset",
      "leading_edge_dataset", "remainder_dataset"
    ),
    stringsAsFactors = FALSE
  )
}

analysis_view_folder_name <- function(track_key, dataset_key) {
  if (identical(dataset_key, "raw_dataset")) return("Original_No_EVS")
  paste0(
    unname(track_short[track_key]),
    "_",
    ifelse(dataset_key == "leading_edge_dataset", "Lead", "Rem")
  )
}

analysis_view_paths <- function(comparison_name, track_key, dataset_key) {
  base <- file.path(
    output_dir,
    comparison_name,
    analysis_view_folder_name(track_key, dataset_key)
  )
  fig <- file.path(base, "figures")
  tab <- file.path(base, "tables")
  dir.create(fig, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab, recursive = TRUE, showWarnings = FALSE)
  list(base = base, figures = fig, tables = tab)
}

export_comparison_manuscript_tables <- function(comparison_name) {
  views <- comparison_analysis_views()
  comparison_counts <- list()
  comparison_sig_rows <- list()

  for (i in seq_len(nrow(views))) {
    track_key <- views$track_key[i]
    dataset_key <- views$dataset_key[i]
    df <- get_registered_result(comparison_name, track_key, dataset_key)
    if (is.null(df) || !nrow(df)) next

    paths <- analysis_view_paths(comparison_name, track_key, dataset_key)
    analysis_label <- analysis_view_label(track_key, dataset_key)
    hc <- suppressWarnings(as.numeric(df$hc_p_threshold_dataset[1]))
    htau <- suppressWarnings(as.numeric(df$hbfss_threshold_dataset[1]))

    method_counts <- data.frame(
      Comparison = comparison_name,
      Analysis = analysis_label,
      EVS_k_per_condition = if (identical(dataset_key, "raw_dataset")) NA_integer_ else get_empirical_evs_cutoff(comparison_name),
      PAS_tested = nrow(df),
      HC_alpha0 = HC_ALPHA0,
      HCp = hc,
      Htau = htau,
      Std = sum(df$standard_flag, na.rm = TRUE),
      Strong = sum(df$strong_cnh_flag, na.rm = TRUE),
      Weak = sum(df$weak_significant_flag, na.rm = TRUE),
      HBFSS = sum(df$hbfss_flag, na.rm = TRUE),
      Ovlp = sum(df$any_overlap, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    expected_overlap <- sum(
      df$hbfss_flag & (df$standard_flag | df$strong_cnh_flag | df$weak_cnh_flag),
      na.rm = TRUE
    )
    if (method_counts$Ovlp[1] != expected_overlap) {
      stop("Overlap-count export mismatch: ", comparison_name, " / ", analysis_label)
    }

    save_csv(method_counts, file.path(paths$tables, "Table_Method_Counts.csv"))
    comparison_counts[[length(comparison_counts) + 1L]] <- method_counts

    sig <- df[
      !is.na(df$final_significant_flag) & df$final_significant_flag,
      ,
      drop = FALSE
    ]

    expected_sig_n <- sum(df$final_significant_flag, na.rm = TRUE)
    if (nrow(sig) != expected_sig_n) {
      stop("Significant-PAS export mismatch: ", comparison_name, " / ", analysis_label)
    }

    if (nrow(sig)) {
      pas <- if ("orig_id" %in% names(sig)) as.character(sig$orig_id) else as.character(sig$feature_id)
      bad_pas <- is.na(pas) | !nzchar(trimws(pas))
      pas[bad_pas] <- as.character(sig$feature_id[bad_pas])

      gene <- if ("gene_symbol" %in% names(sig)) as.character(sig$gene_symbol) else rep(NA_character_, nrow(sig))
      gene[is.na(gene) | !nzchar(trimws(gene))] <- NA_character_

      sig_out <- data.frame(
        Comparison = comparison_name,
        Analysis = analysis_label,
        EVS_k_per_condition = if (identical(dataset_key, "raw_dataset")) NA_integer_ else get_empirical_evs_cutoff(comparison_name),
        PAS = pas,
        Gene = gene,
        Direction = as.character(sig$regulation_direction),
        Apeglm_LFC = as.numeric(sig$lfc_shrunk),
        Std_p = as.numeric(sig$pvalue),
        Std_BH = as.numeric(sig$padj),
        Strong_p = as.numeric(sig$resGA_pvalue),
        Strong_BH = as.numeric(sig$resGA_padj),
        Weak_p = as.numeric(sig$resLA_pvalue),
        Weak_BH = as.numeric(sig$resLA_padj),
        EmpP = as.numeric(sig$empirical_p),
        HBFSS = as.numeric(sig$HBFSS),
        HC_alpha0 = HC_ALPHA0,
        HCp = hc,
        Htau = htau,
        Std = as.logical(sig$standard_flag),
        Strong = as.logical(sig$strong_cnh_flag),
        Weak = as.logical(sig$weak_significant_flag),
        HBFSS_sig = as.logical(sig$hbfss_flag),
        Ovlp = as.logical(sig$any_overlap),
        Support = build_support_label(sig),
        stringsAsFactors = FALSE
      )

      sig_out <- sig_out[
        order(sig_out$EmpP, -sig_out$HBFSS, -abs(sig_out$Apeglm_LFC), na.last = TRUE),
        ,
        drop = FALSE
      ]
    } else {
      sig_out <- data.frame()
    }

    save_csv(sig_out, file.path(paths$tables, "Table_Significant_Sites.csv"))
    if (nrow(sig_out)) comparison_sig_rows[[length(comparison_sig_rows) + 1L]] <- sig_out

    # Individual view figures are optional. Journal-ready review uses the
    # comparison-level panels generated after all five analysis views are fit.
    if (isTRUE(EXPORT_INDIVIDUAL_VIEW_FIGURES)) {
      dataset_name <- paste(comparison_name, unname(track_short[track_key]), dataset_key, sep = "_")
      volcano <- plot_final_volcano(
        df,
        dataset_name = dataset_name,
        short_title = paste0(comparison_name, " | ", analysis_label),
        label_genes = TRUE
      )
      save_grob(volcano, file.path(paths$figures, "Volcano.png"), width = 8.2, height = 6.4)
      save_grob(volcano, file.path(paths$figures, "Volcano.pdf"), width = 8.2, height = 6.4)
    }
  }

  comparison_counts <- if (length(comparison_counts)) dplyr::bind_rows(comparison_counts) else data.frame()
  if (nrow(comparison_counts)) {
    save_csv(
      comparison_counts,
      file.path(output_dir, comparison_name, "Table_Method_Counts_All_Views.csv")
    )
  }

  comparison_sig <- if (length(comparison_sig_rows)) {
    dplyr::bind_rows(comparison_sig_rows)
  } else {
    data.frame()
  }
  save_csv(
    comparison_sig,
    file.path(
      output_dir,
      comparison_name,
      paste0("Table_", comparison_name, "_Significant_Sites_All_Views.csv")
    )
  )

  invisible(comparison_counts)
}

build_overall_manuscript_summary <- function() {
  rows <- list()
  views <- comparison_analysis_views()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (i in seq_len(nrow(views))) {
      track_key <- views$track_key[i]
      dataset_key <- views$dataset_key[i]
      df <- get_registered_result(comparison_name, track_key, dataset_key)
      if (is.null(df) || !nrow(df)) next

      rows[[length(rows) + 1L]] <- data.frame(
        Comparison = comparison_name,
        Analysis = analysis_view_label(track_key, dataset_key),
        EVS_k_per_condition = if (identical(dataset_key, "raw_dataset")) NA_integer_ else get_empirical_evs_cutoff(comparison_name),
        PAS_tested = nrow(df),
        Std = sum(df$standard_flag, na.rm = TRUE),
        Strong = sum(df$strong_cnh_flag, na.rm = TRUE),
        Weak = sum(df$weak_significant_flag, na.rm = TRUE),
        HBFSS = sum(df$hbfss_flag, na.rm = TRUE),
        Ovlp = sum(df$any_overlap, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) data.frame() else dplyr::bind_rows(rows)
}

write_methods_note <- function() {
  methods_text <- c(
    "# Manuscript Methods",
    "",
    "## Study unit and comparisons",
    "The analytical unit was the polyadenylation-site (PAS) feature defined by OrigID in the WTTS-Seq raw-count matrix. PASs were retained as separate observations throughout differential testing and gene symbols were retained as annotation. Four comparisons were analyzed: RT0 versus ZT6, RT2 versus ZT8, RT4 versus ZT10, and RT8 versus ZT14. Within each pairwise comparison, PASs with zero counts across every treatment and control sample were removed before EVS and differential-expression analysis; all remaining nonzero PASs were retained.",
    "",
    "## Eigenvector splitting",
    paste0(
      "This run used the ", CUTOFF_METHOD_LABEL, " EVS selection rule. ",
      "The active per-comparison cutoffs were: ",
      paste(
        paste0(comparison_table$comparison_name, " k=", vapply(comparison_table$comparison_name, get_active_evs_cutoff, integer(1))),
        collapse = ", "
      ),
      ". These cutoffs were locked before downstream differential-expression testing. ",
      "The same active cutoff was applied to NormEVS and RawEVS within each comparison so the two EVS tracks differed only in the matrix used for PC1 ranking, not in the number of PASs selected per condition. ",
      "NormEVS used DESeq2 median-of-ratios normalized counts before PCA, whereas RawEVS used raw counts before PCA. ",
      "Within each condition, prcomp was applied to the transposed feature-by-sample matrix with centering and without scaling. ",
      "PASs were ranked by absolute PC1 loading, and the comparison-specific k* highest-loading PASs were selected independently from the RT and ZT condition eigenvectors."
    ),
    "PASs present in both condition-specific top-k* sets were classified as Joint; PASs present in only one set were classified as Disjoint for that condition. Joint plus both Disjoint sets formed the Leading Edge, and all other PASs formed the Remainder. EVS defined PAS membership only. Downstream DESeq2 analyses received the corresponding raw-count subset and estimated normalization and dispersion parameters within that analysis view. The five reported views were Original (No EVS), NormEVS Lead, NormEVS Rem, RawEVS Lead, and RawEVS Rem.",
    "",
    "## DESeq2 model and log2 fold-change shrinkage",
    "Each analysis view was modeled independently with DESeq2 using a negative-binomial generalized linear model with design ~ condition and ZT/untrt as the reference. Median-of-ratios size factors were estimated once per comparison on all ten RT and ZT samples and supplied to every view; within each view DESeq2 estimated gene-wise dispersions, the mean-dispersion relationship, final dispersions, and Wald statistics for the RT/trt coefficient. Log2 fold changes were then shrunken with apeglm. The apeglm-shrunken log2 fold change was the reported effect estimate, the effect term in HBFSS, the direction indicator, and the x-coordinate for all significance figures.",
    "",
    "## Standard effect",
    sprintf("The Standard effect used the ordinary two-sided DESeq2 Wald p-value with Benjamini-Hochberg adjustment at FDR %.2f and the prespecified manuscript effect boundary on the apeglm-shrunken estimate. A PAS was reported as Standard when padj < %.2f and |apeglm LFC| >= %.1f. The ordinary DESeq2 Wald p-value/padj were retained separately from the empirical-null p-value used for HBFSS.", BH_FDR_STANDARD, BH_FDR_STANDARD, lfc_boundary),
    "",
    "## Strong (decoupled) composite-null effect",
    sprintf("Strong (decoupled) effects were tested with DESeq2 results(..., lfcThreshold=%.1f, altHypothesis='greaterAbs'), followed by Benjamini-Hochberg adjustment at FDR %.2f. Reported Strong PASs additionally had |apeglm LFC| >= %.1f so the plotted/reported shrunken effect remained outside the stated effect boundary.", lfc_boundary, BH_FDR_STRONG, lfc_boundary),
    "",
    "## Weak composite-null effect",
    sprintf("Weak-effect support was tested with DESeq2 results(..., lfcThreshold=%.1f, altHypothesis='lessAbs'). Weak-CNH support required Benjamini-Hochberg adjusted p-value < %.2f and |apeglm LFC| < %.1f. A lessAbs rejection alone was not reported as differential expression. A PAS was reported as final Weak only when the same sub-boundary PAS also passed HBFSS.", lfc_boundary, BH_FDR_WEAK, lfc_boundary),
    "",
    "## Empirical-null calibration, higher criticism, and HBFSS",
    sprintf(
      paste0(
        "Finite ordinary DESeq2 Wald statistics were supplied directly to fdrtool with statistic='normal' and cutoff.method='fndr' to estimate a zero-centered normal empirical null and the corresponding empirical-null p-values. ",
        "These empirical p-values were distinct from the ordinary DESeq2 Wald p-values and from the greaterAbs/lessAbs p-values. ",
        "Higher-Criticism scores were then calculated from the ordered empirical-null p-values with fdrtool::hc.score. ",
        "The HC maximization region was prespecified as the lowest %.0f%% of ordered empirical p-values (alpha0=%.2f), matching the lower-tail restriction provided by fdrtool::hc.thresh(alpha0=...). ",
        "Within that region, HCp was the empirical p-value at the maximum HC score. A dataset/view was assigned an HCp threshold only when that maximum HC score was strictly positive; if the maximum HC score was non-positive or no valid search point existed, HCp and Htau were left undefined and no PAS in that dataset/view was called HBFSS-significant."
      ),
      100 * HC_ALPHA0,
      HC_ALPHA0
    ),
    sprintf("For each PAS, HBFSS = |apeglm-shrunken LFC| x [-log10(empirical p)]. With c=%.1f, Htau = -log10(HCp) x c when a valid positive-HC threshold existed. A PAS was HBFSS-significant when HBFSS > Htau. HCp served only to derive Htau and was not imposed as an additional per-PAS significance gate. HBFSS significance did not use a DESeq2 adjusted-p-value gate.", lfc_boundary),
    "",
    "## Overlap",
    "Overlap was the number of unique HBFSS-significant PASs that also satisfied at least one DESeq2 criterion (Standard, Strong, or Weak-CNH). Overlap was a comparison quantity, not a separate significance test. Because final Weak required HBFSS, every final Weak PAS contributed to the overlap set.",
    "",
    "## PCA comparison after eigenvector splitting",
    "PCA comparison figures used the same matrices that defined EVS: full-comparison median-of-ratios normalized counts for NormEVS and raw counts for RawEVS. No additional log transformation or post-split re-normalization was applied. Original, Leading Edge, and Remainder were compared within the same preprocessing scale. For each view, the full PCA eigenspectrum was used to calculate the percentage of total variance explained by PC1 and PC2. PC1 eigenvalue retention was calculated as the first eigenvalue of the view divided by the first eigenvalue of the corresponding Original matrix. Original-PC1 loading energy captured by a feature subset was calculated as the sum of squared Original PC1 loadings for that subset divided by the sum of squared Original PC1 loadings for all PASs. Leading Edge and Remainder loading-energy fractions were required to sum to 100%. Treatment-control separation was summarized in the PC1-PC2 plane as centroid distance divided by pooled within-group root-mean-square distance. PCA panels reported PAS count, PC1 explained variance, PC1 eigenvalue retention, Original-PC1 loading-energy fraction, and the separation ratio.",
    "",
    "## Volcano figures and tables",
    sprintf("All significance volcanoes used the HBFSS plotting coordinates: apeglm-shrunken log2 fold change on the x-axis and -log10(empirical-null p) on the y-axis. Standard DESeq2, Strong greaterAbs, and Weak lessAbs significance were determined from their own DESeq2 p-values and BH-adjusted p-values, then their markers were projected onto these common HBFSS coordinates only for visual comparison. The HBFSS boundary y=Htau/|LFC| and HCp reference were drawn on the same axes. The same key was used throughout: Std = green diamond, Strong = red square, Weak = blue triangle, HBFSS = purple star. Multiple method markers were superimposed at the same PAS coordinate. LessAbs-only PASs remained background. Each volcano reported the number of significant PASs for Std, Strong, Weak, HBFSS, and their HBFSS/DESeq2 overlap and labeled at most %d final significant PASs.", n_top_labels_volcano),
    "Each comparison/view folder contained one significant-PAS table and one method-count table. Significant-PAS tables contained only the union of final Standard, Strong, Weak, and HBFSS discoveries and reported the DESeq2 p-values/BH-adjusted p-values used by the applicable DESeq2 tests together with empirical p, HBFSS, HC alpha0, HCp, Htau, apeglm LFC, and explicit method indicators. Method-count tables also reported HC alpha0, HCp, and Htau for view-level calibration auditing.",
    "",
    "## 3'aTWAS ortholog overlap",
    "Human 3'aTWAS gene symbols were mapped to rat gene symbols by combining database-supported babelgene human-to-rat ortholog mappings (top=FALSE) with direct case-insensitive symbol-equivalent matches present in the WTTS annotation. For every comparison and analysis view, mapped TWAS orthologs were intersected with Standard, Strong, Weak, and HBFSS WTTS discoveries. Concise TWAS tables reported the human TWAS symbol, rat ortholog, total TWAS records, distinct TWAS transcript count, TWAS multi-transcript/APA status, total WTTS PAS count, WTTS multi-PAS/APA status, significant WTTS PAS identifiers, significant multi-PAS status, and the method(s) and analysis view in which significance was observed.",
    "",
    "## Output organization",
    "Original (No EVS), NormEVS Lead, NormEVS Rem, RawEVS Lead, and RawEVS Rem were written to separate folders within each comparison for view-specific tables. Manuscript figures were organized as comparison-level panels so related volcano, PCA/EVS, and 3'aTWAS results could be reviewed side by side without redundant individual images. The curated figure archive SEQUENCE_ALL_FIGURES.zip contained manuscript PNG panels only. The separate archive SEQUENCE_ALL_TABLES.zip contained the locked empirical EVS cutoff table, EVS split audit table, overall differential-expression method-count table, quantitative EVS/PCA evidence table, one all-view significant-PAS table and method-count table for each comparison, and the combined 3'aTWAS PAS-, gene-, and method-summary tables."
  )

  writeLines(methods_text, file.path(output_dir, "Methods_Manuscript.md"), useBytes = TRUE)
  invisible(TRUE)
}


# =============================================================================
# DESeq2 DISPERSION AFTER EVS, AND ONE FIGURE PER CUTOFF METHOD
# =============================================================================
# Everything here is read from the DESeq2 objects fitted after the split
# (Original, NormEVS Lead/Rem, RawEVS Lead/Rem). Nothing is simulated.
#
# The point of the split: DESeq2 shrinks every gene-wise dispersion toward ONE
# mean-dispersion trend with ONE prior width for all PASs, i.e. it assumes PASs
# at the same mean share the same expected dispersion. When they do not, the
# unsplit prior is centred wrongly and is too wide. EVS separates the two
# populations; each subset gets its own trend and a narrower prior.
#
# Per view (all from the fitted DESeq2 object):
#   log_disp_resid_var  varLogDispEsts = (MAD of log(dispGeneEst/dispFit))^2
#   exp_var_log_disp    sampling variance of a log dispersion estimate with this
#                       design, trigamma((m - p) / 2)  (m samples, p coefficients)
#   prior_width         varLogDispEsts - exp_var_log_disp (not floored): the
#                       width of DESeq2's dispersion prior around the trend.
#                       DESeq2's own dispPriorVar is the same quantity floored
#                       at 0.25 and is reported alongside.
#
# Split score (selection criterion, NormEVS track):
#   W = (n_LE w_LE + n_REM w_REM) / (n_LE + n_REM), lower = narrower prior.

DISP_MIN_GENE_EST <- 1e-6
# Every class here is inside the Leading Edge, so one blue ramp (darkest = Joint).
REGIME_COLS <- c("Joint" = "#08306B", "Disjoint, opposite LE" = "#2171B5",
                 "Disjoint, opposite Divergence" = "#4EB3D3",
                 "Disjoint, opposite Remainder" = "#CFD3E6")

dispersion_metrics_rows <- list()
dispersion_trend_rows <- list()
dispersion_pas_rows <- list()
regime_composition_rows <- list()

cutoff_data_path <- function(name) {
  file.path(Sys.getenv("SEQUENCE_MULTI_ROOT", unset = ""), "Cutoff_Data", name)
}
read_cutoff_data <- function(name) {
  f <- cutoff_data_path(name)
  if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) else NULL
}

dispersion_view_metrics <- function(dds, orig_fun) {
  mc <- S4Vectors::mcols(dds)
  fin <- DESeq2::dispersions(dds)
  ok <- !mc$allZero & is.finite(mc$dispGeneEst) & mc$dispGeneEst >= DISP_MIN_GENE_EST &
        is.finite(mc$dispFit) & mc$dispFit > 0 & is.finite(fin) & mc$baseMean > 0
  lr <- log(mc$dispGeneEst[ok] / mc$dispFit[ok])
  dfun <- DESeq2::dispersionFunction(dds)
  v_deseq <- attr(dfun, "varLogDispEsts")
  prior_var <- attr(dfun, "dispPriorVar")
  fit_type <- attr(dfun, "fitType")
  v_recomp <- stats::mad(lr)^2
  v <- if (length(v_deseq) == 1L && is.finite(v_deseq)) as.numeric(v_deseq) else v_recomp
  X <- stats::model.matrix(DESeq2::design(dds), as.data.frame(SummarizedExperiment::colData(dds)))
  exp_var <- trigamma((ncol(dds) - ncol(X)) / 2)
  outl <- if (is.null(mc$dispOutlier)) rep(FALSE, length(ok)) else (mc$dispOutlier %in% TRUE)
  data.frame(
    n_PAS = sum(ok),
    log_disp_resid_var = v,
    log_disp_resid_var_recomputed = v_recomp,
    exp_var_log_disp = exp_var,
    prior_width = max(v - exp_var, 0),
    disp_prior_var = if (is.null(prior_var)) NA_real_ else as.numeric(prior_var),
    fit_type = if (is.null(fit_type)) NA_character_ else as.character(fit_type),
    dispersion_outlier_pct = 100 * mean(outl[ok]),
    median_dispersion = stats::median(fin[ok]),
    median_log_vs_original_trend = stats::median(log(fin[ok] / orig_fun(mc$baseMean[ok]))),
    median_abs_shrinkage = stats::median(abs(log(fin[ok] / mc$dispGeneEst[ok]))),
    stringsAsFactors = FALSE
  )
}

# Regime classes of the actual child split (NormEVS), using the shared
# normalized-count knots from the master. Knots are converted to depth from
# the top of each PC1 ranking (rank 1 = largest |PC1|), so they apply to the
# comparison-level ranking used by the split.
evs_regime_composition <- function(evs, knots) {
  N_m <- knots$N[1]; c1 <- knots$c1[1]; c2 <- knots$c2[1]
  d_le <- N_m - c2            # opposite depth <= d_le           : opposite arm in its LE regime
  d_div <- N_m - c1 + 1       # d_le < opposite depth <= d_div   : Divergence regime
  rt <- stats::setNames(evs$fit_trt$loading_table$rank, evs$fit_trt$loading_table$feature_id)
  ru <- stats::setNames(evs$fit_untrt$loading_table$rank, evs$fit_untrt$loading_table$feature_id)
  opp <- c(ru[evs$disjoint_trt_ids], rt[evs$disjoint_untrt_ids])
  c("Joint" = length(evs$joint_ids),
    "Disjoint, opposite LE" = sum(opp <= d_le),
    "Disjoint, opposite Divergence" = sum(opp > d_le & opp <= d_div),
    "Disjoint, opposite Remainder" = sum(opp > d_div))
}

collect_dispersion_after_evs <- function(comparison_name) {
  views <- comparison_analysis_views()
  norm_reg <- paper_registry[[registry_key(comparison_name, "normalized_evs")]]
  orig <- norm_reg$analysis_results$raw_dataset$dds
  if (is.null(orig)) return(invisible(NULL))
  orig_fun <- DESeq2::dispersionFunction(orig)
  for (i in seq_len(nrow(views))) {
    obj <- paper_registry[[registry_key(comparison_name, views$track_key[i])]]
    dds <- obj$analysis_results[[views$dataset_key[i]]]$dds
    if (is.null(dds)) next
    label <- analysis_view_label(views$track_key[i], views$dataset_key[i])
    m <- tryCatch(dispersion_view_metrics(dds, orig_fun), error = function(e) NULL)
    if (!is.null(m)) {
      dispersion_metrics_rows[[length(dispersion_metrics_rows) + 1L]] <<-
        cbind(data.frame(Comparison = comparison_name, View = label, stringsAsFactors = FALSE), m)
    }
    mc <- S4Vectors::mcols(dds)
    keep <- !mc$allZero & mc$baseMean > 0
    outl <- if (is.null(mc$dispOutlier)) rep(FALSE, nrow(dds)) else (mc$dispOutlier %in% TRUE)
    dispersion_pas_rows[[length(dispersion_pas_rows) + 1L]] <<- data.frame(
      Comparison = comparison_name, View = label, feature_id = rownames(dds)[keep],
      baseMean = mc$baseMean[keep], dispGeneEst = mc$dispGeneEst[keep], dispFit = mc$dispFit[keep],
      dispersion = DESeq2::dispersions(dds)[keep], dispOutlier = outl[keep], stringsAsFactors = FALSE)
    bm <- mc$baseMean[keep]
    rng <- stats::quantile(bm, c(0.005, 0.995))
    grid_bm <- 10^seq(log10(rng[[1]]), log10(rng[[2]]), length.out = 80)
    dispersion_trend_rows[[length(dispersion_trend_rows) + 1L]] <<- data.frame(
      Comparison = comparison_name, View = label, baseMean = grid_bm,
      trend = DESeq2::dispersionFunction(dds)(grid_bm), stringsAsFactors = FALSE)
  }
  knots <- read_cutoff_data("NormEVS_Regime_Knots.csv")
  if (!is.null(knots) && !is.null(norm_reg$evs)) {
    comp <- evs_regime_composition(norm_reg$evs, knots)
    regime_composition_rows[[length(regime_composition_rows) + 1L]] <<- data.frame(
      Comparison = comparison_name, k = norm_reg$evs$empirical_evs_k,
      LE_n = length(norm_reg$evs$leading_edge_ids), class = names(comp), n = as.integer(comp),
      c1 = knots$c1[1], c2 = knots$c2[1], N_knots = knots$N[1], stringsAsFactors = FALSE)
  }
  invisible(TRUE)
}

# Final (MAP) dispersions: each PAS in the NormEVS Leading Edge / Remainder
# against the same PAS in the Original fit. Size factors are shared, so the
# Wald SE of the LFC scales with sqrt(1/mu + alpha).
final_dispersion_pairs <- function(pas) {
  o <- pas[pas$View == "Original (No EVS)", c("Comparison", "feature_id", "baseMean", "dispersion")]
  names(o)[3:4] <- c("baseMean_Original", "disp_Original")
  s <- pas[pas$View %in% c("NormEVS Lead", "NormEVS Rem"), c("Comparison", "View", "feature_id", "dispersion")]
  names(s)[4] <- "disp_split"
  d <- merge(s, o, by = c("Comparison", "feature_id"))
  d <- d[is.finite(d$disp_split) & d$disp_split > 0 & is.finite(d$disp_Original) & d$disp_Original > 0 &
         is.finite(d$baseMean_Original) & d$baseMean_Original > 0, ]
  d$Set <- ifelse(d$View == "NormEVS Lead", "Leading Edge", "Remainder")
  d$log2_ratio <- log2(d$disp_split / d$disp_Original)
  inv_mu <- 1 / d$baseMean_Original
  d$log_SE_ratio <- 0.5 * log((inv_mu + d$disp_split) / (inv_mu + d$disp_Original))
  d[, c("Comparison", "Set", "feature_id", "baseMean_Original", "disp_Original", "disp_split",
        "log2_ratio", "log_SE_ratio")]
}

final_dispersion_score <- function(pairs) {
  pct <- function(x) 100 * (exp(mean(x)) - 1)
  do.call(rbind, lapply(split(pairs, pairs$Comparison), function(d) {
    le <- d$Set == "Leading Edge"
    data.frame(Comparison = d$Comparison[1], Cutoff_Method = CUTOFF_METHOD_LABEL,
               LE_n = sum(le), REM_n = sum(!le),
               mean_log_SE_ratio = mean(d$log_SE_ratio), SE_change_pct = pct(d$log_SE_ratio),
               SE_change_LE_pct = pct(d$log_SE_ratio[le]), SE_change_REM_pct = pct(d$log_SE_ratio[!le]),
               median_log2_final_dispersion_ratio_LE = stats::median(d$log2_ratio[le]),
               median_log2_final_dispersion_ratio_REM = stats::median(d$log2_ratio[!le]),
               pct_PAS_lower_final_dispersion = 100 * mean(d$log2_ratio < 0),
               stringsAsFactors = FALSE)
  }))
}

split_dispersion_score <- function(met) {
  out <- lapply(unique(met$Comparison), function(cn) {
    le <- met[met$Comparison == cn & met$View == "NormEVS Lead", ]
    re <- met[met$Comparison == cn & met$View == "NormEVS Rem", ]
    o <- met[met$Comparison == cn & met$View == "Original (No EVS)", ]
    if (!nrow(le) || !nrow(re)) return(NULL)
    n <- le$n_PAS + re$n_PAS
    W <- (le$n_PAS * le$prior_width + re$n_PAS * re$prior_width) / n
    V <- (le$n_PAS * le$log_disp_resid_var + re$n_PAS * re$log_disp_resid_var) / n
    wo <- if (nrow(o)) o$prior_width[1] else NA_real_
    vo <- if (nrow(o)) o$log_disp_resid_var[1] else NA_real_
    data.frame(Comparison = cn, Cutoff_Method = CUTOFF_METHOD_LABEL, LE_n = le$n_PAS, REM_n = re$n_PAS,
               w_LE = le$prior_width, w_REM = re$prior_width, W_split = W, w_Original = wo,
               prior_width_reduction_pct = 100 * (1 - W / wo),
               v_LE = le$log_disp_resid_var, v_REM = re$log_disp_resid_var, V_split = V, v_Original = vo,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# ---- figure helpers ---------------------------------------------------------

method_colour <- function() {
  lab <- names(METHOD_SLUGS)[match(CUTOFF_METHOD_SLUG, METHOD_SLUGS)]
  if (length(lab) && !is.na(lab)) unname(METHOD_COLS[[lab]]) else "#1A1A1A"
}
comp_levels <- function() as.character(comparison_table$comparison_name)
comp_lab <- function() stats::setNames(sub("_", " / ", comp_levels()), comp_levels())
log_lab <- function() scales::trans_format("log10", scales::math_format(10^.x))

save_fig_mm <- function(p, stem, w_mm, h_mm) {
  w <- w_mm / 25.4; h <- h_mm / 25.4
  draw <- function() if (inherits(p, "ggplot")) print(p) else { grid::grid.newpage(); grid::grid.draw(p) }
  pdf_path <- file.path(paper_fig_dir, paste0(stem, ".pdf"))
  if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf(pdf_path, width = w, height = h)
  else grDevices::pdf(pdf_path, width = w, height = h, useDingbats = FALSE)
  tryCatch(draw(), finally = grDevices::dev.off())
  png_type <- if (isTRUE(capabilities("cairo"))) "cairo" else getOption("bitmapType")
  grDevices::png(file.path(paper_fig_dir, paste0(stem, ".png")), width = w, height = h,
                 units = "in", res = 600, bg = "white", type = png_type)
  tryCatch(draw(), finally = grDevices::dev.off())
  invisible(TRUE)
}

# rows: list; each element is a ggplot or a list(plots = list(...), widths = c(...))
assemble_rows <- function(rows, heights, title) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    built <- lapply(rows, function(r) if (inherits(r, "ggplot")) r else
      patchwork::wrap_plots(r$plots, nrow = 1, widths = r$widths))
    return(patchwork::wrap_plots(built, ncol = 1, heights = heights) +
      patchwork::plot_annotation(title = title, theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 9, face = "bold", margin = ggplot2::margin(b = 4)))))
  }
  built <- lapply(rows, function(r) if (inherits(r, "ggplot")) ggplot2::ggplotGrob(r) else
    gridExtra::arrangeGrob(grobs = lapply(r$plots, ggplot2::ggplotGrob), nrow = 1, widths = r$widths))
  gridExtra::arrangeGrob(grobs = built, ncol = 1, heights = heights,
                         top = grid::textGrob(title, x = 0.01, hjust = 0,
                                              gp = grid::gpar(fontsize = 9, fontface = "bold")))
}

# ---- panels -----------------------------------------------------------------

panel_scan <- function(scan, weighted, tag) {
  cmp <- comp_levels()
  ks <- data.frame(comparison = cmp, k = vapply(cmp, get_active_evs_cutoff, integer(1)))
  long <- rbind(data.frame(comparison = scan$comparison, k = scan$k, value = scan$G, Q = "G"),
                data.frame(comparison = scan$comparison, k = scan$k, value = scan$R, Q = "R"))
  long$comparison <- factor(long$comparison, levels = cmp)
  ks$comparison <- factor(ks$comparison, levels = cmp)
  q_lab <- c(G = "G: joint + opposite-LE + opposite-Divergence entries", R = "R: opposite-Remainder entries")
  ggplot(long, aes(k, value, colour = Q)) +
    geom_vline(data = ks, aes(xintercept = k), linetype = "22", linewidth = 0.35, colour = "grey25") +
    geom_line(linewidth = 0.6) +
    # G and R both rise from the origin, so the top-left corner is always free
    geom_text(data = ks, aes(x = -Inf, y = Inf, label = paste0("k* = ", scales::comma(k))),
              inherit.aes = FALSE, hjust = -0.12, vjust = 1.4, size = 2.1, fontface = "bold") +
    facet_wrap(~ comparison, nrow = 1, scales = "free", labeller = labeller(comparison = comp_lab())) +
    scale_colour_manual(values = c(G = "#1A1A1A", R = "#9A9A9A"), labels = q_lab) +
    scale_x_continuous(labels = scales::label_comma(), breaks = scales::breaks_pretty(3)) +
    scale_y_continuous(labels = if (weighted) scales::label_percent(accuracy = 1) else scales::label_comma(),
                       expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = "k (top-|PC1| PASs taken from each condition)",
         y = if (weighted) "Share of total weight" else "PASs",
         title = "Candidate cutoffs: what each k adds to the Leading Edge", tag = tag) +
    master_theme() + theme(legend.position = "top", legend.justification = "left",
                           legend.margin = margin(0, 0, -4, 0))
}

panel_frontier <- function(fr, tag) {
  cmp <- comp_levels(); col <- method_colour()
  fr <- fr[order(fr$comparison, fr$remainder_norm, fr$good_norm), ]
  ends <- do.call(rbind, lapply(cmp, function(cn) {
    f <- fr[fr$comparison == cn, ]
    data.frame(comparison = cn, x1 = f$remainder_norm[1], y1 = f$good_norm[1],
               x2 = f$remainder_norm[nrow(f)], y2 = f$good_norm[nrow(f)])
  }))
  sel <- do.call(rbind, lapply(cmp, function(cn) {
    f <- fr[fr$comparison == cn, ]; k <- get_active_evs_cutoff(cn)
    s <- f[f$k == k, , drop = FALSE]
    if (!nrow(s)) s <- f[which.max(f$endpoint_deviation), , drop = FALSE]
    e <- ends[ends$comparison == cn, ]
    dx <- e$x2 - e$x1; dy <- e$y2 - e$y1
    t <- ((s$remainder_norm - e$x1) * dx + (s$good_norm - e$y1) * dy) / (dx^2 + dy^2)
    data.frame(comparison = cn, x = s$remainder_norm[1], y = s$good_norm[1], k = s$k[1],
               fx = e$x1 + t[1] * dx, fy = e$y1 + t[1] * dy, dev = s$endpoint_deviation[1])
  }))
  fr$comparison <- factor(fr$comparison, levels = cmp)
  ends$comparison <- factor(ends$comparison, levels = cmp)
  sel$comparison <- factor(sel$comparison, levels = cmp)
  ggplot(fr, aes(remainder_norm, good_norm)) +
    geom_segment(data = ends, aes(x = x1, y = y1, xend = x2, yend = y2), linetype = "22",
                 linewidth = 0.35, colour = "grey45") +
    geom_segment(data = sel, aes(x = x, y = y, xend = fx, yend = fy), linetype = "11",
                 linewidth = 0.4, colour = "grey20") +
    geom_path(colour = col, linewidth = 0.55) +
    geom_point(colour = col, size = 0.45) +
    geom_point(data = sel, aes(x, y), shape = 21, size = 2.6, stroke = 0.7, fill = col, colour = "black") +
    geom_text(data = sel, aes(x = -Inf, y = Inf, label = sprintf("k* = %s\nd = %.2f", scales::comma(k), dev)),
              hjust = -0.1, vjust = 1.25, size = 2.0, lineheight = 0.9, fontface = "bold") +
    facet_wrap(~ comparison, nrow = 1, labeller = labeller(comparison = comp_lab())) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = TRUE) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = "R (normalized on the Pareto frontier)", y = "G (normalized)",
         title = "Pareto frontier and endpoint chord: k* = maximum deviation d from the chord",
         tag = tag) +
    master_theme() + theme(panel.grid.major.y = element_blank())
}

# Two bars per comparison when a cutoff scan exists:
#   "Cutoff-scan ranking": regime classes at k in the ranking the cutoff was chosen
#                          on (CPM, VST or normalized-count arms, own knots);
#   "Actual NormEVS split": the Leading Edge this run actually fitted with DESeq2
#                          (condition PC1 on the 10-sample comparison), classed
#                          with the normalized-count knots.
# Showing both makes any disagreement between the two rankings visible.
panel_composition <- function(comp, tag, subtitle = NULL, method_by_comparison = NULL, scan_comp = NULL) {
  cmp <- comp_levels()
  act_lab <- "Actual NormEVS split"
  a <- data.frame(Comparison = comp$Comparison, k = comp$k, class = comp$class, n = comp$n,
                  Source = act_lab, stringsAsFactors = FALSE)
  if (!is.null(scan_comp) && nrow(scan_comp)) a <- rbind(a, scan_comp[, names(a)])
  a$class <- factor(a$class, levels = names(REGIME_COLS))
  a$Source <- factor(a$Source, levels = c(act_lab, "Cutoff-scan ranking"))
  a$Comparison <- factor(a$Comparison, levels = cmp)
  tot <- stats::aggregate(n ~ Comparison + Source, data = a, FUN = sum)
  jn <- a[a$class == "Joint", c("Comparison", "Source", "n")]; names(jn)[3] <- "joint"
  tot <- merge(tot, jn, all.x = TRUE)
  tot$lab <- sprintf("%s  (%.0f%% joint)", scales::comma(tot$n), 100 * tot$joint / tot$n)
  kk <- unique(a[a$Source == act_lab, c("Comparison", "k")])
  strip <- stats::setNames(sprintf("%s\nk = %s", comp_lab()[as.character(kk$Comparison)], scales::comma(kk$k)),
                           as.character(kk$Comparison))
  if (!is.null(method_by_comparison)) {
    strip <- stats::setNames(paste0(strip, "\n", method_by_comparison[names(strip)]), names(strip))
  }
  ggplot(a, aes(n, Source, fill = class)) +
    geom_col(width = 0.72, position = position_stack(reverse = TRUE), colour = "white", linewidth = 0.2) +
    geom_text(data = tot, aes(x = n, y = Source, label = lab), inherit.aes = FALSE, hjust = -0.06, size = 1.9) +
    facet_grid(Comparison ~ ., switch = "y", labeller = labeller(Comparison = strip)) +
    scale_fill_manual(values = REGIME_COLS, drop = FALSE) +
    scale_x_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.48))) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(x = "PASs in the Leading Edge", y = NULL, tag = tag,
         title = "What the Leading Edge is made of",
         subtitle = if (is.null(subtitle) && !is.null(scan_comp) && nrow(scan_comp))
           "Cutoff-scan ranking vs the actual NormEVS split fitted by DESeq2" else subtitle) +
    master_theme() + theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
                           panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey92"),
                           strip.placement = "outside", panel.spacing.y = unit(1.2, "mm"),
                           strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6.2,
                                                            face = "bold", lineheight = 0.9),
                           axis.text.y = element_text(size = 5.6, colour = "grey30"),
                           legend.text = element_text(size = 6.2))
}

# Regime counts at the applied k from the cutoff scan of the method that set k.
scan_composition <- function(slug_by_comparison) {
  cache <- list(); out <- list()
  for (cn in names(slug_by_comparison)) {
    sl <- slug_by_comparison[[cn]]
    if (is.null(cache[[sl]])) cache[[sl]] <- read_cutoff_data(paste0(sl, "_scan.csv"))
    sc <- cache[[sl]]
    if (is.null(sc) || !all(c("joint", "opp_LE", "opp_div", "opp_rem") %in% names(sc))) next
    k <- get_active_evs_cutoff(cn)
    r <- sc[sc$comparison == cn & sc$k == k, , drop = FALSE]
    if (!nrow(r)) next
    out[[cn]] <- data.frame(Comparison = cn, k = k, class = names(REGIME_COLS),
                            n = as.integer(c(r$joint[1], r$opp_LE[1], r$opp_div[1], r$opp_rem[1])),
                            Source = "Cutoff-scan ranking", stringsAsFactors = FALSE)
  }
  if (length(out)) do.call(rbind, out) else NULL
}

panel_dispersion <- function(pas, trd, tag, max_pts = 6000L) {
  cmp <- comp_levels()
  d <- pas[pas$View %in% c("NormEVS Lead", "NormEVS Rem") & is.finite(pas$dispersion) &
           pas$dispersion > 0, c("Comparison", "View", "baseMean", "dispersion")]
  set.seed(20260927)
  d <- do.call(rbind, lapply(split(d, list(d$Comparison, d$View), drop = TRUE), function(x)
    if (nrow(x) > max_pts) x[sample.int(nrow(x), max_pts), ] else x))
  d$Set <- factor(ifelse(d$View == "NormEVS Lead", "LE", "REM"), levels = c("REM", "LE"))
  d <- d[order(d$Set), ]
  t <- trd[trd$View %in% c("Original (No EVS)", "NormEVS Lead", "NormEVS Rem"), ]
  t$Line <- factor(c("Original (No EVS)" = "Original trend (no EVS)", "NormEVS Lead" = "Leading Edge trend",
                     "NormEVS Rem" = "Remainder trend")[t$View],
                   levels = c("Leading Edge trend", "Remainder trend", "Original trend (no EVS)"))
  d$Comparison <- factor(d$Comparison, levels = cmp); t$Comparison <- factor(t$Comparison, levels = cmp)
  pts <- if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::geom_point_rast(aes(colour = Set), size = 0.22, alpha = 0.35, stroke = 0,
                             raster.dpi = 600, show.legend = FALSE)
  } else geom_point(aes(colour = Set), size = 0.22, alpha = 0.35, stroke = 0, show.legend = FALSE)
  ggplot(d, aes(baseMean, dispersion)) + pts +
    scale_colour_manual(values = c(LE = unname(SPLIT_POINT_COLS["Leading Edge"]),
                                   REM = unname(SPLIT_POINT_COLS["Remainder"])), guide = "none") +
    trend_layers(t) +
    facet_wrap(~ Comparison, nrow = 1, labeller = labeller(Comparison = comp_lab())) +
    scale_x_log10(labels = log_lab()) + scale_y_log10(labels = log_lab()) +
    annotation_logticks(sides = "bl", short = unit(0.6, "mm"),
                        mid = unit(1, "mm"), long = unit(1.4, "mm"), colour = "grey40") +
    labs(x = "Mean of normalized counts", y = "DESeq2 final (MAP) dispersion", tag = tag,
         title = "DESeq2 final dispersions after the split (NormEVS)",
         subtitle = "Leading Edge (blue) and Remainder (red) were fitted separately; points are the final dispersions used in the Wald test, lines each fit's trend.") +
    master_theme() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
}

# Trend lines drawn with linetype + a fixed dark palette (no second colour scale needed).
trend_layers <- function(t) {
  list(
    geom_line(data = t[t$Line == "Remainder trend", ], aes(baseMean, trend, linetype = Line),
              colour = unname(SPLIT_TREND_COLS["Remainder"]), linewidth = 0.6),
    geom_line(data = t[t$Line == "Leading Edge trend", ], aes(baseMean, trend, linetype = Line),
              colour = unname(SPLIT_TREND_COLS["Leading Edge"]), linewidth = 0.6),
    geom_line(data = t[t$Line == "Original trend (no EVS)", ], aes(baseMean, trend, linetype = Line),
              colour = "black", linewidth = 0.5),
    scale_linetype_manual(values = c("Leading Edge trend" = "solid", "Remainder trend" = "solid",
                                     "Original trend (no EVS)" = "11"), drop = FALSE,
                          guide = guide_legend(override.aes = list(
                            colour = unname(SPLIT_TREND_COLS), linewidth = 0.7)))
  )
}

panel_prior_width <- function(met, score, tag) {
  cmp <- comp_levels(); col <- method_colour()
  m <- met[met$View %in% c("Original (No EVS)", "NormEVS Lead", "NormEVS Rem"), ]
  m$Set <- factor(c("Original (No EVS)" = "Original (no EVS)", "NormEVS Lead" = "Leading Edge",
                    "NormEVS Rem" = "Remainder")[m$View],
                  levels = c("Original (no EVS)", "Leading Edge", "Remainder"))
  m$Comparison <- factor(m$Comparison, levels = cmp)
  s <- score; s$Comparison <- factor(s$Comparison, levels = cmp); s$x <- as.numeric(s$Comparison)
  s$lab <- sprintf("W = %.2f\n%s%.0f%%", s$W_split, ifelse(s$prior_width_reduction_pct >= 0, "-", "+"),
                   abs(s$prior_width_reduction_pct))
  top <- tapply(m$prior_width, m$Comparison, max)
  s$ytxt <- pmax(top[as.character(s$Comparison)], s$W_split)
  ggplot(m, aes(Comparison, prior_width, fill = Set)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72) +
    geom_segment(data = s, aes(x = x - 0.42, xend = x + 0.42, y = W_split, yend = W_split),
                 inherit.aes = FALSE, colour = "white", linewidth = 2.0, lineend = "round") +
    geom_segment(data = s, aes(x = x - 0.4, xend = x + 0.4, y = W_split, yend = W_split),
                 inherit.aes = FALSE, colour = col, linewidth = 1.0, lineend = "round") +
    geom_text(data = s, aes(x = x, y = ytxt, label = lab), inherit.aes = FALSE,
              vjust = -0.35, size = 2.0, lineheight = 0.9, colour = "grey10") +
    scale_fill_manual(values = SPLIT_COLS) +
    scale_x_discrete(labels = comp_lab()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
    labs(x = NULL, y = "Dispersion prior width\n(variance of log dispersion)", tag = tag,
         title = "The split narrows DESeq2's dispersion prior",
         subtitle = "Line: pooled LE + REM prior width W; % narrower than Original.") +
    master_theme() + theme(legend.position = "bottom")
}

panel_final_change <- function(pairs, fscore, tag) {
  cmp <- comp_levels()
  d <- pairs; d$Comparison <- factor(d$Comparison, levels = cmp)
  d$Set <- factor(d$Set, levels = c("Leading Edge", "Remainder"))
  yl <- stats::quantile(d$log2_ratio, c(0.01, 0.99), na.rm = TRUE)
  f <- fscore; f$Comparison <- factor(f$Comparison, levels = cmp)
  f$lab <- sprintf("Wald SE  LE %+.1f%%  |  REM %+.1f%%", f$SE_change_LE_pct, f$SE_change_REM_pct)
  ggplot(d, aes(Comparison, log2_ratio, fill = Set)) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.35, colour = "grey30") +
    geom_violin(position = position_dodge(width = 0.8), width = 0.75, linewidth = 0.2,
                colour = "white", scale = "width", alpha = 0.85) +
    geom_boxplot(aes(group = interaction(Comparison, Set)), position = position_dodge(width = 0.8),
                 width = 0.14, outlier.shape = NA,
                 linewidth = 0.3, colour = "grey10", fill = "white") +
    geom_text(data = f, aes(x = Comparison, y = Inf, label = lab), inherit.aes = FALSE,
              vjust = 1.2, size = 1.9, lineheight = 0.9) +
    coord_cartesian(ylim = c(yl[[1]], yl[[2]] + 0.35 * diff(yl))) +
    scale_fill_manual(values = SPLIT_COLS[c("Leading Edge", "Remainder")]) +
    scale_x_discrete(labels = comp_lab()) +
    labs(x = NULL, y = "log2(final dispersion, split /\nfinal dispersion, Original)", tag = tag,
         title = "Consequence: each PAS's final dispersion moves to its own population's level",
         subtitle = "Leading Edge rises (no longer shrunk toward too low a trend); Remainder falls (sharper tests). Labels: mean change in Wald standard error.") +
    master_theme() + theme(legend.position = "bottom")
}

panel_winner <- function(win_tab, tag) {
  cmp <- comp_levels()
  w <- win_tab
  w$Cutoff_Method <- factor(w$Cutoff_Method, levels = METHOD_LEVELS)
  w$Comparison <- factor(w$Comparison, levels = cmp)
  o <- unique(w[, c("Comparison", "w_Original")])
  ww <- w[w$Winner %in% TRUE, ]
  ggplot(w, aes(Cutoff_Method, W_split)) +
    geom_hline(data = o, aes(yintercept = w_Original), linetype = "22", linewidth = 0.4, colour = "grey25") +
    geom_segment(aes(xend = Cutoff_Method, y = w_Original, yend = W_split, colour = Cutoff_Method),
                 linewidth = 0.6, alpha = 0.6) +
    geom_point(aes(colour = Cutoff_Method), size = 2.3) +
    geom_point(data = ww, shape = 21, size = 4.4, stroke = 0.7, colour = "black", fill = NA) +
    geom_text(data = ww, aes(label = sprintf("selected\nk = %s", scales::comma(k))),
              vjust = 1.7, size = 2.0, lineheight = 0.9, fontface = "bold") +
    scale_x_discrete(expand = expansion(add = 0.8)) +
    facet_wrap(~ Comparison, nrow = 1, scales = "free_y", labeller = labeller(Comparison = comp_lab())) +
    scale_colour_manual(values = METHOD_COLS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.45, 0.12))) +
    labs(x = NULL, y = "Pooled dispersion prior\nwidth W (LE + REM)", tag = tag,
         title = "Dynamic selection: the split that gives DESeq2 the narrowest dispersion prior wins",
         subtitle = "Dashed line: Original (no EVS) prior width. No p-values or discovery counts enter the rule.") +
    master_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

save_method_figure <- function(met, pas, trd, comp, score, pairs, fscore) {
  slug <- CUTOFF_METHOD_SLUG
  stem <- paste0("Figure_Method_", slug)
  comp_row <- NULL
  if (!is.null(comp) && nrow(comp)) comp_row <- comp
  pd <- function(tag) panel_dispersion(pas, trd, tag)
  pe <- function(tag) panel_prior_width(met, score, tag)
  pf <- function(tag) panel_final_change(pairs, fscore, tag)
  slug_by_cmp <- stats::setNames(rep(slug, length(comp_levels())), comp_levels())
  win0 <- read_cutoff_data("Dynamic_Winner.csv")
  if (identical(slug, "Dynamic_Best") && !is.null(win0)) {
    w0 <- win0[win0$Winner %in% TRUE, ]
    slug_by_cmp <- stats::setNames(unname(METHOD_SLUGS[as.character(w0$Cutoff_Method)]), w0$Comparison)
  }
  scan_comp <- tryCatch(scan_composition(slug_by_cmp), error = function(e) NULL)
  pc <- function(tag, sub = NULL, mbc = NULL) panel_composition(comp_row, tag, sub, mbc, scan_comp)

  # Story in every figure: dispersions after the split (trends) -> the prior
  # narrows (headline) and what the Leading Edge is -> consequence for the
  # final dispersions each test uses.
  add_core <- function(rows, hts, mbc = NULL, sub = NULL) {
    L <- function() letters[length(rows) + 1L]
    rows[[length(rows) + 1L]] <- pd(L()); hts <- c(hts, 1.05)
    a <- L(); b <- letters[length(rows) + 2L]
    rows[[length(rows) + 1L]] <- if (!is.null(comp_row))
      list(plots = list(pe(a), pc(b, sub, mbc)), widths = c(1.1, 1)) else pe(a)
    hts <- c(hts, 1.45)
    rows[[length(rows) + 1L]] <- pf(letters[length(rows) + if (!is.null(comp_row)) 2L else 1L]); hts <- c(hts, 1)
    list(rows = rows, hts = hts)
  }
  if (slug %in% c("CPM_Empirical", "VST_Empirical")) {
    scan <- read_cutoff_data(paste0(slug, "_scan.csv"))
    fr <- read_cutoff_data(paste0(slug, "_frontier.csv"))
    rows <- list(); hts <- numeric(0)
    if (!is.null(scan)) { rows[[length(rows) + 1L]] <- panel_scan(scan, FALSE, letters[length(rows) + 1L]); hts <- c(hts, 0.95) }
    if (!is.null(fr)) { rows[[length(rows) + 1L]] <- panel_frontier(fr, letters[length(rows) + 1L]); hts <- c(hts, 0.95) }
    core <- add_core(rows, hts)
    fig <- assemble_rows(core$rows, core$hts, paste0(CUTOFF_METHOD_LABEL, ": cutoff selection and DESeq2 dispersion prior after EVS"))
    save_fig_mm(fig, stem, 180, 240)
  } else if (identical(slug, "Likelihood_Width") && length(likelihood_scan_rows)) {
    core <- add_core(list(panel_loglik("a"), panel_width_vs_k("b")), c(1.05, 0.95))
    fig <- assemble_rows(core$rows, core$hts,
                         "Likelihood k*: cutoff chosen by maximum likelihood (Leading Edge / Remainder dispersion model)")
    save_fig_mm(fig, stem, 180, 240)
  } else if (identical(slug, "Fixed_5000")) {
    core <- add_core(list(), numeric(0), sub = "k is prespecified (no cutoff search).")
    fig <- assemble_rows(core$rows, core$hts, "Fixed 5,000 (prespecified comparator): DESeq2 dispersion prior after EVS")
    save_fig_mm(fig, stem, 180, 205)
  } else {
    win <- read_cutoff_data("Dynamic_Winner.csv")
    rows <- list(); hts <- numeric(0)
    if (!is.null(win)) { rows[[1]] <- panel_winner(win, "a"); hts <- 1 }
    mbc <- NULL
    if (!is.null(win)) {
      ww <- win[win$Winner %in% TRUE, ]
      mbc <- stats::setNames(as.character(ww$Cutoff_Method), ww$Comparison)
    }
    core <- add_core(rows, hts, mbc = mbc)
    fig <- assemble_rows(core$rows, core$hts, "Dynamic best estimator: selection and DESeq2 dispersion prior after EVS")
    save_fig_mm(fig, stem, 180, 240)
  }
  invisible(TRUE)
}

save_dispersion_after_evs_outputs <- function() {
  if (!length(dispersion_metrics_rows)) return(invisible(NULL))
  met <- dplyr::bind_rows(dispersion_metrics_rows)
  trd <- dplyr::bind_rows(dispersion_trend_rows)
  pas <- dplyr::bind_rows(dispersion_pas_rows)
  comp <- if (length(regime_composition_rows)) dplyr::bind_rows(regime_composition_rows) else NULL
  score <- split_dispersion_score(met)
  pairs <- final_dispersion_pairs(pas)
  fscore <- final_dispersion_score(pairs)
  save_csv(pairs, file.path(summary_table_dir, "Table_Final_Dispersion_Change_By_PAS.csv"))
  save_csv(fscore, file.path(summary_table_dir, "Table_Final_Dispersion_Score.csv"))
  save_csv(met, file.path(summary_table_dir, "Table_Dispersion_After_EVS.csv"))
  save_csv(score, file.path(summary_table_dir, "Table_Dispersion_Split_Score.csv"))
  save_csv(trd, file.path(summary_table_dir, "Table_Dispersion_Trends_By_View.csv"))
  save_csv(pas, file.path(summary_table_dir, "Table_Dispersion_Estimates_By_View.csv"))
  if (!is.null(comp)) save_csv(comp, file.path(summary_table_dir, "Table_EVS_Regime_Composition.csv"))
  tryCatch(save_method_figure(met, pas, trd, comp, score, pairs, fscore),
           error = function(e) warning("Method figure skipped: ", conditionMessage(e), call. = FALSE))
  invisible(TRUE)
}

# =============================================================================
# LIKELIHOOD CUTOFF: the k whose Leading Edge / Remainder dispersion model fits best
# =============================================================================
# Varying-dispersion framework (after Tzougas 2020, Risks 8(3):97): dispersion is
# modelled as a function of covariates and the model is chosen by maximum
# likelihood. Here the covariate is EVS membership at cutoff k.
#
#   z_i(k) = 1 if PAS i is in the top k by |PC1 loading| in either condition
#            (the NormEVS Leading Edge), else 0
#   log alpha_i = log alpha_g(mu_i) + eps_i,  eps_i ~ N(0, tau_g^2),  g = z_i(k)
#   log alphahat_i = log alpha_i + e_i,       e_i ~ N(0, s^2),  s^2 = trigamma((m-p)/2)
#
# alpha_g is DESeq2's own trend fitted separately in each group; alphahat_i is
# DESeq2's gene-wise estimate. Gene-wise estimates depend only on each PAS's own
# counts, the design and the shared size factors, so they are computed ONCE per
# comparison; only the two trends change with k. With the variance profiled out
# (robust scale, as DESeq2 itself uses: sigma_g^2 = MAD^2 of the log residuals),
#
#   l(k) = -1/2 * sum_g n_g * [ log(2*pi*sigma_g^2) + 1 ],   sigma_g^2 = tau_g^2 + s^2
#
# k* = argmax l(k). Every k uses the same PASs and the same number of
# parameters, so l(k) is directly comparable across k. The support interval is
# the set of k with l(k) >= l(k*) - 2 (likelihood-ratio interval). No p-value or
# discovery count is used.

LIKELIHOOD_K <- list()
likelihood_scan_rows <- list()
likelihood_selected_rows <- list()
LIKELIHOOD_MIN_GROUP <- 300L

likelihood_cutoff_scan <- function(comparison_name, count_matrix, coldata) {
  message(comparison_name, ": likelihood cutoff scan (dispersion-only DESeq2 fits; no testing).")
  raw <- coerce_raw_count_matrix_for_deseq2(count_matrix, context = paste0(comparison_name, " likelihood cutoff scan"))
  rank_obj <- make_rank_matrix_for_track(raw, coldata, "normalized_evs")
  ids <- colnames(raw)
  trt <- ids[coldata$condition == "trt"]; untrt <- ids[coldata$condition == "untrt"]
  lt <- compute_pc1_loading_table(rank_obj$rank_matrix, trt, top_n = 1L)$loading_table
  lu <- compute_pc1_loading_table(rank_obj$rank_matrix, untrt, top_n = 1L)$loading_table
  depth <- pmin(stats::setNames(lt$rank, lt$feature_id)[rownames(raw)],
                stats::setNames(lu$rank, lu$feature_id)[rownames(raw)])
  names(depth) <- rownames(raw)

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw, colData = coldata, design = make_design_formula(coldata))
  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
  DESeq2::sizeFactors(dds) <- rank_obj$size_factors[colnames(dds)]
  dds <- DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE)
  d <- as.numeric(depth[rownames(dds)])
  X <- stats::model.matrix(DESeq2::design(dds), as.data.frame(SummarizedExperiment::colData(dds)))
  s2 <- trigamma((ncol(dds) - ncol(X)) / 2)

  group_scale <- function(sel) {
    g <- dds[sel, ]
    g <- tryCatch(suppressMessages(DESeq2::estimateDispersionsFit(g, quiet = TRUE)),
                  error = function(e) suppressMessages(DESeq2::estimateDispersionsFit(g, fitType = "local", quiet = TRUE)))
    mc <- S4Vectors::mcols(g)
    ok <- !mc$allZero & is.finite(mc$dispGeneEst) & mc$dispGeneEst >= 1e-6 & is.finite(mc$dispFit) & mc$dispFit > 0
    r <- log(mc$dispGeneEst[ok] / mc$dispFit[ok])
    c(n = sum(ok), sigma2 = stats::mad(r)^2)
  }
  loglik <- function(parts) -0.5 * sum(vapply(parts, function(p) p[["n"]] * (log(2 * pi * p[["sigma2"]]) + 1), numeric(1)))

  orig <- group_scale(rep(TRUE, nrow(dds)))
  l0 <- loglik(list(orig))
  eval_k <- function(k) {
    le <- d <= k
    if (sum(le) < LIKELIHOOD_MIN_GROUP || sum(!le) < LIKELIHOOD_MIN_GROUP) return(NULL)
    a <- group_scale(le); b <- group_scale(!le)
    l <- loglik(list(a, b))
    data.frame(Comparison = comparison_name, k = as.integer(k), LE_n = as.integer(a[["n"]]), REM_n = as.integer(b[["n"]]),
               sigma2_LE = a[["sigma2"]], sigma2_REM = b[["sigma2"]],
               tau2_LE = max(a[["sigma2"]] - s2, 0), tau2_REM = max(b[["sigma2"]] - s2, 0),
               tau2_Original = max(orig[["sigma2"]] - s2, 0), sampling_var = s2,
               loglik = l, delta_loglik_vs_Original = l - l0, stringsAsFactors = FALSE)
  }

  k_max <- min(15000L, as.integer(floor(0.6 * nrow(dds))))
  coarse <- sort(unique(c(seq(250L, k_max, by = 250L), 5000L)))
  rows <- do.call(rbind, lapply(coarse, eval_k))
  kb <- rows$k[which.max(rows$loglik)]
  fine <- setdiff(seq(max(25L, kb - 250L), min(k_max, kb + 250L), by = 25L), rows$k)
  if (length(fine)) rows <- rbind(rows, do.call(rbind, lapply(fine, eval_k)))
  rows <- rows[order(rows$k), ]
  rows$W_split <- (rows$LE_n * rows$tau2_LE + rows$REM_n * rows$tau2_REM) / (rows$LE_n + rows$REM_n)

  best <- rows[which.max(rows$loglik), ]
  sup <- rows$k[rows$loglik >= best$loglik - 2]
  at5k <- rows[rows$k == 5000L, ]
  likelihood_scan_rows[[comparison_name]] <<- rows
  likelihood_selected_rows[[comparison_name]] <<- data.frame(
    Comparison = comparison_name, k_star = best$k, support_lo = min(sup), support_hi = max(sup),
    loglik_k_star = best$loglik, delta_loglik_vs_Original = best$delta_loglik_vs_Original,
    delta_loglik_k_star_vs_5000 = if (nrow(at5k)) best$loglik - at5k$loglik else NA_real_,
    tau2_LE = best$tau2_LE, tau2_REM = best$tau2_REM, tau2_Original = best$tau2_Original,
    W_split = best$W_split, prior_width_reduction_pct = 100 * (1 - best$W_split / best$tau2_Original),
    stringsAsFactors = FALSE)
  message(sprintf("  %s: k* = %d (support %d-%d); delta log-lik vs unsplit = %.1f; vs k = 5,000 = %.1f",
                  comparison_name, best$k, min(sup), max(sup), best$delta_loglik_vs_Original,
                  if (nrow(at5k)) best$loglik - at5k$loglik else NA_real_))
  LIKELIHOOD_K[[comparison_name]] <<- as.integer(best$k)
  invisible(best$k)
}

save_likelihood_tables <- function() {
  if (!length(likelihood_scan_rows)) return(invisible(NULL))
  save_csv(do.call(rbind, likelihood_scan_rows), file.path(summary_table_dir, "Table_Likelihood_Cutoff_Scan.csv"))
  save_csv(do.call(rbind, likelihood_selected_rows), file.path(summary_table_dir, "Table_Likelihood_Cutoff_Selected.csv"))
  invisible(TRUE)
}

# ---- figure panels ------------------------------------------------------------

panel_loglik <- function(tag) {
  cmp <- comp_levels()
  sc <- do.call(rbind, likelihood_scan_rows); sel <- do.call(rbind, likelihood_selected_rows)
  sc$Comparison <- factor(sc$Comparison, levels = cmp); sel$Comparison <- factor(sel$Comparison, levels = cmp)
  col <- method_colour()
  sel$lab <- sprintf("k* = %s\nsupport %s-%s\nlog-lik gain vs 5,000: %+.0f", scales::comma(sel$k_star),
                     scales::comma(sel$support_lo), scales::comma(sel$support_hi), sel$delta_loglik_k_star_vs_5000)
  pk <- sc[match(paste(sel$Comparison, sel$k_star), paste(sc$Comparison, sc$k)), ]
  p5 <- sc[sc$k == 5000L, ]
  ggplot(sc, aes(k, delta_loglik_vs_Original)) +
    geom_rect(data = sel, aes(xmin = support_lo, xmax = support_hi, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = col, alpha = 0.12) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.35, colour = "grey30") +
    geom_vline(xintercept = 5000, linetype = "22", linewidth = 0.35, colour = METHOD_COLS[["Fixed 5,000"]]) +
    geom_line(linewidth = 0.6, colour = "grey15") +
    geom_point(data = p5, colour = METHOD_COLS[["Fixed 5,000"]], size = 1.8) +
    geom_point(data = pk, shape = 21, fill = col, colour = "black", size = 2.6, stroke = 0.6) +
    geom_text(data = sel, aes(x = Inf, y = -Inf, label = lab), inherit.aes = FALSE,
              hjust = 1.04, vjust = -0.25, size = 1.9, lineheight = 0.9) +
    facet_wrap(~ Comparison, nrow = 1, scales = "free_y", labeller = labeller(Comparison = comp_lab())) +
    scale_x_continuous(labels = scales::label_comma(), breaks = scales::breaks_pretty(3)) +
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0.45, 0.08))) +
    labs(x = "k (top-|PC1| PASs taken from each condition)",
         y = "Log-likelihood gain over\nthe unsplit DESeq2 model", tag = tag,
         title = "Choosing k by maximum likelihood: how well two dispersion priors (LE, REM) fit the data",
         subtitle = "Filled point: k*. Shaded: likelihood support interval (within 2 log-lik units). Grey dashed line and point: Fixed 5,000.") +
    master_theme()
}

panel_width_vs_k <- function(tag) {
  cmp <- comp_levels()
  sc <- do.call(rbind, likelihood_scan_rows); sel <- do.call(rbind, likelihood_selected_rows)
  long <- rbind(data.frame(Comparison = sc$Comparison, k = sc$k, w = sc$tau2_LE, Set = "Leading Edge"),
                data.frame(Comparison = sc$Comparison, k = sc$k, w = sc$tau2_REM, Set = "Remainder"),
                data.frame(Comparison = sc$Comparison, k = sc$k, w = sc$W_split, Set = "Pooled W"))
  long$Set <- factor(long$Set, levels = c("Leading Edge", "Remainder", "Pooled W"))
  long$Comparison <- factor(long$Comparison, levels = cmp)
  o <- unique(sc[, c("Comparison", "tau2_Original")]); o$Comparison <- factor(o$Comparison, levels = cmp)
  sel$Comparison <- factor(sel$Comparison, levels = cmp)
  ggplot(long, aes(k, w, colour = Set, linetype = Set)) +
    geom_hline(data = o, aes(yintercept = tau2_Original), linetype = "11", linewidth = 0.45, colour = "black") +
    geom_vline(data = sel, aes(xintercept = k_star), linewidth = 0.4, colour = method_colour()) +
    geom_vline(xintercept = 5000, linetype = "22", linewidth = 0.35, colour = METHOD_COLS[["Fixed 5,000"]]) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ Comparison, nrow = 1, labeller = labeller(Comparison = comp_lab())) +
    scale_colour_manual(values = c("Leading Edge" = unname(SPLIT_COLS["Leading Edge"]),
                                   "Remainder" = unname(SPLIT_COLS["Remainder"]), "Pooled W" = "grey20")) +
    scale_linetype_manual(values = c("Leading Edge" = "solid", "Remainder" = "solid", "Pooled W" = "42")) +
    scale_x_continuous(labels = scales::label_comma(), breaks = scales::breaks_pretty(3)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(x = "k (top-|PC1| PASs taken from each condition)", y = "Dispersion prior width", tag = tag,
         title = "Prior widths across cutoffs: one for the Leading Edge, one for the Remainder",
         subtitle = "Dotted line: unsplit (Original) prior width. Solid vertical: k*. Dashed vertical: Fixed 5,000.") +
    master_theme() + theme(legend.position = "bottom")
}

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  dir.create(file.path(output_dir, comparison_name), recursive = TRUE, showWarnings = FALSE)

  # Likelihood k* is estimated here, on this comparison's own counts, before
  # any split or test (dispersion-only DESeq2 fits; no p-values).
  if (identical(CUTOFF_METHOD_KEY, "ml_width")) {
    likelihood_cutoff_scan(comparison_name, count_matrix, coldata)
  }

  for (track_key in c("normalized_evs", "raw_evs")) {
    message("Running ", comparison_name, " ", unname(track_short[track_key]))
    ok <- run_one_evs_track(
      comparison_name = comparison_name,
      track_key = track_key,
      count_matrix = count_matrix,
      coldata = coldata,
      annot_df = annot_df
    )
    if (!isTRUE(ok)) stop("Analysis track failed: ", comparison_name, " / ", track_key)
  }

  tryCatch(collect_dispersion_after_evs(comparison_name),
           error = function(e) warning("Dispersion diagnostics skipped for ", comparison_name, ": ", conditionMessage(e)))
  export_comparison_manuscript_tables(comparison_name)
  TRUE
}


# -----------------------------------------------------------------------------
# Paper-ready multi-comparison volcano panels
# -----------------------------------------------------------------------------

read_result_table_for_panel <- function(comparison_name, track_key, dataset_key) {
  df <- get_registered_result(
    comparison_name = comparison_name,
    track_key = track_key,
    dataset_key = dataset_key
  )

  if (is.null(df)) {
    warning(
      "Missing registered result for manuscript panel: ",
      comparison_name, " / ", track_key, " / ", dataset_key
    )
    return(NULL)
  }

  df
}


compute_global_volcano_y_limit <- function() {
  views <- comparison_analysis_views()
  running_max <- NA_real_

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (i in seq_len(nrow(views))) {
      df <- read_result_table_for_panel(
        comparison_name = comparison_name,
        track_key = views$track_key[i],
        dataset_key = views$dataset_key[i]
      )
      if (is.null(df)) next

      plot_df <- tryCatch(
        build_final_volcano_df(df, y_col = "neglog10_empirical_p"),
        error = function(e) NULL
      )
      if (is.null(plot_df) || !nrow(plot_df)) next

      this_max <- suppressWarnings(max(plot_df$neglog10_empirical_p, na.rm = TRUE))
      if (is.finite(this_max)) {
        running_max <- if (is.na(running_max)) this_max else max(running_max, this_max)
      }
    }
  }

  if (is.na(running_max)) return(NULL)
  running_max * 1.05
}

save_paper_volcano_panels <- function() {
  dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)
  views <- comparison_analysis_views()
  short_view <- c(
    "Original (No EVS)" = "Orig",
    "NormEVS Lead" = "N-Lead",
    "NormEVS Rem" = "N-Rem",
    "RawEVS Lead" = "R-Lead",
    "RawEVS Rem" = "R-Rem"
  )

  # Computed once so every comparison/view volcano in the manuscript shares
  # the same y-axis scale (see plot_final_volcano's y_limit_override).
  shared_y_limit <- compute_global_volcano_y_limit()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    plots <- list()
    for (i in seq_len(nrow(views))) {
      track_key <- views$track_key[i]
      dataset_key <- views$dataset_key[i]
      df <- read_result_table_for_panel(comparison_name, track_key, dataset_key)
      if (is.null(df)) next
      analysis_label <- analysis_view_label(track_key, dataset_key)
      dataset_name <- paste(comparison_name, unname(track_short[track_key]), dataset_key, sep = "_")

      plots[[length(plots) + 1L]] <- plot_final_volcano(
        df = df,
        dataset_name = dataset_name,
        short_title = unname(short_view[analysis_label]),
        label_genes = TRUE,
        y_limit_override = shared_y_limit,
        n_labels = 10L
      ) +
        theme(
          plot.title = element_text(size = base_theme_size, face = "bold"),
          axis.title = element_text(size = base_theme_size - 1),
          axis.text = element_text(size = base_theme_size - 2),
          plot.caption = element_text(size = base_theme_size - 3)
        )
    }

    plots <- Filter(Negate(is.null), plots)
    if (!length(plots)) next

    panel <- assemble_one_legend_panel(
      plots,
      panel_title = paste0(comparison_name, " | ", CUTOFF_METHOD_LABEL, " | significance across ", length(plots), " analysis view", ifelse(length(plots) == 1L, "", "s")),
      ncol = length(plots)
    )

    panel_dir <- file.path(output_dir, comparison_name, "Panels")
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
    base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_Volcano_", length(plots), "Views"))
    panel_width <- max(10.0, 5.0 * length(plots))
    save_grob(panel, paste0(base, ".png"), width = panel_width, height = 7.6)
    save_grob(panel, paste0(base, ".pdf"), width = panel_width, height = 7.6)
  }

  invisible(TRUE)
}

pca_variance_stats <- function(value_df, coldata) {
  x <- as.matrix(value_df)
  storage.mode(x) <- "numeric"

  common_samples <- intersect(colnames(x), rownames(coldata))
  x <- x[, common_samples, drop = FALSE]
  coldata_use <- coldata[common_samples, , drop = FALSE]

  keep <- rowSums(is.finite(x)) == ncol(x) & apply(x, 1, stats::var) > 0
  x_use <- x[keep, , drop = FALSE]

  if (nrow(x_use) < 2L || ncol(x_use) < 3L) return(NULL)

  # Full PCA is used here so PC1/PC2 explained-variance percentages use the
  # complete non-zero eigenspectrum rather than only the first two components.
  fit <- stats::prcomp(
    t(x_use),
    center = TRUE,
    scale. = FALSE
  )

  eig <- fit$sdev^2
  total_var <- sum(eig)
  if (!is.finite(total_var) || total_var <= 0) return(NULL)

  score_df <- data.frame(
    Sample = rownames(fit$x),
    PC1 = fit$x[, 1],
    PC2 = fit$x[, 2],
    Condition = factor(
      as.character(coldata_use[rownames(fit$x), "condition"]),
      levels = c("untrt", "trt")
    ),
    stringsAsFactors = FALSE
  )

  centroids <- stats::aggregate(cbind(PC1, PC2) ~ Condition, data = score_df, FUN = mean)
  centroid_distance <- NA_real_
  if (nrow(centroids) == 2L) {
    centroid_distance <- sqrt(
      (centroids$PC1[1] - centroids$PC1[2])^2 +
        (centroids$PC2[1] - centroids$PC2[2])^2
    )
  }

  centroid_lookup <- merge(
    score_df,
    centroids,
    by = "Condition",
    suffixes = c("", "_centroid"),
    sort = FALSE
  )
  within_distance <- sqrt(
    (centroid_lookup$PC1 - centroid_lookup$PC1_centroid)^2 +
      (centroid_lookup$PC2 - centroid_lookup$PC2_centroid)^2
  )
  within_group_rms <- sqrt(mean(within_distance^2, na.rm = TRUE))
  separation_ratio <- if (
    is.finite(centroid_distance) && is.finite(within_group_rms) && within_group_rms > 0
  ) {
    centroid_distance / within_group_rms
  } else {
    NA_real_
  }

  list(
    fit = fit,
    coldata = coldata_use,
    score_df = score_df,
    centroids = centroids,
    n_features_input = nrow(x),
    n_features_pca = nrow(x_use),
    pc1_var = eig[1],
    pc2_var = if (length(eig) >= 2L) eig[2] else NA_real_,
    total_var = total_var,
    pc1_fraction = eig[1] / total_var,
    pc2_fraction = if (length(eig) >= 2L) eig[2] / total_var else NA_real_,
    centroid_distance = centroid_distance,
    within_group_rms = within_group_rms,
    separation_ratio = separation_ratio
  )
}

pc1_loading_energy_pct <- function(original_fit, feature_ids) {
  if (is.null(original_fit) || is.null(original_fit$rotation) || !ncol(original_fit$rotation)) {
    return(NA_real_)
  }

  load <- original_fit$rotation[, 1]
  denom <- sum(load^2, na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)

  ids <- intersect(as.character(feature_ids), names(load))
  100 * sum(load[ids]^2, na.rm = TRUE) / denom
}

pca_track_matrix <- function(registry_obj, dataset_key) {
  if (is.null(registry_obj) || is.null(registry_obj$evs)) return(NULL)
  evs <- registry_obj$evs
  mat <- as.data.frame(evs$rank_matrix)

  if (identical(dataset_key, "raw_dataset")) return(mat)
  if (identical(dataset_key, "leading_edge_dataset")) {
    return(mat[evs$leading_edge_ids, , drop = FALSE])
  }
  if (identical(dataset_key, "remainder_dataset")) {
    return(mat[evs$remainder_ids, , drop = FALSE])
  }
  stop("Unknown PCA dataset key: ", dataset_key)
}

compute_pca_support_plot <- function(value_df, coldata, short_title,
                                     original_stats = NULL) {
  stats_obj <- pca_variance_stats(value_df, coldata)
  if (is.null(stats_obj)) return(NULL)

  pca_df <- stats_obj$score_df
  centroid_df <- stats_obj$centroids

  retained <- if (
    !is.null(original_stats) &&
      is.finite(original_stats$pc1_var) &&
      original_stats$pc1_var > 0
  ) {
    100 * stats_obj$pc1_var / original_stats$pc1_var
  } else {
    100
  }

  energy <- if (!is.null(original_stats)) {
    pc1_loading_energy_pct(original_stats$fit, rownames(value_df))
  } else {
    100
  }

  # Segments from each sample to its condition centroid visualize within-group
  # dispersion; centroid crosses summarize treatment/control separation.
  segment_df <- merge(
    pca_df,
    centroid_df,
    by = "Condition",
    suffixes = c("", "_centroid"),
    sort = FALSE
  )

  cap <- paste0(
    "n=", stats_obj$n_features_input,
    " | PC1=", round(100 * stats_obj$pc1_fraction, 1), "%",
    " | λ1=", round(retained, 1), "% Orig",
    " | E1=", round(energy, 1), "%",
    " | Sep=", ifelse(is.finite(stats_obj$separation_ratio),
                       format(round(stats_obj$separation_ratio, 2), trim = TRUE), "NA")
  )

  ggplot(
    pca_df,
    aes(PC1, PC2, label = Sample, shape = Condition, fill = Condition)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey78") +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey78") +
    geom_segment(
      data = segment_df,
      aes(
        x = PC1,
        y = PC2,
        xend = PC1_centroid,
        yend = PC2_centroid,
        color = Condition
      ),
      inherit.aes = FALSE,
      linewidth = 0.32,
      alpha = 0.42,
      show.legend = FALSE
    ) +
    geom_point(size = 2.9, colour = "white", stroke = 0.55) +
    geom_point(
      data = centroid_df,
      aes(PC1, PC2, color = Condition),
      inherit.aes = FALSE,
      shape = 4,
      stroke = 1.15,
      size = 3.5,
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      size = 1.65,
      max.overlaps = 10,
      force = 0.9,
      box.padding = 0.16,
      point.padding = 0.08,
      min.segment.length = 0,
      segment.alpha = 0.45,
      segment.size = 0.16
    ) +
    scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition") +
    scale_fill_manual(values = condition_fills, labels = condition_labels, name = "Condition") +
    scale_color_manual(values = condition_fills, guide = "none") +
    labs(
      title = short_title,
      x = "PC1",
      y = "PC2",
      caption = cap
    ) +
    coord_cartesian(clip = "off") +
    manuscript_theme() +
    theme(
      legend.position = "bottom",
      plot.caption = element_text(size = base_theme_size - 3, hjust = 0.5),
      plot.margin = margin(8, 10, 8, 10)
    )
}

build_evs_split_audit_table <- function() {
  rows <- list()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (track_key in c("normalized_evs", "raw_evs")) {
      obj <- paper_registry[[registry_key(comparison_name, track_key)]]
      if (is.null(obj) || is.null(obj$evs)) next

      evs <- obj$evs
      total_n <- nrow(evs$raw_dataset)
      lead_n <- length(evs$leading_edge_ids)
      rem_n <- length(evs$remainder_ids)

      if (length(evs$trt_top_ids) != evs$empirical_evs_k ||
          length(evs$untrt_top_ids) != evs$empirical_evs_k ||
          lead_n + rem_n != total_n) {
        stop(comparison_name, " ", unname(track_short[track_key]),
             ": empirical EVS split audit failed.")
      }

      rows[[length(rows) + 1L]] <- data.frame(
        Comparison = comparison_name,
        Track = unname(track_short[track_key]),
        Empirical_EVS_k_per_condition = evs$empirical_evs_k,
        TRT_selected_n = length(evs$trt_top_ids),
        UNTRT_selected_n = length(evs$untrt_top_ids),
        Joint_n = length(evs$joint_ids),
        TRT_disjoint_n = length(evs$disjoint_trt_ids),
        UNTRT_disjoint_n = length(evs$disjoint_untrt_ids),
        Leading_edge_union_n = lead_n,
        Remainder_n = rem_n,
        Total_nonzero_PAS_n = total_n,
        Leading_edge_pct = 100 * lead_n / total_n,
        Remainder_pct = 100 * rem_n / total_n,
        TRT_PC1_loading_cutoff = as.numeric(evs$fit_trt$cutoff),
        UNTRT_PC1_loading_cutoff = as.numeric(evs$fit_untrt$cutoff),
        Cutoff_basis = evs$cutoff_basis,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) data.frame() else dplyr::bind_rows(rows)
}


build_evs_pca_evidence_table <- function(comparison_name) {
  rows <- list()

  for (track_key in c("normalized_evs", "raw_evs")) {
    obj <- paper_registry[[registry_key(comparison_name, track_key)]]
    if (is.null(obj) || is.null(obj$evs)) next

    original_mat <- pca_track_matrix(obj, "raw_dataset")
    original_stats <- pca_variance_stats(original_mat, obj$coldata)
    if (is.null(original_stats)) next

    track_rows <- list()

    for (dataset_key in c("raw_dataset", "leading_edge_dataset", "remainder_dataset")) {
      mat <- pca_track_matrix(obj, dataset_key)
      st <- pca_variance_stats(mat, obj$coldata)
      if (is.null(st)) next

      energy <- pc1_loading_energy_pct(original_stats$fit, rownames(mat))
      dataset_label <- c(
        raw_dataset = "Orig",
        leading_edge_dataset = "Lead",
        remainder_dataset = "Rem"
      )[[dataset_key]]

      track_rows[[dataset_key]] <- data.frame(
        Comparison = comparison_name,
        Track = unname(track_short[track_key]),
        Dataset = dataset_label,
        Empirical_EVS_k_per_condition = obj$evs$empirical_evs_k,
        Leading_edge_union_n = length(obj$evs$leading_edge_ids),
        Remainder_n = length(obj$evs$remainder_ids),
        TRT_PC1_loading_cutoff = as.numeric(obj$evs$fit_trt$cutoff),
        UNTRT_PC1_loading_cutoff = as.numeric(obj$evs$fit_untrt$cutoff),
        PAS_n = nrow(mat),
        PCA_PAS_n = st$n_features_pca,
        PC1_explained_pct = 100 * st$pc1_fraction,
        PC2_explained_pct = 100 * st$pc2_fraction,
        PC1_eigenvalue = st$pc1_var,
        PC1_eigenvalue_retained_pct = 100 * st$pc1_var / original_stats$pc1_var,
        Original_PC1_loading_energy_pct = energy,
        Centroid_distance_PC1_PC2 = st$centroid_distance,
        Within_group_RMS_PC1_PC2 = st$within_group_rms,
        Separation_ratio = st$separation_ratio,
        stringsAsFactors = FALSE
      )
    }

    track_tbl <- dplyr::bind_rows(track_rows)

    # Lead and Rem form a partition of the Original feature set. Their shares of
    # Original PC1 squared-loading energy must therefore sum to 100% apart from
    # floating-point tolerance. This assertion catches membership or ID errors.
    lr <- track_tbl[track_tbl$Dataset %in% c("Lead", "Rem"), , drop = FALSE]
    if (nrow(lr) == 2L && all(is.finite(lr$Original_PC1_loading_energy_pct))) {
      energy_sum <- sum(lr$Original_PC1_loading_energy_pct)
      if (abs(energy_sum - 100) > 1e-6) {
        stop(
          comparison_name, " ", unname(track_short[track_key]),
          ": Lead + Rem Original-PC1 loading energy did not sum to 100%."
        )
      }
    }

    rows[[length(rows) + 1L]] <- track_tbl
  }

  if (!length(rows)) data.frame() else dplyr::bind_rows(rows)
}

plot_evs_pca_evidence <- function(comparison_name) {
  tbl <- build_evs_pca_evidence_table(comparison_name)
  if (!nrow(tbl)) return(NULL)

  tbl$Dataset <- factor(tbl$Dataset, levels = c("Orig", "Lead", "Rem"))

  # Three directly interpretable quantities are shown together: PC1 explained
  # variance, PC1 eigenvalue retained relative to Orig, and the fraction of the
  # Original PC1 squared-loading energy contained in each feature subset.
  long <- tbl %>%
    dplyr::select(
      Comparison, Track, Dataset,
      PC1_explained_pct,
      PC1_eigenvalue_retained_pct,
      Original_PC1_loading_energy_pct
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        PC1_explained_pct,
        PC1_eigenvalue_retained_pct,
        Original_PC1_loading_energy_pct
      ),
      names_to = "Metric",
      values_to = "Percent"
    ) %>%
    dplyr::mutate(
      Metric = factor(
        Metric,
        levels = c(
          "PC1_explained_pct",
          "PC1_eigenvalue_retained_pct",
          "Original_PC1_loading_energy_pct"
        ),
        labels = c("PC1 explained", "λ1 vs Orig", "Orig PC1 energy")
      )
    )

  ggplot(long, aes(Dataset, Percent, group = Metric, shape = Metric)) +
    geom_hline(yintercept = 100, linewidth = 0.3, linetype = "dashed", colour = "grey70") +
    geom_line(aes(linetype = Metric), linewidth = 0.55, position = position_dodge(width = 0.08)) +
    geom_point(size = 2.6, position = position_dodge(width = 0.08)) +
    geom_text(
      aes(label = paste0(round(Percent, 1), "%")),
      position = position_dodge(width = 0.08),
      vjust = -0.7,
      size = 2.5,
      check_overlap = TRUE
    ) +
    facet_wrap(~ Track, nrow = 1) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    labs(
      title = paste0(comparison_name, " | quantitative EVS/PC1 evidence"),
      x = NULL,
      y = "Percent",
      shape = NULL,
      linetype = NULL
    ) +
    manuscript_theme() +
    theme(legend.position = "bottom")
}

plot_empirical_hbfss_support <- function(df, short_title) {
  req <- c("empirical_p", "HBFSS", "hc_p_threshold_dataset", "hbfss_threshold_dataset")
  if (length(setdiff(req, names(df))) > 0L) return(NULL)

  plot_df <- df[
    is.finite(df$empirical_p) & !is.na(df$empirical_p) &
      is.finite(df$HBFSS) & !is.na(df$HBFSS),
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(NULL)

  plot_df$empirical_p <- pmax(plot_df$empirical_p, 1e-300)
  method_df <- build_significance_plot_long(plot_df)

  hc_raw <- suppressWarnings(as.numeric(plot_df$hc_p_threshold_dataset[1]))
  htau <- suppressWarnings(as.numeric(plot_df$hbfss_threshold_dataset[1]))

  p <- ggplot() +
    geom_point(
      data = plot_df,
      aes(empirical_p, HBFSS),
      color = plot_palette$background,
      size = 0.45,
      alpha = 0.20
    ) +
    geom_point(
      data = method_df,
      aes(empirical_p, HBFSS, color = Method, shape = Method),
      size = 2.25,
      alpha = 0.98,
      stroke = 0.80
    ) +
    scale_color_manual(
      values = significance_method_colors,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method",
      guide = guide_legend(override.aes = list(size = 3.0, stroke = 0.95))
    ) +
    scale_shape_manual(
      values = significance_method_shapes,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method"
    ) +
    scale_x_log10(labels = scales::label_scientific()) +
    labs(
      title = short_title,
      x = "Empirical p",
      y = "HBFSS",
      caption = paste0(
        if (is.finite(hc_raw) && !is.na(hc_raw)) paste0("HCp=", signif(hc_raw, 3)) else "HCp=NA",
        "  ",
        if (is.finite(htau) && !is.na(htau)) paste0("Hτ=", signif(htau, 3)) else "Hτ=NA"
      )
    ) +
    manuscript_theme() +
    theme(legend.position = "bottom")

  if (is.finite(hc_raw) && !is.na(hc_raw) && hc_raw > 0 && hc_raw <= 1) {
    p <- p + geom_vline(xintercept = hc_raw, color = plot_palette$hc, linewidth = 0.55, linetype = "dotted")
  }
  if (is.finite(htau) && !is.na(htau)) {
    p <- p + geom_hline(yintercept = htau, color = plot_palette$hbfss_line, linewidth = 0.55, linetype = "dashed")
  }
  p
}

save_paper_pca_panels <- function() {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) return(invisible(FALSE))

  all_evidence <- list()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    plots <- list()

    for (track_key in c("normalized_evs", "raw_evs")) {
      obj <- paper_registry[[registry_key(comparison_name, track_key)]]
      if (is.null(obj) || is.null(obj$evs)) next

      original_mat <- pca_track_matrix(obj, "raw_dataset")
      original_stats <- pca_variance_stats(original_mat, obj$coldata)
      if (is.null(original_stats)) next

      track_abbr <- if (track_key == "normalized_evs") "Norm" else "Raw"

      for (dataset_key in c("raw_dataset", "leading_edge_dataset", "remainder_dataset")) {
        mat <- pca_track_matrix(obj, dataset_key)
        ds_abbr <- c(
          raw_dataset = "Orig",
          leading_edge_dataset = "Lead",
          remainder_dataset = "Rem"
        )[[dataset_key]]

        plots[[length(plots) + 1L]] <- compute_pca_support_plot(
          value_df = mat,
          coldata = obj$coldata,
          short_title = paste0(track_abbr, "-", ds_abbr),
          original_stats = original_stats
        )
      }
    }

    plots <- Filter(Negate(is.null), plots)
    panel_dir <- file.path(output_dir, comparison_name, "Panels")
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

    if (length(plots)) {
      panel <- assemble_one_legend_panel(
        plots,
        panel_title = paste0(comparison_name, " | ", CUTOFF_METHOD_LABEL, " | PCA structure before and after EVS"),
        ncol = 3
      )
      base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_PCA_EVS_2x3"))
      save_grob(panel, paste0(base, ".png"), width = 16.5, height = 10.8)
      save_grob(panel, paste0(base, ".pdf"), width = 16.5, height = 10.8)
    }

    evidence_tbl <- build_evs_pca_evidence_table(comparison_name)
    if (nrow(evidence_tbl)) {
      all_evidence[[comparison_name]] <- evidence_tbl
      save_csv(
        evidence_tbl,
        file.path(output_dir, comparison_name, paste0("Table_", comparison_name, "_EVS_PCA_Evidence.csv"))
      )
    }

    pevidence <- plot_evs_pca_evidence(comparison_name)
    if (!is.null(pevidence)) {
      base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_EVS_PC1_Evidence"))
      save_grob(pevidence, paste0(base, ".png"), width = 11.0, height = 5.2)
      save_grob(pevidence, paste0(base, ".pdf"), width = 11.0, height = 5.2)
    }
  }

  evidence_all <- if (length(all_evidence)) dplyr::bind_rows(all_evidence) else data.frame()
  if (nrow(evidence_all)) {
    save_csv(evidence_all, file.path(summary_table_dir, "Table_EVS_PCA_Evidence.csv"))
  }

  invisible(evidence_all)
}

save_paper_empirical_hbfss_panels <- function() {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) {
    return(invisible(FALSE))
  }

  comparison_order <- as.character(comparison_table$comparison_name)

  for (dataset_key in dataset_key_order) {
    plots <- list()

    if (identical(dataset_key, "raw_dataset")) {
      track_order_use <- "normalized_evs"
      panel_title <- paste0("HBFSS calibration | Orig | ", CUTOFF_METHOD_LABEL)
      output_suffix <- "Raw_AllComparisons"
      panel_height <- 5.8
    } else {
      track_order_use <- c("normalized_evs", "raw_evs")
      panel_title <- paste0("HBFSS calibration | ", unname(dataset_short[dataset_key]), " (Norm vs Raw) | ", CUTOFF_METHOD_LABEL)
      output_suffix <- paste0(unname(dataset_short[dataset_key]), "_AllComparisons_NormEVS_vs_RawEVS")
      panel_height <- 11.0
    }

    for (track_key in track_order_use) {
      for (comparison_name in comparison_order) {
        df <- read_result_table_for_panel(
          comparison_name = comparison_name,
          track_key = track_key,
          dataset_key = dataset_key
        )

        if (is.null(df)) {
          next
        }

        short_title <- if (identical(dataset_key, "raw_dataset")) {
          comparison_name
        } else {
          paste0(comparison_name, "\n", unname(track_short[track_key]))
        }

        plots[[length(plots) + 1L]] <- plot_empirical_hbfss_support(
          df = df,
          short_title = short_title
        )
      }
    }

    plots <- Filter(Negate(is.null), plots)

    if (length(plots) == 0L) {
      next
    }

    panel <- assemble_one_legend_panel(
      plots,
      panel_title = panel_title,
      ncol = length(comparison_order)
    )

    save_grob(
      panel,
      file.path(
        paper_fig_dir,
        paste0("Figure_Manuscript_Empirical_HBFSS_", output_suffix, ".png")
      ),
      width = 18.0,
      height = panel_height
    )
  }

  invisible(TRUE)
}

build_discovery_long_table <- function(summary_df) {
  if (!is.data.frame(summary_df) || nrow(summary_df) == 0L) return(data.frame())

  rows <- lapply(seq_len(nrow(summary_df)), function(i) {
    sm <- summary_df[i, , drop = FALSE]
    data.frame(
      Comparison = sm$Comparison,
      Analysis = sm$Analysis,
      Method = factor(
        c("Standard", "Strong", "Weak", "HBFSS"),
        levels = significance_method_levels
      ),
      Count = as.numeric(c(sm$Std, sm$Strong, sm$Weak, sm$HBFSS)),
      Ovlp = as.numeric(sm$Ovlp),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}


save_discovery_count_panel <- function(summary_df) {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) return(invisible(FALSE))
  long_df <- build_discovery_long_table(summary_df)
  if (!nrow(long_df)) return(invisible(FALSE))

  long_df$Method <- factor(long_df$Method, levels = significance_method_levels)

  p <- ggplot(
    long_df,
    aes(Comparison, Count, color = Method, shape = Method)
  ) +
    geom_point(
      position = position_dodge(width = 0.58),
      size = 3.2,
      alpha = 0.98,
      stroke = 0.90
    ) +
    geom_text(
      aes(label = Count),
      position = position_dodge(width = 0.58),
      vjust = -0.65,
      size = 2.7,
      color = "black",
      show.legend = FALSE
    ) +
    facet_wrap(~ Analysis, scales = "free_y", nrow = 1) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    scale_color_manual(
      values = significance_method_colors,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method",
      guide = guide_legend(override.aes = list(size = 3.0, stroke = 0.95))
    ) +
    scale_shape_manual(
      values = significance_method_shapes,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method"
    ) +
    labs(
      title = paste0("Significant PAS counts by method | ", CUTOFF_METHOD_LABEL),
      x = NULL,
      y = "Significant PASs",
      caption = "Ovlp is reported separately; method totals are not mutually exclusive."
    ) +
    manuscript_theme() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(face = "bold")
    )

  save_grob(p, file.path(paper_fig_dir, "Figure_Manuscript_Discovery_Counts.png"), width = 21.0, height = 6.5)
  save_grob(p, file.path(paper_fig_dir, "Figure_Manuscript_Discovery_Counts.pdf"), width = 21.0, height = 6.5)
  invisible(TRUE)
}

save_paper_support_figures <- function(summary_df) {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) {
    return(invisible(FALSE))
  }

  support_steps <- list(
    PCA = function() save_paper_pca_panels(),
    Empirical_HBFSS = function() save_paper_empirical_hbfss_panels(),
    Counts = function() save_discovery_count_panel(summary_df)
  )

  for (nm in names(support_steps)) {
    tryCatch(
      support_steps[[nm]](),
      error = function(e) {
        warning("Support figure generation failed at ", nm, ": ", conditionMessage(e))
      }
    )
  }

  invisible(TRUE)
}


# =============================================================================
# FINAL 3'aTWAS ORTHOLOG COMPARISON
# =============================================================================

ensure_babelgene <- function() {
  if (requireNamespace("babelgene", quietly = TRUE)) {
    return(TRUE)
  }

  message("Installing CRAN package 'babelgene' for database-supported human-to-rat ortholog mapping...")
  suppressWarnings(
    try(
      utils::install.packages(
        "babelgene",
        repos = "https://cloud.r-project.org",
        quiet = TRUE
      ),
      silent = TRUE
    )
  )

  if (!requireNamespace("babelgene", quietly = TRUE)) {
    warning(
      "babelgene could not be installed. TWAS overlap will still run using direct case-insensitive human/rat symbol matches, but database-supported non-identical ortholog symbols cannot be added in this run."
    )
    return(FALSE)
  }

  TRUE
}

read_twas_study <- function(path) {
  twas <- utils::read.csv(
    path,
    skip = 1,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  assert_required_columns(
    twas,
    c("Disease", "PANEL", "Transcript ID", "Gene symbol", "3'aTWAS.Z", "3'aTWAS.P"),
    object_name = "3'aTWAS file"
  )

  twas %>%
    dplyr::transmute(
      Disease = trimws(as.character(Disease)),
      Panel = trimws(as.character(PANEL)),
      TWAS_Transcript = trimws(as.character(`Transcript ID`)),
      Human_TWAS_Gene = trimws(as.character(`Gene symbol`)),
      TWAS_Z = suppressWarnings(as.numeric(`3'aTWAS.Z`)),
      TWAS_P = suppressWarnings(as.numeric(`3'aTWAS.P`)),
      COLOC_PP4 = if ("COLOC.PP4" %in% names(twas)) suppressWarnings(as.numeric(COLOC.PP4)) else NA_real_
    ) %>%
    dplyr::filter(valid_gene_symbol(Human_TWAS_Gene))
}

build_twas_ortholog_map <- function(twas) {
  human_genes <- sort(unique(twas$Human_TWAS_Gene))

  # Direct symbol equivalence is retained for genes whose human symbol differs
  # from the rat WTTS symbol only by capitalization. This preserves obvious
  # one-to-one symbol matches and is combined with database-supported orthology.
  rat_symbols <- sort(unique(OrigID_Symbol$gene_symbol[valid_gene_symbol(OrigID_Symbol$gene_symbol)]))
  rat_lookup <- data.frame(
    gene_key = gene_key(rat_symbols),
    Rat_Ortholog = rat_symbols,
    stringsAsFactors = FALSE
  ) %>% dplyr::filter(!is.na(gene_key)) %>% dplyr::distinct(gene_key, .keep_all = TRUE)

  direct <- data.frame(
    Human_TWAS_Gene = human_genes,
    gene_key = gene_key(human_genes),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::inner_join(rat_lookup, by = "gene_key") %>%
    dplyr::transmute(
      Human_TWAS_Gene,
      Rat_Ortholog,
      Ortholog_support_n = NA_integer_,
      Ortholog_support = "case-insensitive symbol match",
      Mapping_source = "symbol"
    )

  babel <- data.frame()
  if (ensure_babelgene()) {
    orth <- tryCatch(
      babelgene::orthologs(
        genes = human_genes,
        species = TWAS_TARGET_SPECIES,
        human = TRUE,
        min_support = TWAS_ORTHOLOG_MIN_SUPPORT,
        top = FALSE
      ),
      error = function(e) {
        warning("babelgene ortholog mapping failed; continuing with direct symbol matches: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(orth) && nrow(orth)) {
      orth <- as.data.frame(orth, stringsAsFactors = FALSE)
      assert_required_columns(
        orth,
        c("human_symbol", "symbol", "support_n"),
        object_name = "babelgene ortholog output"
      )

      babel <- data.frame(
        Human_TWAS_Gene = as.character(orth$human_symbol),
        Rat_Ortholog = as.character(orth$symbol),
        Ortholog_support_n = suppressWarnings(as.integer(orth$support_n)),
        Ortholog_support = if ("support" %in% names(orth)) as.character(orth$support) else NA_character_,
        Mapping_source = "babelgene",
        stringsAsFactors = FALSE
      ) %>%
        dplyr::filter(
          Human_TWAS_Gene %in% human_genes,
          valid_gene_symbol(Rat_Ortholog)
        )
    }
  }

  dplyr::bind_rows(babel, direct) %>%
    dplyr::mutate(gene_key = gene_key(Rat_Ortholog)) %>%
    dplyr::filter(!is.na(gene_key)) %>%
    dplyr::arrange(Human_TWAS_Gene, dplyr::desc(Ortholog_support_n), Mapping_source) %>%
    dplyr::distinct(Human_TWAS_Gene, Rat_Ortholog, .keep_all = TRUE)
}

build_twas_rat_metadata <- function(twas, mapping) {
  mapped <- twas %>%
    dplyr::left_join(mapping, by = "Human_TWAS_Gene") %>%
    dplyr::filter(!is.na(gene_key))

  mapped %>%
    dplyr::group_by(gene_key) %>%
    dplyr::summarise(
      Rat_Ortholog = dplyr::first(Rat_Ortholog[valid_gene_symbol(Rat_Ortholog)]),
      Human_TWAS_Genes = collapse_unique(Human_TWAS_Gene),
      TWAS_Diseases = collapse_unique(Disease),
      TWAS_n = dplyr::n(),
      TWAS_Panel_n = dplyr::n_distinct(Panel),
      TWAS_Transcript_n = dplyr::n_distinct(TWAS_Transcript),
      TWAS_APA = dplyr::n_distinct(TWAS_Transcript) >= 2L,
      TWAS_min_P = safe_min_numeric(TWAS_P),
      TWAS_max_abs_Z = safe_max_numeric(abs(TWAS_Z)),
      TWAS_max_COLOC_PP4 = safe_max_numeric(COLOC_PP4),
      Ortholog_support_n = safe_max_numeric(Ortholog_support_n),
      Ortholog_mapping_source = collapse_unique(Mapping_source),
      .groups = "drop"
    ) %>%
    dplyr::left_join(WTTS_gene_pas_summary, by = "gene_key")
}

collect_twas_analysis_results <- function() {
  views <- comparison_analysis_views()
  rows <- list()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (i in seq_len(nrow(views))) {
      track_key <- views$track_key[i]
      dataset_key <- views$dataset_key[i]
      df <- get_registered_result(comparison_name, track_key, dataset_key)
      if (is.null(df) || !nrow(df)) next

      pas <- if ("orig_id" %in% names(df)) as.character(df$orig_id) else as.character(df$feature_id)
      bad_pas <- is.na(pas) | !nzchar(trimws(pas))
      pas[bad_pas] <- as.character(df$feature_id[bad_pas])

      rows[[length(rows) + 1L]] <- data.frame(
        Comparison = comparison_name,
        Analysis = analysis_view_label(track_key, dataset_key),
        feature_id = as.character(df$feature_id),
        PAS = pas,
        WTTS_Gene = as.character(df$gene_symbol),
        gene_key = gene_key(df$gene_symbol),
        Apeglm_LFC = suppressWarnings(as.numeric(df$lfc_shrunk)),
        Std_p = suppressWarnings(as.numeric(df$pvalue)),
        Std_BH = suppressWarnings(as.numeric(df$padj)),
        Strong_p = suppressWarnings(as.numeric(df$resGA_pvalue)),
        Strong_BH = suppressWarnings(as.numeric(df$resGA_padj)),
        Weak_p = suppressWarnings(as.numeric(df$resLA_pvalue)),
        Weak_BH = suppressWarnings(as.numeric(df$resLA_padj)),
        EmpP = suppressWarnings(as.numeric(df$empirical_p)),
        HBFSS_score = suppressWarnings(as.numeric(df$HBFSS)),
        HC_alpha0 = HC_ALPHA0,
        HCp = suppressWarnings(as.numeric(df$hc_p_threshold_dataset)),
        Htau = suppressWarnings(as.numeric(df$hbfss_threshold_dataset)),
        Std = !is.na(df$standard_flag) & df$standard_flag,
        Strong = !is.na(df$strong_cnh_flag) & df$strong_cnh_flag,
        Weak = !is.na(df$weak_significant_flag) & df$weak_significant_flag,
        HBFSS = !is.na(df$hbfss_flag) & df$hbfss_flag,
        Ovlp = !is.na(df$any_overlap) & df$any_overlap,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) return(data.frame())
  dplyr::bind_rows(rows) %>% dplyr::filter(!is.na(gene_key))
}

build_twas_overlap_pas_table <- function(all_results, twas_rat_meta) {
  if (!nrow(all_results) || !nrow(twas_rat_meta)) return(data.frame())

  all_results %>%
    dplyr::filter(Std | Strong | Weak | HBFSS) %>%
    dplyr::inner_join(twas_rat_meta, by = "gene_key") %>%
    dplyr::mutate(
      Support = vapply(seq_len(dplyr::n()), function(i) {
        tags <- c(
          if (Std[i]) "Std" else NULL,
          if (Strong[i]) "Str" else NULL,
          if (Weak[i]) "Weak" else NULL,
          if (HBFSS[i]) "HBFSS" else NULL
        )
        paste(tags, collapse = "+")
      }, character(1)),
      Analysis = factor(Analysis, levels = c(
        "Original (No EVS)",
        "NormEVS Lead", "NormEVS Rem",
        "RawEVS Lead", "RawEVS Rem"
      ))
    ) %>%
    dplyr::transmute(
      Comparison,
      Analysis = as.character(Analysis),
      Human_TWAS = Human_TWAS_Genes,
      Rat_Ortholog,
      TWAS_n,
      TWAS_Transcript_n,
      TWAS_APA,
      WTTS_PAS_n,
      WTTS_APA = APA_multi_PAS,
      WTTS_DE_PAS = PAS,
      Apeglm_LFC,
      Std, Strong, Weak, HBFSS, Ovlp,
      Support,
      Std_BH,
      Strong_BH,
      Weak_BH,
      EmpP,
      HBFSS_score
    ) %>%
    dplyr::arrange(
      Comparison,
      factor(Analysis, levels = c(
        "Original (No EVS)",
        "NormEVS Lead", "NormEVS Rem",
        "RawEVS Lead", "RawEVS Rem"
      )),
      Rat_Ortholog,
      WTTS_DE_PAS
    )
}

build_twas_overlap_gene_table <- function(pas_table) {
  if (!nrow(pas_table)) return(data.frame())

  pas_table %>%
    dplyr::mutate(
      Std_flag = Std,
      Strong_flag = Strong,
      Weak_flag = Weak,
      HBFSS_flag = HBFSS,
      Ovlp_flag = Ovlp
    ) %>%
    dplyr::group_by(
      Comparison,
      Analysis,
      Human_TWAS,
      Rat_Ortholog,
      TWAS_n,
      TWAS_Transcript_n,
      TWAS_APA,
      WTTS_PAS_n,
      WTTS_APA
    ) %>%
    dplyr::summarise(
      Std = any(Std_flag, na.rm = TRUE),
      Strong = any(Strong_flag, na.rm = TRUE),
      Weak = any(Weak_flag, na.rm = TRUE),
      HBFSS = any(HBFSS_flag, na.rm = TRUE),
      Ovlp = any(Ovlp_flag, na.rm = TRUE),
      Std_PAS_n = dplyr::n_distinct(WTTS_DE_PAS[Std_flag]),
      Strong_PAS_n = dplyr::n_distinct(WTTS_DE_PAS[Strong_flag]),
      Weak_PAS_n = dplyr::n_distinct(WTTS_DE_PAS[Weak_flag]),
      HBFSS_PAS_n = dplyr::n_distinct(WTTS_DE_PAS[HBFSS_flag]),
      WTTS_DE_PAS_n = dplyr::n_distinct(WTTS_DE_PAS),
      WTTS_DE_APA = dplyr::n_distinct(WTTS_DE_PAS) >= 2L,
      WTTS_DE_PAS_IDs = collapse_unique(WTTS_DE_PAS),
      Methods = paste(
        c(
          if (any(Std_flag, na.rm = TRUE)) "Std" else NULL,
          if (any(Strong_flag, na.rm = TRUE)) "Str" else NULL,
          if (any(Weak_flag, na.rm = TRUE)) "Weak" else NULL,
          if (any(HBFSS_flag, na.rm = TRUE)) "HBFSS" else NULL
        ),
        collapse = "+"
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      Comparison,
      factor(Analysis, levels = c(
        "Original (No EVS)",
        "NormEVS Lead", "NormEVS Rem",
        "RawEVS Lead", "RawEVS Rem"
      )),
      Rat_Ortholog
    )
}

build_twas_overlap_summary <- function(gene_table) {
  if (!nrow(gene_table)) return(data.frame())

  gene_table %>%
    dplyr::group_by(Comparison, Analysis) %>%
    dplyr::summarise(
      Std = sum(Std, na.rm = TRUE),
      Strong = sum(Strong, na.rm = TRUE),
      Weak = sum(Weak, na.rm = TRUE),
      HBFSS = sum(HBFSS, na.rm = TRUE),
      Ovlp = sum(Ovlp, na.rm = TRUE),
      TWAS_APA = sum(TWAS_APA, na.rm = TRUE),
      WTTS_APA = sum(WTTS_APA, na.rm = TRUE),
      WTTS_DE_APA = sum(WTTS_DE_APA, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_twas_gene_support <- function(gene_table, figure_dir) {
  if (!nrow(gene_table)) return(invisible(NULL))

  analysis_map <- c(
    "Original (No EVS)"="Orig",
    "NormEVS Lead"="N-Lead",
    "NormEVS Rem"="N-Rem",
    "RawEVS Lead"="R-Lead",
    "RawEVS Rem"="R-Rem"
  )
  analysis_levels <- unname(analysis_map)

  method_rows <- list()
  for (method in c("Std", "Strong", "Weak", "HBFSS")) {
    keep <- !is.na(gene_table[[method]]) & gene_table[[method]]
    if (!any(keep)) next
    tmp <- gene_table[keep, , drop = FALSE]
    tmp$Method <- c(Std="Standard", Strong="Strong", Weak="Weak", HBFSS="HBFSS")[[method]]
    method_rows[[method]] <- tmp
  }
  if (!length(method_rows)) return(invisible(NULL))

  long <- dplyr::bind_rows(method_rows) %>%
    dplyr::mutate(
      Analysis = factor(unname(analysis_map[Analysis]), levels = analysis_levels),
      Method = factor(Method, levels = significance_method_levels),
      Gene_label = ifelse(
        toupper(Rat_Ortholog) == toupper(Human_TWAS),
        Rat_Ortholog,
        paste0(Rat_Ortholog, " [", Human_TWAS, "]")
      )
    )

  long$Gene_label <- factor(long$Gene_label, levels = rev(sort(unique(long$Gene_label))))

  p <- ggplot(long, aes(Analysis, Gene_label, color = Method, shape = Method)) +
    geom_point(position = position_dodge(width = 0.42), size = 2.7, alpha = 0.98, stroke = 0.85) +
    facet_wrap(~ Comparison, ncol = 2, scales = "free_y") +
    scale_color_manual(
      values = significance_method_colors,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method",
      guide = guide_legend(override.aes = list(size = 3.2, stroke = 0.95))
    ) +
    scale_shape_manual(
      values = significance_method_shapes,
      breaks = significance_method_levels,
      labels = unname(significance_method_labels[significance_method_levels]),
      drop = FALSE,
      name = "Method"
    ) +
    labs(
      title = paste0("3'aTWAS ortholog genes identified in WTTS-Seq | ", CUTOFF_METHOD_LABEL),
      x = NULL,
      y = "Rat ortholog [human TWAS]"
    ) +
    manuscript_theme() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      axis.text.y = element_text(size = 6.8),
      legend.position = "bottom",
      plot.margin = margin(10, 18, 10, 18)
    )

  max_genes <- max(table(long$Comparison))
  height <- min(24, max(8.0, 5.0 + 0.17 * max_genes))
  save_grob(p, file.path(figure_dir, "Figure_TWAS_Gene_Support.png"), width = 15.5, height = height)
  save_grob(p, file.path(figure_dir, "Figure_TWAS_Gene_Support.pdf"), width = 15.5, height = height)
  invisible(p)
}

plot_twas_view_counts <- function(summary_row, title) {
  if (!nrow(summary_row)) return(NULL)
  df <- data.frame(
    Method = factor(c("Standard", "Strong", "Weak", "HBFSS"), levels = significance_method_levels),
    n = c(summary_row$Std[1], summary_row$Strong[1], summary_row$Weak[1], summary_row$HBFSS[1]),
    stringsAsFactors = FALSE
  )

  ggplot(df, aes(Method, n, color = Method, shape = Method)) +
    geom_point(size = 3.2, stroke = 0.90) +
    geom_text(aes(label = n), vjust = -0.8, size = 3.0, show.legend = FALSE) +
    scale_color_manual(values = significance_method_colors, breaks = significance_method_levels,
                       labels = unname(significance_method_labels[significance_method_levels]), name = "Method") +
    scale_shape_manual(values = significance_method_shapes, breaks = significance_method_levels,
                       labels = unname(significance_method_labels[significance_method_levels]), name = "Method") +
    scale_x_discrete(labels = unname(significance_method_labels[significance_method_levels])) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
    labs(
      title = compact_title(title, width = 42),
      subtitle = paste0("Ovlp=", summary_row$Ovlp[1],
                        "  TWAS-APA=", summary_row$TWAS_APA[1],
                        "  WTTS-DE-APA=", summary_row$WTTS_DE_APA[1]),
      x = NULL,
      y = "TWAS genes"
    ) +
    manuscript_theme() +
    theme(legend.position = "none", plot.margin = margin(10, 14, 10, 14))
}

export_twas_by_view <- function(pas_table, gene_table, summary_table) {
  views <- comparison_analysis_views()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (i in seq_len(nrow(views))) {
      track_key <- views$track_key[i]
      dataset_key <- views$dataset_key[i]
      analysis_label <- analysis_view_label(track_key, dataset_key)
      paths <- analysis_view_paths(comparison_name, track_key, dataset_key)

      pas_sub <- if (nrow(pas_table)) {
        pas_table[
          pas_table$Comparison == comparison_name & pas_table$Analysis == analysis_label,
          ,
          drop = FALSE
        ]
      } else pas_table

      gene_sub <- if (nrow(gene_table)) {
        gene_table[
          gene_table$Comparison == comparison_name & gene_table$Analysis == analysis_label,
          ,
          drop = FALSE
        ]
      } else gene_table

      summary_sub <- summary_table[
        summary_table$Comparison == comparison_name & summary_table$Analysis == analysis_label,
        ,
        drop = FALSE
      ]

      save_csv(pas_sub, file.path(paths$tables, "Table_TWAS_PAS.csv"))
      save_csv(gene_sub, file.path(paths$tables, "Table_TWAS_Genes.csv"))

      if (isTRUE(EXPORT_INDIVIDUAL_VIEW_FIGURES)) {
        p <- plot_twas_view_counts(
          summary_sub,
          paste0(comparison_name, " | ", analysis_label, " | 3'aTWAS")
        )
        if (!is.null(p)) {
          save_grob(p, file.path(paths$figures, "TWAS_Overlap.png"), width = 7.2, height = 5.4)
          save_grob(p, file.path(paths$figures, "TWAS_Overlap.pdf"), width = 7.2, height = 5.4)
        }
      }
    }
  }

  invisible(TRUE)
}

save_twas_comparison_panels <- function(summary_table) {
  for (comparison_name in as.character(comparison_table$comparison_name)) {
    sub <- summary_table[summary_table$Comparison == comparison_name, , drop = FALSE]
    if (!nrow(sub)) next

    preferred_views <- c("Original (No EVS)", "NormEVS Lead", "NormEVS Rem", "RawEVS Lead", "RawEVS Rem")
    views <- preferred_views[preferred_views %in% unique(as.character(sub$Analysis))]
    plots <- lapply(views, function(v) {
      one <- sub[sub$Analysis == v, , drop = FALSE]
      if (!nrow(one)) return(NULL)
      short <- c(
        "Original (No EVS)"="Orig",
        "NormEVS Lead"="N-Lead",
        "NormEVS Rem"="N-Rem",
        "RawEVS Lead"="R-Lead",
        "RawEVS Rem"="R-Rem"
      )[[v]]
      plot_twas_view_counts(one, short)
    })
    plots <- Filter(Negate(is.null), plots)
    n_views <- length(plots)
    if (!n_views) next

    panel <- assemble_one_legend_panel(
      lapply(plots, function(p) p + theme(legend.position = "bottom")),
      panel_title = paste0(
        comparison_name, " | ", CUTOFF_METHOD_LABEL, " | 3'aTWAS overlap across ",
        n_views, " view", ifelse(n_views == 1L, "", "s")
      ),
      ncol = n_views
    )
    panel_dir <- file.path(output_dir, comparison_name, "Panels")
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
    base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_TWAS_", n_views, "Views"))
    save_grob(panel, paste0(base, ".png"), width = 20.0, height = 6.0)
    save_grob(panel, paste0(base, ".pdf"), width = 20.0, height = 6.0)
  }
  invisible(TRUE)
}

run_twas_overlap_analysis <- function() {
  twas_file <- resolve_existing_file(twas_file_candidates, "3'aTWAS file")
  message("Running final 3'aTWAS ortholog overlap: ", twas_file)

  twas <- read_twas_study(twas_file)
  mapping <- build_twas_ortholog_map(twas)
  rat_meta <- build_twas_rat_metadata(twas, mapping)
  all_results <- collect_twas_analysis_results()

  pas_table <- build_twas_overlap_pas_table(all_results, rat_meta)
  gene_table <- build_twas_overlap_gene_table(pas_table)
  summary_table <- build_twas_overlap_summary(gene_table)

  twas_dir <- file.path(output_dir, "TWAS")
  twas_fig_dir <- file.path(twas_dir, "figures")
  twas_tab_dir <- file.path(twas_dir, "tables")
  dir.create(twas_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(twas_tab_dir, recursive = TRUE, showWarnings = FALSE)

  save_csv(pas_table, file.path(twas_tab_dir, "Table_TWAS_Overlap_PAS.csv"))
  save_csv(gene_table, file.path(twas_tab_dir, "Table_TWAS_Overlap_Genes.csv"))
  save_csv(summary_table, file.path(twas_tab_dir, "Table_TWAS_Method_Counts.csv"))

  plot_twas_gene_support(gene_table, twas_fig_dir)
  export_twas_by_view(pas_table, gene_table, summary_table)
  save_twas_comparison_panels(summary_table)

  invisible(list(mapping = mapping, pas = pas_table, genes = gene_table, summary = summary_table))
}

comparison_inputs <- lapply(seq_len(nrow(comparison_table)), function(i) {
  prepare_comparison_data(
    comparison_name = comparison_table$comparison_name[i],
    group1_prefix = comparison_table$group1_prefix[i],
    group2_prefix = comparison_table$group2_prefix[i],
    WTTS_Seq = WTTS_Seq,
    meta_all = meta_all
  )
})
names(comparison_inputs) <- comparison_table$comparison_name

failed_comparisons <- list()

for (cmp in names(comparison_inputs)) {
  message("\n=====================================================")
  message("Running comparison: ", cmp)
  message("=====================================================")

  input_obj <- comparison_inputs[[cmp]]
  out <- tryCatch(
    run_full_comparison_pipeline(
      comparison_name = input_obj$comparison_name,
      count_matrix = input_obj$count_matrix,
      coldata = input_obj$coldata,
      annot_df = OrigID_Symbol
    ),
    error = function(e) {
      failed_comparisons[[cmp]] <<- data.frame(
        comparison_name = cmp,
        error_message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      NULL
    }
  )

  if (is.null(out) && is.null(failed_comparisons[[cmp]])) {
    failed_comparisons[[cmp]] <- data.frame(
      comparison_name = cmp,
      error_message = "Comparison did not complete.",
      stringsAsFactors = FALSE
    )
  }
}

if (length(failed_comparisons) > 0L) {
  failed_df <- dplyr::bind_rows(failed_comparisons)
  writeLines(
    apply(failed_df, 1, function(x) paste(x, collapse = " | ")),
    file.path(output_dir, "Failed_Comparisons.txt")
  )
  stop("One or more comparisons failed. See Failed_Comparisons.txt.")
}

overall_summary <- build_overall_manuscript_summary()
if (nrow(overall_summary) > 0L) {
  save_csv(overall_summary, file.path(summary_table_dir, "Table_DE_Method_Counts.csv"))
}

empirical_cutoff_export <- comparison_table[, c("comparison_name", "group1_prefix", "group2_prefix", "fixed_5000_k"), drop = FALSE]
names(empirical_cutoff_export) <- c("Comparison", "RT_prefix", "ZT_prefix", "Fixed_5000_k")
empirical_cutoff_export$Active_Cutoff_Method <- CUTOFF_METHOD_LABEL
empirical_cutoff_export$Active_k_per_condition <- vapply(
  empirical_cutoff_export$Comparison,
  function(cn) tryCatch(get_active_evs_cutoff(cn), error = function(e) NA_integer_), integer(1)
)
empirical_cutoff_export$Basis <- EMPIRICAL_EVS_CUTOFF_BASIS
empirical_cutoff_export$Applied_to <- "NormEVS and RawEVS"
save_csv(
  empirical_cutoff_export,
  file.path(summary_table_dir, "Table_EVS_Cutoffs.csv")
)
# Compatibility filename retained for the existing table-ZIP collector.
save_csv(
  empirical_cutoff_export,
  file.path(summary_table_dir, "Table_EVS_Empirical_Cutoffs.csv")
)

evs_split_audit <- build_evs_split_audit_table()
if (nrow(evs_split_audit) > 0L) {
  save_csv(
    evs_split_audit,
    file.path(summary_table_dir, "Table_EVS_Split_Audit.csv")
  )
}

write_methods_note()
save_paper_volcano_panels()
save_paper_support_figures(overall_summary)
tryCatch(save_likelihood_tables(),
         error = function(e) warning("Likelihood cutoff tables skipped: ", conditionMessage(e)))
tryCatch(save_dispersion_after_evs_outputs(),
         error = function(e) warning("Dispersion-after-EVS outputs skipped: ", conditionMessage(e)))

# Final analytical stage: 3'aTWAS ortholog overlap with every reported WTTS method.
twas_results <- run_twas_overlap_analysis()

figure_zip <- create_all_figures_zip()
table_zip <- create_all_tables_zip()

cat("\n=====================================================\n")
cat("Pipeline complete.\n")
cat("Build: ", PIPELINE_BUILD, "\n", sep = "")
cat("Cutoff method: ", CUTOFF_METHOD_LABEL, "\n", sep = "")
cat("Repository root:\n", repo_root, "\n", sep = "")
cat("Output directory:\n", output_dir, "\n", sep = "")
cat("3'aTWAS overlap: complete\n")
cat("Figure ZIP:\n", figure_zip, "\n", sep = "")
cat("Table ZIP:\n", table_zip, "\n", sep = "")
cat("=====================================================\n\n")

if (nrow(overall_summary) > 0L) print(overall_summary)
