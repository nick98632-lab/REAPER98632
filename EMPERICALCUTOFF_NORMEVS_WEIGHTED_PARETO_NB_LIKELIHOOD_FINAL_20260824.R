#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# CPM-EVS vs DESeq2-EVS — FINAL MANUSCRIPT VERSION
# =============================================================================
#
# DESIGN
# ------
# 1) One experiment-wide PAS universe.
# 2) Two EVS preprocessing methods:
#      CPM-EVS    = log1p(CPM) -> arm-specific PCA
#      DESeq2-EVS = log1p(DESeq2-normalized counts) -> arm-specific PCA
# 3) For each method, all 8 arm-specific D_g(r) curves are fit jointly to
#    estimate one shared c1/c2 regime system.
# 4) Each RT/ZT comparison gets its own independent cutoff k*.
# 5) Cutoff scoring is cross-arm rank-distance weighted and identical across preprocessing methods:
#      - Joint sites receive full unit benefit.
#      - For a disjoint PAS, Delta=|r_selected-r_opposite|/(N-1).
#      - If the opposite-arm rank is in Leading Edge or Divergence, benefit is 1-Delta.
#      - If the opposite-arm rank is in Remainder, contamination cost is Delta.
# 6) k* maximizes equal-weight normalized Pareto utility U=G_norm-R_norm on the
#    comparison-specific Pareto frontier.
# 7) Remainder-crossing sites are NOT silently discarded from the top-k union.
#    They are flagged as penalized/low-confidence. A separate high-confidence
#    subset contains Joint + opposite-LE + opposite-Divergence sites.
# 8) Exactly four composite manuscript figures are generated, at final
#    print size, by the publication figure layer defined below.
#
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT  <- "/root/REAPER98632/exports/cpm_vs_deseq2_evs_manuscript_final"

GROUP_PATTERNS <- c(
  RT0="^R0_", ZT6="^ZT6_", RT2="^R2_", ZT8="^ZT8_",
  RT4="^R4_", ZT10="^ZT10_", RT8="^R8_", ZT14="^ZT14_"
)

COMPARISONS <- list(
  RT0_ZT6=c(control="RT0", treatment="ZT6"),
  RT2_ZT8=c(control="RT2", treatment="ZT8"),
  RT4_ZT10=c(control="RT4", treatment="ZT10"),
  RT8_ZT14=c(control="RT8", treatment="ZT14")
)

METHODS <- c("CPM_EVS", "DESeq2_EVS")

dir.create(OUT_ROOT, recursive=TRUE, showWarnings=FALSE)
FIG_DIR <- file.path(OUT_ROOT, "Figures")
TAB_DIR <- file.path(OUT_ROOT, "Tables")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TAB_DIR, recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# INPUT / NORMALIZATION
# =============================================================================

read_counts <- function(path) {
  if (!file.exists(path)) stop("Count file not found: ", path)

  x <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  idx <- sort(unique(unlist(lapply(GROUP_PATTERNS, \(p) grep(p, names(x))))))

  if (!length(idx)) stop("No sample columns matched GROUP_PATTERNS.")

  ids <- trimws(as.character(x[[1]]))
  blank <- is.na(ids) | ids==""
  ids[blank] <- paste0("feature_", which(blank))
  ids <- make.unique(ids, sep="__dup_")

  m <- do.call(cbind, lapply(x[,idx,drop=FALSE], \(z)
    suppressWarnings(as.numeric(as.character(z)))
  ))

  rownames(m) <- ids
  colnames(m) <- names(x)[idx]
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  m <- pmax(m, 0)

  # One experiment-wide PAS universe.
  m <- m[rowSums(m) > 0, , drop=FALSE]

  if (nrow(m) < 20L) stop("Too few nonzero PASs.")
  m
}

assign_groups <- function(samples) {
  g <- rep(NA_character_, length(samples))
  names(g) <- samples

  for (nm in names(GROUP_PATTERNS)) {
    idx <- grep(GROUP_PATTERNS[[nm]], samples)
    if (any(!is.na(g[idx]))) stop("A sample matched >1 group.")
    g[idx] <- nm
  }

  if (anyNA(g)) {
    stop("Unassigned samples: ", paste(names(g)[is.na(g)], collapse=", "))
  }

  factor(g, levels=names(GROUP_PATTERNS))
}

deseq2_normalize <- function(counts, groups) {
  if (!requireNamespace("DESeq2", quietly=TRUE)) stop("DESeq2 is required.")

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=round(counts),
    colData=data.frame(group=groups, row.names=colnames(counts)),
    design=~group
  )

  dds <- tryCatch(
    DESeq2::estimateSizeFactors(dds),
    error=\(e) DESeq2::estimateSizeFactors(dds, type="poscounts")
  )

  list(
    counts=DESeq2::counts(dds, normalized=TRUE),
    size_factors=DESeq2::sizeFactors(dds)
  )
}

cpm_log1p <- function(x) {
  lib <- colSums(x)
  lib[!is.finite(lib) | lib<=0] <- 1
  log1p(sweep(x, 2, lib/1e6, "/"))
}

pooled_within_group_var <- function(norm, groups) {
  sse <- rep(0, nrow(norm))
  df <- 0L

  for (g in levels(groups)) {
    z <- norm[,groups==g,drop=FALSE]
    mu <- rowMeans(z)
    sse <- sse + rowSums((z-mu)^2)
    df <- df + ncol(z)-1L
  }

  if (df < 2L) stop("Insufficient pooled residual degrees of freedom.")
  pmax(sse/df, 0)
}

# =============================================================================
# PC1 RANKING + CUMULATIVE DIVERGENCE
# =============================================================================

pc1_rank <- function(x) {
  p <- prcomp(t(x), center=TRUE, scale.=FALSE, rank.=1)

  loading <- p$rotation[,1]
  loading[!is.finite(loading)] <- 0

  list(
    order=order(abs(loading), decreasing=FALSE),
    loading=loading,
    contribution=(p$sdev[1]^2)*loading^2
  )
}

build_arm <- function(method, arm, raw, norm, pooled_var) {
  rank_matrix <- switch(
    method,
    CPM_EVS=cpm_log1p(raw),
    DESeq2_EVS=log1p(norm),
    stop("Unknown method: ", method)
  )

  pc <- pc1_rank(rank_matrix)
  ord <- pc$order

  P <- pc$contribution[ord]
  mu <- rowMeans(norm)[ord]
  V <- pooled_var[ord]
  E <- pmax(V-mu, 0)

  if (sum(P)<=0 || sum(E)<=0) {
    stop("Undefined PC1/NB mass for ", method, " / ", arm)
  }

  p <- P/sum(P)
  q <- E/sum(E)

  data.frame(
    method=method,
    arm=arm,
    rank=seq_along(ord),
    feature_id=rownames(raw)[ord],
    abs_pc1_loading=abs(pc$loading[ord]),
    pc1_contribution=P,
    pc1_mass=p,
    excess_variance=E,
    excess_mass=q,
    F_P=cumsum(p),
    F_E=cumsum(q),
    D=cumsum(q)-cumsum(p),
    stringsAsFactors=FALSE
  )
}

# =============================================================================
# SHARED c1/c2 FIT
# =============================================================================

piecewise_basis <- function(x,c1,c2) {
  cbind(1, x, pmax(x-c1,0), pmax(x-c2,0))
}

fit_shared_knots <- function(arms) {
  N <- nrow(arms[[1]])
  if (length(unique(vapply(arms,nrow,integer(1)))) != 1L) {
    stop("All arms must use the same PAS universe.")
  }

  x_full <- (seq_len(N)-1)/(N-1)
  D_full <- do.call(cbind,lapply(arms,\(z) z$D))
  colnames(D_full) <- names(arms)

  # Efficient optimization on a coarse grid; full fit performed once afterward.
  opt_n <- min(4000L,N)
  idx <- unique(as.integer(round(seq(1,N,length.out=opt_n))))
  x <- x_full[idx]
  D <- D_full[idx,,drop=FALSE]

  min_gap <- max(4/(N-1),1e-4)

  objective <- function(par) {
    c1 <- par[1]
    c2 <- par[2]

    if (!is.finite(c1) || !is.finite(c2) ||
        c1<=0 || c2>=1 || c2-c1<=min_gap) return(1e100)

    X <- piecewise_basis(x,c1,c2)
    b <- tryCatch(qr.coef(qr(X),D),error=\(e) NULL)

    if (is.null(b) || any(!is.finite(b))) return(1e100)
    sum((D-X%*%b)^2)
  }

  starts <- list(c(.05,.95),c(.12,.88),c(.22,.78),c(.32,.68))

  fits <- lapply(starts,\(st)
    optim(
      st,objective,
      method="Nelder-Mead",
      control=list(maxit=500,reltol=1e-10)
    )
  )

  best <- fits[[which.min(vapply(fits,`[[`,numeric(1),"value"))]]

  c1 <- as.integer(round(1+best$par[1]*(N-1)))
  c2 <- as.integer(round(1+best$par[2]*(N-1)))

  c1 <- max(2L,min(N-2L,c1))
  c2 <- max(c1+1L,min(N-1L,c2))

  X_full <- piecewise_basis(
    x_full,
    (c1-1)/(N-1),
    (c2-1)/(N-1)
  )

  coef <- qr.coef(qr(X_full),D_full)
  fitted <- X_full%*%coef

  list(
    c1=c1,
    c2=c2,
    N=N,
    groups=names(arms),
    fitted=fitted,
    SSE=sum((D_full-fitted)^2)
  )
}

# =============================================================================
# FAST COUNT-BASED TOP-k SCORING
# =============================================================================

rank_map <- function(df) setNames(df$rank,df$feature_id)

add_events <- function(pos,weight,L) {
  out <- numeric(L)
  ok <- is.finite(pos) & is.finite(weight) & pos>=1 & pos<=L
  if (!any(ok)) return(out)
  z <- rowsum(weight[ok], group=as.integer(pos[ok]), reorder=FALSE)
  out[as.integer(rownames(z))] <- z[,1]
  out
}

activate_from <- function(depth,weight,K) {
  cumsum(add_events(depth,weight,K))
}

active_interval <- function(start,end,weight,K) {
  # Active for start <= k < end.
  delta <- numeric(K+1L)
  ok <- is.finite(start) & is.finite(end) & is.finite(weight) &
        start>=1 & start<=K & end>start
  if (!any(ok)) return(numeric(K))
  ss <- as.integer(start[ok])
  ee <- pmin(as.integer(end[ok]),K+1L)
  w <- weight[ok]
  delta <- delta + add_events(ss,w,K+1L)
  delta <- delta - add_events(ee,w,K+1L)
  cumsum(delta)[seq_len(K)]
}

scan_pair_fast <- function(control_df,treatment_df,c1,c2) {
  N <- nrow(control_df)
  K <- N-c2
  if (K<1L) stop("No candidate k exists beyond c2.")

  rC_map <- rank_map(control_df)
  rT_map <- rank_map(treatment_df)
  ids <- control_df$feature_id
  rC <- as.integer(rC_map[ids])
  rT <- as.integer(rT_map[ids])
  dC <- N-rC+1L
  dT <- N-rT+1L

  # Joint sites receive full unit benefit once both arms contain the PAS.
  joint_depth <- pmax(dC,dT)
  joint_n <- activate_from(joint_depth,rep(1,N),K)

  # Cross-arm rank distance is normalized to [0,1].
  # Supported disjoint PAS: benefit = 1 - delta.
  # Opposite-arm Remainder PAS: contamination = delta.
  rank_delta <- abs(rC-rT)/max(N-1L,1L)

  # Control-only interval.
  idxC <- which(dC<dT & dC<=K)
  startC <- dC[idxC]
  endC <- pmin(dT[idxC],K+1L)
  oppC <- rT[idxC]
  deltaC <- rank_delta[idxC]
  supportC <- ifelse(oppC>=c1, 1-deltaC, 0)
  costC <- ifelse(oppC<c1, deltaC, 0)

  disC_n       <- active_interval(startC,endC,rep(1,length(idxC)),K)
  disC_rem     <- active_interval(startC,endC,as.numeric(oppC<c1),K)
  disC_div     <- active_interval(startC,endC,as.numeric(oppC>=c1 & oppC<=c2),K)
  disC_le      <- active_interval(startC,endC,as.numeric(oppC>c2),K)
  disC_support <- active_interval(startC,endC,supportC,K)
  disC_cost    <- active_interval(startC,endC,costC,K)

  # Treatment-only interval (symmetric rule).
  idxT <- which(dT<dC & dT<=K)
  startT <- dT[idxT]
  endT <- pmin(dC[idxT],K+1L)
  oppT <- rC[idxT]
  deltaT <- rank_delta[idxT]
  supportT <- ifelse(oppT>=c1, 1-deltaT, 0)
  costT <- ifelse(oppT<c1, deltaT, 0)

  disT_n       <- active_interval(startT,endT,rep(1,length(idxT)),K)
  disT_rem     <- active_interval(startT,endT,as.numeric(oppT<c1),K)
  disT_div     <- active_interval(startT,endT,as.numeric(oppT>=c1 & oppT<=c2),K)
  disT_le      <- active_interval(startT,endT,as.numeric(oppT>c2),K)
  disT_support <- active_interval(startT,endT,supportT,K)
  disT_cost    <- active_interval(startT,endT,costT,K)

  disjoint_n <- disC_n+disT_n
  opposite_le_n <- disC_le+disT_le
  opposite_divergence_n <- disC_div+disT_div
  remainder_cross_n <- disC_rem+disT_rem

  # Count-based quantities are retained as diagnostics.
  good_n <- joint_n+opposite_le_n+opposite_divergence_n
  union_n <- 2*seq_len(K)-joint_n
  if (!all(good_n+remainder_cross_n == union_n)) {
    stop("Internal union-count mismatch in distance-weighted top-k scan.")
  }

  # Primary continuous Pareto coordinates.
  supported_disjoint_score <- disC_support+disT_support
  good_score <- joint_n+supported_disjoint_score
  remainder_distance_cost <- disC_cost+disT_cost

  data.frame(
    k=seq_len(K),
    joint_n=joint_n,
    disjoint_n=disjoint_n,
    opposite_le_n=opposite_le_n,
    opposite_divergence_n=opposite_divergence_n,
    remainder_cross_n=remainder_cross_n,
    good_n=good_n,
    union_n=union_n,
    supported_disjoint_score=supported_disjoint_score,
    good_score=good_score,
    remainder_distance_cost=remainder_distance_cost,
    retained_fraction=ifelse(union_n>0,good_n/union_n,NA_real_),
    remainder_cross_fraction=ifelse(union_n>0,remainder_cross_n/union_n,NA_real_),
    stringsAsFactors=FALSE
  )
}

# =============================================================================
# PARETO FRONTIER + EQUAL-WEIGHT NORMALIZED UTILITY
# =============================================================================

mark_pareto_frontier <- function(scan_df, good_col="good_n", cost_col="remainder_cross_n") {
  tmp <- scan_df %>%
    transmute(row_id=row_number(), k=k, good=.data[[good_col]], cost=.data[[cost_col]]) %>%
    arrange(cost,desc(good),desc(k)) %>%
    group_by(cost) %>% slice(1L) %>% ungroup() %>%
    arrange(cost,desc(good))

  running_best_before <- c(-Inf, head(cummax(tmp$good),-1L))
  tmp$is_frontier_coord <- tmp$good > running_best_before
  frontier <- tmp %>% filter(is_frontier_coord) %>% arrange(cost,good,k)

  key_all <- paste(scan_df[[cost_col]],scan_df[[good_col]],sep="::")
  key_frontier <- paste(frontier$cost,frontier$good,sep="::")
  out <- scan_df
  out$is_pareto <- key_all %in% key_frontier
  list(scan=out,frontier=frontier)
}

select_weighted_pareto_optimum <- function(scan_df, good_col="good_n", cost_col="remainder_cross_n",
                                           benefit_weight=1, contamination_weight=1) {
  marked <- mark_pareto_frontier(scan_df,good_col,cost_col)
  frontier <- marked$frontier %>% arrange(cost,good,k)
  if (nrow(frontier)<1L) stop("No Pareto-optimal cutoff points were identified.")

  gr <- range(frontier$good,na.rm=TRUE)
  cr <- range(frontier$cost,na.rm=TRUE)
  norm_good <- function(x) if (diff(gr)==0) rep(1,length(x)) else (x-gr[1])/diff(gr)
  norm_cost <- function(x) if (diff(cr)==0) rep(0,length(x)) else (x-cr[1])/diff(cr)

  frontier$good_norm <- norm_good(frontier$good)
  frontier$remainder_norm <- norm_cost(frontier$cost)
  frontier$weighted_utility <- benefit_weight*frontier$good_norm - contamination_weight*frontier$remainder_norm

  best <- max(frontier$weighted_utility,na.rm=TRUE)
  chosen <- frontier %>%
    filter(abs(weighted_utility-best)<1e-12) %>%
    arrange(desc(good),cost,desc(k)) %>% slice(1L)

  out <- marked$scan
  out$good_norm <- norm_good(out[[good_col]])
  out$remainder_norm <- norm_cost(out[[cost_col]])
  out$weighted_utility <- benefit_weight*out$good_norm - contamination_weight*out$remainder_norm
  out$is_selected_weighted <- out$k==chosen$k[1]

  list(scan=out,frontier=frontier,selected_k=as.integer(chosen$k[1]),
       selected_good=chosen$good[1],selected_cost=chosen$cost[1],
       selected_good_norm=chosen$good_norm[1],selected_remainder_norm=chosen$remainder_norm[1],
       selected_utility=chosen$weighted_utility[1],benefit_weight=benefit_weight,
       contamination_weight=contamination_weight)
}

# =============================================================================
# CLASSIFICATION AT k*
# =============================================================================

classify_at_k <- function(k,comparison,control_df,treatment_df,c1,c2) {
  N <- nrow(control_df)

  rC <- rank_map(control_df)
  rT <- rank_map(treatment_df)

  topC <- tail(control_df$feature_id,k)
  topT <- tail(treatment_df$feature_id,k)
  ids <- union(topC,topT)

  a <- as.integer(rC[ids])
  b <- as.integer(rT[ids])

  inC <- ids%in%topC
  inT <- ids%in%topT
  joint <- inC & inT
  disC <- inC & !inT
  disT <- inT & !inC

  opp <- ifelse(disC,b,ifelse(disT,a,NA_integer_))
  selected_rank <- ifelse(disC,a,ifelse(disT,b,NA_integer_))

  opp_region <- ifelse(
    joint,"Joint",
    ifelse(
      opp>c2,"Opposite Leading Edge",
      ifelse(opp>=c1,"Opposite Divergence","Opposite Remainder")
    )
  )

  # Continuous cross-arm rank concordance. Joint PASs receive full support.
  rank_distance_norm <- ifelse(joint,0,abs(selected_rank-opp)/max(N-1L,1L))
  support <- ifelse(joint,1,ifelse(opp>=c1,1-rank_distance_norm,0))
  penalty <- ifelse(joint,0,ifelse(opp<c1,rank_distance_norm,0))

  data.frame(
    comparison=comparison,
    feature_id=ids,
    selected_k=k,
    control_rank=a,
    treatment_rank=b,
    class=ifelse(joint,"Joint",ifelse(disC,"Control only","Treatment only")),
    opposite_region=opp_region,
    cross_arm_rank_distance=ifelse(joint,0,abs(selected_rank-opp)),
    cross_arm_rank_distance_norm=rank_distance_norm,
    support_score=support,
    remainder_penalty=penalty,
    selected_topk_union=TRUE,
    high_confidence=joint | opp_region%in%c(
      "Opposite Leading Edge","Opposite Divergence"
    ),
    penalized_remainder_crossing=opp_region=="Opposite Remainder",
    stringsAsFactors=FALSE
  )
}

# =============================================================================
# FIGURE SYSTEM (publication layer)
# =============================================================================
#
# Figures are authored at 180 mm final width with 7 pt body text at 100%
# reproduction, in a colourblind-safe palette, and exported as vector PDF plus
# 600 dpi PNG. Figure_Legends.md is generated from the fitted objects so the
# numbers in the legends cannot drift from the numbers in the figures.
#
# Optional: patchwork is used for panel alignment if installed. Without it a
# gtable aligner is used instead and everything still runs.
# =============================================================================

# =============================================================================
# 1. DESIGN SYSTEM
# =============================================================================
# Figures are authored at final print size. A figure drawn 15 in wide and then
# reduced to a 180 mm journal column shrinks all type by ~60%, which is the
# single most common reason draft figures fail production. Everything below is
# sized for 180 mm (double column) reproduced at 100%.

EVS_FIG_W_MM   <- 180    # double-column width
EVS_BASE_PT    <- 7      # body text, in points, at final size
EVS_PNG_DPI    <- 600
EVS_WRITE_TIFF <- FALSE  # set TRUE if the journal requires TIFF

mm2in <- function(mm) mm / 25.4
pt2mm <- function(pt) pt / .pt   # ggplot geom text size is in mm, not points

# Okabe-Ito: colourblind-safe and still separable in greyscale.
EVS_COL <- list(
  cpm       = "#0072B2",  # blue
  deseq     = "#D55E00",  # vermillion
  median    = "#CC79A7",  # reddish purple
  fit       = "#000000",
  arms      = "#9AA3AD",
  knee      = "#CC79A7",
  chord     = "#B9C0C7",
  shared    = "#7A7F86",
  joint     = "#009E73",
  opp_le    = "#56B4E9",
  opp_div   = "#F0E442",
  opp_rem   = "#E69F00",
  band_rem  = "#E8ECF2",  # pale regime fills ...
  band_div  = "#FBF0D0",
  band_lead = "#DCEFE4",
  edge_rem  = "#5B6672",  # ... with saturated partners for rules and labels
  edge_div  = "#8A6D12",
  edge_lead = "#2E8B62"
)

EVS_METHOD_COL <- c("CPM-EVS" = EVS_COL$cpm, "DESeq2-EVS" = EVS_COL$deseq)

evs_method_label <- function(m) ifelse(m == "CPM_EVS", "CPM-EVS", "DESeq2-EVS")

# "RT0_ZT6" -> "RT0 vs ZT6". Underscores do not belong on a published axis.
evs_comparison_label <- function(x) sub("_", " vs ", x, fixed = TRUE)

evs_num <- function(x, digits = 0) {
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

evs_comma <- function(v) {
  format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
}

# legend.position.inside arrived in ggplot2 3.5.0; older versions take a
# numeric legend.position. Keeps the module working on both.
evs_legend_inside <- function(x, y, just, direction = "vertical") {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside",
          legend.position.inside = c(x, y),
          legend.justification = just,
          legend.direction = direction)
  } else {
    theme(legend.position = c(x, y),
          legend.justification = just,
          legend.direction = direction)
  }
}

evs_theme <- function(base_pt = EVS_BASE_PT) {
  theme_classic(base_size = base_pt) +
    theme(
      plot.title        = element_text(size = base_pt + 0.5, face = "bold",
                                       hjust = 0, margin = margin(b = 2.5)),
      plot.tag          = element_text(size = base_pt + 2, face = "bold"),
      plot.tag.position = "topleft",
      axis.title        = element_text(size = base_pt),
      axis.title.x      = element_text(margin = margin(t = 2.5)),
      axis.title.y      = element_text(margin = margin(r = 2.5)),
      axis.text         = element_text(size = base_pt - 0.5, colour = "grey15"),
      axis.line         = element_line(linewidth = 0.25, colour = "grey20"),
      axis.ticks        = element_line(linewidth = 0.25, colour = "grey20"),
      axis.ticks.length = unit(1.1, "pt"),
      legend.position   = "none",
      legend.title      = element_blank(),
      legend.text       = element_text(size = base_pt - 1),
      legend.key.size   = unit(6, "pt"),
      legend.spacing.x  = unit(2, "pt"),
      legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
      legend.margin     = margin(1, 2, 1, 2),
      panel.background  = element_blank(),
      plot.background   = element_blank(),
      plot.margin       = margin(3, 4, 3, 3)
    )
}

evs_theme_inset <- function() {
  theme_classic(base_size = 5) +
    theme(
      axis.title        = element_text(size = 5, colour = "grey20"),
      axis.title.x      = element_text(margin = margin(t = 0.5)),
      axis.title.y      = element_text(margin = margin(r = 0.5)),
      axis.text         = element_text(size = 4.4, colour = "grey25"),
      axis.line         = element_line(linewidth = 0.2, colour = "grey30"),
      axis.ticks        = element_line(linewidth = 0.2, colour = "grey30"),
      axis.ticks.length = unit(0.8, "pt"),
      legend.position   = "none",
      plot.background   = element_rect(fill = "white", colour = "grey75",
                                       linewidth = 0.2),
      plot.margin       = margin(2, 3, 1, 1)
    )
}

# An 8-arm D(r) plot over ~32,000 ranks carries ~256,000 line vertices.
# Thinning is visually identical at print size and keeps the PDF small enough
# for a submission portal.
evs_thin <- function(n, out = 2000L) {
  if (n <= out) return(seq_len(n))
  unique(c(1L, as.integer(round(seq(1, n, length.out = out))), n))
}

# =============================================================================
# 2. DEVICES AND PANEL COMPOSITION
# =============================================================================

evs_open_device <- function(file, width_in, height_in, kind) {
  cairo_ok <- isTRUE(capabilities("cairo"))
  if (kind == "pdf") {
    if (cairo_ok) grDevices::cairo_pdf(file, width = width_in, height = height_in)
    else          grDevices::pdf(file, width = width_in, height = height_in)
  } else if (kind == "png") {
    if (cairo_ok) {
      grDevices::png(file, width = width_in, height = height_in, units = "in",
                     res = EVS_PNG_DPI, bg = "white", type = "cairo")
    } else {
      grDevices::png(file, width = width_in, height = height_in, units = "in",
                     res = EVS_PNG_DPI, bg = "white")
    }
  } else if (kind == "tiff") {
    if (cairo_ok) {
      grDevices::tiff(file, width = width_in, height = height_in, units = "in",
                      res = EVS_PNG_DPI, bg = "white", compression = "lzw",
                      type = "cairo")
    } else {
      grDevices::tiff(file, width = width_in, height = height_in, units = "in",
                      res = EVS_PNG_DPI, bg = "white", compression = "lzw")
    }
  } else {
    stop("Unknown device: ", kind)
  }
}

# Align panel edges across the grid so axes line up column-wise and row-wise.
# Uses gtable only (ships with ggplot2), so no patchwork/cowplot dependency.
evs_align <- function(plots, nrow, ncol) {
  gs <- lapply(plots, ggplot2::ggplotGrob)

  same_w <- length(unique(vapply(gs, function(g) length(g$widths), integer(1)))) == 1L
  same_h <- length(unique(vapply(gs, function(g) length(g$heights), integer(1)))) == 1L

  if (same_w) {
    for (j in seq_len(ncol)) {
      idx <- which(((seq_along(gs) - 1L) %% ncol) + 1L == j)
      if (length(idx) > 1L) {
        w <- do.call(grid::unit.pmax, lapply(gs[idx], function(g) g$widths))
        for (i in idx) gs[[i]]$widths <- w
      }
    }
  }
  if (same_h) {
    for (r in seq_len(nrow)) {
      idx <- which(((seq_along(gs) - 1L) %/% ncol) + 1L == r)
      if (length(idx) > 1L) {
        h <- do.call(grid::unit.pmax, lapply(gs[idx], function(g) g$heights))
        for (i in idx) gs[[i]]$heights <- h
      }
    }
  }
  if (!same_w || !same_h) {
    message("  note: heterogeneous panel layouts; drawing without full alignment")
  }
  gs
}

evs_draw_grid <- function(grobs, nrow, ncol) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow, ncol)))
  for (i in seq_along(grobs)) {
    r  <- ((i - 1L) %/% ncol) + 1L
    cc <- ((i - 1L) %% ncol) + 1L
    pushViewport(viewport(layout.pos.row = r, layout.pos.col = cc))
    grid.draw(grobs[[i]])
    popViewport()
  }
  popViewport()
}

evs_save <- function(plots, file_base, nrow, ncol,
                     width_mm = EVS_FIG_W_MM, height_mm = 150) {
  # patchwork aligns panels correctly even when some are faceted and some are
  # not (Figure 4C is faceted, A/B/D are not). The gtable path below only
  # aligns gtables of identical structure, so it is the fallback.
  use_pw <- requireNamespace("patchwork", quietly = TRUE)
  grobs  <- if (use_pw) NULL else evs_align(plots, nrow, ncol)

  w <- mm2in(width_mm)
  h <- mm2in(height_mm)
  kinds <- c("pdf", "png", if (isTRUE(EVS_WRITE_TIFF)) "tiff")

  for (k in kinds) {
    f <- paste0(file_base, ".", k)
    evs_open_device(f, w, h, k)
    ok <- tryCatch({
      if (use_pw) {
        print(patchwork::wrap_plots(plots, nrow = nrow, ncol = ncol))
      } else {
        evs_draw_grid(grobs, nrow, ncol)
      }
      TRUE
    }, error = function(e) {
      message("  ERROR drawing ", basename(f), ": ", conditionMessage(e))
      FALSE
    })
    grDevices::dev.off()
    if (ok) message("  wrote ", basename(f))
  }
  invisible(file_base)
}

# =============================================================================
# 3. FIGURE 1 PANELS — regime definition
# =============================================================================

evs_arm_matrix <- function(arms) {
  m <- do.call(cbind, lapply(arms, function(z) z$D))
  if (is.null(colnames(m))) colnames(m) <- paste0("arm", seq_len(ncol(m)))
  m
}

panel_regime <- function(method, arms, knot, tag, show_legend = FALSE) {
  N   <- knot$N
  c1  <- knot$c1
  c2  <- knot$c2
  idx <- evs_thin(N, 1800L)

  Dm <- evs_arm_matrix(arms)[idx, , drop = FALSE]
  rr <- idx

  n_arms    <- ncol(Dm)
  lab_arms  <- sprintf("Individual arms (n = %d)", n_arms)
  lab_med   <- "Median observed"
  lab_fit   <- "Shared two-knot fit"
  key_levels <- c(lab_arms, lab_med, lab_fit)

  arms_long <- data.frame(
    rank  = rep(rr, times = n_arms),
    D     = as.vector(Dm),
    arm   = rep(colnames(Dm), each = length(rr)),
    key   = lab_arms,
    stringsAsFactors = FALSE
  )
  med <- data.frame(rank = rr, D = apply(Dm, 1, median), key = lab_med)
  fit <- data.frame(rank = rr,
                    D = apply(knot$fitted[idx, , drop = FALSE], 1, median),
                    key = lab_fit)

  yr   <- range(c(arms_long$D, fit$D), na.rm = TRUE)
  ypad <- diff(yr)
  if (!is.finite(ypad) || ypad <= 0) ypad <- 1
  ylo <- yr[1] - 0.04 * ypad
  yhi <- yr[2] + 0.32 * ypad          # headroom for the regime labels

  bands <- data.frame(
    xmin = c(1, c1, c2),
    xmax = c(c1, c2, N),
    fill = c(EVS_COL$band_rem, EVS_COL$band_div, EVS_COL$band_lead),
    lab  = c("Remainder", "Divergence", "Leading edge"),
    col  = c(EVS_COL$edge_rem, EVS_COL$edge_div, EVS_COL$edge_lead),
    stringsAsFactors = FALSE
  )
  bands$xmid <- (bands$xmin + bands$xmax) / 2
  bands$y    <- yr[2] + 0.265 * ypad

  knots <- data.frame(
    x   = c(c1, c2),
    y   = yr[2] + 0.135 * ypad,
    col = c(EVS_COL$edge_rem, EVS_COL$edge_lead),
    hj  = c(1.06, -0.06),
    lab = c(sprintf("italic(c)[1] * \" = \" * \"%s\"", evs_num(c1)),
            sprintf("italic(c)[2] * \" = \" * \"%s\"", evs_num(c2))),
    stringsAsFactors = FALSE
  )

  p <- ggplot() +
    annotate("rect", xmin = bands$xmin, xmax = bands$xmax,
             ymin = -Inf, ymax = Inf, fill = bands$fill) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dotted",
               linewidth = 0.2) +
    geom_line(data = arms_long,
              aes(rank, D, group = arm, colour = key, linetype = key,
                  linewidth = key)) +
    geom_line(data = med,
              aes(rank, D, colour = key, linetype = key, linewidth = key)) +
    geom_line(data = fit,
              aes(rank, D, colour = key, linetype = key, linewidth = key)) +
    geom_vline(xintercept = c1, colour = EVS_COL$edge_rem,
               linetype = "22", linewidth = 0.3) +
    geom_vline(xintercept = c2, colour = EVS_COL$edge_lead,
               linetype = "22", linewidth = 0.3) +
    geom_text(data = bands, aes(x = xmid, y = y, label = lab),
              colour = bands$col, size = pt2mm(5.8), vjust = 1,
              inherit.aes = FALSE) +
    geom_text(data = knots, aes(x = x, y = y, label = lab),
              colour = knots$col, size = pt2mm(5.6), hjust = knots$hj,
              vjust = 1, fontface = "bold", parse = TRUE, inherit.aes = FALSE) +
    scale_colour_manual(breaks = key_levels,
                        values = setNames(c(EVS_COL$arms, EVS_COL$median,
                                            EVS_COL$fit), key_levels)) +
    scale_linetype_manual(breaks = key_levels,
                          values = setNames(c("solid", "solid", "22"), key_levels)) +
    scale_linewidth_manual(breaks = key_levels,
                           values = setNames(c(0.15, 0.55, 0.42), key_levels)) +
    scale_x_continuous(labels = evs_comma) +
    coord_cartesian(xlim = c(1, N), ylim = c(ylo, yhi), expand = FALSE) +
    labs(
      tag   = tag,
      title = evs_method_label(method),
      x     = "PAS rank by |PC1 loading| (ascending)",
      y     = expression(italic(D)(italic(r)) ==
                           italic(F)[E](italic(r)) - italic(F)[P](italic(r)))
    ) +
    evs_theme()

  if (show_legend) {
    p <- p +
      guides(colour    = guide_legend(order = 1),
             linetype  = guide_legend(order = 1),
             linewidth = guide_legend(order = 1)) +
      evs_legend_inside(0.02, 0.60, c(0, 1))
  }
  p
}

panel_residuals <- function(method_objects, tag, n_bins = 180L) {
  res <- bind_rows(lapply(names(method_objects), function(m) {
    ob <- method_objects[[m]]
    D  <- evs_arm_matrix(ob$arms)
    R  <- D - ob$knot$fitted
    N  <- nrow(R)
    bin <- cut(seq_len(N), breaks = n_bins, labels = FALSE)
    data.frame(
      rank   = rep(seq_len(N), times = ncol(R)),
      bin    = rep(bin, times = ncol(R)),
      resid  = as.vector(R),
      method = evs_method_label(m),
      stringsAsFactors = FALSE
    ) %>%
      group_by(method, bin) %>%
      summarise(rank = mean(rank),
                lo   = stats::quantile(resid, 0.25, names = FALSE),
                mid  = stats::median(resid),
                hi   = stats::quantile(resid, 0.75, names = FALSE),
                .groups = "drop")
  }))
  res$method <- factor(res$method, levels = names(EVS_METHOD_COL))

  ggplot(res, aes(rank, mid, colour = method, fill = method)) +
    geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.25) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, colour = NA) +
    geom_line(aes(linetype = method), linewidth = 0.4) +
    scale_colour_manual(values = EVS_METHOD_COL) +
    scale_fill_manual(values = EVS_METHOD_COL) +
    scale_linetype_manual(values = c("CPM-EVS" = "solid", "DESeq2-EVS" = "22")) +
    scale_x_continuous(labels = evs_comma) +
    labs(
      tag   = tag,
      title = "Shared-knot fit residuals",
      x     = "PAS rank by |PC1 loading| (ascending)",
      y     = expression(paste("Residual, ",
                               italic(D)[italic(g)](italic(r)) -
                                 hat(italic(D))[italic(g)](italic(r))))
    ) +
    evs_theme() +
    evs_legend_inside(0.98, 0.98, c(1, 1), direction = "horizontal")
}

panel_regime_widths <- function(method_objects, tag) {
  d <- bind_rows(lapply(seq_along(method_objects), function(i) {
    m <- names(method_objects)[i]
    k <- method_objects[[m]]$knot
    data.frame(
      method = evs_method_label(m),
      y      = length(method_objects) - i + 1,
      xmin   = c(1, k$c1, k$c2),
      xmax   = c(k$c1, k$c2, k$N),
      regime = c("Remainder", "Divergence", "Leading edge"),
      n      = c(k$c1 - 1L, k$c2 - k$c1 + 1L, k$N - k$c2),
      c1     = k$c1, c2 = k$c2, N = k$N,
      stringsAsFactors = FALSE
    )
  }))
  d$regime <- factor(d$regime, levels = c("Remainder", "Divergence", "Leading edge"))
  d$xmid   <- (d$xmin + d$xmax) / 2

  kn <- d %>%
    distinct(method, y, c1, c2) %>%
    tidyr::pivot_longer(c("c1", "c2"), names_to = "which", values_to = "x") %>%
    mutate(
      lab = ifelse(which == "c1",
                   sprintf("italic(c)[1] * \" = \" * \"%s\"", evs_num(x)),
                   sprintf("italic(c)[2] * \" = \" * \"%s\"", evs_num(x))),
      col = ifelse(which == "c1", EVS_COL$edge_rem, EVS_COL$edge_lead)
    )

  Nmax <- max(d$N)
  h    <- 0.24
  lab_df <- distinct(d, method, y)

  ggplot(d) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - h, ymax = y + h,
                  fill = regime), colour = "grey45", linewidth = 0.2) +
    geom_text(aes(x = xmid, y = y, label = evs_comma(n)), size = pt2mm(5.6)) +
    geom_text(data = kn, aes(x = x, y = y + h + 0.09, label = lab),
              colour = kn$col, size = pt2mm(5.2), vjust = 0,
              parse = TRUE, inherit.aes = FALSE) +
    geom_text(data = lab_df, aes(x = -0.015 * Nmax, y = y, label = method),
              hjust = 1, size = pt2mm(6.2), inherit.aes = FALSE) +
    scale_fill_manual(values = c("Remainder"    = EVS_COL$band_rem,
                                 "Divergence"   = EVS_COL$band_div,
                                 "Leading edge" = EVS_COL$band_lead)) +
    scale_x_continuous(labels = evs_comma, breaks = scales::pretty_breaks(4)) +
    coord_cartesian(xlim = c(-0.40 * Nmax, Nmax),
                    ylim = c(0.30, max(d$y) + 0.58),
                    expand = FALSE, clip = "off") +
    labs(tag = tag, title = "Regime widths (PAS counts)",
         x = "PAS rank by |PC1 loading| (ascending)", y = NULL) +
    evs_theme() +
    theme(
      axis.line.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.text.y       = element_blank(),
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.background = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 1))
}

# =============================================================================
# 4. FIGURES 2-3 PANELS — Pareto frontier and weighted-utility selection
# =============================================================================

panel_pareto <- function(method, comparison, scan, frontier, kstar, tag,
                         method_col) {
  fr <- frontier[order(frontier$cost,frontier$good), , drop = FALSE]
  sel <- fr[fr$k==kstar,,drop=FALSE]
  if (nrow(sel)!=1L) stop("k* not found on the frontier for ", comparison)

  fr_plot <- fr[evs_thin(nrow(fr),1500L),,drop=FALSE]
  xr <- range(fr$cost,na.rm=TRUE); yr <- range(fr$good,na.rm=TRUE)
  xd <- if (diff(xr)>0) diff(xr) else 1; yd <- if (diff(yr)>0) diff(yr) else 1

  ins <- ggplot(fr_plot,aes(k,weighted_utility)) +
    geom_line(colour="grey25",linewidth=0.25) +
    geom_vline(xintercept=sel$k,colour=EVS_COL$knee,linetype="22",linewidth=0.3) +
    labs(x=expression(italic(k)),y="U(k)") + evs_theme_inset()

  ggplot() +
    geom_line(data=fr_plot,aes(cost,good),colour=method_col,linewidth=0.5) +
    annotation_custom(ggplotGrob(ins),
      xmin=xr[1]+0.50*(xr[2]-xr[1]+1e-9),xmax=xr[2]+0.05*xd,
      ymin=yr[1]+0.04*(yr[2]-yr[1]+1e-9),ymax=yr[1]+0.46*(yr[2]-yr[1]+1e-9)) +
    geom_point(data=sel,aes(cost,good),shape=23,size=1.5,stroke=0.3,
               fill=EVS_COL$knee,colour="white") +
    annotate("text",x=sel$cost-0.015*xd,y=sel$good+0.035*yd,
             label=sprintf("italic(k)^"*" * " = " * "%s"",evs_num(sel$k)),
             parse=TRUE,hjust=1,vjust=0,size=pt2mm(6.2),
             colour=EVS_COL$knee,fontface="bold") +
    scale_x_continuous(labels=evs_comma) + scale_y_continuous(labels=evs_comma) +
    labs(tag=tag,title=evs_comparison_label(comparison),
         x=expression(paste("Distance-weighted Remainder disagreement, ",italic(R)[Delta](k))),
         y=expression(paste("Concordance-weighted retained support, ",italic(G)[Delta](k)))) +
    evs_theme()
}

# =============================================================================
# 5. FIGURE 4 PANELS — cross-method summary
# =============================================================================

evs_prep_summary <- function(summary_df, comparisons) {
  summary_df %>%
    mutate(
      method     = factor(evs_method_label(method), levels = names(EVS_METHOD_COL)),
      comparison = factor(evs_comparison_label(comparison),
                          levels = evs_comparison_label(names(comparisons)))
    )
}

panel_kstar <- function(summary_df, comparisons, tag) {
  d <- evs_prep_summary(summary_df, comparisons)
  ggplot(d, aes(comparison, selected_k, fill = method)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62) +
    geom_text(aes(label = evs_comma(selected_k)),
              position = position_dodge(width = 0.72),
              vjust = -0.45, size = pt2mm(5.4)) +
    scale_fill_manual(values = EVS_METHOD_COL) +
    scale_y_continuous(labels = evs_comma,
                       expand = expansion(mult = c(0, 0.16))) +
    labs(tag = tag, title = "Dataset-specific cutoffs", x = NULL,
         y = expression(paste("Selected cutoff, ", italic(k)^"*", " (PAS)"))) +
    evs_theme() +
    evs_legend_inside(0.99, 0.99, c(1, 1), direction = "horizontal")
}

panel_hc_membership <- function(overlap_df, comparisons, tag) {
  d <- overlap_df %>%
    mutate(
      cpm_only   = CPM_high_confidence    - overlap_high_confidence,
      shared     = overlap_high_confidence,
      deseq_only = DESeq2_high_confidence - overlap_high_confidence,
      comparison = factor(evs_comparison_label(comparison),
                          levels = evs_comparison_label(names(comparisons)))
    )

  tot <- d %>%
    transmute(comparison,
              total   = cpm_only + shared + deseq_only,
              jaccard = jaccard_high_confidence)

  long <- d %>%
    select(comparison, cpm_only, shared, deseq_only) %>%
    tidyr::pivot_longer(-comparison, names_to = "set", values_to = "n") %>%
    mutate(set = factor(recode(set,
                               cpm_only   = "CPM-EVS only",
                               shared     = "Shared",
                               deseq_only = "DESeq2-EVS only"),
                        levels = c("CPM-EVS only", "Shared", "DESeq2-EVS only")))

  ggplot(long, aes(comparison, n, fill = set)) +
    geom_col(width = 0.56, colour = "white", linewidth = 0.2) +
    geom_text(aes(label = evs_comma(n)),
              position = position_stack(vjust = 0.5),
              size = pt2mm(5.0), colour = "white") +
    geom_text(data = tot, inherit.aes = FALSE,
              aes(x = comparison, y = total,
                  label = sprintf("italic(J) * \" = \" * \"%.2f\"", jaccard)),
              parse = TRUE, vjust = -0.6, size = pt2mm(5.4), fontface = "bold") +
    scale_fill_manual(values = c("CPM-EVS only"    = EVS_COL$cpm,
                                 "Shared"          = EVS_COL$shared,
                                 "DESeq2-EVS only" = EVS_COL$deseq)) +
    scale_y_continuous(labels = evs_comma,
                       expand = expansion(mult = c(0, 0.20))) +
    labs(tag = tag, title = "High-confidence EVS membership", x = NULL,
         y = "High-confidence PAS (count)") +
    evs_theme() +
    evs_legend_inside(0.5, 1.0, c(0.5, 1), direction = "horizontal")
}

panel_composition <- function(summary_df, comparisons, tag) {
  d <- evs_prep_summary(summary_df, comparisons) %>%
    select(method, comparison, joint_n, opposite_le_n,
           opposite_divergence_n, remainder_cross_n) %>%
    tidyr::pivot_longer(-c("method", "comparison"),
                        names_to = "cls", values_to = "n") %>%
    mutate(cls = factor(recode(cls,
                               joint_n               = "Joint",
                               opposite_le_n         = "Opposite LE",
                               opposite_divergence_n = "Opposite Divergence",
                               remainder_cross_n     = "Opposite Remainder"),
                        levels = c("Joint", "Opposite LE",
                                   "Opposite Divergence", "Opposite Remainder"))) %>%
    group_by(method, comparison) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  ggplot(d, aes(method, pct, fill = cls)) +
    geom_col(width = 0.62, colour = "white", linewidth = 0.2) +
    geom_text(aes(label = ifelse(pct >= 8, sprintf("%.0f", pct), "")),
              position = position_stack(vjust = 0.5), size = pt2mm(4.8)) +
    facet_wrap(~ comparison, nrow = 1, strip.position = "bottom") +
    scale_fill_manual(values = c("Joint"               = EVS_COL$joint,
                                 "Opposite LE"         = EVS_COL$opp_le,
                                 "Opposite Divergence" = EVS_COL$opp_div,
                                 "Opposite Remainder"  = EVS_COL$opp_rem)) +
    scale_x_discrete(labels = c("CPM-EVS" = "CPM", "DESeq2-EVS" = "DESeq2")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(tag = tag,
         title = expression(paste("Site composition at ", italic(k)^"*")),
         x = NULL,
         y = expression(paste("Composition of top-", italic(k)^"*", " union (%)"))) +
    evs_theme() +
    theme(
      strip.background  = element_blank(),
      strip.placement   = "outside",
      strip.text        = element_text(size = EVS_BASE_PT - 0.5,
                                       margin = margin(t = 1.5)),
      panel.spacing.x   = unit(3, "pt"),
      axis.text.x       = element_text(size = EVS_BASE_PT - 1.6, angle = 90,
                                       hjust = 1, vjust = 0.5),
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.background = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}

panel_le_fraction <- function(summary_df, comparisons, tag) {
  d <- evs_prep_summary(summary_df, comparisons) %>%
    mutate(pct = 100 * k_over_leading_edge)
  ggplot(d, aes(comparison, pct, fill = method)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62) +
    geom_text(aes(label = sprintf("%.1f", pct)),
              position = position_dodge(width = 0.72),
              vjust = -0.45, size = pt2mm(5.4)) +
    scale_fill_manual(values = EVS_METHOD_COL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(tag = tag, title = "Fraction of candidate leading edge selected",
         x = NULL,
         y = expression(paste(italic(k)^"*", " as % of leading-edge domain"))) +
    evs_theme() +
    evs_legend_inside(0.99, 0.99, c(1, 1), direction = "horizontal")
}

# =============================================================================
# 6. AUTO-GENERATED FIGURE LEGENDS
# =============================================================================
# Legends are written with the run's real numbers, so the manuscript text
# cannot drift out of sync with the figures.

evs_write_legends <- function(method_objects, summary_df, overlap_df,
                              comparisons, out_root) {

  kn <- lapply(method_objects, function(o) o$knot)
  N  <- kn[[1]]$N

  knot_txt <- paste(vapply(names(kn), function(m) {
    k <- kn[[m]]
    sprintf("%s, c1 = %s and c2 = %s", evs_method_label(m),
            evs_num(k$c1), evs_num(k$c2))
  }, character(1)), collapse = "; ")

  # What fraction of candidate k are actually non-dominated? If the frontier
  # contains nearly all of them, the normalized Pareto utility is doing
  # the work, and the legend should say so rather than imply a sharp trade-off.
  fr_frac <- vapply(names(method_objects), function(m) {
    ob <- method_objects[[m]]
    mean(vapply(names(comparisons), function(nm) {
      nrow(ob$frontiers[[nm]]) / nrow(ob$scans[[nm]])
    }, numeric(1)))
  }, numeric(1))

  kstar_of <- function(m) {
    vapply(names(comparisons),
           function(nm) method_objects[[m]]$pairs[[nm]]$kstar, numeric(1))
  }
  kstar_txt <- paste(vapply(names(method_objects), function(m) {
    sprintf("%s, k* = %s", evs_method_label(m),
            paste(evs_num(kstar_of(m)), collapse = ", "))
  }, character(1)), collapse = "; ")

  jac <- sprintf("%.2f-%.2f",
                 min(overlap_df$jaccard_high_confidence, na.rm = TRUE),
                 max(overlap_df$jaccard_high_confidence, na.rm = TRUE))

  L <- c(
    "# Figure legends",
    "",
    "<!-- Auto-generated by EVS_Manuscript_Figures.R. Do not hand-edit: edit the",
    "     generator instead, so the numbers can never drift from the figures. -->",
    "",
    "## Abbreviations",
    "",
    "CPM, counts per million; c1 and c2, shared lower and upper rank knots;",
    "D(r), cumulative PC1-excess-variance divergence at rank r; DESeq2,",
    "median-of-ratios normalisation; EVS, excess-variance selection; F_E(r),",
    "cumulative excess-variance mass; F_P(r), cumulative PC1 variance-contribution",
    "mass; G_Delta(k), concordance-weighted retained support; HC, high confidence; J, Jaccard",
    "index; k*, selected per-comparison cutoff; LE, leading edge; N, size of the",
    "experiment-wide PAS universe; PAS, polyadenylation site; PC1, first principal",
    "component; R_Delta(k), distance-weighted opposite-arm Remainder disagreement; RT, [AUTHORS: define];",
    "ZT, zeitgeber time.",
    "",
    "> Two abbreviations could not be recovered from the analysis code and must be",
    "> set by the authors before submission: **EVS** (written above as",
    "> \"excess-variance selection\") and **RT**. Edit `evs_write_legends()` once",
    "> the intended expansions are fixed, and they will propagate to every run.",
    "",
    "---",
    "",
    "**Figure 1. Definition of the Remainder, Divergence and Leading-Edge rank",
    "regimes under two normalisation strategies.**",
    sprintf(paste0(
      "(A, B) Cumulative PC1-excess-variance divergence, D(r) = F_E(r) - F_P(r), ",
      "against PAS rank by ascending absolute PC1 loading, for (A) CPM-EVS and ",
      "(B) DESeq2-EVS. Grey lines show all eight experimental arms, the coloured ",
      "line the across-arm median, and the dashed black line the shared two-knot ",
      "piecewise-linear fit estimated jointly across arms by minimising the summed ",
      "squared residual. Vertical dashed lines mark the shared knots (%s), which ",
      "partition the universe of N = %s PAS into the Remainder (rank < c1), ",
      "Divergence (c1 <= rank <= c2) and Leading Edge (rank > c2), shaded left to ",
      "right. (C) Residuals of the shared-knot fit, D_g(r) - Dhat_g(r), binned by ",
      "rank; lines show the across-arm median and shaded envelopes the ",
      "interquartile range. (D) Widths of the three regimes for each method, ",
      "annotated with PAS counts. c1 and c2 define rank regimes common to all arms ",
      "and are not themselves the selected cutoff k*."),
      knot_txt, evs_num(N)),
    "",
    "---",
    "",
    "**Figure 2. Per-comparison cutoff selection under CPM-EVS.**",
    sprintf(paste0(
      "(A-D) Pareto frontier of cross-arm distance-weighted retained support G_Delta(k) against ",
      "distance-weighted opposite-arm Remainder disagreement R_Delta(k) over all candidate cutoffs ",
      "k, shown separately for each RT/ZT comparison; each comparison is optimised ",
      "independently over 1 <= k <= N - c2. The filled diamond marks k*, the ",
      "Pareto-optimal candidate maximizing U(k)=G_norm(k)-R_norm(k) after ",
      "frontier-specific min-max normalization of benefit and contamination. ",
      "Insets plot U(k) across the Pareto frontier. Non-dominated candidates make ",
      "up %.0f%% of the k values scanned under CPM-EVS. Selected values: ",
      "k* = %s for %s respectively; full diagnostics are in Table 1."),
      100 * fr_frac[["CPM_EVS"]],
      paste(evs_num(kstar_of("CPM_EVS")), collapse = ", "),
      paste(evs_comparison_label(names(comparisons)), collapse = ", ")),
    "",
    "---",
    "",
    "**Figure 3. Per-comparison cutoff selection under DESeq2-EVS.**",
    sprintf(paste0(
      "Panels, axes and annotations are as in Figure 2, computed from ",
      "log1p-transformed DESeq2 median-of-ratios normalised counts. Non-dominated ",
      "candidates make up %.0f%% of the k values scanned. Selected values: ",
      "k* = %s."),
      100 * fr_frac[["DESeq2_EVS"]],
      paste(evs_num(kstar_of("DESeq2_EVS")), collapse = ", ")),
    "",
    "---",
    "",
    "**Figure 4. CPM-EVS and DESeq2-EVS compared across all four RT/ZT",
    "comparisons.**",
    sprintf(paste0(
      "(A) Selected cutoff k* for each comparison and method (%s). (B) ",
      "High-confidence PAS membership, partitioned into sites recovered only by ",
      "CPM-EVS, only by DESeq2-EVS, or by both; J is the Jaccard index of the two ",
      "high-confidence sets (range across comparisons, %s). High-confidence sites ",
      "are Joint sites together with disjoint sites whose opposite-arm rank falls ",
      "in the Leading Edge or Divergence regime. (C) Composition of the selected ",
      "top-k* union at k*, as a percentage of union size; segment labels are ",
      "omitted below 8%%. Opposite-Remainder sites are retained in the exported ",
      "union with a penalty flag rather than discarded. (D) k* as a percentage of ",
      "the candidate Leading-Edge domain, N - c2. Bars in A and D are annotated ",
      "with their values; exact counts are given in Tables 1 and 2."),
      kstar_txt, jac),
    "",
    "---",
    "",
    "*Figures were authored at 180 mm final width and are supplied as vector PDF",
    "and 600 dpi raster. Colours follow the Okabe-Ito palette and remain",
    "distinguishable under common forms of colour-vision deficiency and in",
    "greyscale.*"
  )

  writeLines(L, file.path(out_root, "Figure_Legends.md"))
  message("  wrote Figure_Legends.md")
  invisible(TRUE)
}

# =============================================================================
# 7. TOP-LEVEL RENDER
# =============================================================================

evs_render_all <- function(method_objects, summary_df, overlap_df,
                           comparisons, fig_dir, out_root) {

  dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  tags <- LETTERS[1:4]

  message("Figure 1 ...")
  f1 <- list(
    panel_regime("CPM_EVS", method_objects[["CPM_EVS"]]$arms,
                 method_objects[["CPM_EVS"]]$knot, "A", show_legend = TRUE),
    panel_regime("DESeq2_EVS", method_objects[["DESeq2_EVS"]]$arms,
                 method_objects[["DESeq2_EVS"]]$knot, "B"),
    panel_residuals(method_objects, "C"),
    panel_regime_widths(method_objects, "D")
  )
  evs_save(f1, file.path(fig_dir, "Figure_1_Regime_Definition"),
           nrow = 2, ncol = 2, height_mm = 150)

  for (m in names(method_objects)) {
    fig_no <- if (m == "CPM_EVS") 2L else 3L
    message("Figure ", fig_no, " ...")
    ob  <- method_objects[[m]]
    col <- unname(EVS_METHOD_COL[evs_method_label(m)])
    pl  <- lapply(seq_along(comparisons), function(i) {
      nm <- names(comparisons)[i]
      panel_pareto(m, nm, ob$scans[[nm]], ob$frontiers[[nm]],
                   ob$pairs[[nm]]$kstar, tags[i], col)
    })
    evs_save(pl, file.path(fig_dir,
                           sprintf("Figure_%d_%s_Pareto_Cutoffs", fig_no, m)),
             nrow = 2, ncol = 2, height_mm = 145)
  }

  message("Figure 4 ...")
  f4 <- list(
    panel_kstar(summary_df, comparisons, "A"),
    panel_hc_membership(overlap_df, comparisons, "B"),
    panel_composition(summary_df, comparisons, "C"),
    panel_le_fraction(summary_df, comparisons, "D")
  )
  evs_save(f4, file.path(fig_dir, "Figure_4_CPM_vs_DESeq2_EVS_Summary"),
           nrow = 2, ncol = 2, height_mm = 150)

  evs_write_legends(method_objects, summary_df, overlap_df, comparisons, out_root)
  message("Figures complete: ", fig_dir)
  invisible(TRUE)
}

# =============================================================================
# RUN ANALYSIS
# =============================================================================

message("Reading count matrix...")
counts <- read_counts(COUNT_FILE)
groups <- assign_groups(colnames(counts))

message(
  "Features: ",nrow(counts),
  " | samples: ",ncol(counts),
  " | arms: ",length(levels(groups))
)

message("Global DESeq2 normalization...")
norm_obj <- deseq2_normalize(counts,groups)
norm <- norm_obj$counts

message("Pooled normalized within-group variance...")
pooled_var <- pooled_within_group_var(norm,groups)

method_objects <- list()
summary_rows <- list()
site_rows <- list()
scan_rows <- list()
frontier_rows <- list()

for (method in METHODS) {

  message("============================================================")
  message("METHOD: ",method)

  arms <- list()

  for (g in levels(groups)) {
    idx <- which(groups==g)

    arms[[g]] <- build_arm(
      method=method,
      arm=g,
      raw=counts[,idx,drop=FALSE],
      norm=norm[,idx,drop=FALSE],
      pooled_var=pooled_var
    )
  }

  knot <- fit_shared_knots(arms)
  K <- knot$N-knot$c2

  message(
    "Shared c1=",knot$c1,
    " | c2=",knot$c2,
    " | candidate k=1..",K
  )

  scans <- list()
  frontiers <- list()
  pair_results <- list()

  for (nm in names(COMPARISONS)) {
    mp <- COMPARISONS[[nm]]

    scan <- scan_pair_fast(
      arms[[mp[["control"]]]],
      arms[[mp[["treatment"]]]],
      knot$c1,
      knot$c2
    )

    opt <- select_weighted_pareto_optimum(
      scan, good_col="good_score", cost_col="remainder_distance_cost",
      benefit_weight=1, contamination_weight=1
    )
    scan <- opt$scan
    frontier <- opt$frontier
    kstar <- opt$selected_k

    selected <- scan[scan$k==kstar,,drop=FALSE]

    cls <- classify_at_k(
      k=kstar,
      comparison=nm,
      control_df=arms[[mp[["control"]]]],
      treatment_df=arms[[mp[["treatment"]]]],
      c1=knot$c1,
      c2=knot$c2
    )

    cls$method <- method

    high_conf_n <- sum(cls$high_confidence)
    rem_n <- sum(cls$penalized_remainder_crossing)

    row <- data.frame(
      method=method,
      comparison=nm,
      N=knot$N,
      c1=knot$c1,
      c2=knot$c2,
      remainder_size=knot$c1-1L,
      divergence_size=knot$c2-knot$c1+1L,
      leading_edge_candidate_size=K,
      selected_k=kstar,
      cutoff_rank=knot$N-kstar+1L,
      k_over_N=kstar/knot$N,
      k_over_leading_edge=kstar/K,
      good_score=selected$good_score,
      remainder_distance_cost=selected$remainder_distance_cost,
      good_n=selected$good_n,
      remainder_cross_n=selected$remainder_cross_n,
      good_norm=selected$good_norm,
      remainder_norm=selected$remainder_norm,
      weighted_utility=selected$weighted_utility,
      joint_n=selected$joint_n,
      opposite_le_n=selected$opposite_le_n,
      opposite_divergence_n=selected$opposite_divergence_n,
      selected_union_n=selected$union_n,
      high_confidence_n=high_conf_n,
      penalized_remainder_n=rem_n,
      knot_SSE=knot$SSE,
      stringsAsFactors=FALSE
    )

    scans[[nm]] <- scan
    frontiers[[nm]] <- frontier
    pair_results[[nm]] <- list(
      kstar=kstar,
      sites=cls,
      summary=row
    )

    key <- paste(method,nm,sep="__")

    summary_rows[[key]] <- row
    site_rows[[key]] <- cls
    scan_rows[[key]] <- mutate(
      scan,method=method,comparison=nm
    )
    frontier_rows[[key]] <- mutate(
      frontier,method=method,comparison=nm
    )

    message(
      "  ",nm,
      " | k*=",kstar,
      " | Gdist=",round(selected$good_score,3),
      " | Rdist=",round(selected$remainder_distance_cost,3),
      " | U=",round(selected$weighted_utility,4),
      " | Joint=",selected$joint_n,
      " | Rem-cross=",selected$remainder_cross_n
    )
  }

  method_objects[[method]] <- list(
    arms=arms,
    knot=knot,
    scans=scans,
    frontiers=frontiers,
    pairs=pair_results
  )
}

summary_df <- bind_rows(summary_rows)
sites_df <- bind_rows(site_rows)
scans_df <- bind_rows(scan_rows)
frontiers_df <- bind_rows(frontier_rows)

# =============================================================================
# SUMMARY / OVERLAP TABLES
# =============================================================================

overlap_df <- bind_rows(lapply(names(COMPARISONS),\(nm) {

  cpm <- sites_df %>%
    filter(
      method=="CPM_EVS",
      comparison==nm,
      high_confidence
    ) %>%
    pull(feature_id) %>%
    unique()

  des <- sites_df %>%
    filter(
      method=="DESeq2_EVS",
      comparison==nm,
      high_confidence
    ) %>%
    pull(feature_id) %>%
    unique()

  u <- union(cpm,des)

  data.frame(
    comparison=nm,
    CPM_high_confidence=length(cpm),
    DESeq2_high_confidence=length(des),
    overlap_high_confidence=length(intersect(cpm,des)),
    union_high_confidence=length(u),
    jaccard_high_confidence=ifelse(
      length(u)>0,
      length(intersect(cpm,des))/length(u),
      NA_real_
    ),
    stringsAsFactors=FALSE
  )
}))

write.csv(
  summary_df,
  file.path(TAB_DIR,"Table_1_EVS_Cutoff_Summary.csv"),
  row.names=FALSE
)

write.csv(
  overlap_df,
  file.path(TAB_DIR,"Table_2_CPM_vs_DESeq2_HighConfidence_Overlap.csv"),
  row.names=FALSE
)

write.csv(
  sites_df,
  file.path(TAB_DIR,"Table_3_All_Selected_TopK_Union_Sites.csv"),
  row.names=FALSE
)

write.csv(
  sites_df %>% filter(high_confidence),
  file.path(TAB_DIR,"Table_4_HighConfidence_EVS_Sites.csv"),
  row.names=FALSE
)

write.csv(
  sites_df %>% filter(penalized_remainder_crossing),
  file.path(TAB_DIR,"Table_5_Penalized_Remainder_Crossing_Sites.csv"),
  row.names=FALSE
)

write.csv(
  scans_df,
  file.path(TAB_DIR,"Table_S1_All_Cutoff_Scans.csv"),
  row.names=FALSE
)

write.csv(
  frontiers_df,
  file.path(TAB_DIR,"Table_S2_All_Pareto_Frontiers.csv"),
  row.names=FALSE
)

# =============================================================================
# MANUSCRIPT FIGURES
# =============================================================================

# Clear figures from earlier runs. The manifest and zip steps below glob every
# PNG/PDF in FIG_DIR, so stale files from a previous naming scheme would
# otherwise be packaged alongside the current ones.
old_figure_files <- list.files(
  FIG_DIR,
  pattern="\\.(png|pdf|tif|tiff)$",
  full.names=TRUE,
  ignore.case=TRUE
)

if (length(old_figure_files)) unlink(old_figure_files)

# Cache the fitted objects so figures can be re-rendered later without re-running
# the whole pipeline:
#   z <- readRDS(file.path(OUT_ROOT,"analysis_objects.rds"))
#   evs_render_all(z$method_objects, z$summary_df, z$overlap_df,
#                  COMPARISONS, FIG_DIR, OUT_ROOT)
saveRDS(
  list(
    method_objects=method_objects,
    summary_df=summary_df,
    overlap_df=overlap_df,
    sites_df=sites_df
  ),
  file.path(OUT_ROOT,"analysis_objects.rds")
)

evs_render_all(
  method_objects = method_objects,
  summary_df     = summary_df,
  overlap_df     = overlap_df,
  comparisons    = COMPARISONS,
  fig_dir        = FIG_DIR,
  out_root       = OUT_ROOT
)

# =============================================================================
# MANUSCRIPT METHODS
# (Figure legends are generated by evs_write_legends(), above.)
# =============================================================================

methods_lines <- c(
  "# Methods",
  "",
  "Two EVS preprocessing strategies were evaluated using a common experiment-wide PAS universe. CPM-EVS used log1p-transformed counts per million, whereas DESeq2-EVS used log1p-transformed DESeq2 median-of-ratios normalized counts. PCA was performed independently within each experimental arm, and PASs were ranked from lowest to highest absolute PC1 loading.",
  "",
  "For each PAS, PC1 variance contribution was P_i=lambda_1*v_i1^2. A pooled within-group variance was calculated from globally DESeq2-normalized counts across the eight experimental arms, and arm-specific excess variance was E_ig=max(V_pool,i-mu_ig,0). PC1 contribution and excess variance were converted to rank-wise probability masses and cumulative distributions. Cumulative divergence was D_g(r)=F_E,g(r)-F_P,g(r).",
  "",
  "For each EVS method separately, the eight arm-specific D_g(r) curves were jointly fit using a continuous two-knot piecewise-linear model. The shared c1 and c2 values minimized the summed squared residual error across all arms. Rank<c1 defined the Remainder, c1<=rank<=c2 the Divergence interval, and rank>c2 the Leading Edge.",
  "",
  "Each RT/ZT comparison was optimized independently over 1<=k<=N-c2. For each disjoint PAS, cross-arm rank distance was Delta=|r_selected-r_opposite|/(N-1). Joint PASs received full unit benefit. Disjoint PASs whose opposite-arm rank remained in the Leading Edge or Divergence interval received support 1-Delta; disjoint PASs whose opposite-arm rank crossed below c1 into the Remainder received contamination cost Delta. Thus the regime boundaries determine whether a disjoint PAS is supportive or contaminating, while its observed cross-arm rank distance determines the magnitude.",
  "",
  "For each candidate k, G_Delta(k)=N_Joint+sum_supported(1-Delta_i) and R_Delta(k)=sum_remainder(Delta_i). Nondominated [R_Delta(k),G_Delta(k)] candidates defined the Pareto frontier. On that frontier, benefit and contamination were min-max normalized and, with equal weights, U(k)=G_norm(k)-R_norm(k). The empirical k* was the Pareto-optimal candidate maximizing U(k); ties were resolved by greater G_Delta, lower R_Delta, then larger k. Raw Joint/Leading-Edge/Divergence/Remainder counts were retained as diagnostics and all site-level distances and scores were exported."
)

writeLines(
  methods_lines,
  file.path(OUT_ROOT,"Methods_Manuscript.md")
)

# =============================================================================
# MANIFESTS + ZIP FILES
# =============================================================================

figure_files <- list.files(
  FIG_DIR,
  pattern="\\.(png|pdf)$",
  full.names=TRUE
)

table_files <- list.files(
  TAB_DIR,
  pattern="\\.csv$",
  full.names=TRUE
)

fig_manifest <- data.frame(
  file=basename(figure_files),
  type=ifelse(grepl("\\.pdf$",figure_files,ignore.case=TRUE),"PDF","PNG"),
  stringsAsFactors=FALSE
)

tab_manifest <- data.frame(
  file=basename(table_files),
  stringsAsFactors=FALSE
)

write.csv(
  fig_manifest,
  file.path(OUT_ROOT,"Manifest_Manuscript_Figures.csv"),
  row.names=FALSE
)

write.csv(
  tab_manifest,
  file.path(OUT_ROOT,"Manifest_Manuscript_Tables.csv"),
  row.names=FALSE
)

zip_folder <- function(zipfile,files,root) {
  files <- files[file.exists(files)]
  if (!length(files)) return(FALSE)

  if (file.exists(zipfile)) unlink(zipfile)

  old <- getwd()
  on.exit(setwd(old),add=TRUE)
  setwd(root)

  utils::zip(
    zipfile=zipfile,
    files=basename(files)
  )

  file.exists(zipfile)
}

FIG_ZIP <- file.path(OUT_ROOT,"Manuscript_Ready_Figures.zip")
TAB_ZIP <- file.path(OUT_ROOT,"Manuscript_Ready_Tables.zip")
ALL_ZIP <- file.path(OUT_ROOT,"CPM_vs_DESeq2_EVS_Complete_Outputs.zip")

# Figure zip.
old <- getwd()
setwd(FIG_DIR)
if (file.exists(FIG_ZIP)) unlink(FIG_ZIP)
utils::zip(FIG_ZIP,files=basename(figure_files))
setwd(old)

# Table zip.
old <- getwd()
setwd(TAB_DIR)
if (file.exists(TAB_ZIP)) unlink(TAB_ZIP)
utils::zip(TAB_ZIP,files=basename(table_files))
setwd(old)

# Complete package zip from OUT_ROOT using recursive relative paths.
all_package_files <- c(
  list.files("Figures",recursive=TRUE,full.names=TRUE),
  list.files("Tables",recursive=TRUE,full.names=TRUE),
  "Methods_Manuscript.md",
  "Figure_Legends.md",
  "Manifest_Manuscript_Figures.csv",
  "Manifest_Manuscript_Tables.csv"
)

old <- getwd()
setwd(OUT_ROOT)
if (file.exists(ALL_ZIP)) unlink(ALL_ZIP)
utils::zip(ALL_ZIP,files=all_package_files)
setwd(old)

message("============================================================")
message("ANALYSIS COMPLETE")
message("Output root: ",OUT_ROOT)
message("Figures ZIP: ",FIG_ZIP)
message("Tables ZIP: ",TAB_ZIP)
message("Complete ZIP: ",ALL_ZIP)
message("============================================================")
print(summary_df)
