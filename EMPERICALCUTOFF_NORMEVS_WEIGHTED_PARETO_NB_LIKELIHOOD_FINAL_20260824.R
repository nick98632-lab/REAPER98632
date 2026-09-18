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
# 5) Cutoff scoring is continuous:
#      - Joint sites receive full benefit.
#      - Disjoint opposite-arm Leading-Edge sites receive graded LE support.
#      - Disjoint opposite-arm Divergence sites receive graded divergence support.
#      - Disjoint opposite-arm Remainder sites receive a distance-weighted penalty.
# 6) k* is the geometric knee of the comparison-specific weighted Pareto frontier.
# 7) Remainder-crossing sites are NOT silently discarded from the top-k union.
#    They are flagged as penalized/low-confidence. A separate high-confidence
#    subset contains Joint + opposite-LE + opposite-Divergence sites.
# 8) Exactly four composite manuscript figures are generated.
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
PNG_DPI <- 360

COL <- list(
  control="#0072B2",
  treatment="#009E73",
  observed="#6A3D9A",
  fit="#222222",
  c1="#D73027",
  c2="#1A9850",
  selected="#B5179E",
  remainder="#DCE6F2",
  divergence="#FFF0B3",
  leading="#D8F3E7",
  cpm="#3B82F6",
  deseq="#10B981"
)

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
# FAST CONTINUOUS TOP-k SCORING
# =============================================================================

rank_map <- function(df) setNames(df$rank,df$feature_id)

div_support <- function(r,c1,c2) {
  pmin(pmax((r-c1)/(c2-c1),0),1)
}

lead_support <- function(r,c2,N) {
  pmin(pmax((r-c2)/(N-c2),0),1)
}

rem_depth <- function(r,c1) {
  pmin(pmax((c1-r)/(c1-1),0),1)
}

rank_gap <- function(a,b,N) {
  abs(a-b)/(N-1)
}

add_events <- function(pos,weight,L) {
  out <- numeric(L)
  ok <- is.finite(pos) & is.finite(weight) & pos>=1 & pos<=L

  if (!any(ok)) return(out)

  z <- rowsum(
    weight[ok],
    group=as.integer(pos[ok]),
    reorder=FALSE
  )

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

  s <- as.integer(start[ok])
  e <- pmin(as.integer(end[ok]),K+1L)
  w <- weight[ok]

  delta <- delta + add_events(s,w,K+1L)
  delta <- delta - add_events(e,w,K+1L)

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

  # Joint sites: full unit benefit after both arms include the PAS.
  joint_depth <- pmax(dC,dT)
  joint_n <- activate_from(joint_depth,rep(1,N),K)

  # Control-only interval.
  idxC <- which(dC<dT & dC<=K)
  startC <- dC[idxC]
  endC <- pmin(dT[idxC],K+1L)
  oppC <- rT[idxC]
  selC <- rC[idxC]

  supportC <- ifelse(
    oppC>c2,
    lead_support(oppC,c2,N),
    ifelse(oppC>=c1,div_support(oppC,c1,c2),0)
  )

  penaltyC <- ifelse(
    oppC<c1,
    rem_depth(oppC,c1)^2 * rank_gap(selC,oppC,N),
    0
  )

  disC_n       <- active_interval(startC,endC,rep(1,length(idxC)),K)
  disC_support <- active_interval(startC,endC,supportC,K)
  disC_penalty <- active_interval(startC,endC,penaltyC,K)
  disC_rem     <- active_interval(startC,endC,as.numeric(oppC<c1),K)
  disC_div     <- active_interval(startC,endC,as.numeric(oppC>=c1 & oppC<=c2),K)
  disC_le      <- active_interval(startC,endC,as.numeric(oppC>c2),K)

  # Treatment-only interval.
  idxT <- which(dT<dC & dT<=K)
  startT <- dT[idxT]
  endT <- pmin(dC[idxT],K+1L)
  oppT <- rC[idxT]
  selT <- rT[idxT]

  supportT <- ifelse(
    oppT>c2,
    lead_support(oppT,c2,N),
    ifelse(oppT>=c1,div_support(oppT,c1,c2),0)
  )

  penaltyT <- ifelse(
    oppT<c1,
    rem_depth(oppT,c1)^2 * rank_gap(selT,oppT,N),
    0
  )

  disT_n       <- active_interval(startT,endT,rep(1,length(idxT)),K)
  disT_support <- active_interval(startT,endT,supportT,K)
  disT_penalty <- active_interval(startT,endT,penaltyT,K)
  disT_rem     <- active_interval(startT,endT,as.numeric(oppT<c1),K)
  disT_div     <- active_interval(startT,endT,as.numeric(oppT>=c1 & oppT<=c2),K)
  disT_le      <- active_interval(startT,endT,as.numeric(oppT>c2),K)

  disjoint_n <- disC_n+disT_n

  data.frame(
    k=seq_len(K),
    joint_n=joint_n,
    disjoint_n=disjoint_n,
    opposite_le_n=disC_le+disT_le,
    opposite_divergence_n=disC_div+disT_div,
    remainder_cross_n=disC_rem+disT_rem,
    union_n=2*seq_len(K)-joint_n,
    weighted_benefit=joint_n+disC_support+disT_support,
    weighted_penalty=disC_penalty+disT_penalty,
    stringsAsFactors=FALSE
  )
}

# =============================================================================
# PARETO FRONTIER + GEOMETRIC KNEE
# =============================================================================

norm01 <- function(x) {
  r <- range(x,na.rm=TRUE)
  if (!is.finite(diff(r)) || diff(r)==0) return(rep(0,length(x)))
  (x-r[1])/diff(r)
}

pareto_frontier <- function(df) {
  x <- df %>%
    arrange(weighted_penalty,desc(weighted_benefit),desc(k))

  keep <- logical(nrow(x))
  best <- -Inf

  for (i in seq_len(nrow(x))) {
    if (x$weighted_benefit[i] > best) {
      keep[i] <- TRUE
      best <- x$weighted_benefit[i]
    }
  }

  x[keep,,drop=FALSE] %>% arrange(weighted_penalty)
}

pareto_knee <- function(frontier) {
  if (nrow(frontier)==1L) return(frontier$k[1])

  x <- norm01(frontier$weighted_penalty)
  y <- norm01(frontier$weighted_benefit)

  x1 <- x[1]
  y1 <- y[1]
  x2 <- x[length(x)]
  y2 <- y[length(y)]

  den <- sqrt((y2-y1)^2+(x2-x1)^2)

  if (!is.finite(den) || den==0) {
    return(frontier$k[which.max(y-x)])
  }

  distance <- abs(
    (y2-y1)*x - (x2-x1)*y + x2*y1 - y2*x1
  )/den

  frontier$k[which.max(distance)]
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

  penalty <- ifelse(
    joint,0,
    ifelse(
      opp<c1,
      rem_depth(opp,c1)^2 * rank_gap(selected_rank,opp,N),
      0
    )
  )

  support <- ifelse(
    joint,1,
    ifelse(
      opp>c2,
      lead_support(opp,c2,N),
      ifelse(opp>=c1,div_support(opp,c1,c2),0)
    )
  )

  data.frame(
    comparison=comparison,
    feature_id=ids,
    selected_k=k,
    control_rank=a,
    treatment_rank=b,
    class=ifelse(joint,"Joint",ifelse(disC,"Control only","Treatment only")),
    opposite_region=opp_region,
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
# FIGURE HELPERS
# =============================================================================

theme_pub <- function(base=11.5) {
  theme_classic(base_size=base) +
    theme(
      plot.title=element_text(face="bold",size=base+1.2),
      plot.subtitle=element_text(size=base-.8,color="grey30"),
      axis.title=element_text(face="bold"),
      axis.text=element_text(color="grey20"),
      panel.border=element_rect(fill=NA,color="grey78",linewidth=.4),
      legend.position="top",
      legend.title=element_blank(),
      plot.margin=margin(8,10,8,10)
    )
}

save_composite <- function(plots,file,nrow,ncol,width,height) {
  png(
    file,width=width,height=height,units="in",
    res=PNG_DPI,bg="white"
  )

  grid.newpage()
  pushViewport(
    viewport(
      layout=grid.layout(
        nrow=nrow,
        ncol=ncol,
        widths=unit(rep(1,ncol),"null"),
        heights=unit(rep(1,nrow),"null")
      )
    )
  )

  for (i in seq_along(plots)) {
    r <- ((i-1)%/%ncol)+1
    c <- ((i-1)%%ncol)+1
    print(plots[[i]],vp=viewport(layout.pos.row=r,layout.pos.col=c))
  }

  dev.off()

  pdf(
    sub("\\.png$",".pdf",file),
    width=width,height=height,useDingbats=FALSE
  )

  grid.newpage()
  pushViewport(
    viewport(
      layout=grid.layout(
        nrow=nrow,
        ncol=ncol,
        widths=unit(rep(1,ncol),"null"),
        heights=unit(rep(1,nrow),"null")
      )
    )
  )

  for (i in seq_along(plots)) {
    r <- ((i-1)%/%ncol)+1
    c <- ((i-1)%%ncol)+1
    print(plots[[i]],vp=viewport(layout.pos.row=r,layout.pos.col=c))
  }

  dev.off()
}

make_regime_panel <- function(method,arms,knot,panel_label) {
  obs <- bind_rows(lapply(arms,\(z) select(z,rank,D)))

  med <- obs %>%
    group_by(rank) %>%
    summarise(D=median(D),.groups="drop")

  fit <- data.frame(
    rank=seq_len(knot$N),
    D=apply(knot$fitted,1,median)
  )

  method_name <- ifelse(method=="CPM_EVS","CPM-EVS","DESeq2-EVS")

  ggplot() +
    annotate(
      "rect",xmin=1,xmax=knot$c1,ymin=-Inf,ymax=Inf,
      fill=COL$remainder,alpha=.34
    ) +
    annotate(
      "rect",xmin=knot$c1,xmax=knot$c2,ymin=-Inf,ymax=Inf,
      fill=COL$divergence,alpha=.34
    ) +
    annotate(
      "rect",xmin=knot$c2,xmax=knot$N,ymin=-Inf,ymax=Inf,
      fill=COL$leading,alpha=.34
    ) +
    geom_hline(yintercept=0,color="grey70",linetype="dotted") +
    geom_line(
      data=med,aes(rank,D),
      color=COL$observed,linewidth=1.05
    ) +
    geom_line(
      data=fit,aes(rank,D),
      color=COL$fit,linewidth=1.15
    ) +
    geom_vline(
      xintercept=knot$c1,
      color=COL$c1,linetype="dashed",linewidth=.75
    ) +
    geom_vline(
      xintercept=knot$c2,
      color=COL$c2,linetype="longdash",linewidth=.75
    ) +
    annotate(
      "label",
      x=knot$c1,y=Inf,
      label=paste0("c1 = ",knot$c1),
      hjust=1.03,vjust=1.18,size=3,
      fill="white",label.size=.15
    ) +
    annotate(
      "label",
      x=knot$c2,y=Inf,
      label=paste0("c2 = ",knot$c2),
      hjust=-.03,vjust=1.18,size=3,
      fill="white",label.size=.15
    ) +
    labs(
      title=paste0(panel_label,"  ",method_name),
      subtitle=paste0(
        "Shared 8-arm fit | Divergence width = ",
        knot$c2-knot$c1+1,
        " | Candidate LE = ",knot$N-knot$c2
      ),
      x="PC1 rank",
      y="D(r) = F_E(r) - F_P(r)"
    ) +
    theme_pub() +
    theme(legend.position="none")
}

make_regime_explanation_panel <- function() {
  df <- data.frame(
    x=c(.17,.50,.83),
    y=1,
    regime=c("Remainder","Divergence","Leading edge"),
    definition=c(
      "rank < c1",
      "c1 <= rank <= c2",
      "rank > c2"
    )
  )

  ggplot(df,aes(x,y)) +
    annotate("rect",xmin=.02,xmax=.32,ymin=.55,ymax=1.45,
             fill=COL$remainder) +
    annotate("rect",xmin=.34,xmax=.66,ymin=.55,ymax=1.45,
             fill=COL$divergence) +
    annotate("rect",xmin=.68,xmax=.98,ymin=.55,ymax=1.45,
             fill=COL$leading) +
    geom_text(aes(label=regime),fontface="bold",size=4) +
    geom_text(aes(y=.78,label=definition),size=3.3) +
    annotate(
      "text",x=.5,y=.26,
      label="c1/c2 define the common rank regimes; they are not the final k*.",
      fontface="bold",size=3.6
    ) +
    coord_cartesian(xlim=c(0,1),ylim=c(0,1.55),clip="off") +
    labs(title="C  Regime interpretation") +
    theme_void(base_size=12) +
    theme(
      plot.title=element_text(face="bold",size=13),
      plot.margin=margin(15,15,15,15)
    )
}

make_equation_panel <- function() {
  txt <- paste(
    "Shared-knot model:",
    "D_g(x) = beta_0g + beta_1g x + gamma_1g(x-c1)+ + gamma_2g(x-c2)+",
    "",
    "x = (rank - 1)/(N - 1)",
    "",
    "Choose c1,c2 to minimize:",
    "SUM_g SUM_r [D_g(r) - Dhat_g(r; c1,c2)]^2",
    "",
    "Then each RT/ZT pair independently chooses k* within rank > c2.",
    sep="\n"
  )

  ggplot() +
    annotate(
      "label",x=.5,y=.53,label=txt,
      hjust=.5,vjust=.5,size=3.6,
      lineheight=1.18,fill="white",label.size=.3
    ) +
    coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip="off") +
    labs(title="D  How c1 and c2 are chosen") +
    theme_void(base_size=12) +
    theme(
      plot.title=element_text(face="bold",size=13),
      plot.margin=margin(15,15,15,15)
    )
}

make_pareto_panel <- function(
    method,comparison,scan,frontier,kstar,knot,panel_label) {

  selected <- scan[scan$k==kstar,,drop=FALSE]
  if (nrow(selected)!=1L) stop("Selected k not found in scan.")

  method_name <- ifelse(method=="CPM_EVS","CPM-EVS","DESeq2-EVS")

  high_conf <- selected$joint_n +
    selected$opposite_le_n +
    selected$opposite_divergence_n

  label <- paste0(
    "c1=",knot$c1,"   c2=",knot$c2,
    "\nk*=",kstar,
    "\nG_w=",formatC(selected$weighted_benefit,digits=1,format="f"),
    "   R_w=",formatC(selected$weighted_penalty,digits=2,format="f"),
    "\nJoint=",selected$joint_n,
    "   HC=",high_conf,
    "\nRem-cross=",selected$remainder_cross_n
  )

  ggplot() +
    geom_path(
      data=scan,
      aes(weighted_penalty,weighted_benefit),
      color="grey83",linewidth=.5
    ) +
    geom_path(
      data=frontier,
      aes(weighted_penalty,weighted_benefit),
      color=COL$observed,linewidth=1.25
    ) +
    geom_point(
      data=selected,
      aes(weighted_penalty,weighted_benefit),
      shape=23,size=4.5,
      fill=COL$selected,color=COL$selected
    ) +
    annotate(
      "label",
      x=selected$weighted_penalty,
      y=selected$weighted_benefit,
      label=label,
      hjust=ifelse(
        selected$weighted_penalty > median(scan$weighted_penalty),
        1.03,-.03
      ),
      vjust=1.05,
      size=2.9,
      fill="white",
      label.size=.15
    ) +
    labs(
      title=paste0(panel_label,"  ",comparison),
      subtitle=paste0(method_name," | geometric Pareto knee"),
      x="Distance-weighted Remainder penalty",
      y="Weighted retained-site benefit"
    ) +
    theme_pub(base=10.5) +
    theme(legend.position="none")
}

make_cutoff_bar_panel <- function(summary_df) {
  dat <- summary_df %>%
    mutate(
      method=factor(
        method,
        levels=c("CPM_EVS","DESeq2_EVS"),
        labels=c("CPM-EVS","DESeq2-EVS")
      ),
      comparison=factor(comparison,levels=names(COMPARISONS))
    )

  ggplot(dat,aes(comparison,selected_k,fill=method)) +
    geom_col(
      position=position_dodge(width=.74),
      width=.64
    ) +
    geom_text(
      aes(label=selected_k),
      position=position_dodge(width=.74),
      vjust=-.4,
      fontface="bold",
      size=3.3
    ) +
    labs(
      title="A  Dataset-specific empirical cutoffs",
      subtitle="Eight independent k* values: 4 CPM-EVS + 4 DESeq2-EVS",
      x=NULL,
      y="Selected k*",
      fill=NULL
    ) +
    theme_pub() +
    theme(legend.position="top") +
    expand_limits(y=max(dat$selected_k)*1.15)
}

make_overlap_bar_panel <- function(overlap_df) {
  dat <- overlap_df %>%
    select(
      comparison,
      CPM_high_confidence,
      DESeq2_high_confidence,
      overlap_high_confidence
    ) %>%
    pivot_longer(
      -comparison,
      names_to="metric",
      values_to="n"
    ) %>%
    mutate(
      metric=recode(
        metric,
        CPM_high_confidence="CPM high-confidence",
        DESeq2_high_confidence="DESeq2 high-confidence",
        overlap_high_confidence="Shared high-confidence"
      ),
      comparison=factor(comparison,levels=names(COMPARISONS))
    )

  ggplot(dat,aes(comparison,n,fill=metric)) +
    geom_col(
      position=position_dodge(width=.76),
      width=.68
    ) +
    labs(
      title="B  High-confidence EVS membership",
      subtitle="Joint + opposite-Leading-Edge + opposite-Divergence sites",
      x=NULL,
      y="PAS count",
      fill=NULL
    ) +
    theme_pub(base=10.8) +
    theme(
      legend.position="top",
      legend.text=element_text(size=8.8)
    )
}

make_remainder_bar_panel <- function(summary_df) {
  dat <- summary_df %>%
    mutate(
      method=factor(
        method,
        levels=c("CPM_EVS","DESeq2_EVS"),
        labels=c("CPM-EVS","DESeq2-EVS")
      ),
      comparison=factor(comparison,levels=names(COMPARISONS))
    )

  ggplot(dat,aes(comparison,remainder_cross_n,fill=method)) +
    geom_col(
      position=position_dodge(width=.74),
      width=.64
    ) +
    geom_text(
      aes(label=remainder_cross_n),
      position=position_dodge(width=.74),
      vjust=-.35,size=3
    ) +
    labs(
      title="C  Penalized Remainder crossings at k*",
      subtitle="These remain flagged in the selected top-k union; they are not silently discarded",
      x=NULL,
      y="Remainder-crossing PASs",
      fill=NULL
    ) +
    theme_pub(base=10.8) +
    theme(legend.position="top") +
    expand_limits(y=max(dat$remainder_cross_n)*1.15)
}

make_method_fraction_panel <- function(summary_df) {
  dat <- summary_df %>%
    mutate(
      method=factor(
        method,
        levels=c("CPM_EVS","DESeq2_EVS"),
        labels=c("CPM-EVS","DESeq2-EVS")
      ),
      comparison=factor(comparison,levels=names(COMPARISONS)),
      percent=100*k_over_leading_edge
    )

  ggplot(dat,aes(comparison,percent,fill=method)) +
    geom_col(
      position=position_dodge(width=.74),
      width=.64
    ) +
    geom_text(
      aes(label=paste0(round(percent,1),"%")),
      position=position_dodge(width=.74),
      vjust=-.35,size=3
    ) +
    labs(
      title="D  Fraction of candidate Leading Edge selected",
      subtitle="k* / (N - c2)",
      x=NULL,
      y="% of Leading-edge candidate domain",
      fill=NULL
    ) +
    theme_pub(base=10.8) +
    theme(legend.position="top") +
    expand_limits(y=max(dat$percent)*1.15)
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

    frontier <- pareto_frontier(scan)
    kstar <- pareto_knee(frontier)

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
      weighted_benefit=selected$weighted_benefit,
      weighted_penalty=selected$weighted_penalty,
      joint_n=selected$joint_n,
      opposite_le_n=selected$opposite_le_n,
      opposite_divergence_n=selected$opposite_divergence_n,
      remainder_cross_n=selected$remainder_cross_n,
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
      " | G_w=",round(selected$weighted_benefit,2),
      " | R_w=",round(selected$weighted_penalty,3),
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
# FOUR COMPOSITE MANUSCRIPT FIGURES
# =============================================================================

# Figure 1: Method/regime definition.
fig1_plots <- list(
  make_regime_panel(
    "CPM_EVS",
    method_objects[["CPM_EVS"]]$arms,
    method_objects[["CPM_EVS"]]$knot,
    "A"
  ),
  make_regime_panel(
    "DESeq2_EVS",
    method_objects[["DESeq2_EVS"]]$arms,
    method_objects[["DESeq2_EVS"]]$knot,
    "B"
  ),
  make_regime_explanation_panel(),
  make_equation_panel()
)

save_composite(
  fig1_plots,
  file.path(FIG_DIR,"Figure_1_Regime_Definition_and_Shared_Knots.png"),
  nrow=2,ncol=2,width=15,height=10.5
)

# Figure 2: CPM Pareto panels.
letters4 <- LETTERS[1:4]

fig2_plots <- lapply(seq_along(COMPARISONS),\(i) {
  nm <- names(COMPARISONS)[i]
  obj <- method_objects[["CPM_EVS"]]

  make_pareto_panel(
    method="CPM_EVS",
    comparison=nm,
    scan=obj$scans[[nm]],
    frontier=obj$frontiers[[nm]],
    kstar=obj$pairs[[nm]]$kstar,
    knot=obj$knot,
    panel_label=letters4[i]
  )
})

save_composite(
  fig2_plots,
  file.path(FIG_DIR,"Figure_2_CPM_EVS_Dataset_Pareto_Cutoffs.png"),
  nrow=2,ncol=2,width=15,height=11
)

# Figure 3: DESeq2 Pareto panels.
fig3_plots <- lapply(seq_along(COMPARISONS),\(i) {
  nm <- names(COMPARISONS)[i]
  obj <- method_objects[["DESeq2_EVS"]]

  make_pareto_panel(
    method="DESeq2_EVS",
    comparison=nm,
    scan=obj$scans[[nm]],
    frontier=obj$frontiers[[nm]],
    kstar=obj$pairs[[nm]]$kstar,
    knot=obj$knot,
    panel_label=letters4[i]
  )
})

save_composite(
  fig3_plots,
  file.path(FIG_DIR,"Figure_3_DESeq2_EVS_Dataset_Pareto_Cutoffs.png"),
  nrow=2,ncol=2,width=15,height=11
)

# Figure 4: Final quantitative summary bars.
fig4_plots <- list(
  make_cutoff_bar_panel(summary_df),
  make_overlap_bar_panel(overlap_df),
  make_remainder_bar_panel(summary_df),
  make_method_fraction_panel(summary_df)
)

save_composite(
  fig4_plots,
  file.path(FIG_DIR,"Figure_4_CPM_vs_DESeq2_EVS_Summary.png"),
  nrow=2,ncol=2,width=15,height=10.5
)

# =============================================================================
# MANUSCRIPT METHODS + FIGURE LEGENDS
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
  "Each RT/ZT comparison was then optimized independently over 1<=k<=N-c2. Joint top-k PASs received full benefit. For disjoint PASs, opposite-arm Leading-Edge support increased with (r-c2)/(N-c2), and opposite-arm Divergence support increased with (r-c1)/(c2-c1). Opposite-arm Remainder crossings were penalized by [(c1-r_opp)/(c1-1)]^2 * |r_C-r_T|/(N-1).",
  "",
  "For every candidate k, weighted retained-site benefit and weighted Remainder disagreement were calculated. Nondominated candidates defined the Pareto frontier, and the empirical k* was selected as the geometric knee of that frontier. Remainder-crossing PASs were retained in the exported top-k union with an explicit penalty flag; a separate high-confidence subset contained Joint, opposite-Leading-Edge, and opposite-Divergence PASs."
)

writeLines(
  methods_lines,
  file.path(OUT_ROOT,"Methods_Manuscript.md")
)

legend_lines <- c(
  "# Figure legends",
  "",
  "Figure 1. Regime definition and shared knot estimation. Panels A-B show the median observed cumulative PC1-NB divergence and shared two-knot fit for CPM-EVS and DESeq2-EVS. Numeric c1 and c2 values are printed directly on the panels. Panels C-D summarize regime definitions and the shared-knot fitting objective.",
  "",
  "Figure 2. CPM-EVS dataset-specific Pareto cutoffs. Each panel corresponds to one RT/ZT comparison. The faint trajectory shows all candidate k values, the purple curve is the Pareto frontier, and the highlighted diamond is the geometric-knee k*. Each panel reports c1, c2, k*, weighted benefit, weighted penalty, Joint sites, high-confidence sites, and Remainder crossings.",
  "",
  "Figure 3. DESeq2-EVS dataset-specific Pareto cutoffs. Layout and annotations are identical to Figure 2.",
  "",
  "Figure 4. Cross-method summary. Panel A compares all eight dataset-specific k* values. Panel B compares high-confidence EVS membership and overlap. Panel C reports penalized Remainder crossings at each k*. Panel D reports k* as a percentage of the available Leading-Edge candidate domain."
)

writeLines(
  legend_lines,
  file.path(OUT_ROOT,"Figure_Legends.md")
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
