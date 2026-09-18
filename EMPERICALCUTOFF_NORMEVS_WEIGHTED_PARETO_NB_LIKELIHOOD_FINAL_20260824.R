#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# CPM-EVS vs DESeq2-EVS — FAST FINAL
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT  <- "/root/REAPER98632/exports/cpm_vs_deseq2_evs_fast_final"

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
  control="#0072B2", treatment="#009E73",
  observed="#6A3D9A", fit="#222222",
  c1="#D73027", c2="#1A9850", selected="#B5179E",
  remainder="#DCE6F2", divergence="#FFF0B3", leading="#D8F3E7"
)

dir.create(OUT_ROOT, recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# INPUT / NORMALIZATION
# =============================================================================

read_counts <- function(path) {
  x <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)

  idx <- sort(unique(unlist(lapply(GROUP_PATTERNS, \(p) grep(p, names(x))))))
  if (!length(idx)) stop("No sample columns matched GROUP_PATTERNS.")

  ids <- trimws(as.character(x[[1]]))
  ids[is.na(ids) | ids==""] <- paste0("feature_", which(is.na(ids) | ids==""))
  ids <- make.unique(ids, sep="__dup_")

  m <- do.call(cbind, lapply(x[, idx, drop=FALSE], \(z)
    suppressWarnings(as.numeric(as.character(z)))
  ))

  rownames(m) <- ids
  colnames(m) <- names(x)[idx]
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  m <- pmax(m, 0)

  # One experiment-wide PAS universe.
  m[rowSums(m) > 0, , drop=FALSE]
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
  lib[!is.finite(lib) | lib <= 0] <- 1
  log1p(sweep(x, 2, lib/1e6, "/"))
}

pooled_within_group_var <- function(norm, groups) {
  sse <- rep(0, nrow(norm))
  df <- 0L

  for (g in levels(groups)) {
    z <- norm[, groups==g, drop=FALSE]
    mu <- rowMeans(z)
    sse <- sse + rowSums((z-mu)^2)
    df <- df + ncol(z)-1L
  }

  pmax(sse/df, 0)
}

# =============================================================================
# PC1 / DIVERGENCE
# =============================================================================

pc1_rank <- function(x) {
  p <- prcomp(t(x), center=TRUE, scale.=FALSE, rank.=1)
  loading <- p$rotation[,1]
  loading[!is.finite(loading)] <- 0

  list(
    order=order(abs(loading), decreasing=FALSE),
    loading=loading,
    contribution=(p$sdev[1]^2) * loading^2
  )
}

build_arm <- function(method, arm, raw, norm, pooled_var) {
  x <- if (method=="CPM_EVS") cpm_log1p(raw) else log1p(norm)

  pc <- pc1_rank(x)
  ord <- pc$order

  P <- pc$contribution[ord]
  mu <- rowMeans(norm)[ord]
  E <- pmax(pooled_var[ord] - mu, 0)

  if (sum(P) <= 0 || sum(E) <= 0) {
    stop("Undefined variance mass: ", method, " / ", arm)
  }

  p <- P/sum(P)
  q <- E/sum(E)

  data.frame(
    method=method,
    arm=arm,
    rank=seq_along(ord),
    feature_id=rownames(raw)[ord],
    abs_pc1_loading=abs(pc$loading[ord]),
    pc1_mass=p,
    excess_mass=q,
    F_P=cumsum(p),
    F_E=cumsum(q),
    D=cumsum(q)-cumsum(p),
    stringsAsFactors=FALSE
  )
}

piecewise_basis <- function(x, c1, c2) {
  cbind(1, x, pmax(x-c1,0), pmax(x-c2,0))
}

fit_shared_knots <- function(arms) {
  N <- nrow(arms[[1]])
  x_full <- (seq_len(N)-1)/(N-1)
  D_full <- do.call(cbind, lapply(arms, \(z) z$D))

  # FAST: optimize on at most 4,000 evenly spaced ranks.
  opt_n <- min(4000L, N)
  idx <- unique(round(seq(1, N, length.out=opt_n)))
  x <- x_full[idx]
  D <- D_full[idx,,drop=FALSE]

  min_gap <- max(4/(N-1), 1e-4)

  objective <- function(par) {
    c1 <- par[1]; c2 <- par[2]

    if (!is.finite(c1) || !is.finite(c2) ||
        c1 <= 0 || c2 >= 1 || c2-c1 <= min_gap) return(1e100)

    X <- piecewise_basis(x, c1, c2)
    b <- tryCatch(qr.coef(qr(X), D), error=\(e) NULL)

    if (is.null(b) || any(!is.finite(b))) return(1e100)
    sum((D - X %*% b)^2)
  }

  starts <- list(c(.05,.95), c(.12,.88), c(.22,.78), c(.32,.68))

  fits <- lapply(starts, \(st)
    optim(st, objective, method="Nelder-Mead",
          control=list(maxit=500, reltol=1e-10))
  )

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]

  c1 <- as.integer(round(1 + best$par[1]*(N-1)))
  c2 <- as.integer(round(1 + best$par[2]*(N-1)))

  c1 <- max(2L, min(N-2L, c1))
  c2 <- max(c1+1L, min(N-1L, c2))

  # Full matrix is fit ONCE after knot selection.
  X_full <- piecewise_basis(
    x_full,
    (c1-1)/(N-1),
    (c2-1)/(N-1)
  )

  coef <- qr.coef(qr(X_full), D_full)
  fitted <- X_full %*% coef

  list(
    c1=c1, c2=c2, N=N,
    fitted=fitted,
    arms=names(arms),
    SSE=sum((D_full-fitted)^2)
  )
}

# =============================================================================
# FAST CONTINUOUS WEIGHTED TOP-k SCAN
# =============================================================================

rank_map <- function(df) setNames(df$rank, df$feature_id)

div_support <- function(r, c1, c2) {
  pmin(pmax((r-c1)/(c2-c1), 0), 1)
}

lead_support <- function(r, c2, N) {
  pmin(pmax((r-c2)/(N-c2), 0), 1)
}

rem_depth <- function(r, c1) {
  pmin(pmax((c1-r)/(c1-1), 0), 1)
}

rank_gap <- function(a, b, N) abs(a-b)/(N-1)

add_events <- function(pos, weight, L) {
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

activate_from <- function(depth, weight, K) {
  cumsum(add_events(depth, weight, K))
}

active_interval <- function(start, end, weight, K) {
  # Active for start <= k < end. end may equal K+1.
  delta <- numeric(K+1L)

  ok <- is.finite(start) & is.finite(end) & is.finite(weight) &
        start>=1 & start<=K & end>start

  if (!any(ok)) return(numeric(K))

  s <- as.integer(start[ok])
  e <- pmin(as.integer(end[ok]), K+1L)
  w <- weight[ok]

  delta <- delta + add_events(s, w, K+1L)
  delta <- delta - add_events(e, w, K+1L)

  cumsum(delta)[seq_len(K)]
}

scan_pair_fast <- function(control_df, treatment_df, c1, c2) {
  N <- nrow(control_df)
  K <- N-c2

  if (K < 1L) stop("No candidate k beyond c2.")

  rC_map <- rank_map(control_df)
  rT_map <- rank_map(treatment_df)

  ids <- control_df$feature_id
  rC <- as.integer(rC_map[ids])
  rT <- as.integer(rT_map[ids])

  dC <- N-rC+1L
  dT <- N-rT+1L

  # Joint starts once both top-k sets contain the PAS.
  joint_depth <- pmax(dC,dT)
  joint_n <- activate_from(joint_depth, rep(1,N), K)

  # Control-only intervals.
  idxC <- which(dC < dT & dC <= K)
  startC <- dC[idxC]
  endC <- pmin(dT[idxC], K+1L)
  oppC <- rT[idxC]
  selC <- rC[idxC]

  supportC <- ifelse(
    oppC > c2,
    lead_support(oppC,c2,N),
    ifelse(oppC >= c1, div_support(oppC,c1,c2), 0)
  )

  penaltyC <- ifelse(
    oppC < c1,
    rem_depth(oppC,c1)^2 * rank_gap(selC,oppC,N),
    0
  )

  disC_n <- active_interval(startC,endC,rep(1,length(idxC)),K)
  disC_support <- active_interval(startC,endC,supportC,K)
  disC_penalty <- active_interval(startC,endC,penaltyC,K)
  disC_rem <- active_interval(startC,endC,as.numeric(oppC<c1),K)
  disC_gap <- active_interval(startC,endC,rank_gap(selC,oppC,N),K)

  # Treatment-only intervals.
  idxT <- which(dT < dC & dT <= K)
  startT <- dT[idxT]
  endT <- pmin(dC[idxT], K+1L)
  oppT <- rC[idxT]
  selT <- rT[idxT]

  supportT <- ifelse(
    oppT > c2,
    lead_support(oppT,c2,N),
    ifelse(oppT >= c1, div_support(oppT,c1,c2), 0)
  )

  penaltyT <- ifelse(
    oppT < c1,
    rem_depth(oppT,c1)^2 * rank_gap(selT,oppT,N),
    0
  )

  disT_n <- active_interval(startT,endT,rep(1,length(idxT)),K)
  disT_support <- active_interval(startT,endT,supportT,K)
  disT_penalty <- active_interval(startT,endT,penaltyT,K)
  disT_rem <- active_interval(startT,endT,as.numeric(oppT<c1),K)
  disT_gap <- active_interval(startT,endT,rank_gap(selT,oppT,N),K)

  disjoint_n <- disC_n + disT_n
  weighted_benefit <- joint_n + disC_support + disT_support
  weighted_penalty <- disC_penalty + disT_penalty
  remainder_cross_n <- disC_rem + disT_rem
  mean_rank_gap <- ifelse(
    disjoint_n > 0,
    (disC_gap + disT_gap)/disjoint_n,
    0
  )

  k <- seq_len(K)

  data.frame(
    k=k,
    joint_n=joint_n,
    disjoint_n=disjoint_n,
    union_n=2*k-joint_n,
    weighted_benefit=weighted_benefit,
    weighted_penalty=weighted_penalty,
    remainder_cross_n=remainder_cross_n,
    mean_rank_gap=mean_rank_gap,
    stringsAsFactors=FALSE
  )
}

# =============================================================================
# GLOBAL METHOD-SPECIFIC CUTOFF
# =============================================================================

norm01 <- function(x) {
  r <- range(x, na.rm=TRUE)
  if (!is.finite(diff(r)) || diff(r)==0) return(rep(0,length(x)))
  (x-r[1])/diff(r)
}

pareto_frontier <- function(df) {
  x <- df %>%
    arrange(weighted_penalty, desc(weighted_benefit), desc(k))

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

  x1 <- x[1]; y1 <- y[1]
  x2 <- x[length(x)]; y2 <- y[length(y)]

  den <- sqrt((y2-y1)^2 + (x2-x1)^2)

  if (!is.finite(den) || den==0) {
    return(frontier$k[which.max(y-x)])
  }

  d <- abs((y2-y1)*x - (x2-x1)*y + x2*y1 - y2*x1)/den
  frontier$k[which.max(d)]
}

aggregate_scans <- function(scans) {
  bind_rows(lapply(names(scans), \(nm)
    mutate(scans[[nm]], comparison=nm)
  )) %>%
    group_by(k) %>%
    summarise(
      weighted_benefit=sum(weighted_benefit),
      weighted_penalty=sum(weighted_penalty),
      joint_n=sum(joint_n),
      disjoint_n=sum(disjoint_n),
      union_n=sum(union_n),
      remainder_cross_n=sum(remainder_cross_n),
      mean_rank_gap=mean(mean_rank_gap),
      .groups="drop"
    )
}

# =============================================================================
# CLASSIFICATION AT SHARED METHOD-SPECIFIC k*
# =============================================================================

classify_at_k <- function(k, comparison, control_df, treatment_df, c1, c2) {
  rC <- rank_map(control_df)
  rT <- rank_map(treatment_df)

  topC <- tail(control_df$feature_id,k)
  topT <- tail(treatment_df$feature_id,k)
  ids <- union(topC,topT)

  a <- as.integer(rC[ids])
  b <- as.integer(rT[ids])

  inC <- ids %in% topC
  inT <- ids %in% topT
  joint <- inC & inT
  disC <- inC & !inT
  disT <- inT & !inC

  opp <- ifelse(disC,b,ifelse(disT,a,NA_integer_))

  opp_region <- ifelse(
    joint, "Joint",
    ifelse(opp>c2, "Opposite Leading Edge",
           ifelse(opp>=c1, "Opposite Divergence", "Opposite Remainder"))
  )

  retained <- joint | opp_region %in%
    c("Opposite Leading Edge","Opposite Divergence")

  data.frame(
    comparison=comparison,
    feature_id=ids,
    selected_k=k,
    control_rank=a,
    treatment_rank=b,
    class=ifelse(joint,"Joint",ifelse(disC,"Control only","Treatment only")),
    opposite_region=opp_region,
    retained=retained,
    stringsAsFactors=FALSE
  )
}

classification_counts <- function(x) {
  x %>%
    mutate(category=case_when(
      opposite_region=="Joint" ~ "Joint",
      opposite_region=="Opposite Leading Edge" ~ "Disjoint: opposite LE",
      opposite_region=="Opposite Divergence" ~ "Disjoint: opposite Divergence",
      TRUE ~ "Excluded: opposite Remainder"
    )) %>%
    count(category, name="n") %>%
    complete(
      category=c(
        "Joint",
        "Disjoint: opposite LE",
        "Disjoint: opposite Divergence",
        "Excluded: opposite Remainder"
      ),
      fill=list(n=0)
    )
}

# =============================================================================
# FIGURES — MINIMAL MANUSCRIPT SET
# =============================================================================

theme_pub <- function() {
  theme_classic(base_size=12) +
    theme(
      plot.title=element_text(face="bold", size=13),
      plot.subtitle=element_text(size=10.5, color="grey30"),
      axis.title=element_text(face="bold"),
      axis.text=element_text(color="grey20"),
      panel.border=element_rect(fill=NA, color="grey78", linewidth=.4),
      legend.position="none",
      plot.margin=margin(8,10,8,10)
    )
}

save_plot <- function(p,file,w=10,h=6.5) {
  ggsave(file,p,width=w,height=h,dpi=PNG_DPI,bg="white")
  ggsave(sub("\\.png$",".pdf",file),p,width=w,height=h,bg="white")
}

plot_regime <- function(method, arms, knot, file) {
  obs <- bind_rows(lapply(arms, \(z) select(z,rank,D)))

  med <- obs %>%
    group_by(rank) %>%
    summarise(D=median(D),.groups="drop")

  fit <- data.frame(
    rank=seq_len(knot$N),
    D=apply(knot$fitted,1,median)
  )

  p <- ggplot() +
    annotate(
      "rect", xmin=1, xmax=knot$c1, ymin=-Inf, ymax=Inf,
      fill=COL$remainder, alpha=.35
    ) +
    annotate(
      "rect", xmin=knot$c1, xmax=knot$c2, ymin=-Inf, ymax=Inf,
      fill=COL$divergence, alpha=.35
    ) +
    annotate(
      "rect", xmin=knot$c2, xmax=knot$N, ymin=-Inf, ymax=Inf,
      fill=COL$leading, alpha=.35
    ) +
    geom_hline(yintercept=0,color="grey70",linetype="dotted") +
    geom_line(data=med,aes(rank,D),color=COL$observed,linewidth=1.15) +
    geom_line(data=fit,aes(rank,D),color=COL$fit,linewidth=1.15) +
    geom_vline(
      xintercept=knot$c1,
      color=COL$c1,linetype="dashed",linewidth=.8
    ) +
    geom_vline(
      xintercept=knot$c2,
      color=COL$c2,linetype="longdash",linewidth=.8
    ) +
    annotate(
      "label",x=knot$c1,y=Inf,
      label=paste0("c1 = ",knot$c1),
      hjust=1.05,vjust=1.2,size=3,
      fill="white",label.size=.15
    ) +
    annotate(
      "label",x=knot$c2,y=Inf,
      label=paste0("c2 = ",knot$c2),
      hjust=-.05,vjust=1.2,size=3,
      fill="white",label.size=.15
    ) +
    labs(
      title=paste0(method, ": shared divergence geometry"),
      subtitle="Median observed D(r) and shared two-knot fit across all eight arms",
      x="PC1 rank",
      y="D(r) = F_E(r) - F_P(r)"
    ) +
    theme_pub()

  save_plot(p,file,11.5,6.3)
}

plot_pareto <- function(method, comparison, scan, frontier, kstar, file) {
  selected <- scan[scan$k==kstar,,drop=FALSE]

  # Full candidate path is faint; Pareto frontier is the emphasized curve.
  p <- ggplot() +
    geom_path(
      data=scan,
      aes(weighted_penalty, weighted_benefit),
      color="grey82",
      linewidth=.55
    ) +
    geom_path(
      data=frontier,
      aes(weighted_penalty, weighted_benefit),
      color=COL$observed,
      linewidth=1.35
    ) +
    geom_point(
      data=selected,
      aes(weighted_penalty, weighted_benefit),
      shape=23,
      size=5,
      fill=COL$selected,
      color=COL$selected
    ) +
    annotate(
      "label",
      x=selected$weighted_penalty,
      y=selected$weighted_benefit,
      label=paste0(
        "k* = ",kstar,
        "\nJoint = ",round(selected$joint_n),
        "\nRemainder crossings = ",round(selected$remainder_cross_n),
        "\nWeighted penalty = ",formatC(selected$weighted_penalty,digits=2,format="f")
      ),
      hjust=-.05,
      vjust=1.05,
      size=3.1,
      fill="white",
      label.size=.15
    ) +
    labs(
      title=paste0(method," — ",comparison,": weighted Pareto cutoff"),
      subtitle="Purple: Pareto frontier; diamond: selected geometric knee",
      x="Distance-weighted Remainder disagreement",
      y="Weighted retained-site benefit"
    ) +
    theme_pub()

  save_plot(p,file,10.5,7)
}

plot_method_comparison <- function(summary,file) {
  summary <- summary %>%
    mutate(
      method=factor(
        method,
        levels=c("CPM_EVS","DESeq2_EVS"),
        labels=c("CPM-EVS","DESeq2-EVS")
      ),
      comparison=factor(comparison,levels=names(COMPARISONS))
    )

  p <- ggplot(
    summary,
    aes(comparison,selected_k,fill=method)
  ) +
    geom_col(
      position=position_dodge(width=.74),
      width=.64
    ) +
    geom_text(
      aes(label=selected_k),
      position=position_dodge(width=.74),
      vjust=-.4,
      fontface="bold",
      size=3.4
    ) +
    labs(
      title="Empirical EVS cutoffs by dataset and normalization method",
      subtitle="Each RT/ZT dataset has an independent Pareto-knee cutoff under CPM-EVS and DESeq2-EVS",
      x=NULL,
      y="Selected k*",
      fill=NULL
    ) +
    theme_pub() +
    theme(legend.position="top") +
    expand_limits(y=max(summary$selected_k,na.rm=TRUE)*1.14)

  save_plot(p,file,11.5,6.5)
}

# =============================================================================
# RUN
# =============================================================================

message("Reading counts...")
counts <- read_counts(COUNT_FILE)
groups <- assign_groups(colnames(counts))

message("Features: ",nrow(counts)," | samples: ",ncol(counts))
message("Global DESeq2 normalization...")
norm_obj <- deseq2_normalize(counts,groups)
norm <- norm_obj$counts

message("Pooled normalized within-group variance...")
pooled_var <- pooled_within_group_var(norm,groups)

summary_list <- list()
sites_list <- list()

for (method in METHODS) {

  message("============================================================")
  message("METHOD: ",method)

  mdir <- file.path(OUT_ROOT,method)
  fdir <- file.path(mdir,"Figures")
  tdir <- file.path(mdir,"Tables")
  dir.create(fdir,recursive=TRUE,showWarnings=FALSE)
  dir.create(tdir,recursive=TRUE,showWarnings=FALSE)

  message("Building 8 arm-specific PC1 divergence profiles...")
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

  message("Fast shared-knot fit...")
  knot <- fit_shared_knots(arms)
  K <- knot$N-knot$c2

  message(
    "c1=",knot$c1,
    " | c2=",knot$c2,
    " | candidate k=1..",K
  )

  plot_regime(
    method,arms,knot,
    file.path(fdir,paste0("Figure_",method,"_Shared_Divergence.png"))
  )

  message("Vectorized pairwise cutoff scans with independent k* per dataset...")
  scans <- list()
  cls_method <- list()
  pair_summary <- list()

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

    scans[[nm]] <- scan

    message("  ", nm, " k* = ", kstar)

    write.csv(
      scan,
      file.path(tdir,paste0("Table_",method,"_",nm,"_Cutoff_Scan.csv")),
      row.names=FALSE
    )

    write.csv(
      frontier,
      file.path(tdir,paste0("Table_",method,"_",nm,"_Pareto_Frontier.csv")),
      row.names=FALSE
    )

    plot_pareto(
      method=method,
      comparison=nm,
      scan=scan,
      frontier=frontier,
      kstar=kstar,
      file=file.path(
        fdir,
        paste0("Figure_",method,"_",nm,"_Pareto_Cutoff.png")
      )
    )

    cls <- classify_at_k(
      kstar,
      nm,
      arms[[mp[["control"]]]],
      arms[[mp[["treatment"]]]],
      knot$c1,
      knot$c2
    )

    cls$method <- method
    cls_method[[nm]] <- cls

    write.csv(
      cls,
      file.path(tdir,paste0("Table_",method,"_",nm,"_Sites.csv")),
      row.names=FALSE
    )

    sel <- scan[scan$k==kstar,,drop=FALSE]

    pair_summary[[nm]] <- data.frame(
      method=method,
      comparison=nm,
      N=knot$N,
      c1=knot$c1,
      c2=knot$c2,
      divergence_size=knot$c2-knot$c1+1L,
      leading_edge_candidate_size=K,
      selected_k=kstar,
      k_over_leading_edge=kstar/K,
      weighted_benefit=sel$weighted_benefit,
      weighted_penalty=sel$weighted_penalty,
      joint_n=sel$joint_n,
      remainder_cross_n=sel$remainder_cross_n,
      knot_SSE=knot$SSE,
      stringsAsFactors=FALSE
    )
  }

  pair_summary_df <- bind_rows(pair_summary)

  write.csv(
    pair_summary_df,
    file.path(tdir,paste0("Table_",method,"_Dataset_Cutoffs.csv")),
    row.names=FALSE
  )

  summary_list[[method]] <- pair_summary_df
  sites_list[[method]] <- bind_rows(cls_method)
}

summary <- bind_rows(summary_list)
all_sites <- bind_rows(sites_list)

write.csv(
  summary,
  file.path(OUT_ROOT,"Table_CPM_vs_DESeq2_EVS_Dataset_Cutoffs.csv"),
  row.names=FALSE
)

write.csv(
  all_sites,
  file.path(OUT_ROOT,"Table_CPM_vs_DESeq2_EVS_All_Sites.csv"),
  row.names=FALSE
)

overlap <- bind_rows(lapply(names(COMPARISONS), \(nm) {
  a <- all_sites %>%
    filter(method=="CPM_EVS",comparison==nm,retained) %>%
    pull(feature_id) %>% unique()

  b <- all_sites %>%
    filter(method=="DESeq2_EVS",comparison==nm,retained) %>%
    pull(feature_id) %>% unique()

  u <- union(a,b)

  data.frame(
    comparison=nm,
    CPM_retained=length(a),
    DESeq2_retained=length(b),
    overlap=length(intersect(a,b)),
    union=length(u),
    jaccard=ifelse(length(u),length(intersect(a,b))/length(u),NA_real_)
  )
}))

write.csv(
  overlap,
  file.path(OUT_ROOT,"Table_CPM_vs_DESeq2_EVS_Overlap.csv"),
  row.names=FALSE
)

plot_method_comparison(
  summary,
  file.path(OUT_ROOT,"Figure_CPM_vs_DESeq2_EVS_Cutoff_Comparison_Bar.png")
)


# =============================================================================
# FINAL MANUSCRIPT EXPORTS
# =============================================================================

zip_directory_files <- function(zipfile, files, root_dir) {
  files <- unique(files[file.exists(files)])
  if (!length(files)) {
    warning("No files found for zip: ", zipfile)
    return(FALSE)
  }

  if (file.exists(zipfile)) unlink(zipfile)

  old <- getwd()
  on.exit(setwd(old), add=TRUE)
  setwd(root_dir)

  rel <- sub(
    paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", normalizePath(root_dir)), "/?"),
    "",
    normalizePath(files)
  )

  utils::zip(zipfile=zipfile, files=rel)
  file.exists(zipfile)
}

collect_files <- function(root, pattern) {
  list.files(
    root,
    pattern=pattern,
    recursive=TRUE,
    full.names=TRUE,
    ignore.case=TRUE
  )
}

write_manifest <- function(files, type, out_file) {
  if (!length(files)) {
    write.csv(
      data.frame(type=character(), file=character()),
      out_file,
      row.names=FALSE
    )
    return(invisible(NULL))
  }

  root_norm <- normalizePath(OUT_ROOT)
  f_norm <- normalizePath(files)

  rel <- sub(
    paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", root_norm), "/?"),
    "",
    f_norm
  )

  out <- data.frame(
    type=type,
    file=rel,
    stringsAsFactors=FALSE
  )

  write.csv(out, out_file, row.names=FALSE)
}

# Collect final manuscript-ready outputs.
all_figure_files <- collect_files(
  OUT_ROOT,
  "\\.(png|pdf)$"
)

all_table_files <- collect_files(
  OUT_ROOT,
  "\\.csv$"
)

# Exclude manifests themselves from recursive inclusion until written.
all_table_files <- all_table_files[
  !grepl("Manifest", basename(all_table_files), ignore.case=TRUE)
]

figure_manifest <- file.path(
  OUT_ROOT,
  "Manifest_Manuscript_Figures.csv"
)

table_manifest <- file.path(
  OUT_ROOT,
  "Manifest_Manuscript_Tables.csv"
)

write_manifest(
  all_figure_files,
  "figure",
  figure_manifest
)

write_manifest(
  all_table_files,
  "table",
  table_manifest
)

# Add manifests to table bundle.
all_table_files <- c(
  all_table_files,
  figure_manifest,
  table_manifest
)

FIGURE_ZIP <- file.path(
  OUT_ROOT,
  "Manuscript_Ready_Figures.zip"
)

TABLE_ZIP <- file.path(
  OUT_ROOT,
  "Manuscript_Ready_Tables.zip"
)

COMPLETE_ZIP <- file.path(
  OUT_ROOT,
  "CPM_vs_DESeq2_EVS_Complete_Outputs.zip"
)

message("Creating final manuscript-ready ZIP files...")

fig_zip_ok <- zip_directory_files(
  FIGURE_ZIP,
  all_figure_files,
  OUT_ROOT
)

tab_zip_ok <- zip_directory_files(
  TABLE_ZIP,
  all_table_files,
  OUT_ROOT
)

complete_files <- unique(c(
  all_figure_files,
  all_table_files
))

complete_zip_ok <- zip_directory_files(
  COMPLETE_ZIP,
  complete_files,
  OUT_ROOT
)

if (!fig_zip_ok) warning("Figure ZIP was not created.")
if (!tab_zip_ok) warning("Table ZIP was not created.")
if (!complete_zip_ok) warning("Complete output ZIP was not created.")

message("------------------------------------------------------------")
message("FINAL EXPORTS")
message("Figures ZIP: ", FIGURE_ZIP)
message("Tables ZIP: ", TABLE_ZIP)
message("Complete ZIP: ", COMPLETE_ZIP)
message("------------------------------------------------------------")

message("============================================================")
message("COMPLETE")
message("Output: ",OUT_ROOT)
message("Primary figures are shared divergence plots, dataset-specific Pareto curves, and one cutoff comparison bar graph;")
message("shared divergence geometry is retained only where a curve is required.")
message("Figures ZIP: ", FIGURE_ZIP)
message("Tables ZIP: ", TABLE_ZIP)
message("Complete ZIP: ", COMPLETE_ZIP)
print(summary)
message("============================================================")
