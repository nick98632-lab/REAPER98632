#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# CPM-EVS vs DESeq2-EVS
# Shared 8-arm regime geometry + continuous rank-weighted cutoff selection
# =============================================================================

COUNT_FILE <- "/root/REAPER98632/data/WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
OUT_ROOT <- "/root/REAPER98632/exports/cpm_vs_deseq2_evs_refined"

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
  benefit="#0072B2", cost="#D55E00",
  remainder="#EAF0F6", divergence="#FFF3BF", leading="#E4F5EC"
)

dir.create(OUT_ROOT, recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# I/O + NORMALIZATION
# =============================================================================

read_counts <- function(path) {
  x <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  sample_idx <- sort(unique(unlist(lapply(GROUP_PATTERNS, \(p) grep(p, names(x))))))
  if (!length(sample_idx)) stop("No sample columns matched GROUP_PATTERNS.")
  ids <- make.unique(trimws(as.character(x[[1]])), sep="__dup_")
  m <- as.matrix(data.frame(lapply(x[, sample_idx, drop=FALSE], as.numeric)))
  rownames(m) <- ids
  colnames(m) <- names(x)[sample_idx]
  m[!is.finite(m)] <- 0
  m <- pmax(m, 0)
  m[rowSums(m) > 0, , drop=FALSE]
}

assign_groups <- function(samples) {
  g <- rep(NA_character_, length(samples))
  names(g) <- samples
  for (nm in names(GROUP_PATTERNS)) {
    idx <- grep(GROUP_PATTERNS[[nm]], samples)
    if (any(!is.na(g[idx]))) stop("Sample matched >1 group.")
    g[idx] <- nm
  }
  if (anyNA(g)) stop("Unassigned samples: ", paste(names(g)[is.na(g)], collapse=", "))
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
  log1p(sweep(x, 2, lib / 1e6, "/"))
}

pooled_within_group_var <- function(norm, groups) {
  sse <- rep(0, nrow(norm))
  df <- 0L
  for (g in levels(groups)) {
    idx <- which(groups == g)
    z <- norm[, idx, drop=FALSE]
    mu <- rowMeans(z)
    sse <- sse + rowSums((z - mu)^2)
    df <- df + ncol(z) - 1L
  }
  pmax(sse / df, 0)
}

# =============================================================================
# PC1 + DIVERGENCE
# =============================================================================

pc1_rank <- function(x) {
  p <- prcomp(t(x), center=TRUE, scale.=FALSE, rank.=1)
  loading <- p$rotation[,1]
  loading[!is.finite(loading)] <- 0
  ord <- order(abs(loading), decreasing=FALSE)
  lambda1 <- p$sdev[1]^2
  list(
    order=ord,
    loading=loading,
    contribution=lambda1 * loading^2
  )
}

build_arm <- function(method, arm, raw, norm, pooled_var) {
  x <- switch(
    method,
    CPM_EVS=cpm_log1p(raw),
    DESeq2_EVS=log1p(norm),
    stop("Unknown method")
  )

  pc <- pc1_rank(x)
  ord <- pc$order

  mu <- rowMeans(norm)
  P <- pc$contribution[ord]
  E <- pmax(pooled_var[ord] - mu[ord], 0)

  if (sum(P) <= 0 || sum(E) <= 0) stop("Undefined variance mass: ", method, " / ", arm)

  p <- P / sum(P)
  q <- E / sum(E)

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
    D=cumsum(q) - cumsum(p),
    stringsAsFactors=FALSE
  )
}

piecewise_basis <- function(x, c1, c2) {
  cbind(1, x, pmax(x-c1,0), pmax(x-c2,0))
}

fit_shared_knots <- function(arms) {
  N <- nrow(arms[[1]])
  x <- (seq_len(N)-1)/(N-1)
  D <- do.call(cbind, lapply(arms, \(z) z$D))

  sse <- function(par, x, D) {
    c1 <- par[1]; c2 <- par[2]
    if (!is.finite(c1) || !is.finite(c2) || c1 <= 0 || c2 >= 1 || c2 <= c1)
      return(1e100)
    X <- piecewise_basis(x, c1, c2)
    b <- tryCatch(qr.coef(qr(X), D), error=\(e) NULL)
    if (is.null(b) || any(!is.finite(b))) return(1e100)
    sum((D - X %*% b)^2)
  }

  starts <- list(c(.08,.92), c(.15,.85), c(.25,.75), c(.35,.65))
  fits <- lapply(starts, \(st) optim(st, sse, x=x, D=D, method="Nelder-Mead"))
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  fit <- optim(best$par, sse, x=x, D=D, method="Nelder-Mead",
               control=list(maxit=1500, reltol=1e-12))

  c1 <- round(1 + fit$par[1]*(N-1))
  c2 <- round(1 + fit$par[2]*(N-1))
  c1 <- max(2L, min(N-2L, c1))
  c2 <- max(c1+1L, min(N-1L, c2))

  X <- piecewise_basis(x, (c1-1)/(N-1), (c2-1)/(N-1))
  coef <- qr.coef(qr(X), D)

  list(
    c1=c1, c2=c2, N=N,
    fitted=X %*% coef,
    arms=names(arms),
    SSE=sum((D - X %*% coef)^2)
  )
}

# =============================================================================
# CONTINUOUS AGREEMENT / DISAGREEMENT SCORING
# =============================================================================

rank_map <- function(df) setNames(df$rank, df$feature_id)

divergence_support <- function(r, c1, c2) {
  pmin(pmax((r-c1)/(c2-c1), 0), 1)
}

leading_support <- function(r, c2, N) {
  pmin(pmax((r-c2)/(N-c2), 0), 1)
}

remainder_depth <- function(r, c1) {
  pmin(pmax((c1-r)/(c1-1), 0), 1)
}

rank_gap <- function(r1, r2, N) abs(r1-r2)/(N-1)

score_pair_at_k <- function(k, control_df, treatment_df, c1, c2) {
  N <- nrow(control_df)
  rC <- rank_map(control_df)
  rT <- rank_map(treatment_df)

  topC <- tail(control_df$feature_id, k)
  topT <- tail(treatment_df$feature_id, k)
  ids <- union(topC, topT)

  a <- as.integer(rC[ids])
  b <- as.integer(rT[ids])
  inC <- ids %in% topC
  inT <- ids %in% topT
  joint <- inC & inT
  disC <- inC & !inT
  disT <- inT & !inC

  # Joint sites receive full support.
  joint_score <- sum(joint)

  # Disjoint sites receive continuous support according to opposite-arm rank.
  opp_rank <- ifelse(disC, b, ifelse(disT, a, NA_integer_))
  selected_rank <- ifelse(disC, a, ifelse(disT, b, NA_integer_))
  is_dis <- disC | disT

  opp_support <- ifelse(
    opp_rank > c2,
    leading_support(opp_rank, c2, N),
    ifelse(
      opp_rank >= c1,
      divergence_support(opp_rank, c1, c2),
      0
    )
  )

  disjoint_score <- sum(opp_support[is_dis], na.rm=TRUE)

  # Deep opposite-arm Remainder crossings are penalized more strongly,
  # especially when the two arm ranks are far apart.
  penalty <- ifelse(
    is_dis & opp_rank < c1,
    remainder_depth(opp_rank, c1)^2 *
      rank_gap(selected_rank, opp_rank, N),
    0
  )

  data.frame(
    k=k,
    joint_n=sum(joint),
    union_n=length(ids),
    disjoint_n=sum(is_dis),
    weighted_benefit=joint_score + disjoint_score,
    weighted_penalty=sum(penalty, na.rm=TRUE),
    remainder_cross_n=sum(is_dis & opp_rank < c1, na.rm=TRUE),
    mean_rank_gap=if (any(is_dis)) mean(rank_gap(selected_rank[is_dis], opp_rank[is_dis], N)) else 0,
    stringsAsFactors=FALSE
  )
}

scan_pair <- function(control_df, treatment_df, c1, c2) {
  N <- nrow(control_df)
  K <- N-c2
  bind_rows(lapply(seq_len(K), \(k) score_pair_at_k(k, control_df, treatment_df, c1, c2)))
}

# =============================================================================
# GLOBAL METHOD-SPECIFIC CUTOFF
# =============================================================================

normalize01 <- function(x) {
  r <- range(x, na.rm=TRUE)
  if (!is.finite(diff(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x-r[1])/diff(r)
}

pareto_frontier <- function(df, benefit="weighted_benefit", cost="weighted_penalty") {
  x <- df %>%
    arrange(.data[[cost]], desc(.data[[benefit]]), desc(k))

  keep <- logical(nrow(x))
  best <- -Inf
  for (i in seq_len(nrow(x))) {
    if (x[[benefit]][i] > best) {
      keep[i] <- TRUE
      best <- x[[benefit]][i]
    }
  }
  x[keep, , drop=FALSE] %>% arrange(.data[[cost]])
}

pareto_knee <- function(frontier) {
  if (nrow(frontier) == 1) return(frontier$k[1])

  x <- normalize01(frontier$weighted_penalty)
  y <- normalize01(frontier$weighted_benefit)

  # Distance from the chord joining frontier endpoints.
  x1 <- x[1]; y1 <- y[1]
  x2 <- x[length(x)]; y2 <- y[length(y)]
  den <- sqrt((y2-y1)^2 + (x2-x1)^2)

  if (!is.finite(den) || den == 0) {
    return(frontier$k[which.max(y-x)])
  }

  d <- abs((y2-y1)*x - (x2-x1)*y + x2*y1 - y2*x1) / den
  frontier$k[which.max(d)]
}

aggregate_scans <- function(scans) {
  bind_rows(lapply(names(scans), \(nm) mutate(scans[[nm]], comparison=nm))) %>%
    group_by(k) %>%
    summarise(
      weighted_benefit=sum(weighted_benefit),
      weighted_penalty=sum(weighted_penalty),
      joint_n=sum(joint_n),
      union_n=sum(union_n),
      disjoint_n=sum(disjoint_n),
      remainder_cross_n=sum(remainder_cross_n),
      mean_rank_gap=mean(mean_rank_gap),
      .groups="drop"
    )
}

# =============================================================================
# FINAL SITE CLASSIFICATION
# =============================================================================

classify_at_k <- function(k, comparison, control_df, treatment_df, c1, c2) {
  N <- nrow(control_df)
  rC <- rank_map(control_df)
  rT <- rank_map(treatment_df)

  topC <- tail(control_df$feature_id, k)
  topT <- tail(treatment_df$feature_id, k)
  ids <- union(topC, topT)

  a <- as.integer(rC[ids])
  b <- as.integer(rT[ids])
  inC <- ids %in% topC
  inT <- ids %in% topT
  joint <- inC & inT
  disC <- inC & !inT
  disT <- inT & !inC
  opp <- ifelse(disC, b, ifelse(disT, a, NA_integer_))

  opp_region <- ifelse(
    joint, "Joint",
    ifelse(opp > c2, "LeadingEdge",
           ifelse(opp >= c1, "Divergence", "Remainder"))
  )

  retained <- joint | opp_region %in% c("LeadingEdge","Divergence")

  data.frame(
    comparison=comparison,
    feature_id=ids,
    selected_k=k,
    control_rank=a,
    treatment_rank=b,
    class=ifelse(joint, "Joint", ifelse(disC, "Control_only", "Treatment_only")),
    opposite_region=opp_region,
    retained=retained,
    stringsAsFactors=FALSE
  )
}

# =============================================================================
# FIGURES
# =============================================================================

theme_pub <- function() {
  theme_classic(base_size=12) +
    theme(
      plot.title=element_text(face="bold", size=13),
      plot.subtitle=element_text(size=10.5, color="grey30"),
      axis.title=element_text(face="bold"),
      axis.text=element_text(color="grey20"),
      panel.border=element_rect(fill=NA, color="grey75", linewidth=.4),
      legend.position="top",
      legend.title=element_blank(),
      plot.margin=margin(8,10,8,10)
    )
}

shade_regions <- function(p, c1, c2, N) {
  p +
    annotate("rect", xmin=1, xmax=c1, ymin=-Inf, ymax=Inf,
             fill=COL$remainder, alpha=.28) +
    annotate("rect", xmin=c1, xmax=c2, ymin=-Inf, ymax=Inf,
             fill=COL$divergence, alpha=.28) +
    annotate("rect", xmin=c2, xmax=N, ymin=-Inf, ymax=Inf,
             fill=COL$leading, alpha=.28) +
    geom_vline(xintercept=c1, color=COL$c1, linetype="dashed", linewidth=.7) +
    geom_vline(xintercept=c2, color=COL$c2, linetype="longdash", linewidth=.7)
}

save_plot <- function(p, file, w=12, h=7) {
  ggsave(file, p, width=w, height=h, dpi=PNG_DPI, bg="white")
  ggsave(sub("\\.png$", ".pdf", file), p, width=w, height=h, bg="white")
}

plot_global_regime <- function(method, arms, knot, file) {
  N <- knot$N
  obs <- bind_rows(lapply(names(arms), \(g)
    transmute(arms[[g]], rank, D, arm=g)
  ))

  med <- obs %>%
    group_by(rank) %>%
    summarise(D=median(D), .groups="drop")

  fit_med <- data.frame(
    rank=seq_len(N),
    fit=apply(knot$fitted, 1, median)
  )

  p <- shade_regions(
    ggplot() +
      geom_line(data=med, aes(rank, D), color=COL$observed, linewidth=1.15) +
      geom_line(data=fit_med, aes(rank, fit), color=COL$fit, linewidth=1.15),
    knot$c1, knot$c2, N
  ) +
    annotate("label", x=knot$c1, y=Inf, label=paste0("c1 = ",knot$c1),
             hjust=1.05, vjust=1.2, size=3, fill="white", label.size=.15) +
    annotate("label", x=knot$c2, y=Inf, label=paste0("c2 = ",knot$c2),
             hjust=-.05, vjust=1.2, size=3, fill="white", label.size=.15) +
    labs(
      title=paste0(method, ": shared experiment-wide divergence geometry"),
      subtitle="Median observed D(r) and shared two-knot fit across all 8 arms",
      x="PC1 rank",
      y="D(r) = F_E(r) - F_P(r)"
    ) +
    theme_pub() +
    theme(legend.position="none")

  save_plot(p, file, 12.5, 6.8)
}

plot_global_cutoff <- function(method, global_scan, frontier, kstar, file) {
  sel <- global_scan %>% filter(k == kstar) %>% slice(1)

  p <- ggplot() +
    geom_path(
      data=frontier,
      aes(weighted_penalty, weighted_benefit),
      color=COL$observed, linewidth=1.3
    ) +
    geom_point(
      data=sel,
      aes(weighted_penalty, weighted_benefit),
      shape=23, size=5, fill=COL$selected, color=COL$selected
    ) +
    annotate(
      "label",
      x=sel$weighted_penalty,
      y=sel$weighted_benefit,
      label=paste0(
        "k* = ", kstar,
        "\nJoint = ", sel$joint_n,
        "\nRemainder crossings = ", sel$remainder_cross_n
      ),
      hjust=-.05, vjust=1.05, size=3.2,
      fill="white", label.size=.15
    ) +
    labs(
      title=paste0(method, ": global weighted Pareto knee"),
      subtitle="Benefit rewards overlap and opposite-arm support; penalty increases with Remainder depth and rank disagreement",
      x="Distance-weighted Remainder penalty",
      y="Weighted agreement benefit"
    ) +
    theme_pub() +
    theme(legend.position="none")

  save_plot(p, file, 10.5, 7)
}

plot_comparison <- function(method, comparison, control, treatment,
                            knot, pair_scan, kstar, file) {
  N <- knot$N
  pdat <- bind_rows(
    transmute(control, rank, D, arm="Control"),
    transmute(treatment, rank, D, arm="Treatment")
  )

  pA <- shade_regions(
    ggplot(pdat, aes(rank, D, linetype=arm)) +
      geom_line(linewidth=1.0),
    knot$c1, knot$c2, N
  ) +
    scale_linetype_manual(values=c(Control="solid", Treatment="longdash")) +
    labs(
      title=paste0("A  ", comparison, " divergence"),
      subtitle="Only the two arm-specific D(r) curves are shown",
      x=NULL, y="D(r)"
    ) +
    theme_pub()

  pair_scaled <- pair_scan %>%
    mutate(
      benefit_scaled=normalize01(weighted_benefit),
      penalty_scaled=normalize01(weighted_penalty)
    ) %>%
    select(k, benefit_scaled, penalty_scaled) %>%
    pivot_longer(-k, names_to="metric", values_to="value")

  pB <- ggplot(pair_scaled, aes(k, value, linetype=metric)) +
    geom_line(linewidth=1.05) +
    geom_vline(xintercept=kstar, color=COL$selected, linetype="dotdash", linewidth=.85) +
    annotate(
      "label", x=kstar, y=1,
      label=paste0("shared k* = ",kstar),
      hjust=-.05, vjust=1.2, size=3,
      fill="white", label.size=.15
    ) +
    scale_linetype_manual(
      values=c(benefit_scaled="solid", penalty_scaled="longdash"),
      labels=c(benefit_scaled="Weighted benefit", penalty_scaled="Weighted penalty")
    ) +
    labs(
      title="B  Pair-specific behavior at the shared method cutoff",
      subtitle="Two normalized curves only: agreement benefit and disagreement penalty",
      x="Candidate top-k",
      y="Scaled score"
    ) +
    theme_pub()

  png(file, width=13, height=9.5, units="in", res=PNG_DPI, bg="white")
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2,1, heights=unit(c(1,1),"null"))))
  print(pA, vp=viewport(layout.pos.row=1))
  print(pB, vp=viewport(layout.pos.row=2))
  dev.off()

  pdf(sub("\\.png$", ".pdf", file), width=13, height=9.5, useDingbats=FALSE)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2,1, heights=unit(c(1,1),"null"))))
  print(pA, vp=viewport(layout.pos.row=1))
  print(pB, vp=viewport(layout.pos.row=2))
  dev.off()
}

plot_method_comparison <- function(summary, file) {
  p <- ggplot(summary, aes(method, selected_k)) +
    geom_col(width=.55) +
    geom_text(aes(label=selected_k), vjust=-.5, fontface="bold", size=4) +
    labs(
      title="Empirical EVS cutoff by normalization method",
      subtitle="One experiment-wide cutoff for CPM-EVS and one for DESeq2-EVS",
      x=NULL, y="Selected k*"
    ) +
    theme_pub() +
    theme(legend.position="none")
  save_plot(p, file, 8.5, 6.2)
}

# =============================================================================
# RUN
# =============================================================================

counts <- read_counts(COUNT_FILE)
groups <- assign_groups(colnames(counts))
norm_obj <- deseq2_normalize(counts, groups)
norm <- norm_obj$counts
pooled_var <- pooled_within_group_var(norm, groups)

all_summary <- list()
all_sites <- list()

for (method in METHODS) {

  method_dir <- file.path(OUT_ROOT, method)
  fig_dir <- file.path(method_dir, "Figures")
  tab_dir <- file.path(method_dir, "Tables")
  dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
  dir.create(tab_dir, recursive=TRUE, showWarnings=FALSE)

  # Eight arm-specific D(r) curves under one method.
  arms <- list()

  for (g in levels(groups)) {
    idx <- which(groups == g)
    arms[[g]] <- build_arm(
      method=method,
      arm=g,
      raw=counts[,idx,drop=FALSE],
      norm=norm[,idx,drop=FALSE],
      pooled_var=pooled_var
    )
  }

  # One shared regime system per normalization method.
  knot <- fit_shared_knots(arms)
  K <- knot$N - knot$c2
  if (K < 1) stop("No Leading-edge candidate domain for ", method)

  # Pair-specific scans are diagnostics; the selected cutoff is global.
  pair_scans <- list()

  for (nm in names(COMPARISONS)) {
    mp <- COMPARISONS[[nm]]
    pair_scans[[nm]] <- scan_pair(
      arms[[mp[["control"]]]],
      arms[[mp[["treatment"]]]],
      knot$c1,
      knot$c2
    ) %>% mutate(comparison=nm)
  }

  global_scan <- aggregate_scans(pair_scans)
  frontier <- pareto_frontier(global_scan)
  kstar <- pareto_knee(frontier)

  # Export global method results.
  write.csv(global_scan,
            file.path(tab_dir, paste0("Table_",method,"_Global_Cutoff_Scan.csv")),
            row.names=FALSE)
  write.csv(frontier,
            file.path(tab_dir, paste0("Table_",method,"_Global_Pareto_Frontier.csv")),
            row.names=FALSE)

  regime <- data.frame(
    method=method,
    N=knot$N,
    c1=knot$c1,
    c2=knot$c2,
    divergence_size=knot$c2-knot$c1+1L,
    leading_edge_candidate_size=K,
    selected_k=kstar,
    k_over_leading_edge=kstar/K,
    knot_SSE=knot$SSE
  )

  write.csv(regime,
            file.path(tab_dir, paste0("Table_",method,"_Key_Results.csv")),
            row.names=FALSE)

  plot_global_regime(
    method, arms, knot,
    file.path(fig_dir, paste0("Figure_",method,"_Shared_Regime.png"))
  )

  plot_global_cutoff(
    method, global_scan, frontier, kstar,
    file.path(fig_dir, paste0("Figure_",method,"_Global_Cutoff.png"))
  )

  # Apply the one method-specific cutoff to every RT/ZT comparison.
  method_sites <- list()

  for (nm in names(COMPARISONS)) {
    mp <- COMPARISONS[[nm]]
    control <- arms[[mp[["control"]]]]
    treatment <- arms[[mp[["treatment"]]]]

    cls <- classify_at_k(
      kstar, nm, control, treatment, knot$c1, knot$c2
    ) %>% mutate(method=method)

    method_sites[[nm]] <- cls

    write.csv(
      cls,
      file.path(tab_dir, paste0("Table_",method,"_",nm,"_Selected_Sites.csv")),
      row.names=FALSE
    )

    write.csv(
      pair_scans[[nm]],
      file.path(tab_dir, paste0("Table_",method,"_",nm,"_Pair_Scan.csv")),
      row.names=FALSE
    )

    plot_comparison(
      method=method,
      comparison=nm,
      control=control,
      treatment=treatment,
      knot=knot,
      pair_scan=pair_scans[[nm]],
      kstar=kstar,
      file=file.path(fig_dir, paste0("Figure_",method,"_",nm,".png"))
    )
  }

  all_summary[[method]] <- regime
  all_sites[[method]] <- bind_rows(method_sites)
}

summary <- bind_rows(all_summary)
sites <- bind_rows(all_sites)

write.csv(summary,
          file.path(OUT_ROOT, "Table_CPM_vs_DESeq2_EVS_Summary.csv"),
          row.names=FALSE)

write.csv(sites,
          file.path(OUT_ROOT, "Table_CPM_vs_DESeq2_EVS_All_Selected_Sites.csv"),
          row.names=FALSE)

# Cross-method overlap by comparison.
overlap <- bind_rows(lapply(names(COMPARISONS), \(nm) {
  a <- sites %>% filter(method=="CPM_EVS", comparison==nm, retained) %>% pull(feature_id)
  b <- sites %>% filter(method=="DESeq2_EVS", comparison==nm, retained) %>% pull(feature_id)
  data.frame(
    comparison=nm,
    CPM_retained=length(unique(a)),
    DESeq2_retained=length(unique(b)),
    overlap=length(intersect(a,b)),
    union=length(union(a,b)),
    jaccard=ifelse(length(union(a,b))>0, length(intersect(a,b))/length(union(a,b)), NA_real_)
  )
}))

write.csv(overlap,
          file.path(OUT_ROOT, "Table_CPM_vs_DESeq2_EVS_Overlap.csv"),
          row.names=FALSE)

plot_method_comparison(
  summary,
  file.path(OUT_ROOT, "Figure_CPM_vs_DESeq2_EVS_Cutoff_Comparison.png")
)

message("============================================================")
message("CPM-EVS vs DESeq2-EVS analysis complete")
message("Output: ", OUT_ROOT)
print(summary)
message("============================================================")
