# =============================================================================
# SEQUENCE - UNIFIED MANUSCRIPT PIPELINE
# Build: SEQUENCE_UNIFIED v1.0.0
#
# SCOPE
# Derives the empirical EVS cutoff from PC1 excess-variance geometry, applies it
# to split each comparison into leading-edge and remainder strata, tests both
# strata with DESeq2 under four decision rules, and writes the manuscript
# figures and tables. One invocation reproduces the complete study.
#
# -----------------------------------------------------------------------------
# FIGURE AND OUTPUT CONVENTIONS
# -----------------------------------------------------------------------------
# Figure canvases are specified in millimetres at final print size and capped at
# journal page geometry: 180 mm double column, 88 mm single column, 240 mm
# maximum height. Base text size and font family are set once and applied to
# every panel, so no figure is reduced during production.
#
# Rendering is routed through cairo devices (cairo_pdf, png type="cairo", LZW
# TIFF) so PDF and PNG renditions of a figure agree, PDF fonts are embedded and
# raster output is antialiased. Dense non-significant volcano layers are
# rasterised when ggrastr is installed, keeping vector files small while axes,
# rules and labelled points stay vector. Assembled panels are aligned with
# patchwork when available and carry A/B/C panel tags.
#
# The method palette is the Okabe-Ito colour-blind-safe set. Volcano captions
# are ASCII and report how many PASs fall above the shared y-axis limit rather
# than clipping silently.
#
# Provenance: a fixed RNG seed, sessionInfo() and package-version manifests are
# written into every output tree, together with a Figure_Manifest.csv recording
# the canvas size and resolution of every figure. Deleting an output root
# recursive deletion requires the pipeline sentinel file in the output root.
# =============================================================================

#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Global reproducibility settings. Applied in BOTH the master process and every
# child run, before any stochastic or order-dependent step.
# -----------------------------------------------------------------------------
SEQUENCE_BUILD_VERSION <- "SEQUENCE_UNIFIED v1.0.0"
SEQUENCE_SEED <- 20260823L
set.seed(SEQUENCE_SEED)

# Millimetre-based figure geometry shared by the empirical-cutoff module and
# the downstream manuscript figures. Sizes are FINAL PRINT sizes.
FIG_SINGLE_COL_MM <- 88     # single-column width
FIG_DOUBLE_COL_MM <- 180    # double-column width
FIG_MAX_HEIGHT_MM <- 240    # maximum printable height on a journal page
fig_mm2in <- function(mm) mm / 25.4

# Set to "Helvetica" or "Arial" if that family is installed and resolvable by
# the graphics device. Empty string means "use the device default", which is
# always safe.
FIG_FONT_FAMILY <- ""

# Output raster formats. TIFF is off by default because it is slow and large;
# enable it when the target journal requires TIFF.
FIG_WRITE_TIFF <- FALSE

# Maximum number of panel columns before an assembled figure wraps onto a new
# row. Five volcano panels in a single 180 mm row are unreadable in print.
PANEL_MAX_COLS <- 3L

PIPELINE_BUILD <- "SEQUENCE_REFINED_EMPIRICAL_EVS_HC10_2026-08-23"

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
#   CPMEVS Lead / Rem   [CPM empirical cutoff child]
#   VSTEVS Lead / Rem   [VST empirical cutoff child]
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
# the selected PAS subset and estimates its own size factors and dispersions.
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
# PROVENANCE
# =============================================================================
# Written into every output tree so a reviewer can reconstruct the exact
# software state that produced the submitted figures and tables.
write_sequence_provenance <- function(dir_path, scope = "run") {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }

  info_path <- file.path(dir_path, paste0("Provenance_SessionInfo_", scope, ".txt"))

  header <- c(
    SEQUENCE_BUILD_VERSION,
    paste0("scope: ", scope),
    paste0("timestamp: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    paste0("RNG seed: ", SEQUENCE_SEED),
    paste0("RNG kind: ", paste(RNGkind(), collapse = " / ")),
    paste0("platform: ", R.version$platform),
    paste0("R version: ", as.character(getRversion())),
    paste0("cairo available: ", isTRUE(capabilities("cairo"))),
    paste0("working directory: ", getwd()),
    "",
    "---- sessionInfo() ----"
  )

  session_lines <- tryCatch(
    utils::capture.output(utils::sessionInfo()),
    error = function(e) paste("sessionInfo() failed:", conditionMessage(e))
  )

  writeLines(c(header, session_lines), con = info_path)

  pkgs <- c(
    "DESeq2", "apeglm", "fdrtool", "ggplot2", "ggrepel", "dplyr", "tidyr",
    "gridExtra", "scales", "S4Vectors", "patchwork", "ggrastr"
  )
  versions <- vapply(
    pkgs,
    function(p) {
      if (requireNamespace(p, quietly = TRUE)) {
        as.character(utils::packageVersion(p))
      } else {
        NA_character_
      }
    },
    character(1)
  )
  utils::write.csv(
    data.frame(
      package = pkgs,
      version = unname(versions),
      installed = !is.na(versions),
      stringsAsFactors = FALSE
    ),
    file.path(dir_path, paste0("Provenance_PackageVersions_", scope, ".csv")),
    row.names = FALSE
  )

  invisible(info_path)
}

# =============================================================================
# DOWNSTREAM COMPARISON SWITCHES
# =============================================================================
# All four prespecified RT/ZT comparisons are enabled for downstream analysis.
# Empirical cutoff derivation uses all eight experimental arms.
RUN_COMPARISON_FLAGS <- c(
  RT0_ZT6  = TRUE,
  RT2_ZT8  = TRUE,
  RT4_ZT10 = TRUE,
  RT8_ZT14 = TRUE
)

if (!any(RUN_COMPARISON_FLAGS)) {
  stop("At least one RUN_COMPARISON_FLAGS entry must be TRUE.")
}

# Active EVS tracks by cutoff method:
# Fixed 5,000: NormEVS and RawEVS. CPM empirical: NormEVS, RawEVS, and CPMEVS.
# VST empirical: NormEVS, RawEVS, and VSTEVS. All DESeq2 tests use raw counts.
RUN_MATCHED_TRANSFORM_TRACKS_ONLY <- TRUE

# Dispersion diagnostics re-estimate DESeq2 dispersions across the specified k grid
# and include the active comparison-specific cutoff as an exactly evaluated point.
ABBREV_EVS <- "excess-variance selection"
ABBREV_RT  <- ""

EXPORT_DISPERSION_TRADEOFF <- TRUE
DISPERSION_SWEEP_ENABLED   <- TRUE
DISPERSION_SWEEP_TRACK     <- "normalized_evs"
DISPERSION_SWEEP_GRID      <- seq(2000L, 9000L, by = 1000L)


# =============================================================================
# INTEGRATED EMPIRICAL CUTOFF + MANUSCRIPT FIGURE ENGINE
# =============================================================================
# Runs once in the master process. The same fitted CPM-EVS/VST-EVS objects are
# used both to render the cutoff manuscript figures/tables and to provide the
# comparison-specific k* values consumed by the downstream SEQUENCE child runs.
# No empirical cutoff is hard-coded.
# =============================================================================

run_sequence_empirical_cutoff_module <- function(count_path, out_root) {

  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(grid)
  })

  options(stringsAsFactors = FALSE)

  # =============================================================================
  # CPM-EVS vs VST-EVS - FINAL MANUSCRIPT VERSION
  # =============================================================================
  #
  # DESIGN
  # ------
  # 1) One experiment-wide PAS universe.
  # 2) Two EVS preprocessing methods:
  #      CPM-EVS    = log1p(CPM) -> arm-specific PCA
  #      VST-EVS     = DESeq2 variance-stabilizing transformation -> arm-specific PCA
  #    VST changes ONLY the EVS/PCA geometry; the pooled NB excess-variance
  #    reference remains based on globally DESeq2-normalized counts for both methods.
  # 3) For each method, all 8 arm-specific D_g(r) curves are fit jointly to
  #    estimate one shared c1/c2 regime system.
  # 4) Each RT/ZT comparison gets its own independent cutoff k*.
  # 5) Cutoff scoring is count-based and identical across preprocessing methods:
  #      - Joint sites receive full unit benefit.
  #      - Disjoint sites are permissible when the opposite-arm rank is in the
  #        Leading Edge or Divergence interval; each receives full unit benefit.
  #      - Only disjoint sites crossing into the opposite-arm Remainder count as
  #        contamination.
  # 6) On the comparison-specific Pareto frontier, G(k) and R(k) are min-max
  #    normalized to [0,1]. The two normalized frontier endpoints define a chord.
  #    k* is the Pareto-optimal candidate with MAXIMUM perpendicular distance
  #    from that endpoint chord (maximum endpoint deviation).
  # 7) Remainder-crossing disjoint sites remain in the top-k union and are flagged
  #    as penalized/low-confidence. The high-confidence subset contains Joint,
  #    opposite-Leading-Edge, and opposite-Divergence sites.
  # 8) Exactly four composite manuscript figures are generated, at final
  #    print size, by the publication figure layer defined below.
  #
  # =============================================================================

  COUNT_FILE <- count_path
  OUT_ROOT <- out_root

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

  METHODS <- c("CPM_EVS", "VST_EVS")

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
    if (!requireNamespace("SummarizedExperiment", quietly=TRUE)) {
      stop("SummarizedExperiment is required for the VST assay.")
    }

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData=round(counts),
      colData=data.frame(group=groups, row.names=colnames(counts)),
      design=~group
    )

    dds <- tryCatch(
      DESeq2::estimateSizeFactors(dds),
      error=\(e) DESeq2::estimateSizeFactors(dds, type="poscounts")
    )

    # VST is used ONLY for the EVS/PCA geometry.  blind=TRUE prevents the
    # experimental design from being used to preserve group differences during
    # the transformation, keeping the PCA comparison unsupervised.
    vst_obj <- DESeq2::varianceStabilizingTransformation(dds, blind=TRUE)
    vst_mat <- SummarizedExperiment::assay(vst_obj)

    if (!identical(dim(vst_mat), dim(counts))) {
      stop("VST matrix dimensions do not match the count matrix.")
    }
    if (!identical(rownames(vst_mat), rownames(counts)) ||
        !identical(colnames(vst_mat), colnames(counts))) {
      stop("VST matrix feature/sample ordering does not match the count matrix.")
    }
    if (any(!is.finite(vst_mat))) stop("VST matrix contains non-finite values.")

    list(
      counts=DESeq2::counts(dds, normalized=TRUE),
      size_factors=DESeq2::sizeFactors(dds),
      vst=vst_mat
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

  build_arm <- function(method, arm, raw, norm, vst, pooled_var) {
    rank_matrix <- switch(
      method,
      CPM_EVS=cpm_log1p(raw),
      VST_EVS=vst,
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

    # Joint sites become active when both arms contain the PAS.
    joint_depth <- pmax(dC,dT)
    joint_n <- activate_from(joint_depth,rep(1,N),K)

    # Control-only interval.
    idxC <- which(dC<dT & dC<=K)
    startC <- dC[idxC]
    endC <- pmin(dT[idxC],K+1L)
    oppC <- rT[idxC]
    disC_n   <- active_interval(startC,endC,rep(1,length(idxC)),K)
    disC_rem <- active_interval(startC,endC,as.numeric(oppC<c1),K)
    disC_div <- active_interval(startC,endC,as.numeric(oppC>=c1 & oppC<=c2),K)
    disC_le  <- active_interval(startC,endC,as.numeric(oppC>c2),K)

    # Treatment-only interval.
    idxT <- which(dT<dC & dT<=K)
    startT <- dT[idxT]
    endT <- pmin(dC[idxT],K+1L)
    oppT <- rC[idxT]
    disT_n   <- active_interval(startT,endT,rep(1,length(idxT)),K)
    disT_rem <- active_interval(startT,endT,as.numeric(oppT<c1),K)
    disT_div <- active_interval(startT,endT,as.numeric(oppT>=c1 & oppT<=c2),K)
    disT_le  <- active_interval(startT,endT,as.numeric(oppT>c2),K)

    disjoint_n <- disC_n+disT_n
    opposite_le_n <- disC_le+disT_le
    opposite_divergence_n <- disC_div+disT_div
    remainder_cross_n <- disC_rem+disT_rem
    good_n <- joint_n+opposite_le_n+opposite_divergence_n
    union_n <- 2*seq_len(K)-joint_n

    if (!all(good_n+remainder_cross_n == union_n)) {
      stop("Internal union-count mismatch in count-based top-k scan.")
    }

    data.frame(
      k=seq_len(K),
      joint_n=joint_n,
      disjoint_n=disjoint_n,
      opposite_le_n=opposite_le_n,
      opposite_divergence_n=opposite_divergence_n,
      remainder_cross_n=remainder_cross_n,
      good_n=good_n,
      union_n=union_n,
      retained_fraction=ifelse(union_n>0,good_n/union_n,NA_real_),
      remainder_cross_fraction=ifelse(union_n>0,remainder_cross_n/union_n,NA_real_),
      stringsAsFactors=FALSE
    )
  }

  # =============================================================================
  # PARETO FRONTIER + MAXIMUM ENDPOINT-CHORD DEVIATION
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

  select_max_endpoint_deviation <- function(scan_df,
                                            good_col="good_n",
                                            cost_col="remainder_cross_n") {
    marked <- mark_pareto_frontier(scan_df, good_col, cost_col)
    frontier <- marked$frontier %>% arrange(cost,good,k)
    if (nrow(frontier) < 2L) {
      stop("At least two Pareto-optimal points are required to define the endpoint chord.")
    }

    # Normalize ONLY on the Pareto frontier.
    gr <- range(frontier$good, na.rm=TRUE)
    cr <- range(frontier$cost, na.rm=TRUE)
    norm_good <- function(x) {
      if (diff(gr)==0) rep(0,length(x)) else (x-gr[1])/diff(gr)
    }
    norm_cost <- function(x) {
      if (diff(cr)==0) rep(0,length(x)) else (x-cr[1])/diff(cr)
    }

    frontier$good_norm <- norm_good(frontier$good)
    frontier$remainder_norm <- norm_cost(frontier$cost)

    # Endpoint chord in normalized Pareto space.
    x1 <- frontier$remainder_norm[1]
    y1 <- frontier$good_norm[1]
    x2 <- frontier$remainder_norm[nrow(frontier)]
    y2 <- frontier$good_norm[nrow(frontier)]

    denom <- sqrt((y2-y1)^2 + (x2-x1)^2)
    if (!is.finite(denom) || denom <= .Machine$double.eps) {
      stop("Normalized Pareto endpoints are coincident; endpoint chord is undefined.")
    }

    frontier$endpoint_deviation <- abs(
      (y2-y1)*frontier$remainder_norm -
        (x2-x1)*frontier$good_norm +
        x2*y1 - y2*x1
    ) / denom

    best <- max(frontier$endpoint_deviation, na.rm=TRUE)
    chosen <- frontier %>%
      filter(abs(endpoint_deviation-best) < 1e-12) %>%
      arrange(desc(good), cost, desc(k)) %>%
      slice(1L)

    out <- marked$scan
    out$good_norm <- norm_good(out[[good_col]])
    out$remainder_norm <- norm_cost(out[[cost_col]])
    out$endpoint_deviation <- abs(
      (y2-y1)*out$remainder_norm -
        (x2-x1)*out$good_norm +
        x2*y1 - y2*x1
    ) / denom
    out$is_selected_endpoint_max <- out$k == chosen$k[1]

    chord <- data.frame(
      remainder_norm=c(x1,x2),
      good_norm=c(y1,y2),
      stringsAsFactors=FALSE
    )

    list(
      scan=out,
      frontier=frontier,
      chord=chord,
      selected_k=as.integer(chosen$k[1]),
      selected_good=chosen$good[1],
      selected_cost=chosen$cost[1],
      selected_good_norm=chosen$good_norm[1],
      selected_remainder_norm=chosen$remainder_norm[1],
      selected_endpoint_deviation=chosen$endpoint_deviation[1]
    )
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

    # Binary bookkeeping consistent with the validated count-based objective.
    penalty <- ifelse(joint,0,as.numeric(opp<c1))
    support <- ifelse(joint,1,as.numeric(opp>=c1))

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

  EVS_METHOD_COL <- c("CPM-EVS" = EVS_COL$cpm, "VST-EVS" = EVS_COL$deseq)

  evs_method_label <- function(m) ifelse(m == "CPM_EVS", "CPM-EVS", "VST-EVS")

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
  # 3. FIGURE 1 PANELS - regime definition
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
      scale_linetype_manual(values = c("CPM-EVS" = "solid", "VST-EVS" = "22")) +
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
  # 4. FIGURES 2-3 PANELS - candidate path, Pareto frontier, endpoint chord
  # =============================================================================

  panel_pareto <- function(method, comparison, scan, frontier, kstar, tag,
                           method_col, show_legend = FALSE) {
    fr <- frontier[order(frontier$cost,frontier$good), , drop=FALSE]
    sel <- fr[fr$k==kstar,,drop=FALSE]
    if (nrow(sel)!=1L) stop("k* not found on the Pareto frontier for ", comparison)
    if (nrow(fr)<2L) stop("Need at least two Pareto points for endpoint chord in ", comparison)

    # Full candidate trajectory gives context; Pareto frontier is emphasized.
    sc <- scan[order(scan$k), , drop=FALSE]
    sc_plot <- sc[evs_thin(nrow(sc),2500L),,drop=FALSE]
    fr_plot <- fr[evs_thin(nrow(fr),1500L),,drop=FALSE]

    xr <- range(c(sc$remainder_cross_n,fr$cost),na.rm=TRUE)
    yr <- range(c(sc$good_n,fr$good),na.rm=TRUE)
    xd <- if (diff(xr)>0) diff(xr) else 1
    yd <- if (diff(yr)>0) diff(yr) else 1

    # Chord endpoints are the first and last points of the normalized frontier,
    # mapped back to the raw R(k), G(k) coordinates for display.
    chord_raw <- data.frame(
      cost=c(fr$cost[1],fr$cost[nrow(fr)]),
      good=c(fr$good[1],fr$good[nrow(fr)])
    )

    # The candidate trajectory, the Pareto frontier and the endpoint chord share
    # one key. The selected point carries its own k* label and needs no entry.
    lab_scan   <- "All candidate k"
    lab_front  <- "Pareto frontier"
    lab_chord  <- "Endpoint chord"
    key_levels <- c(lab_scan, lab_front, lab_chord)
    sc_plot$key    <- lab_scan
    fr_plot$key    <- lab_front
    chord_raw$key  <- lab_chord

    ins <- ggplot(fr_plot,aes(k,endpoint_deviation)) +
      geom_line(colour="grey25",linewidth=0.28) +
      geom_vline(xintercept=sel$k,colour=EVS_COL$knee,
                 linetype="22",linewidth=0.32) +
      geom_point(data=sel,aes(k,endpoint_deviation),
                 inherit.aes=FALSE,shape=21,size=1.0,stroke=0.25,
                 fill=EVS_COL$knee,colour="white") +
      labs(x=expression(italic(k)),y="Endpoint deviation") +
      evs_theme_inset()

    ggplot() +
      geom_path(data=sc_plot,aes(remainder_cross_n,good_n,colour=key,linetype=key),
                linewidth=0.32,alpha=0.9) +
      geom_line(data=fr_plot,aes(cost,good,colour=key,linetype=key),
                linewidth=0.65) +
      geom_line(data=chord_raw,aes(cost,good,colour=key,linetype=key),
                linewidth=0.42) +
      annotation_custom(
        ggplotGrob(ins),
        xmin=xr[1]+0.50*(xr[2]-xr[1]+1e-9),
        xmax=xr[2]+0.05*xd,
        ymin=yr[1]+0.04*(yr[2]-yr[1]+1e-9),
        ymax=yr[1]+0.46*(yr[2]-yr[1]+1e-9)
      ) +
      geom_point(data=sel,aes(cost,good),shape=23,size=1.6,stroke=0.3,
                 fill=EVS_COL$knee,colour="white") +
      annotate(
        "text",
        x=sel$cost-0.015*xd,
        y=sel$good+0.035*yd,
        label=sprintf("italic(k)^\"*\" * \" = \" * \"%s\"",evs_num(sel$k)),
        parse=TRUE,
        hjust=1,vjust=0,size=pt2mm(6.2),
        colour=EVS_COL$knee,fontface="bold"
      ) +
      scale_colour_manual(
        breaks=key_levels,
        values=setNames(c("grey78",method_col,"grey45"),key_levels)
      ) +
      scale_linetype_manual(
        breaks=key_levels,
        values=setNames(c("solid","solid","22"),key_levels)
      ) +
      scale_x_continuous(labels=evs_comma) +
      scale_y_continuous(labels=evs_comma) +
      labs(
        tag=tag,
        title=evs_comparison_label(comparison),
        x=expression(paste("Opposite-arm Remainder crossings, ",italic(R)(k))),
        y=expression(paste("Joint + permissible Disjoint, ",italic(G)(k)))
      ) +
      evs_theme() +
      (if (isTRUE(show_legend)) {
        list(
          guides(colour=guide_legend(order=1), linetype=guide_legend(order=1)),
          evs_legend_inside(0.02, 0.985, c(0, 1))
        )
      } else NULL)
  }

  # =============================================================================
  # 5. FIGURE 4 PANELS - cross-method summary
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
        vst_only = VST_high_confidence - overlap_high_confidence,
        comparison = factor(evs_comparison_label(comparison),
                            levels = evs_comparison_label(names(comparisons)))
      )

    tot <- d %>%
      transmute(comparison,
                total   = cpm_only + shared + vst_only,
                jaccard = jaccard_high_confidence)

    long <- d %>%
      select(comparison, cpm_only, shared, vst_only) %>%
      tidyr::pivot_longer(-comparison, names_to = "set", values_to = "n") %>%
      mutate(set = factor(recode(set,
                                 cpm_only   = "CPM-EVS only",
                                 shared     = "Shared",
                                 vst_only = "VST-EVS only"),
                          levels = c("CPM-EVS only", "Shared", "VST-EVS only")))

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
                                   "VST-EVS only" = EVS_COL$deseq)) +
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
      scale_x_discrete(labels = c("CPM-EVS" = "CPM", "VST-EVS" = "VST")) +
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
      paste0("median-of-ratios normalisation; EVS, ", ABBREV_EVS, "; F_E(r),"),
      "cumulative excess-variance mass; F_P(r), cumulative PC1 variance-contribution",
      "mass; G(k), Joint plus permissible Disjoint benefit; HC, high confidence; J, Jaccard",
      "index; k*, selected per-comparison cutoff; LE, leading edge; N, size of the",
      "experiment-wide PAS universe; PAS, polyadenylation site; PC1, first principal",
      paste0("component; R(k), opposite-arm Remainder-crossing contamination count; ",
             if (nzchar(ABBREV_RT)) paste0("RT, ", ABBREV_RT, "; ") else "",
             "ZT, zeitgeber time."),
      "",
      "---",
      "",
      "**Figure 1. Definition of the Remainder, Divergence and Leading-Edge rank",
      "regimes under two normalisation strategies.**",
      sprintf(paste0(
        "(A, B) Cumulative PC1-excess-variance divergence, D(r) = F_E(r) - F_P(r), ",
        "against PAS rank by ascending absolute PC1 loading, for (A) CPM-EVS and ",
        "(B) VST-EVS. Grey lines show all eight experimental arms, the coloured ",
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
        "(A-D) Pareto frontier of count-based retained-site benefit G(k) against ",
        "opposite-arm Remainder-crossing contamination R(k) over all candidate cutoffs ",
        "k, shown separately for each RT/ZT comparison; each comparison is optimised ",
        "independently over 1 <= k <= N - c2. The filled diamond marks k*, the ",
        "Pareto-optimal candidate with the maximum perpendicular deviation from ",
        "the chord joining the normalized Pareto endpoints after frontier-specific ",
        "min-max normalization of benefit and contamination. The dashed line is the ",
        "endpoint chord; insets plot endpoint deviation across the Pareto frontier. ",
        "Non-dominated candidates make ",
        "up %.0f%% of the k values scanned under CPM-EVS. Selected values: ",
        "k* = %s for %s respectively; full diagnostics are in Table 1."),
        100 * fr_frac[["CPM_EVS"]],
        paste(evs_num(kstar_of("CPM_EVS")), collapse = ", "),
        paste(evs_comparison_label(names(comparisons)), collapse = ", ")),
      "",
      "---",
      "",
      "**Figure 3. Per-comparison cutoff selection under VST-EVS.**",
      sprintf(paste0(
        "Panels, axes and annotations are as in Figure 2, computed from ",
        "DESeq2 variance-stabilized counts (varianceStabilizingTransformation, blind=TRUE). Non-dominated ",
        "candidates make up %.0f%% of the k values scanned. Selected values: ",
        "k* = %s."),
        100 * fr_frac[["VST_EVS"]],
        paste(evs_num(kstar_of("VST_EVS")), collapse = ", ")),
      "",
      "---",
      "",
      "**Figure 4. CPM-EVS and VST-EVS compared across all four RT/ZT",
      "comparisons.**",
      sprintf(paste0(
        "(A) Selected cutoff k* for each comparison and method (%s). (B) ",
        "High-confidence PAS membership, partitioned into sites recovered only by ",
        "CPM-EVS, only by VST-EVS, or by both; J is the Jaccard index of the two ",
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
      panel_regime("VST_EVS", method_objects[["VST_EVS"]]$arms,
                   method_objects[["VST_EVS"]]$knot, "B"),
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
                     ob$pairs[[nm]]$kstar, tags[i], col,
                     show_legend = (i == 1L))
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
    evs_save(f4, file.path(fig_dir, "Figure_4_CPM_vs_VST_EVS_Summary"),
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
  vst  <- norm_obj$vst

  message("VST matrix computed globally for the EVS/PCA branch (blind=TRUE).")
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
        vst=vst[,idx,drop=FALSE],
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

      opt <- select_max_endpoint_deviation(
        scan, good_col="good_n", cost_col="remainder_cross_n"
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
        good_n=selected$good_n,
        remainder_cross_n=selected$remainder_cross_n,
        good_norm=selected$good_norm,
        remainder_norm=selected$remainder_norm,
        endpoint_deviation=selected$endpoint_deviation,
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
        " | G=",selected$good_n,
        " | R=",selected$remainder_cross_n,
        " | endpoint deviation=",round(selected$endpoint_deviation,4),
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

    vst_sites <- sites_df %>%
      filter(
        method=="VST_EVS",
        comparison==nm,
        high_confidence
      ) %>%
      pull(feature_id) %>%
      unique()

    u <- union(cpm,vst_sites)

    data.frame(
      comparison=nm,
      CPM_high_confidence=length(cpm),
      VST_high_confidence=length(vst_sites),
      overlap_high_confidence=length(intersect(cpm,vst_sites)),
      union_high_confidence=length(u),
      jaccard_high_confidence=ifelse(
        length(u)>0,
        length(intersect(cpm,vst_sites))/length(u),
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
    file.path(TAB_DIR,"Table_2_CPM_vs_VST_HighConfidence_Overlap.csv"),
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

  # Clear prior figure files so manifests and archives contain only the current run.
  old_figure_files <- list.files(
    FIG_DIR,
    pattern="\\.(png|pdf|tif|tiff)$",
    full.names=TRUE,
    ignore.case=TRUE
  )

  if (length(old_figure_files)) unlink(old_figure_files)

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
  # =============================================================================

  methods_lines <- c(
    "# Methods",
    "",
    "CPM-EVS and VST-EVS were evaluated on the common experiment-wide PAS count matrix. CPM-EVS used log1p-transformed counts per million. VST-EVS used the DESeq2 variance-stabilizing transformation with blind=TRUE after median-of-ratios size-factor estimation. For each transformed matrix, PCA was performed separately within each experimental arm and PASs were ranked by absolute PC1 loading.",
    "",
    "For PAS i, PC1 variance contribution was defined as P_i=lambda_1*v_i1^2. Pooled within-arm variance was calculated from experiment-wide DESeq2 median-of-ratios normalized counts, and arm-specific excess variance was E_ig=max(V_pool,i-mu_ig,0). PC1 contribution and excess variance were converted to rank-wise probability masses and cumulative distributions, and cumulative divergence was D_g(r)=F_E,g(r)-F_P,g(r).",
    "",
    "For each EVS preprocessing method, the eight arm-specific divergence curves were jointly fit with a continuous two-knot piecewise-linear model. Shared knots c1 and c2 minimized the summed squared residual error across all arms. Ranks below c1 were classified as Remainder, ranks from c1 through c2 as the Divergence interval, and ranks above c2 as Leading Edge.",
    "",
    "Each RT/ZT comparison was optimized independently over candidate top-k values from 1 through N-c2. At each k, benefit G(k) counted Joint top-k PASs plus Disjoint PASs whose opposite-arm rank remained in the Leading Edge or Divergence interval. Contamination R(k) counted Disjoint PASs whose opposite-arm rank fell below c1 into the Remainder.",
    "",
    "Nondominated [R(k),G(k)] candidates defined the Pareto frontier. Benefit and contamination were min-max normalized on the frontier. The empirical cutoff k* was the Pareto candidate with the greatest perpendicular distance from the chord joining the first and last normalized frontier points. Ties were resolved by greater G, then lower R, then larger k. The selected k* was calculated independently for each comparison and separately for CPM-EVS and VST-EVS."
  )

  writeLines(
    methods_lines,
    file.path(OUT_ROOT,"Methods_Manuscript.md")
  )

  # =============================================================================

  cpm_rows <- summary_df[summary_df$method == "CPM_EVS", , drop=FALSE]
  vst_rows <- summary_df[summary_df$method == "VST_EVS", , drop=FALSE]
  cpm_k <- setNames(as.integer(cpm_rows$selected_k), cpm_rows$comparison)
  vst_k <- setNames(as.integer(vst_rows$selected_k), vst_rows$comparison)

  audit <- summary_df
  names(audit)[names(audit) == "method"] <- "Method"
  names(audit)[names(audit) == "comparison"] <- "Comparison"
  names(audit)[names(audit) == "selected_k"] <- "k_star"

  list(
    cpm = cpm_k,
    vst = vst_k,
    audit = audit,
    method_objects = method_objects,
    summary_df = summary_df,
    overlap_df = overlap_df,
    sites_df = sites_df,
    scans_df = scans_df,
    frontiers_df = frontiers_df,
    output_root = OUT_ROOT,
    figure_dir = FIG_DIR,
    table_dir = TAB_DIR
  )
}

# =============================================================================
# CROSS-CUTOFF SENSITIVITY MODULE
# =============================================================================
# The three cutoff methods differ only in the number of PASs admitted to the
# leading edge: the PC1 ranking, the DESeq2 model and the four decision rules
# are identical across runs, so the comparison is a one-variable sensitivity
# analysis in k.
#
# Because the rankings are fixed, the three top-k sets are nested and overlap
# statistics between them are determined by the ratios of the k values alone.
# What is informative is what happens to the calls, reported in two strata:
#
#   Common stratum   PASs admitted under the smallest cutoff and therefore
#                    tested in all three runs. A call that changes here cannot
#                    be explained by a feature entering or leaving the set; it
#                    reflects the dispersion trend, independent filtering, the
#                    BH denominator and the empirical null being re-fit on a
#                    different surrounding set. Its flip rate is the headline
#                    sensitivity number.
#
#   Margin stratum   PASs admitted only at larger k, which can only gain calls
#                    as k grows. Differences there are expected rather than
#                    evidence of instability.
#
# This module runs in the master process after all child runs complete, since
# no child can see another child's output. It is self-contained for the same
# reason as the empirical cutoff module: the master reaches it before the
# downstream figure helpers are defined.
# =============================================================================

run_cutoff_sensitivity_module <- function(multi_root, methods, cutoff_manifest) {

  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(grid)
  })

  MM_PER_IN <- 25.4
  sens_pretty <- function(x) sub("_", " vs ", x, fixed = TRUE)

  sens_theme <- function(base_pt = 7) {
    theme_classic(base_size = base_pt) +
      theme(
        plot.title   = element_text(size = base_pt + 0.5, face = "bold", hjust = 0,
                                    margin = margin(b = 2.5)),
        plot.caption = element_text(size = base_pt - 1.5, colour = "grey35",
                                    hjust = 0, margin = margin(t = 3)),
        plot.tag     = element_text(size = base_pt + 2, face = "bold"),
        plot.tag.position = "topleft",
        axis.title   = element_text(size = base_pt),
        axis.text    = element_text(size = base_pt - 0.5, colour = "grey15"),
        axis.line    = element_line(linewidth = 0.25, colour = "grey20"),
        axis.ticks   = element_line(linewidth = 0.25, colour = "grey20"),
        strip.background = element_blank(),
        strip.text   = element_text(size = base_pt - 0.5, face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text  = element_text(size = base_pt - 1),
        legend.key.size = unit(7, "pt"),
        panel.background = element_blank(),
        plot.margin  = margin(3, 4, 3, 3)
      )
  }

  sens_save <- function(g, path_base, width_mm = 180, height_mm = 200) {
    w <- width_mm / MM_PER_IN
    h <- height_mm / MM_PER_IN
    cairo_ok <- isTRUE(capabilities("cairo"))
    if (cairo_ok) grDevices::cairo_pdf(paste0(path_base, ".pdf"), width = w, height = h)
    else          grDevices::pdf(paste0(path_base, ".pdf"), width = w, height = h)
    grid::grid.draw(g); grDevices::dev.off()
    if (cairo_ok) {
      grDevices::png(paste0(path_base, ".png"), width = w, height = h, units = "in",
                     res = 600, bg = "white", type = "cairo")
    } else {
      grDevices::png(paste0(path_base, ".png"), width = w, height = h, units = "in",
                     res = 600, bg = "white")
    }
    grid::grid.draw(g); grDevices::dev.off()
  }


  sensitivity_method_levels <- function() {
    c("Fixed_5000", "CPM_Empirical", "VST_Empirical")
  }

  sensitivity_method_labels <- function() {
    c(Fixed_5000 = "Fixed 5,000",
      CPM_Empirical = "CPM k*",
      VST_Empirical = "VST k*")
  }

  # Per-PAS leading-edge membership, written by each child run. Only leading-edge
  # identifiers are stored; the remainder is the complement within a comparison.
  read_membership_tables <- function(multi_root, methods) {
    rows <- list()
    for (m in methods) {
      f <- file.path(multi_root, m, "Summary_Tables",
                     "Table_EVS_Leading_Edge_Membership.csv")
      if (!file.exists(f)) next
      d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
      if (!nrow(d)) next
      d$Cutoff_Method <- m
      rows[[length(rows) + 1L]] <- d
    }
    if (!length(rows)) return(NULL)
    dplyr::bind_rows(rows)
  }

  read_significance_tables <- function(multi_root, methods) {
    rows <- list()
    for (m in methods) {
      root <- file.path(multi_root, m)
      if (!dir.exists(root)) next
      files <- list.files(root, pattern = "^Table_Significant_Sites\\.csv$",
                          recursive = TRUE, full.names = TRUE)
      for (f in files) {
        d <- tryCatch(
          utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
          error = function(e) NULL
        )
        if (is.null(d) || !nrow(d)) next
        d$Cutoff_Method <- m
        rows[[length(rows) + 1L]] <- d
      }
    }
    if (!length(rows)) return(NULL)

    out <- dplyr::bind_rows(rows)
    needed <- c("Comparison", "Analysis", "PAS", "Std", "Strong", "Weak", "HBFSS_sig")
    if (!all(needed %in% names(out))) {
      warning("Significance tables are missing expected columns: ",
              paste(setdiff(needed, names(out)), collapse = ", "))
      return(NULL)
    }
    for (v in c("Std", "Strong", "Weak", "HBFSS_sig")) {
      out[[v]] <- as.logical(out[[v]])
      out[[v]][is.na(out[[v]])] <- FALSE
    }
    out
  }

  # Long form: one row per PAS x method x decision rule that was called.
  sensitivity_long_calls <- function(sig_df) {
    rules <- c(Standard = "Std", Strong = "Strong", Weak = "Weak", HBFSS = "HBFSS_sig")
    dplyr::bind_rows(lapply(names(rules), function(rule) {
      col <- rules[[rule]]
      d <- sig_df[sig_df[[col]], c("Comparison", "Analysis", "PAS", "Cutoff_Method"),
                  drop = FALSE]
      if (!nrow(d)) return(NULL)
      d$Rule <- rule
      d
    }))
  }

  # Agreement across cutoffs, over the union of PASs called by any of them.
  sensitivity_agreement <- function(calls_df, n_methods) {
    if (is.null(calls_df) || !nrow(calls_df)) return(NULL)
    calls_df %>%
      dplyr::distinct(Comparison, Analysis, Rule, PAS, Cutoff_Method) %>%
      dplyr::group_by(Comparison, Analysis, Rule, PAS) %>%
      dplyr::summarise(n_methods_calling = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(Comparison, Analysis, Rule) %>%
      dplyr::summarise(
        union_calls = dplyr::n(),
        called_by_all = sum(n_methods_calling == n_methods),
        called_by_two = sum(n_methods_calling == 2L),
        called_by_one = sum(n_methods_calling == 1L),
        pct_called_by_all = 100 * sum(n_methods_calling == n_methods) / dplyr::n(),
        .groups = "drop"
      )
  }

  # Flip rate within the common stratum. A PAS counts as flipped when it is
  # called by at least one cutoff but not by all of them, restricted to PASs that
  # every cutoff admitted and therefore tested.
  sensitivity_flip_rate <- function(calls_df, member_df, methods) {
    if (is.null(calls_df) || !nrow(calls_df)) return(NULL)

    if (is.null(member_df) || !nrow(member_df)) {
      warning("Leading-edge membership tables were not found; the common-stratum ",
              "flip rate cannot be computed and is omitted.")
      return(NULL)
    }

    # Nested rankings mean the smallest cutoff's leading edge is contained in the
    # others, so its membership defines the commonly tested stratum.
    common <- member_df %>%
      dplyr::distinct(Comparison, Track, PAS, Cutoff_Method) %>%
      dplyr::group_by(Comparison, Track, PAS) %>%
      dplyr::summarise(n_admitting = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(n_admitting == length(methods)) %>%
      dplyr::select(Comparison, Track, PAS)

    if (!nrow(common)) return(NULL)

    # Analysis labels carry the track, e.g. "NormEVS Lead"; remainder views are
    # outside the leading-edge stratum and are excluded here.
    calls <- calls_df[grepl("Lead", calls_df$Analysis, fixed = TRUE), , drop = FALSE]
    if (!nrow(calls)) return(NULL)
    calls$Track <- ifelse(grepl("NormEVS", calls$Analysis, fixed = TRUE),
                          "NormEVS", "RawEVS")

    calls <- dplyr::inner_join(calls, common, by = c("Comparison", "Track", "PAS"))
    if (!nrow(calls)) return(NULL)

    calls %>%
      dplyr::distinct(Comparison, Analysis, Rule, PAS, Cutoff_Method) %>%
      dplyr::group_by(Comparison, Analysis, Rule, PAS) %>%
      dplyr::summarise(n_methods_calling = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(Comparison, Analysis, Rule) %>%
      dplyr::summarise(
        n_called_in_common_stratum = dplyr::n(),
        n_flipped = sum(n_methods_calling < length(methods)),
        pct_flipped = 100 * sum(n_methods_calling < length(methods)) / dplyr::n(),
        .groups = "drop"
      )
  }

  sensitivity_counts <- function(calls_df, cutoff_manifest) {
    if (is.null(calls_df) || !nrow(calls_df)) return(NULL)
    k_long <- tidyr::pivot_longer(
      cutoff_manifest,
      cols = dplyr::any_of(c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14")),
      names_to = "Comparison", values_to = "k"
    )
    k_long$k <- as.integer(k_long$k)

    calls_df %>%
      dplyr::distinct(Comparison, Analysis, Rule, PAS, Cutoff_Method) %>%
      dplyr::group_by(Comparison, Analysis, Rule, Cutoff_Method) %>%
      dplyr::summarise(n_significant = dplyr::n(), .groups = "drop") %>%
      dplyr::left_join(
        k_long[, c("Cutoff_Method", "Comparison", "k")],
        by = c("Cutoff_Method", "Comparison")
      )
  }

  # -----------------------------------------------------------------------------
  # Panels
  # -----------------------------------------------------------------------------

  SENSITIVITY_AGREE_COLORS <- c(
    "Called by all cutoffs" = "#009E73",
    "Called by two"         = "#9AA3AD",
    "Called by one only"    = "#D55E00"
  )

  plot_sensitivity_counts <- function(cnt_df, stratum_regex, panel_title) {
    if (is.null(cnt_df) || !nrow(cnt_df)) return(NULL)
    d <- cnt_df[grepl(stratum_regex, cnt_df$Analysis), , drop = FALSE]
    d <- d[!is.na(d$k), , drop = FALSE]
    if (!nrow(d)) return(NULL)

    d <- d %>%
      dplyr::group_by(Comparison, Analysis, Cutoff_Method, k) %>%
      dplyr::summarise(n_significant = sum(n_significant), .groups = "drop")
    d$Track <- ifelse(grepl("NormEVS", d$Analysis, fixed = TRUE), "NormEVS", "RawEVS")
    d$Comparison <- sens_pretty(d$Comparison)

    ggplot(d, aes(x = k, y = n_significant, colour = Comparison, linetype = Track)) +
      geom_line(linewidth = 0.45) +
      geom_point(size = 0.9) +
      scale_x_continuous(labels = function(v) format(v, big.mark = ",", trim = TRUE)) +
      labs(
        title = panel_title,
        x = "Leading-edge size, k (PAS per condition)",
        y = "Significant PAS (any decision rule)"
      ) +
      sens_theme() +
      theme(legend.position = "bottom", legend.title = element_blank())
  }

  plot_sensitivity_agreement <- function(agree_df) {
    if (is.null(agree_df) || !nrow(agree_df)) return(NULL)
    d <- agree_df[grepl("Lead", agree_df$Analysis, fixed = TRUE), , drop = FALSE]
    if (!nrow(d)) return(NULL)

    long <- d %>%
      dplyr::select(Comparison, Analysis, Rule, called_by_all, called_by_two, called_by_one) %>%
      tidyr::pivot_longer(c("called_by_all", "called_by_two", "called_by_one"),
                          names_to = "agreement", values_to = "n") %>%
      dplyr::mutate(
        agreement = factor(
          dplyr::recode(agreement,
                        called_by_all = "Called by all cutoffs",
                        called_by_two = "Called by two",
                        called_by_one = "Called by one only"),
          levels = names(SENSITIVITY_AGREE_COLORS)
        ),
        Comparison = sens_pretty(Comparison)
      ) %>%
      dplyr::group_by(Comparison, Rule, agreement) %>%
      dplyr::summarise(n = sum(n), .groups = "drop") %>%
      dplyr::group_by(Comparison, Rule) %>%
      dplyr::mutate(pct = 100 * n / sum(n)) %>%
      dplyr::ungroup()

    ggplot(long, aes(x = Rule, y = pct, fill = agreement)) +
      geom_col(width = 0.68, colour = "white", linewidth = 0.2) +
      facet_wrap(~ Comparison, nrow = 1) +
      scale_fill_manual(values = SENSITIVITY_AGREE_COLORS) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
      labs(
        title = "Agreement of calls across the three cutoffs",
        x = NULL, y = "% of union of calls"
      ) +
      sens_theme() +
      theme(
        legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }

  plot_sensitivity_flips <- function(flip_df) {
    if (is.null(flip_df) || !nrow(flip_df)) return(NULL)
    d <- flip_df %>%
      dplyr::group_by(Comparison, Rule) %>%
      dplyr::summarise(
        n_flipped = sum(n_flipped),
        n_total = sum(n_called_in_common_stratum),
        .groups = "drop"
      ) %>%
      dplyr::mutate(pct = 100 * n_flipped / n_total,
                    Comparison = sens_pretty(Comparison))

    ggplot(d, aes(x = Comparison, y = pct, fill = Rule)) +
      geom_col(position = position_dodge(width = 0.78), width = 0.70) +
      geom_text(aes(label = sprintf("%.0f", pct)),
                position = position_dodge(width = 0.78),
                vjust = -0.4, size = 1.9) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
      labs(
        title = "Calls that change within the commonly tested stratum",
        x = NULL, y = "% of calls that flip",
        caption = paste0(
          "Restricted to PASs admitted by every cutoff, so set membership is held ",
          "constant and a flip reflects re-fitting rather than a feature entering ",
          "or leaving the analysis."
        )
      ) +
      sens_theme() +
      theme(legend.position = "bottom", legend.title = element_blank())
  }

  # -----------------------------------------------------------------------------
  # Driver
  # -----------------------------------------------------------------------------

  run_cutoff_sensitivity_analysis <- function(multi_root, methods, cutoff_manifest) {
    out_dir <- file.path(multi_root, "Cutoff_Sensitivity")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    sig_df <- read_significance_tables(multi_root, methods)
    if (is.null(sig_df)) {
      warning("Cutoff sensitivity analysis skipped: no significance tables found.")
      return(invisible(FALSE))
    }

    calls_df <- sensitivity_long_calls(sig_df)
    if (is.null(calls_df) || !nrow(calls_df)) {
      warning("Cutoff sensitivity analysis skipped: no significant calls found.")
      return(invisible(FALSE))
    }

    member_df <- read_membership_tables(multi_root, methods)

    cnt   <- sensitivity_counts(calls_df, cutoff_manifest)
    agree <- sensitivity_agreement(calls_df, length(methods))
    flips <- sensitivity_flip_rate(calls_df, member_df, methods)

    if (!is.null(cnt))   utils::write.csv(cnt,   file.path(out_dir, "Table_Sensitivity_Discovery_Counts.csv"), row.names = FALSE)
    if (!is.null(agree)) utils::write.csv(agree, file.path(out_dir, "Table_Sensitivity_Call_Agreement.csv"), row.names = FALSE)
    if (!is.null(flips)) utils::write.csv(flips, file.path(out_dir, "Table_Sensitivity_Flip_Rate.csv"), row.names = FALSE)

    panels <- list(
      plot_sensitivity_counts(cnt, "Lead", "Discoveries in the leading edge"),
      plot_sensitivity_counts(cnt, "Rem",  "Discoveries left in the remainder"),
      plot_sensitivity_agreement(agree),
      plot_sensitivity_flips(flips)
    )
    panels <- Filter(Negate(is.null), panels)
    if (!length(panels)) {
      warning("Cutoff sensitivity figure skipped: no panel could be built.")
      return(invisible(FALSE))
    }

    tagged <- lapply(seq_along(panels), function(i) {
      panels[[i]] +
        labs(tag = LETTERS[i]) +
        theme(plot.tag = element_text(face = "bold", size = 7 + 1,
                                      family = NULL),
              plot.tag.position = "topleft")
    })

    body <- if (requireNamespace("patchwork", quietly = TRUE)) {
      patchwork::patchworkGrob(
        patchwork::wrap_plots(tagged, ncol = if (length(tagged) > 1L) 2L else 1L)
      )
    } else {
      do.call(gridExtra::arrangeGrob,
              c(tagged, list(ncol = if (length(tagged) > 1L) 2L else 1L)))
    }

    g <- gridExtra::arrangeGrob(
      body, ncol = 1,
      top = grid::textGrob(
        "Sensitivity of SEQUENCE results to the leading-edge cutoff",
        gp = grid::gpar(
          fontface = "bold",
          fontsize = 7 + 2,
          fontfamily = ""
        )
      )
    )

    sens_save(g, file.path(out_dir, "Figure_Cutoff_Sensitivity"),
              width_mm = FIG_DOUBLE_COL_MM, height_mm = 200)

    message("Cutoff sensitivity analysis written to: ", out_dir)
    invisible(TRUE)
  }
}


# =============================================================================
# MULTI-CUTOFF ORCHESTRATION
# Runs this script once for each cutoff method in isolated output folders, then creates one
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
    unset = file.path(repo_guess, "exports", "sequence_final_three_cutoffs")
  )
  # Output roots are deleted only when marked by this pipeline's sentinel file.
  multi_root_sentinel <- ".sequence_output_root"
  if (dir.exists(multi_root)) {
    if (!file.exists(file.path(multi_root, multi_root_sentinel))) {
      # Preserve an unmarked directory and write to a timestamped sibling root.
      original_multi_root <- multi_root
      multi_root <- paste0(
        original_multi_root,
        "_run_",
        format(Sys.time(), "%Y%m%d_%H%M%S")
      )
      warning(
        "Existing output directory has no ", multi_root_sentinel,
        " sentinel and will NOT be deleted: ", original_multi_root,
        ". Writing this run to: ", multi_root
      )
    } else {
      unlink(multi_root, recursive = TRUE, force = TRUE)
    }
  }
  dir.create(multi_root, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      SEQUENCE_BUILD_VERSION,
      paste0("created=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
      "This file marks the directory as a SEQUENCE output root and permits the",
      "pipeline to clear it on the next run. Do not place other data here."
    ),
    file.path(multi_root, multi_root_sentinel)
  )

  # Provenance for the run as a whole.
  write_sequence_provenance(multi_root, scope = "master")

  # Build banner. Each entry is detected from this file, so a stripped or older
  # script reports FALSE here rather than silently producing fewer figures.
  cat("\n--- SEQUENCE build ---\n")
  cat(sprintf("  script                : %s\n", this_script))
  cat(sprintf("  cutoff figures 1-4    : %s\n", exists("run_sequence_empirical_cutoff_module")))
  cat(sprintf("  dispersion trade-off  : %s\n", exists("save_dispersion_tradeoff_panels")))
  cat(sprintf("  cutoff sensitivity    : %s\n", exists("run_cutoff_sensitivity_module")))
  cat("----------------------\n\n")

  count_path <- if (file.exists(file.path(repo_guess, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"))) {
    file.path(repo_guess, "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
  } else {
    file.path(repo_guess, "data", "WTTS-Seq_2022.2_DE_raw_read_numbers.csv")
  }

  empirical_output_root <- file.path(multi_root, "Empirical_Cutoff_Manuscript")
  cutoff_fit <- run_sequence_empirical_cutoff_module(
    count_path = count_path,
    out_root = empirical_output_root
  )
  cutoff_manifest <- data.frame(
    Cutoff_Method = c("Fixed_5000", "CPM_Empirical", "VST_Empirical"),
    RT0_ZT6 = c(5000L, cutoff_fit$cpm[["RT0_ZT6"]], cutoff_fit$vst[["RT0_ZT6"]]),
    RT2_ZT8 = c(5000L, cutoff_fit$cpm[["RT2_ZT8"]], cutoff_fit$vst[["RT2_ZT8"]]),
    RT4_ZT10 = c(5000L, cutoff_fit$cpm[["RT4_ZT10"]], cutoff_fit$vst[["RT4_ZT10"]]),
    RT8_ZT14 = c(5000L, cutoff_fit$cpm[["RT8_ZT14"]], cutoff_fit$vst[["RT8_ZT14"]]),
    Source = c(
      "Prespecified fixed comparator",
      "Computed internally: log1p(CPM)-EVS + shared c1/c2 + count Pareto + maximum endpoint-chord deviation",
      "Computed internally: DESeq2 VST-EVS + shared c1/c2 + count Pareto + maximum endpoint-chord deviation"
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(cutoff_manifest, file.path(multi_root, "Cutoff_Method_Manifest.csv"), row.names = FALSE)
  utils::write.csv(cutoff_fit$audit, file.path(multi_root, "Empirical_Cutoff_Derivation_Audit.csv"), row.names = FALSE)
  utils::write.csv(
    cutoff_fit$summary_df,
    file.path(multi_root, "Empirical_Cutoff_Manuscript_Summary.csv"),
    row.names = FALSE
  )
  message("Empirical cutoff manuscript figures written to: ", cutoff_fit$figure_dir)
  message("Empirical cutoff manuscript tables written to: ", cutoff_fit$table_dir)

  # Cutoff derivation package. The derivation is complete at this point: it
  # depends only on the count matrix, not on any downstream DESeq2 run, so its
  # figures and tables are archived now rather than waiting for the child runs.
  create_cutoff_derivation_package <- function(multi_root, emp_root) {
    if (!dir.exists(emp_root)) {
      warning("Cutoff derivation package skipped: ", emp_root, " does not exist.")
      return(invisible(NULL))
    }

    # The three derivation CSVs written beside the tree belong with it.
    for (f in c("Cutoff_Method_Manifest.csv",
                "Empirical_Cutoff_Derivation_Audit.csv",
                "Empirical_Cutoff_Manuscript_Summary.csv")) {
      src <- file.path(multi_root, f)
      if (file.exists(src)) {
        file.copy(src, file.path(emp_root, f), overwrite = TRUE)
      }
    }

    n_fig <- length(list.files(file.path(emp_root, "Figures"),
                               pattern = "\\.(pdf|png|tif|tiff)$", ignore.case = TRUE))
    n_tab <- length(list.files(file.path(emp_root, "Tables"),
                               pattern = "\\.csv$", ignore.case = TRUE))
    cat(sprintf("\n  cutoff derivation: %d figure file(s), %d table file(s)\n",
                n_fig, n_tab))
    if (n_fig == 0L) {
      warning("Cutoff derivation produced no figure files in ",
              file.path(emp_root, "Figures"))
    }

    zip_path <- file.path(dirname(multi_root), "SEQUENCE_CUTOFF_DERIVATION.zip")
    if (file.exists(zip_path)) unlink(zip_path, force = TRUE)

    oldwd <- getwd()
    tryCatch({
      setwd(multi_root)
      rel <- list.files(basename(emp_root), recursive = TRUE, all.files = FALSE)
      items <- file.path(basename(emp_root), rel)
      items <- items[file.exists(items)]
      if (length(items)) {
        utils::zip(zipfile = zip_path, files = items, flags = "-q")
      }
    }, finally = {
      setwd(oldwd)
    })

    if (file.exists(zip_path)) {
      cat("  cutoff derivation ZIP: ",
          normalizePath(zip_path, winslash = "/", mustWork = FALSE), "\n\n", sep = "")
    } else {
      warning("Cutoff derivation ZIP creation failed: ", zip_path)
    }
    invisible(zip_path)
  }
  create_cutoff_derivation_package(multi_root, empirical_output_root)

  log_dir <- file.path(multi_root, "Logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  # system2(env = ) is POSIX-only. On Windows the variables are set in the
  # parent process and restored afterwards, which is equivalent for a
  # sequential launcher.
  launch_child_run <- function(method, log_path) {
    child_env <- c(
      SEQUENCE_CHILD_RUN = "1",
      SEQUENCE_CUTOFF_METHOD = method,
      SEQUENCE_MULTI_ROOT = multi_root,
      SEQUENCE_REPO_ROOT = repo_guess
    )

    if (identical(.Platform$OS.type, "windows")) {
      previous <- Sys.getenv(names(child_env), unset = NA_character_, names = TRUE)
      do.call(Sys.setenv, as.list(child_env))
      on.exit({
        restore <- previous[!is.na(previous)]
        if (length(restore)) do.call(Sys.setenv, as.list(restore))
        drop_names <- names(previous)[is.na(previous)]
        if (length(drop_names)) Sys.unsetenv(drop_names)
      }, add = TRUE)

      return(system2(
        command = file.path(R.home("bin"), "Rscript"),
        args = shQuote(this_script),
        stdout = log_path,
        stderr = log_path
      ))
    }

    system2(
      command = file.path(R.home("bin"), "Rscript"),
      args = shQuote(this_script),
      env = paste0(names(child_env), "=", unname(child_env)),
      stdout = log_path,
      stderr = log_path
    )
  }

  methods <- c("fixed5000", "cpm_empirical", "vst_empirical")
  for (method in methods) {
    message("=====================================================")
    message("Starting cutoff method: ", method)
    log_path <- file.path(log_dir, paste0("child_", method, ".log"))
    started_at <- Sys.time()
    status <- launch_child_run(method, log_path)
    message(
      "  finished in ",
      format(round(difftime(Sys.time(), started_at, units = "mins"), 2)),
      "; log: ", log_path
    )

    if (!identical(as.integer(status), 0L)) {
      tail_lines <- tryCatch(utils::tail(readLines(log_path, warn = FALSE), 60L),
                             error = function(e) character(0))
      if (length(tail_lines)) {
        message("----- last 60 lines of ", basename(log_path), " -----")
        message(paste(tail_lines, collapse = "\n"))
        message("--------------------------------------------------")
      }
      stop("Cutoff-method run failed: ", method, " (exit status ", status,
           "). Full log: ", log_path)
    }
  }

  # Cross-cutoff sensitivity: every child has now written its tables, so the
  # three runs can be compared against each other.
  tryCatch(
    run_cutoff_sensitivity_module(multi_root, methods, cutoff_manifest),
    error = function(e) {
      warning("Cutoff sensitivity analysis failed: ", conditionMessage(e))
    }
  )


  # Verify the expected figure set actually reached disk. A missing entry means
  # the stage failed or is absent from this build; either way it is reported
  # here rather than discovered later in the archive.
  verify_expected_outputs <- function(multi_root, methods) {
    checks <- list()
    add <- function(label, path) {
      hit <- length(Sys.glob(path)) > 0L
      checks[[length(checks) + 1L]] <<- data.frame(
        Output = label, Found = hit, stringsAsFactors = FALSE
      )
    }
    emp <- file.path(multi_root, "Empirical_Cutoff_Manuscript", "Figures")
    add("Figure_1_Regime_Definition",       file.path(emp, "Figure_1_Regime_Definition.pdf"))
    add("Figure_2_CPM_EVS_Pareto_Cutoffs",  file.path(emp, "Figure_2_CPM_EVS_Pareto_Cutoffs.pdf"))
    add("Figure_3_VST_EVS_Pareto_Cutoffs",  file.path(emp, "Figure_3_VST_EVS_Pareto_Cutoffs.pdf"))
    add("Figure_4_CPM_vs_VST_EVS_Summary",  file.path(emp, "Figure_4_CPM_vs_VST_EVS_Summary.pdf"))
    add("Figure_Cutoff_Sensitivity",
        file.path(multi_root, "Cutoff_Sensitivity", "Figure_Cutoff_Sensitivity.pdf"))
    add("SEQUENCE_CUTOFF_DERIVATION.zip",
        file.path(dirname(multi_root), "SEQUENCE_CUTOFF_DERIVATION.zip"))
    for (m in methods) {
      add(paste0(m, ": Figure_Dispersion_Tradeoff"),
          file.path(multi_root, m, "Combined_Figures", "Figure_Dispersion_Tradeoff_*.pdf"))
    }
    out <- do.call(rbind, checks)
    cat("\n--- expected figure outputs ---\n")
    for (i in seq_len(nrow(out))) {
      cat(sprintf("  [%s] %s\n", if (out$Found[i]) "ok " else "MISS", out$Output[i]))
    }
    cat("-------------------------------\n\n")
    if (any(!out$Found)) {
      warning("Some expected figures were not produced: ",
              paste(out$Output[!out$Found], collapse = "; "))
    }
    utils::write.csv(out, file.path(multi_root, "Figure_Output_Check.csv"),
                     row.names = FALSE)
    invisible(out)
  }
  verify_expected_outputs(multi_root, methods)

  # Manuscript figure package: the derivation figures, the cross-cutoff
  # sensitivity figure and the per-cutoff dispersion figures collected into one
  # archive, so the figures that go into the paper are not scattered across
  # three run trees.
  create_manuscript_figure_package <- function(multi_root, methods) {
    stage <- file.path(multi_root, "Manuscript_Figures")
    unlink(stage, recursive = TRUE, force = TRUE)
    dir.create(file.path(stage, "Cutoff_Derivation"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(stage, "Cutoff_Sensitivity"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(stage, "Dispersion_Tradeoff"), recursive = TRUE, showWarnings = FALSE)

    copy_glob <- function(pattern, dest, prefix = "") {
      hits <- Sys.glob(pattern)
      if (!length(hits)) return(0L)
      ok <- vapply(hits, function(f) {
        file.copy(f, file.path(dest, paste0(prefix, basename(f))), overwrite = TRUE)
      }, logical(1))
      sum(ok)
    }

    emp <- file.path(multi_root, "Empirical_Cutoff_Manuscript")
    copy_glob(file.path(emp, "Figures", "Figure_*.pdf"), file.path(stage, "Cutoff_Derivation"))
    copy_glob(file.path(emp, "Figures", "Figure_*.png"), file.path(stage, "Cutoff_Derivation"))
    copy_glob(file.path(emp, "Figure_Legends.md"),       file.path(stage, "Cutoff_Derivation"))
    copy_glob(file.path(emp, "Methods_Manuscript.md"),   file.path(stage, "Cutoff_Derivation"))

    sens <- file.path(multi_root, "Cutoff_Sensitivity")
    copy_glob(file.path(sens, "Figure_*.pdf"), file.path(stage, "Cutoff_Sensitivity"))
    copy_glob(file.path(sens, "Figure_*.png"), file.path(stage, "Cutoff_Sensitivity"))
    copy_glob(file.path(sens, "Table_*.csv"),  file.path(stage, "Cutoff_Sensitivity"))

    # Dispersion figures carry the cutoff method in their name, since one
    # exists per run and the file names are otherwise identical.
    for (m in methods) {
      copy_glob(
        file.path(multi_root, m, "Combined_Figures", "Figure_Dispersion_Tradeoff_*"),
        file.path(stage, "Dispersion_Tradeoff"),
        prefix = paste0(m, "_")
      )
    }

    copy_glob(file.path(multi_root, "Cutoff_Method_Manifest.csv"), stage)
    copy_glob(file.path(multi_root, "Figure_Output_Check.csv"),    stage)

    n_files <- length(list.files(stage, recursive = TRUE))
    if (!n_files) {
      warning("Manuscript figure package is empty; nothing was copied.")
      return(invisible(NULL))
    }

    fig_zip <- file.path(dirname(multi_root), "SEQUENCE_MANUSCRIPT_FIGURES.zip")
    if (file.exists(fig_zip)) unlink(fig_zip, force = TRUE)
    oldwd <- getwd()
    tryCatch({
      setwd(dirname(multi_root))
      rel <- list.files(stage, recursive = TRUE)
      items <- file.path(basename(multi_root), "Manuscript_Figures", rel)
      items <- items[file.exists(items)]
      utils::zip(zipfile = basename(fig_zip), files = items, flags = "-q")
    }, finally = {
      setwd(oldwd)
    })

    if (file.exists(fig_zip)) {
      cat(sprintf("Manuscript figure package: %s (%d files)\n",
                  normalizePath(fig_zip, winslash = "/", mustWork = FALSE), n_files))
    } else {
      warning("Manuscript figure ZIP creation failed: ", fig_zip)
    }
    invisible(fig_zip)
  }
  create_manuscript_figure_package(multi_root, methods)

  final_zip <- file.path(dirname(multi_root), "SEQUENCE_FINAL_ALL_CUTOFF_METHODS.zip")
  if (file.exists(final_zip)) unlink(final_zip, force = TRUE)
  oldwd <- getwd()
  tryCatch({
    setwd(dirname(multi_root))
    zip_rel <- list.files(multi_root, recursive = TRUE, all.files = FALSE)
    zip_rel <- zip_rel[!grepl("SEQUENCE_(ALL_(FIGURES|TABLES)|MANUSCRIPT_FIGURES|CUTOFF_DERIVATION)\\.zip$", zip_rel)]
    zip_items <- file.path(basename(multi_root), zip_rel)
    zip_items <- zip_items[file.exists(zip_items)]
    utils::zip(zipfile = basename(final_zip), files = zip_items, flags = "-q")
  }, finally = {
    setwd(oldwd)
  })
  if (!file.exists(final_zip)) stop("Final multi-cutoff ZIP creation failed: ", final_zip)

  cat("\n=====================================================\n")
  cat("All three cutoff-method runs complete (Fixed 5,000, CPM-EVS empirical, VST-EVS empirical).\n")
  cat("Final combined ZIP:\n", normalizePath(final_zip, winslash = "/", mustWork = TRUE), "\n", sep = "")
  cut_zip <- file.path(dirname(multi_root), "SEQUENCE_CUTOFF_DERIVATION.zip")
  if (file.exists(cut_zip)) {
    cat("Cutoff derivation ZIP:\n", normalizePath(cut_zip, winslash = "/", mustWork = FALSE), "\n", sep = "")
  }
  man_zip <- file.path(dirname(multi_root), "SEQUENCE_MANUSCRIPT_FIGURES.zip")
  if (file.exists(man_zip)) {
    cat("Manuscript figures ZIP:\n", normalizePath(man_zip, winslash = "/", mustWork = FALSE), "\n", sep = "")
  }
  cat("=====================================================\n\n")
  quit(save = "no", status = 0L)
}

# Active child-run cutoff method.
CUTOFF_METHOD_KEY <- Sys.getenv("SEQUENCE_CUTOFF_METHOD", unset = "fixed5000")
CUTOFF_METHOD_INFO <- switch(
  CUTOFF_METHOD_KEY,
  fixed5000 = list(slug = "Fixed_5000", label = "Fixed 5,000", basis = "Prespecified fixed top-5,000 PASs per condition"),
  cpm_empirical = list(slug = "CPM_Empirical", label = "CPM empirical k*", basis = "Internally computed log1p(CPM)-EVS count-Pareto maximum endpoint-chord cutoff"),
  vst_empirical = list(slug = "VST_Empirical", label = "VST empirical k*", basis = "Internally computed DESeq2 VST-EVS count-Pareto maximum endpoint-chord cutoff"),
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

# Re-seed inside the child process. The master seeds before relaunching, but a
# child is a fresh R session and must set its own RNG state.
set.seed(SEQUENCE_SEED)

message("=====================================================")
message(SEQUENCE_BUILD_VERSION)
message("Child run: cutoff method = ", CUTOFF_METHOD_LABEL,
        " (", CUTOFF_METHOD_SLUG, ")")
message("R ", getRversion(), " | cairo: ", isTRUE(capabilities("cairo")),
        " | patchwork: ", requireNamespace("patchwork", quietly = TRUE),
        " | ggrastr: ", requireNamespace("ggrastr", quietly = TRUE))
message("=====================================================")

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
EMPIRICAL_EVS_CUTOFF_BASIS <- CUTOFF_METHOD_INFO$basis

# 600 dpi is the usual minimum for combined line-art/raster figures. Because
# canvases are now specified at final print size, this is the true output
# resolution rather than a value that will be re-scaled in production.
figure_dpi <- 600

# Base text size in POINTS AT FINAL PRINT SIZE. Derived sizes elsewhere in the
# script are expressed as base_theme_size - 1, - 2 and - 3, so 8 pt keeps the
# smallest annotation at 5 pt, which is the practical legibility floor.
base_theme_size <- 8

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
# Scope of the empirical null
# -----------------------------------------------------------------------------
# fdrtool is fit to every finite DESeq2 Wald statistic returned for the view,
# including PASs with NA BH-adjusted p-values after DESeq2 independent filtering.

# -----------------------------------------------------------------------------
# Palette
# -----------------------------------------------------------------------------

# Okabe-Ito colour-blind-safe palette used consistently across significance figures.
plot_palette <- list(
  background = "#9E9E9E",   # non-significant PASs, darker so they survive print
  threshold  = "#8C6D31",   # +/- LFC boundary rules
  hc         = "#B8860B",   # higher-criticism p threshold rule
  hbfss_line = "#000000",   # HBFSS boundary curve, deliberately not a point colour
  standard   = "#0072B2",   # Okabe-Ito blue
  strong     = "#D55E00",   # Okabe-Ito vermillion
  weak       = "#009E73",   # Okabe-Ito bluish green
  hbfss      = "#CC79A7",   # Okabe-Ito reddish purple
  overlap    = "#56B4E9",   # Okabe-Ito sky blue
  control    = "#4D4D4D",
  treatment  = "#0072B2",
  histogram  = "#969696"
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
  raw_evs = "RawEVS",
  cpm_evs = "CPMEVS",
  vst_evs = "VSTEVS"
)

active_evs_track_keys <- function() {
  base_tracks <- c("normalized_evs", "raw_evs")

  if (!isTRUE(RUN_MATCHED_TRANSFORM_TRACKS_ONLY)) {
    return(c(base_tracks, "cpm_evs", "vst_evs"))
  }

  if (identical(CUTOFF_METHOD_KEY, "cpm_empirical")) {
    return(c(base_tracks, "cpm_evs"))
  }

  if (identical(CUTOFF_METHOD_KEY, "vst_empirical")) {
    return(c(base_tracks, "vst_evs"))
  }

  base_tracks
}

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
required_manifest_methods <- c("Fixed_5000", "CPM_Empirical", "VST_Empirical")
if (!all(required_manifest_methods %in% cutoff_manifest_runtime$Cutoff_Method)) {
  stop("Cutoff manifest does not contain all required methods: ", paste(required_manifest_methods, collapse=", "))
}
manifest_row <- function(method) cutoff_manifest_runtime[match(method, cutoff_manifest_runtime$Cutoff_Method),,drop=FALSE]
fixed_row <- manifest_row("Fixed_5000")
cpm_row <- manifest_row("CPM_Empirical")
vst_row <- manifest_row("VST_Empirical")

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix   = c("R0", "R2", "R4", "R8"),
  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),
  fixed_5000_k    = as.integer(unlist(fixed_row[c("RT0_ZT6","RT2_ZT8","RT4_ZT10","RT8_ZT14")], use.names=FALSE)),
  cpm_empirical_k = as.integer(unlist(cpm_row[c("RT0_ZT6","RT2_ZT8","RT4_ZT10","RT8_ZT14")], use.names=FALSE)),
  vst_empirical_k = as.integer(unlist(vst_row[c("RT0_ZT6","RT2_ZT8","RT4_ZT10","RT8_ZT14")], use.names=FALSE)),
  stringsAsFactors = FALSE
)

enabled_comparisons <- names(RUN_COMPARISON_FLAGS)[RUN_COMPARISON_FLAGS]
comparison_table <- comparison_table[
  comparison_table$comparison_name %in% enabled_comparisons,
  ,
  drop = FALSE
]

if (!nrow(comparison_table)) {
  stop("No downstream comparisons are enabled by RUN_COMPARISON_FLAGS.")
}

message(
  "Downstream comparison switch: running ",
  paste(comparison_table$comparison_name, collapse = ", "),
  "; disabled: ",
  paste(setdiff(names(RUN_COMPARISON_FLAGS), comparison_table$comparison_name), collapse = ", ")
)

get_active_evs_cutoff <- function(comparison_name) {
  idx <- match(as.character(comparison_name), comparison_table$comparison_name)
  if (is.na(idx)) stop("No EVS cutoff is defined for comparison: ", comparison_name)
  k <- switch(
    CUTOFF_METHOD_KEY,
    fixed5000 = comparison_table$fixed_5000_k[idx],
    cpm_empirical = comparison_table$cpm_empirical_k[idx],
    vst_empirical = comparison_table$vst_empirical_k[idx]
  )
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || !is.finite(k) || k < 1L) {
    stop("Invalid EVS cutoff for ", comparison_name, " under ", CUTOFF_METHOD_LABEL)
  }
  k
}

if (anyDuplicated(comparison_table$comparison_name)) stop("comparison_table contains duplicated comparison names.")
cutoff_cols <- c("fixed_5000_k", "cpm_empirical_k", "vst_empirical_k")
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
  unset = file.path(repo_root, "exports", "sequence_final_three_cutoffs")
)
output_dir <- file.path(multi_output_root, CUTOFF_METHOD_SLUG)

if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE, force = TRUE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paper_fig_dir <- file.path(output_dir, "Combined_Figures")
summary_table_dir <- file.path(output_dir, "Summary_Tables")
dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_table_dir, recursive = TRUE, showWarnings = FALSE)

# Per-child provenance, so each cutoff-method output tree is self-describing.
write_sequence_provenance(output_dir, scope = CUTOFF_METHOD_SLUG)

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

# -----------------------------------------------------------------------------
# Figure device layer
# -----------------------------------------------------------------------------
# All figure output is routed through cairo devices. This gives antialiased
# raster output, embedded PDF fonts, and PNG
# and PDF renditions of the same figure that actually match.
figure_open_device <- function(path, width_in, height_in, dpi = figure_dpi, bg = "white") {
  ext <- tolower(tools::file_ext(path))
  cairo_ok <- isTRUE(capabilities("cairo"))

  if (ext == "pdf") {
    if (cairo_ok) {
      grDevices::cairo_pdf(path, width = width_in, height = height_in, bg = bg)
    } else {
      grDevices::pdf(path, width = width_in, height = height_in, bg = bg,
                     useDingbats = FALSE)
    }
  } else if (ext == "png") {
    if (cairo_ok) {
      grDevices::png(path, width = width_in, height = height_in, units = "in",
                     res = dpi, bg = bg, type = "cairo")
    } else {
      grDevices::png(path, width = width_in, height = height_in, units = "in",
                     res = dpi, bg = bg)
    }
  } else if (ext %in% c("tif", "tiff")) {
    if (cairo_ok) {
      grDevices::tiff(path, width = width_in, height = height_in, units = "in",
                      res = dpi, bg = bg, compression = "lzw", type = "cairo")
    } else {
      grDevices::tiff(path, width = width_in, height = height_in, units = "in",
                      res = dpi, bg = bg, compression = "lzw")
    }
  } else {
    stop("Unsupported figure extension: ", ext)
  }

  invisible(TRUE)
}

# Every figure written is recorded, so the submission package carries a machine
# readable inventory of canvas geometry and resolution.
figure_manifest_path <- function() {
  file.path(output_dir, "Figure_Manifest.csv")
}

register_written_figure <- function(path, width_in, height_in, dpi) {
  row <- data.frame(
    file = normalizePath(path, winslash = "/", mustWork = FALSE),
    format = tolower(tools::file_ext(path)),
    width_mm = round(width_in * 25.4, 1),
    height_mm = round(height_in * 25.4, 1),
    width_in = round(width_in, 3),
    height_in = round(height_in, 3),
    dpi = dpi,
    build = SEQUENCE_BUILD_VERSION,
    cutoff_method = if (exists("CUTOFF_METHOD_SLUG")) CUTOFF_METHOD_SLUG else NA_character_,
    written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )

  manifest <- figure_manifest_path()
  utils::write.table(
    row,
    file = manifest,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(manifest),
    append = file.exists(manifest),
    qmethod = "double"
  )

  invisible(row)
}

# Canvas dimensions are given at final print size; `units` may be "mm" or "in".
save_grob <- function(g, path, width = FIG_DOUBLE_COL_MM, height = 150,
                      dpi = figure_dpi, bg = "white", units = c("mm", "in")) {
  if (is.null(g)) {
    return(invisible(NULL))
  }
  if (!should_write_figure(path)) {
    return(invisible(NULL))
  }

  units <- match.arg(units)
  width_in <- if (identical(units, "mm")) fig_mm2in(width) else width
  height_in <- if (identical(units, "mm")) fig_mm2in(height) else height

  if (!is.finite(width_in) || !is.finite(height_in) ||
      width_in <= 0 || height_in <= 0) {
    stop("Invalid figure canvas for ", path)
  }

  # Guard against the oversized canvases that made production down-scaling
  # illegible. This is a warning, not an error, so an intentionally tall
  # supplementary figure still gets written.
  if (width_in > fig_mm2in(FIG_DOUBLE_COL_MM) + 1e-6) {
    warning(
      sprintf(
        "Figure wider than the %.0f mm double-column limit (%.0f mm): %s",
        FIG_DOUBLE_COL_MM, width_in * 25.4, basename(path)
      ),
      call. = FALSE
    )
  }
  if (height_in > fig_mm2in(FIG_MAX_HEIGHT_MM) + 1e-6) {
    warning(
      sprintf(
        "Figure taller than the %.0f mm page limit (%.0f mm): %s",
        FIG_MAX_HEIGHT_MM, height_in * 25.4, basename(path)
      ),
      call. = FALSE
    )
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  figure_open_device(path, width_in, height_in, dpi = dpi, bg = bg)
  drawn <- tryCatch({
    if (inherits(g, "ggplot")) {
      print(g)
    } else {
      grid::grid.newpage()
      grid::grid.draw(g)
    }
    TRUE
  }, error = function(e) {
    message("  ERROR drawing ", basename(path), ": ", conditionMessage(e))
    FALSE
  })
  grDevices::dev.off()

  if (!drawn) {
    unlink(path, force = TRUE)
    return(invisible(NULL))
  }

  register_written_figure(path, width_in, height_in, dpi)
  invisible(path)
}

# Write each figure to every configured output format from one canvas specification.
save_figure <- function(g, path_base, width, height, dpi = figure_dpi,
                        bg = "white", units = c("mm", "in")) {
  units <- match.arg(units)
  formats <- c("pdf", "png", if (isTRUE(FIG_WRITE_TIFF)) "tiff")

  for (fmt in formats) {
    save_grob(
      g,
      paste0(path_base, ".", fmt),
      width = width,
      height = height,
      dpi = dpi,
      bg = bg,
      units = units
    )
  }

  invisible(path_base)
}

# Rasterises a dense point layer when ggrastr is available. A volcano with tens
# of thousands of non-significant PASs produces a vector PDF that is enormous
# and effectively uneditable; rasterising only that layer keeps axes, rules,
# significant markers and text fully vector.
maybe_rasterise <- function(layer, dpi = figure_dpi) {
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    return(ggrastr::rasterise(layer, dpi = dpi))
  }
  layer
}

# Display form of a comparison key: "RT0_ZT6" -> "RT0 vs ZT6". The underscored
# key is used for file paths and list names; figures show the readable form.
pretty_comparison <- function(x) sub("_", " vs ", x, fixed = TRUE)

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
    return(paste(pretty_comparison(comparison_name), track_label,
                 pretty_dataset_type(dataset_key), sep = " | "))
  }

  dataset_key <- paste(parts[3:length(parts)], collapse = "_")

  paste(pretty_comparison(comparison_name), pretty_dataset_type(dataset_key), sep = " | ")
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
  n_search <- max(1L, min(n_total, floor(HC_ALPHA0 * n_total)))
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
    stop(
      sprintf(
        "[%s] Internal HC audit failed: explicit lower-tail HCp %.17g != fdrtool::hc.thresh(alpha0=%.3f) %.17g.",
        dataset_name,
        hc_p,
        HC_ALPHA0,
        package_hc_p
      )
    )
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

# Marker sizes used consistently across significance figures.
significance_method_sizes <- c(
  "Standard" = 1.55,
  "Strong" = 1.45,
  "Weak" = 1.55,
  "HBFSS" = 1.75
)

significance_method_colors <- c(
  "Standard" = plot_palette$standard,
  "Strong" = plot_palette$strong,
  "Weak" = plot_palette$weak,
  "HBFSS" = plot_palette$hbfss
)

# Distinct marker shapes for Standard, Strong, Weak, and HBFSS.
significance_method_shapes <- c(
  "Standard" = 16,
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
  theme_bw(base_size = base_theme_size, base_family = FIG_FONT_FAMILY) +
    theme(
      text = element_text(family = FIG_FONT_FAMILY),
      plot.tag = element_text(
        face = "bold",
        size = base_theme_size + 1,
        family = FIG_FONT_FAMILY
      ),
      plot.tag.position = "topleft",
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
      # Horizontal reference lines only. A full grid competes with the point
      # cloud in a dense volcano once the panel is printed at 60-90 mm wide.
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"),
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.3, colour = "grey35"),
      plot.margin = margin(3, 4, 3, 3)
    )
}

# Extracts the shared method legend from an assembled panel. ggplot2 < 3.5.0
# emits one gtable grob named "guide-box"; ggplot2 >= 3.5.0 emits
# position-specific boxes ("guide-box-bottom", "guide-box-right", and so on)
# either ggplot2 generation.
shared_panel_legend <- function(plot_obj) {
  g <- ggplotGrob(plot_obj + theme(legend.position = "bottom"))

  grob_names <- vapply(
    g$grobs,
    function(x) if (is.null(x$name)) "" else as.character(x$name)[1],
    character(1)
  )

  candidates <- which(startsWith(grob_names, "guide-box"))
  if (!length(candidates)) {
    warning("No guide-box grob found; panel legend will be omitted.", call. = FALSE)
    return(NULL)
  }

  has_content <- vapply(
    candidates,
    function(i) {
      gb <- g$grobs[[i]]
      if (inherits(gb, "zeroGrob")) return(FALSE)
      kids <- tryCatch(length(gb$grobs), error = function(e) 0L)
      isTRUE(kids > 0L)
    },
    logical(1)
  )

  keep <- candidates[has_content]
  if (!length(keep)) {
    warning("Only empty guide-box placeholders found; panel legend will be omitted.",
            call. = FALSE)
    return(NULL)
  }

  g$grobs[[keep[1]]]
}

# Build the shared legend from one synthetic row per method so all four markers are shown.
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

# Returns the number of columns and rows a panel set will actually occupy.
# Exported so figure canvases can be sized from the real grid rather than from
# the requested column count.
panel_grid_dim <- function(n_plots, ncol_requested = n_plots) {
  n_plots <- max(1L, as.integer(n_plots))
  ncol_use <- max(1L, min(as.integer(ncol_requested), PANEL_MAX_COLS, n_plots))
  c(ncol = ncol_use, nrow = as.integer(ceiling(n_plots / ncol_use)))
}

assemble_one_legend_panel <- function(plot_list, panel_title, ncol = length(plot_list), width_legend = TRUE) {
  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0L) {
    return(NULL)
  }

  dims <- panel_grid_dim(length(plot_list), ncol)
  ncol_use <- unname(dims[["ncol"]])

  legend <- build_full_method_legend()

  # Panel tags are applied per sub-plot rather than by the compositor, so they
  # survive both the patchwork and the gridExtra assembly paths.
  tagged <- lapply(
    seq_along(plot_list),
    function(i) {
      strip_legend(plot_list[[i]]) +
        labs(tag = LETTERS[((i - 1L) %% 26L) + 1L]) +
        theme(
          plot.tag = element_text(
            face = "bold",
            size = base_theme_size + 1,
            family = FIG_FONT_FAMILY
          ),
          plot.tag.position = "topleft"
        )
    }
  )

  # patchwork aligns panel interiors across rows and columns even when the
  # sub-plots have axis labels of different widths; gridExtra does not. The
  # patchwork object is converted straight to a gtable so the rest of the
  # assembly (legend strip, title) is identical on both paths.
  row <- if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::patchworkGrob(
      patchwork::wrap_plots(tagged, ncol = ncol_use)
    )
  } else {
    do.call(
      gridExtra::arrangeGrob,
      c(tagged, list(ncol = ncol_use))
    )
  }

  if (is.null(legend)) {
    return(
      gridExtra::arrangeGrob(
        row,
        ncol = 1,
        top = grid::textGrob(
          panel_title,
          gp = grid::gpar(
            fontface = "bold",
            fontsize = base_theme_size + 2,
            fontfamily = if (nzchar(FIG_FONT_FAMILY)) FIG_FONT_FAMILY else ""
          )
        )
      )
    )
  }

  gridExtra::arrangeGrob(
    row,
    legend,
    ncol = 1,
    heights = grid::unit.c(grid::unit(1, "null"), grid::unit(9, "mm")),
    top = grid::textGrob(
      panel_title,
      gp = grid::gpar(
        fontface = "bold",
        fontsize = base_theme_size + 2,
        fontfamily = if (nzchar(FIG_FONT_FAMILY)) FIG_FONT_FAMILY else ""
      )
    )
  )
}

artifact_relative_path <- function(path, root) {
  root_norm <- normalizePath(root, winslash = "/", mustWork = TRUE)
  path_norm <- normalizePath(path, winslash = "/", mustWork = TRUE)
  substring(path_norm, nchar(root_norm) + 2L)
}

copy_artifact_preserving_path <- function(source_path, source_root, destination_root) {
  rel <- artifact_relative_path(source_path, source_root)
  destination <- file.path(destination_root, rel)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(source_path, destination, overwrite = TRUE, copy.date = TRUE)) {
    stop("Failed to copy artifact: ", source_path)
  }
  destination
}

create_comparison_artifact_package <- function(comparison_name) {
  comparison_dir <- file.path(output_dir, comparison_name)
  if (!dir.exists(comparison_dir)) stop("Comparison output directory is missing: ", comparison_dir)

  figures_dir <- file.path(comparison_dir, "Figures")
  tables_dir <- file.path(comparison_dir, "Tables")
  if (dir.exists(figures_dir)) unlink(figures_dir, recursive = TRUE, force = TRUE)
  if (dir.exists(tables_dir)) unlink(tables_dir, recursive = TRUE, force = TRUE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  source_files <- list.files(comparison_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  source_files <- source_files[file.info(source_files)$isdir %in% FALSE]
  normalized <- gsub("\\\\", "/", source_files)
  source_files <- source_files[!grepl("/(Figures|Tables)/", normalized)]
  source_files <- source_files[!grepl("\\.zip$", source_files, ignore.case = TRUE)]

  fig_files <- source_files[grepl("\\.(png|pdf|tif|tiff)$", source_files, ignore.case = TRUE)]
  tab_files <- source_files[grepl("\\.(csv|tsv|txt)$", source_files, ignore.case = TRUE)]

  if (length(fig_files)) {
    invisible(vapply(fig_files, copy_artifact_preserving_path, character(1), source_root = comparison_dir, destination_root = figures_dir))
  }
  if (length(tab_files)) {
    invisible(vapply(tab_files, copy_artifact_preserving_path, character(1), source_root = comparison_dir, destination_root = tables_dir))
  }
  if (!length(fig_files) && !length(tab_files)) stop("No tables or figures found for ", comparison_name)

  zip_path <- file.path(comparison_dir, paste0(comparison_name, ".zip"))
  package_items <- c()
  if (length(list.files(figures_dir, recursive = TRUE))) package_items <- c(package_items, "Figures")
  if (length(list.files(tables_dir, recursive = TRUE))) package_items <- c(package_items, "Tables")
  if (file.exists(zip_path)) unlink(zip_path, force = TRUE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(comparison_dir)
  utils::zip(zipfile = basename(zip_path), files = package_items, flags = "-rq")
  if (!file.exists(zip_path)) stop("Comparison ZIP creation failed: ", zip_path)
  normalizePath(zip_path, winslash = "/", mustWork = TRUE)
}

create_cutoff_artifact_package <- function() {
  figures_dir <- file.path(output_dir, "Figures")
  tables_dir <- file.path(output_dir, "Tables")
  if (dir.exists(figures_dir)) unlink(figures_dir, recursive = TRUE, force = TRUE)
  if (dir.exists(tables_dir)) unlink(tables_dir, recursive = TRUE, force = TRUE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  artifact_roots <- c(
    file.path(output_dir, "Combined_Figures"),
    file.path(output_dir, "Summary_Tables"),
    file.path(output_dir, "TWAS"),
    file.path(output_dir, as.character(comparison_table$comparison_name))
  )
  artifact_roots <- artifact_roots[dir.exists(artifact_roots)]

  source_files <- unlist(lapply(
    artifact_roots,
    function(root) list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  ), use.names = FALSE)
  source_files <- unique(source_files[file.info(source_files)$isdir %in% FALSE])
  normalized <- gsub("\\\\", "/", source_files)
  source_files <- source_files[!grepl("/(Figures|Tables)/", normalized)]
  source_files <- source_files[!grepl("\\.zip$", source_files, ignore.case = TRUE)]

  fig_files <- source_files[grepl("\\.(png|pdf|tif|tiff)$", source_files, ignore.case = TRUE)]
  tab_files <- source_files[grepl("\\.(csv|tsv|txt)$", source_files, ignore.case = TRUE)]

  if (length(fig_files)) {
    invisible(vapply(fig_files, copy_artifact_preserving_path, character(1), source_root = output_dir, destination_root = figures_dir))
  }
  if (length(tab_files)) {
    invisible(vapply(tab_files, copy_artifact_preserving_path, character(1), source_root = output_dir, destination_root = tables_dir))
  }
  if (!length(fig_files) && !length(tab_files)) stop("No cutoff-level tables or figures found for ", CUTOFF_METHOD_SLUG)

  cutoff_zip <- file.path(output_dir, paste0(CUTOFF_METHOD_SLUG, ".zip"))
  if (file.exists(cutoff_zip)) unlink(cutoff_zip, force = TRUE)
  package_items <- c("Figures", "Tables")
  if (file.exists(file.path(output_dir, "Methods_Manuscript.md"))) {
    package_items <- c(package_items, "Methods_Manuscript.md")
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)
  utils::zip(zipfile = basename(cutoff_zip), files = package_items, flags = "-rq")
  if (!file.exists(cutoff_zip)) stop("Cutoff ZIP creation failed: ", cutoff_zip)
  normalizePath(cutoff_zip, winslash = "/", mustWork = TRUE)
}

create_all_figures_zip <- function() {
  # Collect manuscript PNG panels for the figure archive.
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
  # Collect manuscript summary tables for the table archive.
  wanted <- character(0)

  add_if_exists <- function(path) {
    if (file.exists(path)) wanted <<- c(wanted, path)
  }

  add_if_exists(file.path(summary_table_dir, "Table_DE_Method_Counts.csv"))
  add_if_exists(file.path(summary_table_dir, "Table_EVS_Cutoffs.csv"))
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
# Lazy experiment-wide CPM/VST matrices for the added matched EVS pathways
# -----------------------------------------------------------------------------
# These reproduce the preprocessing used by the empirical cutoff engine.
# They are used ONLY for PCA/PC1 ranking and Lead/Rem membership. The selected
# feature IDs are always applied back to the raw integer count matrix before
# downstream DESeq2/HBFSS/TWAS testing.
.evs_transform_cache <- new.env(parent = emptyenv())

full_wtts_raw_counts <- function() {
  x <- as.matrix(WTTS_Seq[, meta_all$id, drop = FALSE])
  storage.mode(x) <- "numeric"
  rownames(x) <- rownames(WTTS_Seq)
  x <- x[rowSums(x) > 0, , drop = FALSE]
  x
}

get_experimentwide_cpm_evs_matrix <- function() {
  if (exists("cpm", envir = .evs_transform_cache, inherits = FALSE)) {
    return(get("cpm", envir = .evs_transform_cache, inherits = FALSE))
  }

  raw <- full_wtts_raw_counts()
  lib <- colSums(raw)
  lib[!is.finite(lib) | lib <= 0] <- 1
  cpm_mat <- log1p(sweep(raw, 2, lib / 1e6, "/"))

  if (any(!is.finite(cpm_mat))) {
    stop("Experiment-wide CPM-EVS matrix contains non-finite values.")
  }

  assign("cpm", cpm_mat, envir = .evs_transform_cache)
  cpm_mat
}

get_experimentwide_vst_evs_matrix <- function() {
  if (exists("vst", envir = .evs_transform_cache, inherits = FALSE)) {
    return(get("vst", envir = .evs_transform_cache, inherits = FALSE))
  }

  raw <- full_wtts_raw_counts()

  time_group <- rep(NA_character_, ncol(raw))
  names(time_group) <- colnames(raw)
  group_patterns <- c(
    RT0 = "^R0_", ZT6 = "^ZT6_",
    RT2 = "^R2_", ZT8 = "^ZT8_",
    RT4 = "^R4_", ZT10 = "^ZT10_",
    RT8 = "^R8_", ZT14 = "^ZT14_"
  )

  for (nm in names(group_patterns)) {
    idx <- grep(group_patterns[[nm]], colnames(raw))
    if (any(!is.na(time_group[idx]))) {
      stop("A sample matched more than one VST time-group pattern.")
    }
    time_group[idx] <- nm
  }

  if (anyNA(time_group)) {
    stop(
      "Unassigned samples in experiment-wide VST-EVS matrix: ",
      paste(names(time_group)[is.na(time_group)], collapse = ", ")
    )
  }

  time_group <- factor(time_group, levels = names(group_patterns))

  dds_vst <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(raw),
    colData = data.frame(group = time_group, row.names = colnames(raw)),
    design = ~ group
  )

  dds_vst <- tryCatch(
    DESeq2::estimateSizeFactors(dds_vst),
    error = function(e) DESeq2::estimateSizeFactors(dds_vst, type = "poscounts")
  )

  vst_obj <- DESeq2::varianceStabilizingTransformation(dds_vst, blind = TRUE)
  vst_mat <- SummarizedExperiment::assay(vst_obj)

  if (!identical(dim(vst_mat), dim(raw)) ||
      !identical(rownames(vst_mat), rownames(raw)) ||
      !identical(colnames(vst_mat), colnames(raw))) {
    stop("Experiment-wide VST-EVS matrix does not match the raw count matrix.")
  }

  if (any(!is.finite(vst_mat))) {
    stop("Experiment-wide VST-EVS matrix contains non-finite values.")
  }

  assign("vst", vst_mat, envir = .evs_transform_cache)
  vst_mat
}


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

make_rank_matrix_for_track <- function(count_matrix, coldata, track_key) {
  track_key <- match.arg(
    track_key,
    c("normalized_evs", "raw_evs", "cpm_evs", "vst_evs")
  )

  if (track_key == "raw_evs") {
    return(
      list(
        rank_matrix = as.data.frame(count_matrix),
        preprocessing_label = "Raw counts prior to EVS; PCA is performed directly on raw counts"
      )
    )
  }

  if (track_key == "cpm_evs") {
    full_cpm <- get_experimentwide_cpm_evs_matrix()
    missing_rows <- setdiff(rownames(count_matrix), rownames(full_cpm))
    missing_cols <- setdiff(colnames(count_matrix), colnames(full_cpm))
    if (length(missing_rows) || length(missing_cols)) {
      stop("CPM-EVS rank matrix is missing comparison features or samples.")
    }

    cpm_sub <- full_cpm[
      rownames(count_matrix),
      colnames(count_matrix),
      drop = FALSE
    ]

    return(
      list(
        rank_matrix = as.data.frame(cpm_sub),
        preprocessing_label = "Experiment-wide log1p(CPM) prior to EVS; PCA is performed directly on CPM-transformed values"
      )
    )
  }

  if (track_key == "vst_evs") {
    full_vst <- get_experimentwide_vst_evs_matrix()
    missing_rows <- setdiff(rownames(count_matrix), rownames(full_vst))
    missing_cols <- setdiff(colnames(count_matrix), colnames(full_vst))
    if (length(missing_rows) || length(missing_cols)) {
      stop("VST-EVS rank matrix is missing comparison features or samples.")
    }

    vst_sub <- full_vst[
      rownames(count_matrix),
      colnames(count_matrix),
      drop = FALSE
    ]

    return(
      list(
        rank_matrix = as.data.frame(vst_sub),
        preprocessing_label = "Experiment-wide DESeq2 VST (blind=TRUE) prior to EVS; PCA is performed directly on VST values"
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
    preprocessing_label = "Median-of-ratios normalized counts prior to EVS; PCA is performed directly on normalized counts"
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
  track_key <- match.arg(track_key, c("normalized_evs", "raw_evs", "cpm_evs", "vst_evs"))
  empirical_k <- get_active_evs_cutoff(comparison_name)

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
    raw_dataset = raw_count_matrix,
    leading_edge_dataset = raw_count_matrix[leading_edge_ids, , drop = FALSE],
    remainder_dataset = raw_count_matrix[remainder_ids, , drop = FALSE]
  )
}

run_core_analysis <- function(count_mat, coldata, dataset_name, annot_df) {
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

  # The y-axis limit is shared across every panel, so some PASs can fall above
  # it. coord_cartesian() clipped them silently; the count is now reported.
  n_above_axis <- sum(
    is.finite(plot_df$neglog10_empirical_p) &
      plot_df$neglog10_empirical_p > y_limit,
    na.rm = TRUE
  )

  # ASCII only. A literal Unicode tau in a caption string is encoding-fragile
  # across graphics devices and locales, and rendered as a missing glyph on
  # non-cairo PDF devices.
  count_text <- paste0(
    "Std=", std_n,
    "  Str=", str_n,
    "  Wk=", wk_n,
    "  HBFSS=", hbfss_n,
    "  Ovlp=", ovlp_n,
    if (is.finite(hc_raw) && !is.na(hc_raw)) paste0("  HCp=", signif(hc_raw, 3)) else "",
    if (is.finite(htau) && !is.na(htau)) paste0("  Htau=", signif(htau, 3)) else "",
    if (n_above_axis > 0L) paste0("  above y-limit=", n_above_axis) else ""
  )

  plot_title <- if (is.null(short_title)) {
    compact_title(pretty_dataset_label(dataset_name), width = 42)
  } else short_title

  p <- ggplot() +
    maybe_rasterise(
      geom_point(
        data = plot_df,
        aes(x = lfc_shrunk, y = neglog10_empirical_p),
        color = plot_palette$background,
        shape = 16,
        size = 0.30,
        alpha = 0.35
      )
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

  # Original is independent of the EVS preprocessing path and is analyzed once,
  # under the first/original NormEVS registry entry only. Every added EVS track
  # contributes Lead/Rem views without redundantly re-running Original.
  analysis_dataset_list <- dataset_list
  if (!identical(track_key, "normalized_evs")) analysis_dataset_list$raw_dataset <- NULL

  analysis_results <- list()
  for (nm in names(analysis_dataset_list)) {
    dataset_name <- paste(comparison_name, unname(track_short[track_key]), nm, sep = "_")
    fit <- run_core_analysis(
      count_mat = analysis_dataset_list[[nm]],
      coldata = coldata,
      dataset_name = dataset_name,
      annot_df = annot_df
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
  tracks <- active_evs_track_keys()

  rows <- list(
    data.frame(
      track_key = "normalized_evs",
      dataset_key = "raw_dataset",
      stringsAsFactors = FALSE
    )
  )

  for (track_key in tracks) {
    rows[[length(rows) + 1L]] <- data.frame(
      track_key = track_key,
      dataset_key = "leading_edge_dataset",
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- data.frame(
      track_key = track_key,
      dataset_key = "remainder_dataset",
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(rows)
}

analysis_view_levels <- function() {
  views <- comparison_analysis_views()
  unique(vapply(seq_len(nrow(views)), function(i) {
    analysis_view_label(views$track_key[i], views$dataset_key[i])
  }, character(1)))
}

analysis_view_short_map <- function() {
  lev <- analysis_view_levels()
  out <- setNames(lev, lev)
  out["Original (No EVS)"] <- "Orig"

  for (nm in names(track_short)) {
    lead <- paste0(unname(track_short[nm]), " Lead")
    rem  <- paste0(unname(track_short[nm]), " Rem")
    prefix <- switch(
      nm,
      normalized_evs = "N",
      raw_evs = "R",
      cpm_evs = "C",
      vst_evs = "V",
      unname(track_short[nm])
    )
    if (lead %in% names(out)) out[lead] <- paste0(prefix, "-Lead")
    if (rem %in% names(out)) out[rem] <- paste0(prefix, "-Rem")
  }

  out
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
      EVS_k_per_condition = if (identical(dataset_key, "raw_dataset")) NA_integer_ else get_active_evs_cutoff(comparison_name),
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
        EVS_k_per_condition = if (identical(dataset_key, "raw_dataset")) NA_integer_ else get_active_evs_cutoff(comparison_name),
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
    # comparison-level panels generated after all active analysis views are fit.
    if (isTRUE(EXPORT_INDIVIDUAL_VIEW_FIGURES)) {
      dataset_name <- paste(comparison_name, unname(track_short[track_key]), dataset_key, sep = "_")
      volcano <- plot_final_volcano(
        df,
        dataset_name = dataset_name,
        short_title = paste0(pretty_comparison(comparison_name), " | ", analysis_label),
        label_genes = TRUE
      )
      # Single-column final print size.
      save_figure(
        volcano,
        file.path(paths$figures, "Volcano"),
        width = FIG_SINGLE_COL_MM,
        height = 82
      )
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
        EVS_k_per_condition = if (identical(dataset_key, "raw_dataset")) NA_integer_ else get_active_evs_cutoff(comparison_name),
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
    "## PAS filtering and pairwise comparisons",
    paste0("The analytical unit was the polyadenylation-site (PAS) feature identified by OrigID in the WTTS-Seq raw-count matrix. The analyzed pairwise comparisons were ", paste(comparison_table$comparison_name, collapse = ", "), ". Within each comparison, PASs with zero counts across all samples in the two compared groups were removed; all remaining PASs were retained for analysis."),
    "",
    "## Eigenvector splitting",
    paste0("This analysis used the ", CUTOFF_METHOD_LABEL, " cutoff rule. The comparison-specific top-k values were ", paste(paste0(comparison_table$comparison_name, " k=", vapply(comparison_table$comparison_name, get_active_evs_cutoff, integer(1))), collapse = ", "), ". NormEVS ranked PASs using PCA of DESeq2 median-of-ratios normalized counts, and RawEVS ranked PASs using PCA of raw counts. PCA was performed separately within each condition by applying prcomp to the transposed PAS-by-sample matrix with centering and without feature scaling. PASs were ranked by absolute PC1 loading, and the k highest-ranking PASs were selected independently in the RT and ZT arms."),
    "PASs selected in both arms were classified as Joint. PASs selected in only one arm were classified as Disjoint for that arm. The union of Joint and Disjoint PASs defined the Leading Edge; all unselected PASs defined the Remainder. EVS determined PAS membership only. Differential-expression testing for every view used the corresponding subset of raw integer counts.",
    paste0(if (identical(CUTOFF_METHOD_KEY, "cpm_empirical")) "For the CPM empirical analysis, an additional CPMEVS track ranked PASs using experiment-wide log1p(CPM) values. " else "", if (identical(CUTOFF_METHOD_KEY, "vst_empirical")) "For the VST empirical analysis, an additional VSTEVS track ranked PASs using experiment-wide DESeq2 variance-stabilized values with blind=TRUE. " else "", "The Original view contained all nonzero PASs without EVS splitting."),
    "",
    "## Differential-expression analysis",
    "Each analysis view was fit independently with DESeq2 using a negative-binomial generalized linear model with design ~ condition and the ZT/untrt group as the reference. DESeq2 estimated median-of-ratios size factors, dispersions, and Wald statistics. Log2 fold changes for the RT/trt coefficient were shrunken with apeglm and the shrunken estimate was used as the reported effect size.",
    sprintf("Standard significance required the ordinary two-sided DESeq2 Wald Benjamini-Hochberg adjusted p-value < %.2f and |apeglm-shrunken log2 fold change| >= %.1f.", BH_FDR_STANDARD, lfc_boundary),
    sprintf("Strong composite-null significance used DESeq2 results with lfcThreshold=%.1f and altHypothesis='greaterAbs', followed by Benjamini-Hochberg adjustment at FDR %.2f; reported Strong PASs also required |apeglm-shrunken log2 fold change| >= %.1f.", lfc_boundary, BH_FDR_STRONG, lfc_boundary),
    sprintf("Weak composite-null support used DESeq2 results with lfcThreshold=%.1f and altHypothesis='lessAbs', followed by Benjamini-Hochberg adjustment at FDR %.2f and |apeglm-shrunken log2 fold change| < %.1f. Final Weak significance additionally required HBFSS significance.", lfc_boundary, BH_FDR_WEAK, lfc_boundary),
    "",
    "## Empirical-null calibration, higher criticism, and HBFSS",
    sprintf("Finite ordinary DESeq2 Wald statistics were supplied to fdrtool with statistic='normal' and cutoff.method='fndr' to estimate empirical-null p-values. Higher-Criticism scores were calculated from ordered empirical-null p-values with fdrtool::hc.score. The HC search was restricted to the lowest %.0f%% of ordered empirical p-values (alpha0=%.2f). When the maximum HC score was positive, HCp was the empirical p-value at that maximum; otherwise HCp and Htau were undefined and no PAS in that view was classified as HBFSS-significant.", 100 * HC_ALPHA0, HC_ALPHA0),
    sprintf("For each PAS, HBFSS was |apeglm-shrunken log2 fold change| multiplied by -log10(empirical-null p). Htau was -log10(HCp) multiplied by %.1f. A PAS was HBFSS-significant when HBFSS > Htau.", lfc_boundary),
    "",
    "## Method overlap",
    "Method overlap was defined as the set of HBFSS-significant PASs that also satisfied at least one DESeq2-based criterion: Standard, Strong, or Weak-CNH. Final Weak significance was the intersection of Weak-CNH and HBFSS.",
    "",
    "## PCA summaries after EVS",
    "PCA summaries used the same preprocessing matrix used to define each EVS track. Original, Leading Edge, and Remainder subsets were evaluated on that common scale. Reported metrics included PC1 and PC2 variance explained, PC1 eigenvalue retention relative to the Original matrix, the fraction of Original-PC1 loading energy contained in each subset, and RT-ZT separation in the PC1-PC2 plane defined as centroid distance divided by pooled within-group root-mean-square distance.",
    "",
    "## 3'aTWAS overlap",
    "Human 3'aTWAS gene symbols were mapped to rat genes using babelgene human-to-rat ortholog mappings together with direct case-insensitive symbol-equivalent matches present in the WTTS annotation. For each comparison and analysis view, mapped TWAS orthologs were intersected with Standard, Strong, final Weak, and HBFSS WTTS discoveries."
  )

  writeLines(methods_text, file.path(output_dir, "Methods_Manuscript.md"), useBytes = TRUE)
  invisible(TRUE)
}

run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  dir.create(file.path(output_dir, comparison_name), recursive = TRUE, showWarnings = FALSE)

  for (track_key in active_evs_track_keys()) {
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
  short_view <- analysis_view_short_map()

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
      panel_title = paste0(pretty_comparison(comparison_name), " | ", CUTOFF_METHOD_LABEL, " | significance across ", length(plots), " analysis view", ifelse(length(plots) == 1L, "", "s")),
      ncol = length(plots)
    )

    panel_dir <- file.path(output_dir, comparison_name, "Panels")
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
    base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_Volcano_", length(plots), "Views"))

    # Size the canvas from the grid the compositor will use: panels wrap to at
    # most PANEL_MAX_COLS columns inside a double-column width.
    grid_dim <- panel_grid_dim(length(plots), length(plots))
    panel_width_mm <- min(
      FIG_DOUBLE_COL_MM,
      max(FIG_SINGLE_COL_MM, 60 * unname(grid_dim[["ncol"]]))
    )
    panel_height_mm <- min(
      FIG_MAX_HEIGHT_MM,
      18 + 72 * unname(grid_dim[["nrow"]])
    )
    save_figure(panel, base, width = panel_width_mm, height = panel_height_mm)
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
    " | lambda1=", round(retained, 1), "% Orig",
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
    for (track_key in active_evs_track_keys()) {
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

  for (track_key in active_evs_track_keys()) {
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
        labels = c("PC1 explained", "lambda1 vs Orig", "Orig PC1 energy")
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
      title = paste0(pretty_comparison(comparison_name), " | quantitative EVS/PC1 evidence"),
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
        if (is.finite(htau) && !is.na(htau)) paste0("Htau=", signif(htau, 3)) else "Htau=NA"
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

    for (track_key in active_evs_track_keys()) {
      obj <- paper_registry[[registry_key(comparison_name, track_key)]]
      if (is.null(obj) || is.null(obj$evs)) next

      original_mat <- pca_track_matrix(obj, "raw_dataset")
      original_stats <- pca_variance_stats(original_mat, obj$coldata)
      if (is.null(original_stats)) next

      track_abbr <- switch(
        track_key,
        normalized_evs = "Norm",
        raw_evs = "Raw",
        cpm_evs = "CPM",
        vst_evs = "VST",
        unname(track_short[track_key])
      )

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
        panel_title = paste0(pretty_comparison(comparison_name), " | ", CUTOFF_METHOD_LABEL, " | PCA structure before and after EVS"),
        ncol = 3
      )
      base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_PCA_EVS_2x3"))
      save_figure(panel, base, width = FIG_DOUBLE_COL_MM, height = 165)
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
      save_figure(pevidence, base, width = FIG_DOUBLE_COL_MM, height = 92)
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
      track_order_use <- active_evs_track_keys()
      panel_title <- paste0("HBFSS calibration | ", unname(dataset_short[dataset_key]), " | active EVS tracks | ", CUTOFF_METHOD_LABEL)
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

    # panel_height is computed upstream in inches; rescale it to the printable
    # page while preserving the intended aspect ratio.
    hbfss_height_mm <- min(
      FIG_MAX_HEIGHT_MM,
      max(70, (panel_height / 18.0) * FIG_DOUBLE_COL_MM)
    )
    save_figure(
      panel,
      file.path(
        paper_fig_dir,
        paste0("Figure_Manuscript_Empirical_HBFSS_", output_suffix)
      ),
      width = FIG_DOUBLE_COL_MM,
      height = hbfss_height_mm
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

  save_figure(
    p,
    file.path(paper_fig_dir, "Figure_Manuscript_Discovery_Counts"),
    width = FIG_DOUBLE_COL_MM,
    height = 115
  )
  invisible(TRUE)
}
# =============================================================================
# DISPERSION TRADE-OFF FIGURE
# =============================================================================
# DESeq2 dispersion fits are summarized for Leading Edge and Remainder strata.
# The sweep evaluates residual dispersion spread across k values and marks the
# active comparison-specific cutoff. No DE p-values are used in this module.
# =============================================================================

dispersion_frame_from_dds <- function(dds) {
  if (is.null(dds)) return(NULL)

  md <- tryCatch(SummarizedExperiment::mcols(dds), error = function(e) NULL)
  if (is.null(md)) return(NULL)

  needed <- c("baseMean", "dispGeneEst", "dispFit")
  if (!all(needed %in% colnames(md))) return(NULL)

  out <- data.frame(
    feature_id  = as.character(rownames(dds)),
    baseMean    = as.numeric(md$baseMean),
    dispGeneEst = as.numeric(md$dispGeneEst),
    dispFit     = as.numeric(md$dispFit),
    stringsAsFactors = FALSE
  )

  keep <- is.finite(out$baseMean) & out$baseMean > 0 &
    is.finite(out$dispGeneEst) & out$dispGeneEst > 0 &
    is.finite(out$dispFit) & out$dispFit > 0
  out <- out[keep, , drop = FALSE]
  if (!nrow(out)) return(NULL)

  out$log_resid <- log(out$dispGeneEst / out$dispFit)
  out
}

dispersion_prior_var <- function(dds) {
  f <- tryCatch(DESeq2::dispersionFunction(dds), error = function(e) NULL)
  if (is.null(f)) return(NA_real_)
  v <- attr(f, "dispPriorVar")
  if (is.null(v) || !is.finite(v)) NA_real_ else as.numeric(v)
}

# Residual spread around the fitted trend. MAD is used rather than SD so a
# handful of dispersion outliers cannot drive the comparison between cutoffs.
dispersion_residual_mad <- function(df) {
  if (is.null(df) || !nrow(df)) return(NA_real_)
  stats::mad(df$log_resid, constant = 1.4826, na.rm = TRUE)
}

collect_dispersion_fits <- function() {
  rows <- list()
  views <- comparison_analysis_views()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (i in seq_len(nrow(views))) {
      track_key   <- views$track_key[i]
      dataset_key <- views$dataset_key[i]

      obj <- paper_registry[[registry_key(comparison_name, track_key)]]
      if (is.null(obj) || is.null(obj$analysis_results)) next
      entry <- obj$analysis_results[[dataset_key]]
      if (is.null(entry) || is.null(entry$dds)) next

      df <- dispersion_frame_from_dds(entry$dds)
      if (is.null(df)) next

      df$comparison <- pretty_comparison(comparison_name)
      df$track      <- unname(track_short[track_key])
      df$view       <- analysis_view_label(track_key, dataset_key)
      df$stratum    <- switch(
        dataset_key,
        raw_dataset          = "Original (unsplit)",
        leading_edge_dataset = "Leading edge",
        remainder_dataset    = "Remainder",
        dataset_key
      )
      df$prior_var  <- dispersion_prior_var(entry$dds)
      rows[[length(rows) + 1L]] <- df
    }
  }

  if (!length(rows)) return(NULL)
  out <- dplyr::bind_rows(rows)
  out$stratum <- factor(
    out$stratum,
    levels = c("Original (unsplit)", "Leading edge", "Remainder")
  )
  out
}

DISPERSION_STRATUM_COLORS <- c(
  "Original (unsplit)" = "#4A5560",
  "Leading edge"       = "#009E73",
  "Remainder"          = "#E69F00"
)

# -----------------------------------------------------------------------------
# Panels
# -----------------------------------------------------------------------------

plot_dispersion_trends <- function(disp_df, track_label) {
  d <- disp_df[disp_df$track == track_label |
                 disp_df$stratum == "Original (unsplit)", , drop = FALSE]
  if (!nrow(d)) return(NULL)

  # Scatter is thinned per facet: a full leading edge plus remainder is tens of
  # thousands of points per comparison and would dominate the vector file.
  pts <- dplyr::bind_rows(lapply(
    split(d, list(d$comparison, d$stratum), drop = TRUE),
    function(g) if (nrow(g) > 2500L) g[sample.int(nrow(g), 2500L), , drop = FALSE] else g
  ))

  trend <- d[order(d$comparison, d$stratum, d$baseMean), , drop = FALSE]

  p <- ggplot() +
    geom_point(
      data = pts,
      aes(x = baseMean, y = dispGeneEst),
      colour = "grey78", size = 0.25, stroke = 0, alpha = 0.55
    ) +
    geom_line(
      data = trend,
      aes(x = baseMean, y = dispFit, colour = stratum, linetype = stratum),
      linewidth = 0.45
    ) +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_manual(values = DISPERSION_STRATUM_COLORS, drop = FALSE) +
    scale_linetype_manual(
      values = c("Original (unsplit)" = "22", "Leading edge" = "solid",
                 "Remainder" = "solid"),
      drop = FALSE
    ) +
    facet_wrap(~ comparison, nrow = 2) +
    labs(
      title = "Fitted dispersion trends",
      x = "Mean of normalized counts",
      y = "Dispersion"
    ) +
    manuscript_theme() +
    theme(legend.position = "bottom", legend.title = element_blank())
  p
}

# Ratio of the two fitted trends where both strata actually have features. If
# that overlap is narrow the comparison is not supported and the panel says so
# rather than extrapolating either trend beyond its data.
plot_dispersion_trend_ratio <- function(disp_df, track_label, n_grid = 120L) {
  d <- disp_df[disp_df$track == track_label &
                 disp_df$stratum %in% c("Leading edge", "Remainder"), , drop = FALSE]
  if (!nrow(d)) return(NULL)

  rows <- list()
  for (cmp in unique(d$comparison)) {
    le  <- d[d$comparison == cmp & d$stratum == "Leading edge", , drop = FALSE]
    rem <- d[d$comparison == cmp & d$stratum == "Remainder", , drop = FALSE]
    if (nrow(le) < 50L || nrow(rem) < 50L) next

    lo <- max(stats::quantile(le$baseMean, 0.01), stats::quantile(rem$baseMean, 0.01))
    hi <- min(stats::quantile(le$baseMean, 0.99), stats::quantile(rem$baseMean, 0.99))
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo * 1.2) next

    grid <- exp(seq(log(lo), log(hi), length.out = n_grid))
    f_le  <- stats::approx(le$baseMean,  le$dispFit,  xout = grid, rule = 2, ties = mean)$y
    f_rem <- stats::approx(rem$baseMean, rem$dispFit, xout = grid, rule = 2, ties = mean)$y

    rows[[length(rows) + 1L]] <- data.frame(
      comparison = cmp,
      baseMean   = grid,
      ratio      = f_le / f_rem,
      n_le       = sum(le$baseMean  >= lo & le$baseMean  <= hi),
      n_rem      = sum(rem$baseMean >= lo & rem$baseMean <= hi),
      stringsAsFactors = FALSE
    )
  }

  if (!length(rows)) return(NULL)
  r <- dplyr::bind_rows(rows)
  r <- r[is.finite(r$ratio) & r$ratio > 0, , drop = FALSE]
  if (!nrow(r)) return(NULL)

  lab <- r %>%
    dplyr::group_by(comparison) %>%
    dplyr::summarise(
      baseMean = min(baseMean, na.rm = TRUE),
      ratio    = min(ratio, na.rm = TRUE),
      txt      = paste0(format(dplyr::first(n_le), big.mark = ","), " LE / ",
                        format(dplyr::first(n_rem), big.mark = ","), " Rem"),
      .groups  = "drop"
    )

  ggplot(r, aes(x = baseMean, y = ratio)) +
    geom_hline(yintercept = 1, colour = "grey35", linetype = "22", linewidth = 0.35) +
    geom_line(colour = "#0072B2", linewidth = 0.5) +
    geom_text(
      data = lab, aes(label = txt),
      hjust = 0, vjust = 0, size = 1.9, colour = "grey35"
    ) +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(~ comparison, nrow = 2) +
    labs(
      title = "Leading-edge vs remainder trend, matched mean",
      x = "Mean of normalized counts (overlapping range only)",
      y = "Dispersion ratio, LE / Rem"
    ) +
    manuscript_theme()
}

plot_dispersion_residual_spread <- function(disp_df, track_label) {
  d <- disp_df[disp_df$track == track_label |
                 disp_df$stratum == "Original (unsplit)", , drop = FALSE]
  if (!nrow(d)) return(NULL)

  s <- d %>%
    dplyr::group_by(comparison, stratum) %>%
    dplyr::summarise(
      mad = stats::mad(log_resid, constant = 1.4826, na.rm = TRUE),
      .groups = "drop"
    )
  s <- s[is.finite(s$mad), , drop = FALSE]
  if (!nrow(s)) return(NULL)

  ggplot(s, aes(x = comparison, y = mad, fill = stratum)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.66) +
    geom_text(
      aes(label = sprintf("%.3f", mad)),
      position = position_dodge(width = 0.75),
      vjust = -0.4, size = 1.9
    ) +
    scale_fill_manual(values = DISPERSION_STRATUM_COLORS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "Residual spread around the fitted trend",
      x = NULL,
      y = "MAD of log(gene-wise / fitted) dispersion"
    ) +
    manuscript_theme() +
    theme(legend.position = "bottom", legend.title = element_blank())
}

# -----------------------------------------------------------------------------
# Cutoff sweep
# -----------------------------------------------------------------------------
# Dispersions only: no model fitting, no testing, no shrinkage. The split is
# re-derived from the stored PC1 loading tables, so no PCA is recomputed either.

fit_dispersions_only <- function(count_mat, coldata) {
  if (is.null(count_mat) || !nrow(count_mat)) return(NULL)

  dds <- tryCatch(
    DESeq2::DESeqDataSetFromMatrix(
      countData = coerce_raw_count_matrix_for_deseq2(
        count_mat,
        context = "dispersion sweep subset"
      ),
      colData = coldata,
      design  = make_design_formula(coldata)
    ),
    error = function(e) NULL
  )
  if (is.null(dds)) return(NULL)

  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
  if (nrow(dds) < 50L) return(NULL)

  dds <- tryCatch(
    DESeq2::estimateSizeFactors(dds),
    error = function(e) tryCatch(
      DESeq2::estimateSizeFactors(dds, type = "poscounts"),
      error = function(e2) NULL
    )
  )
  if (is.null(dds)) return(NULL)

  dds <- tryCatch(
    DESeq2::estimateDispersions(dds, quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(dds)) return(NULL)

  dispersion_frame_from_dds(dds)
}

dispersion_sweep_one <- function(comparison_name, track_key, k_grid) {
  obj <- paper_registry[[registry_key(comparison_name, track_key)]]
  if (is.null(obj) || is.null(obj$evs)) return(NULL)

  evs <- obj$evs
  raw <- evs$raw_dataset
  cd  <- obj$coldata
  lt_t <- evs$fit_trt$loading_table
  lt_u <- evs$fit_untrt$loading_table
  if (is.null(lt_t) || is.null(lt_u)) return(NULL)

  all_ids <- rownames(raw)
  n_total <- length(all_ids)
  # Include the active cutoff in the evaluated dispersion sweep.
  k_active <- as.integer(get_active_evs_cutoff(comparison_name))
  k_grid <- c(as.integer(k_grid), k_active)
  k_grid <- sort(unique(k_grid[k_grid >= 500L & k_grid < n_total / 2]))
  if (!length(k_grid)) return(NULL)

  rows <- list()
  for (k in k_grid) {
    top_t <- as.character(lt_t$feature_id[lt_t$rank <= k])
    top_u <- as.character(lt_u$feature_id[lt_u$rank <= k])
    le_ids  <- union(top_t, top_u)
    rem_ids <- setdiff(all_ids, le_ids)
    if (!length(le_ids) || length(rem_ids) < 50L) next

    d_le  <- fit_dispersions_only(raw[le_ids, , drop = FALSE], cd)
    d_rem <- fit_dispersions_only(raw[rem_ids, , drop = FALSE], cd)
    if (is.null(d_le) || is.null(d_rem)) next

    m_le  <- dispersion_residual_mad(d_le)
    m_rem <- dispersion_residual_mad(d_rem)
    if (!is.finite(m_le) || !is.finite(m_rem)) next

    n_le <- nrow(d_le); n_rem <- nrow(d_rem)
    rows[[length(rows) + 1L]] <- data.frame(
      comparison       = pretty_comparison(comparison_name),
      track            = unname(track_short[track_key]),
      k                = k,
      mad_leading_edge = m_le,
      mad_remainder    = m_rem,
      n_leading_edge   = n_le,
      n_remainder      = n_rem,
      weighted_mad     = (m_le * n_le + m_rem * n_rem) / (n_le + n_rem),
      stringsAsFactors = FALSE
    )
  }

  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

run_dispersion_sweep <- function(track_key, k_grid) {
  rows <- list()
  for (comparison_name in as.character(comparison_table$comparison_name)) {
    message("  dispersion sweep: ", comparison_name, " ",
            unname(track_short[track_key]), " (", length(k_grid), " values of k)")
    r <- tryCatch(
      dispersion_sweep_one(comparison_name, track_key, k_grid),
      error = function(e) {
        warning("Dispersion sweep failed for ", comparison_name, ": ",
                conditionMessage(e))
        NULL
      }
    )
    if (!is.null(r)) rows[[length(rows) + 1L]] <- r
  }
  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

plot_dispersion_sweep <- function(sweep_df, unsplit_mad = NULL) {
  if (is.null(sweep_df) || !nrow(sweep_df)) return(NULL)

  active <- data.frame(
    comparison = pretty_comparison(as.character(comparison_table$comparison_name)),
    k = vapply(
      as.character(comparison_table$comparison_name),
      function(cmp) as.integer(get_active_evs_cutoff(cmp)),
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
  active <- merge(active, sweep_df[, c("comparison", "k", "weighted_mad")],
                  by = c("comparison", "k"), all.x = FALSE)

  minima <- sweep_df %>%
    dplyr::group_by(comparison) %>%
    dplyr::slice_min(weighted_mad, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  p <- ggplot(sweep_df, aes(x = k, y = weighted_mad, colour = comparison)) +
    geom_line(linewidth = 0.5) +
    geom_point(data = minima, shape = 21, fill = "white", size = 1.3, stroke = 0.5) +
    scale_x_continuous(labels = function(v) format(v, big.mark = ",", trim = TRUE)) +
    labs(
      title = paste0("Cutoff sweep | ", CUTOFF_METHOD_LABEL),
      x = "Leading-edge size, k (PAS per condition)",
      y = "Size-weighted residual spread",
      caption = paste0(
        "Open circles mark the minimum; diamonds mark the cutoff used in this run. ",
        "Computed from dispersion fits only, before any testing."
      )
    ) +
    manuscript_theme() +
    theme(legend.position = "bottom", legend.title = element_blank())

  if (nrow(active)) {
    p <- p + geom_point(
      data = active,
      aes(x = k, y = weighted_mad, colour = comparison),
      shape = 23, fill = "white", size = 1.7, stroke = 0.6,
      inherit.aes = FALSE, show.legend = FALSE
    )
  }

  if (!is.null(unsplit_mad) && is.finite(unsplit_mad)) {
    p <- p + geom_hline(
      yintercept = unsplit_mad,
      colour = "grey35", linetype = "22", linewidth = 0.35
    )
  }

  p
}

# -----------------------------------------------------------------------------
# Assembly
# -----------------------------------------------------------------------------

assemble_dispersion_panel <- function(plot_list, panel_title) {
  plot_list <- Filter(Negate(is.null), plot_list)
  if (!length(plot_list)) return(NULL)

  tagged <- lapply(seq_along(plot_list), function(i) {
    plot_list[[i]] +
      labs(tag = LETTERS[i]) +
      theme(
        plot.tag = element_text(face = "bold", size = base_theme_size + 1,
                                family = FIG_FONT_FAMILY),
        plot.tag.position = "topleft"
      )
  })

  ncol_use <- if (length(tagged) > 1L) 2L else 1L

  body <- if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::patchworkGrob(patchwork::wrap_plots(tagged, ncol = ncol_use))
  } else {
    do.call(gridExtra::arrangeGrob, c(tagged, list(ncol = ncol_use)))
  }

  gridExtra::arrangeGrob(
    body,
    ncol = 1,
    top = grid::textGrob(
      panel_title,
      gp = grid::gpar(
        fontface = "bold",
        fontsize = base_theme_size + 2,
        fontfamily = if (nzchar(FIG_FONT_FAMILY)) FIG_FONT_FAMILY else ""
      )
    )
  )
}

save_dispersion_tradeoff_panels <- function() {
  if (!isTRUE(EXPORT_DISPERSION_TRADEOFF)) return(invisible(FALSE))

  disp_df <- collect_dispersion_fits()
  if (is.null(disp_df) || !nrow(disp_df)) {
    warning("Dispersion trade-off figure skipped: no dispersion data collected.")
    return(invisible(FALSE))
  }

  dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)

  unsplit <- disp_df[disp_df$stratum == "Original (unsplit)", , drop = FALSE]
  unsplit_mad <- if (nrow(unsplit)) dispersion_residual_mad(unsplit) else NA_real_

  # One figure per active EVS track: the tracks differ in how the ranking was
  # built, so their dispersion geometry is not interchangeable.
  for (track_key in active_evs_track_keys()) {
    track_label <- unname(track_short[track_key])

    sweep_df <- if (isTRUE(DISPERSION_SWEEP_ENABLED) &&
                    identical(track_key, DISPERSION_SWEEP_TRACK)) {
      run_dispersion_sweep(track_key, DISPERSION_SWEEP_GRID)
    } else {
      NULL
    }

    if (!is.null(sweep_df) && nrow(sweep_df)) {
      save_csv(
        sweep_df,
        file.path(
          summary_table_dir,
          paste0("Table_Dispersion_Sweep_", track_label, ".csv")
        )
      )
    }

    panels <- list(
      plot_dispersion_trends(disp_df, track_label),
      plot_dispersion_trend_ratio(disp_df, track_label),
      plot_dispersion_residual_spread(disp_df, track_label),
      plot_dispersion_sweep(sweep_df, unsplit_mad)
    )

    g <- assemble_dispersion_panel(
      panels,
      paste0("Mean-variance basis for the EVS split | ", track_label,
             " | ", CUTOFF_METHOD_LABEL)
    )
    if (is.null(g)) next

    save_figure(
      g,
      file.path(paper_fig_dir, paste0("Figure_Dispersion_Tradeoff_", track_label)),
      width = FIG_DOUBLE_COL_MM,
      height = 200
    )
  }

  # Per-view residual spread and prior variance, as a table the Methods can cite.
  summary_tbl <- disp_df %>%
    dplyr::group_by(comparison, track, view, stratum) %>%
    dplyr::summarise(
      n_features        = dplyr::n(),
      median_baseMean   = stats::median(baseMean),
      residual_mad      = stats::mad(log_resid, constant = 1.4826, na.rm = TRUE),
      dispersion_prior_var = dplyr::first(prior_var),
      .groups = "drop"
    )

  save_csv(
    summary_tbl,
    file.path(summary_table_dir, "Table_Dispersion_Residual_Spread.csv")
  )

  invisible(TRUE)
}



save_paper_support_figures <- function(summary_df) {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) {
    return(invisible(FALSE))
  }

  support_steps <- list(
    PCA = function() save_paper_pca_panels(),
    Empirical_HBFSS = function() save_paper_empirical_hbfss_panels(),
    Counts = function() save_discovery_count_panel(summary_df),
    Dispersion = function() save_dispersion_tradeoff_panels()
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
      Analysis = factor(Analysis, levels = analysis_view_levels())
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
      factor(Analysis, levels = analysis_view_levels()),
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
      factor(Analysis, levels = analysis_view_levels()),
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

  analysis_map <- analysis_view_short_map()
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
  # Height grows with the number of genes but is capped at the printable page.
  height_mm <- min(FIG_MAX_HEIGHT_MM, max(90, 55 + 2.6 * max_genes))
  save_figure(
    p,
    file.path(figure_dir, "Figure_TWAS_Gene_Support"),
    width = FIG_DOUBLE_COL_MM,
    height = height_mm
  )
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
          paste0(pretty_comparison(comparison_name), " | ", analysis_label, " | 3'aTWAS")
        )
        if (!is.null(p)) {
          save_figure(
            p,
            file.path(paths$figures, "TWAS_Overlap"),
            width = FIG_SINGLE_COL_MM,
            height = 70
          )
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

    preferred_views <- analysis_view_levels()
    views <- preferred_views[preferred_views %in% unique(as.character(sub$Analysis))]
    plots <- lapply(views, function(v) {
      one <- sub[sub$Analysis == v, , drop = FALSE]
      if (!nrow(one)) return(NULL)
      short <- unname(analysis_view_short_map()[[v]])
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
    twas_dim <- panel_grid_dim(n_views, n_views)
    save_figure(
      panel,
      base,
      width = min(FIG_DOUBLE_COL_MM,
                  max(FIG_SINGLE_COL_MM, 60 * unname(twas_dim[["ncol"]]))),
      height = min(FIG_MAX_HEIGHT_MM, 16 + 62 * unname(twas_dim[["nrow"]]))
    )
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

empirical_cutoff_export <- comparison_table[, c(
  "comparison_name", "group1_prefix", "group2_prefix",
  "fixed_5000_k", "cpm_empirical_k", "vst_empirical_k"
), drop = FALSE]
names(empirical_cutoff_export) <- c(
  "Comparison", "RT_prefix", "ZT_prefix",
  "Fixed_5000_k", "CPM_Empirical_k", "VST_Empirical_k"
)
empirical_cutoff_export$Active_Cutoff_Method <- CUTOFF_METHOD_LABEL
empirical_cutoff_export$Active_k_per_condition <- vapply(
  empirical_cutoff_export$Comparison, get_active_evs_cutoff, integer(1)
)
empirical_cutoff_export$Basis <- EMPIRICAL_EVS_CUTOFF_BASIS
empirical_cutoff_export$Applied_to <- "NormEVS and RawEVS"
save_csv(
  empirical_cutoff_export,
  file.path(summary_table_dir, "Table_EVS_Cutoffs.csv")
)
# Per-PAS leading-edge membership, consumed by the master's cross-cutoff
# sensitivity stage to identify PASs that every cutoff admitted and tested.
build_leading_edge_membership_table <- function() {
  rows <- list()
  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (track_key in active_evs_track_keys()) {
      obj <- paper_registry[[registry_key(comparison_name, track_key)]]
      if (is.null(obj) || is.null(obj$evs)) next
      ids <- obj$evs$leading_edge_ids
      if (!length(ids)) next
      rows[[length(rows) + 1L]] <- data.frame(
        Comparison = comparison_name,
        Track = unname(track_short[track_key]),
        EVS_k_per_condition = obj$evs$empirical_evs_k,
        PAS = as.character(ids),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(data.frame())
  dplyr::bind_rows(rows)
}

leading_edge_membership <- build_leading_edge_membership_table()
if (nrow(leading_edge_membership) > 0L) {
  save_csv(
    leading_edge_membership,
    file.path(summary_table_dir, "Table_EVS_Leading_Edge_Membership.csv")
  )
}


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

# Final analytical stage: 3'aTWAS ortholog overlap with every reported WTTS method.
twas_results <- run_twas_overlap_analysis()

comparison_zips <- setNames(
  vapply(as.character(comparison_table$comparison_name), create_comparison_artifact_package, character(1)),
  as.character(comparison_table$comparison_name)
)
figure_zip <- create_all_figures_zip()
table_zip <- create_all_tables_zip()
cutoff_zip <- create_cutoff_artifact_package()

cat("\n=====================================================\n")
cat("Pipeline complete.\n")
cat("Build: ", PIPELINE_BUILD, "\n", sep = "")
cat("Cutoff method: ", CUTOFF_METHOD_LABEL, "\n", sep = "")
cat("Repository root:\n", repo_root, "\n", sep = "")
cat("Output directory:\n", output_dir, "\n", sep = "")
cat("3'aTWAS overlap: complete\n")
cat("Figure ZIP:\n", figure_zip, "\n", sep = "")
cat("Table ZIP:\n", table_zip, "\n", sep = "")
cat("Cutoff ZIP:\n", cutoff_zip, "\n", sep = "")
cat("Per-comparison ZIPs:\n", paste(unname(comparison_zips), collapse = "\n"), "\n", sep = "")
cat("=====================================================\n\n")

if (nrow(overall_summary) > 0L) print(overall_summary)
