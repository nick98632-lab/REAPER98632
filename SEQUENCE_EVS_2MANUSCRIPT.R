#!/usr/bin/env Rscript
# =============================================================================
# SEQUENCE EVS manuscript pipeline — hierarchical dispersion model (V2)
# =============================================================================
# Run: Rscript SEQUENCE_EVS_FINAL.R
# Set SEQUENCE_REPO_ROOT to the directory containing the original WTTS count
# matrix and 3'aTWAS input. SEQUENCE_MULTI_ROOT controls the output directory.
# A previous run in the output directory is moved aside with a timestamp.
# Primary endpoint: ordinary two-sided Wald/BH per comparison and analysis.
# Secondary: strong and nonzero-small effect tests; HBFSS is exploratory.
# EVS is a novel empirical-Bayes extension, not a validated FDR guarantee.
# Cutoff selection uses a discretized Cox–Reid adjusted plug-in likelihood
# and a log-normal dispersion prior. The cutoff is selected without DE tests,
# but arm-specific ranks and eligibility use tested counts; an independent
# reference ranking is run as a sensitivity analysis. Environment options:
# SEQUENCE_CORES (parallel comparisons, default up to 4), SEQUENCE_RANDOM_SPLITS
# (permutation searches, default 5), SEQUENCE_DISP_GRID_STEP (0.25),
# SEQUENCE_PRIOR_VAR_MIN (0.01), SEQUENCE_NULL_METHOD (fndr).
# Calibration: SEQUENCE_EVS_CALIBRATION.R. Methods text: MANUSCRIPT_METHODS.
# =============================================================================

MANUSCRIPT_METHODS <- c(
  "# SEQUENCE EVS | Statistical methods",
  "",
  "**Hierarchical dispersion model with eigenvector splitting**  ",
  "Methods accompanying `SEQUENCE_EVS_FINAL.R` (build SEQUENCE_EVS_MANUSCRIPT_HIERARCHICAL_V2) · 25 November 2026",
  "",
  "This document describes the algorithm as implemented. Numerical results are written by the pipeline to the Summary_Tables folder and are not part of this text. The novel components are the eigenvector-guided partition and its marginal adjusted-likelihood cutoff; the cited literature supports the constituent statistical methods, not EVS itself.",
  "",
  "## Study design and input",
  "",
  "The analytical unit is a WTTS-Seq polyadenylation site (PAS), identified by OrigID; gene symbols are annotations, never units of differential testing. The study is a 4 × 2 factorial time course: rats in recovery after sleep deprivation (R arms, treatment) sampled at 0, 2, 4 and 8 h of recovery, and undisturbed controls sampled at the same clock times (ZT6, ZT8, ZT10 and ZT14), five samples per arm (40 samples, eight arms). Four prespecified pairwise comparisons are analysed: RT0 vs ZT6, RT2 vs ZT8, RT4 vs ZT10 and RT8 vs ZT14; the full design is also analysed jointly (Time-course model). For each comparison a PAS is eligible if it has at least one count in at least one of the comparison's samples (MIN_PAS_COUNT = 1, MIN_PAS_SAMPLES = 1, prespecified); the same universe is used for the Original and both partitions. No result-dependent filtering is performed. Rows with missing or nonnumeric counts are removed at import; negative counts stop the run. Samples are assumed independent; no batch, pairing or longitudinal structure is modelled. [1, 5]",
  "",
  "## Normalization and eigenvector ranking",
  "",
  "DESeq2 median-of-ratios size factors are estimated once per comparison on the full eligible comparison matrix (poscounts if the standard estimator fails). The same factors and the same gene-wise dispersion estimates are supplied to the Original and to both partitions, so the views differ only in their dispersion moderation. Raw counts, never VST values, enter the negative-binomial models. As a check for compositional bias, the size factors are compared with edgeR TMM and upper-quartile effective library sizes (each rescaled to geometric mean 1); agreement within 20% for every sample is recorded as a validation check, with per-sample library sizes and zero fractions. [1, 11]",
  "",
  "For ranking only, a blind variance-stabilizing transformation (VST) is computed across all 40 samples. Within each arm, PCA (centred, unscaled) is applied to the PAS-by-sample VST matrix and PASs are ranked by |PC1 loading|. For PAS i, the union depth dᵢ is the smaller of its two positions in the comparison's arms (1 = largest loading). At cutoff k, the Leading Edge (LE) is every eligible PAS with dᵢ ≤ k (the union of the two arms' top-k sets) and the Remainder (REM) is the complement; the partition is exhaustive and disjoint, and ties follow a fixed feature order. The ranking uses the tested arms (without condition labels or test statistics) and the VST uses all samples, so the split is data-adaptive; an independent reference ranking built from held-out samples only is analysed as a sensitivity analysis (below). [1, 5]",
  "",
  "## Hierarchical dispersion likelihood",
  "",
  "For PAS i in sample j the negative-binomial model has mean μᵢⱼ = sⱼ exp(xⱼᵀβᵢ) and variance μᵢⱼ + αᵢμᵢⱼ². The log dispersion zᵢ = log αᵢ has a population-specific normal prior with centre mᵢ = log(a₀ + a₁ μ̄ᵢ⁻ᵇ) and variance τ²; a₀ > 0, a₁ > 0 and b are estimated. The trend extends DESeq2's parametric trend (b = 1). [1–3]",
  "",
  "DESeq2 gene-wise dispersion estimation supplies fitted means μ̂ᵢⱼ. Conditional on these, each PAS's Cox–Reid adjusted log likelihood is ℓᵢ(z) = Σⱼ log NB(yᵢⱼ; μ̂ᵢⱼ, eᶻ) − ½ log det(XᵀWᵢX), with Wᵢ = diag(μ̂ᵢⱼ/(1 + eᶻμ̂ᵢⱼ)); for the saturated categorical design the determinant is the product of the within-group weight sums. This is a plug-in adjusted likelihood (means held fixed), not a fully profiled or integrated GLM likelihood. [1, 2]",
  "",
  "ℓᵢ is evaluated on a fixed grid of log dispersions from log(10⁻⁶) to log(10³) with spacing 0.25. Below 10⁻⁶ the likelihood is close to its Poisson limit at the count levels in these data, so assigning all prior mass below the grid to the lowest cell has negligible effect; above 10³ the likelihood is negligible. Normal CDF differences assign the prior mass to grid cells, and the marginal objective is Σᵢ log Σₕ Lᵢ(zₕ) P(zᵢ ∈ cell h | mᵢ, τ²), with Lᵢ the relative likelihood. PASs with weakly informative likelihoods therefore contribute less, boundary gene-wise estimates are retained, and no uniform digamma bias correction is needed.",
  "",
  "The parameters are estimated by bound-constrained L-BFGS-B with analytic gradients (at most 1,000 iterations), with a₀ ∈ [10⁻⁸, 10³], a₁ ∈ [10⁻⁸, 10⁸], b ∈ [0.5, 3.5] and τ² ∈ [0.01, 16]. The lower bound 0.01 on τ² regularizes against collapse of the prior variance to zero and is identical in the cutoff search and in the final moderation; a population whose estimate reaches the bound is flagged. The estimated τ² is the prior variance used in DESeq2's MAP shrinkage. DESeq2's conventional 0.25 floor on the prior variance is run as a reference analysis. Pooled and final fits use three starting values; nonconverged fits are excluded. The bounds, the grid spacing and the 0.01 regularization bound are prespecified implementation choices. Their numerical adequacy is checked for every comparison at k*: the Leading-Edge and Remainder fits are repeated from several starting values, with the grid spacing halved, and with the τ² bound reduced tenfold, and the changes in log likelihood, trend parameters and τ² are reported and checked automatically. [1–3]",
  "",
  "## Cutoff selection and model comparison",
  "",
  "For each comparison the other six arms (30 samples) supply the held-out likelihood curves (design ~ arm; size factors estimated on those 30 samples alone). For each admissible k, separate hierarchical marginal likelihoods are fitted to the LE and REM PASs, and the profile log likelihood is the sum of the two maximized objectives. Admissible cutoffs leave at least 200 PASs in each partition, in both the held-out fitting set and the comparison universe, with k at most half the comparison's PASs. The search evaluates the admissible endpoints and a grid of step 250, then refines at steps 25, 5 and 1 around the three best candidates of each stage (ties prefer smaller k). Each fit starts from the parameters of the nearest evaluated cutoff and falls back to the pooled start, with multiple starts, if that fit does not converge. The selected k* is the largest converged objective on the evaluated grid: a grid-refined maximum, not a proven global optimum. No Wald statistic, p-value, discovery count or HBFSS score enters selection, and a counter of significance computations verifies this. The four comparisons are independent and are searched in parallel with fixed per-comparison seeds.",
  "",
  "A BIC-type criterion compares one pooled trend (4 parameters) with the split (8 parameters plus the cutoff), with the number of held-out PASs in the penalty. The split is used for the primary results only if it improves this criterion; otherwise the Original is the selected model and EVS is reported as exploratory. Because the cutoff is not identified under a single trend, the likelihood-ratio statistic does not have its usual distribution, and the criteria are comparative diagnostics, not p-values.",
  "",
  "Four comparators are fitted on the same held-out likelihood: (i) a continuous rank effect, a centred natural cubic spline (3 df) of log(1 + dᵢ) added to the prior centre, testing whether dispersion changes gradually rather than in two populations; (ii) a binary split ranked by within-arm VST variance; (iii) a binary split ranked by each PAS's contribution to the PC1–PC2 subspace, λ₁v²ᵢ₁ + λ₂v²ᵢ₂ (using all nonzero components would reduce this to the variance); and (iv) permutation searches, in which the whole ranking is permuted and searched on the coarse grid at the same resolution as the observed coarse-grid maximum (five permutations by default; SEQUENCE_RANDOM_SPLITS). Comparators (ii) and (iii) use the same search with stages 250 and 25. The model with the lowest criterion is reported for each comparison; the comparators are descriptive and do not replace the prespecified pooled-versus-EVS selection. [1–3, 5]",
  "",
  "## Differential PAS abundance",
  "",
  "For each comparison the Original and the two partitions use the same raw counts, shared size factors, condition design and gene-wise dispersions. The hierarchical model is refitted to the PASs of each view; its trend and τ² are supplied to DESeq2 (dispersionFunction, estimateDispersionsMAP), and DESeq2's dispersion-outlier rule and Wald test are retained. The Original refitted with DESeq2's default parametric trend is reported as a reference. [1, 2]",
  "",
  "One Cook's-distance exclusion list is derived from the Original (DESeq2's default rule) and applied to every view; excluded PASs enter neither the tests nor the empirical null, higher criticism or HBFSS. Benjamini–Hochberg adjustment is applied per comparison and test, across the Original's PASs or across the LE and REM PASs combined, without independent filtering. The primary endpoint is the ordinary two-sided Wald test of β = 0 at BH < 0.10. The secondary strong-effect test is DESeq2's greaterAbs test of |β| > 1 at BH < 0.10. The secondary weak-effect test is an intersection–union test of 0 < |β| < 1, using the larger of the ordinary p-value and the lessAbs equivalence p-value, at BH < 0.20; a prespecified minimum effect δ (δ < |β| < 1, replacing the ordinary test by greaterAbs at δ) is available but set to 0. No fold-change gate is added to any test, and no FDR claim is made for the union of the test families. For comparison with earlier analyses, Standard calls with |apeglm log₂FC| ≥ 1 are also counted, as a descriptive column only (a fold-change filter applied after BH is not FDR-controlled; the Strong test is the formal counterpart). [1, 4, 5]",
  "",
  "Fold changes are shrunk with apeglm for interpretation and visualization, with the prior scale adapted within each view, consistent with treating the views as separate populations; shrunken values do not define rejections. The selected-model tables report one analysis per comparison; all views are retained in the all-PAS and secondary/exploratory tables. [1, 6]",
  "",
  "## Exploratory higher criticism and HBFSS",
  "",
  "For finite Wald statistics of non-excluded PASs, fdrtool fits a normal empirical null within each view (cutoff method fndr). Higher criticism is maximized over the lowest 10% of the Original's ordered empirical p-values and must be positive; its threshold p_HC is applied to every view of the comparison. HBFSS = |apeglm log₂ fold change| × [−log₁₀(empirical p)], and a score above −log₁₀(p_HC) flags a PAS for prioritization. HBFSS calls are exploratory, are not part of the primary discoveries, and have no claimed FDR control. [7, 8]",
  "",
  "## Sensitivity analyses",
  "",
  "For each comparison, under the primary testing rules (the Original's Cook's list and HBFSS threshold, one BH per analysis): (i) the Original and EVS at k* with DESeq2's 0.25 floor on the MAP prior variance; (ii) EVS with an independent reference ranking, in which the VST is computed on the 30 held-out samples alone, PASs are positioned by |PC1 loading| within each held-out arm, a condition's reference position is the mean over its three held-out arms, the union depth is the smaller of the treatment and control positions, and the cutoff is selected by the same search, so that no tested count enters the split; the overlap of its Leading Edge with the primary one is reported (Jaccard index); (iii) independent hypothesis weighting (IHW) of the Original's p-values with the PC1 union depth as covariate, the established covariate-weighting alternative to splitting, run when the IHW package is installed; and (iv) negative binomial count splitting, in which each tested count is split into two independent negative binomial parts by beta-binomial thinning (ε = 0.5, size 1/α from the Original's final dispersion), the split is ranked on one half (VST and per-arm PC1 of the ten training samples, primary k*), and the Original and EVS are both tested on the other half, so the ranking is independent of the tested counts. [9, 10, 15]",
  "",
  "## Time-course model",
  "",
  "The factorial design is analysed jointly in one DESeq2 model across all 40 samples, ~ time + condition + time:condition, with experiment-wide size factors and the hierarchical dispersion trend (one population). Likelihood-ratio tests against ~ time (a treatment effect at any recovery time) and against ~ time + condition (a treatment effect that changes with recovery time) are adjusted by BH within each test. The treatment log₂ fold change and Wald test at each recovery time are obtained from the condition and interaction coefficients. This tests time explicitly instead of comparing separate significance lists. [1, 16]",
  "",
  "## Calibration",
  "",
  "Error rates are assessed by simulation (SEQUENCE_EVS_CALIBRATION.R), using the pipeline's own functions. For each comparison, 10,000 real PASs are subsampled; the truth's means, size factors and dispersion model (trend and τ², from the same marginal likelihood) are fitted on the real counts, either as one population or as the Leading Edge and Remainder at k* scaled to the subsample, optionally with the real samples' PC1 structure. Ten percent of PASs receive true log₂ fold changes; null scenarios have none. The Original and EVS are run with the estimated prior variance and with the 0.25 floor. Each call type is scored against its own truth: Standard and HBFSS against any true change, Strong against |true log₂FC| > 1, and Weak against 0 < |true log₂FC| < 1. FDR, false non-discovery rate, power, false calls under the null, the empirical-null spread, and the error of the final dispersions and shrunken fold changes are reported.",
  "",
  "## Alternative polyadenylation",
  "",
  "Differential abundance does not establish a change in which PAS a gene uses, so APA is analysed separately for each pairwise comparison with edgeR quasi-likelihood models (robust dispersion) on raw counts, with offsets from the shared DESeq2 size factors. Gene units are defined by gene symbol, chromosome and strand, and PASs are ordered along the direction of transcription by their peak coordinate. (i) Usage change: for genes with at least two eligible PASs, diffSpliceDGE tests each PAS against the other PASs of its gene; gene-level p-values are Simes combinations, with BH across genes, and, when the stageR package is installed, a stage-wise analysis screens genes by the Simes p-value and confirms individual PASs within significant genes with overall FDR control at 0.10. (ii) Direction: for each gene, the two most used PASs (mean normalized count ≥ 5 each) form a proximal/distal pair, and the relative expression difference RED is the treatment log₂ fold change of the distal PAS minus that of the proximal PAS, i.e. the change in log₂ distal/proximal usage; positive values indicate 3′ UTR lengthening (more distal use in recovery) and negative values shortening. The gene-level diffSpliceDGE test of the pair gives the p-value, with BH across genes (FDR < 0.10). Whether the pair's PASs lie in the EVS Leading Edge is reported. [11, 12, 17]",
  "",
  "3′aTWAS genes are mapped from human to rat with database-supported babelgene orthologs (at least one supporting source); symbol similarity alone is not orthology. Gene-level enrichment uses the eligible, finite-p-value, ortholog-mapped rat genes of the selected model as the background, each gene counted once. Logistic regression relates any primary PAS discovery in a gene to TWAS membership, adjusting for log(1 + number of tested PASs), log(1 + mean normalized count) and log(number of supporting ortholog databases); a sensitivity analysis restricts to one-to-one mappings. Odds ratios, Wald 95% intervals and p-values are reported with BH adjustment across the enrichment analyses. The same model relates significant APA genes (shortening or lengthening) to 3′aTWAS membership, with the ortholog-mapped tested genes of the APA pair analysis as background, adjusting for log(1 + mean normalized count of the pair) and ortholog support; because 3′aTWAS identifies genetically regulated APA, this is the most direct cross-species comparison. Overlap is corroboration, not causal evidence. [13, 14]",
  "",
  "## Figures",
  "",
  "Figure 1, the EVS cutoff (held-out dispersions with the Leading-Edge, Remainder and single trends at k*; profile likelihood over k). Figure 2, DESeq2 dispersions after the split (dispersion plots, trend exponent, prior variance used). Figure 3, model comparison and sensitivity (BIC of each alternative relative to EVS; observed versus permuted rankings; Standard discoveries in the primary and sensitivity analyses). Figure 4, the time-course model (significant PASs per recovery time; treatment log₂FC over recovery for the strongest time-varying PASs). Figure 5, alternative polyadenylation (genes with a usage change and with significant shortening or lengthening per recovery time; RED volcano plots; 3′aTWAS enrichment of APA genes). Per comparison: volcano plots, PCA before and after EVS, PC1 evidence and 3′aTWAS overlap; discovery counts and HBFSS calibration across comparisons; 3′aTWAS gene support; calibration (FDR, false non-discovery rate, power, null false calls, empirical-null spread, estimation error) and the HBFSS null choice.",
  "",
  "## Reproducibility and limitations",
  "",
  "The pipeline records the cutoff profile, fitted coefficients and prior variances, the selected model, all comparator criteria, the sensitivity analyses, full PAS-level results, input MD5 checksums, sessionInfo(), settings and seeds, and a table of automated validation checks. The time-course model uses one dispersion population; EVS is applied to the pairwise comparisons. Count splitting halves the counts available to each step, so its discovery counts are compared only with the Original on the same half. The nominal 10% and 20% FDR levels are prespecified, not optimized. No discovery count was used to tune any setting.",
  "",
  "# References",
  "",
  "1. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15:550. https://doi.org/10.1186/s13059-014-0550-8",
  "2. McCarthy DJ, Chen Y, Smyth GK. Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research. 2012;40:4288–4297. https://doi.org/10.1093/nar/gks042",
  "3. Wu H, Wang C, Wu Z. A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data. Biostatistics. 2013;14:232–243. https://doi.org/10.1093/biostatistics/kxs033",
  "4. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society B. 1995;57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
  "5. Bourgon R, Gentleman R, Huber W. Independent filtering increases detection power for high-throughput experiments. PNAS. 2010;107:9546–9551. https://doi.org/10.1073/pnas.0914005107",
  "6. Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. 2019;35:2084–2092. https://doi.org/10.1093/bioinformatics/bty895",
  "7. Strimmer K. fdrtool: a versatile R package for estimating local and tail area-based false discovery rates. Bioinformatics. 2008;24:1461–1462. https://doi.org/10.1093/bioinformatics/btn209",
  "8. Donoho D, Jin J. Higher criticism thresholding: optimal feature selection when useful features are rare and weak. PNAS. 2008;105:14790–14795. https://doi.org/10.1073/pnas.0807471105",
  "9. Ignatiadis N, Klaus B, Zaugg JB, Huber W. Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. Nature Methods. 2016;13:577–580. https://doi.org/10.1038/nmeth.3885",
  "10. Efron B. Simultaneous inference: when should hypothesis testing problems be combined? Annals of Applied Statistics. 2008;2:197–223. https://doi.org/10.1214/07-AOAS141",
  "11. Chen Y, Lun ATL, Smyth GK. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. F1000Research. 2016;5:1438. https://doi.org/10.12688/f1000research.8987.2 (edgeR version as recorded in SessionInfo.txt.)",
  "12. Van den Berge K, Soneson C, Robinson MD, Clement L. stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology. 2017;18:151. https://doi.org/10.1186/s13059-017-1277-0",
  "13. Alternative polyadenylation transcriptome-wide association study identifies APA-linked susceptibility genes in brain disorders. Nature Communications. 2023;14. https://doi.org/10.1038/s41467-023-36311-8",
  "14. babelgene: gene orthologs for model organisms in a tidy data format. https://igordot.github.io/babelgene/",
  "15. Neufeld A, Popp J, Gao LL, Battle A, Witten D. Negative binomial count splitting for single-cell RNA sequencing data. arXiv:2307.12985 (2023). https://arxiv.org/abs/2307.12985",
  "16. Love MI, Anders S, Huber W. Analyzing RNA-seq data with DESeq2 (Bioconductor vignette; sections on likelihood ratio tests and interactions). https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html",
  "17. Hoque M, Ji Z, Zheng D, Luo W, Li W, You B, Park JY, Yehia G, Tian B. Analysis of alternative cleavage and polyadenylation by 3′ region extraction and deep sequencing. Nature Methods. 2013;10:133–139. https://doi.org/10.1038/nmeth.2288"
)

PIPELINE_BUILD <- "SEQUENCE_EVS_MANUSCRIPT_HIERARCHICAL_V2"


# =============================================================================
# CONSTANTS AND SHARED FIGURE STYLE
# =============================================================================

COMPARISON_NAMES <- c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14")
METHOD_SLUGS <- c("EVS k*" = "EVS_k")
CUTOFF_METHOD_LABEL <- names(METHOD_SLUGS)
CUTOFF_METHOD_BASIS <- paste0("Experiment-wide VST; grid-refined marginal adjusted-likelihood split of ",
 "held-out counts with log-normal dispersion priors; tested-arm ranks are data-adaptive; ",
 "no DE statistics used; pooled model eligible by information criterion")


# EVS PCA input: the experiment-wide DESeq2 VST (see EIGENVECTOR SPLITTING).
PCA_METHOD_LEVELS <- "VST"
PCA_METHOD_LABELS <- c(VST = "experiment-wide DESeq2 VST (blind = TRUE)")
EVS_PCA_METHOD <- "VST"
ANALYSIS_FOLDER <- "EVS_Analysis"
ORIGINAL_LABEL <- "Original (No EVS)"
SET_COLS <- c("Original (No EVS)" = "#2B2B2B", "Leading Edge" = "#1F6FB2", "Remainder" = "#C0392B")

# Purpose: Comp lab.
comp_lab <- function(cmp = COMPARISON_NAMES) stats::setNames(sub("_", " / ", cmp), cmp)

# Purpose: Fig theme.
fig_theme <- function(base = 7.5) {
  theme_minimal(base_size = base, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
      axis.line = element_line(linewidth = 0.3, colour = "grey35"),
      axis.ticks = element_line(linewidth = 0.3, colour = "grey35"),
      axis.text = element_text(size = base - 1.3, colour = "grey20"),
      axis.title = element_text(size = base - 0.4, colour = "grey10"),
      strip.background = element_rect(fill = "grey94", colour = NA),
      strip.text = element_text(face = "bold", size = base - 0.3, colour = "grey10",
                                margin = margin(2.5, 2, 2.5, 2)),
      panel.spacing = unit(4, "mm"),
      legend.title = element_blank(), legend.text = element_text(size = base - 1),
      legend.key.size = unit(3, "mm"), legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
      plot.title = element_text(size = base + 0.6, face = "bold", margin = margin(b = 1.5)),
      plot.subtitle = element_text(size = base - 1, colour = "grey35", margin = margin(b = 4)),
      plot.title.position = "plot", plot.tag = element_text(size = base + 4, face = "bold"),
      plot.margin = margin(4, 5, 4, 4))
}

# PDF (vector) and 600-dpi PNG of a ggplot, patchwork or grob.
save_figure <- function(p, out_dir, stem, w_mm, h_mm) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  w <- w_mm / 25.4; h <- h_mm / 25.4
  draw <- function() if (inherits(p, "ggplot")) print(p) else { grid::grid.newpage(); grid::grid.draw(p) }
  cairo <- isTRUE(capabilities("cairo"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  if (cairo) grDevices::cairo_pdf(pdf_path, width = w, height = h) else
    grDevices::pdf(pdf_path, width = w, height = h, useDingbats = FALSE)
  tryCatch(draw(), finally = grDevices::dev.off())
  grDevices::png(file.path(out_dir, paste0(stem, ".png")), width = w, height = h, units = "in", res = 600,
                 bg = "white", type = if (cairo) "cairo" else getOption("bitmapType"))
  tryCatch(draw(), finally = grDevices::dev.off())
  invisible(TRUE)
}

# rows: list; each element is a ggplot or list(plots = list(...), widths = c(...)).
assemble_rows <- function(rows, heights, title) {
  title_theme <- theme(plot.title = element_text(size = 10, face = "bold", margin = margin(b = 6)))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    built <- lapply(rows, function(r) if (inherits(r, "ggplot")) r else
      patchwork::wrap_plots(r$plots, nrow = 1, widths = r$widths))
    return(patchwork::wrap_plots(built, ncol = 1, heights = heights) +
             patchwork::plot_annotation(title = title, theme = title_theme))
  }
  built <- lapply(rows, function(r) if (inherits(r, "ggplot")) ggplot2::ggplotGrob(r) else
    gridExtra::arrangeGrob(grobs = lapply(r$plots, ggplot2::ggplotGrob), nrow = 1, widths = r$widths))
  gridExtra::arrangeGrob(grobs = built, ncol = 1, heights = heights,
                         top = grid::textGrob(title, x = 0.01, hjust = 0,
                                              gp = grid::gpar(fontsize = 10, fontface = "bold")))
}

# Purpose: Try fig.
try_fig <- function(expr, what) {
  tryCatch(expr, error = function(e) {
    warning("Figure skipped (", what, "): ", conditionMessage(e), call. = FALSE)
    invisible(NULL)
  })
}

# =============================================================================
# MASTER PROCESS: run the analysis once (single-threaded), then one ZIP
# =============================================================================

SEQUENCE_CHILD_RUN <- identical(Sys.getenv("SEQUENCE_CHILD_RUN", unset = "0"), "1")

if (!SEQUENCE_CHILD_RUN) {
  script_args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", script_args, value = TRUE)
  if (!length(script_arg)) stop("Run this pipeline with Rscript so the master process can relaunch itself.")
  this_script <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)

  count_name <- "WTTS-Seq_2022.2_DE_raw_read_numbers.csv"
  root_candidates <- unique(c(Sys.getenv("SEQUENCE_REPO_ROOT", unset = ""), getwd(), dirname(getwd()),
                              "/root/REAPER98632", dirname(this_script)))
  root_candidates <- root_candidates[nzchar(root_candidates)]
  repo_guess <- NULL
  for (cand in root_candidates) {
    if (file.exists(file.path(cand, count_name)) || file.exists(file.path(cand, "data", count_name))) {
      repo_guess <- normalizePath(cand, winslash = "/", mustWork = TRUE)
      break
    }
  }
  if (is.null(repo_guess)) repo_guess <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

  multi_root <- Sys.getenv("SEQUENCE_MULTI_ROOT",
                           unset = file.path(repo_guess, "exports", "sequence_evs_hierarchical_candidate"))
  # A previous run is never overwritten: it is moved aside with a timestamp.
  if (dir.exists(multi_root) && length(list.files(multi_root, all.files=TRUE, no..=TRUE))) {
    moved <- paste0(multi_root, "_previous_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    if (!file.rename(multi_root, moved)) stop("Output directory is not empty and could not be moved: ", multi_root)
    message("Previous run moved to ", moved)
  }
  dir.create(multi_root, recursive = TRUE, showWarnings = FALSE)

  count_path <- if (file.exists(file.path(repo_guess, count_name))) file.path(repo_guess, count_name) else
    file.path(repo_guess, "data", count_name)
  if (!file.exists(count_path)) stop("Cannot find ", count_name, ". Set SEQUENCE_REPO_ROOT.")

  # One BLAS / OpenMP thread (unless set by the user).
  thread_vars <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
  thread_env <- paste0(thread_vars[!nzchar(Sys.getenv(thread_vars))], "=1")
  message("=====================================================")
  message("Starting the analysis run")
  message("=====================================================")
  status <- system2(
    command = file.path(R.home("bin"), "Rscript"), args = shQuote(this_script),
    env = c("SEQUENCE_CHILD_RUN=1", paste0("SEQUENCE_MULTI_ROOT=", multi_root),
            paste0("SEQUENCE_REPO_ROOT=", repo_guess), thread_env))
  if (!identical(as.integer(status), 0L)) stop("The analysis run failed (exit status ", status, ")")

  final_zip <- file.path(dirname(multi_root), paste0(basename(multi_root), ".zip"))
  if (file.exists(final_zip)) unlink(final_zip, force = TRUE)
  oldwd <- getwd()
  setwd(dirname(multi_root))
  zip_rel <- list.files(multi_root, recursive = TRUE)
  zip_rel <- zip_rel[!grepl("SEQUENCE_ALL_(FIGURES|TABLES)\\.zip$", zip_rel)]
  utils::zip(zipfile = basename(final_zip), files = file.path(basename(multi_root), zip_rel), flags = "-q")
  setwd(oldwd)
  if (!file.exists(final_zip)) stop("Final ZIP creation failed: ", final_zip)

  rd <- function(f) {
    p <- file.path(multi_root, ANALYSIS_FOLDER, "Summary_Tables", f)
    if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) else NULL
  }
  cat("\n=====================================================\n")
  cat("Analysis complete.\n")
  act <- rd("Table_Active_Cutoff.csv"); cnt <- rd("Table_DE_Method_Counts.csv")
  if (!is.null(act) && !is.null(cnt)) {
    for (i in seq_len(nrow(act))) {
      cmp <- act$Comparison[i]
      tot <- function(an, col) sum(cnt[[col]][cnt$Comparison == cmp & cnt$Analysis %in% an], na.rm = TRUE)
      o_std <- tot(ORIGINAL_LABEL, "Std"); e_std <- tot(c("NormEVS Lead", "NormEVS Rem"), "Std")
      o_hb <- tot(ORIGINAL_LABEL, "HBFSS"); e_hb <- tot(c("NormEVS Lead", "NormEVS Rem"), "HBFSS")
      cat(sprintf("  %-9s k* = %5d  LE %5d  EVS Std %5d (%+5d)  HBFSS %5d (%+5d)\n", cmp, act$k[i],
                  act$Leading_Edge[i], e_std, e_std - o_std, e_hb, e_hb - o_hb))
    }
  }
  chk <- rd("Table_Validation_Checks.csv")
  if (!is.null(chk)) {
    cat(sprintf("  Validation checks: %d PASS, %d FAIL\n", sum(chk$Status == "PASS"), sum(chk$Status == "FAIL")))
    bad <- chk[chk$Status == "FAIL", , drop = FALSE]
    for (i in seq_len(nrow(bad))) cat("    FAIL:", bad$Comparison[i], "-", bad$Check[i], "-", bad$Detail[i], "\n")
  }
  cat("Final ZIP:\n", normalizePath(final_zip, winslash = "/", mustWork = TRUE), "\n", sep = "")
  cat("=====================================================\n\n")
  quit(save = "no", status = 0L)
}

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
  "S4Vectors",
  "edgeR",
  "babelgene",
  "splines",
  "parallel"
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
# Optional minimum meaningful effect for Weak (delta < |LFC| < 1). Default 0:
# Weak tests 0 < |LFC| < 1. With delta > 0 the nonzero test is replaced by
# DESeq2's greaterAbs test at lfcThreshold = delta. Prespecify before analysis.
WEAK_DELTA <- as.numeric(Sys.getenv("SEQUENCE_WEAK_DELTA", "0"))
stopifnot(is.finite(WEAK_DELTA), WEAK_DELTA >= 0, WEAK_DELTA < 1)
# Prespecified abundance eligibility: default 1 count in 1 sample preserves the
# original universe. Override BEFORE looking at results for a sensitivity run.
MIN_PAS_COUNT <- 1L
MIN_PAS_SAMPLES <- 1L
# Purpose: Eligible pas.
eligible_pas <- function(x) rowSums(x >= MIN_PAS_COUNT) >= MIN_PAS_SAMPLES

# Higher Criticism searches only the lower tail of the ordered empirical-null
# p-value distribution. alpha0 = 0.10 means that the maximum HC score is sought
# only among the lowest 10% of empirical p-values. A non-positive maximum HC
# score yields no HC threshold and therefore no HBFSS candidates in that
# analysis. HC is computed once per analysis, over the empirical p-values of all
# its views (call_significance).
HC_ALPHA0 <- 0.10

if (!is.numeric(HC_ALPHA0) || length(HC_ALPHA0) != 1L ||
    !is.finite(HC_ALPHA0) || HC_ALPHA0 <= 0 || HC_ALPHA0 > 1) {
  stop("HC_ALPHA0 must be a single finite number in (0, 1].")
}

# EVS selection size k* is comparison-specific, chosen before any test
# (EVS CUTOFF k*).
EMPIRICAL_EVS_CUTOFF_BASIS <- paste0(
  CUTOFF_METHOD_BASIS,
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
# Statistical decision rules (separate families, per comparison and analysis).
# Standard: ordinary Wald BH < 0.10 (primary). Strong: greaterAbs BH < 0.10.
# Weak: BH < 0.20 for max(ordinary p, lessAbs p), testing 0 < |LFC| < 1.
# HBFSS: exploratory |apeglm LFC| x -log10(empirical p) > Htau, no FDR claim.
# All calls exclude common Cook's outliers. No post-hoc LFC rejection gates.
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

track_short <- c(normalized_evs = "NormEVS")

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

if (!nzchar(Sys.getenv("SEQUENCE_MULTI_ROOT", unset = ""))) {
  stop("Run this script through the master process (Rscript SEQUENCE_EVS_FINAL.R); do not launch the analysis run directly.")
}

comparison_table <- data.frame(
  comparison_name = c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14"),
  group1_prefix   = c("R0", "R2", "R4", "R8"),
  group2_prefix   = c("ZT6", "ZT8", "ZT10", "ZT14"),
  stringsAsFactors = FALSE
)

# k* per comparison, set by compute_cutoffs() before any comparison runs.
OPTIMAL_K <- list()
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(as.integer(Sys.getenv("OPENBLAS_NUM_THREADS", unset = "1")))
  RhpcBLASctl::omp_set_num_threads(as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1")))
}
# Run-wide validation checks, and the significance-call counter used to show
# that the cutoff rule performs no differential-expression test.
# The Original's Cook's outlier PASs per comparison: applied in every view.
original_cooks_cache <- list()
# The Original refitted with DESeq2's default parametric trend (reference only).
default_trend_rows <- list()
validation_rows <- list()
pca_ranking_rows <- list()   # per comparison: PC1 ranks and union depth
# Per comparison: the raw matrix, the arm ranks and union depth, shared size
# factors and gene-wise dispersions (prepare_comparison_cache), reused by the
# analysis.
comparison_cache <- list()
# wald: DESeq2::nbinomWaldTest fits; null: call_significance (fdrtool / HC / HBFSS)
SIGNIFICANCE_CALLS <- list2env(list(wald = 0L, null = 0L))

# Purpose: Get active evs cutoff.
get_active_evs_cutoff <- function(comparison_name) {
  idx <- match(as.character(comparison_name), comparison_table$comparison_name)
  if (is.na(idx)) stop("No EVS cutoff is defined for comparison: ", comparison_name)
  k <- OPTIMAL_K[[as.character(comparison_name)]] %||%
    stop("The cutoff k has not been chosen for ", comparison_name)
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || !is.finite(k) || k < 1L) {
    stop("Invalid EVS cutoff for ", comparison_name, " under ", CUTOFF_METHOD_LABEL)
  }
  k
}


if (anyDuplicated(comparison_table$comparison_name)) stop("comparison_table contains duplicated comparison names.")

# -----------------------------------------------------------------------------
# Repository and file helpers
# -----------------------------------------------------------------------------

# Purpose: Resolve existing file.
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

# Purpose: Get script path.
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

# Purpose: Find repo root.
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
  unset = file.path(repo_root, "exports", "sequence_evs_hierarchical_candidate")
)
output_dir <- file.path(multi_output_root, ANALYSIS_FOLDER)

if (dir.exists(output_dir) && length(list.files(output_dir, all.files=TRUE, no..=TRUE)))
  stop("Analysis directory is not empty: ",output_dir)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paper_fig_dir <- file.path(output_dir, "Combined_Figures")
summary_table_dir <- file.path(output_dir, "Summary_Tables")
dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_table_dir, recursive = TRUE, showWarnings = FALSE)

paper_registry <- list()

# Purpose: Registry key.
registry_key <- function(comparison_name, track_key) {
  paste(comparison_name, track_key, sep = "__")
}

# Purpose: Should write figure.
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

# Purpose: Assert required columns.
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

# Purpose: Safe neglog10.
safe_neglog10 <- function(x, pseudocount = 1e-12) {
  -log10(pmax(x, pseudocount))
}

# Purpose: Clip probabilities.
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

# Purpose: Finite plot df.
finite_plot_df <- function(df, x_col, y_col) {
  keep <- is.finite(df[[x_col]]) &
    !is.na(df[[x_col]]) &
    is.finite(df[[y_col]]) &
    !is.na(df[[y_col]])

  df[keep, , drop = FALSE]
}

# Purpose: Compact title.
compact_title <- function(x, width = 54) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

# Purpose: Compact caption.
compact_caption <- function(x, width = 118) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

# Purpose: Save csv.
save_csv <- function(df, path) {
  write.csv(df, file = path, row.names = FALSE)
}

# Purpose: Save grob.
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

# Purpose: Pretty dataset type.
pretty_dataset_type <- function(dataset_key) {
  switch(
    dataset_key,
    raw_dataset = "Original dataset",
    leading_edge_dataset = "Leading-edge dataset",
    remainder_dataset = "Remainder dataset",
    dataset_key
  )
}

# Purpose: Pretty dataset label.
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


# Purpose: Make design formula.
make_design_formula <- function(coldata) {
  ~ condition
}

# Purpose: Get condition coef.
get_condition_coef <- function(dds) {
  rn <- DESeq2::resultsNames(dds)
  idx <- grep("^condition_", rn)

  if (length(idx) == 0L) {
    stop("Could not identify condition coefficient in resultsNames(dds).")
  }

  rn[idx[1]]
}



# Purpose: Resolve top n cutoff.
resolve_top_n_cutoff <- function(sorted_values_desc, top_n) {
  n_total <- length(sorted_values_desc)
  top_n <- as.integer(top_n)

  if (n_total == 0L) {
    stop("resolve_top_n_cutoff() received an empty vector.")
  }

  if (!is.finite(top_n) || is.na(top_n) || top_n < 1L) {
    stop("top_n must be a positive integer.")
  }

  # Use the comparison-specific k exactly. Do not
  # silently shrink k or replace the rank rule with a loading-value cutoff.
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

# Null for the HBFSS empirical p-values (SEQUENCE_NULL_METHOD):
#   fndr         fdrtool's censored empirical null, central region chosen by the
#                false-nondiscovery-rate rule (fdrtool default; the paper's choice)
#   pct0         fdrtool's censored empirical null fitted to the central
#                EMPIRICAL_NULL_PCT0 fraction of the statistics (narrower centre,
#                less affected by many small real effects)
#   theoretical  N(0, 1), no empirical correction (reference only)
EMPIRICAL_NULL_METHOD <- Sys.getenv("SEQUENCE_NULL_METHOD", unset = "fndr")
EMPIRICAL_NULL_PCT0 <- 0.5
if (!EMPIRICAL_NULL_METHOD %in% c("fndr", "pct0", "theoretical")) {
  stop("SEQUENCE_NULL_METHOD must be fndr, pct0 or theoretical.")
}

# Purpose: Run empirical null fdrtool.
run_empirical_null_fdrtool <- function(stat_vec, dataset_name) {
  stat_vec <- as.numeric(stat_vec)
  stat_vec <- stat_vec[is.finite(stat_vec) & !is.na(stat_vec)]
  stat_vec <- unname(stat_vec)

  if (length(stat_vec) < 5L) {
    stop(sprintf("[%s] Fewer than 5 finite Wald statistics were available for fdrtool.", dataset_name))
  }

  if (identical(EMPIRICAL_NULL_METHOD, "theoretical")) {
    return(list(pval = clip_probabilities(2 * stats::pnorm(-abs(stat_vec))),
                param = matrix(c(1, NA_real_), nrow = 1, dimnames = list(NULL, c("sd", "eta0")))))
  }
  fit <- fdrtool::fdrtool(
    stat_vec,
    statistic = "normal",
    plot = FALSE,
    verbose = FALSE,
    cutoff.method = EMPIRICAL_NULL_METHOD,
    pct0 = EMPIRICAL_NULL_PCT0
  )
  fit$pval <- clip_probabilities(fit$pval)
  fit
}

# Purpose: Safe hc thresh.
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
  # asserted for the analysis.
  if (!is.finite(best_hc) || is.na(best_hc) || best_hc <= 0) {
    message(
      sprintf(
        "[%s] No positive Higher-Criticism signal within the lowest %.1f%% of empirical p-values; HBFSS disabled for this analysis.",
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

# Purpose: Build significance plot long.
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

# Purpose: Plot expand xy.
plot_expand_xy <- function() {
  list(
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.10))),
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))
  )
}

# Purpose: Manuscript theme.
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

# Purpose: Shared panel legend.
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



# Purpose: Strip legend.
strip_legend <- function(p) {
  p + theme(legend.position = "none")
}

# Purpose: Assemble one legend panel.
assemble_one_legend_panel <- function(plot_list, panel_title, ncol = length(plot_list),
                                      legend = build_full_method_legend()) {
  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0L) {
    return(NULL)
  }

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

# Purpose: Create all figures zip.
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
  # per-comparison three-view TWAS panels already contain the same count summary.
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

# Purpose: Create all tables zip.
create_all_tables_zip <- function() {
  # The manuscript table archive contains the concise tables used to interpret
  # significance, EVS/PCA evidence, and 3'aTWAS overlap. Per-view working tables
  # remain in their comparison folders but are not duplicated in the archive.
  wanted <- character(0)

  add_if_exists <- function(path) {
    if (file.exists(path)) wanted <<- c(wanted, path)
  }

  add_if_exists(file.path(summary_table_dir, "Table_DE_Method_Counts.csv"))
  for (tab in c("Table_EVS_Cutoffs.csv", "Table_Empirical_Null_By_Analysis.csv",
                "Table_Dispersion_After_EVS.csv", "Table_Active_Cutoff.csv",
                "Table_EVS_Composition.csv", "Table_Cutoff_Summary.csv", "Table_Cutoff_Profile_Likelihood.csv",
                "Table_Cutoff_Random_Splits.csv", "Table_Sensitivity_Analyses.csv",
                "Table_Selected_Model_Primary_Discoveries.csv", "Table_TWAS_Enrichment.csv",
                "Table_APA_Summary.csv", "Table_APA_Shortening_Lengthening.csv", "Table_APA_Gene_Usage_Simes.csv",
                "Table_APA_TWAS_Enrichment.csv", "Table_TimeCourse_Summary.csv",
                "Table_TimeCourse_PAS.csv", "Table_Normalization_Check.csv",
                "Table_Validation_Checks.csv", "Table_Original_DESeq2_Default_Trend.csv")) {
    add_if_exists(file.path(summary_table_dir, tab))
  }
  add_if_exists(file.path(summary_table_dir, "Table_EVS_Split_Audit.csv"))
  add_if_exists(file.path(summary_table_dir, "Table_EVS_PCA_Evidence.csv"))
  add_if_exists(file.path(summary_table_dir, "Table_EVS_PCA_Ranking.csv"))

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

# Purpose: Coerce count column.
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

# Purpose: Valid gene symbol.
valid_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & nzchar(x) & x != "-" & grepl("[A-Za-z0-9]", x)
}

# Purpose: Gene key.
gene_key <- function(x) {
  x <- trimws(as.character(x))
  x[!valid_gene_symbol(x)] <- NA_character_
  toupper(x)
}

# Purpose: Collapse unique.
collapse_unique <- function(x, sep = "; ") {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  paste(sort(x), collapse = sep)
}

# Purpose: Safe min numeric.
safe_min_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else min(x)
}

# Purpose: Safe max numeric.
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

# Purpose: Prepare comparison data.
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
  count_sub <- count_sub[eligible_pas(as.matrix(count_sub)), , drop = FALSE]

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


# Purpose: Coerce raw count matrix for deseq2.
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

# Per comparison: the comparison's slice of the experiment-wide
# PCA input and the per-condition PC1 fits, cached (keyed "<comparison>|<method>").
PCA_RANK_CACHE <- list()

# The matrix the EVS PCA is run on: the experiment-wide VST matrix
# (EVS CUTOFF k*) restricted to the comparison's PASs and samples, the same
# values the cutoff rule ranked. Used only to rank PASs; no
# transformed matrix is ever passed to DESeq2. (size_factors is unused: the
# ranking uses experiment-wide normalization.)
make_rank_matrix_for_track <- function(count_matrix, coldata, track_key, method = EVS_PCA_METHOD,
                                       size_factors = NULL, cache_key = NULL) {
  track_key <- match.arg(track_key, "normalized_evs")
  method <- match.arg(method, PCA_METHOD_LEVELS)
  key <- if (is.null(cache_key)) NULL else paste(cache_key, method, sep = "|")
  hit <- if (is.null(key)) NULL else PCA_RANK_CACHE[[key]]
  if (!is.null(hit) && identical(hit$ids, rownames(count_matrix)) && identical(hit$samples, colnames(count_matrix))) {
    return(hit$rank_obj)
  }
  full <- experiment_rank_matrix(method)
  miss <- c(setdiff(rownames(count_matrix), rownames(full)), setdiff(colnames(count_matrix), colnames(full)))
  if (length(miss)) stop("The experiment-wide ", method, " matrix lacks comparison PASs or samples.")
  x <- full[rownames(count_matrix), colnames(count_matrix), drop = FALSE]
  rank_obj <- list(
    rank_matrix = as.data.frame(x),
    method = method,
    preprocessing_label = paste0(PCA_METHOD_LABELS[[method]], " prior to EVS; ",
                                 "PCA is performed directly on this matrix (no further transformation)")
  )
  if (!is.null(key)) PCA_RANK_CACHE[[key]] <<- list(ids = rownames(count_matrix), samples = colnames(count_matrix),
                                                   rank_obj = rank_obj, fits = NULL)
  rank_obj
}

# Per-condition PC1 fits with the full ranking (cached with the matrix).
pc1_condition_fits <- function(rank_obj, coldata, cache_key = NULL) {
  key <- if (is.null(cache_key)) NULL else paste(cache_key, rank_obj$method, sep = "|")
  if (!is.null(key) && !is.null(PCA_RANK_CACHE[[key]]$fits)) return(PCA_RANK_CACHE[[key]]$fits)
  ids <- colnames(rank_obj$rank_matrix)
  fits <- lapply(c(trt = "trt", untrt = "untrt"), function(cond)
    compute_pc1_loading_table(rank_obj$rank_matrix, ids[coldata$condition == cond], top_n = 1L,
                              preprocessing_label = rank_obj$preprocessing_label))
  if (!is.null(key) && !is.null(PCA_RANK_CACHE[[key]])) PCA_RANK_CACHE[[key]]$fits <<- fits
  fits
}

# A cached PC1 fit at cutoff k: identical to compute_pc1_loading_table(top_n = k).
pc1_fit_at_k <- function(fit, k) {
  lt <- fit$loading_table
  ci <- resolve_top_n_cutoff(lt$pc1_loading_abs, top_n = k)
  lt$split_class <- ifelse(lt$rank <= ci$top_n_actual, "high_loading", "background_loading")
  fit$loading_table <- lt
  fit$cutoff <- ci$cutoff_value
  fit$top_n_used <- ci$top_n_actual
  fit
}

# Purpose: Compute pc1 loading table.
compute_pc1_loading_table <- function(value_df, sample_names, top_n, preprocessing_label = "Normalized prior to EVS") {
  x <- as.matrix(value_df[, sample_names, drop = FALSE])
  storage.mode(x) <- "numeric"

  if (ncol(x) < 2L) {
    stop("EVS PCA requires at least two samples in each condition.")
  }

  # Refined EVS method: run PCA directly on the pre-split matrix for the
  # condition, extract the PC1 feature eigenvector/loading, take absolute values,
  # rank descending, and select the comparison-specific k. The matrix is used
  # as supplied (make_rank_matrix_for_track); nothing is log-transformed again.
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

# Purpose: Build eigenvector split.
build_eigenvector_split <- function(comparison_name, count_matrix, coldata, track_key) {
  track_key <- match.arg(track_key, "normalized_evs")
  empirical_k <- get_active_evs_cutoff(comparison_name)

  # The EVS matrix is used only to choose feature IDs. NormEVS
  # follows the supplied workflow: median-of-ratios normalization of the full
  # RT/ZT comparison matrix, followed by condition-specific PCA/PC1 eigenvectors.
  # Downstream DESeq2 always receives the corresponding raw-count subset.
  raw_count_matrix <- coerce_raw_count_matrix_for_deseq2(
    count_matrix,
    context = paste0(track_key, " full comparison matrix before EVS")
  )

  if (nrow(raw_count_matrix) < empirical_k) {
    stop(
      comparison_name, " ", track_key,
      ": EVS k=", empirical_k,
      " exceeds the ", nrow(raw_count_matrix),
      " PASs available after zero-row filtering."
    )
  }

  # The comparison slice of the experiment-wide PCA input and its PC1 rankings (cached).
  rank_obj <- make_rank_matrix_for_track(
    count_matrix = raw_count_matrix,
    coldata = coldata,
    track_key = track_key,
    size_factors = comparison_cache[[comparison_name]]$sf,
    cache_key = comparison_name
  )

  rank_matrix <- rank_obj$rank_matrix
  preprocessing_label <- rank_obj$preprocessing_label

  fits <- pc1_condition_fits(rank_obj, coldata, comparison_name)
  fit_trt <- pc1_fit_at_k(fits$trt, empirical_k)
  fit_untrt <- pc1_fit_at_k(fits$untrt, empirical_k)

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
      ": EVS failed to select exactly k=", empirical_k,
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
      ": remainder dataset is empty after EVS k=", empirical_k, "."
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
    pca_method = rank_obj$method,
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

# =============================================================================
# HIERARCHICAL DISPERSION MODEL (Love et al. 2014; McCarthy et al. 2012;
# Wu et al. 2013). The EVS marginal adjusted likelihood is a novel extension.
# =============================================================================
# Each PAS has a latent log dispersion z ~ N(log(a0+a1*mean^-b), tau^2).
# The negative-binomial likelihood is evaluated on a fixed grid at the fitted
# means from estimateDispersionsGeneEst, with the Cox-Reid adjustment.
# This is a plug-in adjusted likelihood, not a fully integrated GLM likelihood.
# Boundary gene-wise estimates are retained; no digamma correction is needed.
TREND_B_RANGE <- c(0.5, 3.5)
TREND_MIN_COEF <- 1e-8
TREND_MAX_COEF <- c(a0 = 1e3, a1 = 1e8)
# Latent prior variance tau^2 (natural-log dispersion scale). It is estimated
# by the marginal likelihood for every population, with a small lower bound
# (DISP_PRIOR_VAR_MIN) that regularizes against collapse to zero; this bound
# applies identically to the cutoff search and to the final moderation.
# DESeq2's conventional 0.25 floor on the MAP prior variance is run as a
# reference analysis (SENSITIVITY ANALYSES), not imposed on the primary fit.
DESEQ2_PRIOR_VAR_FLOOR <- 0.25
DISP_PRIOR_VAR_MIN <- as.numeric(Sys.getenv("SEQUENCE_PRIOR_VAR_MIN", "0.01"))
stopifnot(is.finite(DISP_PRIOR_VAR_MIN), DISP_PRIOR_VAR_MIN > 0, DISP_PRIOR_VAR_MIN <= DESEQ2_PRIOR_VAR_FLOOR)
DISP_PRIOR_VAR_MAX <- 16
TREND_MAXIT <- 1000L
# Grid width 0.25 on natural-log scale; endpoint cells absorb tail mass. The
# grid spans 1e-6 to 1e3: below 1e-6 the NB likelihood has reached its Poisson
# plateau, so the lowest cell (which absorbs all prior mass below it) is exact
# to grid accuracy; above 1e3 the likelihood is negligible for any PAS.
# Rerun with SEQUENCE_DISP_GRID_STEP=0.125 as a numerical sensitivity analysis.
DISP_GRID_STEP <- as.numeric(Sys.getenv("SEQUENCE_DISP_GRID_STEP", "0.25"))
stopifnot(is.finite(DISP_GRID_STEP), DISP_GRID_STEP > 0, DISP_GRID_STEP <= 0.5)
DISP_LOG_GRID <- seq(log(1e-6), log(1e3), by = DISP_GRID_STEP)
DISP_LOG_EDGES <- c(-Inf, head(DISP_LOG_GRID, -1) + diff(DISP_LOG_GRID)/2, Inf)

# Return the positive mean-dispersion function expected by DESeq2, with metadata.
make_trend_function <- function(a0, a1, b) {
  force(a0); force(a1); force(b)
  trend <- function(means) a0 + a1 * means^(-b)
  attr(trend, "fitType") <- "hierarchical power trend"
  attr(trend, "coefficients") <- c(a0 = a0, a1 = a1, b = b)
  trend
}

# Compute per-PAS adjusted likelihood curves once. The design must be a full
# rank saturated group-indicator model (one fitted mean per biological arm).
# Under that design det(X'WX) equals the product of the group weight sums,
# up to a constant independent of dispersion for a change of design basis.
# Relative likelihoods remove arbitrary feature-specific log constants.
make_dispersion_likelihood <- function(dds) {
  X <- stats::model.matrix(DESeq2::design(dds), as.data.frame(SummarizedExperiment::colData(dds)))
  if (qr(X)$rank != ncol(X) || nrow(X) <= ncol(X)) stop("Invalid dispersion model design.")
  if ("weights" %in% SummarizedExperiment::assayNames(dds)) stop("Weighted counts are unsupported by the custom likelihood.")
  grp <- factor(apply(X, 1L, paste, collapse = ":"))
  if (nlevels(grp) != ncol(X)) stop("Custom adjusted likelihood requires a saturated categorical arm design.")
  ids <- rownames(dds)
  y <- DESeq2::counts(dds)
  mu <- SummarizedExperiment::assays(dds)[["mu"]]
  if (is.null(mu)) stop("DESeq2 did not retain the fitted means; check the installed version.")
  mc <- S4Vectors::mcols(dds)
  use <- !mc$allZero & is.finite(mc$baseMean) & mc$baseMean > 0 &
    apply(mu, 1L, function(z) all(is.finite(z) & z > 0))
  use[is.na(use)] <- FALSE
  y <- y[use, , drop = FALSE]; mu <- mu[use, , drop = FALSE]
  if (!nrow(y)) stop("No informative PASs for the adjusted likelihood.")
  ll <- vapply(exp(DISP_LOG_GRID), function(a) {
    z <- rowSums(matrix(stats::dnbinom(as.vector(y), mu = as.vector(mu), size = 1/a, log = TRUE), nrow(y)))
    w <- mu/(1 + a*mu)
    for (g in levels(grp)) z <- z - 0.5 * log(rowSums(w[, grp == g, drop = FALSE]))
    z
  }, numeric(nrow(y)))
  ll <- matrix(ll, nrow = nrow(y), dimnames = list(ids[use], NULL))
  maxll <- apply(ll, 1L, max)
  if (any(!is.finite(maxll))) stop("Nonfinite adjusted likelihood curve.")
  exp(ll - maxll)
}

# Integrate the relative adjusted likelihood against a discretized log-normal
# prior. Normal CDF bin probabilities sum to one, including boundary cells.
# A continuous rank covariate can be supplied for an exploratory comparator.
# Optimize in log coefficients and log prior variance, enforcing the DESeq2
# variance floor. Analytic derivatives reduce the cost of repeated cutoff fits.
fit_marginal_trend <- function(lx, likelihood, start = NULL, basis = NULL, multistart = FALSE, fallback = NULL) {
  if (length(lx) != nrow(likelihood) || !length(lx)) stop("Likelihood/mean mismatch.")
  if (is.null(basis)) basis <- matrix(numeric(), nrow = length(lx), ncol = 0L)
  d <- ncol(basis)
  lo <- c(log(TREND_MIN_COEF), log(TREND_MIN_COEF), TREND_B_RANGE[1], log(DISP_PRIOR_VAR_MIN), rep(-5, d))
  hi <- c(log(TREND_MAX_COEF), TREND_B_RANGE[2], log(DISP_PRIOR_VAR_MAX), rep(5, d))
  cache <- new.env(parent = emptyenv())
  calc <- function(q) {
    if (!is.null(cache$q) && identical(q, cache$q)) return(cache$out)
    a <- exp(q[1]); v <- exp(q[2] - q[3]*lx); m <- log(a+v)
    if (d) m <- m + drop(basis %*% q[4L + seq_len(d)])
    sd <- exp(q[4]/2)
    # Boundary derivatives are zero because the CDF at +/-Inf is constant.
    z <- outer(m, DISP_LOG_EDGES, function(m, e) (e-m)/sd)
    F <- stats::pnorm(z); D <- stats::dnorm(z); ZD <- z*D; ZD[!is.finite(ZD)] <- 0
    j <- seq_len(ncol(likelihood)); k <- j+1L
    w <- F[, k, drop=FALSE]-F[, j, drop=FALSE]
    den <- pmax(rowSums(likelihood*w), 1e-300)
    dm <- rowSums(likelihood*(D[, j, drop=FALSE]-D[, k, drop=FALSE]))/(sd*den)
    dv <- rowSums(likelihood*(ZD[, j, drop=FALSE]-ZD[, k, drop=FALSE]))/(2*den)
    g <- c(sum(dm*a/(a+v)), sum(dm*v/(a+v)), sum(dm*(-lx)*v/(a+v)), sum(dv))
    if (d) g <- c(g, drop(crossprod(basis, dm)))
    out <- list(value = -sum(log(den)), grad = -g)
    cache$q <- q; cache$out <- out; out
  }
  st <- if (is.null(start)) c(log(0.1), 0, 1, log(0.5)) else start[1:4]
  starts <- list(c(st, rep(0,d)))
  if (multistart || is.null(start)) starts <- c(starts, list(c(log(0.02), log(2), 1.5, log(1), rep(0,d)),
                                                              c(log(0.5), log(0.1), 0.7, log(0.25), rep(0,d))))
  fits <- lapply(starts, function(st) tryCatch(stats::optim(pmin(hi,pmax(lo,st)),
    function(q) calc(q)$value, function(q) calc(q)$grad, method="L-BFGS-B", lower=lo, upper=hi,
    control=list(maxit=TREND_MAXIT, factr=1e7)), error=function(e) NULL))
  fits <- Filter(function(o) !is.null(o) && o$convergence==0L && is.finite(o$value), fits)
  # A warm start that fails to converge is retried from the fallback (pooled) start.
  if (!length(fits) && !is.null(fallback))
    return(fit_marginal_trend(lx, likelihood, start = fallback, basis = basis, multistart = TRUE))
  if (!length(fits)) return(list(converged=FALSE, loglik=-Inf, par=rep(NA_real_,4+d), n=length(lx), tau2=NA_real_))
  o <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  list(converged=TRUE, loglik=-o$value, par=o$par, n=length(lx), tau2=exp(o$par[4]))
}

# Fit the identical hierarchical objective to the comparison counts before MAP
# shrinkage. The prior variance is estimated directly, not by subtracting a
# common trigamma sampling variance.
set_dispersion_trend <- function(dds) {
  lk <- make_dispersion_likelihood(dds)
  means <- S4Vectors::mcols(dds)$baseMean[match(rownames(lk), rownames(dds))]
  f <- fit_marginal_trend(log(means), lk, multistart=TRUE)
  if (!f$converged) stop("Hierarchical dispersion trend failed to converge.")
  p <- f$par
  trend <- make_trend_function(exp(p[1]), exp(p[2]), p[3])
  attr(trend, "hierarchical_prior_var") <- f$tau2
  attr(trend, "b_at_range_edge") <- min(abs(p[3]-TREND_B_RANGE)) < 1e-5
  attr(trend, "coef_at_bound") <- any(exp(p[1:2]) <= TREND_MIN_COEF*1.001 | exp(p[1:2]) >= TREND_MAX_COEF*.999)
  attr(trend, "n_exponents_fitted") <- NA_integer_
  attr(trend, "log_sampling_bias") <- 0
  DESeq2::`dispersionFunction<-`(dds, value=trend)
}

# DESeq2's dispersion steps with the trend above: size factors, gene-wise
# dispersions, trend, MAP shrinkage. These are the steps DESeq() runs, with the
# trend fit replaced. Size factors: the comparison's shared factors, estimated
# once on the full comparison matrix, so a PAS has identical normalized counts
# in the Original and in its EVS partition. Supplied factors are never
# re-estimated downstream (NULL, used only for a whole comparison matrix,
# estimates them from that matrix).
apply_size_factors <- function(dds, size_factors = NULL) {
  if (is.null(size_factors)) return(DESeq2::estimateSizeFactors(dds))
  if (!setequal(names(size_factors), colnames(dds)) || any(!is.finite(size_factors)) || any(size_factors <= 0)) {
    stop("Shared size factors do not match the samples of this dataset.")
  }
  sizeFactors(dds) <- unname(size_factors[colnames(dds)])
  dds
}

# Shared size factors of a comparison: DESeq2 median-of-ratios on the full
# comparison matrix, used by the Original and every EVS view.
comparison_size_factors <- function(raw, coldata, formula = make_design_formula(coldata)) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    coerce_raw_count_matrix_for_deseq2(raw, context = "shared size factors"),
    colData = coldata, design = formula)
  sf <- DESeq2::sizeFactors(tryCatch(DESeq2::estimateSizeFactors(dds),
    error = function(e) { warning("Ratio normalization failed; using poscounts: ", conditionMessage(e)); DESeq2::estimateSizeFactors(dds, type = "poscounts") }))
  stats::setNames(as.numeric(sf), colnames(raw))
}

# Gene-wise dispersions computed once per comparison.
# A PAS's gene-wise (Cox-Reid) dispersion depends only on its own counts, the
# design and the size factors, so with shared size factors it is the same in
# the Original and in any EVS partition: it is estimated once on the full
# comparison matrix and reused for every split. Only the trend, the prior and
# the MAP dispersions depend on the split. Every run checks the equality on its
# final analysis (validate_track_normalization).
comparison_gene_dds <- function(raw, coldata, size_factors) {
  if (is.null(size_factors)) return(NULL)
  dds <- DESeq2::DESeqDataSetFromMatrix(coerce_raw_count_matrix_for_deseq2(raw, context = "gene-wise dispersions"),
                                        colData = coldata, design = make_design_formula(coldata))
  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
  dds <- apply_size_factors(dds, size_factors)
  suppressMessages(DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE))
}

# Purpose: Estimate dispersions with trend.
estimate_dispersions_with_trend <- function(dds, size_factors = NULL, gene_dds = NULL, prior_floor = DISP_PRIOR_VAR_MIN) {
  if (!is.null(gene_dds)) {
    if (!all(rownames(dds) %in% rownames(gene_dds))) stop("Precomputed gene-wise dispersions lack some PASs.")
    dds <- gene_dds[rownames(dds), ]
  } else {
    dds <- apply_size_factors(dds, size_factors)
    dds <- suppressMessages(DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE))
  }
  dds <- set_dispersion_trend(dds)
  # Store the prior variance on the trend (as DESeq2 does internally), then
  # pass it to the MAP step so DESeq2 does not re-estimate it.
  prior_var <- dispersion_prior_variance(dds, prior_floor)
  f <- DESeq2::dispersionFunction(dds)
  attr(f, "dispPriorVar") <- prior_var
  attr(f, "dispPriorVar_floor") <- prior_floor
  dds <- DESeq2::`dispersionFunction<-`(dds, value = f, estimateVar = FALSE)
  suppressMessages(DESeq2::estimateDispersionsMAP(dds, dispPriorVar = prior_var, quiet = TRUE))
}

# The fitted latent variance tau^2, bounded below by prior_floor
# (DISP_PRIOR_VAR_MIN in the primary analysis; 0.25 in the DESeq2-floor reference).
dispersion_prior_variance <- function(dds, prior_floor = DISP_PRIOR_VAR_MIN) {
  v <- attr(DESeq2::dispersionFunction(dds), "hierarchical_prior_var")
  if (length(v)!=1L || !is.finite(v)) stop("Missing hierarchical prior variance.")
  max(v, prior_floor)
}

# =============================================================================
# DIFFERENTIAL EXPRESSION: MODEL STAGE (per analysis view) AND SIGNIFICANCE
# STAGE (per analysis view)
# =============================================================================
# Each analysis view (Original; Leading Edge; Remainder) is its own DESeq2
# dataset with its own dispersion trend, prior and shrinkage (size factors are
# shared from the full comparison matrix)
# (model stage, fit_view_model), and its own significance calibration
# (significance stage, call_significance): fdrtool fits an empirical null to
# the view's own Wald statistics, and higher criticism sets one HBFSS threshold
# per comparison from the Original's empirical p-values, applied to every view. The Leading Edge and Remainder are
# different populations (different expression ranges, dispersion laws and
# fractions of changed PASs), and a null shared by the two would be too strict
# for one and too loose for the other (Efron 2008, "Simultaneous inference:
# when should hypothesis testing problems be combined?").

# Model stage: DESeq2 with the power trend, the three Wald tests, apeglm.
fit_view_model <- function(count_mat, coldata, dataset_name, size_factors = NULL, gene_dds = NULL,
                           cooks_outliers = NULL, trend = c("power", "deseq2_default"),
                           prior_floor = DISP_PRIOR_VAR_MIN) {
  trend <- match.arg(trend)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = coerce_raw_count_matrix_for_deseq2(
      count_mat,
      context = paste0(dataset_name, " DESeq2 input after EVS split")
    ),
    colData = coldata,
    design = make_design_formula(coldata)
  )
  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]

  # DESeq2 with the power trend fitted above, then the Wald test (no beta prior;
  # with 5 replicates per group DESeq() would not replace outliers either).
  # trend = "deseq2_default" (reference baseline only): DESeq2's own
  # parametric gamma-GLM trend, with DESeq2's automatic fall-back to a local fit.
  dds <- if (identical(trend, "power")) estimate_dispersions_with_trend(dds, size_factors, gene_dds, prior_floor) else
    suppressMessages(DESeq2::estimateDispersions(apply_size_factors(dds, size_factors), fitType = "parametric", quiet = TRUE))
  dds <- DESeq2::nbinomWaldTest(dds, quiet = TRUE)
  SIGNIFICANCE_CALLS$wald <- SIGNIFICANCE_CALLS$wald + 1L

  # The three DESeq2 tests with independentFiltering = FALSE and
  # cooksCutoff = FALSE. The Cook's outlier rule is applied below with ONE
  # outlier list per comparison (the Original's), and the BH correction is
  # applied later, once per analysis (call_significance), so the Original and
  # the EVS views test the same PASs under the same rules.
  contrast <- c("condition", "trt", "untrt")
  tests <- c(list(std = list(alpha = BH_FDR_STANDARD)),
             if (WEAK_DELTA > 0) list(delta = list(lfcThreshold = WEAK_DELTA, altHypothesis = "greaterAbs", alpha = BH_FDR_WEAK)),
             list(strong = list(lfcThreshold = lfc_boundary, altHypothesis = "greaterAbs", alpha = BH_FDR_STRONG),
                  weak = list(lfcThreshold = lfc_boundary, altHypothesis = "lessAbs", alpha = BH_FDR_WEAK)))
  prim <- lapply(tests, function(a) do.call(DESeq2::results, c(list(dds, contrast = contrast, pAdjustMethod = "BH",
                                                                    independentFiltering = FALSE, cooksCutoff = FALSE), a)))

  res_df <- as.data.frame(prim$std)
  res_df$feature_id <- as.character(rownames(res_df))
  # Cook's outliers by DESeq2's default rule (Cook's distance above the 0.99
  # quantile of F(p, m - p)), detected in the Original's fit (cooks_outliers =
  # NULL) and supplied as that list to every EVS view of the comparison.
  if (is.null(cooks_outliers)) {
    cooks_p <- DESeq2::results(dds, contrast = contrast, independentFiltering = FALSE)$pvalue
    res_df$cooks_outlier <- is.na(cooks_p) & !is.na(prim$std$pvalue)
  } else {
    res_df$cooks_outlier <- res_df$feature_id %in% cooks_outliers
  }
  bh_view <- function(p) { q <- rep(NA_real_, length(p)); ok <- is.finite(p); q[ok] <- stats::p.adjust(p[ok], "BH"); q }
  res_df$pvalue[res_df$cooks_outlier] <- NA_real_
  res_df$padj <- bh_view(res_df$pvalue)      # replaced by the analysis-wide BH in call_significance

  shr <- DESeq2::lfcShrink(dds, coef = get_condition_coef(dds), type = "apeglm")
  shr_df <- data.frame(feature_id = rownames(shr), lfc_shrunk = shr$log2FoldChange)
  res_df <- dplyr::left_join(res_df, shr_df, by = "feature_id")

  res_df$regulation_direction <- ifelse(
    is.na(res_df$lfc_shrunk), NA_character_,
    ifelse(res_df$lfc_shrunk > 0, "upregulated",
           ifelse(res_df$lfc_shrunk < 0, "downregulated", "no_change")))

  i_ga <- match(res_df$feature_id, rownames(prim$strong))
  i_la <- match(res_df$feature_id, rownames(prim$weak))
  res_df$resGA_pvalue <- ifelse(res_df$cooks_outlier, NA_real_, prim$strong$pvalue[i_ga])
  # Intersection-union p-value: reject both zero effect and >= boundary.
  res_df$equivalence_pvalue <- ifelse(res_df$cooks_outlier, NA_real_, prim$weak$pvalue[i_la])
  # Lower side of Weak: |LFC| > 0 (ordinary test) or |LFC| > delta when set.
  lower_p <- if (WEAK_DELTA > 0) ifelse(res_df$cooks_outlier, NA_real_, prim$delta$pvalue[match(res_df$feature_id, rownames(prim$delta))]) else res_df$pvalue
  res_df$resLA_pvalue <- pmax(lower_p, res_df$equivalence_pvalue)
  res_df$resGA_padj <- bh_view(res_df$resGA_pvalue)
  res_df$resLA_padj <- bh_view(res_df$resLA_pvalue)

  list(dds = dds, res_df = res_df, dataset_name = dataset_name)
}

# Significance stage for one analysis (a named list of model-stage fits):
#   1. one fdrtool empirical null per view, fitted to that view's own Wald
#      statistics, so each view's empirical p-values are calibrated to its own
#      population (the Leading Edge and the Remainder differ);
#   2. one higher-criticism threshold per analysis, from the empirical p-values
#      of all its views together (the Original alone; the Leading Edge and the
#      Remainder combined), as the BH correction is one per analysis. Every
#      view of an analysis therefore has the same HCp and Htau, and the Original
#      and the split each get one threshold from one family of all PASs;
#   3. each view's HBFSS scores and calls.
# Returns the per-view results and one summary row per view.
# hc_p_fixed: when given (the EVS analysis), the HC threshold of the
# comparison's Original is used instead of a new one, so the Original and the
# Leading Edge + Remainder share one HBFSS threshold per comparison.
call_significance <- function(fits, analysis_label, annot_df, lean = FALSE, hc_p_fixed = NULL) {
  SIGNIFICANCE_CALLS$null <- SIGNIFICANCE_CALLS$null + 1L
  fits <- combine_bh_across_views(fits)
  null_rows <- list()
  for (i in seq_along(fits)) {
    f <- fits[[i]]
    view_label <- paste(analysis_label, names(fits)[i])
    valid <- is.finite(f$res_df$stat) & is.finite(f$res_df$pvalue) & !f$res_df$cooks_outlier
    if (sum(valid) < 5L) {
      fdr_fit <- NULL
    } else fdr_fit <- tryCatch(run_empirical_null_fdrtool(f$res_df$stat[valid], dataset_name = view_label),
      error=function(e) {warning(view_label, ": exploratory null fit failed: ",conditionMessage(e)); NULL})
    if (is.null(fdr_fit)) {
      fits[[i]]$res_df$empirical_p <- NA_real_
      null_rows[[i]] <- data.frame(View=names(fits)[i],n_wald_statistics=sum(valid),null_sd=NA_real_,null_fraction_eta0=NA_real_)
      next
    }
    if (length(fdr_fit$pval) != sum(valid)) {
      stop(sprintf("[%s] fdrtool empirical-p output length mismatch.", view_label))
    }
    param <- fdr_fit$param
    null_rows[[i]] <- data.frame(
      View = names(fits)[i],
      n_wald_statistics = sum(valid),
      null_sd = if ("sd" %in% colnames(param)) unname(param[1, "sd"]) else NA_real_,
      null_fraction_eta0 = if ("eta0" %in% colnames(param)) unname(param[1, "eta0"]) else NA_real_)
    fits[[i]]$res_df$empirical_p <- NA_real_
    fits[[i]]$res_df$empirical_p[valid] <- as.numeric(fdr_fit$pval)
  }
  all_p <- unlist(lapply(fits, function(f) f$res_df$empirical_p), use.names = FALSE)
  all_p <- all_p[is.finite(all_p)]
  hc_source <- if (is.null(hc_p_fixed)) "this analysis" else "the comparison's Original"
  hc_p <- if (is.null(hc_p_fixed)) safe_hc_thresh(all_p, dataset_name = paste(analysis_label, "(all views)")) else hc_p_fixed
  hc_ok <- !is.na(hc_p) && is.finite(hc_p) && hc_p > 0 && hc_p <= 1
  htau <- if (hc_ok) abs(log10(hc_p)) * lfc_boundary else NA_real_
  null_rows <- do.call(rbind, null_rows)
  null_rows$hc_p_threshold <- hc_p
  null_rows$hc_threshold_source <- hc_source
  null_rows$hbfss_threshold <- htau
  null_rows$n_hc_empirical_p <- length(all_p)
  views <- list()
  for (i in seq_along(fits)) views[[names(fits)[i]]] <- classify_view(fits[[i]], hc_p, htau, annot_df, lean = lean)
  list(views = views, null = null_rows)
}

# One BH correction per analysis for each DESeq2 test: over the Original's raw
# p-values, or over the Leading Edge and Remainder raw p-values combined (one
# hypothesis family; the views are disjoint). No independent filtering.
combine_bh_across_views <- function(fits) {
  n <- vapply(fits, function(f) nrow(f$res_df), integer(1))
  idx <- rep(seq_along(fits), n)
  for (col in c("pvalue", "resGA_pvalue", "resLA_pvalue")) {
    p <- unlist(lapply(fits, function(f) f$res_df[[col]]), use.names = FALSE)
    q <- rep(NA_real_, length(p))
    ok <- is.finite(p)
    q[ok] <- stats::p.adjust(p[ok], method = "BH")
    for (i in seq_along(fits)) fits[[i]]$res_df[[sub("pvalue", "padj", col)]] <- q[idx == i]
  }
  fits
}

# Reference baseline: the Original with DESeq2's default parametric trend
# (gamma GLM, b = 1, DESeq2's automatic local fall-back), all other steps and
# testing rules unchanged, next to the Original with the power trend. Shows
# how much of any EVS gain over standard DESeq2 practice comes from the trend.
record_default_trend_baseline <- function(comparison_name, raw, coldata, annot_df, sf, power_out) {
  fit <- fit_view_model(raw, coldata, paste(comparison_name, "Original (DESeq2 default trend)"),
                        size_factors = sf, cooks_outliers = original_cooks_cache[[comparison_name]], trend = "deseq2_default")
  sig <- call_significance(list(raw_dataset = fit), paste(comparison_name, "Original (DESeq2 default trend)"),
                           annot_df, lean = TRUE)
  row <- function(label, res, dds, nl, fit_type) {
    cnt <- count_discoveries(discovery_flags(list(raw_dataset = res)))
    f <- DESeq2::dispersionFunction(dds)
    X <- stats::model.matrix(DESeq2::design(dds), as.data.frame(SummarizedExperiment::colData(dds)))
    data.frame(Comparison = comparison_name, Analysis = label, Trend = fit_type,
               Median_Final_Dispersion = stats::median(S4Vectors::mcols(dds)$dispersion, na.rm = TRUE),
               Share_PAS_Below_Trend = mean(S4Vectors::mcols(dds)$dispGeneEst < S4Vectors::mcols(dds)$dispFit, na.rm = TRUE),
               Prior_Variance = as.numeric(attr(f, "dispPriorVar")),
               Null_SD = nl$null_sd[1], Htau = nl$hbfss_threshold[1],
               Standard = cnt[["Standard"]], Strong = cnt[["Strong"]], Weak = cnt[["Weak"]], HBFSS = cnt[["HBFSS"]])
  }
  ft <- attr(DESeq2::dispersionFunction(sig$views$raw_dataset$dds), "fitType")
  default_trend_rows[[comparison_name]] <<- rbind(
    row(ORIGINAL_LABEL, power_out$views$raw_dataset$results, power_out$views$raw_dataset$dds,
        power_out$null, "power trend (this pipeline)"),
    row("Original (DESeq2 default trend)", sig$views$raw_dataset$results, sig$views$raw_dataset$dds, sig$null,
        paste("DESeq2", if (is.null(ft)) "parametric" else ft)))
  invisible(TRUE)
}

# The Original's Cook's outlier list of a comparison (fitted once; the EVS views
# of that comparison use it instead of their own).
remember_original_cooks <- function(comparison_name, results) {
  original_cooks_cache[[comparison_name]] <<- as.character(results$feature_id[results$cooks_outlier %in% TRUE])
  invisible(original_cooks_cache[[comparison_name]])
}
# Purpose: Original cooks outliers.
original_cooks_outliers <- function(comparison_name, raw, coldata, size_factors) {
  if (!is.null(original_cooks_cache[[comparison_name]])) return(original_cooks_cache[[comparison_name]])
  fit <- fit_view_model(raw, coldata, paste(comparison_name, "Original (Cook's outliers)"), size_factors = size_factors)
  remember_original_cooks(comparison_name, fit$res_df)
}

# Purpose: Run core analysis.
run_core_analysis <- function(count_mats, coldata, dataset_names, annot_df, analysis_label, size_factors = NULL,
                              gene_dds = NULL, cooks_outliers = NULL, hc_p_fixed = NULL) {
  fits <- Map(function(m, nm) fit_view_model(m, coldata, nm, size_factors = size_factors, gene_dds = gene_dds,
                                             cooks_outliers = cooks_outliers), count_mats, dataset_names)
  names(fits) <- names(count_mats)
  call_significance(fits, analysis_label, annot_df, hc_p_fixed = hc_p_fixed)
}

# HBFSS scores and the Standard / Strong / Weak / HBFSS calls for one view,
# given its empirical p-values and the analysis's HC threshold.
# lean = TRUE (used only by the default-trend baseline) returns the decision table
# without normalized counts, annotation and dispersion columns; every flag and
# runtime check is identical. The final analysis always uses lean = FALSE.
classify_view <- function(fit, hc_p_threshold, hbfss_threshold, annot_df, lean = FALSE) {
  dds <- fit$dds
  res_df <- fit$res_df
  dataset_name <- fit$dataset_name

  res_df$HBFSS <- abs(res_df$lfc_shrunk) * (-log10(pmax(res_df$empirical_p, 1e-300)))
  res_df$HBFSS_core_pass <- if (is.na(hbfss_threshold)) FALSE else
    !is.na(res_df$HBFSS) & is.finite(res_df$HBFSS) & res_df$HBFSS > hbfss_threshold

  # Formal BH families are separate. The primary family is ordinary Wald/BH;
  # strong and weak effect tests are secondary. Shrunken LFCs are estimates,
  # not extra rejection gates. HBFSS never defines primary significance.
  finite_lfc <- is.finite(res_df$lfc_shrunk)
  abs_shrunk_lfc <- abs(res_df$lfc_shrunk)
  eligible <- !res_df$cooks_outlier & is.finite(res_df$pvalue)
  res_df$standard_flag <- eligible & !is.na(res_df$padj) & res_df$padj < BH_FDR_STANDARD
  # Descriptive only (not FDR-controlled): Standard calls with |apeglm LFC| >= 1,
  # the definition used in earlier versions. "Strong" is the formal test of |LFC| > 1.
  res_df$standard_lfc1_descriptive <- res_df$standard_flag & finite_lfc & abs_shrunk_lfc >= lfc_boundary
  res_df$strong_cnh_flag <- eligible & !is.na(res_df$resGA_padj) & res_df$resGA_padj < BH_FDR_STRONG
  res_df$weak_cnh_flag <- eligible & !is.na(res_df$resLA_padj) & res_df$resLA_padj < BH_FDR_WEAK
  res_df$weak_significant_flag <- res_df$weak_cnh_flag
  res_df$hbfss_flag <- eligible & finite_lfc & res_df$HBFSS_core_pass
  res_df$standard_hbfss_overlap <- res_df$standard_flag & res_df$hbfss_flag
  res_df$strong_hbfss_overlap <- res_df$strong_cnh_flag & res_df$hbfss_flag
  res_df$weak_hbfss_overlap <- res_df$weak_cnh_flag & res_df$hbfss_flag
  res_df$any_overlap <- res_df$hbfss_flag & (res_df$standard_flag | res_df$strong_cnh_flag | res_df$weak_cnh_flag)
  res_df$final_significant_flag <- res_df$standard_flag
  res_df$secondary_or_exploratory_flag <- res_df$strong_cnh_flag | res_df$weak_significant_flag | res_df$hbfss_flag
  stopifnot(!any(res_df$cooks_outlier & (res_df$final_significant_flag | res_df$secondary_or_exploratory_flag)),
            all(is.na(res_df$empirical_p[res_df$cooks_outlier])))
  # Independent thresholds and effect families need not yield nested sets. No
  # FDR claim is made for the union of the three formal families or HBFSS.


  if (lean) {
    res_df$dataset_name <- dataset_name
    res_df$hc_p_threshold_dataset <- hc_p_threshold
    res_df$hbfss_threshold_dataset <- hbfss_threshold
    return(list(dds = dds, results = res_df))
  }

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
  final_df$hc_p_threshold_dataset <- hc_p_threshold
  final_df$hbfss_threshold_dataset <- hbfss_threshold

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
    "standard_lfc1_descriptive",
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

  list(dds = dds, results = final_df)
}


# Purpose: Build final volcano df.
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

# Purpose: Select final volcano labels.
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


# Purpose: Make hbfss boundary df.
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

# Purpose: Plot final volcano.
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

# Purpose: Run one evs track.
run_one_evs_track <- function(comparison_name, track_key, count_matrix, coldata, annot_df) {
  evs <- build_eigenvector_split(
    comparison_name = comparison_name,
    count_matrix = count_matrix,
    coldata = coldata,
    track_key = track_key
  )

  message(
    comparison_name, " ", unname(track_short[track_key]),
    ": EVS k=", evs$empirical_evs_k,
    " selected by held-out likelihood using each condition PC1 eigenvector; ",
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

  analysis_dataset_list <- dataset_list

  # Two analyses: the Original, then the Leading Edge with the Remainder. Each
  # view gets its own empirical null; each analysis one BH. One HC / HBFSS
  # threshold per comparison: set by the Original and reused by the EVS
  # analysis (call_significance, hc_p_fixed).
  analyses <- list(original = "raw_dataset", split = c("leading_edge_dataset", "remainder_dataset"))

  # Shared size factors and gene-wise dispersions: estimated once on the full
  # comparison matrix (prepare_comparison_cache) and reused by every view.
  cache <- comparison_cache[[comparison_name]]
  if (is.null(cache) || !identical(rownames(cache$raw), rownames(evs$raw_dataset)))
    stop(comparison_name, ": the comparison cache does not match the analysed matrix")
  sf <- cache$sf

  analysis_results <- list()
  original_hc_p <- NULL
  for (a in names(analyses)) {
    keys <- analyses[[a]]
    analysis_label <- if (a == "original") ORIGINAL_LABEL else paste(unname(track_short[track_key]), "LE + REM")
    # Cook's outliers: the Original uses its own; the EVS views use the Original's.
    cooks <- if (a == "original") NULL else original_cooks_outliers(comparison_name, evs$raw_dataset, coldata, sf)
    out <- run_core_analysis(
      count_mats = analysis_dataset_list[keys],
      coldata = coldata,
      dataset_names = paste(comparison_name, unname(track_short[track_key]), keys, sep = "_"),
      annot_df = annot_df,
      analysis_label = paste(comparison_name, analysis_label),
      size_factors = sf,
      gene_dds = cache$gene_dds,
      cooks_outliers = cooks,
      hc_p_fixed = if (a == "original") NULL else original_hc_p
    )
    if (a == "original") {
      original_hc_p <- out$null$hc_p_threshold[1]
      if (is.null(original_hc_p) || !is.finite(original_hc_p)) original_hc_p <- NA_real_
      remember_original_cooks(comparison_name, out$views$raw_dataset$results)
      # Reference baseline: the Original with DESeq2's default trend (once).
      tryCatch(record_default_trend_baseline(comparison_name, evs$raw_dataset, coldata, annot_df, sf, out),
               error = function(e) warning(comparison_name, ": DESeq2 default-trend baseline failed: ",
                                           conditionMessage(e), call. = FALSE))
    } else {
      add_validation("One HBFSS threshold per comparison (Original and EVS share the Original's HC threshold)",
                     comparison_name, identical(unique(out$null$hc_p_threshold), original_hc_p),
                     sprintf("HC p = %.4g; Htau = %.3f", original_hc_p, out$null$hbfss_threshold[1]))
    }
    for (nm in keys) {
      analysis_results[[nm]] <- list(dds = out$views[[nm]]$dds, results = out$views[[nm]]$results,
                                     dataset_mat = analysis_dataset_list[[nm]])
    }
    null_tab <- out$null
    null_tab$View <- vapply(null_tab$View, function(k) analysis_view_label(track_key, k), character(1))
    empirical_null_rows[[length(empirical_null_rows) + 1L]] <<- cbind(
      data.frame(Comparison = comparison_name, Analysis = analysis_label), null_tab)
  }

  tryCatch(run_sensitivity_analyses(comparison_name, evs, coldata, annot_df, sf, cache$gene_dds,
                                    original_cooks_cache[[comparison_name]], original_hc_p, analysis_results,
                                    cache$depth),
           error = function(e) warning(comparison_name, ": sensitivity analyses failed: ", conditionMessage(e), call. = FALSE))
  reg <- paper_registry[[registry_key(comparison_name, track_key)]]
  reg$analysis_results <- analysis_results
  reg$size_factors <- sf
  paper_registry[[registry_key(comparison_name, track_key)]] <<- reg
  validate_track_normalization(comparison_name, track_key, analysis_results, sf, evs, coldata)
  TRUE
}

# =============================================================================
# SENSITIVITY ANALYSES (per comparison, after the primary analysis)
# =============================================================================
# Counts only, under the primary testing rules (one BH per analysis, the
# Original's Cook's outliers and HBFSS threshold):
#   1. DESeq2 0.25 floor on the MAP prior variance (Original and EVS at k*);
#   2. EVS with an independent reference ranking built from the 30 held-out
#      samples only (no tested count enters the split), at its own cutoff;
#   3. IHW (Ignatiadis et al. 2016) on the Original's p-values with the PC1
#      union depth as covariate: the established covariate-weighting
#      alternative to splitting (run when the IHW package is installed).
sensitivity_rows <- list()

run_sensitivity_analyses <- function(comparison_name, evs, coldata, annot_df, sf, gene_dds, cooks, hc_p,
                                     analysis_results, depth) {
  count_row <- function(label, res_list, note, k = NA_integer_, le = NA_integer_) {
    cnt <- count_discoveries(discovery_flags(res_list))
    data.frame(Comparison = comparison_name, Analysis = label, k = k, Leading_Edge = le,
               PAS_tested = cnt[["n_tested"]], Standard = cnt[["Standard"]], Strong = cnt[["Strong"]],
               Weak = cnt[["Weak"]], HBFSS = cnt[["HBFSS"]], Note = note)
  }
  fit_lean <- function(mats, label, prior_floor) {
    fits <- Map(function(m, nm) fit_view_model(m, coldata, paste(comparison_name, label, nm), size_factors = sf,
                                               gene_dds = gene_dds, cooks_outliers = cooks, prior_floor = prior_floor),
                mats, names(mats))
    sig <- call_significance(fits, paste(comparison_name, label), annot_df, lean = TRUE, hc_p_fixed = hc_p)
    lapply(sig$views, `[[`, "results")
  }
  k <- get_active_evs_cutoff(comparison_name)
  rows <- list(
    count_row(ORIGINAL_LABEL, list(raw_dataset = analysis_results$raw_dataset$results), "primary"),
    count_row("EVS (primary)", list(leading_edge_dataset = analysis_results$leading_edge_dataset$results,
                                    remainder_dataset = analysis_results$remainder_dataset$results), "primary",
              k, length(evs$leading_edge_ids)))
  rows[[3]] <- count_row(paste0(ORIGINAL_LABEL, ", DESeq2 0.25 prior floor"),
                         fit_lean(list(raw_dataset = evs$raw_dataset), "Original floor 0.25", DESEQ2_PRIOR_VAR_FLOOR),
                         "reference: MAP prior variance max(tau2, 0.25)")
  rows[[4]] <- count_row("EVS, DESeq2 0.25 prior floor",
                         fit_lean(list(leading_edge_dataset = evs$leading_edge_dataset,
                                       remainder_dataset = evs$remainder_dataset), "EVS floor 0.25", DESEQ2_PRIOR_VAR_FLOOR),
                         "reference: MAP prior variance max(tau2, 0.25)", k, length(evs$leading_edge_ids))
  k_ref <- REFERENCE_K[[comparison_name]]
  ref <- CUTOFF$ref_depth[[comparison_name]]
  if (!is.null(k_ref) && is.finite(k_ref) && !is.null(ref)) {
    ids <- rownames(evs$raw_dataset)
    le <- ids[ref[ids] <= k_ref]; rem <- setdiff(ids, le)
    if (length(le) >= SPLIT_MIN_SET && length(rem) >= SPLIT_MIN_SET)
      rows[[length(rows) + 1L]] <- count_row("EVS, independent reference ranking",
        fit_lean(list(leading_edge_dataset = evs$raw_dataset[le, , drop = FALSE],
                      remainder_dataset = evs$raw_dataset[rem, , drop = FALSE]), "EVS reference ranking", DISP_PRIOR_VAR_MIN),
        "ranking from the 30 held-out samples only (VST and PC1 on held-out arms)", as.integer(k_ref), length(le))
  }
  o <- analysis_results$raw_dataset$results
  if (requireNamespace("IHW", quietly = TRUE)) {
    ok <- is.finite(o$pvalue)
    ih <- IHW::ihw(o$pvalue[ok], as.numeric(depth[o$feature_id[ok]]), alpha = BH_FDR_STANDARD)
    rows[[length(rows) + 1L]] <- data.frame(Comparison = comparison_name, Analysis = "IHW, Original p-values",
      k = NA_integer_, Leading_Edge = NA_integer_, PAS_tested = sum(ok), Standard = IHW::rejections(ih),
      Strong = NA_integer_, Weak = NA_integer_, HBFSS = NA_integer_,
      Note = "covariate: PC1 union depth; alpha = Standard BH level")
  } else {
    rows[[length(rows) + 1L]] <- data.frame(Comparison = comparison_name, Analysis = "IHW, Original p-values",
      k = NA_integer_, Leading_Edge = NA_integer_, PAS_tested = NA_integer_, Standard = NA_integer_,
      Strong = NA_integer_, Weak = NA_integer_, HBFSS = NA_integer_, Note = "not run: IHW package not installed")
  }
  # 4. Negative binomial count splitting (Neufeld et al.): each tested count is
  #    split into two independent NB parts (beta-binomial thinning with the
  #    Original's final dispersion, epsilon = 0.5). The split is ranked on the
  #    training half (VST and PC1 of the ten training samples, the primary k*)
  #    and both the Original and EVS are tested on the other half, so the
  #    ranking is independent of the tested counts. Both use half the data.
  cs <- tryCatch(count_split_analysis(comparison_name, evs, coldata, annot_df, k,
                                      DESeq2::dispersions(analysis_results$raw_dataset$dds),
                                      rownames(analysis_results$raw_dataset$dds)),
                 error = function(e) { warning(comparison_name, ": count splitting failed: ", conditionMessage(e), call. = FALSE); NULL })
  if (!is.null(cs)) rows <- c(rows, cs)
  sensitivity_rows[[comparison_name]] <<- dplyr::bind_rows(rows)
  invisible(TRUE)
}

COUNT_SPLIT_EPSILON <- 0.5

count_split_analysis <- function(comparison_name, evs, coldata, annot_df, k, disp, disp_ids) {
  set.seed(SPLIT_SEED + match(comparison_name, comparison_table$comparison_name))
  y <- evs$raw_dataset
  alpha <- disp[match(rownames(y), disp_ids)]
  alpha[!is.finite(alpha) | alpha <= 0] <- stats::median(alpha[is.finite(alpha) & alpha > 0])
  size <- matrix(1 / pmax(alpha, 1e-8), nrow(y), ncol(y))
  p <- stats::rbeta(length(y), COUNT_SPLIT_EPSILON * size, (1 - COUNT_SPLIT_EPSILON) * size)
  train <- matrix(stats::rbinom(length(y), as.vector(y), p), nrow(y), dimnames = dimnames(y))
  test <- y - train
  storage.mode(train) <- "integer"; storage.mode(test) <- "integer"
  # Ranking on the training half only.
  tr <- train[rowSums(train) > 0, , drop = FALSE]
  dtr <- DESeq2::DESeqDataSetFromMatrix(tr, colData = coldata, design = make_design_formula(coldata))
  dtr <- tryCatch(DESeq2::estimateSizeFactors(dtr), error = function(e) DESeq2::estimateSizeFactors(dtr, type = "poscounts"))
  v <- as.matrix(SummarizedExperiment::assay(suppressMessages(DESeq2::varianceStabilizingTransformation(dtr, blind = TRUE))))
  pos <- lapply(c("trt", "untrt"), function(cnd) {
    r <- arm_pc1_ranks(v, colnames(v)[coldata[colnames(v), "condition"] == cnd]); nrow(v) - r + 1L })
  depth <- pmin(pos[[1]], pos[[2]])
  te <- test[rowSums(test) > 0, , drop = FALSE]
  d_te <- depth[rownames(te)]; d_te[is.na(d_te)] <- .Machine$integer.max
  le <- rownames(te)[d_te <= k]; rem <- setdiff(rownames(te), le)
  if (length(le) < SPLIT_MIN_SET || length(rem) < SPLIT_MIN_SET) stop("count-split partition too small")
  # Testing on the other half, with the primary rules.
  sf <- comparison_size_factors(te, coldata); gdds <- comparison_gene_dds(te, coldata, sf)
  lab <- paste(comparison_name, "count split")
  o <- fit_view_model(te, coldata, paste(lab, "Original"), size_factors = sf, gene_dds = gdds)
  cooks <- as.character(o$res_df$feature_id[o$res_df$cooks_outlier %in% TRUE])
  so <- call_significance(list(raw_dataset = o), paste(lab, "Original"), annot_df, lean = TRUE)
  hc <- so$null$hc_p_threshold[1]
  fits <- Map(function(m, nm) fit_view_model(m, coldata, paste(lab, nm), size_factors = sf, gene_dds = gdds,
                                             cooks_outliers = cooks),
              list(leading_edge_dataset = te[le, , drop = FALSE], remainder_dataset = te[rem, , drop = FALSE]),
              c("Lead", "Rem"))
  se <- call_significance(fits, paste(lab, "EVS"), annot_df, lean = TRUE, hc_p_fixed = if (is.finite(hc)) hc else NA_real_)
  row <- function(label, res_list, le_n) {
    cnt <- count_discoveries(discovery_flags(res_list))
    data.frame(Comparison = comparison_name, Analysis = label, k = if (is.na(le_n)) NA_integer_ else k,
               Leading_Edge = le_n, PAS_tested = cnt[["n_tested"]], Standard = cnt[["Standard"]],
               Strong = cnt[["Strong"]], Weak = cnt[["Weak"]], HBFSS = cnt[["HBFSS"]],
               Note = sprintf("NB count splitting, epsilon %.1f: ranked on one half, tested on the other; LE Jaccard with primary %.2f",
                              COUNT_SPLIT_EPSILON, length(intersect(le, evs$leading_edge_ids)) / length(union(le, evs$leading_edge_ids))))
  }
  list(row("Count split: Original (test half)", list(raw_dataset = so$views$raw_dataset$results), NA_integer_),
       row("Count split: EVS (test half)", lapply(se$views, `[[`, "results"), length(le)))
}

# Purpose: Analysis view label.
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

# Purpose: Get registered result.
get_registered_result <- function(comparison_name, track_key, dataset_key) {
  obj <- paper_registry[[registry_key(comparison_name, track_key)]]

  if (is.null(obj) || is.null(obj$analysis_results) ||
      is.null(obj$analysis_results[[dataset_key]])) {
    return(NULL)
  }

  obj$analysis_results[[dataset_key]]$results
}

# Purpose: Build support label.
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

# Purpose: Comparison analysis views.
comparison_analysis_views <- function() {
  data.frame(
    track_key = "normalized_evs",
    dataset_key = c("raw_dataset", "leading_edge_dataset", "remainder_dataset"),
    stringsAsFactors = FALSE
  )
}

# Purpose: Analysis view folder name.
analysis_view_folder_name <- function(track_key, dataset_key) {
  if (identical(dataset_key, "raw_dataset")) return("Original_No_EVS")
  paste0(
    unname(track_short[track_key]),
    "_",
    ifelse(dataset_key == "leading_edge_dataset", "Lead", "Rem")
  )
}

# Purpose: Analysis view paths.
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

# Purpose: Export comparison manuscript tables.
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
      Std_LFC1_descriptive = sum(df$standard_lfc1_descriptive, na.rm = TRUE),
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
    df$Selected_Model <- CUTOFF$summary$Selected_Model[match(comparison_name,CUTOFF$summary$Comparison)]
    df$Selected_For_Primary <- if(identical(dataset_key,"raw_dataset")) df$Selected_Model=="Pooled" else df$Selected_Model=="EVS"
    save_csv(df, file.path(paths$tables, "Table_All_PAS_Results.csv"))
    save_csv(df[df$secondary_or_exploratory_flag %in% TRUE,,drop=FALSE],
             file.path(paths$tables,"Table_Secondary_And_Exploratory_PAS.csv"))
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
    # comparison-level panels generated after all three analysis views are fit.
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

# Purpose: Build overall manuscript summary.
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
      Std_LFC1_descriptive = sum(df$standard_lfc1_descriptive, na.rm = TRUE),
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

# Purpose: Write methods note.
write_methods_note <- function() {
  writeLines(MANUSCRIPT_METHODS, file.path(output_dir,"Methods_Manuscript.md"), useBytes=TRUE)
  invisible(TRUE)
}

# =============================================================================
# DESeq2 DISPERSIONS AFTER THE SPLIT: SUMMARIES AND FIGURE 2
# =============================================================================
# Everything here is read from the fitted DESeq2 objects of each view; nothing
# is refitted.

DISP_MIN_GENE_EST <- 1e-8

# ---- per-run state -----------------------------------------------------------

dispersion_metrics_rows <- list()
dispersion_trend_rows <- list()
dispersion_pas_rows <- list()
composition_rows <- list()
active_cutoff_rows <- list()
empirical_null_rows <- list()

# ---- DESeq2 dispersions after the split (read from the fitted objects) ------

VIEW_SET <- c("Original (No EVS)" = "Original (No EVS)", "NormEVS Lead" = "Leading Edge",
              "NormEVS Rem" = "Remainder")

`%||%` <- function(a, b) if (is.null(a)) b else a

# Purpose: Dispersion view metrics.
dispersion_view_metrics <- function(dds) {
  mc <- S4Vectors::mcols(dds)
  fin <- DESeq2::dispersions(dds)
  ok <- !mc$allZero & is.finite(mc$dispGeneEst) & mc$dispGeneEst >= DISP_MIN_GENE_EST &
        is.finite(mc$dispFit) & mc$dispFit > 0 & is.finite(fin) & mc$baseMean > 0
  dfun <- DESeq2::dispersionFunction(dds)
  X <- stats::model.matrix(DESeq2::design(dds), as.data.frame(SummarizedExperiment::colData(dds)))
  v <- attr(dfun, "varLogDispEsts")
  outl <- if (is.null(mc$dispOutlier)) rep(FALSE, length(ok)) else (mc$dispOutlier %in% TRUE)
  data.frame(
    n_PAS = sum(ok),
    var_log_disp_residuals = if (length(v) == 1L) as.numeric(v) else NA_real_,
    sampling_var_log_disp = trigamma((ncol(dds) - ncol(X)) / 2),
    prior_var_used = as.numeric(attr(dfun, "dispPriorVar") %||% NA_real_),
    prior_var_hierarchical = as.numeric(attr(dfun, "hierarchical_prior_var") %||% NA_real_),
    prior_var_lower_bound_active = isTRUE(attr(dfun,"hierarchical_prior_var") <= DISP_PRIOR_VAR_MIN*1.00001),
    prior_var_below_deseq2_floor = isTRUE(attr(dfun,"hierarchical_prior_var") < DESEQ2_PRIOR_VAR_FLOOR),
    trend_a0 = attr(dfun, "coefficients")[["a0"]],
    trend_a1 = attr(dfun, "coefficients")[["a1"]],
    trend_b = attr(dfun, "coefficients")[["b"]],
    dispersion_outlier_pct = 100 * mean(outl[ok]),
    median_final_dispersion = stats::median(fin[ok]))
}

# Purpose: Collect dispersion after evs.
collect_dispersion_after_evs <- function(comparison_name) {
  views <- comparison_analysis_views()
  norm_reg <- paper_registry[[registry_key(comparison_name, "normalized_evs")]]
  for (i in seq_len(nrow(views))) {
    obj <- paper_registry[[registry_key(comparison_name, views$track_key[i])]]
    dds <- obj$analysis_results[[views$dataset_key[i]]]$dds
    if (is.null(dds)) next
    label <- analysis_view_label(views$track_key[i], views$dataset_key[i])
    m <- tryCatch(dispersion_view_metrics(dds), error = function(e) NULL)
    if (!is.null(m)) {
      dispersion_metrics_rows[[length(dispersion_metrics_rows) + 1L]] <<-
        cbind(data.frame(Comparison = comparison_name, View = label), m)
    }
    if (!label %in% names(VIEW_SET)) next
    mc <- S4Vectors::mcols(dds)
    keep <- !mc$allZero & mc$baseMean > 0
    outl <- if (is.null(mc$dispOutlier)) rep(FALSE, nrow(dds)) else (mc$dispOutlier %in% TRUE)
    dispersion_pas_rows[[length(dispersion_pas_rows) + 1L]] <<- data.frame(
      Comparison = comparison_name, View = label, feature_id = rownames(dds)[keep],
      baseMean = mc$baseMean[keep], dispGeneEst = mc$dispGeneEst[keep], dispFit = mc$dispFit[keep],
      dispersion = DESeq2::dispersions(dds)[keep], dispOutlier = outl[keep])
    bm <- mc$baseMean[keep]
    rng <- stats::quantile(bm, c(0.005, 0.995), names = FALSE)
    grid_bm <- 10^seq(log10(rng[1]), log10(rng[2]), length.out = 100)
    dispersion_trend_rows[[length(dispersion_trend_rows) + 1L]] <<- data.frame(
      Comparison = comparison_name, View = label, baseMean = grid_bm,
      trend = pmin(DESeq2::dispersionFunction(dds)(grid_bm), max(10, ncol(dds))))  # DESeq2's dispersion ceiling
  }
  evs <- norm_reg$evs
  if (!is.null(evs)) {
    composition_rows[[comparison_name]] <<- data.frame(
      Comparison = comparison_name, k = evs$empirical_evs_k,
      Class = c("Joint (top-k in both)", "Disjoint (top-k in one)"),
      n = c(length(evs$joint_ids), length(evs$disjoint_trt_ids) + length(evs$disjoint_untrt_ids)))
  }
  invisible(TRUE)
}

# ---- figure helpers ----------------------------------------------------------

# Purpose: Comp levels.
comp_levels <- function() as.character(comparison_table$comparison_name)
# Purpose: Log lab.
log_lab <- function() scales::trans_format("log10", scales::math_format(10^.x))
# Purpose: As comp.
as_comp <- function(x) factor(x, levels = comp_levels())

# ---- Figure 2: DESeq2's dispersions and the width of its prior ---------------

# Purpose: Set bars scale.
set_bars_scale <- function() {
  scale_fill_manual(values = SET_COLS, breaks = names(SET_COLS))
}

# DESeq2's dispersion plot for the Original, Leading Edge and Remainder:
# gene-wise estimates (grey), final MAP estimates (colour) and the fitted
# parametric trend a0 + a1 * mean^(-b) (orange, capped at DESeq2's maximum, 10).
panel_dispersion_grid <- function(pas, trd, met, tag, n_points = 3000L) {
  d <- pas[is.finite(pas$dispGeneEst) & pas$dispGeneEst >= DISP_MIN_GENE_EST & is.finite(pas$dispersion), ]
  set.seed(20261115L)
  d <- do.call(rbind, lapply(split(d, list(d$Comparison, d$View), drop = TRUE), function(x)
    x[sort(sample.int(nrow(x), min(nrow(x), n_points))), , drop = FALSE]))
  d$Set <- factor(unname(VIEW_SET[d$View]), levels = names(SET_COLS)); d$Comparison <- as_comp(d$Comparison)
  t <- trd; t$Set <- factor(unname(VIEW_SET[t$View]), levels = names(SET_COLS)); t$Comparison <- as_comp(t$Comparison)
  m <- met[met$View %in% names(VIEW_SET), ]
  m$Set <- factor(unname(VIEW_SET[m$View]), levels = names(SET_COLS)); m$Comparison <- as_comp(m$Comparison)
  m$label <- sprintf("n = %s\ntrend %.3f + %.2f x mean^-%.2f\nprior variance %.3f", format(m$n_PAS, big.mark = ","),
                     m$trend_a0, m$trend_a1, m$trend_b, m$prior_var_used)
  ggplot(d, aes(baseMean)) +
    geom_point(aes(y = dispGeneEst), colour = "grey72", size = 0.12, alpha = 0.5) +
    geom_point(aes(y = dispersion, colour = Set), size = 0.12, alpha = 0.6) +
    geom_line(data = t, aes(baseMean, trend), colour = "#E8A317", linewidth = 0.6) +
    geom_text(data = m, aes(x = Inf, y = Inf, label = label), hjust = 1.03, vjust = 1.15, size = 1.8, lineheight = 0.95) +
    facet_grid(Set ~ Comparison, labeller = labeller(Comparison = comp_lab())) +
    scale_colour_manual(values = SET_COLS, guide = "none") +
    scale_x_log10(labels = log_lab()) + scale_y_log10(labels = log_lab()) +
    coord_cartesian(ylim = c(1e-3, 20)) +
    labs(x = "Mean of normalized counts", y = "Dispersion", tag = tag,
         title = "DESeq2 dispersion estimates: Original vs Leading Edge and Remainder",
         subtitle = "Grey: gene-wise estimates. Colour: final (shrunken) estimates. Orange: parametric trend a0 + a1 x mean^-b, fitted using marginal adjusted likelihood.") +
    fig_theme() + theme(panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.3))
}

# Exponent b of the fitted trend a0 + a1 * mean^(-b), Original vs each NormEVS
# set: a different b means the Leading Edge and Remainder follow different
# dispersion laws.
panel_trend_exponent <- function(met, tag) {
  d <- met[met$View %in% names(VIEW_SET), c("Comparison", "View", "trend_b")]
  d$Set <- factor(unname(VIEW_SET[d$View]), levels = names(SET_COLS)); d$Comparison <- as_comp(d$Comparison)
  ggplot(d, aes(Comparison, trend_b, fill = Set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.74) +
    geom_text(aes(label = sprintf("%.2f", trend_b)), position = position_dodge(width = 0.8), vjust = -0.4, size = 1.8) +
    set_bars_scale() +
    scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Trend exponent b", tag = tag,
         title = "Shape of the dispersion trend",
         subtitle = "b in a0 + a1 x mean^-b for each set (b and a0 trade off; compare whole curves in panel a).") +
    fig_theme() + theme(panel.grid.major.x = element_blank())
}

# Dispersion prior variance used in MAP shrinkage: the latent variance tau^2
# estimated by the marginal likelihood (bounded below by DISP_PRIOR_VAR_MIN),
# for the Original and each NormEVS set; W = n-weighted pool of the two sets.
panel_prior_width <- function(met, tag) {
  d <- met[met$View %in% names(VIEW_SET), c("Comparison", "View", "n_PAS", "prior_var_used")]
  d$Set <- factor(unname(VIEW_SET[d$View]), levels = names(SET_COLS)); d$Comparison <- as_comp(d$Comparison)
  s <- do.call(rbind, lapply(split(d, d$Comparison, drop = TRUE), function(x) {
    sp <- x[x$Set != ORIGINAL_LABEL, ]
    w <- sum(sp$n_PAS * sp$prior_var_used) / sum(sp$n_PAS)
    w0 <- x$prior_var_used[x$Set == ORIGINAL_LABEL]
    data.frame(Comparison = x$Comparison[1], W = w, top = max(x$prior_var_used),
               pct = if (length(w0)) 100 * (1 - w / w0) else NA_real_)
  }))
  s$x <- as.numeric(s$Comparison)
  s$lab <- ifelse(is.finite(s$pct), sprintf("W %.2f\n%.0f%% %s", s$W, abs(s$pct), ifelse(s$pct >= 0, "narrower", "wider")),
                  sprintf("W %.2f", s$W))
  ggplot(d, aes(Comparison, prior_var_used, fill = Set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.74) +
    geom_hline(yintercept = DESEQ2_PRIOR_VAR_FLOOR, linetype = "22", colour = "grey45", linewidth = 0.35) +
    geom_segment(data = s, aes(x = x - 0.42, xend = x + 0.42, y = W, yend = W),
                 inherit.aes = FALSE, colour = "grey10", linewidth = 0.9) +
    geom_text(data = s, aes(x = x, y = top, label = lab), inherit.aes = FALSE, vjust = -0.3, size = 2.0,
              lineheight = 0.9, fontface = "bold", colour = "grey10") +
    set_bars_scale() +
    scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    labs(x = NULL, y = "Dispersion prior variance\n(used in MAP shrinkage)", tag = tag,
         title = "Width of the dispersion prior",
         subtitle = "Estimated latent variance tau^2. Black line: Leading Edge and Remainder pooled (n-weighted). Dashed: DESeq2's 0.25 floor (reference analysis).") +
    fig_theme() + theme(panel.grid.major.x = element_blank())
}

# ---- outputs ------------------------------------------------------------------

# Purpose: Save dispersion outputs.
save_dispersion_outputs <- function() {
  out <- summary_table_dir
  pas <- if (length(dispersion_pas_rows)) dplyr::bind_rows(dispersion_pas_rows) else NULL
  trd <- if (length(dispersion_trend_rows)) dplyr::bind_rows(dispersion_trend_rows) else NULL
  met <- if (length(dispersion_metrics_rows)) dplyr::bind_rows(dispersion_metrics_rows) else NULL
  comp <- if (length(composition_rows)) dplyr::bind_rows(composition_rows) else NULL
  if (!is.null(met)) save_csv(met, file.path(out, "Table_Dispersion_After_EVS.csv"))
  if (length(empirical_null_rows)) {
    save_csv(dplyr::bind_rows(empirical_null_rows), file.path(out, "Table_Empirical_Null_By_Analysis.csv"))
  }
  if (!is.null(comp)) save_csv(comp, file.path(out, "Table_EVS_Composition.csv"))
  if (!is.null(trd)) save_csv(trd, file.path(out, "Table_Dispersion_Trends_By_View.csv"))
  if (!is.null(pas)) save_csv(pas, file.path(out, "Table_Dispersion_Estimates_By_View.csv"))
  if (length(active_cutoff_rows)) save_csv(dplyr::bind_rows(active_cutoff_rows), file.path(out, "Table_Active_Cutoff.csv"))
  save_cutoff_outputs()

  if (!is.null(pas) && !is.null(trd) && !is.null(met)) try_fig({
    fig2 <- assemble_rows(list(panel_dispersion_grid(pas, trd, met, "a"),
                               list(plots = list(panel_trend_exponent(met, "b"), panel_prior_width(met, "c")),
                                    widths = c(1, 1))),
                          heights = c(2.2, 0.9),
                          title = paste0("Figure 2. DESeq2 dispersions, trend shape and prior width, Original vs split (",
                                         CUTOFF_METHOD_LABEL, ")"))
    save_figure(fig2, paper_fig_dir, "Figure_2_Dispersion", 180, 215)
  }, "Figure 2")
  invisible(TRUE)
}

# =============================================================================
# EVS CUTOFF k*: HELD-OUT SPLIT LIKELIHOOD
# =============================================================================
# The held-out likelihood responses are derived from six non-tested arms.
# VST, eligibility and PC1 ranks still use tested samples: this is not external
# validation. Also scored on the same held-out likelihood: one pooled trend, a
# continuous rank effect, splits ranked by within-arm variance and by the
# PC1-PC2 subspace score, an independent reference ranking built from the
# held-out samples only, and permutation searches.
EXPERIMENT_ARMS <- c("R0", "ZT6", "R2", "ZT8", "R4", "ZT10", "R8", "ZT14")
EXPERIMENT_CACHE <- new.env(parent = emptyenv())
SPLIT_MIN_SET <- 200L
SPLIT_GRID_STEPS <- c(250L, 25L, 5L, 1L)
# Permutation searches (whole ranking permuted, searched on the coarse grid at
# the same resolution as the observed coarse-grid maximum). Descriptive.
SPLIT_RANDOM_SPLITS <- as.integer(Sys.getenv("SEQUENCE_RANDOM_SPLITS", "5"))
stopifnot(!is.na(SPLIT_RANDOM_SPLITS), SPLIT_RANDOM_SPLITS >= 0L)
SPLIT_CORES <- as.integer(Sys.getenv("SEQUENCE_CORES",
  as.character(max(1L, min(4L, parallel::detectCores(logical = FALSE) %||% 1L)))))
if (is.na(SPLIT_CORES) || SPLIT_CORES < 1L) SPLIT_CORES <- 1L
SPLIT_SEED <- 20261125L
CUTOFF <- NULL
ARM_RANKS <- NULL
ARM_COMPARATOR_POS <- NULL
REFERENCE_K <- list()
# FALSE skips the comparators, the reference ranking and the permutations
# (used by the calibration script, which needs k* only).
CUTOFF_COMPARATORS <- TRUE

# Experiment-wide raw counts: every PAS with a nonzero count in any sample.
experiment_counts <- function() {
  if (!is.null(EXPERIMENT_CACHE$raw)) return(EXPERIMENT_CACHE$raw)
  x <- as.matrix(WTTS_Seq[, meta_all$id, drop = FALSE])
  storage.mode(x) <- "numeric"
  rownames(x) <- rownames(WTTS_Seq)
  x[!is.finite(x)] <- 0
  x <- x[rowSums(x) > 0, , drop = FALSE]
  EXPERIMENT_CACHE$raw <- coerce_raw_count_matrix_for_deseq2(x, context = "experiment-wide count matrix")
  EXPERIMENT_CACHE$raw
}

# Purpose: Experiment arm.
experiment_arm <- function(samples) {
  g <- rep(NA_character_, length(samples))
  for (a in EXPERIMENT_ARMS) {
    hit <- grepl(paste0("^", a, "_"), samples)
    if (any(!is.na(g[hit]))) stop("A sample matched more than one arm: ", a)
    g[hit] <- a
  }
  if (anyNA(g)) stop("Samples without an arm: ", paste(samples[is.na(g)], collapse = ", "))
  factor(g, levels = EXPERIMENT_ARMS)
}

# Experiment-wide DESeq2 object with median-of-ratios size factors (used for
# the VST and the held-out dispersions only; every DE test uses the
# comparison's own shared size factors).
experiment_dds <- function() {
  if (!is.null(EXPERIMENT_CACHE$dds)) return(EXPERIMENT_CACHE$dds)
  raw <- experiment_counts()
  dds <- DESeq2::DESeqDataSetFromMatrix(raw, colData = data.frame(arm = experiment_arm(colnames(raw)),
                                                                   row.names = colnames(raw)), design = ~ arm)
  dds <- tryCatch(DESeq2::estimateSizeFactors(dds),
                  error = function(e) DESeq2::estimateSizeFactors(dds, type = "poscounts"))
  EXPERIMENT_CACHE$dds <- dds
  dds
}

# The experiment-wide PCA input: DESeq2 VST (blind = TRUE), cached.
experiment_rank_matrix <- function(method = EVS_PCA_METHOD) {
  method <- match.arg(method, PCA_METHOD_LEVELS)
  if (!is.null(EXPERIMENT_CACHE$vst)) return(EXPERIMENT_CACHE$vst)
  raw <- experiment_counts()
  x <- as.matrix(SummarizedExperiment::assay(
    suppressMessages(DESeq2::varianceStabilizingTransformation(experiment_dds(), blind = TRUE))))
  dimnames(x) <- dimnames(raw)
  if (any(!is.finite(x))) stop("The experiment-wide VST matrix has non-finite values.")
  EXPERIMENT_CACHE$vst <- x
  x
}

# One arm's PAS ranks by |PC1 loading| (rank 1 = smallest, rank N = largest).
# Ties are broken exactly as in the split (order(|loading|, decreasing = TRUE)
# keeps tied PASs in row order), so the top k here is the split's top k.
arm_pc1_ranks <- function(rank_matrix, arm_samples) {
  z <- rank_matrix[, arm_samples, drop = FALSE]
  loading <- stats::prcomp(t(z), center = TRUE, scale. = FALSE, rank. = 1)$rotation[, 1]
  loading[!is.finite(loading)] <- 0
  ord <- rev(order(abs(loading), decreasing = TRUE))
  stats::setNames(match(seq_along(loading), ord), rownames(rank_matrix))
}

# Position from the top (1 = largest |PC1 loading|) of each comparison PAS in
# its two arms, and the union depth: a PAS is in the Leading Edge at k exactly
# when its union depth is <= k.
union_depths <- function(comparison_name, ids) {
  i <- match(comparison_name, comparison_table$comparison_name)
  N <- length(ARM_RANKS[[1]])
  pos <- function(arm) N - as.integer(ARM_RANKS[[arm]][ids]) + 1L
  rt <- pos(comparison_table$group1_prefix[i]); zt <- pos(comparison_table$group2_prefix[i])
  if (anyNA(rt) || anyNA(zt)) stop(comparison_name, ": PASs missing from the experiment-wide ranking")
  list(RT = stats::setNames(rt, ids), ZT = stats::setNames(zt, ids), depth = stats::setNames(pmin(rt, zt), ids))
}

# Comparison inputs computed once and reused by the analysis: raw counts,
# shared size factors, gene-wise dispersions and the union depths.
prepare_comparison_cache <- function(comparison_name, count_matrix, coldata) {
  raw <- coerce_raw_count_matrix_for_deseq2(count_matrix, context = paste0(comparison_name, " comparison matrix"))
  sf <- comparison_size_factors(raw, coldata)
  gdds <- comparison_gene_dds(raw, coldata, sf)
  dp <- union_depths(comparison_name, rownames(raw))
  comparison_cache[[comparison_name]] <<- list(raw = raw, coldata = coldata, depth = dp$depth, rank_RT = dp$RT,
                                               rank_ZT = dp$ZT, sf = sf, gene_dds = gdds)
  invisible(TRUE)
}

# ---- Profile marginal adjusted likelihood over the EVS cutoff --------------
# Each side is fitted with the same likelihood and bounds used in the final
# moderation. Fits are warm-started from the nearest evaluated cutoff and fall
# back to the pooled start (with multiple starts) if the warm start fails.
split_loglik <- function(k, depth, lx, likelihood, start_le, start_rem, fallback) {
  le <- depth <= k
  a <- fit_marginal_trend(lx[le], likelihood[le,,drop=FALSE], start_le, fallback = fallback)
  b <- fit_marginal_trend(lx[!le], likelihood[!le,,drop=FALSE], start_rem, fallback = fallback)
  list(k=k, loglik=a$loglik+b$loglik, le=a, rem=b, converged=a$converged && b$converged)
}

# Build independent-response likelihood curves using only the other six arms.
# Their size factors are estimated on those arms alone (not on target counts).
heldout_dispersions <- function(comparison_name, tested_arms) {
  raw <- experiment_counts(); arm <- experiment_arm(colnames(raw))
  ho <- colnames(raw)[!as.character(arm) %in% tested_arms]
  cd <- data.frame(arm=droplevels(arm[match(ho,colnames(raw))]), row.names=ho)
  dds <- DESeq2::DESeqDataSetFromMatrix(raw[,ho,drop=FALSE],colData=cd,design=~arm)
  dds <- apply_size_factors(dds, comparison_size_factors(raw[,ho,drop=FALSE], cd, formula=~arm))
  dds <- suppressMessages(DESeq2::estimateDispersionsGeneEst(dds,quiet=TRUE))
  mc <- S4Vectors::mcols(dds); lk <- make_dispersion_likelihood(dds)
  list(samples=ho, dds=dds, mean=stats::setNames(mc$baseMean,rownames(dds)),
       disp=stats::setNames(mc$dispGeneEst,rownames(dds)),
       use=stats::setNames(rownames(dds) %in% rownames(lk),rownames(dds)), likelihood=lk)
}

# Grid-refined maximum of the split likelihood over k. The first stage
# evaluates the admissible endpoints and a grid of step steps[1]; each later
# stage refines around the three best candidates. Ties prefer smaller k. This
# is a grid-refined maximum, not a proven global maximum.
search_cutoff <- function(depth, lx, likelihood, one, k_max, depth_all=depth, steps=SPLIT_GRID_STEPS) {
  valid <- function(k) k>=1L && k<=k_max && sum(depth<=k)>=SPLIT_MIN_SET &&
    sum(depth>k)>=SPLIT_MIN_SET && sum(depth_all<=k)>=SPLIT_MIN_SET && sum(depth_all>k)>=SPLIT_MIN_SET
  prof <- list()
  nearest <- function(k) {
    if (!length(prof)) return(NULL)
    ks <- as.integer(names(prof)); ok <- vapply(prof, function(p) isTRUE(p$converged), logical(1))
    if (!any(ok)) return(NULL)
    prof[[as.character(ks[ok][which.min(abs(ks[ok]-k))])]]
  }
  eval_k <- function(ks) for (k in sort(unique(as.integer(ks)))) if (valid(k) && is.null(prof[[as.character(k)]])) {
    nb <- nearest(k)
    prof[[as.character(k)]] <<- split_loglik(k, depth, lx, likelihood,
                                             if (is.null(nb)) one$par else nb$le$par,
                                             if (is.null(nb)) one$par else nb$rem$par, one$par)
  }
  lower <- max(sort(depth)[SPLIT_MIN_SET], sort(depth_all)[SPLIT_MIN_SET])
  upper <- min(k_max, sort(depth)[length(depth)-SPLIT_MIN_SET+1L]-1L,
               sort(depth_all)[length(depth_all)-SPLIT_MIN_SET+1L]-1L)
  if (lower>upper) stop("No admissible two-partition cutoff.")
  coarse_k <- unique(c(lower, seq(lower, upper, by=steps[1]), upper))
  eval_k(coarse_k)
  top <- function(n=3L) {
    v <- vapply(prof, `[[`, numeric(1), "loglik"); k <- as.integer(names(v)); ok <- is.finite(v)
    head(k[ok][order(-v[ok],k[ok])],n)
  }
  if (!length(top())) stop("No converged candidate cutoff.")
  for (j in seq_along(steps)[-1]) {
    w <- steps[j-1L]; centers <- top()
    eval_k(unlist(lapply(centers,function(k) seq(k-w,k+w,by=steps[j]))))
  }
  k <- top(1L)
  coarse_ll <- vapply(prof[as.character(coarse_k[coarse_k %in% as.integer(names(prof))])], `[[`, numeric(1), "loglik")
  list(k=k, best=prof[[as.character(k)]], prof=prof, coarse_max_loglik=max(coarse_ll[is.finite(coarse_ll)]))
}

# Position from the top (1 = largest score) of every PAS for a per-arm score.
score_positions <- function(score) stats::setNames(rank(-score, ties.method = "first"), names(score))

# Comparator rankings of one arm (computed once on the experiment-wide VST):
#   pc12      lambda1 v1^2 + lambda2 v2^2 (weighted contribution to the PC1-PC2
#             subspace; all nonzero components would reduce it to the variance)
#   variance  within-arm VST variance
arm_comparator_positions <- function(rank_matrix, arm_samples) {
  z <- rank_matrix[, arm_samples, drop = FALSE]
  pc <- stats::prcomp(t(z), center = TRUE, scale. = FALSE, rank. = 2)
  lam <- pc$sdev[1:2]^2
  s12 <- lam[1] * pc$rotation[, 1]^2 + lam[2] * pc$rotation[, 2]^2
  list(pc12 = score_positions(stats::setNames(s12, rownames(z))),
       variance = score_positions(stats::setNames(apply(z, 1L, stats::var), rownames(z))))
}

# Independent reference ranking (sensitivity analysis): the ranking is built
# only from the 30 held-out samples, so no tested count enters the split. The
# VST is computed on the held-out samples alone; within each held-out arm PASs
# are positioned by |PC1 loading|; a condition's reference position is the mean
# of its three arms' positions, re-ranked; the union depth is the smaller of
# the treatment and control reference positions.
reference_depths <- function(ho, ids) {
  arm <- experiment_arm(ho$samples)
  vst <- as.matrix(SummarizedExperiment::assay(suppressMessages(
    DESeq2::varianceStabilizingTransformation(ho$dds, blind = TRUE))))
  pos <- lapply(levels(droplevels(arm)), function(a) {
    l <- stats::prcomp(t(vst[, ho$samples[arm == a], drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 1)$rotation[, 1]
    score_positions(stats::setNames(abs(l), rownames(vst)))
  })
  names(pos) <- levels(droplevels(arm))
  cond_pos <- function(prefix) {
    arms <- names(pos)[grepl(prefix, names(pos))]
    m <- rowMeans(do.call(cbind, lapply(pos[arms], function(p) p[rownames(vst)])))
    score_positions(-m)              # smallest mean position = top
  }
  rt <- cond_pos("^R"); zt <- cond_pos("^ZT")
  d <- pmin(rt[ids], zt[ids])
  d[is.na(d)] <- max(d, na.rm = TRUE) + 1L   # PASs absent from the held-out samples: never Leading Edge
  stats::setNames(as.integer(d), ids)
}

# The cutoff and all comparators of one comparison (run in parallel).
cutoff_one <- function(i) {
  set.seed(SPLIT_SEED + i)
  raw <- experiment_counts(); arm <- experiment_arm(colnames(raw))
  cmp <- comparison_table$comparison_name[i]; g <- c(comparison_table$group1_prefix[i],comparison_table$group2_prefix[i])
  tested <- colnames(raw)[as.character(arm)%in%g]
  ids <- rownames(raw)[eligible_pas(raw[,tested,drop=FALSE])]
  depth_all <- union_depths(cmp,ids)$depth; ho <- heldout_dispersions(cmp,g)
  fit_ids <- ids[ho$use[ids]]; depth <- depth_all[fit_ids]; lx <- log(ho$mean[fit_ids]); lk <- ho$likelihood[fit_ids,,drop=FALSE]
  n <- length(lx); if (n < 2*SPLIT_MIN_SET) stop(cmp,": too few held-out PASs.")
  one <- fit_marginal_trend(lx,lk,multistart=TRUE)
  if (!one$converged) stop(cmp,": pooled model failed to converge.")
  l_one <- one$loglik; k_max <- floor(length(ids)/2)
  sr <- search_cutoff(depth,lx,lk,one,k_max,depth_all); k <- sr$k; best <- sr$best; prof <- sr$prof
  lrt <- 2*(best$loglik-l_one); n_le <- sum(depth<=k)
  bic_one <- -2*l_one+4*log(n); bic_split <- -2*best$loglik+9*log(n)
  # Continuous comparator: a smooth rank effect (natural cubic spline, 3 df)
  # in the prior centre instead of a binary split.
  bic_cont <- NA_real_; vs <- ps <- rs <- NULL; ref_depth <- stats::setNames(rep(NA_integer_, length(ids)), ids)
  if (CUTOFF_COMPARATORS) {
  basis <- splines::ns(log1p(depth),df=3); basis <- scale(basis,center=TRUE,scale=FALSE)
  continuous <- fit_marginal_trend(lx,lk,one$par,basis=basis,fallback=one$par)
  bic_cont <- if(continuous$converged) -2*continuous$loglik+7*log(n) else NA_real_
  # Binary comparators with the same search (stages 250 and 25).
  comp_search <- function(dp) tryCatch(search_cutoff(dp[fit_ids],lx,lk,one,k_max,dp,steps=SPLIT_GRID_STEPS[1:2]),
                                      error=function(e) NULL)
  comp_depth <- function(key) { d <- pmin(ARM_COMPARATOR_POS[[g[1]]][[key]][ids], ARM_COMPARATOR_POS[[g[2]]][[key]][ids]); names(d) <- ids; d }
  vs <- comp_search(comp_depth("variance")); ps <- comp_search(comp_depth("pc12"))
  # Independent reference ranking (held-out samples only), same full search.
  ref_depth <- reference_depths(ho, ids)
  rs <- tryCatch(search_cutoff(ref_depth[fit_ids],lx,lk,one,k_max,ref_depth),error=function(e) NULL)
  }
  bic_var <- if(is.null(vs)) NA_real_ else -2*vs$best$loglik+9*log(n)
  bic_pc12 <- if(is.null(ps)) NA_real_ else -2*ps$best$loglik+9*log(n)
  bic_ref <- if(is.null(rs)) NA_real_ else -2*rs$best$loglik+9*log(n)
  le_main <- ids[depth_all<=k]; le_ref <- if(is.null(rs)) character() else ids[ref_depth<=rs$k]
  # Permutations at matched resolution: the entire eligible ranking is
  # permuted and searched on the same coarse grid (step SPLIT_GRID_STEPS[1]);
  # compared with the observed coarse-grid maximum. Descriptive only.
  rnd <- vapply(seq_len(if (CUTOFF_COMPARATORS) SPLIT_RANDOM_SPLITS else 0L),function(r) {
    dp <- stats::setNames(sample(unname(depth_all)),ids)
    z <- tryCatch(search_cutoff(dp[fit_ids],lx,lk,one,k_max,dp,steps=SPLIT_GRID_STEPS[1]),error=function(e) NULL)
    if(is.null(z)) NA_real_ else 2*(z$best$loglik-l_one)
  },numeric(1))
  valid_rnd <- rnd[is.finite(rnd)]; coarse_lrt <- 2*(sr$coarse_max_loglik-l_one)
  pr <- data.frame(Comparison=cmp,k=as.integer(names(prof)),LogLik=vapply(prof,`[[`,numeric(1),"loglik"),
    Leading_Edge_fitted=vapply(as.integer(names(prof)),function(kk) sum(depth<=kk),integer(1)),
    Converged=vapply(prof,function(p) isTRUE(p$converged),logical(1)))
  pr <- pr[order(pr$k),,drop=FALSE]; pr$LogLik[!pr$Converged] <- NA_real_
  pr$LRT_vs_one_trend <- 2*(pr$LogLik-l_one); pr$Delta_from_max <- 2*(pr$LogLik-best$loglik); pr$Chosen <- pr$k==k
  drop_at <- function(d) {j<-which.min(abs(pr$k-(k+d))); if(abs(pr$k[j]-(k+d))<=125) pr$Delta_from_max[j] else NA_real_}
  q <- function(p) c(a0=exp(p[1]),a1=exp(p[2]),b=p[3])
  pl <- q(best$le$par); pm <- q(best$rem$par); po <- q(one$par)
  bics <- c(pooled=bic_one, evs=bic_split, continuous=bic_cont, variance=bic_var, pc12=bic_pc12)
  bics <- bics[is.finite(bics)]
  summary <- data.frame(Comparison=cmp,EVS_PCA_Method=EVS_PCA_METHOD,k=k,
    Leading_Edge=sum(depth_all<=k),Remainder=sum(depth_all>k),Heldout_Samples=paste(ho$samples,collapse=","),
    n_PAS_fitted=n,LE_fitted=n_le,REM_fitted=n-n_le,LogLik_split=best$loglik,LogLik_one_trend=l_one,LRT_vs_one_trend=lrt,
    BIC_pooled=bic_one,BIC_EVS=bic_split,BIC_continuous_rank=bic_cont,BIC_variance_split=bic_var,BIC_pc12_split=bic_pc12,
    Variance_k=if(is.null(vs)) NA_integer_ else vs$k, PC12_k=if(is.null(ps)) NA_integer_ else ps$k,
    Best_Model_By_BIC=names(bics)[which.min(bics)],
    Comparator_Better_Than_EVS=any(c(bic_cont,bic_var,bic_pc12)<bic_split,na.rm=TRUE),
    Selected_Model=if(bic_split<bic_one) "EVS" else "Pooled",
    Reference_k=if(is.null(rs)) NA_integer_ else rs$k, Reference_Leading_Edge=length(le_ref),
    BIC_reference_ranking=bic_ref,
    Reference_LE_Jaccard=if(length(le_ref)) length(intersect(le_main,le_ref))/length(union(le_main,le_ref)) else NA_real_,
    Coarse_grid_LRT_observed=coarse_lrt,
    Permutation_searches=SPLIT_RANDOM_SPLITS, Permutation_searches_converged=length(valid_rnd),
    Random_Split_LRT_Mean=if(length(valid_rnd)) mean(valid_rnd) else NA_real_,
    Random_Split_LRT_Max=if(length(valid_rnd)) max(valid_rnd) else NA_real_,
    Permutations_at_or_above_observed=if(length(valid_rnd)) sum(valid_rnd>=coarse_lrt) else NA_integer_,
    Delta_at_minus500=drop_at(-500),Delta_at_plus500=drop_at(500),
    LE_a0=pl[["a0"]],LE_a1=pl[["a1"]],LE_b=pl[["b"]],LE_resid_sd=sqrt(best$le$tau2)/log(10),
    REM_a0=pm[["a0"]],REM_a1=pm[["a1"]],REM_b=pm[["b"]],REM_resid_sd=sqrt(best$rem$tau2)/log(10),
    One_a0=po[["a0"]],One_a1=po[["a1"]],One_b=po[["b"]],
    LE_prior_var=best$le$tau2,REM_prior_var=best$rem$tau2,One_prior_var=one$tau2,
    LE_prior_var_at_bound=best$le$tau2<=DISP_PRIOR_VAR_MIN*1.00001, REM_prior_var_at_bound=best$rem$tau2<=DISP_PRIOR_VAR_MIN*1.00001,
    Converged=isTRUE(best$converged) && isTRUE(one$converged))
  # Numerical checks at k* (run last; in a parallel child the global changes
  # below affect only this comparison's process, and they are restored):
  # optimizer restarts, a finer dispersion grid (step / 2), and a lower tau^2
  # bound (bound / 10). Reported as changes in log-likelihood and parameters.
  num <- if (!CUTOFF_COMPARATORS) NULL else tryCatch({
    le_k <- depth <= k
    ms <- function(lkm, lxv) list(le = fit_marginal_trend(lxv[le_k], lkm[le_k,,drop=FALSE], multistart = TRUE),
                                  rem = fit_marginal_trend(lxv[!le_k], lkm[!le_k,,drop=FALSE], multistart = TRUE))
    r0 <- ms(lk, lx)
    old <- list(g = DISP_LOG_GRID, e = DISP_LOG_EDGES, m = DISP_PRIOR_VAR_MIN)
    on.exit({ DISP_LOG_GRID <<- old$g; DISP_LOG_EDGES <<- old$e; DISP_PRIOR_VAR_MIN <<- old$m }, add = TRUE)
    DISP_LOG_GRID <<- seq(min(old$g), max(old$g), by = DISP_GRID_STEP / 2)
    DISP_LOG_EDGES <<- c(-Inf, head(DISP_LOG_GRID, -1) + diff(DISP_LOG_GRID) / 2, Inf)
    lk2 <- make_dispersion_likelihood(ho$dds)[fit_ids, , drop = FALSE]
    r1 <- ms(lk2, lx)
    DISP_LOG_GRID <<- old$g; DISP_LOG_EDGES <<- old$e; DISP_PRIOR_VAR_MIN <<- old$m / 10
    r2 <- ms(lk, lx)
    DISP_PRIOR_VAR_MIN <<- old$m
    pd <- function(a, b) max(abs(c(a$le$par[1:3] - b$le$par[1:3], a$rem$par[1:3] - b$rem$par[1:3])))
    tau <- function(a) c(a$le$tau2, a$rem$tau2)
    data.frame(Check_multistart_dLogLik = (r0$le$loglik + r0$rem$loglik) - best$loglik,
               Check_finer_grid_max_par_change = pd(r1, r0),
               Check_finer_grid_max_tau2_rel_change = max(abs(tau(r1) / tau(r0) - 1)),
               Check_lower_tau2_bound_max_par_change = pd(r2, r0),
               Check_lower_tau2_bound_tau2_LE = r2$le$tau2, Check_lower_tau2_bound_tau2_REM = r2$rem$tau2)
  }, error = function(e) data.frame(Check_multistart_dLogLik = NA_real_, Check_finer_grid_max_par_change = NA_real_,
                                    Check_finer_grid_max_tau2_rel_change = NA_real_, Check_lower_tau2_bound_max_par_change = NA_real_,
                                    Check_lower_tau2_bound_tau2_LE = NA_real_, Check_lower_tau2_bound_tau2_REM = NA_real_))
  if (!is.null(num)) summary <- cbind(summary, num)
  pas <- data.frame(Comparison=cmp,feature_id=ids,Union_Depth=unname(depth_all),Reference_Depth=unname(ref_depth[ids]),
    Heldout_Mean=unname(ho$mean[ids]),Heldout_Dispersion=unname(ho$disp[ids]),In_Trend_Fit=unname(ho$use[ids]),
    In_Leading_Edge=unname(depth_all<=k))
  message(sprintf("%s: k* = %d (LE %d); LRT %.0f; selected %s; best BIC %s; reference k %s",
                  cmp,k,summary$Leading_Edge,lrt,summary$Selected_Model,summary$Best_Model_By_BIC,
                  if(is.null(rs)) "NA" else rs$k))
  list(cmp=cmp, k=k, summary=summary, profile=pr,
       random=data.frame(Comparison=rep(cmp,length(rnd)),Search=seq_along(rnd),Coarse_LRT_vs_one_trend=rnd),
       pas=pas, ref_depth=ref_depth, ref_k=if(is.null(rs)) NA_integer_ else rs$k)
}

# Lock k* (and the comparators) for every comparison before any DE test. The
# comparisons are independent and run in parallel (SEQUENCE_CORES, default up
# to 4); each uses its own fixed seed, so results do not depend on the cores.
compute_cutoffs <- function() {
  raw <- experiment_counts(); arm <- experiment_arm(colnames(raw)); rm <- experiment_rank_matrix()
  ARM_RANKS <<- lapply(stats::setNames(EXPERIMENT_ARMS,EXPERIMENT_ARMS),function(a) arm_pc1_ranks(rm,colnames(raw)[arm==a]))
  ARM_COMPARATOR_POS <<- lapply(stats::setNames(EXPERIMENT_ARMS,EXPERIMENT_ARMS),function(a)
    arm_comparator_positions(rm,colnames(raw)[arm==a]))
  cores <- if (.Platform$OS.type == "unix") SPLIT_CORES else 1L
  res <- parallel::mclapply(seq_len(nrow(comparison_table)), function(i)
    tryCatch(cutoff_one(i), error = function(e) structure(conditionMessage(e), class = "cutoff_error")),
    mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) inherits(r, "cutoff_error") || inherits(r, "try-error") || is.null(r), logical(1))
  if (any(bad)) stop("Cutoff selection failed: ", paste(vapply(res[bad], function(r) paste(as.character(r), collapse=" "), ""), collapse = " | "))
  for (r in res) { OPTIMAL_K[[r$cmp]] <<- r$k; REFERENCE_K[[r$cmp]] <<- r$ref_k }
  CUTOFF <<- list(summary=dplyr::bind_rows(lapply(res,`[[`,"summary")),
                  profiles=stats::setNames(lapply(res,`[[`,"profile"), vapply(res,`[[`,"","cmp")),
                  pas=stats::setNames(lapply(res,`[[`,"pas"), vapply(res,`[[`,"","cmp")),
                  ref_depth=stats::setNames(lapply(res,`[[`,"ref_depth"), vapply(res,`[[`,"","cmp")),
                  random=dplyr::bind_rows(lapply(res,`[[`,"random")))
  invisible(CUTOFF)
}

# ---- cutoff outputs: tables and Figure 1 --------------------------------------

# Purpose: Save cutoff outputs.
save_cutoff_outputs <- function() {
  if (is.null(CUTOFF)) return(invisible(FALSE))
  out <- summary_table_dir
  save_csv(CUTOFF$summary, file.path(out, "Table_Cutoff_Summary.csv"))
  save_csv(dplyr::bind_rows(CUTOFF$profiles), file.path(out, "Table_Cutoff_Profile_Likelihood.csv"))
  save_csv(CUTOFF$random, file.path(out, "Table_Cutoff_Random_Splits.csv"))
  save_csv(dplyr::bind_rows(CUTOFF$pas), file.path(out, "Table_Cutoff_PAS_Heldout_Dispersions.csv"))
  le <- dplyr::bind_rows(lapply(CUTOFF$pas, function(p) p[p$In_Leading_Edge, c("Comparison", "feature_id")]))
  save_csv(le, file.path(out, "Table_Cutoff_Leading_Edge_Sites.csv"))
  try_fig(save_figure(plot_cutoff_figure(), paper_fig_dir, "Figure_1_EVS_Cutoff", 180, 140), "Figure 1 (cutoff)")
  invisible(TRUE)
}

# Figure 1: (a) held-out dispersions with the Leading-Edge and Remainder trends
# at k* (and the single trend), (b) the profile likelihood of the split over k.
plot_cutoff_figure <- function() {
  summ <- CUTOFF$summary; summ$Comparison <- as_comp(summ$Comparison)
  pts <- dplyr::bind_rows(CUTOFF$pas)
  pts <- pts[pts$In_Trend_Fit %in% TRUE, , drop = FALSE]
  pts$Comparison <- as_comp(pts$Comparison)
  pts$Set <- factor(ifelse(pts$In_Leading_Edge, "Leading Edge", "Remainder"), levels = c("Leading Edge", "Remainder"))
  set.seed(SPLIT_SEED)
  pts <- do.call(rbind, lapply(split(pts, pts$Comparison, drop = TRUE), function(z)
    z[sort(sample.int(nrow(z), min(nrow(z), 12000L))), , drop = FALSE]))
  pts <- pts[order(pts$Set, decreasing = TRUE), , drop = FALSE]
  curves <- do.call(rbind, lapply(seq_len(nrow(summ)), function(i) {
    xr <- range(pts$Heldout_Mean[pts$Comparison == summ$Comparison[i]])
    xs <- exp(seq(log(max(xr[1], 0.1)), log(xr[2]), length.out = 200))
    tr <- function(a0, a1, b) a0 + a1 * xs^(-b)
    rbind(data.frame(Comparison = summ$Comparison[i], x = xs, Trend = "Leading Edge",
                     y = tr(summ$LE_a0[i], summ$LE_a1[i], summ$LE_b[i])),
          data.frame(Comparison = summ$Comparison[i], x = xs, Trend = "Remainder",
                     y = tr(summ$REM_a0[i], summ$REM_a1[i], summ$REM_b[i])),
          data.frame(Comparison = summ$Comparison[i], x = xs, Trend = "One trend (no split)",
                     y = tr(summ$One_a0[i], summ$One_a1[i], summ$One_b[i])))
  }))
  lab <- data.frame(Comparison = summ$Comparison,
                    lab = sprintf("k* = %s; LE %s PASs\nLE %.3f + %.1f m^-%.2f\nREM %.3f + %.1f m^-%.2f",
                                  format(summ$k, big.mark = ","), format(summ$Leading_Edge, big.mark = ","),
                                  summ$LE_a0, summ$LE_a1, summ$LE_b, summ$REM_a0, summ$REM_a1, summ$REM_b))
  tcols <- c("Leading Edge" = SET_COLS[["Leading Edge"]], "Remainder" = SET_COLS[["Remainder"]],
             "One trend (no split)" = "grey20")
  pa <- ggplot(pts, aes(Heldout_Mean, Heldout_Dispersion)) +
    geom_point(aes(colour = Set), size = 0.25, alpha = 0.35, stroke = 0, show.legend = FALSE) +
    geom_line(data = curves, aes(x, y, colour = Trend, linetype = Trend), linewidth = 0.6) +
    geom_text(data = lab, aes(label = lab), x = Inf, y = Inf, hjust = 1.03, vjust = 1.1, size = 1.75,
              lineheight = 0.9, inherit.aes = FALSE) +
    facet_wrap(~ Comparison, nrow = 1, labeller = labeller(Comparison = comp_lab())) +
    scale_x_log10(labels = scales::label_comma(accuracy = 1)) +
    scale_y_log10(labels = function(v) format(v, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
    scale_colour_manual(values = tcols, breaks = names(tcols), name = NULL) +
    scale_linetype_manual(values = c("Leading Edge" = "solid", "Remainder" = "solid", "One trend (no split)" = "22"),
                          breaks = names(tcols), name = NULL) +
    labs(x = "Mean, held-out samples", y = "Gene-wise dispersion,\nheld-out samples", tag = "a",
         title = "Dispersion trends of the Leading Edge and the Remainder at k* (30 samples not in the comparison)",
         subtitle = "Points: PASs coloured by set at k* (up to 12,000 per panel). Lines: a0 + a1 x mean^-b fitted to each set; dashed: one trend for all PASs.") +
    fig_theme()
  prof <- dplyr::bind_rows(CUTOFF$profiles)
  k_of <- stats::setNames(CUTOFF$summary$k, CUTOFF$summary$Comparison)
  # Retain all evaluated grid points, including shifted coarse-grid points.
  prof$Comparison <- as_comp(prof$Comparison)
  rnd <- data.frame(Comparison = summ$Comparison, y = summ$Random_Split_LRT_Max)
  kl <- data.frame(Comparison = summ$Comparison, k = summ$k, y = summ$LRT_vs_one_trend,
                   lab = sprintf("k* = %s\nLRT = %s", format(summ$k, big.mark = ","),
                                 format(round(summ$LRT_vs_one_trend), big.mark = ",")))
  pb <- ggplot(prof, aes(k, LRT_vs_one_trend)) +
    geom_hline(data = rnd, aes(yintercept = y), colour = "grey45", linetype = "22", linewidth = 0.35) +
    geom_line(linewidth = 0.45, colour = "grey15", na.rm=TRUE) +
    geom_vline(data = kl, aes(xintercept = k), colour = "#D55E00", linewidth = 0.4, linetype = "22") +
    geom_point(data = kl, aes(k, y), colour = "#D55E00", size = 1.4) +
    geom_text(data = kl, aes(k, y, label = lab), hjust = -0.08, vjust = 1.2, size = 1.9, colour = "#D55E00",
              fontface = "bold", lineheight = 0.9) +
    facet_wrap(~ Comparison, nrow = 1, scales = "free_y", labeller = labeller(Comparison = comp_lab())) +
    scale_x_continuous(labels = scales::label_comma()) + scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "k (top-|PC1| PASs per arm)", y = "LRT, split vs one trend", tag = "b",
         title = "Profile adjusted likelihood: grid-refined candidate cutoff",
         subtitle = "LRT = 2 x [marginal adjusted log-likelihood difference]; dashed horizontal: largest coarse-grid LRT of permutation searches (same resolution).") +
    fig_theme()
  assemble_rows(list(pa, pb), heights = c(1.1, 1),
                title = "Figure 1. EVS candidate cutoff: grid-refined marginal adjusted likelihood (VST ranking)")
}

# =============================================================================
# VALIDATION CHECKS AND DISCOVERY COUNTS
# =============================================================================

NORMALIZATION_MODE <- "shared"   # size factors from the full comparison matrix, in every view

# Purpose: Add validation.
add_validation <- function(check, comparison, passed, detail = "") {
  validation_rows[[length(validation_rows) + 1L]] <<- data.frame(
    Check = check, Comparison = comparison, Cutoff_Method = CUTOFF_METHOD_LABEL, EVS_PCA_Method = EVS_PCA_METHOD,
    Status = if (isTRUE(passed)) "PASS" else if (is.na(passed)) "NOT APPLICABLE" else "FAIL",
    Detail = detail)
  if (identical(passed, FALSE)) warning("Validation failed (", check, ", ", comparison, "): ", detail, call. = FALSE)
  invisible(passed)
}

# Purpose: Size factors match.
size_factors_match <- function(dds, size_factors) {
  if (is.null(size_factors)) return(NA)
  isTRUE(all.equal(unname(sizeFactors(dds)), unname(size_factors[colnames(dds)]), tolerance = 1e-12))
}

# Per-PAS calls of one analysis (a named list of view results): Standard,
# Strong and Weak from one BH over all its views' raw p-values (as in
# call_significance), HBFSS as calibrated per view.
discovery_flags <- function(res_list) {
  cols <- c("feature_id", "pvalue", "resGA_pvalue", "resLA_pvalue", "lfc_shrunk",
            "standard_flag", "strong_cnh_flag", "weak_significant_flag", "hbfss_flag")
  d <- do.call(rbind, lapply(names(res_list), function(v) {
    r <- res_list[[v]]
    miss <- setdiff(cols, names(r))
    if (length(miss)) stop("Result table ", v, " lacks: ", paste(miss, collapse = ", "))
    data.frame(View = v, r[, cols, drop = FALSE], stringsAsFactors = FALSE)
  }))
  if (anyDuplicated(d$feature_id)) stop("A PAS appears in more than one view of the same analysis.")
  bh <- function(p) {
    q <- rep(NA_real_, length(p))
    ok <- is.finite(p)
    q[ok] <- stats::p.adjust(p[ok], method = "BH")
    q
  }
  lfc_ok <- is.finite(d$lfc_shrunk)
  big <- lfc_ok & abs(d$lfc_shrunk) >= lfc_boundary
  small <- lfc_ok & abs(d$lfc_shrunk) < lfc_boundary
  d$padj_global <- bh(d$pvalue)
  d$padjGA_global <- bh(d$resGA_pvalue)
  d$padjLA_global <- bh(d$resLA_pvalue)
  d$HBFSS <- d$hbfss_flag %in% TRUE
  d$Standard <- !is.na(d$padj_global) & d$padj_global < BH_FDR_STANDARD
  d$Strong <- !is.na(d$padjGA_global) & d$padjGA_global < BH_FDR_STRONG
  d$Weak <- !is.na(d$padjLA_global) & d$padjLA_global < BH_FDR_WEAK
  d
}

# Purpose: Count discoveries.
count_discoveries <- function(flags) {
  c(Standard = sum(flags$Standard), Strong = sum(flags$Strong), Weak = sum(flags$Weak), HBFSS = sum(flags$HBFSS),
    n_tested = nrow(flags), n_with_pvalue = sum(is.finite(flags$pvalue)))
}

# Shared normalization reached every view, the samples match, the partition is
# exact, and a PAS has identical normalized counts in the
# Original and in its EVS partition.
validate_track_normalization <- function(comparison_name, track_key, analysis_results, sf, evs, coldata) {
  tr <- unname(track_short[track_key])
  ids <- rownames(evs$raw_dataset)
  add_validation("Every PAS is assigned exactly once (LE / REM)", comparison_name,
                 !length(intersect(evs$leading_edge_ids, evs$remainder_ids)) &&
                   setequal(c(evs$leading_edge_ids, evs$remainder_ids), ids) &&
                   length(evs$leading_edge_ids) + length(evs$remainder_ids) == length(ids), tr)
  same_samples <- all(vapply(analysis_results, function(a) identical(colnames(a$dds), rownames(coldata)), logical(1)))
  add_validation("Original and EVS sample identities match", comparison_name, same_samples, tr)
  ev <- c("leading_edge_dataset", "remainder_dataset")
  # The Original's own flagged PASs, read from its results (not from the cache),
  # so that a comparison with no outliers is checked as an empty list.
  r0 <- analysis_results[["raw_dataset"]]$results
  oc <- if (is.null(r0)) NULL else as.character(r0$feature_id[r0$cooks_outlier %in% TRUE])
  cooks_ok <- !is.null(oc) && all(vapply(ev, function(v) {
    r <- analysis_results[[v]]$results
    out <- r$cooks_outlier %in% TRUE
    setequal(r$feature_id[out], intersect(oc, r$feature_id)) && all(is.na(r$pvalue[out]))
  }, logical(1)))
  add_validation("Cook's outliers are the Original's in every view", comparison_name, cooks_ok,
                 sprintf("%s; %d outlier PASs in the Original", tr, length(oc)))
  ab <- do.call(rbind, lapply(ev, function(v) analysis_results[[v]]$results[, c("pvalue", "padj")]))
  q <- rep(NA_real_, nrow(ab))
  okp <- is.finite(ab$pvalue)
  q[okp] <- stats::p.adjust(ab$pvalue[okp], method = "BH")
  add_validation("Standard padj is one BH over LE + REM, without independent filtering", comparison_name,
                 isTRUE(all.equal(q, ab$padj)), tr)
  add_validation("Shared normalization is propagated", comparison_name,
                 all(vapply(analysis_results, function(a) isTRUE(size_factors_match(a$dds, sf)), logical(1))), tr)
  o0 <- analysis_results$raw_dataset %||% {
    reg0 <- paper_registry[[registry_key(comparison_name, "normalized_evs")]]
    if (is.null(reg0)) NULL else reg0$analysis_results$raw_dataset
  }
  if (!is.null(o0)) {
    og <- S4Vectors::mcols(o0$dds)$dispGeneEst
    names(og) <- rownames(o0$dds)
    same_gw <- all(vapply(c("leading_edge_dataset", "remainder_dataset"), function(v) {
      d <- analysis_results[[v]]$dds
      isTRUE(all.equal(unname(S4Vectors::mcols(d)$dispGeneEst), unname(og[rownames(d)]), tolerance = 1e-8))
    }, logical(1)))
    add_validation("Gene-wise dispersions are identical in the Original and the partition (computed once per comparison)",
                   comparison_name, same_gw, tr)
  }
  o <- analysis_results$raw_dataset %||% {
    reg <- paper_registry[[registry_key(comparison_name, "normalized_evs")]]
    if (is.null(reg)) NULL else reg$analysis_results$raw_dataset
  }
  if (!is.null(o)) {
    on <- DESeq2::counts(o$dds, normalized = TRUE)
    same <- all(vapply(c("leading_edge_dataset", "remainder_dataset"), function(v) {
      vn <- DESeq2::counts(analysis_results[[v]]$dds, normalized = TRUE)
      isTRUE(all.equal(unname(vn), unname(on[rownames(vn), colnames(vn), drop = FALSE]), tolerance = 1e-10))
    }, logical(1)))
    add_validation("Normalized counts are identical in the Original and the EVS partition", comparison_name, same, tr)
  }
  invisible(TRUE)
}

# The chosen k* is the one used: the split has k = k*, and its Leading Edge is
# exactly the set scored at k*.
validate_selected_k_used <- function(comparison_name) {
  k <- OPTIMAL_K[[comparison_name]]
  reg <- paper_registry[[registry_key(comparison_name, "normalized_evs")]]
  cache <- comparison_cache[[comparison_name]]
  srow <- if (is.null(CUTOFF)) NULL else CUTOFF$summary[CUTOFF$summary$Comparison == comparison_name, , drop = FALSE]
  pas <- if (is.null(CUTOFF)) NULL else CUTOFF$pas[[comparison_name]]
  sites <- if (is.null(pas)) NULL else pas$feature_id[pas$In_Leading_Edge]
  if (is.null(reg) || is.null(cache) || is.null(srow) || nrow(srow) != 1L) {
    add_validation("The selected k* is used throughout the analysis", comparison_name, FALSE, "missing cutoff fit or analysis")
    return(invisible(FALSE))
  }
  evs <- reg$evs
  ks <- c(selected = k, split = as.integer(evs$empirical_evs_k), active = get_active_evs_cutoff(comparison_name),
          table = as.integer(srow$k))
  add_validation("The selected k* is used throughout the analysis", comparison_name,
                 length(unique(ks)) == 1L, paste(names(ks), ks, sep = " = ", collapse = "; "))
  add_validation("The analysed Leading Edge is the set scored at k*", comparison_name,
                 setequal(evs$leading_edge_ids, sites) &&
                   setequal(sites, rownames(cache$raw)[cache$depth[rownames(cache$raw)] <= k]),
                 sprintf("Leading Edge %d; scored set %d", length(evs$leading_edge_ids), length(sites)))
  invisible(TRUE)
}

# PCA preprocessing checks for this run's method, plus the ranking table.
validate_pca_preprocessing <- function(comparison_name) {
  cc <- comparison_cache[[comparison_name]]
  raw <- cc$raw; coldata <- cc$coldata
  ro <- make_rank_matrix_for_track(raw, coldata, "normalized_evs", cache_key = comparison_name)
  x <- as.matrix(ro$rank_matrix)
  full <- experiment_rank_matrix(ro$method)
  add_validation("PCA input is finite, matches the raw PAS IDs and samples, and equals the experiment-wide matrix",
                 comparison_name,
                 all(is.finite(x)) && identical(rownames(x), rownames(raw)) && identical(colnames(x), colnames(raw)) &&
                   identical(unname(x), unname(full[rownames(raw), colnames(raw), drop = FALSE])),
                 sprintf("%s: %d PASs x %d samples; range %.3g to %.3g", ro$method, nrow(x), ncol(x),
                         min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
  fits <- pc1_condition_fits(ro, coldata, comparison_name)
  for (cond in c("trt", "untrt")) {
    v <- fits[[cond]]$pca_fit$sdev^2
    tot <- sum(v)
    ok <- is.finite(tot) && tot > 1e-12 && v[1] > 0 && (length(v) < 2L || v[1] > v[2] * (1 + 1e-8))
    add_validation("PCA has sufficient variance (PC1 nonzero and distinct from PC2)", comparison_name, ok,
                   sprintf("%s %s: PC1 %.1f%% of variance, PC1/PC2 variance ratio %.2f", ro$method,
                           if (cond == "trt") "RT" else "ZT", 100 * v[1] / tot,
                           if (length(v) > 1L) v[1] / v[2] else NA_real_))
  }
  rk <- lapply(fits, function(f) stats::setNames(f$loading_table$rank, f$loading_table$feature_id)[rownames(raw)])
  k <- OPTIMAL_K[[comparison_name]]
  same_top <- setequal(names(rk$trt)[rk$trt <= k], names(cc$rank_RT)[cc$rank_RT <= k]) &&
    setequal(names(rk$untrt)[rk$untrt <= k], names(cc$rank_ZT)[cc$rank_ZT <= k])
  add_validation("Comparison PC1 top-k* sets equal the experiment-wide arm top-k* sets", comparison_name, same_top,
                 sprintf("%s, k* = %d", ro$method, k))
  pca_ranking_rows[[comparison_name]] <<- data.frame(
    Comparison = comparison_name, EVS_PCA_Method = ro$method, feature_id = rownames(raw),
    Rank_RT = unname(rk$trt), Rank_ZT = unname(rk$untrt), Union_Depth = unname(cc$depth[rownames(raw)]))
  invisible(TRUE)
}

# Validation table, the DESeq2 default-trend baseline and the PC1 rankings
# (end of run).
save_run_checks <- function() {
  out <- summary_table_dir
  if (length(validation_rows)) save_csv(dplyr::bind_rows(validation_rows), file.path(out, "Table_Validation_Checks.csv"))
  if (length(default_trend_rows)) save_csv(dplyr::bind_rows(default_trend_rows), file.path(out, "Table_Original_DESeq2_Default_Trend.csv"))
  if (length(pca_ranking_rows)) save_csv(dplyr::bind_rows(pca_ranking_rows), file.path(out, "Table_EVS_PCA_Ranking.csv"))
  if (length(sensitivity_rows)) save_csv(dplyr::bind_rows(sensitivity_rows), file.path(out, "Table_Sensitivity_Analyses.csv"))
  if (length(normalization_rows)) save_csv(dplyr::bind_rows(normalization_rows), file.path(out, "Table_Normalization_Check.csv"))
  if (!is.null(CUTOFF)) try_fig(save_figure(plot_model_sensitivity_figure(), paper_fig_dir,
                                            "Figure_3_Model_Comparison_Sensitivity", 180, 170), "Figure 3")
  invisible(TRUE)
}


# Purpose: Run full comparison pipeline.
run_full_comparison_pipeline <- function(comparison_name, count_matrix, coldata, annot_df) {
  dir.create(file.path(output_dir, comparison_name), recursive = TRUE, showWarnings = FALSE)

  # Stage 1: k* was chosen by compute_cutoffs() (no test of any kind);
  # comparison inputs reused by the analysis are prepared here.
  calls_before <- c(SIGNIFICANCE_CALLS$wald, SIGNIFICANCE_CALLS$null)
  prepare_comparison_cache(comparison_name, count_matrix, coldata)
  calls_after <- c(SIGNIFICANCE_CALLS$wald, SIGNIFICANCE_CALLS$null)
  tryCatch(normalization_check(comparison_name, comparison_cache[[comparison_name]]$raw, comparison_cache[[comparison_name]]$sf),
           error = function(e) warning(comparison_name, ": normalization check failed: ", conditionMessage(e), call. = FALSE))
  add_validation("Comparison preparation performed no significance testing", comparison_name,
                 identical(calls_before, calls_after),
                 sprintf("Wald fits %d -> %d; significance calls %d -> %d",
                         calls_before[1], calls_after[1], calls_before[2], calls_after[2]))
  cs <- CUTOFF$summary[CUTOFF$summary$Comparison == comparison_name, , drop = FALSE]
  le_n <- sum(comparison_cache[[comparison_name]]$depth <= get_active_evs_cutoff(comparison_name))
  active_cutoff_rows[[comparison_name]] <<- data.frame(
    Comparison = comparison_name, Cutoff_Method = CUTOFF_METHOD_LABEL, k = get_active_evs_cutoff(comparison_name),
    Leading_Edge = le_n, Remainder = nrow(comparison_cache[[comparison_name]]$raw) - le_n,
    LRT_vs_one_trend = cs$LRT_vs_one_trend, Random_Split_LRT_Max = cs$Random_Split_LRT_Max, Selected_Model = cs$Selected_Model,
    Best_Model_By_BIC = cs$Best_Model_By_BIC, Reference_k = cs$Reference_k, Reference_LE_Jaccard = cs$Reference_LE_Jaccard,
    LE_prior_var = cs$LE_prior_var, REM_prior_var = cs$REM_prior_var,
    LE_a0 = cs$LE_a0, LE_a1 = cs$LE_a1, LE_b = cs$LE_b, REM_a0 = cs$REM_a0, REM_a1 = cs$REM_a1, REM_b = cs$REM_b,
    Trend_fits_converged = cs$Converged,
    EVS_PCA_Method = EVS_PCA_METHOD, PCA_Preprocessing = PCA_METHOD_LABELS[[EVS_PCA_METHOD]],
    Normalization_Mode = NORMALIZATION_MODE)
  ho <- strsplit(cs$Heldout_Samples, ",")[[1]]
  add_validation("Held-out likelihood response excludes tested samples (rank does not)", comparison_name,
                 !any(ho %in% colnames(comparison_cache[[comparison_name]]$raw)) && length(ho) == 30L,
                 sprintf("held-out samples: %s", cs$Heldout_Samples))
  add_validation("Dispersion trend fits at k* converged", comparison_name, isTRUE(cs$Converged),
                 sprintf("LE %.3f + %.2f x mean^-%.2f; REM %.3f + %.2f x mean^-%.2f",
                         cs$LE_a0, cs$LE_a1, cs$LE_b, cs$REM_a0, cs$REM_a1, cs$REM_b))
  add_validation("Observed coarse-grid LRT exceeds every permutation search at the same resolution (descriptive)",
                 comparison_name, if (cs$Permutation_searches_converged > 0) cs$Permutations_at_or_above_observed == 0 else NA,
                 sprintf("observed %.1f; permutations max %.1f (%d of %d converged)", cs$Coarse_grid_LRT_observed,
                         cs$Random_Split_LRT_Max, cs$Permutation_searches_converged, cs$Permutation_searches))
  add_validation("Cutoff fit: restarts from several starting values do not improve the likelihood (< 0.5)", comparison_name,
                 if (is.null(cs$Check_multistart_dLogLik) || !is.finite(cs$Check_multistart_dLogLik)) NA else cs$Check_multistart_dLogLik < 0.5,
                 sprintf("log-likelihood change %.3g", cs$Check_multistart_dLogLik %||% NA_real_))
  add_validation("Cutoff fit: halving the dispersion grid step changes the trend parameters by < 0.05", comparison_name,
                 if (is.null(cs$Check_finer_grid_max_par_change) || !is.finite(cs$Check_finer_grid_max_par_change)) NA else cs$Check_finer_grid_max_par_change < 0.05,
                 sprintf("largest parameter change %.3g (log a0, log a1, b); tau2 relative change %.3g",
                         cs$Check_finer_grid_max_par_change %||% NA_real_, cs$Check_finer_grid_max_tau2_rel_change %||% NA_real_))
  add_validation("Prior variances not at their lower bound", comparison_name,
                 !isTRUE(cs$LE_prior_var_at_bound) && !isTRUE(cs$REM_prior_var_at_bound),
                 sprintf("LE tau2 %.4f; REM tau2 %.4f; bound %.3g", cs$LE_prior_var, cs$REM_prior_var, DISP_PRIOR_VAR_MIN))
  validate_pca_preprocessing(comparison_name)

  # Stage 2: one significance analysis at the selected k (Original; LE + REM).
  message("Running ", comparison_name, " NormEVS at k* = ", OPTIMAL_K[[comparison_name]])
  ok <- run_one_evs_track(comparison_name = comparison_name, track_key = "normalized_evs",
                          count_matrix = count_matrix, coldata = coldata, annot_df = annot_df)
  if (!isTRUE(ok)) stop("Analysis failed: ", comparison_name)
  validate_selected_k_used(comparison_name)
  add_validation("Exactly one k* selected for this comparison", comparison_name,
                 sum(CUTOFF$profiles[[comparison_name]]$Chosen) == 1L && length(OPTIMAL_K[[comparison_name]]) == 1L,
                 sprintf("%s: k* = %d", EVS_PCA_METHOD, OPTIMAL_K[[comparison_name]]))

  tryCatch(collect_dispersion_after_evs(comparison_name),
           error = function(e) warning("Dispersion summaries skipped for ", comparison_name, ": ",
                                       conditionMessage(e), call. = FALSE))
  export_comparison_manuscript_tables(comparison_name)
  TRUE
}



# -----------------------------------------------------------------------------
# Paper-ready multi-comparison volcano panels
# -----------------------------------------------------------------------------

# Purpose: Read result table for panel.
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


# Purpose: Compute global volcano y limit.
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

# Purpose: Save paper volcano panels.
save_paper_volcano_panels <- function() {
  dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)
  views <- comparison_analysis_views()
  short_view <- c(
    "Original (No EVS)" = "Orig",
    "NormEVS Lead" = "N-Lead",
    "NormEVS Rem" = "N-Rem"
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

# Purpose: Pca variance stats.
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

# Purpose: Pc1 loading energy pct.
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

# Purpose: Pca track matrix.
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

# Purpose: Compute pca support plot.
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

# Purpose: Build evs split audit table.
build_evs_split_audit_table <- function() {
  rows <- list()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    for (track_key in "normalized_evs") {
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


# Purpose: Build evs pca evidence table.
build_evs_pca_evidence_table <- function(comparison_name) {
  rows <- list()

  for (track_key in "normalized_evs") {
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

# Purpose: Plot evs pca evidence.
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

# Purpose: Plot empirical hbfss support.
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

# Purpose: Save paper pca panels.
save_paper_pca_panels <- function() {
  if (!isTRUE(EXPORT_SUPPORT_FIGURES)) return(invisible(FALSE))

  all_evidence <- list()

  for (comparison_name in as.character(comparison_table$comparison_name)) {
    plots <- list()

    for (track_key in "normalized_evs") {
      obj <- paper_registry[[registry_key(comparison_name, track_key)]]
      if (is.null(obj) || is.null(obj$evs)) next

      original_mat <- pca_track_matrix(obj, "raw_dataset")
      original_stats <- pca_variance_stats(original_mat, obj$coldata)
      if (is.null(original_stats)) next


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
          short_title = ds_abbr,
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
        ncol = 3,
        legend = grid::textGrob(
          "Blue triangles: RT samples. Grey circles: ZT samples. X: condition centroid; lines join each sample to its centroid.",
          gp = grid::gpar(fontsize = 11, col = "grey20"))
      )
      base <- file.path(panel_dir, paste0("Figure_", comparison_name, "_PCA_EVS"))
      save_grob(panel, paste0(base, ".png"), width = 16.5, height = 6.0)
      save_grob(panel, paste0(base, ".pdf"), width = 16.5, height = 6.0)
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

# Purpose: Save paper empirical hbfss panels.
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
      track_order_use <- "normalized_evs"
      panel_title <- paste0("HBFSS calibration | ", unname(dataset_short[dataset_key]), " | ", CUTOFF_METHOD_LABEL)
      output_suffix <- paste0(unname(dataset_short[dataset_key]), "_AllComparisons")
      panel_height <- 5.8
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

# Purpose: Build discovery long table.
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


# Purpose: Save discovery count panel.
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

# Purpose: Save paper support figures.
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
# NORMALIZATION CHECK (compositional bias)
# =============================================================================
# DESeq2 median-of-ratios size factors are compared with edgeR TMM and
# upper-quartile effective library sizes (each rescaled to geometric mean 1).
# Large disagreement indicates composition effects (a few highly abundant or
# strongly changing PASs) to which the normalization is sensitive. Zero-count
# fractions and library sizes are reported per sample.
NORMALIZATION_CHECK_TOL <- log2(1.2)
normalization_rows <- list()

normalization_check <- function(comparison_name, raw, sf) {
  gm <- function(x) x / exp(mean(log(x)))
  y <- edgeR::DGEList(counts = raw)
  eff <- function(method) gm(y$samples$lib.size * edgeR::calcNormFactors(y, method = method)$samples$norm.factors)
  d <- data.frame(Comparison = comparison_name, Sample = colnames(raw), Library_size = colSums(raw),
                  Zero_fraction = colMeans(raw == 0), SF_DESeq2 = gm(sf[colnames(raw)]),
                  SF_TMM = eff("TMM"), SF_upper_quartile = eff("upperquartile"))
  d$log2_DESeq2_vs_TMM <- log2(d$SF_DESeq2 / d$SF_TMM)
  d$log2_DESeq2_vs_UQ <- log2(d$SF_DESeq2 / d$SF_upper_quartile)
  normalization_rows[[comparison_name]] <<- d
  m <- max(abs(c(d$log2_DESeq2_vs_TMM, d$log2_DESeq2_vs_UQ)))
  add_validation("Size factors agree with TMM and upper-quartile within 20% (compositional bias)", comparison_name,
                 m <= NORMALIZATION_CHECK_TOL, sprintf("largest |log2 ratio| %.3f", m))
  invisible(d)
}

# =============================================================================
# TIME-COURSE MODEL: RECOVERY TIME x CONDITION (all 40 samples)
# =============================================================================
# The design is a 4 x 2 factorial: recovery after sleep deprivation (R arms,
# treatment) and undisturbed controls at the same clock time (ZT arms), at
# recovery 0, 2, 4 and 8 h (ZT6, ZT8, ZT10, ZT14). One DESeq2 model,
# ~ time + condition + time:condition, with the hierarchical dispersion trend
# and experiment-wide size factors, gives:
#   any-time effect   LRT against ~ time (treatment changes the PAS at any time)
#   time-varying      LRT against ~ time + condition (the treatment effect
#                     differs between recovery times)
#   per-time effects  Wald log2 fold change and test of trt vs untrt at each time
# BH within each test, no independent filtering. This complements the four
# pairwise comparisons by testing time explicitly; it uses one (pooled) trend.
TIME_LEVELS <- c(T0 = "0 h (ZT6)", T2 = "2 h (ZT8)", T4 = "4 h (ZT10)", T8 = "8 h (ZT14)")
TIME_OF_ARM <- c(R0 = "T0", ZT6 = "T0", R2 = "T2", ZT8 = "T2", R4 = "T4", ZT10 = "T4", R8 = "T8", ZT14 = "T8")
TIME_COURSE <- NULL

run_time_course <- function() {
  raw <- experiment_counts()
  arm <- as.character(experiment_arm(colnames(raw)))
  cd <- data.frame(time = factor(unname(TIME_OF_ARM[arm]), levels = names(TIME_LEVELS)),
                   condition = meta_all[colnames(raw), "condition"], row.names = colnames(raw))
  if (anyNA(cd$time) || anyNA(cd$condition)) stop("Time-course design: unmapped samples.")
  full <- ~ time + condition + time:condition
  dds <- DESeq2::DESeqDataSetFromMatrix(raw, colData = cd, design = full)
  sf <- DESeq2::sizeFactors(experiment_dds())[colnames(raw)]
  dds <- estimate_dispersions_with_trend(dds, size_factors = sf)
  lrt <- function(reduced) {
    d <- DESeq2::nbinomLRT(dds, reduced = reduced, quiet = TRUE)
    as.data.frame(DESeq2::results(d, independentFiltering = FALSE))[, c("stat", "pvalue", "padj")]
  }
  any_t <- lrt(~ time); int_t <- lrt(~ time + condition)
  w <- DESeq2::nbinomWaldTest(dds, quiet = TRUE)
  SIGNIFICANCE_CALLS$wald <- SIGNIFICANCE_CALLS$wald + 3L
  rn <- DESeq2::resultsNames(w)
  main <- grep("^condition_", rn, value = TRUE)
  per <- lapply(names(TIME_LEVELS), function(tl) {
    inter <- grep(paste0("^time", tl, "\\.condition"), rn, value = TRUE)
    r <- DESeq2::results(w, contrast = list(c(main, inter)), independentFiltering = FALSE)
    data.frame(lfc = r$log2FoldChange, padj = r$padj)
  })
  names(per) <- names(TIME_LEVELS)
  out <- data.frame(feature_id = rownames(raw),
                    gene_symbol = OrigID_Symbol$gene_symbol[match(rownames(raw), OrigID_Symbol$feature_id)],
                    baseMean = S4Vectors::mcols(w)$baseMean,
                    LRT_any_time_stat = any_t$stat, LRT_any_time_p = any_t$pvalue, LRT_any_time_padj = any_t$padj,
                    LRT_time_varying_stat = int_t$stat, LRT_time_varying_p = int_t$pvalue,
                    LRT_time_varying_padj = int_t$padj)
  for (tl in names(per)) {
    out[[paste0("LFC_", tl)]] <- per[[tl]]$lfc
    out[[paste0("padj_", tl)]] <- per[[tl]]$padj
  }
  sig <- function(q) sum(!is.na(q) & q < BH_FDR_STANDARD)
  summ <- data.frame(Test = c("Treatment effect at any time (LRT)", "Time-varying treatment effect (LRT)",
                              paste("Treatment effect at", unname(TIME_LEVELS), "(Wald)")),
                     Significant_BH_0.10 = c(sig(out$LRT_any_time_padj), sig(out$LRT_time_varying_padj),
                                             vapply(names(TIME_LEVELS), function(tl) sig(out[[paste0("padj_", tl)]]), integer(1))),
                     PAS_tested = nrow(out),
                     Prior_variance = attr(DESeq2::dispersionFunction(dds), "dispPriorVar"))
  TIME_COURSE <<- list(pas = out, summary = summ)
  save_csv(out, file.path(summary_table_dir, "Table_TimeCourse_PAS.csv"))
  save_csv(summ, file.path(summary_table_dir, "Table_TimeCourse_Summary.csv"))
  try_fig(save_figure(plot_time_course_figure(), paper_fig_dir, "Figure_4_Time_Course", 180, 120), "Figure 4 (time course)")
  invisible(TIME_COURSE)
}

# Figure 4: (a) significant PASs per recovery time in the time-course model,
# next to the pairwise Original and EVS Standard counts; (b) treatment log2
# fold change over recovery time for the 40 PASs with the strongest
# time-varying effect.
plot_time_course_figure <- function() {
  s <- TIME_COURSE$summary; pas <- TIME_COURSE$pas
  per <- data.frame(Time = factor(unname(TIME_LEVELS), levels = unname(TIME_LEVELS)),
                    Analysis = "Time-course model (Wald, per time)",
                    n = s$Significant_BH_0.10[3:6])
  cnt <- build_overall_manuscript_summary()
  pw <- do.call(rbind, lapply(seq_len(nrow(comparison_table)), function(i) {
    cmp <- comparison_table$comparison_name[i]; z <- cnt[cnt$Comparison == cmp, , drop = FALSE]
    data.frame(Time = unname(TIME_LEVELS)[i], Analysis = c("Pairwise Original (Standard)", "Pairwise EVS (Standard)"),
               n = c(sum(z$Std[z$Analysis == ORIGINAL_LABEL]), sum(z$Std[z$Analysis %in% c("NormEVS Lead", "NormEVS Rem")])))
  }))
  d <- rbind(per, pw); d$Time <- factor(d$Time, levels = unname(TIME_LEVELS))
  d$Analysis <- factor(d$Analysis, levels = c("Time-course model (Wald, per time)", "Pairwise Original (Standard)",
                                              "Pairwise EVS (Standard)"))
  pa <- ggplot(d, aes(Time, n, fill = Analysis)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.74) +
    geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.3, size = 1.9) +
    scale_fill_manual(values = c("grey35", SET_COLS[["Original (No EVS)"]], SET_COLS[["Leading Edge"]]), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Recovery time", y = "Significant PASs (BH < 0.10)", tag = "a",
         title = "Treatment effect at each recovery time",
         subtitle = sprintf("Time-course model: %s PASs with an effect at any time; %s with a time-varying effect.",
                            format(s$Significant_BH_0.10[1], big.mark = ","), format(s$Significant_BH_0.10[2], big.mark = ","))) +
    fig_theme() + theme(panel.grid.major.x = element_blank())
  top <- pas[is.finite(pas$LRT_time_varying_p), , drop = FALSE]
  top <- head(top[order(top$LRT_time_varying_p), , drop = FALSE], 40L)
  long <- do.call(rbind, lapply(names(TIME_LEVELS), function(tl)
    data.frame(feature_id = top$feature_id, Time = unname(TIME_LEVELS[tl]), lfc = top[[paste0("LFC_", tl)]])))
  long$Time <- factor(long$Time, levels = unname(TIME_LEVELS))
  med <- stats::aggregate(lfc ~ Time, data = long, FUN = stats::median)
  pb <- ggplot(long, aes(Time, lfc, group = feature_id)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(colour = SET_COLS[["Leading Edge"]], alpha = 0.35, linewidth = 0.35) +
    geom_point(colour = SET_COLS[["Leading Edge"]], alpha = 0.5, size = 0.6) +
    labs(x = "Recovery time", y = "Treatment log2 fold change", tag = "b",
         title = "Strongest time-varying PASs",
         subtitle = "The 40 PASs with the smallest time-varying LRT p-value; one line per PAS.") +
    fig_theme()
  assemble_rows(list(list(plots = list(pa, pb), widths = c(1.2, 1))), heights = 1,
                title = "Figure 4. Recovery time x condition model across all 40 samples")
}

# =============================================================================
# FIGURE 3: MODEL COMPARISON AND SENSITIVITY
# =============================================================================
# (a) BIC of each alternative minus BIC of EVS (negative = better than EVS);
# (b) observed coarse-grid LRT vs the permutation searches; (c) Standard
# discoveries across the primary and sensitivity analyses.
plot_model_sensitivity_figure <- function() {
  s <- CUTOFF$summary
  alt <- c(BIC_pooled = "One trend (no split)", BIC_continuous_rank = "Continuous rank effect",
           BIC_variance_split = "Variance-ranked split", BIC_pc12_split = "PC1+PC2-ranked split",
           BIC_reference_ranking = "Independent reference ranking")
  da <- do.call(rbind, lapply(names(alt), function(v)
    data.frame(Comparison = s$Comparison, Model = unname(alt[v]), dBIC = s[[v]] - s$BIC_EVS)))
  da <- da[is.finite(da$dBIC), , drop = FALSE]
  da$Comparison <- as_comp(da$Comparison); da$Model <- factor(da$Model, levels = unname(alt))
  pa <- ggplot(da, aes(Comparison, dBIC, fill = Model)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_col(position = position_dodge(width = 0.82), width = 0.76) +
    scale_fill_manual(values = c("#2B2B2B", "#7F7F7F", "#E69F00", "#009E73", "#56B4E9"), name = NULL) +
    scale_x_discrete(labels = function(x) sub("_", "\n", x)) + scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "BIC minus BIC of EVS", tag = "a", title = "Model comparison on held-out samples",
         subtitle = "Positive: EVS fits the held-out dispersions better.") +
    fig_theme() + theme(panel.grid.major.x = element_blank())
  rnd <- CUTOFF$random; rnd <- rnd[is.finite(rnd$Coarse_LRT_vs_one_trend), , drop = FALSE]
  ob <- data.frame(Comparison = as_comp(s$Comparison), LRT = s$Coarse_grid_LRT_observed)
  pb <- ggplot() +
    geom_point(data = ob, aes(Comparison, LRT), colour = "#D55E00", size = 2.2, shape = 18) +
    scale_x_discrete(labels = function(x) sub("_", "\n", x), drop = FALSE) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "LRT, split vs one trend (coarse grid)", tag = "b", title = "Observed split vs permuted rankings",
         subtitle = "Diamond: observed PC1 ranking. Grey: permutation searches at the same resolution.") +
    fig_theme() + theme(panel.grid.major.x = element_blank())
  if (nrow(rnd)) pb <- pb + geom_jitter(data = transform(rnd, Comparison = as_comp(Comparison)),
                                        aes(Comparison, Coarse_LRT_vs_one_trend), colour = "grey45", width = 0.12, height = 0, size = 1)
  sens <- if (length(sensitivity_rows)) dplyr::bind_rows(sensitivity_rows) else NULL
  if (is.null(sens) || !nrow(sens)) return(assemble_rows(list(list(plots = list(pa, pb), widths = c(1.3, 1))), 1,
                                                       "Figure 3. Model comparison and sensitivity"))
  lv <- c(ORIGINAL_LABEL, "EVS (primary)", paste0(ORIGINAL_LABEL, ", DESeq2 0.25 prior floor"), "EVS, DESeq2 0.25 prior floor",
          "EVS, independent reference ranking", "IHW, Original p-values",
          "Count split: Original (test half)", "Count split: EVS (test half)")
  sens <- sens[sens$Analysis %in% lv, , drop = FALSE]
  sens$Analysis <- factor(sens$Analysis, levels = lv); sens$Comparison <- as_comp(sens$Comparison)
  cols <- c("#2B2B2B", "#1F6FB2", "#8C8C8C", "#8FB8DE", "#56B4E9", "#E69F00", "#B07AA1", "#7B3294")
  names(cols) <- lv
  bar <- function(var, tag, title, keep) {
    d <- sens[sens$Analysis %in% keep & is.finite(sens[[var]]), , drop = FALSE]
    ggplot(d, aes(Comparison, .data[[var]], fill = Analysis)) +
      geom_col(position = position_dodge(width = 0.86), width = 0.8) +
      scale_fill_manual(values = cols, name = NULL, drop = TRUE) +
      scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
      scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.08))) +
      labs(x = NULL, y = paste(var, "discoveries"), tag = tag, title = title) +
      fig_theme() + theme(panel.grid.major.x = element_blank(), legend.position = "right")
  }
  pc <- bar("Standard", "c", "Standard discoveries: primary and sensitivity analyses", lv)
  assemble_rows(list(list(plots = list(pa, pb), widths = c(1.3, 1)), pc), heights = c(1, 1),
                title = "Figure 3. Model comparison and sensitivity analyses")
}

# =============================================================================
# FINAL 3'aTWAS ORTHOLOG COMPARISON
# =============================================================================

# Orthology is a database-supported annotation, never inferred from case alone.
ensure_babelgene <- function() {
  if (!requireNamespace("babelgene",quietly=TRUE)) stop("Install babelgene before running the ortholog analysis.")
  TRUE
}

# Purpose: Read twas study.
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

# Purpose: Build twas ortholog map.
build_twas_ortholog_map <- function(twas) {
  human_genes <- sort(unique(twas$Human_TWAS_Gene))

  # Only database-supported mappings enter the biological overlap.


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
        stop("babelgene ortholog mapping failed: ", conditionMessage(e))
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

  if (!nrow(babel)) stop("No database-supported human-to-rat TWAS mappings.")
  babel %>%
    dplyr::mutate(gene_key = gene_key(Rat_Ortholog)) %>%
    dplyr::filter(!is.na(gene_key)) %>%
    dplyr::arrange(Human_TWAS_Gene, dplyr::desc(Ortholog_support_n), Mapping_source) %>%
    dplyr::distinct(Human_TWAS_Gene, Rat_Ortholog, .keep_all = TRUE)
}

# Purpose: Build twas rat metadata.
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

# Purpose: Collect twas analysis results.
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

# Purpose: Build twas overlap pas table.
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
        "NormEVS Lead", "NormEVS Rem"
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
        "NormEVS Lead", "NormEVS Rem"
      )),
      Rat_Ortholog,
      WTTS_DE_PAS
    )
}

# Purpose: Build twas overlap gene table.
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
        "NormEVS Lead", "NormEVS Rem"
      )),
      Rat_Ortholog
    )
}

# Purpose: Build twas overlap summary.
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

# Purpose: Plot twas gene support.
plot_twas_gene_support <- function(gene_table, figure_dir) {
  if (!nrow(gene_table)) return(invisible(NULL))

  analysis_map <- c(
    "Original (No EVS)"="Orig",
    "NormEVS Lead"="N-Lead",
    "NormEVS Rem"="N-Rem"
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

# Purpose: Plot twas view counts.
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
                        "  TWAS multi-transcript=", summary_row$TWAS_APA[1],
                        "  WTTS-DE-APA=", summary_row$WTTS_DE_APA[1]),
      x = NULL,
      y = "TWAS genes"
    ) +
    manuscript_theme() +
    theme(legend.position = "none", plot.margin = margin(10, 14, 10, 14))
}

# Purpose: Export twas by view.
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

# Purpose: Save twas comparison panels.
save_twas_comparison_panels <- function(summary_table) {
  for (comparison_name in as.character(comparison_table$comparison_name)) {
    sub <- summary_table[summary_table$Comparison == comparison_name, , drop = FALSE]
    if (!nrow(sub)) next

    preferred_views <- c("Original (No EVS)", "NormEVS Lead", "NormEVS Rem")
    views <- preferred_views[preferred_views %in% unique(as.character(sub$Analysis))]
    plots <- lapply(views, function(v) {
      one <- sub[sub$Analysis == v, , drop = FALSE]
      if (!nrow(one)) return(NULL)
      short <- c(
        "Original (No EVS)"="Orig",
        "NormEVS Lead"="N-Lead",
        "NormEVS Rem"="N-Rem"
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

# Purpose: Run twas overlap analysis.
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

# =============================================================================
# COMPLEMENTARY ANALYSES AND PRIMARY MODEL OUTPUT
# =============================================================================
# Export exactly one selected model per contrast; when the pooled criterion
# wins, EVS views are still retained as exploratory diagnostics. BH is per
# contrast, not jointly across contrasts or across test families.
export_selected_primary <- function() {
  rows <- lapply(COMPARISON_NAMES,function(cmp) {
    cs <- CUTOFF$summary[CUTOFF$summary$Comparison==cmp,,drop=FALSE]
    reg <- paper_registry[[registry_key(cmp,"normalized_evs")]]$analysis_results
    keys <- if(cs$Selected_Model=="EVS") c("leading_edge_dataset","remainder_dataset") else "raw_dataset"
    dplyr::bind_rows(lapply(keys,function(k) {
      z <- reg[[k]]$results; z$Comparison <- cmp; z$Selected_Model <- cs$Selected_Model; z$View <- k; z
    }))
  })
  d <- dplyr::bind_rows(rows)
  save_csv(d,file.path(summary_table_dir,"Table_Selected_Model_All_PAS.csv"))
  save_csv(d[d$standard_flag %in% TRUE,,drop=FALSE],file.path(summary_table_dir,"Table_Selected_Model_Primary_Discoveries.csv"))
  invisible(d)
}

# =============================================================================
# ALTERNATIVE POLYADENYLATION: RELATIVE PAS USAGE, 3' UTR SHORTENING /
# LENGTHENING, AND 3'aTWAS
# =============================================================================
# Differential abundance (EVS) does not establish a change in which PAS a gene
# uses. APA is analysed per pairwise comparison with edgeR quasi-likelihood
# models on raw counts with offsets from the shared DESeq2 size factors.
# Gene units are gene symbol + chromosome + strand; PASs are ordered along the
# direction of transcription by their peak coordinate (proximal = most 5').
#   1. Any usage change: diffSpliceDGE tests each PAS against the other PASs of
#      its gene (genes with >= 2 eligible PASs); gene level by Simes, BH across
#      genes.
#   2. Direction: for each gene, the two most used PASs (mean normalized count
#      >= APA_MIN_MEAN each) form a proximal/distal pair. RED = treatment log2
#      fold change of the distal PAS minus that of the proximal PAS (the change
#      in log2 distal/proximal usage; Hoque et al. 2013); positive =
#      lengthening (more distal use in recovery), negative = shortening. The
#      gene-level diffSpliceDGE test of the pair gives the p-value; BH across
#      genes, FDR < 0.10.
#   3. 3'aTWAS: enrichment of APA genes among 3'aTWAS genes (logistic, the
#      ortholog-mapped tested multi-PAS genes as background, adjusted for
#      expression and ortholog support).
APA_MIN_MEAN <- 5
APA_FDR <- 0.10
APA_RESULTS <- NULL

pas_annotation <- function() {
  need <- c("Chromosome", "Strand", "Peak")
  if (!all(need %in% names(WTTS_Seq))) stop("The count file lacks PAS coordinates (Chromosome, Strand, Peak).")
  d <- data.frame(feature_id = WTTS_Seq$feature_id, gene = trimws(as.character(WTTS_Seq$Symbol)),
                  chr = as.character(WTTS_Seq$Chromosome), strand = as.character(WTTS_Seq$Strand),
                  pos = suppressWarnings(as.numeric(gsub(",", "", as.character(WTTS_Seq$Peak)))))
  d <- d[valid_gene_symbol(d$gene) & d$strand %in% c("+", "-") & is.finite(d$pos), , drop = FALSE]
  d$gene_unit <- paste(d$gene, d$chr, d$strand, sep = "|")
  d$tx_pos <- ifelse(d$strand == "+", d$pos, -d$pos)       # increases 5' -> 3'
  d
}

apa_fit <- function(counts, genes, cache) {
  X <- stats::model.matrix(make_design_formula(cache$coldata), cache$coldata)
  yy <- edgeR::DGEList(counts = counts, genes = genes)
  yy$offset <- matrix(log(cache$sf[colnames(counts)]), nrow = nrow(yy), ncol = ncol(yy), byrow = TRUE)
  yy <- edgeR::estimateDisp(yy, X, robust = TRUE)
  fit <- edgeR::glmQLFit(yy, X, robust = TRUE)
  list(fit = fit, coef = ncol(X),
       ds = edgeR::diffSpliceDGE(fit, coef = ncol(X), geneid = "gene_unit", exonid = "PAS", verbose = FALSE))
}

run_apa_analysis <- function() {
  ann <- pas_annotation()
  pas_out <- list(); gene_out <- list(); red_out <- list(); stage_out <- list()
  for (cmp in COMPARISON_NAMES) {
    cache <- comparison_cache[[cmp]]; y <- cache$raw
    a <- ann[match(rownames(y), ann$feature_id), , drop = FALSE]
    ok <- !is.na(a$gene_unit)
    n_per <- table(a$gene_unit[ok]); multi <- ok & a$gene_unit %in% names(n_per[n_per >= 2])
    le <- CUTOFF$pas[[cmp]]; le_ids <- if (is.null(le)) character() else le$feature_id[le$In_Leading_Edge]
    # 1. Any usage change (all PASs of multi-PAS genes).
    if (sum(multi) >= 2L) {
      f <- apa_fit(y[multi, , drop = FALSE],
                   data.frame(gene_unit = a$gene_unit[multi], Gene = a$gene[multi], PAS = rownames(y)[multi]), cache)
      pt <- edgeR::topSpliceDGE(f$ds, test = "exon", number = Inf, FDR = 1); pt$Comparison <- cmp
      gt <- edgeR::topSpliceDGE(f$ds, test = "Simes", number = Inf, FDR = 1); gt$Comparison <- cmp
      # Stage-wise testing (stageR, Van den Berge et al. 2017): genes screened
      # by the Simes p-value, then PASs confirmed within significant genes, with
      # overall FDR control at APA_FDR. Run when the stageR package is installed.
      if (requireNamespace("stageR", quietly = TRUE)) {
        st <- tryCatch({
          s <- stageR::stageRTx(pScreen = stats::setNames(gt$P.Value, gt$gene_unit),
                                pConfirmation = matrix(pt$P.Value, ncol = 1, dimnames = list(pt$PAS, "usage")),
                                pScreenAdjusted = FALSE, tx2gene = data.frame(pt$PAS, pt$gene_unit))
          s <- stageR::stageWiseAdjustment(s, method = "dtu", alpha = APA_FDR)
          out_st <- stageR::getAdjustedPValues(s, onlySignificantGenes = FALSE, order = FALSE)
          out_st$Comparison <- cmp; out_st
        }, error = function(e) { warning(cmp, ": stageR failed: ", conditionMessage(e), call. = FALSE); NULL })
        if (!is.null(st)) stage_out[[cmp]] <- st
      }
      pas_out[[cmp]] <- pt; gene_out[[cmp]] <- gt
    }
    # 2. Proximal / distal pair of the two most used PASs.
    nm <- rowMeans(sweep(y, 2, cache$sf[colnames(y)], "/"))
    cand <- data.frame(feature_id = rownames(y)[ok], unit = a$gene_unit[ok], gene = a$gene[ok],
                       tx = a$tx_pos[ok], pos = a$pos[ok], mean = nm[ok])
    cand <- cand[cand$mean >= APA_MIN_MEAN, , drop = FALSE]
    cand <- cand[order(cand$unit, -cand$mean, cand$tx), , drop = FALSE]
    top2 <- do.call(rbind, lapply(split(cand, cand$unit), function(z) if (nrow(z) >= 2L) z[1:2, , drop = FALSE] else NULL))
    if (is.null(top2) || nrow(top2) < 4L) next
    top2 <- top2[order(top2$unit, top2$tx), , drop = FALSE]
    top2$role <- rep(c("proximal", "distal"), nrow(top2) / 2)
    f2 <- apa_fit(y[top2$feature_id, , drop = FALSE],
                  data.frame(gene_unit = top2$unit, Gene = top2$gene, PAS = top2$feature_id), cache)
    lfc <- f2$fit$coefficients[, f2$coef] / log(2)
    g2 <- edgeR::topSpliceDGE(f2$ds, test = "gene", number = Inf, FDR = 1)
    prox <- top2[top2$role == "proximal", , drop = FALSE]; dist <- top2[top2$role == "distal", , drop = FALSE]
    red <- data.frame(Comparison = cmp, Gene = prox$gene, Gene_unit = prox$unit,
                      Proximal_PAS = prox$feature_id, Distal_PAS = dist$feature_id,
                      Distance_bp = abs(dist$pos - prox$pos), Mean_proximal = prox$mean, Mean_distal = dist$mean,
                      RED_log2 = unname(lfc[match(dist$feature_id, top2$feature_id)] - lfc[match(prox$feature_id, top2$feature_id)]))
    red$P <- g2$P.Value[match(red$Gene_unit, g2$gene_unit)]
    red$FDR <- stats::p.adjust(red$P, "BH")
    red$Significant <- is.finite(red$FDR) & red$FDR < APA_FDR
    red$Direction <- ifelse(red$RED_log2 > 0, "Lengthening", "Shortening")
    red$Proximal_in_LE <- red$Proximal_PAS %in% le_ids; red$Distal_in_LE <- red$Distal_PAS %in% le_ids
    red_out[[cmp]] <- red
  }
  pas_tab <- dplyr::bind_rows(pas_out); gene_tab <- dplyr::bind_rows(gene_out); red_tab <- dplyr::bind_rows(red_out)
  save_csv(pas_tab, file.path(summary_table_dir, "Table_APA_PAS_Usage.csv"))
  save_csv(gene_tab, file.path(summary_table_dir, "Table_APA_Gene_Usage_Simes.csv"))
  save_csv(red_tab, file.path(summary_table_dir, "Table_APA_Shortening_Lengthening.csv"))
  save_csv(if (length(stage_out)) dplyr::bind_rows(stage_out) else data.frame(Note = "stageR not installed or not run"),
           file.path(summary_table_dir, "Table_APA_StageR.csv"))
  summ <- do.call(rbind, lapply(COMPARISON_NAMES, function(cmp) {
    r <- red_tab[red_tab$Comparison == cmp, , drop = FALSE]; g <- gene_tab[gene_tab$Comparison == cmp, , drop = FALSE]
    data.frame(Comparison = cmp, Genes_usage_tested = nrow(g),
               Genes_usage_changed = if (nrow(g)) sum(g$FDR < APA_FDR, na.rm = TRUE) else 0L,
               Genes_pair_tested = nrow(r), Shortening = sum(r$Significant & r$Direction == "Shortening"),
               Lengthening = sum(r$Significant & r$Direction == "Lengthening"),
               Significant_with_PAS_in_LE = sum(r$Significant & (r$Proximal_in_LE | r$Distal_in_LE)))
  }))
  save_csv(summ, file.path(summary_table_dir, "Table_APA_Summary.csv"))
  APA_RESULTS <<- list(pas = pas_tab, gene = gene_tab, red = red_tab, summary = summ, enrichment = NULL)
  invisible(APA_RESULTS)
}

# 3'aTWAS enrichment among APA genes (significant shortening or lengthening).
run_apa_twas_enrichment <- function(mapping) {
  red <- APA_RESULTS$red
  if (is.null(red) || !nrow(red)) stop("No APA results.")
  rat <- unique(red$Gene[valid_gene_symbol(red$Gene)])
  bg <- as.data.frame(babelgene::orthologs(genes = rat, species = TWAS_TARGET_SPECIES, human = FALSE,
                                           min_support = TWAS_ORTHOLOG_MIN_SUPPORT, top = FALSE))
  assert_required_columns(bg, c("symbol", "human_symbol"), "rat background orthology (APA)")
  bgkey <- unique(gene_key(bg$symbol)); hit <- unique(mapping$gene_key)
  sup <- if ("support_n" %in% names(bg)) tapply(suppressWarnings(as.numeric(bg$support_n)), gene_key(bg$symbol), max) else NULL
  rows <- lapply(split(red, red$Comparison), function(z) {
    z$gene_key <- gene_key(z$Gene)
    g <- z |> dplyr::filter(gene_key %in% bgkey, is.finite(P)) |> dplyr::group_by(gene_key) |>
      dplyr::summarise(APA = any(Significant), mean_count = mean(Mean_proximal + Mean_distal), .groups = "drop")
    g$TWAS <- g$gene_key %in% hit
    g$support <- if (is.null(sup)) 1 else as.numeric(sup[g$gene_key]); g$support[!is.finite(g$support)] <- 1
    if (nrow(g) < 10L || length(unique(g$TWAS)) < 2L || length(unique(g$APA)) < 2L) return(NULL)
    form <- if (length(unique(g$support)) > 1L) APA ~ TWAS + log1p(mean_count) + log(support) else APA ~ TWAS + log1p(mean_count)
    f <- stats::glm(form, data = g, family = stats::binomial())
    cc <- summary(f)$coefficients
    if (!"TWASTRUE" %in% rownames(cc) || !f$converged) return(NULL)
    b <- cc["TWASTRUE", 1]; se <- cc["TWASTRUE", 2]
    data.frame(Comparison = z$Comparison[1], Genes = nrow(g), TWAS_genes = sum(g$TWAS), APA_genes = sum(g$APA),
               Overlap = sum(g$APA & g$TWAS), Overlap_genes = paste(sort(unique(z$Gene[gene_key(z$Gene) %in% g$gene_key[g$APA & g$TWAS]])), collapse = ";"),
               OR = exp(b), Lower95 = exp(b - 1.96 * se), Upper95 = exp(b + 1.96 * se), P = cc["TWASTRUE", 4])
  })
  d <- dplyr::bind_rows(rows); if (nrow(d)) d$BH <- stats::p.adjust(d$P, "BH")
  save_csv(d, file.path(summary_table_dir, "Table_APA_TWAS_Enrichment.csv"))
  APA_RESULTS$enrichment <<- d
  invisible(d)
}

# Figure 5: (a) genes with a usage change and with significant shortening or
# lengthening per recovery time; (b) RED volcano per comparison; (c) 3'aTWAS
# enrichment of APA genes (odds ratio and 95% interval), when available.
plot_apa_figure <- function() {
  s <- APA_RESULTS$summary; red <- APA_RESULTS$red
  tl <- stats::setNames(unname(TIME_LEVELS), COMPARISON_NAMES)
  d <- rbind(data.frame(Comparison = s$Comparison, Class = "Any usage change (gene, Simes)", n = s$Genes_usage_changed),
             data.frame(Comparison = s$Comparison, Class = "Shortening", n = s$Shortening),
             data.frame(Comparison = s$Comparison, Class = "Lengthening", n = s$Lengthening))
  d$Time <- factor(tl[d$Comparison], levels = unname(TIME_LEVELS))
  d$Class <- factor(d$Class, levels = c("Any usage change (gene, Simes)", "Shortening", "Lengthening"))
  cols <- c("Any usage change (gene, Simes)" = "grey40", "Shortening" = "#C0392B", "Lengthening" = "#1F6FB2")
  pa <- ggplot(d, aes(Time, n, fill = Class)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.74) +
    geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.3, size = 1.9) +
    scale_fill_manual(values = cols, name = NULL) + scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Recovery time", y = sprintf("Genes (FDR < %.2f)", APA_FDR), tag = "a",
         title = "Alternative polyadenylation in recovery",
         subtitle = "Shortening / lengthening: two most used PASs per gene (RED). Any change: all PASs of the gene.") +
    fig_theme() + theme(panel.grid.major.x = element_blank())
  red$Time <- factor(tl[red$Comparison], levels = unname(TIME_LEVELS))
  red$Class <- ifelse(red$Significant, red$Direction, "Not significant")
  red$Class <- factor(red$Class, levels = c("Shortening", "Lengthening", "Not significant"))
  red <- red[is.finite(red$RED_log2) & is.finite(red$FDR), , drop = FALSE]
  red <- red[order(red$Class, decreasing = TRUE), , drop = FALSE]
  pb <- ggplot(red, aes(RED_log2, -log10(pmax(FDR, 1e-300)), colour = Class)) +
    geom_point(size = 0.5, alpha = 0.6, stroke = 0) +
    geom_hline(yintercept = -log10(APA_FDR), linetype = "22", colour = "grey45", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    facet_wrap(~ Time, nrow = 1) +
    scale_colour_manual(values = c("Shortening" = "#C0392B", "Lengthening" = "#1F6FB2", "Not significant" = "grey75"), name = NULL) +
    labs(x = "RED: change in log2 distal / proximal usage (recovery vs control)", y = "-log10 FDR", tag = "b",
         title = "Direction of APA changes",
         subtitle = "Positive: more distal PAS use (3' UTR lengthening); negative: more proximal (shortening).") +
    fig_theme()
  rows <- list(pa, pb)
  en <- APA_RESULTS$enrichment
  if (!is.null(en) && nrow(en)) {
    en$Time <- factor(tl[en$Comparison], levels = unname(TIME_LEVELS))
    rows[[3]] <- ggplot(en, aes(Time, OR)) +
      geom_hline(yintercept = 1, linetype = "22", colour = "grey45", linewidth = 0.3) +
      geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0, linewidth = 0.4) +
      geom_point(size = 2, colour = "#7B3294") +
      geom_text(aes(label = sprintf("%d genes", Overlap)), vjust = -1, size = 1.9) +
      scale_y_log10() +
      labs(x = "Recovery time", y = "Odds ratio (95% CI)", tag = "c",
           title = "APA genes among 3'aTWAS brain-disorder genes",
           subtitle = "Logistic regression with ortholog-mapped tested genes as background, adjusted for expression and ortholog support.") +
      fig_theme()
  }
  assemble_rows(rows, heights = rep(1, length(rows)), title = "Figure 5. Alternative polyadenylation after sleep deprivation")
}

# Gene-level 3'aTWAS overlap enrichment. Map the entire tested rat background
# to human genes, not just TWAS hits. Collapse mapping to one rat gene per row.
# Logistic adjustment for mean abundance and number of tested PASs reduces
# differential detectability bias; one-to-many mappings are reported separately.
# This is exploratory enrichment, not causal or cross-species validation.
run_twas_enrichment <- function(selected, mapping) {
  rat <- unique(selected$gene_symbol[valid_gene_symbol(selected$gene_symbol)])
  bg <- babelgene::orthologs(genes=rat,species=TWAS_TARGET_SPECIES,human=FALSE,
                              min_support=TWAS_ORTHOLOG_MIN_SUPPORT,top=FALSE)
  bg <- as.data.frame(bg); assert_required_columns(bg,c("symbol","human_symbol"),"rat background orthology")
  bgkey <- unique(gene_key(bg$symbol)); hit <- unique(mapping$gene_key)
  rows <- lapply(split(selected,selected$Comparison),function(z) {
    z$gene_key <- gene_key(z$gene_symbol)
    g <- z |>
      dplyr::filter(gene_key %in% bgkey,is.finite(pvalue)) |>
      dplyr::group_by(gene_key) |>
      dplyr::summarise(DE=any(standard_flag),nPAS=dplyr::n(),mean_count=mean(baseMean),.groups="drop")
    g$TWAS <- g$gene_key %in% hit
    # Ortholog mapping confidence: the largest number of supporting databases
    # for the gene's human ortholog(s).
    sup <- if ("support_n" %in% names(bg)) tapply(suppressWarnings(as.numeric(bg$support_n)), gene_key(bg$symbol), max) else NULL
    g$support <- if (is.null(sup)) 1 else as.numeric(sup[g$gene_key])
    g$support[!is.finite(g$support)] <- 1
    g$ortholog_count <- vapply(g$gene_key,function(k) length(unique(bg$human_symbol[gene_key(bg$symbol)==k])),integer(1))
    # Restrict to one-to-one mapping sensitivity as well as all supported mappings.
    estimate <- function(g,scope) {
      if(nrow(g)<10L || length(unique(g$TWAS))<2L || length(unique(g$DE))<2L) return(NULL)
      form <- if (length(unique(g$support)) > 1L) DE~TWAS+log1p(nPAS)+log1p(mean_count)+log(support) else
        DE~TWAS+log1p(nPAS)+log1p(mean_count)
      f <- stats::glm(form,data=g,family=stats::binomial())
      cc <- summary(f)$coefficients; if(!"TWASTRUE" %in% rownames(cc) || !f$converged) return(NULL)
      b<-cc["TWASTRUE",1]; se<-cc["TWASTRUE",2]
      data.frame(Comparison=z$Comparison[1],Scope=scope,Genes=nrow(g),TWAS_genes=sum(g$TWAS),
                 DE_genes=sum(g$DE),Overlap=sum(g$DE & g$TWAS),OR=exp(b),Lower95=exp(b-1.96*se),Upper95=exp(b+1.96*se),P=cc["TWASTRUE",4])
    }
    # A strict reciprocal 1:1 map within the eligible ortholog background.
    unique_pairs <- unique(bg[,c("symbol","human_symbol")]); hcount <- table(unique_pairs$human_symbol)
    single_h <- unique_pairs$human_symbol %in% names(hcount[hcount==1])
    one_keys <- gene_key(unique_pairs$symbol[single_h])
    dplyr::bind_rows(estimate(g,"all_supported"),estimate(g[g$ortholog_count==1 & g$gene_key %in% one_keys,,drop=FALSE],"one_to_one_within_background"))
  })
  d <- dplyr::bind_rows(rows); if(nrow(d)) d$BH <- stats::p.adjust(d$P,"BH")
  save_csv(d,file.path(summary_table_dir,"Table_TWAS_Enrichment.csv")); invisible(d)
}

# Reproducibility manifest: exact input hashes, software versions, settings.
write_provenance <- function() {
  writeLines(capture.output(sessionInfo()),file.path(output_dir,"SessionInfo.txt"))
  files <- c(count_file,resolve_existing_file(twas_file_candidates,"3'aTWAS file"),get_script_path())
  save_csv(data.frame(Path=files,MD5=unname(tools::md5sum(files))),file.path(output_dir,"Input_Checksums.csv"))
  writeLines(c(paste("Build",PIPELINE_BUILD),paste("Grid step",DISP_GRID_STEP),paste("Minimum prior variance",DISP_PRIOR_VAR_MIN),
               paste("Cores",SPLIT_CORES),paste("Optimizer max iterations",TREND_MAXIT),
               paste("Random searches",SPLIT_RANDOM_SPLITS),paste("Random seed",SPLIT_SEED),paste("Minimum PAS count",MIN_PAS_COUNT),
               paste("Minimum PAS samples",MIN_PAS_SAMPLES)),file.path(output_dir,"Analysis_Settings.txt"))
}

# The cutoff k* of every comparison, before any DE test (EVS CUTOFF k*).
cutoff_calls_before <- c(SIGNIFICANCE_CALLS$wald, SIGNIFICANCE_CALLS$null)
compute_cutoffs()
cutoff_calls_after <- c(SIGNIFICANCE_CALLS$wald, SIGNIFICANCE_CALLS$null)
add_validation("Cutoff selection performed no significance testing", "All comparisons",
               identical(cutoff_calls_before, cutoff_calls_after),
               sprintf("Wald fits %d -> %d; significance calls %d -> %d", cutoff_calls_before[1],
                       cutoff_calls_after[1], cutoff_calls_before[2], cutoff_calls_after[2]))

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
  overall_summary$EVS_PCA_Method <- EVS_PCA_METHOD
  save_csv(overall_summary, file.path(summary_table_dir, "Table_DE_Method_Counts.csv"))
}

empirical_cutoff_export <- comparison_table[, c("comparison_name", "group1_prefix", "group2_prefix"), drop = FALSE]
names(empirical_cutoff_export) <- c("Comparison", "RT_prefix", "ZT_prefix")
empirical_cutoff_export$Active_Cutoff_Method <- CUTOFF_METHOD_LABEL
empirical_cutoff_export$Active_k_per_condition <- vapply(
  empirical_cutoff_export$Comparison,
  function(cn) tryCatch(get_active_evs_cutoff(cn), error = function(e) NA_integer_), integer(1)
)
empirical_cutoff_export$Basis <- EMPIRICAL_EVS_CUTOFF_BASIS
empirical_cutoff_export$Applied_to <- "NormEVS"
save_csv(
  empirical_cutoff_export,
  file.path(summary_table_dir, "Table_EVS_Cutoffs.csv")
)

evs_split_audit <- build_evs_split_audit_table()
if (nrow(evs_split_audit) > 0L) {
  save_csv(
    evs_split_audit,
    file.path(summary_table_dir, "Table_EVS_Split_Audit.csv")
  )
}

write_methods_note()
selected_primary <- export_selected_primary()
save_paper_volcano_panels()
save_paper_support_figures(overall_summary)
tryCatch(save_dispersion_outputs(),
         error = function(e) warning("Dispersion outputs skipped: ", conditionMessage(e), call. = FALSE))
tryCatch(save_run_checks(),
         error = function(e) warning("Validation outputs skipped: ", conditionMessage(e), call. = FALSE))

# Final analytical stage: 3'aTWAS ortholog overlap with every reported WTTS method.
twas_results <- tryCatch(run_twas_overlap_analysis(), error = function(e) {
  warning("3'aTWAS overlap failed: ", conditionMessage(e), call. = FALSE); NULL })
# Record complementary-stage errors explicitly; primary abundance results are
# still useful, but the manuscript must not claim a missing analysis ran.
complementary_errors <- character()
for (nm in c("Time_course","APA","TWAS_enrichment","APA_TWAS_enrichment","APA_figure")) tryCatch({
  if (nm=="Time_course") run_time_course() else
  if (nm=="APA") run_apa_analysis() else
  if (nm=="APA_figure") { if (is.null(APA_RESULTS)) stop("APA analysis did not run")
    save_figure(plot_apa_figure(), paper_fig_dir, "Figure_5_APA", 180, 170) } else
  if (is.null(twas_results)) stop("3'aTWAS overlap did not run") else
  if (nm=="APA_TWAS_enrichment") { if (is.null(APA_RESULTS)) stop("APA analysis did not run"); run_apa_twas_enrichment(twas_results$mapping) } else
    run_twas_enrichment(selected_primary,twas_results$mapping)
},error=function(e) {complementary_errors <<- c(complementary_errors,paste(nm,conditionMessage(e),sep=": ")); warning(tail(complementary_errors,1))})
writeLines(if(length(complementary_errors)) complementary_errors else "All complementary stages completed (time course, APA, TWAS enrichment, APA-TWAS enrichment, Figure 5).",file.path(output_dir,"Complementary_Status.txt"))
write_provenance()

figure_zip <- create_all_figures_zip()
table_zip <- create_all_tables_zip()

cat("\n=====================================================\n")
cat("Pipeline complete.\n")
cat("Build: ", PIPELINE_BUILD, "\n", sep = "")
cat("Cutoff method: ", CUTOFF_METHOD_LABEL, "\n", sep = "")
cat("Repository root:\n", repo_root, "\n", sep = "")
cat("Output directory:\n", output_dir, "\n", sep = "")
cat("3'aTWAS overlap: ", if (is.null(twas_results)) "FAILED (see warnings)" else "complete", "\n", sep = "")
cat("Figure ZIP:\n", figure_zip, "\n", sep = "")
cat("Table ZIP:\n", table_zip, "\n", sep = "")
cat("=====================================================\n\n")

if (nrow(overall_summary) > 0L) print(overall_summary)
