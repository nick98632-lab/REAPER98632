#!/usr/bin/env Rscript
# =============================================================================
# SEQUENCE EVS manuscript pipeline — hierarchical dispersion model (V9, final)
# =============================================================================
# Run: SEQUENCE_CORES=1 Rscript SEQUENCE_EVS_2MANUSCRIPT.R
# Set SEQUENCE_REPO_ROOT to the directory containing the original WTTS count
# matrix and 3'aTWAS input. Output: exports/sequence_evs_final (and .zip),
# written to a staging directory and swapped in only after a successful run.
# Order of the analysis:
#   1. design checks; per comparison (ten samples): abundance filter (>= 10
#      reads in >= 5 samples), shared size factors, profiled Cox-Reid
#      likelihoods, comparison VST, per-condition PC1 rankings, k* by the
#      in-sample hierarchical split likelihood, pooled-vs-split selection by
#      parametric bootstrap, numerical checks;
#   2. per comparison: pooled hierarchical model and EVS (Leading Edge +
#      Remainder, separate moderation; reported whether or not selected),
#      Standard / Strong / Weak tests with one BH per analysis, exploratory
#      HBFSS; sensitivity analyses; DESeq2 default, harmonized DESeq2, edgeR
#      robust QL, TMM and alternative partitions;
#   3. experiment level: multiplicity, the 40-sample time-course model, APA,
#      3'aTWAS enrichment;
#   4. validation: negative binomial simulations (7 scenarios), cutoff
#      stability, filter sensitivity;
#   5. external validation on GSE100349 (rat WTTS-seq) and GSE151173 (mouse
#      liver 3' mRNA-seq, WT): discoveries, real-data null splits, replication;
#   6. figures 1-8, supplementary figures S1-S5, tables, lockfile.
# Primary endpoint: ordinary two-sided Wald/BH per comparison and analysis.
# Secondary: Strong (|LFC| > 1) and Weak (composite null) tests; HBFSS is
# exploratory. EVS's FDR is not guaranteed by theory; it is measured by the
# validation stage. No setting is tuned on discoveries.
# Environment options (defaults: development / publication profile):
#   SEQUENCE_PROFILE (development | publication), SEQUENCE_CORES (1 on a 4 GB
#   server), SEQUENCE_SELECTION_BOOT (19 / 99), SEQUENCE_SIM_REPS (2 / 1000),
#   SEQUENCE_SIM_NPAS (3000 / 10000), SEQUENCE_SIM_SCENARIOS (all 7),
#   SEQUENCE_SIM_MODELS, SEQUENCE_STABILITY_SUBSETS (3 / 25),
#   SEQUENCE_STABILITY_METHOD (subsample | bootstrap),
#   SEQUENCE_FILTER_THRESHOLDS ("1x1,20x5" / "1x1,5x3,20x5"),
#   SEQUENCE_MIN_PAS_COUNT (10), SEQUENCE_MIN_PAS_SAMPLES (5),
#   SEQUENCE_PROFILED_LIKELIHOOD (1), SEQUENCE_VALIDATION (1),
#   SEQUENCE_RANDOM_PARTITIONS (3), SEQUENCE_GTF (optional rat GTF),
#   SEQUENCE_WEAK_DELTA (0), SEQUENCE_DISP_GRID_STEP (0.25),
#   SEQUENCE_PRIOR_VAR_MIN (0.01), SEQUENCE_NULL_METHOD (fndr),
#   SEQUENCE_EXTERNAL (1), SEQUENCE_EXTERNAL_WTTS / SEQUENCE_EXTERNAL_LIVER (paths
#   to GSE100349 / GSE151173 count tables; default data/<GEO file name>),
#   SEQUENCE_EXTERNAL_NULL_SPLITS
#   (3 / 60), SEQUENCE_EXTERNAL_REPLICATION_SPLITS (2 / 10),
#   SEQUENCE_FDR_TOLERANCE (0.025),
#   SEQUENCE_SELFTEST (1 = unit tests on synthetic data only).
# Methods text: MANUSCRIPT_METHODS.
# =============================================================================

MANUSCRIPT_METHODS <- c(
  "# SEQUENCE EVS | Statistical methods",
  "",
  "**Hierarchical dispersion model with eigenvector splitting**  ",
  "Methods accompanying `SEQUENCE_EVS_2MANUSCRIPT.R` (build SEQUENCE_EVS_MANUSCRIPT_HIERARCHICAL_V9_FINAL) · 27 November 2026",
  "",
  "This document describes the algorithm as implemented. Numerical results are written by the pipeline to the Summary_Tables folder and are not part of this text. The novel components are the eigenvector-guided partition and its marginal adjusted-likelihood cutoff; the cited literature supports the constituent statistical methods, not EVS itself. Because EVS chooses its partition from the counts it tests, its false discovery rate (FDR) is not guaranteed by theory; it is assessed empirically by full-pipeline simulation and on real null data from an independent dataset, and nominal and empirical FDR are reported separately.",
  "",
  "## Study design and input",
  "",
  "The analytical unit is a WTTS-Seq polyadenylation site (PAS), identified by OrigID; gene symbols are annotations, never units of differential testing. The study is a 4 × 2 factorial time course: rats in recovery after sleep deprivation (R arms, treatment) sampled at 0, 2, 4 and 8 h of recovery, and undisturbed controls sampled at the same clock times (ZT6, ZT8, ZT10 and ZT14), five animals per arm (40 samples, eight arms). Four prespecified pairwise comparisons are analysed: RT0 vs ZT6, RT2 vs ZT8, RT4 vs ZT10 and RT8 vs ZT14; the full design is analysed jointly in a separate time-course model. For each comparison a PAS is eligible if it has at least 10 reads in at least 5 of the comparison's ten samples (prespecified; the size of one group; designs with other group sizes, in the external datasets, require 10 reads in as many samples as one group holds). The rule ignores condition labels, so it cannot favour PASs that differ between groups, and it removes PASs too sparse to carry information on dispersion. The same universe is used for every model and for the 40-sample time-course model (applied across all 40 samples). Sensitivity to the filter is assessed by repeating the pooled and forced EVS analyses (filter, ranking, cutoff search, moderation and tests) at ≥1 read in ≥1 sample and ≥20 reads in ≥5 samples (and ≥5 reads in ≥3 samples in the publication profile). Before any analysis the pipeline checks that every arm has five distinct animals assigned to the correct condition, that counts are non-negative integers, and that the design matrix has full rank. Each animal is sampled once and the metadata record no pairing or batch, so none is modelled; the pipeline screens for duplicate and near-identical sequencing libraries, which detects duplicated sequencing but does not by itself establish independence; independence rests on the design (one animal per sample). [1, 5]",
  "",
  "## Normalization and eigenvector ranking",
  "",
  "All steps of a comparison use only its ten samples. DESeq2 median-of-ratios size factors are estimated once per comparison on the eligible comparison matrix (poscounts if the standard estimator fails), and the same factors and gene-wise dispersion estimates are supplied to every model, so the models differ only in their dispersion moderation. Raw counts, never transformed values, enter the negative binomial models. As a check for compositional bias, the size factors are compared with edgeR TMM and upper-quartile effective library sizes (each rescaled to geometric mean 1), with agreement within 20% for every sample recorded as a validation check. [1, 11, 21]",
  "",
  "For ranking only, a blind variance-stabilizing transformation (VST) is computed on the ten samples of the comparison with the shared size factors. Within each condition, PCA (centred, unscaled) is applied to the PAS-by-sample VST matrix and PASs are positioned by |PC1 loading| (1 = largest). The union depth dᵢ of PAS i is the smaller of its two positions. At cutoff k, the Leading Edge (LE) is every eligible PAS with dᵢ ≤ k (the union of the two conditions' top-k sets) and the Remainder (REM) is the complement; the partition is exhaustive and disjoint, and ties follow a fixed feature order. Samples from other time points are not used, because their dispersion structure need not match the tested time point. The ranking uses no condition contrast, test statistic or p-value, but it is computed from the counts that are tested, which is why error rates are assessed empirically (Validation). [1, 5]",
  "",
  "## Hierarchical dispersion likelihood",
  "",
  "For PAS i in sample j the negative binomial model has mean μᵢⱼ = sⱼ exp(xⱼᵀβᵢ) and variance μᵢⱼ + αᵢμᵢⱼ², with design ~ condition. The log dispersion zᵢ = log αᵢ has a population-specific normal prior with centre mᵢ = log(a₀ + a₁ μ̄ᵢ⁻ᵇ) and variance τ²; a₀ > 0, a₁ > 0 and b are estimated. The trend extends DESeq2's parametric trend (b = 1). [1–3]",
  "",
  "Each PAS's Cox–Reid adjusted profile log likelihood is ℓᵢ(z) = Σⱼ log NB(yᵢⱼ; μ̂ᵢⱼ(z), eᶻ) − ½ log det(XᵀWᵢX), with Wᵢ = diag(μ̂ᵢⱼ/(1 + eᶻμ̂ᵢⱼ)); for the saturated two-group design the determinant is the product of the within-group weight sums. The group means μ̂ᵢⱼ(z) are re-estimated by maximum likelihood at every grid dispersion (damped Newton steps on the log mean, warm-started from the neighbouring grid point), so the likelihood is properly profiled over the means. The plug-in version, which holds the means at DESeq2's fitted values, is computed at k* only as a check, with the median change in log dispersion, the resulting k* and the Jaccard similarity of the Leading Edges reported. [1, 2, 17]",
  "",
  "ℓᵢ is evaluated on a fixed grid of log dispersions from log(10⁻⁶) to log(10³) with spacing 0.25. Normal CDF differences assign the prior mass to grid cells, and the marginal objective is Σᵢ log Σₕ Lᵢ(zₕ) P(zᵢ ∈ cell h | mᵢ, τ²), with Lᵢ the relative likelihood. PASs with weakly informative likelihoods therefore contribute less, boundary gene-wise estimates are retained, and no uniform bias correction is needed. [18]",
  "",
  "The parameters are estimated by bound-constrained L-BFGS-B with analytic gradients (at most 1,000 iterations; the gradient is checked against central differences at every fitted optimum, with a maximum relative error below 10⁻³ required), with a₀ ∈ [10⁻⁸, 10³], a₁ ∈ [10⁻⁸, 10⁸], b ∈ [0.5, 3.5] and τ² ∈ [0.01, 16]. The lower bound 0.01 regularizes against collapse of τ² to zero and is identical in the cutoff search and in the final moderation; an estimate at the bound is flagged. The estimated τ² is the prior variance used in DESeq2's MAP shrinkage and is compared with DESeq2's conventional 0.25 floor in a reference analysis. That DESeq2 uses the intended prior variance is verified for every fitted view: re-running the MAP step with the stored prior variance reproduces the final dispersions, and re-running it with four times that variance changes them. Numerical adequacy is checked for every comparison at k*: the LE and REM fits are repeated from several starting values, with the grid spacing halved and with the τ² bound reduced tenfold, and the changes in log likelihood, trend parameters and τ² are checked automatically. [1–3]",
  "",
  "## Cutoff selection",
  "",
  "For each admissible k, separate hierarchical marginal likelihoods are fitted to the LE and REM PASs of the comparison, and the profile log likelihood is the sum of the two maximized objectives. Admissible cutoffs leave at least 200 PASs in each partition, with k at most half the comparison's PASs. The search evaluates the admissible endpoints and a grid of step 250, then refines at steps 25, 5 and 1 around the three best candidates of each stage (ties prefer smaller k); each fit starts from the parameters of the nearest evaluated cutoff and falls back to the pooled start, with multiple starts, if it does not converge. The selected k* is the largest converged objective on the evaluated grid: a grid-refined maximum, not a proven global optimum. No Wald statistic, p-value, fold change, discovery count or HBFSS score enters selection, and a counter of significance computations verifies this. Because the ranking, the likelihood and the tests use the same counts and the cutoff is chosen adaptively, the likelihood-ratio statistic and BIC are reported as descriptive measures of model fit, not as tests. EVS ranked by PC1 is the prespecified primary split.",
  "",
  "Pooled versus split. Whether the data support two dispersion populations at all is decided by a parametric bootstrap of the complete adaptive procedure (McLachlan's bootstrap of the likelihood-ratio statistic, which has no standard reference distribution here because k is chosen and the ranking is learned from the data). B data sets (19 in the development profile, 99 in the publication profile) are simulated from the fitted pooled model of the comparison: DESeq2's fitted means and dispersions drawn from the pooled trend and τ². Each is filtered, re-ranked (VST, per-condition PC1), its likelihoods are recomputed and its best split is found on the coarse grid. The statistic is the coarse-grid likelihood ratio, split versus pooled, compared with the observed coarse-grid statistic at the same resolution; p = (1 + #{LRT* ≥ LRT_obs})/(B + 1), computed from the finite bootstrap statistics, and EVS is selected when p ≤ 0.05; at least 80% of the bootstrap replicates must succeed, otherwise the decision is reported as undetermined and no adaptive result is given. The same rule is used for the observed comparisons and for every validation and external analysis, and every one of those analyses generates a new bootstrap distribution from its own counts: each simulated data set, stability subset, external comparison, null split and replication half is filtered, normalized, ranked and searched, and B bootstrap data sets are simulated from its own fitted pooled model and analysed in the same way. No bootstrap distribution is reused between analyses. The adaptive procedure reports EVS when selected and the pooled model otherwise; forced EVS is always reported as well, and the two are labelled separately in every table. [19, 28]",
  "",
  "## Differential PAS abundance",
  "",
  "For each comparison the Original (pooled hierarchical moderation, one population) and the two EVS partitions use the same raw counts, size factors, design and gene-wise dispersions. The hierarchical model is fitted to the PASs of each view; its trend and τ² are supplied to DESeq2 (dispersionFunction, estimateDispersionsMAP), and DESeq2's dispersion-outlier rule and Wald test are retained. [1, 2]",
  "",
  "One Cook's-distance exclusion list is derived from the Original (DESeq2's default rule) and applied to every view and model; excluded PASs enter neither the tests nor the empirical null, higher criticism or HBFSS. Benjamini–Hochberg (BH) adjustment is applied per comparison and test, across the Original's PASs or across the LE and REM PASs combined, without independent filtering; experiment-wide BH across the PASs of all four comparisons is reported in addition. The primary endpoint is the ordinary two-sided Wald test of β = 0 at BH < 0.10. The secondary strong-effect test is DESeq2's greaterAbs test of |β| > 1 at BH < 0.10. The secondary weak-effect test addresses the composite null hypothesis β = 0 or |β| ≥ 1 against 0 < |β| < 1 with the intersection–union test, whose p-value is the larger of the ordinary and lessAbs p-values, across all eligible PASs, followed by BH at 0.20; with five animals per group this test has little power, and a result of no Weak calls is reported as such. A prespecified minimum effect δ (δ < |β| < 1) is available but set to 0. No fold-change gate is added to any test, and no FDR claim is made for the union of the test families. Standard calls with |apeglm log₂FC| ≥ 1 are counted as a descriptive column only. [1, 4, 5, 15]",
  "",
  "Fold changes are shrunk with apeglm for interpretation and visualization, with the prior scale adapted within each view; shrunken values do not define rejections. [1, 6]",
  "",
  "## Exploratory higher criticism and HBFSS",
  "",
  "For finite Wald statistics of non-excluded PASs, fdrtool fits a normal empirical null within each view (cutoff method fndr). Higher criticism is maximized over the lowest 10% of the Original's ordered empirical p-values and must be positive; its threshold p_HC is applied to every view and model of the comparison. HBFSS = |apeglm log₂FC| × [−log₁₀(empirical p)], and a score above −log₁₀(p_HC) flags a PAS for prioritization. Empirical-null p-values, higher criticism and HBFSS are exploratory, have no claimed FDR control and are not part of the primary discoveries; their overlap with the Standard, Strong and Weak calls is reported. [7, 8, 10]",
  "",
  "## Model comparison",
  "",
  "A complete differential abundance analysis is run for: DESeq2 with its default settings (own size factors, parametric trend, its prior variance, Cook's flagging and independent filtering); edgeR robust quasi-likelihood (TMM normalization, robust empirical Bayes dispersions and QL F-tests) [11, 20, 21]; and, under the primary testing rules (the Original's Cook's list and HBFSS threshold, one BH per analysis, no independent filtering), harmonized conventional DESeq2 (parametric trend and DESeq2's prior variance); the pooled hierarchical model with TMM instead of median-of-ratios size factors; pooled hierarchical moderation (the Original); PC1-based EVS (primary); a split ranked by within-condition VST variance; a split ranked by each PAS's contribution to the PC1–PC2 subspace, λ₁v²ᵢ₁ + λ₂v²ᵢ₂; continuous-rank moderation, in which a centred natural cubic spline (3 df) of log(1 + dᵢ) shifts each PAS's prior centre (applied to DESeq2's dispFit before the MAP step) with one prior variance; random partitions (three by default), in which the union depths are permuted; the adaptive procedure (pooled versus EVS, above); and EVS at the cutoff of the published analysis, k = 5,000 per condition (Leading Edge = union of the two top-5,000 sets), with the current ranking, moderation and tests, so that only the cutoff differs from PC1 EVS at k*. The published analysis is also reproduced as published: PASs with at least one read in the comparison; PCA (centred, unscaled) on DESeq2 median-of-ratios normalized counts within each condition; the 5,000 PASs with the largest |PC1 loading| per condition, their union forming the Leading Edge and all other PASs the Remainder; standard DESeq2 (own size factors and dispersions) run separately on each subset with BH within each subset; Standard = BH < 0.10 and |apeglm log₂FC| ≥ 1, Strong = greaterAbs(1) BH < 0.10 and |apeglm log₂FC| ≥ 1, as defined in the published analysis (its Weak and HBFSS calls are not reproduced). Its Standard calls, Leading Edge size and overlap with PC1 EVS are reported side by side (the published Standard definition includes a fold-change gate that the current Standard test does not). Every split model is searched for its own k with the same likelihood, admissibility rules and grid resolution as PC1 EVS. Discovery counts and BIC are descriptive; whether a model controls the FDR and gains power is judged from the simulations. [1–3, 5]",
  "",
  "## Sensitivity analyses",
  "",
  "For each comparison: (i) the Original and EVS at k* with DESeq2's 0.25 floor on the MAP prior variance; and (ii) independent hypothesis weighting (IHW) of the Original's p-values with the mean of normalized counts as covariate, a covariate that is independent of the p-values under the null, run when the IHW package is installed. [5, 9]",
  "",
  "## Validation of error rates, power and stability",
  "",
  "Every replicate below re-runs the complete adaptive analysis on its own counts: filter, normalization, VST, PC1 ranking, profiled likelihood, k*, the pooled-versus-split decision (with a new bootstrap generated from that replicate, above), moderation, testing and BH. No observed-data ranking, cutoff or dispersion is reused, and the pipeline verifies that the observed-data caches are unchanged after validation.",
  "",
  "Label permutation of the observed data is not used: permuting animals within a comparison dilutes but does not remove real treatment effects, so it provides no true null. Error rates are estimated from simulations with known truth and from real-data null splits of an independent dataset (External validation).",
  "",
  "Negative binomial simulation. Counts are simulated from one comparison's fitted parameters (RT4 vs ZT10 by default): PAS means, size factors, the fitted pooled or LE/REM trends and prior variances, and PAS dependence through a sample-level effect shared by the PASs of a gene (SD 0.10 on the log scale). Seven scenarios: a complete null and effects (2% of PASs at each of |log₂FC| = 0.5, 1 and 2) under homogeneous dispersion; null and effects with two dispersion populations as fitted (the only scenarios generated from the EVS populations); effects with continuous heterogeneity (a PAS-level shift of log dispersion, SD 0.8); and null and effects with PC1 unrelated to dispersion (two populations assigned at random with respect to the ranking, plus a latent sample factor on 10% of PASs that dominates PC1). All scenarios use the base comparison's design (five animals per group) and sequencing depth. All models are applied to every data set: DESeq2 default, harmonized DESeq2, edgeR robust QL, pooled hierarchical (median-of-ratios and TMM), PC1 EVS, variance split, random partition and the adaptive procedure (plus the PC1–PC2 split and continuous-rank moderation in the publication profile). Reported per model and scenario as the Monte Carlo mean with a 95% interval (mean ± 1.96 SE; for the empirical FDR an empirical-Bernstein interval, below): empirical FDR (mean false discovery proportion, zero when there are no discoveries), the probability of any false discovery, power overall and by effect size, false-positive rate, specificity, precision, false non-discovery rate, bias and root mean squared error of the final dispersions (log scale) and of the log₂ fold changes, coverage of Wald 95% intervals, the empirical FDR and power of the Strong and Weak calls against their own truths and of the exploratory HBFSS flags against all non-null PASs, runtime, memory, failures and how often the adaptive procedure selects EVS. Because the false discovery proportion is bounded and is often zero in every replicate, its interval is the empirical-Bernstein 95% interval [31], which keeps a non-zero width in that case, instead of a normal approximation. With a prespecified tolerance of 2.5 percentage points, a method's FDR is classified as acceptable when the upper bound is at most 12.5%, unacceptable when the lower bound exceeds 12.5%, and inconclusive otherwise; with the development settings almost every result is inconclusive. A simulated data set whose analysis fails is re-analysed up to twice; a replicate that still fails is counted as failed for every model (attempted, successful and failed numbers are reported), never as zero discoveries. EVS is judged to gain only if its power is higher while its FDR is classified as acceptable. The design and reporting follow Morris et al. [27]",
  "",
  "Cutoff stability. Biological replicates are resampled within each arm by leaving one animal out per arm (25 subsets per comparison). Resampling with replacement is available but duplicates animals, which DESeq2 would treat as independent replicates. For each subset, k*, the PC1 ranking, LE membership and size, the trends and prior variances, the pooled-versus-split decision and the Standard discoveries are recomputed; reported are the median k* and its 2.5–97.5% range, LE size, the LE Jaccard similarity and the Spearman correlation of union depths with the full data, the variability of τ² and of the trend exponents, the frequency with which EVS is selected, and the proportion of full-data discoveries reproduced. No cutoff is chosen retrospectively.",
  "",
  "Run profiles. A development profile (19 bootstrap replicates, 2 simulations per scenario of 3,000 PASs, 3 stability subsets, two filter thresholds) checks that every stage runs and reveals gross departures only. The publication profile (99 bootstrap replicates, 1,000 simulations per scenario of 10,000 PASs, 25 stability subsets, three filter thresholds) provides the precision required for claims of FDR control and power; only its results support such claims.",
  "",
  "## External validation",
  "",
  "The frozen procedure (filter rule, ranking, likelihood, cutoff search, pooled-versus-split bootstrap, moderation and tests; nothing changed after these data were seen) is applied to two independent public 3′-end datasets chosen by prespecified criteria: public, PAS-level raw integer counts, independent biological replicates, at least three per group and a defined two-group contrast. Each analysis is compared with DESeq2 default, edgeR robust QL and the pooled hierarchical model, and every adaptive analysis generates its own bootstrap. Before analysis, PAS identifiers are checked for duplicates (identifiers are never altered), and sample identities, group assignments and raw count integrity (non-negative integers, no missing values) are verified; a failed check stops that dataset's analysis.",
  "",
  "GEO GSE100349 (rat WTTS-seq, the same 3′-end protocol and species as the primary study; Control versus OB, four animals per group, PAS table as deposited; one unlabelled total row removed) has no known truth. Reported are each method's discoveries, k* and the pooled-versus-split decision.",
  "",
  "GEO GSE151173 (mouse liver 3′ mRNA-seq every 4 h, ZT02–ZT22, six animals per time point; only the wild-type animals are used, each sampled once) provides two benchmarks. A real-data null: at each time point the six animals are split three versus three (10 distinct splits per time point, 60 in all); genotype, time and tissue are identical, so every discovery is false and no simulation assumption is involved. Replication with independent animals: for ZT14 versus ZT02 (antiphase, prespecified) the animals are divided into two disjoint halves of three plus three; each method is applied to half A, and half B supplies a reference defined independently of EVS, the PASs called (BH < 0.10) by both DESeq2 default and edgeR robust QL. Reported are the share of the reference recovered in half A and the share of half-A discoveries with the same sign and p < 0.05 in half B under both conventional methods. Because the reference is built from conventional methods it can favour them; it is not ground truth, and observed discovery comparisons are not validated power.",
  "",
  "## Time-course model",
  "",
  "The factorial design is analysed jointly in one DESeq2 model across all 40 samples, ~ time + condition + time:condition, with experiment-wide size factors and the hierarchical dispersion trend (one population). Likelihood-ratio tests against ~ time (a treatment effect at any recovery time) and against ~ time + condition (a treatment effect that changes with recovery time) are adjusted by BH within each test, and the treatment log₂ fold change and Wald test at each recovery time are obtained from the condition and interaction coefficients. This is a separate biological analysis; it is not used to choose any EVS setting. [1, 16]",
  "",
  "## Alternative polyadenylation",
  "",
  "APA is analysed separately from PAS abundance, for every gene with at least two eligible PASs, without requiring any abundance change. Gene units are defined by gene symbol, chromosome and strand, and PASs are ordered along the direction of transcription by their peak coordinate. edgeR quasi-likelihood models (robust dispersion) are fitted to raw counts with offsets from the shared size factors. (i) Usage change: diffSpliceDGE tests each PAS against the other PASs of its gene; gene-level p-values are Simes combinations with BH across genes, and, when stageR is installed, a stage-wise analysis screens genes and confirms individual PASs with overall FDR control at 0.10. (ii) Relative usage: the share of each gene's normalized reads at each PAS is computed per animal, with the absolute change in percentage points between conditions and its Welch 95% interval; a positional index, the usage-weighted mean relative position of the gene's PASs (0 most proximal, 1 most distal), captures switching among more than two PASs and is compared between conditions by Welch's t-test (with a 95% interval for the change) with BH across genes. (iii) Direction: the two most used PASs (mean normalized count ≥ 5 each) form a proximal/distal pair, and RED is the treatment log₂ fold change of the distal PAS minus that of the proximal PAS, tested by the gene-level diffSpliceDGE test of the pair with BH (FDR < 0.10); positive values are distal shifts and negative values proximal shifts. When a rat transcript annotation (GTF) of the same genome build is supplied, each PAS is placed in a 3′UTR (from the end of a transcript's CDS to the transcript end), elsewhere in an annotated gene, or outside annotated genes, annotated genes provide stable gene identifiers, and only shifts between two 3′UTR PASs of the same gene are labelled 3′UTR lengthening or shortening; otherwise all shifts are reported as distal/proximal switching. Whether the pair's PASs lie in the EVS Leading Edge is reported. [11, 12, 20, 22, 23]",
  "",
  "## 3′aTWAS enrichment",
  "",
  "3′aTWAS genes are mapped from human to rat with database-supported babelgene orthologs (at least one supporting source); symbol similarity alone is not orthology. Because 3′aTWAS identifies genetically regulated APA, the primary endpoint is differential APA (gene-level usage change); differential PAS abundance (any Standard call in the gene) is a secondary endpoint, analysed for the Original and EVS on one identical background (ortholog-mapped genes with a finite p-value in both). Gene-level enrichment uses Firth's penalized-likelihood logistic regression, which gives finite estimates for sparse and zero overlaps, adjusting for log(1 + number of tested PASs), log(1 + mean normalized count) and log(number of supporting ortholog databases), with profile penalized-likelihood 95% intervals (a failed or non-finite Firth fit is reported as not estimable, never replaced by ordinary logistic regression); a sensitivity analysis restricts to one-to-one orthologs. A threshold-free competitive test (limma cameraPR) compares the gene statistic z = Φ⁻¹(1 − p/2) of 3′aTWAS genes with all other background genes, with p the gene-level p-value (Simes over the gene's PASs for abundance; the Simes usage p-value for APA). All comparisons are reported, including zero overlaps, with BH across rows within each method. TWAS results are not used to choose any EVS setting. Overlap is corroboration, not causal evidence. [13, 14, 24–26]",
  "",
  "## Figures",
  "",
  "Figure 1, EVS: per-condition PC1 ranking and the selected cutoff (gene-wise dispersions with the LE, REM and pooled trends at k*; profile likelihood over k). Figure 2, dispersion moderation (DESeq2 dispersion plots, trend exponent, prior variance used). Figure 3, empirical FDR and power of every method across the simulation scenarios. Figure 4, differential PAS abundance and effect sizes (discoveries of every model, apeglm fold changes of the discoveries, discoveries gained and lost by EVS relative to the pooled model). Figure 5, stability and model selection (k* across resamples, the pooled-versus-split bootstrap). Figure 6, alternative polyadenylation. Figure 7, 3′aTWAS enrichment (Firth odds ratios with 95% intervals; non-estimable rows listed). Figure 8, the time-course model. Supplementary Figure S1, model fit and sensitivity analyses (descriptive BIC, discoveries of every model, prior-floor and IHW analyses); S2, APA detail; S3, further simulation endpoints (false non-discovery rate, dispersion and fold-change accuracy, interval coverage); S4, abundance filter sensitivity; S5, external validation (discoveries, false discoveries on real null data, replication with independent animals). Per comparison: volcano plots, PCA, PC1 evidence and 3′aTWAS overlap.",
  "",
  "## Reproducibility and limitations",
  "",
  "The run writes into a staging directory that replaces the previous output only after the run completes, so a failed run leaves earlier results intact; the result is one output directory and one ZIP archive of it, and raw data, annotations and code are never modified. Simulation results are written once, as aggregate tables, at the end of the simulation stage. Random numbers use L'Ecuyer-CMRG streams with fixed seeds. Unit tests on synthetic data (analytic gradients, the profile score equations, the filter, partition disjointness, joint BH, Simes, recovery of a known two-population split and determinism) are run with SEQUENCE_SELFTEST=1 and written with the synthetic dataset. The pipeline records the cutoff profile, the pooled-versus-split bootstrap, fitted coefficients and prior variances, the model comparison, the sensitivity and validation analyses, full PAS-level results, input MD5 checksums, sessionInfo(), an renv lockfile (renv is a required package; renv.lock: R, Bioconductor and every package with its version and installation source, restorable with renv::restore()), settings and seeds, stage runtimes and memory, and a table of automated validation checks. EVS selects its partition from the tested counts, so its nominal FDR is a claim that rests on the empirical validation, whose precision depends on the number of simulations and external null splits analysed. The true-null splits of the external liver data are the real-data null benchmark. The simulations use parameters fitted to one comparison. External validation uses one rat dataset with only four animals per group and unknown truth, and one mouse liver dataset in which the true-null splits have three animals per group; neither is a sleep-deprivation study, so they test the method's error control and reproducibility rather than the biological findings. The time-course model uses one dispersion population. Without a transcript annotation, APA is reported as distal/proximal switching. The nominal 10% and 20% FDR levels are prespecified, not optimized, and no discovery count or TWAS result was used to tune any setting.",
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
  "11. Chen Y, Lun ATL, Smyth GK. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. F1000Research. 2016;5:1438. https://doi.org/10.12688/f1000research.8987.2",
  "12. Van den Berge K, Soneson C, Robinson MD, Clement L. stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology. 2017;18:151. https://doi.org/10.1186/s13059-017-1277-0",
  "13. Alternative polyadenylation transcriptome-wide association study identifies APA-linked susceptibility genes in brain disorders. Nature Communications. 2023;14. https://doi.org/10.1038/s41467-023-36311-8",
  "14. babelgene: gene orthologs for model organisms in a tidy data format. https://igordot.github.io/babelgene/",
  "15. Berger RL, Hsu JC. Bioequivalence trials, intersection–union tests and equivalence confidence sets. Statistical Science. 1996;11:283–319. https://doi.org/10.1214/ss/1032280304",
  "16. Love MI, Anders S, Huber W. Analyzing RNA-seq data with DESeq2 (Bioconductor vignette; likelihood ratio tests and interactions). https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html",
  "17. Cox DR, Reid N. Parameter orthogonality and approximate conditional inference. Journal of the Royal Statistical Society B. 1987;49:1–39. https://doi.org/10.1111/j.2517-6161.1987.tb01422.x",
  "18. Robinson MD, Smyth GK. Moderated statistical tests for assessing differences in tag abundance. Bioinformatics. 2007;23:2881–2887. https://doi.org/10.1093/bioinformatics/btm453",
  "19. Phipson B, Smyth GK. Permutation p-values should never be zero: calculating exact p-values when permutations are randomly drawn. Statistical Applications in Genetics and Molecular Biology. 2010;9:39. https://doi.org/10.2202/1544-6115.1585",
  "20. Lund SP, Nettleton D, McCarthy DJ, Smyth GK. Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates. Statistical Applications in Genetics and Molecular Biology. 2012;11(5). https://doi.org/10.1515/1544-6115.1826",
  "21. Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology. 2010;11:R25. https://doi.org/10.1186/gb-2010-11-3-r25",
  "22. Hoque M, Ji Z, Zheng D, Luo W, Li W, You B, Park JY, Yehia G, Tian B. Analysis of alternative cleavage and polyadenylation by 3′ region extraction and deep sequencing. Nature Methods. 2013;10:133–139. https://doi.org/10.1038/nmeth.2288",
  "23. Simes RJ. An improved Bonferroni procedure for multiple tests of significance. Biometrika. 1986;73:751–754. https://doi.org/10.1093/biomet/73.3.751",
  "24. Firth D. Bias reduction of maximum likelihood estimates. Biometrika. 1993;80:27–38. https://doi.org/10.1093/biomet/80.1.27",
  "25. Heinze G, Schemper M. A solution to the problem of separation in logistic regression. Statistics in Medicine. 2002;21:2409–2419. https://doi.org/10.1002/sim.1047",
  "26. Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic Acids Research. 2012;40:e133. https://doi.org/10.1093/nar/gks461",
  "27. Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. Statistics in Medicine. 2019;38:2074–2102. https://doi.org/10.1002/sim.8086",
  "28. McLachlan GJ. On bootstrapping the likelihood ratio test statistic for the number of components in a normal mixture. Journal of the Royal Statistical Society Series C (Applied Statistics). 1987;36:318–324. https://doi.org/10.2307/2347790",
  "29. NCBI Gene Expression Omnibus GSE100349. WTTS-seq expression table, rat, Control and OB. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100349",
  "30. NCBI Gene Expression Omnibus GSE151173. Investigation of diurnal polyadenylation site usage reveals differential regulation in the transcription of gene isoforms (mouse liver 3′ mRNA-seq). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151173",
  "31. Maurer A, Pontil M. Empirical Bernstein bounds and sample variance penalization. Proceedings of the 22nd Conference on Learning Theory (COLT). 2009. https://arxiv.org/abs/0907.3740"
)

PIPELINE_BUILD <- "SEQUENCE_EVS_MANUSCRIPT_HIERARCHICAL_V9_FINAL"

# Run profile. development (default): small replicate numbers, enough to see
# gross departures and to check that every stage runs; publication: the
# replicate numbers the manuscript's FDR and power claims require
# (SEQUENCE_PROFILE=publication). Each number can be overridden individually.
env_int <- function(name, default) { v <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (is.na(v) || v < 0L) default else v }
RUN_PROFILE <- match.arg(Sys.getenv("SEQUENCE_PROFILE", "development"), c("development", "publication"))
PUB <- RUN_PROFILE == "publication"


# =============================================================================
# CONSTANTS AND SHARED FIGURE STYLE
# =============================================================================

COMPARISON_NAMES <- c("RT0_ZT6", "RT2_ZT8", "RT4_ZT10", "RT8_ZT14")
METHOD_SLUGS <- c("EVS k*" = "EVS_k")
CUTOFF_METHOD_LABEL <- names(METHOD_SLUGS)
CUTOFF_METHOD_BASIS <- paste0("Comparison-specific VST (ten samples); per-condition |PC1| union ranking; ",
 "grid-refined maximum of the combined hierarchical marginal adjusted likelihood of the two partitions ",
 "(design ~ condition); no p-value, fold change or discovery count used")


# EVS PCA input: the comparison's own DESeq2 VST (see EIGENVECTOR SPLITTING).
PCA_METHOD_LEVELS <- "VST"
PCA_METHOD_LABELS <- c(VST = "comparison-specific DESeq2 VST (ten samples, blind = TRUE)")
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

  # Output handling: the run writes into a staging directory. Only after it
  # finishes successfully is the previous final output (directory and ZIP,
  # identified by its marker file) replaced, so a failed run never destroys
  # earlier results and a successful run leaves exactly one output directory
  # and one ZIP. Earlier timestamped "_previous_" copies made by older versions
  # of this pipeline are removed at that point. Raw data, TWAS inputs,
  # annotations and code are never touched.
  final_root <- Sys.getenv("SEQUENCE_MULTI_ROOT",
                           unset = file.path(repo_guess, "exports", "sequence_evs_final"))
  final_root <- normalizePath(final_root, winslash = "/", mustWork = FALSE)
  if (identical(Sys.getenv("SEQUENCE_SELFTEST", "0"), "1")) final_root <- paste0(final_root, "_selftest")
  if (!grepl("/exports/", final_root) && !nzchar(Sys.getenv("SEQUENCE_ALLOW_OUTPUT_ANYWHERE")))
    stop("The output directory must lie under an 'exports' folder: ", final_root)
  OUTPUT_MARKER <- ".sequence_evs_output"
  multi_root <- paste0(final_root, "_staging")
  if (dir.exists(multi_root)) unlink(multi_root, recursive = TRUE, force = TRUE)
  dir.create(multi_root, recursive = TRUE, showWarnings = FALSE)
  writeLines("SEQUENCE EVS output directory", file.path(multi_root, OUTPUT_MARKER))

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
  if (!identical(as.integer(status), 0L))
    stop("The analysis run failed (exit status ", status, "). Previous final results are unchanged; ",
         "the incomplete run is in ", multi_root)

  # Replace the previous final output only now.
  if (dir.exists(final_root)) {
    if (!file.exists(file.path(final_root, OUTPUT_MARKER)) && length(list.files(final_root, all.files = TRUE, no.. = TRUE)))
      stop("Refusing to replace ", final_root, ": it is not a SEQUENCE EVS output directory (no marker file).")
    unlink(final_root, recursive = TRUE, force = TRUE)
  }
  old <- list.files(dirname(final_root), full.names = TRUE,
                    pattern = paste0("^", basename(final_root), "_previous_[0-9]{8}_[0-9]{6}$"))
  for (d in old) unlink(d, recursive = TRUE, force = TRUE)
  # Superseded outputs of earlier versions of this pipeline (same exports folder).
  for (legacy in c("sequence_evs_hierarchical_candidate", "sequence_evs_manuscript")) {
    hits <- list.files(dirname(final_root), full.names = TRUE,
                       pattern = paste0("^", legacy, "(_previous_[0-9]{8}_[0-9]{6})?(\\.zip)?$"))
    for (d in hits) unlink(d, recursive = TRUE, force = TRUE)
  }
  if (!file.rename(multi_root, final_root)) stop("Could not move the staging directory to ", final_root)
  multi_root <- final_root
  final_zip <- file.path(dirname(multi_root), paste0(basename(multi_root), ".zip"))
  if (file.exists(final_zip)) unlink(final_zip, force = TRUE)
  oldwd <- getwd()
  setwd(dirname(multi_root))
  zip_rel <- list.files(multi_root, recursive = TRUE)
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
  "limma",
  "GenomicRanges",
  "IRanges",
  "SummarizedExperiment",
  "splines",
  "parallel",
  "logistf",
  "renv"
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
# Prespecified, label-independent abundance filter: at least 10 reads in at
# least 5 of the comparison's 10 samples (5 = the size of one group, so a PAS
# expressed in only one condition can still be eligible). PASs below this have
# too little information for a dispersion estimate (Bourgon et al. 2010; Chen
# et al. 2016). Sensitivity to other thresholds is assessed in the validation
# stage (FILTER_SENSITIVITY). The same eligible population is used by every
# model of a comparison.
MIN_PAS_COUNT <- as.integer(Sys.getenv("SEQUENCE_MIN_PAS_COUNT", "10"))
MIN_PAS_SAMPLES <- as.integer(Sys.getenv("SEQUENCE_MIN_PAS_SAMPLES", "5"))
stopifnot(is.finite(MIN_PAS_COUNT), MIN_PAS_COUNT >= 1L, is.finite(MIN_PAS_SAMPLES), MIN_PAS_SAMPLES >= 1L)
# Purpose: Eligible pas.
eligible_pas <- function(x, min_count = MIN_PAS_COUNT, min_samples = MIN_PAS_SAMPLES)
  rowSums(x >= min_count) >= min_samples
# Designs with other group sizes (external datasets) scale the sample requirement with the group size (5 of 5
# here), so the rule keeps its meaning: >= 10 reads in as many samples as one
# group holds.
filter_rule <- function(coldata) {
  n <- min(table(droplevels(factor(coldata$condition))))
  c(MIN_PAS_COUNT, max(1L, as.integer(round(MIN_PAS_SAMPLES * n / 5))))
}

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
# Weak: composite null (beta = 0 or |beta| >= 1) tested by the intersection-union test
# max(ordinary p, lessAbs p) across all eligible PASs, then BH < 0.20.
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

# Parallel-safe random number streams (reproducible with mclapply).
RNGkind("L'Ecuyer-CMRG")

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
      "[%s] HC threshold: alpha0=%.3f, HCmax=%.6f, HCp=%.8g, search=%d/%d empirical p-values.",
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

# Per comparison: the comparison's PCA input (its own VST) and the
# per-condition PC1 fits, cached (keyed "<comparison>|<method>").
PCA_RANK_CACHE <- list()

# The matrix the EVS PCA is run on: the comparison's own VST (ten samples,
# shared size factors), computed once in the comparison cache, the same values
# the cutoff rule ranked. Used only to rank PASs; no transformed matrix is ever
# passed to DESeq2.
make_rank_matrix_for_track <- function(count_matrix, coldata, track_key, method = EVS_PCA_METHOD,
                                       size_factors = NULL, cache_key = NULL) {
  track_key <- match.arg(track_key, "normalized_evs")
  method <- match.arg(method, PCA_METHOD_LEVELS)
  key <- if (is.null(cache_key)) NULL else paste(cache_key, method, sep = "|")
  hit <- if (is.null(key)) NULL else PCA_RANK_CACHE[[key]]
  if (!is.null(hit) && identical(hit$ids, rownames(count_matrix)) && identical(hit$samples, colnames(count_matrix))) {
    return(hit$rank_obj)
  }
  full <- if (is.null(cache_key)) NULL else comparison_cache[[cache_key]]$vst
  if (is.null(full)) stop("The comparison VST is missing; prepare the comparison cache first.")
  miss <- c(setdiff(rownames(count_matrix), rownames(full)), setdiff(colnames(count_matrix), colnames(full)))
  if (length(miss)) stop("The comparison ", method, " matrix lacks PASs or samples.")
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

  # The comparison's PCA input (its own VST) and its PC1 rankings (cached).
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

# ---- Per-PAS Cox-Reid adjusted dispersion likelihood -------------------------
# Purpose. For every PAS, the Cox-Reid adjusted log likelihood of its dispersion
# on a fixed grid of log dispersions: the input of the hierarchical marginal
# likelihood (cutoff search and final moderation use the same curves).
# Rationale. The adjustment -1/2 log det(X'WX) removes the first-order bias of
# the dispersion MLE caused by estimating the means (Cox & Reid 1987; McCarthy
# et al. 2012; DESeq2 uses the same adjusted likelihood).
# Assumptions. Negative binomial counts, Var = mu + alpha mu^2; a saturated
# categorical design (one mean per group: ~ condition for a comparison, eight
# arms for the time course), for which det(X'WX) is the product of the group
# weight sums W_g = sum_j mu_j / (1 + alpha mu_j) (up to a constant).
# Profiled (default) vs fixed-mean. With profile = TRUE the group means are
# re-estimated by maximum likelihood at every grid dispersion (Newton steps on
# the log mean; the NB log likelihood is concave in the log mean), giving the
# properly profiled adjusted likelihood. profile = FALSE holds the means fixed
# at DESeq2's fitted values (a plug-in approximation); it is computed only to
# report how much the approximation would change the dispersion estimates and
# k* (cutoff_numerical_checks). SEQUENCE_PROFILED_LIKELIHOOD=0 reverts.
# Output. A PAS x grid matrix of relative likelihoods (maximum 1 per PAS), with
# the informative PASs as row names.
DISP_PROFILE <- !identical(Sys.getenv("SEQUENCE_PROFILED_LIKELIHOOD", "1"), "0")

design_groups <- function(dds) {
  X <- stats::model.matrix(DESeq2::design(dds), as.data.frame(SummarizedExperiment::colData(dds)))
  if (qr(X)$rank != ncol(X) || nrow(X) <= ncol(X)) stop("Invalid dispersion model design.")
  grp <- factor(apply(X, 1L, paste, collapse = ":"))
  if (nlevels(grp) != ncol(X)) stop("The adjusted likelihood requires a saturated categorical design.")
  grp
}

# Maximum-likelihood group means of every PAS at a fixed dispersion a, by
# Newton steps on the log mean. start: the log means at the previous grid
# dispersion (warm start); NULL starts from the Poisson MLE.
profile_group_means <- function(y, s, grp, a, start = NULL, iter = 30L) {
  mu <- matrix(0, nrow(y), ncol(y))
  E <- matrix(0, nrow(y), nlevels(grp))
  for (gi in seq_along(levels(grp))) {
    g <- levels(grp)[gi]
    j <- which(grp == g); yg <- y[, j, drop = FALSE]; sg <- s[j]
    eta <- if (is.null(start)) log(pmax(rowSums(yg) / sum(sg), 1e-8)) else start[, gi]
    for (it in seq_len(iter)) {
      M <- outer(exp(eta), sg)
      f <- rowSums((yg - M) / (1 + a * M))
      fp <- -rowSums(M * (1 + a * yg) / (1 + a * M)^2)
      step <- f / fp
      step[!is.finite(step)] <- 0
      step <- pmax(pmin(step, 1), -1)            # damped Newton step
      eta <- pmin(pmax(eta - step, log(1e-8)), log(1e9))
      if (max(abs(step)) < 1e-8) break
    }
    mu[, j] <- outer(exp(eta), sg)
    E[, gi] <- eta
  }
  list(mu = mu, eta = E)
}

make_dispersion_likelihood <- function(dds, profile = DISP_PROFILE) {
  grp <- design_groups(dds)
  if ("weights" %in% SummarizedExperiment::assayNames(dds)) stop("Weighted counts are unsupported by the custom likelihood.")
  ids <- rownames(dds)
  y <- DESeq2::counts(dds)
  s <- DESeq2::sizeFactors(dds)
  if (is.null(s)) stop("Size factors are required for the adjusted likelihood.")
  mu0 <- SummarizedExperiment::assays(dds)[["mu"]]
  if (is.null(mu0)) stop("DESeq2 did not retain the fitted means; check the installed version.")
  mc <- S4Vectors::mcols(dds)
  use <- !mc$allZero & is.finite(mc$baseMean) & mc$baseMean > 0 &
    apply(mu0, 1L, function(z) all(is.finite(z) & z > 0))
  use[is.na(use)] <- FALSE
  y <- y[use, , drop = FALSE]; mu0 <- mu0[use, , drop = FALSE]
  if (!nrow(y)) stop("No informative PASs for the adjusted likelihood.")
  ll <- matrix(NA_real_, nrow(y), length(DISP_LOG_GRID), dimnames = list(ids[use], NULL))
  eta <- NULL
  for (h in seq_along(DISP_LOG_GRID)) {
    a <- exp(DISP_LOG_GRID[h])
    mu <- if (profile) { pm <- profile_group_means(y, s, grp, a, eta); eta <- pm$eta; pm$mu } else mu0
    z <- rowSums(matrix(stats::dnbinom(as.vector(y), mu = as.vector(mu), size = 1/a, log = TRUE), nrow(y)))
    w <- mu/(1 + a*mu)
    for (g in levels(grp)) z <- z - 0.5 * log(pmax(rowSums(w[, grp == g, drop = FALSE]), 1e-300))
    ll[, h] <- z
  }
  maxll <- apply(ll, 1L, max)
  if (any(!is.finite(maxll))) stop("Nonfinite adjusted likelihood curve.")
  exp(ll - maxll)
}

# ---- Hierarchical marginal likelihood ------------------------------------------
# Purpose. The log marginal likelihood of a population of PASs under the
# hierarchical prior log alpha_i ~ N(log(a0 + a1 mean_i^-b) + basis_i' c, tau^2),
# and its analytic gradient in q = (log a0, log a1, b, log tau^2, c).
# Rationale. Integrating each PAS's adjusted likelihood over the prior weights
# PASs by how informative their likelihoods are and retains boundary estimates
# (empirical Bayes moderation; Robinson & Smyth 2007; Love et al. 2014).
# Numerical integration. Normal CDF differences assign the prior mass to the
# grid cells (cell edges midway between grid points, outermost cells open), so
# the cell probabilities sum to one. Accuracy is checked by halving the grid
# step; the analytic gradient is checked against central differences.
marginal_objective <- function(q, lx, likelihood, basis) {
  d <- ncol(basis)
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
  list(value = -sum(log(den)), grad = -g)
}

# Largest relative difference between the analytic gradient and central
# differences (step 1e-5) at q.
gradient_check <- function(q, lx, likelihood, basis = NULL) {
  if (is.null(basis)) basis <- matrix(numeric(), nrow = length(lx), ncol = 0L)
  ga <- marginal_objective(q, lx, likelihood, basis)$grad
  h <- 1e-5
  gn <- vapply(seq_along(q), function(i) {
    e <- rep(0, length(q)); e[i] <- h
    (marginal_objective(q + e, lx, likelihood, basis)$value - marginal_objective(q - e, lx, likelihood, basis)$value) / (2 * h)
  }, numeric(1))
  max(abs(ga - gn)) / max(1, max(abs(gn)))
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
  hi <- unname(c(log(TREND_MAX_COEF), TREND_B_RANGE[2], log(DISP_PRIOR_VAR_MAX), rep(5, d)))
  cache <- new.env(parent = emptyenv())
  calc <- function(q) {
    if (!is.null(cache$q) && identical(q, cache$q)) return(cache$out)
    out <- marginal_objective(q, lx, likelihood, basis)
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
  par <- unname(o$par)
  list(converged=TRUE, loglik=-o$value, par=par, n=length(lx), tau2=exp(par[4]))
}

# Fit the identical hierarchical objective to the comparison counts before MAP
# shrinkage. The prior variance is estimated directly, not by subtracting a
# common trigamma sampling variance.
# basis (continuous-rank model only): a matrix with one row per PAS whose
# columns shift the prior centre (log scale); NULL for the power trend alone.
set_dispersion_trend <- function(dds, basis = NULL) {
  lk <- make_dispersion_likelihood(dds)
  means <- S4Vectors::mcols(dds)$baseMean[match(rownames(lk), rownames(dds))]
  B <- if (is.null(basis)) NULL else basis[rownames(lk), , drop = FALSE]
  f <- fit_marginal_trend(log(means), lk, basis = B, multistart=TRUE)
  if (!f$converged) stop("Hierarchical dispersion trend failed to converge.")
  p <- f$par
  trend <- make_trend_function(exp(p[1]), exp(p[2]), p[3])
  attr(trend, "basis_coef") <- if (length(p) > 4L) p[-(1:4)] else numeric()
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
estimate_dispersions_with_trend <- function(dds, size_factors = NULL, gene_dds = NULL, prior_floor = DISP_PRIOR_VAR_MIN,
                                            basis = NULL) {
  if (!is.null(gene_dds)) {
    if (!all(rownames(dds) %in% rownames(gene_dds))) stop("Precomputed gene-wise dispersions lack some PASs.")
    dds <- gene_dds[rownames(dds), ]
  } else {
    dds <- apply_size_factors(dds, size_factors)
    dds <- suppressMessages(DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE))
  }
  dds <- set_dispersion_trend(dds, basis)
  # Store the prior variance on the trend (as DESeq2 does internally), then
  # pass it to the MAP step so DESeq2 does not re-estimate it.
  prior_var <- dispersion_prior_variance(dds, prior_floor)
  f <- DESeq2::dispersionFunction(dds)
  attr(f, "dispPriorVar") <- prior_var
  attr(f, "dispPriorVar_floor") <- prior_floor
  dds <- DESeq2::`dispersionFunction<-`(dds, value = f, estimateVar = FALSE)
  # Continuous-rank model: DESeq2's MAP step centres each PAS's prior on its
  # dispFit, so the rank shift is applied to dispFit after the trend is set.
  bc <- attr(f, "basis_coef")
  if (!is.null(basis) && length(bc)) {
    off <- drop(basis[rownames(dds), , drop = FALSE] %*% bc)
    S4Vectors::mcols(dds)$dispFit <- S4Vectors::mcols(dds)$dispFit * exp(off)
  }
  suppressMessages(DESeq2::estimateDispersionsMAP(dds, dispPriorVar = prior_var, quiet = TRUE))
}

# Verification that DESeq2's MAP step uses the supplied prior variance (item 2):
# re-running the MAP step with the stored prior variance reproduces the final
# dispersions, and re-running it with four times that variance changes them.
verify_prior_variance_used <- function(dds) {
  f <- DESeq2::dispersionFunction(dds)
  pv <- attr(f, "dispPriorVar")
  if (is.null(pv) || !is.finite(pv)) return(list(ok = FALSE, detail = "no stored prior variance"))
  d0 <- S4Vectors::mcols(dds)$dispersion
  same <- suppressMessages(DESeq2::estimateDispersionsMAP(dds, dispPriorVar = pv, quiet = TRUE))
  wide <- suppressMessages(DESeq2::estimateDispersionsMAP(dds, dispPriorVar = 4 * pv, quiet = TRUE))
  d1 <- S4Vectors::mcols(same)$dispersion; d4 <- S4Vectors::mcols(wide)$dispersion
  ok_same <- isTRUE(all.equal(d0, d1, tolerance = 1e-6))
  ch <- mean(abs(log(d4) - log(d0)) > 1e-4, na.rm = TRUE)
  list(ok = ok_same && ch > 0, detail = sprintf("prior variance %.4f: same-variance rerun reproduces the dispersions (%s); 4x variance changes %.1f%% of them",
                                               pv, if (ok_same) "yes" else "no", 100 * ch))
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
# (model stage, fit_view_model), and its own significance threshold
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
                           prior_floor = DISP_PRIOR_VAR_MIN, basis = NULL) {
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
  dds <- if (identical(trend, "power")) estimate_dispersions_with_trend(dds, size_factors, gene_dds, prior_floor, basis) else
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
  res_df$dispersion_final <- S4Vectors::mcols(dds)$dispersion[match(res_df$feature_id, rownames(dds))]
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
#      statistics, so each view's empirical p-values are referenced to its own
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
    " selected by the in-sample split likelihood using each condition's PC1 eigenvector; ",
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
                                    original_cooks_cache[[comparison_name]], original_hc_p, analysis_results),
           error = function(e) warning(comparison_name, ": sensitivity analyses failed: ", conditionMessage(e), call. = FALSE))
  tryCatch(run_model_comparison(comparison_name, evs, coldata, annot_df, sf, cache$gene_dds,
                                original_cooks_cache[[comparison_name]], original_hc_p, analysis_results),
           error = function(e) warning(comparison_name, ": model comparison failed: ", conditionMessage(e), call. = FALSE))
  for (nm in c("raw_dataset", "leading_edge_dataset", "remainder_dataset")) {
    pv <- tryCatch(verify_prior_variance_used(analysis_results[[nm]]$dds),
                   error = function(e) list(ok = FALSE, detail = conditionMessage(e)))
    add_validation("DESeq2's MAP step uses the intended custom prior variance", paste(comparison_name, nm), pv$ok, pv$detail)
  }
  reg <- paper_registry[[registry_key(comparison_name, track_key)]]
  reg$analysis_results <- analysis_results
  reg$size_factors <- sf
  paper_registry[[registry_key(comparison_name, track_key)]] <<- reg
  validate_track_normalization(comparison_name, track_key, analysis_results, sf, evs, coldata)
  TRUE
}

# =============================================================================
# SENSITIVITY ANALYSES AND MODEL COMPARISON (per comparison, after the primary
# analysis)
# =============================================================================
# All under the primary testing rules: the Original's Cook's outlier list, one
# BH correction per analysis (across all views of the analysis), the
# Original's HBFSS threshold, no fold-change gate.
# Sensitivity analyses:
#   1. DESeq2's conventional 0.25 floor on the MAP prior variance (Original and
#      EVS at k*), to compare the estimated tau^2 with the DESeq2 default;
#   2. IHW (Ignatiadis et al. 2016) on the Original's p-values with the mean of
#      normalized counts as covariate. This covariate is independent of the
#      p-values under the null (Bourgon et al. 2010; the covariate used in the
#      DESeq2 documentation); the PC1 depth is not used because it is computed
#      from the same counts within each condition.
# Model comparison: complete differential abundance analysis with DESeq2
#   default settings, edgeR robust quasi-likelihood, DESeq2 harmonized to the
#   primary rules ("Conventional DESeq2": parametric trend and DESeq2 prior
#   variance, shared filter and size factors, Original's Cook's list, no
#   independent filtering), pooled hierarchical moderation (the Original) with
#   median-of-ratios and with TMM size factors, PC1 EVS, variance-ranked and
#   PC1-PC2-ranked splits (each at its own k*, same search), continuous-rank
#   moderation and random partitions (each at its own k*). The DE counts are
#   descriptive; whether a model controls the FDR and gains power is assessed
#   by the simulations (VALIDATION section), not by these counts or by BIC.
sensitivity_rows <- list()
model_comparison_rows <- list()
MODEL_FLAGS <- list()    # per comparison and model: per-PAS calls (experiment-wide BH, TWAS)

# Fit one analysis (a named list of count matrices, one per view) under the
# primary rules and return its per-view result tables (lean).
lean_analysis <- function(mats, coldata, label, sf, gene_dds, cooks, hc_p, annot_df,
                          prior_floor = DISP_PRIOR_VAR_MIN, trend = "power", basis = NULL) {
  fits <- Map(function(m, nm) fit_view_model(m, coldata, paste(label, nm), size_factors = sf, gene_dds = gene_dds,
                                             cooks_outliers = cooks, trend = trend, prior_floor = prior_floor,
                                             basis = basis),
              mats, names(mats))
  sig <- call_significance(fits, label, annot_df, lean = TRUE, hc_p_fixed = hc_p)
  lapply(sig$views, `[[`, "results")
}

# Leading Edge / Remainder matrices of a ranking at its cutoff (NULL if either
# partition is smaller than SPLIT_MIN_SET).
split_matrices <- function(y, depth, k) {
  ids <- rownames(y); le <- ids[depth[ids] <= k]; rem <- setdiff(ids, le)
  if (length(le) < SPLIT_MIN_SET || length(rem) < SPLIT_MIN_SET) return(NULL)
  list(leading_edge_dataset = y[le, , drop = FALSE], remainder_dataset = y[rem, , drop = FALSE])
}

# ---- External comparator models ---------------------------------------------
# Per-PAS calls in the format of discovery_flags for models fitted outside the
# primary machinery (Standard only; Strong, Weak and HBFSS are not defined).
flags_from_p <- function(feature_id, pvalue, lfc, se = NA_real_, padj = NULL, dispersion = NA_real_) {
  q <- if (is.null(padj)) { qq <- rep(NA_real_, length(pvalue)); ok <- is.finite(pvalue)
    qq[ok] <- stats::p.adjust(pvalue[ok], "BH"); qq } else padj
  data.frame(View = "all", feature_id = as.character(feature_id), pvalue = pvalue, padj_global = q,
             Standard = !is.na(q) & q < BH_FDR_STANDARD, Strong = NA, Weak = NA, HBFSS = NA,
             log2FoldChange = lfc, lfcSE = se, dispersion_final = dispersion, stringsAsFactors = FALSE)
}

# DESeq2 with its default settings (own size factors, parametric trend,
# DESeq2 prior variance, Cook's flagging and independent filtering in
# results()); Love et al. 2014.
fit_deseq2_default <- function(y, coldata) {
  dds <- DESeq2::DESeqDataSetFromMatrix(y, colData = coldata, design = make_design_formula(coldata))
  dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
  r <- DESeq2::results(dds, contrast = c("condition", "trt", "untrt"), alpha = BH_FDR_STANDARD)
  flags_from_p(rownames(r), r$pvalue, r$log2FoldChange, r$lfcSE, padj = r$padj,
               dispersion = DESeq2::dispersions(dds))
}

# edgeR robust quasi-likelihood (TMM normalization, robust empirical Bayes
# dispersion and QL F-test; Robinson & Oshlack 2010; Lund et al. 2012;
# Chen et al. 2016). Fold changes have no Wald standard error.
fit_edger_ql <- function(y, coldata) {
  X <- stats::model.matrix(make_design_formula(coldata), coldata)
  d <- edgeR::calcNormFactors(edgeR::DGEList(counts = y), method = "TMM")
  d <- edgeR::estimateDisp(d, X, robust = TRUE)
  f <- edgeR::glmQLFit(d, X, robust = TRUE)
  tt <- edgeR::topTags(edgeR::glmQLFTest(f, coef = ncol(X)), n = Inf, sort.by = "none")$table
  flags_from_p(rownames(tt), tt$PValue, tt$logFC, dispersion = d$tagwise.dispersion[match(rownames(tt), rownames(d))])
}

# TMM size factors (effective library sizes rescaled to geometric mean 1), for
# the pooled hierarchical model with TMM normalization.
tmm_size_factors <- function(y) {
  d <- edgeR::calcNormFactors(edgeR::DGEList(counts = y), method = "TMM")
  e <- d$samples$lib.size * d$samples$norm.factors
  stats::setNames(e / exp(mean(log(e))), colnames(y))
}

run_sensitivity_analyses <- function(comparison_name, evs, coldata, annot_df, sf, gene_dds, cooks, hc_p,
                                     analysis_results) {
  count_row <- function(label, res_list, note, k = NA_integer_, le = NA_integer_) {
    cnt <- count_discoveries(discovery_flags(res_list))
    data.frame(Comparison = comparison_name, Analysis = label, k = k, Leading_Edge = le,
               PAS_tested = cnt[["n_tested"]], Standard = cnt[["Standard"]], Strong = cnt[["Strong"]],
               Weak = cnt[["Weak"]], HBFSS = cnt[["HBFSS"]], Note = note)
  }
  k <- get_active_evs_cutoff(comparison_name)
  rows <- list(
    count_row(ORIGINAL_LABEL, list(raw_dataset = analysis_results$raw_dataset$results), "primary"),
    count_row("EVS (primary)", list(leading_edge_dataset = analysis_results$leading_edge_dataset$results,
                                    remainder_dataset = analysis_results$remainder_dataset$results), "primary",
              k, length(evs$leading_edge_ids)))
  rows[[3]] <- count_row(paste0(ORIGINAL_LABEL, ", DESeq2 0.25 prior floor"),
                         lean_analysis(list(raw_dataset = evs$raw_dataset), coldata,
                                       paste(comparison_name, "Original floor 0.25"), sf, gene_dds, cooks, hc_p,
                                       annot_df, prior_floor = DESEQ2_PRIOR_VAR_FLOOR),
                         "reference: MAP prior variance max(tau2, 0.25)")
  rows[[4]] <- count_row("EVS, DESeq2 0.25 prior floor",
                         lean_analysis(list(leading_edge_dataset = evs$leading_edge_dataset,
                                            remainder_dataset = evs$remainder_dataset), coldata,
                                       paste(comparison_name, "EVS floor 0.25"), sf, gene_dds, cooks, hc_p,
                                       annot_df, prior_floor = DESEQ2_PRIOR_VAR_FLOOR),
                         "reference: MAP prior variance max(tau2, 0.25)", k, length(evs$leading_edge_ids))
  o <- analysis_results$raw_dataset$results
  if (requireNamespace("IHW", quietly = TRUE)) {
    ok <- is.finite(o$pvalue) & is.finite(o$baseMean)
    ih <- IHW::ihw(o$pvalue[ok], o$baseMean[ok], alpha = BH_FDR_STANDARD)
    rows[[length(rows) + 1L]] <- data.frame(Comparison = comparison_name, Analysis = "IHW, Original p-values",
      k = NA_integer_, Leading_Edge = NA_integer_, PAS_tested = sum(ok), Standard = IHW::rejections(ih),
      Strong = NA_integer_, Weak = NA_integer_, HBFSS = NA_integer_,
      Note = "covariate: mean of normalized counts (null-independent); alpha = Standard BH level")
  } else {
    rows[[length(rows) + 1L]] <- data.frame(Comparison = comparison_name, Analysis = "IHW, Original p-values",
      k = NA_integer_, Leading_Edge = NA_integer_, PAS_tested = NA_integer_, Standard = NA_integer_,
      Strong = NA_integer_, Weak = NA_integer_, HBFSS = NA_integer_, Note = "not run: IHW package not installed")
  }
  sensitivity_rows[[comparison_name]] <<- dplyr::bind_rows(rows)
  invisible(TRUE)
}

# Complete DE analysis of every model of the comparison (item 4).
# Published analysis (SEQUENCE.R, build SEQUENCE_JOURNAL_READY_2026-08-22),
# reproduced as published so its results can be compared directly: PASs with
# at least one read in the comparison; PCA (centred, unscaled) on DESeq2
# median-of-ratios normalized counts within each condition; the 5,000 PASs with
# the largest |PC1 loading| per condition; Leading Edge = the union of the two
# top-5,000 sets (Joint + Disjoint), Remainder = all other PASs; standard DESeq2
# (own size factors and dispersions, default settings) run separately on each
# subset with BH within each subset; Standard = BH < 0.10 and |apeglm log2FC|
# >= 1, Strong = greaterAbs(1) BH < 0.10 and |apeglm log2FC| >= 1 (the
# published definitions). The published Weak and HBFSS calls depend on its
# empirical-null step and are not reproduced (NA).
PUBLISHED_TOP_N <- 5000L
fit_published_5k <- function(raw, coldata) {
  raw <- raw[, rownames(coldata), drop = FALSE]
  dds0 <- DESeq2::estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(raw, colData = coldata, design = make_design_formula(coldata)))
  nc <- DESeq2::counts(dds0, normalized = TRUE)
  top <- function(cond) {
    smp <- rownames(coldata)[coldata$condition == cond]
    ld <- abs(stats::prcomp(t(nc[, smp, drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 1)$rotation[, 1])
    names(ld)[order(ld, decreasing = TRUE)][seq_len(min(PUBLISHED_TOP_N, length(ld)))]
  }
  le <- union(top("trt"), top("untrt")); rem <- setdiff(rownames(raw), le)
  one <- function(ids, view) {
    dds <- DESeq2::DESeqDataSetFromMatrix(raw[ids, , drop = FALSE], colData = coldata, design = make_design_formula(coldata))
    dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
    dds <- suppressMessages(DESeq2::DESeq(dds, betaPrior = FALSE, quiet = TRUE))
    ct <- c("condition", "trt", "untrt")
    r <- DESeq2::results(dds, contrast = ct, alpha = BH_FDR_STANDARD, pAdjustMethod = "BH")
    rs <- DESeq2::results(dds, contrast = ct, lfcThreshold = lfc_boundary, altHypothesis = "greaterAbs",
                          alpha = BH_FDR_STRONG, pAdjustMethod = "BH")
    shr <- suppressMessages(DESeq2::lfcShrink(dds, coef = get_condition_coef(dds), type = "apeglm", quiet = TRUE))$log2FoldChange
    big <- is.finite(shr) & abs(shr) >= lfc_boundary
    data.frame(View = view, feature_id = rownames(r), pvalue = r$pvalue, padj_global = r$padj,
               Standard = !is.na(r$padj) & r$padj < BH_FDR_STANDARD & big,
               Strong = !is.na(rs$padj) & rs$padj < BH_FDR_STRONG & big, Weak = NA, HBFSS = NA,
               log2FoldChange = r$log2FoldChange, lfcSE = r$lfcSE, lfc_shrunk = shr,
               dispersion_final = DESeq2::dispersions(dds), stringsAsFactors = FALSE)
  }
  f <- rbind(one(le, "leading_edge_dataset"), one(rem, "remainder_dataset"))
  attr(f, "leading_edge") <- length(le)
  f
}

run_model_comparison <- function(comparison_name, evs, coldata, annot_df, sf, gene_dds, cooks, hc_p,
                                 analysis_results) {
  y <- evs$raw_dataset; ids <- rownames(y)
  lab <- function(x) paste(comparison_name, x)
  res <- list()
  res[["Conventional DESeq2"]] <- lean_analysis(list(raw_dataset = y), coldata, lab("conventional DESeq2"), sf,
                                                gene_dds, cooks, hc_p, annot_df, trend = "deseq2_default")
  res[["Pooled hierarchical"]] <- list(raw_dataset = analysis_results$raw_dataset$results)
  sf_tmm <- tmm_size_factors(y)
  res[["Pooled hierarchical, TMM"]] <- tryCatch(lean_analysis(list(raw_dataset = y), coldata, lab("pooled TMM"), sf_tmm,
                                                              comparison_gene_dds(y, coldata, sf_tmm), cooks, hc_p, annot_df),
                                                error = function(e) { warning(comparison_name, ": TMM model failed: ", conditionMessage(e), call. = FALSE); NULL })
  res[["PC1 EVS"]] <- list(leading_edge_dataset = analysis_results$leading_edge_dataset$results,
                           remainder_dataset = analysis_results$remainder_dataset$results)
  mk <- MODEL_K[[comparison_name]]; md <- MODEL_DEPTHS[[comparison_name]]
  label_of <- c(variance = "Variance split", pc12 = "PC1-PC2 split")
  for (m in names(mk)) {
    if (m == "pc1") next
    mats <- split_matrices(y, md[[m]], mk[[m]])
    if (is.null(mats)) next
    nm <- if (m %in% names(label_of)) label_of[[m]] else paste("Random partition", sub("^random", "", m))
    res[[nm]] <- tryCatch(lean_analysis(mats, coldata, lab(nm), sf, gene_dds, cooks, hc_p, annot_df),
                          error = function(e) { warning(comparison_name, ": ", nm, " failed: ", conditionMessage(e), call. = FALSE); NULL })
  }
  if (CUTOFF_COMPARATORS) {
    res[["Continuous rank"]] <- tryCatch(lean_analysis(list(raw_dataset = y), coldata, lab("continuous rank"), sf,
                                                       gene_dds, cooks, hc_p, annot_df,
                                                       basis = continuous_basis(md$pc1, ids)),
                                         error = function(e) { warning(comparison_name, ": continuous-rank model failed: ",
                                                                       conditionMessage(e), call. = FALSE); NULL })
  }
  # The paper's cutoff (k = 5,000 per condition) on the current ranking and
  # framework: only the cutoff differs from PC1 EVS.
  mats5 <- split_matrices(y, md$pc1, PUBLISHED_TOP_N)
  if (!is.null(mats5))
    res[["EVS, k = 5,000 (paper cutoff)"]] <- tryCatch(lean_analysis(mats5, coldata, lab("EVS k = 5000"), sf, gene_dds, cooks, hc_p, annot_df),
                                                       error = function(e) { warning(comparison_name, ": k = 5,000 model failed: ", conditionMessage(e), call. = FALSE); NULL })
  res <- Filter(Negate(is.null), res)
  cm <- CUTOFF$models[CUTOFF$models$Comparison == comparison_name, , drop = FALSE]
  flags <- lapply(res, discovery_flags)
  ext <- list("DESeq2 default" = tryCatch(fit_deseq2_default(y, coldata), error = function(e) NULL),
              "edgeR robust QL" = tryCatch(fit_edger_ql(y, coldata), error = function(e) NULL),
              "Published analysis (k = 5,000)" = tryCatch(fit_published_5k(comparison_raw_counts(comparison_name), coldata),
                                                          error = function(e) { warning(comparison_name, ": published analysis failed: ",
                                                                                        conditionMessage(e), call. = FALSE); NULL }))
  flags <- c(flags, Filter(Negate(is.null), ext))
  sel_model <- CUTOFF$summary$Selected_Model[CUTOFF$summary$Comparison == comparison_name]
  if (isTRUE(sel_model %in% c("EVS", "Pooled")))
    flags[["Adaptive (pooled vs EVS)"]] <- if (identical(sel_model, "EVS")) flags[["PC1 EVS"]] else flags[["Pooled hierarchical"]]
  MODEL_FLAGS[[comparison_name]] <<- flags
  model_comparison_rows[[comparison_name]] <<- dplyr::bind_rows(lapply(names(flags), function(nm) {
    cnt <- count_discoveries(flags[[nm]])
    j <- match(nm, cm$Model)
    k5 <- grepl("5,000", nm)
    le5 <- if (nm == "Published analysis (k = 5,000)") attr(flags[[nm]], "leading_edge") %||% NA_integer_ else
      if (k5) sum(md$pc1[ids] <= PUBLISHED_TOP_N) else NA_integer_
    data.frame(Comparison = comparison_name, Model = nm,
               k = if (k5) PUBLISHED_TOP_N else if (is.na(j)) NA_integer_ else cm$k[j],
               Leading_Edge = if (k5) as.integer(le5) else if (is.na(j)) NA_integer_ else cm$Leading_Edge[j],
               PAS_tested = cnt[["n_tested"]], Standard = cnt[["Standard"]], Strong = cnt[["Strong"]],
               Weak = cnt[["Weak"]], HBFSS = cnt[["HBFSS"]],
               Selected_by_adaptive_procedure = (nm == "PC1 EVS" && identical(sel_model, "EVS")) ||
                 (nm == "Pooled hierarchical" && identical(sel_model, "Pooled")),
               BIC_descriptive = if (is.na(j)) NA_real_ else cm$BIC[j],
               Delta_BIC_vs_PC1_EVS = if (is.na(j)) NA_real_ else cm$Delta_BIC_vs_PC1_EVS[j])
  }))
  invisible(TRUE)
}

# Comparison-specific and experiment-wide multiplicity (item 8): Standard calls
# with BH within each comparison (primary) and with one BH across the PASs of
# all four comparisons, for every model.
save_multiplicity_table <- function() {
  if (!length(MODEL_FLAGS)) return(invisible(NULL))
  mods <- unique(unlist(lapply(MODEL_FLAGS, names)))
  rows <- lapply(mods, function(m) {
    d <- dplyr::bind_rows(lapply(names(MODEL_FLAGS), function(cmp) {
      f <- MODEL_FLAGS[[cmp]][[m]]
      if (is.null(f)) NULL else data.frame(Comparison = cmp, pvalue = f$pvalue, Standard = f$Standard,
                                           lfc_gate = if (m == "Published analysis (k = 5,000)") is.finite(f$lfc_shrunk) & abs(f$lfc_shrunk) >= lfc_boundary else TRUE)
    }))
    if (!nrow(d)) return(NULL)
    q <- rep(NA_real_, nrow(d)); ok <- is.finite(d$pvalue); q[ok] <- stats::p.adjust(d$pvalue[ok], "BH")
    d$Experiment <- !is.na(q) & q < BH_FDR_STANDARD & d$lfc_gate
    dplyr::bind_rows(lapply(split(d, d$Comparison), function(z)
      data.frame(Comparison = z$Comparison[1], Model = m, Standard_BH_within_comparison = sum(z$Standard),
                 Standard_BH_experiment_wide = sum(z$Experiment))))
  })
  out <- dplyr::bind_rows(rows)
  save_csv(out, file.path(summary_table_dir, "Table_Multiplicity.csv"))
  invisible(out)
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
         subtitle = "b in a0 + a1 x mean^-b (b and a0 trade off; see panel a).") +
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
         subtitle = "Estimated tau^2. Black: LE + REM pooled. Dashed: DESeq2 0.25 floor.") +
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
# EVS CUTOFF k*: IN-SAMPLE HIERARCHICAL SPLIT LIKELIHOOD
# =============================================================================
# Statistical purpose. EVS partitions the PASs of one comparison into a Leading
# Edge (PASs among the top-k |PC1| loadings of either condition: the union) and
# a Remainder, and moderates the dispersions of each partition with its own
# mean-dispersion trend and prior variance. k* is the number of top-|PC1| PASs
# per condition that maximizes the combined hierarchical marginal likelihood of
# the two partitions.
# Data used. Only the ten samples of the comparison: the VST used for ranking,
# the per-PAS Cox-Reid adjusted dispersion likelihoods (design ~ condition) and
# the cutoff search. Samples from other time points are not used, because their
# dispersion structure need not match the tested time point.
# Assumptions. Within a partition, gene-wise log dispersions are normal around
# a power trend a0 + a1 mean^-b (the parametric form of DESeq2; Love et al.
# 2014) with a partition-specific prior variance tau^2; each PAS's adjusted
# profile likelihood (Cox & Reid 1987; McCarthy et al. 2012) is integrated
# over that prior on a grid, giving a marginal likelihood (cf. Robinson & Smyth
# 2007). No p-value, fold change or discovery count enters the choice of k*.
# Because the ranking, the likelihood and the tests use the same counts, the
# likelihood-ratio statistic and BIC are descriptive model-fit measures (the
# adaptive cutoff invalidates their reference distributions); error rates are
# assessed by full-pipeline simulation and the external null splits (VALIDATION).
EXPERIMENT_ARMS <- c("R0", "ZT6", "R2", "ZT8", "R4", "ZT10", "R8", "ZT14")
EXPERIMENT_CACHE <- new.env(parent = emptyenv())
SPLIT_MIN_SET <- 200L
SPLIT_GRID_STEPS <- c(250L, 25L, 5L, 1L)
SPLIT_CORES <- as.integer(Sys.getenv("SEQUENCE_CORES",
  as.character(max(1L, min(4L, parallel::detectCores(logical = FALSE) %||% 1L)))))
if (is.na(SPLIT_CORES) || SPLIT_CORES < 1L) SPLIT_CORES <- 1L
SPLIT_SEED <- 20261125L
CUTOFF <- NULL
# Model comparison (item 4): the alternative splits and the continuous model are
# searched/fitted with the same likelihood and grid resolution as PC1 EVS.
# SEQUENCE_CUTOFF_COMPARATORS=0 skips them (k* and the primary analysis are
# unchanged). Random partitions: SEQUENCE_RANDOM_PARTITIONS (default 3).
CUTOFF_COMPARATORS <- !identical(Sys.getenv("SEQUENCE_CUTOFF_COMPARATORS", "1"), "0")
MODEL_RANDOM_PARTITIONS <- as.integer(Sys.getenv("SEQUENCE_RANDOM_PARTITIONS", "3"))
if (is.na(MODEL_RANDOM_PARTITIONS) || MODEL_RANDOM_PARTITIONS < 0L) MODEL_RANDOM_PARTITIONS <- 3L
MODEL_K <- list()        # per comparison: k of each split model
MODEL_DEPTHS <- list()   # per comparison: union depths of each ranking (incl. random)

# Experiment-wide raw counts (every PAS with a nonzero count in any sample).
# Used only by the separate 40-sample time-course model, never by EVS.
experiment_counts <- function() {
  if (!is.null(EXPERIMENT_CACHE$raw)) return(EXPERIMENT_CACHE$raw)
  x <- as.matrix(WTTS_Seq[, meta_all$id, drop = FALSE])
  storage.mode(x) <- "numeric"
  rownames(x) <- rownames(WTTS_Seq)
  x[!is.finite(x)] <- 0
  x <- x[eligible_pas(x), , drop = FALSE]     # same prespecified abundance filter (any 5 of 40 samples)
  EXPERIMENT_CACHE$raw <- coerce_raw_count_matrix_for_deseq2(x, context = "experiment-wide count matrix")
  EXPERIMENT_CACHE$raw
}

# Arm of each sample (R0 ... ZT14), from the sample-name prefix.
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

# Experiment-wide DESeq2 object with median-of-ratios size factors (time-course
# model only).
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

# ---- Rankings ---------------------------------------------------------------
# One condition's PAS ranks by |PC1 loading| of the VST (rank 1 = smallest,
# rank N = largest). Ties are broken exactly as in the split
# (order(|loading|, decreasing = TRUE) keeps tied PASs in row order), so the top
# k here is the split's top k.
arm_pc1_ranks <- function(rank_matrix, arm_samples) {
  z <- rank_matrix[, arm_samples, drop = FALSE]
  loading <- stats::prcomp(t(z), center = TRUE, scale. = FALSE, rank. = 1)$rotation[, 1]
  loading[!is.finite(loading)] <- 0
  ord <- rev(order(abs(loading), decreasing = TRUE))
  stats::setNames(match(seq_along(loading), ord), rownames(rank_matrix))
}

# Position from the top (1 = largest score) of every PAS for a per-arm score.
score_positions <- function(score) stats::setNames(rank(-score, ties.method = "first"), names(score))

# Alternative rankings of one condition (model comparison):
#   pc12      lambda1 v1^2 + lambda2 v2^2, the PAS's weighted contribution to
#             the PC1-PC2 subspace (summing over all components would give the
#             variance);
#   variance  within-condition VST variance.
arm_comparator_positions <- function(rank_matrix, arm_samples) {
  z <- rank_matrix[, arm_samples, drop = FALSE]
  pc <- stats::prcomp(t(z), center = TRUE, scale. = FALSE, rank. = 2)
  lam <- pc$sdev[1:2]^2
  s12 <- lam[1] * pc$rotation[, 1]^2 + lam[2] * pc$rotation[, 2]^2
  list(pc12 = score_positions(stats::setNames(s12, rownames(z))),
       variance = score_positions(stats::setNames(apply(z, 1L, stats::var), rownames(z))))
}

# Comparison-specific ranking (item 1). The VST (DESeq2, blind = TRUE) is
# computed on the ten samples of the comparison alone, with the comparison's
# shared size factors. Within each condition PASs are positioned by |PC1
# loading| (and, for the model comparison, by PC1-PC2 contribution and by
# variance). The union depth of a PAS is the smaller of its two positions: the
# PAS is in the Leading Edge at k exactly when its union depth is <= k.
comparison_rankings <- function(raw, coldata, sf) {
  dds <- DESeq2::DESeqDataSetFromMatrix(raw, colData = coldata, design = ~ 1)
  dds <- apply_size_factors(dds, sf)
  v <- as.matrix(SummarizedExperiment::assay(suppressMessages(
    DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))))
  dimnames(v) <- dimnames(raw)
  if (any(!is.finite(v))) stop("The comparison VST matrix has non-finite values.")
  cond <- as.character(coldata[colnames(raw), "condition"])
  N <- nrow(v)
  arm_pos <- function(samples) {
    cp <- arm_comparator_positions(v, samples)
    list(pc1 = N - arm_pc1_ranks(v, samples) + 1L, pc12 = cp$pc12, variance = cp$variance)
  }
  tr <- arm_pos(colnames(v)[cond == "trt"]); un <- arm_pos(colnames(v)[cond == "untrt"])
  depths <- lapply(c(pc1 = "pc1", pc12 = "pc12", variance = "variance"), function(m)
    stats::setNames(as.integer(pmin(tr[[m]][rownames(v)], un[[m]][rownames(v)])), rownames(v)))
  list(vst = v, depths = depths, rank_RT = tr$pc1[rownames(v)], rank_ZT = un$pc1[rownames(v)])
}

# Everything the analysis of one comparison reuses: raw counts, shared size
# factors (median of ratios on the comparison matrix), gene-wise dispersions
# (design ~ condition), the comparison VST and the union depths.
build_comparison_cache <- function(raw, coldata) {
  sf <- comparison_size_factors(raw, coldata)
  gdds <- comparison_gene_dds(raw, coldata, sf)
  rk <- comparison_rankings(raw, coldata, sf)
  list(raw = raw, coldata = coldata, sf = sf, gene_dds = gdds, vst = rk$vst,
       depth = rk$depths$pc1, depths = rk$depths, rank_RT = rk$rank_RT, rank_ZT = rk$rank_ZT)
}

# Store the cache of a comparison once (later calls with the same matrix reuse it).
prepare_comparison_cache <- function(comparison_name, count_matrix, coldata) {
  raw <- coerce_raw_count_matrix_for_deseq2(count_matrix, context = paste0(comparison_name, " comparison matrix"))
  hit <- comparison_cache[[comparison_name]]
  if (!is.null(hit) && identical(dimnames(hit$raw), dimnames(raw))) return(invisible(TRUE))
  comparison_cache[[comparison_name]] <<- build_comparison_cache(raw, coldata)
  invisible(TRUE)
}

# ---- Split likelihood and grid-refined search -------------------------------
# Inputs of the split likelihood: informative PASs, their log means and their
# relative adjusted likelihood curves over the dispersion grid.
split_inputs <- function(cache) {
  lk <- make_dispersion_likelihood(cache$gene_dds)
  mc <- S4Vectors::mcols(cache$gene_dds)
  ids <- rownames(lk)
  list(ids = ids, lx = log(mc$baseMean[match(ids, rownames(cache$gene_dds))]), lk = lk,
       mean = stats::setNames(mc$baseMean, rownames(cache$gene_dds)),
       disp = stats::setNames(mc$dispGeneEst, rownames(cache$gene_dds)))
}

# Combined marginal log-likelihood of the two partitions at cutoff k. Each side
# is fitted with the same likelihood and bounds as the final moderation; fits
# are warm-started from the nearest evaluated cutoff and fall back to the
# pooled start (with multiple starts) if the warm start fails.
split_loglik <- function(k, depth, lx, likelihood, start_le, start_rem, fallback) {
  le <- depth <= k
  a <- fit_marginal_trend(lx[le], likelihood[le,,drop=FALSE], start_le, fallback = fallback)
  b <- fit_marginal_trend(lx[!le], likelihood[!le,,drop=FALSE], start_rem, fallback = fallback)
  list(k=k, loglik=a$loglik+b$loglik, le=a, rem=b, converged=a$converged && b$converged)
}

# Grid-refined maximum of the split likelihood over k (item 1: the efficient
# search is preserved). The first stage evaluates the admissible endpoints and a
# grid of step steps[1]; each later stage refines around the three best
# candidates with the next step (250, 25, 5, 1). Ties prefer smaller k. Each
# partition must hold at least SPLIT_MIN_SET PASs. This is a grid-refined
# maximum, not a proven global maximum.
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
  cv <- vapply(prof[as.character(coarse_k[as.character(coarse_k) %in% names(prof)])], `[[`, numeric(1), "loglik")
  coarse_max <- if (any(is.finite(cv))) max(cv[is.finite(cv)]) else NA_real_
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
  list(k=k, best=prof[[as.character(k)]], prof=prof, coarse_max_loglik=coarse_max)
}

# k* of one ranking for a prepared cache (used by the primary analysis and by
# every simulation, subsample and external analysis, so all use the same rule).
select_cutoff <- function(si, depth_all, one) {
  k_max <- floor(length(depth_all) / 2)
  search_cutoff(depth_all[si$ids], si$lx, si$lk, one, k_max, depth_all)
}

# Continuous-rank model: a natural cubic spline (3 df) of log(1 + PC1 union
# depth), centred, shifts the prior centre of each PAS; one prior variance.
continuous_basis <- function(depth, ids) {
  B <- splines::ns(log1p(as.numeric(depth[ids])), df = 3)
  B <- scale(unclass(B)[, , drop = FALSE], center = TRUE, scale = FALSE)
  attr(B, "scaled:center") <- NULL
  rownames(B) <- ids
  B
}

# ---- Pooled versus split selection ---------------------------------------------
# Purpose. Allow the no-split outcome: EVS is used for the primary results only
# if splitting improves the fit by more than the complete adaptive procedure
# (ranking plus k search) achieves on data without two dispersion populations.
# Rationale. Under one population k is not identified, so the likelihood-ratio
# statistic has no chi-square reference and BIC ignores the search; its null
# distribution is therefore obtained by parametric bootstrap from the fitted
# pooled model (the reference distribution of the criterion under homogeneous dispersion).
# Each replicate simulates counts with the comparison's fitted means (DESeq2
# mu, so true condition effects are kept), size factors and design, and
# dispersions drawn from the fitted pooled trend and tau^2; it is filtered,
# re-ranked (VST, per-condition PC1), its likelihood recomputed and its split
# searched on the coarse grid. The statistic is the coarse-grid LRT, compared
# with the observed coarse-grid LRT at the same resolution, and
# p = (1 + #{LRT* >= LRT_obs}) / (B + 1) (Phipson & Smyth 2010). EVS is
# selected if p <= SELECTION_ALPHA; otherwise the pooled model is selected and
# EVS is reported as a forced analysis. The selection rule's behaviour under
# heterogeneous dispersion is reported in the simulations.
SELECTION_BOOT <- env_int("SEQUENCE_SELECTION_BOOT", if (PUB) 99L else 19L)
if (is.na(SELECTION_BOOT) || SELECTION_BOOT < 0L) SELECTION_BOOT <- 19L
SELECTION_ALPHA <- 0.05
SELECTION <- list()      # per comparison: bootstrap statistics, observed statistic, p, critical value
# A decision requires at least 80% of the bootstrap replicates to succeed;
# otherwise it is "Undetermined" and no adaptive result is reported.
SELECTION_MIN_FRACTION <- 0.8

# The selection rule, shared by cutoff_one() and reanalyze_matrix(): finite
# bootstrap statistics only, p = (1 + #{LRT* >= LRT}) / (B_ok + 1), EVS when
# p <= SELECTION_ALPHA.
selection_decision <- function(observed, boot) {
  ok <- boot[is.finite(boot)]
  if (!length(boot) || length(ok) < ceiling(SELECTION_MIN_FRACTION * length(boot)) || !is.finite(observed))
    return(list(p = NA_real_, crit = NA_real_, model = "Undetermined", n_ok = length(ok), n = length(boot)))
  p <- (1 + sum(ok >= observed)) / (length(ok) + 1)
  list(p = p, crit = unname(stats::quantile(ok, 1 - SELECTION_ALPHA, type = 1)),
       model = if (p <= SELECTION_ALPHA) "EVS" else "Pooled", n_ok = length(ok), n = length(boot))
}

coarse_lrt <- function(si, depth_all, one) {
  k_max <- floor(length(depth_all) / 2)
  s <- search_cutoff(depth_all[si$ids], si$lx, si$lk, one, k_max, depth_all, steps = SPLIT_GRID_STEPS[1])
  2 * (s$coarse_max_loglik - one$loglik)
}

selection_bootstrap <- function(cache, one, seed, filter = filter_rule(cache$coldata)) {
  set.seed(seed)
  gd <- cache$gene_dds
  mu <- as.matrix(SummarizedExperiment::assays(gd)[["mu"]]); dimnames(mu) <- dimnames(gd)
  bm <- S4Vectors::mcols(gd)$baseMean
  ok <- rowSums(is.finite(mu) & mu > 0) == ncol(mu) & is.finite(bm) & bm > 0
  mu <- mu[ok, , drop = FALSE]; bm <- bm[ok]
  p <- unname(one$par); f <- make_trend_function(exp(p[1]), exp(p[2]), p[3])
  vapply(seq_len(SELECTION_BOOT), function(b) {
    alpha <- pmin(pmax(exp(log(f(bm)) + stats::rnorm(length(bm), 0, sqrt(exp(p[4])))), 1e-8), 1e3)
    y <- matrix(stats::rnbinom(length(mu), mu = as.vector(mu), size = rep(1 / alpha, times = ncol(mu))),
                nrow(mu), dimnames = dimnames(mu))
    y <- coerce_raw_count_matrix_for_deseq2(y[eligible_pas(y, filter[1], filter[2]), , drop = FALSE], context = "selection bootstrap")
    tryCatch({
      cb <- build_comparison_cache(y, cache$coldata)
      sb <- split_inputs(cb)
      ob <- fit_marginal_trend(sb$lx, sb$lk, multistart = TRUE)
      if (!ob$converged) NA_real_ else coarse_lrt(sb, cb$depths$pc1, ob)
    }, error = function(e) NA_real_)
  }, numeric(1))
}

# Numerical checks at k* (item 2): optimizer restarts, a finer dispersion grid
# (step / 2) and a lower tau^2 bound (bound / 10). The globals changed here are
# restored before returning.
cutoff_numerical_checks <- function(cache, si, depth, k, best, depth_full) {
  le_k <- depth <= k
  ms <- function(lkm) list(le = fit_marginal_trend(si$lx[le_k], lkm[le_k,,drop=FALSE], multistart = TRUE),
                           rem = fit_marginal_trend(si$lx[!le_k], lkm[!le_k,,drop=FALSE], multistart = TRUE))
  old <- list(g = DISP_LOG_GRID, e = DISP_LOG_EDGES, m = DISP_PRIOR_VAR_MIN)
  restore <- function() { DISP_LOG_GRID <<- old$g; DISP_LOG_EDGES <<- old$e; DISP_PRIOR_VAR_MIN <<- old$m }
  on.exit(restore(), add = TRUE)
  r0 <- ms(si$lk)
  DISP_LOG_GRID <<- seq(min(old$g), max(old$g), by = DISP_GRID_STEP / 2)
  DISP_LOG_EDGES <<- c(-Inf, head(DISP_LOG_GRID, -1) + diff(DISP_LOG_GRID) / 2, Inf)
  lk2 <- make_dispersion_likelihood(cache$gene_dds)[si$ids, , drop = FALSE]
  r1 <- ms(lk2)
  restore(); DISP_PRIOR_VAR_MIN <<- old$m / 10
  r2 <- ms(si$lk)
  restore()
  pd <- function(a, b) max(abs(c(a$le$par[1:3] - b$le$par[1:3], a$rem$par[1:3] - b$rem$par[1:3])))
  tau <- function(a) c(a$le$tau2, a$rem$tau2)
  # Analytic gradient vs central differences at the fitted LE and REM parameters.
  gerr <- max(gradient_check(best$le$par, si$lx[le_k], si$lk[le_k, , drop = FALSE]),
              gradient_check(best$rem$par, si$lx[!le_k], si$lk[!le_k, , drop = FALSE]))
  # The other likelihood (fixed-mean if the primary is profiled, and vice
  # versa): change in gene-wise dispersion estimates and in k*.
  lk_alt <- make_dispersion_likelihood(cache$gene_dds, profile = !DISP_PROFILE)
  common <- intersect(si$ids, rownames(lk_alt))
  dla <- abs(DISP_LOG_GRID[max.col(si$lk[common, , drop = FALSE], ties.method = "first")] -
             DISP_LOG_GRID[max.col(lk_alt[common, , drop = FALSE], ties.method = "first")])
  si_alt <- list(ids = common, lx = si$lx[match(common, si$ids)], lk = lk_alt[common, , drop = FALSE])
  one_alt <- fit_marginal_trend(si_alt$lx, si_alt$lk, multistart = TRUE)
  k_alt <- tryCatch(select_cutoff(si_alt, depth_full, one_alt)$k, error = function(e) NA_integer_)
  le_full <- names(depth_full)[depth_full <= k]
  le_alt <- if (is.na(k_alt)) character() else names(depth_full)[depth_full <= k_alt]
  data.frame(Check_gradient_max_rel_error = gerr,
             Check_other_likelihood = if (DISP_PROFILE) "fixed-mean" else "profiled",
             Check_other_likelihood_median_abs_dlogalpha = stats::median(dla),
             Check_other_likelihood_k = k_alt,
             Check_other_likelihood_LE_Jaccard = if (length(le_alt)) length(intersect(le_full, le_alt)) / length(union(le_full, le_alt)) else NA_real_,
             Check_multistart_dLogLik = (r0$le$loglik + r0$rem$loglik) - best$loglik,
             Check_finer_grid_max_par_change = pd(r1, r0),
             Check_finer_grid_max_tau2_rel_change = max(abs(tau(r1) / tau(r0) - 1)),
             Check_lower_tau2_bound_max_par_change = pd(r2, r0),
             Check_lower_tau2_bound_tau2_LE = r2$le$tau2, Check_lower_tau2_bound_tau2_REM = r2$rem$tau2)
}

# k*, the model comparison and the numerical checks of one comparison.
cutoff_one <- function(i) {
  set.seed(SPLIT_SEED + i)
  cmp <- comparison_table$comparison_name[i]
  cache <- comparison_cache[[cmp]]
  if (is.null(cache)) stop(cmp, ": comparison inputs were not prepared")
  si <- split_inputs(cache)
  ids <- rownames(cache$raw); n <- length(si$ids)
  if (n < 2L * SPLIT_MIN_SET) stop(cmp, ": too few informative PASs")
  one <- fit_marginal_trend(si$lx, si$lk, multistart = TRUE)
  if (!one$converged) stop(cmp, ": pooled model failed to converge")
  l_one <- one$loglik
  depth_all <- cache$depths$pc1
  sr <- select_cutoff(si, depth_all, one); k <- sr$k; best <- sr$best; prof <- sr$prof
  depth <- depth_all[si$ids]
  lrt <- 2 * (best$loglik - l_one); n_le <- sum(depth <= k)
  bic <- function(ll, p) -2 * ll + p * log(n)
  # Pooled versus split selection (parametric bootstrap of the adaptive search).
  sel_obs <- 2 * (sr$coarse_max_loglik - l_one)
  sel_boot <- if (SELECTION_BOOT > 0L) selection_bootstrap(cache, one, SPLIT_SEED + 7000L * i) else numeric()
  sel_dec <- selection_decision(sel_obs, sel_boot)
  sel_ok <- sel_boot[is.finite(sel_boot)]
  sel_p <- sel_dec$p; sel_crit <- sel_dec$crit; sel_model <- sel_dec$model
  models <- list(data.frame(Comparison = cmp, Model = "Pooled hierarchical", k = NA_integer_, Leading_Edge = NA_integer_,
                            LogLik = l_one, Parameters = 4L, BIC = bic(l_one, 4)),
                 data.frame(Comparison = cmp, Model = "PC1 EVS", k = k, Leading_Edge = sum(depth_all <= k),
                            LogLik = best$loglik, Parameters = 9L, BIC = bic(best$loglik, 9)))
  model_k <- list(pc1 = k)
  model_depths <- list(pc1 = depth_all)
  if (CUTOFF_COMPARATORS) {
    for (m in c("variance", "pc12")) {
      s <- tryCatch(select_cutoff(si, cache$depths[[m]], one), error = function(e) NULL)
      if (is.null(s)) next
      model_k[[m]] <- s$k; model_depths[[m]] <- cache$depths[[m]]
      models[[length(models) + 1L]] <- data.frame(Comparison = cmp,
        Model = if (m == "variance") "Variance split" else "PC1-PC2 split", k = s$k,
        Leading_Edge = sum(cache$depths[[m]] <= s$k), LogLik = s$best$loglik, Parameters = 9L, BIC = bic(s$best$loglik, 9))
    }
    B <- continuous_basis(depth_all, ids)
    cf <- fit_marginal_trend(si$lx, si$lk, one$par, basis = B[si$ids, , drop = FALSE], fallback = one$par)
    if (cf$converged) models[[length(models) + 1L]] <- data.frame(Comparison = cmp, Model = "Continuous rank",
      k = NA_integer_, Leading_Edge = NA_integer_, LogLik = cf$loglik, Parameters = 7L, BIC = bic(cf$loglik, 7))
    for (r in seq_len(MODEL_RANDOM_PARTITIONS)) {
      dp <- stats::setNames(sample(unname(depth_all)), ids)
      s <- tryCatch(select_cutoff(si, dp, one), error = function(e) NULL)
      if (is.null(s)) next
      nm <- paste0("random", r)
      model_k[[nm]] <- s$k; model_depths[[nm]] <- dp
      models[[length(models) + 1L]] <- data.frame(Comparison = cmp, Model = paste("Random partition", r), k = s$k,
        Leading_Edge = sum(dp <= s$k), LogLik = s$best$loglik, Parameters = 9L, BIC = bic(s$best$loglik, 9))
    }
  }
  models <- dplyr::bind_rows(models)
  models$Delta_BIC_vs_PC1_EVS <- models$BIC - models$BIC[models$Model == "PC1 EVS"]
  pr <- data.frame(Comparison=cmp,k=as.integer(names(prof)),LogLik=vapply(prof,`[[`,numeric(1),"loglik"),
    Leading_Edge_fitted=vapply(as.integer(names(prof)),function(kk) sum(depth<=kk),integer(1)),
    Converged=vapply(prof,function(p) isTRUE(p$converged),logical(1)))
  pr <- pr[order(pr$k),,drop=FALSE]; pr$LogLik[!pr$Converged] <- NA_real_
  pr$LRT_vs_one_trend <- 2*(pr$LogLik-l_one); pr$Delta_from_max <- 2*(pr$LogLik-best$loglik); pr$Chosen <- pr$k==k
  drop_at <- function(d) {j<-which.min(abs(pr$k-(k+d))); if(abs(pr$k[j]-(k+d))<=125) pr$Delta_from_max[j] else NA_real_}
  q <- function(p) { p <- unname(p); c(a0=exp(p[1]),a1=exp(p[2]),b=p[3]) }
  pl <- q(best$le$par); pm <- q(best$rem$par); po <- q(one$par)
  bm <- models[is.finite(models$BIC), , drop = FALSE]
  summary <- data.frame(Comparison=cmp,EVS_PCA_Method=EVS_PCA_METHOD,k=k,
    Leading_Edge=sum(depth_all<=k),Remainder=sum(depth_all>k),Samples=paste(colnames(cache$raw),collapse=","),
    n_PAS_fitted=n,LE_fitted=n_le,REM_fitted=n-n_le,LogLik_split=best$loglik,LogLik_one_trend=l_one,LRT_vs_one_trend=lrt,
    BIC_pooled=bic(l_one,4),BIC_EVS=bic(best$loglik,9),
    Lowest_BIC_Model=bm$Model[which.min(bm$BIC)],
    Coarse_LRT_observed=sel_obs, Selection_bootstrap_replicates=length(sel_ok),
    Selection_bootstrap_max=if (length(sel_ok)) max(sel_ok) else NA_real_,
    Selection_critical_value=sel_crit, Selection_p=sel_p,
    Selection_rule=sprintf("parametric bootstrap, alpha %.2f, >= %.0f%% successful replicates (%d of %d)",
                           SELECTION_ALPHA, 100 * SELECTION_MIN_FRACTION, sel_dec$n_ok, sel_dec$n),
    Selected_Model=sel_model,
    Delta_at_minus500=drop_at(-500),Delta_at_plus500=drop_at(500),
    LE_a0=pl[["a0"]],LE_a1=pl[["a1"]],LE_b=pl[["b"]],LE_resid_sd=sqrt(best$le$tau2)/log(10),
    REM_a0=pm[["a0"]],REM_a1=pm[["a1"]],REM_b=pm[["b"]],REM_resid_sd=sqrt(best$rem$tau2)/log(10),
    One_a0=po[["a0"]],One_a1=po[["a1"]],One_b=po[["b"]],
    LE_prior_var=best$le$tau2,REM_prior_var=best$rem$tau2,One_prior_var=one$tau2,
    LE_prior_var_at_bound=best$le$tau2<=DISP_PRIOR_VAR_MIN*1.00001, REM_prior_var_at_bound=best$rem$tau2<=DISP_PRIOR_VAR_MIN*1.00001,
    Converged=isTRUE(best$converged) && isTRUE(one$converged))
  num <- tryCatch(cutoff_numerical_checks(cache, si, depth, k, best, depth_all),
                  error = function(e) data.frame(Check_gradient_max_rel_error = NA_real_, Check_multistart_dLogLik = NA_real_, Check_finer_grid_max_par_change = NA_real_,
                                                 Check_finer_grid_max_tau2_rel_change = NA_real_, Check_lower_tau2_bound_max_par_change = NA_real_,
                                                 Check_lower_tau2_bound_tau2_LE = NA_real_, Check_lower_tau2_bound_tau2_REM = NA_real_))
  summary <- cbind(summary, num)
  pas <- data.frame(Comparison=cmp,feature_id=ids,Union_Depth=unname(depth_all[ids]),
    Mean=unname(si$mean[ids]),Gene_Dispersion=unname(si$disp[ids]),In_Trend_Fit=ids %in% si$ids,
    In_Leading_Edge=unname(depth_all[ids]<=k))
  message(sprintf("%s: k* = %d (LE %d); LRT %.0f; selection p = %s -> %s; lowest BIC: %s", cmp, k, summary$Leading_Edge, lrt,
                  format(signif(sel_p, 3)), sel_model, summary$Lowest_BIC_Model))
  list(cmp=cmp, k=k, summary=summary, profile=pr, models=models, pas=pas, model_k=model_k, model_depths=model_depths,
       fits=list(one=one$par, le=best$le$par, rem=best$rem$par),
       selection=list(boot=sel_boot, observed=sel_obs, p=sel_p, crit=sel_crit, model=sel_model))
}

# Prepare every comparison (sequentially) and lock k* and the model comparison
# for all of them before any DE test. Comparisons are independent and run in
# parallel (SEQUENCE_CORES); each uses its own fixed seed, so results do not
# depend on the number of cores.
compute_cutoffs <- function() {
  for (i in seq_len(nrow(comparison_table))) {
    inp <- prepare_comparison_data(comparison_table$comparison_name[i], comparison_table$group1_prefix[i],
                                   comparison_table$group2_prefix[i], WTTS_Seq, meta_all)
    prepare_comparison_cache(inp$comparison_name, inp$count_matrix, inp$coldata)
  }
  cores <- if (.Platform$OS.type == "unix") SPLIT_CORES else 1L
  res <- parallel::mclapply(seq_len(nrow(comparison_table)), function(i)
    tryCatch(cutoff_one(i), error = function(e) structure(conditionMessage(e), class = "cutoff_error")),
    mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) inherits(r, "cutoff_error") || inherits(r, "try-error") || is.null(r), logical(1))
  if (any(bad)) stop("Cutoff selection failed: ", paste(vapply(res[bad], function(r) paste(as.character(r), collapse=" "), ""), collapse = " | "))
  nms <- vapply(res, `[[`, "", "cmp")
  for (r in res) { OPTIMAL_K[[r$cmp]] <<- r$k; MODEL_K[[r$cmp]] <<- r$model_k; MODEL_DEPTHS[[r$cmp]] <<- r$model_depths
                   SELECTION[[r$cmp]] <<- r$selection }
  CUTOFF <<- list(summary=dplyr::bind_rows(lapply(res,`[[`,"summary")),
                  profiles=stats::setNames(lapply(res,`[[`,"profile"), nms),
                  pas=stats::setNames(lapply(res,`[[`,"pas"), nms),
                  models=dplyr::bind_rows(lapply(res,`[[`,"models")),
                  fits=stats::setNames(lapply(res,`[[`,"fits"), nms))
  invisible(CUTOFF)
}

# ---- cutoff outputs: tables and Figure 1 --------------------------------------

# Purpose: Save cutoff outputs.
save_cutoff_outputs <- function() {
  if (is.null(CUTOFF)) return(invisible(FALSE))
  out <- summary_table_dir
  save_csv(CUTOFF$summary, file.path(out, "Table_Cutoff_Summary.csv"))
  save_csv(dplyr::bind_rows(CUTOFF$profiles), file.path(out, "Table_Cutoff_Profile_Likelihood.csv"))
  save_csv(CUTOFF$models, file.path(out, "Table_Model_Comparison_BIC.csv"))
  save_csv(dplyr::bind_rows(CUTOFF$pas), file.path(out, "Table_Cutoff_PAS_Dispersions.csv"))
  le <- dplyr::bind_rows(lapply(CUTOFF$pas, function(p) p[p$In_Leading_Edge, c("Comparison", "feature_id")]))
  save_csv(le, file.path(out, "Table_Cutoff_Leading_Edge_Sites.csv"))
  try_fig(save_figure(plot_cutoff_figure(), paper_fig_dir, "Figure_1_EVS_Cutoff", 180, 140), "Figure 1 (cutoff)")
  invisible(TRUE)
}

# Figure 1: (a) gene-wise dispersions of the comparison with the Leading-Edge and Remainder trends
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
    xr <- range(pts$Mean[pts$Comparison == summ$Comparison[i]])
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
  pa <- ggplot(pts[is.finite(pts$Gene_Dispersion) & pts$Gene_Dispersion > 1e-6, , drop = FALSE], aes(Mean, Gene_Dispersion)) +
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
    labs(x = "Mean of normalized counts", y = "Gene-wise dispersion", tag = "a",
         title = "Dispersion trends of the Leading Edge and the Remainder at k* (the ten samples of each comparison)",
         subtitle = "Points: PASs coloured by set at k* (up to 12,000 per panel). Lines: a0 + a1 x mean^-b fitted to each set; dashed: one trend for all PASs.") +
    fig_theme()
  prof <- dplyr::bind_rows(CUTOFF$profiles)
  k_of <- stats::setNames(CUTOFF$summary$k, CUTOFF$summary$Comparison)
  # Retain all evaluated grid points, including shifted coarse-grid points.
  prof$Comparison <- as_comp(prof$Comparison)
  kl <- data.frame(Comparison = summ$Comparison, k = summ$k, y = summ$LRT_vs_one_trend,
                   lab = sprintf("k* = %s\nLRT = %s", format(summ$k, big.mark = ","),
                                 format(round(summ$LRT_vs_one_trend), big.mark = ",")))
  pb <- ggplot(prof, aes(k, LRT_vs_one_trend)) +
    geom_line(linewidth = 0.45, colour = "grey15", na.rm=TRUE) +
    geom_vline(data = kl, aes(xintercept = k), colour = "#D55E00", linewidth = 0.4, linetype = "22") +
    geom_point(data = kl, aes(k, y), colour = "#D55E00", size = 1.4) +
    geom_text(data = kl, aes(k, y, label = lab), hjust = -0.08, vjust = 1.2, size = 1.9, colour = "#D55E00",
              fontface = "bold", lineheight = 0.9) +
    facet_wrap(~ Comparison, nrow = 1, scales = "free_y", labeller = labeller(Comparison = comp_lab())) +
    scale_x_continuous(labels = scales::label_comma()) + scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "k (top-|PC1| PASs per arm)", y = "LRT, split vs one trend", tag = "b",
         title = "Profile adjusted likelihood: grid-refined candidate cutoff",
         subtitle = "LRT = 2 x [marginal adjusted log-likelihood difference], split vs one trend (descriptive: k* is chosen on the same counts).") +
    fig_theme()
  assemble_rows(list(pa, pb), heights = c(1.1, 1),
                title = "Figure 1. EVS: per-condition PC1 ranking and the selected cutoff (grid-refined marginal adjusted likelihood)")
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
# call_significance), HBFSS thresholded per view.
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
  n <- function(x) if (is.null(x) || all(is.na(x))) NA_integer_ else sum(x %in% TRUE)
  c(Standard = n(flags$Standard), Strong = n(flags$Strong), Weak = n(flags$Weak), HBFSS = n(flags$HBFSS),
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
  full <- cc$vst
  add_validation("PCA input is finite, matches the raw PAS IDs and samples, and is the comparison's own VST",
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
  add_validation("PC1 loadings index PASs, not samples", comparison_name,
                 all(vapply(fits, function(f) nrow(f$loading_table) == nrow(raw) && setequal(f$loading_table$feature_id, rownames(raw)), logical(1))),
                 sprintf("%d PASs; %d samples per arm", nrow(raw), ncol(raw) / 2))
  add_validation("Deterministic execution: recomputing the ranking gives identical union depths", comparison_name,
                 identical(comparison_rankings(raw, coldata, cc$sf)$depths$pc1, cc$depth))
  rk <- lapply(fits, function(f) stats::setNames(f$loading_table$rank, f$loading_table$feature_id)[rownames(raw)])
  k <- OPTIMAL_K[[comparison_name]]
  same_top <- setequal(names(rk$trt)[rk$trt <= k], names(cc$rank_RT)[cc$rank_RT <= k]) &&
    setequal(names(rk$untrt)[rk$untrt <= k], names(cc$rank_ZT)[cc$rank_ZT <= k])
  add_validation("PC1 top-k* sets of the split equal those used by the cutoff search", comparison_name, same_top,
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
  if (length(model_comparison_rows)) save_csv(dplyr::bind_rows(model_comparison_rows), file.path(out, "Table_Model_Comparison_DE.csv"))
  if (length(normalization_rows)) save_csv(dplyr::bind_rows(normalization_rows), file.path(out, "Table_Normalization_Check.csv"))
  if (!is.null(CUTOFF)) try_fig(save_figure(plot_model_sensitivity_figure(), paper_fig_dir,
                                            "Supplementary_Figure_S1_Model_Fit_Sensitivity", 180, 170), "Supplementary Figure S1")
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
    LRT_vs_one_trend_descriptive = cs$LRT_vs_one_trend, Selected_Model = cs$Selected_Model,
    Lowest_BIC_Model_descriptive = cs$Lowest_BIC_Model,
    LE_prior_var = cs$LE_prior_var, REM_prior_var = cs$REM_prior_var,
    LE_a0 = cs$LE_a0, LE_a1 = cs$LE_a1, LE_b = cs$LE_b, REM_a0 = cs$REM_a0, REM_a1 = cs$REM_a1, REM_b = cs$REM_b,
    Trend_fits_converged = cs$Converged,
    EVS_PCA_Method = EVS_PCA_METHOD, PCA_Preprocessing = PCA_METHOD_LABELS[[EVS_PCA_METHOD]],
    Normalization_Mode = NORMALIZATION_MODE)
  cc <- comparison_cache[[comparison_name]]
  add_validation("Ranking, dispersion likelihood and cutoff use only the comparison's ten samples", comparison_name,
                 identical(colnames(cc$vst), colnames(cc$raw)) && identical(colnames(cc$gene_dds), colnames(cc$raw)) &&
                   identical(strsplit(cs$Samples, ",")[[1]], colnames(cc$raw)) && ncol(cc$raw) == 10L,
                 sprintf("samples: %s", cs$Samples))
  add_validation("Dispersion trend fits at k* converged", comparison_name, isTRUE(cs$Converged),
                 sprintf("LE %.3f + %.2f x mean^-%.2f; REM %.3f + %.2f x mean^-%.2f",
                         cs$LE_a0, cs$LE_a1, cs$LE_b, cs$REM_a0, cs$REM_a1, cs$REM_b))
  add_validation("Marginal likelihood: analytic gradient matches central differences (relative error < 1e-3)", comparison_name,
                 if (is.null(cs$Check_gradient_max_rel_error) || !is.finite(cs$Check_gradient_max_rel_error)) NA else cs$Check_gradient_max_rel_error < 1e-3,
                 sprintf("largest relative error %.2g", cs$Check_gradient_max_rel_error %||% NA_real_))
  sel <- SELECTION[[comparison_name]]
  sel_ok <- if (is.null(sel)) numeric() else sel$boot[is.finite(sel$boot)]
  rule_model <- if (is.null(sel)) NA_character_ else selection_decision(sel$observed, sel$boot)$model
  add_validation("Pooled-vs-split decision follows the prespecified bootstrap rule (pooled model retained when the split is unsupported)",
                 comparison_name, identical(rule_model, cs$Selected_Model),
                 sprintf("observed coarse LRT %.1f; %d bootstrap replicates (max %.1f); p = %s; selected %s",
                         cs$Coarse_LRT_observed, length(sel_ok), if (length(sel_ok)) max(sel_ok) else NA_real_,
                         format(signif(cs$Selection_p, 3)), cs$Selected_Model))
  add_validation("Cutoff fit: restarts from several starting values do not improve the likelihood (< 0.5)", comparison_name,
                 if (is.null(cs$Check_multistart_dLogLik) || !is.finite(cs$Check_multistart_dLogLik)) NA else cs$Check_multistart_dLogLik < 0.5,
                 sprintf("log-likelihood change %.3g", cs$Check_multistart_dLogLik %||% NA_real_))
  add_validation("Cutoff fit: halving the dispersion grid step changes the trend parameters by < 0.05", comparison_name,
                 if (is.null(cs$Check_finer_grid_max_par_change) || !is.finite(cs$Check_finer_grid_max_par_change)) NA else cs$Check_finer_grid_max_par_change < 0.05,
                 sprintf("largest parameter change %.3g (log a0, log a1, b); tau2 relative change %.3g",
                         cs$Check_finer_grid_max_par_change %||% NA_real_, cs$Check_finer_grid_max_tau2_rel_change %||% NA_real_))
  add_validation("Cutoff fit: lowering the prior-variance bound tenfold leaves the trend parameters unchanged (< 0.05)", comparison_name,
                 if (is.null(cs$Check_lower_tau2_bound_max_par_change) || !is.finite(cs$Check_lower_tau2_bound_max_par_change)) NA else cs$Check_lower_tau2_bound_max_par_change < 0.05,
                 sprintf("largest parameter change %.3g", cs$Check_lower_tau2_bound_max_par_change %||% NA_real_))
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
      panel_title <- paste0("HBFSS threshold | Orig | ", CUTOFF_METHOD_LABEL)
      output_suffix <- "Raw_AllComparisons"
      panel_height <- 5.8
    } else {
      track_order_use <- "normalized_evs"
      panel_title <- paste0("HBFSS threshold | ", unname(dataset_short[dataset_key]), " | ", CUTOFF_METHOD_LABEL)
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
  try_fig(save_figure(plot_time_course_figure(), paper_fig_dir, "Figure_8_Time_Course", 180, 120), "Figure 8 (time course)")
  invisible(TIME_COURSE)
}

# Figure 8: (a) significant PASs per recovery time in the time-course model,
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
                title = "Figure 8. Recovery time x condition model across all 40 samples")
}

# =============================================================================
# FIGURE 3: MODEL COMPARISON AND SENSITIVITY
# =============================================================================
# (a) BIC of each alternative minus BIC of EVS (negative = better than EVS);
# (b) observed coarse-grid LRT vs the bootstrap searches; (c) Standard
# discoveries across the primary and sensitivity analyses.
# Supplementary Figure S1: (a) BIC of each model relative to PC1 EVS (descriptive model fit;
# the adaptive cutoff invalidates formal BIC inference), (b) Standard
# discoveries of the complete DE analysis of every model (descriptive; FDR and
# power are judged from the simulations, Figure 3), (c) sensitivity analyses.
plot_model_sensitivity_figure <- function() {
  mods <- CUTOFF$models
  mods <- mods[mods$Model != "PC1 EVS" & is.finite(mods$Delta_BIC_vs_PC1_EVS), , drop = FALSE]
  mods$Model <- ifelse(grepl("^Random partition", mods$Model), "Random partition", mods$Model)
  da <- mods |> dplyr::group_by(Comparison, Model) |>
    dplyr::summarise(dBIC = stats::median(Delta_BIC_vs_PC1_EVS), .groups = "drop")
  lvA <- c("Pooled hierarchical", "Continuous rank", "Variance split", "PC1-PC2 split", "Random partition")
  colA <- c("#2B2B2B", "#CC79A7", "#E69F00", "#009E73", "#56B4E9"); names(colA) <- lvA
  da <- da[da$Model %in% lvA, , drop = FALSE]
  da$Model <- factor(da$Model, levels = lvA); da$Comparison <- as_comp(da$Comparison)
  pa <- ggplot(da, aes(Comparison, dBIC, fill = Model)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_col(position = position_dodge(width = 0.82), width = 0.76) +
    scale_fill_manual(values = colA, name = NULL, drop = TRUE) +
    scale_x_discrete(labels = function(x) sub("_", "\n", x)) + scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "BIC minus BIC of PC1 EVS", tag = "a", title = "Model fit (descriptive)",
         subtitle = "Positive: PC1 EVS fits the dispersions better. Random partition: median.") +
    fig_theme() + theme(panel.grid.major.x = element_blank())
  rows <- list(pa)
  mc <- if (length(model_comparison_rows)) dplyr::bind_rows(model_comparison_rows) else NULL
  if (!is.null(mc) && nrow(mc)) {
    mc$Model2 <- ifelse(grepl("^Random partition", mc$Model), "Random partition", mc$Model)
    lvB <- c("Conventional DESeq2", "Pooled hierarchical", "PC1 EVS", "Variance split", "PC1-PC2 split",
             "Continuous rank", "Random partition")
    colB <- c("#7F7F7F", "#2B2B2B", "#1F6FB2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9"); names(colB) <- lvB
    db <- mc |> dplyr::group_by(Comparison, Model2) |> dplyr::summarise(Standard = stats::median(Standard), .groups = "drop")
    db$Model2 <- factor(db$Model2, levels = lvB); db$Comparison <- as_comp(db$Comparison)
    rows[[2]] <- ggplot(db, aes(Comparison, Standard, fill = Model2)) +
      geom_col(position = position_dodge(width = 0.86), width = 0.8) +
      scale_fill_manual(values = colB, name = NULL, drop = TRUE) +
      scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
      scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.08))) +
      labs(x = NULL, y = "Standard discoveries", tag = "b", title = "Complete DE analysis of every model",
           subtitle = "Counts are not error rates: FDR and power are assessed by simulation (Figure 3).") +
      fig_theme() + theme(panel.grid.major.x = element_blank())
  }
  sens <- if (length(sensitivity_rows)) dplyr::bind_rows(sensitivity_rows) else NULL
  if (!is.null(sens) && nrow(sens)) {
    lv <- c(ORIGINAL_LABEL, "EVS (primary)", paste0(ORIGINAL_LABEL, ", DESeq2 0.25 prior floor"),
            "EVS, DESeq2 0.25 prior floor", "IHW, Original p-values")
    sens <- sens[sens$Analysis %in% lv & is.finite(sens$Standard), , drop = FALSE]
    sens$Analysis <- factor(sens$Analysis, levels = lv); sens$Comparison <- as_comp(sens$Comparison)
    cols <- c("#2B2B2B", "#1F6FB2", "#8C8C8C", "#8FB8DE", "#E69F00"); names(cols) <- lv
    rows[[length(rows) + 1L]] <- ggplot(sens, aes(Comparison, Standard, fill = Analysis)) +
      geom_col(position = position_dodge(width = 0.86), width = 0.8) +
      scale_fill_manual(values = cols, name = NULL, drop = TRUE) +
      scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
      scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.08))) +
      labs(x = NULL, y = "Standard discoveries", tag = letters[length(rows) + 1L],
           title = "Sensitivity: DESeq2 0.25 prior floor and IHW (mean-count covariate)") +
      fig_theme() + theme(panel.grid.major.x = element_blank())
  }
  assemble_rows(rows, heights = rep(1, length(rows)), title = "Supplementary Figure S1. Model fit (descriptive BIC), discoveries of every model and sensitivity analyses")
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
# ALTERNATIVE POLYADENYLATION (APA): RELATIVE PAS USAGE, POSITIONAL SHIFTS
# AND PROXIMAL / DISTAL SWITCHING
# =============================================================================
# Differential PAS abundance does not establish a change in which PAS a gene
# uses, so APA is analysed as a separate outcome, per pairwise comparison, for
# all eligible multi-PAS genes. Gene units are gene symbol + chromosome +
# strand; PASs are ordered along the direction of transcription by their peak
# coordinate (proximal = most 5'). Details of each analysis are given above
# run_apa_analysis. APA_MIN_MEAN: minimum mean normalized count of each PAS of
# a proximal/distal pair; APA_FDR: BH level of all APA tests.
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

# ---- Transcript annotation (optional, item 6) --------------------------------
# A rat GTF matching the genome build of the PAS coordinates can be supplied
# with SEQUENCE_GTF (plain or .gz; Ensembl or NCBI). Each coding transcript's
# 3'UTR is taken as the interval from the end of its CDS to the end of the
# transcript (strand aware). A PAS is "3'UTR" when it lies in the 3'UTR of any
# transcript, "Gene body" when it lies elsewhere within an annotated gene, and
# "Outside annotated genes" otherwise; the gene_id of the overlapping gene is
# its stable identifier. Without a GTF, APA is reported as distal/proximal
# switching only, never as confirmed 3'UTR lengthening or shortening.
GTF_FILE <- Sys.getenv("SEQUENCE_GTF", "")

load_gtf_annotation <- function() {
  if (!nzchar(GTF_FILE)) return(NULL)
  if (!file.exists(GTF_FILE)) stop("SEQUENCE_GTF does not exist: ", GTF_FILE)
  con <- if (grepl("\\.gz$", GTF_FILE)) gzfile(GTF_FILE, "r") else file(GTF_FILE, "r")
  on.exit(close(con), add = TRUE)
  keep <- list()
  repeat {
    x <- readLines(con, n = 200000L, warn = FALSE)
    if (!length(x)) break
    x <- x[!startsWith(x, "#") & grepl("\t(transcript|CDS)\t", x)]
    if (length(x)) keep[[length(keep) + 1L]] <- x
  }
  x <- unlist(keep, use.names = FALSE)
  if (!length(x)) stop("No transcript or CDS records in the GTF.")
  f <- strsplit(x, "\t", fixed = TRUE)
  fld <- function(i) vapply(f, function(z) if (length(z) >= i) z[[i]] else NA_character_, "")
  g <- data.frame(chr = sub("^chr", "", fld(1L)), type = fld(3L), start = as.numeric(fld(4L)), end = as.numeric(fld(5L)),
                  strand = fld(7L), attr = fld(9L), stringsAsFactors = FALSE)
  ga <- function(key) {
    re <- paste0(key, ' "([^"]*)"'); hit <- regexpr(re, g$attr)
    out <- rep(NA_character_, nrow(g)); out[hit > 0] <- sub(re, "\\1", regmatches(g$attr, hit)); out
  }
  g$gene_id <- ga("gene_id"); g$transcript_id <- ga("transcript_id"); g$gene_name <- ga("gene_name")
  g <- g[g$strand %in% c("+", "-") & is.finite(g$start) & is.finite(g$end) & !is.na(g$transcript_id), , drop = FALSE]
  tx <- g[g$type == "transcript", , drop = FALSE]; cds <- g[g$type == "CDS", , drop = FALSE]
  if (!nrow(tx)) stop("The GTF has no transcript records.")
  cmin <- tapply(cds$start, cds$transcript_id, min); cmax <- tapply(cds$end, cds$transcript_id, max)
  tx$cds_min <- as.numeric(cmin[tx$transcript_id]); tx$cds_max <- as.numeric(cmax[tx$transcript_id])
  u <- tx[is.finite(tx$cds_min), , drop = FALSE]; plus <- u$strand == "+"
  utr <- data.frame(chr = u$chr, strand = u$strand, gene_id = u$gene_id,
                    start = ifelse(plus, u$cds_max + 1, u$start), end = ifelse(plus, u$end, u$cds_min - 1))
  utr <- utr[utr$end >= utr$start, , drop = FALSE]
  key <- paste(tx$gene_id, tx$chr, tx$strand, sep = "\r")
  gs <- tapply(tx$start, key, min); ge <- tapply(tx$end, key, max)
  parts <- strsplit(names(gs), "\r", fixed = TRUE)
  genes <- data.frame(gene_id = vapply(parts, `[`, "", 1L), chr = vapply(parts, `[`, "", 2L),
                      strand = vapply(parts, `[`, "", 3L), start = as.numeric(gs), end = as.numeric(ge[names(gs)]))
  list(utr = utr, genes = genes, file = GTF_FILE)
}

annotate_pas_regions <- function(ann, gtf) {
  gr <- GenomicRanges::GRanges(sub("^chr", "", ann$chr), IRanges::IRanges(as.integer(round(ann$pos)), width = 1L),
                               strand = ann$strand)
  ur <- GenomicRanges::GRanges(gtf$utr$chr, IRanges::IRanges(as.integer(gtf$utr$start), as.integer(gtf$utr$end)),
                               strand = gtf$utr$strand)
  gn <- GenomicRanges::GRanges(gtf$genes$chr, IRanges::IRanges(as.integer(gtf$genes$start), as.integer(gtf$genes$end)),
                               strand = gtf$genes$strand)
  og <- suppressWarnings(IRanges::findOverlaps(gr, gn))
  ou <- suppressWarnings(IRanges::findOverlaps(gr, ur))
  gq <- S4Vectors::queryHits(og); gh <- S4Vectors::subjectHits(og)
  uq <- S4Vectors::queryHits(ou); uh <- S4Vectors::subjectHits(ou)
  ann$gene_id <- NA_character_; ann$utr_gene_id <- NA_character_; ann$region <- "Outside annotated genes"
  ann$gene_id[gq[!duplicated(gq)]] <- gtf$genes$gene_id[gh[!duplicated(gq)]]
  ann$region[unique(gq)] <- "Gene body"
  ann$utr_gene_id[uq[!duplicated(uq)]] <- gtf$utr$gene_id[uh[!duplicated(uq)]]
  ann$region[unique(uq)] <- "3'UTR"
  ann$gene_id <- ifelse(is.na(ann$utr_gene_id), ann$gene_id, ann$utr_gene_id)
  ann
}

# ---- Relative PAS usage and the all-PAS positional index (item 6) ---------------
# For every gene with at least two eligible PASs: the share of the gene's
# normalized reads at each PAS in each animal, its mean in each condition and the
# absolute change in percentage points; and a positional index, the
# usage-weighted mean relative position of the gene's PASs (0 = most proximal,
# 1 = most distal), which captures switching among more than two PASs. The
# index is compared between conditions with Welch's t-test (BH across genes).
# Descriptive measures of the same usage that diffSpliceDGE tests.
apa_usage_metrics <- function(cmp, y, a, cache) {
  cond <- as.character(cache$coldata[colnames(y), "condition"]); trt <- cond == "trt"
  nc <- sweep(y, 2, cache$sf[colnames(y)], "/")
  unit <- a$gene_unit
  tot <- rowsum(nc, unit)
  use <- nc / tot[unit, , drop = FALSE]
  use[!is.finite(use)] <- NA_real_
  ord <- stats::ave(a$tx_pos, unit, FUN = function(v) rank(v, ties.method = "first"))
  npas <- stats::ave(a$tx_pos, unit, FUN = length)
  rel <- (ord - 1) / pmax(npas - 1, 1)
  pas <- data.frame(Comparison = cmp, feature_id = rownames(y), Gene = a$gene, Gene_unit = unit,
                    PAS_order_5to3 = as.integer(ord), PAS_in_gene = as.integer(npas), Relative_position = rel,
                    Usage_control = rowMeans(use[, !trt, drop = FALSE], na.rm = TRUE),
                    Usage_recovery = rowMeans(use[, trt, drop = FALSE], na.rm = TRUE))
  pas$Delta_usage_pp <- 100 * (pas$Usage_recovery - pas$Usage_control)
  # Welch 95% interval of the usage change (percentage points), per PAS.
  welch_ci <- function(a1, a0) {
    a1 <- a1[is.finite(a1)]; a0 <- a0[is.finite(a0)]
    if (length(a1) < 2L || length(a0) < 2L || stats::var(c(a1, a0)) == 0) return(c(NA_real_, NA_real_))
    tryCatch(as.numeric(stats::t.test(a1, a0)$conf.int), error = function(e) c(NA_real_, NA_real_))
  }
  ci <- t(vapply(seq_len(nrow(use)), function(i) 100 * welch_ci(use[i, trt], use[i, !trt]), numeric(2)))
  pas$Delta_usage_pp_lo95 <- ci[, 1]; pas$Delta_usage_pp_hi95 <- ci[, 2]
  idx <- rowsum(use * rel, unit)
  gene <- do.call(rbind, lapply(rownames(idx), function(g) {
    v <- idx[g, ]; a1 <- v[trt]; a0 <- v[!trt]; a1 <- a1[is.finite(a1)]; a0 <- a0[is.finite(a0)]
    tt <- if (length(a1) >= 2L && length(a0) >= 2L && stats::var(c(a1, a0)) > 0)
      tryCatch(stats::t.test(a1, a0), error = function(e) NULL) else NULL
    data.frame(Gene_unit = g, Index_control = if (length(a0)) mean(a0) else NA_real_,
               Index_recovery = if (length(a1)) mean(a1) else NA_real_,
               Delta_index_lo95 = if (is.null(tt)) NA_real_ else tt$conf.int[1],
               Delta_index_hi95 = if (is.null(tt)) NA_real_ else tt$conf.int[2],
               P = if (is.null(tt)) NA_real_ else tt$p.value)
  }))
  gene$Comparison <- cmp
  gene$Gene <- sub("\\|.*$", "", gene$Gene_unit)
  gene$Delta_index <- gene$Index_recovery - gene$Index_control
  gene$FDR <- stats::p.adjust(gene$P, "BH")
  mx <- tapply(abs(pas$Delta_usage_pp), pas$Gene_unit, function(v) if (all(!is.finite(v))) NA_real_ else max(v, na.rm = TRUE))
  gene$Max_abs_delta_usage_pp <- as.numeric(mx[gene$Gene_unit])
  list(pas = pas, gene = gene[, c("Comparison", "Gene", "Gene_unit", "Index_control", "Index_recovery", "Delta_index",
                                  "Delta_index_lo95", "Delta_index_hi95", "Max_abs_delta_usage_pp", "P", "FDR")])
}

# ---- APA analysis (independent of differential PAS abundance) -------------------
# Every gene with at least two eligible PASs is analysed, whether or not any of
# its PASs changed in abundance; APA and abundance are separate outcomes.
#  (i) usage change: edgeR quasi-likelihood GLMs (robust dispersion; Lund et al.
#      2012; Chen et al. 2016) with offsets from the shared DESeq2 size
#      factors; diffSpliceDGE tests each PAS against the gene's other PASs;
#      gene-level Simes p-values with BH; stage-wise PAS confirmation with
#      stageR (Van den Berge et al. 2017) when installed;
# (ii) direction: the two most used PASs (mean normalized count >= 5 each) form
#      a proximal/distal pair; RED = treatment log2FC of the distal minus that of
#      the proximal PAS (Hoque et al. 2013), tested by the gene-level
#      diffSpliceDGE test of the pair, BH across genes. RED > 0 is a distal
#      shift, RED < 0 a proximal shift; with a transcript annotation, shifts
#      between two 3'UTR PASs of the same gene are labelled 3'UTR lengthening
#      or shortening;
#(iii) relative usage (percentage points) and the positional index (above).
run_apa_analysis <- function() {
  ann <- pas_annotation()
  gtf <- tryCatch(load_gtf_annotation(), error = function(e) {
    warning("GTF annotation not used: ", conditionMessage(e), call. = FALSE); NULL })
  if (!is.null(gtf)) ann <- annotate_pas_regions(ann, gtf)
  pas_out <- list(); gene_out <- list(); red_out <- list(); stage_out <- list(); use_out <- list(); idx_out <- list()
  for (cmp in COMPARISON_NAMES) {
    cache <- comparison_cache[[cmp]]; y <- cache$raw
    a <- ann[match(rownames(y), as.character(ann$feature_id)), , drop = FALSE]
    ok <- !is.na(a$gene_unit)
    n_per <- table(a$gene_unit[ok]); multi <- ok & a$gene_unit %in% names(n_per[n_per >= 2])
    le <- CUTOFF$pas[[cmp]]; le_ids <- if (is.null(le)) character() else le$feature_id[le$In_Leading_Edge]
    if (sum(multi) >= 2L) {
      f <- apa_fit(y[multi, , drop = FALSE],
                   data.frame(gene_unit = a$gene_unit[multi], Gene = a$gene[multi], PAS = rownames(y)[multi]), cache)
      pt <- edgeR::topSpliceDGE(f$ds, test = "exon", number = Inf, FDR = 1); pt$Comparison <- cmp
      gt <- edgeR::topSpliceDGE(f$ds, test = "Simes", number = Inf, FDR = 1); gt$Comparison <- cmp
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
      um <- apa_usage_metrics(cmp, y[multi, , drop = FALSE], a[multi, , drop = FALSE], cache)
      pt$Delta_usage_pp <- um$pas$Delta_usage_pp[match(as.character(pt$PAS), um$pas$feature_id)]
      use_out[[cmp]] <- um$pas; idx_out[[cmp]] <- um$gene
      pas_out[[cmp]] <- pt; gene_out[[cmp]] <- gt
    }
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
    red$Direction <- ifelse(red$RED_log2 > 0, "Distal shift", "Proximal shift")
    if (!is.null(gtf)) {
      pr <- ann[match(red$Proximal_PAS, as.character(ann$feature_id)), , drop = FALSE]
      dr <- ann[match(red$Distal_PAS, as.character(ann$feature_id)), , drop = FALSE]
      red$Gene_id <- pr$gene_id; red$Proximal_region <- pr$region; red$Distal_region <- dr$region
      utr_pair <- pr$region %in% "3'UTR" & dr$region %in% "3'UTR" & !is.na(pr$utr_gene_id) & pr$utr_gene_id == dr$utr_gene_id
      red$Event <- ifelse(utr_pair, ifelse(red$RED_log2 > 0, "3'UTR lengthening", "3'UTR shortening"),
                          paste(red$Direction, "(not a confirmed 3'UTR pair)"))
    } else {
      red$Gene_id <- NA_character_; red$Proximal_region <- NA_character_; red$Distal_region <- NA_character_
      red$Event <- paste(red$Direction, "(no transcript annotation)")
    }
    red$Proximal_in_LE <- red$Proximal_PAS %in% le_ids; red$Distal_in_LE <- red$Distal_PAS %in% le_ids
    red_out[[cmp]] <- red
  }
  pas_tab <- dplyr::bind_rows(pas_out); gene_tab <- dplyr::bind_rows(gene_out); red_tab <- dplyr::bind_rows(red_out)
  use_tab <- dplyr::bind_rows(use_out); idx_tab <- dplyr::bind_rows(idx_out)
  save_csv(pas_tab, file.path(summary_table_dir, "Table_APA_PAS_Usage.csv"))
  save_csv(gene_tab, file.path(summary_table_dir, "Table_APA_Gene_Usage_Simes.csv"))
  save_csv(red_tab, file.path(summary_table_dir, "Table_APA_Proximal_Distal.csv"))
  save_csv(use_tab, file.path(summary_table_dir, "Table_APA_PAS_Usage_Shares.csv"))
  save_csv(idx_tab, file.path(summary_table_dir, "Table_APA_Positional_Index.csv"))
  save_csv(if (length(stage_out)) dplyr::bind_rows(stage_out) else data.frame(Note = "stageR not installed or not run"),
           file.path(summary_table_dir, "Table_APA_StageR.csv"))
  summ <- do.call(rbind, lapply(COMPARISON_NAMES, function(cmp) {
    r <- red_tab[red_tab$Comparison == cmp, , drop = FALSE]; g <- gene_tab[gene_tab$Comparison == cmp, , drop = FALSE]
    ix <- idx_tab[idx_tab$Comparison == cmp, , drop = FALSE]
    data.frame(Comparison = cmp, Genes_usage_tested = nrow(g),
               Genes_usage_changed = if (nrow(g)) sum(g$FDR < APA_FDR, na.rm = TRUE) else 0L,
               Genes_positional_shift = if (nrow(ix)) sum(ix$FDR < APA_FDR, na.rm = TRUE) else 0L,
               Genes_pair_tested = nrow(r),
               Proximal_shift = sum(r$Significant & r$Direction == "Proximal shift"),
               Distal_shift = sum(r$Significant & r$Direction == "Distal shift"),
               UTR_shortening_confirmed = if (!is.null(gtf)) sum(r$Significant & r$Event == "3'UTR shortening") else NA_integer_,
               UTR_lengthening_confirmed = if (!is.null(gtf)) sum(r$Significant & r$Event == "3'UTR lengthening") else NA_integer_,
               Significant_with_PAS_in_LE = sum(r$Significant & (r$Proximal_in_LE | r$Distal_in_LE)),
               Transcript_annotation = if (is.null(gtf)) "none" else basename(gtf$file))
  }))
  save_csv(summ, file.path(summary_table_dir, "Table_APA_Summary.csv"))
  APA_RESULTS <<- list(pas = pas_tab, gene = gene_tab, red = red_tab, usage = use_tab, index = idx_tab,
                       summary = summ, enrichment = NULL)
  invisible(APA_RESULTS)
}

# Figure 6: (a) genes with a usage change, a positional shift and a proximal or
# distal shift per recovery time; (b) RED volcano per comparison; (c) the
# primary 3'aTWAS endpoint (Firth odds ratio and 95% interval per comparison).
plot_apa_figure <- function() {
  s <- APA_RESULTS$summary; red <- APA_RESULTS$red
  tl <- stats::setNames(unname(TIME_LEVELS), COMPARISON_NAMES)
  d <- rbind(data.frame(Comparison = s$Comparison, Class = "Any usage change (gene, Simes)", n = s$Genes_usage_changed),
             data.frame(Comparison = s$Comparison, Class = "Positional shift (all PASs)", n = s$Genes_positional_shift),
             data.frame(Comparison = s$Comparison, Class = "Proximal shift", n = s$Proximal_shift),
             data.frame(Comparison = s$Comparison, Class = "Distal shift", n = s$Distal_shift))
  d$Time <- factor(tl[d$Comparison], levels = unname(TIME_LEVELS))
  d$Class <- factor(d$Class, levels = c("Any usage change (gene, Simes)", "Positional shift (all PASs)", "Proximal shift", "Distal shift"))
  cols <- c("Any usage change (gene, Simes)" = "grey40", "Positional shift (all PASs)" = "#7B3294",
            "Proximal shift" = "#C0392B", "Distal shift" = "#1F6FB2")
  pa <- ggplot(d, aes(Time, n, fill = Class)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.74) +
    geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.3, size = 1.9) +
    scale_fill_manual(values = cols, name = NULL) + scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Recovery time", y = sprintf("Genes (FDR < %.2f)", APA_FDR), tag = "a",
         title = "Alternative polyadenylation in recovery",
         subtitle = "Any change and positional shift: all PASs of the gene. Proximal / distal: two most used PASs (RED).") +
    fig_theme() + theme(panel.grid.major.x = element_blank())
  red$Time <- factor(tl[red$Comparison], levels = unname(TIME_LEVELS))
  red$Class <- ifelse(red$Significant, red$Direction, "Not significant")
  red$Class <- factor(red$Class, levels = c("Proximal shift", "Distal shift", "Not significant"))
  red <- red[is.finite(red$RED_log2) & is.finite(red$FDR), , drop = FALSE]
  red <- red[order(red$Class, decreasing = TRUE), , drop = FALSE]
  pb <- ggplot(red, aes(RED_log2, -log10(pmax(FDR, 1e-300)), colour = Class)) +
    geom_point(size = 0.5, alpha = 0.6, stroke = 0) +
    geom_hline(yintercept = -log10(APA_FDR), linetype = "22", colour = "grey45", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    facet_wrap(~ Time, nrow = 1) +
    scale_colour_manual(values = c("Proximal shift" = "#C0392B", "Distal shift" = "#1F6FB2", "Not significant" = "grey75"), name = NULL) +
    labs(x = "RED: change in log2 distal / proximal usage (recovery vs control)", y = "-log10 FDR", tag = "b",
         title = "Direction of APA changes",
         subtitle = "Positive: more distal PAS use; negative: more proximal. 3'UTR lengthening/shortening requires a transcript annotation.") +
    fig_theme()
  rows <- list(pa, pb)
  en <- APA_RESULTS$enrichment
  # Every estimable comparison is plotted (Firth estimates are finite for zero overlaps).
  if (!is.null(en)) en <- en[is.finite(en$OR) & is.finite(en$Lower95) & is.finite(en$Upper95), , drop = FALSE]
  if (!is.null(en) && nrow(en)) {
    en$Time <- factor(tl[en$Comparison], levels = unname(TIME_LEVELS))
    rows[[3]] <- ggplot(en, aes(Time, OR)) +
      geom_hline(yintercept = 1, linetype = "22", colour = "grey45", linewidth = 0.3) +
      geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0, linewidth = 0.4) +
      geom_point(size = 2, colour = "#7B3294") +
      geom_text(aes(label = sprintf("%d genes", Overlap)), vjust = -1, size = 1.9) +
      scale_y_log10() +
      labs(x = "Recovery time", y = "Odds ratio (95% CI)", tag = "c",
           title = "Primary 3'aTWAS endpoint: genes with differential APA",
           subtitle = "Firth logistic regression; background: ortholog-mapped genes in the usage test; adjusted for PASs, expression, ortholog support.") +
      fig_theme()
  }
  assemble_rows(rows, heights = rep(1, length(rows)), title = "Figure 6. Alternative polyadenylation and relative PAS usage")
}

# =============================================================================
# FIGURE 6: APA DETAIL. (a) PAS usage of example genes, (b) RED over recovery
# time, (c) APA vs PAS abundance change. Descriptive; uses APA_RESULTS, the
# comparison caches and the selected primary abundance calls.
# =============================================================================
apa_gene_abundance <- function(ann) {
  sp <- selected_primary[selected_primary$standard_flag %in% TRUE, c("Comparison", "feature_id"), drop = FALSE]
  sp$unit <- ann$gene_unit[match(as.character(sp$feature_id), as.character(ann$feature_id))]
  out <- list()
  for (cmp in COMPARISON_NAMES) {
    cache <- comparison_cache[[cmp]]; y <- cache$raw
    a <- ann[match(rownames(y), as.character(ann$feature_id)), , drop = FALSE]
    ok <- !is.na(a$gene_unit)
    nc <- sweep(y[ok, , drop = FALSE], 2, cache$sf[colnames(y)], "/")
    g <- rowsum(nc, a$gene_unit[ok])
    trt <- as.character(cache$coldata[colnames(y), "condition"]) == "trt"
    lfc <- log2((rowMeans(g[, trt, drop = FALSE]) + 0.5) / (rowMeans(g[, !trt, drop = FALSE]) + 0.5))
    out[[cmp]] <- data.frame(Comparison = cmp, Gene_unit = rownames(g), Gene_log2FC = unname(lfc),
                             Abundance_changed = rownames(g) %in% sp$unit[sp$Comparison == cmp])
  }
  dplyr::bind_rows(out)
}

apa_usage_long <- function(ex, ann) {
  tl <- stats::setNames(unname(TIME_LEVELS), COMPARISON_NAMES)
  out <- list()
  for (i in seq_len(nrow(ex))) {
    cmp <- ex$Comparison[i]; cache <- comparison_cache[[cmp]]; y <- cache$raw
    a <- ann[ann$gene_unit == ex$Gene_unit[i] & as.character(ann$feature_id) %in% rownames(y), , drop = FALSE]
    if (nrow(a) < 2L) next
    a <- a[order(a$tx_pos), , drop = FALSE]
    ids <- as.character(a$feature_id)
    nc <- sweep(y[ids, , drop = FALSE], 2, cache$sf[colnames(y)], "/")
    tot <- colSums(nc)
    use <- sweep(nc, 2, ifelse(tot > 0, tot, NA_real_), "/")
    role <- ifelse(ids == as.character(ex$Proximal_PAS[i]), " prox", ifelse(ids == as.character(ex$Distal_PAS[i]), " dist", ""))
    pas_lab <- paste0("P", seq_along(ids), role)
    cond <- ifelse(as.character(cache$coldata[colnames(y), "condition"]) == "trt", "Recovery", "Control")
    out[[i]] <- data.frame(
      Label = sprintf("%s, %s\nRED %+.2f, FDR %.2g", ex$Gene[i], tl[[cmp]], ex$RED_log2[i], ex$FDR[i]),
      PAS = factor(rep(pas_lab, times = ncol(use)), levels = pas_lab),
      Condition = factor(rep(cond, each = length(ids)), levels = c("Control", "Recovery")),
      Usage = as.vector(use))
  }
  dplyr::bind_rows(out)
}

plot_apa_detail_figure <- function(n_examples = 6L, n_heat = 60L) {
  red <- APA_RESULTS$red
  if (is.null(red) || !nrow(red) || !any(red$Significant %in% TRUE)) stop("no significant APA genes")
  ann <- pas_annotation()
  tl <- stats::setNames(unname(TIME_LEVELS), COMPARISON_NAMES)
  ccols <- c("Control" = "grey45", "Recovery" = "#D55E00")
  sig <- red[red$Significant %in% TRUE, , drop = FALSE]

  # (a) Example genes: APA genes that are 3'aTWAS genes first, then smallest FDR.
  en <- APA_RESULTS$enrichment
  twas_genes <- if (!is.null(en) && "Overlap_genes" %in% names(en))
    unique(unlist(strsplit(as.character(stats::na.omit(en$Overlap_genes)), ";", fixed = TRUE))) else character()
  ex <- sig[order(!(sig$Gene %in% twas_genes), sig$FDR), , drop = FALSE]
  ex <- utils::head(ex[!duplicated(ex$Gene_unit), , drop = FALSE], n_examples)
  u <- apa_usage_long(ex, ann)
  u <- u[is.finite(u$Usage), , drop = FALSE]
  if (!nrow(u)) stop("no example genes with two or more PASs")
  u$PAS <- as.character(u$PAS)
  u$PAS <- factor(u$PAS, levels = unique(u$PAS[order(as.integer(sub("^P([0-9]+).*$", "\\1", u$PAS)), u$PAS)]))
  pd <- position_dodge(width = 0.6)
  pa <- ggplot(u, aes(PAS, Usage, colour = Condition)) +
    stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.5,
                 linewidth = 0.3, position = pd) +
    geom_point(position = pd, size = 0.8, alpha = 0.8) +
    facet_wrap(~ Label, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = ccols, name = NULL) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    labs(x = "PASs of the gene, 5' to 3' (prox / dist: pair used for RED)", y = "Share of gene reads", tag = "a",
         title = "PAS usage of example APA genes",
         subtitle = "Points: rats; bars: group means. 3'aTWAS genes shown first, then smallest FDR.") +
    fig_theme()

  # (b) RED over recovery time for genes significant at any time.
  best <- sig[order(sig$FDR), , drop = FALSE]
  units <- utils::head(unique(best$Gene_unit), n_heat)
  h <- red[red$Gene_unit %in% units & is.finite(red$RED_log2), c("Comparison", "Gene_unit", "RED_log2", "Significant"), drop = FALSE]
  first <- sig[sig$Gene_unit %in% units, , drop = FALSE]
  first <- first[order(match(first$Comparison, COMPARISON_NAMES), -first$RED_log2), , drop = FALSE]
  lev <- unique(first$Gene_unit)
  h$Time <- factor(tl[h$Comparison], levels = unname(TIME_LEVELS))
  h$Gene <- factor(h$Gene_unit, levels = rev(lev))
  h$RED_clip <- pmax(pmin(h$RED_log2, 4), -4)
  pb <- ggplot(h, aes(Time, Gene, fill = RED_clip)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_point(data = h[h$Significant %in% TRUE, , drop = FALSE], size = 0.5, colour = "black") +
    scale_fill_gradient2(low = "#C0392B", mid = "grey95", high = "#1F6FB2", midpoint = 0, limits = c(-4, 4),
                         name = "RED\n(log2)") +
    scale_y_discrete(labels = function(x) sub("\\|.*", "", x)) +
    labs(x = "Recovery time", y = NULL, tag = "b", title = "APA over recovery time",
         subtitle = sprintf("Genes significant at any time (up to %d). Dot: FDR < %.2f. Blank: not tested.", n_heat, APA_FDR)) +
    fig_theme() + theme(axis.text.y = element_text(size = 4), panel.grid = element_blank())

  # (c) APA vs abundance change of the gene.
  ab <- apa_gene_abundance(ann)
  r <- merge(red[is.finite(red$RED_log2), , drop = FALSE], ab, by = c("Comparison", "Gene_unit"), all.x = TRUE)
  r$Class <- ifelse(!(r$Significant %in% TRUE), "Not significant",
                    ifelse(r$Abundance_changed %in% TRUE, "APA + PAS abundance change", "APA only"))
  r$Class <- factor(r$Class, levels = c("Not significant", "APA only", "APA + PAS abundance change"))
  r$Time <- factor(tl[r$Comparison], levels = unname(TIME_LEVELS))
  r <- r[order(r$Class), , drop = FALSE]
  save_csv(r[r$Significant %in% TRUE, c("Comparison", "Gene", "Gene_unit", "RED_log2", "FDR", "Direction",
                                         "Gene_log2FC", "Abundance_changed"), drop = FALSE],
           file.path(summary_table_dir, "Table_APA_Abundance_Overlap.csv"))
  lab <- do.call(rbind, lapply(split(r, r$Time), function(z) if (nrow(z)) data.frame(
    Time = z$Time[1], txt = sprintf("%d APA genes\n%d also abundance", sum(z$Significant %in% TRUE),
                                     sum(z$Class == "APA + PAS abundance change")))))
  pc <- ggplot(r, aes(Gene_log2FC, RED_log2, colour = Class)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_point(size = 0.6, alpha = 0.7, stroke = 0, na.rm = TRUE) +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
              hjust = -0.05, vjust = 1.1, size = 1.8, lineheight = 0.9) +
    facet_wrap(~ Time, nrow = 1, drop = FALSE) +
    scale_colour_manual(values = c("Not significant" = "grey80", "APA only" = "#7B3294",
                                   "APA + PAS abundance change" = "#E69F00"), name = NULL) +
    labs(x = "Gene log2 fold change (all PASs, recovery vs control)", y = "RED (log2)", tag = "c",
         title = "APA vs abundance change",
         subtitle = "Abundance change: any PAS of the gene significant (Standard, selected model).") +
    fig_theme()

  assemble_rows(list(pa, list(plots = list(pb, pc), widths = c(0.8, 2))), heights = c(1, 1.4),
                title = "Supplementary Figure S2. APA detail: example genes, time course and relation to abundance")
}

# Gene-level 3'aTWAS overlap enrichment. Map the entire tested rat background
# to human genes, not just TWAS hits. Collapse mapping to one rat gene per row.
# Logistic adjustment for mean abundance and number of tested PASs reduces
# differential detectability bias; one-to-many mappings are reported separately.
# This is exploratory enrichment, not causal or cross-species validation.
# =============================================================================
# 3'aTWAS GENE-SET ENRICHMENT (item 7)
# =============================================================================
# Primary endpoint: differential APA (gene-level usage change, Simes FDR <
# APA_FDR), because 3'aTWAS genes are defined by genetically regulated APA.
# Secondary endpoint: differential PAS abundance (any Standard PAS in the gene),
# for the Original and for EVS on one identical gene background.
# Background: database-supported human-to-rat orthologs (babelgene, at least
# TWAS_ORTHOLOG_MIN_SUPPORT source) among the genes tested for the endpoint: for
# abundance, genes with a finite p-value in both the Original and EVS; for APA,
# the genes of the usage test. Covariates: log(1 + PASs tested), log(1 + mean
# normalized count) and log(number of supporting ortholog databases).
# Estimation: Firth's penalized-likelihood logistic regression (Firth 1993;
# Heinze & Schemper 2002; R package logistf), which gives finite estimates for
# sparse and zero overlaps, with profile penalized-likelihood 95% intervals. If
# the Firth fit fails or returns non-finite estimates, the row is reported as
# not estimable (NA estimates). Sensitivity: one-to-one orthologs. Threshold-free test: limma
# cameraPR (Wu & Smyth 2012) of the gene statistic z = qnorm(1 - p/2), with p
# the gene-level p-value (Simes over the gene's PASs for abundance; the Simes
# usage p-value for APA). Every comparison is reported, including zero
# overlaps, with BH across rows within each method. TWAS results never inform
# any EVS setting.
simes_p <- function(p) { p <- sort(p[is.finite(p)]); if (!length(p)) NA_real_ else min(1, min(p * length(p) / seq_along(p))) }

twas_logistic <- function(g, base) {
  out <- function(OR = NA_real_, lo = NA_real_, hi = NA_real_, P = NA_real_, method = NA_character_, note = "")
    cbind(base, data.frame(Method = method, OR = OR, Lower95 = lo, Upper95 = hi, P = P, Direction = NA_character_, Note = note))
  if (nrow(g) < 10L) return(out(note = "not estimable: fewer than 10 background genes"))
  if (!any(g$TWAS) || all(g$TWAS)) return(out(note = "not estimable: no TWAS variation in the background"))
  if (!any(g$Y) || all(g$Y)) return(out(note = "not estimable: no outcome variation"))
  g$TWAS01 <- as.numeric(g$TWAS); g$Y01 <- as.numeric(g$Y)
  form <- if (length(unique(g$support)) > 1L) Y01 ~ TWAS01 + log1p(nPAS) + log1p(mean_count) + log(support) else
    Y01 ~ TWAS01 + log1p(nPAS) + log1p(mean_count)
  f <- tryCatch(logistf::logistf(form, data = g, plconf = 2L), error = function(e) NULL)
  if (is.null(f)) return(out(method = "Firth logistic (profile CI)", note = "not estimable: Firth fit failed"))
  b <- unname(f$coefficients["TWAS01"]); lo <- unname(f$ci.lower["TWAS01"]); hi <- unname(f$ci.upper["TWAS01"])
  pv <- unname(f$prob["TWAS01"])
  if (!all(is.finite(c(b, lo, hi, pv))))
    return(out(method = "Firth logistic (profile CI)", note = "not estimable: non-finite Firth estimates"))
  out(exp(b), exp(lo), exp(hi), pv, "Firth logistic (profile CI)")
}

twas_camera <- function(g, base) {
  z <- stats::qnorm(pmax(g$gene_p, 1e-300) / 2, lower.tail = FALSE)
  idx <- which(g$TWAS)
  r <- if (length(idx) >= 2L && length(idx) < nrow(g)) tryCatch(
    limma::cameraPR(stats::setNames(z, g$gene_key), list(TWAS = idx)), error = function(e) NULL) else NULL
  cbind(base, data.frame(Method = "cameraPR (threshold-free)", OR = NA_real_, Lower95 = NA_real_, Upper95 = NA_real_,
                         P = if (is.null(r)) NA_real_ else r$PValue[1],
                         Direction = if (is.null(r)) NA_character_ else as.character(r$Direction[1]),
                         Note = if (is.null(r)) "not estimable" else ""))
}

run_twas_enrichment <- function(mapping) {
  sym <- stats::setNames(OrigID_Symbol$gene_symbol, OrigID_Symbol$feature_id)
  # Abundance gene tables (Original and EVS on the identical background).
  ab <- lapply(COMPARISON_NAMES, function(cmp) {
    fl <- MODEL_FLAGS[[cmp]]; o <- fl[["Pooled hierarchical"]]; e <- fl[["PC1 EVS"]]
    if (is.null(o) || is.null(e)) return(NULL)
    ids <- intersect(o$feature_id[is.finite(o$pvalue)], e$feature_id[is.finite(e$pvalue)])
    mc <- S4Vectors::mcols(comparison_cache[[cmp]]$gene_dds)
    bm <- stats::setNames(mc$baseMean, rownames(comparison_cache[[cmp]]$gene_dds))
    d <- data.frame(feature_id = ids, gene = unname(sym[ids]), pas_mean = unname(bm[ids]),
                    p_o = o$pvalue[match(ids, o$feature_id)], s_o = o$Standard[match(ids, o$feature_id)],
                    p_e = e$pvalue[match(ids, e$feature_id)], s_e = e$Standard[match(ids, e$feature_id)])
    d <- d[valid_gene_symbol(d$gene), , drop = FALSE]
    d$gene_key <- gene_key(d$gene)
    list(cmp = cmp, d = d)
  })
  ab <- Filter(Negate(is.null), ab)
  apa <- if (is.null(APA_RESULTS) || is.null(APA_RESULTS$gene) || !nrow(APA_RESULTS$gene)) NULL else APA_RESULTS$gene
  rat <- unique(c(unlist(lapply(ab, function(x) x$d$gene)), if (!is.null(apa)) sub("\\|.*$", "", apa$gene_unit)))
  rat <- rat[valid_gene_symbol(rat)]
  bg <- as.data.frame(babelgene::orthologs(genes = rat, species = TWAS_TARGET_SPECIES, human = FALSE,
                                           min_support = TWAS_ORTHOLOG_MIN_SUPPORT, top = FALSE))
  assert_required_columns(bg, c("symbol", "human_symbol"), "rat background orthology")
  bgkey <- unique(gene_key(bg$symbol)); hit <- unique(mapping$gene_key)
  sup <- if ("support_n" %in% names(bg)) tapply(suppressWarnings(as.numeric(bg$support_n)), gene_key(bg$symbol), max) else NULL
  n_orth <- tapply(bg$human_symbol, gene_key(bg$symbol), function(v) length(unique(v)))
  up <- unique(bg[, c("symbol", "human_symbol")]); hc <- table(up$human_symbol)
  one_keys <- unique(gene_key(up$symbol[up$human_symbol %in% names(hc[hc == 1])]))
  one_keys <- one_keys[as.numeric(n_orth[one_keys]) == 1]
  finish <- function(g) {
    g <- g[g$gene_key %in% bgkey, , drop = FALSE]
    g$TWAS <- g$gene_key %in% hit
    g$support <- if (is.null(sup)) 1 else as.numeric(sup[g$gene_key]); g$support[!is.finite(g$support)] <- 1
    g
  }
  rows <- list()
  add <- function(g, cmp, endpoint, model, scope, genes_label) {
    base <- data.frame(Comparison = cmp, Endpoint = endpoint, Model = model, Scope = scope, Genes = nrow(g),
                       TWAS_genes = sum(g$TWAS), Positive_genes = sum(g$Y), Overlap = sum(g$Y & g$TWAS),
                       Overlap_genes = paste(sort(unique(genes_label[g$gene_key[g$Y & g$TWAS]])), collapse = ";"))
    rows[[length(rows) + 1L]] <<- twas_logistic(g, base)
    if (scope == "all_supported") rows[[length(rows) + 1L]] <<- twas_camera(g, base)
  }
  for (x in ab) {
    d <- x$d
    for (m in c("Original", "EVS")) {
      pc <- if (m == "Original") "p_o" else "p_e"; sc <- if (m == "Original") "s_o" else "s_e"
      g <- d |> dplyr::group_by(gene_key) |>
        dplyr::summarise(gene = dplyr::first(gene), Y = any(.data[[sc]] %in% TRUE), nPAS = dplyr::n(),
                         mean_count = mean(pas_mean), gene_p = simes_p(.data[[pc]]), .groups = "drop")
      g <- finish(as.data.frame(g)); lab <- stats::setNames(g$gene, g$gene_key)
      add(g, x$cmp, "Secondary: PAS abundance", m, "all_supported", lab)
      add(g[g$gene_key %in% one_keys, , drop = FALSE], x$cmp, "Secondary: PAS abundance", m, "one_to_one", lab)
    }
  }
  if (!is.null(apa)) for (cmp in unique(apa$Comparison)) {
    z <- apa[apa$Comparison == cmp, , drop = FALSE]
    cache <- comparison_cache[[cmp]]
    nm <- rowMeans(sweep(cache$raw, 2, cache$sf[colnames(cache$raw)], "/"))
    pas <- APA_RESULTS$pas[APA_RESULTS$pas$Comparison == cmp, , drop = FALSE]
    unit_mean <- tapply(nm[as.character(pas$PAS)], pas$gene_unit, mean)
    unit_n <- table(pas$gene_unit)
    z$gene <- sub("\\|.*$", "", z$gene_unit); z <- z[valid_gene_symbol(z$gene), , drop = FALSE]
    z$gene_key <- gene_key(z$gene)
    z$pas_mean <- as.numeric(unit_mean[z$gene_unit]); z$n <- as.numeric(unit_n[z$gene_unit])
    g <- z |> dplyr::group_by(gene_key) |>
      dplyr::summarise(gene = dplyr::first(gene), Y = any(FDR < APA_FDR, na.rm = TRUE), nPAS = sum(n),
                       mean_count = mean(pas_mean, na.rm = TRUE), gene_p = simes_p(P.Value), .groups = "drop")
    g <- finish(as.data.frame(g)); lab <- stats::setNames(g$gene, g$gene_key)
    add(g, cmp, "Primary: differential APA (usage)", "APA (model-independent)", "all_supported", lab)
    add(g[g$gene_key %in% one_keys, , drop = FALSE], cmp, "Primary: differential APA (usage)",
        "APA (model-independent)", "one_to_one", lab)
  }
  d <- dplyr::bind_rows(rows)
  if (nrow(d)) {
    d$BH <- NA_real_
    for (mm in unique(d$Method)) { w <- which(d$Method %in% mm & is.finite(d$P)); if (length(w)) d$BH[w] <- stats::p.adjust(d$P[w], "BH") }
    d <- d[order(d$Endpoint != "Primary: differential APA (usage)", d$Comparison, d$Model, d$Scope, d$Method), , drop = FALSE]
  }
  if (nrow(d)) for (cmp in unique(d$Comparison)) for (sc in c("all_supported", "one_to_one")) {
    w <- d$Comparison == cmp & d$Scope == sc & d$Endpoint == "Secondary: PAS abundance" & !grepl("cameraPR", d$Method)
    o <- d[w & d$Model == "Original", , drop = FALSE]; e <- d[w & d$Model == "EVS", , drop = FALSE]
    if (nrow(o) && nrow(e)) add_validation("3'aTWAS: Original and EVS use identical gene backgrounds", paste(cmp, sc),
                                           identical(o$Genes, e$Genes) && identical(o$TWAS_genes, e$TWAS_genes),
                                           sprintf("genes %d vs %d", o$Genes[1], e$Genes[1]))
  }
  save_csv(d, file.path(summary_table_dir, "Table_TWAS_Enrichment.csv"))
  if (!is.null(APA_RESULTS)) APA_RESULTS$enrichment <<- d[d$Endpoint == "Primary: differential APA (usage)" &
                                                           d$Scope == "all_supported" & !grepl("cameraPR", d$Method), , drop = FALSE]
  invisible(d)
}

# =============================================================================
# VALIDATION: ERROR RATES, POWER, STABILITY AND FILTER SENSITIVITY
# =============================================================================
# Purpose. EVS chooses its partition from the counts it tests, so whether its
# p-values are valid and whether it gains power cannot be settled by theory;
# both are measured here by re-running the complete adaptive analysis (filter,
# normalization, VST, PC1 ranking, likelihood, k* search, pooled-vs-split
# decision, moderation, testing, BH) on data whose truth is known or controlled.
# Every replicate builds its own cache from its own counts; no observed-data
# ranking, cutoff or dispersion is reused (checked at the end of the stage).
#  * Label permutation of the observed data is not used: permuting animals
#    dilutes but does not remove real treatment effects, so it gives no true
#    null. Error rates come from the simulations (known truth) and from the
#    real-data null splits of the external liver data.
#  * Negative binomial simulation with known truth, from the fitted parameters
#    of one comparison (abundances, size factors, mean-dispersion trends and
#    prior variances, gene-level sample effects for PAS dependence). Scenarios
#    include a complete null, mixtures of null and non-null PASs with weak
#    (|log2FC| 0.5), moderate (1) and strong (2) effects, homogeneous, two-
#    population and continuous dispersion heterogeneity, and PC1 unrelated to
#    dispersion (seven scenarios). Only two scenarios are generated from the
#    EVS populations; the others do not favour EVS by construction.
#  * Cutoff stability by resampling animals within arms (default: leave one
#    animal out per arm; bootstrap optional).
#  * Sensitivity to the abundance filter.
# Endpoints (Monte Carlo mean and 95% interval mean +/- 1.96 SE): empirical FDR
# (mean false discovery proportion, zero for replicates without discoveries),
# probability of any false discovery, power (true-positive rate) overall and by
# effect size, false-positive rate, specificity, precision, false
# non-discovery rate, dispersion bias and RMSE (log scale), log2FC bias and
# RMSE, Wald 95% interval coverage, runtime, memory and failures. Forced EVS
# and the complete adaptive pooled-vs-EVS procedure are reported separately.
# Settings. SEQUENCE_PROFILE=development (default) or publication sets the
# replicate numbers; each can be overridden (SEQUENCE_PERMUTATIONS,
# SEQUENCE_SIM_REPS, SEQUENCE_SIM_NPAS, SEQUENCE_STABILITY_SUBSETS,
# SEQUENCE_FILTER_THRESHOLDS, SEQUENCE_SIM_SCENARIOS, SEQUENCE_SIM_MODELS).
# Development settings reveal gross departures only; claims of FDR control and
# power require the publication settings.
# Adaptive procedure. Every replicate (simulation, stability subset, external
# analysis) generates its own pooled-vs-split bootstrap from
# its own counts, normalization, ranking, likelihood and cutoff search, and
# applies the same rule (selection_decision) as the observed comparisons. No
# bootstrap distribution is reused between analyses.
# FDR assessment. The empirical FDR (mean false discovery proportion, a
# bounded [0, 1] variable) gets an empirical-Bernstein 95% interval (Maurer &
# Pontil 2009), which has non-zero width even when every replicate has zero
# false discoveries. With the prespecified tolerance (SEQUENCE_FDR_TOLERANCE,
# default 0.025) the threshold is 0.10 + 0.025: Acceptable if the upper bound
# is at most the threshold, Unacceptable if the lower bound exceeds it,
# Inconclusive otherwise.
VALIDATION_ON <- !identical(Sys.getenv("SEQUENCE_VALIDATION", "1"), "0")
SIM_REPS <- env_int("SEQUENCE_SIM_REPS", if (PUB) 1000L else 2L)
SIM_NPAS <- env_int("SEQUENCE_SIM_NPAS", if (PUB) 10000L else 3000L)
SIM_BASE <- Sys.getenv("SEQUENCE_SIM_BASE", "RT4_ZT10")
SIM_MODELS <- trimws(strsplit(Sys.getenv("SEQUENCE_SIM_MODELS",
  paste0("deseq2_default,conventional,edger,pooled,pooled_tmm,pc1,variance,random,adaptive",
         if (PUB) ",pc12,continuous" else "")), ",")[[1]])
ALL_SCENARIOS <- c("null_homogeneous", "effects_homogeneous", "null_two_populations", "effects_two_populations",
                   "effects_continuous", "null_pc1_unrelated", "effects_pc1_unrelated")
SIM_SCENARIOS <- {
  v <- trimws(strsplit(Sys.getenv("SEQUENCE_SIM_SCENARIOS", paste(ALL_SCENARIOS, collapse = ",")), ",")[[1]])
  v[v %in% ALL_SCENARIOS]
}
SIM_EFFECTS <- c(weak = 0.5, moderate = 1, strong = 2)
SIM_RETRIES <- 2L                # re-analyses of the same simulated matrix after a failure
FDR_TOLERANCE <- suppressWarnings(as.numeric(Sys.getenv("SEQUENCE_FDR_TOLERANCE", "0.025")))
if (!is.finite(FDR_TOLERANCE) || FDR_TOLERANCE < 0) FDR_TOLERANCE <- 0.025
SIM_EFFECT_FRACTION <- 0.02      # of PASs per effect class (6% non-null in effect scenarios)
SIM_GENE_SD <- 0.10              # sample-level effect shared by the PASs of a gene (PAS dependence)
SIM_CONT_SD <- 0.8               # continuous heterogeneity: SD of a PAS-level shift of log dispersion
SIM_FACTOR_SD <- 0.60            # loadings of the latent factor that dominates PC1 (pc1_unrelated)
SIM_FACTOR_FRACTION <- 0.10
STAB_N <- env_int("SEQUENCE_STABILITY_SUBSETS", if (PUB) 25L else 3L)
STAB_METHOD <- match.arg(Sys.getenv("SEQUENCE_STABILITY_METHOD", "subsample"), c("subsample", "bootstrap"))
FILTER_THRESHOLDS <- trimws(strsplit(Sys.getenv("SEQUENCE_FILTER_THRESHOLDS", if (PUB) "1x1,5x3,20x5" else "1x1,20x5"), ",")[[1]])
MODEL_LABELS <- c(deseq2_default = "DESeq2 default", conventional = "Conventional DESeq2", edger = "edgeR robust QL",
                  pooled = "Pooled hierarchical", pooled_tmm = "Pooled hierarchical, TMM", pc1 = "PC1 EVS",
                  variance = "Variance split", pc12 = "PC1-PC2 split", continuous = "Continuous rank",
                  random = "Random partition", adaptive = "Adaptive (pooled vs EVS)")

# Complete re-analysis of one count matrix with the primary rules.
# filter: c(min_count, min_samples) (default: the design's rule). When the
# adaptive model is requested, a new pooled-vs-split bootstrap is generated
# from this matrix (its own normalization, ranking, likelihood and cutoff
# search, with the same filter) and the shared rule selection_decision() is
# applied. Returns per-model calls, k*, split fits, the PC1 depths, the coarse
# LRT, the bootstrap statistics and selection, runtimes, failures and memory.
reanalyze_matrix <- function(y, coldata, label, models, annot_df, seed, filter = NULL) {
  set.seed(seed)
  if (is.null(filter)) filter <- filter_rule(coldata)
  invisible(gc(reset = TRUE))
  t0 <- proc.time()[["elapsed"]]
  y <- y[eligible_pas(y, filter[1], filter[2]), , drop = FALSE]
  y <- coerce_raw_count_matrix_for_deseq2(y, context = label)
  cache <- build_comparison_cache(y, coldata)
  y <- cache$raw; ids <- rownames(y)
  si <- split_inputs(cache)
  one <- fit_marginal_trend(si$lx, si$lk, multistart = TRUE)
  times <- list(); failed <- character(); out <- list(); ext <- list(); ks <- list(); fits <- list(); lrt <- NA_real_
  timed <- function(nm, expr) {
    st <- proc.time()[["elapsed"]]
    v <- tryCatch(expr, error = function(e) { failed <<- c(failed, nm); NULL })
    times[[nm]] <<- proc.time()[["elapsed"]] - st
    if (is.null(v) && !nm %in% failed) failed <<- c(failed, nm)
    v
  }
  o <- fit_view_model(y, coldata, paste(label, "pooled"), size_factors = cache$sf, gene_dds = cache$gene_dds)
  cooks <- as.character(o$res_df$feature_id[o$res_df$cooks_outlier %in% TRUE])
  so <- call_significance(list(raw_dataset = o), paste(label, "pooled"), annot_df, lean = TRUE)
  hc <- so$null$hc_p_threshold[1]; if (is.null(hc) || !is.finite(hc)) hc <- NA_real_
  out$pooled <- list(raw_dataset = so$views$raw_dataset$results)
  times$pooled <- proc.time()[["elapsed"]] - t0
  if ("conventional" %in% models)
    out$conventional <- timed("conventional", lean_analysis(list(raw_dataset = y), coldata, paste(label, "conventional"),
                                                            cache$sf, cache$gene_dds, cooks, hc, annot_df, trend = "deseq2_default"))
  if ("pooled_tmm" %in% models) out$pooled_tmm <- timed("pooled_tmm", {
    sft <- tmm_size_factors(y)
    lean_analysis(list(raw_dataset = y), coldata, paste(label, "pooled TMM"), sft, comparison_gene_dds(y, coldata, sft),
                  cooks, hc, annot_df) })
  if ("deseq2_default" %in% models) ext$deseq2_default <- timed("deseq2_default", fit_deseq2_default(y, coldata))
  if ("edger" %in% models) ext$edger <- timed("edger", fit_edger_ql(y, coldata))
  split_model <- function(dp, nm) {
    s <- tryCatch(select_cutoff(si, dp, one), error = function(e) NULL)
    if (is.null(s)) return(NULL)
    if (nm == "pc1") lrt <<- 2 * (s$coarse_max_loglik - one$loglik)
    mats <- split_matrices(y, dp, s$k)
    if (is.null(mats)) return(NULL)
    ks[[nm]] <<- s$k
    fits[[nm]] <<- list(le_tau2 = s$best$le$tau2, rem_tau2 = s$best$rem$tau2, le_par = s$best$le$par,
                        rem_par = s$best$rem$par, le_ids = rownames(mats$leading_edge_dataset))
    lean_analysis(mats, coldata, paste(label, nm), cache$sf, cache$gene_dds, cooks, hc, annot_df)
  }
  need_pc1 <- any(c("pc1", "adaptive") %in% models)
  for (m in c(if (need_pc1) "pc1", intersect(models, c("variance", "pc12"))))
    out[[m]] <- timed(m, split_model(cache$depths[[m]], m))
  if ("random" %in% models)
    out$random <- timed("random", split_model(stats::setNames(sample(unname(cache$depths$pc1)), ids), "random"))
  if ("continuous" %in% models)
    out$continuous <- timed("continuous", lean_analysis(list(raw_dataset = y), coldata, paste(label, "continuous"), cache$sf,
                                                        cache$gene_dds, cooks, hc, annot_df,
                                                        basis = continuous_basis(cache$depths$pc1, ids)))
  out <- Filter(Negate(is.null), out)
  calls <- lapply(out, function(res) {
    f <- discovery_flags(res)
    est <- dplyr::bind_rows(lapply(res, function(r) r[, c("feature_id", "log2FoldChange", "lfcSE", "dispersion_final"), drop = FALSE]))
    i <- match(f$feature_id, est$feature_id)
    f$log2FoldChange <- est$log2FoldChange[i]; f$lfcSE <- est$lfcSE[i]; f$dispersion_final <- est$dispersion_final[i]
    f
  })
  calls <- c(calls, Filter(Negate(is.null), ext))
  selected <- NA_character_; sel_p <- NA_real_; boot <- numeric()
  if ("adaptive" %in% models && !is.null(calls$pc1)) {
    boot <- timed("selection bootstrap", if (SELECTION_BOOT > 0L) selection_bootstrap(cache, one, seed + 7919L, filter) else numeric())
    if (is.null(boot)) boot <- numeric()
    d <- selection_decision(lrt, boot); selected <- d$model; sel_p <- d$p
    if (selected %in% c("EVS", "Pooled")) calls$adaptive <- if (selected == "EVS") calls$pc1 else calls$pooled
  }
  if (!"pooled" %in% models) calls$pooled <- NULL
  if (!"pc1" %in% models) calls$pc1 <- NULL
  g <- gc()
  list(calls = calls, k = ks, fits = fits, depth = cache$depths$pc1, lrt = lrt, selected = selected, selection_p = sel_p,
       boot_n = length(boot), boot_ok = sum(is.finite(boot)),
       times = times, failed = failed, memory_mb = sum(g[, ncol(g)]),
       total_seconds = proc.time()[["elapsed"]] - t0, n_pas = nrow(y))
}

# Null p-value summaries of one model's calls (type I error, uniformity).
null_metrics <- function(f) {
  p <- f$pvalue[is.finite(f$pvalue)]
  data.frame(PAS_tested = nrow(f), Standard = sum(f$Standard %in% TRUE), Any_Standard = any(f$Standard %in% TRUE),
             TypeI_p05 = mean(p < 0.05), TypeI_p01 = mean(p < 0.01),
             KS_D = if (length(p) > 1L) unname(suppressWarnings(stats::ks.test(p, "punif"))$statistic) else NA_real_)
}

# ---- Simulations --------------------------------------------------------------
simulation_base <- function() {
  if (!SIM_BASE %in% names(comparison_cache)) stop("Unknown simulation base comparison: ", SIM_BASE)
  cache <- comparison_cache[[SIM_BASE]]; fit <- CUTOFF$fits[[SIM_BASE]]
  mc <- S4Vectors::mcols(cache$gene_dds)
  ids <- rownames(cache$gene_dds)[is.finite(mc$baseMean) & mc$baseMean > 0]
  ann <- pas_annotation()
  unit <- ann$gene_unit[match(ids, as.character(ann$feature_id))]; unit[is.na(unit)] <- ids[is.na(unit)]
  k <- OPTIMAL_K[[SIM_BASE]]
  tr <- function(p) { p <- unname(p); list(f = make_trend_function(exp(p[1]), exp(p[2]), p[3]), tau2 = exp(p[4])) }
  list(ids = ids, mean = stats::setNames(mc$baseMean[match(ids, rownames(cache$gene_dds))], ids),
       unit = stats::setNames(unit, ids), in_le = stats::setNames(cache$depth[ids] <= k, ids),
       sf = cache$sf, coldata = cache$coldata, one = tr(fit$one), le = tr(fit$le), rem = tr(fit$rem),
)
}

simulate_counts <- function(base, scenario, seed) {
  set.seed(seed)
  n <- min(SIM_NPAS, length(base$ids))
  ids <- sample(base$ids, n)
  mu <- base$mean[ids]; unit <- base$unit[ids]
  # Design and depth of the base comparison (five animals per group).
  sf <- base$sf; coldata <- base$coldata
  smp <- names(sf)
  la <- if (grepl("two_populations", scenario)) {
    le <- base$in_le[ids]
    ifelse(le, log(base$le$f(mu)) + stats::rnorm(n, 0, sqrt(base$le$tau2)),
           log(base$rem$f(mu)) + stats::rnorm(n, 0, sqrt(base$rem$tau2)))
  } else if (grepl("pc1_unrelated", scenario)) {
    le <- sample(base$in_le[ids])                  # two populations assigned at random to PASs
    ifelse(le, log(base$le$f(mu)) + stats::rnorm(n, 0, sqrt(base$le$tau2)),
           log(base$rem$f(mu)) + stats::rnorm(n, 0, sqrt(base$rem$tau2)))
  } else if (grepl("continuous", scenario)) {
    log(base$one$f(mu)) + SIM_CONT_SD * stats::rnorm(n) + stats::rnorm(n, 0, sqrt(base$one$tau2))
  } else log(base$one$f(mu)) + stats::rnorm(n, 0, sqrt(base$one$tau2))
  alpha <- pmin(pmax(exp(la), 1e-8), 1e3)
  beta <- rep(0, n); cls <- rep("null", n)
  if (grepl("^effects", scenario)) {
    m <- round(n * SIM_EFFECT_FRACTION); idx <- sample.int(n, 3L * m)
    for (j in seq_along(SIM_EFFECTS)) {
      w <- idx[((j - 1L) * m + 1L):(j * m)]
      beta[w] <- SIM_EFFECTS[[j]] * sample(c(-1, 1), m, replace = TRUE); cls[w] <- names(SIM_EFFECTS)[j]
    }
  }
  gu <- unique(unit)
  G <- matrix(stats::rnorm(length(gu) * length(smp), 0, SIM_GENE_SD), length(gu))[match(unit, gu), , drop = FALSE]
  if (grepl("pc1_unrelated", scenario)) {
    load <- ifelse(stats::runif(n) < SIM_FACTOR_FRACTION, stats::rnorm(n, 0, SIM_FACTOR_SD), 0)
    G <- G + outer(load, stats::rnorm(length(smp)))
  }
  trt <- as.character(coldata[smp, "condition"]) == "trt"
  M <- outer(mu, sf[smp]) * exp(G) * 2^outer(beta, as.numeric(trt))
  y <- matrix(stats::rnbinom(length(M), mu = as.vector(M), size = rep(1 / alpha, times = ncol(M))),
              nrow = n, dimnames = list(ids, smp))
  list(y = y, coldata = coldata, beta = stats::setNames(beta, ids), class = stats::setNames(cls, ids),
       alpha = stats::setNames(alpha, ids))
}

simulation_metrics <- function(f, sim) {
  b <- sim$beta[f$feature_id]; cl <- sim$class[f$feature_id]; at <- sim$alpha[f$feature_id]
  call <- f$Standard %in% TRUE
  nz <- b != 0; small <- b != 0 & abs(b) < lfc_boundary; big <- abs(b) > lfc_boundary
  fdp <- function(cl_, tr) { if (all(is.na(cl_))) return(NA_real_); cl_ <- cl_ %in% TRUE
    if (sum(cl_)) sum(cl_ & !tr) / sum(cl_) else 0 }
  tpr <- function(cl_, tr) { if (all(is.na(cl_))) return(NA_real_); cl_ <- cl_ %in% TRUE
    if (sum(tr)) sum(cl_ & tr) / sum(tr) else NA_real_ }
  est <- f$log2FoldChange; se <- f$lfcSE; ok <- is.finite(est)
  okse <- ok & is.finite(se)
  ad <- f$dispersion_final; okd <- is.finite(ad) & ad > 0
  out <- data.frame(PAS = nrow(f), Nonzero = sum(nz), Discoveries = sum(call),
    FDP = fdp(call, nz), Any_false = sum(call & !nz) > 0, TPR = tpr(call, nz),
    FPR = if (sum(!nz)) sum(call & !nz) / sum(!nz) else NA_real_,
    Specificity = if (sum(!nz)) 1 - sum(call & !nz) / sum(!nz) else NA_real_,
    Precision = if (sum(call)) sum(call & nz) / sum(call) else NA_real_,
    FNDR = if (sum(!call)) sum(!call & nz) / sum(!call) else NA_real_,
    FDP_Strong = fdp(f$Strong, big), TPR_Strong = tpr(f$Strong, big),
    FDP_Weak = fdp(f$Weak, small), TPR_Weak = tpr(f$Weak, small),
    FDP_HBFSS = if (is.null(f$HBFSS)) NA_real_ else fdp(f$HBFSS, nz), TPR_HBFSS = if (is.null(f$HBFSS)) NA_real_ else tpr(f$HBFSS, nz),
    Dispersion_log_bias = if (any(okd)) mean(log(ad[okd]) - log(at[okd])) else NA_real_,
    Dispersion_log_RMSE = if (any(okd)) sqrt(mean((log(ad[okd]) - log(at[okd]))^2)) else NA_real_,
    LFC_bias_nonzero = if (any(nz & ok)) mean(est[nz & ok] - b[nz & ok]) else NA_real_,
    LFC_RMSE = if (any(ok)) sqrt(mean((est[ok] - b[ok])^2)) else NA_real_,
    LFC_RMSE_nonzero = if (any(nz & ok)) sqrt(mean((est[nz & ok] - b[nz & ok])^2)) else NA_real_,
    CI95_coverage = if (any(okse)) mean(abs(est[okse] - b[okse]) <= stats::qnorm(0.975) * se[okse]) else NA_real_)
  for (ec in names(SIM_EFFECTS)) out[[paste0("TPR_", ec)]] <- tpr(call, cl == ec)
  out
}

# Empirical-Bernstein 95% interval for the mean of a [0, 1] variable (Maurer &
# Pontil 2009, Theorem 4, applied to each side at level 0.025).
eb_interval <- function(x, level = 0.95) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 2L) return(c(0, 1))
  l <- log(4 / (1 - level)); m <- mean(x); v <- stats::var(x)
  h <- sqrt(2 * v * l / n) + 7 * l / (3 * (n - 1))
  c(max(0, m - h), min(1, m + h))
}

mc_summary <- function(d, groups) {
  se <- function(x) { x <- x[is.finite(x)]; if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_ }
  mn <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
  vars <- c("FDP", "TPR", "FPR", "Specificity", "Precision", "FNDR", "FDP_Strong", "TPR_Strong", "FDP_Weak", "TPR_Weak", "FDP_HBFSS", "TPR_HBFSS",
            "Dispersion_log_bias", "Dispersion_log_RMSE", "LFC_bias_nonzero", "LFC_RMSE", "LFC_RMSE_nonzero",
            "CI95_coverage", "TPR_weak", "TPR_moderate", "TPR_strong", "Discoveries", "Seconds")
  out <- lapply(split(d, d[groups], drop = TRUE), function(z) {
    r <- z[1, groups, drop = FALSE]; r$Attempted <- nrow(z); r$Failed <- sum(z$Failed)
    z <- z[!z$Failed, , drop = FALSE]; r$Successful <- nrow(z)
    r$Prob_any_false_discovery <- mean(z$Any_false %in% TRUE)
    for (v in intersect(vars, names(z))) {
      m <- mn(z[[v]]); e <- se(z[[v]])
      r[[v]] <- m; r[[paste0(v, "_lo95")]] <- m - 1.96 * e; r[[paste0(v, "_hi95")]] <- m + 1.96 * e
    }
    if (nrow(z) && "FDP" %in% names(z)) { b <- eb_interval(z$FDP); r$FDP_lo95 <- b[1]; r$FDP_hi95 <- b[2] }
    r
  })
  out <- dplyr::bind_rows(out)
  out$Nominal_FDR <- BH_FDR_STANDARD
  out$Empirical_FDR <- out$FDP
  out$FDR_interval <- "empirical Bernstein 95%"
  out$FDR_tolerance <- FDR_TOLERANCE
  thr <- BH_FDR_STANDARD + FDR_TOLERANCE
  out$FDR_assessment <- ifelse(!is.finite(out$FDP_hi95) | !is.finite(out$FDP_lo95), "Inconclusive",
                               ifelse(out$FDP_hi95 <= thr, "Acceptable",
                                      ifelse(out$FDP_lo95 > thr, "Unacceptable", "Inconclusive")))
  out
}

run_simulation_validation <- function(annot_df) {
  base <- simulation_base(); rows <- list(); perf <- list()
  for (si in seq_along(SIM_SCENARIOS)) for (r in seq_len(SIM_REPS)) {
    scn <- SIM_SCENARIOS[si]; seed <- SPLIT_SEED + 10000L * si + r
    sim <- simulate_counts(base, scn, seed)
    # Bounded retry on the same simulated matrix (a failure is never counted as zero discoveries).
    res <- NULL; attempt <- 0L
    while (is.null(res) && attempt <= SIM_RETRIES) {
      res <- tryCatch(reanalyze_matrix(sim$y, sim$coldata, sprintf("simulation %s %d", scn, r), SIM_MODELS, annot_df,
                                       seed + 100000L * attempt), error = function(e) NULL)
      attempt <- attempt + 1L
    }
    if (is.null(res)) {
      for (m in SIM_MODELS) rows[[length(rows) + 1L]] <- data.frame(Scenario = scn, Replicate = r, Model = MODEL_LABELS[[m]], Failed = TRUE)
      next
    }
    for (m in setdiff(SIM_MODELS, names(res$calls)))
      rows[[length(rows) + 1L]] <- data.frame(Scenario = scn, Replicate = r, Model = MODEL_LABELS[[m]], Failed = TRUE)
    for (m in names(res$calls)) rows[[length(rows) + 1L]] <- cbind(
      data.frame(Scenario = scn, Replicate = r, Model = MODEL_LABELS[[m]], Failed = FALSE,
                 k = if (is.null(res$k[[m]])) NA_integer_ else res$k[[m]],
                 Selected = if (m == "adaptive") res$selected else NA_character_,
                 Seconds = if (is.null(res$times[[m]])) NA_real_ else res$times[[m]]),
      simulation_metrics(res$calls[[m]], sim))
    perf[[length(perf) + 1L]] <- data.frame(Scenario = scn, Replicate = r, PAS = res$n_pas, Total_seconds = res$total_seconds,
                                            Max_memory_MB = res$memory_mb, Failed_models = paste(res$failed, collapse = ";"),
                                            PC1_coarse_LRT = res$lrt, Adaptive_selected = res$selected,
                                            Selection_p = res$selection_p, Bootstrap_successful = res$boot_ok,
                                            Attempts = attempt)
  }
  d <- dplyr::bind_rows(rows)
  if (!nrow(d)) return(invisible(NULL))
  # Aggregate results only (no per-replicate files).
  save_csv(mc_summary(d, c("Scenario", "Model")), file.path(summary_table_dir, "Table_Validation_Simulation_Summary.csv"))
  pf <- dplyr::bind_rows(perf)
  save_csv(as.data.frame(pf |> dplyr::group_by(Scenario) |>
    dplyr::summarise(Replicates_analysed = dplyr::n(), PAS_median = stats::median(PAS),
                     Seconds_mean = mean(Total_seconds), Seconds_max = max(Total_seconds),
                     Max_memory_MB = max(Max_memory_MB), Mean_attempts = mean(Attempts),
                     Replicates_with_failed_models = sum(nzchar(Failed_models)), .groups = "drop")),
    file.path(summary_table_dir, "Table_Validation_Simulation_Performance.csv"))
  sel <- pf |> dplyr::group_by(Scenario) |>
    dplyr::summarise(Replicates = dplyr::n(), EVS_selected_fraction = mean(Adaptive_selected == "EVS", na.rm = TRUE),
                     Undetermined_fraction = mean(Adaptive_selected == "Undetermined", na.rm = TRUE),
                     Median_PC1_coarse_LRT = stats::median(PC1_coarse_LRT, na.rm = TRUE), .groups = "drop")
  save_csv(as.data.frame(sel), file.path(summary_table_dir, "Table_Validation_Model_Selection.csv"))
  invisible(d)
}

# ---- Cutoff stability ---------------------------------------------------------
stability_sets <- function(trt, ctl) {
  if (STAB_METHOD == "subsample") {
    g <- expand.grid(a = trt, b = ctl, stringsAsFactors = FALSE)
    lapply(seq_len(nrow(g)), function(r) setdiff(c(trt, ctl), c(g$a[r], g$b[r])))
  } else lapply(seq_len(max(STAB_N, 1L)), function(r) c(sample(trt, replace = TRUE), sample(ctl, replace = TRUE)))
}

run_stability_analysis <- function(annot_df) {
  rows <- list()
  for (i in seq_len(nrow(comparison_table))) {
    cmp <- comparison_table$comparison_name[i]; cache <- comparison_cache[[cmp]]
    cd <- cache$coldata; raw <- comparison_raw_counts(cmp); smp <- colnames(raw)
    trt <- smp[cd[smp, "condition"] == "trt"]; ctl <- smp[cd[smp, "condition"] == "untrt"]
    set.seed(SPLIT_SEED + 500L * i)
    sets <- stability_sets(trt, ctl)
    use <- if (STAB_N >= length(sets)) seq_along(sets) else sort(sample.int(length(sets), STAB_N))
    full_le <- names(cache$depth)[cache$depth <= OPTIMAL_K[[cmp]]]
    full_std <- MODEL_FLAGS[[cmp]][["PC1 EVS"]]; full_std <- if (is.null(full_std)) character() else full_std$feature_id[full_std$Standard]
    cs <- CUTOFF$summary[CUTOFF$summary$Comparison == cmp, , drop = FALSE]
    for (j in use) {
      s <- sets[[j]]
      y <- raw[, s, drop = FALSE]
      if (anyDuplicated(s)) colnames(y) <- make.unique(s, sep = "_rep")
      cdj <- data.frame(condition = factor(as.character(cd[s, "condition"]), levels = levels(cd$condition)),
                        row.names = colnames(y))
      r <- tryCatch(reanalyze_matrix(y, cdj, sprintf("%s subset %d", cmp, j), c("pooled", "pc1", "adaptive"), annot_df,
                                     SPLIT_SEED + 2000L * i + j),
                    error = function(e) { warning(cmp, " stability subset ", j, " failed: ", conditionMessage(e), call. = FALSE); NULL })
      if (is.null(r) || is.null(r$k$pc1)) next
      le <- r$fits$pc1$le_ids
      common <- intersect(names(r$depth), names(cache$depth))
      std <- r$calls$pc1$feature_id[r$calls$pc1$Standard]
      p <- unname(r$fits$pc1$le_par); q <- unname(r$fits$pc1$rem_par)
      rows[[length(rows) + 1L]] <- data.frame(Comparison = cmp, Method = STAB_METHOD, Subset = j,
        Samples = paste(s, collapse = ","), k = r$k$pc1, k_full = OPTIMAL_K[[cmp]], Leading_Edge = length(le),
        Leading_Edge_full = length(full_le),
        LE_Jaccard_with_full = length(intersect(le, full_le)) / length(union(le, full_le)),
        Depth_Spearman_with_full = suppressWarnings(stats::cor(r$depth[common], cache$depth[common], method = "spearman")),
        LE_prior_var = r$fits$pc1$le_tau2, REM_prior_var = r$fits$pc1$rem_tau2,
        LE_prior_var_full = cs$LE_prior_var, REM_prior_var_full = cs$REM_prior_var,
        LE_a0 = exp(p[1]), LE_a1 = exp(p[2]), LE_b = p[3], REM_a0 = exp(q[1]), REM_a1 = exp(q[2]), REM_b = q[3],
        Selected = r$selected, Standard = length(std), Standard_full = length(full_std),
        Full_discoveries_recovered = if (length(full_std)) length(intersect(std, full_std)) / length(full_std) else NA_real_,
        Discovery_Jaccard_with_full = if (length(union(std, full_std))) length(intersect(std, full_std)) / length(union(std, full_std)) else NA_real_)
      save_csv(dplyr::bind_rows(rows), file.path(summary_table_dir, "Table_Cutoff_Stability.csv"))
      message(sprintf("[validation] %s stability subset %d done", cmp, j))
    }
  }
  d <- dplyr::bind_rows(rows)
  if (!nrow(d)) return(invisible(NULL))
  s <- d |> dplyr::group_by(Comparison, Method) |>
    dplyr::summarise(Subsets = dplyr::n(), k_full = dplyr::first(k_full), k_median = stats::median(k),
                     k_q025 = stats::quantile(k, 0.025), k_q975 = stats::quantile(k, 0.975),
                     LE_size_median = stats::median(Leading_Edge),
                     LE_Jaccard_mean = mean(LE_Jaccard_with_full), Depth_Spearman_mean = mean(Depth_Spearman_with_full),
                     LE_prior_var_median = stats::median(LE_prior_var), LE_prior_var_IQR = stats::IQR(LE_prior_var),
                     REM_prior_var_median = stats::median(REM_prior_var), REM_prior_var_IQR = stats::IQR(REM_prior_var),
                     LE_b_IQR = stats::IQR(LE_b), REM_b_IQR = stats::IQR(REM_b),
                     EVS_selected_fraction = mean(Selected == "EVS", na.rm = TRUE),
                     Full_discoveries_recovered_mean = mean(Full_discoveries_recovered, na.rm = TRUE),
                     Discovery_Jaccard_mean = mean(Discovery_Jaccard_with_full, na.rm = TRUE), .groups = "drop")
  save_csv(as.data.frame(s), file.path(summary_table_dir, "Table_Cutoff_Stability_Summary.csv"))
  invisible(d)
}

# ---- Filter sensitivity -------------------------------------------------------
run_filter_sensitivity <- function(annot_df) {
  rows <- list()
  for (thr in FILTER_THRESHOLDS) {
    v <- suppressWarnings(as.integer(strsplit(thr, "x", fixed = TRUE)[[1]]))
    if (length(v) != 2L || anyNA(v)) next
    for (i in seq_len(nrow(comparison_table))) {
      cmp <- comparison_table$comparison_name[i]; cache <- comparison_cache[[cmp]]
      r <- tryCatch(reanalyze_matrix(comparison_raw_counts(cmp), cache$coldata, sprintf("%s filter %s", cmp, thr),
                                     c("pooled", "pc1"), annot_df, SPLIT_SEED + 3000L * i, filter = v),
                    error = function(e) { warning(cmp, " filter ", thr, " failed: ", conditionMessage(e), call. = FALSE); NULL })
      if (is.null(r)) next
      rows[[length(rows) + 1L]] <- data.frame(Comparison = cmp, Filter = sprintf(">= %d reads in >= %d samples", v[1], v[2]),
        PAS_eligible = r$n_pas, k = r$k$pc1 %||% NA_integer_, Leading_Edge = length(r$fits$pc1$le_ids %||% character()),
        Standard_pooled = sum(r$calls$pooled$Standard), Standard_PC1_EVS = if (is.null(r$calls$pc1)) NA_integer_ else sum(r$calls$pc1$Standard))
      save_csv(dplyr::bind_rows(rows), file.path(summary_table_dir, "Table_Filter_Sensitivity.csv"))
      message(sprintf("[validation] %s filter %s done", cmp, thr))
    }
  }
  invisible(dplyr::bind_rows(rows))
}

# Unfiltered raw counts of a comparison's ten samples (validation re-applies the filter).
comparison_raw_counts <- function(cmp) {
  i <- match(cmp, comparison_table$comparison_name)
  smp <- rownames(comparison_cache[[cmp]]$coldata)
  x <- as.matrix(WTTS_Seq[, smp, drop = FALSE]); storage.mode(x) <- "numeric"
  rownames(x) <- rownames(WTTS_Seq); x[!is.finite(x)] <- 0
  x <- x[rowSums(x) > 0, , drop = FALSE]
  coerce_raw_count_matrix_for_deseq2(x, context = paste(cmp, "unfiltered counts"))
}

cache_fingerprint <- function() vapply(comparison_cache, function(c) paste(nrow(c$raw), signif(sum(c$sf), 10),
                                                                         sum(c$depth), sep = "|"), "")

run_validation <- function(annot_df) {
  status <- character(); before <- cache_fingerprint()
  stage <- function(nm, expr) {
    st <- proc.time()[["elapsed"]]
    v <- tryCatch(expr, error = function(e) { status <<- c(status, paste(nm, conditionMessage(e), sep = ": "))
      warning(nm, " failed: ", conditionMessage(e), call. = FALSE); NULL })
    record_stage_time(paste("Validation:", nm), proc.time()[["elapsed"]] - st)
    v
  }
  VALIDATION_RESULTS$sim <<- if (SIM_REPS > 0L) stage("Simulations", run_simulation_validation(annot_df)) else NULL
  VALIDATION_RESULTS$stab <<- if (STAB_N > 0L) stage("Stability", run_stability_analysis(annot_df)) else NULL
  VALIDATION_RESULTS$filter <<- if (length(FILTER_THRESHOLDS)) stage("Filter sensitivity", run_filter_sensitivity(annot_df)) else NULL
  add_validation("Validation replicates did not alter the observed-data caches (no leakage)", "All comparisons",
                 identical(before, cache_fingerprint()), "fingerprints of size factors, PAS sets and rankings")
  writeLines(c(sprintf("Profile %s. Simulations: %d per scenario x %d scenarios (%d PASs, base %s); stability: %d subsets (%s); filters: %s",
                       RUN_PROFILE, SIM_REPS, length(SIM_SCENARIOS), SIM_NPAS, SIM_BASE, STAB_N, STAB_METHOD,
                       paste(FILTER_THRESHOLDS, collapse = ", ")),
               if (length(status)) status else "All validation stages completed."),
             file.path(output_dir, "Validation_Status.txt"))
  invisible(TRUE)
}
VALIDATION_RESULTS <- list()

# =============================================================================
# MANUSCRIPT FIGURES 3, 4, 5 AND 7, SUPPLEMENTARY FIGURES S3-S4, AND SUMMARY
# TABLES (gained/lost discoveries, discovery totals, design and sequencing
# quality, computational performance, package lockfile)
# =============================================================================
# All panels are drawn from the tables written by the analysis; nothing is
# refitted here. Counts are descriptive; error rates and power come only from
# the simulations.
MODEL_COLS <- c("DESeq2 default" = "#9E9E9E", "Conventional DESeq2" = "#7F7F7F", "edgeR robust QL" = "#8C564B",
                "Pooled hierarchical" = "#2B2B2B", "Pooled hierarchical, TMM" = "#5E5E5E", "PC1 EVS" = "#1F6FB2",
                "Variance split" = "#E69F00", "PC1-PC2 split" = "#009E73", "Continuous rank" = "#CC79A7",
                "Random partition" = "#56B4E9", "Adaptive (pooled vs EVS)" = "#D55E00",
                "EVS, k = 5,000 (paper cutoff)" = "#4C9BD6", "Published analysis (k = 5,000)" = "#000000")

STAGE_TIMES <- list()
record_stage_time <- function(stage, seconds) {
  STAGE_TIMES[[length(STAGE_TIMES) + 1L]] <<- data.frame(Stage = stage, Seconds = round(seconds, 1),
                                                          Max_memory_MB = round(sum(gc()[, 6]), 1))
}

read_summary <- function(f) {
  p <- file.path(summary_table_dir, f)
  if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) else NULL
}

# Figure 3: empirical FDR and power of every method across simulation scenarios.
plot_fdr_power_figure <- function() {
  s <- read_summary("Table_Validation_Simulation_Summary.csv")
  if (is.null(s) || !nrow(s)) stop("no simulation results")
  s$Scenario <- factor(gsub("_", " ", s$Scenario), levels = gsub("_", " ", ALL_SCENARIOS))
  s$Model <- factor(s$Model, levels = names(MODEL_COLS))
  pa <- ggplot(s, aes(Scenario, FDP, colour = Model)) +
    geom_hline(yintercept = BH_FDR_STANDARD, linetype = "22", colour = "grey40") +
    geom_pointrange(aes(ymin = pmax(FDP_lo95, 0), ymax = FDP_hi95), position = position_dodge(width = 0.75), size = 0.15) +
    scale_colour_manual(values = MODEL_COLS, name = NULL, drop = TRUE) +
    labs(x = NULL, y = "Empirical FDR (Standard)", tag = "a", title = "Empirical FDR against the nominal 10%",
         subtitle = sprintf("Mean false discovery proportion (0 when there are no discoveries); empirical-Bernstein 95%% interval; dotted: %.3f.", BH_FDR_STANDARD + FDR_TOLERANCE)) +
    geom_hline(yintercept = BH_FDR_STANDARD + FDR_TOLERANCE, linetype = "13", colour = "grey55") +
    fig_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  pw <- s[grepl("^effects", s$Scenario), , drop = FALSE]
  pb <- ggplot(pw, aes(Scenario, TPR, colour = Model)) +
    geom_pointrange(aes(ymin = pmax(TPR_lo95, 0), ymax = pmin(TPR_hi95, 1)), position = position_dodge(width = 0.75), size = 0.15) +
    scale_colour_manual(values = MODEL_COLS, name = NULL, drop = TRUE) +
    labs(x = NULL, y = "Power (true-positive rate)", tag = "b", title = "Power at the nominal 10% FDR",
         subtitle = "EVS is better only if its power is higher with comparable empirical FDR.") +
    fig_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  assemble_rows(list(pa, pb), heights = c(1, 1), title = "Figure 3. Empirical FDR and power of competing methods (simulation)")
}

# Figure 4: differential PAS abundance: discoveries of every model, effect
# sizes of the Standard calls, and discoveries gained and lost by EVS.
plot_abundance_figure <- function() {
  mc <- if (length(model_comparison_rows)) dplyr::bind_rows(model_comparison_rows) else NULL
  if (is.null(mc) || !nrow(mc)) stop("no model comparison")
  mc$Model2 <- ifelse(grepl("^Random partition", mc$Model), "Random partition", mc$Model)
  db <- mc |> dplyr::group_by(Comparison, Model2) |> dplyr::summarise(Standard = stats::median(Standard), .groups = "drop")
  db$Model2 <- factor(db$Model2, levels = names(MODEL_COLS)); db$Comparison <- as_comp(db$Comparison)
  pa <- ggplot(db, aes(Comparison, Standard, fill = Model2)) +
    geom_col(position = position_dodge(width = 0.88), width = 0.82) +
    scale_fill_manual(values = MODEL_COLS, name = NULL, drop = TRUE) +
    scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Standard discoveries (BH < 0.10)", tag = "a", title = "Discoveries of every model",
         subtitle = "Counts are not error rates (see Figure 3). Random partitions: median.") +
    fig_theme() + theme(panel.grid.major.x = element_blank())
  eff <- dplyr::bind_rows(lapply(names(MODEL_FLAGS), function(cmp) dplyr::bind_rows(lapply(c("Pooled hierarchical", "PC1 EVS"), function(m) {
    f <- MODEL_FLAGS[[cmp]][[m]]
    if (is.null(f) || !"lfc_shrunk" %in% names(f)) return(NULL)
    f <- f[f$Standard %in% TRUE & is.finite(f$lfc_shrunk), , drop = FALSE]
    if (!nrow(f)) NULL else data.frame(Comparison = cmp, Model = m, LFC = f$lfc_shrunk)
  }))))
  rows <- list(pa)
  if (nrow(eff)) {
    eff$Comparison <- as_comp(eff$Comparison)
    rows[[2]] <- ggplot(eff, aes(LFC, fill = Model)) +
      geom_histogram(bins = 40, position = "identity", alpha = 0.55) +
      geom_vline(xintercept = c(-1, 1), linetype = "22", colour = "grey40", linewidth = 0.3) +
      facet_wrap(~ Comparison, nrow = 1, scales = "free_y", labeller = labeller(Comparison = comp_lab())) +
      scale_fill_manual(values = MODEL_COLS, name = NULL) +
      labs(x = "apeglm log2 fold change of Standard calls", y = "PASs", tag = "b", title = "Effect sizes of the discoveries") +
      fig_theme()
  }
  gl <- read_summary("Table_Gained_Lost_Discoveries.csv")
  if (!is.null(gl) && nrow(gl)) {
    g <- as.data.frame(table(Comparison = gl$Comparison, Change = gl$Change))
    g$Comparison <- as_comp(as.character(g$Comparison))
    rows[[length(rows) + 1L]] <- ggplot(g, aes(Comparison, Freq, fill = Change)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.75) +
      scale_fill_manual(values = c("Gained by EVS" = "#1F6FB2", "Lost by EVS" = "#C0392B"), name = NULL) +
      scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
      labs(x = NULL, y = "PASs", tag = letters[length(rows) + 1L], title = "Discoveries gained and lost by EVS relative to the pooled model") +
      fig_theme() + theme(panel.grid.major.x = element_blank())
  }
  assemble_rows(rows, heights = rep(1, length(rows)), title = "Figure 4. Differential PAS abundance and effect sizes")
}

# Figure 5: cutoff stability and model selection: k* across resamples, the
# pooled-vs-split bootstrap.
plot_stability_figure <- function() {
  rows <- list()
  st <- read_summary("Table_Cutoff_Stability.csv")
  if (!is.null(st) && nrow(st)) {
    st$Comparison <- as_comp(st$Comparison)
    rows[[length(rows) + 1L]] <- ggplot(st, aes(Comparison, k)) +
      geom_jitter(width = 0.12, height = 0, size = 1, colour = "grey35") +
      geom_point(aes(y = k_full), shape = 18, size = 2.6, colour = "#D55E00") +
      scale_x_discrete(labels = function(x) sub("_", "\n", x)) + scale_y_continuous(labels = scales::label_comma()) +
      labs(x = NULL, y = "k*", tag = letters[length(rows) + 1L], title = "Cutoff stability",
           subtitle = sprintf("Grey: k* after resampling animals within arms (%s); diamond: full data.", STAB_METHOD)) +
      fig_theme()
  }
  if (length(SELECTION)) {
    bs <- dplyr::bind_rows(lapply(names(SELECTION), function(cmp) data.frame(Comparison = cmp, LRT = SELECTION[[cmp]]$boot)))
    ob <- dplyr::bind_rows(lapply(names(SELECTION), function(cmp) data.frame(Comparison = cmp, LRT = SELECTION[[cmp]]$observed,
                                                                            Model = SELECTION[[cmp]]$model)))
    bs$Comparison <- as_comp(bs$Comparison); ob$Comparison <- as_comp(ob$Comparison)
    rows[[length(rows) + 1L]] <- ggplot() +
      geom_jitter(data = bs[is.finite(bs$LRT), , drop = FALSE], aes(Comparison, LRT), width = 0.12, height = 0, colour = "grey45", size = 1) +
      geom_point(data = ob, aes(Comparison, LRT, shape = Model), colour = "#D55E00", size = 2.6) +
      scale_shape_manual(values = c(EVS = 18, Pooled = 4), name = "Selected") +
      scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
      scale_y_continuous(trans = "log1p", labels = scales::label_comma()) +
      labs(x = NULL, y = "Coarse-grid LRT, split vs pooled", tag = letters[length(rows) + 1L], title = "Pooled versus split selection",
           subtitle = "Grey: complete adaptive search on data simulated from the fitted pooled model; orange: observed.") +
      fig_theme()
  }
  if (!length(rows)) stop("no stability or selection results")
  assemble_rows(rows, heights = rep(1, length(rows)), title = "Figure 5. Cutoff stability and model-selection robustness")
}

# Figure 7: 3'aTWAS enrichment (Firth odds ratios and 95% intervals) for every
# comparison, endpoint, model and ortholog scope; non-estimable rows are listed.
plot_twas_figure <- function() {
  d <- read_summary("Table_TWAS_Enrichment.csv")
  if (is.null(d) || !nrow(d)) stop("no TWAS enrichment results")
  d <- d[!grepl("cameraPR", d$Method), , drop = FALSE]
  d$Label <- paste(d$Comparison, d$Model, d$Scope, sep = " | ")
  d$Estimable <- is.finite(d$OR) & is.finite(d$Lower95) & is.finite(d$Upper95)
  d$Label <- factor(d$Label, levels = rev(unique(d$Label)))
  e <- d[d$Estimable, , drop = FALSE]
  p <- ggplot(e, aes(OR, Label, colour = Scope)) +
    geom_vline(xintercept = 1, linetype = "22", colour = "grey45") +
    geom_errorbarh(aes(xmin = Lower95, xmax = Upper95), height = 0, linewidth = 0.4) +
    geom_point(size = 1.6) +
    geom_text(aes(label = sprintf("%d/%d", Overlap, TWAS_genes)), x = Inf, hjust = 1.05, size = 1.8, colour = "grey20") +
    facet_wrap(~ Endpoint, ncol = 1, scales = "free_y") +
    scale_x_log10() + scale_colour_manual(values = c(all_supported = "#7B3294", one_to_one = "#56B4E9"), name = "Orthologs") +
    labs(x = "Odds ratio (95% CI)", y = NULL, tag = "a", title = "3'aTWAS enrichment",
         subtitle = paste0("Firth logistic regression; label: overlap / 3'aTWAS genes in the background.",
                           if (any(!d$Estimable)) sprintf(" Not estimable: %d rows (see table).", sum(!d$Estimable)) else "")) +
    fig_theme() + theme(axis.text.y = element_text(size = 5))
  assemble_rows(list(p), heights = 1, title = "Figure 7. 3'aTWAS enrichment: differential APA (primary) and PAS abundance (secondary)")
}

# Supplementary Figure S3: further simulation endpoints.
plot_simulation_supplement <- function() {
  s <- read_summary("Table_Validation_Simulation_Summary.csv")
  if (is.null(s) || !nrow(s)) stop("no simulation results")
  s$Scenario <- factor(gsub("_", " ", s$Scenario), levels = gsub("_", " ", ALL_SCENARIOS))
  s$Model <- factor(s$Model, levels = names(MODEL_COLS))
  one <- function(v, ylab, tag) ggplot(s[is.finite(s[[v]]), , drop = FALSE], aes(Scenario, .data[[v]], colour = Model)) +
    geom_point(position = position_dodge(width = 0.75), size = 1.1) +
    scale_colour_manual(values = MODEL_COLS, name = NULL, drop = TRUE) +
    labs(x = NULL, y = ylab, tag = tag) + fig_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  assemble_rows(list(list(plots = list(one("FNDR", "False non-discovery rate", "a"), one("Dispersion_log_RMSE", "Dispersion RMSE (log)", "b")), widths = c(1, 1)),
                     list(plots = list(one("LFC_RMSE", "log2FC RMSE", "c"), one("CI95_coverage", "Wald 95% CI coverage", "d")), widths = c(1, 1))),
                heights = c(1, 1), title = "Supplementary Figure S3. Simulation: non-discovery, dispersion and effect-size accuracy")
}

# Supplementary Figure S4: sensitivity to the abundance filter.
plot_filter_supplement <- function() {
  d <- read_summary("Table_Filter_Sensitivity.csv")
  if (is.null(d) || !nrow(d)) stop("no filter sensitivity results")
  primary <- sprintf(">= %d reads in >= %d samples", MIN_PAS_COUNT, MIN_PAS_SAMPLES)
  d <- rbind(d[, c("Comparison", "Filter", "Standard_pooled", "Standard_PC1_EVS")],
             dplyr::bind_rows(lapply(names(MODEL_FLAGS), function(cmp) data.frame(Comparison = cmp, Filter = paste(primary, "(primary)"),
               Standard_pooled = sum(MODEL_FLAGS[[cmp]][["Pooled hierarchical"]]$Standard),
               Standard_PC1_EVS = sum(MODEL_FLAGS[[cmp]][["PC1 EVS"]]$Standard)))))
  l <- rbind(data.frame(Comparison = d$Comparison, Filter = d$Filter, Model = "Pooled hierarchical", Standard = d$Standard_pooled),
             data.frame(Comparison = d$Comparison, Filter = d$Filter, Model = "PC1 EVS", Standard = d$Standard_PC1_EVS))
  l$Comparison <- as_comp(l$Comparison)
  p <- ggplot(l, aes(Filter, Standard, fill = Model)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    facet_wrap(~ Comparison, nrow = 1, scales = "free_y", labeller = labeller(Comparison = comp_lab())) +
    scale_fill_manual(values = MODEL_COLS, name = NULL) +
    labs(x = NULL, y = "Standard discoveries", tag = "a", title = "Sensitivity to the abundance filter") +
    fig_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  assemble_rows(list(p), heights = 1, title = "Supplementary Figure S4. Abundance filter sensitivity")
}

# Discoveries gained and lost by EVS relative to the pooled model, and totals
# (total calls, unique PASs across comparisons, per comparison).
save_discovery_tables <- function() {
  gl <- dplyr::bind_rows(lapply(names(MODEL_FLAGS), function(cmp) {
    o <- MODEL_FLAGS[[cmp]][["Pooled hierarchical"]]; e <- MODEL_FLAGS[[cmp]][["PC1 EVS"]]
    if (is.null(o) || is.null(e)) return(NULL)
    ids <- union(o$feature_id[o$Standard], e$feature_id[e$Standard])
    if (!length(ids)) return(NULL)
    d <- data.frame(Comparison = cmp, feature_id = ids,
                    Pooled_padj = o$padj_global[match(ids, o$feature_id)], EVS_padj = e$padj_global[match(ids, e$feature_id)],
                    Pooled_Standard = ids %in% o$feature_id[o$Standard], EVS_Standard = ids %in% e$feature_id[e$Standard],
                    EVS_view = e$View[match(ids, e$feature_id)], lfc_shrunk_EVS = e$lfc_shrunk[match(ids, e$feature_id)])
    d$Change <- ifelse(d$EVS_Standard & !d$Pooled_Standard, "Gained by EVS",
                       ifelse(!d$EVS_Standard & d$Pooled_Standard, "Lost by EVS", "Both"))
    d[d$Change != "Both", , drop = FALSE]
  }))
  gl$gene_symbol <- OrigID_Symbol$gene_symbol[match(gl$feature_id, OrigID_Symbol$feature_id)]
  save_csv(gl, file.path(summary_table_dir, "Table_Gained_Lost_Discoveries.csv"))
  mods <- unique(unlist(lapply(MODEL_FLAGS, names)))
  tot <- dplyr::bind_rows(lapply(mods, function(m) {
    per <- vapply(names(MODEL_FLAGS), function(cmp) { f <- MODEL_FLAGS[[cmp]][[m]]; if (is.null(f)) NA_integer_ else sum(f$Standard) }, integer(1))
    ids <- unique(unlist(lapply(MODEL_FLAGS, function(x) { f <- x[[m]]; if (is.null(f)) NULL else f$feature_id[f$Standard] })))
    cbind(data.frame(Model = m, Total_discoveries = sum(per, na.rm = TRUE), Unique_PASs = length(ids)), as.data.frame(as.list(per)))
  }))
  save_csv(tot, file.path(summary_table_dir, "Table_Discovery_Totals.csv"))
  # Direct comparison with the published analysis (k = 5,000).
  pap <- dplyr::bind_rows(lapply(names(MODEL_FLAGS), function(cmp) {
    fl <- MODEL_FLAGS[[cmp]]; pub <- fl[["Published analysis (k = 5,000)"]]; e <- fl[["PC1 EVS"]]; e5 <- fl[["EVS, k = 5,000 (paper cutoff)"]]
    if (is.null(pub)) return(NULL)
    std <- function(f) if (is.null(f)) character() else f$feature_id[f$Standard %in% TRUE]
    cs <- CUTOFF$summary[CUTOFF$summary$Comparison == cmp, , drop = FALSE]
    data.frame(Comparison = cmp, k_star = OPTIMAL_K[[cmp]] %||% NA_integer_,
               Leading_Edge_k_star = if (nrow(cs) && "Leading_Edge" %in% names(cs)) cs$Leading_Edge[1] else NA_integer_,
               Published_PAS_tested = nrow(pub), Published_Leading_Edge = attr(pub, "leading_edge") %||% NA_integer_,
               Published_Standard = length(std(pub)),
               Published_BH_only_no_LFC_gate = sum(!is.na(pub$padj_global) & pub$padj_global < BH_FDR_STANDARD),
               Published_Strong = sum(pub$Strong %in% TRUE),
               EVS_k5000_Standard = length(std(e5)), PC1_EVS_k_star_Standard = length(std(e)),
               PC1_EVS_and_Published = length(intersect(std(e), std(pub))),
               PC1_EVS_only = length(setdiff(std(e), std(pub))), Published_only = length(setdiff(std(pub), std(e))),
               Note = "Published Standard requires |apeglm log2FC| >= 1 (published definition); current Standard has no fold-change gate")
  }))
  if (nrow(pap)) save_csv(pap, file.path(summary_table_dir, "Table_Published_5K_Comparison.csv"))
  invisible(TRUE)
}

# Experimental design and sequencing quality (one row per sample and comparison).
save_design_quality_table <- function() {
  d <- dplyr::bind_rows(lapply(names(comparison_cache), function(cmp) {
    c <- comparison_cache[[cmp]]; raw <- comparison_raw_counts(cmp)
    data.frame(Comparison = cmp, Sample = colnames(c$raw), Arm = as.character(experiment_arm(colnames(c$raw))),
               Condition = as.character(c$coldata[colnames(c$raw), "condition"]),
               Total_reads = colSums(raw[, colnames(c$raw), drop = FALSE]), PAS_detected = colSums(raw[, colnames(c$raw), drop = FALSE] > 0),
               PAS_eligible = nrow(c$raw), Reads_in_eligible = colSums(c$raw), Zero_fraction_eligible = colMeans(c$raw == 0),
               Size_factor = unname(c$sf[colnames(c$raw)]), TMM_size_factor = unname(tmm_size_factors(c$raw)[colnames(c$raw)]))
  }))
  save_csv(d, file.path(summary_table_dir, "Table_Design_Sequencing_Quality.csv"))
}

# Computational performance and the renv lockfile. renv.lock records R, the
# Bioconductor release and every installed package with its version and
# source (CRAN, Bioconductor, GitHub), so renv::restore() can rebuild the
# library. renv is a required package. The lockfile is written into the output
# directory (and so into the final ZIP); a snapshot failure is reported in
# renv_status.txt and as a warning.
save_reproducibility_tables <- function() {
  if (length(STAGE_TIMES)) save_csv(dplyr::bind_rows(STAGE_TIMES), file.path(summary_table_dir, "Table_Computational_Performance.csv"))
  lock <- file.path(output_dir, "renv.lock")
  msg <- tryCatch({
      proj <- file.path(tempdir(), "sequence_evs_renv"); dir.create(proj, showWarnings = FALSE)
      if (requireNamespace("BiocManager", quietly = TRUE)) options(repos = BiocManager::repositories())
      renv::snapshot(project = proj, library = .libPaths(), lockfile = lock, type = "all", prompt = FALSE, force = TRUE)
      if (file.exists(lock)) "renv.lock written (restore with renv::restore(lockfile = 'renv.lock'))." else "renv::snapshot did not write renv.lock."
    }, error = function(e) paste("renv::snapshot failed:", conditionMessage(e)))
  if (!file.exists(lock)) warning(msg, call. = FALSE)
  writeLines(msg, file.path(output_dir, "renv_status.txt"))
}

# =============================================================================
# EXTERNAL VALIDATION ON INDEPENDENT PUBLIC PAS-LEVEL DATA
# =============================================================================
# Purpose. Apply the frozen EVS procedure (the same filter rule, ranking,
# likelihood, cutoff search, pooled-vs-split bootstrap, moderation and tests;
# no setting changed after these data were seen) to two independent public
# 3'-end datasets, and compare it with the conventional methods.
# Datasets (prespecified eligibility: public, PAS-level raw integer counts,
# biological replicates, at least three per group, a defined two-group
# contrast):
#  * GSE100349: rat WTTS-seq (the same 3'-end protocol and species as the
#    primary study), Control vs OB, four animals per group. Truth is unknown:
#    reported are discoveries of every method, k* and the pooled-vs-split
#    decision.
#  * GSE151173: mouse liver 3' mRNA-seq every 4 h (ZT02-ZT22), six animals per
#    time point per group; only the WT animals are used (each animal sampled
#    once). (i) Real-data null: at each time point the six WT animals are split
#    3 vs 3 (10 distinct splits per time point, 60 in all); same genotype, time
#    and tissue, so every discovery is false, with no simulation assumptions.
#    (ii) Replication with independent animals: for WT ZT02 vs ZT14
#    (antiphase, prespecified), the animals are split into two disjoint halves
#    (3 + 3 per time point each). Each method is applied to half A; half B
#    defines a method-independent reference, the PASs called (BH < 0.10) by
#    both DESeq2 default and edgeR robust QL. Reported: the share of the
#    reference recovered in A (recovery), and the share of A's discoveries with
#    the same sign and p < 0.05 in B under both conventional methods
#    (replication). The reference is built from conventional methods, so it
#    can favour them; it is not ground truth.
# Every adaptive analysis (each full comparison, null split and half) generates its own pooled-vs-split bootstrap inside reanalyze_matrix().
# Duplicate PAS identifiers, sample identities, group assignments and count
# integrity are checked first; a failed check stops that dataset.
# Files: SEQUENCE_EXTERNAL_WTTS (default data/GSE100349_Expression_table_WTTS-seq_rat.txt.gz)
# and SEQUENCE_EXTERNAL_LIVER (default data/GSE151173_PAS_RAW_Table.txt.gz);
# missing files are reported, not silently skipped.
EXTERNAL_ON <- !identical(Sys.getenv("SEQUENCE_EXTERNAL", "1"), "0")
EXT_MODELS <- c("deseq2_default", "edger", "pooled", "pc1", "adaptive")
EXT_NULL_N <- env_int("SEQUENCE_EXTERNAL_NULL_SPLITS", if (PUB) 60L else 3L)
EXT_REPL_N <- env_int("SEQUENCE_EXTERNAL_REPLICATION_SPLITS", if (PUB) 10L else 2L)
EXT_ROWS <- list(summary = list(), null = list(), repl = list())

external_path <- function(env, name) {
  cand <- c(Sys.getenv(env, ""), file.path(repo_root, "data", name), file.path(repo_root, name))
  cand <- cand[nzchar(cand) & file.exists(cand)]
  if (length(cand)) cand[1] else NA_character_
}

read_external_counts <- function(path, sample_pattern) {
  con <- gzfile(path); hdr <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]; close(con)
  if (!nzchar(hdr[1])) hdr[1] <- "ID"          # GSE151173: unnamed identifier column
  d <- utils::read.delim(path, header = FALSE, skip = 1L, quote = "", comment.char = "",
                         check.names = FALSE, stringsAsFactors = FALSE)
  d <- d[, !vapply(d, function(v) all(is.na(v) | v == ""), logical(1)), drop = FALSE]   # trailing tabs
  if (ncol(d) != length(hdr)) stop("Header and body widths differ in ", basename(path))
  names(d) <- hdr
  id <- d[[1]]
  smp <- grep(sample_pattern, names(d), value = TRUE)
  x <- as.matrix(d[, smp, drop = FALSE]); storage.mode(x) <- "numeric"
  keep <- !is.na(id) & nzchar(as.character(id))          # GSE100349 ends with an unlabelled total row
  x <- x[keep, , drop = FALSE]; id <- as.character(id[keep])
  if (anyDuplicated(id)) stop("Duplicate PAS identifiers in ", basename(path), ": ",
                              paste(utils::head(unique(id[duplicated(id)]), 5), collapse = ", "))
  if (anyDuplicated(smp)) stop("Duplicate sample columns in ", basename(path))
  if (anyNA(x) || any(x < 0) || any(x != round(x))) stop("External counts are not non-negative integers: ", basename(path))
  rownames(x) <- id
  x
}

ext_coldata <- function(samples, trt) data.frame(condition = factor(ifelse(trt, "trt", "untrt"), levels = c("untrt", "trt")),
                                                 row.names = samples)

# Metadata checks that stop the affected external analysis when they fail.
ext_require <- function(check, dataset, ok, detail = "") {
  add_validation(check, dataset, ok, detail)
  if (!isTRUE(ok)) stop(check, " failed: ", detail)
}

ext_record <- function(dataset, contrast, r) {
  for (m in names(r$calls)) EXT_ROWS$summary[[length(EXT_ROWS$summary) + 1L]] <<- data.frame(
    Dataset = dataset, Contrast = contrast, Model = MODEL_LABELS[[m]], PAS_eligible = r$n_pas,
    k = if (is.null(r$k[[m]])) NA_integer_ else r$k[[m]], Standard = sum(r$calls[[m]]$Standard %in% TRUE),
    Strong = if (all(is.na(r$calls[[m]]$Strong))) NA_integer_ else sum(r$calls[[m]]$Strong %in% TRUE),
    Selected = if (m == "adaptive") r$selected else NA_character_,
    Coarse_LRT = r$lrt, Bootstrap_replicates = r$boot_n, Bootstrap_successful = r$boot_ok,
    Selection_p = r$selection_p,
    Seconds = r$total_seconds)
}

run_external_validation <- function(annot_df) {
  status <- character()
  ann <- data.frame(feature_id = character(), gene_symbol = character())
  # ---- GSE100349: rat WTTS-seq, Control vs OB ----
  p1 <- external_path("SEQUENCE_EXTERNAL_WTTS", "GSE100349_Expression_table_WTTS-seq_rat.txt.gz")
  if (is.na(p1)) status <- c(status, "GSE100349: file not found") else tryCatch({
    x <- read_external_counts(p1, "mapped reads$")
    expected <- c(sprintf("Rat-Control%d-WTTS mapped reads", 1:4), sprintf("Rat-OB%d-WTTS mapped reads", 1:4))
    ext_require("External GSE100349: sample identities and count columns (Control1-4, OB1-4, one column each)", "GSE100349",
                setequal(colnames(x), expected) && !anyDuplicated(colnames(x)), paste(colnames(x), collapse = ", "))
    trt <- grepl("-OB[0-9]", colnames(x)); cd <- ext_coldata(colnames(x), trt)
    ext_require("External GSE100349: four animals per group", "GSE100349", sum(trt) == 4L && sum(!trt) == 4L)
    ext_require("External GSE100349: PAS identifiers unique", "GSE100349", !anyDuplicated(rownames(x)))
    r <- reanalyze_matrix(x, cd, "GSE100349 OB vs Control", EXT_MODELS, ann, SPLIT_SEED + 70002L)
    ext_record("GSE100349 (rat WTTS-seq)", "OB vs Control", r)
  }, error = function(e) status <<- c(status, paste("GSE100349:", conditionMessage(e))))
  # ---- GSE151173: mouse liver, WT only ----
  p2 <- external_path("SEQUENCE_EXTERNAL_LIVER", "GSE151173_PAS_RAW_Table.txt.gz")
  if (is.na(p2)) status <- c(status, "GSE151173: file not found") else tryCatch({
    x <- read_external_counts(p2, "^WT[0-9]+ZT[0-9]+$")
    zt <- sub("^WT[0-9]+", "", colnames(x)); an <- sub("ZT[0-9]+$", "", colnames(x))
    expected <- as.vector(outer(sprintf("WT%d", 1:6), sprintf("ZT%02d", seq(2, 22, by = 4)), paste0))
    ext_require("External GSE151173: WT sample identities and count columns (WT1-6 x ZT02-ZT22, one column each)", "GSE151173",
                setequal(colnames(x), expected) && !anyDuplicated(colnames(x)), sprintf("%d WT columns", ncol(x)))
    ext_require("External GSE151173: six WT animals at each of six time points", "GSE151173",
                all(table(zt) == 6L) && length(unique(zt)) == 6L, paste(names(table(zt)), table(zt), collapse = "; "))
    ext_require("External GSE151173: PAS identifiers unique", "GSE151173", !anyDuplicated(rownames(x)))
    # (i) Real-data null splits (3 vs 3 within a time point).
    a <- colnames(x)[zt == "ZT02"]; b <- colnames(x)[zt == "ZT14"]
    cdp <- ext_coldata(c(b, a), c(b, a) %in% b)
    rpos <- reanalyze_matrix(x[, c(b, a)], cdp, "GSE151173 WT ZT14 vs ZT02", EXT_MODELS, ann, SPLIT_SEED + 80002L)
    ext_record("GSE151173 (mouse liver 3' mRNA-seq)", "WT ZT14 vs ZT02 (6 vs 6)", rpos)
    splits <- do.call(c, lapply(sort(unique(zt)), function(t) {
      s <- colnames(x)[zt == t]
      lapply(Filter(function(g) s[1] %in% g, utils::combn(s, 3L, simplify = FALSE)), function(g) list(t = t, s = s, g = g))
    }))
    set.seed(SPLIT_SEED + 80003L)
    use <- if (EXT_NULL_N >= length(splits)) seq_along(splits) else sort(sample.int(length(splits), EXT_NULL_N))
    for (j in use) {
      sp <- splits[[j]]
      cdn <- ext_coldata(sp$s, sp$s %in% sp$g)
      rn <- tryCatch(reanalyze_matrix(x[, sp$s], cdn, sprintf("GSE151173 null %s %d", sp$t, j),
                                      EXT_MODELS, ann, SPLIT_SEED + 80100L + j), error = function(e) NULL)
      if (is.null(rn)) { status <- c(status, sprintf("GSE151173 null split %d failed", j)); next }
      for (m in names(rn$calls)) EXT_ROWS$null[[length(EXT_ROWS$null) + 1L]] <<- cbind(
        data.frame(Dataset = "GSE151173 (mouse liver 3' mRNA-seq)", Null_type = sprintf("WT %s, 3 vs 3 (true null)", sp$t),
                   Split = j, Group1 = paste(sp$g, collapse = ","), Model = MODEL_LABELS[[m]],
                   Selected = if (m == "adaptive") rn$selected else NA_character_), null_metrics(rn$calls[[m]]))
      message(sprintf("[external] GSE151173 null split %d done", j))
    }
    # (ii) Replication across disjoint halves of the animals.
    halves <- lapply(Filter(function(g) a[1] %in% g, utils::combn(a, 3L, simplify = FALSE)), function(ga)
      lapply(utils::combn(b, 3L, simplify = FALSE), function(gb) list(a = ga, b = gb)))
    halves <- do.call(c, halves)
    set.seed(SPLIT_SEED + 80004L)
    use <- if (EXT_REPL_N >= length(halves)) seq_along(halves) else sort(sample.int(length(halves), EXT_REPL_N))
    for (j in use) {
      hA <- c(halves[[j]]$b, halves[[j]]$a); hB <- c(setdiff(b, halves[[j]]$b), setdiff(a, halves[[j]]$a))
      cdA <- ext_coldata(hA, hA %in% b)
      rA <- tryCatch(reanalyze_matrix(x[, hA], cdA, sprintf("GSE151173 half A %d", j), EXT_MODELS, ann,
                                      SPLIT_SEED + 80200L + j), error = function(e) NULL)
      rB <- tryCatch(reanalyze_matrix(x[, hB], ext_coldata(hB, hB %in% b), sprintf("GSE151173 half B %d", j),
                                      c("deseq2_default", "edger"), ann, SPLIT_SEED + 80300L + j), error = function(e) NULL)
      if (is.null(rA) || is.null(rB)) { status <- c(status, sprintf("GSE151173 replication split %d failed", j)); next }
      d1 <- rB$calls$deseq2_default; e1 <- rB$calls$edger
      ref <- intersect(d1$feature_id[d1$Standard %in% TRUE], e1$feature_id[e1$Standard %in% TRUE])
      for (m in names(rA$calls)) {
        f <- rA$calls[[m]]; disc <- f$feature_id[f$Standard %in% TRUE]
        sgn <- sign(f$log2FoldChange[match(disc, f$feature_id)])
        rep_ok <- function(g) { i <- match(disc, g$feature_id)
          !is.na(i) & is.finite(g$pvalue[i]) & g$pvalue[i] < 0.05 & sign(g$log2FoldChange[i]) == sgn }
        rp <- rep_ok(d1) & rep_ok(e1)
        EXT_ROWS$repl[[length(EXT_ROWS$repl) + 1L]] <<- data.frame(
          Dataset = "GSE151173 (mouse liver 3' mRNA-seq)", Contrast = "WT ZT14 vs ZT02, 3 vs 3 per half", Split = j,
          Model = MODEL_LABELS[[m]], Selected = if (m == "adaptive") rA$selected else NA_character_,
          Discoveries_A = length(disc), Reference_B = length(ref),
          Recovered = length(intersect(disc, ref)), Recovery = if (length(ref)) length(intersect(disc, ref)) / length(ref) else NA_real_,
          Replicated = sum(rp), Replication_rate = if (length(disc)) mean(rp) else NA_real_)
      }
      message(sprintf("[external] GSE151173 replication split %d done", j))
    }
  }, error = function(e) status <<- c(status, paste("GSE151173:", conditionMessage(e))))
  if (length(EXT_ROWS$summary)) save_csv(dplyr::bind_rows(EXT_ROWS$summary), file.path(summary_table_dir, "Table_External_Validation_Discoveries.csv"))
  if (length(EXT_ROWS$null)) {
    d <- dplyr::bind_rows(EXT_ROWS$null)
    save_csv(d, file.path(summary_table_dir, "Table_External_Validation_Null.csv"))
    s <- d |> dplyr::group_by(Dataset, Null_type, Model) |>
      dplyr::summarise(Splits = dplyr::n(), Mean_Standard = mean(Standard), Max_Standard = max(Standard),
                       Prob_any_discovery = mean(Any_Standard), Mean_TypeI_p05 = mean(TypeI_p05),
                       Mean_TypeI_p01 = mean(TypeI_p01), Median_KS_D = stats::median(KS_D), .groups = "drop")
    save_csv(as.data.frame(s), file.path(summary_table_dir, "Table_External_Validation_Null_Summary.csv"))
  }
  if (length(EXT_ROWS$repl)) {
    d <- dplyr::bind_rows(EXT_ROWS$repl)
    save_csv(d, file.path(summary_table_dir, "Table_External_Validation_Replication.csv"))
  }
  writeLines(c(sprintf("Profile %s. GSE151173 null splits %d, replication splits %d.",
                       RUN_PROFILE, EXT_NULL_N, EXT_REPL_N),
               if (length(status)) status else "All external validation stages completed."),
             file.path(output_dir, "External_Validation_Status.txt"))
  invisible(TRUE)
}

# Supplementary Figure S5: external validation.
plot_external_supplement <- function() {
  rows <- list()
  d <- read_summary("Table_External_Validation_Discoveries.csv")
  if (!is.null(d) && nrow(d)) {
    d$Model <- factor(d$Model, levels = names(MODEL_COLS))
    rows[[length(rows) + 1L]] <- ggplot(d, aes(Contrast, Standard, fill = Model)) +
      geom_col(position = position_dodge(width = 0.85), width = 0.8) +
      facet_wrap(~ Dataset, scales = "free") + scale_fill_manual(values = MODEL_COLS, name = NULL, drop = TRUE) +
      scale_y_continuous(labels = scales::label_comma()) +
      labs(x = NULL, y = "Standard discoveries", tag = letters[length(rows) + 1L], title = "Discoveries on independent data",
           subtitle = "Truth unknown; counts are not error rates.") + fig_theme()
  }
  n <- read_summary("Table_External_Validation_Null.csv")
  if (!is.null(n) && nrow(n)) {
    n$Model <- factor(n$Model, levels = names(MODEL_COLS))
    rows[[length(rows) + 1L]] <- ggplot(n, aes(Model, Standard + 1, colour = Model)) +
      geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.8) +
      scale_y_log10() + scale_colour_manual(values = MODEL_COLS, guide = "none") +
      labs(x = NULL, y = "Standard discoveries + 1", tag = letters[length(rows) + 1L], title = "False discoveries on real null data (GSE151173, WT 3 vs 3 within a time point)") +
      fig_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }
  r <- read_summary("Table_External_Validation_Replication.csv")
  if (!is.null(r) && nrow(r)) {
    r$Model <- factor(r$Model, levels = names(MODEL_COLS))
    rl <- rbind(data.frame(Model = r$Model, Measure = "Recovery of half-B reference", Value = r$Recovery),
                data.frame(Model = r$Model, Measure = "Replication in half B", Value = r$Replication_rate))
    rows[[length(rows) + 1L]] <- ggplot(rl[is.finite(rl$Value), , drop = FALSE], aes(Model, Value, colour = Model)) +
      geom_jitter(width = 0.15, height = 0, size = 1.2) + facet_wrap(~ Measure) +
      scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
      scale_colour_manual(values = MODEL_COLS, guide = "none") +
      labs(x = NULL, y = NULL, tag = letters[length(rows) + 1L], title = "Replication with independent animals (GSE151173, WT ZT14 vs ZT02)",
           subtitle = "Reference: PASs called by both DESeq2 default and edgeR in the other half (not ground truth).") +
      fig_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }
  if (!length(rows)) stop("no external validation results")
  assemble_rows(rows, heights = rep(1, length(rows)), title = "Supplementary Figure S5. External validation on independent public PAS-level data")
}

# =============================================================================
# UNIT TESTS ON SYNTHETIC DATA (SEQUENCE_SELFTEST=1)
# =============================================================================
# Run only on request (SEQUENCE_SELFTEST=1): they use small simulated data, not
# the WTTS-Seq counts, write Table_Unit_Tests.csv and the synthetic dataset
# (Synthetic_Test_Counts.csv) into a separate "_selftest" output directory,
# and stop. The analysis run itself does not execute them.
SELFTEST <- identical(Sys.getenv("SEQUENCE_SELFTEST", "0"), "1")

synthetic_dataset <- function(n = 1200L, n_per = 5L, seed = 1L) {
  set.seed(seed)
  smp <- c(sprintf("RTX_%d", seq_len(n_per)), sprintf("ZTX_%d", seq_len(n_per)))
  coldata <- data.frame(condition = factor(rep(c("trt", "untrt"), each = n_per), levels = c("untrt", "trt")), row.names = smp)
  mu <- exp(stats::runif(n, log(20), log(2000)))
  hi <- seq_len(n) <= n %/% 3                      # two dispersion populations
  alpha <- ifelse(hi, 0.25, 0.02) * exp(stats::rnorm(n, 0, 0.3))
  beta <- ifelse(stats::runif(n) < 0.05, sample(c(-1, 1), n, replace = TRUE), 0)
  sf <- stats::setNames(exp(stats::rnorm(length(smp), 0, 0.15)), smp)
  M <- outer(mu, sf) * 2^outer(beta, as.numeric(coldata$condition == "trt"))
  y <- matrix(stats::rnbinom(length(M), mu = as.vector(M), size = rep(1 / alpha, times = ncol(M))), n,
              dimnames = list(sprintf("PAS%05d", seq_len(n)), smp))
  list(y = y, coldata = coldata, alpha = stats::setNames(alpha, rownames(y)), beta = stats::setNames(beta, rownames(y)),
       high_dispersion = stats::setNames(hi, rownames(y)))
}

run_unit_tests <- function(out_dir) {
  res <- list()
  test <- function(name, expr) {
    st <- proc.time()[["elapsed"]]
    v <- tryCatch(expr, error = function(e) structure(FALSE, detail = conditionMessage(e)))
    res[[length(res) + 1L]] <<- data.frame(Test = name, Status = if (isTRUE(v)) "PASS" else "FAIL",
                                           Detail = as.character(attr(v, "detail") %||% ""),
                                           Seconds = round(proc.time()[["elapsed"]] - st, 2))
  }
  syn <- synthetic_dataset()
  save_csv(data.frame(feature_id = rownames(syn$y), syn$y, true_alpha = syn$alpha, true_log2FC = syn$beta,
                      high_dispersion_population = syn$high_dispersion, check.names = FALSE),
           file.path(out_dir, "Synthetic_Test_Counts.csv"))
  # Synthetic likelihood curves from known dispersions.
  set.seed(2)
  n <- 400L; lx <- stats::runif(n, log(10), log(5000))
  la <- log(0.02 + 0.5 * exp(-0.8 * lx)) + stats::rnorm(n, 0, 0.5)
  lk <- t(vapply(la, function(a) exp(-0.5 * ((DISP_LOG_GRID - a) / 0.6)^2), numeric(length(DISP_LOG_GRID))))
  test("Analytic gradient of the marginal likelihood matches central differences", {
    e <- gradient_check(c(log(0.02), log(0.5), 0.8, log(0.3)), lx, lk)
    structure(e < 1e-4, detail = sprintf("max relative error %.2e", e)) })
  test("Analytic gradient with a spline basis matches central differences", {
    B <- cbind(scale(lx)[, 1], scale(lx^2)[, 1])
    e <- gradient_check(c(log(0.02), log(0.5), 0.8, log(0.3), 0.1, -0.1), lx, lk, B)
    structure(e < 1e-4, detail = sprintf("max relative error %.2e", e)) })
  test("Marginal fit recovers the prior variance (true 0.25)", {
    f <- fit_marginal_trend(lx, lk, multistart = TRUE); t2 <- f$tau2
    structure(is.finite(t2) && abs(t2 - 0.25) < 0.15, detail = sprintf("tau2 %.3f", t2)) })
  test("Profiled group means solve the score equations", {
    grp <- factor(as.character(syn$coldata$condition)); s <- rep(1, ncol(syn$y)); a <- 0.1
    pm <- profile_group_means(syn$y[1:200, ], s, grp, a)
    sc <- vapply(levels(grp), function(g) { j <- grp == g
      max(abs(rowSums((syn$y[1:200, j] - pm$mu[, j]) / (1 + a * pm$mu[, j])))) }, numeric(1))
    structure(max(sc) < 1e-4, detail = sprintf("max |score| %.2e", max(sc))) })
  test("Abundance filter keeps PASs with >= 10 reads in >= 5 samples", {
    x <- rbind(a = c(10, 10, 10, 10, 10, 0, 0, 0, 0, 0), b = c(10, 10, 10, 10, 9, 0, 0, 0, 0, 0))
    identical(unname(eligible_pas(x, 10, 5)), c(TRUE, FALSE)) })
  test("Leading Edge and Remainder are disjoint and complete", {
    d <- stats::setNames(sample.int(nrow(syn$y)), rownames(syn$y))
    m <- split_matrices(syn$y, d, 400L)
    a <- rownames(m$leading_edge_dataset); b <- rownames(m$remainder_dataset)
    !length(intersect(a, b)) && setequal(c(a, b), rownames(syn$y)) })
  test("BH is applied jointly across the two views", {
    mk <- function(ids, p) data.frame(feature_id = ids, pvalue = p, resGA_pvalue = p, resLA_pvalue = p, lfc_shrunk = 0.5,
                                      standard_flag = FALSE, strong_cnh_flag = FALSE, weak_significant_flag = FALSE, hbfss_flag = FALSE)
    p1 <- stats::runif(50)^3; p2 <- stats::runif(80)
    f <- discovery_flags(list(leading_edge_dataset = mk(sprintf("a%d", 1:50), p1), remainder_dataset = mk(sprintf("b%d", 1:80), p2)))
    isTRUE(all.equal(f$padj_global, stats::p.adjust(c(p1, p2), "BH"))) })
  test("Simes combination", isTRUE(all.equal(simes_p(c(0.01, 0.04, 0.5)), min(c(0.01 * 3, 0.04 * 3 / 2, 0.5)))))
  test("Cutoff search on two simulated dispersion populations finds a split near the truth", {
    cache <- build_comparison_cache(syn$y[eligible_pas(syn$y), ], syn$coldata)
    si <- split_inputs(cache)
    one <- fit_marginal_trend(si$lx, si$lk, multistart = TRUE)
    true_rank <- stats::setNames(ifelse(syn$high_dispersion[names(cache$depth)], 1L, 2L) * 1e4 + rank(cache$depth), names(cache$depth))
    d <- stats::setNames(rank(true_rank, ties.method = "first"), names(true_rank))
    s <- select_cutoff(si, d, one)
    nh <- sum(syn$high_dispersion[names(d)])
    structure(abs(s$k - nh) / nh < 0.25 && s$best$le$tau2 > 0, detail = sprintf("k %d, true %d", s$k, nh)) })
  test("Deterministic: identical seeds give identical simulated counts", {
    b <- list(ids = rownames(syn$y), mean = stats::setNames(rowMeans(syn$y) + 1, rownames(syn$y)),
              unit = stats::setNames(rownames(syn$y), rownames(syn$y)), in_le = syn$high_dispersion,
              sf = stats::setNames(rep(1, 10), colnames(syn$y)), coldata = syn$coldata,
              one = list(f = function(m) 0.05 + 1 / m, tau2 = 0.2), le = list(f = function(m) 0.2 + 1 / m, tau2 = 0.2),
              rem = list(f = function(m) 0.02 + 1 / m, tau2 = 0.2))
    identical(simulate_counts(b, "effects_two_populations", 7L)$y, simulate_counts(b, "effects_two_populations", 7L)$y) })
  d <- dplyr::bind_rows(res)
  save_csv(d, file.path(out_dir, "Table_Unit_Tests.csv"))
  d
}

# Reproducibility manifest: exact input hashes, software versions, settings.
write_provenance <- function() {
  writeLines(capture.output(sessionInfo()),file.path(output_dir,"SessionInfo.txt"))
  files <- c(count_file,resolve_existing_file(twas_file_candidates,"3'aTWAS file"),get_script_path())
  ext <- c(external_path("SEQUENCE_EXTERNAL_WTTS", "GSE100349_Expression_table_WTTS-seq_rat.txt.gz"),
           external_path("SEQUENCE_EXTERNAL_LIVER", "GSE151173_PAS_RAW_Table.txt.gz"))
  files <- c(files, ext[!is.na(ext)])
  save_csv(data.frame(Path=files,MD5=unname(tools::md5sum(files))),file.path(output_dir,"Input_Checksums.csv"))
  writeLines(c(paste("Build",PIPELINE_BUILD),paste("Grid step",DISP_GRID_STEP),paste("Minimum prior variance",DISP_PRIOR_VAR_MIN),
               paste("Cores",SPLIT_CORES),paste("Optimizer max iterations",TREND_MAXIT),
               paste("Random partitions",MODEL_RANDOM_PARTITIONS),paste("Random seed",SPLIT_SEED),paste("Minimum PAS count",MIN_PAS_COUNT),
               paste("Minimum PAS samples",MIN_PAS_SAMPLES)),file.path(output_dir,"Analysis_Settings.txt"))
}

if (SELFTEST) {
  ut <- run_unit_tests(output_dir)
  print(ut)
  message(sprintf("Unit tests: %d of %d passed.", sum(ut$Status == "PASS"), nrow(ut)))
  quit(save = "no", status = if (all(ut$Status == "PASS")) 0L else 1L)
}

# Experimental design checks before any analysis.
run_design_checks <- function() {
  for (i in seq_len(nrow(comparison_table))) {
    cmp <- comparison_table$comparison_name[i]
    ids <- meta_all$id[grepl(paste0("^(", comparison_table$group1_prefix[i], "|", comparison_table$group2_prefix[i], ")_"), meta_all$id)]
    cd <- meta_all[meta_all$id %in% ids, , drop = FALSE]
    tab <- table(cd$condition)
    add_validation("Sample identities: five animals per arm, arms matched to conditions", cmp,
                   all(tab == 5L) && length(tab) == 2L &&
                     all(grepl(paste0("^", comparison_table$group1_prefix[i], "_"), cd$id[cd$condition == "trt"])),
                   paste(names(tab), tab, collapse = "; "))
    add_validation("No duplicated samples", cmp, !anyDuplicated(cd$id) && all(cd$id %in% colnames(WTTS_Seq)))
    x <- as.matrix(WTTS_Seq[, rownames(cd), drop = FALSE])
    r <- suppressWarnings(stats::cor(log1p(x))); diag(r) <- NA
    add_validation("Duplicate and near-identical sequencing library screening", cmp,
                   !anyDuplicated(t(x)) && max(r, na.rm = TRUE) < 0.9999,
                   sprintf("max between-sample correlation of log counts %.4f (independence itself rests on the design: one animal per sample)", max(r, na.rm = TRUE)))
    add_validation("Counts are non-negative integers", cmp, all(is.finite(x)) && all(x >= 0) && all(x == round(x)))
    X <- stats::model.matrix(~ condition, data = droplevels(cd))
    add_validation("Design matrix has full rank", cmp, qr(X)$rank == ncol(X), sprintf("rank %d of %d", qr(X)$rank, ncol(X)))
  }
}
.st <- proc.time()[["elapsed"]]
run_design_checks()
record_stage_time("Data loading and design checks", proc.time()[["elapsed"]] - .st)

# The cutoff k* of every comparison, before any DE test (EVS CUTOFF k*).
cutoff_calls_before <- c(SIGNIFICANCE_CALLS$wald, SIGNIFICANCE_CALLS$null)
.st <- proc.time()[["elapsed"]]
compute_cutoffs()
record_stage_time("Cutoff selection and pooled-vs-split bootstrap", proc.time()[["elapsed"]] - .st)
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
  .st <- proc.time()[["elapsed"]]

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

  record_stage_time(paste("Comparison", cmp, "(fits, tests, model comparison, sensitivity)"), proc.time()[["elapsed"]] - .st)
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
.st <- proc.time()[["elapsed"]]
for (nm in c("Time_course","APA","TWAS_enrichment","APA_figure","APA_detail_figure")) tryCatch({
  if (nm=="Time_course") run_time_course() else
  if (nm=="APA") run_apa_analysis() else
  if (nm=="APA_figure") { if (is.null(APA_RESULTS)) stop("APA analysis did not run")
    save_figure(plot_apa_figure(), paper_fig_dir, "Figure_6_APA", 180, 170) } else
  if (nm=="APA_detail_figure") { if (is.null(APA_RESULTS)) stop("APA analysis did not run")
    save_figure(plot_apa_detail_figure(), paper_fig_dir, "Supplementary_Figure_S2_APA_Detail", 180, 170) } else
  if (is.null(twas_results)) stop("3'aTWAS overlap did not run") else
    run_twas_enrichment(twas_results$mapping)
},error=function(e) {complementary_errors <<- c(complementary_errors,paste(nm,conditionMessage(e),sep=": ")); warning(tail(complementary_errors,1))})
writeLines(if(length(complementary_errors)) complementary_errors else "All complementary stages completed (time course, APA, 3'aTWAS enrichment (APA primary, abundance secondary), Figure 6 and Supplementary Figure S2).",file.path(output_dir,"Complementary_Status.txt"))
record_stage_time("Time course, APA and 3'aTWAS", proc.time()[["elapsed"]] - .st)
tryCatch(save_multiplicity_table(), error = function(e) warning("Multiplicity table skipped: ", conditionMessage(e), call. = FALSE))
for (f in list(save_discovery_tables, save_design_quality_table))
  tryCatch(f(), error = function(e) warning("Summary table skipped: ", conditionMessage(e), call. = FALSE))
if (length(SELECTION)) save_csv(dplyr::bind_rows(lapply(names(SELECTION), function(cmp) data.frame(
  Comparison = cmp, Coarse_LRT_observed = SELECTION[[cmp]]$observed, Bootstrap_replicates = sum(is.finite(SELECTION[[cmp]]$boot)),
  Critical_value = SELECTION[[cmp]]$crit, Selection_p = SELECTION[[cmp]]$p, Selected_Model = SELECTION[[cmp]]$model,
  Forced_EVS_reported = TRUE))), file.path(summary_table_dir, "Table_Model_Selection.csv"))
if (VALIDATION_ON) run_validation(OrigID_Symbol) else
  writeLines("Validation skipped (SEQUENCE_VALIDATION=0).", file.path(output_dir, "Validation_Status.txt"))
if (EXTERNAL_ON) {
  .st <- proc.time()[["elapsed"]]
  tryCatch(run_external_validation(OrigID_Symbol), error = function(e) warning("External validation failed: ", conditionMessage(e), call. = FALSE))
  record_stage_time("External validation (GSE100349, GSE151173)", proc.time()[["elapsed"]] - .st)
} else writeLines("External validation skipped (SEQUENCE_EXTERNAL=0).", file.path(output_dir, "External_Validation_Status.txt"))
try_fig(save_figure(plot_external_supplement(), paper_fig_dir, "Supplementary_Figure_S5_External_Validation", 180, 200), "Supplementary Figure S5")
try_fig(save_figure(plot_fdr_power_figure(), paper_fig_dir, "Figure_3_FDR_Power", 180, 170), "Figure 3 (FDR and power)")
try_fig(save_figure(plot_abundance_figure(), paper_fig_dir, "Figure_4_Differential_Abundance", 180, 200), "Figure 4 (differential abundance)")
try_fig(save_figure(plot_stability_figure(), paper_fig_dir, "Figure_5_Stability_Selection", 180, 200), "Figure 5 (stability and selection)")
try_fig(save_figure(plot_twas_figure(), paper_fig_dir, "Figure_7_TWAS", 180, 200), "Figure 7 (3'aTWAS)")
try_fig(save_figure(plot_simulation_supplement(), paper_fig_dir, "Supplementary_Figure_S3_Simulation_Details", 180, 160), "Supplementary Figure S3")
try_fig(save_figure(plot_filter_supplement(), paper_fig_dir, "Supplementary_Figure_S4_Filter_Sensitivity", 180, 90), "Supplementary Figure S4")
tryCatch(save_run_checks(), error = function(e) warning("Validation outputs skipped: ", conditionMessage(e), call. = FALSE))
save_reproducibility_tables()
write_provenance()


cat("\n=====================================================\n")
cat("Pipeline complete.\n")
cat("Build: ", PIPELINE_BUILD, "\n", sep = "")
cat("Cutoff method: ", CUTOFF_METHOD_LABEL, "\n", sep = "")
cat("Repository root:\n", repo_root, "\n", sep = "")
cat("Output directory:\n", output_dir, "\n", sep = "")
cat("3'aTWAS overlap: ", if (is.null(twas_results)) "FAILED (see warnings)" else "complete", "\n", sep = "")
cat("=====================================================\n\n")

if (nrow(overall_summary) > 0L) print(overall_summary)
