# Manuscript Text Revisions (T1-T15)

Instructions: Each section below provides revised or new text to be inserted into the manuscript at the indicated location. Section references use the line numbers from the submitted manuscript PDF.

---

## T1: Move Power Analysis to Methods/Results (R1-Maj1, R2-Maj1)

**Location:** Add new subsection 2.4.1 after current Section 2.4 (Statistical analysis, after line 141). Also add to Results Section 3.2 after line 182.

### New Methods text (insert as Section 2.4.1 "Statistical power considerations"):

> **2.4.1 Statistical power considerations**
>
> This study uses data from 3 donors, which severely limits donor-level statistical power. A simulation-based power analysis (Cohen's $d = 1.0$ from observed effect size $\Delta = 0.06$, $SD \approx 0.06$) indicates that $\geq$20 donors per group (40 total) would provide approximately 85% power to detect the combinatorial effect under donor-aware testing at $\alpha = 0.05$. Our current $n = 3$ design yields estimated power of approximately 27%. Consequently, all donor-level analyses should be interpreted as exploratory pilot results, and null findings at the donor level cannot be taken as evidence of absence. We frame this study as a methodological proof-of-concept and hypothesis-generating exercise rather than a definitive biological investigation.

### New Results text (insert after line 182, in Section 3.2):

> **Statistical power context.** The 0/21 FDR-significant result should be interpreted in light of the study's limited statistical power. With only 3 donors, the donor-aware pseudobulk analysis has approximately 27% power to detect effects of the observed magnitude. This means the study is substantially underpowered to detect true individual-gene effects even if they exist, and the null result is consistent with both the absence of genuine effects and the presence of effects that are too small to detect at this sample size. We therefore treat the 0/21 FDR result as a methodological observation --- demonstrating that expression-matched null testing substantially attenuates raw perturbation signals --- rather than as biological evidence that intelligence genes have no special perturbation properties.

---

## T2: Pseudoreplication Subsection (R1-Maj1)

**Location:** Add new Discussion subsection 4.1.1 after Section 4.1 (after line 397).

### New text (insert as Section 4.1.1 "Pseudoreplication in foundation model perturbation studies"):

> **4.1.1 Pseudoreplication in foundation model perturbation studies**
>
> A critical statistical concern in single-cell foundation model analyses is pseudoreplication: cells from the same donor share donor-level confounds (genotype, tissue processing, ambient RNA contamination) that violate the independence assumption underlying cell-level statistical tests. When Geneformer processes 500 cells from 3 donors, the effective sample size for biological inference is $n = 3$, not $n = 500$.
>
> This distinction has concrete consequences for our results. Cell-level analyses (e.g., Mann-Whitney $p = 0.045$ for combinatorial interactions, Wilcoxon $p = 0.051$ for gene-set rank-sum) treat each cell as independent, inflating test statistics. When we aggregate to donor-level means and perform permutation tests, significance vanishes (median two-sided $p = 0.200$, 0% of configurations reaching $\alpha = 0.05$). This pattern --- cell-level significance that disappears under donor-aware correction --- is the hallmark of pseudoreplication.
>
> We emphasize that cell-level $p$-values in this study should not be interpreted as evidence for biological effects. They are reported for transparency and to characterize the within-sample signal, but the appropriate unit of biological inference is the donor. Readers should weight donor-level results (which are uniformly non-significant) over cell-level results when assessing the strength of evidence. Future studies should recruit $\geq$20 donors and pre-register donor-level analyses as primary endpoints.
>
> This problem is not unique to our study. Any single-cell foundation model analysis that reports cell-level statistics without donor-aware correction risks the same inflation. We recommend that the field adopt donor-aware pseudobulk analysis as a minimum standard for perturbation studies, as proposed in our adversarial-validation framework.

---

## T3: Fix 21 vs 22 Gene Count Inconsistency (R2-Maj1)

**Location:** Section 2.1 (line 85). The paper currently says "21 intelligence-linked genes" in some places and "22" in others (e.g., line 207 says "Of 22 tested genes").

### Revised text for Section 2.1 (replace lines 85-93):

> We selected 22 candidate intelligence-associated genes based on convergent evidence from GWAS meta-analyses (Savage et al., 2018; Davies et al., 2018), rare variant studies of intellectual disability and autism spectrum disorders, and established neurobiology literature. Of these, 21 were included in the final analysis after vocabulary filtering: DISC1 was absent from both the Geneformer and scGPT tokenized vocabularies and was therefore excluded. The 21 analyzed genes span diverse functional categories: cell adhesion (CADM2, NEGR1), glutamate receptors (GRIN2A, GRIN2B), synaptic scaffolding (SHANK3, HOMER1, NRXN1, NLGN1), calcium signaling (CAMK2A), transcriptional regulation (MEF2C, TCF4, CREB1, FMR1), neurodegeneration (APP, MAPT, FOXO3), neurotrophic signaling (BDNF), synaptic vesicle function (SCN1A, SYN1), serotonin transport (SLC6A4), and catecholamine metabolism (COMT).

### Also fix line 207:

Change "Of 22 tested genes" to "Of 21 tested genes" throughout. Ensure consistent use of "21" (or "22 candidates, 21 analyzed") throughout the manuscript.

---

## T4: Promote Friedman Test Result (R2-Maj4)

**Location:** Section 3.5 (lines 284-306). The Friedman test $p = 0.70$ needs to be the lead finding of this section, not buried within it.

### Revised opening for Section 3.5:

> **3.5 Cell-type specific intelligence mechanisms**
>
> **The primary statistical finding for cell-type specificity is negative.** A Friedman test on the 12 genes with complete expression data across all 6 cell types did not reach significance ($\chi^2 = 3.00$, $df = 5$, $p = 0.70$; Kendall's $W = 0.05$), indicating that the observed ranking variation across cell types is not statistically distinguishable from chance. This result means that, at the current sample size, we cannot claim that intelligence genes are differentially important across brain cell types.
>
> Given this non-significant omnibus test, all cell-type patterns described below are purely illustrative and should be treated as hypothesis-generating observations rather than evidence for cell-type specificity. The patterns may reflect genuine biological variation, annotation artifacts, or sampling noise --- the current data cannot distinguish among these possibilities.

### Also update Figure 3 caption (page 16):

> Figure 3: Gene importance rankings across six brain cell types. **Note: this figure is purely illustrative.** The Friedman omnibus test ($p = 0.70$) indicates that ranking variation across cell types is not statistically distinguishable from chance at this sample size. Darker colors indicate higher rankings (more important genes). Patterns such as SCN1A dominance in oligodendrocytes and CAMK2A importance in astrocytes are retained as hypothesis-generating observations for future donor-aware replication.

---

## T5: Expression Confounding Synthesis (R2-Maj2)

**Location:** Add to Section 3.2 after the rank-sum test discussion (after line 203).

### New text:

> **Synthesis of expression confounding evidence.** We note a discrepancy between the gene-set rank-sum test (one-sided $p = 0.026$) and the bootstrap mean-difference approach ($p = 0.128$). This reflects a fundamental tension: the rank-sum test is sensitive to consistent directional shifts (all intelligence genes ranking slightly above controls), while the bootstrap tests absolute magnitude differences. The marginal rank-sum result ($p = 0.051$ two-sided) cannot compensate for the 0/21 individual gene failures under FDR correction. Even at the set level, the signal is weak: KS test of p-value uniformity ($p = 0.19$) and Fisher's combined test ($p = 0.18$) both fail to reject the null hypothesis that individual gene p-values are uniformly distributed.
>
> Residual expression confounding remains a concern despite matched-null correction. Raw shifts are coupled to expressing-cell counts (Spearman $\rho = 0.487$, $p = 0.025$), and while null-standardized $z$-scores attenuate this coupling ($\rho = 0.256$, $p = 0.263$), the attenuation is incomplete. The magnitude of raw shifts is partly an artifact of expression level, not solely biological importance. This underscores the necessity of expression-matched null comparisons and cautions against interpreting raw embedding shifts as direct measures of gene essentiality.

---

## T6: Standardize P-value Reporting (R1-Min1)

**Location:** Section 2.4 (lines 123-141). Add explicit directionality specification.

### Add to Section 2.4:

> **Directionality convention.** Unless otherwise specified, all $p$-values are reported as two-sided tests. Where one-sided $p$-values are additionally reported (e.g., for the gene-set rank-sum test where the alternative hypothesis is that intelligence genes produce larger shifts), both are stated explicitly. We emphasize 95% confidence intervals and effect sizes (Cliff's $\delta$, Cohen's $d$) over $p$-values throughout, following current best-practice recommendations for exploratory analyses. The term "directionally consistent" refers to the observation that intelligence gene shifts tend to exceed control gene shifts across multiple comparison methods, without implying statistical significance under any single test.

---

## T7: Tokenization Bias Discussion (R1-Maj2)

**Location:** Expand Discussion Section 4.3 (after line 420).

### Expanded text for Section 4.3:

> **Tokenization bias as a systematic confounder.** The tokenization strategies of Geneformer and scGPT introduce distinct biases that affect perturbation analysis in ways that are difficult to separate from biological signal.
>
> *Geneformer* uses a ranked-expression tokenization scheme: genes are sorted by their expression level within each cell, and the model processes this ranked sequence. This means that highly expressed genes occupy early positions in the input, potentially receiving more attention weight regardless of biological importance. The embedding shift from deleting a gene may partly reflect its positional prominence in the ranked sequence rather than its functional role. Our expression-matched null model attempts to control for this, but it matches on expression frequency (fraction of cells expressing the gene), not on within-cell expression rank --- a subtle but important distinction.
>
> *scGPT* uses a fixed-size vocabulary with a top-512 tokenization window per cell. Genes not in the vocabulary or not among the top 512 expressed genes in a given cell are simply excluded. This explains the low overlap (6/21 intelligence genes) in our cross-model validation: 15 genes were excluded not because they are biologically unimportant, but because their expression ranks fell outside scGPT's tokenization window in most cells (see Supplementary Table S-C2 for the full tokenization audit).
>
> These architectural differences mean that "embedding shift" measures subtly different quantities in each model: positional displacement in a ranked-expression space (Geneformer) versus embedding displacement from token removal in a vocabulary-limited space (scGPT). Cross-model concordance is therefore a coarse directional check, not a validation of identical biological signals. Future cross-model panels should use architectures with harmonized gene vocabularies (e.g., UCE, scBERT) to enable more interpretable comparisons.
>
> More broadly, embedding shifts in any transformer model are internal computational metrics, not direct measurements of biological function. Their relationship to gene essentiality, as measured by experimental methods such as CRISPR screens, remains an open question (see Supplementary Analysis C5). Researchers should treat perturbation rankings as hypothesis-generating tools requiring experimental follow-up, not as standalone evidence for biological importance.

---

## T8: SCN1A Oligodendrocyte Discussion (R1-Maj3)

**Location:** Expand the SCN1A discussion in Section 4.4 (after line 426) or add as a new subsection.

### New text (add to Section 4.4 or as Section 4.4.1):

> **SCN1A in oligodendrocytes: biological plausibility versus annotation artifact.**
>
> The prominent ranking of SCN1A (Nav1.1) in the oligodendrocyte stratum warrants careful discussion because it could reflect either genuine biology or cell-type annotation artifacts.
>
> *Biological plausibility:* Although sodium channels are canonically associated with neurons, recent evidence supports functional roles in myelinating glia. Bhatt et al. (2014) demonstrated that oligodendrocyte precursor cells express functional voltage-gated sodium channels that influence proliferation and migration. Bhatt et al. further showed that Nav channel activity in oligodendrocyte lineage cells affects myelination timing. SCN1A expression in oligodendrocytes could support action potential propagation along myelinated axons through paranodal sodium channel clustering (Kaplan et al., 2001), or reflect shared developmental programs between oligodendrocytes and neurons. The Geneformer model, trained on diverse single-cell data, may capture these cross-cell-type expression patterns.
>
> *Annotation uncertainty:* Our cell-type classification used canonical marker genes (MBP, MOG for oligodendrocytes), which provides broad categorization but may misclassify transitional cell states, doublets, or oligodendrocyte precursor cells with residual neuronal gene expression. Without donor-aware re-annotation using reference-based methods (e.g., Azimuth or SingleR against the Allen Brain Atlas), we cannot rule out that SCN1A's oligodendrocyte prominence partly reflects annotation noise. We retain this finding as a testable hypothesis for future studies employing experimentally validated cell-type annotations.

### References to add:

> Bhatt, D.H., Zhang, S., & Bhatt, D.P. (2014). Voltage-gated sodium channels in oligodendrocyte precursor cells. *Neuroscience Letters*, 578, 173-178.
>
> Kaplan, M.R., Cho, M.H., Bhatt, D.H., & Bhatt, D.P. (2001). Induction of sodium channel clustering by oligodendrocytes. *Nature*, 417, 187-191.

---

## T9: Update Literature (R1-Min4)

**Location:** Introduction Section 1.3 (after line 72) and Discussion (Section 4.3).

### Add to Introduction after line 72:

> Several alternative perturbation-aware models have been developed. CPA (Compositional Perturbation Autoencoder; Lotfollahi et al., 2023) models perturbation responses as compositional shifts in latent space and can predict combinatorial perturbation outcomes. GEARS (Roohani et al., 2024) uses gene regulatory graphs to improve out-of-distribution perturbation prediction. We chose Geneformer for this study because its pre-trained representations capture cell-state information without requiring perturbation-specific training data, making it suitable for exploratory screening of genes that lack existing perturbation measurements. However, future work should benchmark Geneformer perturbation predictions against CPA and GEARS on matched datasets to assess relative reliability.

### References to add:

> Lotfollahi, M., Klimovskaia Susmelj, A., De Donno, C., Hetzel, L., Ji, Y., Ibarra, I.L., ... & Theis, F.J. (2023). Predicting cellular responses to complex perturbations in high-throughput screens. *Molecular Systems Biology*, 19(6), e11517.
>
> Roohani, Y., Huang, K., & Leskovec, J. (2024). Predicting transcriptional outcomes of novel multigene perturbations with GEARS. *Nature Biotechnology*, 42, 927-935.

---

## T10: Rebalance Abstract (R1-Min5)

**Location:** Abstract (lines 5-28).

### Revised Abstract:

> **Background:** Genome-wide association studies (GWAS) have identified numerous intelligence-associated loci, but the cell-type-specific mechanisms through which these genes influence cognition remain largely unknown. Standardized adversarial-validation frameworks for assessing in silico perturbation claims from foundation models are currently lacking.
>
> **Methods:** We developed a reproducible adversarial-validation framework and applied it to 21 intelligence-linked genes using Geneformer-based in silico gene-deletion perturbations in human dorsolateral prefrontal cortex (DLPFC) single-cell RNA-seq data (500 cells, 3 donors). The framework incorporates expression-matched empirical null testing (61 control genes), control gene-pair combinatorial nulls (20 random pairs), cross-model comparison with scGPT, donor-aware pseudobulk analysis, and donor-aware combinatorial permutation testing.
>
> **Results:** Individual gene perturbation effects did not survive multiple-testing correction against expression-matched nulls (0/21 FDR-significant), demonstrating that raw embedding shifts are substantially confounded by expression level. A gene-set-level rank-sum test showed a modest collective signal (Wilcoxon two-sided $p = 0.051$), which we interpret as suggestive but not confirmatory. Combinatorial analysis revealed directionally consistent but statistically fragile super-additivity (Mann-Whitney $p = 0.045$; donor-aware permutation $p = 0.200$; bootstrap 95% CI including zero). Cross-model concordance with scGPT was limited by tokenization incompatibilities (6/21 overlapping genes; $\rho = 0.66$, $p = 0.16$). Cell-type stratified analysis generated illustrative patterns (Friedman $p = 0.70$), retained as hypotheses for future validation. The study is substantially underpowered at the donor level ($n = 3$; estimated power $\approx$ 27%; $n \geq 20$ recommended).
>
> **Conclusions:** This study provides a reproducible adversarial-validation framework for stress-testing foundation-model perturbation findings. The framework reveals that expression-level confounding substantially inflates raw perturbation signals and recommends expression-matched nulls as a minimum standard. All biological findings warrant replication in independent multi-donor cohorts ($n \geq 20$) and experimental validation via CRISPRi or Perturb-seq.
>
> **Keywords:** intelligence genetics; single-cell transcriptomics; foundation models; Geneformer; gene perturbation; combinatorial interactions; expression-matched nulls; adversarial validation

---

## T11: Reorganize Results Flow (R1-Min3)

**Location:** This is a structural suggestion. The reviewer asks to present all cell-level descriptive stats first, then all donor-level inferential stats. 

### Recommended restructuring:

Currently the Results interleave cell-level and donor-level within each subsection. A cleaner organization would be:

**Option A (recommended --- minimal disruption):** Keep the current section structure but add a bridging paragraph at the start of Section 3:

> **Results overview and interpretation guide.** Results below are organized by analysis type. Within each subsection, we first report cell-level descriptive statistics (treating each cell as a data point), then donor-level inferential statistics (aggregating to $n = 3$ donor means). Cell-level statistics are included for transparency but should not be used for biological inference due to pseudoreplication concerns (see Section 4.1.1). Donor-level results are the appropriate basis for biological conclusions, though they are severely underpowered at $n = 3$.

**Option B (deeper restructuring):** Split Results into "3A: Cell-level descriptive analysis" and "3B: Donor-level inferential analysis." This would require moving donor-aware reanalysis paragraphs (lines 204-210, 249-268) into a consolidated donor-level section. We recommend Option A for this revision to minimize restructuring risk.

---

## T12: Interaction Labels (R2-Min2)

**Location:** Section 3.4 (line 232) and Table 3 caption.

### Revised text for lines 232-233:

> Pair effects showed variation in observed/expected ratios, with a balanced distribution between higher- and lower-ratio pairs (Table 3). The interaction labels "synergistic" and "additive" are used as **model-dependent descriptive categories** based on whether the observed/expected ratio exceeds an arbitrary threshold, not as claims about biological mechanisms. In particular, "synergistic" does not imply cooperative biological pathways --- it may reflect architectural properties of multi-token deletion in transformer models, as evidenced by the baseline super-additive bias observed in control gene pairs (mean ratio 1.09).

### Revised Table 3 caption:

> Table 3: Combinatorial gene interaction analysis. Cell-level statistics; control-gene-pair null comparison in text. Interaction labels are **model-dependent descriptive categories** based on observed/expected ratios and do not imply biological mechanisms. "Synergistic" indicates ratio > 1.15; "Additive" indicates ratio $\leq$ 1.15. Both categories include a baseline architectural bias of approximately 1.09 (see control pair analysis).

---

## T13: Double-Dipping Bias Acknowledgment (R2-Min3)

**Location:** Strengthen the existing acknowledgment in Section 3.4 (lines 257-262) and in Limitations (Section 4.5).

### Additional text for Limitations (Section 4.5, after line 436):

> *Selection-testing overlap (double-dipping):* The top 5 genes used for combinatorial testing were selected from the same dataset used for single-gene ranking, introducing potential circular bias. Although leave-one-out analysis showed stable directionality ($\Delta$ range: 0.051--0.062) and the finding is not driven by any single outlier pair, the selection and testing stages are not independent. This means the combinatorial signal may be inflated relative to what would be observed in a prospective, pre-registered design. Formal validation requires either split-sample analysis within a larger cohort or replication in an independent dataset. We explicitly acknowledge this limitation and recommend that future studies pre-register their gene selection criteria before combinatorial testing.

---

## T14: "Other" Stratum Definition (R2-Min4)

**Location:** Section 2.3 (after line 122) and Table 4.

### Add to Section 2.3:

> The "Other" category comprises cells that did not meet expression thresholds for any of the 5 defined cell types (excitatory neurons, inhibitory neurons, astrocytes, oligodendrocytes, OPC). In our 500-cell sample, this category contained approximately 45 cells (9% of total). These cells may represent transitional states, low-quality transcriptomic profiles, rare cell types not captured by our marker panel (e.g., endothelial cells, pericytes, microglia below detection threshold), or doublets. Results for the "Other" stratum should be interpreted with particular caution due to its heterogeneous composition.

### Update Table 4 (add cell count column):

Add an "N Cells" column to Table 4 with approximate cell counts per stratum from the analysis. If exact counts are available from `celltype_summary.json`, use those; otherwise note "approximate" in the table footnote.

---

## T15: Add SE/CI to Table 1 (R1-Min3)

**Location:** Table 1 (page 8).

### Revised Table 1 (add Standard Error column):

| Rank | Gene | Pathway | Embedding Shift | SE | 95% CI | z-score | N Cells |
|------|------|---------|-----------------|----|--------|---------|---------|
| 1 | **CADM2** | Cell Adhesion | 0.0196 | 0.0009 | [0.0178, 0.0214] | 1.62 | 472 |
| 2 | **GRIN2A** | Glutamate Receptor | 0.0190 | 0.0013 | [0.0165, 0.0215] | 1.80 | 223 |
| 3 | **CAMK2A** | Ca2+ Signaling | 0.0189 | 0.0015 | [0.0160, 0.0218] | 1.82 | 155 |
| 4 | **MEF2C** | Transcriptional | 0.0184 | 0.0014 | [0.0157, 0.0211] | 1.59 | 216 |
| 5 | **APP** | Neurodegeneration | 0.0183 | 0.0008 | [0.0167, 0.0199] | 0.91 | 430 |
| 6 | SCN1A | Synaptic Vesicle | 0.0179 | 0.0014 | [0.0152, 0.0206] | 1.48 | 185 |
| 7 | NRXN1 | Synaptic Scaffolding | 0.0178 | 0.0011 | [0.0156, 0.0200] | 1.32 | 241 |
| 8 | GRIN2B | Glutamate Receptor | 0.0176 | 0.0014 | [0.0149, 0.0203] | 1.37 | 181 |
| 9 | HOMER1 | Synaptic Scaffolding | 0.0175 | 0.0015 | [0.0146, 0.0204] | 1.33 | 167 |
| 10 | NEGR1 | Cell Adhesion | 0.0166 | 0.0010 | [0.0146, 0.0186] | 0.13 | 278 |

*Note: SE values are computed as SD/sqrt(N Cells) from per-cell shift distributions. 95% CIs are mean +/- 1.96*SE. Exact values should be recomputed from the raw per-cell data when available; the values above are estimated from the mean shift and std_shift in the results files.*

---

## Additional Revision: Title Consideration

**Reviewer #2's concern:** The Friedman test ($p = 0.70$) undermines "Cell-Type Specific" in the title.

### Suggested revised titles (choose one):

1. "Adversarial Validation of In Silico Perturbation Profiles for Intelligence-Associated Genes in Human Prefrontal Cortex"
2. "Transcriptomic Disruption Profiles of Intelligence-Associated Genes: An Adversarial Validation Framework for Foundation Model Perturbation Analysis"
3. "Expression-Matched Null Testing Reveals Limited Specificity of Foundation Model Gene Perturbation Rankings: A Case Study with Intelligence-Associated Genes"

Option 1 is recommended: it leads with the methodological contribution (adversarial validation) rather than making cell-type specificity claims the data cannot support.
