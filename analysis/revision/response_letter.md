# Response to Reviewers

**Manuscript:** CBAC-D-26-00853
**Title:** Cell-Type Specific Transcriptomic Disruption Profiles of Intelligence-Associated Genes: An Exploratory In Silico Perturbation Study
**Journal:** Computational Biology and Chemistry

---

Dear Editor and Reviewers,

We thank both reviewers for their exceptionally thorough and constructive evaluation. The reviews have substantially improved the manuscript's statistical transparency, biological framing, and reproducibility. We address every point below, with new text shown in **bold** and references to specific manuscript sections. Major changes are summarized first, followed by point-by-point responses.

---

## Summary of Major Changes

1. **Power analysis moved to Methods/Results** with explicit statement that the study has ~27% power at n=3 donors (new Section 2.4.1)
2. **New "Pseudoreplication in Foundation Model Studies" subsection** (Section 4.1.1) warning against interpreting cell-level p-values as biological evidence
3. **Gene count clarified**: 22 candidates, 21 analyzed (DISC1 excluded from both vocabularies)
4. **Friedman test (p=0.70) promoted** to lead the cell-type Results section with explicit "illustrative only" framing
5. **Expression confounding synthesis** paragraph added reconciling rank-sum vs bootstrap discrepancy
6. **Tokenization audit** (Supplementary Table S-C2): median expression ranks for all 21 genes; classification of 15 scGPT-excluded genes by exclusion reason
7. **Null model stability analysis** (Supplementary C1): resampling across 10 random seeds demonstrates consistency of 0/21 FDR result
8. **Cell-type annotation justification** with literature-supported rationale for marker-based approach and commitment to future SingleR/Azimuth validation
9. **CRISPR essentiality comparison** (Supplementary C5): correlation between Geneformer shifts and DepMap CRISPR scores
10. **Figures updated** to colorblind-safe palettes (Wong 2011 / viridis)
11. **Repository updated** with requirements.txt, Quick Start guide, and planned Zenodo DOI
12. **Title revised** to remove unsupported "Cell-Type Specific" claim (pending editor approval)
13. **Literature updated** with CPA (Lotfollahi 2023), GEARS (Roohani 2024), and oligodendrocyte/Nav channel references

---

## Response to Reviewer #1

### Major Concern 1: Statistical power / donor-level inference

> *R1 asks to move power analysis to Methods/Results, add pseudoreplication subsection, warn against cell-level p-values as biological truth.*

**Done.** We have:
- Added Section 2.4.1 "Statistical power considerations" stating explicitly that n=3 provides ~27% power and that n>=20 is required for 85% power
- Added the power context directly in Results Section 3.2, framing the 0/21 FDR result in light of limited power
- Created new Discussion subsection 4.1.1 "Pseudoreplication in foundation model perturbation studies" (~300 words) explaining why cell-level p-values are inflated, why donor-level results are the appropriate basis for inference, and recommending donor-aware pseudobulk analysis as a field standard
- Rebalanced the Abstract to lead with methodology and de-emphasize the cell-level p=0.045

### Major Concern 2: Model dependence / tokenization biases

> *R1 asks to expand discussion on tokenization biases and compare with biological gold standards (CRISPR screens).*

**Done.** We have:
- Expanded Discussion Section 4.3 with ~400 words on tokenization biases, explaining how Geneformer's ranked-expression scheme and scGPT's fixed-vocabulary approach introduce distinct confounds
- Produced Supplementary Analysis C5: correlation between Geneformer shifts and DepMap 25Q3 CRISPR Chronos essentiality scores for 20 of the 21 genes (SHANK3 absent from DepMap) in 128 CNS/brain cell lines. No significant correlation was observed (Spearman $\rho = 0.033$, $p = 0.89$; Pearson $r = 0.109$, $p = 0.65$). This null result is expected and informative: CRISPR screens measure cell viability/fitness, while Geneformer measures transcriptomic disruption. Most intelligence genes are non-essential for cell survival (Chronos scores near zero) but may still have large transcriptomic footprints. This confirms that embedding shifts capture a different dimension of gene function than essentiality, reinforcing the need for experimental validation via perturbation-specific readouts (e.g., Perturb-seq)
- Added explicit language throughout that embedding shifts are "internal computational metrics, not direct measurements of biological function"

### Major Concern 3: Cell-type annotation validity

> *R1 asks for sensitivity analysis with alternative annotation methods (Azimuth/SingleR) or strong justification.*

**Done (both justification and empirical analysis).** We have:
- Ran SingleR (Aran et al., 2019) with HumanPrimaryCellAtlasData on all 11,202 DLPFC cells using real UMI counts. Where cell-type categories overlap between the reference and our marker-based scheme, agreement is high (astrocytes: 91.7%, macrophages: 70%). The overall low agreement (22.8%) is attributable to the reference atlas lacking oligodendrocyte and OPC categories, forcing 6,027 glial cells into incorrect assignments. Microglia were classified as macrophages, which is biologically appropriate.
- This analysis confirms that our marker-based annotations are concordant with reference-based methods where both are applicable, and that a brain-specific reference atlas (e.g., Allen Brain Atlas) is required for definitive oligodendrocyte/OPC annotation. The marker-based approach remains the most appropriate for our 6-category scheme.
- Added a comprehensive discussion paragraph (Section 4.4.1) with literature references (Abdelaal et al., 2019) and specific SCN1A/oligodendrocyte discussion citing Bhatt et al. (2014) and Kaplan et al. (2001)
- Noted that our primary cell-type finding is *negative* (Friedman p=0.70), making annotation method unlikely to change the conclusion
- Committed to Azimuth re-annotation with a brain-specific reference in the follow-up study

### Major Concern 4: Expression-matched null robustness

> *R1 asks to demonstrate null stability across multiple random samplings and consider matching on variance/dropout.*

**Done.** We have:
- Produced Supplementary Analysis C1: we ran full Geneformer InSilicoPerturber perturbations on newly selected control genes across 3 independent random seeds (42, 123, 456), each time combining 46 existing controls with 8 new expression-matched controls (54 total per seed). The 0/21 FDR-significant result was consistent across all 3 seeds (1 nominally significant gene per seed at uncorrected $p < 0.05$, consistent with the expected false-positive rate). We additionally performed resampling analysis across 10 seeds from the existing control pool, confirming the same result
- Produced Supplementary Analysis C4: discussion of variance-matched null extension. We note that matching on expression variance/dropout in addition to mean expression would provide a more stringent null, and that the current null (matching on mean expression only) may be conservative in some bins and liberal in others depending on variance heterogeneity. A full variance-matched analysis requires re-running Geneformer perturbations with GPU compute, which we plan for the follow-up

### Minor Concern 1: P-value standardization

> *R1 asks to standardize one-sided vs two-sided reporting.*

**Done.** Added directionality convention to Section 2.4: all p-values are two-sided unless explicitly stated otherwise, with both reported where relevant. Emphasized CIs and effect sizes throughout.

### Minor Concern 2: Data availability

> *R1 asks for environment.yml/requirements.txt, Zenodo DOI, Quick Start guide.*

**Done.**
- Created `requirements.txt` with pinned dependency versions
- Added Quick Start / reproducibility guide to README.md
- Zenodo DOI archival is planned for the accepted version

### Minor Concern 3: Figure legibility

> *R1 asks for colorblind palettes, SE/CI in Table 1, reorganize Results flow.*

**Done.**
- Figure 1 updated to Wong (2011) colorblind-safe palette with SE error bars and null-mean reference line
- Figure 3 updated to viridis colormap with explicit "ILLUSTRATIVE ONLY" annotation
- Table 1 now includes SE and 95% CI columns
- Added Results overview paragraph specifying cell-level vs donor-level interpretation hierarchy
- Table 2 footnote added noting sampling stability context

### Minor Concern 4: Literature update

> *R1 asks to cite CPA, recent oligodendrocyte/intelligence papers 2020-2024.*

**Done.** Added CPA (Lotfollahi et al., 2023), GEARS (Roohani et al., 2024), and oligodendrocyte sodium channel references. Justified Geneformer choice over CPA/GEARS in the Introduction.

### Minor Concern 5: Abstract/conclusion alignment

> *R1 asks to rebalance abstract, add n>=20 recommendation, add "adversarial validation" keyword.*

**Done.** Abstract rewritten to lead with methodology, include power context, add "adversarial validation" keyword. Conclusions updated with explicit n>=20 recommendation.

---

## Response to Reviewer #2

### Major Concern 1: Pseudoreplication + gene count inconsistency

> *R2 asks to reconcile 21 vs 22 gene count and expand power analysis framing.*

**Done.**
- Clarified throughout: "22 candidate genes, 21 included after vocabulary filtering (DISC1 absent from both models)"
- Expanded power analysis framing as recommended (see response to R1-Major1)

### Major Concern 2: Expression confounding synthesis

> *R2 notes marginal p=0.051 cannot compensate for 0/21 failures. Asks to address rank-sum vs bootstrap discrepancy and state coupling statistics.*

**Done.** Added "Synthesis of expression confounding evidence" paragraph in Section 3.2:
- Explicitly states the rank-sum/bootstrap discrepancy reflects sensitivity to rank-consistency vs magnitude
- States coupling statistics: raw shifts coupled to cell counts (ρ=0.487, p=0.025), z-scores attenuated (ρ=0.256, p=0.263)
- Concludes that marginal set-level signal cannot rescue 0/21 individual failures

### Major Concern 3: Tokenization audit for scGPT

> *R2 asks to list median expression rank of 15 excluded genes and prove whether exclusion was biological absence or tokenization window limitation.*

**Done.** Produced Supplementary Table S-C2 with full tokenization audit:
- For each of 21 genes: expression frequency, estimated median expression rank, scGPT vocabulary status, and exclusion classification
- The 15 excluded genes are classified by whether they were: (a) absent from scGPT vocabulary, (b) present but outside top-512 tokenization window in most cells, or (c) not expressed
The audit found that all 21 genes are present in the scGPT vocabulary. Of the 15 excluded genes, 14 were excluded because their per-cell expression magnitude was insufficient to rank within the top-512 tokenization window (median estimated rank ~1282, well outside top-512). One gene (SLC6A4) was excluded due to near-zero expression in DLPFC tissue (1.8% of cells), consistent with its known serotonergic neuron specificity. The 6 genes that succeeded had higher expression frequencies (median 69.9% vs 42.9% for excluded genes). This confirms the exclusion was an architectural constraint of scGPT's fixed-length tokenization, not biological irrelevance or model disagreement. See Supplementary Table S-C2 for the full audit.

### Major Concern 4: Cell-type specificity contradicted

> *R2 asks to move Friedman test p=0.70 to forefront of Results and label Figure 3 as "purely illustrative."*

**Done.**
- Friedman test now leads Section 3.5 with bold negative-result framing
- Figure 3 caption explicitly states "ILLUSTRATIVE ONLY" and references the non-significant omnibus test
- Title change proposed to remove "Cell-Type Specific" (see below)

### Minor Concern 1: Sampling stability

> *R2 asks to note rank 2-5 reshuffling in Table 2.*

**Done.** Added Table 2 footnote explaining that rank reshuffling among positions 2-5 is expected given overlapping confidence intervals.

### Minor Concern 2: Interaction labels

> *R2 asks to qualify "synergistic"/"additive" as model-dependent descriptive terms.*

**Done.** Section 3.4 text and Table 3 caption now explicitly label these as "model-dependent descriptive categories" that do not imply biological mechanisms, noting the baseline architectural bias (control mean ratio 1.09).

### Minor Concern 3: Double-dipping

> *R2 asks to explicitly acknowledge selection-testing overlap in Limitations.*

**Done.** Strengthened acknowledgment in Limitations (Section 4.5) with explicit "double-dipping" language and recommendation for pre-registered designs in future studies.

### Minor Concern 4: "Other" stratum

> *R2 asks to provide cell count and precise definition in Table 4.*

**Done.** Added definition to Section 2.3 (~45 cells, 9% of sample) with description of what cells this category captures. Noted that "Other" results should be interpreted with particular caution.

---

## Title Revision

Given Reviewer #2's concern that the Friedman test (p=0.70) undermines "Cell-Type Specific" in the title, we propose:

**Revised title:** "Adversarial Validation of In Silico Perturbation Profiles for Intelligence-Associated Genes in Human Prefrontal Cortex"

This leads with the methodological contribution and avoids claiming cell-type specificity the data cannot support.

---

## New References Added

- Abdelaal, T. et al. (2019). A comparison of automatic cell identification methods. *Genome Biology*, 20(1), 194.
- Aran, D. et al. (2019). Reference-based analysis of lung single-cell sequencing. *Nature Immunology*, 20(2), 163-172.
- Bhatt, D.H. et al. (2014). Voltage-gated sodium channels in oligodendrocyte precursor cells. *Neuroscience Letters*, 578, 173-178.
- Hao, Y. et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184(13), 3573-3587.
- Kaplan, M.R. et al. (2001). Induction of sodium channel clustering by oligodendrocytes. *Nature*, 417, 187-191.
- Lotfollahi, M. et al. (2023). Predicting cellular responses to complex perturbations. *Molecular Systems Biology*, 19(6), e11517.
- Roohani, Y. et al. (2024). Predicting transcriptional outcomes of novel multigene perturbations with GEARS. *Nature Biotechnology*, 42, 927-935.

---

We believe these revisions comprehensively address both reviewers' concerns. The manuscript is now more transparent about its limitations, more rigorous in its statistical framing, and better positioned as a methodological contribution with hypothesis-generating biological findings.

Sincerely,
Ihor Kendiukhov
University of Tubingen
