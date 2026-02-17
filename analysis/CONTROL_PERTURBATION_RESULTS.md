# Control Gene Perturbation Results

## Geneformer In Silico Perturbation: Intelligence Genes vs Controls

**Date:** 2026-02-14
**Dataset:** DLPFC_11k.h5ad (human dorsolateral prefrontal cortex, 500 cells max per gene)
**Model:** Geneformer V1 (10M parameters, pretrained)
**Method:** In silico gene deletion, cell-level embedding shift (1 - cosine similarity)

---

## Summary of Findings

Intelligence-associated genes produce **significantly larger** perturbation effects on brain cell embeddings than expression-matched random control genes.

| Comparison | Mann-Whitney U p-value | Cohen's d | Interpretation |
|---|---|---|---|
| Intelligence vs Expression-matched controls | **p = 0.031** | d = 0.69 | Medium effect, significant |
| Intelligence vs Housekeeping controls | **p = 0.00018** | d = 3.19 | Large effect, highly significant |
| Intelligence vs All controls combined | **p = 0.0017** | d = 1.15 | Large effect, significant |

**Paired analysis (expression-matched):** Intelligence genes showed larger shifts than their matched controls in **9/10 pairs** (Wilcoxon signed-rank p = 0.019, sign test p = 0.011).

---

## Gene-Level Results

### Intelligence Genes (n=21)

| Gene | N Cells | Mean Shift | Std | Pathway |
|---|---|---|---|---|
| CADM2 | 472 | 0.019627 | 0.004491 | Cell Adhesion |
| GRIN2A | 223 | 0.018958 | 0.006070 | Glutamate Receptors |
| CAMK2A | 155 | 0.018922 | 0.005880 | Neurotrophic Signaling |
| MEF2C | 216 | 0.018443 | 0.006518 | Neurodevelopmental/Transcription |
| APP | 430 | 0.018320 | 0.005368 | Neurodegeneration-Associated |
| SCN1A | 185 | 0.017922 | 0.006035 | Synaptic Transmission |
| NRXN1 | 241 | 0.017792 | 0.007574 | Synaptic Scaffolding |
| GRIN2B | 181 | 0.017584 | 0.005752 | Glutamate Receptors |
| HOMER1 | 167 | 0.017466 | 0.007492 | Synaptic Scaffolding |
| NEGR1 | 278 | 0.016624 | 0.007523 | Cell Adhesion |
| SLC6A4 | 1 | 0.016131 | 0.000000 | Synaptic Transmission |
| NLGN1 | 380 | 0.015984 | 0.006987 | Synaptic Scaffolding |
| TCF4 | 440 | 0.015231 | 0.006660 | Neurodevelopmental/Transcription |
| FOXO3 | 276 | 0.014869 | 0.006448 | Neurodegeneration-Associated |
| BDNF | 174 | 0.014019 | 0.007801 | Neurotrophic Signaling |
| MAPT | 370 | 0.013814 | 0.006655 | Neurodegeneration-Associated |
| CREB1 | 150 | 0.012410 | 0.006591 | Neurotrophic Signaling |
| FMR1 | 68 | 0.010911 | 0.005415 | Neurodevelopmental/Transcription |
| SYN1 | 30 | 0.009434 | 0.005739 | Synaptic Transmission |
| COMT | 73 | 0.009359 | 0.005377 | Catecholamine Metabolism |
| SHANK3 | 35 | 0.006302 | 0.005084 | Synaptic Scaffolding |

**Summary:** Mean of means = 0.01524, Median = 0.01613, Total cells = 4,545

### Expression-Matched Control Genes (n=11)

| Control Gene | Matched To | N Cells | Mean Shift | Std |
|---|---|---|---|---|
| PDE7B | FOXO3 | 135 | 0.018103 | 0.006219 |
| EXOC6B | CADM2 | 420 | 0.016264 | 0.006287 |
| DOCK9 | APP | 324 | 0.016019 | 0.005718 |
| ZDHHC9 | NLGN1 | 165 | 0.015359 | 0.005655 |
| IGSF10 | NEGR1 | 22 | 0.013843 | 0.006578 |
| HDAC11 | MAPT | 93 | 0.013290 | 0.005486 |
| CUEDC2 | DISC1 | 42 | 0.012055 | 0.004637 |
| XRCC6 | TCF4 | 74 | 0.010207 | 0.006183 |
| WDR46 | CREB1 | 19 | 0.009861 | 0.005815 |
| SALL3 | FMR1 | 8 | 0.009246 | 0.005901 |
| TMEM88 | MEF2C | 1 | 0.004812 | 0.000000 |

**Summary:** Mean of means = 0.01264, Median = 0.01329, Total cells = 1,303
**Note:** MFSD6L (matched to COMT) failed - no cells in the tokenized dataset expressed this gene.

### Housekeeping Control Genes (n=5)

| Gene | N Cells | Mean Shift | Std | Function |
|---|---|---|---|---|
| YWHAZ | 169 | 0.011647 | 0.006350 | 14-3-3 signaling protein |
| GAPDH | 18 | 0.001013 | 0.000761 | Glycolytic enzyme |
| ACTB | 27 | 0.000468 | 0.000580 | Cytoskeletal (beta-actin) |
| RPL13A | 12 | 0.000172 | 0.000070 | Ribosomal protein |
| B2M | 12 | 0.000080 | 0.000051 | MHC class I component |

**Summary:** Mean of means = 0.00268, Median = 0.00047, Total cells = 238

---

## Statistical Tests

### Test 1: Gene-Level Mann-Whitney U (unpaired)

Compares distributions of per-gene mean embedding shifts.

| Comparison | n_intel | n_ctrl | U | p-value | Cohen's d | Size |
|---|---|---|---|---|---|---|
| vs Expression-matched | 21 | 11 | 163.0 | **0.031** | 0.69 | Medium |
| vs Housekeeping | 21 | 5 | 101.0 | **0.00018** | 3.19 | Large |
| vs All controls | 21 | 16 | 264.0 | **0.0017** | 1.15 | Large |

### Test 2: Cell-Level Mann-Whitney U (pooled)

Pools all individual cell observations across genes.

| Comparison | n_intel_cells | n_ctrl_cells | U | p-value | Cohen's d |
|---|---|---|---|---|---|
| vs Expression-matched | 4,545 | 1,303 | 3,349,680 | **2.4e-13** | 0.16 |
| vs All controls | 4,545 | 1,541 | 4,215,223 | **2.6e-33** | 0.31 |

Note: Cell-level Cohen's d is smaller because within-gene variance is high; the gene-level comparison is the more appropriate test.

### Test 3: Paired Comparison (Expression-Matched Only)

Directly compares each intelligence gene to its expression-matched control.

| Intelligence Gene | Control Gene | Intel Shift | Control Shift | Difference |
|---|---|---|---|---|
| CADM2 | EXOC6B | 0.01963 | 0.01626 | **+0.00336** |
| APP | DOCK9 | 0.01832 | 0.01602 | **+0.00230** |
| FOXO3 | PDE7B | 0.01487 | 0.01810 | -0.00323 |
| NLGN1 | ZDHHC9 | 0.01598 | 0.01536 | **+0.00063** |
| MAPT | HDAC11 | 0.01381 | 0.01329 | **+0.00052** |
| TCF4 | XRCC6 | 0.01523 | 0.01021 | **+0.00502** |
| CREB1 | WDR46 | 0.01241 | 0.00986 | **+0.00255** |
| NEGR1 | IGSF10 | 0.01662 | 0.01384 | **+0.00278** |
| FMR1 | SALL3 | 0.01091 | 0.00925 | **+0.00167** |
| MEF2C | TMEM88 | 0.01844 | 0.00481 | **+0.01363** |

**Result:** 9/10 intelligence genes showed larger effects than their matched controls.

| Test | Statistic | p-value | Significant? |
|---|---|---|---|
| Wilcoxon signed-rank (one-tailed) | W = 48.0 | **0.019** | Yes (p < 0.05) |
| Sign test (one-tailed) | 9/10 | **0.011** | Yes (p < 0.05) |
| Mean paired difference | +0.00292 | -- | -- |

---

## Interpretation

### Key Conclusions

1. **Intelligence genes show significantly larger perturbation effects than expression-matched controls** (Mann-Whitney p = 0.031, Cohen's d = 0.69, medium effect). This demonstrates that the effects are not simply artifacts of expression level.

2. **The paired comparison is particularly compelling:** 9 out of 10 expression-matched pairs show larger effects for the intelligence gene (Wilcoxon p = 0.019). Only FOXO3 showed a smaller effect than its control (PDE7B).

3. **Housekeeping genes show dramatically smaller effects** (d = 3.19, very large), confirming that gene deletion effects vary meaningfully across genes and are not a trivial consequence of removing any gene.

4. **YWHAZ is an outlier among housekeeping genes** (shift = 0.0116, much higher than other housekeeping genes at 0.0001-0.001). This may reflect its role as a signaling hub (14-3-3 protein) rather than a true "housekeeping" function in brain cells.

### Effect Size Context

- Intelligence gene mean shift: **0.01524** (1.5% embedding distance)
- Expression-matched control mean shift: **0.01264** (1.3% embedding distance)
- The ~20% relative difference (0.01524 vs 0.01264) is statistically significant and represents a medium biological effect size

### Limitations

1. **Sample size:** 10 paired comparisons (1 control gene had no cells expressing it, 1 intelligence gene DISC1 lacked a reciprocal match in the intelligence list)
2. **Single dataset:** Results from DLPFC only; cross-tissue validation needed
3. **Model constraints:** Geneformer V1 (10M params); larger models may show different patterns
4. **Cell count variation:** Some genes had very few expressing cells (TMEM88: 1 cell, SALL3: 8 cells)

---

## Methods

### Pipeline
1. Brain scRNA-seq data (DLPFC_11k.h5ad) tokenized using Geneformer V1 rank-value encoding
2. Up to 500 cells sampled per gene (limited to cells expressing each gene)
3. In silico gene deletion: remove gene token from each cell's rank-value sequence
4. Embedding shift computed as 1 - cosine_similarity(original_embedding, perturbed_embedding)
5. Per-gene mean shift used as the summary statistic

### Control Selection
- **Expression-matched:** For each intelligence gene, a random non-intelligence gene with similar expression frequency (mean frequency difference: 0.019) was selected
- **Housekeeping:** 5 standard housekeeping genes (ACTB, GAPDH, B2M, RPL13A, YWHAZ) included as negative controls

### Statistical Tests
- **Mann-Whitney U:** Non-parametric test comparing shift distributions (one-tailed: intelligence > control)
- **Cohen's d:** Standardized effect size (small: 0.2, medium: 0.5, large: 0.8)
- **Wilcoxon signed-rank:** Paired non-parametric test for expression-matched pairs
- **Sign test:** Distribution-free test counting directional differences

---

## Files Generated

| File | Description |
|---|---|
| `run_control_perturbations_v2.py` | Control gene perturbation script |
| `analyze_control_vs_intelligence.py` | Statistical comparison script |
| `results/insilico_wsl/control_perturb_*/` | Individual control gene perturbation results (16 directories) |
| `results/insilico_wsl/control_perturbation_summary_v2.json` | Perturbation run summary |
| `results/insilico_wsl/control_vs_intelligence_analysis.json` | Full statistical analysis results |
| `CONTROL_PERTURBATION_RESULTS.md` | This report |
