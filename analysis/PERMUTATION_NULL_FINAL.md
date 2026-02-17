# Permutation Null Experiment: Intelligence Gene Perturbation Effects

**Date:** 2026-02-17 01:26:24
**Dataset:** DLPFC_11k.h5ad (human dorsolateral prefrontal cortex)
**Model:** Geneformer V1 (10M parameters, pretrained)
**Method:** In silico gene deletion, 500 cells per gene, cell-level embedding shift (1 - cosine similarity)

---

## Summary

This analysis builds an empirical null distribution by running Geneformer
delete-perturbation on 50 random expression-matched control genes (10 per
expression quintile bin), plus 11 previously-run expression-matched controls.
Each control gene was run in a fresh subprocess to avoid OOM and broken-pipe
errors. Empirical p-values are computed per expression bin, with Benjamini-
Hochberg FDR correction across all 21 intelligence genes.

- **Intelligence genes tested:** 21
- **New control genes (permutation):** 44
- **Existing controls reused:** 11
- **Total controls in null:** 54

- **Nominally significant (p < 0.05):** 1/21 genes
- **FDR < 0.10:** 0/21 genes
- **FDR < 0.05:** 0/21 genes
- **FDR < 0.01:** 0/21 genes

- **Mean intelligence gene shift:** 0.015244
- **Mean null (control) shift:** 0.014011
- **Null std:** 0.003294
- **Null median:** 0.014172
- **Null 95th percentile:** 0.018500

---

## Method

1. **Expression frequency computation:** For each gene in DLPFC_11k.h5ad,
   compute the fraction of cells expressing it (count > 0).
2. **Expression matching:** For each intelligence gene, identify all genes
   within +/- 20% expression frequency that are in Geneformer's token dictionary
   (excluding intelligence genes, existing controls, housekeeping genes).
3. **Stratified sampling:** Divide expression range into 5 equal-width bins.
   Sample 50 control genes proportionally across bins (approx 10 per bin).
4. **Perturbation:** Run Geneformer delete-perturbation on each control gene
   in a **fresh subprocess** (500 cells max, batch size 32). Each gene's result
   is saved incrementally to JSON before the next gene starts.
5. **Null augmentation:** Include 11 previously-run expression-matched controls
   (from initial paired control analysis) in the null distribution.
6. **Empirical p-value:** For each intelligence gene:
   p = (n_controls_in_same_bin_with_shift >= intel_shift + 1) / (n_controls_in_bin + 1)
7. **FDR correction:** Benjamini-Hochberg across all 21 genes.

---

## Results

| Rank | Gene | Pathway | Shift | Null Mean | Null N | z-score | p (emp) | FDR q | Sig |
|------|------|---------|-------|-----------|--------|---------|---------|-------|-----|
| 1 | **CADM2** | Cell Adhesion | 0.019627 | 0.016409 | 11 | 1.62 | 0.1667 | 0.4375 | ns |
| 2 | **GRIN2A** | Glutamate Receptors | 0.018958 | 0.014567 | 14 | 1.80 | 0.0667 | 0.4286 | ns |
| 3 | **CAMK2A** | Neurotrophic Signaling | 0.018922 | 0.013495 | 23 | 1.82 | 0.0417 | 0.4286 | ns |
| 4 | **MEF2C** | Neurodevelopmental/Transcription | 0.018443 | 0.014567 | 14 | 1.59 | 0.0667 | 0.4286 | ns |
| 5 | **APP** | Neurodegeneration-Associated | 0.018320 | 0.016336 | 9 | 0.91 | 0.3000 | 0.6300 | ns |
| 6 | **SCN1A** | Synaptic Vesicle/Transmission | 0.017922 | 0.013495 | 23 | 1.48 | 0.1250 | 0.4286 | ns |
| 7 | **NRXN1** | Synaptic Scaffolding | 0.017792 | 0.014567 | 14 | 1.32 | 0.2000 | 0.4667 | ns |
| 8 | **GRIN2B** | Glutamate Receptors | 0.017584 | 0.013495 | 23 | 1.37 | 0.1250 | 0.4286 | ns |
| 9 | **HOMER1** | Synaptic Scaffolding | 0.017466 | 0.013495 | 23 | 1.33 | 0.1250 | 0.4286 | ns |
| 10 | **NEGR1** | Cell Adhesion | 0.016624 | 0.016336 | 9 | 0.13 | 0.4000 | 0.7636 | ns |
| 11 | **SLC6A4** | Synaptic Vesicle/Transmission | 0.016131 | 0.010294 | 6 | 1.45 | 0.1429 | 0.4286 | ns |
| 12 | **NLGN1** | Synaptic Scaffolding | 0.015984 | 0.016336 | 9 | -0.16 | 0.6000 | 0.9000 | ns |
| 13 | **TCF4** | Neurodevelopmental/Transcription | 0.015231 | 0.016409 | 11 | -0.59 | 0.7500 | 0.9583 | ns |
| 14 | **FOXO3** | Neurodegeneration-Associated | 0.014869 | 0.014567 | 14 | 0.12 | 0.4667 | 0.8167 | ns |
| 15 | **BDNF** | Neurotrophic Signaling | 0.014019 | 0.014567 | 14 | -0.22 | 0.5333 | 0.8615 | ns |
| 16 | **MAPT** | Neurodegeneration-Associated | 0.013814 | 0.016336 | 9 | -1.15 | 0.9000 | 0.9583 | ns |
| 17 | **CREB1** | Neurotrophic Signaling | 0.012410 | 0.014567 | 14 | -0.88 | 0.8000 | 0.9583 | ns |
| 18 | **FMR1** | Neurodevelopmental/Transcription | 0.010911 | 0.013495 | 23 | -0.87 | 0.8750 | 0.9583 | ns |
| 19 | **SYN1** | Synaptic Vesicle/Transmission | 0.009434 | 0.013495 | 23 | -1.36 | 0.9167 | 0.9583 | ns |
| 20 | **COMT** | Catecholamine Metabolism | 0.009359 | 0.013495 | 23 | -1.39 | 0.9167 | 0.9583 | ns |
| 21 | **SHANK3** | Synaptic Scaffolding | 0.006302 | 0.013495 | 23 | -2.41 | 0.9583 | 0.9583 | ns |

*Significance: \*\*\* FDR<0.001, \*\* FDR<0.01, \* FDR<0.05, + FDR<0.10, ns = not significant*

---

## Expression Bin Details

| Bin | Freq Range | Intel Genes | Controls in Null | Intel Gene Names |
|-----|-----------|-------------|-----------------|------------------|
| 0 | 0.0141-0.2472 | 1 | 6 | SLC6A4 |
| 1 | 0.2472-0.4803 | 8 | 23 | GRIN2B, CAMK2A, SHANK3, HOMER1, FMR1, SYN1, SCN1A, COMT |
| 2 | 0.4803-0.7135 | 6 | 14 | MEF2C, BDNF, NRXN1, GRIN2A, FOXO3, CREB1 |
| 3 | 0.7135-0.9466 | 4 | 9 | APP, NEGR1, NLGN1, MAPT |
| 4 | 0.9466-1.1798 | 2 | 2 | CADM2, TCF4 |

---

## Null Distribution Statistics

- **N total controls:** 54
- **Mean:** 0.014011
- **Std:** 0.003294
- **5th percentile:** 0.008849
- **25th percentile:** 0.012889
- **Median:** 0.014172
- **75th percentile:** 0.015962
- **95th percentile:** 0.018500

### Per-Bin Null Statistics

| Bin | N Controls | Mean Shift | Std | Min | Max |
|-----|-----------|-----------|-----|-----|-----|
| 0 | 6 | 0.010294 | 0.004033 | 0.002326 | 0.014434 |
| 1 | 23 | 0.013495 | 0.002982 | 0.004207 | 0.018831 |
| 2 | 14 | 0.014567 | 0.002443 | 0.010396 | 0.018321 |
| 3 | 9 | 0.016336 | 0.002184 | 0.012866 | 0.020316 |
| 4 | 2 | 0.016739 | 0.000475 | 0.016264 | 0.017214 |

---

## Interpretation

While 1 intelligence genes show nominally significant effects (p < 0.05),
none survive FDR correction. The intelligence genes as a group show a trend toward
larger perturbation effects, but individual gene significance is limited by the
size of the null distribution.

### Comparison with Previous Paired Analysis

The original paired analysis (11 expression-matched controls) found:
- Mann-Whitney U: p = 0.031 (intelligence > controls)
- Cohen's d = 0.69 (medium effect)
- 9/10 intelligence genes > matched control (Wilcoxon p = 0.019)

This expanded permutation null (54 controls) provides a more
robust empirical null for gene-level inference, complementing the group-level
paired analysis.

---

## Technical Details

- Each control gene was perturbed in a **fresh Python subprocess** to avoid
  memory leaks and broken-pipe errors from Geneformer's internal multiprocessing.
- Results were saved incrementally after each gene completed.
- Failed genes were logged and skipped (not retried).
- Total controls attempted: 50
- Successful: 44
- Failed: 6

### Failed Control Genes

| Gene | Error |
|------|-------|
| NRDC | Timeout (600s) |
| TRIM2 | CUDA error: out of memory
Search for `cudaErrorMemoryAllocation' in https://docs |
| DNM3 | Timeout (600s) |
| ANK2 | Timeout (600s) |
| TCF12 | Timeout (600s) |
| PTPRD | Timeout (600s) |

---

*Generated: 2026-02-17 01:26:24*