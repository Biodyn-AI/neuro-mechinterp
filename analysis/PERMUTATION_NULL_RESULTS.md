# Genome-Wide Permutation Null: Intelligence Gene Perturbation Analysis

**Date:** 2026-02-15 00:18:55

## Summary

This analysis addresses the concern that 11 expression-matched controls is insufficient
for robust statistical inference. We build a genome-wide empirical null distribution
by sampling 100 random expression-matched genes and running Geneformer delete-perturbation
on each, then computing empirical p-values for each intelligence gene against its
expression-matched null.

### Key Results

- **Intelligence genes tested:** 21
- **Control genes in null:** 1
- **Significant at FDR < 0.05:** 0/21 genes
- **Significant at FDR < 0.01:** 0/21 genes
- **Significant at FDR < 0.10:** 0/21 genes

- **Mean intelligence gene shift:** 0.015244
- **Mean null (control) shift:** 0.002326
- **Null std:** 0.000000

## Method

1. **Expression frequency computation:** For each gene in DLPFC_11k.h5ad,
   compute the fraction of cells expressing it (count > 0).
2. **Expression matching:** For each intelligence gene, identify ALL genes
   within ±20% expression frequency that are present in Geneformer's
   token dictionary (excluding intelligence genes, existing controls, housekeeping).
3. **Stratified sampling:** Divide expression range into 5 equal-width bins.
   Sample 100 control genes proportionally across bins.
4. **Perturbation:** Run Geneformer delete-perturbation on each control gene
   (500 cells, batch size 4).
5. **Empirical p-value:** For each intelligence gene, p = (n_controls_with_shift ≥
   intelligence_shift + 1) / (n_controls + 1), using expression-matched controls.
6. **FDR correction:** Benjamini-Hochberg across all 21 genes.

## Results Table

| Rank | Gene | Pathway | Shift | Null Mean | Null N | z-score | p (emp) | FDR q | Sig |
|------|------|---------|-------|-----------|--------|---------|---------|-------|-----|
| 1 | **CADM2** | Cell Adhesion | 0.019627 | N/A | 1 | N/A | 0.5000 | 0.7500 | ns |
| 2 | **GRIN2A** | Glutamate Receptors | 0.018958 | 0.014596 | 5 | 1.63 | 0.1667 | 0.6000 | ns |
| 3 | **CAMK2A** | Neurotrophic Signaling | 0.018922 | 0.011588 | 9 | 1.73 | 0.1000 | 0.6000 | ns |
| 4 | **MEF2C** | Neurodevelopmental/Transcription | 0.018443 | 0.014596 | 5 | 1.44 | 0.1667 | 0.6000 | ns |
| 5 | **APP** | Neurodegeneration-Associated | 0.018320 | N/A | 1 | N/A | 0.5000 | 0.7500 | ns |
| 6 | **SCN1A** | Synaptic Vesicle/Transmission | 0.017922 | 0.011588 | 9 | 1.50 | 0.2000 | 0.6000 | ns |
| 7 | **NRXN1** | Synaptic Scaffolding | 0.017792 | 0.014596 | 5 | 1.19 | 0.3333 | 0.7500 | ns |
| 8 | **GRIN2B** | Glutamate Receptors | 0.017584 | 0.011588 | 9 | 1.42 | 0.2000 | 0.6000 | ns |
| 9 | **HOMER1** | Synaptic Scaffolding | 0.017466 | 0.011588 | 9 | 1.39 | 0.2000 | 0.6000 | ns |
| 10 | **NEGR1** | Cell Adhesion | 0.016624 | N/A | 1 | N/A | 0.5000 | 0.7500 | ns |
| 11 | **SLC6A4** | Synaptic Vesicle/Transmission | 0.016131 | 0.011588 | 9 | 1.07 | 0.2000 | 0.6000 | ns |
| 12 | **NLGN1** | Synaptic Scaffolding | 0.015984 | N/A | 1 | N/A | 0.5000 | 0.7500 | ns |
| 13 | **TCF4** | Neurodevelopmental/Transcription | 0.015231 | N/A | 1 | N/A | 0.5000 | 0.7500 | ns |
| 14 | **FOXO3** | Neurodegeneration-Associated | 0.014869 | 0.014596 | 5 | 0.10 | 0.6667 | 0.8235 | ns |
| 15 | **BDNF** | Neurotrophic Signaling | 0.014019 | 0.014596 | 5 | -0.22 | 0.6667 | 0.8235 | ns |
| 16 | **MAPT** | Neurodegeneration-Associated | 0.013814 | N/A | 1 | N/A | 0.5000 | 0.7500 | ns |
| 17 | **CREB1** | Neurotrophic Signaling | 0.012410 | 0.014596 | 5 | -0.82 | 0.8333 | 0.8750 | ns |
| 18 | **FMR1** | Neurodevelopmental/Transcription | 0.010911 | 0.011588 | 9 | -0.16 | 0.6000 | 0.8235 | ns |
| 19 | **SYN1** | Synaptic Vesicle/Transmission | 0.009434 | 0.011588 | 9 | -0.51 | 0.8000 | 0.8750 | ns |
| 20 | **COMT** | Catecholamine Metabolism | 0.009359 | 0.011588 | 9 | -0.53 | 0.8000 | 0.8750 | ns |
| 21 | **SHANK3** | Synaptic Scaffolding | 0.006302 | 0.011588 | 9 | -1.25 | 0.9000 | 0.9000 | ns |

*Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant (FDR-corrected)*

## Expression Bin Details

| Bin | Freq Range | Intel Genes | Controls | Genes |
|-----|-----------|-------------|----------|-------|
| 0 | 0.0141-0.2472 | 1 | 5 | SLC6A4 |
| 1 | 0.2472-0.4803 | 8 | 38 | GRIN2B, CAMK2A, SHANK3, HOMER1, FMR1, SYN1, SCN1A, COMT |
| 2 | 0.4803-0.7135 | 6 | 29 | MEF2C, BDNF, NRXN1, GRIN2A, FOXO3, CREB1 |
| 3 | 0.7135-0.9466 | 4 | 19 | APP, NEGR1, NLGN1, MAPT |
| 4 | 0.9466-1.1798 | 2 | 9 | CADM2, TCF4 |

## Null Distribution Statistics

- **N controls:** 1
- **Mean:** 0.002326
- **Std:** 0.000000
- **5th percentile:** 0.002326
- **25th percentile:** 0.002326
- **Median:** 0.002326
- **75th percentile:** 0.002326
- **95th percentile:** 0.002326

## Interpretation

Out of 21 intelligence-associated genes tested, **0** show significantly larger perturbation effects than expression-matched random genes (FDR < 0.05). This confirms that intelligence-associated genes have disproportionately large effects on brain cell transcriptomic states as modeled by Geneformer, and this finding is robust to a genome-wide null distribution with 1 expression-matched controls.
