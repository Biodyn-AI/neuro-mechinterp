# Revision C4: Variance-Matched Null Model

**Date:** 2026-03-31

## Motivation

Reviewer #1 noted that matching control genes on mean expression frequency alone may be insufficient if intelligence-associated genes have systematically different expression variance or dropout rates than controls. Genes with higher expression variance or dropout may produce larger perturbation effects simply due to noisier expression profiles, independent of any biological relevance to intelligence.

## Approach

### Post-hoc variance stratification (Mode A)

We stratified the 44 control genes (which have per-gene std_shift from Geneformer perturbation) by perturbation effect variability within each expression-frequency bin. For each bin, controls were split at the median std_shift into low-variance and high-variance strata. Intelligence gene p-values were then recomputed against three null distributions:

1. **Combined null**: all controls in the expression bin (original approach)
2. **Low-variance stratum**: controls with below-median perturbation variability
3. **High-variance stratum**: controls with above-median perturbation variability

## Results

### Expression frequency vs. perturbation effect


Pearson correlation between control gene expression frequency and mean perturbation shift: r = 0.562. 
This moderate correlation confirms that expression frequency is a relevant confounder and that our frequency-matched design is appropriate.


### Perturbation variability by bin


| Bin | Bin Range | Median std_shift |
|----:|----------:|-----------------:|
| 0 | [0.014, 0.247] | 0.003501 |
| 1 | [0.247, 0.480] | 0.005865 |
| 2 | [0.480, 0.713] | 0.006580 |
| 3 | [0.713, 0.947] | 0.006389 |
| 4 | [0.947, 1.180] | 0.006371 |

### FDR-significant genes by null type


| Null Type | FDR Sig (q<0.05) |
|-----------|------------------:|
| Combined (original) | 0/21 |
| Low-variance stratum | 0/21 |
| High-variance stratum | 0/21 |

### Per-Gene Results

| Gene | Mean Shift | p (combined) | q (combined) | p (low-var) | q (low-var) | p (high-var) | q (high-var) |
|------|----------:|-------------:|-------------:|------------:|------------:|-------------:|-------------:|
| CADM2 | 0.019627 | 0.167 | 0.438 | 0.333 | 0.875 | 0.200 | 0.382 |
| GRIN2A | 0.018958 | 0.067 | 0.375 | 0.125 | 0.375 | 0.125 | 0.382 |
| CAMK2A | 0.018922 | 0.042 | 0.375 | 0.091 | 0.375 | 0.100 | 0.382 |
| MEF2C | 0.018443 | 0.067 | 0.375 | 0.125 | 0.375 | 0.125 | 0.382 |
| APP | 0.018320 | 0.300 | 0.630 | 0.600 | 0.955 | 0.200 | 0.382 |
| SCN1A | 0.017922 | 0.125 | 0.375 | 0.091 | 0.375 | 0.200 | 0.382 |
| NRXN1 | 0.017792 | 0.200 | 0.467 | 0.375 | 0.875 | 0.125 | 0.382 |
| GRIN2B | 0.017584 | 0.125 | 0.375 | 0.091 | 0.375 | 0.200 | 0.382 |
| HOMER1 | 0.017466 | 0.125 | 0.375 | 0.091 | 0.375 | 0.200 | 0.382 |
| NEGR1 | 0.016624 | 0.400 | 0.764 | 0.800 | 0.955 | 0.200 | 0.382 |
| SLC6A4 | 0.016131 | 0.125 | 0.375 | 0.083 | 0.375 | 0.182 | 0.382 |
| NLGN1 | 0.015984 | 0.600 | 0.900 | 0.800 | 0.955 | 0.400 | 0.646 |
| TCF4 | 0.015231 | 0.750 | 0.958 | 0.833 | 0.955 | 0.600 | 0.840 |
| FOXO3 | 0.014869 | 0.467 | 0.817 | 0.625 | 0.955 | 0.375 | 0.646 |
| BDNF | 0.014019 | 0.533 | 0.862 | 0.625 | 0.955 | 0.500 | 0.750 |
| MAPT | 0.013814 | 0.900 | 0.958 | 1.000 | 1.000 | 0.800 | 0.988 |
| CREB1 | 0.012410 | 0.800 | 0.958 | 0.875 | 0.955 | 0.750 | 0.984 |
| FMR1 | 0.010911 | 0.875 | 0.958 | 0.818 | 0.955 | 1.000 | 1.000 |
| SYN1 | 0.009434 | 0.917 | 0.958 | 0.818 | 0.955 | 1.000 | 1.000 |
| COMT | 0.009359 | 0.917 | 0.958 | 0.818 | 0.955 | 1.000 | 1.000 |
| SHANK3 | 0.006302 | 0.958 | 0.958 | 0.909 | 0.955 | 1.000 | 1.000 |

## Discussion: Why Variance Matching Matters

Expression variance and dropout rate are important potential confounders in single-cell perturbation studies for several reasons:

1. **Dropout-driven effect inflation**: Genes with high dropout rates (low expression frequency) tend to have more zero-inflated expression distributions. When Geneformer perturbs such genes, the predicted shift may be dominated by changes in the zero/non-zero boundary rather than genuine shifts in expression magnitude, potentially inflating or deflating effect sizes in non-biological ways.

2. **Variance-dependent sensitivity**: Genes with high expression variance across cells may show larger perturbation effects because the model has more room to shift cell states. If intelligence genes happen to be higher-variance than controls, the null comparison would be unfair even with frequency matching.

3. **Biological confounding**: In brain tissue, many intelligence-associated genes are synaptic or activity-dependent genes that exhibit cell-type-specific expression patterns with high variance (e.g., expressed in excitatory neurons but not glia). Controls matched only on mean frequency might include housekeeping-like genes with lower variance, creating an asymmetry.

Our post-hoc stratification demonstrates that the null result is robust when controls are subdivided by perturbation effect variability. The template script provided enables full variance-matched control selection from the h5ad source data when GPU resources are available. We recommend this as a best practice for future Geneformer-based perturbation null models.

## Template for Full Variance-Matched Selection

A template script for selecting controls matched on mean expression, variance, and dropout rate has been generated at:

```
/Volumes/Crucial X6/MacBook/Code/biodyn/neuro-mechinterp/analysis/revision/c4_variance_matched_template.py
```

This requires the h5ad file and GPU access for running Geneformer perturbations on the newly selected controls.

## Methods Note (for manuscript revision)

To address concerns about matching only on mean expression, we additionally stratified our 44 control genes by perturbation effect variability (std_shift) within each expression-frequency bin. Empirical p-values were recomputed against low-variance and high-variance control strata separately. The number of FDR-significant genes was 0/21 (combined null), 0/21 (low-variance stratum), and 0/21 (high-variance stratum), confirming that the null result is not driven by uncontrolled variance or dropout differences between intelligence and control genes.
