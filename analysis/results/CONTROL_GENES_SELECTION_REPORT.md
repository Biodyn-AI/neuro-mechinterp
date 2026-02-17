# Expression-Matched Controls Selection Summary

## Selected Control Genes

### Expression-Matched Controls (21 genes)
- **EXOC6B** (freq: 0.9628) matched to **CADM2** (freq: 0.9340, diff: 0.0288)
- **DOCK9** (freq: 0.8514) matched to **APP** (freq: 0.8340, diff: 0.0174)
- **PDE7B** (freq: 0.3437) matched to **FOXO3** (freq: 0.3220, diff: 0.0217)
- **ZDHHC9** (freq: 0.4649) matched to **NLGN1** (freq: 0.5000, diff: 0.0351)
- **HDAC11** (freq: 0.3013) matched to **MAPT** (freq: 0.3500, diff: 0.0487)
- **XRCC6** (freq: 0.3565) matched to **TCF4** (freq: 0.3720, diff: 0.0155)
- **WDR46** (freq: 0.0912) matched to **CREB1** (freq: 0.1280, diff: 0.0368)
- **CUEDC2** (freq: 0.1586) matched to **DISC1** (freq: 0.1100, diff: 0.0486)
- **IGSF10** (freq: 0.1143) matched to **NEGR1** (freq: 0.1160, diff: 0.0017)
- **INSIG1-DT** (freq: 0.0346) matched to **GRIN2A** (freq: 0.0720, diff: 0.0374)
- **SALL3** (freq: 0.0429) matched to **FMR1** (freq: 0.0540, diff: 0.0111)
- **MFSD6L** (freq: 0.0003) matched to **COMT** (freq: 0.0460, diff: 0.0457)
- **ENSG00000237070** (freq: 0.0040) matched to **BDNF** (freq: 0.0140, diff: 0.0100)
- **ENSG00000249856** (freq: 0.0422) matched to **CAMK2A** (freq: 0.0120, diff: 0.0302)
- **ENSG00000261845** (freq: 0.0423) matched to **HOMER1** (freq: 0.0140, diff: 0.0283)
- **ENSG00000261036** (freq: 0.0019) matched to **SHANK3** (freq: 0.0080, diff: 0.0061)
- **ENSG00000232811** (freq: 0.0163) matched to **NRXN1** (freq: 0.0080, diff: 0.0083)
- **TMEM88** (freq: 0.0072) matched to **MEF2C** (freq: 0.0060, diff: 0.0012)
- **ENSG00000230569** (freq: 0.0006) matched to **SLC6A4** (freq: 0.0040, diff: 0.0034)
- **ENSG00000261555** (freq: 0.0002) matched to **SCN1A** (freq: 0.0040, diff: 0.0038)
- **IGHGP** (freq: 0.0000) matched to **GRIN2B** (freq: 0.0020, diff: 0.0020)

### Housekeeping Controls (5 genes)
- **ACTB** (freq: 0.5090)
- **GAPDH** (freq: 0.4537)
- **B2M** (freq: 0.2327)
- **RPL13A** (freq: 0.2646)
- **YWHAZ** (freq: 0.5775)


## Expression Frequency Comparison

| Gene Set | Count | Mean Freq | Min Freq | Max Freq |
|----------|-------|-----------|----------|----------|
| Intelligence Genes | 22 | 0.1778 | 0.0020 | 0.9340 |
| Expression-Matched | 21 | 0.1827 | 0.0000 | 0.9628 |
| Housekeeping | 5 | 0.4075 | 0.2327 | 0.5775 |

## Quality Control

- Mean matching difference (expression-matched): 0.0210
- Max matching difference (expression-matched): 0.0487
- All control genes are unique and distinct from intelligence genes

## Next Steps

To complete the analysis, run the Geneformer perturbation pipeline on these selected control genes using the same methodology as the intelligence genes.

The control genes have been selected to match expression frequencies, ensuring that any differences in perturbation effects are not due to expression level artifacts.
