# Virtual Perturbation Screening for Cognitive Enhancement Targets

## Project Goal
Use single-cell foundation models (Geneformer/scGPT) to perform in-silico perturbation screening 
on brain tissue data, identifying novel gene targets that shift cell states toward high-cognition 
expression profiles. Cross-reference with drug databases for repurposing candidates.

## Approach 
1. Obtain brain single-cell RNA-seq data (DLPFC pyramidal neurons, astrocytes, oligodendrocytes)
2. Obtain intelligence GWAS summary statistics (Savage et al. 2018, ~1000 loci)
3. Use Geneformer or scGPT for virtual gene perturbation (knockout/overexpression of all ~20K genes)
4. Score perturbation effects against intelligence-associated gene signatures
5. Identify master regulators and hub targets
6. Cross-reference with DrugBank/ChEMBL for repurposing
7. Validate findings against known intelligence biology (BDNF, CREB1, Wnt, etc.)

## Author
Mr. Ihor Kendiukhov, University of Tuebingen, Computer Science Department
Email: kendiukhov@gmail.com

## Key References
- Savage et al. 2018 — Intelligence GWAS meta-analysis
- Theodoris et al. 2023 — Geneformer
- Cui et al. 2024 — scGPT

## Folders
- literature/ — related work analysis, novelty assessment
- data/ — downloaded datasets (GWAS, scRNA-seq, drug databases)
- analysis/ — computational pipeline and results
- paper/ — manuscript drafts
