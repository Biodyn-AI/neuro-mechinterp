# Datasets & Tools for Virtual Perturbation Screening Project

**Date:** 2026-02-13

---

## 1. Foundation Models

### Geneformer V2 (RECOMMENDED — Primary)
- **HuggingFace:** https://huggingface.co/ctheodoris/Geneformer
- **Models available:**
  - Geneformer-V1-10M (original, ~30M cells, 10M params)
  - **Geneformer-V2-104M** (~104M cells, 104M params)
  - **Geneformer-V2-316M** (~104M cells, 316M params) ← default/recommended
  - Geneformer-V2-104M_CLcancer (cancer continual learning variant)
- **Input:** Rank value encoding of single-cell transcriptomes (max 4096 genes for V2)
- **In-silico perturbation:** Built-in zero-shot capability; delete/overexpress genes and measure embedding shift
- **Documentation:** https://geneformer.readthedocs.io
- **Installation:** `git clone https://huggingface.co/ctheodoris/Geneformer && cd Geneformer && pip install .`
- **Paper:** Theodoris et al. 2023, Nature. PMID: 37258680
- **V2 preprint:** https://www.biorxiv.org/content/10.1101/2024.08.16.608180v1

### scGPT (Alternative backbone)
- **GitHub:** https://github.com/bowang-lab/scGPT
- **HuggingFace:** https://huggingface.co/bowang-lab (check for latest)
- **Paper:** Cui et al. 2024, Nat Methods. PMID: 38409223
- **Perturbation:** Has perturbation prediction module

### Mouse-Geneformer (for validation in mouse models)
- **Paper:** Ito et al. 2025, PLoS Genet. PMID: 40106407
- Mouse-specific version for cross-species validation

---

## 2. Brain Single-Cell RNA-seq Datasets

### Siletti et al. 2023 — Adult Human Brain Atlas (PRIMARY)
- **Paper:** Science 2023, 382:eadf7044
- **Data portal:** https://portal.brain-map.org/
- **CellxGene:** https://cellxgene.cziscience.com/ (search "Siletti" or "adult human brain")
- **Content:** 3M+ nuclei from entire adult human brain, 461 superclusters, >3000 clusters
- **Key cell types for us:** UL-IT (upper-layer intratelencephalic) = L2/3 pyramidal neurons of DLPFC
- **Format:** h5ad (AnnData)

### Allen Brain Cell Atlas (BICCN)
- **Portal:** https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas
- **SEA-AD (Aging/Dementia):** https://portal.brain-map.org/explore/seattle-alzheimers-disease
  - DLPFC snRNA-seq from aged donors with varying cognitive status
  - Could be useful for comparing intelligence-associated signatures in aging
- **Human MTG (Middle Temporal Gyrus):** https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-smart-seq
- **Human M1 (Primary Motor Cortex):** https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x

### PsychENCODE
- **Portal:** http://resource.psychencode.org/
- **Content:** Multi-region brain snRNA-seq, DLPFC prominently featured
- **Relevance:** Psychiatric phenotype annotations; some donors have cognitive data

### Human Cell Atlas — Brain
- **Portal:** https://www.humancellatlas.org/
- **CellxGene:** https://cellxgene.cziscience.com/collections (filter by brain)

---

## 3. Intelligence / Cognition GWAS Summary Statistics

### Savage et al. 2018 (PRIMARY)
- **Paper:** Nat Genet 50:912-919. PMID: 29942086
- **Summary stats:** https://ctg.cncr.nl/software/summary_statistics
- **Direct:** Search "Intelligence" at https://ctg.cncr.nl/software/summary_statistics
- **N:** 269,867 individuals
- **Loci:** 205 genomic loci, 1,016 genes

### Davies et al. 2018
- **Paper:** Nat Commun 9:767
- **N:** 300,486
- **Summary stats:** https://www.ed.ac.uk/psychology/research/research-themes/individual-differences

### Hill et al. 2019
- **Paper:** Mol Psychiatry 24:169-181
- **N:** 248,482
- **Summary stats:** https://ctg.cncr.nl/software/summary_statistics

### UK Biobank Cognitive Phenotypes
- **Portal:** https://www.ukbiobank.ac.uk/
- **GWAS catalog:** https://www.ebi.ac.uk/gwas/
- Search for: "cognitive performance", "fluid intelligence", "reaction time"
- **Neale lab UKB GWAS:** http://www.nealelab.is/uk-biobank/ (Field 20016 = fluid intelligence)

### MAGMA Gene-Level Results
- Use MAGMA (https://ctg.cncr.nl/software/magma) to convert SNP-level to gene-level statistics
- Alternatively, FUMA (https://fuma.ctglab.nl/) for web-based gene mapping

---

## 4. Drug Databases

### DrugBank
- **URL:** https://go.drugbank.com/
- **Download:** https://go.drugbank.com/releases/latest (requires academic license)
- **Content:** ~15,000 drug entries with targets, mechanisms, clinical status
- **API:** Available for programmatic access

### ChEMBL
- **URL:** https://www.ebi.ac.uk/chembl/
- **Download:** https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/
- **Content:** Bioactivity data, drug-target relationships
- **Python:** `pip install chembl_webresource_client`

### DGIdb (Drug-Gene Interaction Database)
- **URL:** https://www.dgidb.org/
- **API:** https://dgidb.org/api
- **Content:** Curated drug-gene interactions from multiple sources
- **Best for:** Quick lookup of druggable genes

### OpenTargets
- **URL:** https://www.opentargets.org/
- **Content:** Target-disease associations, drug tractability assessments
- **GraphQL API:** https://api.platform.opentargets.org/

### STITCH (Chemical-Protein Interactions)
- **URL:** http://stitch.embl.de/
- **Content:** Chemical-protein interaction network

### Connectivity Map (CMap) / LINCS L1000
- **URL:** https://clue.io/
- **Content:** Gene expression signatures of drug perturbations
- **Relevance:** Could validate whether predicted targets are affected by known nootropics

---

## 5. Gene Network & Pathway Resources

### STRING
- **URL:** https://string-db.org/
- Protein-protein interaction networks

### Enrichr / Gene Ontology
- **URL:** https://maayanlab.cloud/Enrichr/
- Pathway enrichment analysis of identified targets

### GWAS Catalog
- **URL:** https://www.ebi.ac.uk/gwas/
- Cross-reference targets with other trait associations

---

## 6. Software / Analysis Tools

### Geneformer In-Silico Perturbation
```python
# Key classes from Geneformer package:
from geneformer import InSilicoPerturber, InSilicoPerturberStats
# See examples at: https://huggingface.co/ctheodoris/Geneformer/tree/main/examples
```

### Scanpy (Single-cell analysis)
```bash
pip install scanpy
```

### MAGMA (Gene-level GWAS analysis)
- https://ctg.cncr.nl/software/magma
- Convert SNP-level GWAS to gene-level scores

### LDSC (LD Score Regression)
- https://github.com/bulik/ldsc
- Cell-type-specific heritability enrichment

### FUMA
- https://fuma.ctglab.nl/
- Web-based GWAS functional annotation

---

## 7. Compute Requirements (Estimated)

| Step | GPU Memory | Time (est.) |
|------|-----------|-------------|
| Tokenize brain atlas (~100K DLPFC cells) | CPU only | ~1 hour |
| Geneformer V2-316M perturbation (20K genes × 100K cells) | 40-80 GB GPU | ~24-72 hours |
| Geneformer V2-104M perturbation (lighter) | 16-24 GB GPU | ~8-24 hours |
| Downstream analysis (scoring, ranking, drug mapping) | CPU | ~2-4 hours |

**Recommended:** A100 80GB or H100 GPU. Cloud options: Google Colab Pro+ (A100), AWS p4d, GCP a2-highgpu.
