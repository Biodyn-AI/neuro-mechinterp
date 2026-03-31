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

---

## Quick Start

### 1. Environment Setup

```bash
# Create a fresh conda environment (Python 3.10 recommended for torch/geneformer compatibility)
conda create -n neuro-mechinterp python=3.10 -y
conda activate neuro-mechinterp

# Install PyTorch with CUDA (adjust cuda version to match your driver)
pip install torch --index-url https://download.pytorch.org/whl/cu118

# Install all other dependencies
pip install -r requirements.txt
```

**Note on Geneformer:** Install from the official Hugging Face repo if the pip package is unavailable:
```bash
pip install git+https://huggingface.co/ctheodoris/Geneformer
```

**Note on scGPT:** Install from the official repo:
```bash
pip install git+https://github.com/bowang-lab/scGPT.git
```

### 2. Data

The pipeline requires the following input data placed under `data/`:

| File | Description | Source |
|------|-------------|--------|
| `data/brain_scrna/DLPFC_11k.h5ad` | Human DLPFC single-cell RNA-seq (~11k cells) | [SEA-AD / Allen Brain Map](https://portal.brain-map.org/atlases-and-data/rnaseq) |
| `data/gwas/SavageJansen_2018_intelligence_metaanalysis.txt.gz` | Intelligence GWAS summary statistics | [Savage et al. 2018](https://doi.org/10.1038/s41588-018-0152-6) |
| `data/drugs/interactions.tsv` | Gene-drug interactions | [DGIdb](https://www.dgidb.org/downloads) |

**Model checkpoints** (downloaded automatically or manually):
- Geneformer V1-10M from [Hugging Face](https://huggingface.co/ctheodoris/Geneformer)
- scGPT brain checkpoint from [scGPT releases](https://github.com/bowang-lab/scGPT)

### 3. Running the Analysis Pipeline

All scripts are in `analysis/`. Run them in this order:

```bash
cd analysis/

# Step 1 — Pilot run (200 cells, 50 genes, fast sanity check)
python pilot_perturbation.py

# Step 2 — Full Geneformer perturbation pipeline (500 cells, 21 intelligence genes)
python full_perturbation_pipeline.py

# Step 3 — Expression-matched random-gene controls
python expression_matched_controls.py

# Step 4 — Permutation null distribution and empirical p-values
python permutation_null.py

# Step 5 — Cross-model validation with scGPT
python scgpt_perturbation_crossmodel.py

# Step 6 — Cell-type-specific analysis
python celltype_pathway_analysis.py

# Step 7 — Combinatorial perturbation analysis
python combo_analysis.py

# Step 8 — Figures and summary
python regenerate_figures.py
```

### 4. Expected Outputs

- `analysis/results/` — perturbation scores (CSV/JSON), gene rankings, drug candidates
- `analysis/figures/` — publication-ready figures (gene ranking bar plots, pathway boxplots, heatmaps, violin plots)
- `analysis/FINAL_REPORT.md` — full results narrative with tables and statistics

Key result files:
| File | Contents |
|------|----------|
| `results/perturbation_results_streamlined.csv` | Per-gene embedding shift scores |
| `results/drug_candidates_streamlined.csv` | Drug repurposing candidates |
| `figures/fig1_gene_ranking.png` | Ranked perturbation effects |
| `figures/fig2_pathway_boxplot.png` | Pathway-level aggregation |

### 5. Hardware Requirements

- **GPU:** NVIDIA GPU with >= 6 GB VRAM (developed on RTX 2060 6 GB). CUDA 11.8+ required.
- **RAM:** >= 16 GB system RAM.
- **Disk:** ~5 GB for model checkpoints + data.
- **Runtime:** Full pipeline takes approximately 4-8 hours on a single GPU. The pilot script (`pilot_perturbation.py`) finishes in ~15 minutes and is useful for verifying the setup.
