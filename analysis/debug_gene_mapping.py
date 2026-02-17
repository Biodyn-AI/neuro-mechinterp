#!/usr/bin/env python3
"""
Debug script to examine gene name formats and fix mapping issues
"""
import sys
import os
import pickle
from pathlib import Path

# Project paths
PROJ_ROOT = Path(r"D:\openclaw\intelligence-augmentation")
DATA_DIR = PROJ_ROOT / "data"
MODEL_DIR = PROJ_ROOT / "models" / "Geneformer"
BRAIN_DATA = DATA_DIR / "brain_scrna" / "DLPFC_11k.h5ad"
GENE_DICT_DIR = MODEL_DIR / "geneformer" / "gene_dictionaries_30m"

def examine_gene_formats():
    """Compare gene names between data and dictionaries"""
    print("=== Gene Format Investigation ===")
    
    # Load brain data
    print("\n1. Loading brain data...")
    import scanpy as sc
    adata = sc.read_h5ad(BRAIN_DATA)
    
    brain_genes = list(adata.var_names)
    print(f"Brain data genes: {len(brain_genes)}")
    print(f"Sample brain genes: {brain_genes[:10]}")
    
    # Check gene format
    is_ensembl = brain_genes[0].startswith("ENSG") if brain_genes else False
    print(f"Gene format appears to be: {'Ensembl' if is_ensembl else 'Gene symbols'}")
    
    # Load Geneformer dictionaries
    print("\n2. Loading Geneformer gene dictionaries...")
    with open(GENE_DICT_DIR / "gene_name_id_dict_gc30M.pkl", "rb") as f:
        gene_name_to_id = pickle.load(f)
    
    with open(GENE_DICT_DIR / "token_dictionary_gc30M.pkl", "rb") as f:
        token_dict = pickle.load(f)
    
    print(f"Gene name->ID mappings: {len(gene_name_to_id)}")
    print(f"Token dictionary: {len(token_dict)}")
    
    # Sample gene names from dictionary
    sample_gene_names = list(gene_name_to_id.keys())[:10]
    print(f"Sample dictionary gene names: {sample_gene_names}")
    
    sample_ensembl_ids = list(token_dict.keys())[:10]
    print(f"Sample Ensembl IDs in token dict: {sample_ensembl_ids}")
    
    # Check overlap
    print("\n3. Checking gene overlap...")
    if is_ensembl:
        # Brain genes are Ensembl IDs
        direct_overlap = set(brain_genes) & set(token_dict.keys())
        print(f"Direct Ensembl ID overlap: {len(direct_overlap)} genes")
        print(f"Sample overlapping genes: {list(direct_overlap)[:5]}")
    else:
        # Brain genes are symbols - need to map via gene_name_to_id
        mapped_count = 0
        mapped_examples = []
        
        for gene in brain_genes[:1000]:  # Check first 1000
            # Try different variants
            variants = [gene, gene.upper(), gene.lower()]
            for variant in variants:
                if variant in gene_name_to_id:
                    ensembl_id = gene_name_to_id[variant]
                    if ensembl_id in token_dict:
                        mapped_count += 1
                        if len(mapped_examples) < 5:
                            mapped_examples.append(f"{gene} -> {ensembl_id}")
                        break
        
        print(f"Symbol->Ensembl mapping success (first 1000): {mapped_count}")
        print(f"Sample mappings: {mapped_examples}")
    
    # Check for specific intelligence genes
    print("\n4. Checking intelligence genes...")
    INTEL_GENES = ["FOXO3", "MEF2C", "COMT", "BDNF", "CADM2"]
    
    for gene in INTEL_GENES:
        if is_ensembl:
            # Look for gene in name_to_id mapping
            if gene in gene_name_to_id:
                ensembl_id = gene_name_to_id[gene]
                in_tokens = ensembl_id in token_dict
                in_brain = ensembl_id in brain_genes
                print(f"  {gene} -> {ensembl_id}: in_tokens={in_tokens}, in_brain={in_brain}")
            else:
                print(f"  {gene}: not in gene_name_to_id")
        else:
            in_brain = gene in brain_genes
            if gene in gene_name_to_id:
                ensembl_id = gene_name_to_id[gene]
                in_tokens = ensembl_id in token_dict
                print(f"  {gene}: in_brain={in_brain}, maps_to={ensembl_id}, in_tokens={in_tokens}")
            else:
                print(f"  {gene}: in_brain={in_brain}, not_in_mapping=True")

if __name__ == "__main__":
    examine_gene_formats()