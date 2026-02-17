#!/usr/bin/env python3
"""Investigate the exact structure of Geneformer output"""
import pickle
import numpy as np
from pathlib import Path

def investigate_file(file_path):
    """Investigate structure of a pickle file"""
    print(f"🔍 Investigating: {file_path.name}")
    
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    
    print(f"Top level type: {type(data)}")
    print(f"Top level keys: {list(data.keys())}")
    
    # Get the data for the gene
    key = list(data.keys())[0]
    embeddings = data[key]
    
    print(f"\nEmbedding data:")
    print(f"  Type: {type(embeddings)}")
    print(f"  Length: {len(embeddings)}")
    
    if len(embeddings) > 0:
        first_item = embeddings[0]
        print(f"  First item type: {type(first_item)}")
        
        if hasattr(first_item, 'shape'):
            print(f"  First item shape: {first_item.shape}")
        elif isinstance(first_item, (list, tuple)):
            print(f"  First item length: {len(first_item)}")
            if len(first_item) > 0:
                print(f"  First sub-item type: {type(first_item[0])}")
        
        # Look at a few more items
        if len(embeddings) > 1:
            print(f"  Item types: {[type(item) for item in embeddings[:5]]}")
        
        # If items are arrays, check shapes
        if hasattr(first_item, 'shape'):
            shapes = [item.shape for item in embeddings[:5] if hasattr(item, 'shape')]
            print(f"  First 5 shapes: {shapes}")
            
    # The issue might be that this is the output format for cell embeddings
    # but not the perturbation comparison format
    print(f"\nThis appears to be cell embeddings for the perturbed condition.")
    print(f"We may need to compare against a baseline/control run or")
    print(f"the perturbation analysis needs to be run differently.")

def main():
    """Investigate all pickle files"""
    results_dir = Path("/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled")
    
    # Find all pickle files
    pickle_files = []
    for gene in ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]:
        gene_dir = results_dir / f"perturb_{gene}_1000"
        if gene_dir.exists():
            files = list(gene_dir.glob("*.pickle"))
            if files:
                pickle_files.append(files[0])
    
    print("Found pickle files:")
    for f in pickle_files:
        print(f"  {f}")
    
    # Investigate the first one in detail
    if pickle_files:
        investigate_file(pickle_files[0])
        
        print(f"\n" + "="*60)
        print("ANALYSIS:")
        print("="*60)
        print("The Geneformer InSilicoPerturber seems to output only the")
        print("perturbed cell embeddings, not a comparison structure.")
        print("")
        print("To compute perturbation effects, we need:")
        print("1. Control/baseline embeddings (unperturbed)")
        print("2. Perturbed embeddings (what we have)")
        print("")
        print("The comparison might need to be done against:")
        print("- A separate control run")
        print("- The original dataset embeddings")
        print("- Or use a different perturbation analysis approach")

if __name__ == "__main__":
    main()