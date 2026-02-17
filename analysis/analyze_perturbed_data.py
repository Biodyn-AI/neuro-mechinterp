#!/usr/bin/env python3
"""Analyze the perturbed data structure more thoroughly"""
import pickle
import numpy as np
from pathlib import Path

def analyze_pickle_data(gene_symbol):
    """Thoroughly analyze what's in the perturbed pickle files"""
    pickle_path = f"/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/perturb_{gene_symbol}_1000"
    pickle_dir = Path(pickle_path)
    
    if not pickle_dir.exists():
        print(f"❌ Directory not found: {pickle_dir}")
        return None
    
    pickle_files = list(pickle_dir.glob("*.pickle"))
    if not pickle_files:
        print(f"❌ No pickle files found for {gene_symbol}")
        return None
    
    pickle_file = pickle_files[0]
    
    print(f"\n{'='*60}")
    print(f"Analyzing {gene_symbol}")
    print(f"{'='*60}")
    print(f"File: {pickle_file}")
    
    try:
        with open(pickle_file, 'rb') as f:
            data = pickle.load(f)
        
        print(f"\nData type: {type(data)}")
        print(f"Data keys: {list(data.keys())}")
        
        for key, value in data.items():
            print(f"\nKey: {key}")
            print(f"  Type: {type(value)}")
            
            if isinstance(value, list):
                print(f"  List length: {len(value)}")
                if len(value) > 0:
                    print(f"  First element type: {type(value[0])}")
                    
                    # Check if the list contains embeddings
                    if hasattr(value[0], '__len__') and not isinstance(value[0], str):
                        if hasattr(value[0], 'shape'):
                            print(f"  First element shape: {value[0].shape}")
                        else:
                            print(f"  First element length: {len(value[0])}")
                    
                    # Sample some values
                    if isinstance(value[0], (int, float, np.number)):
                        sample_size = min(10, len(value))
                        print(f"  Sample values: {value[:sample_size]}")
                    elif hasattr(value[0], '__getitem__'):
                        try:
                            elem = value[0]
                            if hasattr(elem, 'shape'):
                                print(f"  Element 0 shape: {elem.shape}")
                                print(f"  Element 0 sample: {elem[:5]}")
                            else:
                                sample = elem[:5] if len(elem) > 5 else elem
                                print(f"  Element 0 sample: {sample}")
                        except Exception as e:
                            print(f"  Could not sample: {e}")
            
            elif hasattr(value, 'shape'):
                print(f"  Shape: {value.shape}")
                print(f"  Sample: {value[:5] if len(value) > 5 else value}")
            
            elif hasattr(value, '__len__') and not isinstance(value, str):
                print(f"  Length: {len(value)}")
                if hasattr(value, '__getitem__'):
                    try:
                        sample = value[:3] if len(value) > 3 else value
                        print(f"  Sample: {sample}")
                    except:
                        pass
        
        return data
    
    except Exception as e:
        print(f"❌ Error loading {gene_symbol}: {e}")
        import traceback
        traceback.print_exc()
        return None

def investigate_embedding_dimensions():
    """Check if all embeddings have the same dimensions"""
    genes = ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]
    
    print("\n" + "="*80)
    print("EMBEDDING DIMENSIONS ANALYSIS")
    print("="*80)
    
    dimensions = {}
    
    for gene in genes:
        data = analyze_pickle_data(gene)
        
        if data:
            for key, value in data.items():
                if isinstance(value, list) and len(value) > 0:
                    first_elem = value[0]
                    if hasattr(first_elem, '__len__') and not isinstance(first_elem, str):
                        if hasattr(first_elem, 'shape'):
                            dim = first_elem.shape[0] if len(first_elem.shape) == 1 else first_elem.shape
                        else:
                            dim = len(first_elem)
                        dimensions[gene] = dim
                        print(f"{gene}: dimension = {dim}")
                        break
    
    print(f"\nDimension summary: {dimensions}")
    
    # Check if dimensions are consistent
    if len(set(str(d) for d in dimensions.values())) == 1:
        print("✅ All genes have consistent embedding dimensions")
    else:
        print("⚠️ Inconsistent embedding dimensions across genes")
    
    return dimensions

def test_embedding_structure():
    """Test what the embedding structure actually contains"""
    print("\n" + "="*80)
    print("EMBEDDING STRUCTURE TEST")
    print("="*80)
    
    # Just test one gene for now
    test_gene = "CADM2"
    data = analyze_pickle_data(test_gene)
    
    if not data:
        return False
    
    # Get the embedding data
    key = list(data.keys())[0]
    embeddings = data[key]
    
    print(f"\nTesting with {test_gene} embeddings:")
    print(f"Number of embeddings: {len(embeddings)}")
    
    if len(embeddings) >= 2:
        emb1 = np.array(embeddings[0])
        emb2 = np.array(embeddings[1])
        
        print(f"Embedding 1 shape: {emb1.shape}")
        print(f"Embedding 2 shape: {emb2.shape}")
        
        # Compute cosine similarity between two embeddings
        from sklearn.metrics.pairwise import cosine_similarity
        
        sim = cosine_similarity(emb1.reshape(1, -1), emb2.reshape(1, -1))[0, 0]
        print(f"Cosine similarity between embeddings 1 and 2: {sim:.4f}")
        
        # Check if these are reasonable embedding values
        print(f"Embedding 1 range: [{emb1.min():.4f}, {emb1.max():.4f}]")
        print(f"Embedding 1 mean: {emb1.mean():.4f}")
        print(f"Embedding 1 std: {emb1.std():.4f}")
        
        # These are likely perturbed embeddings, so let's see if we can use them differently
        print(f"\n🤔 Interpretation: These appear to be perturbed embeddings only.")
        print(f"   Without control embeddings, we can't compute shifts directly.")
        print(f"   However, we might be able to analyze relative differences between genes.")
    
    return True

if __name__ == "__main__":
    print("🔍 COMPREHENSIVE PERTURBED DATA ANALYSIS")
    print("="*80)
    
    # Analyze dimensions across all genes
    dimensions = investigate_embedding_dimensions()
    
    # Test the actual structure
    test_embedding_structure()
    
    print(f"\n{'='*80}")
    print("SUMMARY")
    print("="*80)
    print("✅ Successfully analyzed perturbed embedding files")
    print("❌ No control embeddings found in the files")
    print("💡 Next step: Need to generate control embeddings or use alternative analysis")