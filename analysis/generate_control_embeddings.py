#!/usr/bin/env python3
"""Generate control (unperturbed) embeddings for 1000-cell dataset"""
import os
import sys
import time
import torch
import pickle
import numpy as np
from pathlib import Path

def generate_control_embeddings(n_cells=1000):
    """Generate control embeddings using Geneformer without perturbations"""
    print(f"\n{'='*60}")
    print(f"Generating control embeddings for {n_cells:,} cells")
    print(f"{'='*60}")
    
    # Memory settings for RTX 2060 (6GB)
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:256"
    
    # Clear GPU memory
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"Available memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    
    # Paths
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    token_dict = "/mnt/d/openclaw/intelligence-augmentation/models/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
    dataset_path = f"/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/brain_{n_cells}.dataset"
    output_dir = f"/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/control_{n_cells}"
    
    print(f"Dataset: {dataset_path}")
    print(f"Output: {output_dir}")
    
    # Check if dataset exists
    if not os.path.exists(dataset_path):
        print(f"❌ Dataset not found: {dataset_path}")
        return False
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        print("Loading Geneformer...")
        start_time = time.time()
        
        # Import Geneformer
        sys.path.append('/mnt/d/openclaw/intelligence-augmentation/models/geneformer')
        from geneformer import EmbExtractor
        
        # Conservative settings for 6GB GPU
        batch_size = 4  # Small batch for memory efficiency
        
        print(f"Initializing EmbExtractor...")
        extractor = EmbExtractor(
            model_type="Pretrained",
            num_classes=0,
            emb_mode="cell",
            cell_emb_style="mean_pool",
            filter_data=None,
            max_ncells=n_cells,
            emb_layer=-1,
            forward_batch_size=batch_size,
            nproc=1,
            model_version="V1",
            token_dictionary_file=token_dict,
        )
        
        init_time = time.time() - start_time
        print(f"Initialization time: {init_time:.1f}s")
        
        print("Extracting embeddings...")
        embed_start = time.time()
        
        extractor.extract_embs(
            model_directory=model_path,
            input_data_file=dataset_path,
            output_directory=output_dir,
            output_prefix=f"control_{n_cells}",
        )
        
        embed_time = time.time() - embed_start
        total_time = time.time() - start_time
        
        print(f"Embedding time: {embed_time:.1f}s")
        print(f"Total time: {total_time:.1f}s")
        
        # Check output files
        import glob
        output_files = glob.glob(os.path.join(output_dir, "*.pkl"))
        
        if output_files:
            print(f"✓ Success! Generated {len(output_files)} output files:")
            for f in output_files:
                size_mb = os.path.getsize(f) / 1e6
                print(f"  {os.path.basename(f)}: {size_mb:.1f} MB")
            return True
        else:
            print("✗ No output files generated")
            return False
            
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

def alternative_control_generation(n_cells=1000):
    """Alternative approach using InSilicoPerturber with empty gene list"""
    print(f"\n{'='*60}")
    print(f"Alternative: Using InSilicoPerturber with no genes (control)")
    print(f"{'='*60}")
    
    # Memory settings for RTX 2060 (6GB)
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:256"
    
    # Clear GPU memory
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    
    # Paths
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    token_dict = "/mnt/d/openclaw/intelligence-augmentation/models/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
    dataset_path = f"/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/brain_{n_cells}.dataset"
    output_dir = f"/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/control_{n_cells}_alt"
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Import Geneformer
        sys.path.append('/mnt/d/openclaw/intelligence-augmentation/models/geneformer')
        from geneformer import InSilicoPerturber
        
        # Use InSilicoPerturber with no genes to perturb (should give us original embeddings)
        print(f"Initializing InSilicoPerturber with no target genes...")
        isp = InSilicoPerturber(
            perturb_type="delete",
            genes_to_perturb=[],  # Empty list - no genes to perturb
            model_type="Pretrained", 
            num_classes=0,
            emb_mode="cell",
            cell_emb_style="mean_pool",
            filter_data=None,
            max_ncells=n_cells,
            emb_layer=-1,
            forward_batch_size=4,
            nproc=1,
            model_version="V1",
            token_dictionary_file=token_dict,
        )
        
        print("Running 'perturbation' with no genes (control)...")
        isp.perturb_data(
            model_directory=model_path,
            input_data_file=dataset_path,
            output_directory=output_dir,
            output_prefix=f"control_{n_cells}",
        )
        
        # Check output files
        import glob
        output_files = glob.glob(os.path.join(output_dir, "*.pickle"))
        
        if output_files:
            print(f"✓ Control generated using alternative method!")
            for f in output_files:
                size_mb = os.path.getsize(f) / 1e6
                print(f"  {os.path.basename(f)}: {size_mb:.1f} MB")
            return True
        else:
            print("✗ Alternative method failed")
            return False
            
    except Exception as e:
        print(f"✗ Alternative error: {e}")
        return False

def manual_embedding_extraction(n_cells=1000):
    """Manual embedding extraction using Geneformer model directly"""
    print(f"\n{'='*60}")
    print(f"Manual embedding extraction for {n_cells:,} cells")
    print(f"{'='*60}")
    
    try:
        # Import required modules
        sys.path.append('/mnt/d/openclaw/intelligence-augmentation/models/geneformer')
        from transformers import BertForMaskedLM, BertConfig
        import datasets
        
        # Load tokenized dataset
        dataset_path = f"/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/brain_{n_cells}.dataset"
        print(f"Loading dataset: {dataset_path}")
        
        dataset = datasets.load_from_disk(dataset_path)
        print(f"Dataset loaded: {len(dataset)} cells")
        
        # Load model
        model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
        print(f"Loading model: {model_path}")
        
        config = BertConfig.from_pretrained(model_path)
        model = BertForMaskedLM.from_pretrained(model_path, config=config)
        
        if torch.cuda.is_available():
            model = model.cuda()
            print("Model moved to GPU")
        
        model.eval()
        
        # Extract embeddings
        embeddings = []
        batch_size = 4
        
        print("Extracting embeddings...")
        with torch.no_grad():
            for i in range(0, len(dataset), batch_size):
                batch_end = min(i + batch_size, len(dataset))
                batch_data = dataset[i:batch_end]
                
                # Convert to tensors
                input_ids = torch.tensor([item['input_ids'] for item in batch_data])
                attention_mask = torch.tensor([item['attention_mask'] for item in batch_data])
                
                if torch.cuda.is_available():
                    input_ids = input_ids.cuda()
                    attention_mask = attention_mask.cuda()
                
                # Get hidden states
                outputs = model(input_ids=input_ids, attention_mask=attention_mask, output_hidden_states=True)
                hidden_states = outputs.hidden_states[-1]  # Last layer
                
                # Mean pooling
                for j, mask in enumerate(attention_mask):
                    valid_tokens = hidden_states[j][mask.bool()]
                    cell_embedding = valid_tokens.mean(dim=0).cpu().numpy()
                    embeddings.append(cell_embedding)
                
                if (i // batch_size + 1) % 10 == 0:
                    print(f"  Processed {i + batch_size}/{len(dataset)} cells...")
        
        embeddings = np.array(embeddings)
        print(f"Generated embeddings shape: {embeddings.shape}")
        
        # Save embeddings
        output_dir = f"/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/control_{n_cells}_manual"
        os.makedirs(output_dir, exist_ok=True)
        
        output_file = os.path.join(output_dir, f"control_{n_cells}_embeddings.pkl")
        with open(output_file, 'wb') as f:
            pickle.dump({
                'embeddings': embeddings,
                'n_cells': len(embeddings),
                'embedding_dim': embeddings.shape[1] if len(embeddings) > 0 else 0,
                'method': 'manual_extraction'
            }, f)
        
        print(f"✓ Control embeddings saved: {output_file}")
        return True
        
    except Exception as e:
        print(f"✗ Manual extraction error: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Generate control embeddings using multiple approaches"""
    print("Generating control embeddings for 1000-cell dataset")
    print("="*60)
    
    n_cells = 1000
    success = False
    
    # Try EmbExtractor first
    print("\n🔹 Method 1: EmbExtractor")
    try:
        success = generate_control_embeddings(n_cells)
    except Exception as e:
        print(f"EmbExtractor failed: {e}")
    
    # If that fails, try InSilicoPerturber with empty gene list  
    if not success:
        print("\n🔹 Method 2: InSilicoPerturber with empty gene list")
        try:
            success = alternative_control_generation(n_cells)
        except Exception as e:
            print(f"Alternative method failed: {e}")
    
    # If both fail, try manual extraction
    if not success:
        print("\n🔹 Method 3: Manual embedding extraction")
        try:
            success = manual_embedding_extraction(n_cells)
        except Exception as e:
            print(f"Manual extraction failed: {e}")
    
    if success:
        print(f"\n✅ Control embeddings generated successfully!")
    else:
        print(f"\n❌ All methods failed to generate control embeddings")
    
    return success

if __name__ == "__main__":
    main()