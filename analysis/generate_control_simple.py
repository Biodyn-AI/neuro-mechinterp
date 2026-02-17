#!/usr/bin/env python3
"""Simple control embedding generation for 1000-cell dataset"""
import os
import sys
import torch
import pickle
import numpy as np
from pathlib import Path

def generate_control_embeddings_simple():
    """Generate control embeddings using direct Geneformer model approach"""
    print("🔹 Generating control embeddings (simple approach)")
    print("=" * 60)
    
    # Set memory constraints
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:256"
    
    try:
        # Import required modules
        sys.path.append('/mnt/d/openclaw/intelligence-augmentation/models/geneformer')
        from transformers import BertForMaskedLM, BertConfig
        import datasets
        import torch
        
        # Clear GPU memory
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            print(f"GPU: {torch.cuda.get_device_name(0)}")
        
        # Paths
        dataset_path = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/brain_1000.dataset"
        model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
        output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/control_1000_simple"
        
        print(f"Dataset: {dataset_path}")
        print(f"Model: {model_path}")
        print(f"Output: {output_dir}")
        
        # Load tokenized dataset
        print("Loading dataset...")
        dataset = datasets.load_from_disk(dataset_path)
        print(f"Dataset loaded: {len(dataset)} cells")
        
        # Show dataset structure
        print("Dataset structure:")
        if len(dataset) > 0:
            sample = dataset[0]
            print(f"  Keys: {sample.keys()}")
            print(f"  Input IDs shape: {len(sample['input_ids']) if 'input_ids' in sample else 'N/A'}")
        
        # Load model
        print("Loading Geneformer model...")
        config = BertConfig.from_pretrained(model_path)
        model = BertForMaskedLM.from_pretrained(model_path, config=config)
        
        if torch.cuda.is_available():
            model = model.cuda()
            print("Model moved to GPU")
        
        model.eval()
        print("Model ready for inference")
        
        # Extract embeddings in batches
        print("Extracting embeddings...")
        embeddings = []
        batch_size = 8  # Conservative for 6GB GPU
        
        with torch.no_grad():
            for i in range(0, len(dataset), batch_size):
                batch_end = min(i + batch_size, len(dataset))
                print(f"  Processing batch {i//batch_size + 1}: cells {i+1}-{batch_end}")
                
                batch_input_ids = []
                batch_attention_masks = []
                
                # Prepare batch
                for j in range(i, batch_end):
                    cell_data = dataset[j]
                    
                    # Handle different possible formats
                    if isinstance(cell_data, dict):
                        input_ids = cell_data['input_ids']
                        attention_mask = cell_data.get('attention_mask')
                    else:
                        # If it's not a dict, assume it's the input_ids directly
                        input_ids = cell_data
                        attention_mask = None
                    
                    # Convert to list if needed
                    if isinstance(input_ids, torch.Tensor):
                        input_ids = input_ids.tolist()
                    
                    # Create attention mask if not provided
                    if attention_mask is None:
                        attention_mask = [1] * len(input_ids)
                    elif isinstance(attention_mask, torch.Tensor):
                        attention_mask = attention_mask.tolist()
                    
                    batch_input_ids.append(input_ids)
                    batch_attention_masks.append(attention_mask)
                
                # Convert to tensors
                # Pad sequences to the same length within batch
                max_len = max(len(seq) for seq in batch_input_ids)
                
                padded_input_ids = []
                padded_attention_masks = []
                
                for input_ids, attention_mask in zip(batch_input_ids, batch_attention_masks):
                    # Pad with zeros (or appropriate pad token)
                    padding_length = max_len - len(input_ids)
                    padded_input_ids.append(input_ids + [0] * padding_length)
                    padded_attention_masks.append(attention_mask + [0] * padding_length)
                
                # Convert to tensors
                input_tensor = torch.tensor(padded_input_ids, dtype=torch.long)
                attention_tensor = torch.tensor(padded_attention_masks, dtype=torch.long)
                
                if torch.cuda.is_available():
                    input_tensor = input_tensor.cuda()
                    attention_tensor = attention_tensor.cuda()
                
                # Get embeddings
                outputs = model(input_ids=input_tensor, attention_mask=attention_tensor, output_hidden_states=True)
                hidden_states = outputs.hidden_states[-1]  # Last layer
                
                # Mean pooling over valid tokens for each cell
                for k, mask in enumerate(attention_tensor):
                    valid_positions = mask.bool()
                    if valid_positions.sum() > 0:
                        cell_embedding = hidden_states[k][valid_positions].mean(dim=0)
                        embeddings.append(cell_embedding.cpu().numpy())
                    else:
                        # Fallback: use the first token embedding if no valid tokens
                        cell_embedding = hidden_states[k][0]
                        embeddings.append(cell_embedding.cpu().numpy())
                
                # Memory cleanup
                del outputs, hidden_states, input_tensor, attention_tensor
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
        
        # Convert to numpy array
        embeddings_array = np.array(embeddings)
        print(f"Generated embeddings: {embeddings_array.shape}")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Save embeddings
        output_file = os.path.join(output_dir, "control_1000_embeddings.pkl")
        with open(output_file, 'wb') as f:
            pickle.dump({
                'embeddings': embeddings_array,
                'n_cells': len(embeddings_array),
                'embedding_dim': embeddings_array.shape[1] if len(embeddings_array) > 0 else 0,
                'method': 'direct_extraction'
            }, f)
        
        print(f"✅ Control embeddings saved: {output_file}")
        print(f"   Shape: {embeddings_array.shape}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = generate_control_embeddings_simple()
    if success:
        print("\n🎉 Control embedding generation completed successfully!")
    else:
        print("\n💥 Control embedding generation failed!")