#!/usr/bin/env python3
"""
Run 20 control gene pair combo perturbations through Geneformer,
then compare super-additivity rates with intelligence combos.
"""
import os, sys, pickle, json, torch, gc, random, itertools
import numpy as np
from pathlib import Path

def get_control_singles(results_dir):
    """Load all control single-gene perturbation results."""
    import glob
    controls = {}
    for d in sorted(glob.glob(os.path.join(results_dir, 'control_perturb_*'))):
        gene = os.path.basename(d).replace('control_perturb_', '')
        pkls = glob.glob(os.path.join(d, '*.pickle'))
        if not pkls: continue
        with open(pkls[0], 'rb') as f:
            data = pickle.load(f)
        for k, v in data.items():
            if isinstance(k, tuple) and 'cell_emb' in str(k):
                arr = np.array(v)
                shift = 1.0 - np.mean(arr)
                token_id = k[0]
                if len(arr) >= 20:  # need enough cells
                    controls[gene] = {'shift': shift, 'n': len(arr), 'token_id': int(token_id)}
                    print(f"  {gene}: shift={shift:.4f}, n={len(arr)}, token={token_id}")
                break
    return controls

def get_gene_token_map(controls):
    """Get gene->ensembl mapping from token dict."""
    token_dict_path = os.path.expanduser(
        "~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
    )
    with open(token_dict_path, 'rb') as f:
        token_dict = pickle.load(f)
    # Reverse: token_id -> ensembl
    inv = {v: k for k, v in token_dict.items()}
    
    gene_ensembl = {}
    for gene, info in controls.items():
        tid = info['token_id']
        if tid in inv:
            gene_ensembl[gene] = inv[tid]
            print(f"  {gene}: token {tid} -> {inv[tid]}")
    return gene_ensembl

def main():
    results_dir = '/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl'
    model_path = '/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M'
    tokenized_data = os.path.join(results_dir, 'brain.dataset')
    token_dict_path = os.path.expanduser(
        "~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
    )
    
    print("=== Loading control single-gene results ===")
    controls = get_control_singles(results_dir)
    print(f"\nFound {len(controls)} control genes with >=20 cells")
    
    # Pick 20 random pairs from controls with sufficient cells
    genes = sorted(controls.keys())
    all_pairs = list(itertools.combinations(genes, 2))
    random.seed(42)
    random.shuffle(all_pairs)
    pairs = all_pairs[:20]
    
    print(f"\n=== Will run {len(pairs)} control combo pairs ===")
    for i, (g1, g2) in enumerate(pairs):
        print(f"  {i+1}. {g1} + {g2}")
    
    # Get ensembl IDs
    print("\n=== Resolving Ensembl IDs ===")
    gene_ensembl = get_gene_token_map(controls)
    
    # Run each pair
    from geneformer import InSilicoPerturber
    
    combo_results = {}
    for i, (g1, g2) in enumerate(pairs):
        pair_name = f"{g1}_{g2}"
        pair_output = os.path.join(results_dir, f"control_combo_{pair_name}")
        
        # Skip if already done
        if os.path.exists(pair_output):
            existing_pkls = [f for f in os.listdir(pair_output) if f.endswith('.pickle')]
            if existing_pkls:
                print(f"\n--- Pair {i+1}/20: {g1} + {g2} --- ALREADY DONE, loading")
                with open(os.path.join(pair_output, existing_pkls[0]), 'rb') as f:
                    data = pickle.load(f)
                for k, v in data.items():
                    if isinstance(k, tuple) and 'cell_emb' in str(k):
                        arr = np.array(v)
                        shift = 1.0 - np.mean(arr)
                        expected = controls[g1]['shift'] + controls[g2]['shift']
                        ratio = shift / expected if expected > 0 else 1.0
                        combo_results[pair_name] = {
                            'g1': g1, 'g2': g2, 'shift': shift, 'n': len(arr),
                            'expected': expected, 'ratio': ratio
                        }
                        print(f"  shift={shift:.4f}, exp={expected:.4f}, ratio={ratio:.3f}")
                        break
                continue
        
        if g1 not in gene_ensembl or g2 not in gene_ensembl:
            print(f"\n--- Pair {i+1}/20: {g1} + {g2} --- SKIP (no ensembl ID)")
            continue
        
        print(f"\n--- Pair {i+1}/20: {g1} + {g2} ---")
        os.makedirs(pair_output, exist_ok=True)
        
        torch.cuda.empty_cache()
        gc.collect()
        
        try:
            isp = InSilicoPerturber(
                perturb_type="delete",
                genes_to_perturb=[gene_ensembl[g1], gene_ensembl[g2]],
                model_type="Pretrained",
                num_classes=0,
                emb_mode="cell",
                filter_data=None,
                max_ncells=500,
                emb_layer=-1,
                forward_batch_size=32,
                nproc=1,
                model_version="V1",
                token_dictionary_file=token_dict_path,
            )
            
            isp.perturb_data(
                model_directory=model_path,
                input_data_file=tokenized_data,
                output_directory=pair_output,
                output_prefix=f"control_combo_{pair_name}",
            )
            
            # Extract result
            pkl_files = [f for f in os.listdir(pair_output) if f.endswith('.pickle')]
            for pf in pkl_files:
                with open(os.path.join(pair_output, pf), 'rb') as f:
                    data = pickle.load(f)
                for k, v in data.items():
                    if isinstance(k, tuple) and 'cell_emb' in str(k):
                        arr = np.array(v)
                        shift = 1.0 - np.mean(arr)
                        expected = controls[g1]['shift'] + controls[g2]['shift']
                        ratio = shift / expected if expected > 0 else 1.0
                        combo_results[pair_name] = {
                            'g1': g1, 'g2': g2, 'shift': shift, 'n': len(arr),
                            'expected': expected, 'ratio': ratio
                        }
                        print(f"  RESULT: shift={shift:.4f}, exp={expected:.4f}, ratio={ratio:.3f}, n={len(arr)}")
                        break
            
            print(f"  DONE ({len(combo_results)}/{len(pairs)} complete)")
            
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
    
    # Summary
    print(f"\n\n=== SUMMARY: {len(combo_results)}/{len(pairs)} control combos completed ===")
    ratios = [r['ratio'] for r in combo_results.values()]
    super_additive = sum(1 for r in ratios if r > 1.0)
    print(f"Super-additive (ratio>1): {super_additive}/{len(ratios)}")
    print(f"Mean ratio: {np.mean(ratios):.3f} ± {np.std(ratios):.3f}")
    print(f"Range: [{min(ratios):.3f}, {max(ratios):.3f}]")
    
    # Compare with intelligence combos
    print("\n=== Intelligence combo ratios (from paper) ===")
    intel_ratios = []
    import glob
    for d in sorted(glob.glob(os.path.join(results_dir, 'perturb_combo_*'))):
        name = os.path.basename(d).replace('perturb_combo_', '')
        pkls = glob.glob(os.path.join(d, '*.pickle'))
        if not pkls: continue
        with open(pkls[0], 'rb') as f:
            data = pickle.load(f)
        for k, v in data.items():
            if isinstance(k, tuple) and 'cell_emb' in str(k):
                arr = np.array(v)
                shift = 1.0 - np.mean(arr)
                # Get singles
                genes = name.split('_')
                intel_singles = {
                    'CADM2': 0.01963, 'GRIN2A': 0.01896, 'CAMK2A': 0.01892,
                    'MEF2C': 0.01832, 'APP': 0.01832
                }
                expected = intel_singles.get(genes[0], 0) + intel_singles.get(genes[1], 0)
                ratio = shift / expected if expected > 0 else 1.0
                intel_ratios.append(ratio)
                print(f"  {name}: ratio={ratio:.3f}")
                break
    
    print(f"\nIntelligence combos: mean ratio={np.mean(intel_ratios):.3f}")
    print(f"Control combos: mean ratio={np.mean(ratios):.3f}")
    
    from scipy import stats
    if len(ratios) >= 3 and len(intel_ratios) >= 3:
        u, p = stats.mannwhitneyu(intel_ratios, ratios, alternative='greater')
        print(f"Mann-Whitney U (intel > control): U={u}, p={p:.6f}")
    
    # Save
    output = {
        'control_combos': combo_results,
        'control_ratios': ratios,
        'intel_ratios': intel_ratios,
        'control_mean_ratio': float(np.mean(ratios)),
        'intel_mean_ratio': float(np.mean(intel_ratios)),
        'n_super_additive_control': super_additive,
        'n_super_additive_intel': sum(1 for r in intel_ratios if r > 1.0),
    }
    out_path = '/mnt/d/openclaw/intelligence-augmentation/analysis/CONTROL_COMBO_RESULTS.json'
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved to {out_path}")

if __name__ == '__main__':
    main()
