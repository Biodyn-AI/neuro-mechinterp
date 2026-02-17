"""
Phase 5: Drug Repurposing Analysis
Cross-reference top perturbation genes with DGIdb for therapeutic targets
Focus on MEF2C, NEGR1, SHANK3, GRIN2A, CAMK2A, HOMER1, NRXN1, GRIN2B, APP, SCN1A
"""
import os, sys
import pandas as pd
import numpy as np
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

INPUT_DIR = r"D:\openclaw\intelligence-augmentation\analysis\results\phase5"
DGIDB_PATH = r"D:\openclaw\intelligence-augmentation\data\drugs\interactions.tsv"
OUTPUT_DIR = INPUT_DIR

# Priority target genes based on Phase 4/5 results and autism relevance
PRIORITY_GENES = [
    'MEF2C', 'NEGR1', 'SHANK3', 'GRIN2A', 'CAMK2A', 
    'HOMER1', 'NRXN1', 'GRIN2B', 'APP', 'SCN1A'
]

def load_perturbation_results():
    """Load Phase 5 gene statistics"""
    gene_stats_path = os.path.join(INPUT_DIR, 'gene_statistics_1000cells.csv')
    if not os.path.exists(gene_stats_path):
        print(f"Warning: {gene_stats_path} not found. Using Phase 4 results.")
        gene_stats_path = r"D:\openclaw\intelligence-augmentation\analysis\results\phase4\gene_statistics_combined.csv"
    
    gene_stats = pd.read_csv(gene_stats_path)
    print(f"Loaded perturbation results for {len(gene_stats)} genes")
    return gene_stats

def load_drug_interactions():
    """Load DGIdb drug-gene interactions"""
    if not os.path.exists(DGIDB_PATH):
        print(f"Error: DGIdb file not found at {DGIDB_PATH}")
        return None
    
    # Read DGIdb interactions
    interactions = pd.read_csv(DGIDB_PATH, sep='\t', low_memory=False)
    print(f"Loaded {len(interactions)} drug-gene interactions from DGIdb")
    print(f"Columns: {list(interactions.columns)}")
    
    # Show sample of data to understand structure
    print("\nSample interactions:")
    print(interactions.head())
    
    return interactions

def find_gene_drugs(gene_name, interactions_df):
    """Find drugs targeting a specific gene"""
    if interactions_df is None:
        return pd.DataFrame()
    
    # Handle different possible column names for genes
    gene_cols = [col for col in interactions_df.columns if 'gene' in col.lower()]
    if not gene_cols:
        print("Warning: No gene columns found in interactions data")
        return pd.DataFrame()
    
    gene_col = gene_cols[0]  # Use first gene column
    
    # Search for gene (case-insensitive)
    gene_matches = interactions_df[interactions_df[gene_col].str.upper() == gene_name.upper()]
    
    if len(gene_matches) == 0:
        # Try alternative gene symbols (common variants)
        alt_symbols = {
            'GRIN2A': ['GRINA', 'NR2A'],
            'GRIN2B': ['GRINB', 'NR2B'],
            'CAMK2A': ['CAMK2', 'CAMKA'],
            'SCN1A': ['SCN1'],
        }
        
        if gene_name in alt_symbols:
            for alt in alt_symbols[gene_name]:
                gene_matches = interactions_df[interactions_df[gene_col].str.upper() == alt.upper()]
                if len(gene_matches) > 0:
                    break
    
    return gene_matches

def categorize_drugs_by_mechanism(drug_data):
    """Categorize drugs by interaction type/mechanism"""
    if len(drug_data) == 0:
        return {}
    
    # Handle different possible column names
    interaction_cols = [col for col in drug_data.columns if 'interaction' in col.lower() or 'type' in col.lower()]
    
    if not interaction_cols:
        return {'unknown': drug_data}
    
    interaction_col = interaction_cols[0]
    
    categories = {
        'agonist': [],
        'antagonist': [],
        'inhibitor': [],
        'modulator': [],
        'other': []
    }
    
    for _, row in drug_data.iterrows():
        interaction_type = str(row[interaction_col]).lower()
        
        if 'agonist' in interaction_type:
            categories['agonist'].append(row)
        elif 'antagonist' in interaction_type or 'blocker' in interaction_type:
            categories['antagonist'].append(row)
        elif 'inhibitor' in interaction_type:
            categories['inhibitor'].append(row)
        elif 'modulator' in interaction_type:
            categories['modulator'].append(row)
        else:
            categories['other'].append(row)
    
    # Convert back to DataFrames
    for cat in categories:
        if categories[cat]:
            categories[cat] = pd.DataFrame(categories[cat])
        else:
            categories[cat] = pd.DataFrame()
    
    return categories

def analyze_drug_repurposing():
    """Main drug repurposing analysis"""
    print("=== DRUG REPURPOSING ANALYSIS ===\n")
    
    # Load data
    gene_stats = load_perturbation_results()
    interactions = load_drug_interactions()
    
    if interactions is None:
        print("Cannot proceed without DGIdb data")
        return
    
    # Get top genes from perturbation analysis
    top_genes = gene_stats.sort_values('effect_size', ascending=False).head(15)
    
    # Focus on priority genes + top perturbation genes
    target_genes = set(PRIORITY_GENES + top_genes['gene'].tolist())
    
    print(f"Analyzing drug targets for {len(target_genes)} genes...")
    
    # Drug discovery results
    all_drug_candidates = []
    gene_drug_summary = []
    
    for gene in target_genes:
        print(f"\n--- Analyzing {gene} ---")
        
        # Get gene stats
        gene_data = gene_stats[gene_stats['gene'] == gene]
        if len(gene_data) > 0:
            effect_size = gene_data['effect_size'].iloc[0]
            p_value = gene_data.get('median_empirical_p', pd.Series([np.nan])).iloc[0]
            fdr = gene_data.get('fdr', pd.Series([np.nan])).iloc[0]
        else:
            effect_size = np.nan
            p_value = np.nan
            fdr = np.nan
        
        # Find drugs targeting this gene
        gene_drugs = find_gene_drugs(gene, interactions)
        
        if len(gene_drugs) == 0:
            print(f"  No drugs found for {gene}")
            gene_drug_summary.append({
                'gene': gene,
                'effect_size': effect_size,
                'p_value': p_value,
                'fdr': fdr,
                'n_drugs': 0,
                'priority_gene': gene in PRIORITY_GENES
            })
            continue
        
        print(f"  Found {len(gene_drugs)} drug interactions")
        
        # Categorize by mechanism
        drug_categories = categorize_drugs_by_mechanism(gene_drugs)
        
        # Extract drug information
        drug_cols = [col for col in gene_drugs.columns if 'drug' in col.lower() or 'compound' in col.lower()]
        if drug_cols:
            drug_col = drug_cols[0]
            drugs = gene_drugs[drug_col].unique()
            
            for drug in drugs:
                drug_rows = gene_drugs[gene_drugs[drug_col] == drug]
                
                # Get interaction types
                interaction_cols = [col for col in drug_rows.columns if 'interaction' in col.lower() or 'type' in col.lower()]
                if interaction_cols:
                    interactions_list = drug_rows[interaction_cols[0]].unique()
                else:
                    interactions_list = ['unknown']
                
                all_drug_candidates.append({
                    'gene': gene,
                    'drug_name': drug,
                    'interaction_types': ', '.join([str(x) for x in interactions_list]),
                    'effect_size': effect_size,
                    'p_value': p_value,
                    'fdr': fdr,
                    'priority_gene': gene in PRIORITY_GENES,
                    'n_interactions': len(drug_rows)
                })
        
        gene_drug_summary.append({
            'gene': gene,
            'effect_size': effect_size,
            'p_value': p_value,
            'fdr': fdr,
            'n_drugs': len(gene_drugs[drug_cols[0]].unique()) if drug_cols else 0,
            'priority_gene': gene in PRIORITY_GENES
        })
    
    # Save results
    if all_drug_candidates:
        drug_candidates_df = pd.DataFrame(all_drug_candidates)
        drug_candidates_df = drug_candidates_df.sort_values(['priority_gene', 'effect_size'], ascending=[False, False])
        
        output_path = os.path.join(OUTPUT_DIR, 'drug_repurposing_candidates.csv')
        drug_candidates_df.to_csv(output_path, index=False)
        print(f"\nSaved {len(drug_candidates_df)} drug candidates to: {output_path}")
        
        # Summary statistics
        print(f"\n=== DRUG REPURPOSING SUMMARY ===")
        print(f"Total drug candidates identified: {len(drug_candidates_df)}")
        print(f"Unique drugs: {drug_candidates_df['drug_name'].nunique()}")
        print(f"Genes with drug targets: {drug_candidates_df['gene'].nunique()}")
        
        # Priority gene drugs
        priority_drugs = drug_candidates_df[drug_candidates_df['priority_gene'] == True]
        print(f"Drugs targeting priority autism genes: {len(priority_drugs)}")
        
        # Top drug candidates
        print(f"\n=== TOP DRUG CANDIDATES ===")
        top_candidates = drug_candidates_df.head(20)
        print(top_candidates[['gene', 'drug_name', 'interaction_types', 'effect_size']].to_string(index=False))
        
    else:
        print("No drug candidates found!")
    
    # Save gene summary
    gene_summary_df = pd.DataFrame(gene_drug_summary)
    gene_summary_path = os.path.join(OUTPUT_DIR, 'gene_drug_target_summary.csv')
    gene_summary_df.to_csv(gene_summary_path, index=False)
    print(f"\nGene-drug summary saved to: {gene_summary_path}")
    
    print("\nTASK 3 COMPLETE: Drug repurposing analysis finished")

def analyze_mechanism_classes():
    """Analyze drug mechanisms for therapeutic insights"""
    print("\n=== MECHANISM ANALYSIS ===")
    
    # Load drug candidates
    candidates_path = os.path.join(OUTPUT_DIR, 'drug_repurposing_candidates.csv')
    if not os.path.exists(candidates_path):
        print("No drug candidates file found")
        return
    
    candidates = pd.read_csv(candidates_path)
    
    # Analyze interaction types
    all_interactions = []
    for interactions_str in candidates['interaction_types'].dropna():
        interactions = [x.strip() for x in str(interactions_str).split(',')]
        all_interactions.extend(interactions)
    
    interaction_counts = pd.Series(all_interactions).value_counts()
    print("\nMost common interaction types:")
    print(interaction_counts.head(10))
    
    # Focus on priority genes
    priority_candidates = candidates[candidates['priority_gene'] == True]
    print(f"\nPriority gene drug targets ({len(priority_candidates)} candidates):")
    for gene in PRIORITY_GENES:
        gene_drugs = priority_candidates[priority_candidates['gene'] == gene]
        if len(gene_drugs) > 0:
            print(f"{gene}: {gene_drugs['drug_name'].nunique()} drugs")

if __name__ == '__main__':
    analyze_drug_repurposing()
    analyze_mechanism_classes()