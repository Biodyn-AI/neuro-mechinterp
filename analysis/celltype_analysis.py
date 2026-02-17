#!/usr/bin/env python3
"""
Cell-type Stratified Analysis of Intelligence Gene Perturbations
================================================================

This script performs comprehensive cell-type stratified re-analysis of the 
single-gene perturbation results, as recommended in FINAL_REPORT.md.

Key analyses:
1. Cell-type specific gene rankings comparison
2. Consistency analysis across cell types
3. Cell-type specific pathway analysis
4. Statistical significance testing
5. Visualization generation
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr, pearsonr
from scipy.cluster.hierarchy import dendrogram, linkage
import pickle
import os
from pathlib import Path

# Set up plotting style
plt.style.use('default')
sns.set_palette("husl")

def load_celltype_data():
    """Load cell-type stratified perturbation results"""
    with open('results/insilico_wsl/celltype_summary.json', 'r') as f:
        celltype_data = json.load(f)
    
    # Load original rankings for comparison
    with open('results/insilico_wsl/final_rankings.json', 'r') as f:
        original_rankings = json.load(f)
    
    return celltype_data, original_rankings

def analyze_gene_rankings(celltype_data):
    """Analyze and compare gene rankings across cell types"""
    
    # Extract rankings for each cell type
    cell_types = list(celltype_data['rankings_by_type'].keys())
    genes = [item[0] for item in celltype_data['rankings_by_type']['Excitatory']]  # Use excitatory as reference
    
    # Create ranking matrix
    ranking_matrix = []
    effect_matrix = []
    
    for cell_type in cell_types:
        rankings = celltype_data['rankings_by_type'][cell_type]
        
        # Create dict for easy lookup
        rank_dict = {gene: i+1 for i, (gene, effect, count) in enumerate(rankings)}
        effect_dict = {gene: effect for gene, effect, count in rankings}
        
        # Get rankings for all genes (fill with max+1 for missing genes)
        max_rank = len(rankings) + 1
        cell_rankings = [rank_dict.get(gene, max_rank) for gene in genes]
        cell_effects = [effect_dict.get(gene, 0) for gene in genes]
        
        ranking_matrix.append(cell_rankings)
        effect_matrix.append(cell_effects)
    
    ranking_df = pd.DataFrame(ranking_matrix, 
                             index=cell_types, 
                             columns=genes)
    
    effect_df = pd.DataFrame(effect_matrix,
                            index=cell_types,
                            columns=genes)
    
    return ranking_df, effect_df, genes, cell_types

def compute_ranking_correlations(ranking_df):
    """Compute pairwise correlations between cell type rankings"""
    
    # Spearman correlation (rank-based)
    spearman_corr = ranking_df.T.corr(method='spearman')
    
    # Pearson correlation 
    pearson_corr = ranking_df.T.corr(method='pearson')
    
    return spearman_corr, pearson_corr

def identify_consistent_genes(ranking_df, top_n=5):
    """Identify genes that rank consistently high across cell types"""
    
    # Count how often each gene appears in top N for each cell type
    consistency_scores = {}
    
    for gene in ranking_df.columns:
        # Count appearances in top N across cell types
        top_appearances = sum(ranking_df[gene] <= top_n)
        consistency_scores[gene] = top_appearances
    
    # Sort by consistency
    consistent_genes = sorted(consistency_scores.items(), 
                            key=lambda x: x[1], 
                            reverse=True)
    
    return consistent_genes, consistency_scores

def analyze_cell_type_specificity(ranking_df, effect_df):
    """Identify cell-type specific effects"""
    
    specificity_analysis = {}
    
    for cell_type in ranking_df.index:
        # Find genes that rank much better in this cell type
        cell_ranks = ranking_df.loc[cell_type]
        other_ranks = ranking_df.drop(cell_type).mean()
        
        # Calculate rank difference (negative = better in this cell type)
        rank_diff = cell_ranks - other_ranks
        
        # Get top cell-type specific genes
        specific_genes = rank_diff.nsmallest(5).index.tolist()
        
        specificity_analysis[cell_type] = {
            'specific_genes': specific_genes,
            'rank_differences': rank_diff.to_dict(),
            'top_effects': effect_df.loc[cell_type].nlargest(5).to_dict()
        }
    
    return specificity_analysis

def pathway_analysis_by_celltype(celltype_data):
    """Analyze pathway-level effects by cell type"""
    
    # Define pathway groupings based on FINAL_REPORT.md
    pathways = {
        'Cell Adhesion': ['CADM2', 'NEGR1'],
        'Glutamate Receptors': ['GRIN2A', 'GRIN2B'], 
        'Synaptic Scaffolding': ['NRXN1', 'HOMER1', 'NLGN1', 'SHANK3'],
        'Neurotrophic Signaling': ['CAMK2A', 'BDNF', 'CREB1'],
        'Neurodegeneration': ['APP', 'FOXO3', 'MAPT'],
        'Neurodevelopmental/TF': ['MEF2C', 'TCF4', 'FMR1'],
        'Synaptic Vesicle': ['SCN1A', 'SLC6A4', 'SYN1'],
        'Catecholamine': ['COMT']
    }
    
    pathway_effects = {}
    
    for cell_type, rankings in celltype_data['rankings_by_type'].items():
        # Convert to dict for easy lookup
        effect_dict = {gene: effect for gene, effect, count in rankings}
        
        cell_pathway_effects = {}
        for pathway_name, pathway_genes in pathways.items():
            # Get effects for genes in this pathway
            pathway_vals = [effect_dict.get(gene, 0) for gene in pathway_genes]
            pathway_vals = [v for v in pathway_vals if v > 0]  # Remove missing
            
            if pathway_vals:
                cell_pathway_effects[pathway_name] = {
                    'mean_effect': np.mean(pathway_vals),
                    'max_effect': np.max(pathway_vals),
                    'gene_count': len(pathway_vals),
                    'genes': pathway_genes
                }
        
        pathway_effects[cell_type] = cell_pathway_effects
    
    return pathway_effects

def generate_visualizations(ranking_df, effect_df, spearman_corr, 
                          consistent_genes, specificity_analysis, 
                          pathway_effects):
    """Generate comprehensive visualizations"""
    
    # Create figures directory
    os.makedirs('figures/celltype', exist_ok=True)
    
    # 1. Heatmap of gene rankings across cell types
    plt.figure(figsize=(14, 8))
    sns.heatmap(ranking_df.iloc[:, :15],  # Top 15 genes only for readability
                annot=True, fmt='d', cmap='RdYlBu_r', 
                cbar_kws={'label': 'Rank (lower is better)'})
    plt.title('Gene Rankings Across Cell Types\n(Top 15 Genes by Overall Performance)')
    plt.xlabel('Intelligence Genes')
    plt.ylabel('Cell Types')
    plt.tight_layout()
    plt.savefig('figures/celltype/fig1_ranking_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Effect size heatmap
    plt.figure(figsize=(14, 8))
    sns.heatmap(effect_df.iloc[:, :15], 
                annot=True, fmt='.4f', cmap='viridis',
                cbar_kws={'label': 'Embedding Shift (cosine distance)'})
    plt.title('Gene Perturbation Effects Across Cell Types\n(Top 15 Genes)')
    plt.xlabel('Intelligence Genes')
    plt.ylabel('Cell Types')
    plt.tight_layout()
    plt.savefig('figures/celltype/fig2_effects_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Cell type correlation matrix
    plt.figure(figsize=(8, 6))
    sns.heatmap(spearman_corr, annot=True, fmt='.3f', cmap='coolwarm',
                center=0, square=True)
    plt.title('Spearman Correlation Between Cell Type Rankings')
    plt.tight_layout()
    plt.savefig('figures/celltype/fig3_celltype_correlation.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Consistency score bar plot
    genes_list, scores_list = zip(*consistent_genes[:15])
    plt.figure(figsize=(12, 6))
    bars = plt.bar(range(len(genes_list)), scores_list, color='steelblue')
    plt.xlabel('Intelligence Genes')
    plt.ylabel('Consistency Score (Top-5 appearances across cell types)')
    plt.title('Gene Ranking Consistency Across Cell Types')
    plt.xticks(range(len(genes_list)), genes_list, rotation=45, ha='right')
    
    # Add value labels on bars
    for i, bar in enumerate(bars):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                f'{int(height)}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig('figures/celltype/fig4_consistency_scores.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. Cell-type specific top genes
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.ravel()
    
    for i, (cell_type, data) in enumerate(specificity_analysis.items()):
        if i < len(axes):
            specific_genes_list = data['specific_genes']
            effects = [data['top_effects'].get(gene, 0) for gene in specific_genes_list]
            
            axes[i].bar(range(len(specific_genes_list)), effects, color=f'C{i}')
            axes[i].set_title(f'{cell_type} - Top Specific Genes')
            axes[i].set_xlabel('Genes')
            axes[i].set_ylabel('Embedding Shift')
            axes[i].set_xticks(range(len(specific_genes_list)))
            axes[i].set_xticklabels(specific_genes_list, rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig('figures/celltype/fig5_celltype_specificity.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 6. Pathway analysis across cell types
    # Prepare data for pathway heatmap
    pathway_names = list(next(iter(pathway_effects.values())).keys())
    cell_types = list(pathway_effects.keys())
    
    pathway_matrix = []
    for cell_type in cell_types:
        row = []
        for pathway in pathway_names:
            if pathway in pathway_effects[cell_type]:
                row.append(pathway_effects[cell_type][pathway]['mean_effect'])
            else:
                row.append(0)
        pathway_matrix.append(row)
    
    plt.figure(figsize=(12, 8))
    sns.heatmap(pathway_matrix, 
                xticklabels=pathway_names,
                yticklabels=cell_types,
                annot=True, fmt='.4f', cmap='plasma',
                cbar_kws={'label': 'Mean Pathway Effect'})
    plt.title('Pathway-Level Effects Across Cell Types')
    plt.xlabel('Functional Pathways')
    plt.ylabel('Cell Types')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('figures/celltype/fig6_pathway_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

def statistical_analysis(ranking_df, effect_df):
    """Perform statistical analyses"""
    
    stats_results = {}
    
    # 1. Test for significant differences in rankings between cell types
    # Using Friedman test (non-parametric ANOVA for rankings)
    gene_rankings = ranking_df.values
    friedman_stat, friedman_p = stats.friedmanchisquare(*gene_rankings)
    
    stats_results['friedman_test'] = {
        'statistic': friedman_stat,
        'p_value': friedman_p,
        'interpretation': 'Significant differences in rankings across cell types' if friedman_p < 0.05 else 'No significant differences'
    }
    
    # 2. Coefficient of variation for each gene across cell types
    cv_scores = {}
    for gene in ranking_df.columns:
        ranks = ranking_df[gene].values
        cv = np.std(ranks) / np.mean(ranks) if np.mean(ranks) > 0 else 0
        cv_scores[gene] = cv
    
    # Most and least variable genes
    most_variable = sorted(cv_scores.items(), key=lambda x: x[1], reverse=True)[:5]
    least_variable = sorted(cv_scores.items(), key=lambda x: x[1])[:5]
    
    stats_results['ranking_variability'] = {
        'most_variable_genes': most_variable,
        'least_variable_genes': least_variable,
        'all_cv_scores': cv_scores
    }
    
    # 3. Cell type similarity analysis
    cell_similarities = {}
    cell_types = ranking_df.index.tolist()
    
    for i, ct1 in enumerate(cell_types):
        for j, ct2 in enumerate(cell_types[i+1:], i+1):
            rho, p_val = spearmanr(ranking_df.loc[ct1], ranking_df.loc[ct2])
            cell_similarities[f"{ct1}_vs_{ct2}"] = {
                'spearman_rho': rho,
                'p_value': p_val
            }
    
    stats_results['cell_type_similarities'] = cell_similarities
    
    return stats_results

def main():
    """Main analysis pipeline"""
    
    print("🧠 Starting Cell-Type Stratified Analysis...")
    
    # Load data
    print("📊 Loading cell-type data...")
    celltype_data, original_rankings = load_celltype_data()
    
    # Analyze rankings
    print("🔍 Analyzing gene rankings across cell types...")
    ranking_df, effect_df, genes, cell_types = analyze_gene_rankings(celltype_data)
    
    # Compute correlations
    print("📈 Computing ranking correlations...")
    spearman_corr, pearson_corr = compute_ranking_correlations(ranking_df)
    
    # Identify consistent genes
    print("🎯 Identifying consistently high-ranking genes...")
    consistent_genes, consistency_scores = identify_consistent_genes(ranking_df)
    
    # Analyze specificity
    print("🔬 Analyzing cell-type specificity...")
    specificity_analysis = analyze_cell_type_specificity(ranking_df, effect_df)
    
    # Pathway analysis
    print("🛤️ Performing pathway analysis...")
    pathway_effects = pathway_analysis_by_celltype(celltype_data)
    
    # Statistical analysis
    print("📊 Running statistical tests...")
    stats_results = statistical_analysis(ranking_df, effect_df)
    
    # Generate visualizations
    print("📊 Generating visualizations...")
    generate_visualizations(ranking_df, effect_df, spearman_corr, 
                           consistent_genes, specificity_analysis, 
                           pathway_effects)
    
    # Save results
    results = {
        'ranking_df': ranking_df.to_dict(),
        'effect_df': effect_df.to_dict(),
        'spearman_correlations': spearman_corr.to_dict(),
        'consistent_genes': consistent_genes,
        'specificity_analysis': specificity_analysis,
        'pathway_effects': pathway_effects,
        'statistical_results': stats_results,
        'summary': {
            'n_cell_types': len(cell_types),
            'n_genes_analyzed': len(genes),
            'top_consistent_gene': consistent_genes[0][0],
            'most_correlated_cell_types': spearman_corr.unstack().nlargest(2).index.tolist()[1]  # Exclude self-correlation
        }
    }
    
    with open('results/insilico_wsl/celltype_analysis_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print("✅ Analysis complete! Results saved.")
    print(f"📁 Figures saved to: figures/celltype/")
    print(f"📄 Results saved to: results/insilico_wsl/celltype_analysis_results.json")
    
    return results

if __name__ == "__main__":
    # Change to analysis directory
    os.chdir('/mnt/d/openclaw/intelligence-augmentation/analysis')
    results = main()