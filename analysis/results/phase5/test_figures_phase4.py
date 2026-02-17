"""
Test publication figures using Phase 4 data (200 cells)
This will verify the figure generation pipeline while waiting for 1000-cell results
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform
import warnings
warnings.filterwarnings('ignore')

# Set matplotlib parameters for high-quality output
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'axes.linewidth': 1.0,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'legend.frameon': False,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1
})

PHASE4_DIR = r"D:\openclaw\intelligence-augmentation\analysis\results\phase4"
OUTPUT_DIR = r"D:\openclaw\intelligence-augmentation\analysis\results\phase5\figures"

def load_phase4_results():
    """Load Phase 4 results for testing"""
    # Gene-level statistics
    gene_stats = pd.read_csv(os.path.join(PHASE4_DIR, 'gene_statistics_combined.csv'))
    
    # Cell-type specific results
    celltype_stats = pd.read_csv(os.path.join(PHASE4_DIR, 'gene_celltype_stats_phase4.csv'))
    
    # Raw cell-level results
    cell_results = pd.read_csv(os.path.join(PHASE4_DIR, 'cell_level_results_combined.csv'))
    
    print(f"Loaded Phase 4 results: {len(gene_stats)} genes, {len(celltype_stats)} gene-celltype pairs")
    return gene_stats, celltype_stats, cell_results

def create_test_heatmap(gene_stats, celltype_stats):
    """Test heatmap creation"""
    # Create pivot table for heatmap
    pivot_data = celltype_stats.pivot(index='gene', columns='cell_type', values='mean_distance')
    
    # Fill NaN with 0 and filter to genes with sufficient data
    pivot_data = pivot_data.fillna(0)
    gene_counts = celltype_stats.groupby('gene')['n_cells'].sum()
    valid_genes = gene_counts[gene_counts >= 5].index  # Genes with at least 5 cells total
    pivot_data = pivot_data.loc[valid_genes]
    
    # Create clustered heatmap
    try:
        g = sns.clustermap(pivot_data, 
                           cmap='RdBu_r', center=0,
                           figsize=(10, 8),
                           cbar_kws={'label': 'Mean Perturbation Effect'},
                           xticklabels=True, yticklabels=True,
                           method='ward')
        
        g.fig.suptitle('Gene × Cell Type Perturbation Effects\n(Phase 4: 200 cells, test)', 
                       fontsize=14, fontweight='bold', y=1.02)
        
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
        
        output_path = os.path.join(OUTPUT_DIR, 'test_heatmap_phase4')
        g.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Test heatmap saved")
        return True
    except Exception as e:
        print(f"Heatmap creation failed: {e}")
        return False

def create_test_ranked_bars(gene_stats):
    """Test ranked bar chart"""
    try:
        top_genes = gene_stats.sort_values('effect_size', ascending=False).head(10)
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Simple bar plot
        colors = ['#d73027' if p < 0.3 else '#4575b4' 
                  for p in top_genes.get('median_p', top_genes.get('fdr', [1.0]*len(top_genes)))]
        
        bars = ax.bar(range(len(top_genes)), top_genes['effect_size'], 
                      color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
        
        ax.set_xlabel('Gene', fontsize=12, fontweight='bold')
        ax.set_ylabel('Perturbation Effect Size', fontsize=12, fontweight='bold')
        ax.set_title('Ranked Gene Perturbation Effects (Phase 4 Test)', 
                     fontsize=14, fontweight='bold')
        
        ax.set_xticks(range(len(top_genes)))
        ax.set_xticklabels(top_genes['gene'], rotation=45, ha='right')
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        
        plt.tight_layout()
        
        output_path = os.path.join(OUTPUT_DIR, 'test_ranked_bars_phase4')
        plt.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Test ranked bars saved")
        return True
    except Exception as e:
        print(f"Bar chart creation failed: {e}")
        return False

def create_test_volcano(gene_stats):
    """Test volcano plot"""
    try:
        plot_data = gene_stats.copy()
        
        # Use available p-value columns
        p_col = 'fdr' if 'fdr' in plot_data.columns else 'median_p'
        if p_col not in plot_data.columns:
            print("No p-value column found for volcano plot")
            return False
        
        plot_data['neg_log10_p'] = -np.log10(plot_data[p_col].fillna(1.0))
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Color by significance
        colors = ['#d73027' if p < 0.3 else '#4575b4' 
                  for p in plot_data[p_col].fillna(1.0)]
        
        scatter = ax.scatter(plot_data['effect_size'], plot_data['neg_log10_p'], 
                            c=colors, alpha=0.7, s=50, edgecolor='black', linewidth=0.5)
        
        # Label top genes
        top_genes = plot_data.nlargest(5, 'effect_size')
        for _, gene_data in top_genes.iterrows():
            ax.annotate(gene_data['gene'], 
                       (gene_data['effect_size'], gene_data['neg_log10_p']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=9, ha='left', va='bottom',
                       bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.7))
        
        ax.set_xlabel('Perturbation Effect Size', fontsize=12, fontweight='bold')
        ax.set_ylabel(f'-log₁₀({p_col})', fontsize=12, fontweight='bold')
        ax.set_title('Volcano Plot: Effect Size vs Significance (Phase 4 Test)', 
                     fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        
        output_path = os.path.join(OUTPUT_DIR, 'test_volcano_phase4')
        plt.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Test volcano plot saved")
        return True
    except Exception as e:
        print(f"Volcano plot creation failed: {e}")
        return False

def main():
    """Test figure generation pipeline"""
    print("Testing publication figure generation with Phase 4 data...")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load results
    gene_stats, celltype_stats, cell_results = load_phase4_results()
    
    # Test figures
    print("\nTesting heatmap...")
    success1 = create_test_heatmap(gene_stats, celltype_stats)
    
    print("\nTesting ranked bars...")
    success2 = create_test_ranked_bars(gene_stats)
    
    print("\nTesting volcano plot...")
    success3 = create_test_volcano(gene_stats)
    
    if all([success1, success2, success3]):
        print("\nAll test figures generated successfully!")
        print("Figure generation pipeline is working correctly.")
    else:
        print("\nSome figures failed - check error messages above.")
    
    print(f"Test figures saved to: {OUTPUT_DIR}")

if __name__ == '__main__':
    main()