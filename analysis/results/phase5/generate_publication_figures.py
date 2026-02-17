"""
Phase 5: Publication-Quality Figures Generation
Creates professional figures from 1000-cell perturbation analysis results
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

INPUT_DIR = r"D:\openclaw\intelligence-augmentation\analysis\results\phase5"
OUTPUT_DIR = os.path.join(INPUT_DIR, "figures")

def load_results():
    """Load analysis results"""
    # Gene-level statistics
    gene_stats_path = os.path.join(INPUT_DIR, 'gene_statistics_1000cells.csv')
    gene_stats = pd.read_csv(gene_stats_path)
    
    # Cell-type specific results
    celltype_stats_path = os.path.join(INPUT_DIR, 'gene_celltype_stats_1000cells.csv')
    celltype_stats = pd.read_csv(celltype_stats_path)
    
    # Raw cell-level results
    cell_results_path = os.path.join(INPUT_DIR, 'cell_level_results_1000.csv')
    cell_results = pd.read_csv(cell_results_path)
    
    print(f"Loaded results: {len(gene_stats)} genes, {len(celltype_stats)} gene-celltype pairs")
    return gene_stats, celltype_stats, cell_results

def create_figure_a_heatmap(gene_stats, celltype_stats):
    """
    Figure A: Heatmap of gene x cell_type mean perturbation effect with clustering
    """
    # Create pivot table for heatmap
    pivot_data = celltype_stats.pivot(index='gene', columns='cell_type', values='mean_distance')
    
    # Fill NaN with 0 and filter to genes with sufficient data
    pivot_data = pivot_data.fillna(0)
    gene_counts = celltype_stats.groupby('gene')['n_cells'].sum()
    valid_genes = gene_counts[gene_counts >= 20].index  # Genes with at least 20 cells total
    pivot_data = pivot_data.loc[valid_genes]
    
    # Hierarchical clustering
    gene_linkage = linkage(pivot_data.values, method='ward', metric='euclidean')
    celltype_linkage = linkage(pivot_data.T.values, method='ward', metric='euclidean')
    
    # Create clustered heatmap
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Use seaborn clustermap for hierarchical clustering
    g = sns.clustermap(pivot_data, 
                       cmap='RdBu_r', center=0,
                       figsize=(12, 10),
                       cbar_kws={'label': 'Mean Perturbation Effect'},
                       xticklabels=True, yticklabels=True,
                       method='ward')
    
    g.fig.suptitle('Gene × Cell Type Perturbation Effects\n(1000 cells, hierarchically clustered)', 
                   fontsize=16, fontweight='bold', y=1.02)
    
    # Rotate x-axis labels
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    
    # Save
    output_path = os.path.join(OUTPUT_DIR, 'figure_a_heatmap')
    g.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
    g.savefig(f'{output_path}.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Figure A (Heatmap) saved")

def create_figure_b_ranked_bar(gene_stats):
    """
    Figure B: Ranked bar chart of genes by effect size with 95% CI and significance stars
    """
    # Sort by effect size and take top genes
    top_genes = gene_stats.sort_values('effect_size', ascending=False).head(15)
    
    # Calculate 95% CI (using std_cosine_distance as proxy for SE)
    top_genes = top_genes.copy()
    top_genes['se'] = top_genes['std_cosine_distance'] / np.sqrt(top_genes['n_cells'])
    top_genes['ci_lower'] = top_genes['effect_size'] - 1.96 * top_genes['se']
    top_genes['ci_upper'] = top_genes['effect_size'] + 1.96 * top_genes['se']
    
    # Significance levels
    def get_sig_stars(p_val):
        if pd.isna(p_val):
            return ''
        elif p_val < 0.001:
            return '***'
        elif p_val < 0.01:
            return '**'
        elif p_val < 0.05:
            return '*'
        else:
            return ''
    
    top_genes['sig_stars'] = top_genes['median_empirical_p'].apply(get_sig_stars)
    
    # Create bar plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Colors based on significance
    colors = ['#d73027' if p < 0.05 else '#4575b4' 
              for p in top_genes['median_empirical_p'].fillna(1.0)]
    
    bars = ax.bar(range(len(top_genes)), top_genes['effect_size'], 
                  yerr=[top_genes['effect_size'] - top_genes['ci_lower'],
                        top_genes['ci_upper'] - top_genes['effect_size']],
                  capsize=5, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # Add significance stars
    for i, (bar, stars) in enumerate(zip(bars, top_genes['sig_stars'])):
        if stars:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + top_genes.iloc[i]['se'] * 2,
                   stars, ha='center', va='bottom', fontsize=14, fontweight='bold')
    
    # Formatting
    ax.set_xlabel('Gene', fontsize=14, fontweight='bold')
    ax.set_ylabel('Perturbation Effect Size', fontsize=14, fontweight='bold')
    ax.set_title('Ranked Gene Perturbation Effects\n(1000 cells, with 95% confidence intervals)', 
                 fontsize=16, fontweight='bold', pad=20)
    
    ax.set_xticks(range(len(top_genes)))
    ax.set_xticklabels(top_genes['gene'], rotation=45, ha='right')
    
    # Add horizontal line at 0
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#d73027', label='p < 0.05'),
                      Patch(facecolor='#4575b4', label='p ≥ 0.05')]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    
    # Save
    output_path = os.path.join(OUTPUT_DIR, 'figure_b_ranked_bars')
    plt.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_path}.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Figure B (Ranked bars) saved")

def create_figure_c_volcano(gene_stats):
    """
    Figure C: Volcano plot of effect size vs -log10(FDR p-value)
    """
    # Prepare data
    plot_data = gene_stats.copy()
    plot_data['neg_log10_fdr'] = -np.log10(plot_data['fdr'].fillna(1.0))
    
    # Significance thresholds
    fdr_threshold = 0.1  # 10% FDR
    effect_threshold = 0.001  # Minimum effect size
    
    # Color coding
    plot_data['significant'] = (plot_data['fdr'] < fdr_threshold) & (plot_data['effect_size'] > effect_threshold)
    colors = plot_data['significant'].map({True: '#d73027', False: '#4575b4'})
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 9))
    
    scatter = ax.scatter(plot_data['effect_size'], plot_data['neg_log10_fdr'], 
                        c=colors, alpha=0.7, s=60, edgecolor='black', linewidth=0.5)
    
    # Add threshold lines
    ax.axhline(y=-np.log10(fdr_threshold), color='red', linestyle='--', alpha=0.6, 
               label=f'FDR = {fdr_threshold}')
    ax.axvline(x=effect_threshold, color='red', linestyle='--', alpha=0.6,
               label=f'Effect = {effect_threshold}')
    
    # Label top genes
    top_genes = plot_data.nlargest(10, 'effect_size')
    for _, gene_data in top_genes.iterrows():
        ax.annotate(gene_data['gene'], 
                   (gene_data['effect_size'], gene_data['neg_log10_fdr']),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=10, ha='left', va='bottom',
                   bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.7))
    
    # Formatting
    ax.set_xlabel('Perturbation Effect Size', fontsize=14, fontweight='bold')
    ax.set_ylabel('-log₁₀(FDR p-value)', fontsize=14, fontweight='bold')
    ax.set_title('Volcano Plot: Effect Size vs Statistical Significance\n(1000 cells)', 
                 fontsize=16, fontweight='bold', pad=20)
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#d73027', label=f'Significant (FDR < {fdr_threshold})'),
                      Patch(facecolor='#4575b4', label=f'Not significant (FDR ≥ {fdr_threshold})')]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    
    # Save
    output_path = os.path.join(OUTPUT_DIR, 'figure_c_volcano')
    plt.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_path}.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Figure C (Volcano plot) saved")

def create_figure_d_top_genes_celltype(gene_stats, celltype_stats):
    """
    Figure D: Top 5 genes cell-type profiles (grouped bar chart)
    """
    # Get top 5 genes by effect size
    top_5_genes = gene_stats.nlargest(5, 'effect_size')['gene'].tolist()
    
    # Filter cell-type data for top 5 genes
    top_data = celltype_stats[celltype_stats['gene'].isin(top_5_genes)].copy()
    
    # Get common cell types (appear for multiple genes)
    celltype_counts = top_data['cell_type'].value_counts()
    common_celltypes = celltype_counts[celltype_counts >= 3].index.tolist()
    top_data = top_data[top_data['cell_type'].isin(common_celltypes)]
    
    # Create grouped bar chart
    fig, ax = plt.subplots(figsize=(15, 8))
    
    # Prepare data for plotting
    plot_data = top_data.pivot(index='cell_type', columns='gene', values='mean_distance')
    plot_data = plot_data.fillna(0)
    
    # Create grouped bar plot
    plot_data.plot(kind='bar', ax=ax, width=0.8, 
                   colormap='Set2', alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # Formatting
    ax.set_xlabel('Cell Type', fontsize=14, fontweight='bold')
    ax.set_ylabel('Mean Perturbation Effect', fontsize=14, fontweight='bold')
    ax.set_title('Cell-Type Specific Effects for Top 5 Genes\n(1000 cells)', 
                 fontsize=16, fontweight='bold', pad=20)
    
    ax.legend(title='Gene', bbox_to_anchor=(1.05, 1), loc='upper left', 
              title_fontsize=12, fontsize=11)
    
    # Rotate x-axis labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    
    # Add horizontal line at 0
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    
    # Save
    output_path = os.path.join(OUTPUT_DIR, 'figure_d_celltype_profiles')
    plt.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_path}.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Figure D (Cell-type profiles) saved")

def create_summary_statistics_table(gene_stats):
    """Create a summary table for publication"""
    summary = gene_stats.head(10)[['gene', 'n_cells', 'mean_cosine_distance', 'effect_size', 
                                   'median_empirical_p', 'fdr', 'pct_sig']].copy()
    
    # Format for publication
    summary['effect_size'] = summary['effect_size'].round(4)
    summary['mean_cosine_distance'] = summary['mean_cosine_distance'].round(4)
    summary['median_empirical_p'] = summary['median_empirical_p'].round(4)
    summary['fdr'] = summary['fdr'].round(4)
    summary['pct_sig'] = summary['pct_sig'].round(1)
    
    summary.to_csv(os.path.join(OUTPUT_DIR, 'summary_statistics_table.csv'), index=False)
    print("Summary statistics table saved")

def main():
    """Generate all publication figures"""
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Figures will be saved to: {OUTPUT_DIR}")
    
    # Load results
    gene_stats, celltype_stats, cell_results = load_results()
    
    # Generate figures
    print("\nGenerating Figure A: Heatmap...")
    create_figure_a_heatmap(gene_stats, celltype_stats)
    
    print("\nGenerating Figure B: Ranked bar chart...")
    create_figure_b_ranked_bar(gene_stats)
    
    print("\nGenerating Figure C: Volcano plot...")
    create_figure_c_volcano(gene_stats)
    
    print("\nGenerating Figure D: Cell-type profiles...")
    create_figure_d_top_genes_celltype(gene_stats, celltype_stats)
    
    print("\nGenerating summary table...")
    create_summary_statistics_table(gene_stats)
    
    print(f"\nAll figures saved to: {OUTPUT_DIR}")
    print("TASK 2 COMPLETE: Publication-quality figures generated")

if __name__ == '__main__':
    main()