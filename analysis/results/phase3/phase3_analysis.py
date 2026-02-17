#!/usr/bin/env python3
"""
Phase 3: Statistical Framework and Publication Figures
Virtual Perturbation Screening Pipeline

This script performs comprehensive statistical analysis and generates
publication-quality figures for the perturbation screening results.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import false_discovery_control
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality matplotlib parameters
plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1
})

# Color palette
colors = sns.color_palette("viridis", 8)
sns.set_palette(colors)

class Phase3Analysis:
    def __init__(self, data_dir):
        self.data_dir = data_dir
        self.output_dir = data_dir + r"\phase3"
        self.load_data()
        
    def load_data(self):
        """Load all Phase 2 results"""
        print("Loading Phase 2 data...")
        
        # Main results
        self.perturbation_results = pd.read_csv(
            f"{self.data_dir}/phase2/perturbation_results_phase2_final.csv"
        )
        
        # Summary by cell type
        self.cell_type_summary = pd.read_csv(
            f"{self.data_dir}/phase2/summary_by_celltype_phase2_final.csv"
        )
        
        # Drug candidates
        self.drug_candidates = pd.read_csv(
            f"{self.data_dir}/phase2/drug_candidates_phase2_final.csv"
        )
        
        # Gene ranking
        self.gene_ranking = pd.read_csv(
            f"{self.data_dir}/phase2/gene_ranking_by_mean_score_phase2.csv"
        )
        
        print(f"Loaded data: {len(self.perturbation_results)} cell-gene pairs")
        print(f"Genes: {self.perturbation_results['gene'].nunique()}")
        print(f"Cell types: {self.perturbation_results['cell_type'].nunique()}")
        
    def permutation_test(self, scores, n_permutations=1000):
        """Perform permutation test for empirical p-values"""
        observed_mean = np.mean(scores)
        
        # Generate null distribution by shuffling
        null_distribution = []
        for _ in range(n_permutations):
            shuffled = np.random.permutation(scores)
            null_distribution.append(np.mean(shuffled))
        
        null_distribution = np.array(null_distribution)
        
        # Calculate empirical p-value (two-tailed)
        p_value = np.mean(np.abs(null_distribution) >= np.abs(observed_mean))
        
        return p_value, null_distribution
    
    def cohens_d(self, scores, population_scores):
        """Calculate Cohen's d effect size"""
        mean_diff = np.mean(scores) - np.mean(population_scores)
        pooled_std = np.sqrt(
            ((len(scores) - 1) * np.var(scores, ddof=1) + 
             (len(population_scores) - 1) * np.var(population_scores, ddof=1)) /
            (len(scores) + len(population_scores) - 2)
        )
        
        if pooled_std == 0:
            return 0
        
        return mean_diff / pooled_std
    
    def cell_type_specificity_index(self, gene_data):
        """Calculate entropy-based cell-type specificity index"""
        cell_type_means = gene_data.groupby('cell_type')['perturbation_score'].mean()
        
        # Normalize to probabilities
        probs = cell_type_means / cell_type_means.sum()
        probs = probs[probs > 0]  # Remove zero probabilities
        
        # Calculate entropy
        entropy = -np.sum(probs * np.log2(probs))
        
        # Normalize by max entropy (log2 of number of cell types)
        max_entropy = np.log2(len(cell_type_means))
        specificity = 1 - (entropy / max_entropy) if max_entropy > 0 else 0
        
        return specificity
    
    def statistical_analysis(self):
        """Perform comprehensive statistical analysis"""
        print("Performing statistical analysis...")
        
        # Initialize results dataframe
        genes = self.perturbation_results['gene'].unique()
        cell_types = self.perturbation_results['cell_type'].unique()
        
        results = []
        
        # Get population scores for effect size calculation
        population_scores = self.perturbation_results['perturbation_score'].values
        
        # Analysis by gene
        gene_stats = []
        for gene in genes:
            gene_data = self.perturbation_results[
                self.perturbation_results['gene'] == gene
            ]
            gene_scores = gene_data['perturbation_score'].values
            
            # Permutation test
            p_value, null_dist = self.permutation_test(gene_scores)
            
            # Effect size
            effect_size = self.cohens_d(gene_scores, population_scores)
            
            # Cell-type specificity
            specificity = self.cell_type_specificity_index(gene_data)
            
            # Confidence intervals
            mean_score = np.mean(gene_scores)
            se_score = stats.sem(gene_scores)
            ci_95 = stats.t.interval(
                0.95, len(gene_scores) - 1, loc=mean_score, scale=se_score
            )
            
            gene_stats.append({
                'gene': gene,
                'mean_score': mean_score,
                'se_score': se_score,
                'ci_95_low': ci_95[0],
                'ci_95_high': ci_95[1],
                'p_value': p_value,
                'effect_size': effect_size,
                'specificity_index': specificity,
                'n_cells': len(gene_scores)
            })
            
            # Analysis by gene x cell type
            for cell_type in cell_types:
                cell_gene_data = gene_data[gene_data['cell_type'] == cell_type]
                
                if len(cell_gene_data) == 0:
                    continue
                    
                scores = cell_gene_data['perturbation_score'].values
                p_val, _ = self.permutation_test(scores)
                effect = self.cohens_d(scores, population_scores)
                
                results.append({
                    'gene': gene,
                    'cell_type': cell_type,
                    'mean_score': np.mean(scores),
                    'se_score': stats.sem(scores) if len(scores) > 1 else 0,
                    'p_value': p_val,
                    'effect_size': effect,
                    'n_cells': len(scores)
                })
        
        # Convert to DataFrames
        self.gene_statistics = pd.DataFrame(gene_stats)
        self.gene_celltype_stats = pd.DataFrame(results)
        
        # FDR correction
        self.gene_statistics['p_value_fdr'] = false_discovery_control(
            self.gene_statistics['p_value']
        )
        self.gene_celltype_stats['p_value_fdr'] = false_discovery_control(
            self.gene_celltype_stats['p_value']
        )
        
        # Add significance stars
        def add_significance_stars(p_values):
            stars = []
            for p in p_values:
                if p < 0.001:
                    stars.append('***')
                elif p < 0.01:
                    stars.append('**')
                elif p < 0.05:
                    stars.append('*')
                else:
                    stars.append('')
            return stars
        
        self.gene_statistics['significance'] = add_significance_stars(
            self.gene_statistics['p_value_fdr']
        )
        
        print(f"Statistical analysis complete. {len(self.gene_statistics)} genes analyzed.")
        
        # Correlation analysis: expression vs perturbation
        self.correlation_analysis()
        
    def correlation_analysis(self):
        """Analyze correlation between expression level and perturbation effect"""
        print("Performing correlation analysis...")
        
        # Calculate correlation for each gene
        correlation_results = []
        
        for gene in self.perturbation_results['gene'].unique():
            gene_data = self.perturbation_results[
                self.perturbation_results['gene'] == gene
            ]
            
            # Correlation between original expression and perturbation score
            corr_coef, corr_p = stats.pearsonr(
                gene_data['original_expression'],
                gene_data['perturbation_score']
            )
            
            correlation_results.append({
                'gene': gene,
                'correlation': corr_coef,
                'correlation_p': corr_p,
                'r_squared': corr_coef ** 2
            })
        
        self.correlation_stats = pd.DataFrame(correlation_results)
        
        # Overall correlation
        overall_corr, overall_p = stats.pearsonr(
            self.perturbation_results['original_expression'],
            self.perturbation_results['perturbation_score']
        )
        
        print(f"Overall correlation: r = {overall_corr:.4f}, p = {overall_p:.2e}")
        self.overall_correlation = {
            'correlation': overall_corr,
            'p_value': overall_p,
            'r_squared': overall_corr ** 2
        }
        
    def create_heatmap_figure(self):
        """Create Main Figure 1: Heatmap with dendrograms"""
        print("Creating Main Figure 1: Heatmap...")
        
        # Create pivot table for heatmap
        heatmap_data = self.gene_celltype_stats.pivot(
            index='gene', columns='cell_type', values='mean_score'
        ).fillna(0)
        
        # Create figure with subplots for dendrograms
        fig = plt.figure(figsize=(12, 10))
        
        # Define grid
        gs = fig.add_gridspec(
            3, 3, 
            width_ratios=[0.2, 0.8, 0.1], 
            height_ratios=[0.2, 0.8, 0.1],
            hspace=0.05, wspace=0.05
        )
        
        # Row dendrogram
        ax_row_dend = fig.add_subplot(gs[1, 0])
        row_linkage = linkage(pdist(heatmap_data.values), method='ward')
        row_dend = dendrogram(
            row_linkage, orientation='left', ax=ax_row_dend,
            color_threshold=0.3*max(row_linkage[:, 2])
        )
        ax_row_dend.set_xticks([])
        ax_row_dend.set_yticks([])
        ax_row_dend.spines['top'].set_visible(False)
        ax_row_dend.spines['right'].set_visible(False)
        ax_row_dend.spines['bottom'].set_visible(False)
        ax_row_dend.spines['left'].set_visible(False)
        
        # Column dendrogram
        ax_col_dend = fig.add_subplot(gs[0, 1])
        col_linkage = linkage(pdist(heatmap_data.T.values), method='ward')
        col_dend = dendrogram(
            col_linkage, orientation='top', ax=ax_col_dend,
            color_threshold=0.3*max(col_linkage[:, 2])
        )
        ax_col_dend.set_xticks([])
        ax_col_dend.set_yticks([])
        ax_col_dend.spines['top'].set_visible(False)
        ax_col_dend.spines['right'].set_visible(False)
        ax_col_dend.spines['bottom'].set_visible(False)
        ax_col_dend.spines['left'].set_visible(False)
        
        # Reorder data according to dendrograms
        row_order = row_dend['leaves']
        col_order = col_dend['leaves']
        
        ordered_data = heatmap_data.iloc[row_order, col_order]
        
        # Main heatmap
        ax_heatmap = fig.add_subplot(gs[1, 1])
        im = ax_heatmap.imshow(
            ordered_data.values, cmap='viridis', aspect='auto',
            interpolation='nearest'
        )
        
        # Set labels
        ax_heatmap.set_xticks(range(len(ordered_data.columns)))
        ax_heatmap.set_yticks(range(len(ordered_data.index)))
        ax_heatmap.set_xticklabels(
            ordered_data.columns, rotation=45, ha='right'
        )
        ax_heatmap.set_yticklabels(ordered_data.index)
        
        ax_heatmap.set_xlabel('Cell Type')
        ax_heatmap.set_ylabel('Gene')
        
        # Colorbar
        ax_colorbar = fig.add_subplot(gs[1, 2])
        cbar = plt.colorbar(im, cax=ax_colorbar)
        cbar.set_label('Mean Perturbation Score')
        
        # Add panel letter
        fig.text(0.02, 0.98, 'A', fontsize=16, fontweight='bold', 
                transform=fig.transFigure)
        
        plt.suptitle('Gene × Cell Type Perturbation Effects', y=0.95)
        
        # Save in multiple formats
        plt.savefig(f"{self.output_dir}/Figure1_Heatmap.png", dpi=300)
        plt.savefig(f"{self.output_dir}/Figure1_Heatmap.pdf")
        plt.close()
        
    def create_ranked_bar_chart(self):
        """Create Main Figure 2: Ranked bar chart with error bars"""
        print("Creating Main Figure 2: Ranked bar chart...")
        
        # Sort genes by mean score
        sorted_genes = self.gene_statistics.sort_values(
            'mean_score', ascending=False
        )
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create bar chart
        x_pos = range(len(sorted_genes))
        bars = ax.bar(
            x_pos, 
            sorted_genes['mean_score'],
            yerr=sorted_genes['se_score'] * 1.96,  # 95% CI
            capsize=3,
            color=colors[0],
            alpha=0.8,
            edgecolor='black',
            linewidth=0.5
        )
        
        # Add significance stars
        for i, (_, row) in enumerate(sorted_genes.iterrows()):
            if row['significance']:
                ax.text(
                    i, row['mean_score'] + row['se_score'] * 1.96 + 
                    max(sorted_genes['mean_score']) * 0.02,
                    row['significance'],
                    ha='center', va='bottom', fontweight='bold'
                )
        
        # Formatting
        ax.set_xlabel('Gene')
        ax.set_ylabel('Mean Perturbation Score ± 95% CI')
        ax.set_title('B')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(
            sorted_genes['gene'], rotation=45, ha='right'
        )
        
        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        
        # Save in multiple formats
        plt.savefig(f"{self.output_dir}/Figure2_RankedBars.png", dpi=300)
        plt.savefig(f"{self.output_dir}/Figure2_RankedBars.pdf")
        plt.close()
        
    def create_volcano_plot(self):
        """Create Main Figure 3: Volcano plot"""
        print("Creating Main Figure 3: Volcano plot...")
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Calculate -log10(p-value)
        neg_log_p = -np.log10(self.gene_celltype_stats['p_value_fdr'] + 1e-300)
        
        # Create scatter plot
        scatter = ax.scatter(
            self.gene_celltype_stats['effect_size'],
            neg_log_p,
            c=self.gene_celltype_stats['mean_score'],
            cmap='viridis',
            alpha=0.6,
            s=30
        )
        
        # Add significance thresholds
        ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.7)
        ax.axhline(y=-np.log10(0.01), color='red', linestyle='--', alpha=0.7)
        ax.axvline(x=0, color='gray', linestyle='-', alpha=0.5)
        
        # Label significant points
        significant = self.gene_celltype_stats[
            (self.gene_celltype_stats['p_value_fdr'] < 0.05) & 
            (np.abs(self.gene_celltype_stats['effect_size']) > 0.5)
        ]
        
        for _, row in significant.iterrows():
            ax.annotate(
                f"{row['gene']}\n{row['cell_type'][:10]}",
                (row['effect_size'], -np.log10(row['p_value_fdr'] + 1e-300)),
                xytext=(5, 5), textcoords='offset points',
                fontsize=8, alpha=0.8
            )
        
        # Formatting
        ax.set_xlabel('Effect Size (Cohen\'s d)')
        ax.set_ylabel('-log₁₀(FDR-corrected p-value)')
        ax.set_title('C')
        
        # Colorbar
        cbar = plt.colorbar(scatter)
        cbar.set_label('Mean Perturbation Score')
        
        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        
        # Save in multiple formats
        plt.savefig(f"{self.output_dir}/Figure3_Volcano.png", dpi=300)
        plt.savefig(f"{self.output_dir}/Figure3_Volcano.pdf")
        plt.close()
        
    def create_drug_network(self):
        """Create Main Figure 4: Drug repurposing network"""
        print("Creating Main Figure 4: Drug network...")
        
        # Get top genes and their drug interactions
        top_genes = self.gene_ranking.head(10)['gene'].tolist()
        drug_subset = self.drug_candidates[
            self.drug_candidates['gene'].isin(top_genes) & 
            (self.drug_candidates['approved'] == True)
        ].head(30)  # Top 30 approved drugs
        
        # Count interactions by type
        interaction_counts = drug_subset.groupby(['gene', 'interaction_type']).size().reset_index(name='count')
        
        # Create Sankey-like visualization using stacked bars
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))
        
        # Left panel: Genes and their drug counts
        gene_drug_counts = drug_subset.groupby('gene').size().sort_values(ascending=True)
        
        bars1 = ax1.barh(range(len(gene_drug_counts)), gene_drug_counts.values, 
                        color=colors[1], alpha=0.8)
        ax1.set_yticks(range(len(gene_drug_counts)))
        ax1.set_yticklabels(gene_drug_counts.index)
        ax1.set_xlabel('Number of Drug Interactions')
        ax1.set_ylabel('Gene')
        ax1.set_title('Gene → Drug Interactions')
        
        # Right panel: Drug interaction types
        interaction_type_counts = drug_subset['interaction_type'].value_counts().head(8)
        
        bars2 = ax2.bar(range(len(interaction_type_counts)), interaction_type_counts.values,
                       color=colors[2], alpha=0.8)
        ax2.set_xticks(range(len(interaction_type_counts)))
        ax2.set_xticklabels(interaction_type_counts.index, rotation=45, ha='right')
        ax2.set_ylabel('Number of Interactions')
        ax2.set_title('Interaction Types')
        
        # Remove spines
        for ax in [ax1, ax2]:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        
        # Add panel letter
        fig.text(0.02, 0.95, 'D', fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        # Save in multiple formats
        plt.savefig(f"{self.output_dir}/Figure4_DrugNetwork.png", dpi=300)
        plt.savefig(f"{self.output_dir}/Figure4_DrugNetwork.pdf")
        plt.close()
        
    def create_expression_scatter(self):
        """Create Supplementary Figure 1: Expression vs Perturbation scatter"""
        print("Creating Supplementary Figure 1: Expression scatter...")
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Sample data for plotting (too many points otherwise)
        sample_size = min(5000, len(self.perturbation_results))
        sample_data = self.perturbation_results.sample(n=sample_size, random_state=42)
        
        # Create scatter plot
        scatter = ax.scatter(
            sample_data['original_expression'],
            sample_data['perturbation_score'],
            alpha=0.5,
            s=10,
            c=colors[3]
        )
        
        # Add regression line
        x = sample_data['original_expression']
        y = sample_data['perturbation_score']
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        ax.plot(x, p(x), color='red', linewidth=2, alpha=0.8)
        
        # Add statistics text
        corr = self.overall_correlation['correlation']
        p_val = self.overall_correlation['p_value']
        r_squared = self.overall_correlation['r_squared']
        
        ax.text(0.05, 0.95, 
                f'R² = {r_squared:.4f}\nr = {corr:.4f}\np = {p_val:.2e}',
                transform=ax.transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Formatting
        ax.set_xlabel('Original Expression Level')
        ax.set_ylabel('Perturbation Score')
        ax.set_title('Expression Level vs Perturbation Effect')
        
        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        
        # Save in multiple formats
        plt.savefig(f"{self.output_dir}/SupplementaryFigure1_ExpressionScatter.png", dpi=300)
        plt.savefig(f"{self.output_dir}/SupplementaryFigure1_ExpressionScatter.pdf")
        plt.close()
        
    def create_cell_type_profiles(self):
        """Create Supplementary Figure 2: Cell-type profiles for top 5 genes"""
        print("Creating Supplementary Figure 2: Cell-type profiles...")
        
        # Get top 5 genes
        top_5_genes = self.gene_ranking.head(5)['gene'].tolist()
        
        # Create radar plot data
        fig, axes = plt.subplots(1, 5, figsize=(20, 4), subplot_kw=dict(projection='polar'))
        if not isinstance(axes, np.ndarray):
            axes = [axes]
        
        # Get cell types
        cell_types = sorted(self.perturbation_results['cell_type'].unique())
        angles = np.linspace(0, 2 * np.pi, len(cell_types), endpoint=False)
        
        for i, gene in enumerate(top_5_genes):
            ax = axes[i]
            
            # Get data for this gene
            gene_data = self.cell_type_summary[
                self.cell_type_summary['gene_'] == gene
            ].set_index('cell_type_')
            
            # Ensure all cell types are represented
            values = []
            for ct in cell_types:
                if ct in gene_data.index:
                    values.append(gene_data.loc[ct, 'perturbation_score_mean'])
                else:
                    values.append(0)
            
            # Close the plot
            values += values[:1]
            angles_plot = list(angles) + [angles[0]]
            
            # Plot
            ax.plot(angles_plot, values, 'o-', linewidth=2, color=colors[i])
            ax.fill(angles_plot, values, alpha=0.25, color=colors[i])
            
            # Labels
            ax.set_xticks(angles)
            ax.set_xticklabels([ct.replace(' ', '\n') for ct in cell_types], fontsize=8)
            ax.set_title(gene, pad=20)
            ax.grid(True)
            
        plt.tight_layout()
        
        # Save in multiple formats
        plt.savefig(f"{self.output_dir}/SupplementaryFigure2_CellTypeProfiles.png", dpi=300)
        plt.savefig(f"{self.output_dir}/SupplementaryFigure2_CellTypeProfiles.pdf")
        plt.close()
        
    def generate_results_table(self):
        """Generate LaTeX results table"""
        print("Generating LaTeX results table...")
        
        # Merge gene statistics with drug information
        drug_counts = self.drug_candidates.groupby('gene')['drug_name'].count().reset_index()
        drug_counts.rename(columns={'drug_name': 'drug_count'}, inplace=True)
        
        top_drugs = self.drug_candidates.loc[
            self.drug_candidates.groupby('gene')['perturbation_score'].idxmax()
        ][['gene', 'drug_name']].rename(columns={'drug_name': 'top_drug'})
        
        # Find top cell type for each gene
        top_cell_types = self.gene_celltype_stats.loc[
            self.gene_celltype_stats.groupby('gene')['mean_score'].idxmax()
        ][['gene', 'cell_type']].rename(columns={'cell_type': 'top_cell_type'})
        
        # Merge all data
        table_data = self.gene_statistics.merge(
            drug_counts.reset_index(), on='gene', how='left'
        ).merge(
            top_drugs, on='gene', how='left'
        ).merge(
            top_cell_types, on='gene', how='left'
        )
        
        # Sort by mean score (descending)
        table_data = table_data.sort_values('mean_score', ascending=False)
        
        # Generate LaTeX table
        latex_table = """
\\begin{table}[htbp]
\\centering
\\caption{Statistical Results for Virtual Perturbation Screening}
\\label{tab:perturbation_results}
\\begin{tabular}{lcccccc}
\\toprule
Gene & Mean Effect ± SE & FDR p-value & Top Cell Type & Drug Interactions & Top Drug \\\\
\\midrule
"""
        
        for _, row in table_data.head(15).iterrows():  # Top 15 genes
            mean_se = f"{row['mean_score']:.4f} ± {row['se_score']:.4f}"
            p_val = f"{row['p_value_fdr']:.2e}" if row['p_value_fdr'] < 0.01 else f"{row['p_value_fdr']:.3f}"
            drug_count = int(row['drug_count']) if pd.notna(row['drug_count']) else 0
            top_drug = str(row['top_drug'])[:20] if pd.notna(row['top_drug']) else 'N/A'
            top_cell = str(row['top_cell_type'])[:15] if pd.notna(row['top_cell_type']) else 'N/A'
            
            latex_table += f"{row['gene']} & {mean_se} & {p_val} & {top_cell} & {drug_count} & {top_drug} \\\\\n"
        
        latex_table += """
\\bottomrule
\\end{tabular}
\\end{table}
"""
        
        # Save table
        with open(f"{self.output_dir}/results_table.tex", 'w') as f:
            f.write(latex_table)
            
    def save_statistical_results(self):
        """Save all statistical results to CSV files"""
        print("Saving statistical results...")
        
        # Gene-level statistics
        self.gene_statistics.to_csv(
            f"{self.output_dir}/gene_statistics.csv", index=False
        )
        
        # Gene x cell type statistics  
        self.gene_celltype_stats.to_csv(
            f"{self.output_dir}/gene_celltype_statistics.csv", index=False
        )
        
        # Correlation results
        self.correlation_stats.to_csv(
            f"{self.output_dir}/correlation_statistics.csv", index=False
        )
        
        print("Statistical results saved.")
        
    def run_full_analysis(self):
        """Run the complete Phase 3 analysis"""
        print("Starting Phase 3 Analysis...")
        
        # Statistical analysis
        self.statistical_analysis()
        
        # Generate all figures
        self.create_heatmap_figure()
        self.create_ranked_bar_chart()
        self.create_volcano_plot()
        self.create_drug_network()
        self.create_expression_scatter()
        self.create_cell_type_profiles()
        
        # Generate results table
        self.generate_results_table()
        
        # Save statistical results
        self.save_statistical_results()
        
        print("Phase 3 Analysis Complete!")
        print(f"Results saved to: {self.output_dir}")
        

if __name__ == "__main__":
    # Run analysis
    data_dir = r"D:\openclaw\intelligence-augmentation\analysis\results"
    analysis = Phase3Analysis(data_dir)
    analysis.run_full_analysis()