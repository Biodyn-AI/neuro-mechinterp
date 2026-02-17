#!/usr/bin/env python3
"""
Cross-Model Validation: Correlation-Based Analysis
Alternative to scGPT perturbation for cross-validating Geneformer results
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
from scipy.sparse import load_npz
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json

class CorrelationCentralityAnalyzer:
    def __init__(self, results_dir="results"):
        self.results_dir = Path(results_dir)
        self.intelligence_genes = [
            'CADM2', 'GRIN2A', 'CAMK2A', 'MEF2C', 'APP', 'SCN1A', 'NRXN1', 
            'GRIN2B', 'HOMER1', 'NEGR1', 'SLC6A4', 'NLGN1', 'TCF4', 'FOXO3',
            'BDNF', 'MAPT', 'CREB1', 'FMR1', 'SYN1', 'COMT', 'SHANK3'
        ]
        
    def load_geneformer_results(self):
        """Load Geneformer perturbation results"""
        results_file = self.results_dir / "perturbation_results_streamlined.csv"
        if results_file.exists():
            df = pd.read_csv(results_file)
            return df
        else:
            print(f"❌ Geneformer results not found at {results_file}")
            return None
            
    def generate_mock_expression_data(self, n_cells=500, n_genes=1000):
        """
        Generate mock single-cell expression data for demonstration
        In real analysis, this would load actual brain scRNA-seq data
        """
        print("📊 Generating mock brain cell expression data...")
        
        # Create correlated expression patterns that might exist in real brain cells
        np.random.seed(42)
        
        # Generate base expression levels
        base_expr = np.random.lognormal(mean=1, sigma=1.5, size=(n_cells, n_genes))
        
        # Add correlation structure for intelligence genes
        gene_names = [f"GENE_{i}" for i in range(n_genes)]
        
        # Replace some genes with our intelligence genes
        for i, gene in enumerate(self.intelligence_genes):
            if i < len(gene_names):
                gene_names[i] = gene
                
        # Add some biological correlation structure
        # Synaptic genes tend to be co-expressed
        synaptic_genes = ['SHANK3', 'HOMER1', 'NRXN1', 'NLGN1', 'SYN1', 'CAMK2A']
        synaptic_indices = [i for i, gene in enumerate(gene_names) if gene in synaptic_genes]
        
        # Make synaptic genes correlated
        if len(synaptic_indices) > 1:
            synaptic_factor = np.random.normal(0, 1, n_cells).reshape(-1, 1)
            for idx in synaptic_indices:
                base_expr[:, idx] += 0.3 * synaptic_factor.flatten() * base_expr[:, idx]
                
        # NMDA receptor subunits should be correlated
        nmda_genes = ['GRIN2A', 'GRIN2B']
        nmda_indices = [i for i, gene in enumerate(gene_names) if gene in nmda_genes]
        if len(nmda_indices) > 1:
            nmda_factor = np.random.normal(0, 1, n_cells).reshape(-1, 1)
            for idx in nmda_indices:
                base_expr[:, idx] += 0.4 * nmda_factor.flatten() * base_expr[:, idx]
                
        return pd.DataFrame(base_expr, columns=gene_names)
        
    def compute_gene_correlations(self, expression_df):
        """Compute gene-gene correlation matrix"""
        print("🔗 Computing gene-gene correlations...")
        
        # Focus on intelligence genes that are present
        present_genes = [g for g in self.intelligence_genes if g in expression_df.columns]
        expr_subset = expression_df[present_genes]
        
        # Compute correlation matrix
        corr_matrix = expr_subset.corr(method='spearman')
        
        return corr_matrix, present_genes
        
    def compute_centrality_metrics(self, corr_matrix):
        """Compute various centrality metrics from correlation matrix"""
        print("📐 Computing correlation centrality metrics...")
        
        # Convert correlations to adjacency matrix (absolute values)
        adj_matrix = np.abs(corr_matrix.values)
        np.fill_diagonal(adj_matrix, 0)  # Remove self-loops
        
        centrality_metrics = {}
        gene_names = corr_matrix.columns.tolist()
        
        # Degree centrality (sum of absolute correlations)
        degree_centrality = np.sum(adj_matrix, axis=1)
        
        # Weighted degree centrality (sum of squared correlations)
        weighted_centrality = np.sum(adj_matrix**2, axis=1)
        
        # Average correlation centrality
        avg_correlation = np.mean(adj_matrix, axis=1)
        
        # Maximum correlation centrality
        max_correlation = np.max(adj_matrix, axis=1)
        
        for i, gene in enumerate(gene_names):
            centrality_metrics[gene] = {
                'degree_centrality': degree_centrality[i],
                'weighted_centrality': weighted_centrality[i],
                'avg_correlation': avg_correlation[i],
                'max_correlation': max_correlation[i]
            }
            
        return centrality_metrics
        
    def rank_genes_by_centrality(self, centrality_metrics):
        """Rank genes by different centrality metrics"""
        
        rankings = {}
        metrics = ['degree_centrality', 'weighted_centrality', 'avg_correlation', 'max_correlation']
        
        for metric in metrics:
            # Sort genes by metric value (descending)
            sorted_genes = sorted(centrality_metrics.keys(), 
                                key=lambda g: centrality_metrics[g][metric], 
                                reverse=True)
            
            # Create rank mapping (1 = highest centrality)
            rankings[metric] = {gene: rank+1 for rank, gene in enumerate(sorted_genes)}
            
        return rankings
        
    def compare_with_geneformer(self, geneformer_df, centrality_rankings):
        """Compare correlation centrality rankings with Geneformer perturbation effects"""
        print("⚖️  Comparing correlation centrality with Geneformer rankings...")
        
        # Create Geneformer ranking (by perturbation impact)
        geneformer_df = geneformer_df.sort_values('perturbation_impact', ascending=False)
        geneformer_ranking = {gene: rank+1 for rank, gene in enumerate(geneformer_df['gene_symbol'])}
        
        comparison_results = {}
        
        # Compare each centrality metric with Geneformer
        for metric, rankings in centrality_rankings.items():
            # Get common genes
            common_genes = set(geneformer_ranking.keys()) & set(rankings.keys())
            
            if len(common_genes) < 3:
                print(f"⚠️  Too few common genes for {metric}")
                continue
                
            # Extract rankings for common genes
            gf_ranks = [geneformer_ranking[gene] for gene in common_genes]
            cent_ranks = [rankings[gene] for gene in common_genes]
            
            # Compute correlations
            spearman_corr, spearman_p = spearmanr(gf_ranks, cent_ranks)
            pearson_corr, pearson_p = pearsonr(gf_ranks, cent_ranks)
            
            comparison_results[metric] = {
                'spearman_correlation': spearman_corr,
                'spearman_p_value': spearman_p,
                'pearson_correlation': pearson_corr,
                'pearson_p_value': pearson_p,
                'n_genes': len(common_genes),
                'common_genes': list(common_genes)
            }
            
        return comparison_results, geneformer_ranking
        
    def create_comparison_plots(self, geneformer_ranking, centrality_rankings, comparison_results):
        """Create visualization plots"""
        print("📊 Creating comparison plots...")
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Correlation Centrality vs. Geneformer Perturbation Rankings', fontsize=16)
        
        metrics = ['degree_centrality', 'weighted_centrality', 'avg_correlation', 'max_correlation']
        metric_labels = ['Degree Centrality', 'Weighted Centrality', 'Average Correlation', 'Max Correlation']
        
        for i, (metric, label) in enumerate(zip(metrics, metric_labels)):
            ax = axes[i//2, i%2]
            
            if metric not in comparison_results:
                ax.text(0.5, 0.5, f'No data for {label}', ha='center', va='center', transform=ax.transAxes)
                continue
                
            # Get data for scatter plot
            common_genes = comparison_results[metric]['common_genes']
            gf_ranks = [geneformer_ranking[gene] for gene in common_genes]
            cent_ranks = [centrality_rankings[metric][gene] for gene in common_genes]
            
            # Scatter plot
            ax.scatter(gf_ranks, cent_ranks, alpha=0.7)
            
            # Add gene labels
            for j, gene in enumerate(common_genes):
                ax.annotate(gene, (gf_ranks[j], cent_ranks[j]), 
                          xytext=(5, 5), textcoords='offset points', fontsize=8)
                
            # Add correlation info
            corr = comparison_results[metric]['spearman_correlation']
            p_val = comparison_results[metric]['spearman_p_value']
            ax.text(0.05, 0.95, f'ρ = {corr:.3f}\np = {p_val:.3f}', 
                   transform=ax.transAxes, verticalalignment='top', fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                   
            ax.set_xlabel('Geneformer Rank (1 = highest impact)')
            ax.set_ylabel(f'{label} Rank (1 = highest centrality)')
            ax.set_title(f'{label} vs. Geneformer')
            
        plt.tight_layout()
        plt.savefig('correlation_centrality_comparison.png', dpi=300, bbox_inches='tight')
        print("📸 Saved comparison plot as correlation_centrality_comparison.png")
        
        return fig
        
    def run_full_analysis(self):
        """Run the complete correlation-based cross-validation analysis"""
        print("🧬 Starting Correlation-Based Cross-Validation Analysis")
        print("="*60)
        
        # Load Geneformer results
        geneformer_df = self.load_geneformer_results()
        if geneformer_df is None:
            return None
            
        # Generate/load expression data
        expression_df = self.generate_mock_expression_data()
        print(f"📋 Expression data shape: {expression_df.shape}")
        
        # Compute correlations
        corr_matrix, present_genes = self.compute_gene_correlations(expression_df)
        print(f"🧮 Correlation matrix computed for {len(present_genes)} intelligence genes")
        
        # Compute centrality metrics
        centrality_metrics = self.compute_centrality_metrics(corr_matrix)
        
        # Rank genes by centrality
        centrality_rankings = self.rank_genes_by_centrality(centrality_metrics)
        
        # Compare with Geneformer
        comparison_results, geneformer_ranking = self.compare_with_geneformer(geneformer_df, centrality_rankings)
        
        # Create plots
        fig = self.create_comparison_plots(geneformer_ranking, centrality_rankings, comparison_results)
        
        # Compile final results
        results = {
            'geneformer_ranking': geneformer_ranking,
            'centrality_rankings': centrality_rankings,
            'comparison_results': comparison_results,
            'correlation_matrix': corr_matrix.to_dict(),
            'present_genes': present_genes
        }
        
        return results

if __name__ == "__main__":
    analyzer = CorrelationCentralityAnalyzer()
    results = analyzer.run_full_analysis()
    
    if results:
        # Save results
        with open('correlation_analysis_results.json', 'w') as f:
            # Convert numpy types to Python types for JSON serialization
            json_results = {}
            for key, value in results.items():
                if key == 'correlation_matrix':
                    json_results[key] = value
                elif key == 'comparison_results':
                    json_results[key] = {
                        k: {
                            sub_k: (float(sub_v) if isinstance(sub_v, (np.float64, np.float32)) else sub_v)
                            for sub_k, sub_v in v.items()
                        }
                        for k, v in value.items()
                    }
                else:
                    json_results[key] = value
            
            json.dump(json_results, f, indent=2)
            
        print("💾 Results saved to correlation_analysis_results.json")
        print("✅ Analysis complete!")
    else:
        print("❌ Analysis failed")