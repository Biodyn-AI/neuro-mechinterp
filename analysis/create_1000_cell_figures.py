#!/usr/bin/env python3
"""Generate comparison figures for 1000-cell vs 500-cell analysis"""
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style
plt.style.use('default')
sns.set_palette("husl")

def load_results():
    """Load both 500-cell and 1000-cell results"""
    
    # 500-cell baseline (from FINAL_REPORT.md)
    baseline_500 = {
        "CADM2": {"shift_mean": 0.0196, "shift_std": 0.0045, "n_cells": 472},
        "GRIN2A": {"shift_mean": 0.0190, "shift_std": 0.0061, "n_cells": 223},
        "CAMK2A": {"shift_mean": 0.0189, "shift_std": 0.0059, "n_cells": 155},
        "MEF2C": {"shift_mean": 0.0184, "shift_std": 0.0065, "n_cells": 216},
        "APP": {"shift_mean": 0.0183, "shift_std": 0.0054, "n_cells": 430}
    }
    
    # 1000-cell results
    results_1000_file = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/results_1000.json"
    try:
        with open(results_1000_file, 'r') as f:
            results_1000 = json.load(f)
    except:
        print("❌ Could not load 1000-cell results")
        return None, None
    
    return baseline_500, results_1000

def create_comparison_barplot(baseline_500, results_1000):
    """Create side-by-side bar plot comparing effect sizes"""
    genes = ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]
    
    shifts_500 = [baseline_500[gene]["shift_mean"] for gene in genes]
    shifts_1000 = [results_1000[gene]["shift_mean"] for gene in genes]
    errors_500 = [baseline_500[gene]["shift_std"] for gene in genes]
    errors_1000 = [results_1000[gene]["shift_std"] for gene in genes]
    
    x = np.arange(len(genes))
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    bars1 = ax.bar(x - width/2, shifts_500, width, yerr=errors_500, 
                   label='500 cells', alpha=0.8, capsize=5)
    bars2 = ax.bar(x + width/2, shifts_1000, width, yerr=errors_1000,
                   label='1000 cells', alpha=0.8, capsize=5)
    
    ax.set_xlabel('Intelligence Genes', fontsize=12)
    ax.set_ylabel('Perturbation Shift (1 - cosine similarity)', fontsize=12)
    ax.set_title('Gene Perturbation Effects: 500 vs 1000 Cells', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(genes)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add pathway annotations
    pathways = {
        'CADM2': 'Cell Adhesion',
        'GRIN2A': 'Glutamate',
        'CAMK2A': 'Neurotrophic',
        'MEF2C': 'Neurodevelopmental',
        'APP': 'Neurodegeneration'
    }
    
    for i, gene in enumerate(genes):
        ax.text(i, max(shifts_500[i], shifts_1000[i]) + max(errors_500[i], errors_1000[i]) + 0.001,
                pathways[gene], ha='center', va='bottom', fontsize=9, rotation=45)
    
    plt.tight_layout()
    return fig

def create_ranking_comparison(baseline_500, results_1000):
    """Create ranking comparison visualization"""
    genes = ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]
    
    # Get rankings
    ranking_500 = sorted(genes, key=lambda x: baseline_500[x]["shift_mean"], reverse=True)
    ranking_1000 = sorted(genes, key=lambda x: results_1000[x]["shift_mean"], reverse=True)
    
    # Create rank positions
    ranks_500 = {gene: i+1 for i, gene in enumerate(ranking_500)}
    ranks_1000 = {gene: i+1 for i, gene in enumerate(ranking_1000)}
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot lines connecting rankings
    for gene in genes:
        x_vals = [0, 1]
        y_vals = [ranks_500[gene], ranks_1000[gene]]
        
        # Color by rank change
        if ranks_1000[gene] < ranks_500[gene]:
            color = 'green'  # Moved up
            alpha = 0.8
        elif ranks_1000[gene] > ranks_500[gene]:
            color = 'red'    # Moved down
            alpha = 0.8
        else:
            color = 'gray'   # No change
            alpha = 0.6
        
        ax.plot(x_vals, y_vals, 'o-', color=color, alpha=alpha, linewidth=2, markersize=8)
        
        # Add gene labels
        ax.text(-0.1, ranks_500[gene], gene, ha='right', va='center', fontweight='bold')
        ax.text(1.1, ranks_1000[gene], gene, ha='left', va='center', fontweight='bold')
    
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(0.5, 5.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['500 cells', '1000 cells'], fontsize=12)
    ax.set_ylabel('Rank (1 = highest effect)', fontsize=12)
    ax.set_title('Gene Ranking Changes: 500 vs 1000 Cells', fontsize=14, fontweight='bold')
    ax.invert_yaxis()  # Rank 1 at top
    ax.grid(True, alpha=0.3)
    
    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='green', lw=2, label='Rank improved'),
        Line2D([0], [0], color='red', lw=2, label='Rank declined'), 
        Line2D([0], [0], color='gray', lw=2, label='No change')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    return fig

def create_effect_ratio_plot(baseline_500, results_1000):
    """Create plot showing effect size ratios"""
    genes = ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]
    
    ratios = []
    for gene in genes:
        ratio = results_1000[gene]["shift_mean"] / baseline_500[gene]["shift_mean"]
        ratios.append(ratio)
    
    # Sort by ratio for better visualization
    sorted_data = sorted(zip(genes, ratios), key=lambda x: x[1])
    sorted_genes, sorted_ratios = zip(*sorted_data)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    colors = ['red' if r < 1.0 else 'green' for r in sorted_ratios]
    bars = ax.barh(range(len(sorted_genes)), sorted_ratios, color=colors, alpha=0.7)
    
    # Add 1.0 reference line
    ax.axvline(1.0, color='black', linestyle='--', alpha=0.8, linewidth=2)
    ax.text(1.02, len(sorted_genes)/2, 'Equal effect', rotation=90, va='center', fontweight='bold')
    
    ax.set_yticks(range(len(sorted_genes)))
    ax.set_yticklabels(sorted_genes)
    ax.set_xlabel('Effect Size Ratio (1000-cell / 500-cell)', fontsize=12)
    ax.set_title('Effect Size Scaling: 1000-cell vs 500-cell', fontsize=14, fontweight='bold')
    
    # Add value labels on bars
    for i, (gene, ratio) in enumerate(zip(sorted_genes, sorted_ratios)):
        ax.text(ratio + 0.005, i, f'{ratio:.3f}', va='center', fontsize=10)
    
    ax.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    return fig

def create_cell_count_comparison(baseline_500, results_1000):
    """Create plot comparing cell counts and statistical power"""
    genes = ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]
    
    counts_500 = [baseline_500[gene]["n_cells"] for gene in genes]
    counts_1000 = [results_1000[gene]["n_cells_processed"] for gene in genes]
    
    x = np.arange(len(genes))
    width = 0.35
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Cell count comparison
    bars1 = ax1.bar(x - width/2, counts_500, width, label='500-cell dataset', alpha=0.8)
    bars2 = ax1.bar(x + width/2, counts_1000, width, label='1000-cell dataset', alpha=0.8)
    
    ax1.set_xlabel('Genes', fontsize=12)
    ax1.set_ylabel('Number of Expressing Cells', fontsize=12)
    ax1.set_title('Cell Count Comparison by Gene', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(genes)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Standard error comparison (inversely related to statistical power)
    se_500 = [baseline_500[gene]["shift_std"] / np.sqrt(baseline_500[gene]["n_cells"]) for gene in genes]
    se_1000 = [results_1000[gene]["shift_std"] / np.sqrt(results_1000[gene]["n_cells_processed"]) for gene in genes]
    
    bars3 = ax2.bar(x - width/2, se_500, width, label='500-cell dataset', alpha=0.8)
    bars4 = ax2.bar(x + width/2, se_1000, width, label='1000-cell dataset', alpha=0.8)
    
    ax2.set_xlabel('Genes', fontsize=12)
    ax2.set_ylabel('Standard Error of Mean Shift', fontsize=12)
    ax2.set_title('Statistical Precision Comparison', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(genes)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig

def main():
    """Generate all comparison figures"""
    print("🎨 Creating 1000-cell vs 500-cell comparison figures")
    print("=" * 60)
    
    # Load results
    baseline_500, results_1000 = load_results()
    if baseline_500 is None or results_1000 is None:
        return
    
    # Create output directory
    output_dir = Path("/mnt/d/openclaw/intelligence-augmentation/analysis/figures")
    output_dir.mkdir(exist_ok=True)
    
    # Create figures
    figures = []
    
    print("📊 Creating effect size comparison...")
    fig1 = create_comparison_barplot(baseline_500, results_1000)
    fig1.savefig(output_dir / "1000_vs_500_effect_sizes.png", dpi=300, bbox_inches='tight')
    figures.append("1000_vs_500_effect_sizes.png")
    
    print("📈 Creating ranking comparison...")
    fig2 = create_ranking_comparison(baseline_500, results_1000)
    fig2.savefig(output_dir / "1000_vs_500_rankings.png", dpi=300, bbox_inches='tight')
    figures.append("1000_vs_500_rankings.png")
    
    print("🔄 Creating effect ratio plot...")
    fig3 = create_effect_ratio_plot(baseline_500, results_1000)
    fig3.savefig(output_dir / "1000_vs_500_effect_ratios.png", dpi=300, bbox_inches='tight')
    figures.append("1000_vs_500_effect_ratios.png")
    
    print("📋 Creating statistical power comparison...")
    fig4 = create_cell_count_comparison(baseline_500, results_1000)
    fig4.savefig(output_dir / "1000_vs_500_statistical_power.png", dpi=300, bbox_inches='tight')
    figures.append("1000_vs_500_statistical_power.png")
    
    plt.close('all')  # Clean up
    
    print(f"\n✅ Generated {len(figures)} comparison figures:")
    for fig in figures:
        print(f"   📁 {output_dir / fig}")
    
    print(f"\n🎉 Figure generation complete!")

if __name__ == "__main__":
    main()