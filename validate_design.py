#!/usr/bin/env python3
"""Validate design matrix output."""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    design_path = "test_data/design.csv"
else:
    design_path = sys.argv[1]

# Load design matrix
try:
    design = pd.read_csv(design_path)
except FileNotFoundError:
    print(f"Error: Could not find {design_path}")
    sys.exit(1)

print("=== Design Matrix Summary ===")
print(f"Shape: {design.shape}")
print(f"Genes: {design['gene'].nunique()}")
print(f"Cell lines: {design['cell_line'].nunique()}")
print(f"Chromosomes: {design['chrom'].nunique()}")

print("\n=== CN Statistics ===")
print(design['cn'].describe())

print("\n=== Breakpoint Distance Statistics ===")
print(design['bp_dist'].describe())

print("\n=== Proximity Kernel Statistics ===")
for col in design.columns:
    if col.startswith('prox_'):
        non_zero_pct = (design[col] > 1e-10).mean() * 100
        print(f"{col}: mean={design[col].mean():.6f}, non-zero={non_zero_pct:.1f}%")

print("\n=== CN-Proximity Correlations ===")
prox_cols = [c for c in design.columns if c.startswith('prox_')]
if 'cn' in design.columns and prox_cols:
    corr = design[['cn'] + prox_cols].corr().loc['cn']
    for col in prox_cols:
        print(f"{col}: r={corr[col]:.3f}")
    
    max_corr = corr[prox_cols].abs().max()
    if max_corr > 0.4:
        print(f"\n⚠️  WARNING: High CN-proximity correlation ({max_corr:.3f})")
        print("This suggests entanglement - SV breakpoints may still be confounded with CN")
    else:
        print(f"\n✓ CN-proximity correlations are acceptable (max={max_corr:.3f})")

# Plot distributions
try:
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    axes[0, 0].hist(design['cn'], bins=50, edgecolor='black')
    axes[0, 0].set_title('Copy Number Distribution')
    axes[0, 0].set_xlabel('CN')
    
    axes[0, 1].hist(np.log10(design['bp_dist'] + 1), bins=50, edgecolor='black')
    axes[0, 1].set_title('Breakpoint Distance (log10)')
    axes[0, 1].set_xlabel('log10(bp_dist + 1)')
    
    if 'prox_exp_100k' in design.columns:
        axes[0, 2].hist(design['prox_exp_100k'], bins=50, edgecolor='black')
        axes[0, 2].set_title('Proximity Kernel (exp, λ=100k)')
        axes[0, 2].set_xlabel('prox_exp_100k')
    
    if 'cn' in design.columns and 'prox_exp_100k' in design.columns:
        axes[1, 0].scatter(design['cn'], design['prox_exp_100k'], alpha=0.1, s=1)
        if 'cn' in design.columns and prox_cols:
            corr_val = corr.get('prox_exp_100k', 0)
            axes[1, 0].set_title(f'CN vs Proximity (r={corr_val:.3f})')
        axes[1, 0].set_xlabel('CN')
        axes[1, 0].set_ylabel('prox_exp_100k')
    
    if 'prox_exp_100k' in design.columns:
        axes[1, 1].scatter(np.log10(design['bp_dist'] + 1), design['prox_exp_100k'], alpha=0.1, s=1)
        axes[1, 1].set_title('Distance vs Proximity (should be monotonic)')
        axes[1, 1].set_xlabel('log10(bp_dist + 1)')
        axes[1, 1].set_ylabel('prox_exp_100k')
    
    # Genes per cell line
    genes_per_cell = design.groupby('cell_line')['gene'].nunique()
    axes[1, 2].hist(genes_per_cell, bins=30, edgecolor='black')
    axes[1, 2].set_title('Genes per Cell Line')
    axes[1, 2].set_xlabel('Number of genes')
    
    plt.tight_layout()
    plt.savefig('logs/design_matrix_validation.png', dpi=150)
    print("\n✓ Saved validation plots to logs/design_matrix_validation.png")
except Exception as e:
    print(f"\n⚠️  Could not generate plots: {e}")

# Check specific genes if this is synthetic test data
if 'GENE_A' in design['gene'].values:
    print("\n=== Synthetic Test Validation ===")
    gene_a = design[design['gene'] == 'GENE_A']
    print(f"GENE_A bp_dist range: {gene_a['bp_dist'].min()} - {gene_a['bp_dist'].max()}")
    print("  Expected: ~100,000 bp (SV at 1.1Mb, gene at 1.0Mb)")
    
    gene_h = design[design['gene'] == 'GENE_H']
    print(f"GENE_H bp_dist range: {gene_h['bp_dist'].min()} - {gene_h['bp_dist'].max()}")
    print("  Expected: >> 1Mb (no nearby SVs on chr2)")

