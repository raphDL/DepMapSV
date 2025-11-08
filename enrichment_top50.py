#!/usr/bin/env python3
"""Re-run enrichment with top 50 genes only."""
import pandas as pd
import numpy as np
from pathlib import Path

# Load candidates
candidates = pd.read_csv("out_v3/candidate_genes_fdr01.csv")

# Take top 50 by |prox_coef|
top50 = candidates.head(50)

print("=== Top 50 Genes for Enrichment ===")
print(f"Selected {len(top50)} genes (top by |proximity coefficient|)")
print(f"Median |prox_coef|: {top50['abs_prox_coef'].median():.4f}")
print(f"Range: {top50['abs_prox_coef'].min():.4f} to {top50['abs_prox_coef'].max():.4f}")

print("\n=== Gene List (for g:Profiler/Enrichr) ===")
gene_list = top50['gene'].tolist()
print("\n".join(gene_list))

# Save gene list
Path("out_v3/enrichment").mkdir(parents=True, exist_ok=True)
with open("out_v3/enrichment/gene_list_top50.txt", "w") as f:
    f.write("\n".join(gene_list))
print(f"\n✓ Saved to out_v3/enrichment/gene_list_top50.txt")

# Save full details
top50.to_csv("out_v3/candidate_genes_top50.csv", index=False)
print(f"✓ Saved full details to out_v3/candidate_genes_top50.csv")

# Summary statistics
print("\n=== Top 50 Summary ===")
print(f"Same sign (CN & prox): {(np.sign(top50['coef_cn']) == np.sign(top50['coef_prox_exp_100k'])).sum()}")
print(f"Mean R²: {top50['r2'].mean():.4f}")
print(f"Genes with R² > 0.05: {(top50['r2'] > 0.05).sum()}")

# Try Enrichr if available
print("\n=== Running Enrichr (if available) ===")
try:
    from gseapy import enrichr
    
    # Run enrichment
    enr = enrichr(
        gene_list=gene_list,
        gene_sets=['GO_Biological_Process_2021', 'KEGG_2021_Human', 'WikiPathway_2021_Human'],
        organism='human',
        outdir='out_v3/enrichment_top50',
        cutoff=0.05
    )
    
    print("\n✓ Enrichment results saved to out_v3/enrichment_top50/")
    
    # Print top results
    try:
        results = pd.read_csv("out_v3/enrichment_top50/GO_Biological_Process_2021.txt", sep="\t")
        print("\nTop 10 enriched pathways:")
        print(results[['Term', 'Adjusted P-value', 'Genes']].head(10))
    except FileNotFoundError:
        print("  (Results file not found - check out_v3/enrichment_top50/)")
    
except ImportError:
    print("⚠️  gseapy not installed. Install with: pip install gseapy")
    print("Or use web version: https://maayanlab.cloud/Enrichr/")
    print("\n→ Copy the gene list above and paste into Enrichr")

print("\n=== Expected Enrichment (if hypothesis correct) ===")
print("  - DNA repair (GO:0006281)")
print("  - Response to DNA damage (GO:0006974)")
print("  - Double-strand break repair (GO:0006302)")
print("  - Chromatin remodeling (GO:0006338)")
print("  - Cell cycle checkpoint (GO:0000075)")

