#!/usr/bin/env python3
"""Pathway enrichment for candidate genes."""
import pandas as pd
from pathlib import Path

# Load candidates
candidates = pd.read_csv("out_v3/candidate_genes_fdr01.csv")
gene_list = candidates['gene'].tolist()

print(f"Analyzing {len(gene_list)} candidate genes")

# Option 1: Use g:Profiler (web-based)
print("\n=== Option 1: g:Profiler (Web-based) ===")
print("Gene list for g:Profiler:")
print("\n".join(gene_list))
print("\n→ Copy genes above and paste into https://biit.cs.ut.ee/gprofiler/gost")

# Option 2: Use Enrichr (Python API)
print("\n=== Option 2: Enrichr (Python API) ===")
try:
    from gseapy import enrichr
    
    # Run enrichment
    enr = enrichr(
        gene_list=gene_list,
        gene_sets=['GO_Biological_Process_2021', 'KEGG_2021_Human', 'WikiPathway_2021_Human'],
        organism='human',
        outdir='out_v3/enrichment',
        cutoff=0.05
    )
    
    print("\n✓ Enrichment results saved to out_v3/enrichment/")
    
    # Print top results
    results = pd.read_csv("out_v3/enrichment/GO_Biological_Process_2021.txt", sep="\t")
    print("\nTop 10 enriched pathways:")
    print(results[['Term', 'Adjusted P-value', 'Genes']].head(10))
    
except ImportError:
    print("\n⚠️  gseapy not installed. Install with: pip install gseapy")
    print("Or use web version: https://maayanlab.cloud/Enrichr/")
    
    # Save gene list for manual upload
    Path("out_v3/enrichment").mkdir(parents=True, exist_ok=True)
    with open("out_v3/enrichment/gene_list.txt", "w") as f:
        f.write("\n".join(gene_list))
    print(f"\n✓ Saved gene list to out_v3/enrichment/gene_list.txt")
    print("   Upload this file to Enrichr web interface")

# Option 3: Simple GO term lookup (if available)
print("\n=== Option 3: Manual Pathway Check ===")
print("Expected enrichment (if hypothesis is correct):")
print("  - DNA repair (GO:0006281)")
print("  - Response to DNA damage (GO:0006974)")
print("  - Double-strand break repair (GO:0006302)")
print("  - Chromatin remodeling (GO:0006338)")
print("  - Cell cycle checkpoint (GO:0000075)")
print("\nIf you see enrichment → Strong support for mechanism")
print("If no enrichment → Effect may be non-specific spatial artifact")

