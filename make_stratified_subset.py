#!/usr/bin/env python3
"""
Create a stratified gene subset based on CN variance × proximity prevalence.

Forms a 10×10 grid of deciles and samples ~20 genes per cell (up to 2k total).
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path


def compute_gene_stats(dep_long: pd.DataFrame, 
                      genes_bed: pd.DataFrame,
                      sv_bedpe: pd.DataFrame,
                      cnv_bed: pd.DataFrame,
                      features_cache_dir: Path = None) -> pd.DataFrame:
    """
    Compute per-gene CN variance and proximity prevalence.
    
    Args:
        dep_long: Long-format dependency data (gene, cell_line, dependency)
        genes_bed: Gene BED file
        sv_bedpe: SV breakpoints
        cnv_bed: CNV segments
        features_cache_dir: Optional directory with cached features
        
    Returns:
        DataFrame with gene, cn_variance, proximity_prevalence
    """
    # Try to load from cache first
    if features_cache_dir:
        cache_dir = Path(features_cache_dir)
        cn_cache = cache_dir / "features_cn.parquet"
        bp_cache = cache_dir / "features_bp.parquet"
        
        if cn_cache.exists() and bp_cache.exists():
            print(f"Loading features from cache: {cache_dir}")
            gene_cn = pd.read_parquet(cn_cache)
            bp_dist = pd.read_parquet(bp_cache)
        else:
            print(f"Cache not found, computing features...")
            gene_cn, bp_dist = compute_features(dep_long, genes_bed, sv_bedpe, cnv_bed)
    else:
        gene_cn, bp_dist = compute_features(dep_long, genes_bed, sv_bedpe, cnv_bed)
    
    # Merge with dependency to get full data
    dep_merged = dep_long.merge(gene_cn, on=["gene", "cell_line"], how="left")
    dep_merged = dep_merged.merge(bp_dist, on=["gene", "cell_line"], how="left")
    
    # Compute per-gene CN variance
    cn_var = dep_merged.groupby("gene")["cn"].var()
    
    # Compute proximity prevalence (p<=1Mb)
    dep_merged["prox_active"] = (dep_merged["bp_dist"] <= 1_000_000).astype(float)
    prox_prev = dep_merged.groupby("gene")["prox_active"].mean()
    
    # Combine
    stats = pd.DataFrame({
        "gene": cn_var.index,
        "cn_variance": cn_var.values,
        "proximity_prevalence": prox_prev.values
    })
    
    # Fill NaN with 0
    stats = stats.fillna(0)
    
    return stats


def compute_features(dep_long: pd.DataFrame,
                    genes_bed: pd.DataFrame,
                    sv_bedpe: pd.DataFrame,
                    cnv_bed: pd.DataFrame):
    """
    Compute features (simplified version - can use pipeline's FeatureComputer if available).
    For now, we'll do a basic computation.
    """
    # This is a simplified version - in practice, you'd use the pipeline's FeatureComputer
    # For now, we'll compute basic features
    
    # Gene CN: overlap CNV segments with genes
    gene_cn_list = []
    for _, gene in genes_bed.iterrows():
        chrom = gene["chrom"]
        g_start = gene["start"]
        g_end = gene["end"]
        g_name = gene["gene"]
        
        # Find overlapping CNV segments
        cnv_chr = cnv_bed[cnv_bed["chrom"] == chrom]
        for _, cnv in cnv_chr.iterrows():
            if not (cnv["end"] < g_start or cnv["start"] > g_end):
                overlap_start = max(cnv["start"], g_start)
                overlap_end = min(cnv["end"], g_end)
                overlap_len = max(0, overlap_end - overlap_start)
                gene_len = max(1, g_end - g_start)
                
                weighted_cn = (overlap_len / gene_len) * cnv["cn"]
                gene_cn_list.append({
                    "gene": g_name,
                    "cell_line": cnv["cell_line"],
                    "cn": weighted_cn
                })
    
    gene_cn = pd.DataFrame(gene_cn_list)
    if len(gene_cn) > 0:
        gene_cn = gene_cn.groupby(["gene", "cell_line"])["cn"].sum().reset_index()
    
    # BP distances: minimum distance to breakpoint
    bp_list = []
    for _, sv in sv_bedpe.iterrows():
        for side in ["1", "2"]:
            chrom_col = f"chrom{side}"
            start_col = f"start{side}"
            end_col = f"end{side}"
            bp_pos = (sv[start_col] + sv[end_col]) // 2
            
            genes_chr = genes_bed[genes_bed["chrom"] == sv[chrom_col]]
            for _, gene in genes_chr.iterrows():
                gene_mid = (gene["start"] + gene["end"]) // 2
                dist = abs(gene_mid - bp_pos)
                bp_list.append({
                    "gene": gene["gene"],
                    "cell_line": sv["cell_line"],
                    "bp_dist": min(dist, 5_000_000)  # Cap at 5Mb
                })
    
    bp_dist = pd.DataFrame(bp_list)
    if len(bp_dist) > 0:
        bp_dist = bp_dist.groupby(["gene", "cell_line"])["bp_dist"].min().reset_index()
    
    return gene_cn, bp_dist


def create_stratified_subset(stats: pd.DataFrame, 
                             target_size: int = 2000,
                             genes_per_cell: int = 20) -> list:
    """
    Create stratified subset using 10×10 grid of deciles.
    
    Args:
        stats: DataFrame with gene, cn_variance, proximity_prevalence
        target_size: Target total number of genes
        genes_per_cell: Target genes per grid cell
        
    Returns:
        List of gene symbols
    """
    # Remove genes with NaN or zero variance/prevalence (if needed)
    stats = stats[(stats["cn_variance"] > 0) | (stats["proximity_prevalence"] > 0)].copy()
    
    # Compute deciles
    cn_deciles = pd.qcut(stats["cn_variance"], q=10, labels=False, duplicates="drop")
    prox_deciles = pd.qcut(stats["proximity_prevalence"], q=10, labels=False, duplicates="drop")
    
    stats["cn_decile"] = cn_deciles
    stats["prox_decile"] = prox_deciles
    
    # Sample from each cell
    selected_genes = []
    
    for cn_d in range(10):
        for prox_d in range(10):
            cell = stats[(stats["cn_decile"] == cn_d) & (stats["prox_decile"] == prox_d)]
            if len(cell) == 0:
                continue
            
            n_sample = min(genes_per_cell, len(cell), target_size - len(selected_genes))
            if n_sample <= 0:
                break
            
            sampled = cell.sample(n=n_sample, random_state=42)
            selected_genes.extend(sampled["gene"].tolist())
            
            if len(selected_genes) >= target_size:
                break
        if len(selected_genes) >= target_size:
            break
    
    return selected_genes[:target_size]


def main():
    parser = argparse.ArgumentParser(
        description="Create stratified gene subset based on CN variance × proximity prevalence",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--dependency-long", required=True,
                       help="Long-format dependency CSV (gene, cell_line, dependency)")
    parser.add_argument("--genes-bed", required=True,
                       help="Gene BED file (TSV: chrom, start, end, gene)")
    parser.add_argument("--sv-bedpe", required=True,
                       help="SV breakpoints (BEDPE-like TSV)")
    parser.add_argument("--cnv-bed", required=True,
                       help="CNV segments (BED-like TSV)")
    parser.add_argument("--features-cache", type=str, default=None,
                       help="Optional directory with cached features (features_cn.parquet, features_bp.parquet)")
    parser.add_argument("--out", default="genes_stratified_2k.txt",
                       help="Output file with gene list")
    parser.add_argument("--target-size", type=int, default=2000,
                       help="Target number of genes")
    parser.add_argument("--genes-per-cell", type=int, default=20,
                       help="Target genes per grid cell")
    
    args = parser.parse_args()
    
    # Load data
    print("Loading data...")
    dep_long = pd.read_csv(args.dependency_long)
    genes_bed = pd.read_csv(args.genes_bed, sep="\t", header=None,
                           names=["chrom", "start", "end", "gene", "strand"])
    sv_bedpe = pd.read_csv(args.sv_bedpe, sep="\t")
    cnv_bed = pd.read_csv(args.cnv_bed, sep="\t")
    
    # Compute gene statistics
    print("Computing gene statistics...")
    stats = compute_gene_stats(
        dep_long, genes_bed, sv_bedpe, cnv_bed,
        features_cache_dir=args.features_cache
    )
    
    print(f"Computed stats for {len(stats):,} genes")
    print(f"  CN variance range: [{stats['cn_variance'].min():.4f}, {stats['cn_variance'].max():.4f}]")
    print(f"  Proximity prevalence range: [{stats['proximity_prevalence'].min():.4f}, {stats['proximity_prevalence'].max():.4f}]")
    
    # Create stratified subset
    print(f"Creating stratified subset (target: {args.target_size} genes)...")
    selected = create_stratified_subset(stats, target_size=args.target_size, 
                                       genes_per_cell=args.genes_per_cell)
    
    print(f"Selected {len(selected):,} genes")
    
    # Save
    with open(args.out, 'w') as f:
        for gene in sorted(selected):
            f.write(f"{gene}\n")
    
    print(f"Saved to {args.out}")


if __name__ == "__main__":
    main()

