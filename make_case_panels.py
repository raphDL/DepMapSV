#!/usr/bin/env python3
"""
Make Case Panels - Generate 4-panel figures for each case study gene.

Panels:
A) Structural context (CN segments + WGS SV breakpoints ±1Mb)
B) CN vs dependency scatter (before/after correction) with Δ|corr|
C) Distance vs dependency (log-x) stratified by CN tertiles
D) Orthogonal validation mini-panel (RNAi/Drug/pathway)
"""
import argparse
import logging
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Set matplotlib backend
plt.switch_backend('Agg')


def setup_logging(verbosity: int = 1):
    """Setup logging with verbosity control."""
    level = logging.WARNING if verbosity <= 0 else (logging.INFO if verbosity == 1 else logging.DEBUG)
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S",
        handlers=[logging.StreamHandler(sys.stderr)],
        force=True,
    )


def load_gene_data(gene: str, true_dir: str, sv_dir: str, cnv_path: str, genes_bed: str = None) -> tuple:
    """Load all data needed for a gene panel."""
    # Load design matrix
    design_path = os.path.join(true_dir, "design_matrix_wgs.csv")
    design = pd.read_csv(design_path)
    gene_design = design[design["gene"] == gene].copy()
    
    # Load corrected dependencies
    corr_path = os.path.join(true_dir, "dependency_corrected.csv")
    corrected = pd.read_csv(corr_path)
    gene_corr = corrected[corrected["gene"] == gene].copy()
    
    # Merge
    gene_data = gene_design.merge(
        gene_corr[["cell_line", "dependency", "prox_only_corrected"]],
        on="cell_line", how="left"
    )
    
    # Load CNV segments for this gene's region
    cnv = pd.read_csv(cnv_path, sep="\t")
    # Try to get gene coordinates from design matrix or load from genes BED
    gene_info = gene_design.iloc[0] if len(gene_design) > 0 else None
    chrom = None
    gene_start = None
    gene_end = None
    
    if gene_info is not None:
        # Try design matrix columns first
        if "chrom" in gene_design.columns:
            chrom = gene_info["chrom"]
        if "start" in gene_design.columns:
            gene_start = gene_info["start"]
        if "end" in gene_design.columns:
            gene_end = gene_info["end"]
    
    # If not in design matrix, try loading from genes BED file
    if (chrom is None or gene_start is None or gene_end is None) and genes_bed:
        try:
            genes_df = pd.read_csv(genes_bed, sep="\t", header=None,
                                  names=["chrom", "start", "end", "gene", "strand"])
            gene_row = genes_df[genes_df["gene"] == gene]
            if len(gene_row) > 0:
                chrom = gene_row["chrom"].iloc[0]
                gene_start = int(gene_row["start"].iloc[0])
                gene_end = int(gene_row["end"].iloc[0])
        except Exception as e:
            logging.debug(f"Could not load gene coordinates from {genes_bed}: {e}")
    
    # Get CNV region if we have coordinates
    if chrom and gene_start and gene_end:
        # Get region ±1Mb
        region_start = max(0, gene_start - 1_000_000)
        region_end = gene_end + 1_000_000
        cnv_region = cnv[
            (cnv["chrom"] == chrom) &
            (cnv["start"] < region_end) &
            (cnv["end"] > region_start)
        ].copy()
    else:
        cnv_region = pd.DataFrame()
    
    # Load SV breakpoints for this gene's region
    sv_region = pd.DataFrame()
    if chrom and gene_start and gene_end:
        region_start = max(0, gene_start - 1_000_000)
        region_end = gene_end + 1_000_000
        sv_path = Path(sv_dir)
        if sv_path.is_file():
            sv_files = [sv_path]
        else:
            sv_files = list(sv_path.glob("*.bedpe"))
        
        for fpath in sv_files:
            try:
                sv = pd.read_csv(fpath, sep="\t")
                if "chrom1" in sv.columns:
                    sv_sub = sv[
                        ((sv["chrom1"] == chrom) & (sv["start1"] < region_end) & (sv["end1"] > region_start)) |
                        ((sv["chrom2"] == chrom) & (sv["start2"] < region_end) & (sv["end2"] > region_start))
                    ]
                    sv_region = pd.concat([sv_region, sv_sub], ignore_index=True)
            except Exception as e:
                logging.debug(f"Error loading {fpath}: {e}")
    
    return gene_data, cnv_region, sv_region


def panel_a_structural_context(gene_data: pd.DataFrame, cnv_region: pd.DataFrame,
                               sv_region: pd.DataFrame, ax):
    """Panel A: Structural context (CN segments + WGS SV breakpoints ±1Mb)."""
    if len(gene_data) == 0:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("A) Structural Context")
        return
    
    # Get gene coordinates
    if "start" in gene_data.columns and "end" in gene_data.columns:
        gene_start = gene_data["start"].iloc[0]
        gene_end = gene_data["end"].iloc[0]
        region_start = max(0, gene_start - 1_000_000)
        region_end = gene_end + 1_000_000
    else:
        ax.text(0.5, 0.5, "No coordinates", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("A) Structural Context")
        return
    
    # Plot CNV segments
    if len(cnv_region) > 0:
        for _, seg in cnv_region.iterrows():
            ax.barh(0, seg["end"] - seg["start"], left=seg["start"], 
                   height=0.3, alpha=0.5, color="blue")
    
    # Plot SV breakpoints
    if len(sv_region) > 0:
        for _, sv in sv_region.iterrows():
            if "start1" in sv:
                ax.scatter(sv["start1"], 0.5, marker="|", s=100, color="red", alpha=0.7)
            if "start2" in sv:
                ax.scatter(sv["start2"], 0.5, marker="|", s=100, color="red", alpha=0.7)
    
    # Highlight gene region
    ax.axvspan(gene_start, gene_end, alpha=0.2, color="green", label="Gene")
    
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(-0.5, 1.0)
    ax.set_xlabel("Genomic Position (bp)")
    ax.set_ylabel("")
    ax.set_title("A) Structural Context")
    ax.legend()


def panel_b_cn_vs_dependency(gene_data: pd.DataFrame, ax):
    """Panel B: CN vs dependency scatter (before/after correction) with Δ|corr|."""
    if len(gene_data) == 0:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("B) CN vs Dependency")
        return
    
    # Before correction
    before = gene_data[["cn", "dependency"]].dropna()
    if len(before) > 0:
        ax.scatter(before["cn"], before["dependency"], alpha=0.5, s=20, 
                  label="Before", color="gray")
        r_before = before.corr().iloc[0, 1]
    
    # After correction
    after = gene_data[["cn", "prox_only_corrected"]].dropna()
    if len(after) > 0:
        ax.scatter(after["cn"], after["prox_only_corrected"], alpha=0.5, s=20,
                  label="After prox-only", color="green")
        r_after = after.corr().iloc[0, 1]
        
        delta_corr = abs(r_before) - abs(r_after) if len(before) > 0 else 0.0
        ax.text(0.05, 0.95, f"Δ|corr| = {delta_corr:.3f}", 
               transform=ax.transAxes, va="top", ha="left",
               bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    
    ax.set_xlabel("Copy Number")
    ax.set_ylabel("Dependency")
    ax.set_title("B) CN vs Dependency")
    ax.legend()


def panel_c_distance_vs_dependency(gene_data: pd.DataFrame, ax):
    """Panel C: Distance vs dependency (log-x) stratified by CN tertiles."""
    if len(gene_data) == 0:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("C) Distance vs Dependency")
        return
    
    # Stratify by CN tertiles
    cn_vals = gene_data["cn"].dropna()
    if len(cn_vals) < 3:
        ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("C) Distance vs Dependency")
        return
    
    tertiles = cn_vals.quantile([0.33, 0.67])
    gene_data["cn_tertile"] = pd.cut(gene_data["cn"], 
                                     bins=[-np.inf, tertiles.iloc[0], tertiles.iloc[1], np.inf],
                                     labels=["Low", "Mid", "High"])
    
    colors = {"Low": "blue", "Mid": "orange", "High": "red"}
    for tertile in ["Low", "Mid", "High"]:
        sub = gene_data[gene_data["cn_tertile"] == tertile]
        if len(sub) > 0:
            ax.scatter(sub["bp_dist"], sub["dependency"], 
                      alpha=0.5, s=20, label=f"CN {tertile}", color=colors[tertile])
    
    ax.set_xscale("log")
    ax.set_xlabel("Breakpoint Distance (bp)")
    ax.set_ylabel("Dependency")
    ax.set_title("C) Distance vs Dependency")
    ax.legend()


def panel_d_orthogonal_validation(gene: str, validated_data: pd.DataFrame, ax):
    """Panel D: Orthogonal validation mini-panel (RNAi/Drug/pathway)."""
    if gene not in validated_data["gene"].values:
        ax.text(0.5, 0.5, "No validation data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("D) Orthogonal Validation")
        return
    
    row = validated_data[validated_data["gene"] == gene].iloc[0]
    
    # Create simple bar chart of validation status
    validations = {
        "RNAi": row.get("rnai_validated", False),
        "Drug": row.get("drug_validated", False),
        "Pathway": row.get("pathway_validated", False)
    }
    
    colors = ["green" if v else "red" for v in validations.values()]
    ax.barh(list(validations.keys()), [1 if v else 0 for v in validations.values()],
           color=colors, alpha=0.7)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Validated")
    ax.set_title("D) Orthogonal Validation")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["No", "Yes"])


def make_panel(gene: str, true_dir: str, sv_dir: str, cnv_path: str,
               validated_data: pd.DataFrame, out_path: str, genes_bed: str = None):
    """Create 4-panel figure for a gene."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    # Load data
    gene_data, cnv_region, sv_region = load_gene_data(gene, true_dir, sv_dir, cnv_path, genes_bed)
    
    # Create panels
    panel_a_structural_context(gene_data, cnv_region, sv_region, axes[0])
    panel_b_cn_vs_dependency(gene_data, axes[1])
    panel_c_distance_vs_dependency(gene_data, axes[2])
    panel_d_orthogonal_validation(gene, validated_data, axes[3])
    
    plt.suptitle(f"Case Study: {gene}", fontsize=14, fontweight="bold")
    plt.tight_layout()
    
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    logging.info(f"Saved panel for {gene} to {out_path}")


def main():
    ap = argparse.ArgumentParser(
        description="Generate 4-panel figures for case studies"
    )
    ap.add_argument("--cases", required=True,
                    help="Case studies CSV (with validation columns)")
    ap.add_argument("--true", required=True,
                    help="TRUE output directory")
    ap.add_argument("--sv-dir", required=True,
                    help="WGS SV directory or file")
    ap.add_argument("--cnv", required=True,
                    help="CNV segments BED file")
    ap.add_argument("--genes", default=None,
                    help="Genes BED file (for gene coordinates if not in design matrix)")
    ap.add_argument("--out", required=True,
                    help="Output directory for panels")
    ap.add_argument("--max-genes", type=int, default=None,
                    help="Maximum number of genes to plot (default: all)")
    ap.add_argument("-v", "--verbose", action="count", default=1,
                    help="Increase verbosity")
    args = ap.parse_args()
    
    setup_logging(args.verbose)
    
    os.makedirs(args.out, exist_ok=True)
    
    # Load case studies
    logging.info("Loading case studies...")
    case_studies = pd.read_csv(args.cases)
    
    if args.max_genes:
        case_studies = case_studies.head(args.max_genes)
    
    logging.info(f"Generating panels for {len(case_studies)} genes...")
    
    for _, row in case_studies.iterrows():
        gene = row["gene"]
        out_path = os.path.join(args.out, f"{gene}.png")
        
        try:
            make_panel(gene, args.true, args.sv_dir, args.cnv, case_studies, out_path, args.genes)
        except Exception as e:
            logging.error(f"Error creating panel for {gene}: {e}")
    
    logging.info(f"Panels saved to {args.out}")


if __name__ == "__main__":
    main()

