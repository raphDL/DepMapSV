#!/usr/bin/env python3
"""
Coefficient→effect consistency per line.

For each shortlisted line, test Spearman correlation between:
- per-gene |β_prox|max (max absolute proximity coefficient)
- per-gene corr-drop (Δ|corr|) within that line

TRUE should show a more positive ρ than ROTATE.
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path

def load_data(base_path):
    """Load design matrix and coefficients."""
    dm_path = Path(base_path) / "design_matrix.csv"
    coef_path = Path(base_path) / "models_coefficients.csv"
    
    dm = pd.read_csv(dm_path, usecols=["gene","cell_line","dependency","dependency_corrected","cn"])
    coef = pd.read_csv(coef_path, usecols=["gene","bp_within_100000","bp_within_1000000"])
    
    # Compute max absolute proximity coefficient per gene
    coef["abs_coef_max"] = coef[["bp_within_100000","bp_within_1000000"]].abs().max(axis=1)
    
    # Merge
    df = dm.merge(coef[["gene","abs_coef_max"]], on="gene", how="left")
    
    return df

def compute_per_gene_corr_drop(df, cell_line):
    """Compute corr-drop per gene within a cell line."""
    sub = df[df["cell_line"] == cell_line].copy()
    if len(sub) < 10:
        return None
    
    rows = []
    for gene, gsub in sub.groupby("gene"):
        if len(gsub) < 1:  # Usually 1 row per gene×line, but check
            continue
        
        # For per-gene, we need multiple cell lines, so compute across all lines
        # Actually, within a line, we have one value per gene
        # So we compute across genes in that line
        # But we need the gene-level corr-drop which is computed across all lines
        # Let's compute it properly: corr-drop is computed per gene across all lines
        pass
    
    # Actually, for per-line analysis, we compute:
    # For each gene in this line, what's its corr-drop across all lines?
    # Then correlate with |coef|max for that gene
    return None

def compute_line_coef_effect_rho(df, cell_line):
    """
    Compute Spearman ρ between |coef|max and corr-drop for genes in this line.
    
    For each gene, we need:
    1. |coef|max (already in df)
    2. corr-drop = |corr(dep, CN)| - |corr(dep_corr, CN)| computed across all cell lines
    """
    # Get genes in this line
    line_genes = set(df[df["cell_line"] == cell_line]["gene"].unique())
    
    # Compute corr-drop per gene (across all lines)
    gene_corr_drops = []
    for gene in line_genes:
        gsub = df[df["gene"] == gene][["dependency","dependency_corrected","cn"]].dropna()
        if len(gsub) < 20:
            continue
        
        try:
            r_before, _ = spearmanr(gsub["dependency"], gsub["cn"])
            r_after, _ = spearmanr(gsub["dependency_corrected"], gsub["cn"])
            if np.isfinite(r_before) and np.isfinite(r_after):
                corr_drop = abs(r_before) - abs(r_after)
                coef_max = df[df["gene"] == gene]["abs_coef_max"].iloc[0] if len(df[df["gene"] == gene]) > 0 else np.nan
                if np.isfinite(coef_max):
                    gene_corr_drops.append({
                        "gene": gene,
                        "abs_coef_max": coef_max,
                        "corr_drop": corr_drop
                    })
        except Exception:
            continue
    
    if len(gene_corr_drops) < 10:
        return None
    
    gdf = pd.DataFrame(gene_corr_drops)
    rho, pval = spearmanr(gdf["abs_coef_max"], gdf["corr_drop"])
    
    return {
        "cell_line": cell_line,
        "spearman_rho": rho,
        "p_value": pval,
        "n_genes": len(gdf)
    }

def main():
    base_true = "out_focus_true/unstable"
    base_rot = "out_focus_rotate/unstable"
    
    # Load shortlist
    shortlist_path = Path("figs_full/psi_line_shortlist_for_doc.csv")
    if not shortlist_path.exists():
        print("Shortlist not found. Using top 5 lines from EB results...")
        eb = pd.read_csv("figs_full/psi_line_eb_full.csv")
        shortlist = eb.nlargest(5, "post_diff_mean")["cell_line"].tolist()
    else:
        shortlist = pd.read_csv(shortlist_path)["cell_line"].head(5).tolist()
    
    print(f"Analyzing {len(shortlist)} shortlisted lines...")
    
    # Load data
    print("Loading TRUE data...")
    df_true = load_data(base_true)
    
    print("Loading ROTATE data...")
    df_rot = load_data(base_rot)
    
    # Compute for each line
    results = []
    for cl in shortlist:
        print(f"  Processing {cl}...")
        true_result = compute_line_coef_effect_rho(df_true, cl)
        rot_result = compute_line_coef_effect_rho(df_rot, cl)
        
        if true_result and rot_result:
            results.append({
                "cell_line": cl,
                "rho_TRUE": true_result["spearman_rho"],
                "pval_TRUE": true_result["p_value"],
                "n_genes_TRUE": true_result["n_genes"],
                "rho_ROTATE": rot_result["spearman_rho"],
                "pval_ROTATE": rot_result["p_value"],
                "n_genes_ROTATE": rot_result["n_genes"],
                "rho_diff": true_result["spearman_rho"] - rot_result["spearman_rho"]
            })
    
    if not results:
        print("No results computed.")
        return
    
    out = pd.DataFrame(results)
    outdir = Path("figs_full")
    outdir.mkdir(exist_ok=True)
    out.to_csv(outdir / "psi_line_coef_effect.csv", index=False)
    
    print("\n" + "="*60)
    print("COEFFICIENT→EFFECT CONSISTENCY")
    print("="*60)
    print(out.to_string(index=False))
    print(f"\nMedian ρ (TRUE): {out['rho_TRUE'].median():.4f}")
    print(f"Median ρ (ROTATE): {out['rho_ROTATE'].median():.4f}")
    print(f"Median ρ difference: {out['rho_diff'].median():.4f}")
    
    if out['rho_diff'].median() > 0:
        print("\n✓ TRUE shows stronger coefficient→effect consistency")
    else:
        print("\n✗ No evidence TRUE > ROTATE")
    
    print(f"\nSaved: {outdir}/psi_line_coef_effect.csv")

if __name__ == "__main__":
    import sys
    main()

