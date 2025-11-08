#!/usr/bin/env python3
"""
Per-line Δ|corr| across genes (CN-trimmed), TRUE vs ROTATE.

For each cell line:
- Take all genes from dependency_corrected.csv
- Keep genes with |CN-2| ≤ 0.5 (or trim to 10th-90th CN percentiles)
- Compute Spearman |corr(dep, CN)| before and after prox-only correction
- Record Δ|corr| = |r_before| - |r_after|
- Compare TRUE vs ROTATE with paired bootstrap CI
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path
from numpy.random import default_rng

def load_dependency_corrected(base_path):
    """Load dependency_corrected.csv and compute prox-only correction."""
    dep_path = Path(base_path) / "dependency_corrected.csv"
    dm_path = Path(base_path) / "design_matrix.csv"
    
    # Try design_matrix.csv first (has all we need)
    if dm_path.exists():
        df = pd.read_csv(dm_path)
        # Check what columns we have
        needed = ["gene","cell_line","dependency"]
        if not all(c in df.columns for c in needed):
            raise ValueError(f"Missing required columns in design_matrix.csv: {needed}")
        
        # Get CN - might be in design_matrix or need to merge
        if "cn" not in df.columns:
            # Try to get from dependency_corrected if it exists
            if dep_path.exists():
                dep_df = pd.read_csv(dep_path, usecols=["gene","cell_line","cn"])
                df = df.merge(dep_df, on=["gene","cell_line"], how="left")
            else:
                raise ValueError("No 'cn' column found and dependency_corrected.csv not available")
        
        # Compute prox_contrib from coefficients
        coef_path = Path(base_path) / "models_coefficients.csv"
        if coef_path.exists():
            coef = pd.read_csv(coef_path, usecols=["gene","bp_within_100000","bp_within_1000000"])
            coef = coef.rename(columns={"bp_within_100000":"coef_bp_near", "bp_within_1000000":"coef_bp_far"})
            if "bp_dist" in df.columns:
                w1, w2 = 250_000, 2_000_000
                df["bp_near"] = (df["bp_dist"].fillna(5_000_000) <= w1).astype(float)
                df["bp_far"] = ((df["bp_dist"].fillna(5_000_000) > w1) & (df["bp_dist"] <= w2)).astype(float)
                df = df.merge(coef, on="gene", how="left")
                df["prox_contrib"] = df["coef_bp_near"].fillna(0)*df["bp_near"] + df["coef_bp_far"].fillna(0)*df["bp_far"]
                df["prox_only_corrected"] = df["dependency"] - df["prox_contrib"].fillna(0)
            else:
                # Use dependency_corrected as proxy
                if "dependency_corrected" in df.columns:
                    df["prox_only_corrected"] = df["dependency_corrected"]
                else:
                    df["prox_only_corrected"] = df["dependency"]  # fallback
        else:
            if "dependency_corrected" in df.columns:
                df["prox_only_corrected"] = df["dependency_corrected"]
            else:
                df["prox_only_corrected"] = df["dependency"]  # fallback
        return df[["gene","cell_line","dependency","prox_only_corrected","cn"]].copy()
    
    # Fallback to dependency_corrected.csv
    if dep_path.exists():
        df = pd.read_csv(dep_path)
        # Compute prox_only_corrected if needed
        if "prox_only_corrected" not in df.columns:
            coef_path = Path(base_path) / "models_coefficients.csv"
            if coef_path.exists() and "bp_dist" in df.columns:
                coef = pd.read_csv(coef_path, usecols=["gene","bp_within_100000","bp_within_1000000"])
                coef = coef.rename(columns={"bp_within_100000":"coef_bp_near", "bp_within_1000000":"coef_bp_far"})
                w1, w2 = 250_000, 2_000_000
                df["bp_near"] = (df["bp_dist"].fillna(5_000_000) <= w1).astype(float)
                df["bp_far"] = ((df["bp_dist"].fillna(5_000_000) > w1) & (df["bp_dist"] <= w2)).astype(float)
                df = df.merge(coef, on="gene", how="left")
                df["prox_contrib"] = df["coef_bp_near"].fillna(0)*df["bp_near"] + df["coef_bp_far"].fillna(0)*df["bp_far"]
                df["prox_only_corrected"] = df["dependency"] - df["prox_contrib"]
            else:
                df["prox_only_corrected"] = df.get("dependency_corrected", df["dependency"])
        return df[["gene","cell_line","dependency","prox_only_corrected","cn"]].copy()
    
    raise FileNotFoundError(f"Neither dependency_corrected.csv nor design_matrix.csv found in {base_path}")

def compute_line_corr_delta(df, cn_trim_method="abs", cn_thresh=0.5):
    """
    Compute per-line Δ|corr| across genes.
    
    Args:
        df: DataFrame with gene, cell_line, dependency, prox_only_corrected, cn
        cn_trim_method: "abs" (|CN-2| ≤ thresh) or "percentile" (10th-90th)
        cn_thresh: threshold for |CN-2| or None for percentile
    """
    rows = []
    for cl, sub in df.groupby("cell_line"):
        # CN trimming
        if cn_trim_method == "abs":
            sub = sub[sub["cn"].sub(2.0).abs() <= cn_thresh].copy()
        elif cn_trim_method == "percentile":
            p10, p90 = sub["cn"].quantile([0.1, 0.9])
            sub = sub[sub["cn"].between(p10, p90)].copy()
        
        if len(sub) < 20:  # Need enough genes
            continue
        
        # Compute correlations across genes (not within genes)
        # This is per-line: corr(dependency, CN) across all genes in that line
        before = sub[["dependency", "cn"]].dropna()
        after = sub[["prox_only_corrected", "cn"]].dropna()
        
        if len(before) >= 20 and len(after) >= 20:
            try:
                r_before, _ = spearmanr(before["dependency"], before["cn"])
                r_after, _ = spearmanr(after["prox_only_corrected"], after["cn"])
                
                if np.isfinite(r_before) and np.isfinite(r_after):
                    delta = abs(r_before) - abs(r_after)
                    rows.append({
                        "cell_line": cl,
                        "r_before": r_before,
                        "r_after": r_after,
                        "abs_r_before": abs(r_before),
                        "abs_r_after": abs(r_after),
                        "delta_abs_corr": delta,
                        "n_genes": len(sub)
                    })
            except Exception:
                continue
    
    return pd.DataFrame(rows)

def bootstrap_paired_diff(true_deltas, rot_deltas, n_reps=5000, seed=0):
    """Bootstrap CI for paired difference (TRUE - ROTATE)."""
    rng = default_rng(seed)
    
    # Align by cell_line
    merged = pd.merge(
        true_deltas[["cell_line", "delta_abs_corr"]].rename(columns={"delta_abs_corr": "true"}),
        rot_deltas[["cell_line", "delta_abs_corr"]].rename(columns={"delta_abs_corr": "rot"}),
        on="cell_line",
        how="inner"
    )
    
    if len(merged) == 0:
        return {"median_diff": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "n_pairs": 0}
    
    diffs = merged["true"] - merged["rot"]
    median_diff = diffs.median()
    
    # Bootstrap
    boot_diffs = []
    for _ in range(n_reps):
        idx = rng.choice(len(merged), size=len(merged), replace=True)
        boot_diffs.append(diffs.iloc[idx].median())
    
    boot_diffs = np.array(boot_diffs)
    ci_lo, ci_hi = np.percentile(boot_diffs, [2.5, 97.5])
    
    return {
        "median_diff": median_diff,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "n_pairs": len(merged),
        "median_true": merged["true"].median(),
        "median_rot": merged["rot"].median()
    }

def main():
    base_true = "out_focus_true/unstable"
    base_rot = "out_focus_rotate/unstable"
    outdir = Path("figs_full")
    outdir.mkdir(exist_ok=True)
    
    print("Loading TRUE data...")
    df_true = load_dependency_corrected(base_true)
    
    print("Loading ROTATE data...")
    df_rot = load_dependency_corrected(base_rot)
    
    print("Computing per-line Δ|corr| (TRUE)...")
    true_deltas = compute_line_corr_delta(df_true, cn_trim_method="abs", cn_thresh=0.5)
    
    print("Computing per-line Δ|corr| (ROTATE)...")
    rot_deltas = compute_line_corr_delta(df_rot, cn_trim_method="abs", cn_thresh=0.5)
    
    # Bootstrap comparison
    print("Computing paired bootstrap CI...")
    boot_result = bootstrap_paired_diff(true_deltas, rot_deltas)
    
    # Save results
    true_deltas["condition"] = "TRUE"
    rot_deltas["condition"] = "ROTATE"
    combined = pd.concat([true_deltas, rot_deltas], ignore_index=True)
    combined.to_csv(outdir / "psi_line_corr_delta.csv", index=False)
    
    # Summary
    summary = pd.DataFrame([{
        "metric": "median_delta_TRUE",
        "value": true_deltas["delta_abs_corr"].median()
    }, {
        "metric": "median_delta_ROTATE",
        "value": rot_deltas["delta_abs_corr"].median()
    }, {
        "metric": "paired_median_diff",
        "value": boot_result["median_diff"]
    }, {
        "metric": "paired_ci_lo",
        "value": boot_result["ci_lo"]
    }, {
        "metric": "paired_ci_hi",
        "value": boot_result["ci_hi"]
    }, {
        "metric": "n_pairs",
        "value": boot_result["n_pairs"]
    }])
    summary.to_csv(outdir / "psi_line_corr_delta_summary.csv", index=False)
    
    print("\n" + "="*60)
    print("PER-LINE Δ|corr| RESULTS (CN-trimmed, |CN-2| ≤ 0.5)")
    print("="*60)
    print(f"TRUE median Δ|corr|: {true_deltas['delta_abs_corr'].median():.4f}")
    print(f"ROTATE median Δ|corr|: {rot_deltas['delta_abs_corr'].median():.4f}")
    print(f"\nPaired difference (TRUE - ROTATE):")
    print(f"  Median: {boot_result['median_diff']:.4f}")
    print(f"  95% CI: [{boot_result['ci_lo']:.4f}, {boot_result['ci_hi']:.4f}]")
    print(f"  N pairs: {boot_result['n_pairs']}")
    
    if boot_result['ci_lo'] > 0:
        print("\n✓ STRONG: CI excludes 0, TRUE > ROTATE")
    elif boot_result['median_diff'] > 0:
        print("\n⚠ MODERATE: Median positive but CI includes 0")
    else:
        print("\n✗ WEAK: No evidence TRUE > ROTATE")
    
    print(f"\nSaved: {outdir}/psi_line_corr_delta.csv")
    print(f"Saved: {outdir}/psi_line_corr_delta_summary.csv")

if __name__ == "__main__":
    main()

