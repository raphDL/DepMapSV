#!/usr/bin/env python3
"""
Compute partial correlation metric from existing pipeline outputs.

Computes Δ|partial corr(dep, CN | proximity)| using OLS residualization:
- dep_resid = dep - Proj(dep | bp_near, bp_far)
- dep_corr_resid = dep_corrected - Proj(dep_corrected | bp_near, bp_far)
- cn_resid = cn - Proj(cn | bp_near, bp_far)
- r0 = corr(dep_resid, cn_resid), r1 = corr(dep_corr_resid, cn_resid)
- per-gene metric = |r0| - |r1|

Bootstrap median (2,000 reps) and write summary CSV.
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LinearRegression

def compute_partial_corr_delta(base_folder: Path, min_rows: int = 20) -> dict:
    """
    Compute partial correlation delta from design matrix.
    
    Args:
        base_folder: Folder containing dependency_corrected.csv
        min_rows: Minimum rows per gene to include
        
    Returns:
        Dictionary with summary statistics
    """
    # Load dependency_corrected.csv (should have gene, cell_line, dependency, dependency_corrected, cn, bp_dist)
    dep_path = base_folder / "dependency_corrected.csv"
    if not dep_path.exists():
        raise FileNotFoundError(f"dependency_corrected.csv not found in {base_folder}")
    
    df = pd.read_csv(dep_path)
    
    # Check required columns
    required = ["gene", "cell_line", "dependency", "dependency_corrected", "cn"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    # Check for bp_dist or bp_within flags
    has_bp_dist = "bp_dist" in df.columns
    has_bp_near = "bp_within_100000" in df.columns or "bp_within_250000" in df.columns
    has_bp_far = "bp_within_1000000" in df.columns or "bp_within_2000000" in df.columns
    
    if not (has_bp_dist or (has_bp_near and has_bp_far)):
        raise ValueError("Need either bp_dist or bp_within_* columns for proximity")
    
    # Create bp_near and bp_far if needed
    if has_bp_dist:
        # Use 250k and 2Mb windows (adjustable, but match common defaults)
        w1, w2 = 250_000, 2_000_000
        df["bp_near"] = (df["bp_dist"] <= w1).astype(float)
        df["bp_far"] = ((df["bp_dist"] > w1) & (df["bp_dist"] <= w2)).astype(float)
    else:
        # Use existing flags
        if "bp_within_250000" in df.columns:
            df["bp_near"] = df["bp_within_250000"].astype(float)
        elif "bp_within_100000" in df.columns:
            df["bp_near"] = df["bp_within_100000"].astype(float)
        else:
            raise ValueError("Cannot determine bp_near")
        
        if "bp_within_2000000" in df.columns:
            df["bp_far"] = (df["bp_within_2000000"] & ~df["bp_near"]).astype(float)
        elif "bp_within_1000000" in df.columns:
            df["bp_far"] = (df["bp_within_1000000"] & ~df["bp_near"]).astype(float)
        else:
            df["bp_far"] = 0.0
    
    # Fill NaN proximity indicators
    df["bp_near"] = df["bp_near"].fillna(0.0)
    df["bp_far"] = df["bp_far"].fillna(0.0)
    
    # Per-gene computation
    results = []
    rng = np.random.default_rng(42)
    
    for gene, gdf in df.groupby("gene"):
        # Drop rows with missing critical data
        valid = gdf[["dependency", "dependency_corrected", "cn", "bp_near", "bp_far"]].dropna()
        if len(valid) < min_rows:
            continue
        
        try:
            # Extract features
            X_prox = valid[["bp_near", "bp_far"]].values
            y_dep = valid["dependency"].values
            y_dep_corr = valid["dependency_corrected"].values
            y_cn = valid["cn"].values
            
            # Residualize: remove projection onto proximity
            # dep_resid = dep - Proj(dep | bp_near, bp_far)
            lr_dep = LinearRegression()
            lr_dep.fit(X_prox, y_dep)
            dep_resid = y_dep - lr_dep.predict(X_prox)
            
            # dep_corr_resid = dep_corrected - Proj(dep_corrected | bp_near, bp_far)
            lr_dep_corr = LinearRegression()
            lr_dep_corr.fit(X_prox, y_dep_corr)
            dep_corr_resid = y_dep_corr - lr_dep_corr.predict(X_prox)
            
            # cn_resid = cn - Proj(cn | bp_near, bp_far)
            lr_cn = LinearRegression()
            lr_cn.fit(X_prox, y_cn)
            cn_resid = y_cn - lr_cn.predict(X_prox)
            
            # Compute partial correlations
            r0 = np.corrcoef(dep_resid, cn_resid)[0, 1]
            r1 = np.corrcoef(dep_corr_resid, cn_resid)[0, 1]
            
            if np.isnan(r0) or np.isnan(r1):
                continue
            
            delta = abs(r0) - abs(r1)
            results.append({
                "gene": gene,
                "partial_corr_delta_abs": delta,
                "r0_abs": abs(r0),
                "r1_abs": abs(r1),
                "n_obs": len(valid)
            })
        except Exception as e:
            continue
    
    if not results:
        return {
            "n_genes": 0,
            "median_partial_delta_abs_corr": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan
        }
    
    results_df = pd.DataFrame(results)
    vals = results_df["partial_corr_delta_abs"].values
    
    # Bootstrap median (2,000 reps)
    n = len(vals)
    boot_medians = []
    for _ in range(2000):
        idx = rng.integers(0, n, size=n)
        boot_medians.append(np.median(vals[idx]))
    
    boot_medians = np.array(boot_medians)
    
    summary = {
        "n_genes": int(len(results_df)),
        "median_partial_delta_abs_corr": float(np.median(vals)),
        "ci_low": float(np.percentile(boot_medians, 2.5)),
        "ci_high": float(np.percentile(boot_medians, 97.5)),
        "mean_partial_delta_abs_corr": float(np.mean(vals)),
        "frac_positive": float((vals > 0).mean())
    }
    
    # Save per-gene results
    results_df.to_csv(base_folder / "partial_corr_per_gene.csv", index=False)
    
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Compute partial correlation metric from pipeline outputs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("base_folder", type=str, help="Folder containing dependency_corrected.csv")
    parser.add_argument("--min-rows", type=int, default=20, help="Minimum rows per gene")
    
    args = parser.parse_args()
    base_folder = Path(args.base_folder)
    
    if not base_folder.exists():
        raise FileNotFoundError(f"Base folder not found: {base_folder}")
    
    summary = compute_partial_corr_delta(base_folder, min_rows=args.min_rows)
    
    # Save summary
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(base_folder / "partial_corr_summary.csv", index=False)
    
    # Print one-line summary
    print(f"{base_folder.name}: n={summary['n_genes']}, "
          f"median={summary['median_partial_delta_abs_corr']:.4f} "
          f"[{summary['ci_low']:.4f}, {summary['ci_high']:.4f}]")


if __name__ == "__main__":
    main()

