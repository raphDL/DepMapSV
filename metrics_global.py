#!/usr/bin/env python3
"""
Global Metrics - Prevalence-matched directional, coef→effect consistency, excess signal.

Pre-registered metrics:
- Directional prevalence-matched (deciles/tertiles): TRUE vs ROTATE, equal-weight average, 10k bootstrap CI
- Coef→effect consistency: Spearman ρ between mean |prox_contrib| and Δ|corr(dep,CN)| after prox-only correction
- Excess signal: coefficient-excess and Δ|corr|-excess (TRUE − control) with CIs
"""
import argparse
import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import HuberRegressor

# Pre-registered proximity kernels
KERNEL_COLS = [
    "prox_exp_50k", "prox_exp_100k", "prox_exp_250k", "prox_exp_500k",
    "prox_inv", "prox_spline"
]


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


def load_data(base_dir: str, kernel_col: str) -> tuple:
    """
    Load design matrix, coefficients, and corrected dependencies.
    
    Returns:
        (design_df, coef_df, corrected_df)
    """
    design_path = os.path.join(base_dir, "design_matrix_wgs.csv")
    coef_path = os.path.join(base_dir, "model_coefficients.csv")
    corr_path = os.path.join(base_dir, "dependency_corrected.csv")
    
    design = pd.read_csv(design_path)
    coef = pd.read_csv(coef_path)
    corrected = pd.read_csv(corr_path)
    
    # Merge to get proximity contribution
    # Need to compute prox_contrib from design and coefficients
    coef_kernel_col = f"coef_{kernel_col}"
    if coef_kernel_col not in coef.columns:
        raise ValueError(f"Kernel coefficient column '{coef_kernel_col}' not found in coefficients")
    
    # Merge coefficients
    design = design.merge(coef[["gene", coef_kernel_col]], on="gene", how="left")
    design["prox_contrib"] = design[coef_kernel_col].fillna(0.0) * design[kernel_col].fillna(0.0)
    
    # Merge corrected dependencies
    design = design.merge(corrected[["gene", "cell_line", "dependency", "prox_only_corrected"]], 
                          on=["gene", "cell_line"], how="left")
    
    return design, coef, corrected


def compute_per_gene_metrics(df: pd.DataFrame, contrib_thresh: float = 0.01, 
                              cell_frac: float = 0.10) -> pd.DataFrame:
    """
    Compute per-gene metrics: active_frac, corr_drop.
    
    active_frac: fraction of cells with |prox_contrib| >= contrib_thresh
    corr_drop: Δ|corr(dep,CN)| = |corr(dep,CN)| - |corr(prox_only_corrected,CN)|
    """
    rows = []
    
    for gene, sub in df.groupby("gene"):
        # Active fraction
        active_frac = (sub["prox_contrib"].abs() >= contrib_thresh).mean()
        
        # Correlation drop
        a = sub[["dependency", "cn"]].dropna()
        b = sub[["prox_only_corrected", "cn"]].dropna()
        
        if len(a) > 20 and len(b) > 20:
            try:
                r0 = a.corr().iloc[0, 1]
                r1 = b.corr().iloc[0, 1]
                if pd.notna(r0) and pd.notna(r1):
                    corr_drop = abs(r0) - abs(r1)
                    rows.append({
                        "gene": gene,
                        "active_frac": float(active_frac),
                        "corr_drop": float(corr_drop),
                        "mean_prox_contrib_abs": float(sub["prox_contrib"].abs().mean()),
                        "n": len(sub)
                    })
            except Exception:
                continue
    
    return pd.DataFrame(rows)


def compute_prevalence_matched_directional(true_df: pd.DataFrame, rotate_df: pd.DataFrame,
                                          n_bins: int = 10, n_bootstrap: int = 10000) -> dict:
    """
    Compute prevalence-matched directional metric.
    
    Stratify genes by proximity prevalence (deciles or tertiles), compare TRUE vs ROTATE
    within each stratum, equal-weight average, bootstrap CI.
    """
    logging.info(f"Computing prevalence-matched directional (n_bins={n_bins})...")
    
    # Compute prevalence (fraction of cells with breakpoint within 1Mb)
    prev_true = true_df.assign(
        p1m=(true_df["bp_dist"].fillna(5_000_000) <= 1_000_000)
    ).groupby("gene")["p1m"].mean().rename("prev")
    
    prev_rotate = rotate_df.assign(
        p1m=(rotate_df["bp_dist"].fillna(5_000_000) <= 1_000_000)
    ).groupby("gene")["p1m"].mean().rename("prev")
    
    # Per-gene metrics
    true_metrics = compute_per_gene_metrics(true_df)
    rotate_metrics = compute_per_gene_metrics(rotate_df)
    
    # Merge prevalence
    true_metrics = true_metrics.merge(prev_true, left_on="gene", right_index=True, how="left")
    rotate_metrics = rotate_metrics.merge(prev_rotate, left_on="gene", right_index=True, how="left")
    
    # Create bins on TRUE prevalence
    valid_prev = true_metrics["prev"].dropna()
    if len(valid_prev) < n_bins:
        logging.warning(f"Not enough genes for {n_bins} bins, using {len(valid_prev)} bins")
        n_bins = len(valid_prev)
    
    bins = np.quantile(valid_prev.values, np.linspace(0, 1, n_bins + 1))
    bins[0] -= 1e-9
    bins[-1] += 1e-9
    
    true_metrics["bin"] = pd.cut(true_metrics["prev"], bins=bins, labels=False)
    rotate_metrics["bin"] = pd.cut(rotate_metrics["prev"], bins=bins, labels=False)
    
    # Directional flag: active AND corr_drop > 0
    true_metrics["dir"] = (true_metrics["active_frac"] >= 0.10) & (true_metrics["corr_drop"] > 0)
    rotate_metrics["dir"] = (rotate_metrics["active_frac"] >= 0.10) & (rotate_metrics["corr_drop"] > 0)
    
    # Bin means
    bt = true_metrics.groupby("bin")["dir"].mean()
    br = rotate_metrics.groupby("bin")["dir"].mean()
    common_bins = sorted(set(bt.index).intersection(br.index))
    
    if len(common_bins) == 0:
        logging.warning("No common bins between TRUE and ROTATE")
        return {"error": "No common bins"}
    
    # Per-bin table
    bin_table = pd.DataFrame({
        "bin": common_bins,
        "true_dir_frac": [bt.loc[b] for b in common_bins],
        "rotate_dir_frac": [br.loc[b] for b in common_bins],
        "excess": [bt.loc[b] - br.loc[b] for b in common_bins],
        "n_genes_true": [true_metrics[true_metrics["bin"] == b].shape[0] for b in common_bins],
        "n_genes_rotate": [rotate_metrics[rotate_metrics["bin"] == b].shape[0] for b in common_bins]
    })
    
    # Equal-weight average
    dir_true_eq = float(np.mean([bt.loc[b] for b in common_bins]))
    dir_rot_eq = float(np.mean([br.loc[b] for b in common_bins]))
    excess_eq = dir_true_eq - dir_rot_eq
    
    # Bootstrap CI (resample bins)
    rng = np.random.default_rng(0)
    arr_t = np.array([bt.loc[b] for b in common_bins], dtype=float)
    arr_r = np.array([br.loc[b] for b in common_bins], dtype=float)
    n = len(common_bins)
    
    boots = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        boots.append(np.mean(arr_t[idx]) - np.mean(arr_r[idx]))
    
    ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])
    
    return {
        "dir_true": dir_true_eq,
        "dir_rotate": dir_rot_eq,
        "excess": excess_eq,
        "ci_lower": float(ci_lo),
        "ci_upper": float(ci_hi),
        "n_bins": len(common_bins),
        "n_bootstrap": n_bootstrap,
        "bin_table": bin_table.to_dict("records")
    }


def compute_coef_effect_consistency(true_df: pd.DataFrame, rotate_df: pd.DataFrame) -> dict:
    """
    Compute coef→effect consistency.
    
    Spearman ρ between mean |prox_contrib| and Δ|corr(dep,CN)| after prox-only correction.
    Also compute robust slope adjusted for prevalence and n.
    """
    logging.info("Computing coef→effect consistency...")
    
    true_metrics = compute_per_gene_metrics(true_df)
    rotate_metrics = compute_per_gene_metrics(rotate_df)
    
    # Filter to genes with positive corr_drop (directional)
    true_dir = true_metrics[true_metrics["corr_drop"] > 0].copy()
    rotate_dir = rotate_metrics[rotate_metrics["corr_drop"] > 0].copy()
    
    results = {}
    
    # All genes
    for label, metrics in [("true_all", true_metrics), ("rotate_all", rotate_metrics)]:
        if len(metrics) < 10:
            results[label] = {"rho": np.nan, "p": np.nan, "n": len(metrics)}
            continue
        
        rho, p = spearmanr(metrics["mean_prox_contrib_abs"], metrics["corr_drop"])
        results[label] = {
            "rho": float(rho),
            "p": float(p),
            "n": int(len(metrics))
        }
    
    # Directional genes only (Δ>0)
    for label, metrics in [("true_directional", true_dir), ("rotate_directional", rotate_dir)]:
        if len(metrics) < 10:
            results[label] = {"rho": np.nan, "p": np.nan, "n": len(metrics)}
            continue
        
        rho, p = spearmanr(metrics["mean_prox_contrib_abs"], metrics["corr_drop"])
        results[label] = {
            "rho": float(rho),
            "p": float(p),
            "n": int(len(metrics))
        }
    
    # Robust slope (Huber regression) for TRUE directional genes
    if len(true_dir) >= 20:
        try:
            X = true_dir[["mean_prox_contrib_abs", "n"]].values
            y = true_dir["corr_drop"].values
            
            model = HuberRegressor(epsilon=1.35, alpha=1e-4)
            model.fit(X, y)
            
            results["true_directional_robust_slope"] = {
                "slope_prox": float(model.coef_[0]),
                "slope_n": float(model.coef_[1]),
                "intercept": float(model.intercept_)
            }
        except Exception as e:
            logging.warning(f"Robust slope fitting failed: {e}")
            results["true_directional_robust_slope"] = None
    
    return results


def compute_excess_signal(true_coef: pd.DataFrame, rotate_coef: pd.DataFrame,
                         true_metrics: pd.DataFrame, rotate_metrics: pd.DataFrame) -> dict:
    """
    Compute excess signal: coefficient-excess and Δ|corr|-excess (TRUE − ROTATE).
    """
    logging.info("Computing excess signal...")
    
    # Merge on gene
    merged = true_coef.merge(rotate_coef, on="gene", suffixes=("_true", "_rotate"), how="inner")
    
    # Find kernel coefficient columns
    kernel_coef_cols = [c for c in merged.columns if c.startswith("coef_prox_")]
    if not kernel_coef_cols:
        # Try alternative naming
        kernel_coef_cols = [c for c in merged.columns if "coef" in c and "prox" in c]
    
    excess_coef = {}
    for col in kernel_coef_cols:
        if col.endswith("_true") and col.replace("_true", "_rotate") in merged.columns:
            col_rot = col.replace("_true", "_rotate")
            excess = merged[col] - merged[col_rot]
            excess_coef[col.replace("_true", "")] = {
                "mean": float(excess.mean()),
                "median": float(excess.median()),
                "q25": float(excess.quantile(0.25)),
                "q75": float(excess.quantile(0.75))
            }
    
    # Δ|corr| excess
    merged_metrics = true_metrics.merge(rotate_metrics, on="gene", suffixes=("_true", "_rotate"), how="inner")
    if "corr_drop_true" in merged_metrics.columns and "corr_drop_rotate" in merged_metrics.columns:
        corr_excess = merged_metrics["corr_drop_true"] - merged_metrics["corr_drop_rotate"]
        excess_corr = {
            "mean": float(corr_excess.mean()),
            "median": float(corr_excess.median()),
            "q25": float(corr_excess.quantile(0.25)),
            "q75": float(corr_excess.quantile(0.75)),
            "n": int(len(corr_excess))
        }
    else:
        excess_corr = None
    
    return {
        "coefficient_excess": excess_coef,
        "corr_drop_excess": excess_corr
    }


def main():
    ap = argparse.ArgumentParser(
        description="Compute global metrics: prevalence-matched directional, coef→effect, excess"
    )
    ap.add_argument("--true", required=True,
                    help="TRUE output directory")
    ap.add_argument("--rotate", required=True,
                    help="ROTATE output directory")
    ap.add_argument("--shuffle", default=None,
                    help="SHUFFLE output directory (optional)")
    ap.add_argument("--out", required=True,
                    help="Output directory")
    ap.add_argument("--kernel", default="prox_exp_100k",
                    help="Kernel column name (default: prox_exp_100k)")
    ap.add_argument("--n-bins", type=int, default=10,
                    help="Number of prevalence bins (default: 10 for deciles)")
    ap.add_argument("--n-bootstrap", type=int, default=10000,
                    help="Bootstrap iterations (default: 10000)")
    ap.add_argument("-v", "--verbose", action="count", default=1,
                    help="Increase verbosity")
    args = ap.parse_args()
    
    setup_logging(args.verbose)
    
    os.makedirs(args.out, exist_ok=True)
    
    # Load data
    logging.info("Loading TRUE data...")
    true_design, true_coef, true_corrected = load_data(args.true, args.kernel)
    
    logging.info("Loading ROTATE data...")
    rotate_design, rotate_coef, rotate_corrected = load_data(args.rotate, args.kernel)
    
    # Prevalence-matched directional
    logging.info("Computing prevalence-matched directional metric...")
    dir_metrics = compute_prevalence_matched_directional(
        true_design, rotate_design, n_bins=args.n_bins, n_bootstrap=args.n_bootstrap
    )
    
    dir_table = pd.DataFrame(dir_metrics.pop("bin_table"))
    dir_table_path = os.path.join(args.out, "prevalence_matched_directional.csv")
    dir_table.to_csv(dir_table_path, index=False)
    logging.info(f"Saved per-bin table to {dir_table_path}")
    
    # Coef→effect consistency
    logging.info("Computing coef→effect consistency...")
    true_metrics = compute_per_gene_metrics(true_design)
    rotate_metrics = compute_per_gene_metrics(rotate_design)
    
    coef_effect = compute_coef_effect_consistency(true_design, rotate_design)
    
    coef_effect_path = os.path.join(args.out, "coef_effect_summary.csv")
    coef_effect_df = pd.DataFrame([
        {"metric": k, **v} for k, v in coef_effect.items() if isinstance(v, dict)
    ])
    coef_effect_df.to_csv(coef_effect_path, index=False)
    logging.info(f"Saved coef→effect summary to {coef_effect_path}")
    
    # Excess signal
    logging.info("Computing excess signal...")
    excess = compute_excess_signal(true_coef, rotate_coef, true_metrics, rotate_metrics)
    
    # Save all metrics
    all_metrics = {
        "prevalence_matched_directional": dir_metrics,
        "coef_effect_consistency": coef_effect,
        "excess_signal": excess
    }
    
    metrics_path = os.path.join(args.out, "excess_signal.json")
    with open(metrics_path, "w") as f:
        json.dump(all_metrics, f, indent=2)
    logging.info(f"Saved all metrics to {metrics_path}")
    
    # Summary report
    print("\n=== Global Metrics Summary ===")
    print(f"\nPrevalence-matched directional (n_bins={args.n_bins}):")
    print(f"  TRUE:   {dir_metrics['dir_true']:.4f}")
    print(f"  ROTATE: {dir_metrics['dir_rotate']:.4f}")
    print(f"  EXCESS: {dir_metrics['excess']:+.4f} [{dir_metrics['ci_lower']:+.4f}, {dir_metrics['ci_upper']:+.4f}]")
    
    print(f"\nCoef→effect consistency (Spearman ρ):")
    for key, val in coef_effect.items():
        if isinstance(val, dict) and "rho" in val:
            print(f"  {key}: ρ={val['rho']:.4f}, p={val['p']:.4e}, n={val['n']}")


if __name__ == "__main__":
    main()

