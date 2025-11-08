#!/usr/bin/env python3
"""Evaluate TRUE vs ROTATE and select case-study genes."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


def prevalence_matched_directional(df_true: pd.DataFrame, df_rot: pd.DataFrame) -> dict:
    """Compute prevalence-matched directional metric."""
    # df_* must have columns: gene, coef_cn, coef_prox (or coef_prox_*), n
    # Tertiles by prevalence (n)
    out = {}
    
    for name, df in {"TRUE": df_true, "ROTATE": df_rot}.items():
        # Find proximity coefficient column
        prox_col = None
        for col in df.columns:
            if col.startswith("coef_prox") or col == "coef_prox":
                prox_col = col
                break
        
        if prox_col is None or prox_col not in df.columns:
            out[name] = {}
            continue
        
        d = df.dropna(subset=[prox_col, "n"]).copy()
        if len(d) == 0:
            out[name] = {}
            continue
        d["tertile"] = pd.qcut(d["n"], 3, labels=["low", "mid", "high"], duplicates="drop")
        frac = d.groupby("tertile")[prox_col].apply(lambda s: (s > 0).mean()).to_dict()
        out[name] = {k: float(v) for k, v in frac.items()}
    
    return out


def coef_effect_consistency(df_true: pd.DataFrame, dep: pd.DataFrame) -> dict:
    """Compute coef→effect consistency (Spearman correlation)."""
    # Expect per-gene residualized effect and coef_prox (or coef_prox_*)
    if dep.empty or "effect_residual" not in dep.columns:
        return {}
    
    # Find proximity coefficient column
    prox_col = None
    for col in df_true.columns:
        if col.startswith("coef_prox") or col == "coef_prox":
            prox_col = col
            break
    
    if prox_col is None:
        return {}
    
    merged = df_true.merge(dep, on="gene", how="inner")
    if len(merged) < 10:
        return {}
    
    rho, p = spearmanr(merged[prox_col], merged["effect_residual"], nan_policy="omit")
    return {"rho": float(rho) if pd.notna(rho) else np.nan, 
             "p": float(p) if pd.notna(p) else np.nan, 
             "n": int(len(merged))}


def apply_preregistered(df_true: pd.DataFrame, df_rot: pd.DataFrame, alpha_q: float) -> pd.DataFrame:
    """Apply pre-registered inclusion/exclusion criteria."""
    # Simplified placeholders (assumes inputs include needed metrics)
    t = df_true.copy()
    
    # Find proximity coefficient column (coef_prox or coef_prox_*)
    prox_col = None
    for col in t.columns:
        if col.startswith("coef_prox") or col == "coef_prox":
            prox_col = col
            break
    
    # 1) near-proximity coefficient percentile ≥95th (TRUE)
    if prox_col and prox_col in t.columns:
        thr = np.nanpercentile(np.abs(t[prox_col]), 95)
        t["pass_coef_pct"] = np.abs(t[prox_col]) >= thr
        t["coef_prox"] = t[prox_col]  # Standardize name
    else:
        t["pass_coef_pct"] = False
        t["coef_prox"] = 0.0
    
    # 2) Δ|corr(dep,CN)| > 0 after prox-only correction
    if "delta_abs_corr" in t.columns:
        t["pass_delta_corr"] = (t["delta_abs_corr"] > 0)
    else:
        t["pass_delta_corr"] = False
    
    # 3) coef→effect residual (prevalence-adjusted) in top quartile
    if "coef_effect_resid" in t.columns:
        q = np.nanpercentile(t["coef_effect_resid"], 75)
        t["pass_resid_q"] = t["coef_effect_resid"] >= q
    else:
        t["pass_resid_q"] = False
    
    # 4) robustness across kernels/windows (expect precomputed boolean)
    t["pass_robust"] = t.get("robust_across_kernels", True)
    
    # 5) ROTATE residual < TRUE residual
    if "coef_effect_resid" in t.columns and "coef_effect_resid" in df_rot.columns:
        rot = df_rot[["gene", "coef_effect_resid"]].rename(columns={"coef_effect_resid": "rot_resid"})
        t = t.merge(rot, on="gene", how="left")
        t["pass_true_gt_rotate"] = (t["coef_effect_resid"] > t["rot_resid"])
    else:
        t["pass_true_gt_rotate"] = False
    
    # Exclusions
    exclude_flags = []
    if "non_na_dep" in t.columns:
        exclude_flags.append(t["non_na_dep"] < 500)
    if "cn_var_rank" in t.columns:
        exclude_flags.append(t["cn_var_rank"] <= 5)
    if "depmap_qc_flag" in t.columns:
        exclude_flags.append(t["depmap_qc_flag"].fillna(False))
    
    if exclude_flags:
        t["exclude"] = pd.concat(exclude_flags, axis=1).any(axis=1)
    else:
        t["exclude"] = False
    
    # p→q
    if "pval_prox" in t.columns:
        pvals = t["pval_prox"].fillna(1.0).to_numpy()
        _, qvals, _, _ = multipletests(pvals, alpha=alpha_q, method="fdr_bh")
        t["qval_prox"] = qvals
    else:
        t["qval_prox"] = 1.0
    
    # Final keep
    keep = t[
        ~t["exclude"] & 
        t["pass_coef_pct"] & 
        t["pass_delta_corr"] & 
        t["pass_resid_q"] & 
        t["pass_robust"] & 
        t["pass_true_gt_rotate"] & 
        (t["qval_prox"] <= alpha_q)
    ].copy()
    
    return keep.sort_values("qval_prox")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--true-coefs", required=True)
    ap.add_argument("--rotate-coefs", required=True)
    ap.add_argument("--dep-residuals", required=False, default="")  # optional per-gene effect residuals
    ap.add_argument("--alpha-q", type=float, default=0.10)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    dfT = pd.read_csv(args.true_coefs)
    dfR = pd.read_csv(args.rotate_coefs)
    dep = pd.read_csv(args.dep_residuals) if args.dep_residuals and Path(args.dep_residuals).exists() else pd.DataFrame(columns=["gene", "effect_residual"])
    
    metrics = {
        "directional": prevalence_matched_directional(dfT, dfR),
        "coef_effect_consistency": coef_effect_consistency(dfT, dep) if not dep.empty else {},
    }
    (Path(args.out_dir) / "evaluation_metrics.json").write_text(json.dumps(metrics, indent=2))
    
    keep = apply_preregistered(dfT, dfR, alpha_q=args.alpha_q)
    keep.to_csv(Path(args.out_dir) / "case_study_genes.csv", index=False)
    
    print(f"[evaluate] Selected {len(keep)} case-study genes")
    print(f"[evaluate] Metrics saved to {Path(args.out_dir) / 'evaluation_metrics.json'}")


if __name__ == "__main__":
    main()

