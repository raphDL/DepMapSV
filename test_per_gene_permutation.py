#!/usr/bin/env python3
"""
Per-gene permutation p-values for the directional effect (FDR across genes).
"""

import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

def load(base):
    dm = pd.read_csv(f"{base}/design_matrix.csv")
    co = pd.read_csv(f"{base}/models_coefficients.csv")
    w1, w2 = 250000, 2000000
    if "bp_dist" in dm.columns:
        dm["bp_near"] = (dm["bp_dist"].fillna(5_000_000) <= w1).astype(int)
        dm["bp_far"] = ((dm["bp_dist"].fillna(5_000_000) > w1) & (dm["bp_dist"] <= w2)).astype(int)
    else:
        # Use bp_within flags if bp_dist not available
        dm["bp_near"] = dm.get("bp_within_100000", 0).astype(int)
        dm["bp_far"] = (dm.get("bp_within_1000000", 0) & ~dm["bp_near"]).astype(int)
    # Merge coefficients (bp_within_100000, bp_within_1000000, inv_bp)
    coef_cols = ["gene", "bp_within_100000", "bp_within_1000000", "inv_bp"]
    coef_cols = [c for c in coef_cols if c in co.columns]
    dm = dm.merge(co[coef_cols], on="gene", how="left", suffixes=("", "_coef"))
    
    # Find cn column
    cn_col = None
    for col in ["cn", "cn_x"]:
        if col in dm.columns:
            cn_col = col
            break
    
    # Compute proximity contribution - use _coef suffixed versions if available
    prox_contrib = 0
    if "bp_within_100000_coef" in dm.columns:
        prox_contrib += dm["bp_within_100000_coef"].fillna(0) * dm["bp_within_100000"]
    elif "bp_within_100000" in dm.columns:
        prox_contrib += dm.get("bp_within_100000", 0) * dm["bp_within_100000"]
    
    if "bp_within_1000000_coef" in dm.columns:
        prox_contrib += dm["bp_within_1000000_coef"].fillna(0) * dm["bp_within_1000000"]
    elif "bp_within_1000000" in dm.columns:
        prox_contrib += dm.get("bp_within_1000000", 0) * dm["bp_within_1000000"]
    
    if "inv_bp_coef" in dm.columns:
        prox_contrib += dm["inv_bp_coef"].fillna(0) * dm.get("inv_bp", 0)
    elif "inv_bp" in dm.columns:
        prox_contrib += dm.get("inv_bp", 0) * dm.get("inv_bp", 0)
    
    dm["prox_contrib"] = prox_contrib
    dm["prox_only_corrected"] = dm["dependency"] - dm["prox_contrib"]
    
    # Rename cn column for consistency
    if cn_col and cn_col != "cn":
        dm = dm.rename(columns={cn_col: "cn"})
    
    return dm

def gene_pvals(dm, B=400, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for g, sub in dm.groupby("gene"):
        if "cn" not in sub.columns or "dependency" not in sub.columns:
            continue
        a = sub[["dependency", "cn"]].dropna()
        b = sub[["prox_only_corrected", "cn"]].dropna()
        if len(a) > 40 and len(b) > 40:
            r0 = a.corr().iloc[0, 1]
            r1 = b.corr().iloc[0, 1]
            obs = abs(r0) - abs(r1)
            if not np.isfinite(obs):
                continue
            
            # shuffle proximity assignment across cells (preserve counts)
            # Use bp_within_100000 (more variation than 1Mb which is often all cells)
            if "bp_within_100000" in sub.columns:
                k = int(sub["bp_within_100000"].sum())
            elif "bp_within_1000000" in sub.columns:
                k = int(sub["bp_within_1000000"].sum())
            else:
                k = int(sub.get("bp_near", 0).sum() + sub.get("bp_far", 0).sum())
            n = len(sub)
            if k == 0 or k >= n:
                continue
            
            vals = []
            cn = sub["cn"].values
            dep = sub["dependency"].values
            for _ in range(B):
                mask = np.zeros(n, dtype=int)
                mask[rng.choice(n, size=int(k), replace=False)] = 1
                # recompute prox contrib from mask with same coef signs
                # Use bp_within_100000 coefficient as proxy for proximity effect
                pc = mask * (sub["bp_within_100000_coef"].fillna(0).values)  # treat all mask as "proximity"
                dep_corr = dep - pc
                # corr drop under shuffle
                try:
                    r0s = np.corrcoef(dep, cn)[0, 1]
                    r1s = np.corrcoef(dep_corr, cn)[0, 1]
                    vals.append(abs(r0s) - abs(r1s))
                except:
                    pass
            
            if len(vals) > 50:
                p = (np.sum(np.array(vals) >= obs) + 1) / (len(vals) + 1)
                rows.append((g, obs, p))
    
    out = pd.DataFrame(rows, columns=["gene", "corr_drop", "p"])
    if len(out):
        out["q"] = multipletests(out["p"].values, method="fdr_bh")[1]
    return out

BASE_TRUE = "out_focus_true/unstable"
BASE_ROT = "out_focus_rotate/unstable"

dt = load(BASE_TRUE)
dr = load(BASE_ROT)

pt = gene_pvals(dt)
pr = gene_pvals(dr)

qt = (pt["q"] < 0.1).mean() if len(pt) else 0
qr = (pr["q"] < 0.1).mean() if len(pr) else 0

print(f"Directional prox-only (perm p): TRUE q<0.1 frac={qt:.3f}  n={len(pt)}")
print(f"Directional prox-only (perm p): ROTATE q<0.1 frac={qr:.3f}  n={len(pr)}")

