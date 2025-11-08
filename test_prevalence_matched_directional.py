#!/usr/bin/env python3
"""
Prevalence-matched directional metric (fixes the "more breakpoints → more chances" bias).
"""

import pandas as pd
import numpy as np

def load_set(base):
    dm = pd.read_csv(f"{base}/design_matrix.csv")
    co = pd.read_csv(f"{base}/models_coefficients.csv")
    # Merge coefficients (they're named bp_within_100000, bp_within_1000000, inv_bp)
    df = dm.merge(co[["gene", "bp_within_100000", "bp_within_1000000", "inv_bp"]],
                  on="gene", how="left", suffixes=("", "_coef"))
    w1, w2 = 250000, 2000000
    df["bp_near"] = (df["bp_dist"].fillna(5_000_000) <= w1).astype(float)
    df["bp_far"] = ((df["bp_dist"].fillna(5_000_000) > w1) & (df["bp_dist"] <= w2)).astype(float)
    # Use coefficient columns from models_coefficients (bp_within_100000, bp_within_1000000, inv_bp)
    df["prox_contrib"] = (df["bp_within_100000_coef"].fillna(0) * df["bp_within_100000"] + 
                          df["bp_within_1000000_coef"].fillna(0) * df["bp_within_1000000"] +
                          df["inv_bp_coef"].fillna(0) * df.get("inv_bp", 0))
    return df

BASE_TRUE = "out_focus_true/unstable"
BASE_ROTATE = "out_focus_rotate/unstable"

t = load_set(BASE_TRUE)
r = load_set(BASE_ROTATE)

def per_gene_flags(df, contrib_thresh=0.01, cell_frac=0.10):
    out = []
    for g, sub in df.groupby("gene"):
        active_frac = (sub["prox_contrib"].abs() >= contrib_thresh).mean()
        a = sub[["dependency", "cn"]].dropna()
        b = sub[["dependency_corrected", "cn"]].dropna()
        if len(a) > 20 and len(b) > 20:
            r0 = a.corr().iloc[0, 1]
            r1 = b.corr().iloc[0, 1]
            if pd.notna(r0) and pd.notna(r1):
                drop = abs(r0) - abs(r1)
                out.append((g, active_frac, drop))
    return pd.DataFrame(out, columns=["gene", "active_frac", "corr_drop"])

pt = per_gene_flags(t)
pr = per_gene_flags(r)

# prevalence by gene (using <=1Mb)
prev_t = t.assign(p1m=(t["bp_dist"] <= 1_000_000)).groupby("gene")["p1m"].mean().rename("prev")
prev_r = r.assign(p1m=(r["bp_dist"] <= 1_000_000)).groupby("gene")["p1m"].mean().rename("prev")

gt = pt.merge(prev_t, left_on="gene", right_index=True)
gr = pr.merge(prev_r, left_on="gene", right_index=True)

# Create deciles on TRUE prevalence, cut ROTATE to same bin edges (to compare like with like)
bins = np.quantile(gt["prev"].values, np.linspace(0, 1, 11))
bins[0] -= 1e-9
bins[-1] += 1e-9
gt["bin"] = pd.cut(gt["prev"], bins=bins, labels=False)
gr["bin"] = pd.cut(gr["prev"], bins=bins, labels=False)

# Directional flag: active AND corr_drop>0
gt["dir"] = (gt["active_frac"] >= 0.10) & (gt["corr_drop"] > 0)
gr["dir"] = (gr["active_frac"] >= 0.10) & (gr["corr_drop"] > 0)

# Bin means
bt = gt.groupby("bin")["dir"].mean()
br = gr.groupby("bin")["dir"].mean()
common_bins = sorted(set(bt.index).intersection(br.index))

# Equal-weight average across bins
dir_true_eq = float(np.mean([bt.loc[b] for b in common_bins]))
dir_rot_eq = float(np.mean([br.loc[b] for b in common_bins]))
excess_eq = dir_true_eq - dir_rot_eq

# Bootstrap CI (resample bins)
B = 10000
rng = np.random.default_rng(0)
boots = []
arr_t = np.array([bt.loc[b] for b in common_bins], float)
arr_r = np.array([br.loc[b] for b in common_bins], float)
n = len(common_bins)
for _ in range(B):
    idx = rng.integers(0, n, size=n)
    boots.append(np.mean(arr_t[idx]) - np.mean(arr_r[idx]))
ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])

print("Prevalence-matched directional:")
print(f"TRUE={dir_true_eq:.3f}  ROTATE={dir_rot_eq:.3f}  EXCESS={excess_eq:+.3f}  CI[{ci_lo:+.3f},{ci_hi:+.3f}]")

