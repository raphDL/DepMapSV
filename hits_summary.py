#!/usr/bin/env python3
"""
Rank genes by proximity-specific signal using design_matrix.csv and models_coefficients.csv.

Usage: python hits_summary.py out_pilot_final/unstable
"""

import sys, pandas as pd, numpy as np, os

SETDIR = sys.argv[1]  # e.g., out_pilot_final/unstable

dm   = pd.read_csv(os.path.join(SETDIR, "design_matrix.csv"))
coef = pd.read_csv(os.path.join(SETDIR, "models_coefficients.csv")).rename(
    columns={"cn":"b_cn","bp_within_100000":"b_w100k","bp_within_1000000":"b_w1m","inv_bp":"b_inv"}
)

df = dm.merge(coef[["gene","b_cn","b_w100k","b_w1m","b_inv"]], on="gene", how="left")

# Proximity contribution on the same scale the model used
df["prox_contrib"] = (df["b_w100k"]*df["bp_within_100000"].astype(float) +
                      df["b_w1m"]  *df["bp_within_1000000"].astype(float) +
                      df["b_inv"]  *df.get("inv_bp", 0.0).astype(float))

df["prox_only_corrected"] = df["dependency"] - df["prox_contrib"]

rows=[]
for g, sub in df.groupby("gene"):
    sub2 = sub[["cn","dependency","prox_only_corrected","prox_contrib"]].dropna()
    if len(sub2) < 12: 
        continue
    r0 = sub2[["cn","dependency"]].corr().iloc[0,1]
    r1 = sub2[["cn","prox_only_corrected"]].corr().iloc[0,1]
    drop = abs(r0) - abs(r1)
    active_frac = (sub2["prox_contrib"].abs() >= 0.01).mean()    # tune threshold if needed
    strength    = sub2.loc[sub2["prox_contrib"].abs()>=0.01, "prox_contrib"].abs().median() if (sub2["prox_contrib"].abs()>=0.01).any() else 0.0
    cn_var      = sub2["cn"].var()
    rows.append({
        "gene": g,
        "n_cells": len(sub2),
        "prox_active_frac": active_frac,
        "prox_only_abs_corr_drop": drop,
        "prox_contrib_median_if_active": strength,
        "coef_abs_w100k": sub["b_w100k"].abs().median(),
        "coef_abs_w1m":   sub["b_w1m"].abs().median(),
        "coef_abs_inv":   sub["b_inv"].abs().median(),
        "cn_var": cn_var,
    })

out = pd.DataFrame(rows).sort_values(
    by=["prox_active_frac","prox_only_abs_corr_drop","prox_contrib_median_if_active"],
    ascending=[False, False, False]
)

out.to_csv(os.path.join(SETDIR, "pilot_hits_summary.csv"), index=False)
print("Top 10 hits by proximity signal:")
print(out.head(10).to_string(index=False))
print(f"\nSaved to: {os.path.join(SETDIR, 'pilot_hits_summary.csv')}")

