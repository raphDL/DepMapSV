#!/usr/bin/env python3
"""
Coefficient→effect consistency: are big proximity coefs actually fixing CN coupling?
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr

def load(base):
    dm = pd.read_csv(f"{base}/design_matrix.csv")
    co = pd.read_csv(f"{base}/models_coefficients.csv")
    # Only merge coefficient columns, avoid conflicts
    coef_cols = ["gene", "bp_within_100000", "bp_within_1000000", "inv_bp"]
    coef_cols = [c for c in coef_cols if c in co.columns]
    df = dm.merge(co[coef_cols], on="gene", how="left", suffixes=("", "_coef"))
    
    # Find cn column (might be cn, cn_x, etc.)
    cn_col = None
    for col in ["cn", "cn_x"]:
        if col in df.columns:
            cn_col = col
            break
    
    if cn_col is None:
        return pd.DataFrame(columns=["gene", "abs_coef_max", "corr_drop"])
    
    info = []
    for g, sub in df.groupby("gene"):
        a = sub[["dependency", cn_col]].dropna()
        b = sub[["dependency_corrected", cn_col]].dropna()
        if len(a) > 20 and len(b) > 20:
            try:
                r0 = a.corr().iloc[0, 1]
                r1 = b.corr().iloc[0, 1]
                if pd.notna(r0) and pd.notna(r1):
                    drop = abs(r0) - abs(r1)
                    # Get max abs coefficient - check both original and _coef suffixed versions
                    coef_vals = []
                    for col_base in ["bp_within_100000", "bp_within_1000000", "inv_bp"]:
                        for col in [col_base, f"{col_base}_coef"]:
                            if col in sub.columns:
                                val = sub[col].iloc[0]
                                if pd.notna(val):
                                    coef_vals.append(abs(val))
                    if coef_vals:
                        coef = np.nanmax(coef_vals)
                        info.append((g, coef, drop))
            except:
                continue
    return pd.DataFrame(info, columns=["gene", "abs_coef_max", "corr_drop"])

BASE_TRUE = "out_focus_true/unstable"
BASE_ROT = "out_focus_rotate/unstable"

dt = load(BASE_TRUE)
dr = load(BASE_ROT)

rt = spearmanr(dt["abs_coef_max"], dt["corr_drop"], nan_policy="omit")
rr = spearmanr(dr["abs_coef_max"], dr["corr_drop"], nan_policy="omit")

print("Spearman(abs_coef_max, corr_drop)")
print(f"TRUE:   rho={rt.correlation:.3f}  p={rt.pvalue:.2e}  n={len(dt)}")
print(f"ROTATE: rho={rr.correlation:.3f}  p={rr.pvalue:.2e}  n={len(dr)}")

