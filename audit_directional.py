#!/usr/bin/env python3
import os, re
import pandas as pd
import numpy as np

def _read_contrib_thresh(setdir: str, dm: pd.DataFrame) -> float:
    # default: 0.01 if standardized columns exist, else 0.1
    default = 0.01 if any(c.endswith("_std") for c in dm.columns) else 0.1
    p = os.path.join(setdir, "pilot_summary.txt")
    if not os.path.exists(p):
        return default
    with open(p) as f:
        for ln in f:
            if "contrib_thresh=" in ln:
                try:
                    return float(ln.strip().split("=")[1])
                except:
                    return default
    return default

def compute_directional_flags(setdir: str, genes_list, active_frac_thresh: float = 0.10):
    dm_path = os.path.join(setdir, "design_matrix.csv")
    coef_path = os.path.join(setdir, "models_coefficients.csv")
    dm = pd.read_csv(dm_path)
    coef = pd.read_csv(coef_path).rename(columns={
        "cn":"b_cn", "bp_within_100000":"b_w100k", "bp_within_1000000":"b_w1m", "inv_bp":"b_inv"
    })
    x_w100k = "bp_within_100000_std" if "bp_within_100000_std" in dm.columns else "bp_within_100000"
    x_w1m   = "bp_within_1000000_std" if "bp_within_1000000_std" in dm.columns else "bp_within_1000000"
    x_inv   = "inv_bp_std" if "inv_bp_std" in dm.columns else ("inv_bp" if "inv_bp" in dm.columns else None)
    df = dm.merge(coef[["gene","b_w100k","b_w1m","b_inv"]], on="gene", how="left")
    df["__x_w100k"] = df.get(x_w100k, 0.0).astype(float)
    df["__x_w1m"]   = df.get(x_w1m, 0.0).astype(float)
    df["__x_inv"]   = df.get(x_inv, 0.0).astype(float) if x_inv else 0.0
    df["prox_contrib"] = df["b_w100k"]*df["__x_w100k"] + df["b_w1m"]*df["__x_w1m"] + df["b_inv"]*df["__x_inv"]
    df["prox_only_corr"] = df["dependency"] - df["prox_contrib"]
    THR = _read_contrib_thresh(setdir, dm)
    ACTIVE_FRAC = active_frac_thresh
    rows=[]
    # lock to the provided 50 genes and do not drop any gene unless no usable rows
    gene_set = set(genes_list)
    for g, sub in df[df["gene"].isin(gene_set)].groupby("gene", sort=False):
        sub = sub.dropna(subset=["cn","dependency","prox_only_corr"])
        if len(sub) == 0:
            rows.append((g, False, 0.0, np.nan))
            continue
        active_frac = (sub["prox_contrib"].abs() >= THR).mean()
        r0 = sub[["cn","dependency"]].corr().iloc[0,1]
        r1 = sub[["cn","prox_only_corr"]].corr().iloc[0,1]
        drop = abs(r0) - abs(r1)
        dir_flag = (active_frac >= ACTIVE_FRAC) and (drop > 0)
        rows.append((g, dir_flag, active_frac, drop))
    out = pd.DataFrame(rows, columns=["gene","dir_flag","active_frac","drop"])
    # Ensure all genes in genes_list present
    missing = [g for g in genes_list if g not in set(out["gene"])]
    for g in missing:
        out = pd.concat([out, pd.DataFrame([{"gene":g,"dir_flag":False,"active_frac":0.0,"drop":np.nan}])], ignore_index=True)
    return out

