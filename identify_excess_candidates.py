#!/usr/bin/env python3
"""
Identify top candidate genes with excess signal (TRUE > shuffles).
"""

import os
import numpy as np
import pandas as pd

BASE = "out_v2"
RUNS = ["true", "shuffle_within_chrom", "shuffle_across_chrom"]


def load_coefs(tag):
    df = pd.read_csv(f"{BASE}/{tag}/model_coefficients.csv")
    df["abs_near"] = df["coef_bp_near"].abs()
    df["abs_far"] = df["coef_bp_far"].abs()
    df["abs_max"] = df[["abs_near", "abs_far"]].max(axis=1)
    return df[["gene", "abs_near", "abs_far", "abs_max"]]


def load_corrs(tag):
    p = f"{BASE}/{tag}/gene_level_correlations.csv"
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    return df[["gene", "delta_abs_corr"]]


coefs = {k: load_coefs(k) for k in RUNS}
corrs = {k: load_corrs(k) for k in RUNS}

# Merge TRUE vs shuffles, gene-aligned
m = coefs["true"].merge(coefs["shuffle_within_chrom"], on="gene", suffixes=("_true", "_shufW"))
m = m.merge(coefs["shuffle_across_chrom"], on="gene", suffixes=("", "_shufA"))

# Rename columns from third merge to have _shufA suffix
for col in ["abs_near", "abs_far", "abs_max"]:
    if col in m.columns and f"{col}_shufA" not in m.columns:
        m = m.rename(columns={col: f"{col}_shufA"})

# Bring Δ|corr|
dc = corrs["true"].merge(corrs["shuffle_within_chrom"], on="gene", suffixes=("_true", "_shufW"))
dc = dc.merge(corrs["shuffle_across_chrom"], on="gene", suffixes=("", "_shufA"))

# Rename delta_abs_corr column
if "delta_abs_corr" in dc.columns and "delta_abs_corr_shufA" not in dc.columns:
    dc = dc.rename(columns={"delta_abs_corr": "delta_abs_corr_shufA"})

m = m.merge(dc, on="gene", how="inner")

# Excess metrics per gene
m["near_excess_W"] = m["abs_near_true"] - m["abs_near_shufW"]
m["near_excess_A"] = m["abs_near_true"] - m["abs_near_shufA"]
m["max_excess_W"] = m["abs_max_true"] - m["abs_max_shufW"]
m["max_excess_A"] = m["abs_max_true"] - m["abs_max_shufA"]

m["deltacorr_excess_W"] = m["delta_abs_corr_true"] - m["delta_abs_corr_shufW"]
m["deltacorr_excess_A"] = m["delta_abs_corr_true"] - m["delta_abs_corr_shufA"]

# Primary interest flag: TRUE beats at least one shuffle on Δ|corr| AND has positive near/max excess
m["flag_excess_corr"] = (m["deltacorr_excess_W"] > 0) | (m["deltacorr_excess_A"] > 0)
m["flag_excess_coef"] = (m["near_excess_W"] > 0) | (m["near_excess_A"] > 0) | (m["max_excess_W"] > 0) | (m["max_excess_A"] > 0)
m["flag_interest"] = m["flag_excess_corr"] & m["flag_excess_coef"]

cands = m[m["flag_interest"]].copy()

# If very few, fall back: top by near-excess with positive Δ|corr| in TRUE (even if not exceeding shuffles)
if len(cands) < 50:
    fallback = m[(m["near_excess_W"] > 0) | (m["near_excess_A"] > 0)].copy()
    fallback = fallback[fallback["delta_abs_corr_true"] > 0]
    # score: combine rank of near-excess and delta_corr_true
    r1 = (-fallback[["near_excess_W", "near_excess_A"]].max(axis=1)).rank(method="average")
    r2 = (-fallback["delta_abs_corr_true"]).rank(method="average")
    fallback["score"] = r1 + r2
    fallback = fallback.sort_values("score")
    cands = pd.concat([cands, fallback.head(200)], ignore_index=True).drop_duplicates("gene")

# Rank candidates (prefer strong near_excess and Δ|corr| excess)
cands["score"] = (
    2.0 * cands[["near_excess_W", "near_excess_A"]].max(axis=1).rank(ascending=False)
    + 1.0 * cands[["deltacorr_excess_W", "deltacorr_excess_A"]].max(axis=1).rank(ascending=False)
    + 1.0 * cands["delta_abs_corr_true"].rank(ascending=False)
)
cands = cands.sort_values("score", ascending=False)

out_csv = "figs_v2/top_excess_candidates.csv"
cands_cols = [
    "gene", "abs_near_true", "abs_max_true",
    "near_excess_W", "near_excess_A", "max_excess_W", "max_excess_A",
    "delta_abs_corr_true", "deltacorr_excess_W", "deltacorr_excess_A", "score"
]
cands[cands_cols].to_csv(out_csv, index=False)

topN = min(150, len(cands))
with open("figs_v2/pilot_genes.txt", "w") as f:
    for g in cands.head(topN)["gene"]:
        f.write(str(g) + "\n")

print(f"Wrote {out_csv}  (n={len(cands)})")
print(f"Wrote figs_v2/pilot_genes.txt (n={topN})")

# Print a tiny preview
print("\nPreview:")
print(cands[cands_cols].head(12).to_string(index=False))

