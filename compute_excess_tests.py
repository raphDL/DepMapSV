#!/usr/bin/env python3
"""
Compute excess signal tests comparing TRUE vs shuffles.
"""

import json
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = "out_v2"
RUNS = ["true", "shuffle_within_chrom", "shuffle_across_chrom"]
rng = np.random.default_rng(123)


def boot_ci(vals, fn=np.median, B=2000):
    vals = np.asarray(vals)
    n = len(vals)
    boots = [fn(vals[rng.integers(0, n, n)]) for _ in range(B)]
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def load_coefs(tag):
    df = pd.read_csv(f"{BASE}/{tag}/model_coefficients.csv")
    # expected cols: gene, coef_cn, coef_bp_near, coef_bp_far, r2, n_obs, rmse, intercept
    df["abs_near"] = df["coef_bp_near"].abs()
    df["abs_far"] = df["coef_bp_far"].abs()
    df["abs_max"] = df[["abs_near", "abs_far"]].max(axis=1)
    return df[["gene", "abs_near", "abs_far", "abs_max"]]


def load_corrs(tag):
    p = f"{BASE}/{tag}/gene_level_correlations.csv"
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    # expected cols: gene, r_before, r_after, delta_abs_corr, n_obs
    return df[["gene", "delta_abs_corr"]]


coefs = {k: load_coefs(k) for k in RUNS}
corrs = {k: load_corrs(k) for k in RUNS}

# 1) Coefficient-excess: TRUE minus each shuffle (gene-aligned)
rows = []
for metric in ["abs_near", "abs_far", "abs_max"]:
    for ctrl in ["shuffle_within_chrom", "shuffle_across_chrom"]:
        merged = coefs["true"].merge(coefs[ctrl], on="gene", suffixes=("_true", "_shuf"))
        merged["diff"] = merged[f"{metric}_true"] - merged[f"{metric}_shuf"]
        med = np.median(merged["diff"])
        lo, hi = boot_ci(merged["diff"], fn=np.median, B=2000)
        rows.append(dict(
            test="coef_excess",
            metric=metric,
            versus=ctrl,
            n_genes=len(merged),
            median_excess=med,
            ci_low=lo,
            ci_high=hi
        ))

coef_excess = pd.DataFrame(rows)

# 2) Δ|corr|-excess: TRUE minus each shuffle (gene-aligned)
rows = []
for ctrl in ["shuffle_within_chrom", "shuffle_across_chrom"]:
    if corrs["true"] is None or corrs[ctrl] is None:
        continue
    merged = corrs["true"].merge(corrs[ctrl], on="gene", suffixes=("_true", "_shuf"))
    merged["diff"] = merged["delta_abs_corr_true"] - merged["delta_abs_corr_shuf"]
    med = np.median(merged["diff"])
    lo, hi = boot_ci(merged["diff"], fn=np.median, B=2000)
    rows.append(dict(
        test="delta_corr_excess",
        metric="delta_abs_corr",
        versus=ctrl,
        n_genes=len(merged),
        median_excess=med,
        ci_low=lo,
        ci_high=hi
    ))

corr_excess = pd.DataFrame(rows)

# 3) Simple table of medians per run for reference
summary = []
for k in RUNS:
    # coef medians
    m = coefs[k][["abs_near", "abs_far", "abs_max"]].median().to_dict()
    # delta_abs_corr median
    dc = np.nan
    if corrs[k] is not None and len(corrs[k]):
        dc = float(corrs[k]["delta_abs_corr"].median())
    summary.append(dict(
        run=k,
        med_abs_near=float(m["abs_near"]),
        med_abs_far=float(m["abs_far"]),
        med_abs_max=float(m["abs_max"]),
        med_delta_abs_corr=dc
    ))

summary = pd.DataFrame(summary)

# Save CSVs
os.makedirs("figs_v2", exist_ok=True)
coef_excess.to_csv("figs_v2/coef_excess.csv", index=False)
corr_excess.to_csv("figs_v2/delta_corr_excess.csv", index=False)
summary.to_csv("figs_v2/run_medians.csv", index=False)

# Quick plots
plt.figure()
for metric in ["abs_near", "abs_far", "abs_max"]:
    vals = [summary.loc[summary["run"] == k, f"med_{metric}"].values[0] for k in RUNS]
    plt.plot(RUNS, vals, marker='o', label=metric)
plt.xticks(rotation=20)
plt.ylabel("Median |coef|")
plt.title("Proximity coefficient medians")
plt.legend()
plt.tight_layout()
plt.savefig("figs_v2/coef_medians.png", dpi=150)
plt.close()

if len(corr_excess):
    plt.figure()
    xs = np.arange(len(corr_excess))
    plt.bar(xs, corr_excess["median_excess"])
    plt.errorbar(
        xs, corr_excess["median_excess"],
        yerr=[
            corr_excess["median_excess"] - corr_excess["ci_low"],
            corr_excess["ci_high"] - corr_excess["median_excess"]
        ],
        fmt='none', capsize=3
    )
    plt.xticks(xs, [f"TRUE - {v}" for v in corr_excess["versus"]], rotation=10)
    plt.ylabel("Median excess Δ|corr|")
    plt.title("Δ|corr| excess (TRUE minus shuffle)")
    plt.tight_layout()
    plt.savefig("figs_v2/delta_corr_excess.png", dpi=150)
    plt.close()

# Write a one-pager
with open("figs_v2/EXCESS_SUMMARY.md", "w") as f:
    f.write("# Excess tests (full out_v2)\n\n")
    f.write("## Coefficient-excess (TRUE minus shuffle)\n")
    f.write(coef_excess.to_string(index=False))
    f.write("\n\n")
    f.write("## Δ|corr|-excess (TRUE minus shuffle)\n")
    f.write(corr_excess.to_string(index=False))
    f.write("\n\n")
    f.write("## Run medians\n")
    f.write(summary.to_string(index=False))
    f.write("\n\n")
    
    # Decision rubric
    ok_coef = (coef_excess["ci_low"] > 0).sum()
    ok_corr = (corr_excess["ci_low"] > 0).sum() if len(corr_excess) else 0
    
    f.write("## Decision rubric\n")
    f.write(f"- Coefficient-excess CI>0 in {ok_coef}/{len(coef_excess)} contrasts.\n")
    f.write(f"- Δ|corr|-excess CI>0 in {ok_corr}/{len(corr_excess)} contrasts.\n")
    
    if ok_coef >= 2 or ok_corr >= 1:
        f.write("\n**Interpretation:** Some position-specific signal exceeds shuffles → keep exploring (with restraint).\n")
    else:
        f.write("\n**Interpretation:** No excess over shuffles → consider pausing or pivoting.\n")

print("Wrote:")
print("  figs_v2/coef_excess.csv")
print("  figs_v2/delta_corr_excess.csv")
print("  figs_v2/run_medians.csv")
print("  figs_v2/EXCESS_SUMMARY.md")
print("  figs_v2/*.png")

