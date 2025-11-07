#!/usr/bin/env python3
"""
Bootstrap confidence intervals and permutation p-values for pilot metrics.

Usage:
    python pilot_bootstrap_ci.py out_pilot_final/unstable out_pilot_final/stable
"""

import pandas as pd
import numpy as np
import sys
from scipy import stats

def prox_active_frac(dm, coef, contrib_thresh=0.01, cell_frac_thresh=0.10):
    """Compute proximity-active genes fraction."""
    df = dm.merge(coef[["gene","b_w100k","b_w1m","b_inv"]], on="gene", how="left")
    df["prox_contrib"] = (
        df["b_w100k"]*df["bp_within_100000"] +
        df["b_w1m"]  *df["bp_within_1000000"] +
        df["b_inv"]  *df.get("inv_bp",0.0)
    )
    frac = (df.assign(active=(df["prox_contrib"].abs()>=contrib_thresh))
              .groupby("gene")["active"].mean())
    return (frac>=cell_frac_thresh).mean()

def per_gene_drop_median(dm, coef):
    """Compute median |Δr| across genes."""
    df = dm.merge(coef[["gene","b_w100k","b_w1m","b_inv"]], on="gene", how="left")
    df["prox_contrib"] = (
        df["b_w100k"]*df["bp_within_100000"] +
        df["b_w1m"]  *df["bp_within_1000000"] +
        df["b_inv"]  *df.get("inv_bp",0.0)
    )
    df["prox_only_corrected"] = df["dependency"] - df["prox_contrib"]
    drops = []
    for g, sub in df.groupby("gene"):
        a = sub[["cn","dependency"]].dropna()
        b = sub[["cn","prox_only_corrected"]].dropna()
        if len(a) > 10 and len(b) > 10:
            r0 = a.corr().iloc[0,1]
            r1 = b.corr().iloc[0,1]
            if not (np.isnan(r0) or np.isnan(r1)):
                drops.append(abs(r0) - abs(r1))
    return np.median(drops) if drops else np.nan

def bootstrap_ci(values, n_boot=1000, ci=0.95):
    """Bootstrap confidence interval."""
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    boot_medians = []
    for _ in range(n_boot):
        boot_sample = np.random.choice(values, size=len(values), replace=True)
        boot_medians.append(np.median(boot_sample))
    boot_medians = np.array(boot_medians)
    alpha = 1 - ci
    lower = np.percentile(boot_medians, 100 * alpha/2)
    upper = np.percentile(boot_medians, 100 * (1 - alpha/2))
    return np.median(values), lower, upper

def permutation_test(unstable_values, stable_values, n_perm=10000):
    """Permutation test for difference in medians."""
    if len(unstable_values) == 0 or len(stable_values) == 0:
        return np.nan
    observed_diff = np.median(unstable_values) - np.median(stable_values)
    combined = np.concatenate([unstable_values, stable_values])
    n_unstable = len(unstable_values)
    perm_diffs = []
    for _ in range(n_perm):
        np.random.shuffle(combined)
        perm_unstable = combined[:n_unstable]
        perm_stable = combined[n_unstable:]
        perm_diffs.append(np.median(perm_unstable) - np.median(perm_stable))
    p_value = np.mean(np.array(perm_diffs) >= observed_diff)
    return observed_diff, p_value

def load_set(setdir):
    """Load design matrix and coefficients for a set."""
    dm = pd.read_csv(f"{setdir}/design_matrix.csv")
    coef = pd.read_csv(f"{setdir}/models_coefficients.csv")
    coef = coef.rename(columns={
        "cn":"b_cn",
        "bp_within_100000":"b_w100k",
        "bp_within_1000000":"b_w1m",
        "inv_bp":"b_inv"
    })
    return dm, coef

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python pilot_bootstrap_ci.py <unstable_dir> <stable_dir>")
        sys.exit(1)
    
    unstable_dir = sys.argv[1]
    stable_dir = sys.argv[2]
    
    dm_u, coef_u = load_set(unstable_dir)
    dm_s, coef_s = load_set(stable_dir)
    
    # Metric 1: Prox-active genes fraction
    fa_u = prox_active_frac(dm_u, coef_u)
    fa_s = prox_active_frac(dm_s, coef_s)
    diff_fa = fa_u - fa_s
    
    # Bootstrap CI for difference (treat genes as units)
    # For this, we'd need per-gene binary indicators, then bootstrap genes
    # Simplified: report point estimates
    print("=== Proximity-active genes fraction ===")
    print(f"Unstable: {fa_u:.3f}")
    print(f"Stable:   {fa_s:.3f}")
    print(f"Difference: {diff_fa:.3f}")
    print("(Bootstrap CI on difference: implement by bootstrapping genes)")
    print()
    
    # Metric 2: Prox-only |Δr| median
    d_u = per_gene_drop_median(dm_u, coef_u)
    d_s = per_gene_drop_median(dm_s, coef_s)
    
    # Get per-gene drops for bootstrap/permutation
    def get_per_gene_drops(dm, coef):
        df = dm.merge(coef[["gene","b_w100k","b_w1m","b_inv"]], on="gene", how="left")
        df["prox_contrib"] = (
            df["b_w100k"]*df["bp_within_100000"] +
            df["b_w1m"]  *df["bp_within_1000000"] +
            df["b_inv"]  *df.get("inv_bp",0.0)
        )
        df["prox_only_corrected"] = df["dependency"] - df["prox_contrib"]
        drops = []
        for g, sub in df.groupby("gene"):
            a = sub[["cn","dependency"]].dropna()
            b = sub[["cn","prox_only_corrected"]].dropna()
            if len(a) > 10 and len(b) > 10:
                r0 = a.corr().iloc[0,1]
                r1 = b.corr().iloc[0,1]
                if not (np.isnan(r0) or np.isnan(r1)):
                    drops.append(abs(r0) - abs(r1))
        return np.array(drops)
    
    drops_u = get_per_gene_drops(dm_u, coef_u)
    drops_s = get_per_gene_drops(dm_s, coef_s)
    
    print("=== Prox-only median |Δr| ===")
    print(f"Unstable median: {d_u:.4f} (n={len(drops_u)} genes)")
    print(f"Stable median:   {d_s:.4f} (n={len(drops_s)} genes)")
    
    # Bootstrap CI for unstable median
    if len(drops_u) > 0:
        med_u, ci_lower_u, ci_upper_u = bootstrap_ci(drops_u)
        print(f"Unstable 95% CI: [{ci_lower_u:.4f}, {ci_upper_u:.4f}]")
    
    # Bootstrap CI for stable median
    if len(drops_s) > 0:
        med_s, ci_lower_s, ci_upper_s = bootstrap_ci(drops_s)
        print(f"Stable 95% CI:   [{ci_lower_s:.4f}, {ci_upper_s:.4f}]")
    
    # Permutation test for difference
    if len(drops_u) > 0 and len(drops_s) > 0:
        obs_diff, p_value = permutation_test(drops_u, drops_s)
        print(f"Difference (unstable - stable): {obs_diff:.4f}")
        print(f"Permutation p-value: {p_value:.4f}")

