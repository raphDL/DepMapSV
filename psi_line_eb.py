import os, pandas as pd, numpy as np
from pathlib import Path
from scipy.stats import beta
from numpy.random import default_rng

BASE_T = "out_focus_true/unstable"
BASE_R = "out_focus_rotate/unstable"
outdir = Path("figs_full"); outdir.mkdir(exist_ok=True, parents=True)
rng = default_rng(0)

def load_dm(base):
    dm = pd.read_csv(f"{base}/design_matrix.csv")
    co = pd.read_csv(f"{base}/models_coefficients.csv")
    # windows = 250k/2Mb (matches your focused pilot)
    w1, w2 = 250_000, 2_000_000
    dm["bp_near"] = (dm["bp_dist"].fillna(5_000_000) <= w1).astype(int)
    dm["bp_far"]  = ((dm["bp_dist"].fillna(5_000_000) > w1) & (dm["bp_dist"] <= w2)).astype(int)
    # Use actual column names from models_coefficients.csv
    # Rename columns for clarity
    co_renamed = co[["gene","bp_within_100000","bp_within_1000000"]].rename(
        columns={"bp_within_100000": "coef_bp_near", "bp_within_1000000": "coef_bp_far"}
    )
    dm = dm.merge(co_renamed, on="gene", how="left")
    dm["prox_contrib"] = dm["coef_bp_near"].fillna(0)*dm["bp_near"] + dm["coef_bp_far"].fillna(0)*dm["bp_far"]
    dm["prox_only_corr"] = dm["dependency"] - dm["prox_contrib"]
    return dm

def per_line_directional(dm, contrib_thresh=0.01, active_frac=0.15, min_genes=50):
    rows=[]
    for cl, sub in dm.groupby("cell_line"):
        # gene-level flags in this line
        grows=[]
        for g,sg in sub.groupby("gene"):
            a = sg[["dependency","cn"]].dropna()
            b = sg[["prox_only_corr","cn"]].dropna()
            if len(a)>20 and len(b)>20:
                r0 = a.corr().iloc[0,1]; r1 = b.corr().iloc[0,1]
                if pd.notna(r0) and pd.notna(r1):
                    drop = abs(r0) - abs(r1)
                    is_active = (sg["prox_contrib"].abs() >= contrib_thresh).mean() >= active_frac
                    is_dir = is_active and (drop > 0)
                    grows.append((g, is_dir))
        if len(grows) >= min_genes:
            arr = np.array([x[1] for x in grows], int)
            rows.append((cl, int(arr.sum()), int(arr.size)))
    return pd.DataFrame(rows, columns=["cell_line","k_dir","n_genes"])

def eb_shrink(k, n, a0=None, b0=None):
    # Empirical-Bayes beta prior via method-of-moments across lines
    p = (k / n).astype(float)
    m, v = p.mean(), p.var(ddof=1) if len(p)>1 else 1e-6
    if a0 is None or b0 is None:
        # Avoid degenerate prior
        if v <= 1e-8:
            a0, b0 = 1.0, 1.0
        else:
            # Beta variance v = m(1-m)/(a+b+1)
            s = max((m*(1-m)/v - 1), 1e-3)
            a0 = max(m*s, 1e-3)
            b0 = max((1-m)*s, 1e-3)
    post_mean = (k + a0) / (n + a0 + b0)
    return post_mean, a0, b0

def split_half_reliability(dm, reps=200, contrib_thresh=0.01, active_frac=0.15):
    # Spearman across splits of genes per line; returns per-line mean rho
    out=[]
    for cl, sub in dm.groupby("cell_line"):
        gene_list = list(sub["gene"].unique())
        if len(gene_list) < 120:  # need ~60 per half
            continue
        rhos=[]
        for _ in range(reps):
            rng.shuffle(gene_list)
            g1 = set(gene_list[:len(gene_list)//2]); g2 = set(gene_list[len(gene_list)//2:])
            d1 = sub[sub["gene"].isin(g1)]
            d2 = sub[sub["gene"].isin(g2)]
            s1 = per_line_directional(d1, contrib_thresh, active_frac, min_genes=40)
            s2 = per_line_directional(d2, contrib_thresh, active_frac, min_genes=40)
            if len(s1)==0 or len(s2)==0: 
                continue
            # PSI rate per half
            s1r = (s1["k_dir"]/s1["n_genes"]).values
            s2r = (s2["k_dir"]/s2["n_genes"]).values
            if len(s1r)==len(s2r):
                rho = pd.Series(s1r).corr(pd.Series(s2r), method="spearman")
                if np.isfinite(rho):
                    rhos.append(rho)
        if rhos:
            out.append((cl, float(np.mean(rhos)), len(rhos)))
    return pd.DataFrame(out, columns=["cell_line","split_half_rho","n_reps"])

# Load
dmT = load_dm(BASE_T)
dmR = load_dm(BASE_R)

# Per-line directional counts/rates
lt = per_line_directional(dmT)
lr = per_line_directional(dmR)
res = lt.merge(lr, on="cell_line", how="inner", suffixes=("_true","_rot"))

# EB shrinkage for TRUE and ROTATE separately
res["rate_true_raw"] = res["k_dir_true"] / res["n_genes_true"]
res["rate_rot_raw"]  = res["k_dir_rot"]  / res["n_genes_rot"]
res["rate_true_eb"], aT, bT = eb_shrink(res["k_dir_true"], res["n_genes_true"])
res["rate_rot_eb"],  aR, bR = eb_shrink(res["k_dir_rot"],  res["n_genes_rot"])

# Posterior difference with credible interval via beta posterior draws
D = 5000
post_diffs=[]
for _ in range(D):
    s_true = beta.rvs(res["k_dir_true"]+aT, res["n_genes_true"]+bT, random_state=rng)
    s_rot  = beta.rvs(res["k_dir_rot"] +aR, res["n_genes_rot"] +bR, random_state=rng)
    post_diffs.append(s_true - s_rot)
post_diffs = np.stack(post_diffs, axis=1)  # lines x draws
res["post_diff_mean"] = post_diffs.mean(axis=1)
res["post_diff_lo"]   = np.percentile(post_diffs, 2.5, axis=1)
res["post_diff_hi"]   = np.percentile(post_diffs, 97.5, axis=1)
res["post_prob_gt0"]  = (post_diffs>0).mean(axis=1)

# Split-half reliability (optional: a bit slow but uses current data only)
print("Computing split-half reliability (this may take a few minutes)...")
rel = split_half_reliability(dmT, reps=200)
res = res.merge(rel, on="cell_line", how="left")

# Shortlist: high posterior probability and decent reliability
short = res[(res["post_prob_gt0"]>=0.9) & (res["split_half_rho"]>=0.2)].copy().sort_values("post_diff_mean", ascending=False)

res.to_csv(outdir/"psi_line_eb_full.csv", index=False)
short.to_csv(outdir/"psi_line_eb_shortlist.csv", index=False)

print("\nSaved:")
print(" - figs_full/psi_line_eb_full.csv")
print(" - figs_full/psi_line_eb_shortlist.csv")
print(f"EB priors: TRUE a={aT:.2f}, b={bT:.2f} | ROT a={aR:.2f}, b={bR:.2f}")
print(f"Shortlist lines (n={len(short)}):", list(short["cell_line"].head(10)))

