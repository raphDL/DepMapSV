import os, argparse
import pandas as pd, numpy as np
from pathlib import Path
from scipy.stats import beta
from numpy.random import default_rng

def load_dm(base):
    # Only the columns we need -> faster IO
    dm = pd.read_csv(f"{base}/design_matrix.csv",
                     usecols=["gene","cell_line","dependency","cn","bp_dist"])
    co = pd.read_csv(f"{base}/models_coefficients.csv",
                     usecols=["gene","bp_within_100000","bp_within_1000000"])
    co = co.rename(columns={"bp_within_100000":"coef_bp_near",
                            "bp_within_1000000":"coef_bp_far"})
    w1,w2 = 250_000, 2_000_000
    dm["bp_dist"] = dm["bp_dist"].fillna(5_000_000)
    dm["bp_near"] = (dm["bp_dist"] <= w1).astype(np.int8)
    dm["bp_far"]  = ((dm["bp_dist"] > w1) & (dm["bp_dist"] <= w2)).astype(np.int8)
    dm = dm.merge(co, on="gene", how="left")
    # prox-only correction using global coefs
    dm["prox_contrib"] = dm["coef_bp_near"].fillna(0)*dm["bp_near"] + dm["coef_bp_far"].fillna(0)*dm["bp_far"]
    dm["prox_only_corr"] = dm["dependency"] - dm["prox_contrib"]
    return dm

def gene_level_corr_drop(dm: pd.DataFrame) -> pd.Series:
    # One correlation drop per gene across all lines (enough points to be meaningful)
    drops=[]
    for g, sub in dm.groupby("gene"):
        a = sub[["dependency","cn"]].dropna()
        b = sub[["prox_only_corr","cn"]].dropna()
        if len(a) > 40 and len(b) > 40:
            r0 = a.corr().iloc[0,1]; r1 = b.corr().iloc[0,1]
            if pd.notna(r0) and pd.notna(r1):
                drops.append((g, abs(r0)-abs(r1)))
    return pd.Series({g:d for g,d in drops}, name="corr_drop")

def per_line_directional_fast(dm: pd.DataFrame,
                              corr_drop_by_gene: pd.Series,
                              contrib_thresh=0.01,
                              active_frac=0.15,
                              min_genes=50) -> pd.DataFrame:
    # Directional in a line = gene has global corr_drop>0 AND this line is "active" for that gene
    # "Active" evaluated per gene×line using prox_contrib in that line
    dm = dm.copy()
    dm["is_active_point"] = (dm["prox_contrib"].abs() >= contrib_thresh).astype(np.int8)
    # per gene×line activity (share across cells; usually 1 row/line, but keep generic)
    act = (dm.groupby(["cell_line","gene"])["is_active_point"].mean()
             .rename("active_frac")).reset_index()
    # attach global corr drop
    act = act.merge(corr_drop_by_gene.reset_index().rename(columns={"index":"gene"}),
                    on="gene", how="left")
    act["is_dir"] = (act["active_frac"] >= active_frac) & (act["corr_drop"] > 0)
    # aggregate per line
    per_line = (act.groupby("cell_line")["is_dir"]
                   .agg(k_dir=lambda s: int(s.sum()),
                        n_genes=lambda s: int(s.size))
                   .reset_index())
    # keep lines with enough genes
    per_line = per_line[per_line["n_genes"] >= min_genes].copy()
    return per_line

def eb_shrink(k, n):
    # Empirical-Bayes beta prior via MoM across lines
    p = (k / n).astype(float)
    m = float(p.mean())
    v = float(p.var(ddof=1)) if len(p)>1 else 1e-6
    if v <= 1e-8:
        a0, b0 = 1.0, 1.0
    else:
        s = max((m*(1-m)/v - 1), 1e-3)
        a0 = max(m*s, 1e-3)
        b0 = max((1-m)*s, 1e-3)
    post_mean = (k + a0) / (n + a0 + b0)
    return post_mean, a0, b0

def split_half(dm, reps=50, contrib_thresh=0.01, active_frac=0.15):
    # Quick reliability proxy: split genes within each line and correlate rates
    rng = default_rng(0)
    out=[]
    gene_by_line = dm.groupby("cell_line")["gene"].nunique()
    lines = gene_by_line[gene_by_line >= 120].index.tolist()
    if not lines:
        return pd.DataFrame(columns=["cell_line","split_half_rho","n_reps"])
    # Precompute global drops once
    corr_drop = gene_level_corr_drop(dm)
    for cl in lines:
        sub = dm[dm["cell_line"]==cl]
        genes = sub["gene"].unique().tolist()
        rhos=[]
        for _ in range(reps):
            rng.shuffle(genes)
            g1 = set(genes[:len(genes)//2]); g2 = set(genes[len(genes)//2:])
            s1 = per_line_directional_fast(sub[sub["gene"].isin(g1)], corr_drop,
                                           contrib_thresh, active_frac, min_genes=40)
            s2 = per_line_directional_fast(sub[sub["gene"].isin(g2)], corr_drop,
                                           contrib_thresh, active_frac, min_genes=40)
            if len(s1)==0 or len(s2)==0: 
                continue
            r1 = (s1["k_dir"]/s1["n_genes"]).values
            r2 = (s2["k_dir"]/s2["n_genes"]).values
            if len(r1)==len(r2):
                rhos.append(pd.Series(r1).corr(pd.Series(r2), method="spearman"))
        if rhos:
            rhos = [x for x in rhos if np.isfinite(x)]
            if rhos:
                out.append((cl, float(np.mean(rhos)), len(rhos)))
    return pd.DataFrame(out, columns=["cell_line","split_half_rho","n_reps"])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--true", default="out_focus_true/unstable")
    ap.add_argument("--rotate", default="out_focus_rotate/unstable")
    ap.add_argument("--contrib-thresh", type=float, default=0.01)
    ap.add_argument("--active-frac", type=float, default=0.15)
    ap.add_argument("--min-genes", type=int, default=50)
    ap.add_argument("--draws", type=int, default=1000, help="posterior draws per line")
    ap.add_argument("--reps", type=int, default=50, help="split-half repetitions")
    ap.add_argument("--no-reliability", action="store_true", help="skip split-half")
    ap.add_argument("--outdir", default="figs_full")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    rng = default_rng(0)

    print("Loading TRUE…")
    dmT = load_dm(args.true)
    print("Loading ROTATE…")
    dmR = load_dm(args.rotate)

    # Global gene-level directional signal (drop>0) from TRUE once
    print("Computing global gene-level correlation drops (TRUE)…")
    corr_drop_T = gene_level_corr_drop(dmT)

    # Per-line directional rates (TRUE vs ROTATE)
    print("Per-line directional rates…")
    lt = per_line_directional_fast(dmT, corr_drop_T,
                                   args.contrib_thresh, args.active_frac, args.min_genes)
    lr = per_line_directional_fast(dmR, corr_drop_T,   # same drop thresholding
                                   args.contrib_thresh, args.active_frac, args.min_genes)
    res = lt.merge(lr, on="cell_line", suffixes=("_true","_rot"), how="inner")
    if len(res)==0:
        print("No overlapping lines with sufficient genes. Exiting.")
        return

    # EB shrinkage
    res["rate_true_raw"] = res["k_dir_true"]/res["n_genes_true"]
    res["rate_rot_raw"]  = res["k_dir_rot"]/res["n_genes_rot"]
    res["rate_true_eb"], aT, bT = eb_shrink(res["k_dir_true"], res["n_genes_true"])
    res["rate_rot_eb"],  aR, bR = eb_shrink(res["k_dir_rot"],  res["n_genes_rot"])

    # Posterior difference via beta draws (vectorized)
    D = args.draws
    s_true = beta.rvs(res["k_dir_true"]+aT, res["n_genes_true"]+bT, size=(D,len(res)), random_state=rng)
    s_rot  = beta.rvs(res["k_dir_rot"] +aR, res["n_genes_rot"] +bR, size=(D,len(res)), random_state=rng)
    diffs = s_true - s_rot
    res["post_diff_mean"] = diffs.mean(axis=0)
    res["post_diff_lo"]   = np.percentile(diffs, 2.5, axis=0)
    res["post_diff_hi"]   = np.percentile(diffs, 97.5, axis=0)
    res["post_prob_gt0"]  = (diffs > 0).mean(axis=0)

    # Optional split-half reliability (quick)
    if not args.no_reliability:
        print("Split-half reliability…")
        rel = split_half(dmT, reps=args.reps,
                         contrib_thresh=args.contrib_thresh,
                         active_frac=args.active_frac)
        res = res.merge(rel, on="cell_line", how="left")

    # Shortlist
    short = res[(res["post_prob_gt0"]>=0.9) & (res.get("split_half_rho", pd.Series(dtype=float))>=0.2)].copy()
    res.to_csv(outdir/"psi_line_eb_full.csv", index=False)
    short.to_csv(outdir/"psi_line_eb_shortlist.csv", index=False)
    print("\nSaved:")
    print(" -", outdir/"psi_line_eb_full.csv")
    print(" -", outdir/"psi_line_eb_shortlist.csv")
    print(f"EB priors: TRUE a={aT:.2f}, b={bT:.2f} | ROT a={aR:.2f}, b={bR:.2f}")
    print(f"Shortlist lines (n={len(short)}):", list(short["cell_line"].head(10)))

if __name__ == "__main__":
    main()

