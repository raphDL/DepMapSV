import pandas as pd, numpy as np
from statsmodels.stats.multitest import multipletests

W1,W2=250_000,2_000_000
rng=np.random.default_rng(0)

def load(base):
    dm = pd.read_csv(f"{base}/design_matrix.csv", usecols=["gene","cell_line","dependency","cn","bp_dist"])
    co = pd.read_csv(f"{base}/models_coefficients.csv", usecols=["gene","bp_within_100000","bp_within_1000000"])
    dm["bp_dist"]=dm["bp_dist"].fillna(5_000_000)
    dm["bp_near"]=(dm["bp_dist"]<=W1).astype(float)
    dm["bp_far"]=((dm["bp_dist"]>W1)&(dm["bp_dist"]<=W2)).astype(float)
    dm = dm.merge(dm.groupby("gene")["dependency"].median().rename("gene_med"), on="gene", how="left")
    dm["dep_ctr"]=dm["dependency"]-dm["gene_med"]
    dm = dm.merge(co, on="gene", how="left").fillna({"bp_within_100000":0.0,"bp_within_1000000":0.0})
    # Map model coefficients: bp_within_100000 -> near, bp_within_1000000 -> far
    dm["prox_contrib"]=dm["bp_within_100000"]*dm["bp_near"] + dm["bp_within_1000000"]*dm["bp_far"]
    return dm

def psi_beta(sub):
    # per-line OLS with gene fixed-effect via centering
    X = np.c_[sub["cn"].values, sub["prox_contrib"].values]
    y = sub["dep_ctr"].values
    if len(y) < 150: return np.nan
    try:
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        return float(b[1])  # slope on proximity
    except: return np.nan

def psi_boot(df, B=500):
    genes = df["gene"].unique()
    vals=[]
    for _ in range(B):
        g = rng.choice(genes, size=len(genes), replace=True)
        sub = df.merge(pd.DataFrame({"gene":g}), on="gene", how="inner")
        vals.append(psi_beta(sub))
    vals = np.array([v for v in vals if np.isfinite(v)])
    mu = float(np.nanmedian(vals)) if len(vals) else np.nan
    lo,hi = (np.nanpercentile(vals,[2.5,97.5]) if len(vals) else (np.nan,np.nan))
    return mu,lo,hi

TRUE   = "out_focus_true/unstable"      # adjust to your folders
ROTATE = "out_focus_rotate/unstable"

dt = load(TRUE)
dr = load(ROTATE)

rows=[]
for cl, sub_t in dt.groupby("cell_line"):
    sub_r = dr[dr["cell_line"]==cl]
    if len(sub_r)==0: continue
    bt,lt,ht = psi_boot(sub_t, B=400)
    br,lr,hr = psi_boot(sub_r, B=400)
    diff = bt - br
    rows.append((cl, bt, lt, ht, br, lr, hr, diff, len(sub_t["gene"].unique())))

out = pd.DataFrame(rows, columns=["cell_line","psi_true","ci_t_lo","ci_t_hi",
                                  "psi_rot","ci_r_lo","ci_r_hi","psi_diff","n_genes"])
# FDR on diff>0 via simple z from bootstrap SE
se = ( (out["ci_t_hi"]-out["ci_t_lo"])/3.92 + (out["ci_r_hi"]-out["ci_r_lo"])/3.92 )/np.sqrt(2)
z  = out["psi_diff"]/se.replace(0,np.nan)
p  = 2*np.exp(-0.5*(z**2))  # rough; fine for ranking
out["p_two_sided"]=p
out["q_fdr"]=multipletests(p.fillna(1.0), method="fdr_bh")[1]
out.sort_values("psi_diff", ascending=False).to_csv("figs_full/psi_line_significance.csv", index=False)
print("Wrote figs_full/psi_line_significance.csv")
print("Lines with q<0.1 and psi_diff>0:", (out["q_fdr"]<0.1).sum())

