import pandas as pd, numpy as np
rng=np.random.default_rng(1)

def split_half_reliability(dm):
    rel=[]
    for cl, sub in dm.groupby("cell_line"):
        if sub["gene"].nunique()<100: continue
        genes = sub["gene"].unique()
        A = set(rng.choice(genes, size=len(genes)//2, replace=False))
        B = set(genes) - A
        def psi_beta(sub):
            X = np.c_[sub["cn"].values, sub["prox_contrib"].values]
            y = sub["dep_ctr"].values
            if len(y)<50: return np.nan  # Lower threshold for split-half
            try: return float(np.linalg.lstsq(X, y, rcond=None)[0][1])
            except: return np.nan
        pa = psi_beta(sub[sub["gene"].isin(A)])
        pb = psi_beta(sub[sub["gene"].isin(B)])
        rel.append((cl, pa, pb))
    r = pd.DataFrame(rel, columns=["cell_line","psi_A","psi_B"]).dropna()
    r["sh_correlation"]=np.nan if len(r)==0 else np.corrcoef(r["psi_A"], r["psi_B"])[0,1]
    return r

# reuse cached design+prox from previous snippet if you kept it; else re-create quickly
W1,W2=250_000,2_000_000
dm = pd.read_csv("out_focus_true/unstable/design_matrix.csv",
                 usecols=["gene","cell_line","dependency","cn","bp_dist"])
co = pd.read_csv("out_focus_true/unstable/models_coefficients.csv",
                 usecols=["gene","bp_within_100000","bp_within_1000000"])
dm["bp_dist"]=dm["bp_dist"].fillna(5_000_000)
dm["bp_near"]=(dm["bp_dist"]<=W1).astype(float)
dm["bp_far"]=((dm["bp_dist"]>W1)&(dm["bp_dist"]<=W2)).astype(float)
dm = dm.merge(dm.groupby("gene")["dependency"].median().rename("gene_med"), on="gene", how="left")
dm["dep_ctr"]=dm["dependency"]-dm["gene_med"]
dm = dm.merge(co, on="gene", how="left").fillna({"bp_within_100000":0.0,"bp_within_1000000":0.0})
# Map model coefficients: bp_within_100000 -> near, bp_within_1000000 -> far
dm["prox_contrib"]=dm["bp_within_100000"]*dm["bp_near"] + dm["bp_within_1000000"]*dm["bp_far"]

r = split_half_reliability(dm)
r.to_csv("figs_full/psi_line_split_half.csv", index=False)
print("Wrote figs_full/psi_line_split_half.csv; split-half r:",
      (np.corrcoef(r['psi_A'], r['psi_B'])[0,1] if len(r)>2 else np.nan))

