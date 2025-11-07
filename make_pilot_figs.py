import pandas as pd, matplotlib.pyplot as plt, numpy as np, os

def load(setdir):
    dm = pd.read_csv(f"{setdir}/design_matrix.csv")
    coef = pd.read_csv(f"{setdir}/models_coefficients.csv")
    coef = coef.rename(columns={"cn":"b_cn","bp_within_100000":"b_w100k",
                               "bp_within_1000000":"b_w1m","inv_bp":"b_inv"})
    return dm, coef

def per_gene_drop(dm, coef):
    df = dm.merge(coef[["gene","b_w100k","b_w1m","b_inv"]], on="gene", how="left")
    df["prox_contrib"] = (df["b_w100k"]*df["bp_within_100000"] +
                          df["b_w1m"]  *df["bp_within_1000000"] +
                          df["b_inv"]  *df.get("inv_bp",0.0))
    df["prox_only_corrected"] = df["dependency"] - df["prox_contrib"]
    drops=[]
    for g,sub in df.groupby("gene"):
        a=sub[["cn","dependency"]].dropna()
        b=sub[["cn","prox_only_corrected"]].dropna()
        if len(a)>10 and len(b)>10:
            r0=a.corr().iloc[0,1]; r1=b.corr().iloc[0,1]
            if not (np.isnan(r0) or np.isnan(r1)):
                drops.append(abs(r0)-abs(r1))
    return pd.Series(drops)

def prox_active_frac(dm, coef, contrib_thresh=0.01, cell_frac_thresh=0.10):
    df = dm.merge(coef[["gene","b_w100k","b_w1m","b_inv"]], on="gene", how="left")
    df["prox_contrib"] = (df["b_w100k"]*df["bp_within_100000"] +
                          df["b_w1m"]  *df["bp_within_1000000"] +
                          df["b_inv"]  *df.get("inv_bp",0.0))
    frac = (df.assign(active=(df["prox_contrib"].abs()>=contrib_thresh))
              .groupby("gene")["active"].mean())
    return (frac>=cell_frac_thresh).mean()

os.makedirs("figs_pilot", exist_ok=True)

sets={"unstable":"out_pilot_final/unstable","stable":"out_pilot_final/stable"}
dm_u, coef_u = load(sets["unstable"])
dm_s, coef_s = load(sets["stable"])

# Fig1: Proximity-active genes fraction (robust metric)
fa_u = prox_active_frac(dm_u, coef_u)
fa_s = prox_active_frac(dm_s, coef_s)
plt.figure(); plt.bar(["unstable","stable"], [fa_u, fa_s])
plt.ylabel("Genes with ≥10% cells |prox contrib|≥0.01"); plt.tight_layout()
plt.savefig("figs_pilot/prox_active_frac.pdf")
plt.close()

# Fig2: |Δ corr(dep,CN)| after removing ONLY proximity
d_u = per_gene_drop(dm_u, coef_u); d_s = per_gene_drop(dm_s, coef_s)
plt.figure(); plt.boxplot([d_u.dropna(), d_s.dropna()], labels=["unstable","stable"])
plt.ylabel("Median |Δr| (CN vs dep), prox-only"); plt.tight_layout()
plt.savefig("figs_pilot/prox_only_corr_drop.pdf")
plt.close()

# Fig3: Top-5 genes by |Δr| (unstable)
top = d_u.sort_values(ascending=False).head(5)
pd.DataFrame({"gene":top.index, "abs_corr_drop":top.values}).to_csv("figs_pilot/top5_unstable.csv", index=False)
print("Saved: figs_pilot/*.png + figs_pilot/top5_unstable.csv")

