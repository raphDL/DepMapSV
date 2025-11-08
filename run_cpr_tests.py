import os, json

import numpy as np, pandas as pd

BASE="out_v2"

def read_corr(label):
    p = os.path.join(BASE, label, "gene_level_correlations.csv")
    df = pd.read_csv(p)
    # guard: compute delta_abs_corr if missing
    if "delta_abs_corr" not in df.columns and {"r_before","r_after"}.issubset(df.columns):
        df["delta_abs_corr"] = df["r_before"].abs() - df["r_after"].abs()
    df["label"] = label
    return df[["gene","r_before","r_after","delta_abs_corr","label"]]

def read_coef(label):
    p = os.path.join(BASE, label, "model_coefficients.csv")
    df = pd.read_csv(p)
    df["label"] = label
    keep = ["gene","coef_bp_near","coef_bp_far","coef_cn","r2","rmse","n_obs","label"]
    return df[keep]

def read_dep_corr(label):
    # dependency + corrected medians per gene
    p = os.path.join(BASE, label, "dependency_corrected.csv")
    df = pd.read_csv(p, usecols=["gene","dependency","dependency_corrected"])
    g = df.groupby("gene").agg(dep_median=("dependency","median"),
                               depcorr_median=("dependency_corrected","median")).reset_index()
    g["delta_gene_median"] = g["depcorr_median"] - g["dep_median"]
    g["label"] = label
    return g

# 1) Δ|corr| stratified by baseline |r|
true_corr = read_corr("true")
within_corr = read_corr("shuffle_within_chrom")
across_corr = read_corr("shuffle_across_chrom")

# quartiles on baseline |r| from TRUE only, align all sets
true_corr["abs_r_before"] = true_corr["r_before"].abs()
q = true_corr["abs_r_before"].quantile([0.25,0.5,0.75]).to_dict()

def bin_absr(x):
    if x <= q[0.25]: return "Q1 (lowest |r|)"
    if x <= q[0.5]: return "Q2"
    if x <= q[0.75]: return "Q3"
    return "Q4 (highest |r|)"

true_corr["absr_bin"] = true_corr["abs_r_before"].apply(bin_absr)
within_corr = within_corr.merge(true_corr[["gene","absr_bin"]], on="gene", how="inner")
across_corr = across_corr.merge(true_corr[["gene","absr_bin"]], on="gene", how="inner")

def summarize_corr(df,label):
    out=[]
    for b, sub in df.groupby("absr_bin"):
        out.append(dict(label=label, bin=b,
                        n=len(sub),
                        median_delta=sub["delta_abs_corr"].median(),
                        mean_delta=sub["delta_abs_corr"].mean()))
    return pd.DataFrame(out)

corr_strat = pd.concat([
    summarize_corr(true_corr,"true"),
    summarize_corr(within_corr,"shuffle_within"),
    summarize_corr(across_corr,"shuffle_across")
], ignore_index=True)

corr_strat = corr_strat.sort_values(["bin","label"])
corr_strat.to_csv("figs_full/corr_drop_stratified_by_baseline.csv", index=False)

# 2) Proximity coefficient tails/medians
true_coef = read_coef("true")
within_coef = read_coef("shuffle_within_chrom")
across_coef = read_coef("shuffle_across_chrom")

def coef_summary(df,label):
    v = df[["coef_bp_near","coef_bp_far"]].abs().max(axis=1)
    return pd.DataFrame([{
        "label":label,
        "n_genes":len(df),
        "median_abs_prox_coef":float(v.median()),
        "p90_abs_prox_coef":float(v.quantile(0.90)),
        "p99_abs_prox_coef":float(v.quantile(0.99))
    }])

coef_stats = pd.concat([
    coef_summary(true_coef,"true"),
    coef_summary(within_coef,"shuffle_within"),
    coef_summary(across_coef,"shuffle_across")
], ignore_index=True)

coef_stats.to_csv("figs_full/proximity_coef_stats.csv", index=False)

# 3) Nonessential "de-essentialization" test (gene-median shift)
# Expect: for nonessentials, delta_gene_median (after - before) > shuffled (i.e., more + shift)
true_dep = read_dep_corr("true")
within_dep = read_dep_corr("shuffle_within_chrom")
across_dep = read_dep_corr("shuffle_across_chrom")

# Load gene sets if available; otherwise fallback to STRINGENT set from DepMap repo if present
noness_path = "AchillesNonessentialControls.csv"
ess_path = "CRISPRInferredCommonEssentials.csv"

def extract_gene_symbol(s):
    """Extract gene symbol from formats like 'ABCG8 (64241)' or 'ABCG8'"""
    if pd.isna(s):
        return None
    s = str(s).strip()
    # Split on space and take first part
    return s.split()[0] if s else None

noness = set()
ess = set()

if os.path.exists(noness_path):
    df_noness = pd.read_csv(noness_path)
    # Get first column (skip header)
    gene_col = df_noness.columns[0]
    noness = set(df_noness[gene_col].apply(extract_gene_symbol).dropna())

if os.path.exists(ess_path):
    df_ess = pd.read_csv(ess_path)
    # Get first column (skip header)
    gene_col = df_ess.columns[0]
    ess = set(df_ess[gene_col].apply(extract_gene_symbol).dropna())

def tag_set(df):
    if noness:
        df["is_nonessential"] = df["gene"].isin(noness).astype(int)
    else:
        df["is_nonessential"] = 0
    if ess:
        df["is_essential"] = df["gene"].isin(ess).astype(int)
    else:
        df["is_essential"] = 0
    return df

true_dep = tag_set(true_dep)
within_dep = tag_set(within_dep)
across_dep = tag_set(across_dep)

def summarize_shift(df,label):
    rows=[]
    for name, col in [("all", df["delta_gene_median"]),
                      ("nonessential_only", df.loc[df["is_nonessential"]==1,"delta_gene_median"]),
                      ("essential_only", df.loc[df["is_essential"]==1,"delta_gene_median"])]:
        if len(col)>0:
            rows.append(dict(label=label, subset=name,
                             n=len(col),
                             median_shift=float(col.median()),
                             mean_shift=float(col.mean()),
                             p25=float(col.quantile(0.25)),
                             p75=float(col.quantile(0.75))))
    return pd.DataFrame(rows)

shift_stats = pd.concat([
    summarize_shift(true_dep,"true"),
    summarize_shift(within_dep,"shuffle_within"),
    summarize_shift(across_dep,"shuffle_across")
], ignore_index=True)

shift_stats.to_csv("figs_full/gene_median_shifts.csv", index=False)

# 4) Simple GO/NO-GO scoreboard (no heavy plots)
# Criteria:
#  A) Corr-drop grows with baseline |r| and is larger in TRUE than SHUFFLES at Q4
#  B) Prox coef tails (p90/p99) larger in TRUE than SHUFFLES
#  C) Nonessential median shift more positive in TRUE than SHUFFLES

def pick_row(df, label, bin_name):
    sub = df[(df["label"]==label)&(df["bin"]==bin_name)]
    return None if sub.empty else float(sub.iloc[0]["median_delta"])

q4_true = pick_row(corr_strat,"true","Q4 (highest |r|)")
q4_within = pick_row(corr_strat,"shuffle_within","Q4 (highest |r|)")
q4_across = pick_row(corr_strat,"shuffle_across","Q4 (highest |r|)")

p90_true = float(coef_stats.loc[coef_stats["label"]=="true","p90_abs_prox_coef"].values[0])
p90_within = float(coef_stats.loc[coef_stats["label"]=="shuffle_within","p90_abs_prox_coef"].values[0])
p90_across = float(coef_stats.loc[coef_stats["label"]=="shuffle_across","p90_abs_prox_coef"].values[0])

ns_true = shift_stats[(shift_stats["label"]=="true")&(shift_stats["subset"]=="nonessential_only")]
ns_within = shift_stats[(shift_stats["label"]=="shuffle_within")&(shift_stats["subset"]=="nonessential_only")]
ns_across = shift_stats[(shift_stats["label"]=="shuffle_across")&(shift_stats["subset"]=="nonessential_only")]

ns_true_med = float(ns_true["median_shift"].values[0]) if not ns_true.empty else np.nan
ns_within_med = float(ns_within["median_shift"].values[0]) if not ns_within.empty else np.nan
ns_across_med = float(ns_across["median_shift"].values[0]) if not ns_across.empty else np.nan

score = pd.DataFrame([
    dict(check="A_Q4_corrdrop_true>shuffles",
         true=q4_true, within=q4_within, across=q4_across,
         pass_= (q4_true is not None and q4_within is not None and q4_across is not None and q4_true > q4_within and q4_true > q4_across)),
    dict(check="B_coef_p90_true>shuffles",
         true=p90_true, within=p90_within, across=p90_across,
         pass_= (p90_true > p90_within and p90_true > p90_across)),
    dict(check="C_nonessential_median_shift_true>shuffles",
         true=ns_true_med, within=ns_within_med, across=ns_across_med,
         pass_= (not np.isnan(ns_true_med) and not np.isnan(ns_within_med) and not np.isnan(ns_across_med)
                 and ns_true_med > ns_within_med and ns_true_med > ns_across_med))
])

score.to_csv("figs_full/go_nogo_scoreboard.csv", index=False)

print("\n=== Δ|corr(dep,CN)| stratified by baseline |r| (TRUE vs SHUFFLES) ===")
print(corr_strat.to_string(index=False))

print("\n=== Proximity coefficient distribution (TRUE vs SHUFFLES) ===")
print(coef_stats.to_string(index=False))

print("\n=== Gene-median shift after correction (TRUE vs SHUFFLES) ===")
print(shift_stats.to_string(index=False))

print("\n=== GO/NO-GO scoreboard ===")
print(score.to_string(index=False))

