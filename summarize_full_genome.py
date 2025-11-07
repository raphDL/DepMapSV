#!/usr/bin/env python3
import os, glob, numpy as np, pandas as pd

ACTIVE_FRAC = 0.15  # <- change to 0.20 for stricter

def contrib_thresh(dm_cols):
    # 0.01 when standardized columns exist
    return 0.01 if any(c.endswith("_std") for c in dm_cols) else 0.1

def directional_flags(setdir, active_frac=ACTIVE_FRAC):
    dm = pd.read_csv(os.path.join(setdir,"design_matrix.csv"))
    coef = pd.read_csv(os.path.join(setdir,"models_coefficients.csv")).rename(
        columns={"cn":"b_cn","bp_within_100000":"b_w100k","bp_within_1000000":"b_w1m","inv_bp":"b_inv"}
    )
    x_w100k = "bp_within_100000_std" if "bp_within_100000_std" in dm.columns else "bp_within_100000"
    x_w1m   = "bp_within_1000000_std" if "bp_within_1000000_std" in dm.columns else "bp_within_1000000"
    x_inv   = "inv_bp_std" if "inv_bp_std" in dm.columns else ("inv_bp" if "inv_bp" in dm.columns else None)

    df = dm.merge(coef[["gene","b_w100k","b_w1m","b_inv"]], on="gene", how="left")
    df["__x_w100k"] = df.get(x_w100k, 0.0).astype(float)
    df["__x_w1m"]   = df.get(x_w1m, 0.0).astype(float)
    df["__x_inv"]   = df.get(x_inv, 0.0).astype(float) if x_inv else 0.0
    df["prox_contrib"] = df["b_w100k"]*df["__x_w100k"] + df["b_w1m"]*df["__x_w1m"] + df["b_inv"]*df["__x_inv"]
    df["prox_only_corr"] = df["dependency"] - df["prox_contrib"]

    THR = contrib_thresh(dm.columns)
    out=[]
    for g, sub in df.groupby("gene", sort=False):
        sub = sub.dropna(subset=["cn","dependency","prox_only_corr"])
        if len(sub)==0:
            out.append((g, 0.0, np.nan, False))
            continue
        active = (sub["prox_contrib"].abs() >= THR).mean()
        r0 = sub[["cn","dependency"]].corr().iloc[0,1]
        r1 = sub[["cn","prox_only_corr"]].corr().iloc[0,1]
        drop = abs(r0) - abs(r1)
        flag = (active >= active_frac) and (drop > 0)
        out.append((g, active, drop, flag))
    return pd.DataFrame(out, columns=["gene","active_frac","drop","dir_flag"])

# Collect per-chunk flags
def collect(root):
    rows=[]
    for sub in sorted(glob.glob(os.path.join(root,"fullgenes_*"))):
        d=os.path.join(sub,"unstable")  # pilot path layout
        if not os.path.exists(os.path.join(d,"design_matrix.csv")): 
            continue
        f = directional_flags(d)
        f["chunk"]=os.path.basename(sub)
        rows.append(f)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["gene","active_frac","drop","dir_flag","chunk"])

print("Collecting directional flags from non-rotated chunks...")
lin   = collect("out_full_linear")

print("Collecting directional flags from rotate shuffle chunks...")
rot   = collect("out_full_linear_rotate")

# Pair and summarize
df = lin.merge(rot[["gene","dir_flag"]].rename(columns={"dir_flag":"dir_flag_rotate"}), on="gene", how="inner")
df["pairwise_excess"] = df["dir_flag"].astype(int) - df["dir_flag_rotate"].astype(int)

# Save full genome flags
os.makedirs("figs_full", exist_ok=True)
df.to_csv("figs_full/full_directional_flags.csv", index=False)

# Topline: mean directional and pairwise excess with bootstrap CIs
rng = np.random.default_rng(1)
def ci_mean(vec, B=10000):
    vec = np.array(vec)
    bs = np.mean(vec[rng.integers(0, len(vec), size=(B, len(vec)))], axis=1)
    return float(np.mean(vec)), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))

m_dir, lo_dir, hi_dir = ci_mean(df["dir_flag"].astype(int).values)
m_rot, lo_rot, hi_rot = ci_mean(df["dir_flag_rotate"].astype(int).values)
m_exc, lo_exc, hi_exc = ci_mean(df["pairwise_excess"].astype(int).values)

pd.DataFrame([
    ["linear_pairwise", "directional", m_dir, lo_dir, hi_dir, len(df)],
    ["linear_pairwise", "rotate",      m_rot, lo_rot, hi_rot, len(df)],
    ["linear_pairwise", "excess",      m_exc, lo_exc, hi_exc, len(df)],
], columns=["setting","metric","mean","ci_low","ci_high","n_genes"]).to_csv("figs_full/full_directional_summary.csv", index=False)

print("Wrote figs_full/full_directional_flags.csv and figs_full/full_directional_summary.csv")
print("\nSummary:")
print(pd.read_csv("figs_full/full_directional_summary.csv"))

