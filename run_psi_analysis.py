import os, numpy as np, pandas as pd
from pathlib import Path

BASE = os.environ.get("BASE", "out_pilot_baseline_fixed/unstable")
CNV  = os.environ.get("CNV",  "prep_out/cnv_segments.bed")
SV   = os.environ.get("SV",   "prep_out/sv_from_cnv.bedpe")
W1   = int(os.environ.get("W1", "250000"))
W2   = int(os.environ.get("W2", "2000000"))

outdir = Path(BASE).parent / "psi_line"
outdir.mkdir(parents=True, exist_ok=True)

# --- Load training design (has dependency, cn, bp_dist) ---
dm_path = Path(BASE) / "design_matrix.csv"
if not dm_path.exists():
    raise SystemExit(f"Missing {dm_path}. Use a pilot/focus run that writes design_matrix.csv.")

dm = pd.read_csv(dm_path)
# Safety: keep minimal columns
need = {"gene","cell_line","dependency","cn","bp_dist"}
missing = need - set(dm.columns)
if missing:
    raise SystemExit(f"design_matrix.csv is missing columns: {missing}")

# Flags and gene baseline
dm["bp_near"] = (dm["bp_dist"].fillna(5_000_000) <= W1).astype(int)
dm["bp_far"]  = ((dm["bp_dist"].fillna(5_000_000) > W1) & (dm["bp_dist"] <= W2)).astype(int)
gene_median = dm.groupby("gene")["dependency"].median()
dm = dm.merge(gene_median.rename("gene_med_dep"), left_on="gene", right_index=True, how="left")
dm["dep_centered"] = dm["dependency"] - dm["gene_med_dep"]

# --- Per-line PSI: regression beta and partial correlations ---
def partial_corr(y, x, z):
    # residualize y~z and x~z; corr(res_y, res_x)
    y = np.asarray(y); x=np.asarray(x); Z = np.c_[np.ones(len(z)), np.asarray(z)]
    # Ridge guard (tiny) for numeric stability
    lam = 1e-8
    ZtZ = Z.T@Z + lam*np.eye(Z.shape[1])
    beta_y = np.linalg.solve(ZtZ, Z.T@y)
    beta_x = np.linalg.solve(ZtZ, Z.T@x)
    ry = y - Z@beta_y
    rx = x - Z@beta_x
    denom = (np.linalg.norm(ry)*np.linalg.norm(rx))
    return float((ry@rx)/denom) if denom>0 else np.nan

rows=[]
for cl, sub in dm.groupby("cell_line"):
    # Drop NAs and require enough genes
    s=sub.dropna(subset=["dep_centered","cn"])
    if len(s) < 100:  # need enough genes per line for stable stats
        continue

    # Design: [CN, bp_near, bp_far]
    X = s[["cn","bp_near","bp_far"]].values.astype(float)
    y = s["dep_centered"].values.astype(float)
    # Add intercept
    X1 = np.c_[np.ones(len(s)), X]
    # OLS closed-form (robust SEs skipped; we bootstrap below)
    try:
        beta = np.linalg.lstsq(X1, y, rcond=None)[0]
    except Exception:
        continue
    b_cn, b_near, b_far = beta[1], beta[2], beta[3]

    # Partial correlations (dep_centered ⟂ bp_* | CN)
    pr_near = partial_corr(s["dep_centered"], s["bp_near"], s["cn"])
    pr_far  = partial_corr(s["dep_centered"], s["bp_far"],  s["cn"])

    rows.append((cl, len(s), b_cn, b_near, b_far, pr_near, pr_far))

psi = pd.DataFrame(rows, columns=["cell_line","n_genes",
                                  "beta_cn","beta_bp_near","beta_bp_far",
                                  "partial_r_near","partial_r_far"])
psi["psi_beta_abs"] = psi[["beta_bp_near","beta_bp_far"]].abs().max(axis=1)
psi["psi_part_abs"] = psi[["partial_r_near","partial_r_far"]].abs().max(axis=1)
psi.sort_values("psi_beta_abs", ascending=False).to_csv(outdir/"psi_by_line.csv", index=False)
print("Wrote", outdir/"psi_by_line.csv")

# --- Instability covariates: breakpoint count and FGA ---
# Breakpoints per line
sv = pd.read_csv(SV, sep="\t")
bp_counts = sv.groupby("cell_line").size().rename("breakpoints")

# FGA per line (fraction of genome altered; crude: |CN-2|>0.3 weighted by segment length)
cnv = pd.read_csv(CNV, sep="\t")
cnv = cnv.assign(seg_len=(cnv["end"]-cnv["start"]).clip(lower=1))
cnv["is_alt"] = (cnv["cn"].astype(float)-2.0).abs() > 0.3
fga = (cnv.groupby("cell_line")
          .apply(lambda d: float(d.loc[d["is_alt"],"seg_len"].sum())/float(d["seg_len"].sum()) if d["seg_len"].sum()>0 else np.nan)
          .rename("FGA"))

cov = pd.concat([bp_counts, fga], axis=1).reset_index().rename(columns={"index":"cell_line"})
res = psi.merge(cov, on="cell_line", how="left")

# Simple associations
def safecorr(a,b):
    a=pd.to_numeric(a, errors='coerce'); b=pd.to_numeric(b, errors='coerce')
    a=np.asarray(a); b=np.asarray(b)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum()<10: return np.nan
    return float(np.corrcoef(a[m],b[m])[0,1])

c1 = safecorr(res["psi_beta_abs"], res["breakpoints"])
c2 = safecorr(res["psi_beta_abs"], res["FGA"])
c3 = safecorr(res["psi_part_abs"], res["breakpoints"])
c4 = safecorr(res["psi_part_abs"], res["FGA"])

res.to_csv(outdir/"psi_with_instability.csv", index=False)
print("Wrote", outdir/"psi_with_instability.csv")
print(f"corr(psi_beta_abs, breakpoints) = {c1:.3f}")
print(f"corr(psi_beta_abs, FGA)         = {c2:.3f}")
print(f"corr(psi_part_abs, breakpoints) = {c3:.3f}")
print(f"corr(psi_part_abs, FGA)         = {c4:.3f}")

# Top lines to inspect
top = res.sort_values("psi_beta_abs", ascending=False).head(12)["cell_line"].tolist()
print("Top PSI lines:", top)

