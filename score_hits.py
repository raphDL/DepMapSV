#!/usr/bin/env python3
# Rank pilot hits by proximity-specific signal and assign GREEN/AMBER/RED labels.
import sys, os, argparse
import numpy as np
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("setdir", help="Path to set folder, e.g. out_pilot_final/unstable")
    ap.add_argument("--contrib-thresh", type=float, default=0.01,
                    help="Threshold for |proximity contribution| to count a cell as 'active' (default: 0.01)")
    ap.add_argument("--active-cell-frac", type=float, default=0.10,
                    help="Fraction of cells that must be active for a gene to be considered 'active' (default: 0.10)")
    ap.add_argument("--cn2-low", type=float, default=1.8, help="CN≈2 lower bound (default: 1.8)")
    ap.add_argument("--cn2-high", type=float, default=2.2, help="CN≈2 upper bound (default: 2.2)")
    ap.add_argument("--min-cells", type=int, default=30, help="Min cells per gene (default: 30)")
    ap.add_argument("--min-cells-lineage", type=int, default=25, help="Min cells per gene×lineage (default: 25)")
    ap.add_argument("--min-cells-cn2", type=int, default=25, help="Min cells required for CN≈2 check (default: 25)")
    ap.add_argument("--ignore-lineage", action="store_true", help="Ignore lineage breadth in scoring even if lineage column exists")
    ap.add_argument("--lineage-col", default="lineage", help="Column name for lineage in design_matrix.csv (default: lineage)")
    args = ap.parse_args()

    setdir = args.setdir.rstrip("/")
    dm_path = os.path.join(setdir, "design_matrix.csv")
    coef_path = os.path.join(setdir, "models_coefficients.csv")
    if not os.path.exists(dm_path) or not os.path.exists(coef_path):
        sys.exit(f"ERROR: missing inputs. Need {dm_path} and {coef_path}")

    dm = pd.read_csv(dm_path)
    coef = pd.read_csv(coef_path).rename(columns={
        "cn":"b_cn",
        "bp_within_100000":"b_w100k",
        "bp_within_1000000":"b_w1m",
        "inv_bp":"b_inv"
    })

    # Choose feature columns that match the model scale (std if present; else raw)
    has_std = all(c in dm.columns for c in ["bp_within_100000_std","bp_within_1000000_std"]) or ("inv_bp_std" in dm.columns)
    x_w100k = "bp_within_100000_std" if "bp_within_100000_std" in dm.columns else "bp_within_100000"
    x_w1m   = "bp_within_1000000_std" if "bp_within_1000000_std" in dm.columns else "bp_within_1000000"
    x_inv   = "inv_bp_std" if "inv_bp_std" in dm.columns else ("inv_bp" if "inv_bp" in dm.columns else None)

    # Merge coefficients
    coef_keep = ["gene","b_w100k","b_w1m","b_inv"]
    for c in ["b_w100k","b_w1m","b_inv"]:
        if c not in coef.columns:
            coef[c] = 0.0
    coef = coef[coef_keep]

    df = dm.merge(coef, on="gene", how="left")

    # Proximity contribution in dependency units
    df["__x_w100k"] = df.get(x_w100k, 0.0).astype(float)
    df["__x_w1m"]   = df.get(x_w1m, 0.0).astype(float)
    df["__x_inv"]   = df.get(x_inv, 0.0).astype(float) if x_inv else 0.0
    df["prox_contrib"] = (
        df["b_w100k"] * df["__x_w100k"] +
        df["b_w1m"]   * df["__x_w1m"]   +
        df["b_inv"]   * df["__x_inv"]
    )
    df["prox_only_corrected"] = df["dependency"] - df["prox_contrib"]

    # Per-gene metrics
    rows = []
    have_lineage = (args.lineage_col in df.columns) and (not args.ignore_lineage)
    for g, sub in df.groupby("gene"):
        sub = sub.dropna(subset=["cn","dependency","prox_only_corrected"])
        if len(sub) < args.min_cells:
            continue
        r0 = sub[["cn","dependency"]].corr().iloc[0,1]
        r1 = sub[["cn","prox_only_corrected"]].corr().iloc[0,1]
        drop = abs(r0) - abs(r1)
        # activity
        act_mask = (sub["prox_contrib"].abs() >= args.contrib_thresh)
        active_frac = float(act_mask.mean())
        med_contrib_if_active = float(sub.loc[act_mask, "prox_contrib"].abs().median()) if act_mask.any() else 0.0
        # CN≈2 check
        cn2 = sub[(sub["cn"] >= args.cn2_low) & (sub["cn"] <= args.cn2_high)]
        drop_cn2 = np.nan
        if len(cn2) >= max(args.min_cells, args.min_cells_cn2):
            r0c = cn2[["cn","dependency"]].corr().iloc[0,1]
            r1c = cn2[["cn","prox_only_corrected"]].corr().iloc[0,1]
            drop_cn2 = abs(r0c) - abs(r1c)
        # lineage breadth
        breadth = 0
        if have_lineage:
            for lin, ss in sub.groupby(args.lineage_col):
                if len(ss) < args.min_cells_lineage:
                    continue
                r0l = ss[["cn","dependency"]].corr().iloc[0,1]
                r1l = ss[["cn","prox_only_corrected"]].corr().iloc[0,1]
                if (abs(r0l) - abs(r1l)) > 0:
                    breadth += 1
        rows.append(dict(
            gene=g,
            n_cells=int(len(sub)),
            prox_active_frac=active_frac,
            prox_contrib_median_if_active=med_contrib_if_active,
            prox_only_abs_corr_drop=float(drop),
            prox_only_abs_corr_drop_cn2=float(drop_cn2) if not np.isnan(drop_cn2) else np.nan,
            lineage_breadth=int(breadth),
        ))

    res = pd.DataFrame(rows)
    if res.empty:
        sys.exit("No genes passed min cell thresholds; nothing to score.")

    # scoring
    def score_row(r, have_lineage_flag: bool):
        score = 0
        if r.prox_active_frac >= 0.20:
            score += 1
        if r.prox_contrib_median_if_active >= 0.01:
            score += 1
        if r.prox_only_abs_corr_drop > 0:
            score += 1
        if pd.notna(r.prox_only_abs_corr_drop_cn2) and r.prox_only_abs_corr_drop_cn2 > 0:
            score += 1
        if have_lineage_flag and r.lineage_breadth >= 2:
            score += 1

        # If lineage is missing, don't penalize: GREEN at >=3 instead of >=4
        green_cut = 4 if have_lineage_flag else 3
        label = "GREEN" if score >= green_cut else ("AMBER" if score >= 2 else "RED")
        return pd.Series({"interest_score": score, "label": label})

    scored = pd.concat([res, res.apply(lambda r: score_row(r, have_lineage), axis=1)], axis=1)
    scored = scored.sort_values(
        ["interest_score","prox_active_frac","prox_only_abs_corr_drop","prox_contrib_median_if_active"],
        ascending=[False, False, False, False]
    )

    out_csv = os.path.join(setdir, "hits_interest_score.csv")
    scored.to_csv(out_csv, index=False)

    cols = ["gene","interest_score","label","prox_active_frac",
            "prox_only_abs_corr_drop","prox_only_abs_corr_drop_cn2",
            "lineage_breadth","n_cells"]
    print(scored.head(12)[cols])
    print(f"Wrote {out_csv}")
    if not have_lineage:
        print("Note: lineage column not found or ignored; GREEN threshold adapted (>=3). You can merge a mapping and rerun for lineage breadth scoring.")

if __name__ == "__main__":
    main()
