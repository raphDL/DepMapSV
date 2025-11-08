#!/usr/bin/env python3
"""Generate ROTATE and optionally WITHIN shuffles."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def recompute_kernels(df: pd.DataFrame) -> pd.DataFrame:
    """Recompute proximity kernels from bp_dist."""
    d = df["bp_dist"].astype("int64")
    df["prox_exp_50k"] = np.exp(-d / 5e4)
    df["prox_exp_100k"] = np.exp(-d / 1e5)
    df["prox_exp_250k"] = np.exp(-d / 2.5e5)
    df["prox_exp_500k"] = np.exp(-d / 5e5)
    df["prox_inv"] = 1.0 / (d + 1e3)
    return df


def rotate_by_cell(df: pd.DataFrame, seed: int = 1337) -> pd.DataFrame:
    """Rotate bp_dist within each cell_line (circular shift)."""
    rng = np.random.default_rng(seed)
    
    def _rot(g: pd.DataFrame) -> pd.DataFrame:
        n = len(g)
        if n <= 1 or g["bp_dist"].std() < 1e-9:
            return g
        shift = int(rng.integers(1, n))
        g = g.copy()
        g["bp_dist"] = np.roll(g["bp_dist"].to_numpy(), shift)
        return g
    
    return df.groupby("cell_line", group_keys=False).apply(_rot)


def shuffle_within_cell(df: pd.DataFrame, seed: int = 1337) -> pd.DataFrame:
    """Shuffle bp_dist within each cell_line (permutation)."""
    rng = np.random.default_rng(seed)
    
    def _perm(g: pd.DataFrame) -> pd.DataFrame:
        g = g.copy()
        g["bp_dist"] = rng.permutation(g["bp_dist"].to_numpy())
        return g
    
    return df.groupby("cell_line", group_keys=False).apply(_perm)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)      # TRUE design (parquet)
    ap.add_argument("--out-rotate", required=True)
    ap.add_argument("--out-within", default="")
    ap.add_argument("--do-rotate", action="store_true")
    ap.add_argument("--do-within", action="store_true")
    ap.add_argument("--seed", type=int, default=1337)
    args = ap.parse_args()
    
    # Read design (Parquet or CSV fallback)
    if args.design.endswith(".parquet"):
        try:
            df = pd.read_parquet(args.design)
        except ImportError:
            csv_path = args.design.replace(".parquet", ".csv")
            df = pd.read_csv(csv_path)
            print(f"[shuffle] Parquet not available, using CSV: {csv_path}")
    else:
        df = pd.read_csv(args.design)
    base_cols = ["gene", "cell_line", "cn", "bp_dist"]
    keep = [c for c in df.columns if c in base_cols]
    df = df[keep].copy()
    
    if args.do_rotate:
        r = rotate_by_cell(df, seed=args.seed)
        r = recompute_kernels(r)
        Path(Path(args.out_rotate).parent).mkdir(parents=True, exist_ok=True)
        try:
            r.to_parquet(args.out_rotate, index=False)
            print(f"[shuffle] Wrote ROTATE design to {args.out_rotate}")
        except ImportError:
            csv_out = args.out_rotate.replace(".parquet", ".csv")
            r.to_csv(csv_out, index=False)
            print(f"[shuffle] Parquet not available, wrote ROTATE design to {csv_out}")
    
    if args.do_within and args.out_within:
        w = shuffle_within_cell(df, seed=args.seed + 1)
        w = recompute_kernels(w)
        Path(Path(args.out_within).parent).mkdir(parents=True, exist_ok=True)
        try:
            w.to_parquet(args.out_within, index=False)
            print(f"[shuffle] Wrote WITHIN design to {args.out_within}")
        except ImportError:
            csv_out = args.out_within.replace(".parquet", ".csv")
            w.to_csv(csv_out, index=False)
            print(f"[shuffle] Parquet not available, wrote WITHIN design to {csv_out}")


if __name__ == "__main__":
    main()

