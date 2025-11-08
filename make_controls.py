#!/usr/bin/env python3
"""
Make Controls - Generate rotate and across-chrom shuffle controls.

Controls:
- rotate: Rotate SV positions within chromosome (preserve per-cell breakpoint counts)
- shuffle: Shuffle across chromosomes (preserve per-cell breakpoint counts and segment lengths)

Pre-registered: Preserve per-cell SV counts and segment lengths where applicable.
"""
import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# hg38 chromosome sizes (autosomes)
CHR_SIZES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468
}


def setup_logging(verbosity: int = 1):
    """Setup logging with verbosity control."""
    level = logging.WARNING if verbosity <= 0 else (logging.INFO if verbosity == 1 else logging.DEBUG)
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S",
        handlers=[logging.StreamHandler(sys.stderr)],
        force=True,
    )


def norm_chr(s: str) -> str:
    """Normalize chromosome names."""
    s = str(s).strip()
    if s.startswith("chr"):
        return s
    s_clean = s.replace("chr", "").replace("Chr", "").replace("CHR", "")
    return "chr" + s_clean


def load_wgs_sv(sv_dir: str) -> pd.DataFrame:
    """Load WGS SV BEDPE files."""
    sv_path = Path(sv_dir)
    if sv_path.is_file():
        files = [sv_path]
    else:
        files = list(sv_path.glob("*.bedpe"))
        if not files:
            raise ValueError(f"No BEDPE files found in {sv_dir}")
    
    all_sv = []
    for fpath in files:
        df = pd.read_csv(fpath, sep="\t")
        required = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "cell_line"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns in {fpath}: {missing}")
        all_sv.append(df)
    
    sv = pd.concat(all_sv, ignore_index=True)
    sv["cell_line"] = sv["cell_line"].astype(str)
    sv["chrom1"] = sv["chrom1"].map(norm_chr)
    sv["chrom2"] = sv["chrom2"].map(norm_chr)
    return sv


def rotate_within_chrom(sv: pd.DataFrame, seed: int = 1) -> pd.DataFrame:
    """
    Rotate SV positions within chromosome (preserve per-cell breakpoint counts).
    
    For each cell × chromosome, add a random offset modulo chromosome length.
    """
    rng = np.random.default_rng(seed)
    sv_rot = sv.copy()
    
    logging.info("Rotating SV positions within chromosomes...")
    
    # Process chrom1
    for (cell, chrom), idx in sv_rot.groupby(["cell_line", "chrom1"]).groups.items():
        L = CHR_SIZES.get(str(chrom), None)
        if L is None or len(idx) == 0:
            continue
        off = int(rng.integers(0, L))
        s = (sv_rot.loc[idx, "start1"].astype(int).to_numpy() + off) % L
        e = (sv_rot.loc[idx, "end1"].astype(int).to_numpy() + off) % L
        sv_rot.loc[idx, "start1"] = np.minimum(s, e)
        sv_rot.loc[idx, "end1"] = np.maximum(s, e)
    
    # Process chrom2
    for (cell, chrom), idx in sv_rot.groupby(["cell_line", "chrom2"]).groups.items():
        L = CHR_SIZES.get(str(chrom), None)
        if L is None or len(idx) == 0:
            continue
        off = int(rng.integers(0, L))
        s = (sv_rot.loc[idx, "start2"].astype(int).to_numpy() + off) % L
        e = (sv_rot.loc[idx, "end2"].astype(int).to_numpy() + off) % L
        sv_rot.loc[idx, "start2"] = np.minimum(s, e)
        sv_rot.loc[idx, "end2"] = np.maximum(s, e)
    
    logging.info(f"Rotated {len(sv_rot)} SV breakpoints")
    return sv_rot


def shuffle_across_chrom(sv: pd.DataFrame, seed: int = 1) -> pd.DataFrame:
    """
    Shuffle SV breakpoints across chromosomes (preserve per-cell breakpoint counts and segment lengths).
    
    For each cell, shuffle breakpoints across all chromosomes while preserving:
    - Per-cell total breakpoint count
    - Segment lengths (end - start)
    """
    rng = np.random.default_rng(seed)
    sv_shuf = sv.copy()
    
    logging.info("Shuffling SV breakpoints across chromosomes...")
    
    # Process chrom1
    for cell, cell_sv in sv_shuf.groupby("cell_line"):
        # Collect all chrom1 breakpoints for this cell
        chrom1_data = []
        for _, row in cell_sv.iterrows():
            chrom1_data.append({
                "idx": row.name,
                "chrom": row["chrom1"],
                "start": int(row["start1"]),
                "end": int(row["end1"]),
                "seg_len": int(row["end1"]) - int(row["start1"])
            })
        
        if not chrom1_data:
            continue
        
        # Shuffle chromosomes, preserving segment lengths
        chroms = np.array([d["chrom"] for d in chrom1_data])
        seg_lens = np.array([d["seg_len"] for d in chrom1_data])
        
        # Randomly assign to chromosomes (weighted by chromosome size)
        chrom_names = list(CHR_SIZES.keys())
        chrom_sizes = np.array([CHR_SIZES[c] for c in chrom_names])
        probs = chrom_sizes / chrom_sizes.sum()
        
        new_chroms = rng.choice(chrom_names, size=len(chrom1_data), p=probs)
        
        # Assign positions within new chromosomes
        for i, (d, new_chrom) in enumerate(zip(chrom1_data, new_chroms)):
            csize = CHR_SIZES[new_chrom]
            max_start = max(1, csize - seg_lens[i])
            new_start = rng.integers(0, max_start)
            new_end = min(new_start + seg_lens[i], csize)
            
            sv_shuf.loc[d["idx"], "chrom1"] = new_chrom
            sv_shuf.loc[d["idx"], "start1"] = new_start
            sv_shuf.loc[d["idx"], "end1"] = new_end
    
    # Process chrom2 similarly
    for cell, cell_sv in sv_shuf.groupby("cell_line"):
        chrom2_data = []
        for _, row in cell_sv.iterrows():
            chrom2_data.append({
                "idx": row.name,
                "chrom": row["chrom2"],
                "start": int(row["start2"]),
                "end": int(row["end2"]),
                "seg_len": int(row["end2"]) - int(row["start2"])
            })
        
        if not chrom2_data:
            continue
        
        chroms = np.array([d["chrom"] for d in chrom2_data])
        seg_lens = np.array([d["seg_len"] for d in chrom2_data])
        
        chrom_names = list(CHR_SIZES.keys())
        chrom_sizes = np.array([CHR_SIZES[c] for c in chrom_names])
        probs = chrom_sizes / chrom_sizes.sum()
        
        new_chroms = rng.choice(chrom_names, size=len(chrom2_data), p=probs)
        
        for i, (d, new_chrom) in enumerate(zip(chrom2_data, new_chroms)):
            csize = CHR_SIZES[new_chrom]
            max_start = max(1, csize - seg_lens[i])
            new_start = rng.integers(0, max_start)
            new_end = min(new_start + seg_lens[i], csize)
            
            sv_shuf.loc[d["idx"], "chrom2"] = new_chrom
            sv_shuf.loc[d["idx"], "start2"] = new_start
            sv_shuf.loc[d["idx"], "end2"] = new_end
    
    logging.info(f"Shuffled {len(sv_shuf)} SV breakpoints across chromosomes")
    return sv_shuf


def save_bedpe(sv: pd.DataFrame, out_path: str):
    """Save SV as BEDPE file."""
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    sv.to_csv(out_path, sep="\t", index=False)
    logging.info(f"Saved {len(sv)} SVs to {out_path}")


def main():
    ap = argparse.ArgumentParser(
        description="Generate rotate and shuffle control SV datasets"
    )
    ap.add_argument("--sv-dir", required=True,
                    help="Directory with WGS SV BEDPE files or single file path")
    ap.add_argument("--out", required=True,
                    help="Output directory (will create rotate/ and shuffle/ subdirs)")
    ap.add_argument("--seed", type=int, default=1,
                    help="Random seed (default: 1)")
    ap.add_argument("-v", "--verbose", action="count", default=1,
                    help="Increase verbosity")
    args = ap.parse_args()
    
    setup_logging(args.verbose)
    
    # Load original SV data
    logging.info("Loading WGS SV data...")
    sv_original = load_wgs_sv(args.sv_dir)
    logging.info(f"Loaded {len(sv_original)} SVs from {sv_original['cell_line'].nunique()} cell lines")
    
    # Generate rotate control
    logging.info("Generating rotate control...")
    sv_rotate = rotate_within_chrom(sv_original.copy(), seed=args.seed)
    rotate_dir = os.path.join(args.out, "rotate")
    os.makedirs(rotate_dir, exist_ok=True)
    rotate_path = os.path.join(rotate_dir, "sv_rotated.bedpe")
    save_bedpe(sv_rotate, rotate_path)
    
    # Generate shuffle control
    logging.info("Generating shuffle control...")
    sv_shuffle = shuffle_across_chrom(sv_original.copy(), seed=args.seed)
    shuffle_dir = os.path.join(args.out, "shuffle")
    os.makedirs(shuffle_dir, exist_ok=True)
    shuffle_path = os.path.join(shuffle_dir, "sv_shuffled.bedpe")
    save_bedpe(sv_shuffle, shuffle_path)
    
    logging.info("Control generation complete!")
    logging.info(f"Rotate control: {rotate_path}")
    logging.info(f"Shuffle control: {shuffle_path}")


if __name__ == "__main__":
    main()

