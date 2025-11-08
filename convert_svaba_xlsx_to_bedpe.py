#!/usr/bin/env python3
"""
Convert CCLE SvABA translocations Excel to BEDPE format.

Converts CCLE_translocations_SvABA_20181221.xlsx to BEDPE format with:
- 0-based, half-open coordinates
- Point breakends: [pos-1, pos)
- Filters for TRA-like and INV-like events
- Outputs combined and per-sample BEDPE files
"""
import argparse
import logging
import os
import re
import sys
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
    handlers=[logging.StreamHandler(sys.stderr)],
)


def parse_breakpoint(bp_str: str) -> Optional[Tuple[str, int, int, str]]:
    """
    Parse breakpoint string like "13:33777850-33777850+" or "chrX:123456-123456-".
    
    Returns: (chrom, start, end, strand) or None if unparsable
    - chrom: chr-prefixed (e.g., chr13, chrX)
    - start: pos-1 (0-based)
    - end: pos (half-open, so end = start+1 for point breakends)
    - strand: + or -
    """
    if pd.isna(bp_str) or not isinstance(bp_str, str):
        return None
    
    bp_str = bp_str.strip()
    
    # Pattern: chrom:start-end(strand) or chrom:pos(strand)
    # Handle both with and without chr prefix
    # Standard format: chrom:pos-pos(strand)
    # Chromosome pattern: digits, X, Y, M, or MT
    match = re.match(r"^\s*(?:chr)?(MT|[0-9XYM]+):(\d+)-(\d+)([+-])\s*$", bp_str, re.IGNORECASE)
    if match:
        chrom_raw = match.group(1).upper()
        pos1 = int(match.group(2))
        pos2 = int(match.group(3))
        strand = match.group(4)
        # Use first position (they should be the same for point breakends)
        pos = pos1
    else:
        # Try single position format: chrom:pos(strand)
        match = re.match(r"^\s*(?:chr)?(MT|[0-9XYM]+):(\d+)([+-])\s*$", bp_str, re.IGNORECASE)
        if match:
            chrom_raw = match.group(1).upper()
            pos = int(match.group(2))
            strand = match.group(3)
        else:
            return None
    
    # Convert to chr-prefixed chromosome
    if chrom_raw.isdigit():
        chrom = f"chr{chrom_raw}"
    elif chrom_raw in ["X", "Y", "M", "MT"]:
        chrom = f"chr{chrom_raw}"
    else:
        chrom = f"chr{chrom_raw}"
    
    # BEDPE: 0-based, half-open, point breakend [pos-1, pos)
    start = max(0, pos - 1)
    end = pos
    
    return (chrom, start, end, strand)


def load_sample_map(map_csv: Optional[str]) -> dict:
    """Load sample name mapping from CSV: CCLE_name -> DepMap_ID."""
    if not map_csv or not os.path.exists(map_csv):
        return {}
    
    try:
        df = pd.read_csv(map_csv)
        if "CCLE_name" not in df.columns or "DepMap_ID" not in df.columns:
            logging.warning(f"Mapping CSV missing required columns (CCLE_name, DepMap_ID)")
            return {}
        
        return dict(zip(df["CCLE_name"].astype(str), df["DepMap_ID"].astype(str)))
    except Exception as e:
        logging.warning(f"Error loading mapping CSV: {e}")
        return {}


def convert_excel_to_bedpe(
    infile: str,
    map_csv: Optional[str] = None,
    keep_classes: list = None,
    out_combined: str = "sv_wgs_all.bedpe",
    out_dir: str = "sv_wgs_bedpe"
):
    """
    Convert SvABA Excel to BEDPE format.
    
    Args:
        infile: Path to Excel file
        map_csv: Optional CSV with CCLE_name,DepMap_ID mapping
        keep_classes: List of classes to keep (default: ['TRA-like', 'INV-like'])
        out_combined: Output path for combined BEDPE
        out_dir: Output directory for per-sample BEDPEs
    """
    if keep_classes is None:
        keep_classes = ["TRA-like", "INV-like"]
    
    logging.info(f"Reading Excel file: {infile}")
    df = pd.read_excel(infile, engine="openpyxl")
    
    total_rows = len(df)
    logging.info(f"Total rows read: {total_rows}")
    
    # Check required columns
    required_cols = ["CCLE_name", "bp1", "bp2", "class"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    # Load sample mapping
    sample_map = load_sample_map(map_csv)
    if sample_map:
        logging.info(f"Loaded {len(sample_map)} sample mappings")
    
    # Filter by class
    if "class" in df.columns:
        df_filtered = df[df["class"].isin(keep_classes)].copy()
        kept_by_class = len(df_filtered)
        logging.info(f"Rows kept by class filter ({keep_classes}): {kept_by_class}")
    else:
        df_filtered = df.copy()
        kept_by_class = len(df_filtered)
        logging.warning("No 'class' column found, keeping all rows")
    
    # Parse breakpoints
    logging.info("Parsing breakpoints...")
    bp1_parsed = df_filtered["bp1"].apply(parse_breakpoint)
    bp2_parsed = df_filtered["bp2"].apply(parse_breakpoint)
    
    # Count unparsable
    unparsable_bp1 = bp1_parsed.isna().sum()
    unparsable_bp2 = bp2_parsed.isna().sum()
    if unparsable_bp1 > 0:
        logging.warning(f"Unparsable bp1: {unparsable_bp1} rows")
    if unparsable_bp2 > 0:
        logging.warning(f"Unparsable bp2: {unparsable_bp2} rows")
    
    # Filter out unparsable
    valid_mask = bp1_parsed.notna() & bp2_parsed.notna()
    df_valid = df_filtered[valid_mask].copy()
    logging.info(f"Rows with valid breakpoints: {len(df_valid)}")
    
    # Extract parsed breakpoint components
    df_valid["chrom1"] = bp1_parsed[valid_mask].apply(lambda x: x[0] if x else None)
    df_valid["start1"] = bp1_parsed[valid_mask].apply(lambda x: x[1] if x else None)
    df_valid["end1"] = bp1_parsed[valid_mask].apply(lambda x: x[2] if x else None)
    df_valid["strand1"] = bp1_parsed[valid_mask].apply(lambda x: x[3] if x else None)
    
    df_valid["chrom2"] = bp2_parsed[valid_mask].apply(lambda x: x[0] if x else None)
    df_valid["start2"] = bp2_parsed[valid_mask].apply(lambda x: x[1] if x else None)
    df_valid["end2"] = bp2_parsed[valid_mask].apply(lambda x: x[2] if x else None)
    df_valid["strand2"] = bp2_parsed[valid_mask].apply(lambda x: x[3] if x else None)
    
    # Map sample names
    df_valid["sample_id"] = df_valid["CCLE_name"].astype(str).map(
        lambda x: sample_map.get(x, x) if sample_map else x
    )
    
    # Get SV class
    df_valid["svclass"] = df_valid["class"].astype(str)
    
    # Build info field (key=value pairs)
    info_cols = ["svclass", "gene1", "gene2", "fusion", "multi_sv_fusion", "cosmic_fus", "site1", "site2"]
    available_info_cols = [c for c in info_cols if c in df_valid.columns]
    
    def build_info(row):
        parts = []
        for col in available_info_cols:
            val = row[col]
            if pd.notna(val) and str(val).strip():
                parts.append(f"{col}={str(val).strip()}")
        return ";".join(parts) if parts else "."
    
    df_valid["info"] = df_valid.apply(build_info, axis=1)
    
    # Create BEDPE columns
    bedpe = pd.DataFrame({
        "chrom1": df_valid["chrom1"],
        "start1": df_valid["start1"],
        "end1": df_valid["end1"],
        "chrom2": df_valid["chrom2"],
        "start2": df_valid["start2"],
        "end2": df_valid["end2"],
        "name": df_valid["svclass"],  # Use class as name
        "score": 0,  # Placeholder
        "strand1": df_valid["strand1"],
        "strand2": df_valid["strand2"],
        "sample": df_valid["sample_id"],
        "info": df_valid["info"]
    })
    
    # Deduplicate by (sample_id, chrom1, start1, end1, chrom2, start2, end2, svclass)
    before_dedup = len(bedpe)
    bedpe = bedpe.drop_duplicates(
        subset=["sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "name"]
    )
    after_dedup = len(bedpe)
    logging.info(f"Deduplication: {before_dedup} -> {after_dedup} rows")
    
    # Sort
    bedpe = bedpe.sort_values(["sample", "chrom1", "start1", "chrom2", "start2"])
    
    # Sanity checks
    sanity_check(bedpe)
    
    # Write combined BEDPE
    os.makedirs(os.path.dirname(out_combined) or ".", exist_ok=True)
    bedpe.to_csv(out_combined, sep="\t", index=False, header=False)
    logging.info(f"Wrote combined BEDPE: {out_combined} ({len(bedpe)} rows)")
    
    # Write per-sample BEDPEs
    os.makedirs(out_dir, exist_ok=True)
    unique_samples = bedpe["sample"].nunique()
    logging.info(f"Writing per-sample BEDPEs to {out_dir}/...")
    
    for sample_id, sample_bedpe in bedpe.groupby("sample"):
        sample_file = os.path.join(out_dir, f"{sample_id}.bedpe")
        sample_bedpe[["chrom1", "start1", "end1", "chrom2", "start2", "end2", 
                     "name", "score", "strand1", "strand2", "info"]].to_csv(
            sample_file, sep="\t", index=False, header=False
        )
    
    logging.info(f"Wrote {unique_samples} per-sample BEDPE files")
    
    # Summary statistics
    print("\n" + "=" * 60)
    print("CONVERSION SUMMARY")
    print("=" * 60)
    print(f"Total rows read: {total_rows:,}")
    print(f"Rows kept by class filter: {kept_by_class:,}")
    print(f"Rows with valid breakpoints: {len(df_valid):,}")
    print(f"Rows after deduplication: {after_dedup:,}")
    print(f"Unique samples: {unique_samples:,}")
    print(f"Total breakpoints written: {after_dedup:,}")
    
    # Top 3 samples with most events
    top_samples = bedpe["sample"].value_counts().head(3)
    print(f"\nTop 3 samples by event count:")
    for sample, count in top_samples.items():
        print(f"  {sample}: {count} events")
    
    print("=" * 60)


def sanity_check(bedpe: pd.DataFrame):
    """Run sanity checks on BEDPE data."""
    logging.info("Running sanity checks...")
    
    # Check chromosome format
    all_chroms = pd.concat([bedpe["chrom1"], bedpe["chrom2"]]).unique()
    non_chr_prefixed = [c for c in all_chroms if not str(c).startswith("chr")]
    if non_chr_prefixed:
        raise ValueError(f"Found chromosomes without 'chr' prefix: {non_chr_prefixed[:5]}")
    logging.info("✓ All chromosomes have 'chr' prefix")
    
    # Check coordinates: 0 ≤ start < end and end = start+1 (point breakends)
    for side in ["1", "2"]:
        start_col = f"start{side}"
        end_col = f"end{side}"
        
        invalid_start = (bedpe[start_col] < 0).sum()
        if invalid_start > 0:
            raise ValueError(f"Found {invalid_start} breakends with start < 0")
        
        invalid_end = (bedpe[end_col] <= bedpe[start_col]).sum()
        if invalid_end > 0:
            raise ValueError(f"Found {invalid_end} breakends with end <= start")
        
        # Point breakends: end should equal start+1
        point_breakends = (bedpe[end_col] == bedpe[start_col] + 1).sum()
        total = len(bedpe)
        if point_breakends != total:
            logging.warning(f"Only {point_breakends}/{total} breakends are point intervals (end=start+1)")
        else:
            logging.info(f"✓ All breakends are point intervals (end=start+1)")
    
    logging.info("✓ All sanity checks passed")


def main():
    ap = argparse.ArgumentParser(
        description="Convert CCLE SvABA translocations Excel to BEDPE format"
    )
    ap.add_argument("--infile", required=True,
                    help="Input Excel file (CCLE_translocations_SvABA_20181221.xlsx)")
    ap.add_argument("--map-csv", default=None,
                    help="Optional CSV with CCLE_name,DepMap_ID mapping")
    ap.add_argument("--keep-classes", nargs="+", default=["TRA-like", "INV-like"],
                    help="SV classes to keep (default: TRA-like INV-like)")
    ap.add_argument("--out-combined", default="sv_wgs_all.bedpe",
                    help="Output combined BEDPE file (default: sv_wgs_all.bedpe)")
    ap.add_argument("--out-dir", default="sv_wgs_bedpe",
                    help="Output directory for per-sample BEDPEs (default: sv_wgs_bedpe)")
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="Verbose logging")
    
    args = ap.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    convert_excel_to_bedpe(
        infile=args.infile,
        map_csv=args.map_csv,
        keep_classes=args.keep_classes,
        out_combined=args.out_combined,
        out_dir=args.out_dir
    )


if __name__ == "__main__":
    main()

