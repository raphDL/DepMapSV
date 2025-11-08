#!/usr/bin/env python3
"""
Generate SV→CNV sample mapping CSV automatically.

This script collects sample names from SV BEDPE files and CNV segments,
optionally uses sample_info.csv for ACH↔CCLE mapping, and outputs
a 2-column mapping CSV (source_name, cell_line).
"""
import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

import pandas as pd

# Import functions from sv_ingest_wgs
sys.path.insert(0, str(Path(__file__).parent.parent))
from sv_ingest_wgs import (
    canonicalize_cell,
    infer_sample_from_filename,
    autodetect_sample_column,
    load_sample_info,
    _canon_cols,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)


def collect_sv_names(sv_dir: str) -> Set[str]:
    """Collect SV sample names from BEDPE filenames (canonicalized)."""
    sv_path = Path(sv_dir)
    if sv_path.is_file():
        files = [sv_path]
    else:
        files = list(sv_path.glob("*.bedpe"))
    
    sv_names: Set[str] = set()
    for fpath in files:
        inferred = infer_sample_from_filename(fpath)
        sv_names.add(inferred)
    
    logging.info(f"Collected {len(sv_names)} unique SV sample names from {len(files)} files")
    return sv_names


def collect_cnv_names(cnv_bed: str) -> Set[str]:
    """Collect CNV sample names (autodetect column, canonicalized)."""
    cnv = pd.read_csv(cnv_bed, sep="\t", header=0)
    sample_col_idx = autodetect_sample_column(cnv)
    
    if sample_col_idx is None:
        raise ValueError(
            f"Could not autodetect sample column in CNV file {cnv_bed}. "
            f"Columns: {list(cnv.columns)}"
        )
    
    sample_col = cnv.columns[sample_col_idx]
    cnv_names = set(cnv[sample_col].astype(str).map(canonicalize_cell).unique())
    
    logging.info(f"Collected {len(cnv_names)} unique CNV sample names")
    return cnv_names


def build_mapping(
    sv_names: Set[str],
    cnv_names: Set[str],
    depmap_to_ccle: Dict[str, str],
    ccle_to_depmap: Dict[str, str],
) -> Dict[str, str]:
    """
    Build SV→CNV mapping with fallback strategies.
    
    Order:
    1. If sample_info available: SV(CCLE) → DepMap_ID if that DepMap_ID exists in CNV
    2. Exact name match SV == CNV
    3. Prefix match: SV.split('_',1)[0] equals any CNV or its prefix
    """
    mapping: Dict[str, str] = {}
    
    # Strategy 1: Use sample_info if available
    if ccle_to_depmap:
        for sv_name in sv_names:
            # Try to map SV (CCLE) → DepMap_ID
            depmap_id = ccle_to_depmap.get(sv_name)
            if depmap_id and depmap_id in cnv_names:
                mapping[sv_name] = depmap_id
                continue
    
    # Strategy 2: Exact match
    for sv_name in sv_names:
        if sv_name not in mapping and sv_name in cnv_names:
            mapping[sv_name] = sv_name
    
    # Strategy 3: Prefix match
    for sv_name in sv_names:
        if sv_name not in mapping:
            sv_prefix = sv_name.split("_", 1)[0] if "_" in sv_name else sv_name
            for cnv_name in cnv_names:
                cnv_prefix = cnv_name.split("_", 1)[0] if "_" in cnv_name else cnv_name
                if sv_prefix == cnv_prefix or sv_prefix in cnv_name or cnv_name.startswith(sv_prefix):
                    mapping[sv_name] = cnv_name
                    break
    
    return mapping


def main():
    ap = argparse.ArgumentParser(
        description="Generate SV→CNV sample mapping CSV automatically"
    )
    ap.add_argument("--sv-dir", required=True,
                    help="Directory with WGS SV BEDPE files")
    ap.add_argument("--cnv", required=True,
                    help="CNV segments BED file")
    ap.add_argument("--out", required=True,
                    help="Output mapping CSV path")
    ap.add_argument("--sample-info", type=str, default=None,
                    help="Optional DepMap sample_info.csv for ACH↔CCLE mapping")
    args = ap.parse_args()
    
    # Collect names
    logging.info("Collecting SV sample names...")
    sv_names = collect_sv_names(args.sv_dir)
    
    logging.info("Collecting CNV sample names...")
    cnv_names = collect_cnv_names(args.cnv)
    
    # Load sample_info if provided
    depmap_to_ccle = {}
    ccle_to_depmap = {}
    if args.sample_info:
        depmap_to_ccle, ccle_to_depmap = load_sample_info(args.sample_info)
    
    # Build mapping
    logging.info("Building SV→CNV mapping...")
    mapping = build_mapping(sv_names, cnv_names, depmap_to_ccle, ccle_to_depmap)
    
    # Write CSV
    if mapping:
        df = pd.DataFrame({
            "source_name": list(mapping.keys()),
            "cell_line": list(mapping.values())
        })
        df = df.sort_values("source_name")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out, index=False)
        logging.info(f"Wrote {len(mapping)} mappings to {args.out}")
    else:
        logging.warning("No mappings found; writing empty CSV")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=["source_name", "cell_line"]).to_csv(args.out, index=False)
    
    # Print summary
    mapped_sv = set(mapping.keys())
    mapped_cnv = set(mapping.values())
    predicted_overlap = len(mapped_cnv & cnv_names)
    
    print("\n=== Mapping Summary ===")
    print(f"SV unique samples: {len(sv_names)}")
    print(f"CNV unique samples: {len(cnv_names)}")
    print(f"Mapped SV→CNV pairs: {len(mapping)}")
    print(f"Predicted overlap: {predicted_overlap}")
    print(f"Unmapped SV samples: {len(sv_names - mapped_sv)}")
    if len(sv_names - mapped_sv) > 0 and len(sv_names - mapped_sv) <= 10:
        print(f"  Examples: {sorted(list(sv_names - mapped_sv))[:5]}")


if __name__ == "__main__":
    main()

