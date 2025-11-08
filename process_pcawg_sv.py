#!/usr/bin/env python3
"""
Process PCAWG consensus SV data from compressed archives.

Extracts .tgz files, decompresses .bedpe.gz files, and converts them
to the format expected by sv_ingest_wgs.py by adding sample IDs from filenames.
"""
import argparse
import gzip
import logging
import os
import sys
import tarfile
from pathlib import Path
from typing import Optional

import pandas as pd

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


def extract_sample_id(filename: str) -> str:
    """
    Extract sample ID from PCAWG filename.
    
    Format: {uuid}.pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz
    Returns: {uuid} as sample ID
    """
    # Remove path and extension
    basename = Path(filename).stem  # Remove .gz
    basename = Path(basename).stem  # Remove .bedpe
    # Extract UUID (first part before first dot)
    sample_id = basename.split('.')[0]
    return sample_id


def process_pcawg_archive(tgz_path: str, out_dir: str, svtypes: Optional[list] = None) -> int:
    """
    Extract and process a PCAWG .tgz archive.
    
    Args:
        tgz_path: Path to .tgz file
        out_dir: Output directory for converted BEDPE files
        svtypes: Optional list of SV types to filter (e.g., ['TRA', 'INV'])
    
    Returns:
        Number of files processed
    """
    if svtypes is None:
        svtypes = ['TRA', 'INV']
    
    tgz_path = Path(tgz_path)
    if not tgz_path.exists():
        raise ValueError(f"Archive not found: {tgz_path}")
    
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    logging.info(f"Processing archive: {tgz_path.name}")
    
    # Extract archive to temp directory
    temp_dir = out_dir / f"temp_{tgz_path.stem}"
    temp_dir.mkdir(exist_ok=True)
    
    try:
        # Extract all .bedpe.gz files
        with tarfile.open(tgz_path, 'r:gz') as tar:
            bedpe_files = [m for m in tar.getmembers() if m.name.endswith('.bedpe.gz')]
            logging.info(f"Found {len(bedpe_files)} BEDPE files in archive")
            
            processed = 0
            for member in bedpe_files:
                try:
                    # Extract to temp location
                    tar.extract(member, path=temp_dir)
                    extracted_path = temp_dir / member.name
                    
                    # Get sample ID from filename
                    sample_id = extract_sample_id(member.name)
                    
                    # Read and process the BEDPE file
                    with gzip.open(extracted_path, 'rt') as f:
                        df = pd.read_csv(f, sep='\t')
                    
                    # Check required columns
                    required = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
                    missing = [c for c in required if c not in df.columns]
                    if missing:
                        logging.warning(f"Skipping {member.name}: missing columns {missing}")
                        continue
                    
                    # Add sample ID column
                    df['cell_line'] = sample_id
                    
                    # Handle svtype/svclass column
                    if 'svclass' in df.columns:
                        df['svtype'] = df['svclass']
                    elif 'svtype' not in df.columns:
                        df['svtype'] = 'UNKNOWN'
                    
                    # Filter by svtype if specified
                    if svtypes and 'svtype' in df.columns:
                        before = len(df)
                        df = df[df['svtype'].isin(svtypes)].copy()
                        after = len(df)
                        if before > after:
                            logging.debug(f"Filtered {member.name}: {before} -> {after} rows ({svtypes})")
                    
                    if len(df) == 0:
                        logging.debug(f"Skipping {member.name}: no SVs after filtering")
                        continue
                    
                    # Select columns for output (standard BEDPE format with header)
                    output_cols = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 
                                  'cell_line', 'svtype']
                    # Keep additional columns if present
                    for col in ['strand1', 'strand2', 'pe_support', 'sv_id']:
                        if col in df.columns:
                            output_cols.append(col)
                    
                    df_out = df[output_cols].copy()
                    
                    # Write output BEDPE file
                    out_file = out_dir / f"{sample_id}.bedpe"
                    df_out.to_csv(out_file, sep='\t', index=False)
                    
                    processed += 1
                    if processed % 100 == 0:
                        logging.info(f"Processed {processed} files...")
                    
                    # Clean up extracted file
                    extracted_path.unlink()
                    
                except Exception as e:
                    logging.error(f"Error processing {member.name}: {e}")
                    continue
        
        # Clean up temp directory
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)
        
        logging.info(f"Processed {processed} files from {tgz_path.name}")
        return processed
        
    except Exception as e:
        logging.error(f"Error processing archive {tgz_path}: {e}")
        raise


def main():
    ap = argparse.ArgumentParser(
        description="Extract and convert PCAWG consensus SV data to BEDPE format"
    )
    ap.add_argument("--tgz-dir", required=True,
                    help="Directory containing PCAWG .tgz archives")
    ap.add_argument("--out-dir", required=True,
                    help="Output directory for converted BEDPE files")
    ap.add_argument("--svtypes", nargs="+", default=["TRA", "INV"],
                    help="SV types to include (default: TRA INV)")
    ap.add_argument("-v", "--verbose", action="count", default=1,
                    help="Increase verbosity")
    args = ap.parse_args()
    
    setup_logging(args.verbose)
    
    tgz_dir = Path(args.tgz_dir)
    if not tgz_dir.exists():
        raise ValueError(f"Directory not found: {tgz_dir}")
    
    # Find all .tgz files
    tgz_files = list(tgz_dir.glob("*.tgz"))
    if not tgz_files:
        raise ValueError(f"No .tgz files found in {tgz_dir}")
    
    logging.info(f"Found {len(tgz_files)} archive(s) to process")
    
    total_processed = 0
    for tgz_file in tgz_files:
        try:
            n = process_pcawg_archive(tgz_file, args.out_dir, svtypes=args.svtypes)
            total_processed += n
        except Exception as e:
            logging.error(f"Failed to process {tgz_file}: {e}")
            continue
    
    logging.info(f"Total: Processed {total_processed} BEDPE files")
    logging.info(f"Output directory: {args.out_dir}")


if __name__ == "__main__":
    main()


