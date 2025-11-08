#!/usr/bin/env python3
"""
SV Ingest WGS - Load WGS structural variant breakpoints and compute proximity kernels.

This module loads WGS BEDPE files, filters for TRA/INV (and optional complex),
computes nearest breakpoint distances, and generates all proximity kernels:
- Exponential: exp(-distance / λ) for λ ∈ {50k, 100k, 250k, 500k}
- Inverse: 1 / (distance + 1e3)
- Cubic spline on log10(distance) with df=3

Pre-registered: NO CN-derived breakpoints used.
"""
import argparse
import logging
import os
import re
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple, Iterable

import numpy as np
import pandas as pd

# Chromosome normalization
_CHR_REMAP = {
    "M": "chrM", "MT": "chrM", "X": "chrX", "Y": "chrY"
}


def normalize_chrom(ch: str, scheme: str = "auto") -> str:
    """
    Normalize chromosome names to chr1, chr2, ..., chrX, chrY, chrM format.
    Drops alt/random scaffolds.
    """
    if ch is None:
        return ""
    s = str(ch).strip()
    s = s.replace("chr", "")
    s = s.upper()
    if s in _CHR_REMAP:
        s = _CHR_REMAP[s]
    elif s.isdigit():
        s = f"chr{s}"
    elif s.startswith("CHR"):
        s = s  # already fine (e.g., CHRX)
    else:
        # drop alt/random scaffolds
        if "_" in s or "RANDOM" in s or "ALT" in s:
            return ""
        s = f"chr{s}"
    return s.replace("CHR", "chr")


# Cell line canonicalization
ACH_ID_RE = re.compile(r'^\s*ACH[\s_-]?(\d{6})\s*$', re.I)


def is_ach_id(s: str) -> bool:
    """Check if string is an ACH ID (ACH-######)."""
    return bool(ACH_ID_RE.match(str(s)))


def normalize_ach(s: str) -> str:
    """
    Normalize ACH IDs: "ACH000001", "ACH 000001", "ACH_000001" → "ACH-000001".
    Only normalizes when the value is truly an ACH ID. Leaves other strings untouched.
    """
    m = ACH_ID_RE.match(str(s))
    return f"ACH-{m.group(1)}".upper() if m else str(s).strip()


def normalize_ccle(s: str) -> str:
    """
    Normalize CCLE-style names: treat -, _, ., and space as equivalent.
    Keep hyphens only for ACH IDs.
    """
    s = str(s).strip()
    if is_ach_id(s):
        # Leave to normalize_ach()
        return normalize_ach(s)
    # Collapse separators to one underscore
    s = re.sub(r'[\s\.-]+', '_', s)
    return s.upper()


def key_for_prefix(s: str) -> str:
    """
    Key used for comparisons (SV vs CNV), makes - and _ interchangeable.
    """
    return re.sub(r'[-_]+', '_', str(s).upper())


def infer_sample_from_filename(p: Path) -> str:
    """Infer sample name from filename: take stem up to first dot."""
    stem = p.name
    if "." in stem:
        stem = stem.split(".", 1)[0]
    return normalize_ccle(stem)


def cell_from_filename(fp: str) -> str:
    """Extract cell line from filename."""
    stem = Path(fp).name  # keep only basename
    if stem.lower().endswith(".bedpe"):
        stem = stem[:-6]
    return normalize_ccle(stem)


_CCLE_LIKE = re.compile(r"^[A-Z0-9]+(?:_[A-Z0-9]+)+$")
_ACH_LIKE = re.compile(r"^ACH-\d{6}$")


def is_valid_sample_value(v: str) -> bool:
    """Check if a value looks like a valid sample name (CCLE or ACH pattern)."""
    if v is None:
        return False
    s = str(v)
    # Reject obvious INFO-like fields
    if ("=" in s) or (";" in s):
        return False
    s = s.strip()
    # Check if it's an ACH ID or CCLE-like pattern
    return bool(is_ach_id(s) or _CCLE_LIKE.match(s))


def pick_sample_value_from_df(df: pd.DataFrame) -> Optional[str]:
    """
    Try to find a constant, valid sample value in a column.
    Must be the same for ~all rows and pass pattern checks.
    """
    best = None
    for col in df.columns:
        vals = df[col].dropna().astype(str)
        if vals.empty:
            continue
        # Ignore very long text columns (INFO)
        if vals.str.len().median() > 40:
            continue
        # Candidate must be (almost) constant
        top = vals.value_counts(dropna=False).head(1)
        if top.empty:
            continue
        value = top.index[0]
        frac = top.iloc[0] / len(vals)
        if frac >= 0.95 and is_valid_sample_value(value):
            best = value
            break
    return normalize_ccle(best) if best else None


def _canon_cols(cols) -> Dict[str, str]:
    """Lower, strip, collapse spaces/underscores."""
    out: Dict[str, str] = {}
    for c in cols:
        k = str(c).strip().lower().replace(" ", "_")
        out[c] = k
    return out


def _pick_source_target(colmap: Dict[str, str]) -> Tuple[Optional[str], Optional[str]]:
    """Accept many aliases including CCLE_name/DepMap_ID."""
    src_aliases = {
        "source_name", "source", "sample", "from", "raw", "input", "orig", "original",
        "ccle_name", "ccle", "ccleid"
    }
    tgt_aliases = {
        "cell_line", "target", "to", "depmap", "depmap_id", "model", "name", "model_id"
    }
    inv = {v: k for k, v in colmap.items()}
    src = next((inv[a] for a in src_aliases if a in inv), None)
    tgt = next((inv[a] for a in tgt_aliases if a in inv), None)
    return src, tgt


def load_sample_map(path: Optional[str]) -> Dict[str, str]:
    """
    Load optional sample mapping CSV.
    - Accepts flexible headers (e.g., CCLE_name, DepMap_ID).
    - Ignores comment lines starting with '#'.
    - If exactly 2 columns and headers unrecognized, uses col0->col1.
    - If no usable rows, returns {} with a warning.
    """
    if not path:
        return {}
    
    try:
        df = pd.read_csv(path, comment="#")
    except Exception as e:
        logging.error("Failed to read sample map %s: %s", path, e)
        return {}
    
    if df.empty or df.shape[1] == 0:
        logging.warning("Sample map %s is empty; continuing without remap.", path)
        return {}
    
    colmap = _canon_cols(df.columns)
    df = df.rename(columns=colmap)
    
    src_col, tgt_col = _pick_source_target(colmap)
    if src_col and tgt_col:
        src_c = colmap[src_col]
        tgt_c = colmap[tgt_col]
    else:
        if df.shape[1] == 2:
            src_c, tgt_c = df.columns[0], df.columns[1]
            logging.warning(
                "Sample map headers unrecognized; assuming first column is source and second is target: %s -> %s",
                src_c, tgt_c
            )
        else:
            logging.error(
                "Sample map must have two columns (source,target). Got columns: %s",
                list(df.columns)
            )
            return {}
    
    # Drop fully missing rows
    df = df[[src_c, tgt_c]].dropna(how="any")
    if df.empty:
        logging.warning("Sample map %s has no data rows; continuing without remap.", path)
        return {}
    
    src_vals = df[src_c].astype(str).apply(normalize_ccle)
    tgt_vals = df[tgt_c].astype(str).apply(normalize_ccle)
    mapping = dict(zip(src_vals, tgt_vals))
    
    if not mapping:
        logging.warning("Sample map produced an empty mapping; continuing without remap.")
    else:
        # Log a small example for sanity
        k0 = next(iter(mapping.keys()))
        logging.info("Loaded sample map with %d entries. Example: %s -> %s", len(mapping), k0, mapping[k0])
    return mapping

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




def _read_one_bedpe(path: Path, sample_map: Dict[str, str], svtypes: list,
                    force_cell_from_filename: bool = False) -> pd.DataFrame:
    """Read one BEDPE file with robust cell_line assignment."""
    # Try reading without header first to check column count
    df_test = pd.read_csv(path, sep="\t", header=None, nrows=1, comment="#")
    n_cols = len(df_test.columns)
    
    # Read full file
    if n_cols >= 11:
        # Standard BEDPE format: no header
        df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)
        if len(df.columns) >= 11:
            # Standard BEDPE: chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2 sample [info]
            df.columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", 
                         "name", "score", "strand1", "strand2", "sample"] + \
                        ([f"col{i}" for i in range(11, len(df.columns))] if len(df.columns) > 11 else [])
            # Map extra columns to expected names if present
            if len(df.columns) > 11:
                if "col11" in df.columns:
                    df = df.rename(columns={"col11": "info"})
        else:
            raise ValueError(f"BEDPE file {path} has unexpected number of columns: {len(df.columns)}")
    else:
        # Header format
        df = pd.read_csv(path, sep="\t", comment="#")
        required = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns in {path}: {missing}")
        # Handle cell_line vs sample column
        if "cell_line" not in df.columns and "sample" in df.columns:
            df = df.rename(columns={"sample": "cell_line"})
        # Handle svtype
        if "svtype" not in df.columns and "name" in df.columns:
            df["svtype"] = df["name"].astype(str).str.replace("-like", "", regex=False)
        elif "svtype" not in df.columns:
            df["svtype"] = "UNKNOWN"
    
    # Determine cell_line for this file robustly
    if force_cell_from_filename:
        # Always use filename when forced - never attempt autodetection
        cell_fn = cell_from_filename(str(path))
        # IMPORTANT: set and log the real filename-derived value
        logging.warning(
            "%s: multiple 'sample' values detected in file; overriding with filename-derived '%s'",
            path.name, cell_fn
        )
        # Apply sample map if provided
        cell = sample_map.get(cell_fn, cell_fn)
    else:
        # Try to autodetect from dataframe
        sample_guess = pick_sample_value_from_df(df)
        cell_fn = cell_from_filename(str(path))
        if sample_guess is None:
            logging.info("%s: no reliable 'sample' column; using filename '%s'", path.name, cell_fn)
            cell = sample_map.get(cell_fn, cell_fn)
        else:
            # Use autodetected value, but apply sample map
            cell = sample_map.get(sample_guess, sample_guess)
            # Guard: if this file contains many distinct "sample" values, override with filename
            if ("sample" in df.columns) and (df["sample"].nunique() > 3):
                cell_fn = cell_from_filename(str(path))
                logging.warning(
                    "%s: multiple 'sample' values detected in file; overriding with filename-derived '%s'",
                    path.name, cell_fn
                )
                cell = sample_map.get(cell_fn, cell_fn)
    
    df["cell_line"] = cell
    df["cell_line"] = df["cell_line"].astype("string")
    
    # Extract svtype (if not already set from header format)
    if "svtype" not in df.columns:
        if "name" in df.columns:
            df["svtype"] = df["name"].astype(str).str.replace("-like", "", regex=False)
        elif "svclass" in df.columns:
            df["svtype"] = df["svclass"].astype(str).str.replace("-like", "", regex=False)
        elif "info" in df.columns:
            df["svtype"] = df["info"].astype(str).str.extract(r"svclass=([^;]+)", expand=False)
            df["svtype"] = df["svtype"].str.replace("-like", "", regex=False)
        else:
            df["svtype"] = "UNKNOWN"
    
    # Filter by svtype
    if svtypes:
        normalized_svtypes = [s.replace("-like", "") for s in svtypes]
        df = df[df["svtype"].isin(normalized_svtypes)].copy()
    
    # Normalize chromosomes
    df["chrom1"] = df["chrom1"].map(lambda c: normalize_chrom(c))
    df["chrom2"] = df["chrom2"].map(lambda c: normalize_chrom(c))
    
    return df


def load_wgs_sv(sv_dir: str, svtypes: list = None, sample_map: Dict[str, str] = None,
                force_cell_from_filename: bool = False) -> pd.DataFrame:
    """
    Load WGS SV BEDPE files from directory with robust sample name handling.
    
    Args:
        sv_dir: Directory containing BEDPE files or single file path
        svtypes: List of SV types to keep (default: ['TRA', 'INV'])
        sample_map: Optional mapping from source_name to cell_line
        force_cell_from_filename: If True, always use filename for cell_line
    
    Returns:
        DataFrame with columns: chrom1,start1,end1,chrom2,start2,end2,cell_line,svtype
    """
    if svtypes is None:
        svtypes = ['TRA', 'INV']
    if sample_map is None:
        sample_map = {}
    
    sv_path = Path(sv_dir)
    if sv_path.is_file():
        files = [sv_path]
    else:
        if not sv_path.exists():
            raise ValueError(f"SV directory does not exist: {sv_dir}")
        if not sv_path.is_dir():
            raise ValueError(f"SV path is not a directory: {sv_dir}")
        files = list(sv_path.glob("*.bedpe"))
        if not files:
            raise ValueError(f"No BEDPE files found in {sv_dir}. "
                           f"Expected files matching *.bedpe in this directory.")
    
    logging.info(f"Loading WGS SV from {len(files)} file(s)...")
    all_sv = []
    
    for fpath in files:
        try:
            df = _read_one_bedpe(fpath, sample_map, svtypes, force_cell_from_filename)
            all_sv.append(df)
        except Exception as e:
            logging.error(f"Error loading {fpath}: {e}")
            raise
    
    if not all_sv:
        raise ValueError("No SV data loaded")
    
    sv = pd.concat(all_sv, ignore_index=True)
    
    # Ensure required columns exist
    if "cell_line" not in sv.columns:
        raise ValueError("No cell_line column found in SV data")
    
    # Log SV-only counts (before any merge)
    n_sv_events = len(sv)
    n_sv_lines = sv["cell_line"].nunique()
    n_files = len(files)
    logging.info(f"Loaded {n_sv_events} WGS SV breakpoints from {n_sv_lines} SV samples (from {n_files} BEDPE files)")
    
    if n_sv_lines > n_files * 2:
        logging.warning(
            "Suspicious: SV sample count (%d) >> number of files (%d). "
            "Check sample mapping / filename inference.",
            n_sv_lines, n_files
        )
        # Print top offenders
        logging.warning("Top SV samples:\n%s", sv["cell_line"].value_counts().head(10).to_string())
    
    return sv


def load_genes_bed(genes_bed: str) -> pd.DataFrame:
    """Load gene BED file: chrom,start,end,gene,strand"""
    g = pd.read_csv(genes_bed, sep="\t", header=None, 
                    names=["chrom", "start", "end", "gene", "strand"])
    g["chrom"] = g["chrom"].map(lambda c: normalize_chrom(c))
    # Canonicalize gene names (remove trailing numbers in parentheses)
    g["gene"] = g["gene"].astype(str).str.replace(r"\s*\(\d+\)\s*$", "", regex=True).str.strip()
    return g


def autodetect_sample_column(df: pd.DataFrame) -> Optional[int]:
    """
    Autodetect sample column in a CNV table. Prefer column-name heuristics, then values.
    Returns column index or None if not found.
    """
    name_priority = ["modelid", "depmap", "cell_line", "sample", "ach"]
    for i, col in enumerate(df.columns):
        low = str(col).lower()
        if any(k in low for k in name_priority):
            return i
    # Fallback by value patterns
    ach_pattern = re.compile(r"^ACH-\d{6}$")
    for i, col in enumerate(df.columns):
        vals = df[col].astype(str).head(100)
        if (vals.str.fullmatch(ach_pattern).sum() > 50):
            return i
    return None


def autodetect_cn_column(df: pd.DataFrame) -> Tuple[int, str]:
    """
    Autodetect CN column. Returns (index, mode) where mode in {"abs","segmean"}.
    """
    # Step 1: Name-based detection (most reliable)
    candidates_abs = ["copynumber", "copy_number", "cn", "abs_cn", "total_cn"]
    candidates_seg = ["segment_mean", "segmean", "log2ratio", "log2_ratio", "log2"]
    
    for i, c in enumerate(df.columns):
        col_lower = str(c).lower().replace("_", "").replace(" ", "")
        # Exact match on known absolute CN names
        if col_lower in [s.replace("_", "") for s in candidates_abs]:
            return i, "abs"
    
    for i, c in enumerate(df.columns):
        col_lower = str(c).lower().replace("_", "").replace(" ", "")
        # Exact match on known log2/segment mean names
        if col_lower in [s.replace("_", "") for s in candidates_seg]:
            return i, "segmean"
    
    # Step 2: Numeric heuristic (fallback)
    numeric_cols = []
    for i, c in enumerate(df.columns):
        if pd.api.types.is_numeric_dtype(df[c]):
            numeric_cols.append(i)
    
    if not numeric_cols:
        raise ValueError("No numeric columns found for CN")
    
    # Compute means and choose based on range
    means = []
    for i in numeric_cols:
        vals = pd.to_numeric(df.iloc[:min(500, len(df)), i], errors='coerce')
        m = float(vals.mean()) if vals.notna().sum() > 0 else np.nan
        means.append((i, m))
    
    # Pick column with mean closest to 2.0
    i_best, m_best = min(means, key=lambda x: abs((x[1] or 0) - 2.0))
    
    # Decision rule: if mean is in [0.5, 6.0], assume absolute CN
    # Otherwise (mean near 0, or very negative/positive), assume log2
    if 0.5 <= m_best <= 6.0:
        mode = "abs"
        logging.info(f"CN autodetected as absolute (mean={m_best:.2f} ~ 2.0)")
    else:
        mode = "segmean"
        logging.info(f"CN autodetected as log2/segmean (mean={m_best:.2f} != 2.0)")
    
    return i_best, mode


def build_prefix_mapping(sv_cells: set, cnv_cells: set) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Build prefix-based mapping between CCLE-like (SV) and CNV names.
    Uses key_for_prefix() to make matching hyphen/underscore agnostic.
    """
    sv_to_cnv: Dict[str, str] = {}
    cnv_to_sv: Dict[str, str] = {}
    
    # Build lookup by normalized key
    sv_by_key: Dict[str, list] = {}
    for s in sv_cells:
        k = key_for_prefix(s)
        sv_by_key.setdefault(k, []).append(s)
    
    # Match CNV to SV
    for cn in cnv_cells:
        k = key_for_prefix(cn)
        # Match exact or "prefix_" where separators are normalized
        matches = []
        for k_sv, originals in sv_by_key.items():
            if k_sv == k or k_sv.startswith(k + "_"):
                matches.extend(originals)
        
        if len(matches) == 1:
            cnv_to_sv[cn] = matches[0]
            sv_to_cnv[matches[0]] = cn
        elif len(matches) > 1:
            # Pick shortest original (most specific), but log it
            best = min(matches, key=len)
            cnv_to_sv[cn] = best
            sv_to_cnv[best] = cn
            logging.warning("Multiple SV matches for %s → picked %s from %s", cn, best, matches)
    
    return sv_to_cnv, cnv_to_sv


def apply_name_mapping(
    sv: pd.DataFrame,
    cnv: pd.DataFrame,
    map_cnv_to: str,  # 'sv' or 'cnv'
    sample_info_df: Optional[pd.DataFrame] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    """
    Apply name mapping with *chaining*:
    - Normalize: CNV as ACH only, SV as CCLE
    - If sample_info provided (DepMap_ID <-> CCLE_name), use it
    - Then apply prefix mapping as a refinement, without overwriting good mappings
    """
    sv = sv.copy()
    cnv = cnv.copy()
    
    # Normalize SV sample names (from filename or BEDPE column) → normalize_ccle
    sv["cell_line"] = sv["cell_line"].apply(normalize_ccle)
    
    # Normalize CNV: ACH only first
    cnv["cell_line"] = cnv["cell_line"].apply(normalize_ach)
    
    # Step 1: if sample_info present, try ACH <-> CCLE mapping
    ach_to_ccle_dict = {}
    ccle_to_ach_dict = {}
    if sample_info_df is not None:
        # Try to find DepMap_ID and CCLE_name columns (flexible)
        colmap = _canon_cols(sample_info_df.columns)
        df_renamed = sample_info_df.rename(columns=colmap)
        
        depmap_col = None
        ccle_col = None
        for alias in ["depmap_id", "depmap", "model_id", "modelid", "model", "ach_id"]:
            if alias in colmap.values():
                depmap_col = alias
                break
        for alias in ["ccle_name", "ccle", "ccleid", "stripped_cell_line_name", "strippedcelllinename"]:
            if alias in colmap.values():
                ccle_col = alias
                break
        
        if depmap_col and ccle_col:
            for _, r in df_renamed[[depmap_col, ccle_col]].dropna().iterrows():
                dep = normalize_ach(r[depmap_col])  # Normalize ACH ID
                ccl = normalize_ccle(r[ccle_col])  # Normalize CCLE name
                if dep and ccl:
                    ach_to_ccle_dict[dep] = ccl
                    ccle_to_ach_dict[ccl] = dep
    
    if map_cnv_to == "sv":
        # Map CNV (ACH) → CCLE using sample_info, if possible
        if ach_to_ccle_dict:
            cnv["cell_line"] = cnv["cell_line"].map(ach_to_ccle_dict).fillna(cnv["cell_line"])
            # Now normalize the mapped CCLE names
            cnv["cell_line"] = cnv["cell_line"].apply(normalize_ccle)
        
        # Step 2: prefix refinement (only fill where still unmatched)
        sv_to_cnv, cnv_to_sv = build_prefix_mapping(
            set(sv["cell_line"].unique()),
            set(cnv["cell_line"].unique())
        )
        # Only map CNV names that aren't already matched
        cnv["cell_line"] = cnv["cell_line"].map(lambda x: cnv_to_sv.get(x, x))
        
        # Write mapping used
        mapping_used = cnv[["cell_line"]].drop_duplicates()
        mapping_used.to_csv("logs/name_mapping_used.csv", index=False)
        mapping_stats = {
            "direction": "cnv->sv",
            "n_prefix_mapped": int((cnv["cell_line"].isin(set(sv["cell_line"].unique()))).sum())
        }
    else:
        # Map SV (CCLE) → ACH using sample_info, if possible
        if ccle_to_ach_dict:
            sv["cell_line"] = sv["cell_line"].map(ccle_to_ach_dict).fillna(sv["cell_line"])
            # Normalize ACH IDs
            sv["cell_line"] = sv["cell_line"].apply(normalize_ach)
        
        # Step 2: prefix refinement (now in the ACH namespace)
        sv_to_cnv, cnv_to_sv = build_prefix_mapping(
            set(sv["cell_line"].unique()),
            set(cnv["cell_line"].unique())
        )
        sv["cell_line"] = sv["cell_line"].map(lambda x: sv_to_cnv.get(x, x))
        mapping_stats = {
            "direction": "sv->cnv",
            "n_prefix_mapped": int((sv["cell_line"].isin(set(cnv["cell_line"].unique()))).sum())
        }
    
    return sv, cnv, mapping_stats


def load_sample_info(sample_info_path: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Load DepMap sample_info.csv and build ACH↔CCLE mappings.
    
    Returns:
        (depmap_to_ccle, ccle_to_depmap) dictionaries
    """
    try:
        df = pd.read_csv(sample_info_path, comment="#")
    except Exception as e:
        logging.error("Failed to read sample_info.csv %s: %s", sample_info_path, e)
        return {}, {}
    
    if df.empty:
        logging.warning("sample_info.csv is empty")
        return {}, {}
    
    # Normalize column names
    colmap = _canon_cols(df.columns)
    df = df.rename(columns=colmap)
    
    # Find DepMap_ID and CCLE_name columns (use canonicalized names after rename)
    depmap_col = None
    ccle_col = None
    
    # Try exact match first, then substring match
    depmap_aliases = ["depmap_id", "depmap", "model_id", "modelid", "model", "ach_id"]
    for alias in depmap_aliases:
        if alias in colmap.values():
            # Get the canonicalized column name (which is now the actual column name after rename)
            depmap_col = alias
            break
    # If no exact match, try substring (e.g., "modelid" contains "model")
    if depmap_col is None:
        for canon_col in df.columns:
            if any(alias in canon_col or canon_col in alias for alias in depmap_aliases):
                depmap_col = canon_col
                break
    
    ccle_aliases = ["ccle_name", "ccle", "ccleid", "stripped_cell_line_name", "strippedcelllinename"]
    for alias in ccle_aliases:
        if alias in colmap.values():
            ccle_col = alias
            break
    # If no exact match, try substring
    if ccle_col is None:
        for canon_col in df.columns:
            if any(alias in canon_col or canon_col in alias for alias in ccle_aliases):
                ccle_col = canon_col
                break
    
    if depmap_col is None or ccle_col is None:
        logging.warning(
            "Could not find DepMap_ID and/or CCLE_name in sample_info.csv. "
            f"Columns: {list(df.columns)}"
        )
        return {}, {}
    
    # Build mappings with proper normalization
    depmap_vals = df[depmap_col].astype(str).apply(normalize_ach)  # ACH IDs
    ccle_vals = df[ccle_col].astype(str).apply(normalize_ccle)  # CCLE names
    
    # Drop rows with missing values
    valid_mask = depmap_vals.notna() & ccle_vals.notna()
    depmap_to_ccle = dict(zip(depmap_vals[valid_mask], ccle_vals[valid_mask]))
    ccle_to_depmap = dict(zip(ccle_vals[valid_mask], depmap_vals[valid_mask]))
    
    logging.info(
        f"Loaded sample_info: {len(depmap_to_ccle)} DepMap_ID↔CCLE_name mappings"
    )
    
    return depmap_to_ccle, ccle_to_depmap


def load_cnv(cnv_bed: str) -> Tuple[pd.DataFrame, int, int]:
    """
    Load CNV segments with autodetection of sample and CN columns.
    
    Returns:
        (cnv_df, sample_col_idx, cn_col_idx)
    """
    # Read as TSV with no header assumption
    cnv = pd.read_csv(cnv_bed, sep="\t", header=0)
    
    # Autodetect sample column
    sample_col_idx = autodetect_sample_column(cnv)
    if sample_col_idx is None:
        raise ValueError(
            f"Could not autodetect sample column in CNV file {cnv_bed}. "
            f"Expected ACH-000000 pattern or CCLE-like pattern (e.g., A549_LUNG). "
            f"Columns: {list(cnv.columns)}"
        )
    
    sample_col = cnv.columns[sample_col_idx]
    # CNV usually starts as ACH → normalize ACH only first
    cnv["cell_line"] = cnv[sample_col].astype(str).apply(normalize_ach)
    
    # Autodetect CN column
    try:
        cn_col_idx, cn_mode = autodetect_cn_column(cnv)
    except ValueError:
        raise ValueError(f"Could not autodetect CN column in {cnv_bed}")
    
    cn_col = cnv.columns[cn_col_idx]
    if cn_mode == "segmean":
        # Convert Segment_Mean to approximate CN: cn = 2.0 * 2^segmean
        segmean = pd.to_numeric(cnv[cn_col], errors='coerce')
        cnv["cn"] = 2.0 * np.power(2.0, segmean)
        logging.info(f"CN column detected as Segment_Mean (mode=segmean), converting to absolute CN")
    else:
        cnv["cn"] = pd.to_numeric(cnv[cn_col], errors='coerce')
    
    # Autodetect chrom, start, end columns
    chrom_col = None
    start_col = None
    end_col = None
    
    for col in cnv.columns:
        col_lower = str(col).lower()
        if chrom_col is None and ("chrom" in col_lower or "chr" in col_lower):
            chrom_col = col
        if start_col is None and "start" in col_lower:
            start_col = col
        if end_col is None and "end" in col_lower:
            end_col = col
    
    if chrom_col is None or start_col is None or end_col is None:
        raise ValueError(
            f"Could not autodetect chrom/start/end columns in {cnv_bed}. "
            f"Found: chrom={chrom_col}, start={start_col}, end={end_col}"
        )
    
    cnv["chrom"] = cnv[chrom_col].astype(str).map(lambda c: normalize_chrom(c))
    cnv["start"] = pd.to_numeric(cnv[start_col], errors='coerce')
    cnv["end"] = pd.to_numeric(cnv[end_col], errors='coerce')
    
    # Store CN column name for logging (before filtering)
    cn_col_name = "computed"
    if cn_col_idx is not None:
        cn_col_name = cnv.columns[cn_col_idx] if cn_col_idx < len(cnv.columns) else "computed"
    
    # Keep only required columns
    cnv = cnv[["cell_line", "chrom", "start", "end", "cn"]].copy()
    
    logging.info(
        f"Autodetected CNV columns: sample={sample_col} (idx {sample_col_idx}), "
        f"CN={cn_col_name} (idx {cn_col_idx})"
    )
    
    return cnv, sample_col_idx, cn_col_idx


def explode_bedpe_to_breakpoints(sv_df: pd.DataFrame, cell_col: str = "cell_line") -> pd.DataFrame:
    """
    Explode BEDPE to 1D breakpoints (one row per end).
    Expect columns: chrom1,start1,end1,chrom2,start2,end2, and cell_col.
    Returns: cell_line, chrom, pos
    """
    bp1 = sv_df[[cell_col, "chrom1", "start1", "end1"]].copy()
    bp1["chrom"] = bp1["chrom1"].map(lambda c: normalize_chrom(c))
    bp1["pos"] = ((bp1["start1"].astype(np.int64) + bp1["end1"].astype(np.int64)) // 2).astype(np.int64)
    bp1 = bp1[[cell_col, "chrom", "pos"]]
    
    bp2 = sv_df[[cell_col, "chrom2", "start2", "end2"]].copy()
    bp2["chrom"] = bp2["chrom2"].map(lambda c: normalize_chrom(c))
    bp2["pos"] = ((bp2["start2"].astype(np.int64) + bp2["end2"].astype(np.int64)) // 2).astype(np.int64)
    bp2 = bp2[[cell_col, "chrom", "pos"]]
    
    bps = pd.concat([bp1, bp2], ignore_index=True)
    bps = bps[bps["chrom"] != ""]
    return bps


def dedupe_breakpoints(bps: pd.DataFrame, tol: int = 1) -> pd.DataFrame:
    """
    Deduplicate nearby breakpoints within tolerance (bp) per (cell_line, chrom).
    Returns one row per "cluster" within 'tol' bp.
    """
    out = []
    for (cell, chrom), g in bps.sort_values("pos").groupby(["cell_line", "chrom"], sort=False):
        pos = g["pos"].to_numpy()
        if pos.size == 0:
            continue
        keep = np.r_[True, np.diff(pos) > tol]
        out.append(g.loc[g.index[keep]])
    if not out:
        return bps.iloc[0:0]
    return pd.concat(out, ignore_index=True)


def nearest_dist_vectorized(bps_sorted: np.ndarray, mids: np.ndarray, cap: int) -> np.ndarray:
    """
    Vectorized nearest distance computation.
    Returns distances from mids to nearest breakpoint in bps_sorted, capped at cap.
    """
    if bps_sorted.size == 0 or mids.size == 0:
        return np.full(mids.shape, cap, dtype=np.int64)
    
    idx = np.searchsorted(bps_sorted, mids)
    left_idx = np.clip(idx - 1, 0, bps_sorted.size - 1)
    right_idx = np.clip(idx, 0, bps_sorted.size - 1)
    
    left = np.abs(mids - bps_sorted[left_idx]).astype(np.float64)
    right = np.abs(bps_sorted[right_idx] - mids).astype(np.float64)
    
    left[idx == 0] = np.inf
    right[idx == bps_sorted.size] = np.inf
    
    dist = np.minimum(left, right)
    return np.minimum(dist, cap).astype(np.int64)


def compute_nearest_breakpoint_distance(
    sv: pd.DataFrame,
    genes: pd.DataFrame,
    cap_bp: int,
    dedupe_tol: int = 1
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute nearest breakpoint distance per (gene, cell_line, chrom) using vectorized search.
    Returns (bp_dist_df, deduped_bps_df) where:
    - bp_dist_df: columns ['gene','cell_line','chrom','bp_dist']
    - deduped_bps_df: columns ['cell_line','chrom','pos'] (for QC tables)
    """
    # Explode & normalize
    bps = explode_bedpe_to_breakpoints(sv, cell_col="cell_line")
    n_before = len(bps)
    bps_deduped = dedupe_breakpoints(bps, tol=dedupe_tol)
    logging.info("Breakpoints: %d → %d after dedup (tol=%d bp)", n_before, len(bps_deduped), dedupe_tol)
    
    # Precompute gene midpoints per chrom
    genes = genes.copy()
    genes["mid"] = ((genes["start"].astype(np.int64) + genes["end"].astype(np.int64)) // 2).astype(np.int64)
    
    rows = []
    # Cache mids by chrom
    mids_by_chrom = {c: g["mid"].to_numpy() for c, g in genes.groupby("chrom", sort=False)}
    gene_by_chrom = {c: g["gene"].to_numpy() for c, g in genes.groupby("chrom", sort=False)}
    
    # Sort bps per (cell, chrom)
    for (cell, chrom), g in bps_deduped.groupby(["cell_line", "chrom"], sort=False):
        if chrom not in mids_by_chrom:
            continue
        mids = mids_by_chrom[chrom]
        genes_chr = gene_by_chrom[chrom]
        bps_sorted = np.sort(g["pos"].to_numpy())
        d = nearest_dist_vectorized(bps_sorted, mids, cap=cap_bp)
        # Append rows
        rows.append(pd.DataFrame({
            "gene": genes_chr,
            "cell_line": cell,
            "chrom": chrom,
            "bp_dist": d
        }))
    
    if not rows:
        return pd.DataFrame(columns=["gene", "cell_line", "chrom", "bp_dist"]), bps_deduped
    out = pd.concat(rows, ignore_index=True)
    return out, bps_deduped


def save_design_matrix(df: pd.DataFrame, out_path: str, also_csv: bool = False) -> None:
    """
    Save df to Parquet if pyarrow/fastparquet available; otherwise to CSV.
    Log exactly what was written and where. Never claim Parquet on fallback.
    Also write a CSV sidecar when Parquet succeeds, only if also_csv=True.
    """
    out_path_obj = Path(out_path)
    os.makedirs(out_path_obj.parent or Path("."), exist_ok=True)
    
    # Try Parquet first (pyarrow, then fastparquet)
    if out_path.endswith(".parquet") or out_path.endswith(".parq"):
        parquet_success = False
        try:
            import pyarrow
            df.to_parquet(out_path, index=False, engine="pyarrow")
            logging.info("Saved Parquet to %s", out_path)
            parquet_success = True
        except (ImportError, Exception):
            try:
                import fastparquet
                df.to_parquet(out_path, index=False, engine="fastparquet")
                logging.info("Saved Parquet to %s", out_path)
                parquet_success = True
            except (ImportError, Exception):
                pass
        
        if not parquet_success:
            # Fallback to CSV
            csv_path = str(out_path_obj.with_suffix(".csv"))
            df.to_csv(csv_path, index=False)
            logging.warning("Parquet engine not available; saved CSV to %s", csv_path)
        elif also_csv:
            # Write CSV sidecar
            csv_path = str(out_path_obj.with_suffix(".csv"))
            df.to_csv(csv_path, index=False)
            logging.info("Also wrote CSV sidecar to %s", csv_path)
    else:
        # CSV output
        df.to_csv(out_path, index=False)
        logging.info("Saved CSV to %s", out_path)


def compute_proximity_kernels(bp: pd.Series) -> pd.DataFrame:
    """
    Build proximity features from base-pair distances.
    - Exponential RBFs at 50k/100k/250k/500k
    - Inverse distance (1 / (1 + bp_dist))
    - Boolean windows (≤100k/≤250k/≤500k)
    """
    x = bp.to_numpy().astype(float)
    eps = 1e-12
    k = pd.DataFrame({
        "prox_exp_50k": np.exp(-x / 50_000.0),
        "prox_exp_100k": np.exp(-x / 100_000.0),
        "prox_exp_250k": np.exp(-x / 250_000.0),
        "prox_exp_500k": np.exp(-x / 500_000.0),
        "prox_inv": 1.0 / (x + 1.0)
    }, index=bp.index)
    # Boolean windows (interpretable)
    k["prox_any_100k"] = (x <= 100_000).astype(int)
    k["prox_any_250k"] = (x <= 250_000).astype(int)
    k["prox_any_500k"] = (x <= 500_000).astype(int)
    return k


def main():
    ap = argparse.ArgumentParser(
        description="Ingest WGS SV breakpoints and compute proximity kernels"
    )
    ap.add_argument("--sv-dir", required=True,
                    help="Directory with WGS SV BEDPE files or single file path")
    ap.add_argument("--genes", required=True,
                    help="Gene BED file: chrom,start,end,gene,strand")
    ap.add_argument("--cnv", required=True,
                    help="CNV segments BED: cell_line,chrom,start,end,cn")
    ap.add_argument("--out", required=True,
                    help="Output design matrix path (CSV or Parquet)")
    ap.add_argument("--svtypes", nargs="+", default=["TRA", "INV"],
                    help="SV types to include (default: TRA INV)")
    ap.add_argument("--cap-bp", type=int, default=5_000_000,
                    help="Cap breakpoint distance at this value (default: 5Mb)")
    ap.add_argument("--dedupe-bp-tol", type=int, default=1,
                    help="Tolerance (bp) for de-duplicating nearby breakpoints (default: 1)")
    ap.add_argument("--chrom-map", choices=["auto", "hg19", "hg38"], default="auto",
                    help="Chromosome normalization scheme (default: auto)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Load data, harmonize names, report stats, then exit before heavy compute")
    ap.add_argument("--also-csv", action="store_true",
                    help="When saving Parquet, also write CSV sidecar")
    ap.add_argument("--nonzero-eps", type=float, default=1e-8,
                    help="Epsilon threshold for 'non-trivial' kernel values (default: 1e-8)")
    save_qc_group = ap.add_mutually_exclusive_group()
    save_qc_group.add_argument("--save-qc-tables", action="store_true", default=True,
                               help="Save QC tables (gene proximity rates, cell breakpoint counts) [default]")
    save_qc_group.add_argument("--no-save-qc-tables", action="store_false", dest="save_qc_tables",
                               help="Skip QC table generation")
    ap.add_argument("--join", choices=["inner", "left", "outer"], default="inner",
                    help="Join type for merging CN and SV data (default: inner)")
    drop_group = ap.add_mutually_exclusive_group()
    drop_group.add_argument("--drop-no-sv-lines", action="store_true",
                           help="Drop rows without SV coverage (bp_dist NA) [default]")
    drop_group.add_argument("--keep-no-sv-lines", action="store_true",
                           help="Keep rows without SV coverage")
    ap.add_argument(
        "--sample-map",
        type=str,
        default=None,
        help="CSV mapping from raw SV sample names to DepMap cell_line (e.g., 'CCLE_name,DepMap_ID'). "
             "Flexible headers allowed; '#' comment lines ignored."
    )
    ap.add_argument("--force-cell-from-filename", action="store_true",
                    help="Force cell_line from filename instead of sample column")
    ap.add_argument("--sample-info", type=str, default=None,
                    help="Optional DepMap sample_info.csv for ACH↔CCLE mapping")
    ap.add_argument("--map-cnv-to", choices=["sv", "cnv"], default="sv",
                    help="Map CNV names to match SV (sv) or SV names to match CNV (cnv). Default: sv")
    ap.add_argument("-v", "--verbose", action="count", default=1,
                    help="Increase verbosity")
    args = ap.parse_args()
    
    # Handle mutually exclusive drop flags (default: drop)
    if args.keep_no_sv_lines:
        args.drop_no_sv_lines = False
    else:
        args.drop_no_sv_lines = True  # Default
    
    setup_logging(args.verbose)
    
    # Create logs directory
    Path("logs").mkdir(exist_ok=True)
    
    # Load sample map if provided
    sample_map = load_sample_map(args.sample_map)
    if sample_map:
        logging.info(f"Loaded {len(sample_map)} sample name mappings")
    
    # Load sample_info if provided
    sample_info_df = None
    if args.sample_info:
        try:
            sample_info_df = pd.read_csv(args.sample_info, comment="#")
            logging.info(f"Loaded sample_info from {args.sample_info}")
        except Exception as e:
            logging.error(f"Failed to read sample_info.csv {args.sample_info}: {e}")
            sample_info_df = None
    
    # Load data
    logging.info("Loading WGS SV data...")
    sv = load_wgs_sv(
        args.sv_dir, 
        svtypes=args.svtypes,
        sample_map=sample_map,
        force_cell_from_filename=args.force_cell_from_filename
    )
    
    logging.info("Loading gene annotations...")
    genes = load_genes_bed(args.genes)
    
    logging.info("Loading CNV segments...")
    cnv, sample_col_idx, cn_col_idx = load_cnv(args.cnv)
    
    # Apply name mapping with chaining (sample_info then prefix)
    sv, cnv, mapping_stats = apply_name_mapping(sv, cnv, args.map_cnv_to, sample_info_df)
    logging.info(
        f"Name mapping: direction={mapping_stats['direction']}, "
        f"prefix_mapped={mapping_stats['n_prefix_mapped']}"
    )
    
    # Early overlap check (before heavy computation)
    sv_cells = set(sv["cell_line"].unique())
    cnv_cells = set(cnv["cell_line"].unique())
    overlap = sv_cells & cnv_cells
    
    # Save audit artifacts
    pd.Series(sorted(sv_cells)).to_csv("logs/sv_cells.txt", index=False, header=False)
    pd.Series(sorted(cnv_cells)).to_csv("logs/cnv_cells.txt", index=False, header=False)
    pd.Series(sorted(overlap)).to_csv("logs/overlap_cells.txt", index=False, header=False)
    
    logging.info(f"CNV table covers {len(cnv_cells)} cell lines")
    logging.info(
        "Pre-merge overlap check — SV samples: %d | CNV lines: %d | Overlap: %d",
        len(sv_cells), len(cnv_cells), len(overlap)
    )
    
    # --dry-run: exit after reporting counts
    if args.dry_run:
        n_genes = genes["gene"].nunique()
        n_sv_cells = len(sv_cells)
        n_cnv_cells = len(cnv_cells)
        n_overlap = len(overlap)
        logging.info("=== DRY RUN ===")
        logging.info("Genes=%d | SV cells=%d | CNV cells=%d | Overlap=%d",
                     n_genes, n_sv_cells, n_cnv_cells, n_overlap)
        logging.info("Would compute ~ %d gene–cell pairs (upper bound).", n_genes * max(n_overlap, 1))
        sys.exit(0)
    
    if args.join == "inner" and len(overlap) == 0:
        # Write debug report with 30 examples from each side
        os.makedirs("logs", exist_ok=True)
        sv_examples = sorted(list(sv_cells))[:30]
        cnv_examples = sorted(list(cnv_cells))[:30]
        max_len = max(len(sv_examples), len(cnv_examples))
        dbg = pd.DataFrame({
            "sv_example": sv_examples + [""] * (max_len - len(sv_examples)),
            "cnv_example": cnv_examples + [""] * (max_len - len(cnv_examples))
        })
        dbg.to_csv("logs/overlap_debug_samples.csv", index=False)
        if args.sample_map:
            logging.warning(
                "Sample map provided: %s. If overlap=0, verify headers (e.g., CCLE_name,DepMap_ID) "
                "and that the CSV has data rows.",
                args.sample_map
            )
        if args.sample_info:
            logging.warning(
                "Sample info provided: %s. If overlap=0, verify --map-cnv-to setting and that "
                "sample_info.csv contains matching DepMap_ID and CCLE_name columns.",
                args.sample_info
            )
        logging.error(
            "No cell-line overlap between SV and CNV after canonicalization. "
            "Wrote logs/overlap_debug_samples.csv with 30 examples from each side."
        )
        raise SystemExit(2)
    
    # Compute CN for genes (length-weighted overlap) - preserve chrom
    logging.info("Computing gene-level CN...")
    def length_weighted_gene_cn(cnv_df, genes_df):
        """Compute length-weighted mean CN for each (gene, cell_line, chrom) triple."""
        result_rows = []
        cnv_by_cell = {k: v.sort_values(["chrom", "start"]) for k, v in cnv_df.groupby("cell_line")}
        genes_by_chrom = {k: v.sort_values(["start"]) for k, v in genes_df.groupby("chrom")}
        
        for chrom in sorted(genes_by_chrom.keys()):
            chrom_norm = normalize_chrom(str(chrom))  # Ensure normalized
            gsub = genes_by_chrom[chrom]
            for _, grow in gsub.iterrows():
                gstart, gend, gene = int(grow["start"]), int(grow["end"]), grow["gene"]
                glen = max(1, gend - gstart)
                for cell_line, csub in cnv_by_cell.items():
                    cchr = csub[csub["chrom"] == chrom_norm]
                    if cchr.empty:
                        continue
                    mask = (cchr["end"] > gstart) & (cchr["start"] < gend)
                    ov = cchr.loc[mask, ["start", "end", "cn"]]
                    if ov.empty:
                        continue
                    ov_len = (np.minimum(ov["end"].values, gend) - np.maximum(ov["start"].values, gstart)).clip(min=0)
                    weights = ov_len / glen
                    lw_cn = float(np.sum(weights * ov["cn"].values))
                    result_rows.append([gene, cell_line, chrom_norm, lw_cn])
        
        if not result_rows:
            return pd.DataFrame(columns=["gene", "cell_line", "chrom", "cn"])
        return pd.DataFrame(result_rows, columns=["gene", "cell_line", "chrom", "cn"])
    
    gene_cn = length_weighted_gene_cn(cnv, genes)
    
    # Compute breakpoint distances (returns gene, cell_line, chrom, bp_dist) and deduped breakpoints
    logging.info("Computing nearest breakpoint distances...")
    bp_dist, bps_deduped = compute_nearest_breakpoint_distance(sv, genes, cap_bp=args.cap_bp, dedupe_tol=args.dedupe_bp_tol)
    
    # Merge CN and breakpoint distances on (gene, cell_line, chrom)
    design = gene_cn.merge(bp_dist, on=["gene", "cell_line", "chrom"], how=args.join)
    
    # Diagnostics before drop
    logging.info("Design (pre-drop): rows=%d genes=%d cells=%d",
                 len(design), design['gene'].nunique(), design['cell_line'].nunique())
    n_with = design["bp_dist"].notna().sum()
    n_without = design["bp_dist"].isna().sum()
    logging.info("Rows with SV=%d, without SV=%d", n_with, n_without)
    
    # Save genes missing SV (before drop)
    if len(design) > 0:
        genes_no_sv = design.groupby("gene")["bp_dist"].apply(lambda x: x.isna().mean()).sort_values(ascending=False)
        genes_no_sv.head(100).to_csv("logs/genes_missing_sv.txt", header=["frac_missing"])
    
    # Optional post-filter: drop rows without SV coverage
    if args.drop_no_sv_lines:
        design = design[design["bp_dist"].notna()].copy()
        logging.info("After --drop-no-sv-lines: rows=%d", len(design))
    
    # Add has_same_chr_bp feature (before kernels)
    design["has_same_chr_bp"] = (design["bp_dist"] < args.cap_bp).astype(int)
    num_true = design["has_same_chr_bp"].sum()
    total = len(design)
    pct = float(num_true / total * 100) if total > 0 else 0.0
    logging.info("same-chr breakpoint present: %.1f%% (%d/%d)", pct, num_true, total)
    
    # Compute proximity kernels
    logging.info("Computing proximity kernels...")
    kernels = compute_proximity_kernels(design["bp_dist"])
    design = pd.concat([design, kernels], axis=1)
    
    # Ensure bp_dist is present (fill NaN with cap)
    design["bp_dist"] = design["bp_dist"].fillna(args.cap_bp)
    
    # Kernel statistics (with both >0 and >eps metrics)
    logging.info("=== Proximity Kernel Statistics ===")
    eps = args.nonzero_eps
    for col in kernels.columns:
        nz_gt0 = float((kernels[col] > 0).mean() * 100)
        nz_gt_eps = float((kernels[col].abs() > eps).mean() * 100)
        logging.info("%s: mean=%.6f std=%.6f >0=%.1f%% >%.0e=%.1f%%",
                     col, float(kernels[col].mean()), float(kernels[col].std()), nz_gt0, eps, nz_gt_eps)
    
    cap = args.cap_bp
    pct_capped = float((design["bp_dist"] >= cap).mean() * 100)
    logging.info("bp_dist capped at %d for %.1f%% of rows", cap, pct_capped)
    
    # CN-proximity correlations (quick confounding check)
    if "cn" in design.columns and len(design) > 100:
        cols = ["cn"] + list(kernels.columns)
        corr = design[cols].corr().loc["cn"]
        logging.info("=== CN–Proximity Correlations ===")
        max_abs = 0.0
        for c in kernels.columns:
            r = float(corr[c])
            logging.info("%s: r=%.3f", c, r)
            max_abs = max(max_abs, abs(r))
        if max_abs > 0.4:
            logging.warning("HIGH CN–PROX CORRELATION (%.3f) → potential entanglement", max_abs)
        else:
            logging.info("✓ CN–proximity correlations are low (max=%.3f)", max_abs)
    
    # Save QC tables if requested
    if args.save_qc_tables:
        # A. Per-gene proximity rates
        g = design.groupby("gene")
        df_rates = pd.DataFrame({
            "n_cells": g["cell_line"].nunique(),
            "within_100k": g["prox_any_100k"].sum(),
            "within_250k": g["prox_any_250k"].sum(),
            "within_500k": g["prox_any_500k"].sum(),
        })
        for w in (100, 250, 500):
            df_rates[f"frac_{w}k"] = df_rates[f"within_{w}k"] / df_rates["n_cells"]
        df_rates.to_csv("logs/gene_proximity_rates.csv")
        logging.info("Saved gene proximity rates to logs/gene_proximity_rates.csv")
        
        # B. Per-cell breakpoint density (after dedupe)
        cell_counts = bps_deduped.groupby("cell_line").size().reset_index(name="n_breakpoints_total")
        cell_counts = cell_counts.sort_values("n_breakpoints_total", ascending=False)
        cell_counts.to_csv("logs/cell_breakpoint_counts.csv", index=False)
        logging.info("Saved cell breakpoint counts to logs/cell_breakpoint_counts.csv")
        if len(cell_counts) > 0:
            top5 = cell_counts.head(5)["n_breakpoints_total"].tolist()
            bottom5 = cell_counts.tail(5)["n_breakpoints_total"].tolist()
            logging.info("Top 5 breakpoint counts: %s", top5)
            logging.info("Bottom 5 breakpoint counts: %s", bottom5)
    
    # Sort and save
    design = design.sort_values(["gene", "cell_line", "chrom"]).reset_index(drop=True)
    
    save_design_matrix(design, args.out, also_csv=args.also_csv)
    logging.info("Design matrix: %d rows, %d columns", len(design), len(design.columns))
    logging.info("Columns: %s", list(design.columns))


if __name__ == "__main__":
    main()

