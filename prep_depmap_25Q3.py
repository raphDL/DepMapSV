#!/usr/bin/env python3
import argparse
import logging
import os
import sys
import time
from contextlib import contextmanager
from typing import Optional, List

import numpy as np
import pandas as pd

try:
    # Use tqdm if available; otherwise no-op
    from tqdm import tqdm as _tqdm
except Exception:  # noqa: BLE001
    def _tqdm(x, **kwargs):  # type: ignore
        return x


def setup_logging(verbosity: int = 1, log_file: str | None = None):
    """
    verbosity: 0=WARNING, 1=INFO, 2+=DEBUG
    """
    level = logging.WARNING if verbosity <= 0 else (logging.INFO if verbosity == 1 else logging.DEBUG)
    handlers = [logging.StreamHandler(sys.stderr)]
    if log_file:
        handlers.append(logging.FileHandler(log_file, mode="w", encoding="utf-8"))
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S",
        handlers=handlers,
        force=True,
    )
    logging.debug("Logger initialized (level=%s, file=%s)", logging.getLevelName(level), log_file)


@contextmanager
def step(name: str):
    t0 = time.perf_counter()
    logging.info("▶ %s ...", name)
    try:
        yield
    finally:
        dt = time.perf_counter() - t0
        logging.info("✓ %s (%.2fs)", name, dt)


def maybe_tqdm(iterable, enable: bool, **kwargs):
    return _tqdm(iterable, **kwargs) if enable else iterable


def norm_chr(s: str) -> str:
    """Normalize chromosome names: ensure consistent 'chr' prefix, handle MT/M"""
    s = str(s).strip()
    # Handle MT -> M normalization
    if s.upper() in ("MT", "CHRMT"):
        s = "M"
    elif s.upper().startswith("CHR") and s.upper()[3:] == "MT":
        s = "M"
    
    if s.startswith("chr"):
        return s
    # Remove any existing chr prefix first, then add it
    s_clean = s.replace("chr", "").replace("Chr", "").replace("CHR", "")
    return "chr" + s_clean


def read_csv_fast(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, engine="pyarrow")
    except Exception:
        return pd.read_csv(path)


def pick_first_existing(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def compute_absolute_cn(segments: pd.DataFrame, assume_diploid: bool, global_signatures: Optional[pd.DataFrame], cn_is_log2: bool = False) -> pd.Series:
    # Default ploidy per model
    if assume_diploid or global_signatures is None:
        model_to_ploidy = {}
    else:
        # Try to discover ploidy column
        model_col = pick_first_existing(global_signatures, ["ModelID", "DepMapModelIdentifier", "model_id", "cell_line"])
        ploidy_col = pick_first_existing(global_signatures, ["Ploidy", "PLOIDY", "ploidy"])  # heuristic
        if model_col is None or ploidy_col is None:
            model_to_ploidy = {}
        else:
            model_to_ploidy = dict(zip(global_signatures[model_col].astype(str), global_signatures[ploidy_col].astype(float)))

    scn_col = pick_first_existing(segments, ["SEGMENT_COPY_NUMBER", "SegmentCopyNumber", "segment_copy_number", "Log2CopyNumber", "CopyNumber", "copy_number"])
    if scn_col is None:
        raise ValueError("CN segments file must contain a CN column (SEGMENT_COPY_NUMBER, CopyNumber, etc.)")

    model_col = pick_first_existing(segments, ["ModelID", "DepMapModelIdentifier", "model_id", "cell_line"])
    if model_col is None:
        raise ValueError("CN segments file must contain a ModelID-like column")

    # Guardrail: detect if data looks like absolute CN when --cn_is_log2 is set
    scn_series = pd.to_numeric(segments[scn_col], errors="coerce")
    p95 = float(scn_series.quantile(0.95)) if len(scn_series) > 0 else float("nan")
    p99 = float(scn_series.quantile(0.99)) if len(scn_series) > 0 else float("nan")
    mx = float(scn_series.max()) if len(scn_series) > 0 else float("nan")
    md = float(scn_series.median()) if len(scn_series) > 0 else float("nan")
    mn = float(scn_series.min()) if len(scn_series) > 0 else float("nan")
    logging.debug("SEGMENT_COPY_NUMBER stats: min=%.3f med=%.3f p95=%.3f p99=%.3f max=%.3f", mn, md, p95, p99, mx)
    
    # Override cn_is_log2 if data appears to be absolute CN
    cn_is_log2_effective = cn_is_log2
    if cn_is_log2 and (p99 > 20 or mx > 50):
        logging.warning(
            "Input CN appears to be ABSOLUTE (p99=%.2f, max=%.2f). Overriding --cn_is_log2 → False to avoid exponentiating absolute CN.",
            p99, mx
        )
        cn_is_log2_effective = False

    # cn = ploidy * 2 ** scn if log2, else already absolute
    def to_abs(row):
        ploidy = model_to_ploidy.get(str(row[model_col]), 2.0)
        scn = float(row[scn_col])
        if cn_is_log2_effective:
            return float(ploidy) * (2.0 ** scn)
        else:
            return scn  # already absolute

    abs_cn = segments.apply(to_abs, axis=1)
    logging.info("CN computed: min=%.2f, mean=%.2f, max=%.2f", float(np.nanmin(abs_cn)), float(np.nanmean(abs_cn)), float(np.nanmax(abs_cn)))
    return abs_cn


def write_cnv_bed(segments: pd.DataFrame, out_path: str):
    contig_col = pick_first_existing(segments, ["CONTIG", "Chromosome", "chrom", "chr"]) or "CONTIG"
    start_col = pick_first_existing(segments, ["START", "start", "Start" ]) or "START"
    end_col = pick_first_existing(segments, ["END", "end", "End"]) or "END"
    model_col = pick_first_existing(segments, ["ModelID", "DepMapModelIdentifier", "model_id", "cell_line"]) or "ModelID"

    df = segments[[model_col, contig_col, start_col, end_col, "cn"]].copy()
    df.columns = ["cell_line", "chrom", "start", "end", "cn"]
    # Normalize chromosome names
    df["chrom"] = df["chrom"].map(norm_chr)
    # Ensure integer coordinates
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df.to_csv(out_path, sep="\t", index=False)


def write_sv_from_cnv_boundaries(segments: pd.DataFrame, out_path: str, show_progress: bool = False):
    contig_col = pick_first_existing(segments, ["CONTIG", "Chromosome", "chrom", "chr"]) or "CONTIG"
    start_col = pick_first_existing(segments, ["START", "start", "Start" ]) or "START"
    end_col = pick_first_existing(segments, ["END", "end", "End"]) or "END"
    model_col = pick_first_existing(segments, ["ModelID", "DepMapModelIdentifier", "model_id", "cell_line"]) or "ModelID"

    # For each segment, generate two boundary breakpoints (start and end) as 1bp intervals
    rows = []
    for model_id, sub in maybe_tqdm(segments.groupby(model_col), show_progress, desc="Models"):
        for _, r in sub.iterrows():
            chrom = norm_chr(str(r[contig_col]))
            s = int(r[start_col])
            e = int(r[end_col])
            # start boundary
            s1 = max(0, s - 1)
            rows.append([chrom, s1, s, chrom, s1, s, str(model_id), "BOUNDARY"])
            # end boundary
            e1 = max(0, e - 1)
            rows.append([chrom, e1, e, chrom, e1, e, str(model_id), "BOUNDARY"])

    out_df = pd.DataFrame(rows, columns=["chrom1","start1","end1","chrom2","start2","end2","cell_line","svtype"])
    out_df.to_csv(out_path, sep="\t", index=False)
    logging.debug("Generated %d pseudo-SV breakpoints from %d segments", len(rows), len(segments))


def prepare_expression_long(expr: pd.DataFrame) -> pd.DataFrame:
    # Try to detect default-entry flags
    default_cols = [
        "is_default_entry", "IsDefaultEntryForModel", "IsDefaultEntryForMC",
        "isDefault", "is_model_default"
    ]
    default_col = pick_first_existing(expr, default_cols)

    gene_col = pick_first_existing(expr, ["HugoSymbol", "gene", "Gene", "GeneSymbol"]) or "HugoSymbol"
    model_col = pick_first_existing(expr, ["ModelID", "DepMapModelIdentifier", "model_id", "cell_line", "Model" ])
    value_col = pick_first_existing(expr, ["TPM", "TPMLogp1", "logTPM", "Expression", "expression"])  # heuristic

    if model_col is not None and value_col is not None and gene_col in expr.columns:
        df = expr[[gene_col, model_col, value_col] + ([default_col] if default_col else [])].copy()
        if default_col:
            df = df[df[default_col] == True]
        df = df[[gene_col, model_col, value_col]].copy()
        df.columns = ["gene", "cell_line", "expression"]
        return df

    # Otherwise, assume a matrix with rows as genes and columns as models
    if gene_col in expr.columns:
        expr = expr.set_index(gene_col)
    expr_long = expr.melt(ignore_index=False, var_name="cell_line", value_name="expression").reset_index()
    # After reset_index(), the gene column is named 'index' unless original index had a name
    if "gene" not in expr_long.columns:
        if "index" in expr_long.columns:
            expr_long.rename(columns={"index": "gene"}, inplace=True)
        else:
            # fallback: if the original index kept its name
            expr_long.rename(columns={gene_col: "gene"}, inplace=True)
    if default_col and default_col in expr_long.columns:
        expr_long = expr_long[expr_long[default_col] == True]
        expr_long = expr_long[["gene", "cell_line", "expression"]]
    return expr_long[["gene", "cell_line", "expression"]]


def main():
    ap = argparse.ArgumentParser(description="Prepare DepMap inputs: CNV segments (absolute), pseudo-SV from CNV boundaries, expression long format.")
    ap.add_argument("--cnv_segments_wgs", required=True, help="OmicsCNSegmentsWGS.csv")
    ap.add_argument("--expression_tpm", required=True, help="OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv")
    ap.add_argument("--global_signatures", required=False, help="OmicsGlobalSignatures.csv (optional for per-model ploidy)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--cn_assume_diploid", action="store_true", help="If set, use ploidy=2.0 regardless of global signatures")
    ap.add_argument("--cn_is_log2", action="store_true", help="If set, treat SEGMENT_COPY_NUMBER as log2 and exponentiate; otherwise assume absolute CN")
    ap.add_argument("-v", "--verbose", action="count", default=1, help="Increase verbosity (-v=INFO, -vv=DEBUG).")
    ap.add_argument("-q", "--quiet", action="store_true", help="Quiet mode (only warnings/errors).")
    ap.add_argument("--log-file", type=str, default=None, help="Optional path to write logs.")
    ap.add_argument("--progress", action="store_true", help="Show progress bars for long-running loops.")
    args = ap.parse_args()

    # Setup logging
    verbosity = args.verbose if not args.quiet else 0
    setup_logging(verbosity=verbosity, log_file=args.log_file)
    show_progress = args.progress

    # Echo arguments
    try:
        logging.info("Args: %s", vars(args))
    except Exception:
        pass

    os.makedirs(args.outdir, exist_ok=True)

    with step("Load CNV segments"):
        segments = read_csv_fast(args.cnv_segments_wgs)
        logging.info("CNV segments: %d rows, columns=%s", len(segments), list(segments.columns)[:8])
        logging.debug("Example CNV row: %s", segments.iloc[0].to_dict() if len(segments) > 0 else "N/A")

    globals_df = None
    if args.global_signatures:
        with step("Load global signatures"):
            globals_df = read_csv_fast(args.global_signatures)
            logging.info("Global signatures: %d rows", len(globals_df))

    with step("Compute absolute copy number"):
        segments = segments.copy()
        segments["cn"] = compute_absolute_cn(segments, args.cn_assume_diploid, globals_df, args.cn_is_log2)
        logging.info("CN computed: min=%.2f, max=%.2f, mean=%.2f", segments["cn"].min(), segments["cn"].max(), segments["cn"].mean())
        logging.debug("CN scaling: assume_diploid=%s, cn_is_log2=%s", args.cn_assume_diploid, args.cn_is_log2)

    with step("Write CNV BED"):
        cnv_bed_path = os.path.join(args.outdir, "cnv_segments.bed")
        write_cnv_bed(segments, cnv_bed_path)
        logging.info("Wrote: %s (%d segments)", cnv_bed_path, len(segments))

    with step("Write pseudo-SV from CNV boundaries"):
        sv_bedpe_path = os.path.join(args.outdir, "sv_from_cnv.bedpe")
        write_sv_from_cnv_boundaries(segments, sv_bedpe_path, show_progress)
        logging.info("Wrote: %s", sv_bedpe_path)

    with step("Prepare expression long format"):
        expr_df = read_csv_fast(args.expression_tpm)
        logging.info("Expression input: %d rows, columns=%s", len(expr_df), list(expr_df.columns)[:8])
        expr_long = prepare_expression_long(expr_df)
        expr_out = os.path.join(args.outdir, "expression_long.csv")
        expr_long.to_csv(expr_out, index=False)
        logging.info("Wrote: %s (%d gene-cell pairs)", expr_out, len(expr_long))

    logging.info("Prep complete: wrote %s, %s, %s", cnv_bed_path, sv_bedpe_path, expr_out)


if __name__ == "__main__":
    main()


