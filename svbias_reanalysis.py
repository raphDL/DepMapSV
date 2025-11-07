#!/usr/bin/env python3
import argparse
import json
import logging
import os
import sys
import time
from contextlib import contextmanager
from typing import Optional, Tuple, Dict, List

import numpy as np
import pandas as pd
from sklearn.linear_model import HuberRegressor, LinearRegression
from sklearn.metrics import roc_auc_score, average_precision_score

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


def read_csv_fast(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, engine="pyarrow")
    except Exception:
        return pd.read_csv(path)


def load_dependency(dep_path: str) -> pd.DataFrame:
    df = read_csv_fast(dep_path)
    cols = {c.lower() for c in df.columns}
    if {"gene", "cell_line", "dependency"}.issubset(cols):
        # Normalize column names
        colmap = {c: c.lower() for c in df.columns}
        df = df.rename(columns=colmap)
        out = df[["gene", "cell_line", "dependency"]].copy()
        out["cell_line"] = out["cell_line"].astype(str)
        out["gene"] = out["gene"].astype(str)
        # Normalize DepMap gene names like "A1BG (1)" -> "A1BG" for merging
        out["gene"] = out["gene"].str.replace(r"\s*\(\d+\)\s*$", "", regex=True).str.strip()
        logging.debug("Dependency loaded as long format")
        return out
    # Assume wide: rows=ModelID, columns=genes
    # Prefer specific ID column names
    id_candidates = ["ModelID", "DepMap_ID", "DepMapModelIdentifier", "cell_line"]
    id_col = None
    for cand in id_candidates:
        if cand in df.columns:
            id_col = cand
            break
    if id_col is None:
        # Fall back to first column
        id_col = df.columns[0]
        logging.debug("Using first column as ID: %s", id_col)
    else:
        logging.debug("Using ID column: %s", id_col)
    df[id_col] = df[id_col].astype(str)
    long = df.melt(id_vars=[id_col], var_name="gene", value_name="dependency")
    long.rename(columns={id_col: "cell_line"}, inplace=True)
    # Normalize DepMap gene names like "A1BG (1)" -> "A1BG" for merging
    long["gene"] = long["gene"].astype(str).str.replace(r"\s*\(\d+\)\s*$", "", regex=True).str.strip()
    logging.debug("Dependency loaded as wide format, melted to long")
    return long[["gene", "cell_line", "dependency"]]


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


def load_cnv(cnv_bed_path: str) -> pd.DataFrame:
    cnv = pd.read_csv(cnv_bed_path, sep="\t")
    # Expect: cell_line, chrom, start, end, cn
    required = {"cell_line", "chrom", "start", "end", "cn"}
    if not required.issubset(set(cnv.columns)):
        raise ValueError(f"CNV file must contain columns: {required}")
    cnv["cell_line"] = cnv["cell_line"].astype(str)
    cnv["chrom"] = cnv["chrom"].map(norm_chr)
    return cnv


def load_sv(sv_bedpe_path: str) -> pd.DataFrame:
    sv = pd.read_csv(sv_bedpe_path, sep="\t")
    # Expect: chrom1,start1,end1,chrom2,start2,end2,cell_line,svtype
    needed = ["chrom1","start1","end1","chrom2","start2","end2","cell_line"]
    for c in needed:
        if c not in sv.columns:
            raise ValueError("SV BEDPE missing required columns")
    sv["cell_line"] = sv["cell_line"].astype(str)
    sv["chrom1"] = sv["chrom1"].map(norm_chr)
    sv["chrom2"] = sv["chrom2"].map(norm_chr)
    return sv
# hg38 chrom sizes (autosomes)
CHR_SIZES = {
    "chr1":248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,"chr5":181538259,
    "chr6":170805979,"chr7":159345973,"chr8":145138636,"chr9":138394717,"chr10":133797422,
    "chr11":135086622,"chr12":133275309,"chr13":114364328,"chr14":107043718,"chr15":101991189,
    "chr16":90338345,"chr17":83257441,"chr18":80373285,"chr19":58617616,"chr20":64444167,
    "chr21":46709983,"chr22":50818468
}

def _shuffle_uniform_on(sv: pd.DataFrame, chrom_col: str, s_col: str, e_col: str, rng: np.random.Generator) -> pd.DataFrame:
    sv = sv.copy()
    groups = sv.groupby(["cell_line", chrom_col]).groups
    for (cell, chrom), idx in groups.items():
        L = CHR_SIZES.get(str(chrom), None)
        if L is None or len(idx) == 0:
            continue
        n = len(idx)
        seglen = (sv.loc[idx, e_col].astype(int) - sv.loc[idx, s_col].astype(int)).clip(lower=1).to_numpy()
        mids = rng.integers(0, L, size=n)
        starts = np.maximum(0, mids - seglen//2)
        ends   = np.minimum(L, starts + seglen)
        sv.loc[idx, s_col] = np.minimum(starts, ends)
        sv.loc[idx, e_col] = np.maximum(starts, ends)
    return sv

def _shuffle_rotate_on(sv: pd.DataFrame, chrom_col: str, s_col: str, e_col: str, rng: np.random.Generator) -> pd.DataFrame:
    sv = sv.copy()
    groups = sv.groupby(["cell_line", chrom_col]).groups
    for (cell, chrom), idx in groups.items():
        L = CHR_SIZES.get(str(chrom), None)
        if L is None or len(idx) == 0:
            continue
        off = int(rng.integers(0, L))
        s = (sv.loc[idx, s_col].astype(int).to_numpy() + off) % L
        e = (sv.loc[idx, e_col].astype(int).to_numpy() + off) % L
        sv.loc[idx, s_col] = np.minimum(s, e)
        sv.loc[idx, e_col] = np.maximum(s, e)
    return sv

def apply_shuffle_modes(sv: pd.DataFrame, shuffle_uniform: bool = False, shuffle_rotate: bool = False, seed: int = 1) -> pd.DataFrame:
    if not (shuffle_uniform or shuffle_rotate):
        return sv
    rng = np.random.default_rng(seed)
    out = sv.copy()
    if shuffle_uniform:
        out = _shuffle_uniform_on(out, "chrom1", "start1", "end1", rng)
        out = _shuffle_uniform_on(out, "chrom2", "start2", "end2", rng)
    if shuffle_rotate:
        out = _shuffle_rotate_on(out, "chrom1", "start1", "end1", rng)
        out = _shuffle_rotate_on(out, "chrom2", "start2", "end2", rng)
    return out


def load_genes_bed(genes_bed: str) -> pd.DataFrame:
    g = pd.read_csv(genes_bed, sep="\t", header=None, names=["chrom","start","end","gene","strand"])
    g["chrom"] = g["chrom"].map(norm_chr)
    return g


def load_expression_long(expr_path: Optional[str]) -> Optional[pd.DataFrame]:
    if not expr_path:
        return None
    df = read_csv_fast(expr_path)
    cols = set(df.columns)
    if {"gene", "cell_line", "expression"}.issubset(cols):
        df["cell_line"] = df["cell_line"].astype(str)
        df["gene"] = df["gene"].astype(str)
        return df[["gene","cell_line","expression"]]
    raise ValueError("Expression must be long format with gene,cell_line,expression")


def length_weighted_gene_cn(cnv: pd.DataFrame, genes: pd.DataFrame, show_progress: bool = False) -> pd.DataFrame:
    # For each gene and cell_line, compute length-weighted mean CN of overlapping segments
    # naive interval overlap implementation
    result_rows = []
    cnv_by_cell = {k: v.sort_values(["chrom","start"]) for k, v in cnv.groupby("cell_line")}
    genes_by_chrom = {k: v.sort_values(["start"]) for k, v in genes.groupby("chrom")}

    chroms = sorted(genes_by_chrom.keys())
    for chrom in maybe_tqdm(chroms, show_progress, desc="Chromosomes"):
        gsub = genes_by_chrom[chrom]
        # Pre-index CNV by chrom across cell lines on demand
        for _, grow in gsub.iterrows():
            gstart, gend, gene = int(grow["start"]), int(grow["end"]), grow["gene"]
            glen = max(1, gend - gstart)
            for cell_line, csub in cnv_by_cell.items():
                cchr = csub[csub["chrom"] == chrom]
                if cchr.empty:
                    continue
                # find overlaps
                mask = (cchr["end"] > gstart) & (cchr["start"] < gend)
                ov = cchr.loc[mask, ["start","end","cn"]]
                if ov.empty:
                    continue
                ov_len = (np.minimum(ov["end"].values, gend) - np.maximum(ov["start"].values, gstart)).clip(min=0)
                weights = ov_len / glen
                lw_cn = float(np.sum(weights * ov["cn"].values))
                result_rows.append([gene, cell_line, lw_cn])
    if not result_rows:
        return pd.DataFrame(columns=["gene","cell_line","cn"])
    out = pd.DataFrame(result_rows, columns=["gene","cell_line","cn"])
    logging.debug("Computed CN for %d gene-cell pairs", len(out))
    return out


def nearest_breakpoint_distance(sv: pd.DataFrame, genes: pd.DataFrame, cap_bp: int = 2_000_000, show_progress: bool = False) -> pd.DataFrame:
    # Compute, for each gene x cell_line, distance to nearest breakpoint on same chromosome
    # Bedpe breakpoints: (chrom1,endpoints) and (chrom2,endpoints)
    # Build per cell_line per chrom breakpoint list
    logging.debug("Building breakpoint index...")
    bp_records: Dict[Tuple[str,str], np.ndarray] = {}
    for (cell, chrom), sub in pd.concat([
        sv[["cell_line","chrom1","start1","end1"]].rename(columns={"chrom1":"chrom","start1":"start","end1":"end"}),
        sv[["cell_line","chrom2","start2","end2"]].rename(columns={"chrom2":"chrom","start2":"start","end2":"end"}),
    ], ignore_index=True).groupby(["cell_line","chrom"]):
        # Use midpoints of 1bp intervals
        mids = ((sub["start"].astype(int) + sub["end"].astype(int)) // 2).values
        bp_records[(str(cell), str(chrom))] = np.unique(mids)

    logging.debug("Computing distances for %d genes...", len(genes))
    rows = []
    # Use iterrows for simpler attribute access
    genes_iter = genes.iterrows()
    if show_progress:
        genes_iter = maybe_tqdm(genes_iter, True, desc="Genes", total=len(genes))
    for _, g in genes_iter:
        chrom = str(g["chrom"])
        gstart = int(g["start"])
        gend = int(g["end"])
        gene = str(g["gene"])
        gmid = int((gstart + gend) // 2)
        for (cell, cchrom), bps in bp_records.items():
            if cchrom != chrom:
                continue
            if bps.size == 0:
                continue
            # nearest distance
            idx_search = np.searchsorted(bps, gmid)
            candidates = []
            if idx_search > 0:
                candidates.append(abs(gmid - bps[idx_search-1]))
            if idx_search < bps.size:
                candidates.append(abs(gmid - bps[idx_search]))
            dist = min(candidates) if candidates else cap_bp
            rows.append([gene, cell, int(min(dist, cap_bp))])
    if not rows:
        return pd.DataFrame(columns=["gene","cell_line","bp_dist"])
    out = pd.DataFrame(rows, columns=["gene","cell_line","bp_dist"])
    logging.debug("Computed breakpoint distances for %d gene-cell pairs", len(out))
    return out 


def fit_models_and_predict(dep_long: pd.DataFrame,
                           features: pd.DataFrame,
                           model_type: str,
                           alpha: float,
                           epsilon: float,
                           use_expr: bool,
                           bp_windows: Tuple[int,int],
                           default_ploidy: float = 2.0,
                           show_progress: bool = False,
                           add_continuous_proximity: bool = False,
                           standardize_predictors: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # Merge features
    df = dep_long.merge(features, on=["gene","cell_line"], how="left")
    # Build window indicators - use large distance for missing (far away)
    w100k, w1m = bp_windows
    df["bp_dist"] = df["bp_dist"].fillna(10_000_000)  # "far away" default
    df["bp_within_100000"] = (df["bp_dist"] <= w100k).astype(int)
    df["bp_within_1000000"] = (df["bp_dist"] <= w1m).astype(int)
    
    if add_continuous_proximity:
        df["inv_bp"] = 1.0 / (df["bp_dist"] + 1e3)
        logging.info("Added continuous proximity term (1/(bp_dist+1e3))")

    # Prepare results containers
    pred_struct = []
    model_rows = []

    genes = sorted(df["gene"].unique())
    logging.info("Fitting models for %d genes", len(genes))
    if standardize_predictors:
        logging.info("Standardizing predictors per gene")
    
    for gene in maybe_tqdm(genes, show_progress, desc="Per-gene robust fits"):
        gdf = df[df["gene"] == gene].copy()
        y = gdf["dependency"].values
        # Features: intercept + CN + windows + optional expression + optional continuous proximity
        X_cols = ["cn", "bp_within_100000", "bp_within_1000000"]
        if add_continuous_proximity:
            X_cols.append("inv_bp")
        if use_expr and "expression" in gdf.columns:
            X_cols.append("expression")
        # Fill missing CN with ploidy, not 0
        gdf_filled = gdf[X_cols].copy()
        gdf_filled["cn"] = gdf_filled["cn"].fillna(default_ploidy)
        if "expression" in gdf_filled.columns:
            gdf_filled["expression"] = gdf_filled["expression"].fillna(0.0)  # expression can be 0 if missing
        
        # Optional standardization per gene
        scaling_info = {}
        if standardize_predictors:
            for col in ["cn", "bp_within_100000", "bp_within_1000000"]:
                if col in gdf_filled.columns:
                    mean_val = gdf_filled[col].mean()
                    std_val = gdf_filled[col].std()
                    if std_val > 0:
                        gdf_filled[col] = (gdf_filled[col] - mean_val) / std_val
                        scaling_info[col] = (mean_val, std_val)
        
        X = gdf_filled.astype(float).values

        # guard against empty
        if len(gdf) < 5:
            pred = np.full(len(gdf), np.nan)
            r2 = np.nan
            coef = [np.nan] * len(X_cols)
            intercept = np.nan
        else:
            if model_type == "huber":
                model = HuberRegressor(alpha=alpha, epsilon=epsilon)
            else:
                model = LinearRegression()
            try:
                model.fit(X, y)
                # r2
                ybar = np.mean(y)
                pred_full = model.predict(X)
                ss_res = np.sum((y - pred_full) ** 2)
                ss_tot = np.sum((y - ybar) ** 2)
                r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
                coef = list(model.coef_.ravel())
                intercept = float(model.intercept_)
                
                # CRITICAL: Only subtract structural effects (CN + proximity), NOT expression
                # Build coefficient map
                coef_map = {c: float(v) for c, v in zip(X_cols, coef)}
                cn_val = gdf["cn"].astype(float).fillna(default_ploidy).values
                bp1 = gdf["bp_within_100000"].astype(float).values
                bp2 = gdf["bp_within_1000000"].astype(float).values
                # Include inv_bp if continuous proximity was used
                inv = gdf["inv_bp"].astype(float).values if add_continuous_proximity and "inv_bp" in gdf.columns else None
                
                pred = intercept \
                    + coef_map.get("cn", 0.0) * cn_val \
                    + coef_map.get("bp_within_100000", 0.0) * bp1 \
                    + coef_map.get("bp_within_1000000", 0.0) * bp2 \
                    + (coef_map.get("inv_bp", 0.0) * inv if inv is not None else 0.0)
            except Exception:
                pred = np.full(len(gdf), np.nan)
                r2 = np.nan
                coef = [np.nan] * len(X_cols)
                intercept = np.nan

        gdf["pred_structural"] = pred
        pred_struct.append(gdf[["gene","cell_line","pred_structural"]])
        model_row = {"gene": gene, "intercept": intercept, "n": len(gdf), "r2": r2}
        for c, v in zip(X_cols, coef):
            model_row[c] = v
        model_rows.append(model_row)

    pred_df = pd.concat(pred_struct, ignore_index=True) if pred_struct else pd.DataFrame(columns=["gene","cell_line","pred_structural"])
    models_df = pd.DataFrame(model_rows)
    # Combine predictions back into dependency table
    out = dep_long.merge(pred_df, on=["gene","cell_line"], how="left")
    return out, models_df


def correct_dependencies(pred_df: pd.DataFrame) -> pd.DataFrame:
    # dependency_corrected = dependency - (pred_structural - median_pred_structural_per_gene)
    med = pred_df.groupby("gene")["pred_structural"].median().rename("gene_median_pred")
    out = pred_df.merge(med, on="gene", how="left")
    out["dependency_corrected"] = out["dependency"] - (out["pred_structural"] - out["gene_median_pred"])
    return out[["gene","cell_line","dependency","dependency_corrected","pred_structural"]]


def qc_flags(features: pd.DataFrame,
             default_ploidy: float,
             cn_abs_threshold: float,
             bp_close: int,
             bp_very_close: int,
             require_two: bool) -> pd.DataFrame:
    f = features.copy()
    f["flag_cn"] = (np.abs(f["cn"] - default_ploidy) >= cn_abs_threshold).astype(int)
    f["flag_bp_close"] = (f["bp_dist"].fillna(np.inf) <= bp_close).astype(int)
    f["flag_bp_very_close"] = (f["bp_dist"].fillna(np.inf) <= bp_very_close).astype(int)
    if require_two:
        f["confounded"] = ((f[["flag_cn","flag_bp_close","flag_bp_very_close"]].sum(axis=1)) >= 2).astype(int)
    else:
        f["confounded"] = ((f[["flag_cn","flag_bp_close","flag_bp_very_close"]].sum(axis=1)) >= 1).astype(int)
    return f[["gene","cell_line","flag_cn","flag_bp_close","flag_bp_very_close","confounded"]]


def evaluate(dep: pd.DataFrame, essentials_path: Optional[str], nonessentials_path: Optional[str], corrected_col: str) -> Optional[pd.DataFrame]:
    if not essentials_path or not nonessentials_path:
        return None
    ess = set(read_csv_fast(essentials_path).iloc[:,0].astype(str))
    noness = set(read_csv_fast(nonessentials_path).iloc[:,0].astype(str))
    keep = dep[dep["gene"].isin(ess.union(noness))].copy()
    if keep.empty:
        return None
    # Drop NA before metrics
    keep = keep.dropna(subset=["dependency", corrected_col])
    if keep.empty:
        return None
    y = keep["gene"].isin(ess).astype(int).values  # 1 = essential
    # Scores: more negative = more essential; invert sign for AUROC/AUPRC to align higher = more essential
    s_orig = -keep["dependency"].values
    s_corr = -keep[corrected_col].values
    def safe_auc(y_true, scores):
        try:
            return float(roc_auc_score(y_true, scores))
        except Exception:
            return np.nan
    def safe_aupr(y_true, scores):
        try:
            return float(average_precision_score(y_true, scores))
        except Exception:
            return np.nan
    rows = [
        ["original", safe_auc(y, s_orig), safe_aupr(y, s_orig)],
        ["corrected", safe_auc(y, s_corr), safe_aupr(y, s_corr)]
    ]
    return pd.DataFrame(rows, columns=["version","auroc","auprc"])


def main():
    ap = argparse.ArgumentParser(description="Reanalysis pipeline to detect/correct SV/CN bias in DepMap CRISPR gene-effect scores")
    ap.add_argument("--dependency", required=True, help="Dependency file (wide ModelID x Gene or long gene,cell_line,dependency)")
    ap.add_argument("--cnv", required=True, help="CNV segments BED-like: cell_line,chrom,start,end,cn")
    ap.add_argument("--sv", required=True, help="SV BEDPE: chrom1,start1,end1,chrom2,start2,end2,cell_line")
    ap.add_argument("--genes", required=True, help="Gene BED: chrom,start,end,gene,strand")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--expression", required=False, help="Expression long: gene,cell_line,expression")
    ap.add_argument("--essential", required=False, help="Common essentials (one gene per line)")
    ap.add_argument("--nonessential", required=False, help="Nonessentials (one gene per line)")
    ap.add_argument("--model", choices=["huber","linear"], default="huber")
    ap.add_argument("--alpha", type=float, default=1e-4)
    ap.add_argument("--epsilon", type=float, default=1.35)
    ap.add_argument("--bp_windows", nargs=2, type=int, default=[100000, 1000000])
    ap.add_argument("--default_ploidy", type=float, default=2.0)
    ap.add_argument("--cap_cn", type=float, default=8.0)
    ap.add_argument("--qc_cn_abs_threshold", type=float, default=1.0)
    ap.add_argument("--qc_close_bp", type=int, default=100000)
    ap.add_argument("--qc_very_close_bp", type=int, default=20000)
    ap.add_argument("--qc_require_two", action="store_true")
    # Pilot mode flags
    ap.add_argument("--pilot", action="store_true", help="Enable fast pilot mode on a small gene subset")
    ap.add_argument("--pilot-size", type=int, default=50, help="Number of genes per set (unstable and optional stable)")
    ap.add_argument("--pilot-unstable-only", action="store_true", help="Only run unstable set (skip stable control)")
    ap.add_argument("--pilot-genes-file", default=None, help="Optional file with one gene per line; if set in --pilot, use exactly these genes.")
    ap.add_argument("--shuffle-sv", action="store_true", help="Negative control: shuffle breakpoint positions within each chr×cell before analysis")
    ap.add_argument("-v", "--verbose", action="count", default=1, help="Increase verbosity (-v=INFO, -vv=DEBUG).")
    ap.add_argument("-q", "--quiet", action="store_true", help="Quiet mode (only warnings/errors).")
    ap.add_argument("--log-file", type=str, default=None, help="Optional path to write logs.")
    ap.add_argument("--progress", action="store_true", help="Show progress bars for long-running loops.")
    ap.add_argument("--add-continuous-proximity", action="store_true", help="Add continuous proximity term (1/(bp_dist+1e3)) as predictor.")
    ap.add_argument("--standardize-predictors", action="store_true", help="Z-score CN and proximity features per gene before fitting.")
    ap.add_argument("--shuffle-uniform", action="store_true", help="Shuffle SV breakpoints uniformly within chromosome lengths (per cell × chrom).")
    ap.add_argument("--shuffle-rotate", action="store_true", help="Rotate SV positions by a random offset per cell × chrom (mod chromosome length).")
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

    os.makedirs(args.out, exist_ok=True)

    with step("Load dependency matrix"):
        dep = load_dependency(args.dependency)
        logging.info("Dependencies: genes=%d, cells=%d, pairs=%d", dep["gene"].nunique(), dep["cell_line"].nunique(), len(dep))

    with step("Load CNV segments"):
        cnv = load_cnv(args.cnv)
        logging.info("CNV segments: %d rows, %d cells", len(cnv), cnv["cell_line"].nunique())

    with step("Load SV breakpoints"):
        sv = load_sv(args.sv)
        # Apply shuffle modes if requested
        if args.shuffle_uniform or args.shuffle_rotate or args.shuffle_sv:
            if args.shuffle_uniform or args.shuffle_rotate:
                logging.info("Applying strong shuffle: uniform=%s rotate=%s", args.shuffle_uniform, args.shuffle_rotate)
                sv = apply_shuffle_modes(sv, shuffle_uniform=args.shuffle_uniform, shuffle_rotate=args.shuffle_rotate, seed=1)
            elif args.shuffle_sv:
                logging.info("Applying legacy shuffle (--shuffle-sv)")
                # Legacy: permute start/end within groups (cell×chrom) separately for sides 1 and 2
                def _legacy_shuffle_block(df: pd.DataFrame, start_col: str, end_col: str) -> pd.DataFrame:
                    idx = np.random.permutation(df.index)
                    df[start_col] = df.loc[idx, start_col].values
                    df[end_col] = df.loc[idx, end_col].values
                    return df
                sv = sv.groupby(["cell_line","chrom1"], group_keys=False).apply(lambda d: _legacy_shuffle_block(d.copy(), "start1","end1")).reset_index(drop=True)
                sv = sv.groupby(["cell_line","chrom2"], group_keys=False).apply(lambda d: _legacy_shuffle_block(d.copy(), "start2","end2")).reset_index(drop=True)
        logging.info("SV breakpoints: %d rows, %d cells", len(sv), sv["cell_line"].nunique())

    with step("Load gene annotations"):
        genes = load_genes_bed(args.genes)
        logging.info("Genes: %d", len(genes))

    # Guardrails: overlaps
    dep_genes = set(dep["gene"].astype(str).unique())
    genes_bed_genes = set(genes["gene"].astype(str).unique())
    gene_overlap = len(dep_genes & genes_bed_genes)
    logging.info("Gene overlap between dependency and BED: %d", gene_overlap)
    if gene_overlap == 0:
        logging.error("No overlapping genes between dependency and genes BED. Check --genes path (use prep_out/genes.depmap.unique.bed) and gene normalization.")
        # Continue, but likely nothing will be fitted
    dep_cells = set(dep["cell_line"].astype(str).unique())
    cnv_cells = set(cnv["cell_line"].astype(str).unique())
    sv_cells = set(sv["cell_line"].astype(str).unique())
    cell_overlap = len(dep_cells & cnv_cells & sv_cells)
    logging.info("Cell overlap across dep/CNV/SV: %d", cell_overlap)
    if cell_overlap < 500:
        logging.warning("Low cell overlap (<500). Pilot recommended or check identifiers.")

    # Helper: select genes by per-gene breakpoint prevalence across intersecting cells
    def select_genes_by_prevalence(dep_df: pd.DataFrame,
                                   sv_df: pd.DataFrame,
                                   genes_bed_df: pd.DataFrame,
                                   pilot_size: int,
                                   unstable_only: bool,
                                   random_state: int = 1) -> tuple[list[str], list[str], pd.DataFrame]:
        # Exclude sex chromosomes to reduce idiosyncrasies in pilot
        genes_auto = genes_bed_df[~genes_bed_df["chrom"].isin(["chrX","chrY"])].copy()
        # Compute bp distances for autosomal genes
        bp_all = nearest_breakpoint_distance(sv_df, genes_auto, cap_bp=2_000_000, show_progress=False)
        # Limit to (gene,cell) that exist in dependency
        pairs = dep_df[["gene","cell_line"]].drop_duplicates()
        bp = pairs.merge(bp_all, on=["gene","cell_line"], how="left")
        bp["bp_dist"] = bp["bp_dist"].fillna(10_000_000)
        prev = (
            bp.groupby("gene")["bp_dist"]
              .agg(p100k=lambda x: float((x <= 100_000).mean()),
                   p1m=lambda x: float((x <= 1_000_000).mean()))
              .reset_index()
        )
        # Thresholds
        unstable_df = prev.query("p100k >= 0.05 or p1m >= 0.20").copy()
        stable_df = prev.query("p100k == 0.0 and p1m <= 0.005").copy()
        # Sample up to pilot_size
        if len(unstable_df) == 0:
            logging.warning("Pilot selection: no genes met unstable thresholds; falling back to top by p1m")
            unstable_df = prev.sort_values("p1m", ascending=False).head(pilot_size)
        unstable_df = unstable_df.sample(min(len(unstable_df), int(pilot_size)), random_state=random_state) if len(unstable_df) > int(pilot_size) else unstable_df
        if not unstable_only:
            stable_df = stable_df[~stable_df["gene"].isin(unstable_df["gene"])].copy()
            if len(stable_df) == 0:
                logging.warning("Pilot selection: no genes met stable thresholds; falling back to bottom by p1m")
                stable_df = prev[~prev["gene"].isin(unstable_df["gene"])].sort_values("p1m", ascending=True).head(pilot_size)
            stable_df = stable_df.sample(min(len(stable_df), int(pilot_size)), random_state=random_state) if len(stable_df) > int(pilot_size) else stable_df
        # Log medians
        if len(unstable_df):
            logging.info("Pilot unstable selection medians: p100k=%.4f p1m=%.4f", float(unstable_df["p100k"].median()), float(unstable_df["p1m"].median()))
        if not unstable_only and len(stable_df):
            logging.info("Pilot stable selection medians:   p100k=%.4f p1m=%.4f", float(stable_df["p100k"].median()), float(stable_df["p1m"].median()))
            try:
                assert float(unstable_df["p1m"].median()) > float(stable_df["p1m"].median()), "Stable set is not actually more stable."
            except Exception as e:
                logging.warning("Pilot selection sanity: %s", e)
        return list(unstable_df["gene"].astype(str)), (list(stable_df["gene"].astype(str)) if not unstable_only else []), prev

    if args.pilot:
        with step("Pilot: subset reanalysis"):
            # Restrict to overlapping cell lines
            keep_cells = sorted(dep_cells & cnv_cells & sv_cells)
            logging.info("Pilot: restricting to %d overlapping cell lines", len(keep_cells))
            dep_sub = dep[dep["cell_line"].isin(keep_cells)].copy()
            cnv_sub = cnv[cnv["cell_line"].isin(keep_cells)].copy()
            sv_sub = sv[sv["cell_line"].isin(keep_cells)].copy()
            if args.shuffle_sv:
                with step("Pilot: shuffle SV breakpoints (negative control)"):
                    # Shuffle start/end within each cell_line × chrom for both sides independently
                    def _shuffle_block(df: pd.DataFrame, start_col: str, end_col: str) -> pd.DataFrame:
                        idx = np.random.permutation(df.index)
                        df[start_col] = df.loc[idx, start_col].values
                        df[end_col] = df.loc[idx, end_col].values
                        return df
                    # Shuffle side 1 by cell_line×chrom1
                    sv_sub = sv_sub.groupby(["cell_line","chrom1"], group_keys=False).apply(lambda d: _shuffle_block(d, "start1","end1")).reset_index(drop=True)
                    # Shuffle side 2 by cell_line×chrom2
                    sv_sub = sv_sub.groupby(["cell_line","chrom2"], group_keys=False).apply(lambda d: _shuffle_block(d, "start2","end2")).reset_index(drop=True)
                    logging.info("Shuffled SV breakpoints within chr×cell for negative control")

            # Optional forced gene list
            forced_genes = None
            if args.pilot_genes_file:
                try:
                    with open(args.pilot_genes_file, "r", encoding="utf-8") as f:
                        forced_genes = [ln.strip().split()[0] for ln in f if ln.strip()]
                    logging.info("Pilot: forcing %d genes from --pilot-genes-file", len(forced_genes))
                except Exception as e:
                    logging.warning("Pilot: failed to read --pilot-genes-file: %s", e)

            # Choose genes by prevalence across intersecting cells (autosomes only), unless forced
            if forced_genes is not None:
                unstable_genes = sorted(set(str(g) for g in forced_genes))
                stable_genes = []
                # Build prev_table for logging using all genes for reference
                _, _, prev_table = select_genes_by_prevalence(dep_sub, sv_sub, genes, pilot_size=int(args.pilot_size), unstable_only=True)
            else:
                unstable_genes, stable_genes, prev_table = select_genes_by_prevalence(
                    dep_sub, sv_sub, genes, pilot_size=int(args.pilot_size), unstable_only=bool(args.pilot_unstable_only)
                )
            logging.info("Pilot unstable genes selected: %d", len(unstable_genes))
            if not args.pilot_unstable_only:
                logging.info("Pilot stable genes selected: %d", len(stable_genes))

            # Run sets and collect metrics for contrast checks
            unstable_metrics = {}
            stable_metrics = {}
            
            def run_set_with_metrics(tag: str, gene_list: list[str], prev_table: pd.DataFrame, metrics_dict: dict):
                if not gene_list:
                    logging.warning("Pilot %s: empty gene list; skipping", tag)
                    return
                out_dir = os.path.join(args.out, tag)
                os.makedirs(out_dir, exist_ok=True)
                with step(f"Pilot: build features [{tag}]"):
                    genes_subset = genes[genes["gene"].astype(str).isin(set(gene_list))].copy()
                    # CN cap
                    cnv_local = cnv_sub.copy()
                    cnv_local["cn"] = np.minimum(cnv_local["cn"].astype(float), float(args.cap_cn))
                    gene_cn = length_weighted_gene_cn(cnv_local, genes_subset, show_progress)
                    bp = nearest_breakpoint_distance(sv_sub, genes_subset, cap_bp=max(args.bp_windows[1], 2_000_000), show_progress=show_progress)
                    feats = gene_cn.merge(bp, on=["gene","cell_line"], how="outer")
                    # Join expression if present
                    expr = None
                    if args.expression:
                        expr = load_expression_long(args.expression)
                        feats = feats.merge(expr, on=["gene","cell_line"], how="left")
                    feats["cn"] = feats["cn"].fillna(args.default_ploidy)
                    feats["bp_dist"] = feats["bp_dist"].fillna(10_000_000)
                    # Prevalence logs over training rows (inner join: only rows present in both dep and feats)
                    train = dep_sub[dep_sub["gene"].isin(gene_list)].merge(
                        feats[["gene","cell_line","bp_dist","cn"] + (["expression"] if (args.expression and "expression" in feats.columns) else [])],
                        on=["gene","cell_line"], how="inner"
                    )
                    train["bp_dist"] = train["bp_dist"].fillna(10_000_000)
                    # Diagnostic: compare selection vs training universes
                    cells_sel = set(dep_sub["cell_line"].unique()) & set(sv_sub["cell_line"].unique())
                    sel_bp = bp[bp["gene"].isin(gene_list) & bp["cell_line"].isin(cells_sel)].copy()
                    if len(sel_bp):
                        sel_bp["bp_dist"] = sel_bp["bp_dist"].fillna(10_000_000)
                    if len(sel_bp) and len(train):
                        prev_sel = sel_bp.groupby("gene")["bp_dist"].agg(
                            p100k=lambda x: float((x <= 100_000).mean()),
                            p1m=lambda x: float((x <= 1_000_000).mean())
                        )
                        prev_trn = train.groupby("gene")["bp_dist"].agg(
                            p100k=lambda x: float((x <= 100_000).mean()),
                            p1m=lambda x: float((x <= 1_000_000).mean())
                        )
                        logging.debug("Pilot %s: selection medians p100k=%.4f p1m=%.4f (n=%d pairs)", tag,
                                     float(prev_sel["p100k"].median()), float(prev_sel["p1m"].median()), len(sel_bp))
                        logging.debug("Pilot %s: training medians p100k=%.4f p1m=%.4f (n=%d pairs)", tag,
                                     float(prev_trn["p100k"].median()), float(prev_trn["p1m"].median()), len(train))
                    if len(train):
                        p100k = float((train["bp_dist"] <= args.qc_close_bp).mean())
                        p1m = float((train["bp_dist"] <= args.bp_windows[1]).mean())
                        # Diagnostic: per-gene prevalences in training
                        pg = train.groupby("gene")["bp_dist"].agg(
                            p100k=lambda x: float((x <= 100_000).mean()),
                            p1m=lambda x: float((x <= 1_000_000).mean())
                        )
                        logging.info("Pilot %s: mean(bp<=100k)=%.4f mean(bp<=1Mb)=%.4f", tag, p100k, p1m)
                        logging.debug("Pilot %s: training per-gene medians p100k=%.4f p1m=%.4f", tag,
                                     float(pg["p100k"].median()), float(pg["p1m"].median()))
                        if args.expression:
                            expr_rate = float(train["expression"].notna().mean()) if "expression" in train.columns else 0.0
                            logging.info("Pilot %s: expression non-NA rate=%.4f", tag, expr_rate)
                            if expr_rate < 0.05:
                                logging.warning("Pilot %s: Very low expression non-NA rate; check IDs and --expression file.", tag)
                    merged = train  # Use train for downstream (inner join universe)
                with step(f"Pilot: fit models [{tag}]"):
                    pred_df, models_df = fit_models_and_predict(
                        dep_sub[dep_sub["gene"].isin(gene_list)],
                        feats,
                        args.model,
                        args.alpha,
                        args.epsilon,
                        use_expr=bool(args.expression),
                        bp_windows=tuple(args.bp_windows),
                        default_ploidy=args.default_ploidy,
                        show_progress=show_progress,
                        add_continuous_proximity=args.add_continuous_proximity,
                        standardize_predictors=args.standardize_predictors,
                    )
                with step(f"Pilot: write outputs [{tag}]"):
                    corrected = correct_dependencies(pred_df)
                    corrected_path = os.path.join(out_dir, "dependency_corrected.csv")
                    corrected.to_csv(corrected_path, index=False)
                    models_path = os.path.join(out_dir, "models_coefficients.csv")
                    models_df.to_csv(models_path, index=False)
                    # QC flags on pilot feats
                    qc = qc_flags(feats, args.default_ploidy, args.qc_cn_abs_threshold, args.qc_close_bp, args.qc_very_close_bp, args.qc_require_two)
                    qc_path = os.path.join(out_dir, "qc_flags.csv")
                    qc.to_csv(qc_path, index=False)
                    # Optional evaluation
                    eval_written = False
                    if args.essential and args.nonessential:
                        eval_df = evaluate(corrected, args.essential, args.nonessential, corrected_col="dependency_corrected")
                        if eval_df is not None:
                            eval_path = os.path.join(out_dir, "evaluation_summary.csv")
                            eval_df.to_csv(eval_path, index=False)
                            eval_written = True
                    # Tiny pilot summary
                    prox_cols = [c for c in ["bp_within_100000","bp_within_1000000","inv_bp"] if c in models_df.columns]
                    frac_nonzero = 0.0
                    if prox_cols:
                        any_nonzero = (models_df[prox_cols].abs().sum(axis=1) > 0).mean()
                        frac_nonzero = float(any_nonzero)
                    # Median |Δ corr(dep, CN)| and top 5 genes by |Δcorr|
                    med_delta_corr = np.nan
                    top5: list[tuple[str, float]] = []
                    try:
                        joined = dep_sub.merge(feats[["gene","cell_line","cn"]], on=["gene","cell_line"], how="left")
                        joined = joined.merge(corrected[["gene","cell_line","dependency_corrected"]], on=["gene","cell_line"], how="left")
                        rows = []
                        for g, sub in joined.groupby("gene"):
                            if sub["cn"].notna().sum() >= 5:
                                c1 = sub[["dependency","cn"]].dropna()
                                c2 = sub[["dependency_corrected","cn"]].dropna()
                                if len(c1) >= 5 and len(c2) >= 5:
                                    corr_before = float(np.corrcoef(c1["dependency"], c1["cn"])[0,1]) if c1["dependency"].std() > 0 and c1["cn"].std() > 0 else np.nan
                                    corr_after = float(np.corrcoef(c2["dependency_corrected"], c2["cn"])[0,1]) if c2["dependency_corrected"].std() > 0 and c2["cn"].std() > 0 else np.nan
                                    if not np.isnan(corr_before) and not np.isnan(corr_after):
                                        rows.append((g, abs(corr_before) - abs(corr_after)))
                        if rows:
                            med_delta_corr = float(np.median([r[1] for r in rows]))
                            top5 = sorted(rows, key=lambda x: x[1], reverse=True)[:5]
                    except Exception:
                        pass
                    # Selection medians from prev_table filtered to gene_list
                    sel_prev = prev_table[prev_table["gene"].isin(gene_list)].copy()
                    median_p100k = float(sel_prev["p100k"].median()) if len(sel_prev) else float("nan")
                    median_p1m = float(sel_prev["p1m"].median()) if len(sel_prev) else float("nan")
                    # Training design prevalence
                    train_mean_p100k = float((merged["bp_dist"] <= 100_000).mean()) if len(merged) else float("nan")
                    train_mean_p1m = float((merged["bp_dist"] <= 1_000_000).mean()) if len(merged) else float("nan")
                    # n_cells_training: median and IQR across genes
                    n_cells_per_gene = merged.groupby("gene")["cell_line"].nunique() if len(merged) else pd.Series(dtype=float)
                    n_cells_median = float(n_cells_per_gene.median()) if len(n_cells_per_gene) else float("nan")
                    n_cells_q25 = float(n_cells_per_gene.quantile(0.25)) if len(n_cells_per_gene) else float("nan")
                    n_cells_q75 = float(n_cells_per_gene.quantile(0.75)) if len(n_cells_per_gene) else float("nan")
                    # flag_bp_close rate over feats
                    flag_rate = float(qc["flag_bp_close"].mean()) if "flag_bp_close" in qc.columns and len(qc) else float("nan")
                    # Build design matrix (training rows only) and compute robust metrics
                    design = merged.copy()
                    logging.debug("Pilot %s: merged columns: %s", tag, list(merged.columns))
                    # Ensure dependency is present (should be from dep_sub merge, but double-check)
                    if "dependency" not in design.columns:
                        logging.debug("Pilot %s: dependency not in merged, attempting merge", tag)
                        design = design.merge(dep_sub[["gene","cell_line","dependency"]], on=["gene","cell_line"], how="left")
                    if "dependency" not in design.columns:
                        logging.warning("Pilot %s: dependency column missing from design matrix; skipping design_matrix.csv", tag)
                        design = None
                    if design is not None:
                        design["bp_within_100000"] = (design["bp_dist"] <= args.qc_close_bp).astype(float)
                        design["bp_within_1000000"] = (design["bp_dist"] <= args.bp_windows[1]).astype(float)
                        if args.add_continuous_proximity:
                            design["inv_bp"] = 1.0 / (design["bp_dist"] + 1e3)
                        else:
                            design["inv_bp"] = 0.0
                        # If standardization was used, recompute standardized columns per-gene
                        if args.standardize_predictors:
                            for col in ["cn", "bp_within_100000", "bp_within_1000000"]:
                                if col in design.columns:
                                    std_col = f"{col}_std"
                                    design[std_col] = np.nan
                                    for g, sub in design.groupby("gene"):
                                        mean_val = sub[col].mean()
                                        std_val = sub[col].std()
                                        if std_val > 0:
                                            design.loc[sub.index, std_col] = (sub[col] - mean_val) / std_val
                                        else:
                                            design.loc[sub.index, std_col] = 0.0
                        # Add dependency_corrected for compatibility with analysis scripts
                        design = design.merge(corrected[["gene","cell_line","dependency_corrected"]], on=["gene","cell_line"], how="left")
                        # Write design matrix
                        design_cols = ["gene","cell_line","dependency","dependency_corrected","cn","bp_dist","bp_within_100000","bp_within_1000000","inv_bp"]
                        if args.standardize_predictors:
                            design_cols.extend(["cn_std","bp_within_100000_std","bp_within_1000000_std"])
                            if args.add_continuous_proximity:
                                design["inv_bp_std"] = np.nan
                                for g, sub in design.groupby("gene"):
                                    mean_val = sub["inv_bp"].mean()
                                    std_val = sub["inv_bp"].std()
                                    if std_val > 0:
                                        design.loc[sub.index, "inv_bp_std"] = (sub["inv_bp"] - mean_val) / std_val
                                    else:
                                        design.loc[sub.index, "inv_bp_std"] = 0.0
                                design_cols.append("inv_bp_std")
                        design_matrix_path = os.path.join(out_dir, "design_matrix.csv")
                        design[design_cols].to_csv(design_matrix_path, index=False)
                        logging.debug("Pilot %s: wrote design_matrix.csv with %d rows, columns: %s", tag, len(design), design_cols)
                    # Compute robust proximity metrics
                    # Use lower threshold when standardization is used (contributions tend to be smaller)
                    THRESH_CONTRIB = 0.01 if args.standardize_predictors else 0.1
                    THRESH_CELL_FRAC = 0.10
                    prox_only_medabs_drop = np.nan
                    prox_active_genes_frac = np.nan
                    cn_var_median = np.nan
                    if design is not None:
                        try:
                            # Always use raw columns for contributions (convert standardized coeffs to raw scale if needed)
                            x_w100k, x_w1m, x_inv = "bp_within_100000", "bp_within_1000000", "inv_bp" if args.add_continuous_proximity else None
                            # Merge coefficients
                            keep_coef = models_df[[c for c in ["gene","bp_within_100000","bp_within_1000000","inv_bp"] if c in models_df.columns]].copy()
                            for c in ["bp_within_100000","bp_within_1000000","inv_bp"]:
                                if c not in keep_coef.columns:
                                    keep_coef[c] = 0.0
                            keep_coef = keep_coef.rename(columns={
                                "bp_within_100000": "b_w100k",
                                "bp_within_1000000": "b_w1m",
                                "inv_bp": "b_inv"
                            })
                            df_robust = design.merge(keep_coef, on="gene", how="left")
                            # If standardization was used, convert coefficients from standardized to raw scale
                            if args.standardize_predictors:
                                for col, b_col in [("bp_within_100000", "b_w100k"), ("bp_within_1000000", "b_w1m")]:
                                    if col in df_robust.columns and b_col in df_robust.columns:
                                        for g, sub in df_robust.groupby("gene"):
                                            std_val = sub[col].std()
                                            if std_val > 0:
                                                df_robust.loc[sub.index, b_col] = df_robust.loc[sub.index, b_col] / std_val
                                if x_inv and "inv_bp" in df_robust.columns and "b_inv" in df_robust.columns:
                                    for g, sub in df_robust.groupby("gene"):
                                        std_val = sub["inv_bp"].std()
                                        if std_val > 0:
                                            df_robust.loc[sub.index, "b_inv"] = df_robust.loc[sub.index, "b_inv"] / std_val
                            # Per-row proximity contribution (in dependency units, using raw columns)
                            df_robust["prox_contrib"] = (
                                df_robust["b_w100k"].astype(float) * df_robust[x_w100k].astype(float)
                                + df_robust["b_w1m"].astype(float) * df_robust[x_w1m].astype(float)
                                + (df_robust["b_inv"].astype(float) * df_robust[x_inv].astype(float) if x_inv and x_inv in df_robust.columns else 0.0)
                            )
                            # Proximity-only corrected dependency
                            df_robust["prox_only_corrected"] = df_robust["dependency"].astype(float) - df_robust["prox_contrib"].astype(float)
                            # Gene-level CN coupling before/after (proximity-only) and directional activity
                            drops = []
                            cn_vars = []
                            directional_flags = []
                            ACTIVE_FRAC_THR = THRESH_CELL_FRAC
                            for g, sub in df_robust.groupby("gene"):
                                ok1 = sub[["cn","dependency"]].dropna()
                                ok2 = sub[["cn","prox_only_corrected"]].dropna()
                                # per-gene activity (cells with |contrib| >= THRESH_CONTRIB)
                                active_frac_g = float((sub["prox_contrib"].abs() >= THRESH_CONTRIB).mean()) if len(sub) else 0.0
                                if len(ok1) > 10 and len(ok2) > 10:
                                    r_before = ok1.corr().iloc[0,1]
                                    r_after = ok2.corr().iloc[0,1]
                                    if not (np.isnan(r_before) or np.isnan(r_after)):
                                        d = abs(r_before) - abs(r_after)
                                        drops.append(d)
                                        # directional: active AND positive drop
                                        directional_flags.append(1 if (active_frac_g >= ACTIVE_FRAC_THR and d > 0) else 0)
                                    cn_vars.append(sub["cn"].var())
                            if drops:
                                prox_only_medabs_drop = float(np.median(drops))
                            if cn_vars:
                                cn_var_median = float(np.median(cn_vars))
                            # Proximity activity: genes with >=10% cells where |contrib| >= THRESH_CONTRIB
                            gene_active_frac = df_robust.assign(active=(np.abs(df_robust["prox_contrib"]) >= THRESH_CONTRIB)).groupby("gene")["active"].mean()
                            prox_active_genes_frac = float((gene_active_frac >= THRESH_CELL_FRAC).mean()) if len(gene_active_frac) else np.nan
                            # Directional variant: require positive drop too
                            if directional_flags:
                                prox_active_genes_frac_dir = float(np.mean(directional_flags))
                            else:
                                prox_active_genes_frac_dir = np.nan
                        except Exception as e:
                            logging.debug("Pilot %s: robust metrics computation failed: %s", tag, e)
                    # Store metrics for contrast checks
                    metrics_dict.update({
                        "median_p1m": median_p1m,
                        "train_mean_p1m": train_mean_p1m,
                        "flag_rate": flag_rate,
                        "frac_nonzero": frac_nonzero,
                        "med_delta_corr": med_delta_corr,
                        "prox_active_genes_frac": prox_active_genes_frac,
                        "prox_only_medabs_drop": prox_only_medabs_drop,
                        "cn_var_median": cn_var_median,
                        "prox_active_genes_frac_dir": prox_active_genes_frac_dir,
                    })
                    with open(os.path.join(out_dir, "pilot_summary.txt"), "w", encoding="utf-8") as fsum:
                        fsum.write(f"genes_fitted={len(models_df)}\n")
                        fsum.write(f"fraction_genes_with_nonzero_proximity_coeffs={frac_nonzero:.4f}\n")
                        fsum.write(f"median_p100k={median_p100k:.4f}\n")
                        fsum.write(f"median_p1m={median_p1m:.4f}\n")
                        fsum.write(f"train_mean_p100k={train_mean_p100k:.4f}\n")
                        fsum.write(f"train_mean_p1m={train_mean_p1m:.4f}\n")
                        fsum.write(f"n_cells_training_median={n_cells_median:.0f}\n")
                        fsum.write(f"n_cells_training_iqr=[{n_cells_q25:.0f}, {n_cells_q75:.0f}]\n")
                        fsum.write(f"flag_bp_close_rate={flag_rate:.4f}\n")
                        if not np.isnan(med_delta_corr):
                            fsum.write(f"median_abs_corr_dep_cn_delta_before_to_after={med_delta_corr:.4f}\n")
                            if top5:
                                fsum.write("top5_genes_by_abs_corr_drop=\n")
                                for g, d in top5:
                                    fsum.write(f"  {g}\t{d:.4f}\n")
                        # Robust proximity metrics
                        fsum.write(f"contrib_thresh={THRESH_CONTRIB:.2f}\n")
                        fsum.write(f"cell_frac_thresh={THRESH_CELL_FRAC:.2f}\n")
                        if not np.isnan(prox_only_medabs_drop):
                            fsum.write(f"prox_only_median_abs_corr_drop={prox_only_medabs_drop:.4f}\n")
                        if not np.isnan(prox_active_genes_frac):
                            fsum.write(f"prox_active_genes_frac_(|contrib|>={THRESH_CONTRIB:.2f}_in>={int(THRESH_CELL_FRAC*100)}%cells)={prox_active_genes_frac:.3f}\n")
                        if not np.isnan(prox_active_genes_frac_dir):
                            fsum.write(f"prox_active_genes_frac_dirpos_(|contrib|>={THRESH_CONTRIB:.2f}_in>={int(THRESH_CELL_FRAC*100)}%cells_&_drop>0)={prox_active_genes_frac_dir:.3f}\n")
                        if not np.isnan(cn_var_median):
                            fsum.write(f"cn_var_median={cn_var_median:.4f}\n")
                        if not eval_written:
                            fsum.write("note=evaluation_summary.csv omitted (no overlap between pilot genes and essential/nonessential sets)\n")
                    # Write manifest.json
                    manifest = {
                        "pilot_set": tag,
                        "n_genes": len(gene_list),
                        "n_genes_fitted": len(models_df),
                        "n_cells_training_median": int(n_cells_median) if not np.isnan(n_cells_median) else None,
                        "n_cells_training_iqr": [int(n_cells_q25), int(n_cells_q75)] if not (np.isnan(n_cells_q25) or np.isnan(n_cells_q75)) else None,
                        "args": {
                            "model": args.model,
                            "alpha": args.alpha,
                            "epsilon": args.epsilon,
                            "bp_windows": args.bp_windows,
                            "default_ploidy": args.default_ploidy,
                            "cap_cn": args.cap_cn,
                            "add_continuous_proximity": args.add_continuous_proximity,
                            "standardize_predictors": args.standardize_predictors,
                        },
                        "thresholds": {
                            "contrib_thresh": THRESH_CONTRIB,
                            "cell_frac_thresh": THRESH_CELL_FRAC,
                            "qc_close_bp": args.qc_close_bp,
                            "qc_very_close_bp": args.qc_very_close_bp,
                        },
                        "evaluation_included": eval_written,
                    }
                    manifest_path = os.path.join(out_dir, "manifest.json")
                    with open(manifest_path, "w", encoding="utf-8") as f:
                        json.dump(manifest, f, indent=2)
                    logging.info("Pilot %s: nonzero proximity coeffs fraction=%.3f", tag, frac_nonzero)
                    if not np.isnan(prox_active_genes_frac):
                        logging.info("Pilot %s: prox_active_genes_frac=%.3f", tag, prox_active_genes_frac)
                    if not np.isnan(prox_only_medabs_drop):
                        logging.info("Pilot %s: prox_only_median_abs_corr_drop=%.4f", tag, prox_only_medabs_drop)

            # Run sets
            run_set_with_metrics("unstable", unstable_genes, prev_table, unstable_metrics)
            if not args.pilot_unstable_only:
                run_set_with_metrics("stable", stable_genes, prev_table, stable_metrics)
                # Contrast checks
                if unstable_metrics and stable_metrics:
                    # Label uses the same thresholds as computed in pilot summary
                    label_prox = f"Prox active genes frac (|contrib|>={THRESH_CONTRIB:.2f} in >={int(THRESH_CELL_FRAC*100)}% cells)"
                    label_prox_dir = f"Directional prox-active (above + drop>0)"
                    checks = [
                        ("median_p1m", "Selection median p1m"),
                        ("train_mean_p1m", "Training mean p1m"),
                        ("flag_rate", "flag_bp_close rate"),
                        ("prox_active_genes_frac", label_prox),
                        ("prox_only_medabs_drop", "Prox-only median |Δcorr(dep,CN)|"),
                        ("prox_active_genes_frac_dir", label_prox_dir),
                    ]
                    for key, label in checks:
                        u_val = unstable_metrics.get(key, float("nan"))
                        s_val = stable_metrics.get(key, float("nan"))
                        if not (np.isnan(u_val) or np.isnan(s_val)):
                            if u_val <= s_val:
                                logging.warning("Pilot contrast check FAILED: %s unstable=%.4f <= stable=%.4f", label, u_val, s_val)
                            else:
                                logging.info("Pilot contrast check PASSED: %s unstable=%.4f > stable=%.4f", label, u_val, s_val)
                    # Context: CN variance (should be similar or lower in unstable)
                    u_cn_var = unstable_metrics.get("cn_var_median", float("nan"))
                    s_cn_var = stable_metrics.get("cn_var_median", float("nan"))
                    if not (np.isnan(u_cn_var) or np.isnan(s_cn_var)):
                        logging.info("Pilot CN variance context: unstable=%.4f, stable=%.4f", u_cn_var, s_cn_var)

        # Pilot mode ends here; skip full-genome path
        return

    expr = None
    if args.expression:
        with step("Load expression"):
            expr = load_expression_long(args.expression)
            logging.info("Expression: %d gene-cell pairs", len(expr))

    with step("Cap CN values"):
        cnv = cnv.copy()
        cnv["cn"] = np.minimum(cnv["cn"].astype(float), float(args.cap_cn))
        logging.debug("CN capped at %.2f", args.cap_cn)

    with step("Compute gene CN features"):
        gene_cn = length_weighted_gene_cn(cnv, genes, show_progress)
        logging.info("Gene CN: %d gene-cell pairs", len(gene_cn))

    with step("Compute breakpoint distances"):
        bp = nearest_breakpoint_distance(sv, genes, cap_bp=max(args.bp_windows[1], 2_000_000), show_progress=show_progress)
        logging.info("Breakpoint distances: %d gene-cell pairs", len(bp))

    with step("Assemble features"):
        feats = gene_cn.merge(bp, on=["gene","cell_line"], how="outer")
        if expr is not None:
            feats = feats.merge(expr, on=["gene","cell_line"], how="left")
        logging.info("Features: genes=%d, cells=%d, pairs=%d", feats["gene"].nunique(), feats["cell_line"].nunique(), len(feats))

    with step("Preprocess features"):
        # Sensible NA handling before fit
        feats["cn"] = feats["cn"].fillna(args.default_ploidy)
        feats["bp_dist"] = feats["bp_dist"].fillna(10_000_000)  # "far away"
        logging.debug("Filled missing CN with ploidy=%.2f, missing bp_dist with 10Mb", args.default_ploidy)

    with step("Fit per-gene robust regressions"):
        pred_df, models_df = fit_models_and_predict(
            dep, feats, args.model, args.alpha, args.epsilon, expr is not None,
            tuple(args.bp_windows), args.default_ploidy, show_progress,
            args.add_continuous_proximity, args.standardize_predictors
        )
        logging.info("Fitted models: %d genes", len(models_df))
        logging.debug("Model R2: mean=%.3f, median=%.3f", models_df["r2"].mean(), models_df["r2"].median())

    with step("Apply corrections"):
        corrected = correct_dependencies(pred_df)
        corrected_path = os.path.join(args.out, "dependency_corrected.csv")
        corrected.to_csv(corrected_path, index=False)
        logging.info("Wrote corrected matrix → %s", corrected_path)

    with step("Compute QC flags"):
        qc = qc_flags(feats, args.default_ploidy, args.qc_cn_abs_threshold, args.qc_close_bp, args.qc_very_close_bp, args.qc_require_two)
        qc_path = os.path.join(args.out, "qc_flags.csv")
        qc.to_csv(qc_path, index=False)
        logging.info("Wrote QC flags → %s", qc_path)
        logging.debug("Confounded pairs: %d (%.1f%%)", qc["confounded"].sum(), 100*qc["confounded"].mean())

    with step("Write model coefficients"):
        models_path = os.path.join(args.out, "models_coefficients.csv")
        models_df.to_csv(models_path, index=False)
        logging.info("Wrote model coefficients → %s", models_path)

    if args.essential and args.nonessential:
        with step("Evaluate vs essential/nonessential sets"):
            eval_df = evaluate(corrected, args.essential, args.nonessential, corrected_col="dependency_corrected")
            if eval_df is not None:
                eval_path = os.path.join(args.out, "evaluation_summary.csv")
                eval_df.to_csv(eval_path, index=False)
                logging.info("Evaluation summary → %s", eval_path)
                for _, row in eval_df.iterrows():
                    logging.info("  %s: AUROC=%.3f, AUPRC=%.3f", row["version"], row["auroc"], row["auprc"])


if __name__ == "__main__":
    main()


