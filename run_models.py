#!/usr/bin/env python3
"""
Run Models - Fit Huber and ElasticNetCV models per gene.

Main model: Huber regression: dependency ~ CN + proximity_kernel
Sensitivity: ElasticNetCV per gene

Outputs:
- model_coefficients.csv: coef_cn, kernel coefs, r2, n
- dependency_corrected.csv: full correction and prox-only correction
- metrics.json: summary statistics
"""
import argparse
import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import HuberRegressor, ElasticNetCV
from sklearn.preprocessing import StandardScaler

# Pre-registered proximity kernels
KERNEL_COLS = [
    "prox_exp_50k", "prox_exp_100k", "prox_exp_250k", "prox_exp_500k",
    "prox_inv", "prox_any_100k", "prox_any_250k", "prox_any_500k"
]


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


def load_table(pth: str) -> pd.DataFrame:
    """Load table from Parquet (with CSV fallback) or CSV."""
    from pathlib import Path
    p = Path(pth)
    if p.suffix == ".parquet" or p.suffix == ".parq":
        try:
            return pd.read_parquet(p)
        except (ImportError, Exception):
            # Fall back to CSV with same basename
            csv_path = str(p).replace(".parquet", ".csv").replace(".parq", ".csv")
            return pd.read_csv(csv_path)
    return pd.read_csv(p)


def load_sample_info_mapping(sample_info_path: str) -> dict:
    """Load ACH→CCLE mapping from sample_info.csv. Returns dict mapping normalized ACH to normalized CCLE."""
    try:
        df = pd.read_csv(sample_info_path, comment="#")
    except Exception:
        return {}
    
    if df.empty:
        return {}
    
    # Normalize column names
    cols_lower = {c.lower().replace(" ", "_") for c in df.columns}
    depmap_col = None
    ccle_col = None
    
    for c in df.columns:
        cl = c.lower().replace(" ", "_")
        if cl in ["depmap_id", "depmap", "model_id", "modelid", "model", "ach_id"]:
            depmap_col = c
        if cl in ["ccle_name", "ccle", "ccleid", "stripped_cell_line_name"]:
            ccle_col = c
    
    if not depmap_col or not ccle_col:
        return {}
    
    # Build mapping: normalize ACH IDs and CCLE names
    mapping = {}
    for _, row in df[[depmap_col, ccle_col]].dropna().iterrows():
        ach = str(row[depmap_col]).strip().upper()
        # Normalize ACH: ACH-000001 format
        if ach.startswith("ACH"):
            ach = ach.replace("_", "-").replace(" ", "-")
            if not "-" in ach and len(ach) > 3:
                ach = f"ACH-{ach[3:]}"
        ccle = str(row[ccle_col]).strip().upper()
        import re
        ccle = re.sub(r"[\s\.-]+", "_", ccle)
        if ach and ccle:
            mapping[ach] = ccle
    
    return mapping


def load_dependency(dep_path: str, sample_info_path: str | None = None) -> pd.DataFrame:
    """Load dependency data (long format: gene,cell_line,dependency)."""
    ext = os.path.splitext(dep_path)[1].lower()
    if ext in (".parquet", ".parq"):
        try:
            df = pd.read_parquet(dep_path)
        except (ImportError, Exception) as e:
            # Fall back to CSV with same basename if it exists
            csv_path = dep_path.replace(".parquet", ".csv").replace(".parq", ".csv")
            if os.path.exists(csv_path):
                logging.warning(f"Parquet engine not available; falling back to CSV: {csv_path}")
                try:
                    df = pd.read_csv(csv_path, encoding="utf-8")
                except UnicodeDecodeError:
                    df = pd.read_csv(csv_path, encoding="latin1")
            else:
                raise ImportError(
                    f"Parquet engine not available and CSV fallback not found: {csv_path}. "
                    f"Please install pyarrow or fastparquet, or provide a CSV file."
                ) from e
    else:
        # Be lenient with encodings
        try:
            df = pd.read_csv(dep_path, encoding="utf-8")
        except UnicodeDecodeError:
            df = pd.read_csv(dep_path, encoding="latin1")
    
    # Expect columns already standardized to ['gene','cell_line','dependency'].
    # If not, add column inference / renaming here.
    cols = {c.lower() for c in df.columns}
    if {"gene", "cell_line", "dependency"}.issubset(cols):
        colmap = {c: c.lower() for c in df.columns}
        df = df.rename(columns=colmap)
        out = df[["gene", "cell_line", "dependency"]].copy()
    elif {"gene", "cell_line", "dep"}.issubset(cols):
        # Alternative: 'dep' instead of 'dependency'
        colmap = {c: c.lower() for c in df.columns}
        df = df.rename(columns=colmap)
        df = df.rename(columns={"dep": "dependency"})
        out = df[["gene", "cell_line", "dependency"]].copy()
    else:
        raise ValueError(f"Dependency file must have columns: gene, cell_line, dependency (or dep). Got: {list(df.columns)}")
    
    out["cell_line"] = out["cell_line"].astype(str).str.strip()
    out["gene"] = out["gene"].astype(str)
    # Normalize DepMap gene names
    out["gene"] = out["gene"].str.replace(r"\s*\(\d+\)\s*$", "", regex=True).str.strip()
    
    # Map ACH IDs to CCLE names if sample_info provided
    if sample_info_path:
        ach_to_ccle = load_sample_info_mapping(sample_info_path)
        if ach_to_ccle:
            # Normalize ACH IDs in dependency
            out["cell_line_ach"] = out["cell_line"].str.upper().str.replace("_", "-").str.replace(" ", "-")
            out["cell_line"] = out["cell_line_ach"].map(ach_to_ccle).fillna(out["cell_line"])
            out = out.drop(columns=["cell_line_ach"])
            logging.info(f"Mapped {out['cell_line'].isin(ach_to_ccle.values()).sum()} dependency cell lines via sample_info")
    
    return out


def fit_huber_per_gene(df: pd.DataFrame, kernel_col: str, default_ploidy: float = 2.0,
                       min_cells: int = 200, alpha: float = 1e-4, epsilon: float = 1.35) -> tuple:
    """
    Fit Huber regression per gene: dependency ~ CN + proximity_kernel
    
    Returns:
        (coefficients_df, predictions_df)
    """
    coef_rows = []
    pred_rows = []
    
    genes = sorted(df["gene"].unique())
    logging.info(f"Fitting Huber models for {len(genes)} genes (min_cells={min_cells})...")
    
    for gene in genes:
        gdf = df[df["gene"] == gene].copy()
        
        # Drop rows with missing dependency
        gdf = gdf.dropna(subset=["dependency"])
        
        if len(gdf) < min_cells:
            # Skip but record
            coef_rows.append({
                "gene": gene,
                "n": len(gdf),
                "r2": np.nan,
                "coef_cn": np.nan,
                f"coef_{kernel_col}": np.nan,
                "intercept": np.nan
            })
            continue
        
        # Prepare features
        X = gdf[["cn", kernel_col]].copy()
        X["cn"] = X["cn"].fillna(default_ploidy)
        X[kernel_col] = X[kernel_col].fillna(0.0)
        
        y = gdf["dependency"].values
        X_vals = X.values
        
        try:
            model = HuberRegressor(alpha=alpha, epsilon=epsilon, max_iter=200)
            model.fit(X_vals, y)
            
            # Compute R²
            y_pred = model.predict(X_vals)
            y_mean = np.mean(y)
            ss_res = np.sum((y - y_pred) ** 2)
            ss_tot = np.sum((y - y_mean) ** 2)
            r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else np.nan
            
            coef_rows.append({
                "gene": gene,
                "n": len(gdf),
                "r2": float(r2),
                "coef_cn": float(model.coef_[0]),
                f"coef_{kernel_col}": float(model.coef_[1]),
                "intercept": float(model.intercept_)
            })
            
            # Store predictions for correction
            gdf["pred_structural"] = y_pred
            gdf["prox_contrib"] = model.coef_[1] * X[kernel_col].values
            pred_rows.append(gdf[["gene", "cell_line", "dependency", "pred_structural", "prox_contrib"]])
            
        except Exception as e:
            logging.warning(f"Error fitting {gene}: {e}")
            coef_rows.append({
                "gene": gene,
                "n": len(gdf),
                "r2": np.nan,
                "coef_cn": np.nan,
                f"coef_{kernel_col}": np.nan,
                "intercept": np.nan
            })
    
    coef_df = pd.DataFrame(coef_rows)
    pred_df = pd.concat(pred_rows, ignore_index=True) if pred_rows else pd.DataFrame()
    
    return coef_df, pred_df


def fit_elasticnet_per_gene(df: pd.DataFrame, kernel_col: str, default_ploidy: float = 2.0,
                            min_cells: int = 200) -> pd.DataFrame:
    """
    Fit ElasticNetCV per gene for sensitivity analysis.
    
    Returns:
        DataFrame with ElasticNet coefficients
    """
    coef_rows = []
    
    genes = sorted(df["gene"].unique())
    logging.info(f"Fitting ElasticNetCV models for {len(genes)} genes (sensitivity)...")
    
    for gene in genes:
        gdf = df[df["gene"] == gene].copy()
        gdf = gdf.dropna(subset=["dependency"])
        
        if len(gdf) < min_cells:
            continue
        
        X = gdf[["cn", kernel_col]].copy()
        X["cn"] = X["cn"].fillna(default_ploidy)
        X[kernel_col] = X[kernel_col].fillna(0.0)
        
        y = gdf["dependency"].values
        X_vals = X.values
        
        try:
            # Standardize for ElasticNet
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X_vals)
            
            model = ElasticNetCV(cv=5, random_state=42, max_iter=1000)
            model.fit(X_scaled, y)
            
            # Transform coefficients back to original scale
            coef_original = model.coef_ / scaler.scale_
            intercept_original = model.intercept_ - np.sum(coef_original * scaler.mean_)
            
            coef_rows.append({
                "gene": gene,
                "n": len(gdf),
                "coef_cn_elasticnet": float(coef_original[0]),
                f"coef_{kernel_col}_elasticnet": float(coef_original[1]),
                "intercept_elasticnet": float(intercept_original),
                "alpha_elasticnet": float(model.alpha_),
                "l1_ratio_elasticnet": float(model.l1_ratio_)
            })
        except Exception as e:
            logging.debug(f"ElasticNetCV failed for {gene}: {e}")
    
    return pd.DataFrame(coef_rows)


def compute_corrected_dependencies(pred_df: pd.DataFrame, coef_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute dependency_corrected and prox_only_corrected.
    
    dependency_corrected = dependency - (pred_structural - median_pred_per_gene)
    prox_only_corrected = dependency - prox_contrib
    """
    if len(pred_df) == 0:
        # Return empty DataFrame with expected columns
        return pd.DataFrame(columns=["gene", "cell_line", "dependency", "dependency_corrected", 
                                     "prox_only_corrected", "pred_structural", "prox_contrib"])
    
    # Merge predictions with coefficients
    if len(coef_df) > 0 and "n" in coef_df.columns:
        df = pred_df.merge(coef_df[["gene", "n"]], on="gene", how="left")
    else:
        df = pred_df.copy()
        df["n"] = np.nan
    
    # Full correction: subtract structural prediction (centered per gene)
    gene_medians = df.groupby("gene")["pred_structural"].median().rename("gene_median_pred")
    df = df.merge(gene_medians, on="gene", how="left")
    df["dependency_corrected"] = df["dependency"] - (df["pred_structural"] - df["gene_median_pred"])
    
    # Prox-only correction
    df["prox_only_corrected"] = df["dependency"] - df["prox_contrib"]
    
    return df[["gene", "cell_line", "dependency", "dependency_corrected", 
               "prox_only_corrected", "pred_structural", "prox_contrib"]]


def main():
    ap = argparse.ArgumentParser(
        description="Fit Huber and ElasticNetCV models per gene"
    )
    ap.add_argument("--dep", required=True,
                    help="Dependency file (long: gene,cell_line,dependency)")
    ap.add_argument("--design", required=True,
                    help="Design matrix CSV/Parquet from sv_ingest_wgs.py")
    ap.add_argument("--sample-info", type=str, default=None,
                    help="Optional sample_info.csv to map ACH IDs (dependency) to CCLE names (design)")
    ap.add_argument("--out", required=True,
                    help="Output directory")
    ap.add_argument("--model", choices=["huber", "both"], default="huber",
                    help="Model type: huber (main) or both (huber + elasticnet)")
    ap.add_argument("--kernel", default="prox_exp_100k",
                    help="Primary proximity kernel to use (default: prox_exp_100k)")
    ap.add_argument("--min-cells", type=int, default=200,
                    help="Minimum cells per gene to fit (default: 200)")
    ap.add_argument("--alpha", type=float, default=1e-4,
                    help="Huber regularization alpha (default: 1e-4)")
    ap.add_argument("--epsilon", type=float, default=1.35,
                    help="Huber epsilon (default: 1.35)")
    ap.add_argument("--default-ploidy", type=float, default=2.0,
                    help="Default ploidy for missing CN (default: 2.0)")
    ap.add_argument("-v", "--verbose", action="count", default=1,
                    help="Increase verbosity")
    args = ap.parse_args()
    
    setup_logging(args.verbose)
    
    os.makedirs(args.out, exist_ok=True)
    
    # Load data
    logging.info("Loading dependency data...")
    dep = load_dependency(args.dep, sample_info_path=args.sample_info)
    
    logging.info("Loading design matrix...")
    design = load_table(args.design)
    
    # Check kernel column exists
    if args.kernel not in design.columns:
        available = [c for c in design.columns if c.startswith("prox_")]
        raise ValueError(f"Kernel '{args.kernel}' not found. Available: {available}")
    
    # Aggregate design matrix across chromosomes (if chrom column exists)
    # For proximity kernels: take max (any breakpoint nearby)
    # For CN: take mean (average copy number across chromosomes)
    if "chrom" in design.columns:
        logging.info("Aggregating design matrix across chromosomes...")
        agg_dict = {}
        
        # CN: mean across chromosomes
        if "cn" in design.columns:
            agg_dict["cn"] = "mean"
        
        # Proximity kernels: max (any breakpoint nearby)
        prox_cols = [c for c in design.columns if c.startswith("prox_")]
        for col in prox_cols:
            agg_dict[col] = "max"
        
        # bp_dist: min (nearest breakpoint)
        if "bp_dist" in design.columns:
            agg_dict["bp_dist"] = "min"
        
        # has_same_chr_bp: max (any chromosome has breakpoint)
        if "has_same_chr_bp" in design.columns:
            agg_dict["has_same_chr_bp"] = "max"
        
        # Keep other columns as first (should be constant per gene,cell_line)
        other_cols = [c for c in design.columns if c not in agg_dict and c not in ["gene", "cell_line", "chrom"]]
        for col in other_cols:
            agg_dict[col] = "first"
        
        design_agg = design.groupby(["gene", "cell_line"], as_index=False).agg(agg_dict)
        logging.info(f"Aggregated from {len(design)} (gene,cell,chrom) rows to {len(design_agg)} (gene,cell) rows")
        design = design_agg
    else:
        logging.info("No 'chrom' column found; using design matrix as-is")
    
    # Normalize cell_line names for merge (uppercase, consistent separators)
    dep["cell_line"] = dep["cell_line"].str.upper().str.replace(r"[\s\.-]+", "_", regex=True)
    design["cell_line"] = design["cell_line"].astype(str).str.upper().str.replace(r"[\s\.-]+", "_", regex=True)
    
    # Merge dependency with design
    df = dep.merge(design, on=["gene", "cell_line"], how="inner")
    logging.info(f"Merged data: {len(df)} gene-cell pairs, {df['gene'].nunique()} genes, {df['cell_line'].nunique()} cells")
    
    if len(df) == 0:
        logging.warning("No overlapping gene-cell pairs after merge!")
        logging.warning(f"  Dependency: {dep['gene'].nunique()} genes, {dep['cell_line'].nunique()} cells")
        logging.warning(f"  Design: {design['gene'].nunique()} genes, {design['cell_line'].nunique()} cells")
        logging.warning(f"  Sample dep cells: {sorted(dep['cell_line'].unique())[:5]}")
        logging.warning(f"  Sample design cells: {sorted(design['cell_line'].unique())[:5]}")
        raise ValueError("No overlapping gene-cell pairs. Check cell_line name normalization.")
    
    # Fit Huber models
    logging.info(f"Fitting Huber models with kernel: {args.kernel}")
    coef_huber, pred_huber = fit_huber_per_gene(
        df, args.kernel, default_ploidy=args.default_ploidy,
        min_cells=args.min_cells, alpha=args.alpha, epsilon=args.epsilon
    )
    
    # Compute corrected dependencies
    logging.info("Computing corrected dependencies...")
    corrected = compute_corrected_dependencies(pred_huber, coef_huber)
    
    # Fit ElasticNet if requested
    coef_elasticnet = None
    if args.model == "both":
        logging.info("Fitting ElasticNetCV models (sensitivity)...")
        coef_elasticnet = fit_elasticnet_per_gene(
            df, args.kernel, default_ploidy=args.default_ploidy, min_cells=args.min_cells
        )
        # Merge ElasticNet coefficients
        coef_huber = coef_huber.merge(coef_elasticnet, on="gene", how="left")
    
    # Save outputs
    coef_path = os.path.join(args.out, "model_coefficients.csv")
    coef_huber.to_csv(coef_path, index=False)
    logging.info(f"Saved coefficients to {coef_path}")
    
    corr_path = os.path.join(args.out, "dependency_corrected.csv")
    corrected.to_csv(corr_path, index=False)
    logging.info(f"Saved corrected dependencies to {corr_path}")
    
    # Save metrics
    metrics = {
        "n_genes": int(df["gene"].nunique()),
        "n_cells": int(df["cell_line"].nunique()),
        "n_pairs": int(len(df)),
        "kernel_used": args.kernel,
        "min_cells": args.min_cells,
        "n_genes_fitted": int(coef_huber["n"].notna().sum()),
        "mean_r2": float(coef_huber["r2"].mean()) if coef_huber["r2"].notna().any() else None,
        "model_type": args.model
    }
    
    metrics_path = os.path.join(args.out, "metrics.json")
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)
    logging.info(f"Saved metrics to {metrics_path}")


if __name__ == "__main__":
    main()

