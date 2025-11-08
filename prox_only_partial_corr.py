#!/usr/bin/env python3
"""
Compute proximity-only partial correlation metric.

Measures Δ|partial r(dep, CN | prox)| after removing only the proximity component
(not the CN component). This isolates the positional artifact from CN-driven effects.
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LinearRegression


def resid(y, X):
    """Residualize y against X using OLS."""
    X_with_intercept = np.c_[np.ones(len(X)), X]
    beta = np.linalg.lstsq(X_with_intercept, y, rcond=None)[0]
    return y - X_with_intercept @ beta


def compute_prox_only_correction(dep_corr_path: str, coef_path: str, 
                                 bp_windows: tuple = (250000, 2000000)) -> pd.DataFrame:
    """
    Compute proximity-only correction from dependency_corrected.csv and model_coefficients.csv.
    
    Args:
        dep_corr_path: Path to dependency_corrected.csv
        coef_path: Path to model_coefficients.csv
        bp_windows: Tuple of (near_window, far_window) in bp
        
    Returns:
        DataFrame with prox_only_corrected column
    """
    w1, w2 = bp_windows
    
    # Load coefficients
    coef = pd.read_csv(coef_path)
    
    # Load dependency data (chunked if too large)
    print(f"Loading dependency data from {dep_corr_path}...")
    try:
        df = pd.read_csv(dep_corr_path, usecols=["gene", "cell_line", "dependency", "cn", "bp_dist"])
    except MemoryError:
        # If too large, process in chunks
        print("File too large, processing in chunks...")
        chunks = []
        for chunk in pd.read_csv(dep_corr_path, chunksize=100000, 
                                usecols=["gene", "cell_line", "dependency", "cn", "bp_dist"]):
            chunks.append(chunk)
        df = pd.concat(chunks, ignore_index=True)
    
    # Create proximity flags
    df["bp_dist"] = df["bp_dist"].fillna(5_000_000)  # Cap at max distance
    df["bp_near"] = (df["bp_dist"] <= w1).astype(float)
    df["bp_far"] = ((df["bp_dist"] > w1) & (df["bp_dist"] <= w2)).astype(float)
    
    # Merge coefficients
    coef_subset = coef[["gene", "coef_bp_near", "coef_bp_far"]].copy()
    df = df.merge(coef_subset, on="gene", how="left")
    
    # Fill missing coefficients with 0
    df["coef_bp_near"] = df["coef_bp_near"].fillna(0.0)
    df["coef_bp_far"] = df["coef_bp_far"].fillna(0.0)
    
    # Compute proximity contribution
    df["prox_contrib"] = (df["coef_bp_near"] * df["bp_near"] + 
                          df["coef_bp_far"] * df["bp_far"])
    
    # Gene-level median proximity contribution
    med_prox = df.groupby("gene")["prox_contrib"].median().rename("med_prox")
    df = df.merge(med_prox, left_on="gene", right_index=True, how="left")
    
    # Proximity-only correction
    df["prox_only_corrected"] = df["dependency"] - (df["prox_contrib"] - df["med_prox"])
    
    return df


def compute_prox_only_partial_corr(df: pd.DataFrame, min_rows: int = 20) -> dict:
    """
    Compute Δ|partial r(dep, CN | prox)| using proximity-only corrected dependency.
    
    Args:
        df: DataFrame with dependency, prox_only_corrected, cn, bp_near, bp_far
        min_rows: Minimum rows per gene
        
    Returns:
        Dictionary with summary statistics
    """
    results = []
    rng = np.random.default_rng(42)
    
    for gene, sub in df.groupby("gene"):
        valid = sub[["dependency", "prox_only_corrected", "cn", "bp_near", "bp_far"]].dropna()
        if len(valid) < min_rows:
            continue
        
        try:
            Xprox = valid[["bp_near", "bp_far"]].values
            
            # Partial correlation: residualize both dep and cn against proximity
            dep_resid = resid(valid["dependency"].values, Xprox)
            dep_prox_corr_resid = resid(valid["prox_only_corrected"].values, Xprox)
            cn_resid = resid(valid["cn"].values, Xprox)
            
            # Compute partial correlations
            r0 = np.corrcoef(dep_resid, cn_resid)[0, 1]
            r1 = np.corrcoef(dep_prox_corr_resid, cn_resid)[0, 1]
            
            if np.isnan(r0) or np.isnan(r1):
                continue
            
            delta = abs(r0) - abs(r1)
            results.append({
                "gene": gene,
                "prox_only_partial_delta": delta,
                "r0_abs": abs(r0),
                "r1_abs": abs(r1),
                "n_obs": len(valid)
            })
        except Exception as e:
            continue
    
    if not results:
        return {
            "n_genes": 0,
            "median_prox_only_delta": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "frac_positive": np.nan
        }
    
    results_df = pd.DataFrame(results)
    vals = results_df["prox_only_partial_delta"].values
    
    # Bootstrap median (2,000 reps)
    n = len(vals)
    boot_medians = []
    for _ in range(2000):
        idx = rng.integers(0, n, size=n)
        boot_medians.append(np.median(vals[idx]))
    
    boot_medians = np.array(boot_medians)
    
    summary = {
        "n_genes": int(len(results_df)),
        "median_prox_only_delta": float(np.median(vals)),
        "mean_prox_only_delta": float(np.mean(vals)),
        "frac_positive": float((vals > 0).mean()),
        "ci_low": float(np.percentile(boot_medians, 2.5)),
        "ci_high": float(np.percentile(boot_medians, 97.5))
    }
    
    return summary, results_df


def compute_prevalence_match(df_true: pd.DataFrame, df_shuffle: pd.DataFrame,
                            tolerance: float = 0.02) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """
    Filter genes by proximity prevalence match.
    
    Args:
        df_true: True data DataFrame
        df_shuffle: Shuffled data DataFrame
        tolerance: Maximum difference in prevalence
        
    Returns:
        Filtered dataframes and count of common genes
    """
    def compute_prevalence(df):
        df = df.copy()
        df["prox_active"] = ((df["bp_near"] > 0) | (df["bp_far"] > 0)).astype(float)
        return df.groupby("gene")["prox_active"].mean()
    
    true_prev = compute_prevalence(df_true)
    shuf_prev = compute_prevalence(df_shuffle)
    
    # Find genes with matching prevalence
    common_genes = set(true_prev.index) & set(shuf_prev.index)
    matched_genes = []
    
    for gene in common_genes:
        diff = abs(true_prev[gene] - shuf_prev[gene])
        if diff <= tolerance:
            matched_genes.append(gene)
    
    # Filter dataframes
    true_filtered = df_true[df_true["gene"].isin(matched_genes)].copy()
    shuf_filtered = df_shuffle[df_shuffle["gene"].isin(matched_genes)].copy()
    
    return true_filtered, shuf_filtered, len(matched_genes)


def compute_cn_variance_stratum(df: pd.DataFrame, quartile: int = 2) -> pd.DataFrame:
    """
    Filter to genes in lowest quartiles of CN variance.
    
    Args:
        df: DataFrame with gene, cn columns
        quartile: Number of quartiles to keep (1=Q1, 2=Q1+Q2, etc.)
        
    Returns:
        Filtered DataFrame
    """
    cn_var = df.groupby("gene")["cn"].var()
    thresholds = cn_var.quantile([0.25, 0.5, 0.75])
    
    if quartile == 1:
        threshold = thresholds[0.25]
        keep_genes = cn_var[cn_var <= threshold].index
    elif quartile == 2:
        threshold = thresholds[0.5]
        keep_genes = cn_var[cn_var <= threshold].index
    elif quartile == 3:
        threshold = thresholds[0.75]
        keep_genes = cn_var[cn_var <= threshold].index
    else:
        keep_genes = cn_var.index  # All genes
    
    return df[df["gene"].isin(keep_genes)].copy()


def main():
    parser = argparse.ArgumentParser(
        description="Compute proximity-only partial correlation metric",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("base_folder", type=str, help="Base folder with dependency_corrected.csv and model_coefficients.csv")
    parser.add_argument("--bp-windows", nargs=2, type=int, default=[250000, 2000000],
                       metavar=("W_NEAR", "W_FAR"), help="Breakpoint proximity windows (bp)")
    parser.add_argument("--min-rows", type=int, default=20, help="Minimum rows per gene")
    parser.add_argument("--prevalence-match", type=str, default=None,
                       help="Path to shuffled data folder for prevalence matching")
    parser.add_argument("--prevalence-tolerance", type=float, default=0.02,
                       help="Tolerance for prevalence matching")
    parser.add_argument("--low-cn-variance", type=int, default=None, choices=[1, 2, 3],
                       help="Filter to lowest N quartiles of CN variance")
    parser.add_argument("--output", type=str, default=None,
                       help="Output directory (default: base_folder)")
    
    args = parser.parse_args()
    base_folder = Path(args.base_folder)
    output_dir = Path(args.output) if args.output else base_folder
    
    bp_windows = tuple(args.bp_windows)
    
    print(f"Processing: {base_folder.name}")
    print(f"BP windows: {bp_windows[0]:,}bp, {bp_windows[1]:,}bp")
    
    # Load and compute proximity-only correction
    dep_corr_path = base_folder / "dependency_corrected.csv"
    coef_path = base_folder / "model_coefficients.csv"
    
    if not dep_corr_path.exists() or not coef_path.exists():
        raise FileNotFoundError(f"Required files not found in {base_folder}")
    
    df = compute_prox_only_correction(str(dep_corr_path), str(coef_path), bp_windows)
    print(f"Loaded {len(df):,} rows for {df['gene'].nunique():,} genes")
    
    # Apply filters
    if args.low_cn_variance:
        df_before = len(df)
        df = compute_cn_variance_stratum(df, quartile=args.low_cn_variance)
        print(f"Low CN variance filter (Q1-Q{args.low_cn_variance}): "
              f"{len(df):,} rows ({df['gene'].nunique():,} genes) from {df_before:,}")
    
    if args.prevalence_match:
        shuffle_folder = Path(args.prevalence_match)
        shuffle_dep = shuffle_folder / "dependency_corrected.csv"
        shuffle_coef = shuffle_folder / "model_coefficients.csv"
        
        if not shuffle_dep.exists() or not shuffle_coef.exists():
            raise FileNotFoundError(f"Required files not found in {shuffle_folder}")
        
        df_shuffle = compute_prox_only_correction(str(shuffle_dep), str(shuffle_coef), bp_windows)
        
        if args.low_cn_variance:
            df_shuffle = compute_cn_variance_stratum(df_shuffle, quartile=args.low_cn_variance)
        
        df, df_shuffle, n_matched = compute_prevalence_match(
            df, df_shuffle, tolerance=args.prevalence_tolerance
        )
        print(f"Prevalence matching (tolerance ±{args.prevalence_tolerance}): "
              f"{n_matched:,} genes kept")
    
    # Compute partial correlation
    summary, per_gene = compute_prox_only_partial_corr(df, min_rows=args.min_rows)
    
    # Print summary
    print(f"\n{'='*80}")
    print(f"PROXIMITY-ONLY PARTIAL CORRELATION: {base_folder.name}")
    print(f"{'='*80}")
    print(f"n_genes: {summary['n_genes']:,}")
    print(f"median Δ|partial r|: {summary['median_prox_only_delta']:.4f}")
    print(f"95% CI: [{summary['ci_low']:.4f}, {summary['ci_high']:.4f}]")
    print(f"mean Δ|partial r|: {summary['mean_prox_only_delta']:.4f}")
    print(f"fraction > 0: {summary['frac_positive']:.3f}")
    print(f"{'='*80}\n")
    
    # Save results
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(output_dir / "prox_only_partial_corr_summary.csv", index=False)
    per_gene.to_csv(output_dir / "prox_only_partial_corr_per_gene.csv", index=False)
    
    print(f"Results saved to {output_dir}")


if __name__ == "__main__":
    main()

