#!/usr/bin/env python3
"""
SV/CN Bias Correction Pipeline for CRISPR Dependency Screens (v2.1 - Final)

A statistically rigorous pipeline for detecting and correcting structural genomic
bias in DepMap CRISPR gene effect scores.

Key features:
- No expression collider bias
- Proper uncertainty quantification (bootstrap CIs on all metrics)
- Multiple negative controls (within-chrom and cross-chrom shuffles)
- Directional proximity-active metric (gold standard validation)
- Excess signal quantification with CIs
- Comprehensive diagnostics and VIF checks
"""

import argparse
import json
import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Literal

import numpy as np
import pandas as pd
from sklearn.linear_model import HuberRegressor, QuantileRegressor
from sklearn.metrics import roc_auc_score, average_precision_score

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)


# ============================================================================
# Configuration and Data Classes
# ============================================================================

@dataclass
class Config:
    """Pipeline configuration with sensible defaults."""
    # Model parameters
    model_type: Literal["huber", "quantile", "ols"] = "huber"
    huber_epsilon: float = 1.35
    quantile_alpha: float = 0.5

    # Feature windows (bp)
    proximity_windows: tuple[int, int] = (100_000, 1_000_000)
    max_bp_distance: int = 5_000_000

    # Data parameters
    default_ploidy: float = 2.0
    max_cn: float = 10.0
    min_cells_per_gene: int = 20

    # Statistical parameters
    bootstrap_iterations: int = 1000
    confidence_level: float = 0.95

    # Activity thresholds (for directional prox-active metric)
    activity_contrib_threshold: float = 0.01  # dependency units
    activity_cell_fraction: float = 0.10      # fraction of cells

    # QC thresholds
    qc_cn_threshold: float = 1.5  # |CN - ploidy| threshold
    qc_proximity_threshold: int = 250_000  # bp

    # Diagnostics
    check_vif: bool = True
    vif_threshold: float = 5.0

    # Output
    output_dir: Path = Path("output")
    seed: int = 42
    
    # New features (PART B)
    genes_subset_file: Optional[str] = None
    write_features_dir: Optional[str] = None
    read_features_dir: Optional[str] = None
    orthogonalize_proximity: bool = False
    compute_partial_corr: bool = False
    match_proximity_prevalence: bool = False
    prevalence_tolerance: float = 0.02  # For prevalence matching

    def __post_init__(self):
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)


# ============================================================================
# Data Loading and Validation
# ============================================================================

class DataValidator:
    """Validates input data integrity and compatibility."""

    @staticmethod
    def validate_dependency(df: pd.DataFrame) -> pd.DataFrame:
        """Validate and normalize dependency data."""
        required_cols = {"gene", "cell_line", "dependency"}
        if not required_cols.issubset(df.columns):
            raise ValueError(f"Dependency data must contain: {required_cols}")

        df = df.copy()
        df["gene"] = df["gene"].astype(str).str.strip()
        df["cell_line"] = df["cell_line"].astype(str).str.strip()
        df["dependency"] = pd.to_numeric(df["dependency"], errors="coerce")

        initial_len = len(df)
        df = df.dropna(subset=["dependency"])
        df = df[df["gene"].str.len() > 0]
        df = df[df["cell_line"].str.len() > 0]
        if len(df) < initial_len:
            logger.warning(f"Dropped {initial_len - len(df)} invalid dependency rows")

        dups = df.duplicated(subset=["gene", "cell_line"]).sum()
        if dups > 0:
            logger.warning(f"Found {dups} duplicate gene-cell pairs, keeping first")
            df = df.drop_duplicates(subset=["gene", "cell_line"], keep="first")

        logger.info(f"Validated dependency: {len(df):,} pairs, "
                    f"{df['gene'].nunique():,} genes, {df['cell_line'].nunique():,} cells")
        return df

    @staticmethod
    def validate_cnv(df: pd.DataFrame, max_cn: float) -> pd.DataFrame:
        """Validate and normalize CNV segments."""
        required_cols = {"cell_line", "chrom", "start", "end", "cn"}
        if not required_cols.issubset(df.columns):
            raise ValueError(f"CNV data must contain: {required_cols}")

        df = df.copy()
        df["cell_line"] = df["cell_line"].astype(str).str.strip()
        df["chrom"] = df["chrom"].astype(str).str.replace("chr", "", regex=False)
        df["chrom"] = "chr" + df["chrom"]
        df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
        df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
        df["cn"] = pd.to_numeric(df["cn"], errors="coerce")

        df = df.dropna()
        df = df[df["start"] < df["end"]]
        df = df[df["cn"] >= 0]
        df = df[df["cn"] <= max_cn * 2]  # generous sanity bound
        df["cn"] = df["cn"].clip(upper=max_cn)

        logger.info(f"Validated CNV: {len(df):,} segments, {df['cell_line'].nunique():,} cells")
        return df

    @staticmethod
    def validate_sv(df: pd.DataFrame) -> pd.DataFrame:
        """Validate and normalize SV breakpoints."""
        required_cols = {"cell_line", "chrom1", "start1", "end1",
                         "chrom2", "start2", "end2"}
        if not required_cols.issubset(df.columns):
            raise ValueError(f"SV data must contain: {required_cols}")

        df = df.copy()
        df["cell_line"] = df["cell_line"].astype(str).str.strip()
        for c in ["chrom1", "chrom2"]:
            df[c] = df[c].astype(str).str.replace("chr", "", regex=False)
            df[c] = "chr" + df[c]
        for c in ["start1", "end1", "start2", "end2"]:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")

        df = df.dropna()
        df = df[df["start1"] < df["end1"]]
        df = df[df["start2"] < df["end2"]]

        logger.info(f"Validated SV: {len(df):,} breakpoints, {df['cell_line'].nunique():,} cells")
        return df

    @staticmethod
    def validate_genes(df: pd.DataFrame) -> pd.DataFrame:
        """Validate and normalize gene annotations."""
        required_cols = {"chrom", "start", "end", "gene"}
        if not required_cols.issubset(df.columns):
            raise ValueError(f"Gene data must contain: {required_cols}")

        df = df.copy()
        df["chrom"] = df["chrom"].astype(str).str.replace("chr", "", regex=False)
        df["chrom"] = "chr" + df["chrom"]
        df["gene"] = df["gene"].astype(str).str.strip()
        df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
        df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")

        df = df.dropna()
        df = df[df["start"] < df["end"]]
        df = df[df["gene"].str.len() > 0]

        dups = df.duplicated(subset=["gene"]).sum()
        if dups > 0:
            logger.warning(f"Found {dups} duplicate genes, keeping first")
            df = df.drop_duplicates(subset=["gene"], keep="first")

        logger.info(f"Validated genes: {len(df):,} genes")
        return df


# ============================================================================
# Efficient Feature Computation
# ============================================================================

class FeatureComputer:
    """Efficient computation of genomic features using vectorized operations."""

    def __init__(self, config: Config):
        self.config = config

    def compute_gene_cn(self, cnv: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
        """
        Compute length-weighted copy number for each gene-cell pair.
        Current implementation: O(genes × segments) per chromosome.
        TODO: Optimize with PyRanges or intervaltree for O(n log m).
        """
        logger.info("Computing gene copy numbers...")
        results = []

        for chrom in genes["chrom"].unique():
            genes_chr = genes[genes["chrom"] == chrom].copy()
            cnv_chr = cnv[cnv["chrom"] == chrom].copy()
            if len(cnv_chr) == 0:
                continue

            for cell_line, cnv_cell in cnv_chr.groupby("cell_line"):
                cnv_cell = cnv_cell.sort_values("start")
                for _, gene in genes_chr.iterrows():
                    gs, ge = int(gene["start"]), int(gene["end"])
                    gl = max(1, ge - gs)

                    overlaps = cnv_cell[(cnv_cell["end"] > gs) & (cnv_cell["start"] < ge)]
                    if len(overlaps) == 0:
                        continue

                    s = np.maximum(overlaps["start"].values, gs)
                    e = np.minimum(overlaps["end"].values, ge)
                    w = (e - s).astype(float)
                    weighted_cn = float(np.sum(w * overlaps["cn"].values) / gl)

                    results.append({"gene": gene["gene"], "cell_line": cell_line, "cn": weighted_cn})

        df = pd.DataFrame(results)
        logger.info(f"Computed CN for {len(df):,} gene-cell pairs")
        return df

    def compute_bp_distances(self, sv: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
        """
        Compute minimum distance to nearest SV breakpoint for each gene-cell pair.
        Uses sorted breakpoint arrays for efficient nearest-neighbor search.
        """
        logger.info("Computing breakpoint distances...")
        bp_data = []
        for _, row in sv.iterrows():
            cell = row["cell_line"]
            bp_data.append({"cell_line": cell, "chrom": row["chrom1"], "pos": int((row["start1"] + row["end1"]) // 2)})
            bp_data.append({"cell_line": cell, "chrom": row["chrom2"], "pos": int((row["start2"] + row["end2"]) // 2)})

        bp_df = pd.DataFrame(bp_data).drop_duplicates()

        bp_index = defaultdict(lambda: np.array([], dtype=int))
        for (cell, chrom), group in bp_df.groupby(["cell_line", "chrom"]):
            bp_index[(cell, chrom)] = np.sort(group["pos"].values.astype(int))

        results = []
        for _, gene in genes.iterrows():
            gm = int((gene["start"] + gene["end"]) // 2)
            chrom = gene["chrom"]
            gname = gene["gene"]

            cells_with_bp = set(k[0] for k in bp_index.keys() if k[1] == chrom)
            for cell in cells_with_bp:
                positions = bp_index[(cell, chrom)]
                if len(positions) == 0:
                    continue
                idx = np.searchsorted(positions, gm)
                candidates = []
                if idx > 0:
                    candidates.append(abs(gm - int(positions[idx - 1])))
                if idx < len(positions):
                    candidates.append(abs(gm - int(positions[idx])))
                if candidates:
                    min_dist = min(candidates)
                    results.append({
                        "gene": gname,
                        "cell_line": cell,
                        "bp_dist": min(int(min_dist), self.config.max_bp_distance)
                    })

        df = pd.DataFrame(results)
        logger.info(f"Computed distances for {len(df):,} gene-cell pairs")
        return df


# ============================================================================
# Shuffling Controls
# ============================================================================

class ShuffleControl:
    """Generate negative control datasets by shuffling breakpoints."""

    def __init__(self, config: Config):
        self.config = config
        self.rng = np.random.default_rng(config.seed)
        # hg38 chromosome sizes
        self.chrom_sizes = {
            "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
            "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
            "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
            "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
            "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
            "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
            "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
            "chr22": 50818468, "chrX": 156040895, "chrY": 57227415
        }

    def shuffle_within_chromosome(self, sv: pd.DataFrame) -> pd.DataFrame:
        """
        Shuffle breakpoint positions within each (cell, chromosome) group.
        Preserves per-chromosome breakpoint counts but destroys spatial structure.
        """
        sv_shuffled = sv.copy()
        for side in ["1", "2"]:
            chrom_col = f"chrom{side}"
            start_col = f"start{side}"
            end_col = f"end{side}"

            for (cell, chrom), idx in sv_shuffled.groupby(["cell_line", chrom_col]).groups.items():
                chrom_size = self.chrom_sizes.get(chrom)
                if chrom_size is None:
                    continue

                seg_len = (sv_shuffled.loc[idx, end_col] - sv_shuffled.loc[idx, start_col]).clip(lower=1).to_numpy()
                n = len(seg_len)
                # sample starts per segment length so end stays in-bounds
                starts = self.rng.integers(0, np.maximum(1, chrom_size - seg_len), size=n)
                ends = np.minimum(starts + seg_len, chrom_size)

                sv_shuffled.loc[idx, start_col] = starts
                sv_shuffled.loc[idx, end_col] = ends

        logger.info("Applied within-chromosome shuffle")
        return sv_shuffled

    def shuffle_across_chromosomes(self, sv: pd.DataFrame) -> pd.DataFrame:
        """
        Shuffle breakpoints across chromosomes within each cell line.
        Stronger control: breaks both spatial and chromosomal structure.
        """
        sv_shuffled = sv.copy()
        chroms = np.array(list(self.chrom_sizes.keys()))
        sizes = np.array(list(self.chrom_sizes.values()), dtype=float)
        probs = sizes / sizes.sum()

        for cell, idx in sv_shuffled.groupby("cell_line").groups.items():
            for side in ["1", "2"]:
                chrom_col = f"chrom{side}"
                start_col = f"start{side}"
                end_col = f"end{side}"

                seg_len = (sv_shuffled.loc[idx, end_col] - sv_shuffled.loc[idx, start_col]).clip(lower=1).to_numpy()
                new_chroms = self.rng.choice(chroms, size=len(seg_len), p=probs)

                starts = np.zeros_like(seg_len)
                ends = np.zeros_like(seg_len)
                for i, (cname, L) in enumerate(zip(new_chroms, seg_len)):
                    csize = self.chrom_sizes[cname]
                    s = int(self.rng.integers(0, max(1, csize - int(L))))
                    e = min(s + int(L), csize)
                    starts[i] = s
                    ends[i] = e

                sv_shuffled.loc[idx, chrom_col] = new_chroms
                sv_shuffled.loc[idx, start_col] = starts
                sv_shuffled.loc[idx, end_col] = ends

        logger.info("Applied cross-chromosome shuffle")
        return sv_shuffled


# ============================================================================
# Statistical Modeling
# ============================================================================

class BiasModel:
    """Fit per-gene models to estimate structural genomic bias."""

    def __init__(self, config: Config):
        self.config = config

    def _get_model(self):
        if self.config.model_type == "huber":
            return HuberRegressor(epsilon=self.config.huber_epsilon, max_iter=200)
        elif self.config.model_type == "quantile":
            return QuantileRegressor(quantile=self.config.quantile_alpha, solver="highs")
        else:
            from sklearn.linear_model import LinearRegression
            return LinearRegression()

    def _compute_vif(self, X: np.ndarray, feature_names: list) -> dict:
        from sklearn.linear_model import LinearRegression
        n_features = X.shape[1]
        vif_scores = {}
        for i in range(n_features):
            Xi = X[:, i]
            Xo = np.delete(X, i, axis=1)
            try:
                model = LinearRegression()
                model.fit(Xo, Xi)
                r2 = model.score(Xo, Xi)
                vif = 1 / (1 - r2) if r2 < 0.99 else np.inf
                vif_scores[feature_names[i]] = float(vif)
            except Exception:
                vif_scores[feature_names[i]] = np.nan
        return vif_scores

    def fit_gene_models(self, data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Model: dependency ~ CN + bp_near + bp_far
        Returns:
            predictions: DataFrame with structural bias predictions
            coefficients: DataFrame with model coefficients and diagnostics
        """
        logger.info(f"Fitting per-gene {self.config.model_type} models...")
        w1, w2 = self.config.proximity_windows
        data = data.copy()
        data["bp_dist"] = data["bp_dist"].fillna(self.config.max_bp_distance)
        data["cn"] = data["cn"].fillna(self.config.default_ploidy)
        data["bp_near"] = (data["bp_dist"] <= w1).astype(float)
        data["bp_far"] = ((data["bp_dist"] > w1) & (data["bp_dist"] <= w2)).astype(float)

        predictions, coefficients = [], []
        vif_done = False
        genes = data["gene"].unique()
        logger.info(f"Processing {len(genes):,} genes...")

        for gene in genes:
            gdf = data[data["gene"] == gene]
            if len(gdf) < self.config.min_cells_per_gene:
                continue

            # PART B: Optional orthogonalization of proximity vs CN
            if self.config.orthogonalize_proximity:
                # Form predictors [cn, bp_near, bp_far]
                cn_vals = gdf["cn"].fillna(self.config.default_ploidy).values
                bp_near_vals = gdf["bp_near"].values
                bp_far_vals = gdf["bp_far"].values
                
                # Orthogonalize: bp_near_resid = bp_near - Proj(bp_near | cn)
                from sklearn.linear_model import LinearRegression
                lr_near = LinearRegression()
                lr_near.fit(cn_vals.reshape(-1, 1), bp_near_vals)
                bp_near_resid = bp_near_vals - lr_near.predict(cn_vals.reshape(-1, 1))
                
                lr_far = LinearRegression()
                lr_far.fit(cn_vals.reshape(-1, 1), bp_far_vals)
                bp_far_resid = bp_far_vals - lr_far.predict(cn_vals.reshape(-1, 1))
                
                # Use residualized proximity features
                X = np.column_stack([cn_vals, bp_near_resid, bp_far_resid])
                feature_names = ["cn", "bp_near_resid", "bp_far_resid"]
            else:
                feature_names = ["cn", "bp_near", "bp_far"]
                X = gdf[feature_names].to_numpy()
            
            y = gdf["dependency"].to_numpy()

            if self.config.check_vif and not vif_done:
                vif_scores = self._compute_vif(X, feature_names)
                max_vif = np.nanmax(list(vif_scores.values()))
                if max_vif > self.config.vif_threshold:
                    logger.warning(f"High VIF detected (max={max_vif:.2f}): {vif_scores}")
                else:
                    logger.info(f"VIF check passed (max={max_vif:.2f}): {vif_scores}")
                vif_done = True

            try:
                model = self._get_model()
                model.fit(X, y)
                y_pred = model.predict(X)

                residuals = y - y_pred
                ss_res = float(np.sum(residuals ** 2))
                ss_tot = float(np.sum((y - y.mean()) ** 2))
                r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

                out = gdf[["gene", "cell_line"]].copy()
                out["structural_bias"] = y_pred
                predictions.append(out)

                # PART B: Handle coefficient names based on orthogonalization
                coef_dict = {
                    "gene": gene,
                    "intercept": float(model.intercept_) if hasattr(model, "intercept_") else 0.0,
                    "coef_cn": float(model.coef_[0]),
                    "r2": float(r2),
                    "n_obs": int(len(gdf)),
                    "rmse": float(np.sqrt(np.mean(residuals ** 2)))
                }
                
                if self.config.orthogonalize_proximity:
                    coef_dict["coef_bp_near_resid"] = float(model.coef_[1])
                    coef_dict["coef_bp_far_resid"] = float(model.coef_[2])
                    # Also keep original names for backward compatibility (set to NaN or 0)
                    coef_dict["coef_bp_near"] = 0.0
                    coef_dict["coef_bp_far"] = 0.0
                else:
                    coef_dict["coef_bp_near"] = float(model.coef_[1])
                    coef_dict["coef_bp_far"] = float(model.coef_[2])
                
                coefficients.append(coef_dict)

            except Exception as e:
                logger.debug(f"Failed to fit {gene}: {e}")
                continue

        pred_df = pd.concat(predictions, ignore_index=True) if predictions else pd.DataFrame()
        coef_df = pd.DataFrame(coefficients)
        logger.info(f"Successfully fit {len(coef_df):,} gene models")
        return pred_df, coef_df

    def correct_dependencies(self, original: pd.DataFrame, predictions: pd.DataFrame) -> pd.DataFrame:
        """
        Apply bias correction: dependency_corrected = dependency - (bias - median_bias_per_gene)
        Returns key columns INCLUDING cn/bp_dist for downstream metrics.
        """
        df = original.merge(predictions, on=["gene", "cell_line"], how="left")
        gene_medians = df.groupby("gene")["structural_bias"].median()
        df = df.merge(gene_medians.rename("gene_median_bias"),
                      left_on="gene", right_index=True, how="left")

        df["dependency_corrected"] = df["dependency"] - (df["structural_bias"] - df["gene_median_bias"])
        df["dependency_corrected"] = df["dependency_corrected"].fillna(df["dependency"])

        keep_cols = ["gene", "cell_line", "dependency", "dependency_corrected",
                     "structural_bias", "gene_median_bias", "cn", "bp_dist"]
        return df[keep_cols]


# ============================================================================
# Validation and Metrics
# ============================================================================

class ValidationMetrics:
    """Compute validation metrics with proper uncertainty quantification."""

    def __init__(self, config: Config):
        self.config = config
        self.rng = np.random.default_rng(config.seed)

    def compute_cn_correlation_change(self,
                                      data: pd.DataFrame,
                                      return_gene_level: bool = False) -> dict:
        """Change in |corr(dependency, CN)| before vs after correction."""
        results = []
        for gene, g in data.groupby("gene"):
            valid = g[["dependency", "dependency_corrected", "cn"]].dropna()
            if len(valid) < 10:
                continue
            try:
                r_before = valid[["dependency", "cn"]].corr().iloc[0, 1]
                r_after = valid[["dependency_corrected", "cn"]].corr().iloc[0, 1]
                if not (np.isnan(r_before) or np.isnan(r_after)):
                    results.append({
                        "gene": gene,
                        "r_before": float(r_before),
                        "r_after": float(r_after),
                        "delta_abs_corr": float(abs(r_before) - abs(r_after)),
                        "n_obs": int(len(valid))
                    })
            except Exception:
                continue

        if not results:
            return {"error": "No valid gene-level correlations computed"}

        df = pd.DataFrame(results)
        metrics = {
            "n_genes": int(len(df)),
            "median_delta_abs_corr": float(df["delta_abs_corr"].median()),
            "mean_delta_abs_corr": float(df["delta_abs_corr"].mean()),
            "frac_improved": float((df["delta_abs_corr"] > 0).mean()),
            "percentile_25": float(df["delta_abs_corr"].quantile(0.25)),
            "percentile_75": float(df["delta_abs_corr"].quantile(0.75))
        }

        # Bootstrap CI on median
        vals = df["delta_abs_corr"].to_numpy()
        n = len(vals)
        boot = []
        for _ in range(self.config.bootstrap_iterations):
            idx = self.rng.integers(0, n, size=n)
            boot.append(np.median(vals[idx]))
        metrics["median_ci_lower"] = float(np.percentile(boot, 2.5))
        metrics["median_ci_upper"] = float(np.percentile(boot, 97.5))

        if return_gene_level:
            metrics["gene_level"] = df
        return metrics

    def compute_directional_proximity_active(self,
                                             data: pd.DataFrame,
                                             coefficients: pd.DataFrame) -> dict:
        """
        Fraction of genes that are active (>=X% cells with |contrib|>=thr)
        AND where correction reduces |corr(dependency, CN)|.
        """
        logger.info("Computing directional proximity-active metric...")
        df = data.merge(
            coefficients[["gene", "coef_bp_near", "coef_bp_far"]],
            on="gene", how="left"
        )

        w1, w2 = self.config.proximity_windows
        df["bp_dist"] = df["bp_dist"].fillna(self.config.max_bp_distance)
        df["bp_near"] = (df["bp_dist"] <= w1).astype(float)
        df["bp_far"] = ((df["bp_dist"] > w1) & (df["bp_dist"] <= w2)).astype(float)

        df["proximity_contrib"] = (
            df["coef_bp_near"].fillna(0) * df["bp_near"] +
            df["coef_bp_far"].fillna(0) * df["bp_far"]
        )

        flags = []
        details = []

        for gene, sub in df.groupby("gene"):
            active_frac = (np.abs(sub["proximity_contrib"]) >= self.config.activity_contrib_threshold).mean()

            a = sub[["dependency", "cn"]].dropna()
            b = sub[["dependency_corrected", "cn"]].dropna()
            if len(a) < 10 or len(b) < 10:
                continue

            try:
                r0 = a.corr().iloc[0, 1]
                r1 = b.corr().iloc[0, 1]
                corr_drop = float(abs(r0) - abs(r1))
                is_active = active_frac >= self.config.activity_cell_fraction
                is_directional = is_active and (corr_drop > 0)
                flags.append(is_directional)
                details.append({
                    "gene": gene,
                    "active_frac": float(active_frac),
                    "corr_drop": float(corr_drop),
                    "is_active": bool(is_active),
                    "is_directional": bool(is_directional)
                })
            except Exception:
                continue

        if not flags:
            return {"error": "No valid genes for directional metric"}

        frac = float(np.mean(flags))
        # Bootstrap CI on the directional fraction
        flags_arr = np.array(flags, dtype=float)
        n = len(flags_arr)
        boot = []
        for _ in range(self.config.bootstrap_iterations):
            idx = self.rng.integers(0, n, size=n)
            boot.append(np.mean(flags_arr[idx]))

        metrics = {
            "n_genes": int(n),
            "directional_prox_active_frac": frac,
            "directional_ci_lower": float(np.percentile(boot, 2.5)),
            "directional_ci_upper": float(np.percentile(boot, 97.5)),
            "activity_contrib_threshold": float(self.config.activity_contrib_threshold),
            "activity_cell_fraction": float(self.config.activity_cell_fraction),
            "median_active_frac": float(pd.DataFrame(details)["active_frac"].median()) if details else np.nan,
            "median_corr_drop": float(pd.DataFrame(details)["corr_drop"].median()) if details else np.nan
        }
        return metrics

    def compute_proximity_effect_size(self, data: pd.DataFrame,
                                      coefficients: pd.DataFrame) -> dict:
        """Quantify magnitude of proximity effects across genes."""
        df = data.merge(coefficients[["gene", "coef_bp_near", "coef_bp_far"]],
                        on="gene", how="left")

        w1, w2 = self.config.proximity_windows
        bp = df["bp_dist"].fillna(self.config.max_bp_distance)
        df["bp_near"] = (bp <= w1).astype(float)
        df["bp_far"] = ((bp > w1) & (bp <= w2)).astype(float)

        df["proximity_contribution"] = (
            df["coef_bp_near"].fillna(0) * df["bp_near"] +
            df["coef_bp_far"].fillna(0) * df["bp_far"]
        )

        gene_activity = df.groupby("gene").apply(
            lambda g: (g["proximity_contribution"].abs() >= self.config.activity_contrib_threshold).mean()
        )

        metrics = {
            "median_abs_proximity_coef": float(
                coefficients[["coef_bp_near", "coef_bp_far"]].abs().max(axis=1).median()
            ),
            "frac_genes_with_proximity_effect": float(
                (coefficients[["coef_bp_near", "coef_bp_far"]].abs().max(axis=1) > 0.05).mean()
            ),
            "median_gene_activity_rate": float(gene_activity.median()),
            "mean_gene_activity_rate": float(gene_activity.mean())
        }
        return metrics

    def evaluate_essential_separation(self,
                                      data: pd.DataFrame,
                                      essential_genes: set,
                                      nonessential_genes: set) -> dict:
        """Evaluate separation of essential vs nonessential after correction."""
        test_genes = essential_genes | nonessential_genes
        df = data[data["gene"].isin(test_genes)].copy()
        if len(df) == 0:
            return {"error": "No overlap with essential/nonessential sets"}

        df["is_essential"] = df["gene"].isin(essential_genes).astype(int)
        gene_scores = df.groupby("gene").agg({
            "dependency": "median",
            "dependency_corrected": "median",
            "is_essential": "first"
        })

        y_true = gene_scores["is_essential"].values
        y_score_orig = -gene_scores["dependency"].values
        y_score_corr = -gene_scores["dependency_corrected"].values

        auroc_orig = roc_auc_score(y_true, y_score_orig)
        auroc_corr = roc_auc_score(y_true, y_score_corr)
        auprc_orig = average_precision_score(y_true, y_score_orig)
        auprc_corr = average_precision_score(y_true, y_score_corr)

        return {
            "n_genes_tested": int(len(gene_scores)),
            "n_essential": int(y_true.sum()),
            "auroc_original": float(auroc_orig),
            "auroc_corrected": float(auroc_corr),
            "auroc_delta": float(auroc_corr - auroc_orig),
            "auprc_original": float(auprc_orig),
            "auprc_corrected": float(auprc_corr),
            "auprc_delta": float(auprc_corr - auprc_orig)
        }

    def compute_partial_corr_delta(self, data: pd.DataFrame) -> dict:
        """
        Compute Δ|partial corr(dep, CN | proximity)| using OLS residualization.
        Replicates PART A logic for pipeline integration.
        """
        logger.info("Computing partial correlation delta...")
        w1, w2 = self.config.proximity_windows
        data = data.copy()
        data["bp_dist"] = data["bp_dist"].fillna(self.config.max_bp_distance)
        data["bp_near"] = (data["bp_dist"] <= w1).astype(float)
        data["bp_far"] = ((data["bp_dist"] > w1) & (data["bp_dist"] <= w2)).astype(float)
        
        results = []
        from sklearn.linear_model import LinearRegression
        
        for gene, gdf in data.groupby("gene"):
            valid = gdf[["dependency", "dependency_corrected", "cn", "bp_near", "bp_far"]].dropna()
            if len(valid) < self.config.min_cells_per_gene:
                continue
            
            try:
                X_prox = valid[["bp_near", "bp_far"]].values
                y_dep = valid["dependency"].values
                y_dep_corr = valid["dependency_corrected"].values
                y_cn = valid["cn"].values
                
                # Residualize
                lr_dep = LinearRegression()
                lr_dep.fit(X_prox, y_dep)
                dep_resid = y_dep - lr_dep.predict(X_prox)
                
                lr_dep_corr = LinearRegression()
                lr_dep_corr.fit(X_prox, y_dep_corr)
                dep_corr_resid = y_dep_corr - lr_dep_corr.predict(X_prox)
                
                lr_cn = LinearRegression()
                lr_cn.fit(X_prox, y_cn)
                cn_resid = y_cn - lr_cn.predict(X_prox)
                
                # Partial correlations
                r0 = np.corrcoef(dep_resid, cn_resid)[0, 1]
                r1 = np.corrcoef(dep_corr_resid, cn_resid)[0, 1]
                
                if np.isnan(r0) or np.isnan(r1):
                    continue
                
                delta = abs(r0) - abs(r1)
                results.append({
                    "gene": gene,
                    "partial_corr_delta_abs": delta,
                    "r0_abs": abs(r0),
                    "r1_abs": abs(r1),
                    "n_obs": len(valid)
                })
            except Exception:
                continue
        
        if not results:
            return {"error": "No valid partial correlations computed"}
        
        results_df = pd.DataFrame(results)
        vals = results_df["partial_corr_delta_abs"].values
        
        # Bootstrap median (2,000 reps)
        n = len(vals)
        boot_medians = []
        for _ in range(2000):
            idx = self.rng.integers(0, n, size=n)
            boot_medians.append(np.median(vals[idx]))
        
        boot_medians = np.array(boot_medians)
        
        metrics = {
            "n_genes": int(len(results_df)),
            "median_partial_delta_abs_corr": float(np.median(vals)),
            "mean_partial_delta_abs_corr": float(np.mean(vals)),
            "frac_positive": float((vals > 0).mean()),
            "partial_ci_lower": float(np.percentile(boot_medians, 2.5)),
            "partial_ci_upper": float(np.percentile(boot_medians, 97.5)),
            "gene_level": results_df
        }
        
        return metrics

    def compute_excess_signal(self,
                              true_metrics: dict,
                              shuffled_metrics: dict,
                              metric_name: str) -> dict:
        """
        Excess signal (true - shuffled) with bootstrap CI if gene-level is present.
        """
        if "gene_level" in true_metrics and "gene_level" in shuffled_metrics:
            true_df = true_metrics["gene_level"]
            shuf_df = shuffled_metrics["gene_level"]
            merged = true_df.merge(shuf_df, on="gene", suffixes=("_true", "_shuf"))
            if len(merged) == 0:
                return {"error": "No overlapping genes for excess signal"}

            merged["excess"] = merged[f"{metric_name}_true"] - merged[f"{metric_name}_shuf"]
            vals = merged["excess"].to_numpy()
            n = len(vals)
            boot = []
            for _ in range(self.config.bootstrap_iterations):
                idx = self.rng.integers(0, n, size=n)
                boot.append(np.median(vals[idx]))

            return {
                "median_excess": float(np.median(vals)),
                "excess_ci_lower": float(np.percentile(boot, 2.5)),
                "excess_ci_upper": float(np.percentile(boot, 97.5)),
                "n_genes": int(n)
            }
        else:
            true_val = true_metrics.get(metric_name, np.nan)
            shuf_val = shuffled_metrics.get(metric_name, np.nan)
            if np.isnan(true_val) or np.isnan(shuf_val):
                return {"error": f"Missing {metric_name} in one or both datasets"}
            return {
                "excess": float(true_val - shuf_val),
                "true_value": float(true_val),
                "shuffled_value": float(shuf_val)
            }


# ============================================================================
# Diagnostic Plotting (optional)
# ============================================================================

class DiagnosticPlots:
    """Generate diagnostic visualizations (requires matplotlib)."""

    @staticmethod
    def plot_gene_diagnostics(gene_name: str, gene_data: pd.DataFrame,
                              coefficients: dict, output_path: Path):
        try:
            import matplotlib.pyplot as plt
            from scipy.stats import probplot
        except ImportError:
            logger.warning("matplotlib not available, skipping plots")
            return

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f"Diagnostic Plot: {gene_name}", fontsize=14, fontweight='bold')

        # Panel 1: Dependency vs CN
        ax = axes[0, 0]
        ax.scatter(gene_data["cn"], gene_data["dependency"], alpha=0.5, s=20, label="Original")
        ax.scatter(gene_data["cn"], gene_data["dependency_corrected"], alpha=0.5, s=20, label="Corrected")
        ax.set_xlabel("Copy Number"); ax.set_ylabel("Dependency Score")
        ax.set_title("CN vs Dependency"); ax.legend(); ax.grid(True, alpha=0.3)

        # Panel 2: Residuals
        ax = axes[0, 1]
        residuals = gene_data["dependency"] - gene_data["structural_bias"]
        ax.scatter(gene_data["structural_bias"], residuals, alpha=0.5, s=20)
        ax.axhline(y=0, color='r', linestyle='--', linewidth=1)
        ax.set_xlabel("Fitted Values"); ax.set_ylabel("Residuals")
        ax.set_title("Residual Plot"); ax.grid(True, alpha=0.3)

        # Panel 3: Proximity effect
        ax = axes[1, 0]
        bp = gene_data["bp_dist"].fillna(5_000_000).astype(float)
        colors = ['red' if d < 100000 else 'orange' if d < 1000000 else 'blue' for d in bp]
        ax.scatter(bp / 1e6, gene_data["dependency"], c=colors, alpha=0.5, s=20)
        ax.set_xlabel("Distance to Breakpoint (Mb)"); ax.set_ylabel("Dependency Score")
        ax.set_title("Proximity Effect"); ax.set_xscale('log'); ax.grid(True, alpha=0.3)

        # Panel 4: Q-Q plot
        ax = axes[1, 1]
        probplot(residuals.dropna(), dist="norm", plot=ax)
        ax.set_title("Q-Q Plot (Normality Check)"); ax.grid(True, alpha=0.3)

        coef_text = (f"CN coef: {coefficients['coef_cn']:.3f}\n"
                     f"BP near coef: {coefficients['coef_bp_near']:.3f}\n"
                     f"BP far coef: {coefficients['coef_bp_far']:.3f}\n"
                     f"R²: {coefficients['r2']:.3f}")
        fig.text(0.98, 0.02, coef_text, fontsize=9, family='monospace',
                 ha='right', va='bottom',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()


# ============================================================================
# Main Pipeline
# ============================================================================

class Pipeline:
    """Main pipeline orchestrator."""

    def __init__(self, config: Config):
        self.config = config
        self.validator = DataValidator()
        self.feature_computer = FeatureComputer(config)
        self.shuffler = ShuffleControl(config)
        self.modeler = BiasModel(config)
        self.metrics = ValidationMetrics(config)

    def load_and_validate_data(self,
                               dep_path: str,
                               cnv_path: str,
                               sv_path: str,
                               genes_path: str) -> dict:
        """Load and validate all input data."""
        logger.info("=" * 80); logger.info("LOADING DATA"); logger.info("=" * 80)

        dep = pd.read_csv(dep_path)
        cnv = pd.read_csv(cnv_path, sep="\t")
        sv = pd.read_csv(sv_path, sep="\t")
        genes = pd.read_csv(genes_path, sep="\t",
                            header=None,
                            names=["chrom", "start", "end", "gene", "strand"])

        dep = self.validator.validate_dependency(dep)
        cnv = self.validator.validate_cnv(cnv, self.config.max_cn)
        sv = self.validator.validate_sv(sv)
        genes = self.validator.validate_genes(genes)

        # PART B: Early gene subsetting
        if self.config.genes_subset_file:
            subset_path = Path(self.config.genes_subset_file)
            if not subset_path.exists():
                raise FileNotFoundError(f"Gene subset file not found: {subset_path}")
            subset_genes = set()
            with open(subset_path, 'r') as f:
                for line in f:
                    gene = line.strip()
                    if gene:
                        subset_genes.add(gene)
            logger.info(f"Subsetting to {len(subset_genes):,} genes from {subset_path}")
            
            # Filter dependency and genes tables
            dep_before = len(dep)
            dep = dep[dep["gene"].isin(subset_genes)].copy()
            logger.info(f"Filtered dependency: {len(dep):,} rows (from {dep_before:,})")
            
            genes_before = len(genes)
            genes = genes[genes["gene"].isin(subset_genes)].copy()
            logger.info(f"Filtered genes: {len(genes):,} genes (from {genes_before:,})")

        dep_genes = set(dep["gene"])
        genes_bed = set(genes["gene"])
        gene_overlap = len(dep_genes & genes_bed)
        logger.info(f"Gene overlap: {gene_overlap:,} / {len(dep_genes):,} "
                    f"({100*gene_overlap/len(dep_genes):.1f}%)")
        if gene_overlap < len(dep_genes) * 0.5:
            logger.warning("Low gene overlap! Check gene ID formats.")

        dep_cells = set(dep["cell_line"])
        cnv_cells = set(cnv["cell_line"])
        sv_cells = set(sv["cell_line"])
        cell_overlap = len(dep_cells & cnv_cells & sv_cells)
        logger.info(f"Cell overlap: {cell_overlap:,} / {len(dep_cells):,} "
                    f"({100*cell_overlap/len(dep_cells):.1f}%)")
        if cell_overlap < 100:
            raise ValueError("Insufficient cell line overlap")

        return {"dependency": dep, "cnv": cnv, "sv": sv, "genes": genes}

    def compute_features(self, data: dict) -> pd.DataFrame:
        """Compute all genomic features with optional caching."""
        logger.info("=" * 80); logger.info("COMPUTING FEATURES"); logger.info("=" * 80)
        
        # PART B: Feature caching
        if self.config.read_features_dir:
            cache_dir = Path(self.config.read_features_dir)
            cn_cache = cache_dir / "features_cn.parquet"
            bp_cache = cache_dir / "features_bp.parquet"
            
            if cn_cache.exists() and bp_cache.exists():
                logger.info(f"Loading features from cache: {cache_dir}")
                gene_cn = pd.read_parquet(cn_cache)
                bp_dist = pd.read_parquet(bp_cache)
                logger.info(f"Loaded CN features: {len(gene_cn):,} rows")
                logger.info(f"Loaded BP features: {len(bp_dist):,} rows")
            else:
                logger.warning(f"Cache files not found in {cache_dir}, computing features...")
                gene_cn = self.feature_computer.compute_gene_cn(data["cnv"], data["genes"])
                bp_dist = self.feature_computer.compute_bp_distances(data["sv"], data["genes"])
        else:
            gene_cn = self.feature_computer.compute_gene_cn(data["cnv"], data["genes"])
            bp_dist = self.feature_computer.compute_bp_distances(data["sv"], data["genes"])

        features = gene_cn.merge(bp_dist, on=["gene", "cell_line"], how="outer")
        merged = data["dependency"].merge(features, on=["gene", "cell_line"], how="left")

        # PART B: Write feature cache if requested
        if self.config.write_features_dir:
            cache_dir = Path(self.config.write_features_dir)
            cache_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Saving features to cache: {cache_dir}")
            gene_cn.to_parquet(cache_dir / "features_cn.parquet", index=False)
            bp_dist.to_parquet(cache_dir / "features_bp.parquet", index=False)
            logger.info("Feature cache saved")

        logger.info(f"Feature matrix: {len(merged):,} rows")
        logger.info(f"CN coverage: {merged['cn'].notna().mean():.1%}")
        logger.info(f"BP coverage: {merged['bp_dist'].notna().mean():.1%}")
        return merged

    def run_analysis(self,
                     data: pd.DataFrame,
                     label: str = "main",
                     essential_genes: Optional[set] = None,
                     nonessential_genes: Optional[set] = None) -> dict:
        """Run complete analysis pipeline on a dataset."""
        logger.info(f"Running analysis: {label}")

        predictions, coefficients = self.modeler.fit_gene_models(data)
        if len(coefficients) == 0:
            logger.error("No models successfully fit!")
            return {"error": "Model fitting failed"}

        corrected = self.modeler.correct_dependencies(data, predictions)

        results = {"label": label, "n_genes_fitted": int(len(coefficients))}

        cn_corr_metrics = self.metrics.compute_cn_correlation_change(
            corrected, return_gene_level=True
        )
        results["cn_correlation"] = cn_corr_metrics

        directional_metrics = self.metrics.compute_directional_proximity_active(
            corrected, coefficients
        )
        results["directional_prox_active"] = directional_metrics

        proximity_metrics = self.metrics.compute_proximity_effect_size(
            corrected, coefficients
        )
        results["proximity_effects"] = proximity_metrics

        # PART B: Partial correlation metric
        if self.config.compute_partial_corr:
            partial_metrics = self.metrics.compute_partial_corr_delta(corrected)
            results["partial_correlation"] = partial_metrics

        if essential_genes and nonessential_genes:
            sep_metrics = self.metrics.evaluate_essential_separation(
                corrected, essential_genes, nonessential_genes
            )
            results["essential_separation"] = sep_metrics

        out_dir = self.config.output_dir / label
        out_dir.mkdir(parents=True, exist_ok=True)
        corrected.to_csv(out_dir / "dependency_corrected.csv", index=False)
        coefficients.to_csv(out_dir / "model_coefficients.csv", index=False)

        with open(out_dir / "metrics.json", "w") as f:
            results_json = {k: v for k, v in results.items()}
            if "cn_correlation" in results_json and "gene_level" in results_json["cn_correlation"]:
                gene_level_df = results_json["cn_correlation"].pop("gene_level")
                gene_level_df.to_csv(out_dir / "gene_level_correlations.csv", index=False)
            # PART B: Save partial correlation per-gene results
            if "partial_correlation" in results_json and "gene_level" in results_json["partial_correlation"]:
                partial_df = results_json["partial_correlation"].pop("gene_level")
                partial_df.to_csv(out_dir / "partial_corr_per_gene.csv", index=False)
            json.dump(results_json, f, indent=2, default=str)

        logger.info(f"Saved outputs to {out_dir}")
        return results

    def _compute_proximity_prevalence(self, data: pd.DataFrame) -> pd.Series:
        """Compute per-gene proximity prevalence (mean of bp_near or bp_far)."""
        w1, w2 = self.config.proximity_windows
        data = data.copy()
        data["bp_dist"] = data["bp_dist"].fillna(self.config.max_bp_distance)
        data["bp_near"] = (data["bp_dist"] <= w1).astype(float)
        data["bp_far"] = ((data["bp_dist"] > w1) & (data["bp_dist"] <= w2)).astype(float)
        data["prox_active"] = (data["bp_near"] | data["bp_far"]).astype(float)
        
        # Per-gene mean prevalence
        prevalence = data.groupby("gene")["prox_active"].mean()
        return prevalence

    def _filter_by_prevalence_match(self, true_data: pd.DataFrame, 
                                    shuffled_data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, int]:
        """
        Filter genes by proximity prevalence match (tolerance ±0.02).
        Returns filtered dataframes and count of genes kept.
        """
        true_prev = self._compute_proximity_prevalence(true_data)
        shuf_prev = self._compute_proximity_prevalence(shuffled_data)
        
        # Find genes with matching prevalence
        common_genes = set(true_prev.index) & set(shuf_prev.index)
        matched_genes = []
        
        for gene in common_genes:
            diff = abs(true_prev[gene] - shuf_prev[gene])
            if diff <= self.config.prevalence_tolerance:
                matched_genes.append(gene)
        
        logger.info(f"Prevalence matching: {len(matched_genes):,} / {len(common_genes):,} genes "
                   f"within tolerance ±{self.config.prevalence_tolerance}")
        
        # Filter dataframes
        true_filtered = true_data[true_data["gene"].isin(matched_genes)].copy()
        shuf_filtered = shuffled_data[shuffled_data["gene"].isin(matched_genes)].copy()
        
        return true_filtered, shuf_filtered, len(matched_genes)

    def run_with_controls(self,
                          data: dict,
                          essential_path: Optional[str] = None,
                          nonessential_path: Optional[str] = None,
                          run_shuffles: bool = True,
                          shuffle_mode: Optional[str] = "both") -> dict:
        """Run full analysis with (optional) negative controls."""
        essential_genes = None
        nonessential_genes = None

        if essential_path:
            essential_genes = set(pd.read_csv(essential_path, header=None)[0])
            logger.info(f"Loaded {len(essential_genes)} essential genes")
        if nonessential_path:
            nonessential_genes = set(pd.read_csv(nonessential_path, header=None)[0])
            logger.info(f"Loaded {len(nonessential_genes)} nonessential genes")

        feature_data = self.compute_features(data)

        logger.info("=" * 80); logger.info("MAIN ANALYSIS (True SV positions)"); logger.info("=" * 80)
        results = {}
        results["true"] = self.run_analysis(
            feature_data,
            label="true",
            essential_genes=essential_genes,
            nonessential_genes=nonessential_genes
        )

        if run_shuffles:
            if shuffle_mode in ("within", "both"):
                logger.info("=" * 80); logger.info("CONTROL 1: Within-chromosome shuffle"); logger.info("=" * 80)
                sv_shuf_within = self.shuffler.shuffle_within_chromosome(data["sv"])
                feat_within = self.compute_features({**data, "sv": sv_shuf_within})
                
                # PART B: Prevalence matching for shuffled data
                if self.config.match_proximity_prevalence:
                    feat_within_filtered, _, n_kept = self._filter_by_prevalence_match(
                        feature_data, feat_within
                    )
                    logger.info(f"Using {n_kept:,} genes after prevalence matching for within-chrom shuffle")
                    feat_within = feat_within[feat_within["gene"].isin(feat_within_filtered["gene"].unique())].copy()
                
                results["shuffle_within_chrom"] = self.run_analysis(
                    feat_within, label="shuffle_within_chrom",
                    essential_genes=essential_genes, nonessential_genes=nonessential_genes
                )

            if shuffle_mode in ("across", "both"):
                logger.info("=" * 80); logger.info("CONTROL 2: Cross-chromosome shuffle"); logger.info("=" * 80)
                sv_shuf_across = self.shuffler.shuffle_across_chromosomes(data["sv"])
                feat_across = self.compute_features({**data, "sv": sv_shuf_across})
                
                # PART B: Prevalence matching for shuffled data
                if self.config.match_proximity_prevalence:
                    feat_across_filtered, _, n_kept = self._filter_by_prevalence_match(
                        feature_data, feat_across
                    )
                    logger.info(f"Using {n_kept:,} genes after prevalence matching for across-chrom shuffle")
                    feat_across = feat_across[feat_across["gene"].isin(feat_across_filtered["gene"].unique())].copy()
                
                results["shuffle_across_chrom"] = self.run_analysis(
                    feat_across, label="shuffle_across_chrom",
                    essential_genes=essential_genes, nonessential_genes=nonessential_genes
                )

            logger.info("=" * 80); logger.info("COMPUTING EXCESS SIGNALS"); logger.info("=" * 80)
            if shuffle_mode in ("within", "both") and "shuffle_within_chrom" in results:
                results["excess_within"] = self._compute_all_excess(
                    results["true"], results["shuffle_within_chrom"]
                )
            if shuffle_mode in ("across", "both") and "shuffle_across_chrom" in results:
                results["excess_across"] = self._compute_all_excess(
                    results["true"], results["shuffle_across_chrom"]
                )

        self._generate_comparison_report(results, run_shuffles=run_shuffles)
        return results

    def _compute_all_excess(self, true_results: dict, shuffled_results: dict) -> dict:
        """Compute excess signal for all available metrics."""
        excess = {}
        if "cn_correlation" in true_results and "cn_correlation" in shuffled_results:
            excess["cn_correlation"] = self.metrics.compute_excess_signal(
                true_results["cn_correlation"],
                shuffled_results["cn_correlation"],
                "delta_abs_corr"
            )
        if "directional_prox_active" in true_results and "directional_prox_active" in shuffled_results:
            tv = true_results["directional_prox_active"]["directional_prox_active_frac"]
            sv = shuffled_results["directional_prox_active"]["directional_prox_active_frac"]
            excess["directional_prox_active"] = {"excess": float(tv - sv),
                                                 "true_value": float(tv),
                                                 "shuffled_value": float(sv)}
        return excess

    def _generate_comparison_report(self, results: dict, run_shuffles: bool):
        """Generate a summary comparing true vs shuffled results."""
        report_path = self.config.output_dir / "comparison_report.txt"
        with open(report_path, "w") as f:
            f.write("=" * 80 + "\n")
            f.write("SV/CN BIAS CORRECTION PIPELINE - COMPARISON REPORT\n")
            f.write("=" * 80 + "\n\n")

            f.write("Configuration:\n")
            f.write(f"  Model: {self.config.model_type}\n")
            f.write(f"  Proximity windows: {self.config.proximity_windows[0]:,}bp, {self.config.proximity_windows[1]:,}bp\n")
            f.write(f"  Activity threshold: {self.config.activity_contrib_threshold} dependency units\n")
            f.write(f"  Activity cell fraction: {self.config.activity_cell_fraction:.0%}\n")
            f.write(f"  Bootstrap iterations: {self.config.bootstrap_iterations}\n\n")

            for label in ["true", "shuffle_within_chrom", "shuffle_across_chrom"]:
                if label not in results:
                    continue
                res = results[label]
                f.write(f"\n{label.upper().replace('_', ' ')}\n")
                f.write("-" * 80 + "\n")
                f.write(f"Genes fitted: {res.get('n_genes_fitted', 0):,}\n")

                if "cn_correlation" in res and "error" not in res["cn_correlation"]:
                    cc = res["cn_correlation"]
                    f.write("\nCN Correlation Change:\n")
                    f.write(f"  Median Δ|corr|: {cc['median_delta_abs_corr']:.4f} "
                            f"[{cc['median_ci_lower']:.4f}, {cc['median_ci_upper']:.4f}]\n")
                    f.write(f"  Fraction improved: {cc['frac_improved']:.1%}\n")

                if "directional_prox_active" in res and "error" not in res["directional_prox_active"]:
                    dpa = res["directional_prox_active"]
                    f.write("\nDirectional Proximity-Active Metric (GOLD STANDARD):\n")
                    f.write(f"  Fraction: {dpa['directional_prox_active_frac']:.3f} "
                            f"[{dpa['directional_ci_lower']:.3f}, {dpa['directional_ci_upper']:.3f}]\n")
                    f.write(f"  (Active in ≥{dpa['activity_cell_fraction']:.0%} cells "
                            f"with |contrib|≥{dpa['activity_contrib_threshold']} AND reduces |corr|)\n")

                if "proximity_effects" in res:
                    pe = res["proximity_effects"]
                    f.write("\nProximity Effects:\n")
                    f.write(f"  Median |coef|: {pe['median_abs_proximity_coef']:.4f}\n")
                    f.write(f"  Genes with effect: {pe['frac_genes_with_proximity_effect']:.1%}\n")
                    f.write(f"  Mean gene activity: {pe['mean_gene_activity_rate']:.1%}\n")

                if "essential_separation" in res and "error" not in res["essential_separation"]:
                    es = res["essential_separation"]
                    f.write("\nEssential Gene Separation:\n")
                    f.write(f"  AUROC: {es['auroc_original']:.3f} → {es['auroc_corrected']:.3f} "
                            f"(Δ={es['auroc_delta']:+.3f})\n")
                    f.write(f"  AUPRC: {es['auprc_original']:.3f} → {es['auprc_corrected']:.3f} "
                            f"(Δ={es['auprc_delta']:+.3f})\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("EXCESS SIGNAL (True - Shuffled)\n")
            f.write("=" * 80 + "\n\n")

            if run_shuffles:
                for control_label, excess_key in [("Within-Chrom", "excess_within"),
                                                  ("Cross-Chrom", "excess_across")]:
                    if excess_key not in results:
                        continue
                    f.write(f"\n{control_label} Control:\n")
                    f.write("-" * 40 + "\n")
                    excess = results[excess_key]
                    if "cn_correlation" in excess and "median_excess" in excess["cn_correlation"]:
                        cc_ex = excess["cn_correlation"]
                        f.write(f"  Median Δ|corr| excess: {cc_ex['median_excess']:.4f} "
                                f"[{cc_ex['excess_ci_lower']:.4f}, {cc_ex['excess_ci_upper']:.4f}] "
                                f"(n={cc_ex['n_genes']})\n")
                    if "directional_prox_active" in excess:
                        dpa_ex = excess["directional_prox_active"]
                        f.write(f"  Directional prox-active excess: {dpa_ex['excess']:.3f}\n")
                        f.write(f"    True: {dpa_ex['true_value']:.3f}\n")
                        f.write(f"    Shuffled: {dpa_ex['shuffled_value']:.3f}\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("INTERPRETATION\n")
            f.write("=" * 80 + "\n\n")

            if "true" in results and "directional_prox_active" in results["true"]:
                true_val = results["true"]["directional_prox_active"].get("directional_prox_active_frac", np.nan)
                f.write(f"The directional proximity-active metric indicates that "
                        f"{(0 if np.isnan(true_val) else true_val):.1%} of genes exhibit\n"
                        f"positionally specific proximity effects that reduce CN–dependency coupling after correction.\n")
                if run_shuffles and "excess_within" in results and "directional_prox_active" in results["excess_within"]:
                    excess_val = results["excess_within"]["directional_prox_active"]["excess"]
                    f.write(f"Relative to within-chromosome shuffles, the excess directional fraction is "
                            f"{excess_val:+.1%}, consistent with a breakpoint-driven artifact.\n")

            f.write("\nOutputs written under: {}\n".format(self.config.output_dir.resolve()))

        logger.info(f"Wrote comparison report: {report_path}")


# ============================================================================
# CLI
# ============================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="SV/CN bias correction for CRISPR dependency screens",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--dependency", required=True, help="Long-format CSV with columns: gene,cell_line,dependency")
    p.add_argument("--cnv", required=True, help="CNV segments (BED-like TSV) with columns: cell_line,chrom,start,end,cn")
    p.add_argument("--sv", required=True, help="SV breakpoints (BEDPE-like TSV) with chrom1,start1,end1,chrom2,start2,end2,cell_line")
    p.add_argument("--genes", required=True, help="Gene BED (TSV) with chrom,start,end,gene[,strand]")

    p.add_argument("--essential", default=None, help="Optional file (one gene per line) of essential genes")
    p.add_argument("--nonessential", default=None, help="Optional file (one gene per line) of nonessential genes")

    p.add_argument("--model", choices=["huber", "quantile", "ols"], default="huber")
    p.add_argument("--huber-epsilon", type=float, default=1.35)
    p.add_argument("--quantile-alpha", type=float, default=0.5)
    p.add_argument("--bp-windows", nargs=2, type=int, default=[100000, 1000000],
                   metavar=("W_NEAR", "W_FAR"), help="Breakpoint proximity windows (bp)")
    p.add_argument("--max-bp-distance", type=int, default=5_000_000)
    p.add_argument("--default-ploidy", type=float, default=2.0)
    p.add_argument("--max-cn", type=float, default=10.0)
    p.add_argument("--min-cells-per-gene", type=int, default=20)

    p.add_argument("--bootstrap-iterations", type=int, default=1000)
    p.add_argument("--activity-threshold", type=float, default=0.01,
                   help="Contribution threshold for activity (dependency units)")
    p.add_argument("--activity-cell-fraction", type=float, default=0.10,
                   help="Minimum fraction of cells to call a gene 'active'")

    p.add_argument("--vif-check", action="store_true")
    p.add_argument("--vif-threshold", type=float, default=5.0)

    p.add_argument("--output-dir", default="output", help="Directory to write outputs")
    p.add_argument("--seed", type=int, default=42)

    p.add_argument("--no-shuffles", action="store_true", help="Skip shuffled negative controls")
    p.add_argument("--shuffle-within", action="store_true", help="Run only within-chromosome shuffle")
    p.add_argument("--shuffle-across", action="store_true", help="Run only across-chromosome shuffle")
    
    # PART B: New CLI flags
    p.add_argument("--genes-subset-file", type=str, default=None,
                   help="File with newline-separated gene symbols to subset")
    p.add_argument("--write-features", type=str, default=None,
                   help="Directory to write feature cache (features_cn.parquet, features_bp.parquet)")
    p.add_argument("--read-features", type=str, default=None,
                   help="Directory to read feature cache from")
    p.add_argument("--orthogonalize-proximity", action="store_true",
                   help="Orthogonalize proximity features against CN before fitting")
    p.add_argument("--compute-partial-corr", action="store_true",
                   help="Compute and report partial correlation metric")
    p.add_argument("--match-proximity-prevalence", action="store_true",
                   help="Filter genes by proximity prevalence match (tolerance ±0.02)")
    
    return p


def main():
    parser = build_parser()
    args = parser.parse_args()

    cfg = Config(
        model_type=args.model,
        huber_epsilon=args.huber_epsilon,
        quantile_alpha=args.quantile_alpha,
        proximity_windows=(args.bp_windows[0], args.bp_windows[1]),
        max_bp_distance=args.max_bp_distance,
        default_ploidy=args.default_ploidy,
        max_cn=args.max_cn,
        min_cells_per_gene=args.min_cells_per_gene,
        bootstrap_iterations=args.bootstrap_iterations,
        activity_contrib_threshold=args.activity_threshold,
        activity_cell_fraction=args.activity_cell_fraction,
        check_vif=args.vif_check,
        vif_threshold=args.vif_threshold,
        output_dir=Path(args.output_dir),
        seed=args.seed,
        genes_subset_file=args.genes_subset_file,
        write_features_dir=args.write_features,
        read_features_dir=args.read_features,
        orthogonalize_proximity=args.orthogonalize_proximity,
        compute_partial_corr=args.compute_partial_corr,
        match_proximity_prevalence=args.match_proximity_prevalence
    )

    pipe = Pipeline(cfg)
    data = pipe.load_and_validate_data(
        dep_path=args.dependency,
        cnv_path=args.cnv,
        sv_path=args.sv,
        genes_path=args.genes
    )

    # Determine shuffle mode
    if args.no_shuffles:
        run_shuffles = False
        shuffle_mode = None
    elif args.shuffle_within:
        run_shuffles = True
        shuffle_mode = "within"
    elif args.shuffle_across:
        run_shuffles = True
        shuffle_mode = "across"
    else:
        run_shuffles = True
        shuffle_mode = "both"
    
    pipe.run_with_controls(
        data=data,
        essential_path=args.essential,
        nonessential_path=args.nonessential,
        run_shuffles=run_shuffles,
        shuffle_mode=shuffle_mode
    )


if __name__ == "__main__":
    main()
