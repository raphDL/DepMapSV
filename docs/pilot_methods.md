# Pilot Methods: Fast Subset Reanalysis for Structural-Variant Proximity Bias

## Overview

We implemented a fast pilot mode to validate the structural-bias reanalysis pipeline on a targeted subset of genes before running the full ~18,000-gene analysis. The pilot selects genes based on per-gene breakpoint prevalence and fits models with proximity features enabled.

**Methods: directional metric & strong shuffles.** We use a directional proximity-active metric (genes with ≥10% active cells AND prox-only |Δr| > 0) as the primary falsification metric, validated against strong negative controls (rotate and uniform shuffles of breakpoint positions within chromosome×cell) to ensure signal is breakpoint-driven and robust.

## Gene Selection

**Selection universe:** We restrict to autosomal genes (excluding chrX/chrY) and cell lines present in both dependency data and SV/CNV data.

**Prevalence calculation:** For each gene, we compute breakpoint distances for all (gene, cell) pairs in the dependency × SV intersection. We then calculate per-gene prevalence metrics:
- `p100k`: fraction of cells where the gene is within 100 kb of a breakpoint
- `p1m`: fraction of cells where the gene is within 1 Mb of a breakpoint

**Unstable set:** Genes with `p100k ≥ 0.05` OR `p1m ≥ 0.20` (default: 50 genes, sampled randomly if more candidates exist).

**Stable set:** Genes with `p100k = 0.0` AND `p1m ≤ 0.005` (default: 50 genes, sampled randomly from non-overlapping candidates).

**Disjointness:** Stable genes are explicitly excluded from the unstable set to ensure clear contrast.

## Feature Engineering

**Copy number:** Length-weighted gene-level CN from CNV segments, capped at 8.0.

**Breakpoint proximity:**
- Binary windows: `bp_within_100000` (≤100 kb), `bp_within_1000000` (≤1 Mb)
- Continuous term: `inv_bp = 1 / (bp_dist + 1e3)` (when `--add-continuous-proximity` is enabled)

**Expression:** Optional (omitted in pilot if ID mismatches cause low non-NA rates).

## Model Fitting

**Regression:** Per-gene robust regression (HuberRegressor by default, `alpha=1e-4`, `epsilon=1.35`).

**Standardization:** Optional per-gene z-scoring of predictors (`--standardize-predictors`). When enabled, coefficients are converted back to raw scale for contribution computation.

**Training data:** Inner join of dependency × features (only rows present in both), ensuring training prevalence matches the selection universe.

## Robust Metrics

**Directional proximity-active (primary falsification metric):** Fraction of genes that satisfy BOTH:
1) ≥10% of cells have `|proximity_contribution| ≥ threshold` and 2) prox-only `|Δr| = |corr(CN, dep)| - |corr(CN, dep - prox)| > 0`.

Thresholds: `contrib_thresh = 0.01` when `--standardize-predictors` (predictors z-scored), else `0.1`.

**Non-directional proximity-active (supplemental):** Fraction of genes with ≥10% active cells (ignores `|Δr|` sign). Reported for context only.

**Proximity-only correction:** Remove only proximity effects: `dependency_prox_only_corrected = dependency - proximity_contribution`. Compute `|corr(dep, CN)|` before and after; report median `|Δr|` across genes.

**Contribution computation:** Uses raw predictor columns with coefficients converted from standardized scale if needed, ensuring contributions are in dependency units.

## Guardrails

**Gene overlap:** Assert > 0 genes overlap between dependency and gene BED.

**Cell overlap:** Warn if < 500 cells overlap across dependency/CNV/SV.

**Selection sanity:** Assert `unstable_median_p1m > stable_median_p1m` after selection.

**Contrast checks:** After both sets run, verify unstable > stable on:
- Selection median p1m
- Training mean p1m
- flag_bp_close rate
- Directional prox-active genes fraction
- Prox-only median |Δr|

**Expression:** Warn if expression non-NA rate < 5% (likely ID mismatch).

## Outputs

**Per set (`out/<set>/`):**
- `design_matrix.csv`: Training rows with all features (raw + standardized if used)
- `dependency_corrected.csv`: Full correction (CN + proximity)
- `models_coefficients.csv`: Per-gene model coefficients
- `qc_flags.csv`: QC flags (flag_bp_close, flag_cn, confounded)
- `pilot_summary.txt`: Key metrics (selection medians, training means, robust metrics)
- `manifest.json`: Args, thresholds, n_genes, n_cells

**Evaluation:** `evaluation_summary.csv` omitted if pilot genes don't overlap essential/nonessential sets (noted in summary).

## Data Sources

- **Dependency:** DepMap 25Q3 CRISPR gene-effect scores (`CRISPRGeneEffect.csv`)
- **SV:** Structural variants inferred from CNV segments (`prep_out/sv_from_cnv.bedpe`)
- **CNV:** Copy-number segments (`prep_out/cnv_segments.bed`)
- **Genes:** DepMap-matched gene annotations (`prep_out/genes.depmap.unique.bed`)

## Runtime
## Negative Controls (Falsification)

We implement two strong shuffles of SV breakpoints applied per cell × chromosome:

- **Rotate:** Add a random offset modulo chromosome length to all SV positions (preserves relative structure, destroys alignment to genes).
- **Uniform:** Resample midpoints uniformly within chromosome lengths while preserving segment lengths.

We expect the directional prox-active metric to drop toward low values under these shuffles, while remaining high under robust settings (e.g., linear + wider windows).

Pilot completes in ~100 seconds (both sets, 50 genes each, 830 cells) vs hours for full genome, enabling rapid iteration on thresholds and model parameters.

