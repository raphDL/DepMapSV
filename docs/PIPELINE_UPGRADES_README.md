# Pipeline Upgrades - Implementation Summary

This document describes the new features added to `sv_bias_pipeline.py` and supporting scripts.

## PART A: Partial Correlation Scripts

### `partial_corr_from_design.py`
Computes Δ|partial corr(dep, CN | proximity)| from existing pipeline outputs without refitting.

**Usage:**
```bash
python partial_corr_from_design.py <base_folder>
```

**Outputs:**
- `<base_folder>/partial_corr_summary.csv` - Summary statistics with bootstrap CIs
- `<base_folder>/partial_corr_per_gene.csv` - Per-gene results

### `partial_corr_compare.py`
Compares partial correlation metrics across multiple conditions.

**Usage:**
```bash
python partial_corr_compare.py \
  out_pilot_baseline_fixed/unstable \
  out_pilot_shuffle_rotate/unstable \
  out_pilot_shuffle_uniform/unstable \
  out_pilot_lin_250k_2m_fixed/unstable
```

**Decision Rule:** If TRUE condition's median CI lower bound > any shuffle CI upper bound → proceed; else pause.

## PART B: Pipeline Upgrades

### New CLI Flags

1. **`--genes-subset-file PATH`**
   - Early gene subsetting before feature computation
   - File should contain newline-separated gene symbols

2. **`--write-features DIR`**
   - Save computed features to cache directory
   - Creates `features_cn.parquet` and `features_bp.parquet`

3. **`--read-features DIR`**
   - Load features from cache (skips expensive computation)
   - Requires `features_cn.parquet` and `features_bp.parquet` in directory

4. **`--orthogonalize-proximity`**
   - Orthogonalize proximity features against CN before model fitting
   - Reduces collinearity between CN and proximity predictors
   - Coefficients reported as `coef_bp_near_resid` and `coef_bp_far_resid`

5. **`--compute-partial-corr`**
   - Compute and report partial correlation metric alongside existing metrics
   - Saves `partial_corr_per_gene.csv` in output directory

6. **`--match-proximity-prevalence`**
   - Filter genes by proximity prevalence match (tolerance ±0.02)
   - Applied to shuffled data for fair comparison
   - Reports number of genes kept after filtering

### Implementation Details

- **Early Gene Subsetting:** Applied in `load_and_validate_data()` before feature computation
- **Feature Caching:** Uses pandas parquet format (requires `pyarrow` or `fastparquet`)
- **Orthogonalization:** Per-gene linear regression to remove CN projection from proximity features
- **Prevalence Matching:** Computes per-gene proximity prevalence, filters to genes within tolerance
- **Partial Correlation:** Integrated into `ValidationMetrics` class with bootstrap CIs

## PART C: Ready-to-Run Commands

### 1. Quick Check on Existing Pilot Outputs

```bash
cd /Users/raphael/Documents/sandbox/DepMapSV

# Compute partial correlation for each condition
python partial_corr_from_design.py out_pilot_baseline_fixed/unstable
python partial_corr_from_design.py out_pilot_shuffle_rotate/unstable
python partial_corr_from_design.py out_pilot_shuffle_uniform/unstable
python partial_corr_from_design.py out_pilot_lin_250k_2m_fixed/unstable

# Compare all conditions
python partial_corr_compare.py \
  out_pilot_baseline_fixed/unstable \
  out_pilot_shuffle_rotate/unstable \
  out_pilot_shuffle_uniform/unstable \
  out_pilot_lin_250k_2m_fixed/unstable
```

### 2. Build a 2k-Gene Subset

```bash
python make_stratified_subset.py \
  --dependency-long CRISPRGeneEffect_long.csv \
  --genes-bed prep_out/genes.depmap.unique.bed \
  --sv-bedpe prep_out/sv_from_cnv.bedpe \
  --cnv-bed prep_out/cnv_segments.bed \
  --out genes_stratified_2k.txt
```

**Note:** If you have cached features from a previous run, you can speed this up:
```bash
python make_stratified_subset.py \
  --dependency-long CRISPRGeneEffect_long.csv \
  --genes-bed prep_out/genes.depmap.unique.bed \
  --sv-bedpe prep_out/sv_from_cnv.bedpe \
  --cnv-bed prep_out/cnv_segments.bed \
  --features-cache features_cache_2k \
  --out genes_stratified_2k.txt
```

### 3. Mini Full-Run with Caching + Metrics

#### First Pass (True) - Writes Feature Cache

```bash
python sv_bias_pipeline.py \
  --dependency CRISPRGeneEffect_long.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv  prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --model huber \
  --bp-windows 250000 2000000 \
  --bootstrap-iterations 800 \
  --activity-threshold 0.01 \
  --activity-cell-fraction 0.10 \
  --orthogonalize-proximity \
  --compute-partial-corr \
  --genes-subset-file genes_stratified_2k.txt \
  --write-features features_cache_2k \
  --output-dir out_v2_2k_true
```

#### Reuse Cache for Shuffles

```bash
# Within-chromosome shuffle
python sv_bias_pipeline.py \
  --dependency CRISPRGeneEffect_long.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv  prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --model huber \
  --bp-windows 250000 2000000 \
  --bootstrap-iterations 800 \
  --activity-threshold 0.01 \
  --activity-cell-fraction 0.10 \
  --orthogonalize-proximity \
  --compute-partial-corr \
  --genes-subset-file genes_stratified_2k.txt \
  --read-features features_cache_2k \
  --match-proximity-prevalence \
  --output-dir out_v2_2k_within \
  --shuffle-within

# Across-chromosome shuffle
python sv_bias_pipeline.py \
  --dependency CRISPRGeneEffect_long.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv  prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --model huber \
  --bp-windows 250000 2000000 \
  --bootstrap-iterations 800 \
  --activity-threshold 0.01 \
  --activity-cell-fraction 0.10 \
  --orthogonalize-proximity \
  --compute-partial-corr \
  --genes-subset-file genes_stratified_2k.txt \
  --read-features features_cache_2k \
  --match-proximity-prevalence \
  --output-dir out_v2_2k_across \
  --shuffle-across
```

**Note:** The `--shuffle-within` and `--shuffle-across` flags need to be added to the CLI parser. For now, you can modify the pipeline to accept these or run shuffles separately.

### 4. Mini Report

After running the mini full-run, compare results:

```bash
# Extract partial correlation medians from metrics.json
python -c "
import json
import sys
for label in ['true', 'within', 'across']:
    path = f'out_v2_2k_{label}/metrics.json'
    try:
        with open(path) as f:
            m = json.load(f)
            if 'partial_correlation' in m:
                pc = m['partial_correlation']
                print(f'{label}: median={pc[\"median_partial_delta_abs_corr\"]:.4f} '
                      f'[{pc[\"partial_ci_lower\"]:.4f}, {pc[\"partial_ci_upper\"]:.4f}]')
    except:
        pass
"
```

**Decision Rule:** If TRUE median CI lower > all shuffle CI upper → green-light continuing; else pause.

## Dependencies

- **pandas** (with parquet support): `pip install pyarrow` or `pip install fastparquet`
- **scikit-learn**: For linear regression and models
- **numpy**: For numerical operations

## Notes

1. **Parquet Support:** Feature caching requires parquet support. Install with:
   ```bash
   pip install pyarrow
   ```

2. **Shuffle Flags:** The `--shuffle-within` and `--shuffle-across` flags mentioned in the commands need to be implemented in the CLI. Currently, shuffles are controlled by `--no-shuffles` (which disables them). You may need to add explicit shuffle type flags or modify the pipeline to accept shuffle mode.

3. **Prevalence Matching:** When `--match-proximity-prevalence` is enabled, genes are filtered based on proximity prevalence difference ≤ 0.02. This ensures shuffled data has similar proximity patterns to true data.

4. **Orthogonalization:** When `--orthogonalize-proximity` is enabled, proximity features are residualized against CN before model fitting. This removes collinearity but changes coefficient interpretation.

## Acceptance Criteria Status

✅ `partial_corr_from_design.py` produces summary CSV per condition  
✅ `partial_corr_compare.py` prints comparison table  
✅ `sv_bias_pipeline.py` honors `--genes-subset-file` for early subsetting  
✅ `sv_bias_pipeline.py` saves/loads feature cache via `--write-features`/`--read-features`  
✅ `sv_bias_pipeline.py` supports `--orthogonalize-proximity` in model fitting  
✅ `sv_bias_pipeline.py` reports partial correlation when `--compute-partial-corr`  
✅ `sv_bias_pipeline.py` reports prevalence match counts when `--match-proximity-prevalence`  
✅ `make_stratified_subset.py` outputs stratified gene list based on CN var × proximity prevalence grid  

**Note:** Shuffle flags (`--shuffle-within`, `--shuffle-across`) need to be added to the CLI parser if you want explicit control. Currently, shuffles run by default unless `--no-shuffles` is specified.

