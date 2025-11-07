# Implementation Summary: Pilot Mode Validation & Hit Scoring

## What We Built

### 1. Enhanced Pilot Outputs

**Problem:** Pilot mode needed more detailed outputs to validate that hits are real.

**Solution:**
- **`design_matrix.csv`** per set: Contains all training rows (gene×cell pairs) with:
  - Raw features: `dependency`, `cn`, `bp_dist`, `bp_within_100000`, `bp_within_1000000`, `inv_bp`
  - Standardized features (if `--standardize-predictors` used): `cn_std`, `bp_within_100000_std`, etc.
  - Corrected dependency: `dependency_corrected`
- **`pilot_summary.txt`** enhanced with robust metrics:
  - `prox_active_genes_frac`: Fraction of genes where ≥10% cells have meaningful proximity contributions
  - `prox_only_median_abs_corr_drop`: Median reduction in CN–dependency coupling when removing only proximity
  - `cn_var_median`: Median CN variance (context)
  - `n_cells_training_median/iqr`: Training sample sizes
  - Top 5 genes by |Δcorr|
- **`manifest.json`**: Records all arguments, thresholds, and metadata for reproducibility

### 2. Hit Scoring System

**Problem:** Need to identify which genes are truly "interesting" (real proximity artifacts vs noise).

**Solution:**
- **`score_hits.py`**: Scores each gene on 5 criteria:
  1. Activity: `prox_active_frac >= 0.20` (proximity matters in ≥20% of cells)
  2. Strength: `prox_contrib_median_if_active >= 0.01` (when active, contribution is meaningful)
  3. Direction: `prox_only_abs_corr_drop > 0` (proximity ADDS CN coupling = artifact)
  4. CN≈2 robustness: Effect persists in CN-neutral cells (1.8 ≤ CN ≤ 2.2)
  5. Lineage breadth: Effect seen in ≥2 lineages (if lineage data available)
- **Labels:**
  - **GREEN**: Score ≥4 (or ≥3 if lineage missing) = strong artifact candidate
  - **AMBER**: Score 2-3 = moderate signal
  - **RED**: Score 0-1 = weak/no signal

### 3. Case Study Visualization

**Problem:** Need to visualize why specific genes are interesting.

**Solution:**
- **`gene_panel.py`**: Creates 3-panel figure per gene:
  - Panel 1: Dependency vs distance-to-breakpoint (before correction)
  - Panel 2: CN vs dependency (original vs prox-only corrected)
  - Panel 3: Histogram of proximity contributions
- **`make_case_studies.py`**: Batch generates panels for top hits
- **`make_pilot_figs.py`**: Summary figures (prox-active fraction bar chart, prox-only |Δr| boxplot)
- All figures saved as **PDF** (not PNG)

### 4. Negative Control (Shuffled Breakpoints)

**Problem:** Need to verify signal isn't spurious.

**Solution:**
- Added `--shuffle-sv` flag to `svbias_reanalysis.py`
- Shuffles breakpoint positions within each chromosome×cell combination
- **Expected:** If signal is real, shuffled run should show:
  - prox_active_frac in "unstable" collapses toward "stable" levels
  - Fewer/no GREEN hits in unstable set
- **Ran:** `out_pilot_shuffle/` with shuffled SV

### 5. Robustness Check (Different Windows/Model)

**Problem:** Need to verify signal is robust to parameter choices.

**Solution:**
- Ran pilot with different parameters:
  - Windows: `--bp-windows 250000 2000000` (instead of 100k/1Mb)
  - Model: `--model linear` (instead of huber)
- **Expected:** 
  - Unstable > stable on prox-active fraction
  - Some GREEN hits persist (rank stability)
- **Ran:** `out_pilot_lin_250k_2m/` with linear model and wider windows

### 6. Optional Lineage Analysis

**Problem:** Want to check if effects generalize across cell lineages.

**Solution:**
- **`add_lineage_to_design.py`**: Merges cell_line→lineage mapping into design_matrix.csv
- **`score_hits.py`** computes lineage breadth (how many lineages show the effect)
- Can be run later when lineage data is available

## Results Summary

### Original Run (`out_pilot_final/`)
- **Unstable:** 4 GREEN, 11 AMBER, 29 RED
- **Stable:** 0 GREEN, 16 AMBER, 28 RED
- **Top GREEN hits:** DCANP1, ZNF524, KIAA1671, LRRC51
- **Contrast:** Clear unstable > stable on all metrics ✅

### Shuffled Run (`out_pilot_shuffle/`)
- **Status:** Completed, needs comparison to original
- **Expected:** Should show collapse of signal (fewer GREEN hits)

### Robustness Run (`out_pilot_lin_250k_2m/`)
- **Unstable:** 4 GREEN hits (different genes: INTS3, NFS1, PACS1, RPS6KB2)
- **Stable:** 3 GREEN hits (ZFX, AMMECR1, STARD8)
- **Note:** Some overlap with original (e.g., PACS1 was in original top-5)

## Key Files Created

### Scripts
- `score_hits.py` - Score and label hits GREEN/AMBER/RED
- `hits_summary.py` - Rank genes by proximity signal
- `gene_panel.py` - Generate 3-panel case study per gene
- `make_case_studies.py` - Batch generate case study panels
- `make_pilot_figs.py` - Generate summary figures
- `add_lineage_to_design.py` - Merge lineage data (optional)
- `pilot_bootstrap_ci.py` - Bootstrap CIs and permutation tests

### Documentation
- `pilot_methods.md` - 1-page methods section
- `pilot_figure_captions.md` - Figure captions
- `HITS_VALIDATION_SUMMARY.md` - Quick reference for hit validation
- `DATA_VALIDATION_REPORT.md` - Detailed analysis of data realism
- `PILOT_SUMMARY.md` - Overall pilot summary

## What's Next (Per ChatGPT Suggestions)

1. **Compare shuffled vs original:** Check if prox_active_frac collapsed in shuffled run
2. **Compare robustness runs:** Verify rank stability of GREEN hits across parameter settings
3. **Add lineage (when available):** Merge lineage data and rescore with stricter GREEN threshold
4. **CN≈2 check:** Verify effects persist in CN-neutral cells (already computed in score_hits.py)

## Quick Commands

```bash
# Score hits (original run)
python score_hits.py out_pilot_final/unstable --ignore-lineage
python score_hits.py out_pilot_final/stable --ignore-lineage

# Generate figures
python make_pilot_figs.py
python make_case_studies.py

# Negative control (already run)
# python svbias_reanalysis.py ... --shuffle-sv --out out_pilot_shuffle

# Robustness (already run)
# python svbias_reanalysis.py ... --bp-windows 250000 2000000 --model linear --out out_pilot_lin_250k_2m
```

