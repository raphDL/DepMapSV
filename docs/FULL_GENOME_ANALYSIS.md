# Full-Genome Analysis Setup

## Overview

This setup processes all ~17,697 autosomal genes in chunks of 1,000 to validate the proximity effect at full scale using a pairwise rotate shuffle control.

## Files Created

1. **`all_genes_autosomes.txt`**: List of 17,697 autosomal genes from `prep_out/genes.depmap.unique.bed`
2. **`chunks/fullgenes_*`**: 18 chunk files (1,000 genes each, last chunk has 697)
3. **`run_full_genome.sh`**: Script to run all chunks (non-rotated and rotate shuffle)
4. **`summarize_full_genome.py`**: Script to compute directional flags and summary statistics

## Settings

- **Model**: Linear regression
- **Windows**: 250k / 1.5Mb (recommended) or 250k / 2Mb
- **Predictors**: Standardized (contrib_thresh = 0.01)
- **Active fraction**: 0.15 (can be changed to 0.20 in `summarize_full_genome.py`)

## Usage

### Step 1: Run full-genome analysis

```bash
export MPLCONFIGDIR="$PWD/.mplconfig" XDG_CACHE_HOME="$PWD/.cache"
export PYTHONHASHSEED=1 OMP_NUM_THREADS=1
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

bash run_full_genome.sh
```

This will:
- Process 18 chunks × 2 conditions = 36 runs
- Output to `out_full_linear/fullgenes_*/unstable/` (non-rotated)
- Output to `out_full_linear_rotate/fullgenes_*/unstable/` (rotate shuffle)

**Expected runtime**: Several hours (depends on system)

### Step 2: Summarize results

```bash
python summarize_full_genome.py
```

This generates:
- `figs_full/full_directional_flags.csv`: Per-gene directional flags (non-rotated, rotate, pairwise_excess)
- `figs_full/full_directional_summary.csv`: Summary statistics with bootstrap CIs

## Success Criteria

1. **directional (all genes) >> rotate**: CIs should be clearly separated
2. **pairwise_excess mean > 0**: 95% CI should not cross 0
3. **Optional**: Stratify by breakpoint prevalence deciles to see monotonic increase

## Adjusting Active Fraction

To use a stricter threshold (0.20 instead of 0.15), edit `summarize_full_genome.py`:

```python
ACTIVE_FRAC = 0.20  # Change from 0.15
```

No need to re-fit models; just re-run the summary script.

## Output Structure

```
out_full_linear/
  fullgenes_000/
    unstable/
      design_matrix.csv
      models_coefficients.csv
      pilot_summary.txt
      ...
  fullgenes_001/
    unstable/
      ...
  ...

out_full_linear_rotate/
  fullgenes_000/
    unstable/
      ...
  ...
```

## Notes

- Uses `--pilot --pilot-genes-file` to reuse the pilot pipeline without code changes
- Keeps memory bounded by processing in chunks
- All conditions use the same model settings for fair comparison
- Pairwise comparison (gene-by-gene) controls for gene-specific effects


