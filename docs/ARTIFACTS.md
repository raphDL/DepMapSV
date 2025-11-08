# Project Artifacts and Results

This document catalogs the key outputs and results from the SV/CN proximity correction project.

## Final Results (out_v2/)

The final pipeline run using `sv_bias_pipeline.py` (v2.1) with bootstrap CIs and negative controls:

- `out_v2/comparison_report.txt` - Summary comparing true vs shuffled controls
- `out_v2/true/` - Results from true SV breakpoint positions
  - `dependency_corrected.csv` - Corrected dependency scores
  - `model_coefficients.csv` - Per-gene model coefficients
  - `gene_level_correlations.csv` - Gene-level CN correlation changes
  - `metrics.json` - Summary metrics with bootstrap CIs
- `out_v2/shuffle_within_chrom/` - Within-chromosome shuffle control
- `out_v2/shuffle_across_chrom/` - Cross-chromosome shuffle control

### Key Findings from out_v2/

- **Directional prox-active fraction**: True = 31.3%, Shuffled = 52.9-57.2% (negative excess)
- **CN correlation improvement**: True Δ|corr| = 0.0125, Shuffled = 0.0182-0.0203
- **Essential/non-essential separation**: ΔAUROC ≈ 0 (no improvement)
- **ΔΔAUROC 95% CI**: [0.000000, 0.000005] (overlaps zero)

## Figures

### Pilot Figures (figs_pilot/)
- `final_directional_bar.png/pdf` - Directional metric comparison
- `final_directional_excess_bar.png/pdf` - Excess signal vs controls
- `directional_flag_audit.csv` - Audit of directional flags
- `final_directional_stats.csv` - Summary statistics
- Case study panels: `DCANP1_panel.png`, `KIAA1671_panel.png`, `ZNF524_panel.png`, etc.

### Full Genome Figures (figs_full/)
- `full_directional_linear_only.csv` - Full dataset directional metrics
- `full_directional_pairwise_allgenes.csv` - Pairwise comparisons
- `full_directional_summary.csv` - Summary statistics

## Documentation

- `PILOT_SUMMARY.md` - Pilot analysis summary
- `PILOT_GO_NOGO.md` - Go/no-go decision document
- `FULL_GENOME_ANALYSIS.md` - Full genome analysis results
- `HITS_VALIDATION_SUMMARY.md` - Hit validation summary
- `IMPLEMENTATION_SUMMARY.md` - Implementation details
- `DATA_VALIDATION_REPORT.md` - Data validation results

## Intermediate Outputs

Multiple pilot runs in `out_pilot_*/` directories testing different:
- Window sizes (250k, 1.5M, 2M)
- Model types (linear, robust)
- Shuffling strategies
- Baseline corrections

These can be archived or deleted as they were exploratory.

## Code

All analysis code is in the root directory:
- `sv_bias_pipeline.py` - Final v2.1 pipeline (recommended)
- `svbias_reanalysis.py` - Original pipeline with pilot mode
- Supporting scripts for data prep, visualization, and analysis

