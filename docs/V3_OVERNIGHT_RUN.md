# V3 Overnight Run - One-Command Pipeline

## Overview

The V3 overnight run pipeline automates the complete analysis workflow:
1. Ingests WGS SV data
2. Fits TRUE + ROTATE models across multiple kernels
3. Evaluates TRUE vs ROTATE
4. Selects case-study genes with diagnostics

## Quick Start

```bash
cd /Users/raphael/Documents/sandbox/DepMapSV
./scripts/run_full_v3.sh
```

## Configuration

Edit `configs/v3.yaml` to customize:
- Input paths (dependency, genes, CNV, SV directory)
- Ingest settings (join type, drop flags, spline kernel)
- Model settings (kernels, min_cells, model type)
- Shuffle settings (rotate, within-chrom)
- Evaluation and selection thresholds

## Pipeline Steps

1. **Optional annotation**: Annotate SVs with PCAWG if index available
2. **Ingest TRUE design**: Build design matrix from WGS SVs (Parquet output)
3. **Generate shuffles**: Create ROTATE (and optionally WITHIN) controls
4. **Fit TRUE models**: Run models across all kernels (prox_exp_50k, prox_exp_100k, prox_exp_250k, prox_exp_500k, prox_inv)
5. **Fit ROTATE models**: Run models on shuffled data across all kernels
6. **Evaluate and select**: Compare TRUE vs ROTATE, apply pre-registered criteria, select case studies

## Output Structure

```
out_v3/
├── true/
│   ├── design_matrix_wgs.parquet
│   ├── prox_exp_50k/
│   │   ├── model_coefficients.csv
│   │   ├── dependency_corrected.csv
│   │   └── metrics.json
│   ├── prox_exp_100k/
│   └── ...
├── rotate/
│   ├── design_matrix_wgs_rotate.parquet
│   ├── prox_exp_50k/
│   └── ...
└── summary/
    ├── evaluation_metrics.json
    └── case_study_genes.csv
```

## Logs

All steps log to `logs/`:
- `01_annotate.log` - Optional PCAWG annotation
- `02_ingest_true.log` - TRUE design matrix ingestion
- `03_shuffles.log` - Shuffle generation
- `04_fit_true.log` - TRUE model fitting
- `05_fit_rotate.log` - ROTATE model fitting
- `06_evaluate_select.log` - Evaluation and case study selection

## Defaults

- **Join**: `inner` (prevents inflating design with CN-only lines)
- **Drop no-SV lines**: `true` (keeps only rows with proximity info)
- **Kernels**: `prox_exp_50k, prox_exp_100k, prox_exp_250k, prox_exp_500k, prox_inv`
- **Min cells**: `50` (lower threshold for faster iteration)
- **Model**: `huber` (robust regression)
- **FDR threshold**: `0.10`

## Requirements

- Python packages: `pandas`, `numpy`, `scipy`, `scikit-learn`, `statsmodels`, `pyyaml`
- For Parquet: `pyarrow` or `fastparquet`
- For optional annotation: `annotate_sv_catalogs.py` (if using PCAWG index)

## Notes

- Design matrices are saved as Parquet for efficiency
- All kernels are tested in parallel (sequential execution)
- Case study selection uses pre-registered criteria
- Missing columns in coefficient files are handled gracefully (criteria fail safely)

