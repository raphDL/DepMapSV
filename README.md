# DepMapSV: Structural Variant Proximity Analysis Pipeline

Pipeline for analyzing structural variant (SV) proximity effects on gene dependency using DepMap CRISPR data and WGS SV breakpoints.

## Quick Start

### 1. Install dependencies

```bash
conda env create -f environment.yml
conda activate svbias
```

### 2. Prepare test data (small example)

```bash
# Generate synthetic test data
python test_synthetic.py

# Run ingestion on test data
python sv_ingest_wgs.py \
  --sv-dir data/test \
  --genes data/test/genes.bed \
  --cnv data/test/cnv.csv \
  --out out_test/design_matrix.parquet
```

### 3. Run models

```bash
python run_models.py \
  --design out_test/design_matrix.parquet \
  --dep data/test/dependency.csv \
  --kernel prox_exp_100k \
  --min-cells 10 \
  --model huber \
  --out out_test/models
```

## Full Pipeline

For production runs with real data, see:
- **`docs/QUICK_START_WGS.md`** - Complete workflow with WGS SV data
- **`docs/V3_OVERNIGHT_RUN.md`** - One-command overnight run
- **`docs/WGS_SV_INGEST.md`** - Detailed ingestion guide

## Data Requirements

Large/third-party datasets (DepMap, CCLE, PCAWG) should be placed in `data/external/` (gitignored). See `data/external/README.md` for details.

**Required inputs:**
- WGS SV BEDPE files (one per sample)
- Gene annotations (BED format)
- CNV segments (BED format)
- CRISPR dependency data (DepMap format)

**Optional:**
- `sample_info.csv` for ACH↔CCLE name mapping
- RNAi data for orthogonal validation

## Key Features

- **Robust data ingestion**: Flexible sample name mapping, CNV autodetection, chromosome normalization
- **Vectorized distance computation**: Fast, scalable breakpoint-to-gene distance calculations
- **Multiple proximity kernels**: Exponential RBFs, inverse distance, boolean windows
- **Negative controls**: ROTATE and WITHIN shuffles for statistical validation
- **Comprehensive diagnostics**: QC tables, overlap checks, correlation monitoring

## Project Structure

```
DepMapSV/
├── README.md              # This file
├── LICENSE                # MIT License
├── CITATION.cff           # Citation metadata
├── environment.yml        # Conda environment
├── Makefile               # Build automation
├── sv_ingest_wgs.py       # Main ingestion script
├── run_models.py          # Model fitting
├── scripts/               # Helper scripts
├── configs/               # Configuration files
├── tests/                 # Unit tests
├── data/
│   ├── test/              # Small test data (tracked)
│   └── external/          # Large/licensed data (gitignored)
└── docs/                  # Documentation
    ├── figs/              # Figures
    └── *.md               # Project docs
```

## Documentation

All documentation is in `docs/`:
- **QUICK_START_WGS.md** - Full workflow guide
- **WGS_SV_INGEST.md** - Ingestion details
- **TESTING.md** - Testing guide
- **PROJECT_STATUS.md** - Current status and findings

## Common Issues

### No overlap between SV and CNV

- Ensure sample names match (use `--sample-info` for ACH↔CCLE mapping)
- Check `logs/overlap_debug_samples.csv` for examples
- See `docs/WGS_SV_INGEST.md` for troubleshooting

### Parquet engine missing

- Install `pyarrow`: `mamba install -c conda-forge pyarrow`
- Pipeline automatically falls back to CSV if Parquet unavailable

## Citation

If you use this software, please cite it. See `CITATION.cff` for details.

## License

MIT License - see `LICENSE` for details.

## Contributing

This is a research pipeline. For questions or issues, please open a GitHub issue.
