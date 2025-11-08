# External data (not checked in)

Place licensed or large inputs here, e.g. DepMap, CCLE, PCAWG and WGS BEDPE.

This directory is gitignored. See `docs/QUICK_START_WGS.md` for how to obtain them.

## Expected files

- `CRISPRGeneEffect_long.csv` or `.parquet` - DepMap CRISPR dependency data
- `cnv_segments.bed` - Copy number variation segments
- `sv_wgs_bedpe/` - Directory of WGS SV BEDPE files (one per sample)
- `genes.bed` - Gene annotations in BED format
- `sample_info.csv` - DepMap sample metadata (optional, for name mapping)

## Data sources

- **DepMap**: https://depmap.org/portal/download/
- **CCLE**: https://depmap.org/portal/download/
- **PCAWG**: https://dcc.icgc.org/releases/PCAWG

## License

These datasets are subject to their respective licenses. Please review:
- DepMap: https://depmap.org/portal/download/
- CCLE: https://depmap.org/portal/download/
- PCAWG: https://dcc.icgc.org/releases/PCAWG

