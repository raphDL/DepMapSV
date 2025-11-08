# WGS SV Data Ingestion

This document describes how to convert CCLE SvABA translocation data from Excel format to BEDPE format for use with the SV bias pipeline.

## Overview

The CCLE SvABA translocations dataset (`CCLE_translocations_SvABA_20181221.xlsx`) contains structural variant calls from whole-genome sequencing. This script converts it to BEDPE format (0-based, half-open coordinates) suitable for the SV bias analysis pipeline.

## Prerequisites

```bash
pip install pandas openpyxl
```

## Convert Excel to BEDPE

### Basic Usage

```bash
cd /Users/raphael/Documents/sandbox/DepMapSV

python convert_svaba_xlsx_to_bedpe.py \
  --infile CCLE_translocations_SvABA_20181221.xlsx \
  --out-combined sv_wgs_all.bedpe \
  --out-dir sv_wgs_bedpe
```

This will:
- Read the Excel file
- Filter for TRA-like and INV-like events (default)
- Parse breakpoint coordinates
- Convert to BEDPE format (0-based, half-open, point breakends)
- Output combined BEDPE: `sv_wgs_all.bedpe`
- Output per-sample BEDPEs: `sv_wgs_bedpe/{sample_id}.bedpe`

### With Sample Name Mapping

If you have a mapping file to convert CCLE names to DepMap IDs:

```bash
python convert_svaba_xlsx_to_bedpe.py \
  --infile CCLE_translocations_SvABA_20181221.xlsx \
  --map-csv templates/depmap_model_map.csv \
  --out-combined sv_wgs_all.bedpe \
  --out-dir sv_wgs_bedpe
```

The mapping CSV should have columns: `CCLE_name,DepMap_ID`

To create the mapping file:
1. Get DepMap sample metadata (e.g., from `sample_info.csv`)
2. Create a CSV with columns `CCLE_name,DepMap_ID`
3. Map CCLE names (as they appear in the SvABA Excel) to DepMap model IDs

See `templates/depmap_model_map.csv` for a template.

### Custom Class Filtering

To keep different SV classes:

```bash
python convert_svaba_xlsx_to_bedpe.py \
  --infile CCLE_translocations_SvABA_20181221.xlsx \
  --keep-classes TRA-like INV-like DEL-like \
  --out-combined sv_wgs_all.bedpe \
  --out-dir sv_wgs_bedpe
```

## Spot-Check Outputs

After conversion, verify the outputs:

```bash
# Check combined BEDPE
head -5 sv_wgs_all.bedpe
wc -l sv_wgs_all.bedpe

# Check per-sample files
ls -1 sv_wgs_bedpe | head
wc -l sv_wgs_bedpe/*.bedpe | tail -5
```

Expected BEDPE format (tab-delimited):
```
chrom1  start1  end1  chrom2  start2  end2  name  score  strand1  strand2  sample  info
```

Where:
- `start1/start2`: pos-1 (0-based)
- `end1/end2`: pos (half-open, so end = start+1 for point breakends)
- `info`: semicolon-separated key=value pairs with metadata

## Run Pipeline with Real WGS SVs

Once you have the BEDPE files, run the SV bias pipeline:

```bash
python sv_bias_pipeline.py \
  --dependency CRISPRGeneEffect_long.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv sv_wgs_all.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --model huber \
  --bp-windows 250000 2000000 \
  --bootstrap-iterations 2000 \
  --activity-threshold 0.01 \
  --activity-cell-fraction 0.10 \
  --output-dir out_v3_wgs
```

Or using the V3 pipeline modules:

```bash
# Step 1: Ingest WGS SV data
python sv_ingest_wgs.py \
  --sv-dir sv_wgs_bedpe/ \
  --genes prep_out/genes.depmap.unique.bed \
  --cnv prep_out/cnv_segments.bed \
  --out out_v3/true/design_matrix_wgs.csv \
  --svtypes TRA INV

# Step 2: Fit models
python run_models.py \
  --dep CRISPRGeneEffect_long.csv \
  --design out_v3/true/design_matrix_wgs.csv \
  --out out_v3/true \
  --kernel prox_exp_100k \
  --min-cells 200
```

## BEDPE Format Details

The converter produces BEDPE files with:
- **0-based, half-open coordinates**: `[start, end)`
- **Point breakends**: For a breakpoint at position `pos`, we use `start=pos-1, end=pos`
- **Chromosome format**: All chromosomes prefixed with `chr` (e.g., `chr1`, `chrX`, `chrM`)
- **Strand information**: Preserved from original breakpoint strings (`+` or `-`)

## Validation

The script includes automatic sanity checks:
- All chromosomes have `chr` prefix
- All coordinates satisfy `0 ≤ start < end` and `end = start+1` (point breakends)
- Logs counts and top samples

Run unit tests:

```bash
pytest tests/test_bp_parse.py -v
```

## Troubleshooting

### Unparsable Breakpoints

If you see warnings about unparsable breakpoints, check the format in the Excel file. The parser expects:
- `chrom:pos-pos(strand)` or `chrchrom:pos-pos(strand)`
- Examples: `13:33777850-33777850+`, `chrX:123456-123456-`

### Missing Sample Mappings

If using `--map-csv` but some samples aren't mapped, the script will fall back to using the CCLE_name as the sample_id.

### Coordinate Issues

If sanity checks fail, verify:
- Breakpoint positions are non-negative integers
- The Excel file format matches expected columns

