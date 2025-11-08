# Quick Start: Using WGS SV Data

## Step-by-Step Workflow

### 1. Convert Excel to BEDPE

First, run the converter to create BEDPE files from the Excel file:

```bash
cd /Users/raphael/Documents/sandbox/DepMapSV

python convert_svaba_xlsx_to_bedpe.py \
  --infile CCLE_translocations_SvABA_20181221.xlsx \
  --out-combined sv_wgs_all.bedpe \
  --out-dir sv_wgs_bedpe
```

This will create:
- `sv_wgs_all.bedpe` - Combined BEDPE file
- `sv_wgs_bedpe/` - Directory with per-sample BEDPE files

### 2. Verify Outputs

Check that files were created:

```bash
# Check combined file
head -5 sv_wgs_all.bedpe
wc -l sv_wgs_all.bedpe

# Check per-sample directory
ls -1 sv_wgs_bedpe | head -10
ls -1 sv_wgs_bedpe | wc -l
```

### 3. Run V3 Pipeline

Once BEDPE files exist, run the pipeline:

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

## Troubleshooting

### Error: "No BEDPE files found in sv_wgs_bedpe/"

**Solution:** Run the converter first (Step 1 above). The directory and files are created by the converter.

### Error: "SV directory does not exist"

**Solution:** Make sure you've run the converter and the `sv_wgs_bedpe/` directory exists.

### Error: File format issues

The `sv_ingest_wgs.py` script now supports:
- Standard BEDPE format (no header, 11-12 columns)
- Header format (with column names)

The converter outputs standard BEDPE format, which is automatically detected.

## Alternative: Use Combined BEDPE File

You can also use the combined BEDPE file directly:

```bash
python sv_ingest_wgs.py \
  --sv-dir sv_wgs_all.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --cnv prep_out/cnv_segments.bed \
  --out out_v3/true/design_matrix_wgs.csv \
  --svtypes TRA INV
```

Note: When using a single file, specify the file path directly (not a directory).

