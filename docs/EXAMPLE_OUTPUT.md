# Example Outputs from WGS SV Converter

## Console Summary

When running the converter, you'll see output like this:

```
2024-01-15 10:30:45 | INFO | Reading Excel file: CCLE_translocations_SvABA_20181221.xlsx
2024-01-15 10:30:47 | INFO | Total rows read: 15234
2024-01-15 10:30:47 | INFO | Rows kept by class filter (['TRA-like', 'INV-like']): 12456
2024-01-15 10:30:47 | INFO | Parsing breakpoints...
2024-01-15 10:30:48 | WARNING | Unparsable bp1: 23 rows
2024-01-15 10:30:48 | WARNING | Unparsable bp2: 18 rows
2024-01-15 10:30:48 | INFO | Rows with valid breakpoints: 12415
2024-01-15 10:30:48 | INFO | Running sanity checks...
2024-01-15 10:30:48 | INFO | ✓ All chromosomes have 'chr' prefix
2024-01-15 10:30:48 | INFO | ✓ All breakends are point intervals (end=start+1)
2024-01-15 10:30:48 | INFO | ✓ All sanity checks passed
2024-01-15 10:30:48 | INFO | Deduplication: 12415 -> 11892 rows
2024-01-15 10:30:48 | INFO | Wrote combined BEDPE: sv_wgs_all.bedpe (11892 rows)
2024-01-15 10:30:48 | INFO | Writing per-sample BEDPEs to sv_wgs_bedpe/...
2024-01-15 10:30:49 | INFO | Wrote 847 per-sample BEDPE files

============================================================
CONVERSION SUMMARY
============================================================
Total rows read: 15,234
Rows kept by class filter: 12,456
Rows with valid breakpoints: 12,415
Rows after deduplication: 11,892
Unique samples: 847
Total breakpoints written: 11,892

Top 3 samples by event count:
  A549_LUNG: 45 events
  HT29_LARGE_INTESTINE: 38 events
  MCF7_BREAST: 32 events
============================================================
```

## First 5 Lines of sv_wgs_all.bedpe

```bash
head -5 sv_wgs_all.bedpe
```

Expected output (tab-delimited):

```
chr1	33777849	33777850	chr13	45678900	45678901	TRA-like	0	+	-	A549_LUNG	svclass=TRA-like;gene1=GENE1;gene2=GENE2;fusion=GENE1-GENE2
chr2	12345678	12345679	chr5	98765432	98765433	INV-like	0	-	+	HT29_LARGE_INTESTINE	svclass=INV-like;gene1=GENE3;gene2=GENE4
chr3	50000000	50000001	chr8	20000000	20000001	TRA-like	0	+	+	MCF7_BREAST	svclass=TRA-like;fusion=GENE5-GENE6
chr4	15000000	15000001	chr10	30000000	30000001	INV-like	0	-	-	A549_LUNG	svclass=INV-like;gene1=GENE7
chr5	25000000	25000001	chr12	40000000	40000001	TRA-like	0	+	-	HT29_LARGE_INTESTINE	svclass=TRA-like;cosmic_fus=COSMIC123
```

Note:
- Coordinates are 0-based, half-open: `[start, end)`
- Point breakends: `end = start + 1`
- All chromosomes prefixed with `chr`
- Sample names are CCLE names (or DepMap IDs if mapping provided)

## File Counts

```bash
# Combined BEDPE
wc -l sv_wgs_all.bedpe
# Expected: ~11,893 lines (header + data)

# Per-sample files
ls -1 sv_wgs_bedpe | wc -l
# Expected: ~847 files

# Sample with most events
ls -1 sv_wgs_bedpe | head -5
# Example: A549_LUNG.bedpe, HT29_LARGE_INTESTINE.bedpe, ...
```

## Running the Pipeline

### Using sv_bias_pipeline.py (V2)

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

### Using V3 Pipeline Modules

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

## Validation

After conversion, verify:

1. **Coordinate format**: All breakends are point intervals
   ```bash
   awk '{if ($3 != $2+1 || $6 != $5+1) print NR, $0}' sv_wgs_all.bedpe | head
   # Should be empty (no output)
   ```

2. **Chromosome format**: All chromosomes have `chr` prefix
   ```bash
   cut -f1,4 sv_wgs_all.bedpe | grep -v "^chr" | head
   # Should be empty (no output)
   ```

3. **Sample names**: Check sample column
   ```bash
   cut -f11 sv_wgs_all.bedpe | sort -u | head -10
   ```

