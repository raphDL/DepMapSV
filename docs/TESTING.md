# Testing Guide for sv_ingest_wgs.py

## Quick Test on Synthetic Data

### Step 1: Generate Test Data
```bash
python3 test_synthetic.py
```

### Step 2: Run Pipeline on Synthetic Data
```bash
python3 sv_ingest_wgs.py \
  --sv-dir test_data/ \
  --genes test_data/genes.bed \
  --cnv test_data/cnv.csv \
  --out test_data/design.csv \
  --svtypes TRA INV \
  --cap-bp 5000000 \
  --dedupe-bp-tol 100 \
  -vv
```

### Step 3: Validate Output
```bash
python3 validate_design.py test_data/design.csv
```

## Expected Results (Synthetic Test)

- **Overlap**: 2 cell lines
- **Design matrix**: 8 genes × 2 cells = 16 rows
- **GENE_A**: bp_dist ~ 100kb (SV at 1.1Mb, gene at 1.0Mb)
- **GENE_B**: bp_dist ~ 50kb (SV at 2.05Mb, gene at 2.0Mb)
- **GENE_H**: bp_dist >> 1Mb (no nearby SVs on chr2)

## Real Data Testing

### Step 1: Dry Run
```bash
python3 sv_ingest_wgs.py \
  --sv-dir /path/to/depmap/wgs_sv_bedpe/ \
  --genes gencode.v49.gene.bed \
  --cnv OmicsCNSegmentsWGS.csv \
  --sample-info Model.csv \
  --map-cnv-to sv \
  --out design_matrix.csv \
  --dry-run \
  -vv 2>&1 | tee logs/dry_run.log
```

**Check logs for:**
- Overlap: Should be 50-200 (realistic WGS SV availability)
- If overlap = 0, check `logs/overlap_debug_samples.csv`

### Step 2: Full Run (only after dry-run passes!)
```bash
mkdir -p logs

python3 sv_ingest_wgs.py \
  --sv-dir /path/to/depmap/wgs_sv_bedpe/ \
  --genes gencode.v49.gene.bed \
  --cnv OmicsCNSegmentsWGS.csv \
  --sample-info Model.csv \
  --map-cnv-to sv \
  --out design_matrix.csv \
  --svtypes TRA INV \
  --cap-bp 5000000 \
  --dedupe-bp-tol 100 \
  --drop-no-sv-lines \
  -vv 2>&1 | tee logs/sv_ingest_full.log
```

### Step 3: Validate Real Data Output
```bash
python3 validate_design.py design_matrix.csv
```

## Stop/Go Decision Matrix

### ✅ GREEN LIGHT - Proceed with Modeling
- Overlap: 50-200 cell lines ✓
- CN-proximity max |r| < 0.3 ✓
- Kernel non-zero %: 70-99% ✓
- Gene coverage: >70% genes have SV data ✓
- CN distribution: Mean 1.5-2.5, range 0.5-6.0 ✓

### ⚠️ YELLOW LIGHT - Investigate but Maybe Proceed
- Overlap: 30-50 cell lines (limited power)
- CN-proximity max |r|: 0.3-0.4 (some entanglement)
- Kernel non-zero %: 50-70% (sparse breakpoints)
- Gene coverage: 50-70% (many genes missing)

**Action**: Document limitations, proceed with case-study approach only

### 🛑 RED LIGHT - Do Not Proceed
- Overlap: <30 cell lines (insufficient data)
- CN-proximity max |r| > 0.4 (still entangled)
- Kernel non-zero % < 50% (too sparse)
- Gene coverage < 50% (inadequate)
- CN distribution: Mean <1.0 or >4.0 (wrong scaling)

**Action**: Fix data issues, revisit approach, or abandon analysis

