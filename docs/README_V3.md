# SV/CN Bias V3 - WGS-based, Pre-registered Case-Study Pipeline

## Overview

V3 replaces CN-derived pseudo-SVs with **real WGS SV breakpoints**, implements **pre-registered** metrics and selection criteria, and produces a **two-part story**: (1) global upper bounds (negative), (2) focused gene case studies with orthogonal validation.

## Pre-registered Ground Rules

### Data Requirements
- **Do not** use CN-derived breakpoints for proximity
- Use real WGS SV breakpoints from BEDPE files (TRA, INV, optional complex)

### Proximity Kernels
All of the following kernels are tested:
- Exponential: `exp(-distance / λ)` with λ ∈ {50k, 100k, 250k, 500k}
- Inverse: `1 / (distance + 1e3)`
- Cubic spline on log10(distance), df=3 (knots at 10th/50th/90th percentiles)

### Models
- **Main model:** Huber regression per gene: `dependency ~ CN + proximity_kernel`
  - Requires n≥200 cells/gene to fit
- **Sensitivity:** ElasticNetCV per gene

### Controls
- **Rotate:** Rotate SV positions within chromosome (preserve per-cell breakpoint counts)
- **Shuffle:** Shuffle across chromosomes (preserve per-cell breakpoint counts and segment lengths)

### Metrics
- **Directional metric:** Prevalence-matched (tertiles/deciles) TRUE vs ROTATE with equal weight per stratum, 10k bootstrap CIs
- **Coef→effect consistency:** Spearman ρ between **mean |prox_contrib|** and **Δ|corr(dep,CN)|** after **prox-only** correction; add robust slope adjusted for prevalence and n
- **Multiple testing:** BH-FDR at q<0.10 for gene-level tests

### Case-Study Inclusion Criteria (all must pass)
1. near-coef percentile ≥95 in TRUE
2. Δ|corr(dep,CN)| > 0 after prox-only correction
3. prevalence-adjusted coef→effect residual in top quartile
4. effect persists in ≥2/3 kernels
5. TRUE residual > ROTATE residual

### Exclusions
- <500 non-NA dep obs
- top 5% CN variance
- DepMap QC blacklist

### Validation (need ≥2/3)
- RNAi concordance improves
- Drug sensitivity improves (matched target)
- Pathway membership (DDR/replication stress)

### Stop/Go Gates
- **Go (Patterns):** ≥20 genes pass, ≥10 at FDR<0.10, ≥3–5 with orthogonal validation
- **Pivot (methods/null):** fewer than above → report upper bounds + framework

## Installation

```bash
# Install dependencies
pip install pandas numpy scipy scikit-learn matplotlib statsmodels

# Optional: for pathway enrichment
pip install gseapy
```

## Data Preparation

### Required Inputs

1. **Dependency data:** `depmap/CRISPRGeneEffect_long.csv`
   - Columns: `gene,cell_line,dependency`

2. **CNV segments:** `depmap/cnv_segments.bed`
   - Columns: `cell_line,chrom,start,end,cn`

3. **WGS SV breakpoints:** `wgs_sv/*.bedpe`
   - Columns: `cell_line,chrom1,start1,end1,chrom2,start2,end2,svtype`
   - Filter for `svtype ∈ {TRA, INV}` (and optional complex)

4. **Gene annotations:** `genes.bed`
   - Columns: `chrom,start,end,gene,strand`

### Optional Inputs

- `rnai_long.csv`: RNAi dependency data (gene,cell_line,dependency)
- `drug_sensitivity.csv`: Drug sensitivity data (compound,cell_line,sensitivity,target)
- Pathway GMT files for enrichment

## Usage

### Step 1: Build WGS Design Matrix

```bash
python sv_ingest_wgs.py \
  --sv-dir wgs_sv/ \
  --genes genes.bed \
  --cnv depmap/cnv_segments.bed \
  --out out_v3/true/design_matrix_wgs.csv
```

This generates `design_matrix_wgs.csv` with columns:
- `gene,cell_line,cn,bp_dist`
- `prox_exp_50k,prox_exp_100k,prox_exp_250k,prox_exp_500k`
- `prox_inv,prox_spline`

### Step 2: Fit Models (TRUE)

```bash
python run_models.py \
  --dep depmap/CRISPRGeneEffect_long.csv \
  --design out_v3/true/design_matrix_wgs.csv \
  --out out_v3/true \
  --model huber \
  --kernel prox_exp_100k \
  --min-cells 200
```

Outputs:
- `model_coefficients.csv`: Per-gene coefficients
- `dependency_corrected.csv`: Full and prox-only corrections
- `metrics.json`: Summary statistics

### Step 3: Generate Controls

```bash
python make_controls.py \
  --sv-dir wgs_sv/ \
  --out out_v3/controls \
  --seed 1
```

This creates:
- `out_v3/controls/rotate/sv_rotated.bedpe`
- `out_v3/controls/shuffle/sv_shuffled.bedpe`

Then re-run Steps 1-2 for each control:
```bash
# Rotate control
python sv_ingest_wgs.py --sv-dir out_v3/controls/rotate/ --genes genes.bed --cnv depmap/cnv_segments.bed --out out_v3/controls/rotate/design_matrix_wgs.csv
python run_models.py --dep depmap/CRISPRGeneEffect_long.csv --design out_v3/controls/rotate/design_matrix_wgs.csv --out out_v3/controls/rotate --kernel prox_exp_100k

# Shuffle control
python sv_ingest_wgs.py --sv-dir out_v3/controls/shuffle/ --genes genes.bed --cnv depmap/cnv_segments.bed --out out_v3/controls/shuffle/design_matrix_wgs.csv
python run_models.py --dep depmap/CRISPRGeneEffect_long.csv --design out_v3/controls/shuffle/design_matrix_wgs.csv --out out_v3/controls/shuffle --kernel prox_exp_100k
```

### Step 4: Global Metrics

```bash
python metrics_global.py \
  --true out_v3/true \
  --rotate out_v3/controls/rotate \
  --shuffle out_v3/controls/shuffle \
  --out out_v3/ \
  --kernel prox_exp_100k \
  --n-bins 10 \
  --n-bootstrap 10000
```

Outputs:
- `prevalence_matched_directional.csv`: Per-bin table
- `coef_effect_summary.csv`: Coef→effect consistency
- `excess_signal.json`: All global metrics

### Step 5: Case Study Selection

```bash
python select_case_studies.py \
  --true out_v3/true \
  --rotate out_v3/controls/rotate \
  --out out_v3/case_studies \
  --kernel prox_exp_100k \
  --fdr-thresh 0.10
```

Outputs:
- `case_studies.csv`: Selected genes with FDR, inclusion flags, validation flags

### Step 6: Orthogonal Validation

```bash
python validate_orthogonal.py \
  --cases out_v3/case_studies/case_studies.csv \
  --true out_v3/true \
  --rnai rnai_long.csv \
  --drug drug_sensitivity.csv \
  --out out_v3/case_studies
```

Outputs:
- `case_studies_validated.csv`: Case studies with validation flags

### Step 7: Generate Panels

```bash
python make_case_panels.py \
  --cases out_v3/case_studies/case_studies_validated.csv \
  --true out_v3/true \
  --sv-dir wgs_sv/ \
  --cnv depmap/cnv_segments.bed \
  --out out_v3/panels
```

Generates 4-panel figures for each case study gene:
- Panel A: Structural context (CN segments + WGS SV breakpoints ±1Mb)
- Panel B: CN vs dependency scatter (before/after correction)
- Panel C: Distance vs dependency (log-x) stratified by CN tertiles
- Panel D: Orthogonal validation mini-panel

### Step 8: Final Report

```bash
python comparison_report.py \
  --root out_v3/ \
  --out out_v3/comparison_report.txt
```

Generates final summary report with:
- Configuration and sample sizes
- Global metrics with CIs
- Upper bounds on global effect
- Case study table (all selected genes)
- Stop/Go gates assessment

## Output Structure

```
out_v3/
├── true/
│   ├── design_matrix_wgs.csv
│   ├── model_coefficients.csv
│   ├── dependency_corrected.csv
│   └── metrics.json
├── controls/
│   ├── rotate/
│   │   ├── design_matrix_wgs.csv
│   │   ├── model_coefficients.csv
│   │   └── dependency_corrected.csv
│   └── shuffle/
│       ├── design_matrix_wgs.csv
│       ├── model_coefficients.csv
│       └── dependency_corrected.csv
├── prevalence_matched_directional.csv
├── coef_effect_summary.csv
├── excess_signal.json
├── case_studies/
│   ├── case_studies.csv
│   └── case_studies_validated.csv
├── panels/
│   ├── GENE1.png
│   ├── GENE2.png
│   └── ...
└── comparison_report.txt
```

## Acceptance Checks

The pipeline includes assertions to ensure:

- ✅ **No** CN-derived breakpoints used in proximity features
- ✅ Directional prevalence-matched reports **bin counts** and **equal-weight average** with CI
- ✅ Coef→effect reports **ρ**, p, **robust slope** with CI, and analysis **restricted to directional genes (Δ>0)**
- ✅ Case-study CSV includes **FDR**, kernel support count, **TRUE>ROTATE** residual check, and validation flags
- ✅ Final report states **upper bound** on global effect and lists **all** selected case genes (not just 3 prettiest)

## Notes

- If WGS SV data source is ambiguous, the pipeline will propose the cleanest public dataset that overlaps with DepMap
- All proximity kernels are computed and stored, but one primary kernel (default: `prox_exp_100k`) is used for main analysis
- Bootstrap CIs use 10,000 iterations by default (configurable)
- FDR threshold is q<0.10 by default (configurable)

## Citation

If using this pipeline, please cite the pre-registered analysis plan and any resulting publications.

