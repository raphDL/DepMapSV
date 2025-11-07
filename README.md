# SV/CN Proximity Correction – Project Postmortem

## Decision

**Stop**. Multiple falsification tests show no positive, positionally-specific signal attributable to true SV proximity; shuffled controls outperform the true data on the directional metric. Correction yields no measurable gain on essential/non-essential separation.

## What we attempted

- Per-gene robust regression (CN + proximity windows) on DepMap 25Q3.
- Directional prox-active metric, pilot → full-dataset.
- Two negative controls: within-chrom and cross-chrom shuffles.
- Full-dataset run with bootstrap CIs; coefficient prevalence, corr-drop distributions.
- Essential vs non-essential AUROC/AUPRC (gene-level medians).

## Key results

- **Directional prox-active (true)** < **shuffled** (within & across chrom): negative excess.
- **Essential/non-essential**: AUROC ≈ 1.0 pre-correction; **ΔAUROC ≈ 0** after correction; **ΔΔAUROC CI overlaps 0**.
- Full-dataset directional fraction (true) ~0.31; shuffled ~0.53–0.57 (numbers from run logs).
- Pilot and fixed-gene audits replicated the pattern; widening windows boosts all conditions, not just true.

## Interpretation

- Shuffling breaks structure yet increases the directional metric → our metric/model captures variance that is **easier** (more "artifact-like") in shuffled data than in real SV geography.
- With DepMap-derived essential sets, baseline already perfectly separates labels (partly circular), so there's **no headroom** for improvement.

## What's salvageable

- Clean, documented pipeline pieces:
  - Data validation + manifesting
  - Feature builders (CN, SV proximity)
  - Shuffling frameworks (within & cross-chrom)
  - Pilot mode + design matrices
  - Figure and audit scripts (directional flags, CIs)
- These are reusable for other assays/questions (e.g., enhancer proximity, CNA-aware QC).

## If revisited later (out of scope now)

- Truly orthogonal validation (e.g., external screens or pathway/complex recovery not derived from DepMap).
- Alternative causal designs (local randomization around "natural experiments" like focal events).
- Different targets (e.g., **trans** contacts from Hi-C instead of linear proximity).

## Artifacts saved

- `out_v2/comparison_report.txt` summarizing true vs shuffles
- Final directional bar/excess plots and stats (`figs_pilot/`, `figs_full/`)
- Full-run metrics JSON + per-gene coefficients (`out_v2/`)
- Pilot summaries and case panels

## Code structure

### Main pipelines
- `sv_bias_pipeline.py` - Final v2.1 pipeline with bootstrap CIs and negative controls
- `svbias_reanalysis.py` - Original reanalysis pipeline with pilot mode

### Data preparation
- `prep_depmap_25Q3.py` - Converts DepMap CNV segments and expression to standardized formats
- `gtf_to_gene_bed.py` - Converts GTF annotations to BED format

### Analysis and visualization
- `make_pilot_figs.py` - Generate pilot figures
- `make_final_directional_fig.py` - Generate directional metric figures
- `make_case_studies.py` - Generate case study panels
- `audit_directional.py` - Audit directional flags
- `summarize_full_genome.py` - Summarize full genome results
- `score_hits.py` - Score and analyze hits

### Utilities
- `pilot_bootstrap_ci.py` - Bootstrap confidence intervals for pilot metrics
- `hits_summary.py` - Summarize hit analysis
- `add_lineage_to_design.py` - Add lineage information to design matrices

## Installation

```bash
conda env create -f environment.yml
conda activate svbias
```

## Usage example (for reference)

```bash
# Prepare data
python prep_depmap_25Q3.py \
  --cnv_segments_wgs OmicsCNSegmentsWGS.csv \
  --expression_tpm OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv \
  --outdir prep_out \
  --cn_assume_diploid

python gtf_to_gene_bed.py gencode.v49.basic.annotation.gtf prep_out/gencode.v49.gene.bed

# Run pipeline (final v2.1)
python sv_bias_pipeline.py \
  --dependency CRISPRGeneEffect_long.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --model huber --bp-windows 250000 2000000 \
  --bootstrap-iterations 2000 \
  --activity-threshold 0.01 --activity-cell-fraction 0.10 \
  --output-dir out_v2
```

## Notes

- Large data files (CSV inputs, outputs) are excluded via `.gitignore`
- Key results are in `out_v2/` (final run) and `figs_pilot/`, `figs_full/` (visualizations)
- See individual script `--help` for detailed options
