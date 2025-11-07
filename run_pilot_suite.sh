#!/usr/bin/env bash
set -euo pipefail

# Silence matplotlib cache warnings
export MPLCONFIGDIR="$PWD/.mplconfig"
export XDG_CACHE_HOME="$PWD/.cache"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME" figs_pilot

# Lock unstable genes from baseline
BASE_UNSTABLE="out_pilot_final/unstable"
cut -d, -f1 "$BASE_UNSTABLE/design_matrix.csv" | tail -n +2 | sort -u > pilot_unstable_genes.txt

# Baseline fixed (Huber, 100k/1Mb)
mkdir -p out_pilot_baseline_fixed
python -u svbias_reanalysis.py \
  --pilot --pilot-size 50 --pilot-unstable-only \
  --pilot-genes-file pilot_unstable_genes.txt \
  --dependency CRISPRGeneEffect.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv  prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --add-continuous-proximity --standardize-predictors \
  --out out_pilot_baseline_fixed \
  --progress -q --log-file out_pilot_baseline_fixed/reanalysis.log
python score_hits.py out_pilot_baseline_fixed/unstable --ignore-lineage

# Rotate shuffle fixed
mkdir -p out_pilot_shuffle_rotate
python -u svbias_reanalysis.py \
  --pilot --pilot-size 50 --pilot-unstable-only \
  --pilot-genes-file pilot_unstable_genes.txt \
  --dependency CRISPRGeneEffect.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv  prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --add-continuous-proximity --standardize-predictors \
  --shuffle-rotate \
  --out out_pilot_shuffle_rotate \
  --progress -q --log-file out_pilot_shuffle_rotate/reanalysis.log
python score_hits.py out_pilot_shuffle_rotate/unstable --ignore-lineage

# Uniform shuffle fixed
mkdir -p out_pilot_shuffle_uniform
python -u svbias_reanalysis.py \
  --pilot --pilot-size 50 --pilot-unstable-only \
  --pilot-genes-file pilot_unstable_genes.txt \
  --dependency CRISPRGeneEffect.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv  prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --add-continuous-proximity --standardize-predictors \
  --shuffle-uniform \
  --out out_pilot_shuffle_uniform \
  --progress -q --log-file out_pilot_shuffle_uniform/reanalysis.log
python score_hits.py out_pilot_shuffle_uniform/unstable --ignore-lineage

# Robustness: linear + (250k,2Mb) fixed
mkdir -p out_pilot_lin_250k_2m_fixed
python -u svbias_reanalysis.py \
  --pilot --pilot-size 50 --pilot-unstable-only \
  --pilot-genes-file pilot_unstable_genes.txt \
  --dependency CRISPRGeneEffect.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv  prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --bp_windows 250000 2000000 --model linear \
  --add-continuous-proximity --standardize-predictors \
  --out out_pilot_lin_250k_2m_fixed \
  --progress -q --log-file out_pilot_lin_250k_2m_fixed/reanalysis.log
python score_hits.py out_pilot_lin_250k_2m_fixed/unstable --ignore-lineage

# Optional: linear/250k–2Mb + rotate shuffle on fixed genes
mkdir -p out_pilot_lin_250k_2m_rotate
python -u svbias_reanalysis.py \
  --pilot --pilot-size 50 --pilot-unstable-only \
  --pilot-genes-file pilot_unstable_genes.txt \
  --dependency CRISPRGeneEffect.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv  prep_out/sv_from_cnv.bedpe \
  --bp_windows 250000 2000000 --model linear \
  --genes prep_out/genes.depmap.unique.bed \
  --add-continuous-proximity --standardize-predictors \
  --shuffle-rotate \
  --out out_pilot_lin_250k_2m_rotate \
  --progress -q --log-file out_pilot_lin_250k_2m_rotate/reanalysis.log
python score_hits.py out_pilot_lin_250k_2m_rotate/unstable --ignore-lineage

# Linear + 250k/1.5Mb (fixed genes)
mkdir -p out_pilot_lin_250k_1p5m_fixed
python -u svbias_reanalysis.py \
  --pilot --pilot-size 50 --pilot-unstable-only \
  --pilot-genes-file pilot_unstable_genes.txt \
  --dependency CRISPRGeneEffect.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv  prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --bp_windows 250000 1500000 --model linear \
  --add-continuous-proximity --standardize-predictors \
  --out out_pilot_lin_250k_1p5m_fixed \
  --progress -q --log-file out_pilot_lin_250k_1p5m_fixed/reanalysis.log
python score_hits.py out_pilot_lin_250k_1p5m_fixed/unstable --ignore-lineage

# Linear + 250k/1.5Mb + rotate (fixed genes)
mkdir -p out_pilot_lin_250k_1p5m_rotate
python -u svbias_reanalysis.py \
  --pilot --pilot-size 50 --pilot-unstable-only \
  --pilot-genes-file pilot_unstable_genes.txt \
  --dependency CRISPRGeneEffect.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv  prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --bp_windows 250000 1500000 --model linear \
  --add-continuous-proximity --standardize-predictors \
  --shuffle-rotate \
  --out out_pilot_lin_250k_1p5m_rotate \
  --progress -q --log-file out_pilot_lin_250k_1p5m_rotate/reanalysis.log
python score_hits.py out_pilot_lin_250k_1p5m_rotate/unstable --ignore-lineage

# Final figure and stats
python make_final_directional_fig.py

# Echo summary
python - << 'PY'
import re
paths = {
  "baseline_fixed": "out_pilot_baseline_fixed/unstable/pilot_summary.txt",
  "rotate":         "out_pilot_shuffle_rotate/unstable/pilot_summary.txt",
  "uniform":        "out_pilot_shuffle_uniform/unstable/pilot_summary.txt",
  "linear250k2m":   "out_pilot_lin_250k_2m_fixed/unstable/pilot_summary.txt",
}

def grab(p):
    d={}
    with open(p) as f:
        for ln in f:
            if "prox_active_genes_frac_(" in ln and "_&_drop>0" not in ln:
                d["prox_active"] = float(re.findall(r"=([0-9.]+)$", ln)[0])
            if "prox_active_genes_frac_dirpos_(" in ln:
                d["prox_active_dir"] = float(re.findall(r"=([0-9.]+)$", ln)[0])
    return d

for k,v in paths.items():
    d = grab(v)
    print(k, "→ prox_active=", d.get("prox_active"), " prox_active_dir=", d.get("prox_active_dir"))
PY
