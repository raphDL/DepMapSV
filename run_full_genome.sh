#!/usr/bin/env bash
set -euo pipefail

# Set environment
export MPLCONFIGDIR="$PWD/.mplconfig"
export XDG_CACHE_HOME="$PWD/.cache"
export PYTHONHASHSEED=1
export OMP_NUM_THREADS=1
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

# Choose robust setting (recommended: 250k/1.5Mb)
WIDE="250000 1500000"

# Create output directories
mkdir -p out_full_linear out_full_linear_rotate

# Non-rotated
echo "Running non-rotated chunks..."
for f in chunks/fullgenes_*; do
  tag=$(basename "$f")
  outdir="out_full_linear/$tag"
  mkdir -p "$outdir"
  echo "Processing $tag (non-rotated)..."
  python -u svbias_reanalysis.py \
    --pilot --pilot-size 1000 --pilot-unstable-only \
    --pilot-genes-file "$f" \
    --dependency CRISPRGeneEffect.csv \
    --cnv prep_out/cnv_segments.bed \
    --sv  prep_out/sv_from_cnv.bedpe \
    --genes prep_out/genes.depmap.unique.bed \
    --model linear --bp_windows $WIDE \
    --add-continuous-proximity --standardize-predictors \
    --out "$outdir" --progress -q --log-file "$outdir/reanalysis.log"
done

# Rotate shuffle (pairwise null)
echo "Running rotate shuffle chunks..."
for f in chunks/fullgenes_*; do
  tag=$(basename "$f")
  outdir="out_full_linear_rotate/$tag"
  mkdir -p "$outdir"
  echo "Processing $tag (rotate shuffle)..."
  python -u svbias_reanalysis.py \
    --pilot --pilot-size 1000 --pilot-unstable-only \
    --pilot-genes-file "$f" \
    --dependency CRISPRGeneEffect.csv \
    --cnv prep_out/cnv_segments.bed \
    --sv  prep_out/sv_from_cnv.bedpe \
    --genes prep_out/genes.depmap.unique.bed \
    --model linear --bp_windows $WIDE \
    --add-continuous-proximity --standardize-predictors \
    --shuffle-rotate \
    --out "$outdir" --progress -q --log-file "$outdir/reanalysis.log"
done

echo "Full-genome analysis complete. Run summarize_full_genome.py to compute directional flags and summary."


