#!/usr/bin/env bash
set -euo pipefail

CFG="configs/v3.yaml"
ts() { date +"%Y-%m-%d %H:%M:%S"; }

echo "[$(ts)] V3 overnight run starting with $CFG"

# Parse a few keys from YAML (python is more portable than yq)
read_var() { python - "$CFG" "$1" <<'PY'
import sys, yaml, pathlib
p = sys.argv[1]; key = sys.argv[2]
cfg = yaml.safe_load(open(p))
def get(d, path):
    for k in path.split("."):
        d = d.get(k, {})
    return d if d != {} else ""
v = get(cfg, key)
if v is None:
    print("")
elif isinstance(v, bool):
    print("True" if v else "False")
elif isinstance(v, list):
    print(" ".join(str(x) for x in v))
else:
    print(str(v))
PY
}

DEP=$(read_var paths.dep)
GENES=$(read_var paths.genes_bed)
CNV=$(read_var paths.cnv_bed)
SVDIR=$(read_var paths.sv_dir)
OUTDIR=$(read_var paths.out_dir)
PCAWG=$(read_var paths.pcawg_index)
OUTPARQ=$(read_var ingest.out_parquet)
SVTYPES=$(read_var ingest.svtypes)
JOIN=$(read_var ingest.join)
DROP=$(read_var ingest.drop_no_sv_lines)
SPLINE=$(read_var ingest.spline)

KERNS=$(read_var models.kernels)
MINCELLS=$(read_var models.min_cells)
MODEL=$(read_var models.model)
SEED=$(read_var models.seed)

DO_ROT=$(read_var shuffle.make_rotate)
DO_WIN=$(read_var shuffle.make_within)
ROT_OUT=$(read_var shuffle.rotate_out)
WIN_OUT=$(read_var shuffle.within_out)

BOOT=$(read_var evaluation.bootstrap_iters)
ALPHA=$(read_var evaluation.alpha)
Q=$(read_var selection.fdr_q)

mkdir -p "$OUTDIR" logs out_v3/rotate out_v3/within

echo "[$(ts)] Using join=${JOIN}, drop_no_sv_lines=${DROP}"

echo "[$(ts)] Optional: annotate CCLE SVs with PCAWG"
if [[ -n "$PCAWG" && -f "$PCAWG" ]]; then
  python scripts/annotate_if_available.py --sv sv_wgs_all.bedpe --pcawg-index "$PCAWG" --out sv_wgs_all_annot.bedpe | tee -a logs/01_annotate.log
  SV_INPUT="sv_wgs_all_annot.bedpe"
else
  SV_INPUT="sv_wgs_all.bedpe"
fi

echo "[$(ts)] Ingest TRUE design"
python sv_ingest_wgs.py \
  --sv-dir "$SVDIR" \
  --genes "$GENES" \
  --cnv "$CNV" \
  --out "$OUTPARQ" \
  --svtypes $(python - <<PY
import yaml,sys
cfg=yaml.safe_load(open("$CFG"))
print(" ".join(cfg["ingest"]["svtypes"]))
PY
) \
  --join "$JOIN" $( [[ "$DROP" == "True" || "$DROP" == "true" ]] && echo "--drop-no-sv-lines" || echo "--keep-no-sv-lines" ) \
  $( [[ "$SPLINE" == "True" || "$SPLINE" == "true" ]] && echo "--spline" ) \
  | tee logs/02_ingest_true.log

echo "[$(ts)] Generate shuffles"
ROTATE_FLAG=""
WITHIN_FLAG=""
[[ "$DO_ROT" == "True" || "$DO_ROT" == "true" ]] && ROTATE_FLAG="--do-rotate"
[[ "$DO_WIN" == "True" || "$DO_WIN" == "true" ]] && WITHIN_FLAG="--do-within"

python scripts/generate_shuffles.py \
  --design "$OUTPARQ" \
  --out-rotate "$ROT_OUT" \
  --out-within "$WIN_OUT" \
  $ROTATE_FLAG $WITHIN_FLAG \
  --seed "$SEED" | tee logs/03_shuffles.log || true

echo "[$(ts)] Fit TRUE models across kernels"
python - <<'PY' "$CFG" "$MODEL" "$MINCELLS"
import sys, yaml, subprocess, os, pathlib
cfg = yaml.safe_load(open(sys.argv[1]))
model = sys.argv[2]
min_cells = sys.argv[3]
kernels = cfg["models"]["kernels"]
design = cfg["ingest"]["out_parquet"]
dep = cfg["paths"]["dep"]
out_root = cfg["paths"]["out_dir"]
for k in kernels:
    outdir = pathlib.Path(out_root) / k
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "run_models.py",
           "--dep", dep, "--design", design, "--out", str(outdir),
           "--kernel", k, "--min-cells", str(min_cells)]
    print("[TRUE]", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)
PY | tee logs/04_fit_true.log

echo "[$(ts)] Fit ROTATE models across kernels"
python - <<'PY' "$CFG" "$MODEL" "$MINCELLS"
import sys, yaml, subprocess, os, pathlib
cfg = yaml.safe_load(open(sys.argv[1]))
model = sys.argv[2]
min_cells = sys.argv[3]
kernels = cfg["models"]["kernels"]
design = cfg["shuffle"]["rotate_out"]
dep = cfg["paths"]["dep"]
out_root = "out_v3/rotate"
for k in kernels:
    outdir = pathlib.Path(out_root) / k
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "run_models.py",
           "--dep", dep, "--design", design, "--out", str(outdir),
           "--kernel", k, "--min-cells", str(min_cells)]
    print("[ROTATE]", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)
PY | tee logs/05_fit_rotate.log

echo "[$(ts)] Evaluate TRUE vs ROTATE and select case studies (prox_exp_100k)"
python scripts/evaluate_and_select.py \
  --true-coefs out_v3/true/prox_exp_100k/model_coefficients.csv \
  --rotate-coefs out_v3/rotate/prox_exp_100k/model_coefficients.csv \
  --alpha-q "$Q" \
  --out-dir out_v3/summary | tee logs/06_evaluate_select.log

echo "[$(ts)] DONE"

