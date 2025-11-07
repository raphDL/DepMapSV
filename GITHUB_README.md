# Repository Status: Ready for GitHub

This repository has been cleaned and organized for GitHub publication.

## What's Included

### Code
- ✅ All Python scripts (`.py` files)
- ✅ `environment.yml` for conda environment
- ✅ `Makefile` and shell scripts
- ✅ Documentation (`.md` files)

### Key Results (kept)
- ✅ `out_v2/comparison_report.txt` - Final comparison report
- ✅ `out_v2/*/metrics.json` - Summary metrics with bootstrap CIs
- ✅ `figs_pilot/` - All pilot figures (PNG/PDF)
- ✅ `figs_full/` - Full genome analysis figures

### Excluded (via .gitignore)
- ❌ Large data files (CSV inputs, BED files, GTF files)
- ❌ Large output CSVs (dependency_corrected.csv, model_coefficients.csv)
- ❌ Intermediate outputs (out_pilot_*, out_full_*)
- ❌ Log files
- ❌ Python cache (__pycache__)

## Optional Cleanup (if needed)

If you want to reduce repository size further, you could delete:

```bash
# Remove intermediate pilot runs (keep out_v2/ and figs/)
rm -rf out_pilot_*
rm -rf out_full_*
rm -rf out/
rm -rf chunks/
```

These are already excluded by `.gitignore`, but deleting them locally will reduce local disk usage.

## Initializing Git Repository

If not already initialized:

```bash
git init
git add .gitignore .gitattributes README.md ARTIFACTS.md *.py *.yml *.md Makefile *.sh
git add figs_pilot/ figs_full/
git add out_v2/comparison_report.txt out_v2/*/metrics.json
git commit -m "Initial commit: SV/CN proximity correction project postmortem"
```

## License

Consider adding a LICENSE file if you plan to share this code.

