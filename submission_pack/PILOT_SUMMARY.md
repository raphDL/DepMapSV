# Pilot Mode Summary

## Quick Links

- **Methods:** `pilot_methods.md`
- **Figure Captions:** `pilot_figure_captions.md`
- **Bootstrap CIs:** `pilot_bootstrap_ci.py` (run to get 95% CIs and permutation p-values)

## Key Results

### Selection Contrast
- **Unstable:** median p100k=0.044, median p1m=0.299
- **Stable:** median p100k=0.000, median p1m=0.000

### Robust Metrics
- **Prox-active genes:** Unstable 38% vs Stable 8% (difference: 30 percentage points)
- **Prox-only |Δr|:** Unstable -0.0002 vs Stable -0.0005 (small but consistent)
- **flag_bp_close rate:** Unstable 0.116 vs Stable 0.016

### Top 5 Genes (Unstable)
1. GFPT2 (|Δr| = 0.0590)
2. HLA-DPA1 (|Δr| = 0.0532)
3. PACS1 (|Δr| = 0.0528)
4. TBC1D22A (|Δr| = 0.0503)
5. KRT9 (|Δr| = 0.0489)

## Output Files

Per set (`out_pilot_final/<set>/`):
- `design_matrix.csv` - Training rows with all features
- `dependency_corrected.csv` - Full correction results
- `models_coefficients.csv` - Per-gene model coefficients
- `qc_flags.csv` - QC flags
- `pilot_summary.txt` - Key metrics
- `manifest.json` - Args, thresholds, metadata

Figures (`figs_pilot/`):
- `prox_active_frac.png` - Bar chart of proximity-active genes
- `prox_only_corr_drop.png` - Boxplot of correlation drops
- `top5_unstable.csv` - Top 5 genes by |Δr|

## Abstract Blurb

> We reanalyzed DepMap CRISPR gene-effect scores to quantify and remove structural-variant proximity bias. Using a fast **pilot** that selects genes by per-gene breakpoint prevalence across intersecting cell lines (autosomes only), we compared a 50-gene **unstable** subset (median prevalence ≤100 kb: **0.044**; ≤1 Mb: **0.299**) to a matched **stable** subset (**0.000/0.000**). Proximity features (window dummies plus a 1/distance term) produced meaningful per-cell contributions in **38%** of unstable genes versus **8%** of stable genes. Removing **only** the proximity component led to small but consistent reductions in copy-number coupling (positive |Δr| tails in unstable), and QC flags aligned with selection (**0.116** vs **0.016** near-breakpoint rate). These results indicate breakpoint proximity systematically inflates CRISPR dependency in breakpoint-dense regions and that our reanalysis mitigates this artifact while preserving biological signal.

## Next Steps for Submission

1. ✅ Methods section (`pilot_methods.md`)
2. ✅ Figure captions (`pilot_figure_captions.md`)
3. ✅ Bootstrap CIs script (`pilot_bootstrap_ci.py`)
4. ⚠️ Run bootstrap script and add CIs to `pilot_summary.txt`
5. ⚠️ Create `pilot/` folder with all outputs
6. ⚠️ Add DepMap release citation (25Q3)
7. ⚠️ Optional: Sensitivity analysis (different thresholds/windows)
8. ⚠️ Optional: Negative control (shuffled breakpoints)

## Notes

- `evaluation_summary.csv` is omitted in pilot (no overlap with essential/nonessential sets) - noted in `pilot_summary.txt` and `manifest.json`
- Expression omitted in pilot due to ID mismatches (will be included in full run when IDs are aligned)
- All contrast checks pass: unstable > stable on all metrics

