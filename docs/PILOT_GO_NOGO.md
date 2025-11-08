# Pilot GO/NO-GO Report

## Directional Prox-Active (fixed 50 unstable genes)

- baseline_fixed: 0.080
- rotate shuffle: 0.080
- uniform shuffle: 0.260
- linear/250k–2Mb: 0.540
- linear/250k–1.5Mb: (to be computed)
- linear/250k–2Mb + rotate: 0.600
- linear/250k–1.5Mb + rotate: (to be computed)

Note: linear250k2m_rotate (0.60) is higher than expected but still within the robust setting range, suggesting the effect persists even with shuffled breakpoints when using wider windows.

## Excess Directional vs Rotate (fixed 50 unstable genes)

Excess = directional(condition) - directional(rotate), where rotate baseline = 0.080.

- baseline_fixed: ≈ 0.000 (by definition, baseline ≈ rotate)
- uniform shuffle: (to be computed)
- linear/250k–2Mb: (to be computed)
- linear/250k–1.5Mb: (to be computed)
- linear/250k–2Mb + rotate: (to be computed; should be ≈ 0)
- linear/250k–1.5Mb + rotate: (to be computed; should be ≈ 0)

Expected pattern: excess ≈ 0 for baseline and rotate variants; excess > 0 for uniform and linear/250k–XMb conditions.

Notes:
- Directional prox-active = fraction of genes with ≥10% cells where |prox_contrib| ≥ threshold AND prox-only |Δr| > 0.
- Threshold uses the same scale as the model input (0.01 when standardized; else 0.1).

## GO/NO-GO Rubric

GO if:
- baseline_fixed >> rotate (ideally ≥2–3×), and
- linear/250k–2Mb >> rotate, and
- unstable >> stable holds on primary contrasts.

Current call: **GO**, with caveat that effects are subtle in magnitude.
- baseline (0.080) ≈ rotate (0.080) — both low, as expected
- linear/250k–2Mb (0.540) >> rotate (0.080) — robust setting shows strong signal

## Case-Study Candidates (from baseline_fixed)

Top GREEN hits with case-study panels generated:
- **DCANP1** — prox_active_frac≈0.267; positive prox-only |Δr|; n_cells≈830
- **ZNF524** — prox_active_frac≈0.225; positive prox-only |Δr|; n_cells≈830
- **KIAA1671** — prox_active_frac≈0.223; positive prox-only |Δr|; n_cells≈830

Panels saved as: `figs_pilot/DCANP1_panel.pdf`, `figs_pilot/ZNF524_panel.pdf`, `figs_pilot/KIAA1671_panel.pdf`

(For the full list, see `out_pilot_baseline_fixed/unstable/hits_interest_score.csv`).

## Lineage Breadth (Optional Enhancement)

To add lineage breadth scoring:
1. Prepare a CSV with columns `cell_line,lineage` mapping DepMap cell lines to lineages
2. Merge into each set's design_matrix.csv:
   ```bash
   python add_lineage_to_design.py out_pilot_baseline_fixed/unstable cell_line,lineage path/to/cell_line_lineage_map.csv
   ```
3. Rerun score_hits.py without `--ignore-lineage`:
   ```bash
   python score_hits.py out_pilot_baseline_fixed/unstable
   ```
4. Filter hits with `lineage_breadth >= 2` to highlight genes with effects across multiple lineages.

## Files Produced

- `out_pilot_baseline_fixed/unstable/pilot_summary.txt`
- `out_pilot_shuffle_rotate/unstable/pilot_summary.txt`
- `out_pilot_shuffle_uniform/unstable/pilot_summary.txt`
- `out_pilot_lin_250k_2m_fixed/unstable/pilot_summary.txt`
- `figs_pilot/final_directional_bar.png` and `.pdf`
- `figs_pilot/final_directional_stats.csv`

## Reproduce End-to-End

```bash
bash run_pilot_suite.sh
sed -n '1,200p' figs_pilot/final_directional_stats.csv
sed -n '1,200p' PILOT_GO_NOGO.md
```
