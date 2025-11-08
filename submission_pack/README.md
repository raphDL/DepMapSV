# Submission Pack: SV/CN Proximity Bias Analysis

## Part 1: Pipeline + Falsification
- Robust pipeline with extensive negative controls (rotate/uniform shuffles)
- Gene-level and line-level metrics with EB shrinkage
- Key finding: Proximity coefficients show TRUE > shuffle for near windows, but global Δ|corr| is modest
- Clear limits: Shuffles often look "easier to fix" than real genomes

## Part 2: Focused Signals
- Coefficient→effect consistency: TRUE shows positive ρ (0.17) vs ROTATE negative ρ (-0.07)
- Top-loci permutation: [See top200_perm_summary.csv]
- PSI-by-line analysis identifies genomically unstable lines with elevated proximity sensitivity

## Figures
- `final_directional_bar.*` - Directional correction summary
- `final_directional_excess_bar.*` - TRUE vs ROTATE excess
- `psi_line_compact_figure.png` - PSI by cell line
- `panel_*.png` - Gene×line exemplar panels

## Data Files
- `psi_line_shortlist_for_doc.csv` - Shortlisted high-PSI lines
- `part2_prevalence_matched_win2m.csv` - Prevalence-matched directional (tertiles)
- `coef_effect_param_sweep.csv` - Coefficient→effect robustness
- `top200_perm_summary.csv` - Top-loci permutation results
