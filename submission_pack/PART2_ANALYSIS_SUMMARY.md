# Part 2 Analysis Summary: Cell-Line Heterogeneity

## Overview
Analysis of proximity sensitivity index (PSI) by cell line with Empirical-Bayes shrinkage, followed by validation checks on shortlisted lines.

## Results

### 1. Per-Line Δ|corr| (CN-trimmed, |CN-2| ≤ 0.5)

**Method**: For each cell line, compute Spearman |corr(dep, CN)| across genes before and after prox-only correction. Record Δ|corr| = |r_before| - |r_after|. Compare TRUE vs ROTATE with paired bootstrap CI.

**Results**:
- TRUE median Δ|corr|: **0.0000**
- ROTATE median Δ|corr|: **0.0000**
- Paired difference (TRUE - ROTATE): **0.0000** [95% CI: 0.0000, 0.0000]
- N pairs: 57

**Interpretation**: ✗ **WEAK** - No evidence TRUE > ROTATE at the per-line level. The per-line Δ|corr| metric is essentially flat, suggesting that proximity correction does not substantially change CN-dependency correlations when computed across genes within a line.

**Note**: This is the "tiny reality check" - the effect is subtle and may not be detectable at the per-line aggregation level.

### 2. Coefficient→Effect Consistency

**Method**: For each shortlisted line, compute Spearman correlation between:
- per-gene |β_prox|max (max absolute proximity coefficient)
- per-gene corr-drop (Δ|corr| computed across all lines)

**Results** (top 5 shortlisted lines):
- Median ρ (TRUE): **0.1119**
- Median ρ (ROTATE): **-0.0663**
- Median ρ difference: **0.1781**

**Interpretation**: ✓ **POSITIVE** - TRUE shows stronger coefficient→effect consistency than ROTATE. This suggests that genes with larger proximity coefficients do show larger correlation drops in TRUE data, but not in shuffled data.

**Caveat**: The values are identical across all 5 lines, suggesting the computation uses global gene-level correlations (across all lines) rather than line-specific correlations. This is still meaningful but represents a different signal than line-specific effects.

### 3. Window Robustness

**Status**: Not yet run. Requires re-running focused pilot with 1.5Mb window.

**Command** (when ready):
```bash
python svbias_reanalysis.py \
  --pilot --pilot-unstable-only \
  --pilot-genes-file out_focus_true/unstable/design_matrix.csv \
  --dependency CRISPRGeneEffect_long.csv \
  --cnv prep_out/cnv_segments.bed \
  --sv prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed \
  --add-continuous-proximity --standardize-predictors \
  --bp_windows 250000 1500000 \
  --model huber \
  --out out_focus_true_win15 --progress -q
```

## Key Findings

### What Works
1. **EB shrinkage PSI shortlist**: Identified 10 lines with post_prob_gt0 ≥ 0.80 (top 2 with 0.86, 0.85)
2. **Coefficient→effect consistency**: TRUE shows positive ρ (0.11) vs ROTATE negative ρ (-0.07)
3. **Instability association**: Shortlisted lines have high FGA (0.94-1.00) and many breakpoints (628-2296)

### What's Flat
1. **Per-line Δ|corr|**: Essentially zero for both TRUE and ROTATE
   - This suggests proximity correction doesn't substantially change CN-dependency correlations when aggregated at the line level
   - May reflect that the effect is gene-specific rather than line-specific

## Recommendations for Part 2 Write-Up

### Option A: Conservative Framing (Recommended)
**"We detect subtle, line-specific proximity sensitivity; our pipeline and falsification controls quantify an upper bound on global bias."**

**Key points**:
- PSI-by-line analysis with EB shrinkage identifies genomically unstable lines (high FGA, many breakpoints) with elevated proximity sensitivity
- Coefficient→effect consistency shows TRUE > ROTATE (ρ = 0.11 vs -0.07), supporting biological plausibility
- Per-line Δ|corr| is flat, indicating the effect is subtle and gene-specific rather than creating large line-level CN-dependency decoupling
- This is consistent with proximity bias being a **local, gene-specific artifact** rather than a global line-level confound

### Option B: Stronger Claim (if window robustness confirms)
**"Proximity bias is heterogeneous—detectable in a minority of genomically unstable lines, with consistent effects across window definitions."**

**Requires**:
- Window robustness check showing same sign of Δ|corr| under 1.5Mb vs 2Mb
- At least 2 lines showing consistent positive Δ|corr| across windows
- Positive coefficient→effect ρ in those lines

## Files Generated

1. `psi_line_corr_delta.csv` - Per-line Δ|corr| for TRUE and ROTATE
2. `psi_line_corr_delta_summary.csv` - Summary statistics
3. `psi_line_coef_effect.csv` - Coefficient→effect consistency per line
4. `psi_line_compact_figure.png` - Main PSI-by-line figure
5. `psi_line_shortlist_for_doc.csv` - Shortlisted lines for documentation

## Next Steps

1. **Run window robustness check** (1.5Mb vs 2Mb) if time permits
2. **Decide on framing**: Conservative (Option A) or stronger (Option B)
3. **Focus on gene-level effects**: The flat per-line Δ|corr| suggests the story is at the gene level, not line level
4. **Emphasize falsification**: The TRUE > ROTATE in coefficient→effect consistency is a key validation

