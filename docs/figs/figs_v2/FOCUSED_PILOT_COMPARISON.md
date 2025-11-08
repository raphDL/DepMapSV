# Focused Pilot Comparison: TRUE vs Rotate-Shuffle

## Overview

Focused pilot analysis on top 150 excess candidate genes using richer model specification:
- Linear model
- Continuous proximity term (1/(bp_dist+1e3))
- Standardized predictors
- BP windows: 250k, 2Mb

## Key Metrics Comparison

| Metric | TRUE | ROTATE | Winner |
|--------|------|--------|--------|
| **Median Δ|corr|** | **0.0054** | 0.0053 | ✅ TRUE (barely) |
| **Prox-only median Δ|corr|** | -0.0164 | **-0.0035** | ❌ ROTATE |
| **Directional prox-active frac** | 0.293 | **0.457** | ❌ ROTATE |
| Prox-active genes frac | 0.933 | 0.933 | Tie |
| Nonzero proximity coeffs frac | 0.933 | 0.933 | Tie |

## Detailed Results

### TRUE Condition
- **Genes fitted:** 150
- **Median Δ|corr|:** 0.0054
- **Prox-only median Δ|corr|:** -0.0164 (negative = correction increases correlation)
- **Directional prox-active frac:** 0.293 (29.3% of genes show directional correction)
- **Proximity prevalence:** p100k=0.0351, p1m=0.1793
- **Top genes by corr drop:** CCN3 (0.1046), TLR2 (0.0429), BCAP31 (0.0401), HLA-DQB1 (0.0363), MROH2B (0.0342)

### ROTATE Condition
- **Genes fitted:** 150
- **Median Δ|corr|:** 0.0053
- **Prox-only median Δ|corr|:** -0.0035 (less negative = better)
- **Directional prox-active frac:** 0.457 (45.7% of genes show directional correction)
- **Proximity prevalence:** p100k=0.0628, p1m=0.3062 (higher than TRUE)
- **Top genes by corr drop:** CCN3 (0.1021), TLR2 (0.0424), BCAP31 (0.0401), MROH2B (0.0354), SCFD2 (0.0315)

## Interpretation

### Mixed Results

1. **Median Δ|corr|:** TRUE wins by a tiny margin (0.0054 vs 0.0053)
   - This is essentially a tie
   - Suggests the richer model spec doesn't help TRUE more than shuffles

2. **Prox-only Δ|corr|:** ROTATE wins decisively (-0.0035 vs -0.0164)
   - TRUE shows more negative value (correction increases correlation)
   - ROTATE shows less negative value (correction works better)
   - This suggests proximity-only correction works better on shuffled data

3. **Directional prox-active frac:** ROTATE wins decisively (0.457 vs 0.293)
   - ROTATE shows 45.7% of genes with directional correction
   - TRUE shows only 29.3%
   - This suggests shuffles are easier to "fix" than real data

## Conclusion

**ROTATE ≥ TRUE on 2/3 key metrics**

The focused pilot on top excess candidate genes shows that:
- Even with richer model specification (linear, continuous proximity, standardized predictors)
- Even on genes selected for showing excess signal
- **Rotate-shuffle still performs as well or better than TRUE**

This strongly suggests:
1. **Shuffled worlds are easier to "fix"** - Shuffles create simpler artifacts that are easier to correct
2. **Real-world complexity** - TRUE data has more complex structure that is harder to model
3. **Model limitations** - Current model specification may not capture the true artifact structure

## Recommendations

### Option 1: Narrow-Scope Paper
- Document subtle proximity effects where they matter
- Focus on specific genes/contexts where TRUE > shuffles
- Acknowledge limitations and mixed results

### Option 2: Pivot
- Consider alternative approaches:
  - Independent SV calls (WGS-based) instead of CN-derived breakpoints
  - Trans controls (genes on different chromosomes from CN events)
  - Different proximity parameterizations (100k/1Mb, inverse distance)
- Focus on understanding why shuffles outperform TRUE

### Option 3: Deep Investigation
- Analyze why specific genes (CCN3, TLR2, BCAP31) show large correlation drops in both conditions
- Investigate the structure of artifacts in TRUE vs shuffles
- Consider whether the correction is removing real signal rather than artifacts

## Files Generated

- `out_focus_true/unstable/` - TRUE condition outputs
- `out_focus_rotate/unstable/` - ROTATE condition outputs
- Both contain: `design_matrix.csv`, `models_coefficients.csv`, `pilot_summary.txt`, `dependency_corrected.csv`

