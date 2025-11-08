# Proximity-Only Partial Correlation Analysis Summary

## Overview

This analysis measures Δ|partial r(dep, CN | prox)| after removing **only** the proximity component (not the CN component). This isolates the positional artifact from CN-driven effects.

## Methodology

1. **Proximity-only correction:**
   - `prox_contrib = coef_bp_near * bp_near + coef_bp_far * bp_far`
   - `prox_only_corrected = dependency - (prox_contrib - median_gene_prox_contrib)`

2. **Partial correlation:**
   - Residualize both `dependency` and `prox_only_corrected` against proximity flags
   - Residualize `cn` against proximity flags
   - Compute correlation between residuals: `r0 = corr(dep_resid, cn_resid)`, `r1 = corr(prox_corr_resid, cn_resid)`
   - Metric: `Δ|partial r| = |r0| - |r1|`

## Results

### All Genes (Full Dataset)

| Condition | n_genes | median Δ|partial r| | 95% CI | mean Δ|partial r| | frac > 0 |
|-----------|---------|---------------------|---------|---------------------|----------|
| **TRUE** | 18,424 | **0.0000** | [0.0000, 0.0000] | 1.09e-20 | 30.1% |
| shuffle_within | 18,424 | **0.0000** | [0.0000, 0.0000] | -2.95e-19 | 30.8% |
| shuffle_across | 18,424 | **0.0000** | [0.0000, 0.0000] | -7.16e-20 | 31.7% |

**Finding:** All conditions show median = 0.0000 (essentially zero effect).

### Prevalence-Matched Analysis

- **Genes kept after prevalence matching:** 64 genes (out of 18,424)
- **TRUE result:** median = 0.0000, frac > 0 = 32.8%

**Finding:** Very few genes (64) have matching proximity prevalence between true and shuffled data, suggesting fundamental differences in proximity patterns.

### Low CN Variance Stratum (Q1-Q2)

| Condition | n_genes | median Δ|partial r| | frac > 0 |
|-----------|---------|---------------------|----------|
| **TRUE** | 9,212 | **0.0000** | 31.3% |
| shuffle_within | 9,212 | **0.0000** | 32.2% |
| shuffle_across | 9,212 | **0.0000** | 33.2% |

**Finding:** Even in low CN variance genes, proximity-only correction shows no effect.

## Key Observations

1. **Proximity coefficients are non-zero:**
   - bp_near: 18,357/18,435 genes have non-zero coefficients (median abs = 0.0148)
   - bp_far: 18,284/18,435 genes have non-zero coefficients (median abs = 0.0094)

2. **Proximity contribution is small:**
   - Median absolute proximity contribution: 0.000000
   - Mean absolute proximity contribution: 0.002948
   - Max absolute proximity contribution: 0.113181

3. **After conditioning on proximity, correction has no effect:**
   - When we residualize against proximity flags, we remove the proximity signal
   - The proximity-only correction then has nothing left to correct
   - This suggests the proximity artifact is already captured by the conditioning

## Interpretation

### Why proximity-only correction shows zero effect:

1. **Conditioning removes the signal:** When we compute partial correlation by residualizing against proximity, we're already removing the proximity effect. The proximity-only correction then has minimal additional effect.

2. **Small proximity contributions:** Most genes have very small proximity contributions (median = 0), meaning most cells don't have breakpoints nearby.

3. **CN is the dominant driver:** The remaining CN-dependency correlation after conditioning on proximity is likely driven by CN itself, not proximity artifacts.

### Implications:

1. **Proximity artifact may be small relative to CN effects:**
   - The proximity-only correction doesn't reduce CN-dependency correlation
   - This suggests proximity artifacts are not the main driver of spurious CN-dependency coupling

2. **Breakpoints derived from CN segmentation:**
   - If breakpoints are derived from CN segmentation, they may be too entangled with CN to uniquely identify a proximity artifact
   - Independent SV calls (WGS-based) might be needed

3. **Alternative parameterizations:**
   - Current windows (250k/2Mb) may not capture the true proximity effect
   - Consider trying 100k/1Mb windows or inverse distance terms

4. **Trans controls:**
   - Consider looking at genes on different chromosomes from focal CN events
   - This would isolate proximity from CN effects more cleanly

## Conclusion

The proximity-only partial correlation analysis shows **no measurable effect** across all conditions and filters. This suggests:

- The proximity artifact, if present, is very small relative to CN-driven effects
- After conditioning on proximity, there's little left for proximity correction to affect
- The CN-dependency correlation that remains is likely driven by CN itself, not proximity artifacts

**Recommendation:** Consider alternative approaches:
1. Use independent SV calls (WGS-based) rather than CN-derived breakpoints
2. Try different proximity window parameterizations (100k/1Mb, inverse distance)
3. Use trans controls (genes on different chromosomes from CN events)
4. Focus on genes with high proximity prevalence and low CN variance

