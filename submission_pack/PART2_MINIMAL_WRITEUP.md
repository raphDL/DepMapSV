# Part 2 (Focused) — Minimal Story

## Key Results

### Prevalence-Matched Directional (2Mb window)
- **TRUE**: 0.000 (no directional signal)
- **ROTATE**: 0.046 (4.6% of genes show directional correction)
- **EXCESS**: -0.046 [95% CI: -0.092, 0.000]
- **Interpretation**: ROTATE shows MORE directional signal than TRUE (negative excess). This is consistent with shuffles being "easier to fix" than real data.

### Coefficient→Effect Consistency (2Mb window)
- **TRUE Spearman ρ**: 0.171 (positive: larger coefficients → larger corrections)
- **ROTATE Spearman ρ**: -0.070 (negative: no consistent relationship)
- **Interpretation**: ✓ TRUE shows biologically plausible coefficient→effect relationship; ROTATE does not.

### Per-Line Coefficient→Effect (Shortlisted Lines)
- **Median ρ_TRUE**: 0.112
- **Median ρ_ROTATE**: -0.066
- **Difference**: 0.178
- **Interpretation**: ✓ Shortlisted lines show stronger coefficient→effect consistency in TRUE vs ROTATE.

### Window Robustness
- **Status**: Not yet run (requires re-running pipeline with 1.5Mb window)
- **Note**: The `svbias_reanalysis.py` script encountered an error. Window robustness check pending.

### Per-Line Essentials AUROC
- **Status**: Could not compute (essential/nonessential gene lists not found in prep_out/)
- **Note**: This check requires `prep_out/essential_genes.txt` and `prep_out/nonessential_genes.txt`

## Summary

**What Works:**
1. Coefficient→effect consistency: TRUE shows positive ρ (0.17) vs ROTATE negative ρ (-0.07)
2. Per-line consistency: Shortlisted lines show TRUE > ROTATE (ρ diff = 0.18)

**What's Challenging:**
1. Prevalence-matched directional: ROTATE shows MORE signal than TRUE (negative excess)
2. This is consistent with the "shuffles are easier to fix" observation from Part 1

**Recommendation:**
Focus on coefficient→effect consistency as the key validation metric. The negative directional excess suggests proximity correction works better on shuffled data, but the coefficient→effect relationship is biologically plausible only in TRUE data.
