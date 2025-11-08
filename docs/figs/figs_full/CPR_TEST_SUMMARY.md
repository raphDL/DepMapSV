# CPR Test Results Summary
**Date:** Generated from `out_v2/` outputs  
**Purpose:** Decisive, low-effort tests to assess actionable signal beyond pilot

---

## Executive Summary

**Score: 2/3 tests PASSED** ⚠️

- ✅ **Test B:** Proximity coefficient signal (TRUE > shuffles)
- ✅ **Test C:** Nonessential de-essentialization (TRUE > shuffles)
- ❌ **Test A:** Correlation drop at Q4 (TRUE < shuffles) — **CRITICAL FAILURE**

---

## Detailed Results

### Test A: Δ|corr(dep,CN)| Stratified by Baseline |r| Quartiles

**Hypothesis:** Correction should show larger improvement where baseline coupling is strongest (Q4).

**Result:** ❌ **FAILED**

| Quartile | TRUE median Δ|corr| | shuffle_within | shuffle_across |
|----------|--------------|----------------|----------------|
| Q1 (lowest \|r\|) | -0.0040 | -0.0035 | -0.0043 |
| Q2 | 0.0099 | 0.0146 | 0.0159 |
| Q3 | 0.0190 | 0.0283 | 0.0311 |
| **Q4 (highest \|r\|)** | **0.0422** | **0.0589** | **0.0640** |

**Interpretation:** 
- TRUE shows **smaller** correlation drop at Q4 than both shuffles
- This suggests the correction is **underperforming** compared to random for high baseline correlation genes
- This is the most concerning finding — where we expect the biggest impact, we see the worst relative performance

---

### Test B: Proximity Coefficient Distribution Tails

**Hypothesis:** Real positional signal should produce larger proximity coefficients in TRUE than shuffles.

**Result:** ✅ **PASSED**

| Metric | TRUE | shuffle_within | shuffle_across |
|--------|------|----------------|----------------|
| median \|prox coef\| | 0.0186 | 0.0142 | 0.0158 |
| **p90 \|prox coef\|** | **0.0548** | **0.0405** | **0.0440** |
| p99 \|prox coef\| | 0.1635 | 0.1111 | 0.1200 |

**Interpretation:**
- TRUE shows consistently larger proximity coefficients across all percentiles
- p90 and p99 are substantially higher in TRUE, suggesting real positional signal
- This indicates the model is detecting genuine proximity effects

---

### Test C: Nonessential "De-essentialization" Test

**Hypothesis:** Nonessential genes should show more positive median shift (dependency_corrected − dependency) in TRUE than shuffles, indicating removal of spurious essentiality.

**Result:** ✅ **PASSED**

| Subset | TRUE median shift | shuffle_within | shuffle_across |
|--------|-------------------|----------------|----------------|
| all genes | 0.000186 | -0.000242 | -0.000284 |
| **nonessential_only** | **0.000284** | **-0.000207** | **-0.000329** |
| essential_only | -0.000572 | -0.000649 | -0.000819 |

**Interpretation:**
- Nonessential genes show **positive** shift in TRUE (0.00028) vs **negative** in shuffles
- This suggests the correction is successfully pulling nonessential genes up (reducing spurious essentiality) more in real data
- Essential genes show negative shifts in all conditions (expected, as they're being corrected downward)

---

## GO/NO-GO Scoreboard

| Check | TRUE value | shuffle_within | shuffle_across | Pass? |
|-------|------------|----------------|----------------|-------|
| A: Q4 corr-drop TRUE > shuffles | 0.0422 | 0.0589 | 0.0640 | ❌ **NO** |
| B: p90 prox coef TRUE > shuffles | 0.0548 | 0.0405 | 0.0440 | ✅ **YES** |
| C: Nonessential shift TRUE > shuffles | 0.00028 | -0.00021 | -0.00033 | ✅ **YES** |

**Overall: 2/3 PASSED**

---

## Critical Findings

### 🚨 Major Concern: Test A Failure

The failure of Test A is particularly concerning because:
1. **Q4 genes are where we expect the biggest impact** — these are genes with the strongest baseline CNV-dependency coupling
2. **Shuffles outperform TRUE** — this suggests the correction may be removing signal rather than artifact
3. **Pattern is consistent across quartiles** — TRUE consistently shows smaller drops than shuffles at Q2, Q3, and Q4

### Possible Explanations:
- The correction may be **over-correcting** and removing real biological signal
- The proximity model may not be capturing the true structure of the artifact
- Baseline correlations in TRUE may include more real signal than assumed
- The shuffles may be creating artifacts that are easier to "correct" than real artifacts

### ✅ Positive Signals:
- **Test B** confirms real positional signal exists (proximity coefficients are meaningful)
- **Test C** shows the correction is working directionally for nonessential genes

---

## Recommendations

### If proceeding (2/3 pass threshold):
1. **Investigate Test A failure deeply** — why does TRUE underperform shuffles at Q4?
2. **Consider alternative correction approaches** — current method may be too aggressive
3. **Focus manuscript on Tests B & C** — proximity signal and nonessential de-essentialization
4. **Acknowledge limitations** — correlation drop results are mixed

### If pausing/re-scoping:
- The Q4 failure suggests fundamental issues with the correction approach
- May need to revisit the underlying model assumptions
- Consider whether the artifact structure is more complex than proximity alone

---

## Files Generated

- `corr_drop_stratified_by_baseline.csv` — Detailed quartile breakdown
- `proximity_coef_stats.csv` — Coefficient distribution statistics
- `gene_median_shifts.csv` — Gene-level median shifts by subset
- `go_nogo_scoreboard.csv` — Binary pass/fail results

