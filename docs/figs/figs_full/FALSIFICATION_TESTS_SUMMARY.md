# Falsification Tests Summary
**Date:** Generated from `out_v2/` outputs  
**Purpose:** Additional decisive tests to validate correction beyond initial CPR tests

---

## Executive Summary

**Score: 1/4 tests PASSED** ⚠️⚠️

- ❌ **Test 1:** Held-out falsification (TRUE < shuffles)
- ✅ **Test 2:** CN-neutral falsification (TRUE > shuffles) — **BARE PASS**
- ❌ **Test 3:** Effect-size normalization (TRUE < shuffles)
- ❌ **Test 4:** Proximity-enriched set (TRUE < shuffles)
- ℹ️ **Test 5:** Alternative AUROC metrics (all near zero, no clear winner)

**Combined with CPR tests: 3/7 total tests passed**

---

## Detailed Results

### Test 1: Held-Out Falsification (80/20 Split)

**Hypothesis:** Training on 80% of cells, testing on 20% should show TRUE > shuffles if correction generalizes.

**Result:** ❌ **FAILED**

| Condition | n_genes | median_te_drop | frac_positive |
|-----------|---------|----------------|---------------|
| **TRUE** | 17,922 | **0.0081** | 0.630 |
| shuffle_within | 17,922 | **0.0118** | 0.658 |
| shuffle_across | 17,922 | **0.0133** | 0.663 |

**Interpretation:**
- TRUE shows **smaller** held-out correlation drop than both shuffles
- Fraction positive is also lower in TRUE (0.63 vs 0.66-0.66)
- This is a **serious strike** — the correction doesn't generalize better than random on held-out data
- Suggests the model may be overfitting or the correction is removing signal rather than artifact

---

### Test 2: CN-Neutral Falsification (Ploidy ~2)

**Hypothesis:** Isolating proximity from CN amplitude by focusing on cells with CN ~2 should show TRUE > shuffles.

**Result:** ✅ **PASSED (barely)**

| Condition | n_genes | median_cn2_drop | frac_positive |
|-----------|---------|------------------|---------------|
| **TRUE** | 540 | **0.00129** | 0.539 |
| shuffle_within | 540 | **0.00129** | 0.546 |
| shuffle_across | 540 | **0.00041** | 0.513 |

**Interpretation:**
- TRUE shows **slightly larger** drop than shuffle_within (0.00129 vs 0.00129 — essentially tied)
- TRUE shows **larger** drop than shuffle_across (0.00129 vs 0.00041)
- This is a **marginal pass** — the effect is very small and barely distinguishable from shuffle_within
- Only 540 genes have sufficient data at CN ~2, so power is limited
- Fraction positive is actually slightly lower in TRUE (0.539 vs 0.546 for within)

---

### Test 3: Effect-Size Normalization

**Hypothesis:** Scaling Δ|r| by dependency variability should show TRUE > shuffles when comparing apples to apples.

**Result:** ❌ **FAILED**

| Condition | n_genes | median_scaled_drop | frac_pos |
|-----------|---------|-------------------|----------|
| **TRUE** | 17,922 | **0.141** | 0.697 |
| shuffle_within | 17,922 | **0.208** | 0.752 |
| shuffle_across | 17,922 | **0.236** | 0.763 |

**Interpretation:**
- TRUE shows **much smaller** scaled drop than both shuffles
- This suggests that when normalized by variability, the correction is **less effective** in TRUE
- The shuffles are showing larger normalized effects, which is concerning
- Fraction positive is also lower in TRUE (0.70 vs 0.75-0.76)

---

### Test 4: Proximity-Enriched Set (Unstable Genes)

**Hypothesis:** Focusing on unstable genes (top decile by p<=1Mb) should show TRUE > shuffles where proximity actually occurs.

**Result:** ❌ **FAILED**

| Condition | dir_frac |
|-----------|----------|
| **TRUE** | **0.549** |
| shuffle_within | **0.596** |
| shuffle_across | **0.599** |

**Interpretation:**
- TRUE shows **lower** directional fraction than both shuffles
- This is particularly concerning because unstable genes are where we expect the biggest effect
- The shuffles are showing **more** directional correction on proximity-enriched genes
- This suggests the correction may not be capturing the true structure of proximity artifacts

---

### Test 5: Alternative AUROC Metrics

**Hypothesis:** Alternative metrics (AUPRC on bottom 50%, simplified AUROC) should show improvement.

**Result:** ℹ️ **INCONCLUSIVE**

| Condition | n_genes | auprc_bottom50_delta | auroc_delta |
|-----------|---------|---------------------|-------------|
| **TRUE** | 2,254 | 4.83e-06 | 2.07e-05 |
| shuffle_within | 2,254 | 4.12e-06 | 1.89e-05 |
| shuffle_across | 2,254 | 9.46e-06 | 3.69e-05 |

**Interpretation:**
- All values are **extremely small** (near zero)
- TRUE is slightly better than shuffle_within on both metrics
- shuffle_across actually shows the largest improvement
- These differences are so small they're likely not meaningful
- The original AUROC saturation issue persists even with alternative metrics

---

## Combined Analysis: CPR + Falsification Tests

### Overall Scorecard

| Test Category | Test | Result |
|---------------|------|--------|
| **CPR Tests** | A: Q4 corr-drop | ❌ FAILED |
| | B: Prox coef p90 | ✅ PASSED |
| | C: Nonessential shift | ✅ PASSED |
| **Falsification** | 1: Held-out | ❌ FAILED |
| | 2: CN-neutral | ✅ PASSED (barely) |
| | 3: Effect-size | ❌ FAILED |
| | 4: Unstable genes | ❌ FAILED |
| | 5: AUROC alt | ℹ️ INCONCLUSIVE |

**Total: 3/7 tests passed (43%)**

---

## Critical Patterns

### 🚨 Consistent Failure Pattern

**TRUE consistently underperforms shuffles on:**
1. Correlation drop metrics (Q4, held-out, effect-size)
2. Directional fraction on unstable genes
3. Fraction positive metrics

**This suggests:**
- The correction may be **removing real signal** rather than artifact
- The proximity model may not capture the true artifact structure
- Shuffles may be creating artifacts that are easier to "correct" than real artifacts
- The correction may be **over-aggressive** and removing legitimate CNV-dependency relationships

### ✅ Limited Positive Signals

**TRUE outperforms shuffles on:**
1. Proximity coefficient magnitudes (p90, p99) — suggests real positional signal exists
2. Nonessential gene median shift — suggests directional correction for nonessentials
3. CN-neutral test (barely) — very small effect, essentially tied with shuffle_within

**These are weak signals compared to the consistent failures.**

---

## Recommendations

### 🛑 Strong Signal to Pause/Re-scope

The consistent pattern of TRUE underperforming shuffles across multiple independent tests is a **strong signal** that:

1. **The correction approach may be fundamentally flawed**
   - It's not generalizing better than random
   - It's not working better on proximity-enriched genes
   - It's not showing larger effects when normalized appropriately

2. **The model assumptions may be incorrect**
   - The proximity-based artifact model may not match reality
   - The correction may be removing real biological signal
   - The shuffles may not be appropriate null models

3. **Alternative approaches should be considered**
   - Different artifact models
   - Different correction strategies
   - Different validation approaches

### If Proceeding (Not Recommended)

If you choose to proceed despite these results:

1. **Acknowledge limitations prominently**
   - Most falsification tests fail
   - TRUE underperforms shuffles on key metrics
   - Effects are small and inconsistent

2. **Focus on what works**
   - Proximity coefficients show real signal
   - Nonessential de-essentialization shows directional effect
   - Frame as "proof of concept" rather than validated method

3. **Investigate the failures**
   - Why do shuffles outperform TRUE?
   - What is the true structure of the artifact?
   - Are we over-correcting?

---

## Files Generated

- `falsification_tests_summary.csv` — Detailed results table
- `falsification_tests_detailed.json` — Full results in JSON format
- `falsification_scoreboard.csv` — Binary pass/fail results
- `unstable_genes.txt` — List of unstable genes (top decile by p<=1Mb)

---

## Conclusion

The falsification tests provide **strong evidence against** the current correction approach. The consistent pattern of TRUE underperforming shuffles across multiple independent tests suggests fundamental issues with either:

1. The correction method itself
2. The underlying model assumptions
3. The validation approach

**Recommendation: Pause and re-scope** before investing more time in this approach.

