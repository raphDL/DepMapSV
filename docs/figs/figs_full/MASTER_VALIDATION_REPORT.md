# Master Validation Report: SV Bias Correction
**Date:** Generated from `out_v2/` outputs  
**Purpose:** Comprehensive validation of SV bias correction approach

---

## Executive Summary

**Overall Score: 3/7 tests PASSED (43%)**

This report combines:
- **CPR Tests** (3 tests): Initial validation on correlation drops, proximity coefficients, and gene shifts
- **Falsification Tests** (4 tests): Additional validation on held-out data, CN-neutral effects, effect sizes, and unstable genes

### Quick Scorecard

| Test | Result | Notes |
|------|--------|-------|
| **CPR-A:** Q4 corr-drop | ❌ FAILED | TRUE < shuffles |
| **CPR-B:** Prox coef p90 | ✅ PASSED | TRUE > shuffles |
| **CPR-C:** Nonessential shift | ✅ PASSED | TRUE > shuffles |
| **FALS-1:** Held-out | ❌ FAILED | TRUE < shuffles |
| **FALS-2:** CN-neutral | ✅ PASSED | TRUE > shuffles (barely) |
| **FALS-3:** Effect-size | ❌ FAILED | TRUE < shuffles |
| **FALS-4:** Unstable genes | ❌ FAILED | TRUE < shuffles |

---

## Critical Finding: Consistent Failure Pattern

### 🚨 The Core Problem

**TRUE consistently underperforms shuffles on correlation drop metrics:**

| Metric | TRUE | shuffle_within | shuffle_across | TRUE wins? |
|--------|------|----------------|----------------|------------|
| Q4 corr-drop (CPR) | 0.042 | 0.059 | 0.064 | ❌ |
| Held-out drop | 0.008 | 0.012 | 0.013 | ❌ |
| Effect-size drop | 0.141 | 0.208 | 0.236 | ❌ |
| Unstable dir_frac | 0.549 | 0.596 | 0.599 | ❌ |

**This pattern appears across:**
- Different gene subsets (Q4, unstable, all genes)
- Different validation approaches (full data, held-out, CN-neutral)
- Different normalization schemes (raw, scaled)

### Possible Explanations

1. **Over-correction:** The model is removing real biological signal, not just artifacts
2. **Wrong artifact model:** The proximity-based model doesn't match the true artifact structure
3. **Shuffle advantage:** Shuffles create simpler artifacts that are easier to "correct"
4. **Signal vs noise:** Real CNV-dependency relationships are being removed as "artifacts"

---

## Positive Signals (Limited)

### ✅ What Works

1. **Proximity Coefficients (CPR-B)**
   - TRUE shows larger p90/p99 proximity coefficients than shuffles
   - Suggests real positional signal is being detected
   - **But:** This doesn't translate to better correction performance

2. **Nonessential De-essentialization (CPR-C)**
   - Nonessential genes show positive shift in TRUE vs negative in shuffles
   - Suggests directional correction is working for nonessentials
   - **But:** Effect is very small (0.00028 vs -0.00021)

3. **CN-Neutral Test (FALS-2)**
   - TRUE shows slightly larger drop than shuffle_across at CN ~2
   - **But:** Essentially tied with shuffle_within, very small effect, limited power (540 genes)

### ⚠️ Weak Signals

- All positive signals are **small in magnitude**
- Most are **barely distinguishable** from shuffles
- None show **consistent superiority** across multiple tests

---

## Detailed Test Results

### CPR Tests

#### Test A: Q4 Correlation Drop
- **Result:** ❌ FAILED
- **Finding:** TRUE (0.042) < shuffle_within (0.059) < shuffle_across (0.064)
- **Interpretation:** Correction underperforms where we expect biggest impact

#### Test B: Proximity Coefficients
- **Result:** ✅ PASSED
- **Finding:** TRUE p90 (0.055) > shuffle_within (0.040) > shuffle_across (0.044)
- **Interpretation:** Real positional signal detected, but doesn't improve correction

#### Test C: Nonessential Shift
- **Result:** ✅ PASSED
- **Finding:** TRUE (+0.00028) > shuffle_within (-0.00021) > shuffle_across (-0.00033)
- **Interpretation:** Directional correction works for nonessentials, but effect is tiny

### Falsification Tests

#### Test 1: Held-Out Falsification
- **Result:** ❌ FAILED
- **Finding:** TRUE (0.008) < shuffle_within (0.012) < shuffle_across (0.013)
- **Interpretation:** Correction doesn't generalize better than random on held-out data

#### Test 2: CN-Neutral Falsification
- **Result:** ✅ PASSED (barely)
- **Finding:** TRUE (0.00129) ≈ shuffle_within (0.00129) > shuffle_across (0.00041)
- **Interpretation:** Marginal pass, essentially tied with shuffle_within

#### Test 3: Effect-Size Normalization
- **Result:** ❌ FAILED
- **Finding:** TRUE (0.141) < shuffle_within (0.208) < shuffle_across (0.236)
- **Interpretation:** Normalized effects are smaller in TRUE, suggesting less effective correction

#### Test 4: Unstable Genes
- **Result:** ❌ FAILED
- **Finding:** TRUE (0.549) < shuffle_within (0.596) < shuffle_across (0.599)
- **Interpretation:** Correction underperforms on proximity-enriched genes where we expect biggest effect

---

## Statistical Power Considerations

### Tests with Limited Power
- **CN-neutral test:** Only 540 genes have sufficient data at CN ~2
- **Unstable genes:** 1,843 genes (top decile), but still shows failure
- **Held-out test:** 17,922 genes, good power, still fails

### Tests with Good Power
- **Q4 correlation drop:** 4,606 genes per quartile
- **Effect-size normalization:** 17,922 genes
- **Proximity coefficients:** 18,435 genes

**Conclusion:** Power is not the issue — the failures are consistent across tests with varying power.

---

## Recommendations

### 🛑 Primary Recommendation: Pause and Re-scope

The evidence is **strong and consistent** that the current correction approach has fundamental issues:

1. **Multiple independent tests fail** — not a single fluke
2. **Consistent pattern** — TRUE underperforms shuffles across metrics
3. **Fails where it should work best** — Q4, unstable genes, held-out data
4. **Limited positive signals** — only 3/7 tests pass, and those are weak

### Alternative Actions

#### Option 1: Deep Investigation (If Continuing)
- Investigate why shuffles outperform TRUE
- Analyze what the shuffles are actually "correcting"
- Compare correction patterns between TRUE and shuffles
- Consider if shuffles are appropriate null models

#### Option 2: Method Revision
- Revisit the proximity-based artifact model
- Consider alternative artifact structures
- Test different correction approaches
- Validate on known artifacts if available

#### Option 3: Scope Reduction
- Frame as "proof of concept" rather than validated method
- Focus on proximity signal detection (Test B) rather than correction
- Acknowledge limitations prominently
- Publish as exploratory analysis

### If Proceeding Despite Results

**Required disclosures:**
1. Most falsification tests fail
2. TRUE consistently underperforms shuffles on key metrics
3. Effects are small and inconsistent
4. Correction may be removing real signal
5. Results should be interpreted with caution

---

## Files and Data

### Generated Files
- `corr_drop_stratified_by_baseline.csv` — CPR Test A results
- `proximity_coef_stats.csv` — CPR Test B results
- `gene_median_shifts.csv` — CPR Test C results
- `go_nogo_scoreboard.csv` — CPR binary results
- `falsification_tests_summary.csv` — Falsification detailed results
- `falsification_scoreboard.csv` — Falsification binary results
- `unstable_genes.txt` — List of unstable genes
- `CPR_TEST_SUMMARY.md` — CPR test interpretation
- `FALSIFICATION_TESTS_SUMMARY.md` — Falsification test interpretation
- `MASTER_VALIDATION_REPORT.md` — This document

### Source Data
- `out_v2/true/` — True SV data results
- `out_v2/shuffle_within_chrom/` — Within-chromosome shuffle results
- `out_v2/shuffle_across_chrom/` — Across-chromosome shuffle results

---

## Conclusion

The comprehensive validation suite provides **strong evidence against** the current SV bias correction approach. The consistent pattern of TRUE underperforming shuffles across multiple independent tests suggests fundamental issues that should be addressed before proceeding.

**The data speaks clearly: the correction is not working as intended.**

---

## Next Steps

1. **Review this report** with collaborators
2. **Decide on path forward:**
   - Pause and re-scope (recommended)
   - Deep investigation of failures
   - Method revision
   - Scope reduction
3. **If continuing:** Address the consistent failure pattern
4. **If pausing:** Document learnings for future work

---

*Generated by automated validation pipeline*

