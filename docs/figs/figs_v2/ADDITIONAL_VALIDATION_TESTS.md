# Additional Validation Tests on Focused Pilot

## Overview

Three additional validation tests on the focused pilot (top 150 excess candidate genes) to assess whether TRUE shows genuine positional signal beyond prevalence artifacts.

## Test 1: Prevalence-Matched Directional Metric

**Goal:** Fix the "more breakpoints → more chances" bias by stratifying genes by proximity prevalence and comparing within each stratum.

**Method:**
- Stratify genes by p1Mb deciles (based on TRUE prevalence)
- Compare TRUE vs ROTATE within each decile
- Equal-weight average across deciles
- Bootstrap CI (10,000 reps)

**Results:**
```
Prevalence-matched directional:
TRUE=0.018  ROTATE=0.000  EXCESS=+0.018  CI[+0.000,+0.054]
```

**Interpretation:**
- ✅ TRUE shows higher directional fraction (0.018 vs 0.000)
- ⚠️ **CI includes 0** - Excess is positive but not statistically significant
- The CI lower bound is exactly 0.000, suggesting marginal evidence

**Conclusion:** Weak positive signal, but not statistically significant.

---

## Test 2: Coefficient→Effect Consistency

**Goal:** Check if larger proximity coefficients correlate with larger correlation drops (Spearman correlation).

**Method:**
- Compute max absolute proximity coefficient per gene
- Compute correlation drop (Δ|corr|) per gene
- Spearman correlation between coefficient magnitude and correlation drop

**Results:**
```
Spearman(abs_coef_max, corr_drop)
TRUE:   rho=0.019  p=8.22e-01  n=140
ROTATE: rho=-0.029  p=7.29e-01  n=140
```

**Interpretation:**
- Both correlations are **essentially zero** (not significant)
- TRUE: rho=0.019 (p=0.82) - no correlation
- ROTATE: rho=-0.029 (p=0.73) - no correlation
- Neither condition shows that larger coefficients → larger drops

**Conclusion:** No evidence that proximity coefficients predict correlation improvement in either condition.

---

## Test 3: Per-Gene Permutation P-Values

**Goal:** Within each gene, shuffle proximity flags and get p-value for directional effect (FDR across genes).

**Method:**
- For each gene, shuffle proximity assignment (preserve counts)
- Compute correlation drop under shuffle (400 permutations)
- P-value: fraction of shuffles with drop ≥ observed
- FDR correction (q<0.1 threshold)

**Results:**
```
Directional prox-only (perm p): TRUE q<0.1 frac=0.329  n=143
Directional prox-only (perm p): ROTATE q<0.1 frac=0.203  n=143
```

**Interpretation:**
- ✅ **TRUE wins:** 32.9% of genes significant (q<0.1) vs 20.3% for ROTATE
- This is the **strongest positive signal** for TRUE
- Suggests proximity effects are more consistent in TRUE than ROTATE

**Conclusion:** TRUE shows more genes with significant proximity-only correction effects.

---

## Combined Assessment

| Test | TRUE | ROTATE | Winner | Significance |
|------|------|--------|--------|--------------|
| **1. Prevalence-matched directional** | 0.018 | 0.000 | ✅ TRUE | ⚠️ CI includes 0 |
| **2. Coef→effect consistency** | rho=0.019 | rho=-0.029 | Tie | ❌ Both non-significant |
| **3. Permutation p-values (q<0.1)** | 32.9% | 20.3% | ✅ TRUE | ✅ TRUE > ROTATE |

**Overall:** 2/3 tests favor TRUE, but signals are weak/marginal.

---

## Key Findings

### Positive Signals for TRUE:
1. **Prevalence-matched directional:** TRUE shows excess (0.018 vs 0.000), but CI includes 0
2. **Permutation p-values:** TRUE has more significant genes (32.9% vs 20.3%)

### Negative/Neutral Signals:
1. **Coefficient→effect consistency:** No correlation in either condition
2. **Prevalence-matched CI:** Includes 0, not statistically significant

### Why Results Are Mixed:

1. **Prevalence artifact:** ROTATE had higher proximity prevalence (p100k=0.0628 vs 0.0351), which can inflate threshold-based metrics
2. **CN-derived breakpoints:** If breakpoints are derived from CN segmentation, they're entangled with CN, making it hard to isolate proximity effects
3. **Real-world complexity:** TRUE data has more complex structure than shuffles, making correction harder

---

## Recommendations

### If Proceeding (Weak Evidence):
- Focus on **Test 3 (permutation p-values)** as the strongest signal
- Acknowledge that prevalence-matched directional metric is marginal (CI includes 0)
- Note that coefficient→effect consistency is absent in both conditions
- Frame as "subtle proximity effects detected in subset of genes"

### If Pausing/Pivoting:
- The mixed results (2/3 favor TRUE, but all are weak) suggest:
  - Positional signal exists but is subtle
  - Current correction approach doesn't reliably reduce CN-coupling
  - Shuffles create simpler artifacts that are easier to "fix"
- Consider:
  - Narrow-scope paper documenting where proximity effects matter
  - Alternative approaches (independent SV calls, trans controls)
  - Focus on understanding why shuffles outperform TRUE

---

## Conclusion

The additional validation tests provide **weak, mixed evidence** for TRUE > ROTATE:
- Test 1: Marginal (CI includes 0)
- Test 2: No difference (both show no correlation)
- Test 3: TRUE wins (32.9% vs 20.3% significant genes)

**Overall assessment:** There is some positional signal, but it's subtle and the correction approach has limitations. The evidence is not strong enough to confidently claim the correction works better on real data than shuffles.

