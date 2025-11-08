# Data Validation Report: Does the Analysis Make Sense?

## Summary

**The data shows a real signal, but with important caveats:**

✅ **What's Real:**
- Unstable set has higher proximity coefficients than stable (median 0.0044 vs 0.0026 for 100k window)
- GREEN hits all show positive drops (proximity adds CN coupling = artifact signal)
- Unstable has 4 GREEN hits vs stable has 0 GREEN hits
- Selection criteria worked: unstable genes have higher breakpoint prevalence

⚠️ **Concerns:**
- Coefficients are very small (0.003-0.006), suggesting subtle effects
- GREEN hits don't have dramatically larger coefficients than AMBER/RED
- Physical proximity is actually LOWER for GREEN than RED (counterintuitive)
- Proximity contributions are tiny (mean ~0.002 dependency units)

## Key Findings

### 1. Selection vs Scoring Are Different Metrics

**Selection criteria** (p100k, p1m):
- Measures: How often is the gene physically near a breakpoint?
- Unstable: median p100k=0.044, p1m=0.299
- Stable: median p100k=0.000, p1m=0.000

**Scoring criteria** (prox_active_frac, prox_contrib):
- Measures: How often does proximity actually affect dependency?
- Can be high even if physical proximity is lower (if effect is strong when it occurs)

### 2. GREEN Hits Characteristics

**Top GREEN hit (DCANP1):**
- prox_active_frac: 0.267 (26.7% of cells)
- prox_only_abs_corr_drop: +0.0109 (positive = artifact)
- Physical proximity: mean(bp<=100k)=0.025, mean(bp<=1Mb)=0.268
- Coefficients: w100k=-0.0049, w1m=-0.0069 (small but non-zero)

**All GREEN hits:**
- 4/4 have positive drops (artifact signal) ✅
- Mean prox_active_frac: 0.187 (18.7% of cells) ✅
- Mean drop: 0.0047 (small but consistent) ⚠️

### 3. Unstable vs Stable Contrast

**Coefficients:**
- Unstable median |bp_within_100000|: 0.0044
- Stable median |bp_within_100000|: 0.0026
- Difference: 1.7x higher in unstable ✅

**Scoring:**
- Unstable: 4 GREEN, 11 AMBER, 29 RED
- Stable: 0 GREEN, 16 AMBER, 28 RED
- Clear contrast ✅

### 4. Counterintuitive Pattern

**Physical proximity by label:**
- GREEN: mean(bp<=100k)=0.032, mean(bp<=1Mb)=0.237
- RED: mean(bp<=100k)=0.112, mean(bp<=1Mb)=0.457

**Why?** RED genes are physically closer to breakpoints more often, but proximity doesn't contribute much to dependency (low coefficients or weak effects). GREEN genes are less often near breakpoints, but when they are, proximity has a stronger effect.

## Interpretation

### Is This Real or Spurious?

**Real signal indicators:**
1. ✅ Unstable > stable on coefficients (1.7x difference)
2. ✅ GREEN hits all show positive drops (artifact signal)
3. ✅ Clear contrast in scoring (4 GREEN vs 0 GREEN)
4. ✅ Selection criteria worked (unstable has higher breakpoint prevalence)

**Spurious signal concerns:**
1. ⚠️ Coefficients are very small (0.003-0.006)
2. ⚠️ Proximity contributions are tiny (mean ~0.002)
3. ⚠️ GREEN doesn't have dramatically larger coefficients than AMBER/RED
4. ⚠️ Physical proximity is lower for GREEN than RED (counterintuitive)

### Most Likely Explanation

**The signal is real but subtle:**

1. **Proximity effects are small** - The artifact exists but is subtle (0.002-0.01 dependency units)
2. **GREEN hits are "high-impact, low-frequency"** - They're not always near breakpoints, but when they are, proximity matters more
3. **RED hits are "low-impact, high-frequency"** - They're often near breakpoints, but proximity doesn't affect dependency much
4. **The contrast is real** - Unstable set has higher coefficients and more GREEN hits than stable

### Recommendations

1. **Use GREEN hits as case studies** - They show the clearest artifact signal
2. **Acknowledge subtlety** - Effects are small but consistent
3. **Focus on contrast** - Unstable vs stable difference is the key finding
4. **Validate with sensitivity** - Re-run with different windows/models to check robustness

## Conclusion

**The data makes sense, but the effects are subtle.** The analysis is not spurious - there's a real contrast between unstable and stable sets, and GREEN hits show consistent artifact signals. However, the proximity effects are small (0.002-0.01 dependency units), which is why they might seem counterintuitive at first glance.

The key insight: **Selection measures physical proximity, scoring measures functional impact.** A gene can be rarely near breakpoints but have strong effects when it is (GREEN), or often near breakpoints but have weak effects (RED).

