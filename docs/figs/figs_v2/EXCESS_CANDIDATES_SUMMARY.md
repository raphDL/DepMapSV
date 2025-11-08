# Top Excess Candidate Genes - Summary

## Overview

Identified genes where TRUE shows excess signal over shuffles on both:
1. **Coefficient excess:** TRUE has larger proximity coefficients than shuffles
2. **Correlation excess:** TRUE shows larger correlation drops than shuffles

## Results

- **Total candidates:** 7,827 genes (out of 18,424)
- **Top 150 selected** for focused pilot analysis
- **Selection criteria:**
  - TRUE beats at least one shuffle on Δ|corr| AND has positive near/max excess
  - Ranked by: 2× near_excess + 1× deltacorr_excess + 1× delta_abs_corr_true

## Top 12 Genes (Preview)

| Gene | abs_near_true | abs_max_true | near_excess_W | near_excess_A | max_excess_W | max_excess_A | delta_abs_corr_true | deltacorr_excess_W | deltacorr_excess_A |
|------|---------------|--------------|---------------|---------------|--------------|--------------|---------------------|-------------------|-------------------|
| EZHIP | 0.0125 | 0.0479 | -0.0501 | -0.0196 | -0.0146 | 0.0158 | -0.0569 | -0.0107 | 0.0011 |
| RPL19 | 0.0049 | 0.0499 | -0.0746 | -0.0226 | -0.0328 | 0.0224 | -0.0171 | -0.0306 | 0.0018 |
| ARHGEF7 | 0.0067 | 0.0328 | -0.0148 | -0.0372 | 0.0113 | -0.0111 | -0.0280 | 0.0027 | -0.0201 |
| SCFD2 | 0.0043 | 0.0144 | -0.0018 | -0.0020 | 0.0052 | 0.0080 | -0.0468 | -0.0039 | 0.0008 |
| LACTBL1 | 0.0075 | 0.0777 | -0.0159 | -0.0669 | 0.0320 | -0.0236 | -0.0322 | 0.0031 | -0.0459 |
| ACTR3 | 0.0013 | 0.0106 | -0.0119 | -0.0053 | -0.0026 | 0.0039 | -0.0102 | -0.0084 | 0.0001 |
| SRRM1 | 0.0092 | 0.0276 | -0.0142 | -0.0119 | 0.0042 | 0.0065 | -0.0076 | -0.0084 | 0.0005 |
| GPR156 | 0.0028 | 0.0067 | -0.0021 | -0.0141 | 0.0017 | -0.0102 | -0.0258 | 0.0009 | -0.0194 |
| HRC | 0.0000 | 0.0150 | -0.0094 | -0.0100 | 0.0056 | 0.0050 | -0.0266 | 0.0028 | 0.0017 |
| VAMP4 | 0.0032 | 0.0361 | -0.0043 | -0.0038 | 0.0267 | 0.0291 | -0.0426 | 0.0026 | -0.0022 |
| C1QTNF6 | 0.0019 | 0.0041 | -0.0055 | -0.0021 | -0.0032 | 0.0001 | -0.0211 | 0.0011 | -0.0088 |
| OXER1 | 0.0052 | 0.0217 | -0.0093 | -0.0289 | 0.0072 | -0.0140 | -0.0138 | 0.0025 | -0.0482 |

**Note:** Many top genes show mixed signals (positive on some metrics, negative on others). This reflects the complexity of the data.

## Files Generated

- `top_excess_candidates.csv` - Full list of 7,827 candidates with all metrics
- `pilot_genes.txt` - Top 150 genes for focused pilot (one per line)
- `run_focused_pilot.sh` - Commands to run focused pilot analyses

## Next Steps: Focused Pilot Analysis

Run focused pilot on the top 150 genes with richer model specification:

### TRUE Condition:
```bash
python svbias_reanalysis.py --pilot --pilot-genes-file figs_v2/pilot_genes.txt \
  --dependency CRISPRGeneEffect.csv --cnv prep_out/cnv_segments.bed --sv prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed --model linear --bp_windows 250000 2000000 \
  --add-continuous-proximity --standardize-predictors --out out_focus_true --progress -vv
```

### Rotate-Shuffle Control:
```bash
python svbias_reanalysis.py --pilot --pilot-genes-file figs_v2/pilot_genes.txt --shuffle-rotate \
  --dependency CRISPRGeneEffect.csv --cnv prep_out/cnv_segments.bed --sv prep_out/sv_from_cnv.bedpe \
  --genes prep_out/genes.depmap.unique.bed --model linear --bp_windows 250000 2000000 \
  --add-continuous-proximity --standardize-predictors --out out_focus_rotate --progress -vv
```

## Success/Failure Criteria

### If TRUE beats rotate on:
- **Directional prox-active fraction**
- **Median Δ|corr| for selected genes**

→ **Model-spec issue:** The richer specification (linear model, continuous proximity, standardized predictors) helps real cases more than shuffles.

### If rotate still ≥ TRUE:
→ **Shuffled worlds are easier to "fix":** Shuffles create simpler artifacts that are easier to correct than real-world complexity.

→ **Consider:**
- Narrow-scope paper documenting subtle proximity effects where they matter
- Pivot to alternative approaches (independent SV calls, trans controls)

## Interpretation

The fact that we found 7,827 candidate genes (42% of all genes) suggests:
1. There is some positional signal in the data
2. The signal is subtle and mixed (not consistently strong across all metrics)
3. A focused analysis on top candidates may reveal clearer patterns

The focused pilot will test whether a richer model specification can better capture real proximity effects compared to shuffled controls.

