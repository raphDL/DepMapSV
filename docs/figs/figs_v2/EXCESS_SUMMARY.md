# Excess tests (full out_v2)

## Coefficient-excess (TRUE minus shuffle)
       test   metric               versus  n_genes  median_excess    ci_low   ci_high
coef_excess abs_near shuffle_within_chrom    18435       0.003381  0.003137  0.003608
coef_excess abs_near shuffle_across_chrom    18435       0.002118  0.001870  0.002405
coef_excess  abs_far shuffle_within_chrom    18435      -0.000412 -0.000588 -0.000244
coef_excess  abs_far shuffle_across_chrom    18435      -0.001613 -0.001792 -0.001431
coef_excess  abs_max shuffle_within_chrom    18435       0.003973  0.003744  0.004204
coef_excess  abs_max shuffle_across_chrom    18435       0.002732  0.002466  0.002991

## Δ|corr|-excess (TRUE minus shuffle)
             test         metric               versus  n_genes  median_excess    ci_low   ci_high
delta_corr_excess delta_abs_corr shuffle_within_chrom    18424      -0.006125 -0.006450 -0.005826
delta_corr_excess delta_abs_corr shuffle_across_chrom    18424      -0.008158 -0.008593 -0.007845

## Run medians
                 run  med_abs_near  med_abs_far  med_abs_max  med_delta_abs_corr
                true      0.014849     0.009368     0.018636            0.012463
shuffle_within_chrom      0.010868     0.010131     0.014232            0.018236
shuffle_across_chrom      0.012525     0.011496     0.015811            0.020317

## Decision rubric
- Coefficient-excess CI>0 in 4/6 contrasts.
- Δ|corr|-excess CI>0 in 0/2 contrasts.

**Interpretation:** Some position-specific signal exceeds shuffles → keep exploring (with restraint).
