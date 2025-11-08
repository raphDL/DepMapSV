# Biology-focused proximity tests (focused pilot)

Sets: TRUE=out_focus_true/unstable  ROTATE=out_focus_rotate/unstable

Params: W1/W2=250,000–2,000,000 bp, contrib_thr=0.01, active_frac=0.1

## Prevalence-matched directional (decile-weighted)
TRUE=0.000  ROTATE=0.046  EXCESS=-0.046  CI[-0.092,+0.000]

Per-bin summary (equal-weight across bins):

Saved: pivot_bio_bins.csv

## Coefficient → effect consistency (Spearman)
TRUE:   rho=0.171, p=4.36e-02, n=140
ROTATE: rho=-0.070, p=4.09e-01, n=140

Top genes by |coef| (TRUE):

   gene  active_frac  corr_drop  n_cells  abs_coef_max
  PPWD1     0.048193   0.004237      830      0.018228
TUBGCP4     0.057831   0.012337      830      0.014059
    MTR     0.055422  -0.000455      830      0.011585
 SOWAHB     0.037349   0.002839      830      0.011352
    MT4     0.049398   0.004799      830      0.011095
  CSE1L     0.062651   0.002449      830      0.010812
  VAMP4     0.053012   0.001074      830      0.010700
  NEDD1     0.063855   0.005325      830      0.010299
DENND4A     0.000000   0.017426      830      0.009546
  KCNH4     0.000000   0.014791      830      0.009245
  APOL5     0.000000   0.000779      830      0.008894
   LY6D     0.000000   0.011191      830      0.008807
  RAB18     0.000000   0.007655      830      0.008748
  ITGB5     0.000000  -0.000249      830      0.008489
   RTF2     0.000000   0.004337      830      0.007922
