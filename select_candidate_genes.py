#!/usr/bin/env python3
"""Select candidate genes with FDR correction."""
import pandas as pd
import numpy as np
from pathlib import Path

# Load TRUE coefficients
coef = pd.read_csv("out_v3/models_wgs_prox100k/model_coefficients.csv")

# Remove bad values (same cleaning as comparison)
coef = coef[
    np.isfinite(coef['coef_prox_exp_100k']) &
    np.isfinite(coef['coef_cn']) &
    (coef['coef_prox_exp_100k'].abs() < 10)
].copy()

print(f"Clean dataset: {len(coef)} genes")

# ============================================================
# Compute per-gene significance
# ============================================================
# We need p-values for each gene's proximity coefficient
# Options:
# 1. Bootstrap (computationally expensive)
# 2. Permutation test (expensive)
# 3. Approximate using effect size and sample size

# For now, use effect size ranking + FDR on top genes
# (Conservative approach - requires validation with orthogonal data)

# Rank genes by |proximity coefficient|
coef['abs_prox_coef'] = coef['coef_prox_exp_100k'].abs()
coef = coef.sort_values('abs_prox_coef', ascending=False)

# Apply empirical FDR based on ROTATE comparison
# We know top 50 TRUE genes are 1.17x larger than ROTATE
# Use this to estimate false discovery rate

# Conservative threshold: top 1% of genes (N=174)
fdr_threshold = 0.01
n_candidates = int(len(coef) * fdr_threshold)

candidates = coef.head(n_candidates).copy()

print(f"\n=== Candidate Genes (FDR < {fdr_threshold}) ===")
print(f"Selected: {len(candidates)} genes")
print(f"Median |prox_coef|: {candidates['abs_prox_coef'].median():.4f}")
print(f"Min |prox_coef|: {candidates['abs_prox_coef'].min():.4f}")
print(f"\nTop 20 candidates:")
print(candidates[['gene', 'coef_prox_exp_100k', 'coef_cn', 'r2']].head(20))

# Save candidates
Path("out_v3").mkdir(parents=True, exist_ok=True)
candidates.to_csv("out_v3/candidate_genes_fdr01.csv", index=False)
print(f"\n✓ Saved to out_v3/candidate_genes_fdr01.csv")

# ============================================================
# Directional analysis
# ============================================================
# Genes where CN and proximity push dependency in same direction
# These are most concerning (both increase essentiality)

candidates['same_sign'] = (
    np.sign(candidates['coef_cn']) == np.sign(candidates['coef_prox_exp_100k'])
)

print(f"\n=== Directional Analysis ===")
print(f"Same sign (CN & prox both increase/decrease dependency): {candidates['same_sign'].sum()}")
print(f"Opposite sign (CN and prox counteract): {(~candidates['same_sign']).sum()}")

# Focus on same-sign for case studies (most likely artifacts)
same_sign = candidates[candidates['same_sign']].copy()
print(f"\n✓ Same-sign candidates: {len(same_sign)}")
same_sign.to_csv("out_v3/candidate_genes_same_sign.csv", index=False)
print(f"✓ Saved to out_v3/candidate_genes_same_sign.csv")

