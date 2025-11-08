#!/usr/bin/env python3
"""Select case study genes for detailed validation."""
import pandas as pd
import numpy as np
from pathlib import Path

# Load candidates
candidates = pd.read_csv("out_v3/candidate_genes_fdr01.csv")

print(f"Loaded {len(candidates)} candidate genes")

# Selection criteria:
# 1. Same sign with CN (most likely artifacts)
# 2. High R² (model explains variance well)
# 3. Top proximity coefficient (effect size)

# Add same_sign flag
candidates['same_sign'] = (
    np.sign(candidates['coef_cn']) == np.sign(candidates['coef_prox_exp_100k'])
)

# Select case studies (prioritize same-sign, high effect size)
# Use top genes by absolute proximity coefficient among same-sign candidates
same_sign_candidates = candidates[candidates['same_sign']].copy()
same_sign_candidates = same_sign_candidates.sort_values('abs_prox_coef', ascending=False)

# Take top 30 same-sign genes, optionally filter by R²
case_studies = same_sign_candidates[
    (same_sign_candidates['r2'] > 0.01)  # Very relaxed R² threshold
].head(30)

print(f"\n=== Case Study Genes ({len(case_studies)}) ===")
print("Selection criteria:")
print("  - Same sign (CN & proximity)")
print("  - R² > 0.01")
print("  - Top 30 by |proximity coefficient|")
print(f"\nSelected {len(case_studies)} genes:\n")

for idx, (_, row) in enumerate(case_studies.iterrows(), 1):
    print(f"{idx:2d}. {row['gene']:15s} prox={row['coef_prox_exp_100k']:7.4f}  "
          f"cn={row['coef_cn']:7.4f}  r2={row['r2']:.4f}")

# Save case studies
Path("out_v3").mkdir(parents=True, exist_ok=True)
case_studies.to_csv("out_v3/case_study_genes.csv", index=False)
print(f"\n✓ Saved to out_v3/case_study_genes.csv")

# Summary statistics
print(f"\n=== Case Study Summary ===")
print(f"Mean |prox_coef|: {case_studies['abs_prox_coef'].mean():.4f}")
print(f"Mean R²: {case_studies['r2'].mean():.4f}")
print(f"Mean |CN coef|: {case_studies['coef_cn'].abs().mean():.4f}")
