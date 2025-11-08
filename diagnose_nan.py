#!/usr/bin/env python3
"""Diagnose NaN issue in comparison."""
import pandas as pd
import numpy as np

# Load coefficients
true_coef = pd.read_csv("out_v3/models_wgs_prox100k/model_coefficients.csv")
rotate_coef = pd.read_csv("out_v3/models_wgs_rotate_prox100k/model_coefficients.csv")

print("=== TRUE Coefficients Diagnostics ===")
print(f"Total genes: {len(true_coef)}")
print(f"NaN in coef_prox: {true_coef['coef_prox_exp_100k'].isna().sum()}")
print(f"Inf in coef_prox: {np.isinf(true_coef['coef_prox_exp_100k']).sum()}")
print(f"\nCoefficient range:")
print(f"  Min: {true_coef['coef_prox_exp_100k'].min()}")
print(f"  Max: {true_coef['coef_prox_exp_100k'].max()}")
print(f"  Mean: {true_coef['coef_prox_exp_100k'].mean()}")
print(f"  Std: {true_coef['coef_prox_exp_100k'].std()}")

# Find extreme outliers
extreme = true_coef[true_coef['coef_prox_exp_100k'].abs() > 10]
print(f"\nExtreme coefficients (|coef| > 10): {len(extreme)}")
if len(extreme) > 0:
    print(extreme[['gene', 'coef_prox_exp_100k', 'coef_cn', 'r2']].head(20))

print("\n=== ROTATE Coefficients Diagnostics ===")
print(f"Total genes: {len(rotate_coef)}")
print(f"NaN in coef_prox: {rotate_coef['coef_prox_exp_100k'].isna().sum()}")
print(f"Inf in coef_prox: {np.isinf(rotate_coef['coef_prox_exp_100k']).sum()}")
print(f"\nCoefficient range:")
print(f"  Min: {rotate_coef['coef_prox_exp_100k'].min()}")
print(f"  Max: {rotate_coef['coef_prox_exp_100k'].max()}")

# Check CN coefficients too
print("\n=== CN Coefficients (TRUE) ===")
print(f"NaN in coef_cn: {true_coef['coef_cn'].isna().sum()}")
print(f"Inf in coef_cn: {np.isinf(true_coef['coef_cn']).sum()}")
print(f"Range: {true_coef['coef_cn'].min()} to {true_coef['coef_cn'].max()}")

# Check R² values
print("\n=== R² Values ===")
print(f"TRUE R² - NaN: {true_coef['r2'].isna().sum()}")
print(f"TRUE R² - negative: {(true_coef['r2'] < 0).sum()}")
print(f"TRUE R² range: {true_coef['r2'].min()} to {true_coef['r2'].max()}")

