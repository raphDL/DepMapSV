#!/usr/bin/env python3
"""Compare TRUE vs ROTATE shuffle - ROBUST VERSION."""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path

# Load coefficients
true_coef = pd.read_csv("out_v3/models_wgs_prox100k/model_coefficients.csv")
rotate_coef = pd.read_csv("out_v3/models_wgs_rotate_prox100k/model_coefficients.csv")

print(f"Loaded {len(true_coef)} TRUE genes, {len(rotate_coef)} ROTATE genes")

# ============================================================
# CRITICAL: Clean data before comparison
# ============================================================
print("\n=== Data Cleaning ===")

# Remove NaN and Inf values
def clean_coefficients(df, label):
    original = len(df)
    
    # Remove NaN
    df = df.dropna(subset=['coef_prox_exp_100k', 'coef_cn', 'r2'])
    after_nan = len(df)
    
    # Remove Inf
    mask_finite = (
        np.isfinite(df['coef_prox_exp_100k']) & 
        np.isfinite(df['coef_cn']) & 
        np.isfinite(df['r2'])
    )
    df = df[mask_finite].copy()
    after_inf = len(df)
    
    # Remove extreme outliers (likely numerical errors)
    # Keep only |coef| < 10 (dependency is [-1, 1], so coef > 10 is absurd)
    mask_reasonable = (
        (df['coef_prox_exp_100k'].abs() < 10) &
        (df['coef_cn'].abs() < 10)
    )
    df = df[mask_reasonable].copy()
    after_outlier = len(df)
    
    print(f"{label}:")
    print(f"  Original: {original}")
    print(f"  After removing NaN: {after_nan} (dropped {original - after_nan})")
    print(f"  After removing Inf: {after_inf} (dropped {after_nan - after_inf})")
    print(f"  After removing outliers: {after_outlier} (dropped {after_inf - after_outlier})")
    
    return df

true_coef = clean_coefficients(true_coef, "TRUE")
rotate_coef = clean_coefficients(rotate_coef, "ROTATE")

# Merge on gene (inner join to keep only genes present in both)
# Include coef_cn from TRUE only (it's the same for both)
df = true_coef[['gene', 'coef_prox_exp_100k', 'coef_cn', 'r2']].merge(
    rotate_coef[['gene', 'coef_prox_exp_100k', 'r2']],
    on='gene',
    suffixes=('_true', '_rotate'),
    how='inner'  # Only genes in both datasets
)
# Rename coef_cn to coef_cn_true for consistency
df = df.rename(columns={'coef_cn': 'coef_cn_true'})

print(f"\nFinal comparison dataset: {len(df)} genes")

if len(df) < 100:
    print("\n⚠️  WARNING: Very few genes remained after cleaning!")
    print("This suggests serious data quality issues.")
    print("Check your model fitting - Huber regression may have failed.")

# ============================================================
# TEST 1: Global coefficient distribution
# ============================================================
print("\n=== TEST 1: Coefficient Distributions ===")
print(f"TRUE:   mean={df['coef_prox_exp_100k_true'].mean():.6f}, "
      f"std={df['coef_prox_exp_100k_true'].std():.6f}")
print(f"ROTATE: mean={df['coef_prox_exp_100k_rotate'].mean():.6f}, "
      f"std={df['coef_prox_exp_100k_rotate'].std():.6f}")

# Two-sample KS test
try:
    ks_stat, ks_p = stats.ks_2samp(
        df['coef_prox_exp_100k_true'].abs(),
        df['coef_prox_exp_100k_rotate'].abs()
    )
    print(f"\nKS test: D={ks_stat:.4f}, p={ks_p:.4e}")
except Exception as e:
    print(f"\n⚠️  KS test failed: {e}")
    ks_p = np.nan

# Wilcoxon signed-rank test (paired)
try:
    w_stat, w_p = stats.wilcoxon(
        df['coef_prox_exp_100k_true'].abs(),
        df['coef_prox_exp_100k_rotate'].abs(),
        nan_policy='raise'  # Fail if NaN present
    )
    print(f"Wilcoxon: W={w_stat:.0f}, p={w_p:.4e}")
except Exception as e:
    print(f"⚠️  Wilcoxon test failed: {e}")
    w_p = np.nan

# ============================================================
# TEST 2: R² comparison
# ============================================================
print("\n=== TEST 2: Model Fit (R²) ===")
print(f"TRUE:   mean R²={df['r2_true'].mean():.6f}")
print(f"ROTATE: mean R²={df['r2_rotate'].mean():.6f}")

r2_diff = df['r2_true'] - df['r2_rotate']
print(f"Difference: mean={r2_diff.mean():.6f}, median={r2_diff.median():.6f}")

# Test if TRUE R² > ROTATE R² (one-sided)
try:
    w_r2, p_r2 = stats.wilcoxon(r2_diff, alternative='greater', nan_policy='raise')
    print(f"Paired test (TRUE > ROTATE): W={w_r2:.0f}, p={p_r2:.4e}")
except Exception as e:
    print(f"⚠️  R² test failed: {e}")
    p_r2 = np.nan

# ============================================================
# TEST 3: Top genes consistency
# ============================================================
print("\n=== TEST 3: Top Genes ===")

# Define "directional" genes
df['cn_prox_same_sign'] = (
    np.sign(df['coef_cn_true']) == np.sign(df['coef_prox_exp_100k_true'])
)
directional_genes = df[df['cn_prox_same_sign']]['gene'].values
print(f"Directional genes: {len(directional_genes)}")

# Top genes comparison
top_true_mean = df.nlargest(50, 'coef_prox_exp_100k_true')['coef_prox_exp_100k_true'].abs().mean()
top_rotate_mean = df.nlargest(50, 'coef_prox_exp_100k_rotate')['coef_prox_exp_100k_rotate'].abs().mean()
top_ratio = top_true_mean / top_rotate_mean if top_rotate_mean > 0 else np.nan

print(f"Top 50 TRUE mean: {top_true_mean:.6f}")
print(f"Top 50 ROTATE mean: {top_rotate_mean:.6f}")
print(f"Ratio: {top_ratio:.2f}x")

# Rank correlation
try:
    top_n = min(100, len(df))
    df_top = df.nlargest(top_n, 'coef_prox_exp_100k_true')
    rho, rho_p = stats.spearmanr(
        df_top['coef_prox_exp_100k_true'].abs(),
        df_top['coef_prox_exp_100k_rotate'].abs()
    )
    print(f"Top {top_n} genes rank correlation: ρ={rho:.3f}, p={rho_p:.4e}")
except Exception as e:
    print(f"⚠️  Rank correlation failed: {e}")
    rho_p = np.nan

# ============================================================
# VISUALIZATIONS
# ============================================================
Path("out_v3/plots").mkdir(parents=True, exist_ok=True)

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Panel A: Coefficient scatter
ax = axes[0, 0]
ax.scatter(df['coef_prox_exp_100k_rotate'], df['coef_prox_exp_100k_true'],
           alpha=0.3, s=1)
ax.axline((0, 0), slope=1, color='red', linestyle='--', label='y=x')
ax.set_xlabel('ROTATE coefficient')
ax.set_ylabel('TRUE coefficient')
ax.set_title('Proximity Coefficients')
ax.legend()

# Panel B: |Coefficient| distributions
ax = axes[0, 1]
ax.hist(df['coef_prox_exp_100k_true'].abs(), bins=50, alpha=0.5, label='TRUE', density=True)
ax.hist(df['coef_prox_exp_100k_rotate'].abs(), bins=50, alpha=0.5, label='ROTATE', density=True)
ax.set_xlabel('|Proximity coefficient|')
ax.set_ylabel('Density')
ks_p_str = f"{ks_p:.2e}" if not np.isnan(ks_p) else "NaN"
ax.set_title(f'KS test: p={ks_p_str}')
ax.legend()
ax.set_xlim(0, 0.05)  # Focus on main distribution

# Panel C: R² scatter
ax = axes[0, 2]
ax.scatter(df['r2_rotate'], df['r2_true'], alpha=0.3, s=1)
ax.axline((0, 0), slope=1, color='red', linestyle='--', label='y=x')
ax.set_xlabel('ROTATE R²')
ax.set_ylabel('TRUE R²')
ax.set_title('Model Fit Comparison')
ax.legend()

# Panel D: R² difference histogram
ax = axes[1, 0]
ax.hist(r2_diff, bins=50, edgecolor='black')
ax.axvline(0, color='red', linestyle='--', label='No difference')
ax.set_xlabel('R² (TRUE - ROTATE)')
ax.set_ylabel('Count')
p_r2_str = f"{p_r2:.2e}" if not np.isnan(p_r2) else "NaN"
ax.set_title(f'Mean diff={r2_diff.mean():.4f}, p={p_r2_str}')
ax.legend()

# Panel E: Top genes comparison
ax = axes[1, 1]
top_n = 50
df_plot = df.nlargest(top_n, 'coef_prox_exp_100k_true')
x = np.arange(len(df_plot))
ax.bar(x - 0.2, df_plot['coef_prox_exp_100k_true'].abs(), width=0.4, 
       label='TRUE', alpha=0.7)
ax.bar(x + 0.2, df_plot['coef_prox_exp_100k_rotate'].abs(), width=0.4,
       label='ROTATE', alpha=0.7)
ax.set_xlabel(f'Gene rank (top {top_n})')
ax.set_ylabel('|Proximity coefficient|')
ax.set_title('Top Genes: TRUE vs ROTATE')
ax.legend()

# Panel F: Cumulative distribution
ax = axes[1, 2]
sorted_true = np.sort(df['coef_prox_exp_100k_true'].abs())
sorted_rotate = np.sort(df['coef_prox_exp_100k_rotate'].abs())
ax.plot(sorted_true, np.linspace(0, 1, len(sorted_true)), label='TRUE')
ax.plot(sorted_rotate, np.linspace(0, 1, len(sorted_rotate)), label='ROTATE')
ax.set_xlabel('|Proximity coefficient|')
ax.set_ylabel('Cumulative probability')
ax.set_title(f'CDF Comparison (KS: p={ks_p_str})')
ax.legend()
ax.set_xlim(0, 0.1)

plt.tight_layout()
plt.savefig('out_v3/plots/true_vs_rotate_comparison.png', dpi=300)
print("\n✓ Saved plot to out_v3/plots/true_vs_rotate_comparison.png")

# ============================================================
# DECISION SUMMARY
# ============================================================
print("\n" + "="*60)
print("DECISION SUMMARY")
print("="*60)

go_criteria = []

# 1. Global distribution test
if not np.isnan(ks_p) and ks_p < 0.05:
    go_criteria.append("✓ TRUE ≠ ROTATE (KS test, p<0.05)")
    print(f"✓ TRUE coefficients differ from ROTATE (p={ks_p:.2e})")
else:
    go_criteria.append("✗ TRUE = ROTATE (KS test, p≥0.05)")
    print(f"✗ TRUE coefficients not distinguishable from ROTATE (p={ks_p})")

# 2. R² improvement
if not np.isnan(p_r2) and p_r2 < 0.05 and r2_diff.mean() > 0:
    go_criteria.append("✓ TRUE R² > ROTATE R² (p<0.05)")
    print(f"✓ TRUE models fit better than ROTATE (Δ R²={r2_diff.mean():.4f}, p={p_r2:.2e})")
else:
    go_criteria.append("✗ TRUE R² not better than ROTATE")
    print(f"✗ TRUE models don't fit better than ROTATE (p={p_r2})")

# 3. Top genes criterion
if not np.isnan(top_ratio) and top_ratio > 1.2:
    go_criteria.append("✓ Top 50 TRUE genes > ROTATE by 20%")
    print(f"✓ Top TRUE genes have 20%+ larger coefficients than ROTATE ({top_ratio:.2f}x)")
else:
    go_criteria.append("✗ Top TRUE genes not substantially larger")
    print(f"✗ Top TRUE genes not substantially larger than ROTATE ({top_ratio:.2f}x)")

# Final decision
n_pass = sum('✓' in c for c in go_criteria)
print(f"\nCriteria passed: {n_pass}/3")

if n_pass >= 2:
    print("\n🟢 GO: Spatial signal detected. Proceed with case studies.")
    decision = "GO"
elif n_pass == 1:
    print("\n🟡 INVESTIGATE: Weak/ambiguous signal. Run follow-up analyses.")
    decision = "INVESTIGATE"
else:
    print("\n🔴 STOP: No clear spatial signal. Effect may be noise.")
    decision = "STOP"

# Save decision
with open("out_v3/decision.txt", "w") as f:
    f.write(f"Decision: {decision}\n")
    f.write(f"Criteria passed: {n_pass}/3\n\n")
    for c in go_criteria:
        f.write(f"{c}\n")
    f.write(f"\nData quality:\n")
    f.write(f"  Genes compared: {len(df)}\n")
    f.write(f"  KS p-value: {ks_p}\n")
    f.write(f"  R² test p-value: {p_r2}\n")
    f.write(f"  Top 50 ratio: {top_ratio:.2f}\n")

print("\n✓ Decision saved to out_v3/decision.txt")
