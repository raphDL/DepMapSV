#!/usr/bin/env python3
"""
Create a 3-panel case-study figure for a single gene.

Usage: python gene_panel.py out_pilot_final/unstable GFPT2 figs_pilot/GFPT2_panel.png
"""

import sys, pandas as pd, numpy as np, matplotlib.pyplot as plt, os

SETDIR, GENE, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

dm   = pd.read_csv(os.path.join(SETDIR,"design_matrix.csv"))
coef = pd.read_csv(os.path.join(SETDIR,"models_coefficients.csv")).rename(
    columns={"cn":"b_cn","bp_within_100000":"b_w100k","bp_within_1000000":"b_w1m","inv_bp":"b_inv"}
)

sub = dm.query("gene == @GENE").merge(coef[coef["gene"]==GENE][["gene","b_w100k","b_w1m","b_inv"]], on="gene", how="left")
# Use standardized columns if present (same scale as model)
x_w100k = "bp_within_100000_std" if "bp_within_100000_std" in sub.columns else "bp_within_100000"
x_w1m   = "bp_within_1000000_std" if "bp_within_1000000_std" in sub.columns else "bp_within_1000000"
x_inv   = "inv_bp_std" if "inv_bp_std" in sub.columns else ("inv_bp" if "inv_bp" in sub.columns else None)
sub["prox_contrib"] = (sub["b_w100k"]*sub.get(x_w100k, 0.0).astype(float) +
                       sub["b_w1m"]  *sub.get(x_w1m, 0.0).astype(float) +
                       (sub["b_inv"]*sub.get(x_inv, 0.0).astype(float) if x_inv and x_inv in sub.columns else 0.0))
sub["prox_only_corrected"] = sub["dependency"] - sub["prox_contrib"]

plt.figure(figsize=(9,3.4))
plt.subplot(1,3,1)
plt.scatter(sub["bp_dist"], sub["dependency"], s=6, alpha=0.6)
plt.xscale("symlog", linthresh=1e5); plt.xlabel("Nearest BP distance (bp)"); plt.ylabel("Dependency"); plt.title(GENE)

plt.subplot(1,3,2)
plt.scatter(sub["cn"], sub["dependency"], s=6, label="original", alpha=0.6)
plt.scatter(sub["cn"], sub["prox_only_corrected"], s=6, alpha=0.6, label="prox-only removed", color='#2a9d8f')
plt.xlabel("Copy number"); plt.ylabel("Dependency"); plt.legend(loc="best", fontsize=8)

plt.subplot(1,3,3)
plt.hist(sub["prox_contrib"].abs(), bins=30, edgecolor='black', alpha=0.7)
plt.xlabel("|proximity contribution|"); plt.ylabel("cells")
plt.tight_layout(); os.makedirs(os.path.dirname(OUT), exist_ok=True)
plt.savefig(OUT, dpi=150)
if OUT.endswith(".png"):
    plt.savefig(OUT.replace(".png", ".pdf"))
elif OUT.endswith(".pdf"):
    plt.savefig(OUT.replace(".pdf", ".png"), dpi=150)
plt.close()
print("wrote", OUT)

