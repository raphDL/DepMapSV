import pandas as pd, subprocess

sig = pd.read_csv("figs_full/psi_line_significance.csv")
rel = pd.read_csv("figs_full/psi_line_split_half.csv")

# Get top lines by psi_diff (TRUE > ROTATE)
keep = sig[sig["psi_diff"] > 0].sort_values("psi_diff", ascending=False).head(5)
print("Top 5 lines by psi_diff:")
print(keep[["cell_line", "psi_diff", "psi_true", "psi_rot", "q_fdr"]])

# Merge with reliability data
keep = keep.merge(rel, on="cell_line", how="left")
top_lines = keep["cell_line"].tolist()[:3]  # Top 3 for panels
print(f"\nGenerating panels for top 3 lines: {top_lines}")

# pick a few genes you already know show proximity patterns in panels
genes = ["GFPT2","RPS6KB2","KCNQ1","AAMDC","ALK"]

# Generate panels (gene_panel.py doesn't support --focus-cell, so we'll generate standard panels)
# and note which cell lines are high-PSI in the filename
for cl in top_lines:
    for g in genes:
        outfile = f"figs_full/panel_{g}_{cl}.png"
        subprocess.run(["python","gene_panel.py","out_focus_true/unstable",g,outfile], 
                      check=False)
        print(f"Generated {outfile} (high-PSI line: {cl})")

print("\nPanels written to figs_full/")
print("Note: Panels show all cell lines; high-PSI lines are listed in filenames")

