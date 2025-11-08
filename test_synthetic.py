#!/usr/bin/env python3
"""Generate minimal synthetic test data."""
import pandas as pd
from pathlib import Path

# Create test directory
Path("test_data").mkdir(exist_ok=True)

# 1. Genes BED (5 genes on chr1, 3 on chr2)
genes = pd.DataFrame({
    "chrom": ["chr1"]*5 + ["chr2"]*3,
    "start": [1_000_000, 2_000_000, 3_000_000, 4_000_000, 5_000_000,
              1_000_000, 2_000_000, 3_000_000],
    "end":   [1_010_000, 2_010_000, 3_010_000, 4_010_000, 5_010_000,
              1_010_000, 2_010_000, 3_010_000],
    "gene": ["GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E",
             "GENE_F", "GENE_G", "GENE_H"],
    "strand": ["+"]*8
})
genes.to_csv("test_data/genes.bed", sep="\t", index=False, header=False)

# 2. CNV (2 cell lines, 3 segments each)
cnv = pd.DataFrame({
    "ModelID": ["ACH-000001"]*3 + ["ACH-000002"]*3,
    "Chromosome": ["chr1"]*6,
    "Start": [0, 1_500_000, 3_500_000, 0, 2_500_000, 4_500_000],
    "End": [1_500_000, 3_500_000, 6_000_000, 2_500_000, 4_500_000, 6_000_000],
    "CopyNumber": [2.0, 4.5, 1.2, 2.1, 5.2, 0.8]
})
cnv.to_csv("test_data/cnv.csv", sep="\t", index=False)

# 3. SV BEDPE (2 cell lines, separate files - standard BEDPE format)
# Cell 1: SVs near GENE_A (1Mb) and GENE_D (4Mb)
# Cell 2: SVs near GENE_B (2Mb) and GENE_E (5Mb)

# Create separate BEDPE files per cell line (more realistic)
sv_data_1 = [
    ["chr1", 1100000, 1100001, "chr2", 5000000, 5000001, "TRA", ".", "+", "-", "ACH-000001"],
    ["chr1", 4100000, 4100001, "chr3", 6000000, 6000001, "TRA", ".", "+", "-", "ACH-000001"],
    ["chr2", 1200000, 1200001, "chr2", 1300000, 1300001, "INV", ".", "+", "+", "ACH-000001"],
]
sv1 = pd.DataFrame(sv_data_1)
sv1.to_csv("test_data/ACH-000001.bedpe", sep="\t", index=False, header=False)

sv_data_2 = [
    ["chr1", 2050000, 2050001, "chr4", 7000000, 7000001, "TRA", ".", "+", "-", "ACH-000002"],
    ["chr1", 5200000, 5200001, "chr5", 8000000, 8000001, "TRA", ".", "+", "-", "ACH-000002"],
    ["chr2", 2100000, 2100001, "chr2", 2200000, 2200001, "INV", ".", "+", "+", "ACH-000002"],
]
sv2 = pd.DataFrame(sv_data_2)
sv2.to_csv("test_data/ACH-000002.bedpe", sep="\t", index=False, header=False)

print("✓ Created test_data/")
print("  - genes.bed: 8 genes (5 on chr1, 3 on chr2)")
print("  - cnv.csv: 2 cell lines × 3 segments = 6 rows")
print("  - ACH-000001.bedpe: 3 SVs")
print("  - ACH-000002.bedpe: 3 SVs")
print("  - Total: 2 cell lines × 3 SVs = 6 events")
print("\nExpected results:")
print("  - Overlap: 2 cell lines")
print("  - Design matrix: 8 genes × 2 cells = 16 rows")
print("  - GENE_A: bp_dist ~ 100kb (SV at 1.1Mb, gene at 1.0Mb)")
print("  - GENE_B: bp_dist ~ 50kb (SV at 2.05Mb, gene at 2.0Mb)")
print("  - GENE_H: bp_dist >> 1Mb (no nearby SVs on chr2)")

