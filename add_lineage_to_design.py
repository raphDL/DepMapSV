#!/usr/bin/env python3
# Merge a cell_line -> lineage CSV into design_matrix.csv
# Usage:
#   python add_lineage_to_design.py out_pilot_final/unstable cell_line,lineage path/to/cell_line_lineage_map.csv
#   python add_lineage_to_design.py out_pilot_final/stable   cell_line,lineage path/to/cell_line_lineage_map.csv

import sys, os, pandas as pd

def main():
    if len(sys.argv) != 4:
        print("Usage: python add_lineage_to_design.py <setdir> <key_cols> <mapping_csv>")
        print("Example key_cols: cell_line,lineage")
        sys.exit(1)
    setdir, key_cols, mapping_csv = sys.argv[1], sys.argv[2], sys.argv[3]
    key_src, key_dst = key_cols.split(",")[0], key_cols.split(",")[1]

    dm_path = os.path.join(setdir, "design_matrix.csv")
    if not os.path.exists(dm_path):
        sys.exit(f"ERROR: {dm_path} not found")

    dm = pd.read_csv(dm_path)
    m = pd.read_csv(mapping_csv)
    if key_src not in m.columns or key_dst not in m.columns:
        sys.exit(f"ERROR: mapping CSV must contain columns {key_src} and {key_dst}")

    out = dm.merge(m[[key_src, key_dst]], left_on="cell_line", right_on=key_src, how="left")
    out = out.drop(columns=[key_src]).rename(columns={key_dst: "lineage"})
    out.to_csv(dm_path, index=False)
    print(f"Appended 'lineage' to {dm_path}")

if __name__ == "__main__":
    main()
