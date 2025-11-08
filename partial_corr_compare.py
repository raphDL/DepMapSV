#!/usr/bin/env python3
"""
Compare partial correlation metrics across multiple conditions.

Reads partial_corr_summary.csv from each base folder and prints a comparison table.
"""

import argparse
import pandas as pd
from pathlib import Path
from partial_corr_from_design import compute_partial_corr_delta


def main():
    parser = argparse.ArgumentParser(
        description="Compare partial correlation metrics across conditions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("base_folders", nargs="+", help="Base folders to compare")
    parser.add_argument("--min-rows", type=int, default=20, help="Minimum rows per gene")
    parser.add_argument("--recompute", action="store_true", help="Recompute even if summary exists")
    
    args = parser.parse_args()
    
    results = []
    
    for folder_str in args.base_folders:
        folder = Path(folder_str)
        if not folder.exists():
            print(f"Warning: {folder} does not exist, skipping")
            continue
        
        summary_path = folder / "partial_corr_summary.csv"
        
        # Compute if needed
        if args.recompute or not summary_path.exists():
            try:
                compute_partial_corr_delta(folder, min_rows=args.min_rows)
            except Exception as e:
                print(f"Error processing {folder}: {e}")
                continue
        
        # Load summary
        if summary_path.exists():
            try:
                summary = pd.read_csv(summary_path).iloc[0].to_dict()
                summary["condition"] = folder.name
                results.append(summary)
            except Exception as e:
                print(f"Error reading {summary_path}: {e}")
                continue
    
    if not results:
        print("No valid results found")
        return
    
    # Create comparison table
    df = pd.DataFrame(results)
    
    # Reorder columns for readability
    cols = ["condition", "n_genes", "median_partial_delta_abs_corr", 
            "ci_low", "ci_high", "mean_partial_delta_abs_corr", "frac_positive"]
    cols = [c for c in cols if c in df.columns]
    df = df[cols]
    
    # Print table
    print("\n" + "=" * 100)
    print("PARTIAL CORRELATION COMPARISON")
    print("=" * 100)
    print(df.to_string(index=False))
    print("=" * 100)
    
    # Decision rule check
    if len(df) > 1:
        true_row = df[df["condition"].str.contains("true|baseline|main", case=False, na=False)]
        if len(true_row) == 0:
            true_row = df.iloc[0:1]  # Use first as "true" if no match
        
        if len(true_row) > 0:
            true_val = true_row.iloc[0]
            true_lower = true_val["ci_low"]
            
            # Check if TRUE lower bound > any shuffle upper bound
            shuffles = df[~df.index.isin(true_row.index)]
            if len(shuffles) > 0:
                max_shuffle_upper = shuffles["ci_high"].max()
                
                print(f"\nDecision Rule:")
                print(f"  TRUE lower bound: {true_lower:.4f}")
                print(f"  Max shuffle upper bound: {max_shuffle_upper:.4f}")
                if true_lower > max_shuffle_upper:
                    print(f"  ✓ PROCEED: TRUE CI lower > all shuffle CI upper")
                else:
                    print(f"  ✗ PAUSE: TRUE CI lower ≤ at least one shuffle CI upper")


if __name__ == "__main__":
    main()

