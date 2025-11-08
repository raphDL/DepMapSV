#!/usr/bin/env python3
"""
Compare proximity-only partial correlation metrics across conditions.
"""

import argparse
import pandas as pd
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Compare proximity-only partial correlation metrics",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("folders", nargs="+", help="Base folders to compare")
    
    args = parser.parse_args()
    
    results = []
    
    for folder_str in args.folders:
        folder = Path(folder_str)
        summary_path = folder / "prox_only_partial_corr_summary.csv"
        
        if not summary_path.exists():
            print(f"Warning: {summary_path} not found, skipping {folder.name}")
            continue
        
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
    
    # Reorder columns
    cols = ["condition", "n_genes", "median_prox_only_delta", 
            "ci_low", "ci_high", "mean_prox_only_delta", "frac_positive"]
    cols = [c for c in cols if c in df.columns]
    df = df[cols]
    
    # Print table
    print("\n" + "=" * 100)
    print("PROXIMITY-ONLY PARTIAL CORRELATION COMPARISON")
    print("=" * 100)
    print(df.to_string(index=False))
    print("=" * 100)
    
    # Decision rule
    if len(df) > 1:
        true_row = df[df["condition"].str.contains("true|baseline|main", case=False, na=False)]
        if len(true_row) == 0:
            true_row = df.iloc[0:1]
        
        if len(true_row) > 0:
            true_val = true_row.iloc[0]
            true_lower = true_val["ci_low"]
            true_median = true_val["median_prox_only_delta"]
            
            shuffles = df[~df.index.isin(true_row.index)]
            if len(shuffles) > 0:
                max_shuffle_upper = shuffles["ci_high"].max()
                max_shuffle_median = shuffles["median_prox_only_delta"].max()
                
                print(f"\nDecision Rule:")
                print(f"  TRUE median: {true_median:.4f} [CI: {true_lower:.4f}, {true_val['ci_high']:.4f}]")
                print(f"  Max shuffle median: {max_shuffle_median:.4f} [CI upper: {max_shuffle_upper:.4f}]")
                
                if true_lower > max_shuffle_upper:
                    print(f"  ✓ PROCEED: TRUE CI lower > all shuffle CI upper")
                elif true_median > max_shuffle_median:
                    print(f"  ⚠ MARGINAL: TRUE median > max shuffle median, but CIs overlap")
                else:
                    print(f"  ✗ PAUSE: TRUE ≤ shuffles on proximity-only metric")


if __name__ == "__main__":
    main()

