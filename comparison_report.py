#!/usr/bin/env python3
"""
Comparison Report - Generate final summary report with upper bounds and case study table.

Reports:
- Config and sample sizes
- Global metrics with CIs
- Excess signal
- Case study table (all selected genes, not just prettiest)
- Stop/Go gates assessment
"""
import argparse
import json
import logging
import os
import sys
from pathlib import Path

import pandas as pd


def setup_logging(verbosity: int = 1):
    """Setup logging with verbosity control."""
    level = logging.WARNING if verbosity <= 0 else (logging.INFO if verbosity == 1 else logging.DEBUG)
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S",
        handlers=[logging.StreamHandler(sys.stderr)],
        force=True,
    )


def load_metrics(root_dir: str) -> dict:
    """Load all metrics from output directories."""
    metrics = {}
    
    # Global metrics
    excess_path = os.path.join(root_dir, "excess_signal.json")
    if os.path.exists(excess_path):
        with open(excess_path, "r") as f:
            metrics["global"] = json.load(f)
    
    # Prevalence-matched directional
    dir_path = os.path.join(root_dir, "prevalence_matched_directional.csv")
    if os.path.exists(dir_path):
        metrics["directional"] = pd.read_csv(dir_path)
    
    # Case studies
    case_path = os.path.join(root_dir, "case_studies", "case_studies_validated.csv")
    if not os.path.exists(case_path):
        case_path = os.path.join(root_dir, "case_studies", "case_studies.csv")
    if os.path.exists(case_path):
        metrics["case_studies"] = pd.read_csv(case_path)
    
    # TRUE metrics
    true_metrics_path = os.path.join(root_dir, "true", "metrics.json")
    if os.path.exists(true_metrics_path):
        with open(true_metrics_path, "r") as f:
            metrics["true"] = json.load(f)
    
    return metrics


def assess_stop_go(metrics: dict) -> dict:
    """
    Assess Stop/Go gates.
    
    Go (Patterns): ≥20 genes pass, ≥10 at FDR<0.10, ≥3–5 with orthogonal validation.
    Pivot (methods/null): fewer than above → report upper bounds + framework.
    """
    case_studies = metrics.get("case_studies")
    if case_studies is None or len(case_studies) == 0:
        return {
            "status": "PIVOT",
            "reason": "No case studies selected",
            "n_genes_pass": 0,
            "n_fdr_10": 0,
            "n_validated": 0
        }
    
    n_genes_pass = len(case_studies)
    n_fdr_10 = (case_studies["q_fdr"] < 0.10).sum() if "q_fdr" in case_studies.columns else 0
    n_validated = (case_studies["passes_validation"]).sum() if "passes_validation" in case_studies.columns else 0
    
    if n_genes_pass >= 20 and n_fdr_10 >= 10 and n_validated >= 3:
        status = "GO"
        reason = "All criteria met"
    else:
        status = "PIVOT"
        reason = f"Criteria not met: n_pass={n_genes_pass} (need ≥20), n_fdr={n_fdr_10} (need ≥10), n_validated={n_validated} (need ≥3)"
    
    return {
        "status": status,
        "reason": reason,
        "n_genes_pass": int(n_genes_pass),
        "n_fdr_10": int(n_fdr_10),
        "n_validated": int(n_validated)
    }


def generate_report(metrics: dict, root_dir: str, out_path: str):
    """Generate final comparison report."""
    lines = []
    
    lines.append("=" * 80)
    lines.append("SV/CN Bias V3 - Final Comparison Report")
    lines.append("=" * 80)
    lines.append("")
    
    # Config and sample sizes
    lines.append("## Configuration and Sample Sizes")
    lines.append("")
    if "true" in metrics:
        true_metrics = metrics["true"]
        lines.append(f"- Genes analyzed: {true_metrics.get('n_genes', 'N/A')}")
        lines.append(f"- Cell lines: {true_metrics.get('n_cells', 'N/A')}")
        lines.append(f"- Gene-cell pairs: {true_metrics.get('n_pairs', 'N/A')}")
        lines.append(f"- Genes fitted: {true_metrics.get('n_genes_fitted', 'N/A')}")
        lines.append(f"- Kernel used: {true_metrics.get('kernel_used', 'N/A')}")
    lines.append("")
    
    # Global metrics
    lines.append("## Global Metrics")
    lines.append("")
    if "global" in metrics:
        global_metrics = metrics["global"]
        
        # Prevalence-matched directional
        if "prevalence_matched_directional" in global_metrics:
            dir_metrics = global_metrics["prevalence_matched_directional"]
            lines.append("### Prevalence-Matched Directional")
            lines.append(f"- TRUE:   {dir_metrics.get('dir_true', 'N/A'):.4f}")
            lines.append(f"- ROTATE: {dir_metrics.get('dir_rotate', 'N/A'):.4f}")
            lines.append(f"- EXCESS: {dir_metrics.get('excess', 'N/A'):+.4f}")
            lines.append(f"- 95% CI: [{dir_metrics.get('ci_lower', 'N/A'):+.4f}, {dir_metrics.get('ci_upper', 'N/A'):+.4f}]")
            lines.append("")
        
        # Coef→effect consistency
        if "coef_effect_consistency" in global_metrics:
            coef_effect = global_metrics["coef_effect_consistency"]
            lines.append("### Coef→Effect Consistency")
            if "true_directional" in coef_effect:
                td = coef_effect["true_directional"]
                lines.append(f"- TRUE directional: ρ={td.get('rho', 'N/A'):.4f}, p={td.get('p', 'N/A'):.4e}, n={td.get('n', 'N/A')}")
            lines.append("")
        
        # Excess signal
        if "excess_signal" in global_metrics:
            excess = global_metrics["excess_signal"]
            lines.append("### Excess Signal")
            if "corr_drop_excess" in excess and excess["corr_drop_excess"]:
                corr_excess = excess["corr_drop_excess"]
                lines.append(f"- Δ|corr| excess (median): {corr_excess.get('median', 'N/A'):.4f}")
                lines.append(f"- Δ|corr| excess (Q25-Q75): [{corr_excess.get('q25', 'N/A'):.4f}, {corr_excess.get('q75', 'N/A'):.4f}]")
            lines.append("")
    
    # Upper bounds
    lines.append("## Upper Bounds on Global Effect")
    lines.append("")
    if "global" in metrics and "prevalence_matched_directional" in metrics["global"]:
        dir_metrics = metrics["global"]["prevalence_matched_directional"]
        excess = dir_metrics.get("excess", 0)
        ci_upper = dir_metrics.get("ci_upper", 0)
        lines.append(f"Prevalence-matched directional excess: {excess:.4f}")
        lines.append(f"Upper bound (95% CI): {ci_upper:.4f}")
        lines.append("")
        lines.append("**Interpretation:** The global effect of SV proximity on dependency-CN")
        lines.append("correlation is bounded above by the excess directional fraction.")
    lines.append("")
    
    # Case studies
    lines.append("## Case Studies")
    lines.append("")
    case_studies = metrics.get("case_studies")
    if case_studies is not None and len(case_studies) > 0:
        lines.append(f"Total case studies selected: {len(case_studies)}")
        lines.append("")
        lines.append("### All Selected Case Studies")
        lines.append("")
        
        # Create table
        cols_to_show = ["gene", "corr_drop_true", "q_fdr"]
        if "passes_validation" in case_studies.columns:
            cols_to_show.append("passes_validation")
        if "n_validations" in case_studies.columns:
            cols_to_show.append("n_validations")
        
        available_cols = [c for c in cols_to_show if c in case_studies.columns]
        table_df = case_studies[available_cols].copy()
        
        # Format for display
        if "corr_drop_true" in table_df.columns:
            table_df["corr_drop_true"] = table_df["corr_drop_true"].apply(lambda x: f"{x:.4f}")
        if "q_fdr" in table_df.columns:
            table_df["q_fdr"] = table_df["q_fdr"].apply(lambda x: f"{x:.4f}")
        
        lines.append(table_df.to_string(index=False))
        lines.append("")
    else:
        lines.append("No case studies selected.")
        lines.append("")
    
    # Stop/Go assessment
    lines.append("## Stop/Go Gates Assessment")
    lines.append("")
    stop_go = assess_stop_go(metrics)
    lines.append(f"**Status:** {stop_go['status']}")
    lines.append(f"**Reason:** {stop_go['reason']}")
    lines.append("")
    lines.append("Criteria:")
    lines.append(f"- Genes passing all criteria: {stop_go['n_genes_pass']} (need ≥20)")
    lines.append(f"- Genes at FDR < 0.10: {stop_go['n_fdr_10']} (need ≥10)")
    lines.append(f"- Genes with orthogonal validation: {stop_go['n_validated']} (need ≥3-5)")
    lines.append("")
    
    if stop_go["status"] == "GO":
        lines.append("**Conclusion:** Patterns detected. Proceed with case study analysis.")
    else:
        lines.append("**Conclusion:** Pivot to methods/null framework. Report upper bounds.")
    lines.append("")
    
    # Write report
    report_text = "\n".join(lines)
    
    with open(out_path, "w") as f:
        f.write(report_text)
    
    logging.info(f"Saved report to {out_path}")
    
    # Also print to stdout
    print(report_text)


def main():
    ap = argparse.ArgumentParser(
        description="Generate final comparison report"
    )
    ap.add_argument("--root", required=True,
                    help="Root output directory (out_v3/)")
    ap.add_argument("--out", required=True,
                    help="Output report file path")
    ap.add_argument("-v", "--verbose", action="count", default=1,
                    help="Increase verbosity")
    args = ap.parse_args()
    
    setup_logging(args.verbose)
    
    # Load metrics
    logging.info("Loading metrics...")
    metrics = load_metrics(args.root)
    
    # Generate report
    logging.info("Generating report...")
    generate_report(metrics, args.root, args.out)


if __name__ == "__main__":
    main()

