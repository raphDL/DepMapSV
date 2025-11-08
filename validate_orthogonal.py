#!/usr/bin/env python3
"""
Orthogonal Validation - RNAi concordance, drug sensitivity, pathway enrichment.

For case studies, validate with:
- RNAi: Δcorr/ΔAUROC (CRISPR vs RNAi) before/after prox-only correction
- Drug: if target-matched compound exists, association improves
- Pathway: GO/KEGG enrichment for DDR/replication stress

Need ≥2/3 validations to pass.
"""
import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr
from sklearn.metrics import roc_auc_score

try:
    from gseapy import enrichr
except ImportError:
    enrichr = None
    logging.warning("gseapy not available, pathway enrichment will be skipped")


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


def validate_rnai(case_studies: pd.DataFrame, rnai_path: str, 
                  true_dir: str, kernel_col: str) -> pd.DataFrame:
    """
    Validate with RNAi: Δcorr/ΔAUROC (CRISPR vs RNAi) before/after prox-only correction.
    
    Returns DataFrame with RNAi validation columns.
    """
    if not os.path.exists(rnai_path):
        logging.warning(f"RNAi file not found: {rnai_path}")
        return case_studies.assign(
            rnai_delta_corr=np.nan,
            rnai_delta_auroc=np.nan,
            rnai_validated=False
        )
    
    logging.info("Loading RNAi data...")
    rnai = pd.read_csv(rnai_path)
    required = ["gene", "cell_line", "dependency"]
    if not all(c in rnai.columns for c in required):
        logging.warning(f"RNAi file missing required columns: {required}")
        return case_studies.assign(
            rnai_delta_corr=np.nan,
            rnai_delta_auroc=np.nan,
            rnai_validated=False
        )
    
    # Load CRISPR data
    corr_path = os.path.join(true_dir, "dependency_corrected.csv")
    crispr = pd.read_csv(corr_path)
    
    # Merge RNAi with CRISPR
    merged = crispr.merge(rnai, on=["gene", "cell_line"], 
                          suffixes=("_crispr", "_rnai"), how="inner")
    
    # Per-gene validation
    rnai_results = []
    for gene in case_studies["gene"]:
        g_data = merged[merged["gene"] == gene].dropna(
            subset=["dependency_crispr", "dependency_rnai", "prox_only_corrected"]
        )
        
        if len(g_data) < 20:
            rnai_results.append({
                "gene": gene,
                "rnai_delta_corr": np.nan,
                "rnai_delta_auroc": np.nan,
                "rnai_validated": False
            })
            continue
        
        # Before correction: corr(CRISPR, RNAi)
        r_before = g_data[["dependency_crispr", "dependency_rnai"]].corr().iloc[0, 1]
        
        # After prox-only correction: corr(prox_only_corrected, RNAi)
        r_after = g_data[["prox_only_corrected", "dependency_rnai"]].corr().iloc[0, 1]
        
        delta_corr = abs(r_after) - abs(r_before) if pd.notna(r_after) and pd.notna(r_before) else np.nan
        
        # AUROC improvement (if we have essential/nonessential labels)
        # For now, skip AUROC
        delta_auroc = np.nan
        
        rnai_results.append({
            "gene": gene,
            "rnai_delta_corr": float(delta_corr) if pd.notna(delta_corr) else np.nan,
            "rnai_delta_auroc": delta_auroc,
            "rnai_validated": pd.notna(delta_corr) and delta_corr > 0
        })
    
    rnai_df = pd.DataFrame(rnai_results)
    return case_studies.merge(rnai_df, on="gene", how="left")


def validate_drug(case_studies: pd.DataFrame, drug_path: str,
                  true_dir: str) -> pd.DataFrame:
    """
    Validate with drug sensitivity: if target-matched compound exists, association improves.
    
    Returns DataFrame with drug validation columns.
    """
    if not os.path.exists(drug_path):
        logging.warning(f"Drug sensitivity file not found: {drug_path}")
        return case_studies.assign(
            drug_effect_size=np.nan,
            drug_p_value=np.nan,
            drug_validated=False
        )
    
    logging.info("Loading drug sensitivity data...")
    drug = pd.read_csv(drug_path)
    required = ["compound", "cell_line", "sensitivity", "target"]
    if not all(c in drug.columns for c in required):
        logging.warning(f"Drug file missing required columns: {required}")
        return case_studies.assign(
            drug_effect_size=np.nan,
            drug_p_value=np.nan,
            drug_validated=False
        )
    
    # Load CRISPR data
    corr_path = os.path.join(true_dir, "dependency_corrected.csv")
    crispr = pd.read_csv(corr_path)
    
    # Per-gene validation
    drug_results = []
    for gene in case_studies["gene"]:
        # Find drugs targeting this gene
        target_drugs = drug[drug["target"].str.contains(gene, case=False, na=False)]
        
        if len(target_drugs) == 0:
            drug_results.append({
                "gene": gene,
                "drug_effect_size": np.nan,
                "drug_p_value": np.nan,
                "drug_validated": False
            })
            continue
        
        # Use first matching drug
        drug_name = target_drugs["compound"].iloc[0]
        drug_data = target_drugs[target_drugs["compound"] == drug_name].merge(
            crispr[crispr["gene"] == gene][["cell_line", "prox_only_corrected"]],
            on="cell_line", how="inner"
        ).dropna(subset=["sensitivity", "prox_only_corrected"])
        
        if len(drug_data) < 20:
            drug_results.append({
                "gene": gene,
                "drug_effect_size": np.nan,
                "drug_p_value": np.nan,
                "drug_validated": False
            })
            continue
        
        # Correlation between drug sensitivity and prox-only corrected dependency
        rho, p = spearmanr(drug_data["sensitivity"], drug_data["prox_only_corrected"])
        
        drug_results.append({
            "gene": gene,
            "drug_effect_size": float(rho) if pd.notna(rho) else np.nan,
            "drug_p_value": float(p) if pd.notna(p) else np.nan,
            "drug_validated": pd.notna(p) and p < 0.05 and pd.notna(rho) and rho < 0
        })
    
    drug_df = pd.DataFrame(drug_results)
    return case_studies.merge(drug_df, on="gene", how="left")


def validate_pathway(case_studies: pd.DataFrame, 
                    ddr_genes: set = None, replication_stress_genes: set = None) -> pd.DataFrame:
    """
    Validate with pathway membership: DDR/replication stress.
    
    Returns DataFrame with pathway validation columns.
    """
    if ddr_genes is None:
        # Default DDR genes (simplified list)
        ddr_genes = {
            "ATM", "ATR", "CHEK1", "CHEK2", "BRCA1", "BRCA2", "RAD51", "RAD52",
            "PARP1", "PARP2", "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE"
        }
    
    if replication_stress_genes is None:
        # Default replication stress genes
        replication_stress_genes = {
            "ATR", "CHEK1", "RPA1", "RPA2", "RPA3", "PCNA", "RFC1", "RFC2"
        }
    
    # Check pathway membership
    case_studies["in_ddr"] = case_studies["gene"].isin(ddr_genes)
    case_studies["in_replication_stress"] = case_studies["gene"].isin(replication_stress_genes)
    case_studies["pathway_validated"] = case_studies["in_ddr"] | case_studies["in_replication_stress"]
    
    # Optional: use enrichr for GO/KEGG enrichment
    if enrichr and len(case_studies) > 0:
        try:
            gene_list = case_studies["gene"].tolist()
            # Run enrichr (this requires internet connection)
            # For now, skip automated enrichment
            pass
        except Exception as e:
            logging.debug(f"Enrichr failed: {e}")
    
    return case_studies


def main():
    ap = argparse.ArgumentParser(
        description="Validate case studies with orthogonal data"
    )
    ap.add_argument("--cases", required=True,
                    help="Case studies CSV from select_case_studies.py")
    ap.add_argument("--true", required=True,
                    help="TRUE output directory")
    ap.add_argument("--rnai", default=None,
                    help="RNAi dependency file (long: gene,cell_line,dependency)")
    ap.add_argument("--drug", default=None,
                    help="Drug sensitivity file (compound,cell_line,sensitivity,target)")
    ap.add_argument("--ddr-genes", default=None,
                    help="File with DDR genes (one per line)")
    ap.add_argument("--replication-stress-genes", default=None,
                    help="File with replication stress genes (one per line)")
    ap.add_argument("--kernel", default="prox_exp_100k",
                    help="Kernel column name (default: prox_exp_100k)")
    ap.add_argument("--out", required=True,
                    help="Output directory")
    ap.add_argument("-v", "--verbose", action="count", default=1,
                    help="Increase verbosity")
    args = ap.parse_args()
    
    setup_logging(args.verbose)
    
    os.makedirs(args.out, exist_ok=True)
    
    # Load case studies
    logging.info("Loading case studies...")
    case_studies = pd.read_csv(args.cases)
    logging.info(f"Loaded {len(case_studies)} case studies")
    
    # RNAi validation
    if args.rnai:
        logging.info("Validating with RNAi...")
        case_studies = validate_rnai(case_studies, args.rnai, args.true, args.kernel)
    else:
        case_studies = case_studies.assign(
            rnai_delta_corr=np.nan,
            rnai_delta_auroc=np.nan,
            rnai_validated=False
        )
    
    # Drug validation
    if args.drug:
        logging.info("Validating with drug sensitivity...")
        case_studies = validate_drug(case_studies, args.drug, args.true)
    else:
        case_studies = case_studies.assign(
            drug_effect_size=np.nan,
            drug_p_value=np.nan,
            drug_validated=False
        )
    
    # Pathway validation
    ddr_genes = None
    if args.ddr_genes and os.path.exists(args.ddr_genes):
        ddr_genes = set(pd.read_csv(args.ddr_genes, header=None)[0].astype(str))
    
    replication_stress_genes = None
    if args.replication_stress_genes and os.path.exists(args.replication_stress_genes):
        replication_stress_genes = set(pd.read_csv(args.replication_stress_genes, header=None)[0].astype(str))
    
    logging.info("Validating with pathway membership...")
    case_studies = validate_pathway(case_studies, ddr_genes, replication_stress_genes)
    
    # Count validations (need ≥2/3)
    validation_cols = ["rnai_validated", "drug_validated", "pathway_validated"]
    case_studies["n_validations"] = case_studies[validation_cols].sum(axis=1)
    case_studies["passes_validation"] = case_studies["n_validations"] >= 2
    
    # Save
    out_path = os.path.join(args.out, "case_studies_validated.csv")
    case_studies.to_csv(out_path, index=False)
    logging.info(f"Saved validated case studies to {out_path}")
    
    # Summary
    print(f"\n=== Orthogonal Validation Summary ===")
    print(f"Total case studies: {len(case_studies)}")
    print(f"RNAi validated: {case_studies['rnai_validated'].sum()}")
    print(f"Drug validated: {case_studies['drug_validated'].sum()}")
    print(f"Pathway validated: {case_studies['pathway_validated'].sum()}")
    print(f"Pass validation (≥2/3): {case_studies['passes_validation'].sum()}")


if __name__ == "__main__":
    main()

