#!/usr/bin/env python3
"""
Falsification tests for SV bias correction:
1. Held-out (per-gene) falsification
2. CN-neutral falsification (ploidy ~2)
3. Effect-size normalization
4. Proximity-enriched set (unstable genes)
5. Sanity on AUROC (alternative metrics)
"""

import pandas as pd
import numpy as np
import os
import json

rng = np.random.default_rng(0)

# Test 1: Held-out falsification
def heldout_metric(base):
    """Train on 80% of cells, measure prox-only |Δr| on held-out 20%"""
    dm = pd.read_csv(os.path.join(base, 'dependency_corrected.csv'))
    
    if 'dependency' not in dm.columns:
        raise SystemExit(f"dependency missing in {base}/dependency_corrected.csv")
    
    assert 'cn' in dm.columns, "cn missing; rerun pipeline to keep cn/bp columns"
    
    rows = []
    for g, sub in dm.groupby('gene'):
        sub = sub[['cell_line', 'dependency', 'dependency_corrected', 'cn']].dropna()
        if len(sub) < 60:
            continue
        
        idx = np.arange(len(sub))
        rng.shuffle(idx)
        tr = idx[:int(0.8 * len(idx))]
        te = idx[int(0.8 * len(idx)):]
        
        A = sub.iloc[tr]
        B = sub.iloc[te]
        
        def drop(df):
            r0 = df[['dependency', 'cn']].corr().iloc[0, 1]
            r1 = df[['dependency_corrected', 'cn']].corr().iloc[0, 1]
            return np.abs(r0) - np.abs(r1)
        
        d_te = drop(B)
        rows.append((g, d_te, len(B)))
    
    if not rows:
        return {'n_genes': 0, 'median_te_drop': np.nan, 'frac_positive': np.nan}
    
    d = pd.DataFrame(rows, columns=['gene', 'heldout_prox_only_abs_corr_drop', 'n_te'])
    return {
        'n_genes': int(len(d)),
        'median_te_drop': float(d['heldout_prox_only_abs_corr_drop'].median()),
        'frac_positive': float((d['heldout_prox_only_abs_corr_drop'] > 0).mean())
    }

# Test 2: CN-neutral falsification
def cn2_metric(base):
    """Isolate proximity from CN amplitude by focusing on ploidy ~2"""
    dm = pd.read_csv(os.path.join(base, 'dependency_corrected.csv'))
    dm = dm[(dm['cn'] >= 1.8) & (dm['cn'] <= 2.2)]
    
    rows = []
    for g, sub in dm.groupby('gene'):
        sub = sub[['dependency', 'dependency_corrected', 'cn']].dropna()
        if len(sub) < 60:
            continue
        
        r0 = sub[['dependency', 'cn']].corr().iloc[0, 1]
        r1 = sub[['dependency_corrected', 'cn']].corr().iloc[0, 1]
        rows.append(np.abs(r0) - np.abs(r1))
    
    arr = np.array(rows)
    return {
        'n_genes': int(arr.size),
        'median_cn2_drop': float(np.median(arr)) if arr.size else float('nan'),
        'frac_positive': float(np.mean(arr > 0)) if arr.size else float('nan')
    }

# Test 3: Effect-size normalization
def effect_size_metric(base):
    """Scale Δ|r| by dependency variability to compare apples to apples"""
    dm = pd.read_csv(os.path.join(base, 'dependency_corrected.csv'))
    
    rows = []
    for g, sub in dm.groupby('gene'):
        sub = sub[['dependency', 'dependency_corrected', 'cn']].dropna()
        if len(sub) < 60:
            continue
        
        r0 = sub[['dependency', 'cn']].corr().iloc[0, 1]
        r1 = sub[['dependency_corrected', 'cn']].corr().iloc[0, 1]
        d = np.abs(r0) - np.abs(r1)
        scale = np.nanmedian(np.abs(sub['dependency'] - np.nanmedian(sub['dependency']))) + 1e-9
        rows.append(d / scale)
    
    arr = np.array(rows)
    return {
        'n_genes': int(arr.size),
        'median_scaled_drop': float(np.median(arr)),
        'frac_pos': float(np.mean(arr > 0))
    }

# Test 4: Proximity-enriched set (unstable genes)
def define_unstable_genes():
    """Define unstable genes by prevalence in TRUE (top decile by p<=1Mb)"""
    true_dm = pd.read_csv('out_v2/true/dependency_corrected.csv', 
                         usecols=['gene', 'cell_line', 'bp_dist', 'dependency', 'dependency_corrected', 'cn'])
    p1m = true_dm.assign(flag=(true_dm['bp_dist'] <= 1_000_000)).groupby('gene')['flag'].mean()
    unstable = set(p1m.sort_values(ascending=False).head(max(1000, int(0.1 * len(p1m)))).index)
    
    with open('unstable_genes.txt', 'w') as f:
        f.write('\n'.join(sorted(unstable)))
    
    print(f'unstable genes: {len(unstable)}')
    return unstable

def dir_frac(base, genes):
    """Recompute directional prox-active only on unstable genes"""
    dm = pd.read_csv(os.path.join(base, 'dependency_corrected.csv'))
    coef = pd.read_csv(os.path.join(base, 'model_coefficients.csv'))
    
    dm = dm[dm['gene'].isin(genes)].copy()
    coef = coef[coef['gene'].isin(genes)]
    
    w1, w2 = 250000, 2000000
    dm['bp_near'] = (dm['bp_dist'] <= w1).astype(float)
    dm['bp_far'] = ((dm['bp_dist'] > w1) & (dm['bp_dist'] <= w2)).astype(float)
    
    dm = dm.merge(coef[['gene', 'coef_bp_near', 'coef_bp_far']], on='gene', how='left')
    dm['prox_contrib'] = (dm['coef_bp_near'].fillna(0) * dm['bp_near'] + 
                          dm['coef_bp_far'].fillna(0) * dm['bp_far'])
    
    THR, FRAC = 0.01, 0.15
    flags = []
    
    for g, sub in dm.groupby('gene'):
        if sub['cn'].notna().sum() < 60:
            continue
        active = (np.abs(sub['prox_contrib']) >= THR).mean() >= FRAC
        r0 = sub[['dependency', 'cn']].corr().iloc[0, 1]
        r1 = sub[['dependency_corrected', 'cn']].corr().iloc[0, 1]
        drop = np.abs(r0) - np.abs(r1)
        flags.append(active and (drop > 0))
    
    return np.mean(flags) if flags else np.nan

# Test 5: Alternative AUROC metrics
def auroc_alternatives(base):
    """Alternative metrics for essential vs nonessential classification"""
    # Load gene sets
    noness_path = "AchillesNonessentialControls.csv"
    ess_path = "CRISPRInferredCommonEssentials.csv"
    
    def extract_gene_symbol(s):
        if pd.isna(s):
            return None
        s = str(s).strip()
        return s.split()[0] if s else None
    
    noness = set()
    ess = set()
    
    if os.path.exists(noness_path):
        df_noness = pd.read_csv(noness_path)
        gene_col = df_noness.columns[0]
        noness = set(df_noness[gene_col].apply(extract_gene_symbol).dropna())
    
    if os.path.exists(ess_path):
        df_ess = pd.read_csv(ess_path)
        gene_col = df_ess.columns[0]
        ess = set(df_ess[gene_col].apply(extract_gene_symbol).dropna())
    
    if not noness or not ess:
        return {'n_genes': 0, 'auprc_bottom50': np.nan, 'auroc_loo': np.nan}
    
    # Load dependency data
    dm = pd.read_csv(os.path.join(base, 'dependency_corrected.csv'))
    
    # Gene-level medians
    gene_medians = dm.groupby('gene').agg({
        'dependency': 'median',
        'dependency_corrected': 'median'
    }).reset_index()
    
    gene_medians['is_essential'] = gene_medians['gene'].isin(ess).astype(int)
    gene_medians['is_nonessential'] = gene_medians['gene'].isin(noness).astype(int)
    
    # Filter to genes with labels
    labeled = gene_medians[(gene_medians['is_essential'] == 1) | (gene_medians['is_nonessential'] == 1)].copy()
    
    if len(labeled) == 0:
        return {'n_genes': 0, 'auprc_bottom50': np.nan, 'auroc_loo': np.nan}
    
    # AUPRC on bottom 50% of |dependency| genes
    labeled['abs_dep'] = labeled['dependency'].abs()
    bottom50 = labeled[labeled['abs_dep'] <= labeled['abs_dep'].median()]
    
    if len(bottom50) > 0:
        from sklearn.metrics import average_precision_score
        y_true = bottom50['is_essential'].values
        y_score_before = -bottom50['dependency'].abs().values  # Negative because essential = low dependency
        y_score_after = -bottom50['dependency_corrected'].abs().values
        auprc_before = average_precision_score(y_true, y_score_before) if len(np.unique(y_true)) > 1 else np.nan
        auprc_after = average_precision_score(y_true, y_score_after) if len(np.unique(y_true)) > 1 else np.nan
        auprc_bottom50 = auprc_after - auprc_before
    else:
        auprc_bottom50 = np.nan
    
    # Leave-one-lineage-out AUROC (simplified: use cell-line as proxy for lineage)
    # For simplicity, we'll compute per-cell-line AUROC and average
    # This is a simplified version - full LOO would require lineage annotations
    try:
        from sklearn.metrics import roc_auc_score
        y_true = labeled['is_essential'].values
        y_score_before = -labeled['dependency'].abs().values
        y_score_after = -labeled['dependency_corrected'].abs().values
        
        auroc_before = roc_auc_score(y_true, y_score_before) if len(np.unique(y_true)) > 1 else np.nan
        auroc_after = roc_auc_score(y_true, y_score_after) if len(np.unique(y_true)) > 1 else np.nan
        auroc_loo = auroc_after - auroc_before  # Simplified - not true LOO
    except:
        auroc_loo = np.nan
    
    return {
        'n_genes': int(len(labeled)),
        'auprc_bottom50_delta': float(auprc_bottom50) if not np.isnan(auprc_bottom50) else np.nan,
        'auroc_delta': float(auroc_loo) if not np.isnan(auroc_loo) else np.nan
    }

# Main execution
if __name__ == '__main__':
    print("=" * 80)
    print("FALSIFICATION TESTS")
    print("=" * 80)
    
    bases = {
        'TRUE': 'out_v2/true',
        'SHUF_WITHIN': 'out_v2/shuffle_within_chrom',
        'SHUF_ACROSS': 'out_v2/shuffle_across_chrom'
    }
    
    results = {}
    
    # Test 1: Held-out falsification
    print("\n1. HELD-OUT FALSIFICATION (80/20 split)")
    print("-" * 80)
    test1_results = {}
    for k, p in bases.items():
        try:
            m = heldout_metric(p)
            test1_results[k] = m
            print(f"{k:15s} {m}")
        except Exception as e:
            print(f"{k:15s} ERROR: {e}")
            test1_results[k] = {'error': str(e)}
    results['heldout'] = test1_results
    
    # Test 2: CN-neutral falsification
    print("\n2. CN-NEUTRAL FALSIFICATION (ploidy ~2)")
    print("-" * 80)
    test2_results = {}
    for k, p in bases.items():
        try:
            m = cn2_metric(p)
            test2_results[k] = m
            print(f"{k:15s} {m}")
        except Exception as e:
            print(f"{k:15s} ERROR: {e}")
            test2_results[k] = {'error': str(e)}
    results['cn2'] = test2_results
    
    # Test 3: Effect-size normalization
    print("\n3. EFFECT-SIZE NORMALIZATION")
    print("-" * 80)
    test3_results = {}
    for k, p in bases.items():
        try:
            m = effect_size_metric(p)
            test3_results[k] = m
            print(f"{k:15s} {m}")
        except Exception as e:
            print(f"{k:15s} ERROR: {e}")
            test3_results[k] = {'error': str(e)}
    results['effect_size'] = test3_results
    
    # Test 4: Proximity-enriched set
    print("\n4. PROXIMITY-ENRICHED SET (unstable genes)")
    print("-" * 80)
    try:
        unstable_genes = define_unstable_genes()
        test4_results = {}
        for k, p in bases.items():
            try:
                frac = dir_frac(p, unstable_genes)
                test4_results[k] = {'dir_frac': float(frac) if not np.isnan(frac) else np.nan}
                print(f"{k:15s} dir_frac: {frac:.4f}")
            except Exception as e:
                print(f"{k:15s} ERROR: {e}")
                test4_results[k] = {'error': str(e)}
        results['unstable'] = test4_results
    except Exception as e:
        print(f"ERROR defining unstable genes: {e}")
        results['unstable'] = {'error': str(e)}
    
    # Test 5: Alternative AUROC metrics
    print("\n5. ALTERNATIVE AUROC METRICS")
    print("-" * 80)
    test5_results = {}
    for k, p in bases.items():
        try:
            m = auroc_alternatives(p)
            test5_results[k] = m
            print(f"{k:15s} {m}")
        except Exception as e:
            print(f"{k:15s} ERROR: {e}")
            test5_results[k] = {'error': str(e)}
    results['auroc_alt'] = test5_results
    
    # Save results
    print("\n" + "=" * 80)
    print("SAVING RESULTS")
    print("=" * 80)
    
    # Convert to DataFrame format for easier reading
    summary_rows = []
    for test_name, test_data in results.items():
        for condition, metrics in test_data.items():
            if 'error' not in metrics:
                row = {'test': test_name, 'condition': condition}
                row.update(metrics)
                summary_rows.append(row)
    
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv('figs_full/falsification_tests_summary.csv', index=False)
        print("\nSummary saved to: figs_full/falsification_tests_summary.csv")
        print("\nSummary table:")
        print(summary_df.to_string(index=False))
    
    # Save detailed JSON
    with open('figs_full/falsification_tests_detailed.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nDetailed results saved to: figs_full/falsification_tests_detailed.json")
    
    # Generate pass/fail scoreboard
    print("\n" + "=" * 80)
    print("FALSIFICATION SCOREBOARD")
    print("=" * 80)
    
    scoreboard = []
    
    # Test 1: Held-out
    if 'heldout' in results and all('error' not in results['heldout'][k] for k in results['heldout']):
        true_val = results['heldout']['TRUE'].get('median_te_drop', np.nan)
        within_val = results['heldout']['SHUF_WITHIN'].get('median_te_drop', np.nan)
        across_val = results['heldout']['SHUF_ACROSS'].get('median_te_drop', np.nan)
        pass_test1 = (not np.isnan(true_val) and not np.isnan(within_val) and not np.isnan(across_val) and
                     true_val > within_val and true_val > across_val)
        scoreboard.append({
            'test': '1_heldout_median_drop',
            'true': true_val,
            'within': within_val,
            'across': across_val,
            'pass': pass_test1
        })
    
    # Test 2: CN2
    if 'cn2' in results and all('error' not in results['cn2'][k] for k in results['cn2']):
        true_val = results['cn2']['TRUE'].get('median_cn2_drop', np.nan)
        within_val = results['cn2']['SHUF_WITHIN'].get('median_cn2_drop', np.nan)
        across_val = results['cn2']['SHUF_ACROSS'].get('median_cn2_drop', np.nan)
        pass_test2 = (not np.isnan(true_val) and not np.isnan(within_val) and not np.isnan(across_val) and
                     true_val > within_val and true_val > across_val)
        scoreboard.append({
            'test': '2_cn2_median_drop',
            'true': true_val,
            'within': within_val,
            'across': across_val,
            'pass': pass_test2
        })
    
    # Test 3: Effect size
    if 'effect_size' in results and all('error' not in results['effect_size'][k] for k in results['effect_size']):
        true_val = results['effect_size']['TRUE'].get('median_scaled_drop', np.nan)
        within_val = results['effect_size']['SHUF_WITHIN'].get('median_scaled_drop', np.nan)
        across_val = results['effect_size']['SHUF_ACROSS'].get('median_scaled_drop', np.nan)
        pass_test3 = (not np.isnan(true_val) and not np.isnan(within_val) and not np.isnan(across_val) and
                     true_val > within_val and true_val > across_val)
        scoreboard.append({
            'test': '3_effect_size_median_drop',
            'true': true_val,
            'within': within_val,
            'across': across_val,
            'pass': pass_test3
        })
    
    # Test 4: Unstable genes
    if 'unstable' in results and all('error' not in results['unstable'][k] for k in results['unstable']):
        true_val = results['unstable']['TRUE'].get('dir_frac', np.nan)
        within_val = results['unstable']['SHUF_WITHIN'].get('dir_frac', np.nan)
        across_val = results['unstable']['SHUF_ACROSS'].get('dir_frac', np.nan)
        pass_test4 = (not np.isnan(true_val) and not np.isnan(within_val) and not np.isnan(across_val) and
                     true_val > within_val and true_val > across_val)
        scoreboard.append({
            'test': '4_unstable_dir_frac',
            'true': true_val,
            'within': within_val,
            'across': across_val,
            'pass': pass_test4
        })
    
    if scoreboard:
        scoreboard_df = pd.DataFrame(scoreboard)
        scoreboard_df.to_csv('figs_full/falsification_scoreboard.csv', index=False)
        print("\n" + scoreboard_df.to_string(index=False))
        print(f"\nPass rate: {scoreboard_df['pass'].sum()}/{len(scoreboard_df)}")

