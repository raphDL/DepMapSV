import os, pandas as pd, numpy as np, matplotlib.pyplot as plt

UNSTABLE_DIR = os.environ.get('UNSTABLE_DIR', 'out_pilot_final/unstable')
OUT_DIR = 'figs_pilot'
os.makedirs(OUT_DIR, exist_ok=True)

dm = pd.read_csv(os.path.join(UNSTABLE_DIR, 'design_matrix.csv'))
coef = pd.read_csv(os.path.join(UNSTABLE_DIR, 'models_coefficients.csv'))
coef = coef.rename(columns={
    'cn':'b_cn', 'bp_within_100000':'b_w100k',
    'bp_within_1000000':'b_w1m', 'inv_bp':'b_inv'
})

# choose genes: prefer top5_unstable.csv if present
cands = []
top5_csv = os.path.join('figs_pilot', 'top5_unstable.csv')
if os.path.exists(top5_csv):
    top = pd.read_csv(top5_csv)
    cands = [str(g) for g in top['gene'].head(3).tolist()]
else:
    # fallback to earlier surfaced examples
    cands = [g for g in ['GFPT2','TBC1D22A','TLR3','KRT9','PACS1'] if g in dm['gene'].unique()][:3]

if not cands:
    # final fallback: pick 3 genes with highest |b_inv|+windows
    tmp = coef.copy()
    for c in ['b_w100k','b_w1m','b_inv']:
        if c not in tmp: tmp[c] = 0.0
    tmp['prox_mag'] = tmp[['b_w100k','b_w1m','b_inv']].abs().sum(axis=1)
    cands = tmp.sort_values('prox_mag', ascending=False)['gene'].head(3).astype(str).tolist()

for gene in cands:
    d = dm[dm['gene'].astype(str)==gene].copy()
    cf = coef[coef['gene'].astype(str)==gene][['gene','b_cn','b_w100k','b_w1m','b_inv']].copy()
    if d.empty or cf.empty:
        continue
    b = cf.iloc[0]
    # compute proximity contribution on raw features
    for col in ['bp_within_100000','bp_within_1000000','inv_bp']:
        if col not in d.columns: d[col] = 0.0
    d['prox_contrib'] = (
        float(b.get('b_w100k',0.0))*d['bp_within_100000'] +
        float(b.get('b_w1m',0.0))*d['bp_within_1000000'] +
        float(b.get('b_inv',0.0))*d['inv_bp']
    )
    d['prox_only_corrected'] = d['dependency'] - d['prox_contrib']
    # derive distance from bp_dist (already in bp units)
    dplot = d[['cell_line','dependency','prox_only_corrected','bp_dist','cn']].dropna(subset=['dependency','bp_dist'])

    fig, ax = plt.subplots(1, 2, figsize=(9,3.2), sharey=True)
    ax[0].scatter(dplot['bp_dist'], dplot['dependency'], s=8, alpha=0.6, edgecolors='none')
    ax[0].set_xscale('symlog', linthresh=1e4)
    ax[0].set_xlabel('nearest breakpoint distance (bp)')
    ax[0].set_ylabel('dependency')
    ax[0].set_title(f'{gene}: before')

    ax[1].scatter(dplot['bp_dist'], dplot['prox_only_corrected'], s=8, alpha=0.6, edgecolors='none', color='#2a9d8f')
    ax[1].set_xscale('symlog', linthresh=1e4)
    ax[1].set_xlabel('nearest breakpoint distance (bp)')
    ax[1].set_title('after prox-only')

    # annotate proximity coefficients
    txt = f"b_w100k={float(b.get('b_w100k',0.0)):.3f}\n" \
          f"b_w1m={float(b.get('b_w1m',0.0)):.3f}\n" \
          f"b_inv={float(b.get('b_inv',0.0)):.3f}"
    fig.text(0.99, 0.05, txt, ha='right', va='bottom', fontsize=9, family='monospace')
    fig.tight_layout()
    outpath = os.path.join(OUT_DIR, f'case_{gene}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'Wrote {outpath}')
