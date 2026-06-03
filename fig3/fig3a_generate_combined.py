#!/usr/bin/env python3
"""
Generate merged Figure 3A: CCR8+ and MKI67+ TF bar plots side by side.
Single (A) label, clean publication layout.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import os, warnings
warnings.filterwarnings('ignore')

OUT_DIR = 'D:/research/cucumber/fig3'
DATA_DIR = 'D:/research/tomato'
TF_ACT_PATH = os.path.join(DATA_DIR, 'results/scenic_input/dorothea_tf_activity.csv')
META_PATH = os.path.join(DATA_DIR, 'results/scenic_input/treg_no_foxp3_metadata.csv')

COLOR_HIGHER_R = '#2E5AAC'
COLOR_HIGHER_NR = '#C0392B'

# ── Load ──
tf_df = pd.read_csv(TF_ACT_PATH, index_col=0)
meta = pd.read_csv(META_PATH)
meta = meta[meta['CellID'].isin(tf_df.index)]
data = tf_df.loc[meta['CellID']].copy()
data['sub_cell_type'] = meta['sub_cell_type'].values
data['response'] = meta['response'].values

def diff_analysis(data, subtype):
    sub = data[data['sub_cell_type'] == subtype]
    tf_names = [c for c in sub.columns if c not in ('sub_cell_type', 'response')]
    results = []
    for tf in tf_names:
        rv = sub.loc[sub['response'] == 'Responder', tf].dropna()
        nrv = sub.loc[sub['response'] == 'Non-responder', tf].dropna()
        if len(rv) < 3 or len(nrv) < 3:
            continue
        _, p = mannwhitneyu(rv, nrv, alternative='two-sided')
        results.append({'TF': tf, 'effect_size': rv.median() - nrv.median(),
                        'p_value': p, 'R_median': rv.median(), 'NR_median': nrv.median()})
    res = pd.DataFrame(results)
    _, padj, _, _ = multipletests(res['p_value'], method='fdr_bh')
    res['padj'] = padj
    return res.sort_values('effect_size', ascending=False).reset_index(drop=True)

def star(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    return 'ns'

def make_subpanel(ax, res, title, show_ylabel=True):
    top = res.nlargest(10, 'effect_size').iloc[::-1].copy()
    top['star'] = top['padj'].apply(star)
    top['label'] = top['TF'] + '  ' + top['star']
    colors = [COLOR_HIGHER_R if e > 0 else COLOR_HIGHER_NR for e in top['effect_size']]
    
    ax.barh(range(len(top)), top['effect_size'].values, color=colors,
            edgecolor='white', linewidth=0.3, height=0.7)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top['label'].values, fontsize=8)
    ax.axvline(0, color='grey', linewidth=0.6, alpha=0.5)
    ax.set_title(title, fontsize=10, fontweight='bold')
    
    for i, (_, row) in enumerate(top.iterrows()):
        val = row['effect_size']
        st = row['star']
        if st != 'ns':
            xr = top['effect_size'].max() - top['effect_size'].min()
            offset = 0.02 * xr
            c = COLOR_HIGHER_R if val > 0 else COLOR_HIGHER_NR
            ha = 'left' if val > 0 else 'right'
            ax.text(val + (offset if val > 0 else -offset), i, st,
                    va='center', fontsize=7, color=c, fontweight='bold', ha=ha)
    
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    
    xlim = max(abs(top['effect_size'].min()), abs(top['effect_size'].max()))
    ax.set_xlim(-xlim * 1.25, xlim * 1.25)
    
    if show_ylabel:
        ax.set_ylabel('Transcription Factor', fontsize=9)
    else:
        ax.set_yticklabels([])
        ax.set_ylabel('')
    
    return top

# ── Compute ──
res_ccr8 = diff_analysis(data, 'CD4T_Treg_CCR8')
res_mki67 = diff_analysis(data, 'CD4T_Treg_MKI67')

# ── Figure ──
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), gridspec_kw={'wspace': 0.3})

top_a = make_subpanel(axes[0], res_ccr8, 'CCR8+ Treg', show_ylabel=True)
top_b = make_subpanel(axes[1], res_mki67, 'MKI67+ Treg', show_ylabel=False)

# Shared x-axis label
fig.text(0.5, 0.01, 'Effect size (Responder − Non-responder)', ha='center', fontsize=9)

# Single (A) label
fig.text(0.01, 0.98, '(A)', fontsize=14, fontweight='bold', transform=fig.transFigure)

# Shared legend
legend_elements = [
    mpatches.Patch(facecolor=COLOR_HIGHER_R, edgecolor='white', label='Higher in Responder'),
    mpatches.Patch(facecolor=COLOR_HIGHER_NR, edgecolor='white', label='Higher in Non-responder'),
]
fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 1.02),
           ncol=2, fontsize=8, frameon=False)

fig.tight_layout(rect=[0, 0.03, 1, 0.97])
fig.savefig(os.path.join(OUT_DIR, 'fig3a_combined_bar.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'fig3a_combined_bar.pdf'), dpi=300, bbox_inches='tight')
plt.close(fig)

# Save results
combined = pd.concat([
    res_ccr8.assign(subtype='CCR8'),
    res_mki67.assign(subtype='MKI67')
], ignore_index=True)
combined.to_csv(os.path.join(OUT_DIR, 'fig3a_tf_differential_results.csv'), index=False)

print(f"Saved: fig3a_combined_bar.png/pdf")
print(f"  CCR8 sig: {res_ccr8['padj'].lt(0.05).sum()} / {len(res_ccr8)} TFs")
print(f"  MKI67 sig: {res_mki67['padj'].lt(0.05).sum()} / {len(res_mki67)} TFs")
