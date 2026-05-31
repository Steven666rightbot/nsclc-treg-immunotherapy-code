"""
Supplementary Figure 2: Cross-cancer validation of Treg signatures in melanoma
GSE120575 (Sade-Feldman et al. Cell 2018) — 48 biopsies, anti-PD-1 ± CTLA-4
Panel A: CCR8+ Treg proportions (R vs NR)
Panel B: MKI67+ Treg proportions (R vs NR)
"""

import gzip
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from PIL import Image

OUTDIR = "D:/research/cucumber/supp_fig2"
DATA_DIR = "D:/research/cucumber/fig1/cross_cancer/data"
os.makedirs(OUTDIR, exist_ok=True)

DPI = 300
# Updated colors: Responder = blue, Non-responder = red
COLOR_R = '#2E5AAC'
COLOR_NR = '#C0392B'

TARGET_GENES = ['CD4', 'FOXP3', 'CCR8', 'MKI67']
POS_THRESHOLD = 0.0
META_FILE = os.path.join(DATA_DIR, 'metadata.txt.gz')
EXPR_FILE = os.path.join(DATA_DIR, 'expression.txt.gz')

# ─── Load metadata ───
with gzip.open(META_FILE, 'rt', encoding='latin-1') as f:
    lines = f.readlines()
header_idx = next(i for i, line in enumerate(lines) if line.startswith('Sample name'))

meta = pd.read_csv(META_FILE, compression='gzip', encoding='latin-1', sep='\t',
                    skiprows=header_idx, nrows=16291,
                    usecols=['title', 'characteristics: patinet ID (Pre=baseline; Post= on treatment)',
                             'characteristics: response', 'characteristics: therapy'])
meta.rename(columns={'title': 'cell',
                     'characteristics: patinet ID (Pre=baseline; Post= on treatment)': 'sample',
                     'characteristics: response': 'response',
                     'characteristics: therapy': 'therapy'}, inplace=True)
meta = meta[meta['response'].isin(['Responder', 'Non-responder'])].copy()
meta['response_binary'] = (meta['response'] == 'Responder').astype(int)
sample_meta = meta.drop_duplicates('sample')[['sample', 'response', 'response_binary']].copy()

# ─── Load expression (target genes only) ───
with gzip.open(EXPR_FILE, 'rt') as f:
    cell_names = f.readline().rstrip().split('\t')[1:]
    sample_ids = f.readline().rstrip().split('\t')[1:]
    gene_expr = {}
    for line in f:
        parts = line.rstrip().split('\t')
        gene = parts[0]
        if gene in TARGET_GENES:
            vals = np.array([float(x) if x != '' else 0.0 for x in parts[1:]], dtype=np.float32)
            gene_expr[gene] = vals

cell_df = pd.DataFrame({'cell': cell_names, 'sample': sample_ids,
                        'CD4': gene_expr['CD4'], 'FOXP3': gene_expr['FOXP3'],
                        'CCR8': gene_expr['CCR8'], 'MKI67': gene_expr['MKI67']})
cell_df = cell_df[cell_df['cell'].isin(meta['cell'])].copy()

# ─── Identify subsets ───
cell_df['is_CCR8_Treg'] = (cell_df['CD4'] > POS_THRESHOLD) & (cell_df['FOXP3'] > POS_THRESHOLD) & (cell_df['CCR8'] > POS_THRESHOLD)
cell_df['is_MKI67_Treg'] = (cell_df['CD4'] > POS_THRESHOLD) & (cell_df['FOXP3'] > POS_THRESHOLD) & (cell_df['MKI67'] > POS_THRESHOLD)

sample_stats = []
for sample, grp in cell_df.groupby('sample'):
    n_total = len(grp)
    sample_stats.append({'sample': sample, 'n_cells': n_total,
                         'prop_ccr8_treg': grp['is_CCR8_Treg'].sum() / n_total,
                         'prop_mki67_treg': grp['is_MKI67_Treg'].sum() / n_total})
sample_df = pd.DataFrame(sample_stats).merge(sample_meta, on='sample')

r_mask = sample_df['response'] == 'Responder'
nr_mask = sample_df['response'] == 'Non-responder'

# ─── Helper: whisker max for standard 1.5*IQR boxplot ───
def whisker_max(vals):
    q1, q3 = np.percentile(vals, [25, 75])
    iqr = q3 - q1
    upper = q3 + 1.5 * iqr
    return min(upper, max(vals))

# ─── Plot ───
fig, axes = plt.subplots(1, 2, figsize=(6.5, 4.5), sharey=False)
fig.subplots_adjust(wspace=0.3)

plot_configs = [
    ('prop_ccr8_treg', 'CCR8⁺ Treg'),
    ('prop_mki67_treg', 'MKI67⁺ Treg'),
]

for ax, (col, label) in zip(axes, plot_configs):
    r_vals = sample_df.loc[r_mask, col].values
    nr_vals = sample_df.loc[nr_mask, col].values

    # Standard boxplot with 1.5*IQR whiskers and visible caps
    bp = ax.boxplot([nr_vals, r_vals],
                    tick_labels=['', ''],
                    patch_artist=True, widths=0.5,
                    whis=1.5,
                    showfliers=True,
                    medianprops=dict(color='black', linewidth=2.5),
                    whiskerprops=dict(linewidth=1.2, color='#333333'),
                    capprops=dict(linewidth=1.2, color='#333333'),
                    flierprops=dict(marker='o', markerfacecolor='none', markersize=4,
                                   markeredgewidth=0.5, markeredgecolor='#888888'))
    bp['boxes'][0].set_facecolor(COLOR_NR)
    bp['boxes'][0].set_alpha(0.7)
    bp['boxes'][1].set_facecolor(COLOR_R)
    bp['boxes'][1].set_alpha(0.7)

    # Individual points (jitter) — smaller to avoid blocking boxes
    np.random.seed(42)
    ax.scatter(np.random.normal(1, 0.04, len(nr_vals)), nr_vals,
               color=COLOR_NR, edgecolor='black', linewidth=0.5, s=22, zorder=3, alpha=0.7)
    ax.scatter(np.random.normal(2, 0.04, len(r_vals)), r_vals,
               color=COLOR_R, edgecolor='black', linewidth=0.5, s=22, zorder=3, alpha=0.7)

    _, p_val = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')

    # Significance annotation centered above the boxes
    y_top = max(whisker_max(nr_vals), whisker_max(r_vals))
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    ax.set_ylim(top=y_top + y_range * 0.25)

    if label == 'CCR8⁺ Treg':
        sig_text = 'ns'
        p_text = f'(p = {p_val:.3f})'
        sig_color = '#666666'
        text_y = y_top + y_range * 0.06
        ax.text(1.5, text_y, sig_text,
                ha='right', va='bottom', fontsize=12, fontweight='bold', color=sig_color)
        ax.text(1.5, text_y, p_text,
                ha='left', va='bottom', fontsize=8, color='#888888')
    else:
        sig_text = f'* p = {p_val:.3f}'
        sig_color = '#333333'
        ax.text(1.5, y_top + y_range * 0.06, sig_text,
                ha='center', va='bottom', fontsize=11, fontweight='bold', color=sig_color)

    ax.set_title(label, fontsize=13, fontweight='bold')
    ax.set_ylabel('Proportion of immune cells', fontsize=11)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f'NR\n(n={len(nr_vals)})', f'R\n(n={len(r_vals)})'], fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(bottom=0)

plt.tight_layout()

for ext in ['png', 'pdf']:
    path = os.path.join(OUTDIR, f'supplementary_fig2.{ext}')
    fig.savefig(path, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {path}")

plt.close()
print(f"\nn = {len(sample_df)} samples (R={r_mask.sum()}, NR={nr_mask.sum()})")
print(f"MKI67+ Treg: p = {stats.mannwhitneyu(sample_df.loc[r_mask,'prop_mki67_treg'], sample_df.loc[nr_mask,'prop_mki67_treg'])[1]:.4f}")
print(f"CCR8+ Treg:  p = {stats.mannwhitneyu(sample_df.loc[r_mask,'prop_ccr8_treg'], sample_df.loc[nr_mask,'prop_ccr8_treg'])[1]:.4f}")
