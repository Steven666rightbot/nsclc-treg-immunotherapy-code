#!/usr/bin/env python3
"""
Cross-cancer validation: GSE120575 melanoma scRNA-seq (Sade-Feldman Cell 2018)

Validate whether CCR8+ Treg and MKI67+ Treg signatures distinguish
immunotherapy responders from non-responders across cancer types.

Strategy B — Quick signature validation:
  1. Load GSE120575 processed TPM data + metadata
  2. Identify CD4+ T cells, then CCR8+FOXP3+ and MKI67+FOXP3+ Treg subsets
  3. Compute per-sample proportions among immune cells
  4. Compare responder (R) vs non-responder (NR) with Mann-Whitney U
  5. Export boxplot (PNG/PDF) and print key statistics

Data source:
  GEO GSE120575 — 48 tumor biopsies from 32 melanoma patients,
  anti-PD-1 ± CTLA-4, CD45+ Smart-seq2 scRNA-seq.

Output:
  cross_cancer_validation.png / .pdf
"""

import os
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# ─── Configuration ──────────────────────────────────────────────────────────
DATA_DIR = r'D:\research\cucumber\fig1\cross_cancer\data'
OUT_DIR  = r'D:\research\cucumber\fig1\cross_cancer'
os.makedirs(OUT_DIR, exist_ok=True)

EXPR_FILE = os.path.join(DATA_DIR, 'expression.txt.gz')
META_FILE = os.path.join(DATA_DIR, 'metadata.txt.gz')

TARGET_GENES = ['CD4', 'FOXP3', 'CCR8', 'MKI67']
POS_THRESHOLD = 0.0          # TPM > 0 considered positive (consistent with GSE207422)

# Figure style (match Figure 1)
COLOR_R  = '#8bc98b'         # soft green — Responder
COLOR_NR = '#b08ebf'         # soft purple — Non-responder
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 12,
    'axes.facecolor': 'white',
    'figure.facecolor': 'white',
})

# ─── 1. Load metadata ───────────────────────────────────────────────────────
print("=" * 60)
print("GSE120575 Cross-Cancer Validation (Melanoma)")
print("=" * 60)

with gzip.open(META_FILE, 'rt', encoding='latin-1') as f:
    lines = f.readlines()

# Find header row
header_idx = None
for i, line in enumerate(lines):
    if line.startswith('Sample name'):
        header_idx = i
        break

meta = pd.read_csv(
    os.path.join(DATA_DIR, 'metadata.txt.gz'),
    compression='gzip',
    encoding='latin-1',
    sep='\t',
    skiprows=header_idx,
    nrows=16291,              # exactly the number of cells
    usecols=['title',
             'characteristics: patinet ID (Pre=baseline; Post= on treatment)',
             'characteristics: response',
             'characteristics: therapy'],
)
meta.rename(columns={
    'title': 'cell',
    'characteristics: patinet ID (Pre=baseline; Post= on treatment)': 'sample',
    'characteristics: response': 'response',
    'characteristics: therapy': 'therapy',
}, inplace=True)

# Drop any malformed rows
meta = meta[meta['response'].isin(['Responder', 'Non-responder'])].copy()
meta['response_binary'] = (meta['response'] == 'Responder').astype(int)

print(f"Metadata loaded: {len(meta)} cells")
print(f"  Responders    : {(meta['response']=='Responder').sum()}")
print(f"  Non-responders: {(meta['response']=='Non-responder').sum()}")
print(f"  Unique samples: {meta['sample'].nunique()}")
print(f"  Therapy breakdown:\n{meta['therapy'].value_counts()}")

# Sample-level response map (each sample has a single response)
sample_meta = meta.drop_duplicates('sample')[['sample', 'response', 'response_binary', 'therapy']].copy()

# ─── 2. Load expression matrix (target genes only) ──────────────────────────
print("\nLoading expression matrix (target genes only)...")

with gzip.open(EXPR_FILE, 'rt') as f:
    cell_names = f.readline().rstrip().split('\t')[1:]   # skip leading empty field
    sample_ids = f.readline().rstrip().split('\t')[1:]

    gene_expr = {}
    for line in f:
        parts = line.rstrip().split('\t')
        gene = parts[0]
        if gene in TARGET_GENES:
            # parse numeric values; empty strings → 0.0
            vals = np.array([float(x) if x != '' else 0.0 for x in parts[1:]], dtype=np.float32)
            gene_expr[gene] = vals
            if len(gene_expr) == len(TARGET_GENES):
                break

# Build cell-level DataFrame
cell_df = pd.DataFrame({
    'cell': cell_names,
    'sample': sample_ids,
    'CD4': gene_expr['CD4'],
    'FOXP3': gene_expr['FOXP3'],
    'CCR8': gene_expr['CCR8'],
    'MKI67': gene_expr['MKI67'],
})

# Restrict to cells present in metadata
cell_df = cell_df[cell_df['cell'].isin(meta['cell'])].copy()
print(f"Expression loaded: {len(cell_df)} cells × {len(TARGET_GENES)} genes")

# ─── 3. Identify cell subsets ───────────────────────────────────────────────
print("\nIdentifying cell subsets...")

# Boolean masks (positive = expression > threshold)
cell_df['is_CD4']    = cell_df['CD4'] > POS_THRESHOLD
cell_df['is_FOXP3']  = cell_df['FOXP3'] > POS_THRESHOLD
cell_df['is_CCR8']   = cell_df['CCR8'] > POS_THRESHOLD
cell_df['is_MKI67']  = cell_df['MKI67'] > POS_THRESHOLD

# Define subsets as requested:
#   CCR8+ Treg = CD4+ & FOXP3+ & CCR8+
#   MKI67+ Treg = CD4+ & FOXP3+ & MKI67+
cell_df['is_CCR8_Treg']  = cell_df['is_CD4'] & cell_df['is_FOXP3'] & cell_df['is_CCR8']
cell_df['is_MKI67_Treg'] = cell_df['is_CD4'] & cell_df['is_FOXP3'] & cell_df['is_MKI67']

# Also compute plain Treg (CD4+ FOXP3+) for reference
cell_df['is_Treg'] = cell_df['is_CD4'] & cell_df['is_FOXP3']

print(f"  CD4+ cells      : {cell_df['is_CD4'].sum():,}")
print(f"  Treg (CD4+FOXP3+): {cell_df['is_Treg'].sum():,}")
print(f"  CCR8+ Treg      : {cell_df['is_CCR8_Treg'].sum():,}")
print(f"  MKI67+ Treg     : {cell_df['is_MKI67_Treg'].sum():,}")

# ─── 4. Aggregate per-sample proportions ────────────────────────────────────
print("\nAggregating per-sample proportions (among immune cells)...")

sample_stats = []
for sample, grp in cell_df.groupby('sample'):
    n_total = len(grp)                      # all cells = immune cells (CD45+ sorted)
    n_ccr8_treg  = grp['is_CCR8_Treg'].sum()
    n_mki67_treg = grp['is_MKI67_Treg'].sum()
    n_treg       = grp['is_Treg'].sum()

    sample_stats.append({
        'sample': sample,
        'n_cells': n_total,
        'n_treg': n_treg,
        'n_ccr8_treg': n_ccr8_treg,
        'n_mki67_treg': n_mki67_treg,
        'prop_ccr8_treg': n_ccr8_treg / n_total if n_total > 0 else 0,
        'prop_mki67_treg': n_mki67_treg / n_total if n_total > 0 else 0,
        'prop_treg': n_treg / n_total if n_total > 0 else 0,
    })

sample_df = pd.DataFrame(sample_stats)

# Merge response labels
sample_df = sample_df.merge(sample_meta[['sample', 'response', 'response_binary', 'therapy']],
                            on='sample', how='left')

# Drop samples with missing response info (should not happen)
sample_df = sample_df.dropna(subset=['response']).copy()

print(f"Samples aggregated: {len(sample_df)}")
print(f"  Responder samples    : {(sample_df['response']=='Responder').sum()}")
print(f"  Non-responder samples: {(sample_df['response']=='Non-responder').sum()}")

# ─── 5. Statistical comparison ──────────────────────────────────────────────
print("\n" + "=" * 60)
print("Statistical Results (Mann-Whitney U)")
print("=" * 60)

r_mask = sample_df['response'] == 'Responder'
nr_mask = sample_df['response'] == 'Non-responder'

for label, col in [('CCR8+ Treg', 'prop_ccr8_treg'),
                   ('MKI67+ Treg', 'prop_mki67_treg')]:
    r_vals = sample_df.loc[r_mask, col]
    nr_vals = sample_df.loc[nr_mask, col]

    # Mann-Whitney U (two-sided)
    u_stat, p_val = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')

    print(f"\n{label}:")
    print(f"  Responder    : n={len(r_vals)}, mean={r_vals.mean():.4f}, median={r_vals.median():.4f}")
    print(f"  Non-responder: n={len(nr_vals)}, mean={nr_vals.mean():.4f}, median={nr_vals.median():.4f}")
    print(f"  Mann-Whitney U = {u_stat:.1f}, p = {p_val:.2e}")

# ─── 6. Visualization ───────────────────────────────────────────────────────
print("\nGenerating figure...")

fig, axes = plt.subplots(1, 2, figsize=(10, 5.5), sharey=False)
fig.subplots_adjust(wspace=0.35)

plot_configs = [
    ('prop_ccr8_treg', 'CCR8+ Treg', 'Proportion of immune cells'),
    ('prop_mki67_treg', 'MKI67+ Treg', 'Proportion of immune cells'),
]

for ax, (col, title, ylabel) in zip(axes, plot_configs):
    r_vals = sample_df.loc[r_mask, col].values
    nr_vals = sample_df.loc[nr_mask, col].values

    # Boxplot
    bp = ax.boxplot([nr_vals, r_vals],
                    tick_labels=[f'NR\nn={len(nr_vals)}', f'R\nn={len(r_vals)}'],
                    patch_artist=True,
                    widths=0.5,
                    medianprops=dict(color='black', linewidth=2),
                    whiskerprops=dict(linewidth=1.2),
                    capprops=dict(linewidth=1.2))

    bp['boxes'][0].set_facecolor(COLOR_NR)
    bp['boxes'][0].set_alpha(0.7)
    bp['boxes'][1].set_facecolor(COLOR_R)
    bp['boxes'][1].set_alpha(0.7)

    # Strip plot (individual samples)
    np.random.seed(42)
    jitter_nr = np.random.normal(1, 0.04, size=len(nr_vals))
    jitter_r  = np.random.normal(2, 0.04, size=len(r_vals))
    ax.scatter(jitter_nr, nr_vals, color=COLOR_NR, edgecolor='black', linewidth=0.5,
               s=45, zorder=3, alpha=0.8)
    ax.scatter(jitter_r,  r_vals,  color=COLOR_R,  edgecolor='black', linewidth=0.5,
               s=45, zorder=3, alpha=0.8)

    # P-value annotation (top-right inside plot)
    _, p_val = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
    if p_val < 1e-10:
        p_text = 'p < 1e-10'
    elif p_val < 0.001:
        p_text = f'p = {p_val:.2e}'
    else:
        p_text = f'p = {p_val:.3f}'

    ax.text(0.97, 0.97, p_text, transform=ax.transAxes,
            ha='right', va='top', fontsize=12, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.9))

    # Styling
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=11)
    ax.set_ylim(bottom=0)

    # Add n annotations
    ylim_min, ylim_max = ax.get_ylim()
    ax.text(1, ylim_min - (ylim_max - ylim_min) * 0.02, f"n={len(nr_vals)}",
            ha='center', va='top', fontsize=10, color='#555555', transform=ax.transData)
    ax.text(2, ylim_min - (ylim_max - ylim_min) * 0.02, f"n={len(r_vals)}",
            ha='center', va='top', fontsize=10, color='#555555', transform=ax.transData)

# Overall title
fig.suptitle('Cross-cancer Validation — GSE120575 Melanoma (n=48 biopsies)',
             fontsize=15, fontweight='bold', y=1.02)

plt.tight_layout()

png_path = os.path.join(OUT_DIR, 'cross_cancer_validation.png')
pdf_path = os.path.join(OUT_DIR, 'cross_cancer_validation.pdf')
fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(pdf_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"\nFigure saved:")
print(f"  {png_path}")
print(f"  {pdf_path}")

# ─── 7. Save per-sample data ────────────────────────────────────────────────
csv_path = os.path.join(OUT_DIR, 'cross_cancer_sample_stats.csv')
sample_df.to_csv(csv_path, index=False)
print(f"\nPer-sample data saved: {csv_path}")

# ─── 7. Pre-treatment only supplementary analysis ───────────────────────────
print("\n" + "=" * 60)
print("Supplementary: Pre-treatment only analysis")
print("=" * 60)

pre_df = sample_df[sample_df['sample'].str.startswith('Pre')].copy()
print(f"Pre-treatment samples: {len(pre_df)}")
print(f"  Responder    : {(pre_df['response']=='Responder').sum()}")
print(f"  Non-responder: {(pre_df['response']=='Non-responder').sum()}")

r_mask_pre = pre_df['response'] == 'Responder'
nr_mask_pre = pre_df['response'] == 'Non-responder'

for label, col in [('CCR8+ Treg', 'prop_ccr8_treg'),
                   ('MKI67+ Treg', 'prop_mki67_treg')]:
    r_vals = pre_df.loc[r_mask_pre, col]
    nr_vals = pre_df.loc[nr_mask_pre, col]
    if len(r_vals) > 0 and len(nr_vals) > 0:
        u_stat, p_val = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
        print(f"\n{label} (Pre only):")
        print(f"  Responder    : n={len(r_vals)}, mean={r_vals.mean():.4f}, median={r_vals.median():.4f}")
        print(f"  Non-responder: n={len(nr_vals)}, mean={nr_vals.mean():.4f}, median={nr_vals.median():.4f}")
        print(f"  Mann-Whitney U = {u_stat:.1f}, p = {p_val:.2e}")

# Plot pre-treatment only
fig2, axes2 = plt.subplots(1, 2, figsize=(10, 5.5), sharey=False)
fig2.subplots_adjust(wspace=0.35)

for ax, (col, title, ylabel) in zip(axes2, plot_configs):
    r_vals = pre_df.loc[r_mask_pre, col].values
    nr_vals = pre_df.loc[nr_mask_pre, col].values

    bp = ax.boxplot([nr_vals, r_vals],
                    tick_labels=[f'NR\nn={len(nr_vals)}', f'R\nn={len(r_vals)}'],
                    patch_artist=True,
                    widths=0.5,
                    medianprops=dict(color='black', linewidth=2),
                    whiskerprops=dict(linewidth=1.2),
                    capprops=dict(linewidth=1.2))

    bp['boxes'][0].set_facecolor(COLOR_NR)
    bp['boxes'][0].set_alpha(0.7)
    bp['boxes'][1].set_facecolor(COLOR_R)
    bp['boxes'][1].set_alpha(0.7)

    np.random.seed(42)
    jitter_nr = np.random.normal(1, 0.04, size=len(nr_vals))
    jitter_r  = np.random.normal(2, 0.04, size=len(r_vals))
    ax.scatter(jitter_nr, nr_vals, color=COLOR_NR, edgecolor='black', linewidth=0.5,
               s=45, zorder=3, alpha=0.8)
    ax.scatter(jitter_r,  r_vals,  color=COLOR_R,  edgecolor='black', linewidth=0.5,
               s=45, zorder=3, alpha=0.8)

    _, p_val = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
    if p_val < 1e-10:
        p_text = 'p < 1e-10'
    elif p_val < 0.001:
        p_text = f'p = {p_val:.2e}'
    else:
        p_text = f'p = {p_val:.3f}'

    ax.text(0.97, 0.97, p_text, transform=ax.transAxes,
            ha='right', va='top', fontsize=12, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.9))

    ax.set_title(title, fontsize=14, fontweight='bold', pad=10)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=11)
    ax.set_ylim(bottom=0)

fig2.suptitle('Cross-cancer Validation — GSE120575 Pre-treatment Melanoma (n=19 biopsies)',
              fontsize=15, fontweight='bold', y=1.02)

plt.tight_layout()
png_path2 = os.path.join(OUT_DIR, 'cross_cancer_validation_pre_only.png')
pdf_path2 = os.path.join(OUT_DIR, 'cross_cancer_validation_pre_only.pdf')
fig2.savefig(png_path2, dpi=300, bbox_inches='tight', facecolor='white')
fig2.savefig(pdf_path2, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"\nPre-treatment figure saved:")
print(f"  {png_path2}")
print(f"  {pdf_path2}")

# ─── 8. Save per-sample data ────────────────────────────────────────────────
csv_path = os.path.join(OUT_DIR, 'cross_cancer_sample_stats.csv')
sample_df.to_csv(csv_path, index=False)
print(f"\nPer-sample data saved: {csv_path}")

print("\n" + "=" * 60)
print("Validation complete.")
print("=" * 60)
