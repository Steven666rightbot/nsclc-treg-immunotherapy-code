#!/usr/bin/env python3
"""
Figure 6E (New) — Pseudotime Key Gene Expression Heatmap

Shows smooth expression of 22 key genes along pseudotime, split by response.
Rows = genes ordered by peak expression time.
Columns = pseudotime bins.
Color = mean log10(CPM+1) expression.

Data source: fig6/data/slingshot_gene_expr.csv (26,725 cells, 22 genes)

Output: fig6/fig6e_pseudotime_heatmap.{png,pdf}
"""

import os, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.ndimage import gaussian_filter1d

warnings.filterwarnings('ignore')

# ─── Paths ───────────────────────────────────────────────────────────────────
DATA_FILE = r'D:/research/cucumber/fig6/data/slingshot_gene_expr.csv'
OUT_DIR = r'D:/research/cucumber/fig6'
os.makedirs(OUT_DIR, exist_ok=True)

# ─── Publication style ──────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 7,
    'font.sans-serif': ['Arial'],
    'axes.labelsize': 7.5,
    'axes.titlesize': 8,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
})

# ─── Colors ──────────────────────────────────────────────────────────────────
C_R = '#2E5AAC'  # Responder blue
C_NR = '#B83A3A'  # Non-responder red
CMAP_R = 'Blues'
CMAP_NR = 'Reds'

# ─── Gene groups for annotation ──────────────────────────────────────────────
GENE_GROUPS = {
    'Treg core': ['CCR8', 'FOXP3', 'IL2RA', 'CTLA4'],
    'Proliferation': ['MKI67'],
    'MHC-I': ['HLA-A', 'B2M'],
    'MHC-II / MIF': ['HLA-DRB5', 'CD74'],
    'Effector / Cytokine': ['GZMB', 'PRF1', 'IL10', 'TGFB1', 'CXCL13'],
    'Cytoskeleton': ['ACTA2', 'MYL9', 'TLN1', 'VCL'],
    'Other': ['TIGIT', 'CCL3', 'CCL4', 'HSF1'],
}
ALL_GENES = []
for group, genes in GENE_GROUPS.items():
    ALL_GENES.extend(genes)

# ════════════════════════════════════════════════════════════════════════════════
# 1. Load data
# ════════════════════════════════════════════════════════════════════════════════
print("Loading data...")
df = pd.read_csv(DATA_FILE)
print(f"  Cells: {df['cell_id'].nunique()}, Genes: {df['gene'].nunique()}")

# Filter to our gene list
df = df[df['gene'].isin(ALL_GENES)].copy()
print(f"  After filtering: {len(df)} rows")

# ════════════════════════════════════════════════════════════════════════════════
# 2. Bin pseudotime and compute mean expression per gene × bin
# ════════════════════════════════════════════════════════════════════════════════
print("Binning pseudotime...")
pt_min, pt_max = df['pseudotime'].min(), df['pseudotime'].max()
print(f"  Pseudotime range: {pt_min:.2f} – {pt_max:.2f}")

N_BINS = 25
bins = np.linspace(pt_min, pt_max, N_BINS + 1)
bin_centers = (bins[:-1] + bins[1:]) / 2

def build_heatmap(sub_df, label):
    """Build a genes × bins matrix of mean expression."""
    heat = np.zeros((len(ALL_GENES), N_BINS))
    for gi, gene in enumerate(ALL_GENES):
        gdf = sub_df[sub_df['gene'] == gene]
        if len(gdf) == 0:
            continue
        for bi in range(N_BINS):
            mask = (gdf['pseudotime'] >= bins[bi]) & (gdf['pseudotime'] < bins[bi + 1])
            if mask.sum() > 0:
                heat[gi, bi] = gdf.loc[mask, 'expression'].mean()
    return heat

# Build R and NR heatmaps
df_r = df[df['response_binary'] == 'Responder']
df_nr = df[df['response_binary'] == 'Non-responder']
print(f"  R cells: {df_r['cell_id'].nunique()}, NR cells: {df_nr['cell_id'].nunique()}")

heat_r = build_heatmap(df_r, 'R')
heat_nr = build_heatmap(df_nr, 'NR')

# ════════════════════════════════════════════════════════════════════════════════
# 3. Determine gene order by peak expression time (average of R and NR)
# ════════════════════════════════════════════════════════════════════════════════
print("Ordering genes by peak expression...")
# Smooth each row with Gaussian to find robust peak
sigma = 1.2
peak_bins = []
for gi in range(len(ALL_GENES)):
    # Use average of R and NR for peak finding
    row_avg = (heat_r[gi] + heat_nr[gi]) / 2
    row_smooth = gaussian_filter1d(row_avg, sigma=sigma)
    peak_bins.append(np.argmax(row_smooth))

# Sort by peak position
gene_order = np.argsort(peak_bins)
sorted_genes = [ALL_GENES[i] for i in gene_order]
print(f"  Gene order (early→late):")
for i, gi in enumerate(gene_order):
    print(f"    {i+1:2d}. {ALL_GENES[gi]:10s}  peak bin={peak_bins[gi]:2d}")

# Reorder heatmaps
heat_r_ordered = heat_r[gene_order]
heat_nr_ordered = heat_nr[gene_order]

# ════════════════════════════════════════════════════════════════════════════════
# 4. Build gene group color annotations
# ════════════════════════════════════════════════════════════════════════════════
group_colors = {
    'Treg core': '#8bc98b',
    'Proliferation': '#e6a8d7',
    'MHC-I': '#2E5AAC',
    'MHC-II / MIF': '#56B4E9',
    'Effector / Cytokine': '#E69F00',
    'Cytoskeleton': '#D55E00',
    'Other': '#999999',
}

# Build annotation bar (gene groups in sorted order)
gene_to_group = {}
for group, genes in GENE_GROUPS.items():
    for g in genes:
        gene_to_group[g] = group

annot_colors = []
for gene in sorted_genes:
    annot_colors.append(group_colors[gene_to_group[gene]])

# ════════════════════════════════════════════════════════════════════════════════
# 5. Plot
# ════════════════════════════════════════════════════════════════════════════════
print("Plotting...")

# Global color scale (shared between R and NR)
vmin = min(heat_r.min(), heat_nr.min())
vmax = max(heat_r.max(), heat_nr.max())

fig = plt.figure(figsize=(8, 6.5))

# Grid layout
# [0] = annotation bar | [1] = heatmap R | [2] = heatmap NR | [3] = colorbar
gs = fig.add_gridspec(1, 5, width_ratios=[0.15, 1.8, 0.15, 1.8, 0.05],
                       wspace=0.02, left=0.08, right=0.95, bottom=0.08, top=0.95)

# ── (a) Gene group annotation bar ──
ax_annot = fig.add_subplot(gs[0])
for gi, color in enumerate(annot_colors):
    ax_annot.barh(gi, 1, height=0.9, color=color, edgecolor='none')
ax_annot.set_xlim(0, 1)
ax_annot.set_ylim(-0.5, len(sorted_genes) - 0.5)
ax_annot.invert_yaxis()
ax_annot.axis('off')

# ── (b) Heatmap R ──
ax_r = fig.add_subplot(gs[1])
im_r = ax_r.imshow(heat_r_ordered, aspect='auto', cmap=CMAP_R,
                   vmin=vmin, vmax=vmax,
                   extent=[bins[0], bins[-1], len(sorted_genes) - 0.5, -0.5])
ax_r.set_xlabel('Pseudotime', fontsize=7.5)
ax_r.set_ylabel('Gene', fontsize=7.5)
ax_r.set_title('Responder', fontsize=8, fontweight='bold', color=C_R, pad=4)
ax_r.set_xticks(np.linspace(bins[0], bins[-1], 5))
ax_r.tick_params(labelsize=6)
ax_r.set_yticks(range(len(sorted_genes)))
ax_r.set_yticklabels(sorted_genes, fontsize=6.5)
ax_r.xaxis.set_ticks_position('bottom')

# ── (c) Separator ──
ax_sep = fig.add_subplot(gs[2])
ax_sep.axis('off')

# ── (d) Heatmap NR ──
ax_nr = fig.add_subplot(gs[3])
im_nr = ax_nr.imshow(heat_nr_ordered, aspect='auto', cmap=CMAP_NR,
                     vmin=vmin, vmax=vmax,
                     extent=[bins[0], bins[-1], len(sorted_genes) - 0.5, -0.5])
ax_nr.set_xlabel('Pseudotime', fontsize=7.5)
ax_nr.set_ylabel('')
ax_nr.set_title('Non-Responder', fontsize=8, fontweight='bold', color=C_NR, pad=4)
ax_nr.set_xticks(np.linspace(bins[0], bins[-1], 5))
ax_nr.tick_params(labelsize=6)
ax_nr.set_yticks([])
ax_nr.xaxis.set_ticks_position('bottom')

# ── (e) Color bar ──
cbar_ax = fig.add_subplot(gs[4])
cbar = fig.colorbar(im_r, cax=cbar_ax, ticks=[vmin, (vmin+vmax)/2, vmax])
cbar.set_label('Mean Expression\n(log₁₀ CPM+1)', fontsize=6.5, labelpad=2)
cbar.ax.tick_params(labelsize=5.5)

# ── Legend for gene groups ──
from matplotlib.patches import Patch
legend_handles = [Patch(facecolor=color, edgecolor='none', label=group)
                  for group, color in group_colors.items()]
# Place legend below the heatmaps
fig.legend(handles=legend_handles, loc='lower center',
           ncol=4, fontsize=6, frameon=True, framealpha=0.9,
           bbox_to_anchor=(0.55, 0.01), handlelength=1)

# Adjust y-limits of heatmaps to leave room for legend at bottom
for ax in [ax_r, ax_nr, ax_annot]:
    ylim = ax.get_ylim()
    ax.set_ylim(ylim[0], ylim[1] + 0.5)

# ── Panel label ──
ax_r.text(-0.18, 1.02, 'E', transform=ax_r.transAxes, fontsize=14,
          fontweight='bold', va='bottom', ha='right')

# ── Save ──
for ext in ['png', 'pdf']:
    path = os.path.join(OUT_DIR, f'fig6e_pseudotime_heatmap.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {path}")

plt.close(fig)
print("\nDone!")
