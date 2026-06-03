#!/usr/bin/env python3
"""
Figure 3B — HSF1 Target Gene Expression Heatmap (Pseudobulk)
=============================================================
Cells: MKI67+ Treg only (Responder + Non-responder)
Genes: HSF1 target genes from literature (not in pySCENIC regulons)
Layout: Pseudobulk heatmap — columns = R/NR groups
         Shows mean expression per group, Z-score normalized per gene

Fixed from single-cell heatmap: now aggregates to pseudobulk for
clean, publication-ready visualization.
"""

import os, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler
import loompy

warnings.filterwarnings('ignore')

# ── Paths ──
DATA_DIR = r'D:\\research\\tomato'
OUT_DIR  = r'D:\\research\\cucumber\\fig3'
LOOM_PATH     = os.path.join(DATA_DIR, 'results/scenic_input/treg_no_foxp3_8k.loom')
META_PATH     = os.path.join(DATA_DIR, 'results/scenic_input/treg_no_foxp3_metadata.csv')

os.makedirs(OUT_DIR, exist_ok=True)

# ── HSF1 target genes (curated, 29 targets) ──
HSF1_TARGETS = [
    'HSPA1A', 'HSPA1B', 'HSP90AA1', 'HSP90AB1', 'HSPB1', 'HSPE1',
    'DNAJB1', 'BAG3', 'HSPH1', 'HSPA8', 'HSPA4', 'SERPINH1',
    'CRYAB', 'ABCF2',
    'HSPA6', 'HSPA1L', 'DNAJA1', 'DNAJA4', 'DNAJB4', 'DNAJB6',
    'HSPBP1', 'STIP1', 'PTGES3', 'FKBP4',
    'HSP90B1', 'SIL1', 'HYOU1',
    'HSPA2', 'CLU',
]

# ── 1. Load metadata, filter MKI67+ ──
meta = pd.read_csv(META_PATH)
meta['CellID'] = meta['CellID'].astype(str)
mkidf = meta[meta['sub_cell_type'].str.contains('MKI67', case=False, na=False)].copy()
print(f"MKI67+ Treg cells in metadata: {len(mkidf)}")

# ── 2. Load expression from loom ──
ds = loompy.connect(LOOM_PATH, 'r')
gene_names = ds.ra['Gene'][:].astype(str)
gene_idx_map = {g.upper(): i for i, g in enumerate(gene_names)}
loom_cell_ids = ds.ca['CellID'][:].astype(str)
loom_response = ds.ca['response'][:].astype(str)

# Map MKI67+ cells to loom columns
loom_cell_id_list = list(loom_cell_ids)
col_indices = []
cell_responses = []
for _, row in mkidf.iterrows():
    try:
        idx = loom_cell_id_list.index(row['CellID'])
        col_indices.append(idx)
        cell_responses.append(row['response'])
    except ValueError:
        pass

print(f"MKI67+ cells found in loom: {len(col_indices)}")
cell_responses = np.array(cell_responses)

# ── 3. Get HSF1 target expression ──
available_genes = []
for g in HSF1_TARGETS:
    if g.upper() in gene_idx_map:
        available_genes.append((g, gene_idx_map[g.upper()]))
    else:
        print(f"  Missing: {g}")

print(f"Available HSF1 targets: {len(available_genes)}")

# Extract expression
n_genes = len(available_genes)
n_cells = len(col_indices)
expr = np.zeros((n_genes, n_cells))
gene_labels = []
for i, (gname, gidx) in enumerate(available_genes):
    expr[i, :] = ds[[gidx], :][0, col_indices]
    gene_labels.append(gname)
ds.close()

# ── 4. Pseudobulk: average per group ──
r_mask = cell_responses == 'Responder'
nr_mask = cell_responses == 'Non-responder'

r_mean = expr[:, r_mask].mean(axis=1)
nr_mean = expr[:, nr_mask].mean(axis=1)

# Also compute % expressing and FC
r_pct = (expr[:, r_mask] > 0).mean(axis=1) * 100
nr_pct = (expr[:, nr_mask] > 0).mean(axis=1) * 100

pb_df = pd.DataFrame({
    'gene': gene_labels,
    'R_mean': r_mean,
    'NR_mean': nr_mean,
    'R_pct': r_pct,
    'NR_pct': nr_pct,
    'log2FC_R_vs_NR': np.log2((r_mean + 0.01) / (nr_mean + 0.01))
})
pb_df.to_csv(os.path.join(OUT_DIR, 'fig3b_hsf1_pseudobulk.csv'), index=False)
print(f"\nPseudobulk saved. Top genes by log2FC:")
print(pb_df.sort_values('log2FC_R_vs_NR').head(10)[['gene', 'R_mean', 'NR_mean', 'log2FC_R_vs_NR']].to_string(index=False))

# ── 5. Create two-panel figure ──
# Panel B1: Pseudobulk heatmap (2 columns x N genes)
# Panel B2: Dot plot (% expressing vs mean expression)

# Z-score normalize the pseudobulk matrix (rows)
pb_matrix = np.column_stack([r_mean, nr_mean])
scaler = StandardScaler()
pb_z = scaler.fit_transform(pb_matrix)
pb_z = np.clip(pb_z, -2, 2)

# Row clustering
row_dist = pdist(pb_matrix, metric='euclidean')
row_link = linkage(row_dist, method='complete')
row_order = leaves_list(row_link)

pb_z_clustered = pb_z[row_order, :]
gene_labels_clustered = [gene_labels[i] for i in row_order]

# ── 6. Figure layout: side-by-side heatmap + dotplot ──
fig, axes = plt.subplots(1, 2, figsize=(12, max(5, n_genes * 0.32)),
                         gridspec_kw={'width_ratios': [1.8, 2.5], 'wspace': 0.3})

# ── Panel B1: Pseudobulk heatmap ──
ax = axes[0]
im = ax.imshow(pb_z_clustered, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)

ax.set_yticks(range(n_genes))
ax.set_yticklabels(gene_labels_clustered, fontsize=6.5)
ax.set_xticks([0, 1])
ax.set_xticklabels(['Responder', 'Non-responder'], fontsize=9, rotation=0)

# Annotate values in cells
for i in range(n_genes):
    for j in range(2):
        val = pb_z_clustered[i, j]
        color = 'white' if abs(val) > 1.2 else 'black'
        ax.text(j, i, f'{pb_matrix[row_order[i], j]:.2f}',
                ha='center', va='center', fontsize=7, color=color, fontweight='bold')

ax.set_title('Mean Expression', fontsize=10, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Colorbar
cbar = plt.colorbar(im, ax=ax, fraction=0.15, pad=0.04)
cbar.set_label('Z-score', fontsize=8)

# ── Panel B2: Dot plot (% expressing vs expression level) ──
ax = axes[1]
# Separate R and NR
for i, gene in enumerate(gene_labels_clustered):
    # Find original index
    orig_idx = gene_labels.index(gene)
    
    # Responder dot
    ax.scatter(r_mean[orig_idx], i, s=r_pct[orig_idx] * 3, 
               c='#2E5AAC', edgecolors='#555555', linewidth=0.3, alpha=0.8, zorder=3)
    # Non-responder dot
    ax.scatter(nr_mean[orig_idx], i, s=nr_pct[orig_idx] * 3,
               c='#C0392B', edgecolors='#555555', linewidth=0.3, alpha=0.8, zorder=3)

# Connect paired dots with lines
for i, gene in enumerate(gene_labels_clustered):
    orig_idx = gene_labels.index(gene)
    ax.plot([r_mean[orig_idx], nr_mean[orig_idx]], [i, i],
            color='#cccccc', linewidth=0.4, zorder=1)

ax.set_yticks(range(n_genes))
ax.set_yticklabels(gene_labels_clustered, fontsize=6.5)
ax.set_xlabel('Mean Expression', fontsize=9)
ax.set_title('Expression × % Expressing', fontsize=10, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.invert_yaxis()

# Legend
legend_elements = [
    plt.scatter([], [], s=30, c='#2E5AAC', edgecolors='#555555', linewidth=0.3, label='Responder'),
    plt.scatter([], [], s=30, c='#C0392B', edgecolors='#555555', linewidth=0.3, label='Non-responder'),
]
# Size legend
size_handles = []
for pct in [20, 50, 80]:
    size_handles.append(plt.scatter([], [], s=pct*3, c='gray', edgecolors='#555555', linewidth=0.3, label=f'{pct}%'))
legend2 = ax.legend(handles=size_handles, title='% Expressing', loc='lower right',
                     frameon=False, fontsize=7, title_fontsize=7.5)
ax.add_artist(ax.legend(handles=legend_elements, loc='upper right', frameon=False, fontsize=7))

# Panel label
fig.text(0.01, 0.98, 'B', fontsize=14, fontweight='bold', transform=fig.transFigure)

fig.savefig(os.path.join(OUT_DIR, 'fig3b_hsf1_heatmap.png'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(OUT_DIR, 'fig3b_hsf1_heatmap.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)

print(f"\n[DONE] Figure 3B saved to {OUT_DIR}/")
print(f"  - fig3b_hsf1_heatmap.png")
print(f"  - fig3b_hsf1_heatmap.pdf")
