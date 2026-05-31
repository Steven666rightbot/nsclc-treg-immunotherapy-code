#!/usr/bin/env python3
"""
Figure 3B — HSF1 Target Gene Expression Heatmap
================================================

Cells: MKI67+ Treg only (Responder + Non-responder)
Genes: HSF1 target genes from literature (HSF1 not found in pySCENIC regulons)
Layout: Columns = cells (Responder first, then Non-responder), Rows = HSF1 target genes
Color: Blue-white-red diverging (RdBu_r), Z-score per gene
Row clustering: hierarchical (complete linkage, Euclidean distance)
Column order: manual (no clustering), R then NR

Outputs:
  - figures/fig3b_hsf1_heatmap.png
  - figures/fig3b_hsf1_heatmap.pdf
  - results/fig3b_hsf1_targets.csv
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = r'D:\research\tomato'
LOOM_PATH     = os.path.join(BASE_DIR, r'results\scenic_input\treg_no_foxp3_8k.loom')
REGULONS_PATH = os.path.join(BASE_DIR, r'results\scenic_input\scenic_treg_8k_regulons.csv')
META_PATH     = os.path.join(BASE_DIR, r'results\scenic_input\treg_no_foxp3_metadata.csv')
OUT_FIG_PNG   = os.path.join(BASE_DIR, r'figures\fig3b_hsf1_heatmap.png')
OUT_FIG_PDF   = os.path.join(BASE_DIR, r'figures\fig3b_hsf1_heatmap.pdf')
OUT_TARGETS   = os.path.join(BASE_DIR, r'results\fig3b_hsf1_targets.csv')

os.makedirs(os.path.dirname(OUT_FIG_PNG), exist_ok=True)
os.makedirs(os.path.dirname(OUT_TARGETS), exist_ok=True)

# ── HSF1 target genes (literature-based, since not in pySCENIC regulons) ───────
HSF1_TARGETS_LIT = [
    'HSPA1A', 'HSPA1B', 'HSP90AA1', 'HSP90AB1', 'HSPB1', 'HSPE1',
    'DNAJB1', 'BAG3', 'HSPH1', 'HSPA8', 'HSPA4', 'SERPINH1',
    'CLU', 'CRYAB', 'ABCF2', 'HSF1',           # include HSF1 itself
    'HSPA6', 'HSPA7', 'HSPA2', 'HSPA1L',
    'DNAJA1', 'DNAJA4', 'DNAJB4', 'DNAJB6',
    'HSPBP1', 'STIP1', 'PTGES3', 'FKBP4',
    'HSP90B1', 'TRA1', 'SIL1', 'HYOU1',
]

# ── 1. Load metadata ─────────────────────────────────────────────────────────
meta = pd.read_csv(META_PATH)
meta['CellID'] = meta['CellID'].astype(str)

# Identify MKI67+ Treg cells
mkidf = meta[meta['sub_cell_type'].str.contains('MKI67', case=False, na=False)].copy()
print(f"[INFO] MKI67+ Treg cells: {len(mkidf)}")

# ── 2. Load expression matrix from loom ──────────────────────────────────────
import loompy

print("[INFO] Loading loom file ...")
ds = loompy.connect(LOOM_PATH, 'r')

# Gene names (row axis)
gene_names = ds.ra['Gene'][:].astype(str)
gene_idx_map = {g.upper(): i for i, g in enumerate(gene_names)}

# Cell barcodes and metadata from loom
loom_cell_ids = ds.ca['CellID'][:].astype(str)
loom_sub_ct    = ds.ca['sub_cell_type'][:].astype(str)
loom_response  = ds.ca['response'][:].astype(str)

# Map metadata to loom column indices
# Use the metadata file cell order for filtering MKI67+ cells
mkidf_sorted = mkidf.sort_values('response', ascending=False).copy()
# Responder first (ascending=False means 'Responder' < 'Non-responder' string-sorted)
# But we want Responder first then Non-responder
# Let's sort: R first, then NR
order_map = {'Responder': 0, 'Non-responder': 1}
mkidf_sorted['_order'] = mkidf_sorted['response'].map(order_map)
mkidf_sorted = mkidf_sorted.sort_values('_order').reset_index(drop=True)

mki_cell_ids = mkidf_sorted['CellID'].tolist()
mki_responses = mkidf_sorted['response'].tolist()
mki_subtypes   = mkidf_sorted['sub_cell_type'].tolist()

# Find loom column indices for MKI67+ cells
loom_cell_id_list = list(loom_cell_ids)
col_indices = []
valid_cell_ids = []
valid_responses = []
valid_subtypes = []
for cid, resp, sub in zip(mki_cell_ids, mki_responses, mki_subtypes):
    try:
        idx = loom_cell_id_list.index(cid)
        col_indices.append(idx)
        valid_cell_ids.append(cid)
        valid_responses.append(resp)
        valid_subtypes.append(sub)
    except ValueError:
        print(f"[WARN] CellID {cid} not found in loom; skipping.")

print(f"[INFO] MKI67+ cells found in loom: {len(col_indices)}")

# ── 3. Select HSF1 target genes available in dataset ─────────────────────────
# First, check if HSF1 regulon exists in the pySCENIC output
hsf1_regulon_genes = None
if os.path.exists(REGULONS_PATH):
    regulons_df = pd.read_csv(REGULONS_PATH, skiprows=2)
    hsf1_rows = regulons_df[regulons_df['TF'].str.strip().str.upper() == 'HSF1']
    if len(hsf1_rows) > 0:
        print("[INFO] HSF1 regulon found in pySCENIC output!")
        # Parse TargetGenes column
        import ast
        target_str = hsf1_rows.iloc[0]['TargetGenes']
        # TargetGenes is string representation of list of tuples
        target_list = ast.literal_eval(target_str)
        hsf1_regulon_genes = [t[0] for t in target_list]
        print(f"[INFO] HSF1 regulon target genes ({len(hsf1_regulon_genes)}): {hsf1_regulon_genes[:10]}...")
    else:
        print("[INFO] HSF1 not found in pySCENIC regulons. Using literature-based target genes.")

if hsf1_regulon_genes is not None:
    target_gene_symbols = hsf1_regulon_genes
else:
    target_gene_symbols = HSF1_TARGETS_LIT

# Map target genes to loom indices
available_genes = []
missing_genes = []
for g in target_gene_symbols:
    if g.upper() in gene_idx_map:
        available_genes.append((g, gene_idx_map[g.upper()]))
    else:
        missing_genes.append(g)

print(f"[INFO] Available HSF1 target genes: {len(available_genes)}")
if missing_genes:
    print(f"[INFO] Missing from dataset: {missing_genes}")

# If more than 30 target genes, keep top 30 by mean expression across MKI67+ cells
if len(available_genes) > 30:
    print(f"[INFO] >30 target genes available; selecting top 30 by mean expression ...")
    # Compute mean expression for each gene across MKI67+ cells
    means = []
    for gname, gidx in available_genes:
        expr = ds[[gidx], :][0, :]  # 1D array over all cells
        mean_expr = np.mean(expr[col_indices])
        means.append((gname, mean_expr, gidx))
    means.sort(key=lambda x: x[1], reverse=True)
    available_genes = [(g, idx) for g, m, idx in means[:30]]
    print(f"[INFO] Top 30: {[g for g, idx in available_genes]}")

# Extract expression matrix for selected genes and cells
n_genes = len(available_genes)
n_cells = len(col_indices)
expr_matrix = np.zeros((n_genes, n_cells), dtype=np.float64)

gene_labels = []
for i, (gname, gidx) in enumerate(available_genes):
    expr_matrix[i, :] = ds[[gidx], :][0, col_indices]
    gene_labels.append(gname)

ds.close()

print(f"[INFO] Expression matrix shape: {expr_matrix.shape}")

# ── 4. Z-score normalize per gene (StandardScaler) ───────────────────────────
scaler = StandardScaler()
expr_z = scaler.fit_transform(expr_matrix)  # shape (n_genes, n_cells)

# Clip extreme Z-scores for better visualization
expr_z = np.clip(expr_z, -3, 3)

# ── 5. Row clustering (hierarchical, complete linkage, Euclidean) ─────────────
row_dist = pdist(expr_z, metric='euclidean')
row_link = linkage(row_dist, method='complete')
row_order = leaves_list(row_link)

expr_z_clustered = expr_z[row_order, :]
gene_labels_clustered = [gene_labels[i] for i in row_order]

print(f"[INFO] Row clustering complete. Order: {gene_labels_clustered}")

# ── 6. Build annotation colors ───────────────────────────────────────────────
# Responder = green, Non-responder = purple
resp_colors = {'Responder': '#4CAF50', 'Non-responder': '#9C27B0'}
annotation_colors = [resp_colors[r] for r in valid_responses]

# ── 7. Create figure ─────────────────────────────────────────────────────────
# Dimensions: publication quality
# Height depends on number of genes
row_height = 0.45  # inches per gene row
heatmap_width = 8.0
label_width = 2.0        # space for gene labels on the right
annotation_height = 0.35  # top annotation bar
colorbar_width = 0.3

fig_height = max(3.0, n_genes * row_height + annotation_height + 1.0)
fig_width = heatmap_width + label_width + colorbar_width + 0.5

fig = plt.figure(figsize=(fig_width, fig_height))

# Layout: [heatmap | gene_labels | colorbar]
# Top bar for annotations
ax_annot = fig.add_axes([0.05, 1 - annotation_height/fig_height,
                         heatmap_width/fig_width, annotation_height/fig_height])
ax_heat = fig.add_axes([0.05, 0.12,
                        heatmap_width/fig_width, 1 - 0.12 - annotation_height/fig_height - 0.02])
ax_cbar = fig.add_axes([0.05 + heatmap_width/fig_width + 0.01, 0.12,
                        colorbar_width/fig_width, 1 - 0.12 - annotation_height/fig_height - 0.02])

# ── 7a. Top annotation bar ──────────────────────────────────────────────────
# Response bar (top half)
bar_height = annotation_height / fig_height * 0.4
bar_y = 1 - bar_height - 0.01

for ci in range(n_cells):
    color = annotation_colors[ci]
    ax_annot.add_patch(Rectangle(
        (ci / n_cells, 0.55), 1 / n_cells, 0.45,
        facecolor=color, edgecolor='none', linewidth=0
    ))
    # Subtype marker (bottom half) - use a lighter shade of same color
    sub_color = '#E8F5E9' if valid_responses[ci] == 'Responder' else '#F3E5F5'
    ax_annot.add_patch(Rectangle(
        (ci / n_cells, 0.05), 1 / n_cells, 0.45,
        facecolor=sub_color, edgecolor='none', linewidth=0
    ))

ax_annot.set_xlim(0, 1)
ax_annot.set_ylim(0, 1)
ax_annot.axis('off')

# Add annotation labels on the left
ax_annot.text(-0.02, 0.75, 'Response', transform=ax_annot.transAxes,
              fontsize=7, fontname='Arial', va='center', ha='right')
ax_annot.text(-0.02, 0.25, 'Subtype', transform=ax_annot.transAxes,
              fontsize=7, fontname='Arial', va='center', ha='right')

# Legend for response groups
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#4CAF50', label='Responder'),
    Patch(facecolor='#9C27B0', label='Non-responder'),
]
ax_annot.legend(handles=legend_elements, loc='upper right', fontsize=6,
                frameon=True, framealpha=0.9, edgecolor='#cccccc',
                title='Response', title_fontsize=7)

# ── 7b. Heatmap ─────────────────────────────────────────────────────────────
im = ax_heat.imshow(expr_z_clustered, aspect='auto', cmap='RdBu_r',
                    interpolation='nearest', vmin=-3, vmax=3)

ax_heat.set_xticks([])  # no column labels
# Gene labels on right side
ax_heat.set_yticks(range(n_genes))
ax_heat.set_yticklabels(gene_labels_clustered, fontsize=8, fontname='Arial')
ax_heat.tick_params(axis='y', left=False)  # remove tick marks

# Horizontal lines between genes
for i in range(n_genes - 1):
    ax_heat.axhline(i + 0.5, color='#eeeeee', linewidth=0.3)

# Subtle vertical line between Responder and Non-responder groups
n_resp = sum(1 for r in valid_responses if r == 'Responder')
if 0 < n_resp < n_cells:
    ax_heat.axvline(n_resp - 0.5, color='#333333', linewidth=1.2, linestyle='-', alpha=0.6)

# Add group labels on top of heatmap
ax_heat.text(n_resp / 2 - 0.5, -1.5, 'Responder',
             fontsize=9, fontname='Arial', ha='center', va='bottom',
             color='#4CAF50', fontweight='bold')
ax_heat.text(n_resp + (n_cells - n_resp) / 2 - 0.5, -1.5, 'Non-responder',
             fontsize=9, fontname='Arial', ha='center', va='bottom',
             color='#9C27B0', fontweight='bold')

# ── 7c. Colorbar ─────────────────────────────────────────────────────────────
cbar = fig.colorbar(im, cax=ax_cbar)
cbar.set_label('Z-score (normalized expression)', fontsize=8, fontname='Arial')
cbar.ax.tick_params(labelsize=7)

# ── 7d. Panel label ──────────────────────────────────────────────────────────
fig.text(0.01, 0.98, 'B', fontsize=14, fontname='Arial', fontweight='bold',
         transform=fig.transFigure)

# ── 8. Save ──────────────────────────────────────────────────────────────────
fig.savefig(OUT_FIG_PNG, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_FIG_PDF, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)

print(f"[DONE] Heatmap saved to:")
print(f"       {OUT_FIG_PNG}")
print(f"       {OUT_FIG_PDF}")

# ── 9. Save target gene list ─────────────────────────────────────────────────
pd.DataFrame({
    'gene': gene_labels,
    'symbol': gene_labels,
    'source': 'literature' if hsf1_regulon_genes is None else 'pySCENIC_regulon'
}).to_csv(OUT_TARGETS, index=False)

print(f"[DONE] Target gene list saved to: {OUT_TARGETS}")
print(f"[INFO] Used {n_genes} HSF1 target genes (HSF1 "
      + ("found in pySCENIC regulons)" if hsf1_regulon_genes else "not in pySCENIC; literature-based)"))
print("[INFO] All done.")
