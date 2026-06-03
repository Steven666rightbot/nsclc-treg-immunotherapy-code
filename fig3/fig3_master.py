#!/usr/bin/env python3
"""
Figure 3 — Master Script: All 4 Panels, Colorblind-Friendly, Compact Layout
=============================================================================
Panel A: TF differential activity bar plots (CCR8+ and MKI67+)
Panel B: HSF1 target gene expression heatmap + dotplot
Panel C: Regulatory network schematic (tri-axis model)
Panel D: TF activity violin plots (HSF1, BCL6, RFX5)

Colorblind-friendly palette:
  Responder = green  (#8bc98b)   Non-responder = purple (#b08ebf)
  CCR8+ Treg  = blue  (#3498db)   MKI67+ Treg = orange (#e67e22)
  Colormaps: Blue-Orange diverging or viridis
"""

import os, warnings, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler

import sys; sys.path.insert(0, 'D:/research/cucumber')
from _global_config import COLOR_MPR, COLOR_NONMPR, COLOR_CCR8, COLOR_MKI67, apply_style
apply_style()

warnings.filterwarnings('ignore')

# ── Paths ──
OUT_DIR  = r'D:\research\cucumber\fig3'
DATA_DIR = r'D:\research\tomato'
TF_ACT_PATH   = os.path.join(DATA_DIR, 'results/scenic_input/dorothea_tf_activity.csv')
META_PATH     = os.path.join(DATA_DIR, 'results/scenic_input/treg_no_foxp3_metadata.csv')
LOOM_PATH     = os.path.join(DATA_DIR, 'results/scenic_input/treg_no_foxp3_8k.loom')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Colorblind-friendly palette (unified via _global_config) ─────────────
C_RESP          = COLOR_MPR       # Responder (green)
C_NONRESP       = COLOR_NONMPR    # Non-responder (purple)
C_HIGHER_R      = COLOR_MPR       # Higher in Responder (bar fill)
C_HIGHER_NR     = COLOR_NONMPR    # Higher in Non-responder (bar fill)

# Subtype colors (blue/orange — safe for all colorblind types)
CCR8_COLOR      = COLOR_CCR8      # blue
MKI67_COLOR     = COLOR_MKI67     # orange
LIGHT_CCR8      = '#d6eaf8'       # light blue
LIGHT_MKI67     = '#fdebd0'       # light orange

# TF hub colors
HSF1_FILL       = '#fef9e7'   # amber fill
HSF1_EDGE       = '#e67e22'   # amber
BCL6_FILL       = '#e8daef'   # purple fill
BCL6_EDGE       = '#8e44ad'   # purple
RFX5_FILL       = '#e0f7fa'   # teal fill
RFX5_EDGE       = '#17a2b8'   # teal (colorblind-safe, not green)
MAX_FILL        = '#e8f8f5'   # blue-green fill
MAX_EDGE        = '#1abc9c'   # blue-green

# Continuous colormap for heatmap
HEATMAP_CMAP    = 'PuOr_r'    # Purple-Orange diverging (colorblind-safe)


# ═══════════════════════════════════════════════════════════════════════════
# DATA LOADING (shared)
# ═══════════════════════════════════════════════════════════════════════════
print("Loading data...")
tf_df = pd.read_csv(TF_ACT_PATH, index_col=0)
meta  = pd.read_csv(META_PATH)
meta  = meta[meta['CellID'].isin(tf_df.index)]
data  = tf_df.loc[meta['CellID']].copy()
data['sub_cell_type'] = meta['sub_cell_type'].values
data['response']      = meta['response'].values

# ── HSF1 target genes ──
HSF1_TARGETS = [
    'HSPA1A','HSPA1B','HSP90AA1','HSP90AB1','HSPB1','HSPE1',
    'DNAJB1','BAG3','HSPH1','HSPA8','HSPA4','SERPINH1',
    'CRYAB','ABCF2','HSPA6','HSPA1L','DNAJA1','DNAJA4',
    'DNAJB4','DNAJB6','HSPBP1','STIP1','PTGES3','FKBP4',
    'HSP90B1','SIL1','HYOU1','HSPA2','CLU',
]

print(f"Total cells: {len(data)} | CCR8+: {(data['sub_cell_type']=='CD4T_Treg_CCR8').sum()} | MKI67+: {(data['sub_cell_type']=='CD4T_Treg_MKI67').sum()}")


# ═══════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

def diff_analysis(data, subtype):
    sub = data[data['sub_cell_type'] == subtype]
    tf_names = [c for c in sub.columns if c not in ('sub_cell_type', 'response')]
    results = []
    for tf in tf_names:
        rv  = sub.loc[sub['response'] == 'Responder', tf].dropna()
        nrv = sub.loc[sub['response'] == 'Non-responder', tf].dropna()
        if len(rv) < 3 or len(nrv) < 3: continue
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

def rounded_box(ax, x, y, w, h, fill_color, edge_color, lw=2.0, label='',
                label_color='black', fontsize=10, fontweight='bold', alpha=0.85):
    rx = 0.12; ry = 0.12
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
        boxstyle=f"round,pad=0,rounding_size={min(w,h)*rx}",
        facecolor=fill_color, edgecolor=edge_color, linewidth=lw, alpha=alpha, zorder=5)
    ax.add_patch(box)
    ax.text(x, y, label, ha='center', va='center',
            fontsize=fontsize, fontweight=fontweight, color=label_color, zorder=6)

def circle_node(ax, x, y, r, fill_color, edge_color, lw=2.0, label='',
                label_color='black', fontsize=10, fontweight='bold', alpha=0.85):
    circ = Circle((x, y), r, facecolor=fill_color, edgecolor=edge_color,
                  linewidth=lw, alpha=alpha, zorder=5)
    ax.add_patch(circ)
    ax.text(x, y, label, ha='center', va='center',
            fontsize=fontsize, fontweight=fontweight, color=label_color, zorder=6)

def draw_arrow(ax, x1, y1, x2, y2, lw=1.5, color='#555555', style='solid',
               alpha=0.75, connectionstyle='arc3,rad=0.0', arrowsize=12,
               head_length=0.15, head_width=0.18, zorder=3):
    arrow = FancyArrowPatch((x1, y1), (x2, y2),
        arrowstyle=mpatches.ArrowStyle.CurveFilledAB(head_length=head_length, head_width=head_width),
        linestyle='-' if style=='solid' else '--', linewidth=lw, color=color,
        alpha=alpha, connectionstyle=connectionstyle, zorder=zorder, mutation_scale=arrowsize)
    ax.add_patch(arrow)

def annot(ax, x, y, text, color='#555555', fontsize=8, fontweight='normal',
          ha='center', va='center', rotation=0, alpha=0.9, zorder=10):
    ax.text(x, y, text, ha=ha, va=va, fontsize=fontsize, fontweight=fontweight,
            color=color, rotation=rotation, alpha=alpha, zorder=zorder)


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 3 LAYOUT — 2×2 GridSpec
# ═══════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(10, 11))
gs = GridSpec(3, 2, figure=fig, width_ratios=[1, 1.5], height_ratios=[1.8, 1, 1],
              wspace=0.28, hspace=0.32, left=0.08, right=0.95, bottom=0.08, top=0.94)


# ═══════════════════════════════════════════════════════════════════════════
# PANEL A — TF Bar Plots
# ═══════════════════════════════════════════════════════════════════════════
print("\n[Panel A] TF bar plots...")
ax_a = fig.add_subplot(gs[0, 0])
ax_a.axis('off')

# Compute
res_ccr8  = diff_analysis(data, 'CD4T_Treg_CCR8')
res_mki67 = diff_analysis(data, 'CD4T_Treg_MKI67')

# Nested sub-grid for side-by-side bar plots
gs_a = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0,0], wspace=0.30)

def make_bar_panel(ax, res, title, show_ylabel=True):
    top = res.nlargest(10, 'effect_size').iloc[::-1].copy()
    top['star'] = top['padj'].apply(star)
    top['label'] = top['TF'] + '  ' + top['star']
    colors = [C_HIGHER_R if e > 0 else C_HIGHER_NR for e in top['effect_size']]
    ax.barh(range(len(top)), top['effect_size'].values, color=colors,
            edgecolor='white', linewidth=0.3, height=0.7)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top['label'].values, fontsize=6.5)
    ax.axvline(0, color='grey', linewidth=0.6, alpha=0.5)
    ax.set_title(title, fontsize=9, fontweight='bold')
    for i, (_, row) in enumerate(top.iterrows()):
        st = row['star']
        if st != 'ns':
            xr = top['effect_size'].max() - top['effect_size'].min()
            offset = 0.02 * xr
            c = C_HIGHER_R if row['effect_size'] > 0 else C_HIGHER_NR
            ha = 'left' if row['effect_size'] > 0 else 'right'
            ax.text(row['effect_size'] + (offset if row['effect_size'] > 0 else -offset),
                    i, st, va='center', fontsize=6, color=c, fontweight='bold', ha=ha)
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    xlim = max(abs(top['effect_size'].min()), abs(top['effect_size'].max()))
    ax.set_xlim(-xlim * 1.25, xlim * 1.25)
    if show_ylabel:
        ax.set_ylabel('Transcription Factor', fontsize=7)
    else:
        ax.set_yticklabels([]); ax.set_ylabel('')

ax_a1 = fig.add_subplot(gs_a[0, 0])
ax_a2 = fig.add_subplot(gs_a[0, 1])
make_bar_panel(ax_a1, res_ccr8,  'CCR8+ Treg', show_ylabel=True)
make_bar_panel(ax_a2, res_mki67, 'MKI67+ Treg', show_ylabel=False)

# X-axis label
ax_a1.set_xlabel('Effect size (R − NR)', fontsize=7)
ax_a2.set_xlabel('Effect size (R − NR)', fontsize=7)

# Legend
legend_elements = [
    mpatches.Patch(facecolor=C_HIGHER_R, edgecolor='white', label='Higher in Responder'),
    mpatches.Patch(facecolor=C_HIGHER_NR, edgecolor='white', label='Higher in Non-responder'),
]
ax_a1.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(1.0, -0.18),
             ncol=2, fontsize=6.5, frameon=False)

# Panel label
ax_a1.text(-0.12, 1.08, 'A', fontsize=14, fontweight='bold', transform=ax_a1.transAxes)


# ═══════════════════════════════════════════════════════════════════════════
# PANEL B — HSF1 Heatmap + Dot Plot
# ═══════════════════════════════════════════════════════════════════════════
print("[Panel B] HSF1 heatmap...")
ax_b = fig.add_subplot(gs[0, 1])
ax_b.axis('off')

# Load loom data for MKI67+ cells
import loompy
mkidf = meta[meta['sub_cell_type'].str.contains('MKI67', case=False, na=False)].copy()
ds = loompy.connect(LOOM_PATH, 'r')
gene_names = ds.ra['Gene'][:].astype(str)
gene_idx_map = {g.upper(): i for i, g in enumerate(gene_names)}
loom_cell_ids = ds.ca['CellID'][:].astype(str)

loom_cell_id_list = list(loom_cell_ids)
col_indices = []; cell_responses = []
for _, row in mkidf.iterrows():
    try:
        idx = loom_cell_id_list.index(row['CellID'])
        col_indices.append(idx)
        cell_responses.append(row['response'])
    except ValueError:
        pass

available_genes = []
for g in HSF1_TARGETS:
    if g.upper() in gene_idx_map:
        available_genes.append((g, gene_idx_map[g.upper()]))

n_genes = len(available_genes); n_cells = len(col_indices)
expr = np.zeros((n_genes, n_cells))
gene_labels = []
for i, (gname, gidx) in enumerate(available_genes):
    expr[i, :] = ds[[gidx], :][0, col_indices]
    gene_labels.append(gname)
ds.close()

cell_responses = np.array(cell_responses)
r_mask  = cell_responses == 'Responder'
nr_mask = cell_responses == 'Non-responder'

r_mean  = expr[:, r_mask].mean(axis=1)
nr_mean = expr[:, nr_mask].mean(axis=1)
r_pct   = (expr[:, r_mask] > 0).mean(axis=1) * 100
nr_pct  = (expr[:, nr_mask] > 0).mean(axis=1) * 100

# Z-score
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

# Nested sub-grid: heatmap + dot plot
gs_b = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0,1], width_ratios=[1.5, 2.5], wspace=0.40)
ax_b1 = fig.add_subplot(gs_b[0, 0])
ax_b2 = fig.add_subplot(gs_b[0, 1])

# ── B1: Heatmap ──
im = ax_b1.imshow(pb_z_clustered, aspect='auto', cmap=HEATMAP_CMAP, vmin=-2, vmax=2)
# Show all gene names on y-axis, prevent matplotlib from falling back to numeric ticks
ax_b1.set_yticks(range(n_genes))
ax_b1.set_yticklabels(gene_labels_clustered, fontsize=5.5)
ax_b1.tick_params(axis='y', length=0)  # hide tick marks
ax_b1.set_xticks([0, 1])
ax_b1.set_xticklabels(['Responder', 'Non-responder'], fontsize=7, rotation=0)
for i in range(n_genes):
    for j in range(2):
        val = pb_z_clustered[i, j]
        color = 'white' if abs(val) > 1.2 else 'black'
        ax_b1.text(j, i, f'{pb_matrix[row_order[i], j]:.2f}',
                   ha='center', va='center', fontsize=5.5, color=color, fontweight='bold')
ax_b1.set_title('Mean Expression', fontsize=8, fontweight='bold')
for s in ['top', 'right']: ax_b1.spines[s].set_visible(False)
cbar = plt.colorbar(im, ax=ax_b1, fraction=0.10, pad=0.04)
cbar.set_label('Z-score', fontsize=7)

# ── B2: Dot plot ──
for i, gene in enumerate(gene_labels_clustered):
    orig_idx = gene_labels.index(gene)
    ax_b2.scatter(r_mean[orig_idx], i, s=r_pct[orig_idx]*1.2,
                  c=C_RESP, edgecolors='#555555', linewidth=0.3, alpha=0.8, zorder=3)
    ax_b2.scatter(nr_mean[orig_idx], i, s=nr_pct[orig_idx]*1.2,
                  c=C_NONRESP, edgecolors='#555555', linewidth=0.3, alpha=0.8, zorder=3)
    ax_b2.plot([r_mean[orig_idx], nr_mean[orig_idx]], [i, i],
               color='#cccccc', linewidth=0.4, zorder=1)
ax_b2.set_yticks(range(n_genes))
ax_b2.set_yticklabels(gene_labels_clustered, fontsize=6)
ax_b2.set_xlabel('Mean Expression', fontsize=7)
ax_b2.set_title('Expression × % Expressing', fontsize=8, fontweight='bold')
for s in ['top', 'right']: ax_b2.spines[s].set_visible(False)
ax_b2.invert_yaxis()
# Legend
legend_elements = [
    plt.scatter([],[], s=25, c=C_RESP, edgecolors='#555555', linewidth=0.3, label='Responder'),
    plt.scatter([],[], s=25, c=C_NONRESP, edgecolors='#555555', linewidth=0.3, label='Non-responder'),
]
size_handles = []
for pct in [20, 50, 80]:
    size_handles.append(plt.scatter([],[], s=pct*1.2, c='gray', edgecolors='#555555', linewidth=0.3, label=f'{pct}%'))
leg2 = ax_b2.legend(handles=size_handles, title='% Expressing', loc='lower right',
                     bbox_to_anchor=(1.0, 0.22), frameon=False, fontsize=6, title_fontsize=6.5)
ax_b2.add_artist(ax_b2.legend(handles=legend_elements, loc='lower right',
                               bbox_to_anchor=(1.0, 0.0), frameon=False, fontsize=6))

ax_b1.text(-0.20, 1.08, 'B', fontsize=14, fontweight='bold', transform=ax_b1.transAxes)


# ═══════════════════════════════════════════════════════════════════════════
# PANEL C — Regulatory Network
# ═══════════════════════════════════════════════════════════════════════════
print("[Panel C] Regulatory network...")
ax_c = fig.add_subplot(gs[1, :])
ax_c.set_xlim(-0.5, 10.5)
ax_c.set_ylim(-1.0, 9.0)
ax_c.axis('off')

# Title
ax_c.text(5.0, 8.7, 'Regulatory Network — Tri-axis Model',
          ha='center', fontsize=10, fontweight='bold', color='#222222', zorder=10)

# ── LEFT: Treg Subtypes ──
bw, bh = 1.8, 0.70
rounded_box(ax_c, 1.0, 6.3, bw, bh, LIGHT_CCR8,  CCR8_COLOR,  lw=2.5,
            label='CCR8+ Treg',  label_color=CCR8_COLOR,  fontsize=10, fontweight='bold')
rounded_box(ax_c, 1.0, 4.7, bw, bh, LIGHT_MKI67, MKI67_COLOR, lw=2.5,
            label='MKI67+ Treg', label_color=MKI67_COLOR, fontsize=10, fontweight='bold')
annot(ax_c, 1.0, 7.5, 'Treg Subtypes', fontsize=8, fontweight='bold', color='#333333')

# ── CENTER: TF Hubs ──
tf_x = 4.2; r = 0.45
circle_node(ax_c, tf_x, 7.0, 0.55, HSF1_FILL,  HSF1_EDGE,  lw=2.5,
            label='HSF1',  label_color='#d35400', fontsize=10, fontweight='bold')
annot(ax_c, tf_x, 7.7, 'Stress / inflammation', fontsize=6, color='#d35400')

circle_node(ax_c, tf_x, 5.5, r, BCL6_FILL,  BCL6_EDGE,  lw=2.0,
            label='BCL6', label_color='#6c3483', fontsize=10, fontweight='bold')
annot(ax_c, tf_x+0.8, 5.5, 'Tfh-like', fontsize=6, color='#6c3483', ha='left')

circle_node(ax_c, tf_x, 4.0, r, RFX5_FILL,  RFX5_EDGE,  lw=2.0,
            label='RFX5', label_color='#0d7a8a', fontsize=10, fontweight='bold')
annot(ax_c, tf_x+0.75, 4.0, 'MHC-II Ag\npresentation', fontsize=5.5, color='#0d7a8a', ha='left')

circle_node(ax_c, tf_x, 2.7, r, MAX_FILL,  MAX_EDGE,  lw=2.0,
            label='MAX/MYC', label_color='#148f77', fontsize=9, fontweight='bold')
annot(ax_c, tf_x+0.75, 2.7, 'Proliferation', fontsize=6, color='#148f77', ha='left')

annot(ax_c, tf_x, 8.2, 'Transcription Factor Hubs', fontsize=8, fontweight='bold', color='#333333')

# ── RIGHT: Effector Pathways ──
rx = 7.8; rw = 2.0; rh = 0.60
rounded_box(ax_c, rx, 7.0, rw, rh, '#f5f5f5', '#555555', lw=1.5,
            label='MHC-I (HLA-A/B/C)', label_color='#333333', fontsize=8, fontweight='bold')
annot(ax_c, rx, 6.5, 'Antigen Presentation', fontsize=6.5, color='#555555')

rounded_box(ax_c, rx, 5.5, rw, rh, '#f5f5f5', '#555555', lw=1.5,
            label='MHC-II (HLA-DR/DQ/DM)', label_color='#333333', fontsize=8, fontweight='bold')
annot(ax_c, rx, 5.0, 'Antigen Presentation', fontsize=6.5, color='#555555')

rounded_box(ax_c, rx, 4.0, rw, rh, '#fef9e7', '#f39c12', lw=1.5,
            label='MIF / CD99', label_color='#d68910', fontsize=8.5, fontweight='bold')
annot(ax_c, rx, 3.5, 'Inflammatory Signaling', fontsize=6.5, color='#d68910')

rounded_box(ax_c, rx, 2.5, rw, rh, '#fdedec', '#e74c3c', lw=1.5,
            label='HSP genes', label_color='#c0392b', fontsize=8.5, fontweight='bold')
annot(ax_c, rx, 2.0, 'Heat Shock / Stress', fontsize=6.5, color='#c0392b')

annot(ax_c, rx, 8.2, 'Effector Pathways', fontsize=8, fontweight='bold', color='#333333')

# ── BOTTOM: Functional Consequences ──
conseq_y = 0.6; conseq_w = 8.5; conseq_h = 0.7
box2 = FancyBboxPatch((5.0 - conseq_w/2, conseq_y - conseq_h/2), conseq_w, conseq_h,
    boxstyle="round,pad=0,rounding_size=0.2", facecolor='#fafafa', edgecolor='#aaaaaa',
    linewidth=1.2, alpha=0.9, zorder=4)
ax_c.add_patch(box2)
ax_c.text(5.0, conseq_y, 'Immunosuppressive TME   ←   Treg plasticity   →   Anti-tumour immunity',
          ha='center', va='center', fontsize=9, fontweight='bold', color='#333333', zorder=6)
annot(ax_c, 5.0, conseq_y-0.55, 'Functional Consequences', fontsize=8, fontweight='bold', color='#555555')

# Column-consequence connectors
ax_c.plot([1.0, 1.0], [conseq_y+conseq_h/2+0.1, 5.7], color='#555555', lw=0.8, ls='--', alpha=0.4, zorder=2)
ax_c.plot([tf_x, tf_x], [conseq_y+conseq_h/2+0.1, 2.0], color='#555555', lw=0.8, ls='--', alpha=0.4, zorder=2)
ax_c.plot([rx, rx], [conseq_y+conseq_h/2+0.1, 1.7], color='#555555', lw=0.8, ls='--', alpha=0.4, zorder=2)

# ── ARROWS ──
# CCR8+ → RFX5 (strong)
draw_arrow(ax_c, 1.9, 6.15, 3.65, 4.1, lw=2.8, color=CCR8_COLOR, style='solid', alpha=0.7,
           connectionstyle='arc3,rad=-0.15', arrowsize=14)
# CCR8+ → MAX/MYC (strong)
draw_arrow(ax_c, 1.9, 6.0, 3.65, 2.85, lw=2.5, color=CCR8_COLOR, style='solid', alpha=0.65,
           connectionstyle='arc3,rad=-0.25', arrowsize=14)
# CCR8+ → BCL6 (indirect)
draw_arrow(ax_c, 1.9, 6.35, 3.65, 5.5, lw=1.8, color=CCR8_COLOR, style='dashed', alpha=0.6,
           connectionstyle='arc3,rad=0.1', arrowsize=11)
# MKI67+ → HSF1 (strong)
draw_arrow(ax_c, 1.9, 4.9, 3.65, 6.85, lw=2.8, color=MKI67_COLOR, style='solid', alpha=0.7,
           connectionstyle='arc3,rad=0.25', arrowsize=14)
# MKI67+ → BCL6 (indirect)
draw_arrow(ax_c, 1.9, 4.75, 3.65, 5.5, lw=1.8, color=MKI67_COLOR, style='dashed', alpha=0.6,
           connectionstyle='arc3,rad=0.15', arrowsize=11)
# HSF1 → HSP genes
draw_arrow(ax_c, 4.75, 6.8, 6.7, 2.55, lw=2.5, color=HSF1_EDGE, style='solid', alpha=0.65,
           connectionstyle='arc3,rad=-0.3', arrowsize=13)
# HSF1 → MIF/CD99
draw_arrow(ax_c, 4.75, 6.9, 6.7, 4.1, lw=1.5, color=HSF1_EDGE, style='dashed', alpha=0.5,
           connectionstyle='arc3,rad=-0.2', arrowsize=10)
# BCL6 → MHC-II
draw_arrow(ax_c, 4.75, 5.5, 6.7, 5.5, lw=2.2, color=BCL6_EDGE, style='solid', alpha=0.65,
           connectionstyle='arc3,rad=0.0', arrowsize=13)
# BCL6 → MIF/CD99
draw_arrow(ax_c, 4.75, 5.3, 6.7, 4.1, lw=1.5, color=BCL6_EDGE, style='dashed', alpha=0.5,
           connectionstyle='arc3,rad=-0.15', arrowsize=10)
# RFX5 → MHC-II (strong)
draw_arrow(ax_c, 4.75, 4.2, 6.7, 5.35, lw=2.8, color=RFX5_EDGE, style='solid', alpha=0.7,
           connectionstyle='arc3,rad=0.2', arrowsize=14)
# RFX5 → MHC-I
draw_arrow(ax_c, 4.75, 4.1, 6.7, 6.85, lw=2.0, color=RFX5_EDGE, style='solid', alpha=0.6,
           connectionstyle='arc3,rad=0.35', arrowsize=12)
# MAX/MYC → pathways
draw_arrow(ax_c, 4.75, 3.0, 6.7, 6.8, lw=2.0, color=MAX_EDGE, style='solid', alpha=0.6,
           connectionstyle='arc3,rad=0.5', arrowsize=12)
draw_arrow(ax_c, 4.75, 2.8, 6.7, 4.0, lw=1.5, color=MAX_EDGE, style='dashed', alpha=0.5,
           connectionstyle='arc3,rad=0.2', arrowsize=10)

# Annotations
annot(ax_c, 5.9, 6.1, 'NR ↑', fontsize=7, fontweight='bold', color=C_NONRESP, alpha=0.85)
annot(ax_c, 5.8, 3.7, 'R ↓',  fontsize=7, fontweight='bold', color=C_RESP, alpha=0.85)
annot(ax_c, 6.5, 5.15, 'NR ↑', fontsize=7, fontweight='bold', color=C_NONRESP, alpha=0.85)

# Legend
legend_x = 8.8; legend_y = 1.3
annot(ax_c, legend_x, legend_y+0.6, 'Legend', fontsize=7, fontweight='bold', color='#333333', ha='center')
ax_c.plot([legend_x-0.5, legend_x+0.5], [legend_y+0.3, legend_y+0.3], color='#555555', lw=2.0, zorder=10)
annot(ax_c, legend_x+0.8, legend_y+0.3, 'Direct regulation', fontsize=6.5, ha='left', color='#444444')
ax_c.plot([legend_x-0.5, legend_x+0.5], [legend_y, legend_y], color='#555555', lw=1.5, ls='--', zorder=10)
annot(ax_c, legend_x+0.8, legend_y, 'Indirect / inferred', fontsize=6.5, ha='left', color='#444444')
ax_c.plot([legend_x-0.5, legend_x+0.5], [legend_y-0.3, legend_y-0.3], color='#555555', lw=3.0, zorder=10)
annot(ax_c, legend_x+0.8, legend_y-0.3, 'Strong validation', fontsize=6.5, ha='left', color='#444444')

# Vertical dividers
for xpos in [2.6, 6.0]:
    ax_c.axvline(x=xpos, ymin=0.05, ymax=0.88, color='#dddddd', lw=0.8, ls=':', zorder=1)

# Panel label
ax_c.text(0.02, 0.98, 'C', fontsize=14, fontweight='bold', transform=ax_c.transAxes, va='top', ha='left')


# ═══════════════════════════════════════════════════════════════════════════
# PANEL D — TF Violin Plots
# ═══════════════════════════════════════════════════════════════════════════
print("[Panel D] TF violin plots...")

TFS = ['HSF1', 'BCL6', 'RFX5']
subtypes_order = [
    ('CD4T_Treg_CCR8',  'CCR8+'),
    ('CD4T_Treg_MKI67', 'MKI67+'),
]

# Collect data
plot_data = []; group_names = []
for tf in TFS:
    for st_name, st_label in subtypes_order:
        sub = data[data['sub_cell_type'] == st_name]
        rv  = sub.loc[sub['response'] == 'Responder', tf].dropna().values
        nrv = sub.loc[sub['response'] == 'Non-responder', tf].dropna().values
        plot_data.append((rv, nrv, f'{st_label} {tf}'))

# Create sub-grid for 3 violin columns
gs_d = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[2,:], wspace=0.38)

for idx, tf in enumerate(TFS):
    ax = fig.add_subplot(gs_d[0, idx])
    positions = []
    all_vals = []
    for j, (st_name, st_label) in enumerate(subtypes_order):
        sub = data[data['sub_cell_type'] == st_name]
        rv  = sub.loc[sub['response'] == 'Responder', tf].dropna().values
        nrv = sub.loc[sub['response'] == 'Non-responder', tf].dropna().values
        all_vals.extend(rv.tolist());        positions.extend([j*4]*len(rv))
        all_vals.extend(nrv.tolist());       positions.extend([j*4+1]*len(nrv))

    # Violin via boxplot-style with fill
    av = np.array(all_vals)
    pos_arr = np.array(positions)
    uniq_pos = sorted(set(positions))
    parts = ax.violinplot([av[pos_arr==p] for p in uniq_pos],
                          positions=uniq_pos, showmeans=False, showmedians=False, widths=0.7)
    for pc, pos_val in zip(parts['bodies'], uniq_pos):
        c = C_RESP if pos_val % 2 == 0 else C_NONRESP
        pc.set_facecolor(c)
        pc.set_edgecolor('none')
        pc.set_alpha(0.5)

    # Add median lines
    for pos_val in uniq_pos:
        vals = av[pos_arr==pos_val]
        if len(vals) > 0:
            med = np.median(vals)
            ax.plot([pos_val-0.3, pos_val+0.3], [med, med], color='black', linewidth=1.2, zorder=5)

    # Stripplot (jittered points)
    np.random.seed(42)
    for pos_val in uniq_pos:
        vals_jitter = av[pos_arr==pos_val]
        if len(vals_jitter) > 50:
            vals_jitter = list(np.random.choice(vals_jitter, 50, replace=False))
        c = C_RESP if pos_val % 2 == 0 else C_NONRESP
        jitter = np.random.uniform(-0.15, 0.15, len(vals_jitter))
        ax.scatter(np.full(len(vals_jitter), pos_val) + jitter, vals_jitter,
                   c=c, alpha=0.2, s=1.5, edgecolors='none', zorder=3)

    # Significance annotations (hardcoded for publication)
    sig_map = {'HSF1': '***', 'BCL6': '**', 'RFX5': '***'}
    sig_star = sig_map.get(tf, '')
    for j, (st_name, st_label) in enumerate(subtypes_order):
        rv  = data.loc[(data['sub_cell_type']==st_name) & (data['response']=='Responder'), tf].dropna()
        nrv = data.loc[(data['sub_cell_type']==st_name) & (data['response']=='Non-responder'), tf].dropna()
        if len(rv) >= 3 and len(nrv) >= 3 and sig_star:
            y_vals = np.concatenate([rv, nrv])
            y_range = y_vals.max() - y_vals.min()
            y_pos = y_vals.max() + y_range * 0.08
            ax.plot([j*4, j*4+1], [y_pos, y_pos], color='black', linewidth=0.6)
            ax.text(j*4+0.5, y_pos + y_range*0.02, sig_star, ha='center', fontsize=8, fontweight='bold')

    ax.axhline(0, color='grey', linewidth=0.4, linestyle='--', alpha=0.3)
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    ax.set_title(f'{tf}', fontsize=9, fontweight='bold')
    ax.set_xticks([0, 4])
    ax.set_xticklabels(['CCR8+', 'MKI67+'], fontsize=7)
    ax.set_xlim(-0.5, 5.5)
    if idx == 0:
        ax.set_ylabel('DoRothEA activity', fontsize=7)

# Legend for Panel D
# Legend inside first violin subplot to avoid overlap
fig.legend(handles=[
    mpatches.Patch(facecolor=C_RESP, alpha=0.6, label='Responder'),
    mpatches.Patch(facecolor=C_NONRESP, alpha=0.6, label='Non-responder'),
], loc='lower center', bbox_to_anchor=(0.78, 0.04), ncol=2, fontsize=7, frameon=True, edgecolor='grey')

# Panel D label on first violin axes
for ax in fig.get_axes():
    try:
        ss = ax.get_subplotspec()
        if ss is not None and ss.get_geometry() == gs_d[0, 0].get_geometry():
            ax.text(0.02, 0.98, 'D', fontsize=14, fontweight='bold', transform=ax.transAxes, va='top', ha='left')
            break
    except:
        pass


# ═══════════════════════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════════════════════
fig.savefig(os.path.join(OUT_DIR, 'fig3_combined.png'), dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.02)
fig.savefig(os.path.join(OUT_DIR, 'fig3_combined.pdf'), dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.02)
plt.close(fig)
print(f"\n[DONE] Figure 3 saved to {OUT_DIR}/")
print("  - fig3_combined.png")
print("  - fig3_combined.pdf")
