#!/usr/bin/env python3
"""Figure 6 — Slingshot Pseudotime Trajectory Analysis of CCR8+ / MKI67+ Treg.

Generates 4 panels as PNG + PDF, plus a 2×2 combined layout.
Data source: run_slingshot.R → CSV exports.

Panels:
  6A — UMAP colored by subtype with Slingshot curve overlay
  6B — UMAP colored by pseudotime with curve overlay  
  6C — Pseudotime density split by subtype × response
  6D — Key gene LOESS trends along pseudotime (CCR8, MKI67, CTLA4, CD74)
"""

import os, sys, re
sys.path.insert(0, 'D:/research/cucumber')
from _global_config import COLOR_MPR, COLOR_NONMPR, COLOR_CCR8, COLOR_MKI67, apply_style

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from scipy.interpolate import interp1d

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DATA_DIR = "D:/Research/cucumber/fig6/data"
OUT_DIR  = "D:/Research/cucumber/fig6"

os.makedirs(OUT_DIR, exist_ok=True)

DPI = 300
apply_style()

# Colors (from _global_config)
CCR8_COLOR  = COLOR_CCR8       # blue
MKI67_COLOR = COLOR_MKI67      # orange
RESP_COLOR  = COLOR_MPR        # responder green
NONRESP_COLOR = COLOR_NONMPR   # non-responder purple
PT_CMAP     = 'plasma'         # pseudotime heat colormap

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cells = pd.read_csv(os.path.join(DATA_DIR, "slingshot_pseudotime_umap.csv"))
curves = pd.read_csv(os.path.join(DATA_DIR, "slingshot_curves.csv"))
gexpr = pd.read_csv(os.path.join(DATA_DIR, "slingshot_gene_expr.csv"))

# Filter main trajectory (Lineage 1 & 2)
main = cells[cells['primary_lineage'].isin([1, 2])].copy()
main['subtype_short'] = main['sub_cell_type'].map({
    'CD4T_Treg_CCR8': 'CCR8+ Treg',
    'CD4T_Treg_MKI67': 'MKI67+ Treg'
})

# Curve for lineage 1 only (the main CCR8→MKI67 curve)
curve1 = curves[curves['curve'] == 1].copy()

# ---------------------------------------------------------------------------
# Helper: save figure
# ---------------------------------------------------------------------------
def save_fig(fig, name):
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(OUT_DIR, f"{name}.{ext}"), dpi=DPI, bbox_inches='tight')
    print(f"  Saved: {name}")

# ===========================================================================
# 6A — UMAP colored by subtype with Slingshot curve overlay
# ===========================================================================
def panel_6a():
    fig, ax = plt.subplots(figsize=(4.2, 4))

    # Scatter: colored by subtype
    for st, color, label in [('CD4T_Treg_CCR8', CCR8_COLOR, 'CCR8+ Treg'),
                              ('CD4T_Treg_MKI67', MKI67_COLOR, 'MKI67+ Treg')]:
        sub = main[main['sub_cell_type'] == st]
        ax.scatter(sub['UMAP_1'], sub['UMAP_2'], c=color, s=0.4, alpha=0.4,
                   rasterized=True, label=label)

    # Slingshot curve overlay (Lineage 1 — main trajectory)
    ax.plot(curve1['x'], curve1['y'], color='black', linewidth=1.2, linestyle='-', alpha=0.8)

    # Arrow to indicate direction
    mid_idx = len(curve1) // 2
    dx = curve1['x'].iloc[mid_idx + 5] - curve1['x'].iloc[mid_idx]
    dy = curve1['y'].iloc[mid_idx + 5] - curve1['y'].iloc[mid_idx]
    norm = np.sqrt(dx**2 + dy**2)
    if norm > 0:
        ax.annotate('', xy=(curve1['x'].iloc[-10], curve1['y'].iloc[-10]),
                    xytext=(curve1['x'].iloc[-30], curve1['y'].iloc[-30]),
                    arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

    ax.set_xlabel('UMAP 1', fontsize=10)
    ax.set_ylabel('UMAP 2', fontsize=10)
    ax.set_title('Trajectory by Subtype', fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend
    legend_elements = [
        mpatches.Patch(color=CCR8_COLOR, label='CCR8+ Treg'),
        mpatches.Patch(color=MKI67_COLOR, label='MKI67+ Treg'),
        Line2D([0], [0], color='black', linewidth=1.2, label='Slingshot trajectory'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, frameon=False)

    # Subtype counts annotation
    n_ccr8 = (main['sub_cell_type'] == 'CD4T_Treg_CCR8').sum()
    n_mki67 = (main['sub_cell_type'] == 'CD4T_Treg_MKI67').sum()
    ax.text(0.02, 0.02, f'n(CCR8+) = {n_ccr8:,}\nn(MKI67+) = {n_mki67:,}',
            transform=ax.transAxes, fontsize=7, va='bottom', ha='left',
            color='gray')

    save_fig(fig, 'fig6a_trajectory_subtype')
    plt.close(fig)

# ===========================================================================
# 6B — UMAP colored by pseudotime with curve overlay
# ===========================================================================
def panel_6b():
    fig, ax = plt.subplots(figsize=(4.2, 4))

    valid = main[main['pseudotime'].notna()]
    sc = ax.scatter(valid['UMAP_1'], valid['UMAP_2'], c=valid['pseudotime'],
                    cmap=PT_CMAP, s=0.4, alpha=0.5, rasterized=True, vmin=0, vmax=16)

    # Curve
    ax.plot(curve1['x'], curve1['y'], color='black', linewidth=1.2, alpha=0.8)

    ax.set_xlabel('UMAP 1', fontsize=10)
    ax.set_ylabel('UMAP 2', fontsize=10)
    ax.set_title('Pseudotime Trajectory', fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    cbar = plt.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label('Pseudotime', fontsize=9)

    fig.tight_layout()
    save_fig(fig, 'fig6b_trajectory_pseudotime')
    plt.close(fig)

# ===========================================================================
# 6C — Pseudotime density by subtype × response
# ===========================================================================
def panel_6c():
    fig, ax = plt.subplots(figsize=(4.2, 3.2))

    valid = main[main['pseudotime'].notna()].copy()
    valid['resp_group'] = valid['response_binary'].map({
        'Responder': 'Responder', 'Non-responder': 'Non-responder'
    })

    color_map = {
        ('CD4T_Treg_CCR8', 'Responder'): '#2E5AAC',
        ('CD4T_Treg_CCR8', 'Non-responder'): '#C0392B',
        ('CD4T_Treg_MKI67', 'Responder'): '#5B9BD5',
        ('CD4T_Treg_MKI67', 'Non-responder'): '#E74C3C',
    }

    for (st, resp), grp in valid.groupby(['sub_cell_type', 'response_binary']):
        color = color_map.get((st, resp), 'gray')
        label_short = st.replace('CD4T_Treg_', '') + ' ' + resp
        # Use stat in kde to show density
        grp['pseudotime'].plot.kde(ax=ax, color=color, linewidth=1.5,
                                    label=label_short, alpha=0.8)

    ax.set_xlabel('Pseudotime', fontsize=10)
    ax.set_ylabel('Density', fontsize=10)
    ax.set_title('Pseudotime Distribution by Group', fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=7, frameon=False, loc='upper left', bbox_to_anchor=(0.02, 0.98))

    fig.tight_layout()
    save_fig(fig, 'fig6c_pseudotime_density')
    plt.close(fig)

# ===========================================================================
# 6D — Key gene LOESS trends along pseudotime
# ===========================================================================
def panel_6d():
    """4-panel grid: CCR8, MKI67, CTLA4, CD74 expression trends."""
    plot_genes = ['CCR8', 'MKI67', 'CTLA4', 'CD74']
    gene_labels = ['CCR8', 'MKI67', 'CTLA4', 'CD74 (MIF receptor)']

    fig, axes = plt.subplots(2, 2, figsize=(4.5, 3.5))
    fig.suptitle('Gene Expression Along Pseudotime', fontsize=10, fontweight='bold', y=1.02)

    axes_flat = axes.flatten() if hasattr(axes, 'flatten') else axes
    for idx, (gene, label) in enumerate(zip(plot_genes, gene_labels)):
        ax = axes_flat[idx]
        sub = gexpr[gexpr['gene'] == gene].copy()

        for resp_grp, color, resp_label in [
            ('Responder', RESP_COLOR, 'Responder'),
            ('Non-responder', NONRESP_COLOR, 'Non-responder')
        ]:
            grp = sub[sub['response_binary'] == resp_grp].dropna(subset=['pseudotime', 'expression'])
            if len(grp) < 50:
                continue

            # Sort by pseudotime
            grp = grp.sort_values('pseudotime')

            # Binned LOESS-style: rolling mean with ±1 SEM
            bins = np.linspace(0, 16, 33)
            grp['pt_bin'] = pd.cut(grp['pseudotime'], bins=bins)
            binned = grp.groupby('pt_bin', observed=False).agg(
                mean_expr=('expression', 'mean'),
                sem_expr=('expression', lambda x: x.std() / np.sqrt(len(x))),
                n=('expression', 'count')
            ).reset_index()
            binned['pt_mid'] = [(interval.left + interval.right) / 2 for interval in binned['pt_bin'].cat.categories]

            # Filter bins with few cells
            binned = binned[binned['n'] >= 10]

            if len(binned) < 3:
                continue

            ax.plot(binned['pt_mid'], binned['mean_expr'], color=color, linewidth=1.5,
                    label=resp_label, alpha=0.9)
            ax.fill_between(binned['pt_mid'],
                            binned['mean_expr'] - binned['sem_expr'],
                            binned['mean_expr'] + binned['sem_expr'],
                            color=color, alpha=0.1)

        ax.set_xlabel('Pseudotime', fontsize=9)
        ax.set_ylabel('log10(CPM+1)', fontsize=9)
        ax.set_title(label, fontsize=10, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(0, 16)
        if idx == 0:
            ax.legend(fontsize=7, frameon=False, loc='upper right')

    fig.tight_layout(rect=[0, 0, 1, 0.98])
    save_fig(fig, 'fig6d_gene_trends')
    plt.close(fig)

# ===========================================================================
# Combined layout
# ===========================================================================
def combined_figure():
    """Assemble 6A+6B+6C+6D into a 2x2 multi-panel figure using matplotlib add_axes.
    Preserves original figure proportions and natural whitespace."""
    from PIL import Image
    import numpy as np

    panels = {
        'A': 'fig6a_trajectory_subtype.png',
        'B': 'fig6b_trajectory_pseudotime.png',
        'C': 'fig6c_pseudotime_density.png',
        'D': 'fig6d_gene_trends.png',
    }

    # Load images
    imgs = {}
    aspects = {}
    for label, filename in panels.items():
        img_path = os.path.join(OUT_DIR, filename)
        if os.path.exists(img_path):
            img = Image.open(img_path)
            imgs[label] = np.array(img)
            aspects[label] = img.width / img.height
        else:
            raise FileNotFoundError(f"Panel image not found: {img_path}")

    # Layout: 2x2 grid with natural spacing
    # Compute display sizes that fit within a reasonable figure width
    TARGET_WIDTH_IN = 8.0  # inches (fits JTM full-width ~2400px @ 300 DPI)
    GAP_IN = 0.15          # inches between panels
    MARGIN_IN = 0.25       # inches outer margin

    avail_w = TARGET_WIDTH_IN - 2 * MARGIN_IN

    # Column widths: proportional to average aspect ratio
    # Left column (A, C): fit to narrower aspect
    # Right column (B, D): fit to wider aspect
    w_left = avail_w * 0.48
    w_right = avail_w * 0.52

    # Display heights
    h_a = w_left / aspects['A']
    h_b = w_right / aspects['B']
    h_c = w_left / aspects['C']
    h_d = w_right / aspects['D']

    # Row heights
    h_row0 = max(h_a, h_b)
    h_row1 = max(h_c, h_d)

    fig_h = 2 * MARGIN_IN + h_row0 + h_row1 + GAP_IN

    # Compute positions (left, bottom, width, height) in figure coordinates
    left_a = MARGIN_IN / TARGET_WIDTH_IN
    left_b = (MARGIN_IN + w_left + GAP_IN) / TARGET_WIDTH_IN
    left_c = MARGIN_IN / TARGET_WIDTH_IN
    left_d = (MARGIN_IN + w_left + GAP_IN) / TARGET_WIDTH_IN

    width_a = w_left / TARGET_WIDTH_IN
    width_b = w_right / TARGET_WIDTH_IN
    width_c = w_left / TARGET_WIDTH_IN
    width_d = w_right / TARGET_WIDTH_IN

    bottom_a = (MARGIN_IN + h_row1 + GAP_IN) / fig_h
    bottom_b = (MARGIN_IN + h_row1 + GAP_IN) / fig_h
    bottom_c = MARGIN_IN / fig_h
    bottom_d = MARGIN_IN / fig_h

    height_a = h_a / fig_h
    height_b = h_b / fig_h
    height_c = h_c / fig_h
    height_d = h_d / fig_h

    fig = plt.figure(figsize=(TARGET_WIDTH_IN, fig_h))

    positions = {
        'A': [left_a, bottom_a, width_a, height_a],
        'B': [left_b, bottom_b, width_b, height_b],
        'C': [left_c, bottom_c, width_c, height_c],
        'D': [left_d, bottom_d, width_d, height_d],
    }

    for label, pos in positions.items():
        ax = fig.add_axes(pos)
        ax.imshow(imgs[label])
        ax.axis('off')
        # Panel label
        ax.text(0.015, 0.98, label, transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='top', ha='left',
                color='black',
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                         alpha=0.9, edgecolor='none', linewidth=0))

    # Save
    out_png = os.path.join(OUT_DIR, 'fig6_combined.png')
    out_pdf = os.path.join(OUT_DIR, 'fig6_combined.pdf')
    fig.savefig(out_png, dpi=DPI, bbox_inches='tight', pad_inches=0.02)
    fig.savefig(out_pdf, dpi=DPI, bbox_inches='tight', pad_inches=0.02)
    print("  Saved: fig6_combined")
    plt.close(fig)

# ===========================================================================
# Main
# ===========================================================================
if __name__ == '__main__':
    print("Generating Figure 6 panels...")
    panel_6a()
    panel_6b()
    panel_6c()
    panel_6d()
    # Panel 6E is generated by fig6e_pseudotime_heatmap.py
    combined_figure()
    print("\nDone! All figures in:", OUT_DIR)
