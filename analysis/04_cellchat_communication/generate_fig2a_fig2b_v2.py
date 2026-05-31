#!/usr/bin/env python3
"""
Fig 2A + Fig 2B v2 — Two-panel layout with side-by-side comparison
- Fig 2A: Weak (left) vs Strong (right) source×target communication heatmaps
- Fig 2B: Grouped bar chart with all 5 KO conditions, both Weak and Strong shown
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import matplotlib.patheffects as pe

BASE_DIR = r'D:\Research\tomato'
DECENT_DIR = r'D:\Research\decent'
FIGURES_DIR = os.path.join(BASE_DIR, 'figures', 'demo', '2')
os.makedirs(FIGURES_DIR, exist_ok=True)

# Colors
COLOR_RESPONDER = '#8bc98b'
COLOR_NONRESPONDER = '#b08ebf'
COLOR_MID = '#ffffff'

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 12,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
    'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
})

def savefig_both(fig, name):
    fig.savefig(os.path.join(FIGURES_DIR, f'{name}.png'), dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(FIGURES_DIR, f'{name}.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved {name}.png/pdf")

# ==============================================================================
# FIG 2A: Side-by-side Weak vs Strong communication heatmaps
# ==============================================================================
print("Generating Fig 2A...")

df_weak = pd.read_csv(os.path.join(DECENT_DIR, 'cellchat', 'interactions_weak.csv'))
df_strong = pd.read_csv(os.path.join(DECENT_DIR, 'cellchat', 'interactions_strong.csv'))

def build_matrix(df):
    st = df.groupby(['source', 'target'])['prob'].sum().reset_index()
    return st.pivot(index='source', columns='target', values='prob')

w_mat = build_matrix(df_weak)
s_mat = build_matrix(df_strong)

cell_types = ['Fibroblast', 'Epithelial', 'Endothelial',
              'Myeloid', 'DC', 'Mast',
              'NK', 'CD8_T', 'CD4_T', 'T_cells', 'Treg', 'B_cells', 'Plasma']

w_aligned = w_mat.reindex(index=cell_types, columns=cell_types).fillna(0)
s_aligned = s_mat.reindex(index=cell_types, columns=cell_types).fillna(0)

# Shared log scale for both panels
prob_max = max(w_aligned.values.max(), s_aligned.values.max())
prob_min = min(w_aligned.values.min(), s_aligned.values.min())
prob_min = max(prob_min, 1e-12)

fig, axes = plt.subplots(1, 3, figsize=(14, 6.5),
                         gridspec_kw={'width_ratios': [1, 1, 0.05]})

# Purple sequential colormap for communication probability
colors_seq = ['#f7f4f9', '#e7e1ef', '#d4b9da', '#c994c7', '#df65b0', '#e7298a', '#ce1256', '#91003f']
cmap_seq = LinearSegmentedColormap.from_list('seq_purp', colors_seq)

for ax, mat, title, label_color in [
    (axes[0], w_aligned, 'Weak Responder', COLOR_NONRESPONDER),
    (axes[1], s_aligned, 'Strong Responder', COLOR_RESPONDER)
]:
    im = ax.imshow(np.log10(mat.values + 1e-12), cmap=cmap_seq,
                   vmin=np.log10(1e-12), vmax=np.log10(prob_max), aspect='auto')
    
    ax.set_xticks(np.arange(len(cell_types)))
    ax.set_yticks(np.arange(len(cell_types)))
    ax.set_xticklabels(cell_types, rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(cell_types, fontsize=9)
    
    # Bold highlight for Fibroblast and Treg
    for i, lbl in enumerate(ax.get_yticklabels()):
        if cell_types[i] == 'Fibroblast':
            lbl.set_fontweight('bold'); lbl.set_color('#222222')
    for i, lbl in enumerate(ax.get_xticklabels()):
        if cell_types[i] == 'Treg':
            lbl.set_fontweight('bold'); lbl.set_color('#222222')
    
    ax.set_xlabel('Target', fontsize=11)
    if ax == axes[0]:
        ax.set_ylabel('Source', fontsize=11)
    
    # Title with colored background strip
    ax.set_title(title, fontsize=12, fontweight='bold', color=label_color, pad=8)
    
    # Minor grid
    ax.set_xticks(np.arange(-0.5, len(cell_types), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(cell_types), 1), minor=True)
    ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.3, alpha=0.4)
    
    # Box around Fibroblast->Treg
    fib_idx = cell_types.index('Fibroblast')
    treg_idx = cell_types.index('Treg')
    rect = Rectangle((treg_idx-0.45, fib_idx-0.45), 0.9, 0.9,
                     fill=False, edgecolor='black', linewidth=2)
    ax.add_patch(rect)

# Shared colorbar
cbar = plt.colorbar(im, cax=axes[2])
cbar.set_label('log10(Communication Probability)', fontsize=10)
cbar.ax.tick_params(labelsize=8)

plt.suptitle('A. CellChat Communication Landscape', fontsize=13, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])
savefig_both(fig, 'fig2a_cellchat_sidebyside')
plt.close()

# ==============================================================================
# FIG 2B: KO grouped bar chart with all conditions + both groups
# ==============================================================================
print("Generating Fig 2B...")

ko = pd.read_csv(os.path.join(BASE_DIR, 'results', 'cellchat_ko_corrected', 'ko_summary.csv'))

# All 5 KO conditions, ordered by severity
show_kos = ['KO_COL1A1', 'KO_COL1A2', 'KO_FN1', 'KO_COL1A1_COL1A2', 'KO_ALL_ECM']
ko_labels = {
    'KO_COL1A1': 'COL1A1',
    'KO_COL1A2': 'COL1A2',
    'KO_FN1': 'FN1',
    'KO_COL1A1_COL1A2': 'COL1A1+2',
    'KO_ALL_ECM': 'All ECM',
}

ko_plot = ko[ko['KO'].isin(show_kos)].copy()
ko_plot['Label'] = ko_plot['KO'].map(ko_labels)
ko_plot['order'] = ko_plot['KO'].map({k: i for i, k in enumerate(show_kos)})
ko_plot = ko_plot.sort_values('order')

fig, ax = plt.subplots(figsize=(6.5, 5))

x = np.arange(len(show_kos))
width = 0.35

# Weak = Signal Loss % (positive, meaningful)
weak_df = ko_plot[ko_plot['Group']=='Weak'].set_index('KO').loc[show_kos]
weak_vals = weak_df['Lost_Pct'].clip(lower=0).values
bars_w = ax.bar(x - width/2, weak_vals, width, label='Weak (Signal Lost %)',
                color=COLOR_NONRESPONDER, edgecolor='black', linewidth=0.6, zorder=3)

# Strong = absolute KO prob (right y-axis) — more honest than "loss %" for a group with no ECM baseline
ax2 = ax.twinx()
strong_df = ko_plot[ko_plot['Group']=='Strong'].set_index('KO').loc[show_kos]
strong_vals = strong_df['KO_Prob'].values * 1e6  # scale to ×10^-6
bars_s = ax2.bar(x + width/2, strong_vals, width, label='Strong (KO Prob ×10⁻⁶)',
                 color=COLOR_RESPONDER, edgecolor='black', linewidth=0.6, alpha=0.7, zorder=3)

# Labels on bars
for bar, val in zip(bars_w, weak_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
            f'{val:.1f}%', ha='center', va='bottom', fontsize=8, fontweight='bold')

for bar, val in zip(bars_s, strong_vals):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
             f'{val:.2f}', ha='center', va='bottom', fontsize=8, fontweight='bold', color='#2d5a2d')

# Axis settings
ax.set_ylabel('Weak: Signal Lost (%)', fontsize=11, color=COLOR_NONRESPONDER)
ax.set_ylim(0, 105)
ax.tick_params(axis='y', labelcolor=COLOR_NONRESPONDER)

ax2.set_ylabel('Strong: KO Prob (×10⁻⁶)', fontsize=11, color='#2d5a2d')
ax2.set_ylim(0, max(strong_vals) * 1.3)
ax2.tick_params(axis='y', labelcolor='#2d5a2d')

ax.set_xticks(x)
ax.set_xticklabels([ko_labels[k] for k in show_kos], fontsize=9)
ax.set_xlabel('KO Condition', fontsize=11)

# Title
ax.set_title('B. Virtual Knockout of ECM Ligands (Fibroblast → Treg)', 
             fontsize=12, fontweight='bold', pad=10)

# Combined legend
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left', frameon=False, fontsize=8.5)

# Annotation
ax.annotate('Strong group has negligible ECM baseline;\nvalues shown are post-KO absolute probabilities',
            xy=(0.98, 0.02), xycoords='axes fraction', fontsize=7.5,
            ha='right', va='bottom', color='#555555',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#f8f8f8', edgecolor='#cccccc'))

ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.yaxis.grid(True, linestyle='--', alpha=0.3, zorder=0)
ax.set_axisbelow(True)

plt.tight_layout()
savefig_both(fig, 'fig2b_virtual_ko_with_strong')
plt.close()

print("\nDone! Both panels generated.")
