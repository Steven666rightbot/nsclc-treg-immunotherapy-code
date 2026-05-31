#!/usr/bin/env python3
"""
Fig 2B — Virtual Knockout Signal Loss (CORRECTED annotation)
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 12,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
    'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
})

OUT_DIR = r"D:\Research\tomato\figures\demo\2"
os.makedirs(OUT_DIR, exist_ok=True)

# Load corrected KO summary
ko = pd.read_csv(r"D:\Research\tomato\results\cellchat_ko_corrected\ko_summary.csv")

# Filter to the 4 conditions we want to show (exclude single COL1A1 since COL1A1+COL1A2 covers it better)
# Actually, show: COL1A2, FN1, COL1A1+COL1A2, ALL_ECM
show_kos = ['KO_COL1A2', 'KO_FN1', 'KO_COL1A1_COL1A2', 'KO_ALL_ECM']
ko_plot = ko[ko['KO'].isin(show_kos)].copy()

# Clip negative loss to 0 for visualization (Strong group has no ECM baseline)
ko_plot['Lost_Pct_clip'] = ko_plot['Lost_Pct'].clip(lower=0)

# Map KO names to display labels
label_map = {
    'KO_COL1A2': 'COL1A2',
    'KO_FN1': 'FN1',
    'KO_COL1A1_COL1A2': 'COL1A1+COL1A2',
    'KO_ALL_ECM': 'All ECM',
}
ko_plot['Label'] = ko_plot['KO'].map(label_map)

# Colors
COLOR_WEAK = '#b08ebf'   # Purple = non-responder / weak
COLOR_STRONG = '#8bc98b' # Green = responder / strong

fig, ax = plt.subplots(figsize=(5, 4.5))

x = np.arange(len(show_kos))
width = 0.35

for i, grp in enumerate(['Weak', 'Strong']):
    df_g = ko_plot[ko_plot['Group'] == grp].set_index('KO').loc[show_kos]
    vals = df_g['Lost_Pct_clip'].values
    colors = [COLOR_WEAK if grp == 'Weak' else COLOR_STRONG] * len(vals)
    bars = ax.bar(x + (i - 0.5) * width, vals, width, label=grp, color=colors, 
                  edgecolor='black', linewidth=0.6, zorder=3)
    
    # Add value labels on bars
    for bar, val, raw_val in zip(bars, vals, df_g['Lost_Pct'].values):
        if raw_val < 0:
            text = 'n/a*'
            ypos = 2
        else:
            text = f'{val:.1f}%'
            ypos = bar.get_height() + 1.5
        ax.text(bar.get_x() + bar.get_width()/2, ypos, text,
                ha='center', va='bottom', fontsize=8, fontweight='bold')

ax.set_ylabel('Signal Loss (%)', fontsize=11)
ax.set_xticks(x)
ax.set_xticklabels([label_map[k] for k in show_kos], fontsize=9)
ax.set_ylim(0, 105)
ax.set_title('B. Virtual Knockout of ECM Ligands\n(Fibroblast → Treg)', 
             fontsize=12, fontweight='bold')

# Legend
weak_patch = mpatches.Patch(facecolor=COLOR_WEAK, edgecolor='black', linewidth=0.6, label='Weak responder')
strong_patch = mpatches.Patch(facecolor=COLOR_STRONG, edgecolor='black', linewidth=0.6, label='Strong responder')
ax.legend(handles=[weak_patch, strong_patch], loc='upper right', frameon=False)

# Annotation for Strong n/a
ax.annotate('*Strong group has no baseline ECM signal\n(MHC-II only), KO effect not applicable',
            xy=(0.02, 0.02), xycoords='axes fraction', fontsize=7.5,
            ha='left', va='bottom', color='#555555',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#f8f8f8', edgecolor='#cccccc'))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.yaxis.grid(True, linestyle='--', alpha=0.4, zorder=0)
ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'fig2b_virtual_ko_signal_loss.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUT_DIR, 'fig2b_virtual_ko_signal_loss.pdf'), dpi=300, bbox_inches='tight')
plt.close()

print(f"Saved: {OUT_DIR}/fig2b_virtual_ko_signal_loss.png")
print(f"Saved: {OUT_DIR}/fig2b_virtual_ko_signal_loss.pdf")

# Print data table for verification
print("\nData used:")
print(ko_plot[['KO','Group','Lost_Pct','Lost_Pct_clip']].to_string(index=False))
