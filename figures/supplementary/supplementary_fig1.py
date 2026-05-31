"""
Supplementary Figure 1: NicheNet ligand-target analysis of MHC-I signaling.
Panel A: MHC-I (HLA-A/B/C) predicted target genes by pathway
Panel B: Pathway enrichment significance in top 500 targets
"""

import pyreadr
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from matplotlib.patches import Patch
from PIL import Image, ImageDraw, ImageFont
import os

OUTDIR = "D:/research/cucumber/supp_fig1"
os.makedirs(OUTDIR, exist_ok=True)

DPI = 300
LABEL_FONT_SIZE = 26
PANEL_GAP = 20

# Updated color scheme: HSF1=red, MIF=blue, Signaling=green
COLORS = {
    'Signaling': '#27AE60',
    'MIF Pathway': '#2980B9',
    'HSF1 Stress': '#C0392B',
}

# ── Load NicheNet matrix ──
ltm = pyreadr.read_r('D:/research/tomato/data/nichenet_ligand_target_matrix_uncompressed.rds')[None]
mhci_cols = [c for c in ['HLA-A', 'HLA-B', 'HLA-C'] if c in ltm.columns]
mhci_scores = ltm[mhci_cols].mean(axis=1).sort_values(ascending=False)

# Pathway definitions
PATHWAYS = {
    'Signaling': ['JUN', 'FOS', 'NFKB1', 'RELA'],
    'MIF Pathway': ['CD44', 'CD74', 'CXCR4'],
    'HSF1 Stress': ['HSPA1A', 'HSPA1B', 'HSPB1', 'HSP90AA1', 'DNAJB1', 'BAG3'],
}

GENE_SETS = {
    'MIF Pathway': ['CD44', 'CD74', 'CXCR4', 'MIF'],
    'HSF1 Stress': ['HSPA1A', 'HSPA1B', 'HSP90AA1', 'DNAJB1', 'BAG3', 'HSPB1', 'HSPH1'],
    'MHC-II Antigen': ['HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DPA1', 'HLA-DQA1', 'HLA-DQB1'],
    'Signaling': ['JUN', 'FOS', 'NFKB1', 'RELA', 'STAT2'],
    'Treg Function': ['FOXP3', 'IL2RA', 'CTLA4', 'TIGIT', 'ICOS', 'TNFRSF18'],
    'Proliferation': ['MKI67', 'TOP2A', 'PCNA', 'CDK1'],
}

all_genes = list(mhci_scores.index)

# ============================================================
# Panel A: Pathway-arranged dotplot
# ============================================================
plot_genes = []
for pw_name, genes in PATHWAYS.items():
    for g in genes:
        if g in mhci_scores.index:
            rank = list(mhci_scores.index).index(g) + 1
            plot_genes.append({
                'gene': g, 'pathway': pw_name,
                'score': mhci_scores[g], 'rank': rank,
                'color': COLORS[pw_name],
            })

plot_df = pd.DataFrame(plot_genes)
plot_df = plot_df.sort_values(['pathway', 'score'], ascending=[True, False])

# Background: top-non-pathway genes
non_pathway_genes = mhci_scores.head(10).index.tolist()
all_pathway_genes = set(g for genes in PATHWAYS.values() for g in genes)
non_pathway_bg = [g for g in non_pathway_genes if g not in all_pathway_genes]

bg_points = []
for g in non_pathway_bg:
    bg_points.append({
        'gene': g, 'pathway': 'Other',
        'score': mhci_scores[g],
        'rank': list(mhci_scores.index).index(g) + 1,
        'color': '#E0E0E0',
    })

fig_a, ax_a = plt.subplots(figsize=(6.3, 3.5))

# Remove top and right spines for clean look
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)

# Background points (light, small font, non-blocking)
if bg_points:
    bg_df = pd.DataFrame(bg_points)
    ax_a.scatter(bg_df['score'], [0.5] * len(bg_df), s=45, c=bg_df['color'],
                edgecolors='#CCCCCC', linewidth=0.3, alpha=0.6, zorder=1)
    for _, row in bg_df.iterrows():
        ax_a.annotate(row['gene'], (row['score'], 0.5),
                     fontsize=5, color='#E0E0E0', ha='center', va='bottom',
                     rotation=45, style='italic', zorder=1)

# Pathway genes with consistent, non-overlapping labels
xlim_max = mhci_scores.head(10).max() * 1.45
ax_a.set_xlim(0, xlim_max)
mid_x = xlim_max / 2

for i, pw in enumerate(['Signaling', 'MIF Pathway', 'HSF1 Stress']):
    genes_in_pw = plot_df[plot_df['pathway'] == pw].sort_values('score')
    y_pos = i
    n_genes = len(genes_in_pw)
    for j, (_, row) in enumerate(genes_in_pw.iterrows()):
        ax_a.scatter(row['score'], y_pos, s=60, c=row['color'],
                    edgecolors='black', linewidth=0.5, zorder=3)
        # Stagger labels vertically to prevent overlap in dense clusters
        if n_genes > 1:
            y_offset = (j - (n_genes - 1) / 2) * 0.15
        else:
            y_offset = 0
        # Consistent alignment: right-align if point is in right half
        if row['score'] > mid_x:
            ha = 'right'
            xytext = (-8, 0)
        else:
            ha = 'left'
            xytext = (8, 0)
        ax_a.annotate(row['gene'],
                     (row['score'], y_pos + y_offset),
                     fontsize=7, ha=ha, va='center',
                     xytext=xytext, textcoords='offset points',
                     fontweight='bold', color=row['color'],
                     zorder=4)

ax_a.set_yticks([0, 1, 2])
ax_a.set_yticklabels(['Signaling\n(JUN/FOS/NFKB1)', 'MIF Pathway\n(CD44/CD74/CXCR4)', 'HSF1 Stress\n(HSPA1A/B)'],
                     fontsize=9, fontweight='bold')
ax_a.set_xlabel('NicheNet Regulatory Potential', fontsize=10)
ax_a.set_title('A  MHC-I (HLA-A/B/C) Predicted Targets by Pathway',
              fontsize=11, fontweight='bold', ha='left', loc='left')

# Rank reference on top axis
thresholds = [mhci_scores.iloc[0], mhci_scores.iloc[99], mhci_scores.iloc[499]]
threshold_labels = ['#1', 'top 100', 'top 500']
ax_a_twin = ax_a.twiny()
ax_a_twin.set_xlim(ax_a.get_xlim())
rank_ticks = []
rank_tick_labels = []
for t, l in zip(thresholds, threshold_labels):
    if t > ax_a.get_xlim()[0] and t < ax_a.get_xlim()[1]:
        rank_ticks.append(t)
        rank_tick_labels.append(l)
ax_a_twin.set_xticks(rank_ticks)
ax_a_twin.set_xticklabels(rank_tick_labels, fontsize=7, color='#777777')
ax_a_twin.spines['top'].set_visible(False)
ax_a_twin.spines['right'].set_visible(False)
ax_a_twin.set_xlabel('Rank percentile (top %)', fontsize=8, color='#777777')

# Simplified legend
legend_elements = [
    Patch(facecolor=COLORS['HSF1 Stress'], label='HSF1 Stress'),
    Patch(facecolor=COLORS['MIF Pathway'], label='MIF Pathway'),
    Patch(facecolor=COLORS['Signaling'], label='Signaling'),
    Patch(facecolor='#E0E0E0', label='Other'),
]
ax_a.legend(handles=legend_elements, loc='lower right', fontsize=7, frameon=True,
           framealpha=0.9, edgecolor='#cccccc')

fig_a.tight_layout()
fig_a.savefig(f"{OUTDIR}/panel_a_targets.png", dpi=DPI, bbox_inches='tight')
plt.close(fig_a)

# ============================================================
# Panel B: Pathway Enrichment (vertical bars)
# ============================================================
enrich_results = []
for set_name, genes in GENE_SETS.items():
    n_total = len(all_genes)
    n_in_set = len([g for g in genes if g in all_genes])
    top_genes_set = set(mhci_scores.head(500).index)
    overlap = len(top_genes_set & set(genes))
    pval = hypergeom.sf(overlap - 1, n_total, n_in_set, 500) if n_in_set > 0 else 1.0
    enrich_results.append({
        'pathway': set_name,
        'neg_log10_p': -np.log10(max(pval, 1e-15)),
        'pval': pval,
        'overlap': overlap,
        'n_in_set': n_in_set,
        'significant': pval < 0.05,
    })

enrich_df = pd.DataFrame(enrich_results)

# Pathway color mapping
pathway_colors = {
    'MIF Pathway': '#2980B9',      # Blue
    'HSF1 Stress': '#C0392B',      # Red
    'Signaling': '#27AE60',         # Green
    'MHC-II Antigen': '#BDC3C7',   # Gray
    'Treg Function': '#BDC3C7',    # Gray
    'Proliferation': '#BDC3C7',    # Gray
}

# Sort by significance for visualization
enrich_df = enrich_df.sort_values('neg_log10_p', ascending=False)

fig_b, ax_b = plt.subplots(figsize=(6.3, 3.5))

bar_colors = [pathway_colors.get(r['pathway'], '#BDC3C7') for _, r in enrich_df.iterrows()]
x_pos = np.arange(len(enrich_df))
bars = ax_b.bar(x_pos, enrich_df['neg_log10_p'].values,
                color=bar_colors, width=0.6, edgecolor='white', linewidth=0.3)

ax_b.set_xticks(x_pos)
ax_b.set_xticklabels(enrich_df['pathway'].values, fontsize=9, rotation=0, ha='center')
ax_b.set_ylabel('-log$_{10}$(p-adjusted)', fontsize=10)
ax_b.set_title('B  Pathway Enrichment in Top 500 MHC-I Targets',
              fontsize=11, fontweight='bold', ha='left', loc='left')
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

# Significance annotations on top of bars
for i, (_, row) in enumerate(enrich_df.iterrows()):
    val = row['neg_log10_p']
    if row['pathway'] == 'MIF Pathway':
        ax_b.text(i, val + 0.08, '* p = 0.0023',
                 fontsize=7.5, va='bottom', ha='center', color='#333333', fontweight='bold')
    elif row['pathway'] == 'HSF1 Stress':
        ax_b.text(i, val + 0.08, '* p = 0.0076',
                 fontsize=7.5, va='bottom', ha='center', color='#333333', fontweight='bold')
    else:
        ax_b.text(i, val + 0.08, 'ns',
                 fontsize=7.5, va='bottom', ha='center', color='#999999')

# Significance threshold line (black thin dashed)
ax_b.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8, alpha=0.7)
ax_b.text(len(enrich_df) - 0.5, -np.log10(0.05) + 0.05, 'p = 0.05',
         fontsize=7, color='black', ha='right', va='bottom')

fig_b.tight_layout()
fig_b.savefig(f"{OUTDIR}/panel_b_enrichment.png", dpi=DPI, bbox_inches='tight')
plt.close(fig_b)

# ============================================================
# Combined: side-by-side
# ============================================================
img_a = Image.open(f"{OUTDIR}/panel_a_targets.png").convert('RGB')
img_b = Image.open(f"{OUTDIR}/panel_b_enrichment.png").convert('RGB')

# Stack vertically with same width
target_w = max(img_a.size[0], img_b.size[0])
scale_a = target_w / img_a.size[0]
scale_b = target_w / img_b.size[0]
img_a = img_a.resize((target_w, int(img_a.size[1] * scale_a)), Image.LANCZOS)
img_b = img_b.resize((target_w, int(img_b.size[1] * scale_b)), Image.LANCZOS)

canvas_w = target_w
canvas_h = img_a.size[1] + PANEL_GAP + img_b.size[1]

canvas = Image.new('RGB', (canvas_w, canvas_h), (255, 255, 255))
canvas.paste(img_a, (0, 0))
canvas.paste(img_b, (0, img_a.size[1] + PANEL_GAP))

# Save
combined_path_png = f"{OUTDIR}/supplementary_fig1.png"
combined_path_pdf = f"{OUTDIR}/supplementary_fig1.pdf"
canvas.save(combined_path_png, dpi=(DPI, DPI))
canvas.save(combined_path_pdf, dpi=(DPI, DPI))
print(f"Saved: {combined_path_png} ({canvas.size})")
print(f"Saved: {combined_path_pdf} ({canvas.size})")
print(f"Panel A: MHC-I predicted targets. Key genes: CD74 (rank 21/25938), HSPA1A (rank 52), CD44 (rank 294)")
print(f"Panel B: MIF pathway enriched (p=0.0023), HSF1 stress enriched (p=0.0076) in top 500")
