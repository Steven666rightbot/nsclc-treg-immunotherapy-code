"""
Figure 5 (Redesigned): NicheNet-based validation of MHC-I signaling targets.
Panel A: Pathway-grouped dotplot of MHC-I predicted targets
Panel B: Pathway enrichment significance
"""
import pyreadr
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from matplotlib.patches import Patch
import os

OUTDIR = "D:/research/cucumber/fig5"
os.makedirs(OUTDIR, exist_ok=True)

COLORS = {
    'Signaling': '#2E5AAC',
    'MIF Pathway': '#C0392B',
    'HSF1 Stress': '#C0392B',
}

# Load matrix
ltm = pyreadr.read_r('D:/research/tomato/data/nichenet_ligand_target_matrix_uncompressed.rds')[None]
mhci_cols = [c for c in ['HLA-A', 'HLA-B', 'HLA-C'] if c in ltm.columns]
mhci_scores = ltm[mhci_cols].mean(axis=1).sort_values(ascending=False)

# Gene → pathway mapping
PATHWAYS = {
    'Signaling': ['JUN', 'FOS', 'NFKB1', 'RELA'],
    'MIF Pathway': ['CD44', 'CD74', 'CXCR4'],
    'HSF1 Stress': ['HSPA1A', 'HSPA1B', 'HSPB1', 'HSP90AA1', 'DNAJB1', 'BAG3'],
}

# Gene sets for enrichment (same as before)
GENE_SETS = {
    'Signaling': ['JUN', 'FOS', 'NFKB1', 'RELA', 'STAT2'],
    'MIF Pathway': ['CD44', 'CD74', 'CXCR4', 'MIF'],
    'HSF1 Stress': ['HSPA1A', 'HSPA1B', 'HSP90AA1', 'DNAJB1', 'BAG3', 'HSPB1', 'HSPH1'],
    'MHC-II Antigen': ['HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DPA1', 'HLA-DQA1', 'HLA-DQB1'],
    'Treg Function': ['FOXP3', 'IL2RA', 'CTLA4', 'TIGIT', 'ICOS', 'TNFRSF18'],
}

# ============================================================
# Panel A: Pathway-arranged dotplot
# ============================================================
# Build the plot data
plot_genes = []
for pw_name, genes in PATHWAYS.items():
    for g in genes:
        if g in mhci_scores.index:
            rank = list(mhci_scores.index).index(g) + 1
            plot_genes.append({
                'gene': g,
                'pathway': pw_name,
                'score': mhci_scores[g],
                'rank': rank,
                'color': COLORS[pw_name],
            })

plot_df = pd.DataFrame(plot_genes)
plot_df = plot_df.sort_values(['pathway', 'score'], ascending=[True, False])

# Also add top non-pathway genes as background
non_pathway_genes = mhci_scores.head(10).index.tolist()
all_pathway_genes = set(g for genes in PATHWAYS.values() for g in genes)
non_pathway_bg = [g for g in non_pathway_genes if g not in all_pathway_genes]

bg_points = []
for g in non_pathway_bg:
    bg_points.append({
        'gene': g,
        'pathway': 'Other',
        'score': mhci_scores[g],
        'rank': list(mhci_scores.index).index(g) + 1,
        'color': '#BDC3C7',
    })

fig_a, ax_a = plt.subplots(figsize=(7, 4))

# Draw background points
if bg_points:
    bg_df = pd.DataFrame(bg_points)
    ax_a.scatter(bg_df['score'], [0.5] * len(bg_df), s=40, c=bg_df['color'],
                edgecolors='#999999', linewidth=0.3, alpha=0.6, zorder=1)
    for _, row in bg_df.iterrows():
        ax_a.annotate(row['gene'], (row['score'], 0.5),
                     fontsize=7, color='#999999', ha='center', va='bottom',
                     rotation=45, style='italic')

# Draw pathway genes - each pathway on its own row
for i, pw in enumerate(['Signaling', 'MIF Pathway', 'HSF1 Stress']):
    genes_in_pw = plot_df[plot_df['pathway'] == pw]
    y_pos = i  # 0, 1, 2
    
    # Individual genes as scatter points
    for _, row in genes_in_pw.iterrows():
        ax_a.scatter(row['score'], y_pos, s=120, c=row['color'],
                    edgecolors='black', linewidth=0.5, zorder=3)
        # Gene label
        ax_a.annotate(f"{row['gene']} (rank {row['rank']})", 
                     (row['score'], y_pos),
                     fontsize=8, ha='left', va='center',
                     xytext=(5, 0), textcoords='offset points',
                     fontweight='bold', color=row['color'])

# Pathway labels on y-axis
ax_a.set_yticks([0, 1, 2])
ax_a.set_yticklabels(['Signaling\n(JUN/FOS/NFKB1)', 'MIF Pathway\n(CD44/CD74/CXCR4)', 'HSF1 Stress\n(HSPA1A/B)'], 
                     fontsize=9, fontweight='bold')
ax_a.set_xlabel('NicheNet Regulatory Potential', fontsize=10)

# Vertical reference lines at key ranks
# Add a second x-axis showing percentile
ax_a.set_title('A  MHC-I (HLA-A/B/C): Predicted Target Genes by Pathway', 
              fontsize=11, fontweight='bold', ha='left', loc='left')

ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.set_xlim(0, mhci_scores.head(10).max() * 1.3)

# Add rank labels on top
ax_a_twin = ax_a.twiny()
ax_a_twin.set_xlim(ax_a.get_xlim())
# Map score to rank
score_at_rank = {mhci_scores.iloc[r]: r+1 for r in [0, 99, 499, 999]}
rank_positions = list(score_at_rank.keys())
rank_labels = [f'rank {score_at_rank[s]}' for s in rank_positions]
# Only keep a few
rank_ticks = []
rank_tick_labels = []
thresholds = [mhci_scores.iloc[0], mhci_scores.iloc[99], mhci_scores.iloc[499]]
threshold_labels = ['#1', 'top 100', 'top 500']
# Find these on the axis
for t, l in zip(thresholds, threshold_labels):
    if t > ax_a.get_xlim()[0] and t < ax_a.get_xlim()[1]:
        rank_ticks.append(t)
        rank_tick_labels.append(l)
rank_ticks.append(mhci_scores.iloc[999])
rank_tick_labels.append('top 1000')

ax_a_twin.set_xticks(rank_ticks)
ax_a_twin.set_xticklabels(rank_tick_labels, fontsize=7, color='#777777')
ax_a_twin.spines['top'].set_visible(False)

# Legend
legend_elements = [
    Patch(facecolor=COLORS['Signaling'], label='Signaling module'),
    Patch(facecolor=COLORS['MIF Pathway'], label='MIF pathway'),
    Patch(facecolor=COLORS['HSF1 Stress'], label='HSF1 stress program'),
    Patch(facecolor='#BDC3C7', label='Other (top 10)'),
]
ax_a.legend(handles=legend_elements, loc='lower right', fontsize=8, frameon=True,
           framealpha=0.9, edgecolor='#cccccc')

fig_a.tight_layout()
fig_a.savefig(f"{OUTDIR}/fig5a_mhci_targets.png", dpi=300, bbox_inches='tight')
fig_a.savefig(f"{OUTDIR}/fig5a_mhci_targets.pdf", dpi=300, bbox_inches='tight')
plt.close(fig_a)
print("Panel A saved")

# ============================================================
# Panel B: Pathway Enrichment (keep same design, maybe cleaner)
# ============================================================
all_genes = list(mhci_scores.index)
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

fig_b, ax_b = plt.subplots(figsize=(6, 3.5))

colors_b = ['#2E5AAC' if r['significant'] else '#CCCCCC' for _, r in enrich_df.iterrows()]
bars = ax_b.barh(range(len(enrich_df)), enrich_df['neg_log10_p'].values,
                 color=colors_b, height=0.55, edgecolor='white', linewidth=0.3)
ax_b.set_yticks(range(len(enrich_df)))
ax_b.set_yticklabels(enrich_df['pathway'].values, fontsize=8.5)
ax_b.invert_yaxis()
ax_b.set_xlabel('−log₁₀(p-value)', fontsize=10)
ax_b.set_title('B  Pathway Enrichment in Top 500 MHC-I Targets', 
              fontsize=11, fontweight='bold', ha='left', loc='left')
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

# Annotate
for i, (_, row) in enumerate(enrich_df.iterrows()):
    val = row['neg_log10_p']
    if row['significant']:
        ax_b.text(val + 0.05, i, f'p = {row["pval"]:.4f}\n{row["overlap"]}/{row["n_in_set"]} genes',
                 fontsize=7.5, va='center', color='#333333')
    else:
        ax_b.text(val + 0.05, i, f'ns',
                 fontsize=7.5, va='center', color='#999999')

ax_b.axvline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.6, alpha=0.5)
ax_b.text(-np.log10(0.05) + 0.02, len(enrich_df) - 0.3, 'p = 0.05', 
         fontsize=7, color='gray', style='italic')

fig_b.tight_layout()
fig_b.savefig(f"{OUTDIR}/fig5b_pathway_enrichment.png", dpi=300, bbox_inches='tight')
fig_b.savefig(f"{OUTDIR}/fig5b_pathway_enrichment.pdf", dpi=300, bbox_inches='tight')
plt.close(fig_b)
print("Panel B saved")

# Combined
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw={'width_ratios': [1.3, 1]})

# Recreate Panel A in combined
ax = axes[0]
# BG points
if bg_points:
    bg_df = pd.DataFrame(bg_points)
    ax.scatter(bg_df['score'], [0.5] * len(bg_df), s=50, c=bg_df['color'],
              edgecolors='#999999', linewidth=0.3, alpha=0.5, zorder=1)
    for _, row in bg_df.iterrows():
        ax.annotate(row['gene'], (row['score'], 0.5),
                   fontsize=6.5, color='#999999', ha='center', va='bottom',
                   rotation=45, style='italic')

for i, pw in enumerate(['Signaling', 'MIF Pathway', 'HSF1 Stress']):
    genes_in_pw = plot_df[plot_df['pathway'] == pw]
    y_pos = i
    for _, row in genes_in_pw.iterrows():
        ax.scatter(row['score'], y_pos, s=130, c=row['color'],
                  edgecolors='black', linewidth=0.5, zorder=3)
        ax.annotate(f"{row['gene']} (rank {row['rank']})", 
                   (row['score'], y_pos),
                   fontsize=7.5, ha='left', va='center',
                   xytext=(5, 0), textcoords='offset points',
                   fontweight='bold', color=row['color'])

ax.set_yticks([0, 1, 2])
ax.set_yticklabels(['Signaling\n(JUN/FOS/NFKB1)', 'MIF Pathway\n(CD44/CD74/CXCR4)', 'HSF1 Stress\n(HSPA1A/B)'], 
                   fontsize=8.5, fontweight='bold')
ax.set_xlabel('Regulatory Potential', fontsize=9)
ax.set_title('A  MHC-I Predicted Targets by Pathway', fontsize=10, fontweight='bold', ha='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim(0, mhci_scores.head(10).max() * 1.25)

# Recreate Panel B in combined
ax = axes[1]
colors_b = ['#2E5AAC' if r['significant'] else '#CCCCCC' for _, r in enrich_df.iterrows()]
ax.barh(range(len(enrich_df)), enrich_df['neg_log10_p'].values,
        color=colors_b, height=0.55, edgecolor='white', linewidth=0.3)
ax.set_yticks(range(len(enrich_df)))
ax.set_yticklabels(enrich_df['pathway'].values, fontsize=8.5)
ax.invert_yaxis()
ax.set_xlabel('−log₁₀(p)', fontsize=9)
ax.set_title('B  Enrichment in Top 500', fontsize=10, fontweight='bold', ha='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for i, (_, row) in enumerate(enrich_df.iterrows()):
    val = row['neg_log10_p']
    if row['significant']:
        ax.text(val + 0.03, i, f'p={row["pval"]:.4f}', fontsize=7, va='center', color='#333333')

ax.axvline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.6, alpha=0.5)

fig.tight_layout()
fig.savefig(f"{OUTDIR}/fig5_combined_nichenet.png", dpi=300, bbox_inches='tight')
fig.savefig(f"{OUTDIR}/fig5_combined_nichenet.pdf", dpi=300, bbox_inches='tight')
plt.close(fig)
print("Combined saved")

# Copy to fig5_combined for pipeline
import shutil
shutil.copy(f"{OUTDIR}/fig5_combined_nichenet.png", f"{OUTDIR}/fig5_combined.png")
shutil.copy(f"{OUTDIR}/fig5_combined_nichenet.pdf", f"{OUTDIR}/fig5_combined.pdf")

print("\nDONE!")
