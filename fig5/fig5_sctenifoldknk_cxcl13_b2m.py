"""
Figure 5: scTenifoldKnk virtual knockout analysis
Panels:
  5A: Volcano plot - CXCL13 KO in CCR8+ Treg
  5B: Volcano plot - B2M KO in MKI67+ Treg
  5C: Top 15 perturbed genes (combined KOs)
  5D: GO BP enrichment for B2M KO (MKI67+)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import sys; sys.path.insert(0, 'D:/research/cucumber')
from _global_config import apply_style; apply_style()

OUTDIR = 'D:/research/cucumber/fig5'
os.makedirs(OUTDIR, exist_ok=True)

# KO-specific colors (not MPR/NR — these are different KO targets)
COLOR_CXCL13 = '#D55E00'   # Bang Wong vermillion (colorblind-safe red alt)
COLOR_B2M = '#0072B2'      # Bang Wong blue (exact)
COLOR_NS = '#BBBBBB'       # non-significant grey

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
df_cxcl13 = pd.read_csv(
    'D:/research/tomato/results/scTenifoldKnk_ko_ccr8_cxcl13/ko_CCR8_CXCL13.csv'
)
df_b2m = pd.read_csv(
    'D:/research/tomato/results/scTenifoldKnk_ko_mki67_b2m/ko_MKI67_B2M_lowres.csv'
)

# Guard against zero p.adj for log transform
for df in (df_cxcl13, df_b2m):
    df['log10_padj'] = -np.log10(df['p.adj'].replace(0, np.nextafter(0, 1)))
    df['absZ'] = df['Z'].abs()

# ------------------------------------------------------------------
# Stats to stdout
# ------------------------------------------------------------------
sig_cxcl13 = df_cxcl13[df_cxcl13['p.adj'] < 0.05]
sig_b2m = df_b2m[df_b2m['p.adj'] < 0.05]

print("=" * 60)
print("scTenifoldKnk Virtual KO Summary")
print("=" * 60)
print(f"\nCXCL13 KO in CCR8+ Treg:")
print(f"  Total genes: {len(df_cxcl13)}")
print(f"  Significant genes (p.adj < 0.05): {len(sig_cxcl13)}")
print(f"  CXCL13 self-KO Z = {df_cxcl13.loc[df_cxcl13['gene'] == 'CXCL13', 'Z'].values[0]:.2f}")
print("  Top 5 genes by |Z|:")
for _, row in df_cxcl13.nlargest(5, 'absZ').iterrows():
    print(f"    {row['gene']:12s}  Z = {row['Z']:7.3f}  p.adj = {row['p.adj']:.2e}")

print(f"\nB2M KO in MKI67+ Treg:")
print(f"  Total genes: {len(df_b2m)}")
print(f"  Significant genes (p.adj < 0.05): {len(sig_b2m)}")
print(f"  B2M self-KO Z = {df_b2m.loc[df_b2m['gene'] == 'B2M', 'Z'].values[0]:.2f}")
print("  Top 5 genes by |Z|:")
for _, row in df_b2m.nlargest(5, 'absZ').iterrows():
    print(f"    {row['gene']:12s}  Z = {row['Z']:7.3f}  p.adj = {row['p.adj']:.2e}")
print("=" * 60)

# ------------------------------------------------------------------
# Helper: volcano plot
# ------------------------------------------------------------------
def plot_volcano(ax, df, title, exclude_gene=None, top_n=10):
    """
    Volcano-style scatter on axis `ax`.
    Non-significant in grey; significant colored by Z direction.
    """
    sig = df['p.adj'] < 0.05
    z_pos = sig & (df['Z'] > 0)
    z_neg = sig & (df['Z'] < 0)

    # Non-sig grey
    ax.scatter(
        df.loc[~sig, 'Z'],
        df.loc[~sig, 'log10_padj'],
        c=COLOR_NS, s=8, alpha=0.6, edgecolors='none', rasterized=True,
        label='n.s.'
    )
    # Sig positive Z (reddish)
    ax.scatter(
        df.loc[z_pos, 'Z'],
        df.loc[z_pos, 'log10_padj'],
        c=COLOR_CXCL13 if 'CXCL13' in title else COLOR_B2M,
        s=18, alpha=0.85, edgecolors='none',
        label='p.adj < 0.05'
    )
    # Sig negative Z (blueish) – same palette but we keep it simple:
    # the prompt asks red-blue gradient for Z>0 vs Z<0.
    neg_color = '#3C5488'  # dark blue
    ax.scatter(
        df.loc[z_neg, 'Z'],
        df.loc[z_neg, 'log10_padj'],
        c=neg_color, s=18, alpha=0.85, edgecolors='none',
    )

    # Annotate top N by |Z| (excluding self-KO if requested)
    annot_df = df.copy()
    if exclude_gene:
        annot_df = annot_df[annot_df['gene'] != exclude_gene]
    top_genes = annot_df.nlargest(top_n, 'absZ')

    for _, row in top_genes.iterrows():
        ax.annotate(
            row['gene'],
            xy=(row['Z'], row['log10_padj']),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=7,
            fontfamily='Arial',
            color='black',
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.7),
        )

    # Threshold lines
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)

    ax.set_xlabel('Z-score', fontfamily='Arial')
    ax.set_ylabel('-log$_{10}$(adjusted p-value)', fontfamily='Arial')
    ax.set_title(title, fontfamily='Arial', fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # legend
    legend_elements = [
        Patch(facecolor=COLOR_NS, edgecolor='none', label='n.s.'),
        Patch(facecolor=COLOR_CXCL13 if 'CXCL13' in title else COLOR_B2M,
              edgecolor='none', label='Z > 0, p.adj < 0.05'),
        Patch(facecolor=neg_color, edgecolor='none', label='Z < 0, p.adj < 0.05'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=False, fontsize=7)

# ------------------------------------------------------------------
# Figure 5A & 5B (volcano plots)
# ------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(5.3, 3))

plot_volcano(axes[0], df_cxcl13, 'CXCL13 KO in CCR8+ Treg', exclude_gene='CXCL13', top_n=10)
plot_volcano(axes[1], df_b2m, 'B2M KO in MKI67+ Treg', exclude_gene='B2M', top_n=10)

plt.tight_layout()
fig.savefig(os.path.join(OUTDIR, 'fig5A_volcano_cxcl13.pdf'), format='pdf', bbox_inches='tight')
fig.savefig(os.path.join(OUTDIR, 'fig5A_volcano_cxcl13.png'), format='png', bbox_inches='tight')
fig.savefig(os.path.join(OUTDIR, 'fig5B_volcano_b2m.pdf'), format='pdf', bbox_inches='tight')
fig.savefig(os.path.join(OUTDIR, 'fig5B_volcano_b2m.png'), format='png', bbox_inches='tight')
plt.close(fig)
print("Saved fig5A and fig5B (volcano plots).")

# ------------------------------------------------------------------
# Figure 5C: Top 15 significant genes (combined KOs)
# ------------------------------------------------------------------
sig_cxcl13 = df_cxcl13[df_cxcl13['p.adj'] < 0.05].copy()
sig_b2m = df_b2m[df_b2m['p.adj'] < 0.05].copy()

sig_cxcl13['KO'] = 'CXCL13 KO'
sig_b2m['KO'] = 'B2M KO'

combined = pd.concat([sig_cxcl13, sig_b2m], ignore_index=True)
combined_top15 = combined.nlargest(15, 'absZ').sort_values('absZ', ascending=True)

# Horizontal bar plot
fig, ax = plt.subplots(figsize=(3.5, 4.5))
colors = [COLOR_CXCL13 if ko == 'CXCL13 KO' else COLOR_B2M for ko in combined_top15['KO']]

bars = ax.barh(
    combined_top15['gene'],
    combined_top15['absZ'],
    color=colors,
    edgecolor='white',
    height=0.65,
)

ax.set_xlabel('|Z-score|', fontfamily='Arial')
ax.set_title('Top Perturbed Genes by Virtual KO', fontfamily='Arial', fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend
legend_elements = [
    Patch(facecolor=COLOR_CXCL13, edgecolor='none', label='CXCL13 KO (CCR8+ Treg)'),
    Patch(facecolor=COLOR_B2M, edgecolor='none', label='B2M KO (MKI67+ Treg)'),
]
ax.legend(handles=legend_elements, loc='lower right', frameon=False, fontsize=8)

plt.tight_layout()
fig.savefig(os.path.join(OUTDIR, 'fig5C_top15_perturbed_genes.pdf'), format='pdf', bbox_inches='tight')
fig.savefig(os.path.join(OUTDIR, 'fig5C_top15_perturbed_genes.png'), format='png', bbox_inches='tight')
plt.close(fig)
print("Saved fig5C (top 15 perturbed genes).")

# ------------------------------------------------------------------
# Figure 5D: GO BP enrichment (manual mapping) for B2M KO
# ------------------------------------------------------------------
# Curated pathway gene sets (case-insensitive)
PATHWAY_MAP = {
    'Antigen Processing & Presentation': [
        'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB',
        'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DQB2',
        'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5',
        'HLA-E', 'HLA-F', 'HLA-G', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L',
        'HLA-P', 'HLA-S', 'HLA-V', 'HLA-W', 'HLA-Z',
        'B2M', 'TAP1', 'TAP2', 'PSMB8', 'PSMB9', 'CALR', 'CANX', 'PDIA3',
        'HSP90AA1', 'HSP90AB1', 'CD74',
    ],
    'MHC Class I Signaling': [
        'HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G',
        'B2M', 'TAP1', 'TAP2', 'NLRC5', 'CIITA',
    ],
    'MHC Class II Signaling': [
        'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5',
        'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DQB2',
        'HLA-DPA1', 'HLA-DPB1', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB',
        'CD74',
    ],
    'Immune Response': [
        'CD8A', 'CD8B', 'CD4', 'IFNG', 'IL2', 'IL2RA', 'FOXP3',
        'CTLA4', 'PDCD1', 'LAG3', 'TIGIT', 'GZMA', 'GZMB', 'PRF1',
        'TNFSF4', 'TNFRSF4', 'ICOS', 'ICOSLG', 'CD28', 'CD80', 'CD86',
        'CCR8', 'CXCL13', 'CCL22', 'CCL17', 'CCL4', 'CCL5',
        'STAT1', 'STAT3', 'STAT5A', 'STAT5B', 'IRF1', 'IRF4', 'IRF8',
        'NFKB1', 'NFKB2', 'RELA', 'RELB', 'HAVCR1',
    ],
    'Cell Cycle': [
        'MKI67', 'TOP2A', 'CDK1', 'CCNB1', 'CCNB2', 'CCNA2', 'CDK2',
        'CCNE1', 'CCNE2', 'MCM2', 'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MCM7',
        'PCNA', 'AURKA', 'AURKB', 'PLK1', 'BUB1', 'BUB1B', 'CDC20',
        'CDCA8', 'NDC80', 'KPNA2',
    ],
    'Apoptosis': [
        'BCL2', 'BAX', 'BAK1', 'CASP3', 'CASP8', 'CASP9', 'CASP7',
        'FAS', 'FASLG', 'TNF', 'TNFSF10', 'TNFRSF10A', 'TNFRSF10B',
        'BID', 'BAD', 'BAK1', 'BCL2L1', 'MCL1', 'BCL2A1', 'PMAIP1',
        'BBC3', 'BCL2L11', 'BIM',
    ],
}

# Normalise gene names to upper case for matching
sig_b2m_genes = set(g.upper() for g in sig_b2m['gene'].tolist())

pathway_counts = []
for pathway_name, genes in PATHWAY_MAP.items():
    matched = [g for g in genes if g.upper() in sig_b2m_genes]
    count = len(matched)
    if count > 0:
        pathway_counts.append({
            'pathway': pathway_name,
            'count': count,
            'genes': ', '.join(matched[:10]) + ('...' if len(matched) > 10 else ''),
        })

pathway_df = pd.DataFrame(pathway_counts).sort_values('count', ascending=True)

print(f"\nB2M KO significant genes mapped to curated pathways:")
for _, row in pathway_df.iterrows():
    print(f"  {row['pathway']:35s} : {row['count']:2d} genes")

fig, ax = plt.subplots(figsize=(3.5, 3))
bars = ax.barh(
    pathway_df['pathway'],
    pathway_df['count'],
    color=COLOR_B2M,
    edgecolor='white',
    height=0.55,
)

ax.set_xlabel('Number of significant genes', fontfamily='Arial')
ax.set_title('Pathway Enrichment: B2M KO in MKI67+ Treg', fontfamily='Arial', fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add count labels
for bar, count in zip(bars, pathway_df['count']):
    ax.text(
        bar.get_width() + 0.3,
        bar.get_y() + bar.get_height() / 2,
        str(int(count)),
        va='center',
        ha='left',
        fontsize=9,
        fontfamily='Arial',
    )

plt.tight_layout()
fig.savefig(os.path.join(OUTDIR, 'fig5D_pathway_enrichment_b2m.pdf'), format='pdf', bbox_inches='tight')
fig.savefig(os.path.join(OUTDIR, 'fig5D_pathway_enrichment_b2m.png'), format='png', bbox_inches='tight')
plt.close(fig)
print("Saved fig5D (pathway enrichment).")

# ------------------------------------------------------------------
# Combined 4-panel figure — 2×2 GridSpec with adequate spacing
fig = plt.figure(figsize=(8.0, 7.0))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1],
                      wspace=0.40, hspace=0.30,
                      left=0.10, right=0.95, bottom=0.10, top=0.95)

ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])

plot_volcano(ax_a, df_cxcl13, 'CXCL13 KO in CCR8+ Treg', exclude_gene='CXCL13', top_n=10)
plot_volcano(ax_b, df_b2m, 'B2M KO in MKI67+ Treg', exclude_gene='B2M', top_n=10)

# 5C inside combined figure
sig_cxcl13 = df_cxcl13[df_cxcl13['p.adj'] < 0.05].copy()
sig_b2m = df_b2m[df_b2m['p.adj'] < 0.05].copy()
sig_cxcl13['KO'] = 'CXCL13 KO'
sig_b2m['KO'] = 'B2M KO'
combined = pd.concat([sig_cxcl13, sig_b2m], ignore_index=True)
combined_top15 = combined.nlargest(15, 'absZ').sort_values('absZ', ascending=True)
colors = [COLOR_CXCL13 if ko == 'CXCL13 KO' else COLOR_B2M for ko in combined_top15['KO']]
ax_c.barh(combined_top15['gene'], combined_top15['absZ'], color=colors, edgecolor='white', height=0.65)
ax_c.set_xlabel('|Z-score|')
ax_c.set_title('Top Perturbed Genes by Virtual KO', fontweight='bold')
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
legend_elements = [
    Patch(facecolor=COLOR_CXCL13, edgecolor='none', label='CXCL13 KO (CCR8+ Treg)'),
    Patch(facecolor=COLOR_B2M, edgecolor='none', label='B2M KO (MKI67+ Treg)'),
]
ax_c.legend(handles=legend_elements, loc='lower right', frameon=False, fontsize=7)

# 5D inside combined figure
ax_d.barh(pathway_df['pathway'], pathway_df['count'], color=COLOR_B2M, edgecolor='white', height=0.55)
ax_d.set_xlabel('Number of significant genes')
ax_d.set_title('Pathway Enrichment: B2M KO in MKI67+ Treg', fontweight='bold')
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)
for bar, count in zip(ax_d.patches, pathway_df['count']):
    ax_d.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
              str(int(count)), va='center', ha='left', fontsize=8, fontfamily='Arial')

# Add panel labels
for ax, label in zip([ax_a, ax_b, ax_c, ax_d], ['A', 'B', 'C', 'D']):
    ax.text(-0.12, 1.05, label, transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top', ha='right', fontfamily='Arial')

fig.savefig(os.path.join(OUTDIR, 'fig5_combined.pdf'), format='pdf', bbox_inches='tight', pad_inches=0.15)
fig.savefig(os.path.join(OUTDIR, 'fig5_combined.png'), format='png', bbox_inches='tight', pad_inches=0.15)
plt.close(fig)
print("Saved fig5_combined (4-panel layout).")
print("\nAll outputs written to:", OUTDIR)
