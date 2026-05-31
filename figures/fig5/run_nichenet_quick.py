"""
Quick NicheNet analysis using pre-downloaded matrix.
Predict downstream targets of MHC-I (HLA-A/B/C) and MIF ligands.
"""
import pyreadr
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUTDIR = "D:/research/cucumber/fig5/nichenet_output"
os.makedirs(OUTDIR, exist_ok=True)

print("Loading NicheNet ligand-target matrix...")
ltm = pyreadr.read_r('D:/Research/tomato/data/nichenet_ligand_target_matrix_uncompressed.rds')[None]
print(f"  Matrix: {ltm.shape}")

# ============================================================
# MHC-I ligands: HLA-A, HLA-B, HLA-C
# ============================================================
print("\n=== MHC-I Ligands (HLA-A/B/C) ===")
mhci_ligands = ['HLA-A', 'HLA-B', 'HLA-C']
mhci_ligands = [l for l in mhci_ligands if l in ltm.columns]
print(f"  Found: {mhci_ligands}")

# Average regulatory potential across the three HLA ligands
mhci_scores = ltm[mhci_ligands].mean(axis=1).sort_values(ascending=False)
mhci_scores.name = 'MHC-I_regulatory_potential'

print(f"  Top 20 predicted targets:")
for gene, score in mhci_scores.head(20).items():
    print(f"    {gene:15s}: {score:.6f}")

# Check our key pathway genes
print(f"\n  Key pathway genes in MHC-I targets:")
key_genes = ['CD74', 'HLA-DRB5', 'HLA-DRA', 'HLA-DRB1', 'CD44', 'FOXP3', 'CTLA4',
             'HSPA1A', 'HSPA1B', 'HSP90AA1', 'DNAJB1', 'BAG3', 'NFKB1', 'JUN', 'FOS',
             'IL2RA', 'TNFRSF18', 'MKI67', 'TOP2A', 'IRF4', 'PRDM1', 'BCL6', 'RFX5', 'STAT2', 'ELK1']
for g in key_genes:
    if g in mhci_scores.index:
        rank = list(mhci_scores.index).index(g) + 1
        print(f"    {g:15s}: rank={rank}, score={mhci_scores[g]:.6f}")

mhci_scores.to_csv(f"{OUTDIR}/mhci_target_scores.csv", header=True)

# ============================================================
# MIF ligand
# ============================================================
print(f"\n=== MIF Ligand ===")
if 'MIF' in ltm.columns:
    mif_scores = ltm['MIF'].sort_values(ascending=False)
    mif_scores.name = 'MIF_regulatory_potential'
    
    print(f"  Top 20 predicted targets:")
    for gene, score in mif_scores.head(20).items():
        print(f"    {gene:15s}: {score:.6f}")
    
    print(f"\n  Key pathway genes in MIF targets:")
    for g in key_genes:
        if g in mif_scores.index:
            rank = list(mif_scores.index).index(g) + 1
            print(f"    {g:15s}: rank={rank}, score={mif_scores[g]:.6f}")
    
    mif_scores.to_csv(f"{OUTDIR}/mif_target_scores.csv", header=True)
else:
    print("  MIF NOT FOUND in matrix")

# ============================================================
# CXCL13 for comparison (our original scTenifoldKnk target)
# ============================================================
print(f"\n=== CXCL13 ===")
if 'CXCL13' in ltm.columns:
    cxcl13_scores = ltm['CXCL13'].sort_values(ascending=False)
    print(f"  Top 10:")
    for gene, score in cxcl13_scores.head(10).items():
        print(f"    {gene:15s}: {score:.6f}")

# ============================================================
# Generate comparison figure
# ============================================================
print("\n=== Generating Figure ===")

fig, axes = plt.subplots(1, 3, figsize=(16, 6))

# MHC-I targets
ax = axes[0]
top20 = mhci_scores.head(20).iloc[::-1]
colors = ['#2E5AAC' if g in ['CD74', 'HLA-DRB5', 'HLA-DRA', 'HLA-DRB1', 'CD44', 'B2M'] 
          else '#C0392B' if g in ['HSPA1A', 'HSPA1B', 'HSP90AA1', 'DNAJB1', 'BAG3', 'HSPH1']
          else '#7F8C8D' for g in top20.index]
ax.barh(range(len(top20)), top20.values, color=colors, height=0.7, edgecolor='white')
ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20.index, fontsize=9)
ax.invert_yaxis()
ax.set_xlabel('Regulatory Potential', fontsize=10)
ax.set_title('MHC-I (HLA-A/B/C): Predicted Targets', fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# MIF targets
ax = axes[1]
top20 = mif_scores.head(20).iloc[::-1]
colors = ['#8E44AD' if g in ['CD74', 'CD44'] 
          else '#27AE60' if g in ['NFKB1', 'JUN', 'FOS', 'RELA', 'NFKB2']
          else '#7F8C8D' for g in top20.index]
ax.barh(range(len(top20)), top20.values, color=colors, height=0.7, edgecolor='white')
ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20.index, fontsize=9)
ax.invert_yaxis()
ax.set_xlabel('Regulatory Potential', fontsize=10)
ax.set_title('MIF: Predicted Targets', fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# MHC-I + MIF: comparison of shared targets
ax = axes[2]
# Combine into a dataframe
combined = pd.DataFrame({
    'MHC-I': mhci_scores,
    'MIF': mif_scores,
})
combined = combined.fillna(0)

# Select genes that are in top 100 for either
top100_mhci = set(mhci_scores.head(100).index)
top100_mif = set(mif_scores.head(100).index)
union = top100_mhci | top100_mif

# Key pathway genes to highlight
highlight_df = combined.loc[combined.index.intersection(key_genes)]
highlight_df = highlight_df.sort_values('MHC-I', ascending=False)

if len(highlight_df) > 0:
    for i, (gene, row) in enumerate(highlight_df.iterrows()):
        ax.scatter(row['MHC-I'], row['MIF'], s=60, 
                   c='#E74C3C' if gene in ['CD74', 'HLA-DRB5', 'HLA-DRA', 'HLA-DRB1'] 
                   else '#C0392B' if gene in ['HSPA1A', 'HSPA1B', 'HSP90AA1', 'DNAJB1', 'BAG3']
                   else '#95A5A6',
                   edgecolors='black', linewidth=0.3, zorder=3)
        ax.annotate(gene, (row['MHC-I'], row['MIF']), fontsize=7,
                    xytext=(3, 3), textcoords='offset points')

# Background
bg_genes = [g for g in combined.index if g not in highlight_df.index]
ax.scatter(combined.loc[bg_genes, 'MHC-I'], combined.loc[bg_genes, 'MIF'],
           s=5, c='#ECF0F1', edgecolors='#BDC3C7', linewidth=0.2, alpha=0.5)

ax.set_xlabel('MHC-I Regulatory Potential', fontsize=10)
ax.set_ylabel('MIF Regulatory Potential', fontsize=10)
ax.set_title('MHC-I vs MIF: Target Comparison', fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#E74C3C', label='MHC-II / MIF pathway'),
    Patch(facecolor='#C0392B', label='HSF1 stress program'),
    Patch(facecolor='#95A5A6', label='Other key genes'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=7, frameon=False)

plt.tight_layout()
plt.savefig(f"{OUTDIR}/nichenet_combined.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{OUTDIR}/nichenet_combined.pdf", dpi=300, bbox_inches='tight')
plt.close()

print(f"\nSaved: {OUTDIR}/nichenet_combined.png/pdf")
print("\nDONE!")
