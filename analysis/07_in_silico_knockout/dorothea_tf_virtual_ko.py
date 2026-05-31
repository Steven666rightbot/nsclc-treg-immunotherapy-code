#!/usr/bin/env python3
"""
Dorothea TF Virtual Knockout Analysis
Simulate TF knockout by reducing target gene expression,
then measure impact on Treg subtype signatures.
"""

import pandas as pd
import numpy as np
import loompy as lp
import matplotlib.pyplot as plt
import seaborn as sns
import decoupler as dc
from pathlib import Path
import os

os.makedirs('figures', exist_ok=True)
os.makedirs('results/tf_ko', exist_ok=True)

# ======================== 配置 ========================

LOOM_PATH = 'results/scenic_input/treg_no_foxp3_8k.loom'
OUTPUT_DIR = Path('results/tf_ko')

# Target TFs: mechanosensitive + subtype-differentiated
TARGET_TFS = [
    'MEF2C',   # integrin-FAK-HDAC, CCR8+ enriched
    'TEAD1',   # YAP/TAZ partner, MKI67+ enriched
    'NFKB1',   # NF-kB, mechanical stress
    'RELA',    # NF-kB p65
    'STAT2',   # MKI67+ enriched
    'IRF8',    # MKI67+ enriched
    'JUN',     # AP-1, mechanosensitive
    'SP1',     # MKI67+ vs CCR8+ differential
    'MAX',     # Top CCR8+ TF
    'ATF3',    # CCR8+ related
]

# Signature gene sets
SIGNATURES = {
    'CCR8_Treg': [
        'DIRAS3', 'ACTG2', 'TMPRSS6', 'LINC00519', 'CD177', 'XIRP1',
        'CACNB2', 'GNG8', 'CPE', 'LINC02099', 'ECEL1', 'LAYN', 'EEF1A2',
        'TNFRSF4', 'AC133644.2', 'C2CD4A', 'EPHX2', 'MAGEH1', 'ADTRP',
        'LRRC32', 'TNFRSF18', 'BATF', 'CTLA4'
    ],
    'MKI67_Treg': [
        'LINC00589', 'FANK1', 'PTN', 'FOXP3', 'MEOX1', 'AC017002.3',
        'RTKN2', 'TRAV36DV7', 'HLF', 'F5', 'IL12RB2', 'TLE2', 'LINC01229',
        'TSHR', 'IL2RA', 'CPA5', 'ZBED2', 'LAIR2', 'GPR19', 'IKZF4',
        'HACD1', 'CEP55', 'DDX11-AS1', 'NEIL3', 'FAM184A', 'BTNL8',
        'LINC01281', 'XXYLT1-AS2', 'MKI67', 'TOP2A', 'STMN1', 'TYMS', 'PCNA'
    ],
    'Contractile': [
        'ACTA2', 'MYL9', 'VCL', 'TLN1', 'TRPV4', 'MOB1B',
        'ITGA1', 'DDR1', 'DDR2', 'SRC', 'ITGB1', 'CD44', 'PTK2', 'ZYX', 'FLNA',
        'FOXP3', 'IL2RA', 'CTLA4'
    ],
    'Mechanoscore': [
        'YAP1', 'WWTR1', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4',
        'ITGB1', 'ITGA1', 'ITGA5', 'PTK2', 'SRC', 'VCL', 'TLN1'
    ]
}

# Clean signatures: remove empty/None
for k, v in SIGNATURES.items():
    SIGNATURES[k] = [g for g in v if g and isinstance(g, str)]

# ======================== 加载数据 ========================

print("=" * 60)
print("Loading expression matrix from loom...")
print("=" * 60)

with lp.connect(LOOM_PATH, mode='r') as ds:
    genes = ds.ra.Gene.tolist()
    cells = ds.ca.CellID.tolist()
    cell_types = ds.ca.sub_cell_type.tolist()
    responses = ds.ca.response.tolist()
    # decoupler expects cells x genes
    exp_mtx = pd.DataFrame(ds[:, :].T, index=cells, columns=genes)

print(f"Matrix: {exp_mtx.shape}")
print(f"Responders: {sum(1 for r in responses if r == 'Responder')}")
print(f"Non-responders: {sum(1 for r in responses if r == 'Non-responder')}")

# ======================== 加载 Dorothea regulons ========================

print("\n" + "=" * 60)
print("Loading Dorothea regulons (A+B)...")
print("=" * 60)

net = dc.get_dorothea(organism='human', levels=['A', 'B'])
print(f"Regulons: {net['source'].nunique()} TFs, {len(net)} interactions")

# Filter to expressed genes
net = net[net['target'].isin(exp_mtx.columns)]
print(f"After filtering: {net['source'].nunique()} TFs, {len(net)} interactions")

# ======================== 计算签名得分函数 ========================

def calc_signature_scores(mtx, sig_dict):
    """Calculate mean log1p expression of signature genes per cell."""
    scores = pd.DataFrame(index=mtx.index)
    for name, gene_list in sig_dict.items():
        available = [g for g in gene_list if g in mtx.columns]
        if available:
            scores[name] = np.log1p(mtx[available]).mean(axis=1)
        else:
            scores[name] = 0
            print(f"  WARNING: No genes available for {name}")
    return scores

# Base scores
print("\nCalculating baseline signature scores...")
base_scores = calc_signature_scores(exp_mtx, SIGNATURES)
base_scores['response'] = responses
print(base_scores.head())

# ======================== TF 虚拟敲除 ========================

print("\n" + "=" * 60)
print("Running TF virtual knockout (50% target reduction)")
print("=" * 60)

ko_results = []

for tf in TARGET_TFS:
    print(f"\n--- KO: {tf} ---")
    
    # Get target genes for this TF
    tf_net = net[net['source'] == tf]
    if len(tf_net) == 0:
        print(f"  No regulon found for {tf}, skipping")
        continue
    
    targets = tf_net['target'].unique().tolist()
    targets_in_data = [t for t in targets if t in exp_mtx.columns]
    print(f"  Target genes: {len(targets)} total, {len(targets_in_data)} in data")
    
    if len(targets_in_data) == 0:
        print(f"  No targets in expression matrix, skipping")
        continue
    
    # Create KO matrix: reduce target expression by 50%
    ko_mtx = exp_mtx.copy()
    ko_mtx[targets_in_data] = ko_mtx[targets_in_data] * 0.5
    
    # Calculate KO signature scores
    ko_scores = calc_signature_scores(ko_mtx, SIGNATURES)
    
    # Compare: mean change per signature
    for sig_name in SIGNATURES.keys():
        base_mean = base_scores[sig_name].mean()
        ko_mean = ko_scores[sig_name].mean()
        change = ko_mean - base_mean
        pct_change = (change / (base_mean + 1e-10)) * 100
        
        # Also compare by response group
        for resp in ['Responder', 'Non-responder']:
            mask = [r == resp for r in responses]
            base_group = base_scores.loc[mask, sig_name].mean()
            ko_group = ko_scores.loc[mask, sig_name].mean()
            group_change = ko_group - base_group
            
            ko_results.append({
                'TF': tf,
                'Signature': sig_name,
                'Response': resp,
                'Base_Mean': base_group,
                'KO_Mean': ko_group,
                'Abs_Change': group_change,
                'Pct_Change': (group_change / (base_group + 1e-10)) * 100,
                'N_Targets': len(targets_in_data)
            })
    
    print(f"  CCR8 change: {ko_results[-4]['Abs_Change']:.4f}")
    print(f"  MKI67 change: {ko_results[-3]['Abs_Change']:.4f}")

# Save results
ko_df = pd.DataFrame(ko_results)
ko_df.to_csv(OUTPUT_DIR / 'tf_ko_signature_changes.csv', index=False)
print(f"\nSaved: {OUTPUT_DIR / 'tf_ko_signature_changes.csv'}")

# ======================== 可视化 ========================

print("\n" + "=" * 60)
print("Generating visualizations...")
print("=" * 60)

# Pivot for heatmap
pivot_abs = ko_df.pivot_table(
    index='TF', columns=['Signature', 'Response'],
    values='Abs_Change'
)

pivot_pct = ko_df.pivot_table(
    index='TF', columns=['Signature', 'Response'],
    values='Pct_Change'
)

# Plot 1: Absolute change heatmap
fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(pivot_abs, cmap='RdBu_r', center=0,
            annot=True, fmt='.3f', linewidths=0.5,
            cbar_kws={'label': 'Signature Score Change (KO - Base)'},
            ax=ax)
ax.set_title('TF Virtual KO: Absolute Signature Score Changes', fontsize=14, fontweight='bold')
ax.set_xlabel('Signature × Response Group', fontsize=12)
ax.set_ylabel('Transcription Factor', fontsize=12)
plt.tight_layout()
fig.savefig('figures/tf_ko_absolute_changes.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig('figures/tf_ko_absolute_changes.pdf', dpi=300, bbox_inches='tight', facecolor='white')
print("Saved: figures/tf_ko_absolute_changes.png")

# Plot 2: Percentage change heatmap
fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(pivot_pct, cmap='RdBu_r', center=0,
            annot=True, fmt='.1f', linewidths=0.5,
            cbar_kws={'label': 'Signature Score Change (%)'},
            ax=ax)
ax.set_title('TF Virtual KO: Percentage Signature Score Changes', fontsize=14, fontweight='bold')
ax.set_xlabel('Signature × Response Group', fontsize=12)
ax.set_ylabel('Transcription Factor', fontsize=12)
plt.tight_layout()
fig.savefig('figures/tf_ko_percentage_changes.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig('figures/tf_ko_percentage_changes.pdf', dpi=300, bbox_inches='tight', facecolor='white')
print("Saved: figures/tf_ko_percentage_changes.png")

# Plot 3: Bar plot of key TFs vs key signatures
key_sigs = ['CCR8_Treg', 'MKI67_Treg', 'Contractile']
plot_df = ko_df[ko_df['Signature'].isin(key_sigs)]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for idx, sig in enumerate(key_sigs):
    ax = axes[idx]
    sig_df = plot_df[plot_df['Signature'] == sig]
    pivot = sig_df.pivot(index='TF', columns='Response', values='Abs_Change')
    pivot.plot(kind='barh', ax=ax, color=['#3498DB', '#E74C3C'], width=0.7)
    ax.set_title(f'{sig}', fontsize=12, fontweight='bold')
    ax.set_xlabel('Score Change (KO - Base)', fontsize=10)
    ax.axvline(0, color='black', linewidth=0.8)
    ax.legend(title='Response', loc='lower right')

plt.suptitle('TF Virtual KO Impact on Treg Subtype Signatures', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig('figures/tf_ko_subtype_comparison.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig('figures/tf_ko_subtype_comparison.pdf', dpi=300, bbox_inches='tight', facecolor='white')
print("Saved: figures/tf_ko_subtype_comparison.png")

# ======================== 关键发现摘要 ========================

print("\n" + "=" * 60)
print("KEY FINDINGS")
print("=" * 60)

for sig in key_sigs:
    print(f"\n{sig}:")
    sig_df = ko_df[ko_df['Signature'] == sig]
    top_down = sig_df.nsmallest(3, 'Abs_Change')
    top_up = sig_df.nlargest(3, 'Abs_Change')
    
    print("  Most decreased by KO:")
    for _, row in top_down.iterrows():
        print(f"    {row['TF']} ({row['Response']}): {row['Abs_Change']:.4f} ({row['Pct_Change']:.1f}%)")
    
    print("  Most increased by KO:")
    for _, row in top_up.iterrows():
        print(f"    {row['TF']} ({row['Response']}): {row['Abs_Change']:.4f} ({row['Pct_Change']:.1f}%)")

print("\n" + "=" * 60)
print("Analysis complete!")
print("=" * 60)
