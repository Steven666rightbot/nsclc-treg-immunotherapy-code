"""
NicheNet-style ligand-target analysis
Validates Fibroblast ECM ligands -> Treg mechanosensing targets
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Try to read RDS files
try:
    import pyreadr
    HAS_PYREADR = True
except ImportError:
    HAS_PYREADR = False
    print("pyreadr not installed. Attempting to install...")
    import subprocess
    subprocess.check_call(["python", "-m", "pip", "install", "pyreadr", "-q"])
    import pyreadr
    HAS_PYREADR = True

# Paths
LT_MATRIX = 'D:/Research/tomato/data/nichenet_ligand_target_matrix.rds'
LR_NETWORK = 'D:/Research/tomato/data/nichenet_lr_network.rds'
OUT_DIR = 'D:/Research/tomato/figures'
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("NicheNet Ligand-Target Validation")
print("=" * 60)

# 1. Load ligand-target matrix
print("\n[1] Loading ligand-target matrix...")
lt_result = pyreadr.read_r(LT_MATRIX)
lt_matrix = list(lt_result.values())[0]
print(f"    Matrix shape: {lt_matrix.shape}")
print(f"    Index (ligands): {lt_matrix.index[:5].tolist()}...")
print(f"    Columns (targets): {lt_matrix.columns[:5].tolist()}...")

# 2. Load ligand-receptor network
print("\n[2] Loading ligand-receptor network...")
lr_result = pyreadr.read_r(LR_NETWORK)
lr_network = list(lr_result.values())[0]
print(f"    Network shape: {lr_network.shape}")
print(f"    Columns: {lr_network.columns.tolist()}")

# 3. Define our ligands of interest (from CellChat results)
# Fibroblast ECM ligands active in Weak responder
cellchat_ligands = [
    'COL1A1', 'COL1A2', 'COL6A1', 'COL6A2', 'COL6A3', 'COL4A2',
    'FN1', 'LAMA4', 'LAMB1', 'LAMC1', 'MDK', 'LGALS9', 'THBS1', 'THBS2'
]

# Receptors expressed in Treg (from our analysis)
treg_receptors = ['CD44', 'ITGA4', 'ITGB1', 'ITGAV']

# Downstream targets (contractile + mechanotransduction genes)
target_genes = [
    'ACTA2', 'MYL9', 'VCL', 'TLN1', 'TRPV4', 'MOB1B',
    'ITGA1', 'DDR1', 'DDR2', 'SRC', 'RHOA', 'WWTR1',
    'YAP1', 'TEAD1', 'TEAD3', 'TEAD4', 'CTGF'
]

# 4. Filter ligands present in matrix
available_ligands = [l for l in cellchat_ligands if l in lt_matrix.index]
missing_ligands = [l for l in cellchat_ligands if l not in lt_matrix.index]
print(f"\n[3] Ligands in NicheNet matrix: {len(available_ligands)}/{len(cellchat_ligands)}")
print(f"    Available: {available_ligands}")
if missing_ligands:
    print(f"    Missing: {missing_ligands}")

# 5. Filter targets present in matrix
available_targets = [t for t in target_genes if t in lt_matrix.columns]
missing_targets = [t for t in target_genes if t not in lt_matrix.columns]
print(f"\n[4] Targets in NicheNet matrix: {len(available_targets)}/{len(target_genes)}")
print(f"    Available: {available_targets}")
if missing_targets:
    print(f"    Missing: {missing_targets}")

# 6. Extract ligand-target scores
if available_ligands and available_targets:
    sub_matrix = lt_matrix.loc[available_ligands, available_targets]
    print(f"\n[5] Sub-matrix shape: {sub_matrix.shape}")
    
    # For each ligand, find top targets
    print("\n[6] Top target predictions per ligand:")
    for ligand in available_ligands:
        scores = sub_matrix.loc[ligand].sort_values(ascending=False)
        top5 = scores.head(5)
        print(f"\n    {ligand}:")
        for target, score in top5.items():
            marker = " ***" if target in ['ACTA2', 'MYL9', 'VCL', 'TLN1', 'ITGA1', 'DDR1', 'DDR2'] else ""
            print(f"      {target}: {score:.4f}{marker}")
    
    # For each target, find top ligands
    print("\n[7] Top ligand predictions per target:")
    for target in available_targets:
        scores = sub_matrix[target].sort_values(ascending=False)
        top5 = scores.head(5)
        print(f"\n    {target}:")
        for ligand, score in top5.items():
            marker = " ***" if ligand in ['COL1A1', 'COL1A2', 'FN1'] else ""
            print(f"      {ligand}: {score:.4f}{marker}")
    
    # 7. Aggregate scores: sum across all ligands for each target
    aggregated = sub_matrix.sum(axis=0).sort_values(ascending=False)
    print("\n[8] Aggregated ligand-target scores (sum across all ECM ligands):")
    for target, score in aggregated.items():
        print(f"    {target}: {score:.4f}")
    
    # 8. Visualization
    print("\n[9] Generating heatmap...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))
    
    # Heatmap of ligand-target scores
    ax1 = axes[0]
    sns.heatmap(sub_matrix, cmap='YlOrRd', annot=True, fmt='.3f', 
                linewidths=0.5, ax=ax1, cbar_kws={'label': 'Regulatory Potential'})
    ax1.set_title('NicheNet: Ligand-Target Regulatory Potential\n(Fibroblast ECM -> Treg Mechanosensing)', 
                  fontsize=11, fontweight='bold')
    ax1.set_xlabel('Target Genes', fontsize=10)
    ax1.set_ylabel('Ligands (from CellChat)', fontsize=10)
    
    # Barplot of aggregated scores
    ax2 = axes[1]
    colors = ['#E74C3C' if t in ['ACTA2', 'MYL9', 'VCL', 'TLN1'] 
              else '#3498DB' if t in ['ITGA1', 'DDR1', 'DDR2'] 
              else '#95A5A6' for t in aggregated.index]
    bars = ax2.barh(range(len(aggregated)), aggregated.values, color=colors, alpha=0.8, edgecolor='black')
    ax2.set_yticks(range(len(aggregated)))
    ax2.set_yticklabels(aggregated.index, fontsize=9)
    ax2.set_xlabel('Aggregated Regulatory Potential', fontsize=10)
    ax2.set_title('Aggregated ECM Ligand Potential\nper Target Gene', fontsize=11, fontweight='bold')
    ax2.invert_yaxis()
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E74C3C', label='Cytoskeleton (ACTA2/MYL9/VCL/TLN1)', alpha=0.8),
        Patch(facecolor='#3498DB', label='ECM Receptor (ITGA1/DDR1/DDR2)', alpha=0.8),
        Patch(facecolor='#95A5A6', label='Other', alpha=0.8)
    ]
    ax2.legend(handles=legend_elements, loc='lower right', fontsize=8)
    
    plt.tight_layout()
    out_path = os.path.join(OUT_DIR, 'nichenet_ligand_target_heatmap.png')
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"    Saved: {out_path}")
    plt.close()
    
    # 9. Save results table
    sub_matrix.to_csv(os.path.join(OUT_DIR, 'nichenet_ligand_target_matrix_subset.csv'))
    aggregated.to_csv(os.path.join(OUT_DIR, 'nichenet_aggregated_scores.csv'))
    print(f"    Tables saved to {OUT_DIR}")

else:
    print("ERROR: No overlapping ligands or targets found!")

print("\n" + "=" * 60)
print("NicheNet analysis complete!")
print("=" * 60)
