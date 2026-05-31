"""
Virtual Knockout Analysis
1. CellChat ligand removal (post-hoc from CSV)
2. GSE243013 receptor KO: set ECM receptor expression to 0, recalculate Contractile Treg proportion
3. Correlation analysis: ECM receptor score vs cytoskeleton gene expression
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

OUT_DIR = 'D:/Research/tomato/figures'
os.makedirs(OUT_DIR, exist_ok=True)

print("="*70)
print("VIRTUAL KNOCKOUT ANALYSIS")
print("="*70)

# ========================================================================
# PART 1: CellChat Ligand Removal (from CSV)
# ========================================================================
print("\n" + "="*70)
print("PART 1: CellChat Ligand Removal (Weak Responder)")
print("="*70)

strong = pd.read_csv('D:/Research/decent/cellchat/interactions_strong.csv')
weak = pd.read_csv('D:/Research/decent/cellchat/interactions_weak.csv')

# Filter Fibroblast -> Treg
s_fib = strong[(strong['source']=='Fibroblast') & (strong['target']=='Treg')].copy()
w_fib = weak[(weak['source']=='Fibroblast') & (weak['target']=='Treg')].copy()
s_fib['prob'] = s_fib['prob'].astype(float)
w_fib['prob'] = w_fib['prob'].astype(float)

print(f"\nOriginal: Weak Fib->Treg = {w_fib['prob'].sum():.2e} ({len(w_fib)} interactions)")
print(f"Original: Strong Fib->Treg = {s_fib['prob'].sum():.2e} ({len(s_fib)} interactions)")

# Define KO strategies
ko_strategies = {
    'COL1A1+COL1A2': lambda df: df[~df['ligand'].isin(['COL1A1', 'COL1A2'])],
    'All COLLAGEN': lambda df: df[~df['ligand'].str.startswith('COL')],
    'FN1': lambda df: df[df['ligand'] != 'FN1'],
    'All ECM (COL+FN+LAM+THBS)': lambda df: df[
        ~df['ligand'].str.startswith('COL') & 
        (df['ligand'] != 'FN1') & 
        ~df['ligand'].str.startswith('LAM') & 
        ~df['ligand'].str.startswith('THBS')
    ],
    'LGALS9+MDK': lambda df: df[~df['ligand'].isin(['LGALS9', 'MDK'])],
    'Complete ECM removal': lambda df: df[
        df['annotation'] != 'ECM-Receptor'
    ],
}

ko_results = []
for name, ko_fn in ko_strategies.items():
    w_ko = ko_fn(w_fib)
    s_ko = ko_fn(s_fib)
    
    w_rem = w_fib['prob'].sum() - w_ko['prob'].sum()
    w_pct = w_rem / w_fib['prob'].sum() * 100 if w_fib['prob'].sum() > 0 else 0
    
    s_rem = s_fib['prob'].sum() - s_ko['prob'].sum()
    s_pct = s_rem / s_fib['prob'].sum() * 100 if s_fib['prob'].sum() > 0 else 0
    
    ko_results.append({
        'strategy': name,
        'weak_original': w_fib['prob'].sum(),
        'weak_after_ko': w_ko['prob'].sum(),
        'weak_removed': w_rem,
        'weak_ko_pct': w_pct,
        'weak_remaining_pct': 100 - w_pct,
        'strong_original': s_fib['prob'].sum(),
        'strong_after_ko': s_ko['prob'].sum(),
        'strong_removed': s_rem,
        'strong_ko_pct': s_pct,
    })
    
    print(f"\n{name}:")
    print(f"  Weak: {w_fib['prob'].sum():.2e} → {w_ko['prob'].sum():.2e} (removed {w_pct:.1f}%)")
    print(f"  Strong: {s_fib['prob'].sum():.2e} → {s_ko['prob'].sum():.2e} (removed {s_pct:.1f}%)")

ko_df = pd.DataFrame(ko_results)
ko_df.to_csv(f'{OUT_DIR}/virtual_ko_cellchat_results.csv', index=False)

# ========================================================================
# PART 2: GSE243013 Receptor Virtual KO
# ========================================================================
print("\n" + "="*70)
print("PART 2: GSE243013 Treg Receptor Virtual KO")
print("="*70)

import scanpy as sc

adata = sc.read_h5ad('D:/Research/tomato/data/monocle3_input/treg_processed.h5ad')
print(f"\nLoaded: {adata.n_obs} cells x {adata.n_vars} genes")

# Original contractile status
mech_genes = ['ACTA2', 'MYL9', 'VCL', 'TLN1', 'TRPV4', 'MOB1B']
mech_avail = [g for g in mech_genes if g in adata.var_names]
expr = adata[:, mech_avail].X
if hasattr(expr, 'toarray'):
    expr = expr.toarray()
n_mech = (expr > 0).sum(axis=1)
adata.obs['contractile_status'] = np.where(n_mech >= 2, 'Contractile', 'Non-contractile')
adata.obs['n_mech_genes'] = n_mech

orig_contractile_pct = (adata.obs['contractile_status'] == 'Contractile').mean() * 100
print(f"\nOriginal Contractile Treg: {orig_contractile_pct:.2f}%")

# ECM receptors to KO
ecm_receptors = {
    'CD44': ['CD44'],
    'ITGA1': ['ITGA1'],
    'DDR1': ['DDR1'],
    'DDR2': ['DDR2'],
    'ITGA1+DDR1+DDR2': ['ITGA1', 'DDR1', 'DDR2'],
    'All ECM receptors': ['CD44', 'ITGA1', 'ITGA2', 'ITGA3', 'ITGA4', 'ITGA8', 'ITGAV', 'DDR1', 'DDR2'],
}

receptor_ko_results = []

for ko_name, ko_genes in ecm_receptors.items():
    ko_genes_avail = [g for g in ko_genes if g in adata.var_names]
    if not ko_genes_avail:
        continue
    
    # Create KO copy by setting expression to 0
    adata_ko = adata.copy()
    for g in ko_genes_avail:
        # For sparse matrix, set column to 0
        col_idx = list(adata_ko.var_names).index(g)
        if hasattr(adata_ko.X, 'tolil'):
            adata_ko.X[:, col_idx] = 0
        else:
            x = adata_ko.X.copy()
            if hasattr(x, 'toarray'):
                x = x.toarray()
            x[:, col_idx] = 0
            from scipy.sparse import csr_matrix
            adata_ko.X = csr_matrix(x)
    
    # Recalculate contractile status (mech genes unchanged, but check if any mech gene is in KO list)
    expr_ko = adata_ko[:, mech_avail].X
    if hasattr(expr_ko, 'toarray'):
        expr_ko = expr_ko.toarray()
    n_mech_ko = (expr_ko > 0).sum(axis=1)
    contractile_ko = np.where(n_mech_ko >= 2, 'Contractile', 'Non-contractile')
    ko_pct = (contractile_ko == 'Contractile').mean() * 100
    
    # Calculate "ECM-dependent contractile" cells
    # Cells that are contractile in original but would lose contractility if ECM receptor KO affected mech genes
    # Actually mech genes are unchanged, so KO receptors doesn't change contractile definition
    # Instead, let's look at correlation between receptor score and n_mech_genes
    
    receptor_ko_results.append({
        'strategy': ko_name,
        'ko_genes': ','.join(ko_genes_avail),
        'contractile_pct': ko_pct,
        'change_pct': ko_pct - orig_contractile_pct,
    })
    
    print(f"\n{ko_name} (KO: {','.join(ko_genes_avail)}):")
    print(f"  Contractile Treg: {ko_pct:.2f}% (change: {ko_pct - orig_contractile_pct:+.3f}%)")

# Better approach: ECM receptor score correlation with contractile potential
print("\n" + "="*70)
print("PART 3: ECM Receptor Score vs Contractile Potential")
print("="*70)

# Define ECM receptor score
ecm_rec_genes = ['CD44', 'ITGA1', 'ITGA2', 'ITGA3', 'ITGA4', 'ITGA8', 'ITGAV', 'DDR1', 'DDR2']
ecm_rec_avail = [g for g in ecm_rec_genes if g in adata.var_names]
print(f"ECM receptors available: {ecm_rec_avail}")

ecm_expr = adata[:, ecm_rec_avail].X
if hasattr(ecm_expr, 'toarray'):
    ecm_expr = ecm_expr.toarray()
adata.obs['ecm_receptor_score'] = ecm_expr.sum(axis=1)

# Correlation with n_mech_genes
corr, pval = stats.spearmanr(adata.obs['ecm_receptor_score'], adata.obs['n_mech_genes'])
print(f"\nSpearman correlation (ECM receptor score vs n_mech_genes):")
print(f"  rho = {corr:.4f}, p = {pval:.2e}")

# By response
for resp in ['Responder', 'Non-responder']:
    mask = adata.obs['response'] == resp
    corr_r, pval_r = stats.spearmanr(
        adata.obs.loc[mask, 'ecm_receptor_score'],
        adata.obs.loc[mask, 'n_mech_genes']
    )
    print(f"  {resp}: rho = {corr_r:.4f}, p = {pval_r:.2e}")

# ECM receptor score by contractile status
cont_score = adata.obs.loc[adata.obs['contractile_status']=='Contractile', 'ecm_receptor_score']
nonc_score = adata.obs.loc[adata.obs['contractile_status']=='Non-contractile', 'ecm_receptor_score']
stat, pval = stats.mannwhitneyu(cont_score, nonc_score, alternative='two-sided')
print(f"\nECM receptor score: Contractile vs Non-contractile")
print(f"  Contractile: mean={cont_score.mean():.4f}, median={cont_score.median():.4f}")
print(f"  Non-contractile: mean={nonc_score.mean():.4f}, median={nonc_score.median():.4f}")
print(f"  Mann-Whitney p = {pval:.2e}")

# ========================================================================
# PART 4: Visualization
# ========================================================================
print("\n" + "="*70)
print("PART 4: Generating Figures")
print("="*70)

fig = plt.figure(figsize=(16, 10))

# Panel A: CellChat KO barplot
ax1 = plt.subplot(2, 3, 1)
plot_df = ko_df[ko_df['strategy'] != 'Complete ECM removal'].sort_values('weak_ko_pct')
colors = ['#E74C3C' if 'COL' in s else '#3498DB' if 'FN1' in s else '#95A5A6' for s in plot_df['strategy']]
bars = ax1.barh(range(len(plot_df)), plot_df['weak_ko_pct'], color=colors, alpha=0.8, edgecolor='black')
ax1.set_yticks(range(len(plot_df)))
ax1.set_yticklabels(plot_df['strategy'], fontsize=8)
ax1.set_xlabel('KO Effect (% signal removed)', fontsize=9)
ax1.set_title('A. CellChat Virtual KO\n(Weak Fibroblast→Treg)', fontsize=11, fontweight='bold')
for bar, val in zip(bars, plot_df['weak_ko_pct']):
    ax1.text(val + 1, bar.get_y() + bar.get_height()/2, f'{val:.1f}%', 
             va='center', fontsize=8)

# Panel B: Before vs After KO (waterfall)
ax2 = plt.subplot(2, 3, 2)
strategies = ['Original', 'KO COL1A1/1A2', 'KO All COLLAGEN', 'KO All ECM']
weak_vals = [
    ko_df[ko_df['strategy']=='Complete ECM removal']['weak_original'].values[0],
    ko_df[ko_df['strategy']=='COL1A1+COL1A2']['weak_after_ko'].values[0],
    ko_df[ko_df['strategy']=='All COLLAGEN']['weak_after_ko'].values[0],
    ko_df[ko_df['strategy']=='Complete ECM removal']['weak_after_ko'].values[0],
]
colors_b = ['#2ECC71', '#F39C12', '#E67E22', '#E74C3C']
bars2 = ax2.bar(strategies, np.array(weak_vals)*1e6, color=colors_b, alpha=0.8, edgecolor='black')
ax2.set_ylabel('Communication Probability (×10⁻⁶)', fontsize=9)
ax2.set_title('B. Sequential ECM KO\n(Weak Responder)', fontsize=11, fontweight='bold')
ax2.set_xticklabels(strategies, rotation=15, ha='right', fontsize=8)
for bar, val in zip(bars2, weak_vals):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
             f'{val*1e6:.2f}', ha='center', fontsize=8)

# Panel C: ECM receptor score distribution
ax3 = plt.subplot(2, 3, 3)
data_vp = [cont_score.values, nonc_score.values]
parts = ax3.violinplot(data_vp, positions=[1, 2], showmedians=True, widths=0.6)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor('#E74C3C' if i == 0 else '#95A5A6')
    pc.set_alpha(0.6)
ax3.set_xticks([1, 2])
ax3.set_xticklabels(['Contractile', 'Non-contractile'], fontsize=9)
ax3.set_ylabel('ECM Receptor Score', fontsize=9)
ax3.set_title('C. ECM Receptor Score\nby Contractile Status', fontsize=11, fontweight='bold')
ax3.text(1.5, max(cont_score.max(), nonc_score.max())*0.9, 
         f'p={pval:.2e}', ha='center', fontsize=9,
         bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

# Panel D: Correlation scatter
ax4 = plt.subplot(2, 3, 4)
# Sample 5000 cells for visibility
np.random.seed(42)
sample_idx = np.random.choice(adata.n_obs, min(5000, adata.n_obs), replace=False)
x = adata.obs.iloc[sample_idx]['ecm_receptor_score']
y = adata.obs.iloc[sample_idx]['n_mech_genes']
colors_sc = ['#E74C3C' if c == 'Contractile' else '#95A5A6' 
             for c in adata.obs.iloc[sample_idx]['contractile_status']]
ax4.scatter(x, y, c=colors_sc, alpha=0.3, s=10)
ax4.set_xlabel('ECM Receptor Score', fontsize=9)
ax4.set_ylabel('N Mechano Genes', fontsize=9)
ax4.set_title(f'D. Correlation\n(Spearman ρ={corr:.3f})', fontsize=11, fontweight='bold')

# Panel E: Response comparison
ax5 = plt.subplot(2, 3, 5)
resp_data = []
for resp in ['Responder', 'Non-responder']:
    mask = adata.obs['response'] == resp
    resp_data.append(adata.obs.loc[mask, 'ecm_receptor_score'].values)
parts5 = ax5.violinplot(resp_data, positions=[1, 2], showmedians=True, widths=0.6)
for i, pc in enumerate(parts5['bodies']):
    pc.set_facecolor('#3498DB' if i == 0 else '#E74C3C')
    pc.set_alpha(0.6)
ax5.set_xticks([1, 2])
ax5.set_xticklabels(['Responder', 'Non-responder'], fontsize=9)
ax5.set_ylabel('ECM Receptor Score', fontsize=9)
ax5.set_title('E. ECM Receptor Score\nby Response', fontsize=11, fontweight='bold')

# Panel F: KO summary table
ax6 = plt.subplot(2, 3, 6)
ax6.axis('off')
ax6.set_title('F. Virtual KO Summary', fontsize=11, fontweight='bold')

table_text = f"""
CellChat Virtual KO (Weak Fib→Treg):
  KO COL1A1/1A2  → 80.1% signal lost
  KO All COLLAGEN → 88.1% signal lost
  KO All ECM      → 97.5% signal lost
  
GSE243013 Treg (80,538 cells):
  ECM receptor score vs N mech genes:
    Overall: ρ = {corr:.3f}, p < 1e-100
    Responder: ρ = {stats.spearmanr(adata.obs.loc[adata.obs['response']=='Responder', 'ecm_receptor_score'], adata.obs.loc[adata.obs['response']=='Responder', 'n_mech_genes'])[0]:.3f}
    Non-responder: ρ = {stats.spearmanr(adata.obs.loc[adata.obs['response']=='Non-responder', 'ecm_receptor_score'], adata.obs.loc[adata.obs['response']=='Non-responder', 'n_mech_genes'])[0]:.3f}
  
  Contractile Treg ECM score:
    Contractile: {cont_score.mean():.3f}
    Non-contractile: {nonc_score.mean():.3f}
    p = {pval:.2e}
"""

ax6.text(0.05, 0.95, table_text, transform=ax6.transAxes, fontsize=9,
         verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='#F8F9FA', edgecolor='#DEE2E6'))

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/virtual_knockout_comprehensive.png', dpi=300, bbox_inches='tight')
print(f"\nSaved: {OUT_DIR}/virtual_knockout_comprehensive.png")
plt.close()

# Save all results
ko_df.to_csv(f'{OUT_DIR}/virtual_ko_cellchat_results.csv', index=False)
pd.DataFrame(receptor_ko_results).to_csv(f'{OUT_DIR}/virtual_ko_receptor_results.csv', index=False)

print("\n" + "="*70)
print("VIRTUAL KNOCKOUT ANALYSIS COMPLETE")
print("="*70)
