"""
Dorothea TF analysis: CCR8+ vs MKI67+ subtypes + SHAP cross-validation
"""
import pandas as pd
import numpy as np
import loompy as lp
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import decoupler as dc
import os

os.makedirs('figures', exist_ok=True)
os.makedirs('results/scenic_input', exist_ok=True)

# =============================================================================
# Load data
# =============================================================================
print("=" * 60)
print("Loading data")
print("=" * 60)

loom_path = 'results/scenic_input/treg_no_foxp3_8k.loom'
with lp.connect(loom_path, mode='r') as ds:
    genes = ds.ra.Gene
    cells = ds.ca.CellID
    cell_types = ds.ca.sub_cell_type
    responses = ds.ca.response
    exp_mtx = pd.DataFrame(ds[:, :].T, index=cells, columns=genes)

# Load Dorothea regulons
net = dc.op.dorothea(organism='human', levels=['A', 'B'])
net_filt = net[net['target'].isin(exp_mtx.columns)]
print(f"Dorothea: {net_filt['source'].nunique()} TFs, {len(net_filt)} interactions")

# Load pre-computed TF activities (from previous run)
tf_acts = pd.read_csv('results/scenic_input/dorothea_tf_activity.csv', index_col=0)
print(f"TF activity matrix: {tf_acts.shape}")

# Load SHAP importance
shap_df = pd.read_csv('data/shap_importance.csv')
print(f"SHAP features: {len(shap_df)}")

# =============================================================================
# Part 1: CCR8+ vs MKI67+ TF activity differences
# =============================================================================
print("\n" + "=" * 60)
print("Part 1: Differential TF activity (CCR8+ vs MKI67+)")
print("=" * 60)

tf_acts['sub_cell_type'] = list(cell_types)

results_subtype = []
for tf in tf_acts.columns:
    if tf == 'sub_cell_type':
        continue
    ccr8_vals = tf_acts[tf_acts['sub_cell_type'] == 'CD4T_Treg_CCR8'][tf]
    mki67_vals = tf_acts[tf_acts['sub_cell_type'] == 'CD4T_Treg_MKI67'][tf]
    
    _, pval = stats.mannwhitneyu(ccr8_vals, mki67_vals, alternative='two-sided')
    effect = ccr8_vals.mean() - mki67_vals.mean()
    
    results_subtype.append({
        'TF': tf,
        'CCR8_mean': ccr8_vals.mean(),
        'MKI67_mean': mki67_vals.mean(),
        'effect_size': effect,
        'p_value': pval,
        '-log10_p': -np.log10(pval) if pval > 0 else 50
    })

subtype_df = pd.DataFrame(results_subtype)
subtype_df['padj'] = subtype_df['p_value'] * len(subtype_df)
subtype_df = subtype_df.sort_values('p_value')
print(subtype_df.head(20).to_string(index=False))

n_sig = (subtype_df['p_value'] < 0.05).sum()
print(f"\nSignificant TFs (p < 0.05): {n_sig}")
subtype_df.to_csv('results/scenic_input/dorothea_differential_ccr8_vs_mki67.csv', index=False)

# Volcano plot
fig, ax = plt.subplots(figsize=(8, 8))
subtype_df['sig'] = 'NS'
subtype_df.loc[(subtype_df['p_value'] < 0.05) & (subtype_df['effect_size'] > 0), 'sig'] = 'Up in CCR8+'
subtype_df.loc[(subtype_df['p_value'] < 0.05) & (subtype_df['effect_size'] < 0), 'sig'] = 'Up in MKI67+'

colors = {'NS': '#999999', 'Up in CCR8+': '#2E86AB', 'Up in MKI67+': '#A23B72'}
for label, grp in subtype_df.groupby('sig'):
    ax.scatter(grp['effect_size'], -np.log10(grp['p_value']),
               c=colors[label], s=30, alpha=0.7, edgecolors='none', label=label)

# Label top
top_to_label = subtype_df.nsmallest(15, 'p_value')
for _, row in top_to_label.iterrows():
    ax.annotate(row['TF'], (row['effect_size'], -np.log10(row['p_value'])),
                fontsize=8, xytext=(5, 5), textcoords='offset points')

ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8, alpha=0.5)
ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
ax.set_xlabel('Effect size (CCR8+ mean - MKI67+ mean)', fontsize=12)
ax.set_ylabel('-log10(p-value)', fontsize=12)
ax.set_title('Differential TF Activity: CCR8+ vs MKI67+ Treg', fontsize=14, fontweight='bold')
ax.legend(loc='upper right')
plt.tight_layout()
fig.savefig('figures/dorothea_volcano_ccr8_vs_mki67.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig('figures/dorothea_volcano_ccr8_vs_mki67.pdf', dpi=300, bbox_inches='tight', facecolor='white')
print("Saved: figures/dorothea_volcano_ccr8_vs_mki67.png")

# Heatmap
n_top = min(30, len(subtype_df))
top_tfs = subtype_df.nsmallest(n_top, 'p_value')['TF'].tolist()
heatmap_data = tf_acts[top_tfs + ['sub_cell_type']].copy()
heatmap_data = heatmap_data.sort_values('sub_cell_type')

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
heatmap_z = pd.DataFrame(
    scaler.fit_transform(heatmap_data[top_tfs]),
    index=heatmap_data.index,
    columns=top_tfs
)

# Downsample
ccr8_cells = heatmap_data[heatmap_data['sub_cell_type'] == 'CD4T_Treg_CCR8'].index
mki67_cells = heatmap_data[heatmap_data['sub_cell_type'] == 'CD4T_Treg_MKI67'].index
if len(ccr8_cells) > 1000:
    ccr8_cells = np.random.choice(ccr8_cells, 1000, replace=False)
if len(mki67_cells) > 1000:
    mki67_cells = np.random.choice(mki67_cells, 1000, replace=False)
plot_cells = list(ccr8_cells) + list(mki67_cells)

heatmap_z_plot = heatmap_z.loc[plot_cells]
subtype_plot = heatmap_data.loc[plot_cells, 'sub_cell_type']

fig, ax = plt.subplots(figsize=(12, 10))
subtype_colors = {'CD4T_Treg_CCR8': '#2E86AB', 'CD4T_Treg_MKI67': '#A23B72'}
row_colors = subtype_plot.map(subtype_colors)

sns.heatmap(heatmap_z_plot, cmap='RdBu_r', center=0,
            xticklabels=True, yticklabels=False,
            cbar_kws={'label': 'Z-score (TF activity)'},
            ax=ax)
ax.set_title('Top Differential TF Activities: CCR8+ vs MKI67+', fontsize=14, fontweight='bold')
ax.set_xlabel('Transcription Factors', fontsize=12)
ax.set_ylabel('Cells', fontsize=12)

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#2E86AB', label='CCR8+'),
                   Patch(facecolor='#A23B72', label='MKI67+')]
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.15, 1))

plt.tight_layout()
fig.savefig('figures/dorothea_heatmap_ccr8_vs_mki67.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig('figures/dorothea_heatmap_ccr8_vs_mki67.pdf', dpi=300, bbox_inches='tight', facecolor='white')
print("Saved: figures/dorothea_heatmap_ccr8_vs_mki67.png")

# =============================================================================
# Part 2: SHAP cross-validation
# =============================================================================
print("\n" + "=" * 60)
print("Part 2: SHAP x Dorothea cross-validation")
print("=" * 60)

# Extract gene symbols from SHAP features
# Format: "CellType_Gene" or "CellType_Gene+Gene2" or just "CellType"
def extract_genes(feature_name):
    """Extract gene symbols from SHAP feature names."""
    # Skip features without underscore (whole cell types)
    if '_' not in feature_name:
        return []
    # Take the last part after the last underscore (gene symbol)
    gene_part = feature_name.rsplit('_', 1)[-1]
    # Handle cases like "GZMK+GZMH+" -> ['GZMK', 'GZMH']
    genes = []
    for g in gene_part.split('+'):
        g = g.strip()
        if g and g not in ['cell', 'Cell']:
            genes.append(g)
    return genes

shap_genes = set()
for feat in shap_df['feature']:
    for g in extract_genes(feat):
        shap_genes.add(g)

print(f"Extracted {len(shap_genes)} unique genes from SHAP features")
print(f"Examples: {list(shap_genes)[:10]}")

# For each TF in Dorothea, check overlap with SHAP genes
crossval_results = []
for tf in net_filt['source'].unique():
    tf_targets = set(net_filt[net_filt['source'] == tf]['target'])
    overlap = tf_targets & shap_genes
    overlap_frac = len(overlap) / len(tf_targets) if len(tf_targets) > 0 else 0
    
    crossval_results.append({
        'TF': tf,
        'n_targets': len(tf_targets),
        'n_shap_overlap': len(overlap),
        'shap_overlap_frac': overlap_frac,
        'shap_genes': ';'.join(sorted(overlap)) if overlap else ''
    })

crossval_df = pd.DataFrame(crossval_results)

# Merge with R vs NR differential results
diff_rnr = pd.read_csv('results/scenic_input/dorothea_differential_tfs.csv')
crossval_df = crossval_df.merge(diff_rnr[['TF', 'p_value', 'effect_size']], on='TF', how='left')
crossval_df = crossval_df.rename(columns={'p_value': 'rnr_p', 'effect_size': 'rnr_effect'})

# Merge with CCR8 vs MKI67 differential results
crossval_df = crossval_df.merge(subtype_df[['TF', 'p_value', 'effect_size']], on='TF', how='left')
crossval_df = crossval_df.rename(columns={'p_value': 'subtype_p', 'effect_size': 'subtype_effect'})

crossval_df = crossval_df.sort_values('n_shap_overlap', ascending=False)
print("\nTop TFs by SHAP gene overlap:")
print(crossval_df.head(15)[['TF', 'n_targets', 'n_shap_overlap', 'shap_overlap_frac', 'rnr_p', 'subtype_p']].to_string(index=False))

# Significant TFs with SHAP overlap
sig_rnr = crossval_df[crossval_df['rnr_p'] < 0.05]
sig_subtype = crossval_df[crossval_df['subtype_p'] < 0.05]
print(f"\nSignificant R vs NR TFs with SHAP overlap: {(sig_rnr['n_shap_overlap'] > 0).sum()}/{len(sig_rnr)}")
print(f"Significant CCR8 vs MKI67 TFs with SHAP overlap: {(sig_subtype['n_shap_overlap'] > 0).sum()}/{len(sig_subtype)}")

# Top overlapping TFs that are also significant
sig_with_overlap = crossval_df[(crossval_df['rnr_p'] < 0.05) & (crossval_df['n_shap_overlap'] > 0)]
sig_with_overlap = sig_with_overlap.sort_values('n_shap_overlap', ascending=False)
print("\nSignificant R/NR TFs with SHAP overlap (top 10):")
print(sig_with_overlap.head(10)[['TF', 'n_shap_overlap', 'shap_overlap_frac', 'rnr_effect', 'rnr_p']].to_string(index=False))

# Save
crossval_df.to_csv('results/scenic_input/dorothea_shap_crossvalidation.csv', index=False)
print("\nSaved: results/scenic_input/dorothea_shap_crossvalidation.csv")

# =============================================================================
# Cross-validation visualization
# =============================================================================
print("\n" + "=" * 60)
print("Part 3: Cross-validation plots")
print("=" * 60)

# Plot 1: Scatter -log10(p) vs SHAP overlap
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# R vs NR
ax = axes[0]
colors = np.where(crossval_df['n_shap_overlap'] > 0, '#E74C3C', '#999999')
ax.scatter(crossval_df['n_shap_overlap'], -np.log10(crossval_df['rnr_p'].clip(lower=1e-300)),
           c=colors, s=30, alpha=0.7, edgecolors='none')
ax.set_xlabel('Number of SHAP gene overlaps', fontsize=12)
ax.set_ylabel('-log10(p-value) R vs NR', fontsize=12)
ax.set_title('Dorothea TFs: R/NR significance vs SHAP overlap', fontsize=13, fontweight='bold')
ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8, alpha=0.5)

# Label top points
top_label = crossval_df[(crossval_df['rnr_p'] < 0.05) & (crossval_df['n_shap_overlap'] > 0)].nsmallest(8, 'rnr_p')
for _, row in top_label.iterrows():
    ax.annotate(row['TF'], (row['n_shap_overlap'], -np.log10(row['rnr_p'])),
                fontsize=8, xytext=(5, 5), textcoords='offset points')

# CCR8 vs MKI67
ax = axes[1]
colors = np.where(crossval_df['n_shap_overlap'] > 0, '#A23B72', '#999999')
ax.scatter(crossval_df['n_shap_overlap'], -np.log10(crossval_df['subtype_p'].clip(lower=1e-300)),
           c=colors, s=30, alpha=0.7, edgecolors='none')
ax.set_xlabel('Number of SHAP gene overlaps', fontsize=12)
ax.set_ylabel('-log10(p-value) CCR8+ vs MKI67+', fontsize=12)
ax.set_title('Dorothea TFs: Subtype significance vs SHAP overlap', fontsize=13, fontweight='bold')
ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8, alpha=0.5)

top_label = crossval_df[(crossval_df['subtype_p'] < 0.05) & (crossval_df['n_shap_overlap'] > 0)].nsmallest(8, 'subtype_p')
for _, row in top_label.iterrows():
    ax.annotate(row['TF'], (row['n_shap_overlap'], -np.log10(row['subtype_p'])),
                fontsize=8, xytext=(5, 5), textcoords='offset points')

plt.tight_layout()
fig.savefig('figures/dorothea_shap_crossval_scatter.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig('figures/dorothea_shap_crossval_scatter.pdf', dpi=300, bbox_inches='tight', facecolor='white')
print("Saved: figures/dorothea_shap_crossval_scatter.png")

# Plot 2: Bar plot of top TFs with SHAP overlap
fig, axes = plt.subplots(1, 2, figsize=(14, 8))

# Top R/NR TFs with SHAP overlap
top_rnr = sig_with_overlap.head(15)
ax = axes[0]
colors = ['#E74C3C' if e > 0 else '#3498DB' for e in top_rnr['rnr_effect']]
ax.barh(range(len(top_rnr)), top_rnr['n_shap_overlap'], color=colors, alpha=0.7)
ax.set_yticks(range(len(top_rnr)))
ax.set_yticklabels(top_rnr['TF'], fontsize=10)
ax.set_xlabel('Number of SHAP gene overlaps', fontsize=12)
ax.set_title('Top R/NR TFs with SHAP overlap', fontsize=13, fontweight='bold')
ax.invert_yaxis()

# Top Subtype TFs with SHAP overlap
sig_sub_overlap = crossval_df[(crossval_df['subtype_p'] < 0.05) & (crossval_df['n_shap_overlap'] > 0)]
sig_sub_overlap = sig_sub_overlap.sort_values('n_shap_overlap', ascending=False).head(15)
ax = axes[1]
colors = ['#2E86AB' if e > 0 else '#A23B72' for e in sig_sub_overlap['subtype_effect']]
ax.barh(range(len(sig_sub_overlap)), sig_sub_overlap['n_shap_overlap'], color=colors, alpha=0.7)
ax.set_yticks(range(len(sig_sub_overlap)))
ax.set_yticklabels(sig_sub_overlap['TF'], fontsize=10)
ax.set_xlabel('Number of SHAP gene overlaps', fontsize=12)
ax.set_title('Top Subtype TFs with SHAP overlap', fontsize=13, fontweight='bold')
ax.invert_yaxis()

plt.tight_layout()
fig.savefig('figures/dorothea_shap_crossval_bar.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig('figures/dorothea_shap_crossval_bar.pdf', dpi=300, bbox_inches='tight', facecolor='white')
print("Saved: figures/dorothea_shap_crossval_bar.png")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"CCR8+ vs MKI67+ significant TFs: {n_sig}")
print(f"Top CCR8+ TF: {subtype_df.iloc[0]['TF']} (effect={subtype_df.iloc[0]['effect_size']:.3f}, p={subtype_df.iloc[0]['p_value']:.2e})")
print(f"Top MKI67+ TF: {subtype_df.sort_values('effect_size').iloc[0]['TF']} (effect={subtype_df.sort_values('effect_size').iloc[0]['effect_size']:.3f})")
print()
print(f"SHAP genes extracted: {len(shap_genes)}")
print(f"Dorothea TFs with SHAP overlap: {(crossval_df['n_shap_overlap'] > 0).sum()}/{len(crossval_df)}")
print(f"Significant R/NR TFs with SHAP overlap: {len(sig_with_overlap)}")
print(f"Significant Subtype TFs with SHAP overlap: {len(sig_sub_overlap)}")
print("\nAll outputs saved to results/scenic_input/ and figures/")
print("=" * 60)
