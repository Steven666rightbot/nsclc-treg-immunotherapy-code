"""
SCENIC post-analysis: AUCell + differential analysis + heatmap + top TF
"""
import pandas as pd
import numpy as np
import loompy as lp
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
import os

print("=" * 60)
print("Step 1: Loading data")
print("=" * 60)

# Load expression matrix from loom
loom_path = 'results/scenic_input/treg_no_foxp3_8k.loom'
with lp.connect(loom_path, mode='r') as ds:
    genes = ds.ra.Gene
    cells = ds.ca.CellID
    cell_types = ds.ca.sub_cell_type
    responses = ds.ca.response
    exp_mtx = pd.DataFrame(ds[:, :], index=genes, columns=cells)
    
print(f"Expression matrix: {exp_mtx.shape}")
print(f"Cell types: {pd.Series(cell_types).value_counts().to_dict()}")
print(f"Response: {pd.Series(responses).value_counts().to_dict()}")

# Load regulons
print("\n" + "=" * 60)
print("Step 2: Loading regulons")
print("=" * 60)

df_motifs = load_motifs('results/scenic_input/scenic_treg_8k_regulons.csv')
regulons = df2regulons(df_motifs)
print(f"Regulons: {len(regulons)}")
for i, reg in enumerate(regulons[:5]):
    print(f"  {reg.name}: {len(reg.genes)} targets")

# Calculate AUCell
print("\n" + "=" * 60)
print("Step 3: Calculating AUCell")
print("=" * 60)

auc_mtx = aucell(exp_mtx.T, regulons, num_workers=1)
print(f"AUCell matrix: {auc_mtx.shape}")
print(f"AUCell range: {auc_mtx.min().min():.3f} - {auc_mtx.max().max():.3f}")

# Save AUCell matrix
auc_mtx.to_csv('results/scenic_input/aucell_matrix.csv')
print("Saved: results/scenic_input/aucell_matrix.csv")

# Add response label
auc_mtx['response'] = list(responses)

# Differential analysis
print("\n" + "=" * 60)
print("Step 4: Differential analysis (R vs NR)")
print("=" * 60)

results = []
for tf in auc_mtx.columns:
    if tf == 'response':
        continue
    r_vals = auc_mtx[auc_mtx['response'] == 'Responder'][tf]
    nr_vals = auc_mtx[auc_mtx['response'] == 'Non-responder'][tf]
    
    # Mann-Whitney U test
    _, pval = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
    
    # Effect size (mean difference)
    effect = r_vals.mean() - nr_vals.mean()
    
    results.append({
        'TF': tf,
        'R_mean': r_vals.mean(),
        'NR_mean': nr_vals.mean(),
        'effect_size': effect,
        'p_value': pval,
        '-log10_p': -np.log10(pval) if pval > 0 else 50
    })

diff_df = pd.DataFrame(results)
diff_df['padj'] = diff_df['p_value'] * len(diff_df)  # Bonferroni
diff_df = diff_df.sort_values('p_value')
print(diff_df.to_string(index=False))

# Top TF筛选
print("\n" + "=" * 60)
print("Step 5: Top TF regulons")
print("=" * 60)

top_tfs = diff_df[diff_df['p_value'] < 0.05].sort_values('effect_size', key=abs, ascending=False)
print(f"Significant regulons (p < 0.05): {len(top_tfs)}")
print(top_tfs[['TF', 'R_mean', 'NR_mean', 'effect_size', 'p_value']].to_string(index=False))

# Save results
diff_df.to_csv('results/scenic_input/regulon_differential_analysis.csv', index=False)
print("\nSaved: results/scenic_input/regulon_differential_analysis.csv")

# Heatmap
print("\n" + "=" * 60)
print("Step 6: Generating heatmap")
print("=" * 60)

# Select top significant regulons for heatmap
n_top = min(12, len(top_tfs))
if n_top > 0:
    top_tf_names = top_tfs['TF'].head(n_top).tolist()
    
    # Prepare data for heatmap
    heatmap_data = auc_mtx[top_tf_names + ['response']].copy()
    heatmap_data = heatmap_data.sort_values('response')
    
    # Z-score normalization per regulon
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    heatmap_z = pd.DataFrame(
        scaler.fit_transform(heatmap_data[top_tf_names]),
        index=heatmap_data.index,
        columns=top_tf_names
    )
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Color bar for response
    response_colors = {'Responder': '#3498DB', 'Non-responder': '#E74C3C'}
    row_colors = heatmap_data['response'].map(response_colors)
    
    sns.heatmap(heatmap_z, cmap='RdBu_r', center=0, 
                xticklabels=True, yticklabels=False,
                cbar_kws={'label': 'Z-score (AUCell)'},
                ax=ax)
    ax.set_title('Regulon Activity Heatmap (Top Differential TFs)', fontsize=14, fontweight='bold')
    ax.set_xlabel('Transcription Factor Regulons', fontsize=12)
    ax.set_ylabel('Cells', fontsize=12)
    
    # Add response legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#3498DB', label='Responder'),
                       Patch(facecolor='#E74C3C', label='Non-responder')]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.15, 1))
    
    plt.tight_layout()
    fig.savefig('figures/scenic_regulon_heatmap.png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig('figures/scenic_regulon_heatmap.pdf', dpi=300, bbox_inches='tight', facecolor='white')
    print("Saved: figures/scenic_regulon_heatmap.png")
else:
    print("No significant regulons found for heatmap.")

print("\n" + "=" * 60)
print("Analysis complete!")
print("=" * 60)
