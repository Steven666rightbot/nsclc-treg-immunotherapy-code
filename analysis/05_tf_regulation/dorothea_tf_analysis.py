"""
Dorothea-based TF activity analysis (alternative to SCENIC)
Uses curated literature-based regulons instead of de novo inference.
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

print("=" * 60)
print("Step 1: Loading expression matrix")
print("=" * 60)

loom_path = 'results/scenic_input/treg_no_foxp3_8k.loom'
with lp.connect(loom_path, mode='r') as ds:
    genes = ds.ra.Gene
    cells = ds.ca.CellID
    cell_types = ds.ca.sub_cell_type
    responses = ds.ca.response
    # decoupler expects cells x genes
    exp_mtx = pd.DataFrame(ds[:, :].T, index=cells, columns=genes)

print(f"Expression matrix: {exp_mtx.shape}")
print(f"Cell types: {pd.Series(cell_types).value_counts().to_dict()}")
print(f"Response: {pd.Series(responses).value_counts().to_dict()}")

print("\n" + "=" * 60)
print("Step 2: Loading Dorothea regulons (A+B confidence)")
print("=" * 60)

net = dc.op.dorothea(organism='human', levels=['A', 'B'])
print(f"Regulons: {net['source'].nunique()} TFs, {len(net)} interactions")
print(f"Confidence: A={sum(net['confidence']=='A')}, B={sum(net['confidence']=='B')}")

# Filter to TFs present in expression matrix
net_filt = net[net['target'].isin(exp_mtx.columns)]
print(f"After filtering to expressed genes: {net_filt['source'].nunique()} TFs, {len(net_filt)} interactions")

print("\n" + "=" * 60)
print("Step 3: Running ULM (Univariate Linear Model) for TF activity")
print("=" * 60)

# ULM estimates TF activity per cell
from decoupler.mt._ulm import ulm
tf_acts, tf_pvals = ulm(exp_mtx, net_filt, verbose=True)
print(f"TF activity matrix: {tf_acts.shape}")
print(f"Activity range: {tf_acts.min().min():.3f} - {tf_acts.max().max():.3f}")

# Save
tf_acts.to_csv('results/scenic_input/dorothea_tf_activity.csv')
print("Saved: results/scenic_input/dorothea_tf_activity.csv")

print("\n" + "=" * 60)
print("Step 4: Differential TF activity (Responder vs Non-responder)")
print("=" * 60)

# Add response label
tf_acts['response'] = list(responses)

results = []
for tf in tf_acts.columns:
    if tf == 'response':
        continue
    r_vals = tf_acts[tf_acts['response'] == 'Responder'][tf]
    nr_vals = tf_acts[tf_acts['response'] == 'Non-responder'][tf]
    
    # Mann-Whitney U test
    _, pval = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
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
print(diff_df.head(20).to_string(index=False))

n_sig = (diff_df['p_value'] < 0.05).sum()
print(f"\nSignificant TFs (p < 0.05): {n_sig}")

# Save
diff_df.to_csv('results/scenic_input/dorothea_differential_tfs.csv', index=False)
print("Saved: results/scenic_input/dorothea_differential_tfs.csv")

print("\n" + "=" * 60)
print("Step 5: Volcano plot")
print("=" * 60)

fig, ax = plt.subplots(figsize=(8, 8))

# Color by significance
diff_df['sig'] = 'NS'
diff_df.loc[(diff_df['p_value'] < 0.05) & (diff_df['effect_size'] > 0), 'sig'] = 'Up in R'
diff_df.loc[(diff_df['p_value'] < 0.05) & (diff_df['effect_size'] < 0), 'sig'] = 'Up in NR'

colors = {'NS': '#999999', 'Up in R': '#E74C3C', 'Up in NR': '#3498DB'}
for label, grp in diff_df.groupby('sig'):
    ax.scatter(grp['effect_size'], -np.log10(grp['p_value']),
               c=colors[label], s=30, alpha=0.7, edgecolors='none', label=label)

# Label top TFs
top_to_label = diff_df.nsmallest(10, 'p_value')
for _, row in top_to_label.iterrows():
    ax.annotate(row['TF'], (row['effect_size'], -np.log10(row['p_value'])),
                fontsize=8, xytext=(5, 5), textcoords='offset points')

ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8, alpha=0.5)
ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
ax.set_xlabel('Effect size (R mean - NR mean)', fontsize=12)
ax.set_ylabel('-log10(p-value)', fontsize=12)
ax.set_title('Differential TF Activity (Dorothea A+B)', fontsize=14, fontweight='bold')
ax.legend(loc='upper right')
plt.tight_layout()
fig.savefig('figures/dorothea_volcano.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig('figures/dorothea_volcano.pdf', dpi=300, bbox_inches='tight', facecolor='white')
print("Saved: figures/dorothea_volcano.png")

print("\n" + "=" * 60)
print("Step 6: Heatmap of top differential TFs")
print("=" * 60)

n_top = min(30, len(diff_df))
if n_top > 0:
    top_tfs = diff_df.nsmallest(n_top, 'p_value')['TF'].tolist()
    
    # Prepare data
    heatmap_data = tf_acts[top_tfs + ['response']].copy()
    heatmap_data = heatmap_data.sort_values('response')
    
    # Z-score per TF
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    heatmap_z = pd.DataFrame(
        scaler.fit_transform(heatmap_data[top_tfs]),
        index=heatmap_data.index,
        columns=top_tfs
    )
    
    # Downsample cells for visibility (max 1000 per group)
    r_cells = heatmap_data[heatmap_data['response'] == 'Responder'].index
    nr_cells = heatmap_data[heatmap_data['response'] == 'Non-responder'].index
    if len(r_cells) > 1000:
        r_cells = np.random.choice(r_cells, 1000, replace=False)
    if len(nr_cells) > 1000:
        nr_cells = np.random.choice(nr_cells, 1000, replace=False)
    plot_cells = list(r_cells) + list(nr_cells)
    
    heatmap_z_plot = heatmap_z.loc[plot_cells]
    response_plot = heatmap_data.loc[plot_cells, 'response']
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 10))
    response_colors = {'Responder': '#3498DB', 'Non-responder': '#E74C3C'}
    row_colors = response_plot.map(response_colors)
    
    sns.heatmap(heatmap_z_plot, cmap='RdBu_r', center=0,
                xticklabels=True, yticklabels=False,
                cbar_kws={'label': 'Z-score (TF activity)'},
                ax=ax)
    ax.set_title('Top Differential TF Activities (Dorothea A+B)', fontsize=14, fontweight='bold')
    ax.set_xlabel('Transcription Factors', fontsize=12)
    ax.set_ylabel('Cells', fontsize=12)
    
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#3498DB', label='Responder'),
                       Patch(facecolor='#E74C3C', label='Non-responder')]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.15, 1))
    
    plt.tight_layout()
    fig.savefig('figures/dorothea_heatmap.png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig('figures/dorothea_heatmap.pdf', dpi=300, bbox_inches='tight', facecolor='white')
    print("Saved: figures/dorothea_heatmap.png")
else:
    print("No significant TFs for heatmap.")

print("\n" + "=" * 60)
print("Analysis complete!")
print("=" * 60)
