import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Load metadata
obs = pd.read_csv(os.path.join(BASE_DIR, 'results', 'spatial', 'spatial_spot_metadata.csv'), index_col=0)

# Reload a single sample for spatial plots
sample_id = 'PDAC_1'
adata = sc.read_visium(os.path.join(BASE_DIR, 'data', 'spatial', 'pdac_visium', 'extracted', sample_id))
adata.var_names_make_unique()
adata = adata[adata.obs['in_tissue'] == 1].copy()

sample_mask = obs['sample_id'] == sample_id
obs_sample = obs[sample_mask].copy()
common_idx = list(set(adata.obs.index) & set(obs_sample.index))
adata = adata[common_idx].copy()

for col in ['stromal_score', 'tumor_score', 'region', 'ccr8_score', 'mki67_score', 'ccr8_score_z', 'mki67_score_z']:
    if col in obs_sample.columns:
        adata.obs[col] = obs_sample.loc[adata.obs.index, col]

print(f"Reloaded {sample_id}: {adata.shape}")
print(f"Region counts: {adata.obs['region'].value_counts().to_dict()}")

# Plot 1: Region definition
print("Plotting regions...")
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.spatial(adata, color='stromal_score', ax=axes[0], show=False, title='Stromal Score', cmap='RdYlBu_r', spot_size=50)
sc.pl.spatial(adata, color='tumor_score', ax=axes[1], show=False, title='Tumor Score', cmap='RdYlBu_r', spot_size=50)
sc.pl.spatial(adata, color='region', ax=axes[2], show=False, title='Region', palette={'hard_stroma': '#E41A1C', 'tumor_core': '#377EB8', 'other': '#999999'}, spot_size=50)
plt.tight_layout()
plt.savefig(os.path.join(BASE_DIR, 'figures', 'spatial_region_definition.png'), dpi=300, bbox_inches='tight')
plt.close()

# Plot 2: Treg scores
print("Plotting Treg scores...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
sc.pl.spatial(adata, color='ccr8_score', ax=axes[0], show=False, title='CCR8+ Treg Signature', cmap='YlOrRd', spot_size=50)
sc.pl.spatial(adata, color='mki67_score', ax=axes[1], show=False, title='MKI67+ Treg Signature', cmap='YlOrRd', spot_size=50)
plt.tight_layout()
plt.savefig(os.path.join(BASE_DIR, 'figures', 'spatial_treg_scores.png'), dpi=300, bbox_inches='tight')
plt.close()

# Plot 3: Boxplots (all samples)
print("Plotting boxplots...")
plot_data = []
for score_col, label in [('ccr8_score_z', 'CCR8+ Treg'), ('mki67_score_z', 'MKI67+ Treg')]:
    if score_col not in obs.columns:
        continue
    for region in ['hard_stroma', 'tumor_core', 'other']:
        mask = obs['region'] == region
        vals = obs.loc[mask, score_col].dropna().values
        for v in vals:
            plot_data.append({'Score': label, 'Region': region, 'Value': v})

if plot_data:
    df_plot = pd.DataFrame(plot_data)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.boxplot(data=df_plot, x='Score', y='Value', hue='Region',
                palette={'hard_stroma': '#E41A1C', 'tumor_core': '#377EB8', 'other': '#999999'}, ax=ax)
    ax.set_ylabel('Z-score (per-sample normalized)')
    ax.set_title('Treg Subtype Signature by Region')
    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(BASE_DIR, 'figures', 'spatial_treg_boxplot.png'), dpi=300, bbox_inches='tight')
    plt.close()

# Plot 4: Scatter
print("Plotting scatter...")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for idx, (score_col, title) in enumerate([('ccr8_score_z', 'CCR8+ vs Stroma'), ('mki67_score_z', 'MKI67+ vs Stroma')]):
    if score_col in obs.columns:
        x = obs['stromal_score'].values
        y = obs[score_col].values
        m = ~(np.isnan(x) | np.isnan(y))
        x, y = x[m], y[m]
        if len(x) > 0:
            r, p = stats.pearsonr(x, y)
            axes[idx].scatter(x, y, alpha=0.2, s=3)
            axes[idx].set_xlabel('Stromal Score')
            axes[idx].set_ylabel('Z-score')
            axes[idx].set_title(f'{title}\nr={r:.3f}, p={p:.2e}')
plt.tight_layout()
plt.savefig(os.path.join(BASE_DIR, 'figures', 'spatial_stroma_treg_scatter.png'), dpi=300, bbox_inches='tight')
plt.close()

print("All plots done!")
