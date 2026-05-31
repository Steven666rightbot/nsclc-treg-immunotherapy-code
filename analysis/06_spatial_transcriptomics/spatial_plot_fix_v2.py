import os
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'data', 'spatial', 'pdac_visium', 'extracted')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

# Load metadata
obs = pd.read_csv(os.path.join(BASE_DIR, 'results', 'spatial', 'spatial_spot_metadata.csv'), index_col=0)

# Reload single sample
sample_id = 'PDAC_1'
print(f"Reloading {sample_id}...")
adata = sc.read_visium(os.path.join(DATA_DIR, sample_id))
adata.var_names_make_unique()
adata = adata[adata.obs['in_tissue'] == 1].copy()

# Match metadata: strip suffix from merged obs indices
obs_sample = obs[obs['sample_id'] == sample_id].copy()
# Merged index is like 'AAACAAGTATCTCCCA-1-PDAC_1', original is 'AAACAAGTATCTCCCA-1'
# Remove the '-PDAC_1' suffix
obs_sample['orig_index'] = obs_sample.index.str.replace(f'-{sample_id}$', '', regex=True)
obs_sample = obs_sample.set_index('orig_index')

common_idx = list(set(adata.obs.index) & set(obs_sample.index))
adata = adata[common_idx].copy()
print(f"Matched {len(common_idx)} spots")

for col in ['stromal_score', 'tumor_score', 'region', 'ccr8_score', 'mki67_score', 'ccr8_score_z', 'mki67_score_z']:
    if col in obs_sample.columns:
        adata.obs[col] = obs_sample.loc[adata.obs.index, col]

print(f"Region counts: {adata.obs['region'].value_counts().to_dict()}")

# Plot 1: Region definition
print("Plotting regions...")
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.spatial(adata, color='stromal_score', ax=axes[0], show=False, title='Stromal Score', cmap='RdYlBu_r', spot_size=50)
sc.pl.spatial(adata, color='tumor_score', ax=axes[1], show=False, title='Tumor Score', cmap='RdYlBu_r', spot_size=50)
sc.pl.spatial(adata, color='region', ax=axes[2], show=False, title='Region', 
              palette={'hard_stroma': '#E41A1C', 'tumor_core': '#377EB8', 'other': '#999999'}, spot_size=50)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'spatial_region_definition.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved spatial_region_definition.png")

# Plot 2: Treg scores
print("Plotting Treg scores...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
sc.pl.spatial(adata, color='ccr8_score', ax=axes[0], show=False, title='CCR8+ Treg Signature', cmap='YlOrRd', spot_size=50)
sc.pl.spatial(adata, color='mki67_score', ax=axes[1], show=False, title='MKI67+ Treg Signature', cmap='YlOrRd', spot_size=50)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'spatial_treg_scores.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved spatial_treg_scores.png")

print("Done!")
