"""
Spatial Treg subtype mapping on PDAC Visium data (optimized version)
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import pickle
from scipy import stats
from itertools import combinations

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'data', 'spatial', 'pdac_visium', 'extracted')
RESULTS_DIR = os.path.join(BASE_DIR, 'results', 'spatial')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

SAMPLES = {
    'LG_1': 'Low Grade', 'LG_2': 'Low Grade', 'LG_3': 'Low Grade',
    'LG_4': 'Low Grade', 'LG_5': 'Low Grade', 'LG_6': 'Low Grade', 'LG_7': 'Low Grade',
    'HG_1': 'High Grade', 'HG_2': 'High Grade', 'HG_3': 'High Grade',
    'PDAC_1': 'PDAC', 'PDAC_2': 'PDAC', 'PDAC_3': 'PDAC',
}

# ==============================================================================
# STEP 1: Load all Visium samples
# ==============================================================================
print("=" * 60)
print("STEP 1: Loading Visium samples")
print("=" * 60)

adatas = []
for sample_id, grade in SAMPLES.items():
    sample_path = os.path.join(DATA_DIR, sample_id)
    print(f"  {sample_id}...", end=" ", flush=True)
    adata = sc.read_visium(sample_path)
    adata.var_names_make_unique()
    adata.obs['sample_id'] = sample_id
    adata.obs['grade'] = grade
    adata = adata[adata.obs['in_tissue'] == 1].copy()
    print(f"{adata.n_obs} spots")
    adatas.append(adata)

print(f"\nMerging {len(adatas)} samples...")
adata_spatial = anndata.concat(
    adatas, label='sample_id', keys=list(SAMPLES.keys()),
    index_unique='-', join='outer', merge='same'
)
# Restore spatial coords
for adata in adatas:
    sid = adata.obs['sample_id'].iloc[0]
    mask = adata_spatial.obs['sample_id'] == sid
    for col in ['array_row', 'array_col']:
        if col in adata.obs.columns:
            adata_spatial.obs.loc[mask, col] = adata.obs[col].values

print(f"Merged: {adata_spatial.shape}")

# ==============================================================================
# STEP 2: Preprocessing
# ==============================================================================
print("\n" + "=" * 60)
print("STEP 2: Preprocessing")
print("=" * 60)

sc.pp.filter_genes(adata_spatial, min_cells=10)
print(f"After filter: {adata_spatial.shape}")
sc.pp.normalize_total(adata_spatial, target_sum=1e4)
sc.pp.log1p(adata_spatial)

# ==============================================================================
# STEP 3: Define regions
# ==============================================================================
print("\n" + "=" * 60)
print("STEP 3: Defining regions")
print("=" * 60)

STROMAL_MARKERS = ['POSTN', 'COL1A1', 'COL1A2', 'ACTA2', 'TAGLN', 'FAP']
TUMOR_MARKERS = ['EPCAM', 'KRT19', 'KRT8', 'KRT18']

avail_stromal = [g for g in STROMAL_MARKERS if g in adata_spatial.var_names]
avail_tumor = [g for g in TUMOR_MARKERS if g in adata_spatial.var_names]
print(f"  Stromal markers available: {avail_stromal}")
print(f"  Tumor markers available: {avail_tumor}")

sc.tl.score_genes(adata_spatial, gene_list=avail_stromal, score_name='stromal_score')
sc.tl.score_genes(adata_spatial, gene_list=avail_tumor, score_name='tumor_score')

# Define regions
adata_spatial.obs['region'] = 'other'
for sample_id in SAMPLES.keys():
    mask = adata_spatial.obs['sample_id'] == sample_id
    if mask.sum() == 0:
        continue
    s_vals = adata_spatial.obs.loc[mask, 'stromal_score'].values
    t_vals = adata_spatial.obs.loc[mask, 'tumor_score'].values
    s_med, s_std = np.median(s_vals), np.std(s_vals)
    t_med, t_std = np.median(t_vals), np.std(t_vals)
    s_thresh = s_med + 0.5 * s_std
    t_thresh = t_med + 0.5 * t_std
    regions = []
    for s, t in zip(s_vals, t_vals):
        if s > s_thresh and s > t:
            regions.append('hard_stroma')
        elif t > t_thresh and t > s:
            regions.append('tumor_core')
        else:
            regions.append('other')
    adata_spatial.obs.loc[mask, 'region'] = regions
    counts = pd.Series(regions).value_counts().to_dict()
    print(f"    {sample_id}: {counts}")

# ==============================================================================
# STEP 4: Load precomputed Treg signatures
# ==============================================================================
print("\n" + "=" * 60)
print("STEP 4: Loading Treg signatures")
print("=" * 60)

sig_path = os.path.join(RESULTS_DIR, 'treg_signatures.pkl')
with open(sig_path, 'rb') as f:
    signatures = pickle.load(f)

# Filter to available genes and add canonical markers
for subtype in list(signatures.keys()):
    genes = signatures[subtype]
    available = [g for g in genes if g in adata_spatial.var_names]
    signatures[subtype] = available
    print(f"  {subtype}: {len(available)}/{len(genes)} markers in spatial data")

for gene in ['CCR8', 'LAYN', 'TNFRSF18', 'BATF', 'CTLA4']:
    if gene in adata_spatial.var_names and gene not in signatures['CD4T_Treg_CCR8']:
        signatures['CD4T_Treg_CCR8'].append(gene)
for gene in ['MKI67', 'TOP2A', 'STMN1', 'TYMS', 'PCNA']:
    if gene in adata_spatial.var_names and gene not in signatures['CD4T_Treg_MKI67']:
        signatures['CD4T_Treg_MKI67'].append(gene)

# ==============================================================================
# STEP 5: Score spots
# ==============================================================================
print("\n" + "=" * 60)
print("STEP 5: Scoring spots")
print("=" * 60)

for subtype, genes in signatures.items():
    score_name = subtype.replace('CD4T_Treg_', '').lower() + '_score'
    if len(genes) >= 3:
        sc.tl.score_genes(adata_spatial, gene_list=genes, score_name=score_name)
        print(f"  {subtype} -> {score_name} ({len(genes)} genes)")
    else:
        print(f"  WARNING: {subtype} skipped (only {len(genes)} genes)")

# Z-score per sample
for score_col in ['ccr8_score', 'mki67_score']:
    if score_col in adata_spatial.obs.columns:
        z_all = np.full(adata_spatial.n_obs, np.nan)
        for sample_id in SAMPLES.keys():
            mask = adata_spatial.obs['sample_id'] == sample_id
            vals = adata_spatial.obs.loc[mask, score_col].values
            if len(vals) > 0 and np.std(vals) > 0:
                z_all[mask.values] = (vals - np.mean(vals)) / np.std(vals)
        adata_spatial.obs[score_col + '_z'] = z_all

# ==============================================================================
# STEP 6: Statistics
# ==============================================================================
print("\n" + "=" * 60)
print("STEP 6: Statistics")
print("=" * 60)

results = []
for score_col in ['ccr8_score', 'mki67_score', 'ccr8_score_z', 'mki67_score_z']:
    if score_col not in adata_spatial.obs.columns:
        continue
    print(f"\n  {score_col}:")
    regions = ['hard_stroma', 'tumor_core', 'other']
    vals_by_region = {r: adata_spatial.obs.loc[adata_spatial.obs['region'] == r, score_col].dropna().values for r in regions}
    for r in regions:
        v = vals_by_region[r]
        print(f"    {r}: mean={np.mean(v):.4f}, n={len(v)}")
    for r1, r2 in combinations(regions, 2):
        v1, v2 = vals_by_region[r1], vals_by_region[r2]
        if len(v1) > 0 and len(v2) > 0:
            stat, pval = stats.mannwhitneyu(v1, v2, alternative='two-sided')
            print(f"    {r1} vs {r2}: p={pval:.2e}")
            results.append({
                'score': score_col, 'comparison': f'{r1}_vs_{r2}',
                'mean1': np.mean(v1), 'mean2': np.mean(v2),
                'n1': len(v1), 'n2': len(v2), 'pvalue': pval
            })

pd.DataFrame(results).to_csv(os.path.join(RESULTS_DIR, 'spatial_treg_stats.csv'), index=False)

# ==============================================================================
# STEP 7: Save data BEFORE plotting
# ==============================================================================
print("\n" + "=" * 60)
print("STEP 7: Saving data")
print("=" * 60)

adata_spatial.obs.to_csv(os.path.join(RESULTS_DIR, 'spatial_spot_metadata.csv'))
print("  Metadata saved")

# Save lighter h5ad
adata_light = adata_spatial.copy()
if 'spatial' in adata_light.uns:
    for sid in list(adata_light.uns.get('spatial', {}).keys()):
        if 'images' in adata_light.uns['spatial'].get(sid, {}):
            adata_light.uns['spatial'][sid]['images'] = {}
try:
    adata_light.write(os.path.join(RESULTS_DIR, 'spatial_processed.h5ad'))
    print("  h5ad saved")
except Exception as e:
    print(f"  h5ad save failed: {e}")

# Save signatures
sig_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in signatures.items()]))
sig_df.to_csv(os.path.join(RESULTS_DIR, 'treg_subtype_signatures.csv'), index=False)
print("  Signatures saved")

# ==============================================================================
# STEP 8: Plotting (reload single sample for spatial plots)
# ==============================================================================
print("\n" + "=" * 60)
print("STEP 8: Plotting")
print("=" * 60)

import seaborn as sns

# -- 8A. Spatial plots from single sample --
rep_sample = 'PDAC_1'
print(f"  8A: Reloading {rep_sample} for spatial plots...")
adata_rep = sc.read_visium(os.path.join(DATA_DIR, rep_sample))
adata_rep.var_names_make_unique()
adata_rep = adata_rep[adata_rep.obs['in_tissue'] == 1].copy()

# Match metadata
obs_rep = adata_spatial.obs[adata_spatial.obs['sample_id'] == rep_sample].copy()
common_idx = list(set(adata_rep.obs.index) & set(obs_rep.index))
adata_rep = adata_rep[common_idx].copy()
for col in ['stromal_score', 'tumor_score', 'region', 'ccr8_score', 'mki67_score', 'ccr8_score_z', 'mki67_score_z']:
    if col in obs_rep.columns:
        adata_rep.obs[col] = obs_rep.loc[adata_rep.obs.index, col]

print(f"    Matched {len(common_idx)} spots")

# Region definition plot
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.spatial(adata_rep, color='stromal_score', ax=axes[0], show=False, title='Stromal Score', cmap='RdYlBu_r', spot_size=50)
sc.pl.spatial(adata_rep, color='tumor_score', ax=axes[1], show=False, title='Tumor Score', cmap='RdYlBu_r', spot_size=50)
sc.pl.spatial(adata_rep, color='region', ax=axes[2], show=False, title='Region', palette={'hard_stroma': '#E41A1C', 'tumor_core': '#377EB8', 'other': '#999999'}, spot_size=50)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'spatial_region_definition.png'), dpi=300, bbox_inches='tight')
plt.close()
print("    Saved spatial_region_definition.png")

# Treg score plots
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
sc.pl.spatial(adata_rep, color='ccr8_score', ax=axes[0], show=False, title='CCR8+ Treg Signature', cmap='YlOrRd', spot_size=50)
sc.pl.spatial(adata_rep, color='mki67_score', ax=axes[1], show=False, title='MKI67+ Treg Signature', cmap='YlOrRd', spot_size=50)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'spatial_treg_scores.png'), dpi=300, bbox_inches='tight')
plt.close()
print("    Saved spatial_treg_scores.png")

# -- 8B. Boxplot (all samples) --
print("  8B: Boxplot...")
plot_data = []
for score_col, label in [('ccr8_score_z', 'CCR8+ Treg'), ('mki67_score_z', 'MKI67+ Treg')]:
    if score_col not in adata_spatial.obs.columns:
        continue
    for region in ['hard_stroma', 'tumor_core', 'other']:
        mask = adata_spatial.obs['region'] == region
        vals = adata_spatial.obs.loc[mask, score_col].dropna().values
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
    plt.savefig(os.path.join(FIGURES_DIR, 'spatial_treg_boxplot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print("    Saved spatial_treg_boxplot.png")

# -- 8C. Per-sample bar chart --
print("  8C: Per-sample summary...")
sample_summary = []
for sample_id in SAMPLES.keys():
    mask = adata_spatial.obs['sample_id'] == sample_id
    if mask.sum() == 0:
        continue
    for region in ['hard_stroma', 'tumor_core', 'other']:
        rmask = mask & (adata_spatial.obs['region'] == region)
        if rmask.sum() == 0:
            continue
        sample_summary.append({
            'sample': sample_id,
            'region': region,
            'ccr8_z': adata_spatial.obs.loc[rmask, 'ccr8_score_z'].mean(),
            'mki67_z': adata_spatial.obs.loc[rmask, 'mki67_score_z'].mean(),
            'n_spots': rmask.sum()
        })

df_summary = pd.DataFrame(sample_summary)
df_summary.to_csv(os.path.join(RESULTS_DIR, 'spatial_per_sample_summary.csv'), index=False)

# Plot per-sample heatmap
if len(df_summary) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    pivot_ccr8 = df_summary.pivot(index='sample', columns='region', values='ccr8_z')
    pivot_mki67 = df_summary.pivot(index='sample', columns='region', values='mki67_z')
    
    # Reorder columns
    col_order = ['hard_stroma', 'tumor_core', 'other']
    pivot_ccr8 = pivot_ccr8[[c for c in col_order if c in pivot_ccr8.columns]]
    pivot_mki67 = pivot_mki67[[c for c in col_order if c in pivot_mki67.columns]]
    
    sns.heatmap(pivot_ccr8, annot=True, fmt='.2f', cmap='RdYlGn', center=0, ax=axes[0], cbar_kws={'label': 'CCR8+ Z-score'})
    axes[0].set_title('CCR8+ Treg Signature by Sample')
    sns.heatmap(pivot_mki67, annot=True, fmt='.2f', cmap='RdYlGn', center=0, ax=axes[1], cbar_kws={'label': 'MKI67+ Z-score'})
    axes[1].set_title('MKI67+ Treg Signature by Sample')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'spatial_per_sample_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print("    Saved spatial_per_sample_heatmap.png")

print("\n" + "=" * 60)
print("DONE!")
print("=" * 60)
