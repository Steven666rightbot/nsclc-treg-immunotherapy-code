"""
NSCLC Spatial Validation (GSE189487)
Same-cancer validation using NSCLC Visium + NSCLC Treg signatures
"""

import os
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy.io import mmread
from scipy import sparse
from scipy import stats
from itertools import combinations
import pickle
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Paths
BASE_DIR = r'D:\Research\tomato'
SPATIAL_DIR = r'D:\Research\potato\data\space\GSE189487_RAW'
RESULTS_DIR = os.path.join(BASE_DIR, 'results', 'spatial_nsclc')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')
os.makedirs(RESULTS_DIR, exist_ok=True)

SAMPLES = ['GSM5702473_TD1', 'GSM5702474_TD2', 'GSM5702475_TD3',
           'GSM5702476_TD5', 'GSM5702477_TD6', 'GSM5702478_TD8']

# ==============================================================================
# STEP 1: Load all NSCLC Visium samples
# ==============================================================================
print("=" * 60)
print("STEP 1: Loading NSCLC Visium samples (GSE189487)")
print("=" * 60)

adatas = []
for sample in SAMPLES:
    print(f"  {sample}...", end=" ", flush=True)
    
    # Read barcodes
    with gzip.open(os.path.join(SPATIAL_DIR, f'{sample}_barcodes.tsv.gz'), 'rt') as f:
        barcodes = [line.strip() for line in f]
    
    # Read features
    with gzip.open(os.path.join(SPATIAL_DIR, f'{sample}_features.tsv.gz'), 'rt') as f:
        features = [line.strip().split('\t')[1] for line in f]
    
    # Read matrix
    X = mmread(os.path.join(SPATIAL_DIR, f'{sample}_matrix.mtx.gz'))
    X = sparse.csr_matrix(X.T)
    
    # Read positions
    pos = pd.read_csv(os.path.join(SPATIAL_DIR, f'{sample}_tissue_positions_list.csv.gz'), header=None)
    pos.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_col', 'pxl_row']
    pos = pos.set_index('barcode')
    
    # Create AnnData
    adata = anndata.AnnData(X=X)
    adata.obs_names = barcodes
    adata.var_names = features
    adata.var_names_make_unique()  # Fix duplicate gene names
    
    # Add spatial coords
    adata.obs = adata.obs.join(pos)
    adata.obs['sample_id'] = sample
    
    # Filter to tissue spots
    adata = adata[adata.obs['in_tissue'] == 1].copy()
    
    print(f"{adata.n_obs} spots")
    adatas.append(adata)

print(f"\nMerging {len(adatas)} samples...")
adata_spatial = anndata.concat(
    adatas, label='sample_id', keys=SAMPLES,
    index_unique='-', join='outer', merge='same'
)
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
print(f"  Stromal markers: {avail_stromal}")
print(f"  Tumor markers: {avail_tumor}")

sc.tl.score_genes(adata_spatial, gene_list=avail_stromal, score_name='stromal_score')
sc.tl.score_genes(adata_spatial, gene_list=avail_tumor, score_name='tumor_score')

adata_spatial.obs['region'] = 'other'
for sample in SAMPLES:
    mask = adata_spatial.obs['sample_id'] == sample
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
    print(f"    {sample}: {counts}")

# ==============================================================================
# STEP 4: Load Treg signatures
# ==============================================================================
print("\n" + "=" * 60)
print("STEP 4: Loading Treg signatures")
print("=" * 60)

sig_path = os.path.join(BASE_DIR, 'results', 'spatial', 'treg_signatures.pkl')
with open(sig_path, 'rb') as f:
    signatures = pickle.load(f)

for subtype in list(signatures.keys()):
    genes = signatures[subtype]
    available = [g for g in genes if g in adata_spatial.var_names]
    signatures[subtype] = available
    print(f"  {subtype}: {len(available)}/{len(genes)} markers in spatial data")

# Add canonical markers
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
        for sample in SAMPLES:
            mask = adata_spatial.obs['sample_id'] == sample
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

pd.DataFrame(results).to_csv(os.path.join(RESULTS_DIR, 'nsclc_spatial_stats.csv'), index=False)

# ==============================================================================
# STEP 7: Save
# ==============================================================================
print("\n" + "=" * 60)
print("STEP 7: Saving")
print("=" * 60)

adata_spatial.obs.to_csv(os.path.join(RESULTS_DIR, 'nsclc_spot_metadata.csv'))
print("  Metadata saved")

# ==============================================================================
# STEP 8: Plots
# ==============================================================================
print("\n" + "=" * 60)
print("STEP 8: Plots")
print("=" * 60)

# 8A. Scatter plot of spatial coordinates colored by region (for one sample)
print("  8A: Spatial scatter (TD1)...")
sample = 'GSM5702473_TD1'
mask = adata_spatial.obs['sample_id'] == sample
adata_rep = adata_spatial[mask].copy()

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, color, title in zip(axes, 
    ['stromal_score', 'tumor_score', 'region'],
    ['Stromal Score', 'Tumor Score', 'Region']):
    if color == 'region':
        palette = {'hard_stroma': '#E41A1C', 'tumor_core': '#377EB8', 'other': '#999999'}
        for region, color_val in palette.items():
            rmask = adata_rep.obs['region'] == region
            ax.scatter(adata_rep.obs.loc[rmask, 'pxl_col'], 
                      adata_rep.obs.loc[rmask, 'pxl_row'],
                      c=color_val, label=region, s=5, alpha=0.6)
        ax.legend()
    else:
        sc = ax.scatter(adata_rep.obs['pxl_col'], adata_rep.obs['pxl_row'],
                       c=adata_rep.obs[color], cmap='RdYlBu_r', s=5, alpha=0.6)
        plt.colorbar(sc, ax=ax)
    ax.set_title(f'{sample}: {title}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.invert_yaxis()
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'nsclc_spatial_region_definition.png'), dpi=300, bbox_inches='tight')
plt.close()
print("    Saved nsclc_spatial_region_definition.png")

# 8B. Treg scores scatter
print("  8B: Treg score scatter...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, color, title in zip(axes, ['ccr8_score', 'mki67_score'],
    ['CCR8+ Treg Signature', 'MKI67+ Treg Signature']):
    if color in adata_rep.obs.columns:
        sc = ax.scatter(adata_rep.obs['pxl_col'], adata_rep.obs['pxl_row'],
                       c=adata_rep.obs[color], cmap='YlOrRd', s=5, alpha=0.6)
        plt.colorbar(sc, ax=ax)
    ax.set_title(f'{sample}: {title}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.invert_yaxis()
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'nsclc_spatial_treg_scores.png'), dpi=300, bbox_inches='tight')
plt.close()
print("    Saved nsclc_spatial_treg_scores.png")

# 8C. Boxplot
print("  8C: Boxplot...")
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
    ax.set_title('NSCLC: Treg Subtype Signature by Region')
    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'nsclc_spatial_treg_boxplot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print("    Saved nsclc_spatial_treg_boxplot.png")

# 8D. Per-sample heatmap
print("  8D: Per-sample heatmap...")
sample_summary = []
for sample in SAMPLES:
    mask = adata_spatial.obs['sample_id'] == sample
    if mask.sum() == 0:
        continue
    for region in ['hard_stroma', 'tumor_core', 'other']:
        rmask = mask & (adata_spatial.obs['region'] == region)
        if rmask.sum() == 0:
            continue
        sample_summary.append({
            'sample': sample,
            'region': region,
            'ccr8_z': adata_spatial.obs.loc[rmask, 'ccr8_score_z'].mean(),
            'mki67_z': adata_spatial.obs.loc[rmask, 'mki67_score_z'].mean(),
            'n_spots': rmask.sum()
        })

df_summary = pd.DataFrame(sample_summary)
df_summary.to_csv(os.path.join(RESULTS_DIR, 'nsclc_per_sample_summary.csv'), index=False)

if len(df_summary) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    pivot_ccr8 = df_summary.pivot(index='sample', columns='region', values='ccr8_z')
    pivot_mki67 = df_summary.pivot(index='sample', columns='region', values='mki67_z')
    col_order = ['hard_stroma', 'tumor_core', 'other']
    pivot_ccr8 = pivot_ccr8[[c for c in col_order if c in pivot_ccr8.columns]]
    pivot_mki67 = pivot_mki67[[c for c in col_order if c in pivot_mki67.columns]]
    
    sns.heatmap(pivot_ccr8, annot=True, fmt='.2f', cmap='RdYlGn', center=0, ax=axes[0])
    axes[0].set_title('CCR8+ Treg Signature by Sample')
    sns.heatmap(pivot_mki67, annot=True, fmt='.2f', cmap='RdYlGn', center=0, ax=axes[1])
    axes[1].set_title('MKI67+ Treg Signature by Sample')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'nsclc_per_sample_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print("    Saved nsclc_per_sample_heatmap.png")

print("\n" + "=" * 60)
print("DONE!")
print("=" * 60)
