"""
Two-step Treg subtype annotation on NSCLC spatial data:
Step 1: Identify Treg-enriched spots using pan-Treg markers
Step 2: Within Treg spots, score CCR8+ vs MKI67+ subtypes
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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = r'D:\Research\tomato'
RAW_DIR = r'D:\Research\potato\data\space\GSE189487_RAW'
META_PATH = os.path.join(BASE_DIR, 'results', 'spatial_nsclc', 'nsclc_spot_metadata.csv')
FIG_DIR = os.path.join(BASE_DIR, 'figures', 'nsclc_treg_2step')
OUT_DIR = os.path.join(BASE_DIR, 'results', 'spatial_nsclc')
os.makedirs(FIG_DIR, exist_ok=True)

SAMPLES = ['GSM5702473_TD1', 'GSM5702474_TD2', 'GSM5702475_TD3',
           'GSM5702476_TD5', 'GSM5702477_TD6', 'GSM5702478_TD8']

# ================================================================
# Step 0: Find Treg subtype markers from scRNA-seq
# ================================================================
print("=" * 60)
print("STEP 0: Finding Treg subtype markers from scRNA-seq")
print("=" * 60)

# Load Treg h5ad (already normalized)
treg = sc.read(os.path.join(BASE_DIR, 'data', 'monocle3_input', 'treg_processed.h5ad'))
print(f"Treg data: {treg.n_obs} cells x {treg.n_vars} genes")

# Use scanpy's rank_genes_groups for fast DE
# CCR8+ vs FOXP3+
print("  Running DE: CCR8+ vs FOXP3+...")
treg_ccr8_foxp3 = treg[treg.obs['sub_cell_type'].isin(['CD4T_Treg_CCR8', 'CD4T_Treg_FOXP3'])].copy()
sc.tl.rank_genes_groups(treg_ccr8_foxp3, groupby='sub_cell_type', groups=['CD4T_Treg_CCR8'], reference='CD4T_Treg_FOXP3', method='t-test')
ccr8_de = sc.get.rank_genes_groups_df(treg_ccr8_foxp3, group='CD4T_Treg_CCR8')
ccr8_markers = ccr8_de[(ccr8_de['logfoldchanges'] > 1.0) & (ccr8_de['pvals_adj'] < 0.01)]['names'].tolist()[:50]
print(f"\nCCR8+ markers (vs FOXP3+): {len(ccr8_markers)} genes")
print(ccr8_markers[:20])

# MKI67+ vs FOXP3+
print("  Running DE: MKI67+ vs FOXP3+...")
treg_mki_foxp3 = treg[treg.obs['sub_cell_type'].isin(['CD4T_Treg_MKI67', 'CD4T_Treg_FOXP3'])].copy()
sc.tl.rank_genes_groups(treg_mki_foxp3, groupby='sub_cell_type', groups=['CD4T_Treg_MKI67'], reference='CD4T_Treg_FOXP3', method='t-test')
mki_de = sc.get.rank_genes_groups_df(treg_mki_foxp3, group='CD4T_Treg_MKI67')
mki_markers = mki_de[(mki_de['logfoldchanges'] > 1.0) & (mki_de['pvals_adj'] < 0.01)]['names'].tolist()[:50]
print(f"\nMKI67+ markers (vs FOXP3+): {len(mki_markers)} genes")
print(mki_markers[:20])

# Pan-Treg markers (genes high in all Treg vs general expectation)
# Use known canonical markers + top genes from CCR8+ and MKI67+ that are also high in FOXP3+
pan_treg_candidates = ['FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF18', 'BATF', 'IKZF2', 'TIGIT', 'LAG3']
# Add genes that are highly expressed in FOXP3+ as well
foxp3_mean = np.asarray(treg[treg.obs['sub_cell_type']=='CD4T_Treg_FOXP3'].X.mean(axis=0)).ravel()
top_foxp3 = treg.var_names[np.argsort(foxp3_mean)[-30:]].tolist()
pan_treg_markers = list(set(pan_treg_candidates + [g for g in top_foxp3 if g in treg.var_names]))[:25]
pan_treg_markers = [g for g in pan_treg_markers if g in treg.var_names]
print(f"\nPan-Treg markers: {len(pan_treg_markers)} genes")
print(pan_treg_markers)

# ================================================================
# Step 1: Load NSCLC spatial and score Treg enrichment
# ================================================================
print("\n" + "=" * 60)
print("STEP 1: Loading NSCLC spatial data")
print("=" * 60)

adatas = []
for sample in SAMPLES:
    with gzip.open(os.path.join(RAW_DIR, f'{sample}_barcodes.tsv.gz'), 'rt') as f:
        barcodes = [line.strip() for line in f]
    with gzip.open(os.path.join(RAW_DIR, f'{sample}_features.tsv.gz'), 'rt') as f:
        features = [line.strip().split('\t')[1] for line in f]
    X = mmread(os.path.join(RAW_DIR, f'{sample}_matrix.mtx.gz'))
    X = sparse.csr_matrix(X.T)
    ad = anndata.AnnData(X=X)
    ad.obs_names = barcodes
    ad.var_names = features
    ad.var_names_make_unique()
    p = pd.read_csv(os.path.join(RAW_DIR, f'{sample}_tissue_positions_list.csv.gz'), header=None)
    p.columns = ['barcode','in_tissue','array_row','array_col','pxl_col','pxl_row']
    p = p.set_index('barcode')
    ad.obs = ad.obs.join(p)
    ad.obs['sample_id'] = sample
    ad = ad[ad.obs['in_tissue']==1].copy()
    adatas.append(ad)

adata = anndata.concat(adatas, label='sample_id', keys=SAMPLES,
                       index_unique='-', join='outer', merge='same')
print(f"Merged: {adata.n_obs} spots x {adata.n_vars} genes")

# Preprocess
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Merge region metadata
meta = pd.read_csv(META_PATH, index_col=0)
common = list(set(adata.obs.index) & set(meta.index))
adata = adata[common].copy()
adata.obs['region'] = meta.loc[adata.obs.index, 'region']
print(f"Matched metadata: {len(common)} spots")

# ================================================================
# Step 1b: Score Treg enrichment
# ================================================================
print("\n" + "=" * 60)
print("STEP 1b: Scoring Treg enrichment")
print("=" * 60)

avail_pan = [g for g in pan_treg_markers if g in adata.var_names]
print(f"Pan-Treg markers in spatial: {len(avail_pan)}/{len(pan_treg_markers)}")
print(avail_pan)

sc.tl.score_genes(adata, gene_list=avail_pan, score_name='treg_score')
print(f"Treg score range: {adata.obs['treg_score'].min():.3f} to {adata.obs['treg_score'].max():.3f}")

# Define Treg-enriched spots (top 20% per sample)
adata.obs['treg_enriched'] = False
for sample in SAMPLES:
    mask = adata.obs['sample_id'] == sample
    if mask.sum() == 0:
        continue
    vals = adata.obs.loc[mask, 'treg_score']
    thresh = vals.quantile(0.8)
    adata.obs.loc[mask & (adata.obs['treg_score'] >= thresh), 'treg_enriched'] = True

n_treg = adata.obs['treg_enriched'].sum()
print(f"Treg-enriched spots: {n_treg} / {adata.n_obs} ({100*n_treg/adata.n_obs:.1f}%)")

# ================================================================
# Step 2: Within Treg spots, score CCR8+ vs MKI67+
# ================================================================
print("\n" + "=" * 60)
print("STEP 2: Scoring CCR8+ vs MKI67+ in Treg spots")
print("=" * 60)

# Only score markers available in spatial
avail_ccr8 = [g for g in ccr8_markers if g in adata.var_names]
avail_mki = [g for g in mki_markers if g in adata.var_names]

print(f"CCR8+ markers in spatial: {len(avail_ccr8)}/{len(ccr8_markers)}")
print(f"MKI67+ markers in spatial: {len(avail_mki)}/{len(mki_markers)}")

# Score all spots (will only use Treg-enriched ones for classification)
if len(avail_ccr8) >= 3:
    sc.tl.score_genes(adata, gene_list=avail_ccr8, score_name='ccr8_score_new')
else:
    print("WARNING: Too few CCR8 markers, using single gene if available")
    if 'CCR8' in adata.var_names:
        adata.obs['ccr8_score_new'] = adata[:, 'CCR8'].X.toarray().ravel()
    else:
        adata.obs['ccr8_score_new'] = 0

if len(avail_mki) >= 3:
    sc.tl.score_genes(adata, gene_list=avail_mki, score_name='mki67_score_new')
else:
    print("WARNING: Too few MKI67 markers")
    adata.obs['mki67_score_new'] = 0

# Classify Treg-enriched spots
treg_mask = adata.obs['treg_enriched']
adata.obs['treg_subtype'] = 'non_Treg'

# Within Treg spots, classify based on relative scores
ccr8_vals = adata.obs.loc[treg_mask, 'ccr8_score_new']
mki_vals = adata.obs.loc[treg_mask, 'mki67_score_new']

# Z-score within Treg spots
ccr8_z = (ccr8_vals - ccr8_vals.mean()) / ccr8_vals.std() if ccr8_vals.std() > 0 else pd.Series(0, index=ccr8_vals.index)
mki_z = (mki_vals - mki_vals.mean()) / mki_vals.std() if mki_vals.std() > 0 else pd.Series(0, index=mki_vals.index)

for idx in treg_mask[treg_mask].index:
    c = ccr8_z.loc[idx]
    m = mki_z.loc[idx]
    if c > 0.5 and c > m:
        adata.obs.loc[idx, 'treg_subtype'] = 'CCR8+_Treg'
    elif m > 0.5 and m > c:
        adata.obs.loc[idx, 'treg_subtype'] = 'MKI67+_Treg'
    else:
        adata.obs.loc[idx, 'treg_subtype'] = 'unclassified_Treg'

print("\nTreg subtype distribution:")
print(adata.obs['treg_subtype'].value_counts())

# ================================================================
# Step 3: Validate by region
# ================================================================
print("\n" + "=" * 60)
print("STEP 3: Regional distribution validation")
print("=" * 60)

crosstab = pd.crosstab(adata.obs['region'], adata.obs['treg_subtype'], normalize='index') * 100
print("\n% within each region:")
print(crosstab.round(2).to_string())

# Statistical test: CCR8+ vs MKI67+ in hard_stroma vs tumor_core
hs_mask = adata.obs['region'] == 'hard_stroma'
tc_mask = adata.obs['region'] == 'tumor_core'

for score_col in ['ccr8_score_new', 'mki67_score_new']:
    print(f"\n{score_col}:")
    hs_vals = adata.obs.loc[hs_mask, score_col].dropna()
    tc_vals = adata.obs.loc[tc_mask, score_col].dropna()
    stat, pval = stats.mannwhitneyu(hs_vals, tc_vals, alternative='two-sided')
    print(f"  hard_stroma: mean={hs_vals.mean():.4f}, n={len(hs_vals)}")
    print(f"  tumor_core:  mean={tc_vals.mean():.4f}, n={len(tc_vals)}")
    print(f"  p={pval:.2e}")

# ================================================================
# Step 4: Plotting
# ================================================================
print("\n" + "=" * 60)
print("STEP 4: Plotting")
print("=" * 60)

# 4A. Treg score distribution by region
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
sns.boxplot(data=adata.obs, x='region', y='treg_score', ax=ax, showfliers=False,
            order=['hard_stroma', 'tumor_core', 'other'])
ax.set_title('Pan-Treg Score by Region', fontweight='bold')
ax.set_ylabel('Treg Score')

ax = axes[1]
sns.boxplot(data=adata.obs[adata.obs['treg_enriched']], x='region', y='ccr8_score_new', ax=ax, showfliers=False,
            order=['hard_stroma', 'tumor_core', 'other'])
ax.set_title('CCR8+ Score (Treg spots only)', fontweight='bold')
ax.set_ylabel('CCR8+ Score')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'treg_scores_by_region.png'), dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: treg_scores_by_region.png")

# 4B. Subtype composition by region
fig, ax = plt.subplots(figsize=(8, 5))
crosstab_counts = pd.crosstab(adata.obs['region'], adata.obs['treg_subtype'])
crosstab_counts = crosstab_counts[['CCR8+_Treg', 'MKI67+_Treg', 'unclassified_Treg']]
crosstab_counts.plot(kind='bar', stacked=True, ax=ax, 
                     color=['#e74c3c', '#3498db', '#95a5a6'])
ax.set_title('Treg Subtype Composition by Region', fontweight='bold')
ax.set_ylabel('Number of Spots')
ax.set_xlabel('Region')
ax.legend(title='Treg Subtype')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'treg_subtype_by_region.png'), dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: treg_subtype_by_region.png")

# 4C. Scatter: CCR8+ vs MKI67+ scores in Treg spots
fig, ax = plt.subplots(figsize=(7, 7))
treg_spots = adata.obs[adata.obs['treg_enriched']]
colors = {'CCR8+_Treg': '#e74c3c', 'MKI67+_Treg': '#3498db', 'unclassified_Treg': '#95a5a6'}
for subtype in ['CCR8+_Treg', 'MKI67+_Treg', 'unclassified_Treg']:
    sub = treg_spots[treg_spots['treg_subtype'] == subtype]
    ax.scatter(sub['ccr8_score_new'], sub['mki67_score_new'], 
               c=colors[subtype], label=subtype, alpha=0.5, s=10)
ax.set_xlabel('CCR8+ Score', fontsize=11)
ax.set_ylabel('MKI67+ Score', fontsize=11)
ax.set_title('Treg Spots: CCR8+ vs MKI67+', fontweight='bold')
ax.legend()
ax.axhline(y=treg_spots['mki67_score_new'].mean(), color='gray', linestyle='--', alpha=0.3)
ax.axvline(x=treg_spots['ccr8_score_new'].mean(), color='gray', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'treg_subtype_scatter.png'), dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: treg_subtype_scatter.png")

# ================================================================
# Save results
# ================================================================
print("\n" + "=" * 60)
print("Saving results")
print("=" * 60)

# Save marker lists
marker_info = pd.DataFrame({
    'pan_treg': pd.Series(pan_treg_markers),
    'ccr8_specific': pd.Series(ccr8_markers),
    'mki67_specific': pd.Series(mki_markers),
})
marker_info.to_csv(os.path.join(OUT_DIR, 'treg_2step_markers.csv'), index=False)

# Save spot annotations
adata.obs.to_csv(os.path.join(OUT_DIR, 'treg_2step_spot_annotations.csv'))

# Save summary
summary = pd.DataFrame({
    'metric': ['total_spots', 'treg_enriched_spots', 'ccr8_treg_spots', 'mki67_treg_spots',
               'ccr8_markers_found', 'mki67_markers_found', 'pan_treg_markers'],
    'value': [adata.n_obs, n_treg, 
              (adata.obs['treg_subtype']=='CCR8+_Treg').sum(),
              (adata.obs['treg_subtype']=='MKI67+_Treg').sum(),
              len(avail_ccr8), len(avail_mki), len(avail_pan)]
})
summary.to_csv(os.path.join(OUT_DIR, 'treg_2step_summary.csv'), index=False)

print(f"\nAll results saved to: {OUT_DIR}")
print(f"All figures saved to: {FIG_DIR}")
