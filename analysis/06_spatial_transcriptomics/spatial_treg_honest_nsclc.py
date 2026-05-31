"""
Honest Treg spatial analysis for NSCLC Visium data.
Does NOT claim CCR8+ or MKI67+ Treg subtypes.
Reports: pan-Treg enrichment, proliferation signal, activation signal.
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
FIG_DIR = os.path.join(BASE_DIR, 'figures', 'nsclc_treg_honest')
OUT_DIR = os.path.join(BASE_DIR, 'results', 'spatial_nsclc')
os.makedirs(FIG_DIR, exist_ok=True)

SAMPLES = ['GSM5702473_TD1', 'GSM5702474_TD2', 'GSM5702475_TD3',
           'GSM5702476_TD5', 'GSM5702477_TD6', 'GSM5702478_TD8']

# ================================================================
# Marker definitions (evidence-based, not label-name-based)
# ================================================================
PAN_TREG_MARKERS = ['FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF18', 'TIGIT', 'LAG3', 'IKZF2']
PROLIF_MARKERS = ['CCNB1', 'AURKA', 'PCNA', 'TOP2A', 'STMN1', 'UBE2C']
ACTIVATION_MARKERS = ['TNFRSF18', 'HLA-DRB1', 'CX3CR1', 'CD28', 'TNFRSF8']

print("=" * 60)
print("NSCLC: Honest Treg Spatial Analysis")
print("=" * 60)

# ================================================================
# Step 1: Load NSCLC spatial data
# ================================================================
print("\nSTEP 1: Loading NSCLC spatial data")
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
# Step 2: Score pan-Treg enrichment
# ================================================================
print("\nSTEP 2: Scoring pan-Treg enrichment")
avail_pan = [g for g in PAN_TREG_MARKERS if g in adata.var_names]
print(f"Pan-Treg markers available: {len(avail_pan)}/{len(PAN_TREG_MARKERS)}")
print(avail_pan)

sc.tl.score_genes(adata, gene_list=avail_pan, score_name='treg_score')
print(f"Treg score range: {adata.obs['treg_score'].min():.3f} to {adata.obs['treg_score'].max():.3f}")

# Define Treg-enriched spots (top 20% per sample)
adata.obs['treg_enriched'] = False
for sample in SAMPLES:
    mask = adata.obs['sample_id'] == sample
    if mask.sum() == 0: continue
    thresh = adata.obs.loc[mask, 'treg_score'].quantile(0.8)
    adata.obs.loc[mask & (adata.obs['treg_score'] >= thresh), 'treg_enriched'] = True

n_treg = adata.obs['treg_enriched'].sum()
print(f"Treg-enriched spots: {n_treg} / {adata.n_obs} ({100*n_treg/adata.n_obs:.1f}%)")

# ================================================================
# Step 3: Score proliferation and activation signals
# ================================================================
print("\nSTEP 3: Scoring biological signals in Treg spots")

avail_prolif = [g for g in PROLIF_MARKERS if g in adata.var_names]
avail_act = [g for g in ACTIVATION_MARKERS if g in adata.var_names]
print(f"Proliferation markers available: {len(avail_prolif)}/{len(PROLIF_MARKERS)}")
print(avail_prolif)
print(f"Activation markers available: {len(avail_act)}/{len(ACTIVATION_MARKERS)}")
print(avail_act)

if len(avail_prolif) >= 3:
    sc.tl.score_genes(adata, gene_list=avail_prolif, score_name='prolif_score')
else:
    adata.obs['prolif_score'] = 0.0

if len(avail_act) >= 3:
    sc.tl.score_genes(adata, gene_list=avail_act, score_name='act_score')
else:
    adata.obs['act_score'] = 0.0

# ================================================================
# Step 4: Classify Treg spots by signal dominance
# ================================================================
print("\nSTEP 4: Classifying Treg-enriched spots by dominant signal")
treg_mask = adata.obs['treg_enriched']
adata.obs['treg_class'] = 'non_Treg'

# Z-score within Treg spots
for score_col in ['prolif_score', 'act_score']:
    vals = adata.obs.loc[treg_mask, score_col]
    adata.obs.loc[treg_mask, f'{score_col}_z'] = (vals - vals.mean()) / vals.std() if vals.std() > 0 else 0

for idx in treg_mask[treg_mask].index:
    p_z = adata.obs.loc[idx, 'prolif_score_z']
    a_z = adata.obs.loc[idx, 'act_score_z']
    if p_z > 0.5 and p_z > a_z:
        adata.obs.loc[idx, 'treg_class'] = 'prolif_dominant'
    elif a_z > 0.5 and a_z > p_z:
        adata.obs.loc[idx, 'treg_class'] = 'act_dominant'
    else:
        adata.obs.loc[idx, 'treg_class'] = 'unclassified_Treg'

print("\nTreg spot classification:")
print(adata.obs['treg_class'].value_counts())

# ================================================================
# Step 5: Regional validation
# ================================================================
print("\n" + "=" * 60)
print("STEP 5: Regional distribution")
print("=" * 60)

crosstab = pd.crosstab(adata.obs['region'], adata.obs['treg_class'], normalize='index') * 100
print("\n% within each region:")
print(crosstab.round(2).to_string())

# Statistical tests: hard_stroma vs tumor_core
hs_mask = adata.obs['region'] == 'hard_stroma'
tc_mask = adata.obs['region'] == 'tumor_core'

for score_col in ['treg_score', 'prolif_score', 'act_score']:
    hs_vals = adata.obs.loc[hs_mask, score_col].dropna()
    tc_vals = adata.obs.loc[tc_mask, score_col].dropna()
    stat, pval = stats.mannwhitneyu(hs_vals, tc_vals, alternative='two-sided')
    print(f"\n{score_col}:")
    print(f"  hard_stroma: mean={hs_vals.mean():.4f}, n={len(hs_vals)}")
    print(f"  tumor_core:  mean={tc_vals.mean():.4f}, n={len(tc_vals)}")
    print(f"  p={pval:.2e}")

# ================================================================
# Step 6: Plotting
# ================================================================
print("\n" + "=" * 60)
print("STEP 6: Plotting")
print("=" * 60)

# 6A. Pan-Treg score by region
fig, ax = plt.subplots(figsize=(7, 5))
sns.boxplot(data=adata.obs, x='region', y='treg_score', ax=ax, showfliers=False,
            order=['hard_stroma', 'tumor_core', 'other'])
ax.set_title('Pan-Treg Score by Region (NSCLC)', fontweight='bold')
ax.set_ylabel('Pan-Treg Score')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'nsclc_treg_score_by_region.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: nsclc_treg_score_by_region.png")

# 6B. Proliferation vs Activation in Treg spots
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
sns.boxplot(data=adata.obs[adata.obs['treg_enriched']], x='region', y='prolif_score',
            ax=ax, showfliers=False, order=['hard_stroma', 'tumor_core', 'other'])
ax.set_title('Proliferation Signal (Treg spots)', fontweight='bold')
ax.set_ylabel('Proliferation Score')

ax = axes[1]
sns.boxplot(data=adata.obs[adata.obs['treg_enriched']], x='region', y='act_score',
            ax=ax, showfliers=False, order=['hard_stroma', 'tumor_core', 'other'])
ax.set_title('Activation Signal (Treg spots)', fontweight='bold')
ax.set_ylabel('Activation Score')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'nsclc_treg_signals_by_region.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: nsclc_treg_signals_by_region.png")

# 6C. Treg spot composition by region
fig, ax = plt.subplots(figsize=(8, 5))
crosstab_counts = pd.crosstab(adata.obs['region'], adata.obs['treg_class'])
crosstab_counts = crosstab_counts[['prolif_dominant', 'act_dominant', 'unclassified_Treg']]
crosstab_counts.plot(kind='bar', stacked=True, ax=ax,
                     color=['#e74c3c', '#3498db', '#95a5a6'])
ax.set_title('Treg Spot Composition by Region (NSCLC)', fontweight='bold')
ax.set_ylabel('Number of Spots')
ax.set_xlabel('Region')
ax.legend(title='Signal Type', labels=['Proliferation-dominant', 'Activation-dominant', 'Unclassified'])
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'nsclc_treg_composition_by_region.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: nsclc_treg_composition_by_region.png")

# 6D. Scatter: proliferation vs activation in Treg spots
fig, ax = plt.subplots(figsize=(7, 7))
treg_spots = adata.obs[adata.obs['treg_enriched']]
colors = {'prolif_dominant': '#e74c3c', 'act_dominant': '#3498db', 'unclassified_Treg': '#95a5a6'}
for subtype in ['prolif_dominant', 'act_dominant', 'unclassified_Treg']:
    sub = treg_spots[treg_spots['treg_class'] == subtype]
    ax.scatter(sub['prolif_score'], sub['act_score'],
               c=colors[subtype], label=subtype, alpha=0.5, s=10)
ax.set_xlabel('Proliferation Score', fontsize=11)
ax.set_ylabel('Activation Score', fontsize=11)
ax.set_title('Treg Spots: Proliferation vs Activation (NSCLC)', fontweight='bold')
ax.legend()
ax.axhline(y=treg_spots['act_score'].mean(), color='gray', linestyle='--', alpha=0.3)
ax.axvline(x=treg_spots['prolif_score'].mean(), color='gray', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'nsclc_treg_signal_scatter.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: nsclc_treg_signal_scatter.png")

# ================================================================
# Save results
# ================================================================
print("\n" + "=" * 60)
print("Saving results")
print("=" * 60)

# Save marker lists
marker_info = pd.DataFrame({
    'pan_treg': pd.Series(PAN_TREG_MARKERS),
    'proliferation': pd.Series(PROLIF_MARKERS),
    'activation': pd.Series(ACTIVATION_MARKERS),
})
marker_info.to_csv(os.path.join(OUT_DIR, 'treg_honest_markers.csv'), index=False)

# Save spot annotations
adata.obs.to_csv(os.path.join(OUT_DIR, 'treg_honest_spot_annotations.csv'))

# Save summary
summary = pd.DataFrame({
    'metric': ['total_spots', 'treg_enriched_spots', 'prolif_dominant_spots', 'act_dominant_spots',
               'prolif_markers_found', 'act_markers_found', 'pan_treg_markers_found'],
    'value': [adata.n_obs, n_treg,
              (adata.obs['treg_class']=='prolif_dominant').sum(),
              (adata.obs['treg_class']=='act_dominant').sum(),
              len(avail_prolif), len(avail_act), len(avail_pan)]
})
summary.to_csv(os.path.join(OUT_DIR, 'treg_honest_summary.csv'), index=False)

print(f"\nAll results saved to: {OUT_DIR}")
print(f"All figures saved to: {FIG_DIR}")
print("\nNOTE: This analysis reports Treg-related BIOLOGICAL SIGNALS,")
print("      NOT specific Treg subtypes. CCR8+ and MKI67+ labels are NOT used.")
