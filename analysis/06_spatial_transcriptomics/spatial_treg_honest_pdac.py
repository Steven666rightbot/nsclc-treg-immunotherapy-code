"""
Honest Treg spatial analysis for PDAC Visium data.
Does NOT claim CCR8+ or MKI67+ Treg subtypes.
Reports: pan-Treg enrichment, proliferation signal, activation signal.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = r'D:\Research\tomato'
DATA_DIR = os.path.join(BASE_DIR, 'data', 'spatial', 'pdac_visium', 'extracted')
META_PATH = os.path.join(BASE_DIR, 'results', 'spatial_pdac', 'pdac_spot_metadata.csv')
FIG_DIR = os.path.join(BASE_DIR, 'figures', 'pdac_treg_honest')
OUT_DIR = os.path.join(BASE_DIR, 'results', 'spatial_pdac')
os.makedirs(FIG_DIR, exist_ok=True)

SAMPLES = ['HG_1', 'HG_2', 'HG_3',
           'LG_1', 'LG_2', 'LG_3', 'LG_4', 'LG_5', 'LG_6', 'LG_7',
           'PDAC_1', 'PDAC_2', 'PDAC_3']

# ================================================================
# Marker definitions
# ================================================================
PAN_TREG_MARKERS = ['FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF18', 'TIGIT', 'LAG3', 'IKZF2']
PROLIF_MARKERS = ['CCNB1', 'AURKA', 'PCNA', 'TOP2A', 'STMN1', 'UBE2C']
ACTIVATION_MARKERS = ['TNFRSF18', 'HLA-DRB1', 'CX3CR1', 'CD28', 'TNFRSF8']

print("=" * 60)
print("PDAC: Honest Treg Spatial Analysis")
print("=" * 60)

# ================================================================
# Step 1: Load PDAC spatial data
# ================================================================
print("\nSTEP 1: Loading PDAC spatial data")
adatas = []
for sample in SAMPLES:
    sample_path = os.path.join(DATA_DIR, sample)
    if not os.path.exists(sample_path):
        print(f"  WARNING: {sample} not found, skipping")
        continue
    ad = sc.read_visium(sample_path, count_file='filtered_feature_bc_matrix.h5', load_images=False)
    ad.var_names_make_unique()
    # filtered_feature_bc_matrix already contains only in-tissue spots
    ad.obs['sample_id'] = sample
    # Add sample suffix for unique barcodes across samples
    ad.obs['barcode_raw'] = ad.obs_names
    ad.obs_names = [f"{b}_{sample}" for b in ad.obs_names]
    adatas.append(ad)
    print(f"  {sample}: {ad.n_obs} spots x {ad.n_vars} genes")

adata = anndata.concat(adatas, label='sample_id', keys=[ad.obs['sample_id'].iloc[0] for ad in adatas],
                       index_unique=None, join='outer', merge='same')
print(f"Merged: {adata.n_obs} spots x {adata.n_vars} genes")

# Preprocess
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Merge metadata
meta = pd.read_csv(META_PATH, index_col=0)
print(f"\nMetadata barcodes example: {meta.index[0]}")
print(f"AnnData barcodes example: {adata.obs_names[0]} (raw={adata.obs['barcode_raw'].iloc[0]})")

# Metadata barcodes have format like ACTATTTCCGGGCCCA-1_5 (suffix is sample index)
# H5 barcodes have format like ACTATTTCCGGGCCCA-1
# Need to strip the trailing _N from metadata barcodes
meta['barcode_clean'] = meta.index.str.replace(r'_\d+$', '', regex=True)

# For each sample, match barcodes using barcode_raw
matched_barcodes = []
for sample in SAMPLES:
    if sample not in meta['Code'].values:
        continue
    meta_sample = meta[meta['Code'] == sample]
    adata_sample = adata[adata.obs['sample_id'] == sample]
    
    meta_bcs = set(meta_sample['barcode_clean'])
    adata_bcs = set(adata_sample.obs['barcode_raw'])
    overlap = list(meta_bcs & adata_bcs)
    
    # Map adata barcodes (with suffix) to metadata index
    meta_map = dict(zip(meta_sample['barcode_clean'], meta_sample.index))
    adata_raw_map = dict(zip(adata_sample.obs['barcode_raw'], adata_sample.obs_names))
    for bc in overlap:
        matched_barcodes.append((adata_raw_map[bc], sample, meta_map[bc]))
    print(f"  {sample}: matched={len(overlap)}, meta={len(meta_sample)}, adata={len(adata_sample)}")

if len(matched_barcodes) > 0:
    match_df = pd.DataFrame(matched_barcodes, columns=['adata_name', 'sample_id', 'meta_barcode'])
    adata = adata[match_df['adata_name'].tolist()].copy()
    match_df = match_df.set_index('adata_name')
    adata.obs['region'] = meta.loc[match_df.loc[adata.obs_names, 'meta_barcode'].values, 'Region'].values
    adata.obs['type'] = meta.loc[match_df.loc[adata.obs_names, 'meta_barcode'].values, 'Type'].values
    print(f"Matched metadata: {len(matched_barcodes)} spots")
else:
    print("WARNING: Could not match metadata, proceeding without region labels")
    adata.obs['region'] = 'unknown'
    adata.obs['type'] = 'unknown'

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
# Step 5: Regional and type validation
# ================================================================
print("\n" + "=" * 60)
print("STEP 5: Regional distribution")
print("=" * 60)

crosstab = pd.crosstab(adata.obs['region'], adata.obs['treg_class'], normalize='index') * 100
print("\n% within each region:")
print(crosstab.round(2).to_string())

# By Type (LG vs HG vs PDAC)
type_tab = pd.crosstab(adata.obs['type'], adata.obs['treg_class'], normalize='index') * 100
print("\n% within each type:")
print(type_tab.round(2).to_string())

# Statistical tests: Epi vs Peri
epi_mask = adata.obs['region'] == 'Epi'
peri_mask = adata.obs['region'] == 'Peri'

for score_col in ['treg_score', 'prolif_score', 'act_score']:
    epi_vals = adata.obs.loc[epi_mask, score_col].dropna()
    peri_vals = adata.obs.loc[peri_mask, score_col].dropna()
    if len(epi_vals) > 0 and len(peri_vals) > 0:
        stat, pval = stats.mannwhitneyu(epi_vals, peri_vals, alternative='two-sided')
        print(f"\n{score_col} (Epi vs Peri):")
        print(f"  Epi:  mean={epi_vals.mean():.4f}, n={len(epi_vals)}")
        print(f"  Peri: mean={peri_vals.mean():.4f}, n={len(peri_vals)}")
        print(f"  p={pval:.2e}")

# LG vs HG vs PDAC
for score_col in ['treg_score', 'prolif_score', 'act_score']:
    groups = []
    for t in ['LG', 'HG', 'IPMNPDAC']:
        vals = adata.obs.loc[adata.obs['type'] == t, score_col].dropna()
        if len(vals) > 0:
            groups.append(vals)
    if len(groups) >= 2:
        stat, pval = stats.kruskal(*groups)
        print(f"\n{score_col} (LG vs HG vs PDAC): Kruskal p={pval:.2e}")
        for t in ['LG', 'HG', 'IPMNPDAC']:
            vals = adata.obs.loc[adata.obs['type'] == t, score_col].dropna()
            if len(vals) > 0:
                print(f"  {t}: mean={vals.mean():.4f}, n={len(vals)}")

# ================================================================
# Step 6: Plotting
# ================================================================
print("\n" + "=" * 60)
print("STEP 6: Plotting")
print("=" * 60)

# 6A. Pan-Treg score by region
fig, ax = plt.subplots(figsize=(7, 5))
sns.boxplot(data=adata.obs, x='region', y='treg_score', ax=ax, showfliers=False,
            order=['Epi', 'Juxta', 'Peri'])
ax.set_title('Pan-Treg Score by Region (PDAC)', fontweight='bold')
ax.set_ylabel('Pan-Treg Score')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'pdac_treg_score_by_region.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: pdac_treg_score_by_region.png")

# 6B. By type
fig, ax = plt.subplots(figsize=(7, 5))
sns.boxplot(data=adata.obs, x='type', y='treg_score', ax=ax, showfliers=False,
            order=['LG', 'HG', 'IPMNPDAC'])
ax.set_title('Pan-Treg Score by Lesion Type (PDAC)', fontweight='bold')
ax.set_ylabel('Pan-Treg Score')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'pdac_treg_score_by_type.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: pdac_treg_score_by_type.png")

# 6C. Proliferation vs Activation in Treg spots
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
sns.boxplot(data=adata.obs[adata.obs['treg_enriched']], x='region', y='prolif_score',
            ax=ax, showfliers=False, order=['Epi', 'Juxta', 'Peri'])
ax.set_title('Proliferation Signal (Treg spots)', fontweight='bold')
ax.set_ylabel('Proliferation Score')

ax = axes[1]
sns.boxplot(data=adata.obs[adata.obs['treg_enriched']], x='region', y='act_score',
            ax=ax, showfliers=False, order=['Epi', 'Juxta', 'Peri'])
ax.set_title('Activation Signal (Treg spots)', fontweight='bold')
ax.set_ylabel('Activation Score')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'pdac_treg_signals_by_region.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: pdac_treg_signals_by_region.png")

# 6D. Composition by region
fig, ax = plt.subplots(figsize=(8, 5))
crosstab_counts = pd.crosstab(adata.obs['region'], adata.obs['treg_class'])
cols = [c for c in ['prolif_dominant', 'act_dominant', 'unclassified_Treg'] if c in crosstab_counts.columns]
crosstab_counts = crosstab_counts[cols]
crosstab_counts.plot(kind='bar', stacked=True, ax=ax, color=['#e74c3c', '#3498db', '#95a5a6'])
ax.set_title('Treg Spot Composition by Region (PDAC)', fontweight='bold')
ax.set_ylabel('Number of Spots')
ax.set_xlabel('Region')
ax.legend(title='Signal Type')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'pdac_treg_composition_by_region.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: pdac_treg_composition_by_region.png")

# 6E. Composition by type
fig, ax = plt.subplots(figsize=(8, 5))
crosstab_counts = pd.crosstab(adata.obs['type'], adata.obs['treg_class'])
cols = [c for c in ['prolif_dominant', 'act_dominant', 'unclassified_Treg'] if c in crosstab_counts.columns]
crosstab_counts = crosstab_counts[cols]
crosstab_counts.plot(kind='bar', stacked=True, ax=ax, color=['#e74c3c', '#3498db', '#95a5a6'])
ax.set_title('Treg Spot Composition by Lesion Type (PDAC)', fontweight='bold')
ax.set_ylabel('Number of Spots')
ax.set_xlabel('Type')
ax.legend(title='Signal Type')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'pdac_treg_composition_by_type.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: pdac_treg_composition_by_type.png")

# 6F. Scatter
fig, ax = plt.subplots(figsize=(7, 7))
treg_spots = adata.obs[adata.obs['treg_enriched']]
colors = {'prolif_dominant': '#e74c3c', 'act_dominant': '#3498db', 'unclassified_Treg': '#95a5a6'}
for subtype in ['prolif_dominant', 'act_dominant', 'unclassified_Treg']:
    sub = treg_spots[treg_spots['treg_class'] == subtype]
    if len(sub) > 0:
        ax.scatter(sub['prolif_score'], sub['act_score'],
                   c=colors[subtype], label=subtype, alpha=0.5, s=10)
ax.set_xlabel('Proliferation Score', fontsize=11)
ax.set_ylabel('Activation Score', fontsize=11)
ax.set_title('Treg Spots: Proliferation vs Activation (PDAC)', fontweight='bold')
ax.legend()
ax.axhline(y=treg_spots['act_score'].mean(), color='gray', linestyle='--', alpha=0.3)
ax.axvline(x=treg_spots['prolif_score'].mean(), color='gray', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'pdac_treg_signal_scatter.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: pdac_treg_signal_scatter.png")

# ================================================================
# Save results
# ================================================================
print("\n" + "=" * 60)
print("Saving results")
print("=" * 60)

marker_info = pd.DataFrame({
    'pan_treg': pd.Series(PAN_TREG_MARKERS),
    'proliferation': pd.Series(PROLIF_MARKERS),
    'activation': pd.Series(ACTIVATION_MARKERS),
})
marker_info.to_csv(os.path.join(OUT_DIR, 'treg_honest_markers.csv'), index=False)

adata.obs.to_csv(os.path.join(OUT_DIR, 'treg_honest_spot_annotations.csv'))

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
