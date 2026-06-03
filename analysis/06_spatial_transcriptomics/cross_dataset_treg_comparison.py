"""
Cross-dataset comparison: NSCLC vs PDAC Treg spatial signals.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = r'D:\Research\tomato'
FIG_DIR = os.path.join(BASE_DIR, 'figures', 'cross_dataset_comparison')
OUT_DIR = os.path.join(BASE_DIR, 'results', 'cross_dataset')
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("Cross-Dataset Comparison: NSCLC vs PDAC")
print("=" * 60)

# ================================================================
# Load annotations
# ================================================================
nsclc = pd.read_csv(os.path.join(BASE_DIR, 'results', 'spatial_nsclc', 'treg_honest_spot_annotations.csv'), index_col=0)
pdac = pd.read_csv(os.path.join(BASE_DIR, 'results', 'spatial_pdac', 'treg_honest_spot_annotations.csv'), index_col=0)

nsclc['dataset'] = 'NSCLC'
pdac['dataset'] = 'PDAC'

# Add simplified region labels
region_map_nsclc = {'hard_stroma': 'Stroma', 'tumor_core': 'Tumor', 'other': 'Other'}
region_map_pdac = {'Epi': 'Tumor', 'Juxta': 'Juxta', 'Peri': 'Stroma'}
nsclc['region_simple'] = nsclc['region'].map(region_map_nsclc)
pdac['region_simple'] = pdac['region'].map(region_map_pdac)

# Combine
combined = pd.concat([nsclc[['dataset', 'region', 'region_simple', 'treg_score', 'prolif_score', 'act_score', 'treg_enriched', 'treg_class']],
                      pdac[['dataset', 'region', 'region_simple', 'treg_score', 'prolif_score', 'act_score', 'treg_enriched', 'treg_class']]],
                     ignore_index=True)

print(f"NSCLC: {len(nsclc)} spots")
print(f"PDAC: {len(pdac)} spots")

# ================================================================
# Summary statistics
# ================================================================
print("\n" + "=" * 60)
print("Summary Statistics")
print("=" * 60)

for metric in ['treg_score', 'prolif_score', 'act_score']:
    print(f"\n{metric}:")
    for ds in ['NSCLC', 'PDAC']:
        vals = combined[combined['dataset']==ds][metric].dropna()
        print(f"  {ds}: mean={vals.mean():.4f}, median={vals.median():.4f}, std={vals.std():.4f}, n={len(vals)}")
    # Mann-Whitney U test
    nsclc_vals = combined[combined['dataset']=='NSCLC'][metric].dropna()
    pdac_vals = combined[combined['dataset']=='PDAC'][metric].dropna()
    stat, pval = stats.mannwhitneyu(nsclc_vals, pdac_vals, alternative='two-sided')
    print(f"  Mann-Whitney p={pval:.2e}")

# ================================================================
# Regional comparisons
# ================================================================
print("\n" + "=" * 60)
print("Regional Comparisons")
print("=" * 60)

for ds in ['NSCLC', 'PDAC']:
    print(f"\n{ds}:")
    sub = combined[combined['dataset']==ds]
    for metric in ['treg_score', 'prolif_score', 'act_score']:
        print(f"  {metric} by region:")
        groups = []
        for region in sub['region_simple'].unique():
            vals = sub[sub['region_simple']==region][metric].dropna()
            if len(vals) > 0:
                groups.append(vals)
                print(f"    {region}: mean={vals.mean():.4f}, n={len(vals)}")
        if len(groups) >= 2:
            stat, pval = stats.kruskal(*groups)
            print(f"    Kruskal p={pval:.2e}")

# ================================================================
# Plotting
# ================================================================
print("\n" + "=" * 60)
print("Plotting")
print("=" * 60)

# 1A. Overall score distributions
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
metrics = ['treg_score', 'prolif_score', 'act_score']
titles = ['Pan-Treg Score', 'Proliferation Signal', 'Activation Signal']
for ax, metric, title in zip(axes, metrics, titles):
    sns.boxplot(data=combined, x='dataset', y=metric, ax=ax, showfliers=False)
    ax.set_title(title, fontweight='bold')
    ax.set_ylabel('Score')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'cross_dataset_score_comparison.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: cross_dataset_score_comparison.png")

# 1B. By region
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, metric, title in zip(axes, metrics, titles):
    # Only Stroma vs Tumor for fair comparison
    sub = combined[combined['region_simple'].isin(['Stroma', 'Tumor'])]
    sns.boxplot(data=sub, x='dataset', y=metric, hue='region_simple', ax=ax, showfliers=False)
    ax.set_title(title, fontweight='bold')
    ax.set_ylabel('Score')
    ax.legend(title='Region')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'cross_dataset_region_comparison.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: cross_dataset_region_comparison.png")

# 1C. Treg spot classification composition
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, ds in zip(axes, ['NSCLC', 'PDAC']):
    sub = combined[combined['dataset']==ds]
    ct = pd.crosstab(sub['region_simple'], sub['treg_class'], normalize='index') * 100
    cols = [c for c in ['prolif_dominant', 'act_dominant', 'unclassified_Treg'] if c in ct.columns]
    ct = ct[cols]
    ct.plot(kind='bar', stacked=True, ax=ax, color=['#e74c3c', '#3498db', '#95a5a6'])
    ax.set_title(f'{ds}', fontweight='bold')
    ax.set_ylabel('% of Spots')
    ax.set_xlabel('Region')
    ax.legend(title='Signal', labels=['Proliferation', 'Activation', 'Unclassified'])
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'cross_dataset_treg_composition.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: cross_dataset_treg_composition.png")

# 1D. Treg enrichment rate by region
fig, ax = plt.subplots(figsize=(8, 6))
rate_data = []
for ds in ['NSCLC', 'PDAC']:
    sub = combined[combined['dataset']==ds]
    for region in sub['region_simple'].unique():
        rsub = sub[sub['region_simple']==region]
        rate = rsub['treg_enriched'].mean() * 100
        rate_data.append({'dataset': ds, 'region': region, 'rate': rate})
rate_df = pd.DataFrame(rate_data)
rate_pivot = rate_df.pivot(index='region', columns='dataset', values='rate')
rate_pivot.plot(kind='bar', ax=ax, color=['#e67e22', '#9b59b6'])
ax.set_title('Treg Enrichment Rate by Region', fontweight='bold')
ax.set_ylabel('% Treg-enriched spots')
ax.set_xlabel('Region')
ax.legend(title='Dataset')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'cross_dataset_treg_enrichment_rate.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: cross_dataset_treg_enrichment_rate.png")

# ================================================================
# Save summary
# ================================================================
print("\n" + "=" * 60)
print("Saving results")
print("=" * 60)

# Summary table
summary_rows = []
for ds in ['NSCLC', 'PDAC']:
    sub = combined[combined['dataset']==ds]
    summary_rows.append({
        'dataset': ds,
        'total_spots': len(sub),
        'treg_enriched': sub['treg_enriched'].sum(),
        'treg_enriched_pct': sub['treg_enriched'].mean() * 100,
        'prolif_dominant': (sub['treg_class']=='prolif_dominant').sum(),
        'act_dominant': (sub['treg_class']=='act_dominant').sum(),
        'treg_score_mean': sub['treg_score'].mean(),
        'prolif_score_mean': sub['prolif_score'].mean(),
        'act_score_mean': sub['act_score'].mean(),
    })
summary = pd.DataFrame(summary_rows)
summary.to_csv(os.path.join(OUT_DIR, 'cross_dataset_summary.csv'), index=False)
print(summary.to_string(index=False))

# Detailed regional summary
regional_summary = []
for ds in ['NSCLC', 'PDAC']:
    sub = combined[combined['dataset']==ds]
    for region in sub['region_simple'].unique():
        rsub = sub[sub['region_simple']==region]
        regional_summary.append({
            'dataset': ds,
            'region': region,
            'n_spots': len(rsub),
            'treg_enriched_pct': rsub['treg_enriched'].mean() * 100,
            'treg_score_mean': rsub['treg_score'].mean(),
            'prolif_score_mean': rsub['prolif_score'].mean(),
            'act_score_mean': rsub['act_score'].mean(),
            'prolif_dominant_pct': (rsub['treg_class']=='prolif_dominant').mean() * 100,
            'act_dominant_pct': (rsub['treg_class']=='act_dominant').mean() * 100,
        })
regional_df = pd.DataFrame(regional_summary)
regional_df.to_csv(os.path.join(OUT_DIR, 'cross_dataset_regional_summary.csv'), index=False)

print(f"\nAll results saved to: {OUT_DIR}")
print(f"All figures saved to: {FIG_DIR}")
