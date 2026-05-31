"""
Validate NSCLC cell annotation: check if Treg subtype marker genes
are differentially expressed in hard_stroma vs tumor_core regions.
"""

import os
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy.io import mmread
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ================================================================
# Config
# ================================================================
BASE_DIR = r'D:\Research\tomato'
RAW_DIR = r'D:\Research\potato\data\space\GSE189487_RAW'
META_PATH = os.path.join(BASE_DIR, 'results', 'spatial_nsclc', 'nsclc_spot_metadata.csv')
FIG_DIR = os.path.join(BASE_DIR, 'figures', 'nsclc_validation')
os.makedirs(FIG_DIR, exist_ok=True)

# Marker genes to validate
MARKERS = {
    'CCR8+ Treg': ['CCR8', 'LAYN', 'BATF', 'CTLA4', 'TNFRSF18', 'FOXP3'],
    'MKI67+ Treg': ['MKI67', 'TOP2A', 'STMN1', 'TYMS', 'PCNA'],
    'Mechanosensing': ['ITGB1', 'CD44', 'YAP1', 'POSTN', 'COL1A2', 'ACTA2'],
}

SAMPLES = ['GSM5702473_TD1', 'GSM5702474_TD2', 'GSM5702475_TD3',
           'GSM5702476_TD5', 'GSM5702477_TD6', 'GSM5702478_TD8']

# ================================================================
# 1. Load all NSCLC samples
# ================================================================
print("Loading NSCLC samples...")
adatas = []
for sample in SAMPLES:
    print(f"  {sample}...", end=" ", flush=True)
    with gzip.open(os.path.join(RAW_DIR, f'{sample}_barcodes.tsv.gz'), 'rt') as f:
        barcodes = [line.strip() for line in f]
    with gzip.open(os.path.join(RAW_DIR, f'{sample}_features.tsv.gz'), 'rt') as f:
        features = [line.strip().split('\t')[1] for line in f]
    X = mmread(os.path.join(RAW_DIR, f'{sample}_matrix.mtx.gz'))
    X = sparse.csr_matrix(X.T)
    ad = anndata.AnnData(X=X)
    ad.obs_names = [f"{b}-{sample}" for b in barcodes]
    ad.var_names = features
    ad.var_names_make_unique()
    ad.obs['sample_id'] = sample
    adatas.append(ad)
    print(f"{ad.n_obs} spots")

adata = anndata.concat(adatas, label='sample_id', keys=SAMPLES,
                       index_unique='-', join='outer', merge='same')
print(f"\nMerged: {adata.n_obs} spots x {adata.n_vars} genes")

# ================================================================
# 2. Merge region metadata
# ================================================================
print("\nMerging region metadata...")
meta = pd.read_csv(META_PATH, index_col=0)
# NSCLC metadata index format: spot_id (e.g., AAACAAGTATCTCCCA-1-GSM5702473_TD1)
common = list(set(adata.obs.index) & set(meta.index))
adata = adata[common].copy()
adata.obs['region'] = meta.loc[adata.obs.index, 'region']
adata.obs['ccr8_score'] = meta.loc[adata.obs.index, 'ccr8_score']
adata.obs['mki67_score'] = meta.loc[adata.obs.index, 'mki67_score']
print(f"  Matched {len(common)} spots with metadata")

# ================================================================
# 3. Preprocess
# ================================================================
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ================================================================
# 4. Check marker gene availability
# ================================================================
print("\nMarker gene availability:")
available_markers = {}
for group, genes in MARKERS.items():
    avail = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    available_markers[group] = avail
    print(f"  {group}: {avail} (missing: {missing})")

# ================================================================
# 5. Statistical test: hard_stroma vs tumor_core
# ================================================================
print("\n" + "="*60)
print("Statistical tests (Mann-Whitney U): hard_stroma vs tumor_core")
print("="*60)

results = []
hs_mask = adata.obs['region'] == 'hard_stroma'
tc_mask = adata.obs['region'] == 'tumor_core'

for group, genes in available_markers.items():
    for gene in genes:
        gidx = adata.var_names.get_loc(gene)
        expr = adata.X[:, gidx].toarray().ravel()
        hs_expr = expr[hs_mask]
        tc_expr = expr[tc_mask]
        stat, pval = stats.mannwhitneyu(hs_expr, tc_expr, alternative='two-sided')
        
        results.append({
            'group': group, 'gene': gene,
            'hs_mean': hs_expr.mean(), 'tc_mean': tc_expr.mean(),
            'hs_median': np.median(hs_expr), 'tc_median': np.median(tc_expr),
            'log2fc': np.log2((tc_expr.mean() + 0.01) / (hs_expr.mean() + 0.01)),
            'pvalue': pval,
        })
        
        sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
        direction = 'TC>' if tc_expr.mean() > hs_expr.mean() else 'HS>'
        print(f"  {gene:12s} | HS={hs_expr.mean():.4f}  TC={tc_expr.mean():.4f}  "
              f"log2FC={results[-1]['log2fc']:+.3f}  p={pval:.2e} {sig} [{direction}]")

df_res = pd.DataFrame(results)
df_res.to_csv(os.path.join(FIG_DIR, 'marker_expression_stats.csv'), index=False)

# ================================================================
# 6. Boxplot comparison
# ================================================================
print("\nGenerating boxplots...")

# Flatten for seaborn
plot_data = []
for _, row in df_res.iterrows():
    gene = row['gene']
    gidx = adata.var_names.get_loc(gene)
    expr = adata.X[:, gidx].toarray().ravel()
    for region, mask in [('hard_stroma', hs_mask), ('tumor_core', tc_mask)]:
        for val in expr[mask]:
            plot_data.append({
                'gene': gene, 'region': region,
                'group': row['group'], 'expression': val
            })
df_plot = pd.DataFrame(plot_data)

# One figure per marker group
for group, genes in available_markers.items():
    if len(genes) == 0:
        continue
    df_g = df_plot[df_plot['group'] == group]
    
    fig, ax = plt.subplots(figsize=(max(6, len(genes)*1.2), 5))
    sns.boxplot(data=df_g, x='gene', y='expression', hue='region',
                palette={'hard_stroma': '#e74c3c', 'tumor_core': '#3498db'},
                ax=ax, showfliers=False)
    ax.set_title(f'{group} Marker Expression by Region', fontsize=13, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('log1p(norm counts)', fontsize=11)
    ax.legend(title='Region')
    
    # Add significance stars
    for i, gene in enumerate(genes):
        row = df_res[(df_res['group']==group) & (df_res['gene']==gene)].iloc[0]
        p = row['pvalue']
        stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        y_max = df_g[df_g['gene']==gene]['expression'].max()
        ax.text(i, y_max*1.05, stars, ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    out = os.path.join(FIG_DIR, f'{group.replace("+", "plus").replace(" ", "_")}_boxplot.png')
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")

# ================================================================
# 7. Single-sample spatial expression map (TD1)
# ================================================================
print("\nGenerating spatial expression maps (TD1)...")

td1_mask = adata.obs['sample_id'] == 'GSM5702473_TD1'
adata_td1 = adata[td1_mask].copy()

# Load coords for TD1
pos = pd.read_csv(os.path.join(RAW_DIR, 'GSM5702473_TD1_tissue_positions_list.csv.gz'), header=None)
pos.columns = ['barcode','in_tissue','array_row','array_col','pxl_col','pxl_row']
pos = pos.set_index('barcode')
# Match barcodes
bc_map = {b.split('-')[0] if '-' in b else b: b for b in adata_td1.obs.index}
pos_td1 = pos.loc[[b.split('-')[0] if '-' in b else b for b in adata_td1.obs.index]]
coords = pos_td1[['pxl_col', 'pxl_row']].values

# Plot key genes
key_genes = ['CCR8', 'MKI67', 'ITGB1', 'COL1A2']
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
axes = axes.flatten()

for idx, gene in enumerate(key_genes):
    ax = axes[idx]
    if gene not in adata_td1.var_names:
        ax.set_title(f'{gene} NOT FOUND')
        ax.axis('off')
        continue
    
    gidx = adata_td1.var_names.get_loc(gene)
    expr = adata_td1.X[:, gidx].toarray().ravel()
    
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=expr, cmap='YlOrRd',
                    s=15, alpha=0.8, edgecolors='none')
    ax.set_title(f'{gene}', fontsize=14, fontweight='bold')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.axis('off')
    plt.colorbar(sc, ax=ax, shrink=0.6, label='Expression')

plt.suptitle('NSCLC TD1: Spatial Expression of Key Markers', fontsize=16, fontweight='bold')
plt.tight_layout()
out = os.path.join(FIG_DIR, 'TD1_spatial_expression_map.png')
plt.savefig(out, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {out}")

# ================================================================
# 8. Scatter: CCR8 score vs actual CCR8 expression
# ================================================================
print("\nValidating CCR8/MKI67 scores against actual gene expression...")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# CCR8
ax = axes[0]
gidx = adata.var_names.get_loc('CCR8') if 'CCR8' in adata.var_names else None
if gidx is not None:
    ccr8_expr = adata.X[:, gidx].toarray().ravel()
    ax.scatter(adata.obs['ccr8_score'], ccr8_expr, alpha=0.1, s=1)
    ax.set_xlabel('CCR8 Score (signature-based)', fontsize=11)
    ax.set_ylabel('CCR8 Actual Expression', fontsize=11)
    ax.set_title('CCR8 Score vs Actual CCR8', fontsize=12, fontweight='bold')
    # Correlation
    from scipy.stats import spearmanr
    r, p = spearmanr(adata.obs['ccr8_score'], ccr8_expr)
    ax.text(0.05, 0.95, f'Spearman r={r:.3f}\np={p:.2e}', transform=ax.transAxes,
            fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# MKI67
ax = axes[1]
gidx = adata.var_names.get_loc('MKI67') if 'MKI67' in adata.var_names else None
if gidx is not None:
    mki67_expr = adata.X[:, gidx].toarray().ravel()
    ax.scatter(adata.obs['mki67_score'], mki67_expr, alpha=0.1, s=1)
    ax.set_xlabel('MKI67 Score (signature-based)', fontsize=11)
    ax.set_ylabel('MKI67 Actual Expression', fontsize=11)
    ax.set_title('MKI67 Score vs Actual MKI67', fontsize=12, fontweight='bold')
    r, p = spearmanr(adata.obs['mki67_score'], mki67_expr)
    ax.text(0.05, 0.95, f'Spearman r={r:.3f}\np={p:.2e}', transform=ax.transAxes,
            fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
out = os.path.join(FIG_DIR, 'score_vs_actual_expression.png')
plt.savefig(out, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {out}")

# ================================================================
# Done
# ================================================================
print("\n" + "="*60)
print("VALIDATION COMPLETE")
print("="*60)
print(f"Results saved to: {FIG_DIR}")
print(f"Stats: {os.path.join(FIG_DIR, 'marker_expression_stats.csv')}")
