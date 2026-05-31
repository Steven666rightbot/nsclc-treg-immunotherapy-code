"""Validate NSCLC annotation v2 - simplified working version."""
import os, gzip, numpy as np, pandas as pd, anndata
from scipy.io import mmread
from scipy import sparse
import scanpy as sc
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = r'D:\Research\tomato'
RAW_DIR = r'D:\Research\potato\data\space\GSE189487_RAW'
META_PATH = os.path.join(BASE_DIR, 'results', 'spatial_nsclc', 'nsclc_spot_metadata.csv')
FIG_DIR = os.path.join(BASE_DIR, 'figures', 'nsclc_validation')
os.makedirs(FIG_DIR, exist_ok=True)

MARKERS = {
    'CCR8+ Treg': ['CCR8', 'LAYN', 'BATF', 'CTLA4', 'TNFRSF18', 'FOXP3'],
    'MKI67+ Treg': ['MKI67', 'TOP2A', 'STMN1', 'TYMS', 'PCNA'],
    'Mechanosensing': ['ITGB1', 'CD44', 'YAP1', 'POSTN', 'COL1A2', 'ACTA2'],
}
SAMPLES = ['GSM5702473_TD1', 'GSM5702474_TD2', 'GSM5702475_TD3',
           'GSM5702476_TD5', 'GSM5702477_TD6', 'GSM5702478_TD8']

# 1. Load all samples
print("Loading NSCLC samples...")
adatas = []
for s in SAMPLES:
    print(f"  {s}...", end=" ", flush=True)
    with gzip.open(os.path.join(RAW_DIR, f'{s}_barcodes.tsv.gz'), 'rt') as f:
        bc = [line.strip() for line in f]
    with gzip.open(os.path.join(RAW_DIR, f'{s}_features.tsv.gz'), 'rt') as f:
        ft = [line.strip().split('\t')[1] for line in f]
    X = mmread(os.path.join(RAW_DIR, f'{s}_matrix.mtx.gz'))
    X = sparse.csr_matrix(X.T)
    a = anndata.AnnData(X=X)
    a.obs_names = bc
    a.var_names = ft
    a.var_names_make_unique()
    p = pd.read_csv(os.path.join(RAW_DIR, f'{s}_tissue_positions_list.csv.gz'), header=None)
    p.columns = ['barcode','in_tissue','array_row','array_col','pxl_col','pxl_row']
    p = p.set_index('barcode')
    a.obs = a.obs.join(p)
    a.obs['sample_id'] = s
    a = a[a.obs['in_tissue']==1].copy()
    print(f"{a.n_obs} spots")
    adatas.append(a)

adata = anndata.concat(adatas, label='sample_id', keys=SAMPLES,
                       index_unique='-', join='outer', merge='same')
print(f"\nMerged: {adata.n_obs} spots x {adata.n_vars} genes")

# 2. Merge metadata
meta = pd.read_csv(META_PATH, index_col=0)
common = list(set(adata.obs.index) & set(meta.index))
adata = adata[common].copy()
adata.obs['region'] = meta.loc[adata.obs.index, 'region']
adata.obs['ccr8_score'] = meta.loc[adata.obs.index, 'ccr8_score']
adata.obs['mki67_score'] = meta.loc[adata.obs.index, 'mki67_score']
print(f"Matched metadata: {len(common)} spots")

# 3. Preprocess
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 4. Check marker availability
print("\nMarker gene availability:")
available_markers = {}
for group, genes in MARKERS.items():
    avail = [g for g in genes if g in adata.var_names]
    available_markers[group] = avail
    print(f"  {group}: {avail}")

# 5. Statistical tests
print("\n" + "="*60)
print("Mann-Whitney U: hard_stroma vs tumor_core")
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
        log2fc = np.log2((tc_expr.mean() + 0.01) / (hs_expr.mean() + 0.01))
        results.append({'group': group, 'gene': gene,
                        'hs_mean': hs_expr.mean(), 'tc_mean': tc_expr.mean(),
                        'log2fc': log2fc, 'pvalue': pval})
        sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
        direction = 'TC>' if tc_expr.mean() > hs_expr.mean() else 'HS>'
        print(f"  {gene:12s} | HS={hs_expr.mean():.4f} TC={tc_expr.mean():.4f} "
              f"log2FC={log2fc:+.3f} p={pval:.2e} {sig} [{direction}]")

df_res = pd.DataFrame(results)
df_res.to_csv(os.path.join(FIG_DIR, 'marker_expression_stats.csv'), index=False)

# 6. Boxplots
print("\nGenerating boxplots...")
plot_data = []
for _, row in df_res.iterrows():
    gene = row['gene']
    gidx = adata.var_names.get_loc(gene)
    expr = adata.X[:, gidx].toarray().ravel()
    for region, mask in [('hard_stroma', hs_mask), ('tumor_core', tc_mask)]:
        for val in expr[mask]:
            plot_data.append({'gene': gene, 'region': region,
                              'group': row['group'], 'expression': val})
df_plot = pd.DataFrame(plot_data)

for group, genes in available_markers.items():
    if len(genes) == 0:
        continue
    df_g = df_plot[df_plot['group'] == group]
    fig, ax = plt.subplots(figsize=(max(6, len(genes)*1.2), 5))
    sns.boxplot(data=df_g, x='gene', y='expression', hue='region',
                palette={'hard_stroma': '#e74c3c', 'tumor_core': '#3498db'},
                ax=ax, showfliers=False)
    ax.set_title(f'{group} Marker Expression by Region', fontsize=13, fontweight='bold')
    ax.set_ylabel('log1p(norm counts)', fontsize=11)
    ax.legend(title='Region')
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

# 7. Score vs actual expression correlation
print("\nScore vs actual expression...")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for idx, (gene, score_col) in enumerate([('CCR8', 'ccr8_score'), ('MKI67', 'mki67_score')]):
    ax = axes[idx]
    if gene in adata.var_names:
        gidx = adata.var_names.get_loc(gene)
        expr = adata.X[:, gidx].toarray().ravel()
        ax.scatter(adata.obs[score_col], expr, alpha=0.05, s=1)
        ax.set_xlabel(f'{score_col} (signature)', fontsize=11)
        ax.set_ylabel(f'{gene} actual expression', fontsize=11)
        ax.set_title(f'{gene} Score vs Actual', fontsize=12, fontweight='bold')
        r, p = stats.spearmanr(adata.obs[score_col], expr)
        ax.text(0.05, 0.95, f'Spearman r={r:.3f}\np={p:.2e}', transform=ax.transAxes,
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
plt.tight_layout()
out = os.path.join(FIG_DIR, 'score_vs_actual_expression.png')
plt.savefig(out, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {out}")

print("\n" + "="*60)
print("VALIDATION COMPLETE")
print("="*60)
print(f"Results: {FIG_DIR}")
