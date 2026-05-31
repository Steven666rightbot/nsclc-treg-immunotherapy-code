#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE243013 - Treg细胞拟时序分析 (Pseudotime Analysis) v2
使用pandas分块读取22GB mtx，只保留Treg细胞
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import sparse
from scipy.stats import spearmanr, mannwhitneyu

warnings.filterwarnings('ignore')
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, facecolor='white', fontsize=10)

# ===================== 路径配置 =====================
DATA_DIR = r"D:\Research\土豆\数据\raw\GSE243013"
OUT_DIR = r"D:\Research\tomato\figures\pseudotime_treg"
RESULT_DIR = r"D:\Research\tomato\data\pseudotime_treg"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(RESULT_DIR, exist_ok=True)

# 进度日志
LOG_FILE = r"D:\Research\tomato\data\pseudotime_treg\progress.log"
class TeeLogger:
    def __init__(self, filepath):
        self.terminal = sys.stdout
        self.log = open(filepath, 'w', encoding='utf-8')
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()
    def flush(self):
        self.terminal.flush()
        self.log.flush()
sys.stdout = TeeLogger(LOG_FILE)
sys.stderr = sys.stdout

# ===================== 1. 读取Treg metadata =====================
print("=" * 60)
print("Step 1: 读取Treg metadata")
print("=" * 60)

treg_meta = pd.read_csv(
    os.path.join(DATA_DIR, "GSE243013_Treg_cells_metadata.csv"),
    low_memory=False
)
print(f"Treg cells: {len(treg_meta)}")
print(f"Subtypes: {treg_meta['sub_cell_type'].value_counts().to_dict()}")

# 过滤unknown响应
treg_meta = treg_meta[treg_meta['pathological_response'] != 'unknowm'].copy()
print(f"过滤unknown后: {len(treg_meta)}")

# 定义响应标签
treg_meta['response'] = treg_meta['pathological_response'].map({
    'MPR': 'Responder', 'pCR': 'Responder', 'non-MPR': 'Non-responder'
})
treg_meta['response_binary'] = treg_meta['pathological_response'].map({
    'MPR': 1, 'pCR': 1, 'non-MPR': 0
})

# ===================== 2. 读取表达矩阵并筛选Treg =====================
print("\n" + "=" * 60)
print("Step 2: 读取表达矩阵 (22GB mtx) - 分块读取")
print("=" * 60)

# 读取barcodes
barcodes = pd.read_csv(
    os.path.join(DATA_DIR, "GSE243013_barcodes.csv.gz"),
    compression='gzip', header=None, names=['barcode']
)
barcodes_list = barcodes['barcode'].tolist()
print(f"Total barcodes: {len(barcodes_list)}")

# 读取features (genes)
features = pd.read_csv(
    os.path.join(DATA_DIR, "features.csv.gz"),
    compression='gzip', header=None, names=['geneSymbol']
)
gene_names = features['geneSymbol'].tolist()
print(f"Total genes: {len(gene_names)}")

# 读取预存的Treg行掩码
is_treg = np.load(r"D:\Research\tomato\data\treg_row_mask.npy")
print(f"Treg rows: {is_treg.sum()}")

# 分块读取mtx
mtx_path = os.path.join(DATA_DIR, "GSE243013_NSCLC_immune_scRNA_counts.mtx")
print(f"Reading mtx: {mtx_path}")
print("This will take ~15-20 minutes...")

start_time = time.time()

chunk_rows = []
chunk_cols = []
chunk_data = []
chunk_count = 0
total_kept = 0

for chunk in pd.read_csv(
    mtx_path, skiprows=2, header=None, sep=r'\s+',
    dtype={0: 'int32', 1: 'int32', 2: 'float32'},
    chunksize=100_000_000, engine='c'
):
    chunk_start = time.time()
    mask = is_treg[chunk[0].values]
    kept = mask.sum()
    if kept > 0:
        chunk_rows.append(chunk[0].values[mask] - 1)
        chunk_cols.append(chunk[1].values[mask] - 1)
        chunk_data.append(chunk[2].values[mask])
        total_kept += kept
    chunk_count += 1
    elapsed = time.time() - start_time
    progress = chunk_count / 20  # approx 20 chunks total
    eta = elapsed / progress - elapsed if progress > 0 else 0
    print(f"  Chunk {chunk_count:02d}/~20: kept {kept:>8,}/{len(chunk):>10,} "
          f"({kept/len(chunk)*100:5.2f}%) | elapsed: {elapsed/60:5.1f}min | ETA: {eta/60:5.1f}min")
    sys.stdout.flush()

print(f"\nTotal chunks: {chunk_count}")
print(f"Total kept entries: {total_kept}")
print(f"Read time: {(time.time()-start_time)/60:.1f} minutes")
sys.stdout.flush()

# 合并数据
print("\nMerging chunks...")
all_rows = np.concatenate(chunk_rows)
all_cols = np.concatenate(chunk_cols)
all_data = np.concatenate(chunk_data)

# 构建稀疏矩阵 (genes x cells)
# mtx是 cells x genes 的坐标格式，但读取后我需要转置为 cells x genes
# 原始格式: row=cell_idx, col=gene_idx
X = sparse.coo_matrix(
    (all_data, (all_rows, all_cols)),
    shape=(len(barcodes_list) - 1, len(gene_names))  # exclude header row
)
print(f"Sparse matrix shape (cells x genes): {X.shape}")
print("Converting to CSR format for fast row slicing...")
X = X.tocsr()
print("CSR conversion done")

# 筛选Treg细胞
treg_indices = np.where(is_treg[1:])[0]  # exclude header row (index 0)
print(f"Treg cell indices: {len(treg_indices)}")

# 注意：is_treg[1:] 排除了表头行，treg_indices对应稀疏矩阵的行号
# 但treg_indices是相对于原始mtx的0-based行号
# 我们需要用这些行号来子集化X

# 筛选并重新排序，与metadata中的cellID一致
treg_cell_ids = treg_meta['cellID'].tolist()
barcodes_array = np.array(barcodes_list[1:])  # exclude header

# 创建行号到cellID的映射
row_to_cellid = {i: barcodes_array[i] for i in treg_indices}

# 快速构建barcode到索引的映射（O(1)查找）
print("Building barcode index map...")
sys.stdout.flush()
barcode_to_idx = {bc: i for i, bc in enumerate(barcodes_array)}

ordered_indices = []
ordered_cell_ids = []
for cid in treg_cell_ids:
    if cid in barcode_to_idx:
        ordered_indices.append(barcode_to_idx[cid])
        ordered_cell_ids.append(cid)

ordered_indices = np.array(ordered_indices)
print(f"Ordered Treg cells: {len(ordered_indices)}")
sys.stdout.flush()

# 子集化稀疏矩阵
X_treg = X[ordered_indices, :]
print(f"Treg matrix shape: {X_treg.shape}")

# 构建AnnData
print("Building AnnData...")
sys.stdout.flush()
adata_treg = sc.AnnData(X_treg.tocsr())
adata_treg.obs_names = ordered_cell_ids
adata_treg.var_names = gene_names
print(f"AnnData built: {adata_treg.shape}")

# 合并metadata
print("Merging metadata...")
sys.stdout.flush()
obs_df = treg_meta.set_index('cellID')
adata_treg.obs = adata_treg.obs.join(obs_df)
print(f"Final Treg adata: {adata_treg.shape}")
sys.stdout.flush()

# 降采样：每亚型最多5000个细胞，加速DPT计算
print("Subsampling for DPT speed...")
sys.stdout.flush()
n_per_subtype = 5000
sampled_indices = []
for subtype in ['CD4T_Treg_FOXP3', 'CD4T_Treg_CCR8', 'CD4T_Treg_MKI67']:
    subtype_idx = np.where(adata_treg.obs['sub_cell_type'] == subtype)[0]
    if len(subtype_idx) > n_per_subtype:
        np.random.seed(42)
        sampled = np.random.choice(subtype_idx, n_per_subtype, replace=False)
    else:
        sampled = subtype_idx
    sampled_indices.extend(sampled)
    print(f"  {subtype}: {len(sampled)} / {len(subtype_idx)}")

adata_treg = adata_treg[sampled_indices].copy()
print(f"After subsampling: {adata_treg.shape}")
sys.stdout.flush()

# 释放内存
del X, X_treg, chunk_rows, chunk_cols, chunk_data, all_rows, all_cols, all_data

# ===================== 3. 预处理 =====================
print("\n" + "=" * 60)
print("Step 3: 预处理")
print("=" * 60)

sc.pp.filter_genes(adata_treg, min_cells=10)
print(f"After filtering genes (min_cells=10): {adata_treg.shape}")

sc.pp.normalize_total(adata_treg, target_sum=1e4)
sc.pp.log1p(adata_treg)

sc.pp.highly_variable_genes(adata_treg, n_top_genes=2000, subset=True)
print(f"After HVG selection: {adata_treg.shape}")

print("Scaling (zero_center=False to keep sparse)...")
sys.stdout.flush()
sc.pp.scale(adata_treg, max_value=10, zero_center=False)
print("Scale done")
sys.stdout.flush()

# ===================== 4. 降维与聚类 =====================
print("\n" + "=" * 60)
print("Step 4: PCA + 邻居图 + UMAP + 聚类")
print("=" * 60)

print("Running PCA...")
sys.stdout.flush()
sc.tl.pca(adata_treg, svd_solver='arpack')
print("PCA done")

print("Running neighbors...")
sys.stdout.flush()
sc.pp.neighbors(adata_treg, n_neighbors=15, n_pcs=30)
print("Neighbors done")

print("Running UMAP...")
sys.stdout.flush()
sc.tl.umap(adata_treg)
print("UMAP done")

print("Running leiden clustering...")
sys.stdout.flush()
sc.tl.leiden(adata_treg, resolution=0.5)
print("Leiden done")

print(f"Leiden clusters: {adata_treg.obs['leiden'].nunique()}")
sys.stdout.flush()

# ===================== 5. PAGA + 轨迹推断 =====================
print("\n" + "=" * 60)
print("Step 5: PAGA + Diffusion Map + DPT")
print("=" * 60)

print("Running PAGA...")
sys.stdout.flush()
sc.tl.paga(adata_treg, groups='sub_cell_type')
print("PAGA done")

print("Running diffmap...")
sys.stdout.flush()
sc.tl.diffmap(adata_treg, n_comps=15)
print("Diffmap done")
sys.stdout.flush()

# 以FOXP3+ Treg作为根节点
foxp3_idx = np.where(adata_treg.obs['sub_cell_type'] == 'CD4T_Treg_FOXP3')[0]
if len(foxp3_idx) > 0:
    root_idx = foxp3_idx[np.argmin(adata_treg.obsm['X_diffmap'][foxp3_idx, 1])]
    adata_treg.uns['iroot'] = root_idx
    print(f"Root set to FOXP3+ Treg cell index {root_idx}")

print("Running DPT...")
sys.stdout.flush()
sc.tl.dpt(adata_treg, n_dcs=10)
print(f"DPT done. Range: {adata_treg.obs['dpt_pseudotime'].min():.3f} - {adata_treg.obs['dpt_pseudotime'].max():.3f}")
sys.stdout.flush()

# ===================== 6. 可视化 =====================
print("\n" + "=" * 60)
print("Step 6: 可视化")
print("=" * 60)

# 6.1 UMAP: Treg亚型
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sc.pl.umap(adata_treg, color='sub_cell_type', ax=axes[0], show=False, title='Treg Subtypes',
           legend_loc='right margin', palette={'CD4T_Treg_FOXP3': '#1f77b4',
                                               'CD4T_Treg_CCR8': '#ff7f0e',
                                               'CD4T_Treg_MKI67': '#2ca02c'})
axes[0].set_xlabel('UMAP1')
axes[0].set_ylabel('UMAP2')

sc.pl.umap(adata_treg, color='response', ax=axes[1], show=False, title='Response',
           legend_loc='right margin', palette={'Responder': '#d62728', 'Non-responder': '#9467bd'})
axes[1].set_xlabel('UMAP1')
axes[1].set_ylabel('UMAP2')

sc.pl.umap(adata_treg, color='dpt_pseudotime', ax=axes[2], show=False, title='Pseudotime',
           color_map='viridis')
axes[2].set_xlabel('UMAP1')
axes[2].set_ylabel('UMAP2')

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '01_Treg_UMAP_overview.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 01_Treg_UMAP_overview.png")

# 6.2 PAGA图
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.paga(adata_treg, color='sub_cell_type', ax=ax, show=False,
           title='PAGA: Treg Trajectory',
           node_size_scale=3, edge_width_scale=1)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '02_Treg_PAGA.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 02_Treg_PAGA.png")

# 6.3 UMAP + PAGA
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.paga_compare(adata_treg, color='sub_cell_type', ax=ax, show=False,
                   title='Treg Trajectory (PAGA + UMAP)')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '03_Treg_PAGA_UMAP.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 03_Treg_PAGA_UMAP.png")

# ===================== 7. 拟时序分析 =====================
print("\n" + "=" * 60)
print("Step 7: 拟时序分析")
print("=" * 60)

# 7.1 各亚型在拟时序上的分布
fig, ax = plt.subplots(figsize=(10, 6))
for subtype, color in [('CD4T_Treg_FOXP3', '#1f77b4'),
                        ('CD4T_Treg_CCR8', '#ff7f0e'),
                        ('CD4T_Treg_MKI67', '#2ca02c')]:
    subset = adata_treg.obs[adata_treg.obs['sub_cell_type'] == subtype]['dpt_pseudotime']
    sns.kdeplot(subset, label=subtype, color=color, ax=ax, fill=True, alpha=0.3)
ax.set_xlabel('Pseudotime')
ax.set_ylabel('Density')
ax.set_title('Treg Subtype Distribution along Pseudotime')
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '04_Subtype_pseudotime_density.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 04_Subtype_pseudotime_density.png")

# 7.2 响应状态在拟时序上的分布
fig, ax = plt.subplots(figsize=(10, 6))
for resp, color in [('Responder', '#d62728'), ('Non-responder', '#9467bd')]:
    subset = adata_treg.obs[adata_treg.obs['response'] == resp]['dpt_pseudotime']
    sns.kdeplot(subset, label=resp, color=color, ax=ax, fill=True, alpha=0.3)
ax.set_xlabel('Pseudotime')
ax.set_ylabel('Density')
ax.set_title('Response Distribution along Pseudotime')
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '05_Response_pseudotime_density.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 05_Response_pseudotime_density.png")

# 7.3 箱线图
fig, ax = plt.subplots(figsize=(10, 6))
df_plot = adata_treg.obs[['sub_cell_type', 'response', 'dpt_pseudotime']].copy()
sns.boxplot(data=df_plot, x='sub_cell_type', y='dpt_pseudotime', hue='response',
            palette={'Responder': '#d62728', 'Non-responder': '#9467bd'}, ax=ax)
ax.set_xlabel('Treg Subtype')
ax.set_ylabel('Pseudotime')
ax.set_title('Pseudotime by Subtype and Response')
plt.xticks(rotation=15)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '06_Pseudotime_boxplot.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 06_Pseudotime_boxplot.png")

# ===================== 8. 拟时序相关基因 =====================
print("\n" + "=" * 60)
print("Step 8: 拟时序相关基因")
print("=" * 60)

adata_treg.raw = adata_treg

X_raw = adata_treg.raw.X
if sparse.issparse(X_raw):
    X_raw = X_raw.toarray()

genes = adata_treg.raw.var_names.tolist()
pseudotime = adata_treg.obs['dpt_pseudotime'].values

print(f"Computing correlations for {len(genes)} genes...")
corrs = []
pvals = []
for i, gene in enumerate(genes):
    expr = X_raw[:, i]
    mask = expr > 0
    if mask.sum() < 50:
        corrs.append(np.nan)
        pvals.append(np.nan)
        continue
    corr, pval = spearmanr(expr[mask], pseudotime[mask])
    corrs.append(corr)
    pvals.append(pval)

gene_corr_df = pd.DataFrame({
    'gene': genes,
    'spearman_r': corrs,
    'p_value': pvals
})
gene_corr_df = gene_corr_df.dropna()
gene_corr_df['abs_r'] = gene_corr_df['spearman_r'].abs()
gene_corr_df = gene_corr_df.sort_values('abs_r', ascending=False)

print(f"Top positive: {gene_corr_df.iloc[0]['gene']} (r={gene_corr_df.iloc[0]['spearman_r']:.3f})")
print(f"Top negative: {gene_corr_df.iloc[-1]['gene']} (r={gene_corr_df.iloc[-1]['spearman_r']:.3f})")

gene_corr_df.to_csv(os.path.join(RESULT_DIR, 'pseudotime_gene_correlations.csv'), index=False)
print(f"Saved: pseudotime_gene_correlations.csv")

# ===================== 9. 关键基因可视化 =====================
print("\n" + "=" * 60)
print("Step 9: 关键基因可视化")
print("=" * 60)

key_genes = ['FOXP3', 'CCR8', 'MKI67', 'IL2RA', 'CTLA4', 'TIGIT', 'LAG3',
             'PDCD1', 'IL10', 'TGFB1', 'GZMB', 'PRF1', 'CXCR3', 'CCR4',
             'BATF', 'IRF4', 'STAT5A', 'STAT3', 'NFKB1', 'RELA',
             'TOX', 'NR4A1', 'EOMES', 'TBX21', 'RORC', 'GATA3']

available_genes = [g for g in key_genes if g in adata_treg.raw.var_names]
print(f"Available key genes: {available_genes}")

# 9.1 UMAP上展示关键基因表达
n_genes = min(len(available_genes), 12)
fig, axes = plt.subplots(3, 4, figsize=(20, 15))
axes = axes.flatten()
for i, gene in enumerate(available_genes[:n_genes]):
    sc.pl.umap(adata_treg, color=gene, ax=axes[i], show=False,
               title=gene, color_map='Reds', use_raw=True)
    axes[i].set_xlabel('UMAP1')
    axes[i].set_ylabel('UMAP2')
for j in range(i+1, 12):
    axes[j].axis('off')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '07_Key_genes_UMAP.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 07_Key_genes_UMAP.png")

# 9.2 关键基因随拟时序的表达趋势
fig, axes = plt.subplots(3, 4, figsize=(20, 15))
axes = axes.flatten()
for i, gene in enumerate(available_genes[:n_genes]):
    ax = axes[i]
    gene_idx = list(adata_treg.raw.var_names).index(gene)
    expr = X_raw[:, gene_idx]
    mask = expr > 0

    df_gene = pd.DataFrame({
        'pseudotime': pseudotime[mask],
        'expression': expr[mask],
        'subtype': adata_treg.obs['sub_cell_type'].values[mask]
    })

    for subtype, color in [('CD4T_Treg_FOXP3', '#1f77b4'),
                            ('CD4T_Treg_CCR8', '#ff7f0e'),
                            ('CD4T_Treg_MKI67', '#2ca02c')]:
        df_sub = df_gene[df_gene['subtype'] == subtype]
        if len(df_sub) > 0:
            ax.scatter(df_sub['pseudotime'], df_sub['expression'],
                      c=color, alpha=0.1, s=1, label=subtype)

    ax.set_xlabel('Pseudotime')
    ax.set_ylabel('Expression (log1p)')
    ax.set_title(gene)
    ax.legend(markerscale=5, fontsize=6)

for j in range(i+1, 12):
    axes[j].axis('off')

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '08_Genes_vs_pseudotime_trends.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 08_Genes_vs_pseudotime_trends.png")

# ===================== 10. 响应差异分析 =====================
print("\n" + "=" * 60)
print("Step 10: 响应者 vs 非响应者 拟时序差异")
print("=" * 60)

resp_pt = adata_treg.obs[adata_treg.obs['response'] == 'Responder']['dpt_pseudotime']
nonresp_pt = adata_treg.obs[adata_treg.obs['response'] == 'Non-responder']['dpt_pseudotime']

stat, pval = mannwhitneyu(resp_pt, nonresp_pt, alternative='two-sided')
print(f"Mann-Whitney U: p={pval:.2e}")
print(f"Responder median PT: {resp_pt.median():.3f}")
print(f"Non-responder median PT: {nonresp_pt.median():.3f}")

# 各亚型内的响应差异
print("\nSubtype-specific response differences:")
for subtype in ['CD4T_Treg_FOXP3', 'CD4T_Treg_CCR8', 'CD4T_Treg_MKI67']:
    sub_data = adata_treg.obs[adata_treg.obs['sub_cell_type'] == subtype]
    if len(sub_data) < 50:
        continue
    resp_sub = sub_data[sub_data['response'] == 'Responder']['dpt_pseudotime']
    nonresp_sub = sub_data[sub_data['response'] == 'Non-responder']['dpt_pseudotime']
    if len(resp_sub) > 10 and len(nonresp_sub) > 10:
        stat_sub, pval_sub = mannwhitneyu(resp_sub, nonresp_sub, alternative='two-sided')
        print(f"  {subtype}: Resp median={resp_sub.median():.3f}, "
              f"NonResp median={nonresp_sub.median():.3f}, p={pval_sub:.3f}")

# 小提琴图
fig, ax = plt.subplots(figsize=(10, 6))
df_violin = adata_treg.obs[['sub_cell_type', 'response', 'dpt_pseudotime']].copy()
sns.violinplot(data=df_violin, x='sub_cell_type', y='dpt_pseudotime', hue='response',
               palette={'Responder': '#d62728', 'Non-responder': '#9467bd'},
               split=True, ax=ax)
ax.set_xlabel('Treg Subtype')
ax.set_ylabel('Pseudotime')
ax.set_title('Pseudotime Distribution by Subtype and Response')
plt.xticks(rotation=15)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '09_Pseudotime_violin_response.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 09_Pseudotime_violin_response.png")

# ===================== 11. 保存结果 =====================
print("\n" + "=" * 60)
print("Step 11: 保存结果")
print("=" * 60)

dpt_results = adata_treg.obs[['sampleID', 'sub_cell_type', 'response',
                               'dpt_pseudotime', 'leiden', 'pathological_response']].copy()
dpt_results.to_csv(os.path.join(RESULT_DIR, 'Treg_pseudotime_per_cell.csv'))
print(f"Saved: Treg_pseudotime_per_cell.csv")

sample_dpt = adata_treg.obs.groupby('sampleID').agg({
    'dpt_pseudotime': ['mean', 'median', 'std'],
    'sub_cell_type': lambda x: x.value_counts().index[0],
    'response': 'first',
    'pathological_response': 'first'
}).reset_index()
sample_dpt.columns = ['sampleID', 'dpt_mean', 'dpt_median', 'dpt_std',
                      'dominant_subtype', 'response', 'pathological_response']
sample_dpt.to_csv(os.path.join(RESULT_DIR, 'Treg_pseudotime_per_sample.csv'), index=False)
print(f"Saved: Treg_pseudotime_per_sample.csv")

adata_treg.write_h5ad(os.path.join(RESULT_DIR, 'Treg_pseudotime_adata.h5ad'))
print(f"Saved: Treg_pseudotime_adata.h5ad")

# ===================== 12. 总结 =====================
print("\n" + "=" * 60)
print("分析完成!")
print("=" * 60)
print(f"输出目录: {OUT_DIR}")
print(f"结果目录: {RESULT_DIR}")
print(f"\n核心发现:")
print(f"- Treg细胞总数: {adata_treg.n_obs}")
print(f"- 拟时序范围: {adata_treg.obs['dpt_pseudotime'].min():.3f} - {adata_treg.obs['dpt_pseudotime'].max():.3f}")
print(f"- 响应者中位拟时序: {resp_pt.median():.3f}")
print(f"- 非响应者中位拟时序: {nonresp_pt.median():.3f}")
print(f"- 响应差异显著性: p={pval:.2e}")
