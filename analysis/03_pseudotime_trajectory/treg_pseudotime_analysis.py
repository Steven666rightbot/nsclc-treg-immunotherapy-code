#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE243013 - Treg细胞拟时序分析 (Pseudotime Analysis)
聚焦Treg三个亚型：FOXP3 / CCR8 / MKI67
使用 diffusion pseudotime (DPT) + PAGA 推断分化轨迹
比较响应者(MPR/pCR) vs 非响应者(non-MPR)的拟时序差异
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import sparse
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')
sc.settings.verbosity = 3  # verbosity level: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=150, facecolor='white', fontsize=10)

# ===================== 路径配置 =====================
DATA_DIR = r"D:\Research\土豆\数据\raw\GSE243013"
OUT_DIR = r"D:\Research\tomato\figures\pseudotime_treg"
RESULT_DIR = r"D:\Research\tomato\data\pseudotime_treg"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(RESULT_DIR, exist_ok=True)

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
print("Step 2: 读取表达矩阵 (22GB mtx)")
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

# 读取mtx
mtx_path = os.path.join(DATA_DIR, "GSE243013_NSCLC_immune_scRNA_counts.mtx")
print(f"Reading mtx: {mtx_path}")
print("This may take 2-5 minutes...")

adata = sc.read_mtx(mtx_path).T  # transpose: cells x genes
adata.obs_names = barcodes_list
adata.var_names = gene_names
print(f"Full adata: {adata.shape}")

# 筛选Treg细胞
treg_barcodes = treg_meta['cellID'].tolist()
adata_treg = adata[adata.obs_names.isin(treg_barcodes)].copy()
print(f"Treg adata: {adata_treg.shape}")

# 释放内存
del adata

# 合并metadata
obs_df = treg_meta.set_index('cellID')
adata_treg.obs = adata_treg.obs.join(obs_df)
print(f"Final Treg adata: {adata_treg.shape}")

# ===================== 3. 预处理 =====================
print("\n" + "=" * 60)
print("Step 3: 预处理")
print("=" * 60)

# 基础过滤：去掉在所有细胞中表达量都为0的基因
sc.pp.filter_genes(adata_treg, min_cells=10)
print(f"After filtering genes (min_cells=10): {adata_treg.shape}")

# 归一化到10000 counts per cell
sc.pp.normalize_total(adata_treg, target_sum=1e4)
sc.pp.log1p(adata_treg)

# 高变基因选择 (HVG)
sc.pp.highly_variable_genes(adata_treg, n_top_genes=2000, subset=True)
print(f"After HVG selection: {adata_treg.shape}")

# 缩放
sc.pp.scale(adata_treg, max_value=10)

# ===================== 4. 降维与聚类 =====================
print("\n" + "=" * 60)
print("Step 4: PCA + 邻居图 + UMAP + 聚类")
print("=" * 60)

sc.tl.pca(adata_treg, svd_solver='arpack')
sc.pp.neighbors(adata_treg, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_treg)
sc.tl.leiden(adata_treg, resolution=0.5)

print(f"Leiden clusters: {adata_treg.obs['leiden'].nunique()}")

# ===================== 5. PAGA + 轨迹推断 =====================
print("\n" + "=" * 60)
print("Step 5: PAGA + Diffusion Map + DPT")
print("=" * 60)

# PAGA: 基于leiden聚类构建轨迹图
sc.tl.paga(adata_treg, groups='sub_cell_type')

# Diffusion Map
sc.tl.diffmap(adata_treg, n_comps=15)

# DPT (Diffusion Pseudotime)
# 以FOXP3+ Treg作为根节点（生物学上认为是naive/resting状态）
foxp3_cells = adata_treg.obs[adata_treg.obs['sub_cell_type'] == 'CD4T_Treg_FOXP3'].index
if len(foxp3_cells) > 0:
    # 设置根节点为FOXP3+ Treg中diffusion map坐标最小的细胞（最"初始"的）
    adata_treg.uns['iroot'] = np.argmin(adata_treg.obsm['X_diffmap'][:, 1])
    # 或者更精确：用FOXP3+群体中diffusion coordinate最小的
    foxp3_idx = np.where(adata_treg.obs['sub_cell_type'] == 'CD4T_Treg_FOXP3')[0]
    root_idx = foxp3_idx[np.argmin(adata_treg.obsm['X_diffmap'][foxp3_idx, 1])]
    adata_treg.uns['iroot'] = root_idx
    print(f"Root set to FOXP3+ Treg cell index {root_idx}")

sc.tl.dpt(adata_treg, n_dcs=10)
print(f"DPT range: {adata_treg.obs['dpt_pseudotime'].min():.3f} - {adata_treg.obs['dpt_pseudotime'].max():.3f}")

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

# 6.2 UMAP: 响应状态
sc.pl.umap(adata_treg, color='response', ax=axes[1], show=False, title='Response',
           legend_loc='right margin', palette={'Responder': '#d62728', 'Non-responder': '#9467bd'})
axes[1].set_xlabel('UMAP1')
axes[1].set_ylabel('UMAP2')

# 6.3 UMAP: 拟时序
sc.pl.umap(adata_treg, color='dpt_pseudotime', ax=axes[2], show=False, title='Pseudotime',
           color_map='viridis')
axes[2].set_xlabel('UMAP1')
axes[2].set_ylabel('UMAP2')

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '01_Treg_UMAP_overview.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 01_Treg_UMAP_overview.png")

# 6.4 PAGA图
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.paga(adata_treg, color='sub_cell_type', ax=ax, show=False,
           title='PAGA: Treg Trajectory',
           node_size_scale=3, edge_width_scale=1)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '02_Treg_PAGA.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 02_Treg_PAGA.png")

# 6.5 UMAP + PAGA
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

# 7.3 箱线图：拟时序 vs 亚型 × 响应
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

# 计算每个基因与拟时序的Spearman相关性
# 使用原始log1p表达值（在adata_treg.raw中）
# 先保存raw
adata_treg.raw = adata_treg

# 恢复raw数据用于计算
X = adata_treg.raw.X
if sparse.issparse(X):
    X = X.toarray()

genes = adata_treg.raw.var_names.tolist()
pseudotime = adata_treg.obs['dpt_pseudotime'].values

print(f"Computing correlations for {len(genes)} genes...")
corrs = []
pvals = []
for i, gene in enumerate(genes):
    expr = X[:, i]
    # 只在表达>0的细胞中计算
    mask = expr > 0
    if mask.sum() < 50:  # 至少50个细胞表达
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

print(f"Top positive correlation: {gene_corr_df.iloc[0]['gene']} (r={gene_corr_df.iloc[0]['spearman_r']:.3f})")
print(f"Top negative correlation: {gene_corr_df.iloc[-1]['gene']} (r={gene_corr_df.iloc[-1]['spearman_r']:.3f})")

# 保存结果
gene_corr_df.to_csv(os.path.join(RESULT_DIR, 'pseudotime_gene_correlations.csv'), index=False)
print(f"Saved: {RESULT_DIR}\\pseudotime_gene_correlations.csv")

# ===================== 9. 关键基因可视化 =====================
print("\n" + "=" * 60)
print("Step 9: 关键基因可视化")
print("=" * 60)

# 选择与Treg功能相关的关键marker基因
key_genes = ['FOXP3', 'CCR8', 'MKI67', 'IL2RA', 'CTLA4', 'TIGIT', 'LAG3',
             'PDCD1', 'IL10', 'TGFB1', 'GZMB', 'PRF1', 'CXCR3', 'CCR4',
             'BATF', 'IRF4', 'STAT5A', 'STAT3', 'NFKB1', 'RELA',
             'TOX', 'NR4A1', 'EOMES', 'TBX21', 'RORC', 'GATA3']

# 只保留在数据中的基因
available_genes = [g for g in key_genes if g in adata_treg.raw.var_names]
print(f"Available key genes: {available_genes}")

# 9.1 UMAP上展示关键基因表达
n_genes = len(available_genes)
fig, axes = plt.subplots(3, 4, figsize=(20, 15))
axes = axes.flatten()
for i, gene in enumerate(available_genes[:12]):
    sc.pl.umap(adata_treg, color=gene, ax=axes[i], show=False, 
               title=gene, color_map='Reds', use_raw=True)
    axes[i].set_xlabel('UMAP1')
    axes[i].set_ylabel('UMAP2')
for j in range(i+1, 12):
    axes[j].axis('off')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '07_Key_genes_UMAP_1.png'), dpi=300, bbox_inches='tight')
plt.close()

if n_genes > 12:
    fig, axes = plt.subplots(3, 4, figsize=(20, 15))
    axes = axes.flatten()
    for i, gene in enumerate(available_genes[12:24]):
        sc.pl.umap(adata_treg, color=gene, ax=axes[i], show=False, 
                   title=gene, color_map='Reds', use_raw=True)
        axes[i].set_xlabel('UMAP1')
        axes[i].set_ylabel('UMAP2')
    for j in range(i+1, 12):
        axes[j].axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, '07_Key_genes_UMAP_2.png'), dpi=300, bbox_inches='tight')
    plt.close()
print("Saved: 07_Key_genes_UMAP.png")

# 9.2 关键基因随拟时序的表达趋势
fig, axes = plt.subplots(3, 4, figsize=(20, 15))
axes = axes.flatten()
for i, gene in enumerate(available_genes[:12]):
    ax = axes[i]
    # 获取表达值
    gene_idx = list(adata_treg.raw.var_names).index(gene)
    expr = X[:, gene_idx]
    mask = expr > 0
    
    # 散点图 + LOESS平滑
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
    
    # LOESS平滑
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess
        for subtype in ['CD4T_Treg_FOXP3', 'CD4T_Treg_CCR8', 'CD4T_Treg_MKI67']:
            df_sub = df_gene[df_gene['subtype'] == subtype]
            if len(df_sub) > 100:
                smooth = lowess(df_sub['expression'], df_sub['pseudotime'], 
                               frac=0.3, return_sorted=True)
                ax.plot(smooth[:, 0], smooth[:, 1], 
                       color={'CD4T_Treg_FOXP3': '#1f77b4', 
                              'CD4T_Treg_CCR8': '#ff7f0e',
                              'CD4T_Treg_MKI67': '#2ca02c'}[subtype],
                       linewidth=2)
    except ImportError:
        pass
    
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

# 10.1 统计检验
from scipy.stats import mannwhitneyu, kruskal

resp_pt = adata_treg.obs[adata_treg.obs['response'] == 'Responder']['dpt_pseudotime']
nonresp_pt = adata_treg.obs[adata_treg.obs['response'] == 'Non-responder']['dpt_pseudotime']

stat, pval = mannwhitneyu(resp_pt, nonresp_pt, alternative='two-sided')
print(f"Mann-Whitney U: p={pval:.2e}")
print(f"Responder median PT: {resp_pt.median():.3f}")
print(f"Non-responder median PT: {nonresp_pt.median():.3f}")

# 10.2 各亚型内的响应差异
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

# 10.3 小提琴图
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

# 保存细胞级别的拟时序信息
dpt_results = adata_treg.obs[['sampleID', 'sub_cell_type', 'response', 
                               'dpt_pseudotime', 'leiden', 'pathological_response']].copy()
dpt_results.to_csv(os.path.join(RESULT_DIR, 'Treg_pseudotime_per_cell.csv'))
print(f"Saved: {RESULT_DIR}\\Treg_pseudotime_per_cell.csv")

# 保存样本级别的拟时序统计
sample_dpt = adata_treg.obs.groupby('sampleID').agg({
    'dpt_pseudotime': ['mean', 'median', 'std'],
    'sub_cell_type': lambda x: x.value_counts().index[0],
    'response': 'first',
    'pathological_response': 'first'
}).reset_index()
sample_dpt.columns = ['sampleID', 'dpt_mean', 'dpt_median', 'dpt_std', 
                      'dominant_subtype', 'response', 'pathological_response']
sample_dpt.to_csv(os.path.join(RESULT_DIR, 'Treg_pseudotime_per_sample.csv'), index=False)
print(f"Saved: {RESULT_DIR}\\Treg_pseudotime_per_sample.csv")

# 保存AnnData
adata_treg.write_h5ad(os.path.join(RESULT_DIR, 'Treg_pseudotime_adata.h5ad'))
print(f"Saved: {RESULT_DIR}\\Treg_pseudotime_adata.h5ad")

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
