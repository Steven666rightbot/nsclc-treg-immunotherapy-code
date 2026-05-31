#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE243013 - Treg数据预处理并导出为Monocle3输入格式
输出：counts.mtx + metadata.csv + genes.csv
"""

import os
import time
import numpy as np
import pandas as pd
from scipy import sparse

# ===================== 路径配置 =====================
DATA_DIR = r"D:\Research\土豆\数据\raw\GSE243013"
OUT_DIR = r"D:\Research\tomato\data\monocle3_input"
os.makedirs(OUT_DIR, exist_ok=True)

# ===================== 1. 读取Treg metadata =====================
print("=" * 60)
print("Step 1: 读取Treg metadata")
print("=" * 60)

treg_meta = pd.read_csv(
    os.path.join(DATA_DIR, "GSE243013_Treg_cells_metadata.csv"),
    low_memory=False
)
print(f"Treg cells: {len(treg_meta)}")

# 过滤unknown响应
treg_meta = treg_meta[treg_meta['pathological_response'] != 'unknowm'].copy()
print(f"过滤unknown后: {len(treg_meta)}")

# 定义响应标签
treg_meta['response'] = treg_meta['pathological_response'].map({
    'MPR': 'Responder', 'pCR': 'Responder', 'non-MPR': 'Non-responder'
})

# 只保留关键列用于monocle3
cols_keep = ['cellID', 'sampleID', 'sub_cell_type', 'pathological_response', 
             'response', 'n_genes', 'total_counts', 'pct_counts_mt']
meta_out = treg_meta[cols_keep].copy()

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

# 构建Treg行掩码
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
    mask = is_treg[chunk[0].values]
    kept = mask.sum()
    if kept > 0:
        chunk_rows.append(chunk[0].values[mask] - 1)
        chunk_cols.append(chunk[1].values[mask] - 1)
        chunk_data.append(chunk[2].values[mask])
        total_kept += kept
    chunk_count += 1
    elapsed = time.time() - start_time
    progress = chunk_count / 9
    eta = elapsed / progress - elapsed if progress > 0 else 0
    print(f"  Chunk {chunk_count:02d}/~09: kept {kept:>8,}/{len(chunk):>10,} "
          f"({kept/len(chunk)*100:5.2f}%) | elapsed: {elapsed/60:5.1f}min | ETA: {eta/60:5.1f}min")

print(f"\nTotal chunks: {chunk_count}")
print(f"Total kept entries: {total_kept}")
print(f"Read time: {(time.time()-start_time)/60:.1f} minutes")

# 合并数据
print("\nMerging chunks...")
all_rows = np.concatenate(chunk_rows)
all_cols = np.concatenate(chunk_cols)
all_data = np.concatenate(chunk_data)

# 构建稀疏矩阵并转CSR
X = sparse.coo_matrix(
    (all_data, (all_rows, all_cols)),
    shape=(len(barcodes_list) - 1, len(gene_names))
)
print(f"COO matrix: {X.shape}")
X = X.tocsr()
print(f"CSR conversion done")

# 筛选Treg细胞并按metadata顺序排列
treg_cell_ids = meta_out['cellID'].tolist()
barcodes_array = np.array(barcodes_list[1:])

print("Building barcode index map...")
barcode_to_idx = {bc: i for i, bc in enumerate(barcodes_array)}

ordered_indices = []
ordered_cell_ids = []
for cid in treg_cell_ids:
    if cid in barcode_to_idx:
        ordered_indices.append(barcode_to_idx[cid])
        ordered_cell_ids.append(cid)

ordered_indices = np.array(ordered_indices)
print(f"Ordered Treg cells: {len(ordered_indices)}")

X_treg = X[ordered_indices, :]
print(f"Treg matrix shape: {X_treg.shape}")

# ===================== 3. 预处理 =====================
print("\n" + "=" * 60)
print("Step 3: 预处理")
print("=" * 60)

import scanpy as sc
adata = sc.AnnData(X_treg.tocsr())
adata.obs_names = ordered_cell_ids
adata.var_names = gene_names

# 合并metadata
obs_df = meta_out.set_index('cellID')
adata.obs = adata.obs.join(obs_df)

print(f"Raw adata: {adata.shape}")

# 基础过滤：只去掉完全不在任何细胞中表达的基因
sc.pp.filter_genes(adata, min_cells=1)
print(f"After filter_genes (min_cells=1): {adata.shape}")

# 过滤总表达量为0的细胞
print("Filtering zero-read cells...")
n_before = adata.n_obs
adata = adata[adata.X.sum(axis=1).A1 > 0, :].copy()
n_after = adata.n_obs
print(f"Removed {n_before - n_after} zero-read cells. Final: {n_after}")

# 归一化 + log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 不做HVG选择，保留全部过滤后的基因
# Monocle3的preprocess_cds内部会自己选HVG
print(f"Retaining all filtered genes: {adata.shape}")

# ===================== 4. 导出 =====================
print("\n" + "=" * 60)
print("Step 4: 导出为Monocle3格式")
print("=" * 60)

# 导出稀疏矩阵 (MTX格式)
# Monocle3 的 new_cell_data_set 支持 Matrix 包读取的稀疏矩阵
# 我们需要基因 x 细胞的格式 (转置)
print("Exporting expression matrix (genes x cells)...")
X_out = adata.X.T  # genes x cells
if sparse.issparse(X_out):
    X_out = X_out.tocsr()

# 保存为MTX (Matrix Market格式)
# scipy.io.mmwrite 输出的是坐标格式
from scipy.io import mmwrite
mmwrite(os.path.join(OUT_DIR, "counts.mtx"), X_out)
print(f"Saved: counts.mtx ({X_out.shape})")

# 导出metadata
meta_export = adata.obs.reset_index(drop=True)
meta_export['cell_id'] = adata.obs_names
# 确保列名符合monocle3要求
meta_export.to_csv(os.path.join(OUT_DIR, "metadata.csv"), index=False)
print(f"Saved: metadata.csv ({len(meta_export)} cells)")

# 导出基因列表
gene_df = pd.DataFrame({'gene_short_name': adata.var_names})
gene_df.to_csv(os.path.join(OUT_DIR, "genes.csv"), index=False)
print(f"Saved: genes.csv ({len(gene_df)} genes)")

# 额外导出：仅Treg细胞ID列表（用于验证）
pd.DataFrame({'cell_id': adata.obs_names}).to_csv(
    os.path.join(OUT_DIR, "cell_ids.csv"), index=False)
print("Saved: cell_ids.csv")

# ===================== 5. 保存h5ad（备选） =====================
print("\nSaving h5ad for Python backup...")
adata.write_h5ad(os.path.join(OUT_DIR, "treg_processed.h5ad"))
print("Saved: treg_processed.h5ad")

# ===================== 总结 =====================
print("\n" + "=" * 60)
print("导出完成!")
print("=" * 60)
print(f"输出目录: {OUT_DIR}")
print(f"\nMonocle3 输入文件:")
print(f"  1. counts.mtx      - 表达矩阵 (genes x cells, sparse)")
print(f"  2. metadata.csv    - 细胞元数据")
print(f"  3. genes.csv       - 基因列表")
print(f"\n细胞数: {adata.n_obs}")
print(f"基因数: {adata.n_vars}")
print(f"\n亚型分布:")
print(adata.obs['sub_cell_type'].value_counts())
print(f"\n响应分布:")
print(adata.obs['response'].value_counts())
