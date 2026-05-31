#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE207422 Monocle3 输入准备
从完整 92,330 cells × 24,292 genes 矩阵中提取细胞，构建 Monocle3 输入
"""

import os
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

# ===================== 配置 =====================
MATRIX_PATH = r"D:\research\potato\data\raw\GSE207422\GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
ANNOTATION_PATH = r"D:\research\potato\data\raw\GSE207422\annotation_results\all_cells_annotation.csv"
META_WEAK = r"D:\research\tomato\data\gse207422_cellchat_input\weak_metadata.csv"
META_STRONG = r"D:\research\tomato\data\gse207422_cellchat_input\strong_metadata.csv"
OUT_DIR = r"D:\research\tomato\data\gse207422_monocle3_input"

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("GSE207422 Monocle3 Input Preparation")
print("=" * 60)

# ===================== 1. 读取注释和分组信息 =====================
print("\n[Step 1] Reading annotations and metadata...")

ann = pd.read_csv(ANNOTATION_PATH, index_col=0)
ann.index.name = "cell"
print(f"Annotation: {ann.shape[0]} cells")

meta_weak = pd.read_csv(META_WEAK)
meta_strong = pd.read_csv(META_STRONG)

# 合并分组信息
meta_weak['response'] = 'Weak'
meta_strong['response'] = 'Strong'
meta_all = pd.concat([meta_weak, meta_strong], ignore_index=True)
meta_all = meta_all.set_index('cell')

# 合并注释和分组
ann = ann.join(meta_all[['cell_type', 'response_group']], how='left')
ann['response_group'] = ann['response_group'].fillna('Unknown')

# 重命名 Treg
ann['cell_type_major'] = ann['cell_type']
ann.loc[ann['predicted_labels'] == 'Regulatory T cells', 'cell_type_major'] = 'Treg'
ann.loc[ann['majority_voting'] == 'Regulatory T cells', 'cell_type_major'] = 'Treg'

print("\nCell type distribution:")
print(ann['cell_type_major'].value_counts())
print("\nResponse group distribution:")
print(ann['response_group'].value_counts())

# ===================== 2. 筛选细胞 =====================
print("\n[Step 2] Selecting cells for Monocle3...")

# 策略：保留主要细胞类型（数量>500），排除极稀有类型
# 但为了聚焦 Treg 轨迹，我们保留所有 Treg + CD4_T + Fibroblast + 其他主要类型
cell_counts = ann['cell_type_major'].value_counts()
major_types = cell_counts[cell_counts >= 500].index.tolist()

# 确保包含 Fibroblast（即使数量少）
if 'Fibroblast' not in major_types:
    major_types.append('Fibroblast')

selected_cells = ann[ann['cell_type_major'].isin(major_types)].index.tolist()
print(f"Selected {len(selected_cells)} cells from types: {major_types}")

# 保存 metadata
meta_out = ann.loc[selected_cells, ['cell_type_major', 'response_group', 'predicted_labels', 'majority_voting']].copy()
meta_out = meta_out.rename(columns={'cell_type_major': 'cell_type'})
meta_out.to_csv(os.path.join(OUT_DIR, "metadata.csv"))
print(f"Metadata saved: {len(meta_out)} cells")

# ===================== 3. 从原始矩阵提取表达 =====================
print("\n[Step 3] Extracting expression matrix...")
print(f"Source: {MATRIX_PATH}")
print("This may take a few minutes (4.5GB matrix)...")

# 读取 header，找到需要的列索引
with gzip.open(MATRIX_PATH, "rt") as f:
    header = f.readline().strip().split("\t")
    all_cells = header[1:]
    
    # 找到 selected_cells 在 header 中的位置
    cell_to_idx = {c: i for i, c in enumerate(all_cells)}
    selected_indices = []
    selected_in_order = []
    for c in selected_cells:
        if c in cell_to_idx:
            selected_indices.append(cell_to_idx[c])
            selected_in_order.append(c)
    
    print(f"Found {len(selected_indices)} / {len(selected_cells)} cells in matrix")
    
    # 读取基因和表达
    genes = []
    data_list = []
    
    for line_idx, line in enumerate(f):
        parts = line.strip().split("\t")
        gene = parts[0]
        genes.append(gene)
        
        # 只取需要的列
        expr = [float(parts[i+1]) for i in selected_indices]
        data_list.append(expr)
        
        if (line_idx + 1) % 5000 == 0:
            print(f"  Processed {line_idx + 1} / ~24292 genes", end="\r")
    
    print(f"\n  Processed {len(genes)} genes")

# 构建 AnnData
print("\nBuilding AnnData...")
data = np.array(data_list, dtype=np.float32)  # genes x cells
data_T = data.T  # cells x genes

adata = sc.AnnData(data_T)
adata.obs_names = selected_in_order
adata.var_names = genes
adata.obs = meta_out.loc[selected_in_order]

# 保存 h5ad
h5ad_path = os.path.join(OUT_DIR, "gse207422_selected.h5ad")
adata.write(h5ad_path)
print(f"Saved: {h5ad_path}")
print(f"Shape: {adata.shape}")

# ===================== 4. 导出 Monocle3 输入 =====================
print("\n[Step 4] Exporting Monocle3 inputs...")

# 4.1 稀疏矩阵
counts_sparse = csr_matrix(adata.X.T)  # genes x cells
mmwrite(os.path.join(OUT_DIR, "counts.mtx"), counts_sparse)
print("Saved: counts.mtx")

# 4.2 基因列表
genes_df = pd.DataFrame({'gene_short_name': adata.var_names})
genes_df.to_csv(os.path.join(OUT_DIR, "genes.csv"), index=False)
print("Saved: genes.csv")

# 4.3 Metadata（需要 cell_id 列）
meta_monocle = adata.obs.copy()
meta_monocle['cell_id'] = meta_monocle.index
meta_monocle = meta_monocle.reset_index(drop=True)
meta_monocle.to_csv(os.path.join(OUT_DIR, "metadata.csv"), index=False)
print("Saved: metadata.csv")

# ===================== 5. 统计 =====================
print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"Total cells: {adata.n_obs}")
print(f"Total genes: {adata.n_vars}")
print(f"Cell types:")
print(adata.obs['cell_type'].value_counts())
print(f"Response groups:")
print(adata.obs['response_group'].value_counts())
print(f"\nOutput directory: {OUT_DIR}")
