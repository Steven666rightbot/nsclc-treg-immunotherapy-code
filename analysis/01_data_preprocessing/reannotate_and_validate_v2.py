#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE207422 重新注释 + XGBoost验证 (内存优化版)
"""

import os
import gc
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
import xgboost as xgb

warnings.filterwarnings('ignore')
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ========================== 配置 ==========================
REF_META = r"D:\Research\土豆\数据\raw\GSE243013\GSE243013_NSCLC_immune_scRNA_metadata.csv"
REF_BARCODES = r"D:\Research\土豆\数据\raw\GSE243013\GSE243013_barcodes.csv.gz"
REF_FEATURES = r"D:\Research\土豆\数据\raw\GSE243013\features.csv.gz"
REF_MTX = r"D:\Research\土豆\数据\raw\GSE243013\GSE243013_NSCLC_immune_scRNA_counts.mtx"

VAL_UMI = r"D:\Research\potato\data\raw\GSE207422\GSE207422_NSCLC_scRNAseq_UMI_matrix.txt"
VAL_META = r"D:\Research\potato\data\raw\GSE207422\GSE207422_NSCLC_scRNAseq_metadata.xlsx"
VAL_ANNO = r"D:\Research\potato\data\raw\GSE207422\annotation_results\all_cells_annotation.csv"

OUT_DIR = r"D:\Research\tomato\figures"
DATA_DIR = r"D:\Research\tomato\data"

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

N_CELLS_PER_TYPE = 500  # 每类抽样的细胞数（减少以节省内存）
CHUNK_SIZE = 5_000_000  # mtx分块读取行数

print("="*70)
print("Step 1: 构建GSE243013参考表达谱")
print("="*70)

# ---- 1.1 读取metadata ----
print("\n[1.1] 读取GSE243013 metadata...")
meta = pd.read_csv(REF_META, low_memory=False)
meta = meta[meta['pathological_response'] != 'unknowm'].copy()
print(f"  有效细胞数: {len(meta)}")

# ---- 1.2 读取barcodes和features ----
print("\n[1.2] 读取barcodes和features...")
barcodes = pd.read_csv(REF_BARCODES, header=None)[0].tolist()[1:]  # skip header
features = pd.read_csv(REF_FEATURES, header=None)[0].tolist()[1:]  # skip header
print(f"  Barcodes: {len(barcodes)}, Genes: {len(features)}")

# ---- 1.3 分层抽样 ----
print("\n[1.3] 分层抽样...")
np.random.seed(42)
sampled_cell_type = {}  # index → sub_cell_type
type_counts = meta['sub_cell_type'].value_counts()

for cell_type, count in type_counts.items():
    cells = meta[meta['sub_cell_type'] == cell_type]['cellID'].values
    n_sample = min(N_CELLS_PER_TYPE, len(cells))
    selected = np.random.choice(cells, size=n_sample, replace=False)
    for cid in selected:
        idx = barcodes.index(cid) + 1
        sampled_cell_type[idx] = cell_type

print(f"  抽样细胞数: {len(sampled_cell_type)}")

# ---- 1.4 扫描mtx文件 (pandas分块读取) ----
print(f"\n[1.4] 扫描mtx文件 (pandas chunksize={CHUNK_SIZE})...")
print("  这可能需要5-15分钟...")

from collections import defaultdict
cell_type_data = defaultdict(lambda: defaultdict(list))

current_row = None
current_genes = []
current_vals = []

processed_cells = 0
kept_cells = 0
chunk_idx = 0

# 手动跳过头部
with open(REF_MTX, 'r') as f:
    for line in f:
        if line.startswith('%'):
            continue
        parts = line.strip().split()
        if len(parts) == 3:
            n_rows, n_cols, n_nnz = int(parts[0]), int(parts[1]), int(parts[2])
            print(f"  mtx: {n_rows} rows, {n_cols} cols, {n_nnz} nonzeros")
            break

# 用pandas分块读取
reader = pd.read_csv(
    REF_MTX, skiprows=2, sep=' ', header=None,
    names=['row', 'col', 'val'], dtype={'row': np.int32, 'col': np.int32, 'val': np.float32},
    chunksize=CHUNK_SIZE
)

for chunk in reader:
    chunk_idx += 1
    rows = chunk['row'].values
    cols = chunk['col'].values
    vals = chunk['val'].values
    
    # 找到这一chunk中row变化的位置
    row_changes = np.where(np.diff(rows, prepend=rows[0]-1) != 0)[0]
    
    # 处理每个cell
    for i in range(len(row_changes)):
        start = row_changes[i]
        end = row_changes[i+1] if i+1 < len(row_changes) else len(rows)
        row_idx = rows[start]
        
        if row_idx != current_row:
            # 保存上一个cell
            if current_row is not None and current_row in sampled_cell_type:
                ct = sampled_cell_type[current_row]
                for g, v in zip(current_genes, current_vals):
                    cell_type_data[ct][g].append(v)
                kept_cells += 1
            
            current_row = row_idx
            current_genes = list(cols[start:end])
            current_vals = list(vals[start:end])
            processed_cells += 1
        else:
            current_genes.extend(cols[start:end])
            current_vals.extend(vals[start:end])
    
    if chunk_idx % 20 == 0:
        print(f"    chunk {chunk_idx}, 处理cells: {processed_cells}, 保留: {kept_cells}")

# 处理最后一个cell
if current_row is not None and current_row in sampled_cell_type:
    ct = sampled_cell_type[current_row]
    for g, v in zip(current_genes, current_vals):
        cell_type_data[ct][g].append(v)
    kept_cells += 1

print(f"  完成! 处理cells: {processed_cells}, 保留cells: {kept_cells}")

# ---- 1.5 计算pseudo-bulk ----
print("\n[1.5] 计算pseudo-bulk参考谱...")
n_genes = len(features)
cell_types = sorted(cell_type_data.keys())
reference_matrix = np.zeros((len(cell_types), n_genes))

for i, ct in enumerate(cell_types):
    gene_data = cell_type_data[ct]
    n_cells_in_type = len([k for k, v in sampled_cell_type.items() if v == ct])
    for g_idx, vals in gene_data.items():
        reference_matrix[i, g_idx - 1] = sum(vals) / n_cells_in_type

reference_matrix = np.log1p(reference_matrix)
print(f"  参考谱: {reference_matrix.shape}")

ref_df = pd.DataFrame(reference_matrix, index=cell_types, columns=features)
ref_df.to_csv(os.path.join(DATA_DIR, 'GSE243013_reference_pseudobulk.csv'))
print(f"  参考谱已保存")

del cell_type_data, meta, chunk, reader
gc.collect()

# ========================== Step 2: 读取GSE207422 (稀疏格式) ==========================
print("\n" + "="*70)
print("Step 2: 读取GSE207422并转为稀疏矩阵")
print("="*70)

print("\n[2.1] 读取header并识别术后样本...")
val_meta = pd.read_excel(VAL_META)
val_meta = val_meta[val_meta['Sample'].notna()].copy()
val_meta_post = val_meta[val_meta['Resource'] == 'Post-treatment surgery'].copy()
val_meta_post['response'] = val_meta_post['Pathologic Response'].map({'MPR': 1, 'pCR': 1, 'NMPR': 0})
val_meta_post = val_meta_post[val_meta_post['response'].notna()].copy()
post_samples = set(val_meta_post['Sample'].tolist())
print(f"  术后可评估样本: {post_samples}")

# 读取header获取所有细胞
with open(VAL_UMI, 'r') as f:
    header = f.readline().strip().split('\t')
all_cells = header[1:]
print(f"  总细胞数: {len(all_cells)}")

# 识别术后样本的细胞
cell_samples = ['_'.join(c.split('_')[:2]) for c in all_cells]
post_mask = [cs in post_samples for cs in cell_samples]
post_indices = [i for i, m in enumerate(post_mask) if m]
post_cells = [all_cells[i] for i in post_indices]
print(f"  术后样本细胞数: {len(post_cells)}")

# 读取基因名并逐行构建稀疏矩阵
print("\n[2.2] 逐行读取UMI矩阵，构建稀疏矩阵...")
row_indices = []
col_indices = []
data_values = []
gene_names = []

with open(VAL_UMI, 'r') as f:
    f.readline()  # skip header
    row_idx = 0
    for line in f:
        parts = line.strip().split('\t')
        gene = parts[0]
        values = np.array(parts[1:], dtype=np.float32)
        post_vals = values[post_indices]
        
        # 只保留非零值
        nonzero_mask = post_vals > 0
        if nonzero_mask.any():
            nz_cols = np.where(nonzero_mask)[0]
            row_indices.extend([row_idx] * len(nz_cols))
            col_indices.extend(nz_cols.tolist())
            data_values.extend(post_vals[nz_cols].tolist())
        
        gene_names.append(gene)
        row_idx += 1
        if row_idx % 5000 == 0:
            print(f"    已处理 {row_idx} 个基因...")

n_genes_val = len(gene_names)
n_cells_post = len(post_cells)
val_sparse = sp.csr_matrix((data_values, (row_indices, col_indices)), shape=(n_genes_val, n_cells_post))
print(f"  稀疏矩阵: {val_sparse.shape}, nonzeros: {val_sparse.nnz}")
print(f"  稀疏度: {val_sparse.nnz / (val_sparse.shape[0] * val_sparse.shape[1]) * 100:.2f}%")

# 释放列表内存
del row_indices, col_indices, data_values
gc.collect()

# ========================== Step 3: 基因对齐并注释 ==========================
print("\n" + "="*70)
print("Step 3: 基因对齐并映射注释")
print("="*70)

# ---- 3.1 基因对齐 ----
print("\n[3.1] 基因对齐...")
common_genes = list(set(features) & set(gene_names))
print(f"  交集基因: {len(common_genes)}")

ref_idx = [features.index(g) for g in common_genes]
val_idx = [gene_names.index(g) for g in common_genes]

ref_X = reference_matrix[:, ref_idx]  # 51 x common_genes
val_X = val_sparse[val_idx, :].T      # cells x common_genes (CSR -> CSC for column access)
val_X = val_X.tocsr()
print(f"  参考谱: {ref_X.shape}, 验证集: {val_X.shape}")

# ---- 3.2 CPM + log1p 标准化 ----
print("\n[3.2] CPM标准化...")
val_cpm = val_sparse[val_idx, :].T.tocsr()
# 每行(cell)除以总和，乘1e6
row_sums = np.array(val_cpm.sum(axis=1)).flatten()
row_sums[row_sums == 0] = 1  # 避免除0
# 手动做CPM
val_cpm_data = val_cpm.data.copy()
val_cpm_indices = val_cpm.indices.copy()
val_cpm_indptr = val_cpm.indptr.copy()

for i in range(val_cpm.shape[0]):
    start = val_cpm_indptr[i]
    end = val_cpm_indptr[i+1]
    if end > start:
        val_cpm_data[start:end] = np.log1p(val_cpm_data[start:end] / row_sums[i] * 1e6)

val_cpm = sp.csr_matrix((val_cpm_data, val_cpm_indices, val_cpm_indptr), shape=val_cpm.shape)
print(f"  CPM矩阵: {val_cpm.shape}")

# 将稀疏矩阵转为dense（交集基因通常几百到几千个，dense矩阵不大）
val_cpm_dense = val_cpm.toarray()
print(f"  Dense CPM: {val_cpm_dense.shape}, 内存: {val_cpm_dense.nbytes / 1e6:.1f} MB")

# ---- 3.3 余弦相似度映射 ----
print("\n[3.3] 余弦相似度映射...")
similarities = cosine_similarity(val_cpm_dense, ref_X)
predicted_types_idx = np.argmax(similarities, axis=1)
predicted_types = [cell_types[i] for i in predicted_types_idx]
max_similarities = np.max(similarities, axis=1)

print(f"  平均最大相似度: {max_similarities.mean():.4f}")
print(f"  中位最大相似度: {np.median(max_similarities):.4f}")

val_anno = pd.DataFrame({
    'cell': post_cells,
    'sub_cell_type': predicted_types,
    'cosine_similarity': max_similarities
})

print("\n  预测细胞类型分布 (Top 15):")
print(val_anno['sub_cell_type'].value_counts().head(15))

val_anno.to_csv(os.path.join(DATA_DIR, 'GSE207422_reannotated.csv'), index=False)

# ========================== Step 4: 计算比例并验证 ==========================
print("\n" + "="*70)
print("Step 4: 计算细胞比例并验证模型")
print("="*70)

# ---- 4.1 计算细胞比例 ----
print("\n[4.1] 计算细胞比例...")
val_anno['sample'] = val_anno['cell'].str.split('_').str[0] + '_' + val_anno['cell'].str.split('_').str[1]

val_counts = val_anno.groupby(['sample', 'sub_cell_type']).size().unstack(fill_value=0)
val_props = val_counts.div(val_counts.sum(axis=1), axis=0)
print(f"  术后样本细胞比例矩阵: {val_props.shape}")

# 对齐到训练集的51类
train_features = pd.read_csv(os.path.join(DATA_DIR, 'cell_proportions.csv'), index_col=0).columns.tolist()
for f in train_features:
    if f not in val_props.columns:
        val_props[f] = 0.0
val_props = val_props[train_features]
val_props.to_csv(os.path.join(DATA_DIR, 'GSE207422_cell_proportions_post.csv'))

# ---- 4.2 训练最终模型 ----
print("\n[4.2] 训练最终模型...")
train_props = pd.read_csv(os.path.join(DATA_DIR, 'cell_proportions.csv'), index_col=0)
train_labels = pd.read_csv(os.path.join(DATA_DIR, 'sample_labels.csv'), index_col=0)['response']
train_props = train_props.loc[train_labels.index]

X_train = train_props.values
y_train = train_labels.values

scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)

final_model = xgb.XGBClassifier(
    n_estimators=200, max_depth=4, learning_rate=0.05,
    subsample=0.8, colsample_bytree=0.8,
    reg_alpha=0.1, reg_lambda=1.0,
    random_state=42, eval_metric='logloss',
    use_label_encoder=False, n_jobs=4
)
final_model.fit(X_train_s, y_train)
print("  最终模型训练完成")

# ---- 4.3 验证 ----
print("\n[4.3] 验证集预测...")
X_val = val_props.values
sample_list = val_props.index.tolist()
y_val = val_meta_post.set_index('Sample').loc[sample_list, 'response'].values

X_val_s = scaler.transform(X_val)
y_prob = final_model.predict_proba(X_val_s)[:, 1]
y_pred = final_model.predict(X_val_s)

auc = roc_auc_score(y_val, y_prob)
acc = accuracy_score(y_val, y_pred)
print(f"\n  验证集 AUC: {auc:.4f}")
print(f"  验证集 ACC: {acc:.4f}")
print(f"\n  样本预测详情:")
for s, p, t in zip(sample_list, y_prob, y_val):
    resp_label = {1: 'MPR/pCR', 0: 'NMPR'}.get(int(t), 'Unknown')
    print(f"    {s}: prob={p:.3f}, pred={int(p>0.5)}, true={resp_label}")

# ========================== Step 5: 可视化 ==========================
print("\n" + "="*70)
print("Step 5: 生成可视化")
print("="*70)

# ---- 5.1 ROC曲线 ----
fig, ax = plt.subplots(figsize=(8, 8))
train_cv = pd.read_csv(os.path.join(DATA_DIR, 'cv_results.csv'))
mean_auc_train = train_cv['auc'].mean()
std_auc_train = train_cv['auc'].std()

fpr, tpr, _ = roc_curve(y_val, y_prob)
ax.plot(fpr, tpr, color='darkred', lw=2.5,
        label=f'GSE207422 Validation (AUC = {auc:.3f}, n={len(y_val)})')
ax.plot([0, 1], [0, 1], 'k--', lw=1.5, label='Chance')
ax.set_xlim([-0.02, 1.02])
ax.set_ylim([-0.02, 1.02])
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title(f'Model Validation on GSE207422\nTraining AUC = {mean_auc_train:.3f} ± {std_auc_train:.3f}',
             fontsize=14, fontweight='bold')
ax.legend(loc='lower right', fontsize=11)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '08_Validation_ROC_GSE207422.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [已保存] 08_Validation_ROC_GSE207422.png")

# ---- 5.2 细胞比例对比 ----
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
responders = val_meta_post[val_meta_post['response'] == 1]['Sample'].tolist()
non_responders = val_meta_post[val_meta_post['response'] == 0]['Sample'].tolist()

top10 = val_props.mean().sort_values(ascending=False).head(10).index.tolist()
resp_mean = val_props.loc[val_props.index.intersection(responders), top10].mean()
nonresp_mean = val_props.loc[val_props.index.intersection(non_responders), top10].mean()

x = np.arange(len(top10))
width = 0.35
axes[0].bar(x - width/2, resp_mean.values, width, label='Responder (MPR/pCR)', color='coral')
axes[0].bar(x + width/2, nonresp_mean.values, width, label='Non-Responder (NMPR)', color='steelblue')
axes[0].set_xticks(x)
axes[0].set_xticklabels(top10, rotation=45, ha='right', fontsize=9)
axes[0].set_ylabel('Mean Proportion', fontsize=11)
axes[0].set_title('Top 10 Cell Type Proportions (GSE207422)', fontsize=12, fontweight='bold')
axes[0].legend()
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

val_props[top10].T.plot(kind='barh', stacked=True, ax=axes[1], colormap='tab20', legend=False, width=0.7)
axes[1].set_xlabel('Proportion', fontsize=11)
axes[1].set_title('Cell Type Composition per Sample', fontsize=12, fontweight='bold')
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '09_Validation_celltype_composition.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [已保存] 09_Validation_celltype_composition.png")

# ---- 5.3 新旧注释对比 ----
print("\n[5.3] 新旧注释对比...")
old_anno = pd.read_csv(VAL_ANNO)
old_anno['sample'] = old_anno['cell'].str.split('_').str[0] + '_' + old_anno['cell'].str.split('_').str[1]
old_anno = old_anno[old_anno['sample'].isin(post_samples)]

major_map = {
    'T/NK cell': ['CD4T_', 'CD8T_', 'NK_', 'T_gdT_', 'ILC3_'],
    'B cell': ['Bm_', 'Bn_', 'Plasma_cell', 'B_prf_', 'GCB_'],
    'Myeloid cell': ['Mφ_', 'Mast cell', 'cDC', 'mDC', 'pDC', 'Neu_', 'M_']
}

def get_major(cell_type):
    for major, prefixes in major_map.items():
        for p in prefixes:
            if str(cell_type).startswith(p):
                return major
    return 'Other'

val_anno['major_new'] = val_anno['sub_cell_type'].apply(get_major)

new_major = val_anno.groupby(['sample', 'major_new']).size().unstack(fill_value=0)
new_major = new_major.div(new_major.sum(axis=1), axis=0)

old_major_keep = ['Tem/Trm cytotoxic T cells', 'Epithelial cells', 'Regulatory T cells',
                  'Alveolar macrophages', 'Classical monocytes', 'Memory B cells',
                  'Tcm/Naive helper T cells', 'Macrophages', 'Plasma cells', 'NK cells']
old_major = old_anno[old_anno['predicted_labels'].isin(old_major_keep)].groupby(
    ['sample', 'predicted_labels']).size().unstack(fill_value=0)
old_major = old_major.div(old_major.sum(axis=1), axis=0)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
new_major.plot(kind='bar', stacked=True, ax=axes[0], colormap='Set2', width=0.8)
axes[0].set_title('New Annotation (Aligned to GSE243013)', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Proportion', fontsize=11)
axes[0].legend(loc='upper right', fontsize=8)
axes[0].tick_params(axis='x', rotation=45)
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

old_major.plot(kind='bar', stacked=True, ax=axes[1], colormap='tab10', width=0.8)
axes[1].set_title('Original Annotation', fontsize=12, fontweight='bold')
axes[1].set_ylabel('Proportion', fontsize=11)
axes[1].legend(loc='upper right', fontsize=8)
axes[1].tick_params(axis='x', rotation=45)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '11_Old_vs_New_annotation.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [已保存] 11_Old_vs_New_annotation.png")

print("\n" + "="*70)
print("全部完成!")
print(f"结果: {OUT_DIR}")
print(f"数据: {DATA_DIR}")
print("="*70)
