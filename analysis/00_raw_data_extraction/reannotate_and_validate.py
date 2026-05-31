#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE207422 重新注释 + XGBoost验证
基于GSE243013的51类sub_cell_type参考谱对GSE207422做映射注释
"""

import os
import gc
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial.distance import cosine
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import xgboost as xgb
import pickle

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
CODE_DIR = r"D:\Research\tomato\code"

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

N_CELLS_PER_TYPE = 2000  # 每类抽样的细胞数
MIN_CELLS_PER_TYPE = 100  # 最少细胞数

print("="*70)
print("Step 1: 构建GSE243013参考表达谱")
print("="*70)

# ---- 1.1 读取metadata，建立cellID → sub_cell_type映射 ----
print("\n[1.1] 读取GSE243013 metadata...")
meta = pd.read_csv(REF_META, low_memory=False)
meta = meta[meta['pathological_response'] != 'unknowm'].copy()
cell2type = dict(zip(meta['cellID'], meta['sub_cell_type']))
print(f"  有效细胞数: {len(cell2type)}")

# ---- 1.2 读取barcodes，建立index → cellID映射 ----
print("\n[1.2] 读取barcodes...")
barcodes = pd.read_csv(REF_BARCODES, header=None)[0].tolist()
# 第一行是header "barcode"
barcodes = barcodes[1:]
print(f"  Barcodes数: {len(barcodes)}")
assert len(barcodes) == 1254749, f"Expected 1254749, got {len(barcodes)}"

# ---- 1.3 读取features，建立index → gene映射 ----
print("\n[1.3] 读取features...")
features = pd.read_csv(REF_FEATURES, header=None)[0].tolist()
features = features[1:]  # skip header
print(f"  Genes数: {len(features)}")
assert len(features) == 31831, f"Expected 31831, got {len(features)}"

# ---- 1.4 分层抽样 ----
print("\n[1.4] 分层抽样...")
np.random.seed(42)
sampled_indices = set()
sampled_cell_type = {}  # index → sub_cell_type

type_counts = meta['sub_cell_type'].value_counts()
for cell_type, count in type_counts.items():
    cells = meta[meta['sub_cell_type'] == cell_type]['cellID'].values
    n_sample = min(N_CELLS_PER_TYPE, len(cells))
    if n_sample < MIN_CELLS_PER_TYPE:
        print(f"  WARN: {cell_type} only has {len(cells)} cells")
        continue
    selected = np.random.choice(cells, size=n_sample, replace=False)
    for cid in selected:
        # 找到barcode index (1-based)
        idx = barcodes.index(cid) + 1  # mtx uses 1-based indexing
        sampled_indices.add(idx)
        sampled_cell_type[idx] = cell_type

print(f"  抽样细胞数: {len(sampled_indices)}")

# ---- 1.5 扫描mtx文件，提取抽样细胞的表达 ----
print("\n[1.5] 扫描mtx文件 (约22GB, 20亿行)...")
print("  这可能需要10-20分钟，请耐心等待...")

# 按sub_cell_type聚合: {cell_type: {gene_idx: [values]}}
from collections import defaultdict
cell_type_data = defaultdict(lambda: defaultdict(list))

current_row = None
current_genes = []
current_vals = []

processed_cells = 0
kept_cells = 0
line_num = 0

with open(REF_MTX, 'r') as f:
    # 跳过头部
    for line in f:
        line_num += 1
        if line.startswith('%'):
            continue
        parts = line.strip().split()
        if len(parts) == 3:
            n_rows, n_cols, n_nnz = int(parts[0]), int(parts[1]), int(parts[2])
            print(f"  mtx: {n_rows} rows, {n_cols} cols, {n_nnz} nonzeros")
            break
    
    # 读取数据
    for line in f:
        line_num += 1
        parts = line.strip().split()
        if len(parts) != 3:
            continue
        row_idx = int(parts[0])
        col_idx = int(parts[1])
        val = float(parts[2])
        
        if row_idx != current_row:
            # 保存上一个cell的数据
            if current_row is not None and current_row in sampled_indices:
                ct = sampled_cell_type[current_row]
                for g, v in zip(current_genes, current_vals):
                    cell_type_data[ct][g].append(v)
                kept_cells += 1
            
            current_row = row_idx
            current_genes = [col_idx]
            current_vals = [val]
            processed_cells += 1
            
            if processed_cells % 50000 == 0:
                print(f"    已处理 {processed_cells} cells, 保留 {kept_cells}")
        else:
            current_genes.append(col_idx)
            current_vals.append(val)
    
    # 处理最后一个cell
    if current_row is not None and current_row in sampled_indices:
        ct = sampled_cell_type[current_row]
        for g, v in zip(current_genes, current_vals):
            cell_type_data[ct][g].append(v)
        kept_cells += 1

print(f"  完成! 处理cells: {processed_cells}, 保留cells: {kept_cells}")

# ---- 1.6 计算每个sub_cell_type的平均表达 (pseudo-bulk) ----
print("\n[1.6] 计算pseudo-bulk参考谱...")
n_genes = len(features)
cell_types = sorted(cell_type_data.keys())
reference_matrix = np.zeros((len(cell_types), n_genes))

for i, ct in enumerate(cell_types):
    gene_data = cell_type_data[ct]
    # 对每个基因，计算平均表达
    for g_idx, vals in gene_data.items():
        # 注意: vals是多个cell的该基因表达值
        # 需要除以该cell_type的cell数（不是非零cell数）
        n_cells_in_type = len([k for k, v in sampled_cell_type.items() if v == ct])
        reference_matrix[i, g_idx - 1] = sum(vals) / n_cells_in_type

# 转换为log1p
reference_matrix = np.log1p(reference_matrix)

print(f"  参考谱: {reference_matrix.shape}")
print(f"  细胞类型: {cell_types}")

# 清理内存
del cell_type_data, meta, cell2type
gc.collect()

# 保存参考谱
ref_df = pd.DataFrame(reference_matrix, index=cell_types, columns=features)
ref_df.to_csv(os.path.join(DATA_DIR, 'GSE243013_reference_pseudobulk.csv'))
print(f"  参考谱已保存")

# ========================== Step 2: 读取GSE207422 ==========================
print("\n" + "="*70)
print("Step 2: 读取并处理GSE207422")
print("="*70)

# ---- 2.1 读取UMI矩阵 ----
print("\n[2.1] 读取GSE207422 UMI矩阵 (约2GB)...")
val_df = pd.read_csv(VAL_UMI, sep='\t', index_col=0)
print(f"  矩阵: {val_df.shape}")
val_genes = val_df.index.tolist()
val_cells = val_df.columns.tolist()

# 转换为numpy数组
val_X = val_df.values.astype(np.float32)
print(f"  内存占用: {val_X.nbytes / 1e9:.2f} GB")

# ---- 2.2 基因名对齐 ----
print("\n[2.2] 基因对齐...")
common_genes = list(set(features) & set(val_genes))
print(f"  交集基因: {len(common_genes)}")

ref_df_aligned = ref_df[common_genes]
val_df_aligned = val_df.loc[common_genes]

ref_X = ref_df_aligned.values
val_X_aligned = val_df_aligned.values.T  # cells x genes
print(f"  参考谱对齐后: {ref_X.shape}")
print(f"  验证集对齐后: {val_X_aligned.shape}")

# ---- 2.3 Log1p标准化 ----
print("\n[2.3] 标准化...")
# 对验证集每个cell做CPM -> log1p
val_cpm = val_X_aligned / val_X_aligned.sum(axis=1, keepdims=True) * 1e6
val_cpm = np.log1p(val_cpm)

# 参考谱已经是log1p的

# ---- 2.4 余弦相似度映射 ----
print("\n[2.4] 余弦相似度映射注释...")
from sklearn.metrics.pairwise import cosine_similarity

# 计算每个验证集细胞与51个参考谱的相似度
similarities = cosine_similarity(val_cpm, ref_X)
predicted_types_idx = np.argmax(similarities, axis=1)
predicted_types = [cell_types[i] for i in predicted_types_idx]
max_similarities = np.max(similarities, axis=1)

print(f"  平均最大相似度: {max_similarities.mean():.4f}")
print(f"  中位最大相似度: {np.median(max_similarities):.4f}")

# 过滤低置信度注释（可选）
# low_conf_thresh = np.percentile(max_similarities, 10)
# print(f"  低置信度阈值(10th percentile): {low_conf_thresh:.4f}")

# 构建注释DataFrame
val_anno = pd.DataFrame({
    'cell': val_cells,
    'sub_cell_type': predicted_types,
    'cosine_similarity': max_similarities
})

print("\n  预测细胞类型分布:")
print(val_anno['sub_cell_type'].value_counts().head(20))

# 保存注释结果
val_anno.to_csv(os.path.join(DATA_DIR, 'GSE207422_reannotated.csv'), index=False)

# ========================== Step 3: 计算细胞比例并验证 ==========================
print("\n" + "="*70)
print("Step 3: 计算细胞比例并验证模型")
print("="*70)

# ---- 3.1 读取验证集metadata ----
print("\n[3.1] 读取验证集metadata...")
val_meta = pd.read_excel(VAL_META)
# 去掉注释行
val_meta = val_meta[val_meta['Sample'].notna()].copy()
print(val_meta[['Sample', 'Patient', 'Resource', 'Pathologic Response']])

# 只保留Post-treatment surgery样本
val_meta_post = val_meta[val_meta['Resource'] == 'Post-treatment surgery'].copy()
# 映射响应标签
val_meta_post['response'] = val_meta_post['Pathologic Response'].map({
    'MPR': 1, 'pCR': 1, 'NMPR': 0
})
val_meta_post = val_meta_post[val_meta_post['response'].notna()].copy()
print(f"\n  术后可评估样本: {len(val_meta_post)}")
print(val_meta_post[['Sample', 'Pathologic Response', 'response']])

# ---- 3.2 计算细胞比例 ----
print("\n[3.2] 计算细胞比例...")
# 从cell名推断sample
val_anno['sample'] = val_anno['cell'].str.split('_').str[0] + '_' + val_anno['cell'].str.split('_').str[1]

# 计算每个sample的细胞比例
val_counts = val_anno.groupby(['sample', 'sub_cell_type']).size().unstack(fill_value=0)
val_props = val_counts.div(val_counts.sum(axis=1), axis=0)

# 只保留术后样本
val_props_post = val_props.loc[val_props.index.intersection(val_meta_post['Sample'])]
print(f"  术后样本细胞比例矩阵: {val_props_post.shape}")

# 确保包含所有51类细胞（对齐到训练集的51类）
train_features = pd.read_csv(os.path.join(DATA_DIR, 'cell_proportions.csv'), index_col=0).columns.tolist()
for f in train_features:
    if f not in val_props_post.columns:
        val_props_post[f] = 0.0
val_props_post = val_props_post[train_features]

# 保存
val_props_post.to_csv(os.path.join(DATA_DIR, 'GSE207422_cell_proportions_post.csv'))

# ---- 3.3 加载之前训练的模型 ----
print("\n[3.3] 加载训练好的模型...")

# 重新训练一个全量模型（用全部GSE243013数据）
train_props = pd.read_csv(os.path.join(DATA_DIR, 'cell_proportions.csv'), index_col=0)
train_labels = pd.read_csv(os.path.join(DATA_DIR, 'sample_labels.csv'), index_col=0)['response']

# 对齐
train_props = train_props.loc[train_labels.index]
X_train = train_props.values
y_train = train_labels.values

# 标准化
scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)

# 训练最终模型
final_model = xgb.XGBClassifier(
    n_estimators=200,
    max_depth=4,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    reg_alpha=0.1,
    reg_lambda=1.0,
    random_state=42,
    eval_metric='logloss',
    use_label_encoder=False,
    n_jobs=4
)
final_model.fit(X_train_s, y_train)
print("  最终模型训练完成")

# ---- 3.4 验证集预测 ----
print("\n[3.4] 验证集预测...")
X_val = val_props_post.values
y_val = val_meta_post.set_index('Sample').loc[val_props_post.index, 'response'].values

X_val_s = scaler.transform(X_val)
y_prob = final_model.predict_proba(X_val_s)[:, 1]
y_pred = final_model.predict(X_val_s)

# 评估
auc = roc_auc_score(y_val, y_prob)
acc = accuracy_score(y_val, y_pred)
print(f"\n  验证集 AUC: {auc:.4f}")
print(f"  验证集 ACC: {acc:.4f}")
print(f"\n  样本预测概率:")
for s, p, t in zip(val_props_post.index, y_prob, y_val):
    print(f"    {s}: prob={p:.3f}, true={int(t)}")

# ========================== Step 4: 可视化 ==========================
print("\n" + "="*70)
print("Step 4: 生成可视化")
print("="*70)

# ---- 4.1 验证集ROC曲线 ----
fig, ax = plt.subplots(figsize=(8, 8))

# 训练集五折平均ROC（从之前的结果读取）
train_cv = pd.read_csv(os.path.join(DATA_DIR, 'cv_results.csv'))
mean_auc_train = train_cv['auc'].mean()
std_auc_train = train_cv['auc'].std()

# 绘制验证集ROC
fpr, tpr, thresholds = roc_curve(y_val, y_prob)
ax.plot(fpr, tpr, color='darkred', lw=2.5,
        label=f'GSE207422 Validation (AUC = {auc:.3f}, n={len(y_val)})')

# 参考线
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

# ---- 4.2 验证集细胞比例堆叠柱状图 ----
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# 按响应分组的颜色
responders = val_meta_post[val_meta_post['response'] == 1]['Sample'].tolist()
non_responders = val_meta_post[val_meta_post['response'] == 0]['Sample'].tolist()

# Top 10细胞类型的比例比较
top10 = val_props_post.mean().sort_values(ascending=False).head(10).index.tolist()

resp_mean = val_props_post.loc[val_props_post.index.intersection(responders), top10].mean()
nonresp_mean = val_props_post.loc[val_props_post.index.intersection(non_responders), top10].mean()

x = np.arange(len(top10))
width = 0.35
axes[0].bar(x - width/2, resp_mean.values, width, label='Responder (MPR/pCR)', color='coral')
axes[0].bar(x + width/2, nonresp_mean.values, width, label='Non-Responder (NMPR)', color='steelblue')
axes[0].set_xticks(x)
axes[0].set_xticklabels(top10, rotation=45, ha='right')
axes[0].set_ylabel('Mean Proportion', fontsize=11)
axes[0].set_title('Top 10 Cell Type Proportions\n(GSE207422 Validation)', fontsize=12, fontweight='bold')
axes[0].legend()
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

# 每个样本的细胞比例堆叠图
val_props_plot = val_props_post[top10].T
val_props_plot.plot(kind='barh', stacked=True, ax=axes[1], 
                    colormap='tab20', legend=False, width=0.7)
axes[1].set_xlabel('Proportion', fontsize=11)
axes[1].set_title('Cell Type Composition per Sample\n(GSE207422)', fontsize=12, fontweight='bold')
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '09_Validation_celltype_composition.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [已保存] 09_Validation_celltype_composition.png")

# ---- 4.3 相似度分布图 ----
fig, ax = plt.subplots(figsize=(10, 6))
val_anno['similarity_bin'] = pd.cut(val_anno['cosine_similarity'], bins=20)
sim_by_type = val_anno.groupby(['sub_cell_type', 'similarity_bin']).size().unstack(fill_value=0)
# 只画top 15类型
top15_types = val_anno['sub_cell_type'].value_counts().head(15).index
sim_top15 = sim_by_type.loc[top15_types]
sim_top15.div(sim_top15.sum(axis=1), axis=0).plot(kind='barh', stacked=True, 
                                                   colormap='viridis', ax=ax, legend=False)
ax.set_xlabel('Proportion', fontsize=11)
ax.set_title('Cosine Similarity Distribution by Cell Type\n(GSE207422 Reannotation)', 
             fontsize=12, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '10_Reannotation_similarity_distribution.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [已保存] 10_Reannotation_similarity_distribution.png")

# ---- 4.4 对比原有注释和新注释 ----
print("\n[4.4] 注释对比...")
old_anno = pd.read_csv(VAL_ANNO)
old_anno['sample'] = old_anno['cell'].str.split('_').str[0] + '_' + old_anno['cell'].str.split('_').str[1]

# 主要类型映射对比
major_map = {
    'T/NK cell': ['CD4T_', 'CD8T_', 'NK_', 'T_gdT_', 'ILC3_', 'CD4T_'],
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
old_anno['major_old'] = old_anno['predicted_labels']

# 对比每个样本的主要细胞类型比例
new_major_props = val_anno.groupby(['sample', 'major_new']).size().unstack(fill_value=0)
new_major_props = new_major_props.div(new_major_props.sum(axis=1), axis=0)

old_major_keep = ['Tem/Trm cytotoxic T cells', 'Epithelial cells', 'Regulatory T cells',
                  'Alveolar macrophages', 'Classical monocytes', 'Memory B cells',
                  'Tcm/Naive helper T cells', 'Macrophages', 'Plasma cells', 'NK cells']
old_major_props = old_anno[old_anno['predicted_labels'].isin(old_major_keep)].groupby(
    ['sample', 'predicted_labels']).size().unstack(fill_value=0)
old_major_props = old_major_props.div(old_major_props.sum(axis=1), axis=0)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
new_major_props.plot(kind='bar', stacked=True, ax=axes[0], colormap='Set2', width=0.8)
axes[0].set_title('New Annotation (Aligned to GSE243013)', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Proportion', fontsize=11)
axes[0].legend(loc='upper right', fontsize=8)
axes[0].tick_params(axis='x', rotation=45)
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

old_major_props.plot(kind='bar', stacked=True, ax=axes[1], colormap='tab10', width=0.8)
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
print(f"结果保存在: {OUT_DIR}")
print(f"数据保存在: {DATA_DIR}")
print("="*70)
