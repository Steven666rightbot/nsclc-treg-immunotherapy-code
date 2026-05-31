#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE207422 基于Marker基因的重新注释 + XGBoost验证
利用GSE243013细胞类型名称中的marker基因，对GSE207422做两层注释
"""

import os
import gc
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import xgboost as xgb

warnings.filterwarnings('ignore')
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ========================== 配置 ==========================
VAL_UMI = r"D:\Research\potato\data\raw\GSE207422\GSE207422_NSCLC_scRNAseq_UMI_matrix.txt"
VAL_META = r"D:\Research\potato\data\raw\GSE207422\GSE207422_NSCLC_scRNAseq_metadata.xlsx"
VAL_ANNO = r"D:\Research\potato\data\raw\GSE207422\annotation_results\all_cells_annotation.csv"
DATA_DIR = r"D:\Research\tomato\data"
OUT_DIR = r"D:\Research\tomato\figures"
os.makedirs(OUT_DIR, exist_ok=True)

# ========================== 51类Marker定义 ==========================
# 从GSE243013的细胞类型名称 + 已知免疫细胞marker推导
MARKER_DICT = {
    # CD4 T cells
    'CD4T_Tn_CCR7':       {'positive': ['CCR7', 'LEF1', 'SELL'],           'negative': ['FOXP3', 'GZMA', 'CXCL13']},
    'CD4T_Tm_ANXA1':      {'positive': ['ANXA1', 'IL7R', 'CD27'],          'negative': ['CCR7', 'FOXP3']},
    'CD4T_Tem_GZMA':      {'positive': ['GZMA', 'GZMK', 'CCL5'],           'negative': ['CCR7', 'FOXP3']},
    'CD4T_Treg_FOXP3':    {'positive': ['FOXP3', 'IL2RA', 'CTLA4'],        'negative': ['GZMA', 'CXCL13']},
    'CD4T_Treg_CCR8':     {'positive': ['CCR8', 'FOXP3', 'IL2RA'],         'negative': ['MKI67']},
    'CD4T_Treg_MKI67':    {'positive': ['MKI67', 'FOXP3', 'TOP2A'],        'negative': ['CCR8']},
    'CD4T_Tfh_CXCL13':    {'positive': ['CXCL13', 'BCL6', 'PDCD1'],        'negative': ['FOXP3']},
    'CD4T_Th1-like_CXCL13': {'positive': ['CXCL13', 'TBX21', 'IFNG'],      'negative': ['BCL6', 'FOXP3']},
    'CD4T_Tm_XCL1':       {'positive': ['XCL1', 'XCL2', 'CD27'],           'negative': ['CCR7', 'FOXP3']},
    
    # CD8 T cells
    'CD8T_Tm_IL7R':       {'positive': ['IL7R', 'CD27', 'TCF7'],           'negative': ['GZMK', 'GZMH', 'CXCL13']},
    'CD8T_Tem_GZMK+GZMH+': {'positive': ['GZMK', 'GZMH', 'CCL5'],          'negative': ['CXCL13', 'ZNF683']},
    'CD8T_Tem_GZMK+NR4A1+': {'positive': ['GZMK', 'NR4A1', 'CCL5'],        'negative': ['GZMH', 'CXCL13']},
    'CD8T_Tex_CXCL13':    {'positive': ['CXCL13', 'PDCD1', 'HAVCR2'],      'negative': ['ZNF683', 'LAYN']},
    'CD8T_Trm_ZNF683':    {'positive': ['ZNF683', 'ITGAE', 'CD69'],         'negative': ['CXCL13', 'LAYN']},
    'CD8T_terminal_Tex_LAYN': {'positive': ['LAYN', 'CXCL13', 'HAVCR2'],   'negative': ['ZNF683']},
    'CD8T_prf_MKI67':     {'positive': ['MKI67', 'TOP2A', 'STMN1'],        'negative': ['CXCL13']},
    'CD8T_ISG15':         {'positive': ['ISG15', 'IFIT1', 'MX1'],           'negative': ['CXCL13']},
    'CD8T_MAIT_KLRB1':    {'positive': ['KLRB1', 'SLC4A10', 'ZBTB16'],     'negative': ['CXCL13']},
    'CD8T_NK-like_FGFBP2': {'positive': ['FGFBP2', 'FCGR3A', 'SPON2'],     'negative': ['CD4']},
    
    # NK cells
    'NK_CD16hi_FGFBP2':   {'positive': ['FGFBP2', 'FCGR3A', 'SPON2'],     'negative': ['CD3D', 'CD3E']},
    'NK_CD16low_GZMK':    {'positive': ['GZMK', 'NCAM1', 'XCL1'],          'negative': ['CD3D', 'FCGR3A']},
    
    # B cells
    'Bm_PDE4D':           {'positive': ['PDE4D', 'CD27', 'CD79A'],          'negative': ['CD38', 'MKI67']},
    'Bm_TNFSF9':          {'positive': ['TNFSF9', 'CD27', 'CD79A'],         'negative': ['CD38']},
    'Bm_TNF':             {'positive': ['TNF', 'CD27', 'CD79A'],            'negative': ['CD38']},
    'Bm_CD74':            {'positive': ['CD74', 'CD27', 'HLA-DRA'],        'negative': ['CD38']},
    'Bm_MT2A':            {'positive': ['MT2A', 'CD27', 'CD79A'],          'negative': ['CD38']},
    'Bm_FCRL4':           {'positive': ['FCRL4', 'CD27', 'CD79A'],          'negative': ['CD38']},
    'Bn_TCL1A':           {'positive': ['TCL1A', 'IGHM', 'CD79A'],          'negative': ['CD27', 'CD38']},
    'Plasma_cell':        {'positive': ['JCHAIN', 'MZB1', 'CD38'],          'negative': ['CD27', 'MS4A1']},
    'B_prf_MKI67':        {'positive': ['MKI67', 'TOP2A', 'CD79A'],         'negative': ['CD38', 'JCHAIN']},
    'GCB_RGS13':          {'positive': ['RGS13', 'BCL6', 'AICDA'],          'negative': ['CD27']},
    
    # Myeloid
    'Mφ_VCAN':            {'positive': ['VCAN', 'S100A8', 'S100A9'],        'negative': ['CD14']},
    'Mφ_MARCO':           {'positive': ['MARCO', 'MSR1', 'MRC1'],           'negative': ['S100A8']},
    'Mφ_FOLR2':           {'positive': ['FOLR2', 'MRC1', 'CD163'],          'negative': ['S100A8', 'CXCL10']},
    'Mφ_FCGR3A':          {'positive': ['FCGR3A', 'MS4A7', 'CD163'],        'negative': ['S100A8', 'CXCL10']},
    'Mφ_CXCL2':           {'positive': ['CXCL2', 'IL1B', 'CCL3'],           'negative': ['CXCL10']},
    'Mφ_DNAJB1':          {'positive': ['DNAJB1', 'HSPA1A', 'HSPA1B'],      'negative': ['CXCL10']},
    'Mφ_S100A10':         {'positive': ['S100A10', 'S100A11', 'ANXA2'],     'negative': ['CXCL10']},
    'Mφ_CXCL10':          {'positive': ['CXCL10', 'CXCL9', 'ISG15'],        'negative': ['S100A8']},
    'Mφ_MK67':            {'positive': ['MKI67', 'TOP2A', 'STMN1'],         'negative': ['S100A8']},
    'Mφ_ISG15':           {'positive': ['ISG15', 'IFIT1', 'MX1'],           'negative': ['S100A8']},
    'Mφ_MMP9':            {'positive': ['MMP9', 'MMP12', 'CTSD'],           'negative': ['S100A8']},
    'Mast cell':          {'positive': ['TPSAB1', 'KIT', 'CPA3'],           'negative': ['CD14', 'CD68']},
    'cDC1_CLEC9A':        {'positive': ['CLEC9A', 'XCR1', 'BATF3'],         'negative': ['CD14', 'CD163']},
    'cDC2_CD1C':          {'positive': ['CD1C', 'FCER1A', 'CLEC10A'],       'negative': ['CD14', 'CLEC9A']},
    'mDC_LAMP3':          {'positive': ['LAMP3', 'CCR7', 'BIRC3'],          'negative': ['CD14']},
    'pDC_LILRA4':         {'positive': ['LILRA4', 'CLEC4C', 'GZMB'],         'negative': ['CD14']},
    'Neu_FCGR3B':         {'positive': ['FCGR3B', 'CSF3R', 'CXCR2'],        'negative': ['CD14', 'CD68']},
    
    # Other
    'T_gdT_TRDV1':        {'positive': ['TRDV1', 'TRDC', 'TRGC1'],          'negative': ['TRAC']},
    'T_gdT_TRDV2':        {'positive': ['TRDV2', 'TRDC', 'TRGC2'],          'negative': ['TRAC']},
    'ILC3_KIT':           {'positive': ['KIT', 'RORC', 'IL22'],             'negative': ['TRAC', 'TRDC']},
}

# 大类Marker (用于第一层过滤)
LINEAGE_MARKERS = {
    'T_NK':    ['CD3D', 'CD3E', 'CD3G'],
    'CD4T':    ['CD4'],
    'CD8T':    ['CD8A', 'CD8B'],
    'B':       ['CD79A', 'CD79B', 'MS4A1'],
    'Myeloid': ['CD14', 'CD68', 'LYZ', 'CSF1R'],
    'NK':      ['NCAM1', 'NKG7', 'GNLY'],
    'Mast':    ['TPSAB1', 'KIT'],
    'DC':      ['FLT3', 'HLA-DQA1'],
    'Neu':     ['FCGR3B', 'CSF3R'],
}

def score_cell(expr_dict, marker_info):
    """计算一个细胞的marker得分"""
    pos = marker_info['positive']
    neg = marker_info.get('negative', [])
    
    pos_score = np.mean([expr_dict.get(g, 0) for g in pos if g in expr_dict])
    neg_score = np.mean([expr_dict.get(g, 0) for g in neg if g in expr_dict]) if neg else 0
    
    return pos_score - 0.3 * neg_score

# ========================== Step 1: 读取GSE207422 ==========================
print("="*70)
print("Step 1: 读取GSE207422 UMI矩阵")
print("="*70)

print("\n[1.1] 读取metadata...")
val_meta = pd.read_excel(VAL_META)
val_meta = val_meta[val_meta['Sample'].notna()].copy()
val_meta_post = val_meta[val_meta['Resource'] == 'Post-treatment surgery'].copy()
val_meta_post['response'] = val_meta_post['Pathologic Response'].map({'MPR': 1, 'pCR': 1, 'NMPR': 0})
val_meta_post = val_meta_post[val_meta_post['response'].notna()].copy()
post_samples = set(val_meta_post['Sample'].tolist())
print(f"  术后可评估样本: {sorted(post_samples)}")

print("[1.2] 读取UMI矩阵header...")
with open(VAL_UMI, 'r') as f:
    header = f.readline().strip().split('\t')
all_cells = header[1:]
cell_samples = ['_'.join(c.split('_')[:2]) for c in all_cells]
post_mask = [cs in post_samples for cs in cell_samples]
post_indices = [i for i, m in enumerate(post_mask) if m]
post_cells = [all_cells[i] for i in post_indices]
print(f"  总细胞: {len(all_cells)}, 术后细胞: {len(post_cells)}")

print("[1.3] 逐行读取，构建基因表达字典...")
# 由于细胞数多，我们不构建完整矩阵，而是对每个细胞计算大类归属
# 策略: 先读取所有基因表达，存储为 {gene: [values for post cells]}
gene_expr = {}  # gene -> numpy array of values for post cells

total_lines = 0
with open(VAL_UMI, 'r') as f:
    f.readline()
    for line in f:
        parts = line.strip().split('\t')
        gene = parts[0]
        vals = np.array(parts[1:], dtype=np.float32)
        gene_expr[gene] = vals[post_indices]
        total_lines += 1
        if total_lines % 5000 == 0:
            print(f"    {total_lines} genes read...")

print(f"  总共 {total_lines} 个基因")

# 计算每个细胞的大类得分
print("[1.4] 计算大类得分...")
n_post = len(post_cells)
lineage_scores = {}
for lineage, markers in LINEAGE_MARKERS.items():
    score = np.zeros(n_post, dtype=np.float32)
    valid = 0
    for m in markers:
        if m in gene_expr:
            score += gene_expr[m]
            valid += 1
    if valid > 0:
        score /= valid
    lineage_scores[lineage] = score

# 确定每个细胞的大类
lineage_df = pd.DataFrame(lineage_scores)
main_lineage = lineage_df.idxmax(axis=1).values
print("  大类分布:")
print(pd.Series(main_lineage).value_counts())

# ========================== Step 2: 细分类注释 ==========================
print("\n" + "="*70)
print("Step 2: 基于Marker的精细注释")
print("="*70)

print(f"[2.1] 对 {len(MARKER_DICT)} 类细胞计算得分...")
cell_type_scores = {}

for ct, markers in MARKER_DICT.items():
    score = np.zeros(n_post, dtype=np.float32)
    pos_genes = markers['positive']
    neg_genes = markers.get('negative', [])
    
    pos_valid = 0
    for g in pos_genes:
        if g in gene_expr:
            score += gene_expr[g]
            pos_valid += 1
    if pos_valid > 0:
        score /= pos_valid
    
    # 减去negative marker
    neg_score = np.zeros(n_post, dtype=np.float32)
    neg_valid = 0
    for g in neg_genes:
        if g in gene_expr:
            neg_score += gene_expr[g]
            neg_valid += 1
    if neg_valid > 0:
        neg_score /= neg_valid
    
    score -= 0.3 * neg_score
    cell_type_scores[ct] = score

score_df = pd.DataFrame(cell_type_scores)

# 根据大类限制可选类型
LINEAGE_MAP = {
    'T_NK': ['CD4T_', 'CD8T_', 'NK_', 'T_gdT_', 'ILC3_'],
    'CD4T': ['CD4T_'],
    'CD8T': ['CD8T_'],
    'B': ['Bm_', 'Bn_', 'Plasma_cell', 'B_prf_', 'GCB_'],
    'Myeloid': ['Mφ_', 'Mast cell', 'cDC', 'mDC', 'pDC', 'Neu_'],
    'NK': ['NK_'],
    'Mast': ['Mast cell'],
    'DC': ['cDC', 'mDC', 'pDC'],
    'Neu': ['Neu_'],
}

print("[2.2] 根据大类限制，分配最可能的亚型...")
predicted_types = []
max_scores = []

for i in range(n_post):
    lineage = main_lineage[i]
    allowed_prefixes = LINEAGE_MAP.get(lineage, [])
    
    # 如果是T_NK，需要进一步判断CD4/CD8/NK
    if lineage == 'T_NK':
        cd4_score = gene_expr.get('CD4', np.zeros(n_post))[i]
        cd8_score = gene_expr.get('CD8A', np.zeros(n_post))[i] + gene_expr.get('CD8B', np.zeros(n_post))[i]
        nk_score = gene_expr.get('NCAM1', np.zeros(n_post))[i] + gene_expr.get('NKG7', np.zeros(n_post))[i]
        
        if cd4_score > cd8_score and cd4_score > nk_score * 0.5:
            allowed_prefixes = ['CD4T_']
        elif cd8_score > cd4_score and cd8_score > nk_score * 0.5:
            allowed_prefixes = ['CD8T_']
        elif nk_score > cd4_score and nk_score > cd8_score:
            allowed_prefixes = ['NK_', 'T_gdT_', 'ILC3_']
    
    # 在允许的prefix中找最高分
    best_ct = None
    best_score = -1e9
    for ct in score_df.columns:
        if any(ct.startswith(p) for p in allowed_prefixes):
            s = score_df.loc[i, ct]
            if s > best_score:
                best_score = s
                best_ct = ct
    
    predicted_types.append(best_ct)
    max_scores.append(best_score)

# 处理None的情况
predicted_types = [ct if ct is not None else 'Unknown' for ct in predicted_types]

val_anno = pd.DataFrame({
    'cell': post_cells,
    'sub_cell_type': predicted_types,
    'marker_score': max_scores,
    'main_lineage': main_lineage
})

print("\n  预测分布 (Top 15):")
print(val_anno['sub_cell_type'].value_counts().head(15))

# 保存
val_anno.to_csv(os.path.join(DATA_DIR, 'GSE207422_reannotated_markers.csv'), index=False)

# ========================== Step 3: 验证 ==========================
print("\n" + "="*70)
print("Step 3: 计算比例并验证")
print("="*70)

print("[3.1] 计算细胞比例...")
val_anno['sample'] = val_anno['cell'].str.split('_').str[0] + '_' + val_anno['cell'].str.split('_').str[1]
val_counts = val_anno.groupby(['sample', 'sub_cell_type']).size().unstack(fill_value=0)
val_props = val_counts.div(val_counts.sum(axis=1), axis=0)
print(f"  比例矩阵: {val_props.shape}")

# 对齐训练集51类
train_features = pd.read_csv(os.path.join(DATA_DIR, 'cell_proportions.csv'), index_col=0).columns.tolist()
for f in train_features:
    if f not in val_props.columns:
        val_props[f] = 0.0
val_props = val_props[train_features]
val_props.to_csv(os.path.join(DATA_DIR, 'GSE207422_cell_proportions_post.csv'))

print("[3.2] 训练最终模型...")
train_props = pd.read_csv(os.path.join(DATA_DIR, 'cell_proportions.csv'), index_col=0)
train_labels = pd.read_csv(os.path.join(DATA_DIR, 'sample_labels.csv'), index_col=0)['response']
train_props = train_props.loc[train_labels.index]

scaler = StandardScaler()
X_train_s = scaler.fit_transform(train_props.values)

final_model = xgb.XGBClassifier(
    n_estimators=200, max_depth=4, learning_rate=0.05,
    subsample=0.8, colsample_bytree=0.8,
    reg_alpha=0.1, reg_lambda=1.0,
    random_state=42, eval_metric='logloss',
    use_label_encoder=False, n_jobs=4
)
final_model.fit(X_train_s, train_labels.values)
print("  模型训练完成")

print("[3.3] 验证...")
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
print(f"\n  预测详情:")
for s, p, t in zip(sample_list, y_prob, y_val):
    label = 'MPR/pCR' if t==1 else 'NMPR'
    print(f"    {s}: prob={p:.3f}, pred={int(p>0.5)}, true={label}")

# ========================== Step 4: 可视化 ==========================
print("\n" + "="*70)
print("Step 4: 可视化")
print("="*70)

# ROC
train_cv = pd.read_csv(os.path.join(DATA_DIR, 'cv_results.csv'))
mean_train_auc = train_cv['auc'].mean()
std_train_auc = train_cv['auc'].std()

fig, ax = plt.subplots(figsize=(8, 8))
fpr, tpr, _ = roc_curve(y_val, y_prob)
ax.plot(fpr, tpr, color='darkred', lw=2.5,
        label=f'GSE207422 Validation (AUC = {auc:.3f}, n={len(y_val)})')
ax.plot([0, 1], [0, 1], 'k--', lw=1.5, label='Chance')
ax.set_xlim([-0.02, 1.02])
ax.set_ylim([-0.02, 1.02])
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title(f'Validation on GSE207422 (Marker-based Reannotation)\nTraining AUC = {mean_train_auc:.3f} ± {std_train_auc:.3f}',
             fontsize=13, fontweight='bold')
ax.legend(loc='lower right', fontsize=11)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '08_Validation_ROC_GSE207422.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: 08_Validation_ROC_GSE207422.png")

# Composition
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
responders = val_meta_post[val_meta_post['response'] == 1]['Sample'].tolist()
non_responders = val_meta_post[val_meta_post['response'] == 0]['Sample'].tolist()

top10 = val_props.mean().sort_values(ascending=False).head(10).index.tolist()
resp_mean = val_props.loc[val_props.index.intersection(responders), top10].mean()
nonresp_mean = val_props.loc[val_props.index.intersection(non_responders), top10].mean()

x = np.arange(len(top10))
axes[0].bar(x - 0.35/2, resp_mean.values, 0.35, label='Responder', color='coral')
axes[0].bar(x + 0.35/2, nonresp_mean.values, 0.35, label='Non-Responder', color='steelblue')
axes[0].set_xticks(x)
axes[0].set_xticklabels(top10, rotation=45, ha='right', fontsize=9)
axes[0].set_ylabel('Mean Proportion')
axes[0].set_title('Top 10 Cell Types (GSE207422)')
axes[0].legend()
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

val_props[top10].T.plot(kind='barh', stacked=True, ax=axes[1], colormap='tab20', legend=False, width=0.7)
axes[1].set_xlabel('Proportion')
axes[1].set_title('Per-Sample Composition')
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '09_Validation_celltype_composition.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: 09_Validation_celltype_composition.png")

# Old vs New
print("[4.3] Old vs New annotation...")
old_anno = pd.read_csv(VAL_ANNO)
old_anno['sample'] = old_anno['cell'].str.split('_').str[0] + '_' + old_anno['cell'].str.split('_').str[1]
old_anno = old_anno[old_anno['sample'].isin(post_samples)]

major_map = {
    'T/NK cell': ['CD4T_', 'CD8T_', 'NK_', 'T_gdT_', 'ILC3_'],
    'B cell': ['Bm_', 'Bn_', 'Plasma_cell', 'B_prf_', 'GCB_'],
    'Myeloid cell': ['Mφ_', 'Mast cell', 'cDC', 'mDC', 'pDC', 'Neu_', 'M_']
}

def get_major(ct):
    for m, ps in major_map.items():
        for p in ps:
            if str(ct).startswith(p): return m
    return 'Other'

val_anno['major_new'] = val_anno['sub_cell_type'].apply(get_major)
new_maj = val_anno.groupby(['sample', 'major_new']).size().unstack(fill_value=0)
new_maj = new_maj.div(new_maj.sum(axis=1), axis=0)

old_keep = ['Tem/Trm cytotoxic T cells', 'Epithelial cells', 'Regulatory T cells',
            'Alveolar macrophages', 'Classical monocytes', 'Memory B cells',
            'Tcm/Naive helper T cells', 'Macrophages', 'Plasma cells', 'NK cells']
old_maj = old_anno[old_anno['predicted_labels'].isin(old_keep)].groupby(
    ['sample', 'predicted_labels']).size().unstack(fill_value=0)
old_maj = old_maj.div(old_maj.sum(axis=1), axis=0)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
new_maj.plot(kind='bar', stacked=True, ax=axes[0], colormap='Set2', width=0.8)
axes[0].set_title('New Annotation (Marker-based, Aligned to GSE243013)')
axes[0].tick_params(axis='x', rotation=45)
old_maj.plot(kind='bar', stacked=True, ax=axes[1], colormap='tab10', width=0.8)
axes[1].set_title('Original Annotation')
axes[1].tick_params(axis='x', rotation=45)
for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '11_Old_vs_New_annotation.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: 11_Old_vs_New_annotation.png")

# Marker score distribution
fig, ax = plt.subplots(figsize=(10, 6))
for ct in val_anno['sub_cell_type'].value_counts().head(10).index:
    subset = val_anno[val_anno['sub_cell_type'] == ct]['marker_score']
    ax.hist(subset, bins=30, alpha=0.5, label=ct)
ax.set_xlabel('Marker Score', fontsize=11)
ax.set_ylabel('Cell Count', fontsize=11)
ax.set_title('Marker Score Distribution (Top 10 Cell Types)', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, '12_Marker_score_distribution.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: 12_Marker_score_distribution.png")

print("\n" + "="*70)
print("ALL DONE!")
print(f"Figures: {OUT_DIR}")
print(f"Data: {DATA_DIR}")
print("="*70)
