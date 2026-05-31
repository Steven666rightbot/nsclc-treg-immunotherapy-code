#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE243013 - XGBoost + SHAP 细胞比例分析
探究51类细胞亚型比例对免疫治疗响应（MPR+pCR vs non-MPR）的预测价值
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, precision_score, recall_score, f1_score
import xgboost as xgb
import shap
from matplotlib.colors import LinearSegmentedColormap

warnings.filterwarnings('ignore')

# ── Unified green-white-purple color scheme for Fig1 ──
COLOR_HIGH = '#8bc98b'   # soft green
COLOR_LOW  = '#b08ebf'   # soft purple
CMAP_GWP = LinearSegmentedColormap.from_list('gwp',
    [COLOR_LOW, '#d4c2df', '#ffffff', '#c2e4c2', COLOR_HIGH])
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ===================== 1. 数据读取与预处理 =====================
print("=" * 60)
print("Step 1: 读取数据并预处理")
print("=" * 60)

meta_path = r"D:\Research\土豆\数据\raw\GSE243013\GSE243013_NSCLC_immune_scRNA_metadata.csv"
out_dir = r"D:\Research\土豆\数据\raw\GSE243013\XGBoost_SHAP_results"
os.makedirs(out_dir, exist_ok=True)

meta = pd.read_csv(meta_path, low_memory=False)
print(f"原始metadata: {meta.shape}")

# 过滤掉unknown响应的样本
meta = meta[meta['pathological_response'] != 'unknowm'].copy()
print(f"过滤unknown后: {meta.shape}")

# 定义响应标签：MPR + pCR = 1 (响应者), non-MPR = 0 (非响应者)
meta['response'] = meta['pathological_response'].map({'MPR': 1, 'pCR': 1, 'non-MPR': 0})

# ===================== 2. 计算每个样本的细胞比例 =====================
print("\n" + "=" * 60)
print("Step 2: 计算样本水平的51类细胞比例")
print("=" * 60)

# 统计每个样本中各细胞亚型的数量
cell_counts = meta.groupby(['sampleID', 'sub_cell_type']).size().unstack(fill_value=0)
print(f"细胞计数矩阵: {cell_counts.shape}")

# 计算比例
cell_props = cell_counts.div(cell_counts.sum(axis=1), axis=0)
print(f"细胞比例矩阵: {cell_props.shape}")
print(f"细胞类型数: {cell_props.shape[1]}")

# 获取每个样本的响应标签
sample_labels = meta.groupby('sampleID')['response'].first()

# 对齐数据
common_samples = cell_props.index.intersection(sample_labels.index)
cell_props = cell_props.loc[common_samples]
sample_labels = sample_labels.loc[common_samples]

print(f"最终样本数: {len(common_samples)}")
print(f"响应者(1): {(sample_labels == 1).sum()}, 非响应者(0): {(sample_labels == 0).sum()}")

X = cell_props.values
y = sample_labels.values
feature_names = cell_props.columns.tolist()

# ===================== 3. 五折交叉验证 + XGBoost =====================
print("\n" + "=" * 60)
print("Step 3: 五折交叉验证 + XGBoost")
print("=" * 60)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

fold_results = []
all_shap_values = []
all_test_indices = []
models = []

for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), 1):
    print(f"\n--- Fold {fold} ---")
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    
    # 标准化
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)
    
    # XGBoost
    model = xgb.XGBClassifier(
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
    model.fit(X_train_s, y_train)
    
    # 预测
    y_prob = model.predict_proba(X_test_s)[:, 1]
    y_pred = model.predict(X_test_s)
    
    # 评估
    auc = roc_auc_score(y_test, y_prob)
    acc = accuracy_score(y_test, y_pred)
    prec = precision_score(y_test, y_pred, zero_division=0)
    rec = recall_score(y_test, y_pred, zero_division=0)
    f1 = f1_score(y_test, y_pred, zero_division=0)
    
    print(f"  AUC: {auc:.4f}, ACC: {acc:.4f}, Prec: {prec:.4f}, Rec: {rec:.4f}, F1: {f1:.4f}")
    
    fold_results.append({
        'fold': fold,
        'auc': auc, 'acc': acc, 'prec': prec, 'rec': rec, 'f1': f1,
        'y_test': y_test, 'y_prob': y_prob
    })
    
    # SHAP
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test_s)
    all_shap_values.append(shap_values)
    all_test_indices.append(test_idx)
    models.append((model, scaler))

# 汇总结果
results_df = pd.DataFrame([{k: v for k, v in r.items() if k not in ['y_test', 'y_prob']} for r in fold_results])
print("\n" + "=" * 60)
print("五折交叉验证结果汇总")
print("=" * 60)
print(results_df[['fold', 'auc', 'acc', 'prec', 'rec', 'f1']])
print(f"\n平均 AUC: {results_df['auc'].mean():.4f} ± {results_df['auc'].std():.4f}")
print(f"平均 ACC: {results_df['acc'].mean():.4f} ± {results_df['acc'].std():.4f}")
print(f"平均 F1:  {results_df['f1'].mean():.4f} ± {results_df['f1'].std():.4f}")

# 保存结果
results_df.to_csv(os.path.join(out_dir, 'cv_results.csv'), index=False)

# ===================== 4. SHAP 分析 =====================
print("\n" + "=" * 60)
print("Step 4: SHAP 分析")
print("=" * 60)

# 合并所有fold的SHAP值
shap_matrix = np.vstack(all_shap_values)
# 构建对应的X矩阵（标准化后的）
X_all_test = np.vstack([models[i][1].transform(X[all_test_indices[i]]) for i in range(5)])

print(f"SHAP矩阵: {shap_matrix.shape}")
print(f"X矩阵: {X_all_test.shape}")

# 计算平均绝对SHAP值
mean_abs_shap = np.abs(shap_matrix).mean(axis=0)
shap_importance = pd.DataFrame({
    'feature': feature_names,
    'mean_abs_shap': mean_abs_shap
}).sort_values('mean_abs_shap', ascending=False)

print("\nTop 15 重要细胞类型 (按平均|SHAP|):")
print(shap_importance.head(15).to_string(index=False))
shap_importance.to_csv(os.path.join(out_dir, 'shap_importance.csv'), index=False)

# ===================== 5. 可视化 =====================
print("\n" + "=" * 60)
print("Step 5: 生成可视化图表")
print("=" * 60)

# --- Figure 1: SHAP Summary (Beeswarm) --- [Fig1g]
fig, ax = plt.subplots(figsize=(12, 14))
shap.summary_plot(
    shap_matrix,
    X_all_test,
    feature_names=feature_names,
    max_display=20,
    show=False,
    plot_size=None,
    color=CMAP_GWP
)
plt.title("SHAP Summary Plot (Top 20 Cell Types)", fontsize=14, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'fig1g_SHAP_summary_beeswarm.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [已保存] 01_SHAP_summary_beeswarm.png")

# --- Figure 2: SHAP Bar Plot (Mean |SHAP|) --- [Fig1e]
fig, ax = plt.subplots(figsize=(10, 14))
top_n = 20
top_features = shap_importance.head(top_n)
colors = CMAP_GWP(np.linspace(0, 1, top_n))
bars = ax.barh(range(top_n), top_features['mean_abs_shap'].values[::-1], color=colors[::-1])
ax.set_yticks(range(top_n))
ax.set_yticklabels(top_features['feature'].values[::-1], fontsize=11)
ax.set_xlabel('Mean |SHAP value|', fontsize=12)
ax.set_title('SHAP Feature Importance (Top 20 Cell Types)', fontsize=14, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'fig1e_SHAP_bar_importance.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [已保存] 02_SHAP_bar_importance.png")

# --- Figure 3: ROC Curves (5-Fold) --- [Fig1d]
fig, ax = plt.subplots(figsize=(8, 8))
# Green-purple gradient for 5 folds
fold_colors = ['#c2a5cf', '#d4c2df', '#e8e8e8', '#c2e4c2', '#a6dba0']
mean_fpr = np.linspace(0, 1, 100)
tprs = []

for i, res in enumerate(fold_results):
    fpr, tpr, _ = roc_curve(res['y_test'], res['y_prob'])
    tprs.append(np.interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    ax.plot(fpr, tpr, color=fold_colors[i], lw=1.5, alpha=0.7,
            label=f'Fold {i+1} (AUC = {res["auc"]:.3f})')

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = results_df['auc'].mean()
std_auc = results_df['auc'].std()

ax.plot(mean_fpr, mean_tpr, color='#5a9e5a', lw=2.5,
        label=f'Mean (AUC = {mean_auc:.3f} ± {std_auc:.3f})')
ax.fill_between(mean_fpr, np.maximum(mean_tpr - np.std(tprs, axis=0), 0),
                np.minimum(mean_tpr + np.std(tprs, axis=0), 1), color='#c2e4c2', alpha=0.35)

ax.plot([0, 1], [0, 1], 'k--', lw=1.5, label='Chance')
ax.set_xlim([-0.02, 1.02])
ax.set_ylim([-0.02, 1.02])
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC Curves (5-Fold Cross-Validation)', fontsize=14, fontweight='bold')
ax.legend(loc='lower right', fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'fig1d_ROC_curves_5fold.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [已保存] 03_ROC_curves_5fold.png")

# --- Figure 4: Model Performance Boxplot ---
fig, ax = plt.subplots(figsize=(8, 6))
metrics = ['auc', 'acc', 'prec', 'rec', 'f1']
metric_labels = ['AUC', 'Accuracy', 'Precision', 'Recall', 'F1-Score']
data_for_box = [results_df[m].values for m in metrics]
bp = ax.boxplot(data_for_box, labels=metric_labels, patch_artist=True,
                boxprops=dict(facecolor='steelblue', alpha=0.7),
                medianprops=dict(color='darkred', linewidth=2),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                flierprops=dict(marker='o', color='red', alpha=0.5))
ax.set_ylabel('Score', fontsize=12)
ax.set_title('Model Performance Across 5 Folds', fontsize=14, fontweight='bold')
ax.set_ylim([0, 1.05])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.3)
# 添加均值点
for i, m in enumerate(metrics):
    ax.scatter(i+1, results_df[m].mean(), color='gold', s=100, zorder=5, edgecolors='black', linewidth=1)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, '04_performance_boxplot.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [已保存] 04_performance_boxplot.png")

# --- Figure 5: Top 10 Cell Type SHAP Waterfall-style Heatmap ---
fig, ax = plt.subplots(figsize=(12, 8))
top10 = shap_importance.head(10)['feature'].tolist()
top10_idx = [feature_names.index(f) for f in top10]

# 计算每个细胞类型在响应组和非响应组中的平均比例
plot_df = cell_props[top10].copy()
plot_df['response'] = sample_labels

# 计算组间差异
mean_resp = plot_df[plot_df['response'] == 1][top10].mean()
mean_nonresp = plot_df[plot_df['response'] == 0][top10].mean()
diff = (mean_resp - mean_nonresp).sort_values(ascending=True)

bars = ax.barh(range(len(diff)), diff.values, 
               color=['steelblue' if v < 0 else 'coral' for v in diff.values])
ax.set_yticks(range(len(diff)))
ax.set_yticklabels(diff.index, fontsize=11)
ax.axvline(x=0, color='black', linewidth=0.8)
ax.set_xlabel('Mean Proportion Difference (Responder - Non-Responder)', fontsize=12)
ax.set_title('Cell Proportion Difference Between Responders and Non-Responders\n(Top 10 SHAP-important Cell Types)', 
             fontsize=13, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, '05_celltype_proportion_diff.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [已保存] 05_celltype_proportion_diff.png")

# --- Figure 6: Correlation Heatmap of Top 15 Cell Types ---
fig, ax = plt.subplots(figsize=(14, 12))
top15 = shap_importance.head(15)['feature'].tolist()
corr = cell_props[top15].corr()
mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
sns.heatmap(corr, mask=mask, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
            square=True, linewidths=0.5, cbar_kws={"shrink": 0.8}, ax=ax,
            annot_kws={"size": 8})
ax.set_title('Correlation Matrix of Top 15 SHAP-important Cell Types', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, '06_correlation_heatmap.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [已保存] 06_correlation_heatmap.png")

# --- Figure 7: SHAP Dependence Plots for Top 5 Features ---
top5 = shap_importance.head(5)['feature'].tolist()
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes = axes.flatten()

for idx, feat in enumerate(top5):
    fi = feature_names.index(feat)
    shap.dependence_plot(
        fi, shap_matrix, X_all_test, feature_names=feature_names,
        interaction_index=None, show=False, ax=axes[idx]
    )
    axes[idx].set_title(f'SHAP Dependence: {feat}', fontsize=11, fontweight='bold')

axes[5].axis('off')
plt.suptitle('SHAP Dependence Plots (Top 5 Cell Types)', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, '07_SHAP_dependence_top5.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [已保存] 07_SHAP_dependence_top5.png")

# ===================== 6. 保存额外数据 =====================
print("\n" + "=" * 60)
print("Step 6: 保存分析数据")
print("=" * 60)

# 保存细胞比例矩阵
cell_props.to_csv(os.path.join(out_dir, 'cell_proportions.csv'))
print("  [已保存] cell_proportions.csv")

# 保存样本标签
sample_labels.to_frame('response').to_csv(os.path.join(out_dir, 'sample_labels.csv'))
print("  [已保存] sample_labels.csv")

# 保存Top 20特征详细信息
top20_detail = shap_importance.head(20).copy()
top20_detail['responder_mean'] = [cell_props.loc[sample_labels == 1, f].mean() for f in top20_detail['feature']]
top20_detail['nonresponder_mean'] = [cell_props.loc[sample_labels == 0, f].mean() for f in top20_detail['feature']]
top20_detail['proportion_diff'] = top20_detail['responder_mean'] - top20_detail['nonresponder_mean']
top20_detail.to_csv(os.path.join(out_dir, 'top20_features_detail.csv'), index=False)
print("  [已保存] top20_features_detail.csv")

print("\n" + "=" * 60)
print("分析完成！所有结果保存在:")
print(f"  {out_dir}")
print("=" * 60)
