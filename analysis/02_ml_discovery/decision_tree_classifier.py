#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Decision Tree Classifier for Immunotherapy Response Prediction
Compare with XGBoost (same 5-fold CV setup)
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier, plot_tree, export_text
from sklearn.metrics import (roc_auc_score, roc_curve, accuracy_score, 
                             precision_score, recall_score, f1_score)
import xgboost as xgb

warnings.filterwarnings('ignore')
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ===================== Paths =====================
DATA_DIR = r"D:\Research\tomato\data"
FIG_DIR = r"D:\Research\tomato\figures\decision_tree"
RESULT_DIR = r"D:\Research\tomato\results\decision_tree"

os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RESULT_DIR, exist_ok=True)

# ===================== 1. Load Data =====================
print("=" * 60)
print("Step 1: Load Data")
print("=" * 60)

cell_props = pd.read_csv(os.path.join(DATA_DIR, 'cell_proportions.csv'), index_col='sampleID')
sample_labels = pd.read_csv(os.path.join(DATA_DIR, 'sample_labels.csv'), index_col='sampleID')

# Align
common_samples = cell_props.index.intersection(sample_labels.index)
cell_props = cell_props.loc[common_samples]
sample_labels = sample_labels.loc[common_samples]

X = cell_props.values
y = sample_labels['response'].values
feature_names = cell_props.columns.tolist()

print(f"Samples: {len(common_samples)}")
print(f"Features: {len(feature_names)}")
print(f"Responders(1): {(y==1).sum()}, Non-responders(0): {(y==0).sum()}")

# ===================== 2. 5-Fold CV: Decision Tree vs XGBoost =====================
print("\n" + "=" * 60)
print("Step 2: 5-Fold Cross-Validation")
print("=" * 60)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

dt_results = []
xgb_results = []

# For ROC plotting
dt_fold_data = []
xgb_fold_data = []

for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), 1):
    print(f"\n--- Fold {fold} ---")
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    
    # Standardize
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)
    
    # --- Decision Tree (max_depth=3 for interpretability) ---
    dt = DecisionTreeClassifier(
        max_depth=3,
        min_samples_split=10,
        min_samples_leaf=5,
        criterion='gini',
        random_state=42,
        class_weight='balanced'
    )
    dt.fit(X_train_s, y_train)
    
    dt_prob = dt.predict_proba(X_test_s)[:, 1]
    dt_pred = dt.predict(X_test_s)
    
    dt_auc = roc_auc_score(y_test, dt_prob)
    dt_acc = accuracy_score(y_test, dt_pred)
    dt_prec = precision_score(y_test, dt_pred, zero_division=0)
    dt_rec = recall_score(y_test, dt_pred, zero_division=0)
    dt_f1 = f1_score(y_test, dt_pred, zero_division=0)
    
    dt_results.append({
        'fold': fold, 'auc': dt_auc, 'acc': dt_acc, 
        'prec': dt_prec, 'rec': dt_rec, 'f1': dt_f1
    })
    dt_fold_data.append({'y_test': y_test, 'y_prob': dt_prob})
    
    # --- XGBoost (same params as original) ---
    xgb_model = xgb.XGBClassifier(
        n_estimators=200, max_depth=4, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        reg_alpha=0.1, reg_lambda=1.0,
        random_state=42, eval_metric='logloss',
        use_label_encoder=False, n_jobs=4
    )
    xgb_model.fit(X_train_s, y_train)
    
    xgb_prob = xgb_model.predict_proba(X_test_s)[:, 1]
    xgb_pred = xgb_model.predict(X_test_s)
    
    xgb_auc = roc_auc_score(y_test, xgb_prob)
    xgb_acc = accuracy_score(y_test, xgb_pred)
    xgb_prec = precision_score(y_test, xgb_pred, zero_division=0)
    xgb_rec = recall_score(y_test, xgb_pred, zero_division=0)
    xgb_f1 = f1_score(y_test, xgb_pred, zero_division=0)
    
    xgb_results.append({
        'fold': fold, 'auc': xgb_auc, 'acc': xgb_acc,
        'prec': xgb_prec, 'rec': xgb_rec, 'f1': xgb_f1
    })
    xgb_fold_data.append({'y_test': y_test, 'y_prob': xgb_prob})
    
    print(f"  DT  AUC={dt_auc:.4f} ACC={dt_acc:.4f} F1={dt_f1:.4f}")
    print(f"  XGB AUC={xgb_auc:.4f} ACC={xgb_acc:.4f} F1={xgb_f1:.4f}")

# ===================== 3. Results Summary =====================
print("\n" + "=" * 60)
print("Step 3: Results Summary")
print("=" * 60)

dt_df = pd.DataFrame(dt_results)
xgb_df = pd.DataFrame(xgb_results)

print("\nDecision Tree (max_depth=3):")
print(dt_df[['fold', 'auc', 'acc', 'f1']])
print(f"  Mean AUC: {dt_df['auc'].mean():.4f} ± {dt_df['auc'].std():.4f}")
print(f"  Mean ACC: {dt_df['acc'].mean():.4f} ± {dt_df['acc'].std():.4f}")
print(f"  Mean F1:  {dt_df['f1'].mean():.4f} ± {dt_df['f1'].std():.4f}")

print("\nXGBoost:")
print(xgb_df[['fold', 'auc', 'acc', 'f1']])
print(f"  Mean AUC: {xgb_df['auc'].mean():.4f} ± {xgb_df['auc'].std():.4f}")
print(f"  Mean ACC: {xgb_df['acc'].mean():.4f} ± {xgb_df['acc'].std():.4f}")
print(f"  Mean F1:  {xgb_df['f1'].mean():.4f} ± {xgb_df['f1'].std():.4f}")

# Save
summary = pd.DataFrame({
    'metric': ['AUC', 'ACC', 'Precision', 'Recall', 'F1'],
    'DT_mean': [dt_df['auc'].mean(), dt_df['acc'].mean(), dt_df['prec'].mean(), 
                dt_df['rec'].mean(), dt_df['f1'].mean()],
    'DT_std': [dt_df['auc'].std(), dt_df['acc'].std(), dt_df['prec'].std(),
               dt_df['rec'].std(), dt_df['f1'].std()],
    'XGB_mean': [xgb_df['auc'].mean(), xgb_df['acc'].mean(), xgb_df['prec'].mean(),
                 xgb_df['rec'].mean(), xgb_df['f1'].mean()],
    'XGB_std': [xgb_df['auc'].std(), xgb_df['acc'].std(), xgb_df['prec'].std(),
                xgb_df['rec'].std(), xgb_df['f1'].std()]
})
summary.to_csv(os.path.join(RESULT_DIR, 'model_comparison.csv'), index=False)
print(f"\n[Saved] {RESULT_DIR}/model_comparison.csv")

# ===================== 4. Train Final DT on Full Data =====================
print("\n" + "=" * 60)
print("Step 4: Train Final Decision Tree (Full Data)")
print("=" * 60)

scaler_full = StandardScaler()
X_full_s = scaler_full.fit_transform(X)

dt_final = DecisionTreeClassifier(
    max_depth=3,
    min_samples_split=10,
    min_samples_leaf=5,
    criterion='gini',
    random_state=42,
    class_weight='balanced'
)
dt_final.fit(X_full_s, y)

# Feature importance
fi_df = pd.DataFrame({
    'feature': feature_names,
    'importance': dt_final.feature_importances_
}).sort_values('importance', ascending=False)
fi_df = fi_df[fi_df['importance'] > 0]
print("\nDecision Tree Feature Importance (gini):")
print(fi_df.to_string(index=False))
fi_df.to_csv(os.path.join(RESULT_DIR, 'dt_feature_importance.csv'), index=False)

# Export decision rules
tree_rules = export_text(dt_final, feature_names=feature_names, 
                         max_depth=3, spacing=3,
                         decimals=3)
print("\n" + "=" * 60)
print("Decision Rules:")
print("=" * 60)
print(tree_rules)
with open(os.path.join(RESULT_DIR, 'decision_rules.txt'), 'w') as f:
    f.write(tree_rules)

# ===================== 5. Visualizations =====================
print("\n" + "=" * 60)
print("Step 5: Generate Figures")
print("=" * 60)

# --- Figure 1: ROC Comparison ---
fig, ax = plt.subplots(figsize=(8, 8))
mean_fpr = np.linspace(0, 1, 100)

# DT ROC
dt_tprs = []
for fd in dt_fold_data:
    fpr, tpr, _ = roc_curve(fd['y_test'], fd['y_prob'])
    dt_tprs.append(np.interp(mean_fpr, fpr, tpr))
    dt_tprs[-1][0] = 0.0
dt_mean_tpr = np.mean(dt_tprs, axis=0)
dt_mean_tpr[-1] = 1.0
ax.plot(mean_fpr, dt_mean_tpr, color='steelblue', lw=2.5,
        label=f'Decision Tree (AUC = {dt_df["auc"].mean():.3f} ± {dt_df["auc"].std():.3f})')
ax.fill_between(mean_fpr, np.maximum(dt_mean_tpr - np.std(dt_tprs, axis=0), 0),
                np.minimum(dt_mean_tpr + np.std(dt_tprs, axis=0), 1), 
                color='steelblue', alpha=0.15)

# XGB ROC
xgb_tprs = []
for fd in xgb_fold_data:
    fpr, tpr, _ = roc_curve(fd['y_test'], fd['y_prob'])
    xgb_tprs.append(np.interp(mean_fpr, fpr, tpr))
    xgb_tprs[-1][0] = 0.0
xgb_mean_tpr = np.mean(xgb_tprs, axis=0)
xgb_mean_tpr[-1] = 1.0
ax.plot(mean_fpr, xgb_mean_tpr, color='darkred', lw=2.5,
        label=f'XGBoost (AUC = {xgb_df["auc"].mean():.3f} ± {xgb_df["auc"].std():.3f})')
ax.fill_between(mean_fpr, np.maximum(xgb_mean_tpr - np.std(xgb_tprs, axis=0), 0),
                np.minimum(xgb_mean_tpr + np.std(xgb_tprs, axis=0), 1), 
                color='darkred', alpha=0.15)

ax.plot([0, 1], [0, 1], 'k--', lw=1.5, label='Chance')
ax.set_xlim([-0.02, 1.02])
ax.set_ylim([-0.02, 1.02])
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC Curves: Decision Tree vs XGBoost\n(5-Fold Cross-Validation)', 
             fontsize=14, fontweight='bold')
ax.legend(loc='lower right', fontsize=11)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '01_ROC_comparison.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [Saved] 01_ROC_comparison.png")

# --- Figure 2: Performance Comparison Bar Plot ---
fig, ax = plt.subplots(figsize=(9, 6))
metrics = ['auc', 'acc', 'prec', 'rec', 'f1']
metric_labels = ['AUC', 'Accuracy', 'Precision', 'Recall', 'F1-Score']
x = np.arange(len(metrics))
width = 0.35

dt_means = [dt_df[m].mean() for m in metrics]
dt_stds = [dt_df[m].std() for m in metrics]
xgb_means = [xgb_df[m].mean() for m in metrics]
xgb_stds = [xgb_df[m].std() for m in metrics]

bars1 = ax.bar(x - width/2, dt_means, width, yerr=dt_stds, 
               label='Decision Tree', color='steelblue', alpha=0.8, capsize=4)
bars2 = ax.bar(x + width/2, xgb_means, width, yerr=xgb_stds,
               label='XGBoost', color='coral', alpha=0.8, capsize=4)

ax.set_ylabel('Score', fontsize=12)
ax.set_title('Model Performance Comparison (5-Fold CV)', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(metric_labels, fontsize=11)
ax.set_ylim([0, 1.05])
ax.legend(fontsize=11)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.3)

# Add value labels
for bar in bars1:
    height = bar.get_height()
    ax.annotate(f'{height:.3f}', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)
for bar in bars2:
    height = bar.get_height()
    ax.annotate(f'{height:.3f}', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '02_performance_comparison.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [Saved] 02_performance_comparison.png")

# --- Figure 3: Decision Tree Visualization ---
fig, ax = plt.subplots(figsize=(24, 14))
plot_tree(dt_final, 
          feature_names=feature_names,
          class_names=['Non-responder', 'Responder'],
          filled=True,
          rounded=True,
          fontsize=10,
          max_depth=3,
          impurity=False,
          proportion=True,
          ax=ax)
ax.set_title('Decision Tree for Immunotherapy Response Prediction\n(max_depth=3, trained on full data)', 
             fontsize=16, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '03_decision_tree_structure.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [Saved] 03_decision_tree_structure.png")

# --- Figure 4: Feature Importance (DT vs XGBoost) ---
# Re-train XGBoost on full data for fair comparison
xgb_full = xgb.XGBClassifier(
    n_estimators=200, max_depth=4, learning_rate=0.05,
    subsample=0.8, colsample_bytree=0.8,
    reg_alpha=0.1, reg_lambda=1.0,
    random_state=42, eval_metric='logloss',
    use_label_encoder=False, n_jobs=4
)
xgb_full.fit(X_full_s, y)

xgb_fi = pd.DataFrame({
    'feature': feature_names,
    'importance': xgb_full.feature_importances_
}).sort_values('importance', ascending=False)

# Merge top 15 from each
top_dt = fi_df.head(15)['feature'].tolist()
top_xgb = xgb_fi.head(15)['feature'].tolist()
all_top = list(dict.fromkeys(top_dt + top_xgb))  # unique, preserve order

comp_fi = pd.DataFrame({'feature': all_top})
comp_fi = comp_fi.merge(fi_df, on='feature', how='left').rename(columns={'importance': 'DT'})
comp_fi = comp_fi.merge(xgb_fi, on='feature', how='left').rename(columns={'importance': 'XGB'})
comp_fi = comp_fi.fillna(0)

fig, ax = plt.subplots(figsize=(10, 8))
y_pos = np.arange(len(comp_fi))
height = 0.35
ax.barh(y_pos - height/2, comp_fi['DT'], height, label='Decision Tree', color='steelblue', alpha=0.8)
ax.barh(y_pos + height/2, comp_fi['XGB'], height, label='XGBoost', color='coral', alpha=0.8)
ax.set_yticks(y_pos)
ax.set_yticklabels(comp_fi['feature'], fontsize=10)
ax.set_xlabel('Feature Importance', fontsize=12)
ax.set_title('Feature Importance: Decision Tree vs XGBoost', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '04_feature_importance_comparison.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [Saved] 04_feature_importance_comparison.png")

# ===================== 6. Save Final Model =====================
import joblib
joblib.dump({'model': dt_final, 'scaler': scaler_full, 'feature_names': feature_names},
            os.path.join(RESULT_DIR, 'dt_final_model.pkl'))
print(f"\n[Saved] Final model: {RESULT_DIR}/dt_final_model.pkl")

print("\n" + "=" * 60)
print("Analysis Complete!")
print(f"Figures: {FIG_DIR}")
print(f"Results: {RESULT_DIR}")
print("=" * 60)
