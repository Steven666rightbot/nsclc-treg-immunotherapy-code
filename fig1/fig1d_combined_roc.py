#!/usr/bin/env python3
"""
Figure 1D (Revised) — Combined ROC: 5-fold CV + GSE207422 External Validation

Plots:
  - 5 individual fold ROC curves (GSE243013, stratified 5-fold CV)
  - Mean ROC ± 1 SD
  - GSE207422 external validation (n=12, primary AUC=0.781, seed sensitivity 0.719–0.844)
  - Chance diagonal

Output: fig1d_combined_roc.png / .pdf
"""

import os, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve
import xgboost as xgb

warnings.filterwarnings('ignore')

# ─── Paths ─────────────────────────────────────────────────────────────────────
TOM_DATA = r'D:\Research\tomato\data'
OUT_DIR = r'D:\research\cucumber\fig1'
os.makedirs(OUT_DIR, exist_ok=True)

# ─── Colors ────────────────────────────────────────────────────────────────────
COLOR_FOLDS = ['#c2a5cf', '#d4c2df', '#e0d5e0', '#c2e4c2', '#a6dba0']
COLOR_MEAN  = '#2E5AAC'  # blue
COLOR_EXT   = '#B83A3A'  # red
FILL_CV = '#c2e4c2'

# ─── Publication style ──────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 8.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'axes.linewidth': 0.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# ════════════════════════════════════════════════════════════════════════════════
# 1. Load GSE243013 training data
# ════════════════════════════════════════════════════════════════════════════════
print("=" * 60)
print("Step 1: Load GSE243013 training data")
print("=" * 60)

props = pd.read_csv(os.path.join(TOM_DATA, 'cell_proportions.csv'), index_col=0)
labels = pd.read_csv(os.path.join(TOM_DATA, 'sample_labels.csv'), index_col=0)

common = props.index.intersection(labels.index)
props = props.loc[common]
labels = labels.loc[common]
y = labels['response'].values.ravel()
X = props.values
train_features = props.columns.tolist()

n_pos = (y == 1).sum()
n_neg = (y == 0).sum()
print(f"  Samples: {len(y)} (R={n_pos}, NR={n_neg})")
print(f"  Features: {X.shape[1]} cell subtypes")

# ════════════════════════════════════════════════════════════════════════════════
# 2. 5-fold CV — capture per-fold and mean ROC
# ════════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("Step 2: 5-fold Cross-Validation")
print("=" * 60)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
mean_fpr = np.linspace(0, 1, 200)
tprs_list = []
fold_aucs = []
fold_roc_data = []

for i, (train_idx, test_idx) in enumerate(skf.split(X, y)):
    X_tr, X_te = X[train_idx], X[test_idx]
    y_tr, y_te = y[train_idx], y[test_idx]

    scaler = StandardScaler()
    X_tr_s = scaler.fit_transform(X_tr)
    X_te_s = scaler.transform(X_te)

    model = xgb.XGBClassifier(
        n_estimators=200, max_depth=4, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        reg_alpha=0.1, reg_lambda=1.0,
        random_state=42, use_label_encoder=False,
        eval_metric='logloss'
    )
    model.fit(X_tr_s, y_tr)
    y_prob = model.predict_proba(X_te_s)[:, 1]

    fpr, tpr, _ = roc_curve(y_te, y_prob)
    auc_val = roc_auc_score(y_te, y_prob)
    fold_aucs.append(auc_val)

    tpr_interp = np.interp(mean_fpr, fpr, tpr)
    tpr_interp[0] = 0.0
    tprs_list.append(tpr_interp)

    fold_roc_data.append({
        'fold': i + 1, 'fpr': fpr.tolist(), 'tpr': tpr.tolist(),
        'auc': round(auc_val, 3)
    })
    print(f"  Fold {i+1}: AUC = {auc_val:.4f}")

mean_tpr = np.mean(tprs_list, axis=0)
mean_tpr[-1] = 1.0
std_tpr = np.std(tprs_list, axis=0)
mean_auc = np.mean(fold_aucs)
std_auc = np.std(fold_aucs)
print(f"  Mean AUC = {mean_auc:.3f} ± {std_auc:.3f}")

# ════════════════════════════════════════════════════════════════════════════════
# 3. Train on full GSE243013 + predict GSE207422 (seed=1)
# ════════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("Step 3: Full model training + GSE207422 validation")
print("=" * 60)

# Load validation data
val_props = pd.read_csv(os.path.join(TOM_DATA, 'GSE207422_cell_proportions_post.csv'), index_col=0)
for f in train_features:
    if f not in val_props.columns:
        val_props[f] = 0.0
val_props = val_props[train_features]

val_labels = {
    'BD_immune03': 1, 'BD_immune06': 1, 'BD_immune11': 1, 'BD_immune14': 1,
    'BD_immune02': 0, 'BD_immune04': 0, 'BD_immune07': 0, 'BD_immune09': 0,
    'BD_immune10': 0, 'BD_immune12': 0, 'BD_immune13': 0, 'BD_immune15': 0,
}
samples = [s for s in val_props.index if s in val_labels]
X_val = val_props.loc[samples].values
y_val = np.array([val_labels[s] for s in samples])

# Train with fixed seed=1 (user's choice — avoids data peeking on n=12)
SEED_EXT = 1
scaler_val = StandardScaler()
X_train_s = scaler_val.fit_transform(X)
X_val_s = scaler_val.transform(X_val)

model_val = xgb.XGBClassifier(
    n_estimators=200, max_depth=4, learning_rate=0.05,
    subsample=0.8, colsample_bytree=0.8,
    reg_alpha=0.1, reg_lambda=1.0,
    random_state=SEED_EXT, use_label_encoder=False,
    eval_metric='logloss'
)
model_val.fit(X_train_s, y)

y_prob_val = model_val.predict_proba(X_val_s)[:, 1]
y_pred_val = model_val.predict(X_val_s)

auc_val = roc_auc_score(y_val, y_prob_val)
acc_val = (y_pred_val == y_val).mean()
fpr_val, tpr_val, _ = roc_curve(y_val, y_prob_val)

n_val_pos = y_val.sum()
n_val_neg = len(y_val) - n_val_pos

print(f"  Validation: {len(samples)} samples (R={n_val_pos}, NR={n_val_neg})")
print(f"  Seed = {SEED_EXT}")
print(f"  AUC = {auc_val:.4f}")
print(f"  ACC = {acc_val:.3f} ({int(acc_val*len(samples))}/{len(samples)} correct)")
for s, prob, pred, true in sorted(zip(samples, y_prob_val, y_pred_val, y_val),
                                   key=lambda x: -x[1]):
    mark = '✓' if pred == true else '✗'
    print(f"    {mark} {s}: prob={prob:.4f}, pred={'R' if pred else 'NR'}, "
          f"true={'R' if true else 'NR'}")

# ════════════════════════════════════════════════════════════════════════════════
# 5. Combined ROC Plot
# ════════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("Step 5: Plotting combined ROC")
print("=" * 60)

fig, ax = plt.subplots(figsize=(7.5, 7.5))

# (a) 5 individual folds
for i in range(5):
    fd = fold_roc_data[i]
    ax.plot(fd['fpr'], fd['tpr'], color=COLOR_FOLDS[i], lw=1.2, alpha=0.6,
            label=f"Fold {fd['fold']} (AUC = {fd['auc']:.3f})")

# (b) Mean ROC ± 1 SD
ax.plot(mean_fpr, mean_tpr, color=COLOR_MEAN, lw=2.5,
        label=f'Mean CV (AUC = {mean_auc:.3f} ± {std_auc:.3f})')
ax.fill_between(mean_fpr,
                np.maximum(mean_tpr - std_tpr, 0),
                np.minimum(mean_tpr + std_tpr, 1),
                color=FILL_CV, alpha=0.3, label='±1 SD')

# (c) GSE207422 external validation
ax.plot(fpr_val, tpr_val, color=COLOR_EXT, lw=2.5, linestyle='--',
        label=f'GSE207422 NSCLC (AUC = {auc_val:.3f})')
ax.fill_between(fpr_val, tpr_val, alpha=0.12, color=COLOR_EXT)

# (d) Chance diagonal
ax.plot([0, 1], [0, 1], 'k--', lw=1.2, alpha=0.4, label='Chance')

# Labels & styling
ax.set_xlabel('False Positive Rate (1 − Specificity)', fontweight='bold')
ax.set_ylabel('True Positive Rate (Sensitivity)', fontweight='bold')
ax.set_title('ROC Curves: Internal Cross-Validation + External Validation',
             fontweight='bold', pad=12)
ax.set_xlim([-0.02, 1.02])
ax.set_ylim([-0.02, 1.02])
ax.set_aspect('equal')
ax.legend(loc='lower right', fontsize=8.5, framealpha=0.9, edgecolor='#cccccc')

# Info box
info_text = (
    f'Training: GSE243013 NSCLC (n={len(y)}, R={n_pos}, NR={n_neg})\n'
    f'Validation: GSE207422 NSCLC (n={len(samples)}, R={n_val_pos}, NR={n_val_neg})\n'
    f'Model: XGBoost | 51 cell subtypes | 5-fold stratified CV'
)
ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
        ha='left', va='top', fontsize=8, color='#555555',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='#dddddd', alpha=0.9))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(labelsize=10)
ax.grid(True, alpha=0.15)

plt.tight_layout()

# Save
for ext in ['png', 'pdf']:
    path = os.path.join(OUT_DIR, f'fig1d_combined_roc.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {path}")

plt.close()
print("\nDone!")
