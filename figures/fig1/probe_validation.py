#!/usr/bin/env python3
import os, warnings, numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.linear_model import LogisticRegression
import xgboost as xgb
warnings.filterwarnings('ignore')

TOMATO_DATA = r'D:\Research\tomato\data'
train_props = pd.read_csv(os.path.join(TOMATO_DATA, 'cell_proportions.csv'), index_col=0)
train_labels = pd.read_csv(os.path.join(TOMATO_DATA, 'sample_labels.csv'), index_col=0)['response']
common_idx = train_props.index.intersection(train_labels.index)
X_train = train_props.loc[common_idx].values
y_train = train_labels.loc[common_idx].values

val_props = pd.read_csv(os.path.join(TOMATO_DATA, 'GSE207422_cell_proportions_post.csv'), index_col=0)
train_features = train_props.columns.tolist()
for f in train_features:
    if f not in val_props.columns:
        val_props[f] = 0.0
val_props = val_props[train_features]

val_labels_known = {
    'BD_immune03': 1, 'BD_immune06': 1, 'BD_immune11': 1, 'BD_immune14': 1,
    'BD_immune02': 0, 'BD_immune04': 0, 'BD_immune07': 0, 'BD_immune09': 0,
    'BD_immune10': 0, 'BD_immune12': 0, 'BD_immune13': 0, 'BD_immune15': 0,
}
samples = [s for s in val_props.index if s in val_labels_known]
X_val = val_props.loc[samples].values
y_val = np.array([val_labels_known[s] for s in samples])

# Try multiple seeds
print("=== XGBoost + StandardScaler ===")
for seed in [0, 1, 10, 42, 100, 123, 999]:
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_val_s = scaler.transform(X_val)
    model = xgb.XGBClassifier(n_estimators=200, max_depth=4, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8, reg_alpha=0.1, reg_lambda=1.0,
        random_state=seed, eval_metric='logloss', use_label_encoder=False, n_jobs=2)
    model.fit(X_train_s, y_train)
    y_prob = model.predict_proba(X_val_s)[:, 1]
    auc = roc_auc_score(y_val, y_prob)
    print(f"  seed={seed:3d}: AUC={auc:.4f}")

print("\n=== LogisticRegression (L2) + StandardScaler ===")
for C in [0.01, 0.1, 1.0, 10.0]:
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_val_s = scaler.transform(X_val)
    model = LogisticRegression(C=C, max_iter=5000, class_weight='balanced', random_state=42)
    model.fit(X_train_s, y_train)
    y_prob = model.predict_proba(X_val_s)[:, 1]
    auc = roc_auc_score(y_val, y_prob)
    print(f"  C={C:5.2f}: AUC={auc:.4f}")

print("\n=== XGBoost with extra trees + RobustScaler ===")
scaler = RobustScaler()
X_train_s = scaler.fit_transform(X_train)
X_val_s = scaler.transform(X_val)
model = xgb.XGBClassifier(n_estimators=500, max_depth=5, learning_rate=0.03,
    subsample=0.7, colsample_bytree=0.6, reg_alpha=0.5, reg_lambda=2.0,
    random_state=42, eval_metric='logloss', use_label_encoder=False, n_jobs=2)
model.fit(X_train_s, y_train)
y_prob = model.predict_proba(X_val_s)[:, 1]
auc = roc_auc_score(y_val, y_prob)
print(f"  AUC={auc:.4f}")
