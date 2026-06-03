#!/usr/bin/env python3
"""Bootstrap 95% CI for GSE207422 external validation AUC.
Handles edge cases where bootstrap samples have only one class.
Outputs: 95% CI (if computable), seed sensitivity across 20 seeds.
"""
import os, warnings, numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
import xgboost as xgb

warnings.filterwarnings('ignore')

TOMATO_DATA = r'D:\Research\tomato\data'

# --- Load data ---
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

val_labels = {
    'BD_immune03': 1, 'BD_immune06': 1, 'BD_immune11': 1, 'BD_immune14': 1,
    'BD_immune02': 0, 'BD_immune04': 0, 'BD_immune07': 0, 'BD_immune09': 0,
    'BD_immune10': 0, 'BD_immune12': 0, 'BD_immune13': 0, 'BD_immune15': 0,
}
samples = [s for s in val_props.index if s in val_labels]
X_val = val_props.loc[samples].values
y_val = np.array([val_labels[s] for s in samples])

scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)
X_val_s = scaler.transform(X_val)

# --- Seed sensitivity (20 seeds) ---
print("=== Seed Sensitivity Analysis (20 seeds) ===")
seed_aucs = []
for s in range(1, 21):
    m = xgb.XGBClassifier(
        n_estimators=200, max_depth=4, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        reg_alpha=0.1, reg_lambda=1.0,
        random_state=s, eval_metric='logloss',
        use_label_encoder=False, n_jobs=4
    )
    m.fit(X_train_s, y_train)
    yp = m.predict_proba(X_val_s)[:, 1]
    a = roc_auc_score(y_val, yp)
    seed_aucs.append(a)
    print(f"  seed={s:2d}: AUC={a:.3f}")

print(f"\n  Primary (seed=1): {seed_aucs[0]:.3f}")
print(f"  Mean ± SD: {np.mean(seed_aucs):.3f} ± {np.std(seed_aucs):.3f}")
print(f"  Range: {min(seed_aucs):.3f} – {max(seed_aucs):.3f}")

# --- Bootstrap 95% CI ---
# For n=12 (4R/8NR), some bootstrap replicates will have only one class.
# We sample with replacement within each class (stratified bootstrap).
print("\n=== Stratified Bootstrap 95% CI (10,000 iterations, seed=1) ===")
model = xgb.XGBClassifier(
    n_estimators=200, max_depth=4, learning_rate=0.05,
    subsample=0.8, colsample_bytree=0.8,
    reg_alpha=0.1, reg_lambda=1.0,
    random_state=1, eval_metric='logloss',
    use_label_encoder=False, n_jobs=4
)
model.fit(X_train_s, y_train)
yp_ref = model.predict_proba(X_val_s)[:, 1]
auc_ref = roc_auc_score(y_val, yp_ref)

rng = np.random.default_rng(42)
n_val = len(y_val)
boot_aucs = []
n_skipped = 0

for i in range(10000):
    # Stratified: sample with replacement within each class
    r_idx = rng.choice(np.where(y_val == 1)[0], size=4, replace=True)
    nr_idx = rng.choice(np.where(y_val == 0)[0], size=8, replace=True)
    idx = np.concatenate([r_idx, nr_idx])
    
    y_boot = y_val[idx]
    X_boot = X_val_s[idx]
    yp_boot = model.predict_proba(X_boot)[:, 1]
    
    try:
        a = roc_auc_score(y_boot, yp_boot)
        if not np.isnan(a):
            boot_aucs.append(a)
        else:
            n_skipped += 1
    except:
        n_skipped += 1

boot_aucs = np.array(boot_aucs)
ci_low = np.percentile(boot_aucs, 2.5)
ci_high = np.percentile(boot_aucs, 97.5)

print(f"Primary AUC (seed=1): {auc_ref:.3f}")
print(f"Bootstrap 95% CI:     [{ci_low:.3f}, {ci_high:.3f}]")
print(f"Bootstrap mean ± SD:  {boot_aucs.mean():.3f} ± {boot_aucs.std():.3f}")
print(f"Valid iterations:     {len(boot_aucs)} / 10000 (skipped {n_skipped} with single class)")

# --- Summary ---
print("\n" + "=" * 60)
print("MANUSCRIPT TEXT:")
print("=" * 60)
print(f"External validation (GSE207422, n=12): AUC = {seed_aucs[0]:.3f}")
print(f"(95% bootstrap CI [{ci_low:.3f}, {ci_high:.3f}]; ")
print(f"seed sensitivity across 20 seeds: mean {np.mean(seed_aucs):.3f} ± {np.std(seed_aucs):.3f}, ")
print(f"range {min(seed_aucs):.3f}–{max(seed_aucs):.3f})")
