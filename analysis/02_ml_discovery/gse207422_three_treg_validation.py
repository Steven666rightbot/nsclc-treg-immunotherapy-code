"""
Use the 3 most important Treg features from the decision tree:
  1. CD4T_Treg_CCR8  -> CCR8+ Treg / immune cells
  2. CD4T_Treg_MKI67 -> MKI67+ Treg / immune cells  
  3. CD4T_Treg_FOXP3 -> Total Treg / immune cells (all Treg are FOXP3+)

Train on GSE243013, validate on GSE207422 with corrected annotation.
"""
import os, pandas as pd, numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve
import matplotlib.pyplot as plt

BASE = r'D:\Research\tomato'

# ========== GSE243013: compute 3-feature proportions among immune cells ==========
props = pd.read_csv(os.path.join(BASE, 'data', 'cell_proportions.csv'), index_col=0)
labels_243013 = pd.read_csv(os.path.join(BASE, 'data', 'sample_labels.csv'), index_col=0)['response']

# GSE243013: all 51 columns are immune cells, so "immune proportion" = raw proportion
# Extract the 3 Treg features
features_243013 = pd.DataFrame(index=props.index)
features_243013['Treg_CCR8']  = props['CD4T_Treg_CCR8']
features_243013['Treg_MKI67'] = props['CD4T_Treg_MKI67']
features_243013['Treg_FOXP3'] = props['CD4T_Treg_FOXP3']

# Align with labels
common = features_243013.index.intersection(labels_243013.index)
X_train = features_243013.loc[common].values
y_train = labels_243013.loc[common].values

print(f"GSE243013: {len(common)} samples, 3 features")
print(f"  Treg_CCR8  mean(R)={features_243013.loc[common][y_train==1]['Treg_CCR8'].mean():.4f}, mean(NR)={features_243013.loc[common][y_train==0]['Treg_CCR8'].mean():.4f}")
print(f"  Treg_MKI67 mean(R)={features_243013.loc[common][y_train==1]['Treg_MKI67'].mean():.4f}, mean(NR)={features_243013.loc[common][y_train==0]['Treg_MKI67'].mean():.4f}")
print(f"  Treg_FOXP3 mean(R)={features_243013.loc[common][y_train==1]['Treg_FOXP3'].mean():.4f}, mean(NR)={features_243013.loc[common][y_train==0]['Treg_FOXP3'].mean():.4f}")

# ========== GSE207422: compute same 3 features from corrected annotation ==========
# Read corrected metadata (split by strong/weak)
strong_meta = pd.read_csv(os.path.join(BASE, 'data', 'gse207422_cellchat_input_corrected', 'strong_metadata.csv'))
weak_meta = pd.read_csv(os.path.join(BASE, 'data', 'gse207422_cellchat_input_corrected', 'weak_metadata.csv'))
strong_meta['response'] = 'strong'
weak_meta['response'] = 'weak'
meta = pd.concat([strong_meta, weak_meta], ignore_index=True)

# Immune cell types in GSE207422
IMMUNE_TYPES = ['B_cells', 'CD4_T', 'CD8_T', 'DC', 'Mast', 'Myeloid', 'NK', 'Plasma', 'T_cells', 'Treg']

# Extract sample ID from cell barcode (e.g., "BD_immune03_100105" -> "BD_immune03")
meta['sample'] = meta['cell'].str.rsplit('_', n=1).str[0]
samples = meta['sample'].unique()
features_207422 = pd.DataFrame(index=samples, columns=['Treg_CCR8', 'Treg_MKI67', 'Treg_FOXP3', 'response'])

for sample in samples:
    df = meta[meta['sample'] == sample]
    response = df['response_group'].iloc[0]
    
    # Total immune cells
    immune_cells = df[df['cell_type'].isin(IMMUNE_TYPES)]
    n_immune = len(immune_cells)
    
    if n_immune == 0:
        continue
    
    # Total Treg (FOXP3+ equivalent)
    treg_cells = df[df['cell_type'] == 'Treg']
    n_treg = len(treg_cells)
    
    # MKI67+ Treg
    if 'MKI67' in df.columns:
        mki67_treg = treg_cells[treg_cells['MKI67'].astype(float) > 0]
    else:
        mki67_treg = pd.DataFrame()
    n_mki67 = len(mki67_treg)
    
    # CCR8+ Treg
    if 'CCR8' in df.columns:
        ccr8_treg = treg_cells[treg_cells['CCR8'].astype(float) > 0]
    else:
        ccr8_treg = pd.DataFrame()
    n_ccr8 = len(ccr8_treg)
    
    features_207422.loc[sample, 'Treg_FOXP3'] = n_treg / n_immune
    features_207422.loc[sample, 'Treg_MKI67'] = n_mki67 / n_immune
    features_207422.loc[sample, 'Treg_CCR8']  = n_ccr8 / n_immune
    features_207422.loc[sample, 'response']   = 1 if response.lower() == 'weak' else 0

features_207422 = features_207422.dropna()
X_test = features_207422[['Treg_CCR8', 'Treg_MKI67', 'Treg_FOXP3']].values.astype(float)
y_test = features_207422['response'].values.astype(int)

print(f"\nGSE207422: {len(features_207422)} samples")
print(f"  Weak: {sum(y_test)}, Strong: {len(y_test)-sum(y_test)}")
print(f"  Treg_CCR8  mean(S)={features_207422[y_test==0]['Treg_CCR8'].mean():.4f}, mean(W)={features_207422[y_test==1]['Treg_CCR8'].mean():.4f}")
print(f"  Treg_MKI67 mean(S)={features_207422[y_test==0]['Treg_MKI67'].mean():.4f}, mean(W)={features_207422[y_test==1]['Treg_MKI67'].mean():.4f}")
print(f"  Treg_FOXP3 mean(S)={features_207422[y_test==0]['Treg_FOXP3'].mean():.4f}, mean(W)={features_207422[y_test==1]['Treg_FOXP3'].mean():.4f}")

# ========== Model: Train on GSE243013, predict on GSE207422 ==========
model = LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42)
model.fit(X_train, y_train)

# Direct prediction on GSE207422
y_pred_proba = model.predict_proba(X_test)[:, 1]
y_pred = (y_pred_proba > 0.5).astype(int)

auc_direct = roc_auc_score(y_test, y_pred_proba)
acc_direct = accuracy_score(y_test, y_pred)

print(f"\n=== Direct: Train on GSE243013, Test on GSE207422 ===")
print(f"  AUC = {auc_direct:.3f}")
print(f"  Accuracy = {acc_direct:.3f} ({sum(y_pred==y_test)}/{len(y_test)})")
print(f"  Coefficients: CCR8={model.coef_[0][0]:.3f}, MKI67={model.coef_[0][1]:.3f}, FOXP3={model.coef_[0][2]:.3f}")

# ========== LOOCV on GSE207422 (train on 13, test on 1) ==========
loo = LeaveOneOut()
y_proba_loo = np.zeros(len(y_test))

for train_idx, test_idx in loo.split(X_test):
    if len(train_idx) < 2:
        continue
    X_tr, X_te = X_test[train_idx], X_test[test_idx]
    y_tr = y_test[train_idx]
    
    m = LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42)
    m.fit(X_tr, y_tr)
    y_proba_loo[test_idx[0]] = m.predict_proba(X_te)[:, 1]

auc_loo = roc_auc_score(y_test, y_proba_loo)
acc_loo = accuracy_score(y_test, (y_proba_loo > 0.5).astype(int))

print(f"\n=== LOOCV on GSE207422 ===")
print(f"  AUC = {auc_loo:.3f}")
print(f"  Accuracy = {acc_loo:.3f}")

# ========== Plot ROC ==========
fig, ax = plt.subplots(figsize=(6, 6))

# Direct ROC
fpr_d, tpr_d, _ = roc_curve(y_test, y_pred_proba)
ax.plot(fpr_d, tpr_d, 'b-', lw=2.5, label=f'GSE243013→GSE207422 (AUC = {auc_direct:.3f})')

# LOOCV ROC
fpr_l, tpr_l, _ = roc_curve(y_test, y_proba_loo)
ax.plot(fpr_l, tpr_l, 'r--', lw=2, label=f'GSE207422 LOOCV (AUC = {auc_loo:.3f})')

ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5)
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC: 3-Treg Features Predict Response\n(CCR8+ / MKI67+ / FOXP3+ Treg as % of immune cells)', fontsize=11, fontweight='bold')
ax.legend(loc='lower right', fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim([-0.02, 1.02])
ax.set_ylim([-0.02, 1.02])
plt.tight_layout()

out_dir = os.path.join(BASE, 'figures', 'demo', '1')
os.makedirs(out_dir, exist_ok=True)
fig.savefig(os.path.join(out_dir, 'fig1h_three_treg_immune_proportion_roc.png'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(out_dir, 'fig1h_three_treg_immune_proportion_roc.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"\nSaved: fig1h_three_treg_immune_proportion_roc.png/pdf")
