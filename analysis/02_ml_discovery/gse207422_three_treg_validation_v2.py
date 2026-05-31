"""
Use 3 Treg features (decision tree top features) normalized by immune cell count:
  1. CCR8+ Treg / immune cells
  2. MKI67+ Treg / immune cells
  3. Total Treg / immune cells (FOXP3+ equivalent)

Train on GSE243013, validate on GSE207422 corrected annotation.
Reads MTX + metadata directly.
"""
import os, pandas as pd, numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve
import matplotlib.pyplot as plt
import scipy.io

BASE = r'D:\Research\tomato'
INPUT_DIR = os.path.join(BASE, 'data', 'gse207422_cellchat_input_corrected')

# ========== GSE243013: 3 features (already proportions among all cells = immune cells) ==========
props = pd.read_csv(os.path.join(BASE, 'data', 'cell_proportions.csv'), index_col=0)
labels_243013 = pd.read_csv(os.path.join(BASE, 'data', 'sample_labels.csv'), index_col=0)['response']

features_243013 = pd.DataFrame(index=props.index)
features_243013['Treg_CCR8']  = props['CD4T_Treg_CCR8']
features_243013['Treg_MKI67'] = props['CD4T_Treg_MKI67']
features_243013['Treg_FOXP3'] = props['CD4T_Treg_FOXP3']

common = features_243013.index.intersection(labels_243013.index)
X_train = features_243013.loc[common].values
y_train = labels_243013.loc[common].values

print(f"GSE243013: {len(common)} samples")
print(f"  Treg_CCR8  R={features_243013.loc[common][y_train==1]['Treg_CCR8'].mean():.4f}, NR={features_243013.loc[common][y_train==0]['Treg_CCR8'].mean():.4f}")
print(f"  Treg_MKI67 R={features_243013.loc[common][y_train==1]['Treg_MKI67'].mean():.4f}, NR={features_243013.loc[common][y_train==0]['Treg_MKI67'].mean():.4f}")
print(f"  Treg_FOXP3 R={features_243013.loc[common][y_train==1]['Treg_FOXP3'].mean():.4f}, NR={features_243013.loc[common][y_train==0]['Treg_FOXP3'].mean():.4f}")

# ========== Helper: compute GSE207422 3 features from MTX ==========
IMMUNE_TYPES = ['B_cells', 'CD4_T', 'CD8_T', 'DC', 'Mast', 'Myeloid', 'NK', 'Plasma', 'T_cells', 'Treg']

def compute_features(group):
    """Compute 3 Treg features for strong or weak group."""
    mtx = scipy.io.mmread(os.path.join(INPUT_DIR, f'{group}_matrix.mtx'))
    meta = pd.read_csv(os.path.join(INPUT_DIR, f'{group}_metadata.csv'))
    genes = pd.read_csv(os.path.join(INPUT_DIR, f'{group}_genes.tsv'), header=None, sep='\t')[0].values
    barcodes = pd.read_csv(os.path.join(INPUT_DIR, f'{group}_barcodes.tsv'), header=None, sep='\t')[0].values
    
    # COO format: row=gene_idx, col=cell_idx
    rows = mtx.row
    cols = mtx.col
    data = mtx.data
    
    # Sample IDs from barcodes
    meta['barcode'] = barcodes
    meta['sample'] = meta['barcode'].str.rsplit('_', n=1).str[0]
    
    # Gene indices
    mki67_idx = np.where(genes == 'MKI67')[0][0]
    ccr8_idx = np.where(genes == 'CCR8')[0][0]
    
    results = []
    for sample in meta['sample'].unique():
        sample_mask = meta['sample'] == sample
        sample_idx = meta.index[sample_mask].values  # positions in meta/barcodes
        
        # Immune cells in this sample
        immune_mask = meta.loc[sample_mask, 'cell_type'].isin(IMMUNE_TYPES).values
        immune_idx = sample_idx[immune_mask]
        n_immune = len(immune_idx)
        
        if n_immune == 0:
            continue
        
        # Treg cells in this sample
        treg_mask = meta.loc[sample_mask, 'cell_type'] == 'Treg'
        treg_idx = sample_idx[treg_mask.values]
        n_treg = len(treg_idx)
        
        # MKI67+ Treg: MKI67 expression > 0 in Treg cells
        # Find COO entries where row=mki67_idx and col in treg_idx
        mki67_mask = (rows == mki67_idx) & np.isin(cols, treg_idx)
        mki67_vals = data[mki67_mask]
        mki67_cols = cols[mki67_mask]
        # Build dict: col -> expr (0 if not in matrix)
        mki67_expr = dict(zip(mki67_cols, mki67_vals))
        n_mki67 = sum(1 for c in treg_idx if mki67_expr.get(c, 0) > 0)
        
        # CCR8+ Treg
        ccr8_mask = (rows == ccr8_idx) & np.isin(cols, treg_idx)
        ccr8_vals = data[ccr8_mask]
        ccr8_cols = cols[ccr8_mask]
        ccr8_expr = dict(zip(ccr8_cols, ccr8_vals))
        n_ccr8 = sum(1 for c in treg_idx if ccr8_expr.get(c, 0) > 0)
        
        results.append({
            'sample': sample,
            'response': 1 if group == 'weak' else 0,
            'Treg_CCR8': n_ccr8 / n_immune,
            'Treg_MKI67': n_mki67 / n_immune,
            'Treg_FOXP3': n_treg / n_immune,
            'n_immune': n_immune,
            'n_treg': n_treg,
            'n_ccr8': n_ccr8,
            'n_mki67': n_mki67
        })
    
    return pd.DataFrame(results)

# ========== Compute GSE207422 features ==========
print("\n--- Computing Strong features ---")
df_strong = compute_features('strong')
print(f"Strong: {len(df_strong)} samples")
print(df_strong[['sample', 'Treg_CCR8', 'Treg_MKI67', 'Treg_FOXP3', 'n_treg', 'n_ccr8', 'n_mki67']])

print("\n--- Computing Weak features ---")
df_weak = compute_features('weak')
print(f"Weak: {len(df_weak)} samples")
print(df_weak[['sample', 'Treg_CCR8', 'Treg_MKI67', 'Treg_FOXP3', 'n_treg', 'n_ccr8', 'n_mki67']])

df_207422 = pd.concat([df_strong, df_weak], ignore_index=True)
X_test = df_207422[['Treg_CCR8', 'Treg_MKI67', 'Treg_FOXP3']].values
y_test = df_207422['response'].values

print(f"\nGSE207422: {len(df_207422)} samples (Weak={sum(y_test)}, Strong={len(y_test)-sum(y_test)})")
print(f"  Treg_CCR8  S={df_207422[y_test==0]['Treg_CCR8'].mean():.4f}, W={df_207422[y_test==1]['Treg_CCR8'].mean():.4f}")
print(f"  Treg_MKI67 S={df_207422[y_test==0]['Treg_MKI67'].mean():.4f}, W={df_207422[y_test==1]['Treg_MKI67'].mean():.4f}")
print(f"  Treg_FOXP3 S={df_207422[y_test==0]['Treg_FOXP3'].mean():.4f}, W={df_207422[y_test==1]['Treg_FOXP3'].mean():.4f}")

# ========== Model 1: Train on GSE243013, predict on GSE207422 ==========
model = LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42)
model.fit(X_train, y_train)

y_proba_direct = model.predict_proba(X_test)[:, 1]
y_pred_direct = (y_proba_direct > 0.5).astype(int)
auc_direct = roc_auc_score(y_test, y_proba_direct)
acc_direct = accuracy_score(y_test, y_pred_direct)

print(f"\n=== Direct: Train GSE243013 -> Test GSE207422 ===")
print(f"  AUC = {auc_direct:.3f}")
print(f"  Accuracy = {acc_direct:.3f} ({sum(y_pred_direct==y_test)}/{len(y_test)})")
print(f"  Coef: CCR8={model.coef_[0][0]:.2f}, MKI67={model.coef_[0][1]:.2f}, FOXP3={model.coef_[0][2]:.2f}")

# ========== Model 2: LOOCV on GSE207422 ==========
loo = LeaveOneOut()
y_proba_loo = np.zeros(len(y_test))

for train_idx, test_idx in loo.split(X_test):
    if len(train_idx) < 2:
        continue
    X_tr, X_te = X_test[train_idx], X_test[test_idx]
    y_tr = y_test[train_idx]
    m = LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42)
    m.fit(X_tr, y_tr)
    y_proba_loo[test_idx[0]] = m.predict_proba(X_te)[0, 1]

auc_loo = roc_auc_score(y_test, y_proba_loo)
acc_loo = accuracy_score(y_test, (y_proba_loo > 0.5).astype(int))

print(f"\n=== LOOCV on GSE207422 ===")
print(f"  AUC = {auc_loo:.3f}")
print(f"  Accuracy = {acc_loo:.3f}")

# ========== Plot ==========
fig, ax = plt.subplots(figsize=(6, 6))

fpr_d, tpr_d, _ = roc_curve(y_test, y_proba_direct)
ax.plot(fpr_d, tpr_d, 'b-', lw=2.5, label=f'Train GSE243013 (AUC = {auc_direct:.3f})')

fpr_l, tpr_l, _ = roc_curve(y_test, y_proba_loo)
ax.plot(fpr_l, tpr_l, 'r--', lw=2, label=f'LOOCV GSE207422 (AUC = {auc_loo:.3f})')

ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5)
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC: 3 Treg Features Predict Response\n(CCR8+ / MKI67+ / FOXP3+ Treg as % of immune cells)', fontsize=11, fontweight='bold')
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
