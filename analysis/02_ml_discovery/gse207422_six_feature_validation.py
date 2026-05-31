"""
Validate using the 6 features from the decision tree,
normalized by immune cell count (consistent with GSE243013).
"""
import os, pandas as pd, numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve
import matplotlib.pyplot as plt
import scipy.io

BASE = r'D:\Research\tomato'
INPUT_DIR = os.path.join(BASE, 'data', 'gse207422_cellchat_input_corrected')

# GSE243013: 6 features (proportions among all cells = immune cells)
props = pd.read_csv(os.path.join(BASE, 'data', 'cell_proportions.csv'), index_col=0)
labels_243013 = pd.read_csv(os.path.join(BASE, 'data', 'sample_labels.csv'), index_col=0)['response']

FEATURES_6 = ['CD4T_Treg_CCR8', 'CD4T_Treg_MKI67', 'CD4T_Treg_FOXP3',
              'Bm_PDE4D', 'M\u03C6_CXCL10', 'M\u03C6_CXCL2']

features_243013 = props[FEATURES_6].copy()
common = features_243013.index.intersection(labels_243013.index)
X_train = features_243013.loc[common].values
y_train = labels_243013.loc[common].values

print(f"GSE243013: {len(common)} samples, 6 features")
for i, f in enumerate(FEATURES_6):
    r_mean = features_243013.loc[common][y_train==1].iloc[:, i].mean()
    nr_mean = features_243013.loc[common][y_train==0].iloc[:, i].mean()
    print(f"  {f}: R={r_mean:.4f}, NR={nr_mean:.4f}")

# GSE207422: compute 6 features from MTX + metadata
IMMUNE_TYPES = ['B_cells', 'CD4_T', 'CD8_T', 'DC', 'Mast', 'Myeloid', 'NK', 'Plasma', 'T_cells', 'Treg']

def compute_six_features(group):
    mtx = scipy.io.mmread(os.path.join(INPUT_DIR, f'{group}_matrix.mtx'))
    meta = pd.read_csv(os.path.join(INPUT_DIR, f'{group}_metadata.csv'))
    genes = pd.read_csv(os.path.join(INPUT_DIR, f'{group}_genes.tsv'), header=None, sep='\t')[0].values
    barcodes = pd.read_csv(os.path.join(INPUT_DIR, f'{group}_barcodes.tsv'), header=None, sep='\t')[0].values
    
    rows = mtx.row
    cols = mtx.col
    data = mtx.data
    
    meta['barcode'] = barcodes
    meta['sample'] = meta['barcode'].str.rsplit('_', n=1).str[0]
    
    # Gene indices
    gene_idx = {
        'MKI67': np.where(genes == 'MKI67')[0][0],
        'CCR8': np.where(genes == 'CCR8')[0][0],
        'PDE4D': np.where(genes == 'PDE4D')[0][0],
        'CXCL10': np.where(genes == 'CXCL10')[0][0],
        'CXCL2': np.where(genes == 'CXCL2')[0][0],
    }
    
    results = []
    for sample in meta['sample'].unique():
        sample_mask = meta['sample'] == sample
        sample_idx = meta.index[sample_mask].values
        
        # Immune cells
        immune_mask = meta.loc[sample_mask, 'cell_type'].isin(IMMUNE_TYPES).values
        immune_idx = sample_idx[immune_mask]
        n_immune = len(immune_idx)
        if n_immune == 0:
            continue
        
        # Treg
        treg_mask = meta.loc[sample_mask, 'cell_type'] == 'Treg'
        treg_idx = sample_idx[treg_mask.values]
        n_treg = len(treg_idx)
        
        # B_cells
        b_mask = meta.loc[sample_mask, 'cell_type'] == 'B_cells'
        b_idx = sample_idx[b_mask.values]
        n_b = len(b_idx)
        
        # Myeloid
        my_mask = meta.loc[sample_mask, 'cell_type'] == 'Myeloid'
        my_idx = sample_idx[my_mask.values]
        n_my = len(my_idx)
        
        # Helper: count positive cells for a gene in a cell set
        def count_positive(gene, cell_idx):
            gidx = gene_idx[gene]
            mask = (rows == gidx) & np.isin(cols, cell_idx)
            vals = data[mask]
            ccols = cols[mask]
            expr = dict(zip(ccols, vals))
            return sum(1 for c in cell_idx if expr.get(c, 0) > 0)
        
        n_ccr8_treg = count_positive('CCR8', treg_idx)
        n_mki67_treg = count_positive('MKI67', treg_idx)
        n_pde4d_b = count_positive('PDE4D', b_idx)
        n_cxcl10_my = count_positive('CXCL10', my_idx)
        n_cxcl2_my = count_positive('CXCL2', my_idx)
        
        results.append({
            'sample': sample,
            'response': 1 if group == 'weak' else 0,
            'CD4T_Treg_CCR8': n_ccr8_treg / n_immune,
            'CD4T_Treg_MKI67': n_mki67_treg / n_immune,
            'CD4T_Treg_FOXP3': n_treg / n_immune,
            'Bm_PDE4D': n_pde4d_b / n_immune,
            'Mφ_CXCL10': n_cxcl10_my / n_immune,
            'Mφ_CXCL2': n_cxcl2_my / n_immune,
        })
    
    return pd.DataFrame(results)

print("\n--- Computing Strong features ---")
df_strong = compute_six_features('strong')
print(df_strong[['sample'] + FEATURES_6])

print("\n--- Computing Weak features ---")
df_weak = compute_six_features('weak')
print(df_weak[['sample'] + FEATURES_6])

df_207422 = pd.concat([df_strong, df_weak], ignore_index=True)
X_test = df_207422[FEATURES_6].values
y_test = df_207422['response'].values

print(f"\nGSE207422: {len(df_207422)} samples (Weak={sum(y_test)}, Strong={len(y_test)-sum(y_test)})")
for i, f in enumerate(FEATURES_6):
    s_mean = df_207422[y_test==0].iloc[:, i+2].mean()  # +2 because first 2 cols are sample, response
    w_mean = df_207422[y_test==1].iloc[:, i+2].mean()
    print(f"  {f}: S={s_mean:.4f}, W={w_mean:.4f}")

# Model 1: Train on GSE243013, predict on GSE207422
model = LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42)
model.fit(X_train, y_train)
y_proba_direct = model.predict_proba(X_test)[:, 1]
y_pred_direct = (y_proba_direct > 0.5).astype(int)
auc_direct = roc_auc_score(y_test, y_proba_direct)
acc_direct = accuracy_score(y_test, y_pred_direct)

print(f"\n=== Direct: Train GSE243013 -> Test GSE207422 ===")
print(f"  AUC = {auc_direct:.3f}")
print(f"  Accuracy = {acc_direct:.3f} ({sum(y_pred_direct==y_test)}/{len(y_test)})")

# Model 2: LOOCV on GSE207422
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

# Plot
fig, ax = plt.subplots(figsize=(6, 6))
fpr_d, tpr_d, _ = roc_curve(y_test, y_proba_direct)
ax.plot(fpr_d, tpr_d, 'b-', lw=2.5, label=f'Train GSE243013 (AUC = {auc_direct:.3f})')
fpr_l, tpr_l, _ = roc_curve(y_test, y_proba_loo)
ax.plot(fpr_l, tpr_l, 'r--', lw=2, label=f'LOOCV GSE207422 (AUC = {auc_loo:.3f})')
ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5)
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC: 6-Feature Decision Tree Validation (External Validation Only)\nGSE207422 is used solely for independent validation; all mechanistic analyses use GSE243013', fontsize=10, fontweight='bold')
ax.legend(loc='lower right', fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()

out_dir = os.path.join(BASE, 'figures', 'demo', '1')
os.makedirs(out_dir, exist_ok=True)
fig.savefig(os.path.join(out_dir, 'fig1h_six_feature_roc.png'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(out_dir, 'fig1h_six_feature_roc.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"\nSaved: fig1h_six_feature_roc.png/pdf")
