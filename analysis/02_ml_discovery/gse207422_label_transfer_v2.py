"""
Label transfer GSE207422 using GSE243013 pseudobulk reference.
Vectorized Spearman correlation for speed.
"""
import os, sys, pandas as pd, numpy as np
import scipy.io
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve
import matplotlib.pyplot as plt

BASE = r'D:\Research\tomato'
INPUT_DIR = os.path.join(BASE, 'data', 'gse207422_cellchat_input_corrected')

def flush_print(msg):
    print(msg)
    sys.stdout.flush()

# ========== 1. Load reference ==========
flush_print("[1/7] Loading GSE243013 reference pseudobulk...")
ref = pd.read_csv(os.path.join(BASE, 'data', 'GSE243013_reference_pseudobulk.csv'), index_col=0)
flush_print(f"  Reference: {ref.shape[0]} subtypes x {ref.shape[1]} genes")

# ========== 2. Find marker genes ==========
flush_print("[2/7] Finding marker genes...")
all_markers = set()
for subtype in ref.index:
    other_mean = ref.drop(subtype).mean()
    fc = ref.loc[subtype] - other_mean
    top = fc.nlargest(200)
    all_markers.update(top.index)
all_markers = sorted(all_markers)
flush_print(f"  Total unique marker genes: {len(all_markers)}")

# ========== 3. Read GSE207422 genes ==========
flush_print("[3/7] Reading GSE207422 gene lists...")
genes = pd.read_csv(os.path.join(INPUT_DIR, 'strong_genes.tsv'), header=None, sep='\t')[0].values

common_markers = [g for g in all_markers if g in genes]
flush_print(f"  Marker genes in GSE207422: {len(common_markers)}")

# Build reference expression: genes x subtypes
ref_expr = ref[common_markers].T.values  # n_genes x 51

# ========== 4. Helper: process one group ==========
IMMUNE_TYPES = ['B_cells', 'CD4_T', 'CD8_T', 'DC', 'Mast', 'Myeloid', 'NK', 'Plasma', 'T_cells', 'Treg']

def process_group(prefix, group_name):
    flush_print(f"\n--- Processing {group_name} ---")
    
    # Read metadata
    meta = pd.read_csv(os.path.join(INPUT_DIR, f'{prefix}_metadata.csv'))
    meta_filtered = meta[meta['cell_type'].isin(IMMUNE_TYPES)].copy()
    n_immune = len(meta_filtered)
    n_total = len(meta)
    flush_print(f"  Immune cells: {n_immune:,} / {n_total:,}")
    
    if n_immune == 0:
        return None
    
    # Read MTX
    flush_print(f"  Reading MTX...")
    mtx = scipy.io.mmread(os.path.join(INPUT_DIR, f'{prefix}_matrix.mtx'))
    rows, cols, data = mtx.row, mtx.col, mtx.data
    
    # Gene indices
    gene_idx_map = {g: i for i, g in enumerate(genes)}
    target_gene_idx = np.array([gene_idx_map[g] for g in common_markers])
    
    # Cell indices (immune cells)
    cell_mask = meta['cell_type'].isin(IMMUNE_TYPES).values
    target_cells = np.where(cell_mask)[0]
    
    # Filter COO
    flush_print(f"  Filtering to {len(common_markers)} genes x {len(target_cells)} cells...")
    gene_filter = np.isin(rows, target_gene_idx)
    cell_filter = np.isin(cols, target_cells)
    mask = gene_filter & cell_filter
    
    row_map = {old: new for new, old in enumerate(target_gene_idx)}
    col_map = {old: new for new, old in enumerate(target_cells)}
    
    expr = np.zeros((len(common_markers), len(target_cells)))
    for r, c, v in zip(rows[mask], cols[mask], data[mask]):
        expr[row_map[r], col_map[c]] = v
    
    flush_print(f"  Expression matrix: {expr.shape}, density={np.count_nonzero(expr)/expr.size:.3f}")
    
    # Vectorized Spearman correlation
    flush_print(f"  Computing Spearman correlations (vectorized)...")
    
    # Concatenate reference + cells, compute ranks
    combined = np.concatenate([ref_expr, expr], axis=1)  # genes x (51 + n_cells)
    
    # Rank each column (axis=0)
    flush_print(f"  Ranking...")
    ranks = rankdata(combined, axis=0)  # genes x (51 + n_cells)
    
    # Split back
    ref_ranks = ranks[:, :51]  # genes x 51
    cell_ranks = ranks[:, 51:]  # genes x n_cells
    
    # Center
    ref_centered = ref_ranks - ref_ranks.mean(axis=0)
    cell_centered = cell_ranks - cell_ranks.mean(axis=0)
    
    # Pearson correlation of ranks = Spearman correlation
    flush_print(f"  Correlating...")
    numerator = cell_centered.T @ ref_centered  # n_cells x 51
    denom_cells = np.sqrt(np.sum(cell_centered**2, axis=0)).reshape(-1, 1)  # n_cells x 1
    denom_ref = np.sqrt(np.sum(ref_centered**2, axis=0)).reshape(1, -1)  # 1 x 51
    corr = numerator / (denom_cells @ denom_ref)
    
    # Assign labels
    best_idx = np.argmax(corr, axis=1)
    best_score = np.max(corr, axis=1)
    
    labels = [ref.index[i] for i in best_idx]
    flush_print(f"  Mean correlation score: {best_score.mean():.3f}")
    flush_print(f"  Score range: [{best_score.min():.3f}, {best_score.max():.3f}]")
    
    meta_filtered['label'] = labels
    meta_filtered['score'] = best_score
    
    return meta_filtered

# ========== 5. Process Strong and Weak ==========
strong_meta = process_group('strong', 'Strong')
weak_meta = process_group('weak', 'Weak')

# ========== 6. Compute proportions ==========
flush_print("\n[6/7] Computing proportions...")

def compute_props(meta_filtered, group_name):
    meta_filtered['sample'] = meta_filtered['cell'].str.rsplit('_', n=1).str[0]
    samples = meta_filtered['sample'].unique()
    rows = []
    for sample in samples:
        df = meta_filtered[meta_filtered['sample'] == sample]
        n_immune = len(df)
        if n_immune == 0:
            continue
        props = {'sample': sample, 'response': 1 if group_name == 'weak' else 0, 'n_immune': n_immune}
        for subtype in ref.index:
            props[subtype] = (df['label'] == subtype).sum() / n_immune
        rows.append(props)
    return pd.DataFrame(rows)

strong_props = compute_props(strong_meta, 'strong')
weak_props = compute_props(weak_meta, 'weak')
all_props = pd.concat([strong_props, weak_props], ignore_index=True)

flush_print(f"  Strong samples: {len(strong_props)}")
flush_print(f"  Weak samples: {len(weak_props)}")

# ========== 7. Validation ==========
flush_print("\n[7/7] Validation...")

props_243013 = pd.read_csv(os.path.join(BASE, 'data', 'cell_proportions.csv'), index_col=0)
labels_243013 = pd.read_csv(os.path.join(BASE, 'data', 'sample_labels.csv'), index_col=0)['response']

feature_cols = list(ref.index)
X_train = props_243013[feature_cols].values
y_train = labels_243013.values

X_test = all_props[feature_cols].values
y_test = all_props['response'].values

# Direct
model = LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42)
model.fit(X_train, y_train)
y_proba_direct = model.predict_proba(X_test)[:, 1]
auc_direct = roc_auc_score(y_test, y_proba_direct)
acc_direct = accuracy_score(y_test, (y_proba_direct > 0.5).astype(int))
flush_print(f"Direct: AUC={auc_direct:.3f}, Acc={acc_direct:.3f}")

# LOOCV
loo = LeaveOneOut()
y_proba_loo = np.zeros(len(y_test))
for tr_idx, te_idx in loo.split(X_test):
    if len(tr_idx) < 2:
        continue
    m = LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42)
    m.fit(X_test[tr_idx], y_test[tr_idx])
    y_proba_loo[te_idx[0]] = m.predict_proba(X_test[te_idx])[0, 1]

auc_loo = roc_auc_score(y_test, y_proba_loo)
acc_loo = accuracy_score(y_test, (y_proba_loo > 0.5).astype(int))
flush_print(f"LOOCV:  AUC={auc_loo:.3f}, Acc={acc_loo:.3f}")

# Top features
top_idx = np.argsort(np.abs(model.coef_[0]))[::-1][:10]
flush_print("\nTop 10 features by |coef|:")
for i in top_idx:
    flush_print(f"  {feature_cols[i]}: {model.coef_[0][i]:.3f}")

# Plot
fig, ax = plt.subplots(figsize=(6, 6))
fpr_d, tpr_d, _ = roc_curve(y_test, y_proba_direct)
ax.plot(fpr_d, tpr_d, 'b-', lw=2.5, label=f'Train GSE243013 (AUC = {auc_direct:.3f})')
fpr_l, tpr_l, _ = roc_curve(y_test, y_proba_loo)
ax.plot(fpr_l, tpr_l, 'r--', lw=2, label=f'LOOCV GSE207422 (AUC = {auc_loo:.3f})')
ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5)
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC: 51-Subtype Label Transfer\n(GSE243013 ref -> GSE207422)', fontsize=11, fontweight='bold')
ax.legend(loc='lower right', fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()

out_dir = os.path.join(BASE, 'figures', 'demo', '1')
os.makedirs(out_dir, exist_ok=True)
fig.savefig(os.path.join(out_dir, 'fig1h_51subtype_label_transfer_roc.png'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(out_dir, 'fig1h_51subtype_label_transfer_roc.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

# Save
all_props.to_csv(os.path.join(BASE, 'data', 'gse207422_51subtype_proportions.csv'), index=False)
flush_print(f"\nSaved: data/gse207422_51subtype_proportions.csv")
flush_print(f"Saved: figures/demo/1/fig1h_51subtype_label_transfer_roc.png")
