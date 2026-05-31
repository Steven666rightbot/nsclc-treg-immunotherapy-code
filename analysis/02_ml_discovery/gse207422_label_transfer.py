"""
Label transfer GSE207422 using GSE243013 pseudobulk reference.
Annotate immune cells into 51 subtypes using Spearman correlation.
"""
import os, pandas as pd, numpy as np
import scipy.io
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve
import matplotlib.pyplot as plt

BASE = r'D:\Research\tomato'
INPUT_DIR = os.path.join(BASE, 'data', 'gse207422_cellchat_input_corrected')

def find_markers(ref, n_markers=200):
    """Find top marker genes for each subtype from pseudobulk."""
    markers = {}
    for subtype in ref.index:
        other_mean = ref.drop(subtype).mean()
        fc = ref.loc[subtype] - other_mean
        top = fc.nlargest(n_markers)
        markers[subtype] = list(top.index)
    # Union of all markers
    all_markers = set()
    for m in markers.values():
        all_markers.update(m)
    return markers, list(all_markers)

def read_mtx_filtered(prefix, meta, gene_names, target_genes, cell_types=None):
    """Read MTX and filter to target genes and cells."""
    mtx = scipy.io.mmread(os.path.join(INPUT_DIR, f'{prefix}_matrix.mtx'))
    
    # Find gene indices
    gene_idx_map = {g: i for i, g in enumerate(gene_names)}
    target_gene_idx = [gene_idx_map[g] for g in target_genes if g in gene_idx_map]
    
    if cell_types:
        cell_mask = meta['cell_type'].isin(cell_types).values
        target_cells = np.where(cell_mask)[0]
    else:
        target_cells = np.arange(mtx.shape[1])
    
    # COO format
    rows = mtx.row
    cols = mtx.col
    data = mtx.data
    
    # Filter
    gene_filter = np.isin(rows, target_gene_idx)
    cell_filter = np.isin(cols, target_cells)
    mask = gene_filter & cell_filter
    
    # Build dense matrix: genes x cells
    row_map = {old: new for new, old in enumerate(target_gene_idx)}
    col_map = {old: new for new, old in enumerate(target_cells)}
    
    result = np.zeros((len(target_gene_idx), len(target_cells)))
    for r, c, v in zip(rows[mask], cols[mask], data[mask]):
        result[row_map[r], col_map[c]] = v
    
    return result, target_cells

def label_transfer(expr, ref_expr):
    """
    expr: genes x cells (dense)
    ref_expr: genes x 51 subtypes (dense)
    Returns: labels for each cell
    """
    n_cells = expr.shape[1]
    labels = []
    scores = []
    
    for i in range(n_cells):
        cell_expr = expr[:, i]
        # Spearman correlation with each reference
        corrs = []
        for j in range(ref_expr.shape[1]):
            r, _ = spearmanr(cell_expr, ref_expr[:, j])
            corrs.append(r)
        corrs = np.array(corrs)
        best_idx = np.argmax(corrs)
        labels.append(best_idx)
        scores.append(corrs[best_idx])
    
    return np.array(labels), np.array(scores)

# ========== 1. Load reference ==========
print("Loading GSE243013 reference pseudobulk...")
ref = pd.read_csv(os.path.join(BASE, 'data', 'GSE243013_reference_pseudobulk.csv'), index_col=0)
print(f"  Reference: {ref.shape[0]} subtypes x {ref.shape[1]} genes")

# ========== 2. Find marker genes ==========
print("Finding marker genes...")
markers, all_markers = find_markers(ref, n_markers=200)
print(f"  Total unique marker genes: {len(all_markers)}")

# ========== 3. Read GSE207422 genes ==========
print("Reading GSE207422 gene lists...")
genes_s = pd.read_csv(os.path.join(INPUT_DIR, 'strong_genes.tsv'), header=None, sep='\t')[0].values
genes_w = pd.read_csv(os.path.join(INPUT_DIR, 'weak_genes.tsv'), header=None, sep='\t')[0].values
assert list(genes_s) == list(genes_w), "Strong and Weak gene lists differ!"
genes = genes_s

# Common markers
common_markers = [g for g in all_markers if g in genes]
print(f"  Marker genes in GSE207422: {len(common_markers)}")

# Build reference expression matrix for common markers
ref_expr = ref[common_markers].T.values  # genes x subtypes

# ========== 4. Read metadata ==========
print("Reading metadata...")
strong_meta = pd.read_csv(os.path.join(INPUT_DIR, 'strong_metadata.csv'))
weak_meta = pd.read_csv(os.path.join(INPUT_DIR, 'weak_metadata.csv'))
strong_meta['group'] = 'strong'
weak_meta['group'] = 'weak'

# Immune cell types
IMMUNE_TYPES = ['B_cells', 'CD4_T', 'CD8_T', 'DC', 'Mast', 'Myeloid', 'NK', 'Plasma', 'T_cells', 'Treg']

# ========== 5. Label transfer for Strong ==========
print("\n--- Label transfer: Strong ---")
strong_meta_filtered = strong_meta[strong_meta['cell_type'].isin(IMMUNE_TYPES)].copy()
strong_meta_filtered['orig_index'] = strong_meta_filtered.index
print(f"  Immune cells: {len(strong_meta_filtered)} / {len(strong_meta)}")

strong_expr, strong_cells = read_mtx_filtered('strong', strong_meta, genes, common_markers, IMMUNE_TYPES)
print(f"  Expression matrix: {strong_expr.shape}")

strong_labels, strong_scores = label_transfer(strong_expr, ref_expr)
strong_meta_filtered['label'] = [ref.index[i] for i in strong_labels]
strong_meta_filtered['score'] = strong_scores
print(f"  Mean correlation score: {strong_scores.mean():.3f}")

# ========== 6. Label transfer for Weak ==========
print("\n--- Label transfer: Weak ---")
weak_meta_filtered = weak_meta[weak_meta['cell_type'].isin(IMMUNE_TYPES)].copy()
weak_meta_filtered['orig_index'] = weak_meta_filtered.index
print(f"  Immune cells: {len(weak_meta_filtered)} / {len(weak_meta)}")

weak_expr, weak_cells = read_mtx_filtered('weak', weak_meta, genes, common_markers, IMMUNE_TYPES)
print(f"  Expression matrix: {weak_expr.shape}")

weak_labels, weak_scores = label_transfer(weak_expr, ref_expr)
weak_meta_filtered['label'] = [ref.index[i] for i in weak_labels]
weak_meta_filtered['score'] = weak_scores
print(f"  Mean correlation score: {weak_scores.mean():.3f}")

# ========== 7. Compute proportions per sample ==========
print("\n--- Computing proportions ---")

def compute_props(meta_filtered, group_name):
    meta_filtered['sample'] = meta_filtered['cell'].str.rsplit('_', n=1).str[0]
    samples = meta_filtered['sample'].unique()
    
    rows = []
    for sample in samples:
        df = meta_filtered[meta_filtered['sample'] == sample]
        n_immune = len(df)
        if n_immune == 0:
            continue
        
        props = {}
        for subtype in ref.index:
            props[subtype] = (df['label'] == subtype).sum() / n_immune
        
        props['sample'] = sample
        props['response'] = 1 if group_name == 'weak' else 0
        props['n_immune'] = n_immune
        rows.append(props)
    
    return pd.DataFrame(rows)

strong_props = compute_props(strong_meta_filtered, 'strong')
weak_props = compute_props(weak_meta_filtered, 'weak')
all_props = pd.concat([strong_props, weak_props], ignore_index=True)

print(f"  Strong samples: {len(strong_props)}")
print(f"  Weak samples: {len(weak_props)}")

# ========== 8. Validation ==========
print("\n--- Validation ---")

# GSE243013 training data
props_243013 = pd.read_csv(os.path.join(BASE, 'data', 'cell_proportions.csv'), index_col=0)
labels_243013 = pd.read_csv(os.path.join(BASE, 'data', 'sample_labels.csv'), index_col=0)['response']

# Use all 51 features
feature_cols = list(ref.index)
X_train = props_243013[feature_cols].values
y_train = labels_243013.values

# GSE207422 test
X_test = all_props[feature_cols].values
y_test = all_props['response'].values

print(f"GSE243013: {len(y_train)} samples, 51 features")
print(f"GSE207422: {len(y_test)} samples")

# Model 1: Train on GSE243013, test on GSE207422
model = LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42)
model.fit(X_train, y_train)
y_proba_direct = model.predict_proba(X_test)[:, 1]
auc_direct = roc_auc_score(y_test, y_proba_direct)
acc_direct = accuracy_score(y_test, (y_proba_direct > 0.5).astype(int))
print(f"\nDirect: Train GSE243013 -> Test GSE207422: AUC={auc_direct:.3f}, Acc={acc_direct:.3f}")

# Model 2: LOOCV on GSE207422
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
print(f"LOOCV on GSE207422: AUC={auc_loo:.3f}, Acc={acc_loo:.3f}")

# Non-zero coefficients
for i, feat in enumerate(feature_cols):
    if model.coef_[0][i] != 0:
        print(f"  Coef {feat}: {model.coef_[0][i]:.3f}")

# ========== 9. Plot ==========
fig, ax = plt.subplots(figsize=(6, 6))

fpr_d, tpr_d, _ = roc_curve(y_test, y_proba_direct)
ax.plot(fpr_d, tpr_d, 'b-', lw=2.5, label=f'Train GSE243013 (AUC = {auc_direct:.3f})')

fpr_l, tpr_l, _ = roc_curve(y_test, y_proba_loo)
ax.plot(fpr_l, tpr_l, 'r--', lw=2, label=f'LOOCV GSE207422 (AUC = {auc_loo:.3f})')

ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5)
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC: 51-Subtype Label Transfer Validation\n(GSE243013 ref -> GSE207422)', fontsize=11, fontweight='bold')
ax.legend(loc='lower right', fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()

out_dir = os.path.join(BASE, 'figures', 'demo', '1')
os.makedirs(out_dir, exist_ok=True)
fig.savefig(os.path.join(out_dir, 'fig1h_51subtype_label_transfer_roc.png'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(out_dir, 'fig1h_51subtype_label_transfer_roc.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

# Save proportions
all_props.to_csv(os.path.join(BASE, 'data', 'gse207422_51subtype_proportions.csv'), index=False)
print(f"\nSaved proportions to: data/gse207422_51subtype_proportions.csv")
print(f"Saved ROC plot to: figures/demo/1/fig1h_51subtype_label_transfer_roc.png")
