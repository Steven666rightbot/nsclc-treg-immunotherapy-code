#!/usr/bin/env python3
"""
Full XGBoost model validation on GSE120575 (melanoma)
Train: GSE243013 (NSCLC, 242 samples, 51 subtypes)
Test:  GSE120575 (melanoma, 48 biopsies)

Pipeline:
  1. Load GSE243013 reference pseudobulk (51 cell types × genes)
  2. Load GSE120575 expression matrix + metadata
  3. Map each GSE120575 cell to 51-subtype via cosine similarity (shared genes)
  4. Compute per-sample proportions
  5. Train XGBoost on GSE243013 → predict on GSE120575 → AUC

Output:
  - full_model_validation_results.csv
  - full_model_validation_roc.png / .pdf
"""

import os
import gc
import warnings
import gzip

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score, roc_curve, accuracy_score, confusion_matrix
)
from sklearn.metrics.pairwise import cosine_similarity
import xgboost as xgb

warnings.filterwarnings('ignore')

# ─── Paths ───────────────────────────────────────────────────────────────────
TOMATO_DATA = r'D:\Research\tomato\data'
CROSS_DATA = r'D:\research\cucumber\fig1\cross_cancer\data'
OUT_DIR = r'D:\research\cucumber\fig1\cross_cancer'
os.makedirs(OUT_DIR, exist_ok=True)

EXPR_FILE = os.path.join(CROSS_DATA, 'expression.txt.gz')
META_FILE = os.path.join(CROSS_DATA, 'metadata.txt.gz')

# Figure style (match Figure 1H)
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 12,
    'axes.facecolor': 'white',
    'figure.facecolor': 'white',
})

SEED = 42

print("=" * 60)
print("Full XGBoost validation on GSE120575 (melanoma)")
print("=" * 60)

# ═════════════════════════════════════════════════════════════════════════════
# 1. Load GSE243013 reference pseudobulk (51 cell types)
# ═════════════════════════════════════════════════════════════════════════════
print("\n[1] Loading GSE243013 reference pseudobulk...")
ref = pd.read_csv(
    os.path.join(TOMATO_DATA, 'GSE243013_reference_pseudobulk.csv'),
    index_col=0
)
cell_types = ref.index.tolist()
genes_ref = set(ref.columns.tolist())
print(f"  Reference: {len(cell_types)} cell types, {len(genes_ref)} genes")

# ═════════════════════════════════════════════════════════════════════════════
# 2. Parse GSE120575 metadata (same logic as cross_cancer_validation.py)
# ═════════════════════════════════════════════════════════════════════════════
print("\n[2] Parsing GSE120575 metadata...")

with gzip.open(META_FILE, 'rt', encoding='latin-1') as f:
    lines = f.readlines()

# Find header row
header_idx = None
for i, line in enumerate(lines):
    if line.startswith('Sample name'):
        header_idx = i
        break

meta = pd.read_csv(
    META_FILE,
    compression='gzip',
    encoding='latin-1',
    sep='\t',
    skiprows=header_idx,
    nrows=16300,
    usecols=[
        'title',
        'characteristics: patinet ID (Pre=baseline; Post= on treatment)',
        'characteristics: response',
        'characteristics: therapy',
    ],
)
meta.rename(columns={
    'title': 'cell',
    'characteristics: patinet ID (Pre=baseline; Post= on treatment)': 'sample',
    'characteristics: response': 'response',
    'characteristics: therapy': 'therapy',
}, inplace=True)

# Drop any malformed rows
meta = meta[meta['response'].isin(['Responder', 'Non-responder'])].copy()
meta['response_binary'] = (meta['response'] == 'Responder').astype(int)

print(f"  Metadata loaded: {len(meta)} cells")
print(f"    Responders    : {(meta['response'] == 'Responder').sum()}")
print(f"    Non-responders: {(meta['response'] == 'Non-responder').sum()}")
print(f"    Unique samples: {meta['sample'].nunique()}")

# Sample-level response map
sample_meta = meta.drop_duplicates('sample')[
    ['sample', 'response', 'response_binary', 'therapy']
].copy()
print(f"  Sample-level: {len(sample_meta)} samples")

# Cell → sample mapping (use metadata, not expression row 2, to avoid enriched suffixes)
cell_to_sample = dict(zip(meta['cell'], meta['sample']))

# ═════════════════════════════════════════════════════════════════════════════
# 3. Load GSE120575 expression matrix — shared genes only
# ═════════════════════════════════════════════════════════════════════════════
print("\n[3] Loading GSE120575 expression matrix (shared genes only)...")

with gzip.open(EXPR_FILE, 'rt') as f:
    cell_names = f.readline().rstrip().split('\t')[1:]  # skip leading empty field
    _ = f.readline()  # skip sample IDs row

print(f"  Expression file cells: {len(cell_names)}")

# Intersect with metadata cells to ensure 1:1 mapping
valid_cells = [c for c in cell_names if c in cell_to_sample]
print(f"  Cells with metadata: {len(valid_cells)}")

# Build cell index map
name_to_idx = {name: i for i, name in enumerate(cell_names)}
valid_indices = [name_to_idx[c] for c in valid_cells]

# Identify common genes and load only those
genes_expr_set = set()
with gzip.open(EXPR_FILE, 'rt') as f:
    f.readline()
    f.readline()
    for line in f:
        gene = line.rstrip().split('\t')[0]
        genes_expr_set.add(gene)

common_genes = sorted(genes_ref & genes_expr_set)
print(f"  Common genes: {len(common_genes)}")

# Load expression for common genes only (genes × valid_cells)
common_set = set(common_genes)
expr_rows = []
gene_names_loaded = []

with gzip.open(EXPR_FILE, 'rt') as f:
    f.readline()  # skip cell names
    f.readline()  # skip sample IDs
    for line in f:
        parts = line.rstrip().split('\t')
        gene = parts[0]
        if gene in common_set:
            vals = np.array(
                [float(x) if x != '' else 0.0 for x in parts[1:]],
                dtype=np.float32
            )
            # Subset to valid cells
            expr_rows.append(vals[valid_indices])
            gene_names_loaded.append(gene)

# Align gene order to reference
ref = ref[common_genes]  # ensure column order
expr_matrix = np.array(expr_rows)  # common_genes × valid_cells
# Reorder expression rows to match common_genes (sorted) order
gene_order_map = {gene: i for i, gene in enumerate(gene_names_loaded)}
reordered = np.zeros((len(common_genes), expr_matrix.shape[1]), dtype=np.float32)
for i, gene in enumerate(common_genes):
    reordered[i] = expr_matrix[gene_order_map[gene]]
expr_matrix = reordered
gene_names_loaded = common_genes

print(f"  Expression matrix shape: {expr_matrix.shape} (genes × cells)")

# Transpose to cells × genes for cosine_similarity
expr_matrix = expr_matrix.T  # cells × genes
print(f"  Transposed: {expr_matrix.shape} (cells × genes)")

del genes_expr_set, common_set, expr_rows
gc.collect()

# ═════════════════════════════════════════════════════════════════════════════
# 4. Cosine similarity mapping to 51 subtypes
# ═════════════════════════════════════════════════════════════════════════════
print("\n[4] Mapping cells to 51 subtypes (cosine similarity)...")
ref_mat = ref.values  # 51 × common_genes

similarities = cosine_similarity(expr_matrix, ref_mat)  # cells × 51
predicted_idx = np.argmax(similarities, axis=1)
predicted_types = [cell_types[i] for i in predicted_idx]
max_sim = np.max(similarities, axis=1)

print(f"  Mean max similarity: {max_sim.mean():.4f}")
print(f"  Median max similarity: {np.median(max_sim):.4f}")

# Cell type distribution
from collections import Counter
dist = Counter(predicted_types)
print("  Top 10 predicted types:")
for ct, cnt in dist.most_common(10):
    print(f"    {ct}: {cnt} ({cnt / len(predicted_types) * 100:.1f}%)")

del similarities, expr_matrix, ref_mat
gc.collect()

# ═════════════════════════════════════════════════════════════════════════════
# 5. Compute per-sample 51-subtype proportions
# ═════════════════════════════════════════════════════════════════════════════
print("\n[5] Computing per-sample 51-subtype proportions...")

# Map each cell to its sample using metadata
sample_ids = [cell_to_sample[c] for c in valid_cells]

# Build proportion matrix
unique_samples = sorted(set(sample_ids))
type_to_idx = {t: j for j, t in enumerate(cell_types)}
sample_to_idx = {s: i for i, s in enumerate(unique_samples)}

prop_matrix = np.zeros((len(unique_samples), len(cell_types)))
for s, t in zip(sample_ids, predicted_types):
    prop_matrix[sample_to_idx[s], type_to_idx[t]] += 1

# Normalize to proportions
row_sums = prop_matrix.sum(axis=1, keepdims=True)
row_sums[row_sums == 0] = 1
prop_matrix = prop_matrix / row_sums

prop_df = pd.DataFrame(prop_matrix, index=unique_samples, columns=cell_types)
print(f"  Proportion matrix: {prop_df.shape}")
print(f"  Samples: {len(unique_samples)}")

# ═════════════════════════════════════════════════════════════════════════════
# 6. Align to GSE243013 training features
# ═════════════════════════════════════════════════════════════════════════════
print("\n[6] Aligning features to GSE243013 training set...")
train_props = pd.read_csv(
    os.path.join(TOMATO_DATA, 'cell_proportions.csv'),
    index_col='sampleID'
)
train_features = train_props.columns.tolist()
print(f"  Training features: {len(train_features)}")

# Add missing features as 0
for f in train_features:
    if f not in prop_df.columns:
        prop_df[f] = 0.0
prop_df = prop_df[train_features]  # enforce same order

# ═════════════════════════════════════════════════════════════════════════════
# 7. Load training labels and train XGBoost
# ═════════════════════════════════════════════════════════════════════════════
print("\n[7] Training XGBoost on GSE243013...")
train_labels = pd.read_csv(
    os.path.join(TOMATO_DATA, 'sample_labels.csv'),
    index_col='sampleID'
)['response']

common_idx = train_props.index.intersection(train_labels.index)
X_train = train_props.loc[common_idx].values
y_train = train_labels.loc[common_idx].values
print(f"  Training samples: {len(common_idx)}")
print(f"    R={(y_train == 1).sum()}, NR={(y_train == 0).sum()}")

# Standardize
scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)

# XGBoost
model = xgb.XGBClassifier(
    n_estimators=200,
    max_depth=4,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    reg_alpha=0.1,
    reg_lambda=1.0,
    random_state=SEED,
    eval_metric='logloss',
    use_label_encoder=False,
    n_jobs=4
)
model.fit(X_train_s, y_train)

# ═════════════════════════════════════════════════════════════════════════════
# 8. Predict on GSE120575
# ═════════════════════════════════════════════════════════════════════════════
print("\n[8] Predicting on GSE120575...")

# Merge response labels
prop_df = prop_df.merge(
    sample_meta[['sample', 'response', 'response_binary']],
    left_index=True,
    right_on='sample',
    how='left'
)
prop_df.set_index('sample', inplace=True)

# Drop any samples without response info (should not happen)
prop_df = prop_df.dropna(subset=['response_binary']).copy()
prop_df['response_binary'] = prop_df['response_binary'].astype(int)

samples = prop_df.index.tolist()
X_test = prop_df[train_features].values
y_test = prop_df['response_binary'].values

print(f"  Test samples: {len(samples)}")
print(f"    R={(y_test == 1).sum()}, NR={(y_test == 0).sum()}")

X_test_s = scaler.transform(X_test)
y_prob = model.predict_proba(X_test_s)[:, 1]
y_pred = model.predict(X_test_s)

# Metrics
auc = roc_auc_score(y_test, y_prob)
acc = accuracy_score(y_test, y_pred)
cm = confusion_matrix(y_test, y_pred)
tn, fp, fn, tp = cm.ravel()
sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
npv = tn / (tn + fn) if (tn + fn) > 0 else 0

print("\n" + "=" * 60)
print("RESULTS — GSE120575 Full Model Validation")
print("=" * 60)
print(f"  AUC        = {auc:.3f}")
print(f"  ACC        = {acc:.3f} ({(y_pred == y_test).sum()}/{len(y_test)})")
print(f"  Sensitivity = {sensitivity:.3f} ({tp}/{tp + fn})")
print(f"  Specificity = {specificity:.3f} ({tn}/{tn + fp})")
print(f"  PPV         = {ppv:.3f} ({tp}/{tp + fp})")
print(f"  NPV         = {npv:.3f} ({tn}/{tn + fn})")
print(f"  Confusion Matrix:")
print(f"    TN={tn}  FP={fp}")
print(f"    FN={fn}  TP={tp}")

# Per-sample details
print("\nPer-sample predictions (sorted by probability descending):")
results = []
for s, prob, pred, true in zip(samples, y_prob, y_pred, y_test):
    mark = '✓' if pred == true else '✗'
    rl = 'R' if true == 1 else 'NR'
    pl = 'R' if pred == 1 else 'NR'
    print(f"  {mark} {s}: prob={prob:.4f}, pred={pl}, true={rl}")
    results.append({
        'sample': s,
        'true_label': true,
        'true_response': rl,
        'predicted_prob': prob,
        'predicted_label': pred,
        'predicted_response': pl,
        'correct': int(pred == true),
    })

# Save CSV
results_df = pd.DataFrame(results)
results_df = results_df.sort_values('predicted_prob', ascending=False)
results_path = os.path.join(OUT_DIR, 'full_model_validation_results.csv')
results_df.to_csv(results_path, index=False)
print(f"\n  Saved: {results_path}")

# ═════════════════════════════════════════════════════════════════════════════
# 9. Plot ROC curve
# ═════════════════════════════════════════════════════════════════════════════
print("\n[9] Plotting ROC curve...")

fig, ax = plt.subplots(figsize=(7, 7))
fpr, tpr, _ = roc_curve(y_test, y_prob)

# Color-blind-friendly teal
curve_color = '#56B4E9'
ax.plot(fpr, tpr, color=curve_color, lw=3,
        label=f'XGBoost (51 subtypes)  AUC = {auc:.3f}')
ax.fill_between(fpr, tpr, alpha=0.15, color=curve_color)
ax.plot([0, 1], [0, 1], 'k--', lw=1.5, alpha=0.5, label='Chance')

ax.set_xlabel('False Positive Rate', fontsize=14, fontweight='bold')
ax.set_ylabel('True Positive Rate', fontsize=14, fontweight='bold')
ax.set_title(
    'External Validation — GSE120575 Melanoma\n'
    f'(n={len(y_test)}, R={(y_test == 1).sum()}, NR={(y_test == 0).sum()})',
    fontsize=15, fontweight='bold'
)
ax.legend(loc='lower right', fontsize=13)
ax.set_xlim([-0.02, 1.02])
ax.set_ylim([-0.02, 1.02])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(
    os.path.join(OUT_DIR, 'full_model_validation_roc.png'),
    dpi=300, bbox_inches='tight', facecolor='white'
)
fig.savefig(
    os.path.join(OUT_DIR, 'full_model_validation_roc.pdf'),
    dpi=300, bbox_inches='tight', facecolor='white'
)
plt.close()

print(f"  Saved: {os.path.join(OUT_DIR, 'full_model_validation_roc.png')}")
print(f"  Saved: {os.path.join(OUT_DIR, 'full_model_validation_roc.pdf')}")

print("\n" + "=" * 60)
print("Done.")
print("=" * 60)
