#!/usr/bin/env python3
"""
Recalculate two-feature logistic regression AUC for GSE120575
using the updated FOXP3+IL2RA+CTLA4>=2 Treg definition.
"""
import os, gzip, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

# ─── Paths ───
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, '..', 'fig1', 'cross_cancer', 'data'))
OUTDIR = SCRIPT_DIR
os.makedirs(OUTDIR, exist_ok=True)

TARGET_GENES = ['FOXP3', 'IL2RA', 'CTLA4', 'CCR8', 'MKI67']
POS_THRESHOLD = 0.0
META_FILE = os.path.join(DATA_DIR, 'metadata.txt.gz')
EXPR_FILE = os.path.join(DATA_DIR, 'expression.txt.gz')

# ─── Load metadata ───
with gzip.open(META_FILE, 'rt', encoding='latin-1') as f:
    lines = f.readlines()
header_idx = next(i for i, line in enumerate(lines) if line.startswith('Sample name'))

meta = pd.read_csv(META_FILE, compression='gzip', encoding='latin-1', sep='\t',
                    skiprows=header_idx,
                    usecols=['title',
                             'characteristics: patinet ID (Pre=baseline; Post= on treatment)',
                             'characteristics: response', 'characteristics: therapy'])
meta.rename(columns={'title': 'cell',
                     'characteristics: patinet ID (Pre=baseline; Post= on treatment)': 'sample',
                     'characteristics: response': 'response',
                     'characteristics: therapy': 'therapy'}, inplace=True)
meta = meta[meta['response'].isin(['Responder', 'Non-responder'])].copy()
meta['response_binary'] = (meta['response'] == 'Responder').astype(int)
sample_meta = meta.drop_duplicates('sample')[['sample', 'response', 'response_binary']].copy()

# ─── Load expression ───
with gzip.open(EXPR_FILE, 'rt') as f:
    cell_names = f.readline().rstrip().split('\t')[1:]
    sample_ids = f.readline().rstrip().split('\t')[1:]
    gene_expr = {}
    for line in f:
        parts = line.rstrip().split('\t')
        gene = parts[0]
        if gene in TARGET_GENES:
            vals = np.array([float(x) if x != '' else 0.0 for x in parts[1:]], dtype=np.float32)
            gene_expr[gene] = vals

cell_df = pd.DataFrame({'cell': cell_names, 'sample': sample_ids,
                        'FOXP3': gene_expr['FOXP3'], 'IL2RA': gene_expr['IL2RA'],
                        'CTLA4': gene_expr['CTLA4'], 'CCR8': gene_expr['CCR8'],
                        'MKI67': gene_expr['MKI67']})
cell_df = cell_df[cell_df['cell'].isin(meta['cell'])].copy()

# Map enriched samples back to parent biopsy
cell_df['sample'] = cell_df['sample'].str.replace(r'_T_enriched$|_myeloid_enriched$', '', regex=True)

# Treg definition: FOXP3 + IL2RA + CTLA4 >= 2
treg_score = ((cell_df['FOXP3'] > POS_THRESHOLD).astype(int) +
              (cell_df['IL2RA'] > POS_THRESHOLD).astype(int) +
              (cell_df['CTLA4'] > POS_THRESHOLD).astype(int))
cell_df['is_treg'] = treg_score >= 2
cell_df['is_CCR8_Treg'] = cell_df['is_treg'] & (cell_df['CCR8'] > POS_THRESHOLD)
cell_df['is_MKI67_Treg'] = cell_df['is_treg'] & (cell_df['MKI67'] > POS_THRESHOLD)

sample_stats = []
for sample, grp in cell_df.groupby('sample'):
    n_total = len(grp)
    sample_stats.append({'sample': sample, 'n_cells': n_total,
                         'prop_ccr8_treg': grp['is_CCR8_Treg'].sum() / n_total,
                         'prop_mki67_treg': grp['is_MKI67_Treg'].sum() / n_total})
sample_df = pd.DataFrame(sample_stats).merge(sample_meta, on='sample')

# Save updated sample stats
sample_df.to_csv(os.path.join(DATA_DIR, 'cross_cancer_sample_stats.csv'), index=False)
print(f"Saved updated cross_cancer_sample_stats.csv ({len(sample_df)} samples)")

# ─── Two-feature logistic regression ───
X = sample_df[['prop_ccr8_treg', 'prop_mki67_treg']].values
y = sample_df['response_binary'].values

scaler = StandardScaler()
X_s = scaler.fit_transform(X)

model = LogisticRegression(C=1.0, max_iter=5000, class_weight='balanced', random_state=42)
model.fit(X_s, y)
y_prob = model.predict_proba(X_s)[:, 1]
auc = roc_auc_score(y, y_prob)

print(f"\nTwo-feature Logistic Regression (GSE120575)")
print(f"  AUC = {auc:.3f}")
print(f"  Coefficients: CCR8={model.coef_[0][0]:.3f}, MKI67={model.coef_[0][1]:.3f}")
print(f"  Intercept: {model.intercept_[0]:.3f}")

# ─── Bootstrap 95% CI ───
rng = np.random.default_rng(42)
n_boot = 10000
boot_aucs = []

for i in range(n_boot):
    idx = rng.choice(len(y), size=len(y), replace=True)
    y_boot = y[idx]
    # Ensure both classes present
    if len(np.unique(y_boot)) < 2:
        continue
    X_boot = X_s[idx]
    try:
        m_boot = LogisticRegression(C=1.0, max_iter=5000, class_weight='balanced', random_state=42)
        m_boot.fit(X_boot, y_boot)
        yp_boot = m_boot.predict_proba(X_boot)[:, 1]
        a = roc_auc_score(y_boot, yp_boot)
        boot_aucs.append(a)
    except:
        continue

boot_aucs = np.array(boot_aucs)
ci_low = np.percentile(boot_aucs, 2.5)
ci_high = np.percentile(boot_aucs, 97.5)

print(f"\nBootstrap 95% CI: [{ci_low:.3f}, {ci_high:.3f}]")
print(f"Bootstrap mean ± SD: {boot_aucs.mean():.3f} ± {boot_aucs.std():.3f}")
print(f"Valid iterations: {len(boot_aucs)} / {n_boot}")

# ─── Plot ROC ───
fpr, tpr, _ = roc_curve(y, y_prob)
fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(fpr, tpr, color='#C0392B', lw=2.5,
        label=f'Two-feature LR  AUC = {auc:.3f}\n(95% CI [{ci_low:.3f}, {ci_high:.3f}])')
ax.fill_between(fpr, tpr, alpha=0.15, color='#C0392B')
ax.plot([0, 1], [0, 1], 'k--', lw=1.2, alpha=0.4, label='Chance')
ax.set_xlabel('False Positive Rate', fontsize=12, fontweight='bold')
ax.set_ylabel('True Positive Rate', fontsize=12, fontweight='bold')
ax.set_title('GSE120575 Two-Feature Logistic Regression\n(CCR8+ Treg + MKI67+ Treg proportions)',
             fontsize=13, fontweight='bold')
ax.set_xlim([-0.02, 1.02])
ax.set_ylim([-0.02, 1.02])
ax.set_aspect('equal')
ax.legend(loc='lower right', fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, alpha=0.15)
plt.tight_layout()

for ext in ['png', 'pdf']:
    path = os.path.join(OUTDIR, f'gse120575_two_feature_roc.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {path}")
plt.close()

# ─── Print sample details ───
print(f"\nn = {len(sample_df)} samples (R={(y==1).sum()}, NR={(y==0).sum()})")
print(f"CCR8+ Treg:  p = {stats.mannwhitneyu(sample_df.loc[y==1, 'prop_ccr8_treg'], sample_df.loc[y==0, 'prop_ccr8_treg'], alternative='two-sided')[1]:.4f}")
print(f"MKI67+ Treg: p = {stats.mannwhitneyu(sample_df.loc[y==1, 'prop_mki67_treg'], sample_df.loc[y==0, 'prop_mki67_treg'], alternative='two-sided')[1]:.4f}")
