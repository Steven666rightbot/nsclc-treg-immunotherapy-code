#!/usr/bin/env python3
"""
Supplementary Figure — Cross-Cancer Validation ROC: GSE120575 Melanoma

External validation of the 51-subtype XGBoost model on GSE120575
(Sade-Feldman et al. Cell 2018, n=48 melanoma biopsies, CD45+ Smart-seq2).

AUC = 0.666 — moderate cross-cancer generalizability.
Trained on GSE243013 NSCLC, tested on GSE120575 melanoma.

Output: suppl_cross_cancer_roc.png / .pdf
"""

import os, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve

warnings.filterwarnings('ignore')

OUT_DIR = r'D:\research\cucumber\fig1'
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'axes.linewidth': 0.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# ─── Load GSE120575 predictions ────────────────────────────────────────────────
results = pd.read_csv(
    r'D:\research\cucumber\fig1\cross_cancer\full_model_validation_results.csv'
)

y_true = results['true_label'].values
y_prob = results['predicted_prob'].values

auc = roc_auc_score(y_true, y_prob)
fpr, tpr, _ = roc_curve(y_true, y_prob)

n = len(y_true)
n_r = (y_true == 1).sum()
n_nr = n - n_r

# ─── Load sample-level stats for context ──────────────────────────────────────
stats = pd.read_csv(
    r'D:\research\cucumber\fig1\cross_cancer\cross_cancer_sample_stats.csv'
)
ccr8_med = stats['prop_ccr8_treg'].median()
mki67_med = stats['prop_mki67_treg'].median()

print(f"GSE120575 Melanoma External Validation")
print(f"  Samples: {n} (R={n_r}, NR={n_nr})")
print(f"  AUC = {auc:.4f}")
print(f"  CCR8+ Treg median = {ccr8_med:.4f}")
print(f"  MKI67+ Treg median = {mki67_med:.4f}")

# ─── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6.5, 6.5))

# ROC curve
curve_color = '#B83A3A'  # red accent
ax.plot(fpr, tpr, color=curve_color, lw=2.5,
        label=f'XGBoost (51 subtypes)  AUC = {auc:.3f}')
ax.fill_between(fpr, tpr, alpha=0.15, color=curve_color)

# Chance diagonal
ax.plot([0, 1], [0, 1], 'k--', lw=1.2, alpha=0.4, label='Chance')

ax.set_xlabel('False Positive Rate', fontweight='bold')
ax.set_ylabel('True Positive Rate', fontweight='bold')
ax.set_title('Cross-Cancer Validation — GSE120575 Melanoma', fontweight='bold', pad=12)
ax.set_xlim([-0.02, 1.02])
ax.set_ylim([-0.02, 1.02])
ax.set_aspect('equal')
ax.legend(loc='lower right', fontsize=10, framealpha=0.9, edgecolor='#cccccc')

# Info box
info_text = (
    f'Training: GSE243013 NSCLC (n=242, 51 subtypes)\n'
    f'Testing: GSE120575 Melanoma (n={n}, R={n_r}, NR={n_nr})\n'
    f'Platform: Smart-seq2 | CD45+ sorted\n'
    f'CCR8+ Treg: median={ccr8_med:.1%}\n'
    f'MKI67+ Treg: median={mki67_med:.2%}'
)
ax.text(0.98, 0.02, info_text, transform=ax.transAxes,
        ha='right', va='bottom', fontsize=8, color='#555555',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='#dddddd', alpha=0.9))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(labelsize=10)
ax.grid(True, alpha=0.15)

plt.tight_layout()

for ext in ['png', 'pdf']:
    path = os.path.join(OUT_DIR, f'suppl_cross_cancer_roc.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {path}")

plt.close()
print("\nDone!")
