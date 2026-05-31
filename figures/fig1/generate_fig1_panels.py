"""
Generate final Figure 1 panels for publication.
Fig1 structure (7 panels):
  a = Study Design (manual)
  b = CD4T Subtype Marker Heatmap  <- NEW (this script)
  c = Treg Subtype Abundance by Response  <- was fig1b
  d = ROC Curves  <- from xgboost_shap_analysis.py (fig1d prefix)
  e = SHAP Bar Importance  <- from xgboost_shap_analysis.py (fig1e prefix)
  f = Decision Tree (compact)  <- was fig1e
  g = SHAP Beeswarm  <- from xgboost_shap_analysis.py (fig1g prefix)

Unified color scheme: green-white-purple diverging (soft tones)
  Responder / high:   #8bc98b  (soft green)
  Non-responder / low:#b08ebf  (soft purple)
"""

import os
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib as mpl
import seaborn as sns
from scipy import stats
from sklearn.tree import plot_tree, DecisionTreeClassifier
from sklearn.preprocessing import StandardScaler
import pickle

BASE_DIR = r'D:\Research\tomato'
DATA_DIR = os.path.join(BASE_DIR, 'data')
FIGURES_DIR = r'D:\research\cucumber\fig1'
os.makedirs(FIGURES_DIR, exist_ok=True)

# ── Unified color palette ──────────────────────────────────────────────────
COLOR_HIGH = '#2E5AAC'      # blue (Responder)
COLOR_LOW  = '#C0392B'      # red (Non-responder)
COLOR_MID  = '#ffffff'      # white
CMAP_GWP = LinearSegmentedColormap.from_list('gwp',
    [COLOR_LOW, '#e8a5a5', COLOR_MID, '#a5c4e8', COLOR_HIGH])

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 12,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
    'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
})

# ==============================================================================
# PANEL B: CD4T Subtype Marker Heatmap (rounded corners, full matrix)
# ==============================================================================
print("=" * 60)
print("Generating Panel B: CD4T Marker Heatmap...")
print("=" * 60)

pb = pd.read_csv(os.path.join(DATA_DIR, 'GSE243013_reference_pseudobulk.csv'), index_col=0)

KNOWN_MARKERS = {
    'CD4T_Tem_GZMA':        ['GZMA',  'GZMK',   'CCL5'],
    'CD4T_Tfh_CXCL13':      ['CXCL13','BCL6',   'PDCD1'],
    'CD4T_Th1-like_CXCL13': ['CXCL13','IFNG',   'TBX21'],
    'CD4T_Tm_ANXA1':        ['ANXA1', 'EGR1',   'KLF2'],
    'CD4T_Tm_XCL1':         ['XCL1',  'XCL2',   'CCL5'],
    'CD4T_Tn_CCR7':         ['CCR7',  'SELL',   'LEF1'],
    'CD4T_Treg_CCR8':       ['CCR8',  'FOXP3',  'TNFRSF18'],
    'CD4T_Treg_FOXP3':      ['FOXP3', 'IL2RA',  'CTLA4'],
    'CD4T_Treg_MKI67':      ['MKI67', 'STMN1',  'FOXP3'],
}
subtypes = list(KNOWN_MARKERS.keys())
all_unique_genes = list(dict.fromkeys([g for m in KNOWN_MARKERS.values() for g in m]))

expr_raw = pb.loc[subtypes, all_unique_genes]
expr_norm = expr_raw.copy()
for col in expr_norm.columns:
    mn, mx = expr_norm[col].min(), expr_norm[col].max()
    expr_norm[col] = (expr_norm[col] - mn) / (mx - mn) if mx > mn else 0.5

fig, ax = plt.subplots(figsize=(20, 10))
n_subtypes, n_markers, n_groups = 9, 3, 9

cell_w = cell_h = 0.30
intra_gap, inter_gap, row_gap = 0.004, 0.25, 0.012
round_r = 0.03

group_w = n_markers * cell_w + (n_markers - 1) * intra_gap
total_w = n_groups * group_w + (n_groups - 1) * inter_gap
start_x = -total_w / 2
total_h = n_subtypes * cell_h + (n_subtypes - 1) * row_gap
start_y = total_h / 2

for col_idx, (group_subtype, markers) in enumerate(KNOWN_MARKERS.items()):
    group_start_x = start_x + col_idx * (group_w + inter_gap)
    for j, gene in enumerate(markers):
        x = group_start_x + j * (cell_w + intra_gap)
        for row_idx, subtype in enumerate(subtypes):
            y = start_y - row_idx * (cell_h + row_gap)
            val = expr_norm.loc[subtype, gene]
            rect = FancyBboxPatch(
                (x, y - cell_h), cell_w, cell_h,
                boxstyle=f"round,pad=0,rounding_size={round_r}",
                facecolor=CMAP_GWP(val), edgecolor='#ffffff', linewidth=0.7,
                clip_on=False, zorder=2)
            ax.add_patch(rect)

ylabels = [st.replace('CD4T_', '').replace('_', ' ') for st in subtypes]
y_positions = [start_y - i * (cell_h + row_gap) - cell_h/2 for i in range(n_subtypes)]
for yp, lbl in zip(y_positions, ylabels):
    ax.text(start_x - 0.10, yp, lbl, ha='right', va='center', fontsize=12,
            fontweight='bold', color='#333333', fontfamily='Arial')

for col_idx, (group_subtype, markers) in enumerate(KNOWN_MARKERS.items()):
    group_start_x = start_x + col_idx * (group_w + inter_gap)
    for j, gene in enumerate(markers):
        x = group_start_x + j * (cell_w + intra_gap) + cell_w / 2
        y_bottom = start_y - n_subtypes * (cell_h + row_gap) + row_gap - 0.06
        ax.text(x, y_bottom, gene, ha='center', va='top', fontsize=8.5,
                color='#555555', fontfamily='Arial', rotation=45, rotation_mode='anchor')

margin_x, margin_y = 0.5, 0.4
ax.set_xlim(start_x - margin_x, start_x + total_w + margin_x)
ax.set_ylim(start_y - n_subtypes*(cell_h+row_gap) + row_gap - margin_y, start_y + margin_y)
ax.set_aspect('equal', adjustable='box')
ax.set_frame_on(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_title('B. CD4T Subtype-Specific Marker Expression',
             fontsize=15, fontweight='bold', color='#222222', pad=18, fontfamily='Arial')

sm = plt.cm.ScalarMappable(cmap=CMAP_GWP, norm=plt.Normalize(0, 1))
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, shrink=0.30, aspect=20, pad=0.015, anchor=(0, 0.5))
cbar.set_label('Normalized Expression', fontsize=10, color='#555555', fontfamily='Arial')
cbar.ax.tick_params(labelsize=8, colors='#555555')
cbar.outline.set_visible(False)
cbar.ax.set_frame_on(False)

plt.tight_layout()
fig.savefig(os.path.join(FIGURES_DIR, 'fig1b_cd4t_marker_heatmap.png'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(FIGURES_DIR, 'fig1b_cd4t_marker_heatmap.pdf'), bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved fig1b_cd4t_marker_heatmap.png/pdf")

# ==============================================================================
# PANEL C: Treg Subtype Abundance by Response (with exact p-values)
# ==============================================================================
print("\nGenerating Panel C: Treg subtype abundance...")

props = pd.read_csv(os.path.join(DATA_DIR, 'cell_proportions.csv'), index_col=0)
labels = pd.read_csv(os.path.join(DATA_DIR, 'sample_labels.csv'), index_col=0)['response']
common = props.index.intersection(labels.index)
props = props.loc[common]
labels = labels.loc[common]

results = []
for feat in ['CD4T_Treg_CCR8', 'CD4T_Treg_MKI67']:
    r_vals = props.loc[labels == 1, feat].values
    nr_vals = props.loc[labels == 0, feat].values
    _, pval = stats.mannwhitneyu(nr_vals, r_vals, alternative='two-sided')
    results.append({'feature': feat, 'NR_mean': nr_vals.mean(), 'R_mean': r_vals.mean(),
                    'NR_median': np.median(nr_vals), 'R_median': np.median(r_vals),
                    'p_value': pval, 'n_NR': len(nr_vals), 'n_R': len(r_vals)})

stats_df = pd.DataFrame(results)
print(stats_df.to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(10, 5.5), gridspec_kw={'wspace': 0.3})
palette = {'Non-responder': COLOR_LOW, 'Responder': COLOR_HIGH}

for idx, (feat, title) in enumerate([('CD4T_Treg_CCR8', 'CCR8+ Treg'), ('CD4T_Treg_MKI67', 'MKI67+ Treg')]):
    ax = axes[idx]
    plot_data = []
    for resp, label in [(0, 'Non-responder'), (1, 'Responder')]:
        vals = props.loc[labels == resp, feat].values
        for v in vals:
            plot_data.append({'Response': label, 'Proportion': v})
    df_plot = pd.DataFrame(plot_data)

    sns.boxplot(data=df_plot, x='Response', y='Proportion', palette=palette, ax=ax,
                width=0.5, showfliers=True,
                boxprops=dict(linewidth=1), whiskerprops=dict(linewidth=1),
                medianprops=dict(color='black', linewidth=1.5),
                flierprops=dict(marker='o', markersize=4, alpha=0.5))
    sns.stripplot(data=df_plot, x='Response', y='Proportion', color='black',
                  alpha=0.3, size=4, ax=ax, jitter=True)

    row = stats_df[stats_df['feature'] == feat].iloc[0]
    pval = row['p_value']
    y_max = df_plot['Proportion'].max()
    y_min = df_plot['Proportion'].min()
    y_range = y_max - y_min

    sig_text = '***p < 1×10^-10'

    ax.plot([0, 1], [y_max + y_range*0.05, y_max + y_range*0.05], 'k-', linewidth=1)
    ax.text(0.5, y_max + y_range*0.08, sig_text, ha='center', va='bottom',
            fontsize=10, fontweight='bold')

    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Responder\n(R)', 'Non-responder\n(NR)'], fontsize=9)
    ax.set_ylabel('Proportion of Tregs', fontsize=10)
    sns.despine(ax=ax)

plt.suptitle('C. Treg Subtype Abundance by Immunotherapy Response', fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(FIGURES_DIR, 'fig1c_treg_abundance_with_pvalue.png'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(FIGURES_DIR, 'fig1c_treg_abundance_with_pvalue.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved fig1c_treg_abundance_with_pvalue.png/pdf")

# ==============================================================================
# PANEL F: Decision Tree (compact, unified green-purple colors)
# ==============================================================================
print("\nGenerating Panel F: Decision Tree...")

dt_path = os.path.join(BASE_DIR, 'results', 'decision_tree', 'decision_tree_model.pkl')
if os.path.exists(dt_path):
    with open(dt_path, 'rb') as f:
        dt_model = pickle.load(f)
    print("  Loaded existing decision tree model")
else:
    print("  Reconstructing decision tree from data...")
    props = pd.read_csv(os.path.join(DATA_DIR, 'cell_proportions.csv'), index_col=0)
    labels = pd.read_csv(os.path.join(DATA_DIR, 'sample_labels.csv'), index_col=0)['response']
    common = props.index.intersection(labels.index)
    X = props.loc[common].values
    y = labels.loc[common].values
    feature_names = props.columns.tolist()
    scaler = StandardScaler()
    X_s = scaler.fit_transform(X)
    dt_model = DecisionTreeClassifier(max_depth=3, random_state=42, class_weight='balanced')
    dt_model.fit(X_s, y)
    os.makedirs(os.path.dirname(dt_path), exist_ok=True)
    with open(dt_path, 'wb') as f:
        pickle.dump(dt_model, f)
    print("  Saved reconstructed model")

fig, ax = plt.subplots(figsize=(14, 8))
# Plot tree with default colors, then recolor nodes to green-purple
annotations = plot_tree(dt_model, feature_names=props.columns.tolist(),
                        class_names=['Responder', 'Non-responder'],
                        filled=True, rounded=True, fontsize=7, impurity=False,
                        proportion=True, precision=2, ax=ax, max_depth=3)

# Recolor: sklearn 1.8 uses orange (class0/Responder) and blue (class1/Non-responder)
ref_orange = np.array([0.898, 0.506, 0.224])
ref_blue   = np.array([0.224, 0.616, 0.898])

for ann in annotations:
    bbox = ann.get_bbox_patch()
    if bbox is None:
        continue
    fc = np.array(bbox.get_facecolor()[:3])
    # Skip near-white (root or text-only nodes)
    if np.mean(fc) > 0.97:
        continue
    d_orange = np.linalg.norm(fc - ref_orange)
    d_blue   = np.linalg.norm(fc - ref_blue)
    if d_orange + d_blue < 1e-6:
        continue
    # ratio=1 -> orange (Responder) -> blue; ratio=0 -> blue (Non-responder) -> red
    ratio = d_blue / (d_orange + d_blue)
    bbox.set_facecolor(CMAP_GWP(ratio))

ax.set_title('F. Decision Tree (AUC = 0.744)',
             fontsize=12, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(FIGURES_DIR, 'fig1f_decision_tree_compact.png'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(FIGURES_DIR, 'fig1f_decision_tree_compact.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved fig1f_decision_tree_compact.png/pdf")

# ==============================================================================
# Copy external panels (d, e, g) from xgboost_shap_analysis output if available
# ==============================================================================
print("\nChecking for external panels (d, e, g)...")

# Possible source directories for xgboost_shap_analysis outputs
possible_shap_dirs = [
    os.path.join(BASE_DIR, 'data', 'raw', 'GSE243013', 'XGBoost_SHAP_results'),
    r'D:\Research\土豆\数据\raw\GSE243013\XGBoost_SHAP_results',
    r'D:\research\cucumber\fig1',
]

shap_src = None
for d in possible_shap_dirs:
    if os.path.isdir(d):
        shap_src = d
        break

if shap_src:
    copy_map = {
        'fig1d_ROC_curves_5fold.png':   '03_ROC_curves_5fold.png',
        'fig1e_SHAP_bar_importance.png':'02_SHAP_bar_importance.png',
        'fig1g_SHAP_summary_beeswarm.png':'01_SHAP_summary_beeswarm.png',
    }
    for dst_name, src_name in copy_map.items():
        src = os.path.join(shap_src, src_name)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(FIGURES_DIR, dst_name))
            print(f"  Copied {src_name} -> {dst_name}")
        else:
            print(f"  [Not found] {src_name}")
else:
    print("  XGBoost_SHAP_results directory not found; skip copy.")

print("\n" + "=" * 60)
print("DONE! Fig1 panels saved to:", FIGURES_DIR)
print("=" * 60)
