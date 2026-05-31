#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Publication-quality SHAP figure
Style inspired by: XGBoost-SHAP composite plots
Key technique: simple_beeswarm histogram-based layering
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
import xgboost as xgb
import shap

warnings.filterwarnings('ignore')

# ===================== Config =====================
FIG_DIR = r'D:\Research\tomato\figures'
DATA_DIR = r'D:\Research\tomato\data'
RAW_DIR = r'D:\Research\土豆\data\raw\GSE243013'
CACHE_DIR = os.path.join(DATA_DIR, 'shap_cache')
os.makedirs(CACHE_DIR, exist_ok=True)

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

# Feature groups and colors (inspired by WeChat article style)
FEATURE_GROUPS = {
    'Treg': ['CD4T_Treg_CCR8', 'CD4T_Treg_MKI67', 'CD4T_Treg_FOXP3'],
    'Macrophage': ['Mφ_CXCL10', 'Mφ_ISG15', 'Mφ_S100A10', 'Mφ_DNAJB1', 'Mφ_MMP9',
                   'Mφ_MARCO', 'Mφ_VCAN', 'Mφ_FCGR3A', 'Mφ_MKI67', 'Mφ_FOLR2', 'Mφ_CXCL2'],
    'CD8 T': ['CD8T_Tex_CXCL13', 'CD8T_Trm_ZNF683', 'CD8T_terminal_Tex_LAYN',
              'CD8T_prf_MKI67', 'CD8T_Tem_GZMK+GZMH+', 'CD8T_MAIT_KLRB1',
              'CD8T_NK-like_FGFBP2', 'CD8T_Tm_IL7R', 'CD8T_ISG15', 'CD8T_Tem_GZMK+NR4A1+'],
    'CD4 T': ['CD4T_Tfh_CXCL13', 'CD4T_Th1-like_CXCL13', 'CD4T_Tm_XCL1',
              'CD4T_Tm_ANXA1', 'CD4T_Tem_GZMA', 'CD4T_Tn_CCR7'],
    'B cell': ['Bm_CD74', 'Bm_FCRL4', 'Bm_PDE4D', 'Bm_TNF', 'Bm_MT2A',
               'Bm_TNFSF9', 'Plasma_cell', 'B_prf_MKI67', 'Bn_TCL1A'],
    'NK': ['NK_CD16low_GZMK', 'NK_CD16hi_FGFBP2'],
    'DC': ['cDC1_CLEC9A', 'cDC2_CD1C', 'pDC_LILRA4', 'mDC_LAMP3'],
    'Other': ['T_gdT_TRDV2', 'T_gdT_TRDV1', 'ILC3_KIT', 'Neu_FCGR3B', 'Mast cell', 'GCB_RGS13'],
}

# Color scheme: each group has a base color + feature shades
GROUP_COLOR_SCHEME = {
    'Treg':     {'base': '#c62828', 'shades': ['#b71c1c', '#c62828', '#e53935', '#ef5350', '#e57373']},
    'Macrophage': {'base': '#1565c0', 'shades': ['#0d47a1', '#1565c0', '#1976d2', '#1e88e5', '#42a5f5', '#64b5f6', '#90caf9', '#bbdefb', '#e3f2fd', '#0d47a1', '#1565c0']},
    'CD8 T':    {'base': '#2e7d32', 'shades': ['#1b5e20', '#2e7d32', '#388e3c', '#43a047', '#66bb6a', '#81c784', '#a5d6a7', '#c8e6c9', '#e8f5e9', '#1b5e20']},
    'CD4 T':    {'base': '#6a1b9a', 'shades': ['#4a148c', '#6a1b9a', '#7b1fa2', '#8e24aa', '#ab47bc', '#ce93d8']},
    'B cell':   {'base': '#e65100', 'shades': ['#bf360c', '#e65100', '#ef6c00', '#f57c00', '#fb8c00', '#ff9800', '#ffb74d', '#ffcc80', '#ffe0b2']},
    'NK':       {'base': '#00695c', 'shades': ['#004d40', '#00695c', '#00796b', '#00897b']},
    'DC':       {'base': '#f9a825', 'shades': ['#f57f17', '#f9a825', '#fbc02d', '#fdd835']},
    'Other':    {'base': '#757575', 'shades': ['#424242', '#616161', '#757575', '#9e9e9e', '#bdbdbd', '#e0e0e0']},
}

# Unified viridis-based diverging colormap for whole figure
VIRIDIS_NEG = plt.cm.viridis(0.08)   # deep purple
VIRIDIS_MID = '#f0f0f0'              # light gray
VIRIDIS_POS = plt.cm.viridis(0.92)   # yellow-green
VIRIDIS_CMAP = LinearSegmentedColormap.from_list(
    'viridis_div', [VIRIDIS_NEG, VIRIDIS_MID, VIRIDIS_POS]
)


def build_feature_config():
    """Build feature->color and feature->group mapping."""
    feature_color = {}
    feature_group = {}
    group_base_color = {}
    for group, features in FEATURE_GROUPS.items():
        scheme = GROUP_COLOR_SCHEME[group]
        group_base_color[group] = scheme['base']
        for i, feat in enumerate(features):
            feature_group[feat] = group
            shade_idx = min(i, len(scheme['shades']) - 1)
            feature_color[feat] = scheme['shades'][shade_idx]
    return feature_color, feature_group, group_base_color


FEATURE_COLOR, FEATURE_GROUP, GROUP_BASE_COLOR = build_feature_config()


def simple_beeswarm(y_values, nbins=40, width=0.35):
    """
    Histogram-based beeswarm layering.
    Distribute points along y-axis within a band to avoid overlap.
    """
    y_values = np.asarray(y_values, dtype=float)
    if len(y_values) == 0:
        return np.array([])
    hist_range = (np.min(y_values), np.max(y_values))
    if hist_range[0] == hist_range[1]:
        hist_range = (hist_range[0] - 0.1, hist_range[1] + 0.1)
    counts, edges = np.histogram(y_values, bins=nbins, range=hist_range)
    bin_indices = np.digitize(y_values, edges) - 1
    bin_indices = np.clip(bin_indices, 0, len(counts) - 1)

    offsets = np.zeros(len(y_values))
    for bin_idx in range(len(counts)):
        mask = bin_indices == bin_idx
        n_in_bin = np.sum(mask)
        if n_in_bin > 0:
            # Distribute evenly within [-width/2, +width/2]
            ys = np.linspace(-width / 2, width / 2, n_in_bin)
            offsets[mask] = ys
    return offsets


def compute_or_load_shap():
    """Compute SHAP or load from cache."""
    cache_shap = os.path.join(CACHE_DIR, 'shap_matrix.npy')
    cache_X = os.path.join(CACHE_DIR, 'X_all_test.npy')
    cache_labels = os.path.join(CACHE_DIR, 'sample_labels.csv')
    cache_features = os.path.join(CACHE_DIR, 'feature_names.txt')
    cache_props = os.path.join(CACHE_DIR, 'cell_proportions.csv')

    if all(os.path.exists(f) for f in [cache_shap, cache_X, cache_labels, cache_features, cache_props]):
        print("[Cache] Loading SHAP values...")
        shap_matrix = np.load(cache_shap)
        X_all_test = np.load(cache_X)
        sample_labels = pd.read_csv(cache_labels, index_col=0)['response']
        with open(cache_features, 'r', encoding='utf-8') as f:
            feature_names = [line.strip() for line in f]
        cell_props = pd.read_csv(cache_props, index_col=0)
        return shap_matrix, X_all_test, sample_labels, feature_names, cell_props

    print("[Compute] Loading local data + Running XGBoost + SHAP (5-fold CV)...")
    cell_props = pd.read_csv(os.path.join(DATA_DIR, 'cell_proportions.csv'), index_col=0)
    sample_labels = pd.read_csv(os.path.join(DATA_DIR, 'sample_labels.csv'), index_col=0)['response']

    common = cell_props.index.intersection(sample_labels.index)
    cell_props = cell_props.loc[common]
    sample_labels = sample_labels.loc[common]

    X = cell_props.values
    y = sample_labels.values
    feature_names = cell_props.columns.tolist()

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    all_shap = []
    all_X_test = []
    models = []

    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), 1):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)

        model = xgb.XGBClassifier(
            n_estimators=200, max_depth=4, learning_rate=0.05,
            subsample=0.8, colsample_bytree=0.8,
            reg_alpha=0.1, reg_lambda=1.0, random_state=42,
            eval_metric='logloss', n_jobs=4
        )
        model.fit(X_train_s, y_train)
        explainer = shap.TreeExplainer(model)
        shap_vals = explainer.shap_values(X_test_s)
        all_shap.append(shap_vals)
        all_X_test.append(X_test_s)
        models.append(scaler)

        y_prob = model.predict_proba(X_test_s)[:, 1]
        auc = roc_auc_score(y_test, y_prob)
        print(f"  Fold {fold}: AUC = {auc:.4f}")

    shap_matrix = np.vstack(all_shap)
    X_all_test = np.vstack(all_X_test)

    np.save(cache_shap, shap_matrix)
    np.save(cache_X, X_all_test)
    sample_labels.to_frame('response').to_csv(cache_labels)
    with open(cache_features, 'w', encoding='utf-8') as f:
        for fn in feature_names:
            f.write(fn + '\n')
    cell_props.to_csv(cache_props)

    print(f"  Done. SHAP matrix: {shap_matrix.shape}")
    return shap_matrix, X_all_test, sample_labels, feature_names, cell_props


def plot_custom_beeswarm(ax, shap_matrix, X_all_test, feature_names, max_display=20):
    """
    Custom beeswarm with:
    - Histogram-based layering (simple_beeswarm)
    - Feature-value color gradient (blue-white-red)
    - Group color indicators on y-axis labels
    """
    mean_abs = np.abs(shap_matrix).mean(axis=0)
    top_idx = np.argsort(mean_abs)[::-1][:max_display]

    # Custom colormap
    cmap = LinearSegmentedColormap.from_list('shap_bw', BEESWARM_GRADIENT)

    # Plot each feature
    for rank, fi in enumerate(top_idx):
        y_base = max_display - 1 - rank
        svals = shap_matrix[:, fi]
        fvals = X_all_test[:, fi]

        # Sort by SHAP value
        order = np.argsort(svals)
        svals_sorted = svals[order]
        fvals_sorted = fvals[order]

        # Histogram-based layering
        y_offsets = simple_beeswarm(svals_sorted, nbins=30, width=0.65)
        y_positions = y_base + y_offsets

        # Color by feature value
        vmin, vmax = np.percentile(fvals_sorted, [2, 98])
        if vmin == vmax:
            vmin -= 0.1
            vmax += 0.1
        norm = plt.Normalize(vmin, vmax)
        point_colors = cmap(norm(fvals_sorted))

        ax.scatter(svals_sorted, y_positions, c=point_colors, s=12, alpha=0.7,
                   edgecolors='none', rasterized=True, zorder=2)

    # Y-axis with group color indicators
    labels = [feature_names[i] for i in top_idx[::-1]]
    ax.set_yticks(range(max_display))
    ax.set_yticklabels(labels, fontsize=9)

    # Add colored dots next to labels
    for i, label in enumerate(labels):
        group = FEATURE_GROUP.get(label, 'Other')
        color = GROUP_BASE_COLOR.get(group, '#757575')
        ax.scatter(-0.08, i, c=color, s=60, transform=ax.get_yaxis_transform(),
                   clip_on=False, zorder=10, marker='s', edgecolors='white', linewidth=0.5)

    ax.axvline(x=0, color='#333333', linewidth=0.8, linestyle='-', zorder=1)
    ax.set_xlabel('SHAP value (impact on model output)', fontsize=11)
    ax.set_title('(a) SHAP Beeswarm Plot', fontsize=12, fontweight='bold', loc='left')
    ax.set_xlim(shap_matrix[:, top_idx].min() - 0.15, shap_matrix[:, top_idx].max() + 0.15)
    ax.set_ylim(-0.8, max_display - 0.2)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-2.5, 2.5))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.025, pad=0.015, aspect=35)
    cbar.set_label('Feature value (z-score)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    return labels


def plot_mean_shap_bar(ax, shap_importance_df, max_display=20):
    """Horizontal bar chart of mean |SHAP| with viridis continuous colors."""
    top = shap_importance_df.head(max_display).copy()
    # Sort ascending for barh (top feature at top)
    top = top.iloc[::-1]
    
    # Viridis continuous gradient: high importance = yellow-green, low = blue-purple
    viridis_colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(top)))

    y_pos = np.arange(len(top))
    bars = ax.barh(y_pos, top['mean_abs_shap'].values, color=viridis_colors,
                   alpha=0.9, edgecolor='white', linewidth=0.5, height=0.65)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(top['feature'].values, fontsize=9)
    ax.set_xlabel('Mean |SHAP value|', fontsize=11)
    ax.set_title('(b) Feature Importance', fontsize=12, fontweight='bold', loc='left')
    ax.set_xlim(0, top['mean_abs_shap'].max() * 1.15)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.8)

    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, top['mean_abs_shap'].values)):
        ax.text(val + 0.01, i, f'{val:.3f}', va='center', ha='left', fontsize=8, color='#333333')


def plot_group_diff_bar(ax, shap_importance_df, cell_props, sample_labels, max_display=10):
    """Horizontal bar: Responder vs Non-responder proportion difference.
    Viridis-based: negative = blue-purple, positive = yellow-green."""
    top = shap_importance_df.head(max_display).copy()
    features = top['feature'].values

    resp_means = [cell_props.loc[sample_labels == 1, f].mean() for f in features]
    nonresp_means = [cell_props.loc[sample_labels == 0, f].mean() for f in features]
    diffs = np.array(resp_means) - np.array(nonresp_means)
    max_abs = np.max(np.abs(diffs)) if len(diffs) > 0 else 1

    # Viridis-based diverging: negative = purple end, positive = yellow-green end
    colors = []
    for d in diffs:
        if d > 0:
            intensity = 0.5 + 0.4 * (d / max_abs)  # 0.5-0.9 range
            colors.append(plt.cm.viridis(intensity))
        else:
            intensity = 0.5 - 0.4 * (abs(d) / max_abs)  # 0.1-0.5 range
            colors.append(plt.cm.viridis(intensity))

    y_pos = np.arange(len(features))
    bars = ax.barh(y_pos, diffs, color=colors, alpha=0.9, edgecolor='white', linewidth=0.5, height=0.6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(features, fontsize=9)
    ax.axvline(x=0, color='#333333', linewidth=0.8)
    ax.set_xlabel('Proportion difference\n(Responder − Non-responder)', fontsize=10)
    ax.set_title('(c) Cell Proportion Difference', fontsize=12, fontweight='bold', loc='left')
    ax.invert_yaxis()

    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, diffs)):
        if val >= 0:
            ax.text(val + 0.001, i, f'{val:.4f}', va='center', ha='left', fontsize=8, color='#333333')
        else:
            # For very long negative bars (Treg_CCR8), place label on the right side
            if abs(val) > 0.015:
                ax.text(0.001, i, f'{val:.4f}', va='center', ha='left', fontsize=8, color='#333333')
            else:
                ax.text(val - 0.001, i, f'{val:.4f}', va='center', ha='right', fontsize=8, color='#333333')


def plot_group_legend(ax):
    """Group legend."""
    ax.axis('off')
    handles = []
    for group in ['Treg', 'Macrophage', 'CD8 T', 'CD4 T', 'B cell', 'NK', 'DC', 'Other']:
        color = GROUP_BASE_COLOR.get(group, '#757575')
        handles.append(mpatches.Patch(color=color, label=group))
    ax.legend(handles=handles, loc='center', frameon=False, fontsize=9,
              title='Cell Type Groups', title_fontsize=10, ncol=4,
              columnspacing=1.5, handletextpad=0.5)


def main():
    print("=" * 60)
    print("Publication SHAP Figure (XGBoost-SHAP Style)")
    print("=" * 60)

    # 1. Load/compute SHAP
    shap_matrix, X_all_test, sample_labels, feature_names, cell_props = compute_or_load_shap()

    # 2. Importance ranking
    mean_abs = np.abs(shap_matrix).mean(axis=0)
    shap_importance = pd.DataFrame({
        'feature': feature_names,
        'mean_abs_shap': mean_abs
    }).sort_values('mean_abs_shap', ascending=False)

    print("\nTop 10 features:")
    print(shap_importance.head(10).to_string(index=False))

    # 3. Figure layout
    fig = plt.figure(figsize=(20, 14), dpi=300)
    gs = GridSpec(3, 3, figure=fig,
                  width_ratios=[1.2, 0.85, 0.95],
                  height_ratios=[1, 1, 0.12],
                  wspace=0.35, hspace=0.35)

    # (a) Native shap.summary_plot (classic red-blue coolwarm)
    ax_beeswarm = fig.add_subplot(gs[0:2, 0])
    plt.sca(ax_beeswarm)
    shap.summary_plot(
        shap_matrix, X_all_test,
        feature_names=feature_names,
        max_display=20,
        show=False,
        plot_size=None
    )
    ax_beeswarm.set_title('(a) SHAP Beeswarm Plot', fontsize=12, fontweight='bold', loc='left')
    
    # Adjust colorbar
    cbar_ax = fig.axes[-1]
    cbar_ax.set_ylabel('')
    cbar_ax.set_xlabel('Feature value (z-score)', fontsize=9)
    cbar_ax.xaxis.set_label_position('top')
    cbar_ax.tick_params(labelsize=8)

    # (b) Mean |SHAP| bar (top-middle)
    ax_bar = fig.add_subplot(gs[0, 1])
    plot_mean_shap_bar(ax_bar, shap_importance, max_display=20)

    # (c) Group diff bar (bottom-middle)
    ax_diff = fig.add_subplot(gs[1, 1])
    plot_group_diff_bar(ax_diff, shap_importance, cell_props, sample_labels, max_display=10)

    # (d) Rose chart: larger, no pct labels, viridis high=yellow-green
    gs_rose = gs[0:2, 2].subgridspec(2, 1, height_ratios=[1.75, 0.85], hspace=0.02)
    
    ax_rose = fig.add_subplot(gs_rose[0], projection='polar')
    
    top_n = 10
    top = shap_importance.head(top_n).copy()
    values = top['mean_abs_shap'].values
    total = values.sum()
    pct = values / total * 100
    
    # Viridis: high contribution = yellow-green, low = blue-purple
    viridis = plt.cm.viridis(np.linspace(0.85, 0.15, top_n))
    
    N = len(values)
    theta = np.linspace(0, 2 * np.pi, N, endpoint=False)
    width = 2 * np.pi / N * 0.88
    
    ax_rose.set_theta_zero_location('N')
    ax_rose.set_theta_direction(-1)
    
    bars = ax_rose.bar(theta, values, width=width, bottom=0.01,
                       color=viridis, alpha=0.92,
                       edgecolor='white', linewidth=1.5)
    
    # No percentage labels on rose chart
    ax_rose.set_ylim(0, values.max() * 1.08)
    ax_rose.set_yticks([])
    ax_rose.set_xticks([])
    ax_rose.set_title('(d) Feature Contribution %', fontsize=11, fontweight='bold', pad=12)
    ax_rose.spines['polar'].set_visible(False)
    ax_rose.grid(False)
    
    # Legend below: horizontal bars with percentage
    ax_legend = fig.add_subplot(gs_rose[1])
    ax_legend.axis('off')
    
    for i, (_, row) in enumerate(top.iterrows()):
        feat = row['feature'].replace('CD4T_', '').replace('CD8T_', '').replace('Mφ_', 'Mφ')
        color = viridis[i]
        p = pct[i]
        
        if i < 5:
            x = 0.06
            y = 0.95 - i * 0.19
        else:
            x = 0.54
            y = 0.95 - (i - 5) * 0.19
        
        ax_legend.add_patch(plt.Rectangle((x, y), 0.10, 0.08,
                                           transform=ax_legend.transAxes,
                                           facecolor=color, edgecolor='white',
                                           linewidth=0.8, clip_on=False))
        ax_legend.text(x + 0.13, y + 0.04, f'{feat}\n{p:.1f}%',
                       transform=ax_legend.transAxes, fontsize=8.5,
                       va='center', ha='left', linespacing=1.1)
    
    ax_legend.set_xlim(0, 1)
    ax_legend.set_ylim(0, 1)

    # Bottom: model info only (remove Cell Type Groups legend)
    ax_info = fig.add_subplot(gs[2, 0:3])
    ax_info.axis('off')
    n_resp = (sample_labels == 1).sum()
    n_non = (sample_labels == 0).sum()
    cv_df = pd.read_csv(os.path.join(DATA_DIR, 'cv_results.csv'))
    auc_mean = cv_df['auc'].mean()
    info_text = (
        f"XGBoost 5-fold CV | Mean AUC = {auc_mean:.3f} | "
        f"n = {len(sample_labels)} (Responders = {n_resp}, Non-responders = {n_non}) | "
        f"Features = {len(feature_names)} cell subtypes"
    )
    ax_info.text(0.5, 0.5, info_text, transform=ax_info.transAxes, ha='center', va='center',
                 fontsize=10, style='italic', color='#555555',
                 bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8f9fa', edgecolor='#dee2e6', linewidth=1))

    plt.suptitle(
        'XGBoost-SHAP Analysis: Cell Subtype Proportions Predict Anti-PD-1 Response in NSCLC',
        fontsize=15, fontweight='bold', y=0.98
    )

    # Save
    for fmt, dpi in [('png', 600), ('pdf', 600), ('tiff', 300)]:
        path = os.path.join(FIG_DIR, f'pub_fig_shap_composite.{fmt}')
        kwargs = {'bbox_inches': 'tight', 'facecolor': 'white'}
        if fmt == 'tiff':
            kwargs['pil_kwargs'] = {'compression': 'tiff_lzw'}
        plt.savefig(path, dpi=dpi, **kwargs)
        print(f"[Saved] {path}")

    plt.close()
    print("=" * 60)


if __name__ == '__main__':
    main()
