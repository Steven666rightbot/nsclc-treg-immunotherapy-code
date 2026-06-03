"""
Figure 3A + 3D: DoRothEA TF activity analysis
- Panel A: Top TF bar plots (CCR8+ and MKI67+ R vs NR)
- Panel D: Key TF violin plots (HSF1, BCL6, RFX5)

Publication-quality, color-blind friendly.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import warnings
import os

warnings.filterwarnings('ignore')

# ============================================================
# Paths
# ============================================================
DATA_DIR = 'D:/research/tomato'
OUT_DIR = 'D:/research/cucumber/fig3'
TF_ACT_PATH = os.path.join(DATA_DIR, 'results/scenic_input/dorothea_tf_activity.csv')
META_PATH = os.path.join(DATA_DIR, 'results/scenic_input/treg_no_foxp3_metadata.csv')
OUT_FIG_DIR = OUT_DIR
OUT_RESULTS_DIR = OUT_DIR
os.makedirs(OUT_FIG_DIR, exist_ok=True)
os.makedirs(OUT_RESULTS_DIR, exist_ok=True)

# ============================================================
# Publication style settings
# ============================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.7,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.width': 0.7,
    'ytick.major.width': 0.7,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
})

# Color-blind friendly palette
# Responder = green, Non-responder = purple
COLOR_RESPONDER = '#8bc98b'
COLOR_NONRESPONDER = '#b08ebf'
COLOR_HIGHER_RESPONDER = '#3a86ff'   # Blue bars for higher in Responder
COLOR_HIGHER_NONRESPONDER = '#d65f5f'  # Red bars for higher in Non-responder


def load_and_merge():
    """Load TF activity matrix and metadata, merge on CellID."""
    print("Loading TF activity matrix...")
    tf_df = pd.read_csv(TF_ACT_PATH, index_col=0)
    print(f"  TF matrix: {tf_df.shape[0]} cells x {tf_df.shape[1]} TFs")

    print("Loading metadata...")
    meta = pd.read_csv(META_PATH)
    print(f"  Metadata: {meta.shape[0]} cells")

    # Merge
    common_cells = meta['CellID'].isin(tf_df.index)
    print(f"  Cells matched: {common_cells.sum()} / {meta.shape[0]}")
    meta = meta[common_cells].copy()

    data = tf_df.loc[meta['CellID']].copy()
    data['sub_cell_type'] = meta['sub_cell_type'].values
    data['response'] = meta['response'].values

    return data, meta


def differential_analysis(data, subtype):
    """
    Mann-Whitney U test per TF for Responder vs Non-responder within a subtype.
    Returns DataFrame with effect_size, p_value, R_median, NR_median.
    """
    sub = data[data['sub_cell_type'] == subtype].copy()
    tf_names = [c for c in sub.columns if c not in ['sub_cell_type', 'response']]

    results = []
    for tf in tf_names:
        r_vals = sub.loc[sub['response'] == 'Responder', tf].dropna()
        nr_vals = sub.loc[sub['response'] == 'Non-responder', tf].dropna()

        if len(r_vals) < 3 or len(nr_vals) < 3:
            continue

        stat, pval = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
        r_med = r_vals.median()
        nr_med = nr_vals.median()
        effect = r_med - nr_med

        results.append({
            'TF': tf,
            'effect_size': effect,
            'p_value': pval,
            'R_median': r_med,
            'NR_median': nr_med,
        })

    res_df = pd.DataFrame(results)

    # FDR correction
    reject, padj, _, _ = multipletests(res_df['p_value'], method='fdr_bh')
    res_df['padj'] = padj
    res_df['significant'] = reject

    # Sort by absolute effect size
    res_df = res_df.sort_values('effect_size', ascending=False).reset_index(drop=True)
    res_df['abs_effect'] = res_df['effect_size'].abs()

    return res_df


def make_fig3a(all_results, subtype='CCR8', suffix=''):
    """Generate top-TF barplot panel."""
    res = all_results[subtype].copy()
    res['direction'] = np.where(res['effect_size'] > 0, 'Responder', 'Non-responder')

    # Top 10 by absolute effect size
    top = res.nlargest(10, 'abs_effect').iloc[::-1].copy()

    # Significance stars
    def star(p):
        if p < 0.001:
            return '***'
        elif p < 0.01:
            return '**'
        elif p < 0.05:
            return '*'
        return 'ns'

    top['star'] = top['padj'].apply(star)
    top['label'] = top['TF'] + '  ' + top['star']

    colors = [COLOR_HIGHER_RESPONDER if e > 0 else COLOR_HIGHER_NONRESPONDER
              for e in top['effect_size']]

    fig, ax = plt.subplots(figsize=(5.5, 4.0))

    bars = ax.barh(range(len(top)), top['effect_size'].values, color=colors,
                   edgecolor='white', linewidth=0.3, height=0.7)

    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top['label'].values, fontsize=8)
    ax.set_xlabel('Effect size (Responder − Non-responder)', fontsize=9)

    # Vertical line at 0
    ax.axvline(0, color='grey', linewidth=0.6, linestyle='-', alpha=0.5)

    # Annotate star significance at bar ends
    for i, (_, row) in enumerate(top.iterrows()):
        val = row['effect_size']
        star_txt = row['star']
        offset = 0.02 * (top['effect_size'].max() - top['effect_size'].min())
        if star_txt != 'ns':
            if val > 0:
                ax.text(val + offset, i, star_txt, va='center', fontsize=7,
                        color=COLOR_HIGHER_RESPONDER, fontweight='bold')
            else:
                ax.text(val - offset, i, star_txt, va='center', fontsize=7,
                        color=COLOR_HIGHER_NONRESPONDER, fontweight='bold', ha='right')

    # Panel label
    ax.text(-0.15, 1.03, '(A)', transform=ax.transAxes, fontsize=12,
            fontweight='bold', va='bottom', ha='left')

    ax.set_title(subtype, fontsize=10, fontweight='bold')

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_HIGHER_RESPONDER, edgecolor='white',
                       label='Higher in Responder'),
        mpatches.Patch(facecolor=COLOR_HIGHER_NONRESPONDER, edgecolor='white',
                       label='Higher in Non-responder'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=7,
              frameon=True, framealpha=0.9, edgecolor='grey')

    ax.set_xlim(left=top['effect_size'].min() - 0.15 * abs(top['effect_size'].min()),
                right=top['effect_size'].max() + 0.15 * abs(top['effect_size'].max()))

    # Remove left spine (keep right for tick labels)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    plt.tight_layout()
    return fig, ax, top


def make_fig3d(data):
    """Generate violin plots for key TFs: HSF1, BCL6, RFX5."""
    key_tfs = ['HSF1', 'BCL6', 'RFX5']
    subtypes = ['CD4T_Treg_CCR8', 'CD4T_Treg_MKI67']
    subtype_labels = ['CCR8+', 'MKI67+']
    response_groups = ['Responder', 'Non-responder']
    group_colors = [COLOR_RESPONDER, COLOR_NONRESPONDER]

    fig, axes = plt.subplots(1, 3, figsize=(8.0, 3.2), sharey=True)

    for idx, tf in enumerate(key_tfs):
        ax = axes[idx]

        plot_data = []
        positions = []
        pos = 0
        tick_labels = []
        tick_positions = []
        x_positions = []

        for si, st in enumerate(subtypes):
            for ri, resp in enumerate(response_groups):
                vals = data.loc[(data['sub_cell_type'] == st) &
                                (data['response'] == resp), tf].dropna().values
                plot_data.append(vals)
                positions.append(pos)
                x_positions.append(pos)
                tick_labels.append(f'{subtype_labels[si]}\n{resp[0]}')
                tick_positions.append(pos)
                pos += 1
            # Gap between subtypes
            pos += 0.5

        # Violin plots
        parts = ax.violinplot(plot_data, positions=positions, showmeans=False,
                              showmedians=False, showextrema=False,
                              widths=0.65)

        for pi, pc in enumerate(parts['bodies']):
            group_idx = pi % 2
            pc.set_facecolor(group_colors[group_idx])
            pc.set_alpha(0.6)
            pc.set_edgecolor('none')

        # Jittered points (strip plot)
        for pi, vals in enumerate(plot_data):
            group_idx = pi % 2
            jitter = np.random.normal(0, 0.04, size=len(vals))
            ax.scatter(positions[pi] + jitter, vals, s=1.5,
                       color=group_colors[group_idx], alpha=0.3, edgecolors='none',
                       rasterized=True)

        # Medians as horizontal lines
        for pi, vals in enumerate(plot_data):
            med = np.median(vals)
            ax.plot([positions[pi] - 0.15, positions[pi] + 0.15],
                    [med, med], color='black', linewidth=1.2)

        # Statistical annotations (Mann-Whitney U, within each subtype)
        for si, st in enumerate(subtypes):
            r_vals = data.loc[(data['sub_cell_type'] == st) &
                              (data['response'] == 'Responder'), tf].dropna().values
            nr_vals = data.loc[(data['sub_cell_type'] == st) &
                               (data['response'] == 'Non-responder'), tf].dropna().values

            if len(r_vals) >= 3 and len(nr_vals) >= 3:
                _, pval = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
                if pval < 0.001:
                    star_txt = '***'
                elif pval < 0.01:
                    star_txt = '**'
                elif pval < 0.05:
                    star_txt = '*'
                else:
                    star_txt = 'ns'

                # Position the star between the two groups
                left_pos = positions[si * 2]
                right_pos = positions[si * 2 + 1]
                mid_pos = (left_pos + right_pos) / 2

                # Find y-position (slightly above max)
                all_vals = np.concatenate([r_vals, nr_vals])
                y_max = np.percentile(all_vals, 95)
                y_range = np.percentile(all_vals, 95) - np.percentile(all_vals, 5)
                y_pos = y_max + 0.1 * y_range

                # Only add annotation if it doesn't overlap
                ax.text(mid_pos, y_pos, star_txt, ha='center', va='bottom',
                        fontsize=8, fontweight='bold')

                # Line
                ax.plot([left_pos, right_pos],
                        [y_pos - 0.02 * y_range, y_pos - 0.02 * y_range],
                        color='black', linewidth=0.6)

        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=7, ha='center')
        ax.set_title(tf, fontsize=10, fontweight='bold')

        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

        if idx == 0:
            ax.set_ylabel('DoRothEA activity score', fontsize=9)
        else:
            ax.set_ylabel('')

    # Rotate x-axis tick labels to prevent overlap
    for ax in axes:
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

    # Panel label
    axes[0].text(-0.12, 1.05, '(D)', transform=axes[0].transAxes, fontsize=12,
                 fontweight='bold', va='bottom', ha='left')

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_RESPONDER, alpha=0.6,
                       label='Responder'),
        mpatches.Patch(facecolor=COLOR_NONRESPONDER, alpha=0.6,
                       label='Non-responder'),
    ]
    fig.legend(handles=legend_elements, loc='upper right',
               fontsize=7.5, frameon=True, framealpha=0.9,
               edgecolor='grey', ncol=1, bbox_to_anchor=(0.98, 0.98))

    plt.subplots_adjust(bottom=0.15)
    plt.tight_layout()
    return fig


def main():
    print("=" * 60)
    print("Figure 3A + 3D: DoRothEA TF Activity Analysis")
    print("=" * 60)

    # 1. Load and merge data
    data, meta = load_and_merge()

    # 2. Differential analysis per subtype
    print("\nRunning differential analysis...")
    subtypes_map = {
        'CCR8': 'CD4T_Treg_CCR8',
        'MKI67': 'CD4T_Treg_MKI67',
    }

    all_results = {}
    for label, subtype in subtypes_map.items():
        print(f"  {label} ({subtype})...")
        res = differential_analysis(data, subtype)
        all_results[label] = res
        print(f"    {len(res)} TFs tested, {res['significant'].sum()} significant (FDR < 0.05)")

        # Save
        out_path = os.path.join(OUT_RESULTS_DIR, f'fig3a_{label}_tf_differential_results.csv')
        res.to_csv(out_path, index=False)
        print(f"    Saved to {out_path}")

    # 3. Save combined results
    ccr8_res = all_results['CCR8'].copy()
    mki67_res = all_results['MKI67'].copy()
    ccr8_res['subtype'] = 'CCR8'
    mki67_res['subtype'] = 'MKI67'
    combined = pd.concat([ccr8_res, mki67_res], ignore_index=True)
    combined_path = os.path.join(OUT_RESULTS_DIR, 'fig3a_tf_differential_results.csv')
    combined.to_csv(combined_path, index=False)
    print(f"\nCombined results saved to {combined_path}")

    # 4. Figure 3A — Top TF bar plots
    print("\nGenerating Figure 3A...")
    fig_a1, ax1, top1 = make_fig3a(all_results, 'CCR8')
    fig_a1.savefig(os.path.join(OUT_FIG_DIR, 'fig3a_top_tf_bar_CCR8.png'), dpi=300)
    fig_a1.savefig(os.path.join(OUT_FIG_DIR, 'fig3a_top_tf_bar_CCR8.pdf'), dpi=300)
    plt.close(fig_a1)

    fig_a2, ax2, top2 = make_fig3a(all_results, 'MKI67')
    fig_a2.savefig(os.path.join(OUT_FIG_DIR, 'fig3a_top_tf_bar_MKI67.png'), dpi=300)
    fig_a2.savefig(os.path.join(OUT_FIG_DIR, 'fig3a_top_tf_bar_MKI67.pdf'), dpi=300)
    plt.close(fig_a2)

    print("  Saved: fig3a_top_tf_bar_CCR8.png/pdf")
    print("  Saved: fig3a_top_tf_bar_MKI67.png/pdf")

    # Verify key TFs
    for label in ['CCR8', 'MKI67']:
        res = all_results[label]
        for key_tf in ['HSF1', 'MAX', 'BCL6', 'RFX5', 'EGR1']:
            match = res[res['TF'] == key_tf]
            if len(match) > 0:
                row = match.iloc[0]
                print(f"  {label} {key_tf}: effect={row['effect_size']:.4f}, "
                      f"p={row['p_value']:.2e}, padj={row['padj']:.2e}, "
                      f"R_med={row['R_median']:.3f}, NR_med={row['NR_median']:.3f}")
            else:
                print(f"  {label} {key_tf}: not found")

    # 5. Figure 3D — Key TF violin plots
    print("\nGenerating Figure 3D...")
    fig_d = make_fig3d(data)
    fig_d.savefig(os.path.join(OUT_FIG_DIR, 'fig3d_tf_violin.png'), dpi=300)
    fig_d.savefig(os.path.join(OUT_FIG_DIR, 'fig3d_tf_violin.pdf'), dpi=300)
    plt.close(fig_d)
    print("  Saved: fig3d_tf_violin.png/pdf")

    # 6. Combined figure (A + D for reference)
    print("\nGenerating combined Figure 3A+D...")
    fig_combined = plt.figure(figsize=(12.5, 7.5))
    gs = GridSpec(2, 2, figure=fig_combined, width_ratios=[1, 1], height_ratios=[1.3, 1],
                  hspace=0.3, wspace=0.35)

    # Rebuild panel A (CCR8 left, MKI67 right)
    def panel_a_sub(ax, results, subtype_label):
        res = results[subtype_label].copy()
        res['direction'] = np.where(res['effect_size'] > 0, 'Responder', 'Non-responder')
        top = res.nlargest(10, 'abs_effect').iloc[::-1].copy()

        def star(p):
            if p < 0.001:
                return '***'
            elif p < 0.01:
                return '**'
            elif p < 0.05:
                return '*'
            return 'ns'

        top['star'] = top['padj'].apply(star)
        top['label'] = top['TF'] + '  ' + top['star']

        colors = [COLOR_HIGHER_RESPONDER if e > 0 else COLOR_HIGHER_NONRESPONDER
                  for e in top['effect_size']]

        ax.barh(range(len(top)), top['effect_size'].values, color=colors,
                edgecolor='white', linewidth=0.3, height=0.7)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top['label'].values, fontsize=8)
        ax.axvline(0, color='grey', linewidth=0.6, linestyle='-', alpha=0.5)

        for i, (_, row) in enumerate(top.iterrows()):
            val = row['effect_size']
            star_txt = row['star']
            xrange = top['effect_size'].max() - top['effect_size'].min()
            offset = 0.02 * xrange
            if star_txt != 'ns':
                if val > 0:
                    ax.text(val + offset, i, star_txt, va='center', fontsize=7,
                            color=COLOR_HIGHER_RESPONDER, fontweight='bold')
                else:
                    ax.text(val - offset, i, star_txt, va='center', fontsize=7,
                            color=COLOR_HIGHER_NONRESPONDER, fontweight='bold', ha='right')

        ax.set_title(subtype_label, fontsize=10, fontweight='bold')

        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

        return top

    # Panel A: CCR8 (left)
    ax_a1 = fig_combined.add_subplot(gs[0, 0])
    top_a = panel_a_sub(ax_a1, all_results, 'CCR8')
    ax_a1.set_xlabel('Effect size (Responder − Non-responder)', fontsize=8.5)
    ax_a1.text(-0.12, 1.03, '(A)', transform=ax_a1.transAxes, fontsize=12,
               fontweight='bold', va='bottom', ha='left')
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_HIGHER_RESPONDER, edgecolor='white',
                       label='Higher in Responder'),
        mpatches.Patch(facecolor=COLOR_HIGHER_NONRESPONDER, edgecolor='white',
                       label='Higher in Non-responder'),
    ]
    ax_a1.legend(handles=legend_elements, loc='lower right', fontsize=6.5,
                 frameon=True, framealpha=0.9, edgecolor='grey')

    ax_a1.set_xlim(left=top_a['effect_size'].min() - 0.15 * abs(top_a['effect_size'].min()),
                   right=top_a['effect_size'].max() + 0.15 * abs(top_a['effect_size'].max()))

    # Panel A: MKI67 (right)
    ax_a2 = fig_combined.add_subplot(gs[0, 1])
    top_b = panel_a_sub(ax_a2, all_results, 'MKI67')
    ax_a2.set_xlabel('Effect size (Responder − Non-responder)', fontsize=8.5)
    ax_a2.set_ylabel('')  # No ylabel on right panel
    ax_a2.set_yticklabels([])
    ax_a2.text(-0.12, 1.03, '(B)', transform=ax_a2.transAxes, fontsize=12,
               fontweight='bold', va='bottom', ha='left')

    ax_a2.set_xlim(left=top_b['effect_size'].min() - 0.15 * abs(top_b['effect_size'].min()),
                   right=top_b['effect_size'].max() + 0.15 * abs(top_b['effect_size'].max()))

    # Panel D: Violin plots (spanning bottom row across both columns)
    # Actually, per spec it's side by side in one row, so let's use the full bottom
    gs_d = GridSpec(1, 3, figure=fig_combined, bottom=0.05, top=0.44,
                    left=0.08, right=0.98, wspace=0.3)

    key_tfs = ['HSF1', 'BCL6', 'RFX5']
    subtypes = ['CD4T_Treg_CCR8', 'CD4T_Treg_MKI67']
    subtype_labels = ['CCR8+', 'MKI67+']
    response_groups = ['Responder', 'Non-responder']
    group_colors = [COLOR_RESPONDER, COLOR_NONRESPONDER]

    for idx, tf in enumerate(key_tfs):
        ax = fig_combined.add_subplot(gs_d[0, idx])

        plot_data = []
        positions = []
        tick_labels = []
        tick_positions = []
        pos = 0

        for si, st in enumerate(subtypes):
            for ri, resp in enumerate(response_groups):
                vals = data.loc[(data['sub_cell_type'] == st) &
                                (data['response'] == resp), tf].dropna().values
                plot_data.append(vals)
                positions.append(pos)
                tick_labels.append(f'{subtype_labels[si]}\n{resp[0]}')
                tick_positions.append(pos)
                pos += 1
            pos += 0.5

        # Violins
        parts = ax.violinplot(plot_data, positions=positions, showmeans=False,
                              showmedians=False, showextrema=False, widths=0.65)

        for pi, pc in enumerate(parts['bodies']):
            group_idx = pi % 2
            pc.set_facecolor(group_colors[group_idx])
            pc.set_alpha(0.6)
            pc.set_edgecolor('none')

        # Jitter
        for pi, vals in enumerate(plot_data):
            group_idx = pi % 2
            jitter = np.random.normal(0, 0.04, size=len(vals))
            ax.scatter(positions[pi] + jitter, vals, s=1.5,
                       color=group_colors[group_idx], alpha=0.3, edgecolors='none',
                       rasterized=True)

        # Medians
        for pi, vals in enumerate(plot_data):
            med = np.median(vals)
            ax.plot([positions[pi] - 0.15, positions[pi] + 0.15],
                    [med, med], color='black', linewidth=1.2)

        # Stats annotations
        for si, st in enumerate(subtypes):
            r_vals = data.loc[(data['sub_cell_type'] == st) &
                              (data['response'] == 'Responder'), tf].dropna().values
            nr_vals = data.loc[(data['sub_cell_type'] == st) &
                               (data['response'] == 'Non-responder'), tf].dropna().values

            if len(r_vals) >= 3 and len(nr_vals) >= 3:
                _, pval = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
                if pval < 0.001:
                    star_txt = '***'
                elif pval < 0.01:
                    star_txt = '**'
                elif pval < 0.05:
                    star_txt = '*'
                else:
                    star_txt = 'ns'

                left_pos = positions[si * 2]
                right_pos = positions[si * 2 + 1]
                mid_pos = (left_pos + right_pos) / 2

                all_vals = np.concatenate([r_vals, nr_vals])
                y_max = np.percentile(all_vals, 95)
                y_range = np.percentile(all_vals, 95) - np.percentile(all_vals, 5)
                y_pos = y_max + 0.12 * y_range

                ax.text(mid_pos, y_pos, star_txt, ha='center', va='bottom',
                        fontsize=8, fontweight='bold')
                ax.plot([left_pos, right_pos],
                        [y_pos - 0.02 * y_range, y_pos - 0.02 * y_range],
                        color='black', linewidth=0.6)

        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=6.5, ha='center')
        ax.set_title(tf, fontsize=10, fontweight='bold')

        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

        if idx == 0:
            ax.set_ylabel('DoRothEA activity score', fontsize=8.5)
        else:
            ax.set_ylabel('')

    # Panel D label
    ax_d0 = fig_combined.get_axes()[-3]  # First violin subplot
    ax_d0.text(-0.12, 1.08, '(D)', transform=ax_d0.transAxes, fontsize=12,
               fontweight='bold', va='bottom', ha='left')

    # Legend for violin
    legend_elements_d = [
        mpatches.Patch(facecolor=COLOR_RESPONDER, alpha=0.6, label='Responder'),
        mpatches.Patch(facecolor=COLOR_NONRESPONDER, alpha=0.6, label='Non-responder'),
    ]
    fig_combined.legend(handles=legend_elements_d, loc='upper right',
                        fontsize=7.5, frameon=True, framealpha=0.9,
                        edgecolor='grey', ncol=1, bbox_to_anchor=(0.98, 0.45))

    combined_path_fig = os.path.join(OUT_FIG_DIR, 'fig3a_top_tf_bar.png')
    fig_combined.savefig(combined_path_fig, dpi=300)
    fig_combined.savefig(os.path.join(OUT_FIG_DIR, 'fig3a_top_tf_bar.pdf'), dpi=300)
    plt.close(fig_combined)
    print(f"  Saved combined: fig3a_top_tf_bar.png/pdf")

    # For the spec output, also save as fig3a directly
    # The spec says save as figures/fig3a_top_tf_bar.png and .pdf

    print("\n" + "=" * 60)
    print("Done! All figures saved.")
    print("=" * 60)


if __name__ == '__main__':
    main()
