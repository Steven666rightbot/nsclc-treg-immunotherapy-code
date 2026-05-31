"""
Publication-quality spatial distribution figure
Shows CCR8+ Treg enrichment in hard_stroma vs random MKI67+ distribution
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from scipy import stats

BASE_DIR = r'D:\Research\tomato'
FIG_DIR = os.path.join(BASE_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# Publication style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
})

# Colors
REGION_COLORS = {
    'hard_stroma': '#E74C3C',
    'tumor_core': '#3498DB',
    'other': '#95A5A6',
}

# ============================================================================
# Load data
# ============================================================================
print("Loading spatial metadata...")
pdac = pd.read_csv(os.path.join(BASE_DIR, 'results', 'spatial', 'spatial_spot_metadata.csv'))
nsclc = pd.read_csv(os.path.join(BASE_DIR, 'results', 'spatial_nsclc', 'nsclc_spot_metadata.csv'))

# Select representative samples
pdac_sample = 'PDAC_1'
nsclc_sample = 'GSM5702474_TD2'

pdac_sub = pdac[pdac['sample_id'] == pdac_sample].copy()
nsclc_sub = nsclc[nsclc['sample_id'] == nsclc_sample].copy()

print(f"PDAC {pdac_sample}: {len(pdac_sub)} spots")
print(f"NSCLC {nsclc_sample}: {len(nsclc_sub)} spots")

# ============================================================================
# Helper: compute region convex hull for outline
# ============================================================================
def region_outline(df, region_name, alpha=0.15):
    """Return approximate region boundary using convex hull of spot centers."""
    sub = df[df['region'] == region_name].copy()
    if len(sub) < 3:
        return None
    pts = sub[['array_col', 'array_row']].values.astype(float)
    from scipy.spatial import ConvexHull
    try:
        hull = ConvexHull(pts)
        return pts[hull.vertices]
    except Exception:
        return None

# ============================================================================
# Helper: spatial scatter plot
# ============================================================================
def plot_spatial_sample(ax, df, score_col, title, region_outline_dict=None,
                        vmin=-2, vmax=2, cmap='RdYlBu_r'):
    """Plot Visium spots colored by score."""
    # Background: light gray for all spots
    ax.scatter(df['array_col'], df['array_row'], c='#E8E8E8', s=12, zorder=1, marker='h')
    
    # Colored by score
    sc = ax.scatter(df['array_col'], df['array_row'], c=df[score_col],
                    cmap=cmap, vmin=vmin, vmax=vmax, s=18, zorder=2,
                    edgecolors='none', marker='h')
    
    # Region outlines
    if region_outline_dict:
        for region, verts in region_outline_dict.items():
            if verts is not None and len(verts) > 2:
                color = REGION_COLORS.get(region, '#333333')
                poly = Polygon(verts, closed=True, fill=False,
                               edgecolor=color, linewidth=1.5, linestyle='--', zorder=3)
                ax.add_patch(poly)
    
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_title(title, fontweight='bold', fontsize=10)
    ax.set_xlabel('Array column')
    ax.set_ylabel('Array row')
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Remove spines
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    
    return sc

# ============================================================================
# Compute region outlines
# ============================================================================
print("Computing region outlines...")
pdac_outlines = {
    'hard_stroma': region_outline(pdac_sub, 'hard_stroma'),
    'tumor_core': region_outline(pdac_sub, 'tumor_core'),
    'other': region_outline(pdac_sub, 'other'),
}
nsclc_outlines = {
    'hard_stroma': region_outline(nsclc_sub, 'hard_stroma'),
    'tumor_core': region_outline(nsclc_sub, 'tumor_core'),
    'other': region_outline(nsclc_sub, 'other'),
}

# ============================================================================
# Statistical tests for annotation
# ============================================================================
def region_comparison(df, score_col):
    """Compare hard_stroma vs tumor_core."""
    hs = df[df['region'] == 'hard_stroma'][score_col].dropna()
    tc = df[df['region'] == 'tumor_core'][score_col].dropna()
    if len(hs) < 5 or len(tc) < 5:
        return None, None, None, None
    stat, p = stats.mannwhitneyu(hs, tc, alternative='two-sided')
    return hs.median(), tc.median(), stat, p

pdac_ccr8_hs, pdac_ccr8_tc, _, pdac_ccr8_p = region_comparison(pdac_sub, 'ccr8_score_z')
pdac_mki67_hs, pdac_mki67_tc, _, pdac_mki67_p = region_comparison(pdac_sub, 'mki67_score_z')
nsclc_ccr8_hs, nsclc_ccr8_tc, _, nsclc_ccr8_p = region_comparison(nsclc_sub, 'ccr8_score_z')
nsclc_mki67_hs, nsclc_mki67_tc, _, nsclc_mki67_p = region_comparison(nsclc_sub, 'mki67_score_z')

print(f"\nPDAC CCR8+  HS median={pdac_ccr8_hs:.3f}  TC median={pdac_ccr8_tc:.3f}  p={pdac_ccr8_p:.2e}")
print(f"PDAC MKI67+ HS median={pdac_mki67_hs:.3f}  TC median={pdac_mki67_tc:.3f}  p={pdac_mki67_p:.2e}")
print(f"NSCLC CCR8+ HS median={nsclc_ccr8_hs:.3f}  TC median={nsclc_ccr8_tc:.3f}  p={nsclc_ccr8_p:.2e}")
print(f"NSCLC MKI67+ HS median={nsclc_mki67_hs:.3f}  TC median={nsclc_mki67_tc:.3f}  p={nsclc_mki67_p:.2e}")

# ============================================================================
# Build figure
# ============================================================================
print("\nGenerating figure...")
fig = plt.figure(figsize=(14, 10))

# Use GridSpec for better control
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.08], width_ratios=[1, 1, 0.06],
                      hspace=0.25, wspace=0.15)

# Row 0: PDAC
ax_a = fig.add_subplot(gs[0, 0])
sc_a = plot_spatial_sample(ax_a, pdac_sub, 'ccr8_score_z',
    f'A. PDAC ({pdac_sample})\nCCR8+ Treg Signature',
    region_outline_dict=pdac_outlines)

ax_b = fig.add_subplot(gs[0, 1])
sc_b = plot_spatial_sample(ax_b, pdac_sub, 'mki67_score_z',
    f'B. PDAC ({pdac_sample})\nMKI67+ Treg Signature',
    region_outline_dict=pdac_outlines)

# Row 1: NSCLC
ax_c = fig.add_subplot(gs[1, 0])
sc_c = plot_spatial_sample(ax_c, nsclc_sub, 'ccr8_score_z',
    f'C. NSCLC ({nsclc_sample})\nCCR8+ Treg Signature',
    region_outline_dict=nsclc_outlines)

ax_d = fig.add_subplot(gs[1, 1])
sc_d = plot_spatial_sample(ax_d, nsclc_sub, 'mki67_score_z',
    f'D. NSCLC ({nsclc_sample})\nMKI67+ Treg Signature',
    region_outline_dict=nsclc_outlines)

# Shared colorbar
cbar_ax = fig.add_subplot(gs[0:2, 2])
cbar = fig.colorbar(sc_a, cax=cbar_ax, label='Z-score')
cbar.ax.tick_params(labelsize=8)

# Row 2: Legend (spans columns 0-1)
legend_ax = fig.add_subplot(gs[2, 0:2])
legend_ax.axis('off')
legend_patches = [
    plt.Line2D([0], [0], marker='h', color='w', markerfacecolor=REGION_COLORS['hard_stroma'],
               markersize=10, label='Hard stroma', markeredgecolor='black', markeredgewidth=0.5),
    plt.Line2D([0], [0], marker='h', color='w', markerfacecolor=REGION_COLORS['tumor_core'],
               markersize=10, label='Tumor core', markeredgecolor='black', markeredgewidth=0.5),
    plt.Line2D([0], [0], marker='h', color='w', markerfacecolor=REGION_COLORS['other'],
               markersize=10, label='Other', markeredgecolor='black', markeredgewidth=0.5),
]
legend_ax.legend(handles=legend_patches, loc='center', ncol=3, frameon=False,
                 title='Region boundaries (dashed)', title_fontsize=9)

# Add p-value annotations on each panel
def add_p_annotation(ax, pval, hs_med, tc_med):
    if pval is not None:
        sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
        text = f'HS vs TC: p={pval:.2e} {sig}\nHS median={hs_med:.2f}, TC median={tc_med:.2f}'
        ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=7,
                verticalalignment='top', bbox=dict(boxstyle='round,pad=0.3',
                facecolor='white', alpha=0.85, edgecolor='gray', linewidth=0.5))

add_p_annotation(ax_a, pdac_ccr8_p, pdac_ccr8_hs, pdac_ccr8_tc)
add_p_annotation(ax_b, pdac_mki67_p, pdac_mki67_hs, pdac_mki67_tc)
add_p_annotation(ax_c, nsclc_ccr8_p, nsclc_ccr8_hs, nsclc_ccr8_tc)
add_p_annotation(ax_d, nsclc_mki67_p, nsclc_mki67_hs, nsclc_mki67_tc)

plt.suptitle('Spatial Distribution of CCR8+ vs MKI67+ Treg Signatures', fontsize=13, fontweight='bold', y=1.01)

# Save
out_base = os.path.join(FIG_DIR, 'pub_fig_spatial_ccr8_mki67_distribution')
fig.savefig(f'{out_base}.png', dpi=600, bbox_inches='tight', facecolor='white')
fig.savefig(f'{out_base}.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{out_base}.tiff', dpi=600, bbox_inches='tight', facecolor='white', pil_kwargs={'compression': 'lzw'})
print(f"  Saved: {out_base}.png/pdf/tiff")
plt.close()

# ============================================================================
# Summary panel: Boxplot across all samples
# ============================================================================
print("\nGenerating summary boxplot...")

# Aggregate all samples
all_pdac = []
for sample in pdac['sample_id'].unique():
    sub = pdac[pdac['sample_id'] == sample]
    for region in ['hard_stroma', 'tumor_core']:
        rsub = sub[sub['region'] == region]
        if len(rsub) > 10:
            all_pdac.append({'sample': sample, 'region': region, 'ccr8': rsub['ccr8_score_z'].median(),
                            'mki67': rsub['mki67_score_z'].median()})

all_nsclc = []
for sample in nsclc['sample_id'].unique():
    sub = nsclc[nsclc['sample_id'] == sample]
    for region in ['hard_stroma', 'tumor_core']:
        rsub = sub[sub['region'] == region]
        if len(rsub) > 10:
            all_nsclc.append({'sample': sample, 'region': region, 'ccr8': rsub['ccr8_score_z'].median(),
                             'mki67': rsub['mki67_score_z'].median()})

pdac_summary = pd.DataFrame(all_pdac)
nsclc_summary = pd.DataFrame(all_nsclc)

fig, axes = plt.subplots(1, 2, figsize=(8, 4.5))

# CCR8+ boxplot
data_c = []
for _, row in pdac_summary.iterrows():
    data_c.append({'cohort': 'PDAC', 'region': row['region'], 'score': row['ccr8']})
for _, row in nsclc_summary.iterrows():
    data_c.append({'cohort': 'NSCLC', 'region': row['region'], 'score': row['ccr8']})
df_c = pd.DataFrame(data_c)

positions = [1, 2, 4, 5]
bp = axes[0].boxplot([df_c[(df_c['cohort']=='PDAC')&(df_c['region']=='hard_stroma')]['score'].values,
                      df_c[(df_c['cohort']=='PDAC')&(df_c['region']=='tumor_core')]['score'].values,
                      df_c[(df_c['cohort']=='NSCLC')&(df_c['region']=='hard_stroma')]['score'].values,
                      df_c[(df_c['cohort']=='NSCLC')&(df_c['region']=='tumor_core')]['score'].values],
                     positions=positions, widths=0.5, patch_artist=True,
                     boxprops=dict(linewidth=0.8),
                     medianprops=dict(color='black', linewidth=1.2),
                     whiskerprops=dict(linewidth=0.8),
                     capprops=dict(linewidth=0.8))
for patch, pos in zip(bp['boxes'], positions):
    patch.set_facecolor(REGION_COLORS['hard_stroma'] if pos in [1, 4] else REGION_COLORS['tumor_core'])
    patch.set_alpha(0.7)

axes[0].set_xticks([1.5, 4.5])
axes[0].set_xticklabels(['PDAC', 'NSCLC'])
axes[0].set_ylabel('Median CCR8+ Z-score')
axes[0].set_title('A. CCR8+ Treg: Hard Stroma Enrichment', fontweight='bold')
axes[0].axhline(y=0, color='black', linestyle='--', linewidth=0.5, alpha=0.4)
axes[0].set_ylim(-1.5, 2.5)

# Add significance brackets
from matplotlib.patches import FancyBboxPatch
# PDAC
pdac_hs_c = df_c[(df_c['cohort']=='PDAC')&(df_c['region']=='hard_stroma')]['score'].values
pdac_tc_c = df_c[(df_c['cohort']=='PDAC')&(df_c['region']=='tumor_core')]['score'].values
_, p_pdac_c = stats.wilcoxon(pdac_hs_c, pdac_tc_c) if len(pdac_hs_c)==len(pdac_tc_c) else stats.mannwhitneyu(pdac_hs_c, pdac_tc_c, alternative='two-sided')
# NSCLC
nsclc_hs_c = df_c[(df_c['cohort']=='NSCLC')&(df_c['region']=='hard_stroma')]['score'].values
nsclc_tc_c = df_c[(df_c['cohort']=='NSCLC')&(df_c['region']=='tumor_core')]['score'].values
_, p_nsclc_c = stats.wilcoxon(nsclc_hs_c, nsclc_tc_c) if len(nsclc_hs_c)==len(nsclc_tc_c) else stats.mannwhitneyu(nsclc_hs_c, nsclc_tc_c, alternative='two-sided')

for p, x in [(p_pdac_c, 1.5), (p_nsclc_c, 4.5)]:
    if p is not None and not np.isnan(p):
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        y = 2.2
        axes[0].plot([x-0.4, x+0.4], [y, y], 'k-', linewidth=0.8)
        axes[0].text(x, y+0.1, sig, ha='center', va='bottom', fontsize=9, fontweight='bold')

# MKI67+ boxplot
data_m = []
for _, row in pdac_summary.iterrows():
    data_m.append({'cohort': 'PDAC', 'region': row['region'], 'score': row['mki67']})
for _, row in nsclc_summary.iterrows():
    data_m.append({'cohort': 'NSCLC', 'region': row['region'], 'score': row['mki67']})
df_m = pd.DataFrame(data_m)

bp2 = axes[1].boxplot([df_m[(df_m['cohort']=='PDAC')&(df_m['region']=='hard_stroma')]['score'].values,
                       df_m[(df_m['cohort']=='PDAC')&(df_m['region']=='tumor_core')]['score'].values,
                       df_m[(df_m['cohort']=='NSCLC')&(df_m['region']=='hard_stroma')]['score'].values,
                       df_m[(df_m['cohort']=='NSCLC')&(df_m['region']=='tumor_core')]['score'].values],
                      positions=positions, widths=0.5, patch_artist=True,
                      boxprops=dict(linewidth=0.8),
                      medianprops=dict(color='black', linewidth=1.2),
                      whiskerprops=dict(linewidth=0.8),
                      capprops=dict(linewidth=0.8))
for patch, pos in zip(bp2['boxes'], positions):
    patch.set_facecolor(REGION_COLORS['hard_stroma'] if pos in [1, 4] else REGION_COLORS['tumor_core'])
    patch.set_alpha(0.7)

axes[1].set_xticks([1.5, 4.5])
axes[1].set_xticklabels(['PDAC', 'NSCLC'])
axes[1].set_ylabel('Median MKI67+ Z-score')
axes[1].set_title('B. MKI67+ Treg: No Regional Preference', fontweight='bold')
axes[1].axhline(y=0, color='black', linestyle='--', linewidth=0.5, alpha=0.4)
axes[1].set_ylim(-1.5, 2.5)

# Add legend
legend_patches2 = [
    plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=REGION_COLORS['hard_stroma'],
               markersize=10, label='Hard stroma', markeredgecolor='black', markeredgewidth=0.5),
    plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=REGION_COLORS['tumor_core'],
               markersize=10, label='Tumor core', markeredgecolor='black', markeredgewidth=0.5),
]
axes[1].legend(handles=legend_patches2, loc='upper right', frameon=True, fancybox=False,
               edgecolor='black', fontsize=8)

plt.suptitle('Cross-cohort Regional Comparison of Treg Signatures', fontsize=12, fontweight='bold', y=1.02)
plt.tight_layout()

out_base2 = os.path.join(FIG_DIR, 'pub_fig_spatial_region_boxplot')
fig.savefig(f'{out_base2}.png', dpi=600, bbox_inches='tight', facecolor='white')
fig.savefig(f'{out_base2}.pdf', dpi=300, bbox_inches='tight', facecolor='white')
print(f"  Saved: {out_base2}.png/pdf")
plt.close()

print("\nAll spatial figures generated!")
