"""
GSE189487 Visium spatial transcriptomics analysis.
Early-stage LUAD (AIS -> MIA -> IAC) spatial validation of Contractile Treg signature.
"""

import os, gzip, warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import io, stats
import matplotlib.pyplot as plt

DATA_DIR = r'D:\Research\tomato\data'
RESULT_DIR = r'D:\Research\tomato\results\gse189487_spatial'
FIG_DIR = r'D:\Research\tomato\figures'
os.makedirs(RESULT_DIR, exist_ok=True)

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

# Sample info
SAMPLES = {
    'TD1': ('GSM5702473', 'IAC'),
    'TD2': ('GSM5702474', 'IAC'),
    'TD3': ('GSM5702475', 'MIA'),
    'TD5': ('GSM5702476', 'AIS'),
    'TD6': ('GSM5702477', 'MIA'),
    'TD8': ('GSM5702478', 'AIS'),
}

STAGE_ORDER = ['AIS', 'MIA', 'IAC']
STAGE_COLORS = {'AIS': '#3498db', 'MIA': '#f39c12', 'IAC': '#e74c3c'}

# Signature genes
CORE_MARKERS = ['ACTA2', 'MYL9', 'VCL', 'TLN1', 'TRPV4', 'MOB1B']
EXTENDED = ['ITGA1', 'DDR1', 'DDR2', 'SRC', 'ITGB1', 'CD44', 'PTK2', 'ZYX', 'FLNA']
TREG_ANCHOR = ['FOXP3', 'IL2RA', 'CTLA4']
SIGNATURE = CORE_MARKERS + EXTENDED + TREG_ANCHOR

ECM_GENES = ['COL1A1', 'COL1A2', 'COL3A1', 'COL6A1', 'FN1', 'POSTN', 'THBS1', 'VIM']

print("="*60)
print("Loading GSE189487 Visium data...")
print("="*60)

adatas = []
for sample_name, (gsm, stage) in SAMPLES.items():
    print(f"\nLoading {sample_name} ({stage})...")
    prefix = f"{gsm}_{sample_name}"
    
    # Read matrix
    matrix_path = os.path.join(DATA_DIR, f"{prefix}_matrix.mtx.gz")
    features_path = os.path.join(DATA_DIR, f"{prefix}_features.tsv.gz")
    barcodes_path = os.path.join(DATA_DIR, f"{prefix}_barcodes.tsv.gz")
    positions_path = os.path.join(DATA_DIR, f"{prefix}_tissue_positions_list.csv.gz")
    
    # Load count matrix
    with gzip.open(matrix_path, 'rb') as f:
        mat = io.mmread(f).T.tocsr()
    
    # Load features (genes)
    features = pd.read_csv(features_path, sep='\t', header=None)
    if features.shape[1] >= 2:
        gene_names = features.iloc[:, 1].values
    else:
        gene_names = features.iloc[:, 0].values
    
    # Load barcodes (spots)
    barcodes = pd.read_csv(barcodes_path, sep='\t', header=None)
    spot_names = barcodes.iloc[:, 0].values
    
    # Create AnnData
    adata = sc.AnnData(X=mat)
    adata.obs_names = spot_names
    adata.var_names = gene_names
    adata.var_names_make_unique()
    
    # Load spatial coordinates
    positions = pd.read_csv(positions_path, header=None)
    # Format: barcode, in_tissue, array_row, array_col, pxl_col_in_fullres, pxl_row_in_fullres
    positions.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_col', 'pxl_row']
    positions = positions.set_index('barcode')
    
    # Match barcodes
    matched_pos = positions.loc[adata.obs_names]
    adata.obs['array_row'] = pd.to_numeric(matched_pos['array_row'].values, errors='coerce')
    adata.obs['array_col'] = pd.to_numeric(matched_pos['array_col'].values, errors='coerce')
    adata.obs['pxl_row'] = pd.to_numeric(matched_pos['pxl_row'].values, errors='coerce')
    adata.obs['pxl_col'] = pd.to_numeric(matched_pos['pxl_col'].values, errors='coerce')
    adata.obs['in_tissue'] = pd.to_numeric(matched_pos['in_tissue'].values, errors='coerce')
    adata.obs['sample'] = sample_name
    adata.obs['stage'] = stage
    
    # Filter to tissue spots only
    adata = adata[adata.obs['in_tissue'] == 1].copy()
    
    print(f"  Spots: {adata.n_obs}, Genes: {adata.n_vars}")
    adatas.append(adata)

# ===================== MERGE =====================
print("\n" + "="*60)
print("Merging samples...")
print("="*60)

# Use outer join for genes
merged = sc.concat(adatas, join='outer', label='sample', index_unique='-')
merged.X = np.nan_to_num(merged.X, nan=0)
print(f"Merged: {merged.n_obs} spots x {merged.n_vars} genes")

# Fix data types after concat
for col in ['array_row', 'array_col', 'pxl_row', 'pxl_col', 'in_tissue']:
    merged.obs[col] = pd.to_numeric(merged.obs[col], errors='coerce')

# Make var names unique
merged.var_names_make_unique()

# ===================== QC & NORMALIZE =====================
print("\nNormalizing...")
sc.pp.normalize_total(merged, target_sum=1e4)
sc.pp.log1p(merged)

# ===================== COMPUTE SIGNATURE SCORES =====================
print("\n" + "="*60)
print("Computing signature scores...")
print("="*60)

# Check gene availability
sig_available = [g for g in SIGNATURE if g in merged.var_names]
sig_missing = [g for g in SIGNATURE if g not in merged.var_names]
ecm_available = [g for g in ECM_GENES if g in merged.var_names]
ecm_missing = [g for g in ECM_GENES if g not in merged.var_names]

print(f"Signature genes: {len(sig_available)}/{len(SIGNATURE)} available")
print(f"  Available: {sig_available}")
print(f"  Missing: {sig_missing}")
print(f"ECM genes: {len(ecm_available)}/{len(ECM_GENES)} available")
print(f"  Missing: {ecm_missing}")

# Compute scores using scanpy's score_genes (z-score based)
if len(sig_available) >= 3:
    sc.tl.score_genes(merged, sig_available, score_name='contractile_treg_score', use_raw=False)
    print(f"Contractile score computed")

if len(ecm_available) >= 3:
    sc.tl.score_genes(merged, ecm_available, score_name='ecm_score', use_raw=False)
    print(f"ECM score computed")

# Also compute per-gene expression for key markers
for gene in ['COL1A2', 'POSTN', 'COL1A1', 'FN1', 'ACTA2']:
    if gene in merged.var_names:
        merged.obs[f'{gene}_expr'] = merged[:, gene].X.toarray().flatten()

# ===================== SPATIAL ANALYSIS =====================
print("\n" + "="*60)
print("Spatial analysis...")
print("="*60)

# Stage-level summary
stage_summary = merged.obs.groupby('stage').agg({
    'contractile_treg_score': ['mean', 'std', 'count'],
    'ecm_score': ['mean', 'std'],
})
print("\nStage-level summary:")
print(stage_summary)

# Pairwise comparisons
print("\n--- Contractile score by stage ---")
for stage in STAGE_ORDER:
    vals = merged.obs[merged.obs['stage'] == stage]['contractile_treg_score'].values
    print(f"  {stage}: {vals.mean():.4f} ± {vals.std():.4f} (n={len(vals)})")

# ANOVA
stage_vals = [merged.obs[merged.obs['stage'] == s]['contractile_treg_score'].values for s in STAGE_ORDER]
f_stat, p_anova = stats.f_oneway(*stage_vals)
print(f"\nANOVA: F={f_stat:.3f}, p={p_anova:.4e}")

# Post-hoc: AIS vs IAC
ais_vals = merged.obs[merged.obs['stage'] == 'AIS']['contractile_treg_score'].values
iac_vals = merged.obs[merged.obs['stage'] == 'IAC']['contractile_treg_score'].values
t_stat, p_ttest = stats.ttest_ind(ais_vals, iac_vals)
print(f"AIS vs IAC t-test: t={t_stat:.3f}, p={p_ttest:.4e}")

# ECM-Treg correlation
if 'contractile_treg_score' in merged.obs.columns and 'ecm_score' in merged.obs.columns:
    r, p_corr = stats.pearsonr(merged.obs['ecm_score'], merged.obs['contractile_treg_score'])
    print(f"\nECM score vs Contractile score: r={r:.3f}, p={p_corr:.2e}")

# ===================== VISUALIZATION =====================
print("\nGenerating figures...")

fig = plt.figure(figsize=(14, 10), dpi=300)

# Helper to plot spatial
from matplotlib.colors import LinearSegmentedColormap

def plot_spatial(ax, adata_subset, color_key, cmap='viridis', title='', vmin=None, vmax=None):
    x = adata_subset.obs['pxl_col'].values.astype(float)
    y = adata_subset.obs['pxl_row'].values.astype(float)
    c = adata_subset.obs[color_key].values.astype(float)
    
    if len(y) == 0:
        ax.set_title(title + ' (empty)', fontsize=9)
        return None
    
    # Flip y for image-like orientation
    y = y.max() - y
    
    scatter = ax.scatter(x, y, c=c, s=8, cmap=cmap, vmin=vmin, vmax=vmax, edgecolors='none')
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=9, fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])
    return scatter

# Panel 1-6: Contractile score per sample
for i, (sample_name, (gsm, stage)) in enumerate(SAMPLES.items()):
    ax = fig.add_subplot(3, 4, i+1)
    sub = merged[merged.obs['sample'] == sample_name]
    vmax = merged.obs['contractile_treg_score'].quantile(0.98)
    scat = plot_spatial(ax, sub, 'contractile_treg_score', cmap='RdYlBu_r', 
                        title=f'{sample_name} ({stage})', vmax=vmax)
    if scat is not None:
        plt.colorbar(scat, ax=ax, fraction=0.046, pad=0.04)
    else:
        print(f"  Warning: empty plot for {sample_name}")

# Panel 7: Stage comparison boxplot
ax7 = fig.add_subplot(3, 4, 7)
plot_data = [merged.obs[merged.obs['stage'] == s]['contractile_treg_score'].values for s in STAGE_ORDER]
bp = ax7.boxplot(plot_data, labels=STAGE_ORDER, patch_artist=True,
                 boxprops=dict(edgecolor='black', linewidth=0.5),
                 medianprops=dict(color='black', linewidth=1.5),
                 whiskerprops=dict(linewidth=0.5),
                 capprops=dict(linewidth=0.5))
for patch, stage in zip(bp['boxes'], STAGE_ORDER):
    patch.set_facecolor(STAGE_COLORS[stage])
    patch.set_alpha(0.7)
ax7.set_ylabel('Contractile Treg score')
ax7.set_title('By stage', fontweight='bold')

# Panel 8: ECM score per stage
ax8 = fig.add_subplot(3, 4, 8)
plot_data2 = [merged.obs[merged.obs['stage'] == s]['ecm_score'].values for s in STAGE_ORDER]
bp2 = ax8.boxplot(plot_data2, labels=STAGE_ORDER, patch_artist=True,
                  boxprops=dict(edgecolor='black', linewidth=0.5),
                  medianprops=dict(color='black', linewidth=1.5),
                  whiskerprops=dict(linewidth=0.5),
                  capprops=dict(linewidth=0.5))
for patch, stage in zip(bp2['boxes'], STAGE_ORDER):
    patch.set_facecolor(STAGE_COLORS[stage])
    patch.set_alpha(0.7)
ax8.set_ylabel('ECM score')
ax8.set_title('ECM by stage', fontweight='bold')

# Panel 9: ECM vs Contractile scatter
ax9 = fig.add_subplot(3, 4, 9)
for stage in STAGE_ORDER:
    sub = merged.obs[merged.obs['stage'] == stage]
    ax9.scatter(sub['ecm_score'], sub['contractile_treg_score'], 
                c=STAGE_COLORS[stage], alpha=0.3, s=5, label=stage, edgecolors='none')
ax9.set_xlabel('ECM score')
ax9.set_ylabel('Contractile Treg score')
ax9.set_title(f'ECM vs Contractile\nr={r:.3f}, p={p_corr:.2e}', fontweight='bold')
ax9.legend(fontsize=7, frameon=False)

# Panel 10: COL1A2 expression spatial (example: TD2 IAC)
ax10 = fig.add_subplot(3, 4, 10)
sub_iac = merged[merged.obs['sample'] == 'TD2']
if 'COL1A2_expr' in sub_iac.obs.columns:
    vmax_col = sub_iac.obs['COL1A2_expr'].quantile(0.99)
    scat10 = plot_spatial(ax10, sub_iac, 'COL1A2_expr', cmap='YlOrRd', 
                          title='TD2 (IAC): COL1A2', vmax=vmax_col)
    if scat10 is not None:
        plt.colorbar(scat10, ax=ax10, fraction=0.046, pad=0.04)

# Panel 11: POSTN expression spatial
ax11 = fig.add_subplot(3, 4, 11)
if 'POSTN_expr' in sub_iac.obs.columns:
    vmax_postn = sub_iac.obs['POSTN_expr'].quantile(0.99)
    scat11 = plot_spatial(ax11, sub_iac, 'POSTN_expr', cmap='YlOrRd',
                          title='TD2 (IAC): POSTN', vmax=vmax_postn)
    if scat11 is not None:
        plt.colorbar(scat11, ax=ax11, fraction=0.046, pad=0.04)

# Panel 12: Summary stats table
ax12 = fig.add_subplot(3, 4, 12)
ax12.axis('off')
summary_text = f"""
GSE189487 Spatial Analysis Summary

Samples: 6 Visium sections
  AIS: TD5, TD8 (n=2)
  MIA: TD3, TD6 (n=2)
  IAC: TD1, TD2 (n=2)

Signature genes matched: {len(sig_available)}/{len(SIGNATURE)}
ECM genes matched: {len(ecm_available)}/{len(ECM_GENES)}

ANOVA (stage): F={f_stat:.2f}, p={p_anova:.2e}
AIS vs IAC: t={t_stat:.2f}, p={p_ttest:.2e}
ECM-Contractile: r={r:.3f}, p={p_corr:.2e}

Total spots: {merged.n_obs}
"""
ax12.text(0.1, 0.9, summary_text, transform=ax12.transAxes, fontsize=8,
          verticalalignment='top', fontfamily='monospace',
          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'gse189487_spatial_validation.png'), dpi=600, bbox_inches='tight')
fig.savefig(os.path.join(FIG_DIR, 'gse189487_spatial_validation.pdf'), dpi=600, bbox_inches='tight')
plt.close()

# Save results
merged.obs.to_csv(os.path.join(RESULT_DIR, 'spot_level_scores.csv'))
print(f"\nResults saved to {RESULT_DIR}")
print(f"Figure saved: gse189487_spatial_validation.png/pdf")

print("\n" + "="*60)
print("GSE189487 SPATIAL ANALYSIS COMPLETE")
print("="*60)
