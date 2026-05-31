"""
Re-annotate Treg subtypes from scRNA-seq data.
Find robust marker genes for CCR8+ and MKI67+ Treg using DE analysis,
then rebuild signatures and validate on spatial data.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = r'D:\Research\tomato'
OUT_DIR = os.path.join(BASE_DIR, 'results', 'treg_reannotation')
os.makedirs(OUT_DIR, exist_ok=True)

# ================================================================
# 1. Load scRNA-seq Treg data
# ================================================================
print("Loading Treg scRNA-seq data...")
ad = sc.read(os.path.join(BASE_DIR, 'data', 'monocle3_input', 'treg_processed.h5ad'))
print(f"  {ad.n_obs} cells x {ad.n_vars} genes")
print(f"  Subtypes: {ad.obs['sub_cell_type'].value_counts().to_dict()}")

# Basic preprocessing if needed
if 'log1p' not in ad.uns:
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

# ================================================================
# 2. Differential expression: CCR8+ vs FOXP3+ (rest)
# ================================================================
print("\n" + "="*60)
print("DE: CCR8+ vs FOXP3+ Treg")
print("="*60)

# Subset to CCR8+ and FOXP3+
ccr8_mask = ad.obs['sub_cell_type'] == 'CD4T_Treg_CCR8'
foxp3_mask = ad.obs['sub_cell_type'] == 'CD4T_Treg_FOXP3'
ccr8_cells = ad[ccr8_mask].copy()
foxp3_cells = ad[foxp3_mask].copy()

print(f"  CCR8+ cells: {ccr8_cells.n_obs}")
print(f"  FOXP3+ cells: {foxp3_cells.n_obs}")

# Compute per-gene stats
genes = ad.var_names
de_results = []

for gene in genes:
    gidx = ad.var_names.get_loc(gene)
    ccr8_expr = ccr8_cells.X[:, gidx].toarray().ravel()
    foxp3_expr = foxp3_cells.X[:, gidx].toarray().ravel()
    
    # Pct expression
    ccr8_pct = (ccr8_expr > 0).mean()
    foxp3_pct = (foxp3_expr > 0).mean()
    
    # Mean expression
    ccr8_mean = ccr8_expr.mean()
    foxp3_mean = foxp3_expr.mean()
    
    # Wilcoxon test
    if ccr8_mean > 0 or foxp3_mean > 0:
        stat, pval = stats.mannwhitneyu(ccr8_expr, foxp3_expr, alternative='two-sided')
    else:
        pval = 1.0
    
    log2fc = np.log2((ccr8_mean + 0.01) / (foxp3_mean + 0.01))
    
    de_results.append({
        'gene': gene,
        'ccr8_mean': ccr8_mean, 'foxp3_mean': foxp3_mean,
        'ccr8_pct': ccr8_pct, 'foxp3_pct': foxp3_pct,
        'log2fc': log2fc, 'pvalue': pval,
    })

df_ccr8 = pd.DataFrame(de_results)
df_ccr8['p_adj'] = np.minimum(df_ccr8['pvalue'] * len(df_ccr8), 1.0)  # Bonferroni
df_ccr8 = df_ccr8.sort_values('log2fc', ascending=False)

# Strict filtering for CCR8+ markers
# Up in CCR8+: log2FC > 1, expressed in >30% of CCR8+ cells, p_adj < 0.01
ccr8_up = df_ccr8[(df_ccr8['log2fc'] > 1.0) & 
                    (df_ccr8['ccr8_pct'] > 0.3) & 
                    (df_ccr8['p_adj'] < 0.01)]
print(f"\n  CCR8+ UP markers (strict): {len(ccr8_up)} genes")
print("  Top 20:")
print(ccr8_up.head(20)[['gene', 'ccr8_mean', 'foxp3_mean', 'log2fc', 'ccr8_pct', 'p_adj']].to_string(index=False))

# Also check CCR8 itself
if 'CCR8' in df_ccr8['gene'].values:
    print(f"\n  CCR8 gene itself:")
    print(df_ccr8[df_ccr8['gene']=='CCR8'][['gene','ccr8_mean','foxp3_mean','log2fc','ccr8_pct','p_adj']].to_string(index=False))

# ================================================================
# 3. Differential expression: MKI67+ vs FOXP3+ (rest)
# ================================================================
print("\n" + "="*60)
print("DE: MKI67+ vs FOXP3+ Treg")
print("="*60)

mki67_mask = ad.obs['sub_cell_type'] == 'CD4T_Treg_MKI67'
mki67_cells = ad[mki67_mask].copy()
print(f"  MKI67+ cells: {mki67_cells.n_obs}")
print(f"  FOXP3+ cells: {foxp3_cells.n_obs}")

de_results_mki = []
for gene in genes:
    gidx = ad.var_names.get_loc(gene)
    mki_expr = mki67_cells.X[:, gidx].toarray().ravel()
    foxp3_expr = foxp3_cells.X[:, gidx].toarray().ravel()
    
    mki_pct = (mki_expr > 0).mean()
    foxp3_pct = (foxp3_expr > 0).mean()
    mki_mean = mki_expr.mean()
    foxp3_mean = foxp3_expr.mean()
    
    if mki_mean > 0 or foxp3_mean > 0:
        stat, pval = stats.mannwhitneyu(mki_expr, foxp3_expr, alternative='two-sided')
    else:
        pval = 1.0
    
    log2fc = np.log2((mki_mean + 0.01) / (foxp3_mean + 0.01))
    
    de_results_mki.append({
        'gene': gene,
        'mki67_mean': mki_mean, 'foxp3_mean': foxp3_mean,
        'mki67_pct': mki_pct, 'foxp3_pct': foxp3_pct,
        'log2fc': log2fc, 'pvalue': pval,
    })

df_mki = pd.DataFrame(de_results_mki)
df_mki['p_adj'] = np.minimum(df_mki['pvalue'] * len(df_mki), 1.0)
df_mki = df_mki.sort_values('log2fc', ascending=False)

# Strict filtering for MKI67+ markers
mki_up = df_mki[(df_mki['log2fc'] > 1.0) & 
                (df_mki['mki67_pct'] > 0.3) & 
                (df_mki['p_adj'] < 0.01)]
print(f"\n  MKI67+ UP markers (strict): {len(mki_up)} genes")
print("  Top 20:")
print(mki_up.head(20)[['gene', 'mki67_mean', 'foxp3_mean', 'log2fc', 'mki67_pct', 'p_adj']].to_string(index=False))

# Check MKI67 itself
if 'MKI67' in df_mki['gene'].values:
    print(f"\n  MKI67 gene itself:")
    print(df_mki[df_mki['gene']=='MKI67'][['gene','mki67_mean','foxp3_mean','log2fc','mki67_pct','p_adj']].to_string(index=False))

# ================================================================
# 4. Save results
# ================================================================
print("\n" + "="*60)
print("Saving results")
print("="*60)

# Save all DE results
df_ccr8.to_csv(os.path.join(OUT_DIR, 'de_ccr8_vs_foxp3.csv'), index=False)
df_mki.to_csv(os.path.join(OUT_DIR, 'de_mki67_vs_foxp3.csv'), index=False)

# Save strict markers
ccr8_markers = ccr8_up['gene'].tolist()
mki_markers = mki_up['gene'].tolist()

pd.DataFrame({'CCR8+_markers': pd.Series(ccr8_markers), 
              'MKI67+_markers': pd.Series(mki_markers)}).to_csv(
    os.path.join(OUT_DIR, 'strict_markers.csv'), index=False)

print(f"  CCR8+ strict markers: {len(ccr8_markers)} genes -> {os.path.join(OUT_DIR, 'strict_markers.csv')}")
print(f"  MKI67+ strict markers: {len(mki_markers)} genes")

# Print overlap with old signatures
print("\n  Comparison with old signatures:")
old_sig = pd.read_csv(os.path.join(BASE_DIR, 'results', 'spatial', 'treg_subtype_signatures.csv'))
old_ccr8 = set(old_sig['CD4T_Treg_CCR8'].dropna())
old_mki = set(old_sig['CD4T_Treg_MKI67'].dropna())
print(f"    Old CCR8+ signature: {len(old_ccr8)} genes")
print(f"    New CCR8+ markers overlap: {len(set(ccr8_markers) & old_ccr8)} / {len(old_ccr8)}")
print(f"    Old MKI67+ signature: {len(old_mki)} genes")
print(f"    New MKI67+ markers overlap: {len(set(mki_markers) & old_mki)} / {len(old_mki)}")

print(f"\nAll results in: {OUT_DIR}")
