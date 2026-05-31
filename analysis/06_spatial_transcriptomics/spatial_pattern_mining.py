"""
Spatial Transcriptomics Pattern Mining
Continue finding patterns in spatial Treg distribution across cohorts.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import cdist

OUT_DIR = 'D:/Research/tomato/results/spatial_pattern_mining'
FIG_DIR = 'D:/Research/tomato/figures/spatial_pattern_mining'
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 12,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
    'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
})

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading spatial metadata...")

pdac = pd.read_csv('D:/Research/tomato/results/spatial/spatial_spot_metadata.csv')
nsclc = pd.read_csv('D:/Research/tomato/results/spatial_nsclc/nsclc_spot_metadata.csv')

print(f"PDAC: {pdac.shape[0]} spots")
print(f"NSCLC: {nsclc.shape[0]} spots")

# ---------------------------------------------------------------------------
# Analysis 1: Grade-associated trend (PDAC only)
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("Analysis 1: Grade-associated CCR8+/MKI67+ trend")
print("="*60)

grade_order = ['Low Grade', 'High Grade', 'PDAC']
grade_summary = []
for grade in grade_order:
    sub = pdac[pdac['grade'] == grade]
    for region in ['hard_stroma', 'tumor_core', 'other']:
        rsub = sub[sub['region'] == region]
        if len(rsub) > 10:
            grade_summary.append({
                'grade': grade,
                'region': region,
                'n': len(rsub),
                'ccr8_z_mean': rsub['ccr8_score_z'].mean(),
                'ccr8_z_std': rsub['ccr8_score_z'].std(),
                'mki67_z_mean': rsub['mki67_score_z'].mean(),
                'mki67_z_std': rsub['mki67_score_z'].std(),
            })

grade_df = pd.DataFrame(grade_summary)
print(grade_df.to_string(index=False))

# Grade trend within hard_stroma
hs_by_grade = pdac[pdac['region'] == 'hard_stroma'].groupby('grade').agg({
    'ccr8_score_z': ['mean', 'std', 'count'],
    'mki67_score_z': ['mean', 'std', 'count']
}).reindex(grade_order)
print("\nHard stroma by grade:")
print(hs_by_grade)

# ANOVA for grade effect in hard_stroma
hs_data = pdac[pdac['region'] == 'hard_stroma']
lg = hs_data[hs_data['grade'] == 'Low Grade']['ccr8_score_z']
hg = hs_data[hs_data['grade'] == 'High Grade']['ccr8_score_z']
pdac_g = hs_data[hs_data['grade'] == 'PDAC']['ccr8_score_z']
f_stat, p_anova = stats.f_oneway(lg.dropna(), hg.dropna(), pdac_g.dropna())
print(f"\nANOVA (CCR8+ in hard_stroma across grades): F={f_stat:.3f}, p={p_anova:.4f}")

# ---------------------------------------------------------------------------
# Analysis 2: Spot-level CCR8 vs MKI67 correlation by region
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("Analysis 2: Spot-level CCR8 vs MKI67 correlation")
print("="*60)

corr_results = []
for cohort, df in [('PDAC', pdac), ('NSCLC', nsclc)]:
    for region in ['hard_stroma', 'tumor_core', 'other']:
        sub = df[df['region'] == region]
        if len(sub) > 30:
            rho, pval = stats.spearmanr(sub['ccr8_score_z'], sub['mki67_score_z'])
            corr_results.append({
                'cohort': cohort,
                'region': region,
                'n': len(sub),
                'rho': rho,
                'p_value': pval,
            })
            print(f"  {cohort} {region}: n={len(sub)}, rho={rho:.4f}, p={pval:.2e}")

corr_df = pd.DataFrame(corr_results)

# ---------------------------------------------------------------------------
# Analysis 3: Spatial mutual exclusion (per sample)
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("Analysis 3: Spatial mutual exclusion (CCR8+ vs MKI67+ spots)")
print("="*60)

def analyze_mutual_exclusion(df, cohort_name):
    """For each sample, check if high-CCR8 spots and high-MKI67 spots are spatially separated."""
    results = []
    samples = df['sample_id'].unique()
    
    for sample in samples:
        sub = df[df['sample_id'] == sample].copy()
        if len(sub) < 50:
            continue
        
        # Define high-CCR8 and high-MKI67 spots (top 20% per sample)
        ccr8_q80 = sub['ccr8_score_z'].quantile(0.8)
        mki67_q80 = sub['mki67_score_z'].quantile(0.8)
        
        high_ccr8 = sub[sub['ccr8_score_z'] >= ccr8_q80]
        high_mki67 = sub[sub['mki67_score_z'] >= mki67_q80]
        
        if len(high_ccr8) < 5 or len(high_mki67) < 5:
            continue
        
        # Calculate mean distance between high-CCR8 and high-MKI67 spots
        coords_ccr8 = high_ccr8[['array_row', 'array_col']].values
        coords_mki67 = high_mki67[['array_row', 'array_col']].values
        
        dist_matrix = cdist(coords_ccr8, coords_mki67, metric='euclidean')
        mean_dist = dist_matrix.mean()
        min_dist = dist_matrix.min()
        
        # Permutation test: shuffle labels and recalculate
        np.random.seed(42)
        n_perm = 100
        perm_dists = []
        all_coords = sub[['array_row', 'array_col']].values
        n_ccr8 = len(high_ccr8)
        n_mki67 = len(high_mki67)
        
        for _ in range(n_perm):
            idx = np.random.permutation(len(all_coords))
            perm_ccr8 = all_coords[idx[:n_ccr8]]
            perm_mki67 = all_coords[idx[n_ccr8:n_ccr8+n_mki67]]
            perm_dist = cdist(perm_ccr8, perm_mki67).mean()
            perm_dists.append(perm_dist)
        
        # If real mean_dist > perm_dist, they are more separated than random
        z_score = (mean_dist - np.mean(perm_dists)) / np.std(perm_dists)
        p_val = 1 - stats.norm.cdf(z_score)  # one-sided: test if more separated
        
        results.append({
            'cohort': cohort_name,
            'sample': sample,
            'n_spots': len(sub),
            'high_ccr8_n': len(high_ccr8),
            'high_mki67_n': len(high_mki67),
            'mean_dist': mean_dist,
            'min_dist': min_dist,
            'perm_mean': np.mean(perm_dists),
            'perm_std': np.std(perm_dists),
            'z_score': z_score,
            'p_value': p_val,
            'significant': p_val < 0.05,
        })
    
    return pd.DataFrame(results)

me_pdac = analyze_mutual_exclusion(pdac, 'PDAC')
me_nsclc = analyze_mutual_exclusion(nsclc, 'NSCLC')
me_all = pd.concat([me_pdac, me_nsclc], ignore_index=True)

print("\nMutual exclusion results:")
print(me_all[['cohort', 'sample', 'mean_dist', 'perm_mean', 'z_score', 'p_value', 'significant']].to_string(index=False))

# Summary
n_sig = me_all['significant'].sum()
n_total = len(me_all)
print(f"\nSummary: {n_sig}/{n_total} samples show significant spatial separation (p<0.05)")
print(f"  PDAC: {me_pdac['significant'].sum()}/{len(me_pdac)}")
print(f"  NSCLC: {me_nsclc['significant'].sum()}/{len(me_nsclc)}")

# ---------------------------------------------------------------------------
# Analysis 4: Transition zone gradient (PDAC)
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("Analysis 4: Transition zone gradient")
print("="*60)

def analyze_transition_gradient(df, cohort_name):
    """Analyze gradient from tumor_core -> hard_stroma through 'other' region."""
    results = []
    samples = df['sample_id'].unique()
    
    for sample in samples:
        sub = df[df['sample_id'] == sample].copy()
        if len(sub) < 50:
            continue
        
        # Only analyze samples with all three regions
        regions = sub['region'].unique()
        if not all(r in regions for r in ['hard_stroma', 'tumor_core', 'other']):
            continue
        
        # Calculate distance to nearest hard_stroma and tumor_core spot
        hs_coords = sub[sub['region'] == 'hard_stroma'][['array_row', 'array_col']].values
        tc_coords = sub[sub['region'] == 'tumor_core'][['array_row', 'array_col']].values
        
        if len(hs_coords) < 5 or len(tc_coords) < 5:
            continue
        
        all_coords = sub[['array_row', 'array_col']].values
        dist_to_hs = np.min(cdist(all_coords, hs_coords), axis=1)
        dist_to_tc = np.min(cdist(all_coords, tc_coords), axis=1)
        
        sub['dist_to_hs'] = dist_to_hs
        sub['dist_to_tc'] = dist_to_tc
        
        # Transition score: negative = closer to HS, positive = closer to TC
        sub['transition_score'] = dist_to_tc - dist_to_hs
        
        # Correlation with transition score
        rho_ccr8, p_ccr8 = stats.spearmanr(sub['transition_score'], sub['ccr8_score_z'])
        rho_mki67, p_mki67 = stats.spearmanr(sub['transition_score'], sub['mki67_score_z'])
        
        results.append({
            'cohort': cohort_name,
            'sample': sample,
            'n': len(sub),
            'rho_ccr8': rho_ccr8,
            'p_ccr8': p_ccr8,
            'rho_mki67': rho_mki67,
            'p_mki67': p_mki67,
        })
    
    return pd.DataFrame(results)

trans_pdac = analyze_transition_gradient(pdac, 'PDAC')
trans_nsclc = analyze_transition_gradient(nsclc, 'NSCLC')
trans_all = pd.concat([trans_pdac, trans_nsclc], ignore_index=True)

print("\nTransition gradient results:")
print(trans_all.to_string(index=False))

# Count significant gradients
ccr8_grad = (trans_all['p_ccr8'] < 0.05) & (trans_all['rho_ccr8'] < 0)  # CCR8 decreases toward TC
mki67_grad = (trans_all['p_mki67'] < 0.05) & (trans_all['rho_mki67'] > 0)  # MKI67 increases toward TC
print(f"\nCCR8+ gradient (higher in HS, lower in TC): {ccr8_grad.sum()}/{len(trans_all)} samples")
print(f"MKI67+ gradient (higher in TC, lower in HS): {mki67_grad.sum()}/{len(trans_all)} samples")

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("Generating figures...")
print("="*60)

# Fig 1: Grade trend
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

grade_pivot = grade_df.pivot(index='grade', columns='region', values='ccr8_z_mean').reindex(grade_order)
grade_pivot[['hard_stroma', 'tumor_core', 'other']].plot(kind='bar', ax=axes[0], color=['#E74C3C', '#3498DB', '#95A5A6'], alpha=0.8, edgecolor='black')
axes[0].set_title('A. CCR8+ Z-score by Grade & Region', fontweight='bold')
axes[0].set_ylabel('CCR8+ Z-score')
axes[0].set_xlabel('Grade')
axes[0].legend(title='Region')
axes[0].axhline(y=0, color='black', linestyle='--', linewidth=0.5)

grade_pivot2 = grade_df.pivot(index='grade', columns='region', values='mki67_z_mean').reindex(grade_order)
grade_pivot2[['hard_stroma', 'tumor_core', 'other']].plot(kind='bar', ax=axes[1], color=['#E74C3C', '#3498DB', '#95A5A6'], alpha=0.8, edgecolor='black')
axes[1].set_title('B. MKI67+ Z-score by Grade & Region', fontweight='bold')
axes[1].set_ylabel('MKI67+ Z-score')
axes[1].set_xlabel('Grade')
axes[1].legend(title='Region')
axes[1].axhline(y=0, color='black', linestyle='--', linewidth=0.5)

plt.suptitle('PDAC Spatial Treg Signatures by Tumor Grade', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/fig1_grade_trend.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{FIG_DIR}/fig1_grade_trend.pdf', bbox_inches='tight')
plt.close()
print("  Saved: fig1_grade_trend.png")

# Fig 2: Correlation heatmap
fig, ax = plt.subplots(figsize=(8, 5))
pivot_corr = corr_df.pivot(index='region', columns='cohort', values='rho')
pivot_corr = pivot_corr.reindex(['hard_stroma', 'other', 'tumor_core'])
sns.heatmap(pivot_corr, annot=True, fmt='.3f', cmap='RdBu_r', center=0, vmin=-0.5, vmax=0.5,
            square=True, linewidths=0.5, cbar_kws={'shrink': 0.8, 'label': 'Spearman ρ'}, ax=ax)
ax.set_title('CCR8+ vs MKI67+ Spot-level Correlation', fontweight='bold')
ax.set_xlabel('Cohort')
ax.set_ylabel('Region')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/fig2_correlation_heatmap.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{FIG_DIR}/fig2_correlation_heatmap.pdf', bbox_inches='tight')
plt.close()
print("  Saved: fig2_correlation_heatmap.png")

# Fig 3: Mutual exclusion
fig, ax = plt.subplots(figsize=(10, 6))
plot_df = me_all.sort_values('z_score')
colors = ['#E74C3C' if s else '#95A5A6' for s in plot_df['significant']]
bars = ax.barh(range(len(plot_df)), plot_df['z_score'], color=colors, alpha=0.8, edgecolor='black')
ax.set_yticks(range(len(plot_df)))
ax.set_yticklabels([f"{r['cohort']}:{r['sample']}" for _, r in plot_df.iterrows()], fontsize=8)
ax.set_xlabel('Spatial Separation Z-score', fontsize=11)
ax.set_title('CCR8+ vs MKI67+ Spatial Mutual Exclusion\n(Higher = more spatially separated)', fontweight='bold')
ax.axvline(x=0, color='black', linewidth=0.5)
ax.axvline(x=1.645, color='red', linestyle='--', linewidth=0.8, label='p=0.05')
ax.legend()
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/fig3_mutual_exclusion.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{FIG_DIR}/fig3_mutual_exclusion.pdf', bbox_inches='tight')
plt.close()
print("  Saved: fig3_mutual_exclusion.png")

# Fig 4: Transition gradient
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# CCR8
colors_c = ['#8bc98b' if p < 0.05 and r < 0 else '#95A5A6' for r, p in zip(trans_all['rho_ccr8'], trans_all['p_ccr8'])]
axes[0].scatter(trans_all['rho_ccr8'], range(len(trans_all)), c=colors_c, s=80, alpha=0.8, edgecolors='black')
axes[0].axvline(x=0, color='black', linewidth=0.5)
axes[0].set_yticks(range(len(trans_all)))
axes[0].set_yticklabels([f"{r['cohort']}:{r['sample']}" for _, r in trans_all.iterrows()], fontsize=7)
axes[0].set_xlabel('Spearman ρ (Transition score vs CCR8+)', fontsize=10)
axes[0].set_title('A. CCR8+ Spatial Gradient\n(Negative = higher in HS)', fontweight='bold')

# MKI67
colors_m = ['#b08ebf' if p < 0.05 and r > 0 else '#95A5A6' for r, p in zip(trans_all['rho_mki67'], trans_all['p_mki67'])]
axes[1].scatter(trans_all['rho_mki67'], range(len(trans_all)), c=colors_m, s=80, alpha=0.8, edgecolors='black')
axes[1].axvline(x=0, color='black', linewidth=0.5)
axes[1].set_yticks(range(len(trans_all)))
axes[1].set_yticklabels([f"{r['cohort']}:{r['sample']}" for _, r in trans_all.iterrows()], fontsize=7)
axes[1].set_xlabel('Spearman ρ (Transition score vs MKI67+)', fontsize=10)
axes[1].set_title('B. MKI67+ Spatial Gradient\n(Positive = higher in TC)', fontweight='bold')

plt.suptitle('Spatial Gradient: Tumor Core ← → Hard Stroma', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/fig4_transition_gradient.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{FIG_DIR}/fig4_transition_gradient.pdf', bbox_inches='tight')
plt.close()
print("  Saved: fig4_transition_gradient.png")

# Save tables
grade_df.to_csv(f'{OUT_DIR}/grade_trend.csv', index=False)
corr_df.to_csv(f'{OUT_DIR}/spot_correlation.csv', index=False)
me_all.to_csv(f'{OUT_DIR}/mutual_exclusion.csv', index=False)
trans_all.to_csv(f'{OUT_DIR}/transition_gradient.csv', index=False)

print(f"\nAll tables saved to {OUT_DIR}/")
print("Analysis complete!")
