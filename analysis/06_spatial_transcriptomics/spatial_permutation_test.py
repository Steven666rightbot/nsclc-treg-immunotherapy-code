#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spatial Permutation Test for GSE189487
Test whether ECM-Contractile spatial correlation is significant
by permuting spatial coordinates while preserving score distributions.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# ===================== Paths =====================
DATA_PATH = r"D:\Research\tomato\results\gse189487_spatial\spot_level_scores.csv"
FIG_DIR = r"D:\Research\tomato\figures\spatial_permutation"
RESULT_DIR = r"D:\Research\tomato\results\spatial_permutation"

os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RESULT_DIR, exist_ok=True)

# ===================== 1. Load Data =====================
print("=" * 60)
print("Spatial Permutation Test: GSE189487")
print("=" * 60)

df = pd.read_csv(DATA_PATH)
print(f"Total spots: {len(df)}")
print(f"Samples: {df['sample'].unique()}")
print(f"Stages: {df['stage'].unique()}")

# ===================== 2. Observed Correlation =====================
print("\n" + "=" * 60)
print("Step 2: Observed Correlations")
print("=" * 60)

obs_all = pearsonr(df['ecm_score'], df['contractile_treg_score'])
print(f"Overall: r = {obs_all.statistic:.4f}, p = {obs_all.pvalue:.4e}")

obs_by_stage = {}
for stage in sorted(df['stage'].unique()):
    sub = df[df['stage'] == stage]
    r = pearsonr(sub['ecm_score'], sub['contractile_treg_score'])
    obs_by_stage[stage] = r.statistic
    print(f"  {stage}: r = {r.statistic:.4f}, p = {r.pvalue:.4e}, n = {len(sub)}")

# ===================== 3. Permutation Test =====================
print("\n" + "=" * 60)
print("Step 3: Permutation Test (n_perm = 1000)")
print("=" * 60)

def run_permutation(df_source, n_perm=1000, group_by='sample', seed=42):
    """
    Permute scores within each sample group to preserve sample-specific distributions,
    then compute correlation.
    """
    rng = np.random.RandomState(seed)
    perm_stats = []
    
    for i in range(n_perm):
        df_perm = df_source.copy()
        
        # Shuffle contractile score within each sample
        for grp in df_perm[group_by].unique():
            mask = df_perm[group_by] == grp
            idx = df_perm.loc[mask].index
            perm_vals = rng.permutation(df_perm.loc[mask, 'contractile_treg_score'].values)
            df_perm.loc[idx, 'contractile_treg_score_perm'] = perm_vals
        
        r = pearsonr(df_perm['ecm_score'], df_perm['contractile_treg_score_perm']).statistic
        perm_stats.append(r)
        
        if (i + 1) % 100 == 0:
            print(f"  Progress: {i+1}/{n_perm}")
    
    return np.array(perm_stats)

perm_all = run_permutation(df, n_perm=1000, group_by='sample', seed=42)

# Calculate p-value (two-tailed)
p_val_all = (np.sum(np.abs(perm_all) >= np.abs(obs_all.statistic)) + 1) / (len(perm_all) + 1)
print(f"\nPermutation test (overall):")
print(f"  Observed r: {obs_all.statistic:.4f}")
print(f"  Perm mean ± SD: {perm_all.mean():.4f} ± {perm_all.std():.4f}")
print(f"  95% CI of null: [{np.percentile(perm_all, 2.5):.4f}, {np.percentile(perm_all, 97.5):.4f}]")
print(f"  Permutation p-value: {p_val_all:.4e}")

# ===================== 4. Stage-specific Permutation =====================
print("\n" + "=" * 60)
print("Step 4: Stage-specific Permutation Tests")
print("=" * 60)

stage_perm_results = {}
for stage in sorted(df['stage'].unique()):
    sub = df[df['stage'] == stage].copy()
    obs_r = obs_by_stage[stage]
    
    # Within-stage permutation (shuffle within each sample in this stage)
    rng = np.random.RandomState(42)
    perm_stats = []
    for i in range(1000):
        df_perm = sub.copy()
        for grp in df_perm['sample'].unique():
            mask = df_perm['sample'] == grp
            idx = df_perm.loc[mask].index
            perm_vals = rng.permutation(df_perm.loc[mask, 'contractile_treg_score'].values)
            df_perm.loc[idx, 'contractile_treg_score_perm'] = perm_vals
        r = pearsonr(df_perm['ecm_score'], df_perm['contractile_treg_score_perm']).statistic
        perm_stats.append(r)
    
    perm_stats = np.array(perm_stats)
    p_val = (np.sum(np.abs(perm_stats) >= np.abs(obs_r)) + 1) / (len(perm_stats) + 1)
    
    stage_perm_results[stage] = {
        'observed_r': obs_r,
        'perm_mean': perm_stats.mean(),
        'perm_std': perm_stats.std(),
        'perm_p': p_val,
        'perm_95ci_low': np.percentile(perm_stats, 2.5),
        'perm_95ci_high': np.percentile(perm_stats, 97.5),
        'perm_distribution': perm_stats
    }
    
    print(f"  {stage}: obs_r={obs_r:.4f}, perm_p={p_val:.4e}, "
          f"null=[{np.percentile(perm_stats, 2.5):.4f}, {np.percentile(perm_stats, 97.5):.4f}]")

# ===================== 5. Visualizations =====================
print("\n" + "=" * 60)
print("Step 5: Generate Figures")
print("=" * 60)

# --- Figure 1: Overall Permutation Distribution ---
fig, ax = plt.subplots(figsize=(9, 6))
ax.hist(perm_all, bins=50, color='steelblue', alpha=0.7, edgecolor='white', density=True)
ax.axvline(obs_all.statistic, color='darkred', lw=2.5, linestyle='--',
           label=f'Observed r = {obs_all.statistic:.4f}\np = {p_val_all:.4e}')
ax.axvline(0, color='black', lw=1, linestyle='-')
ax.set_xlabel('Permuted Pearson r', fontsize=12)
ax.set_ylabel('Density', fontsize=12)
ax.set_title('Spatial Permutation Test: ECM vs Contractile Treg\n'
             '(Contractile scores permuted within samples, n=1000)', 
             fontsize=13, fontweight='bold')
ax.legend(fontsize=11)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '01_permutation_overall.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [Saved] 01_permutation_overall.png")

# --- Figure 2: Stage-specific Permutation Distributions ---
fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
stages = sorted(stage_perm_results.keys())

for ax, stage in zip(axes, stages):
    res = stage_perm_results[stage]
    perm_dist = res['perm_distribution']
    
    ax.hist(perm_dist, bins=40, color='steelblue', alpha=0.7, edgecolor='white', density=True)
    ax.axvline(res['observed_r'], color='darkred', lw=2.5, linestyle='--',
               label=f'Observed r = {res["observed_r"]:.4f}\np = {res["perm_p"]:.4e}')
    ax.axvline(0, color='black', lw=0.8, linestyle='-')
    ax.set_xlabel('Permuted Pearson r', fontsize=11)
    ax.set_title(f'{stage}\n(n = {len(df[df["stage"] == stage])} spots)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

axes[0].set_ylabel('Density', fontsize=12)
plt.suptitle('Stage-specific Spatial Permutation Tests', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '02_permutation_by_stage.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [Saved] 02_permutation_by_stage.png")

# --- Figure 3: Scatter with regression line (observed) ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Overall scatter
ax = axes[0]
ax.scatter(df['ecm_score'], df['contractile_treg_score'], 
           c=df['sample'].astype('category').cat.codes, cmap='tab10', 
           alpha=0.3, s=8, edgecolors='none')
z = np.polyfit(df['ecm_score'], df['contractile_treg_score'], 1)
p = np.poly1d(z)
ax.plot(df['ecm_score'].sort_values(), p(df['ecm_score'].sort_values()), 
        "r--", alpha=0.8, lw=2, label=f'r = {obs_all.statistic:.4f}, p = {p_val_all:.4e}')
ax.set_xlabel('ECM Score', fontsize=12)
ax.set_ylabel('Contractile Treg Score', fontsize=12)
ax.set_title('Overall Spatial Correlation', fontsize=13, fontweight='bold')
ax.legend(fontsize=11)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Stage boxplot
ax = axes[1]
stage_order = ['AIS', 'MIA', 'IAC']
df['ecm_contractile_product'] = df['ecm_score'] * df['contractile_treg_score']
df.boxplot(column='ecm_contractile_product', by='stage', ax=ax, 
           grid=False, patch_artist=True,
           boxprops=dict(facecolor='lightblue', alpha=0.7),
           medianprops=dict(color='darkred', lw=2))
ax.set_xlabel('Stage', fontsize=12)
ax.set_ylabel('ECM × Contractile Score', fontsize=12)
ax.set_title('ECM-Contractile Co-localization by Stage', fontsize=13, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.suptitle('')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '03_observed_correlation.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  [Saved] 03_observed_correlation.png")

# ===================== 6. Save Results =====================
print("\n" + "=" * 60)
print("Step 6: Save Results")
print("=" * 60)

# Overall result
overall_result = pd.DataFrame({
    'test': ['overall'],
    'n_spots': [len(df)],
    'observed_r': [obs_all.statistic],
    'parametric_p': [obs_all.pvalue],
    'perm_p': [p_val_all],
    'perm_mean': [perm_all.mean()],
    'perm_std': [perm_all.std()],
    'perm_95ci_low': [np.percentile(perm_all, 2.5)],
    'perm_95ci_high': [np.percentile(perm_all, 97.5)]
})

# Stage results
stage_results = []
for stage in stages:
    res = stage_perm_results[stage]
    stage_results.append({
        'test': stage,
        'n_spots': len(df[df['stage'] == stage]),
        'observed_r': res['observed_r'],
        'parametric_p': pearsonr(df[df['stage']==stage]['ecm_score'], 
                                  df[df['stage']==stage]['contractile_treg_score']).pvalue,
        'perm_p': res['perm_p'],
        'perm_mean': res['perm_mean'],
        'perm_std': res['perm_std'],
        'perm_95ci_low': res['perm_95ci_low'],
        'perm_95ci_high': res['perm_95ci_high']
    })

all_results = pd.concat([overall_result, pd.DataFrame(stage_results)], ignore_index=True)
all_results.to_csv(os.path.join(RESULT_DIR, 'permutation_test_results.csv'), index=False)
print(f"  [Saved] {RESULT_DIR}/permutation_test_results.csv")

# Save permutation distributions
np.save(os.path.join(RESULT_DIR, 'perm_distribution_overall.npy'), perm_all)
for stage in stages:
    np.save(os.path.join(RESULT_DIR, f'perm_distribution_{stage}.npy'), 
            stage_perm_results[stage]['perm_distribution'])
print("  [Saved] Permutation distributions (.npy)")

print("\n" + "=" * 60)
print("Analysis Complete!")
print(f"Figures: {FIG_DIR}")
print(f"Results: {RESULT_DIR}")
print("=" * 60)
