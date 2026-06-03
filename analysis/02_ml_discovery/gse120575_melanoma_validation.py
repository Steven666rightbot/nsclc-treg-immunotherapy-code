#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE120575 Melanoma Cross-Cancer Validation
Sade-Feldman et al., Cell 2018

Validates CCR8+ and MKI67+ Treg proportions as predictors of immunotherapy
response across a non-NSCLC malignancy.

Dataset facts (verified from original publication & GEO):
    - 48 tumor samples (19 pre-treatment, 29 post-treatment)
    - From 32 melanoma patients (11 with longitudinal pre/post biopsies)
    - 16,291 CD45+ immune cells (Smart-seq2, TPM)
    - Response defined per lesion (Responder = CR/PR, Non-responder = SD/PD)

Analysis is performed at the lesion (sample) level to match the original
publication's lesion-level response classification.
"""

import os
import sys
import gzip
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
BASE_DIR = r'D:\Research\tomato'
DATA_DIR = os.path.join(BASE_DIR, 'data', 'raw', 'GSE120575')
RESULT_DIR = os.path.join(BASE_DIR, 'results', 'gse120575_validation')
FIG_DIR = os.path.join(BASE_DIR, 'figures', 'pub_figures')

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# Source files (downloaded from GEO)
TPM_FILE = os.path.join(DATA_DIR, 'GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz')
META_FILE = os.path.join(DATA_DIR, 'GSE120575_patient_ID_single_cells.txt.gz')

# Marker genes
TREG_MARKERS = ['FOXP3', 'IL2RA', 'CTLA4']
CCR8_MARKER = 'CCR8'
MKI67_MARKER = 'MKI67'

# ---------------------------------------------------------------------------
# 1. PARSE GEO METADATA
# ---------------------------------------------------------------------------
def parse_metadata():
    """Parse the GEO metadata file to get cell-level annotation."""
    print("\n" + "=" * 70)
    print("Parsing GEO metadata...")
    print("=" * 70)
    
    with gzip.open(META_FILE, 'rb') as f:
        lines = f.read().decode('latin-1').split('\r\n')
    
    header_idx = next(i for i, l in enumerate(lines) if l.startswith('Sample name'))
    header = lines[header_idx].split('\t')
    
    rows = []
    for line in lines[header_idx + 1:]:
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) >= 6 and parts[0].startswith('Sample'):
            rows.append(parts[:7])
    
    meta = pd.DataFrame(rows, columns=header[:7])
    
    # Keep only valid response entries
    meta = meta[meta['characteristics: response'].isin(['Responder', 'Non-responder'])]
    
    # Parse fields
    meta['cell_id'] = meta['title'].str.strip()
    meta['patient_raw'] = meta['characteristics: patinet ID (Pre=baseline; Post= on treatment)']
    meta['pre_post'] = meta['patient_raw'].apply(
        lambda x: 'Pre' if str(x).startswith('Pre') else ('Post' if str(x).startswith('Post') else 'Unknown')
    )
    meta['patient_base'] = meta['patient_raw'].str.replace(r'^(Pre|Post)_', '', regex=True)
    meta['response'] = meta['characteristics: response']
    meta['therapy'] = meta['characteristics: therapy']
    meta['is_enriched'] = meta['cell_id'].str.contains('_enriched', na=False)
    
    print(f"  Total cells: {len(meta)}")
    print(f"  Pre-treatment: {(meta['pre_post'] == 'Pre').sum()}")
    print(f"  Post-treatment: {(meta['pre_post'] == 'Post').sum()}")
    print(f"  Enriched cells: {meta['is_enriched'].sum()}")
    
    # Sample-level summary
    sample_summary = meta.groupby(['patient_raw', 'pre_post']).agg({
        'cell_id': 'count',
        'response': 'first',
        'patient_base': 'first',
    }).rename(columns={'cell_id': 'n_cells'}).reset_index()
    
    print(f"\n  Unique samples (lesions): {len(sample_summary)}")
    print(f"    Pre: {(sample_summary['pre_post'] == 'Pre').sum()}")
    print(f"    Post: {(sample_summary['pre_post'] == 'Post').sum()}")
    print(f"    Responders: {(sample_summary['response'] == 'Responder').sum()}")
    print(f"    Non-responders: {(sample_summary['response'] == 'Non-responder').sum()}")
    
    # Patients with both pre and post
    patient_prepost = sample_summary.pivot_table(
        index='patient_base', columns='pre_post', values='n_cells', aggfunc='sum'
    )
    n_both = patient_prepost.dropna().shape[0]
    print(f"  Patients with both Pre & Post: {n_both}")
    
    return meta


# ---------------------------------------------------------------------------
# 2. STREAM TPM MATRIX (only needed genes)
# ---------------------------------------------------------------------------
def load_expression(meta):
    """Stream the TPM matrix and extract only marker genes."""
    print("\n" + "=" * 70)
    print("Loading expression (streaming marker genes)...")
    print("=" * 70)
    
    needed_genes = TREG_MARKERS + [CCR8_MARKER, MKI67_MARKER]
    gene_expr = {}
    
    with gzip.open(TPM_FILE, 'rt', encoding='latin-1') as f:
        cell_ids = f.readline().strip().split('\t')
        patient_row = f.readline().strip().split('\t')
        
        for line in f:
            parts = line.strip().split('\t')
            gene = parts[0]
            if gene.upper() in [g.upper() for g in needed_genes]:
                gene_expr[gene.upper()] = np.array(parts[1:], dtype=np.float32)
                print(f"  Found {gene}")
                if len(gene_expr) == len(needed_genes):
                    break
    
    # Build cell-level DataFrame
    cell_df = pd.DataFrame({
        'cell_id': cell_ids,
        'patient_tpm': patient_row,
    })
    
    # Merge with metadata
    meta_cols = ['cell_id', 'patient_raw', 'pre_post', 'patient_base',
                 'response', 'therapy', 'is_enriched']
    cell_df = cell_df.merge(meta[meta_cols], on='cell_id', how='inner')
    
    # Add expression
    for gene in needed_genes:
        cell_df[gene] = gene_expr[gene.upper()]
    
    print(f"\n  Cells with metadata + expression: {len(cell_df)}")
    return cell_df


# ---------------------------------------------------------------------------
# 3. IDENTIFY TREG AND SUBTYPE CELLS
# ---------------------------------------------------------------------------
def classify_cells(cell_df):
    """Classify Treg and CCR8+/MKI67+ subsets."""
    print("\n" + "=" * 70)
    print("Classifying Treg cells...")
    print("=" * 70)
    
    # Treg: express at least 2 of 3 canonical markers (TPM > 0)
    treg_score = (cell_df['FOXP3'] > 0).astype(int) + \
                 (cell_df['IL2RA'] > 0).astype(int) + \
                 (cell_df['CTLA4'] > 0).astype(int)
    cell_df['is_treg'] = treg_score >= 2
    
    # Subtypes
    cell_df['is_ccr8_treg'] = cell_df['is_treg'] & (cell_df['CCR8'] > 0)
    cell_df['is_mki67_treg'] = cell_df['is_treg'] & (cell_df['MKI67'] > 0)
    
    n_treg = cell_df['is_treg'].sum()
    n_ccr8 = cell_df['is_ccr8_treg'].sum()
    n_mki67 = cell_df['is_mki67_treg'].sum()
    
    print(f"  Total cells: {len(cell_df)}")
    print(f"  Treg cells: {n_treg} ({100*n_treg/len(cell_df):.2f}%)")
    print(f"  CCR8+ Treg: {n_ccr8} ({100*n_ccr8/len(cell_df):.2f}% of immune)")
    print(f"  MKI67+ Treg: {n_mki67} ({100*n_mki67/len(cell_df):.2f}% of immune)")
    
    return cell_df


# ---------------------------------------------------------------------------
# 4. COMPUTE PER-SAMPLE PROPORTIONS
# ---------------------------------------------------------------------------
def compute_proportions(cell_df):
    """Compute per-sample (lesion-level) proportions."""
    print("\n" + "=" * 70)
    print("Computing per-sample proportions...")
    print("=" * 70)
    
    summary = []
    for sample_id, grp in cell_df.groupby('patient_raw'):
        total = len(grp)
        if total < 5:
            continue
        
        treg = grp['is_treg'].sum()
        ccr8 = grp['is_ccr8_treg'].sum()
        mki67 = grp['is_mki67_treg'].sum()
        
        response = grp['response'].iloc[0]
        pre_post = grp['pre_post'].iloc[0]
        patient_base = grp['patient_base'].iloc[0]
        
        summary.append({
            'sample_id': sample_id,
            'patient': patient_base,
            'pre_post': pre_post,
            'response': response,
            'total_cells': total,
            'treg_cells': treg,
            'ccr8_treg_cells': ccr8,
            'mki67_treg_cells': mki67,
            'ccr8_pct_of_immune': 100 * ccr8 / total,
            'mki67_pct_of_immune': 100 * mki67 / total,
            'ccr8_pct_of_treg': 100 * ccr8 / treg if treg > 0 else 0,
            'mki67_pct_of_treg': 100 * mki67 / treg if treg > 0 else 0,
        })
    
    summary_df = pd.DataFrame(summary)
    print(f"  Samples summarized: {len(summary_df)}")
    print(summary_df[['sample_id', 'response', 'ccr8_pct_of_immune', 'mki67_pct_of_immune']].head().to_string(index=False))
    
    return summary_df


# ---------------------------------------------------------------------------
# 5. STATISTICAL TESTS
# ---------------------------------------------------------------------------
def statistical_tests(summary_df):
    """Perform Mann-Whitney U tests."""
    print("\n" + "=" * 70)
    print("Statistical tests (Mann-Whitney U, two-sided)...")
    print("=" * 70)
    
    r = summary_df[summary_df['response'] == 'Responder']
    nr = summary_df[summary_df['response'] == 'Non-responder']
    
    results = {}
    
    for metric, label in [
        ('ccr8_pct_of_immune', 'CCR8+ Treg (% of immune cells)'),
        ('mki67_pct_of_immune', 'MKI67+ Treg (% of immune cells)'),
        ('ccr8_pct_of_treg', 'CCR8+ Treg (% of Tregs)'),
        ('mki67_pct_of_treg', 'MKI67+ Treg (% of Tregs)'),
    ]:
        r_vals = r[metric].dropna().values
        nr_vals = nr[metric].dropna().values
        
        stat, pval = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
        
        results[metric] = {
            'label': label,
            'responder_mean': np.mean(r_vals),
            'responder_median': np.median(r_vals),
            'responder_n': len(r_vals),
            'nonresponder_mean': np.mean(nr_vals),
            'nonresponder_median': np.median(nr_vals),
            'nonresponder_n': len(nr_vals),
            'mann_whitney_U': stat,
            'p_value': pval,
        }
        
        print(f"\n  {label}:")
        print(f"    Responder:     mean={np.mean(r_vals):.3f}%, median={np.median(r_vals):.3f}%, n={len(r_vals)}")
        print(f"    Non-responder: mean={np.mean(nr_vals):.3f}%, median={np.median(nr_vals):.3f}%, n={len(nr_vals)}")
        print(f"    Mann-Whitney U = {stat:.2f}, p = {pval:.4f}")
    
    results_df = pd.DataFrame(results).T
    results_df.to_csv(os.path.join(RESULT_DIR, 'gse120575_statistical_tests.csv'), index=True)
    print(f"\n  Saved: {RESULT_DIR}/gse120575_statistical_tests.csv")
    
    return results


# ---------------------------------------------------------------------------
# 6. PLOT SUPPLEMENTARY FIGURE 2
# ---------------------------------------------------------------------------
def plot_supplementary_figure_2(summary_df):
    """Generate Supplementary Figure 2."""
    print("\n" + "=" * 70)
    print("Generating Supplementary Figure 2...")
    print("=" * 70)
    
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 9,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.5))
    colors = {'Responder': '#3498db', 'Non-responder': '#e74c3c'}
    
    # Panel A: CCR8+ Treg
    ax = axes[0]
    plot_data = []
    for resp in ['Responder', 'Non-responder']:
        vals = summary_df[summary_df['response'] == resp]['ccr8_pct_of_immune'].dropna().values
        plot_data.append(vals)
    
    bp = ax.boxplot(plot_data, labels=['R', 'NR'], patch_artist=True,
                    widths=0.5,
                    boxprops=dict(edgecolor='black', linewidth=0.8),
                    medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(linewidth=0.8),
                    capprops=dict(linewidth=0.8))
    for patch, color in zip(bp['boxes'], [colors['Responder'], colors['Non-responder']]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Add jittered points
    for i, resp in enumerate(['Responder', 'Non-responder']):
        vals = summary_df[summary_df['response'] == resp]['ccr8_pct_of_immune'].dropna().values
        jitter = np.random.normal(i + 1, 0.04, size=len(vals))
        ax.scatter(jitter, vals, color='black', s=15, alpha=0.5, zorder=3)
    
    _, pval = stats.mannwhitneyu(plot_data[0], plot_data[1], alternative='two-sided')
    sig = f'p = {pval:.3f}' if pval >= 0.001 else f'p = {pval:.2e}'
    y_max = max(np.max(plot_data[0]) if len(plot_data[0]) > 0 else 0,
                np.max(plot_data[1]) if len(plot_data[1]) > 0 else 0)
    ax.text(1.5, y_max * 0.95, sig, ha='center', va='bottom', fontsize=10, fontweight='bold')
    if pval < 0.05:
        ax.text(1.5, y_max * 1.05, '*', ha='center', va='bottom', fontsize=14, fontweight='bold', color='red')
    
    ax.set_ylabel('CCR8+ Treg (% of immune cells)', fontsize=11)
    ax.set_title('(A) CCR8+ Treg', fontweight='bold')
    ax.set_xlabel('Response', fontsize=11)
    
    # Panel B: MKI67+ Treg
    ax = axes[1]
    plot_data = []
    for resp in ['Responder', 'Non-responder']:
        vals = summary_df[summary_df['response'] == resp]['mki67_pct_of_immune'].dropna().values
        plot_data.append(vals)
    
    bp = ax.boxplot(plot_data, labels=['R', 'NR'], patch_artist=True,
                    widths=0.5,
                    boxprops=dict(edgecolor='black', linewidth=0.8),
                    medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(linewidth=0.8),
                    capprops=dict(linewidth=0.8))
    for patch, color in zip(bp['boxes'], [colors['Responder'], colors['Non-responder']]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    for i, resp in enumerate(['Responder', 'Non-responder']):
        vals = summary_df[summary_df['response'] == resp]['mki67_pct_of_immune'].dropna().values
        jitter = np.random.normal(i + 1, 0.04, size=len(vals))
        ax.scatter(jitter, vals, color='black', s=15, alpha=0.5, zorder=3)
    
    _, pval = stats.mannwhitneyu(plot_data[0], plot_data[1], alternative='two-sided')
    sig = f'p = {pval:.3f}' if pval >= 0.001 else f'p = {pval:.2e}'
    y_max = max(np.max(plot_data[0]) if len(plot_data[0]) > 0 else 0,
                np.max(plot_data[1]) if len(plot_data[1]) > 0 else 0)
    ax.text(1.5, y_max * 0.95, sig, ha='center', va='bottom', fontsize=10, fontweight='bold')
    if pval < 0.05:
        ax.text(1.5, y_max * 1.05, '*', ha='center', va='bottom', fontsize=14, fontweight='bold', color='red')
    
    ax.set_ylabel('MKI67+ Treg (% of immune cells)', fontsize=11)
    ax.set_title('(B) MKI67+ Treg', fontweight='bold')
    ax.set_xlabel('Response', fontsize=11)
    
    fig.suptitle('Supplementary Figure 2 | Cross-cancer validation in melanoma\n'
                 '(GSE120575, n = 48 tumor samples)',
                 fontsize=12, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    
    fig.savefig(os.path.join(FIG_DIR, 'supplementary_figure_S2_melanoma_validation.png'),
                dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(FIG_DIR, 'supplementary_figure_S2_melanoma_validation.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"  Saved: {FIG_DIR}/supplementary_figure_S2_melanoma_validation.png")
    print(f"  Saved: {FIG_DIR}/supplementary_figure_S2_melanoma_validation.pdf")


# ---------------------------------------------------------------------------
# 7. SAVE SUMMARY TABLE
# ---------------------------------------------------------------------------
def save_summary(summary_df):
    """Save per-sample summary table."""
    summary_df.to_csv(os.path.join(RESULT_DIR, 'gse120575_per_sample_summary.csv'), index=False)
    print(f"\n  Saved: {RESULT_DIR}/gse120575_per_sample_summary.csv")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("GSE120575 Melanoma Cross-Cancer Validation")
    print("Sade-Feldman et al., Cell 2018")
    print("=" * 70)
    
    # Check data exists
    if not os.path.exists(TPM_FILE) or not os.path.exists(META_FILE):
        print("\nERROR: Required data files not found.")
        print(f"  Expected: {TPM_FILE}")
        print(f"  Expected: {META_FILE}")
        print("\nPlease download from GEO: GSE120575")
        sys.exit(1)
    
    meta = parse_metadata()
    cell_df = load_expression(meta)
    cell_df = classify_cells(cell_df)
    summary_df = compute_proportions(cell_df)
    statistical_tests(summary_df)
    plot_supplementary_figure_2(summary_df)
    save_summary(summary_df)
    
    print("\n" + "=" * 70)
    print("GSE120575 VALIDATION COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
