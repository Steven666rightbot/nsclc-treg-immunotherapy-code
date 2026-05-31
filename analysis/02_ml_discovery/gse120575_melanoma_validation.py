#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE120575 Melanoma Cross-Cancer Validation (n=48 post-treatment biopsies)
Sade-Feldman et al., Cell 2018

Validates CCR8+ and MKI67+ Treg proportions as predictors of immunotherapy
response across a non-NSCLC malignancy.

Expected data format:
    - GSE120575_counts.csv (or .mtx) : gene x cell expression matrix
    - GSE120575_metadata.csv          : cell-level metadata with columns:
                                        'cell', 'patient', 'response'
                                        ('R' = responder, 'NR' = non-responder)

If data are not present, the script prints manual download instructions.
"""

import os
import sys
import gzip
import warnings
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
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

# Expected files
META_FILE = os.path.join(DATA_DIR, 'GSE120575_metadata.csv')
COUNTS_FILE = os.path.join(DATA_DIR, 'GSE120575_counts.csv.gz')
MTX_FILE = os.path.join(DATA_DIR, 'matrix.mtx')
BARCODES_FILE = os.path.join(DATA_DIR, 'barcodes.tsv')
GENES_FILE = os.path.join(DATA_DIR, 'genes.tsv')

# Marker genes
TREG_MARKERS = ['FOXP3', 'IL2RA', 'CTLA4']
CCR8_MARKER = 'CCR8'
MKI67_MARKER = 'MKI67'

# Response mapping
RESPONSE_MAP = {
    'R': 'Responder',
    'NR': 'Non-responder',
    'Responder': 'Responder',
    'Non-responder': 'Non-responder',
    ' responder': 'Responder',
    ' non-responder': 'Non-responder',
}

# ---------------------------------------------------------------------------
# 0. CHECK / DOWNLOAD DATA
# ---------------------------------------------------------------------------
def check_data():
    """Check if required data files exist; if not, print download instructions."""
    has_meta = os.path.exists(META_FILE)
    has_counts = os.path.exists(COUNTS_FILE) or os.path.exists(MTX_FILE)
    
    if has_meta and has_counts:
        return True
    
    print("=" * 70)
    print("GSE120575 data not found. Manual download required.")
    print("=" * 70)
    print("""
The GSE120575 dataset (Sade-Feldman et al., Cell 2018) is not included in
this repository due to size (~1.5 GB). Please download the files and place
them in: {data_dir}

Recommended download approach:
1. GEOquery (R):
   library(GEOquery)
   gse <- getGEO("GSE120575", GSEMatrix = TRUE)
   # Or use the Series Matrix file from GEO

2. Direct GEO download:
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575
   - Download "Series Matrix File(s)" for metadata
   - Download supplementary count matrix (may be .csv, .mtx, or .h5)

3. From the original publication (Sade-Feldman et al., Cell 2018):
   The single-cell count matrix and metadata are available as supplementary
   materials or via the Broad Institute Single Cell Portal.

Required files (choose one count format):
  A. Dense CSV:
     {meta}
     {counts}
  
  B. Sparse MTX (10x-like):
     {meta}
     {mtx}
     {barcodes}
     {genes}

Metadata CSV must contain at least these columns:
  - 'cell'     : cell barcode / identifier
  - 'patient'  : patient/biopsy identifier (48 unique patients)
  - 'response' : 'R' or 'NR' (or 'Responder' / 'Non-responder')
""".format(
        data_dir=DATA_DIR,
        meta=META_FILE,
        counts=COUNTS_FILE,
        mtx=MTX_FILE,
        barcodes=BARCODES_FILE,
        genes=GENES_FILE,
    ))
    return False


# ---------------------------------------------------------------------------
# 1. LOAD METADATA
# ---------------------------------------------------------------------------
def load_metadata():
    """Load and standardize metadata."""
    print("\n" + "=" * 70)
    print("Loading metadata...")
    print("=" * 70)
    
    meta = pd.read_csv(META_FILE)
    print(f"  Metadata shape: {meta.shape}")
    print(f"  Columns: {list(meta.columns)}")
    
    # Standardize column names (case-insensitive matching)
    col_map = {}
    for c in meta.columns:
        c_lower = c.lower().strip()
        if c_lower in ['cell', 'barcode', 'cell_id']:
            col_map[c] = 'cell'
        elif c_lower in ['patient', 'sample', 'patient_id', 'biopsy']:
            col_map[c] = 'patient'
        elif c_lower in ['response', 'response_status', 'responder', 'best_response']:
            col_map[c] = 'response'
    
    if 'cell' not in col_map.values() or 'patient' not in col_map.values() or 'response' not in col_map.values():
        print("  WARNING: Could not auto-map required columns. Please ensure metadata has:")
        print("    - 'cell'     : cell identifier")
        print("    - 'patient'  : patient/biopsy ID")
        print("    - 'response' : 'R'/'NR' or 'Responder'/'Non-responder'")
        print("  Available columns:", list(meta.columns))
        sys.exit(1)
    
    meta = meta.rename(columns={k: v for k, v in col_map.items()})
    
    # Standardize response labels
    meta['response_group'] = meta['response'].astype(str).str.strip().map(
        lambda x: RESPONSE_MAP.get(x, x)
    )
    
    n_patients = meta['patient'].nunique()
    n_cells = len(meta)
    print(f"  Patients: {n_patients}")
    print(f"  Cells: {n_cells}")
    print(f"  Response groups:\n{meta['response_group'].value_counts()}")
    
    return meta


# ---------------------------------------------------------------------------
# 2. LOAD EXPRESSION MATRIX
# ---------------------------------------------------------------------------
def load_expression():
    """Load count matrix in either dense CSV or sparse MTX format."""
    print("\n" + "=" * 70)
    print("Loading expression matrix...")
    print("=" * 70)
    
    if os.path.exists(COUNTS_FILE):
        print(f"  Reading dense CSV: {COUNTS_FILE}")
        counts = pd.read_csv(COUNTS_FILE, index_col=0)
        print(f"  Shape: {counts.shape}")
        return counts
    
    elif os.path.exists(MTX_FILE):
        print(f"  Reading sparse MTX: {MTX_FILE}")
        mtx = scipy.io.mmread(MTX_FILE)
        
        with open(BARCODES_FILE, 'r') as f:
            barcodes = [line.strip() for line in f]
        with open(GENES_FILE, 'r') as f:
            genes = [line.strip().split('\t')[0] for line in f]
        
        counts = pd.DataFrame(mtx.toarray(), index=genes, columns=barcodes)
        print(f"  Shape: {counts.shape}")
        return counts
    
    else:
        print("  ERROR: No expression matrix found.")
        sys.exit(1)


# ---------------------------------------------------------------------------
# 3. IDENTIFY TREG CELLS
# ---------------------------------------------------------------------------
def identify_tregs(counts, meta):
    """Identify Treg cells using canonical markers."""
    print("\n" + "=" * 70)
    print("Identifying Treg cells...")
    print("=" * 70)
    
    genes = counts.index.tolist()
    genes_upper = [g.upper() for g in genes]
    gene_map = {g.upper(): g for g in genes}
    
    # Find available Treg markers
    treg_avail = [gene_map[g] for g in TREG_MARKERS if g in gene_map]
    print(f"  Treg markers available: {treg_avail}")
    
    if len(treg_avail) < 2:
        print("  WARNING: Fewer than 2 Treg markers found. Trying alternative names...")
        alt_map = {'IL2RA': 'CD25', 'CTLA4': 'CD152'}
        for alt, orig in alt_map.items():
            if alt in gene_map and orig not in treg_avail:
                treg_avail.append(gene_map[alt])
        print(f"  After alternative mapping: {treg_avail}")
    
    # Subset to cells in metadata
    common_cells = list(set(counts.columns) & set(meta['cell']))
    counts_sub = counts[common_cells]
    meta_sub = meta[meta['cell'].isin(common_cells)].copy()
    meta_sub = meta_sub.set_index('cell')
    
    # Detect Treg: express at least 2 of 3 canonical markers
    treg_expr = (counts_sub.loc[treg_avail] > 0).sum(axis=0)
    is_treg = treg_expr >= 2
    
    treg_cells = is_treg[is_treg].index.tolist()
    print(f"  Total cells: {len(common_cells)}")
    print(f"  Treg cells (>=2 markers): {len(treg_cells)} ({100*len(treg_cells)/len(common_cells):.2f}%)")
    
    meta_sub['is_treg'] = is_treg
    
    return counts_sub, meta_sub, treg_cells, gene_map


# ---------------------------------------------------------------------------
# 4. COMPUTE CCR8+ AND MKI67+ TREG PROPORTIONS
# ---------------------------------------------------------------------------
def compute_proportions(counts_sub, meta_sub, treg_cells, gene_map):
    """Compute per-patient CCR8+ and MKI67+ Treg proportions."""
    print("\n" + "=" * 70)
    print("Computing CCR8+ and MKI67+ Treg proportions...")
    print("=" * 70)
    
    treg_counts = counts_sub[treg_cells]
    
    # CCR8+
    ccr8_gene = gene_map.get('CCR8')
    if ccr8_gene:
        ccr8_expr = treg_counts.loc[ccr8_gene]
        is_ccr8 = (ccr8_expr > 0).reindex(meta_sub.index, fill_value=False)
    else:
        print("  WARNING: CCR8 not found in gene list")
        is_ccr8 = pd.Series(False, index=meta_sub.index)
    
    # MKI67+
    mki67_gene = gene_map.get('MKI67')
    if mki67_gene:
        mki67_expr = treg_counts.loc[mki67_gene]
        is_mki67 = (mki67_expr > 0).reindex(meta_sub.index, fill_value=False)
    else:
        print("  WARNING: MKI67 not found in gene list")
        is_mki67 = pd.Series(False, index=meta_sub.index)
    
    meta_sub['is_ccr8_treg'] = is_ccr8 & meta_sub['is_treg']
    meta_sub['is_mki67_treg'] = is_mki67 & meta_sub['is_treg']
    
    # Per-patient aggregation
    patient_summary = []
    for patient, grp in meta_sub.groupby('patient'):
        total = len(grp)
        treg = grp['is_treg'].sum()
        ccr8 = grp['is_ccr8_treg'].sum()
        mki67 = grp['is_mki67_treg'].sum()
        
        response = grp['response_group'].iloc[0] if 'response_group' in grp.columns else 'Unknown'
        
        patient_summary.append({
            'patient': patient,
            'response': response,
            'total_cells': total,
            'treg_cells': treg,
            'ccr8_treg': ccr8,
            'mki67_treg': mki67,
            'treg_pct': 100 * treg / total if total > 0 else 0,
            'ccr8_pct_of_immune': 100 * ccr8 / total if total > 0 else 0,
            'mki67_pct_of_immune': 100 * mki67 / total if total > 0 else 0,
            'ccr8_pct_of_treg': 100 * ccr8 / treg if treg > 0 else 0,
            'mki67_pct_of_treg': 100 * mki67 / treg if treg > 0 else 0,
        })
    
    summary_df = pd.DataFrame(patient_summary)
    print(f"  Patients summarized: {len(summary_df)}")
    print(summary_df.head().to_string(index=False))
    
    return summary_df


# ---------------------------------------------------------------------------
# 5. STATISTICAL TESTS
# ---------------------------------------------------------------------------
def statistical_tests(summary_df):
    """Perform Mann-Whitney U tests for CCR8+ and MKI67+ Treg proportions."""
    print("\n" + "=" * 70)
    print("Statistical tests (Mann-Whitney U)...")
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
        
        if len(r_vals) > 0 and len(nr_vals) > 0:
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
            print(f"    Responder:     mean={np.mean(r_vals):.3f}, median={np.median(r_vals):.3f}, n={len(r_vals)}")
            print(f"    Non-responder: mean={np.mean(nr_vals):.3f}, median={np.median(nr_vals):.3f}, n={len(nr_vals)}")
            print(f"    Mann-Whitney U = {stat:.2f}, p = {pval:.4f}")
        else:
            print(f"\n  {label}: INSUFFICIENT DATA")
    
    results_df = pd.DataFrame(results).T
    results_df.to_csv(os.path.join(RESULT_DIR, 'gse120575_statistical_tests.csv'), index=True)
    print(f"\n  Saved: {RESULT_DIR}/gse120575_statistical_tests.csv")
    
    return results


# ---------------------------------------------------------------------------
# 6. PLOT SUPPLEMENTARY FIGURE 2
# ---------------------------------------------------------------------------
def plot_supplementary_figure_2(summary_df):
    """Generate Supplementary Figure 2: cross-cancer melanoma validation."""
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
    
    r_vals = summary_df[summary_df['response'] == 'Responder']['ccr8_pct_of_immune'].dropna().values
    nr_vals = summary_df[summary_df['response'] == 'Non-responder']['ccr8_pct_of_immune'].dropna().values
    if len(r_vals) > 0 and len(nr_vals) > 0:
        _, pval = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
        sig = f'p = {pval:.3f}' if pval >= 0.001 else f'p = {pval:.2e}'
        y_max = max(np.max(r_vals) if len(r_vals) > 0 else 0,
                    np.max(nr_vals) if len(nr_vals) > 0 else 0)
        ax.text(1.5, y_max * 0.95, sig, ha='center', va='bottom', fontsize=10, fontweight='bold')
    
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
    
    r_vals = summary_df[summary_df['response'] == 'Responder']['mki67_pct_of_immune'].dropna().values
    nr_vals = summary_df[summary_df['response'] == 'Non-responder']['mki67_pct_of_immune'].dropna().values
    if len(r_vals) > 0 and len(nr_vals) > 0:
        _, pval = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
        sig = f'p = {pval:.3f}' if pval >= 0.001 else f'p = {pval:.2e}'
        y_max = max(np.max(r_vals) if len(r_vals) > 0 else 0,
                    np.max(nr_vals) if len(nr_vals) > 0 else 0)
        ax.text(1.5, y_max * 0.95, sig, ha='center', va='bottom', fontsize=10, fontweight='bold')
        # Add significance asterisk if p < 0.05
        if pval < 0.05:
            ax.text(1.5, y_max * 1.05, '*', ha='center', va='bottom', fontsize=14, fontweight='bold', color='red')
    
    ax.set_ylabel('MKI67+ Treg (% of immune cells)', fontsize=11)
    ax.set_title('(B) MKI67+ Treg', fontweight='bold')
    ax.set_xlabel('Response', fontsize=11)
    
    fig.suptitle('Supplementary Figure 2 | Cross-cancer validation in melanoma\n'
                 '(GSE120575, n = 48 post-treatment biopsies)',
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
    """Save per-patient summary table."""
    summary_df.to_csv(os.path.join(RESULT_DIR, 'gse120575_per_patient_summary.csv'), index=False)
    print(f"\n  Saved: {RESULT_DIR}/gse120575_per_patient_summary.csv")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("GSE120575 Melanoma Cross-Cancer Validation")
    print("Sade-Feldman et al., Cell 2018 | n = 48 post-treatment biopsies")
    print("=" * 70)
    
    if not check_data():
        sys.exit(0)
    
    meta = load_metadata()
    counts = load_expression()
    counts_sub, meta_sub, treg_cells, gene_map = identify_tregs(counts, meta)
    summary_df = compute_proportions(counts_sub, meta_sub, treg_cells, gene_map)
    statistical_tests(summary_df)
    plot_supplementary_figure_2(summary_df)
    save_summary(summary_df)
    
    print("\n" + "=" * 70)
    print("GSE120575 VALIDATION COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
