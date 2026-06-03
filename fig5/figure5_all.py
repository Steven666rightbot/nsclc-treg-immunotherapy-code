#!/usr/bin/env python3
"""
Figure 5 — scTenifoldKnk virtual knockout validation
Generates all 4 panels (5A, 5B, 5C, 5D) as PNG + PDF.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['axes.facecolor'] = 'white'
matplotlib.rcParams['figure.facecolor'] = 'white'
matplotlib.rcParams['axes.grid'] = False
matplotlib.rcParams['grid.alpha'] = 0

DPI = 300
OUTDIR = "D:/Research/cucumber/fig5"

def save_both(fig, basename):
    """Save PNG and PDF."""
    fig.savefig(os.path.join(OUTDIR, f"{basename}.png"), dpi=DPI, bbox_inches='tight')
    fig.savefig(os.path.join(OUTDIR, f"{basename}.pdf"), dpi=DPI, bbox_inches='tight')

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cxcl13 = pd.read_csv("D:/Research/tomato/results/scTenifoldKnk_ko_ccr8_cxcl13/ko_CCR8_CXCL13.csv")
b2m = pd.read_csv("D:/Research/tomato/results/scTenifoldKnk_ko_mki67_b2m/ko_MKI67_B2M_lowres.csv")
go_bp = pd.read_csv("D:/Research/cucumber/fig5/go_bp_b2m_ko_top20.csv")

# Compute -log10(p.adj) safely
cxcl13['negLog10padj'] = -np.log10(cxcl13['p.adj'].replace(0, np.nan))
cxcl13['negLog10padj'] = cxcl13['negLog10padj'].fillna(-np.log10(1e-300))

b2m['negLog10padj'] = -np.log10(b2m['p.adj'].replace(0, np.nan))
b2m['negLog10padj'] = b2m['negLog10padj'].fillna(-np.log10(1e-300))

# Significance categories
def sig_color(row):
    if row['p.adj'] < 0.05 and abs(row['Z']) > 2:
        return '#E74C3C'      # BH significant red
    elif row['Z'] > 2:
        return '#7F8C8D'      # downstream dark gray
    else:
        return '#BDC3C7'      # not sig light gray

cxcl13['color'] = cxcl13.apply(sig_color, axis=1)
b2m['color'] = b2m.apply(sig_color, axis=1)

# ---------------------------------------------------------------------------
# Panel 5A — CXCL13 KO Volcano
# ---------------------------------------------------------------------------
def panel_5a():
    fig, ax = plt.subplots(figsize=(5, 5))

    ax.scatter(cxcl13['Z'], cxcl13['negLog10padj'],
               c=cxcl13['color'], s=15, alpha=0.8, edgecolors='none', rasterized=False)

    label_genes = ['CXCL13', 'CLDND2', 'AL136018.1', 'AC034238.2', 'CCL3']
    for idx, g in enumerate(label_genes):
        row = cxcl13[cxcl13['gene'] == g]
        if row.empty:
            continue
        # Alternate vertical offsets to avoid overlap
        offset = (6, 4) if idx % 2 == 0 else (6, -12)
        ax.annotate(g, (row['Z'].values[0], row['negLog10padj'].values[0]),
                    textcoords="offset points", xytext=offset,
                    fontsize=8, fontweight='bold', color='black')

    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.axvline(2, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.axvline(-2, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)

    # Threshold annotation
    ax.text(0.98, 0.92, '|Z| > 2, p.adj < 0.05',
            transform=ax.transAxes, fontsize=8, va='top', ha='right',
            color='#333333', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.8))

    ax.set_xlabel('Z-score', fontsize=10)
    ax.set_ylabel('-log$_{10}$(adjusted p-value)', fontsize=10)
    ax.set_title('CXCL13 KO in CCR8+ Treg', fontsize=11, fontweight='bold')

    n_sig = ((cxcl13['p.adj'] < 0.05) & (cxcl13['Z'].abs() > 2)).sum()
    ax.text(0.02, 0.98, f'{n_sig} significant genes (BH, p.adj < 0.05)',
            transform=ax.transAxes, fontsize=8, va='top', ha='left', style='italic')

    # Legend
    red_patch = mpatches.Patch(color='#E74C3C', label='p.adj < 0.05 & |Z| > 2')
    gray_patch = mpatches.Patch(color='#7F8C8D', label='Z > 2 only')
    ltgray_patch = mpatches.Patch(color='#BDC3C7', label='Not significant')
    ax.legend(handles=[red_patch, gray_patch, ltgray_patch], loc='lower right',
              fontsize=7, frameon=False)

    ax.set_xlim(-4, 8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    save_both(fig, 'fig5a_cxcl13_volcano')
    plt.close(fig)

# ---------------------------------------------------------------------------
# Panel 5B — B2M KO Volcano
# ---------------------------------------------------------------------------
def panel_5b():
    fig, ax = plt.subplots(figsize=(5, 5))

    ax.scatter(b2m['Z'], b2m['negLog10padj'],
               c=b2m['color'], s=15, alpha=0.8, edgecolors='none', rasterized=False)

    label_genes = ['B2M', 'PTMS', 'HAVCR1', 'HLA-DRB5', 'CD74']
    for idx, g in enumerate(label_genes):
        row = b2m[b2m['gene'] == g]
        if row.empty:
            continue
        x, y = row['Z'].values[0], row['negLog10padj'].values[0]
        weight = 'bold' if g in ('HLA-DRB5', 'CD74') else 'normal'
        # Alternate vertical offsets to avoid overlap
        offset = (6, 4) if idx % 2 == 0 else (6, -12)
        ax.annotate(g, (x, y), textcoords="offset points", xytext=offset,
                    fontsize=8, fontweight=weight, color='black')

    # Annotation near CD74
    cd74 = b2m[b2m['gene'] == 'CD74']
    if not cd74.empty:
        ax.annotate('CD74 (MIF receptor)',
                    (cd74['Z'].values[0], cd74['negLog10padj'].values[0]),
                    textcoords="offset points", xytext=(6, -24),
                    fontsize=7, color='#2980B9', fontweight='bold')

    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.axvline(2, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.axvline(-2, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)

    # Threshold annotation
    ax.text(0.98, 0.92, '|Z| > 2, p.adj < 0.05',
            transform=ax.transAxes, fontsize=8, va='top', ha='right',
            color='#333333', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.8))

    ax.set_xlabel('Z-score', fontsize=10)
    ax.set_ylabel('-log$_{10}$(adjusted p-value)', fontsize=10)
    ax.set_title('B2M KO in MKI67+ Treg', fontsize=11, fontweight='bold')

    n_sig = ((b2m['p.adj'] < 0.05) & (b2m['Z'].abs() > 2)).sum()
    ax.text(0.02, 0.98, f'{n_sig} significant genes (BH, p.adj < 0.05)',
            transform=ax.transAxes, fontsize=8, va='top', ha='left', style='italic')

    red_patch = mpatches.Patch(color='#E74C3C', label='p.adj < 0.05 & |Z| > 2')
    gray_patch = mpatches.Patch(color='#7F8C8D', label='Z > 2 only')
    ltgray_patch = mpatches.Patch(color='#BDC3C7', label='Not significant')
    ax.legend(handles=[red_patch, gray_patch, ltgray_patch], loc='lower right',
              fontsize=7, frameon=False)

    ax.set_xlim(-2, 8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    save_both(fig, 'fig5b_b2m_volcano')
    plt.close(fig)

# ---------------------------------------------------------------------------
# Panel 5C — Side-by-side horizontal barplot
# ---------------------------------------------------------------------------
def panel_5c():
    cxcl13_top = cxcl13.reindex(cxcl13['Z'].abs().sort_values(ascending=False).index).head(8).copy()
    b2m_top = b2m.reindex(b2m['Z'].abs().sort_values(ascending=False).index).head(10).copy()

    cxcl13_top['bar_color'] = cxcl13_top.apply(
        lambda r: '#E74C3C' if (r['p.adj'] < 0.05 and abs(r['Z']) > 2) else '#95A5A6', axis=1)
    b2m_top['bar_color'] = b2m_top.apply(
        lambda r: '#E74C3C' if (r['p.adj'] < 0.05 and abs(r['Z']) > 2) else '#95A5A6', axis=1)

    # Override HLA-DRB5 and CD74 to blue
    b2m_top.loc[b2m_top['gene'] == 'HLA-DRB5', 'bar_color'] = '#2980B9'
    b2m_top.loc[b2m_top['gene'] == 'CD74', 'bar_color'] = '#2980B9'

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [1, 1]})

    # Left — CXCL13
    ax = axes[0]
    y_pos = np.arange(len(cxcl13_top))
    ax.barh(y_pos, cxcl13_top['Z'], color=cxcl13_top['bar_color'], height=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(cxcl13_top['gene'], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Z-score', fontsize=10)
    ax.set_title('CXCL13 KO → CCR8+ Treg', fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axvline(0, color='black', linewidth=0.5)

    # Right — B2M
    ax = axes[1]
    y_pos = np.arange(len(b2m_top))
    ax.barh(y_pos, b2m_top['Z'], color=b2m_top['bar_color'], height=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(b2m_top['gene'], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Z-score', fontsize=10)
    ax.set_title('B2M KO → MKI67+ Treg', fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axvline(0, color='black', linewidth=0.5)

    # Legend
    red_patch = mpatches.Patch(color='#E74C3C', label='p.adj < 0.05 (BH)')
    gray_patch = mpatches.Patch(color='#95A5A6', label='Not significant')
    blue_patch = mpatches.Patch(color='#2980B9', label='CellChat link')
    axes[1].legend(handles=[red_patch, gray_patch, blue_patch], loc='lower right',
                   fontsize=8, frameon=False)

    fig.tight_layout()
    save_both(fig, 'fig5c_barplot_comparison')
    plt.close(fig)

# ---------------------------------------------------------------------------
# Panel 5D — Pathway Enrichment Bar Plot (GO BP, B2M KO)
# ---------------------------------------------------------------------------
def panel_5d():
    """
    Pathway enrichment bar plot using GO BP results from B2M KO.
    Y-axis: -log10(p.adjust)
    Colors: MHC-I/II = blue, Immune = red, Other = gray
    """
    df = go_bp.copy()
    df = df[df['p.adjust'] < 0.06].copy()
    df['negLog10padj'] = -np.log10(df['p.adjust'].replace(0, np.nan))
    df['negLog10padj'] = df['negLog10padj'].fillna(-np.log10(1e-300))

    # Truncate long descriptions
    df['short_desc'] = df['Description'].apply(lambda x: x[:55] + '...' if len(x) > 55 else x)

    def get_color(desc):
        desc_lower = desc.lower()
        if 'mhc class ii' in desc_lower or 'mhc-ii' in desc_lower or 'mhc class i' in desc_lower or 'mhc-i' in desc_lower:
            return '#2E5AAC'  # MHC-I/II blue
        elif 'mhc' in desc_lower and 'antigen' in desc_lower:
            return '#2E5AAC'  # MHC-related blue
        elif 'immune' in desc_lower or 'leukocyte' in desc_lower or 'lymphocyte' in desc_lower or 'inflammatory' in desc_lower:
            return '#E74C3C'  # Immune red
        else:
            return '#95A5A6'  # Other gray

    df['bar_color'] = df['Description'].apply(get_color)

    fig, ax = plt.subplots(figsize=(8, 5))

    x_pos = np.arange(len(df))
    bars = ax.bar(x_pos, df['negLog10padj'].values,
                  color=df['bar_color'].values, width=0.6,
                  edgecolor='white', linewidth=0.5, zorder=2)

    # Add value labels on top of bars
    for i, (bar, val) in enumerate(zip(bars, df['negLog10padj'].values)):
        ax.text(bar.get_x() + bar.get_width()/2., val + 0.05,
                f'{val:.2f}', ha='center', va='bottom', fontsize=7, color='#333333')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(df['short_desc'].values, rotation=45, ha='right', fontsize=7.5)
    ax.set_ylabel('-log$_{10}$(p.adjust)', fontsize=10)
    ax.set_title('B2M KO in MKI67+ Treg: GO Biological Process Enrichment', fontsize=11, fontweight='bold')

    # Significance threshold line
    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.text(len(df) - 0.5, -np.log10(0.05) + 0.05, 'p.adj = 0.05',
            ha='right', va='bottom', fontsize=7, color='gray', style='italic')

    # Legend
    legend_elements = [
        mpatches.Patch(color='#2E5AAC', label='MHC-I/II'),
        mpatches.Patch(color='#E74C3C', label='Immune'),
        mpatches.Patch(color='#95A5A6', label='Other'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', frameon=False, fontsize=8)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0, df['negLog10padj'].max() * 1.2)

    fig.tight_layout()
    save_both(fig, 'fig5d_pathway_enrichment')
    plt.close(fig)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    os.makedirs(OUTDIR, exist_ok=True)
    panel_5a()
    panel_5b()
    panel_5c()
    panel_5d()
    print("All panels saved to", OUTDIR)
