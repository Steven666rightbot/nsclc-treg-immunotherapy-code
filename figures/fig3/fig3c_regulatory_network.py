#!/usr/bin/env python3
"""
Figure 3C — Regulatory Network Diagram
========================================
Tri-axis model: Two ML-discovered Treg subtypes (CCR8+ and MKI67+)
are controlled by distinct transcriptional programs, leading to different
communication patterns with the tumor microenvironment.

Outputs: figures/fig3c_regulatory_network.{png,pdf} (300 DPI)
Style: Nature/Cell-level publication quality, vector-compatible
"""

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.path import Path
import numpy as np
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
import os

# ── Colour palette ──────────────────────────────────────────────────────────
CCR8_COLOR    = '#E74C3C'   # red
MKI67_COLOR   = '#C0392B'   # red
RESPONDER     = '#8bc98b'   # green
NON_RESPONDER = '#b08ebf'   # purple
LINE_GRAY     = '#555555'
LIGHT_CCR8    = '#fde0dd'   # light pink
LIGHT_MKI67   = '#deebf7'   # light blue
BG_COLOR      = 'white'
FONT_FAMILY   = 'Arial'

# ── Global style ────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': FONT_FAMILY,
    'font.size': 9,
    'axes.facecolor': BG_COLOR,
    'figure.facecolor': BG_COLOR,
    'savefig.facecolor': BG_COLOR,
    'savefig.dpi': 300,
    'pdf.fonttype': 42,       # editable text in PDF
    'ps.fonttype': 42,
})

# ── Figure setup ────────────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(10, 7))
ax.set_xlim(-0.5, 10.5)
ax.set_ylim(-1.0, 9.0)
ax.axis('off')

# ════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ════════════════════════════════════════════════════════════════════════════

def rounded_box(x, y, w, h, fill_color, edge_color, lw=2.0, label='',
                label_color='black', fontsize=10, fontweight='bold',
                alpha=0.85, rx=0.12, ry=0.12):
    """Draw a rounded rectangle with centred text."""
    box = FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle=f"round,pad=0,rounding_size={min(w,h)*rx}",
        facecolor=fill_color, edgecolor=edge_color,
        linewidth=lw, alpha=alpha, zorder=5,
    )
    ax.add_patch(box)
    ax.text(x, y, label, ha='center', va='center',
            fontsize=fontsize, fontweight=fontweight,
            color=label_color, zorder=6)


def circle_node(x, y, r, fill_color, edge_color, lw=2.0, label='',
                label_color='black', fontsize=10, fontweight='bold',
                alpha=0.85):
    """Draw a circular node with centred text."""
    circ = Circle((x, y), r, facecolor=fill_color, edgecolor=edge_color,
                  linewidth=lw, alpha=alpha, zorder=5)
    ax.add_patch(circ)
    ax.text(x, y, label, ha='center', va='center',
            fontsize=fontsize, fontweight=fontweight,
            color=label_color, zorder=6,
            fontdict={'fontfamily': FONT_FAMILY})


def draw_arrow(x1, y1, x2, y2, lw=1.5, color=LINE_GRAY, style='solid',
               alpha=0.75, connectionstyle='arc3,rad=0.0',
               arrowsize=12, head_length=0.15, head_width=0.18,
               zorder=3):
    """Draw a curved or straight arrow between two points."""
    style_dict = {'solid': '-', 'dashed': '--', 'dotted': ':'}
    linestyle = style_dict.get(style, '-')

    arrow = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle=mpatches.ArrowStyle.CurveFilledAB(
            head_length=head_length, head_width=head_width),
        linestyle=linestyle,
        linewidth=lw,
        color=color,
        alpha=alpha,
        connectionstyle=connectionstyle,
        zorder=zorder,
        mutation_scale=arrowsize,
    )
    ax.add_patch(arrow)


def draw_straight_arrow(x1, y1, x2, y2, lw=1.5, color=LINE_GRAY,
                        style='solid', alpha=0.75, arrowsize=12,
                        zorder=3):
    """Straight arrow with precise control."""
    dx = x2 - x1
    dy = y2 - y1
    arrow = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle=mpatches.ArrowStyle.CurveFilledAB(
            head_length=0.15, head_width=0.18),
        linestyle='-' if style == 'solid' else '--',
        linewidth=lw,
        color=color,
        alpha=alpha,
        connectionstyle='arc3,rad=0.0',
        zorder=zorder,
        mutation_scale=arrowsize,
    )
    ax.add_patch(arrow)


def annotation_text(x, y, text, color=LINE_GRAY, fontsize=8,
                    fontweight='normal', ha='center', va='center',
                    rotation=0, alpha=0.9, zorder=10):
    """Add a text annotation."""
    ax.text(x, y, text, ha=ha, va=va, fontsize=fontsize,
            fontweight=fontweight, color=color, rotation=rotation,
            alpha=alpha, zorder=zorder,
            fontdict={'fontfamily': FONT_FAMILY})


def draw_split_arrow(x_from, y_from, x_to_mid, y_split_top, y_split_bot,
                     color_top, color_bot, lw=1.5, style='solid', label=''):
    """Arrow that splits into two branches with different colours."""
    # Stem
    draw_arrow(x_from, y_from, x_to_mid, (y_split_top + y_split_bot)/2,
               lw=lw, color=LINE_GRAY, style=style)
    # Branches
    draw_arrow(x_to_mid, (y_split_top + y_split_bot)/2,
               x_to_mid + 0.4, y_split_top,
               lw=lw-0.3, color=color_top, style=style)
    draw_arrow(x_to_mid, (y_split_top + y_split_bot)/2,
               x_to_mid + 0.4, y_split_bot,
               lw=lw-0.3, color=color_bot, style=style)


# ════════════════════════════════════════════════════════════════════════════
# PANEL TITLE
# ════════════════════════════════════════════════════════════════════════════
ax.text(5.0, 8.6, 'Regulatory Network — Tri-axis Model',
        ha='center', va='center', fontsize=13, fontweight='bold',
        color='#222222', zorder=10)

# ════════════════════════════════════════════════════════════════════════════
# 1  —  LEFT COLUMN: Treg Subtypes
# ════════════════════════════════════════════════════════════════════════════
box_w = 1.8
box_h = 0.70

# CCR8+ Treg
rounded_box(1.0, 6.3, box_w, box_h,
            fill_color=LIGHT_CCR8, edge_color=CCR8_COLOR, lw=2.5,
            label='CCR8+ Treg', label_color=CCR8_COLOR,
            fontsize=11, fontweight='bold')
# MKI67+ Treg
rounded_box(1.0, 4.7, box_w, box_h,
            fill_color=LIGHT_MKI67, edge_color=MKI67_COLOR, lw=2.5,
            label='MKI67+ Treg', label_color=MKI67_COLOR,
            fontsize=11, fontweight='bold')
# Label above left column
annotation_text(1.0, 7.4, 'Treg Subtypes', fontsize=9, fontweight='bold',
                color='#333333')

# ════════════════════════════════════════════════════════════════════════════
# 2  —  CENTER COLUMN: Transcription Factor Hubs
# ════════════════════════════════════════════════════════════════════════════
tf_x = 4.2
r = 0.45

# HSF1 — prominent (larger radius)
big_r = 0.55
circle_node(tf_x, 7.0, big_r, fill_color='#fff3cd', edge_color='#e67e22',
            lw=2.5, label='HSF1', label_color='#d35400',
            fontsize=10, fontweight='bold')
annotation_text(tf_x, 7.7, 'Stress /', fontsize=6.5, color='#d35400')
annotation_text(tf_x, 7.3, 'inflammation', fontsize=6.5, color='#d35400')

# BCL6
circle_node(tf_x, 5.5, r, fill_color='#e8daef', edge_color='#8e44ad',
            lw=2.0, label='BCL6', label_color='#6c3483',
            fontsize=10, fontweight='bold')
annotation_text(tf_x + 0.8, 5.5, 'GC Tfh-like',
                fontsize=6.5, color='#6c3483', ha='left')

# RFX5
circle_node(tf_x, 4.0, r, fill_color='#d5f5e3', edge_color='#27ae60',
            lw=2.0, label='RFX5', label_color='#1e8449',
            fontsize=10, fontweight='bold')
annotation_text(tf_x + 0.75, 4.0, 'MHC-II Ag\npresentation',
                fontsize=6.0, color='#1e8449', ha='left')

# MAX/MYC
circle_node(tf_x, 2.7, r, fill_color='#e8f8f5', edge_color='#1abc9c',
            lw=2.0, label='MAX/MYC', label_color='#148f77',
            fontsize=9, fontweight='bold')
annotation_text(tf_x + 0.75, 2.7, 'Proliferation',
                fontsize=6.5, color='#148f77', ha='left')

# Label above center column
annotation_text(tf_x, 8.2, 'Transcription Factor Hubs',
                fontsize=9, fontweight='bold', color='#333333')

# ════════════════════════════════════════════════════════════════════════════
# 3  —  RIGHT COLUMN: Effector Pathways
# ════════════════════════════════════════════════════════════════════════════
rx = 7.8
rw = 2.0
rh = 0.60

# MHC-I Antigen Presentation
rounded_box(rx, 7.0, rw, rh,
            fill_color='#f5f5f5', edge_color=LINE_GRAY, lw=1.5,
            label='MHC-I (HLA-A/B/C)', label_color='#333333',
            fontsize=8.5, fontweight='bold')
annotation_text(rx, 6.5, 'Antigen Presentation',
                fontsize=7, color=LINE_GRAY)

# MHC-II Antigen Presentation
rounded_box(rx, 5.5, rw, rh,
            fill_color='#f5f5f5', edge_color=LINE_GRAY, lw=1.5,
            label='MHC-II (HLA-DR/DQ/DM)', label_color='#333333',
            fontsize=8.5, fontweight='bold')
annotation_text(rx, 5.0, 'Antigen Presentation',
                fontsize=7, color=LINE_GRAY)

# MIF/CD99 Inflammatory Signaling
rounded_box(rx, 4.0, rw, rh,
            fill_color='#fef9e7', edge_color='#f39c12', lw=1.5,
            label='MIF / CD99', label_color='#d68910',
            fontsize=9, fontweight='bold')
annotation_text(rx, 3.5, 'Inflammatory Signaling',
                fontsize=7, color='#d68910')

# Heat Shock / Stress Response
rounded_box(rx, 2.5, rw, rh,
            fill_color='#fdedec', edge_color='#e74c3c', lw=1.5,
            label='HSP genes', label_color='#c0392b',
            fontsize=9, fontweight='bold')
annotation_text(rx, 2.0, 'Heat Shock / Stress',
                fontsize=7, color='#c0392b')

# Label above right column
annotation_text(rx, 8.2, 'Effector Pathways',
                fontsize=9, fontweight='bold', color='#333333')

# ════════════════════════════════════════════════════════════════════════════
# 4  —  BOTTOM: Functional Consequences
# ════════════════════════════════════════════════════════════════════════════
# Box spanning the full width at the bottom
conseq_y = 0.6
conseq_w = 8.5
conseq_h = 0.7
box2 = FancyBboxPatch(
    (5.0 - conseq_w/2, conseq_y - conseq_h/2), conseq_w, conseq_h,
    boxstyle="round,pad=0,rounding_size=0.2",
    facecolor='#fafafa', edgecolor='#aaaaaa',
    linewidth=1.2, alpha=0.9, zorder=4,
)
ax.add_patch(box2)
ax.text(5.0, conseq_y,
        'Immunosuppressive TME   ←   Treg plasticity   →   Anti-tumour immunity',
        ha='center', va='center',
        fontsize=10, fontweight='bold',
        color='#333333', zorder=6)
annotation_text(5.0, conseq_y - 0.55,
                'Functional Consequences',
                fontsize=9, fontweight='bold', color='#555555')

# ── Column-consequence connector lines ──
ax.plot([1.0, 1.0], [conseq_y + conseq_h/2 + 0.1, 5.7],
        color=LINE_GRAY, lw=0.8, ls='--', alpha=0.4, zorder=2)
ax.plot([tf_x, tf_x], [conseq_y + conseq_h/2 + 0.1, 2.0],
        color=LINE_GRAY, lw=0.8, ls='--', alpha=0.4, zorder=2)
ax.plot([rx, rx], [conseq_y + conseq_h/2 + 0.1, 1.7],
        color=LINE_GRAY, lw=0.8, ls='--', alpha=0.4, zorder=2)


# ════════════════════════════════════════════════════════════════════════════
# 5  —  CONNECTIONS: Arrows between nodes
# ════════════════════════════════════════════════════════════════════════════

# ── CCR8+ Treg → TF hubs ───────────────────────────────────────────────
# CCR8+ (1.0, 6.3) → RFX5 (4.2, 4.0) — strong, solid
draw_arrow(1.9, 6.15, 3.65, 4.1,
           lw=2.8, color=CCR8_COLOR, style='solid', alpha=0.7,
           connectionstyle='arc3,rad=-0.15', arrowsize=14)

# CCR8+ (1.0, 6.3) → MAX/MYC (4.2, 2.7) — strong, solid
draw_arrow(1.9, 6.0, 3.65, 2.85,
           lw=2.5, color=CCR8_COLOR, style='solid', alpha=0.65,
           connectionstyle='arc3,rad=-0.25', arrowsize=14)

# CCR8+ (1.0, 6.3) → BCL6 (4.2, 5.5) — moderate, dashed (indirect)
draw_arrow(1.9, 6.35, 3.65, 5.5,
           lw=1.8, color=CCR8_COLOR, style='dashed', alpha=0.6,
           connectionstyle='arc3,rad=0.1', arrowsize=11)

# ── MKI67+ Treg → TF hubs ─────────────────────────────────────────────
# MKI67+ (1.0, 4.7) → HSF1 (4.2, 7.0) — strong, solid
draw_arrow(1.9, 4.9, 3.65, 6.85,
           lw=2.8, color=MKI67_COLOR, style='solid', alpha=0.7,
           connectionstyle='arc3,rad=0.25', arrowsize=14)

# MKI67+ (1.0, 4.7) → BCL6 (4.2, 5.5) — moderate, dashed
draw_arrow(1.9, 4.75, 3.65, 5.5,
           lw=1.8, color=MKI67_COLOR, style='dashed', alpha=0.6,
           connectionstyle='arc3,rad=0.15', arrowsize=11)

# ── HSF1 → Effector pathways ───────────────────────────────────────────
# HSF1 (4.2, 7.0) → HSP genes / Stress (7.8, 2.5) — strong
draw_arrow(4.75, 6.8, 6.7, 2.55,
           lw=2.5, color='#e67e22', style='solid', alpha=0.65,
           connectionstyle='arc3,rad=-0.3', arrowsize=13)
# p-value annotation
annotation_text(5.6, 4.0, 'p = 0.003',
                fontsize=7, color='#e67e22', fontweight='bold',
                rotation=-28, alpha=0.75)

# HSF1 → MIF/CD99 (7.8, 4.0) — indirect
draw_arrow(4.75, 6.9, 6.7, 4.1,
           lw=1.5, color='#e67e22', style='dashed', alpha=0.5,
           connectionstyle='arc3,rad=-0.2', arrowsize=10)

# ── BCL6 → Effector pathways ───────────────────────────────────────────
# BCL6 (4.2, 5.5) → MHC-II (7.8, 5.5)
draw_arrow(4.75, 5.5, 6.7, 5.5,
           lw=2.2, color='#8e44ad', style='solid', alpha=0.65,
           connectionstyle='arc3,rad=0.0', arrowsize=13)
annotation_text(5.7, 5.2, 'p = 0.01',
                fontsize=7, color='#8e44ad', fontweight='bold',
                alpha=0.75)

# BCL6 → MIF/CD99 (7.8, 4.0)
draw_arrow(4.75, 5.3, 6.7, 4.1,
           lw=1.5, color='#8e44ad', style='dashed', alpha=0.5,
           connectionstyle='arc3,rad=-0.15', arrowsize=10)

# ── RFX5 → Effector pathways ───────────────────────────────────────────
# RFX5 (4.2, 4.0) → MHC-II (7.8, 5.5) — strong
draw_arrow(4.75, 4.2, 6.7, 5.35,
           lw=2.8, color='#27ae60', style='solid', alpha=0.7,
           connectionstyle='arc3,rad=0.2', arrowsize=14)
annotation_text(5.7, 5.8, 'p = 0.001',
                fontsize=7, color='#27ae60', fontweight='bold',
                alpha=0.75)
# RFX5 → MHC-I (7.8, 7.0)
draw_arrow(4.75, 4.1, 6.7, 6.85,
           lw=2.0, color='#27ae60', style='solid', alpha=0.6,
           connectionstyle='arc3,rad=0.35', arrowsize=12)
annotation_text(5.6, 5.95, 'p = 0.008',
                fontsize=7, color='#27ae60', fontweight='bold',
                rotation=25, alpha=0.75)

# ── MAX/MYC → Effector pathways ────────────────────────────────────────
# MAX/MYC (4.2, 2.7) → MHC-I (7.8, 7.0)
draw_arrow(4.75, 3.0, 6.7, 6.8,
           lw=2.0, color='#1abc9c', style='solid', alpha=0.6,
           connectionstyle='arc3,rad=0.5', arrowsize=12)

# MAX/MYC → MIF/CD99 (7.8, 4.0)
draw_arrow(4.75, 2.8, 6.7, 4.0,
           lw=1.5, color='#1abc9c', style='dashed', alpha=0.5,
           connectionstyle='arc3,rad=0.2', arrowsize=10)


# ════════════════════════════════════════════════════════════════════════════
# 6  —  ANNOTATIONS: NR↑ / R↓ labels
# ════════════════════════════════════════════════════════════════════════════

# On key edges
# CCR8+ → MHC-II via RFX5 → NR↑
annotation_text(5.9, 6.1, 'NR ↑', fontsize=8, fontweight='bold',
                color=NON_RESPONDER, alpha=0.85)
# MKI67+ → Stress via HSF1 → R↓
annotation_text(5.8, 3.7, 'R ↓', fontsize=8, fontweight='bold',
                color=RESPONDER, alpha=0.85)
# BCL6 → MHC-II → both NR↑
annotation_text(6.5, 5.15, 'NR ↑', fontsize=8, fontweight='bold',
                color=NON_RESPONDER, alpha=0.85)

# ════════════════════════════════════════════════════════════════════════════
# 7  —  LEGEND
# ════════════════════════════════════════════════════════════════════════════
legend_x = 8.8
legend_y = 1.3

ax.text(legend_x, legend_y + 0.6, 'Legend',
        fontsize=8, fontweight='bold', color='#333333',
        ha='center', zorder=10)

# Solid line
ax.plot([legend_x - 0.5, legend_x + 0.5], [legend_y + 0.3, legend_y + 0.3],
        color=LINE_GRAY, lw=2.0, zorder=10)
annotation_text(legend_x + 0.8, legend_y + 0.3, 'Direct regulation',
                fontsize=7, ha='left', color='#444444')

# Dashed line
ax.plot([legend_x - 0.5, legend_x + 0.5], [legend_y, legend_y],
        color=LINE_GRAY, lw=1.5, ls='--', zorder=10)
annotation_text(legend_x + 0.8, legend_y, 'Indirect / inferred',
                fontsize=7, ha='left', color='#444444')

# Thick arrow
ax.plot([legend_x - 0.5, legend_x + 0.5], [legend_y - 0.3, legend_y - 0.3],
        color=LINE_GRAY, lw=3.0, zorder=10)
annotation_text(legend_x + 0.8, legend_y - 0.3, 'Strong validation',
                fontsize=7, ha='left', color='#444444')

# Responder / Non-responder
circle_node(legend_x, legend_y - 0.8, 0.12,
            fill_color=RESPONDER, edge_color=RESPONDER,
            lw=0, label='', alpha=0.9)
annotation_text(legend_x + 0.5, legend_y - 0.8, 'Responder (R) ↓',
                fontsize=7, ha='left', color=RESPONDER, fontweight='bold')

circle_node(legend_x, legend_y - 1.1, 0.12,
            fill_color=NON_RESPONDER, edge_color=NON_RESPONDER,
            lw=0, label='', alpha=0.9)
annotation_text(legend_x + 0.5, legend_y - 1.1, 'Non-responder (NR) ↑',
                fontsize=7, ha='left', color=NON_RESPONDER, fontweight='bold')


# ════════════════════════════════════════════════════════════════════════════
# 8  —  SUB-COLUMN HEADERS & DIVIDERS
# ════════════════════════════════════════════════════════════════════════════

# Light vertical dividers between columns
xx = [2.6, 6.0]
for xpos in xx:
    ax.axvline(x=xpos, ymin=0.05, ymax=0.88, color='#dddddd',
               lw=0.8, ls=':', zorder=1)


# ════════════════════════════════════════════════════════════════════════════
# SAVE
# ════════════════════════════════════════════════════════════════════════════
OUT_DIR = 'D:/research/cucumber/fig3'
os.makedirs(OUT_DIR, exist_ok=True)

for ext in ['png', 'pdf']:
    fpath = os.path.join(OUT_DIR, f'fig3c_regulatory_network.{ext}')
    fig.savefig(fpath, dpi=300, bbox_inches='tight',
                facecolor=BG_COLOR, edgecolor='none',
                transparent=False)
    print(f'Saved {fpath}')

plt.close(fig)
print('✓ Figure 3C — Regulatory Network Diagram complete.')
