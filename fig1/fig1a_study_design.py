"""
Figure 1A — Study Design Flowchart (Redesigned)
Clean layout with uniform colorblind-safe palette.
Uses project standard colors from _global_config.py.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os
import sys
sys.path.insert(0, 'D:/research/cucumber')
from _global_config import COLOR_MPR, COLOR_NONMPR, COLOR_CCR8, COLOR_MKI67, COLOR_NEUTRAL

OUTDIR = "D:/research/cucumber/fig1"
os.makedirs(OUTDIR, exist_ok=True)

# ── Standard palette (from _global_config) ──
C_DATA   = COLOR_CCR8              # data input
C_METHOD = COLOR_MKI67             # analysis method
C_FIND   = COLOR_MPR               # key finding
C_VALID  = COLOR_NONMPR            # validation
C_DOWN   = COLOR_NEUTRAL   # #555555  gray   — downstream
C_BG     = '#F8F9FA'

fig, ax = plt.subplots(1, 1, figsize=(8.5, 5.2))
ax.set_xlim(0, 10)
ax.set_ylim(0, 6.2)
ax.axis('off')
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

def draw_box(ax, x, y, w, h, text, color, text_color='white', fontsize=9,
             subtext=None, subtext_size=7.5, subtext_color=None, alpha=0.92, radius=0.08):
    """Draw a rounded rectangle with text."""
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle=f"round,pad={radius}",
                         facecolor=color, edgecolor='white', linewidth=1.5, alpha=alpha,
                         mutation_scale=1)
    ax.add_patch(box)
    # Main text
    ax.text(x + w/2, y + h/2 + 0.12, text, ha='center', va='center',
            fontsize=fontsize, fontweight='bold', color=text_color,
            transform=ax.transData)
    if subtext:
        ax.text(x + w/2, y + h/2 - 0.22, subtext, ha='center', va='top',
                fontsize=subtext_size, color=subtext_color or text_color,
                transform=ax.transData, alpha=0.9, linespacing=1.3)

def draw_arrow(ax, x1, y1, x2, y2, color='#555555', lw=1.3):
    """Draw a clean arrow."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                connectionstyle='arc3,rad=0'))

# ============================================================
# Row 1: Data Input (blue)
# ============================================================
draw_box(ax, 0.3, 4.6, 4.2, 1.2,
         'GSE243013', color=C_DATA, fontsize=11,
         subtext='242 post-neoadjuvant chemo-immunotherapy\nNSCLC tumor samples')

draw_box(ax, 5.5, 4.6, 4.2, 1.2,
         'scRNA-seq Profiling', color=C_DATA, fontsize=11,
         subtext='1,749,316 immune cells\n51 annotated subtypes')

draw_arrow(ax, 4.5, 5.2, 5.5, 5.2)

# ============================================================
# Row 2: ML Method (orange)
# ============================================================
draw_box(ax, 1.5, 3.0, 7.0, 1.2,
         'XGBoost Classifier + SHAP Feature Importance',
         color=C_METHOD, fontsize=11,
         subtext='5-fold cross-validation | AUC = 0.805 ± 0.029 | 51 immune subtypes')

draw_arrow(ax, 5.0, 4.6, 5.0, 4.2)

# ============================================================
# Row 3: Key Finding (green)
# ============================================================
draw_box(ax, 2.0, 1.4, 6.0, 1.2,
         'Top 2 Predictive Features', color=C_FIND, fontsize=11,
         subtext='① CCR8⁺ Treg   (mean |SHAP| = 0.80)\n② MKI67⁺ Treg  (mean |SHAP| = 0.67)')

draw_arrow(ax, 5.0, 3.0, 5.0, 2.6)

# ============================================================
# Row 4: Downstream analyses + External validation
# ============================================================
# Four downstream boxes in a row
# Downstream analyses in a row — cleaner connector
downstream = [
    (0.3,  'CellChat\n(Intercellular)',   'Fig. 3'),
    (2.6,  'DoRothEA\n(TF Regulation)',    'Fig. 4'),
    (4.9,  'Visium\n(Spatial)',            'Fig. 5'),
    (7.2,  'Slingshot\n(Trajectory)',      'Fig. 2'),
]
for x, title, figref in downstream:
    draw_box(ax, x, 0.05, 1.8, 0.95,
             title, color=C_DOWN, fontsize=7.5,
             subtext=figref, subtext_size=7, subtext_color='#CCCCCC',
             alpha=0.85, radius=0.06)

# Horizontal connector ABOVE the downstream boxes
connector_y = 1.2
ax.plot([1.2, 8.1], [connector_y, connector_y], color=C_DOWN, lw=1.5, alpha=0.5, zorder=0)

# Main vertical arrow: Key Finding → connector
draw_arrow(ax, 5.0, 1.4, 5.0, connector_y, lw=1.5)

# Downward arrows: connector → each box
for x, _, _ in downstream:
    draw_arrow(ax, x + 0.9, connector_y, x + 0.9, 1.0, lw=1.0)

# ============================================================
# Side: External Validation (purple)
# ============================================================
draw_box(ax, 8.3, 1.55, 1.5, 0.9,
         'External\nValidation', color=C_VALID, fontsize=8,
         subtext='GSE207422\nAUC = 0.781', subtext_size=7)

# Arrow from Key Finding to External Validation (curved to avoid overlap)
ax.annotate('', xy=(8.3, 2.0), xytext=(8.0, 2.0),
            arrowprops=dict(arrowstyle='->', color=C_VALID, lw=1.2,
                            connectionstyle='arc3,rad=0.15'))

# Panel label (disabled — added by combined figure layout)

plt.savefig(os.path.join(OUTDIR, 'fig1a_study_design.png'), dpi=300,
            bbox_inches='tight', facecolor='white', pad_inches=0.05)
plt.savefig(os.path.join(OUTDIR, 'fig1a_study_design.pdf'), dpi=300,
            bbox_inches='tight', facecolor='white', pad_inches=0.05)
plt.close()

print("Saved: fig1a_study_design.png/pdf")
