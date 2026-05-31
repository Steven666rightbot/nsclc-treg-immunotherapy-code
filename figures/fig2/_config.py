"""Shared publication styling config for Figure 2 panels.
Import this at the top of each plotting script before any matplotlib calls.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

# ── Paths ──
RESULTS_DIR = Path("D:/research/tomato/results/cellchat_gse243013")
OUT_DIR = Path("D:/research/cucumber/fig2_revised_merged")
OUT_DIR.mkdir(exist_ok=True)

# ── Colors ──
# Colorblind-friendly: blue-red pair
COLOR_MPR = '#2E5AAC'       # blue (MPR)
COLOR_NONMPR = '#C0392B'    # red (Non-MPR)
COLOR_MPR_EDGE = '#1a3a7a'  # dark blue edge
COLOR_NONMPR_EDGE = '#8b2522'  # dark red edge
COLOR_MPR_ACCENT = '#2E5AAC'    # blue accent for MPR-high
COLOR_NONMPR_ACCENT = '#C0392B'  # red accent for Non-MPR-high
COLOR_BG = "#FAFAFA"

# ── Matplotlib rcParams (publication quality) ──
mpl.rcParams.update({
    # Font
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'legend.title_fontsize': 7.5,

    # Lines
    'lines.linewidth': 0.8,
    'patch.linewidth': 0.5,
    'axes.linewidth': 0.6,

    # Ticks
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.minor.size': 2,
    'ytick.minor.size': 2,

    # Figure
    'figure.dpi': 150,
    'savefig.dpi': 300,
    # bbox='tight' removed from defaults — each script handles layout explicitly
    'savefig.pad_inches': 0.05,

    # Grid
    'axes.grid': False,
})

# ── Figure size presets (Nature Comms: single col ~89mm, 1.5 col ~140mm, full ~178mm) ──
SINGLE_COL = (3.5, 3.0)    # inches
FULL_WIDTH = (7.0, 4.5)    # inches
WIDE_SHORT = (7.0, 2.2)    # for bar charts
WIDE_MEDIUM = (7.0, 3.5)   # for bubble plots


def save_fig(fig, name):
    """Save PNG + PDF with consistent naming."""
    png_path = OUT_DIR / f"{name}.png"
    pdf_path = OUT_DIR / f"{name}.pdf"
    fig.savefig(str(png_path), dpi=300, bbox_inches='tight')
    fig.savefig(str(pdf_path), bbox_inches='tight')
    print(f"  [OK] {name}.png  ({png_path.stat().st_size // 1024} KB)")
    print(f"  [OK] {name}.pdf  ({pdf_path.stat().st_size // 1024} KB)")
