"""
Global figure styling config for all figures.
统一颜色、字体、DPI，贯穿全文。

色盲友好配色权威标准
======================
本项目配色基于以下权威来源：

【权威来源 1】Okabe & Ito (2002) "Color Universal Design (CUD)"
   University of Tokyo / Jikei Medical School
   → https://jfly.uni-koeln.de/color/
   THE  foundational  work  on  colorblind-safe  palettes  for  scientific
   visualization.  Their  8-color  palette  is  the  gold  standard  used
   worldwide.

【权威来源 2】Bang Wong (2011) "Points of view: Color blindness"
   Nature Methods 8, 441  |  DOI: 10.1038/nmeth.1618
   Popularized Okabe & Ito's palette  in  the  biomedical  visualization
   community.  Key stat:  ~8% of men, ~0.5% of women have red-green
   colorblindness.  If your 3 reviewers  are  Northern  European  males,
   there's a 22% chance at least one is colorblind.

【权威来源 3】Nature Portfolio Final Artwork Guide
   "Avoid the use of red and green" — official journal policy.

【权威来源 4】Rougier et al. (2014) "Ten Simple Rules for Better Figures"
   PLOS Comput Biol 10(9): e1003833  |  Rule 6: Use Color Effectively
   "Avoid using red and green together. Use color to highlight, not decorate."

Okabe-Ito / Bang Wong 8-color safe palette:
┌──────────┬─────────┬──────────────────────────────┐
│ Black    │ #000000 │ text, axes, labels            │
│ Orange   │ #E69F00 │ category 1 (alternative to red)│
│ Sky Blue │ #56B4E9 │ category 2                    │
│ Bl.Green │ #009E73 │ category 3 (alternative to grn)│
│ Yellow   │ #F0E442 │ highlight only (not for area)  │
│ Blue     │ #0072B2 │ category 4                    │
│ Vermil.  │ #D55E00 │ category 5 (red alternative)   │
│ Rd.Purple│ #CC79A7 │ category 6                    │
└──────────┴─────────┴──────────────────────────────┘
Note: This project uses a modified Blue-Red scheme (not the full Okabe-Ito palette)
following the two-color contrast standard common in top journals (Nature, Cell, Science).
Blue (#2E5AAC) for the experimental/responder group, Red (#C0392B) for the
control/non-responder group — both are distinguishable by colorblind readers.
Core principle — "redundant coding": never rely SOLELY on color.
Always pair color with shapes, patterns, or text labels.

本项目配色映射（基于上述权威标准）：
  MPR (Responder)   = Blue   (#2E5AAC)  ~ Okabe Blue (darker, adjusted)
  Non-MPR           = Red    (#C0392B)  ~ Okabe Vermillion / Red alternative
  CCR8+ Treg        = Blue   (#2E5AAC)  ~ Okabe Blue (darker, adjusted)
  MKI67+ Treg       = Red    (#C0392B)  ~ Okabe Vermillion / Red alternative
  Highlight         = #D55E00           = Okabe Vermillion (exact)
  Neutral           = #555555           = gray
"""
import matplotlib as mpl

# ── 项目标准色（色盲友好，基于 Okabe-Ito / Bang Wong）──
COLOR_MPR = "#2E5AAC"          # 蓝 — Responder (MPR)
COLOR_NONMPR = "#C0392B"       # 红 — Non-Responder
COLOR_MPR_EDGE = "#1a3a7a"     # MPR 深蓝边
COLOR_NONMPR_EDGE = "#8b2522"  # Non-MPR 深红边
COLOR_CCR8 = "#2E5AAC"         # 蓝 — CCR8+ Treg
COLOR_MKI67 = "#C0392B"        # 红 — MKI67+ Treg
COLOR_HIGHLIGHT = "#D55E00"    # 强调色（Bang Wong vermillion）
COLOR_NEUTRAL = "#555555"      # 中性灰

# Bang Wong 原版 8 色（可选用）
BANG_WONG = [
    "#000000",  # Black
    "#E69F00",  # Orange
    "#56B4E9",  # Sky Blue
    "#009E73",  # Bluish Green
    "#F0E442",  # Yellow
    "#0072B2",  # Blue
    "#D55E00",  # Vermillion
    "#CC79A7",  # Reddish Purple
]

# ── 全局 rcParams ──
def apply_style():
    """Apply publication-level matplotlib styling."""
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        'axes.titlesize': 9,
        'axes.labelsize': 8,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'lines.linewidth': 0.8,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.pad_inches': 0.02,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })
