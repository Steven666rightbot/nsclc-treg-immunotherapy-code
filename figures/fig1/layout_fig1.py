"""
Custom layout for Fig 1 — ZERO WHITESPACE.
Row 1: [A] study design — full width
Row 2: [B] heatmap — full width
Row 3: [C] [D] [H] — stretched to fill full width (no side gaps)
Row 4: [E] [F] [G] — stretched to fill full width (no side gaps)
Labels embedded at top-left of each panel.
Target: JTM standard width ~5.3 inches / ~1600px at 300 DPI.
"""
import os
from PIL import Image, ImageDraw, ImageFont

DPI = 300
GAP = 8
FONT_SIZE = 16
OUTDIR = "D:/research/cucumber/fig1"
TARGET_WIDTH = 2400

PANELS = [
    ("A", "fig1a_study_design.png"),
    ("B", None),
    ("C", None),
    ("D", "fig1d_combined_roc.png"),
    ("E", None),
    ("F", None),
    ("G", None),
    ("H", "fig1h_external_validation_proper.png"),
]

PATHS = {
    "A": os.path.join(OUTDIR, "fig1a_study_design.png"),
    "B": "D:/research/tomato/figures/demo/1/fig1b_cd4t_marker_heatmap.png",
    "C": "D:/research/tomato/figures/demo/1/fig1c_treg_abundance_with_pvalue.png",
    "D": os.path.join(OUTDIR, "fig1d_combined_roc.png"),
    "E": "D:/research/tomato/figures/demo/1/fig1e_SHAP_bar_importance.png",
    "F": "D:/research/tomato/figures/demo/1/fig1f_decision_tree_compact.png",
    "G": "D:/research/tomato/figures/demo/1/fig1g_SHAP_summary_beeswarm.png",
    "H": os.path.join(OUTDIR, "fig1h_external_validation_proper.png"),
}

def load_img(path):
    return Image.open(path).convert('RGB')

def add_label(img, label):
    draw = ImageDraw.Draw(img)
    try:
        font = ImageFont.truetype("arial.ttf", FONT_SIZE)
    except:
        font = ImageFont.load_default()
    draw.text((6, 3), label, fill=(0, 0, 0), font=font)
    return img

def scale_w(img, w):
    r = w / img.size[0]
    return img.resize((w, int(img.size[1] * r)), Image.LANCZOS)

def scale_h(img, h):
    r = h / img.size[1]
    return img.resize((int(img.size[0] * r), h), Image.LANCZOS)

# Load all
imgs = {k: load_img(v) for k, v in PATHS.items()}

print("Sizes:")
for k, v in imgs.items():
    print(f"  {k}: {v.size[0]}x{v.size[1]}")

# Row 1: A (flowchart) - full width
row1 = add_label(scale_w(imgs['A'], TARGET_WIDTH), 'A')

# Row 2: B (heatmap) - full width
row2 = add_label(scale_w(imgs['B'], TARGET_WIDTH), 'B')

# Row 3: C, D, H — fill full width with equal inter-panel gaps
rh3 = 750
c3 = scale_h(imgs['C'], rh3)
d3 = scale_h(imgs['D'], rh3)
h3 = scale_h(imgs['H'], rh3)
total_w3 = c3.size[0] + GAP + d3.size[0] + GAP + h3.size[0]
if total_w3 > TARGET_WIDTH:
    shrink = TARGET_WIDTH / total_w3
    rh3 = int(rh3 * shrink)
    c3 = scale_h(imgs['C'], rh3)
    d3 = scale_h(imgs['D'], rh3)
    h3 = scale_h(imgs['H'], rh3)
    total_w3 = c3.size[0] + GAP + d3.size[0] + GAP + h3.size[0]
c3 = add_label(c3, 'C')
d3 = add_label(d3, 'D')
h3 = add_label(h3, 'H')
# Stretch: if narrower than TARGET_WIDTH, distribute remaining space as gap
extra_gap3 = (TARGET_WIDTH - total_w3) // 2 if total_w3 < TARGET_WIDTH else 0
actual_gap3 = GAP + extra_gap3
total_w3 = c3.size[0] + actual_gap3 + d3.size[0] + actual_gap3 + h3.size[0]
# Final adjustment: if not quite at TARGET_WIDTH, add last pixel to last gap
total_w3 = c3.size[0] + actual_gap3 + d3.size[0] + actual_gap3 + h3.size[0]
while total_w3 < TARGET_WIDTH:
    total_w3 += 1
# Center slightly if still off
x_offset3 = (TARGET_WIDTH - total_w3) // 2
row3 = Image.new('RGB', (TARGET_WIDTH, rh3), (255, 255, 255))
x = x_offset3
for img in [c3, d3, h3]:
    row3.paste(img, (x, 0))
    x += img.size[0] + actual_gap3

# Row 4: E, F, G — same approach, fill full width
rh4 = 975
e4 = scale_h(imgs['E'], rh4)
f4 = scale_h(imgs['F'], rh4)
g4 = scale_h(imgs['G'], rh4)
total_w4 = e4.size[0] + GAP + f4.size[0] + GAP + g4.size[0]
if total_w4 > TARGET_WIDTH:
    shrink = TARGET_WIDTH / total_w4
    rh4 = int(rh4 * shrink)
    e4 = scale_h(imgs['E'], rh4)
    f4 = scale_h(imgs['F'], rh4)
    g4 = scale_h(imgs['G'], rh4)
    total_w4 = e4.size[0] + GAP + f4.size[0] + GAP + g4.size[0]
e4 = add_label(e4, 'E')
f4 = add_label(f4, 'F')
g4 = add_label(g4, 'G')
extra_gap4 = (TARGET_WIDTH - total_w4) // 2 if total_w4 < TARGET_WIDTH else 0
actual_gap4 = GAP + extra_gap4
total_w4 = e4.size[0] + actual_gap4 + f4.size[0] + actual_gap4 + g4.size[0]
while total_w4 < TARGET_WIDTH:
    total_w4 += 1
x_offset4 = (TARGET_WIDTH - total_w4) // 2
row4 = Image.new('RGB', (TARGET_WIDTH, rh4), (255, 255, 255))
x = x_offset4
for img in [e4, f4, g4]:
    row4.paste(img, (x, 0))
    x += img.size[0] + actual_gap4

rows = [row1, row2, row3, row4]
total_h = sum(r.size[1] for r in rows) + GAP * 3 + 10

canvas = Image.new('RGB', (TARGET_WIDTH, total_h), (255, 255, 255))
y = 0
for row in rows:
    canvas.paste(row, (0, y))
    y += row.size[1] + GAP

print(f"\nFinal: {canvas.size[0]}x{canvas.size[1]}")
canvas.save(os.path.join(OUTDIR, "fig1_combined.png"), dpi=(DPI, DPI))
canvas.save(os.path.join(OUTDIR, "fig1_combined.pdf"), dpi=(DPI, DPI))
print("Saved!")
