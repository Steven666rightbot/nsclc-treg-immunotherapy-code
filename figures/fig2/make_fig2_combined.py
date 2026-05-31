"""
Figure 2 (Manuscript Figure 3 — CellChat) combined layout.
3-row × 2-column grid (six-grid layout):
  Row 1: A_left (MPR/Strong) | A_right (Non-MPR/Weak)
  Row 2: B (heatmap) | C (incoming signaling)
  Row 3: D (MHC bar) | E (MIF/CD99 bar)
Panel labels: 28px Arial Bold, matching Figure 2 (Slingshot/fig6) reference style.
"""
from PIL import Image, ImageDraw, ImageFont
from pathlib import Path
import os

OUT_DIR = Path("D:/research/cucumber/fig2_revised_merged")
DPI = 300
GAP = 8
JTM_W = 2400
LABEL_FONT_SIZE = 16

# ── Panel definitions ──
panels = {
    'A_left':  str(OUT_DIR / "fig2a_circle_merged_cellchat_style_strong.png"),
    'A_right': str(OUT_DIR / "fig2a_circle_merged_cellchat_style_weak.png"),
    'B':       str(OUT_DIR / "fig2b_differential_heatmap_treg_focused_12types.png"),
    'C':       str(OUT_DIR / "fig2c_treg_dual_subtype_incoming.png"),
    'D':       str(OUT_DIR / "fig2d_mhc_grouped_bar.png"),
    'E':       str(OUT_DIR / "fig2e_mif_cd99_grouped_bar.png"),
}

loaded = {}
for label, path in panels.items():
    if os.path.exists(path):
        loaded[label] = Image.open(path).convert("RGB")
    else:
        print(f"  ✗ MISSING: {path}")
        loaded[label] = Image.new("RGB", (400, 300), (245, 245, 245))

try:
    label_font = ImageFont.truetype("arial.ttf", LABEL_FONT_SIZE)
except Exception:
    label_font = ImageFont.load_default()


def layout_row(img1, img2, target_h):
    """
    Scale two images to the same height.
    If total width exceeds JTM_W, shrink proportionally.
    Returns (img1_r, img2_r, actual_h, x1, x2) where x1/x2 are centered positions.
    """
    ratio1 = img1.width / img1.height
    ratio2 = img2.width / img2.height
    
    w1 = target_h * ratio1
    w2 = target_h * ratio2
    total_w = w1 + GAP + w2
    
    if total_w > JTM_W:
        scale = JTM_W / total_w
        actual_h = int(target_h * scale)
        w1 = int(w1 * scale)
        w2 = int(w2 * scale)
    else:
        actual_h = target_h
        w1 = int(w1)
        w2 = int(w2)
    
    img1_r = img1.resize((w1, actual_h), Image.LANCZOS)
    img2_r = img2.resize((w2, actual_h), Image.LANCZOS)
    
    # Center horizontally
    x1 = (JTM_W - w1 - w2 - GAP) // 2
    x2 = x1 + w1 + GAP
    
    return img1_r, img2_r, actual_h, x1, x2


def add_label_with_bg(draw, label, x, y):
    """Draw panel label with white background patch."""
    bbox = draw.textbbox((0, 0), label, font=label_font)
    text_w = bbox[2] - bbox[0]
    text_h = bbox[3] - bbox[1]
    padding = 4
    draw.rectangle(
        [x, y, x + text_w + padding * 2, y + text_h + padding * 2],
        fill=(255, 255, 255)
    )
    draw.text((x + padding, y + padding), label, fill=(0, 0, 0), font=label_font)


# ── Row heights (tuned for balanced composition) ──
# 1800px wide / 300 DPI = 6 inches
# Width: 2400 px = 8 inches
# Row 1 (circle plots) needs ~1000px to keep labels readable
# Total height: 1000+8+650+8+550 = 2216 px
ROW_TARGETS = [1000, 650, 550]

# Row 1: A_left | A_right
a_left, a_right, h1, x1_1, x1_2 = layout_row(loaded['A_left'], loaded['A_right'], ROW_TARGETS[0])

# Row 2: B | C
b_img, c_img, h2, x2_1, x2_2 = layout_row(loaded['B'], loaded['C'], ROW_TARGETS[1])

# Row 3: D | E
d_img, e_img, h3, x3_1, x3_2 = layout_row(loaded['D'], loaded['E'], ROW_TARGETS[2])

# ── Composite ──
total_h = h1 + GAP + h2 + GAP + h3
canvas = Image.new("RGB", (JTM_W, total_h), (255, 255, 255))
draw = ImageDraw.Draw(canvas)

# Row 1
y = 0
canvas.paste(a_left, (x1_1, y))
canvas.paste(a_right, (x1_2, y))
add_label_with_bg(draw, "A", x1_1 + 6, y + 6)

# Row 2
y += h1 + GAP
canvas.paste(b_img, (x2_1, y))
canvas.paste(c_img, (x2_2, y))
add_label_with_bg(draw, "B", x2_1 + 6, y + 6)
add_label_with_bg(draw, "C", x2_2 + 6, y + 6)

# Row 3
y += h2 + GAP
canvas.paste(d_img, (x3_1, y))
canvas.paste(e_img, (x3_2, y))
add_label_with_bg(draw, "D", x3_1 + 6, y + 6)
add_label_with_bg(draw, "E", x3_2 + 6, y + 6)

# ── Save ──
png_path = OUT_DIR / "fig2_combined.png"
pdf_path = OUT_DIR / "fig2_combined.pdf"
canvas.save(str(png_path), dpi=(DPI, DPI))
canvas.save(str(pdf_path), dpi=(DPI, DPI))
print(f"OK {png_path.name}  ({canvas.size[0]}x{canvas.size[1]}, {os.path.getsize(png_path)//1024} KB)")
