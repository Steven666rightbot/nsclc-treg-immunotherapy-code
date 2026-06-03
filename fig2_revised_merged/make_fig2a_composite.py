"""
Figure 2A — Composite: MPR (Strong) vs Non-MPR (Weak) CellChat circle plots.
Clean side-by-side composition without borders or panel labels.
"""
from PIL import Image
from pathlib import Path
import os

OUT_DIR = Path("D:/research/cucumber/fig2_revised_merged")
DPI = 300

# ── Load source images ──
strong_path = OUT_DIR / "fig2a_circle_merged_cellchat_style_strong.png"
weak_path = OUT_DIR / "fig2a_circle_merged_cellchat_style_weak.png"

strong = Image.open(strong_path).convert("RGB")
weak = Image.open(weak_path).convert("RGB")

# Scale to same height for alignment
target_h = 1000
if strong.size[1] != target_h:
    ratio = target_h / strong.size[1]
    strong = strong.resize((int(strong.size[0] * ratio), target_h), Image.LANCZOS)
if weak.size[1] != target_h:
    ratio = target_h / weak.size[1]
    weak = weak.resize((int(weak.size[0] * ratio), target_h), Image.LANCZOS)

# ── Layout ──
GAP = 40

panel_w = max(strong.size[0], weak.size[0])
panel_h = max(strong.size[1], weak.size[1])

canvas_w = panel_w * 2 + GAP
canvas_h = panel_h
canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))

# Paste images side by side, centered vertically
paste_y1 = (canvas_h - strong.size[1]) // 2
paste_y2 = (canvas_h - weak.size[1]) // 2
canvas.paste(strong, (0, paste_y1))
canvas.paste(weak, (panel_w + GAP, paste_y2))

out_path = OUT_DIR / "fig2a_combined_circles.png"
canvas.save(str(out_path), dpi=(DPI, DPI))
print(f"OK {out_path.name}  ({canvas.size[0]}x{canvas.size[1]}, {os.path.getsize(out_path)//1024} KB)")
