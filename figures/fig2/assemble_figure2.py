import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
from pathlib import Path

OUT_DIR = Path("D:/research/cucumber/fig2_revised_merged")
IMG_DIR = OUT_DIR

# ======================== Load & Prepare Images ========================

# 2A: Use the R/circlize CellChat composite (fancier than Python version)
img_2a = Image.open(IMG_DIR / "fig2a_combined_circles.png")

# Other panels
img_2b = Image.open(IMG_DIR / "fig2b_differential_heatmap_treg_focused_12types.png")
img_2c = Image.open(IMG_DIR / "fig2c_treg_dual_subtype_incoming.png")
img_2d = Image.open(IMG_DIR / "fig2d_mhc_grouped_bar.png")
img_2e = Image.open(IMG_DIR / "fig2e_mif_cd99_grouped_bar.png")

# Aspect ratios (width / height)
aspects = {
    'A': img_2a.width / img_2a.height,
    'B': img_2b.width / img_2b.height,
    'C': img_2c.width / img_2c.height,
    'D': img_2d.width / img_2d.height,
    'E': img_2e.width / img_2e.height,
}

# ======================== Compute Layout ========================
fig_w = 6.3
margin = 0.4
avail_w = fig_w - 2 * margin
gap = 0.22

# Row 0: 2A spans full width
disp_w_a = avail_w
disp_h_a = disp_w_a / aspects['A']

# Row 1: optimize column widths so 2B and 2C have similar display heights
w_b = avail_w * aspects['B'] / (aspects['B'] + aspects['C'])
w_c = avail_w - w_b
disp_h_b = w_b / aspects['B']
disp_h_c = w_c / aspects['C']
h_bc = max(disp_h_b, disp_h_c)

# Row 2: equal halves
w_d = avail_w / 2
w_e = avail_w / 2
disp_h_d = w_d / aspects['D']
disp_h_e = w_e / aspects['E']
h_de = max(disp_h_d, disp_h_e, 1.8)  # Ensure DE row is tall enough to read labels

# Total figure height
fig_h = margin * 2 + disp_h_a + h_bc + h_de + 2 * gap

# Bottom coordinates (from bottom of figure, in inches)
bottom_e = margin
bottom_d = margin
bottom_c = margin + h_de + gap + (h_bc - disp_h_c) / 2
bottom_b = margin + h_de + gap + (h_bc - disp_h_b) / 2
bottom_a = margin + h_de + gap + h_bc + gap

left_a = margin
left_b = margin
left_c = margin + w_b
left_d = margin
left_e = margin + avail_w / 2

coords = {
    'A': (left_a / fig_w, bottom_a / fig_h, disp_w_a / fig_w, disp_h_a / fig_h),
    'B': (left_b / fig_w, bottom_b / fig_h, w_b / fig_w, disp_h_b / fig_h),
    'C': (left_c / fig_w, bottom_c / fig_h, w_c / fig_w, disp_h_c / fig_h),
    'D': (left_d / fig_w, bottom_d / fig_h, w_d / fig_w, disp_h_d / fig_h),
    'E': (left_e / fig_w, bottom_e / fig_h, w_e / fig_w, disp_h_e / fig_h),
}

images = {'A': img_2a, 'B': img_2b, 'C': img_2c, 'D': img_2d, 'E': img_2e}

# ======================== Assemble Figure ========================
fig = plt.figure(figsize=(fig_w, fig_h))

for label, (left, bottom, width, height) in coords.items():
    ax = fig.add_axes([left, bottom, width, height])
    ax.imshow(np.array(images[label]))
    ax.axis('off')
    
    # Panel label (A-E)
    ax.text(0.015, 0.98, label, transform=ax.transAxes,
            fontsize=14, fontweight='bold', va='top', ha='left',
            color='black',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                     alpha=0.9, edgecolor='none', linewidth=0))

# NOTE: Figure title removed for journal submission — JTM adds figure numbers during production
# fig.suptitle('Figure 2. Cell-Cell Communication Landscape in MPR vs Non-MPR NSCLC',
#              fontsize=13, fontweight='bold', y=0.98)

fig.savefig(OUT_DIR / 'figure2_final_combined.png', dpi=300, bbox_inches='tight', pad_inches=0.08)
fig.savefig(OUT_DIR / 'figure2_final_combined.pdf', bbox_inches='tight', pad_inches=0.08)
print('Saved: figure2_final_combined.png / .pdf')
