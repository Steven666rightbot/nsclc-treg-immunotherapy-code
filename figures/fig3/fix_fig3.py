#!/usr/bin/env python3
"""Fix color scheme and label overlap issues in Figure 3 scripts."""

import os
import re

# ── fig3b_hsf1_heatmap.py ────────────────────────────────────────────────────
path_b = 'fig3b_hsf1_heatmap.py'
with open(path_b, 'r', encoding='utf-8') as f:
    txt_b = f.read()

# Dimensions & spacing
old = """    row_height = 0.45  # inches per gene row
    heatmap_width = 8.0
    label_width = 2.0        # space for gene labels on the right
    annotation_height = 0.35  # top annotation bar"""
new = """    row_height = 0.5  # inches per gene row
    heatmap_width = 8.5
    label_width = 2.5        # space for gene labels on the right
    annotation_height = 0.5  # top annotation bar"""
assert old in txt_b, "fig3b: dimension block not found"
txt_b = txt_b.replace(old, new)

# Annotation bar rectangles (increase bar height)
old = """        ax_annot.add_patch(Rectangle(
            (ci / n_cells, 0.55), 1 / n_cells, 0.45,
            facecolor=color, edgecolor='none', linewidth=0
        ))
        # Subtype marker (bottom half) - use a lighter shade of same color
        sub_color = '#E8F5E9' if valid_responses[ci] == 'Responder' else '#F3E5F5'
        ax_annot.add_patch(Rectangle(
            (ci / n_cells, 0.05), 1 / n_cells, 0.45,
            facecolor=sub_color, edgecolor='none', linewidth=0
        ))"""
new = """        ax_annot.add_patch(Rectangle(
            (ci / n_cells, 0.5), 1 / n_cells, 0.5,
            facecolor=color, edgecolor='none', linewidth=0
        ))
        # Subtype marker (bottom half) - use a lighter shade of same color
        sub_color = '#E8F5E9' if valid_responses[ci] == 'Responder' else '#F3E5F5'
        ax_annot.add_patch(Rectangle(
            (ci / n_cells, 0.0), 1 / n_cells, 0.5,
            facecolor=sub_color, edgecolor='none', linewidth=0
        ))"""
assert old in txt_b, "fig3b: Rectangle block not found"
txt_b = txt_b.replace(old, new)

# Gene label font size
old = "    ax_heat.set_yticklabels(gene_labels_clustered, fontsize=8, fontname='Arial')"
new = "    ax_heat.set_yticklabels(gene_labels_clustered, fontsize=7.5, fontname='Arial')"
assert old in txt_b, "fig3b: yticklabels line not found"
txt_b = txt_b.replace(old, new)

# Separator line visibility
old = "        ax_heat.axvline(n_resp - 0.5, color='#333333', linewidth=1.2, linestyle='-', alpha=0.6)"
new = "        ax_heat.axvline(n_resp - 0.5, color='#333333', linewidth=2.0, linestyle='-', alpha=0.9)"
assert old in txt_b, "fig3b: axvline line not found"
txt_b = txt_b.replace(old, new)

with open(path_b, 'w', encoding='utf-8') as f:
    f.write(txt_b)
print(f"[DONE] Updated {path_b}")

# ── fig3_ad_tf_analysis.py ───────────────────────────────────────────────────
path_ad = 'fig3_ad_tf_analysis.py'
with open(path_ad, 'r', encoding='utf-8') as f:
    txt_ad = f.read()

# Panel 3A standalone: increase figure width
old = "    fig, ax = plt.subplots(figsize=(5.5, 4.0))"
new = "    fig, ax = plt.subplots(figsize=(6.2, 4.0))"
assert old in txt_ad, "fig3ad: figsize line not found"
txt_ad = txt_ad.replace(old, new)

# Panel 3A standalone: increase star offset
old = "        offset = 0.02 * (top['effect_size'].max() - top['effect_size'].min())"
new = "        offset = 0.03 * (top['effect_size'].max() - top['effect_size'].min())"
assert old in txt_ad, "fig3ad: offset line not found"
txt_ad = txt_ad.replace(old, new)

# Panel 3A standalone: increase xlim padding
old = """    ax.set_xlim(left=top['effect_size'].min() - 0.15 * abs(top['effect_size'].min()),
                right=top['effect_size'].max() + 0.15 * abs(top['effect_size'].max()))"""
new = """    ax.set_xlim(left=top['effect_size'].min() - 0.2 * abs(top['effect_size'].min()),
                right=top['effect_size'].max() + 0.2 * abs(top['effect_size'].max()))"""
assert old in txt_ad, "fig3ad: xlim block not found"
txt_ad = txt_ad.replace(old, new)

# Panel 3A standalone: explicit margins instead of tight_layout
old = """    plt.tight_layout()
    return fig, ax, top"""
new = """    fig.subplots_adjust(left=0.22, right=0.95, top=0.92, bottom=0.12)
    return fig, ax, top"""
assert old in txt_ad, "fig3ad: tight_layout block not found"
txt_ad = txt_ad.replace(old, new)

# Panel 3D standalone: increase figure size
old = "    fig, axes = plt.subplots(1, 3, figsize=(8.0, 3.2), sharey=True)"
new = "    fig, axes = plt.subplots(1, 3, figsize=(8.5, 3.5), sharey=True)"
assert old in txt_ad, "fig3ad: violin figsize not found"
txt_ad = txt_ad.replace(old, new)

# Panel 3D standalone: fix x-axis tick labels (NR instead of N)
old = "                tick_labels.append(f'{subtype_labels[si]}\\n{resp[0]}')"
new = """                resp_label = 'R' if resp == 'Responder' else 'NR'
                tick_labels.append(f'{subtype_labels[si]}\\n{resp_label}')"""
assert old in txt_ad, "fig3ad: tick_labels append not found"
txt_ad = txt_ad.replace(old, new)

# Panel 3D standalone: stagger statistical annotation heights
old = """                y_max = np.percentile(all_vals, 95)
                y_range = np.percentile(all_vals, 95) - np.percentile(all_vals, 5)
                y_pos = y_max + 0.1 * y_range"""
new = """                y_max = np.percentile(all_vals, 95)
                y_range = np.percentile(all_vals, 95) - np.percentile(all_vals, 5)
                base_offset = 0.1 if si == 0 else 0.22
                y_pos = y_max + base_offset * y_range"""
assert old in txt_ad, "fig3ad: y_pos block not found"
txt_ad = txt_ad.replace(old, new)

# Combined panel_a_sub: increase star offset
old = """            xrange = top['effect_size'].max() - top['effect_size'].min()
            offset = 0.02 * xrange"""
new = """            xrange = top['effect_size'].max() - top['effect_size'].min()
            offset = 0.03 * xrange"""
assert old in txt_ad, "fig3ad: combined offset not found"
txt_ad = txt_ad.replace(old, new)

# Combined ax_a1 xlim padding
old = """    ax_a1.set_xlim(left=top_a['effect_size'].min() - 0.15 * abs(top_a['effect_size'].min()),
                   right=top_a['effect_size'].max() + 0.15 * abs(top_a['effect_size'].max()))"""
new = """    ax_a1.set_xlim(left=top_a['effect_size'].min() - 0.2 * abs(top_a['effect_size'].min()),
                   right=top_a['effect_size'].max() + 0.2 * abs(top_a['effect_size'].max()))"""
assert old in txt_ad, "fig3ad: ax_a1 xlim not found"
txt_ad = txt_ad.replace(old, new)

# Combined ax_a2 xlim padding
old = """    ax_a2.set_xlim(left=top_b['effect_size'].min() - 0.15 * abs(top_b['effect_size'].min()),
                   right=top_b['effect_size'].max() + 0.15 * abs(top_b['effect_size'].max()))"""
new = """    ax_a2.set_xlim(left=top_b['effect_size'].min() - 0.2 * abs(top_b['effect_size'].min()),
                   right=top_b['effect_size'].max() + 0.2 * abs(top_b['effect_size'].max()))"""
assert old in txt_ad, "fig3ad: ax_a2 xlim not found"
txt_ad = txt_ad.replace(old, new)

# Combined violin: fix x-axis tick labels
old = "                tick_labels.append(f'{subtype_labels[si]}\\n{resp[0]}')"
new = """                resp_label = 'R' if resp == 'Responder' else 'NR'
                tick_labels.append(f'{subtype_labels[si]}\\n{resp_label}')"""
# This appears twice; replace both
assert old in txt_ad, "fig3ad: combined tick_labels not found"
txt_ad = txt_ad.replace(old, new)

# Combined violin: stagger statistical annotation heights
old = """                y_max = np.percentile(all_vals, 95)
                y_range = np.percentile(all_vals, 95) - np.percentile(all_vals, 5)
                y_pos = y_max + 0.12 * y_range"""
new = """                y_max = np.percentile(all_vals, 95)
                y_range = np.percentile(all_vals, 95) - np.percentile(all_vals, 5)
                base_offset = 0.12 if si == 0 else 0.24
                y_pos = y_max + base_offset * y_range"""
assert old in txt_ad, "fig3ad: combined y_pos not found"
txt_ad = txt_ad.replace(old, new)

with open(path_ad, 'w', encoding='utf-8') as f:
    f.write(txt_ad)
print(f"[DONE] Updated {path_ad}")

# ── fig3c_regulatory_network.py ──────────────────────────────────────────────
path_c = 'fig3c_regulatory_network.py'
with open(path_c, 'r', encoding='utf-8') as f:
    txt_c = f.read()

# Extend ylim to accommodate legend
old = "ax.set_ylim(-1.0, 9.0)"
new = "ax.set_ylim(-1.5, 9.0)"
assert old in txt_c, "fig3c: ylim not found"
txt_c = txt_c.replace(old, new)

# Move legend to bottom left
old = "legend_x = 8.8\nlegend_y = 1.3"
new = "legend_x = 0.5\nlegend_y = -0.2"
assert old in txt_c, "fig3c: legend coords not found"
txt_c = txt_c.replace(old, new)

# Shift 'GC Tfh-like' label up to avoid BCL6→MHC-II arrow
old = "annotation_text(tf_x + 0.8, 5.5, 'GC Tfh-like',\n                fontsize=6.5, color='#6c3483', ha='left')"
new = "annotation_text(tf_x + 0.85, 5.75, 'GC Tfh-like',\n                fontsize=6.5, color='#6c3483', ha='left')"
assert old in txt_c, "fig3c: GC Tfh-like not found"
txt_c = txt_c.replace(old, new)

# Shift 'MHC-II Ag\\npresentation' up
old = "annotation_text(tf_x + 0.75, 4.0, 'MHC-II Ag\\npresentation',\n                fontsize=6.0, color='#1e8449', ha='left')"
new = "annotation_text(tf_x + 0.85, 4.2, 'MHC-II Ag\\npresentation',\n                fontsize=6.0, color='#1e8449', ha='left')"
assert old in txt_c, "fig3c: MHC-II Ag not found"
txt_c = txt_c.replace(old, new)

# Shift 'Proliferation' down to avoid MAX/MYC arrow starts
old = "annotation_text(tf_x + 0.75, 2.7, 'Proliferation',\n                fontsize=6.5, color='#148f77', ha='left')"
new = "annotation_text(tf_x + 0.85, 2.45, 'Proliferation',\n                fontsize=6.5, color='#148f77', ha='left')"
assert old in txt_c, "fig3c: Proliferation not found"
txt_c = txt_c.replace(old, new)

with open(path_c, 'w', encoding='utf-8') as f:
    f.write(txt_c)
print(f"[DONE] Updated {path_c}")

print("\nAll files updated successfully.")
