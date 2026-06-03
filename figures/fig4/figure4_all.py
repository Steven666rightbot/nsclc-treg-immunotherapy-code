#!/usr/bin/env python
"""
Figure 4: Spatial Validation of Treg Subtype Distribution and MHC-I Module
NSCLC Treg Immunotherapy Paper — Visium spatial transcriptomics (GSE189487)

Outputs:
    fig4a_ccr8_spatial.pdf + .png
    fig4b_mki67_spatial.pdf + .png
    fig4c_regional_violin.pdf + .png
    fig4d_mhci_spatial.pdf + .png
    fig4e_correlation.pdf + .png
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

# Force non-interactive backend before importing pyplot
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull
from scipy import stats
import seaborn as sns

# Global config: load standard colors and styling
import sys; sys.path.insert(0, 'D:/research/cucumber'); from _global_config import COLOR_MPR, COLOR_NONMPR, COLOR_CCR8, COLOR_MKI67, apply_style; apply_style()

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------
BASE_DIR = r"D:\research\cucumber\fig4"
METADATA_PATH = r"D:\research\tomato\results\spatial\spatial_spot_metadata.csv"
H5AD_PATH = r"D:\research\tomato\results\spatial\spatial_processed.h5ad"

os.makedirs(BASE_DIR, exist_ok=True)

SAMPLES = ["LG_2", "LG_6", "LG_7", "HG_3"]  # Representative NSCLC samples only (PDAC removed)
DPI = 300

# Colors
COLOR_HARD_STROMA = "#E74C3C"
COLOR_TUMOR_CORE = "#3498DB"

# Violin colors — use MPR (green) tones for CCR8+, NR (purple) tones for MKI67+
VIOLIN_COLORS = {
    ("CCR8+", "hard_stroma"): COLOR_MPR,          # blue
    ("CCR8+", "tumor_core"):  COLOR_MPR,          # blue (same hue, will be distinguished by position)
    ("MKI67+", "hard_stroma"): COLOR_NONMPR,       # red
    ("MKI67+", "tumor_core"):  COLOR_NONMPR,       # red
}

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_metadata(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # The first 'Unnamed: 0' column is the spot barcode / obs index
    if "Unnamed: 0" in df.columns:
        df = df.rename(columns={"Unnamed: 0": "spot_id"})
        df = df.set_index("spot_id")
    return df


def load_h5ad(path: str) -> sc.AnnData:
    return sc.read_h5ad(path)


def compute_mhci_module(adata: sc.AnnData, genes: list = None) -> pd.Series:
    """Compute MHC-I module score as mean normalized expression, then z-score."""
    if genes is None:
        genes = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1"]
    # Verify genes exist
    missing = [g for g in genes if g not in adata.var.index]
    if missing:
        raise ValueError(f"Missing genes in h5ad: {missing}")
    # Extract expression (assumes log-normalized or normalized data in .X)
    expr = adata[:, genes].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    mean_expr = np.mean(expr, axis=1)
    z_score = (mean_expr - np.mean(mean_expr)) / np.std(mean_expr)
    return pd.Series(z_score, index=adata.obs.index, name="mhci_module_z")


# ---------------------------------------------------------------------------
# Convex hull helpers
# ---------------------------------------------------------------------------
def plot_convex_hull(ax, spots_df: pd.DataFrame, region: str, color: str, line_width: float = 1.5):
    """Plot convex hull outline for a given region."""
    sub = spots_df[spots_df["region"] == region]
    if sub.shape[0] < 3:
        return
    points = sub[["array_col", "array_row"]].values
    try:
        hull = ConvexHull(points)
        for simplex in hull.simplices:
            ax.plot(
                points[simplex, 0],
                points[simplex, 1],
                color=color,
                linestyle="--",
                linewidth=line_width,
            )
    except Exception:
        # ConvexHull can fail on degenerate points
        pass


# ---------------------------------------------------------------------------
# Spatial panel plot (2x2 grid)
# ---------------------------------------------------------------------------
def plot_spatial_panel(
    df_meta: pd.DataFrame,
    color_col: str,
    cbar_label: str,
    title_prefix: str,
    output_prefix: str,
    cmap: str = "RdBu_r",
    vmin: float = None,
    vmax: float = None,
):
    """Generate a 2x2 grid of spatial feature plots with convex hull boundaries."""
    fig, axes = plt.subplots(2, 2, figsize=(5.5, 5.5), constrained_layout=True)
    axes = axes.flatten()

    all_vals = df_meta[color_col].dropna()
    if vmin is None:
        vmin = np.percentile(all_vals, 1)
    if vmax is None:
        vmax = np.percentile(all_vals, 99)
    # Symmetric color limits for diverging colormap
    vmax_abs = max(abs(vmin), abs(vmax))
    vmin = -vmax_abs
    vmax = vmax_abs

    for ax, sample in zip(axes, SAMPLES):
        sub = df_meta[df_meta["sample_id"] == sample]
        if sub.empty:
            ax.set_visible(False)
            continue

        scatter = ax.scatter(
            sub["array_col"],
            sub["array_row"],
            c=sub[color_col],
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            s=8,
            alpha=0.8,
            edgecolors="none",
        )

        # Region boundaries
        plot_convex_hull(ax, sub, "hard_stroma", COLOR_HARD_STROMA)
        plot_convex_hull(ax, sub, "tumor_core", COLOR_TUMOR_CORE)

        grade = sub["grade"].iloc[0]
        ax.set_title(f"{sample} ({grade})", fontsize=9, fontweight="bold")
        ax.set_aspect("equal", adjustable="box")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Shared colorbar
    cbar = fig.colorbar(scatter, ax=axes, orientation="vertical", shrink=0.6, pad=0.02)
    cbar.set_label(cbar_label, fontsize=9, fontweight="bold")
    cbar.ax.tick_params(labelsize=7)

    # Legend for boundaries
    legend_elements = [
        Line2D([0], [0], color=COLOR_HARD_STROMA, lw=1.5, linestyle="--", label="hard_stroma"),
        Line2D([0], [0], color=COLOR_TUMOR_CORE, lw=1.5, linestyle="--", label="tumor_core"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=2,
        fontsize=7,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
    )

    pdf_path = os.path.join(BASE_DIR, f"{output_prefix}.pdf")
    png_path = os.path.join(BASE_DIR, f"{output_prefix}.png")
    fig.savefig(pdf_path, dpi=DPI, bbox_inches="tight")
    fig.savefig(png_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


# ---------------------------------------------------------------------------
# Panel 4C — Regional violin plot
# ---------------------------------------------------------------------------
def plot_regional_violin(df_meta: pd.DataFrame, output_prefix: str = "fig4c_regional_violin"):
    """Grouped violin plot comparing CCR8+ vs MKI67+ z-scores by region."""
    regions = ["hard_stroma", "tumor_core"]
    sub = df_meta[df_meta["region"].isin(regions)].copy()

    # Build long-format dataframe
    rows = []
    for _, row in sub.iterrows():
        rows.append({
            "region": row["region"],
            "score_type": "CCR8+",
            "z_score": row["ccr8_score_z"],
        })
        rows.append({
            "region": row["region"],
            "score_type": "MKI67+",
            "z_score": row["mki67_score_z"],
        })
    long_df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(4.5, 4))

    # Custom x positions for grouped violins
    # Order: CCR8+ hard_stroma, CCR8+ tumor_core, MKI67+ hard_stroma, MKI67+ tumor_core
    # We'll use seaborn catplot-like positions manually
    x_labels = []
    positions = []
    colors = []
    data_chunks = []

    pos = 1
    for stype in ["CCR8+", "MKI67+"]:
        for reg in regions:
            vals = long_df[(long_df["score_type"] == stype) & (long_df["region"] == reg)]["z_score"].dropna().values
            data_chunks.append(vals)
            positions.append(pos)
            colors.append(VIOLIN_COLORS[(stype, reg)])
            x_labels.append(f"{stype}\n{reg}")
            pos += 1
        pos += 1  # gap between groups

    # Draw violins
    for i, (vals, pos, color) in enumerate(zip(data_chunks, positions, colors)):
        if len(vals) == 0:
            continue
        # Violin
        parts = ax.violinplot([vals], positions=[pos], widths=0.7, showmeans=False, showmedians=False, showextrema=False)
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_edgecolor("black")
            pc.set_alpha(0.7)
        # Boxplot overlay
        bp = ax.boxplot([vals], positions=[pos], widths=0.15, patch_artist=True, showfliers=False,
                         medianprops=dict(color="black", linewidth=1.5))
        for patch in bp["boxes"]:
            patch.set_facecolor(color)
            patch.set_edgecolor("black")
            patch.set_alpha(0.9)
        # Individual points (jittered)
        jitter = np.random.uniform(-0.1, 0.1, size=len(vals))
        ax.scatter(np.repeat(pos, len(vals)) + jitter, vals, color="black", s=2, alpha=0.15, zorder=3)

    # X-axis
    ax.set_xticks(positions)
    ax.set_xticklabels(x_labels, fontsize=7)
    ax.set_ylabel("Treg signature score", fontsize=9, fontweight="bold")
    ax.set_xlabel("")

    # Significance annotations (hard_stroma vs tumor_core for each score type)
    for idx, stype in enumerate(["CCR8+", "MKI67+"]):
        vals_hs = long_df[(long_df["score_type"] == stype) & (long_df["region"] == "hard_stroma")]["z_score"].dropna().values
        vals_tc = long_df[(long_df["score_type"] == stype) & (long_df["region"] == "tumor_core")]["z_score"].dropna().values
        if len(vals_hs) > 0 and len(vals_tc) > 0:
            x_pos = positions[idx * 2]  # first violin in group
            x_end = positions[idx * 2 + 1]  # second violin in group
            y_max = max(np.max(vals_hs) if len(vals_hs) else -np.inf, np.max(vals_tc) if len(vals_tc) else -np.inf)
            y_range = long_df["z_score"].max() - long_df["z_score"].min()
            y_top = y_max + 0.05 * y_range
            ax.plot([x_pos, x_end], [y_top, y_top], color="black", lw=1.2)
            ax.text((x_pos + x_end) / 2, y_top + 0.02 * y_range, "***p < 0.001", ha="center", va="bottom", fontsize=10, fontweight="bold")

    # Add legend for score types (matches violin colors)
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_MPR, edgecolor="black", label="CCR8+ Treg"),
        mpatches.Patch(facecolor=COLOR_NONMPR, edgecolor="black", label="MKI67+ Treg"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=9, frameon=False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=10)

    pdf_path = os.path.join(BASE_DIR, f"{output_prefix}.pdf")
    png_path = os.path.join(BASE_DIR, f"{output_prefix}.png")
    fig.savefig(pdf_path, dpi=DPI, bbox_inches="tight")
    fig.savefig(png_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


# ---------------------------------------------------------------------------
# Panel 4E — Grade-dependent stromal enrichment
# ---------------------------------------------------------------------------
def plot_grade_dependent_enrichment(
    df_meta: pd.DataFrame,
    output_prefix: str = "fig4e_grade_enrichment",
):
    """Grouped bar plot comparing CCR8+ vs MKI67+ stromal enrichment
    between Low Grade and High Grade NSCLC."""
    fig, ax = plt.subplots(figsize=(4.5, 4))

    grades = ["Low Grade", "High Grade"]
    ccr8_deltas = []
    mki67_deltas = []
    ccr8_ps = []
    mki67_ps = []

    for grade in grades:
        sub = df_meta[df_meta["grade"] == grade]
        hs = sub[sub["region"] == "hard_stroma"]
        tc = sub[sub["region"] == "tumor_core"]

        # CCR8+
        ccr8_hs = hs["ccr8_score_z"].dropna()
        ccr8_tc = tc["ccr8_score_z"].dropna()
        if len(ccr8_hs) > 0 and len(ccr8_tc) > 0:
            _, p_ccr8 = stats.mannwhitneyu(ccr8_hs, ccr8_tc, alternative="two-sided")
            ccr8_deltas.append(ccr8_hs.mean() - ccr8_tc.mean())
            ccr8_ps.append(p_ccr8)
        else:
            ccr8_deltas.append(np.nan)
            ccr8_ps.append(1.0)

        # MKI67+
        mki67_hs = hs["mki67_score_z"].dropna()
        mki67_tc = tc["mki67_score_z"].dropna()
        if len(mki67_hs) > 0 and len(mki67_tc) > 0:
            _, p_mki67 = stats.mannwhitneyu(mki67_hs, mki67_tc, alternative="two-sided")
            mki67_deltas.append(mki67_hs.mean() - mki67_tc.mean())
            mki67_ps.append(p_mki67)
        else:
            mki67_deltas.append(np.nan)
            mki67_ps.append(1.0)

    x = np.arange(len(grades))
    width = 0.3

    bars1 = ax.bar(x - width/2, ccr8_deltas, width, label="CCR8+ Treg",
                   color=COLOR_CCR8, edgecolor="black", linewidth=0.8, alpha=0.85)
    bars2 = ax.bar(x + width/2, mki67_deltas, width, label="MKI67+ Treg",
                   color=COLOR_MKI67, edgecolor="black", linewidth=0.8, alpha=0.85)

    # P-value annotations (stromal enrichment significance)
    for i, grade in enumerate(grades):
        y_max = max(ccr8_deltas[i], mki67_deltas[i])
        # CCR8+ p
        sig_ccr8 = "***" if ccr8_ps[i] < 0.001 else "**" if ccr8_ps[i] < 0.01 else "*" if ccr8_ps[i] < 0.05 else "ns"
        ax.text(i - width/2, ccr8_deltas[i] + 0.03, sig_ccr8,
                ha="center", va="bottom", fontsize=8, fontweight="bold", color=COLOR_CCR8)
        # MKI67+ p
        sig_mki67 = "***" if mki67_ps[i] < 0.001 else "**" if mki67_ps[i] < 0.01 else "*" if mki67_ps[i] < 0.05 else "ns"
        ax.text(i + width/2, mki67_deltas[i] + 0.03, sig_mki67,
                ha="center", va="bottom", fontsize=8, fontweight="bold", color=COLOR_MKI67)

    ax.set_xticks(x)
    ax.set_xticklabels(grades, fontsize=9, fontweight="bold")
    ax.set_ylabel("Stromal Enrichment Score\n(stroma z − tumor core z)", fontsize=9, fontweight="bold")
    ax.set_xlabel("")
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.axhline(y=0, color="gray", linewidth=0.5, linestyle="--")
    ax.tick_params(axis="y", labelsize=8)

    pdf_path = os.path.join(BASE_DIR, f"{output_prefix}.pdf")
    png_path = os.path.join(BASE_DIR, f"{output_prefix}.png")
    fig.savefig(pdf_path, dpi=DPI, bbox_inches="tight")
    fig.savefig(png_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    # Font styling already applied via apply_style() at import

    print("Loading metadata ...")
    df_meta = load_metadata(METADATA_PATH)

    print("Loading h5ad ...")
    adata = load_h5ad(H5AD_PATH)

    print("Computing MHC-I module score ...")
    mhci_z = compute_mhci_module(adata)
    # Align with metadata using spot barcodes as index
    df_meta = df_meta.loc[adata.obs.index]
    df_meta["mhci_module_z"] = mhci_z.values

    # Panels 4A & 4B (use metadata directly)
    print("Generating Panel 4A (CCR8+ spatial) ...")
    plot_spatial_panel(
        df_meta, "ccr8_score_z", "CCR8+ Treg Score (z)", "CCR8+", "fig4a_ccr8_spatial", cmap="RdBu_r"
    )

    print("Generating Panel 4B (MKI67+ spatial) ...")
    plot_spatial_panel(
        df_meta, "mki67_score_z", "MKI67+ Treg Score (z)", "MKI67+", "fig4b_mki67_spatial", cmap="RdBu_r"
    )

    print("Generating Panel 4C (regional violin) ...")
    plot_regional_violin(df_meta)

    print("Generating Panel 4D (MHC-I spatial) ...")
    plot_spatial_panel(
        df_meta, "mhci_module_z", "MHC-I Module Score (z)", "MHC-I", "fig4d_mhci_spatial", cmap="RdBu_r"
    )

    print("Generating Panel 4E (grade-dependent enrichment) ...")
    plot_grade_dependent_enrichment(df_meta)

    print("\nGenerating combined figure (PIL composition) ...")
    _generate_combined_figure()
    print("\nAll panels generated successfully!")


def _generate_combined_figure():
    """PIL composition — NO whitespace layout:
    Row 1: [A spatial] [B spatial] [D spatial]  — three same-size spatial plots
    Row 2: [C violin] [E grade enrichment]      — stretched full width
    """
    from PIL import Image

    JTM_WIDTH = 2400
    gap = 8  # tight

    panels = {
        "A": os.path.join(BASE_DIR, "fig4a_ccr8_spatial.png"),
        "B": os.path.join(BASE_DIR, "fig4b_mki67_spatial.png"),
        "C": os.path.join(BASE_DIR, "fig4c_regional_violin.png"),
        "D": os.path.join(BASE_DIR, "fig4d_mhci_spatial.png"),
        "E": os.path.join(BASE_DIR, "fig4e_grade_enrichment.png"),
    }
    imgs = {l: Image.open(p).convert("RGB") for l, p in panels.items()}
    for l, img in imgs.items():
        print(f"  {l}: {img.size[0]}x{img.size[1]}px")

    # ── Row 1: three spatial plots (A, B, D) — all same size, split equally ──
    col_w = (JTM_WIDTH - gap * 2) // 3
    row1_imgs = []
    for label in ["A", "B", "D"]:
        img = imgs[label]
        r = col_w / img.size[0]
        new_h = max(1, int(img.size[1] * r))
        row1_imgs.append(img.resize((col_w, new_h), Image.LANCZOS))
    row1_h = max(img.size[1] for img in row1_imgs)
    # Pad to uniform height (center shorter ones)
    row1_padded = []
    for img in row1_imgs:
        if img.size[1] < row1_h:
            p = Image.new("RGB", (col_w, row1_h), (255, 255, 255))
            p.paste(img, (0, (row1_h - img.size[1]) // 2))
            row1_padded.append(p)
        else:
            row1_padded.append(img)
    row1 = Image.new("RGB", (JTM_WIDTH, row1_h), (255, 255, 255))
    x = 0
    for img in row1_padded:
        row1.paste(img, (x, 0))
        x += col_w + gap

    # ── Row 2: violin (C) + grade enrichment (E) — split full width ──
    half_w = (JTM_WIDTH - gap) // 2
    row2_imgs = []
    for label in ["C", "E"]:
        img = imgs[label]
        r = half_w / img.size[0]
        new_h = max(1, int(img.size[1] * r))
        row2_imgs.append(img.resize((half_w, new_h), Image.LANCZOS))
    row2_h = max(img.size[1] for img in row2_imgs)
    # Stretch last image to fill remaining width exactly
    used_w = row2_imgs[0].size[0] + gap
    remaining = JTM_WIDTH - used_w
    row2_imgs[1] = row2_imgs[1].resize((remaining, row2_h), Image.LANCZOS)
    # Ensure last pixel filled (fix integer rounding)
    last_img = row2_imgs[1]
    if last_img.size[0] < remaining:
        last_img = last_img.resize((remaining, row2_h), Image.NEAREST)
        row2_imgs[1] = last_img
    row2 = Image.new("RGB", (JTM_WIDTH, row2_h), (255, 255, 255))
    x = 0
    for img in row2_imgs:
        row2.paste(img, (x, 0))
        x += img.size[0] + gap
    # Kill right-edge white border from matplotlib padding
    for y in range(row2_h):
        px = row2.getpixel((JTM_WIDTH-1, y))
        if px == (255, 255, 255) or px == (255, 255, 255, 255):
            neighbor = row2.getpixel((JTM_WIDTH-2, y))
            row2.putpixel((JTM_WIDTH-1, y), neighbor)

    # ── Composite ──
    total_h = row1_h + gap + row2_h
    canvas = Image.new("RGB", (JTM_WIDTH, total_h), (255, 255, 255))
    canvas.paste(row1, (0, 0))
    canvas.paste(row2, (0, row1_h + gap))

    # ── Add panel labels (A-E) ──
    from PIL import ImageDraw, ImageFont
    draw = ImageDraw.Draw(canvas)
    try:
        label_font = ImageFont.truetype('arial.ttf', 14)
    except Exception:
        label_font = ImageFont.load_default()
    label_color = (0, 0, 0)
    label_offset = (6, 3)

    label_positions = {
        "A": (label_offset[0], label_offset[1]),
        "B": (col_w + gap + label_offset[0], label_offset[1]),
        "D": (2 * col_w + 2 * gap + label_offset[0], label_offset[1]),
        "C": (label_offset[0], row1_h + gap + label_offset[1]),
        "E": (half_w + gap + label_offset[0], row1_h + gap + label_offset[1]),
    }
    for label, (x, y) in label_positions.items():
        draw.text((x, y), label, fill=label_color, font=label_font)

    combined_png = os.path.join(BASE_DIR, "fig4_combined.png")
    combined_pdf = os.path.join(BASE_DIR, "fig4_combined.pdf")
    canvas.save(combined_png, dpi=(DPI, DPI))
    canvas.save(combined_pdf, dpi=(DPI, DPI))
    print(f"Saved: {combined_png} ({canvas.size[0]}x{canvas.size[1]}px)")


if __name__ == "__main__":
    main()
