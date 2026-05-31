"""
Spatial Gradient Analysis for PDAC Visium data
Simplified alternative to SLOPER using scanpy + scipy + matplotlib.
Computes spatial gradients and streamlines for key mechanosensitive genes.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, SmoothBivariateSpline
from scipy.ndimage import gaussian_filter

# ================================================================
# Configuration
# ================================================================
BASE_DIR = r'D:\Research\tomato'
DATA_DIR = os.path.join(BASE_DIR, 'data', 'spatial', 'pdac_visium', 'extracted')
RESULTS_DIR = os.path.join(BASE_DIR, 'results', 'spatial')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures', 'spatial_gradient')
os.makedirs(FIGURES_DIR, exist_ok=True)

# Representative sample for spatial plots
REP_SAMPLE = 'PDAC_1'

# Target genes for gradient analysis
TARGET_GENES = {
    'ECM': ['COL1A2', 'COL1A1', 'POSTN', 'FN1'],
    'Mechanosensing': ['ITGB1', 'CD44', 'DDR1', 'DDR2'],
    'YAP_TAZ': ['YAP1', 'WWTR1', 'TEAD1', 'CTGF', 'CYR61'],
    'Treg_subtype': ['CCR8', 'LAYN', 'BATF', 'CTLA4', 'MKI67', 'TOP2A'],
    'Contractile': ['ACTA2', 'MYL9', 'MYH9', 'TPM1', 'TAGLN'],
}

# Grid resolution for interpolation
GRID_RES = 100
SIGMA = 1.5  # Gaussian smoothing sigma

# ================================================================
# 1. Load representative sample
# ================================================================
print("=" * 60)
print(f"Loading {REP_SAMPLE}...")
print("=" * 60)

adata = sc.read_visium(os.path.join(DATA_DIR, REP_SAMPLE))
adata.var_names_make_unique()
adata = adata[adata.obs['in_tissue'] == 1].copy()

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print(f"  Spots: {adata.n_obs}, Genes: {adata.n_vars}")

# Get spatial coordinates
spatial = adata.obsm['spatial']  # shape (n_spots, 2)
coords = spatial.copy()

# Center and scale coords for better interpolation
coords_min = coords.min(axis=0)
coords_max = coords.max(axis=0)
coords_norm = (coords - coords_min) / (coords_max - coords_min)

# ================================================================
# 2. Gradient computation function
# ================================================================
def compute_gradient(coords_norm, values, grid_res=GRID_RES, sigma=SIGMA):
    """
    Compute spatial gradient using interpolation + finite differences.
    Returns: grid_x, grid_y, grid_z (smoothed), grad_x, grad_y
    """
    # Create regular grid
    x = np.linspace(0, 1, grid_res)
    y = np.linspace(0, 1, grid_res)
    grid_x, grid_y = np.meshgrid(x, y)
    
    # Interpolate expression to grid
    grid_z = griddata(coords_norm, values, (grid_x, grid_y), method='cubic', fill_value=0)
    grid_z = np.nan_to_num(grid_z)
    
    # Gaussian smoothing (simulate SLOPER's enhancement)
    grid_z_smooth = gaussian_filter(grid_z, sigma=sigma)
    
    # Compute gradient
    grad_y, grad_x = np.gradient(grid_z_smooth)
    
    return grid_x, grid_y, grid_z_smooth, grad_x, grad_y

# ================================================================
# 3. Plotting function
# ================================================================
def plot_gradient_streamlines(grid_x, grid_y, grid_z, grad_x, grad_y, 
                               gene_name, category, sample_id,
                               output_dir):
    """Plot enhanced expression + gradient streamlines."""
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # --- Panel A: Enhanced expression heatmap ---
    ax = axes[0]
    vmax = np.percentile(grid_z, 99)
    vmin = np.percentile(grid_z, 1)
    im = ax.imshow(grid_z, origin='lower', cmap='viridis', 
                   vmin=vmin, vmax=vmax,
                   extent=[grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()])
    ax.set_title(f'{gene_name}\nEnhanced Expression', fontsize=13, fontweight='bold')
    ax.set_xlabel('Normalized X')
    ax.set_ylabel('Normalized Y')
    plt.colorbar(im, ax=ax, shrink=0.7, label='log1p(norm counts)')
    ax.set_aspect('equal')
    
    # --- Panel B: Gradient streamlines ---
    ax = axes[1]
    
    # Background: expression
    vmax2 = np.percentile(grid_z, 95)
    im2 = ax.imshow(grid_z, origin='lower', cmap='YlOrRd', 
                    vmin=0, vmax=vmax2,
                    extent=[grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()],
                    alpha=0.7)
    
    # Streamlines
    # Subsample for cleaner visualization
    step = max(1, grid_x.shape[0] // 25)
    x_sub = grid_x[::step, ::step]
    y_sub = grid_y[::step, ::step]
    u_sub = grad_x[::step, ::step]
    v_sub = grad_y[::step, ::step]
    
    # Mask low-gradient regions
    mag = np.sqrt(u_sub**2 + v_sub**2)
    mask = mag > np.percentile(mag, 30)
    
    # Plot streamlines
    strm = ax.streamplot(x_sub[0, :], y_sub[:, 0], u_sub.T, v_sub.T,
                         color='white', linewidth=1.2, density=1.5,
                         arrowsize=1.5, arrowstyle='->')
    
    # Overlay quiver for strong gradients
    ax.quiver(x_sub[mask], y_sub[mask], u_sub[mask], v_sub[mask],
              color='cyan', scale=20, width=0.003, alpha=0.8)
    
    ax.set_title(f'{gene_name}\nSpatial Gradient (→ high expr)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Normalized X')
    ax.set_ylabel('Normalized Y')
    plt.colorbar(im2, ax=ax, shrink=0.7, label='Expression')
    ax.set_aspect('equal')
    
    plt.tight_layout()
    
    out_path = os.path.join(output_dir, f'{sample_id}_{category}_{gene_name}_gradient.png')
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")

# ================================================================
# 4. Run analysis for all target genes
# ================================================================
print("\n" + "=" * 60)
print("Computing spatial gradients...")
print("=" * 60)

results_summary = []

for category, genes in TARGET_GENES.items():
    print(f"\n--- {category} ---")
    for gene in genes:
        if gene not in adata.var_names:
            print(f"  {gene}: NOT FOUND, skipping")
            continue
        
        expr = adata[:, gene].X.toarray().flatten()
        
        # Skip if too lowly expressed
        if expr.max() < 0.1:
            print(f"  {gene}: too low expression, skipping")
            continue
        
        print(f"  {gene}: computing gradient...", end=" ", flush=True)
        
        grid_x, grid_y, grid_z, grad_x, grad_y = compute_gradient(
            coords_norm, expr, grid_res=GRID_RES, sigma=SIGMA
        )
        
        # Summary statistics
        grad_mag = np.sqrt(grad_x**2 + grad_y**2)
        results_summary.append({
            'gene': gene,
            'category': category,
            'mean_expr': expr.mean(),
            'max_expr': expr.max(),
            'mean_gradient_magnitude': grad_mag.mean(),
            'max_gradient_magnitude': grad_mag.max(),
        })
        
        # Plot
        plot_gradient_streamlines(grid_x, grid_y, grid_z, grad_x, grad_y,
                                   gene, category, REP_SAMPLE, FIGURES_DIR)
        print("done")

# ================================================================
# 5. Combined multi-gene comparison plot
# ================================================================
print("\n" + "=" * 60)
print("Generating combined comparison plot...")
print("=" * 60)

# Pick top 3 genes from each category that have highest mean gradient
df_summary = pd.DataFrame(results_summary)
df_summary.to_csv(os.path.join(RESULTS_DIR, 'spatial_gradient_summary.csv'), index=False)

top_genes = []
for cat in TARGET_GENES.keys():
    cat_df = df_summary[df_summary['category'] == cat]
    if len(cat_df) > 0:
        top = cat_df.nlargest(min(2, len(cat_df)), 'mean_gradient_magnitude')
        top_genes.extend(top['gene'].tolist())

n_genes = len(top_genes)
if n_genes > 0:
    ncols = 3
    nrows = (n_genes + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 5))
    axes = axes.flatten() if nrows > 1 else [axes] if ncols == 1 else axes
    
    for idx, gene in enumerate(top_genes):
        ax = axes[idx]
        expr = adata[:, gene].X.toarray().flatten()
        grid_x, grid_y, grid_z, grad_x, grad_y = compute_gradient(
            coords_norm, expr, grid_res=GRID_RES, sigma=SIGMA
        )
        
        vmax = np.percentile(grid_z, 95)
        im = ax.imshow(grid_z, origin='lower', cmap='YlOrRd', vmin=0, vmax=vmax,
                       extent=[grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()])
        
        step = max(1, grid_x.shape[0] // 20)
        strm = ax.streamplot(grid_x[0, ::step], grid_y[::step, 0], 
                             grad_x[::step, ::step].T, grad_y[::step, ::step].T,
                             color='white', linewidth=1.0, density=1.2, arrowsize=1.2)
        
        ax.set_title(gene, fontsize=12, fontweight='bold')
        ax.set_aspect('equal')
        ax.axis('off')
    
    # Hide unused subplots
    for idx in range(n_genes, len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle('Spatial Gradients of Key Genes\n(arrows point toward higher expression)', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    out_combined = os.path.join(FIGURES_DIR, f'{REP_SAMPLE}_combined_gradient_streamlines.png')
    plt.savefig(out_combined, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_combined}")

# ================================================================
# 6. Gradient magnitude comparison barplot
# ================================================================
print("\nGenerating gradient magnitude comparison...")

fig, ax = plt.subplots(figsize=(10, 6))
df_plot = df_summary.sort_values('mean_gradient_magnitude', ascending=True)
colors = {'ECM': '#e74c3c', 'Mechanosensing': '#3498db', 
          'YAP_TAZ': '#9b59b6', 'Treg_subtype': '#2ecc71', 'Contractile': '#f39c12'}
bar_colors = [colors.get(c, 'gray') for c in df_plot['category']]

bars = ax.barh(df_plot['gene'], df_plot['mean_gradient_magnitude'], color=bar_colors)
ax.set_xlabel('Mean Gradient Magnitude', fontsize=11)
ax.set_title('Spatial Gradient Strength by Gene', fontsize=13, fontweight='bold')

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors[c], label=c) for c in colors if c in df_plot['category'].values]
ax.legend(handles=legend_elements, loc='lower right')

plt.tight_layout()
out_bar = os.path.join(FIGURES_DIR, 'gradient_magnitude_comparison.png')
plt.savefig(out_bar, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {out_bar}")

# ================================================================
# Done
# ================================================================
print("\n" + "=" * 60)
print(f"DONE! All results in: {FIGURES_DIR}")
print("=" * 60)
print(f"\nSummary: {len(results_summary)} genes analyzed")
print(df_summary.to_string(index=False))
