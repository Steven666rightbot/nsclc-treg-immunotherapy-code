# Figure 4: Spatial Validation of Treg Subtype Distribution and MHC-I Module

## Background
- Figure 2 (CellChat) showed MHC-I pathways predominantly target CCR8+ Treg in responders
- MIF/CD99 pathways target MKI67+ Treg in non-responders
- Figure 4 validates these findings spatially using Visium data (GSE189487)

## Data
- **Spatial metadata**: D:/research/tomato/results/spatial/spatial_spot_metadata.csv
  - Columns: in_tissue, array_row, array_col, sample_id, grade, stromal_score, tumor_score, region, ccr8_score, mki67_score, ccr8_score_z, mki67_score_z
  - 12,245 spots across 13 samples (7 Low Grade NSCLC, 3 High Grade NSCLC, 3 PDAC)
- **Full expression**: D:/research/tomato/results/spatial/spatial_processed.h5ad (scanpy AnnData)
  - 20,346 genes × 12,245 spots
  - Includes: HLA-A, HLA-B, HLA-C, HLA-E, B2M, TAP1, TAP2 (MHC-I); MIF, CD74; CXCL13; FOXP3, IL2RA, MKI67
- **Samples with good hard_stroma + tumor_core coverage**: LG_1, LG_2, LG_3, LG_4, LG_5, LG_6, LG_7, HG_3, PDAC_1, PDAC_2, PDAC_3

## Panel Design

### Panel 4A — CCR8+ Treg Spatial Distribution
- **Type**: 2×2 grid of spatial feature plots
- **Sections**: LG_1 (low grade), LG_6 (low grade), HG_3 (high grade), PDAC_1 (PDAC)
- **Color**: ccr8_score_z, diverging colormap (blue=low, white=mid, red=high)
- **Style**: 
  - Each plot: scatter of array_col vs array_row, colored by ccr8_score_z
  - Spot size = 15 (small), alpha = 0.8
  - Show region boundaries as convex hull outlines:
    - hard_stroma: red (#E74C3C) dashed
    - tumor_core: blue (#3498DB) dashed
  - Title per subplot = sample_id + grade
  - Shared colorbar with label "CCR8+ Treg Score (z)"
- **Font**: Arial, 300 DPI output

### Panel 4B — MKI67+ Treg Spatial Distribution
- **Same layout as 4A** (same 4 sections)
- **Color**: mki67_score_z, diverging colormap
- **Title**: per sample
- Same region boundaries

### Panel 4C — Regional Comparison Violin Plot
- **Type**: Grouped violin plot
- **Data**: All spots with region in ['hard_stroma', 'tumor_core']
- **X-axis**: Region (hard_stroma, tumor_core), grouped by score type (CCR8+, MKI67+)
- **Y-axis**: z-score
- **Colors**: 
  - CCR8+ hard_stroma: #E74C3C (red)
  - CCR8+ tumor_core: #C0392B (dark red)
  - MKI67+ hard_stroma: #3498DB (blue)
  - MKI67+ tumor_core: #2980B9 (dark blue)
- **Statistics**: Wilcoxon rank-sum test, asterisks (*** p<0.001)
- **Add** individual points (small jitter, alpha=0.3)

### Panel 4D — MHC-I Module Spatial Distribution
- **Type**: Same 2×2 grid layout as 4A/4B
- **Computation**: 
  1. Load h5ad
  2. Compute MHC-I module score = mean normalized expression of ['HLA-A', 'HLA-B', 'HLA-C', 'B2M', 'TAP1']
  3. z-score normalize across all spots
- **Color**: MHC-I module z-score, diverging colormap
- **Title**: "MHC-I Module Score"

### Panel 4E (Optional) — Correlation Scatter
- **Type**: Scatter plot
- **X-axis**: MHC-I module z-score per spot
- **Y-axis**: ccr8_score_z
- **Color**: Sample (use distinct colors for LG_1, LG_6, HG_3, PDAC_1)
- **Add**: Spearman correlation (rho and p-value) as text annotation
- **Add**: Linear regression line with confidence band

## Style Requirements
- Arial font throughout
- 300 DPI output
- PDF + PNG for all panels
- Color-blind friendly (avoid red-green)
- Consistent axis sizes across panels 4A, 4B, 4D
- Region boundary convex hulls via scipy.spatial.ConvexHull

## Output
Save to D:/research/cucumber/fig4/:
- fig4a_ccr8_spatial.pdf + .png
- fig4b_mki67_spatial.pdf + .png
- fig4c_regional_violin.pdf + .png
- fig4d_mhci_spatial.pdf + .png
- fig4e_correlation.pdf + .png (optional)
- figure4_all.py (combined script)
- FIGURE4_RESULTS.md (results text + figure legends — I will write this)

## Python Environment
- Use: source D:/research/tomato/venv_scenic/Scripts/activate
- Packages: scanpy, pandas, numpy, matplotlib, scipy, seaborn
- Working dir: D:/research/cucumber/fig4
