# Figure 5 — Design Specification for Kimi Code

## Background
- Figure 2 (CellChat) showed MHC-I/II → CCR8+ Treg in MPR, MIF → MKI67+ Treg in Non-MPR
- Figure 5 validates these findings using scTenifoldKnk in-silico KO on GSE243013 subtype-specific Treg data
- Two KOs: **CXCL13** in CCR8+ Treg (archetypal TLS chemokine), **B2M** in MKI67+ Treg (MHC-I light chain)
- B2M KO in MKI67+ Treg downregulates HLA-DRB5 (MHC-II) and CD74 (MIF receptor) → directly links CellChat predictions

## Panel Design

### Panel 5A — CXCL13 KO Volcano in CCR8+ Treg
- **Data**: `D:/Research/tomato/results/scTenifoldKnk_ko_ccr8_cxcl13/ko_CCR8_CXCL13.csv`
- **Type**: Volcano-like scatter (Z-score vs -log10(p.adj))
- **X-axis**: Z-score (from scTenifoldKnk)
- **Y-axis**: -log10(adjusted p-value)
- **Color by significance**: 
  - p.adj < 0.05 & |Z| > 2: red ("BH significant")
  - Z > 2 but p.adj ≥ 0.05: gray ("Downstream effect")
  - others: light gray ("Not significant")
- **Label top genes**: CLDND2, AL136018.1, AC034238.2, CCL3 (Z>2.5), CXCL13 itself (self-KO)
- **Title**: "CXCL13 KO in CCR8+ Treg"
  - Subtitle: "8 significant genes (BH, p.adj < 0.05)"
- **Color palette**:
  - BH significant: `#E74C3C` (red)
  - Downstream (Z>2 only): `#7F8C8D` (gray)
  - Non-significant: `#BDC3C7` (light gray)
- **Style**: Arial, 300 DPI, grid=False, size 5x5 inches

### Panel 5B — B2M KO Volcano in MKI67+ Treg
- **Data**: `D:/Research/tomato/results/scTenifoldKnk_ko_mki67_b2m/ko_MKI67_B2M_lowres.csv`
- **Type**: Same volcano layout as 5A
- **X-axis**: Z-score
- **Y-axis**: -log10(p.adj)
- **Color**: Same scheme as 5A
- **Label key genes**: PTMS, HAVCR1, **HLA-DRB5**, **CD74**, B2M itself
  - Bold the genes that directly link to CellChat findings
- **Title**: "B2M KO in MKI67+ Treg"
  - Subtitle: "19 significant genes (BH, p.adj < 0.05)"
- **Inset/annotation**: "HLA-DRB5↓, CD74↓ → Links to CellChat MHC-II & MIF pathways"
- **Size**: 5x5 inches

### Panel 5C — Cross-Subtype Comparison Barplot
- **Data**: Combine both KO results, select top 8 genes by |Z| from each KO
- **Type**: Horizontal barplot showing Z-scores
- **X-axis**: Z-score (negative = downregulated, positive = upregulated)
- **Y-axis**: Gene names
- **Two panels side-by-side**:
  - Left: CXCL13 KO in CCR8+ (top 8 by |Z|)
  - Right: B2M KO in MKI67+ (top 10 by |Z|)
- **Color bars by**: 
  - BH significant: `#E74C3C` (red)
  - Non-significant: `#95A5A6` (gray)
- **Key highlight**: In B2M KO panel, mark HLA-DRB5 and CD74 with different color (`#2980B9`, blue) or annotation to emphasize CellChat link
- **Title per panel**: "CXCL13 KO → CCR8+ Treg" / "B2M KO → MKI67+ Treg"
- **Size**: 12x6 inches (combined)

### Panel 5D — Gene-Pathway Bridge: Key B2M KO Targets → Pathway Mapping
- **Data**: Hardcoded from B2M KO significant gene list (4 key targets)
- **Type**: Horizontal bar chart showing Z-scores, pathway annotations, and p-values
- **X-axis**: Z-score (all downregulated)
- **Y-axis**: Gene names (left-aligned text) — HLA-DRB5, CD74, HAVCR1, PTMS
- **Bars**: Horizontal, colored by pathway:
  - HLA-DRB5 (blue `#2980B9`): MHC-II Antigen Presentation → links to Fig 2C/D responder axis
  - CD74 (purple `#8E44AD`): MIF Signaling (MIF receptor) → links to Fig 2E non-responder axis
  - HAVCR1 (orange `#E67E22`): Immune Checkpoint → independent
  - PTMS (green `#27AE60`): MHC-I Peptide Loading → independent
- **Annotations**: Z-score value at bar tip, p.adj below each bar, pathway label on right
- **Threshold line**: Dotted line at |Z|=2 (significance threshold)
- **Legend**: 4 pathway colors with "(Fig. 2C/D)" and "(Fig. 2E)" cross-references
- **Title**: "B2M KO in MKI67+ Treg: Predicted Pathway Targets"
- **Size**: 6x3.5 inches (compact, fits next to 5C in 2×2 grid)
- **Rationale**: Replaces old GO BP dotplot. The 19-gene B2M KO set is too small for meaningful GO enrichment (only 12 Entrez IDs, driven by 2-3 genes). A direct gene→pathway mapping is more honest and more directly connects to CellChat narrative.

## Layout
- **Grid**: 2×2
  - Row 1: 5A | 5B (side by side, equal width)
  - Row 2: 5C | 5D (5C wider, 5D narrower)
- **Label panels**: A, B, C, D in top-left of each panel (bold, 14pt)
- **Figure title**: None (will add in final assembly)

## Output Files
Save to `D:/Research/cucumber/fig5/`:
- `fig5a_cxcl13_volcano.png` / `.pdf`
- `fig5b_b2m_volcano.png` / `.pdf`
- `fig5c_barplot_comparison.png` / `.pdf`
- `fig5d_go_dotplot.png` / `.pdf`
- `figure5_all.py` (combined Python script)

## Python Environment
- Use: `D:/Research/tomato/venv_scenic/Scripts/python.exe`
- Packages: pandas, numpy, matplotlib, seaborn, scipy
- Working dir: `D:/Research/cucumber/fig5/`

## Style Requirements
- Arial font throughout
- 300 DPI output
- PDF + PNG for each panel
- Consistent axis font sizes across panels 5A/5B
- Color-blind friendly (avoid red-green where possible; use red-blue-gray)
- Statistical text: "p.adj < 0.05 (BH)"
- No grid lines
- All panels white background
