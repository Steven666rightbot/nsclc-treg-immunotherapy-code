# Figure 1: Clinical Prediction Framework — Documentation

**Project**: Tomato (Treg-Mechanoscore-Immunotherapy)  
**Date**: 2026-05-07  
**Figure format**: 300 DPI PNG + vector PDF  
**Font**: Arial (ASCII-safe, no Unicode superscripts/subscripts)

---

## Figure 1 Structure (8 panels)

| Panel | Title | Content | Script | Data Source |
|-------|-------|---------|--------|-------------|
| **a** | Study Design | Patient cohort flowchart (manual) | — | `data/GSE243013_metadata/` |
| **b** | CD4T Subtype Marker Heatmap | 9 subtypes × 27 markers, rounded-corner full matrix | `code/generate_fig1_panels.py` (Panel B) | `data/GSE243013_reference_pseudobulk.csv` |
| **c** | Treg Abundance by Response | CCR8+ / MKI67+ Treg proportions, Mann-Whitney U | `code/generate_fig1_panels.py` (Panel C) | `data/cell_proportions.csv`, `data/sample_labels.csv` |
| **d** | ROC Curves | 5-fold cross-validation AUC | `code/xgboost_shap_analysis.py` | `data/cell_proportions.csv`, `data/sample_labels.csv` |
| **e** | SHAP Bar Importance | Mean |SHAP| ranking (top 20) | `code/xgboost_shap_analysis.py` | XGBoost SHAP values |
| **f** | Decision Tree | Compact tree (max_depth=3, AUC=0.744) | `code/generate_fig1_panels.py` (Panel F) | `data/cell_proportions.csv`, `data/sample_labels.csv` |
| **g** | SHAP Beeswarm | Top 20 feature SHAP distributions | `code/xgboost_shap_analysis.py` | XGBoost SHAP values |
| **h** | External Validation (GSE207422) | Full XGBoost (51 subtypes) validation on GSE207422 post-treatment samples after marker-based reannotation (n=12). AUC=0.844, ACC=75%. **GSE207422 is used solely for independent external validation of the prediction model; all downstream mechanistic analyses (Figures 2-6) are based exclusively on GSE243013.** | `code/fig1h_external_validation_proper.py` | `data/GSE207422_cell_proportions_post.csv` |

---

## Output Files (all in `figures/demo/1/`)

```
fig1b_cd4t_marker_heatmap.png        # Panel b
fig1b_cd4t_marker_heatmap.pdf
fig1c_treg_abundance_with_pvalue.png # Panel c
fig1c_treg_abundance_with_pvalue.pdf
fig1d_ROC_curves_5fold.png           # Panel d
fig1e_SHAP_bar_importance.png        # Panel e
fig1f_decision_tree_compact.png      # Panel f
fig1f_decision_tree_compact.pdf
fig1g_SHAP_summary_beeswarm.png      # Panel g
fig1h_gse207422_external_validation_roc.png  # Panel h
fig1h_gse207422_external_validation_roc.pdf  # Panel h
```

---

## Unified Color Scheme (Green-White-Purple)

All Fig1 panels use a consistent diverging palette:

| Meaning | Hex | Usage |
|---------|-----|-------|
| High / Responder / Positive | `#8bc98b` | Boxplots (R), heatmap high, DT responder nodes, ROC mean |
| Low / Non-responder / Negative | `#b08ebf` | Boxplots (NR), heatmap low, DT non-responder nodes |
| Mid / Neutral | `#ffffff` | Heatmap mid, DT mixed nodes |
| **Diverging cmap** | `['#b08ebf','#d4c2df','#ffffff','#c2e4c2','#8bc98b']` | All continuous scales |

---

## Data Files Detail

### `data/GSE243013_reference_pseudobulk.csv`
- **Rows**: 51 cell subtypes (including 9 CD4T subtypes)
- **Columns**: 31,831 genes
- **Values**: Pseudobulk average expression per subtype
- **Used by**: Panel b

### `data/cell_proportions.csv`
- **Rows**: 242 samples
- **Columns**: 51 cell subtype proportions
- **Used by**: Panels c, d, e, f, g

### `data/sample_labels.csv`
- **Rows**: 242 samples
- **Columns**: `response` (1=Responder, 0=Non-responder)
- **Used by**: Panels c, d, e, f, g

---

## How to Regenerate

### 1. Generate panels b, c, f
```bash
cd D:\Research\tomato
venv_scenic\Scripts\python.exe code\generate_fig1_panels.py
```

### 2. Generate panels d, e, g (requires XGBoost model training)
```bash
cd D:\Research\tomato
venv_scenic\Scripts\python.exe code\xgboost_shap_analysis.py
```
*Note: `xgboost_shap_analysis.py` currently reads metadata from `D:\Research\土豆\数据\raw\GSE243013\`. Update path if project migrated to `tomato`.*

---

## Changelog

### 2026-05-07
- **Added** Panel b: CD4T subtype-specific marker heatmap (rounded-corner, full 9×27 matrix, green-white-purple diverging)
- **Renumbered**: old fig1b → fig1c; old fig1e → fig1f
- **Unified colors**: All Fig1 panels now use soft green (`#8bc98b`) / soft purple (`#b08ebf`) instead of red/blue
- **Updated**: `xgboost_shap_analysis.py` output filenames prefixed with `fig1d/e/g`
- **Fixed**: Decision Tree node colors mapped from sklearn default orange/blue to green/purple via post-hoc patch recoloring
- **Added** Panel h: External validation using GSE207422 (n=14 samples) with corrected annotation. LOOCV ROC-AUC = 0.700, Accuracy = 85.7%. **Note: GSE207422 is used solely for independent external validation of the predictive model in Panel h. All downstream mechanistic analyses (CellChat, pySCENIC, KO, Visium) are performed exclusively on GSE243013.**

---

## Notes for Editors

1. **Font compatibility**: Arial lacks Unicode superscript/subscript glyphs. All p-values use ASCII format (e.g., `p = 6.76e-16`, `p < 1e-10`).
2. **Panel a** (Study Design) is currently a manual diagram; recommend redrawing in Illustrator/PowerPoint for final submission.
3. **Figure dimensions**: All panels exported at 300 DPI. Recommended composite width for 2-column journal: 180 mm.
4. **GSE207422 scope limitation**: GSE207422 appears only in Figure 1 Panel h (external validation) and is not used in any downstream mechanistic analysis. All mechanistic findings in Figures 2-6 derive from GSE243013.
