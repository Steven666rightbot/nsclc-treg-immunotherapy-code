# NSCLC Treg Immunotherapy Resistance — Code Repository

Companion code for the manuscript:
> **Integrative Single-Cell Dissection Reveals a CCR8-to-MKI67 Treg Continuum Driving Immunotherapy Resistance in NSCLC**

**Journal**: Journal of Translational Medicine (submitted)

**Repository**: https://github.com/yan.xiaohan/nsclc-treg-immunotherapy-code

---

## Repository Structure

| Directory | Contents |
|-----------|----------|
| `fig1/` | Machine learning discovery (XGBoost + SHAP) |
| `fig2_revised_merged/` | CellChat intercellular communication |
| `fig3/` | TF regulation (pySCENIC / DoRothEA) |
| `fig4/` | Spatial transcriptomics (Visium) |
| `fig5/` | In silico KO (scTenifoldKnk) + NicheNet |
| `fig6/` | Slingshot pseudotime trajectory |
| `supp_fig1/` | NicheNet supplementary figure |
| `supp_fig2/` | Cross-cancer validation (melanoma) |
| `_global_config.py` | Shared color palette and matplotlib rcParams |

## Environment

- **Python**: 3.8+ (scanpy 1.9, numpy 1.21, pandas 1.3, scikit-learn 1.0, xgboost 1.6, matplotlib 3.5)
- **R**: 4.5.2 (Seurat 4.3, CellChat 2.1, Slingshot 2.4, clusterProfiler 4.6)

## Data Availability

Public datasets from GEO:
- GSE243013 (discovery cohort, Liu et al. Cell 2025)
- GSE207422 (NSCLC external validation)
- GSE189487 (Visium spatial transcriptomics)
- GSE120575 (melanoma cross-cancer validation)

## Reproducibility

All XGBoost models used fixed random seeds:
- Five-fold CV: `random_state = 42`
- External validation: `seed = 1`
- Seed sensitivity confirmed across 20 seeds (AUC range 0.719–0.844)

## Note

This repository contains analysis scripts and figure generation code. Raw single-cell data are available via GEO (see above). Some intermediate processed files (e.g., `.h5ad`, `.rds`) exceed GitHub file size limits and are not included; they can be regenerated from the provided scripts and public data.
