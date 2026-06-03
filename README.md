# NSCLC Treg Immunotherapy Resistance — Code Repository

Companion code for the manuscript:
> **Integrative Single-Cell Dissection Reveals a CCR8-to-MKI67 Treg Continuum Driving Immunotherapy Resistance in NSCLC**

**Repository**: https://github.com/Steven666rightbot/nsclc-treg-immunotherapy-code

This repository contains all computational analysis and figure-generation scripts used in the study.

---

## Repository Structure

```
.
├── figures/                          # Publication-quality figure generation
│   ├── fig1/                         # ML discovery (XGBoost + SHAP)
│   ├── fig2/                         # Pseudotime trajectory (Slingshot/Monocle3)
│   ├── fig3/                         # CellChat intercellular communication
│   ├── fig4/                         # TF regulation (DoRothEA/pySCENIC)
│   ├── fig5/                         # Spatial transcriptomics (Visium)
│   ├── fig6/                         # In silico knockout (scTenifoldKnk)
│   ├── supplementary/                # Supplementary figures
│   └── _global_config.py             # Shared colors, fonts, styling
│
├── analysis/                         # Core computational analyses
│   ├── 01_data_preprocessing/        # Dataset curation & cell annotation
│   ├── 02_ml_discovery/              # XGBoost, SHAP, external validation
│   ├── 03_pseudotime_trajectory/     # Slingshot & Monocle3 trajectory inference
│   ├── 04_cellchat_communication/    # Cell-cell communication analysis
│   ├── 05_tf_regulation/             # DoRothEA/pySCENIC TF activity inference
│   ├── 06_spatial_transcriptomics/   # Visium spatial analysis
│   ├── 07_in_silico_knockout/        # scTenifoldKnk virtual KO & GO enrichment
│   └── 08_nichenet_validation/       # NicheNet ligand-target validation
│
└── README.md
```

---

## Datasets Used

| Dataset | GEO ID | Purpose |
|---------|--------|---------|
| Discovery cohort | [GSE243013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013) | 242 post-neoadjuvant chemo-immunotherapy NSCLC tumor samples (Liu et al., *Cell* 2025); 1,749,316 immune cells after quality control |
| External validation (NSCLC) | [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) | 14 NSCLC patients treated with neoadjuvant anti-PD-1 ± chemotherapy; 12 evaluable post-treatment samples used for validation |
| Spatial transcriptomics | [GSE189487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189487) | 10x Visium (13 tissue sections from 8 NSCLC patients); 10 NSCLC sections (7 low-grade LUAD, 3 high-grade NSCLC; 8,686 spots) retained for analysis after excluding PDAC samples |
| Cross-cancer validation (melanoma) | [GSE120575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575) | 48 tumor biopsies from 32 melanoma patients (19 pre-treatment, 29 post-treatment; anti-PD-1 ± CTLA-4) |

---

## Key Analyses & Corresponding Scripts

### Figure 1 — Machine Learning Discovery
- **XGBoost classifier + SHAP**: `analysis/02_ml_discovery/xgboost_shap_analysis.py`
- **Decision tree**: `analysis/02_ml_discovery/decision_tree_classifier.py`
- **External validation**: `analysis/02_ml_discovery/gse207422_six_feature_validation.py`
- **Figure generation**: `figures/fig1/`

### Figure 2 — Pseudotime Trajectory
- **Slingshot (R)**: `analysis/03_pseudotime_trajectory/slingshot_ccr8_mki67_treg.R`
- **Monocle3 (R)**: `analysis/03_pseudotime_trajectory/treg_monocle3_ccr8_mki67_only.R`
- **Figure generation**: `figures/fig2/`

### Figure 3 — CellChat Communication
- **CellChat input preparation**: `analysis/01_data_preprocessing/prepare_gse243013_cellchat_input*.py`
- **CellChat analysis (R)**: `analysis/04_cellchat_communication/cellchat_gse243013_strong_weak.R`
- **Virtual knockout (R)**: `analysis/04_cellchat_communication/cellchat_virtual_knockout*.R`
- **Figure generation**: `figures/fig3/`

### Figure 4 — Transcriptional Regulation
- **DoRothEA/pySCENIC**: `analysis/05_tf_regulation/dorothea_tf_analysis.py`
- **TF differential analysis**: `analysis/05_tf_regulation/dorothea_subtypes_and_shap.py`
- **HSF1 heatmap**: `analysis/05_tf_regulation/fig3b_hsf1_heatmap.py`
- **Regulatory network**: `analysis/05_tf_regulation/fig3c_regulatory_network.py`
- **Figure generation**: `figures/fig4/`

### Figure 5 — Spatial Transcriptomics
- **Visium analysis**: `analysis/06_spatial_transcriptomics/gse189487_spatial_analysis.py`
- **Spatial deconvolution**: `analysis/06_spatial_transcriptomics/spatial_treg_deconvolution.py`
- **Figure generation**: `figures/fig5/`

### Figure 6 — In Silico Knockout
- **scTenifoldKnk (R)**: `analysis/07_in_silico_knockout/scTenifoldKnk_merged_treg_ko*.R`
- **KO result analysis**: `analysis/07_in_silico_knockout/analyze_ko_results.py`
- **GO enrichment (R)**: `analysis/07_in_silico_knockout/run_go.R`
- **Figure generation**: `figures/fig6/`

---

## Environment & Dependencies

**Python** (≥3.8)
```
scanpy, anndata, numpy, pandas, matplotlib, seaborn
scikit-learn, xgboost, shap, scipy, statsmodels
loompy, pillow, python-docx
```

**R** (≥4.0)
```
Seurat, Monocle3, Slingshot, CellChat, clusterProfiler
SCENIC/AUCell, NicheNet, scTenifoldKnk
```

---

## Citation

If you use this code, please cite the corresponding manuscript.

---

## Contact

Corresponding author: Zhuchen Yan (zhyan@tmu.edu.cn)
