# Figure 5 — Causal Perturbation Analysis (scTenifoldKnk Virtual KO)

## Results Text (Draft)

To validate the functional relevance of the CellChat-predicted communication axes (Fig. 2), we performed in-silico gene knockout experiments using scTenifoldKnk on CCR8+ and MKI67+ Treg cells from GSE243013. We targeted two key genes bridging the communication and regulatory landscapes: CXCL13 (a TLS-associated chemokine highly expressed in CCR8+ Treg) and B2M (the MHC-I light chain, central to the MHC-I pathway targeting CCR8+ Treg in responders).

Virtual knockout of CXCL13 in CCR8+ Treg (Fig. 5A) identified 8 significantly perturbed genes (BH-adjusted p < 0.05), including CLDND2 (Z = 4.10), AL136018.1 (Z = 3.96), and AC034238.2 (Z = 3.18). CCL3 (MIP-1α, Z = 2.69) showed a trend toward downregulation below the BH threshold, suggesting a potential chemokine network perturbation.

Virtual knockout of B2M in MKI67+ Treg (Fig. 5B) revealed broader perturbation, with 19 significantly perturbed genes. Among these, **HLA-DRB5** (MHC class II, Z = 2.30, p.adj = 0.0006) and **CD74** (MIF receptor, Z = 2.11, p.adj = 0.045) directly link MHC-I disruption to MHC-II and MIF pathway modulation, connecting back to the CellChat findings. Additional perturbed genes included HAVCR1 (TIM-1, an immune checkpoint) and PTMS (involved in MHC-I peptide loading) (Fig. 5C).

To map these perturbations to functional pathways (Fig. 5D), we annotated each key gene by its primary effector pathway. HLA-DRB5 feeds into MHC-II antigen presentation — the axis selectively active in responder CCR8+ Tregs (Fig. 2C/D). CD74 encodes the MIF receptor, central to the non-responder MKI67+ Treg signaling identified by CellChat (Fig. 2E). HAVCR1 (TIM-1) represents an immune checkpoint receptor, and PTMS participates in MHC-I peptide loading. This gene-pathway bridge directly connects the B2M KO perturbation results to the tri-axis model, suggesting that MHC-I disruption in MKI67+ Tregs may simultaneously affect MHC-II and MIF signaling components.

## Figure Legends

**Figure 5. Causal perturbation analysis of key signaling axes via scTenifoldKnk virtual knockout.**

**(A–B)** Volcano plots showing Z-scores versus -log₁₀(adjusted p-value) for differentially regulated genes following virtual knockout of CXCL13 in CCR8+ Treg (A) and B2M in MKI67+ Treg (B). Red points indicate significantly perturbed genes (BH-adjusted p < 0.05, |Z| > 2). Gray points show genes with |Z| > 2 but not meeting BH significance. Select significantly perturbed genes are labeled.

**(C)** Horizontal bar plots comparing top perturbed genes by |Z| score in CXCL13 KO in CCR8+ Treg (left) and B2M KO in MKI67+ Treg (right). Bar color indicates BH significance (red: p.adj < 0.05; gray: not significant). Blue-highlighted genes (HLA-DRB5, CD74) connect B2M KO effects to MHC-II antigen presentation and MIF signaling pathways identified in Figure 2.

**(D)** Gene-pathway bridge plot. Horizontal bars show Z-scores of four key genes perturbed by B2M KO in MKI67+ Treg, colored by their primary pathway: HLA-DRB5 → MHC-II antigen presentation (blue, linking to Fig. 2C/D responder axis), CD74 → MIF signaling (purple, linking to Fig. 2E non-responder axis), HAVCR1 → immune checkpoint (orange), PTMS → MHC-I peptide loading (green). Pathway labels on the right indicate the effector pathway for each gene.

## Data Sources
- CXCL13 KO results: `ko_CCR8_CXCL13.csv` (2,029 genes, 8 BH-significant)
- B2M KO results: `ko_MKI67_B2M_lowres.csv` (2,013 genes, 19 BH-significant)
- GO BP enrichment (archived): `go_bp_b2m_ko_top20.csv` (deprecated — replaced by gene-pathway bridge)
