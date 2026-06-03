# Figure 6 — Slingshot Pseudotime Trajectory Analysis

## Results Text (Draft)

To investigate the developmental relationship between the two machine-learning-discovered Treg subsets, we performed pseudotime trajectory analysis using Slingshot on 26,725 CCR8+ and MKI67+ Treg cells from GSE243013. Cells were clustered into 10 transcriptional states based on PCA-reduced expression, and a minimum spanning tree was constructed to define the principal trajectory path (Fig. 6A).

The Slingshot trajectory revealed a continuous transition from CCR8+ Treg to MKI67+ Treg along a single major lineage (Lineage 1, 7,236 cells). CCR8+ Treg cells were enriched at low pseudotime values (median pseudotime = 6.28), while MKI67+ Treg cells occupied high pseudotime values (median = 14.26, p < 2.2 × 10⁻¹⁶, Wilcoxon rank-sum test), consistent with a progressive activation and proliferation program (Fig. 6B).

Examination of pseudotime distributions stratified by response group (Fig. 6C) revealed that non-responders had a higher density of MKI67+ Treg cells at the terminal pseudotime range (pseudotime > 12), whereas responder MKI67+ Treg cells were more evenly distributed. This suggests that in non-responsive tumors, Treg cells more efficiently traverse the CCR8→MKI67 differentiation trajectory toward a terminally activated state.

Analysis of key gene expression along pseudotime (Fig. 6D) confirmed that CCR8 expression peaked at intermediate pseudotime and declined at the terminal state, while MKI67 expression monotonically increased along the trajectory. Notably, CTLA4, a key immunosuppressive checkpoint receptor, showed sustained elevation throughout the trajectory in non-responders compared to responders. CD74 (the MIF receptor identified in Figure 2) showed an early peak in CCR8+ Treg that was more pronounced in non-responders, consistent with the enhanced MIF signaling in non-responsive tumors.

These findings establish that CCR8+ and MKI67+ Treg represent sequential states along a differentiation continuum, with non-responder tumors exhibiting more complete progression toward the proliferative, immunosuppressive MKI67+ terminal state.

## Figure Legends

**Figure 6. Pseudotime trajectory analysis reveals CCR8+ Treg → MKI67+ Treg differentiation continuum.**

**(A)** UMAP visualization of 26,725 CCR8+ and MKI67+ Treg cells from GSE243013. Cells are colored by subtype (blue: CCR8+ Treg; rose: MKI67+ Treg). The black line indicates the Slingshot principal curve (Lineage 1), with an arrow indicating the direction of pseudotime progression. The arrow indicates the direction of increasing pseudotime from CCR8+ → MKI67+.

**(B)** Same UMAP embedding colored by Slingshot pseudotime value (plasma colormap, range 0–16), demonstrating the continuous transition from CCR8+ (low pseudotime) to MKI67+ (high pseudotime).

**(C)** Kernel density estimates of pseudotime distribution for each subtype stratified by response group. Solid dark lines: non-responder; lighter lines: responder. Non-responder MKI67+ Treg show higher density at terminal pseudotime values (>12).

**(D)** Binned mean expression (± SEM) of CCR8, MKI67, CTLA4, and CD74 along pseudotime, stratified by response group (green: responder; purple: non-responder). Expression values are log₁₀(CPM + 1) normalized.

## Key Data

- **Total cells**: 26,725 (CCR8+: 24,431; MKI67+: 2,294)
- **Main trajectory**: Lineage 1, clusters 8→1→2→4→5→6→10
- **Pseudotime**: CCR8+ median = 6.28, MKI67+ median = 14.26
- **Analysis**: Slingshot on Monocle3-preprocessed UMAP + clusters
