# Figure 2 — Results Text & Figure Legends

## Results Section (for manuscript)

### Intercellular Communication Rewiring Reveals MHC‑I/II and MIF/CD99 Axes Targeting Distinct Treg Subsets

To investigate how CCR8⁺ and MKI67⁺ Tregs differentially engage in intercellular communication between responders and non-responders, we performed CellChat analysis on stratified ~50K cells from each response group (GSE243013). Global network visualization revealed that the responder (Strong) group was dominated by MHC‑I and MHC‑II signaling, while the non-responder (Weak) group showed prominent MIF and CD99 signaling, with distinct network topologies and interaction strengths (Figure 2A).

To systematically quantify these differences, we constructed a focused 12 × 12 differential communication heatmap comparing signaling probabilities between the two Treg subtypes (CCR8⁺, MKI67⁺) and ten partner cell-type groups (Figure 2B). This analysis revealed striking subtype-specific patterns: ILC↔Treg_CCR8 showed the most extreme negative differential in non-responders versus responders (log₂FC = −14.4), suggesting communication between these two cell types is almost entirely absent in non-responders. Conversely, the γδT cell column was uniformly red-shifted, indicating broadly stronger communication from γδT cells to both Treg subtypes in non-responders. Within the Treg-to-partner axis, Treg_MKI67→NK showed the strongest positive differential (log₂FC = +5.7), suggesting enhanced MKI67⁺ Treg-to-NK signaling in non-responders.

Analysis of incoming signaling pathways to each Treg subtype revealed a striking dichotomy (Figure 2C). In responders, MHC‑I pathway signaling (HLA‑A/B/C → CD8A/CD8B) was specifically directed toward CCR8⁺ Tregs, with a total communication probability of 5.8 × 10⁻⁵ compared to essentially zero in non-responders (Figure 2D, left panel). MHC‑II signaling (HLA‑DRB1/DQA1/DQB1/DMA → CD4) also preferentially targeted CCR8⁺ Tregs in responders, showing ~4-fold higher probability than in MKI67⁺ Tregs. In contrast, MIF signaling (MIF → CD74/CD44) was predominantly directed toward MKI67⁺ Tregs and was uniquely active in non-responders (probability 2.1 × 10⁻⁴ vs. 0 in responders; Figure 2E, left panel). CD99 signaling (CD99 → CD99) showed a similar pattern, with non-responder MKI67⁺ Tregs receiving stronger CD99 signals.

Detailed ligand-receptor pair analysis confirmed that multiple MHC‑I pairs (HLA‑A/B/C → CD8A, HLA‑A/B/C → CD8B) and MHC‑II pairs (HLA‑DRB1 → CD4, HLA‑DMA → CD4) were uniquely active in the responder group, with CCR8⁺ Tregs consistently showing 5–6 fold higher probabilities than MKI67⁺ Tregs (Figure 2D). The MIF axis in non-responders was mediated through MIF → CD74 and MIF → CD44 pairs, with MKI67⁺ Tregs showing a ~2.5-fold higher MIF-receiving probability than CCR8⁺ Tregs (Figure 2E).

These results establish that CCR8⁺ and MKI67⁺ Tregs are engaged by fundamentally different intercellular signaling networks: MHC‑I/II antigen presentation signaling operates preferentially in the responder microenvironment and targets CCR8⁺ Tregs, while MIF/CD99 inflammatory signaling operates preferentially in the non-responder microenvironment and targets MKI67⁺ Tregs. This differential communication wiring provides a mechanistic framework linking specific Treg subsets to immunotherapy response status.

---

## Figure Legend

### Figure 2. CellChat Analysis Reveals Treg Subtype-Specific Intercellular Communication Rewiring Between Responders and Non-Responders.

**(A)** Global intercellular communication networks for responder (Strong, MPR/pCR) and non-responder (Weak, non-MPR) groups, visualized as circle plots. Left: Strong group network dominated by MHC‑I (red) and MHC‑II (green) signaling. Right: Weak group network dominated by MIF (purple) and CD99 (orange) signaling. Edge width represents communication probability; node size reflects cell number in each group. Line colors correspond to distinct signaling pathways (visualized using CellChat's hierarchical pathway clustering).

**(B)** Focused 12 × 12 differential communication heatmap showing log₂ fold change (FC) of communication probability (Non-responder vs. Responder) between two Treg subtypes (CCR8⁺, MKI67⁺) and ten partner cell-type groups (B cell, DC, Endothelial, Epithelial, Fibroblast, γδT, ILC, Macrophage, Mast, NK). Color scale: red = higher in non-responders; blue = higher in responders. Key extremes: ILC↔Treg_CCR8 (log₂FC = −14.4), Treg_MKI67→NK (log₂FC = +5.7).

**(C)** Incoming signaling pathways to CCR8⁺ Treg (left) and MKI67⁺ Treg (right), ranked by total communication probability across both response groups. Bar length represents the sum of communication probabilities for each pathway. Pathways are ordered by combined rank (mean rank across both subtypes). MHC‑I dominates CCR8⁺ Treg incoming signaling in responders; MIF dominates MKI67⁺ Treg incoming signaling in non-responders.

**(D)** MHC‑I (HLA‑A/B/C → CD8A/CD8B) and MHC‑II (HLA‑DRB1/DQA1/DQB1/DMA → CD4) ligand-receptor pair probabilities for CCR8⁺ Treg (blue) and MKI67⁺ Treg (rose) in responders (left) and non-responders (right). CCR8⁺ Treg shows ~5–6 fold higher MHC‑I signaling probability than MKI67⁺ Treg in responders. MHC‑II signaling (HLA‑DRB1 → CD4) is exclusively detected in responder CCR8⁺ Treg.

**(E)** MIF (MIF → CD74/CD44) and CD99 (CD99 → CD99) ligand-receptor pair probabilities. MIF signaling is uniquely active in non-responders and preferentially targets MKI67⁺ Treg (~2.5-fold higher than CCR8⁺ Treg). CD99 signaling shows a concordant pattern with stronger non-responder signaling to MKI67⁺ Treg.
