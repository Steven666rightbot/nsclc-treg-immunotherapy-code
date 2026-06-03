# Figure 4 — Spatial Validation of Treg Subtype Distribution and Grade-Dependent Polarization

## Results Text

To investigate how the CCR8⁺→MKI67⁺ Treg differentiation continuum (Fig. 2) manifests in tissue space, we analyzed Visium spatial transcriptomics data from GSE189487, comprising 8,686 spots across 10 NSCLC tissue sections (7 low-grade LUAD, 3 high-grade NSCLC). We assessed the spatial distribution of CCR8⁺ and MKI67⁺ Treg signature scores and computed an MHC-I antigen presentation module score (HLA-A, HLA-B, HLA-C, B2M, TAP1).

CCR8⁺ Treg signature scores exhibited robust stromal enrichment across all NSCLC samples (Fig. 4A), with the strongest demarcation observed in low-grade LUAD. MKI67⁺ Treg scores showed a more nuanced spatial pattern: significant stromal enrichment in low-grade tumors (p < 0.001), but this polarization was selectively lost in high-grade NSCLC (p = 0.87; Fig. 4B–C, E). Quantitatively, the stromal enrichment score (stroma z-score minus tumor core z-score) for CCR8⁺ Tregs remained high in both low-grade (Δ = 0.86) and high-grade (Δ = 0.40) tumors, whereas MKI67⁺ Treg enrichment dropped from Δ = 0.37 in low-grade to Δ = 0.08 (non-significant) in high-grade tumors (Fig. 4E).

These grade-dependent spatial patterns suggest that the terminal differentiation from CCR8⁺ to MKI67⁺ Tregs — inferred from pseudotime trajectory analysis (Fig. 2) — requires an intact stromal microenvironment, which is progressively disrupted in poorly differentiated tumors. MHC-I module scores showed a heterogeneous spatial pattern enriched in tumor cores (Fig. 4D), consistent with the established role of MHC-I in tumor cell antigen presentation. Notably, MHC-I and CCR8⁺ Treg signatures did not exhibit spot-level co-expression (Spearman ρ = −0.12), reflecting their occupancy of distinct tissue compartments (tumor core vs. stroma) rather than simple spatial colocalization.

---

## Figure Legends

**Figure 4. Spatial validation of Treg subtype distribution and grade-dependent polarization.**

**(A–B)** Spatial feature plots showing CCR8⁺ Treg signature scores (A) and MKI67⁺ Treg signature scores (B) across four representative NSCLC tissue sections: LG_2 and LG_6 (low-grade LUAD), and HG_3 (high-grade NSCLC). Spot colors indicate z-score-normalized signature scores (red = high, blue = low). Dashed outlines denote convex hull boundaries of hard stroma (red) and tumor core (blue) regions.

**(C)** Grouped violin plots comparing CCR8⁺ and MKI67⁺ Treg z-scores between hard stroma and tumor core regions across all NSCLC samples. Box plots (inner) show median and interquartile range. Statistical significance was assessed by Mann–Whitney U test (*** p < 0.001).

**(D)** Spatial feature plots of MHC-I module score (mean normalized expression of HLA-A, HLA-B, HLA-C, B2M, TAP1, z-scored) across the same tissue sections as in (A–B). Color scale as in (A–B).

**(E)** Grade-dependent stromal enrichment of CCR8⁺ and MKI67⁺ Treg signatures. Bar heights represent the difference in mean z-score between hard stroma and tumor core regions (Δ = stroma − tumor). CCR8⁺ Tregs maintained significant stromal enrichment in both low-grade and high-grade NSCLC, whereas MKI67⁺ Treg stromal polarization was significant only in low-grade tumors and lost in high-grade disease. Significance markers indicate stroma vs. tumor core comparisons: *** p < 0.001; ns, not significant.

---

## Key Statistics

| Comparison | Metric | Value | p-value |
|-----------|--------|-------|---------|
| CCR8⁺ stroma vs tumor (Low Grade) | Δ (z-score) | 0.864 | 9.24 × 10⁻¹¹² |
| CCR8⁺ stroma vs tumor (High Grade) | Δ (z-score) | 0.404 | 2.37 × 10⁻⁷ |
| MKI67⁺ stroma vs tumor (Low Grade) | Δ (z-score) | 0.375 | 1.80 × 10⁻²¹ |
| MKI67⁺ stroma vs tumor (High Grade) | Δ (z-score) | 0.081 | 0.869 (ns) |
| CCR8⁺ Low vs High Grade (stroma) | mean diff | 0.151 | 0.022 |
| MKI67⁺ Low vs High Grade (stroma) | mean diff | 0.048 | 0.078 |
| MHC-I vs CCR8⁺ (all spots) | Spearman ρ | −0.120 | 1.79 × 10⁻⁴⁰ |
| MHC-I vs MKI67⁺ (all spots) | Spearman ρ | −0.074 | 3.08 × 10⁻¹⁶ |
| CCR8⁺ vs MKI67⁺ overlap (Jaccard) | mean ± SD | 0.127 ± 0.031 | — |

---

## Narrative Integration

This figure serves three narrative functions in the manuscript:

1. **Spatial validation of the differentiation continuum**: CCR8⁺ and MKI67⁺ Tregs occupy distinct but related spatial niches — CCR8⁺ as a stromal-resident "gatekeeper" and MKI67⁺ as a conditionally polarized "proliferator." Their spatial mutual exclusion (Jaccard ~0.13) supports the pseudotime-inferred continuum rather than independent lineages.

2. **Grade-dependent microenvironment quality**: The selective loss of MKI67⁺ stromal polarization in high-grade tumors suggests that terminal Treg differentiation requires an intact stromal niche. As tumors dedifferentiate, the microenvironment becomes unable to support the CCR8⁺→MKI67⁺ transition, leaving CCR8⁺ Tregs as the dominant suppressive population.

3. **MHC-I as a cross-compartment signal**: MHC-I enrichment in tumor cores (Fig. 4D) and CCR8⁺ enrichment in stroma (Fig. 4A) are spatially segregated, consistent with CellChat's prediction of MHC-I→CCR8⁺ Treg communication as a **functional axis** (across tissue regions) rather than a **colocalization** relationship (within single spots).
