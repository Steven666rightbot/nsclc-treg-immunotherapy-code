# Figure 3 — Results Text & Figure Legends

## Results Section (for manuscript)

### Transcriptional Programs Distinguish CCR8+ and MKI67+ Treg Subsets

To investigate the transcriptional programs driving the functional divergence of the two ML-discovered Treg subsets, we performed transcription factor (TF) activity analysis using DoRothEA regulon enrichment on CCR8+ and MKI67+ Treg cells separately, comparing responders (MPR/pCR) versus non-responders (non-MPR).

Within CCR8+ Treg cells, the top TFs with significantly higher activity in non-responders included MAX (effect = +0.010, p = 6.1 × 10⁻⁴⁵), ZBTB33 (effect = +0.054, p = 7.3 × 10⁻⁴⁵), TCF4 (effect = +0.043, p = 1.5 × 10⁻⁴²), and HSF1 (effect = −0.720, p = 2.3 × 10⁻⁴¹). Conversely, TFs enriched in responder CCR8+ Treg included RFX5 (effect = +0.397, p = 1.8 × 10⁻⁷), a master regulator of MHC-II antigen presentation, and NFYB (effect = +0.162, p = 4.2 × 10⁻³⁸), a component of the NF-Y complex (Figure 3A, left panel).

Within MKI67+ Treg cells, heat shock transcription factor 1 (HSF1) emerged as the most differentially activated TF, showing substantially higher activity in non-responders (effect = −0.888, p = 2.8 × 10⁻¹⁷). Additional TFs upregulated in non-responder MKI67+ Treg included STAT2 (effect = −0.275, p = 3.0 × 10⁻⁷), ELK1 (effect = −0.281, p = 1.4 × 10⁻⁵), and MAF (effect = −0.300, p = 5.7 × 10⁻⁶). In contrast, TFs enriched in responder MKI67+ Treg included RFX5 (effect = +0.350, p = 0.024), NR5A1 (effect = +0.256, p = 6.5 × 10⁻⁵), and CTCF (effect = +0.150, p = 2.7 × 10⁻⁴) (Figure 3A, right panel).

Notably, HSF1 was the top-ranked TF by absolute effect size in both subtypes but was particularly prominent in MKI67+ Treg (effect = −0.888 vs −0.720 in CCR8+), and its elevated activity in non-responders suggested a stress-responsive transcriptional program. To explore this, we examined the expression of known HSF1 target genes (29 heat shock protein genes) specifically within MKI67+ Treg cells (Figure 3B). Hierarchical clustering revealed a coordinated upregulation of multiple HSP genes including HSPA1A, HSPA1B, HSP90AA1, DNAJB1, and BAG3 in non-responder MKI67+ Treg cells, consistent with a heightened cellular stress response in this subset.

To integrate these findings into a mechanistic model, we constructed a regulatory network delineating the relationships between the two Treg subtypes, their key transcriptional regulators, and downstream effector pathways (Figure 3C). The tri-axis model highlights: (i) HSF1-driven stress/inflammatory program (HSP genes) predominantly active in MKI67+ Treg from non-responders; (ii) RFX5-dependent MHC-II antigen presentation preferentially active in CCR8+ Treg from responders; and (iii) BCL6/MAX-mediated transcriptional programs linking Treg plasticity to the MIF/CD99 inflammatory axis.

We further validated the differential activity of three key TFs — HSF1, BCL6, and RFX5 — across both subtypes and response groups using DoRothEA activity scores (Figure 3D). HSF1 activity was significantly elevated in non-responders for both CCR8+ (p < 0.001) and MKI67+ (p < 0.001) Treg. BCL6 showed elevated activity in non-responder MKI67+ Treg (p < 0.01), consistent with its role in T follicular helper-like programming. RFX5, a key regulator of MHC-II genes, showed significantly higher activity in responder CCR8+ Treg (p < 0.001), suggesting that antigen presentation competency is a hallmark of responder-associated Treg function.

Collectively, these results reveal that the two Treg subtypes are governed by distinct transcriptional programs: CCR8+ Treg from responders are characterized by MHC-II antigen presentation competence (RFX5), while MKI67+ Treg from non-responders are dominated by a stress-responsive, HSF1-driven heat shock program. This transcriptional divergence provides a mechanistic basis for the differential cell–cell communication patterns observed in Figure 2.

---

## Figure Legend

### Figure 3. Transcriptional Regulation of CCR8+ and MKI67+ Treg Subsets.

**(A)** Top differentially active transcription factors (TFs) between responder (MPR/pCR) and non-responder (non-MPR) within CCR8+ (left) and MKI67+ (right) Treg cells, ranked by absolute effect size (Responder median − Non-responder median). TFs were scored using DoRothEA regulon enrichment. Bar colors indicate direction of activity: red (higher in non-responder) and blue (higher in responder). Asterisks denote Benjamini-Hochberg adjusted p-values (* p < 0.05, ** p < 0.01, *** p < 0.001).

**(B)** Heatmap showing expression of 29 HSF1 target genes (heat shock protein genes) across MKI67+ Treg cells, split by response group (Responder, n = 466; Non-responder, n = 1,050). Expression values are Z-score normalized per gene. Rows are hierarchically clustered (complete linkage, Euclidean distance); columns are ordered by response group. Top annotation bar: green = Responder, purple = Non-responder.

**(C)** Regulatory network schematic integrating the tri-axis model. Left: CCR8+ and MKI67+ Treg subtypes. Center: key transcriptional hubs (HSF1, BCL6, RFX5, MAX/MYC). Right: downstream effector pathways (MHC-I, MHC-II, MIF/CD99, heat shock response). Arrow thickness indicates strength of association; dashed lines represent indirect or inferred regulation. NR ↑ / R ↓ annotations indicate the response group with higher activity.

**(D)** Violin plots of DoRothEA activity scores for three key TFs (HSF1, BCL6, RFX5) across four groups: CCR8+ Responder, CCR8+ Non-responder, MKI67+ Responder, and MKI67+ Non-responder. Individual cell scores are shown as jittered points. Pairwise Mann-Whitney U test p-values: *** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant.
