# GSEA for MKI67 only
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

INPUT_DIR <- "D:/Research/tomato/vB_revised_integrin_fyn_axis"
OUTPUT_DIR <- "D:/Research/tomato/go_exploration/output"

df <- read.csv(file.path(INPUT_DIR, "de_mki67_responder_vs_nonresponder.csv"))

# Rank genes: sign(logFC) * -log10(pval)
df$rank_score <- sign(df$logfoldchanges) * (-log10(df$pvals + 1e-300))
df <- df[order(df$rank_score, decreasing = TRUE), ]

geneList <- df$rank_score
names(geneList) <- df$names

cat(sprintf("Running GSEA for MKI67 (%d genes)...\n", length(geneList)))

gsea <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "SYMBOL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
  out_df <- as.data.frame(gsea)
  cat(sprintf("%d significant gene sets\n", nrow(out_df)))
  write.csv(out_df, file.path(OUTPUT_DIR, "gsea_mki67.csv"), row.names = FALSE)
  
  cat("\nTop UP in R (positive NES):\n")
  up <- out_df[out_df$NES > 0, ]
  if (nrow(up) > 0) {
    print(up[1:min(8, nrow(up)), c("Description", "NES", "pvalue", "qvalue")], row.names = FALSE)
  }
  
  cat("\nTop UP in NR (negative NES):\n")
  down <- out_df[out_df$NES < 0, ]
  if (nrow(down) > 0) {
    print(down[1:min(8, nrow(down)), c("Description", "NES", "pvalue", "qvalue")], row.names = FALSE)
  }
} else {
  cat("No significant enrichment\n")
}
