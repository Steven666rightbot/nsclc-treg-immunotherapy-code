# GO enrichment for B2M KO significant genes in MKI67+ Treg
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

OUTPUT_DIR <- "D:/Research/cucumber/fig5"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# B2M KO significant genes (p.adj < 0.05, excluding B2M itself)
sig_genes <- c("PTMS", "HAVCR1", "HLA-DRB5", "AC092123.1", "MATN2", 
               "SLC18B1", "RPL34-AS1", "AL603910.1", "AL592166.1", "CD74",
               "EMD", "SOWAHD", "AC012358.1", "VPS52", "AC034198.2",
               "AP001972.2", "NDUFA7", "ATG10-AS1", "AP003392.1")

# Convert to Entrez IDs
gene_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat(sprintf("Converted %d/%d genes to Entrez IDs\n", nrow(gene_ids), length(sig_genes)))
print(gene_ids)

# GO BP enrichment
go_bp <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2,
  minGSSize = 5,
  maxGSSize = 500
)

if (!is.null(go_bp) && nrow(go_bp) > 0) {
  cat(sprintf("GO BP: %d significant terms\n", nrow(go_bp)))
  
  # Save full results
  write.csv(as.data.frame(go_bp), file.path(OUTPUT_DIR, "go_bp_b2m_ko.csv"), row.names = FALSE)
  
  # Save top 20 for figure
  go_top <- head(as.data.frame(go_bp), 20)
  write.csv(go_top, file.path(OUTPUT_DIR, "go_bp_b2m_ko_top20.csv"), row.names = FALSE)
  
  # Print top terms
  print(go_top[, c("ID", "Description", "p.adjust", "Count")])
} else {
  cat("No significant GO BP terms found\n")
}

# GO MF enrichment
go_mf <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2,
  minGSSize = 5,
  maxGSSize = 500
)

if (!is.null(go_mf) && nrow(go_mf) > 0) {
  go_mf_top <- head(as.data.frame(go_mf), 20)
  write.csv(go_mf_top, file.path(OUTPUT_DIR, "go_mf_b2m_ko_top20.csv"), row.names = FALSE)
  cat(sprintf("GO MF: %d significant terms\n", nrow(go_mf)))
  print(go_mf_top[, c("ID", "Description", "p.adjust", "Count")])
} else {
  cat("No significant GO MF terms found\n")
}

# KEGG enrichment
kegg <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.1,
  minGSSize = 5,
  maxGSSize = 500
)

if (!is.null(kegg) && nrow(kegg) > 0) {
  write.csv(as.data.frame(kegg), file.path(OUTPUT_DIR, "kegg_b2m_ko.csv"), row.names = FALSE)
  cat(sprintf("KEGG: %d significant terms\n", nrow(kegg)))
  print(kegg[, c("ID", "Description", "p.adjust", "Count")])
} else {
  cat("No significant KEGG terms found\n")
}

cat("DONE!\n")
