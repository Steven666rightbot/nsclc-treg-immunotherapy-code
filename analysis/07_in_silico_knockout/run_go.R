# ============================================================================
# GO BP Enrichment — CCR8+ / MKI67+ Treg: Responder vs Non-responder DE genes
# ============================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

INPUT_DIR <- "D:/Research/tomato/go_exploration/input"
OUTPUT_DIR <- "D:/Research/tomato/go_exploration/output"

# ----------------------------------------------------------------------------
# Load gene lists
# ----------------------------------------------------------------------------
read_genes <- function(fname) {
  path <- file.path(INPUT_DIR, fname)
  scan(path, what = "character", quiet = TRUE)
}

genesets <- list(
  ccr8_up_R = read_genes("ccr8_up_R.txt"),
  ccr8_down_R = read_genes("ccr8_down_R.txt"),
  mki67_up_R = read_genes("mki67_up_R.txt"),
  mki67_down_R = read_genes("mki67_down_R.txt")
)

bg <- read_genes("background_genes.txt")

cat("Gene sets loaded:\n")
for (nm in names(genesets)) {
  cat(sprintf("  %s: %d genes\n", nm, length(genesets[[nm]])))
}
cat(sprintf("  background: %d genes\n", length(bg)))

# ----------------------------------------------------------------------------
# Run GO BP enrichment
# ----------------------------------------------------------------------------
run_enrich <- function(genes, bg_genes, label) {
  cat(sprintf("\nRunning GO BP for %s...\n", label))
  
  ego <- enrichGO(
    gene = genes,
    universe = bg_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    out_df <- as.data.frame(ego)
    cat(sprintf("  %d enriched terms\n", nrow(out_df)))
    write.csv(out_df, file.path(OUTPUT_DIR, sprintf("go_%s.csv", label)), row.names = FALSE)
    
    # Barplot top 15
    if (nrow(out_df) >= 5) {
      p <- barplot(ego, showCategory = 15, title = sprintf("GO BP: %s", label)) +
        theme(axis.text.y = element_text(size = 8))
      ggsave(file.path(OUTPUT_DIR, sprintf("go_%s_barplot.pdf", label)), p, width = 8, height = 10)
      ggsave(file.path(OUTPUT_DIR, sprintf("go_%s_barplot.png", label)), p, width = 8, height = 10, dpi = 300)
    }
    
    # Dotplot top 15
    if (nrow(out_df) >= 5) {
      p2 <- dotplot(ego, showCategory = 15, title = sprintf("GO BP: %s", label)) +
        theme(axis.text.y = element_text(size = 8))
      ggsave(file.path(OUTPUT_DIR, sprintf("go_%s_dotplot.pdf", label)), p2, width = 8, height = 10)
      ggsave(file.path(OUTPUT_DIR, sprintf("go_%s_dotplot.png", label)), p2, width = 8, height = 10, dpi = 300)
    }
    
    return(ego)
  } else {
    cat("  No significant enrichment\n")
    return(NULL)
  }
}

results <- list()
for (nm in names(genesets)) {
  results[[nm]] <- run_enrich(genesets[[nm]], bg, nm)
}

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
cat("\n", rep("=", 60), "\n", sep = "")
cat("SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")
for (nm in names(results)) {
  res <- results[[nm]]
  if (!is.null(res)) {
    df <- as.data.frame(res)
    cat(sprintf("\n%s (%d sig genes): %d GO terms\n", nm, length(genesets[[nm]]), nrow(df)))
    cat("Top 10 terms:\n")
    print(df[1:min(10, nrow(df)), c("Description", "pvalue", "qvalue", "Count")], row.names = FALSE)
  } else {
    cat(sprintf("\n%s: no enrichment\n", nm))
  }
}

cat("\nAll results saved to:", OUTPUT_DIR, "\n")
