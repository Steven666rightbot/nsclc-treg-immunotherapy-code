# ============================================================================
# GSEA — CCR8+ / MKI67+ Treg: Responder vs Non-responder
# No hard threshold, use full ranked list
# ============================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

INPUT_DIR <- "D:/Research/tomato/vB_revised_integrin_fyn_axis"
OUTPUT_DIR <- "D:/Research/tomato/go_exploration/output"

# ----------------------------------------------------------------------------
# Load full DE results (ranked by sign(logFC) * -log10(p))
# ----------------------------------------------------------------------------
run_gsea <- function(csv_file, label) {
  cat(sprintf("\nRunning GSEA for %s...\n", label))
  
  df <- read.csv(file.path(INPUT_DIR, csv_file))
  
  # Rank genes: sign(logFC) * -log10(pval)
  df$rank_score <- sign(df$logfoldchanges) * (-log10(df$pvals + 1e-300))
  df <- df[order(df$rank_score, decreasing = TRUE), ]
  
  geneList <- df$rank_score
  names(geneList) <- df$names
  
  # GSEA GO BP
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
    cat(sprintf("  %d significant gene sets\n", nrow(out_df)))
    write.csv(out_df, file.path(OUTPUT_DIR, sprintf("gsea_%s.csv", label)), row.names = FALSE)
    
    # Top positive (R-enriched) and negative (NR-enriched) NES
    cat("\n  Top UP in R (positive NES):\n")
    up <- out_df[out_df$NES > 0, ]
    if (nrow(up) > 0) {
      print(up[1:min(5, nrow(up)), c("Description", "NES", "pvalue", "qvalues")], row.names = FALSE)
    }
    
    cat("\n  Top UP in NR (negative NES):\n")
    down <- out_df[out_df$NES < 0, ]
    if (nrow(down) > 0) {
      print(down[1:min(5, nrow(down)), c("Description", "NES", "pvalue", "qvalues")], row.names = FALSE)
    }
    
    return(gsea)
  } else {
    cat("  No significant enrichment\n")
    return(NULL)
  }
}

results <- list()
results[["ccr8"]] <- run_gsea("de_ccr8_responder_vs_nonresponder.csv", "ccr8")
results[["mki67"]] <- run_gsea("de_mki67_responder_vs_nonresponder.csv", "mki67")

cat("\nDone. Results saved to:", OUTPUT_DIR, "\n")
