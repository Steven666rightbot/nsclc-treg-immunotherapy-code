# =============================================================================
# scTenifoldKnk Virtual Knockout — Merged Treg (Weak + Strong)
# Strategy: use ALL Treg cells from GSE207422 without response stratification
# Targets: receptors (CD44/ITGB1/DDR1/DDR2/ITGA4) + TFs (MEF2C/TEAD1/YAP1/WWTR1)
# Post-processing: GO/KEGG enrichment + volcano plot + enrichment dotplot
# =============================================================================

library(scTenifoldKnk)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

set.seed(42)

OUTPUT_DIR <- "D:/Research/tomato/results/scTenifoldKnk_ko_merged"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Key gene sets
CONTRACTILE_GENES <- c("ACTA2", "MYL9", "VCL", "TLN1", "TRPV4", "MOB1B")
MECHANOSENSING_GENES <- c("DDR1", "DDR2", "ITGA1", "SRC", "RHOA", "ROCK1", "ROCK2")
TREG_MARKERS <- c("FOXP3", "IL2RA", "CTLA4")
YAP_TAZ_GENES <- c("YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4")
KEY_GENES <- unique(c(CONTRACTILE_GENES, MECHANOSENSING_GENES, TREG_MARKERS, YAP_TAZ_GENES))

TARGET_RECEPTORS <- c("CD44", "ITGB1", "DDR1", "DDR2", "ITGA4")
TARGET_TFS <- c("MEF2C", "TEAD1", "YAP1", "WWTR1")

# ======================== Load data ========================

cat("Loading preprocessed Treg data...\n")
treg_file <- "D:/Research/potato/data/roosted/treg_207422.rds"
treg_all <- readRDS(treg_file)
cat(sprintf("  Total: %d genes x %d cells\n", nrow(treg_all), ncol(treg_all)))

# ======================== Preprocessing ========================

prep_data <- function(mat, label) {
  cat(sprintf("\n[%s] Preprocessing...\n", label))
  
  # Lower min_cells to retain more genes for denser network
  min_cells <- 200
  filtered_genes <- apply(mat, 1, function(x) sum(x > 0) >= min_cells)
  mat <- mat[filtered_genes, ]
  cat(sprintf("  Retained %d genes (min_cells=%d)\n", nrow(mat), min_cells))
  
  # Force-keep key genes that might have been filtered out
  missing_key <- KEY_GENES[!(KEY_GENES %in% rownames(mat))]
  if (length(missing_key) > 0) {
    cat(sprintf("  Force-keeping %d key genes: %s\n", length(missing_key), paste(missing_key, collapse=", ")))
    mat_key <- mat_all[rownames(mat_all) %in% missing_key, colnames(mat), drop=FALSE]
    if (nrow(mat_key) > 0) {
      mat <- rbind(mat, mat_key)
    }
  }
  
  # log1p normalization
  cat("  Applying log1p...\n")
  mat <- log1p(mat)
  
  cat(sprintf("  Final: %d genes x %d cells\n", nrow(mat), ncol(mat)))
  return(mat)
}

# Keep original for force-keeping
mat_all <- treg_all
treg_merged <- prep_data(treg_all, "Merged Treg")

# ======================== KO function ========================

run_ko <- function(count_matrix, ko_gene, category) {
  cat(sprintf("\n========== KO: %s (%s) ==========\n", ko_gene, category))
  
  if (!(ko_gene %in% rownames(count_matrix))) {
    cat(sprintf("  WARNING: %s not found, skipping\n", ko_gene))
    return(NULL)
  }
  
  expr_level <- mean(count_matrix[ko_gene, ])
  cat(sprintf("  %s mean expression: %.4f\n", ko_gene, expr_level))
  
  if (expr_level < 0.01) {
    cat(sprintf("  WARNING: %s too low, skipping\n", ko_gene))
    return(NULL)
  }
  
  cat("  Running scTenifoldKnk (nNet=10, nCells=500, K=3)...\n")
  start_time <- Sys.time()
  
  result <- tryCatch({
    scTenifoldKnk(
      countMatrix = count_matrix,
      gKO = ko_gene,
      qc_mtThreshold = 0.1,
      qc_minLSize = 1000,
      nc_lambda = 0,
      nc_nNet = 10,        # increased for stability
      nc_nCells = 500,     # increased for stability
      nc_nComp = 3,
      nc_scaleScores = TRUE,
      nc_symmetric = FALSE,
      nc_q = 0.9,
      td_K = 3,
      td_maxIter = 1000,
      td_maxError = 1e-05,
      td_nDecimal = 3,
      ma_nDim = 3
    )
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("  Finished in %.1f min\n", as.numeric(elapsed)))
  
  if (!is.null(result)) {
    diff_reg <- result$diffRegulation
    cat(sprintf("  %d perturbed genes\n", nrow(diff_reg)))
    
    diff_reg$KO_Gene <- ko_gene
    diff_reg$Category <- category
    
    n_sig <- sum(diff_reg$p.adj < 0.05, na.rm = TRUE)
    cat(sprintf("  Significant (p.adj<0.05): %d\n", n_sig))
    
    key_hits <- diff_reg[diff_reg$gene %in% KEY_GENES, ]
    if (nrow(key_hits) > 0) {
      cat(sprintf("  KEY HITS (%d):\n", nrow(key_hits)))
      print(key_hits[, c("gene", "Z", "p.value", "p.adj")])
    }
    
    out_file <- file.path(OUTPUT_DIR, sprintf("ko_merged_%s_%s.csv", category, ko_gene))
    write.csv(diff_reg, out_file, row.names = FALSE)
    cat(sprintf("  Saved: %s\n", out_file))
    return(diff_reg)
  }
  
  return(NULL)
}

# ======================== Run all KOs ========================

cat("\n", rep("=", 60), "\n", sep = "")
cat("Starting Receptor KOs (n=5)\n")
cat(rep("=", 60), "\n", sep = "")

all_results <- list()

for (rec in TARGET_RECEPTORS) {
  res <- run_ko(treg_merged, rec, "receptor")
  if (!is.null(res)) {
    all_results[[paste("receptor", rec, sep = "_")]] <- res
  }
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("Starting TF KOs (n=4)\n")
cat(rep("=", 60), "\n", sep = "")

for (tf in TARGET_TFS) {
  res <- run_ko(treg_merged, tf, "TF")
  if (!is.null(res)) {
    all_results[[paste("TF", tf, sep = "_")]] <- res
  }
}

# ======================== Summarize key hits ========================

cat("\n", rep("=", 60), "\n", sep = "")
cat("SUMMARY: Key gene perturbations\n")
cat(rep("=", 60), "\n", sep = "")

key_summary <- data.frame()
for (name in names(all_results)) {
  res <- all_results[[name]]
  hits <- res[res$gene %in% KEY_GENES & res$p.adj < 0.05, ]
  if (nrow(hits) > 0) {
    key_summary <- rbind(key_summary, hits)
  }
}

if (nrow(key_summary) > 0) {
  key_summary <- key_summary %>% arrange(Category, KO_Gene, p.adj)
  print(key_summary[, c("Category", "KO_Gene", "gene", "Z", "p.value", "p.adj")], row.names = FALSE)
  write.csv(key_summary, file.path(OUTPUT_DIR, "key_genes_summary_merged.csv"), row.names = FALSE)
} else {
  cat("No key genes significantly perturbed.\n")
}

# ======================== GO/KEGG enrichment ========================

cat("\n", rep("=", 60), "\n", sep = "")
cat("GO/KEGG Enrichment Analysis\n")
cat(rep("=", 60), "\n", sep = "")

run_enrichment <- function(diff_reg, ko_name, category) {
  sig_genes <- diff_reg$gene[diff_reg$p.adj < 0.05]
  cat(sprintf("\n%s %s: %d significant genes for enrichment\n", category, ko_name, length(sig_genes)))
  
  if (length(sig_genes) < 10) {
    cat("  Too few genes, skipping enrichment.\n")
    return(NULL)
  }
  
  # GO BP
  go_bp <- tryCatch({
    enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  }, error = function(e) NULL)
  
  # KEGG
  kegg <- tryCatch({
    # Convert SYMBOL to ENTREZID
    entrez <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
    enrichKEGG(gene = entrez, organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  }, error = function(e) NULL)
  
  if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
    out <- file.path(OUTPUT_DIR, sprintf("go_bp_%s_%s.csv", category, ko_name))
    write.csv(as.data.frame(go_bp), out, row.names = FALSE)
    cat(sprintf("  GO BP: %d terms -> %s\n", nrow(as.data.frame(go_bp)), out))
  }
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    out <- file.path(OUTPUT_DIR, sprintf("kegg_%s_%s.csv", category, ko_name))
    write.csv(as.data.frame(kegg), out, row.names = FALSE)
    cat(sprintf("  KEGG: %d terms -> %s\n", nrow(as.data.frame(kegg)), out))
  }
  
  return(list(go_bp = go_bp, kegg = kegg))
}

enrich_list <- list()
for (name in names(all_results)) {
  parts <- strsplit(name, "_")[[1]]
  category <- parts[1]
  ko_name <- paste(parts[-1], collapse = "_")
  res <- all_results[[name]]
  enrich_list[[name]] <- run_enrichment(res, ko_name, category)
}

# ======================== Visualization ========================

cat("\n", rep("=", 60), "\n", sep = "")
cat("Generating plots...\n")
cat(rep("=", 60), "\n", sep = "")

# --- Volcano plot: all KOs combined ---
all_df <- do.call(rbind, all_results)
all_df$significant <- ifelse(all_df$p.adj < 0.05, "Sig", "NS")
all_df$highlight <- ifelse(all_df$gene %in% KEY_GENES, "Key", "Other")
all_df$label <- ifelse(all_df$highlight == "Key" & all_df$significant == "Sig", all_df$gene, "")

# Focus on key genes for cleaner plot
key_df <- all_df[all_df$gene %in% KEY_GENES, ]

if (nrow(key_df) > 0) {
  p_volcano <- ggplot(key_df, aes(x = Z, y = -log10(p.adj), color = KO_Gene, shape = Category)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_text(aes(label = gene), vjust = -1, size = 3, show.legend = FALSE) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Virtual KO Perturbation of Key Genes (Merged Treg, n=1492)",
         x = "Z-score", y = "-log10(adjusted p-value)", color = "KO Target", shape = "Category") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  
  ggsave(file.path(OUTPUT_DIR, "volcano_key_genes_merged.png"), p_volcano, width = 10, height = 7, dpi = 300)
  cat("  Saved: volcano_key_genes_merged.png\n")
}

# --- Bar plot: significant perturbed gene counts per KO ---
sig_counts <- data.frame()
for (name in names(all_results)) {
  res <- all_results[[name]]
  n_sig <- sum(res$p.adj < 0.05, na.rm = TRUE)
  sig_counts <- rbind(sig_counts, data.frame(
    KO = res$KO_Gene[1],
    Category = res$Category[1],
    n_sig = n_sig,
    stringsAsFactors = FALSE
  ))
}

if (nrow(sig_counts) > 0) {
  p_counts <- ggplot(sig_counts, aes(x = reorder(KO, n_sig), y = n_sig, fill = Category)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c("receptor" = "#8bc98b", "TF" = "#b08ebf")) +
    labs(title = "Significant Perturbations per Virtual KO", x = "KO Target", y = "Count (p.adj < 0.05)") +
    theme_minimal(base_size = 12) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  
  ggsave(file.path(OUTPUT_DIR, "sig_count_per_ko_merged.png"), p_counts, width = 8, height = 5, dpi = 300)
  cat("  Saved: sig_count_per_ko_merged.png\n")
}

# --- Top GO/KEGG dotplot for ITGB1 KO (most relevant to our story) ---
itgb1_enrich <- enrich_list[["receptor_ITGB1"]]
if (!is.null(itgb1_enrich) && !is.null(itgb1_enrich$go_bp) && nrow(as.data.frame(itgb1_enrich$go_bp)) > 0) {
  p_go <- dotplot(itgb1_enrich$go_bp, showCategory = 15) +
    ggtitle("GO BP: ITGB1 KO (Merged Treg)") +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(OUTPUT_DIR, "go_bp_ITGB1_merged.png"), p_go, width = 9, height = 7, dpi = 300)
  cat("  Saved: go_bp_ITGB1_merged.png\n")
}

# --- Top GO/KEGG dotplot for YAP1 KO ---
yap1_enrich <- enrich_list[["TF_YAP1"]]
if (!is.null(yap1_enrich) && !is.null(yap1_enrich$go_bp) && nrow(as.data.frame(yap1_enrich$go_bp)) > 0) {
  p_go <- dotplot(yap1_enrich$go_bp, showCategory = 15) +
    ggtitle("GO BP: YAP1 KO (Merged Treg)") +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(OUTPUT_DIR, "go_bp_YAP1_merged.png"), p_go, width = 9, height = 7, dpi = 300)
  cat("  Saved: go_bp_YAP1_merged.png\n")
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("DONE! All results in:", OUTPUT_DIR, "\n")
cat(rep("=", 60), "\n", sep = "")
