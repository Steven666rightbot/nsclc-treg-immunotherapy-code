# =============================================================================
# scTenifoldKnk Receptor KO using preprocessed 207422 Treg data
# Strategy: use user's preprocessed RDS + metadata splitting
# Expected: ~5 min per KO, 10 KOs = ~50 min total
# =============================================================================

library(scTenifoldKnk)
library(data.table)
library(ggplot2)
library(dplyr)

set.seed(42)

OUTPUT_DIR <- "D:/Research/tomato/results/scTenifoldKnk_ko"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

RECEPTORS <- c("CD44", "ITGB1", "ITGA4", "DDR1", "DDR2")
KEY_GENES <- c("FOXP3", "CTLA4", "IL2RA", "ACTA2", "MYL9", "VCL", "TLN1", "TRPV4", "MOB1B")

# ======================== 加载预处理数据 ========================

cat("Loading preprocessed Treg data...\n")
treg_file <- list.files("D:/Research/potato/data/roosted", pattern = "Treg\\.rds$", full.names = TRUE)[1]
treg_all <- readRDS(treg_file)
cat(sprintf("  Total: %d genes x %d cells\n", nrow(treg_all), ncol(treg_all)))

# Load metadata to split weak/strong
meta_weak <- fread("D:/Research/tomato/data/gse207422_cellchat_input_corrected/weak_metadata.csv")
meta_strong <- fread("D:/Research/tomato/data/gse207422_cellchat_input_corrected/strong_metadata.csv")

weak_barcodes <- meta_weak[[1]]
strong_barcodes <- meta_strong[[1]]

cat(sprintf("  Weak metadata: %d cells\n", length(weak_barcodes)))
cat(sprintf("  Strong metadata: %d cells\n", length(strong_barcodes)))

# Split
treg_weak <- treg_all[, colnames(treg_all) %in% weak_barcodes]
treg_strong <- treg_all[, colnames(treg_all) %in% strong_barcodes]

cat(sprintf("  Matched Weak: %d cells\n", ncol(treg_weak)))
cat(sprintf("  Matched Strong: %d cells\n", ncol(treg_strong)))

# ======================== 预处理函数（借鉴成功脚本）====================

prep_data <- function(mat, label) {
  cat(sprintf("\n[%s] Preprocessing...\n", label))
  
  # Gene filter: min_cells=500 (from successful script)
  min_cells <- 500
  filtered_genes <- apply(mat, 1, function(x) sum(x > 0) >= min_cells)
  mat <- mat[filtered_genes, ]
  cat(sprintf("  Retained %d genes (min_cells=%d)\n", nrow(mat), min_cells))
  
  # log1p normalization (from successful script)
  cat("  Applying log1p...\n")
  mat <- log1p(mat)
  
  cat(sprintf("  Final: %d genes x %d cells\n", nrow(mat), ncol(mat)))
  return(mat)
}

treg_weak <- prep_data(treg_weak, "Weak")
treg_strong <- prep_data(treg_strong, "Strong")

# ======================== KO 函数 ========================

run_ko <- function(count_matrix, ko_gene, group_name) {
  cat(sprintf("\n========== KO: %s (%s) ==========\n", ko_gene, group_name))
  
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
  
  cat("  Running scTenifoldKnk (nNet=5, nCells=300, K=2)...\n")
  start_time <- Sys.time()
  
  result <- tryCatch({
    scTenifoldKnk(
      countMatrix = count_matrix,
      gKO = ko_gene,
      qc_mtThreshold = 0.1,
      qc_minLSize = 1000,
      nc_lambda = 0,
      nc_nNet = 5,
      nc_nCells = 300,
      nc_nComp = 2,
      nc_scaleScores = TRUE,
      nc_symmetric = FALSE,
      nc_q = 0.9,
      td_K = 2,
      td_maxIter = 1000,
      td_maxError = 1e-05,
      td_nDecimal = 3,
      ma_nDim = 2
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
    diff_reg$Group <- group_name
    
    key_hits <- diff_reg[diff_reg$gene %in% KEY_GENES, ]
    if (nrow(key_hits) > 0) {
      cat(sprintf("  KEY HITS (%d):\n", nrow(key_hits)))
      print(key_hits[, c("gene", "Z", "p.value")])
    }
    
    out_file <- file.path(OUTPUT_DIR, sprintf("ko_%s_%s.csv", group_name, ko_gene))
    write.csv(diff_reg, out_file, row.names = FALSE)
    cat(sprintf("  Saved: %s\n", out_file))
    return(diff_reg)
  }
  
  return(NULL)
}

# ======================== 主流程 ========================

cat("\n", rep("=", 60), "\n", sep = "")
cat("Starting 10 KOs (5 receptors x 2 groups)\n")
cat("Expected total: ~50 minutes\n")
cat(rep("=", 60), "\n", sep = "")

all_results <- list()

for (group_name in c("weak", "strong")) {
  count_matrix <- if (group_name == "weak") treg_weak else treg_strong
  
  for (rec in RECEPTORS) {
    res <- run_ko(count_matrix, rec, group_name)
    if (!is.null(res)) {
      all_results[[paste(group_name, rec, sep = "_")]] <- res
    }
  }
}

# ======================== 汇总 ========================

cat("\n", rep("=", 60), "\n", sep = "")
cat("SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")

summary_df <- data.frame()
for (name in names(all_results)) {
  res <- all_results[[name]]
  key_hits <- res[res$gene %in% KEY_GENES, ]
  if (nrow(key_hits) > 0) {
    summary_df <- rbind(summary_df, key_hits)
  }
}

if (nrow(summary_df) > 0) {
  summary_df <- summary_df %>% arrange(Group, KO_Gene, p.value)
  print(summary_df[, c("Group", "KO_Gene", "gene", "Z", "p.value", "p.adj")], row.names = FALSE)
  write.csv(summary_df, file.path(OUTPUT_DIR, "key_genes_summary.csv"), row.names = FALSE)
} else {
  cat("No key genes perturbed.\n")
}

# Visualization
cat("\nGenerating plots...\n")
if (nrow(summary_df) > 0) {
  p <- ggplot(summary_df, aes(x = gene, y = Z, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
    facet_wrap(~KO_Gene, ncol = 2) +
    scale_fill_manual(values = c("weak" = "#b08ebf", "strong" = "#8bc98b")) +
    labs(title = "Receptor KO Downstream Perturbations", x = "Gene", y = "Z-score", fill = "Group") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          strip.background = element_rect(fill = "gray90", color = "black"),
          strip.text = element_text(face = "bold"))
  
  ggsave(file.path(OUTPUT_DIR, "ko_perturbation_zscores.png"), p, width = 12, height = 8, dpi = 300)
  cat("  Saved: ko_perturbation_zscores.png\n")
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("DONE! Results in:", OUTPUT_DIR, "\n")
cat(rep("=", 60), "\n", sep = "")
