# =============================================================================
# scTenifoldKnk Virtual Knockout of ECM Receptors in Treg Cells
# Dataset: GSE207422 (corrected annotation)
# Aggressively optimized based on successful local script
# =============================================================================

library(scTenifoldKnk)
library(Matrix)
library(data.table)
library(ggplot2)
library(dplyr)

set.seed(42)

# ======================== 参数配置 ========================
INPUT_DIR  <- "D:/Research/tomato/data/gse207422_cellchat_input_corrected"
OUTPUT_DIR <- "D:/Research/tomato/results/scTenifoldKnk_ko"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

RECEPTORS <- c("CD44", "ITGB1", "ITGA4", "DDR1", "DDR2")
KEY_GENES <- c("FOXP3", "CTLA4", "IL2RA", "ACTA2", "MYL9", "VCL", "TLN1", "TRPV4", "MOB1B")

# scTenifoldKnk parameters (from successful script)
N_NET    <- 5
N_CELLS  <- 300
N_COMP   <- 2
QC_MT    <- 0.1
QC_MINSIZE <- 1000

# Data reduction (from successful script + further optimization)
MIN_CELLS_EXPRESSED <- 500      # Same as successful script
MAX_GENES           <- 2000     # Cap genes to 2000 max
MAX_CELLS           <- 1000     # Subsample cells if too many

# ======================== 辅助函数 ========================

read_and_extract_treg <- function(prefix) {
  cat(sprintf("\n[%s] Reading data...\n", prefix))
  
  mtx_file <- file.path(INPUT_DIR, sprintf("%s_matrix.mtx", prefix))
  expr_matrix <- readMM(mtx_file)
  
  genes_file <- file.path(INPUT_DIR, sprintf("%s_genes.tsv", prefix))
  genes <- fread(genes_file, header = FALSE)$V1
  
  barcodes_file <- file.path(INPUT_DIR, sprintf("%s_barcodes.tsv", prefix))
  barcodes <- fread(barcodes_file, header = FALSE)$V1
  
  meta_file <- file.path(INPUT_DIR, sprintf("%s_metadata.csv", prefix))
  meta <- fread(meta_file)
  
  rownames(expr_matrix) <- genes
  colnames(expr_matrix) <- barcodes
  meta <- meta[match(barcodes, meta[[1]]), ]
  
  # Extract Treg cells
  treg_idx <- which(meta$cell_type == "Treg")
  cat(sprintf("  Total cells: %d, Treg cells: %d\n", ncol(expr_matrix), length(treg_idx)))
  
  if (length(treg_idx) < 100) {
    stop(sprintf("Too few Treg cells in %s: %d", prefix, length(treg_idx)))
  }
  
  treg_matrix <- expr_matrix[, treg_idx]
  
  # Convert to dense matrix
  cat("  Converting to dense matrix...\n")
  treg_dense <- as.matrix(treg_matrix)
  
  # KEY: Cell subsampling if too many (further optimization)
  if (ncol(treg_dense) > MAX_CELLS) {
    cat(sprintf("  Subsampling cells: %d -> %d\n", ncol(treg_dense), MAX_CELLS))
    set.seed(42)
    keep_cells <- sample(ncol(treg_dense), MAX_CELLS)
    treg_dense <- treg_dense[, keep_cells]
  }
  
  # KEY: Gene filter from successful script (min_cells_expressed = 500)
  cat(sprintf("  Filtering genes (expressed in >= %d cells)...\n", MIN_CELLS_EXPRESSED))
  filtered_genes <- apply(treg_dense, 1, function(x) sum(x > 0) >= MIN_CELLS_EXPRESSED)
  treg_dense <- treg_dense[filtered_genes, ]
  cat(sprintf("  Retained %d genes (min_cells=%d)\n", nrow(treg_dense), MIN_CELLS_EXPRESSED))
  
  # KEY: Cap to MAX_GENES if still too many
  if (nrow(treg_dense) > MAX_GENES) {
    cat(sprintf("  Capping to top %d genes by total expression...\n", MAX_GENES))
    gene_totals <- rowSums(treg_dense)
    top_idx <- order(gene_totals, decreasing = TRUE)[1:MAX_GENES]
    treg_dense <- treg_dense[top_idx, ]
  }
  
  # KEY: log1p normalization (from successful script)
  cat("  Applying log1p normalization...\n")
  treg_dense <- log1p(treg_dense)
  
  cat(sprintf("  Final matrix: %d genes x %d cells\n", nrow(treg_dense), ncol(treg_dense)))
  return(treg_dense)
}

run_ko <- function(count_matrix, ko_gene, group_name) {
  cat(sprintf("\n========== KO: %s (%s) ==========\n", ko_gene, group_name))
  
  if (!(ko_gene %in% rownames(count_matrix))) {
    cat(sprintf("  WARNING: %s not found in matrix, skipping\n", ko_gene))
    return(NULL)
  }
  
  expr_level <- mean(count_matrix[ko_gene, ])
  cat(sprintf("  %s mean expression: %.4f\n", ko_gene, expr_level))
  
  if (expr_level < 0.01) {
    cat(sprintf("  WARNING: %s expression too low, skipping\n", ko_gene))
    return(NULL)
  }
  
  cat(sprintf("  Running scTenifoldKnk (nNet=%d, nCells=%d, K=%d)...\n", N_NET, N_CELLS, N_COMP))
  
  result <- tryCatch({
    scTenifoldKnk(
      countMatrix = count_matrix,
      gKO = ko_gene,
      qc_mtThreshold = QC_MT,
      qc_minLSize = QC_MINSIZE,
      nc_lambda = 0,
      nc_nNet = N_NET,
      nc_nCells = N_CELLS,
      nc_nComp = N_COMP,
      nc_scaleScores = TRUE,
      nc_symmetric = FALSE,
      nc_q = 0.9,
      td_K = N_COMP,
      td_maxIter = 1000,
      td_maxError = 1e-05,
      td_nDecimal = 3,
      ma_nDim = 2
    )
  }, error = function(e) {
    cat(sprintf("  ERROR in scTenifoldKnk: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  if (!is.null(result)) {
    diff_reg <- result$diffRegulation
    cat(sprintf("  Completed. %d perturbed genes identified.\n", nrow(diff_reg)))
    
    diff_reg$KO_Gene <- ko_gene
    diff_reg$Group <- group_name
    
    key_hits <- diff_reg[diff_reg$gene %in% KEY_GENES, ]
    if (nrow(key_hits) > 0) {
      cat(sprintf("  KEY GENE HITS (%d):\n", nrow(key_hits)))
      print(key_hits[, c("gene", "Z", "p.value")])
    } else {
      cat("  No key genes in top perturbed genes.\n")
    }
    
    return(diff_reg)
  }
  
  return(NULL)
}

# ======================== 主流程 ========================

cat("=" , rep("=", 60), "\n", sep = "")
cat("scTenifoldKnk Receptor Virtual Knockout (AGGRESSIVE)\n")
cat("GSE207422 Treg Cells\n")
cat("=" , rep("=", 60), "\n", sep = "")

treg_weak <- read_and_extract_treg("weak")
treg_strong <- read_and_extract_treg("strong")

all_results <- list()

for (group_name in c("weak", "strong")) {
  count_matrix <- if (group_name == "weak") treg_weak else treg_strong
  
  for (rec in RECEPTORS) {
    res <- run_ko(count_matrix, rec, group_name)
    if (!is.null(res)) {
      all_results[[paste(group_name, rec, sep = "_")]] <- res
      
      out_file <- file.path(OUTPUT_DIR, sprintf("ko_%s_%s.csv", group_name, rec))
      write.csv(res, out_file, row.names = FALSE)
      cat(sprintf("  Saved: %s\n", out_file))
    }
  }
}

# ======================== 汇总分析 ========================

cat("\n", rep("=", 60), "\n", sep = "")
cat("SUMMARY: Key Gene Perturbations\n")
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
  summary_df <- summary_df %>% 
    arrange(Group, KO_Gene, p.value) %>%
    select(Group, KO_Gene, gene, Z, p.value, p.adj)
  
  print(summary_df, row.names = FALSE)
  write.csv(summary_df, file.path(OUTPUT_DIR, "key_genes_summary.csv"), row.names = FALSE)
} else {
  cat("No key genes were significantly perturbed in any KO condition.\n")
}

# ======================== 可视化 ========================

cat("\n", rep("=", 60), "\n", sep = "")
cat("Generating visualizations...\n")
cat(rep("=", 60), "\n", sep = "")

if (nrow(summary_df) > 0) {
  p <- ggplot(summary_df, aes(x = gene, y = Z, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
    facet_wrap(~KO_Gene, ncol = 2) +
    scale_fill_manual(values = c("weak" = "#b08ebf", "strong" = "#8bc98b")) +
    labs(
      title = "scTenifoldKnk: Receptor KO Downstream Perturbations",
      subtitle = "Z-score of key Treg genes after virtual receptor knockout",
      x = "Downstream Gene", y = "Perturbation Z-score",
      fill = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = "gray90", color = "black"),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(file.path(OUTPUT_DIR, "ko_perturbation_zscores.png"), p, width = 12, height = 8, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "ko_perturbation_zscores.pdf"), p, width = 12, height = 8)
  cat("  Saved: ko_perturbation_zscores.png\n")
}

cat("\nGenerating full perturbation heatmaps...\n")
for (name in names(all_results)) {
  res <- all_results[[name]]
  if (nrow(res) >= 10) {
    top_genes <- head(res[order(res$Z, decreasing = TRUE), ], 50)
    
    p <- ggplot(top_genes, aes(x = reorder(gene, Z), y = Z, fill = Z)) +
      geom_bar(stat = "identity", color = "black", alpha = 0.8) +
      coord_flip() +
      scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C", midpoint = 0) +
      labs(
        title = sprintf("%s: Top Perturbed Genes", name),
        x = "Gene", y = "Perturbation Z-score"
      ) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none")
    
    ggsave(file.path(OUTPUT_DIR, sprintf("top_perturbed_%s.png", name)), p, width = 8, height = 10, dpi = 300)
  }
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("scTenifoldKnk analysis complete!\n")
cat(sprintf("Results saved to: %s\n", OUTPUT_DIR))
cat(rep("=", 60), "\n", sep = "")
