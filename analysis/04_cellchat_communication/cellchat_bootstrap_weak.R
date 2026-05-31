# =============================================================================
# CellChat Bootstrap: Downsample Weak to match Strong cell count (23,097)
# Stratified by cell_type to maintain proportions
# Run 10 iterations, report median + 95% CI for Fibroblast->Treg
# =============================================================================

library(CellChat)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(Matrix)
library(data.table)

INPUT_DIR  <- "D:/Research/tomato/data/gse207422_cellchat_input_corrected"
OUTPUT_DIR <- "D:/Research/tomato/results/cellchat_ko_corrected/bootstrap"
STRONG_RDS <- "D:/Research/tomato/results/cellchat_ko_corrected/cellchat_strong.rds"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_N <- 23097  # Match Strong group
N_BOOT   <- 10
set.seed(42)

# ======================== Read full Weak data ========================
cat("Reading Weak group data...\n")
meta_weak <- fread(file.path(INPUT_DIR, "weak_metadata.csv"))
expr_weak <- readMM(file.path(INPUT_DIR, "weak_matrix.mtx"))
genes <- fread(file.path(INPUT_DIR, "weak_genes.tsv"), header = FALSE)$V1
barcodes <- fread(file.path(INPUT_DIR, "weak_barcodes.tsv"), header = FALSE)$V1
rownames(expr_weak) <- genes
colnames(expr_weak) <- barcodes
expr_weak <- as(expr_weak, "dgCMatrix")

cat(sprintf("Full Weak: %d genes x %d cells\n", nrow(expr_weak), ncol(expr_weak)))
cat(sprintf("Cell type distribution:\n"))
print(table(meta_weak$cell_type))

# ======================== Read Strong baseline ========================
cat("\nReading Strong baseline...\n")
strong_base <- readRDS(STRONG_RDS)
ft_strong <- subsetCommunication(strong_base)
ft_strong <- ft_strong[ft_strong$source == "Fibroblast" & ft_strong$target == "Treg", ]
strong_prob <- if (nrow(ft_strong) > 0) sum(ft_strong$prob) else 0
strong_n <- nrow(ft_strong)
cat(sprintf("Strong Fib->Treg: %d interactions, prob = %.6e\n", strong_n, strong_prob))

# ======================== Bootstrap loop ========================
results <- data.frame(
  iter = integer(),
  n_cells = integer(),
  n_fibroblast = integer(),
  n_treg = integer(),
  fib_treg_n = integer(),
  fib_treg_prob = numeric(),
  stringsAsFactors = FALSE
)

cat("\n", rep("=", 60), "\n", sep = "")
cat("Starting bootstrap CellChat analysis (", N_BOOT, " iterations)\n", sep = "")
cat(rep("=", 60), "\n", sep = "")

for (i in 1:N_BOOT) {
  cat(sprintf("\n========== Bootstrap %d / %d ==========\n", i, N_BOOT))
  
  # Stratified downsample: maintain cell_type proportions
  # Sample exactly TARGET_N cells proportionally
  ct_counts <- table(meta_weak$cell_type)
  ct_props <- ct_counts / sum(ct_counts)
  n_per_ct <- round(ct_props * TARGET_N)
  
  # Ensure sum equals TARGET_N (adjust largest group)
  diff <- TARGET_N - sum(n_per_ct)
  if (diff != 0) {
    largest_ct <- names(ct_counts)[which.max(ct_counts)]
    n_per_ct[largest_ct] <- n_per_ct[largest_ct] + diff
  }
  
  sampled_cells <- c()
  for (ct in names(n_per_ct)) {
    ct_cells <- meta_weak$cell[meta_weak$cell_type == ct]
    n_samp <- min(n_per_ct[ct], length(ct_cells))
    if (n_samp > 0) {
      sampled_cells <- c(sampled_cells, sample(ct_cells, n_samp))
    }
  }
  
  # Extract subset
  col_idx <- match(sampled_cells, barcodes)
  expr_sub <- expr_weak[, col_idx]
  meta_sub <- meta_weak[meta_weak$cell %in% sampled_cells, ]
  colnames(expr_sub) <- sampled_cells
  
  n_fib <- sum(meta_sub$cell_type == "Fibroblast")
  n_treg <- sum(meta_sub$cell_type == "Treg")
  cat(sprintf("Sampled: %d cells (Fib=%d, Treg=%d)\n", length(sampled_cells), n_fib, n_treg))
  
  # Run CellChat pipeline
  cat("  -> createCellChat...\n")
  cc <- createCellChat(object = expr_sub, meta = meta_sub, group.by = "cell_type")
  cc@DB <- CellChatDB.human
  
  cat("  -> subsetData...\n")
  cc <- subsetData(cc)
  cat("  -> identifyOverExpressedGenes...\n")
  cc <- identifyOverExpressedGenes(cc)
  cat("  -> identifyOverExpressedInteractions...\n")
  cc <- identifyOverExpressedInteractions(cc)
  cat("  -> computeCommunProb (~10-15 min)...\n")
  cc <- computeCommunProb(cc, raw.use = TRUE)
  cat("  -> filterCommunication...\n")
  cc <- filterCommunication(cc, min.cells = 10)
  cat("  -> computeCommunProbPathway...\n")
  cc <- computeCommunProbPathway(cc)
  cat("  -> aggregateNet...\n")
  cc <- aggregateNet(cc)
  
  # Extract Fib->Treg
  df <- subsetCommunication(cc)
  ft <- df[df$source == "Fibroblast" & df$target == "Treg", ]
  prob <- if (nrow(ft) > 0) sum(ft$prob) else 0
  
  results <- rbind(results, data.frame(
    iter = i,
    n_cells = length(sampled_cells),
    n_fibroblast = n_fib,
    n_treg = n_treg,
    fib_treg_n = nrow(ft),
    fib_treg_prob = prob,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("  Fib->Treg: %d interactions, prob = %.6e\n", nrow(ft), prob))
  
  # Save RDS
  saveRDS(cc, file.path(OUTPUT_DIR, sprintf("cellchat_weak_bootstrap_%02d.rds", i)))
  
  # Clean up to free memory
  rm(cc, expr_sub, meta_sub, df, ft)
  gc()
}

# ======================== Summary ========================
cat("\n", rep("=", 60), "\n", sep = "")
cat("BOOTSTRAP SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")

print(results, row.names = FALSE)

med_prob <- median(results$fib_treg_prob)
q025 <- quantile(results$fib_treg_prob, 0.025)
q975 <- quantile(results$fib_treg_prob, 0.975)

ratio <- med_prob / strong_prob
ratio_lower <- q025 / strong_prob
ratio_upper <- q975 / strong_prob

cat("\nFibroblast -> Treg Communication Probability:\n")
cat(sprintf("  Strong (baseline):     %.6e  (%d interactions)\n", strong_prob, strong_n))
cat(sprintf("  Weak (bootstrap med):  %.6e  (%d interactions)\n", med_prob, round(median(results$fib_treg_n))))
cat(sprintf("  95%% CI:                [%.6e, %.6e]\n", q025, q975))
cat(sprintf("\nWeak / Strong ratio: %.2f×  [95%% CI: %.2f×, %.2f×]\n", 
            ratio, ratio_lower, ratio_upper))

write.csv(results, file.path(OUTPUT_DIR, "bootstrap_results.csv"), row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n", OUTPUT_DIR))
