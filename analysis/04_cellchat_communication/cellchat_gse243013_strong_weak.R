# =============================================================================
# CellChat on GSE243013 (stratified downsampled)
# Strong vs Weak, sub_cell_type grouping
# =============================================================================

library(CellChat)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(Matrix)
library(data.table)

INPUT_DIR  <- "D:/Research/tomato/data/gse243013_cellchat_input"
OUTPUT_DIR <- "D:/Research/tomato/results/cellchat_gse243013"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ======================== Read data ========================
cat("Reading Strong group...\n")
meta_strong <- fread(file.path(INPUT_DIR, "strong_metadata.csv"))
expr_strong <- readMM(file.path(INPUT_DIR, "strong_matrix.mtx"))
expr_strong <- as(expr_strong, "dgCMatrix")
genes <- fread(file.path(INPUT_DIR, "genes.tsv"), header = FALSE)$V1
barcodes_s <- fread(file.path(INPUT_DIR, "strong_barcodes.tsv"), header = FALSE)$V1
rownames(expr_strong) <- genes
colnames(expr_strong) <- barcodes_s
cat(sprintf("Strong: %d genes x %d cells\n", nrow(expr_strong), ncol(expr_strong)))

cat("\nReading Weak group...\n")
meta_weak <- fread(file.path(INPUT_DIR, "weak_metadata.csv"))
expr_weak <- readMM(file.path(INPUT_DIR, "weak_matrix.mtx"))
expr_weak <- as(expr_weak, "dgCMatrix")
barcodes_w <- fread(file.path(INPUT_DIR, "weak_barcodes.tsv"), header = FALSE)$V1
rownames(expr_weak) <- genes
colnames(expr_weak) <- barcodes_w
cat(sprintf("Weak: %d genes x %d cells\n", nrow(expr_weak), ncol(expr_weak)))

# ======================== Run CellChat for each group ========================
run_cellchat <- function(expr, meta, name) {
  cat(sprintf("\n========== Running CellChat: %s ==========\n", name))
  
  cc <- createCellChat(object = expr, meta = meta, group.by = "cell_type")
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
  
  # Save RDS
  saveRDS(cc, file.path(OUTPUT_DIR, sprintf("cellchat_%s.rds", name)))
  cat(sprintf("  Saved: cellchat_%s.rds\n", name))
  
  # Extract and save interactions CSV
  df <- subsetCommunication(cc)
  write.csv(df, file.path(OUTPUT_DIR, sprintf("interactions_%s.csv", name)), row.names = FALSE)
  cat(sprintf("  Saved: interactions_%s.csv (%d rows)\n", name, nrow(df)))
  
  return(cc)
}

cc_strong <- run_cellchat(expr_strong, meta_strong, "strong")
cc_weak <- run_cellchat(expr_weak, meta_weak, "weak")

# ======================== Compare Treg-related pathways ========================
cat("\n========== Treg-related pathway summary ==========\n")

df_s <- subsetCommunication(cc_strong)
df_w <- subsetCommunication(cc_weak)

# Incoming to Treg
s_in <- df_s[df_s$target == "Treg", ]
w_in <- df_w[df_w$target == "Treg", ]

s_path <- aggregate(prob ~ pathway_name, data = s_in, FUN = sum)
w_path <- aggregate(prob ~ pathway_name, data = w_in, FUN = sum)
merged <- merge(s_path, w_path, by = "pathway_name", all = TRUE, suffixes = c("_strong", "_weak"))
merged$ratio <- merged$prob_weak / (merged$prob_strong + 1e-15)
merged$diff <- merged$prob_weak - merged$prob_strong
merged <- merged[order(-merged$diff), ]

cat("\nPathways targeting Treg (sorted by Weak - Strong diff):\n")
print(merged, row.names = FALSE)

write.csv(merged, file.path(OUTPUT_DIR, "treg_incoming_pathways.csv"), row.names = FALSE)

# Outgoing from Treg
s_out <- df_s[df_s$source == "Treg", ]
w_out <- df_w[df_w$source == "Treg", ]

s_path <- aggregate(prob ~ pathway_name, data = s_out, FUN = sum)
w_path <- aggregate(prob ~ pathway_name, data = w_out, FUN = sum)
merged_out <- merge(s_path, w_path, by = "pathway_name", all = TRUE, suffixes = c("_strong", "_weak"))
merged_out$ratio <- merged_out$prob_weak / (merged_out$prob_strong + 1e-15)
merged_out$diff <- merged_out$prob_weak - merged_out$prob_strong
merged_out <- merged_out[order(-merged_out$diff), ]

cat("\nPathways from Treg (sorted by Weak - Strong diff):\n")
print(merged_out, row.names = FALSE)

write.csv(merged_out, file.path(OUTPUT_DIR, "treg_outgoing_pathways.csv"), row.names = FALSE)

cat(sprintf("\nAll results saved to: %s\n", OUTPUT_DIR))
