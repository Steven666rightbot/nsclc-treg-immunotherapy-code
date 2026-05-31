# ============================================================================
# CellChat Analysis — GSE243013 Weak (non-MPR) vs Strong (pCR+MPR)
# Compare cell-cell communication between responders and non-responders
# ============================================================================

library(CellChat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

DATA_DIR <- "D:/Research/tomato/go_exploration/cellchat_input"
OUT_DIR <- "D:/Research/tomato/go_exploration/cellchat_output"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Step 1: Load data
# ----------------------------------------------------------------------------
cat("Loading expression matrix...\n")
# Read MTX (genes x cells)
X <- Matrix::readMM(file.path(DATA_DIR, "matrix.mtx"))

# Read barcodes and genes
barcodes <- readLines(file.path(DATA_DIR, "barcodes.tsv"))
genes <- readLines(file.path(DATA_DIR, "genes.tsv"))

dim(X)
cat(sprintf("Matrix: %d genes x %d cells\n", nrow(X), ncol(X)))

# Create Seurat object
meta <- read.csv(file.path(DATA_DIR, "cell_annotation.csv"))
rownames(meta) <- meta$cellID

# Ensure barcodes match
common <- intersect(barcodes, meta$cellID)
cat(sprintf("Common barcodes: %d\n", length(common)))

# Subset to common
idx <- match(common, barcodes)
X <- X[, idx]
colnames(X) <- common
rownames(X) <- genes

# Subset metadata
meta <- meta[common, ]

# Create Seurat
seu <- CreateSeuratObject(counts = X, meta.data = meta)
cat(sprintf("Seurat: %d cells\n", ncol(seu)))

# Normalize and identify variable features (required by CellChat)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

# ----------------------------------------------------------------------------
# Step 2: Run CellChat for each group
# ----------------------------------------------------------------------------
run_cellchat <- function(seu_obj, group_name, group_cells) {
  cat(sprintf("\n========== Running CellChat for %s (%d cells) ==========\n", 
              group_name, length(group_cells)))
  
  sub <- subset(seu_obj, cells = group_cells)
  
  # Prepare CellChat input
  data.input <- GetAssayData(sub, layer = "data")
  meta.sub <- sub@meta.data
  
  # Create CellChat object
  cellchat <- createCellChat(object = data.input, meta = meta.sub, group.by = "cell_type")
  
  # Set database
  cellchatDB <- CellChatDB.human
  cellchat@DB <- cellchatDB
  
  # Preprocessing
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Compute communication probability
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Aggregate at pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # Save
  saveRDS(cellchat, file.path(OUT_DIR, sprintf("cellchat_%s.rds", group_name)))
  cat(sprintf("Saved: cellchat_%s.rds\n", group_name))
  
  return(cellchat)
}

# Split by response group
strong_cells <- colnames(seu)[seu$response_group == "Strong"]
weak_cells <- colnames(seu)[seu$response_group == "Weak"]

cat(sprintf("Strong cells: %d\n", length(strong_cells)))
cat(sprintf("Weak cells: %d\n", length(weak_cells)))

# Run CellChat for each group
cc_strong <- run_cellchat(seu, "Strong", strong_cells)
cc_weak <- run_cellchat(seu, "Weak", weak_cells)

# ----------------------------------------------------------------------------
# Step 3: Merge and compare
# ----------------------------------------------------------------------------
cat("\n========== Merging and comparing ==========\n")

object.list <- list(Strong = cc_strong, Weak = cc_weak)
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))

# Save merged
saveRDS(cellchat.merged, file.path(OUT_DIR, "cellchat_merged.rds"))

# Compare number of interactions
cat("\nNumber of interactions:\n")
for (name in names(object.list)) {
  net <- object.list[[name]]@net$count
  cat(sprintf("  %s: %d total interactions\n", name, sum(net)))
}

# Compare pathway-level signaling
pathway_comparison <- rankNet(cellchat.merged, mode = "comparison", 
                               stacked = TRUE, do.stat = TRUE)
ggsave(file.path(OUT_DIR, "pathway_comparison.pdf"), pathway_comparison, 
       width = 10, height = 12)
ggsave(file.path(OUT_DIR, "pathway_comparison.png"), pathway_comparison, 
       width = 10, height = 12, dpi = 300)
cat("Saved pathway comparison plots\n")

# Identify signaling changes
cat("\nSignaling changes (lifted in Weak vs Strong):\n")
# Get all pathways
pathways <- unique(c(cc_strong@netP$pathways, cc_weak@netP$pathways))
pathway_df <- data.frame()

for (pw in pathways) {
  prob_strong <- 0
  prob_weak <- 0
  
  if (pw %in% cc_strong@netP$pathways) {
    idx <- which(cc_strong@netP$pathways == pw)
    prob_strong <- sum(cc_strong@netP$prob[, , idx])
  }
  if (pw %in% cc_weak@netP$pathways) {
    idx <- which(cc_weak@netP$pathways == pw)
    prob_weak <- sum(cc_weak@netP$prob[, , idx])
  }
  
  fold <- ifelse(prob_strong > 0, prob_weak / prob_strong, ifelse(prob_weak > 0, Inf, 1))
  
  pathway_df <- rbind(pathway_df, data.frame(
    pathway = pw,
    prob_strong = prob_strong,
    prob_weak = prob_weak,
    fold_weak_vs_strong = fold,
    present_strong = pw %in% cc_strong@netP$pathways,
    present_weak = pw %in% cc_weak@netP$pathways
  ))
}

pathway_df <- pathway_df %>% arrange(desc(fold_weak_vs_strong))
write.csv(pathway_df, file.path(OUT_DIR, "pathway_comparison.csv"), row.names = FALSE)
cat("Saved pathway comparison table\n")

# Print top changes
cat("\nTop 10 pathways UP in Weak (non-responder):\n")
up_weak <- pathway_df[pathway_df$present_weak, ] %>% head(10)
print(up_weak[, c("pathway", "prob_strong", "prob_weak", "fold_weak_vs_strong")], row.names = FALSE)

cat("\nTop 10 pathways UP in Strong (responder):\n")
up_strong <- pathway_df[pathway_df$present_strong, ] %>% arrange(fold_weak_vs_strong) %>% head(10)
print(up_strong[, c("pathway", "prob_strong", "prob_weak", "fold_weak_vs_strong")], row.names = FALSE)

cat("\n========== DONE ==========\n")
cat(sprintf("All results in: %s\n", OUT_DIR))
