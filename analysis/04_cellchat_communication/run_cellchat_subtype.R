# ============================================================================
# CellChat Analysis — GSE243013 Weak vs Strong with Treg subtypes
# Fine cell type annotation: Treg_FOXP3, Treg_CCR8, Treg_MKI67 + others
# ============================================================================

library(CellChat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

DATA_DIR <- "D:/Research/tomato/go_exploration/cellchat_input_subtype"
OUT_DIR <- "D:/Research/tomato/go_exploration/cellchat_output_subtype"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Step 1: Load data
# ----------------------------------------------------------------------------
cat("Loading expression matrix...\n")
X <- Matrix::readMM(file.path(DATA_DIR, "matrix.mtx"))

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
  
  # Create CellChat object — use cell_type (fine annotation)
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

# ----------------------------------------------------------------------------
# Step 4: Treg-specific analysis
# ----------------------------------------------------------------------------
cat("\n========== Treg-specific analysis ==========\n")

# Extract Treg-related interactions
extract_treg_interactions <- function(cc, group_name) {
  groups <- levels(cc@idents)
  treg_groups <- c("Treg_FOXP3", "Treg_CCR8", "Treg_MKI67")
  
  # Source = Treg
  treg_source <- c()
  for (si in 1:length(groups)) {
    for (ti in 1:length(groups)) {
      if (groups[si] %in% treg_groups) {
        net <- cc@net$prob[si, ti, ]
        if (sum(net) > 0) {
          for (lr in 1:length(net)) {
            if (net[lr] > 0) {
              lr_name <- cc@net$LRs[lr]
              pw <- cc@DB$interaction[lr_name, "pathway_name"]
              treg_source <- rbind(treg_source, data.frame(
                group = group_name,
                source = groups[si],
                target = groups[ti],
                lr_pair = lr_name,
                pathway = pw,
                prob = net[lr]
              ))
            }
          }
        }
      }
    }
  }
  return(treg_source)
}

treg_strong <- extract_treg_interactions(cc_strong, "Strong")
treg_weak <- extract_treg_interactions(cc_weak, "Weak")
treg_all <- rbind(treg_strong, treg_weak)

if (nrow(treg_all) > 0) {
  write.csv(treg_all, file.path(OUT_DIR, "treg_interactions.csv"), row.names = FALSE)
  cat(sprintf("Saved %d Treg interactions\n", nrow(treg_all)))
  
  # Summarize by pathway
  treg_pw <- treg_all %>%
    group_by(pathway, source, target, group) %>%
    summarise(total_prob = sum(prob), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = total_prob, values_fill = 0) %>%
    mutate(diff = Weak - Strong, fold = ifelse(Strong > 0, Weak/Strong, ifelse(Weak > 0, Inf, 1))) %>%
    arrange(desc(diff))
  
  write.csv(treg_pw, file.path(OUT_DIR, "treg_pathway_comparison.csv"), row.names = FALSE)
  cat("Saved Treg pathway comparison\n")
  
  print(treg_pw %>% head(20), row.names = FALSE)
}

cat("\n========== DONE ==========\n")
cat(sprintf("All results in: %s\n", OUT_DIR))
