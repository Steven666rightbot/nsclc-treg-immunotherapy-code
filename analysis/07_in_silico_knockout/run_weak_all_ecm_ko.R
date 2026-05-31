library(CellChat)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(Matrix)
library(data.table)

INPUT_DIR  <- "D:/Research/tomato/data/gse207422_cellchat_input_corrected"
OUTPUT_DIR <- "D:/Research/tomato/results/cellchat_ko_corrected"

cat("[weak] Reading input data...\n")
expr_matrix <- readMM(file.path(INPUT_DIR, "weak_matrix.mtx"))
genes <- fread(file.path(INPUT_DIR, "weak_genes.tsv"), header=FALSE)$V1
barcodes <- fread(file.path(INPUT_DIR, "weak_barcodes.tsv"), header=FALSE)$V1
meta <- fread(file.path(INPUT_DIR, "weak_metadata.csv"))

rownames(expr_matrix) <- genes
colnames(expr_matrix) <- barcodes
expr_matrix <- as(expr_matrix, "dgCMatrix")

cat(sprintf("[weak] Matrix: %d genes x %d cells\n", nrow(expr_matrix), ncol(expr_matrix)))

# KO ALL_ECM
ko_genes <- c("COL1A1", "COL1A2", "FN1")
cat("\n========== KO_ALL_ECM: Weak ==========\n")
expr_ko <- expr_matrix
for (g in ko_genes) {
  if (g %in% rownames(expr_ko)) {
    expr_ko[g, ] <- 0
    cat(sprintf("  KO: %s (set to 0)\n", g))
  }
}

cat("  -> createCellChat...\n")
cellchat <- createCellChat(object = expr_ko, meta = meta, group.by = "cell_type")
cellchat@DB <- CellChatDB.human

cat("  -> subsetData...\n")
cellchat <- subsetData(cellchat)
cat("  -> identifyOverExpressedGenes...\n")
cellchat <- identifyOverExpressedGenes(cellchat)
cat("  -> identifyOverExpressedInteractions...\n")
cellchat <- identifyOverExpressedInteractions(cellchat)
cat("  -> computeCommunProb (this takes ~15 min)...\n")
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cat("  -> filterCommunication...\n")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cat("  -> computeCommunProbPathway...\n")
cellchat <- computeCommunProbPathway(cellchat)
cat("  -> aggregateNet...\n")
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file.path(OUTPUT_DIR, "cellchat_weak_KO_ALL_ECM.rds"))
write.csv(subsetCommunication(cellchat), file.path(OUTPUT_DIR, "interactions_weak_KO_ALL_ECM.csv"), row.names=FALSE)
cat("Saved weak_KO_ALL_ECM RDS and CSV.\n")

# Build summary table
cat("\n========== Building summary ==========\n")
ko_list <- list(
  KO_COL1A1 = c('COL1A1'),
  KO_COL1A2 = c('COL1A2'),
  KO_COL1A1_COL1A2 = c('COL1A1','COL1A2'),
  KO_FN1 = c('FN1'),
  KO_ALL_ECM = c('COL1A1','COL1A2','FN1')
)

base_strong <- readRDS(file.path(OUTPUT_DIR, 'cellchat_strong.rds'))
base_weak <- readRDS(file.path(OUTPUT_DIR, 'cellchat_weak.rds'))
ft_strong_base <- subsetCommunication(base_strong)
ft_strong_base <- ft_strong_base[ft_strong_base$source=='Fibroblast' & ft_strong_base$target=='Treg',]
ft_weak_base <- subsetCommunication(base_weak)
ft_weak_base <- ft_weak_base[ft_weak_base$source=='Fibroblast' & ft_weak_base$target=='Treg',]

base_s <- sum(ft_strong_base$prob)
base_w <- sum(ft_weak_base$prob)

summary <- data.frame()
for (ko_name in names(ko_list)) {
  # Strong
  rds_s <- file.path(OUTPUT_DIR, paste0('cellchat_strong_', ko_name, '.rds'))
  if (file.exists(rds_s)) {
    cc_s <- readRDS(rds_s)
    ft_s <- subsetCommunication(cc_s)
    ft_s <- ft_s[ft_s$source=='Fibroblast' & ft_s$target=='Treg',]
    ko_s <- if(nrow(ft_s)>0) sum(ft_s$prob) else 0
    retained_s <- if(base_s>0) round(ko_s/base_s*100, 2) else NA
    lost_s <- if(base_s>0) round(100 - retained_s, 2) else NA
    summary <- rbind(summary, data.frame(
      KO=ko_name, Group='Strong',
      Base_Interactions=nrow(ft_strong_base), KO_Interactions=nrow(ft_s),
      Base_Prob=base_s, KO_Prob=ko_s,
      Retained_Pct=retained_s, Lost_Pct=lost_s
    ))
  }
  # Weak
  rds_w <- file.path(OUTPUT_DIR, paste0('cellchat_weak_', ko_name, '.rds'))
  if (file.exists(rds_w)) {
    cc_w <- readRDS(rds_w)
    ft_w <- subsetCommunication(cc_w)
    ft_w <- ft_w[ft_w$source=='Fibroblast' & ft_w$target=='Treg',]
    ko_w <- if(nrow(ft_w)>0) sum(ft_w$prob) else 0
    retained_w <- if(base_w>0) round(ko_w/base_w*100, 2) else NA
    lost_w <- if(base_w>0) round(100 - retained_w, 2) else NA
    summary <- rbind(summary, data.frame(
      KO=ko_name, Group='Weak',
      Base_Interactions=nrow(ft_weak_base), KO_Interactions=nrow(ft_w),
      Base_Prob=base_w, KO_Prob=ko_w,
      Retained_Pct=retained_w, Lost_Pct=lost_w
    ))
  }
}

cat("\n========== KO SUMMARY ==========\n")
print(summary, row.names=FALSE)
write.csv(summary, file.path(OUTPUT_DIR, 'ko_summary.csv'), row.names=FALSE)
cat('\nAll done. Saved ko_summary.csv\n')
