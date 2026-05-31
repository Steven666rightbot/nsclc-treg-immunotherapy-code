library(CellChat)
library(Seurat)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
SEED <- args[1]

DATA_DIR <- file.path("D:/Research/tomato/go_exploration", paste0("cellchat_input_subtype_seed", SEED))
OUT_DIR  <- file.path("D:/Research/tomato/go_exploration", paste0("cellchat_output_subtype_seed", SEED))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== CellChat for seed %s ===\n", SEED))

X <- Matrix::readMM(file.path(DATA_DIR, "matrix.mtx"))
barcodes <- readLines(file.path(DATA_DIR, "barcodes.tsv"))
genes <- readLines(file.path(DATA_DIR, "genes.tsv"))

dim(X)
cat(sprintf("Matrix: %d genes x %d cells\n", nrow(X), ncol(X)))

meta <- read.csv(file.path(DATA_DIR, "cell_annotation.csv"))
rownames(meta) <- meta$cellID

common <- intersect(barcodes, meta$cellID)
cat(sprintf("Common barcodes: %d\n", length(common)))

idx <- match(common, barcodes)
X <- X[, idx]
colnames(X) <- common
rownames(X) <- genes

meta <- meta[common, ]
seu <- CreateSeuratObject(counts = X, meta.data = meta)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

run_cellchat <- function(seu_obj, group_name, group_cells) {
  cat(sprintf("\n===== Running CellChat for %s (%d cells) =====\n", group_name, length(group_cells)))
  sub <- subset(seu_obj, cells = group_cells)
  data.input <- GetAssayData(sub, layer = "data")
  meta.sub <- sub@meta.data
  cellchat <- createCellChat(object = data.input, meta = meta.sub, group.by = "cell_type")
  cellchatDB <- CellChatDB.human
  cellchat@DB <- cellchatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  saveRDS(cellchat, file.path(OUT_DIR, sprintf("cellchat_%s.rds", group_name)))
  cat(sprintf("Saved: cellchat_%s.rds\n", group_name))
  return(cellchat)
}

strong_cells <- colnames(seu)[seu$response_group == "Strong"]
weak_cells <- colnames(seu)[seu$response_group == "Weak"]

cat(sprintf("Strong cells: %d\n", length(strong_cells)))
cat(sprintf("Weak cells: %d\n", length(weak_cells)))

cc_strong <- run_cellchat(seu, "Strong", strong_cells)
cc_weak <- run_cellchat(seu, "Weak", weak_cells)

cat("\n===== Merging =====\n")
object.list <- list(Strong = cc_strong, Weak = cc_weak)
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat.merged, file.path(OUT_DIR, "cellchat_merged.rds"))

for (name in names(object.list)) {
  net <- object.list[[name]]@net$count
  cat(sprintf("  %s: %d total interactions\n", name, sum(net)))
}

pathway_comparison <- rankNet(cellchat.merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
ggsave(file.path(OUT_DIR, "pathway_comparison.pdf"), pathway_comparison, width = 10, height = 12)
ggsave(file.path(OUT_DIR, "pathway_comparison.png"), pathway_comparison, width = 10, height = 12, dpi = 300)

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
    pathway = pw, prob_strong = prob_strong, prob_weak = prob_weak,
    fold_weak_vs_strong = fold, present_strong = pw %in% cc_strong@netP$pathways,
    present_weak = pw %in% cc_weak@netP$pathways
  ))
}
pathway_df <- pathway_df %>% arrange(desc(fold_weak_vs_strong))
write.csv(pathway_df, file.path(OUT_DIR, "pathway_comparison.csv"), row.names = FALSE)

cat("\n===== Treg-specific (pathway-level) =====\n")
extract_treg_pathways <- function(cc, group_name) {
  groups <- levels(cc@idents)
  treg_groups <- c("Treg_FOXP3", "Treg_CCR8", "Treg_MKI67")
  pathways <- cc@netP$pathways
  out <- list()
  for (si in 1:length(groups)) {
    for (ti in 1:length(groups)) {
      if (groups[si] %in% treg_groups) {
        for (pi in 1:length(pathways)) {
          prob <- cc@netP$prob[si, ti, pi]
          if (!is.na(prob) && prob > 0) {
            out[[length(out) + 1]] <- data.frame(
              group = group_name, source = groups[si], target = groups[ti],
              pathway = pathways[pi], prob = prob, stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  if (length(out) == 0) return(data.frame())
  do.call(rbind, out)
}

treg_strong <- extract_treg_pathways(cc_strong, "Strong")
treg_weak <- extract_treg_pathways(cc_weak, "Weak")
treg_all <- rbind(treg_strong, treg_weak)
cat(sprintf("Total Treg interactions: %d (S:%d W:%d)\n", nrow(treg_all), nrow(treg_strong), nrow(treg_weak)))
write.csv(treg_all, file.path(OUT_DIR, "treg_interactions.csv"), row.names = FALSE)

if (nrow(treg_all) > 0) {
  treg_pw <- treg_all %>%
    group_by(pathway, source, target, group) %>%
    summarise(total_prob = sum(prob), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = total_prob, values_fill = 0) %>%
    mutate(diff = Weak - Strong, fold = ifelse(Strong > 0, Weak / Strong, ifelse(Weak > 0, Inf, 1))) %>%
    arrange(desc(diff))
  write.csv(treg_pw, file.path(OUT_DIR, "treg_pathway_comparison.csv"), row.names = FALSE)
  cat("Saved Treg pathway comparison\n")
}

cat("\n===== DONE =====\n")
