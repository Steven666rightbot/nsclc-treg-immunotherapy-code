# =============================================================================
# Slingshot pseudotime analysis for CCR8+ and MKI67+ Treg
# =============================================================================

library(Seurat)
library(SeuratObject)
library(slingshot)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)

# Paths
CCR8_DIR <- "D:/Research/tomato/data/treg_ccr8_for_ko"
MKI67_DIR <- "D:/Research/tomato/data/treg_mki67_for_ko"
OUT_DIR <- "D:/Research/tomato/results/slingshot_ccr8_mki67"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ======================== Helper: read 10x-like mtx ========================
read_mtx <- function(dir_path, label) {
  cat(sprintf("Reading %s from %s...\n", label, dir_path))
  expr <- readMM(file.path(dir_path, "matrix.mtx"))
  genes <- readLines(file.path(dir_path, "genes.tsv"))
  barcodes <- readLines(file.path(dir_path, "barcodes.tsv"))
  
  # mtx is genes x cells? Check dimensions
  cat(sprintf("  Raw: %d x %d\n", nrow(expr), ncol(expr)))
  
  # If genes < cells, likely cells x genes, need transpose
  if (nrow(expr) < ncol(expr)) {
    expr <- t(expr)
    cat(sprintf("  Transposed to: %d x %d\n", nrow(expr), ncol(expr)))
  }
  
  rownames(expr) <- genes
  colnames(expr) <- barcodes
  expr <- as(expr, "dgCMatrix")
  
  # Create minimal Seurat object
  obj <- CreateSeuratObject(counts = expr, project = label)
  obj$subtype <- label
  return(obj)
}

# ======================== Read data ========================
ccr8 <- read_mtx(CCR8_DIR, "CCR8")
mki67 <- read_mtx(MKI67_DIR, "MKI67")

# Merge
cat("\nMerging CCR8 + MKI67...\n")
combined <- merge(ccr8, mki67)
cat(sprintf("Combined: %d cells\n", ncol(combined)))

# ======================== Preprocessing ========================
cat("\nPreprocessing...\n")
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.8)
combined <- RunUMAP(combined, dims = 1:20)

# Save basic UMAP
ggsave(file.path(OUT_DIR, "umap_by_subtype.png"), 
       DimPlot(combined, reduction = "umap", group.by = "subtype"),
       width = 8, height = 6)
ggsave(file.path(OUT_DIR, "umap_by_cluster.png"),
       DimPlot(combined, reduction = "umap", group.by = "seurat_clusters"),
       width = 8, height = 6)

# ======================== Slingshot ========================
cat("\nRunning Slingshot...\n")

# Convert to SingleCellExperiment (Seurat v5 compatible)
counts_matrix <- GetAssayData(combined, layer = "counts")
logcounts_matrix <- GetAssayData(combined, layer = "data")

sce <- SingleCellExperiment(
  assays = list(counts = counts_matrix, logcounts = logcounts_matrix),
  colData = combined@meta.data
)

# Add reduced dims
reducedDim(sce, "PCA") <- Embeddings(combined, "pca")
reducedDim(sce, "UMAP") <- Embeddings(combined, "umap")

# Run slingshot using UMAP and clusters
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP')

# Extract pseudotime
pseudotime <- slingPseudotime(sce)
lineages <- slingLineages(sce)
curves <- slingCurves(sce)

cat("Lineages:\n")
print(lineages)

cat("\nPseudotime summary:\n")
print(summary(pseudotime))

# Save pseudotime
pt_df <- as.data.frame(pseudotime)
pt_df$cell <- rownames(pt_df)
pt_df$subtype <- combined$subtype
pt_df$cluster <- combined$seurat_clusters
write.csv(pt_df, file.path(OUT_DIR, "pseudotime.csv"), row.names = FALSE)

# ======================== Visualization ========================
cat("\nPlotting...\n")

# Color by pseudotime
umap_coords <- as.data.frame(Embeddings(combined, "umap"))
umap_coords$pseudotime <- pseudotime[, 1]
umap_coords$subtype <- combined$subtype

p1 <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_gradient(low = "gray90", high = "darkblue") +
  theme_minimal() +
  labs(title = "Slingshot Pseudotime", color = "Pseudotime")
ggsave(file.path(OUT_DIR, "slingshot_pseudotime.png"), p1, width = 8, height = 6)

# Split by subtype
p2 <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_gradient(low = "gray90", high = "darkblue") +
  facet_wrap(~subtype) +
  theme_minimal() +
  labs(title = "Pseudotime by Subtype", color = "Pseudotime")
ggsave(file.path(OUT_DIR, "slingshot_pseudotime_by_subtype.png"), p2, width = 10, height = 5)

# Save SingleCellExperiment
saveRDS(sce, file.path(OUT_DIR, "slingshot_sce.rds"))

# Save Seurat
saveRDS(combined, file.path(OUT_DIR, "slingshot_seurat.rds"))

cat(sprintf("\nAll results saved to: %s\n", OUT_DIR))
cat("Done!\n")
