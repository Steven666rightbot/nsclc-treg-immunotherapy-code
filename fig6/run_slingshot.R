## =============================================================================
## Figure 6 — Slingshot Pseudotime Analysis for CCR8+ / MKI67+ Treg
##
## Pipeline:
##   1. Load Monocle3 CDS (preprocessed, has UMAP + clusters)
##   2. Normalize expression, build SCE
##   3. Run Slingshot (MST + principal curves)
##   4. Export pseudotime + UMAP + curve coordinates to CSV
## =============================================================================

library(slingshot)
library(SingleCellExperiment)
library(monocle3)
library(Matrix)

# Paths
CDS_PATH <- "D:/Research/tomato/results/monocle3_ccr8_mki67/cds_ccr8_mki67_full.rds"
OUT_DIR  <- "D:/Research/cucumber/fig6/data"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat("Figure 6 — Slingshot Pseudotime Analysis\n")
cat("========================================\n")

# ======================== 1. Load CDS ========================
cat("\n[1] Loading Monocle3 CDS...\n")
cds <- readRDS(CDS_PATH)
cat("  Cells:", ncol(cds), "x Genes:", nrow(cds), "\n")
print(table(colData(cds)$sub_cell_type))
print(table(colData(cds)$response_binary))

# ======================== 2. Normalize + Build SCE ========================
cat("\n[2] Normalizing expression...\n")

counts_mat <- counts(cds)
cat("  counts dim:", nrow(counts_mat), "x", ncol(counts_mat), "\n")

# Normalize: log(CPM/10 + 1) — Monocle3 default
size_factors <- colData(cds)$Size_Factor
if (is.null(size_factors)) {
  size_factors <- colSums(counts_mat) / median(colSums(counts_mat))
}
logcounts_mat <- t(t(counts_mat) / size_factors)
logcounts_mat <- log10(logcounts_mat + 1)
logcounts_mat <- as(logcounts_mat, "dgCMatrix")
cat("  logcounts computed\n")

# UMAP + PCA
umap_coords <- reducedDims(cds)$UMAP
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
pca_coords <- reducedDims(cds)$PCA
cat("  UMAP:", nrow(umap_coords), "cells\n")

# Metadata
meta <- as.data.frame(colData(cds))
meta$cluster_id <- as.character(clusters(cds))
cat("  Clusters:", length(unique(meta$cluster_id)), "\n")
print(table(meta$cluster_id))

# Build SCE
sce <- SingleCellExperiment(
  assays = list(counts = counts_mat, logcounts = logcounts_mat),
  colData = meta
)
reducedDims(sce) <- SimpleList(PCA = pca_coords, UMAP = umap_coords)

# ======================== 3. Run Slingshot ========================
cat("\n[3] Running Slingshot...\n")

# Use UMAP + cluster labels
# Main argument tweaks for clean linear trajectory:
#   stretch = 0 — don't extrapolate beyond data range
#   omega = TRUE — allow soft assignment to multiple lineages
sce <- slingshot(sce, clusterLabels = "cluster_id", reducedDim = "UMAP",
                 start.clus = NULL,  # auto-detect start cluster
                 stretch = 0, omega = TRUE)

# Results
pst <- slingPseudotime(sce)
curves <- slingCurves(sce)
lineages <- slingLineages(sce)

cat("\n  Lineages:\n")
print(lineages)
cat("\n  Pseudotime range:", range(pst[,1], na.rm=TRUE), "\n")
cat("  Curves:", length(curves), "\n")

# ======================== 4. Export ========================
cat("\n[4] Exporting to CSV...\n")

# Cell-level data
df_cells <- data.frame(
  cell_id          = colnames(sce),
  UMAP_1           = umap_coords[, 1],
  UMAP_2           = umap_coords[, 2],
  PCA_1            = pca_coords[, 1],
  PCA_2            = pca_coords[, 2],
  pseudotime       = pst[, 1],
  sub_cell_type    = sce$sub_cell_type,
  response_binary  = sce$response_binary,
  cluster_id       = sce$cluster_id,
  stringsAsFactors = FALSE
)

# Handle multiple lineages
if (ncol(pst) > 1) {
  df_cells$pseudotime_L2 <- pst[, 2]
  df_cells$primary_lineage <- apply(pst, 1, which.max)
} else {
  df_cells$primary_lineage <- 1
}

write.csv(df_cells, file.path(OUT_DIR, "slingshot_pseudotime_umap.csv"), row.names = FALSE)
cat("  Cells CSV:", nrow(df_cells), "rows\n")

# Curve coordinates (the fitted principal curve, ~150 points per curve)
curve_list <- list()
for (i in seq_along(curves)) {
  pts <- data.frame(
    curve = rep(i, nrow(curves[[i]]$s)),
    x = curves[[i]]$s[, 1],
    y = curves[[i]]$s[, 2],
    pseudotime = curves[[i]]$lambda[curves[[i]]$ord],
    stringsAsFactors = FALSE
  )
  curve_list[[i]] <- pts
}
curves_df <- do.call(rbind, curve_list)
write.csv(curves_df, file.path(OUT_DIR, "slingshot_curves.csv"), row.names = FALSE)
cat("  Curves CSV:", nrow(curves_df), "points\n")

# Save SCE for potential re-use
saveRDS(sce, file.path(OUT_DIR, "slingshot_sce.rds"))
cat("  SCE saved\n")

cat("\n=== DONE ===\n")
cat("Output:", OUT_DIR, "\n")
