# Extract Treg expression matrix from Monocle3 CDS for SCENIC analysis
library(monocle3)

print("Loading CDS...")
cds <- readRDS("data/monocle3_treg/treg_monocle3_cds.rds")
print(paste("CDS dimensions:", nrow(cds), "genes x", ncol(cds), "cells"))

# Check cell types
print("Cell type table:")
print(table(colData(cds)$sub_cell_type))

# Extract expression matrix (sparse matrix to dense for export)
print("Extracting count matrix...")
counts <- as.matrix(assay(cds, "counts"))
print(paste("Count matrix:", nrow(counts), "x", ncol(counts)))

# Save as CSV (may be very large - consider downsampling)
# For SCENIC, we need cells x genes format
# Let's subset to highly variable genes or TF list

# Get TF list from pySCENIC resource (or use known human TFs)
# For now, save the full matrix with cell annotations
meta <- as.data.frame(colData(cds))
write.csv(meta, "results/scenic_input/treg_metadata.csv", row.names = TRUE)

# Save count matrix as loom or CSV
# loom is preferred for pySCENIC
print("Saving to loom...")

# Check if loomR is available
if (!requireNamespace("loomR", quietly = TRUE)) {
  print("loomR not available, saving as sparse CSV instead...")
  # Save as sparse matrix in Matrix Market format
  library(Matrix)
  counts_sparse <- assay(cds, "counts")
  writeMM(counts_sparse, "results/scenic_input/treg_counts.mtx")
  writeLines(rownames(counts_sparse), "results/scenic_input/genes.txt")
  writeLines(colnames(counts_sparse), "results/scenic_input/cells.txt")
} else {
  # Create loom file
  library(loomR)
  lfile <- create(filename = "results/scenic_input/treg_counts.loom",
                  data = counts,
                  cell.attrs = meta)
  lfile$close_all()
}

print("Done!")
