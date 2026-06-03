## Export key gene expression + pseudotime for Python plotting (Fig 6D)
library(monocle3)
library(Matrix)

CDS_PATH <- "D:/Research/tomato/results/monocle3_ccr8_mki67/cds_ccr8_mki67_full.rds"
PT_FILE <- "D:/Research/cucumber/fig6/data/slingshot_pseudotime_umap.csv"
OUT_DIR <- "D:/Research/cucumber/fig6/data"

cds <- readRDS(CDS_PATH)
pt_df <- read.csv(PT_FILE, stringsAsFactors = FALSE)

# Key genes to export
key_genes <- c("CCR8", "MKI67", "FOXP3", "CTLA4", "IL2RA", "TIGIT",
               "CD74", "HLA-DRB5", "HLA-A", "B2M",
               "ACTA2", "MYL9", "VCL", "TLN1",
               "CXCL13", "GZMB", "PRF1", "TGFB1", "IL10",
               "HSF1", "CCL3", "CCL4")

# Check which are available
gene_names <- rowData(cds)$gene_short_name
available <- key_genes[key_genes %in% gene_names]
cat("Available genes:", length(available), "/", length(key_genes), "\n")
cat(paste(available, collapse=", "), "\n")

# Get counts and normalize
counts_mat <- counts(cds)
size_factors <- colData(cds)$Size_Factor
if (is.null(size_factors)) {
  size_factors <- colSums(counts_mat) / median(colSums(counts_mat))
}

# Extract expression for key genes only (sparse)
gene_indices <- which(gene_names %in% available)
expr_mat <- as.matrix(counts_mat[gene_indices, , drop=FALSE])
rownames(expr_mat) <- gene_names[gene_indices]

# Normalize: log10(CPM + 1)
norm_mat <- t(t(expr_mat) / size_factors)
norm_mat <- log10(norm_mat + 1)

# Build output dataframe
gene_list <- list()
for (g in available) {
  idx <- which(rownames(norm_mat) == g)
  df_g <- data.frame(
    cell_id = colnames(cds),
    gene = g,
    expression = norm_mat[idx, ],
    stringsAsFactors = FALSE
  )
  gene_list[[g]] <- df_g
}
gene_df <- do.call(rbind, gene_list)

# Merge with pseudotime
gene_df <- merge(gene_df, pt_df[, c("cell_id", "pseudotime", "sub_cell_type", 
                                     "response_binary", "primary_lineage")],
                 by = "cell_id", all.x = TRUE)

# Keep only main trajectory cells (lineage 1 or 2) with valid pseudotime
gene_df <- gene_df[gene_df$primary_lineage %in% 1:2 & !is.na(gene_df$pseudotime), ]

write.csv(gene_df, file.path(OUT_DIR, "slingshot_gene_expr.csv"), row.names = FALSE)
cat("\nExported:", nrow(gene_df), "rows for", length(available), "genes\n")
cat("Done.\n")
