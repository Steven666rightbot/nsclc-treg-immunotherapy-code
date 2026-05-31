# ============================================================================
# scTenifoldKnk — MKI67+ Treg KO B2M (low-resolution params for small n)
# 2323 cells → QC ~776; reduce nNet/nCells to avoid rank deficiency
# ============================================================================

library(scTenifoldKnk)
library(Matrix)

cat("=== MKI67+ Treg KO B2M (low-res) ===\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Key genes
MHC_GENES <- c("HLA-A","HLA-B","HLA-C","HLA-E","B2M","CD8A","CD8B")
TREG_FUNCTION <- c("FOXP3","IL2RA","CTLA4","TIGIT","ICOS","IL2RB","MKI67","TOP2A")
KEY_GENES <- unique(c(MHC_GENES, TREG_FUNCTION))

# Load
data_dir <- "D:/Research/tomato/data/treg_mki67_for_ko"
out_dir <- "D:/Research/tomato/results/scTenifoldKnk_ko_mki67_b2m"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading...\n")
X <- readMM(file.path(data_dir, "matrix.mtx"))
genes <- readLines(file.path(data_dir, "genes.tsv"))
barcodes <- readLines(file.path(data_dir, "barcodes.tsv"))
rownames(X) <- genes; colnames(X) <- barcodes
cat(sprintf("  Raw: %d x %d\n", nrow(X), ncol(X)))

expr_level <- mean(X["B2M", ])
nz_frac <- sum(X["B2M", ] > 0) / ncol(X)
cat(sprintf("  B2M mean: %.4f | nonzero: %.2f%%\n", expr_level, nz_frac*100))

# HVG
cat("HVG selection...\n")
expr_var <- apply(X, 1, var)
top2000 <- names(sort(expr_var, decreasing = TRUE))[1:2000]
missing <- c("B2M", KEY_GENES)[!(c("B2M", KEY_GENES) %in% top2000)]
if (length(missing) > 0) {
  cat(sprintf("  Adding %d genes\n", length(missing)))
  top2000 <- c(top2000, missing)
}
mat <- X[top2000, ]
mat <- log1p(mat)
cat(sprintf("  Final: %d x %d\n", nrow(mat), ncol(mat)))

# Run with reduced params for small sample
cat("Running scTenifoldKnk (low-res params)...\n")
start_time <- Sys.time()

result <- tryCatch({
  scTenifoldKnk(
    countMatrix = mat,
    gKO = "B2M",
    qc_mtThreshold = 0.1,
    qc_minLSize = 1000,
    nc_lambda = 0,
    nc_nNet = 3,        # reduced from 5
    nc_nCells = 200,    # reduced from 300
    nc_nComp = 2,
    nc_scaleScores = TRUE,
    nc_symmetric = FALSE,
    nc_q = 0.8,         # reduced from 0.9
    td_K = 1,           # reduced from 2
    td_maxIter = 1000,
    td_maxError = 1e-05,
    td_nDecimal = 3,
    ma_nDim = 2         # reduced from 3
  )
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
  return(NULL)
})

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("  Finished in %.1f min\n", as.numeric(elapsed)))

if (!is.null(result)) {
  diff_reg <- result$diffRegulation
  diff_reg$KO_Gene <- "B2M"
  diff_reg$CellType <- "MKI67"
  n_sig <- sum(diff_reg$p.adj < 0.05, na.rm = TRUE)
  cat(sprintf("  Perturbed: %d | Significant: %d\n", nrow(diff_reg), n_sig))
  
  key_hits <- diff_reg[diff_reg$gene %in% KEY_GENES & diff_reg$p.adj < 0.05, ]
  if (nrow(key_hits) > 0) {
    cat(sprintf("  Key hits (%d):\n", nrow(key_hits)))
    print(key_hits[, c("gene", "Z", "p.value", "p.adj")], row.names = FALSE)
  } else {
    cat("  No significant key hits.\n")
  }
  
  out_file <- file.path(out_dir, "ko_MKI67_B2M_lowres.csv")
  write.csv(diff_reg, out_file, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", out_file))
}

cat("\nDONE!\n")
