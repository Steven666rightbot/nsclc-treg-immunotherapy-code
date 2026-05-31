library(Matrix)

# Check CCR8+ Treg
cat("=== CCR8+ Treg (24,492 cells) ===\n")
mat1 <- readMM("D:/Research/tomato/data/treg_ccr8_for_ko/matrix.mtx")
genes1 <- readLines("D:/Research/tomato/data/treg_ccr8_for_ko/genes.tsv")
rownames(mat1) <- genes1
cat("Matrix:", nrow(mat1), "genes x", ncol(mat1), "cells\n")
cat("CXCL13 mean:", mean(mat1["CXCL13", ]), "| nonzero:", sum(mat1["CXCL13", ] > 0) / ncol(mat1), "\n")
cat("B2M mean:", mean(mat1["B2M", ]), "| nonzero:", sum(mat1["B2M", ] > 0) / ncol(mat1), "\n")

# Check MKI67+ Treg
cat("\n=== MKI67+ Treg (2,323 cells) ===\n")
mat2 <- readMM("D:/Research/tomato/data/treg_mki67_for_ko/matrix.mtx")
genes2 <- readLines("D:/Research/tomato/data/treg_mki67_for_ko/genes.tsv")
rownames(mat2) <- genes2
cat("Matrix:", nrow(mat2), "genes x", ncol(mat2), "cells\n")
cat("CXCL13 mean:", mean(mat2["CXCL13", ]), "| nonzero:", sum(mat2["CXCL13", ] > 0) / ncol(mat2), "\n")
cat("B2M mean:", mean(mat2["B2M", ]), "| nonzero:", sum(mat2["B2M", ] > 0) / ncol(mat2), "\n")
