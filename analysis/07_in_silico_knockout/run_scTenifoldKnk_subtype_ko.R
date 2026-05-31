# ============================================================================
# scTenifoldKnk Virtual Knockout — Subtype-specific Treg KO
# CCR8+ Treg: KO CXCL13 | MKI67+ Treg: KO B2M
# ============================================================================

library(scTenifoldKnk)
library(Matrix)

cat("=== scTenifoldKnk Subtype-specific KO ===\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Key gene sets for interpretation
CXCL_GENES <- c("CXCL13","CXCL9","CXCL10","CXCL11","CXCL12","CXCR1","CXCR2","CXCR3","CXCR4","CXCR5","CXCR6")
MHC_GENES <- c("HLA-A","HLA-B","HLA-C","HLA-E","B2M","CD8A","CD8B")
TREG_FUNCTION <- c("FOXP3","IL2RA","CTLA4","TIGIT","ICOS","IL2RB")
INTEGRIN_DOWNSTREAM <- c("RHOA","ROCK1","ROCK2","FYN","LCK","PIK3CD","AKT1","MAPK1","PTPN11")
KEY_GENES <- unique(c(CXCL_GENES, MHC_GENES, TREG_FUNCTION, INTEGRIN_DOWNSTREAM))

# ======================== Helper ========================
run_ko <- function(data_dir, ko_gene, out_dir, label) {
  cat(rep("=", 60), "\n", sep = "")
  cat(sprintf("[%s] KO: %s\n", label, ko_gene))
  cat(sprintf("Data: %s\n", data_dir))
  cat(rep("=", 60), "\n", sep = "")
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load raw counts
  cat("Loading expression matrix...\n")
  X <- readMM(file.path(data_dir, "matrix.mtx"))
  genes <- readLines(file.path(data_dir, "genes.tsv"))
  barcodes <- readLines(file.path(data_dir, "barcodes.tsv"))
  rownames(X) <- genes
  colnames(X) <- barcodes
  cat(sprintf("  Raw: %d genes x %d cells\n", nrow(X), ncol(X)))
  
  # Check target gene
  if (!(ko_gene %in% rownames(X))) {
    cat(sprintf("ERROR: %s not found in matrix\n", ko_gene))
    return(NULL)
  }
  expr_level <- mean(X[ko_gene, ])
  nz_frac <- sum(X[ko_gene, ] > 0) / ncol(X)
  cat(sprintf("  %s mean: %.4f | nonzero: %.2f%%\n", ko_gene, expr_level, nz_frac * 100))
  
  if (expr_level < 0.01) {
    cat(sprintf("WARNING: %s expression too low, skipping\n", ko_gene))
    return(NULL)
  }
  
  # Select top 2000 HVG
  cat("Selecting top 2000 HVG...\n")
  expr_var <- apply(X, 1, var)
  top2000 <- names(sort(expr_var, decreasing = TRUE))[1:2000]
  
  # Force-keep key genes and KO target
  missing <- c(ko_gene, KEY_GENES)[!(c(ko_gene, KEY_GENES) %in% top2000)]
  if (length(missing) > 0) {
    cat(sprintf("  Adding %d genes to HVG set\n", length(missing)))
    top2000 <- c(top2000, missing)
  }
  
  mat <- X[top2000, ]
  mat <- log1p(mat)
  cat(sprintf("  Final: %d genes x %d cells\n", nrow(mat), ncol(mat)))
  
  # Run scTenifoldKnk
  cat("Running scTenifoldKnk...\n")
  start_time <- Sys.time()
  
  result <- tryCatch({
    scTenifoldKnk(
      countMatrix = mat,
      gKO = ko_gene,
      qc_mtThreshold = 0.1,
      qc_minLSize = 1000,
      nc_lambda = 0,
      nc_nNet = 5,
      nc_nCells = 300,
      nc_nComp = 2,
      nc_scaleScores = TRUE,
      nc_symmetric = FALSE,
      nc_q = 0.9,
      td_K = 2,
      td_maxIter = 1000,
      td_maxError = 1e-05,
      td_nDecimal = 3,
      ma_nDim = 3
    )
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("  Finished in %.1f min\n", as.numeric(elapsed)))
  
  if (!is.null(result)) {
    diff_reg <- result$diffRegulation
    diff_reg$KO_Gene <- ko_gene
    diff_reg$CellType <- label
    
    n_sig <- sum(diff_reg$p.adj < 0.05, na.rm = TRUE)
    cat(sprintf("  Perturbed genes: %d | Significant: %d\n", nrow(diff_reg), n_sig))
    
    # Key hits
    key_hits <- diff_reg[diff_reg$gene %in% KEY_GENES & diff_reg$p.adj < 0.05, ]
    if (nrow(key_hits) > 0) {
      cat(sprintf("  Key hits (%d):\n", nrow(key_hits)))
      print(key_hits[, c("gene", "Z", "p.value", "p.adj")], row.names = FALSE)
    }
    
    out_file <- file.path(out_dir, sprintf("ko_%s_%s.csv", label, ko_gene))
    write.csv(diff_reg, out_file, row.names = FALSE)
    cat(sprintf("  Saved: %s\n", out_file))
    return(diff_reg)
  }
  
  return(NULL)
}

# ======================== Run KOs ========================
res1 <- run_ko(
  data_dir = "D:/Research/tomato/data/treg_ccr8_for_ko",
  ko_gene = "CXCL13",
  out_dir = "D:/Research/tomato/results/scTenifoldKnk_ko_ccr8_cxcl13",
  label = "CCR8"
)

res2 <- run_ko(
  data_dir = "D:/Research/tomato/data/treg_mki67_for_ko",
  ko_gene = "B2M",
  out_dir = "D:/Research/tomato/results/scTenifoldKnk_ko_mki67_b2m",
  label = "MKI67"
)

# ======================== Summary ========================
cat("\n", rep("=", 60), "\n", sep = "")
cat("KO SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")

if (!is.null(res1)) {
  sig1 <- res1[res1$p.adj < 0.05, ]
  cat(sprintf("CCR8+ KO CXCL13: %d sig / %d total\n", nrow(sig1), nrow(res1)))
}
if (!is.null(res2)) {
  sig2 <- res2[res2$p.adj < 0.05, ]
  cat(sprintf("MKI67+ KO B2M: %d sig / %d total\n", nrow(sig2), nrow(res2)))
}

cat("\nDONE!\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
