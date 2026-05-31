# ============================================================================
# scTenifoldKnk — New KOs with better targets
# CCR8+ Treg: KO HLA-DRB5 (MHC-II)
# MKI67+ Treg: KO CD74 (MIF receptor)
# ============================================================================

library(scTenifoldKnk)
library(Matrix)

cat("=== New scTenifoldKnk KOs ===\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Key genes for interpretation
MHC_GENES <- c("HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-DRB5",
               "HLA-DPA1","HLA-DQA1","HLA-DQB1","B2M","CD8A","CD8B","CD4")
TREG_GENES <- c("FOXP3","IL2RA","CTLA4","TIGIT","ICOS","IL2RB","CCR8","TNFRSF18")
MIF_GENES <- c("CD74","CD44","MIF","CXCR4","CXCR2")
STRESS_GENES <- c("HSPA1A","HSPA1B","HSP90AA1","DNAJB1","BAG3","HSF1","HSPH1")
PROLIF_GENES <- c("MKI67","TOP2A","PCNA","CDK1","CDK2")
KEY_GENES <- unique(c(MHC_GENES, TREG_GENES, MIF_GENES, STRESS_GENES, PROLIF_GENES))

# ======================== Helper ========================
run_ko <- function(data_dir, ko_gene, out_dir, label, params = "default") {
  cat(rep("=", 60), "\n", sep = "")
  cat(sprintf("[%s] KO: %s (%s params)\n", label, ko_gene, params))
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
  
  # Determine params based on cell count
  n_cells <- ncol(X)
  if (n_cells < 5000) {
    # Small dataset: reduce params
    nc_nNet <- 3
    nc_nCells <- min(200, n_cells)
    nc_q <- 0.8
    td_K <- 1
    ma_nDim <- 2
    cat(sprintf("  Using LOW-RES params (%d cells)\n", n_cells))
  } else {
    nc_nNet <- 5
    nc_nCells <- min(300, n_cells)
    nc_q <- 0.9
    td_K <- 2
    ma_nDim <- 3
    cat(sprintf("  Using DEFAULT params (%d cells)\n", n_cells))
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
      nc_nNet = nc_nNet,
      nc_nCells = nc_nCells,
      nc_nComp = 2,
      nc_scaleScores = TRUE,
      nc_symmetric = FALSE,
      nc_q = nc_q,
      td_K = td_K,
      td_maxIter = 1000,
      td_maxError = 1e-05,
      td_nDecimal = 3,
      ma_nDim = ma_nDim
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
    cat(sprintf("  Perturbed genes: %d | Significant (BH<0.05): %d\n", nrow(diff_reg), n_sig))
    
    # Key hits (our curated gene sets)
    key_hits <- diff_reg[diff_reg$gene %in% KEY_GENES & diff_reg$p.adj < 0.05, ]
    if (nrow(key_hits) > 0) {
      cat(sprintf("  Key hits (%d):\n", nrow(key_hits)))
      print(key_hits[, c("gene", "Z", "p.value", "p.adj")], row.names = FALSE)
    } else {
      cat("  No significant key hits in curated gene sets.\n")
    }
    
    # Print ALL significant genes for inspection
    sig_all <- diff_reg[diff_reg$p.adj < 0.05 & diff_reg$gene != ko_gene, ]
    if (nrow(sig_all) > 0) {
      cat(sprintf("  All significant (excluding %s itself):\n", ko_gene))
      print(sig_all[order(-abs(sig_all$Z)), c("gene", "Z", "p.adj")], row.names = FALSE)
    }
    
    out_file <- file.path(out_dir, sprintf("ko_%s_%s.csv", label, ko_gene))
    write.csv(diff_reg, out_file, row.names = FALSE)
    cat(sprintf("  Saved: %s\n", out_file))
    return(diff_reg)
  }
  
  return(NULL)
}

# ======================== Run KOs ========================

# 1. HLA-DRB5 KO in CCR8+ Treg (24,492 cells — use default params)
res1 <- run_ko(
  data_dir = "D:/Research/tomato/data/treg_ccr8_for_ko",
  ko_gene = "HLA-DRB5",
  out_dir = "D:/Research/tomato/results/scTenifoldKnk_new_kos",
  label = "CCR8"
)

# 2. CD74 KO in MKI67+ Treg (2,323 cells — use low-res params)
res2 <- run_ko(
  data_dir = "D:/Research/tomato/data/treg_mki67_for_ko",
  ko_gene = "CD74",
  out_dir = "D:/Research/tomato/results/scTenifoldKnk_new_kos",
  label = "MKI67"
)

# ======================== Summary ========================
cat("\n", rep("=", 60), "\n", sep = "")
cat("KO SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")

if (!is.null(res1)) {
  sig1 <- res1[res1$p.adj < 0.05 & res1$gene != "HLA-DRB5", ]
  cat(sprintf("CCR8+ KO HLA-DRB5: %d sig genes (BH<0.05, excl. self)\n", nrow(sig1)))
  if (nrow(sig1) > 0) {
    top1 <- sig1[order(-abs(sig1$Z)), ]
    cat("  Top 10:\n")
    print(top1[1:min(10, nrow(top1)), c("gene", "Z", "p.adj")], row.names = FALSE)
  }
}

if (!is.null(res2)) {
  sig2 <- res2[res2$p.adj < 0.05 & res2$gene != "CD74", ]
  cat(sprintf("MKI67+ KO CD74: %d sig genes (BH<0.05, excl. self)\n", nrow(sig2)))
  if (nrow(sig2) > 0) {
    top2 <- sig2[order(-abs(sig2$Z)), ]
    cat("  Top 10:\n")
    print(top2[1:min(10, nrow(top2)), c("gene", "Z", "p.adj")], row.names = FALSE)
  }
}

cat("\nDONE!\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
