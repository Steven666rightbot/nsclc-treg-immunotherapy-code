library(scTenifoldKnk)

set.seed(42)

OUTPUT_DIR <- "D:/Research/tomato/results/scTenifoldKnk_ko_merged"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(OUTPUT_DIR, "log_test.txt")
sink(log_file, split = TRUE, append = FALSE)

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Starting minimal test\n")

treg_all <- readRDS("D:/Research/potato/data/roosted/treg_207422.rds")
cat("Loaded:", nrow(treg_all), "genes x", ncol(treg_all), "cells\n")

# Minimal preprocessing: keep top 1000 most variable genes + CD44
expr_var <- apply(treg_all, 1, var)
top_genes <- names(sort(expr_var, decreasing = TRUE))[1:1000]
if (!("CD44" %in% top_genes)) top_genes <- c(top_genes, "CD44")
mat <- treg_all[top_genes, ]
mat <- log1p(mat)
cat("Filtered to", nrow(mat), "genes x", ncol(mat), "cells\n")

# Ultra-minimal parameters
cat("Running scTenifoldKnk with nNet=1, nCells=50, K=1...\n")
start <- Sys.time()

result <- tryCatch({
  scTenifoldKnk(
    countMatrix = mat,
    gKO = "CD44",
    qc_mtThreshold = 0.1,
    qc_minLSize = 100,
    nc_lambda = 0,
    nc_nNet = 1,
    nc_nCells = 50,
    nc_nComp = 2,
    nc_scaleScores = TRUE,
    nc_symmetric = FALSE,
    nc_q = 0.9,
    td_K = 1,
    td_maxIter = 100,
    td_maxError = 1e-03,
    td_nDecimal = 3,
    ma_nDim = 2
  )
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  return(NULL)
})

elapsed <- difftime(Sys.time(), start, units = "mins")
cat("Finished in", as.numeric(elapsed), "minutes\n")

if (!is.null(result)) {
  cat("Perturbed genes:", nrow(result$diffRegulation), "\n")
  cat("Top hits:\n")
  print(head(result$diffRegulation[, c("gene", "Z", "p.value", "p.adj")], 10))
  write.csv(result$diffRegulation, file.path(OUTPUT_DIR, "test_cd44_minimal.csv"), row.names = FALSE)
} else {
  cat("Result is NULL\n")
}

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Test complete\n")
sink()
