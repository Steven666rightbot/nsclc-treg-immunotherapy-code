# ================================================================
# scTenifoldKnk Virtual TF Knockout for Mechanosensitive TFs
# ================================================================
# Target TFs: MEF2C, TEAD1, NFKB1, STAT2, IRF8
# Strategy: filter genes with min_cells=500, then FORCE-KEEP target TFs
# so they can be knocked out even if lowly expressed.
# ================================================================

library(scTenifoldKnk)
library(data.table)
library(dplyr)

set.seed(42)

OUTPUT_DIR <- "D:/Research/tomato/results/scTenifoldKnk_ko"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_TFS <- c("MEF2C", "TEAD1", "NFKB1", "STAT2", "IRF8")

# ======================== 加载数据 ========================

cat("Loading preprocessed Treg data...\n")
treg_file <- list.files("D:/Research/potato/data/roosted", pattern = "treg_207422\\.rds$", full.names = TRUE)[1]
treg_all <- readRDS(treg_file)
cat(sprintf("  Total: %d genes x %d cells\n", nrow(treg_all), ncol(treg_all)))

# Load metadata to split weak/strong
meta_weak <- fread("D:/Research/tomato/data/gse207422_cellchat_input_corrected/weak_metadata.csv")
meta_strong <- fread("D:/Research/tomato/data/gse207422_cellchat_input_corrected/strong_metadata.csv")

weak_barcodes <- meta_weak[[1]]
strong_barcodes <- meta_strong[[1]]

cat(sprintf("  Weak metadata: %d cells\n", length(weak_barcodes)))
cat(sprintf("  Strong metadata: %d cells\n", length(strong_barcodes)))

# Split
treg_weak <- treg_all[, colnames(treg_all) %in% weak_barcodes]
treg_strong <- treg_all[, colnames(treg_all) %in% strong_barcodes]

cat(sprintf("  Matched Weak: %d cells\n", ncol(treg_weak)))
cat(sprintf("  Matched Strong: %d cells\n", ncol(treg_strong)))

# ======================== 预处理 (force-keep TFs) ========================

prep_data_keep_tfs <- function(mat, label, target_tfs) {
  cat(sprintf("\n[%s] Preprocessing...\n", label))
  
  # Standard filter: min_cells=500
  min_cells <- 500
  filtered_genes <- apply(mat, 1, function(x) sum(x > 0) >= min_cells)
  mat_filt <- mat[filtered_genes, ]
  cat(sprintf("  Retained %d genes (min_cells=%d)\n", nrow(mat_filt), min_cells))
  
  # Force-add target TFs that were filtered out
  missing_tfs <- target_tfs[!(target_tfs %in% rownames(mat_filt))]
  if (length(missing_tfs) > 0) {
    cat(sprintf("  Force-keeping %d target TFs: %s\n", length(missing_tfs), paste(missing_tfs, collapse = ", ")))
    mat_missing <- mat[missing_tfs, , drop = FALSE]
    mat_filt <- rbind(mat_filt, mat_missing)
  }
  
  # log1p normalization
  cat("  Applying log1p...\n")
  mat_filt <- log1p(mat_filt)
  
  cat(sprintf("  Final: %d genes x %d cells\n", nrow(mat_filt), ncol(mat_filt)))
  return(mat_filt)
}

treg_weak <- prep_data_keep_tfs(treg_weak, "Weak", TARGET_TFS)
treg_strong <- prep_data_keep_tfs(treg_strong, "Strong", TARGET_TFS)

# ======================== KO 函数 ========================

run_tf_ko <- function(count_matrix, ko_gene, group_name) {
  cat(sprintf("\n========== KO: %s (%s) ==========\n", ko_gene, group_name))
  
  if (!(ko_gene %in% rownames(count_matrix))) {
    cat(sprintf("  WARNING: %s not found, skipping\n", ko_gene))
    return(NULL)
  }
  
  expr_level <- mean(count_matrix[ko_gene, ])
  cat(sprintf("  %s mean expression: %.4f\n", ko_gene, expr_level))
  
  cat("  Running scTenifoldKnk (nNet=5, nCells=300, K=2)...\n")
  start_time <- Sys.time()
  
  result <- tryCatch({
    scTenifoldKnk(
      countMatrix = count_matrix,
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
      ma_nDim = 2
    )
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("  Finished in %.1f min\n", as.numeric(elapsed)))
  
  if (!is.null(result)) {
    diff_reg <- result$diffRegulation
    cat(sprintf("  %d perturbed genes\n", nrow(diff_reg)))
    
    diff_reg$KO_TF <- ko_gene
    diff_reg$Group <- group_name
    
    # Count significant
    n_sig <- sum(diff_reg$p.adj < 0.05, na.rm = TRUE)
    cat(sprintf("  Significant (p.adj<0.05): %d\n", n_sig))
    
    out_file <- file.path(OUTPUT_DIR, sprintf("ko_%s_%s.csv", group_name, ko_gene))
    write.csv(diff_reg, out_file, row.names = FALSE)
    cat(sprintf("  Saved: %s\n", out_file))
    return(diff_reg)
  }
  
  return(NULL)
}

# ======================== 主流程 ========================

cat("\n", rep("=", 60), "\n", sep = "")
cat("Starting TF KOs (5 TFs x 2 groups)\n")
cat("Expected total: ~50 minutes\n")
cat(rep("=", 60), "\n", sep = "")

all_results <- list()

for (group_name in c("weak", "strong")) {
  count_matrix <- if (group_name == "weak") treg_weak else treg_strong
  
  if (ncol(count_matrix) < 100) {
    cat(sprintf("\nSKIPPING %s: only %d cells\n", group_name, ncol(count_matrix)))
    next
  }
  
  for (tf in TARGET_TFS) {
    res <- run_tf_ko(count_matrix, tf, group_name)
    if (!is.null(res)) {
      all_results[[paste(group_name, tf, sep = "_")]] <- res
    }
  }
}

# ======================== 汇总 ========================

cat("\n", rep("=", 60), "\n", sep = "")
cat("SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")

summary_list <- list()
for (name in names(all_results)) {
  res <- all_results[[name]]
  n_sig <- sum(res$p.adj < 0.05, na.rm = TRUE)
  top_up <- res$gene[which.max(res$Z * (res$p.adj < 0.05))]
  top_down <- res$gene[which.min(res$Z * (res$p.adj < 0.05))]
  summary_list[[name]] <- data.frame(
    Group = unique(res$Group),
    KO_TF = unique(res$KO_TF),
    n_perturbed = n_sig,
    top_up = if (length(top_up) > 0) top_up else NA,
    top_down = if (length(top_down) > 0) top_down else NA,
    stringsAsFactors = FALSE
  )
}

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  print(summary_df, row.names = FALSE)
  write.csv(summary_df, file.path(OUTPUT_DIR, "ko_summary_tf.csv"), row.names = FALSE)
  cat(sprintf("Summary saved to: %s\n", file.path(OUTPUT_DIR, "ko_summary_tf.csv")))
} else {
  cat("No successful KO results.\n")
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("DONE! Results in:", OUTPUT_DIR, "\n")
cat(rep("=", 60), "\n", sep = "")
