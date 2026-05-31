#!/usr/bin/env Rscript
# Extract Fibroblast->Treg without CellChat package

BOOT_DIR <- "D:/Research/tomato/results/cellchat_ko_corrected/bootstrap"
OUT_DIR  <- "D:/Research/tomato/results/cellchat_ko_corrected/bootstrap"
STRONG_RDS <- "D:/Research/tomato/results/cellchat_ko_corrected/cellchat_strong.rds"

# Helper: try to extract communication data from CellChat object
extract_fib_treg <- function(cc) {
  # Try multiple possible slots where CellChat stores communication data
  prob <- 0
  n <- 0
  
  # Method 1: net$prob array
  if (!is.null(cc@net) && !is.null(cc@net$prob)) {
    prob_array <- cc@net$prob
    # prob_array is typically a 3D array: [source, target, pathway]
    # Or a list of matrices
    if (is.array(prob_array) && length(dim(prob_array)) >= 3) {
      # Get dimension names
      dn <- dimnames(prob_array)
      if (!is.null(dn)) {
        src_idx <- which(dn[[1]] == "Fibroblast")
        tgt_idx <- which(dn[[2]] == "Treg")
        if (length(src_idx) > 0 && length(tgt_idx) > 0) {
          vals <- prob_array[src_idx, tgt_idx, ]
          vals <- vals[!is.na(vals) & vals > 0]
          prob <- sum(vals)
          n <- length(vals)
        }
      }
    } else if (is.list(prob_array)) {
      # Check if it's a list of source-target matrices
      for (item in prob_array) {
        if (is.matrix(item) && !is.null(dimnames(item))) {
          if ("Fibroblast" %in% rownames(item) && "Treg" %in% colnames(item)) {
            val <- item["Fibroblast", "Treg"]
            if (!is.na(val) && val > 0) {
              prob <- prob + val
              n <- n + 1
            }
          }
        }
      }
    }
  }
  
  # Method 2: net$count or net$weight
  if (prob == 0 && !is.null(cc@net) && !is.null(cc@net$count)) {
    count_mat <- cc@net$count
    if (is.matrix(count_mat) && !is.null(dimnames(count_mat))) {
      if ("Fibroblast" %in% rownames(count_mat) && "Treg" %in% colnames(count_mat)) {
        n <- count_mat["Fibroblast", "Treg"]
      }
    }
  }
  
  return(list(prob = prob, n = n))
}

# Read Strong baseline
cat("Reading Strong baseline...\n")
cc_strong <- readRDS(STRONG_RDS)
res_strong <- extract_fib_treg(cc_strong)
cat(sprintf("Strong: %d interactions, prob = %.6e\n", res_strong$n, res_strong$prob))

# Read 6 bootstrap results
results <- data.frame(
  iter = integer(),
  n_fib_treg = integer(),
  fib_treg_prob = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:6) {
  rds_file <- file.path(BOOT_DIR, sprintf("cellchat_weak_bootstrap_%02d.rds", i))
  cat(sprintf("Reading %s ...\n", basename(rds_file)))
  
  cc <- readRDS(rds_file)
  res <- extract_fib_treg(cc)
  
  results <- rbind(results, data.frame(
    iter = i,
    n_fib_treg = res$n,
    fib_treg_prob = res$prob,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("  Bootstrap %d: %d interactions, prob = %.6e\n", i, res$n, res$prob))
  rm(cc)
  gc()
}

# Summary
med_prob <- median(results$fib_treg_prob)
q025 <- quantile(results$fib_treg_prob, 0.025)
q975 <- quantile(results$fib_treg_prob, 0.975)
mean_prob <- mean(results$fib_treg_prob)
sd_prob <- sd(results$fib_treg_prob)

ratio <- med_prob / res_strong$prob
ratio_lower <- q025 / res_strong$prob
ratio_upper <- q975 / res_strong$prob

cat("\n============================================================\n")
cat("BOOTSTRAP SUMMARY (n=6)\n")
cat("============================================================\n")
print(results, row.names = FALSE)

cat(sprintf("\nFibroblast -> Treg Communication Probability:\n"))
cat(sprintf("  Strong (baseline):     %.6e  (%d interactions)\n", res_strong$prob, res_strong$n))
cat(sprintf("  Weak (bootstrap med):  %.6e  (mean = %.6e, sd = %.6e)\n", med_prob, mean_prob, sd_prob))
cat(sprintf("  95%% CI:                [%.6e, %.6e]\n", q025, q975))
cat(sprintf("\nWeak / Strong ratio: %.2fx  [95%% CI: %.2fx, %.2fx]\n", 
            ratio, ratio_lower, ratio_upper))

# Save
write.csv(results, file.path(OUT_DIR, "bootstrap_results_extracted.csv"), row.names = FALSE)
cat(sprintf("\nSaved to: %s\n", file.path(OUT_DIR, "bootstrap_results_extracted.csv")))
