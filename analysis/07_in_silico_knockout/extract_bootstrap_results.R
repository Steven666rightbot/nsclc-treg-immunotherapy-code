#!/usr/bin/env Rscript
# Extract Fibroblast->Treg communication from 6 bootstrap CellChat objects

library(CellChat)

BOOT_DIR <- "D:/Research/tomato/results/cellchat_ko_corrected/bootstrap"
OUT_DIR  <- "D:/Research/tomato/results/cellchat_ko_corrected/bootstrap"
STRONG_RDS <- "D:/Research/tomato/results/cellchat_ko_corrected/cellchat_strong.rds"

# Read Strong baseline
strong_cc <- readRDS(STRONG_RDS)
ft_strong <- subsetCommunication(strong_cc)
ft_strong <- ft_strong[ft_strong$source == "Fibroblast" & ft_strong$target == "Treg", ]
strong_prob <- if (nrow(ft_strong) > 0) sum(ft_strong$prob) else 0
strong_n <- nrow(ft_strong)
cat(sprintf("Strong baseline: %d interactions, prob = %.6e\n", strong_n, strong_prob))

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
  df <- subsetCommunication(cc)
  ft <- df[df$source == "Fibroblast" & df$target == "Treg", ]
  prob <- if (nrow(ft) > 0) sum(ft$prob) else 0
  n <- nrow(ft)
  
  results <- rbind(results, data.frame(
    iter = i,
    n_fib_treg = n,
    fib_treg_prob = prob,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("  Bootstrap %d: %d interactions, prob = %.6e\n", i, n, prob))
  rm(cc, df, ft)
  gc()
}

# Summary statistics
med_prob <- median(results$fib_treg_prob)
q025 <- quantile(results$fib_treg_prob, 0.025)
q975 <- quantile(results$fib_treg_prob, 0.975)
mean_prob <- mean(results$fib_treg_prob)
sd_prob <- sd(results$fib_treg_prob)

ratio <- med_prob / strong_prob
ratio_lower <- q025 / strong_prob
ratio_upper <- q975 / strong_prob

cat("\n============================================================\n")
cat("BOOTSTRAP SUMMARY (n=6)\n")
cat("============================================================\n")
print(results, row.names = FALSE)

cat(sprintf("\nFibroblast -> Treg Communication Probability:\n"))
cat(sprintf("  Strong (baseline):     %.6e  (%d interactions)\n", strong_prob, strong_n))
cat(sprintf("  Weak (bootstrap med):  %.6e  (mean = %.6e, sd = %.6e)\n", med_prob, mean_prob, sd_prob))
cat(sprintf("  95%% CI:                [%.6e, %.6e]\n", q025, q975))
cat(sprintf("\nWeak / Strong ratio: %.2fx  [95%% CI: %.2fx, %.2fx]\n", 
            ratio, ratio_lower, ratio_upper))

# Save
write.csv(results, file.path(OUT_DIR, "bootstrap_results_extracted.csv"), row.names = FALSE)
cat(sprintf("\nSaved to: %s\n", file.path(OUT_DIR, "bootstrap_results_extracted.csv")))
