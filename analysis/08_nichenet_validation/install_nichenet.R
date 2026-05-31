# Install NicheNet
library(remotes)
cat("Installing nichenetr...\n")
result <- tryCatch({
  install_github("saeyslab/nichenetr", quiet = FALSE, upgrade = "never")
  cat("Installation succeeded!\n")
}, error = function(e) {
  cat(sprintf("ERROR: %s\n", conditionMessage(e)))
})
cat("Package installed:\n")
cat(paste(installed.packages()["nichenetr", "Version"], collapse = ""))
cat("\n")
