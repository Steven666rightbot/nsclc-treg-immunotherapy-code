# Extract detailed CellChat results for mechanism validation
# Focus on: CD137, NECTIN, VEGF, CCL pathways

library(CellChat)
library(dplyr)

out_dir <- "D:/Research/tomato/go_exploration/cellchat_output"
setwd(out_dir)

# Load merged object
print("Loading merged CellChat object...")
cellchat <- readRDS("cellchat_merged.rds")

# Extract detailed L-R pairs for key pathways
key_pathways <- c("CD137", "NECTIN", "VEGF", "CCL", "THBS", "COLLAGEN", "CXCL", "MHC-I")

results <- list()

for (pathway in key_pathways) {
  cat("\n==========", pathway, "==========\n")
  
  # Get L-R pairs for this pathway
  lr_pairs <- tryCatch({
    # Access netP slot for pathway info
    netP <- cellchat@netP
    
    # Get the L-R pair information
    pairLR <- searchPair(signaling = pathway, pairLR.use = cellchat@LR$LRsig, key = "pathway_name", matching.exact = TRUE)
    
    if (nrow(pairLR) > 0) {
      cat("L-R pairs:", nrow(pairLR), "\n")
      print(pairLR %>% select(interaction_name, pathway_name, ligand, receptor))
      
      # Get probability for each group
      prob <- subsetCommunication(cellchat, signaling = pathway)
      if (!is.null(prob) && nrow(prob) > 0) {
        cat("\nCommunication probabilities:\n")
        print(prob %>% select(source, target, prob, pval) %>% head(20))
        
        # Save to results
        results[[pathway]] <- list(
          lr_pairs = pairLR,
          communication = prob
        )
      }
    } else {
      cat("No L-R pairs found\n")
    }
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    NULL
  })
}

# Extract all differential interactions
cat("\n\n========== ALL DIFFERENTIAL INTERACTIONS ==========\n")

# Access net slot for all interactions
net <- cellchat@net

# Get differential network
diffnet <- netVisual_diffInteraction(cellchat, weight.scale = TRUE, return.data = TRUE)

# Save comprehensive results
saveRDS(results, "detailed_pathway_results.rds")

# Export summary table
summary_df <- data.frame()
for (pw in names(results)) {
  if (!is.null(results[[pw]]$communication)) {
    comm <- results[[pw]]$communication
    if (nrow(comm) > 0) {
      comm$pathway <- pw
      summary_df <- rbind(summary_df, comm)
    }
  }
}

if (nrow(summary_df) > 0) {
  write.csv(summary_df, "detailed_communication_table.csv", row.names = FALSE)
  cat("\nSaved detailed_communication_table.csv with", nrow(summary_df), "interactions\n")
}

# Specific focus: TNK as source or target (Treg-related)
cat("\n\n========== TNK-INVOLVED COMMUNICATIONS ==========\n")
if (nrow(summary_df) > 0) {
  tnk_comm <- summary_df %>% 
    filter(grepl("TNK", source) | grepl("TNK", target)) %>%
    arrange(desc(prob))
  
  cat("TNK-involved interactions:", nrow(tnk_comm), "\n")
  print(tnk_comm %>% select(pathway, source, target, prob, pval) %>% head(30))
  
  write.csv(tnk_comm, "tnk_communication.csv", row.names = FALSE)
}

print("Done!")
