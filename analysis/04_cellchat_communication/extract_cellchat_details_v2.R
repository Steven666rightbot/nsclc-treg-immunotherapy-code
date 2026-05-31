# Extract detailed CellChat results v2
# Direct slot access to avoid function errors

library(CellChat)
library(dplyr)

out_dir <- "D:/Research/tomato/go_exploration/cellchat_output"
setwd(out_dir)

# Load merged object
print("Loading merged CellChat object...")
cellchat <- readRDS("cellchat_merged.rds")

# Access LR database
lr_db <- cellchat@DB$interaction
print(paste("Total L-R pairs in DB:", nrow(lr_db)))

# Key pathways of interest
key_pathways <- c("CD137", "NECTIN", "VEGF", "CCL", "THBS", "COLLAGEN", "CXCL", "MHC-I", "RLN", "RESISTIN")

# Extract L-R pairs for each pathway
cat("\n========== L-R PAIRS BY PATHWAY ==========\n")
for (pw in key_pathways) {
  pairs <- lr_db %>% filter(pathway_name == pw) %>% 
    select(interaction_name, pathway_name, ligand, receptor, annotation, evidence)
  if (nrow(pairs) > 0) {
    cat("\n---", pw, "(", nrow(pairs), "pairs) ---\n")
    print(pairs %>% select(interaction_name, ligand, receptor))
  }
}

# Access net slot for interaction probabilities
cat("\n\n========== INTERACTION NETWORK ==========\n")
net <- cellchat@net

# Get group labels
groups <- cellchat@idents %>% levels()
cat("Cell groups:", paste(groups, collapse=", "), "\n")

# Extract interaction matrix for each pathway from netP
netP <- cellchat@netP

# Check structure
if (length(netP) > 0) {
  cat("\nPathways in netP:", length(netP), "\n")
  
  # Get pathway names from netP
  pw_names <- names(netP)
  cat("Pathway names:", paste(head(pw_names, 20), collapse=", "), "...\n")
  
  # For each key pathway, extract interaction probabilities
  for (pw in key_pathways) {
    if (pw %in% pw_names) {
      cat("\n---", pw, "network ---\n")
      pw_net <- netP[[pw]]
      cat("Type:", class(pw_net), "\n")
      if (is.matrix(pw_net)) {
        cat("Dimensions:", nrow(pw_net), "x", ncol(pw_net), "\n")
        # Print matrix with group names
        rownames(pw_net) <- groups[1:min(length(groups), nrow(pw_net))]
        colnames(pw_net) <- groups[1:min(length(groups), ncol(pw_net))]
        print(round(pw_net, 6))
      } else if (is.list(pw_net)) {
        cat("Names:", paste(names(pw_net), collapse=", "), "\n")
      }
    }
  }
}

# Extract differential network
net <- cellchat@net
cat("\n\n========== DIFFERENTIAL NETWORK STRUCTURE ==========\n")
cat("net slots:", paste(names(net), collapse=", "), "\n")

# Access differential network
if (!is.null(net$diff)) {
  diff_net <- net$diff
  cat("diff slots:", paste(names(diff_net), collapse=", "), "\n")
  
  # Get differential interaction matrix
  if (!is.null(diff_net$prob)) {
    cat("\nDifferential probability matrix (Weak - Strong):\n")
    diff_prob <- diff_net$prob
    if (is.matrix(diff_prob)) {
      dim_names <- groups[1:min(length(groups), nrow(diff_prob))]
      rownames(diff_prob) <- dim_names
      colnames(diff_prob) <- dim_names
      print(round(diff_prob, 6))
      
      # Find strongest differential interactions
      cat("\n--- Top differential interactions (Weak - Strong) ---\n")
      # Convert to data frame
      diff_df <- expand.grid(source = rownames(diff_prob), target = colnames(diff_prob), stringsAsFactors = FALSE)
      diff_df$prob_diff <- as.vector(diff_prob)
      diff_df <- diff_df %>% filter(source != target) %>% arrange(desc(abs(prob_diff)))
      print(head(diff_df, 20))
    }
  }
}

# Access individual group networks for comparison
cat("\n\n========== INDIVIDUAL GROUP NETWORKS ==========\n")
for (grp in c("Strong", "Weak")) {
  if (!is.null(net[[grp]])) {
    grp_net <- net[[grp]]
    cat("\n", grp, "network:\n")
    if (!is.null(grp_net$prob)) {
      prob_mat <- grp_net$prob
      if (is.matrix(prob_mat) && length(dim(prob_mat)) == 3) {
        # 3D array: source x target x pathway
        cat("Dimensions:", dim(prob_mat), "(source x target x pathway)\n")
        # Sum across pathways
        total_prob <- apply(prob_mat, c(1,2), sum)
        dim_names <- groups[1:min(length(groups), nrow(total_prob))]
        rownames(total_prob) <- dim_names
        colnames(total_prob) <- dim_names
        print(round(total_prob, 6))
      }
    }
  }
}

print("\nDone!")
