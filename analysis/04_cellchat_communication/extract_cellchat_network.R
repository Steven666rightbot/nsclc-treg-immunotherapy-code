# Extract network details from individual CellChat objects
library(CellChat)
library(dplyr)

out_dir <- "D:/Research/tomato/go_exploration/cellchat_output"
setwd(out_dir)

# Load individual objects
chat_strong <- readRDS("cellchat_Strong.rds")
chat_weak <- readRDS("cellchat_Weak.rds")

# Get cell groups
groups <- levels(chat_strong@idents)
cat("Cell groups:", paste(groups, collapse=", "), "\n")

# Key pathways
key_pathways <- c("CD137", "NECTIN", "VEGF", "CCL", "THBS", "COLLAGEN", "CXCL", "MHC-I", "RLN")

# Function to extract pathway network
extract_pathway_net <- function(chat, pathway, group_name) {
  netP <- chat@netP
  if (!pathway %in% names(netP)) return(NULL)
  
  pw_net <- netP[[pathway]]
  if (!is.matrix(pw_net)) return(NULL)
  
  # Ensure dimensions match groups
  n <- min(length(groups), nrow(pw_net), ncol(pw_net))
  mat <- pw_net[1:n, 1:n, drop=FALSE]
  rownames(mat) <- groups[1:n]
  colnames(mat) <- groups[1:n]
  
  # Convert to data frame
  df <- expand.grid(source = rownames(mat), target = colnames(mat), stringsAsFactors = FALSE)
  df$prob <- as.vector(mat)
  df$pathway <- pathway
  df$group <- group_name
  return(df)
}

# Extract all pathway networks
all_net <- data.frame()

for (pw in key_pathways) {
  s_df <- extract_pathway_net(chat_strong, pw, "Strong")
  w_df <- extract_pathway_net(chat_weak, pw, "Weak")
  
  if (!is.null(s_df)) all_net <- rbind(all_net, s_df)
  if (!is.null(w_df)) all_net <- rbind(all_net, w_df)
}

cat("Total interactions extracted:", nrow(all_net), "\n")

# Compare Strong vs Weak
if (nrow(all_net) > 0) {
  # Pivot to compare
  library(tidyr)
  comp <- all_net %>%
    select(source, target, pathway, group, prob) %>%
    pivot_wider(names_from = group, values_from = prob, values_fill = 0) %>%
    mutate(diff = Weak - Strong, fold = ifelse(Strong > 0, Weak/Strong, ifelse(Weak > 0, Inf, 1)))
  
  # Top UP in Weak
  cat("\n========== TOP INTERACTIONS UP IN WEAK ==========\n")
  up_weak <- comp %>% filter(diff > 0) %>% arrange(desc(diff)) %>% head(30)
  print(up_weak %>% select(pathway, source, target, Strong, Weak, diff, fold))
  
  # Top UP in Strong
  cat("\n========== TOP INTERACTIONS UP IN STRONG ==========\n")
  up_strong <- comp %>% filter(diff < 0) %>% arrange(diff) %>% head(30)
  print(up_strong %>% select(pathway, source, target, Strong, Weak, diff, fold))
  
  # Filter for TNK-involved
  cat("\n========== TNK-INVOLVED COMMUNICATIONS ==========\n")
  tnk_comm <- comp %>% 
    filter(grepl("TNK", source) | grepl("TNK", target)) %>%
    filter(abs(diff) > 0.001) %>%
    arrange(desc(abs(diff)))
  cat("TNK-involved interactions (|diff| > 0.001):", nrow(tnk_comm), "\n")
  print(tnk_comm %>% select(pathway, source, target, Strong, Weak, diff) %>% head(40))
  
  # Save
  write.csv(comp, "detailed_network_comparison.csv", row.names = FALSE)
  write.csv(tnk_comm, "tnk_network_comparison.csv", row.names = FALSE)
  cat("\nSaved: detailed_network_comparison.csv\n")
  cat("Saved: tnk_network_comparison.csv\n")
}

print("Done!")
