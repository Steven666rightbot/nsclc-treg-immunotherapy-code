# ============================================================================
# NicheNet Analysis: Predict downstream targets of CellChat-identified ligands
# Two analyses:
#   1. MHC-I (HLA-A/B/C) → CCR8+ Treg targets
#   2. MIF → MKI67+ Treg targets
# ============================================================================

library(nichenetr)
library(dplyr)
library(ggplot2)

OUTDIR <- "D:/Research/cucumber/fig5/nichenet_output"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== NicheNet Analysis ===\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================
# Load NicheNet's prior knowledge networks
# ============================================================
cat("Loading NicheNet models...\n")

# Ligand-target model (from the package)
ligand_target_matrix <- readRDS(url("https://zenodo.org/records/7074291/files/ligand_target_matrix_nsga2r_final.rds?download=1"))
lr_network <- readRDS(url("https://zenodo.org/records/7074291/files/lr_network_human_21122021.rds?download=1"))
weighted_networks <- readRDS(url("https://zenodo.org/records/7074291/files/weighted_networks_nsga2r_final.rds?download=1"))

cat(sprintf("  Ligand-target matrix: %d x %d\n", nrow(ligand_target_matrix), ncol(ligand_target_matrix)))
cat(sprintf("  LR network: %d edges\n", nrow(lr_network)))
cat(sprintf("  Weighted networks: %d\n", length(weighted_networks)))

# ============================================================
# Helper function: run NicheNet for one analysis
# ============================================================
run_nichenet <- function(label, expressed_genes_file, de_file, 
                         ligands_of_interest, genes_of_interest_n = 200) {
  cat(rep("=", 60), "\n", sep = "")
  cat(sprintf("[%s] NicheNet Analysis\n", label))
  cat(rep("=", 60), "\n", sep = "")
  
  # Load data
  expressed_genes <- readLines(expressed_genes_file)
  de <- read.csv(de_file)
  
  cat(sprintf("  Expressed genes: %d\n", length(expressed_genes)))
  cat(sprintf("  Total DE genes: %d\n", nrow(de)))
  
  # Define genes of interest: top DE genes
  # Take top genes upregulated in each direction
  up_in_r <- de %>% filter(logFC > 0) %>% slice_min(pval_adj, n = genes_of_interest_n) %>% pull(gene)
  up_in_nr <- de %>% filter(logFC < 0) %>% slice_min(pval_adj, n = genes_of_interest_n) %>% pull(gene)
  
  cat(sprintf("  Genes up in Responder: %d\n", length(up_in_r)))
  cat(sprintf("  Genes up in Non-responder: %d\n", length(up_in_nr)))
  
  # Combine both directions as genes of interest
  genes_of_interest <- unique(c(up_in_r, up_in_nr))
  cat(sprintf("  Total genes of interest: %d\n", length(genes_of_interest)))
  
  # Background: all expressed genes
  background <- expressed_genes
  cat(sprintf("  Background genes: %d\n", length(background)))
  
  # ============================================================
  # 1. Predict ligand activities
  # ============================================================
  cat("\nPredicting ligand activities...\n")
  
  ligand_activities <- predict_ligand_activities(
    geneset = genes_of_interest,
    background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = ligands_of_interest
  )
  
  ligand_activities <- ligand_activities %>% 
    arrange(pearson) %>% 
    mutate(rank = rank(pearson, ties.method = "average"))
  
  cat("  Ligand activity ranking:\n")
  print(ligand_activities[, c("test_ligand", "pearson", "rank")])
  
  write.csv(ligand_activities, file.path(OUTDIR, sprintf("%s_ligand_activities.csv", label)), 
            row.names = FALSE)
  
  # ============================================================
  # 2. Predict target genes for top ligands
  # ============================================================
  cat("\nPredicting target genes...\n")
  
  # Get the best-performing ligand(s)
  best_ligand <- ligand_activities %>% slice_max(pearson, n = 1) %>% pull(test_ligand)
  cat(sprintf("  Best ligand: %s (pearson = %.4f)\n", 
              best_ligand, 
              ligand_activities %>% filter(test_ligand == best_ligand) %>% pull(pearson)))
  
  # Predict target genes for each ligand of interest
  all_targets <- data.frame()
  for (lig in ligands_of_interest) {
    # Get target genes with regulatory potential
    targets <- predict_target_genes(
      ligands = list(lig),
      ligand_target_matrix = ligand_target_matrix,
      background_expressed_genes = background,
      top_n = 50
    )
    
    # Filter to genes that are in our genes of interest
    targets_in_goi <- targets[targets %in% genes_of_interest]
    
    cat(sprintf("  %s: %d predicted targets (top 50), %d in genes of interest\n", 
                lig, length(targets), length(targets_in_goi)))
    
    if (length(targets_in_goi) > 0) {
      cat("    First 10: ", paste(head(targets_in_goi, 10), collapse = ", "), "\n")
    }
    
    # Get regulatory potential scores for these targets
    if (length(targets) > 0) {
      reg_pot <- ligand_target_matrix[targets, lig, drop = FALSE]
      reg_pot_df <- data.frame(
        ligand = lig,
        target = rownames(reg_pot),
        regulatory_potential = as.numeric(reg_pot[, 1]),
        stringsAsFactors = FALSE
      )
      reg_pot_df <- reg_pot_df %>% arrange(desc(regulatory_potential))
      all_targets <- rbind(all_targets, reg_pot_df)
    }
  }
  
  if (nrow(all_targets) > 0) {
    write.csv(all_targets, file.path(OUTDIR, sprintf("%s_target_predictions.csv", label)), 
              row.names = FALSE)
    cat(sprintf("  Saved %d target predictions\n", nrow(all_targets)))
  }
  
  # ============================================================
  # 3. Generate heatmap of ligand → target regulatory potential
  # ============================================================
  if (nrow(all_targets) > 0) {
    # Take top 15 targets per ligand
    top_targets <- all_targets %>% 
      group_by(ligand) %>% 
      slice_max(regulatory_potential, n = 15) %>%
      ungroup()
    
    # Create matrix for heatmap
    ligands_plot <- unique(top_targets$ligand)
    targets_plot <- unique(top_targets$target)
    
    mat <- matrix(0, nrow = length(targets_plot), ncol = length(ligands_plot))
    rownames(mat) <- targets_plot
    colnames(mat) <- ligands_plot
    
    for (i in 1:nrow(top_targets)) {
      mat[top_targets$target[i], top_targets$ligand[i]] <- top_targets$regulatory_potential[i]
    }
    
    # Heatmap
    png(file.path(OUTDIR, sprintf("%s_ligand_target_heatmap.png", label)), 
        width = 8, height = max(6, length(targets_plot) * 0.3), 
        units = "in", res = 300)
    
    # Simple heatmap with pheatmap if available, otherwise use base R
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      library(pheatmap)
      pheatmap(mat, 
               cluster_rows = TRUE, cluster_cols = FALSE,
               color = colorRampPalette(c("white", "#2E5AAC"))(100),
               main = paste0("Ligand → Target Regulatory Potential (", label, ")"),
               fontsize_row = 8, fontsize_col = 10,
               display_numbers = FALSE)
    } else {
      heatmap(mat, Colv = NA, scale = "none",
              col = colorRampPalette(c("white", "#2E5AAC"))(100),
              main = paste0("Ligand → Target (", label, ")"))
    }
    
    dev.off()
    cat(sprintf("  Heatmap saved: %s_ligand_target_heatmap.png\n", label))
  }
  
  return(list(
    ligand_activities = ligand_activities,
    target_predictions = all_targets,
    best_ligand = best_ligand,
    up_in_r = up_in_r,
    up_in_nr = up_in_nr
  ))
}

# ============================================================
# Analysis 1: MHC-I → CCR8+ Treg
# ============================================================
res_ccr8 <- run_nichenet(
  label = "CCR8",
  expressed_genes_file = "D:/Research/cucumber/fig5/nichenet_input/expressed_genes_CCR8.csv",
  de_file = "D:/Research/cucumber/fig5/nichenet_input/de_CCR8_R_vs_NR.csv",
  ligands_of_interest = c("HLA-A", "HLA-B", "HLA-C", "HLA-E"),
  genes_of_interest_n = 200
)

# ============================================================
# Analysis 2: MIF → MKI67+ Treg
# ============================================================
res_mki67 <- run_nichenet(
  label = "MKI67",
  expressed_genes_file = "D:/Research/cucumber/fig5/nichenet_input/expressed_genes_MKI67.csv",
  de_file = "D:/Research/cucumber/fig5/nichenet_input/de_MKI67_R_vs_NR.csv",
  ligands_of_interest = c("MIF", "CD74", "CD44"),
  genes_of_interest_n = 200
)

# ============================================================
# Summary
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("NICHENET SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")

cat("\n=== CCR8+ Treg: MHC-I Ligand Activities ===\n")
print(res_ccr8$ligand_activities[, c("test_ligand", "pearson", "rank")])

cat("\n=== MKI67+ Treg: MIF Ligand Activities ===\n")
print(res_mki67$ligand_activities[, c("test_ligand", "pearson", "rank")])

cat("\nDONE!\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
