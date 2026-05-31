#!/usr/bin/env Rscript
# WGCNA analysis: find co-expression modules associated with Contractile Treg score

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4)

FIG_DIR <- r"(D:\Research\tomato\figures)"
RES_DIR <- r"(D:\Research\tomato\results\wgcna)"
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

# ===================== 1. Load expression data =====================
cat("[1/6] Loading TCGA-LUAD expression matrix...\n")
expr <- read_tsv(r"(D:\Research\tomato\data\tcga_luad_expression.tsv.gz)", show_col_types = FALSE)
# First column is gene ID
gene_ids <- expr[[1]]
expr_mat <- as.matrix(expr[, -1])
rownames(expr_mat) <- gene_ids
dim(expr_mat)

# Transpose: WGCNA expects genes as rows, samples as columns (already is)
cat(sprintf("  Genes: %d, Samples: %d\n", nrow(expr_mat), ncol(expr_mat)))

# ===================== 2. Filter genes =====================
cat("[2/6] Filtering low-variance genes...\n")
# Keep top 5000 genes by MAD (robust variance filter)
mads <- apply(expr_mat, 1, mad, na.rm = TRUE)
keep <- mads > 0.5 & rowMeans(expr_mat > 0, na.rm = TRUE) > 0.1
expr_filt <- expr_mat[keep, ]
# Further reduce to top 5000 by MAD for speed
if (nrow(expr_filt) > 5000) {
  mad_keep <- mads[keep]
  expr_filt <- expr_filt[order(mad_keep, decreasing = TRUE)[1:5000], ]
}
cat(sprintf("  Retained: %d genes\n", nrow(expr_filt)))

# ===================== 3. Load Contractile score (trait) =====================
cat("[3/6] Loading Contractile Treg score...\n")
score_df <- read_csv(r"(D:\Research\tomato\results\tcga_validation\tcga_luad_contractile_scores.csv)", show_col_types = FALSE)
score_df <- score_df %>% rename(sample = 1, contractile_score = 2)

# Match samples
trait <- data.frame(
  contractile_score = score_df$contractile_score[match(colnames(expr_filt), score_df$sample)]
)
rownames(trait) <- colnames(expr_filt)
trait <- trait %>% drop_na()

# Filter to common samples
common_samples <- intersect(colnames(expr_filt), rownames(trait))
expr_filt <- expr_filt[, common_samples]
trait <- trait[common_samples, , drop = FALSE]
cat(sprintf("  Common samples: %d\n", length(common_samples)))

# ===================== 4. Choose soft threshold =====================
cat("[4/6] Choosing soft threshold...\n")
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(t(expr_filt), powerVector = powers, verbose = 0)

# Auto-select beta: first power where scale-free R^2 > 0.85
r2 <- sft$fitIndices$SFT.R.sq
beta <- sft$fitIndices$Power[which(r2 > 0.85)[1]]
if (is.na(beta)) beta <- 6
cat(sprintf("  Selected beta = %d (R^2 = %.3f)\n", beta, r2[which(sft$fitIndices$Power == beta)]))

# Save soft threshold plot
png(file.path(RES_DIR, "wgcna_soft_threshold.png"), width = 1200, height = 600, res = 150)
par(mfrow = c(1, 2))
plot(sft$fitIndices$Power, sft$fitIndices$SFT.R.sq, type = "b", pch = 19, col = "steelblue",
     xlab = "Soft Threshold (beta)", ylab = "Scale Free Topology Model Fit (R^2)",
     main = "Scale Independence")
abline(h = 0.85, col = "red", lty = 2)
plot(sft$fitIndices$Power, sft$fitIndices$mean.k., type = "b", pch = 19, col = "steelblue",
     xlab = "Soft Threshold (beta)", ylab = "Mean Connectivity", main = "Mean Connectivity")
dev.off()

# ===================== 5. Build network & detect modules =====================
cat("[5/6] Building network and detecting modules...\n")
net <- blockwiseModules(t(expr_filt), power = beta, TOMType = "unsigned",
                        minModuleSize = 30, reassignThreshold = 0,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        pamRespectsDendro = FALSE, verbose = 0,
                        maxBlockSize = 5000)

module_colors <- labels2colors(net$colors)
n_modules <- length(unique(module_colors))
cat(sprintf("  Detected %d modules\n", n_modules))

# Module sizes
table(module_colors)

# ===================== 6. Module-trait association =====================
cat("[6/6] Module-trait association...\n")
MEs <- net$MEs
module_trait_cor <- cor(MEs, trait$contractile_score, use = "p")
module_trait_p <- corPvalueStudent(module_trait_cor, nrow(MEs))

# Format for plotting
cor_df <- data.frame(
  Module = rownames(module_trait_cor),
  Correlation = as.vector(module_trait_cor),
  Pvalue = as.vector(module_trait_p)
) %>%
  mutate(Color = gsub("^ME", "", Module)) %>%
  arrange(desc(abs(Correlation)))

write_csv(cor_df, file.path(RES_DIR, "module_trait_correlation.csv"))
cat("\nTop modules associated with Contractile score:\n")
print(as.data.frame(head(cor_df, 10)))

# ===================== 7. Publication-quality heatmap =====================
cat("[Plot] Generating module-trait heatmap...\n")

# Prepare matrix for heatmap (modules × traits)
heatmap_mat <- module_trait_cor
rownames(heatmap_mat) <- gsub("^ME", "", rownames(heatmap_mat))
colnames(heatmap_mat) <- "Contractile Score"

# Order by correlation
ord <- order(heatmap_mat[, 1], decreasing = TRUE)
heatmap_mat <- heatmap_mat[ord, , drop = FALSE]

# Get module colors for row labels
mod_colors <- rownames(heatmap_mat)

png(file.path(FIG_DIR, "pub_fig_wgcna_module_trait.png"), width = 1200, height = 1600, res = 200)
par(mar = c(6, 12, 4, 4))

# Draw heatmap manually for full control
n_rows <- nrow(heatmap_mat)
n_cols <- 1
x_pos <- 1
y_pos <- n_rows:1

# Background
plot(0, 0, type = "n", xlim = c(0.5, 1.5), ylim = c(0.5, n_rows + 0.5),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n",
     main = "WGCNA Module-Trait Association\n(TCGA-LUAD, n = 576)")

# Color scale: blue (negative) -> white (0) -> red (positive)
cor_range <- max(abs(heatmap_mat))
cor_colors <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(101)

for (i in 1:n_rows) {
  val <- heatmap_mat[i, 1]
  idx <- round((val + cor_range) / (2 * cor_range) * 100) + 1
  idx <- max(1, min(101, idx))
  rect(x_pos - 0.4, y_pos[i] - 0.4, x_pos + 0.4, y_pos[i] + 0.4,
       col = cor_colors[idx], border = "white", lwd = 1)
  
  # Correlation value
  text_col <- ifelse(abs(val) > 0.5, "white", "black")
  text(x_pos, y_pos[i], sprintf("%.2f", val), col = text_col, cex = 0.9, font = 2)
  
  # P-value asterisk
  pval <- cor_df$Pvalue[cor_df$Color == rownames(heatmap_mat)[i]]
  stars <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "")))
  text(x_pos + 0.45, y_pos[i], stars, col = "black", cex = 0.8)
}

# Y-axis labels with color dots
for (i in 1:n_rows) {
  mod_name <- rownames(heatmap_mat)[i]
  # Get actual WGCNA color
  actual_color <- mod_name
  if (actual_color == "grey") actual_color <- "#999999"
  
  points(0.15, y_pos[i], pch = 15, col = actual_color, cex = 1.5)
  text(0.25, y_pos[i], mod_name, adj = 0, cex = 0.85, font = 1)
}

# X-axis label
axis(1, at = 1, labels = "Contractile\nTreg Score", las = 1, tick = FALSE, line = 0.5, cex.axis = 0.9)

# Legend
legend("bottomright", legend = c("*** p<0.001", "** p<0.01", "* p<0.05"),
       bty = "n", cex = 0.8, text.col = "#555555")

dev.off()

# ===================== 8. Top module gene export =====================
top_module_label <- as.integer(gsub("^ME", "", cor_df$Module[1]))
top_module_color <- labels2colors(top_module_label)
cat(sprintf("\n[Export] Top module: %s (ME%d, r = %.3f, p = %.2e)\n",
            top_module_color, top_module_label, cor_df$Correlation[1], cor_df$Pvalue[1]))

# Use numeric label to match net$colors (safer than color name matching)
top_genes <- rownames(expr_filt)[net$colors == top_module_label]
writeLines(as.character(top_genes), file.path(RES_DIR, paste0("module_", top_module_color, "_genes.txt")))
cat(sprintf("  Genes in %s module: %d\n", top_module_color, length(top_genes)))

# Export all module gene lists
all_module_df <- data.frame(
  gene = rownames(expr_filt),
  module_label = net$colors,
  module_color = labels2colors(net$colors)
)
write_csv(all_module_df, file.path(RES_DIR, "module_gene_lists.csv"))

# GO enrichment for top module (simple hypergeometric via clusterProfiler if available)
if (length(top_genes) > 0 && requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  ego <- enrichGO(gene = top_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write_csv(as.data.frame(ego), file.path(RES_DIR, paste0("module_", top_module_color, "_GO_BP.csv")))
    
    png(file.path(FIG_DIR, "pub_fig_wgcna_top_module_go.png"), width = 1400, height = 1000, res = 200)
    print(barplot(ego, showCategory = 15, title = paste("GO BP:", top_module_color, "Module")))
    dev.off()
    cat("  GO enrichment plot saved.\n")
  } else {
    cat("  GO enrichment returned no significant results.\n")
  }
} else {
  if (length(top_genes) == 0) {
    cat("  No genes in top module, skipping GO enrichment.\n")
  } else {
    cat("  clusterProfiler not available, skipping GO enrichment.\n")
  }
}

cat("\n[WGCNA Complete] Results in:", RES_DIR, "\n")
cat("[Figures] Saved to:", FIG_DIR, "\n")
