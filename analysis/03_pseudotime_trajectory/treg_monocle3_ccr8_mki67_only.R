#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# GSE243013 - CCR8+ / MKI67+ Treg 专属 Monocle3 分析
# 
# 核心设计：
#   1. 完全排除 FOXP3+ Treg（XGBoost+SHAP 识别的关键亚型只有 CCR8+/MKI67+）
#   2. 在 CCR8+MKI67+ 混合群体上重建轨迹
#   3. 按 Response 分组各自独立跑 Monocle3（各自学图、各自定根、各自算伪时序）
#   4. 机械记忆基因 + Treg关键基因 随各自伪时序变动
#   5. 分组 Moran's I 计算 + 对比
#   6. 假设：Non-responder（弱组）呈现更明显的趋势

library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix)
library(gridExtra)
library(reshape2)

# ===================== 路径配置 =====================
CDS_PATH <- "D:/Research/tomato/data/monocle3_treg/treg_monocle3_cds.rds"
OUT_DIR <- "D:/Research/tomato/figures/monocle3_ccr8_mki67"
RESULT_DIR <- "D:/Research/tomato/results/monocle3_ccr8_mki67"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("Monocle3: CCR8+ / MKI67+ Treg Only (FOXP3+ Excluded)\n")
cat("============================================================\n")

# ===================== 1. 加载CDS并过滤FOXP3+ =====================
cat("\n[Step 1] Loading CDS and filtering out FOXP3+ cells...\n")
cds_full <- readRDS(CDS_PATH)
cat("Original CDS:", ncol(cds_full), "cells x", nrow(cds_full), "genes\n")

# 只保留 CCR8+ 和 MKI67+
keep_cells <- colnames(cds_full)[colData(cds_full)$sub_cell_type %in% c("CD4T_Treg_CCR8", "CD4T_Treg_MKI67")]
cds_cm <- cds_full[, keep_cells]
cat("After excluding FOXP3+:", ncol(cds_cm), "cells\n")
cat("  CCR8+:", sum(colData(cds_cm)$sub_cell_type == "CD4T_Treg_CCR8"), "\n")
cat("  MKI67+:", sum(colData(cds_cm)$sub_cell_type == "CD4T_Treg_MKI67"), "\n")

# 过滤零read细胞
cds_cm <- cds_cm[, Matrix::colSums(SingleCellExperiment::counts(cds_cm)) != 0]
cat("After zero-filter:", ncol(cds_cm), "cells\n")

# ===================== 2. 在CCR8+/MKI67+上重建Monocle3 =====================
cat("\n[Step 2] Rebuilding Monocle3 on CCR8+/MKI67+ cells...\n")

cds_cm <- preprocess_cds(cds_cm, num_dim = 30, method = "PCA")
cat("Preprocessing done\n")

cds_cm <- reduce_dimension(cds_cm, reduction_method = "UMAP", preprocess_method = "PCA",
                           umap.n_neighbors = 30, umap.min_dist = 0.1)
cat("UMAP done\n")

cds_cm <- cluster_cells(cds_cm, reduction_method = "UMAP")
cat("Clustering done:", length(unique(clusters(cds_cm, reduction_method = "UMAP"))), "clusters\n")

cds_cm <- learn_graph(cds_cm, use_partition = FALSE)
cat("Graph learned\n")

# 根节点选择：以CCR8+最密集区域为根（假设CCR8+是更基础的状态，MKI67+是增殖衍生的）
# 也可以尝试以MKI67+为根，但先以CCR8+为根
ccr8_cells <- which(colData(cds_cm)$sub_cell_type == "CD4T_Treg_CCR8")
umap_cm <- reducedDims(cds_cm)$UMAP
ccr8_center <- colMeans(umap_cm[ccr8_cells, ])
dists <- apply(umap_cm, 1, function(x) sum((x - ccr8_center)^2))
root_cell_cm <- names(which.min(dists))
cat("Root cell (CCR8+ center):", root_cell_cm, "\n")
cat("Root subtype:", colData(cds_cm)[root_cell_cm, "sub_cell_type"], "\n")

cds_cm <- order_cells(cds_cm, root_cells = root_cell_cm)
cat("Pseudotime ordered\n")

# 保存整体CDS
saveRDS(cds_cm, file.path(RESULT_DIR, "cds_ccr8_mki67_full.rds"))

# 可视化整体
cat("\n[Step 3] Visualizing full CCR8+/MKI67+ trajectory...\n")

p1 <- plot_cells(cds_cm, color_cells_by = "sub_cell_type",
                 label_cell_groups = FALSE, cell_size = 0.5,
                 label_leaves = FALSE, label_branch_points = FALSE) +
  scale_color_manual(values = c(
    "CD4T_Treg_CCR8" = "#ff7f0e",
    "CD4T_Treg_MKI67" = "#2ca02c"
  )) +
  ggtitle("CCR8+ / MKI67+ Trajectory") +
  theme(legend.position = "right")
ggsave(file.path(OUT_DIR, "01_trajectory_subtype.png"), p1, width = 8, height = 6, dpi = 300)

p2 <- plot_cells(cds_cm, color_cells_by = "pseudotime",
                 label_cell_groups = FALSE, cell_size = 0.5,
                 label_leaves = FALSE, label_branch_points = FALSE) +
  scale_color_viridis_c() +
  ggtitle("Pseudotime (CCR8+ as Root)") +
  theme(legend.position = "right")
ggsave(file.path(OUT_DIR, "02_trajectory_pseudotime.png"), p2, width = 8, height = 6, dpi = 300)

p3 <- plot_cells(cds_cm, color_cells_by = "response_binary",
                 label_cell_groups = FALSE, cell_size = 0.5,
                 label_leaves = FALSE, label_branch_points = FALSE) +
  scale_color_manual(values = c(
    "Responder" = "#d62728",
    "Non-responder" = "#9467bd"
  )) +
  ggtitle("Response Status") +
  theme(legend.position = "right")
ggsave(file.path(OUT_DIR, "03_trajectory_response.png"), p3, width = 8, height = 6, dpi = 300)

# ===================== 4. 定义关键基因集 =====================
cat("\n[Step 4] Defining key gene sets...\n")

mechanical_memory_genes <- c(
  "ACTA2", "MYL9", "VCL", "TLN1", "TRPV4", "MOB1B",
  "ITGA1", "DDR1", "DDR2", "SRC", "ITGB1", "CD44", "PTK2", "ZYX", "FLNA",
  "ACTG1", "RHOA", "ROCK1", "ROCK2", "MYH9", "MYH10", "PFN1", "CFL1"
)

treg_key_genes <- c(
  "IL2RA", "CTLA4", "TIGIT", "LAG3", "PDCD1",
  "MKI67", "CCR8", "TGFB1", "IL10", "GZMB", "PRF1", "IL2", "FOXP3"
)

all_key_genes <- unique(c(mechanical_memory_genes, treg_key_genes))
available_key_genes <- all_key_genes[all_key_genes %in% rowData(cds_cm)$gene_short_name]
mech_available <- intersect(mechanical_memory_genes, available_key_genes)
treg_available <- intersect(treg_key_genes, available_key_genes)

cat("Mechanical memory genes:", length(mech_available), "/", length(mechanical_memory_genes), "\n")
cat("Treg key genes:", length(treg_available), "/", length(treg_key_genes), "\n")

# ===================== 5. 分组独立Monocle3分析函数 =====================
cat("\n[Step 5] Group-specific Monocle3 analysis...\n")

run_group_analysis <- function(cds_source, group_name, response_label, 
                               out_dir, result_dir, available_genes) {
  cat("\n==========", group_name, "==========\n")
  
  # 提取该组细胞
  cells <- colnames(cds_source)[colData(cds_source)$response_binary == response_label]
  cat("Cells:", length(cells), "\n")
  
  cds_g <- cds_source[, cells]
  cds_g <- cds_g[, Matrix::colSums(SingleCellExperiment::counts(cds_g)) != 0]
  cat("After filter:", ncol(cds_g), "cells\n")
  
  # 独立重建Monocle3
  cat("Preprocessing...\n")
  cds_g <- preprocess_cds(cds_g, num_dim = 30, method = "PCA")
  
  cat("UMAP...\n")
  cds_g <- reduce_dimension(cds_g, reduction_method = "UMAP", preprocess_method = "PCA",
                            umap.n_neighbors = 30, umap.min_dist = 0.1)
  
  cat("Clustering...\n")
  cds_g <- cluster_cells(cds_g, reduction_method = "UMAP")
  
  cat("Learning graph...\n")
  cds_g <- learn_graph(cds_g, use_partition = FALSE)
  
  # 根节点选择：以该组内CCR8+最密集区域为根
  ccr8_idx <- which(colData(cds_g)$sub_cell_type == "CD4T_Treg_CCR8")
  if (length(ccr8_idx) > 0) {
    umap_g <- reducedDims(cds_g)$UMAP
    ccr8_center <- colMeans(umap_g[ccr8_idx, , drop = FALSE])
    dists <- apply(umap_g, 1, function(x) sum((x - ccr8_center)^2))
    root_g <- names(which.min(dists))
  } else {
    # 如果没有CCR8+，以UMAP中心为根
    umap_g <- reducedDims(cds_g)$UMAP
    center <- colMeans(umap_g)
    dists <- apply(umap_g, 1, function(x) sum((x - center)^2))
    root_g <- names(which.min(dists))
  }
  cat("Root cell:", root_g, " Subtype:", colData(cds_g)[root_g, "sub_cell_type"], "\n")
  
  cds_g <- order_cells(cds_g, root_cells = root_g)
  cat("Pseudotime ordered. Range:", min(pseudotime(cds_g)), "-", max(pseudotime(cds_g)), "\n")
  
  # 保存
  saveRDS(cds_g, file.path(result_dir, paste0("cds_", group_name, ".rds")))
  
  # 伪时序统计
  pt_df <- data.frame(
    cell_id = colnames(cds_g),
    pseudotime = pseudotime(cds_g),
    sub_cell_type = colData(cds_g)$sub_cell_type,
    stringsAsFactors = FALSE
  )
  write.csv(pt_df, file.path(result_dir, paste0("pseudotime_", group_name, ".csv")), row.names = FALSE)
  
  cat("Pseudotime by subtype:\n")
  print(pt_df %>% group_by(sub_cell_type) %>% summarise(median_pt = median(pseudotime), n = n(), .groups = "drop"))
  
  # Moran's I
  cat("Computing Moran's I...\n")
  cds_g_key <- cds_g[available_genes, ]
  gene_fits <- graph_test(cds_g_key, neighbor_graph = "principal_graph", cores = 1)
  gene_fits <- gene_fits[order(gene_fits$morans_I, decreasing = TRUE), ]
  gene_fits$group <- group_name
  write.csv(gene_fits, file.path(result_dir, paste0("moransI_", group_name, ".csv")))
  
  cat("Top Moran's I:\n")
  print(head(gene_fits[, c("gene_short_name", "morans_I", "q_value")], 8))
  
  # Pseudotime-expression correlation
  cat("Computing pseudotime-expression correlation...\n")
  pt_vec <- pseudotime(cds_g)
  corr_df <- data.frame(
    gene = available_genes,
    spearman_r = NA_real_,
    spearman_p = NA_real_,
    stringsAsFactors = FALSE
  )
  
  expr_mat <- SingleCellExperiment::counts(cds_g)
  for (i in seq_along(available_genes)) {
    g <- available_genes[i]
    expr <- as.numeric(expr_mat[g, ])
    if (sum(expr > 0) > 10) {
      ct <- cor.test(pt_vec, expr, method = "spearman", exact = FALSE)
      corr_df$spearman_r[i] <- ct$estimate
      corr_df$spearman_p[i] <- ct$p.value
    }
  }
  write.csv(corr_df, file.path(result_dir, paste0("spearman_", group_name, ".csv")), row.names = FALSE)
  
  # 可视化
  cat("Generating plots...\n")
  
  # UMAP
  p_umap <- plot_cells(cds_g, color_cells_by = "pseudotime",
                       label_cell_groups = FALSE, cell_size = 0.5,
                       label_leaves = FALSE, label_branch_points = FALSE) +
    scale_color_viridis_c() +
    ggtitle(paste(group_name, "- Pseudotime")) +
    theme(legend.position = "right")
  ggsave(file.path(out_dir, paste0("04_umap_pseudotime_", group_name, ".png")), 
         p_umap, width = 7, height = 6, dpi = 300)
  
  # Subtype
  p_sub <- plot_cells(cds_g, color_cells_by = "sub_cell_type",
                      label_cell_groups = FALSE, cell_size = 0.5,
                      label_leaves = FALSE, label_branch_points = FALSE) +
    scale_color_manual(values = c(
      "CD4T_Treg_CCR8" = "#ff7f0e",
      "CD4T_Treg_MKI67" = "#2ca02c"
    )) +
    ggtitle(paste(group_name, "- Subtype")) +
    theme(legend.position = "right")
  ggsave(file.path(out_dir, paste0("05_umap_subtype_", group_name, ".png")), 
         p_sub, width = 7, height = 6, dpi = 300)
  
  # Density
  pt_df$sub_cell_type <- factor(pt_df$sub_cell_type, levels = c("CD4T_Treg_CCR8", "CD4T_Treg_MKI67"))
  p_den <- ggplot(pt_df, aes(x = pseudotime, fill = sub_cell_type)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("CD4T_Treg_CCR8" = "#ff7f0e", "CD4T_Treg_MKI67" = "#2ca02c")) +
    labs(x = "Pseudotime", y = "Density", fill = "Subtype") +
    ggtitle(paste(group_name, "- Pseudotime Distribution")) +
    theme_bw()
  ggsave(file.path(out_dir, paste0("06_density_", group_name, ".png")), 
         p_den, width = 7, height = 5, dpi = 300)
  
  return(list(cds = cds_g, gene_fits = gene_fits, corr = corr_df, pt = pt_df))
}

# 运行两组分析
resp_res <- run_group_analysis(cds_cm, "Responder", "Responder", OUT_DIR, RESULT_DIR, available_key_genes)
nonresp_res <- run_group_analysis(cds_cm, "Non_responder", "Non-responder", OUT_DIR, RESULT_DIR, available_key_genes)

# ===================== 6. 对比分析 =====================
cat("\n[Step 6] Cross-group comparison...\n")

df_resp <- resp_res$gene_fits[, c("gene_short_name", "morans_I", "q_value")]
df_nonresp <- nonresp_res$gene_fits[, c("gene_short_name", "morans_I", "q_value")]
df_resp_corr <- resp_res$corr
df_nonresp_corr <- nonresp_res$corr

colnames(df_resp) <- c("gene", "moransI_resp", "q_resp")
colnames(df_nonresp) <- c("gene", "moransI_nonresp", "q_nonresp")

comparison <- merge(df_resp, df_nonresp, by = "gene", all = TRUE)
comparison <- merge(comparison, df_resp_corr, by = "gene", all = TRUE)
colnames(comparison)[ncol(comparison)-1] <- "spearman_r_resp"
colnames(comparison)[ncol(comparison)] <- "spearman_p_resp"
comparison <- merge(comparison, df_nonresp_corr, by = "gene", all = TRUE)
colnames(comparison)[ncol(comparison)-1] <- "spearman_r_nonresp"
colnames(comparison)[ncol(comparison)] <- "spearman_p_nonresp"

comparison$moransI_diff <- comparison$moransI_nonresp - comparison$moransI_resp
comparison$spearman_diff <- comparison$spearman_r_nonresp - comparison$spearman_r_resp

# 分类
comparison$category <- ifelse(comparison$gene %in% mechanical_memory_genes, 
                               "Mechanical_Memory", 
                               ifelse(comparison$gene %in% treg_key_genes, 
                                      "Treg_Key", "Other"))

# Fisher Z-test
n_resp <- ncol(resp_res$cds)
n_nonresp <- ncol(nonresp_res$cds)
comparison$fisher_z <- NA
comparison$fisher_p <- NA

for (i in 1:nrow(comparison)) {
  r1 <- comparison$spearman_r_resp[i]
  r2 <- comparison$spearman_r_nonresp[i]
  if (!is.na(r1) && !is.na(r2)) {
    z1 <- atanh(r1)
    z2 <- atanh(r2)
    se <- sqrt(1/(n_resp - 3) + 1/(n_nonresp - 3))
    z_stat <- (z2 - z1) / se
    comparison$fisher_z[i] <- z_stat
    comparison$fisher_p[i] <- 2 * (1 - pnorm(abs(z_stat)))
  }
}
comparison$fisher_fdr <- p.adjust(comparison$fisher_p, method = "fdr")

# 排序
comparison <- comparison[order(comparison$moransI_diff, decreasing = TRUE), ]
write.csv(comparison, file.path(RESULT_DIR, "comparison_moransI_spearman.csv"), row.names = FALSE)

cat("\n--- Top genes: Higher Moran's I in Non-responder ---\n")
print(head(comparison[, c("gene", "category", "moransI_resp", "moransI_nonresp", "moransI_diff", 
                          "spearman_r_resp", "spearman_r_nonresp", "fisher_p")], 15))

cat("\n--- Top genes: Higher Moran's I in Responder ---\n")
print(tail(comparison[, c("gene", "category", "moransI_resp", "moransI_nonresp", "moransI_diff",
                          "spearman_r_resp", "spearman_r_nonresp", "fisher_p")], 10))

# ===================== 7. 类别汇总 =====================
cat("\n[Step 7] Category summary...\n")

cat_sum <- comparison %>%
  group_by(category) %>%
  summarise(
    n = n(),
    mean_moransI_resp = mean(moransI_resp, na.rm = TRUE),
    mean_moransI_nonresp = mean(moransI_nonresp, na.rm = TRUE),
    mean_moransI_diff = mean(moransI_diff, na.rm = TRUE),
    mean_spearman_resp = mean(spearman_r_resp, na.rm = TRUE),
    mean_spearman_nonresp = mean(spearman_r_nonresp, na.rm = TRUE),
    mean_spearman_diff = mean(spearman_diff, na.rm = TRUE),
    .groups = "drop"
  )
print(cat_sum)
write.csv(cat_sum, file.path(RESULT_DIR, "category_summary.csv"), row.names = FALSE)

# ===================== 8. 可视化对比 =====================
cat("\n[Step 8] Comparison plots...\n")

# Moran's I scatter
p_moran <- ggplot(comparison, aes(x = moransI_resp, y = moransI_nonresp, color = category)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_text(data = comparison[abs(comparison$moransI_diff) > 0.02 & !is.na(comparison$moransI_diff), ],
            aes(label = gene), vjust = -1, size = 3) +
  scale_color_manual(values = c("Mechanical_Memory" = "#E41A1C", 
                                 "Treg_Key" = "#377EB8", 
                                 "Other" = "#999999")) +
  labs(x = "Moran's I (Responder)", y = "Moran's I (Non-responder)", color = "Category") +
  ggtitle("Moran's I Comparison") +
  theme_bw()
ggsave(file.path(OUT_DIR, "07_moransI_scatter.png"), p_moran, width = 8, height = 7, dpi = 300)

# Moran's I diff barplot
p_diff <- ggplot(comparison, aes(x = reorder(gene, moransI_diff), y = moransI_diff, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Mechanical_Memory" = "#E41A1C", 
                                "Treg_Key" = "#377EB8", 
                                "Other" = "#999999")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Gene", y = "Δ Moran's I (Non-responder - Responder)") +
  ggtitle("Moran's I Difference") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))
ggsave(file.path(OUT_DIR, "08_moransI_diff_barplot.png"), p_diff, width = 8, height = 10, dpi = 300)

# Spearman scatter
p_spear <- ggplot(comparison, aes(x = spearman_r_resp, y = spearman_r_nonresp, color = category)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_text(data = comparison[abs(comparison$spearman_diff) > 0.02 & !is.na(comparison$spearman_diff), ],
            aes(label = gene), vjust = -1, size = 3) +
  scale_color_manual(values = c("Mechanical_Memory" = "#E41A1C", 
                                 "Treg_Key" = "#377EB8", 
                                 "Other" = "#999999")) +
  labs(x = "Spearman r (Responder)", y = "Spearman r (Non-responder)", color = "Category") +
  ggtitle("Pseudotime-Expression Correlation") +
  theme_bw()
ggsave(file.path(OUT_DIR, "09_spearman_scatter.png"), p_spear, width = 8, height = 7, dpi = 300)

# Category boxplot
comp_long <- rbind(
  data.frame(gene = comparison$gene, category = comparison$category, 
             value = comparison$moransI_resp, group = "Responder", metric = "MoransI"),
  data.frame(gene = comparison$gene, category = comparison$category, 
             value = comparison$moransI_nonresp, group = "Non-responder", metric = "MoransI"),
  data.frame(gene = comparison$gene, category = comparison$category, 
             value = comparison$spearman_r_resp, group = "Responder", metric = "Spearman"),
  data.frame(gene = comparison$gene, category = comparison$category, 
             value = comparison$spearman_r_nonresp, group = "Non-responder", metric = "Spearman")
)

p_box_moran <- ggplot(comp_long[comp_long$metric == "MoransI", ], 
                      aes(x = category, y = value, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Responder" = "#d62728", "Non-responder" = "#9467bd")) +
  labs(x = "", y = "Moran's I", fill = "Group") +
  ggtitle("Moran's I by Category") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(OUT_DIR, "10_moransI_boxplot_category.png"), p_box_moran, width = 7, height = 5, dpi = 300)

# ===================== 9. 关键基因LOESS趋势图 =====================
cat("\n[Step 9] LOESS trend plots for priority genes...\n")

priority_genes <- c("ACTA2", "MYL9", "VCL", "TLN1", "TRPV4", "ITGA1", "DDR1", "DDR2", "SRC",
                    "MKI67", "CCR8", "CTLA4", "IL2RA", "FOXP3")
priority_genes <- priority_genes[priority_genes %in% available_key_genes]

# 合并两组伪时序数据
pt_resp <- resp_res$pt
pt_resp$group <- "Responder"
pt_nonresp <- nonresp_res$pt
pt_nonresp$group <- "Non-responder"

# 获取各自的表达
expr_resp <- as.matrix(SingleCellExperiment::counts(resp_res$cds)[priority_genes, ])
expr_nonresp <- as.matrix(SingleCellExperiment::counts(nonresp_res$cds)[priority_genes, ])

for (g in priority_genes) {
  cat("Plotting", g, "...\n")
  
  df_g <- rbind(
    data.frame(
      pseudotime = pt_resp$pseudotime,
      expression = as.numeric(expr_resp[g, ]),
      sub_cell_type = pt_resp$sub_cell_type,
      group = "Responder",
      stringsAsFactors = FALSE
    ),
    data.frame(
      pseudotime = pt_nonresp$pseudotime,
      expression = as.numeric(expr_nonresp[g, ]),
      sub_cell_type = pt_nonresp$sub_cell_type,
      group = "Non-responder",
      stringsAsFactors = FALSE
    )
  )
  
  # 只保留表达>0的细胞用于绘图（减少视觉噪音）
  df_plot <- df_g[df_g$expression > 0, ]
  
  if (nrow(df_plot) > 50) {
    # 计算各组Spearman r用于标注
    r_resp <- cor(df_g$pseudotime[df_g$group == "Responder"], 
                  df_g$expression[df_g$group == "Responder"], 
                  method = "spearman", use = "pairwise.complete.obs")
    r_non <- cor(df_g$pseudotime[df_g$group == "Non-responder"], 
                 df_g$expression[df_g$group == "Non-responder"], 
                 method = "spearman", use = "pairwise.complete.obs")
    
    p_trend <- ggplot(df_plot, aes(x = pseudotime, y = expression, color = group)) +
      geom_point(alpha = 0.08, size = 0.4) +
      geom_smooth(method = "loess", se = TRUE, aes(fill = group), alpha = 0.15, span = 0.75) +
      scale_color_manual(values = c("Responder" = "#d62728", "Non-responder" = "#9467bd")) +
      scale_fill_manual(values = c("Responder" = "#d62728", "Non-responder" = "#9467bd")) +
      labs(x = "Pseudotime", y = "Expression", color = "Group", fill = "Group") +
      ggtitle(paste0(g, " (r_resp=", round(r_resp, 3), ", r_non=", round(r_non, 3), ")")) +
      theme_bw()
    
    ggsave(file.path(OUT_DIR, sprintf("11_trend_%s.png", g)), p_trend, width = 8, height = 5, dpi = 300)
  }
}

# ===================== 10. 最终总结 =====================
cat("\n============================================================\n")
cat("Analysis Complete!\n")
cat("============================================================\n")
cat("Output:", OUT_DIR, "\n")
cat("Results:", RESULT_DIR, "\n")
cat("\nCells analyzed:\n")
cat("  Full (CCR8+/MKI67+):", ncol(cds_cm), "\n")
cat("  Responder:", ncol(resp_res$cds), "\n")
cat("  Non-responder:", ncol(nonresp_res$cds), "\n")
cat("\nHypothesis test (Non-responder > Responder):\n")
cat("  Mean Moran's I diff:", mean(comparison$moransI_diff, na.rm = TRUE), "\n")
cat("  Mechanical memory mean diff:", 
    mean(comparison$moransI_diff[comparison$category == "Mechanical_Memory"], na.rm = TRUE), "\n")
cat("  Treg key mean diff:", 
    mean(comparison$moransI_diff[comparison$category == "Treg_Key"], na.rm = TRUE), "\n")
cat("  Mean Spearman r diff:", mean(comparison$spearman_diff, na.rm = TRUE), "\n")
cat("  Mechanical memory Spearman diff:", 
    mean(comparison$spearman_diff[comparison$category == "Mechanical_Memory"], na.rm = TRUE), "\n")
