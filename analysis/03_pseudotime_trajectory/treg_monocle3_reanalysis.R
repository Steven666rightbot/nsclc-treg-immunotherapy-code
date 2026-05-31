#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# GSE243013 - Treg Monocle3 重分析
# 核心变更：
#   1. 重新选择根节点（FOXP3不参与CCR8/MKI67分化路径）
#   2. 机械记忆基因 + Treg关键基因 随伪时序动态
#   3. Moran's I 分组计算（Responder vs Non-responder）并对比
#   4. 假设：Non-responder（弱组）呈现更明显的趋势

library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix)
library(gridExtra)
library(reshape2)

# ===================== 路径配置 =====================
DATA_DIR <- "D:/Research/tomato/data/monocle3_input"
CDS_PATH <- "D:/Research/tomato/data/monocle3_treg/treg_monocle3_cds.rds"
OUT_DIR <- "D:/Research/tomato/figures/monocle3_reanalysis"
RESULT_DIR <- "D:/Research/tomato/results/monocle3_reanalysis"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("Monocle3 Reanalysis: Root Re-selection + Grouped Moran's I\n")
cat("============================================================\n")

# ===================== 1. 加载已有CDS =====================
cat("\n[Step 1] Loading existing CDS...\n")
cds <- readRDS(CDS_PATH)
cat("CDS loaded:", ncol(cds), "cells x", nrow(cds), "genes\n")
cat("Subtypes:", paste(unique(colData(cds)$sub_cell_type), collapse=", "), "\n")

# ===================== 2. 重新选择根节点 =====================
cat("\n[Step 2] Root node re-selection\n")
cat("Finding trajectory graph leaves...\n")

# 获取principal graph
g <- principal_graph(cds)[["UMAP"]]
# 找到degree=1的节点（端点/叶子）
node_degrees <- igraph::degree(g)
leaves <- names(node_degrees[node_degrees == 1])
cat("Graph leaves (endpoints):", length(leaves), "\n")

# 获取UMAP坐标和细胞信息
umap_coords <- reducedDims(cds)$UMAP
meta <- as.data.frame(colData(cds))

# 分析每个端点附近细胞的亚型组成
leaf_analysis <- data.frame(
  leaf = character(),
  n_cells = integer(),
  pct_foxp3 = numeric(),
  pct_ccr8 = numeric(),
  pct_mki67 = numeric(),
  stringsAsFactors = FALSE
)

for (leaf in leaves) {
  # 获取该节点在UMAP中的坐标
  leaf_coord <- igraph::V(g)[leaf]$coordinates
  if (is.null(leaf_coord)) {
    # 某些节点可能没有coordinates属性，从pr_graph_cell_proj_closest_vertex找
    closest <- principal_graph_aux(cds)[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    leaf_cells <- rownames(closest)[closest[,1] == leaf]
  } else {
    # 找距离该节点最近的100个细胞
    dists <- apply(umap_coords, 1, function(x) sqrt(sum((x - leaf_coord)^2)))
    leaf_cells <- names(sort(dists))[1:min(100, length(dists))]
  }
  
  if (length(leaf_cells) > 0) {
    subtypes <- meta[leaf_cells, "sub_cell_type"]
    n_total <- length(subtypes)
    leaf_analysis <- rbind(leaf_analysis, data.frame(
      leaf = leaf,
      n_cells = n_total,
      pct_foxp3 = sum(subtypes == "CD4T_Treg_FOXP3") / n_total * 100,
      pct_ccr8 = sum(subtypes == "CD4T_Treg_CCR8") / n_total * 100,
      pct_mki67 = sum(subtypes == "CD4T_Treg_MKI67") / n_total * 100
    ))
  }
}

cat("\nLeaf composition analysis:\n")
print(leaf_analysis)
write.csv(leaf_analysis, file.path(RESULT_DIR, "leaf_composition.csv"), row.names = FALSE)

# 新根节点选择策略：
# 由于FOXP3不参与CCR8/MKI67分化路径，我们不以任何特定亚型为根。
# 策略A：选择连接最多FOXP3+、最少CCR8+/MKI67+的端点（代表resting state）
# 策略B：选择距离当前根节点最远的端点（翻转轨迹）
# 策略C：选择CCR8+比例最高的端点作为根（假设CCR8+是分化起点）

# 这里我们采用策略A：找到FOXP3比例最高且MKI67比例最低的端点
# 但如果所有端点都混杂，则使用更精细的方法

# 更精细的根节点选择：
# 找到FOXP3+亚型中，CCR8和MKI67表达最低的细胞作为根
# 因为这些细胞代表最"resting"的状态

cat("\nRefined root selection: finding resting FOXP3+ cells...\n")

# 获取FOXP3+细胞的表达
foxp3_cells <- which(meta$sub_cell_type == "CD4T_Treg_FOXP3")
cat("FOXP3+ cells:", length(foxp3_cells), "\n")

# 获取CCR8和MKI67的表达（如果它们在数据集中）
expr_matrix <- SingleCellExperiment::counts(cds)
available_genes <- rowData(cds)$gene_short_name

has_ccr8 <- "CCR8" %in% available_genes
has_mki67 <- "MKI67" %in% available_genes

if (has_ccr8 && has_mki67) {
  ccr8_expr <- as.numeric(expr_matrix["CCR8", foxp3_cells])
  mki67_expr <- as.numeric(expr_matrix["MKI67", foxp3_cells])
  
  # 找到CCR8和MKI67都低的细胞（双低的10%）
  ccr8_rank <- rank(ccr8_expr) / length(ccr8_expr)
  mki67_rank <- rank(mki67_expr) / length(mki67_expr)
  combined_rank <- ccr8_rank + mki67_rank
  
  resting_idx <- order(combined_rank)[1:max(50, round(length(foxp3_cells) * 0.05))]
  resting_cells <- colnames(cds)[foxp3_cells[resting_idx]]
  
  cat("Selected", length(resting_cells), "resting FOXP3+ cells (low CCR8 + low MKI67)\n")
} else {
  # 如果没有这些基因，直接用FOXP3+细胞
  resting_cells <- colnames(cds)[foxp3_cells]
  cat("CCR8/MKI67 not in dataset, using all FOXP3+ cells\n")
}

# 以resting细胞的UMAP中心为参考，找最近的实际细胞作为根节点
resting_umap <- umap_coords[resting_cells, ]
resting_center <- colMeans(resting_umap)
dists_to_center <- apply(umap_coords, 1, function(x) sum((x - resting_center)^2))
root_cell <- names(which.min(dists_to_center))

cat("New root cell:", root_cell, "\n")
cat("Root cell subtype:", meta[root_cell, "sub_cell_type"], "\n")

# ===================== 3. 重新计算伪时序 =====================
cat("\n[Step 3] Recomputing pseudotime with new root...\n")
cds <- order_cells(cds, root_cells = root_cell)
cat("Pseudotime recomputed\n")

# 提取新的伪时序数据
pt_data <- data.frame(
  cell_id = colnames(cds),
  pseudotime = pseudotime(cds),
  sub_cell_type = colData(cds)$sub_cell_type,
  response = colData(cds)$response_binary,
  sampleID = colData(cds)$sampleID,
  stringsAsFactors = FALSE
)

# 查看新伪时序下各亚型的分布
cat("\nNew pseudotime distribution by subtype:\n")
pt_summary <- pt_data %>%
  group_by(sub_cell_type) %>%
  summarise(
    median_pt = median(pseudotime),
    mean_pt = mean(pseudotime),
    min_pt = min(pseudotime),
    max_pt = max(pseudotime),
    .groups = "drop"
  )
print(pt_summary)
write.csv(pt_summary, file.path(RESULT_DIR, "pseudotime_subtype_summary.csv"), row.names = FALSE)

# 保存新的伪时序数据
write.csv(pt_data, file.path(RESULT_DIR, "pseudotime_per_cell_new_root.csv"), row.names = FALSE)

# ===================== 4. 可视化新伪时序 =====================
cat("\n[Step 4] Generating visualizations...\n")

# 4.1 UMAP: 新伪时序
p_pt <- plot_cells(cds, color_cells_by = "pseudotime",
                   label_cell_groups = FALSE,
                   cell_size = 0.5,
                   label_leaves = FALSE,
                   label_branch_points = FALSE) +
  scale_color_viridis_c() +
  ggtitle("Pseudotime (New Root)") +
  theme(legend.position = "right")
ggsave(file.path(OUT_DIR, "01_UMAP_pseudotime_newroot.png"), p_pt, width = 8, height = 6, dpi = 300)

# 4.2 轨迹图 + 亚型
p_traj_subtype <- plot_cells(cds, color_cells_by = "sub_cell_type",
                              label_cell_groups = FALSE,
                              cell_size = 0.5,
                              label_leaves = FALSE,
                              label_branch_points = FALSE) +
  scale_color_manual(values = c(
    "CD4T_Treg_FOXP3" = "#1f77b4",
    "CD4T_Treg_CCR8" = "#ff7f0e",
    "CD4T_Treg_MKI67" = "#2ca02c"
  )) +
  ggtitle("Trajectory by Subtype (New Root)") +
  theme(legend.position = "right")
ggsave(file.path(OUT_DIR, "02_trajectory_subtype_newroot.png"), p_traj_subtype, width = 8, height = 6, dpi = 300)

# 4.3 伪时序密度图
p_density <- ggplot(pt_data, aes(x = pseudotime, fill = sub_cell_type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c(
    "CD4T_Treg_FOXP3" = "#1f77b4",
    "CD4T_Treg_CCR8" = "#ff7f0e",
    "CD4T_Treg_MKI67" = "#2ca02c"
  )) +
  labs(x = "Pseudotime", y = "Density", fill = "Subtype") +
  ggtitle("Subtype Distribution along New Pseudotime") +
  theme_bw()
ggsave(file.path(OUT_DIR, "03_subtype_density_newroot.png"), p_density, width = 8, height = 5, dpi = 300)

# 4.4 响应组密度图
p_resp_density <- ggplot(pt_data, aes(x = pseudotime, fill = response)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c(
    "Responder" = "#d62728",
    "Non-responder" = "#9467bd"
  )) +
  labs(x = "Pseudotime", y = "Density", fill = "Response") +
  ggtitle("Response Distribution along New Pseudotime") +
  theme_bw()
ggsave(file.path(OUT_DIR, "04_response_density_newroot.png"), p_resp_density, width = 8, height = 5, dpi = 300)

cat("Visualizations saved to", OUT_DIR, "\n")

# ===================== 5. 定义关键基因集 =====================
cat("\n[Step 5] Defining key gene sets...\n")

# 机械记忆基因（Contractile + Mechanosensing）
mechanical_memory_genes <- c(
  # Contractile核心（细胞骨架收缩）
  "ACTA2", "MYL9", "VCL", "TLN1", "TRPV4", "MOB1B",
  # Extended mechanosensing
  "ITGA1", "DDR1", "DDR2", "SRC", "ITGB1", "CD44", "PTK2", "ZYX", "FLNA",
  # 附加机械力基因
  "ACTG1", "RHOA", "ROCK1", "ROCK2", "MYH9", "MYH10"
)

# Treg关键基因
treg_key_genes <- c(
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "LAG3", "PDCD1",
  "MKI67", "CCR8", "TGFB1", "IL10", "GZMB", "PRF1", "IL2"
)

# 合并
all_key_genes <- unique(c(mechanical_memory_genes, treg_key_genes))
available_key_genes <- all_key_genes[all_key_genes %in% rowData(cds)$gene_short_name]

cat("Mechanical memory genes available:", 
    paste(intersect(mechanical_memory_genes, available_key_genes), collapse=", "), "\n")
cat("Treg key genes available:", 
    paste(intersect(treg_key_genes, available_key_genes), collapse=", "), "\n")
cat("Total key genes to analyze:", length(available_key_genes), "\n")

# ===================== 6. 整体CDS的Moran's I =====================
cat("\n[Step 6] Moran's I on full CDS...\n")

cds_key <- cds[available_key_genes, ]
gene_fits_full <- graph_test(cds_key, 
                              neighbor_graph = "principal_graph", 
                              cores = 1)
gene_fits_full <- gene_fits_full[order(gene_fits_full$morans_I, decreasing = TRUE), ]
gene_fits_full$group <- "Full"
write.csv(gene_fits_full, file.path(RESULT_DIR, "moransI_full_cds.csv"))
cat("Top 10 Moran's I (Full CDS):\n")
print(head(gene_fits_full[, c("gene_short_name", "morans_I", "q_value")], 10))

# ===================== 7. 分组分析：Responder vs Non-responder =====================
cat("\n[Step 7] Grouped analysis: Responder vs Non-responder\n")

# 提取分组细胞
resp_cells <- colnames(cds)[colData(cds)$response_binary == "Responder"]
nonresp_cells <- colnames(cds)[colData(cds)$response_binary == "Non-responder"]
cat("Responder cells:", length(resp_cells), "\n")
cat("Non-responder cells:", length(nonresp_cells), "\n")

# 函数：对子集CDS重新计算graph并做graph_test
analyze_subset <- function(cds_full, cells, group_name, root_cell) {
  cat("\n--- Analyzing", group_name, "---\n")
  
  # 子集化
  cds_sub <- cds_full[, cells]
  cat("Subset:", ncol(cds_sub), "cells\n")
  
  # 过滤零read细胞
  cds_sub <- cds_sub[, Matrix::colSums(SingleCellExperiment::counts(cds_sub)) != 0]
  cat("After zero-filter:", ncol(cds_sub), "cells\n")
  
  # 清除聚类信息并重新聚类（避免子集与全数据集的partitions不匹配）
  cds_sub@clusters[["UMAP"]] <- NULL
  cat("Re-clustering on subset...\n")
  cds_sub <- cluster_cells(cds_sub, reduction_method = "UMAP")
  
  # 重新learn_graph（在子集上）
  cat("Learning graph on subset...\n")
  cds_sub <- learn_graph(cds_sub, use_partition = FALSE)
  
  # 重新order_cells
  # 如果root_cell在子集中，用它；否则找最近的细胞
  if (root_cell %in% colnames(cds_sub)) {
    sub_root <- root_cell
  } else {
    # 找子集中距离原始root_cell UMAP坐标最近的细胞
    umap_full <- reducedDims(cds_full)$UMAP
    root_coord <- umap_full[root_cell, ]
    umap_sub <- reducedDims(cds_sub)$UMAP
    dists <- apply(umap_sub, 1, function(x) sum((x - root_coord)^2))
    sub_root <- names(which.min(dists))
  }
  cat("Root cell for", group_name, ":", sub_root, "\n")
  cds_sub <- order_cells(cds_sub, root_cells = sub_root)
  
  # 计算Moran's I（关键基因）
  cat("Computing Moran's I...\n")
  cds_sub_key <- cds_sub[available_key_genes, ]
  gene_fits_sub <- graph_test(cds_sub_key, 
                              neighbor_graph = "principal_graph", 
                              cores = 1)
  gene_fits_sub <- gene_fits_sub[order(gene_fits_sub$morans_I, decreasing = TRUE), ]
  gene_fits_sub$group <- group_name
  
  # 计算每个基因的pseudotime-expression相关性
  cat("Computing pseudotime-expression correlation...\n")
  pt_sub <- pseudotime(cds_sub)
  corr_results <- data.frame(
    gene = available_key_genes,
    spearman_r = NA,
    spearman_p = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(available_key_genes)) {
    g <- available_key_genes[i]
    expr <- as.numeric(SingleCellExperiment::counts(cds_sub)[g, ])
    # 只在有表达的细胞中计算（避免全0）
    if (sum(expr > 0) > 10) {
      ct <- cor.test(pt_sub, expr, method = "spearman", exact = FALSE)
      corr_results$spearman_r[i] <- ct$estimate
      corr_results$spearman_p[i] <- ct$p.value
    }
  }
  
  # 合并结果
  gene_fits_sub <- merge(gene_fits_sub, corr_results, 
                          by.x = "gene_short_name", by.y = "gene", all = TRUE)
  
  # 保存子集CDS的伪时序
  pt_df <- data.frame(
    cell_id = colnames(cds_sub),
    pseudotime = pt_sub,
    sub_cell_type = colData(cds_sub)$sub_cell_type,
    stringsAsFactors = FALSE
  )
  write.csv(pt_df, file.path(RESULT_DIR, paste0("pseudotime_", group_name, ".csv")), row.names = FALSE)
  
  return(list(gene_fits = gene_fits_sub, cds = cds_sub))
}

# 分别分析两组
resp_result <- analyze_subset(cds, resp_cells, "Responder", root_cell)
nonresp_result <- analyze_subset(cds, nonresp_cells, "Non-responder", root_cell)

# 保存分组Moran's I
write.csv(resp_result$gene_fits, file.path(RESULT_DIR, "moransI_responder.csv"), row.names = FALSE)
write.csv(nonresp_result$gene_fits, file.path(RESULT_DIR, "moransI_nonresponder.csv"), row.names = FALSE)

# ===================== 8. 对比分析 =====================
cat("\n[Step 8] Cross-group comparison\n")

# 合并对比
df_resp <- resp_result$gene_fits[, c("gene_short_name", "morans_I", "q_value", "spearman_r", "spearman_p")]
df_nonresp <- nonresp_result$gene_fits[, c("gene_short_name", "morans_I", "q_value", "spearman_r", "spearman_p")]

colnames(df_resp) <- c("gene", "moransI_resp", "q_resp", "spearman_r_resp", "spearman_p_resp")
colnames(df_nonresp) <- c("gene", "moransI_nonresp", "q_nonresp", "spearman_r_nonresp", "spearman_p_nonresp")

comparison <- merge(df_resp, df_nonresp, by = "gene")

# 计算差异
comparison$moransI_diff <- comparison$moransI_nonresp - comparison$moransI_resp
comparison$spearman_diff <- comparison$spearman_r_nonresp - comparison$spearman_r_resp

# Fisher Z-test for Spearman correlation difference
comparison$fisher_z <- NA
comparison$fisher_p <- NA

n_resp <- length(resp_cells)
n_nonresp <- length(nonresp_cells)

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

# FDR校正
comparison$fisher_fdr <- p.adjust(comparison$fisher_p, method = "fdr")

# 排序：按Moran's I差异（Non-responder - Responder）降序
comparison <- comparison[order(comparison$moransI_diff, decreasing = TRUE), ]

write.csv(comparison, file.path(RESULT_DIR, "moransI_comparison.csv"), row.names = FALSE)

# 打印关键结果
cat("\n--- Genes with higher Moran's I in Non-responder (Top 15) ---\n")
print(head(comparison[, c("gene", "moransI_resp", "moransI_nonresp", "moransI_diff", "fisher_p")], 15))

cat("\n--- Genes with higher Moran's I in Responder (Top 10) ---\n")
print(tail(comparison[, c("gene", "moransI_resp", "moransI_nonresp", "moransI_diff", "fisher_p")], 10))

# ===================== 9. 按基因类别汇总 =====================
cat("\n[Step 9] Category-level summary\n")

comparison$category <- ifelse(comparison$gene %in% mechanical_memory_genes, 
                               "Mechanical_Memory", 
                               ifelse(comparison$gene %in% treg_key_genes, 
                                      "Treg_Key", "Other"))

cat_summary <- comparison %>%
  group_by(category) %>%
  summarise(
    n_genes = n(),
    mean_moransI_resp = mean(moransI_resp, na.rm = TRUE),
    mean_moransI_nonresp = mean(moransI_nonresp, na.rm = TRUE),
    mean_moransI_diff = mean(moransI_diff, na.rm = TRUE),
    mean_spearman_r_resp = mean(spearman_r_resp, na.rm = TRUE),
    mean_spearman_r_nonresp = mean(spearman_r_nonresp, na.rm = TRUE),
    mean_spearman_diff = mean(spearman_diff, na.rm = TRUE),
    .groups = "drop"
  )
print(cat_summary)
write.csv(cat_summary, file.path(RESULT_DIR, "category_summary.csv"), row.names = FALSE)

# ===================== 10. 可视化对比结果 =====================
cat("\n[Step 10] Generating comparison plots...\n")

# 10.1 Moran's I 分组对比散点图
p_moran_scatter <- ggplot(comparison, aes(x = moransI_resp, y = moransI_nonresp, color = category)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_text(data = comparison[abs(comparison$moransI_diff) > 0.02, ],
            aes(label = gene), vjust = -1, size = 3) +
  scale_color_manual(values = c("Mechanical_Memory" = "#E41A1C", 
                                 "Treg_Key" = "#377EB8", 
                                 "Other" = "#999999")) +
  labs(x = "Moran's I (Responder)", 
       y = "Moran's I (Non-responder)",
       color = "Gene Category") +
  ggtitle("Moran's I: Non-responder vs Responder") +
  theme_bw()
ggsave(file.path(OUT_DIR, "05_moransI_scatter.png"), p_moran_scatter, width = 8, height = 7, dpi = 300)

# 10.2 Moran's I 差异条形图
comparison_plot <- comparison %>%
  mutate(gene = reorder(gene, moransI_diff))

p_moran_diff <- ggplot(comparison_plot, aes(x = gene, y = moransI_diff, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Mechanical_Memory" = "#E41A1C", 
                                "Treg_Key" = "#377EB8", 
                                "Other" = "#999999")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Gene", y = "Δ Moran's I (Non-responder - Responder)", fill = "Category") +
  ggtitle("Moran's I Difference: Non-responder vs Responder") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))
ggsave(file.path(OUT_DIR, "06_moransI_diff_barplot.png"), p_moran_diff, width = 8, height = 10, dpi = 300)

# 10.3 Spearman相关分组对比
p_spear_scatter <- ggplot(comparison, aes(x = spearman_r_resp, y = spearman_r_nonresp, color = category)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_text(data = comparison[abs(comparison$spearman_diff) > 0.02 & !is.na(comparison$spearman_diff), ],
            aes(label = gene), vjust = -1, size = 3) +
  scale_color_manual(values = c("Mechanical_Memory" = "#E41A1C", 
                                 "Treg_Key" = "#377EB8", 
                                 "Other" = "#999999")) +
  labs(x = "Spearman r (Responder)", 
       y = "Spearman r (Non-responder)",
       color = "Gene Category") +
  ggtitle("Pseudotime-Expression Correlation: Non-responder vs Responder") +
  theme_bw()
ggsave(file.path(OUT_DIR, "07_spearman_scatter.png"), p_spear_scatter, width = 8, height = 7, dpi = 300)

# 10.4 机械记忆基因专题图
mech_genes <- comparison[comparison$category == "Mechanical_Memory", ]
if (nrow(mech_genes) > 0) {
  mech_long <- melt(mech_genes[, c("gene", "moransI_resp", "moransI_nonresp")], 
                    id.vars = "gene", 
                    variable.name = "group", 
                    value.name = "morans_I")
  mech_long$group <- ifelse(mech_long$group == "moransI_resp", "Responder", "Non-responder")
  
  p_mech_moran <- ggplot(mech_long, aes(x = reorder(gene, morans_I), y = morans_I, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Responder" = "#d62728", "Non-responder" = "#9467bd")) +
    labs(x = "Gene", y = "Moran's I", fill = "Group") +
    ggtitle("Mechanical Memory Genes: Moran's I by Group") +
    coord_flip() +
    theme_bw()
  ggsave(file.path(OUT_DIR, "08_mechanical_memory_moransI.png"), p_mech_moran, width = 8, height = 6, dpi = 300)
}

# 10.5 Treg关键基因专题图
treg_genes <- comparison[comparison$category == "Treg_Key", ]
if (nrow(treg_genes) > 0) {
  treg_long <- melt(treg_genes[, c("gene", "moransI_resp", "moransI_nonresp")], 
                    id.vars = "gene", 
                    variable.name = "group", 
                    value.name = "morans_I")
  treg_long$group <- ifelse(treg_long$group == "moransI_resp", "Responder", "Non-responder")
  
  p_treg_moran <- ggplot(treg_long, aes(x = reorder(gene, morans_I), y = morans_I, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Responder" = "#d62728", "Non-responder" = "#9467bd")) +
    labs(x = "Gene", y = "Moran's I", fill = "Group") +
    ggtitle("Treg Key Genes: Moran's I by Group") +
    coord_flip() +
    theme_bw()
  ggsave(file.path(OUT_DIR, "09_treg_key_moransI.png"), p_treg_moran, width = 8, height = 6, dpi = 300)
}

# ===================== 11. 伪时序基因表达趋势图 =====================
cat("\n[Step 11] Pseudotime gene expression trends...\n")

# 选取最关键的机械记忆基因做LOESS平滑曲线
priority_mech <- c("ACTA2", "MYL9", "VCL", "TLN1", "TRPV4", "ITGA1", "DDR1", "DDR2", "SRC")
priority_mech <- priority_mech[priority_mech %in% available_key_genes]

# 合并两组数据用于绘图
pt_resp <- read.csv(file.path(RESULT_DIR, "pseudotime_Responder.csv"), stringsAsFactors = FALSE)
pt_nonresp <- read.csv(file.path(RESULT_DIR, "pseudotime_Non-responder.csv"), stringsAsFactors = FALSE)
pt_resp$group <- "Responder"
pt_nonresp$group <- "Non-responder"

# 对每个优先基因，提取表达并绘制LOESS曲线
for (g in priority_mech) {
  cat("Plotting", g, "...\n")
  
  # 获取表达量（在整个CDS上）
  expr <- as.numeric(SingleCellExperiment::counts(cds)[g, ])
  expr_df <- data.frame(
    cell_id = colnames(cds),
    expression = expr,
    pseudotime = pseudotime(cds),
    response = colData(cds)$response_binary,
    stringsAsFactors = FALSE
  )
  
  # 只保留有表达的细胞用于绘图（减少点密度）
  expr_df <- expr_df[expr_df$expression > 0, ]
  
  if (nrow(expr_df) > 50) {
    p_trend <- ggplot(expr_df, aes(x = pseudotime, y = expression, color = response)) +
      geom_point(alpha = 0.1, size = 0.5) +
      geom_smooth(method = "loess", se = TRUE, aes(fill = response), alpha = 0.2) +
      scale_color_manual(values = c("Responder" = "#d62728", "Non-responder" = "#9467bd")) +
      scale_fill_manual(values = c("Responder" = "#d62728", "Non-responder" = "#9467bd")) +
      labs(x = "Pseudotime", y = "Expression", color = "Group", fill = "Group") +
      ggtitle(paste(g, "Expression along Pseudotime")) +
      theme_bw()
    ggsave(file.path(OUT_DIR, sprintf("10_trend_%s.png", g)), p_trend, width = 8, height = 5, dpi = 300)
  }
}

# ===================== 12. 保存最终CDS =====================
cat("\n[Step 12] Saving results...\n")
saveRDS(cds, file.path(RESULT_DIR, "treg_monocle3_cds_newroot.rds"))
saveRDS(resp_result$cds, file.path(RESULT_DIR, "treg_monocle3_responder.rds"))
saveRDS(nonresp_result$cds, file.path(RESULT_DIR, "treg_monocle3_nonresponder.rds"))

# ===================== 13. 总结报告 =====================
cat("\n============================================================\n")
cat("Analysis Complete!\n")
cat("============================================================\n")
cat("Output directory:", OUT_DIR, "\n")
cat("Results directory:", RESULT_DIR, "\n")
cat("\nKey files:\n")
cat("  - leaf_composition.csv\n")
cat("  - pseudotime_subtype_summary.csv\n")
cat("  - moransI_full_cds.csv\n")
cat("  - moransI_responder.csv\n")
cat("  - moransI_nonresponder.csv\n")
cat("  - moransI_comparison.csv\n")
cat("  - category_summary.csv\n")
cat("\nNew root cell:", root_cell, "\n")
cat("Root cell subtype:", meta[root_cell, "sub_cell_type"], "\n")
cat("\nPseudotime ranges (new root):\n")
print(pt_summary)
cat("\nHypothesis test (Non-responder > Responder in Moran's I):\n")
cat("  Mean Moran's I diff (Non-resp - Resp):", mean(comparison$moransI_diff, na.rm = TRUE), "\n")
cat("  Mechanical memory genes mean diff:", mean(comparison$moransI_diff[comparison$category == "Mechanical_Memory"], na.rm = TRUE), "\n")
cat("  Treg key genes mean diff:", mean(comparison$moransI_diff[comparison$category == "Treg_Key"], na.rm = TRUE), "\n")
