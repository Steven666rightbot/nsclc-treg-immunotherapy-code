#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# GSE243013 - CCR8+/MKI67+ Treg 统一伪时序分析（不分组）
# 
# 目标：
#   1. 补充 ACTG2
#   2. 强弱响应合在一起，根据拟时序曲线推断基因调控层级/上下游关系
#   3. 计算各基因动态特征（峰值位置、变化方向、拐点）
#   4. 基因间伪时序动态相关性矩阵

library(monocle3)
library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)

# ===================== 路径配置 =====================
CDS_PATH <- "D:/Research/tomato/results/monocle3_ccr8_mki67/cds_ccr8_mki67_full.rds"
OUT_DIR <- "D:/Research/tomato/figures/monocle3_unified"
RESULT_DIR <- "D:/Research/tomato/results/monocle3_unified"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("Unified Pseudotime Analysis: CCR8+/MKI67+ (All Cells)\n")
cat("============================================================\n")

# ===================== 1. 加载完整CDS =====================
cat("\n[Step 1] Loading full CDS...\n")
cds <- readRDS(CDS_PATH)
cat("CDS:", ncol(cds), "cells x", nrow(cds), "genes\n")
cat("Subtypes:", table(colData(cds)$sub_cell_type), "\n")
cat("Response:", table(colData(cds)$response_binary), "\n")

# ===================== 2. 定义基因集（含ACTG2） =====================
cat("\n[Step 2] Defining gene set (including ACTG2)...\n")

mechanical_memory_genes <- c(
  "ACTA2", "MYL9", "VCL", "TLN1", "TRPV4", "MOB1B",
  "ITGA1", "DDR1", "DDR2", "SRC", "ITGB1", "CD44", "PTK2", "ZYX", "FLNA",
  "ACTG1", "ACTG2", "RHOA", "ROCK1", "ROCK2", "MYH9", "MYH10", "PFN1", "CFL1"
)

treg_key_genes <- c(
  "IL2RA", "CTLA4", "TIGIT", "LAG3", "PDCD1",
  "MKI67", "CCR8", "TGFB1", "IL10", "GZMB", "PRF1", "IL2", "FOXP3"
)

all_key_genes <- unique(c(mechanical_memory_genes, treg_key_genes))
available_key_genes <- all_key_genes[all_key_genes %in% rowData(cds)$gene_short_name]
mech_available <- intersect(mechanical_memory_genes, available_key_genes)
treg_available <- intersect(treg_key_genes, available_key_genes)

cat("Mechanical memory genes:", length(mech_available), "/", length(mechanical_memory_genes), "\n")
cat("  ", paste(mech_available, collapse=", "), "\n")
cat("Treg key genes:", length(treg_available), "/", length(treg_key_genes), "\n")
cat("  ", paste(treg_available, collapse=", "), "\n")

# ===================== 3. Moran's I =====================
cat("\n[Step 3] Computing Moran's I...\n")
cds_key <- cds[available_key_genes, ]
gene_fits <- graph_test(cds_key, neighbor_graph = "principal_graph", cores = 1)
gene_fits <- gene_fits[order(gene_fits$morans_I, decreasing = TRUE), ]
write.csv(gene_fits, file.path(RESULT_DIR, "moransI_unified.csv"))
cat("Top Moran's I genes:\n")
print(gene_fits[, c("gene_short_name", "morans_I", "q_value")])

# ===================== 4. 提取伪时序表达矩阵 =====================
cat("\n[Step 4] Extracting pseudotime expression profiles...\n")

pt_vec <- pseudotime(cds)
expr_mat <- as.matrix(SingleCellExperiment::counts(cds)[available_key_genes, ])

# 为每个基因计算伪时序上的 LOESS 拟合曲线（等间距采样）
pt_grid <- seq(min(pt_vec), max(pt_vec), length.out = 200)

loess_profiles <- data.frame(pseudotime = pt_grid)
peak_positions <- data.frame(
  gene = available_key_genes,
  peak_pt = NA,
  valley_pt = NA,
  trend_type = NA,
  early_expr = NA,
  late_expr = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(available_key_genes)) {
  g <- available_key_genes[i]
  expr <- expr_mat[g, ]
  
  # LOESS 拟合
  df <- data.frame(pt = pt_vec, expr = expr)
  fit <- loess(expr ~ pt, data = df, span = 0.5, degree = 2)
  pred <- predict(fit, newdata = data.frame(pt = pt_grid))
  loess_profiles[[g]] <- pred
  
  # 找峰值和谷值
  peak_idx <- which.max(pred)
  valley_idx <- which.min(pred)
  peak_positions$peak_pt[i] <- pt_grid[peak_idx]
  peak_positions$valley_pt[i] <- pt_grid[valley_idx]
  peak_positions$early_expr[i] <- pred[1]
  peak_positions$late_expr[i] <- pred[length(pred)]
  
  # 分类趋势
  early <- mean(pred[1:20])
  late <- mean(pred[181:200])
  mid <- mean(pred[90:110])
  
  if (early > mid && late > mid) {
    trend <- "U-shaped"
  } else if (early < mid && late < mid) {
    trend <- "Inverted-U"
  } else if (early > late) {
    trend <- "Decreasing"
  } else if (early < late) {
    trend <- "Increasing"
  } else {
    trend <- "Flat"
  }
  peak_positions$trend_type[i] <- trend
}

write.csv(loess_profiles, file.path(RESULT_DIR, "loess_profiles.csv"), row.names = FALSE)
write.csv(peak_positions, file.path(RESULT_DIR, "peak_positions.csv"), row.names = FALSE)

cat("\nTrend classification:\n")
print(table(peak_positions$trend_type))

# ===================== 5. 基因间伪时序动态相关性 =====================
cat("\n[Step 5] Gene-gene pseudotime dynamic correlation...\n")

# 用 LOESS 拟合值计算基因间相关性（避免单细胞噪声）
loess_mat <- as.matrix(loess_profiles[, -1])
gene_corr <- cor(loess_mat, method = "pearson")
write.csv(gene_corr, file.path(RESULT_DIR, "gene_gene_loess_correlation.csv"))

# 找出高相关的基因对
corr_pairs <- data.frame(
  gene1 = character(),
  gene2 = character(),
  correlation = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(nrow(gene_corr)-1)) {
  for (j in (i+1):ncol(gene_corr)) {
    r <- gene_corr[i, j]
    if (abs(r) > 0.7) {
      corr_pairs <- rbind(corr_pairs, data.frame(
        gene1 = rownames(gene_corr)[i],
        gene2 = colnames(gene_corr)[j],
        correlation = r
      ))
    }
  }
}
corr_pairs <- corr_pairs[order(abs(corr_pairs$correlation), decreasing = TRUE), ]
write.csv(corr_pairs, file.path(RESULT_DIR, "high_corr_gene_pairs.csv"), row.names = FALSE)

cat("High-correlation gene pairs (|r| > 0.7):\n")
print(head(corr_pairs, 20))

# ===================== 6. 可视化：所有基因的 LOESS 曲线叠合 =====================
cat("\n[Step 6] Generating unified trend plots...\n")

# 6.1 机械记忆基因统一叠合图
mech_long <- melt(loess_profiles, id.vars = "pseudotime", 
                  measure.vars = mech_available,
                  variable.name = "gene", value.name = "expression")

p_mech_all <- ggplot(mech_long, aes(x = pseudotime, y = expression, color = gene)) +
  geom_line(size = 0.8) +
  labs(x = "Pseudotime", y = "LOESS-fitted Expression", color = "Gene") +
  ggtitle("Mechanical Memory Genes - Unified Pseudotime Trends") +
  theme_bw()
ggsave(file.path(OUT_DIR, "01_unified_mech_trends.png"), p_mech_all, width = 10, height = 6, dpi = 300)

# 6.2 Treg关键基因统一叠合图
treg_long <- melt(loess_profiles, id.vars = "pseudotime",
                  measure.vars = treg_available,
                  variable.name = "gene", value.name = "expression")

p_treg_all <- ggplot(treg_long, aes(x = pseudotime, y = expression, color = gene)) +
  geom_line(size = 0.8) +
  labs(x = "Pseudotime", y = "LOESS-fitted Expression", color = "Gene") +
  ggtitle("Treg Key Genes - Unified Pseudotime Trends") +
  theme_bw()
ggsave(file.path(OUT_DIR, "02_unified_treg_trends.png"), p_treg_all, width = 10, height = 6, dpi = 300)

# 6.3 分面展示：每个基因单独小图
mech_long$category <- "Mechanical_Memory"
treg_long$category <- "Treg_Key"
all_long <- rbind(mech_long, treg_long)

p_facet <- ggplot(all_long, aes(x = pseudotime, y = expression)) +
  geom_line(color = "steelblue", size = 0.7) +
  facet_wrap(~ gene, scales = "free_y", ncol = 6) +
  labs(x = "Pseudotime", y = "Expression") +
  ggtitle("Individual Gene Trends along Pseudotime") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 8))
ggsave(file.path(OUT_DIR, "03_facet_all_genes.png"), p_facet, width = 14, height = 10, dpi = 300)

# ===================== 7. 热图：基因-伪时序表达矩阵 =====================
cat("\n[Step 7] Pseudotime expression heatmap...\n")

# 对 LOESS 曲线进行 z-score 标准化，方便比较形状
loess_z <- apply(loess_mat, 2, function(x) (x - mean(x)) / sd(x))
loess_z_df <- as.data.frame(loess_z)
loess_z_df$pseudotime <- pt_grid

# 按峰值位置排序基因
peak_order <- peak_positions$gene[order(peak_positions$peak_pt)]

heatmap_long <- melt(loess_z_df, id.vars = "pseudotime", variable.name = "gene", value.name = "zscore")
heatmap_long$gene <- factor(heatmap_long$gene, levels = peak_order)

p_heatmap <- ggplot(heatmap_long, aes(x = pseudotime, y = gene, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Pseudotime", y = "Gene", fill = "Z-score") +
  ggtitle("Gene Expression Dynamics Heatmap\n(ordered by peak position)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))
ggsave(file.path(OUT_DIR, "04_pseudotime_heatmap.png"), p_heatmap, width = 10, height = 8, dpi = 300)

# ===================== 8. 上下游层级推断 =====================
cat("\n[Step 8] Inferring regulatory hierarchy...\n")

# 策略：
# 1. 早期高表达 → 可能是上游调控因子
# 2. 晚期高表达 → 可能是下游效应分子
# 3. U型 → 可能是双相调控
# 4. 结合 Moran's I：高 Moran's I = 强轨迹约束 = 可能是核心调控节点

hierarchy <- merge(
  peak_positions,
  gene_fits[, c("gene_short_name", "morans_I")],
  by.x = "gene", by.y = "gene_short_name", all.x = TRUE
)

# 计算综合评分：
# - 早期高表达且 Moran's I 高 → 上游候选
# - 晚期高表达 → 下游候选
hierarchy$upstream_score <- hierarchy$early_expr * hierarchy$morans_I
hierarchy$downstream_score <- hierarchy$late_expr * hierarchy$morans_I

hierarchy <- hierarchy[order(hierarchy$morans_I, decreasing = TRUE), ]
write.csv(hierarchy, file.path(RESULT_DIR, "regulatory_hierarchy.csv"), row.names = FALSE)

cat("\nTop upstream candidates (early high + high Moran's I):\n")
up_candidates <- hierarchy[order(hierarchy$upstream_score, decreasing = TRUE), ]
print(head(up_candidates[, c("gene", "trend_type", "early_expr", "late_expr", "morans_I", "upstream_score")], 10))

cat("\nTop downstream candidates (late high + high Moran's I):\n")
down_candidates <- hierarchy[order(hierarchy$downstream_score, decreasing = TRUE), ]
print(head(down_candidates[, c("gene", "trend_type", "early_expr", "late_expr", "morans_I", "downstream_score")], 10))

# ===================== 9. 峰值位置瀑布图 =====================
cat("\n[Step 9] Peak position waterfall plot...\n")

peak_plot <- peak_positions
peak_plot$gene <- factor(peak_plot$gene, levels = peak_plot$gene[order(peak_plot$peak_pt)])
peak_plot$category <- ifelse(peak_plot$gene %in% mech_available, "Mechanical_Memory", "Treg_Key")

p_peak <- ggplot(peak_plot, aes(x = gene, y = peak_pt, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Mechanical_Memory" = "#E41A1C", "Treg_Key" = "#377EB8")) +
  labs(x = "Gene", y = "Peak Pseudotime Position") +
  ggtitle("Gene Expression Peak Position along Pseudotime") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))
ggsave(file.path(OUT_DIR, "05_peak_position_waterfall.png"), p_peak, width = 8, height = 10, dpi = 300)

# ===================== 10. 基因相关性网络图 =====================
cat("\n[Step 10] Gene correlation network...\n")

# 只保留机械记忆基因的相关性网络
mech_corr <- gene_corr[mech_available, mech_available]
mech_corr_pairs <- data.frame(
  gene1 = character(), gene2 = character(), r = numeric(), stringsAsFactors = FALSE
)
for (i in 1:(nrow(mech_corr)-1)) {
  for (j in (i+1):ncol(mech_corr)) {
    r <- mech_corr[i, j]
    if (abs(r) > 0.5) {
      mech_corr_pairs <- rbind(mech_corr_pairs, data.frame(
        gene1 = rownames(mech_corr)[i], gene2 = colnames(mech_corr)[j], r = r
      ))
    }
  }
}

if (nrow(mech_corr_pairs) > 0) {
  mech_corr_pairs <- mech_corr_pairs[order(abs(mech_corr_pairs$r), decreasing = TRUE), ]
  
  p_corr_net <- ggplot(mech_corr_pairs, aes(x = gene1, y = gene2, fill = r)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(r, 2)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-1, 1)) +
    labs(x = "", y = "", fill = "Correlation") +
    ggtitle("Mechanical Memory Genes\nPseudotime Dynamic Correlation") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8))
  ggsave(file.path(OUT_DIR, "06_mech_corr_network.png"), p_corr_net, width = 8, height = 7, dpi = 300)
}

# ===================== 11. 与PPI网络交叉验证 =====================
cat("\n[Step 11] Cross-reference with PPI network...\n")

# 读取PPI边
ppi_edges <- read.csv("D:/Research/tomato/results/ppi/ppi_edges.csv")

# 找同时有高LOESS相关性和高PPI score的基因对
ppi_mech <- ppi_edges[
  (ppi_edges$source %in% mech_available & ppi_edges$target %in% mech_available),
]

if (nrow(ppi_mech) > 0) {
  ppi_mech$loess_r <- NA
  for (i in 1:nrow(ppi_mech)) {
    g1 <- ppi_mech$source[i]
    g2 <- ppi_mech$target[i]
    if (g1 %in% rownames(gene_corr) && g2 %in% rownames(gene_corr)) {
      ppi_mech$loess_r[i] <- gene_corr[g1, g2]
    }
  }
  ppi_mech <- ppi_mech[order(ppi_mech$score, decreasing = TRUE), ]
  write.csv(ppi_mech, file.path(RESULT_DIR, "ppi_vs_pseudotime_correlation.csv"), row.names = FALSE)
  
  cat("PPI edges among mechanical memory genes:\n")
  print(ppi_mech[, c("source", "target", "score", "loess_r")])
}

# ===================== 12. 总结报告 =====================
cat("\n============================================================\n")
cat("Analysis Complete!\n")
cat("============================================================\n")
cat("Output:", OUT_DIR, "\n")
cat("Results:", RESULT_DIR, "\n")

cat("\nKey findings:\n")
cat("- Total cells:", ncol(cds), "\n")
cat("- Genes analyzed:", length(available_key_genes), "\n")
cat("- Trend types:\n")
print(table(peak_positions$trend_type))
cat("\n- Top Moran's I:\n")
print(head(gene_fits[, c("gene_short_name", "morans_I")], 8))
"
