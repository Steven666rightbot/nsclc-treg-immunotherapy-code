# =============================================================================
# CellChat Virtual Knockout Analysis - GSE207422
# 
# 正宗 R 做法：
#   1. 读取完整表达矩阵
#   2. 构建 CellChat 对象（Strong / Weak）
#   3. 虚拟敲除：将目标基因表达设为 0
#   4. 重新 computeCommunProb
#   5. 对比 KO 前后 Fibroblast → Treg 通信概率
#
# 运行前请先执行 prepare_gse207422_for_cellchat.py 生成输入数据
# =============================================================================

library(CellChat)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(Matrix)
library(data.table)

# ======================== 参数配置 ========================
INPUT_DIR  <- "D:/Research/tomato/data/gse207422_cellchat_input"
OUTPUT_DIR <- "D:/Research/tomato/results/cellchat_ko"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 选择要分析的组
RUN_STRONG <- TRUE   # Strong 响应者基础分析
RUN_WEAK   <- TRUE   # Weak 响应者基础分析（核心）

# KO 条件列表
KO_LIST <- list(
  "KO_COL1A1"        = c("COL1A1"),
  "KO_COL1A2"        = c("COL1A2"),
  "KO_COL1A1_COL1A2" = c("COL1A1", "COL1A2"),
  "KO_FN1"           = c("FN1"),
  "KO_ALL_ECM"       = c("COL1A1", "COL1A2", "FN1")
)

# 如果 RDS 已存在则跳过（支持断点续跑）
SKIP_EXISTING_RDS <- TRUE

# ======================== 辅助函数 ========================

#' 读取预处理后的 MTX 输入数据
read_cellchat_input <- function(prefix) {
  cat(sprintf("\n[%s] 读取输入数据...\n", prefix))
  
  # 读取 MTX (Matrix Market format)
  mtx_file <- file.path(INPUT_DIR, sprintf("%s_matrix.mtx", prefix))
  expr_matrix <- readMM(mtx_file)
  
  # 读取基因名
  genes_file <- file.path(INPUT_DIR, sprintf("%s_genes.tsv", prefix))
  genes <- fread(genes_file, header = FALSE)$V1
  
  # 读取细胞名
  barcodes_file <- file.path(INPUT_DIR, sprintf("%s_barcodes.tsv", prefix))
  barcodes <- fread(barcodes_file, header = FALSE)$V1
  
  # 读取 metadata
  meta_file <- file.path(INPUT_DIR, sprintf("%s_metadata.csv", prefix))
  meta <- fread(meta_file)
  
  # 设置行列名
  rownames(expr_matrix) <- genes
  colnames(expr_matrix) <- barcodes
  
  # 转为 dgCMatrix (CellChat 需要)
  expr_matrix <- as(expr_matrix, "dgCMatrix")
  
  cat(sprintf("[%s] 矩阵: %d genes x %d cells\n", prefix, nrow(expr_matrix), ncol(expr_matrix)))
  return(list(expr = expr_matrix, meta = meta))
}

#' 运行完整 CellChat pipeline
run_cellchat_pipeline <- function(expr_matrix, meta, group.by = "cell_type") {
  cat("  -> createCellChat...\n")
  cellchat <- createCellChat(object = expr_matrix, meta = meta, group.by = group.by)
  
  cat("  -> set database...\n")
  cellchat@DB <- CellChatDB.human
  
  cat("  -> subsetData...\n")
  cellchat <- subsetData(cellchat)
  
  cat("  -> identifyOverExpressedGenes...\n")
  cellchat <- identifyOverExpressedGenes(cellchat)
  
  cat("  -> identifyOverExpressedInteractions...\n")
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cat("  -> computeCommunProb (最耗时，请耐心等待)...\n")
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  
  cat("  -> filterCommunication...\n")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  cat("  -> computeCommunProbPathway...\n")
  cellchat <- computeCommunProbPathway(cellchat)
  
  cat("  -> aggregateNet...\n")
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

#' 虚拟敲除：将指定基因表达设为 0
virtual_ko <- function(expr_matrix, ko_genes) {
  expr_ko <- expr_matrix
  for (g in ko_genes) {
    if (g %in% rownames(expr_ko)) {
      expr_ko[g, ] <- 0
      cat(sprintf("    KO: %s (已设为0)\n", g))
    } else {
      cat(sprintf("    警告: %s 不在表达矩阵中，跳过\n", g))
    }
  }
  return(expr_ko)
}

#' 提取 Fibroblast -> Treg 通信
dump_fib_treg <- function(cellchat, label) {
  df <- subsetCommunication(cellchat)
  ft <- df[df$source == "Fibroblast" & df$target == "Treg", ]
  cat(sprintf("    %s Fibroblast->Treg 互作: %d 条\n", label, nrow(ft)))
  if (nrow(ft) > 0) {
    cat(sprintf("    %s 总 prob: %.6e\n", label, sum(ft$prob)))
  }
  return(ft)
}

#' 保存 CellChat 结果（CSV + RDS）
save_cellchat_results <- function(cellchat, prefix) {
  rds_file <- file.path(OUTPUT_DIR, sprintf("cellchat_%s.rds", prefix))
  csv_file <- file.path(OUTPUT_DIR, sprintf("interactions_%s.csv", prefix))
  
  saveRDS(cellchat, rds_file)
  write.csv(subsetCommunication(cellchat), csv_file, row.names = FALSE)
  
  cat(sprintf("  已保存: %s\n", rds_file))
  cat(sprintf("  已保存: %s\n", csv_file))
}

# ======================== 主流程 ========================

cat("=" , rep("=", 58), "\n", sep = "")
cat("CellChat Virtual Knockout Analysis\n")
cat("=" , rep("=", 58), "\n", sep = "")

# --------------------------------------------------------
# 1. 读取数据
# --------------------------------------------------------
if (RUN_STRONG) {
  strong_data <- read_cellchat_input("strong")
}
if (RUN_WEAK) {
  weak_data <- read_cellchat_input("weak")
}

# --------------------------------------------------------
# 2. 基础 CellChat 分析
# --------------------------------------------------------
results <- list()

if (RUN_STRONG) {
  cat("\n--- Strong 响应者 基础分析 ---\n")
  rds_strong <- file.path(OUTPUT_DIR, "cellchat_strong.rds")
  
  if (SKIP_EXISTING_RDS && file.exists(rds_strong)) {
    cat("RDS 已存在，跳过计算，直接读取...\n")
    cellchat_strong <- readRDS(rds_strong)
  } else {
    cellchat_strong <- run_cellchat_pipeline(strong_data$expr, strong_data$meta)
    save_cellchat_results(cellchat_strong, "strong")
  }
  
  results[["strong_base"]] <- list(
    obj = cellchat_strong,
    fib_treg = dump_fib_treg(cellchat_strong, "Strong")
  )
}

if (RUN_WEAK) {
  cat("\n--- Weak 响应者 基础分析 ---\n")
  rds_weak <- file.path(OUTPUT_DIR, "cellchat_weak.rds")
  
  if (SKIP_EXISTING_RDS && file.exists(rds_weak)) {
    cat("RDS 已存在，跳过计算，直接读取...\n")
    cellchat_weak <- readRDS(rds_weak)
  } else {
    cellchat_weak <- run_cellchat_pipeline(weak_data$expr, weak_data$meta)
    save_cellchat_results(cellchat_weak, "weak")
  }
  
  results[["weak_base"]] <- list(
    obj = cellchat_weak,
    fib_treg = dump_fib_treg(cellchat_weak, "Weak")
  )
}

# --------------------------------------------------------
# 3. 虚拟敲除分析
# --------------------------------------------------------
cat("\n", rep("-", 60), "\n", sep = "")
cat("开始虚拟敲除分析\n")
cat(rep("-", 60), "\n", sep = "")

ko_summary <- data.frame(
  KO = character(),
  Group = character(),
  Base_Interactions = integer(),
  KO_Interactions = integer(),
  Base_Prob = numeric(),
  KO_Prob = numeric(),
  Retained_Pct = numeric(),
  Lost_Pct = numeric(),
  stringsAsFactors = FALSE
)

for (ko_name in names(KO_LIST)) {
  ko_genes <- KO_LIST[[ko_name]]
  cat(sprintf("\n========== %s: %s ==========\n", ko_name, paste(ko_genes, collapse = " + ")))
  
  # --- Strong KO ---
  if (RUN_STRONG) {
    rds_ko_strong <- file.path(OUTPUT_DIR, sprintf("cellchat_strong_%s.rds", ko_name))
    
    if (SKIP_EXISTING_RDS && file.exists(rds_ko_strong)) {
      cat("[Strong] RDS 已存在，跳过...\n")
      cellchat_strong_ko <- readRDS(rds_ko_strong)
    } else {
      cat("[Strong] 执行 KO...\n")
      expr_strong_ko <- virtual_ko(strong_data$expr, ko_genes)
      cellchat_strong_ko <- run_cellchat_pipeline(expr_strong_ko, strong_data$meta)
      save_cellchat_results(cellchat_strong_ko, sprintf("strong_%s", ko_name))
    }
    
    ft_strong_ko <- dump_fib_treg(cellchat_strong_ko, "Strong KO")
    ft_strong_base <- results[["strong_base"]]$fib_treg
    
    base_prob <- if (nrow(ft_strong_base) > 0) sum(ft_strong_base$prob) else 0
    ko_prob   <- if (nrow(ft_strong_ko) > 0) sum(ft_strong_ko$prob) else 0
    retained  <- if (base_prob > 0) round(ko_prob / base_prob * 100, 2) else NA
    lost      <- if (base_prob > 0) round(100 - retained, 2) else NA
    
    ko_summary <- rbind(ko_summary, data.frame(
      KO = ko_name,
      Group = "Strong",
      Base_Interactions = nrow(ft_strong_base),
      KO_Interactions = nrow(ft_strong_ko),
      Base_Prob = base_prob,
      KO_Prob = ko_prob,
      Retained_Pct = retained,
      Lost_Pct = lost,
      stringsAsFactors = FALSE
    ))
  }
  
  # --- Weak KO ---
  if (RUN_WEAK) {
    rds_ko_weak <- file.path(OUTPUT_DIR, sprintf("cellchat_weak_%s.rds", ko_name))
    
    if (SKIP_EXISTING_RDS && file.exists(rds_ko_weak)) {
      cat("[Weak] RDS 已存在，跳过...\n")
      cellchat_weak_ko <- readRDS(rds_ko_weak)
    } else {
      cat("[Weak] 执行 KO...\n")
      expr_weak_ko <- virtual_ko(weak_data$expr, ko_genes)
      cellchat_weak_ko <- run_cellchat_pipeline(expr_weak_ko, weak_data$meta)
      save_cellchat_results(cellchat_weak_ko, sprintf("weak_%s", ko_name))
    }
    
    ft_weak_ko <- dump_fib_treg(cellchat_weak_ko, "Weak KO")
    ft_weak_base <- results[["weak_base"]]$fib_treg
    
    base_prob <- if (nrow(ft_weak_base) > 0) sum(ft_weak_base$prob) else 0
    ko_prob   <- if (nrow(ft_weak_ko) > 0) sum(ft_weak_ko$prob) else 0
    retained  <- if (base_prob > 0) round(ko_prob / base_prob * 100, 2) else NA
    lost      <- if (base_prob > 0) round(100 - retained, 2) else NA
    
    ko_summary <- rbind(ko_summary, data.frame(
      KO = ko_name,
      Group = "Weak",
      Base_Interactions = nrow(ft_weak_base),
      KO_Interactions = nrow(ft_weak_ko),
      Base_Prob = base_prob,
      KO_Prob = ko_prob,
      Retained_Pct = retained,
      Lost_Pct = lost,
      stringsAsFactors = FALSE
    ))
  }
}

# --------------------------------------------------------
# 4. 汇总结果
# --------------------------------------------------------
cat("\n", rep("=", 60), "\n", sep = "")
cat("虚拟敲除汇总结果\n")
cat(rep("=", 60), "\n", sep = "")
print(ko_summary, row.names = FALSE)

# 保存汇总
summary_file <- file.path(OUTPUT_DIR, "ko_summary.csv")
write.csv(ko_summary, summary_file, row.names = FALSE)
cat(sprintf("\n汇总已保存: %s\n", summary_file))

# --------------------------------------------------------
# 5. 生成对比可视化
# --------------------------------------------------------
cat("\n生成对比图...\n")

# 绘制 KO 后保留概率的条形图
if (nrow(ko_summary) > 0) {
  ko_summary$KO_Group <- paste(ko_summary$KO, ko_summary$Group, sep = "\n")
  
  p <- ggplot(ko_summary, aes(x = KO, y = Lost_Pct, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", Lost_Pct)), 
              position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("Strong" = "#2E86AB", "Weak" = "#A23B72")) +
    labs(
      title = "CellChat Virtual Knockout: Fibroblast -> Treg Signal Loss",
      subtitle = "KO of ECM ligands (COL1A1, COL1A2, FN1)",
      x = "Knockout Condition",
      y = "Signal Lost (%)",
      fill = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "top"
    ) +
    ylim(0, 105)
  
  plot_file <- file.path(OUTPUT_DIR, "ko_fib_treg_signal_loss.pdf")
  ggsave(plot_file, p, width = 10, height = 6)
  cat(sprintf("对比图已保存: %s\n", plot_file))
  
  # 同时保存为 PNG
  plot_file_png <- file.path(OUTPUT_DIR, "ko_fib_treg_signal_loss.png")
  ggsave(plot_file_png, p, width = 10, height = 6, dpi = 300)
  cat(sprintf("对比图已保存: %s\n", plot_file_png))
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("CellChat Virtual Knockout 分析完成!\n")
cat(sprintf("所有结果保存在: %s\n", OUTPUT_DIR))
cat(rep("=", 60), "\n", sep = "")
