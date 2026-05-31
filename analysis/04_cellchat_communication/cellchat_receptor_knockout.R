# =============================================================================
# CellChat Receptor Virtual Knockout Analysis - GSE207422 (CORRECTED annotation)
# 
# Knockout ECM receptors specifically in Treg cells, then rebuild CellChat
# to measure loss of Fibroblast->Treg communication probability.
# =============================================================================

library(CellChat)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(Matrix)
library(data.table)

# ======================== 参数配置 ========================
INPUT_DIR  <- "D:/Research/tomato/data/gse207422_cellchat_input_corrected"
OUTPUT_DIR <- "D:/Research/tomato/results/cellchat_ko_corrected"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

RUN_STRONG <- TRUE
RUN_WEAK   <- TRUE

# Receptor KO list: target genes to knockout in Treg cells
KO_LIST <- list(
  "KO_REC_CD44"         = c("CD44"),
  "KO_REC_ITGB1"        = c("ITGB1"),
  "KO_REC_ITGA4"        = c("ITGA4"),
  "KO_REC_DDR1"         = c("DDR1"),
  "KO_REC_DDR2"         = c("DDR2"),
  "KO_REC_CD44_ITGB1"   = c("CD44", "ITGB1"),
  "KO_REC_ALL_ECM_REC"  = c("CD44", "ITGB1", "ITGA4", "DDR1", "DDR2", "ITGA1", "ITGA5")
)

SKIP_EXISTING_RDS <- TRUE  # Skip existing RDS files

# ======================== 辅助函数 ========================

read_cellchat_input <- function(prefix) {
  cat(sprintf("\n[%s] 读取输入数据...\n", prefix))
  
  mtx_file <- file.path(INPUT_DIR, sprintf("%s_matrix.mtx", prefix))
  expr_matrix <- readMM(mtx_file)
  
  genes_file <- file.path(INPUT_DIR, sprintf("%s_genes.tsv", prefix))
  genes <- fread(genes_file, header = FALSE)$V1
  
  barcodes_file <- file.path(INPUT_DIR, sprintf("%s_barcodes.tsv", prefix))
  barcodes <- fread(barcodes_file, header = FALSE)$V1
  
  meta_file <- file.path(INPUT_DIR, sprintf("%s_metadata.csv", prefix))
  meta <- fread(meta_file)
  
  rownames(expr_matrix) <- genes
  colnames(expr_matrix) <- barcodes
  expr_matrix <- as(expr_matrix, "dgCMatrix")
  
  cat(sprintf("[%s] 矩阵: %d genes x %d cells\n", prefix, nrow(expr_matrix), ncol(expr_matrix)))
  return(list(expr = expr_matrix, meta = meta, barcodes = barcodes, genes = genes))
}

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

# Virtual KO that zeros out genes ONLY in specific cell types
virtual_ko_celltype <- function(expr_matrix, ko_genes, meta, target_celltypes = c("Treg")) {
  expr_ko <- expr_matrix
  
  # Get barcodes of target cells
  target_barcodes <- meta$barcode[meta$cell_type %in% target_celltypes]
  if (length(target_barcodes) == 0) {
    # Try alternative column name
    target_barcodes <- meta[[1]][meta$cell_type %in% target_celltypes]
  }
  
  target_idx <- which(colnames(expr_ko) %in% target_barcodes)
  cat(sprintf("    目标细胞类型 %s: %d 个细胞将被KO\n", 
              paste(target_celltypes, collapse = "/"), length(target_idx)))
  
  for (g in ko_genes) {
    if (g %in% rownames(expr_ko)) {
      if (length(target_idx) > 0) {
        expr_ko[g, target_idx] <- 0
        cat(sprintf("    KO: %s (仅在目标细胞类型中设为0, 影响 %d 个细胞)\n", g, length(target_idx)))
      } else {
        cat(sprintf("    警告: 未找到目标细胞类型，跳过 %s\n", g))
      }
    } else {
      cat(sprintf("    警告: %s 不在表达矩阵中，跳过\n", g))
    }
  }
  return(expr_ko)
}

dump_fib_treg <- function(cellchat, label) {
  df <- subsetCommunication(cellchat)
  ft <- df[df$source == "Fibroblast" & df$target == "Treg", ]
  cat(sprintf("    %s Fibroblast->Treg 互作: %d 条\n", label, nrow(ft)))
  if (nrow(ft) > 0) {
    cat(sprintf("    %s 总 prob: %.6e\n", label, sum(ft$prob)))
  }
  return(ft)
}

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
cat("CellChat Receptor Virtual Knockout Analysis (CORRECTED)\n")
cat("Target: Treg cells only\n")
cat("=" , rep("=", 58), "\n", sep = "")

# Read input data
if (RUN_STRONG) {
  strong_data <- read_cellchat_input("strong")
}
if (RUN_WEAK) {
  weak_data <- read_cellchat_input("weak")
}

# Load base results (already computed)
results <- list()

if (RUN_STRONG) {
  cat("\n--- Strong 响应者 基础分析 (读取已有RDS) ---\n")
  rds_strong <- file.path(OUTPUT_DIR, "cellchat_strong.rds")
  if (file.exists(rds_strong)) {
    cellchat_strong <- readRDS(rds_strong)
    cat("RDS 已加载\n")
  } else {
    stop("Strong baseline RDS not found!")
  }
  results[["strong_base"]] <- list(
    obj = cellchat_strong,
    fib_treg = dump_fib_treg(cellchat_strong, "Strong")
  )
}

if (RUN_WEAK) {
  cat("\n--- Weak 响应者 基础分析 (读取已有RDS) ---\n")
  rds_weak <- file.path(OUTPUT_DIR, "cellchat_weak.rds")
  if (file.exists(rds_weak)) {
    cellchat_weak <- readRDS(rds_weak)
    cat("RDS 已加载\n")
  } else {
    stop("Weak baseline RDS not found!")
  }
  results[["weak_base"]] <- list(
    obj = cellchat_weak,
    fib_treg = dump_fib_treg(cellchat_weak, "Weak")
  )
}

# Receptor KO analysis
cat("\n", rep("-", 60), "\n", sep = "")
cat("开始受体虚拟敲除分析 (仅在Treg中敲除)\n")
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
  Target_Cells = integer(),
  stringsAsFactors = FALSE
)

for (ko_name in names(KO_LIST)) {
  ko_genes <- KO_LIST[[ko_name]]
  cat(sprintf("\n========== %s: %s ==========\n", ko_name, paste(ko_genes, collapse = " + ")))
  
  if (RUN_STRONG) {
    rds_ko_strong <- file.path(OUTPUT_DIR, sprintf("cellchat_strong_%s.rds", ko_name))
    
    if (SKIP_EXISTING_RDS && file.exists(rds_ko_strong)) {
      cat("[Strong] RDS 已存在，跳过...\n")
      cellchat_strong_ko <- readRDS(rds_ko_strong)
    } else {
      cat("[Strong] 执行受体 KO (仅在Treg中)...\n")
      expr_strong_ko <- virtual_ko_celltype(strong_data$expr, ko_genes, strong_data$meta, target_celltypes = c("Treg"))
      cellchat_strong_ko <- run_cellchat_pipeline(expr_strong_ko, strong_data$meta)
      save_cellchat_results(cellchat_strong_ko, sprintf("strong_%s", ko_name))
    }
    
    ft_strong_ko <- dump_fib_treg(cellchat_strong_ko, "Strong KO")
    ft_strong_base <- results[["strong_base"]]$fib_treg
    
    base_prob <- if (nrow(ft_strong_base) > 0) sum(ft_strong_base$prob) else 0
    ko_prob   <- if (nrow(ft_strong_ko) > 0) sum(ft_strong_ko$prob) else 0
    retained  <- if (base_prob > 0) round(ko_prob / base_prob * 100, 2) else NA
    lost      <- if (base_prob > 0) round(100 - retained, 2) else NA
    
    treg_count <- sum(strong_data$meta$cell_type == "Treg")
    
    ko_summary <- rbind(ko_summary, data.frame(
      KO = ko_name, Group = "Strong",
      Base_Interactions = nrow(ft_strong_base), KO_Interactions = nrow(ft_strong_ko),
      Base_Prob = base_prob, KO_Prob = ko_prob,
      Retained_Pct = retained, Lost_Pct = lost,
      Target_Cells = treg_count,
      stringsAsFactors = FALSE
    ))
  }
  
  if (RUN_WEAK) {
    rds_ko_weak <- file.path(OUTPUT_DIR, sprintf("cellchat_weak_%s.rds", ko_name))
    
    if (SKIP_EXISTING_RDS && file.exists(rds_ko_weak)) {
      cat("[Weak] RDS 已存在，跳过...\n")
      cellchat_weak_ko <- readRDS(rds_ko_weak)
    } else {
      cat("[Weak] 执行受体 KO (仅在Treg中)...\n")
      expr_weak_ko <- virtual_ko_celltype(weak_data$expr, ko_genes, weak_data$meta, target_celltypes = c("Treg"))
      cellchat_weak_ko <- run_cellchat_pipeline(expr_weak_ko, weak_data$meta)
      save_cellchat_results(cellchat_weak_ko, sprintf("weak_%s", ko_name))
    }
    
    ft_weak_ko <- dump_fib_treg(cellchat_weak_ko, "Weak KO")
    ft_weak_base <- results[["weak_base"]]$fib_treg
    
    base_prob <- if (nrow(ft_weak_base) > 0) sum(ft_weak_base$prob) else 0
    ko_prob   <- if (nrow(ft_weak_ko) > 0) sum(ft_weak_ko$prob) else 0
    retained  <- if (base_prob > 0) round(ko_prob / base_prob * 100, 2) else NA
    lost      <- if (base_prob > 0) round(100 - retained, 2) else NA
    
    treg_count <- sum(weak_data$meta$cell_type == "Treg")
    
    ko_summary <- rbind(ko_summary, data.frame(
      KO = ko_name, Group = "Weak",
      Base_Interactions = nrow(ft_weak_base), KO_Interactions = nrow(ft_weak_ko),
      Base_Prob = ko_prob, KO_Prob = ko_prob,
      Retained_Pct = retained, Lost_Pct = lost,
      Target_Cells = treg_count,
      stringsAsFactors = FALSE
    ))
  }
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("受体虚拟敲除汇总结果\n")
cat(rep("=", 60), "\n", sep = "")
print(ko_summary, row.names = FALSE)

# Append to existing ko_summary.csv or create new
summary_file <- file.path(OUTPUT_DIR, "ko_summary_receptor.csv")
write.csv(ko_summary, summary_file, row.names = FALSE)
cat(sprintf("\n汇总已保存: %s\n", summary_file))

cat("\n", rep("=", 60), "\n", sep = "")
cat("CellChat Receptor Virtual Knockout 分析完成!\n")
cat(sprintf("所有结果保存在: %s\n", OUTPUT_DIR))
cat(rep("=", 60), "\n", sep = "")
