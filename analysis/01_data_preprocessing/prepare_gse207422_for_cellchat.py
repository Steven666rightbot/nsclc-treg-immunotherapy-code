#!/usr/bin/env python3
"""
Prepare GSE207422 data for CellChat virtual knockout analysis
Reads the 4.49GB UMI matrix, splits by strong/weak response, saves as sparse MTX
"""

import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix, vstack
from scipy.io import mmwrite
import os
import sys

# ==================== CONFIG ====================
BASE_DIR = r"D:\Research\potato\data\raw\GSE207422"
OUTPUT_DIR = r"D:\Research\tomato\data\gse207422_cellchat_input"
EXPR_FILE = os.path.join(BASE_DIR, "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt")
META_SAMPLE_FILE = os.path.join(BASE_DIR, "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
META_CELL_FILE = os.path.join(BASE_DIR, "annotation_results", "all_cells_annotation.csv")
# ================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)


def classify_response(pr):
    """MPR/pCR = strong, NMPR = weak, NE = unknown"""
    if pd.isna(pr):
        return "unknown"
    pr = str(pr).strip()
    if pr in ["MPR", "pCR"]:
        return "strong"
    elif pr == "NMPR":
        return "weak"
    else:
        return "unknown"


def build_sparse_matrix(expr_file, cells_subset, meta_subset, output_prefix):
    """Read expression matrix in chunks and build sparse matrix for subset of cells"""
    chunk_size = 1000
    sparse_chunks = []
    gene_names = []
    gene_counter = 0

    n_chunks = 24292 // chunk_size + 1
    print(f"\n[{output_prefix}] 构建稀疏矩阵，约 {n_chunks} 个块...")

    for chunk_idx, chunk in enumerate(
        pd.read_csv(expr_file, sep="\t", usecols=["Gene"] + cells_subset, chunksize=chunk_size)
    ):
        genes_chunk = chunk["Gene"].values.tolist()
        gene_names.extend(genes_chunk)

        # Drop Gene column, convert to numpy int32
        chunk_vals = chunk.drop("Gene", axis=1).values.astype(np.int32)

        # Build sparse chunk
        sparse_chunk = csr_matrix(chunk_vals)
        sparse_chunks.append(sparse_chunk)

        gene_counter += len(genes_chunk)
        if (chunk_idx + 1) % 5 == 0 or chunk_idx == 0:
            print(f"  块 {chunk_idx + 1}/{n_chunks}，已处理 {gene_counter} 个基因")

    # Merge all chunks
    print(f"[{output_prefix}] 合并 {len(sparse_chunks)} 个稀疏块...")
    full_matrix = vstack(sparse_chunks)

    # Save as MTX (Matrix Market format)
    mtx_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_matrix.mtx")
    print(f"[{output_prefix}] 保存 MTX: {mtx_path}")
    mmwrite(mtx_path, full_matrix)

    # Save gene names
    genes_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_genes.tsv")
    with open(genes_path, "w", encoding="utf-8") as f:
        for g in gene_names:
            f.write(g + "\n")

    # Save cell barcodes
    barcodes_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_barcodes.tsv")
    with open(barcodes_path, "w", encoding="utf-8") as f:
        for c in cells_subset:
            f.write(c + "\n")

    # Save metadata (reordered to match matrix columns)
    meta_reordered = meta_subset.set_index("cell").loc[cells_subset].reset_index()
    meta_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_metadata.csv")
    meta_reordered.to_csv(meta_path, index=False)

    print(f"[{output_prefix}] 完成! 矩阵: {full_matrix.shape}, 非零元素: {full_matrix.nnz}")
    return full_matrix


def main():
    print("=" * 60)
    print("GSE207422 CellChat 数据预处理")
    print("=" * 60)

    # 1. Read sample metadata
    print("\n[1/5] 读取样本 metadata...")
    meta_sample = pd.read_excel(META_SAMPLE_FILE)
    meta_sample = meta_sample.dropna(subset=["Sample"])
    meta_sample["response_group"] = meta_sample["Pathologic Response"].apply(classify_response)
    print(meta_sample[["Sample", "Pathologic Response", "response_group"]].to_string(index=False))

    # 2. Read cell annotation
    print("\n[2/5] 读取细胞注释...")
    meta_cell = pd.read_csv(META_CELL_FILE)
    meta_cell["sample"] = meta_cell["cell"].apply(lambda x: "_".join(x.split("_")[:2]))
    meta_cell = meta_cell.merge(
        meta_sample[["Sample", "response_group"]],
        left_on="sample",
        right_on="Sample",
        how="left",
    )

    # Cell type mapping to match previous CellChat analysis
    cell_type_map = {
        "Regulatory T cells": "Treg",
        "Treg(diff)": "Treg",
        "Fibroblasts": "Fibroblast",
        "Epithelial cells": "Epithelial",
        "Endothelial cells": "Endothelial",
        "Alveolar macrophages": "Myeloid",
        "Classical monocytes": "Myeloid",
        "Macrophages": "Myeloid",
        "Intermediate macrophages": "Myeloid",
        "Intestinal macrophages": "Myeloid",
        "Monocytes": "Myeloid",
        "Non-classical monocytes": "Myeloid",
        "Mono-mac": "Myeloid",
        "Kupffer cells": "Myeloid",
        "Erythrophagocytic macrophages": "Myeloid",
        "Migratory DCs": "Myeloid",
        "DC2": "DC",
        "DC": "DC",
        "DC1": "DC",
        "DC3": "DC",
        "DC precursor": "DC",
        "Transitional DC": "DC",
        "Cycling DCs": "DC",
        "pDC": "DC",
        "Tem/Trm cytotoxic T cells": "CD8_T",
        "Tem/Temra cytotoxic T cells": "CD8_T",
        "Trm cytotoxic T cells": "CD8_T",
        "Tcm/Naive cytotoxic T cells": "CD8_T",
        "Cycling T cells": "T_cells",
        "Tem/Effector helper T cells": "CD4_T",
        "Tcm/Naive helper T cells": "CD4_T",
        "Type 1 helper T cells": "CD4_T",
        "Type 17 helper T cells": "CD4_T",
        "Follicular helper T cells": "CD4_T",
        "Tem/Effector helper T cells PD1+": "CD4_T",
        "NK cells": "NK",
        "CD16+ NK cells": "NK",
        "CD16- NK cells": "NK",
        "Transitional NK": "NK",
        "Cycling NK cells": "NK",
        "Memory B cells": "B_cells",
        "Naive B cells": "B_cells",
        "Age-associated B cells": "B_cells",
        "Follicular B cells": "B_cells",
        "Transitional B cells": "B_cells",
        "Germinal center B cells": "B_cells",
        "Proliferative germinal center B cells": "B_cells",
        "B cells": "B_cells",
        "Cycling B cells": "B_cells",
        "Plasma cells": "Plasma",
        "Plasmablasts": "Plasma",
        "Mast cells": "Mast",
        "CRTAM+ gamma-delta T cells": "T_cells",
        "gamma-delta T cells": "T_cells",
        "MAIT cells": "T_cells",
        "NKT cells": "T_cells",
        "T(agonist)": "T_cells",
        "CD8a/b(entry)": "CD8_T",
        "CD8a/a": "CD8_T",
        "Memory CD4+ cytotoxic T cells": "CD4_T",
    }

    meta_cell["cell_type"] = meta_cell["predicted_labels"].map(cell_type_map).fillna("Other")

    # 3. Filter strong and weak
    print("\n[3/5] 筛选 strong / weak 细胞...")
    meta_strong = meta_cell[meta_cell["response_group"] == "strong"][
        ["cell", "cell_type", "response_group"]
    ].copy()
    meta_weak = meta_cell[meta_cell["response_group"] == "weak"][
        ["cell", "cell_type", "response_group"]
    ].copy()

    print(f"Strong cells: {len(meta_strong)}")
    print(f"Weak cells: {len(meta_weak)}")
    print("\nStrong 细胞类型分布:")
    print(meta_strong["cell_type"].value_counts().to_string())
    print("\nWeak 细胞类型分布:")
    print(meta_weak["cell_type"].value_counts().to_string())

    # 4. Read expression matrix header
    print("\n[4/5] 读取表达矩阵列名...")
    df_header = pd.read_csv(EXPR_FILE, sep="\t", nrows=0)
    all_cells = df_header.columns.tolist()[1:]  # Exclude 'Gene'
    print(f"表达矩阵总细胞数: {len(all_cells)}")

    strong_cells = [c for c in all_cells if c in set(meta_strong["cell"])]
    weak_cells = [c for c in all_cells if c in set(meta_weak["cell"])]

    print(f"Strong cells in matrix: {len(strong_cells)}")
    print(f"Weak cells in matrix: {len(weak_cells)}")

    # 5. Build and save sparse matrices
    print("\n[5/5] 构建稀疏矩阵并保存...")
    build_sparse_matrix(EXPR_FILE, strong_cells, meta_strong, "strong")
    build_sparse_matrix(EXPR_FILE, weak_cells, meta_weak, "weak")

    print(f"\n{'=' * 60}")
    print(f"全部完成! 输出目录: {OUTPUT_DIR}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
