#!/usr/bin/env python3
"""
Re-prepare GSE207422 CellChat input using CORRECT annotation (cell_metadata_full.csv)
instead of the flawed all_cells_annotation.csv
"""

import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix, vstack
from scipy.io import mmwrite
import os
import sys

# ==================== CONFIG ====================
# Find correct metadata file (handle encoding issues in path)
import glob
META_FILE = None
for p in glob.glob(r"D:\Research\*\cell_metadata_full.csv") + glob.glob(r"D:\Research\*\*\cell_metadata_full.csv"):
    if "Figure5" not in p:
        META_FILE = p
        break

if not META_FILE:
    raise FileNotFoundError("cell_metadata_full.csv not found")

print(f"Using metadata: {META_FILE}")
EXPR_FILE = r"D:\Research\potato\data\raw\GSE207422\GSE207422_NSCLC_scRNAseq_UMI_matrix.txt"
OUTPUT_DIR = r"D:\Research\tomato\data\gse207422_cellchat_input_corrected"
# ================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)


def build_sparse_matrix(expr_file, cells_subset, meta_subset, output_prefix):
    """Read expression matrix in chunks and build sparse matrix for subset of cells"""
    chunk_size = 1000
    sparse_chunks = []
    gene_names = []
    gene_counter = 0

    n_chunks = 24292 // chunk_size + 1
    print(f"\n[{output_prefix}] Building sparse matrix, ~{n_chunks} chunks...")

    for chunk_idx, chunk in enumerate(
        pd.read_csv(expr_file, sep="\t", usecols=["Gene"] + cells_subset, chunksize=chunk_size)
    ):
        genes_chunk = chunk["Gene"].values.tolist()
        gene_names.extend(genes_chunk)

        chunk_vals = chunk.drop("Gene", axis=1).values.astype(np.int32)
        sparse_chunk = csr_matrix(chunk_vals)
        sparse_chunks.append(sparse_chunk)

        gene_counter += len(genes_chunk)
        if (chunk_idx + 1) % 5 == 0 or chunk_idx == 0:
            print(f"  Chunk {chunk_idx + 1}/{n_chunks}, {gene_counter} genes processed")

    print(f"[{output_prefix}] Merging {len(sparse_chunks)} sparse chunks...")
    full_matrix = vstack(sparse_chunks)

    mtx_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_matrix.mtx")
    print(f"[{output_prefix}] Saving MTX: {mtx_path}")
    mmwrite(mtx_path, full_matrix)

    genes_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_genes.tsv")
    with open(genes_path, "w", encoding="utf-8") as f:
        for g in gene_names:
            f.write(g + "\n")

    barcodes_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_barcodes.tsv")
    with open(barcodes_path, "w", encoding="utf-8") as f:
        for c in cells_subset:
            f.write(c + "\n")

    meta_reordered = meta_subset.set_index("cell").loc[cells_subset].reset_index()
    meta_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_metadata.csv")
    meta_reordered.to_csv(meta_path, index=False)

    print(f"[{output_prefix}] Done! Matrix: {full_matrix.shape}, nnz: {full_matrix.nnz}")
    return full_matrix


def main():
    print("=" * 60)
    print("GSE207422 CellChat Data Re-preparation (CORRECTED)")
    print("=" * 60)

    # 1. Read correct metadata
    print("\n[1/4] Reading correct metadata...")
    meta = pd.read_csv(META_FILE)
    meta = meta.rename(columns={"Unnamed: 0": "cell"})
    print(f"Total cells in metadata: {len(meta)}")
    print(f"Response groups: {dict(meta['response_group'].value_counts())}")

    # Keep only strong/weak
    meta = meta[meta["response_group"].isin(["strong", "weak"])].copy()
    print(f"After filtering strong/weak: {len(meta)}")

    # 2. Read expression matrix header
    print("\n[2/4] Reading expression matrix header...")
    with open(EXPR_FILE, "r") as f:
        header = f.readline().strip().split("\t")
    all_cells = header[1:]
    print(f"Total cells in matrix: {len(all_cells)}")

    # 3. Match cells
    print("\n[3/4] Matching cells...")
    meta_cells = set(meta["cell"])
    matrix_cells = set(all_cells)
    common_cells = sorted(meta_cells & matrix_cells)
    print(f"Common cells: {len(common_cells)}")

    meta_strong = meta[meta["response_group"] == "strong"][["cell", "cell_type", "response_group"]].copy()
    meta_weak = meta[meta["response_group"] == "weak"][["cell", "cell_type", "response_group"]].copy()

    strong_cells = [c for c in common_cells if c in meta_strong["cell"].values]
    weak_cells = [c for c in common_cells if c in meta_weak["cell"].values]

    print(f"Strong cells: {len(strong_cells)}")
    print(f"Weak cells: {len(weak_cells)}")

    print("\nStrong cell type distribution:")
    print(meta_strong[meta_strong["cell"].isin(strong_cells)]["cell_type"].value_counts().to_string())
    print("\nWeak cell type distribution:")
    print(meta_weak[meta_weak["cell"].isin(weak_cells)]["cell_type"].value_counts().to_string())

    # 4. Build and save
    print("\n[4/4] Building sparse matrices...")
    build_sparse_matrix(EXPR_FILE, strong_cells, meta_strong, "strong")
    build_sparse_matrix(EXPR_FILE, weak_cells, meta_weak, "weak")

    print(f"\n{'=' * 60}")
    print(f"All done! Output: {OUTPUT_DIR}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
