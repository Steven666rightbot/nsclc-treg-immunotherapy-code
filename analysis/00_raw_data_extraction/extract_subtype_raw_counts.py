#!/usr/bin/env python3
"""
Extract raw counts for CCR8+ and MKI67+ Treg from GSE243013 original MTX.
No normalization, no log1p — pure integer counts for scTenifoldKnk.
Strategy: top 2000 HVG + KO targets + key genes, export as CSV.
"""

import os
import time
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc

DATA_DIR = r"D:\Research\土豆\数据\raw\GSE243013"
OUT_DIR = r"D:\Research\tomato\vB_revised_integrin_fyn_axis\gse243013_ko_input_raw"
os.makedirs(OUT_DIR, exist_ok=True)

# Load Treg metadata to get barcodes + subtypes
print("Loading Treg metadata...")
treg_meta = pd.read_csv(os.path.join(DATA_DIR, "GSE243013_Treg_cells_metadata.csv"), low_memory=False)
treg_meta = treg_meta[treg_meta['pathological_response'] != 'unknowm'].copy()

# Target subtypes
TARGET_SUBTYPES = ['CD4T_Treg_CCR8', 'CD4T_Treg_MKI67']
subtype_masks = {}
for st in TARGET_SUBTYPES:
    mask = treg_meta['sub_cell_type'] == st
    barcodes = treg_meta.loc[mask, 'cellID'].tolist()
    subtype_masks[st] = set(barcodes)
    print(f"  {st}: {len(barcodes)} cells")

# Combine all target barcodes
all_target_barcodes = subtype_masks['CD4T_Treg_CCR8'] | subtype_masks['CD4T_Treg_MKI67']
print(f"  Total target cells: {len(all_target_barcodes)}")

# Read barcodes and features
print("\nReading barcodes and features...")
barcodes_df = pd.read_csv(os.path.join(DATA_DIR, "GSE243013_barcodes.csv.gz"), compression='gzip', header=None, names=['barcode'])
barcodes_list = barcodes_df['barcode'].tolist()
print(f"  Total barcodes: {len(barcodes_list)}")

features_df = pd.read_csv(os.path.join(DATA_DIR, "features.csv.gz"), compression='gzip', header=None, names=['geneSymbol'])
gene_names = features_df['geneSymbol'].tolist()
print(f"  Total genes: {len(gene_names)}")

# Build barcode -> index map (MTX uses 1-based indexing)
barcode_to_idx = {bc: i for i, bc in enumerate(barcodes_list[1:], start=1)}  # 1-based

# Build target cell index sets
target_indices = {st: set() for st in TARGET_SUBTYPES}
for st in TARGET_SUBTYPES:
    for bc in subtype_masks[st]:
        if bc in barcode_to_idx:
            target_indices[st].add(barcode_to_idx[bc])
    print(f"  {st} indices mapped: {len(target_indices[st])}")

# Also build a combined index set for streaming extraction
all_target_indices = target_indices['CD4T_Treg_CCR8'] | target_indices['CD4T_Treg_MKI67']

# Stream-read MTX and extract target cells
mtx_path = os.path.join(DATA_DIR, "GSE243013_NSCLC_immune_scRNA_counts.mtx")
print(f"\nStreaming MTX: {mtx_path}")
print("This will take ~15-20 minutes...")

start_time = time.time()

# We need to collect (row, col, data) for all target cells
# row = cell index (1-based in MTX), col = gene index (1-based)
all_rows = []
all_cols = []
all_data = []
total_read = 0
chunk_num = 0

# MTX has 3 header lines, then data
for chunk in pd.read_csv(mtx_path, skiprows=3, header=None, sep=r'\s+',
                         dtype={0: 'int32', 1: 'int32', 2: 'float32'},
                         chunksize=100_000_000, engine='c'):
    chunk_num += 1
    total_read += len(chunk)
    
    # MTX format: col gene_idx, row cell_idx, value
    # Wait, let's check. Standard Matrix Market: row_idx col_idx value
    # From the export script: chunk[0] is row (cell), chunk[1] is col (gene)
    cell_idx = chunk[0].values
    gene_idx = chunk[1].values
    vals = chunk[2].values
    
    mask = np.isin(cell_idx, list(all_target_indices))
    n_kept = mask.sum()
    
    if n_kept > 0:
        all_rows.append(cell_idx[mask])
        all_cols.append(gene_idx[mask])
        all_data.append(vals[mask])
    
    if chunk_num % 2 == 0:
        elapsed = time.time() - start_time
        print(f"  Chunk {chunk_num:02d}: read {total_read:>12,} | kept {n_kept:>8,} | elapsed: {elapsed/60:5.1f}min")

print(f"\nTotal chunks: {chunk_num}")
print(f"Total entries read: {total_read}")

# Merge
print("Merging chunks...")
all_rows = np.concatenate(all_rows)
all_cols = np.concatenate(all_cols)
all_data = np.concatenate(all_data)
print(f"Total kept entries: {len(all_data)}")

# Build sparse matrix: cells x genes
# Need to remap cell indices to 0-based contiguous indices for the subset
cell_idx_map = {}
current_idx = 0
for st in TARGET_SUBTYPES:
    for idx in sorted(target_indices[st]):
        if idx not in cell_idx_map:
            cell_idx_map[idx] = current_idx
            current_idx += 1

remapped_rows = np.array([cell_idx_map[idx] for idx in all_rows])
remapped_cols = all_cols - 1  # MTX is 1-based, matrix is 0-based

n_cells = len(cell_idx_map)
n_genes = len(gene_names)

X = sparse.coo_matrix((all_data, (remapped_rows, remapped_cols)), shape=(n_cells, n_genes))
X = X.tocsr()
print(f"CSR matrix: {X.shape}")

# Build cell labels
idx_to_cell = {v: k for k, v in cell_idx_map.items()}
cell_labels = [barcodes_list[idx_to_cell[i]] for i in range(n_cells)]
cell_subtypes = []
for bc in cell_labels:
    for st in TARGET_SUBTYPES:
        if bc in subtype_masks[st]:
            cell_subtypes.append(st)
            break

print(f"Cell label check: {len(cell_labels)} labels, {len(cell_subtypes)} subtypes")

# Create AnnData with raw counts
adata = sc.AnnData(X)
adata.obs_names = cell_labels
adata.var_names = gene_names
adata.obs['sub_cell_type'] = cell_subtypes

print(f"Raw adata: {adata.shape}")
print(f"Value check (should be integers): min={adata.X.min()}, max={adata.X.max()}")

# Now split by subtype and compute HVG for each
TARGET_RECEPTORS = ["CD44", "ITGB1", "DDR1", "DDR2", "ITGA4"]
TARGET_TFS = ["MEF2C", "TEAD1", "YAP1", "WWTR1"]
TARGETS = TARGET_RECEPTORS + TARGET_TFS

TREG_FUNCTION = ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "ICOS", "IL2RB"]
INTEGRIN_DOWNSTREAM = ["RHOA", "ROCK1", "ROCK2", "FYN", "LCK", "PIK3CD", "AKT1", "MAPK1", "PTPN11"]
SURVIVAL_APOPTOSIS = ["MCL1", "BCL2", "BAX", "BCL2L1", "PTEN", "CASP3"]
METABOLISM = ["LDHA", "PFKFB3"]
TREG_MARKERS = ["FOXP3", "IL2RA", "CTLA4"]
KEY_GENES = list(set(TREG_FUNCTION + INTEGRIN_DOWNSTREAM + SURVIVAL_APOPTOSIS + METABOLISM + TREG_MARKERS))

for st in TARGET_SUBTYPES:
    st_short = st.replace('CD4T_Treg_', '')
    print(f"\n{'='*60}")
    print(f"Processing: {st} -> {st_short}")
    print(f"{'='*60}")
    
    sub = adata[adata.obs['sub_cell_type'] == st].copy()
    print(f"  Cells: {sub.n_obs}")
    
    # Compute HVG on raw counts (using variance)
    expr = sub.X.toarray() if hasattr(sub.X, 'toarray') else sub.X
    gene_vars = np.var(expr, axis=0)
    gene_names_arr = np.array(sub.var_names)
    
    top2000_idx = np.argsort(gene_vars)[::-1][:2000]
    top2000_genes = set(gene_names_arr[top2000_idx])
    
    # Force-keep targets and key genes
    all_keep = top2000_genes | set(TARGETS) | set(KEY_GENES)
    all_keep = [g for g in gene_names_arr if g in all_keep]
    
    print(f"  Top 2000 HVG: {len(top2000_genes)}")
    print(f"  Targets added: {sum(1 for t in TARGETS if t not in top2000_genes)}")
    print(f"  Key genes added: {sum(1 for k in KEY_GENES if k not in top2000_genes)}")
    print(f"  Final gene set: {len(all_keep)}")
    
    sub = sub[:, all_keep].copy()
    
    # Export as dense CSV for R (genes as rows, cells as columns)
    print("  Exporting dense CSV...")
    mat = sub.X.toarray() if hasattr(sub.X, 'toarray') else sub.X
    
    # Ensure integer type for raw counts
    if not np.issubdtype(mat.dtype, np.integer):
        mat = np.round(mat).astype(np.int32)
    
    mat_df = pd.DataFrame(mat.T, index=sub.var_names, columns=sub.obs_names)
    
    out_csv = os.path.join(OUT_DIR, f'treg_{st_short}_gse243013_raw.csv')
    mat_df.to_csv(out_csv)
    print(f"  Saved: {out_csv} ({mat_df.shape[0]} genes x {mat_df.shape[1]} cells)")
    
    # Save metadata
    meta = sub.obs[['sub_cell_type']].copy()
    meta_out = os.path.join(OUT_DIR, f'meta_{st_short}_gse243013_raw.csv')
    meta.to_csv(meta_out)
    print(f"  Saved metadata: {meta_out}")
    
    # Verify expression of KO targets
    for t in TARGETS:
        if t in sub.var_names:
            mean_expr = sub[:, t].X.mean()
            print(f"    {t} mean: {mean_expr:.4f}")

elapsed = (time.time() - start_time) / 60
print(f"\n[OK] Total time: {elapsed:.1f} minutes")
print(f"Output: {OUT_DIR}")
