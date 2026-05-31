"""
Extract full CCR8+ and MKI67+ Treg from GSE243013 for scTenifoldKnk KO.
Output: raw count matrix (genes x cells) in .mtx format.
"""
import os, sys, gzip, pandas as pd, numpy as np
from scipy import sparse
from scipy.io import mmwrite

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)

MTX_PATH = r"D:\Research\土豆\数据\raw\GSE243013\GSE243013_NSCLC_immune_scRNA_counts.mtx.gz"
BARCODE_PATH = r"D:\Research\土豆\数据\raw\GSE243013\GSE243013_barcodes.csv.gz"
GENE_PATH = r"D:\Research\土豆\数据\raw\GSE243013\features.csv.gz"
META_PATH = r"D:\Research\土豆\数据\raw\GSE243013\GSE243013_NSCLC_immune_scRNA_metadata.csv"

OUT_BASE = r"D:\Research\tomato\data"

# Load metadata
print("Loading metadata...", flush=True)
meta = pd.read_csv(META_PATH, low_memory=False)
print(f"Total cells: {len(meta)}", flush=True)
print(f"Sub types: {meta['sub_cell_type'].value_counts().head(10).to_dict()}", flush=True)

# Filter Treg subtypes
ccr8_cells = meta[meta['sub_cell_type'] == 'CD4T_Treg_CCR8']['cellID'].tolist()
mki67_cells = meta[meta['sub_cell_type'] == 'CD4T_Treg_MKI67']['cellID'].tolist()

print(f"\nCCR8+ Treg: {len(ccr8_cells)} cells", flush=True)
print(f"MKI67+ Treg: {len(mki67_cells)} cells", flush=True)

# Load barcodes/genes
with gzip.open(BARCODE_PATH, "rt") as f:
    barcodes = [line.strip() for line in f]
with gzip.open(GENE_PATH, "rt") as f:
    gene_names = [line.strip() for line in f]

barcode_to_rowidx = {b: i+1 for i, b in enumerate(barcodes)}

def extract_and_save(target_barcodes, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    target_rowidx = {barcode_to_rowidx[bc]: bc for bc in target_barcodes if bc in barcode_to_rowidx}
    target_rows_set = set(target_rowidx.keys())
    print(f"\nExtracting to {out_dir}: {len(target_rows_set)} cells", flush=True)
    
    row_indices, col_indices, data_values = [], [], []
    chunk_n = 0
    for chunk in pd.read_csv(MTX_PATH, skiprows=3, header=None, sep=r'\s+',
                             dtype={0: 'int32', 1: 'int32', 2: 'float32'},
                             chunksize=300_000_000, engine='c', compression='gzip'):
        chunk_n += 1
        cell_idx = chunk[0].values
        gene_idx = chunk[1].values
        vals = chunk[2].values
        mask = np.isin(cell_idx, list(target_rows_set))
        if mask.sum() > 0:
            row_indices.extend(cell_idx[mask].tolist())
            col_indices.extend(gene_idx[mask].tolist())
            data_values.extend(vals[mask].tolist())
        print(f"  Chunk {chunk_n}: kept {len(data_values)} total", flush=True)
    
    # Build matrix (cells x genes), then transpose to genes x cells for R
    row_map = {old: new for new, old in enumerate(sorted(target_rowidx.keys()))}
    cell_idx_to_barcode = {v: target_rowidx[k] for k, v in row_map.items()}
    
    new_rows = [row_map[r] for r in row_indices]
    new_cols = [c - 1 for c in col_indices]
    n_cells = len(target_rowidx)
    n_genes = len(gene_names)
    
    mat = sparse.csr_matrix((data_values, (new_rows, new_cols)), shape=(n_cells, n_genes))
    mat_t = mat.T.tocsr()  # genes x cells
    
    print(f"Matrix shape: {mat_t.shape}", flush=True)
    
    mmwrite(os.path.join(out_dir, "matrix.mtx"), mat_t)
    with open(os.path.join(out_dir, "barcodes.tsv"), "w") as f:
        for i in range(n_cells):
            f.write(cell_idx_to_barcode[i] + "\n")
    with open(os.path.join(out_dir, "genes.tsv"), "w") as f:
        for g in gene_names:
            f.write(g + "\n")
    
    print(f"Saved to {out_dir}", flush=True)

extract_and_save(ccr8_cells, os.path.join(OUT_BASE, "treg_ccr8_for_ko"))
extract_and_save(mki67_cells, os.path.join(OUT_BASE, "treg_mki67_for_ko"))
print("\nAll done!", flush=True)
