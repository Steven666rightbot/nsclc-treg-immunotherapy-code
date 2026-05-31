"""
Prepare GSE243013 CellChat input - CORRECTED.
Original mtx is (cells, genes). Extract cells by row, transpose to (genes, cells).
"""
import os, gzip, random, json
from datetime import datetime
import pandas as pd

META_PATH = r"D:\Research\tomato\data\temp_umap_all\GSE243013_metadata.csv"
BARCODES_PATH = r"D:\Research\tomato\data\geo_cache\GSE243013_barcodes.csv.gz"
MTX_PATH = r"D:\Research\tomato\data\geo_cache\GSE243013_NSCLC_immune_scRNA_counts.mtx.gz"
GENES_PATH = r"D:\Research\tomato\data\geo_cache\GSE243013_genes.csv.gz"
OUT_DIR = r"D:\Research\tomato\data\gse243013_cellchat_input"

SEED = 20250517
N_PER_GROUP = 1000

os.makedirs(OUT_DIR, exist_ok=True)
random.seed(SEED)

# Step 1: Read metadata
print(f"[{datetime.now()}] Step 1: Reading metadata...")
chunks = []
for chunk in pd.read_csv(META_PATH, low_memory=False, chunksize=20000,
                          usecols=['cellID','sampleID','major_cell_type','sub_cell_type','pathological_response']):
    chunks.append(chunk)
meta = pd.concat(chunks)

meta_valid = meta[meta['pathological_response'].isin(['MPR', 'pCR', 'non-MPR'])].copy()
meta_valid['response_group'] = meta_valid['pathological_response'].apply(
    lambda x: 'strong' if x in ['MPR', 'pCR'] else 'weak')

sampled = []
for (ct, rg), group in meta_valid.groupby(['sub_cell_type', 'response_group']):
    n = min(N_PER_GROUP, len(group))
    sel = group.sample(n=n, random_state=SEED)
    sampled.append(sel)

sampled_df = pd.concat(sampled)
print(f"Sampled: {len(sampled_df)} | Strong: {sum(sampled_df['response_group']=='strong')} | Weak: {sum(sampled_df['response_group']=='weak')}")

strong_df = sampled_df[sampled_df['response_group']=='strong'][['cellID','sub_cell_type']].copy()
weak_df = sampled_df[sampled_df['response_group']=='weak'][['cellID','sub_cell_type']].copy()
strong_df.columns = ['cell', 'cell_type']
weak_df.columns = ['cell', 'cell_type']

# Step 2: Map barcodes to ROW indices (original mtx is cells x genes)
print(f"[{datetime.now()}] Step 2: Mapping barcodes to row indices...")
with gzip.open(BARCODES_PATH, 'rt') as f:
    all_barcodes = [l.strip() for l in f][1:]

barcode_to_row = {b: i+1 for i, b in enumerate(all_barcodes)}  # 1-based row index
strong_cells = set(strong_df['cell'].values)
weak_cells = set(weak_df['cell'].values)

strong_rows = sorted([barcode_to_row[c] for c in strong_cells if c in barcode_to_row])
weak_rows = sorted([barcode_to_row[c] for c in weak_cells if c in barcode_to_row])
strong_set = set(strong_rows)
weak_set = set(weak_rows)

# New col index mapping (consecutive 1-based for output mtx)
strong_col_map = {old: new for new, old in enumerate(sorted(strong_rows), 1)}
weak_col_map = {old: new for new, old in enumerate(sorted(weak_rows), 1)}

print(f"Strong rows: {len(strong_rows)} | Weak rows: {len(weak_rows)}")

# Step 3: Single pass - extract cells by ROW, transpose to (genes, cells)
print(f"[{datetime.now()}] Step 3: Extracting by row (cells) and transposing...")
print(f"Scanning 2B lines...")

strong_tmp = os.path.join(OUT_DIR, '_strong_tmp_v4.mtx')
weak_tmp = os.path.join(OUT_DIR, '_weak_tmp_v4.mtx')

strong_nnz = 0
weak_nnz = 0
line_count = 0

with gzip.open(MTX_PATH, 'rt') as f:
    # Skip header
    for line in f:
        if line.startswith('%'):
            continue
        dims = line.strip().split()
        n_rows_orig, n_cols_orig, n_nnz_total = int(dims[0]), int(dims[1]), int(dims[2])
        break
    
    # n_rows_orig = cells, n_cols_orig = genes
    n_genes = n_cols_orig
    
    with open(strong_tmp, 'w', buffering=8*1024*1024) as f_s, open(weak_tmp, 'w', buffering=8*1024*1024) as f_w:
        for line in f:
            line_count += 1
            if line_count % 100000000 == 0:
                print(f"  {line_count//1000000}M lines...")
            
            parts = line.split()
            row, col, val = int(parts[0]), int(parts[1]), float(parts[2])
            
            # row = cell index, col = gene index
            # Output: (gene, new_cell, val)
            if row in strong_set:
                f_s.write(f"{col} {strong_col_map[row]} {val}\n")
                strong_nnz += 1
            if row in weak_set:
                f_w.write(f"{col} {weak_col_map[row]} {val}\n")
                weak_nnz += 1

print(f"  Strong nnz: {strong_nnz} | Weak nnz: {weak_nnz}")

# Step 4: Write final matrices (genes x cells)
print(f"[{datetime.now()}] Step 4: Writing final matrices...")

strong_out = os.path.join(OUT_DIR, 'strong_matrix.mtx')
weak_out = os.path.join(OUT_DIR, 'weak_matrix.mtx')

with open(strong_out, 'w') as f_out:
    f_out.write("%%MatrixMarket matrix coordinate real general\n")
    f_out.write(f"{n_genes} {len(strong_rows)} {strong_nnz}\n")
    with open(strong_tmp, 'r') as f_in:
        f_out.write(f_in.read())

with open(weak_out, 'w') as f_out:
    f_out.write("%%MatrixMarket matrix coordinate real general\n")
    f_out.write(f"{n_genes} {len(weak_rows)} {weak_nnz}\n")
    with open(weak_tmp, 'r') as f_in:
        f_out.write(f_in.read())

os.remove(strong_tmp)
os.remove(weak_tmp)

# Step 5: Save metadata, barcodes, genes
print(f"[{datetime.now()}] Step 5: Saving metadata...")
strong_df.to_csv(os.path.join(OUT_DIR, 'strong_metadata.csv'), index=False)
weak_df.to_csv(os.path.join(OUT_DIR, 'weak_metadata.csv'), index=False)

strong_barcodes_out = [all_barcodes[i-1] for i in sorted(strong_rows)]
weak_barcodes_out = [all_barcodes[i-1] for i in sorted(weak_rows)]

with open(os.path.join(OUT_DIR, 'strong_barcodes.tsv'), 'w') as f:
    f.write('\n'.join(strong_barcodes_out) + '\n')
with open(os.path.join(OUT_DIR, 'weak_barcodes.tsv'), 'w') as f:
    f.write('\n'.join(weak_barcodes_out) + '\n')

with gzip.open(GENES_PATH, 'rt') as f:
    genes = [l.strip() for l in f][1:]
with open(os.path.join(OUT_DIR, 'genes.tsv'), 'w') as f:
    f.write('\n'.join(genes) + '\n')

with open(os.path.join(OUT_DIR, 'sampling_info.json'), 'w') as f:
    json.dump({
        'seed': SEED, 'n_per_group': N_PER_GROUP,
        'total_sampled': len(sampled_df),
        'strong_n': len(strong_rows), 'weak_n': len(weak_rows),
        'strong_nnz': strong_nnz, 'weak_nnz': weak_nnz,
        'n_genes': n_genes,
    }, f, indent=2)

print(f"\n[{datetime.now()}] DONE!")
