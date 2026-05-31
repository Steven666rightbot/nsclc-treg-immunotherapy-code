"""
Prepare GSE243013 CellChat input with stratified downsampling.
Stream-based extraction to minimize memory usage.
"""
import os, gzip, random, json
from datetime import datetime
import pandas as pd

# Paths
META_PATH = r"D:\Research\tomato\data\temp_umap_all\GSE243013_metadata.csv"
BARCODES_PATH = r"D:\Research\tomato\data\geo_cache\GSE243013_barcodes.csv.gz"
MTX_PATH = r"D:\Research\tomato\data\geo_cache\GSE243013_NSCLC_immune_scRNA_counts.mtx.gz"
GENES_PATH = r"D:\Research\tomato\data\geo_cache\GSE243013_genes.csv.gz"
OUT_DIR = r"D:\Research\tomato\data\gse243013_cellchat_input"

SEED = 20250517
N_PER_GROUP = 1000

os.makedirs(OUT_DIR, exist_ok=True)
random.seed(SEED)

# ===================== Step 1: Read metadata, stratified sampling =====================
print(f"[{datetime.now()}] Step 1: Reading metadata...")
chunks = []
for chunk in pd.read_csv(META_PATH, low_memory=False, chunksize=20000,
                          usecols=['cellID','sampleID','major_cell_type','sub_cell_type','pathological_response']):
    chunks.append(chunk)
meta = pd.concat(chunks)

meta_valid = meta[meta['pathological_response'].isin(['MPR', 'pCR', 'non-MPR'])].copy()
meta_valid['response_group'] = meta_valid['pathological_response'].apply(
    lambda x: 'strong' if x in ['MPR', 'pCR'] else 'weak')

print(f"Total valid: {len(meta_valid)} | Strong: {sum(meta_valid['response_group']=='strong')} | Weak: {sum(meta_valid['response_group']=='weak')}")

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

# ===================== Step 2: Map barcodes =====================
print(f"[{datetime.now()}] Step 2: Mapping barcodes...")
with gzip.open(BARCODES_PATH, 'rt') as f:
    all_barcodes = [l.strip() for l in f][1:]

barcode_to_idx = {b: i+1 for i, b in enumerate(all_barcodes)}
strong_cells = set(strong_df['cell'].values)
weak_cells = set(weak_df['cell'].values)

strong_cols = sorted([barcode_to_idx[c] for c in strong_cells if c in barcode_to_idx])
weak_cols = sorted([barcode_to_idx[c] for c in weak_cells if c in barcode_to_idx])

strong_set = set(strong_cols)
weak_set = set(weak_cols)

print(f"Strong cols: {len(strong_cols)} | Weak cols: {len(weak_cols)}")

# Column remapping
strong_col_map = {old: new for new, old in enumerate(sorted(strong_cols), 1)}
weak_col_map = {old: new for new, old in enumerate(sorted(weak_cols), 1)}

# ===================== Step 3: First pass - count nnz =====================
print(f"[{datetime.now()}] Step 3: Pass 1 - counting non-zero entries...")
strong_nnz = 0
weak_nnz = 0
line_count = 0

with gzip.open(MTX_PATH, 'rt') as f:
    for line in f:
        if line.startswith('%'):
            continue
        dims = line.strip().split()
        n_rows, n_cols, n_nnz_total = int(dims[0]), int(dims[1]), int(dims[2])
        break
    
    for line in f:
        line_count += 1
        if line_count % 100000000 == 0:
            print(f"  {line_count//1000000}M lines scanned...")
        parts = line.strip().split()
        col = int(parts[1])
        if col in strong_set:
            strong_nnz += 1
        if col in weak_set:
            weak_nnz += 1

print(f"  Strong nnz: {strong_nnz} | Weak nnz: {weak_nnz}")
print(f"  Total lines: {line_count}")

# ===================== Step 4: Second pass - write matrices =====================
print(f"[{datetime.now()}] Step 4: Pass 2 - writing matrices...")

strong_out = os.path.join(OUT_DIR, 'strong_matrix.mtx')
weak_out = os.path.join(OUT_DIR, 'weak_matrix.mtx')

with gzip.open(MTX_PATH, 'rt') as f:
    for line in f:
        if line.startswith('%'):
            continue
        break
    
    with open(strong_out, 'w') as f_s, open(weak_out, 'w') as f_w:
        f_s.write("%%MatrixMarket matrix coordinate real general\n")
        f_s.write(f"{n_rows} {len(strong_cols)} {strong_nnz}\n")
        f_w.write("%%MatrixMarket matrix coordinate real general\n")
        f_w.write(f"{n_rows} {len(weak_cols)} {weak_nnz}\n")
        
        for i, line in enumerate(f):
            if i % 100000000 == 0 and i > 0:
                print(f"  {i//1000000}M lines written...")
            parts = line.strip().split()
            row, col, val = int(parts[0]), int(parts[1]), float(parts[2])
            if col in strong_set:
                f_s.write(f"{row} {strong_col_map[col]} {val}\n")
            if col in weak_set:
                f_w.write(f"{row} {weak_col_map[col]} {val}\n")

# ===================== Step 5: Save metadata, barcodes, genes =====================
print(f"[{datetime.now()}] Step 5: Saving metadata...")
strong_df.to_csv(os.path.join(OUT_DIR, 'strong_metadata.csv'), index=False)
weak_df.to_csv(os.path.join(OUT_DIR, 'weak_metadata.csv'), index=False)

strong_barcodes_out = [all_barcodes[i-1] for i in sorted(strong_cols)]
weak_barcodes_out = [all_barcodes[i-1] for i in sorted(weak_cols)]

with open(os.path.join(OUT_DIR, 'strong_barcodes.tsv'), 'w') as f:
    f.write('\n'.join(strong_barcodes_out) + '\n')
with open(os.path.join(OUT_DIR, 'weak_barcodes.tsv'), 'w') as f:
    f.write('\n'.join(weak_barcodes_out) + '\n')

with gzip.open(GENES_PATH, 'rt') as f:
    genes = [l.strip() for l in f][1:]
with open(os.path.join(OUT_DIR, 'genes.tsv'), 'w') as f:
    f.write('\n'.join(genes) + '\n')

# Save sampling info
with open(os.path.join(OUT_DIR, 'sampling_info.json'), 'w') as f:
    json.dump({
        'seed': SEED,
        'n_per_group': N_PER_GROUP,
        'total_sampled': len(sampled_df),
        'strong_n': len(strong_cols),
        'weak_n': len(weak_cols),
        'strong_nnz': strong_nnz,
        'weak_nnz': weak_nnz,
        'n_genes': n_rows,
    }, f, indent=2)

print(f"\n[{datetime.now()}] DONE! Output: {OUT_DIR}")
