"""
Prepare GSE243013 CellChat input with stratified downsampling.
Each sub_cell_type x response_group: sample up to 1000 cells.
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
print(f"[{datetime.now()}] Reading metadata...")
chunks = []
for chunk in pd.read_csv(META_PATH, low_memory=False, chunksize=20000,
                          usecols=['cellID','sampleID','major_cell_type','sub_cell_type','pathological_response']):
    chunks.append(chunk)
meta = pd.concat(chunks)

meta_valid = meta[meta['pathological_response'].isin(['MPR', 'pCR', 'non-MPR'])].copy()
meta_valid['response_group'] = meta_valid['pathological_response'].apply(
    lambda x: 'strong' if x in ['MPR', 'pCR'] else 'weak')

print(f"Total valid cells: {len(meta_valid)}")
print(f"Strong: {sum(meta_valid['response_group']=='strong')}")
print(f"Weak: {sum(meta_valid['response_group']=='weak')}")

# Stratified sampling
sampled = []
for (ct, rg), group in meta_valid.groupby(['sub_cell_type', 'response_group']):
    n = min(N_PER_GROUP, len(group))
    sel = group.sample(n=n, random_state=SEED)
    sampled.append(sel)
    print(f"  {ct} / {rg}: {len(group)} -> {n}")

sampled_df = pd.concat(sampled)
print(f"\nTotal sampled: {len(sampled_df)}")
print(f"  Strong: {sum(sampled_df['response_group']=='strong')}")
print(f"  Weak: {sum(sampled_df['response_group']=='weak')}")

# Split by group
strong_df = sampled_df[sampled_df['response_group']=='strong'][['cellID','sub_cell_type','response_group']].copy()
weak_df = sampled_df[sampled_df['response_group']=='weak'][['cellID','sub_cell_type','response_group']].copy()

strong_df.columns = ['cell', 'cell_type', 'response_group']
weak_df.columns = ['cell', 'cell_type', 'response_group']

# ===================== Step 2: Map barcodes to column indices =====================
print(f"\n[{datetime.now()}] Mapping barcodes to column indices...")
with gzip.open(BARCODES_PATH, 'rt') as f:
    all_barcodes = [l.strip() for l in f][1:]  # Skip header

barcode_to_idx = {b: i+1 for i, b in enumerate(all_barcodes)}  # 1-based

strong_cells = set(strong_df['cell'].values)
weak_cells = set(weak_df['cell'].values)

strong_cols = sorted([barcode_to_idx[c] for c in strong_cells if c in barcode_to_idx])
weak_cols = sorted([barcode_to_idx[c] for c in weak_cells if c in barcode_to_idx])

print(f"Strong columns: {len(strong_cols)}")
print(f"Weak columns: {len(weak_cols)}")

# Save mappings
with open(os.path.join(OUT_DIR, 'sampling_info.json'), 'w') as f:
    json.dump({
        'seed': SEED,
        'n_per_group': N_PER_GROUP,
        'total_sampled': len(sampled_df),
        'strong_n': len(strong_cols),
        'weak_n': len(weak_cols),
        'strong_cell_types': strong_df['cell_type'].value_counts().to_dict(),
        'weak_cell_types': weak_df['cell_type'].value_counts().to_dict(),
    }, f, indent=2)

# ===================== Step 3: Extract expression matrix =====================
print(f"\n[{datetime.now()}] Extracting expression matrix...")
print(f"This will take 10-30 minutes (scanning 2B entries)...")

strong_out = os.path.join(OUT_DIR, 'strong_matrix.mtx')
weak_out = os.path.join(OUT_DIR, 'weak_matrix.mtx')

strong_set = set(strong_cols)
weak_set = set(weak_cols)

strong_entries = []
weak_entries = []
strong_nnz = 0
weak_nnz = 0

with gzip.open(MTX_PATH, 'rt') as f:
    # Read header
    header = []
    for line in f:
        if line.startswith('%'):
            header.append(line)
        else:
            header.append(line)
            break
    
    # Parse dimensions
    dims = header[-1].strip().split()
    n_rows, n_cols, n_nnz = int(dims[0]), int(dims[1]), int(dims[2])
    print(f"Original: {n_rows} genes x {n_cols} cells, {n_nnz} non-zero")
    
    # Scan all entries
    for i, line in enumerate(f):
        if i % 50000000 == 0 and i > 0:
            print(f"  Processed {i}M entries...")
        parts = line.strip().split()
        row, col, val = int(parts[0]), int(parts[1]), float(parts[2])
        
        if col in strong_set:
            strong_entries.append((row, col, val))
            strong_nnz += 1
        if col in weak_set:
            weak_entries.append((row, col, val))
            weak_nnz += 1

print(f"\nStrong: {strong_nnz} non-zero entries")
print(f"Weak: {weak_nnz} non-zero entries")

# Remap column indices to consecutive 1-based
strong_col_map = {old: new for new, old in enumerate(sorted(strong_cols), 1)}
weak_col_map = {old: new for new, old in enumerate(sorted(weak_cols), 1)}

# Write strong matrix
with open(strong_out, 'w') as f:
    f.write("%%MatrixMarket matrix coordinate real general\n")
    f.write(f"{n_rows} {len(strong_cols)} {strong_nnz}\n")
    for row, col, val in strong_entries:
        f.write(f"{row} {strong_col_map[col]} {val}\n")

# Write weak matrix
with open(weak_out, 'w') as f:
    f.write("%%MatrixMarket matrix coordinate real general\n")
    f.write(f"{n_rows} {len(weak_cols)} {weak_nnz}\n")
    for row, col, val in weak_entries:
        f.write(f"{row} {weak_col_map[col]} {val}\n")

print(f"\nSaved matrices to {OUT_DIR}")

# ===================== Step 4: Save metadata =====================
strong_df.to_csv(os.path.join(OUT_DIR, 'strong_metadata.csv'), index=False)
weak_df.to_csv(os.path.join(OUT_DIR, 'weak_metadata.csv'), index=False)

# Save barcodes (remapped order)
strong_barcodes_out = [all_barcodes[i-1] for i in sorted(strong_cols)]
weak_barcodes_out = [all_barcodes[i-1] for i in sorted(weak_cols)]

with open(os.path.join(OUT_DIR, 'strong_barcodes.tsv'), 'w') as f:
    f.write('\n'.join(strong_barcodes_out) + '\n')
with open(os.path.join(OUT_DIR, 'weak_barcodes.tsv'), 'w') as f:
    f.write('\n'.join(weak_barcodes_out) + '\n')

# Copy genes
import shutil
with gzip.open(GENES_PATH, 'rt') as f:
    genes = [l.strip() for l in f][1:]
with open(os.path.join(OUT_DIR, 'genes.tsv'), 'w') as f:
    f.write('\n'.join(genes) + '\n')

print(f"\n[{datetime.now()}] All done!")
print(f"Output dir: {OUT_DIR}")
