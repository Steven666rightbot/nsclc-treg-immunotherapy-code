"""
Fast extraction of data for NicheNet analysis.
Uses direct numpy operations instead of scanpy's rank_genes_groups.
"""
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse
from scipy import stats

OUTDIR = "D:/Research/cucumber/fig5/nichenet_input"
import os
os.makedirs(OUTDIR, exist_ok=True)

adata = sc.read("D:/research/tomato/data/monocle3_input/treg_processed.h5ad")
print(f"Full data: {adata.shape}")

# Check data type
print(f"X dtype: {adata.X.dtype}")
print(f"X min/max: {adata.X.min():.4f}/{adata.X.max():.4f}")
print(f"Is log-transformed? {adata.X.max() < 50}")

# ============================================================
# 1. Split by subtype
# ============================================================
for subtype, label in [("CD4T_Treg_CCR8", "CCR8"), ("CD4T_Treg_MKI67", "MKI67")]:
    mask = adata.obs["sub_cell_type"] == subtype
    sub = adata[mask].copy()
    n_cells = sub.shape[0]
    print(f"\n{label}+ Treg: {sub.shape[0]} cells")
    
    # Get raw expression matrix
    X = sub.X
    if issparse(X):
        X = X.toarray()
    
    # Response groups
    resp_mask = (sub.obs["response"] == "Responder").values
    nr_mask = (sub.obs["response"] == "Non-responder").values
    
    print(f"  Responders: {resp_mask.sum()}, Non-responders: {nr_mask.sum()}")
    
    # ============================================================
    # Quick DE: Mann-Whitney U per gene
    # ============================================================
    genes = sub.var_names.tolist()
    de_results = []
    
    # Check a subset first (only protein-coding / well-annotated genes)
    # For speed, look at the top 3000 most variable genes
    gene_var = X.var(axis=0)
    top_var_idx = np.argsort(gene_var)[-3000:]
    
    print(f"  Running MWU on {len(top_var_idx)} HVGs...")
    
    for idx in top_var_idx:
        g = genes[idx]
        r_expr = X[resp_mask, idx]
        nr_expr = X[nr_mask, idx]
        stat, pval = stats.mannwhitneyu(r_expr, nr_expr, alternative='two-sided')
        logfc = np.mean(r_expr) - np.mean(nr_expr)
        de_results.append({"gene": g, "logFC": logfc, "MWU_stat": stat, "pval": pval})
    
    de_df = pd.DataFrame(de_results)
    # BH correction
    from scipy.stats import false_discovery_control
    de_df["pval_adj"] = false_discovery_control(de_df["pval"])
    de_df = de_df.sort_values("pval_adj")
    de_df.to_csv(f"{OUTDIR}/de_{label}_R_vs_NR.csv", index=False)
    
    n_sig = (de_df["pval_adj"] < 0.05).sum()
    print(f"  Significant DE genes (p.adj<0.05): {n_sig}")
    print(f"  Top 10 up in Responder:")
    print(de_df[de_df["logFC"] > 0].head(10)[["gene", "logFC", "pval_adj"]])
    print(f"  Top 10 up in Non-responder:")
    print(de_df[de_df["logFC"] < 0].head(10)[["gene", "logFC", "pval_adj"]])
    
    # ============================================================
    # Expressed genes (>5% cells)
    # ============================================================
    frac_expressed = (X > 0).mean(axis=0)
    expressed_genes = [g for g, frac in zip(genes, frac_expressed) if frac > 0.05]
    pd.Series(expressed_genes).to_csv(f"{OUTDIR}/expressed_genes_{label}.csv", index=False, header=False)
    print(f"  Expressed genes (>5% cells): {len(expressed_genes)}")
    
    # ============================================================
    # Save pseudobulk expression (sum by sample for sender/receiver)
    # ============================================================
    sub.obs["sample_response"] = sub.obs["sampleID"].astype(str) + "_" + sub.obs["response"].astype(str)
    
    # Pseudobulk: sum counts per sample
    samples = sub.obs["sample_response"].unique()
    pb_matrix = np.zeros((len(genes), len(samples)))
    for i, s in enumerate(samples):
        cell_mask = (sub.obs["sample_response"] == s).values
        pb_matrix[:, i] = X[cell_mask].sum(axis=0)
    
    pb_df = pd.DataFrame(pb_matrix, index=genes, columns=samples)
    pb_df.to_csv(f"{OUTDIR}/pseudobulk_{label}.csv.gz", compression="gzip")
    print(f"  Pseudobulk: {pb_df.shape}")
    print(f"  Samples: {len(samples)}")

print("\nDone!")
