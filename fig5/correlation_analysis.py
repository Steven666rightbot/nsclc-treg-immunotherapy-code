"""
Gene expression correlation analysis to validate CellChat-predicted signaling axes.

Approach:
1. In CCR8+ Treg: find genes correlated with CD8A/CD8B (MHC-I receptors)
   → These genes represent the "MHC-I response program"
2. In MKI67+ Treg: find genes correlated with CD74 (MIF receptor)
   → These genes represent the "MIF response program"
3. Compare correlation patterns between Responders and Non-responders
4. Show top correlated genes with pathway enrichment
"""
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os

OUTDIR = "D:/Research/cucumber/fig5/correlation_analysis"
os.makedirs(OUTDIR, exist_ok=True)

# Load data
adata = sc.read("D:/research/tomato/data/monocle3_input/treg_processed.h5ad")
print(f"Full data: {adata.shape}")

# ============================================================
# Analysis 1: CCR8+ Treg — CD8A/CD8B correlation
# ============================================================
print("\n=== CCR8+ Treg: MHC-I Receptor Correlation ===")
mask_ccr8 = adata.obs["sub_cell_type"] == "CD4T_Treg_CCR8"
ccr8 = adata[mask_ccr8].copy()
X_ccr8 = ccr8.X
if issparse(X_ccr8):
    X_ccr8 = X_ccr8.toarray()

# Receptor genes
receptors_ccr8 = ["CD8A", "CD8B", "CD4"]
# Check if receptors exist
receptors_ccr8 = [g for g in receptors_ccr8 if g in ccr8.var_names]
print(f"Receptors found: {receptors_ccr8}")

# For each receptor, find correlated genes
for receptor in receptors_ccr8:
    print(f"\n  Correlating with {receptor}...")
    ridx = list(ccr8.var_names).index(receptor)
    r_expr = X_ccr8[:, ridx]
    
    # Quick check: how many cells express it
    nz = (r_expr > 0).mean()
    print(f"    Expressed in {nz*100:.1f}% of cells")
    
    if nz < 0.05:
        print(f"    SKIP: {receptor} too lowly expressed")
        continue
    
    # Spearman correlation with ALL other genes
    cors = []
    for i in range(X_ccr8.shape[1]):
        if i == ridx:
            continue
        gene = ccr8.var_names[i]
        g_expr = X_ccr8[:, i]
        # Only compute if gene has reasonable expression
        if (g_expr > 0).mean() < 0.05:
            continue
        rho, pval = spearmanr(r_expr, g_expr)
        cors.append({"gene": gene, "rho": rho, "pval": pval})
    
    cors_df = pd.DataFrame(cors)
    cors_df["pval_adj"] = np.minimum(cors_df["pval"] * len(cors_df), 1.0)
    cors_df = cors_df.sort_values("rho", ascending=False)
    
    print(f"    Total correlated genes: {len(cors_df)}")
    print(f"    Top 10 positive:")
    print(cors_df.head(10)[["gene", "rho", "pval_adj"]].to_string(index=False))
    print(f"    Top 10 negative:")
    print(cors_df.tail(10).iloc[::-1][["gene", "rho", "pval_adj"]].to_string(index=False))
    
    cors_df.to_csv(f"{OUTDIR}/corr_CCR8_{receptor}.csv", index=False)

# ============================================================
# Analysis 2: MKI67+ Treg — CD74 correlation
# ============================================================
print("\n=== MKI67+ Treg: MIF Receptor (CD74) Correlation ===")
mask_mki67 = adata.obs["sub_cell_type"] == "CD4T_Treg_MKI67"
mki67 = adata[mask_mki67].copy()
X_mki67 = mki67.X
if issparse(X_mki67):
    X_mki67 = X_mki67.toarray()

receptors_mki67 = ["CD74", "CD44", "CXCR4"]
receptors_mki67 = [g for g in receptors_mki67 if g in mki67.var_names]
print(f"Receptors found: {receptors_mki67}")

for receptor in receptors_mki67:
    print(f"\n  Correlating with {receptor}...")
    ridx = list(mki67.var_names).index(receptor)
    r_expr = X_mki67[:, ridx]
    
    nz = (r_expr > 0).mean()
    print(f"    Expressed in {nz*100:.1f}% of cells")
    
    if nz < 0.05:
        print(f"    SKIP: {receptor} too lowly expressed")
        continue
    
    cors = []
    for i in range(X_mki67.shape[1]):
        if i == ridx:
            continue
        gene = mki67.var_names[i]
        g_expr = X_mki67[:, i]
        if (g_expr > 0).mean() < 0.05:
            continue
        rho, pval = spearmanr(r_expr, g_expr)
        cors.append({"gene": gene, "rho": rho, "pval": pval})
    
    cors_df = pd.DataFrame(cors)
    cors_df["pval_adj"] = np.minimum(cors_df["pval"] * len(cors_df), 1.0)
    cors_df = cors_df.sort_values("rho", ascending=False)
    
    print(f"    Total correlated genes: {len(cors_df)}")
    print(f"    Top 10 positive:")
    print(cors_df.head(10)[["gene", "rho", "pval_adj"]].to_string(index=False))
    print(f"    Top 10 negative:")
    print(cors_df.tail(10).iloc[::-1][["gene", "rho", "pval_adj"]].to_string(index=False))
    
    # Check for key pathway genes
    for key_gene in ["HLA-DRB5", "HLA-DRA", "HLA-DRB1", "HSPA1A", "HSPA1B", 
                      "HSP90AA1", "DNAJB1", "BAG3", "MIF", "CTLA4", "FOXP3"]:
        if key_gene in cors_df["gene"].values:
            row = cors_df[cors_df["gene"] == key_gene]
            print(f"    {key_gene}: rho={row['rho'].values[0]:.4f}, p.adj={row['pval_adj'].values[0]:.2e}")
    
    cors_df.to_csv(f"{OUTDIR}/corr_MKI67_{receptor}.csv", index=False)

# ============================================================
# Plot: top correlated genes for key receptors
# ============================================================
def plot_top_correlations(cors_csv, title, receptor, n_top=20, outfile=None):
    """Plot top positively and negatively correlated genes."""
    df = pd.read_csv(cors_csv)
    df = df.sort_values("rho", ascending=False)
    
    # Take top positive and top negative
    top_pos = df.head(n_top).iloc[::-1]  # reverse for horizontal bar
    top_neg = df.tail(n_top)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, max(6, n_top * 0.3)))
    
    # Positive correlations
    ax = axes[0]
    ax.barh(range(len(top_pos)), top_pos["rho"].values, color="#2E5AAC", height=0.7)
    ax.set_yticks(range(len(top_pos)))
    ax.set_yticklabels(top_pos["gene"].values, fontsize=8)
    ax.set_xlabel("Spearman ρ")
    ax.set_title(f"Top + correlated with {receptor}", fontsize=10)
    ax.axvline(0, color="gray", linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Negative correlations
    ax = axes[1]
    ax.barh(range(len(top_neg)), top_neg["rho"].values, color="#B83A3A", height=0.7)
    ax.set_yticks(range(len(top_neg)))
    ax.set_yticklabels(top_neg["gene"].values, fontsize=8)
    ax.set_xlabel("Spearman ρ")
    ax.set_title(f"Top − correlated with {receptor}", fontsize=10)
    ax.axvline(0, color="gray", linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    fig.suptitle(title, fontsize=12, fontweight="bold")
    fig.tight_layout()
    
    if outfile:
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
        print(f"  Saved: {outfile}")
    plt.close(fig)

# Generate plots
if os.path.exists(f"{OUTDIR}/corr_CCR8_CD8A.csv"):
    plot_top_correlations(
        f"{OUTDIR}/corr_CCR8_CD8A.csv",
        "CCR8+ Treg: Genes Correlated with CD8A (MHC-I Receptor)",
        "CD8A", n_top=20,
        outfile=f"{OUTDIR}/fig_corr_CCR8_CD8A.png"
    )

if os.path.exists(f"{OUTDIR}/corr_MKI67_CD74.csv"):
    plot_top_correlations(
        f"{OUTDIR}/corr_MKI67_CD74.csv",
        "MKI67+ Treg: Genes Correlated with CD74 (MIF Receptor)",
        "CD74", n_top=20,
        outfile=f"{OUTDIR}/fig_corr_MKI67_CD74.png"
    )

print("\nDone!")
