"""
Statistical enrichment of key pathway genes in NicheNet target rankings.
Tests if our CellChat-relevant genes are ranked higher than random expectation.
"""
import pyreadr
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, hypergeom

print("Loading NicheNet matrix...")
ltm = pyreadr.read_r('D:/Research/tomato/data/nichenet_ligand_target_matrix_uncompressed.rds')[None]

# Our key gene sets (from CellChat and Fig 3 findings)
gene_sets = {
    'MHC-II_antigen': ['HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DPA1', 'HLA-DQA1', 'HLA-DQB1', 'CD74'],
    'MIF_pathway': ['CD74', 'CD44', 'MIF', 'CXCR4'],
    'HSF1_stress': ['HSPA1A', 'HSPA1B', 'HSP90AA1', 'DNAJB1', 'BAG3', 'HSPH1', 'HSPB1'],
    'Treg_function': ['FOXP3', 'IL2RA', 'CTLA4', 'TIGIT', 'ICOS', 'TNFRSF18'],
    'Proliferation': ['MKI67', 'TOP2A', 'PCNA', 'CDK1'],
    'Signaling': ['NFKB1', 'JUN', 'FOS', 'RELA', 'STAT2'],
}

def enrichment_test(ranking, gene_set, all_genes, top_n=500):
    """Test if gene_set is enriched in top_n of ranking."""
    n_total = len(all_genes)
    n_in_set = len([g for g in gene_set if g in all_genes])
    n_top = top_n
    
    # Count overlaps
    top_genes = set(ranking.head(top_n).index)
    overlap = len(top_genes & set(gene_set))
    
    # Hypergeometric test
    if n_in_set > 0:
        pval = hypergeom.sf(overlap - 1, n_total, n_in_set, n_top)
    else:
        pval = 1.0
    
    return overlap, n_in_set, pval

for ligand_name, ligand_col in [('MHC-I (HLA-A/B/C mean)', None), ('MIF', 'MIF')]:
    print(f"\n=== {ligand_name} ===")
    
    if ligand_col is None:
        # Average across HLA-A, B, C
        cols = [c for c in ['HLA-A', 'HLA-B', 'HLA-C'] if c in ltm.columns]
        ranking = ltm[cols].mean(axis=1).sort_values(ascending=False)
    else:
        ranking = ltm[ligand_col].sort_values(ascending=False)
    
    all_genes = list(ranking.index)
    
    for set_name, genes in gene_sets.items():
        overlap, n_in_set, pval = enrichment_test(ranking, genes, all_genes, top_n=500)
        print(f"  {set_name:20s}: {overlap}/{n_in_set} in top500, p={pval:.4f}")
        
        # Also show individual gene ranks
        for g in genes:
            if g in ranking.index:
                rank = list(ranking.index).index(g) + 1
                score = ranking[g]
                percentile = rank / len(ranking) * 100
                marker = " ★" if overlap > 0 and g in set(ranking.head(500).index) else ""
                print(f"    {g:15s}: rank={rank}/{len(ranking)} ({percentile:.1f}%), score={score:.6f}{marker}")

print("\n\n=== Summary for manuscript: ===")
print("Looking at whether MIF or MHC-I signaling is predicted to regulate pathways of interest:")
print()

for ligand_name, ligand_col in [('MHC-I (HLA-A/B/C)', None), ('MIF', 'MIF')]:
    print(f"\n{ligand_name}:")
    if ligand_col is None:
        cols = [c for c in ['HLA-A', 'HLA-B', 'HLA-C'] if c in ltm.columns]
        ranking = ltm[cols].mean(axis=1).sort_values(ascending=False)
    else:
        ranking = ltm[ligand_col].sort_values(ascending=False)
    
    all_genes = list(ranking.index)
    
    for set_name, genes in gene_sets.items():
        overlap, n_in_set, pval = enrichment_test(ranking, genes, all_genes, top_n=500)
        if pval < 0.05:
            print(f"  ✓ {set_name}: enriched in top500 (p={pval:.4f})")
        else:
            print(f"  ✗ {set_name}: not enriched (p={pval:.4f})")
