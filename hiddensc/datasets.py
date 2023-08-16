"""Everything related to loading and manipulating datasets."""
from typing import Tuple

import anndata
import numpy as np
import pandas as pd
import scanpy as sc

from . import types


def augment_for_analysis(adata: types.AnnData) -> None:
    """Augment dataset for visual analysis."""
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=1.2)

def preprocess_data(adata: types.AnnData) -> None:
    """Basic preprocessing for dataset."""
    sc.pp.filter_cells(adata, min_genes=20, inplace=True)
    sc.pp.filter_genes(adata, min_cells=1, inplace=True)
    adata.raw = adata

def normalize_and_log(adata: types.AnnData) -> None:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    
def setup_minimal_data(counts: np.ndarray,
                       gene_ids: pd.Series,
                       barcodes: pd.Series,
                       batch: pd.Series,
                       perturbed: pd.Series) -> types.AnnData:
    """Creates a minimal raw dataset in anndata format."""
    adata = sc.AnnData(counts, dtype=np.int32)
    adata.var_names = anndata.utils.make_index_unique(pd.Index(gene_ids))
    adata.var['gene_ids'] = gene_ids
    adata.obs['barcodes'] = pd.Index(barcodes, name=0)
    adata.layers["counts"] = adata.X.copy()  # preserve counts
    adata.obs['batch'] = batch
    adata.obs['perturbed'] = perturbed
    adata.raw = adata
    return adata

def _get_relevant_genes(name: str, df: pd.DataFrame) -> types.StrArray:
    """Helper function to extract relevant genes."""
    mask = np.logical_and(df[f'{name}_p'] < 0.05, df[f'{name}_l'] > 0)
    genes = df[f'{name}_n'][mask].values.astype(str)
    return genes

def get_de_genes(adata: types.AnnData, label: str, prefix: str = '') -> types.StrArrayDict:
    """Get relevant differentially expressed genes based on a label."""
    adata.obs[label] = adata.obs[label].astype('category')
    sc.tl.rank_genes_groups(adata, label, method='wilcoxon')
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    result_df = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
         for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})

    levels = np.unique(adata.obs[label])
    de_genes = {f'{prefix}{i}': _get_relevant_genes(i, result_df) for i in levels}
    return de_genes

def make_adata_with_hybrid_perturbed_cells(adata: types.AnnData, W_m: float) -> Tuple[types.AnnData, types.AnnData]:
    """Create an ablated dataset by perturbing cells.

    To create the Memory B / Naive B hybrids with weights W_m and W_n:=(1-W_m)
    for each Memory B cell indexed by i, denote the total number of UMIs with N_counts_i
    draw N_m_i = ceil(W_m * N_counts_i) from Multinom(gene_props_i), where gene_props_i is the s
    caled count vector across all genes for cell i
    draw N_n_i := N_counts_i = N_m_i from Multinom(NaiveB_centroid)
    sum the two draws.
    """

    print(f'Generating a dataset with {W_m} Memory B {1 - W_m} Naive B hybrids')
    adata_NaiveB = adata[adata.obs['perturbed'].isin(['Naive B'])]
    df_X_NaiveB = pd.DataFrame(adata_NaiveB.layers['counts'])  # raw expression
    df_X_NaiveB.columns = adata_NaiveB.var.index

    NaiveB_centroid = np.mean(df_X_NaiveB.divide(np.sum(df_X_NaiveB, axis=1), axis=0), axis=0)

    adata_MemoryB = adata[adata.obs['perturbed'].isin(['Memory B'])]

    W_n = 1 - W_m
    hybrid_MemoryB_cells = pd.DataFrame(np.nan, index=adata_MemoryB.obs['barcodes'],
                                        columns=adata_MemoryB.var['gene_ids'])
    for i, row in pd.DataFrame(adata_MemoryB.layers['counts']).iterrows():
        N_counts = sum(row)
        N_m = math.ceil(W_m * N_counts)
        N_n = N_counts - N_m
        hybrid_MemoryB_cells.iloc[i] = np.random.multinomial(N_m, row.values / N_counts) + np.random.multinomial(N_n,
                                                                                                                 NaiveB_centroid)

    adata.layers[f'counts_{W_m}'] = np.concatenate((adata_NaiveB.layers['counts'], hybrid_MemoryB_cells.values), axis=0)

    adata_hybrid = sc.AnnData(np.concatenate((adata_NaiveB.layers['counts'], hybrid_MemoryB_cells.values), axis=0))
    adata_hybrid.var_names = anndata.utils.make_index_unique(pd.Index(adata.var['gene_ids']))
    adata_hybrid.var['gene_ids'] = adata.var['gene_ids']
    adata_hybrid.obs['barcodes'] = pd.Index(adata.obs['barcodes'], name=0)
    adata_hybrid.raw = adata_hybrid

    # filter genes
    sc.pp.filter_genes(adata_hybrid, min_cells=1, inplace=True)
    mito_genes = adata_hybrid.var_names.str.startswith('mt-')
    adata_hybrid.obs['percent_mito'] = np.sum(adata_hybrid[:, mito_genes].X, axis=1) / np.sum(adata_hybrid.X, axis=1)
    adata_hybrid.obs['n_counts'] = adata_hybrid.X.sum(axis=1)
    sc.pp.normalize_total(adata_hybrid, target_sum=1e4)
    sc.pp.log1p(adata_hybrid)
    adata_hybrid.raw = adata_hybrid
    # sc.pp.highly_variable_genes(adata, min_mean=0.03, max_mean=5, min_disp=0.9) #min_mean=0.01, max_mean=4, min_disp=0.2
    # if plot_diagnostics:
    #    sc.pl.highly_variable_genes(adata)
    # adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata_hybrid, max_value=10)
    sc.tl.pca(adata_hybrid, svd_solver='arpack', n_comps=50)
    sc.pp.neighbors(adata_hybrid, n_neighbors=10, n_pcs=50)
    sc.tl.umap(adata_hybrid)
    sc.tl.leiden(adata_hybrid, resolution=1.2)

    adata_hybrid.obs['batch'] = adata.obs['batch']
    adata_hybrid.obs['perturbed'] = adata.obs['perturbed']

    return adata, adata_hybrid
