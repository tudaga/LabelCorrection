#!/usr/bin/env python
"""Create minimal versions of memory B datasets."""
import gc
import os
import sys

root_path = os.path.abspath('./..')
sys.path.insert(0, root_path)

import functools
import pandas as pd
from tqdm.auto import tqdm

import hiddensc
from hiddensc import utils, files

import scanpy as sc
import anndata

EXP_IDS = files.MEMORY_B_VALUES


def load_memoryb_data_folder(dirname: str) -> sc.AnnData:
    """Load memoryB dataset in its original structure."""
    at_data_dir = functools.partial(os.path.join, dirname)
    data = sc.read(at_data_dir('citeseq_rna_adata.h5ad'))
    batch = pd.read_csv(at_data_dir('batch.csv'))['x'].values
    perturbed = pd.read_csv(at_data_dir('celltype.csv'))['x'].values
    data.obs['batch'] = ['Case' if int(b) else 'Control' for b in batch]
    data.obs['perturbed'] = perturbed
    raw_counts_allgenes = pd.read_csv(at_data_dir('counts_allgenes.csv'))
    data.layers['lognorm'] = data.X
    data.layers['counts'] = raw_counts_allgenes.values
    return data


if __name__ == "__main__":
    utils.print_module_versions([sc, anndata, hiddensc])
    for exp_id in tqdm(EXP_IDS):
        data_name = f'naiveB_1900_memoryB_{exp_id:d}'
        rawdata_dir = os.path.join(root_path, files.RAW_DIR, data_name)
        data_dir = os.path.join(root_path, files.DATA_DIR)
        fname = os.path.join(data_dir, f'{data_name}_raw.h5ad')
        # if os.path.exists(fname):
        #    continue
        print(f'Generating data in {data_dir}')
        os.makedirs(data_dir, exist_ok=True)
        memb_data = load_memoryb_data_folder(rawdata_dir)
        adata = hiddensc.datasets.setup_minimal_data(counts=memb_data.layers['counts'],
                                                     gene_ids=memb_data.var['gene_ids'],
                                                     barcodes=memb_data.obs['barcodes'],
                                                     batch=memb_data.obs['batch'],
                                                     perturbed=memb_data.obs['perturbed'])
        files.save_compressed_h5ad(fname, adata)
        del adata
        del memb_data
        gc.collect()
