#!/usr/bin/env python
"""Load dim reduced features and apply a predictive model."""

import gc
import os
import sys

root_path = os.path.abspath('./..')
sys.path.insert(0, root_path)
import pandas as pd
import numpy as np
import itertools
import functools
from tqdm import tqdm
import hiddensc
from hiddensc import utils, files

import scanpy as sc
import scvi
import anndata
import _mixscape as mixscape
from sklearn.exceptions import ConvergenceWarning
import warnings

warnings.simplefilter("ignore", category=ConvergenceWarning)

utils.print_module_versions([sc, anndata, scvi, hiddensc])

EXP_IDS = files.MEMORY_B_VALUES
OVERWRITE = False

if __name__ == "__main__":
    print('== Making predictions ==')
    for exp_id in tqdm(EXP_IDS):
        utils.set_random_seed(utils.RANDOM_SEED)
        
        # Check and setup files.
        data_name = f'naiveB_1900_memoryB_{exp_id:d}'
        at_results_dir = functools.partial(os.path.join, root_path, files.RESULT_DIR, data_name)

        if not os.path.exists(at_results_dir('predictions.csv')):
            print(f'Skipping {data_name}, expected already computed predictions')
            #found outputs and overwrite={OVERWRITE}')
            continue

        # Load dataset.
        main_datafile = os.path.join(root_path, files.DATA_DIR, f'{data_name}_raw.h5ad')
        adata = sc.read(main_datafile)
        
        # Mixscape preprocessing.
        adata.X = adata.layers['counts']
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, subset=True)
        adata.raw = adata
        
        # Mixscape predictions
        mixscape_identifier = mixscape.Mixscape()
        mixscape_identifier.perturbation_signature(adata=adata, pert_key='batch', control='Control')
        mixscape_identifier.mixscape(adata=adata, labels='batch', control='Control')
        p_hat = adata.obs['mixscape_class_p_ko'].astype(np.float)
        p_label = adata.obs['mixscape_class'] != 'Control'
        
        # Data storage.
        info = [('PCA', 'Mixscape', 'p_hat'), ('PCA', 'Mixscape', 'p_label')]
        preds = [ p_hat, p_label]
        ids = adata.obs['barcodes'].values
        
        cols = pd.MultiIndex.from_tuples(info)
        pred_df = pd.DataFrame(np.array(preds).T, index=adata.obs['barcodes'], columns=cols)
        pred_df.to_csv(at_results_dir('mixscape_predictions.csv'))
        
        # Cleanup.
        del adata
        gc.collect()
