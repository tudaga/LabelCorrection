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
        if not os.path.exists(at_results_dir('dim_reduced.npz')):
            print(f'Did not find features file for {data_name}')
            continue
        if os.path.exists(at_results_dir('predictions.csv')) and not OVERWRITE:
            print(f'Skipping {data_name}, found outputs and overwrite={OVERWRITE}')
            continue

        # Load dataset.
        main_datafile = os.path.join(root_path, files.DATA_DIR, f'{data_name}_raw.h5ad')
        adata = sc.read(main_datafile)
        hiddensc.datasets.preprocess_data(adata)
        hiddensc.datasets.normalize_and_log(adata)
    
        # Load dim, setup data.
        feats = files.load_npz(at_results_dir('dim_reduced.npz'))
        y = (adata.obs['batch'].values == 'Case').astype(np.int32)
        y_true = (adata.obs['perturbed'].values == 'Memory B').astype(np.int32)
        ids = adata.obs['barcodes'].values
        # Prediction models.
        pred_fns = {'Logistic': hiddensc.models.logistic_predictions,
                    'SVM': hiddensc.models.svm_predictions}
        # Initial ground truth data.
        preds = [y, y_true]
        info = [('batch', '', 'Case'), ('perturbed', '', 'Memory B')]
        # Loop over features and prediction models.
        combos = list(itertools.product(feats.keys(), pred_fns.keys()))
        for feat_name, strat_name in tqdm(combos):
            rand_state = 0
            x = feats[feat_name]
            p_hat, p_labels = pred_fns[strat_name](x, y, 1, rand_state)
            preds.append(p_hat)
            info.append((feat_name, strat_name, 'p_hat'))
            preds.append(p_labels)
            info.append((feat_name, strat_name, 'p_label'))

        # Assemble prediction csv.
        cols = pd.MultiIndex.from_tuples(info)
        pred_df = pd.DataFrame(np.array(preds).T, index=adata.obs['barcodes'], columns=cols)
        pred_df.to_csv(at_results_dir('predictions.csv'))
        # Cleanup.
        del adata
        gc.collect()
