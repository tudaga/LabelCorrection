#!/usr/bin/env python
"""Evaluate predictions as 1) binary classification problem and 2) with de genes."""

import gc
import os
import sys

root_path = os.path.abspath('./..')
sys.path.insert(0, root_path)
import pandas as pd
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
OVERWRITE = True

if __name__ == "__main__":
    print('== Evaluating predictions over several datasets ==')
    for exp_id in tqdm(EXP_IDS):
        utils.set_random_seed(utils.RANDOM_SEED)
        
        # Setup files.
        data_name = f'naiveB_1900_memoryB_{exp_id:d}'
        at_results_dir = functools.partial(os.path.join, root_path, files.RESULT_DIR, data_name)
        if not os.path.exists(at_results_dir('predictions.csv')):
            print(f'Skipping {data_name}, did not find predictions')
            continue
        
        if os.path.exists(at_results_dir('de_genes.csv')) and not OVERWRITE:
            print(f'Skipping {data_name}, found outputs and overwrite={OVERWRITE}')
            continue
        
        # Load dataset.
        main_datafile = os.path.join(root_path, files.DATA_DIR, f'{data_name}_raw.h5ad')
        adata = sc.read(main_datafile)
        hiddensc.datasets.preprocess_data(adata)
        
        # Load predictions.
        batch, perturbed, pred_df = files.load_predictions(at_results_dir())
        # Calculate classification performance.
        results = []
        for (dim_reduce, method), grp_df in pred_df.groupby(axis=1, level=[0, 1]):
            p_hat = grp_df[(dim_reduce, method, 'p_hat')]
            p_label = grp_df[(dim_reduce, method, 'p_label')]
            stat = {'Dataset': data_name,
                    'DimReduction': dim_reduce,
                    'PredictionMethod': method}
            stat.update(hiddensc.metrics.evaluate_classification(perturbed, p_hat, p_label))
            results.append(stat)

        results_df = pd.DataFrame(results)
        results_df.to_csv(at_results_dir('classification_stats.csv'), index=False)

        # Differentially expressed genes (DE).
        # Calculate on reference data.
        de_genes = hiddensc.datasets.get_de_genes(adata, 'perturbed')
        de_genes.update(hiddensc.datasets.get_de_genes(adata, 'batch'))
        results = []
        info = {'Dataset': data_name,
                'True labels': 'Naive B/Memory B',
                'DimReduction': 'NA',
                'PredictionMethod': 'NA',
                'Approach': 'sample'}
        info.update(hiddensc.metrics.evaluate_de_genes(de_genes,
                                                       'Naive B', 'Memory B',
                                                       'Control', 'Case'))
        results.append(info)
        # Calculate on refined labels.
        for (dim_reduce, method), grp_df in pred_df.groupby(axis=1, level=[0, 1]):
            index = (dim_reduce, method, 'p_label')
            new_label = '_'.join(index[:2])
            adata.obs['refined_labels'] = grp_df[index].values.astype(int).astype(str)
            de_genes.update(hiddensc.datasets.get_de_genes(adata, 'refined_labels', f'{new_label}_'))
            info = {'Dataset': data_name,
                    'True labels': 'Naive B/Memory B',
                    'DimReduction': dim_reduce,
                    'PredictionMethod': method,
                    'Approach': new_label}
            info.update(
                hiddensc.metrics.evaluate_de_genes(de_genes,
                                                   'Naive B', 'Memory B',
                                                   f'{new_label}_0', f'{new_label}_1'))
            results.append(info)
        # DE results
        files.save_npz(at_results_dir('de_genes.npz'), de_genes)
        de_df = pd.DataFrame(results)
        de_df.to_csv(at_results_dir('de_genes.csv'), index=False)
        # Cleanup memory.
        del adata
        gc.collect()
