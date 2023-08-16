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


def run_cna_association(data, id_name, assoc_name):
    print(f' Run CNA with {id_name} as "id"', 3)
    print(f'Looking for association with {assoc_name}')
    data.obs['id'] = data.obs[id_name]
    d = MultiAnnData(data, sampleid='id')
    d.obs['ground_truth_labels'] = data.obs['ground_truth_labels'] == 'Memory B'
    d.obs_to_sample(['sample_level_labels', 'ground_truth_labels'])
        
    res = association(d, d.samplem[assoc_name])
    data.obs['CNA_result'] = res.ncorrs
    print(f'global association ({id_name}, {assoc_name}) p-value: {res.p}')
    
    del data.obs['CNA_result']
    del data.obs['id']
    return res

def map_to_binary(values, case_control_labels, case_cond=1):
    
    # simple thresholding works, since the values are either very close to -1 or very close to 1
    # return np.interp(values, [np.min(values),np.max(values)], [0,1])>=0.5
    
    # but for consistency,
    # using kmeans with n_clusters=2, same as with HiDDEN p_hat binarization of the cells in the case_cond
    
    kmeans_values_res = KMeans(n_clusters=2, random_state=0).fit(pd.DataFrame(values))
    mean_values_res_kmeans_label0 = np.mean(values[kmeans_values_res.labels_==0]) 
    mean_values_res_kmeans_label1 = np.mean(values[kmeans_values_res.labels_==1])
    zero_lab_has_lower_mean = mean_values_res_kmeans_label0 < mean_values_res_kmeans_label1

    df_values_clust = pd.DataFrame(values)
    df_values_clust['kmeans'] = 0
    df_values_clust['kmeans'][(case_control_labels==case_cond).values] = [1 if x==int(zero_lab_has_lower_mean) else 0 for x in kmeans_values_res.labels_]
    
    return df_values_clust['kmeans'].values

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
        
        # CNA preprocessing.
        adata.X = adata.layers['counts']
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, subset=True)
        adata.raw = adata
        
        adata.obs['sample_level_labels'] = (adata.obs['batch'].values.to_numpy() == 'Case').astype(int)
        adata.obs['ground_truth_labels'] = (adata.obs['perturbed'].values.to_numpy() == 'Memory B').astype(int)

        res = run_cna_association(adata, 'sample_level_labels','sample_level_labels')
        p_hat = res.ncorrs
        p_label = map_to_binary(res.ncorrs, adata.obs['sample_level_labels'])
    
        # Data storage.
        info = [('CNA', 'CNA', 'p_hat'), ('CNA', 'CNA', 'p_label')]
        preds = [ p_hat, p_label]
        ids = adata.obs['barcodes'].values
        
        cols = pd.MultiIndex.from_tuples(info)
        pred_df = pd.DataFrame(np.array(preds).T, index=adata.obs['barcodes'], columns=cols)
        pred_df.to_csv(at_results_dir('cna_predictions.csv'))
        
        # Cleanup.
        del adata
        gc.collect()
