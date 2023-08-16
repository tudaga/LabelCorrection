#!/usr/bin/env python
"""Generate dimensionality reduced features."""
import os
import sys

root_path = os.path.abspath('./..')
sys.path.insert(0, root_path)

import numpy as np
import itertools
import functools
from tqdm import tqdm
import gc
import matplotlib.pyplot as plt

import hiddensc
from hiddensc import utils, files, vis, types

import scanpy as sc
import scvi
import anndata

utils.print_module_versions([sc, anndata, scvi, hiddensc])
vis.visual_settings()

OVERWRITE = False
EXP_IDS = files.MEMORY_B_VALUES
N_EPOCHS = 250


def create_vae_representations(adata: types.AnnData, n_epochs: int,
                               ks: list[int], plot_dir: str) -> dict[str, types.Array]:
    """Compute several flavors of a VAE over an anndata."""
    model_classes = [scvi.model.LinearSCVI, scvi.model.SCVI]
    combos = list(itertools.product(model_classes, ks))
    feats = {}
    for model_cls, k in tqdm(combos):
        # Reset seed for each model.
        utils.set_random_seed(utils.RANDOM_SEED)
        local_adata = adata.copy()
        name = f'{model_cls.__name__}_{k}'
        model_cls.setup_anndata(local_adata, layer="counts")
        model = model_cls(local_adata, n_latent=k)
        model.train(max_epochs=n_epochs, plan_kwargs={"lr": 5e-3},
                    check_val_every_n_epoch=5)
        train_elbo = model.history["elbo_train"][1:]
        test_elbo = model.history["elbo_validation"]
        ax = train_elbo.plot()
        test_elbo.plot(ax=ax)
        plt.yscale('log')
        plt.savefig(os.path.join(plot_dir, f'{name}.png'))
        plt.title(name)
        feats[name] = model.get_latent_representation()
        plt.clf()  # plt.show()
        del local_adata
    return feats


if __name__ == "__main__":
    print('== Generating features  ==')

    fname = os.path.join(root_path, 'figures', 'ablation', 'optimal_NUM_PCS_KS_dict.npz')
    optimal_ks = files.load_npz(fname)

    for exp_id in tqdm(EXP_IDS):
        utils.set_random_seed(utils.RANDOM_SEED)
        data_name = f'naiveB_1900_memoryB_{exp_id:d}'
        # Utility filename makers.
        at_results_dir = functools.partial(os.path.join, root_path, files.RESULT_DIR, data_name)
        at_train_dir = functools.partial(os.path.join, root_path, files.RESULT_DIR, data_name, 'training')
        os.makedirs(at_results_dir(), exist_ok=True)
        os.makedirs(at_train_dir(), exist_ok=True)
        fname = at_results_dir('dim_reduced.npz')
        if not OVERWRITE and os.path.exists(fname):
            print(f'Skipping {data_name}, found outputs and overwrite={OVERWRITE}')
            continue
            
        # Load dataset.
        main_datafile = os.path.join(root_path, files.DATA_DIR, f'{data_name}_raw.h5ad')
        adata = sc.read(main_datafile)
        hiddensc.datasets.preprocess_data(adata)
        hiddensc.datasets.normalize_and_log(adata)

        # Features.
        feats = {}
        # Reset seed for PCA.
        utils.set_random_seed(utils.RANDOM_SEED)
        x_pca = hiddensc.models.get_pca(adata)
        k = optimal_ks[data_name]
        feats['PCA'] = x_pca[:, :k]
        feats.update(create_vae_representations(adata,
                                                n_epochs=N_EPOCHS,
                                                ks=[10, 25, 50],
                                                plot_dir=at_train_dir()))
        np.savez_compressed(fname, **feats)
        # Cleanup.
        del adata
        gc.collect()
