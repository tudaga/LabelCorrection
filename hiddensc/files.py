from typing import Tuple

import hdf5plugin
import numpy as np
import pandas as pd

from . import types

RAW_DIR = 'raw_data'
DATA_DIR = 'data'
RESULT_DIR = 'results'

MEMORY_B_VALUES = [49, 100, 154, 211, 271, 335, 403, 475,
                   552, 633, 721, 814, 915, 1023, 1140, 1267, 1404, 1555]


def check_ext(fname: str, ext: str) -> None:
    """Check if extension of file is correct."""
    if not fname.endswith(ext):
        raise ValueError(f'Expected extension "{ext}" in {fname}')


def update_npz(fname: str, new_data: types.ArrayDict) -> None:
    """Update values in a array dict."""
    data = load_npz(fname)
    data.update(new_data)
    save_npz(fname, data)


def load_npz(fname: str) -> types.ArrayDict:
    """Load an array dict."""
    check_ext(fname, 'npz')
    return dict(np.load(fname))


def save_npz(fname: str, data: types.ArrayDict) -> None:
    """Save an array dict."""
    check_ext(fname, 'npz')
    np.savez_compressed(fname, **data)


def save_compressed_h5ad(fname: str, adata: types.AnnData) -> None:
    """Conveniece function for compressed h5ad, make sure we import hdf5plugin."""
    check_ext(fname, 'h5ad')
    adata.write(fname,
                compression=hdf5plugin.FILTERS["zstd"],
                compression_opts=hdf5plugin.Zstd(clevel=9).filter_options
                )


def load_predictions(fname: str) -> Tuple[types.Array, types.Array, pd.DataFrame]:
    """Predictions have a special header + index, so we have a helper function."""
    check_ext(fname, 'csv')
    pred_df = pd.read_csv(fname, index_col=[0], header=[0, 1, 2])
    assert pred_df.columns[0][0] == 'batch', 'Expected batch values in first column'
    batch = pred_df[pred_df.columns[0]].values
    assert pred_df.columns[1][0] == 'perturbed', 'Expected perturbed values in first column'
    perturbed = pred_df[pred_df.columns[1]].values
    # Trim true values.
    pred_df = pred_df[pred_df.columns[2:]]
    return batch, perturbed, pred_df
