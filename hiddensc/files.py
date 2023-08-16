from typing import Tuple

import hdf5plugin
import numpy as np
import pandas as pd
import os
import glob

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

def _read_prediction_csv(fname:str) -> pd.DataFrame:
    """Predictions have a special header + index, so we have a helper function."""
    check_ext(fname, 'csv')
    return pd.read_csv(fname, index_col=[0], header=[0, 1, 2])

def load_predictions(dirname: str) -> Tuple[types.Array, types.Array, pd.DataFrame]:
    """Load all predictions in a folder."""
    pattern = os.path.join(dirname, '*predictions.csv')
    dfs = [_read_prediction_csv(fname) for fname in glob.glob(pattern)]
    pred_df = pd.concat(dfs, axis=1)
    cols = pred_df.columns.tolist()
    values = {}
    for label in ['batch', 'perturbed']:
        is_col = np.array([label in c[0] for c in cols])
        assert np.sum(is_col) == 1, f'Expected at least one {label} column, found {np.sum(is_col)}'
        index = np.argmax(is_col)
        values[label] = pred_df[cols[index]].values
        pred_df.drop(cols[index], axis=1, inplace=True)
        
    return values['batch'], values['perturbed'], pred_df

