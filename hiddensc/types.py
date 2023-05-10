from typing import Sequence, Dict

import anndata
import numpy as np

# An annotated data matrix.
AnnData = anndata._core.anndata.AnnData
Array = np.ndarray
StrArray = Sequence[str]
ArrayDict = Dict[str, np.ndarray]
NumberDict = Dict[str, float]
StrArrayDict = Dict[str, StrArray]
