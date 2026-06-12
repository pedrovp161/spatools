import os
import warnings

import numpy as np
import scanpy as sc
from anndata import AnnData, ImplicitModificationWarning


# ── Outlier detection ──────────────────────────────────────────────────────────

_VALID_METHODS = {"high", "low", "both"}


def is_outlier(x, k: float = 4, method: str = "both") -> np.ndarray:
    """
    Identify outliers using Median Absolute Deviation (MAD).

    Parameters
    ----------
    x : array-like
        Input values.
    k : float
        Number of MADs from the median to consider a value an outlier.
    method : {'high', 'low', 'both'}
        Direction of outlier detection.

    Returns
    -------
    np.ndarray
        Boolean mask — True where a value is an outlier.

    Examples
    --------
    >>> is_outlier([1, 2, 3, 4, 5, 100], method="high")
    array([False, False, False, False, False,  True])
    """
    if method not in _VALID_METHODS:
        raise ValueError(f"method must be one of {_VALID_METHODS}, got {method!r}")

    x = np.asarray(x, dtype=float)
    median = np.median(x)
    deviation = x - median
    mad = np.median(np.abs(deviation))

    if mad == 0:
        return np.zeros(len(x), dtype=bool)

    threshold = k * mad
    conditions = {
        "high": deviation > threshold,
        "low":  deviation < -threshold,
        "both": np.abs(deviation) > threshold,
    }
    return conditions[method]


# ── I/O ───────────────────────────────────────────────────────────────────────

def _normalize_output_path(output_dir: str, name: str) -> str:
    """Return the full .h5ad path for a given sample name."""
    filename = name if name.endswith(".h5ad") else f"{name}.h5ad"
    return os.path.join(output_dir, filename)


def save_adatas(output_dir: str, adatas_dict: dict[str, AnnData]) -> None:
    """
    Write each AnnData in *adatas_dict* to *output_dir* as .h5ad files.

    The directory is created if it does not exist. The caller is responsible
    for validating that *output_dir* is appropriate for the context.

    Parameters
    ----------
    output_dir : str
        Destination directory.
    adatas_dict : dict[str, AnnData]
        Mapping of sample ID → AnnData.
    """
    os.makedirs(output_dir, exist_ok=True)

    warnings.filterwarnings(
        "ignore",
        message="Trying to modify attribute",
        category=ImplicitModificationWarning,
    )

    for name, adata in adatas_dict.items():
        path = _normalize_output_path(output_dir, name)
        sc.write(path, adata)  # type: ignore[arg-type]