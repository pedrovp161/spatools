from anndata import ImplicitModificationWarning
import scanpy as sc
import numpy as np
import warnings
import os

# outlier detection
def is_outlier(x, k=4, method='both') -> np.ndarray:
    """
    Identify outliers in an array of values using Median Absolute Deviation (MAD).

    Parameters
    ----------
    x : array-like
        Input array of values.
    k : int, optional
        The number of median absolute deviations from the median to consider a value
        an outlier. Default is 4.
    method : str, optional
        Method to identify outliers:
        - 'high': identify only high outliers (values above median + k*MAD)
        - 'low': identify only low outliers (values below median - k*MAD)  
        - 'both': identify both high and low outliers (default)

    Returns
    -------
    boolean array
        Boolean array indicating which values are outliers.

    Notes
    -----
    This function uses the Median Absolute Deviation (MAD) to identify outliers.
    The MAD is a robust measure of the spread of the data, and outliers are defined
    as values that are more than k times the MAD away from the median.
    
    Examples
    --------
    >>> data = [1, 2, 3, 4, 5, 100]  # 100 is a high outlier
    >>> is_outlier(data, method='high')
    array([False, False, False, False, False, True])
    
    >>> data = [-100, 2, 3, 4, 5, 6]  # -100 is a low outlier
    >>> is_outlier(data, method='low') 
    array([True, False, False, False, False, False])
    
    >>> data = [-100, 2, 3, 4, 5, 100]  # both -100 and 100 are outliers
    >>> is_outlier(data, method='both')
    array([True, False, False, False, False, True])
    """
    x = np.asarray(x)

    median = np.median(x)
    deviation = x - median
    mad = np.median(np.abs(deviation))

    if mad == 0:
        return np.zeros_like(x, dtype=bool)

    threshold = k * mad

    conditions = {
        "high": deviation > threshold,
        "low": deviation < -threshold,
        "both": np.abs(deviation) > threshold,
    }

    try:
        return conditions[method]
    except KeyError:
        raise ValueError("method must be 'high', 'low', or 'both'")

# saving files
def save_spatial_files(output_dir: str, adatas_dict: dict):# TODO verificar se deveria estar aqui
    # supress irrelevant warning
    warnings.filterwarnings("ignore", message="Trying to modify attribute", category=ImplicitModificationWarning)
    # Verifica se o diretório de saída existe, se não, cria
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for name, adata in adatas_dict.items():
        if name.endswith('.h5ad'):
            output_file_path = os.path.join(output_dir, name)
        else:
            output_file_path = os.path.join(output_dir, f"{name}.h5ad")
        sc.write(output_file_path, adata) #type: ignore

