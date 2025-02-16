from anndata._core.views import ImplicitModificationWarning
from typing import Any
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
import os


def save_spatial_files(output_dir: str, adatas_dir: dict):# TODO verificar se deveria estar aqui
    # supress irrelevant warning
    warnings.filterwarnings("ignore", message="Trying to modify attribute", category=ImplicitModificationWarning)
    # Verifica se o diretório de saída existe, se não, cria
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for name, adata in adatas_dir.items():
        if name.endswith('.h5ad'):
            output_file_path = os.path.join(output_dir, name)
        else:
            output_file_path = os.path.join(output_dir, f"{name}.h5ad")
        sc.write(output_file_path, adata) #type: ignore

# Outliers identification
def is_outlier_high(x: Any, k=4):
    """
    Identify high outliers in an array of values.

    Parameters
    ----------
    x : array-like
        Input array of values.
    k : int, optional
        The number of median absolute deviations from the median to consider a value
        an outlier. Default is 4.

    Returns
    -------
    boolean array
        Boolean array indicating which values are low outliers.

    Notes
    -----
    This function uses the Median Absolute Deviation (MAD) to identify outliers.
    The MAD is a robust measure of the spread of the data, and outliers are defined
    as values that are more than k times the MAD away from the median.
    """
    # Compute the median of the data
    median = np.median(x)
    
    # Compute the Median Absolute Deviation (MAD)
    mad = np.median(np.abs(x - median))
    
    # Compute the threshold value for outliers
    threshold = k * mad
    
    # Return boolean array indicating which values are outliers
    return np.abs(x - median) > threshold


def is_outlier_low(x: Any, k=4):
    """
    Identify low outliers in an array of values.

    Parameters
    ----------
    x : array-like
        Input array of values.
    k : int, optional
        The number of median absolute deviations from the median to consider a value
        an outlier. Default is 4.

    Returns
    -------
    boolean array
        Boolean array indicating which values are low outliers.

    Notes
    -----
    This function uses the Median Absolute Deviation (MAD) to identify outliers.
    The MAD is a robust measure of the spread of the data, and outliers are defined
    as values that are more than k times the MAD away from the median.
    """
    # Compute the median of the data
    median = np.median(x)
    
    # Compute the Median Absolute Deviation (MAD)
    mad = np.median(np.abs(x - median))
    
    # Compute the lower threshold value for outliers
    threshold_low = median - k * mad
    
    # Return boolean array indicating which values are low outliers
    return x < threshold_low

# define função para contagem de células em todas as amostras
def check_summary(dicionario: dict):
    """
    Function to calculate the total number of cells and genes in all samples.
    
    Parameters
    ----------
    dicionario : dict
        A dictionary with the Visium 10X Genomics data. The keys are the sample names and the values are the anndata objects.
    
    Returns
    -------
    tuple
        A tuple with two values. The first value is the total number of cells and the second value is the total number of genes.
    
    """
    cell_sum = 0
    gene_sum = 0
    for i in dicionario:
        adata = dicionario[i]
        cells_before_by_sample = adata.n_obs
        cell_sum += cells_before_by_sample

        genes_before_by_sample = adata.n_vars
        gene_sum += genes_before_by_sample

    return cell_sum, gene_sum

# processing steps
def preprocessar(adatas_dir: dict, 
                 save_files = False, 
                 output_dir = None,
                 genes_and_counts_outliers = True,
                 genes_outliers = False, 
                 counts_outliers = False,
                 mt_percentage_outliers = True):
    """
    Function to preprocess spatial transcriptomics data from Visium 10X Genomics.
    
    Parameters
    ----------
    adatas_dir : dict
        A dictionary with the Visium 10X Genomics data. The keys are the sample names and the values are the anndata objects.
    save_files : bool, default=False
        Whether to save the preprocessed data. If True, the function will save the data in the directory specified by output_dir.
    output_dir : str, default=None
        The directory where to save the preprocessed data. If None, the function will not save the data.
    
    Returns
    -------
    None
    
    Notes
    -----
    The function will preprocess the data by marking and removing outliers based on the total counts, genes, and percentage of mitochondrial genes. It will also filter out genes that are not expressed in at least one cell.
    """
    if (genes_and_counts_outliers and (genes_outliers or counts_outliers)) or (genes_outliers and counts_outliers):
        return "Warning: you are filtering outliers twice, which is not recommended."

    for i in adatas_dir:
        adata = adatas_dir[i]
        adata.var_names_make_unique()

        ### Correcting adata
        if adata.var["gene_ids"].str.startswith("ENSG").iloc[0] == True:
            adata.var["mt"] = adata.var_names.str.startswith("MT-")
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        elif adata.var["gene_ids"].str.startswith("ENSG").iloc[0] == False:
            adata.var["mt"] = adata.var["gene_ids"].str.startswith("MT-")
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        # Outliers cuts
        if counts_outliers:
            adata = adata.copy()
            adata.obs['out_counts_low'] = is_outlier_low(adata.obs['log1p_total_counts'])
            true_indexes = adata.obs['out_counts_low']
            true_indexes = true_indexes[~true_indexes].index
            adata = adata[true_indexes, :]
            
        if genes_outliers:
            adata = adata.copy()
            adata.obs['out_genes'] = is_outlier_low(adata.obs['log1p_n_genes_by_counts'])
            true_indexes = adata.obs['out_genes']
            true_indexes = true_indexes[~true_indexes].index
            adata = adata[true_indexes, :]

        # Outliers combinados (contagens ou genes)
        if genes_and_counts_outliers:
            adata = adata.copy()
            adata.obs['out_counts_low'] = is_outlier_low(adata.obs['log1p_total_counts'])
            adata.obs['out_genes'] = is_outlier_low(adata.obs['log1p_n_genes_by_counts'])
            adata.obs['out_combined'] = adata.obs['out_counts_low'] | adata.obs['out_genes']
            true_indexes = ~adata.obs['out_combined']
            adata = adata[true_indexes, :]
        
        if mt_percentage_outliers:
            adata = adata.copy()
            adata.obs["out_counts_mt"] = is_outlier_high(adata.obs["pct_counts_mt"])
            true_indexes = adata.obs['out_counts_mt']
            true_indexes = true_indexes[~true_indexes].index
            adata = adata[true_indexes, :]

        adata = adata.copy()
        sc.pp.filter_genes(adata, min_cells=1)# type:ignore


        adatas_dir[i] = adata

        if save_files == True:
            if output_dir != None:
                if os.path.exists(output_dir) == False:
                    os.makedirs(output_dir)
                save_spatial_files(adatas_dir=adatas_dir, output_dir=output_dir)
            elif output_dir == None:
                print("You should define output_dir!!!")
                

def processar(adatas_dir: dict):

    for i in adatas_dir:
        adata = adatas_dir[i]

        sc.pp.filter_genes(adata, min_cells=1)#type:ignore
        ### scaling data
        
        sc.experimental.pp.normalize_pearson_residuals(adata)
        sc.pp.scale(adata)#type:ignore
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)#type:ignore
        sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
        adatas_dir[i] = adata

# if __name__ == "__main__":
#     import spatools as st
#     from spatools.reading.read import Reading as rd
#     import scanpy as sc
#     from copy import deepcopy
#     import random
#     import os

#     read = rd(dir_path=r"D:\My_decon_package\redo_\data\raw_")

#     adatas_dir = read.list_dict_with_data_h5ad()
#     print(adatas_dir)

#     adatas_dir_raw = deepcopy(adatas_dir)
#     print(adatas_dir_raw)

#     # Random seed for reproducibility
#     random.seed(42)

#     st.pp.preprocessar(adatas_dir=adatas_dir, save_files=True, mt_percentage_outliers=True)

#     # Check summary of data before and after preprocessing
#     spots_raw, genes_raw = st.pp.check_summary(dicionario=adatas_dir_raw)
#     print(f"Número de celulas antes {spots_raw}, numero de genes antes {genes_raw}")

#     spots, genes = st.pp.check_summary(dicionario=adatas_dir)
#     print(f"Número de celulas depois {spots}, numero de genes depois {genes}")
#     print(1-spots/spots_raw)
