from anndata import AnnData
from typing_extensions import Dict
import scanpy as sc
import os

from .utils import save_spatial_files, is_outlier

def anndataFilters(
    adatas_dict: Dict[str, AnnData],
    genes_outliers: bool = False,
    counts_outliers: bool = False,
    mt_percentage_outliers: bool = True,
    genes_and_counts_outliers: bool = True,
    k: float = 4
) -> Dict[str, AnnData]:

    filtered_adatas = {}

    for sample_id, adata in adatas_dict.items():
        # Copy to avoid modifying original
        adata = adata.copy()

        # Setup stats
        stats = {
            "sample_id": sample_id,
            "initial_n_spots": adata.n_obs,
            "initial_n_genes": adata.n_vars
        }

        adata.var_names_make_unique()

        # MT genes (robusto)
        mt_mask = adata.var_names.str.startswith("MT-")
        if "gene_ids" in adata.var.columns:
            mt_mask |= adata.var["gene_ids"].str.startswith("MT-")

        adata.var["mt"] = mt_mask

        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        # --- Filters ---

        if counts_outliers:
            mask = is_outlier(adata.obs['log1p_total_counts'], method="low", k=k)
            adata = adata[~mask, :].copy()
            stats["n_after_counts_filter"] = adata.n_obs

        if genes_outliers:
            mask = is_outlier(adata.obs['log1p_n_genes_by_counts'], method="low", k=k)
            adata = adata[~mask, :].copy()
            stats["n_after_genes_filter"] = adata.n_obs

        if genes_and_counts_outliers:
            out_c = is_outlier(x = adata.obs['log1p_total_counts'], method="low", k=k)
            out_g = is_outlier(x = adata.obs['log1p_n_genes_by_counts'], method="low", k=k)
            mask = out_c | out_g
            adata = adata[~mask, :].copy()
            stats["n_after_combined_filter"] = adata.n_obs

        if mt_percentage_outliers:
            mask = is_outlier(x = adata.obs["pct_counts_mt"], method="high", k=k)
            adata = adata[~mask, :].copy()
            stats["n_after_mt_filter"] = adata.n_obs

        # Remove genes não expressos
        sc.pp.filter_genes(adata, min_cells=1)

        # Final stats
        stats["final_n_spots"] = adata.n_obs
        stats["final_n_genes"] = adata.n_vars

        adata.uns[f"preprocessing_stats_{sample_id}"] = stats

        # Save back
        filtered_adatas[sample_id] = adata

    return filtered_adatas

def preprocessing(adatas_dict: dict, 
                 output_dir: str = "",
                 save_files: bool = False, 
                 genes_outliers: bool = False, 
                 counts_outliers: bool = False,
                 mt_percentage_outliers: bool = True,
                 genes_and_counts_outliers: bool = True,
                 k: float = 4
                 ) -> Dict[str, AnnData]:
    """
    Preprocess Visium data and store unique stats in .uns for later integration.
    """
    
    # Validações iniciais
    if (genes_and_counts_outliers and (genes_outliers or counts_outliers)) or (genes_outliers and counts_outliers):
        raise ValueError("Error: Redundant outlier filters selected. Operation terminated.")

    if save_files and not output_dir:
        raise ValueError("Error: output_dir must be defined to save files.")
    
    adatas_dict = anndataFilters(adatas_dict = adatas_dict,
                   genes_outliers = genes_outliers,
                   counts_outliers = counts_outliers,
                   mt_percentage_outliers = mt_percentage_outliers,
                   genes_and_counts_outliers = genes_and_counts_outliers,
                   k = k)

    # Saving files if needed
    if save_files:
        os.makedirs(output_dir, exist_ok=True)
        save_spatial_files(output_dir, adatas_dict)

    return adatas_dict

