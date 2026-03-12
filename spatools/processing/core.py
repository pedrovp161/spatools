import scanpy as sc
import os

from .utils import save_spatial_files, is_outlier

def preprocessing(adatas_dict: dict, 
                 output_dir: str = "",
                 save_files: bool = False, 
                 genes_outliers: bool = False, 
                 counts_outliers: bool = False,
                 mt_percentage_outliers: bool = True,
                 genes_and_counts_outliers: bool = True
                 ):
    """
    Preprocess Visium data and store unique stats in .uns for later integration.
    """
    
    # Validações iniciais
    if (genes_and_counts_outliers and (genes_outliers or counts_outliers)) or (genes_outliers and counts_outliers):
        print("Error: Redundant outlier filters selected. Operation terminated.")
        return

    if save_files and not output_dir:
        print("Error: output_dir must be defined to save files.")
        return
    
    for i, adata in adatas_dict.items():
        # Setup
        stats = {
            "sample_id": i,
            "initial_n_spots": adata.n_obs,
            "initial_n_genes": adata.n_vars
        }
        
        adata.var_names_make_unique()

        # Identification o MT genes
        mt_mask = adata.var_names.str.startswith("MT-") | \
                  adata.var["gene_ids"].str.startswith("MT-", na=False)
        
        adata.var["mt"] = mt_mask
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        # Applying filters (MAD)
        if counts_outliers:
            mask = is_outlier(adata.obs['log1p_total_counts'], method="low")
            adata = adata[~mask, :].copy()
            stats["n_after_counts_filter"] = adata.n_obs
            
        if genes_outliers:
            mask = is_outlier(adata.obs['log1p_n_genes_by_counts'], method="low")
            adata = adata[~mask, :].copy()
            stats["n_after_genes_filter"] = adata.n_obs

        if genes_and_counts_outliers:
            out_c = is_outlier(adata.obs['log1p_total_counts'], method="low")
            out_g = is_outlier(adata.obs['log1p_n_genes_by_counts'], method="low")
            mask = out_c | out_g
            adata = adata[~mask, :].copy()
            stats["n_after_combined_filter"] = adata.n_obs
        
        if mt_percentage_outliers:
            mask = is_outlier(adata.obs["pct_counts_mt"], method="high")
            adata = adata[~mask, :].copy()
            stats["n_after_mt_filter"] = adata.n_obs

        sc.pp.filter_genes(adata, min_cells=1)
        
        # Adding final stats
        stats["final_n_spots"] = adata.n_obs
        stats["final_n_genes"] = adata.n_vars

        adata.uns[f"preprocessing_stats_{i}"] = stats
        
        adatas_dict[i] = adata

    # Saving files if needed
    if save_files:
        os.makedirs(output_dir, exist_ok=True)
        save_spatial_files(output_dir, adatas_dict)

    return adatas_dict

