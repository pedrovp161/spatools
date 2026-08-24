import os
import mygene
import itertools
import numpy as np
import pandas as pd
from pandas import Series
from pmenu_lib import pmenu
from anndata import AnnData
from time import perf_counter
from .. import constants as con
import matplotlib.pyplot as plt
from matplotlib.path import Path
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist
from typing import Optional, Union, Final, Tuple

from ..reading import read


def spatools_check(adata):
    if "spatools" in adata.uns:
        print("Overwriting old analysis!")
        i = input("Do you want to proceed? [y or n] ").strip().lower()
        if i == "n":
            raise Exception("Operation canceled by the user.")
        elif i != "y":
            raise Exception("Invalid response. Use 'y' for yes or 'n' for no.")

def mesure_distances(adata: AnnData, cluster_col: str):
    data = pd.DataFrame(adata.obsm["spatial"]) #type: ignore
    data[cluster_col] = adata.obs[cluster_col].values
    data.rename(columns={0: "x", 1: "y"}, inplace=True)

    x = np.array(data["x"])
    y = np.array(data["y"])
    colors = np.array(data[cluster_col])

    # Criando a matriz de pontos
    points = np.column_stack((x, y))
    dist_matrix = cdist(points, points)
    np.fill_diagonal(dist_matrix, np.inf)  # Evita distâncias zero consigo mesmo

    # Encontrando a menor distância e calculando o threshold
    min_distance = np.min(dist_matrix)
    threshold_distance = min_distance * 1.1

    # Aplicando threshold: excluindo distâncias maiores que o threshold
    mask = dist_matrix < threshold_distance

    # Criando um dataframe com os pontos mais próximos dentro do threshold
    nearest_points = []
    for i in range(len(points)):
        neighbors = np.where(mask[i])[0]  # Índices dos pontos dentro do threshold
        for j in neighbors:
            nearest_points.append([x[i], y[i],f"{x[i]}_{y[i]}" ,colors[i], x[j], y[j], colors[j], dist_matrix[i, j]])

    nearest_df = pd.DataFrame(nearest_points, columns=["x", "y","point_name" ,"color", "x_neigh", "y_neigh", "color_neigh", "distance"])
    
    return nearest_df

def check_spots_analysed(adata: AnnData,
                    batch_key: str = "batch",
                    spatools_key: str = "spatools"):
    # Inicializa o dicionário se não existir
    if "check_distances" not in adata.uns:
        adata.uns["check_distances"] = {}

    # Itera por cada batch
    if spatools_key in adata.uns:
        if "point_name" in adata.uns[spatools_key]:
            try:
                if len((adata.obs[batch_key]).unique()) != 1:
                    for i in adata.obs[batch_key].unique():
                        subset = adata.uns[spatools_key][adata.uns[spatools_key][batch_key] == i]
                        counts = subset["point_name"].value_counts()
                        adata.uns["check_distances"][i] = counts

                elif len((adata.obs[batch_key]).unique()) == 1:
                    counts = adata.uns[spatools_key]["point_name"].value_counts()
                    adata.uns["check_distances"][adata.obs[batch_key].unique()[0]] = counts
                else: print("Erro inesperado")
            except KeyError:
                counts = adata.uns[spatools_key]["point_name"].value_counts()
                adata.uns["check_distances"]["Sample"] = counts
        else:
            raise KeyError(f"key 'point_name' not found inside any of the subsets.")
    else:
        raise KeyError(f"Dict '{spatools_key}' not found inside adata.uns.") 

    # Now I will check the number of spots
    adata.uns["check_spots"] = {}
    for i in adata.uns["check_distances"]:
        adata.uns["check_spots"][i] = len(adata.uns["check_distances"][i])

    # Now I will check the number of spots analysed and compare with the number of spots in real object
    df = pd.DataFrame()
    for key, value in adata.uns["check_spots"].items():
        df.index = ["spots_analysed"] #type: ignore
        df[key] = value 
    df = df.T
    try:
        df["total_spots_anndata"] = pd.concat([df, pd.DataFrame(adata.obs[batch_key].value_counts())] ,axis=1)["count"]
    except KeyError: 
        df["total_spots_anndata"] = adata.n_obs
    df["percentage"] = df["spots_analysed"] / df["total_spots_anndata"] * 100

    adata.uns["check_spots"] = df

    return adata

def correlate_distances(adata: AnnData, 
                        is_concatenated=False, 
                        cluster_col: str = "cluster", 
                        batch_key: str = "batch"):
    """
    Calculates the distances between spatial points and stores the nearest neighbors within the threshold.

    Parameters

    adata : AnnData
    An AnnData object containing spatial coordinates in obsm["spatial"].
    is_concatenated : bool, optional
    Indicates whether the data has already been concatenated. Default is False.
    cluster_col : str, optional
    Name of the column in adata.obs containing cluster information.

    Returns

    adata : AnnData
    The AnnData object with the nearest neighbors stored in uns["spatools"] and
    the percentage of spots analyzed from the total in the AnnData object in uns["check_spots"].
    """
    
    # verifying if the analysis has already being done
    spatools_check(adata)

    if is_concatenated:
        merged_df = []
        for i in adata.obs[batch_key].unique():
            subset = adata[adata.obs[batch_key] == i].copy()
            nearest_df = mesure_distances(adata=subset, cluster_col=cluster_col)
            nearest_df[batch_key] = i  # add batch column
            try: # two ways of dealing with the same thing ~ minimize errors
                nearest_df["combination"] = nearest_df.apply(lambda row: tuple(sorted((int(row["color"]), int(row["color_neigh"])))), axis=1)
            except ValueError:
                nearest_df["combination"] = nearest_df.apply(lambda row: tuple(sorted((row["color"], row["color_neigh"]))), axis=1)
            merged_df.append(nearest_df)

        # Concatenating individual DataFrames before storing
        adata.uns["spatools"] = pd.concat(merged_df, ignore_index=True)

    else:
        if "spatools" not in adata.uns:
            nearest_df = mesure_distances(adata=adata, cluster_col=cluster_col)  # add batch column
            try: # two ways of dealing with the same thing ~ minimize errors
                nearest_df["combination"] = nearest_df.apply(lambda row: tuple(sorted((int(row["color"]), int(row["color_neigh"])))), axis=1)
            except ValueError:
                nearest_df["combination"] = nearest_df.apply(lambda row: tuple(sorted((row["color"], row["color_neigh"]))), axis=1)
            adata.uns["spatools"] = nearest_df

    # adding a df for result cheking
    adata = check_spots_analysed(adata, batch_key=batch_key, spatools_key="spatools")

    return adata

def remove_random_rows(df: pd.DataFrame, 
                       num_rows: int):
    # Check if the number of rows to remove is greater than the DataFrame's length
    if num_rows >= len(df):
        return pd.DataFrame()  # Return an empty DataFrame if all rows are to be removed

    # Randomly select rows to remove
    remove_indices = np.random.choice(df.index, size=num_rows, replace=False)

    # Remove selected rows
    df_removed = df.drop(index=remove_indices)#type: ignore

    return df_removed

def translate_anndata_genes(
    adata: AnnData,
    col: Optional[str] = None,
    species: str = "human",
    inplace: bool = True
) -> AnnData:
    """
    Convert Ensembl gene IDs to gene symbols in an AnnData object.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing gene information in `adata.var`.
    
    col : Optional[str], default=None
        Column in `adata.var` containing Ensembl IDs.
        If None, uses `adata.var.index`.
    
    species : str, default="human"
        Species used in mygene query.
    
    inplace : bool, default=True
        Whether to modify the object in place.

    Returns
    -------
    AnnData
        Modified AnnData object with:
        - `gene_symbol` column
        - updated `.var.index`
        - `gene_ids` column preserving original IDs
    """

    if not inplace:
        adata = adata.copy()

    # -------------------------
    # Validate inputs
    # -------------------------
    if col is not None and col not in adata.var.columns:
        raise ValueError(f"Column '{col}' not found in adata.var")

    # -------------------------
    # Extract Ensembl IDs
    # -------------------------
    if col is None:
        ensembl_ids = adata.var.index.astype(str).tolist()
    else:
        ensembl_ids = adata.var[col].astype(str).tolist()

    # -------------------------
    # Query mygene
    # -------------------------
    mg = mygene.MyGeneInfo()

    try:
        results = mg.querymany(
            ensembl_ids,
            scopes="ensembl.gene",
            fields="symbol",
            species=species
        )
    except Exception as e:
        raise RuntimeError(f"MyGene query failed: {e}")

    # -------------------------
    # Build mapping
    # -------------------------
    id_to_symbol = {
        r["query"]: r.get("symbol", r["query"])
        for r in results if "query" in r
    }

    # -------------------------
    # Map gene symbols
    # -------------------------
    if col is None:
        adata.var["gene_ids"] = adata.var.index.astype(str)
        adata.var["gene_symbol"] = adata.var.index.map(id_to_symbol)
    else:
        adata.var["gene_ids"]: Series[str] = adata.var[col].astype(str) # type: ignore
        adata.var["gene_symbol"] = adata.var[col].map(id_to_symbol)

    # fallback se algum não foi mapeado
    adata.var["gene_symbol"]: Series[str] = adata.var["gene_symbol"].fillna(adata.var["gene_ids"]) # type: ignore

    # -------------------------
    # Handle duplicates (CRÍTICO)
    # -------------------------
    new_index = adata.var["gene_symbol"].astype(str)

    duplicated = new_index.duplicated()
    if duplicated.any():
        new_index[duplicated] = (# type: ignore
            new_index[duplicated] + "_" + adata.var["gene_ids"][duplicated]
        )

    # -------------------------
    # Apply new index
    # -------------------------
    adata.var.index = new_index# type: ignore
    adata.var.index.name = None

    # garantir unicidade final
    adata.var_names_make_unique()

    return adata

# merge diferent clusters in same resolution
def merge_clusters(adata: AnnData, 
                   clusters_col: str, 
                   rename_dict: dict, 
                   new_clusters_col: str
                   ):
    """
    Merge clusters from different resolutions in the same AnnData object.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing the data.
    clusters_col : str
        Name of the column containing the cluster labels to be merged.
    rename_dict : dict
        Dictionary mapping old cluster names to new ones.
    new_clusters_col : str
        Name of the new column to store the merged cluster labels.

    Returns
    -------
    AnnData
        AnnData object with merged cluster labels.
    """
    if clusters_col not in adata.obs:
        raise KeyError(f"'{clusters_col}' column not found in 'adata.obs'.")

    # Replace old cluster labels with new ones
    # .astype(str) antes do replace: em coluna categórica o pandas deprecou o replace que
    # altera categorias. O resultado é idêntico, pois a coluna vira string logo abaixo.
    adata.obs[new_clusters_col] = adata.obs[clusters_col].astype(str).replace(rename_dict)

    # Ensure values are in the correct order
    unique_values = sorted(adata.obs[new_clusters_col].unique())

    # Create a mapping of old values to sequential integers
    value_mapping = {old_value: new_value for new_value, old_value in enumerate(unique_values)}

    # Map the new values back to the column
    adata.obs[new_clusters_col] = adata.obs[new_clusters_col].map(value_mapping)

    # Attempt to copy colors from the original column, if available
    if f"{clusters_col}_colors" in adata.uns:
        original_colors = adata.uns[f"{clusters_col}_colors"]

        # Map colors to the new cluster order if possible
        new_colors = [original_colors[value_mapping[old_value]] for old_value in unique_values]
        adata.uns[f"{new_clusters_col}_colors"] = new_colors
    else:
        print(f"Warning: '{clusters_col}_colors' not found in 'adata.uns'. Colors not transferred.")

    # Convert the new cluster column to integers, then back to strings
    adata.obs[new_clusters_col] = adata.obs[new_clusters_col].astype(int).astype(str)# type: ignore

    # Return the updated AnnData object
    return adata

def remove_spots(adata: AnnData, 
                 type: str, 
                 n_spots: int = 1, 
                 x_minimo: Optional[int] = None, 
                 x_maximo: Optional[int] = None, 
                 y_minimo: Optional[int] = None, 
                 y_maximo: Optional[int] = None, 
                 invert_x: bool = False, 
                 invert_y: bool = False):
    """
    Removes spots from the AnnData object based on spatial position.

    Parameters:

    - adata: AnnData
    Object containing spatial coordinates in obsm["spatial"].

    - type: str
    Type of removal ('y_max', 'y_min', 'x_max', 'x_min', 'y', 'x', 'all', 'lower', 'upper').

    - n_spots: int, optional (default: 1)
    Number of spots to remove for each criterion.

    - x_minimo: float, optional (default: None)
    Minimum value for the x coordinate.

    - x_maximo: float, optional (default: None)
    Maximum value for the x coordinate.

    - y_minimo: float, optional (default: None)
    Minimum value for the y coordinate.

    - y_maximo: float, optional (default: None)
    Maximum value for the y coordinate.

    - invert_x: bool, optional (default: False)
    If True, inverts the selection for x_minimo or x_maximo.

    - invert_y: bool, optional (default: False)
    If True, inverts the selection for y_minimo or y_maximo.

    Returns:

    - Updated AnnData with spots removed.
    """
    spatial_coords: np.ndarray = adata.obsm["spatial"]#type: ignore
    n_obs = adata.n_obs

    # Identificar os índices dos valores extremos
    sorted_indices = {
        "y_max": np.argsort(spatial_coords[:, 1])[-n_spots:],
        "y_min": np.argsort(spatial_coords[:, 1])[:n_spots],
        "x_max": np.argsort(spatial_coords[:, 0])[-n_spots:],
        "x_min": np.argsort(spatial_coords[:, 0])[:n_spots]
    }

    # Criar máscara para manter todos os spots
    mask = np.ones(n_obs, dtype=bool)

    if type in sorted_indices:
        mask[sorted_indices[type]] = False

    elif type == "y":
        mask[sorted_indices["y_max"]] = False
        mask[sorted_indices["y_min"]] = False

    elif type == "x":
        mask[sorted_indices["x_max"]] = False
        mask[sorted_indices["x_min"]] = False

    elif type == "all":
        for idx in sorted_indices.values():
            mask[idx] = False

    elif type in ("lower", "upper"):
        # Filtrar por x_minimo e x_maximo
        if x_minimo is not None and x_maximo is not None:
            x_filter = (spatial_coords[:, 0] > x_minimo) & (spatial_coords[:, 0] < x_maximo)
            if invert_x:
                x_filter = ~x_filter

        elif x_minimo is not None:
            x_filter = spatial_coords[:, 0] > x_minimo
            if invert_x:
                x_filter = spatial_coords[:, 0] < x_minimo

        elif x_maximo is not None:
            x_filter = spatial_coords[:, 0] < x_maximo
            if invert_x:
                x_filter = spatial_coords[:, 0] > x_maximo

        else:
            x_filter = np.ones(n_obs, dtype=bool)

        # Filtrar por y_minimo e y_maximo
        if y_minimo is not None and y_maximo is not None:
            y_filter = (spatial_coords[:, 1] > y_minimo) & (spatial_coords[:, 1] < y_maximo)
            if invert_y:
                y_filter = ~y_filter

        elif y_minimo is not None:
            y_filter = spatial_coords[:, 1] > y_minimo
            if invert_y:
                y_filter = spatial_coords[:, 1] < y_minimo

        elif y_maximo is not None:
            y_filter = spatial_coords[:, 1] < y_maximo
            if invert_y:
                y_filter = spatial_coords[:, 1] > y_maximo

        else:
            y_filter = np.ones(n_obs, dtype=bool)

        combined_filter = x_filter & y_filter
        filtered_coords = spatial_coords[combined_filter]
        combined_idx = np.where(combined_filter)[0]

        if len(filtered_coords) < n_spots:
            raise ValueError(f"Apenas {len(filtered_coords)} spots disponíveis após filtragem, mas {n_spots} são necessários.")

        y_max_idx_filtered = np.argsort(filtered_coords[:, 1])[-n_spots:] if type == "lower" else np.argsort(filtered_coords[:, 1])[:n_spots]
        y_max_idx = combined_idx[y_max_idx_filtered]
        mask[y_max_idx] = False

    else:
        raise ValueError(f"Tipo '{type}' inválido. Escolha entre {list(sorted_indices.keys()) + ['y', 'x', 'all', 'lower', 'upper']}.")

    return adata[mask, :].copy()

def z_score(adata: AnnData, 
            filter_column: str = "", 
            filter_value: Union[str, int] = "",
            batch_key: str = "batch"):
    # Verifica se a chave "spatools" existe em uns
    if "spatools" not in adata.uns:
        raise KeyError("A chave 'spatools' não foi encontrada em adata.uns")
    
    # Verifica se a coluna existe
    if filter_column:
        if filter_column not in adata.uns["spatools"]:
            raise ValueError(f"A coluna '{filter_column}' não existe em adata.uns['spatools']")

    # Filtra os dados, se necessário
    df = adata.uns["spatools"].copy()
    if filter_value:
        df = df[df[filter_column] == filter_value]

    merges = {}

    if batch_key not in df.columns:
        df[batch_key] = "sample"

    for i in df[batch_key].unique():
        filtro_batch = df[df[batch_key] == i]

        # 1. Count of observations for each combination
        filtro_batch = filtro_batch[filtro_batch["color_neigh"] != filtro_batch["color"]]
        score = pd.DataFrame(filtro_batch["combination"].value_counts()).reset_index()
        score.columns = ["combination", "count"]

        # 2. Calculating the observed proportion
        score["proportion_observed"] = score["count"] / score["count"].sum()

        # 3. Counting each individual cluster
        cluster_counts = filtro_batch["color"].value_counts()

        # 4. Frequency of clusters
        cluster_frequencies = cluster_counts / cluster_counts.sum()

        # 5. Get all possible combinations of clusters
        clusters = cluster_counts.index.tolist()
        combinacoes = list(itertools.combinations(clusters, 2))

        # 6. Calculate the expected proportion for each combination
        try:
            proporcoes_esperadas = {tuple(sorted((int(c1), int(c2)))): 2 * cluster_frequencies[c1] * cluster_frequencies[c2] for c1, c2 in combinacoes}
        except:
            proporcoes_esperadas = {tuple(sorted((c1, c2))): 2 * cluster_frequencies[c1] * cluster_frequencies[c2] for c1, c2 in combinacoes}

        # 7. convert to DataFrame
        proporcoes_esperadas_df = pd.DataFrame(list(proporcoes_esperadas.items()), columns=["combination", "proportion_expected"])

        # 8. Merge between observed and expected counts
        merged_ordered_df = pd.merge(score, proporcoes_esperadas_df, on="combination", how="outer")

        # 9. filling with 0
        merged_ordered_df.fillna(0, inplace=True)

        # x. Ajustando o número de vizinhos #### TODO REMOVED
        # average_neighbors = 6
        # total_connections = len(filtro_batch) * average_neighbors / 2 

        # xx. Contagem esperada #### TODO REMOVED
        # merged_ordered_df["expected_count"] = merged_ordered_df["proportion_expected"] * total_connections
        # merged_ordered_df["proportion_expected"] = merged_ordered_df["expected_count"] / merged_ordered_df["expected_count"].sum()

        # 10. Calculation of standard deviation
        merged_ordered_df["std_dev"] = np.sqrt((merged_ordered_df["proportion_expected"] * (1 - merged_ordered_df["proportion_expected"])) / len(filtro_batch))

        # 11. Calculating the Z-score
        merged_ordered_df["Z_score"] = (merged_ordered_df["proportion_observed"] - merged_ordered_df["proportion_expected"]) / merged_ordered_df["std_dev"]

        # 12. Adding to dictionary by batch
        merges[i] = merged_ordered_df

    # 13. Adding to adata.uns as a dictionary
    adata.uns["z-score"] = merges

    # Creating the list of correlation matrices by batch
    z_list = {}

    for batch, merged_df in merges.items():
        zscore_matrix = merged_df[["combination", 'Z_score']].copy()
        
        # Extract the tuple values for two separate columns (a, b)
        zscore_matrix[['a', 'b']] = pd.DataFrame(zscore_matrix['combination'].tolist(), index=zscore_matrix.index)

        # To treat combinations (a, b) and (b, a) as equivalent
        zscore_matrix['a'], zscore_matrix['b'] = np.minimum(zscore_matrix['a'], zscore_matrix['b']), np.maximum(zscore_matrix['a'], zscore_matrix['b'])

        # Creating the correlation matrix already as float and zero-filled, so that
        # unobserved pairs stay 0 without an object -> float downcast on fillna
        unique_values = sorted(set(zscore_matrix['a']).union(set(zscore_matrix['b'])))
        z_matrix = pd.DataFrame(
            0.0, index=unique_values, columns=unique_values, dtype=float
        )

        # Fill in the correlation matrix with the Z_scores
        for i in unique_values:
            for j in unique_values:
                if i <= j:
                    pair_score = zscore_matrix[((zscore_matrix['a'] == i) & (zscore_matrix['b'] == j)) | ((zscore_matrix['a'] == j) & (zscore_matrix['b'] == i))]['Z_score']
                    if not pair_score.empty:
                        z_matrix.loc[i, j] = pair_score.values[0]
                        z_matrix.loc[j, i] = pair_score.values[0]

        # Save to the correlation dictionary
        z_list[batch] = z_matrix
        
    adata.uns["zscore_matrix"] = z_list

    return adata

# Função para calcular correlação de Spearman entre dois tipos celulares
def spatial_spearman(adata: AnnData, #type: ignore
                     cell1: str, 
                     cell2: str, 
                     subset: str = "", 
                     value: str = "", 
                     saveFig: bool = False,
                     show: bool = False) -> Tuple[float, float]:
    if subset and value:
        adata: AnnData = adata[adata.obs[subset] == value]

    # Acessar diretamente as colunas corretas
    x = adata.obsm["q05_cell_abundance_w_sf"][f"q05cell_abundance_w_sf_{cell1}"].values.copy()#type: ignore
    y = adata.obsm["q05_cell_abundance_w_sf"][f"q05cell_abundance_w_sf_{cell2}"].values.copy()#type: ignore
    
    # Calcular correlação
    corr, p_value = spearmanr(x, y)
    
    # Grafico de dispersao: so vale a pena montar se for exibido ou salvo. Antes a figura
    # era construida e descartada em toda chamada, o que pesava em spearman_correlation_matrix
    # (N x N chamadas). O retorno e identico nos dois casos.
    if saveFig or show:
        plt.figure(figsize=(12, 8))
        plt.scatter(x, y, c='blue', alpha=0.5)#type: ignore
        plt.xlabel(f"Abundância de {cell1}", fontsize=18)
        plt.ylabel(f"Abundância de {cell2}", fontsize=18)
        plt.title(f"Correlação de Spearman entre {cell1} e {cell2}", fontsize=20)
        plt.grid(True)
        plt.text(0.95, 0.95, f'Correlação: {corr:.4f}', transform=plt.gca().transAxes,
                 fontsize=12, color='red', verticalalignment='top')
        plt.text(0.95, 0.90, f'Valor-p: {p_value:.4e}', transform=plt.gca().transAxes,
                 fontsize=12, color='red', verticalalignment='top')

        if saveFig:
            plt.savefig(f'spatial_spearman_{cell1}_{cell2}.png', dpi=300)
        if show:
            plt.show()

        plt.close()

    return (corr, p_value) # type: ignore
# Função para matriz de correlação
def spearman_correlation_matrix(adata: AnnData, 
                                subset: str = "", 
                                value: str = "") -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not subset or not value:
        raise ValueError("Os parâmetros 'subset' e 'value' devem ser fornecidos.")
    
    # Extrair os nomes das colunas sem cortar errado
    names: list[str] = [col.replace("q05cell_abundance_w_sf_", "") for col in adata.obsm["q05_cell_abundance_w_sf"].columns]#type: ignore

    corr_matrix = pd.DataFrame(index=names, columns=names, dtype=float)
    pval_matrix = pd.DataFrame(index=names, columns=names, dtype=float)

    for name1 in names:
        for name2 in names:
            corr, p_value = spatial_spearman(
                adata=adata,
                cell1=name1,
                cell2=name2,
                subset=subset,
                value=value,
                saveFig=False
            )
            corr_matrix.loc[name1, name2] = float(corr)
            pval_matrix.loc[name1, name2] = float(p_value)
    return corr_matrix, pval_matrix

class SelectionTool:
    """
    Interactive spatial spot selection tool for Visium/spatial transcriptomics data.
    
    This tool provides an interactive interface to manually select spots from spatial 
    transcriptomics images using a lasso selection. Selected spots can be stored in the 
    AnnData object for downstream analysis.
    
    Parameters
    ----------
    dir : str
        Path to a .h5ad file or directory containing .h5ad files. If a directory is 
        provided, a menu will be displayed to select the file.
    
    Attributes
    ----------
    adata : AnnData
        The AnnData object containing spatial data and gene expression.
    sample : str
        Name of the sample file being analyzed.
    coords : np.ndarray
        Spatial coordinates of spots, shape (n_spots, 2).
    image : np.ndarray
        High-resolution tissue image.
    scale : float
        Current scaling factor for displaying coordinates on the image.
    scale_mode : str
        Current scale mode ('hires' or 'lowres').
    selected : np.ndarray
        Boolean array indicating selected spots.
    alpha : float
        Transparency level for unselected spots (0.0 - 1.0).
    spot : int
        Size of spot markers on the plot.
    color : str
        Hex color code for unselected spots.
    
    Methods
    -------
    plot()
        Render the current state of the plot with selected and unselected spots.
    on_press(event)
        Handle mouse click to start lasso selection.
    on_move(event)
        Handle mouse movement to draw lasso path.
    on_release(event)
        Handle mouse release to finalize lasso selection.
    on_key(event)
        Handle keyboard shortcuts.
    on_scroll(event)
        Handle mouse scroll to adjust spot transparency.
    main()
        Initialize the interactive plot interface.
    run()
        Execute the interactive selection tool.
    
    Keyboard Shortcuts
    ------------------
    a : Save selected spots to adata.obs["selected_area"]
    c : Clear all selections
    q : Quit and close the plot
    up : Switch to high-resolution (hires) scale
    down : Switch to low-resolution (lowres) scale
    left : Decrease scale (zoom out)
    right : Increase scale (zoom in)
    d : Reset scale to default (0.998)
    w : Write AnnData object to disk
    z : Increase spot size
    x : Decrease spot size
    v : Cycle through available colors
    scroll up : Increase spot transparency
    scroll down : Decrease spot transparency
    
    Mouse Interactions
    ------------------
    click and drag : Draw lasso selection path
    
    Examples
    --------
    >>> import spatools as sp
    >>> tool = sp.tl.SelectionTool("/path/to/sample.h5ad")
    >>> adata_selected = tool.run()
    
    The selected spots will be added to adata.obs["selected_area"] with values 
    "Selected" or "Not selected".
    """

    def __init__(self, dir: str):
        print(os.path.basename(dir))
        self.process_dir(dir)
        self.dir: str = dir

        # Spatial data
        self.coords: np.ndarray = np.asanyarray(self.adata.obsm["spatial"])
        self.image: np.ndarray = self.adata.uns["spatial"][self.sample.replace(".h5ad", "")]["images"]["hires"]
        self.scales = self.adata.uns["spatial"][self.sample.replace(".h5ad", "")]["scalefactors"]

        self.hires_scale: float = self.scales.get("tissue_hires_scalef", 1.0)
        self.lowres_scale: float = self.scales.get("tissue_lowres_scalef", 1.0)
        self.scale: float = self.hires_scale  # default
        self.scale_mode: str = "hires"

        # States
        self.selected: np.ndarray = np.zeros(len(self.coords), dtype=bool)
        self.verts: list = []
        self.line_artist = None
        self.is_drawing: bool = False  # TODO

        # update time
        self.last_update: float = 0.0
        self.min_update_time: Final = 0.01
        self.a: int = 0

        # parameters
        self.alpha: float = 0.4
        self.spot: int = 10
        self.color: str = "lightgray"
        self.b: int = 0

    def process_dir(self, dir):
        if os.path.isfile(dir):
            if dir.endswith(".h5ad"):
                result = read(dir)
                assert isinstance(result, AnnData), "Expected a single AnnData"
                self.adata: AnnData = result
                self.sample: str = os.path.basename(dir)
        elif os.path.isdir(dir):
            adata, sample = self.menu(dir)
            if adata is None or sample is None:
                raise ValueError("No file selected from menu")
            self.adata = adata
            self.sample = sample

    def plot(self):
        """
        Render the current state of the plot with tissue image and spot markers.
        
        Displays the tissue image as background and overlays selected (red) and 
        unselected (configured color) spots with their respective sizes and transparency.
        """
        self.ax.clear()
        self.ax.imshow(self.image)

        # Aplicar escala às coordenadas originais
        scaled_coords = self.coords * self.scale

        # Spots não selecionados
        self.ax.scatter(scaled_coords[:, 0], scaled_coords[:, 1], 
                        c=self.color, s=self.spot-3, alpha=self.alpha, edgecolors='none')

        # Spots selecionados
        if np.any(self.selected):
            sel = scaled_coords[self.selected]
            self.ax.scatter(sel[:, 0], sel[:, 1], c="red", s=self.spot, edgecolors='none')

        self.ax.set_title("Lasso: Hold and drag | [a] Salvar | [c] Limpar | [w] Write anndata")
        self.fig.canvas.draw_idle()

    def on_press(self, event):
        """
        Handle mouse click event to start lasso selection.
        
        Initializes the lasso drawing path when the mouse button is pressed within 
        the plot axes.
        """
        if event.inaxes != self.ax:
            return

        self.is_drawing = True
        self.verts = [(event.xdata, event.ydata)]

        # conecta o movimento
        self.motion_cid = self.fig.canvas.mpl_connect("motion_notify_event", self.on_move)

        if self.line_artist:
            self.line_artist.remove()

        self.line_artist, = self.ax.plot([], [], color="cyan", lw=2)

    def on_move(self, event):
        """
        Handle mouse movement event during lasso drawing.
        
        Continuously updates the lasso path vertices as the mouse moves, with rate 
        limiting to prevent excessive updates.
        """
        # 1. Sai imediatamente se não estiver desenhando ou estiver fora do gráfico
        if not self.is_drawing or event.inaxes != self.ax:
            return

        # (Opcional) Contador de debug
        self.a += 1 

        # 2. Verifica o tempo apenas se o evento for válido
        current_time = perf_counter()
        if current_time - self.last_update < self.min_update_time:
            return

        # 3. Atualiza os vértices e o desenho
        self.verts.append((event.xdata, event.ydata))

        if len(self.verts) > 1:
            x, y = zip(*self.verts)
            if self.line_artist:
                self.line_artist.set_data(x, y)
                self.line_artist.figure.canvas.draw_idle()

        # 4. Atualiza o último tempo registrado
        self.last_update = current_time

    def on_release(self, event):
        """
        Handle mouse release event to finalize lasso selection.
        
        Closes the lasso path and selects all spots contained within, adding them 
        to the existing selection using OR logic.
        """
        if not self.is_drawing:
            return

        self.is_drawing = False

        # desconecta corretamente usando o ID
        if hasattr(self, "motion_cid"):
            self.fig.canvas.mpl_disconnect(self.motion_cid)

        if len(self.verts) > 3:
            self.verts.append(self.verts[0])
            path = Path(self.verts)

            scaled_coords = self.coords * self.scale

            inside = path.contains_points(scaled_coords)
            self.selected |= inside

        if self.line_artist:
            self.line_artist.remove()
            self.line_artist = None

        self.verts = []

        self.plot()
        print(self.a)
        
        # Limpeza visual
        if self.line_artist:
            self.line_artist.remove()
            self.line_artist = None
        
        self.verts = []
        self.plot()

    def on_key(self, event):
        """
        Handle keyboard events for tool control.
        
        Supports multiple keyboard shortcuts for selection management, scale adjustment,
        visualization parameters, and file operations. See class docstring for full 
        list of available shortcuts.
        """
        if event.key == "a":
            self.adata.obs["selected_area"] = self.selected
            print(f"Sucesso! {np.sum(self.selected)} spots selecionados.")
            self.adata.obs["selected_area"] = self.adata.obs["selected_area"].map({True: "Selected", False: "Not selected"})
            return self.adata
        
        elif event.key == "c":
            self.selected[:] = False
            self.plot()

        elif event.key == "q":
            plt.close(self.fig)

        elif event.key == "up":
            self.scale = self.hires_scale
            self.sacale_mode = "Hiers"
            print("Switched to HIERS")
            self.plot()

        elif event.key == "down":
            self.scale = self.lowres_scale
            self.scale_mode = "Lowres"
            print("Switched to LOWRES")
            self.plot()

        elif event.key == "left":
            self.scale -= 0.005
            print(self.scale)
            self.plot()

        elif event.key == "right":
            self.scale += 0.005
            print(self.scale)
            self.plot()

        elif event.key == "d":
            self.scale = 0.998 # default
            self.plot()

        elif event.key == "w":
            self.adata.write_h5ad(filename=os.path.join(self.dir, self.sample))

        elif event.key == "z":
            self.spot += 1
            self.plot()
            
        elif event.key == "x":
            self.spot -= 1
            self.plot()

        elif event.key == "v":
            self.color = con.COLORS_23_HEX[self.b]
            self.b += 1
            if self.b == len(con.COLORS_23_HEX):
                self.b = 0
            self.plot()

    def on_scroll(self, event):
        """
        Handle mouse scroll event to adjust spot transparency.
        
        Increases transparency (alpha) on scroll up and decreases on scroll down,
        with bounds checking to keep alpha between 0.1 and 0.9.
        """
        if 0.9 > self.alpha:
            if event.button == 'up':
                self.alpha += 0.1

        if 0.1 < self.alpha:
            if event.button == "down":
                self.alpha -= 0.1

        self.plot()

    def menu(self, directory) -> Tuple[Optional[AnnData], Optional[str]]:
        run: bool = True
        while run:
            try:
                # Lista arquivos e pastas
                options = os.listdir(directory)
            except PermissionError:
                print(f"without permission to access {directory}")
                return None, None

            sample: Optional[str] = pmenu(options)

            if not sample:
                print("Nenhuma opção selecionada. Saindo.")
                return None, None

            path = os.path.join(directory, sample)

            if not path.endswith(".h5ad"):
                raise ValueError("File must be in h5ad format")

            if os.path.isfile(path):
                # Entrar na pasta
                result = read(path)
                assert isinstance(result, AnnData), "Expected a single AnnData"
                adata: AnnData = result
                run = False
                return adata, sample
            else:
                print("Select a file and not a directory")

    def main(self):
        """
        Initialize and setup the interactive matplotlib interface.
        
        Creates the figure and axes, renders the initial plot, and connects all 
        event handlers for mouse and keyboard interactions.
        """
        self.fig, self.ax = plt.subplots(figsize=(8, 8))
        self.plot()
        # Event functions
        self.fig.canvas.mpl_connect("button_press_event", self.on_press)
        self.fig.canvas.mpl_connect("button_release_event", self.on_release)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.fig.canvas.mpl_connect("scroll_event", self.on_scroll)
        plt.show()

    def run(self) -> AnnData:
        """
        Execute the interactive selection tool.
        
        Launches the interactive interface and returns the AnnData object with 
        selected spots stored in obs["selected_area"].
        
        Returns
        -------
        AnnData
            The modified AnnData object containing selection results.
            Keyboard Shortcuts

        ------------------
        a : Save selected spots to adata.obs["selected_area"]
        c : Clear all selections
        q : Quit and close the plot
        up : Switch to high-resolution (hires) scale
        down : Switch to low-resolution (lowres) scale
        left : Decrease scale (zoom out)
        right : Increase scale (zoom in)
        d : Reset scale to default (0.998)
        w : Write AnnData object to disk
        z : Increase spot size
        x : Decrease spot size
        v : Cycle through available colors
        scroll up : Increase spot transparency
        scroll down : Decrease spot transparency
        
        Mouse Interactions
        ------------------
        click and drag : Draw lasso selection path
        
        Examples
        --------
        >>> import spatools as sp
        >>> tool = sp.tl.SelectionTool("/path/to/sample.h5ad")
        >>> adata_selected = tool.run()
        
        The selected spots will be added to adata.obs["selected_area"] with values 
        "Selected" or "Not selected".
        """
        self.main()
        return self.adata