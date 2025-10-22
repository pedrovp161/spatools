import os
import scipy
import cv2 as cv
import itertools
import numpy as np
import pandas as pd
from PIL import Image
from skimage import io
from anndata import AnnData
from pybiomart import Server
from .. import constants as con
import matplotlib.pyplot as plt
from scipy.spatial import distance
from typing import List, Any, Optional, Union
from scipy.spatial.distance import cdist
from multiprocessing import Pool, cpu_count

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
                    adata.uns["check_distances"][adata.obs[batch_key].unique()[0]]
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
    adata = check_spots_analysed(adata, batch_key="batch", spatools_key="spatools")

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

def convert_df_ens(ens: Any):
    """
    Given a list of Ensembl gene IDs, convert them to external gene names using the Ensembl BioMart API.

    Parameters
    ----------
    ens : List
        List of Ensembl gene IDs

    Returns
    -------
    df : pd.DataFrame
        A DataFrame with the Ensembl gene ID as index and the external gene name as the only column
    """
    if not isinstance(ens, list):
        raise ValueError("Values must be in list format")

    urls = [
        "http://www.ensembl.org",
        "http://useast.ensembl.org",
        "http://asia.ensembl.org"
    ]

    for url in urls:
        try:
            server = Server(host=url)

            dataset = (server.marts["ENSEMBL_MART_ENSEMBL"]
                          .datasets['hsapiens_gene_ensembl'])

            df = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])

            # Verificar se as colunas esperadas estão presentes
            if 'Gene stable ID' not in df.columns or 'Gene name' not in df.columns:
                raise KeyError(f"Expected columns not found in the dataset from {url}")

            result_dict = df.set_index('Gene stable ID')['Gene name'].to_dict()

            result = {i: result_dict.get(i, None) for i in ens}

            df = pd.DataFrame.from_dict(result, orient="index", columns=["Gene Name"])

            df = df.dropna()

            return df

        except KeyError as ke:
            print(f"KeyError with URL {url}: {ke}")
        except Exception as e:
            print(f"Error with URL {url}: {e}")

    raise RuntimeError("All Ensembl URLs failed.")

def convert_anndata_ens(adata: AnnData, 
                        clusters_col: str = "gene_symbol"):
    """
    Convert Ensembl gene IDs in AnnData object to external gene names.
    
    Parameters
    ----------
    adata : AnnData
        Anndata object containing the data.
    clusters_col : str, optional
        Name of the column to store the external gene names (default: "gene_symbol").
    
    Returns
    -------
    AnnData
        Anndata object with Ensembl gene IDs converted to external gene names.
    """
    gene_ids = adata.var.index.to_list()

    converted = convert_df_ens(gene_ids)  # type: ignore
    if converted is None or converted.empty:
        raise RuntimeError("Conversion failed; no valid mappings were returned.")
    
    adata.var[clusters_col] = converted["Gene Name"]
    df = adata.var

    # Convertendo a coluna para object temporariamente
    df[clusters_col] = df[clusters_col].astype('object')

    # Preenchendo os valores np.nan com os valores dos índices correspondentes
    df[clusters_col] = df[clusters_col].fillna(pd.Series(df.index, index=df.index))

    # Convertendo de volta para category
    df[clusters_col] = df[clusters_col].astype('category')

    adata.var = df

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
    adata.obs[new_clusters_col] = adata.obs[clusters_col].replace(rename_dict)

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
    adata.obs[new_clusters_col] = adata.obs[new_clusters_col].astype(int).astype(str)

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

        # Creating the correlation matrix
        unique_values = sorted(set(zscore_matrix['a']).union(set(zscore_matrix['b'])))
        z_matrix = pd.DataFrame(index=unique_values, columns=unique_values)

        # Fill in the correlation matrix with the Z_scores
        for i in unique_values:
            for j in unique_values:
                if i <= j:  
                    z_score = zscore_matrix[((zscore_matrix['a'] == i) & (zscore_matrix['b'] == j)) | ((zscore_matrix['a'] == j) & (zscore_matrix['b'] == i))]['Z_score']
                    if not z_score.empty:
                        z_matrix.loc[i, j] = z_score.values[0]
                        z_matrix.loc[j, i] = z_score.values[0]

        # Filling values to 0
        z_matrix.fillna(0, inplace=True)

        # Save to the correlation dictionary
        z_list[batch] = z_matrix
        
    adata.uns["zscore_matrix"] = z_list

    return adata


# deprecated
def calculate_distances(args):
    """
    Calculate the distances between each pair of points within a given threshold.

    Parameters
    ----------
    args : tuple
        A tuple containing the following elements:
            - centers_colors : array-like
                A 2D array with shape (n_points, 3) containing the coordinates (x, y) and color of each point.
            - idx : int
                The index of the point for which to calculate the distances.
            - threshold_distance : float
                The maximum distance between two points to consider them close.

    Returns
    -------
    data : list
        A list of lists, where each sublist contains the coordinates (x, y) of the center point, its color, the coordinates (x, y) of a neighboring point, and the distance between the two points.
    """
    centers_colors, idx, threshold_distance = args
    x, y, color_center = centers_colors[idx]
    data = []
    for j, (x2, y2, _) in enumerate(centers_colors):
        if idx != j:
            dist = distance.euclidean((x, y), (x2, y2))
            if dist < threshold_distance:
                data.append([x, y, color_center, x2, y2, dist])
    return data

def process_image(input_image_path, 
                  output_dir: str, 
                  minDist=50, 
                  param1=50, 
                  param2=0.2, 
                  minRadius=50, 
                  maxRadius=100):
    """
    Process an input image to detect circles using Hough Transform.

    Parameters
    ----------
    input_image_path : str
        The path to the input image file.
    output_dir : str
        The directory to save the output files.
    minDist : int, default=50
        Minimum distance between detected circles.
    param1 : int, default=50
        First method-specific parameter for the Hough Transform (higher threshold).
    param2 : float, default=0.2
        Second method-specific parameter for the Hough Transform (accumulator threshold).
    minRadius : int, default=50
        Minimum circle radius to be detected.
    maxRadius : int, default=100
        Maximum circle radius to be detected.

    Returns
    -------
    output_image: png
        Image containing the detected circles outlined by lines generated with Matplotlib.
    output_excel : XLSX
        Path to the Excel file in XLSX format containing a dataframe with the following columns:
        - Center_X: X-coordinate of the center point.
        - Center_Y: Y-coordinate of the center point.
        - Center_Color: Color value of the center point.
        - Neighbor_X: X-coordinate of the neighboring point.
        - Neighbor_Y: Y-coordinate of the neighboring point.
        - Distance: Distance between the center point and the neighboring point.
        - Point_Name: Name of the point in the format "Point_X_Y".
        - Color_Code: Mapped color code from the dictionary.
        - Proximity: Categorization of the distance as 'close' or 'far'.
        - Neighbor_Cluster: Cluster of the neighboring point.
        - Combination: Tuple of sorted color codes of center and neighbor points.
    """
    # Aumentar o limite de pixels
    Image.MAX_IMAGE_PIXELS = None

    # Carregar a imagem
    image = io.imread(input_image_path)

    # Converter RGBA para RGB (ignorando o canal alfa)
    if image.shape[2] == 4:
        image_rgb = image[:, :, :3]
    else:
        image_rgb = image

    # Converter a imagem RGB para escala de cinza
    gray_image = cv.cvtColor(image_rgb, cv.COLOR_BGR2GRAY)

    # Detectar círculos usando a Transformada de Hough
    circles = cv.HoughCircles(
        gray_image,
        cv.HOUGH_GRADIENT_ALT,
        dp=1,
        minDist=minDist,
        param1=param1,
        param2=param2,
        minRadius=minRadius,
        maxRadius=maxRadius 
    )

    if circles is not None:
        circles = np.uint16(np.around(circles[0, :])).astype("int")

        # Obter a cor do centro de cada círculo e seus raios
        centers_colors = [(x, y, image_rgb[y, x]) for x, y, _ in circles]
        radii = circles[:, 2]

        # Calcular a média dos raios e definir a distância limite baseada no círculo e seus 6 vizinhos mais próximos
        mean_radii_with_neighbors = []
        for i, (x, y, r) in enumerate(circles):
            # Calcular a distância para todos os outros círculos
            distances = np.array([distance.euclidean((x, y), (x2, y2)) for (x2, y2, _) in circles if (x2, y2) != (x, y)])
            # Obter os índices dos 6 círculos mais próximos
            nearest_indices = np.argsort(distances)[:6]
            # Calcular a média dos raios desses 6 círculos mais o círculo atual
            mean_radius = np.mean(np.append(radii[nearest_indices], r))
            mean_radii_with_neighbors.append(mean_radius)

        # Definir a distância limite baseada na média dos raios com os vizinhos
        threshold_distance = 2 * np.mean(mean_radii_with_neighbors) * np.sqrt(3) * 0.9

        # Preparar argumentos para paralelização
        args = [(centers_colors, i, threshold_distance) for i in range(len(centers_colors))]

        # Usar Pool para paralelizar o cálculo das distâncias
        with Pool(cpu_count()) as pool:
            results = pool.map(calculate_distances, args)

        # Combinar os resultados
        data = [item for sublist in results for item in sublist]

        df = pd.DataFrame(data, columns=['Center_X', 'Center_Y', 'Center_Color', 'Neighbor_X', 'Neighbor_Y', 'Distance'])

        # Adicionar a coluna 'Point_Name'
        df['Point_Name'] = df.apply(lambda row: f"Point_{row['Center_X']}_{row['Center_Y']}", axis=1)

        # Função para mapear a cor do centro para o dicionário
        def map_color_to_dict(color):
            for key, value in con.COLORS_23.items():
                if tuple(color) == value:
                    return key
            return None

        # Adicionar a coluna 'Color_Code'
        df['Color_Code'] = df['Center_Color'].apply(map_color_to_dict)

        # Adicionar a coluna 'proximity'
        df['proximity'] = df['Distance'].apply(lambda d: 'close' if d < threshold_distance else 'far')

        # Criar um dicionário para mapear as coordenadas dos vizinhos para seus clusters
        neighbor_clusters = {f"{x}_{y}": map_color_to_dict(color) for x, y, color in centers_colors}

        # Adicionar a coluna 'Neighbor_Cluster'
        df['Neighbor_Cluster'] = df.apply(lambda row: neighbor_clusters.get(f"{row['Neighbor_X']}_{row['Neighbor_Y']}"), axis=1)#type:ignore

        # Adicionar a coluna 'combination'
        df['combination'] = df.apply(lambda row: tuple(sorted((row['Color_Code'], row['Neighbor_Cluster']))), axis=1)

        # Salvar o dataframe em Excel
        output_excel_path = os.path.join(output_dir, "output_data.xlsx")
        df.to_excel(output_excel_path, index=False)

        # Plotar a imagem e os círculos detectados
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.imshow(image_rgb)

        # Desenhar os círculos
        for (x, y, r) in circles:
            circle = plt.Circle((x, y), r, color='black', fill=False, linewidth=0.2)#type:ignore
            ax.add_patch(circle)

        ax.set_title('Círculos Detectados')
        plt.axis('off')

        # Salvar a imagem
        output_image_path = os.path.join(output_dir, "detected_circles.png")
        plt.savefig(output_image_path, format="png", dpi=1000)
        plt.close()

        return output_image_path, output_excel_path
    else:
        print("Nenhum círculo foi detectado.")
        return None, None
