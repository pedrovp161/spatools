import os
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
from typing import List, Any, Optional
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
    data = pd.DataFrame(adata.obsm["spatial"])
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
            nearest_points.append([x[i], y[i], colors[i], x[j], y[j], colors[j], dist_matrix[i, j]])

    nearest_df = pd.DataFrame(nearest_points, columns=["x", "y", "color", "x_neigh", "y_neigh", "color_neigh", "distance"])
    
    return nearest_df

def correlate_distances(adata: AnnData, 
                        is_concatenated=False, 
                        cluster_col: str = "cluster", 
                        batch_key: str = "batch"):
    """
    Calcula as distâncias entre pontos espaciais e armazena os vizinhos mais próximos dentro do threshold.

    Parameters
    ----------
    adata : AnnData
        Objeto AnnData contendo as coordenadas espaciais em `obsm["spatial"]`.
    is_concatenated : bool, optional
        Indica se os dados já foram concatenados. Default é False.
    cluster_col : str, optional
        Nome da coluna em `adata.obs` contendo a informação do cluster.

    Returns
    -------
    adata : AnnData
        O objeto AnnData com os vizinhos mais próximos armazenados em `uns["spatools"]`.
    """
    
    # Verificando se a análise já foi feita
    spatools_check(adata)

    if is_concatenated:
        merged_df = []
        for i in adata.obs[batch_key].unique():
            subset = adata[adata.obs[batch_key] == i].copy()
            nearest_df = mesure_distances(adata=subset, cluster_col=cluster_col)
            nearest_df[batch_key] = i  # Adiciona a coluna de batch
            nearest_df["combination"] = nearest_df.apply(lambda row: tuple(sorted((int(row["color"]), int(row["color_neigh"])))), axis=1)
            merged_df.append(nearest_df)

        # Concatenando os DataFrames individuais antes de armazenar
        adata.uns["spatools"] = pd.concat(merged_df, ignore_index=True)

    else:
        if "spatools" not in adata.uns:
            nearest_df = mesure_distances(adata=adata, cluster_col=cluster_col)
            nearest_df[batch_key] = i  # Adiciona a coluna de batch
            nearest_df["combination"] = nearest_df.apply(lambda row: tuple(sorted((int(row["color"]), int(row["color_neigh"])))), axis=1)
            adata.uns["spatools"] = nearest_df

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
    Remove spots do objeto AnnData com base na posição espacial.

    Parâmetros:
    - adata: AnnData
        Objeto contendo coordenadas espaciais em obsm["spatial"].
    - type: str
        Tipo de remoção ('y_max', 'y_min', 'x_max', 'x_min', 'y', 'x', 'all', 'lower', 'upper').
    - n_spots: int, opcional (padrão: 1)
        Número de spots a remover em cada critério.
    - x_minimo: float, opcional (padrão: None)
        Valor mínimo para a coordenada x.
    - x_maximo: float, opcional (padrão: None)
        Valor máximo para a coordenada x.
    - y_minimo: float, opcional (padrão: None)
        Valor mínimo para a coordenada y.
    - y_maximo: float, opcional (padrão: None)
        Valor máximo para a coordenada y.
    - invert_x: bool, opcional (padrão: False)
        Se True, inverte a seleção para x_minimo ou x_maximo.
    - invert_y: bool, opcional (padrão: False)
        Se True, inverte a seleção para y_minimo ou y_maximo.

    Retorna:
    - AnnData atualizado com os spots removidos.
    """
    spatial_coords = adata.obsm["spatial"]
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
            filter_column: str = None, 
            filter_value: str = None):
    # Verifica se a chave "spatools" existe em uns
    if "spatools" not in adata.uns:
        raise KeyError("A chave 'spatools' não foi encontrada em adata.uns")
    
    # Verifica se a coluna existe
    if filter_column:
        if filter_column not in adata.uns["spatools"]:
            raise ValueError(f"A coluna '{filter_column}' não existe em adata.uns['spatools']")

    # Filtra os dados, se necessário
    df = adata.uns["spatools"].copy()
    if filter_value is not None:
        df = df[df[filter_column] == filter_value]

    merges = {}

    for i in df["batch"].unique():
        filtro_batch = df[df["batch"] == i]

        # 1. Contagem de observações para cada combinação
        filtro_batch = filtro_batch[filtro_batch["color_neigh"] != filtro_batch["color"]]
        score = pd.DataFrame(filtro_batch["combination"].value_counts()).reset_index()
        score.columns = ["combination", "count"]

        # 2. Calculando a proporção observada
        score["proportion_observed"] = score["count"] / score["count"].sum()

        # 3. Contagem de cada cluster individual
        cluster_counts = filtro_batch["color"].value_counts()

        # 4. Frequência dos clusters
        cluster_frequencies = cluster_counts / cluster_counts.sum()

        # 5. Obter todas as combinações possíveis de clusters
        clusters = cluster_counts.index.tolist()
        combinacoes = list(itertools.combinations(clusters, 2))

        # 6. Calcular a proporção esperada para cada combinação
        proporcoes_esperadas = {tuple(sorted((int(c1), int(c2)))): 2 * cluster_frequencies[c1] * cluster_frequencies[c2] for c1, c2 in combinacoes}

        # 7. Converter para DataFrame
        proporcoes_df = pd.DataFrame(list(proporcoes_esperadas.items()), columns=["combination", "proportion_expected"])

        # 8. Merge entre as contagens observadas e esperadas
        merged_ordered_df = pd.merge(score, proporcoes_df, on="combination", how="inner")

        # 9. Ajustando o número de vizinhos
        average_neighbors = 6
        total_connections = len(filtro_batch) * average_neighbors / 2 

        # 10. Contagem esperada
        merged_ordered_df["expected_count"] = merged_ordered_df["proportion_expected"] * total_connections
        merged_ordered_df["proportion_expected"] = merged_ordered_df["expected_count"] / merged_ordered_df["expected_count"].sum()

        # 11. Cálculo do desvio padrão
        merged_ordered_df["std_dev"] = np.sqrt((merged_ordered_df["proportion_expected"] * (1 - merged_ordered_df["proportion_expected"])) / len(filtro_batch))

        # 12. Calculando o Z-score
        merged_ordered_df["Z_score"] = (merged_ordered_df["proportion_observed"] - merged_ordered_df["proportion_expected"]) / merged_ordered_df["std_dev"]

        # 13. Adicionando ao dicionário por batch
        merges[i] = merged_ordered_df

    # 14. Adicionando ao adata.uns como um dicionário
    adata.uns["z-score"] = merges

    # Criando a lista de matrizes de correlação por batch
    corr_list = {}

    for batch, merged_df in merges.items():
        correlation_matrix = merged_df[["combination", 'Z_score']].copy()
        
        # Extrair os valores das tuplas para duas colunas separadas (a, b)
        correlation_matrix[['a', 'b']] = pd.DataFrame(correlation_matrix['combination'].tolist(), index=correlation_matrix.index)

        # Para tratar as combinações (a, b) e (b, a) como equivalentes
        correlation_matrix['a'], correlation_matrix['b'] = np.minimum(correlation_matrix['a'], correlation_matrix['b']), np.maximum(correlation_matrix['a'], correlation_matrix['b'])

        # Criar a matriz de correlação
        unique_values = sorted(set(correlation_matrix['a']).union(set(correlation_matrix['b'])))
        corr_matrix = pd.DataFrame(index=unique_values, columns=unique_values)

        # Preencher a matriz de correlação com os Z_scores
        for i in unique_values:
            for j in unique_values:
                if i <= j:  
                    z_score = correlation_matrix[((correlation_matrix['a'] == i) & (correlation_matrix['b'] == j)) | ((correlation_matrix['a'] == j) & (correlation_matrix['b'] == i))]['Z_score']
                    if not z_score.empty:
                        corr_matrix.loc[i, j] = z_score.values[0]
                        corr_matrix.loc[j, i] = z_score.values[0]

        # Preencher valores NaN com 0
        corr_matrix.fillna(0, inplace=True)

        # Salvar no dicionário de correlações
        corr_list[batch] = corr_matrix

    # Adicionando ao adata.uns
    adata.uns["correlation_matrix"] = corr_list

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
