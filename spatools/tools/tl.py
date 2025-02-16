import os
import cv2 as cv
import numpy as np
import pandas as pd
from anndata import AnnData
from PIL import Image
from skimage import io
from .. import constants as con
from typing import List, Any
from pybiomart import Server
import matplotlib.pyplot as plt
from scipy.spatial import distance
from multiprocessing import Pool, cpu_count

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

def Z_score_simple(df):# TODO
    ...
    # df = df[df["Neighbor_Cluster"] != df["Color_Code"]]
    
    # score = pd.DataFrame(df["combination"].value_counts()).reset_index()
    # score.columns = ["combination", "count"]

    # # 2. Calculando a proporção observada
    # score["proportion_observed"] = score["count"] / score["count"].sum()

    # # 3. Contagem de cada cluster individual
    # cluster_counts = df["Color_Code"].value_counts()

    # # 4. Frequência dos clusters
    # cluster_frequencies = cluster_counts / cluster_counts.sum()

    # # 5. Obter todas as combinações possíveis de clusters (excluindo pares iguais)
    # clusters = cluster_counts.index.tolist()
    # # Usando itertools.combinations (não combinations_with_replacement) para garantir que (c1, c2) != (c2, c1)
    # combinacoes = list(itertools.combinations(clusters, 2))

    # # 6. Calcular a proporção esperada para cada combinação
    # proporcoes_esperadas = {}
    # for c1, c2 in combinacoes:
    #     chave = tuple(sorted((c1, c2)))
    #     proporcoes_esperadas[chave] = 2 * cluster_frequencies[c1] * cluster_frequencies[c2]  # Sem pares iguais

    # # 7. Converter as combinações e proporções esperadas para um DataFrame
    # proporcoes_df = pd.DataFrame(list(proporcoes_esperadas.items()), columns=["combination", "proportion_expected"])

    # # 8. Merge entre as contagens observadas e esperadas
    # merged_ordered_df = pd.merge(score, proporcoes_df, on="combination", how="inner")

    # # 9. Ajustando o número de vizinhos para o cálculo (assumindo 6 vizinhos em média por bola)
    # average_neighbors = 6
    # total_connections = len(df) * average_neighbors / 2  # Dividido por 2 para não contar conexões duplicadas

    # # 10. Calcular a contagem esperada para cada combinação de acordo com o total de conexões
    # merged_ordered_df["expected_count"] = merged_ordered_df["proportion_expected"] * total_connections
    # merged_ordered_df["proportion_expected"] = merged_ordered_df["expected_count"] / merged_ordered_df["expected_count"].sum()

    # # 11. Calculando o desvio padrão da contagem esperada usando a variável binomial
    # merged_ordered_df["std_dev"] = np.sqrt((merged_ordered_df["proportion_expected"] * (1 - merged_ordered_df["proportion_expected"])) / len(df))

    # # 12. Calculando o Z-score
    # merged_ordered_df["Z_score"] = (merged_ordered_df["proportion_observed"] - merged_ordered_df["proportion_expected"]) / merged_ordered_df["std_dev"]

    # # 13. todos os df em uma lista
    # merges.append(merged_ordered_df)

