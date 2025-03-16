import os
import numpy as np
import scanpy as sc
import pandas as pd
from copy import deepcopy
from anndata import AnnData
import matplotlib.pyplot as plt
from collections import namedtuple
import matplotlib.patches as _mpatches
from matplotlib.lines import Line2D as _Line2D
from matplotlib.colors import LinearSegmentedColormap


### deprecated
def plot_bar_by_batch(adata: AnnData, clusters_col: str) -> None:
    # Verificar se cluster_col está em adata.obs
    if clusters_col not in adata.obs.columns:
        raise ValueError(f"A coluna '{clusters_col}' não está em adata.obs")
    
    # Verificar se as cores estão definidas em adata.uns
    color_key = f"{clusters_col}_colors"
    if color_key not in adata.uns:
        raise ValueError(f"As cores para '{clusters_col}' não estão definidas em adata.uns['{color_key}']")
    
    # Obter as cores dos clusters
    cluster_colors = adata.uns[color_key]

    # Agrupar dados por batch e clusters_col
    count_data = adata.obs.groupby(['batch', clusters_col]).size().unstack(fill_value=0)

    # Calcular a porcentagem
    percentage_data = count_data.div(count_data.sum(axis=1), axis=0) * 100#type: ignore

    # Definir cores usando a paleta de cores do AnnData
    cluster_labels = percentage_data.columns
    colors = [cluster_colors[int(label)] for label in cluster_labels]

    # Plotar gráfico de barras empilhadas
    ax = percentage_data.plot(kind='bar', stacked=True, figsize=(12, 6), color=colors)
    ax.set_title(f'Porcentagem de clusters em {clusters_col} para cada amostra', fontsize=25)
    ax.set_xlabel('Amostras', fontsize=25)
    ax.set_ylabel('Porcentagem (%)', fontsize=25)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=14)
    ax.legend(title="clusters", ncol=2, loc="right", bbox_to_anchor=(1.17, 0.5))
    plt.tight_layout()
    plt.show()


def plot_bar_by_group(adata: AnnData, clusters_col: str = "leiden_0.5") -> None:
    # Verificar se a coluna 'response' está em adata.obs
    if 'response' not in adata.obs.columns:
        raise ValueError("A coluna 'response' não está em bdata.obs")
    
    # Agrupar dados por response e cluster_col
    count_data = adata.obs.groupby(['response', clusters_col]).size().unstack(fill_value=0)

    # Reordenar as colunas na ordem GR, PR, BR
    count_data = count_data.reindex(['GOR', 'PAR', 'POR'])

    # Calcular a porcentagem
    percentage_data = count_data.div(count_data.sum(axis=1), axis=0) * 100#type: ignore

    # Definir as cores para os batches (opcional: pode ser ajustado conforme necessário)
    batch_colors = adata.uns[clusters_col + "_colors"]
    batch_labels = percentage_data.columns
    colors = [batch_colors[i % len(batch_colors)] for i in range(len(batch_labels))]

    # Plotar gráfico de barras empilhadas
    ax = percentage_data.plot(kind='bar', stacked=True, figsize=(12, 6), color=colors)
    ax.set_title('Porcentagem clusters por tipo de resposta', fontsize=25)
    ax.set_xlabel('Tipo de resposta', fontsize=25)
    ax.set_ylabel('Porcentagem (%)', fontsize=25)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=14, rotation=0)
    ax.legend(title='Clusters', ncol=2, loc="right", bbox_to_anchor=(1.17, 0.5))
    plt.tight_layout()
    plt.show()
### deprecated

def plot_bar(
        adata: AnnData, 
        clusters_col: str, 
        group_by: str, 
        group_order: list = None,  # type: ignore
        title: str = '', 
        xlabel: str = '', 
        ylabel: str = 'Porcentagem (%)',
        use_percentage: bool = True  # Novo argumento
    ) -> None:
    """
    Function to plot stacked bar charts with grouping.

    Parameters
    ----------
    adata : AnnData
        AnnData object.
    clusters_col : str
        Name of the column with clusters in adata.obs.
    group_by : str
        Name of the column to group by (e.g., 'batch' or 'response').
    group_order : list, optional
        Custom order of groups for plotting (default: None).
    title : str, optional
        Title of the chart (default: '').
    xlabel : str, optional
        Label for the X-axis (default: '').
    ylabel : str, optional
        Label for the Y-axis (default: 'Porcentagem (%)').
    use_percentage : bool, optional
        If True, plot percentages; otherwise, plot absolute values (default: True).

    Returns
    -------
    None
        The function displays the plot.
    """

    if clusters_col not in adata.obs.columns:
        raise ValueError(f"A coluna '{clusters_col}' não está em adata.obs")
    
    color_key = f"{clusters_col}_colors"
    if color_key not in adata.uns:
        raise ValueError(f"As cores para '{clusters_col}' não estão definidas em adata.uns['{color_key}']")

    cluster_colors = adata.uns[color_key]

    # verify if group_by is in adata.obs
    if group_by not in adata.obs.columns:
        raise ValueError(f"A coluna '{group_by}' não está em adata.obs")

    # group data by group_by e clusters_col
    count_data = adata.obs.groupby([group_by, clusters_col]).size().unstack(fill_value=0)

    # reorder groups in a specific order
    if group_order:
        count_data = count_data.reindex(group_order)

    # decide between percentage and absolute values
    if use_percentage:
        data_to_plot = count_data.div(count_data.sum(axis=1), axis=0) * 100  # type: ignore
        ylabel = 'Porcentagem (%)'
    else:
        data_to_plot = count_data
        ylabel = 'Valores absolutos' 

    # use colors defined in adata.uns
    cluster_labels = data_to_plot.columns
    colors = [cluster_colors[int(label)] for label in cluster_labels]

    # plot stacked bar
    ax = data_to_plot.plot(kind='bar', stacked=True, figsize=(12, 6), color=colors)
    ax.set_title(title, fontsize=25)
    ax.set_xlabel(xlabel, fontsize=25)
    ax.set_ylabel(ylabel, fontsize=25)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=14, rotation=60)
    ax.legend(title="clusters", ncol=2, loc="right", bbox_to_anchor=(1.17, 0.5))
    plt.tight_layout()
    plt.show()

def plot_clusters_quality_violin_boxplot(
        adata: AnnData, 
        clusters_col: str = "leiden_0.5", 
        value_col: str = "pct_counts_mt", 
        figsize: tuple = (12, 8)):
    """
    Plots violin and box plots for the percentage of mitochondrial genes by cluster.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing the data.
    cluster_col : str, optional
        Name of the column containing the clusters (default: "leiden_0.5").
    value_col : str, optional
        Name of the column containing the values to be used (default: "pct_counts_mt").
    figsize : tuple, optional
        Size of the figure (default: (12, 8)).

    Returns
    -------
    None
        The function displays the plot.
    """

    # Extrair as colunas relevantes para análise
    df = adata.obs[[clusters_col, value_col]]

    # Obter a lista de clusters únicos e ordená-los
    clusters = sorted(df[clusters_col].unique().astype(int).tolist())

    # Obter as cores associadas a cada cluster
    colors = adata.uns[f"{clusters_col}_colors"]

    # Criar a figura
    fig, ax = plt.subplots(figsize=figsize)

    # Configurações dos gráficos
    violin_width = 0.8
    boxplot_width = violin_width * 0.3

    # Adicionar gráficos de violino e boxplot para cada cluster
    for i, cluster in enumerate(clusters):
        cluster_str = str(cluster)
        cluster_data = df[df[clusters_col] == cluster_str][value_col]

        # Gráfico de violino
        parts = ax.violinplot(cluster_data, positions=[i], widths=violin_width, showmeans=False, 
                              showmedians=False, showextrema=False)
        for pc in parts['bodies']:#type: ignore
            pc.set_facecolor(colors[i])  # Cor do violino
            pc.set_edgecolor('black')
            pc.set_alpha(1)

        # Gráfico de boxplot
        ax.boxplot(cluster_data, positions=[i], widths=boxplot_width, patch_artist=True,
                   boxprops=dict(facecolor='white', color='black'),
                   medianprops=dict(color='black'),
                   whiskerprops=dict(color='black'),
                   capprops=dict(color='black'),
                   flierprops=dict(markeredgecolor='black', markersize=3))

    # Adicionar título e rótulos
    if value_col == "pct_counts_mt":
        ax.set_title('Porcentagem de genes mitocondriais por cluster')
        ax.set_ylabel('Porcentagem de genes mitocondriais (%)')

    elif value_col == "total_counts":
        ax.set_title('Número de reads por cluster')
        ax.set_ylabel('Número de reads')

    elif value_col == "n_genes_by_counts":
        ax.set_title('Número de genes por cluster')
        ax.set_ylabel('Número de genes')

    ax.set_xlabel('Clusters')

    # Aplicar tamanhos de fonte após os rótulos serem definidos
    ax.title.set_fontsize(25)
    ax.xaxis.label.set_fontsize(20)
    ax.yaxis.label.set_fontsize(20)

    # Configurar ticks do eixo x
    ax.set_xticks(range(len(clusters)))
    ax.set_xticklabels(clusters, rotation=30, fontsize=18)

    # Configurar ticks do eixo y com tamanho de fonte personalizado
    ax.set_yticklabels(ax.get_yticks(), fontsize=18)

    # Ajustar layout e exibir o gráfico
    fig.tight_layout()
    plt.show()

def plot_spatial_clusters(
        adata: AnnData, 
        clusters_col: str = "leiden_0.5", 
        cols: int = 4, 
        scale_factor: int = 3000, 
        output_file: bool = False, 
        dpi: int = 1000, 
        size = 1.5, 
        include_titles: bool = True):
    
    """
    Plots spatial images for each sample into subplots.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing the data.
    clusters_col : str, optional
        Name of the column containing the clusters (default is "leiden_0.5").
    cols : int, optional
        Number of columns for the subplot (default is 4).
    scale_factor : int, optional
        Scale factor for the figure size (default is 3000).
    output_file : bool or str, optional
        File path to save the figure. If False, the figure is not saved (default is False).
    dpi : int, optional
        Figure resolution in dpi (default is 1000).
    size : float, optional
        Size of the points on the graph (default is 1.5).
    include_titles : bool, optional
        If True, sample titles will be included (default is True).

    Returns
    -------
    None
        The function displays the plot.
    """

    keynames = adata.obs["batch"].unique()

    # Mapeamento das cores dos clusters
    clusters_colors = dict(
        zip([str(i) for i in range(len(adata.uns[f"{clusters_col}_colors"]))], adata.uns[f"{clusters_col}_colors"])
    )

    # Determinar o número de subplots baseado no número de amostras
    num_samples = len(keynames)
    rows = (num_samples + cols - 1) // cols  # Calcular número de linhas

    fig, axs = plt.subplots(rows, cols, figsize=(24, 5 * rows))  # Ajustar altura da figura
    axs = axs.flatten()  # Para indexação simples

    # Iterar sobre as amostras para plotar
    for i, library in enumerate(keynames):
        ad = adata[adata.obs['batch'] == library, :].copy()  # Uso correto de pandas slicing
        sc.pl.spatial(
            ad,
            img_key="hires",
            library_id=library,
            color=f"{clusters_col}",
            size=size,
            legend_loc=None,
            show=False,
            scale_factor=scale_factor,
            frameon=False,
            palette=[
                v for k, v in clusters_colors.items() if k in ad.obs[f'{clusters_col}'].unique().tolist()],
            ax=axs[i],  # Uso correto do eixo
            title='' if not include_titles else library  # Define o título como vazio se include_titles for False
        )
        if include_titles:
            axs[i].set_title(library, fontsize=25)  # Tamanho da fonte do título ajustado para 25

    # Remover eixos não utilizados se houver menos amostras que subplots
    for j in range(i + 1, len(axs)):
        fig.delaxes(axs[j])

    plt.tight_layout()
    if output_file:
        if not os.path.exists(os.path.dirname(output_file)):#type: ignore
            os.makedirs(os.path.dirname(output_file))#type: ignore
        plt.savefig(f"{output_file}_{library}.png", format="png", dpi=dpi)
    plt.show()

def plot_single_spatial_image(
        adata: AnnData,
        clusters_col: str = "leiden_0.5",
        scale_factor: int = 3000,
        output_file=None,
        scale: int = 6,
        title: bool = True,
        size=1.5,
        dpi: int = 1000):
    """
    Plots a single spatial image for each sample in the AnnData object.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing the data.
    clusters_col : str, optional
        Name of the column containing the clusters (default: "leiden_0.5").
    scale_factor : int, optional
        Scale factor for the image (default: 3000).
    output_file : str, optional
        Path to the output file (optional).
    scale : int, optional
        Scale factor for the figure (default: 6).
    title : bool, optional
        If True, adds the sample title (default: True).
    size : float, optional
        Size of the points in the plot (default: 1.5).
    dpi : int, optional
        Resolution of the output image (default: 1000).

    Returns
    -------
    None
        The function displays the plot.
    """


    keynames = adata.obs["batch"].unique()

    # Mapeamento das cores dos clusters
    clusters_colors = dict(
        zip([str(i) for i in range(len(adata.uns[f"{clusters_col}_colors"]))], adata.uns[f"{clusters_col}_colors"])
    )

    # Iterar sobre as amostras para plotar uma por vez
    for library in keynames:
        ad = adata[adata.obs['batch'] == library, :].copy()

        # Criar uma nova figura para cada amostra com largura e altura iguais
        plt.figure(figsize=(scale * 2, scale * 2))  # Aumenta o tamanho da figura
        sc.pl.spatial(
            ad,
            img_key="hires",
            library_id=library,
            color=f"{clusters_col}",
            size=size,
            legend_loc=None,
            show=False,
            scale_factor=scale_factor,
            frameon=False,
            palette=[
                v for k, v in clusters_colors.items() if k in ad.obs[f'{clusters_col}'].unique().tolist()]
        )

        # Condição para adicionar o título
        if title:
            plt.title(library, fontsize=25)
        else:
            plt.gca().set_title('')

        if output_file:
            if not os.path.exists(os.path.dirname(output_file)):
                os.makedirs(output_file)
            plt.savefig(f"{output_file}_{library}.png", format="png", dpi=dpi)  # Ajusta a resolução
            
        plt.show()

def z_score_corr(adata: AnnData, cluster_col: str = "clusters", show = True, title: str = "Score Z das conexões"):
    matrix_size = len(adata.obs[cluster_col].unique())
    accumulation_matrix = np.zeros((matrix_size, matrix_size))

    for key in adata.uns["correlation_matrix"].keys():
        matrix = pd.DataFrame(adata.uns["correlation_matrix"][key])
        if isinstance(matrix, pd.DataFrame):
            rows, cols = matrix.index, matrix.columns
            for i in rows:
                for j in cols:
                    accumulation_matrix[i, j] += matrix.loc[i, j] # faco o calculo do acumulado de cada linha e coluna
        else:
            print("FATAL ERROR")
    # Calcular a média
    # Dividir cada valor de accumulation_matrix pelo tamanho de corr_list
    num_matrices = len(adata.uns["correlation_matrix"].keys())
    average_matrix = accumulation_matrix / num_matrices

    max = average_matrix.max().max()
    min = average_matrix.min().min()
    sum = -min + max
    x=average_matrix.max().max()/sum
    x

    # Copiar a matriz de correlação
    average_matrix = pd.DataFrame(average_matrix)
    corr_matrix = deepcopy(average_matrix)

    # Criar um cmap personalizado que tenha branco no zero
    colors = [(0, 'blue'), (1-x, 'white'), (1, 'red')]  # cores azul -> branco -> vermelho
    cmap = LinearSegmentedColormap.from_list('custom_bwr', colors)

    # Converter a matriz de correlação em valores numéricos (caso ainda não estejam)
    corr_matrix_numeric = corr_matrix.astype(float)

    # Criar uma máscara para a parte inferior da matriz (somente parte superior visível)
    mask = ~np.tril(np.ones_like(corr_matrix_numeric, dtype=bool))

    # Criar o heatmap com matplotlib
    plt.figure(figsize=(18, 12))

    # Aplicar a máscara para ocultar a parte inferior da matriz
    plt.imshow(np.ma.masked_where(mask, corr_matrix_numeric), cmap=cmap, interpolation='nearest', vmin=corr_matrix_numeric.min().min(), vmax=corr_matrix_numeric.max().max())

    # Adicionar barra de cores (colorbar)
    plt.colorbar()

    # Adicionar os valores nas células (somente na parte visível)
    for i in range(corr_matrix_numeric.shape[0]):
        for j in range(i):  # Colocar valores apenas na parte inferior
            if i != j:  # Excluir valores de clusters iguais
                plt.text(j, i, f'{corr_matrix_numeric.iloc[i, j]:.2f}', 
                        ha='center', va='center', color='black')

    # Definir os ticks nos eixos x e y
    plt.xticks(ticks=np.arange(len(corr_matrix.columns)), labels=corr_matrix.columns)
    plt.yticks(ticks=np.arange(len(corr_matrix.index)), labels=corr_matrix.index)
    plt.tick_params(axis="x", labelsize = 16)
    plt.tick_params(axis="y", labelsize = 16)

    # Adicionar título e rótulos dos eixos
    plt.title(title, fontsize=25)
    plt.xlabel('clusters eixo x', fontsize=20)
    plt.ylabel('clusters eixo y', fontsize=20)

    # Exibir o heatmap
    if show:
        plt.tight_layout()
        plt.show()

def extract_P_number(file_name: str) -> int:
    """Extracts the number after the 'P' in the file name."""
    import re
    match = re.search(r'P(\d+)', file_name)
    return int(match.group(1)) if match else float('inf')  # Inf se não encontrar # type: ignore

def sample_classifier(file_name: str, classification_dict: dict) -> str:
    """Classifies a sample based on the file name using a classification dictionary."""
    for key, value in classification_dict.items():
        if key in file_name:
            return key
    return "Unknown"

def outlier_quality(
        *,
        path_to_directory: str = r"", 
        clusters_col: str = "pct_counts_mt", 
        group_by: str = "", 
        outlier: int = 4,
        figsize: tuple = (12, 8),
        title: str = "",
        xlabel: str = "",
        ylabel: str = "",
        outlier_type = "upper",        
        legend1_pos: tuple =(0.75, 0.95),
        legend2_pos: tuple =(0.88, 0.95),
        add_line = True,
        add_outliers = True,
        metric_to_show = "median", # mean also possible
        **kwargs
    ):
    """
    Plot a violin plot of the quality metrics of all samples.

    Parameters
    ----------
    path_to_directory : str
        Path to the directory containing the .h5ad files.
    clusters_col : str
        Column name in the AnnData object containing the quality metric to be plotted.
    group_by : str
        Column name in the AnnData object containing the group information.
    outlier : int
        The number of outliers to be shown in the plot.
    figsize : tuple
        The size of the figure.
    title : str
        The title of the figure.
    xlabel : str
        The label for the x-axis.
    ylabel : str
        The label for the y-axis.
    outlier_type : str
        Whether to show the upper or lower bound of the outliers.
    legend1_pos : tuple
        The position of the first legend.
    legend2_pos : tuple
        The position of the second legend.
    add_line : bool
        Whether to add a line to the plot.
    metric_to_show : str
        Which metric to show in the line.
    **kwargs
            This argument allows you to pass a classification dictionary that associates keys (sample nomenclatures)
            to a classification group and a corresponding color. Each key in the dictionary represents an identifier 
            that will be searched for in the file names, while the values must be lists containing the group to be assigned 
            and the color to be used in the visualization.

            Example of use:
                {
                    “GOR“: [”Good”, ‘deepskyblue’],
                    “PAR“: [”Partial”, ‘Khaki’],
                    “POR“: [”Poor”, ‘coral’]
                }

            In this example, “GOR”, “PAR” and “POR” are keys that correspond to different types of samples. The lists 
            associated with each key specify the group (e.g. “Good”) and the color (e.g. “deepskyblue”) that will be used in the 
            will be used in the graphical representation.


    Returns
    -------
    None. Shows a figure.
    """
    classification_dict = {}
    colors = {}
    for key in kwargs:
        classification_dict[key] = kwargs[key][0]
        if len(kwargs[key]) > 1:
            colors[key] = kwargs[key][1]

    # Verificar se existem arquivos .h5ad dentro da pasta
    files = [os.path.join(path_to_directory, i) for i in os.listdir(path_to_directory) if i.endswith(".h5ad")]

    names = []
    for name in files:
        name = os.path.basename(name).replace(".h5ad", "")
        names.append(name)

    # Verifications needed
    if not files:
        print(f"Aviso: Nenhum arquivo .h5ad encontrado em {path_to_directory}")
        return

    if classification_dict is None:
        raise ValueError("O dicionário de classificação 'classification_dict' deve ser fornecido.")

    # read all files and process them with quality metrics provided by scanpy
    adatas = [sc.read(file) for file in files] # type: ignore

    for adata in adatas:
        if adata.var["gene_ids"].str.startswith("ENSG")[0] == True:
            adata.var["mt"] = adata.var_names.str.startswith("MT-")
        else:
            adata.var["mt"] = adata.var["gene_ids"].str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # make sure that all the gene names are unique
    adata.var_names_make_unique()

    # Create lists to store data limits of the outliers and Iterating on the data to collect information
    all_data = []
    outliers = {}
    for i, adata in enumerate(adatas):
        data = pd.DataFrame({
            clusters_col: adata.obs[clusters_col],
            'sample': names[i]
        })
        all_data.append(data)
        
        if add_outliers:
            # Calculating the value of outliers in each sample
            col = adata.obs[clusters_col]
            median = np.median(col)
            mad = np.median(np.abs(col - median))
            k = outlier
            upper_bound = median + k * mad
            lower_bound = median - k * mad
            
            # storing the lower and upper bounds in a dictionary
            outliers[names[i]] = (lower_bound, upper_bound)

    # Concatenating all the dataframes
    percentages = pd.concat(all_data, ignore_index=True)

    # Classifying each sample using the sample_classifier function
    percentages[group_by] = percentages["sample"].apply(lambda x: sample_classifier(x, classification_dict))

    # Creating a custom order
    ClassificationOrder = namedtuple('ClassificationOrder', classification_dict.keys())

    custom_order = {key: idx for idx, key in enumerate(classification_dict)}

    # Creating a new column 'order' based on custom_order
    percentages['order'] = percentages[group_by].map(custom_order)

    # Order the numbers after "P"
    percentages['number'] = percentages["sample"].apply(extract_P_number)

    # Order the dataframe based on the "P" number and the 'order' column
    sorted_names = percentages.sort_values(by=['order', 'number'])['sample'].unique()

    # Create a figure
    fig, ax = plt.subplots(figsize=figsize)

    # Atualizar a plotagem com base na ordem correta
    positions = range(len(sorted_names))
    width = 0.4  # Largura para o violino
    box_width = width * 0.3  # Largura para o boxplot, menor que a largura do violino

    # Add plots to every sample
    for i, sample in enumerate(sorted_names):
        sample_data = percentages[percentages['sample'] == sample][clusters_col]
        group = percentages[percentages['sample'] == sample][group_by].iloc[0]
        if group == "Unknown":
            print(f"{sample} não foi possível de ser identificado com nenhum dos códigos entregues!")

        # Add violinplot
        parts = ax.violinplot(sample_data, positions=[i], widths=width, showmeans=False, showmedians=False, showextrema=False)
        for pc in parts['bodies']:# type: ignore
            pc.set_facecolor(colors[group])
            pc.set_edgecolor('black')
            pc.set_alpha(1)

        # Add boxplot
        ax.boxplot(sample_data, positions=[i], widths=box_width, patch_artist=True,
                boxprops=dict(facecolor='white', color='black'),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                flierprops=dict(color='black', markeredgecolor='black', markersize=3))

    # Configure the possible titles
    if title == "":
        if clusters_col == "pct_counts_mt":
            title = 'Porcentagem de genes mitocondriais por spot em diferentes amostras'
            if ylabel == "":
                ylabel = "Porcentagem de genes mitocondriais expressos"
        elif clusters_col == "log1p_n_genes_by_counts":
            title = "Número de genes por spot em diferentes amostras"
            if ylabel == "":
                ylabel = "Log1p do número de genes"
    if xlabel == "":
        xlabel = "Amostras"

    ax.set_title(title, fontsize=20)
    ax.set_xlabel(xlabel, fontsize=25)
    ax.set_ylabel(ylabel, fontsize=18)

    # Rotate x-axis labels
    ax.set_xticks(positions)
    ax.set_xticklabels(sorted_names, rotation=60, fontsize=16)

    # Add line
    if add_line:
        if metric_to_show == "mean":
            metric = "Mean"
            sample_mean = percentages[clusters_col].mean()
            ax.axhline(y=sample_mean, linestyle='--', color="brown", label=metric)
        elif metric_to_show == "median":
            metric = "Median"
            sample_median = percentages[clusters_col].median()
            ax.axhline(y=sample_median, linestyle='--', color="brown", label=metric)

    if add_outliers:
        for tick, sample in enumerate(sorted_names):
            lower, upper = outliers[sample]  # Pegar os limites de outliers para a amostra atual
            xpos = tick  # Posição da amostra atual no gráfico

            if outlier_type == "upper":
                # Plotar a linha superior (upper bound) para o outlier
                ax.plot([xpos - 0.5, xpos + 0.5], [upper, upper], color="red")
                if tick == 0:
                    ax.plot([xpos - 0.5, xpos + 0.5], [upper, upper], color="red", label='Outlier')
            
            elif outlier_type == "lower":
                # Plotar a linha inferior (lower bound) para o outlier
                ax.plot([xpos - 0.5, xpos + 0.5], [lower, lower], color="red")
                if tick == 0:
                    ax.plot([xpos - 0.5, xpos + 0.5], [lower, lower], color="red", label='Outlier')

    # Add legends for lines and patches
    legend_patches = [
        _mpatches.Patch(color=colors[key], label=value)
        for key, value in classification_dict.items()
    ]

    # Inicializa legend_lines como uma lista vazia
    legend_lines = []

    # Define legend_lines com base nas opções de add_line e add_outliers
    if add_line and add_outliers:
        legend_lines = [
            _Line2D([0], [0], color='red', lw=2, label='Outlier'),
            _Line2D([0], [0], color='brown', lw=2, label=metric)
        ]
    elif add_line:
        legend_lines = [
            _Line2D([0], [0], color='brown', lw=2, label=metric)
        ]
    elif add_outliers:
        legend_lines = [_Line2D([0], [0], color='red', lw=2, label='Outlier')]

    legend1 = ax.legend(handles=legend_patches, title=group_by, loc='upper left', bbox_to_anchor=legend1_pos, fontsize=12, title_fontsize=12)

    # Apenas cria legend2 se legend_lines não estiver vazio
    if legend_lines:
        legend2 = ax.legend(handles=legend_lines, title="Linhas", loc='upper left', bbox_to_anchor=legend2_pos, fontsize=12, title_fontsize=12)
        ax.add_artist(legend1)

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    # example of use:
    # classification_dict: dict[str, list[str]] = {
    #     "PAR": ["Partial", "Khaki"],
    #     "GOR": ["Good", "deepskyblue"],
    #     "POR": ["Poor", "coral"]
    # }

    # # calling the function outlier_quality
    # outlier_quality(
    #     path_to_directory=r"D:\My_decon_package\My_decon_package\data\raw",
    #     # clusters_col="log1p_n_genes_by_counts",
    #     clusters_col= "log1p_n_genes_by_counts",
    #     group_by="response",
    #     outlier_type="lower",
    #     outlier=4,
    #     legend1_pos=(0.75, 0.25),
    #     legend2_pos=(0.88, 0.25),
    #     **classification_dict #type: ignore
    # )
    #  example of use:
    classification_dict: dict[str, list[str]] = {
        "GOR": ["Good", "deepskyblue"],
        "PAR": ["Partial", "Khaki"],
        "POR": ["Poor", "coral"]
    }

    # calling the function outlier_quality
    outlier_quality(
        path_to_directory=r"/home/pedrovideira/Desktop/pack_v1/data/outlier",
        # clusters_col="pct_counts_mt",
        clusters_col= "pct_counts_mt",
        group_by="response",
        outlier_type="upper",
        outlier=4,
        add_line=True,
        add_outliers=False,
        # legend1_pos=(0.75, 0.25),
        # legend2_pos=(0.88, 0.25),
        **classification_dict #type: ignore
    )
