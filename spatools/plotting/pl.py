import os
import warnings
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.colors
from copy import deepcopy
import matplotlib.cm as cm
from anndata import AnnData
from scipy.stats import norm
import matplotlib.pyplot as plt
from typing import Union, Optional
from collections import namedtuple
from matplotlib.patches import Patch
import matplotlib.patches as _mpatches
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D as _Line2D
from statsmodels.stats.multitest import multipletests
from matplotlib.colors import LinearSegmentedColormap

DEFAULT_COLORS = [
    "#cd3f35", "#62b74f", "#9e5ecf", "#c54eae", "#538342",
    "#676bc6", "#db9448", "#40bbce", "#1d19e6", "#e7b4e9",
    "#4ebba0", "#ca6f87", "#615963"
]

def bar(
        adata: AnnData, 
        clusters_col: str, 
        group_by: str, 
        group_order: list = None,  # type: ignore
        title: str = '', 
        xlabel: str = '', 
        ylabel: str = 'Percentage (%)',
        use_percentage: bool = True,
        angle: int  = 90,
        legend: str = 'Clusters',
        custom_colors: Union[dict, None] = None
    ) -> None:
    
    """
    Plot a stacked bar chart showing the distribution of categories across groups.

    This function groups observations from an AnnData object by two categorical 
    variables and visualizes their distribution as a stacked bar plot. It supports 
    both percentage-based and absolute counts, and allows custom color mappings.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing observations in `adata.obs`.

    clusters_col : str
        Column name in `adata.obs` representing the categories to be stacked 
        (e.g., clusters, samples, cell types).

    group_by : str
        Column name in `adata.obs` used to define the groups on the x-axis 
        (e.g., study, batch, condition).

    group_order : list, optional
        Custom order for the groups on the x-axis. If provided, the plot will 
        follow this order. Missing groups will appear as NaN.

    title : str, optional
        Title of the plot.

    xlabel : str, optional
        Label for the x-axis.

    ylabel : str, optional
        Label for the y-axis. Defaults to 'Percentage (%)'.

    use_percentage : bool, optional
        If True, values are normalized to percentages within each group.
        If False, absolute counts are displayed.

    angle : int, optional
        Rotation angle for x-axis tick labels.

    legend : str, optional
        Title of the legend.

    custom_colors : dict, optional
        Dictionary mapping category names (from `clusters_col`) to colors.
        Example:
            {'Sample1': '#1f77b4', 'Sample2': '#ff7f0e'}

        If None, the function will attempt to use colors stored in:
            `adata.uns[f"{clusters_col}_colors"]`

        Notes:
        - The order of colors will automatically match the plotted categories.
        - This parameter is recommended when working with non-numeric categories
          (e.g., sample names, patient IDs).

    Returns
    -------
    None
        Displays the stacked bar plot.

    Raises
    ------
    ValueError
        If `clusters_col` or `group_by` are not found in `adata.obs`.

    ValueError
        If `custom_colors` is None and no color information is found in
        `adata.uns`.

    Examples
    --------
    Basic usage with percentages:

    >>> bar(adata, clusters_col="leiden", group_by="batch")

    Using absolute counts:

    >>> bar(adata, clusters_col="leiden", group_by="batch", use_percentage=False)

    Using custom colors (e.g., for samples):

    >>> colors = dict(zip(adata.obs["sample"].unique(), st.con.BATCH_COLORS))
    >>> bar(
    ...     adata,
    ...     clusters_col="sample",
    ...     group_by="study",
    ...     custom_colors=colors,
    ...     use_percentage=False
    ... )

    Notes
    -----
    - This function is a generalization of Scanpy-style stacked bar plots.
    - It ensures consistent color mapping even when category order changes.
    - Particularly useful for publication-quality figures involving multiple
      categorical annotations.
    """

    import matplotlib.pyplot as plt

    if clusters_col not in adata.obs.columns:
        raise ValueError(f"A coluna '{clusters_col}' não está em adata.obs")

    if group_by not in adata.obs.columns:
        raise ValueError(f"A coluna '{group_by}' não está em adata.obs")

    count_data = adata.obs.groupby([group_by, clusters_col]).size().unstack(fill_value=0)# type: ignore

    if group_order:
        count_data = count_data.reindex(group_order)

    if use_percentage:
        data_to_plot = count_data.div(count_data.sum(axis=1), axis=0) * 100
    else:
        data_to_plot = count_data

    cluster_labels = data_to_plot.columns
    n_clusters = len(cluster_labels)

    if custom_colors is not None:
        # dict: {label: color}
        if isinstance(custom_colors, dict):
            colors = [custom_colors[label] for label in cluster_labels]
        
        # lista: assume ordem dos clusters
        elif isinstance(custom_colors, (list, tuple)):
            if len(custom_colors) < n_clusters:
                raise ValueError("Número de cores menor que número de clusters")
            colors = list(custom_colors[:n_clusters])
        
        else:
            raise TypeError("custom_colors deve ser lista ou dict")

    else:
        color_key = f"{clusters_col}_colors"

        if color_key in adata.uns:
            cluster_colors = adata.uns[color_key]

            # funciona mesmo se label for string
            try:
                colors = [cluster_colors[int(label)] for label in cluster_labels]
            except:
                # fallback seguro (ordem)
                colors = cluster_colors[:n_clusters]

        else:
            # 🎯 fallback automático
            cmap = cm.get_cmap("tab20", n_clusters)
            colors = [cmap(i) for i in range(n_clusters)]

    ax = data_to_plot.plot(
        kind='bar',
        stacked=True,
        figsize=(12, 6),
        color=colors
    )

    ax.set_title(title, fontsize=25)
    ax.set_xlabel(xlabel, fontsize=23)
    ax.set_ylabel(ylabel, fontsize=23)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=14, rotation=angle)

    ax.legend(
        title=legend,
        ncol=2,
        loc="upper left",
        bbox_to_anchor=(1.05, 1)
    )

    plt.tight_layout()
    plt.show()

def clusters_quality_violin_boxplot(
    adata: AnnData,
    clusters_col: str = "",
    value_col=None,
    titles=None, 
    figsize: tuple = (12, 8),
    show: bool = True
) -> None:

    if not clusters_col:
        print("clusters_col is not defined")
        return

    if value_col is None:
        print("value_col is not defined")
        return

    # 🔹 Garante lista
    if isinstance(value_col, str):
        value_col = [value_col]

    n_plots = len(value_col)

    # 🔹 Subplots
    fig, axes = plt.subplots(1, n_plots, figsize=(figsize[0] * n_plots, figsize[1]))

    if n_plots == 1:
        axes = [axes]

    colors = adata.uns[f"{clusters_col}_colors"]

    for ax, val in zip(axes, value_col):

        df = adata.obs[[clusters_col, val]].copy()
        df[val] = pd.to_numeric(pd.Series(df[val]), errors="coerce")

        clusters = sorted(df[clusters_col].astype(int).unique().tolist())

        violin_width = 0.8
        boxplot_width = violin_width * 0.3

        for i, cluster in enumerate(clusters):
            cluster_str = str(cluster)
            cluster_data = df[df[clusters_col] == cluster_str][val].dropna() # type: ignore

            parts = ax.violinplot(
                [cluster_data],
                positions=[i],
                widths=violin_width,
                showmeans=False,
                showmedians=False,
                showextrema=False
            )

            for pc in parts['bodies']:  # type: ignore
                pc.set_facecolor(colors[i])
                pc.set_edgecolor('black')
                pc.set_alpha(1)

            ax.boxplot(
                cluster_data,
                positions=[i],
                widths=boxplot_width,
                patch_artist=True,
                boxprops=dict(facecolor='white', color='black'),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                flierprops=dict(markeredgecolor='black', markersize=3)
            )

        # 🔥 TÍTULO CUSTOMIZADO
        if titles is not None:
            if isinstance(titles, list):
                title = titles[value_col.index(val)]
            elif isinstance(titles, dict):
                title = titles.get(val, val)
            else:
                title = val
        else:
            # fallback automático
            if val == "pct_counts_mt":
                title = "Percentual of mitochondrial genes by spot"
            elif val == "total_counts":
                title = "Number of reads by spot"
            elif val == "n_genes_by_counts":
                title = "Number of genes by spot"
            elif val == "log1p_n_genes_by_counts":
                title = "log1p of the number of genes by spot"
            elif val == "log1p_total_counts":
                title = "log1p of the number of reads by spot"
            else:
                title = val
        global_mean = df[val].mean()

        ax.axhline(
            global_mean,
            color='red',
            linestyle='--',
            linewidth=2,
            label='Global mean'
        )

        ax.set_title(title)
        ax.set_ylabel(" ".join(title.split(" ")[:-2]))
        ax.set_xlabel('Clusters')

        # Estética
        ax.title.set_fontsize(20)
        ax.xaxis.label.set_fontsize(16)
        ax.yaxis.label.set_fontsize(16)

        ax.set_xticks(range(len(clusters)))
        ax.set_xticklabels(clusters, rotation=30, fontsize=14)
        ax.tick_params(axis='y', labelsize=14)

    if show:
        fig.tight_layout()
        plt.show()

def spatial_plot(
    adata: AnnData,
    group: Optional[str] = None,
    highlight: Union[str, int, list, None] = None,
    sample_key: str = "sample",
    scatter_plot: bool = True,
    ncols: int = 7,
    spot_size: int = 5,
    title_fontsize: int = 18,
    custom_colors: Union[dict, list, None] = None,
    show: bool = True,
    dpi: int = 150
) -> None:
    """
    Plot spatial data for each sample or batch using AnnData spatial coordinates and images.

    Parameters
    ----------
    adata : AnnData
        Annotated data object containing spatial coordinates in `.obsm["spatial"]`,
        metadata in `.obs`, and image data in `.uns["spatial"]`.
    group : Optional[str]
        Observation key in `adata.obs` for categorical coloring, or gene name in
        `adata.var_names` for expression-based coloring. If `None`, only the tissue
        image is shown.
    highlight : Union[str, int, list, None], optional
        Value or list of values to highlight when plotting categorical groups. Highlighted
        categories are shown with dedicated colors and other categories are shown in gray.
    sample_key : str, optional
        Column in `adata.obs` that identifies each sample or batch. Default is `"sample"`.
    scatter_plot : bool, optional
        If `True`, overlay spot points on the tissue image. If `False`, show the tissue
        image with spatial axes scaled in millimeters.
    ncols : int, optional
        Maximum number of columns in the subplot grid.
    spot_size : int, optional
        Size of the scatter points.
    title_fontsize : int, optional
        Font size for each subplot title.
    custom_colors : Union[dict, list, None], optional
        Custom color mapping for categorical groups. If a dict, keys are category labels
        and values are colors. If a list, it must match the number of categories.
    show : bool, optional
        If `True`, call `plt.show()` after plotting.
    dpi : int, optional
        Figure resolution in dots per inch.

    Example:
    >>> group = "leiden_0_5"
    >>> clusters = list(adata.obs[group].cat.categories)
    >>> adata.uns[f"{group}_colors"] = NICHES_COLORS
    >>> spatial_plot(adata=adata, group=group, sample_key="library_id", spot_size=5)
    >>> # Highlight specific niche
    >>> spatial_plot(adata=adata, group=group, highlight="Niche1", sample_key="library_id")
    >>> # Gene expression visualization
    >>> spatial_plot(adata=adata, group="ENSG00000000003", sample_key="library_id")
    """
    
    batches = adata.obs[sample_key].unique()
    n = len(batches)
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(ncols * 4, nrows * 4),
        dpi=dpi,
        constrained_layout=True if scatter_plot else False
    )
    axes = np.atleast_1d(axes).flatten()

    # --- 1. Verificação de Modo ---
    do_scatter = bool(scatter_plot and group)
    is_obs = do_scatter and group in adata.obs
    is_gene = do_scatter and group in adata.var_names

    # contínuo vs categórico
    is_continuous = False
    if is_obs:
        import pandas as pd
        if pd.api.types.is_numeric_dtype(adata.obs[group]):
            is_continuous = True
    elif is_gene:
        is_continuous = True

    if do_scatter and not is_obs and not is_gene:
        raise ValueError(
            f"'{group}' não encontrado em obs ou var_names. "
            "Para plotar apenas o tecido, use scatter_plot=False."
        )
    # --- 2. Lógica de Cores (Categorias) ---
    color_map = {}
    clusters = None

    if is_obs and not is_continuous:
        if not hasattr(adata.obs[group], "cat"):
            adata.obs[group] = adata.obs[group].astype("category")  # type: ignore

        clusters = adata.obs[group].cat.categories
        base_colors = adata.uns.get(f"{group}_colors", plt.cm.tab20.colors)  # type: ignore

        if len(base_colors) < len(clusters):
            base_colors = plt.cm.get_cmap('turbo')(np.linspace(0, 1, len(clusters)))

        color_map = {str(cl): base_colors[i] for i, cl in enumerate(clusters)}

        if custom_colors:
            if isinstance(custom_colors, dict):
                for k, v in custom_colors.items():
                    color_map[str(k)] = v
            elif isinstance(custom_colors, list) and len(custom_colors) == len(clusters):
                color_map = {str(cl): custom_colors[i] for i, cl in enumerate(clusters)}

        if highlight is not None:
            hi_cols = ["red", "yellow", "blue"]
            hi_list = [str(highlight)] if not isinstance(highlight, list) else [str(x) for x in highlight]
            color_map = {
                str(cl): (
                    hi_cols[hi_list.index(str(cl))]
                    if str(cl) in hi_list else "#D3D3D3"
                )
                for cl in clusters
            }

    # --- 3. Lógica de Escala (Contínuo: gene OU obs numérico) ---
    vmin, vmax = None, None

    if is_continuous:
        if is_gene:
            expr_global = adata[:, group].X
        else:
            expr_global = adata.obs[group].values

        if hasattr(expr_global, "toarray"):
            expr_global = expr_global.toarray()  # type: ignore

        expr_global = np.array(expr_global).flatten()
        vmin, vmax = float(np.min(expr_global)), float(np.max(expr_global))

    # --- 4. Loop de Plotagem ---
    for idx, batch in enumerate(batches):
        b = adata[adata.obs[sample_key] == batch]
        ax = axes[idx]
        
        spatial_data = b.uns["spatial"][batch]
        res_key = "hires" if "hires" in spatial_data["images"] else "lowres"
        img = spatial_data["images"][res_key]
        scale = spatial_data["scalefactors"][f"tissue_{res_key}_scalef"]
        
        # Ajuste de escala para o modo métrico
        if not scatter_plot:
            # No modo 6.5mm, a escala é relativa ao tamanho da imagem (proporcional)
            # coordenadas_mm = (coords_originais * scale_do_tecido) / tamanho_em_pixels * 6.5
            h, w = img.shape[:2]
            coords_display = (b.obsm["spatial"] * scale) / w * 6.5
            ax.imshow(img, extent=[0, 6.5, 6.5, 0], aspect="auto")
            
            ax.xaxis.set_major_locator(MultipleLocator(2))
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            ax.yaxis.set_major_locator(MultipleLocator(2))
            ax.yaxis.set_minor_locator(MultipleLocator(1))
            ax.set_xlabel("mm")
            ax.set_ylabel("mm")
        else:
            coords_display = b.obsm["spatial"] * scale
            ax.imshow(img)
            ax.axis("off")

        # Plotar pontos apenas se solicitado
        if do_scatter:
            if is_obs:
                c = [color_map[str(v)] for v in b.obs[group]]
                cmap = None
            else:
                expr = b[:, group].X
                c = expr.toarray().flatten() if hasattr(expr, "toarray") else expr.flatten()# type: ignore
                cmap = "plasma"
            
            ax.scatter(coords_display[:, 0], coords_display[:, 1], c=c, s=spot_size, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0)

        ax.set_title(str(batch), fontsize=title_fontsize)

    # Limpar eixos vazios e Adicionar Legendas
    for j in range(idx + 1, len(axes)): axes[j].axis("off")

    if is_obs:
        handles = [Patch(facecolor=color_map[str(cl)], label=str(cl)) for cl in clusters]# type: ignore
        fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1.1, 0.5), title=group)
    elif is_gene:
        sm = cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), cmap="plasma")# type: ignore
        fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.04, label=group)

    if show: plt.show()

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

def z_score_matrixplot(adata: AnnData, 
                       show=True, 
                       title: str = "Z-score of connections",
                       mask_upper=True,
                       return_object=True,
                       fontsize_ticks=16,
                       fontsize_values=12,
                       fontsize_title=24,
                       fontsize_labels=18,
                       fontsize_colorbar=14):

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    # --- Initial checks ---
    if "zscore_matrix" not in adata.uns:
        raise KeyError("'zscore_matrix' key was not found in adata.uns")
    if not isinstance(adata.uns["zscore_matrix"], dict):
        raise ValueError("'zscore_matrix' is not a dictionary")
    
    # --- Collect all unique labels ---
    all_labels = set()
    for key, mat in adata.uns["zscore_matrix"].items():
        df = pd.DataFrame(mat)
        all_labels.update(df.index)
        all_labels.update(df.columns)
    all_labels = sorted(list(all_labels))

    label_to_idx = {label: idx for idx, label in enumerate(all_labels)}

    matrix_size = len(all_labels)
    accumulation_matrix = np.zeros((matrix_size, matrix_size))

    # --- Sum matrices ---
    for key, mat in adata.uns["zscore_matrix"].items():
        matrix = pd.DataFrame(mat)
        for i in matrix.index:
            for j in matrix.columns:
                accumulation_matrix[label_to_idx[i], label_to_idx[j]] += matrix.loc[i, j]

    # --- Average ---
    num_matrices = len(adata.uns["zscore_matrix"])
    average_matrix = accumulation_matrix / num_matrices

    corr_matrix = pd.DataFrame(average_matrix, index=all_labels, columns=all_labels)

    # --- Colormap ---
    vmax = corr_matrix.values.max()
    vmin = corr_matrix.values.min()
    norm_range = vmax - vmin
    zero_pos = (0 - vmin) / norm_range if norm_range != 0 else 0.5
    colors = [(0, 'blue'), (zero_pos, 'white'), (1, 'red')]
    cmap = LinearSegmentedColormap.from_list('custom_bwr', colors)

    # --- Mask ---
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool)) if mask_upper else None
    masked_matrix = np.ma.masked_where(mask, corr_matrix) if mask is not None else corr_matrix

    # --- Plot ---
    plt.figure(figsize=(18, 12))
    im = plt.imshow(masked_matrix, cmap=cmap, interpolation='nearest', 
                    vmin=vmin, vmax=vmax)

    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize=fontsize_colorbar)

    # --- Values ---
    n = corr_matrix.shape[0]
    for i in range(n):
        for j in range(i):
            if i != j:
                plt.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}', 
                         ha='center', va='center', 
                         color='black', fontsize=fontsize_values)

    # --- Axis ---
    plt.xticks(ticks=np.arange(n), labels=list(corr_matrix.columns), rotation=45, ha='right')
    plt.yticks(ticks=np.arange(n), labels=list(corr_matrix.index))

    plt.tick_params(axis="x", labelsize=fontsize_ticks)
    plt.tick_params(axis="y", labelsize=fontsize_ticks)

    plt.title(title, fontsize=fontsize_title)
    plt.xlabel('x-axis clusters', fontsize=fontsize_labels)
    plt.ylabel('y-axis clusters', fontsize=fontsize_labels)

    if show:
        plt.tight_layout()
        plt.show()

    if return_object:
        return corr_matrix

def boxplot_cluster_correlations(adata: AnnData, 
                                 cluster_col: str = "clusters", 
                                 show=True, 
                                 title: str = "Horizontal Boxplot for niche's correlation", 
                                 subset: bool = False, 
                                 value: Union[str, int] = "",
                                 figsize: tuple[int, int] = (12, 16),
                                 title_font: int = 25,
                                 label_font: int = 18,
                                 ticks_font: int = 18,
                                 limits: Union[tuple[float, float], tuple[None, None]] = (None, None)
                                 ):
    """
    Generate a horizontal boxplot based on inter-cluster correlations (avoiding duplicate symmetric pairs).
    
    Parameters
    ----------
    adata : AnnData
        AnnData object containing the z-score correlation matrices in `adata.uns["zscore_matrix"]`.
    cluster_col : str, optional
        Name of the column that stores cluster labels (default: "clusters").
    show : bool, optional
        Whether to display the plot (default: True).
    title : str, optional
        Plot title (default: "Horizontal Boxplot of Cluster Correlations").
    subset : bool, optional
        If True, only include correlations involving a specific cluster (defined by `value`).
    value : str or int, optional
        The cluster index or label to filter by when `subset=True`.
    """

    # Verificar se adata.uns["correlation_matrix"] contém os dados esperados
    if "zscore_matrix" not in adata.uns.keys():
        raise ValueError("adata.uns não contém a chave 'zscore_matrix'.")
    if cluster_col not in adata.obs.columns:
        raise ValueError(f"adata.obs does not have the {cluster_col} column")

    # Prepare data for the boxplot 
    boxplot_data = []
    samples = adata.uns["zscore_matrix"].keys()

    for sample in samples:
        matrix = adata.uns["zscore_matrix"][sample]
        if not isinstance(matrix, pd.DataFrame):
            raise ValueError(f"The z-score matrix for sample '{sample}' is not a pandas DataFrame.")

        # Handle both string and integer indexing
        row_labels = list(matrix.index)
        col_labels = list(matrix.columns)

        for i in range(len(matrix)):
            for j in range(i + 1, len(matrix)):  # Avoid symmetric duplicates
                include = True

                if subset:
                    # Converte value para string para comparação segura
                    value_str = str(value)

                    include = (
                        str(row_labels[i]) == value_str or 
                        str(col_labels[j]) == value_str
                    )

                if include:
                    boxplot_data.append({
                        "sample_key": sample,
                        "Cluster Pair": f"{row_labels[i]}-{col_labels[j]}",
                        "Correlation": matrix.iloc[i, j]
                    })


    # Convert to df
    final_data = pd.DataFrame(boxplot_data)

    # z-score para p-valor bicaudal
    final_data["pval"] = 2 * norm.sf(np.abs(final_data["Correlation"]))

    # Correção de FDR usando Benjamini-Hochberg
    reject, pvals_corr, _, _ = multipletests(final_data["pval"], method='fdr_bh')
    final_data["FDR_pval"] = pvals_corr
    final_data["significant"] = reject

    # significativos = [final_data["FDR_pval"] < 0.05]
    # final_data["significant"] = significativos
    adata.uns["stats"] = final_data

    # Horizontal boxplot 
    plt.figure(figsize=figsize)
    sns.boxplot(x="Correlation", y="Cluster Pair", data=final_data, hue="Cluster Pair", palette="Set3", orient="h", legend=False)

    # Add shaded area for insignificant region: Correction of borrefeni
    num_comparisons = len(final_data["Cluster Pair"].unique())
    bonferroni_threshold = 1 - (0.05 / (num_comparisons))  # Ajuste do limiar

    # q_normal calculus: calculate the critical value of the normal distribution (q_normal)
    q_normal = norm.ppf(bonferroni_threshold)

    plt.axvspan(-q_normal, q_normal, color='gray', alpha=0.2, label='Não Significativo (|z| < 1.96)')#type: ignore

    #  Dark central line
    plt.axvline(x=0, color='black', linestyle='--', linewidth=2, alpha=0.5)

    # plot
    plt.title(title, fontsize=title_font)
    plt.xlabel("Valores de z-score", fontsize=label_font)
    plt.ylabel("Par de Clusters", fontsize=label_font)
    plt.xticks(fontsize=ticks_font)
    plt.yticks(fontsize=ticks_font)
    if limits != (None, None):
        plt.xlim(limits)
    plt.legend()

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
        ax=None,  # New parameter to accept an axis
        show=True,
        **kwargs
    ):
        """
    
        Function to plot outlier quality with additional customization options.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The axis to plot on. If None, a new figure and axis will be created.

        Parameters
        ----------
        path_to_directory : str
            Path to the directory containing the .h5ad files.
        clusters_col : str
            Column name in the AnnData object containing the quality metric to be plotted.
        group_by : str
            Column name in the AnnData object containing the group information.
        outlier : int
            The k value from the MAD outlier detection.
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
        # If no axis is provided, create a new figure and axis
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        classification_dict = {}
        colors = {}
        for key in kwargs:
            # Validate that kwargs[key] is a list with at least one element
            if isinstance(kwargs[key], list) and len(kwargs[key]) > 0:
                classification_dict[key] = kwargs[key][0]
                if len(kwargs[key]) > 1:
                    colors[key] = kwargs[key][1]
            else:
                raise ValueError(f"Invalid format for classification_dict[{key}]: Expected a list with at least one element.")

        # Verificar se existem arquivos .h5ad dentro da pasta
        files = sorted([os.path.join(path_to_directory, i) for i in os.listdir(path_to_directory) if i.endswith(".h5ad")])

        names = [os.path.basename(f).replace(".h5ad", "") for f in files]

        # Verifications needed
        if not files:
            print(f"Aviso: Nenhum arquivo .h5ad encontrado em {path_to_directory}")
            return

        if classification_dict is None:
            raise ValueError("O dicionário de classificação 'classification_dict' deve ser fornecido.")

        # read all files and process them with quality metrics provided by scanpy
        adatas = [sc.read(file) for file in files] # type: ignore

        for adata in adatas:
            if adata.var["gene_ids"].str.startswith("ENSG").iloc[0] == True:
                adata.var["mt"] = adata.var_names.str.startswith("MT-")
            else:
                adata.var["mt"] = adata.var["gene_ids"].str.startswith("MT-")# type:ignore
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
                mad = np.median(np.abs(col - median))#type: ignore
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
        
        # Create a combined sort key: order (0,1,2...) + number (1-31)
        # This ensures proper sorting: GOR comes first, then PAR, then POR
        # Within each group, samples are sorted by number
        percentages['sort_key'] = percentages['order'] * 1000 + percentages['number']

        # Order the dataframe based on the combined sort key and get unique samples
        sorted_names = sorted(percentages['sample'].unique(), key=lambda x: extract_P_number(x))

        # Atualizar a plotagem com base na ordem correta
        positions = range(len(sorted_names))
        width = 0.4  # Largura para o violino
        box_width = width * 0.3  # Largura para o boxplot, menor que a largura do violino

        # Use existing axis if provided, otherwise ensure `fig` is defined
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        # Adicionar os gráficos para cada amostra
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
                title = 'Percentage of mitochondrial genes per spot in different samples'
                if ylabel == "":
                    ylabel = "Percentage of expressed mitochondrial genes"
            elif clusters_col == "log1p_n_genes_by_counts":
                title = "Number of genes per spot in different samples"
                if ylabel == "":
                    ylabel = "log1p of number of genes"
        if xlabel == "":
            xlabel = "Samples"

        ax.set_title(title, fontsize=20)
        ax.set_xlabel(xlabel, fontsize=25)
        ax.set_ylabel(ylabel, fontsize=18)

        # Rotate x-axis labels
        ax.set_xticks(positions)
        ax.set_xticklabels(sorted_names, rotation=90, fontsize=16)

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
                _Line2D([0], 
                        [0], 
                        color='brown', 
                        lw=2, 
                        label=metric
                        )
            ]
        elif add_outliers:
            legend_lines = [_Line2D([0], 
                                    [0], 
                                    color='red', 
                                    lw=2, 
                                    label='Outlier'
                                    )]

        legend1 = ax.legend(handles=legend_patches, title=group_by, loc='upper left', bbox_to_anchor=legend1_pos, fontsize=12, title_fontsize=12)

        # Apenas cria legend2 se legend_lines não estiver vazio
        if legend_lines:
            legend2 = ax.legend(handles=legend_lines, 
                                title="Linhas", 
                                loc='upper left', 
                                bbox_to_anchor=legend2_pos, 
                                fontsize=12, 
                                title_fontsize=12
                                )
            
            ax.add_artist(legend1)

        fig.tight_layout()
        if show:
            plt.show()

def corr_spearman(adata: AnnData, 
                  corr_matrix, 
                  pval_matrix
                  ):
    # Extrair os nomes das colunas sem cortar errado
    obsm_data = adata.obsm["q05_cell_abundance_w_sf"]
    if isinstance(obsm_data, np.ndarray):
        obsm_data = pd.DataFrame(obsm_data)
    names = [col.replace("q05cell_abundance_w_sf_", "") 
            for col in obsm_data.columns]
      
    # Criar clustermap (com clustering automático de linhas e colunas)
    cg = sns.clustermap(
        corr_matrix.astype(float),
        cmap="coolwarm",
        center=0,
        figsize=(12, 10),
        xticklabels=True,
        yticklabels=True,
        row_cluster=True,   # cluster nas linhas
        col_cluster=True,    # cluster nas colunas
        cbar_pos=(0.98, 0.05, 0.1, 0.05),
        dendrogram_ratio=(0, 0)
    )

    plt.title("Correlação de Spearman entre tipos celulares", fontsize=20)

    # Salvar figura
    plt.savefig("heatmap_spearman_clustered.png", dpi=300)
    plt.close()

def preprocessing_quality_metrics(
    adatas: dict,
    title: str = "Quality Metrics: Aggregate (Top) and Individual (Bottom)",
    xlabel: str = "Filtering Stages",
    ylabel_top: str = "Total Sum",
    ylabel_bottom: str = "Spots per Sample",
    x_labels: list = ["A\n(Raw)", "B\n(Counts & Genes)", "C\n(MT Filter)"],
    legend_title_bottom: str = "Samples",
    legend_label_top: str = "Total Sum",
    figsize: tuple = (12, 9),
    save_path = None,
    font_size_base: int = 18,
    colors: Optional[list] = None,
    dodge: float = 0.02  # ← novo parâmetro
):
    """
    Quality plot for preprocessing using st.pp.Preprocessing.run
    example:
    >>> import spatools as st
    >>> adatas: dict = st.read("/path/to/your/directory")
    >>> adatas = st.pp.Preprocessing.run("MAD_combined", adatas_dict=adatas)
    >>> df = st.pl.preprocessing_quality_metrics(adatas)
    """
    
    # 1. Data Extraction
    adatas_list = list(adatas.values()) if isinstance(adatas, dict) else adatas
    
    if colors is None:
        colors = list(plt.get_cmap('tab20')(np.linspace(0, 1, len(adatas_list))))

    all_values = []
    sample_names = []

    for adata in adatas_list:
        stats_key = next((k for k in adata.uns.keys() if k.startswith("preprocessing_stats_")), None)
        if stats_key is None:
            continue
            
        stats = adata.uns[stats_key]
        values = [
            stats.get("initial_n_spots", np.nan),
            stats.get("n_after_combined_filter", np.nan),
            stats.get("n_after_mt_filter", np.nan)
        ]
        all_values.append(values)
        sample_names.append(stats.get("sample_id", "Unknown"))

    all_values = np.array(all_values)
    df = pd.DataFrame(all_values, columns=["Raw", "Counts&Genes", "MT"], index=sample_names)
    step_sums = df.sum(axis=0)

    # 2. Figure Setup
    plt.rcParams.update({'font.size': font_size_base})
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=figsize, gridspec_kw={'height_ratios': [1, 2.5]})
    fig.subplots_adjust(hspace=0.1)

    # ---------------------------------------------------------
    # Upper Axis (ax1): Aggregate Sum + % Loss
    # ---------------------------------------------------------
    ax1.plot(x_labels, step_sums, marker="D", ls="-", color="darkred", lw=3, ms=10, label=legend_label_top)

    for i, val in enumerate(step_sums):
        is_last = i == len(step_sums) - 1
        ax1.text(i, 
                 val + (val * dodge) if is_last else val - (val * dodge),
                 f"{int(val)}", ha="center", va="bottom", 
                 fontweight="bold", 
                 color="darkred", 
                 fontsize=font_size_base
                 )
        
        if i > 0:
            prev_val = step_sums[i-1]
            loss_pct = ((val - prev_val) / prev_val) * 100
            mid_x = i - 0.5
            mid_y = (val + prev_val) / 2
            ax1.text(mid_x, mid_y, f"{loss_pct:.2f}%", ha="center", va="bottom",
                     bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.2'),
                     color="red", fontweight="bold", fontsize=font_size_base)

    ax1.set_title(title, fontsize=font_size_base + 2, pad=20, fontweight="bold")
    ax1.set_ylabel(ylabel_top, fontweight="bold", fontsize=font_size_base + 1)
    ax1.legend(loc="upper right", fontsize=font_size_base)

    # ← padding para o último label não ser cortado pelo corte de eixo
    y_min1, y_max1 = ax1.get_ylim()
    ax1.set_ylim(y_min1, y_max1 * 1.08)

    # ---------------------------------------------------------
    # Lower Axis (ax2): Individual Trends
    # ---------------------------------------------------------      
    for i in range(all_values.shape[0]):
        ax2.plot(x_labels, all_values[i], marker=".", ls="-", alpha=0.6, color=colors[i], label=sample_names[i])

    ax2.set_xlabel(xlabel, fontweight="bold", fontsize=font_size_base + 1)
    ax2.set_ylabel(ylabel_bottom, fontweight="bold", fontsize=font_size_base + 1)
    
    # Legend formatting
    if len(adatas) <= 10:
        ncol = 1
    elif len(adatas) <= 20:
        ncol = 2
    elif len(adatas) <= 30:
        ncol = 3
    else:
        print("Too many data, max number of cols is 3")
        ncol = 3

    ax2.legend(title=legend_title_bottom, ncol=ncol, bbox_to_anchor=(1.02, 1), 
               loc="upper left", fontsize=font_size_base - 2, title_fontsize=font_size_base)

    # ---------------------------------------------------------
    # Aesthetics & Broken Axis
    # ---------------------------------------------------------
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.tick_params(labeltop=False, bottom=False, labelsize=font_size_base)
    ax2.tick_params(labelsize=font_size_base)

    # ← axhlines depois do set_ylim para pegar os limites corretos
    ax1.axhline(ax1.get_ylim()[0], color='black', ls=':', lw=1)
    ax2.axhline(ax2.get_ylim()[1], color='black', ls=':', lw=1)

    d = .012  
    kwargs = dict(transform=ax1.transAxes, color='black', clip_on=False, lw=2)
    ax1.plot((-d, +d), (-d, +d), **kwargs)        
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  
    kwargs.update(transform=ax2.transAxes)  
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  

    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=200)
        print(f"✅ Figure saved to: {save_path}")

    plt.show()
    return df
