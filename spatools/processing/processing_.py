import scanpy as sc

class Processing:

    pipelines = {}

    @classmethod
    def register(cls, name):
        def decorator(func):
            cls.pipelines[name] = func
            return func
        return decorator

    @classmethod
    def run(cls, name, adata, **kwargs):
        """
        Execute a registered downstream processing pipeline.

        This method acts as a central dispatcher for various single-cell or 
        spatial transcriptomics processing workflows (e.g., Pearson residuals, 
        standard log-normalization).

        Parameters
        ----------
        name : str
            The identifier of the pipeline to execute. Default registered 
            options include:
            - 'pearson_PCA': Normalizes using Pearson residuals, followed by 
              PCA, Neighbors, UMAP, and Leiden clustering.
            - 'lognorm': Performs total count normalization, log1p transformation, 
              HVG selection, scaling, and PCA.
            - 'scVI_paper': Placeholder for scVI-based integration/processing.
        adata : sc.AnnData
            The annotated data matrix to be processed.
        **kwargs : dict, optional
            Parameters passed directly to the specific pipeline function.
            Examples: `n_pcs`, `resolution`, `n_hvg`, or `subset_hvg`.

        Returns
        -------
        sc.AnnData
            The processed AnnData object, typically containing calculated 
            embeddings (X_pca, X_umap) and clustering results in .obs.

        Raises
        ------
        KeyError
            If the provided `name` has not been registered via the `@Processing.register` decorator.
        """
        return cls.pipelines[name](adata, **kwargs)
    

@Processing.register("pearson_PCA")
def pearson_pipeline(
    adata,
    n_pcs: int=50,
    resolution: int=1.0,
    n_hvg: int = 3000
):

    sc.experimental.pp.normalize_pearson_residuals(adata)

    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_hvg
    )

    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)

    return adata

@Processing.register("lognorm")
def lognorm_pipeline(
    adata,
    n_hvg=3000,
    subset_hvg=True,
):

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_hvg
    )

    if subset_hvg:
        adata._inplace_subset_var(adata.var.highly_variable)

    sc.pp.scale(adata)
    sc.pp.pca(adata)

    return adata

@Processing.register("scVI_pearson")
def scVI_pipeline(adata):
    ...

