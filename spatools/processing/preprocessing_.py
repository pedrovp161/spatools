from pathlib import PurePath

from .core import preprocessing

class Preprocessing:

    pipelines = {}

    @classmethod
    def register(cls, name):
        def decorator(func):
            cls.pipelines[name] = func
            return func
        return decorator
    

    @classmethod
    def run(cls, name, adatas_dict, **kwargs):
        """
        Executes a previously registered preprocessing pipeline.

        This method acts as a dispatcher, looking up the pipeline by name 
        in the class registry and executing it with the provided data.

        Parameters
        ----------
        name : str
            The name of the preprocessing pipeline to execute. Available default 
            pipelines are:
            - 'MAD_low': Filters low outliers for genes and counts separately, 
              plus high mitochondrial percentage.
            - 'MAD_combined': Filters outliers by combining gene and count filters,
              plus high mitochondrial percentage.
            - 'MAD_no_mt': Filters outliers by combining gene and count filters, 
              but ignores the mitochondrial percentage filter.
        adatas_dict : dict
            A dictionary containing spatial transcriptomics data (Visium 10X).
            Keys are sample names and values are AnnData objects.
        **kwargs : dict, optional
            Additional arguments that will be passed down to the core `preprocessing` 
            function. Common arguments include `output_dir` (str) and `save_files` (bool).

        Returns
        -------
        None or str
            Returns a warning string if there is redundancy in the applied filters 
            (e.g., genes_and_counts_outliers and genes_outliers simultaneously). 
            Otherwise, modifies `adatas_dict` in-place and returns None.

        Raises
        ------
        KeyError
            If the provided `name` is not registered in the `cls.pipelines` dictionary.
        """
        return cls.pipelines[name](adatas_dict, **kwargs)
              

@Preprocessing.register("MAD_low")
def pipeline_MAD_low(adatas_dict, **kwargs):

    return preprocessing(
        adatas_dict,
        genes_outliers=True,
        counts_outliers=True,
        mt_percentage_outliers=True,
        genes_and_counts_outliers=False,
        **kwargs
    )

@Preprocessing.register("MAD_combined")
def pipeline_MAD_combined(adatas_dict, **kwargs):

    return preprocessing(
        adatas_dict,
        genes_and_counts_outliers=True,
        mt_percentage_outliers=True,
        **kwargs
    )

@Preprocessing.register("MAD_no_mt")
def pipeline_MAD_no_mt(adatas_dict, **kwargs):

    return preprocessing(
        adatas_dict,
        genes_and_counts_outliers=True,
        mt_percentage_outliers=False,
        **kwargs
    )
