# 💻 **spatools**

[![DOI](https://zenodo.org/badge/912254487.svg)](https://doi.org/10.5281/zenodo.14611085)

## 🚧 Under construction! Currently experimenting and planning! 🚧

## Developed by Pedro Videira Pinho at the National Cancer Institute (Brazil) (c) 2024

**spatools** is a Python package built to speed up the analysis of spatial transcriptomics
data. It started as a Scientific Initiation project at the National Cancer Institute (INCA,
Brazil) and focuses on practical tooling for preprocessing, analysis and visualization of
multi-sample spatial experiments — with a particular emphasis on **colocalization**: which
clusters end up next to each other in the tissue, and whether that holds across samples.

## 🧬 Repository

[github.com/pedrovp161/spatools](https://github.com/pedrovp161/spatools)

## 🧬 Modules

| module | alias | what it does |
|---|---|---|
| `reading` | `st.read` | one entry point that detects the format and reads `.h5ad`, Visium folders, or a directory of samples |
| `processing` | `st.pp` | `Preprocessing` (QC pipelines) and `Processing` (normalization, PCA/UMAP, clustering) |
| `tools` | `st.tl` | neighborhood distances, colocalization z-scores, Spearman on deconvolution, cluster merging, spot removal, interactive selection |
| `plotting` | `st.pl` | spatial maps, stacked bars, QC violins, colocalization heatmaps |
| `constants` | `st.con` | ready-made color palettes |

Everything is reachable from the top-level namespace:

```python
import spatools as st

st.read(...)                 # reading
st.pp.Preprocessing.run(...) # preprocessing
st.tl.correlate_distances(...)
st.pl.spatial_plot(...)
st.con.COLORS_23_HEX
```

## Installation

```bash
pip install spatools
```

From source:

```bash
git clone https://github.com/pedrovp161/spatools
cd spatools
pip install -e .
```

Requires Python 3.9+.

## Quick start

### Reading

`st.read()` inspects the path and picks the reader for you:

```python
import spatools as st

adata  = st.read("path/to/sample.h5ad")     # AnnData
adata  = st.read("path/to/visium_sample/")  # AnnData (folder containing spatial/)
adatas = st.read("path/to/all_samples/")    # dict[str, AnnData]
```

### Quality control

```python
import spatools as st

adatas = st.read("data/raw")

# See every available pipeline and what it filters
st.pp.PipelineType.help()

filtered = st.pp.Preprocessing.run(
    name=st.pp.PipelineType.MAD_COMBINED_WITH_MT,
    adatas_dict=adatas,
    k=4,
    save_files=True,
    output_dir="data/filtered",
)

# Per-sample stats land in uns["preprocessing_stats_<sample>"]
st.pl.preprocessing_quality_metrics(filtered)
```

Outlier detection uses the **MAD** (median absolute deviation) rather than mean ± standard
deviation, since the median and MAD are not dragged around by the very outliers being
detected.

### Clustering

```python
adata = st.pp.Processing.run("pearson_PCA", adata, n_pcs=30, resolution=0.5)
# or the classic path:
adata = st.pp.Processing.run("lognorm", adata, n_hvg=3000, subset_hvg=True)
```

### Colocalization

```python
adata = st.tl.correlate_distances(
    adata, is_concatenated=True, cluster_col="clusters_0.6", batch_key="batch"
)
adata = st.tl.z_score(adata)

st.pl.z_score_matrixplot(adata, title="Cluster colocalization")
st.pl.boxplot_cluster_correlations(adata, cluster_col="clusters_0.6", subset=True, value=9)
```

Positive z-score means two clusters are neighbors more often than chance would predict;
negative means they exclude each other.

## 📚 Tutorials

| notebook | covers |
|---|---|
| [Pre-processing](./Tutorials/Pre-processing.ipynb) | the MAD criterion step by step, then `Preprocessing` across many samples, `preprocessing_quality_metrics`, `outlier_quality` |
| [processing](./Tutorials/processing.ipynb) | `Processing` pipelines (`lognorm`, `pearson_PCA`), `merge_clusters`, `remove_spots`, `translate_anndata_genes`, `SelectionTool` |
| [ploting](./Tutorials/ploting.ipynb) | `spatial_plot` in every mode, `bar`, `clusters_quality_violin_boxplot`, color palettes |
| [clustering_correlation_analysis](./Tutorials/clustering_correlation_analysis.ipynb) | `correlate_distances`, `z_score`, `z_score_matrixplot`, `boxplot_cluster_correlations`, Spearman on cell2location output |

## Building and publishing

Metadata lives entirely in `pyproject.toml` (PEP 621) — there is no `setup.py`:

```bash
pip install --upgrade build twine
python -m build
twine check dist/*
twine upload dist/*
```

## Licence

[MIT License](./LICENCE).
