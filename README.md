# 💻 **spatools**

## 🚧 Under construction! Not ready for use yet! Currently experimenting and planning! 🚧

## Developed by Pedro Videira Pinho from National Institute of Câncer (Brazil) (c) 2024
spatools is a Python package developed to facilitate and accelerate the analysis of spatial transcriptomics data. This package was created as part of a Scientific Initiation project at the National Cancer Institute (Brazil), focusing on practical tools for preprocessing, analysis, and visualization of spatial data in transcriptomics experiments. It provides an intuitive interface for spatial data manipulation and multiple visualization functions to assist researchers in studying gene expression patterns across different samples.

## 🧬 Github
The github link: [spatools](https://github.com/pedrovp161/pack_v3.git)
## 🧬 Features
The **`spatools`** package is divided into four main modules, each with a specific purpose:

- **`read`** (📂 **input/output**): Functions for reading and saving spatial data files.
- **`pp`** (🧹 **preprocessing**): Functions for data preprocessing before analysis.
- **`tl`** (🛠️ **tools**): Tools for data manipulation and processing the output of pp, including spatial data integration and cluster operations, as well as image analisys.
- **`pl`** (📊 **plotting**): Functions for data visualization, such as cluster plots, spatial images, and cluster quality plots.
- **`constants`** (🎨 **configuration**): Color definitions and parameters used across various functions in the package.

## Package Structure
The package is organized so that methods can be accessed directly from the main spatools namespace (aliased as st). For example:
- **`Plotting`** (**pl**): st.pl.plot_bar(...)
- **`Tools`** (**tl**): st.tl.process_image(...), st.tl.merge_clusters(...)

## Installation
You can install **`spatools`** directly from [pip](https://pypi.org/):
``` bash
pip install spatools
```

## **Usage**
Here’s a simple example of how to use the package to load data, process it, and generate some visualizations:

## Load data and save
``` python
from spatools.read import Reading, save_spatial_files
import spatools as st
import os

# definindo diretório e output
DIR = os.path.dirname(__file__)

read = Reading(DIR)

adatas_dir = read.list_dict_with_data_free()
print(adatas_dir)

save_spatial_files(adatas_dir=adatas_dir, output_dir=r"D:\path\to\directory\output")


```

## preprocess data and save

```python
from spatools.read import Reading, save_spatial_files
import spatools as st
from copy import deepcopy
import random
import os

DIR = os.path.dirname(__file__)
read = Reading(DIR)

adatas_dir = read.list_dict_with_data_h5ad()
print(adatas_dir)

adatas_dir_raw = deepcopy(adatas_dir)
print(adatas_dir_raw)

random.seed(42)

st.pp.preprocessar(adatas_dir=adatas_dir, save_files=True, output_dir=r'D:\path\to\your\directory\of\output')

# Check summary of data before and after preprocessing
spots_raw, genes_raw = st.pp.check_summary(dicionario=adatas_dir_raw)
print(f"Número de celulas antes {spots_raw}, numero de genes antes {genes_raw}")

spots, genes = st.pp.check_summary(dicionario=adatas_dir)
print(f"Número de celulas depois {spots}, numero de genes depois {genes}")
print(1-spots/spots_raw)
```

...

## Licence
This package is licensed under the [MIT License](https://www.mit.edu/~amini/LICENSE.md).

...
[![DOI](https://zenodo.org/badge/912254487.svg)](https://doi.org/10.5281/zenodo.14611084)