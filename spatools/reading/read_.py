import os
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as an
from PIL import Image
from pathlib import Path
from anndata import AnnData

# ==========================================
# 1. Registry System (The Internal "Engine")
# ==========================================
READERS = {}

def register_reader(name: str):
    """Decorator to register new reading pipelines."""
    def decorator(func):
        READERS[name] = func
        return func
    return decorator

# ==========================================
# 2. Helper Functions (Hidden from user)
# ==========================================
def _post_process_adata(adata: sc.AnnData, filter_empty_genes: bool = True) -> sc.AnnData:
    """Internal function to standardize AnnData after reading."""
    adata.var_names_make_unique()
    
    # Removed the rigid 36601 gene limit to make the package universal
    if filter_empty_genes:
        sc.pp.filter_genes(adata, min_cells=1)
    
    return adata

def _read_single_free(dir_path: str) -> sc.AnnData:
    """Raw manual AnnData assembly logic for spatial data (formerly read_free)."""
    pos_path = os.path.join(dir_path, "spatial", "tissue_positions_list.csv")
    matx_path = os.path.join(dir_path, "filtered_feature_bc_matrix")
    json_path = os.path.join(dir_path, "spatial", "scalefactors_json.json")

    # Images
    spatial_path = os.path.join(dir_path, "spatial")
    hier_path = os.path.join(spatial_path, [f for f in os.listdir(spatial_path) if f.endswith("hires_image.png")][0])
    lower_path = os.path.join(spatial_path, [f for f in os.listdir(spatial_path) if f.endswith("lowres_image.png")][0])

    image_hirer = np.array(Image.open(hier_path))
    image_lower = np.array(Image.open(lower_path))

    # Matrix and positions
    adata = sc.read_10x_mtx(matx_path)
    
    try:
        pos_spatial = pd.read_csv(pos_path, header=None)
    except FileNotFoundError:
        pos_spatial = pd.read_csv(f"{'_'.join(pos_path.split('_')[0:-1])}.csv", header=None)

    barcodes = pd.DataFrame(adata.obs.index)
    pos = pd.merge_ordered(barcodes, pos_spatial, how='inner', left_on=0, right_on=0)
    
    # Obs adjustment (array coordinates)
    pos_obs = pos[[0, 2, 3]].copy()
    pos_obs.index = pos_obs[0]
    del pos_obs[0]
    pos_obs = pos_obs.rename(columns={2: "array_row", 3: "array_col"})
    pos_obs.index.name = None
    adata.obs = pos_obs

    # .uns construction
    with open(json_path, "r") as json_file:
        scale_info = json.load(json_file)
        
    adata.uns = {
        'spatial': {
            f'{os.path.basename(dir_path)}': {
                'images': {'hires': image_hirer, 'lowres': image_lower},
                'scalefactors': scale_info,
                'metadata': {'chemistry_description': "Spatial 3' v1", 'software_version': 'spaceranger-1.2.0'}
            }
        }
    }

    # obsm adjustment (pixel coordinates)
    pos_obsm = pos[[0, 5, 4]].copy()
    pos_obsm.index = pos_obsm[0]
    del pos_obsm[0]
    pos_obsm = pos_obsm.rename(columns={4: "array_row", 5: "array_col"})
    adata.obsm["spatial"] = pos_obsm.values

    return adata

# ==========================================
# 3. Registered Pipelines
# ==========================================
@register_reader("visium")
def read_visium_dir(dir_path: str, filter_empty_genes: bool = True) -> dict:
    dictionary = {}
    subfolders = [os.path.join(dir_path, f) for f in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, f))]

    for folder in subfolders:
        name = os.path.basename(folder)
        print(f"Reading Visium subfolder: {name}")
        adata = sc.read_visium(folder)
        dictionary[name] = _post_process_adata(adata, filter_empty_genes)
        
    return dictionary

@register_reader("h5ad")
def read_single_h5ad(item_path: str, filter_empty_genes: bool = True, **kwargs) -> sc.AnnData:
    # Removi a trava de string fixa para ser mais flexível
    if not str(item_path).lower().endswith('.h5ad'):
        raise ValueError("Não é um arquivo .h5ad")
    
    # sc.read_h5ad pode ser sensível a caminhos absolutos no Linux dependendo da versão
    adata = sc.read_h5ad(item_path)
    return _post_process_adata(adata, filter_empty_genes)
            
    return dictionary

@register_reader("free")
def read_free_dir(dir_path: str, filter_empty_genes: bool = True) -> dict:
    dictionary = {}
    subfolders = [os.path.join(dir_path, f) for f in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, f))]

    for folder in subfolders:
        name = os.path.basename(folder)
        print(f"Reading Free format subfolder: {name}")
        adata = _read_single_free(folder)
        dictionary[name] = _post_process_adata(adata, filter_empty_genes)
        
    return dictionary

@register_reader("auto")
def _detect_format(dir_path: str) -> str:
    """
    Inspects the directory to automatically infer the spatial data type.
    """
    if not os.path.exists(dir_path):
        raise FileNotFoundError(f"Directory '{dir_path}' not found.")

    files = os.listdir(dir_path)
    
    # 1. H5AD rule: Is there any file ending with .h5ad?
    if any(f.endswith('.h5ad') for f in files):
        return "h5ad"
        
    # 2. Visium rule: Are there subfolders? Does the first subfolder have a 'spatial' folder?
    subfolders = [os.path.join(dir_path, d) for d in files if os.path.isdir(os.path.join(dir_path, d))]
    
    if subfolders:
        first_folder = subfolders[0]
        # 10x Genomics pattern always creates a 'spatial' folder
        if os.path.exists(os.path.join(first_folder, "spatial")):
            return "visium"
            
    # If we got here, we don't know what it is
    raise ValueError(
        f"Could not automatically detect format in '{dir_path}'. "
        "Check folder structure or manually specify 'format_name'."
    )

# ==========================================
# 4. Public API (Try-Except Cascade)
# ==========================================
def read(dir_path: str = None, **kwargs) -> dict | AnnData:
    dir_path = dir_path or os.getcwd()
    dictionary = {}
    ordem_tentativas = ["h5ad", "visium", "free"]
    
    if not os.path.exists(dir_path):
        print(f"❌ Erro: O diretório {dir_path} não existe.")
        return {}
    
    if dir_path.endswith(".h5ad"):
        return an.read_h5ad(dir_path)
    if os.listdir(dir_path) == ["spatial", "filtered_feature_bc_matrix.h5"] or os.listdir(dir_path) == ["filtered_feature_bc_matrix.h5", "spatial"]:
        return read_free_dir(dir_path)
    

    print(f"🔎 Analisando diretório: {dir_path}")

    for item_name in sorted(os.listdir(dir_path)):
        if item_name.startswith('.'): continue
        
        item_path = os.path.join(dir_path, item_name)
        print(f"  ➜ Processando: '{item_name}' ... ", end="")
        
        sucesso = False
        erros_detalhados = []

        for fmt in ordem_tentativas:
            try:
                adata = READERS[fmt](item_path, **kwargs)
                dictionary[Path(item_name).stem] = adata
                print(f"✅ Lido como [{fmt.upper()}]")
                sucesso = True
                break
            except Exception as e:
                # Guarda o erro do formato que deveria ter funcionado
                if fmt == "h5ad" and item_name.endswith(".h5ad"):
                    erros_detalhados.append(f"Erro no H5AD: {str(e)}")
                continue
                
        if not sucesso:
            msg_erro = " | ".join(erros_detalhados) if erros_detalhados else "Formato incompatível"
            print(f"❌ Falha. ({msg_erro})")

    return dictionary


if __name__ == "__main__":
    import spatools as st
    adata = st.read("/mnt/SATA/spatialPaper/data/newData/P30_GOR_S4")

    import squidpy as sq
    adata.obs["color"] = "#669ce4"
    sq.pl.spatial_scatter(adata, color = "color", alpha = 0.8)