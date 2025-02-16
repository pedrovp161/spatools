import os
import json
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image
from anndata._core.views import ImplicitModificationWarning


def save_spatial_files(output_dir: str, adatas_dir: dict):
    # supress irrelevant warning
    warnings.filterwarnings("ignore", message="Trying to modify attribute", category=ImplicitModificationWarning)
    # Verifica se o diretório de saída existe, se não, cria
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for name, adata in adatas_dir.items():
        output_file_path = os.path.join(output_dir, f"{name}.h5ad")
        sc.write(output_file_path, adata) #type: ignore


class Reading():
    def __init__(self, dir_path: str = "") -> None:
        # If no path is passed, use the default
        self.DIR = dir_path if dir_path else os.getcwd()  # Uses the current working directory if None
                
        # Make sure the directory exists
        os.makedirs(self.DIR, exist_ok=True)

    def list_path_for_archives(self, pasta: str):

        path_files = []

        # List the file names in the folder
        path_names = os.listdir(pasta)

        # Creates the full paths to each file
        for path_name in path_names:
            complete_path = os.path.join(pasta, path_name)
            path_files.append(complete_path)

        return path_files
    
    def list_dict_with_data_free(self):
        dictionary = {}
        
        subfolder_paths = [os.path.join(self.DIR, name) for name in os.listdir(self.DIR) if os.path.isdir(os.path.join(self.DIR, name))]

        for subfolder_path in subfolder_paths:
            subfolder_name = os.path.basename(subfolder_path)

            print("Reading subfolder:", subfolder_name) # Directory verification

            dictionary[subfolder_name] = self.read_free(subfolder_path)
        
        for key in dictionary.keys():
            if dictionary[key].n_vars == 36601:
                sc.pp.filter_genes(dictionary[key], min_cells=1)#type:ignore

        return dictionary
    
    def list_dict_with_data_visium(self):
        dictionary = {}
        
        caminhos_subpastas = [os.path.join(self.DIR, nome) for nome in os.listdir(self.DIR) if os.path.isdir(os.path.join(self.DIR, nome))]

        for caminho_subpasta in caminhos_subpastas:
            nome_subpasta = os.path.basename(caminho_subpasta)

            print("Lendo subpasta:", nome_subpasta) # Directory verification

            dictionary[nome_subpasta] = sc.read_visium(caminho_subpasta)
        
        # remove NaN
        for key in dictionary.keys():
            if dictionary[key].n_vars == 36601:
                sc.pp.filter_genes(dictionary[key], min_cells=1)#type:ignore

        return dictionary
    
    def list_dict_with_data_h5ad(self):
        dictionary = {}
        
        caminhos_arquivos = [os.path.join(self.DIR, nome) for nome in os.listdir(self.DIR) if nome.endswith('.h5ad')]
        
        for caminho_arquivo in caminhos_arquivos:
            nome_arquivo = os.path.basename(caminho_arquivo)
            
            print("Lendo arquivo:", nome_arquivo)  # Verificando o arquivo que está sendo lido

            try:
                adata = sc.read_h5ad(caminho_arquivo)
                dictionary[nome_arquivo] = adata
            except Exception as e:
                print(f"Erro ao ler o arquivo {nome_arquivo}: {e}")
                continue
        
        for key, adata in dictionary.items():
            if adata.n_vars == 36601:
                sc.pp.filter_genes(adata, min_cells=1) #type: ignore

        return dictionary


    def read_free(self, DIR: str):
        # File paths
        pos_path = os.path.join(DIR, "spatial", "tissue_positions_list.csv")
        matx_path = os.path.join(DIR, "raw_feature_bc_matrix")
        json_path = os.path.join(DIR, "spatial", "scalefactors_json.json")

        # Finding image files
        spatial_path = os.path.join(DIR, "spatial")
        hier_files = [file for file in os.listdir(spatial_path) if file.endswith("hires_image.png")]
        lower_files = [file for file in os.listdir(spatial_path) if file.endswith("lowres_image.png")]
        hier_path = os.path.join(spatial_path, hier_files[0])  # Assuming there's only one hires image
        lower_path = os.path.join(spatial_path, lower_files[0])  # Assuming there's only one lowres image

        # transforming image into numpy array
        imh = Image.open(hier_path)
        iml = Image.open(lower_path)
        image_hirer = np.array(imh)
        image_lower = np.array(iml)

        # Reading the matrix file and transposing it
        adata = sc.read_10x_mtx(matx_path)#type:ignore

        # Reading positional information
        pos_spatial = pd.read_csv(pos_path, header=None)

        barcodes = pd.DataFrame(adata.obs.index)

        # Merging positional information with barcodes
        pos = pd.merge_ordered(barcodes, pos_spatial, how='inner', left_on=0, right_on=0)
        
        # Selecting and renaming columns for positional information
        pos = pos[[0, 2, 3]]
        pos.index = pos[0]#type:ignore
        del pos[0]
        pos = pos.rename(columns={2:"array_row", 3:"array_col"})
        pos.index.name = None

        # Setting up the adata.obs
        adata.obs = pos

        # Making the uns model
        with open(json_path, "r") as arquivo_json:
            scale_info = json.load(arquivo_json)
        modelo_uns = {
            'spatial': {
                f'{os.path.basename(DIR)}': {
                    'images': {
                        'hires': image_hirer,
                        'lowres': image_lower
                    },
                    'scalefactors': scale_info,
                    'metadata': {
                        'chemistry_description': "Spatial 3' v1",
                        'software_version': 'spaceranger-1.2.0'
                    }
                }
            }
        }

        adata.uns = modelo_uns  # Assigning the uns model to adata.uns

        # spatial coordinates
        pos = pd.read_csv(pos_path, header=None)
        pos = pd.merge_ordered(barcodes, pos, how='inner', left_on=0, right_on=0)
        # spatial coordinates
        pos = pos[[0, 5, 4]]
        pos.index = pos[0]#type:ignore
        del pos[0]
        pos = pos.rename(columns={4:"array_row", 5:"array_col"})
        pos.index.name = None

        adata.obsm["spatial"] = pos.values

        # Remove NaN values
        adata.var_names_make_unique()
        sc.pp.filter_genes(adata, min_cells=1)#type:ignore

        return adata
    

if __name__ == "__main__":# test
    
    from copy import deepcopy
    import random
    import spatools.preprocessing.pp as pp

    DIR = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__),"..")))
    read = Reading(DIR)

    adatas_dir = read.list_dict_with_data_h5ad()
    print(adatas_dir)

    adatas_dir_raw = deepcopy(adatas_dir)
    print(adatas_dir_raw)

    # Random seed for reproducibility
    r = random.seed(42)
    print(f"seed utilizada: {r}")

    pp.preprocessar(adatas_dir=adatas_dir, save_files=False, output_dir=r'D:\pack_v1\data\filtered')

    # Check summary of data before and after preprocessing
    spots_raw, genes_raw = pp.check_summary(dicionario=adatas_dir_raw)
    print(f"Número de celulas antes {spots_raw}, numero de genes antes {genes_raw}")

    spots, genes = pp.check_summary(dicionario=adatas_dir)
    print(f"Número de celulas depois {spots}, numero de genes depois {genes}")
    print(1-spots/spots_raw)