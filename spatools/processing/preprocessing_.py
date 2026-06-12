import os
from enum import Enum
from typing import Dict, Optional, Union

import scanpy as sc

from .core import preprocessing


class PipelineType(Enum):
    """
    Registra os tipos de pipelines disponíveis, suas strings de identificação
    e descrições detalhadas para o usuário final.
    """
    MAD_LOW_READS_AND_GENES = (
        "MAD_LOW_READS_AND_GENES",
        "Filtra outliers inferiores para contagens de reads (counts) e genes detectados separadamente, "
        "além de remover células com alta porcentagem de expressão mitocondrial (padrão 5 MADs)."
    )
    MAD_COMBINED_WITH_MT = (
        "MAD_COMBINED_WITH_MT",
        "Filtra outliers combinando os critérios de genes e contagens em uma métrica conjunta, "
        "além de aplicar o filtro padrão para alta porcentagem mitocondrial."
    )
    MAD_WITHOUT_MT = (
        "MAD_WITHOUT_MT",
        "Aplica o filtro combinado de genes e contagens, mas ignora completamente a "
        "porcentagem de genes mitocondriais no processo de filtragem."
    )
    MAD_CUSTOM_MT = (
        "MAD_CUSTOM_MT",
        "Pipeline flexível: Aplica o filtro combinado de genes e contagens, e permite ao usuário "
        "definir um limite manual (threshold) exato para a porcentagem mitocondrial máxima."
    )

    def __init__(self, name: str, description: str):
        self._value_ = name
        self.description = description

    @classmethod
    def help(cls):
        """Imprime no terminal uma tabela explicativa de todas as pipelines."""
        print("\n" + "="*90)
        print(f"{'PIPELINE':<30} | {'DESCRIÇÃO'}")
        print("="*90)
        for pipeline in cls:
            # Quebra o texto para não estourar a tela do terminal
            desc = pipeline.description
            print(f"{pipeline.value:<30} | {desc[:55]}...")
            if len(desc) > 55:
                print(f"{'':<30} | {desc[55:]}")
            print("-"*90)


class Preprocessing:
    pipelines = {}

    @classmethod
    def register(cls, name: Union[str, PipelineType]):
        def decorator(func):
            key = name.value if isinstance(name, PipelineType) else name
            cls.pipelines[key] = func
            return func
        return decorator
    
    @classmethod
    def run(
        cls, 
        name: Union[str, PipelineType], 
        adatas_dict: Dict[str, sc.AnnData], 
        save_files: bool = False,
        output_dir: Optional[str] = None,
        **kwargs
    ) -> Dict[str, sc.AnnData]:
        """Executa a pipeline escolhida e retorna o dicionário atualizado."""
        pipeline_key = name.value if isinstance(name, PipelineType) else name

        if pipeline_key not in cls.pipelines:
            available = list(cls.pipelines.keys())
            raise ValueError(
                f"Pipeline '{pipeline_key}' não encontrada.\n"
                f"Use `PipelineType.help()` para ver as opções válidas e suas descrições."
            )

        result = cls.pipelines[pipeline_key](
            adatas_dict, 
            save_files=save_files, 
            output_dir=output_dir, 
            **kwargs
        )
        
        # Se o resultado for um dicionário (dict), assume que é o resultado filtrado
        # Caso contrário, assume que é um aviso e retorna o dicionário original
        if isinstance(result, dict) and set(result.keys()) == set(adatas_dict.keys()):
            return result
        else:
            if result:
                print(f"[AVISO Preprocessing]: \n{result}")
            return adatas_dict


# ==============================================================================
# DEFINIÇÃO DAS PIPELINES
# ==============================================================================

@Preprocessing.register(PipelineType.MAD_LOW_READS_AND_GENES)
def pipeline_low_reads_genes(adatas_dict, **kwargs):
    return preprocessing(
        adatas_dict,
        genes_outliers=True,
        counts_outliers=True,
        mt_percentage_outliers=True,
        genes_and_counts_outliers=False,
        **kwargs
    )

@Preprocessing.register(PipelineType.MAD_COMBINED_WITH_MT)
def pipeline_combined_with_mt(adatas_dict, **kwargs):
    return preprocessing(
        adatas_dict,
        genes_and_counts_outliers=True,
        mt_percentage_outliers=True,
        **kwargs
    )

@Preprocessing.register(PipelineType.MAD_WITHOUT_MT)
def pipeline_without_mt(adatas_dict, **kwargs):
    return preprocessing(
        adatas_dict,
        genes_and_counts_outliers=True,
        mt_percentage_outliers=False,
        **kwargs
    )

@Preprocessing.register(PipelineType.MAD_CUSTOM_MT)
def pipeline_custom_mt(adatas_dict, threshold_mt: float = 15.0, **kwargs):
    """
    Pipeline flexível com limite manual para porcentagem mitocondrial.
    Recebe um parâmetro explícito `threshold_mt` (ex: 15.0 para 15%).
    """
    # Valida que threshold_mt é um número válido
    try:
        threshold_mt = float(threshold_mt)
    except (TypeError, ValueError):
        raise ValueError(f"threshold_mt deve ser um número entre 0 e 100, recebido: {threshold_mt!r}")
    
    return preprocessing(
        adatas_dict,
        genes_and_counts_outliers=True,
        mt_percentage_outliers=True, 
        threshold_mt=threshold_mt,      # Aplica o corte fixo do usuário
        **kwargs
    )


# ==============================================================================
# EXEMPLO DE USO
# ==============================================================================
if __name__ == "__main__":
    
    # 1. O usuário está em dúvida? Ele pode rodar isso no Jupyter/Terminal:
    PipelineType.help()

    # Simulando os dados do usuário
    data_dir = "/mnt/SATA/spatialPaper/data/corrected"
    if os.path.exists(data_dir):
        files = [f for f in os.listdir(data_dir) if f.endswith('.h5ad')]
        meu_dict_adatas = {f.replace('.h5ad', ''): sc.read_h5ad(os.path.join(data_dir, f)) for f in files}

        # 2. Executando a nova pipeline customizada:
        # Note que agora o parâmetro mt_threshold aparece no autocomplete!
        meu_dict_processado = Preprocessing.run(
            name=PipelineType.MAD_CUSTOM_MT, 
            adatas_dict=meu_dict_adatas, 
            mt_threshold=10.5, # Usuário escolheu cortar em 10.5%
            save_files=False
        )