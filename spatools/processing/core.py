from __future__ import annotations

import os
import re
import warnings
from typing import Any, Union

import numpy as np
import scanpy as sc
from anndata import AnnData

from .utils import is_outlier, save_adatas


# ── Types ─────────────────────────────────────────────────────────────────────

AdatasDict = dict[str, AnnData]


# ── QC setup ──────────────────────────────────────────────────────────────────

# Prefixos mitocondriais conhecidos por organismo
_MT_PREFIXES = ("MT-", "mt-")
_ENS_PATTERN  = re.compile(r"^ENS(MUS)?G\d+")


def _detect_mt_prefix(names) -> str | None:
    """
    Retorna o prefixo mitocondrial encontrado em *names*, ou None se ausente.
    Testa MT- (humano) e mt- (camundongo), nessa ordem.
    """
    for prefix in _MT_PREFIXES:
        if names.str.startswith(prefix).any():
            return prefix
    return None


def _is_ensembl_ids(names) -> bool:
    """Retorna True se a maioria dos IDs segue o padrão Ensembl (ENSG.../ENSMUSG...)."""
    sample = names[:min(50, len(names))]
    return sample.str.match(_ENS_PATTERN).mean() > 0.5


def _flag_mt_genes(adata: AnnData) -> None:
    """
    Anota genes mitocondriais em ``adata.var["mt"]``.

    Estratégia (em ordem de prioridade):
    1. Procura prefixo ``MT-`` (humano) em ``var_names`` e ``gene_ids``.
    2. Procura prefixo ``mt-`` (camundongo) em ``var_names`` e ``gene_ids``.
    3. Se nenhum prefixo for encontrado e os IDs forem Ensembl, lança
       ``ValueError`` orientando o usuário a fornecer uma lista manual.

    Parameters
    ----------
    adata : AnnData
        Objeto AnnData a ser anotado in-place.

    Raises
    ------
    ValueError
        Quando os IDs são Ensembl e não há símbolo de gene disponível
        para inferir genes mitocondriais automaticamente.
    """
    # Constrói lista de Series candidatas para busca
    candidates: Any = [adata.var_names]
    if "gene_ids" in adata.var.columns:
        candidates.append(adata.var["gene_ids"].astype(str))

    mt_mask = None

    for names in candidates:
        prefix = _detect_mt_prefix(names)
        if prefix is not None:
            current = names.str.startswith(prefix)
            mt_mask = current if mt_mask is None else (mt_mask | current)

    if mt_mask is not None:
        adata.var["mt"] = mt_mask
        return

    # Nenhum prefixo encontrado — verifica se são Ensembl IDs
    if _is_ensembl_ids(adata.var_names):
        raise ValueError(
            "Genes mitocondriais não puderam ser detectados automaticamente: "
            "os var_names parecem ser Ensembl IDs (ex: ENSG00000...). "
            "Forneça uma lista manual via:\n\n"
            "    adata.var['mt'] = adata.var_names.isin(mt_gene_list)\n\n"
            "e passe o AnnData já anotado para preprocessing()."
        )

    # var_names são símbolos mas não MT/mt — dataset sem genes mitocondriais
    # (ex: dado já filtrado, ou organismo não suportado). Anota tudo como False.
    warnings.warn(
        "Nenhum gene mitocondrial encontrado com prefixos MT- ou mt-. "
        "A coluna 'mt' foi definida como False para todos os genes. "
        "Verifique se o organismo é suportado ou se os genes já foram removidos.",
        UserWarning,
        stacklevel=3,
    )
    adata.var["mt"] = False


def _compute_qc_metrics(adata: AnnData) -> None:
    """Calculate per-spot QC metrics, including mitochondrial percentage."""
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


# ── Individual filters ─────────────────────────────────────────────────────────

def _filter_by_counts(adata: AnnData, k: float) -> tuple[AnnData, int]:
    mask = is_outlier(adata.obs["log1p_total_counts"], method="low", k=k)
    return adata[~mask].copy(), int(mask.sum())


def _filter_by_genes(adata: AnnData, k: float) -> tuple[AnnData, int]:
    mask = is_outlier(adata.obs["log1p_n_genes_by_counts"], method="low", k=k)
    return adata[~mask].copy(), int(mask.sum())


def _filter_by_counts_or_genes(adata: AnnData, k: float) -> tuple[AnnData, int]:
    out_counts = is_outlier(adata.obs["log1p_total_counts"], method="low", k=k)
    out_genes  = is_outlier(adata.obs["log1p_n_genes_by_counts"], method="low", k=k)
    mask = out_counts | out_genes
    return adata[~mask].copy(), int(mask.sum())


def _filter_by_mt_outlier(adata: AnnData, k: float) -> tuple[AnnData, int]:
    mask = is_outlier(adata.obs["pct_counts_mt"], method="high", k=k)
    return adata[~mask].copy(), int(mask.sum())


def _filter_by_mt_threshold(
    adata: AnnData, threshold: float
) -> tuple[AnnData, int]:
    """Remove spots com porcentagem mitocondrial > threshold."""
    if not (0 <= threshold <= 100):
        raise ValueError(f"threshold_mt deve estar entre 0 e 100, recebido: {threshold!r}")
    
    mask = adata.obs["pct_counts_mt"] > threshold
    n_removed = int(mask.sum())
    
    if n_removed > 0:
        print(f"  Removidos {n_removed} spots com pct_counts_mt > {threshold}%")
    
    return adata[~mask].copy(), n_removed


def _drop_unexpressed_genes(adata: AnnData) -> AnnData:
    sc.pp.filter_genes(adata, min_cells=1)
    return adata


# ── Stats ─────────────────────────────────────────────────────────────────────

def _build_initial_stats(sample_id: str, adata: AnnData) -> dict:
    return {
        "sample_id":       sample_id,
        "initial_n_spots": adata.n_obs,
        "initial_n_genes": adata.n_vars,
    }


def _record_filter_step(
    stats: dict, key: str, adata: AnnData, n_removed: int
) -> None:
    stats[f"n_removed_by_{key}"] = n_removed
    stats[f"n_after_{key}"]      = adata.n_obs


# ── Per-sample filtering ───────────────────────────────────────────────────────

def _filter_sample(
    sample_id: str,
    adata: AnnData,
    *,
    counts_outliers: bool,
    genes_outliers: bool,
    genes_and_counts_outliers: bool,
    mt_percentage_outliers: bool,
    threshold_mt: float | None,
    k: float,
) -> AnnData:
    """
    Apply the configured QC filters to a single sample.

    Returns a filtered copy; the original is not modified.
    """
    adata = adata.copy()
    adata.var_names_make_unique()

    _flag_mt_genes(adata)
    _compute_qc_metrics(adata)

    stats = _build_initial_stats(sample_id, adata)

    if counts_outliers:
        adata, n = _filter_by_counts(adata, k)
        _record_filter_step(stats, "counts_outlier", adata, n)

    if genes_outliers:
        adata, n = _filter_by_genes(adata, k)
        _record_filter_step(stats, "genes_outlier", adata, n)

    if genes_and_counts_outliers:
        adata, n = _filter_by_counts_or_genes(adata, k)
        _record_filter_step(stats, "counts_or_genes_outlier", adata, n)

    if mt_percentage_outliers:
        adata, n = _filter_by_mt_outlier(adata, k)
        _record_filter_step(stats, "mt_outlier", adata, n)

    if threshold_mt is not None:
        adata, n = _filter_by_mt_threshold(adata, threshold_mt)
        stats["threshold_mt_value"] = threshold_mt
        _record_filter_step(stats, "mt_threshold", adata, n)

    adata = _drop_unexpressed_genes(adata)

    stats["final_n_spots"] = adata.n_obs
    stats["final_n_genes"] = adata.n_vars
    adata.uns[f"preprocessing_stats_{sample_id}"] = stats

    return adata


# ── Argument validation ────────────────────────────────────────────────────────

def _validate_filter_args(
    genes_outliers: bool,
    counts_outliers: bool,
    genes_and_counts_outliers: bool,
    threshold_mt: float | None,
    save_files: bool,
    output_dir: str,
) -> None:
    if genes_and_counts_outliers and (genes_outliers or counts_outliers):
        raise ValueError(
            "genes_and_counts_outliers already covers genes and counts individually. "
            "Disable genes_outliers and counts_outliers to avoid redundancy."
        )
    if genes_outliers and counts_outliers:
        raise ValueError(
            "Select either genes_outliers or counts_outliers, not both. "
            "Use genes_and_counts_outliers to apply both simultaneously."
        )
    if threshold_mt is not None and not (0 <= threshold_mt <= 100):
        raise ValueError(
            f"threshold_mt must be between 0 and 100, got {threshold_mt!r}."
        )
    if save_files and not output_dir:
        raise ValueError(
            "output_dir must be provided when save_files=True."
        )


# ── Public API ─────────────────────────────────────────────────────────────────

def preprocessing(
    adatas_dict: AdatasDict,
    output_dir: str = "",
    save_files: bool = False,
    genes_outliers: bool = False,
    counts_outliers: bool = False,
    mt_percentage_outliers: bool = True,
    genes_and_counts_outliers: bool = True,
    k: float = 4,
    threshold_mt: Union[float, int, None] = None,
) -> AdatasDict:
    """
    Quality-control filtering for Visium spatial transcriptomics data.

    Applies MAD-based outlier filters and an optional hard mitochondrial
    threshold to each sample independently. Per-sample statistics are stored
    in ``adata.uns["preprocessing_stats_<sample_id>"]``.

    Parameters
    ----------
    adatas_dict : dict[str, AnnData]
        Mapping of sample ID → AnnData (Visium).
    output_dir : str
        Directory for saving .h5ad files. Required when *save_files* is True.
    save_files : bool
        Write filtered AnnData objects to *output_dir*.
    genes_outliers : bool
        Remove spots with low gene counts (MAD-based, lower tail).
    counts_outliers : bool
        Remove spots with low total counts (MAD-based, lower tail).
    genes_and_counts_outliers : bool
        Remove spots that are outliers in genes **or** total counts.
        Mutually exclusive with *genes_outliers* and *counts_outliers*.
    mt_percentage_outliers : bool
        Remove spots with high mitochondrial percentage (MAD-based, upper tail).
    k : float
        MAD multiplier for outlier thresholds (default 4).
    threshold_mt : float or None
        Hard upper cutoff for ``pct_counts_mt`` (0–100). Applied after all
        MAD-based filters. Useful as a fixed biological limit
        (e.g. ``threshold_mt=8`` for T cell datasets).

    Returns
    -------
    dict[str, AnnData]
        Filtered AnnData objects, keyed by sample ID.
    """
    _validate_filter_args(
        genes_outliers=genes_outliers,
        counts_outliers=counts_outliers,
        genes_and_counts_outliers=genes_and_counts_outliers,
        threshold_mt=threshold_mt,
        save_files=save_files,
        output_dir=output_dir,
    )

    filtered = {
        sample_id: _filter_sample(
            sample_id,
            adata,
            counts_outliers=counts_outliers,
            genes_outliers=genes_outliers,
            genes_and_counts_outliers=genes_and_counts_outliers,
            mt_percentage_outliers=mt_percentage_outliers,
            threshold_mt=threshold_mt,
            k=k,
        )
        for sample_id, adata in adatas_dict.items()
    }

    if save_files:
        save_adatas(output_dir, filtered)

    return filtered