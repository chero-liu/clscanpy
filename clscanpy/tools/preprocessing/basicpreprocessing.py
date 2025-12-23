from anndata import AnnData
from typing import Optional, List, Union
import scanpy as sc
import numpy as np
import pandas as pd
import os


def _filter(
    adata: AnnData,
    mingenes: Optional[float] = None,
    minumis: Optional[float] = None,
    maxgenes: Optional[float] = None,
    maxumis: Optional[float] = None,
    mincells: Optional[int] = None,
) -> AnnData:
    sc.pp.filter_cells(adata, min_genes=mingenes)
    sc.pp.filter_cells(adata, min_counts=minumis)

    maxgenes_new = (
        maxgenes
        if maxgenes > 1
        else np.percentile(adata.obs["n_genes"], 100 * maxgenes) + 1
    )
    maxumis_new = (
        maxumis
        if maxumis > 1
        else np.percentile(adata.obs["n_counts"], 100 * maxumis) + 1
    )

    sc.pp.filter_cells(adata, max_genes=maxgenes_new)
    sc.pp.filter_cells(adata, max_counts=maxumis_new)

    sc.pp.filter_genes(adata, min_cells=mincells)

    return adata


def _load_gene_list(adata: AnnData, gene_file: str) -> Union[np.ndarray, bool]:
    """Load gene list from GMT file or return False if file not found.

    Args:
        adata: AnnData object containing gene names in var_names
        gene_file: Path to GMT file containing gene list

    Returns:
        Boolean array indicating which genes are in the list, or False if file not found
    """
    if not os.path.exists(gene_file):
        print(f"WARNING: {gene_file} not found.")
        return False
    gene_list = get_genelist_from_refgenome_gmt(gene_file)
    return adata.var_names.isin(gene_list)


def get_genelist_from_refgenome_gmt(filepath):
    with open(filepath, "r") as f:
        genes = f.readline().strip().split("\t")[2:]
    return genes


def find_mt_threshold(
    adata: AnnData,
    refgenome: str,
    mt_thresholds: List[float],
    mtfilter: str = "default",
) -> float:
    """
    Calculate mt threshold for each sample (generic version)
    """
    try:
        adata.var["mt"] = _load_gene_list(adata, f"{refgenome}/MT_genelist.gmt")

    except FileNotFoundError:
        print(
            f"WARNING: Mitochondrial gene list file not found: {refgenome}/MT_genelist.gmt"
        )
        adata.var["mt"] = False

    except Exception as e:
        print(f"ERROR: Failed to load mitochondrial gene list: {str(e)}")
        print(f"File path: {refgenome}/MT_genelist.gmt")
        adata.var["mt"] = False

    mt_df = sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=False
    )[0]["pct_counts_mt"]

    mt_count = []
    for ts in mt_thresholds:
        mt_count.append(len(mt_df[mt_df > ts]) / len(mt_df))

    mt_loc = mt_count.index(min(mt_count, key=lambda x: abs(x - 0.05)))
    if mt_count[mt_loc] > 0.3:
        print("WARNING: Filter too much cells according to MT threshold")

    if mtfilter != "default":
        return float(mtfilter)
    else:
        mtfilter = mt_thresholds[mt_loc]
        return float(mtfilter)


def calculate_mt_common(mtfilter: str, mt_list: List[int]) -> int:
    if mtfilter != "default":
        try:
            mt_common = int(mtfilter)
        except (ValueError, TypeError):
            print("Warning: mtfilter is not a valid integer. Using the default value.")
            mt_common = "default"  # You can set a default value or handle it as needed
    else:
        df = pd.DataFrame(mt_list)
        mt_common = df.mode().sort_values(by=0, ascending=False)[0].iloc[0]

    return mt_common


def normalize_total(
    adata: AnnData,
    target_sum: Optional[float] = None,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    key_added: Optional[str] = None,
    layer: Optional[str] = None,
    layers: Optional[List[str]] = None,
    layer_norm: Optional[str] = None,
    inplace: bool = True,
    copy: bool = False,
) -> None:
    """
    Normalize total counts per cell.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    target_sum : float, optional (default: None)
        If None, after normalization, each observation (cell) has a total count equal to the median of the
        *counts_per_cell* before normalization.
    exclude_highly_expressed : bool, optional (default: False)
        Exclude genes that are expressed above `max_fraction` in a cell.
    max_fraction : float, optional (default: 0.05)
        If `exclude_highly_expressed` is True, exclude genes that are expressed above `max_fraction` in a cell.
    key_added : str, optional (default: None)
        Name of the field in `adata.obs` where the total counts per cell are stored.
    layer : str, optional (default: None)
        Name of the layer where the normalized data should be stored. By default, the normalized data is stored in
        `adata.X`.
    layers : List[str], optional (default: None)
        Names of the layers where the normalized data should be stored. If `layers` is not None, `layer` should be None.
    layer_norm : str, optional (default: None)
        If not None, renormalize the values of the layer such that the values sum to `layer_norm`.
    inplace : bool, optional (default: True)
        Whether to place the result back into the `AnnData` object or return it.
    copy : bool, optional (default: False)
        If True, return a copy of `adata` instead of updating it in place.

    Returns
    -------
    None
    """
    sc.pp.normalize_total(
        adata,
        target_sum=target_sum,
        exclude_highly_expressed=exclude_highly_expressed,
        max_fraction=max_fraction,
        key_added=key_added,
        layer=layer,
        layers=layers,
        layer_norm=layer_norm,
        inplace=inplace,
        copy=copy,
    )
    print("Normalization step is finished in adata.X")


def log1p(
    adata: AnnData,
    base: Optional[float] = None,
    copy: bool = False,
    chunked: Optional[bool] = None,
    chunk_size: Optional[int] = None,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
) -> None:
    """
    Logarithmize the data matrix.

    Computes X = log(X + 1), where log denotes the natural logarithm unless a different base is given.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    base : float, optional (default: None)
        Base of the logarithm. Natural logarithm is used by default.
    copy : bool, optional (default: False)
        If an AnnData is passed, determines whether a copy is returned.
    chunked : bool, optional (default: None)
        Process the data matrix in chunks, which will save memory. Applies only to AnnData.
    chunk_size : int, optional (default: None)
        n_obs of the chunks to process the data in.
    layer : str, optional (default: None)
        Entry of layers to transform.
    obsm : str, optional (default: None)
        Entry of obsm to transform.

    Returns
    -------
    None
    """
    sc.pp.log1p(
        adata,
        base=base,
        copy=copy,
        chunked=chunked,
        chunk_size=chunk_size,
        layer=layer,
        obsm=obsm,
    )
    print("Log transformation step is finished in adata.X")


def normalize(adata: AnnData, target_sum: float = 1e4) -> None:
    """
    Normalize the data matrix by total counts and logarithmize the result.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    target_sum : float, optional (default: 1e4)
        If None, after normalization, each observation (cell) has a total count equal to the median of the
        *counts_per_cell* before normalization.

    Returns
    -------
    None
    """
    normalize_total(adata, target_sum=target_sum)
    log1p(adata)


def scale(
    adata: Union[AnnData],
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    Wrap function of scanpy.pp.scale

    Scale data to unit variance and zero mean.
    .. note::
        Variables (genes) that do not display any variation (are constant across
        all observations) are retained and set to 0 during this operation. In
        the future, they might be set to NaNs.
    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    zero_center
        If `False`, omit zero-centering variables, which allows to handle sparse
        input efficiently.
    max_value
        Clip (truncate) to this value after scaling. If `None`, do not clip.
    copy
        If an :class:`~anndata.AnnData` is passed,
        determines whether a copy is returned.
    Returns
    -------
    Depending on `copy` returns or updates `adata` with a scaled `adata.X`.
    """

    sc.pp.scale(adata, zero_center=zero_center, max_value=max_value, copy=copy)

    print("Scale step is finished in adata.X")


def filter_genes_for_pyscenic(
    adata: sc.AnnData, feather_files: list, minCountsPerGene: int = 1
) -> sc.AnnData:
    """
    为 pySCENIC 分析预处理基因过滤

    Parameters
    ----------
    adata : AnnData
        包含 raw counts 的 AnnData 对象，需有 layers['raw']
    feather_files : list
        RcisTarget 数据库 feather 文件路径列表
    minCountsPerGene : int, optional
        基因最低总计数阈值，默认为1

    Returns
    -------
    AnnData
        过滤后的 AnnData 对象
    """
    # Load the databases and get the union of all gene names
    file_type = argparse.FileType("rb")
    dbs = _load_dbs([file_type(file) for file in feather_files])

    all_genes_union = set()
    for db in dbs:
        all_genes_union.update(db.genes)

    genes_in_db = [gene.upper() for gene in all_genes_union]

    minSamples = int(adata.n_obs * 0.01)
    exprMat = adata.layers["raw"]

    if issparse(exprMat):
        n_counts_per_gene = exprMat.sum(axis=0).A1
        n_cells_per_gene = exprMat.getnnz(axis=0)
    else:
        n_counts_per_gene = exprMat.sum(axis=0)
        n_cells_per_gene = (exprMat > 0).sum(axis=0)

    print(f"Genes initially: {adata.n_vars}")
    print(f"Maximum count value: {n_counts_per_gene.max()}")
    print(f"Mean counts per gene: {n_counts_per_gene.mean():.1f}")
    detected_ratio = n_cells_per_gene.sum() / (adata.n_obs * adata.n_vars)
    print(f"Detection ratio: {detected_ratio:.2%}")

    count_mask = n_counts_per_gene > minCountsPerGene
    print(f"Genes after count filter: {count_mask.sum()}")

    cell_mask = n_cells_per_gene > minSamples
    combined_mask = count_mask & cell_mask
    print(f"Genes after cell filter: {combined_mask.sum()}")

    var_names_upper = adata.var_names.str.upper()
    db_mask = var_names_upper.isin(genes_in_db)
    final_mask = combined_mask & db_mask

    print(f"Genes in database: {db_mask.sum()}")
    print(f"Final kept genes: {final_mask.sum()}")
    print(f"Cells : {adata.n_obs}")

    return adata[:, final_mask].copy()
