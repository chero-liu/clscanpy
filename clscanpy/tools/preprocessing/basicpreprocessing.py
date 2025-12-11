from anndata import AnnData
from typing import Optional, List, Union, Tuple
import scanpy as sc
import numpy as np
import math
import os
import pandas as pd
from scipy import stats
from scipy.sparse import issparse
from pyscenic.cli.pyscenic import _load_dbs
import argparse
import logging
LOGGER = logging.getLogger(__name__)

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

    LOGGER.info(f"Genes initially: {adata.n_vars}")
    LOGGER.info(f"Maximum count value: {n_counts_per_gene.max()}")
    LOGGER.info(f"Mean counts per gene: {n_counts_per_gene.mean():.1f}")
    detected_ratio = n_cells_per_gene.sum() / (adata.n_obs * adata.n_vars)
    LOGGER.info(f"Detection ratio: {detected_ratio:.2%}")

    count_mask = n_counts_per_gene > minCountsPerGene
    LOGGER.info(f"Genes after count filter: {count_mask.sum()}")

    cell_mask = n_cells_per_gene > minSamples
    combined_mask = count_mask & cell_mask
    LOGGER.info(f"Genes after cell filter: {combined_mask.sum()}")
    
    var_names_upper = adata.var_names.str.upper()
    db_mask = var_names_upper.isin(genes_in_db)
    final_mask = combined_mask & db_mask

    LOGGER.info(f"Genes in database: {db_mask.sum()}")
    LOGGER.info(f"Final kept genes: {final_mask.sum()}")
    LOGGER.info(f"Cells : {adata.n_obs}")

    return adata[:, final_mask].copy()


def _load_gene_list(adata: AnnData, gene_file: str) -> Union[np.ndarray, bool]:
    """Load gene list from GMT file or return False if file not found.

    Args:
        adata: AnnData object containing gene names in var_names
        gene_file: Path to GMT file containing gene list

    Returns:
        Boolean array indicating which genes are in the list, or False if file not found
    """
    if not os.path.exists(gene_file):
        LOGGER.info(f"WARNING: {gene_file} not found.")
        return False
    gene_list = getGenelistFromRefGenomeGmt(gene_file)
    return adata.var_names.isin(gene_list)


def getGenelistFromRefGenomeGmt(filepath):
    with open(filepath, "r") as f:
        genes = f.readline().strip().split("\t")[2:]
    return genes


def _qc(
    adata: AnnData,
    refgenome: Optional[str] = None,
    genelist: List[str] = ["MT_genelist.gmt", "HB_genelist.gmt"],
) -> Tuple[List[str], AnnData]:
    """Perform quality control metrics calculation on single-cell RNA-seq data.

    Calculates standard QC metrics including:
    - Mitochondrial gene percentage
    - Hemoglobin gene percentage (if reference provided)
    - log10(genes per UMI)
    - Total UMI counts (nCount_RNA)
    - Detected genes (nFeature_RNA)

    Args:
        adata: AnnData object containing raw counts data
        refgenome: Path to reference genome directory containing gene lists
        genelist: List of GMT files for mitochondrial and hemoglobin genes

    Returns:
        Tuple containing:
        - List of QC metric names calculated
        - Modified AnnData object with QC metrics added to obs and var
    """
    # Handle mitochondrial genes
    mt_genes = _load_gene_list(adata, f"{refgenome}/{genelist[0]}")
    hb_genes = _load_gene_list(adata, f"{refgenome}/{genelist[1]}")

    # Calculate QC metrics
    metrics = ["nFeature_RNA", "nCount_RNA"]

    if isinstance(mt_genes, np.ndarray):
        adata.var["mt"] = mt_genes
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        adata.obs.rename(columns={"pct_counts_mt": "percent.mito"}, inplace=True)
        adata.obs["percent.mito"] /= 100
        metrics.insert(0, "percent.mito")

    if isinstance(hb_genes, np.ndarray):
        adata.var["hb"] = hb_genes
        sc.pp.calculate_qc_metrics(adata, qc_vars=["hb"], inplace=True)
        adata.obs.rename(columns={"pct_counts_hb": "percent.HB"}, inplace=True)
        adata.obs["percent.HB"] /= 100
        metrics.insert(0, "percent.HB")

    if not isinstance(mt_genes, np.ndarray) and not isinstance(hb_genes, np.ndarray):
        sc.pp.calculate_qc_metrics(adata, inplace=True)

    adata.obs.rename(
        columns={
            "total_counts": "nCount_RNA",
            "n_genes_by_counts": "nFeature_RNA",
        },
        inplace=True,
    )

    return metrics, adata


def _filter(
    adata: AnnData,
    mingenes: Optional[float] = None,
    minumis: Optional[float] = None,
    maxgenes: Optional[float] = None,
    maxumis: Optional[float] = None,
    mincells: Optional[int] = None,
    hbfilter: Optional[float] = None,
) -> AnnData:
    """Filter cells and genes based on quality control metrics.

    Performs sequential filtering of cells and genes based on various QC metrics:
    - Minimum/maximum genes per cell
    - Minimum/maximum UMIs per cell
    - Minimum cells expressing a gene
    - log10(genes/UMI) ratio thresholds
    - Hemoglobin percentage threshold

    Args:
        adata: AnnData object containing QC metrics
        mingenes: Minimum genes per cell (absolute or percentile)
        minumis: Minimum UMIs per cell (absolute or percentile)
        maxgenes: Maximum genes per cell (absolute or percentile)
        maxumis: Maximum UMIs per cell (absolute or percentile)
        mincells: Minimum cells expressing a gene
        hbfilter: Maximum hemoglobin percentage threshold

    Returns:
        Filtered AnnData object

    Notes:
        - For maxgenes/maxumis:
          If value > 1: used as absolute threshold
          If value <= 1: interpreted as percentile
        - Filtering order:
          1. Cell filtering by gene/UMI counts
          2. Gene filtering by cell counts
          3. Additional cell filtering by log10(genes/UMI) and hemoglobin
    """
    sc.pp.filter_cells(adata, min_genes=mingenes)
    sc.pp.filter_cells(adata, min_counts=minumis)

    maxgenes_new = (
        maxgenes
        if maxgenes > 1
        else np.percentile(adata.obs["nFeature_RNA"], 100 * maxgenes) + 1
    )
    maxumis_new = (
        maxumis
        if maxumis > 1
        else np.percentile(adata.obs["nCount_RNA"], 100 * maxumis) + 1
    )

    sc.pp.filter_cells(adata, max_genes=maxgenes_new)
    sc.pp.filter_cells(adata, max_counts=maxumis_new)

    sc.pp.filter_genes(adata, min_cells=mincells)

    if hbfilter is not None:
        adata = adata[adata.obs["percent.HB"] < hbfilter, :]

    return adata


def _mtfilter(
    adata: AnnData,
    mt_thresholds: List[float],
    mtfilter: str = "default",
) -> float:
    """
    Calculate mt threshold for each sample.
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    mt_thresholds : List[float]
        List of mt thresholds for each sample.
    mtfilter : str, optional (default: "default")
        If "default", calculate the mt threshold based on the 90th percentile of the mt thresholds.
        If a float, use the provided value as the mt threshold.
    Returns
    -------
    float
        The mt threshold, adata after filtering.
    """

    if mtfilter == "default":
        max_mito = mt_thresholds.max()
        mtfilter = math.ceil(max_mito / 0.05) * 0.05
        LOGGER.info(f"各样本线粒体比例q90最大值为: {max_mito}")
        if mtfilter > 0.3:
            LOGGER.info("Warning: The 90th percentile of percent.mito exceed 30%!")
            mtfilter = 0.3
        LOGGER.info(f"The percent.mito cut-off is set as {mtfilter}.")

    if mtfilter != "default":
        mtfilter = float(mtfilter)
    else:
        mtfilter = float(mtfilter)

    adata = adata[adata.obs["percent.mito"] < mtfilter, :]
    return mtfilter, adata


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
    LOGGER.info("Normalization step is finished in adata.X")


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
    LOGGER.info("Log transformation step is finished in adata.X")


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

    LOGGER.info("Scale step is finished in adata.X")


def is_outlier(
    metric: np.ndarray,
    n: int = 3,
    type: str = "both",
    cut_1: str = "median",
    cut_2: str = "mad",
    subset: np.ndarray = None,
    batch: np.ndarray = None,
    min_diff: float = None,
) -> np.ndarray:
    """
    Python implementation of the isOutlier function from R

    Parameters:
    metric: array-like, the metric values
    n: number of standard deviations for threshold (default: 3)
    type: "both", "lower", or "higher" (default: "both")
        - "both": detect both lower and upper outliers
        - "lower": detect only lower outliers (values below threshold)
        - "higher": detect only higher outliers (values above threshold)
    cut_1: "median" or "mean" for central tendency (default: "median")
    cut_2: "sd", "median", or "mad" for dispersion (default: "mad")
    subset: subset of data to use for threshold calculation
    batch: batch information for each observation
    min_diff: minimum difference for threshold calculation

    Returns:
    Boolean array indicating outliers
    """
    metric = np.array(metric)
    N = len(metric)

    if batch is None:
        batch = np.array(["1"] * N)
        nobatch = True
    else:
        batch = np.array(batch).astype(str)
        nobatch = False
        if len(batch) != N:
            raise ValueError("length of 'batch' must equal length of 'metric'")

    if subset is not None:
        M = metric[subset]
        B = batch[subset]
    else:
        M = metric
        B = batch

    # Remove NaN values
    nan_mask = np.isnan(M)
    if np.any(nan_mask):
        M = M[~nan_mask]
        B = B[~nan_mask]
        warnings.warn("missing values ignored during outlier detection")

    # Calculate thresholds by batch
    unique_batches = np.unique(B)
    cur_1 = {}
    cur_2 = {}

    for batch_id in unique_batches:
        batch_mask = B == batch_id
        batch_data = M[batch_mask]

        if cut_1 == "median":
            cur_1[batch_id] = np.median(batch_data)
        elif cut_1 == "mean":
            cur_1[batch_id] = np.mean(batch_data)
        else:
            raise ValueError(f"cut_1 must be 'median' or 'mean', got {cut_1}")

        if cut_2 == "sd":
            cur_2[batch_id] = np.std(batch_data)
        elif cut_2 == "median":
            cur_2[batch_id] = np.median(batch_data)
        elif cut_2 == "mad":
            if cut_1 == "median":
                cur_2[batch_id] = stats.median_abs_deviation(batch_data, scale="normal")
            else:
                cur_2[batch_id] = stats.median_abs_deviation(
                    batch_data - cur_1[batch_id], scale="normal"
                )
        else:
            raise ValueError(f"cut_2 must be 'sd', 'median', or 'mad', got {cut_2}")

    # Calculate limits
    lower_limit = {}
    upper_limit = {}

    for batch_id in unique_batches:
        batch_mask = B == batch_id
        batch_data = M[batch_mask]

        diff_val = max(
            min_diff if min_diff is not None else -np.inf, n * cur_2[batch_id]
        )

        # Lower limit
        lower_val = cur_1[batch_id] - diff_val
        batch_min = np.min(batch_data)
        lower_limit[batch_id] = max(lower_val, batch_min)

        # Upper limit
        upper_val = cur_1[batch_id] + diff_val
        batch_max = np.max(batch_data)
        upper_limit[batch_id] = min(upper_val, batch_max)

    # Adjust limits based on type
    if type == "lower":
        for batch_id in unique_batches:
            upper_limit[batch_id] = np.inf
    elif type == "higher":
        for batch_id in unique_batches:
            lower_limit[batch_id] = -np.inf

    LOGGER.info(f"upper_limit: {upper_limit}")
    LOGGER.info(f"lower_limit: {lower_limit}")
    # Identify outliers
    outliers = np.zeros(N, dtype=bool)
    for i in range(N):
        batch_id = batch[i]
        if batch_id in lower_limit and batch_id in upper_limit:
            outliers[i] = (metric[i] < lower_limit[batch_id]) or (
                metric[i] > upper_limit[batch_id]
            )

    return outliers


def find_outliers(
    adata: sc.AnnData,
    vars: list,
    type: list,
    batch: str = None,
    cut_1: str = "median",
    cut_2: str = "mad",
    n: int = 3,
) -> dict:
    """
    Python implementation of the FindOutliers function from R

    Parameters:
    adata: AnnData object
    vars: list of variables to check for outliers
    batch: batch column name in adata.obs
    type: "both", "lower", or "higher"
    cut_1: "median" or "mean"
    cut_2: "sd", "median", or "mad"
    n: number of standard deviations for threshold

    Returns:
    Dictionary with outlier results for each variable
    """
    outliers_dict = {}

    if batch is not None:
        batch_data = adata.obs[batch].values
    else:
        batch_data = None

    for var, current_type in zip(vars, type):
        if var not in adata.obs.columns:
            LOGGER.info(f"Warning: Variable {var} not found in adata.obs")
            continue

        metric_data = adata.obs[var].values

        # Skip if all values are the same
        if np.min(metric_data) == np.max(metric_data):
            LOGGER.info(f"Skipping {var}: all values are identical")
            continue

        LOGGER.info(var)
        outliers = is_outlier(
            metric_data,
            n=n,
            type=current_type,
            cut_1=cut_1,
            cut_2=cut_2,
            batch=batch_data,
        )
        outliers_dict[var] = outliers

    return outliers_dict
