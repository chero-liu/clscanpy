from typing import Optional, Union
import os
import pandas as pd
from pathlib import Path
import scanpy as sc
from scipy import sparse
import sys


def get_single_data_format(filename):
    if os.path.isdir(filename):
        return "10x-format"
    elif filename.endswith(".h5ad"):
        return "h5ad-format"
    else:
        return "normal-format"


def read(
    filename: Union[str, Path],
    var_names: Optional[str] = "gene_symbols",
    make_unique: bool = True,
    cache: bool = False,
    delimiter: str = "\t",
    dtype: str = None,
    prefix: str = None,
):
    if dtype == "10x-format":
        adata = sc.read_10x_mtx(
            filename, var_names=var_names, cache=cache, make_unique=make_unique
        )
    elif dtype == "h5ad-format":
        adata = sc.read(filename)
        if "raw" in adata.layers:
            adata.X = adata.layers["raw"]
        else:
            sys.exit("No raw counts matrix can be found in Anndata object")
    elif dtype == "normal-format":
        try:
            adata = sc.read_csv(filename, delimiter=delimiter).T
        except ValueError:  # BD matrix
            adata = sc.read_csv(filename, delimiter=delimiter, first_column_names=True)
    else:
        raise ValueError("Not a valid matrix format, plz check.")

    gex_rows = list(map(lambda x: x.replace("_", "-"), adata.var.index))
    adata.var.index = gex_rows

    if prefix:
        gex_cols = list(map(lambda x: prefix + "_" + x, adata.obs.index))
        adata.obs.index = gex_cols

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    return adata


def ensure_sparse_matrix(adata):
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    adata.layers["raw"] = adata.X


def validate_and_normalize_adata(
    adata: sc.AnnData, valid_range=(0, 9.22), valid_sum=1e4
) -> sc.AnnData:
    """
    Validate and normalize an AnnData object to ensure it is properly log1p normalized.

    Parameters:
        adata (sc.AnnData): The AnnData object to validate and normalize.
        valid_range (tuple): A tuple specifying the valid range for `adata.X` values (min, max).
        valid_sum (float): The expected sum of `np.expm1(adata.X[0])`.

    Returns:
        sc.AnnData: The validated and normalized AnnData object.

    Raises:
        TypeError: If the data is not properly normalized and no raw data is available for normalization.
    """

    def is_valid_range(adata, valid_range):
        """Check if the data range is within the valid range."""
        return (adata.X[:1000].min() >= valid_range[0]) and (
            adata.X[:1000].max() <= valid_range[1]
        )

    def is_valid_sum(adata, valid_sum):
        """Check if the sum of the first row (after expm1) is close to the expected value."""
        return np.abs(np.expm1(adata.X[0]).sum() - valid_sum) <= 1

    # Perform validation checks
    range_valid = is_valid_range(adata, valid_range)
    sum_valid = is_valid_sum(adata, valid_sum)

    print(f"数据范围检查: {range_valid}")
    print(f"标准化检查: {sum_valid}")

    # If validation fails, attempt normalization
    if not (range_valid and sum_valid):
        print("数据未进行正确的log1p标准化，进行normalize_total和log1p处理...")
        print("检查是否存在raw数据")
        if adata.raw is not None:
            adata = restore_from_raw(adata)
            normalize(adata)
        else:
            raise TypeError("未找到raw数据,请核查数据是否正确")
    else:
        print("数据已经进行了正确的log1p标准化")

    return adata
