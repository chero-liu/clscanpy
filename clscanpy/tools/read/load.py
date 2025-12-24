import scanpy as sc
import numpy as np
from typing import Tuple, Optional, List, Union
from pathlib import Path
from clscanpy.tools.read.utils import (
    validate_and_normalize_adata,
)
from anndata import AnnData
import pandas as pd
import warnings
from clscanpy.tools.color.colors import add_color_to_dataframe
import re
import os
import logging

LOGGER = logging.getLogger(__name__)


def standardize_colors(
    adata: AnnData,
    use_col: List[str],
    palette: Optional[List[str]] = None,
) -> AnnData:
    """
    Add standardized color mappings to an AnnData object

    Parameters:
        adata: AnnData object
        use_col: List of column names to add color mappings to
        palette: Optional list of palette names corresponding to each column
                 If not provided, uses default palettes based on column name

    Returns:
        Updated AnnData object with color columns in obs ({col}_col)

    Usage:
        adata = standardize_colors(
            adata,
            ['cluster', 'batch'],
            palette=['tableau20', 'set1']
        )
    """
    # Input validation
    if not isinstance(adata, AnnData):
        raise TypeError("Input must be an AnnData object")

    if palette is not None:
        if len(palette) != len(use_col):
            raise ValueError(
                "palette list length must match use_col length. "
                f"Got {len(palette)} palettes for {len(use_col)} columns"
            )

    # Default palettes for common columns
    SPECIAL_PALETTES = {
        "group": "group_default",
        "sample": "cellpaper_custom",
        "cluster": "customecol2_light",
        "celltype": "tableau20",
    }

    # Create palette list with defaults if not provided
    palettes = []
    for colname in use_col:
        if palette is not None:
            palettes.append(palette[use_col.index(colname)])
        else:
            palettes.append(SPECIAL_PALETTES.get(colname, "default"))

    # Process each target column
    for colname, current_palette in zip(use_col, palettes):
        if colname not in adata.obs.columns:
            warnings.warn(f"Column {colname} not found in obs, skipping", UserWarning)
            continue

        # Ensure categorical type for color mapping
        adata = standardize_levels(adata, colname)
        color_col = f"{colname}_col"

        # Create or update color mapping
        if color_col not in adata.obs.columns:
            try:
                adata.obs = add_color_to_dataframe(
                    adata.obs,
                    colname,
                    palette=current_palette,
                )
                LOGGER.info(
                    f"Added color mapping for '{colname}' using palette '{current_palette}'"
                )
            except ValueError as e:
                warnings.warn(
                    f"Color mapping failed for '{colname}': {str(e)}", RuntimeWarning
                )
        else:
            # Update if category-color mapping is out of sync
            if not is_two_group_equality(adata.obs, colname, color_col):
                LOGGER.info(
                    f"Updating color mapping for '{colname}' (previous mapping out of sync)"
                )
                del adata.obs[color_col]
                adata.obs = add_color_to_dataframe(
                    adata.obs,
                    colname,
                    palette=current_palette,
                )
    return adata


def subsetH5AD(
    adata,
    sample=None,
    group=None,
    cluster=None,
    celltype=None,
    predicate=None,
):
    """
    使用标准化分隔符处理函数增强的AnnData子集筛选

    Parameters:
        adata : AnnData
            输入的AnnData对象
        sample : str (可选)
            逗号分隔的样本ID，支持!前缀排除，如 "S1,!S2"
        group : str (可选)
            逗号分隔的分组名称，支持!前缀排除
        cluster : str (可选)
            逗号分隔的细胞簇ID，支持!前缀排除
        celltype : str (可选)
            逗号分隔的细胞类型，支持!前缀排除
        predicate : str (可选)
            Pandas风格的逻辑表达式
    """

    def parse_exclusion(values):
        """使用standardize_delimited_vector解析包含/排除值"""
        if values is None:
            return None, None

        # 标准化输入格式
        std_values = standardize_delimited_vector(values, ",")

        # 分离包含和排除值
        include = [v for v in std_values if not v.startswith("!")]
        exclude = [v[1:] for v in std_values if v.startswith("!")]

        return (include or None), (exclude or None)

    # 处理每个过滤条件
    filter_params = [
        ("celltype", celltype),
        ("sample", sample),
        ("group", group),
        ("cluster", cluster),
    ]

    for col, value in filter_params:
        if value is not None:
            include_vals, exclude_vals = parse_exclusion(value)

            if include_vals:
                adata = adata[adata.obs[col].isin(include_vals)].copy()

            if exclude_vals:
                adata = adata[~adata.obs[col].isin(exclude_vals)].copy()

    if predicate:
        try:
            predicate = predicate.replace("== [", "in [")
            mask = adata.obs.eval(predicate)
            adata = adata[mask].copy()
        except Exception as e:
            raise ValueError(f"逻辑表达式错误: {predicate}\n错误详情: {str(e)}")

    # Standardize levels for relevant columns
    use_col = [
        "sample",
        "group",
        "cluster",
        "celltype",
    ]

    variables = get_predicate_vars(predicate)
    if isinstance(variables, list):
        use_col.extend(variables)
    else:
        use_col.append(variables)
    use_col = list(set([col for col in use_col if col is not None]))

    for col in use_col:
        if col in adata.obs.columns:
            adata = standardize_levels(adata, col)

    return adata


def is_two_group_equality(df, group_col1="group", group_col2="group_col"):
    """
    检查两个分组列的唯一值数量及其组合数是否一致

    Args:
        df: 输入数据框 (pd.DataFrame)
        group_col1: 第一个分组列名 (默认 "group")
        group_col2: 第二个分组列名 (默认 "group_col")

    Returns:
        bool: 三个统计量(列1唯一值数/列2唯一值数/组合数)是否全部相等

    Raises:
        ValueError: 输入不是DataFrame或列不存在
    """
    if not isinstance(df, pd.DataFrame):
        raise ValueError("输入必须是pandas DataFrame")

    for col in [group_col1, group_col2]:
        if col not in df.columns:
            raise ValueError(f"列 '{col}' 不存在于数据框中")

    n_unique1 = df[group_col1].nunique()
    n_unique2 = df[group_col2].nunique()
    n_combined = df[[group_col1, group_col2]].drop_duplicates().shape[0]

    msgs = []
    if n_unique1 != n_unique2:
        msgs.append(
            f"{group_col1}({n_unique1}) 和 {group_col2}({n_unique2}) 唯一值数量不等"
        )
    if n_unique2 != n_combined:
        msgs.append(f"{group_col2}({n_unique2}) 与组合数({n_combined}) 不等")
    if n_unique1 != n_combined:
        msgs.append(f"{group_col1}({n_unique1}) 与组合数({n_combined}) 不等")

    if msgs:
        LOGGER.info("组合差异检测结果:")
        for msg in msgs:
            LOGGER.info(f"  ! {msg}")

    return n_unique1 == n_unique2 == n_combined


def standardize_levels(
    adata: AnnData, colname: str, colname_levels: str = None
) -> AnnData:
    """
    Standardize the levels of a categorical column in the AnnData object's metadata.

    Parameters:
        adata (AnnData): The AnnData object to modify.
        colname (str): The column name in `adata.obs` to standardize.
        colname_levels (str or None): Comma-separated levels to enforce, or None to infer levels.

    Returns:
        AnnData: The updated AnnData object with standardized levels.
    """
    if colname not in adata.obs.columns:
        raise ValueError(f"Column '{colname}' not found in adata.obs.")

    # 检查是否已是Categorical且无新顺序要求
    if (
        isinstance(adata.obs[colname].dtype, pd.CategoricalDtype)
        and colname_levels is None
    ):
        return adata
    else:
        # Helper function to sort numeric-like strings
        def sort_char_numbers(values):
            try:
                return sorted(
                    values,
                    key=lambda x: (
                        float(x)
                        if isinstance(x, str) and x.replace(".", "", 1).isdigit()
                        else x
                    ),
                )
            except ValueError:
                return sorted(values)

        # If no levels are provided, infer levels
        if colname_levels is None:
            current_values = adata.obs[colname]
            # Convert to Series for compatibility with .apply()
            if (
                pd.api.types.is_numeric_dtype(current_values)
                or pd.Series(current_values.astype(str))
                .apply(lambda x: x.replace(".", "", 1).isdigit())
                .all()
            ):
                levels = sort_char_numbers(current_values.unique())
            else:
                levels = (
                    current_values.value_counts().index.tolist()
                )  # Order by frequency
        else:
            # Parse provided levels
            levels = [x.strip() for x in colname_levels.split(",") if x.strip()]

            # Check if all current values are in the provided levels
            current_values = adata.obs[colname].unique()
            if not set(current_values).issubset(set(levels)):
                raise ValueError(
                    "Some values in the column are not in the provided levels."
                )

            # Ensure levels cover all unique values in the column
            if len(current_values) > len(levels):
                raise ValueError(
                    "The number of provided levels must be greater than or equal to the unique values in the column."
                )

        # Convert the column to a categorical type with the specified levels
        adata = adata.copy()
        adata.obs[colname] = pd.Categorical(
            adata.obs[colname], categories=levels, ordered=True
        )

        return adata


def standardize_delimited_vector(input_vec, delimiter=","):
    """
    标准化分隔符分隔的输入向量

    >>> standardize_delimited_vector(" C0 , !C1 ", ",")
    ['C0', '!C1']

    >>> standardize_delimited_vector(["C0", "C1"])
    ['C0', 'C1']
    """
    if isinstance(input_vec, str):
        # 分割字符串并清理空白
        processed = [x.strip() for x in input_vec.split(delimiter) if x.strip()]
        return processed
    elif isinstance(input_vec, (list, tuple)):
        # 直接返回列表/元组
        return list(input_vec)
    else:
        # 处理其他类型输入（如None）
        return input_vec


def get_predicate_vars(
    predicate: Union[str, List[str], None] = None,
) -> Optional[List[str]]:
    """
    Extract variable names from predicate expressions containing 'in' clauses.

    Args:
        predicate: Either a string with 'in' expressions, a list of variable names,
                  or None. If None or "all", returns None.

    Returns:
        List of unique variable names found, or None if predicate is None/"all".
        Empty list if no matches found.

    Raises:
        ValueError: If predicate is not str, list, or None
    """
    if predicate is None or predicate == "all":
        return None

    if isinstance(predicate, list):
        return list(set(predicate))  # Return unique variables if input is list

    if not isinstance(predicate, str):
        raise ValueError("Predicate must be string, list of variables, or None")

    # Find all matches of the pattern: word characters (including . and _) followed by ' in '
    matches = re.findall(r"\b([\w.]+)\s+in\s+", predicate)

    variables = list(set(matches))  # Get unique variables

    if not variables:
        warnings.warn("No variables found matching the pattern")
        return []

    return variables


def get_h5ad(
    input: Union[str, Path, sc.AnnData], validate_normalize=False
) -> sc.AnnData:
    """
    Read h5ad file or return AnnData object directly

    Parameters:
    -----------
    input : Union[str, Path, sc.AnnData]
        Can be either:
        - Path to .h5ad file
        - Existing AnnData object

    Returns:
    --------
    sc.AnnData
        Loaded or passed-through AnnData object

    Raises:
    -------
    ValueError
        If input file is not .h5ad format
    TypeError
        If input is neither path nor AnnData object

    Examples:
    --------
    >>> adata = get_h5ad("data.h5ad")  # Read from file
    >>> adata = get_h5ad(adata_object)  # Pass through existing object
    """
    if isinstance(input, (str, Path)):
        if str(input).endswith(".h5ad"):
            adata = sc.read_h5ad(input)
            if validate_normalize:
                adata = validate_and_normalize_adata(adata)
            return adata
        elif str(input).endswith(".h5"):
            adata = sc.read_10x_h5(input)
            return adata
        else:
            raise ValueError("Input file must be .h5ad, format")
    elif isinstance(input, sc.AnnData):
        if validate_normalize:
            input = validate_and_normalize_adata(input)
        return input
    else:
        raise TypeError("Input must be either h5ad file path or AnnData object")


def get_downsample_adata(
    adata: AnnData,
    group_by: str = None,
    target_cells: str = None,
    error_threshold: int = 1000,
    verbose: bool = True,
) -> AnnData:
    if target_cells == None:
        return adata
    else:
        target_cells = int(target_cells)

    np.random.seed(2025)
    if not isinstance(adata, AnnData):
        raise ValueError("Input 'adata' must be an AnnData object.")
    if target_cells <= 0:
        raise ValueError("'target_cells' must be a positive integer.")
    if error_threshold < 0:
        raise ValueError("'error_threshold' must be a non-negative number.")

    if group_by is None:
        group_by = "cluster"

    if group_by not in adata.obs.columns:
        raise ValueError(f"'group_by' must be a valid column name in adata.obs.")

    if adata.n_obs <= (target_cells + error_threshold):
        LOGGER.info(
            f"Dataset already meets the target cell count: {adata.n_obs} cells."
        )
        return adata

    # Create a combo column for grouping
    adata.obs["combo"] = adata.obs[group_by].astype(str)

    group_counts = adata.obs["combo"].value_counts()
    upper_bound = group_counts.max()

    low = 3
    high = upper_bound
    history = []

    max_iter = target_cells // 1000
    for i in range(max_iter):
        if high - low <= 2:
            if verbose:
                LOGGER.info("Reached minimum search interval")
            break

        current_downsample = (low + high) // 2

        np.random.seed(1236)
        sampled_indices = (
            adata.obs.groupby("combo")
            .apply(
                lambda x: x.sample(n=min(len(x), current_downsample), random_state=1236)
            )
            .index.get_level_values(1)
        )

        sub_adata = adata[sampled_indices]

        actual_cells = sub_adata.n_obs
        error = abs(actual_cells - target_cells)

        history.append(
            {
                "iteration": i + 1,
                "downsample_numb": current_downsample,
                "actual_cells": actual_cells,
                "error": error,
            }
        )

        if error <= error_threshold:
            if verbose:
                LOGGER.info(f"Converged at iteration {i + 1}")
            break

        if actual_cells < target_cells:
            low = current_downsample
        else:
            high = current_downsample

    if verbose:
        LOGGER.info(f"{actual_cells} cells in downsampled AnnData")

    history_df = pd.DataFrame(history)
    best_index = history_df["error"].idxmin()

    LOGGER.info("Iteration history:")
    LOGGER.info(history_df)
    LOGGER.info(
        f"Final downsample value: {history_df.loc[best_index, 'downsample_numb']}"
    )
    LOGGER.info(f"Final error: {history_df.loc[best_index, 'error']}")

    sampled_indices = (
        adata.obs.groupby("combo")
        .apply(
            lambda x: x.sample(
                n=min(len(x), history_df.loc[best_index, "downsample_numb"]),
                random_state=1236,
            )
        )
        .index.get_level_values(1)
    )
    del adata.obs["combo"]
    return adata[sampled_indices]


def loadH5AD(
    input: Union[str, Path, sc.AnnData],
    sample: Union[str, list] = None,
    group: Union[str, list] = None,
    cluster: Union[str, list] = None,
    celltype: Union[str, list] = None,
    predicate: Union[str, list] = None,
    downsample: int = None,
    groupby: str = None,
    groupby_levels: Union[str, list] = None,
    palette: str = None,
) -> sc.AnnData:
    """
    Load h5ad file and filter cells based on provided criteria.
    Parameters:
    -----------
    input : Union[str, Path, sc.AnnData]
        Can be either:
        - Path to .h5ad file
        - Existing AnnData object
    sample : Union[str, list], optional
        Sample ID(s) to filter cells by. If None, all samples are included.
    group : Union[str, list], optional
        Group(s) to filter cells by. If None, all groups are included.
    cluster : Union[str, list], optional
        Cluster(s) to filter cells by. If None, all cluster are included.
    celltype : Union[str, list], optional
        New cell type(s) to filter cells by. If None, all cell types are included.
    predicate : Union[str, list], optional
        Predicate(s) to filter cells by. If None, all predicates are included.predicate="(sample in ['STB1', 'STB4']) and (cluster in ['1', '5', '6'])"
    downsample : int, optional
        Number of cells to downsample to. If None, no downsampling is performed.
    groupby : str, optional
        Column name in adata.obs to group by for filtering. If None, no grouping is performed.
    rm_groupby_cells_lessthan : Union[int, list], optional
        Minimum number of cells required in each group to keep the group. If None, no filtering is performed.
    palette : str, optional
        Color palette for plotting. If None, default palette is used.
    Returns:
    --------
    sc.AnnData
        Filtered AnnData object
    """
    # Load the h5ad file or return the AnnData object directly
    adata = get_h5ad(input)
    use_col = [
        "sample",
        "group",
        "cluster",
        "celltype",
        groupby,
    ]

    variables = get_predicate_vars(predicate)
    if isinstance(variables, list):
        use_col.extend(variables)
    else:
        use_col.append(variables)
    use_col = list(set([col for col in use_col if col is not None]))
    # subset the AnnData object
    adata = subsetH5AD(
        adata,
        sample=sample,
        group=group,
        cluster=cluster,
        celltype=celltype,
        predicate=predicate,
    )
    if groupby_levels != None:
        adata = standardize_levels(
            adata,
            colname=groupby,
            colname_levels=groupby_levels,
        )

    # standardize_colors
    adata = standardize_colors(adata, ["sample", "group", "cluster", "celltype"])
    adata = standardize_colors(adata, [groupby], [palette])

    adata = get_downsample_adata(
        adata,
        group_by=groupby,
        target_cells=downsample,
        error_threshold=1000,
        verbose=True,
    )

    for colname in use_col:
        if colname in adata.obs.columns:
            adata = standardize_levels(adata, colname=colname)

    return adata
