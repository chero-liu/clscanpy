import scanpy as sc
import numpy as np
from typing import Tuple, Optional, List, Union
from pathlib import Path
from oescanpy.tools.read.utils import (
    validate_and_normalize_adata,
    convert_seurat_to_anndata,
)
from anndata import AnnData
import pandas as pd
import warnings
from oescanpy.tools.color.colors import add_color_to_dataframe
import re
import os
import logging

LOGGER = logging.getLogger(__name__)


def update_adata_from_metadata(
    adata: AnnData,
    metadata_path: str,
    id_col: str = "rawbc",
) -> AnnData:
    """
    Update AnnData object with metadata from a CSV/TSV file.
    Parameters:
    - adata: AnnData object to update
    - metadata_path: Path to metadata file (CSV or TSV)
    - id_col: Column name used for cell matching (default: "rawbc")

    Returns:
    Updated AnnData object with additional columns from metadata

    Raises:
    - FileNotFoundError: If metadata file does not exist
    - ValueError: If ID column not found
    """

    if not isinstance(adata, AnnData):
        raise TypeError("adata must be an AnnData object")
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    if not isinstance(id_col, str) or not id_col:
        raise ValueError("id_col must be a non-empty string")

    ext = os.path.splitext(metadata_path)[1].lower()
    sep = "\t" if ext == ".tsv" else ","

    try:
        subdata = pd.read_csv(metadata_path, sep=sep, dtype={id_col: str}, engine="c")
        LOGGER.info(
            f"Loaded metadata: {subdata.shape[0]} rows, {subdata.shape[1]} columns"
        )
    except Exception as e:
        LOGGER.exception(f"Error reading metadata file: {metadata_path}")
        raise

    original_obs_cols = len(adata.obs.columns)
    LOGGER.debug(f"Original adata.obs columns: {original_obs_cols}")

    if subdata[id_col].duplicated().any():
        duplicate_ids = subdata[subdata.duplicated(subset=id_col, keep=False)][
            id_col
        ].unique()
        LOGGER.warning(
            f"Found {len(duplicate_ids)} duplicate {id_col} values in metadata. "
            "Keeping first occurrence and dropping duplicates."
        )
        subdata = subdata.drop_duplicates(subset=id_col, keep="first")
        LOGGER.info(
            f"Reduced metadata to {subdata.shape[0]} rows after removing duplicates"
        )

    if id_col not in subdata.columns:
        err_msg = f"ID column '{id_col}' not found in metadata columns: {subdata.columns.tolist()}"
        LOGGER.error(err_msg)
        raise ValueError(err_msg)

    if id_col not in adata.obs.columns:
        err_msg = f"ID column '{id_col}' not found in adata.obs columns: {adata.obs.columns.tolist()}"
        LOGGER.error(err_msg)
        raise ValueError(err_msg)

    adata.obs[id_col] = adata.obs[id_col].astype(str).str.strip()
    subdata[id_col] = subdata[id_col].astype(str).str.strip()

    original_size = adata.n_obs
    valid_ids = set(subdata[id_col])
    cell_mask = adata.obs[id_col].isin(valid_ids)

    if cell_mask.sum() == 0:
        sample_ids = adata.obs[id_col].head(5).tolist()
        err_msg = (
            f"No matching IDs found. Metadata has {len(valid_ids)} unique IDs, "
            f"first 5 AnnData IDs: {sample_ids}"
        )
        LOGGER.error(err_msg)
        raise ValueError(err_msg)

    adata = adata[cell_mask].copy()
    kept = adata.n_obs
    LOGGER.info(
        f"Filtered cells using '{id_col}'. Kept {kept}/{original_size} cells "
        f"({kept/original_size:.1%})"
    )

    metadata_non_id_cols = [col for col in subdata.columns if col != id_col]

    new_cols = [col for col in metadata_non_id_cols if col not in adata.obs.columns]

    update_cols = [col for col in metadata_non_id_cols if col in adata.obs.columns]

    all_cols = update_cols + new_cols

    if not all_cols:
        LOGGER.info("No columns to update or add")
        return adata

    expected_final_cols = original_obs_cols + len(new_cols)
    LOGGER.info(
        f"Will add {len(new_cols)} new columns and update {len(update_cols)} existing columns. "
        f"Expected final columns: {expected_final_cols}"
    )

    mapping_df = subdata.set_index(id_col)[all_cols].copy()

    for col in all_cols:
        col_series = mapping_df[col].astype(str).str.strip()

        if col_series.index.duplicated().any():
            LOGGER.warning(
                f"Column '{col}' has duplicate IDs. Keeping first occurrence."
            )
            col_series = col_series[~col_series.index.duplicated(keep="first")]

        adata.obs[col] = adata.obs[id_col].map(col_series)

        na_count = adata.obs[col].isna().sum()
        if na_count > 0:
            LOGGER.warning(
                f"Column '{col}' has {na_count} missing values after mapping"
            )

    final_cols = len(adata.obs.columns)
    if final_cols != expected_final_cols:
        LOGGER.warning(
            f"Unexpected column count: Expected {expected_final_cols}, Actual {final_cols}"
        )

    if new_cols:
        LOGGER.info(f"Added new columns to adata.obs: {sorted(new_cols)}")

    LOGGER.info(
        f"Successfully updated/added {len(all_cols)} columns. "
        f"Final adata.obs columns: {final_cols} (original: {original_obs_cols})"
    )

    return adata


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
        "sampleid": "cellpaper_custom",
        "clusters": "customecol2_light",
        "new_celltype": "tableau20",
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
    sampleid=None,
    group=None,
    clusters=None,
    new_celltype=None,
    predicate=None,
    metadata=None,
):
    """
    使用标准化分隔符处理函数增强的AnnData子集筛选

    Parameters:
        adata : AnnData
            输入的AnnData对象
        sampleid : str (可选)
            逗号分隔的样本ID，支持!前缀排除，如 "S1,!S2"
        group : str (可选)
            逗号分隔的分组名称，支持!前缀排除
        clusters : str (可选)
            逗号分隔的细胞簇ID，支持!前缀排除
        new_celltype : str (可选)
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
        ("new_celltype", new_celltype),
        ("sampleid", sampleid),
        ("group", group),
        ("clusters", clusters),
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

    if metadata:
        # Determine separator based on file extension
        adata = update_adata_from_metadata(
            adata,
            metadata_path=metadata,
            id_col="rawbc",
        )

    # Standardize levels for relevant columns
    use_col = [
        "sampleid",
        "group",
        "clusters",
        "new_celltype",
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
    if isinstance(adata.obs[colname].dtype, pd.CategoricalDtype) and colname_levels is None:
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
                levels = current_values.value_counts().index.tolist()  # Order by frequency
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
        elif input.endswith(".h5seurat") or input.endswith(".rds"):
            adata = convert_seurat_to_anndata(input, use_raw_counts=False)
            return adata
        else:
            raise ValueError("Input file must be .h5ad, rds, h5seurat format")
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
        group_by = "clusters"

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
    sampleid: Union[str, list] = None,
    group: Union[str, list] = None,
    clusters: Union[str, list] = None,
    new_celltype: Union[str, list] = None,
    predicate: Union[str, list] = None,
    metadata: Union[str] = None,
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
    sampleid : Union[str, list], optional
        Sample ID(s) to filter cells by. If None, all samples are included.
    group : Union[str, list], optional
        Group(s) to filter cells by. If None, all groups are included.
    clusters : Union[str, list], optional
        Cluster(s) to filter cells by. If None, all clusters are included.
    new_celltype : Union[str, list], optional
        New cell type(s) to filter cells by. If None, all cell types are included.
    predicate : Union[str, list], optional
        Predicate(s) to filter cells by. If None, all predicates are included.predicate="(sampleid in ['STB1', 'STB4']) and (clusters in ['1', '5', '6'])"
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
        "sampleid",
        "group",
        "clusters",
        "new_celltype",
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
        sampleid=sampleid,
        group=group,
        clusters=clusters,
        new_celltype=new_celltype,
        predicate=predicate,
        metadata=metadata,
    )
    if groupby_levels != None:
        adata = standardize_levels(
            adata,
            colname=groupby,
            colname_levels=groupby_levels,
        )

    # standardize_colors
    adata = standardize_colors(adata, ["sampleid", "group", "clusters", "new_celltype"])
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
