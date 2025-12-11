import os
import pandas as pd
from anndata import AnnData
from typing import Dict, Any, Literal, Tuple, List
import json
import importlib
from clscanpy.config import ROOT_PATH
import matplotlib.pyplot as plt
import abc
import numpy as np
import scanpy as sc
import logging
import sys
import re
import pathlib
import warnings

import logging

LOGGER = logging.getLogger(__name__)


def find_assay_init(assay):
    init_module = importlib.import_module(f"clscanpy.{assay}.__init__")
    return init_module


def clean_adata(adata, additional_patterns=None):
    """Remove unnecessary data columns from adata.obs

    Args:
        adata (anndata.AnnData): AnnData object to clean
        additional_patterns (list, optional): Additional patterns to remove.
            These will be combined with the default patterns. Defaults to None.

    Returns:
        anndata.AnnData: Cleaned AnnData object
    """
    # Default patterns to remove
    patterns_to_remove = [
        "n_genes",
        "n_counts",
        "doubletdoublet",
        "doublet",
        "total_counts_mt",
        "total_counts_hb",
        "log1p_n_genes_by_counts",
        "log1p_total_counts",
        "pct_counts_in_top_50_genes",
        "pct_counts_in_top_100_genes",
        "pct_counts_in_top_200_genes",
        "pct_counts_in_top_500_genes",
        "log1p_total_counts_hb",
        "datatype",
    ]

    # Add additional patterns if provided
    if additional_patterns is not None:
        # Ensure additional_patterns is a list
        if not isinstance(additional_patterns, list):
            raise ValueError("additional_patterns must be a list of strings")
        # Combine patterns, avoiding duplicates
        patterns_to_remove = list(set(patterns_to_remove + additional_patterns))

    # Get existing columns in adata.obs
    existing_cols = adata.obs.columns.tolist()

    # Find columns matching patterns (case-insensitive)
    cols_to_remove = set()
    for col in existing_cols:
        if any(pattern.lower() in col.lower() for pattern in patterns_to_remove):
            cols_to_remove.add(col)

    # Remove matched columns
    if cols_to_remove:
        LOGGER.info(
            f"Removing {len(cols_to_remove)} columns: {', '.join(cols_to_remove)}"
        )
        adata.obs = adata.obs.drop(columns=list(cols_to_remove))
    else:
        LOGGER.info("No columns found for removal")

    return adata


def find_step_module(assay, step):
    file_path_dict = {
        "assay": f"{ROOT_PATH}/{assay}/{step}.py",
        "tools": f"{ROOT_PATH}/tools/{step}.py",
    }
    init_module = find_assay_init(assay)
    if os.path.exists(file_path_dict["assay"]):
        step_module = importlib.import_module(f"clscanpy.{assay}.{step}")
    elif hasattr(init_module, "IMPORT_DICT") and step in init_module.IMPORT_DICT:
        module_path = init_module.IMPORT_DICT[step]
        step_module = importlib.import_module(f"{module_path}.{step}")
    elif os.path.exists(file_path_dict["tools"]):
        step_module = importlib.import_module(f"clscanpy.tools.{step}")
    else:
        raise ModuleNotFoundError(f"No module found for {assay}.{step}")

    return step_module


def check_file_exist(file_path):
    if not os.path.exists(file_path):
        return False
    return True


def check_mkdir(path):
    """
    Ensure that the specified directory exists. If it doesn't exist, create it.

    Parameters:
    - path (str): The path of the directory to be checked/created.

    Returns:
    - bool: True if the directory exists or was successfully created, False otherwise.
    """
    try:
        # Check if the directory exists
        if os.path.exists(path):
            return True

        # Create the directory if it doesn't exist
        os.makedirs(path)
        LOGGER.info(f"Directory '{path}' created successfully.")
        return True

    except Exception as e:
        LOGGER.info(f"Error: {e}")
        return False


def save_parameters_to_anndata(adata: AnnData, parameters: Dict[str, Any]) -> AnnData:
    """
    Save parameters to the 'uns' attribute of an AnnData object.

    Parameters:
    - adata (AnnData): The AnnData object where parameters will be saved.
    - parameters (Dict[str, Any]): A dictionary containing parameters to be saved.

    Returns:
    AnnData: Updated AnnData object with saved parameters.
    """
    if "preprocess_para" in adata.uns:
        para = json.loads(adata.uns["preprocess_para"])
    else:
        para = {}
    for key, value in parameters.items():
        para[key] = value
    adata.uns["preprocess_para"] = json.dumps(para)

    return adata


def retrieve_layers(
    adata: AnnData,
    data: AnnData,
    layer_keys: Literal["raw", "normalised"] = ["raw", "normalised"],
    use_layer: str = "normalised",
) -> AnnData:
    """
    Retrieve specified layers from 'data' and update 'adata' object.

    Parameters:
    - adata (AnnData): The AnnData object to be updated with retrieved layers.
    - data (AnnData): The AnnData object containing layers to be retrieved.
    - layer_keys (Literal["raw", "normalised"]): List of layer keys to be retrieved. Default is ["raw", "normalised"].
    - use_layer (str): The layer to be used for 'adata.X'. Default is "normalised".

    Returns:
    AnnData: Updated 'adata' object.
    """
    # Create a new AnnData object for raw data
    raw_adata = adata.raw.to_adata()

    # Set raw data to itself
    raw_adata.raw = raw_adata

    # Update layers with specified data
    for layer in layer_keys:
        raw_adata.layers[layer] = data.layers[layer]

    # Update the original 'adata' object with the preprocessed data
    raw_adata.X = data.layers[use_layer]

    return raw_adata


def check_adata(adata: AnnData) -> AnnData:
    """
    Update `adata.X` based on the existence of a valid matrix named `raw` in `adata.layers`.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.

    Returns
    -------
    AnnData
    """
    if "raw" in adata.layers:
        adata.X = adata.layers["raw"].copy()
        return adata
    else:
        LOGGER.info("adata has no valid matrix raw")


def add_counts_to_data(data: pd.DataFrame, column_name: str) -> pd.DataFrame:
    """
    统计指定列中每个类别对应的数量，并将结果添加为新列。

    Parameters:
        data (pd.DataFrame): 输入的 Pandas 数据框。
        column_name (str): 要统计的列名。

    Returns:
        pd.DataFrame: 添加了新列的 Pandas 数据框，新列名为 `<column_name>_counts`。
    """
    if column_name not in data.columns:
        raise ValueError(f"列名 '{column_name}' 不存在于数据框中")

    counts = data[column_name].value_counts()

    count_column_name = f"{column_name}_counts"
    data[count_column_name] = data[column_name].map(counts)

    return data


def save_figure(
    save_path: str, save_formats: Tuple[str, ...] = ("png", "pdf"), dpi: int = 300
) -> None:
    """Save matplotlib figure to specified path with multiple formats.

    Args:
        save_path: Base path without extension
        save_formats: Tuple of file formats to save (e.g. ("png", "pdf"))
        dpi: Image resolution
    """
    if not save_path:
        return

    check_mkdir(os.path.dirname(save_path))
    for fmt in save_formats:
        plt.savefig(f"{save_path}.{fmt}", dpi=dpi, bbox_inches="tight")


def handle_figure_display(show: bool = True) -> None:
    """Handle figure display or closing.

    Args:
        show: Whether to display the figure (if False, closes it)
    """
    if show:
        plt.show()
    else:
        plt.close()


def s_common(parser):
    """subparser common arguments"""
    # parser.add_argument("--prefix", help="Prefix of all output files.", required=True)
    parser.add_argument("--outdir", help="Output directory.", required=True)
    parser.add_argument("--species", help="Species.", required=False)
    parser.add_argument("--sampleid", help="Sample ID.", default=None, required=False)
    parser.add_argument("--group", help="Group name.", default=None, required=False)
    parser.add_argument("--clusters", help="Clusters.", default=None, required=False)
    parser.add_argument("--new_celltype", help="New cell type.", required=False)
    parser.add_argument("--predicate", help="Predicate.", required=False)
    parser.add_argument("--thread", help="Number of threads.", default=10, type=int)
    parser.add_argument("--downsample", help="downsample", default=None, required=False)
    parser.add_argument(
        "--palette", help="palette", default="customecol2", required=False
    )
    parser.add_argument(
        "--metadata",
        help="Metadata file, if not provided, will use adata.obs.",
        default=None,
        required=False,
    )
    return parser


class Step:
    """
    Step class
    """

    def __init__(self, args):
        self.args = args
        self.outdir = args.outdir
        self.species = args.species
        self.assay = args.subparser_assay
        self.downsample = args.downsample
        self.palette = args.palette
        self.sampleid = args.sampleid
        self.group = args.group
        self.clusters = args.clusters
        self.new_celltype = args.new_celltype
        self.predicate = args.predicate
        self.metadata = args.metadata
        self.thread = int(args.thread)

        check_mkdir(self.outdir)

    @abc.abstractmethod
    def run(self):
        sys.exit("Please implement run() method.")

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        LOGGER.info("bey")


def setup_logger():
    """配置日志记录器"""
    logger = logging.getLogger("celltypist_annotate")

    # 如果logger已经有处理程序，先清除所有已有的处理程序
    if logger.handlers:
        logger.handlers.clear()

    logger.setLevel(logging.INFO)

    # 创建控制台处理程序
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)

    # 创建格式化器
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    console_handler.setFormatter(formatter)

    # 添加处理程序到logger
    logger.addHandler(console_handler)

    # 防止日志传递到父logger
    logger.propagate = False

    return logger


def is_valid_path(path_str: str) -> bool:
    """Determine if a string represents a valid filesystem path

    Parameters:
        path_str: Input string to validate

    Returns:
        bool: True if valid path, False otherwise
    """
    # 1. Basic string validation
    if not isinstance(path_str, str) or not path_str:
        return False

    # 2. Path format validation
    try:
        # Create path object (doesn't verify existence)
        path = pathlib.Path(path_str)

        # Windows-specific checks
        if os.name == "nt":
            # Validate Windows path format (C:\ or \\server\share)
            if ":" in path_str and not path_str.strip().endswith(":"):
                drive, rest = path_str.split(":", 1)
                if len(drive) != 1 or not drive.isalpha():
                    return False
        else:
            # Unix-like systems must start with root
            if not path_str.startswith("/"):
                return False

        # 3. Component validation
        for part in path.parts:
            # Check for illegal characters
            if any(char in part for char in '<>:"|?*\\\0'):
                return False

            # Check Windows reserved names
            if os.name == "nt":
                if part.upper() in [
                    "CON",
                    "PRN",
                    "AUX",
                    "NUL",
                    "COM1",
                    "COM2",
                    "COM3",
                    "COM4",
                    "COM5",
                    "COM6",
                    "COM7",
                    "COM8",
                    "COM9",
                    "LPT1",
                    "LPT2",
                    "LPT3",
                    "LPT4",
                    "LPT5",
                    "LPT6",
                    "LPT7",
                    "LPT8",
                    "LPT9",
                ]:
                    return False

        # 4. Path normalization test
        normalized = os.path.normpath(path_str)

        # 5. OS path operation test
        dirname, basename = os.path.split(path_str)

        # All validation passed
        return True

    except (ValueError, TypeError, OSError) as e:
        # Catch any path parsing exceptions
        return False


from typing import Dict, Any, Literal, Tuple, List


def get_group_order(
    df: pd.DataFrame,
    key_col: str,
    value_col: str,
) -> Dict[str, str]:
    """
    获取两列之间的精确映射关系

    参数:
        df: 输入数据框
        key_col: 键列名 (如clusters)
        value_col: 值列名 (如clusters_col)
    """
    missing_cols = [c for c in [key_col, value_col] if c not in df.columns]
    if missing_cols:
        raise ValueError(f"以下列不存在: {', '.join(missing_cols)}")

    unique_pairs = df[[key_col, value_col]].drop_duplicates()

    unique_pairs = unique_pairs.astype(str)

    return unique_pairs


def save_dict_to_txt(data_dict, file_path):
    """
    将字典保存为简单的键值对文本格式

    参数:
        data_dict: 要保存的字典
        file_path: 要保存的文件路径
    """
    with open(file_path, "w") as f:
        for key, value in data_dict.items():
            f.write(f"{key}={value}\n")
