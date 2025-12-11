from typing import Optional, Union, List
import os
import pandas as pd
import numpy as np
from pathlib import Path
import scanpy as sc
from scipy import sparse
import sys
import gzip
import re
from random import sample
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import anndata
from oescanpy.tools.preprocessing.basicpreprocessing import normalize
from typing import List
import logging

LOGGER = logging.getLogger(__name__)


class MultiSpeciesSplitter:
    """
    Split mixed-species single-cell data (e.g., human-mouse)
    for a single CellRanger output sample.

    Expected structure:
    /path/to/sample/outs/
        ├── filtered_feature_bc_matrix/   <-- input_dir
        └── analysis/gem_classification.csv
    """

    def __init__(self, input_dir: str, mix_species: str, gene_prefix_clean: bool = True):
        """
        Parameters
        ----------
        input_dir : str
            Path to the 'filtered_feature_bc_matrix' directory.
        mix_species : str
            Target species to extract (e.g., "GRCh38", "mm10", "Rat").
        gene_prefix_clean : bool, default=True
            Whether to remove the prefix like "GRCh38-" from gene names.
        """
        self.input_dir = os.path.abspath(input_dir)
        self.mix_species = mix_species
        self.gene_prefix_clean = gene_prefix_clean

        self.sample_name = os.path.basename(os.path.dirname(os.path.dirname(self.input_dir)))
        self.logger = LOGGER

    def _read_gem_classification(self) -> pd.DataFrame:
        """Read gem_classification.csv located under outs/analysis."""
        gem_path = os.path.join(os.path.dirname(self.input_dir), "analysis", "gem_classification.csv")
        if not os.path.exists(gem_path):
            raise FileNotFoundError(f"Cannot find gem_classification.csv at: {gem_path}")

        data = pd.read_csv(gem_path)
        if "barcode" not in data.columns or "call" not in data.columns:
            raise ValueError(f"Missing required columns in {gem_path}: {list(data.columns)}")

        data["cells"] = self.sample_name + "-" + data["barcode"].astype(str)
        data["cells"] = data["cells"].str.replace("-1", "", regex=False)
        self.logger.info(f"Loaded {len(data)} barcode records from {gem_path}")
        return data

    def _filter_barcodes(self, key_data: pd.DataFrame) -> List[str]:
        """Filter barcodes belonging to the target species."""
        filtered = key_data.loc[key_data["call"] == self.mix_species, "cells"].tolist()
        self.logger.info(f"Filtered {len(filtered)} barcodes for species '{self.mix_species}'")
        return filtered

    def _process_adata(self, adata: sc.AnnData) -> sc.AnnData:
        """Subset genes by species prefix and clean their names."""
        safe_prefix = re.escape(self.mix_species)
        genes_of_interest = [g for g in adata.var_names if re.search(safe_prefix, g)]
        self.logger.info(f"Found {len(genes_of_interest)} genes with prefix '{self.mix_species}'")

        sub_adata = adata[:, genes_of_interest].copy()

        if self.gene_prefix_clean:
            cleaned = [re.sub(f"^{safe_prefix}-+", "", g) for g in sub_adata.var_names]
            sub_adata.var_names = cleaned

        return sub_adata

    def split(self, adata: sc.AnnData) -> sc.AnnData:
        """Main entry: subset an AnnData object by target species."""
        self.logger.info(f"Starting split for species '{self.mix_species}' in sample '{self.sample_name}'")

        key_data = self._read_gem_classification()
        barcodes = self._filter_barcodes(key_data)

        sub_adata = adata[adata.obs_names.isin(barcodes)].copy()
        self.logger.info(f"Subset AnnData: {sub_adata.n_obs} cells retained")

        sub_adata = self._process_adata(sub_adata)
        self.logger.info(f"Final AnnData shape: {sub_adata.n_obs} cells × {sub_adata.n_vars} genes")

        return sub_adata

def getSingleDataFormat(datatype, filename):
    if not datatype:
        if os.path.isdir(filename):
            return "mtx-10x"
        elif filename.endswith(".h5ad"):
            return "h5ad-10x"
        else:
            return "normal-10x"
    else:
        return datatype


def reset_barcode(
    raw_barcodes: List[str],
    barcode_pool_file: Path = Path("/data/software/cellranger/cellranger-7.0.1/lib/python/cellranger/barcodes/3M-february-2018.txt.gz")
) -> List[str]:
    if barcode_pool_file.is_file():
        LOGGER.info(f"Barcode格式转换: 目标barcode抽样自{barcode_pool_file}。")
        # 顺序读取 .gz file
        with gzip.open(barcode_pool_file, 'rt') as f:
            barcode_pool = [i.strip() for i in f.readlines()]
    else:
        raise FileNotFoundError(
            f"Barcode pool file {barcode_pool_file} does not exist."
        )
    new_barcodes = sample(barcode_pool, len(raw_barcodes))
    if re.match(r"-[0-9]+$", raw_barcodes[0]):
        suffixes = [re.search(r"-([0-9]+)$", bc).group(1) for bc in raw_barcodes]
        new_barcodes = [f"{bc}-{suffix}" for bc, suffix in zip(new_barcodes, suffixes)]
    else:
        new_barcodes = [f"{bc}-1" for bc in new_barcodes]
    return new_barcodes


def read(
    filename: Union[str, Path],
    var_names: Optional[str] = "gene_symbols",
    make_unique: bool = True,
    cache: bool = False,
    delimiter: str = "\t",
    dtype: str = None,
    prefix: str = None,
    index: Optional[int] = None,
    barcode_pool_file: Union[str, Path] = Path("/data/software/cellranger/cellranger-7.0.1/lib/python/cellranger/barcodes/3M-february-2018.txt.gz")
):
    if dtype.startswith("mtx"):
        adata = sc.read_10x_mtx(
            filename, var_names=var_names, cache=cache, make_unique=make_unique
        )
    elif dtype.startswith("h5ad"):
        adata = sc.read(filename)
        if "count" in adata.layers:
            adata.X = adata.layers["count"]
        else:
            sys.exit("No raw counts matrix can be found in Anndata object")   
    elif dtype.startswith("normal"):
        try:
            adata = sc.read_csv(filename, delimiter=delimiter).T
        except ValueError:
            adata = sc.read_csv(filename, delimiter=delimiter, first_column_names=True)
    else:
        raise ValueError("Not a valid matrix format, plz check.")

    gex_rows = list(map(lambda x: x.replace("_", "-"), adata.var.index))
    adata.var.index = gex_rows

    if index:
        if dtype.endswith("mobi") or dtype.endswith("C4") or dtype.endswith("BD"):
            rawbc_cols = reset_barcode(
                adata.obs_names,
                Path(barcode_pool_file),
            )
            rawbc_cols = list(map(lambda x: x.replace("-1", f"-{index}"), rawbc_cols))

            old_rawbc_cols = list(
                map(lambda x: x.replace("-1", f"-{index}"), adata.obs.index)
            )
            adata.obs["old_rawbc"] = old_rawbc_cols

        else:
            rawbc_cols = list(
                map(lambda x: x.replace("-1", f"-{index}"), adata.obs.index)
            )
        adata.obs["rawbc"] = rawbc_cols
    if prefix:
        gex_cols = list(
            map(lambda x: prefix + "-" + x.replace("-1", ""), adata.obs.index)
        )
        adata.obs.index = gex_cols

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    return adata


def ensureSparseMatrix(adata):
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    adata.layers["raw"] = adata.X


def summarizeAdataMetrics(
    adata,
    metrics,
    group_by="sampleid",
    suffix="qc",
):
    """
    统计每个样本的细胞数、基因数，以及指定指标的最大值、最小值、均值和中位值。

    Parameters:
    -----------
    adata : AnnData
        AnnData 对象，包含单细胞数据。
    metrics : list
        要统计的指标列表，例如 ['nFeature_RNA', 'nCount_RNA', 'percent.mito','percent.HB']。
    group_by : str
        分组变量（通常是样本 ID），默认为 'sampleid'。

    Returns:
    --------
    pd.DataFrame
        每个样本的统计结果，包括细胞数、基因数和指定指标的统计值。
    """

    if group_by not in adata.obs.columns:
        raise ValueError(f"'{group_by}' 不在 adata.obs 中，请检查分组变量。")
    results = []
    for sample, group in adata.obs.groupby(group_by):
        cell_count = group.shape[0]
        stats = {"sampleid": sample, f"cell_count_{suffix}": cell_count}
        for metric in metrics:
            if metric in group.columns:
                stats[f"{metric}_max_{suffix}"] = group[metric].max()
                stats[f"{metric}_min_{suffix}"] = group[metric].min()
                stats[f"{metric}_mean_{suffix}"] = group[metric].mean()
                stats[f"{metric}_median_{suffix}"] = group[metric].median()
            else:
                stats[f"{metric}_max_{suffix}"] = None
                stats[f"{metric}_min_{suffix}"] = None
                stats[f"{metric}_mean_{suffix}"] = None
                stats[f"{metric}_median_{suffix}"] = None
        results.append(stats)
    return pd.DataFrame(results)


def convert_seurat_to_anndata(seurat_path, use_raw_counts=False):
    """
    将 Seurat 对象（.h5seurat文件）转换为 AnnData 对象

    Parameters:
    -----------
    seurat_path : str
        Seurat对象文件路径
    use_raw_counts : bool, default False
        如果为True，则adata.X存放原始计数矩阵
        如果为False，则adata.X存放标准化数据
    """
    # 激活自动转换
    pandas2ri.activate()

    # R代码：读取和转换数据
    r_code = """
    convert_seurat <- function(seurat_path) {
        tryCatch({
            source('/home/liuchenglong/script/lclFunc.r')
            # 加载必要的包
            suppressPackageStartupMessages({
                library(Seurat)
                library(SeuratObject)
                library(SeuratDisk)
            })
            
            # 读取Seurat对象
            seurat_obj = getRDS(seurat_path)
            # 获取所有可用的assay
            print("Available assays:")
            print(Assays(seurat_obj))
            
            # 获取默认assay
            default_assay <- DefaultAssay(seurat_obj)
            print(paste("Default assay:", default_assay))
            
            # 获取原始数据（counts）

            raw_counts <- GetAssayData(seurat_obj, slot = "counts", assay = default_assay)
            raw_counts <- as.matrix(raw_counts)

            
            # 获取标准化数据（data）

            norm_data <- GetAssayData(seurat_obj, slot = "data", assay = default_assay)
            norm_data <- as.matrix(norm_data)

            
            # 获取缩放数据（scale.data，如果存在）

            scale_data <- NULL
            scale_features <- NULL
            tryCatch({
                scale_data <- GetAssayData(seurat_obj, slot = "scale.data", assay = default_assay)
                scale_features <- rownames(scale_data)
                scale_data <- as.matrix(scale_data)
            }, error = function(e) {
                print("No scaled data found")
            })
            
            # 获取元数据
            metadata <- seurat_obj@meta.data
            
            # 获取降维结果
            print("Extracting dimensional reductions...")
            dimreduc_data <- list()
            for (reduc in Reductions(seurat_obj)) {
              if(reduc != "pca"){

                reduc_obj <- Embeddings(seurat_obj[[reduc]])
                dimreduc_data[[reduc]] <- as.matrix(reduc_obj)

            }
            }
            
            # 获取高变基因
            var_features <- VariableFeatures(seurat_obj)

            # 检查并获取空间坐标
            has_spatial <- FALSE
            spatial_coords <- NULL
            
            tryCatch({
                spatial_coords <- GetTissueCoordinates(seurat_obj)
                if (!is.null(spatial_coords) && nrow(spatial_coords) > 0) {
                    has_spatial <- TRUE

                }
            }, error = function(e) {
                print("No spatial coordinates found")
            })
            
            # 获取基因信息
            print("Extracting gene information...")
            gene_info <- rownames(raw_counts)
            gene_symbols <- NULL
            
            # 尝试获取基因符号
            tryCatch({
                if ("Symbol" %in% rownames(seurat_obj[[default_assay]]@meta.features)) {
                    gene_symbols <- seurat_obj[[default_assay]]@meta.features["Symbol",]

                } else if ("gene_name" %in% colnames(seurat_obj[[default_assay]]@meta.features)) {
                    gene_symbols <- seurat_obj[[default_assay]]@meta.features$gene_name

                }
            }, error = function(e) {
                print("No gene symbols found in meta.features")
            })
            
            # 获取metadata中因子信息
            fact_list = list()
            for (col in colnames(metadata)) {
                if (is.factor(metadata[[col]])) {
                    fact_list[[col]] = levels(metadata[[col]])} 
            }
            # 返回结果
            result <- list(
                raw_counts = raw_counts,
                norm_data = norm_data,
                scale_data = scale_data,
                scale_features = scale_features,
                metadata = metadata,
                dimreduc_data = dimreduc_data,
                var_features = var_features,
                has_spatial = has_spatial,
                gene_info = gene_info,  # 基因名称（可能是 ENSEMBL ID）
                gene_symbols = gene_symbols,  # 基因符号（如果可用）
                fact_list = fact_list # 因子信息
            )
            print("Data extraction completed.")
            print(names(result))
            if (has_spatial) {
                result$spatial_coords <- spatial_coords
            }
            
            return(result)
        }, error = function(e) {
            print(paste("Error in R:", e$message))
            stop(e$message)
        })
    }
    """
    # 运行R代码
    robjects.r(r_code)
    convert_func = robjects.globalenv["convert_seurat"]

    # 获取数据
    LOGGER.info("调用R函数转换数据...")
    r_data = convert_func(seurat_path)

    # 转换原始计数矩阵
    LOGGER.info("转换原始计数矩阵...")
    with localconverter(robjects.default_converter + pandas2ri.converter):
        raw_counts = pd.DataFrame(r_data[0])

        # 获取基因信息
        gene_info = list(r_data[8])  # 确保转换为Python列表
        LOGGER.info(f"基因信息示例: {gene_info[:5]}")

        # 显式设置索引
        raw_counts.index = gene_info
        LOGGER.info(f"原始计数矩阵维度: {raw_counts.shape}")

    # 转换标准化数据
    LOGGER.info("转换标准化数据...")
    with localconverter(robjects.default_converter + pandas2ri.converter):
        norm_data = pd.DataFrame(r_data[1])
        # 确保使用相同的基因名称
        norm_data.index = raw_counts.index
        LOGGER.info(f"标准化数据维度: {norm_data.shape}")

    # 创建 var DataFrame，确保保留基因名称
    var_df = pd.DataFrame(index=raw_counts.index)

    # 转换 metadata
    LOGGER.info("转换 metadata...")
    with localconverter(robjects.default_converter + pandas2ri.converter):
        try:
            metadata_py = pd.DataFrame(r_data[4])
            # 正确获取因子信息并转换为字典
            fact_dict = dict(r_data[10].items()) if r_data[10] is not None else {}
            # 处理每个因子列
            for col, levels in fact_dict.items():
                if col in metadata_py.columns:
                    # 确保levels是Python列表
                    levels = list(levels) if isinstance(levels, robjects.vectors.ListVector) else levels
                    # 转换为有序分类
                    metadata_py[col] = pd.Categorical(
                        metadata_py[col], 
                        categories=levels, 
                        ordered=True
                    )
        except Exception as e:
            LOGGER.warning(f"转换 metadata 时出错: {str(e)}")
            metadata_py = pd.DataFrame(
                index=norm_data.columns
            )  # 创建空的 metadata DataFrame

    # 转换降维结果
    LOGGER.info("转换降维结果...")
    obsm = {}
    with localconverter(robjects.default_converter + pandas2ri.converter):
        dimreduc_data = r_data[5]
        if isinstance(dimreduc_data, dict):
            for key, value in dimreduc_data.items():
                # 修改键名以符合scanpy约定
                reduc_name = f"X_{str(key).lower()}"  # 添加 'X_' 前缀
                try:
                    reduc_data = pd.DataFrame(value)
                    obsm[reduc_name] = reduc_data.values
                    LOGGER.info(f"{reduc_name} 维度: {obsm[reduc_name].shape}")
                except Exception as e:
                    LOGGER.warning(f"警告：转换{reduc_name}时出错: {str(e)}")
        else:
            LOGGER.warning("没有找到降维结果")

    # 创建AnnData对象
    LOGGER.info("创建AnnData对象...")
    main_matrix = raw_counts.T.values if use_raw_counts else norm_data.T.values
    LOGGER.info(
        f"使用{'原始计数矩阵' if use_raw_counts else '标准化数据'}作为主要数据 (adata.X)"
    )

    adata = anndata.AnnData(
        X=sparse.csr_matrix(main_matrix, dtype=np.float32),
        obs=metadata_py,
        var=var_df,
        obsm=obsm,
        uns={},
    )

    # 验证AnnData对象中的基因名称
    LOGGER.info("AnnData对象验证:")
    LOGGER.info(f"var_names示例: {list(adata.var_names[:5])}")

    # 设置 raw 属性
    adata.raw = anndata.AnnData(
        X=sparse.csr_matrix(raw_counts.T.values, dtype=np.float32),
        var=var_df.copy(),
        obs=adata.obs,
    )

    # 添加标准化数据到layers
    adata.layers["normalised"] = sparse.csr_matrix(norm_data.T.values, dtype=np.float32)
    LOGGER.info("添加标准化数据到layers中的normalised")

    # 添加counts据到layers
    adata.layers["raw"] = sparse.csr_matrix(raw_counts.T.values, dtype=np.float32)
    LOGGER.info("添加原始计数数据到layers中的raw")

    # 添加高变基因信息
    with localconverter(robjects.default_converter + pandas2ri.converter):
        var_features = r_data[6]
    adata.var["highly_variable"] = adata.var.index.isin(var_features)
    LOGGER.info(f"添加高变基因信息：{sum(adata.var['highly_variable'])}个高变基因")

    # 如果有空间坐标，添加到obsm
    if r_data[7][0]:
        with localconverter(robjects.default_converter + pandas2ri.converter):
            spatial_coords = pd.DataFrame(r_data[8])
        adata.obsm["spatial"] = spatial_coords.values

    return adata


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

    LOGGER.info(f"数据范围检查: {range_valid}")
    LOGGER.info(f"标准化检查: {sum_valid}")

    # If validation fails, attempt normalization
    if not (range_valid and sum_valid):
        LOGGER.info("数据未进行正确的log1p标准化，进行normalize_total和log1p处理...")
        LOGGER.info("检查是否存在raw数据")
        if adata.raw is not None:
            adata = restore_from_raw(adata)
            normalize(adata)
        else:
            raise TypeError("未找到raw数据,请核查数据是否正确")
    else:
        LOGGER.info("数据已经进行了正确的log1p标准化")

    return adata


def restore_from_raw(adata: sc.AnnData) -> sc.AnnData:
    """
    Restore the AnnData object from its raw data.

    Parameters:
        adata (sc.AnnData): The AnnData object with raw data.

    Returns:
        sc.AnnData: A new AnnData object restored from the raw data.
    """
    var_names = adata.var_names
    restored_adata = sc.AnnData(
        X=adata.raw.X, obs=adata.obs, var=adata.raw.var, uns=adata.uns, obsm=adata.obsm
    )
    try:
        restored_adata.var_names = var_names
    except Exception as e:
        LOGGER.info(f"无法恢复变量名: {e}")
    return restored_adata
