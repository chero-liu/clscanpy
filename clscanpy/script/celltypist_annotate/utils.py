import scanpy as sc
import pandas as pd
import numpy as np
import logging
import sys
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import anndata
import logging

LOGGER = logging.getLogger(__name__)

def update_h5seurat_metadata(adata, h5seurat_path, output_h5seurat=None, col_names=None):
    """
    将AnnData的obs数据更新到h5seurat文件中
    
    Parameters:
    -----------
    adata : anndata.AnnData
        输入的AnnData对象，包含要更新的metadata
    h5seurat_path : str
        输入的h5seurat文件路径
    output_h5seurat : str, optional
        输出的h5seurat文件路径。如果为None，则直接更新输入文件
    col_names : list, optional
        要更新的列名列表。如果为None，则更新所有列
        
    Returns:
    --------
    None
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    # 如果没有指定输出路径，则使用输入路径
    if output_h5seurat is None:
        LOGGER.info("没有指定输出路径，请核查")

    # R代码
    r_code = """
    update_seurat_metadata <- function(seurat_path, metadata, output_path) {
        tryCatch({
            # 加载必要的包
            suppressPackageStartupMessages({
                library(Seurat)
                library(SeuratObject)
                library(SeuratDisk)
            })
            
            # 读取Seurat对象
            print("Reading Seurat object...")
            if (endsWith(seurat_path, ".h5seurat")) {
                seurat_obj <- LoadH5Seurat(seurat_path)
            } else {
                seurat_obj <- readRDS(seurat_path)
            }
            print("Seurat object loaded successfully")
            
            # 检查细胞名称是否匹配
            if (!all(rownames(metadata) %in% colnames(seurat_obj))) {
                stop("Cell names in metadata do not match with Seurat object")
            }
            
            # 更新metadata
            message("Updating metadata...")
            for (col in colnames(metadata)) {
                if (col %in% colnames(seurat_obj@meta.data)) {
                    message(paste("Updating existing column:", col))
                } else {
                    message(paste("Adding new column:", col))
                }
                seurat_obj[[col]] <- metadata[colnames(seurat_obj), col]
            }
            
            # 保存更新后的对象
            message("Saving updated h5seurat file...")
            SaveH5Seurat(seurat_obj, filename = output_path, overwrite = TRUE)
            
            message("Update completed successfully!")
        }, error = function(e) {
            message(paste("Error occurred:", e$message))
            stop(e$message)
        })
    }
    """

    LOGGER.info("开始更新metadata...")

    # 准备metadata
    metadata_df = adata.obs.copy()
    if col_names is not None:
        # 验证所有请求的列是否存在
        missing_cols = [col for col in col_names if col not in metadata_df.columns]
        if missing_cols:
            raise ValueError(f"以下列在adata.obs中不存在: {missing_cols}")
        # 只选择指定的列
        metadata_df = metadata_df[col_names]
        LOGGER.info(f"选择更新以下列: {col_names}")

    # 转换为R对象
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_metadata = ro.conversion.py2rpy(metadata_df)

    # 运行R代码
    ro.r(r_code)
    update_func = ro.globalenv['update_seurat_metadata']

    try:
        LOGGER.info(f"正在更新h5seurat文件: {h5seurat_path}")
        LOGGER.info(f"输出路径: {output_h5seurat}")
        update_func(h5seurat_path, r_metadata, output_h5seurat)
        LOGGER.info("更新完成!")
    except Exception as e:
        logger.error(f"更新过程中出错: {str(e)}")
        raise
