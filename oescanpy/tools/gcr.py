# get clusters result
import os
import pandas as pd
import oescanpy as oep
import unittest
from oescanpy.tools.utils import (
    check_mkdir,
    add_counts_to_data,
    save_figure,
    get_group_order,
)
from oescanpy.tools.color.utils import get_color_order
from oescanpy.log import log_function_call

import logging
LOGGER = logging.getLogger(__name__)
import json
import numpy as np
import scipy


def build_gene_index(input_file, gene_ids, output_file):
    """
    构建基因在虚拟连续文本中的字节位置索引
    处理标题行问题并正确计算字节位置

    参数:
    input_file: 输入表达式TSV文件路径
    gene_ids: 基因ID列表（长度需与基因行数一致）
    output_file: 输出索引文件路径
    """
    # 检查文件是否存在
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"输入文件不存在: {input_file}")

    # 使用二进制模式精确计算字节位置
    with open(input_file, "rb") as f:
        # 读取整个文件内容
        content = f.read()

    # 分割为行（保留换行符）
    lines = content.splitlines(keepends=True)

    # 验证行数
    if len(lines) < 2:
        raise ValueError("输入文件至少需要一行标题和一行数据")

    # 只处理基因行（跳过标题行）
    gene_lines = lines[1:]

    # 验证基因行数与提供的gene_ids数量匹配
    if len(gene_lines) != len(gene_ids):
        raise ValueError(
            f"基因ID数量({len(gene_ids)})与文件基因行数({len(gene_lines)})不匹配!"
            f"输入文件: {input_file} 有 {len(gene_lines)} 行基因数据"
        )

    # 计算标题行字节长度（包括换行符）
    header_bytes = len(lines[0])

    # 初始化位置列表
    start_n_list = []
    end_n_list = []

    # 当前字节位置（标题行之后开始）
    current_pos = header_bytes

    # 计算每个基因行的位置
    for line in gene_lines:
        # 行字节长度（包括换行符）
        line_bytes = len(line)

        # 起始位置 = 当前字节位置 + 1（1-based索引）
        start_pos = current_pos + 1
        # 结束位置 = 当前字节位置 + 行字节长度
        end_pos = current_pos + line_bytes

        # 添加到列表
        start_n_list.append(start_pos)
        end_n_list.append(end_pos)

        # 更新当前字节位置
        current_pos = end_pos

    # 创建结果数据框
    result_df = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "start": [str(x) for x in start_n_list],
            "end": [str(x) for x in end_n_list],
        }
    )

    # 写入TSV文件
    result_df.to_csv(
        output_file,
        sep="\t",
        index=False,
        escapechar="\\",
    )
    print(f"成功生成索引文件: {output_file}，包含 {len(result_df)} 个基因索引")


def save_dimension_reduction_results(
    adata,
    outdir,
    method="umap",
    color_by="clusters",
    palette=None,
    clust_method="leiden",
    resolution=None,
):
    """
    Save dimension reduction results (e.g., UMAP, t-SNE) and their visualizations.

    Parameters:
        adata: AnnData object containing the data.
        outdir: Output directory to save results.
        method: Dimension reduction method ('umap' or 'tsne').
        color_by: Column in `adata.obs` to color the plots by.
        clust_method: Key in `adata.uns` to extract resolution for file naming.
    """
    # Validate method
    if method not in ['umap', 'tsne']:
        raise ValueError("Method must be 'umap' or 'tsne'")

    # Create output directory
    outdir = os.path.join(outdir, f'{method}_Dimension_Reduction')
    check_mkdir(outdir)

    # Extract dimension reduction coordinates
    coord_key = f'X_{method}'
    if coord_key not in adata.obsm:
        raise KeyError(f"{coord_key} not found in adata.obsm")

    tmp = pd.DataFrame(adata.obsm[coord_key]).rename(columns={0: f'{method.upper()}_1', 1: f'{method.upper()}_2'})
    tmp.insert(0, 'Barcode', adata.obs['rawbc'].values)
    coord_file = os.path.join(outdir, f'{method}_Dimension_Reduction_coordination.csv')
    tmp.to_csv(coord_file, index=False)

    # add counts to the 'clusters' column
    # adata.obs = add_counts_to_data(adata.obs, color_by)

    oep.plot_dim(
        adata,
        method=method,
        color=color_by,
        palette=get_color_order(adata, color_by, palette),
        show=False,
        figsize=(5, 5),
        title=color_by,
    )

    if resolution == None:
        try:
            resolution = json.loads(adata.uns["preprocess_para"])["resolution"]
        except KeyError:
            resolution = "None"

    plot_file = os.path.join(
        outdir, f"{method}_groupby_{color_by}_resolution{resolution}"
    )
    save_figure(plot_file)

    oep.plot_dim(adata, 
                method=method, 
                color=color_by, 
                palette=get_color_order(adata,color_by,palette),
                show=False, 
                figsize=(5, 5), 
                legend_loc='on data')
    plot_file = os.path.join(
        outdir, f"{method}_groupby_{color_by}_resolution{resolution}_on_data"
    )
    save_figure(plot_file)


def save_clusters_results(
    adata,
    outdir,
    method="umap",
    color_by="clusters",
    resolution=None,
    clust_method="leiden",
    cloudcfg=False,
):
    """
    Save clusters results to a CSV file.

    Parameters:
        adata: AnnData object containing the data.
        outdir: Output directory to save results.
    """
    check_mkdir(outdir)
    clusters_file = os.path.join(outdir, f"{color_by}_result.csv")
    data = adata.obs[["rawbc", "sampleid", "group", "batchid", color_by]]
    data = data.rename(
        columns={
            "rawbc": "Barcode",
        }
    )
    data.to_csv(clusters_file, index=False)

    data = get_group_order(
        data,
        key_col="sampleid",
        value_col="batchid",
    )
    data.to_csv(os.path.join(outdir, "sampleid-batchid.xls"), sep="\t", index=False)
    LOGGER.info(f"Saved clusters results to {clusters_file}")

    if resolution == None:
        try:
            resolution = adata.uns[clust_method]["params"]["resolution"]
        except KeyError:
            resolution = "None"

    if cloudcfg:
        LOGGER.info(f"cloudcfg 为True--开始生成云平台所需文件")

        if "barcode" in adata.obs.columns:
            adata.obs = adata.obs.drop(columns=["barcode"])

        metadata_df = adata.obs.reset_index().rename(columns={"index": "barcode"})

        metadata_df.to_csv(
            os.path.join(outdir, "metadata.tsv"),
            sep="\t",
            na_rep="",
            index=False,
            quoting=3,
        )

        # Extract dimension reduction coordinates
        coord_key = f"X_{method}"
        if coord_key not in adata.obsm:
            raise KeyError(f"{coord_key} not found in adata.obsm")

        tmp = pd.DataFrame(adata.obsm[coord_key]).rename(
            columns={0: f"{method.upper()}_1", 1: f"{method.upper()}_2"}
        )
        tmp.insert(0, "Barcode", adata.obs["rawbc"].values)
        coord_file = os.path.join(
            outdir,
            f"{method}_Dimension_Reduction_coordination.tsv",
        )
        tmp.to_csv(
            coord_file,
            index=False,
            sep="\t",
        )

        embedding_df = (
            pd.DataFrame(
                adata.obsm[f"X_{method}"][:, :2],
                columns=[f"{method.upper()}_1", f"{method.upper()}_2"],
                index=adata.obs.index,
            )
            .reset_index()
            .rename(columns={"index": "cellid"})
        )

        metadata_df = adata.obs.reset_index().rename(columns={"index": "cellid"})
        merged_df = embedding_df.merge(
            metadata_df[["cellid", color_by, "sampleid", "group", f"{color_by}_col"]],
            on="cellid",
            how="left",
        )

        output_file = os.path.join(
            outdir,
            f"{method}_groupby_{color_by}_resolution{resolution}_plot.tsv",
        )
        merged_df.to_csv(output_file, sep="\t", na_rep="", index=False, quoting=3)

        adata.X = adata.layers["normalised"]
        expression_df = (
            pd.DataFrame(
                adata.X.T.toarray() if hasattr(adata.X, "toarray") else adata.X.T,
                index=adata.var.index,
                columns=adata.obs.index,
            )
            .reset_index()
            .rename(columns={"index": "gene_id"})
        )

        expression_file = os.path.join(
            outdir,
            f"{method}_groupby_{color_by}_resolution{resolution}_plot_expression.tsv",
        )
        expression_df.to_csv(expression_file, sep="\t", index=False, quoting=3)

        index_file = os.path.join(
            outdir,
            f"{method}_groupby_{color_by}_resolution{resolution}_plot_gene_index.tsv",
        )
        build_gene_index(expression_file, adata.var.index.tolist(), index_file)


def get_clusters_result(
    adata,
    outdir,
    color_by="clusters",
    method="umap",
    cloudcfg=False,
    palette=None,
    resolution=None,
    clust_method="leiden",
):
    """
    Main function to save clusters and dimension reduction results.

    Parameters:
        adata: AnnData object containing the data.
        outdir: Base output directory.
    """
    # Save clusters results
    save_clusters_results(
        adata,
        outdir=outdir,
        method=method,
        color_by=color_by,
        resolution=resolution,
        cloudcfg=cloudcfg,
        clust_method=clust_method,
    )

    # Save UMAP results
    save_dimension_reduction_results(
        adata,
        outdir,
        method=method,
        color_by=color_by,
        palette=palette,
        clust_method=clust_method,
        resolution=resolution,
    )


class GCR:

    def __init__(
        self,
        input,
        outdir,
        groupby,
        groupby_levels,
        sampleid,
        group,
        clusters,
        new_celltype,
        predicate,
        metadata,
        clust_method,
        method,
        palette,
        cloudcfg,
        resolution,
        save_h5ad: bool = False,
    ):
        self.input = input
        self.outdir = outdir
        self.groupby = groupby
        self.groupby_levels = groupby_levels
        self.sampleid = sampleid
        self.group = group
        self.clusters = clusters
        self.new_celltype = new_celltype
        self.metadata = metadata
        self.predicate = predicate
        self.clust_method = clust_method
        self.method = method
        self.palette = palette
        self.cloudcfg = cloudcfg
        self.resolution = resolution
        self.save_h5ad = save_h5ad

    def run(self):
        LOGGER.info("Start GCR ...")
        adata = oep.loadH5AD(
            self.input,
            sampleid=self.sampleid,
            group=self.group,
            clusters=self.clusters,
            new_celltype=self.new_celltype,
            predicate=self.predicate,
            metadata=self.metadata,
            groupby=self.groupby,
            groupby_levels=self.groupby_levels,
            palette=self.palette,
        )

        get_clusters_result(
            adata,
            self.outdir,
            color_by=self.groupby,
            cloudcfg=self.cloudcfg,
            method=self.method,
            palette=self.palette,
            clust_method=self.clust_method,
            resolution=self.resolution,
        )

        if self.save_h5ad:
            adata.write_h5ad(os.path.join(self.outdir, "adata.h5ad"))

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        LOGGER.info("gcr done")

@log_function_call
def gcr(args):
    with GCR(
        input=args.input,
        outdir=args.outdir,
        groupby=args.groupby,
        groupby_levels=args.groupby_levels,
        sampleid=args.sampleid,
        group=args.group,
        clusters=args.clusters,
        new_celltype=args.new_celltype,
        predicate=args.predicate,
        metadata=args.metadata,
        clust_method=args.clust_method,
        method=args.method,
        palette=args.palette,
        cloudcfg=args.cloudcfg,
        resolution=args.resolution,
        save_h5ad=args.save_h5ad,
    ) as runner:
        runner.run()

def get_opts_gcr(parser, sub_program=True):
    parser.add_argument('-i', "--input", type=str, default=None, help="Input file", required=True)
    parser.add_argument('-o', "--outdir", type=str, default=None, help="Output directory", required=True)
    parser.add_argument('-c', "--groupby", type=str, default='clusters', help="Color by column")
    parser.add_argument('-l', "--groupby_levels", type=str, default=None, help="Groupby levels")
    parser.add_argument(
        "--clust_method", type=str, default="leiden", help="Resolution key"
    )
    parser.add_argument("--resolution", type=str, default=None, help="resolution")
    parser.add_argument('-m', "--method", type=str, default='umap', help="Dimension reduction method (umap or tsne)")
    parser.add_argument("-p", "--palette", type=str, default=None, help="")
    parser.add_argument("--cloudcfg", type=bool, default=False, help="")
    parser.add_argument(
        "--sampleid",
        type=str,
        default=None,
        help="sampleid column name in adata.obs",
    )
    parser.add_argument(
        "--group",
        type=str,
        default=None,
        help="group column name in adata.obs",
    )
    parser.add_argument(
        "--clusters",
        type=str,
        default=None,
        help="clusters column name in adata.obs",
    )
    parser.add_argument(
        "--new_celltype",
        type=str,
        default=None,
        help="new celltype column name in adata.obs",
    )
    parser.add_argument(
        "--predicate",
        type=str,
        default=None,
        help="predicate for filtering adata.obs, e.g. (sampleid in ['STB1', 'STB4']) and ~(clusters in ['1', '5', '6'])",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        help="Metadata file to load, if not provided, will use adata.obs",
    )
    parser.add_argument("--save_h5ad", default=True, help="")
    return parser


if __name__ == "__main__":
    unittest.main()
