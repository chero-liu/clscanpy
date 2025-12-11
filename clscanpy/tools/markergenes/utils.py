import os
import scanpy as sc
from clscanpy.tools.utils import check_mkdir, check_file_exist, save_figure
from anndata import AnnData
import pandas as pd
from clscanpy.log import log_function_call
from clscanpy.tools.plotting.heatmap import plot_marker_heatmap
import clscanpy as oep
from clscanpy.tools.utils import is_valid_path
from pathlib import Path
from clscanpy.tools.color.utils import get_color_order
import logging
import numpy as np
from typing import Optional

LOGGER = logging.getLogger(__name__)

COLUMN_MAP = {
    "names": "gene",
    "logfoldchanges": "log2FoldChange",
    "pvals": "p-value",
    "pvals_adj": "q-value",
    "pct_nz_group": "pct.1",
    "pct_nz_reference": "pct.2",
}


def parse_marker_table(file_path):
    """
    将 marker 表格内容转换为字典。

    Args:
        file_path (str): 文件路径。

    Returns:
        dict: 转换后的字典，键为细胞类型，值为对应的 marker 列表。
    """
    marker_dict = {}

    with open(file_path, "r") as f:
        lines = f.readlines()

        cell_types = lines[0].strip().split("\t")

        for line in lines[1:]:
            markers = line.strip().split("\t")
            for i, marker in enumerate(markers):
                if marker:
                    marker_dict.setdefault(cell_types[i], []).append(marker)
    return marker_dict


def filter_and_adjust_markers(marker_dict, adata):
    """
    过滤 marker 表格中的基因，保留在 adata 中的基因，并调整大小写与 adata 中一致。

    Args:
        marker_dict (dict): marker 表格转换后的字典。
        adata (AnnData): AnnData 对象。

    Returns:
        dict: 过滤并调整后的 marker 字典。
    """
    adata_genes_map = {gene.lower(): gene for gene in adata.var_names}
    filtered_marker_dict = {}
    removed_genes = []

    for cell_type, markers in marker_dict.items():
        adjusted_markers = []
        for marker in markers:
            marker_lower = marker.lower()
            if marker_lower in adata_genes_map:
                adjusted_markers.append(adata_genes_map[marker_lower])
            else:
                removed_genes.append(marker)
        filtered_marker_dict[cell_type] = adjusted_markers

    if removed_genes:
        print("Genes not found in adata:")
        print(", ".join(removed_genes))

    return filtered_marker_dict


def write_diff_genes(
    adata: AnnData,
    outdir: str,
    logfc: float = 0.25,
    minpct: float = 0.25,
    pvals: float = 0.05,
    pvals_adj: float = None,
    refgenome: str = None,
    groupby: str = None,
    top_n: int = 10,
):
    """
    Write differential gene expression results to files with annotation and filtering.

    Parameters:
        adata (AnnData): Annotated data matrix
        outdir (str): Output directory path
        logfc (float): Log-fold change threshold (default: 0.25)
        minpct (float): Minimum percentage expressed (default: 0.25)
        pvals (float): P-value threshold (default: 0.05)
        pvals_adj (float): Adjusted p-value threshold (optional)
        refgenome (str): Path to genome annotation (optional)
        groupby (str): Grouping variable name (optional)
        top_n (int): Number of top markers to save per cluster (default: 10)
    """

    # 1. 配置基本参数
    output_path = Path(outdir)
    output_path.mkdir(parents=True, exist_ok=True)
    use_pvals_adj = pvals_adj is not None
    pval_threshold = pvals_adj if use_pvals_adj else pvals
    pval_label = "qval" if use_pvals_adj else "pval"
    pval_col = "pvals_adj" if use_pvals_adj else "pvals"
    pval_short_label = "q-value" if use_pvals_adj else "p-value"

    # 2. 定义内部辅助函数
    def load_annotation(ref_path: Path) -> pd.DataFrame:
        """加载基因注释文件，支持自动尝试UTF-8和GBK编码"""
        if not ref_path.exists():
            LOGGER.info(f"Annotation file not found at {ref_path}")
            return None

        # 先尝试默认编码(通常是UTF-8)
        try:
            return pd.read_csv(ref_path, sep="\t")
        except UnicodeDecodeError:
            # 捕获编码错误时尝试GBK编码
            LOGGER.warning(
                f"UTF-8 decode failed for {ref_path}, trying GBK encoding..."
            )
            try:
                return pd.read_csv(ref_path, sep="\t", encoding="gbk")
            except Exception as e:
                LOGGER.error(f"GBK decode also failed for {ref_path}: {str(e)}")
                return None
        except Exception as e:
            # 处理其他非编码相关错误
            LOGGER.warning(f"Failed to load annotation from {ref_path}: {str(e)}")
            return None

    def process_comparison() -> tuple:
        """确定比较类型和相关参数"""
        rank_params = adata.uns["rank_genes_groups"]
        groups = rank_params["names"].dtype.names
        reference = rank_params["params"]["reference"]
        is_pairwise = (
            reference != "rest"
            and "pts" in rank_params
            and rank_params["pts"].shape[1] == 2
        )
        return groups, reference, is_pairwise

    def create_renamed_df(group: str) -> pd.DataFrame:
        """为每个组创建重命名的差异表达数据框"""
        df = sc.get.rank_genes_groups_df(adata, group=group)
        df["cluster"] = group
        return df

    def filter_and_count(df: pd.DataFrame) -> tuple:
        """应用过滤条件并返回统计信息"""
        # 通用过滤条件
        filtered = df[
            (df[pval_col] <= pval_threshold) & (df["logfoldchanges"].abs() > logfc)
        ]

        # 添加调控方向
        filtered["Regulation"] = np.where(filtered["logfoldchanges"] > 0, "Up", "Down")

        # 统计信息
        counts = filtered["Regulation"].value_counts()
        return filtered, counts.get("Up", 0), counts.get("Down", 0)

    def save_result(df: pd.DataFrame, filename: str, drop_cols: list = None):
        """保存结果到文件"""
        if drop_cols:
            df = df.drop(columns=drop_cols, errors="ignore")
        file_path = output_path / filename
        file_exists = file_path.exists()
        df.to_csv(
            file_path,
            index=False,
            sep="\t",
            encoding="utf-8",
            mode="a" if file_exists else "w",
            header=not file_exists,
        )

    def filter_protein_coding_genes(
        full_df: pd.DataFrame,
        annotation_df: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """
        过滤非蛋白编码基因，保留标注为'protein_coding'的行

        参数:
        full_df -- 待过滤的完整数据集
        annotation_df -- 触发过滤的注释数据集 (默认为None)

        返回:
        过滤后的DataFrame（优先检查'gene_biotype'列，其次'gene_type'列）
        如果未找到有效列或annotation_df为空，返回原始DataFrame
        """
        # 仅当提供annotation_df时执行过滤
        if annotation_df is not None:
            # 按优先级检查列名
            for col in ["gene_biotype", "gene_type"]:
                if col in full_df.columns:
                    # 创建蛋白编码基因掩码
                    protein_mask = full_df[col] == "protein_coding"
                    # 应用过滤条件
                    return full_df[protein_mask]

        return full_df  # 没有过滤条件时返回原始数据

    # 3. 加载基因注释
    annotation_df = load_annotation(Path(refgenome)) if refgenome else None
    # 4. 处理比较类型
    groups, reference, is_pairwise = process_comparison()

    # 5. 核心处理流程
    if is_pairwise:
        # 成对比较处理
        stats_rows = []
        reference_group = adata.uns["rank_genes_groups"]["pts"].columns[1]

        for group in groups:
            # 处理每个组的差异表达
            de_df = create_renamed_df(group)

            # 添加参考组表达百分比
            pts_df = pd.DataFrame(
                {
                    "names": adata.uns["rank_genes_groups"]["pts"].index,
                    "pct_nz_reference": adata.uns["rank_genes_groups"]["pts"].iloc[
                        :, 1
                    ],
                }
            )
            de_df = de_df.merge(pts_df, on="names", how="left")
            de_df[["FoldChange"]] = 2 ** de_df[["logfoldchanges"]]
            # 应用过滤条件
            filtered_df, up_count, down_count = filter_and_count(de_df)

            # 收集统计信息
            stats_rows.append(
                {
                    "Case": group,
                    "Control": reference_group,
                    "Up_diff": up_count,
                    "Down_diff": down_count,
                    f"Total_diff({pval_short_label}<{pval_threshold}&FoldChange>{2 ** logfc})": f"{up_count + down_count}",
                }
            )

            # 保存结果
            suffix = f"{groupby}_{group}-vs-{reference_group}"
            filtered_df = filtered_df.rename(columns=COLUMN_MAP)
            filtered_df = (
                filtered_df.merge(annotation_df, left_on="gene", right_on="id")
                if annotation_df is not None
                else filtered_df
            )

            if len(filtered_df) != (up_count + down_count):
                LOGGER.warning(
                    f"Filtered count mismatch for {suffix}: expected {up_count + down_count}, got {len(filtered_df)}, check annotation file {refgenome}"
                )

            save_result(
                filtered_df,
                f"{suffix}-diff-{pval_label}-{pval_threshold}-FC-{2 ** logfc}_anno.xls",
                drop_cols=["cluster", "scores"],
            )
            de_df = de_df.rename(columns=COLUMN_MAP)
            # 保存完整结果
            de_df = (
                de_df.merge(annotation_df, left_on="gene", right_on="id")
                if annotation_df is not None
                else de_df
            )
            save_result(
                de_df,
                f"{suffix}-all_diffexp_genes_anno.xls",
                drop_cols=["cluster", "scores"],
            )

        # 保存统计信息
        if stats_rows:
            save_result(pd.DataFrame(stats_rows), "diffexp_results_stat.xls")

    else:
        # cluster_vs_rest 比较处理
        filtered_dfs = []

        for group in groups:
            de_df = create_renamed_df(group)

            # 应用cluster_vs_rest特定过滤
            filtered = de_df[
                (de_df[pval_col] <= pval_threshold)
                & (
                    (de_df["pct_nz_group"] >= minpct)
                    | (de_df["pct_nz_reference"] >= minpct)
                )
                & (de_df["logfoldchanges"] >= 0)
            ]
            filtered_dfs.append(filtered)

        if filtered_dfs:
            # 合并所有结果
            full_df = pd.concat([create_renamed_df(g) for g in groups])
            filtered_df = pd.concat(filtered_dfs)
            filtered_df["gene_diff"] = np.round(
                filtered_df["pct_nz_group"] / filtered_df["pct_nz_reference"], 3
            )

            filtered_df = (
                filtered_df.merge(annotation_df, left_on="names", right_on="id")
                if annotation_df is not None
                else filtered_df
            )

            if annotation_df is not None:
                filtered_df = filter_protein_coding_genes(filtered_df, annotation_df)

            cluster_order = filtered_df["cluster"].unique().tolist()

            temp_df = filtered_df.assign(
                cluster_cat=pd.Categorical(
                    filtered_df["cluster"], categories=cluster_order, ordered=True
                )
            )

            top_markers_df = (
                temp_df.sort_values(
                    ["cluster_cat", "gene_diff"], ascending=[True, False]
                )
                .groupby("cluster", sort=False)
                .head(top_n)
                .drop(columns="cluster_cat")
                .reset_index(drop=True)
                .rename(columns=COLUMN_MAP)
            )
            top_markers_df = top_markers_df.rename(columns=COLUMN_MAP)
            top_markers_df = top_markers_df.rename(
                columns={"log2FoldChange": "avg_log2FC"}
            )

            save_result(top_markers_df, f"top{top_n}_markers_for_each_cluster_anno.xls")

            # 保存所有标记基因
            filtered_df = filtered_df.rename(columns=COLUMN_MAP)
            filtered_df = filtered_df.rename(columns={"log2FoldChange": "avg_log2FC"})

            save_result(filtered_df, "all_markers_for_each_cluster_anno.xls")


def marker_topn(
    adata,
    top: int = 3,
):
    if "rank_genes_groups" not in adata.uns:
        raise ValueError(
            "'rank_genes_groups' not in adata.uns, please rank genes first"
        )

    groups = adata.uns["rank_genes_groups"]["names"].dtype.names
    topn = {}

    for i in groups:
        de_i = sc.get.rank_genes_groups_df(adata, group=i)
        topn[i] = de_i["names"][:top].tolist()
    return topn


def vis_markers(
    adata,
    top,
    groupby,
    outdir,
    method="umap",
    palette=None,
):
    if is_valid_path(top):
        marker_dict = parse_marker_table(top)
        marker_dict = filter_and_adjust_markers(marker_dict, adata)
        top = ""
    else:
        top = int(top)
        marker_dict = marker_topn(adata, top=top)
        top = f"top{top}"

        plot_marker_heatmap(
            adata,
            markers=marker_dict,
            groupby=groupby,
            save_path=f"{outdir}/{top}heatmap",
            show=False,
            palette=palette,
        )

    sc.pl.stacked_violin(
        adata,
        marker_dict,
        groupby=groupby,
        dendrogram=False,
        use_raw=False,
        show=False,
    )
    save_figure(f"{outdir}/{top}stacked_violin")

    sc.pl.dotplot(
        adata,
        marker_dict,
        groupby=groupby,
        dendrogram=False,
        use_raw=False,
        swap_axes=False,
        show=False,
    )
    save_figure(f"{outdir}/{top}dotplot")

    for key in marker_dict:
        oep.plot_dim(
            adata, method=method, color=marker_dict[key], show=False, figsize=(5, 5)
        )
        save_figure(f"{outdir}/markers_vis4cluster{key}_featureplot")
