#!/usr/bin/env python3
"""
可视化脚本 - 生成富集分析结果的可视化图表
"""

import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import gseapy as gp
from gseapy.scipalette import SciPalette
import warnings
import os
import sys
import json
from typing import Dict, List, Optional

matplotlib.use("Agg")


def parse_visualization_arguments():
    """解析可视化脚本的命令行参数"""
    parser = argparse.ArgumentParser(
        description="生成富集分析结果的可视化图表",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # 必需参数
    parser.add_argument("--input-dir", required=True, help="分析结果目录路径")
    parser.add_argument("--output-dir", required=True, help="可视化输出目录路径")
    parser.add_argument(
        "--summary-file", required=True, help="分析摘要文件路径 (analysis_summary.json)"
    )

    # 可选参数
    parser.add_argument(
        "--top-n", type=int, default=20, help="每个类别显示的前N个富集结果 (默认: 20)"
    )
    parser.add_argument(
        "--pvalue", type=float, default=0.05, help="p值阈值 (默认: 0.05)"
    )
    parser.add_argument("--verbose", action="store_true", help="显示详细输出信息")

    return parser.parse_args()


def setup_visualization_dirs(output_dir: str) -> Dict[str, str]:
    """创建可视化输出目录结构"""
    dirs = {
        "root": output_dir,
        "figures": os.path.join(output_dir, "figures"),
    }

    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)

    return dirs


def load_analysis_summary(summary_file: str) -> Dict:
    """加载分析摘要信息"""
    try:
        with open(summary_file, "r") as f:
            summary = json.load(f)
        return summary
    except Exception as e:
        print(f"错误: 无法加载分析摘要文件 {summary_file}: {e}")
        sys.exit(1)


def plot_enrichment_results(
    enr_df: pd.DataFrame,
    output_dir: str,
    prefix: str,
    top_n: int = 20,
    pvalue: float = 0.05,
):
    """绘制富集分析结果图"""

    if enr_df.empty:
        return

    # 过滤和排序
    enr_df = enr_df[enr_df["Adjusted P-value"] <= pvalue]
    enr_df = enr_df.sort_values("Adjusted P-value").head(top_n)

    if enr_df.empty:
        return

    sci = SciPalette()
    NbDr = sci.create_colormap()

    # 1. GO条形图 (按类别分组)
    if "category" in enr_df.columns:
        try:
            ax = gp.barplot(
                enr_df,
                figsize=(10, max(6, len(enr_df) * 0.5)),
                group="category",
                color=NbDr.reversed(),
                cutoff=pvalue,
                top_term=top_n,
            )
            ax.tick_params(axis="x", rotation=45)
            plt.tight_layout()
            ax.figure.savefig(
                os.path.join(output_dir, f"{prefix}_GO_barplot.png"),
                dpi=300,
                bbox_inches="tight",
            )
            ax.figure.savefig(
                os.path.join(output_dir, f"{prefix}_GO_barplot.pdf"),
                dpi=300,
                bbox_inches="tight",
            )
            plt.close(ax.figure)
        except Exception as e:
            warnings.warn(f"绘制条形图失败: {str(e)}")

    # 2. GO点图 (按类别分组)
    if "category" in enr_df.columns:
        try:
            ax = gp.dotplot(
                enr_df,
                figsize=(10, max(6, len(enr_df) * 0.5)),
                x="category",
                x_order=["BP", "CC", "MF"],
                cmap=NbDr.reversed(),
                size=5,
                cutoff=pvalue,
                top_term=top_n,
                show_ring=True,
            )
            plt.tight_layout()
            ax.figure.savefig(
                os.path.join(output_dir, f"{prefix}_GO_dotplot.png"),
                dpi=300,
                bbox_inches="tight",
            )
            ax.figure.savefig(
                os.path.join(output_dir, f"{prefix}_GO_dotplot.pdf"),
                dpi=300,
                bbox_inches="tight",
            )
            plt.close(ax.figure)
        except Exception as e:
            warnings.warn(f"绘制点图失败: {str(e)}")

    # 3. 单个点图（不按类别分组）
    try:
        ax = gp.dotplot(
            enr_df,
            figsize=(10, max(6, len(enr_df) * 0.5)),
            cmap=plt.cm.viridis_r,
            cutoff=pvalue,
            top_term=top_n,
        )
        plt.tight_layout()
        ax.figure.savefig(
            os.path.join(output_dir, f"{prefix}_enrichment_dotplot.png"),
            dpi=300,
            bbox_inches="tight",
        )
        ax.figure.savefig(
            os.path.join(output_dir, f"{prefix}_enrichment_dotplot.pdf"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(ax.figure)
    except Exception as e:
        warnings.warn(f"绘制单个点图失败: {str(e)}")


def plot_kegg_results(
    kegg_df: pd.DataFrame,
    output_dir: str,
    prefix: str,
    top_n: int = 20,
    pvalue: float = 0.05,
):
    """绘制KEGG富集分析结果图"""
    if kegg_df.empty:
        return

    # 过滤和排序
    kegg_df = kegg_df[kegg_df["Adjusted P-value"] <= pvalue]
    kegg_df = kegg_df.sort_values("Adjusted P-value").head(top_n)

    if kegg_df.empty:
        return

    # KEGG条形图
    try:
        ax = gp.barplot(
            kegg_df,
            figsize=(10, max(6, len(kegg_df) * 0.3)),
            title=f"{prefix} - KEGG Enrichment",
            color="b",
            cutoff=pvalue,
            top_term=top_n,
        )
        plt.tight_layout()
        ax.figure.savefig(
            os.path.join(output_dir, f"{prefix}_KEGG_barplot.png"),
            dpi=300,
            bbox_inches="tight",
        )
        ax.figure.savefig(
            os.path.join(output_dir, f"{prefix}_KEGG_barplot.pdf"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(ax.figure)
    except Exception as e:
        warnings.warn(f"绘制KEGG条形图失败: {str(e)}")

    # KEGG点图
    try:
        ax = gp.dotplot(
            kegg_df,
            figsize=(10, max(6, len(kegg_df) * 0.3)),
            title=f"{prefix} - KEGG Enrichment",
            cmap=plt.cm.viridis_r,
            cutoff=pvalue,
            top_term=top_n,
        )
        plt.tight_layout()
        ax.figure.savefig(
            os.path.join(output_dir, f"{prefix}_KEGG_dotplot.png"),
            dpi=300,
            bbox_inches="tight",
        )
        ax.figure.savefig(
            os.path.join(output_dir, f"{prefix}_KEGG_dotplot.pdf"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(ax.figure)
    except Exception as e:
        warnings.warn(f"绘制KEGG点图失败: {str(e)}")


def generate_visualizations(args, summary_data: Dict, input_dir: str):
    """生成所有可视化图表"""
    from datetime import datetime

    print("=" * 60)
    print("富集分析脚本 - 可视化模块")
    print("=" * 60)
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"输入目录: {input_dir}")

    output_dirs = setup_visualization_dirs(args.output_dir)

    if args.verbose:
        print(f"创建可视化输出目录: {args.output_dir}")

    tables_dir = os.path.join(input_dir, "tables")
    if not os.path.exists(tables_dir):
        print(f"错误: 分析结果目录不存在: {tables_dir}")
        sys.exit(1)

    print(f"\n开始生成可视化图表...")

    for cluster_id in summary_data.keys():
        print(f"\n为Cluster {cluster_id}生成可视化图表...")

        # 检查是否有GO合并结果
        go_combined_file = os.path.join(
            tables_dir, f"cluster_{cluster_id}_GO_combined.tsv"
        )
        if os.path.exists(go_combined_file):
            try:
                go_df = pd.read_csv(go_combined_file, sep="\t")
                if not go_df.empty:
                    plot_enrichment_results(
                        go_df,
                        output_dirs["figures"],
                        f"cluster_{cluster_id}",
                        args.top_n,
                        args.pvalue,
                    )
                    print(f"  ✓ 生成GO可视化图表")
                else:
                    print(f"  ⚠ GO合并结果为空")
            except Exception as e:
                print(f"  ✗ 处理GO合并结果失败: {e}")
        else:
            print(f"  ⚠ 未找到GO合并结果文件: {go_combined_file}")

        # 检查KEGG结果
        kegg_file = os.path.join(
            tables_dir, f"cluster_{cluster_id}_KEGG_enrichment.tsv"
        )
        if os.path.exists(kegg_file):
            try:
                kegg_df = pd.read_csv(kegg_file, sep="\t")
                if not kegg_df.empty:
                    plot_kegg_results(
                        kegg_df,
                        output_dirs["figures"],
                        f"cluster_{cluster_id}",
                        args.top_n,
                        args.pvalue,
                    )
                    print(f"  ✓ 生成KEGG可视化图表")
                else:
                    print(f"  ⚠ KEGG结果为空")
            except Exception as e:
                print(f"  ✗ 处理KEGG结果失败: {e}")
        else:
            print(f"  ⚠ 未找到KEGG结果文件: {kegg_file}")

    print(f"\n可视化完成！图表保存在: {args.output_dir}")


def run_visualization(args):
    """运行可视化主函数"""
    # 加载分析摘要
    summary_data = load_analysis_summary(args.summary_file)

    if args.verbose:
        print(f"加载分析摘要: {len(summary_data)} 个clusters")

    # 生成可视化图表
    generate_visualizations(args, summary_data, args.input_dir)


if __name__ == "__main__":
    args = parse_visualization_arguments()
    run_visualization(args)
