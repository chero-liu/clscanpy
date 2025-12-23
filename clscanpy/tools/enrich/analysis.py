#!/usr/bin/env python3
"""
富集分析脚本 - 对差异表达基因进行GO和KEGG富集分析
"""

import argparse
import pandas as pd
import gseapy as gp
import warnings
import numpy as np
import os
import sys
from typing import Dict, List, Tuple, Optional, Any


def parse_analysis_arguments():
    """解析分析脚本的命令行参数"""
    parser = argparse.ArgumentParser(
        description="对差异表达基因进行GO和KEGG富集分析",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # 必需参数
    parser.add_argument("--input", required=True, help="输入文件路径 (TSV格式)")
    parser.add_argument("--outdir", required=True, help="输出目录路径")
    parser.add_argument("--bp", required=True, help="GO BP基因集文件路径")
    parser.add_argument("--cc", required=True, help="GO CC基因集文件路径")
    parser.add_argument("--mf", required=True, help="GO MF基因集文件路径")
    parser.add_argument("--kegg", required=True, help="KEGG基因集文件路径")

    # 可选参数
    parser.add_argument(
        "--clusters",
        default="",
        help="指定分析的cluster列表，逗号分隔，如 '0,1,2' (默认: 分析所有cluster)",
    )
    parser.add_argument(
        "--pvalue", type=float, default=0.05, help="p值阈值 (默认: 0.05)"
    )
    parser.add_argument(
        "--gene-column", default="names", help="基因名列名 (默认: 'names')"
    )
    parser.add_argument(
        "--cluster-column", default="clusters", help="cluster列名 (默认: 'clusters')"
    )
    parser.add_argument(
        "--min-genes", type=int, default=3, help="富集分析最小基因数 (默认: 3)"
    )
    parser.add_argument(
        "--max-genes", type=int, default=500, help="富集分析最大基因数 (默认: 500)"
    )
    parser.add_argument(
        "--truncate-length", type=int, default=30, help="术语截断长度 (默认: 30)"
    )
    parser.add_argument("--verbose", action="store_true", help="显示详细输出信息")

    return parser.parse_args()


def setup_output_dirs(outdir: str) -> Dict[str, str]:
    """创建输出目录结构"""
    tables_dir = os.path.join(outdir, "tables")

    dirs = {
        "root": outdir,
        "tables": tables_dir,
    }

    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)

    return dirs


def load_gene_sets(args) -> Dict[str, str]:
    """加载基因集文件并验证"""
    gene_sets = {"BP": args.bp, "CC": args.cc, "MF": args.mf, "KEGG": args.kegg}

    # 验证文件是否存在
    for name, path in gene_sets.items():
        if not os.path.exists(path):
            raise FileNotFoundError(f"基因集文件不存在: {name} -> {path}")
        if args.verbose:
            print(f"✓ 加载基因集: {name} -> {path}")

    return gene_sets


def smart_truncate(
    df, column_name, max_length=20, truncate_indicator="...", inplace=False
):
    """
    智能截断字符串，保留重要信息
    """
    if not inplace:
        df = df.copy()

    if column_name not in df.columns:
        raise ValueError(f"列 '{column_name}' 不存在于数据框中")

    available_length = max_length - len(truncate_indicator)

    def truncate_func(val):
        if pd.isna(val):
            return val
        str_val = str(val)
        if len(str_val) > max_length:
            if " " in str_val[:available_length]:
                last_space = str_val[:available_length].rfind(" ")
                if last_space > max_length * 0.5:
                    return str_val[:last_space] + truncate_indicator

            if "(" in str_val[:available_length]:
                last_paren = str_val[:available_length].rfind("(")
                if last_paren > max_length * 0.5:
                    return str_val[:last_paren] + truncate_indicator

            for sep in [",", ";", ":"]:
                if sep in str_val[:available_length]:
                    last_sep = str_val[:available_length].rfind(sep)
                    if last_sep > max_length * 0.5:
                        return str_val[:last_sep] + truncate_indicator

            return str_val[:available_length] + truncate_indicator
        else:
            return str_val

    df[f"raw_{column_name}"] = df[column_name]
    df[column_name] = df[column_name].apply(truncate_func)

    return df


def run_enrichr(
    gene_list: List[str], gmt_path: str, min_genes: int = 3, max_genes: int = 500
) -> Optional[Any]:
    """
    运行Enrichr富集分析
    """
    try:
        # 过滤基因列表
        filtered_genes = gene_list[:max_genes]
        if len(filtered_genes) < min_genes:
            warnings.warn(f"基因数量({len(filtered_genes)})少于最小要求({min_genes})")
            return None

        enr = gp.enrichr(
            gene_list=filtered_genes,
            gene_sets=[gmt_path],
            outdir=None,
            verbose=False,
            cutoff=1.0,  # 不过滤，后面自己过滤
        )

        if enr.res2d.empty:
            return None

        return enr

    except Exception as e:
        warnings.warn(f"富集分析失败: {gmt_path}, 错误: {str(e)}")
        return None


def save_enrichment_results(enr_result, output_path: str, category: str):
    """保存富集分析结果到文件"""
    if enr_result is not None and not enr_result.res2d.empty:
        enr_result.res2d.to_csv(output_path, sep="\t", index=False)
        return True
    return False


def analyze_cluster(
    cluster_id, cluster_genes, gene_sets, output_dirs, args, summary_info
):
    """分析单个cluster"""
    results = {}

    if args.verbose:
        print(f"  Cluster {cluster_id}: 分析 {len(cluster_genes)} 个基因")

    # 运行富集分析
    for category, gmt_path in gene_sets.items():
        if args.verbose:
            print(f"    分析 {category}...")

        enr = run_enrichr(cluster_genes, gmt_path, args.min_genes, args.max_genes)

        if enr is not None and not enr.res2d.empty:
            # 添加类别标签
            enr.res2d["category"] = category
            # 智能截断Term名称
            enr.res2d = smart_truncate(enr.res2d, "Term", args.truncate_length)

            # 保存结果
            output_file = os.path.join(
                output_dirs["tables"], f"cluster_{cluster_id}_{category}_enrichment.tsv"
            )
            if save_enrichment_results(enr, output_file, category):
                results[category] = {
                    "file": output_file,
                    "count": len(enr.res2d),
                    "significant": len(
                        enr.res2d[enr.res2d["Adjusted P-value"] <= args.pvalue]
                    ),
                }
            else:
                results[category] = None
        else:
            results[category] = None

    # 合并GO结果
    go_results = []
    for category in ["BP", "CC", "MF"]:
        if results.get(category):
            result_file = results[category]["file"]
            df = pd.read_csv(result_file, sep="\t")
            go_results.append(df)

    if go_results:
        go_df = pd.concat(go_results, ignore_index=True)
        go_df = go_df[go_df["Adjusted P-value"] <= args.pvalue]

        if not go_df.empty:
            # 保存合并的GO结果
            go_output = os.path.join(
                output_dirs["tables"], f"cluster_{cluster_id}_GO_combined.tsv"
            )
            go_df.to_csv(go_output, sep="\t", index=False)
            results["GO_combined"] = {"file": go_output, "count": len(go_df)}

    # 收集分析统计信息
    cluster_summary = {"gene_count": len(cluster_genes), "significant_counts": {}}

    for category in ["BP", "CC", "MF", "KEGG"]:
        if results.get(category):
            info = results[category]
            cluster_summary["significant_counts"][category] = info.get("significant", 0)
        else:
            cluster_summary["significant_counts"][category] = 0

    summary_info[cluster_id] = cluster_summary

    return results


def run_analysis(args):
    """运行分析主函数"""
    from datetime import datetime

    print("=" * 60)
    print("富集分析脚本 - 分析模块")
    print("=" * 60)
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"输入文件: {os.path.basename(args.input)}")

    output_dirs = setup_output_dirs(args.outdir)

    if args.verbose:
        print(f"创建输出目录: {args.outdir}")

    try:
        gene_sets = load_gene_sets(args)
    except FileNotFoundError as e:
        print(f"错误: {e}")
        sys.exit(1)

    try:
        print(args.input)
        degs = pd.read_csv(args.input, sep="\t")
        if args.verbose:
            print(f"成功读取输入文件: {args.input}")
            print(f"数据形状: {degs.shape}")
            print(f"列名: {list(degs.columns)}")
    except Exception as e:
        print(f"读取输入文件失败: {e}")
        sys.exit(1)

    required_columns = [args.gene_column, args.cluster_column]
    missing_columns = [col for col in required_columns if col not in degs.columns]
    if missing_columns:
        print(f"错误: 输入文件缺少必要的列: {missing_columns}")
        print(f"可用的列: {list(degs.columns)}")
        sys.exit(1)

    if args.clusters:
        clusters_to_analyze = [int(c.strip()) for c in args.clusters.split(",")]
        clusters_to_analyze = [
            c for c in clusters_to_analyze if c in degs[args.cluster_column].unique()
        ]
    else:
        clusters_to_analyze = sorted(degs[args.cluster_column].unique())

    print(f"将分析以下clusters: {clusters_to_analyze}")

    summary_info = {}

    print("\n开始富集分析...")

    for cluster_id in clusters_to_analyze:
        print(f"\n分析Cluster {cluster_id}...")

        cluster_degs = degs[degs[args.cluster_column] == cluster_id]
        gene_list = cluster_degs[args.gene_column].dropna().unique().tolist()

        if len(gene_list) < args.min_genes:
            print(
                f"  警告: Cluster {cluster_id} 的基因数({len(gene_list)})少于最小要求({args.min_genes})，跳过"
            )
            continue

        analyze_cluster(
            cluster_id, gene_list, gene_sets, output_dirs, args, summary_info
        )
        print(f"  Cluster {cluster_id} 分析完成")

    # 保存分析摘要信息
    summary_file = os.path.join(output_dirs["root"], "summary.json")
    import json

    # 将numpy类型转换为Python原生类型以便JSON序列化
    def convert_to_serializable(obj):
        if isinstance(
            obj,
            (
                np.int8,
                np.int16,
                np.int32,
                np.int64,
                np.uint8,
                np.uint16,
                np.uint32,
                np.uint64,
            ),
        ):
            return int(obj)
        elif isinstance(obj, (np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            # 处理字典键可能是numpy类型的情况
            new_dict = {}
            for k, v in obj.items():
                # 将键转换为字符串（如果它是numpy类型）
                if isinstance(
                    k,
                    (
                        np.int8,
                        np.int16,
                        np.int32,
                        np.int64,
                        np.uint8,
                        np.uint16,
                        np.uint32,
                        np.uint64,
                    ),
                ):
                    new_key = str(int(k))
                else:
                    new_key = str(k)
                new_dict[new_key] = convert_to_serializable(v)
            return new_dict
        elif isinstance(obj, list):
            return [convert_to_serializable(item) for item in obj]
        else:
            return obj

    serializable_summary = convert_to_serializable(summary_info)
    with open(summary_file, "w") as f:
        json.dump(serializable_summary, f, indent=2, ensure_ascii=False)

    print(f"\n分析完成！结果保存在: {args.outdir}")
    print(f"分析摘要保存为: {summary_file}")

    return summary_info, output_dirs


if __name__ == "__main__":
    args = parse_analysis_arguments()
    run_analysis(args)
