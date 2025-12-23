#!/usr/bin/env python3
"""
报告脚本 - 生成富集分析报告
"""

import argparse
import pandas as pd
import os
import sys
import json
import configparser
from typing import Dict, List
from datetime import datetime


def parse_report_arguments():
    """解析报告脚本的命令行参数"""
    parser = argparse.ArgumentParser(
        description="生成富集分析报告",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # 必需参数
    parser.add_argument("--input-dir", required=True, help="分析结果目录路径")
    parser.add_argument("--output-dir", required=True, help="报告输出目录路径")
    parser.add_argument(
        "--summary-file", required=True, help="分析摘要文件路径 (summary.json)"
    )
    parser.add_argument("--bp", required=True, help="GO BP基因集文件路径")
    parser.add_argument("--cc", required=True, help="GO CC基因集文件路径")
    parser.add_argument("--mf", required=True, help="GO MF基因集文件路径")
    parser.add_argument("--kegg", required=True, help="KEGG基因集文件路径")

    # 可选参数
    parser.add_argument(
        "--pvalue", type=float, default=0.05, help="p值阈值 (默认: 0.05)"
    )
    parser.add_argument("--verbose", action="store_true", help="显示详细输出信息")

    return parser.parse_args()


def load_analysis_arguments(input_dir: str) -> Dict:
    """加载分析参数文件"""
    arguments_file = os.path.join(input_dir, "arguments.txt")
    if not os.path.exists(arguments_file):
        print(f"警告: 未找到参数文件 {arguments_file}")
        return {}

    try:
        config = configparser.ConfigParser()
        # 保留原始大小写
        config.optionxform = str
        config.read(arguments_file)

        arguments = {}
        # 从所有section中提取参数
        for section in config.sections():
            for key, value in config.items(section):
                arguments[key] = value

        return arguments
    except Exception as e:
        print(f"警告: 无法解析参数文件 {arguments_file}: {e}")
        return {}


def load_analysis_summary(summary_file: str) -> Dict:
    """加载分析摘要信息"""
    try:
        with open(summary_file, "r") as f:
            summary = json.load(f)
        return summary
    except Exception as e:
        print(f"错误: 无法加载分析摘要文件 {summary_file}: {e}")
        sys.exit(1)


def create_summary_report(
    summary_data: Dict, output_dir: str, args, input_file: str = None
):
    """创建分析概览报告"""
    # 加载分析参数
    analysis_args = load_analysis_arguments(args.input_dir)

    summary_file = os.path.join(output_dir, "report.txt")

    total_clusters = len(summary_data)
    total_genes = sum([info["gene_count"] for info in summary_data.values()])

    category_stats = {"BP": 0, "CC": 0, "MF": 0, "KEGG": 0}
    clusters_with_enrichment = []

    for cluster_id, info in summary_data.items():
        has_enrichment = False
        for category in category_stats.keys():
            if info["significant_counts"].get(category, 0) > 0:
                category_stats[category] += 1
                has_enrichment = True
        if has_enrichment:
            clusters_with_enrichment.append(cluster_id)

    best_clusters = []
    if summary_data:
        max_significant = max(
            [sum(info["significant_counts"].values()) for info in summary_data.values()]
        )
        if max_significant > 0:
            best_clusters = [
                cluster_id
                for cluster_id, info in summary_data.items()
                if sum(info["significant_counts"].values()) == max_significant
            ]

    with open(summary_file, "w", encoding="utf-8") as f:
        f.write(
            "╔══════════════════════════════════════════════════════════════════════════╗\n"
        )
        f.write(
            "║                             富集分析报告                                 ║\n"
        )
        f.write(
            "╚══════════════════════════════════════════════════════════════════════════╝\n\n"
        )

        # 1. 报告基本信息
        f.write("一、 报告基本信息\n")
        f.write("    " + "─" * 70 + "\n")
        f.write(f"    分析类型: GO和KEGG富集分析\n")
        if input_file:
            f.write(f"    输入文件: {os.path.basename(input_file)}\n")
        f.write("\n")

        f.write("二、 分析参数设置\n")
        f.write("    " + "─" * 70 + "\n")
        f.write("    ▪ 显著性阈值: p < {:.3f}\n".format(args.pvalue))

        # 从分析参数中获取额外的参数
        min_genes = analysis_args.get("min_genes", "3")
        max_genes = analysis_args.get("max_genes", "500")
        top_n = analysis_args.get("top_n", "20")
        truncate_length = analysis_args.get("truncate_length", "30")

        f.write(f"    ▪ 基因数范围: {min_genes} ≤ 基因集大小 ≤ {max_genes}\n")
        f.write(f"    ▪ 每类显示top条目数: {top_n}\n")
        f.write(
            f"    ▪ 术语截断长度(便于展示，超过{truncate_length}字符的通路使用...代替, 原始通路名称于raw_Term中): \n"
        )

        f.write("    ▪ 使用基因集:\n")
        f.write("        • GO生物过程 (BP): {}\n".format(os.path.basename(args.bp)))
        f.write("        • GO细胞组分 (CC): {}\n".format(os.path.basename(args.cc)))
        f.write("        • GO分子功能 (MF): {}\n".format(os.path.basename(args.mf)))
        f.write("        • KEGG通路: {}\n".format(os.path.basename(args.kegg)))
        f.write("\n")

        f.write("三、 分析总体概况\n")
        f.write("    " + "─" * 70 + "\n")
        f.write(f"    1. 共分析了 {total_clusters} 个细胞群体 (Clusters)\n")
        f.write(f"    2. 涉及差异表达基因 {total_genes} 个\n")
        f.write(f"    3. {len(clusters_with_enrichment)} 个Clusters有显著富集结果\n")
        if clusters_with_enrichment:
            f.write(
                f"       - 分别为: Cluster "
                + ", ".join(map(str, clusters_with_enrichment))
                + "\n"
            )

        f.write("    4. 各功能类别富集情况:\n")
        for category, count in category_stats.items():
            percentage = (count / total_clusters * 100) if total_clusters > 0 else 0
            category_name = {
                "BP": "GO生物过程",
                "CC": "GO细胞组分",
                "MF": "GO分子功能",
                "KEGG": "KEGG通路",
            }[category]
            f.write(
                f"       • {category_name}: {count} 个Clusters有显著富集 ({percentage:.1f}%)\n"
            )
        f.write("\n")

        f.write("四、 各Cluster富集分析详细结果\n")
        f.write("    " + "─" * 70 + "\n")

        if not summary_data:
            f.write("    本次分析未获得显著富集结果。\n")
        else:
            f.write(
                "    ┌─────────┬────────────┬─────────────────────────────────────────────┐\n"
            )
            f.write(
                "    │ Cluster │ 基因数     │ 显著富集条目数 (p<{:.3f})                    │\n".format(
                    args.pvalue
                )
            )
            f.write(
                "    │         │            ├────────┬────────┬────────┬────────┬─────────┤\n"
            )
            f.write(
                "    │         │            │  BP   │  CC    │  MF    │ KEGG   │  总计   │\n"
            )
            f.write(
                "    ├─────────┼────────────┼────────┼────────┼────────┼────────┼─────────┤\n"
            )

            for cluster_id in sorted(summary_data.keys()):
                info = summary_data[cluster_id]
                bp_count = info["significant_counts"].get("BP", 0)
                cc_count = info["significant_counts"].get("CC", 0)
                mf_count = info["significant_counts"].get("MF", 0)
                kegg_count = info["significant_counts"].get("KEGG", 0)
                total_count = bp_count + cc_count + mf_count + kegg_count

                f.write(
                    "    │ {:^7} │ {:^10} │ {:^6} │ {:^6} │ {:^6} │ {:^6} │ {:^7} │\n".format(
                        cluster_id,
                        info["gene_count"],
                        bp_count,
                        cc_count,
                        mf_count,
                        kegg_count,
                        total_count,
                    )
                )

            f.write(
                "    └─────────┴────────────┴────────┴────────┴────────┴────────┴─────────┘\n"
            )

        f.write("\n")

        f.write("五、 重点发现\n")
        f.write("    " + "─" * 70 + "\n")

        if best_clusters:
            f.write("    1. 富集最显著的Clusters:\n")
            for cluster_id in best_clusters:
                info = summary_data[cluster_id]
                total_sig = sum(info["significant_counts"].values())
                f.write(f"       • Cluster {cluster_id}: {total_sig} 个显著富集条目\n")
                f.write(f"         包含: {info['gene_count']} 个差异表达基因\n")
                for category in ["BP", "CC", "MF", "KEGG"]:
                    if info["significant_counts"].get(category, 0) > 0:
                        result_file = os.path.join(
                            args.input_dir,
                            "tables",
                            f"cluster_{cluster_id}_{category}_enrichment.tsv",
                        )
                        if os.path.exists(result_file):
                            try:
                                df = pd.read_csv(result_file, sep="\t")
                                df_sig = df[df["Adjusted P-value"] <= args.pvalue]
                                df_sig = df_sig.sort_values("Adjusted P-value").head(3)
                                if not df_sig.empty:
                                    category_name = {
                                        "BP": "生物过程",
                                        "CC": "细胞组分",
                                        "MF": "分子功能",
                                        "KEGG": "通路",
                                    }[category]
                                    f.write(f"         {category_name}:\n")
                                    for _, row in df_sig.iterrows():
                                        term = row.get("Term", "")
                                        pval = row.get("Adjusted P-value", "")
                                        if term:
                                            f.write(
                                                f"           - {term} (p={pval:.2e})\n"
                                            )
                            except:
                                pass

        if clusters_with_enrichment:
            avg_genes = sum(
                [summary_data[c]["gene_count"] for c in clusters_with_enrichment]
            ) / len(clusters_with_enrichment)
            f.write(f"\n    2. 富集统计:\n")
            f.write(
                f"       • 平均每个有富集的Cluster包含 {avg_genes:.1f} 个差异表达基因\n"
            )

            max_category = max(category_stats.items(), key=lambda x: x[1])
            if max_category[1] > 0:
                category_name = {
                    "BP": "GO生物过程",
                    "CC": "GO细胞组分",
                    "MF": "GO分子功能",
                    "KEGG": "KEGG通路",
                }[max_category[0]]
                f.write(
                    f"       • {category_name} 在最多Clusters中富集 ({max_category[1]}个Clusters)\n"
                )

        f.write("\n")

        f.write("六、 结果文件说明\n")
        f.write("    " + "─" * 70 + "\n")
        f.write("    本次分析生成了以下结果文件，可用于进一步的数据挖掘和可视化：\n\n")

        f.write("    1. 数据表格文件 (tables/目录):\n")
        f.write("       • cluster_[编号]_[类别]_enrichment.tsv - 原始富集分析结果\n")
        f.write("       • cluster_[编号]_GO_combined.tsv - GO三个类别的合并结果\n")
        f.write("       \n")
        f.write("       表格列名说明:\n")
        f.write("       - Term: 功能术语或通路名称\n")
        f.write("       - P-value: 原始富集显著性p值\n")
        f.write("       - Adjusted P-value: 多重检验校正后的p值 (FDR)\n")
        f.write("       - Odds Ratio: 比值比，表示富集程度\n")
        f.write("       - Combined Score: 综合评分\n")
        f.write("       - Genes: 输入基因中富集到的基因列表\n")
        f.write("       - Gene Ratio: 富集基因比例 (k/K)\n")
        f.write("       - Bg Ratio: 背景基因比例 (n/N)\n")
        f.write("       - Category: 功能类别 (BP/CC/MF/KEGG)\n")
        f.write("       - Gene Count: 富集到的基因数\n")
        f.write("       \n")

        f.write("    2. 可视化图形文件 (figures/目录):\n")
        f.write("       • *_GO_barplot.png/pdf - GO富集条形图（按类别分组）\n")
        f.write("         - X轴: 富集分数 (-log10(p-value))\n")
        f.write("         - Y轴: GO术语名称（按类别分组）\n")
        f.write("         - 颜色: 不同GO类别 (BP/CC/MF)\n")
        f.write("         - 条形长度: 富集显著性程度\n")
        f.write("         - 说明: 展示每个类别top N显著富集结果\n")
        f.write("       \n")

        f.write("       • *_GO_dotplot.png/pdf - GO富集点图（按类别分组）\n")
        f.write("         - X轴: 基因比例 (Gene Ratio)\n")
        f.write("         - Y轴: GO术语名称（按类别分组）\n")
        f.write("         - 点大小: 富集到的基因数\n")
        f.write("         - 点颜色: -log10(p-value)，颜色越深表示越显著\n")
        f.write("         - 说明: 同时展示富集显著性和富集强度\n")
        f.write("       \n")

        f.write("       • *_enrichment_dotplot.png/pdf - 综合富集点图\n")
        f.write("         - X轴: Cluster编号\n")
        f.write("         - Y轴: 功能术语名称\n")
        f.write("         - 点大小: 富集到的基因数\n")
        f.write("         - 点颜色: -log10(p-value)，颜色越深表示越显著\n")
        f.write("         - 说明: 跨Cluster的富集结果比较\n")
        f.write("       \n")

        f.write("       • *_KEGG_barplot.png/pdf - KEGG通路富集条形图\n")
        f.write("         - X轴: 富集分数 (-log10(p-value))\n")
        f.write("         - Y轴: KEGG通路名称\n")
        f.write("         - 条形长度: 富集显著性程度\n")
        f.write("         - 说明: 展示KEGG通路富集的显著性\n")
        f.write("       \n")

        f.write("       • *_KEGG_dotplot.png/pdf - KEGG通路富集点图\n")
        f.write("         - X轴: 基因比例 (Gene Ratio)\n")
        f.write("         - Y轴: KEGG通路名称\n")
        f.write("         - 点大小: 富集到的基因数\n")
        f.write("         - 点颜色: -log10(p-value)，颜色越深表示越显著\n")
        f.write("         - 说明: 展示KEGG通路的富集强度和显著性\n")
        f.write("       \n")

        f.write(
            "       格式说明: 提供PNG和PDF两种格式，PNG适合网页查看，PDF适合出版发表\n"
        )
        f.write("\n")

        f.write("    3. 本报告文件:\n")
        f.write("       • report.txt - 本文件，包含分析概览和重点发现\n")
        f.write("\n")

        f.write("七、 分析建议\n")
        f.write("    " + "─" * 70 + "\n")

        if clusters_with_enrichment:
            f.write("    1. 重点关注以下Clusters的富集结果:\n")
            f.write(f"       • {', '.join([f'Cluster {c}' for c in best_clusters])}\n")
            f.write("         这些Clusters显示了最显著的生物学功能富集\n\n")

            f.write("    2. 后续分析建议:\n")
            f.write("       • 结合富集结果和基因表达模式，深入理解细胞群体功能特性\n")
            f.write("       • 使用KEGG通路富集结果，构建信号通路网络图\n")
            f.write("       • 对于关键生物学过程，可进一步进行通路活性分析\n")
            f.write("       • 可使用表格中的基因列表进行蛋白互作网络分析\n")
        else:
            f.write("    本次分析未发现显著富集的生物学功能或通路。\n")
            f.write("    可能原因：\n")
            f.write("       • 使用的基因集与物种/组织不匹配\n")
            f.write("       • p值阈值设置过于严格\n")
            f.write("       • 差异表达基因数量不足\n")
            f.write("    建议：\n")
            f.write("       • 检查基因集文件是否与实验物种匹配\n")
            f.write("       • 适当放宽p值阈值（如改为0.1）\n")
            f.write("       • 考虑使用更宽松的差异表达筛选标准\n")
            f.write("       • 调整基因集大小范围 (min-genes/max-genes)\n")

        f.write("\n")

        f.write(
            "╔══════════════════════════════════════════════════════════════════════════╗\n"
        )
        f.write(
            "║                              报告结束                                     ║\n"
        )
        f.write(
            "║                                                                          ║\n"
        )
        f.write(
            "║   注：本报告为自动生成，生物学解释需结合具体实验背景和专业知识                ║\n"
        )
        f.write(
            "╚══════════════════════════════════════════════════════════════════════════╝\n"
        )

    return summary_file


def run_report(args, input_file: str = None):
    """运行报告生成主函数"""
    from datetime import datetime

    print("=" * 60)
    print("富集分析脚本 - 报告模块")
    print("=" * 60)
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"输入目录: {args.input_dir}")

    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)

    if args.verbose:
        print(f"创建报告输出目录: {args.output_dir}")

    # 加载分析参数
    analysis_args = load_analysis_arguments(args.input_dir)
    if args.verbose and analysis_args:
        print(f"加载分析参数: {len(analysis_args)} 个参数")

    # 加载分析摘要
    summary_data = load_analysis_summary(args.summary_file)

    if args.verbose:
        print(f"加载分析摘要: {len(summary_data)} 个clusters")

    # 生成报告
    print(f"\n开始生成分析报告...")
    report_file = create_summary_report(summary_data, args.output_dir, args, input_file)
    print(f"报告生成完成！保存为: {report_file}")


if __name__ == "__main__":
    args = parse_report_arguments()
    run_report(args)
