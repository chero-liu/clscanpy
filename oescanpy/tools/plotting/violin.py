from typing import List, Dict, Optional, Tuple, Union
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import scanpy as sc
import oescanpy as oep


def plot_qc_violin(
    adata,
    groupby: str,
    metrics: str,
    palette: Optional[Dict[str, str]] = None,
    figsize: Tuple[int, int] = (25, 10),
    dpi: int = 300,
    save_path: Optional[str] = None,
    save_formats: Tuple[str, ...] = ("png", "pdf"),
    groupby_levels: Union[str, list] = None,
):
    adata = oep.loadH5AD(adata, groupby=groupby, groupby_levels=groupby_levels)
    sorted_groups = adata.obs[groupby].cat.categories.tolist()
    max_len = max(len(group) for group in sorted_groups)
    n_samples = len(sorted_groups)
    # 自动计算子图布局（2行自动分配列数）
    n_metrics = len(metrics)
    n_cols = (n_metrics + 1) // 2  # 向上取整分配列数
    width = 2 + n_samples * 3
    height = max_len*0.4+6
    # 创建画布
    fig, axes = plt.subplots(2, n_cols, figsize=(width, height), dpi=dpi)
    axes = axes.flatten()  # 展平坐标轴数组便于索引
    # 绘制每个指标
    for idx, metric in enumerate(metrics):
        if idx < n_metrics:
            sc.pl.violin(
                adata,
                keys=metric,
                groupby=groupby,
                jitter=0.2,
                rotation=90,
                multi_panel=False,
                ax=axes[idx],
                show=False,
                order=sorted_groups,
                palette=palette,
            )
            axes[idx].set_xlabel("")
            axes[idx].set_ylabel("")
            axes[idx].tick_params(axis="both", which="major", labelsize=16)
            axes[idx].set_title(
                f"{metrics[idx]}",
                fontsize=20,  # 标题字体大小
                fontweight="bold",  # 加粗
                pad=20,  # 标题与图表间距
            )
            plt.setp(
                axes[idx].get_xticklabels(),
                rotation=45,  # 45度旋转
                ha="right",  # 右对齐
                rotation_mode="anchor",
            )  # 防止文字截断
    for ax in axes[n_metrics:]:
        ax.axis("off")
    # 自动调整布局
    plt.tight_layout()
    # 保存结果
    if save_path:
        for fmt in save_formats:
            plt.savefig(f"{save_path}.{fmt}", dpi=dpi, bbox_inches="tight")
    plt.show()
    plt.close(fig)
