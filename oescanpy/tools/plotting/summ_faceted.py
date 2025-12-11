from typing import Optional, Union, List, Any
import scanpy as sc
from anndata import AnnData
from matplotlib import rc_context
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from matplotlib.legend import Legend
from matplotlib.axes._axes import Axes
from oescanpy.tools.color.utils import get_color_order
from oescanpy.tools.utils import save_figure
from pandas.api.types import CategoricalDtype
import os
import math
import oescanpy as oep
def summ_faceted_plot(
    adata ,
    groupby_all: str="group",  # 实际分组列（如 "group"）
    groupby: str = "sampleid",  # 样本 ID 列
    pointsize: int = 5,
    color: str = "clusters",
    palette: str = None,
    figsize: tuple = (15, 10),
    title: Optional[str] = None,
    save: Union[bool, str, None] = None,
    show: bool = True,
    outdir: Optional[str] = None,
    method: str = "umap"
):
    if groupby_all not in adata.obs.columns:
        raise ValueError(f"分组列 '{groupby_all}' 不存在于 adata.obs 中")
    if groupby not in adata.obs.columns:
        raise ValueError(f"样本列 '{groupby}' 不存在于 adata.obs 中")
    if color not in adata.obs.columns and color not in adata.var_names:
        raise ValueError(f"着色依据 '{color}' 不存在于 obs 或 var 中")
    # 判断是否为 groupby_all == groupby 的情况
    if groupby_all == groupby:
        # 直接按 groupby 的唯一值绘制每行2列的分面图
        return _plot_single_level(adata, groupby, pointsize, color, palette, save, show, outdir, method)
    # 获取 group 列的 level 顺序
    if isinstance(adata.obs[groupby_all].dtype, CategoricalDtype):
        groups = adata.obs[groupby_all].cat.categories.tolist()
    else:
        groups = sorted(adata.obs[groupby_all].unique())
    print(groups)
    # 获取 sampleid 的有序列表
    if isinstance(adata.obs[groupby].dtype, CategoricalDtype):
        all_samples = adata.obs[groupby].cat.categories.tolist()
    else:
        all_samples = sorted(adata.obs[groupby].unique())
    v = adata.obs[groupby].unique()
    n_groups = len(groupby)
    # 统计每个 group 的 sample 数量
    group_samples = {
        group: adata[adata.obs[groupby_all] == group].obs[groupby].unique().tolist()
        for group in groups
    }
    print(group_samples)
    # 获取每个 group 的 sample 数量及最大值
    group_sample_counts = {k: len(v) for k, v in group_samples.items()}
    max_samples = max(group_sample_counts.values())
    layout_mode = "group_level" if max_samples <= 6 else "sample_level"
    unique_values = adata.obs[color].unique()
    max_nchar = max(len(str(x)) for x in unique_values)
    nlevel = len(unique_values)
    legend_ncols = math.ceil(nlevel / 15) if nlevel >= 1 else 1
    print(layout_mode)
    color_dict = get_color_order(adata, color, palette)
    # 获取有序的分类标签
    if isinstance(adata.obs[color].dtype, CategoricalDtype):
        ordered_categories = adata.obs[color].cat.categories.tolist()
    else:
        ordered_categories = sorted(adata.obs[color].unique())
    # 计算图例列数
    nlevel = len(ordered_categories)
    legend_ncols = math.ceil(nlevel / 15) if nlevel >= 1 else 1
    # 动态调整图像尺寸
    if layout_mode == "group_level":
        ncols = max_samples
        nrows = len(groups)
        subplot_size = 4
        width = subplot_size * ncols
        height = subplot_size * nrows
        max_label_length = max(len(str(cat)) for cat in ordered_categories)
        legend_col_width = (max_label_length * 0.008) * legend_ncols +0.1
        width = width + legend_col_width * subplot_size
        legend_x = 0.9*(width/(subplot_size * ncols)) + 0.2 - (nrows - 1) *0.01
        if max_label_length >40:
            legend_x = 0.95*(width/(subplot_size * ncols)) + 0.2 - (nrows - 1) *0.01
        print(legend_x)
    all_x = []
    all_y = []
    buffer_ratio = 0.1
    for group in groups:
        samples = group_samples[group]
        for sample in samples:
            sample_data = adata[adata.obs[groupby] == sample]
            if method == "umap":
                x = sample_data.obsm["X_umap"][:, 0]
                y = sample_data.obsm["X_umap"][:, 1]
            elif method == "tsne":
                x = sample_data.obsm["X_tsne"][:, 0]
                y = sample_data.obsm["X_tsne"][:, 1]
            else:
                raise ValueError(f"Unsupported method: {method}")
            all_x.extend(x)
            all_y.extend(y)
    x_span = np.max(all_x) - np.min(all_x)
    y_span = np.max(all_y) - np.min(all_y)
    max_span = max(x_span, y_span)
    max_span_with_buffer = max_span * (1 + buffer_ratio)
    center_x = (np.min(all_x) + np.max(all_x)) / 2
    center_y = (np.min(all_y) + np.max(all_y)) / 2
    x_range = (center_x - max_span_with_buffer / 2, center_x + max_span_with_buffer / 2)
    y_range = (center_y - max_span_with_buffer / 2, center_y + max_span_with_buffer / 2)
    if layout_mode == "sample_level" or max_samples == 1:
        return _plot_single_level(adata, groupby, pointsize, color, palette, save, show, outdir, method,ncols = 6)
    with rc_context({'figure.figsize': (width, height)}):
        fig, axes = plt.subplots(nrows, ncols, squeeze=False)
        axes = axes.flat
        # 绘制 group_level 布局
        if layout_mode == "group_level":
            ax_idx = 0
            all_samples = 0
            for group in groups:
                samples = group_samples[group]
                n_samples = len(samples)
                all_samples = all_samples +n_samples
                group_axes = axes[ax_idx:ax_idx + ncols]
                for ax,sample in zip(group_axes,samples):
                    sample_data = adata[adata.obs[groupby] == sample]
                    oep.plot_dim(sample_data, color=color, palette=color_dict, ax=ax, show=False, legend_loc='none',size=pointsize,method=method)
                    ax.set_xlim(*x_range)
                    ax.set_ylim(*y_range)
                    ax.set_title(f"{groupby}: {sample}",fontsize=8)
                    ax.set_aspect('equal')
                for ax in group_axes[n_samples:]:
                    ax.axis('off')
                ax_idx += ncols
        # 添加公共图例
        ordered_handles = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_dict[cat], markersize=8)
            for cat in ordered_categories
        ]
        right_pad = legend_col_width / width  # 计算右侧空白比例
        fig.subplots_adjust(right=1 - right_pad * 1.5) 
        title = ''
        if color == "clusters":
            title = color
        fig.legend(
                handles=ordered_handles, 
                labels=ordered_categories,
                loc='center right', 
                bbox_to_anchor=(legend_x, 0.5), 
                title=title,
                frameon=False,
                borderaxespad=0.,
               ncol=legend_ncols)
        plot_file = os.path.join(outdir, f'splitby-{groupby}_split_plot')
        save_figure(plot_file)
        return fig
            
def _plot_single_level(
    adata,
    groupby: str,
    pointsize: int = 5,
    color: str = "clusters",
    palette: str = None,
    save: Union[bool, str, None] = None,
    show: bool = True,
    outdir: Optional[str] = None,
    method: str = "umap",
	ncols: int =6,
):
    """
    当 groupby_all == groupby 时，直接按 groupby 的唯一值绘制每行6列的分面图
    """
    # 获取 groupby 的唯一值列表
    if isinstance(adata.obs[groupby].dtype, CategoricalDtype):
        samples = adata.obs[groupby].cat.categories.tolist()
    else:
        samples = sorted(adata.obs[groupby].unique())
    total_samples = len(samples)
    if total_samples < 6:
        ncols = total_samples
    nrows = (total_samples + ncols - 1) // ncols
    unique_values = adata.obs[color].unique()
    max_nchar = max(len(str(x)) for x in unique_values)
    nlevel = len(unique_values)
    legend_ncols = math.ceil(nlevel / 15) if nlevel >= 1 else 1
    subplot_size = 4
    width = subplot_size * ncols 
    height = subplot_size * nrows
    color_dict = get_color_order(adata, color, palette)
    if isinstance(adata.obs[color].dtype, CategoricalDtype):
        ordered_categories = adata.obs[color].cat.categories.tolist()
    else:
        ordered_categories = sorted(adata.obs[color].unique())
        max_label_length = max(len(str(cat)) for cat in ordered_categories)
    max_label_length = max(len(str(cat)) for cat in ordered_categories)
    legend_col_width = (max_label_length * 0.008) * legend_ncols +0.1
    width = width + legend_col_width * subplot_size
    legend_x = 0.9*(width/(subplot_size * ncols)) + 0.2 - (nrows - 1) *0.01
    if max_label_length >40:
        legend_x = 0.95 * (width/(subplot_size * ncols)) + 0.2 - (nrows - 1) *0.01
    print(legend_x)
    with rc_context({'figure.figsize': (width, height)}):
        fig, axes = plt.subplots(nrows, ncols, squeeze=False,constrained_layout=False,gridspec_kw={'wspace': 0.05})
        axes = axes.flat
        # 获取所有样本的坐标范围（用于统一坐标轴）
        buffer_ratio = 0.1
        all_x = []
        all_y = []
        # 第一次遍历：收集所有坐标范围
        for sample in samples:
            sample_data = adata[adata.obs[groupby] == sample]
            if method == "umap":
                x = sample_data.obsm["X_umap"][:, 0]
                y = sample_data.obsm["X_umap"][:, 1]
            elif method == "tsne":
                x = sample_data.obsm["X_tsne"][:, 0]
                y = sample_data.obsm["X_tsne"][:, 1]
            else:
                raise ValueError(f"Unsupported method: {method}")
            all_x.extend(x)
            all_y.extend(y)
        x_span = np.max(all_x) - np.min(all_x)
        y_span = np.max(all_y) - np.min(all_y)
        max_span = max(x_span, y_span)
        max_span_with_buffer = max_span * (1 + buffer_ratio)
        center_x = (np.min(all_x) + np.max(all_x)) / 2
        center_y = (np.min(all_y) + np.max(all_y)) / 2
        x_range = (center_x - max_span_with_buffer / 2, center_x + max_span_with_buffer / 2)
        y_range = (center_y - max_span_with_buffer / 2, center_y + max_span_with_buffer / 2)
        for ax,sample in zip(axes, samples):
            sample_data = adata[adata.obs[groupby] == sample]
            oep.plot_dim(
                sample_data,
                color=color,
                palette=color_dict,
                ax=ax,
                show=False,
                legend_loc='none',
                size=pointsize,
                method=method
            )
            ax.set_title(f"{groupby}: {sample}")
            ax.set_xlim(*x_range)
            ax.set_ylim(*y_range)
            ax.set_aspect('equal')
        # 隐藏多余子图
        for ax in axes[total_samples:]:
            ax.axis('off')
        # 构建图例
        ordered_handles = [
            plt.Line2D([0], [0], marker='o', color='w',
                       markerfacecolor=color_dict[cat], markersize=8)
            for cat in ordered_categories
        ]
        right_pad = legend_col_width / width  # 计算右侧空白比例
        fig.subplots_adjust(right=1 - right_pad * 1.5)
        title = ''
        if color == "clusters":
            title = color
        # 添加图例
        fig.legend(
            handles=ordered_handles,
            labels=ordered_categories,
            loc='center right',
            title=title,
            frameon=False,
            borderaxespad=0.,
            bbox_to_anchor=(legend_x,0.5),
            ncol=legend_ncols
        )
        plot_file = os.path.join(outdir, f'splitby-{groupby}_split_plot')
        save_figure(plot_file)
        return fig