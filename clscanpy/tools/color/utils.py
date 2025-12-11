from typing import Dict
from anndata import AnnData
import warnings
from clscanpy.tools.color.colors import add_color_to_dataframe

import clscanpy as oep


def get_color_order(
    adata: AnnData,
    use_col: str,
    palette: str = None,
    use_col_levels: str = None,
) -> Dict[str, str]:
    """
    获取或创建指定分类列的颜色映射字典

    参数:
        adata: AnnData对象
        use_col: 需要颜色映射的元数据列名
        palette: 调色板名称或文件路径.csv,两列, 列名为adata.obs中的列名

    返回:
        颜色映射字典 {分类值: HEX颜色码}

    异常:
        ValueError: 输入列不存在或调色板无效
        TypeError: 输入对象不是AnnData
    """
    color_col = f"{use_col}_col"
    adata = oep.loadH5AD(
        adata,
        groupby=use_col,
        groupby_levels=use_col_levels,
    )

    if color_col not in adata.obs.columns:
        adata.obs = add_color_to_dataframe(
            adata.obs,
            colname=use_col,
            palette=palette,
        )
    elif color_col in adata.obs.columns and palette != None:
        del adata.obs[color_col]
        adata.obs = add_color_to_dataframe(
            adata.obs,
            colname=use_col,
            palette=palette,
        )
    categories = adata.obs[use_col].dtype.categories
    colors = []
    for cat in categories:
        mask = adata.obs[use_col] == cat
        color_val = adata.obs.loc[mask, color_col].iloc[0]
        colors.append(color_val)

    return dict(zip(categories, colors))
