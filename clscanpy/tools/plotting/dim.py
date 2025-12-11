from typing import Optional, Union, List, Dict, Any
from anndata import AnnData
import scanpy as sc
from matplotlib import rc_context

def plot_dim(
    adata: AnnData,
    method: str = "umap",
    figsize: Optional[tuple] = (5, 5),
    color: Union[str, List[str], None] = None,
    mask_obs: Union[str, List[str], None] = None,
    gene_symbols: Optional[str] = None,
    use_raw: Optional[bool] = None,
    sort_order: bool = True,
    edges: bool = False,
    edges_width: float = 0.1,
    edges_color: str = "grey",
    neighbors_key: Optional[str] = None,
    arrows: bool = False,
    arrows_kwds: Optional[Dict[str, Any]] = None,
    groups: Union[str, List[str], None] = None,
    components: Union[str, List[str], None] = None,
    dimensions: Union[int, List[int], None] = None,
    layer: Optional[str] = None,
    projection: str = "2d",
    scale_factor: Optional[float] = None,
    color_map: Optional[str] = None,
    cmap: Optional[str] = None,
    palette: Optional[str] = None,
    na_color: str = "lightgray",
    na_in_legend: bool = True,
    size: Optional[float] = None,
    frameon: Optional[bool] = None,
    legend_fontsize: Optional[int] = None,
    legend_fontweight: str = "bold",
    legend_loc: str = "right margin",
    legend_fontoutline: Optional[int] = None,
    colorbar_loc: str = "right",
    vmax: Union[float, List[float], None] = None,
    vmin: Union[float, List[float], None] = None,
    vcenter: Union[float, List[float], None] = None,
    norm: Optional[Any] = None,
    add_outline: bool = False,
    outline_width: tuple = (0.3, 0.05),
    outline_color: tuple = ("black", "white"),
    ncols: int = 4,
    hspace: float = 0.25,
    wspace: Optional[float] = None,
    title: Optional[str] = None,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    ax: Optional[Any] = None,
    return_fig: Optional[bool] = None,
    marker: str = ".",
    **kwargs
) -> Optional[Any]:
    """\
    Wrap function for visualizing UMAP/t-SNE embeddings.

    Parameters
    ----------
    method
        Visualization method to use: 'umap' (default) or 'tsne'.
    adata
        Annotated data matrix with embeddings in `adata.obsm['X_umap']` or `adata.obsm['X_tsne']`.
    Returns
    -------
    If `return_fig` is True, returns a :class:`matplotlib.figure.Figure`.
    Otherwise shows the plot and returns None.

    Example
    -------
    >>> # Plot UMAP
    >>> plot_dim(adata, method='umap', color='louvain')
    >>> # Plot t-SNE
    >>> plot_dim(adata, method='tsne', color='CD79A', use_raw=True)
    """

    if method not in ["umap", "tsne"]:
        raise ValueError("method must be 'umap' or 'tsne'")
    plot_func = getattr(sc.pl, method)
    with rc_context(rc={'figure.figsize': figsize}):
        plot_func(
            adata,
            color=color,
            mask_obs=mask_obs,
            gene_symbols=gene_symbols,
            use_raw=use_raw,
            sort_order=sort_order,
            edges=edges,
            edges_width=edges_width,
            edges_color=edges_color,
            neighbors_key=neighbors_key,
            arrows=arrows,
            arrows_kwds=arrows_kwds,
            groups=groups,
            components=components,
            dimensions=dimensions,
            layer=layer,
            projection=projection,
            scale_factor=scale_factor,
            color_map=color_map,
            cmap=cmap,
            palette=palette,
            na_color=na_color,
            na_in_legend=na_in_legend,
            size=size,
            frameon=frameon,
            legend_fontsize=legend_fontsize,
            legend_fontweight=legend_fontweight,
            legend_loc=legend_loc,
            legend_fontoutline=legend_fontoutline,
            colorbar_loc=colorbar_loc,
            vmax=vmax,
            vmin=vmin,
            vcenter=vcenter,
            norm=norm,
            add_outline=add_outline,
            outline_width=outline_width,
            outline_color=outline_color,
            ncols=ncols,
            hspace=hspace,
            wspace=wspace,
            title=title,
            show=show,
            save=save,
            ax=ax,
            return_fig=return_fig,
            marker=marker,
            **kwargs
        )
