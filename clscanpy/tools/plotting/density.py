from typing import Optional, Union, List, Any
import scanpy as sc
from anndata import AnnData
from matplotlib import rc_context

def plot_density(
    adata: AnnData,
    basis: str = 'umap',
    groupby: Optional[str] = None,
    figsize: tuple = (5, 5),
    key_added: Optional[str] = None,
    components: Optional[Union[str, List[str]]] = None,
    group: Union[str, List[str], None] = 'all',
    color_map: str = 'YlOrRd',
    bg_dotsize: int = 80,
    fg_dotsize: int = 180,
    vmin: float = 0,
    vmax: float = 1,
    vcenter: Optional[float] = None,
    norm: Optional[Any] = None,
    ncols: int = 4,
    hspace: float = 0.25,
    wspace: Optional[float] = None,
    title: Optional[str] = None,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    ax: Optional[Any] = None,
    return_fig: Optional[bool] = None,
    **kwargs
) -> Optional[Any]:
    """
    Wrapper function to compute and visualize embedding density.

    Parameters
    ----------
    adata
        Annotated data matrix.
    basis
        Embedding basis to use (default: 'umap').
    groupby
        Key for grouping cells into categories. If None, computes overall density.
    key_added
        Key under which the density estimates will be stored. If None, generates automatically.
    components
        Embedding dimensions to consider (must be 2 components).
    group
        Categories to plot. Use 'all' for all groups, None for overall density.
    color_map
        Color map for density visualization.
    bg_dotsize
        Dot size for background cells.
    fg_dotsize
        Dot size for foreground cells in the group.
    vmin, vmax, vcenter
        Color scale limits and center.
    norm
        Normalization for the color scale.
    ncols
        Number of columns for subplots.
    hspace, wspace
        Spacing between subplots.
    title
        Plot title.
    show
        Whether to display the plot.
    save
        Save the plot to a file.
    ax
        Matplotlib axis to plot on.
    return_fig
        Whether to return the figure object.
    **kwargs
        Additional keyword arguments passed to sc.pl.embedding_density.

    Returns
    -------
    If `return_fig` is True, returns the matplotlib figure, otherwise shows the plot.
    """

    # Check if embedding exists
    if f'X_{basis}' not in adata.obsm:
        raise ValueError(f"Embedding '{basis}' not found in adata.obsm. Compute it first.")

    # Generate default key_added if not provided
    if key_added is None:
        key_added = f"{basis}_density_{groupby}" if groupby is not None else f"{basis}_density"

    # Compute density if not already present
    if key_added not in adata.obs:
        sc.tl.embedding_density(
            adata,
            basis=basis,
            groupby=groupby,
            key_added=key_added,
            components=components
        )

    # Force group to None if no grouping was performed
    if groupby is None:
        group = None

    with rc_context(rc={'figure.figsize': figsize}):
        sc.pl.embedding_density(
            adata,
            basis=basis,
            key=key_added,
            group=group,
            color_map=color_map,
            bg_dotsize=bg_dotsize,
            fg_dotsize=fg_dotsize,
            vmin=vmin,
            vmax=vmax,
            vcenter=vcenter,
            norm=norm,
            ncols=ncols,
            hspace=hspace,
            wspace=wspace,
            title=title,
            show=show if show is not None else True,
            save=save,
            ax=ax,
            return_fig=return_fig,
            **kwargs
        )
