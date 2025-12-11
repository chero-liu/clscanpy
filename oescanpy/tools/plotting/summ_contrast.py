import oescanpy as oep
import anndata as ad
from oescanpy.tools.color.utils import get_color_order
import os
from oescanpy.tools.utils import save_figure
def summ_contrast_plot(
    adata,
    pointsize: int = 5,
    color: str = "sampleid",
    palette: str = None,
    title: str = None,
    show: bool = True,
	outdir: str = None,
	method: str = "umap"
):
    color_use = get_color_order(adata,color,palette)
    oep.plot_dim(adata, 
                    method=method, 
                    color=color, 
                    palette=color_use,
                    show=False, 
                    figsize=(5, 5), 
                    title = color)
    plot_file = os.path.join(outdir, f'groupby-{color}_contrast_plot')
    save_figure(plot_file)