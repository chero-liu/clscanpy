import matplotlib.pyplot as plt
import scanpy as sc
import matplotlib.patches as mpatches
from oescanpy.tools.utils import save_figure
from oescanpy.tools.color.utils import get_color_order
import numpy as np


def plot_marker_heatmap(
    adata,
    markers,
    groupby="clusters",
    save_path=None,
    save_formats=("png", "pdf"),
    show=False,
    palette=None,
):
    all_genes = set(gene for genes in markers.values() for gene in genes)
    bbox_to_anchor_y = ((len(all_genes)+1) // 2) / 10

    colors_dict = get_color_order(adata, groupby, palette)

    adata.uns[f"{groupby}_colors"] = np.array(list(colors_dict.values()), dtype=object)

    fig = sc.pl.heatmap(
        adata,
        markers,
        groupby=groupby,
        swap_axes=True,
        show_gene_labels=True,
        show=False,
        dendrogram=False,
        standard_scale="var",
    )
    fig["groupby_ax"].set_xticklabels("")
    ax = plt.gca()

    legend_patches = [
        mpatches.Patch(color=color, label=label) for label, color in colors_dict.items()
    ]

    ax.legend(
        handles=legend_patches,
        bbox_to_anchor=(len(markers),bbox_to_anchor_y),
        loc='best',
        ncol=(len(markers) + 9) // 10,
        borderaxespad=0.0,
        title=groupby,
        frameon=False,
    )
    if save_path != None:
        save_figure(save_path)
