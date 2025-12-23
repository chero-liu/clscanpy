# get clusters result
import os
import pandas as pd
import clscanpy as oep
import unittest
from clscanpy.tools.utils import (
    check_mkdir,
    add_counts_to_data,
    save_figure,
    get_group_order,
)
from clscanpy.tools.color.utils import get_color_order
from clscanpy.log import log_function_call
from typing import Union, List
import logging

LOGGER = logging.getLogger(__name__)
import json
import numpy as np
import scipy


def save_dimension_reduction_results(
    adata,
    outdir,
    method="umap",
    color_by="clusters",
    palette=None,
    clust_method="leiden",
    resolution=None,
):
    """
    Save dimension reduction results (e.g., UMAP, t-SNE) and their visualizations.

    Parameters:
        adata: AnnData object containing the data.
        outdir: Output directory to save results.
        method: Dimension reduction method ('umap' or 'tsne').
        color_by: Column in `adata.obs` to color the plots by.
        clust_method: Key in `adata.uns` to extract resolution for file naming.
    """
    # Validate method
    if method not in ["umap", "tsne"]:
        raise ValueError("Method must be 'umap' or 'tsne'")

    # Create output directory
    check_mkdir(os.path.join(outdir, "figures"))
    check_mkdir(os.path.join(outdir, "tables"))

    # Extract dimension reduction coordinates
    coord_key = f"X_{method}"
    if coord_key not in adata.obsm:
        raise KeyError(f"{coord_key} not found in adata.obsm")

    tmp = pd.DataFrame(adata.obsm[coord_key]).rename(
        columns={0: f"{method.upper()}_1", 1: f"{method.upper()}_2"}
    )
    tmp.insert(0, "Barcode", adata.obs["rawbc"].values)
    coord_file = os.path.join(
        os.path.join(outdir, "tables"), f"{method}_Dimension_Reduction_coordination.csv"
    )
    tmp.to_csv(coord_file, index=False)

    # add counts to the 'clusters' column
    # adata.obs = add_counts_to_data(adata.obs, color_by)

    oep.plot_dim(
        adata,
        method=method,
        color=color_by,
        palette=get_color_order(adata, color_by, palette),
        show=False,
        figsize=(5, 5),
        title=color_by,
    )

    if resolution == None:
        try:
            resolution = json.loads(adata.uns["preprocess_para"])["resolution"]
        except KeyError:
            resolution = "None"

    plot_file = os.path.join(
        os.path.join(outdir, "figures"),
        f"{method}_res{resolution}",
    )
    save_figure(plot_file)

    oep.plot_dim(
        adata,
        method=method,
        color=color_by,
        palette=get_color_order(adata, color_by, palette),
        show=False,
        figsize=(5, 5),
        legend_loc="on data",
    )
    plot_file = os.path.join(
        os.path.join(outdir, "figures"),
        f"{method}_res{resolution}_on_data",
    )
    save_figure(plot_file)


def save_clusters_results(
    adata,
    outdir,
    method="umap",
    color_by="clusters",
    resolution=None,
    clust_method="leiden",
):
    """
    Save clusters results to a CSV file.

    Parameters:
        adata: AnnData object containing the data.
        outdir: Output directory to save results.
    """
    check_mkdir(outdir)
    clusters_file = os.path.join(os.path.join(outdir, "tables"), f"metadata.csv")
    data = adata.obs[["rawbc", "sampleid", "group", color_by]]
    data = data.rename(
        columns={
            "rawbc": "Barcode",
        }
    )
    data.to_csv(clusters_file, index=False)

    LOGGER.info(f"Saved clusters results to {clusters_file}")

    if resolution == None:
        try:
            resolution = adata.uns[clust_method]["params"]["resolution"]
        except KeyError:
            resolution = "None"


def get_clusters_result(
    adata,
    outdir,
    color_by="clusters",
    method="umap",
    palette=None,
    resolution=None,
    clust_method="leiden",
):
    """
    Main function to save clusters and dimension reduction results.

    Parameters:
        adata: AnnData object containing the data.
        outdir: Base output directory.
    """
    # Save clusters results
    save_clusters_results(
        adata,
        outdir=outdir,
        method=method,
        color_by=color_by,
        resolution=resolution,
        clust_method=clust_method,
    )

    # Save UMAP results
    save_dimension_reduction_results(
        adata,
        outdir,
        method="umap",
        color_by=color_by,
        palette=palette,
        clust_method=clust_method,
        resolution=resolution,
    )
    save_dimension_reduction_results(
        adata,
        outdir,
        method="tsne",
        color_by=color_by,
        palette=palette,
        clust_method=clust_method,
        resolution=resolution,
    )


class GCR:

    def __init__(
        self,
        input,
        outdir: str = None,
        groupby: str = None,
        groupby_levels: str = None,
        sampleid: Union[str, list] = None,
        group: Union[str, list] = None,
        clusters: Union[str, list] = None,
        new_celltype: Union[str, list] = None,
        predicate: Union[str, list] = None,
        clust_method: str = "leiden",
        palette=None,
        resolution=None,
        save_h5ad: bool = False,
    ):
        self.input = input
        self.outdir = outdir
        self.groupby = groupby
        self.groupby_levels = groupby_levels
        self.sampleid = sampleid
        self.group = group
        self.clusters = clusters
        self.new_celltype = new_celltype
        self.predicate = predicate
        self.clust_method = clust_method
        self.palette = palette
        self.resolution = resolution
        self.save_h5ad = save_h5ad

    def run(self):
        LOGGER.info("Start GCR ...")
        adata = oep.loadH5AD(
            self.input,
            sampleid=self.sampleid,
            group=self.group,
            clusters=self.clusters,
            new_celltype=self.new_celltype,
            predicate=self.predicate,
            groupby=self.groupby,
            groupby_levels=self.groupby_levels,
            palette=self.palette,
        )

        get_clusters_result(
            adata,
            self.outdir,
            color_by=self.groupby,
            palette=self.palette,
            clust_method=self.clust_method,
            resolution=self.resolution,
        )

        if self.save_h5ad:
            adata.write_h5ad(os.path.join(self.outdir, "adata.h5ad"))

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        LOGGER.info("gcr done")


@log_function_call
def gcr(args):
    with GCR(
        input=args.input,
        outdir=args.outdir,
        groupby=args.groupby,
        groupby_levels=args.groupby_levels,
        sampleid=args.sampleid,
        group=args.group,
        clusters=args.clusters,
        new_celltype=args.new_celltype,
        predicate=args.predicate,
        clust_method=args.clust_method,
        palette=args.palette,
        resolution=args.resolution,
        save_h5ad=args.save_h5ad,
    ) as runner:
        runner.run()


def get_opts_gcr(parser, sub_program=True):
    parser.add_argument(
        "-i", "--input", type=str, default=None, help="Input file", required=True
    )
    parser.add_argument(
        "-o", "--outdir", type=str, default=None, help="Output directory", required=True
    )
    parser.add_argument(
        "-c", "--groupby", type=str, default="clusters", help="Color by column"
    )
    parser.add_argument(
        "-l", "--groupby_levels", type=str, default=None, help="Groupby levels"
    )
    parser.add_argument(
        "--clust_method", type=str, default="leiden", help="Resolution key"
    )
    parser.add_argument("--resolution", type=str, default=None, help="resolution")
    parser.add_argument("-p", "--palette", type=str, default=None, help="")
    parser.add_argument(
        "--sampleid",
        type=str,
        default=None,
        help="sampleid column name in adata.obs",
    )
    parser.add_argument(
        "--group",
        type=str,
        default=None,
        help="group column name in adata.obs",
    )
    parser.add_argument(
        "--clusters",
        type=str,
        default=None,
        help="clusters column name in adata.obs",
    )
    parser.add_argument(
        "--new_celltype",
        type=str,
        default=None,
        help="new celltype column name in adata.obs",
    )
    parser.add_argument(
        "--predicate",
        type=str,
        default=None,
        help="predicate for filtering adata.obs, e.g. (sampleid in ['STB1', 'STB4']) and ~(clusters in ['1', '5', '6'])",
    )
    parser.add_argument("--save_h5ad", default=True, help="")
    return parser


if __name__ == "__main__":
    unittest.main()
