import clscanpy as csc
import anndata as ad
import numpy as np
import os
from pathlib import Path
from typing import List, Union
from matplotlib.gridspec import GridSpec
from clscanpy.tools.utils import check_mkdir
from clscanpy.tools.plotting.summ_faceted import summ_faceted_plot, _plot_single_level
from clscanpy.tools.plotting.summ_contrast import summ_contrast_plot
from clscanpy.log import log_function_call
import random
import logging

LOGGER = logging.getLogger(__name__)


class Faceted_plot:
    def __init__(
        self,
        adata,
        groupby: Union[str, list] = None,
        facetby: Union[str, list] = None,
        dosummary: str = "TRUE",
        method: str = "umap",
        palette: str = None,
        outdir: str = None,
        sampleid: Union[str, list] = None,
        group: Union[str, list] = None,
        clusters: Union[str, list] = None,
        new_celltype: Union[str, list] = None,
        predicate: Union[str, list] = None,
        groupby_levels: str = None,
        metadata: str = None,
    ):
        self.adata = adata
        self.groupby = groupby
        self.facetby = facetby
        self.dosummary = dosummary
        self.method = method
        self.palette = palette
        self.outdir = outdir
        self.sampleid = sampleid
        self.group = group
        self.clusters = clusters
        self.new_celltype = new_celltype
        self.predicate = predicate
        self.groupby_levels = groupby_levels
        self.metadata = metadata

    def run(self):
        LOGGER.info(f"Start faceted_plot ...")
        random.seed(2025)
        self.adata = csc.loadH5AD(
            self.adata,
            sampleid=self.sampleid,
            group=self.group,
            clusters=self.clusters,
            new_celltype=self.new_celltype,
            predicate=self.predicate,
            metadata=self.metadata,
            palette=self.palette,
        )
        LOGGER.info(self.adata.obs.group)

        check_mkdir(self.outdir)
        for groupfactor in self.groupby:
            for facet_col in self.facetby:
                summ_faceted_plot(
                    self.adata,
                    groupby=facet_col,
                    color=groupfactor,
                    method=self.method,
                    outdir=self.outdir,
                    palette=self.palette,
                )
                summ_contrast_plot(
                    self.adata,
                    method=self.method,
                    color=facet_col,
                    outdir=self.outdir,
                    palette=self.palette,
                )

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        LOGGER.info("faceted plot done")


@log_function_call
def faceted_plot(args):
    with Faceted_plot(
        adata=args.input,
        groupby=args.groupby.split(","),
        facetby=args.facetby.split(","),
        method=args.method,
        palette=args.palette,
        outdir=args.outdir,
        sampleid=args.sampleid,
        group=args.group,
        clusters=args.clusters,
        new_celltype=args.new_celltype,
        predicate=args.predicate,
        groupby_levels=args.groupby_levels,
        metadata=args.metadata,
    ) as runner:
        runner.run()


def get_opts_faceted_plot(parser, sub_program=True):
    parser.add_argument(
        "-i", "--input", type=str, default=None, help="Input file", required=True
    )
    parser.add_argument(
        "-o", "--outdir", type=str, default=None, help="Output directory", required=True
    )
    parser.add_argument(
        "--groupby", type=str, default=None, help="分组变量", required=True
    )
    parser.add_argument(
        "--facetby", type=str, default=None, help="分面分组变量", required=True
    )
    parser.add_argument("--palette", type=str, default=None, help="颜色调色板")
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        default="umap",
        help="Dimension reduction method (umap or tsne)",
    )
    parser.add_argument(
        "-l", "--groupby_levels", type=str, default=None, help="Groupby levels"
    )
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
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        help="Metadata file to load, if not provided, will use adata.obs",
    )
    return parser


if __name__ == "__main__":
    unittest.main()
