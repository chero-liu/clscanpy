# get clustering result
import os
import pandas as pd
import clscanpy as csc
import unittest
from clscanpy.tools.utils import check_mkdir
from clscanpy.log import log_function_call

import logging
from typing import Union, List

from clscanpy.tools.markergenes.utils import (
    vis_markers,
)


LOGGER = logging.getLogger(__name__)


class Vismarker:
    def __init__(
        self,
        input: str = None,
        outdir: str = None,
        groupby: str = "clusters",
        groupby_levels: str = None,
        sampleid: Union[str, list] = None,
        group: Union[str, list] = None,
        clusters: Union[str, list] = None,
        new_celltype: Union[str, list] = None,
        predicate: Union[str, list] = None,
        method: str = "umap",
        palette: str = None,
        top: int = 3,
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
        self.method = method
        self.palette = palette
        self.top = top

    def run(self):
        LOGGER.info("Start vismarker ...")
        check_mkdir(self.outdir)
        adata = csc.loadH5AD(
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

        vis_markers(
            adata,
            top=self.top,
            groupby=self.groupby,
            outdir=self.outdir,
            method=self.method,
            palette=self.palette,
        )

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        LOGGER.info("vismarker done")


@log_function_call
def vismarker(args):
    with Vismarker(
        input=args.input,
        outdir=args.outdir,
        groupby=args.groupby,
        groupby_levels=args.groupby_levels,
        sampleid=args.sampleid,
        group=args.group,
        clusters=args.clusters,
        new_celltype=args.new_celltype,
        predicate=args.predicate,
        method=args.method,
        palette=args.palette,
        top=args.top,
    ) as runner:
        runner.run()


def get_opts_vismarker(parser, sub_program=True):
    parser.add_argument(
        "-i", "--input", type=str, default=None, help="Input file", required=True
    )
    parser.add_argument(
        "-o", "--outdir", type=str, default=None, help="Output directory", required=True
    )
    parser.add_argument(
        "-t", "--top", default=3, help="Top n marker genes to visualize"
    )
    parser.add_argument(
        "-c", "--groupby", type=str, default="clusters", help="Color by column"
    )
    parser.add_argument(
        "-l", "--groupby_levels", type=str, default=None, help="Groupby levels"
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        default="umap",
        help="Dimension reduction method (umap or tsne)",
    )
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
    return parser


if __name__ == "__main__":
    unittest.main()
