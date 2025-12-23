from clscanpy.tools.utils import s_common, Step
from clscanpy.tools.cluster import Cluster
from clscanpy.tools.rankgenes import RankGenes
from clscanpy.tools.gcr import GCR
from clscanpy.tools.enrich import Enrich
from clscanpy.tools.vismarker import Vismarker
import unittest
import matplotlib
from clscanpy.tools.utils import check_mkdir
from clscanpy.tools.markergenes.utils import write_diff_genes

import subprocess
from clscanpy.config import ROOT_PATH

matplotlib.use("Agg")
import random

random.seed(2025)
from clscanpy.log import log_function_call

import logging

LOGGER = logging.getLogger(__name__)


class Clustering(Step):
    def __init__(self, args):
        Step.__init__(self, args)
        self.input = args.input
        self.outdir = args.outdir
        self.hvg_num = args.hvg_num
        self.n_pcs = args.n_pcs
        self.batch_method = args.batch_method
        self.batch_name = args.batch_name
        self.clust_method = args.clust_method
        self.resolution = args.resolution
        self.sampleid = args.sampleid
        self.group = args.group
        self.clusters = args.clusters
        self.new_celltype = args.new_celltype
        self.predicate = args.predicate
        self.palette = args.palette
        self.refgenome = args.refgenome

    @log_function_call
    def cluster(self):
        create_obj = Cluster(
            adata=self.input,
            hvg_num=self.hvg_num,
            n_pcs=self.n_pcs,
            batch_method=self.batch_method,
            batch_name=self.batch_name,
            clust_method=self.clust_method,
            resolution=self.resolution,
            outdir=f"{self.outdir}/Cluster",
            sampleid=self.sampleid,
            group=self.group,
            clusters=self.clusters,
            new_celltype=self.new_celltype,
            predicate=self.predicate,
        )
        adata = create_obj.run()

        return adata

    @log_function_call
    def gcr(self, adata):
        gcr_obj = GCR(
            input=adata,
            outdir=f"{self.outdir}/Cluster",
            groupby="clusters",
            clust_method=self.clust_method,
            palette=self.palette,
            resolution=self.resolution,
        )
        gcr_obj.run()

    @log_function_call
    def rankgenes(self, adata):
        rankgenes_obj = RankGenes(
            adata=adata,
            groupby="clusters",
        )
        adata = rankgenes_obj.run_RankGeneDefault()
        write_diff_genes(
            adata,
            outdir=f"{self.outdir}/Marker/tables",
            refgenome=f"{self.refgenome}/annotation/gene_annotation.xls",
            groupby="clusters",
        )
        return adata

    @log_function_call
    def vismarker(self, adata):
        vismarker_obj = Vismarker(
            input=adata,
            outdir=f"{self.outdir}/Marker/figures/",
        )
        vismarker_obj.run()

    @log_function_call
    def enrich(self):
        enrich_obj = Enrich(
            input=f"{self.outdir}/Marker/tables/all_markers_for_each_cluster_anno.xls",
            outdir=f"{self.outdir}/Enrich",
            bp=f"{self.refgenome}/annotation/gmt/mouse_2024_gene_go_bp.gmt",
            cc=f"{self.refgenome}/annotation/gmt/mouse_2024_gene_go_cc.gmt",
            mf=f"{self.refgenome}/annotation/gmt/mouse_2024_gene_go_mf.gmt",
            kegg=f"{self.refgenome}/annotation/gmt/mouse_2024_kegg.gmt",
            skip_report=False,
        )
        enrich_obj.run()

    def run(self):
        LOGGER.info("Start Clustering pipeline ...")
        check_mkdir(f"{self.outdir}")

        adata = self.cluster()

        self.gcr(adata)

        adata = self.rankgenes(adata)

        self.vismarker(adata)

        adata.write_h5ad(f"{self.outdir}/PRO_diff.h5ad")

        self.enrich()


@log_function_call
def clustering(args):
    with Clustering(args) as runner:
        runner.run()


def get_opts_clustering(parser, sub_program=True):
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input file (AnnData object in .h5ad format)",
    )
    parser.add_argument(
        "--hvg_num",
        type=int,
        default=2000,
        help="Number of highly variable genes (default: 2000)",
    )
    parser.add_argument(
        "--n_pcs",
        type=int,
        default=10,
        help="Number of principal components (default: 10)",
    )
    parser.add_argument(
        "--batch_method",
        type=str,
        default="No",
        help="Batch correction method: harmony/combat/scanorama/No",
    )
    parser.add_argument(
        "--batch_name",
        type=str,
        default="batchid",
        help="Batch column name in metadata",
    )
    parser.add_argument(
        "--clust_method",
        type=str,
        default="leiden",
        help="Clustering method (default: 'leiden')",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=1.2,
        help="Resolution for clustering (default: 1.2)",
    )
    parser.add_argument(
        "--refgenome",
        type=str,
        default=None,
        help="Reference genome for annotation (default: None)",
    )

    if sub_program:
        parser = s_common(parser)

    return parser


if __name__ == "__main__":
    unittest.main()
