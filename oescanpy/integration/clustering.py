from oescanpy.tools.utils import s_common, Step
from oescanpy.tools.cluster import Cluster
from oescanpy.tools.rankgenes import RankGenes
from oescanpy.tools.gcr import GCR
import unittest
import matplotlib
from oescanpy.tools.utils import check_mkdir
from oescanpy.tools.markergenes.utils import write_diff_genes

import subprocess
from oescanpy.config import ROOT_PATH

matplotlib.use("Agg")
import random

random.seed(2025)
from oescanpy.log import log_function_call

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
        self.min_dist = args.min_dist
        self.resolution = args.resolution
        self.target_sum = args.target_sum
        self.sampleid = args.sampleid
        self.group = args.group
        self.clusters = args.clusters
        self.new_celltype = args.new_celltype
        self.predicate = args.predicate
        self.metadata = args.metadata
        self.method = args.method
        self.palette = args.palette
        self.use_raw = args.use_raw
        self.pts = args.pts
        self.corr_method = args.corr_method
        self.tie_correct = args.tie_correct
        self.rankby_abs = args.rankby_abs
        self.key_added = args.key_added
        self.minpct = args.minpct
        self.pvals = args.pvals
        self.refgenome = args.refgenome
        self.top_n = args.top_n
        self.n_comps = args.n_comps
        self.hvg_add = args.hvg_add
        self.hvg_remove = args.hvg_remove
        self.hvg_replace = args.hvg_replace

    @log_function_call
    def cluster(self):
        create_obj = Cluster(
            adata=self.input,
            hvg_num=self.hvg_num,
            n_pcs=self.n_pcs,
            batch_method=self.batch_method,
            batch_name=self.batch_name,
            clust_method=self.clust_method,
            min_dist=self.min_dist,
            resolution=self.resolution,
            target_sum=self.target_sum,
            outdir=self.outdir,
            sampleid=self.sampleid,
            group=self.group,
            clusters=self.clusters,
            new_celltype=self.new_celltype,
            predicate=self.predicate,
            metadata=self.metadata,
            n_comps=self.n_comps,
            hvg_add=self.hvg_add,
            hvg_remove=self.hvg_remove,
            hvg_replace=self.hvg_replace,
        )
        adata = create_obj.run()

        return adata

    @log_function_call
    def gcr(self, adata):
        gcr_obj = GCR(
            input=adata,
            outdir=f"{self.outdir}",
            groupby="clusters",
            sampleid=None,
            group=None,
            clusters=None,
            new_celltype=None,
            predicate=None,
            metadata=None,
            clust_method=self.clust_method,
            method=self.method,
            palette=self.palette,
            cloudcfg=False,
            groupby_levels=None,
            resolution=self.resolution,
        )
        gcr_obj.run()

    @log_function_call
    def rankgenes(self, adata):
        rankgenes_obj = RankGenes(
            adata=adata,
            outdir=f"{self.outdir}/Markers",
            groupby="clusters",
            use_raw=self.use_raw,
            pts=self.pts,
            method="wilcoxon",
            corr_method=self.corr_method,
            tie_correct=self.tie_correct,
            rankby_abs=self.rankby_abs,
            key_added=self.key_added,
            sampleid=None,
            group=None,
            clusters=None,
            new_celltype=None,
            predicate=None,
            save_h5ad=False,
        )
        adata = rankgenes_obj.run_RankGeneDefault()
        write_diff_genes(
            adata,
            outdir=f"{self.outdir}/Markers",
            minpct=self.minpct,
            pvals=self.pvals,
            refgenome=self.refgenome,
            groupby="clusters",
            top_n=self.top_n,
        )
        return adata

    def run(self):
        LOGGER.info("Start Clustering pipeline ...")
        check_mkdir(f"{self.outdir}/Markers")

        adata = self.cluster()

        # self.gcr(adata)

        adata = self.rankgenes(adata)

        adata.write_h5ad(f"{self.outdir}/adata.h5ad")

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
        type=str,
        default=10,
        help="Number of principal components (default: 10)",
    )
    parser.add_argument(
        "--n_comps",
        type=int,
        default=50,
        help="(default: 50)",
    )
    parser.add_argument(
        "--batch_method",
        type=str,
        default="No",
        help="Batch correction method:harmony/combat/scanorama/No",
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
        "--min_dist",
        default="default",
        help="Minimum distance for clustering (default)",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=0.4,
        help="Resolution for clustering (default: 0.4)",
    )
    parser.add_argument(
        "--target_sum", type=float, default=10000, help="Target sum for normalization"
    )
    parser.add_argument(
        "--method",
        type=str,
        default="umap",
        help="Visualization method (default: 'umap')",
    )
    parser.add_argument(
        "--use_raw",
        type=bool,
        default=False,
        help="Whether to use raw data for differential expression analysis (default: False)",
    )
    parser.add_argument(
        "--pts",
        type=bool,
        default=True,
        help="Whether to calculate percentage of cells expressing each gene (default: True)",
    )
    parser.add_argument(
        "--corr_method",
        type=str,
        choices=["benjamini-hochberg", "bonferroni"],
        default="benjamini-hochberg",
        help="Correction method for p-values (default: 'benjamini-hochberg')",
    )
    parser.add_argument(
        "--tie_correct",
        type=bool,
        default=False,
        help="Whether to apply tie correction in the ranking algorithm (default: False)",
    )
    parser.add_argument(
        "--rankby_abs",
        type=bool,
        default=False,
        help="Whether to rank genes by absolute values (default: False)",
    )
    parser.add_argument(
        "--key_added",
        type=str,
        default=None,
        help="Key added to adata.uns for storing results",
    )
    parser.add_argument(
        "--minpct",
        type=float,
        default=0.25,
        help="Minimum percentage of cells expressing a gene (default: 0.25)",
    )
    parser.add_argument(
        "--pvals",
        type=float,
        default=0.05,
        help="P-value threshold for filtering genes (default: 0.05)",
    )
    parser.add_argument(
        "--refgenome",
        type=str,
        default=None,
        help="Reference genome for annotation (default: None)",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=5,
        help=" (default: 5)",
    )
    parser.add_argument(
        "--hvg_add",
        default=None,
        help="(default: None)",
    )
    parser.add_argument(
        "--hvg_remove",
        default=None,
        help="(default: None)",
    )
    parser.add_argument(
        "--hvg_replace",
        default=None,
        help="(default: None)",
    )

    # Common parameters
    if sub_program:
        parser = s_common(parser)

    return parser


if __name__ == "__main__":
    unittest.main()
