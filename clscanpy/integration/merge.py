from clscanpy.tools.utils import s_common, Step
from clscanpy.tools.create import Create
from clscanpy.tools.cluster import Cluster
from clscanpy.tools.rankgenes import RankGenes
from clscanpy.tools.autoanno import Autoanno
import unittest
import matplotlib

matplotlib.use("Agg")
import random

random.seed(2025)


class Merge(Step):
    def __init__(self, args):
        Step.__init__(self, args)
        # create related parameters
        self.info = args.info
        self.refgenome = args.refgenome
        self.mtfilter = args.mtfilter
        self.hbfilter = args.hbfilter
        self.mingene = args.mingene
        self.maxgene = args.maxgene
        self.minumi = args.minumi
        self.maxumi = args.maxumi
        self.minlgGenePerUmi = args.minlgGenePerUmi
        self.maxlgGenePerUmi = args.maxlgGenePerUmi
        self.mincell = args.mincell
        self.rmdouble = args.rmdouble
        self.rmdoublet_method = args.rmdoublet_method

        # cluster related parameters
        self.hvg_num = args.hvg_num
        self.n_pcs = args.n_pcs
        self.batch_method = args.batch_method
        self.batch_name = args.batch_name
        self.clust_method = args.clust_method
        self.min_dist = args.min_dist
        self.resolution = args.resolution
        self.target_sum = args.target_sum

        # rankgenes related parameters
        self.method = args.method
        self.logfc = args.logfc
        self.minpct = args.minpct
        self.groupby = args.groupby
        self.use_raw = args.use_raw
        self.pts = args.pts
        self.top = args.top

        # autoanno related parameters
        self.model_path = args.model_path
        self.mode = args.mode
        self.cluster_key = args.cluster_key
        self.tissue = args.tissue

        # common parameters
        self.groupby = "clusters"

    def create(self):
        create_obj = Create(
            info=self.info,
            outdir=f"{self.outdir}/QC/",
            species=self.species,
            refgenome=self.refgenome,
            mtfilter=self.mtfilter,
            hbfilter=self.hbfilter,
            mingene=self.mingene,
            maxgene=self.maxgene,
            minumi=self.minumi,
            maxumi=self.maxumi,
            minlgGenePerUmi=self.minlgGenePerUmi,
            maxlgGenePerUmi=self.maxlgGenePerUmi,
            mincell=self.mincell,
            rmdouble=self.rmdouble,
            rmdoublet_method=self.rmdoublet_method,
        )
        adata = create_obj.run()

        return adata

    def cluster(self, adata):
        cluster_obj = Cluster(
            adata=adata,
            hvg_num=self.hvg_num,
            n_pcs=self.n_pcs,
            batch_method=self.batch_method,
            batch_name=self.batch_name,
            clust_method=self.clust_method,
            min_dist=self.min_dist,
            resolution=self.resolution,
            target_sum=self.target_sum,
            outdir=f"{self.outdir}/Clustering/",
        )
        adata = cluster_obj.run()
        cluster_obj.get_result(save=False)
        return adata

    def rankgenes(self, adata):
        rank_obj = RankGenes(
            adata=adata,
            outdir=f"{self.outdir}/Markers/",
            method=self.method,
            logfc=self.logfc,
            minpct=self.minpct,
            groupby=self.groupby,
            use_raw=self.use_raw,
            pts=self.pts,
            top=self.top,
        )
        adata = rank_obj.run_RankGeneDefault()
        rank_obj.get_result()
        return adata

    def autoanno(self):
        autoanno_obj = Autoanno(
            input=f"{self.outdir}/PRO_diff.h5ad",
            outdir=f"{self.outdir}/Autoanno/",
            model_path=self.model_path,
            mode=self.mode,
            cluster_key=self.cluster_key,
            species=self.species,
            tissue=self.tissue,
        )
        autoanno_obj.run()

    def run(self):
        print("Start Merge pipeline ...")

        adata = self.create()

        adata = self.cluster(adata)

        adata = self.rankgenes(adata)

        adata.write(f"{self.outdir}/PRO_diff.h5ad")

        self.autoanno()


def merge(args):
    with Merge(args) as runner:
        runner.run()


def get_opts_merge(parser, sub_program=True):
    # Create related parameters
    parser.add_argument(
        "--refgenome", type=str, default=None, help="reference genome name"
    )
    parser.add_argument(
        "--mtfilter",
        type=str,
        default="default",
        help="if default, automatically choose a mt threshold for all samples. Recommend",
    )
    parser.add_argument(
        "--hbfilter", type=float, default=0.05, help="hb gene filter threshold"
    )
    parser.add_argument(
        "--mingene", type=float, default=200, help="minimum gene content cutoff"
    )
    parser.add_argument(
        "--maxgene", type=float, default=1, help="maximum gene content cutoff"
    )
    parser.add_argument(
        "--minumi", type=float, default=1000, help="minimum umi content cutoff"
    )
    parser.add_argument(
        "--maxumi", type=float, default=1, help="maximum umi content cutoff"
    )
    parser.add_argument(
        "--minlgGenePerUmi",
        type=float,
        default=0.7,
        help="minimum log10(genes per UMI) cutoff",
    )
    parser.add_argument(
        "--maxlgGenePerUmi",
        type=float,
        default=None,
        help="maximum log10(genes per UMI) cutoff",
    )
    parser.add_argument(
        "--mincell",
        type=int,
        default=0,
        help="filter genes if exists in less than an exact number of cells",
    )
    parser.add_argument(
        "--rmdouble", type=str, default="True", help="whether remove doublets"
    )
    parser.add_argument(
        "--rmdoublet_method",
        type=str,
        default="doubletdetection",
        help="method to remove doublets",
    )

    # Cluster related parameters
    parser.add_argument("--hvg_num", type=int, default=2000, help="hvg genes number")
    parser.add_argument(
        "--n_pcs",
        type=int,
        default=10,
        help="calculate principle components automatically or use an exact pc number",
    )
    parser.add_argument(
        "--batch_method",
        type=str,
        default="No",
        help="whether remove batch, and which method used ['harmony', 'combat', 'scanorama', 'No']",
    )
    parser.add_argument(
        "--batch_name",
        type=str,
        default="batchid",
        help="batch column name in adata.obs",
    )
    parser.add_argument(
        "--clust_method",
        type=str,
        default="leiden",
        help="cluster method, louvain or leiden",
    )
    parser.add_argument(
        "--min_dist", type=str, default="default", help="a umap argument"
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=0.4,
        help="a parameter value controlling the coarseness of the clustering",
    )
    parser.add_argument(
        "--target_sum", type=float, default=10000, help="target sum for normalization"
    )

    # RankGenes related parameters
    parser.add_argument(
        "--method", type=str, default="wilcoxon", help="rank genes method."
    )
    parser.add_argument(
        "--logfc", type=float, default=0.25, help="logfoldchanges cutoff"
    )
    parser.add_argument("--minpct", type=float, default=0.1, help="min pct cutoff")
    parser.add_argument(
        "--groupby",
        type=str,
        default="clusters",
        help="group by column for marker gene analysis",
    )
    parser.add_argument(
        "--use_raw",
        type=bool,
        default=False,
        help="use raw data for marker gene analysis",
    )
    parser.add_argument(
        "--pts",
        type=bool,
        default=True,
        help="calculate fraction of cells expressing the genes",
    )
    parser.add_argument("--top", type=int, default=3, help="top n marker genes to show")

    # Autoanno related parameters
    parser.add_argument(
        "--model_path",
        type=str,
        default=None,
        help="celltypist model path, if None, will use default model",
    )
    parser.add_argument(
        "--mode",
        type=str,
        default="best_match",
        help="celltypist mode, best_match or majority_voting",
    )
    parser.add_argument(
        "--cluster_key",
        type=str,
        default="clusters",
        help="celltypist cluster key, default is clusters",
    )
    parser.add_argument(
        "--tissue",
        type=str,
        default="default",
        help="celltypist tissue, default is default, human model only support default",
    )
    # Common parameters
    if sub_program:
        parser = s_common(parser)
        parser.add_argument(
            "--info",
            help="a project description config file, include path,rspname,sampleid,datatype,group,batchid",
        )
    return parser


if __name__ == "__main__":
    unittest.main()
