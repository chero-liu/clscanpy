import os
from clscanpy.tools.utils import s_common, Step, check_mkdir, check_adata
import unittest
import clscanpy as csc
from clscanpy.tools.preprocessing.basicpreprocessing import filter_genes_for_pyscenic
from clscanpy.script.pyscenic.utils import Analysis, Visualize
from clscanpy.config import ROOT_PATH
import subprocess
from clscanpy.log import log_function_call

import logging

LOGGER = logging.getLogger(__name__)


class pySCENIC(Step):
    def __init__(
        self,
        args,
    ):
        Step.__init__(self, args)
        self.input = args.input
        self.outdir = args.outdir
        self.sampleid = args.sampleid
        self.group = args.group
        self.clusters = args.clusters
        self.new_celltype = args.new_celltype
        self.predicate = args.predicate
        self.downsample = args.downsample
        self.groupby = args.groupby
        self.groupby_levels = args.groupby_levels
        self.palette = args.palette
        self.database = args.database
        self.species = args.species
        self.method = args.method
        self.tfs = args.tfs
        self.motifs_tbl = args.motifs_tbl
        self.rank_threshold = args.rank_threshold
        self.auc_threshold = args.auc_threshold
        self.nes_threshold = args.nes_threshold
        self.min_orthologous_identity = args.min_orthologous_identity
        self.max_similarity_fdr = args.max_similarity_fdr
        self.chunk_size = args.chunk_size
        self.thresholds = args.thresholds
        self.top_n_targets = args.top_n_targets
        self.min_genes = args.min_genes
        self.all_modules = args.all_modules
        self.threshold = args.threshold
        self.num_workers = args.thread
        self.topGenes = args.topGenes
        self.extended = args.extended
        self.nclust = args.nclust
        self.utils_path = args.utils_path
        self.filter_genes = args.filter_genes

        if self.species == "human":
            self.genes_vs_motifs = [
                "500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
                "10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
            ]

        elif self.species == "mouse":
            self.genes_vs_motifs = [
                "500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather",
                "500bp-upstream-10species.mc9nr.genes_vs_motifs.rankings.feather",
            ]

        elif self.species == "Malus_domestica_gene_MDP":
            self.genes_vs_motifs = [
                "1000bp_up_1000bp_down_Malus_domestica.regions_vs_motifs.rankings.feather",
            ]
        elif self.species == "Aegilops_tauschii_gene_AT":
            self.genes_vs_motifs = [
                "At.genes_vs_motifs.rankings.feather",
            ]

        self.db_paths = [
            os.path.join(self.database, self.species, fname)
            for fname in self.genes_vs_motifs
        ]

    @log_function_call
    def datapreparation(self):
        adata = csc.loadH5AD(
            input=self.input,
            sampleid=self.sampleid,
            group=self.group,
            clusters=self.clusters,
            new_celltype=self.new_celltype,
            predicate=self.predicate,
            downsample=self.downsample,
            groupby=self.groupby,
            groupby_levels=self.groupby_levels,
            palette=self.palette,
        )
        if self.filter_genes.lower() == "true":
            fadata = filter_genes_for_pyscenic(
                adata,
                self.db_paths,
            )
        else:
            fadata = adata
        return fadata

    @log_function_call
    def analysis(self, input, outdir):
        analysis_obj = Analysis(
            input=input,
            outdir=outdir,
            species=self.species,
            method=self.method,
            database=self.database,
            tfs=self.tfs,
            genes_vs_motifs=self.db_paths,
            motifs_tbl=self.motifs_tbl,
            rank_threshold=self.rank_threshold,
            auc_threshold=self.auc_threshold,
            nes_threshold=self.nes_threshold,
            min_orthologous_identity=self.min_orthologous_identity,
            max_similarity_fdr=self.max_similarity_fdr,
            chunk_size=self.chunk_size,
            thresholds=self.thresholds,
            top_n_targets=self.top_n_targets,
            min_genes=self.min_genes,
            all_modules=self.all_modules,
            num_workers=self.num_workers,
        )
        analysis_obj.run()

    @log_function_call
    def visualize(self, input):
        visualize_obj = Visualize(
            input=os.path.join(input, "sce_SCENIC.loom"),
            outdir=os.path.join(self.outdir, "result"),
            species=self.species,
            rds_filepath=self.input,
            groupby=self.groupby,
            threshold=self.threshold,
            regulons_path=os.path.join(self.outdir, "result", "regulons.xls"),
            topGenes=self.topGenes,
            extended=self.extended,
            nclust=self.nclust,
            utils_path=self.utils_path,
        )
        visualize_obj.run()

    def run(self):
        LOGGER.info("Start pySCENIC pipeline ...")
        check_mkdir(f"{self.outdir}/data/")
        check_mkdir(f"{self.outdir}/script/")
        check_mkdir(f"{self.outdir}/result/")

        adata = self.datapreparation()

        adata = check_adata(adata)
        if not hasattr(adata.var, "Gene"):
            adata.var["Gene"] = adata.var_names
        if not hasattr(adata.obs, "CellID"):
            adata.obs["CellID"] = adata.obs_names
        adata.write_loom(f"{self.outdir}/data/for_scenic.loom", write_obsm_varm=True)

        self.analysis(f"{self.outdir}/data/for_scenic.loom", f"{self.outdir}/data")

        subprocess.run(
            [
                "mv",
                f"{self.outdir}/data/regulons.xls",
                f"{self.outdir}/result/regulons.xls",
            ],
            check=True,
        )
        subprocess.run(
            [
                "mv",
                f"{self.outdir}/data/TF_TargetGenes.xls",
                f"{self.outdir}/result/TF_TargetGenes.xls",
            ],
            check=True,
        )

        self.visualize(f"{self.outdir}/data/")


@log_function_call
def pyscenic(args):
    with pySCENIC(args) as runner:
        runner.run()


def get_opts_pyscenic(parser, sub_program=True):
    # Data Preparation arguments
    parser.add_argument("--input", required=True, help="")

    # Analysis arguments
    parser.add_argument(
        "--method",
        default="grnboost2",
        choices=["grnboost2", "genie3"],
        help="GRN inference method",
    )
    parser.add_argument(
        "--database",
        default="/gpfs/oe-database/pySCENIC/",
        help="Database directory path",
    )
    parser.add_argument("--tfs", default="allTFs.txt", help="TF list filename")
    parser.add_argument(
        "--motifs_tbl",
        default="motifs-v10nr_clust-nr-m0.001-o0.0.tbl",
        help="Motif table filename",
    )
    parser.add_argument(
        "--rank_threshold",
        type=int,
        default=5000,
        help="Gene rank threshold (default: 5000)",
    )
    parser.add_argument(
        "--auc_threshold",
        type=float,
        default=0.05,
        help="AUC calculation threshold (default: 0.05)",
    )
    parser.add_argument(
        "--nes_threshold",
        type=float,
        default=3.0,
        help="NES significance threshold (default: 3.0)",
    )
    parser.add_argument(
        "--min_orthologous_identity",
        type=float,
        default=0.0,
        help="Minimum orthologous identity (default: 0.0)",
    )
    parser.add_argument(
        "--max_similarity_fdr",
        type=float,
        default=0.001,
        help="Maximum similarity FDR (default: 0.001)",
    )
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=200,
        help="Chunk size for module partitioning (default: 200)",
    )
    parser.add_argument(
        "--thresholds",
        type=str,
        default="0.75 0.90",
        help="Correlation coefficient thresholds (default: 0.75 0.90)",
    )
    parser.add_argument(
        "--top_n_targets",
        type=int,
        default=50,
        help="Top N target counts (default: 50)",
    )
    parser.add_argument(
        "--min_genes",
        type=int,
        default=20,
        help="Minimum number of genes (default: 20)",
    )
    parser.add_argument(
        "--all_modules",
        default="False",
        help="Include both positive and negative regulatory modules (default: False)",
    )

    # Visualization arguments
    parser.add_argument(
        "--groupby", default="clusters", help="Grouping column for visualization"
    )
    parser.add_argument("--groupby_levels", default=None, help="")
    parser.add_argument("--threshold", type=float, default=0.0, help="AUC threshold")
    parser.add_argument(
        "--topGenes", type=int, default=3, help="Number of top genes to show"
    )
    parser.add_argument("--extended", default="FALSE", help="Extended visualization")
    parser.add_argument("--nclust", type=int, default=4, help="Number of clusters")
    parser.add_argument(
        "--utils_path",
        default=f"{ROOT_PATH}/script/pyscenic/utils.r",
        help="R utils path for visualization",
    )
    parser.add_argument(
        "--filter_genes",
        default="True",
        help="Whether to filter genes not in the database (default: True)",
    )

    if sub_program:
        parser = s_common(parser)
    return parser


if __name__ == "__main__":
    unittest.main()
