from anndata import AnnData
from typing import Optional, Literal, Union, List, Tuple
from clscanpy.tools.markergenes.rank_genes_groups import rank_genes_groups
import scanpy as sc
import clscanpy as oep
from clscanpy.tools.utils import check_mkdir
import os
import numpy as np
from clscanpy.tools.markergenes.utils import write_diff_genes
from clscanpy.log import log_function_call
from distutils.util import strtobool

import logging

LOGGER = logging.getLogger(__name__)


class RankGenes:
    """
    Find diff genes
    """

    def __init__(
        self,
        adata: AnnData,
        outdir: str = None,
        groupby: str = None,
        use_raw: Optional[bool] = False,
        pts: bool = True,
        method: Literal[
            "logreg", "t-test", "wilcoxon", "t-test_overestim_var"
        ] = "wilcoxon",
        corr_method: Literal["benjamini-hochberg", "bonferroni"] = "benjamini-hochberg",
        tie_correct: bool = False,
        rankby_abs: bool = False,
        key_added: Optional[str] = None,
        sampleid: Union[str, list] = None,
        group: Union[str, list] = None,
        clusters: Union[str, list] = None,
        new_celltype: Union[str, list] = None,
        predicate: Union[str, list] = None,
        metadata: Union[str] = None,
        save_h5ad: bool = False,
    ):
        self.adata = adata
        self.outdir = outdir
        self.method = method
        self.groupby = groupby
        self.use_raw = use_raw
        self.pts = pts
        self.corr_method = corr_method
        self.tie_correct = tie_correct
        self.rankby_abs = rankby_abs
        self.key_added = key_added
        self.sampleid = sampleid
        self.group = group
        self.clusters = clusters
        self.new_celltype = new_celltype
        self.predicate = predicate
        self.metadata = metadata
        self.save_h5ad = save_h5ad

    def run_RankGeneDefault(
        self,
        groups: Optional[Union[str, List[str]]] = None,
        reference: str = "rest",
        copy: bool = False,
    ):
        """
        Perform findmarker analysis for all groups or specific groups compared to a reference.

        Parameters:
        -----------
        groups : Optional[Union[str, List[str]]]
            Group(s) to compare. None for findmarker analysis.
        reference : str (default: 'rest')
            Which group to use as reference.
        copy : bool (default: False)
            Whether to return a copy of the AnnData object.

        Returns:
        --------
        AnnData object with differential expression results.
        """
        check_mkdir(self.outdir)
        self.adata = oep.loadH5AD(
            self.adata,
            sampleid=self.sampleid,
            group=self.group,
            clusters=self.clusters,
            new_celltype=self.new_celltype,
            predicate=self.predicate,
            metadata=self.metadata,
        )

        if groups is None:
            LOGGER.info("Performing findmarker analysis for all groups...")
        else:
            if isinstance(groups, str):
                groups = [groups]
            LOGGER.info(f"Performing rank genes between {groups} and {reference}...")

        kwargs = {
            "adata": self.adata,
            "groupby": self.groupby,
            "method": self.method,
            "pts": self.pts,
            "use_raw": self.use_raw,
            "layer": "normalised",
            "corr_method": self.corr_method,
            "tie_correct": self.tie_correct,
            "rankby_abs": self.rankby_abs,
            "key_added": self.key_added,
        }

        if groups is not None:
            kwargs["groups"] = groups
            if reference:
                kwargs["reference"] = reference

        rank_genes_groups(**kwargs)

        if self.save_h5ad:
            self.adata.write_h5ad(os.path.join(self.outdir, "PRO_diff.h5ad"))
        return self.adata.copy() if copy else self.adata

    def run_RankGeneSpecified(
        self,
        group1: List[str],
        group2: List[str],
        copy: bool = False,
    ):
        """
        Perform differential expression analysis between two specified groups.

        Parameters:
        -----------
        group1 : List[str]
            First group of clusters.
        group2 : List[str]
            Second group of clusters.
        copy : bool (default: False)
            Whether to return a copy of the AnnData object.

        Returns:
        --------
        AnnData object with differential expression results.
        """
        check_mkdir(self.outdir)
        self.adata = oep.loadH5AD(
            self.adata,
            sampleid=self.sampleid,
            group=self.group,
            clusters=self.clusters,
            new_celltype=self.new_celltype,
            predicate=self.predicate,
            metadata=self.metadata,
        )

        LOGGER.info(
            f"Performing rank genes between {'_'.join(group1)} vs {'_'.join(group2)}..."
        )

        overlap = set(group1) & set(group2)
        if overlap:
            raise ValueError(f"Group1 and Group2 have overlapping clusters: {overlap}")

        # Create a temporary grouping in adata for the comparison
        self.adata.obs["_comparison_group"] = "rest"
        for g in group1:
            self.adata.obs.loc[
                self.adata.obs[self.groupby].isin([g]), "_comparison_group"
            ] = "_".join(group1)
        for g in group2:
            self.adata.obs.loc[
                self.adata.obs[self.groupby].isin([g]), "_comparison_group"
            ] = "_".join(group2)

        rank_genes_groups(
            self.adata,
            groupby="_comparison_group",
            groups=["_".join(group1)],
            reference="_".join(group2),
            method=self.method,
            pts=self.pts,
            use_raw=self.use_raw,
            layer="normalised",
            corr_method=self.corr_method,
            tie_correct=self.tie_correct,
            rankby_abs=self.rankby_abs,
            key_added=self.key_added,
        )
        # Remove temporary grouping
        del self.adata.obs["_comparison_group"]

        if self.save_h5ad:
            self.adata.write_h5ad(os.path.join(self.outdir, "PRO_diff.h5ad"))
        return self.adata.copy() if copy else self.adata

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        LOGGER.info("rankgenes done")


@log_function_call
def rankgenes(args):
    with RankGenes(
        adata=args.input,
        outdir=args.outdir,
        groupby=args.groupby,
        use_raw=args.use_raw,
        pts=args.pts,
        method=args.method,
        corr_method=args.corr_method,
        tie_correct=args.tie_correct,
        rankby_abs=args.rankby_abs,
        key_added=args.key_added,
        sampleid=args.sampleid,
        group=args.group,
        clusters=args.clusters,
        new_celltype=args.new_celltype,
        predicate=args.predicate,
        metadata=args.metadata,
        save_h5ad=args.save_h5ad,
    ) as runner:

        if args.group1 is not None:
            adata = runner.run_RankGeneSpecified(
                group1=args.group1,
                group2=args.group2,
            )
        else:
            adata = runner.run_RankGeneDefault()
        write_diff_genes(
            adata,
            outdir=args.outdir,
            logfc=np.log2(args.FC),
            minpct=args.minpct,
            pvals=args.pvals,
            refgenome=args.refgenome,
            pvals_adj=args.qvals,
            groupby=args.groupby,
            top_n=args.top_n,
        )


def get_opts_rankgenes(parser, sub_program=False):
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default=None,
        help="Input file (AnnData object in .h5ad format)",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        default=None,
        help="Output directory for results",
        required=True,
    )
    parser.add_argument(
        "--sampleid",
        type=str,
        default=None,
        help="Column name in adata.obs for sample IDs",
    )
    parser.add_argument(
        "--group",
        type=str,
        default=None,
        help="Column name in adata.obs for group information",
    )
    parser.add_argument(
        "--clusters",
        type=str,
        default=None,
        help="Column name in adata.obs for cluster information",
    )
    parser.add_argument(
        "--new_celltype",
        type=str,
        default=None,
        help="Column name in adata.obs for new cell type information",
    )
    parser.add_argument(
        "--predicate",
        type=str,
        default=None,
        help="Predicate for filtering adata.obs, e.g. (sampleid in ['STB1', 'STB4']) and (clusters in ['1', '5', '6'])",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        help="",
    )
    parser.add_argument(
        "--groupby",
        type=str,
        default="clusters",
        help="Column name in adata.obs for grouping cells (required for differential expression analysis)",
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=["logreg", "t-test", "wilcoxon", "t-test_overestim_var"],
        default="wilcoxon",
        help="Method for differential expression analysis (default: 'wilcoxon')",
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
        type=lambda x: bool(strtobool(x)),
        default=False,
        help="Whether to apply tie correction in the ranking algorithm (default: False)",
    )
    parser.add_argument(
        "--rankby_abs",
        type=lambda x: bool(strtobool(x)),
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
        "--group1",
        type=str,
        nargs="+",
        default=None,
        help="First group of clusters for specified comparison",
    )
    parser.add_argument(
        "--group2",
        type=str,
        nargs="+",
        default=None,
        help="Second group of clusters for specified comparison",
    )
    parser.add_argument(
        "--save_h5ad",
        type=lambda x: bool(strtobool(x)),
        default=True,
        help="Whether to save the updated AnnData object to .h5ad format (default: True)",
    )
    parser.add_argument(
        "--pts",
        type=lambda x: bool(strtobool(x)),
        default=True,
        help="Whether to calculate percentage of cells expressing each gene (default: True)",
    )
    parser.add_argument(
        "--use_raw",
        type=lambda x: bool(strtobool(x)),
        default=False,
        help="Whether to use raw data for differential expression analysis (default: False)",
    )
    parser.add_argument(
        "--FC",
        type=float,
        default=1.5,
        help="fold change threshold for filtering genes (default: 1.5)",
    )
    parser.add_argument(
        "--minpct",
        type=float,
        default=0.25,
        help="Minimum percentage of cells expressing a gene (default: 0.25)",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=3,
        help="",
    )
    parser.add_argument(
        "--pvals",
        type=float,
        default=0.05,
        help="P-value threshold for filtering genes (default: 0.05)",
    )
    parser.add_argument(
        "--qvals",
        type=float,
        default=None,
        help="Q-value threshold for filtering genes (default: 0.1)",
    )
    parser.add_argument(
        "--refgenome",
        type=str,
        default=None,
        help="Reference genome for annotation (default: None)",
    )
    return parser


if __name__ == "__main__":
    unittest.main()
