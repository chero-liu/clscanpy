from anndata import AnnData
import unittest
import os
import oescanpy as oep
from oescanpy.tools.utils import (
    check_adata,
    retrieve_layers,
    save_parameters_to_anndata,
    save_dict_to_txt,
    check_mkdir,
)
from oescanpy.tools.preprocessing.highly_variable_genes import variable_filter
from oescanpy.tools.embedding._batch import batch_correction
from oescanpy.tools.preprocessing.graph import neighbors, clustering, find_elbow
from oescanpy.tools.preprocessing.utils import transform_cluster_labels, umap_argue
from oescanpy.tools.embedding.embedding import run_umap, run_tsne
from oescanpy.tools.preprocessing.basicpreprocessing import normalize
from typing import Union
import random
from distutils.util import strtobool
from oescanpy.log import log_function_call

import logging
import pandas as pd
LOGGER = logging.getLogger(__name__)


class Cluster:

    def __init__(
        self,
        adata: AnnData,
        hvg_num: int = 2000,
        n_pcs: int = 10,
        batch_method: str = "No",
        batch_name: str = "batchid",
        clust_method: str = "leiden",
        min_dist: str = "default",
        resolution: float = 0.4,
        target_sum: int = 10000,
        outdir: str = None,
        sampleid: Union[str, list] = None,
        group: Union[str, list] = None,
        clusters: Union[str, list] = None,
        new_celltype: Union[str, list] = None,
        predicate: Union[str, list] = None,
        metadata: Union[str] = None,
        save_h5ad: bool = False,
        hvg_add: str = None,
        hvg_remove: str = None,
        hvg_replace: str = None,
        n_comps: int = 50,
    ):
        self.adata = adata
        self.hvg_num = hvg_num
        self.n_pcs = n_pcs
        self.batch_method = batch_method
        self.batch_name = batch_name
        self.clust_method = clust_method
        self.min_dist = min_dist
        self.resolution = resolution
        self.target_sum = target_sum
        self.outdir = outdir
        self.sampleid = sampleid
        self.group = group
        self.clusters = clusters
        self.new_celltype = new_celltype
        self.predicate = predicate
        self.metadata = metadata
        self.save_h5ad = save_h5ad
        self.hvg_add = hvg_add
        self.hvg_remove = hvg_remove
        self.hvg_replace = hvg_replace
        self.n_comps = n_comps

    def run(self):
        LOGGER.info(f"Start cluster ...")
        check_mkdir(self.outdir)
        random.seed(12345)

        self.adata = oep.loadH5AD(
            self.adata,
            sampleid=self.sampleid,
            group=self.group,
            clusters=self.clusters,
            new_celltype=self.new_celltype,
            predicate=self.predicate,
            metadata=self.metadata,
        )

        self.adata = check_adata(self.adata)

        self.adata.raw = self.adata
        LOGGER.info("Raw adata is already saved in adata.raw")

        normalize(self.adata, target_sum=self.target_sum)

        self.adata.layers["normalised"] = self.adata.X.copy()
        LOGGER.info(
            "Normalization step is finished ,copy to adata.layers['normalised']"
        )

        data = self.adata.copy()

        self.adata = variable_filter(
            self.adata,
            "normalised",
            self.hvg_num,
            flavor="seurat",
        )
        # 处理自定义高变基因列表
        self._process_custom_hvgs()

        variable_filter_result = self.adata.var.copy()
        self.adata = self.adata[:, self.adata.var.highly_variable]

        #如果batch_name不是Batchid根据batch_name去更改batchid这一列
        if self.batch_name != "batchid" and self.batch_method!="No":
            self.adata.obs['batchid'] = pd.Categorical(self.adata.obs[self.batch_name].cat.codes + 1)
        self.adata, use_rep = batch_correction(
            self.adata,
            batch_key=self.batch_name,
            methods=self.batch_method,
            n_pcs=self.n_comps,
        )

        if self.n_pcs == "auto":
            self.n_pcs = find_elbow(
                self.adata,
                reduction="pca",
            )
        else:
            self.n_pcs = int(self.n_pcs)

        neighbors(
            self.adata,
            n_neighbors=15,
            n_pcs=self.n_pcs,
            use_rep=use_rep,
        )

        clustering(
            self.adata,
            self.clust_method,
            self.resolution,
        )

        transform_cluster_labels(
            self.adata,
            self.clust_method,
        )

        if self.min_dist == "default":
            self.min_dist = umap_argue(self.adata)
        else:
            self.min_dist = float(self.min_dist)

        run_umap(
            self.adata,
            min_dist=self.min_dist,
        )

        run_tsne(
            self.adata,
            n_pcs=self.n_pcs,
        )

        self.adata = retrieve_layers(
            self.adata,
            data,
        )
        self.adata.var = variable_filter_result
        argu = {
            "n_pcs": str(self.n_pcs),
            "n_comps": str(self.n_comps),
            "resolution": str(self.resolution),
            "batch_method": self.batch_method,
            "clust_method": self.clust_method,
            "min_dist": str(self.min_dist),
            "hvg_num": str(self.hvg_num),
            "batch_name": self.batch_name,
            "target_sum": str(self.target_sum),
        }
        save_dict_to_txt(
            argu,
            os.path.join(self.outdir, "arguments.txt"),
        )
        self.adata = save_parameters_to_anndata(self.adata, argu)
        self.adata = oep.loadH5AD(self.adata)

        if self.save_h5ad:
            self.adata.write_h5ad(os.path.join(self.outdir, "adata.h5ad"))
        return self.adata

    def _process_custom_hvgs(self):
        """
        处理自定义高变基因列表，支持添加、删除或替换高变基因
        """
        # 创建基因名映射字典，处理大小写不敏感匹配
        gene_map = {gene.lower(): gene for gene in self.adata.var_names}

        def match_genes(target_genes):
            """匹配目标基因到实际基因名（大小写不敏感）"""
            matched_genes = []
            for gene in target_genes:
                if gene in gene_map:
                    matched_genes.append(gene_map[gene])
                elif gene.lower() in gene_map:
                    matched_genes.append(gene_map[gene.lower()])
            return matched_genes

        # 替换高变基因
        if self.hvg_replace:
            with open(self.hvg_replace, 'r') as f:
                target_genes = [line.strip() for line in f if line.strip()]

            matched_genes = match_genes(target_genes)
            # 重置所有高变基因状态
            self.adata.var['highly_variable'] = False
            # 设置指定基因为高变基因
            self.adata.var.loc[matched_genes, 'highly_variable'] = True
            LOGGER.info(f"Replaced highly variable genes with {len(matched_genes)} genes from {self.hvg_replace}")

        # 添加高变基因
        if self.hvg_add:
            with open(self.hvg_add, 'r') as f:
                target_genes = [line.strip() for line in f if line.strip()]

            matched_genes = match_genes(target_genes)
            # 添加基因到现有高变基因集合
            current_hvgs = self.adata.var_names[self.adata.var['highly_variable']].tolist()
            new_hvgs = list(set(current_hvgs + matched_genes))
            self.adata.var['highly_variable'] = self.adata.var_names.isin(new_hvgs)
            LOGGER.info(f"Added {len(matched_genes)} genes to highly variable genes from {self.hvg_add}")

        # 删除高变基因
        if self.hvg_remove:
            with open(self.hvg_remove, 'r') as f:
                target_genes = [line.strip() for line in f if line.strip()]

            matched_genes = match_genes(target_genes)
            # 从现有高变基因中删除指定基因
            current_hvgs = self.adata.var_names[self.adata.var['highly_variable']].tolist()
            remaining_hvgs = [gene for gene in current_hvgs if gene not in matched_genes]
            self.adata.var['highly_variable'] = self.adata.var_names.isin(remaining_hvgs)
            LOGGER.info(f"Removed {len(matched_genes)} genes from highly variable genes using {self.hvg_remove}")
    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        LOGGER.info("cluster done")


@log_function_call
def cluster(args):
    with Cluster(
        adata=args.input,
        hvg_num=args.hvg_num,
        n_pcs=args.n_pcs,
        batch_method=args.batch_method,
        batch_name=args.batch_name,
        clust_method=args.clust_method,
        min_dist=args.min_dist,
        resolution=args.resolution,
        target_sum=args.target_sum,
        outdir=args.outdir,
        sampleid=args.sampleid,
        group=args.group,
        clusters=args.clusters,
        new_celltype=args.new_celltype,
        predicate=args.predicate,
        metadata=args.metadata,
        save_h5ad=args.save_h5ad,
        hvg_add=args.hvg_add,
        hvg_remove=args.hvg_remove,
        hvg_replace=args.hvg_replace,
        n_comps=args.n_comps,
    ) as runner:
        runner.run()


def get_opts_cluster(parser, sub_program=False):
    parser.add_argument(
        "-i", "--input", type=str, default=None, help="Input file", required=True
    )
    parser.add_argument(
        "-o", "--outdir", type=str, default=None, help="Output directory", required=True
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
        help="predicate for filtering adata.obs, e.g. (sampleid in ['STB1', 'STB4']) and (clusters in ['1', '5', '6'])",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        help="",
    )
    parser.add_argument("--hvg_num", type=int, default=2000, help="hvg genes number")
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
    parser.add_argument(
        "--save_h5ad",
        type=lambda x: bool(strtobool(x)),
        default=True,
        help="",
    )
    parser.add_argument(
        "--hvg_add",
        type=str,
        default=None,
        help="Gene list file for adding highly variable genes (one gene per line)",
    )
    parser.add_argument(
        "--hvg_remove",
        type=str,
        default=None,
        help="Gene list file for removing highly variable genes (one gene per line)",
    )
    parser.add_argument(
        "--hvg_replace",
        type=str,
        default=None,
        help="Gene list file for replacing highly variable genes (one gene per line)",
    )
    return parser


if __name__ == "__main__":
    unittest.main()
