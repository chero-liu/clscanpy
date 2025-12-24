import scanpy as sc
import pandas as pd
import numpy as np
from collections import Counter
from clscanpy.tools.read.utils import get_single_data_format, read, ensure_sparse_matrix
from clscanpy.tools.utils import (
    check_mkdir,
    save_parameters_to_anndata,
    save_dict_to_txt,
)
from clscanpy.tools.preprocessing.basicpreprocessing import (
    _filter,
    find_mt_threshold,
    calculate_mt_common,
    _load_gene_list,
)
import os
from clscanpy.log import log_function_call

MT_THRESHOLD = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]


class Create:
    def __init__(
        self,
        info: None,
        outdir: str,
        prefix: str,
        refgenome: str,
        mtfilter: str = "default",
        hbfilter: float = 0.05,
        mingene: float = 200,
        minumi: float = 0,
        maxumi: float = 0.98,
        maxgene: float = 0.98,
        mincell: int = 5,
        dpi: int = 300,
        save_h5ad: bool = False,
    ):
        self.info = info
        self.outdir = outdir
        self.prefix = prefix
        self.refgenome = refgenome
        self.mtfilter = mtfilter
        self.hbfilter = hbfilter
        self.mingene = mingene
        self.minumi = minumi
        self.maxgene = maxgene
        self.maxumi = maxumi
        self.mincell = mincell
        self.dpi = dpi
        self.save_h5ad = save_h5ad

    def run(self):
        self.info = pd.read_csv(self.info, sep=",")
        adatas, mt_list, genes_list = [], [], []
        raw_sample, raw_cells, raw_genes, raw_umi, raw_gene = [], [], [], [], []
        filter_sample, filter_cells, filter_genes, filter_umi, filter_gene = (
            [],
            [],
            [],
            [],
            [],
        )

        for index, line in self.info.iterrows():
            print(f"Start processing {line.spname} ...")
            data_type = get_single_data_format(line.path)
            adata = read(line.path, prefix=line.spname, dtype=data_type)
            ensure_sparse_matrix(adata)

            df = sc.pp.calculate_qc_metrics(
                adata, percent_top=None, log1p=False, inplace=False
            )[0]
            raw_sample.append(line.spname)
            raw_cells.append(adata.shape[0])
            raw_genes.append(adata.shape[1])
            raw_umi.append(np.median(df["total_counts"]))
            raw_gene.append(np.median(df["n_genes_by_counts"]))

            _filter(
                adata,
                mingenes=self.mingene,
                minumis=self.minumi,
                maxgenes=self.maxgene,
                maxumis=self.maxumi,
                mincells=self.mincell,
            )

            add_names = line.index.drop(["path"])
            for tag in add_names:
                adata.obs[tag] = str(line[tag])
            mt_ = find_mt_threshold(
                adata,
                mt_thresholds=MT_THRESHOLD,
                refgenome=self.refgenome,
                mtfilter=self.mtfilter,
            )
            mt_list.append(mt_)
            genes_list.append(adata.shape[1])
            adatas.append(adata)

        adata = sc.concat(adatas, join="outer")
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
        adata.X = np.nan_to_num(adata.X)

        mt_ = find_mt_threshold(
            adata,
            mt_thresholds=MT_THRESHOLD,
            refgenome=self.refgenome,
            mtfilter=self.mtfilter,
        )
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )
        mt_common = calculate_mt_common(self.mtfilter, mt_list)
        mt_dict = {}
        mt_raw = Counter(adata.obs["spname"])
        adata = adata[adata.obs.pct_counts_mt < int(mt_common), :]
        mt_filter = Counter(adata.obs["spname"])
        for key, value in mt_raw.items():
            mt_dict[key] = mt_filter[key] / mt_raw[key]
        mt_df = pd.DataFrame(
            {"SampleID": mt_dict.keys(), "mt_filtered_percent": mt_dict.values()}
        )

        adata.var["hb"] = _load_gene_list(adata, f"{self.refgenome}/HB_genelist.gmt")
        sc.pp.calculate_qc_metrics(
            adata,
            percent_top=None,
            log1p=False,
            qc_vars=["hb"],
            inplace=True,
        )
        adata.obs["pct_counts_hb"] /= 100
        hb_dict = {}
        hb_raw = Counter(adata.obs["spname"])
        adata = adata[adata.obs["pct_counts_hb"] < self.hbfilter, :]
        hb_filter = Counter(adata.obs["spname"])
        for key, value in mt_raw.items():
            hb_dict[key] = hb_filter[key] / hb_raw[key]
        hb_df = pd.DataFrame(
            {"SampleID": hb_dict.keys(), "hb_filtered_percent": hb_dict.values()}
        )

        # re-order obs slot
        reorder_slot = ["rawname", "spname", "gname"]
        for slot in reorder_slot:
            adata.obs[slot] = adata.obs[slot].astype("category")
            if slot == "gname":
                try:
                    order = pd.DataFrame(self.info.gname.tolist()).drop_duplicates()
                except AttributeError:
                    order = pd.DataFrame(self.info.spname.tolist()).drop_duplicates()
                order = order[0].tolist()
            else:
                order = self.info[slot].tolist()

            order = [str(x) for x in order]
            adata.obs[slot] = adata.obs[slot].cat.set_categories(order)

        adata.obs.rename(columns={"spname": "sample"}, inplace=True)

        del adata.obs["n_genes"]
        del adata.obs["n_counts"]

        for sample in set(adata.obs["sample"]):
            tmp = adata[adata.obs["sample"] == sample]
            filter_sample.append(sample)
            filter_cells.append(tmp.shape[0])
            filter_genes.append(tmp.shape[1])
            filter_umi.append(np.median(tmp.obs["total_counts"]))
            filter_gene.append(np.median(tmp.obs["n_genes_by_counts"]))

        raw_df = pd.DataFrame(
            dict(
                SampleID=raw_sample,
                raw_cells=raw_cells,
                raw_genes=raw_genes,
                raw_median_umi=raw_umi,
                raw_median_gene=raw_gene,
            )
        )
        filter_df = pd.DataFrame(
            dict(
                SampleID=filter_sample,
                filtered_cells=filter_cells,
                filtered_genes=filter_genes,
                filtered_median_umi=filter_umi,
                filtered_median_gene=filter_gene,
            )
        )
        qc_df = pd.merge(raw_df, filter_df, on="SampleID")
        qc_df = pd.merge(
            qc_df,
            pd.merge(mt_df, hb_df, on="SampleID"),
            on="SampleID",
        )
        check_mkdir(self.outdir)
        qc_df.to_csv(
            self.outdir + "/" + "{0}_qc_config.xls".format(self.prefix),
            sep="\t",
            index=False,
        )

        argu = {
            "mingene": str(self.mingene),
            "minumi": str(self.minumi),
            "maxgene": str(self.maxgene),
            "maxumi": str(self.maxumi),
            "mincell": str(self.mincell),
            "mtfilter": str(mt_common),
            "hbfilter": str(self.hbfilter),
        }
        save_dict_to_txt(
            argu,
            os.path.join(self.outdir, "arguments.txt"),
        )
        adata = save_parameters_to_anndata(adata, argu)

        if self.save_h5ad:
            adata.write_h5ad(os.path.join(self.outdir, "filtered.h5ad"))
        return adata

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        print("create done")


@log_function_call
def create(args):
    with Create(
        info=args.info,
        outdir=args.outdir,
        prefix=args.prefix,
        refgenome=args.refgenome,
        mtfilter=args.mtfilter,
        hbfilter=args.hbfilter,
        mingene=args.mingene,
        minumi=args.minumi,
        maxumi=args.maxumi,
        maxgene=args.maxgene,
        mincell=args.mincell,
        save_h5ad=args.save_h5ad,
    ) as runner:
        runner.run()


def get_opts_create(parser, sub_program=False):
    parser.add_argument(
        "-i",
        "--info",
        type=str,
        default=None,
        help="Path to the sample information table (CSV format) containing details of each sample",
        required=True,
    )

    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        default=None,
        help="Path to the output directory",
        required=True,
    )

    parser.add_argument(
        "--prefix",
        type=str,
        default="",
        help="Prefix for output files. If not specified, the output directory name will be used",
    )

    parser.add_argument(
        "--refgenome",
        type=str,
        default=None,
        help="",
    )

    parser.add_argument(
        "--mtfilter",
        type=str,
        default="default",
        help="Mitochondrial gene filtering method, default is 'default'",
    )
    parser.add_argument(
        "--hbfilter",
        type=float,
        default=0.05,
        help="Mitochondrial gene filtering method, default is '0.05'",
    )
    parser.add_argument(
        "--mingene",
        type=float,
        default=200,
        help="Minimum number of genes that must be expressed per cell, default is 200",
    )

    parser.add_argument(
        "--maxgene",
        type=float,
        default=0.98,
        help="Maximum number of genes expressed per cell, default is 7500",
    )

    parser.add_argument(
        "--minumi",
        type=float,
        default=0,
        help="Minimum number of UMIs (Unique Molecular Identifiers) per cell, default is 1000",
    )

    parser.add_argument(
        "--maxumi",
        type=float,
        default=0.98,
        help="Maximum number of UMIs per cell, default is 50000",
    )

    parser.add_argument(
        "--mincell",
        type=int,
        default=5,
        help="Minimum number of cells in which each gene must be expressed, default is 0 (no filtering)",
    )
    parser.add_argument("--save_h5ad", default=True, help="")
    return parser


if __name__ == "__main__":
    unittest.main()
