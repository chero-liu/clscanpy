import scanpy as sc
import pandas as pd
import numpy as np
from clscanpy.tools.plotting.violin import plot_qc_violin
from clscanpy.tools.color.utils import get_color_order
from clscanpy.tools.color.colors import add_color_to_dataframe
from clscanpy.tools.preprocessing.doublet import findDoublet
from clscanpy.tools.utils import (
    check_mkdir,
    clean_adata,
    save_parameters_to_anndata,
    save_dict_to_txt,
)
from clscanpy.tools.read.utils import (
    summarizeAdataMetrics,
    getSingleDataFormat,
    read,
    ensureSparseMatrix,
    MultiSpeciesSplitter,
)
from clscanpy.tools.read.load import standardize_colors
from clscanpy.tools.preprocessing.basicpreprocessing import (
    _qc,
    _filter,
    _mtfilter,
    find_outliers,
)
import clscanpy as oep
import random
import os
import unittest
from clscanpy.log import log_function_call

import logging

LOGGER = logging.getLogger(__name__)

random.seed(2025)


class Create:

    def __init__(
        self,
        info: str = None,
        outdir: str = None,
        refgenome: str = None,
        mtfilter: str = "default",
        hbfilter: float = 0.05,
        mingene: float = 200,
        maxgene: float = 7500,
        minumi: float = 1000,
        maxumi: float = 50000,
        mincell: int = 0,
        rmdoublet_method: str = "doubletdetection",
        qc_vars_vis: list = None,
        save_h5ad: bool = False,
        version: str = "v3",
        sampleid_palette: str = "cellpaper_custom",
        group_palette: str = "group_default",
        nfold: int = None,
        nfold_vars: list = ["nFeature_RNA", "nCount_RNA", "percent.mito"],
        nfold_type: list = ["both", "both", "higher"],
        cut_1: str = "median",
        cut_2: str = "mad",
    ):
        self.info = info
        self.outdir = outdir
        self.refgenome = refgenome
        self.mtfilter = mtfilter
        self.hbfilter = hbfilter
        self.mingene = mingene
        self.maxgene = maxgene
        self.minumi = minumi
        self.maxumi = maxumi
        self.mincell = mincell
        self.rmdoublet_method = rmdoublet_method
        self.qc_vars_vis = qc_vars_vis
        self.save_h5ad = save_h5ad
        self.version = version
        self.sampleid_palette = sampleid_palette
        self.group_palette = group_palette
        self.nfold = nfold
        self.nfold_vars = nfold_vars
        self.nfold_type = nfold_type
        self.cut_1 = cut_1
        self.cut_2 = cut_2

    def run(self):
        LOGGER.info("Start Create ...")
        check_mkdir(self.outdir)

        if self.info.endswith(".csv"):
            self.info = pd.read_csv(self.info)
            # 按 sampleid_order 排序信息表格
            if "sampleid_order" in self.info.columns:
                self.info = self.info.sort_values("sampleid_order")

            # 获取排序后的样本顺序
            sampleid_order = self.info["sampleid"].tolist()

            # 获取分组顺序（按样本顺序首次出现的分组）
            group_order = []
            for group in self.info["group"]:
                if group not in group_order:
                    group_order.append(group)

            adatas = []

            for index, line in self.info.iterrows():
                LOGGER.info(f"Start processing {line.sampleid} ...")
                data_type = getSingleDataFormat(
                    getattr(line, "datatype", None), line.path
                )

                adata = read(
                    line.path,
                    prefix=line.sampleid,
                    dtype=data_type,
                    index=(index + 1),
                )

                ensureSparseMatrix(adata)

                add_names = line.index.drop(["path"])
                for tag in add_names:
                    adata.obs[tag] = str(line[tag])

                if (
                    getattr(line, "mix_species", None)
                    and line.mix_species.lower() != "no"
                ):
                    try:
                        splitter = MultiSpeciesSplitter(
                            input_dir=line.path, mix_species=line.mix_species
                        )
                        adata = splitter.split(adata)
                        LOGGER.info(
                            f"Successfully split species '{line.mix_species}' for {line.path}"
                        )
                    except FileNotFoundError as e:
                        LOGGER.warning(f"Missing required files for {line.path}: {e}")
                    except ValueError as e:
                        LOGGER.warning(f"Invalid data format for {line.path}: {e}")
                    except Exception as e:
                        LOGGER.exception(
                            f"[Error] Unexpected error while splitting {line.path}: {e}"
                        )
                else:
                    LOGGER.debug(
                        f"mix_species not specified or set to 'no' for {line.path}"
                    )

                adatas.append(adata)

            adata = sc.concat(adatas, join="outer", keys=sampleid_order)

            adata.var_names_make_unique()
            adata.obs_names_make_unique()
            adata.X = np.nan_to_num(adata.X)

            # 设置分类顺序
            adata.obs["sampleid"] = pd.Categorical(
                adata.obs["sampleid"], categories=sampleid_order, ordered=True
            )

            adata.obs["group"] = pd.Categorical(
                adata.obs["group"], categories=group_order, ordered=True
            )
        else:
            adata = oep.loadH5AD(self.info)
            adata = clean_adata(
                adata, ["nCount_RNA", "nFeature_RNA", "percent.HB", "percent.mito"]
            )
            del adata.var
            del adata.uns

        adata = standardize_colors(
            adata,
            use_col=["sampleid", "group"],
            palette=[self.sampleid_palette, self.group_palette],
        )

        qc_vars, adata = _qc(
            adata,
            refgenome=self.refgenome,
        )
        if "percent.mito" in qc_vars:
            mitolist = (
                adata.obs.groupby("sampleid")["percent.mito"]
                .quantile(0.9)
                .reset_index(name="q90_percent")
            )
        else:
            self.mtfilter = None

        if "percent.HB" not in qc_vars:
            self.hbfilter = None

        filter_bf = summarizeAdataMetrics(
            adata,
            metrics=qc_vars,
            group_by="sampleid",
            suffix="beforeQC",
        )

        if self.qc_vars_vis is None:
            self.qc_vars_vis = qc_vars

        LOGGER.info(self.qc_vars_vis)
        try:
            plot_qc_violin(
                adata,
                groupby="sampleid",
                metrics=self.qc_vars_vis,
                save_path=f"{self.outdir}/QC_metrics_beforeQC",
                palette=get_color_order(
                    adata,
                    use_col="sampleid",
                ),
            )
        except Exception as e:
            LOGGER.info(f"Warning: Failed to plot QC metrics before QC - {str(e)}")

        if self.nfold is not None:
            outliers_dict = find_outliers(
                adata=adata,
                vars=self.nfold_vars,
                batch="sampleid",
                type=self.nfold_type,
                cut_1=self.cut_1,
                cut_2=self.cut_2,
                n=self.nfold,
            )

            combined_outliers = np.zeros_like(
                outliers_dict[self.nfold_vars[0]], dtype=bool
            )
            for key in outliers_dict:
                combined_outliers = np.logical_or(combined_outliers, outliers_dict[key])
            adata = adata[~combined_outliers, :].copy()
            LOGGER.info(f"过滤后adata形状: {adata.shape}")

        else:
            adata = _filter(
                adata,
                mingenes=self.mingene,
                minumis=self.minumi,
                maxgenes=self.maxgene,
                maxumis=self.maxumi,
                mincells=self.mincell,
                hbfilter=self.hbfilter,
            )

        if self.mtfilter is not None:
            mtfilter_threshold, adata = _mtfilter(
                adata,
                mitolist["q90_percent"],
                mtfilter=self.mtfilter,
            )
        else:
            mtfilter_threshold = None

        filter_af = summarizeAdataMetrics(
            adata,
            metrics=qc_vars,
            group_by="sampleid",
            suffix="afterQC",
        )

        if self.rmdoublet_method != "None":
            doublet = findDoublet(
                adata,
                method=self.rmdoublet_method,
                batch_key="sampleid",
                version=self.version,
            )
            adata.obs["doublet"] = doublet
            adata = adata[~adata.obs["doublet"], :].copy()
            rmdoublet_af = summarizeAdataMetrics(
                adata,
                metrics=qc_vars,
                group_by="sampleid",
                suffix="after_rmdoublets",
            )

            pd.merge(
                pd.merge(
                    filter_bf[["sampleid", "cell_count_beforeQC"]],
                    filter_af,
                    on="sampleid",
                )[["sampleid", "cell_count_beforeQC", "cell_count_afterQC"]],
                rmdoublet_af,
                on="sampleid",
            ).to_csv(
                self.outdir + "/" + "cell_count_before_after_QC.xls",
                sep="\t",
                index=False,
            )

            pd.merge(
                pd.merge(
                    filter_bf,
                    filter_af,
                    on="sampleid",
                ),
                rmdoublet_af,
                on="sampleid",
            ).to_csv(
                self.outdir + "/" + "cell_statitics_before_after_QC.xls",
                sep="\t",
                index=False,
            )

        else:

            pd.merge(
                filter_bf[["sampleid", "cell_count_beforeQC"]],
                filter_af,
                on="sampleid",
            ).to_csv(
                self.outdir + "/" + "cell_count_before_after_QC.xls",
                sep="\t",
                index=False,
            )

            pd.merge(
                filter_bf,
                filter_af,
                on="sampleid",
            ).to_csv(
                self.outdir + "/" + "cell_statitics_before_after_QC.xls",
                sep="\t",
                index=False,
            )

        try:
            plot_qc_violin(
                adata,
                groupby="sampleid",
                metrics=self.qc_vars_vis,
                save_path=f"{self.outdir}/QC_metrics_afterQC",
                palette=get_color_order(
                    adata,
                    use_col="sampleid",
                ),
            )
        except Exception as e:
            LOGGER.info(f"Warning: Failed to plot QC metrics after QC - {str(e)}")

        argu = {
            "mingene": str(self.mingene),
            "minumi": str(self.minumi),
            "maxgene": str(self.maxgene),
            "maxumi": str(self.maxumi),
            "mincell": str(self.mincell),
            "mtfilter": str(mtfilter_threshold),
            "hbfilter": str(self.hbfilter),
            "refgenome": str(self.refgenome),
            "nfold": str(self.nfold),
            "nfold_vars": str(self.nfold_vars),
            "nfold_type": str(self.nfold_type),
            "cut_1": str(self.cut_1),
            "cut_2": str(self.cut_2),
            "rmdoublet_method": str(self.rmdoublet_method),
            "version": str(self.version),
        }
        save_dict_to_txt(
            argu,
            os.path.join(self.outdir, "arguments.txt"),
        )
        adata = save_parameters_to_anndata(adata, argu)

        adata = oep.loadH5AD(adata)
        if self.save_h5ad:
            adata = clean_adata(adata)
            adata.write_h5ad(os.path.join(self.outdir, "filtered.h5ad"))
            metadata = adata.obs
            metadata.insert(0, "Barcode", metadata.index)
            metadata.to_csv(f"{self.outdir}/metadata.tsv", index=False, sep="\t")
        return adata

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        LOGGER.info("create done")


@log_function_call
def create(args):
    with Create(
        info=args.info,
        outdir=args.outdir,
        refgenome=args.refgenome,
        mtfilter=args.mtfilter,
        hbfilter=args.hbfilter,
        mingene=args.mingene,
        maxgene=args.maxgene,
        minumi=args.minumi,
        maxumi=args.maxumi,
        mincell=args.mincell,
        rmdoublet_method=args.rmdoublet_method,
        qc_vars_vis=args.qc_vars_vis,
        save_h5ad=args.save_h5ad,
        version=args.version,
        sampleid_palette=args.sampleid_palette,
        group_palette=args.group_palette,
        nfold=args.nfold,
        nfold_vars=args.nfold_vars,
        nfold_type=args.nfold_type,
        cut_1=args.cut_1,
        cut_2=args.cut_2,
    ) as runner:
        runner.run()


def get_opts_create(parser, sub_program=False):
    parser.add_argument(
        "-i",
        "--info",
        type=str,
        default=None,
        help="样本信息表的文件路径（CSV格式），包含各样本的信息",
        required=True,
    )

    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        default=None,
        help="",
        required=True,
    )

    parser.add_argument(
        "--refgenome",
        type=str,
        default=None,
        help="参考基因组文件路径，用于线粒体和血红蛋白基因的过滤，不提供时将不进行相关过滤",
        required=True,
    )

    parser.add_argument(
        "--mtfilter",
        type=str,
        default="default",
        help="",
    )

    parser.add_argument(
        "--hbfilter",
        type=float,
        default=0.05,
        help="默认为0.05",
    )

    parser.add_argument(
        "--mingene",
        type=float,
        default=200,
        help="每个细胞必须表达的最少基因数，默认为200",
    )

    parser.add_argument(
        "--maxgene",
        type=float,
        default=7500,
        help="",
    )

    parser.add_argument(
        "--minumi",
        type=float,
        default=1000,
        help="每个细胞必须表达的最少UMI（唯一分子标识符）数，默认为1000",
    )

    parser.add_argument(
        "--maxumi",
        type=float,
        default=50000,
        help="",
    )

    parser.add_argument(
        "--mincell",
        type=int,
        default=0,
        help="每个基因必须在最少多少个细胞中表达，默认为0（不进行过滤）",
    )

    parser.add_argument(
        "--qc_vars_vis",
        type=str,
        nargs="+",
        default=None,
        help="['percent.HB', 'percent.mito', 'nFeature_RNA', 'nCount_RNA']",
    )

    parser.add_argument(
        "--rmdoublet_method",
        type=str,
        default="doubletdetection",
        choices=["doubletdetection", "scrublet", "None"],
        help="检测双细胞的方法，默认为doubletdetection",
    )
    parser.add_argument("--save_h5ad", default=True, help="")
    parser.add_argument(
        "--version",
        type=str,
        default="v3",
        help="cellranger chemistry的版本",
    )
    parser.add_argument(
        "--sampleid_palette",
        type=str,
        default="cellpaper_custom",
        help="样本ID的调色板，默认为'cellpaper_custom'",
    )
    parser.add_argument(
        "--group_palette",
        type=str,
        default="group_default",
        help="分组的调色板，默认为'group_default'",
    )
    parser.add_argument(
        "--nfold",
        type=int,
        default=None,
        help="Number of folds for cross-validation (if applicable). When set, outlier detection will be performed in a cross-validation manner. Defaults to None (no cross-validation).",
    )
    parser.add_argument(
        "--nfold_vars",
        type=str,
        nargs="+",
        default=["nFeature_RNA", "nCount_RNA"],
        help="List of variables to perform outlier detection on. These variables must exist in the observation matrix (adata.obs) of the single-cell dataset. Defaults to common QC metrics: ['nFeature_RNA', 'nCount_RNA', 'percent.mito'].",
    )
    parser.add_argument(
        "--nfold_type",
        type=str,
        nargs="+",
        default=["both", "both"],
        help="List of outlier detection directions, one-to-one corresponding with --nfold_vars. Allowed values: 'both' (detect both lower and higher outliers), 'lower' (only detect lower outliers), 'higher' (only detect higher outliers). Defaults to ['both', 'both', 'higher'], matching the default --nfold_vars.",
    )
    parser.add_argument(
        "--cut_1",
        type=str,
        default="median",
        help="Type of central tendency statistic used to calculate thresholds. Allowed values: 'median' (robust to outliers, default), 'mean' (sensitive to outliers).",
    )
    parser.add_argument(
        "--cut_2",
        type=str,
        default="mad",
        help="Type of dispersion statistic used to calculate thresholds. Allowed values: 'mad' (median absolute deviation, robust to outliers, default), 'sd' (standard deviation, sensitive to outliers), 'median' (use median as dispersion).",
    )
    return parser


if __name__ == "__main__":
    unittest.main()
