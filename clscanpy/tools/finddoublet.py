import scanpy as sc
import pandas as pd
import numpy as np
import os
from clscanpy.log import log_function_call
from oescanpy.tools.preprocessing.doublet import find_doublet
import clscanpy as csc
from typing import Union
from clscanpy.tools.utils import (
    check_mkdir,
    save_parameters_to_anndata,
    save_dict_to_txt,
)

class FindDoublet:
    def __init__(
        self,
        input: str = None,
        outdir: str = None,
        prefix: str,
        sample: Union[str, list] = None,
        group: Union[str, list] = None,
        cluster: Union[str, list] = None,
        celltype: Union[str, list] = None,
        method: str = "doubletdetection",
        version: str = "v3",
        save_h5ad: bool = False,
    ):
        self.input = input
        self.outdir = outdir
        self.prefix = prefix
        self.sample = sample
        self.group = group
        self.cluster = cluster
        self.celltype = celltype
        self.method = method
        self.version = version
        self.save_h5ad = save_h5ad

    def run(self):
        check_mkdir(self.outdir)
        adata = csc.loadH5AD(
            self.input,
            sample=self.sample,
            group=self.group,
            cluster=self.cluster,
            celltype=self.celltype,
            predicate=self.predicate,
        )
        doublet = find_doublet(
            adata,
            method=self.method,
            batch_key="sample",
            version=self.version,
        )
        adata.obs["doublet"] = doublet
        adata = adata[~adata.obs["doublet"], :].copy()

        argu = {
            "rmdoublet_method": str(self.method),
            "rmdoublet_method_batch_key": "sample",
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
def finddoublet(args):
    with FindDoublet(
        input=args.input,
        outdir=args.outdir,
        prefix=args.prefix,
    ) as runner:
        runner.run()


def get_opts_finddoublet(parser, sub_program=False):
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default=None,
        help="",
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

    parser.add_argument("--save_h5ad", default=True, help="")
    return parser


if __name__ == "__main__":
    unittest.main()
