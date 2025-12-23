import unittest
import os
from clscanpy.tools.read.utils import convert_seurat_to_anndata
import clscanpy as csc
from clscanpy.tools.utils import Step
from clscanpy.log import log_function_call

import logging

LOGGER = logging.getLogger(__name__)


class Rds2H5ad:
    def __init__(self, args):
        self.input = args.input

    def run(self):
        LOGGER.info("Start Rds2H5ad ...")
        adata = convert_seurat_to_anndata(self.input, use_raw_counts=False)
        adata.write_h5ad(os.path.join(os.path.dirname(self.input), "adata.h5ad"))

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        LOGGER.info("done")


@log_function_call
def rds2h5ad(args):
    with Rds2H5ad(args) as runner:
        runner.run()


def get_opts_rds2h5ad(parser, sub_program=True):
    parser.add_argument("-i", "--input", type=str, default=None, help="Input file")
    return parser


if __name__ == "__main__":
    unittest.main()
