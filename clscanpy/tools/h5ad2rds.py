import os
import sys
import subprocess
from clscanpy.config import ROOT_PATH
from clscanpy.log import log_function_call


class H5ad2Rds:

    def __init__(
        self,
        input,
        outdir,
    ):
        self.input = input
        self.outdir = outdir

    @property
    def shell_script(self):
        shell_script_content = f"""
source /opt/miniconda3/etc/profile.d/conda.sh && conda activate clscanpy
Rscript {ROOT_PATH}/script/data_transformation/h5ad2rds.r -i {self.input}
"""
        return shell_script_content

    def run(self):
        subprocess.run(self.shell_script, shell=True)

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        print("done")


@log_function_call
def h5ad2rds(args):
    with H5ad2Rds(
        input=args.input,
        outdir=args.outdir,
    ) as runner:
        runner.run()


def get_opts_h5ad2rds(parser, sub_program=True):
    parser.add_argument("-i", "--input", type=str, default=None, help="Input file")
    parser.add_argument("-o", "--outdir", type=str, default=None, help="Outdir file")
    return parser


if __name__ == "__main__":
    unittest.main()
