import os
import sys
import subprocess
from oescanpy.config import ROOT_PATH
from oescanpy.log import log_function_call
from oescanpy.tools.utils import check_mkdir
import logging
import unittest

LOGGER = logging.getLogger(__name__)

class Autoanno:

    def __init__(
        self,
        input: str,
        model_path: str = None,
        outdir: str = None,
        mode="best_match",
        cluster_key="clusters",
    ):
        self.input = input
        self.model_path = model_path
        self.outdir = outdir
        self.mode = mode
        self.cluster_key = cluster_key

    def subprocess_script(self, shell_script_content):
        print(shell_script_content)
        subprocess.run(shell_script_content, shell=True)

    @property
    def shell_script(self):
        shell_script_content = f"""
module purge
source /gpfs/oe-scrna/liuchenglong/anaconda3/bin/activate /data/software/conda_envs/scrna_envs/oenv

python {ROOT_PATH}/script/celltypist_annotate/step2.run_celltypist_annotate.py \\
       --input {self.input} \\
       --model-path {self.model_path} \\
       --output-path {self.outdir} \\
       --mode {self.mode} \\
       --majority-voting \\
       --plot-results TRUE \\
       --cluster-key {self.cluster_key}

"""
        return shell_script_content

    def run(self):
        check_mkdir(self.outdir)
        self.subprocess_script(self.shell_script)

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        LOGGER.info("autoanno done")


@log_function_call
def autoanno(args):
    with Autoanno(
        input=args.input,
        model_path=args.model_path,
        outdir=args.outdir,
        mode=args.mode,
        cluster_key=args.cluster_key,
    ) as runner:
        runner.run()


def get_opts_autoanno(parser, sub_program=False):
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="输入文件路径 (h5ad,rds, or h5seurat 格式)",
    )

    parser.add_argument(
        "--model_path",
        type=str,
        default=None,
        help="预训练模型路径，若不提供则根据物种自动选择默认模型",
        required=True,
    )
    parser.add_argument(
        "-o", "--outdir", type=str, default=None, help="Output directory", required=True
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["best_match", "prob_match"],
        default="best_match",
        help="注释模式: 'best_match' (最佳匹配) 或 'prob_match' (概率匹配)，默认为 'best_match'",
    )
    parser.add_argument(
        "--cluster_key",
        type=str,
        default="clusters",
        help="包含聚类信息的 obs 列名，默认为 'clusters'",
    )
    return parser


if __name__ == "__main__":
    unittest.main()
