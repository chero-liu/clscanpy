import os
import unittest
from oescanpy.tools.utils import Step, s_common
from oescanpy.config import ROOT_PATH
from oescanpy.log import log_function_call

import logging
LOGGER = logging.getLogger(__name__)

class V4tov3(Step):
    def __init__(self, args):
        super().__init__(args)
        self.input = args.input

    def run(self):
        LOGGER.info("Start V4tov3 ...")
        for input in self.input.split(","):
            os.system(
                f"""
                module purge && module load OESingleCell/3.0.d
                Rscript {ROOT_PATH}/script/v4tov3.r \\
                    -i {input} \\
                    -f h5seurat \\
                    -o {os.path.dirname(input)}
                module purge && module load OESingleCell/2.0.0
                Rscript {ROOT_PATH}/script/v4tov3.r \\
                    -o {os.path.dirname(input)}
                """
            )

@log_function_call
def v4tov3(args):
    with V4tov3(args) as runner:
        runner.run()


def get_opts_v4tov3(parser, sub_program=True):
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default="./",
        help="path of clean",
    )

    if sub_program:
        parser = s_common(parser)

    return parser


if __name__ == "__main__":
    unittest.main()
