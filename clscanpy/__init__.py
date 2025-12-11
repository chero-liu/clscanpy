import os
import logging
from clscanpy.tools.read.load import get_h5ad
from clscanpy.tools.read.load import loadH5AD
from clscanpy.tools.plotting.dim import plot_dim
from clscanpy.tools.plotting.density import plot_density
from clscanpy.tools.gcr import get_clusters_result as gcr
from clscanpy.log import create_logging_handler
from clscanpy.tools.color.utils import get_color_order

__VERSION__ = "0.1.0"
__author__ = "liuchenglong"

ASSAY_LIST = [
    "tools",
    "integration",
    # 'differentiation',
    # 'featureanalysis',
    # 'development',
    # 'heterogeneity',
    # 'communication',
    "functionalmodule",
    # 'multiomics',
]

# Configure root logger only once
logging_debug_opt = False
logging.basicConfig(
    level=logging.DEBUG if logging_debug_opt else logging.INFO,
    format="\n%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[create_logging_handler(logging_debug_opt)],
)
