# setup.py

from setuptools import setup, find_packages

from oescanpy.__init__ import __VERSION__

setup(
    name="oescanpy",
    version=__VERSION__,
    description="A Python package for single-cell RNA-seq data analysis",
    author="liuchenglong",
    author_email="chenglong.liu@oebiotech.com",
    install_requires=[
        # "argparse",
    ],
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "oescanpy = oescanpy.oescanpy:main",
        ],
    },
    include_package_data=True,
)
