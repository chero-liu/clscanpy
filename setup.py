# setup.py

from setuptools import setup, find_packages

from clscanpy.__init__ import __VERSION__

setup(
    name="clscanpy",
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
            "clscanpy = clscanpy.clscanpy:main",
        ],
    },
    include_package_data=True,
)
