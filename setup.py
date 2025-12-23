# setup.py

from setuptools import setup, find_packages

from clscanpy.__init__ import __VERSION__

setup(
    name="clscanpy",
    version=__VERSION__,
    description="A Python package for single-cell RNA-seq data analysis",
    author="liuchenglong",
    author_email="njlcl@outlook.com",
    install_requires=[
        "scanpy>=1.9.0",
        "anndata>=0.10.0",
        "pandas>=1.5.0",
        "numpy>=1.23.0",
        "matplotlib>=3.7.0",
        "scipy>=1.10.0",
        "doubletdetection>=4.2.0",
        "celltypist>=1.6.0",
        "click>=8.1.0",
    ],
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "clscanpy = clscanpy.clscanpy:main",
        ],
    },
    include_package_data=True,
)
