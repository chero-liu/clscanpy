import unittest
from clscanpy.log import log_function_call
from clscanpy.tools.utils import check_mkdir
import argparse
import os
import sys
import subprocess
import tempfile
import json
import random
from datetime import datetime
from typing import Dict, List, Optional

import logging

LOGGER = logging.getLogger(__name__)


class Enrich:

    def __init__(
        self,
        input: str,
        outdir: str,
        bp: str,
        cc: str,
        mf: str,
        kegg: str,
        clusters: str = "",
        top_n: int = 20,
        pvalue: float = 0.05,
        gene_column: str = "names",
        cluster_column: str = "clusters",
        min_genes: int = 3,
        max_genes: int = 500,
        truncate_length: int = 30,
        verbose: bool = False,
        skip_visualization: bool = False,
        skip_report: bool = False,
    ):
        self.input = input
        self.outdir = outdir
        self.bp = bp
        self.cc = cc
        self.mf = mf
        self.kegg = kegg
        self.clusters = clusters
        self.top_n = top_n
        self.pvalue = pvalue
        self.gene_column = gene_column
        self.cluster_column = cluster_column
        self.min_genes = min_genes
        self.max_genes = max_genes
        self.truncate_length = truncate_length
        self.verbose = verbose
        self.skip_visualization = skip_visualization
        self.skip_report = skip_report

    def run_command(self, cmd: List[str], description: str) -> tuple[bool, str]:
        """Run subcommand"""
        if self.verbose:
            print(f"\n{'='*60}")
            print(f"Executing: {description}")
            print(f"Command: {' '.join(cmd)}")
            print(f"{'='*60}")

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            if self.verbose:
                print(f"Output:\n{result.stdout}")
                if result.stderr:
                    print(f"Error output:\n{result.stderr}")
            return True, result.stdout
        except subprocess.CalledProcessError as e:
            print(f"Error: {description} failed")
            print(f"Return code: {e.returncode}")
            print(f"Output:\n{e.stdout}")
            print(f"Error output:\n{e.stderr}")
            return False, e.stderr

    def run_analysis_module(self) -> Optional[str]:
        """Run analysis module"""
        print(f"\n{'='*60}")
        print("Step 1/3: Running analysis module")
        print(f"{'='*60}")

        # Build analysis command
        analysis_cmd = [
            sys.executable,
            os.path.join(os.path.dirname(__file__), "enrich", "analysis.py"),
            "--input",
            self.input,
            "--outdir",
            self.outdir,
            "--bp",
            self.bp,
            "--cc",
            self.cc,
            "--mf",
            self.mf,
            "--kegg",
            self.kegg,
            "--pvalue",
            str(self.pvalue),
            "--gene-column",
            self.gene_column,
            "--cluster-column",
            self.cluster_column,
            "--min-genes",
            str(self.min_genes),
            "--max-genes",
            str(self.max_genes),
            "--truncate-length",
            str(self.truncate_length),
        ]

        if self.clusters:
            analysis_cmd.extend(["--clusters", self.clusters])

        if self.verbose:
            analysis_cmd.append("--verbose")

        success, output = self.run_command(analysis_cmd, "Analysis module")

        if success:
            summary_file = os.path.join(self.outdir, "summary.json")
            if os.path.exists(summary_file):
                print(f"✓ Analysis completed! Results saved in: {self.outdir}")
                print(f"  Analysis summary: {summary_file}")
                return summary_file
            else:
                print(
                    f"✗ Analysis completed but summary file not found: {summary_file}"
                )
                return None
        else:
            print("✗ Analysis module execution failed")
            return None

    def run_visualization_module(self, summary_file: str) -> bool:
        """Run visualization module"""
        if self.skip_visualization:
            print(f"\n{'='*60}")
            print(
                "Step 2/3: Skipping visualization module (user specified --skip-visualization)"
            )
            print(f"{'='*60}")
            return True

        print(f"\n{'='*60}")
        print("Step 2/3: Running visualization module")
        print(f"{'='*60}")

        # Build visualization command
        visualization_cmd = [
            sys.executable,
            os.path.join(os.path.dirname(__file__), "enrich", "visualization.py"),
            "--input-dir",
            self.outdir,
            "--output-dir",
            self.outdir,
            "--summary-file",
            summary_file,
            "--top-n",
            str(self.top_n),
            "--pvalue",
            str(self.pvalue),
        ]

        if self.verbose:
            visualization_cmd.append("--verbose")

        success, output = self.run_command(visualization_cmd, "Visualization module")

        if success:
            print(f"✓ Visualization completed! Charts saved in: {self.outdir}")
            return True
        else:
            print("✗ Visualization module execution failed")
            return False

    def run_report_module(self, summary_file: str) -> bool:
        """Run report module"""
        if self.skip_report:
            print(f"\n{'='*60}")
            print("Step 3/3: Skipping report module (user specified --skip-report)")
            print(f"{'='*60}")
            return True

        print(f"\n{'='*60}")
        print("Step 3/3: Running report module")
        print(f"{'='*60}")

        # Build report command
        report_cmd = [
            sys.executable,
            os.path.join(os.path.dirname(__file__), "enrich", "report.py"),
            "--input-dir",
            self.outdir,
            "--output-dir",
            self.outdir,
            "--summary-file",
            summary_file,
            "--bp",
            self.bp,
            "--cc",
            self.cc,
            "--mf",
            self.mf,
            "--kegg",
            self.kegg,
            "--pvalue",
            str(self.pvalue),
        ]

        if self.verbose:
            report_cmd.append("--verbose")

        success, output = self.run_command(report_cmd, "Report module")

        if success:
            report_file = os.path.join(self.outdir, "report.txt")
            if os.path.exists(report_file):
                print(f"✓ Report generation completed! Saved as: {report_file}")
                return True
            else:
                print(f"✗ Report generation completed but report file not found")
                return False
        else:
            print("✗ Report module execution failed")
            return False

    def run(self):
        LOGGER.info(f"Start enrich ...")
        check_mkdir(self.outdir)
        random.seed(12345)

        print("=" * 80)
        print("Enrichment Analysis Main Workflow")
        print("=" * 80)
        print(f"Input file: {os.path.basename(self.input)}")
        print(f"Output directory: {self.outdir}")
        print(
            f"Analysis parameters: p-value={self.pvalue}, show top {self.top_n} results"
        )
        if self.skip_visualization:
            print(f"Skipping: Visualization step")
        if self.skip_report:
            print(f"Skipping: Report generation step")
        print("=" * 80)

        # Create output directory structure - flat structure
        tables_dir = os.path.join(self.outdir, "tables")
        figures_dir = os.path.join(self.outdir, "figures")

        os.makedirs(self.outdir, exist_ok=True)
        os.makedirs(tables_dir, exist_ok=True)
        os.makedirs(figures_dir, exist_ok=True)

        # Save analysis parameters to arguments.txt
        arguments_file = os.path.join(self.outdir, "arguments.txt")
        with open(arguments_file, "w") as f:
            f.write("clusters={}\n".format(self.clusters))
            f.write("top_n={}\n".format(self.top_n))
            f.write("pvalue={}\n".format(self.pvalue))
            f.write("gene_column={}\n".format(self.gene_column))
            f.write("cluster_column={}\n".format(self.cluster_column))
            f.write("min_genes={}\n".format(self.min_genes))
            f.write("max_genes={}\n".format(self.max_genes))
            f.write("truncate_length={}\n".format(self.truncate_length))
            f.write("input={}\n".format(os.path.basename(self.input)))
            f.write("bp={}\n".format(os.path.basename(self.bp)))
            f.write("cc={}\n".format(os.path.basename(self.cc)))
            f.write("mf={}\n".format(os.path.basename(self.mf)))
            f.write("kegg={}\n".format(os.path.basename(self.kegg)))

        if self.verbose:
            print(f"✓ Analysis parameters saved to: {arguments_file}")

        # Step 1: Run analysis module
        summary_file = self.run_analysis_module()
        if not summary_file:
            print("✗ Analysis module failed, terminating workflow")
            return

        # Step 2: Run visualization module
        visualization_success = self.run_visualization_module(summary_file)
        if not visualization_success and not self.skip_visualization:
            print("⚠ Visualization module failed, but continuing with report module")

        # Step 3: Run report module
        report_success = self.run_report_module(summary_file)
        if not report_success and not self.skip_report:
            print("⚠ Report module failed")

        # Generate final summary
        print(f"\n{'='*80}")
        print("Enrichment Analysis Workflow Completed")
        print(f"{'='*80}")
        print(f"\nResult files:")
        print(f"  Output directory: {self.outdir}")
        print(f"  1. Data tables: {tables_dir}")
        print(f"  2. Analysis summary: {summary_file}")

        if not self.skip_visualization:
            print(f"  3. Visualization charts: {figures_dir}")

        if not self.skip_report:
            report_file = os.path.join(self.outdir, "report.txt")
            print(f"  4. Analysis report: {report_file}")

        print(f"\n{'='*80}")
        print("Tip: You can run individual modules using the following commands:")
        print(f"  Analysis module: python enrich/analysis.py [parameters]")
        print(f"  Visualization module: python enrich/visualization.py [parameters]")
        print(f"  Report module: python enrich/report.py [parameters]")
        print(f"{'='*80}")

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        LOGGER.info("enrich done")


@log_function_call
def enrich(args):
    with Enrich(
        input=args.input,
        outdir=args.outdir,
        bp=args.bp,
        cc=args.cc,
        mf=args.mf,
        kegg=args.kegg,
        clusters=args.clusters,
        top_n=args.top_n,
        pvalue=args.pvalue,
        gene_column=args.gene_column,
        cluster_column=args.cluster_column,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        truncate_length=args.truncate_length,
        verbose=args.verbose,
        skip_visualization=args.skip_visualization,
        skip_report=args.skip_report,
    ) as runner:
        runner.run()


def get_opts_enrich(parser, sub_program=False):
    """Add command line arguments for enrichment analysis"""
    if sub_program:
        # Required parameters
        parser.add_argument(
            "-i",
            "--input",
            type=str,
            required=True,
            help="Input file path (TSV format)",
        )
        parser.add_argument(
            "-o", "--outdir", type=str, required=True, help="Output directory path"
        )
        parser.add_argument(
            "--bp", type=str, required=True, help="GO BP gene set file path"
        )
        parser.add_argument(
            "--cc", type=str, required=True, help="GO CC gene set file path"
        )
        parser.add_argument(
            "--mf", type=str, required=True, help="GO MF gene set file path"
        )
        parser.add_argument(
            "--kegg", type=str, required=True, help="KEGG gene set file path"
        )

        # Optional parameters
        parser.add_argument(
            "--clusters",
            type=str,
            default="",
            help="Specify cluster list to analyze, comma-separated, e.g., '0,1,2' (default: analyze all clusters)",
        )
        parser.add_argument(
            "--top-n",
            type=int,
            default=20,
            help="Top N enrichment results to display for each category (default: 20)",
        )
        parser.add_argument(
            "--pvalue",
            type=float,
            default=0.05,
            help="p-value threshold (default: 0.05)",
        )
        parser.add_argument(
            "--gene-column",
            type=str,
            default="names",
            help="Gene column name (default: 'names')",
        )
        parser.add_argument(
            "--cluster-column",
            type=str,
            default="clusters",
            help="Cluster column name (default: 'clusters')",
        )
        parser.add_argument(
            "--min-genes",
            type=int,
            default=3,
            help="Minimum gene count for enrichment analysis (default: 3)",
        )
        parser.add_argument(
            "--max-genes",
            type=int,
            default=500,
            help="Maximum gene count for enrichment analysis (default: 500)",
        )
        parser.add_argument(
            "--truncate-length",
            type=int,
            default=30,
            help="Term truncation length (default: 30)",
        )
        parser.add_argument(
            "--verbose", action="store_true", help="Show detailed output information"
        )
        parser.add_argument(
            "--skip-visualization", action="store_true", help="Skip visualization step"
        )
        parser.add_argument(
            "--skip-report", action="store_true", help="Skip report generation step"
        )
    else:
        # When not called as a subprogram, only basic parameters are needed
        parser.add_argument(
            "-i",
            "--input",
            type=str,
            required=True,
            help="Input file path (TSV format)",
        )
        parser.add_argument(
            "-o", "--outdir", type=str, required=True, help="Output directory path"
        )

    return parser


if __name__ == "__main__":
    unittest.main()
