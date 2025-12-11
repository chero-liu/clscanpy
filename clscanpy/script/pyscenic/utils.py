import os
import subprocess
from clscanpy.config import ROOT_PATH

VALID_SPECIES = {"human", "mouse"}


class Analysis:
    def __init__(
        self,
        input: str,
        outdir: str,
        species: str,
        method: str = "grnboost2",
        database: str = "/data/database/pySCENIC/",
        tfs: str = "allTFs.txt",
        genes_vs_motifs: list = None,
        motifs_tbl: str = "motifs-v10nr_clust-nr-m0.001-o0.0.tbl",
        rank_threshold: int = 5000,
        auc_threshold: float = 0.05,
        nes_threshold: float = 3.0,
        min_orthologous_identity: float = 0.0,
        max_similarity_fdr: float = 0.001,
        chunk_size: int = 200,
        thresholds: str = "0.75 0.90",
        top_n_targets: int = 50,
        min_genes: int = 20,
        all_modules: str = "False",
        num_workers: int = 10,
    ):
        self.input = input
        self.outdir = outdir
        self.species = species
        self.method = method
        self.database = database
        self.tfs = tfs
        self.genes_vs_motifs = genes_vs_motifs
        self.motifs_tbl = motifs_tbl

        self.rank_threshold = rank_threshold
        self.auc_threshold = auc_threshold
        self.nes_threshold = nes_threshold
        self.min_orthologous_identity = min_orthologous_identity
        self.max_similarity_fdr = max_similarity_fdr
        self.chunk_size = chunk_size
        self.thresholds = thresholds
        self.top_n_targets = top_n_targets
        self.min_genes = min_genes
        self.all_modules = all_modules
        self.num_workers = num_workers

        if self.all_modules != "False":
            self.min_genes = f"{self.min_genes}  --all_modules"

    @property
    def shell_script(self):
        genes_vs_motifs_arg = f"{self.genes_vs_motifs[0]}"
        if len(self.genes_vs_motifs) > 1:
            genes_vs_motifs_arg += f" {self.genes_vs_motifs[1]}"

        shell_script_content = f"""
module purge
source /gpfs/oe-scrna/liuchenglong/anaconda3/bin/activate {os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(ROOT_PATH))))}

database={self.database}
input={self.input}
outdir={self.outdir}
species={self.species}
num_workers={self.num_workers}
method={self.method}

# Step 1: GRN inference
pyscenic grn \\
    --num_workers $num_workers \\
    --output $outdir/expr_mat.adjacencies.tsv \\
    --method $method \\
    $input \\
    $database/$species/{self.tfs}

# Step 2: Context enrichment
pyscenic ctx \\
    $outdir/expr_mat.adjacencies.tsv \\
    {genes_vs_motifs_arg} \\
    --annotations_fname $database/$species/{self.motifs_tbl} \\
    --expression_mtx_fname $input \\
    --mode "custom_multiprocessing" \\
    --output $outdir/regulons.csv \\
    --num_workers $num_workers \\
    --rank_threshold {self.rank_threshold}  \\
    --auc_threshold {self.auc_threshold} \\
    --nes_threshold {self.nes_threshold} \\
    --min_orthologous_identity {self.min_orthologous_identity} \\
    --max_similarity_fdr {self.max_similarity_fdr} \\
    --chunk_size {self.chunk_size} \\
    --thresholds {self.thresholds} \\
    --top_n_targets {self.top_n_targets} \\
    --min_genes {self.min_genes}

# Step 3: AUCell scoring
pyscenic aucell \\
    --num_workers $num_workers \\
    --rank_threshold {self.rank_threshold}  \\
    --auc_threshold {self.auc_threshold} \\
    --nes_threshold {self.nes_threshold} \\
    --output $outdir/sce_SCENIC.loom \\
    $input \\
    $outdir/regulons.csv
"""
        return shell_script_content

    def save_script(self, shell_script_content):
        shell_path = f"{os.path.dirname(self.outdir)}/script/analysis.sh"
        with open(
            shell_path,
            "w",
        ) as file:
            file.write(shell_script_content)

        subprocess.run(shell_script_content, shell=True)

    def run(self):
        self.save_script(self.shell_script)


class Visualize:
    def __init__(
        self,
        input: str,
        outdir: str,
        species: str,
        rds_filepath: str,
        groupby: str,
        threshold: int,
        regulons_path: str,
        topGenes: int,
        extended: str,
        nclust: int,
        utils_path: str,
        groupby_levels: str = "no",
    ):
        self.input = input
        self.outdir = outdir
        self.species = species
        self.rds_filepath = rds_filepath
        self.groupby = groupby
        self.threshold = threshold
        self.regulons_path = regulons_path
        self.topGenes = topGenes
        self.extended = extended
        self.nclust = nclust
        self.utils_path = utils_path
        self.groupby_levels = groupby_levels

        if str(self.rds_filepath).endswith(".h5ad"):
            self.h5ad_filepath = self.rds_filepath
            self.supplemental_text = f"""
module purge
module load OESingleCell/3.0.d
Rscript {ROOT_PATH}/script/data_transformation/h5ad2rds.r -i {self.h5ad_filepath}


            """
            self.rds_filepath = os.path.join(
                os.path.dirname(self.rds_filepath), "seurat.h5seurat"
            )
        else:
            self.supplemental_text = ""

    @property
    def shell_script(self):
        shell_script_content = f"""
{self.supplemental_text}
#!/bin/bash
module purge
module load OESingleCell/3.0.d

Rscript {os.path.dirname(self.utils_path)}/visualize.r \\
  --input {self.input} \\
  --rds_filepath {self.rds_filepath} \\
  --outdir {self.outdir} \\
  --subnew_celltype all \\
  --subsampleid all  \\
  --subgroup all  \\
  --subcluster all \\
  --predicate  all \\
  --groupby {self.groupby} \\
  --groupby_levels {self.groupby_levels} \\
  --threshold {self.threshold} \\
  --regulons_path {self.regulons_path} \\
  --topGenes {self.topGenes} \\
  --extended {self.extended} \\
  --nclust  {self.nclust} \\
  --utils_path {self.utils_path}

"""
        return shell_script_content

    def save_script(self, shell_script_content):
        shell_path = f"{os.path.dirname(self.outdir)}/script/visualize.sh"
        with open(
            shell_path,
            "w",
        ) as file:
            file.write(shell_script_content)

        subprocess.run(shell_script_content, shell=True)

    def run(self):
        self.save_script(self.shell_script)
