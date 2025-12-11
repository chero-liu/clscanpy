#!/usr/bin/env Rscript

library(Seurat)
library(stringr)
library(tools)
library(optparse)
source("/gpfs/oe-scrna/pipeline/scRNA-seq_further_analysis/function/seuratFc.r")
library(RcisTarget,lib.loc = '/gpfs/oe-software/scrna_software_bak/conda_env/Scenic/lib/R/library')

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input file path (RDS file)"),
  make_option(c("--species"), type = "character", help = "species (e.g., 'human', 'mouse')"),
  make_option(c("-o", "--outdir"), type = "character", default = "./", help = "Output directory to save the CSV file"),
  make_option(c("--groupby"), type = "character", default = "clusters", help = "Group by column"),
  make_option(c("-d", "--downsample"), default = 'no', help = "Downsample number"),
  make_option(c("-c", "--subnew_celltype"), type = "character", default = "all", help = "Subnew cell type"),
  make_option(c("-s", "--subsampleid"), type = "character", default = "all", help = "Subsample ID"),
  make_option(c("-g", "--subgroup"), type = "character", default = "all", help = "Subgroup"),
  make_option(c("-p", "--predicate"), type = "character", default = 'all', help = "Filtering predicate (e.g., '')"),
  make_option(c("-l", "--subcluster"), type = "character", default = "all", help = "Subcluster"),
  make_option(c("--utils_path"), type = "character", help = "utils_path"),
  make_option(c("--hvg"), type = "character", default = "auto",help = "e.g. auto,no,yes")
)

# Parse the command line arguments
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

source(args$utils_path)

rds = loadRDS(rds_filepath = args$input,
              subnew_celltype = args$subnew_celltype,
              subsampleid = args$subsampleid,
              subgroup = args$subgroup,
              subcluster = args$subcluster,
              predicate = args$predicate,
              downsample = args$downsample,
              groupby = args$groupby)

print(paste0("Number of cells: ", ncol(rds)))

if (args$species == 'human') {
  dbFilePath = '/data/database/SCENIC/human/hg19-500bp-upstream-7species.mc9nr.feather'
} else if (args$species == 'mouse'){
  dbFilePath = '/data/database/SCENIC/mouse/mm9-500bp-upstream-7species.mc9nr.feather'
}else{
    stop("NO motif annotation for your specified species.")
}

hvg = args$hvg

if(hvg == 'auto'){
  if (ncol(rds) > 70000){
    hvg = 'yes'
  }else{
    hvg = 'no'
  }
}

if (hvg == 'yes'){
    genesKept = VariableFeatures(rds)
}else{
    Rcpp::sourceCpp(code= '
    #include <Rcpp.h>
    using namespace Rcpp;
    // [[Rcpp::export]]
    NumericMatrix asMatrix(NumericVector rp,
                          NumericVector cp,
                          NumericVector z,
                          int nrows,
                          int ncols){
      R_xlen_t k = z.size();
      NumericMatrix mat(nrows, ncols);
      for (R_xlen_t i = 0; i < k; i++){
        mat(rp[i],cp[i]) = z[i];
      }
      return mat;
    }
        ' )
    as_matrix = function(mat){
        row_pos = mat@i
        col_pos = findInterval(seq(mat@x)-1,mat@p[-1])
        tmp = asMatrix(rp = row_pos,cp = col_pos,z = mat@x , nrows = mat@Dim[1],ncols = mat@Dim[2])
        row.names(tmp) = mat@Dimnames[[1]]
        colnames(tmp) = mat@Dimnames[[2]]
        return(tmp)
    }
    exprMat = GetAssayData(rds, slot = "counts")
    genesKept <- geneFiltering(as_matrix(exprMat), dbFilePath=dbFilePath, minCountsPerGene=1)
}

rds=rds[rownames(rds) %in% genesKept,]

print(paste0("Number of genes kept: ", nrow(rds)))

SeuratDisk::as.loom(rds, filename = paste0(args$outdir,"/for_scenic.loom"), overwrite = TRUE)

