suppressMessages({
library(rhdf5)
library(Seurat)
library(Matrix)
library(optparse)
library(jsonlite)
})

source("/home/liuchenglong/script/lclFunc.r")

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="the h5ad file")
)
argv <- parse_args(OptionParser(option_list=option_list))

to_rds <- function(h5ad_file = NULL){
  # h5ls(h5ad_file)
  counts_layer <- h5read(file=h5ad_file,name="layers/raw")
  # normalised_layer <- h5read(file=h5ad_file,name="layers/normalised")
  var_genes <- h5read(file=h5ad_file,name="/var")$'_index'
  obs_cells <- h5read(file=h5ad_file,name="/obs")$'_index'

  counts_mtx <- try(sparseMatrix(i = as.integer(counts_layer$indices),
                      p = as.integer(counts_layer$indptr),
                      x = as.numeric(counts_layer$data),
                      dims = c(length(var_genes), length(obs_cells)),
                      index1 = F,
                      repr = "C"))
  if ("try-error" %in% class(counts_mtx)){
    counts_mtx <- t(sparseMatrix(i = as.integer(counts_layer$indices),
                      p = as.integer(counts_layer$indptr),
                      x = as.numeric(counts_layer$data),
                      dims = c(length(obs_cells), length(var_genes)),
                      index1 = F,
                      repr = "C"))
  }
  rownames(counts_mtx) = var_genes
  colnames(counts_mtx) = obs_cells

  PRO <- CreateSeuratObject(counts = counts_mtx, min.cells = 0)
  PRO <- NormalizeData(PRO)

  obs <- h5read(file = h5ad_file, name = "/obs")
  
  for (col_name in names(obs)) {
    tryCatch({
      valid_col_name <- make.names(col_name)
      
      col_data <- obs[[col_name]]
      
      if ("array" %in% class(col_data)) {
        PRO[[valid_col_name]] <- col_data
        print(paste("Processed column:", col_name, "as array added directly"))
      } else if ("list" %in% class(col_data)) {
        col_values <- col_data$categories[col_data$codes + 1]
        PRO[[valid_col_name]] <- factor(col_values, levels = col_data$categories)
        print(paste("Processed column:", col_name, "as categorical factor"))
      } else {
        warning(paste("Unsupported data type in column:", col_name, "of class:", class(col_data)))
      }
    }, error = function(e) {
      print(paste("Failed to process column:", col_name, "with error:", e$message))
    })
  }


  tryCatch({ umapinfo <- as.data.frame(t(h5read(file=h5ad_file,name='/obsm/X_umap')))
    rownames(umapinfo) <- colnames(PRO)
    colnames(umapinfo) <- c('UMAP_1','UMAP_2')
    PRO[['umap']] <- CreateDimReducObject(embeddings = as.matrix(umapinfo), key = 'UMAP_', assay = 'RNA', global = TRUE)}
  ,error=function(e){print("NO UMAP")})

  tryCatch({tsneinfo <- as.data.frame(t(h5read(file=h5ad_file,name='/obsm/X_tsne')))
           rownames(tsneinfo) <- colnames(PRO)
           colnames(tsneinfo) <- c('TSNE_1','TSNE_2')
           PRO[['tsne']] <- CreateDimReducObject(embeddings = as.matrix(tsneinfo), key = 'TSNE_', assay = 'RNA', global = TRUE)}
  ,error=function(e){print("NO TSNE")})

  tryCatch({pcainfo <- as.data.frame(t(h5read(file=h5ad_file,name='/obsm/X_pca')))
  rownames(pcainfo) <- colnames(PRO)
  colpc <- c()
  for (i in 1:ncol(pcainfo)){
    colpc <- c(colpc, paste0("PCA_",i))
  }
  colnames(pcainfo) <- colpc
  PRO[['pca']] <- CreateDimReducObject(embeddings = as.matrix(pcainfo), key = 'PCA_', assay = 'RNA', global = TRUE)}
  ,error=function(e){print("NO PCA")})

  tryCatch({harmonyinfo <- as.data.frame(t(h5read(file=h5ad_file,name='/obsm/X_pca_harmony')))
  rownames(harmonyinfo) <- colnames(PRO)
  colpc <- c()
  for (i in 1:ncol(harmonyinfo)){
    colpc <- c(colpc, paste0("HARMONY_",i))
  }
  colnames(harmonyinfo) <- colpc
  PRO[['harmony']] <- CreateDimReducObject(embeddings = as.matrix(harmonyinfo), key = 'HARMONY_', assay = 'RNA', global = TRUE)}
  ,error=function(e){print("NO HARMONY")})

  tryCatch({
    highly_variable <- h5read(file=h5ad_file, name="var/highly_variable")
    var_genes <- h5read(file=h5ad_file, name="/var")$'_index'
    high_var_genes <- var_genes[as.logical(highly_variable)]
    PRO@assays$RNA@var.features <- high_var_genes
  }, error=function(e) {
    print("NO highly_variable")
  })

  tryCatch({
  preprocess_para = h5read(file=h5ad_file,name="uns/preprocess_para")
  param_list <- fromJSON(preprocess_para)
  PRO@misc$preprocess_para = param_list
  }
  ,error=function(e){print("NO preprocess_para")})
  return(PRO)
}

input <- argv$input
# oudat <- gsub('.h5ad$', '', input)
oudat = dirname(input)
PRO <- to_rds(h5ad_file = input)
tryCatch({Idents(PRO) <- PRO$clusters},error=function(e){print("NO RAW CLUSTER")})

saveH5seurat(PRO,outdir = oudat)
