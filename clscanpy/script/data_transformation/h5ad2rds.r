suppressMessages({
  library(rhdf5)
  library(Seurat)
  library(Matrix)
  library(optparse)
  library(jsonlite)
  library(dplyr)
})

source("/home/chenglong.liu/RaD/myscript/seurat_fc.r")

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="the h5ad file")
)
argv <- parse_args(OptionParser(option_list=option_list))

to_rds <- function(h5ad_file = NULL){
  # 打印调试信息
  print(paste("Processing file:", h5ad_file))
  
  # h5ls(h5ad_file)
  counts_layer <- h5read(file=h5ad_file,name="layers/raw")
  var_genes <- h5read(file=h5ad_file,name="/var")$'_index'
  obs_cells <- h5read(file=h5ad_file,name="/obs")$'_index'
  
  # 确保基因名是唯一的
  var_genes <- make.unique(var_genes)
  
  # 确保细胞名是唯一的
  obs_cells <- make.unique(obs_cells)
  
  print(paste("Number of genes:", length(var_genes)))
  print(paste("Number of cells:", length(obs_cells)))
  
  # Create sparse matrix for counts
  counts_mtx <- try(sparseMatrix(i = as.integer(counts_layer$indices),
                      p = as.integer(counts_layer$indptr),
                      x = as.numeric(counts_layer$data),
                      dims = c(length(var_genes), length(obs_cells)),
                      index1 = F,
                      repr = "C"))
  
  # 如果维度不正确，尝试转置
  if ("try-error" %in% class(counts_mtx)){
    print("Trying transposed orientation...")
    counts_mtx <- sparseMatrix(i = as.integer(counts_layer$indices),
                      p = as.integer(counts_layer$indptr),
                      x = as.numeric(counts_layer$data),
                      dims = c(length(obs_cells), length(var_genes)),
                      index1 = F,
                      repr = "C")
    # 需要转置以符合基因×细胞的标准格式
    counts_mtx <- t(counts_mtx)
  }
  
  # 设置行名和列名
  rownames(counts_mtx) = var_genes
  colnames(counts_mtx) = obs_cells
  
  print(paste("Matrix dimensions:", nrow(counts_mtx), "x", ncol(counts_mtx)))
  
  # 检查是否有重复的基因名
  if(any(duplicated(rownames(counts_mtx)))) {
    warning("Found duplicate gene names. Making them unique...")
    rownames(counts_mtx) <- make.unique(rownames(counts_mtx))
  }
  
  # 检查是否有重复的细胞名
  if(any(duplicated(colnames(counts_mtx)))) {
    warning("Found duplicate cell names. Making them unique...")
    colnames(counts_mtx) <- make.unique(colnames(counts_mtx))
  }
  
  # 检查是否有空基因名
  if(any(rownames(counts_mtx) == "")) {
    warning("Found empty gene names. Replacing with 'Unknown'...")
    rownames(counts_mtx)[rownames(counts_mtx) == ""] <- paste0("Unknown_", seq_len(sum(rownames(counts_mtx) == "")))
  }
  
  # 检查是否有空细胞名
  if(any(colnames(counts_mtx) == "")) {
    warning("Found empty cell names. Replacing with 'Cell'...")
    colnames(counts_mtx)[colnames(counts_mtx) == ""] <- paste0("Cell_", seq_len(sum(colnames(counts_mtx) == "")))
  }

  # 检查矩阵中是否有NaN或Inf值
  if(any(is.nan(counts_mtx@x))) {
    warning("Found NaN values in matrix. Replacing with 0...")
    counts_mtx@x[is.nan(counts_mtx@x)] <- 0
  }
  
  if(any(is.infinite(counts_mtx@x))) {
    warning("Found infinite values in matrix. Replacing with 0...")
    counts_mtx@x[is.infinite(counts_mtx@x)] <- 0
  }

  print("Creating Seurat object...")
  # Create Seurat object with v5 syntax
  # 使用tryCatch来捕获可能的错误
  PRO <- tryCatch({
    CreateSeuratObject(counts = counts_mtx, 
                       min.cells = 0, 
                       min.features = 0)
  }, error = function(e) {
    print(paste("Error creating Seurat object:", e$message))
    print("Trying alternative approach...")
    
    # 尝试使用矩阵的子集或转换
    counts_mtx <- as(counts_mtx, "dgCMatrix")
    CreateSeuratObject(counts = counts_mtx, 
                       min.cells = 0, 
                       min.features = 0)
  })
  
  print(paste("Seurat object created with", ncol(PRO), "cells and", nrow(PRO), "features"))
  
  # Normalize data
  print("Normalizing data...")
  PRO <- NormalizeData(PRO, verbose = FALSE)
  
  # Add meta data from obs
  print("Adding metadata...")
  obs <- h5read(file = h5ad_file, name = "/obs")
  
  for (col_name in names(obs)) {
    if (col_name == '_index') next  # Skip the index column
    
    tryCatch({
      valid_col_name <- make.names(col_name)
      
      col_data <- obs[[col_name]]
      
      if (is.list(col_data) && !is.null(col_data$categories) && !is.null(col_data$codes)) {
        # Handle categorical data
        col_values <- col_data$categories[col_data$codes + 1]
        # 确保长度匹配
        if(length(col_values) == ncol(PRO)) {
          PRO@meta.data[[valid_col_name]] <- factor(col_values, levels = col_data$categories)
          print(paste("Processed column:", col_name, "as categorical factor"))
        } else {
          print(paste("Skipping column", col_name, ": length mismatch"))
        }
      } else if (is.list(col_data) && length(col_data) == ncol(PRO)) {
        # Handle lists that match cell count
        PRO@meta.data[[valid_col_name]] <- unlist(col_data)
        print(paste("Processed column:", col_name, "as unlisted"))
      } else if (is.vector(col_data) && length(col_data) == ncol(PRO)) {
        # Handle vectors
        PRO@meta.data[[valid_col_name]] <- col_data
        print(paste("Processed column:", col_name, "as vector"))
      } else {
        print(paste("Skipping column", col_name, ": unexpected type or length"))
      }
    }, error = function(e) {
      print(paste("Failed to process column:", col_name, "with error:", e$message))
    })
  }

  # Add dimensional reductions
  print("Adding dimensional reductions...")
  
  tryCatch({
    umapinfo <- as.data.frame(t(h5read(file=h5ad_file,name='/obsm/X_umap')))
    if(nrow(umapinfo) == ncol(PRO) && ncol(umapinfo) >= 2) {
      rownames(umapinfo) <- colnames(PRO)
      colnames(umapinfo)[1:2] <- c('UMAP_1','UMAP_2')
      PRO[['umap']] <- CreateDimReducObject(
        embeddings = as.matrix(umapinfo[, 1:2]), 
        key = 'UMAP_', 
        assay = DefaultAssay(PRO)
      )
      print("Added UMAP")
    } else {
      print(paste("UMAP dimension mismatch. Expected", ncol(PRO), "rows, got", nrow(umapinfo)))
    }
  }, error=function(e){
    print(paste("NO UMAP:", e$message))
  })

  tryCatch({
    tsneinfo <- as.data.frame(t(h5read(file=h5ad_file,name='/obsm/X_tsne')))
    if(nrow(tsneinfo) == ncol(PRO) && ncol(tsneinfo) >= 2) {
      rownames(tsneinfo) <- colnames(PRO)
      colnames(tsneinfo)[1:2] <- c('TSNE_1','TSNE_2')
      PRO[['tsne']] <- CreateDimReducObject(
        embeddings = as.matrix(tsneinfo[, 1:2]), 
        key = 'TSNE_', 
        assay = DefaultAssay(PRO)
      )
      print("Added TSNE")
    } else {
      print(paste("TSNE dimension mismatch. Expected", ncol(PRO), "rows, got", nrow(tsneinfo)))
    }
  }, error=function(e){
    print(paste("NO TSNE:", e$message))
  })

  tryCatch({
    pcainfo <- as.data.frame(t(h5read(file=h5ad_file,name='/obsm/X_pca')))
    if(nrow(pcainfo) == ncol(PRO)) {
      rownames(pcainfo) <- colnames(PRO)
      colpc <- paste0("PC_", 1:ncol(pcainfo))
      colnames(pcainfo) <- colpc
      PRO[['pca']] <- CreateDimReducObject(
        embeddings = as.matrix(pcainfo), 
        key = 'PC_', 
        assay = DefaultAssay(PRO)
      )
      print("Added PCA")
    } else {
      print(paste("PCA dimension mismatch. Expected", ncol(PRO), "rows, got", nrow(pcainfo)))
    }
  }, error=function(e){
    print(paste("NO PCA:", e$message))
  })

  tryCatch({
    harmonyinfo <- as.data.frame(t(h5read(file=h5ad_file,name='/obsm/X_pca_harmony')))
    if(nrow(harmonyinfo) == ncol(PRO)) {
      rownames(harmonyinfo) <- colnames(PRO)
      colpc <- paste0("harmony_", 1:ncol(harmonyinfo))
      colnames(harmonyinfo) <- colpc
      PRO[['harmony']] <- CreateDimReducObject(
        embeddings = as.matrix(harmonyinfo), 
        key = 'harmony_', 
        assay = DefaultAssay(PRO)
      )
      print("Added HARMONY")
    } else {
      print(paste("HARMONY dimension mismatch. Expected", ncol(PRO), "rows, got", nrow(harmonyinfo)))
    }
  }, error=function(e){
    print(paste("NO HARMONY:", e$message))
  })

  # Add variable features
  print("Adding variable features...")
  tryCatch({
    highly_variable <- h5read(file=h5ad_file, name="var/highly_variable")
    var_genes <- h5read(file=h5ad_file, name="/var")$'_index'
    # 确保var_genes也是唯一的
    var_genes <- make.unique(var_genes)
    
    # 检查长度是否匹配
    if(length(highly_variable) == length(var_genes)) {
      high_var_genes <- var_genes[as.logical(highly_variable)]
      # 确保高变基因在数据集中存在
      high_var_genes <- high_var_genes[high_var_genes %in% rownames(PRO)]
      if(length(high_var_genes) > 0) {
        VariableFeatures(PRO) <- high_var_genes
        print(paste("Added", length(high_var_genes), "variable features"))
      } else {
        print("No valid variable features found")
      }
    } else {
      print(paste("Length mismatch between highly_variable (", length(highly_variable), 
                  ") and var_genes (", length(var_genes), ")"))
    }
  }, error=function(e) {
    print(paste("NO highly_variable or error:", e$message))
  })

  # Add additional metadata from var
  print("Adding feature metadata...")
  tryCatch({
    var_data <- h5read(file=h5ad_file, name="/var")
    for (col_name in names(var_data)) {
      if (col_name == '_index' || col_name == 'highly_variable') next
      tryCatch({
        valid_col_name <- make.names(col_name)
        col_data <- var_data[[col_name]]
        
        # 确保长度匹配
        if(length(col_data) == nrow(PRO)) {
          # 创建或更新feature metadata
          if(is.null(PRO[["RNA"]]@meta.features)) {
            PRO[["RNA"]]@meta.features <- data.frame(row.names = rownames(PRO))
          }
          PRO[["RNA"]]@meta.features[[valid_col_name]] <- col_data
          print(paste("Added feature metadata:", col_name))
        } else {
          print(paste("Skipping feature metadata", col_name, 
                      ": length mismatch, expected", nrow(PRO), "got", length(col_data)))
        }
      }, error=function(e) {
        print(paste("Failed to add feature metadata:", col_name, e$message))
      })
    }
  }, error=function(e) {
    print("No additional var metadata")
  })

  # Add preprocess parameters
  print("Adding preprocess parameters...")
  tryCatch({
    preprocess_para = h5read(file=h5ad_file,name="uns/preprocess_para")
    if(is.raw(preprocess_para)) {
      preprocess_para <- rawToChar(preprocess_para)
    }
    param_list <- fromJSON(preprocess_para)
    PRO@misc$preprocess_para = param_list
    print("Added preprocess parameters")
  }, error=function(e){
    print(paste("NO preprocess_para:", e$message))
  })
  
  # Try to set clusters as active identities
  print("Setting active identities...")
  tryCatch({
    if ("clusters" %in% colnames(PRO@meta.data)) {
      Idents(PRO) <- PRO$clusters
      print("Set clusters as active identities")
    }
  }, error=function(e){
    print("NO RAW CLUSTER or cannot set as identity")
  })
  
  print("Conversion complete!")
  return(PRO)
}

input <- argv$input
oudat <- gsub('.h5ad$', '.rds', input)
print(paste("Input file:", input))
print(paste("Output file:", oudat))

PRO <- to_rds(h5ad_file = input)

# Save the Seurat v5 object
saveRDS(PRO, file = oudat)
print(paste("Successfully saved Seurat v5 object to:", oudat))
