source("/gpfs/oe-scrna/pipeline/scRNA-seq_further_analysis/function/seuratFc.r")
source("/gpfs/oe-scrna/guokaiqi/test/report_plot/oe_theme.R")
load_data <- function(args) {
  print(args$combo)
  seurat_obj = loadRDS(rds_filepath = args$rds_filepath,
                        subnew_celltype = args$subnew_celltype,
                        subsampleid = args$subsampleid,
                        subgroup = args$subgroup,
                        subcluster = args$subcluster,
                        predicate = args$predicate,
                        groupby_levels = args$groupby_levels,
                        groupby = args$groupby,
                        combo = args$combo)
  print(colnames(seurat_obj@meta.data))
  scenic_data <- open_loom(args$input)
  # 数据提取
  regulons_incidMat <- get_regulons(scenic_data, column.attr.name = 'Regulons')
  #   regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(scenic_data, column.attr.name = 'RegulonsAUC')
  regulonAUC <- as.data.frame(regulonAUC)

  # 修改regulonAUC行名,并取出前40个TF
  regulons = read.csv(args$regulons_path,sep = '\t')
  regulons_sorted <- regulons[order(-regulons$NES), ]
  top_TF <- unique(regulons_sorted$TF)[1:40]
  tf_counts <- setNames(regulons$TargetGenes_counts, regulons$TF)
  rownames_updated <- sapply(rownames(regulonAUC), function(x) {
    gene <- sub("\\(.*\\)", "", x)  # 提取基因名
    count <- tf_counts[gene]       # 获取对应的 TargetGenes_counts
    if (!is.na(count)) {
      paste0(gene, "(", count, "g)")  # 替换括号内容
    } else {
      x
    }
  })
  rownames(regulonAUC) <- rownames_updated
  regulonAUC <- regulonAUC[rownames(regulonAUC) %in% paste0(top_TF, "(", tf_counts[top_TF], "g)"), ]

  regulonAUC = as.matrix(regulonAUC)
  # 数据对齐
  seurat_obj <- seurat_obj[, colnames(regulonAUC)]
  regulonAUC <- regulonAUC[, rownames(seurat_obj@meta.data)]
  regulonAUC <- regulonAUC[rowSums(regulonAUC) != 0, ]

  list(
    seurat = loadRDS(seurat_obj),
    regulonAUC = regulonAUC
  )
}

#RemoveOutlier 函数的作用是基于中位数和MAD（Median Absolute Deviation，中位数绝对偏差）来检测并移除异常值，支持 按批次处理数据，并提供不同的异常值筛选模式（上限、下限、或两者）
RemoveOutlier <- function(
  metric,
  nmads = 5,
  type = c("both", "lower", "higher"),
  log = FALSE,
  subset = NULL,
  batch = NULL,
  min_diff = NA
) {
    if (log) {
        metric <- log10(metric)
    }
    if (any(is.na(metric))) {
        warning("missing values ignored during outlier detection")
    }

    if (!is.null(batch)) {
        N <- length(metric)
        if (length(batch) != N) {
            stop("length of 'batch' must equal length of 'metric'")
        }

        # Coercing non-NULL subset into a logical vector.
        if (!is.null(subset)) {
            new.subset <- logical(N)
            names(new.subset) <- names(metric)
            new.subset[subset] <- TRUE
            subset <- new.subset
        }

        # Computing QC metrics for each batch.
        by.batch <- split(seq_len(N), batch)
        collected <- logical(N)
        all.threshold <- vector("list", length(by.batch))
        for (b in seq_along(by.batch)) {
            bdx <- by.batch[[b]]
            current <- Recall(metric[bdx], nmads = nmads, type = type, log = FALSE, subset = subset[bdx], batch = NULL, min_diff = min_diff)
            all.threshold[[b]] <- attr(x, "thresholds")
            collected[bdx] <- current
        }

        all.threshold <- do.call(cbind, all.threshold)
        colnames(all.threshold) <- names(by.batch)
        # return(.store_thresholds(collected, all.threshold, logged=log))
        if ( log ){ val <- 10^all.threshold }
        attr(collected, "thresholds") <- val
        return( collected )
    }
    # Computing median/MAD (possibly based on subset of the data).
    if (!is.null(subset)) {
        submetric <- metric[subset]
        if (length(submetric) == 0L) {
            warning("no observations remaining after subsetting")
        }
    } else {
        submetric <- metric
    }
    cur.med <- median(submetric, na.rm = TRUE)
    cur.mad <- mad(submetric, center = cur.med, na.rm = TRUE)

    diff.val <- max(min_diff, nmads * cur.mad, na.rm = TRUE)
    upper.limit <- cur.med + diff.val
    lower.limit <- cur.med - diff.val

    type <- match.arg(type)
    if (type == "lower") {
        upper.limit <- Inf
    } else if (type == "higher") {
        lower.limit <- -Inf
    }

    kx = metric < lower.limit | upper.limit < metric
    val = c(lower=lower.limit, higher=upper.limit)
    if ( log ){
      val <- 10^val
    }
    attr(kx, "thresholds") <- val
    return( kx )
}


geneFiltering <- function (exprMat, dbFilePath, minCountsPerGene = 3 * 0.01 * 
    ncol(exprMat), minSamples = ncol(exprMat) * 0.01) 
{
    if (is.null(dbFilePath)) 
        stop("dbFilePath")
    if (is.data.frame(exprMat)) {
        supportedClasses <- paste(gsub("AUCell_buildRankings,", 
            "", methods("AUCell_buildRankings")), collapse = ", ")
        supportedClasses <- gsub("-method", "", supportedClasses)
        stop("'exprMat' should be one of the following classes: ", 
            supportedClasses, "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
    }
    if (any(table(rownames(exprMat)) > 1)) 
        stop("The rownames (gene id/name) in the expression matrix should be unique.")
    nCountsPerGene <- rowSums(exprMat, na.rm = T)
    nCellsPerGene <- rowSums(exprMat > 0, na.rm = T)
    message("Maximum value in the expression matrix: ", max(exprMat, 
        na.rm = T))
    message("Ratio of detected vs non-detected: ", signif(sum(exprMat > 
        0, na.rm = T)/sum(exprMat == 0, na.rm = T), 2))
    message("Number of counts (in the dataset units) per gene:")
    print(summary(nCountsPerGene))
    message("Number of cells in which each gene is detected:")
    print(summary(nCellsPerGene))
    message("\nNumber of genes left after applying the following filters (sequential):")
    genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > 
        minCountsPerGene)]
    message("\t", length(genesLeft_minReads), "\tgenes with counts per gene > ", 
        minCountsPerGene)
    nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
    genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > 
        minSamples)]
    message("\t", length(genesLeft_minCells), "\tgenes detected in more than ", 
        minSamples, " cells")
    library(RcisTarget)
    motifRankings <- importRankings(dbFilePath)
    genesInDatabase <- colnames(getRanking(motifRankings))
    genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% 
        genesInDatabase)]
    message("\t", length(genesLeft_minCells_inDatabases), "\tgenes available in RcisTarget database")
    genesKept <- genesLeft_minCells_inDatabases

    return(genesKept)
}

process_regulon_activity <- function(cellInfo, regulonAUC_mat, groupby, nmads = 3, scale_by_group_count = TRUE) {
  group_sizes <- table(cellInfo[[groupby]])
  single_cell_groups <- names(group_sizes[group_sizes == 1])

  # Function to calculate mean after outlier removal
  calculate_mean_no_outliers <- function(values) {
    mean(values[!RemoveOutlier(values, nmads = nmads, type = "higher")])
  }

  if (length(single_cell_groups) == 0) {
    # No single-cell groups
    regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo[[groupby]]), function(cells) {
      apply(regulonAUC_mat[, cells, drop = FALSE], 1, calculate_mean_no_outliers)
    })
  } else {
    # Exclude single-cell groups temporarily
    cellInfo_filtered <- cellInfo[!cellInfo[[groupby]] %in% single_cell_groups, ]
    cellInfo_filtered[[groupby]] <- droplevels(cellInfo_filtered[[groupby]])

    regulonActivity_byCellType0 <- sapply(split(rownames(cellInfo_filtered), cellInfo_filtered[[groupby]]), function(cells) {
      apply(regulonAUC_mat[, cells, drop = FALSE], 1, calculate_mean_no_outliers)
    })

    # Add back single-cell groups
    single_cell_matrix <- regulonAUC_mat[, rownames(cellInfo[cellInfo[[groupby]] %in% single_cell_groups, ]), drop = FALSE]
    regulonActivity_byCellType <- cbind(regulonActivity_byCellType0, single_cell_matrix)
    colnames(regulonActivity_byCellType) <- c(colnames(regulonActivity_byCellType0), single_cell_groups)
  }

  # Scale and process the activity matrix
  scaled_data <- t(scale(t(regulonActivity_byCellType), center = TRUE, scale = ifelse(scale_by_group_count && length(unique(cellInfo[[groupby]])) == 2, FALSE, TRUE)))
  regulonActivity_byCellType_processed <- na.omit(scaled_data)

  result <- regulonActivity_byCellType_processed[which(rowSums(abs(regulonActivity_byCellType_processed)) > max(abs(regulonActivity_byCellType_processed))/4),]
  return(result)
}



RSSRanking <- function(
  rss_df,
  group.by,
  ggrepel_force = 2,
  ggrepel_point_padding = 0.5,
  top_genes = 3,
  plot_extended = FALSE,
  size = 4
){
  require(ggrepel)
  require(cowplot)
  require(dplyr)
  require(ggplot2)

  # 筛选是否显示 extended regulons
  if(plot_extended){
    rss_df <- rss_df %>% subset(grepl("extended", regulon))
  } else {
    rss_df <- rss_df %>% subset(!grepl("extended", regulon))
  }

  # 排名
  rss_df_sub <- rss_df %>%
    dplyr::group_by(.data[[group.by]]) %>%
    mutate(rank = order(order(RSS, decreasing = TRUE)))

  # 绘图
  rrs_ranking_plot <- ggplot(rss_df_sub, aes(rank, RSS, label = regulon)) +
    geom_point(color = "grey20", size = 2) +
    geom_point(data = subset(rss_df_sub, rank <= top_genes), color = "red", size = 2) +
    geom_text_repel(
      data = subset(rss_df_sub, rank <= top_genes),
      force = ggrepel_force,
      point.padding = ggrepel_point_padding,
      size = size # 标签字体大小，可调
    ) +
    theme_bw(base_size = 14) +
    theme(
     panel.border = element_rect(color = "black", size = 1.5),
      panel.grid = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),       # facet标题加粗
      axis.text = element_text(face = "bold", size = 12),        # 坐标轴刻度加粗
      axis.title = element_text(face = "bold", size = 12),       # 坐标轴标题加粗
      panel.spacing = unit(1.5, "lines"),
      plot.margin = margin(10, 40, 10, 10)
    ) +
    labs(x = "Rank", y = "RSS", title = group.by) +
    facet_wrap(as.formula(paste("~", group.by)), ncol = 2, scales = "free_y")+
    scale_y_continuous(
      expand = expansion(mult = c(0.1, 0.2))  # 增加顶部留白（20%）
    )
    rrs_ranking_plot = rrs_ranking_plot + oe_theme(base_size = 12,base_family = "sans")

  return(rrs_ranking_plot)
}


calculate_csi <- function(
  regulonAUC = regulonAUC,
  calc_extended = FALSE,
  verbose = FALSE
){
  if (calc_extended == FALSE){
    regulonAUC <- subset(regulonAUC,!grepl("extended",rownames(regulonAUC)))
  }

  pearson_cor <- cor(t(regulonAUC))
  pearson_cor_df <- as.data.frame(pearson_cor)
  pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
  pearson_cor_long <- pearson_cor_df %>%
    tidyr::gather(regulon_2,pcc,-regulon_1) %>%
    dplyr::mutate("regulon_pair" = paste(regulon_1,regulon_2,sep="_"))

  regulon_names <- unique(colnames(pearson_cor))
  num_of_calculations <- length(regulon_names)*length(regulon_names)
  csi_regulons <- data.frame(matrix(nrow=num_of_calculations,ncol = 3))
  colnames(csi_regulons) <- c("regulon_1", "regulon_2", "CSI")
  num_regulons <- length(regulon_names)

  f <- 0
  for(reg in regulon_names){
    for(reg2 in regulon_names){
      f <- f + 1
      # fraction_lower <- calc_csi(reg,reg2,pearson_cor)
      test_cor <- pearson_cor[reg,reg2]
      total_n <- ncol(pearson_cor)
      pearson_cor_sub <- subset(pearson_cor,rownames(pearson_cor) == reg | rownames(pearson_cor) == reg2)

      # sums <- apply(pearson_cor_sub,MARGIN = 2, FUN = compare_pcc, pcc = test_cor)
      sums <- apply(pearson_cor_sub, 2,function(m) ifelse( length(m[m>test_cor]) == length(m), 0, length(m)) )
      fraction_lower <- length(sums[sums == nrow(pearson_cor_sub)]) / total_n
      csi_regulons[f,] <- c(reg,reg2,fraction_lower)
    }
  }
  csi_regulons$CSI <- as.numeric(csi_regulons$CSI)
  return(csi_regulons)
}

plot_csi_modules <- function(
  csi_df,
  groupby,
  nclust = 4,
  row_anno = NULL,
  font_size_regulons = 6,
  vcolors = "viridis"
){
  ## subset csi data frame based on threshold
  csi_test_mat <- csi_df %>% tidyr::spread(regulon_2,CSI)

  future_rownames <- csi_test_mat$regulon_1
  csi_test_mat <- as.matrix(csi_test_mat[,2:ncol(csi_test_mat)])
  rownames(csi_test_mat) <- future_rownames

  color_map = c("#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b",
    "#666666","#1b9e77","#d95f02","#7570b3","#d01b2a","#43acde",
    "#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff",
    "#ecb9e5","#813139","#743fd2","#434b7e","#e6908e","#214a00",
    "#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7",
    "#006874","#d2ad2c","#b5d7a5","#9e8442","#4e1737","#e482a7",
    "#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99",
    "#3fdacb","#bf5b17")[1:nclust]
    
  names(color_map) <- 1:nclust
  color_use = list()

  color_use[[groupby]] <- color_map  ##different with local script 
  
  pheatmap::pheatmap(csi_test_mat,
           show_colnames = TRUE,
           #border_color = NA,
           color = colorRampPalette(continuous_palette[[vcolors]])(10),
           fontsize = 12,
           fontsize_row = font_size_regulons,
           fontsize_col = font_size_regulons,
           angle_col = 90,
           cutree_cols = nclust,
           cutree_rows = nclust,
           annotation_row = row_anno,
           annotation_colors = color_use,
           treeheight_row = 20,
           treeheight_col = 20,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           width = 4500,
           height = 4000,
           fontface = "bold")
}

calc_csi_module_activity <- function(
  clusters_df,
  regulonAUC_mat,
  metadata
){
  cell_types<- unique(metadata$cell_type)
  regulons <- unique(clusters_df$regulon)

  regulonAUC_mat <- regulonAUC_mat[as.character(regulons),]

  csi_activity_matrix_list <- list()
  csi_cluster_activity <- data.frame("csi_module" = c(),
                                    "mean_activity" = c(),
                                    "cell_type" = c())

  cell_type_counter <- 0
  groupby_1 =c()   # select group with only one cell
  for(i in names(table(metadata$cell_type)) ){
	if( table(metadata$cell_type)[i]==1 ){
		groupby_1= i
	}
  }
  cell_types <- setdiff(cell_types, groupby_1)
  
  regulon_counter <-
    for(ct in cell_types) {
      cell_type_counter <- cell_type_counter + 1
      cell_type_aucs <- rowMeans(regulonAUC_mat[,rownames(subset(metadata,cell_type == ct))])
	  
      cell_type_aucs_df <- data.frame("regulon" = names(cell_type_aucs),
                                      "activtiy"= cell_type_aucs,
                                      "cell_type" = ct)
      csi_activity_matrix_list[[ct]] <- cell_type_aucs_df
    }
  if(!is.null(groupby_1)){
    regulonAUC_mat_1 = regulonAUC_mat[,rownames(subset(metadata,cell_type == groupby_1))]
    csi_activity_matrix_list_1=data.frame("regulon" = names(regulonAUC_mat_1),
                                      "activtiy"= as.vector(regulonAUC_mat_1),
                                      "cell_type" = groupby_1)
    csi_activity_matrix_list[[groupby_1]] <- csi_activity_matrix_list_1
  }
  for(ct in names(csi_activity_matrix_list)){
    for(cluster in unique(clusters_df$csi_module)){
      csi_regulon <- subset(clusters_df,csi_module == cluster)
      csi_regulon_activtiy <- subset(csi_activity_matrix_list[[ct]],regulon %in% csi_regulon$regulon)
      csi_activtiy_mean <- mean(csi_regulon_activtiy$activtiy)
      this_cluster_ct_activity <- data.frame("csi_module" = cluster,
                                             "mean_activity" = csi_activtiy_mean,
                                             "cell_type" = ct)
      csi_cluster_activity <- rbind(csi_cluster_activity,this_cluster_ct_activity)
    }
  }

  csi_cluster_activity[is.na(csi_cluster_activity)] <- 0

  csi_cluster_activity_wide <- csi_cluster_activity %>%
    spread(cell_type,mean_activity)

  rownames(csi_cluster_activity_wide) <- csi_cluster_activity_wide$csi_module
  csi_cluster_activity_wide <- as.matrix(csi_cluster_activity_wide[2:ncol(csi_cluster_activity_wide)])

  return(csi_cluster_activity_wide)
}