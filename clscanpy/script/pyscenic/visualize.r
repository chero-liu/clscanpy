#!/usr/bin/env Rscript
source("/gpfs/oe-scrna/pipeline/scRNA-seq_further_analysis/function/seuratFc.r")
#source("/home/liuchenglong/script/Get_colors.R")
source("/public/scRNA_works/pipeline/scRNA-seq_further_analysis/Get_colors.R")

library(Seurat)
library(tidyr)
library(dplyr)
library(pheatmap)
library(optparse)
library(circlize)
library(ggplot2)
library(ComplexHeatmap)
library(SCENIC,lib.loc = '/gpfs/oe-software/scrna_software_bak/conda_env/Scenic/lib/R/library')
library(AUCell,lib.loc = '/gpfs/oe-software/scrna_software_bak/conda_env/Scenic/lib/R/library')
library(SCopeLoomR,lib.loc = '/gpfs/oe-software/scrna_software_bak/conda_env/Scenic/lib/R/library')

# 自定义函数模块 ---------------------------------------------------------------
handle_parameters <- function() {
    option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Path to the pyscenic sce_SCENIC.loom file path", metavar = "FILE"),
    make_option(c("-r", "--rds_filepath"), type = "character", help = "Path to input RDS file", metavar = "FILE"),
    make_option(c("-o", "--outdir"), type = "character", default = 'auto', help = "Output directory to save the CSV file"),
    make_option(c("-c", "--subnew_celltype"), type = "character", default = "all", help = "Subnew cell type"),
    make_option(c("-s", "--subsampleid"), type = "character", default = "all", help = "Subsample ID"),
    make_option(c("-g", "--subgroup"), type = "character", default = "all", help = "Subgroup"),
    make_option(c("-l", "--subcluster"), type = "character", default = "all", help = "Subcluster"),
    make_option(c("-p", "--predicate"), type = "character", default = 'all', help = "Filtering predicate (e.g.,"),
    make_option(c("--groupby_levels"), type = "character", default = 'no', help = "Tcell,Bcell,NK"),
    make_option(c("--groupby"), type = "character", default = "group_new_celltype", help = "Grouping variable"),
    make_option(c("-t", "--threshold"), type = "numeric", default = 0, help = "Threshold value"),
    make_option(c("--regulons_path"), type = "character", help = "Path to the regulons file"),
    make_option(c("--topGenes"), type = "integer", default = 3, help = "Number of top genes to consider"),
    make_option(c("-e", "--extended"), action = "store_true", default = FALSE, help = "Whether to use extended mode"),
    make_option(c("--nclust"), type = "integer", default = 4, help = "the number of csi modules to use"),
    make_option(c("--utils_path"), type = "character", help = ""),
    make_option(c("--combo"), type = "character", help = "combo metadata column"),
    make_option(c("--ccolors"), type = "character", default = "rdwibl",help = "连续性色板选择"),
    make_option(c("--vcolors"), type = "character", default = "viridis",help = "3.2 csi图色板选择"),
    make_option(c("--rowcluster"), type = "character", default="T",help="是否对行进行聚类，布尔值或 hclust 对象"),
    make_option(c("--colcluster"), type = "character", default="T",help="是否对列进行聚类，布尔值或 hclust 对象")
    )
  parse_args(OptionParser(option_list = option_list))
}

analyze_regulon_activity <- function(seurat_obj, regulonAUC, args) {
  # 注释处理
  cell_info <- seurat_obj@meta.data %>% arrange(!!sym(args$groupby))
  cell_info[[args$groupby]]<- factor(cell_info[[args$groupby]],levels = args$group_levels)
  # 热图生成
  annotation_colors <- list(
    setNames(getColorOrder(seurat_obj, args$groupby), 
            args$group_levels)
  )
  names(annotation_colors) <- args$groupby

  # 过滤0个数大于细胞数1/5的行
  regulonAUC_filtered = regulonAUC[rowSums(regulonAUC[, rownames(cell_info)] == 0) <= (ncol(regulonAUC[, rownames(cell_info)]) / 1.2),rownames(cell_info)]
  col_order <- order(cell_info[[args$groupby]])   # 或者按 levels 排序
  cell_info <- cell_info[col_order, ]
  col_list <- setNames(annotation_colors[[args$groupby]], levels(cell_info[[args$groupby]]))
annot_df <- data.frame(
  group = cell_info[[args$groupby]]
)
names(annot_df) <- args$groupby
annot_colors <- list()
annot_colors[[args$groupby]] <- annotation_colors[[args$groupby]]
col_ha <- HeatmapAnnotation(
  df = annot_df,
  col = annot_colors,
  annotation_name_gp = gpar(fontsize = 12),
  annotation_legend_param = list(
    labels_gp = gpar(fontsize = 12),
    title_gp = gpar(fontsize = 12)
  )
)
regulonAUC_scaled <- t(scale(t(regulonAUC_filtered)))
ht <- Heatmap(
  regulonAUC_scaled,
  name = "AUC",
  col = colorRamp2(
  breaks = c(seq(-2.5, 0, length.out = 100), seq(0, 2.5, length.out = 100)),
  colors = colorRampPalette(continuous_palette[[args$ccolors]])(200)
),show_row_dend = FALSE,
  top_annotation = col_ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 12)
  )
)
max_length <- max(nchar(args$group_levels))
width <- if (max_length < 15) {
  8 
} else {
  10 
}
pdf(file.path(args$outdir, "1.1.regulon_activity_heatmap_groupby_cells.pdf"),
    width = width, height = 6)
draw(ht)
dev.off()
detach("package:ComplexHeatmap", unload = TRUE)
  # activity处理
regulon_activity <- process_regulon_activity(
    cell_info, regulonAUC, args$groupby, nmads = 3
  )

  cbind(regulon = rownames(regulon_activity), regulon_activity) %>% write.table(file.path(args$outdir, "1.2.centered_regulon_activity_groupby_design.xls"), sep = '\t', row.names = FALSE, quote = FALSE)
  
  if (dim(regulon_activity)[1]<11) {
    hig <- 5
  } else {
    hig <- dim(regulon_activity)[1]*0.45
  }

  p_mean_pheatmap = pheatmap::pheatmap(regulon_activity,
                    cellwidth = 18,
                    cellheight = 18,
                    color= colorRampPalette(continuous_palette[[args$ccolors]])(200),
                    angle_col = "45",
                    treeheight_col=20, 
                    treeheight_row=20,
                    border_color=NA,
                    fontsize = 15,  
                    fontsize_row = 13,
                    fontsize_col = 15,
                    cluster_rows=as.logical(args$rowcluster),
                    cluster_cols=as.logical(args$colcluster),
                    fontface = "bold")

  save_plot(p_mean_pheatmap, 
              file.path(args$outdir, "1.3.regulon_activity_heatmap"),
              width = dim(regulon_activity)[2]/2+5,
              height = hig,
              dpi = 300,
              file_formats = c("pdf"))
}

#' RSS分析模块
perform_rss_analysis <- function(seurat_obj, regulonAUC, args) {
  rss_matrix <- calcRSS(regulonAUC, seurat_obj@meta.data[[args$groupby]])
  rss_df <- na.omit(rss_matrix) %>%
    as.data.frame()
  rss_df$regulon = rownames(rss_df)
  rss_df = gather(rss_df,key = !!sym(args$groupby), value = "RSS", -regulon)
  rss_df[[args$groupby]] = factor(rss_df[[args$groupby]],levels = args$group_levels)

  # 可视化
  rss_plot <- RSSRanking(
    rss_df, 
    group.by = args$groupby,
    top_genes = args$topGenes,
    plot_extended = args$extended,
    size = 4
  )
  write.table(rss_df, file.path(args$outdir, "2.1.regulon_RSS_annotation.xls"), sep = "\t", col.names = T, row.names =F, quote =F)
  save_plot(rss_plot, 
              file.path(args$outdir, "2.2.RSS_ranking_plot"),
              height = ceiling(length(unique(rss_df[[args$groupby]]))/2) * 4,width = 7,
              file_formats = c("pdf"))

  #筛选出至少有一个值大于 rss_threshold 的所有行
  if ( is.null(args$threshold) ){
    rss_threshold = 0
  }else{
    rss_threshold = as.numeric(args$threshold)
  }
  rss_specific <- na.omit(rss_matrix)[apply(na.omit(rss_matrix),MARGIN = 1 ,FUN = function(x) any(x > rss_threshold)),]
  rss_specific = rss_specific[,args$group_levels]
  p = pheatmap::pheatmap(rss_specific,
    cellwidth = 18,
    cellheight = 18,
    color=colorRampPalette(continuous_palette[[args$ccolors]])(299),
    angle_col = 45,
    treeheight_col=20,
    treeheight_row=20, 
    border_color=NA,
    fontsize = 15,  
    fontsize_row = 15,
    fontsize_col = 15,
    cluster_rows=as.logical(args$rowcluster),
    cluster_cols=as.logical(args$colcluster),
    fontface = "bold")
  save_plot(p,
    file.path(args$outdir, "2.3.RSS_heatmap"),
    width = dim(rss_specific)[2]/2+3,
    height = dim(rss_specific)[1]*0.5,
    dpi = 300,
    file_formats = c("pdf"))
}

#' CSI模块分析
perform_csi_analysis <- function(regulonAUC, args) {
  # CSI计算与聚类
  regulons_csi <- calculate_csi(regulonAUC)
  
  # 矩阵转换
  csi_matrix <- regulons_csi %>% tidyr::spread(regulon_2,CSI)
  future_rownames <- csi_matrix$regulon_1
  csi_matrix <- as.matrix(csi_matrix[,2:ncol(csi_matrix)])
  rownames(csi_matrix) <- future_rownames
  
  # 层次聚类
  hclust_result <- hclust(dist(csi_matrix, method = "euclidean"))
  clusters <- cutree(hclust_result, k = as.numeric(args$nclust))
  
  # 结果保存
  clusters_df <- data.frame(
    regulon = names(clusters),
    csi_module = clusters,
    row.names = NULL
  )
  write.table(
    clusters_df,
    file.path(args$outdir, "3.1.csi_module_annotation.xls"),
    sep = "\t", 
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
  
  # 热图注释数据准备
  row_annotation <- clusters_df
  rownames(row_annotation) <- NULL
  row_annotation$csi_module <- as.character(row_annotation$csi_module)
  rownames(row_annotation) <- row_annotation$regulon
  row_annotation <- row_annotation[, -which(names(row_annotation) == "regulon"),drop = F]

  # 可视化
  csi_plot <- plot_csi_modules(
    regulons_csi,
    args$groupby,
    nclust = as.numeric(args$nclust),
    row_anno = row_annotation,
    font_size_regulons = 9,
    vcolors = args$vcolors)
  save_plot(
    csi_plot,
    file.path(args$outdir, "3.2.regulons_csi_correlation_heatmap"),
    width = length(unique(regulons_csi[, 1])) * 0.15 + 2,
    height = length(unique(regulons_csi[, 1])) * 0.15,
    dpi = 300,
    file_formats = c("pdf")
  )

  return(clusters_df)
}

#' CSI模块活动分析
analyze_csi_activity <- function(clusters_df, seurat_obj, regulonAUC, args) {
  # 数据准备
  cell_metadata = seurat_obj@meta.data
  cell_metadata[,"cell_type"] = cell_metadata[,args$groupby]

  # 模块活动计算
  activity_matrix <- calc_csi_module_activity(
    clusters_df,
    regulonAUC,
    cell_metadata
  )
  activity_matrix = activity_matrix[,args$group_levels]
  # 热图可视化
  activity_heatmap <- pheatmap::pheatmap(
    activity_matrix,
    show_colnames = TRUE,
    color =  colorRampPalette(continuous_palette[[args$vcolors]])(10),
    cellwidth = 24,
    cellheight = 24,
    cluster_rows=as.logical(args$rowcluster),
    cluster_cols=as.logical(args$colcluster),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    fontsize = 15,  
    fontsize_row = 15,
    fontsize_col = 15,
    fontface = "bold"
  )

  if (dim(activity_matrix)[1]<7) {
    hig <- 6
  } else {
    hig <- dim(activity_matrix)[1]*0.45
  }
  save_plot(
    activity_heatmap,
    file.path(args$outdir, "3.3.csi_module_activity_heatmap"),
    width = dim(activity_matrix)[2]/2+2,
    height = hig,
    dpi = 150,
    file_formats = c("pdf")
  )
}

#
main <- function() {
  args <- handle_parameters()
  source(args$utils_path)
  
  dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
  
  data_objects <- load_data(args)
  # 确保是因子型，如果不是，转换为因子型
 if (args$groupby_levels == "no") {
  if (!is.factor(data_objects$seurat@meta.data[[args$groupby]])) {
    data_objects$seurat@meta.data[[args$groupby]] <- factor(data_objects$seurat@meta.data[[args$groupby]])
  }
  args$group_levels <- levels(data_objects$seurat@meta.data[[args$groupby]])
} else {
  args$group_levels <- strsplit(args$groupby_levels, ",")[[1]]
}
print("args$group_levels")
  print('regulon_activity 分析')
  analyze_regulon_activity(data_objects$seurat, data_objects$regulonAUC, args)

  print('rss 分析')
  perform_rss_analysis(data_objects$seurat, data_objects$regulonAUC, args)

  print('CSI分析')
  clusters_df <- perform_csi_analysis(data_objects$regulonAUC, args)
  
  print('CSI activity 分析')
  analyze_csi_activity(
    clusters_df,
    data_objects$seurat, 
    data_objects$regulonAUC, 
    args
  )
  if (!file.exists(file.path(args$outdir, "pySCENIC分析说明.docx"))) {
    file.copy(
      paste0(dirname(args$utils_path),"/pySCENIC分析说明.docx"),
      args$outdir
    )
  }
  setwd(args$outdir)
  print("Convert pdf to png...")
  system("module load ImageMagick && ls *.pdf | xargs -I {} -P 8 bash -c '/data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 200 -trim  ${0}  -quality 100  -flatten  ${0%.pdf}.png' {}")
  #system("find . -type f ! -name '*.pdf' ! -name '*.png' ! -name '*.xls'  ! -name '*.docx' -exec rm -f {} +") #删除多余产生的文件
}

# 执行主函数
if (!interactive()) {
  main()
}
