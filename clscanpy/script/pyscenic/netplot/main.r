#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
source('/gpfs/oe-scrna/liuchenglong/RaD/pyscenic/oepySCENIC/oepySCENIC/script/netplot/tools/utils.r')
source("/gpfs/oe-scrna/pipeline/scRNA-seq_further_analysis/function/seuratFc.r")

# 定义命令行参数
option_list <- list(
  make_option(c("--input"), type="character", 
              help="输入TSV文件路径", metavar="FILE"),
  make_option(c("--interestTF"), type="character", default="",
              help="感兴趣的TF列表（逗号分隔）"),
  make_option(c("--outdir"), type="character",
              help="输出目录 [默认当前目录]"),
  make_option(c("--targetTop"), type="integer", default=100,
              help="目标基因数量阈值 [默认 %default]"),
  make_option(c("--interestTFColor"), type="character", default="#1E90FF",
              help="TF节点颜色 [默认 %default]"),
  make_option(c("--interestTFSize"), type="integer", default=10,
              help="TF节点大小 [默认 %default]"),
  make_option(c("--showTargetTop"), type="double", default=0.3,
              help="显示目标基因比例 [默认 %default]"),
  make_option(c("--topTargetColor"), type="character", default="#E87D72",
              help="顶部目标基因颜色 [默认 %default]"),
  make_option(c("--topTargetSize"), type="integer", default=3,
              help="顶部目标基因大小 [默认 %default]"),
  make_option(c("--highlightTarget"), type="character", default="no",
              help="高亮目标基因列表（逗号分隔）"),
  make_option(c("--highlightTargetColor"), type="character", default="#00aa55",
              help="高亮目标颜色 [默认 %default]"),
  make_option(c("--highlightTargetSize"), type="integer", default=6,
              help="高亮目标大小 [默认 %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

main <- function(opt) {

  metadata <- read.csv(opt$input, sep = '\t')
  interestTF <- standardizeDelimitedVector(opt$interestTF)
  highlightTarget <- standardizeDelimitedVector(opt$highlightTarget)
  createDir(opt$outdir)
  

  edge_df <- process_tf(interestTF, metadata, opt$targetTop, highlightTarget)
  
  edge_df_top <- edge_df %>%
    group_by(TF) %>%
    slice_max(order_by = importance, n = opt$targetTop * opt$showTargetTop,
             with_ties = FALSE) %>%
    ungroup()
  showTarget <- unique(edge_df_top$target)
  
  allGenes <- unique(c(edge_df$target, edge_df$TF))
  vertex_df <- data.frame(name = allGenes, label = "")
  vertex_df$label <- ifelse(vertex_df$name %in% c(interestTF, showTarget, highlightTarget),
                           vertex_df$name, "")
  
  vertex_df <- vertex_df %>%
    mutate(
      size = as.integer(
        case_when(
          name %in% interestTF ~ opt$interestTFSize,
          name %in% showTarget ~ opt$topTargetSize,
          name %in% highlightTarget ~ opt$highlightTargetSize,
          TRUE ~ 1L  # 使用L后缀确保是整数
        )
      ),
      color = case_when(
        name %in% interestTF ~ opt$interestTFColor,
        name %in% showTarget ~ opt$topTargetColor,
        name %in% highlightTarget ~ opt$highlightTargetColor,
        TRUE ~ "#FFFFFF"
      )
    )

  write.csv(edge_df, paste0(opt$outdir,'/edge_df.csv'))
  write.csv(vertex_df, paste0(opt$outdir,'/vertex_df.csv'))

  numb <- ifelse(opt$highlightTarget != "no", 4, 3)
  
  generate_plot <- function(filename, width, height, device) {
    full_path <- file.path(opt$outdir, filename)
    
    if (device == "pdf") {
      pdf(file = full_path, width = width, height = height)
    } else if (device == "png") {
      png(filename = full_path, width = width, height = height)
    }
    
    net_plot(
      edge_df,
      vertex_df,
      layout = 'layout_with_lgl',
      edge.color = '#00CED1',
      vertex.label.cex = 0.8,
      vertex.size = vertex_df$size,
      vertex.color = vertex_df$color,
      legend.vertex.label = c('TF', 
                            paste0('Top ', opt$showTargetTop*100, '% of target'),
                            'Other target',
                            'highlightTarget')[1:numb],
      legend.vertex.color = c(opt$interestTFColor,
                            opt$topTargetColor,
                            "#FFFFFF",
                            opt$highlightTargetColor)[1:numb]
    )
    dev.off()
  }
  
  generate_plot("TFsNetRegulation.pdf", 12, 8, "pdf")
  generate_plot("TFsNetRegulation.png", 1200, 800, "png")
}

if (!interactive()) {
  main(opt)
}