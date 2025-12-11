net_plot = function(edge,vertex,
                    layout = 'layout_components',#'layout_with_lgl,layout_with_kk,layout_with_fr,layout_with_dh,layout_on_grid,layout_nicely,layout_components'
                    edge.color='#00CED1',#边颜色
                    edge.width=1,#边宽
                    edge.arrow.mode = 2,#箭头朝向，1为向后，2为向前，3为双向
                    edge.arrow.size=.5,#箭头大小
                    edge.label = '',#边标签
                    edge.label.cex = 1,#边标签字体大小
                    edge.label.color = 'black',#边标签字体颜色
                    edge.curved=.0,#线段曲度
                    edge.lty = 1,#线条类型，0 and “blank” mean no edges, 1 and “solid” are for solid lines, the other possible values are: 2 (“dashed”), 3 (“dotted”), 4 (“dotdash”), 5 (“longdash”), 6 (“twodash”)
                    
                    vertex.size=1,#节点大小
                    vertex.color='#1E90FF',#节点颜色
                    vertex.frame.color="grey",#节点边框颜色
                    vertex.label.color="black",#标签字体颜色
                    vertex.label.cex = 0.8,#节点标签字体大小
                    vertex.label.dist = 0,#标签与节点的距离
                    
                    legend.edge.label = 'F',#如需展示边得图例，可填写该参数，应该为向量，且向量中的元素个数应与length(unique(edge.color))一致
                    legend.vertex.label = 'F',#同上
                    legend.edge.color = '#00CED1',
                    legend.vertex.color = '#1E90FF',
                    mai = c(0.1,0.1,0.1,2.5),
                    
                    main='',
                    sub='',
                    xlab='',
                    ylab='',
                    
                    frame = FALSE#图是否加外框
                       ){
                library(igraph)
                set.seed(12346)
                if (is.data.frame(edge) != TRUE | is.data.frame(vertex) != TRUE){
                    stop('Make sure the edge/vertex is a data.frame')
                }
                net_pc<-graph_from_data_frame(edge, directed=TRUE, vertices = vertex)
                l <- do.call(layout, list(net_pc))
                set.seed(1235)
                if (legend.edge.label[1] != 'F' | legend.vertex.label[1] != 'F'){
                     par(mai = mai)
                }
                plot(net_pc, layout = l,
                edge.color=edge.color,
                edge.alpha=0.5,#边透明度
                edge.width = edge.width,
                edge.arrow.size=edge.arrow.size,
                edge.arrow.width=1,#箭头宽度
                edge.arrow.mode = edge.arrow.mode,
                edge.label=edge.label,
                edge.label.family='',#边标签字体
                edge.label.font=2,#边标签字体格式
                edge.label.cex = edge.label.cex,
                edge.label.color = 'black',
                edge.curved=edge.curved,
                edge.lty=edge.lty,# accepted as for the standard graphics par, .

                vertex.size=vertex.size,
                vertex.color=vertex.color,
                vertex.frame.color=vertex.frame.color,
                vertex.frame.width=0.5,#节点边框宽度
                vertex.label.color=vertex.label.color,
                vertex.shapes='pie',#节点形状
                vertex.label.family='Helvetica',#字体样式
                vertex.label.font = 2,#字体格式，节点标签使用的字库中的字体。1是纯文本，2是黑体，3是斜体，4是黑体和斜体，5指定符号字体
                vertex.label.cex = vertex.label.cex,#vertex_df$size/10,
                vertex.label.dist = vertex.label.dist,
                vertex.label.degree = pi/4,#它定义了节点标签的位置，相对于顶点的中心。它被解释为一个弧度的角度，零意味着 "向右"，而 "pi "意味着向左，向上是-pi/2，向下是pi/2.It defines the position of the vertex labels, relative to the center of the vertices. It is interpreted as an angle in radian, zero means ‘to the right’, and ‘pi’ means to the left, up is -pi/2 and down is pi/2.
                vertex.alpha = 0.1,

                margin =0,#图形到边缘的距离
                asp = 1,#一个数字常数，它给出了绘图的asp参数，即长宽比。如果你不想给一个长宽比，在这里提供0,默认为1
                frame = frame,#是否给一个框架,默认为FALSE
                main=main,#主图总标题
                sub=sub,#副标题
                xlab=xlab,ylab=ylab)
                if(legend.edge.label[1] != 'F'){
                    legend(1.2,1,legend = legend.edge.label,lwd = 4,lty = 1,pch = c(18,18),merge = FALSE,pt.cex = 1.5,xpd=TRUE,seg.len=1,
                    col = legend.edge.color,
                    border = vertex.frame.color,
                    bty = 'n')
                }
                if(legend.vertex.label[1] != 'F'){
                legend(1.2,0.85,legend = legend.vertex.label,pch = 16,pt.cex = 2,xpd=TRUE,pt.bg=legend.vertex.color,ncol=1,
                col = legend.vertex.color,
                border = vertex.frame.color,
                bty = 'n')
                legend(1.2,0.85,legend = legend.vertex.label,pch = 1,pt.cex = 2,xpd=TRUE,pt.bg=legend.vertex.color,ncol=1,
                col = vertex.frame.color,
                bty = 'n')
                }
}


process_tf <- function(tf, df, top, highlightTarget) {
  subset_df <- df[df$TF %in% tf, ]
  ordered_df = subset_df %>%
                  group_by(TF) %>%
                  slice_max(order_by = importance, n = top,
                          with_ties = FALSE) %>%
                  ungroup()

  highlightTarget_df = subset_df[subset_df$target %in% highlightTarget,]

  return(rbind(ordered_df,highlightTarget_df))
}
