# ------------------------------------------------------------------------------
# -*- coding: utf-8 -*-
# @Time    : 
# @Author  : liangxing.sun
# @FileName: filename
# @Software: Rstudio
# @Email: liangxing.sun@personalbio.cn
# @Description: Brief description of what the script does.
# ------------------------------------------------------------------------------

pacman::p_load(openxlsx,tidyr,ggplot2,ComplexHeatmap,deeptime,
               circlize,stringr,dplyr,ggraph,igraph,data.table,gridBase,clusterProfiler)

#
#' Title数据读入
#'
#' @param enrich_file 富集文件
#'
#' @return 富集文件读入的表格，如果有富集数据则将description按照50字符宽度换行
#' @export
#'
#' @examples
read_enrich_data<-function(enrich_file){
  enrich_data<-read.xlsx(enrich_file)
  enrich_data$Description<-str_wrap(enrich_data$Description,width = 50)
  return(enrich_data)
}

#
#' Title富集气泡图
#'
#' @param enrich_data 富集表读入的数据
#' @param top_n 取富集表的前几个 默认是20
#' @param order_value 按照什么排序，默认是按照Pvalue
#' @param enrich_type 富集类型，有GO、KEGG、InterPro、Reactome、DO、WikiPathway
#' @param output_dir 图形保存路径
#' @param contrast 比较组信息
#' @param p_color pvalue 映射的颜色按照显著到不显著附上颜色默认是红蓝
#' @param base_family 基础图形字体默认是sans
#' @param x_value x轴的值 在富集表中的列
#' @param y_value y轴的值 在富集表中的列
#' @param size_value 气泡值 在富集表中的列
#' @param color_value 颜色值 在富集表中的列
#' @param x_lab x轴标题。默认是Richfactor
#' @param y_lab y轴标题。默认是NULL,返回enrich_type+Description
#' @param title_lab 主标题。默认是enrich_type+Enrichment
#' @param color_lab 颜色图例标题，默认是NULL,返回"P value"
#' @param size_lab 大小图例标题，默认是 number
#'
#' @return
#' @export 在指定路径下输出一张富集因子图片
#'
#' @examples
create_richfactor_plot<-function(enrich_data,top_n=20,order_value="Pvalue",enrich_type,output_dir,contrast,
                                 p_color=c('#B2182B' ,'#2166AC'),base_family="sans",
                                 x_value="RichFactor",y_value="Description",size_value="DEP",color_value="Pvalue",
                                 x_lab="RichFactor",y_lab=NULL,title_lab=NULL,color_lab=NULL,size_lab="number"){
  if(is.null(y_lab)){
    y_lab = paste0(enrich_type," Description")
  }
  if(is.null(title_lab)){
    title_lab = paste0(enrich_type," Enrichment")
  }
  if(is.null(color_lab)){
    if(color_value == "Pvalue"){
      color_lab = "P value"
    }else if(color_value == "adjustPvalue"){
      color_lab = "adjust P"
    }
  }
  enrichment_final <- enrich_data %>% arrange(!!sym(order_value))
  enrichment_top <- enrichment_final[1:min(top_n,nrow(enrichment_final)), ]
  enrichment_top <- enrichment_top %>% arrange(desc(!!sym(order_value)))
  enrichment_top$Description <- factor(enrichment_top$Description,levels = enrichment_top$Description)
  bubbleplot<-ggplot(data = enrichment_top, aes(x = !!sym(x_value),
                                                y = !!sym(y_value), 
                                                size = !!sym(size_value), 
                                                color = !!sym(color_value))) +
    geom_point() +
    scale_color_gradientn(colours = p_color) +
    theme_classic(base_family = base_family)+
    labs(x = x_lab, y = y_lab,size=size_lab,
         title = title_lab,color=color_lab) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(color = "black", size = 13),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.y = element_text(color = "black", size = 13,vjust = 5),
      axis.text.x = element_text(color = "black", size = 10),
      plot.margin = margin(5,5,5,20)
    )+
    theme(
      panel.grid.major = element_blank(), # 主网格线
      panel.grid.minor = element_blank(), # 次网格线
      panel.border = element_rect(color = "grey", linewidth = 0.2, fill = NA),
      axis.line.x = element_blank(),
      axis.line.y = element_blank()
    )+guides(color=guide_colorbar(order=1),size=guide_legend(order=2))  
  ggsave(paste0(output_dir,"/",contrast, "_",enrich_type,"_richfactor.pdf"), bubbleplot, width = 10, height = 8)
}

#

#' Title富集柱状图
#'
#' @param enrich_data 富集表读入的数据
#' @param top_n 取富集表的前几个 默认是20
#' @param order_value 按照什么排序，默认是按照Pvalue
#' @param enrich_type 富集类型，有GO、KEGG、InterPro、Reactome、DO、WikiPathway
#' @param output_dir 图形保存路径
#' @param contrast 比较组信息
#' @param fill_color 柱状图
#' @param base_family 基础图形字体默认是sans
#' @param x_lab x轴标题。默认是"-Log10(P value)"
#' @param y_lab y轴标题。默认是NULL返回enrich_type+Description
#' @param title_lab 主标题。默认是NULL,返回enrich_type+Enrichment
#' @param size_lab 大小图例标题，默认是 number
#'
#' @return
#' @export 输出一张柱状图
#'
#' @examples
creat_barplot<-function(enrich_data,top_n=20,order_value="Pvalue",enrich_type,
                        output_dir,contrast,
                        fill_color,base_family="sans",
                        x_lab="-Log10(P value)",y_lab=NULL,title_lab=NULL,size_lab="number"){
  if(is.null(y_lab)){
    y_lab = paste0(enrich_type," Description")
  }
  if(is.null(title_lab)){
    title_lab = paste0(enrich_type," Enrichment")
  }
  
  enrichment_final <- enrich_data %>% arrange(!!sym(order_value))
  enrichment_top <- enrichment_final[1:min(top_n,nrow(enrichment_final)), ]
  enrichment_top$`log10(p-value)` <- log10(enrichment_top$Pvalue)
  
  
  if(enrich_type == "GO"){
    enrichment_top$classify<-enrichment_top$Category
    enrichment_top<-enrichment_top%>%group_by(classify)%>%arrange(desc(classify),desc(`log10(p-value)`))
  }else if(enrich_type == "KEGG"){
    enrichment_top$classify<-enrichment_top$level1
    enrichment_top<-enrichment_top%>%group_by(classify)%>%arrange(desc(classify),desc(`log10(p-value)`))
  }else{
    enrichment_top$classify<-enrichment_top$Description
    enrichment_top<-enrichment_top%>%group_by(classify)%>%arrange(desc(`log10(p-value)`))
  }
  
  enrichment_top$Description <- factor(enrichment_top$Description, levels = enrichment_top$Description)
  barplot<-ggplot(data = enrichment_top, aes(x = -`log10(p-value)`, y = Description,fill = classify))
  
  if(enrich_type == "GO" || enrich_type == "KEGG"){
    barplot<-barplot+geom_bar(stat = "identity")
  }else{
    barplot<-barplot+geom_bar(stat = "identity", fill = fill_color[1:min(top_n,nrow(enrichment_top))],show.legend = F)
  } 
  barplot <- barplot+theme_classic(base_family=base_family)+
    labs(x = x_lab, y = y_lab,title=title_lab) +
    theme(
      plot.title=element_text(hjust=0.5),
      axis.title.x = element_text(color = "black", size = 13),
      axis.title.y = element_text(color = "black", size = 13,vjust = 5),
      axis.text.y = element_text(color = "black", size = 10),
      axis.text.x = element_text(color = "black", size = 10),
      plot.margin = margin(5,5,5,20),
      legend.title = element_blank()
    ) 
  ggsave(paste0(output_dir, "/",contrast,"_",enrich_type,"_barplot.pdf"), barplot, width = 10, height = 8)
}


#富集网络图

#' Title富集网络图
#'
#' @param outdir 输出路径
#' @param enrich_data 富集数据
#' @param top_n 前多少个默认20
#' @param contrast 比较组信息
#' @param enrich_type 富集类型
#' @param line_color 线颜色默认是"#F2D4D1"
#' @param point_color 点颜色默认是c( '#2166AC','#B2182B')
#' @param geneSet 读取富集表富集到此次表的列名默认是DEP，选取Up_Protein和Dwon_Protein
#'
#' @return
#' @export 输出网络图、gml文件格式、输出xlsx网络图表格
#'
#' @examples
network_plot<-function(outdir,enrich_data,top_n=20,contrast,enrich_type,line_color="#F2D4D1",point_color=c( '#2166AC','#B2182B'),geneSet="DEP"){
  Result_enrichment<-enrich_data
  type=enrich_type
  gene_strsplit <- function(gene_list) {
    if (gene_list != "-") {
      split_gene <- unlist(strsplit(gene_list, ";"))
    } else {
      split_gene <- NULL
    }
    split_gene
  }
  get_change_label <- function(label) {
    label <- ifelse(ifelse(label == "P-value", "#-log10(Pvalue)",
                           ifelse(label == "adjusted-P", "#-log10(adjustP)",
                                  ifelse(label == "DEP", "DEG_value",
                                         ifelse(label == "-log10(Pvalue)", "#-log10(Pvalue)",
                                                ifelse(label == "-log10(adjustP)", "#-log10(adjustP)",
                                                       label
                                                )
                                         )
                                  )
                           )
    ))
    return(label)
  }
  DDR_filter_label <- function(Result_enrichment, top_n = NULL) {
    sub_enrichment <- as.data.frame(Result_enrichment) %>%dplyr::arrange(Pvalue)%>%
      dplyr::select(1, Pvalue, adjustPvalue, RichFactor) 
    if (is.numeric(top_n)) {
      # sub_enrichment <- sub_enrichment %>% slice_min(Pvalue,n=top_n)
      sub_enrichment <- sub_enrichment %>%
        # arrange(desc(RichFactor), Pvalue) %>%
        arrange(Pvalue, desc(RichFactor)) %>%
        head(top_n)
      label_contant <- paste0(sub_enrichment[, 1], collapse = "//")
    }
    
    return(label_contant)
  }
  
  common_gene <- function(x, df, geneSet) {
    if (geneSet == "DEP") {
      GO1_upgene <- gene_strsplit(df[x[1], "Up_protein"])
      GO1_downgene <- gene_strsplit(df[x[1], "Down_protein"])
      GO1_gene <- c(GO1_upgene, GO1_downgene)
      
      GO2_upgene <- gene_strsplit(df[x[2], "Up_protein"])
      GO2_downgene <- gene_strsplit(df[x[2], "Down_protein"])
      GO2_gene <- c(GO2_upgene, GO2_downgene)
      
      common_gene_n <- length(intersect(GO1_gene, GO2_gene))
      common_gene_n_list <- paste0(common_gene_n, "//", paste0(intersect(GO1_gene, GO2_gene), collapse = ";"))
      
      return(common_gene_n_list)
    } else {
      #geneSet <- ifelse(geneSet == "up" | geneSet == "Up", "Up_protein", "Down_protein")
      GO1_gene <- gene_strsplit(df[x[1], geneSet])
      
      GO2_gene <- gene_strsplit(df[x[2], geneSet])
      
      common_gene_n <- length(intersect(GO1_gene, GO2_gene))
      common_gene_n_list <- paste0(common_gene_n, "//", paste0(intersect(GO1_gene, GO2_gene), collapse = ";"))
      return(common_gene_n_list)
    }
  }
  
  network_list<-DDR_filter_label(Result_enrichment,top_n)
  
  ## step1 - 对富集分结果做过滤
  label_list2 <- unlist(strsplit(network_list, "//"))
  NetPlot_df <- as.data.frame(Result_enrichment)
  colnames(NetPlot_df)[1] <- "ID"
  
  NetPlot_df <- NetPlot_df[which(NetPlot_df$ID %in% label_list2), ]
  rownames(NetPlot_df) <- NetPlot_df$ID
  NetPlot_df$Pvalue<-(-log10(NetPlot_df$Pvalue))
  NetPlot_df<-NetPlot_df%>%dplyr::rename(`-log10(Pvalue)`=Pvalue)
  ## step2 - 构建边的数据
  # if(length(NetPlot_df[, 1]) > 2) 
  edge_df <- as.data.frame(t(combn(NetPlot_df[, 1], 2))) ## 构建 from 和 to;即构建边的数据
  i=1
  test<-c()
  for(i in 1:nrow(edge_df)){
    x1=edge_df[i,1]
    x2=edge_df[i,2]
    x=c(x1,x2)
    test<-c(test,common_gene(x,NetPlot_df,geneSet=geneSet))
  }
  edge_df$inter_Count<-test
  edge_df <- edge_df %>% separate(inter_Count, into = c("inter_Count", "common_diffgene"), sep = "//")
  edge_df <- edge_df[edge_df$inter_Count != 0, ]
  ## 保存网络图表格
  out_network_df <- edge_df
  names(out_network_df) <- c("Node1", "Node2", "Common_Diffgene_Count", "Common_Diffgene")
  if (type == "GO") {
    # Description
    out_network_df <- left_join(out_network_df, NetPlot_df[, c("ID", "Description")], by = c("Node1" = "ID"))
    out_network_df <- left_join(out_network_df, NetPlot_df[, c("ID", "Description")], by = c("Node2" = "ID"))
  } else if (type == "KEGG") {
    # Description
    out_network_df <- left_join(out_network_df, NetPlot_df[, c("ID", "Description")], by = c("Node1" = "ID"))
    out_network_df <- left_join(out_network_df, NetPlot_df[, c("ID", "Description")], by = c("Node2" = "ID"))
  } else {
    # Description
    out_network_df <- left_join(out_network_df, NetPlot_df[, c("ID", "Description")], by = c("Node1" = "ID"))
    out_network_df <- left_join(out_network_df, NetPlot_df[, c("ID", "Description")], by = c("Node2" = "ID"))
  }
  names(out_network_df) <- c("Node1", "Node2", "Common_Diffgene_Count", "Common_Diffgene", "Node1_Description", "Node2_Description")
  out_network_df <- out_network_df[, c("Node1", "Node1_Description", "Node2", "Node2_Description", "Common_Diffgene_Count", "Common_Diffgene")]
  
  write.xlsx(out_network_df, paste0(outdir, "/", contrast, "_", type, "_network.xlsx"))
  
  edge_df$common_diffgene <- NULL
  names(edge_df) <- c("from", "to", "inter_Count")
  edge_df$inter_Count <- as.numeric(edge_df$inter_Count)
  
  ## step3 - 构建点的数据
  node_df <- NetPlot_df
  ## step4 点和边的数据构建网络
  net <- graph_from_data_frame(d = edge_df, vertices = node_df, directed = FALSE)
  ## 按画图需求构建 点和边 的属性
  node_color_label<-"-log10(Pvalue)"
  vertex_attr(net)$Degree <- igraph::degree(net) ## 计算 Degree; 每个节点连接的节点数 --> 对应每个节点的公共基因数？
  label_list <- c("RichFactor", "Degree", "DEP_Number", "-log10(Pvalue)", "-log10(adjustP)")
  vertex_attr(net)$node_color_label <- vertex_attr(net)[[node_color_label]]
  vsize <- vertex_attr(net)[["DEP"]] %>%
    abs() %>%
    round(1)
  vertex_attr(net)$size <- vsize
  
  node_size_label<-"DEP_Number"
  p1<-ggraph(net, layout ="nicely") +
    geom_edge_arc(aes(edge_width = inter_Count), color = line_color, strength = 0.2, angle_calc = "along") +
    scale_edge_width(range = c(0.5, 2)) ## 线宽的系数 最大值和最小值
  if(type == "GO"){
    p1<-p1+geom_node_point(aes(size = size, color = node_color_label,shape = Category))
  }else{
    p1<-p1+geom_node_point(aes(size = size, color = node_color_label))
  }
    p1<-p1+
    scale_color_gradientn(colours = point_color)+
    labs(fill = "-Log10(P value)", color = "-Log10(P value)", size = all_of(node_size_label),edge_width="inter_Count") +
    scale_size(range = c(3, 8))+
    geom_node_text(aes(label = Description), color = "black", check_overlap = F, repel = T, min.segment.length = 0, show.legend = F)+
    theme_graph()+
    guides(color=guide_colorbar(order=1),size=guide_legend(order=3),edge_width=guide_legend(order=2))
  
  
  
  # 补充一个gml文件供下载，可以直接导入cytoscape或gephi
  write.graph(net, file = paste0(outdir, "/", contrast, "_", type, "_network.gml"), format = "gml") # 使用该格式保存的话，需要对net的数据做数据的限制，即，哪些数据是需要输出的，哪些是不需要的
  # 保存图片
  ggsave(paste0(outdir,"/",contrast, "_", type, "_network.pdf"),plot=p1,width = 12, height = 10)
}

#' Title 富集圈图
#'
#' @param enrichment_final 富集表格
#' @param top_n 前多少个默认20
#' @param p_color p值颜色c( '#2166AC','#B2182B')
#' @param save_file 保存路径
#' @param updown_color 上下调颜色c("#f5780c","#0586c0")
#'
#' @return
#' @export
#'
#' @examples
sub_circlize_plot <- function(enrichment_final,top_n=20,p_color=c( '#4E9323',"#FFCA80",'#C72281'),save_file,updown_color = c("#FF964A","#5FC6FF")) {
  # 首先获取当前的 circos 参数
  # current_circos_par <- circos.par()
  
  # 将 cell.padding 的第二个和第四个值设定为0
  # circos.par(cell.padding = c(current_circos_par$cell.padding[1], 0, current_circos_par$cell.padding[3], 0))
  
  # 在修改参数后，绘制你的环形图（此处需要你的其他代码来实际绘制环形图）
  
  # 第一个圈：绘制id
  
  circlize_df <- enrichment_final
  circlize_df$Pvalue<-(-log10(circlize_df$Pvalue))
  circlize_df$adjustPvalue<-(-log10(circlize_df$adjustPvalue))
  circlize_df<-circlize_df%>%dplyr::rename(`#-log10(Pvalue)`=Pvalue,`#-log10(adjustP)`=adjustPvalue)
  colnames(circlize_df)[1] <- "ID"
  if ("level1" %in% colnames(circlize_df)) {
    colnames(circlize_df)[which(colnames(circlize_df) %in% "level1")] <- "Category"
  }else if("Category" %in% colnames(circlize_df)) {
    colnames(circlize_df)[which(colnames(circlize_df) %in% "Category")] <- "Category"
  }else if ("Description" %in% colnames(circlize_df)) {
    colnames(circlize_df)[which(colnames(circlize_df) %in% "Description")] <- "Category"
  }
  
  circlize_label <-  circlize_df%>%arrange(desc(`#-log10(Pvalue)`))%>%head(top_n)## 对结果做过滤,使用过滤后的结果进行画图
  circlize_df <- circlize_df[which(circlize_df$ID %in% circlize_label$ID), ]
  circlize_df <- circlize_df[order(circlize_df$Category), ] ## 按照 Category 进行排序，以确保后续的图的颜色连续
  
  ## 构建圈图数据
  circlize_df$gene_num.min <- 0
  circlize_df$gene_num.rich <- as.numeric(unlist(lapply(circlize_df$BgRatio, function(x) strsplit(x, "/")[[1]][1]))) # 该通路中的基因数
  first_max_tick<-NULL
  if (is.null(first_max_tick) || first_max_tick == 0) {
    circlize_df$gene_num.max <- rep(max(circlize_df$gene_num.rich))
  } else if (first_max_tick < max(circlize_df$gene_num.rich)) {
    circlize_df$gene_num.max <- rep(as.numeric(first_max_tick))
  } else if (first_max_tick >= max(circlize_df$gene_num.rich)) {
    circlize_df$gene_num.max <- rep(max(circlize_df$gene_num.rich))
  }
  rownames(circlize_df) <- circlize_df$ID
  circlize_df$ID <- factor(rownames(circlize_df), levels = rownames(circlize_df))
  
  col=color.ls[["Classic"]]
  col_df <- data.table(Category = unique(circlize_df$Category), col = col[1:length(unique(circlize_df$Category))])
  first_color <- col_df[match(circlize_df$Category, col_df$Category), 2]
  pdf(paste0(save_file, "_circos.pdf"),width = 11, height = 7)
  showtext_auto()
  circle_size <- unit(1, "snpc")
  circos.par(gap.degree = 0.5, start.degree = 90)
  plot_data <- circlize_df[, c("ID", "gene_num.min", "gene_num.max")]
  
  plot.new()
  pushViewport(viewport(
    x = 0, y = 0.5, width = circle_size, height = circle_size,
    just = c("left", "center")
  ))
  par(omi = gridOMI(), new = TRUE)
  
  circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)
  circos.track(
    ylim = c(0, 1), track.margin = c(0, 0), track.height = 0.08, bg.border = NA, bg.col = first_color[[1]],
    panel.fun = function(x, y) {
      ylim <- get.cell.meta.data("ycenter")
      xlim <- get.cell.meta.data("xcenter")
      sector.name <- get.cell.meta.data("sector.index")
      circos.axis(h = "top", labels.cex = 0.7, labels.pos.adjust = T, labels.col = "black", labels.niceFacing = FALSE) ## 基因数据标签
      circos.text(xlim, ylim, sector.name, cex = 0.8, col = "black",  niceFacing = FALSE) ## ID 标签
    }
  )
  circlize_second="P-value"
  # 第二圈，绘制富集的基因和富集p值
  if (circlize_second == "P-value") {
    plot_data <- circlize_df[, c("ID", "gene_num.min", "gene_num.rich", "#-log10(Pvalue)")]
    second_data <- "#-log10(Pvalue)"
    title_second_data <- "-Log10(P value)"
  } else {
    # circlize_second=="adjusted-P"
    plot_data <- circlize_df[, c("ID", "gene_num.min", "gene_num.rich", "#-log10(adjustP)")]
    second_data <- "#-log10(adjustP)"
    title_second_data <- "-Log10(adjustP value)"
  }
  
  label_data <- circlize_df[, "gene_num.rich"]
  label_data <- as.data.frame(label_data)
  rownames(label_data) <- circlize_df$ID
  p_max <- round(max(circlize_df[[second_data]])) + 1
  colorsChoice <- colorRampPalette(p_color)
  color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))
  
  circos.genomicTrackPlotRegion(
    plot_data,
    track.margin = c(0, 0), track.height = 0.1, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
      ylim <- get.cell.meta.data("ycenter")
      xlim <- label_data[get.cell.meta.data("sector.index"), 1] / 2
      sector.name <- label_data[get.cell.meta.data("sector.index"), 1]
      circos.text(xlim, ylim, sector.name, cex = 0.8, col = "black",  niceFacing = FALSE)
    }
  )
  
  ## 第三圈：上下调基因数目
  ## 数据准备
  
  if ("DEP" %in% colnames(circlize_df)) {
    ## 数据准备
    ## 上调 + 下调
    circlize_df$up.proportion <- circlize_df$Up / circlize_df$DEP
    circlize_df$down.proportion <- circlize_df$Down / circlize_df$DEP
    circlize_df$up_v <- circlize_df$up.proportion * circlize_df$gene_num.max
    plot_data_up <- circlize_df[, c("ID", "gene_num.min", "up_v")]
    names(plot_data_up) <- c("id", "start", "end")
    plot_data_up$type <- 1
    
    circlize_df$down_v <- circlize_df$down.proportion * circlize_df$gene_num.max + circlize_df$up_v
    plot_data_down <- circlize_df[, c("ID", "up_v", "down_v")]
    names(plot_data_down) <- c("id", "start", "end")
    plot_data_down$type <- 2
    
    plot_data <- rbind(plot_data_up, plot_data_down)
    label_data <- circlize_df[, c("up_v", "down_v", "Up", "Down")]
    label_data <- as.data.frame(label_data)
    rownames(label_data) <- circlize_df$ID
    color_assign <- colorRamp2(breaks = c(1, 2), col = updown_color)
    
    ## 画图
    suppressWarnings(suppressMessages(circos.genomicTrackPlotRegion(
      plot_data,
      track.margin = c(0, 0), track.height = 0.1, bg.border = NA, stack = TRUE,
      panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
        ylim <- get.cell.meta.data("cell.bottom.radius") - 0.5
        xlim <- label_data[get.cell.meta.data("sector.index"), 1] / 2
        sector.name <- ifelse(label_data[get.cell.meta.data("sector.index"), 3] > 0, label_data[get.cell.meta.data("sector.index"), 3], " ")
        circos.text(xlim, ylim, sector.name, cex = 0.6, col = "black", niceFacing = FALSE)
        xlim <- (label_data[get.cell.meta.data("sector.index"), 2] + label_data[get.cell.meta.data("sector.index"), 1]) / 2
        sector.name <- ifelse(label_data[get.cell.meta.data("sector.index"), 4] > 0, label_data[get.cell.meta.data("sector.index"), 4], " ")
        circos.text(xlim, ylim, sector.name, cex = 0.6, col = "black", niceFacing = FALSE)
      }
    )))
  } else {
    ## 上调
    plot_data <- circlize_df[, c("ID", "gene_num.min", "gene_num.rich")]
    label_v <- ifelse("Up" %in% colnames(circlize_df), "Up", "Down")
    label_col <- ifelse("Up" %in% colnames(circlize_df), upcol, downcol)
    label_data <- as.data.frame(cbind(circlize_df[, "gene_num.rich"], circlize_df[[label_v]]))
    rownames(label_data) <- circlize_df$ID
    
    ## 画图
    suppressWarnings(suppressMessages(circos.genomicTrackPlotRegion(
      plot_data,
      track.margin = c(0, 0), track.height = 0.1, bg.border = NA, stack = TRUE,
      panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = upcol, border = NA, ...)
        ylim <- get.cell.meta.data("cell.bottom.radius") - 0.6
        xlim <- label_data[get.cell.meta.data("sector.index"), 1] / 2
        sector.name <- label_data[get.cell.meta.data("sector.index"), 2]
        circos.text(xlim, ylim, sector.name, cex = third_size, col = third_label_col, font = par(third_font), niceFacing = FALSE)
      }
    )))
  }
  
  # 第四圈，绘制富集因子
  plot_data <- circlize_df[, c("ID", "gene_num.min", "gene_num.max", "RichFactor")]
  label_data <- circlize_df[, c("Category", "RichFactor")]
  label_data$RichFactor <- round(label_data$RichFactor, 2)
  label_data <- as.data.frame(label_data)
  rownames(label_data) <- circlize_df$ID
  
  fcol_df <- as.data.frame(col_df[match(circlize_df$Category, col_df$Category), ])
  rownames(fcol_df) <- circlize_df$ID
  show_RichFactor<-F
  if (show_RichFactor) {
    ## 需要显示 RichFactor 的数字
    suppressWarnings(suppressMessages(circos.genomicTrack(
      plot_data,
      track.margin = c(0.01, 0.04), track.height = 0.3, ylim = c(0.05, 0.95), bg.col = "gray95", bg.border = NA, # ylim = c(0, 1)
      panel.fun = function(region, value, ...) {
        sector.name <- get.cell.meta.data("sector.index")
        for (i in seq(0.1, 0.9, by = 0.1)) {
          circos.lines(c(0, max(region)), c(i, i), col = "gray", lwd = 0.3)
        }
        circos.genomicRect(region, value, col = fcol_df[sector.name, 2], border = NA, ytop.column = 1, ybottom = 0, ...)
        ylim <- label_data[get.cell.meta.data("sector.index"), 2] - 0.05
        xlim <- get.cell.meta.data("xcenter")
        sector.name <- label_data[get.cell.meta.data("sector.index"), 2]
        circos.text(xlim, ylim, sector.name, cex = 0.6, col = "black", niceFacing = FALSE)
      }
    )))
  } else {
    ## 默认情况下，不显示 RichFactor 的数字
    circos.genomicTrack(
      plot_data,
      track.margin = c(0.01, 0.04), track.height = 0.3, ylim = c(0, max(circlize_df$RichFactor)*1.3), bg.col = "gray95", bg.border = NA,
      panel.fun = function(region, value, ...) {
        sector.name <- get.cell.meta.data("sector.index")
        for (i in seq(max(circlize_df$RichFactor)*0.13, max(circlize_df$RichFactor)*1.3, by = 0.1)) {
          circos.lines(c(0, max(region)), c(i, i), col = "gray", lwd = 0.3)
        }
        circos.genomicRect(region, value, col = fcol_df[sector.name, 2], border = NA, ytop.column = 1, ybottom = 0, ...)
      }
    )
  }
  circos.clear()
  
  # 图例
  ## 中心图例
  updown_legend <- Legend(
    labels = c("Number of Proteins", "Up-regulated", "Down-regulated", "Rich Factor"),
    graphics = list(
      function(x = 0, y, w, h) {
        grid.draw(gTree(
          children = gList(
            rectGrob(x = 0, y, w * 2.2, h * 1, gp = gpar(fill = "gray80", col = "gray80")),
            textGrob("100", x = 0, y)
          ),
          gp = gpar(col = "black", fontsize = 10, fontface = "bold")
        ))
      },
      function(x = 0, y, w, h) {
        grid.rect(x = 0, y, w * 2.2, h * 1, gp = gpar(fill = updown_color[1], col = updown_color[1]))
      },
      function(x = 0, y, w, h) {
        grid.rect(x = 0, y, w * 2.2, h * 1, gp = gpar(fill = updown_color[2], col = updown_color[2]))
      },
      function(x = 0, y, w, h) {
        grid.polygon(y = c(0.09, -0.05, -0.12, 0.16), x = c(-0.13, -0.13, 0.13, 0.13), gp = gpar(fill = "gray85", col = "gray85"))
      }
    ),
    labels_gp = gpar(fontsize = 10), # grid_height = unit(0.6, 'cm'), grid_width = unit(0.8, 'cm'),
    row_gap = unit(2, "mm")
  )
  
  ## 第四圈图例
  category_legend <- Legend(
    labels = unique(fcol_df)$Category, # 各二级分类的名称
    type = "points", pch = NA, background = unique(fcol_df)$col, # 各二级分类的颜色
    labels_gp = gpar(fontsize = 10), row_gap = unit(1, "mm")
  ) # grid_height = unit(0.6, 'cm'), grid_width = unit(0.8, 'cm'))
  
  ## 第二圈图例
  pvalue_legend <- Legend(
    col_fun = colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1)),
    legend_height = unit(3, "cm"), labels_gp = gpar(fontsize = 10),
    title_gp = gpar(fontsize = 10), title_position = "topleft", title = paste0(title_second_data, "\n"), row_gap = unit(3, "mm")
  )
  
  lgd_list_vertical <- packLegend(updown_legend)
  lgd_list_vertical2 <- packLegend(category_legend, pvalue_legend)
  
  draw(lgd_list_vertical, just = "center")
  upViewport()
  draw(lgd_list_vertical2, x = circle_size, just = "left")
  dev.off()
  # 如果你想在作图后恢复默认的 circos 参数
  # circos.par(current_circos_par)
}

#' Title 和弦图
#'
#' @param enrich_data 富集表格
#' @param dep 差异表
#' @param order_value 按照什么排序 默认Pvalue
#' @param enrich_type 富集类型
#' @param outdir 输出路径
#' @param contrast 比较组
#' @param up_down_color 上下调颜色默认：c("#f5780c","white","#0586c0")
#' @param pathway_color 通路颜色默认color.ls[["Funny"]]
#'
#' @return
#' @export 输出一张和弦图
#'
#' @examples
chordfun<-function(enrich_data,dep,order_value="Pvalue",enrich_type,outdir,
                   contrast,up_down_color = c("#f5780c","white","#0586c0"),pathway_color=color.ls[["Funny"]]){
  plotdata<-enrich_data%>%arrange(!!sym(order_value))
  plotdata<-plotdata[1:min(10,nrow(plotdata)),]
  plotdata<-plotdata%>%select(ID,Description,Pvalue,Up_protein,Down_protein)%>%
    mutate(`-log10(Pvalue)`= -log10(Pvalue))
  
  UP_Pro<-plotdata%>%select(ID,Up_protein)%>%
    separate_rows(Up_protein,sep = ";")%>%
    mutate(Up_protein = gsub("\\(.+","",Up_protein))%>%
    rename(from = ID,to=Up_protein)
  
  DOWN_Pro<-plotdata%>%select(ID,Down_protein)%>%
    separate_rows(Down_protein,sep = ";")%>%
    mutate(Down_protein = gsub("\\(.+","",Down_protein))%>%
    rename(from = ID,to=Down_protein)
  df1<-rbind(UP_Pro,DOWN_Pro)%>%filter(to!="-")
  color_fun<-function(log2FC,updown_color){
    if(min(log2FC)>0){
      col_fun = colorRamp2(c(0, max(log2FC)), updown_color[2:1])
      my_color<-col_fun(log2FC)
    }else if(max(log2FC)<0){
      col_fun = colorRamp2(c(min(log2FC),0), updown_color[3:2])
      my_color<-col_fun(log2FC)
    }else{
      col_fun = colorRamp2(c(-max(abs(log2FC)),0,max(abs(log2FC))),updown_color[3:1] )
      my_color<-col_fun(log2FC)
    }
    return(my_color)
  }
  #log2FC
  dep<-dep%>%select(Protein.Accessions,Gene.Name,contains("Log2FC"))%>%
    arrange(desc(abs(select(.,contains("Log2FC")))))
  
  df2<-left_join(df1,dep,join_by(to==Protein.Accessions),relationship = "many-to-many")%>%
    unique()%>%
    mutate(ABS=abs(select(.,contains("Log2FC"))))%>%
    group_by(from)%>%
    arrange(desc(ABS),.by_group = T)%>%
    do(head(., n = 10))%>%unique()
  colnames(df2)[c(2,4)]<-c("Accession","log2FC")
  df2<-df2[order(df2$log2FC,decreasing = T),]
  gene_color<-color_fun(df2$log2FC,updown_color = up_down_color)
  names(gene_color)<-df2$Accession
  df3<-df2 %>% left_join(.,select(plotdata,ID,!!sym(order_value)),by=join_by(from==ID))%>%
    arrange(!!sym(order_value))
  pathway_color<-pathway_color[1:length(unique(df3$from))]
  names(pathway_color)<-unique(df3$from)
  
  grid.col<-c(pathway_color,gene_color)
  if(length(unique(df2$Accession))<50){
    cex=0.9
  }else if(length(unique(df2$Accession))<70){
    cex=0.8
  }else if(length(unique(df2$Accession))<90){
    cex=0.7
  }else{
    cex=0.6
  }
  filename=paste0(outdir,"/",contrast,"_",enrich_type,"_chord.pdf")
  #绘图
  pdf(filename,width = 16,height = 10)
  par(mar = c(0, 0, 0, 1))
  circos.par(start.degree = 90,circle.margin=c(0.01,0.8,0.01,0.01),gap.after=1)
  chordDiagram(df3[,c(1,2)], annotationTrack = "grid", order = c(unique(df3$from),rev(df2$Accession)),grid.col = grid.col,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df3))))))
  # we go back to the first track and customize sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = T, adj = c(0, 0.5),cex =cex)
  }, bg.border = NA) # here set bg.border to NA is important
  circos.clear()
  
  col_fun1<-function(log2FC,updown_color){
    if(min(log2FC)>0){
      col_fun = colorRamp2(c(0, max(log2FC)), updown_color[2:1])
    }else if(max(log2FC)<0){
      col_fun = colorRamp2(c(min(log2FC),0), updown_color[3:2])
    }else{
      col_fun = colorRamp2(c(-max(abs(log2FC)),0,max(abs(log2FC))), updown_color[3:1])
    }
    return(col_fun)
  }
  
  
  log2_leged<-Legend(col_fun = col_fun1(df2$log2FC,updown_color=up_down_color),
                     title = "Log2(Fold change)", direction = "horizontal")
  descipt<-Legend(labels =str_wrap(paste(plotdata[,1],plotdata[,2]),width = 50) , title = "Description", legend_gp = gpar(fill = pathway_color),labels_gp = gpar(fontsize = 12))
  legendlist<-packLegend(log2_leged, descipt)
  draw(legendlist,x = unit(13, "in"), y = unit(5, "in"))
  
  dev.off()
}

#' Title 上下调图片
#'
#' @param enrichment_final 富集表格
#' @param order_value 按照什么排序 默认Pvalue
#' @param enrich_type 富集类型
#' @param top_n 前多少默认20
#' @param updown_color  上下调颜色 c("#f5780c","#0586c0")
#' @param output_dir 输出路径
#' @param contrast 比较组
#' @param base_family 默认字体"sans"
#'
#' @return
#' @export 输出上下调个数图
#'
#' @examples
create_updownbar<-function(enrichment_final,order_value="Pvalue",enrich_type,
                           top_n=20,updown_color = c("#f5780c","#0586c0"),output_dir,contrast,base_family="sans"){
  enrichment_final<-enrichment_final%>%arrange(!!sym(order_value))%>%head(top_n)
  enrichment_final<-enrichment_final%>%arrange(DEP)
  levels_fc<-enrichment_final$Description
  plotdata<-enrichment_final%>%select(ID,Description,Up,Down)%>%pivot_longer(cols = c("Up","Down"))
  plotdata$Description<-factor(plotdata$Description,levels = levels_fc)
  p1<-ggplot(plotdata,aes(Description,value,fill=name))+
    geom_col(position = "dodge")+
    geom_text(aes(label = value), position = position_dodge(width = 0.8), hjust = -0.8,size = 3, color = "black") +
    scale_fill_manual(values= c("Up" = updown_color[1], "Down" = updown_color[2]))+
    theme_classic(base_family=base_family) +
    theme(axis.title.x = element_text(color = "black", size = 13),
          axis.title.y = element_text(color = "black", size = 13,vjust = 5),
          axis.text.y = element_text(color = "black", size = 10),
          axis.text.x = element_text(color = "black", size = 10),
          legend.title = element_blank(),
          plot.margin = margin(5,5,5,20),
    )+
    scale_y_continuous(limits = c(0,max(plotdata$value)+0.5))+
    coord_flip() +  # 旋转坐标轴
    labs(x = paste0(enrich_type," Description"), y = "Number of proteins")+
    guides(fill=guide_legend(reverse=T))
  ggsave(filename = paste0(output_dir,"/",contrast,"_",enrich_type,"_updownbar.pdf"),p1,width = 10,height =8,bg="white" )
}

#' Title 分类图（GO、KEGG）
#'
#' @param enrich_data 富集表格
#' @param enrich_type 富集类型
#' @param top_n 前多少默认20,KEGG 5
#' @param output 输出路径
#' @param contrast 比较组
#'
#' @return
#' @export
#'
#' @examples
create_classify_plot<-function(enrich_data,enrich_type,output,contrast){
  if(enrich_type == "GO"){
    top_n = 10
    go.bp <- enrich_data %>%
      filter(Category == "BP") %>%
      arrange(Pvalue)
    go.cc <- enrich_data %>%
      filter(Category == "CC") %>%
      arrange(Pvalue)
    go.mf <- enrich_data %>%
      filter(Category == "MF") %>%
      arrange(Pvalue)
    go.bp.top <- go.bp[1:min(top_n,nrow(go.bp)), ] %>% arrange(desc(DEP))
    go.cc.top <- go.cc[1:min(top_n,nrow(go.cc)), ] %>% arrange(desc(DEP))
    go.mf.top <- go.mf[1:min(top_n,nrow(go.mf)), ] %>% arrange(desc(DEP))
    go.bp.top$Description <- factor(go.bp.top$Description,levels = go.bp.top$Description)
    go.cc.top$Description <- factor(go.cc.top$Description,levels = go.cc.top$Description)
    go.mf.top$Description <- factor(go.mf.top$Description,levels = go.mf.top$Description)
    enrich.go <- rbind(go.bp.top,go.cc.top,go.mf.top)
    if(any(is.na(enrich.go$ID))){
      enrich.go<-enrich.go[-which(is.na(enrich.go$ID)),]
    }
    #修整描述的长度
    des_length = sapply(enrich.go$Description,  function(i){nchar(as.vector(i))})
    limit_length = ceiling(unname(quantile(des_length,0.8,na.rm = T)))
    enrich.go$Description <- factor(enrich.go$Description, level = enrich.go$Description)
    p <- ggplot(data = enrich.go, mapping = aes(x = Description, y = DEP, fill = Category)) +
      geom_bar(stat = 'identity') +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box = "horizontal",
        legend.position = 'right',
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20),
        axis.title = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.text.x = element_text(angle = 80, hjust = 1, color = "black", size = 15),
        strip.text = element_text(size = 15)
      ) +
      labs(x = 'GO Description', y = 'Number of proteins')
    if (!dir.exists(output)) {
      dir.create(output,recursive = T)
    }
    output<-paste0(output,"/",contrast)
    ggsave(paste0(output, '_GO_classify.pdf'), p, width  = 650,height = 400,units = 'mm', dpi = 300)
  }else{
    kegg.plot.data.top20<-enrich_data%>%group_by(level1)%>%
      arrange(desc(Pvalue),.by_group = T)%>%
      do(head(., n = 5))%>%
      arrange(DEP,.by_group = T)%>%
      select(ID,Description,level1,DEP)
    
    colourCount <- length(unique(kegg.plot.data.top20$level1))
    color_test <- c("#00468BFF", "#ED0000FF" ,"#42B540FF", "#0099B4FF", "#925E9FFF", "#FDAF91FF" ,"#AD002AFF")
    col <- color_test[1:colourCount]
    
    # Main bar plot
    p <- ggplot(kegg.plot.data.top20, aes(x = reorder(Description,DEP), y = DEP, fill = level1)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = col) +
      labs(x = "KEGG Description", y = "Number of proteins") +
      geom_text(aes(label = DEP), hjust = -0.3, size = 3) + # Adjust hjust for label position
      theme_bw(base_family = "sans") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(color = "black", size = 14),
            axis.text = element_text(color = "black", size = 10),
            axis.title.y = element_text(color = "black", size = 10,vjust=5),
            legend.position = "none", # Remove legend
            plot.margin = margin(t = 1, r = 0.1, b = 1, l = 1, unit = "cm"),
            strip.text = element_blank()) +
      facet_grid(level1 ~ ., scales = "free_y", space = "free") # Adjust scales and space
    
    # Annotations on the side
    df1 <- data.frame(x = "A", y = kegg.plot.data.top20$Description, group = kegg.plot.data.top20$level1)
    df1$group <- factor(df1$group, levels = unique(df1$group))
    
    p1 <- ggplot(df1, aes(x, y, fill = group)) +
      geom_tile(show.legend = FALSE) +
      facet_grid(group ~ ., scales = 'free_y', space = 'free') +
      labs(x = NULL, y = NULL) +
      #scale_y_discrete(expand = c(0,0))+
      theme(panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.margin = margin(0, -0.05, 0, 0, "cm"),
            strip.text.y = element_text(angle = 0, size = 10, color = "black", hjust = 0, margin = margin(b = 3, t = 3)),
            strip.background = element_rect(colour = NULL, fill = 'white')) +
      scale_fill_manual(values = col)
    
    # Combine plots
    p2 <- ggarrange2(p, p1, nrow = 1, widths = c(2, 0.05))
    if (!dir.exists(output)) {
      dir.create(output,recursive = T)
    }
    output<-paste0(output,"/",contrast)
    ggsave(paste0(output, "_KEGG_classify.pdf"),p2,width = 12,height = 7.2,dpi = 300,bg="white")
  }
    
}


#' Title 一个总函数。但有部分绘图参数没加。开放参数目的是为了售后好做。
#'
#' @param enrich_file 富集文件
#' @param diff_prot_file 差异文件
#' @param top_n 前多少默认20,KEGG 5
#' @param order_value 排序值 默认Pvalue
#' @param enrich_type 富集类型
#' @param output_dir 输出路径
#' @param contrast 比较组
#' @param p_color p值颜色c('#B2182B' ,'#2166AC')
#' @param base_family 基础字体默认是"sans"
#' @param fill_color 填充颜色
#'
#' @return
#' @export 输出同一个类型的6张图
#'
#' @examples
main<-function(enrich_file,diff_prot_file,top_n=20,order_value="Pvalue",enrich_type,output_dir,contrast,
               p_color=c('#B2182B' ,'#2166AC'),base_family="sans",fill_color){
  print(paste0(contrast,"_",enrich_type))
  
  dep<-read.xlsx(diff_prot_file)
  if(file.exists(enrich_file)){
    enrich_df<-read_enrich_data(enrich_file)
    if(nrow(enrich_df)>0){
      if(nrow(enrich_df)>1){
         network_plot(enrich_data = enrich_df,enrich_type=enrich_type,outdir=output_dir,contrast=contrast,geneSet="Protein")
      }else{
	 write.table(x = "Network diagrams require more than two entries to be plotted",file = paste0(output_dir,"/",contrast,"_网络图需要两个条目以上才进行绘制.txt"),row.names = F,col.names = F,quote = F)
      }
      create_richfactor_plot(enrich_data = enrich_df,enrich_type=enrich_type,output_dir=output_dir,contrast=contrast)
      creat_barplot(enrich_data = enrich_df,enrich_type=enrich_type,output_dir=output_dir,contrast=contrast,fill_color=fill_color)
      if(enrich_type == "GO" || enrich_type =="KEGG"){
        create_classify_plot(enrich_data = enrich_df,enrich_type=enrich_type,output=output_dir,contrast=contrast)
      }
      if("Up" %in% colnames(enrich_df)){
        sub_circlize_plot(enrichment_final=enrich_df,save_file=paste0(output_dir,"/",contrast,"_",enrich_type))
        chordfun(enrich_data=enrich_df,dep=dep,enrich_type=enrich_type,outdir=output_dir,contrast=contrast)
        create_updownbar(enrichment_final=enrich_df,enrich_type=enrich_type,output_dir=output_dir,contrast=contrast)
      }else{
        write.table(x = "circos, chord, updownbar plot need diff data",file = paste0(output_dir,"/",contrast,"_没有差异信息不绘制图像.txt"),row.names = F,col.names = F,quote = F)
      }
    }else{
      write.table(x="No enrichment results, cannot plot",file=paste0(output_dir,"/",contrast,"_无富集结果不能绘图.txt"),row.names = F,col.names = F,quote = F)
    }
    
    
  }else{
    write.table(x="No enrichment results, cannot plot",file=paste0(output_dir,"/",contrast,"_无富集结果不能绘图.txt"),row.names = F,col.names = F,quote = F)
  }
  
}

#' Title确定富集路径
#'
#' @param enrich_type 富集类型
#'
#' @return
#' @export 输出富集类型所对应的路径
#'
#' @examples
classify_outputdir<-function(enrich_type){
  if(enrich_type == "GO"){
    outdir<-"./result/5_Enrichment/5.1_GO/"
  }else if(enrich_type == "KEGG"){
    outdir<-"./result/5_Enrichment/5.2_KEGG/"
  }else if(enrich_type == "InterPro"){
    outdir<-"./result/5_Enrichment/5.3_domain/"
  }else if(enrich_type == "Reactome"){
    outdir<-"./result/5_Enrichment/5.4_Reactome/"
  }else if(enrich_type == "DO"){
    outdir<-"./result/5_Enrichment/5.5_DO/"
  }else if(enrich_type == "WikiPathway"){
    outdir<-"./result/5_Enrichment/5.6_WikiPathway/"
  }
  return(outdir)
}

source("./config/set_Theme.R")
source("./config/default.conf")
AddFont("config/myFonts")
showtext_auto()
fill_color<-get_color_palette("Summer")
no_cores <- 4
contrast_file<-"./temp/contrast.xlsx"
contrast_df<-read.xlsx(contrast_file)
sapply(contrast_df$contrast, function(contrast){
  library(parallel)
  start=Sys.time()
  cl<-makeCluster(no_cores,type="FORK")
  
  diff_prot_file<-paste0("./result/4_Diff_Expressed/4.1_DiffStats/",contrast,"-DEP_results.xlsx")
  parSapply(cl, c("GO","KEGG","InterPro","Reactome","DO","WikiPathway"),function(x) {
    showtext_auto()
    outdir<-classify_outputdir(enrich_type = x)
    outdir<-paste0(outdir,"/",contrast,"/")
    enrich_file<-paste0(outdir,"/",contrast,"_",x,"_enrichment.xlsx")
    print(enrich_file)
    main(enrich_file=enrich_file,diff_prot_file=diff_prot_file,enrich_type=x,output_dir=outdir,contrast=contrast,fill_color=fill_color)
  })
  stopCluster(cl)
  end<-Sys.time()
  print(end-start)
})




