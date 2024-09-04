rm(list = ls())
suppressWarnings(suppressMessages(library(tidyverse))) # 2.0.0
suppressWarnings(suppressMessages(library(data.table))) # 1.14.8
suppressWarnings(suppressMessages(library(igraph))) # 1.5.0
suppressWarnings(suppressMessages(library(ggraph))) # 2.1.0
suppressWarnings(suppressMessages(library(ggplot2))) # 3.4.3
suppressWarnings(suppressMessages(library(showtext)))
suppressWarnings(suppressMessages(library(ggnewscale)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
## 获取比较组信息
get_contrast <- function(file_name) {
  data <- read.xlsx(file_name)
  return(data)
}

## 获取输入文件: 5.2结果文件
get_data <- function(file_name) {
  data <- fread(file_name, sep = "\t")
  return(data)
}


## 获取蛋白互作网络图数据
get_net_data <- function(diffExp_df, PPI_df, Score_number, diff_gene_type) {
  if (diff_gene_type == "both") {
    diffExp_df <- diffExp_df## 只能包含上下调的基因
  }
  diffExp_df <- diffExp_df %>% dplyr::rename(Protein = Protein.Accessions,Genes = Gene.Name)
  colnames(diffExp_df)[grep(pattern = "Log2fc",x = colnames(diffExp_df))] <- "Log2fc"
  colnames(diffExp_df)[grep(pattern = "pvalue",x = colnames(diffExp_df))] <- "p_value"
  colnames(diffExp_df)[grep(pattern = "Q_value",x = colnames(diffExp_df))] <- "Q_value"
  colnames(diffExp_df)[grep(pattern = "Regulation",x = colnames(diffExp_df))] <- "Regulation"
  
  if(!("SwissProtName" %in% colnames(diffExp_df))) diffExp_df$SwissProtName <- diffExp_df$Protein
  
  # diffExp_df <- diffExp_df %>%
  #     dplyr::filter(Name != "-") %>%
  #     dplyr::filter(SwissProtName != "-")

  
  data_links <- as.data.frame(PPI_df[PPI_df$Node1 %in% diffExp_df$Genes & PPI_df$Node2 %in% diffExp_df$Genes & PPI_df$Score > Score_number,
                                     c("Node1","Node2","Score")])
  col_names <- names(data_links)
  #data_links <- data_links[c(col_names[1:2], col_names[10], col_names[3:9])]
  temp_data_links <- data_links
  match.index1 <- base::match(temp_data_links$Node1,diffExp_df$Genes)
  temp_data_links$Node1 <- diffExp_df$Protein[match.index1]
  match.index2 <- base::match(temp_data_links$Node2,diffExp_df$Genes)
  temp_data_links$Node2 <- diffExp_df$Protein[match.index2]
  # data_links <- left_join(data_links, diffExp_df[, c("id", "Name")], by = c("Node1" = "id")) ## id 和 name 的对应关系不具有唯一性 
  # data_links <- left_join(data_links, diffExp_df[, c("id", "Name")], by = c("Node2" = "id"))
  # data_links <- cbind(data_links[, c("Name.x", "Name.y")], data_links[, c(1:(ncol(data_links) - 2))])
  # temp_data_links <- data_links
  # data_links$Node1 <- NULL
  # data_links$Node2 <- NULL
  
  colnames(data_links)[1:3] <- c("fromN", "toN", "weight")
  
  # 节点数据
  nodes_list <- data_links %>%
    {
      data.frame(Genes = c(.$fromN, .$toN))
    } %>%
    distinct()
  diffExp_sub_df <- diffExp_df[which(diffExp_df$Genes %in% nodes_list$Genes), 
                               c("Genes", "Regulation", "p_value", "Q_value", "Log2fc", "Protein")] ## 需要增加 Gene_ID、SWissprotID
  nodes <- left_join(x = nodes_list, y = diffExp_sub_df, by = "Genes")
  
  nodes <- as.data.table(nodes)
  nodes$Regulation <- sub("Up", "up", nodes$Regulation)
  nodes$Regulation <- sub("Down", "down", nodes$Regulation)
  nodes$"-log10(Pvalue)" <- -log10(nodes$p_value)
  nodes$"-log10(adjustP)" <- -log10(nodes$Q_value)
  nodes$p_value <- NULL
  nodes$Q_value <- NULL
  
  # names(nodes)[c(1, 3:5)] <- c("Name", "log2FC", "Gene ID", "SWissprotID")
  
  # 创建网络图
  # 根据links和nodes创建
  nodes <- nodes[!duplicated(nodes[,1]),]
  net <- igraph::graph_from_data_frame(d = data_links, vertices = nodes, directed = F)
  ## 计算 Degree
  igraph::V(net)$Degree <- igraph::degree(net) # 每个节点连接的节点数
  # return(list(data_links,net))
  return(list(temp_data_links, net, data_links))
}

## 输出网络图的点和线的表格
save_file <- function(dataset_net, diffExp_df, contrast, suboutdir) {
  ## 线上线下均需要
  network <- dataset_net[[1]]
  network_gene <- dataset_net[[3]]
  colnames(network_gene) <- c("Node1","Node2","Score")
  network <- as.data.frame(network)[,c(1,2,3)]
  
  attributes <- as.data.table(vertex_attr(dataset_net[[2]])) %>% dplyr::rename(GeneName = "name")
  if (nrow(network) < 2) {
    print("差异蛋白之间没有相互作用，请调整PPI的阈值分数")
  }else{
    write.xlsx(as.data.table(network), paste0(suboutdir, "/", contrast, "_network_pro.xlsx"))
    write.xlsx(as.data.table(attributes), paste0(suboutdir, "/", contrast, "_attributes.xlsx"))
    write.xlsx(as.data.table(network_gene), paste0(suboutdir, "/", contrast, "_network_gene.xlsx"))
  }
}
main <- function(conf, offline = "null") {
  source(conf)
  ## 加载颜色方案
  set_theme_module <- paste0(utils_root, "/set_Theme.R")
  source(set_theme_module)
  ## 加载字体方案
  AddFont(font_path)
  showtext_auto()
  ## 加载图片保存方案
  source(paste0(utils_root, "/figure_Save.R"))
  ## 创建输出文件夹
  outdir <- "result/9_PPI"
  if (!dir.exists(outdir)) {
    dir.create(outdir,recursive = T)
  }
  contrast_df <- get_contrast(contrast.file)
  PPI_df <- get_data(PPI_file)
  ## 对每个比较组的文件进行循环处理
  for (contrast in contrast_df$contrast) {
    print(contrast)
    suboutdir <- paste0(outdir,"/",contrast)
    if (!dir.exists(suboutdir)) {
      dir.create(suboutdir)
    }
    pattern <- paste0("^", contrast, "-DEP_results.xlsx")
    diffExp_path <- dir(indir, pattern = pattern) %>% paste0(indir, .)
    ## 判断 5.2 的结果文件是否存在
    if (file.exists(diffExp_path)) {
      ## 修改标签值
      # link_color_label <- get_change_label(link_color_label)
      # node_color_label <- get_change_label(node_color_label)
      # node_size_label <- get_change_label(node_size_label)
      # link_width_label <- get_change_label(link_width_label)
      # font_label <- get_change_label(font_label)
      
      if (analysis) {
        # 获取 差异表达基因  文件: 5.2结果文件
        diffExp_df <- read.xlsx(diffExp_path)
        # 获取 物种对应PPI文件
        
        colnames(PPI_df)[c(1,2,10)] <- c("Node1","Node2","Score")
        ## 获取蛋白互作网络图数据
        dataset_net <- get_net_data(diffExp_df, PPI_df, Score_number, diff_gene_type)
        ## 输出网络图的点线属性表格
        save_file(dataset_net, diffExp_df, contrast, suboutdir)
      } else {
        load(RData)
        print("only update plot")
      }
      # get_net_plot(dataset_net, node_color_label, node_color, node_size_label, link_color_label, link_color, link_width_label, label_list, font_label, node_min, node_max, link_min, link_max, label_min, label_max, layout_algorithm, outdir, contrast, up_shape, down_shape, font_family, font_color, legend_size, legend_color, legend_family, up_color, down_color, node_size_label_fixed, link_width_label_fixed, font_label_fixed)
    } else {
      ## 容错信息
      err_info <- paste0(diffExp_path, " 文件不存在")
      fwrite(as.data.table(err_info), file = paste0(outdir, "/err.log"), sep = "\t", quote = F, col.names = FALSE)
    }
  }
}

# debug for ppi analysis
conf <- "config/default.conf"
source(conf)
PPI_file <- paste0(PPI_file_path, organism_code, "/ppi.txt")
main(conf)





