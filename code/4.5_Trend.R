rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(ComplexHeatmap)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(circlize)))
# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/4_Diff_Expressed/4.1_DiffStats/"
total_proteins_path <- "./result/2_Data_Summary/Protein_SummaryStat.xlsx"
group_file <- "./temp/group.xlsx"
contrast_file <- "./temp/contrast.xlsx"
output <- "./result/4_Diff_Expressed/4.5_Trend/"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}


# 读取数据文件
# setwd("../../Kmeans/")
difp.file <- list.files(path = file_path, pattern = "-DEP_results.xlsx$")
all_proteins <- list()

group_names <- read.xlsx(group_file)
sample_names <- group_names$SampleName
Group_names <- unique(group_names$Group)
num_group <- length(Group_names)
# contrast_names <- read.xlsx(contrast_file)
# column_names <- contrast_names[,1]
# column_names_num <- nrow(contrast_names)

# 对比策略小于2，则停止运行脚本
if (num_group <= 2) {
  error_plot_heatmap <- ggplot() +
    theme_void() +
    geom_text(aes(0, 0, label = "num_group must be at least 3,the heatmap can't be plotted"), size = 10) +
    xlab(NULL)
  ggsave(paste0(output, "/error_info.png"), error_plot_heatmap, width = 15,height = 10,dpi = 300, bg = "white")
  
  #error_plot_kMeans <- ggplot() +
  #  theme_void() +
  #  geom_text(aes(0, 0, label = "num_group must be at least 3,the k_Means can't be plotted"), size = 10) +
  #  xlab(NULL)
  #ggsave(paste0(output, "/error_info.png"), error_plot_kMeans, width = 15,height = 10,dpi = 300, bg = "white")
  stop("num_group must be at least 3,  the trend can't be plotted, it's a sample problem, please ignore this error")
}


# 循环读取每个文件
for (i in seq_along(difp.file)) {
  # 读取文件
  # i = 1
  difp.data <- read_excel(file.path(file_path, difp.file[i]))
  # 更换表头，将 "Protein.Accessions" 替换为 "Protein"
  difp.data <- difp.data %>% dplyr::rename(Protein = Protein.Accessions, Genes = Gene.Name)
  # colnames(difp.data) <- gsub("Protein.Accessions", "Protein", colnames(difp.data))
  # colnames(difp.data) <- gsub("Gene.Name", "Genes", colnames(difp.data))
  # 确保需要删除的列存在
  # if (ncol(difp.data) >= 4) {
  #   difp.data <- difp.data[,-c(1,4)]
  # } else {
  #   stop(paste("File", difp.file[i], "does not have enough columns."))
  # }
  # 查找列名中包含 "Protein" 或 "Mean" 的列索引
  key.index <- key.index <- c(which(colnames(difp.data) %in% c("Protein")), grep(pattern = "Mean", colnames(difp.data)))
  # 选择包含这些关键字的列
  if (length(key.index) > 0) {
    difp.data.select <- difp.data[, key.index]
    all_proteins[[i]] <- difp.data.select$Protein
  } else {
    warning(paste("File", difp.file[i], "does not have columns with 'Protein' or 'Mean'."))
    all_proteins[[i]] <- NULL
  }
}

# 找到所有文件中的蛋白质名称的并集
union.prot <- Reduce(union, all_proteins)

#for (i in 1:length(difp.file)) {
  #difp.data <- read_excel(file = paste0("./result/4_Diff_Expressed/",difp.file[i]))
  #difp.data <- difp.data[,-c(1,4)]
  #key.index <- grep(pattern = "Protein|Mean",colnames(difp.data))
  #difp.data.select <- difp.data[,key.index]
  #all_proteins[[i]] <- difp.data.select$Protein
#}
#修改文件读取 此处为总表读取 不需要循环 20240627 by wjj
#total.prot.file <- list.files(path = "./input/",pattern = "-total_proteins.xlsx$")
#difp.data.merge <- list()
#for (j in 1:length(total.prot.file)) {
  #difp.data <- read.csv(file = paste0("./input/",total.prot.file[j]))
  #difp.data <- difp.data[,-c(1,4)]
  #key.index <- grep(pattern = "Protein|Mean",colnames(difp.data))
  #difp.data.select <- difp.data[,key.index]
  #difp.data.filter <- difp.data.select[difp.data.select$Protein %in% union.prot,]
  #difp.data.merge[[j]] <- difp.data.filter
#}
#merged_data <- do.call(cbind, difp.data.merge)

difp.data <- read_excel(total_proteins_path)
# 更换表头，将 "Protein.Accessions" 替换为 "Protein"
difp.data <- difp.data %>% dplyr::rename(Protein = Protein.Accessions, Genes = Gene.Name)
key.index <- key.index <- c(which(colnames(difp.data) %in% c("Protein")), grep(pattern = "Mean", colnames(difp.data)))
difp.data.select <- difp.data[,key.index]
difp.data.filter <- difp.data.select[difp.data.select$Protein %in% union.prot,]
merged_data <- difp.data.filter


duplicated_cols <- duplicated(colnames(merged_data))
# 保留唯一的列
unique_data <- merged_data[, !duplicated_cols]
colnames(unique_data) <- gsub(pattern = "Mean_",replacement = "",colnames(unique_data))
df <- unique_data
out_diff <- df
df = df[,-1]
dat = as.data.frame(t(apply(df, 1, scale)))
colnames(dat) = colnames(df)
df = dat
if(nrow(df) >= 150){
    center = 9
  }else{
    center = 7
  }
set.seed(1234)
fit <- kmeans(df, centers=center, iter.max=200, nstart=25)
withinss <- fit$tot.withinss
fit_cluster = fit$cluster 
df$cluster = fit$cluster
df$proteins = out_diff$Protein
df_new = melt(df,id.vars=c("cluster","proteins"))
cou = table(df$cluster)
for(i in 1:nrow(df_new)){
  df_new$name[i] = paste0("Cluster",df_new$cluster[i],", ",cou[df_new$cluster[i]]," proteins")
}
for(i in 1:nrow(df_new)){
  df_new$cluster[i] = paste0("Cluster",df_new$cluster[i])
}
# output <- "./output/"
#if (!dir.exists(output)) {
#  dir.create(output)
#}
write.xlsx(df,paste0(output,'kmeans.xlsx'))
## heatmap for cluster
cluster <- unique(df$cluster)
# 颜色
# 加载颜色方案
set_theme_module <-  "./config/set_Theme.R"
source(set_theme_module)
# print(group_colors)
if (group_colors %in% names(color.ls)) {
  group_colors <- color.ls[[group_colors]]
  #print(group_colors)
  #print(class(group_colors))
} else if (group_colors %in% names(color.continuous.ls)) {
  group_colors <- color.continuous.ls[[group_colors]]
}



for (i in 1:length(cluster)) {
  #i = 1
  cluster.prot <- df[df$cluster %in% cluster[i],c("proteins")]
  cluster.data.merge <- list()
  #为一个总表 无需再次循环读取
  difp.data <- read_excel(total_proteins_path)
  # 更换表头，将 "Protein.Accessions" 替换为 "Protein"
  difp.data <- difp.data %>% dplyr::rename(Protein = Protein.Accessions, Genes = Gene.Name)
  key.index <- c(which(colnames(difp.data) %in% c("Protein")), grep(pattern = "Mean", colnames(difp.data)))
  difp.data.select <- difp.data[,c(key.index[1]:(key.index[2]-1))]
  difp.data.filter <- difp.data.select[difp.data.select$Protein %in% cluster.prot,]
  merged_data <- difp.data.filter
  
  duplicated_cols <- duplicated(colnames(merged_data))
  
  merged_data$Protein_Gene <- paste0(merged_data$Protein, "_", merged_data$Genes)
  heat_data <- merged_data
  
  heat_data <- heat_data %>% select(all_of(sample_names))
  rownames(heat_data) <- merged_data$Protein_Gene
  
  #heat_data <- heat_data[,-1]
  # group <- names(heat_data)[grepl("^[A-Za-z]", names(heat_data))]
  # group <- sub("\\d+$", "", group)
  # group_counts <- table(group)
  # groupID <- unique(group)
  groupID <- unique(group_names$Group) 
  
  colors <- group_colors[(1:length(groupID))]
  #print(colors)
  #print(class(colors))
  group_col <- list(type = setNames(colors, groupID))
  logo_colors <- colorRampPalette(c(logo_blue, "white", logo_orange))(100)
  
  annotation_col = as.data.frame(group_names$Group)
  rownames(annotation_col) = colnames(heat_data)
  colnames(annotation_col) = "group"
  # heat_data  scale
  data_scale <- apply(heat_data, 1, scale)
  data_scale <- as.data.frame(data_scale)
  rownames(data_scale) <- colnames(heat_data)
  data_scale <- t(data_scale)
  data_scale <- as.data.frame(data_scale)
  rownames(data_scale) <- rownames(heat_data)
  

  annotation_top <- HeatmapAnnotation(df = annotation_col,
                                      name = "Group",
                                      col = list(group = group_col$type),
                                      #col = list(group = colors),
                                      simple_anno_size = unit(0.4, "cm"),
                                      show_legend = FALSE,
                                      show_annotation_name = FALSE
                                      )
  annotation_blank <- HeatmapAnnotation(
    blank = anno_empty(border = FALSE, height = unit(0.2, "cm")))
  
  protein_num <- length(data_scale[,1])
  

  heat_map <- Heatmap(data_scale,
                     name = "Z_Score",
                     #width =  if (protein_num > 50) unit(30, "cm") else unit(15, "cm"), 
                     width = unit(12, "cm"),
                     #height = if (protein_num > 50) unit(50, "cm") else unit(20, "cm"),
                     height = unit(15, "cm"),
                     
                     col = logo_colors,
                     # 注释
                     top_annotation = annotation_top,
                     bottom_annotation = annotation_blank,
                     # row
                     cluster_rows = TRUE,
                     clustering_distance_rows = "euclidean",
                     clustering_method_rows = "complete",
                     row_dend_side = "left",
                     row_dend_width = unit(25, "mm"), 
                     show_row_dend = TRUE,
                     show_row_names = FALSE,
                     row_title = NULL,
                     
                     # column
                     cluster_columns = TRUE,
                     clustering_method_columns = "complete",
                     clustering_distance_columns = "euclidean",
                     column_dend_side = "top",
                     column_dend_height = unit(15, "mm"),
                     show_column_dend = TRUE,
                     show_column_names = TRUE,
                     column_title = NULL,
                     column_names_rot = 45,
                     column_names_centered = TRUE,
                     
                     column_names_gp = gpar(fontsize = 10),
                     row_names_gp = gpar(fontsize = 10),
                     
                     use_raster = TRUE, # 图形渲染
                     
                     show_heatmap_legend = FALSE,
                     # heatmap_legend_param = heatmap_legend_param
                     )
  
  heat_map2 <- Heatmap(data_scale,
                     name = "Z_Score",
                     # width =  if (protein_num > 50) unit(30, "cm") else unit(15, "cm"), 
                     width = unit(12, "cm"),
                     # height = if (protein_num > 50) unit(50, "cm") else unit(20, "cm"),
                     height = unit(15, "cm"),
                     
                     col = logo_colors,
                     # 注释
                     top_annotation = annotation_top,
                     
                     # row
                     cluster_rows = TRUE,
                     clustering_distance_rows = "euclidean",
                     clustering_method_rows = "complete",
                     row_dend_side = "left",
                     row_dend_width = unit(15, "mm"), 
                     show_row_dend = TRUE,
                     show_row_names = FALSE,
                     row_title = NULL,
                     
                     # column
                     cluster_columns = FALSE,   # 不进行列聚类
                     clustering_method_columns = "complete",
                     clustering_distance_columns = "euclidean",
                     column_dend_side = "top",
                     column_dend_height = unit(15, "mm"),
                     show_column_dend = TRUE,
                     show_column_names = TRUE,
                     column_title = NULL,
                     column_names_rot = 45,
                     column_names_centered = TRUE,
                     
                     column_names_gp = gpar(fontsize = 10),
                     row_names_gp = gpar(fontsize = 10),
                     
                     use_raster = TRUE, # 图形渲染
                     
                     show_heatmap_legend = FALSE,
                     # heatmap_legend_param = heatmap_legend_param
                     )
  
  # 创建图例
  legend_annotation <- Legend(labels = groupID,
                            title = "group",
                            legend_gp = gpar(fill = c(colors[])),
                            title_gp = gpar(fontsize = 10),
                            labels_gp = gpar(fontsize = 8),
                            title_gap = unit(0.2, "cm"),
                            grid_width = unit(0.5, "cm"),
                            grid_height = unit(0.5, "cm")
)
  down <- round(min(apply(data_scale, 2, min)),1)
  top <- round(max(apply(data_scale, 2, max)),1)
  mean_c <- round(mean(apply(data_scale, 2, mean)),1)
  # top_middle <- round((top - mean_c)/2, 1)
  # down_middle <- round((down - mean_c)/2, 1)

  at_values <- c(down, mean_c, top)
  col_heatmap <- colorRamp2(breaks = at_values, colors = colorRampPalette(c(logo_blue, "white", logo_orange))(length(at_values)))
  legend_heatmap <- Legend(col_fun = col_heatmap,
                        #  at = at_values,
                         title = "Z_Score",
                         title_gp = gpar(fontsize = 10),
                         title_gap = unit(0.3, "cm"),
                         labels_gp = gpar(fontsize = 8),
                         legend_height = unit(1, "cm"),
                         legend_width = unit(1, "cm"),
                         grid_width = unit(0.5, "cm"),
                         grid_height = unit(0.5, "cm"),
                         tick_length = unit(1.5, "mm")
  )
  combined_legend <- packLegend(legend_annotation, legend_heatmap)

  heatmap_grob <- grid.grabExpr(draw(heat_map))
  heatmap_nocluster_grob <- grid.grabExpr(draw(heat_map2))
  
  
  
  dir <- paste0(output, "cluster",cluster[i])
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  png(paste0(dir,"/heatmap_cluster_sample.png"), width = 2300 ,height =  2600, res = 300 )
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(1, "npc")-unit(4, "cm"), unit(4, "cm")), heights = unit(1, "npc"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(heatmap_grob)
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  draw(combined_legend, x = unit(1, "npc")-unit(2, "cm") , y = unit(0.838, "npc"), just = c("right", "top"))
  popViewport()
  dev.off()
  
  
  pdf(paste0(dir,"/heatmap_cluster_sample.pdf"), width =  8, height = 8)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(1, "npc")-unit(4, "cm"), unit(4, "cm")), heights = unit(1, "npc"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(heatmap_grob)
  # 
  # 
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  draw(combined_legend, x = unit(1, "npc") + unit(1, "cm") , y = unit(0.868, "npc"), just = c("right", "top"))
  popViewport()
  dev.off()
  ## 列不聚类q
  png(paste0(dir,"/heatmap.png"), width =  2300, height =  2600, res = 300)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(1, "npc")-unit(4, "cm"), unit(4, "cm")), heights = unit(1, "npc"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(heatmap_nocluster_grob)
   
   
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  draw(combined_legend, x = unit(1, "npc") + unit(1, "cm") , y = unit(0.87, "npc"), just = c("right", "top"))
  popViewport()
  dev.off()
  
  pdf(paste0(dir,"/heatmap.pdf"), width = 8, height = 8)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(1, "npc")-unit(4, "cm"), unit(4, "cm")), heights = unit(1, "npc"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(heatmap_nocluster_grob)
  # 
  # 
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  draw(combined_legend, x = unit(1, "npc") + unit(1, "cm") , y = unit(0.90, "npc"), just = c("right", "top"))
  popViewport()
  dev.off()
}
#画图
# 设置颜色调色板
color_palette <- brewer.pal(n = center, name = "Set1")  # 或其他调色板如 "Dark2", "Paired"
p <- ggplot(df_new,aes(variable, value, group=proteins)) + 
  geom_line(aes(colour = cluster),size = 2) + 
    stat_summary(aes(group=1),fun.y = mean, geom="line", size = 0.5, color = "black") + 
  facet_wrap(~name,scales="free") +
  labs(x=" ",y="Standerised Intensity")+
  expand_limits(y=c(-2, 2))+
  scale_color_manual(values = color_palette) + 
  theme_minimal() + 
  theme(panel.grid=element_blank())+
  theme(
    strip.background = element_rect(
      color = "white", fill = "white"),
    panel.grid = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",
        axis.text = element_text(size=12, face = "bold"),
        axis.text.x = element_text(angle=0, hjust=1),
        strip.text = element_text(size = 8, face = "bold"),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(0.5, "in"),
        panel.spacing.y = unit(0.5, "in"))
ggsave(  # 输出位置
  paste0(
    output,  
    'K_Means_Plot.pdf'
  ), p,
  width = 9.6, height = 7.2, units = 'in', dpi = 600
)
ggsave(  # 输出位置
  paste0(
    output, 
    'K_Means_Plot.png'
  ), p,
  width = 9.6, height = 7.2, units = 'in', dpi = 300,bg = "white"
)

