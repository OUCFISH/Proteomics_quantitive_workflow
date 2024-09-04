# 清除工作区
rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library("RColorBrewer")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("ComplexHeatmap")))
suppressWarnings(suppressMessages(library("circlize")))
suppressWarnings(suppressMessages(library("openxlsx")))

# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/2_Data_Summary/Protein_SummaryStat.xlsx"
group_file <- "./temp/group.xlsx"
output <- "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.7_Global_Heatmap/"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# 按分组表格进行列名筛选
process_mat <- function(mat, map) {
  # 获取SampleName列对应的索引顺序
  new_order <- match(map, colnames(mat))
  # 按照新的顺序重新排列矩阵的列
  mat <- mat[, new_order]
  return(mat)
}

# 读取数据
data_total <- read.xlsx(file_path)
group_file <- read.xlsx(group_file)
data_sample <- process_mat(data_total,group_file$SampleName)
# data_sample <- data_total[ ,(9:(9+length(group_file$Group)-1))]
rownames(data_sample) <- data_total$Protein.Accessions


# 分组
group_ids <- group_file$Group
group_names <- group_file$SampleName
group_ID <- unique(group_ids)
#colnames(data_sample) <- group_ids


# 加载颜色方案
set_theme_module <-  "./config/set_Theme.R"
source(set_theme_module)

if (group_colors %in% names(color.ls)) {
  group_colors <- color.ls[[group_colors]]
} else if (group_colors %in% names(color.continuous.ls)) {
  group_colors <- color.continuous.ls[[group_colors]]
}

if (sample_colors %in% names(color.ls)) {
  sample_colors <- color.ls[[sample_colors]]
} else if (sample_colors %in% names(color.continuous.ls)) {
  sample_colors <- color.continuous.ls[[sample_colors]]
}

colors <- group_colors[1:length(group_ID)]
# 颜色
#base_colors <- c("#187679","#E6382E","#FFD119","#C27F35","#FE685F","#887A81","#CC525C","#FA9A4B","#44BBA4","#A7A09F","#5C9966","#86BF8B","#418E91","#E7716A","#F3D455","#C78E4F","#FE685F","#A8959F","#CC666E","#FAB276","#5BBDAA","#BDB3B1","#6DB97A","#A0D6A5","#5C6160","#A14C47","#EE6B27","#F69E20","#EAB522","#D69A2B","#D67743","#EA6F51","#D66E6A","#AF7475","#9E6C74","#B55F68","#DB6A56","#EA8250","#BDA568","#80B086","#65B2A2","#86A9A0","#8E9D8C","#759B79","#6AA572","#78B27E","#6FAE8D","#589E8F","#788484","#AF7A77","#EB9263","#EFB35C","#E4BC53","#D5A551","#D98154","#EB7459","#E17774","#C48689","#B4858E","#C0757E","#DB7F70","#EA9873","#C5B587","#90B998","#7BB9AC","#9CB6AE","#A2B59E","#87B78C","#7EC288","#8FCC96")
# group_colors <- base_colors[1:length(group_ID)]
colors = setNames(colors, group_ID)
mycol1 <- colorRampPalette(c(logo_blue, "white", logo_orange))(100)

# 绘制热图
# "ComplexHeatmap"
# data scale
data_scale <- apply(data_sample, 1, scale)
rownames(data_scale) <- colnames(data_sample)
data_scale <- t(data_scale)
data_scale <- as.data.frame(data_scale)

group <- data.frame(group = group_ids)
rownames(group) = colnames(data_scale)

annotation_top <- HeatmapAnnotation(df = group,
                                    name = "Group",
                                    col = list(group = colors),
                                    simple_anno_size = unit(0.3, "cm"),
                                    show_legend = FALSE,
                                    show_annotation_name = FALSE
                                    )
annotation_blank <- HeatmapAnnotation(
  blank = anno_empty(border = FALSE, height = unit(0.1, "cm")))                                    

# heatmap_legend_param = list(
#   title_position = "topcenter",
#   title_gp = gpar(fontsize = 12), 
#   label_gp = gpar(fontsize = 10), 
#   direction = "vertical", # 图例的方向，水平或垂直
#   
#   legend_width = unit(3, "cm"), # 图例宽度
#   legend_height = unit(1, "cm"), # 图例高度
#   legend_margin = unit(c(1, 1, 1, 1), "line"), # 图例边缘的空白
#   legend_space = unit(0.5, "cm"), # 图例项之间的间隔
#   # 通过调整x和y，你可以精细控制图例的位置
#   legend_pos = list("right", unit(12, "npc")), 
#   just = "top")

heatmap_cluster <- Heatmap(data_scale,
                   name = "Z_Score",
                   width = unit(12, "cm"),
                   height = unit(15, "cm"),

                   
                   col = mycol1,
                   # 注释
                   top_annotation = annotation_top,
                   bottom_annotation = annotation_blank,

                   # row
                   cluster_rows = TRUE,
                   clustering_distance_rows = "euclidean",
                   clustering_method_rows = "complete",
                   row_dend_side = "left",
                   row_dend_width = unit(20, "mm"), 
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

heatmap_no_cluster <- Heatmap(data_scale,
                   name = "Z_Score",
                   width = unit(12, "cm"),
                   height = unit(15, "cm"),

                   
                   col = mycol1,
                   # 注释
                   top_annotation = annotation_top,
                   bottom_annotation = annotation_blank,

                   # row
                   cluster_rows = TRUE,
                   clustering_distance_rows = "euclidean",
                   clustering_method_rows = "complete",
                   row_dend_side = "left",
                   row_dend_width = unit(20, "mm"), 
                   show_row_dend = TRUE,
                   show_row_names = FALSE,
                   row_title = NULL,
                   
                   # column
                   cluster_columns = FALSE,
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

#创建图例
legend_annotation <- Legend(labels = group_ID,
                            title = "group",
                            legend_gp = gpar(fill = c(colors[])),
                            title_gp = gpar(fontsize = 12),
                            labels_gp = gpar(fontsize = 10),
                            title_gap = unit(0.2, "cm"),
                            grid_width = unit(0.5, "cm"),
                            grid_height = unit(0.5, "cm")
)
# 图例断点
down <- round(min(apply(data_scale, 2, min)),1)
top <- round(max(apply(data_scale, 2, max)),1)
mean <- round(mean(apply(data_scale, 2, mean)),1)
top_middle <- round((top - mean)/2, 1)
down_middle <- round((down - mean)/2, 1)

at_values <- c(down, down_middle, mean, top_middle, top)
# at_values <- c(-4, -2, 0, 2, 4)
col_heatmap <- colorRamp2(breaks = at_values, colors = colorRampPalette(c(logo_blue, "white", logo_orange))(length(at_values)))
legend_heatmap <- Legend(col_fun = col_heatmap,
                         at = at_values,
                         title = "Z_Score",
                         title_gp = gpar(fontsize = 12),
                         title_gap = unit(0.3, "cm"),
                         labels_gp = gpar(fontsize = 10),
                         legend_height = unit(1, "cm"),
                         legend_width = unit(1, "cm"),
                         grid_width = unit(0.5, "cm"),
                         grid_height = unit(0.5, "cm"),
                         tick_length = unit(1.5, "mm")
)
combined_legend <- packLegend(legend_annotation, legend_heatmap)

heatmap_no_cluster_grob <- grid.grabExpr(draw(heatmap_no_cluster))
heatmap_cluster_grob <- grid.grabExpr(draw(heatmap_cluster))



#保存图片
# cluster
pdf(paste0(output, "heatmap_cluster.pdf"), width = 10, height = 9)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(1, "npc")-unit(4, "cm"), unit(4, "cm")), heights = unit(1, "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(heatmap_cluster_grob)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(combined_legend, x = unit(1, "npc")-unit(5, "cm") , y = unit(0.81, "npc"), just = c("right", "top"))
popViewport()
dev.off()
png(paste0(output, "heatmap_cluster.png"), width = 2400, height = 2200, res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(1, "npc")-unit(4, "cm"), unit(4, "cm")), heights = unit(1, "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(heatmap_cluster_grob)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(combined_legend, x = unit(1, "npc")-unit(2.8, "cm") , y = unit(0.88, "npc"), just = c("right", "top"))
popViewport()
dev.off()


# no _cluster
pdf(paste0(output, "heatmap_no_cluster.pdf"), width = 10, height = 9)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(1, "npc")-unit(4, "cm"), unit(4, "cm")), heights = unit(1, "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(heatmap_no_cluster_grob)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(combined_legend, x = unit(1, "npc")-unit(5, "cm") , y = unit(0.84, "npc"), just = c("right", "top"))
popViewport()
dev.off()
png(paste0(output, "heatmap_no_cluster.png"), width = 2400, height = 2200, res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(1, "npc")-unit(4, "cm"), unit(4, "cm")), heights = unit(1, "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(heatmap_no_cluster_grob)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(combined_legend, x = unit(1, "npc")-unit(2.8, "cm") , y = unit(0.93, "npc"), just = c("right", "top"))
popViewport()
dev.off()

