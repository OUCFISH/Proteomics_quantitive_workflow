rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ComplexHeatmap)))
suppressWarnings(suppressMessages(library(circlize)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(grid)))


# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/4_Diff_Expressed/4.1_DiffStats/"
group_file <- "./temp/group.xlsx"
contrast_file <- "./temp/contrast.xlsx"
output <- "./result/4_Diff_Expressed/4.4_Heatmap/"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}



# 读取数据
group_names <- read.xlsx(group_file)
contrast_names <- read.xlsx(contrast_file)
column_names <- contrast_names[,1]

# 匹配差异表达结果文件
patterns <- paste0(column_names, "-DEP_results.xlsx")
excel_files <- paste0(file_path, patterns)
data_list <- lapply(excel_files, read.xlsx)

sample_names <- group_names$SampleName
#print(sample_names)
selected_data_list <- lapply(data_list, function(df) {
  # 选择固定的列
  fixed_cols <- c("Protein.Accessions", "Gene.Name")
  #print(fixed_cols)
  # 匹配以"pvalue"结尾的列名
  pvalue_cols <- grep("pvalue$", names(df), value = TRUE)
  #print(pvalue_cols)
  # 合并所有需要选择的列名
  cols_to_select <- c(fixed_cols, sample_names, pvalue_cols)
  #print(cols_to_select)
  # 使用select函数选择列
  selected_df <- df %>% select(any_of(cols_to_select))
  #print(head(selected_df))
  return(selected_df)
})
selected_data_list <- setNames(selected_data_list, column_names)

# 总分组
total_group <- group_names
group_ids_total <- total_group$Group
group_ID_total <- unique(group_ids_total)


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

colors <- group_colors[1:length(group_ID_total)] 
group_col <- list(type = setNames(colors, group_ID_total))
logo_colors <- colorRampPalette(c(logo_blue, "white", logo_orange))(100)
# print(group_col)


# 循环绘制不同对比策略下的heatmap
for (contrast_item in contrast_names$contrast) {
  # print(contrast_item)

if (!dir.exists(paste0(output, contrast_item))) {
  dir.create(paste0(output, contrast_item), recursive = TRUE)
}


# 筛选表格内容
data <- selected_data_list[[contrast_item]]
data$names <- paste(data$Protein.Accessions, data$Gene.Name, sep = "_")
rownames(data) <- data$names
data <- data[,-c(1:2,(ncol(data)-1):ncol(data))]



# 筛选p值top50
data_sorted <- selected_data_list[[contrast_item]]
data_sorted$names <- paste(data_sorted$Protein.Accessions, data_sorted$Gene.Name, sep = "_")
rownames(data_sorted) <- data_sorted$names
data_sorted <- data_sorted[,-c(1:2,ncol(data_sorted))]
pvalue_cols <- grep("pvalue$", names(data_sorted), value = TRUE)
data_sorted <- data_sorted[order(data_sorted[,pvalue_cols], decreasing = TRUE), ]
num_top <- nrow(data_sorted)
data_sorted_50 <- head(data_sorted,if (num_top > 50) 50 else num_top,)
data_sorted_50 <- select(data_sorted_50, -all_of(pvalue_cols))


# 分组
filtered_group <- group_names[group_names$SampleName %in% colnames(data), ]
group_ids <- filtered_group$Group
group_ID <- unique(group_ids)


# 选颜色

colors <- group_col$type[group_ID[]]


# "ComplexHeatmap"
# datascale
data_scale <- apply(data, 1, scale)
rownames(data_scale) <- colnames(data)
data_scale <- t(data_scale)
data_scale <- as.data.frame(data_scale)

group_ids <- data.frame(group = group_ids)
rownames(group_ids) = colnames(data_scale)
# print(group_ids)

data_sorted_scale <- apply(data_sorted_50, 1, scale)
rownames(data_sorted_scale) <- colnames(data_sorted_50)
data_sorted_scale <- t(data_sorted_scale)
data_sorted_scale <- as.data.frame(data_sorted_scale)


#绘制heatmap

annotation_top <- HeatmapAnnotation(df = group_ids,
                                    name = "group",
                                    col = list(group = colors),
                                    simple_anno_size = unit(0.8, "cm"),
                                    show_legend = FALSE,
                                    show_annotation_name = FALSE
                                    )

annotation_top_50 <- HeatmapAnnotation(df = group_ids,
                                    name = "group",
                                    col = list(group = colors),
                                    simple_anno_size = unit(0.8, "cm"),
                                    show_annotation_name = FALSE,
                                    show_legend = FALSE                                  
)
annotation_blank <- HeatmapAnnotation(
  blank = anno_empty(border = FALSE, height = unit(0.5, "cm")))

heatmap <- Heatmap(data_scale,
                   name = "Z_Score",
                   width = unit(40, "cm"),
                   height = unit(50, "cm"),
                   
                   
                   col = logo_colors,
                   # 注释
                   top_annotation = annotation_top,
                   bottom_annotation = annotation_blank,
                   
                   # row
                   cluster_rows = TRUE,
                   clustering_distance_rows = "euclidean",
                   clustering_method_rows = "complete",
                   row_dend_side = "left",
                   row_dend_width = unit(30, "mm"), 
                   show_row_dend = TRUE,
                   show_row_names = FALSE,
                   row_title = NULL,
                   
                   # column
                   cluster_columns = TRUE,
                   clustering_method_columns = "complete",
                   clustering_distance_columns = "euclidean",
                   column_dend_side = "top",
                   column_dend_height = unit(30, "mm"),
                   show_column_dend = TRUE,
                   show_column_names = TRUE,
                   column_names_rot = 45,
                   column_names_centered = TRUE,                   
                  
                   
                   column_names_gp = gpar(fontsize = 30),
                   row_names_gp = gpar(fontsize = 20),
                   
                   use_raster = TRUE, # 图形渲染
                   
                   show_heatmap_legend = FALSE,
                   # heatmap_legend_param = heatmap_legend_param
)

heatmap_50 <- Heatmap(data_sorted_scale,
                   name = "Z_Score",
                   width = unit(40, "cm"),
                   height = unit(50, "cm"),
                   
                   
                   col = logo_colors,
                   # 注释
                   top_annotation = annotation_top_50,
                   bottom_annotation = annotation_blank,
                   # row
                   cluster_rows = TRUE,
                   clustering_distance_rows = "euclidean",
                   clustering_method_rows = "complete",
                   row_dend_side = "left",
                   row_dend_width = unit(30, "mm"), 
                   show_row_dend = TRUE,
                   show_row_names =TRUE,
                   row_title = NULL,
                   
                   # column
                   cluster_columns = TRUE,
                   clustering_method_columns = "complete",
                   clustering_distance_columns = "euclidean",
                   column_dend_side = "top",
                   column_dend_height = unit(30, "mm"),
                   show_column_dend = TRUE,
                   show_column_names = TRUE,                   
                   column_names_rot = 45,
                   column_names_centered = TRUE,
                   
                   column_names_gp = gpar(fontsize = 30),
                   row_names_gp = gpar(fontsize = 20),
                   
                   use_raster = TRUE, # 图形渲染
                   
                   show_heatmap_legend = FALSE,
                   
)

# 创建图例
#col_group <- colorRamp2(breaks = c(group_ID), colors = colors[])
legend_annotation <- Legend(labels = c(group_ID),
                            title = "group",
                            legend_gp = gpar(fill = c(colors[])),
                            labels_gp = gpar(fontsize = 22),
                            title_gp = gpar(fontsize = 24),
                            title_gap = unit(0.5, "cm"),
                            legend_height = unit(8, "cm"),
                            legend_width = unit(8, "cm"),
                            grid_width = unit(1, "cm"),
                            grid_height = unit(1, "cm")
)
down <- round(min(apply(data_scale, 2, min)),1)
top <- round(max(apply(data_scale, 2, max)),1)
mean_c <- round(mean(apply(data_scale, 2, mean)),1)
at_values <- c(down, mean_c, top)
  
# at_values <- c(round(min(apply(data_scale, 2, min)), 4) + 0.0001, round(mean(apply(data_scale, 2, mean)), 4), round(max(apply(data_scale, 2, max))), 4)
# at_values <- c(-2, -1, 0, 1, 2)
col_heatmap <- colorRamp2(breaks = at_values, colors = colorRampPalette(c(logo_blue, "white", logo_orange))(length(at_values)))
legend_heatmap <- Legend(col_fun = col_heatmap,
                        #  at = at_values,       
                         title = "Z_Score",
                         title_gp = gpar(fontsize = 22),
                         title_gap = unit(0.8, "cm"),
                         labels_gp = gpar(fontsize = 22),
                         legend_height = unit(8, "cm"),
                         legend_width = unit(8, "cm"),
                         grid_width = unit(1, "cm"),
                         grid_height = unit(1, "cm"),
                         tick_length = unit(1.5, "mm")
                         
)


combined_legend <- packLegend(legend_annotation, legend_heatmap)

heatmap_50_grob <- grid.grabExpr(draw(heatmap_50))
heatmap_grob <- grid.grabExpr(draw(heatmap))

#保存图片
pdf(paste0(output, contrast_item,"/", "heatmap.pdf"), width = 24, height = 26)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(1, "npc")-unit(4, "cm"), unit(4, "cm")), heights = unit(1, "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(heatmap_grob)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(combined_legend, x = unit(1, "npc")-unit(6, "cm") , y = unit(0.88, "npc"), just = c("right", "top"))
popViewport()
dev.off()


png(paste0(output, contrast_item,"/", "heatmap.png"), width = 7400, height = 7600, res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(1, "npc")-unit(4, "cm"), unit(4, "cm")), heights = unit(1, "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(heatmap_grob)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(combined_legend, x = unit(1, "npc")-unit(7, "cm") , y = unit(0.88, "npc"), just = c("right", "top"))
popViewport()
dev.off()

pdf(paste0(output, contrast_item,"/", "heatmap_top50.pdf"), width = 24, height = 26)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(1, "npc")-unit(4, "cm"), unit(4, "cm")), heights = unit(1, "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(heatmap_50_grob)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(combined_legend, x = unit(1, "npc")-unit(2, "cm") , y = unit(0.88, "npc"), just = c("right", "top"))
popViewport()
dev.off()


png(paste0(output, contrast_item,"/", "heatmap_top50.png"), width = 8000, height = 8200, res = 300)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(1, "npc")-unit(4, "cm"), unit(4, "cm")), heights = unit(1, "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(heatmap_50_grob)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(combined_legend, x = unit(1, "npc")-unit(2.5, "cm") , y = unit(0.86, "npc"), just = c("right", "top"))
popViewport()
dev.off()

}





