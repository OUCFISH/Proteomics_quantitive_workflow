rm(list = ls())

# 加载库
library(ggvenn)
library(grid)
library(ggplot2)
library(openxlsx)


conf <- "./config/default.conf"
source(conf)
file_path <- "./result/4_Diff_Expressed/4.1_DiffStats/"
group_file <- "./temp/group.xlsx"
contrast_file <- "./temp/contrast.xlsx"
output <- "./result/4_Diff_Expressed/4.2_Diff_Repeat/"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}


# 读取数据
group_names <- read.xlsx(group_file)
contrast_names <- read.xlsx(contrast_file)
column_names <- contrast_names[,1]
# 对比策略数目必须在[2:4]， 否则停止运行脚本
column_names_num <- nrow(contrast_names)
#print(column_names_num)
if (column_names_num < 2 | column_names_num > 4) {
  error_message <- "column_names_num must be [2:4], cause contrast is not right, the venn can't be plotted, it's a sample problem, please ignore this error"
  writeLines(error_message, paste0(output, "Venn不能绘制的原因.txt"))
  stop("column_names_num must be [2:4], cause contrast is not right, the venn can't be plotted, it's a sample problem, please ignore this error")
}

# 匹配差异表达结果文件
patterns <- paste0(column_names, "-DEP_results.xlsx")
# print(patterns)
excel_files <- paste0(file_path, patterns)
# print(excel_files)
data_list <- lapply(excel_files, read.xlsx)
selected_column_name <- "Protein.Accessions"
selected_data_list <- list()
for (i in seq_along(column_names)) {
  pro_names <- select(data_list[[i]], selected_column_name)
  names <- pro_names$Protein.Accessions
  selected_data_list[[column_names[i]]] <- names
}


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
colors <- group_colors[1:length(column_names)] 

# 绘制Venn图
#combined_proteins_list <- as.list(combined_df)
vene <- ggvenn(selected_data_list,
               show_percentage = FALSE,
               fill_color = colors,
               stroke_color = "white",
               stroke_size = 1,
               set_name_size = 8,
               text_size = 8,
               
               )

# 保存为PNG文件
ggsave(paste0(output, "venn.png"), vene, width = 12, height = 10, dpi = 300)

# 保存为PDF文件
ggsave(paste0(output, "venn.pdf"), vene, width = 12, height = 10)

# 输出表格
# combined_proteins_replaced <- replace(combined_df, is.na(combined_df), "\t")
# write.csv(combined_proteins_replaced, file = paste0(output, "venn_table.csv"), row.names = FALSE)

