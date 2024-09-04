rm(list = ls())


suppressWarnings(suppressMessages(library(UpSetR)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(writexl)))

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
column_names_num <- nrow(contrast_names)
#print(column_names_num)
# 对比策略小于2，则停止运行脚本
if (column_names_num < 2) {
  error_message <- "column_contrast must be at least 2, Because there is only one contrast, the upset can't be plotted, it's a sample problem, please ignore this error"
  writeLines(error_message, paste0(output, "upset不能绘制的原因.txt"))

  stop("column_names_num must be at least 2, Because there is only one contrast, the upset can't be plotted, it's a sample problem, please ignore this error")
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
  

# 统计差异表达结果
max_length <- max(sapply(data_list, nrow)) 
combined_df <- data.frame(matrix(ncol = length(column_names), nrow = max_length))
for (i in seq_along(selected_data_list)) {
  # 提取蛋白名称
  current.df <- unlist(selected_data_list[[i]])
  combined_df[,i] <- c(current.df, rep(NA, max_length - length(current.df)))
  
}
colnames(combined_df) <- column_names

# 所有蛋白的并集
all_protein_names <- unlist(combined_df)
all_protein_names <- all_protein_names[!is.na(all_protein_names)]
all_proteins <- unique(all_protein_names)
group_ID <- colnames(combined_df)

# 初始化0-1矩阵
presence_matrix <- matrix(0, nrow = length(all_proteins), ncol = length(group_ID), dimnames = list(all_proteins, group_ID))

# 遍历每个样本，更新0-1矩阵
for (i in seq_along(group_ID)) {
  current_proteins <- combined_df[!is.na(combined_df[group_ID[i]]),i]
  
  presence_matrix[current_proteins, i] <- 1
}
presence_matrix <- data.frame(Protein_Name = rownames(presence_matrix), presence_matrix)



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



# png
png(paste0(output,"upset_plot.png"), width = 2400, height = 2000 , res = 300)
upset(fromList(selected_data_list),
      nsets = length(selected_data_list),
      keep.order = TRUE,
      point.size = 3,
      line.size = 1.5,
      main.bar.color = logo_blue,
      matrix.color = logo_orange,
      sets.bar.color = logo_blue,
      mb.ratio = c(0.7, 0.3),
      shade.color = "#d3d3d3",
      text.scale = c(1.5, 1.3, 1.5, 1.3, 1.5, 1.3)
)
dev.off()


# pdf
pdf(paste0(output, "upset_plot.pdf"), width = 8, height = 6, onefile = FALSE)
upset(fromList(selected_data_list),
      nsets = length(selected_data_list),
      keep.order = TRUE,
      point.size = 3,
      line.size = 1,
      main.bar.color = logo_blue,
      matrix.color = logo_orange,
      sets.bar.color = logo_blue,
      mb.ratio = c(0.7, 0.3),
      shade.color = "#d3d3d3",
      text.scale = c(1.5, 1.3, 1.5, 1.3, 1.5, 1.3)
)
dev.off()

# 输出表格
# combined_proteins_replaced <- replace(combined_df, is.na(combined_df), "\t")
# write.csv(combined_proteins_replaced, file = paste0(output, "upset_table.csv"), row.names = FALSE)
write.xlsx(presence_matrix, file = paste0(output,"protein_table.xlsx"))
