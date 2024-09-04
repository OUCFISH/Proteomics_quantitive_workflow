# 清除工作区
rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(ggvenn)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(openxlsx)))

# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/2_Data_Summary/ProteinGroup.xlsx"
group_file <- "./temp/group.xlsx"
output <- "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.2_Expression_Repeat/"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# 读取数据
data <- read.xlsx(file_path)
if ("Protein.Group" %in% colnames(data)) {
  protein_names <- data$Protein.Group %>% strsplit(split = ";") %>% sapply(`[`, 1)
  rownames(data) <- protein_names
}else{
  protein_names <- data[,2]
  rownames(data) <- protein_names
}
# 分组
group_file <- read.xlsx(group_file)  
group_ids <- group_file$Group
group_names <- group_file$SampleName
group_ID <- unique(group_ids)

# 按照分组顺序重新排序
process_mat <- function(mat, map) {
  # 获取SampleName列对应的索引顺序
  new_order <- match(map, colnames(mat))
  # 按照新的顺序重新排列矩阵的列
  mat <- mat[, new_order]
  return(mat)
}

# 分组匹配的列号
process_num <- function(data, names) {
  # 获取SampleName列对应的列索引
  col_num <- match(names, colnames(data))
  return(col_num)
}

sorted_data <- process_mat(data,group_names)

# 筛选分组列名,用来绘图
names_list <- list()
for (i in seq_along(group_ID)){
  sample_names <- group_file$SampleName[group_file$Group == group_ID[i]]
  result <- process_num(sorted_data,sample_names)
  sub_data <- sorted_data[,result]
  data_filter <- sub_data %>% filter(rowSums(is.na(.)) < ncol(sub_data))
  pro_names <- rownames(data_filter)
  names_list[[group_ID[i]]] <- pro_names 
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


# 绘制venn图
#combined_proteins_list <- as.list(combined_df)
colors <- group_colors[1:length(group_ID)] 
# print(colors)
vene_with_group <- ggvenn(names_list,
              show_percentage = FALSE,
              fill_color = colors,
              stroke_color = "white",
              stroke_size = 1,
              set_name_size = 8,
              text_size = 8,
)

# print(vene_with_group)

# 保存为PNG文件
ggsave(paste0(output, "venn_group.png"), vene_with_group, width = 12, height = 10, dpi = 300)

# 保存为PDF文件
ggsave(paste0(output, "venn_group.pdf"), vene_with_group, width = 12, height = 10)


# 输出表格
# combined_proteins_replaced <- replace(combined_df, is.na(combined_df), "\t")
# write.csv(combined_proteins_replaced, file = paste0(output, "venn_table.csv"), row.names = FALSE)

