# 清除工作区
rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(UpSetR)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(tidyverse)))

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
  protein_names <- data[, 2]
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

# 重新筛选表格
sorted_data <- process_mat(data,group_names)

# 筛选分组列名,用来绘图
col_list <- list()
names_list <- list()
for (i in seq_along(group_ID)){
  sample_names <- group_file$SampleName[group_file$Group == group_ID[i]]
  result <- process_num(sorted_data,sample_names)
  sub_data <- sorted_data[,result]
  data_filter <- sub_data %>% filter(rowSums(is.na(.)) < ncol(sub_data))
  pro_names <- rownames(data_filter)
  col_list[[group_ID[i]]] <- result
  names_list[[group_ID[i]]] <- pro_names 
}

# combined_df用来输出表格
length_proteins <- nrow(sorted_data)
combined_df <- data.frame(matrix(ncol = length(group_ID), nrow = length_proteins))

# 构建蛋白名表
for (i in seq_along(group_ID)) {
    sub_data <- sorted_data[,col_list[[i]]]
    data_filter <- sub_data %>% filter(rowSums(is.na(.)) < ncol(sub_data))
    current.df <- rownames(data_filter)
    combined_df[,i] <- c(current.df, rep(NA, length_proteins - length(current.df)))
}

colnames(combined_df) <- group_ID



# 所有蛋白的并集
all_protein_names <- unlist(names_list)
all_protein_names <- all_protein_names[!is.na(all_protein_names)]
all_proteins <- unique(all_protein_names)

# 初始化0-1矩阵
presence_matrix <- matrix(0, nrow = length(all_proteins), ncol = length(group_ID), dimnames = list(all_proteins, group_ID))

# 遍历每个样本，更新0-1矩阵
for (i in seq_along(group_ID)) {
  current_proteins <- combined_df[!is.na(combined_df[group_ID[i]]),i]
  
  presence_matrix[current_proteins, i] <- 1
}
presence_matrix <- data.frame(Protein_Name = rownames(presence_matrix), presence_matrix)



# png
png(paste0(output,"upset_plot.png"), width = 2600, height = 2200, res = 300)
upset(fromList(names_list),
      nsets = length(names_list),
      keep.order = TRUE,
      point.size = 3,
      line.size = 1.5,
      main.bar.color = logo_blue,
      matrix.color = logo_orange,
      sets.bar.color = logo_blue,
      mb.ratio = c(0.7, 0.3),
      shade.color = "#d3d3d3",
      text.scale = c(1.5, 1.3, 1.5, 1.3, 1.5, 1.3),
      order.by = c("freq")
      )
dev.off()


# pdf
pdf(paste0(output,"upset_plot.pdf"), width = 8, height = 6, onefile = FALSE)
upset(fromList(names_list),
      nsets = length(names_list),
      keep.order = TRUE,
      point.size = 3,
      line.size = 1,
      main.bar.color = logo_blue,
      matrix.color = logo_orange,
      sets.bar.color = logo_blue,
      mb.ratio = c(0.7, 0.3),
      shade.color = "#d3d3d3",
      text.scale = c(1.5, 1.3, 1.5, 1.3,  1.5, 1.3),
      order.by = c("freq")
)
dev.off()

# 输出表格
# combined_proteins_replaced <- replace(combined_df, is.na(combined_df), "\t")
# write.csv(combined_proteins_replaced, file = paste0(output,"upset_table.csv"), row.names = FALSE)
 
write.xlsx(presence_matrix, file = paste0(output,"protein_table.xlsx"))
