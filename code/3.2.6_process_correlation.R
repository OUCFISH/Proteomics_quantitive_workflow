# 清除工作区
rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(tidyverse)))


# 文件路径和参数
# conf <- "./config/default.conf"
# source(conf)
file_path <- "./result/2_Data_Summary/Protein_SummaryStat.xlsx"
group_file <- "./temp/group.xlsx"
output <- "./temp/"
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
group_file <- read.xlsx(group_file)
input_data <- read.xlsx(file_path)
data1 <- input_data$Protein.Accessions
data2 <- process_mat(input_data,group_file$SampleName)
# data2 <- input_data[ ,(9:(9+length(group_file$Group)-1))]
data <- cbind(data1,data2)
colnames(data) <- c("Protein", colnames(data2))
colnames(group_file) <- c("SampleID","SampleName","Default") 

# 保存文件
write.table(data, file = paste0(output, "fpkm.tsv"), row.names = FALSE, sep = "\t")
write.table(group_file, file = paste0(output, "group.tsv"), row.names = FALSE, sep = "\t")

# 后续分析用到的文件
colnames(data) <- c("Protein.Accessions", colnames(data2))
write.table(data, file = paste0(output, "fpkm_A.tsv"), row.names = FALSE, sep = "\t")
