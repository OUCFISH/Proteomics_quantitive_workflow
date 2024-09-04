# 清除工作区
rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(tidyverse)))


# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/2_Data_Summary/Protein_SummaryStat.xlsx"
group_file <- "./temp/group.xlsx"
output <- "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.3_Expression_Violin/"
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
input.data <- read.xlsx(file_path)
input.data <- process_mat(input.data,group_file$SampleName)

# input.data <- input.data[ ,(9:(8+length(group_file$Group)))]

# 转换为长格式
long_data <- pivot_longer(input.data, cols = everything(), names_to = "Variable", values_to = "Value")
long_data$Value <- log2(long_data$Value) 
#long_data$group <- sapply(strsplit(x = long_data$Variable,split = "(?<=\\D)(?=\\d)"),FUN = function(x)x[1])

# 分组
group_ids <- group_file$Group
group_names <- group_file$SampleName
group_ID <- unique(group_ids)
# print(group_ID)
long_data$group <- rep(group_ids, length.out = nrow(long_data))

# 加载颜色
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
colors = setNames(colors, group_ID)


# print(long_data)
# 绘制小提琴图
p <- ggplot(long_data, aes(x = Variable, y = Value, fill = group)) +
  geom_violin(color = "black") +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  labs(title = NULL,
       x = NULL, y = expression(Log[2]~"(Intensity)"), fill = "group") +
  # scale_fill_discrete(name = "Group") +
  scale_fill_manual(values = colors, name = "group") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(color = "black", size = 18),
    axis.text = element_text(color = "black", size = 12),
    axis.ticks = element_line(linewidth = 0.5),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.3, "cm"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
  ) 



ggsave(paste0(output, "Expression_Violin.png"),p,width = 12,height = 7.2,dpi = 300)
ggsave(paste0(output, "Expression_Violin.pdf"),p,width = 12,height = 7.2,dpi = 300)


























