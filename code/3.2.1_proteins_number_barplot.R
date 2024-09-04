# 清除工作区
rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(impute)))

# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/2_Data_Summary/ProteinGroup.xlsx"
group_path <- "./temp/group.xlsx"
output <- "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.1_Expression_Statistics/"
if (!dir.exists(output)) {
  dir.create(output, recursive = T)
}

# 导入数据
proteins_data <- read.xlsx(file_path)
group_file <- read.xlsx(group_path)
sample_names <- group_file$SampleName
col_numbers <- sapply(sample_names, function(x) which(colnames(proteins_data) == x))
non_empty_counts <- apply(proteins_data[,col_numbers], 2, function(col) sum(!is.na(col)))

# 将统计结果转换为数据框
counts_df <- data.frame(Column = sample_names, Count = non_empty_counts)

# 绘制柱形图
p <- ggplot(counts_df, aes(x = Column, y = Count)) +
  geom_bar(stat = "identity", fill = logo_blue, width = 0.5) +
  labs(x = NULL, y = "Number of Proteins") +
  theme_classic() +
  theme(axis.title = element_text(color = "black", size = 20),
        axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.ticks = element_line(linewidth = 0.5),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        ) +
  geom_text(aes(label = Count), vjust = -0.6, size = 4) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(counts_df$Count) * 1.1))  # 设置y轴范围，使其起始点为0


# 保存图片
ggsave(paste0(output, "proteins_number_barplot.png"),p,width = 12,height = 8,dpi = 300)
ggsave(paste0(output, "proteins_number_barplot.pdf"),p,width = 12,height = 8,dpi = 300)






















