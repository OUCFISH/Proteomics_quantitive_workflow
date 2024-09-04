# 清除工作区
rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(ggplot2)))

# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/2_Data_Summary/Protein_SummaryStat.xlsx"
group_file <- "./temp/group.xlsx"
output <- "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.5_Expression_Scatter/"
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
proteins.data <- read.xlsx(file_path)
group_file <- read.xlsx(group_file)
proteins.data <- process_mat(proteins.data,group_file$SampleName)
# proteins.data <- proteins.data[ ,(9:(9+length(group_file$Group)-1))]

log_data <- log10(proteins.data)
sorted_log_data <- apply(log_data, 2, function(x) sort(x, na.last = TRUE,decreasing = T))
# 转换为数据框并添加序号
sorted_log_df <- data.frame(
  index = 1:nrow(sorted_log_data),
  sorted_log_data
)
# 颜色配置
#color_func <- colorRampPalette(ggsci::pal_locuszoom("default")(7))
n <- length(unique(sorted_log_df))
#colors <- color_func(n)

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
colors <- sample_colors[1:n] 


# 转换为长格式

long_data <- reshape2::melt(sorted_log_df, id.vars = "index", variable.name = "variable", value.name = "value")
# 绘制散点图
p <- ggplot(long_data, aes(x = index, y = value, color = variable)) +
  geom_point() +
  labs(x = "rank", y = expression(Log[10]~"(Intensity)"), color = "Sample") +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(color = "black", size = 18),
    axis.text = element_text(color = "black", size = 14),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.3, "cm"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
  ) 

ggsave(paste0(output, "proteins_intensity_scatterplot.png"),p,width = 12,height = 7.2,dpi = 300)
ggsave(paste0(output, "proteins_intensity_scatterplot.pdf"),p,width = 12,height = 7.2,dpi = 300)














