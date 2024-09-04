# 清除工作区
rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggsci)))
suppressWarnings(suppressMessages(library(tidyverse)))
# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/2_Data_Summary/Protein_SummaryStat.xlsx"
output <- "./result/3_Quality_Control/3.1_Quality_Control/3.1.2_protein_molecular_distribution/"
if (!dir.exists(output)) {
  dir.create(output, recursive = T)
}



input.data <- read.xlsx(file_path)

# 确定区间宽度去除NA值
data <- input.data$`Molecular.Weight[KDa]`
data <- na.omit(data)
bin_width <- 20# 计算每个数据的区间
bins <- cut(data, breaks = seq(0, max(data) + bin_width, by = bin_width), right = FALSE)
frequency <- table(bins)

# 创建数据框
df <- data.frame(Interval = as.character(names(frequency)), frequency = as.numeric(frequency))

rows_number <- grep("^\\[200", df$Interval)
sum_freq_gt_200 <- sum(df$frequency[(rows_number):nrow(df)])
df <- rbind(df[1:(rows_number-1), ], 
              data.frame(Interval = ">200", frequency = sum_freq_gt_200))
df$Interval <- factor(x = df$Interval, levels = df$Interval)

# 创建颜色函数
color_func <- colorRampPalette(c(logo_blue, "lightblue"))
n <- length(unique(df$Interval))
colors <- color_func(n)

# 绘制柱状图
p <- ggplot(df, aes(x = Interval, y = frequency, fill = Interval)) +
  geom_bar(stat = "identity") +
  labs(x = "Molecular Weight(KDa)", y = "Number of Proteins") +
  scale_fill_manual(values = colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 24)) +
  geom_text(aes(label = frequency), vjust = -0.5, size = 10) +
  theme(
    axis.title = element_text(color = "black", size = 36),
    axis.text = element_text(color = "black", size = 28),
    legend.position = "none",
    # legend.text = element_text(size = 28),
    # legend.title = element_text(size = 40),
    # legend.key.size = unit(1, "cm"),
    # legend.key.spacing = unit(1, "cm"),
    # legend.key.spacing.y = unit(0.2, "cm"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.ticks = element_line(linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.3, "cm")
  ) +
  scale_y_continuous(expand = c(0,0),limits = c(0, max(df$frequency) * 1.1)) 



ggsave(paste0(output, "proteins_molecular_distribution.png"), p, width  = 600,height = 400,units = 'mm', dpi = 300)
ggsave(paste0(output, "proteins_molecular_distribution.pdf"), p, width  = 600,height = 400,units = 'mm', dpi = 300)























