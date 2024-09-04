# 清除工作区
rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(dplyr)))

# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/2_Data_Summary/Protein_SummaryStat.xlsx"
output <- "./result/3_Quality_Control/3.1_Quality_Control/3.1.4_proteins_isoelectric_distribution/"
if (!dir.exists(output)) {
  dir.create(output, recursive = T)
}
# 导入数据
input.data <- read.xlsx(file_path)
plot.df <- data.frame(x = input.data$Isoelectric.Point)
plot.df <- na.omit(plot.df)
# 统计结果以及获取最大值
# count.data <- ggplot(plot.df, aes(x = x)) +
#   geom_histogram(binwidth = 0.5, fill = "#0586c0", color = "black")
# max.count <- max(count.data$data$count)

bins <- cut(plot.df$x, breaks = seq(round(min(plot.df$x),0), round(max(plot.df$x),0),by = 0.5))
hist.df <- as.data.frame(table(bins))
max.count <- max(hist.df$Freq)


p <- ggplot(hist.df, aes(x = bins, y = Freq)) +
  geom_col(fill = logo_blue, color = "white") +
  labs(title = NULL, x = "Isoelectric Point", y = "Number of Proteins")+
  geom_text(aes(label = Freq), vjust = -0.5, size = 5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(color = "black", size = 20),
    axis.text.x = element_text(angle = 45, color = "black", size = 14, hjust = 1),
    axis.text.y = element_text(color = "black", size = 14),
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(0.3, "cm"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
  ) + 
  scale_y_continuous(expand = c(0,0),limits = c(0,max.count * 1.1)) 



ggsave(paste0(output, "proteins_isoelectric_distribution.png"),p,width = 14,height = 8,dpi = 300)
ggsave(paste0(output, "proteins_isoelectric_distribution.pdf"),p,width = 14,height = 8,dpi = 300)



























