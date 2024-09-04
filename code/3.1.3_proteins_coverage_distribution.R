
#Author: zhaowu zhang
#Date: 2024-07-19 05:28:55
#FilePath: /PERSONALBIO1/prot2/test/proteomics_workflow/SP231015994/protein_workflow_v1.0.0/code/3.1.3_proteins_coverage_distribution.R
#Description: 
#LastEditTime: 2024-07-24 02:48:40
#Copyright (c) 2024 by <zhaowu.zhang@personalbio.cn>, All Rights Reserved. 

# 清除工作区
rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(ggsci)))

# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/2_Data_Summary/Protein_SummaryStat.xlsx"
output <- "./result/3_Quality_Control/3.1_Quality_Control/3.1.3_proteins_coverage_distribution/"
if (!dir.exists(output)) {
  dir.create(output, recursive = T)
}

input.data <- read.xlsx(file_path)
data <- input.data$Coverage
# 将数据按照给定的区间进行分组
groups <- cut(data, breaks = c(0, 10, 20, 30, 40, 50, Inf), labels = c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%", ">50%"))
# 统计每个区间的数据个数
counts <- table(groups)
# 创建数据框
df <- data.frame(Coverage = names(counts), Count = as.numeric(counts))
  # 绘制饼图并添加数字
p <- ggplot(df, aes(x = "", y = Count, fill = Coverage)) +
  geom_bar(stat = 'identity',width = 1, color = "white") +
  geom_text(aes(x = 1.65, label = Count), size = 8, position = position_stack(vjust = 0.5),) +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "white"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        legend.key.size = unit(1, "cm"),
        legend.key.spacing = unit(1, "cm"),
        legend.key.spacing.y = unit(0.2, "cm"),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        ) +
  
  scale_fill_brewer(palette = "Set1", breaks = c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%", ">50%")) 

# 添加引导线 
plot_info = ggplot_build(p)
y_coord = plot_info$data[[2]]$y
p <- p + geom_segment(aes(x = 1.5, xend = 1.55, y = y_coord, yend = y_coord), color = "black", linewidth = 0.5) 



# total_count <- sum(df$Count)
# df <- df %>%
#   mutate(Percentage = Count / total_count * 100)  
# 
# coverage_factor <- factor(df$Coverage, levels = c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%", ">50%"))
# proportions <- table(coverage_factor) / length(coverage_factor)
# 
# pie(df$Percentage, labels = df$Count, col = brewer.pal(length(proportions), "Set1"))


ggsave(paste0(output, "proteins_coverage_distribution.png"),p,width = 12,height = 10,dpi = 300)
ggsave(paste0(output, "proteins_coverage_distribution.pdf"),p,width = 12,height = 10,dpi = 300)




















