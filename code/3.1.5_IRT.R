#Author: zhaowu zhang
#Date: 2024-07-23 22:37:05
#FilePath: /PERSONALBIO1/prot2/test/proteomics_workflow/SP231015994/protein_workflow_v1.0.0/code/3.1.5_IRT.R
#Description: 
#LastEditTime: 2024-07-23 22:37:49
#Copyright (c) 2024 by <zhaowu.zhang@personalbio.cn>, All Rights Reserved. 

rm(list = ls())

library(ggplot2)
library(openxlsx)
library(dplyr)
library(ggsci)


# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
IRT_file_path <- "./data/IRT.xlsx"
data_report_path <- "./data/report.tsv"
output <- "./result/3_Quality_Control/3.1_Quality_Control/3.1.5_IRT/"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# 读取数据
# IRT_file_path <- "F:/work_file/test/蛋白组流程/202406011/3.1.1IRT/IRT.xlsx"
# data_report_path <- "F:/work_file/test/蛋白组流程/202406011/3.1.1IRT/report.tsv"
IRT_file <- read.xlsx(IRT_file_path)
data_report <- read.csv(data_report_path, sep = "\t")


# 根据肽段帅选表格
process_IRT <- function(peptide,report_file){
  #description: 根据标准肽段筛选IRT数据表格,并绘制IRT时间图
  #param {标准肽段，IRT数据表格}
  #return {可以绘制IRT时间图的plot数据}

  # 提取表格需要的列
  data_IRT <- report_file[,c("Stripped.Sequence", "Run", "RT")]
  colnames(data_IRT) = c("Peptide", "Sample", "RT")
  
  IRT_plot <- data_IRT %>% 
    filter(Peptide %in% peptide)
  
  return(IRT_plot)
}

# 数据处理
IRT_peptide <- unique(IRT_file$standard)
IRT_data <- process_IRT(IRT_peptide,data_report)


# 加载颜色
IRT_colors <- "Summer"
set_theme_module <-  "./config/set_Theme.R"
source(set_theme_module)

if (IRT_colors %in% names(color.ls)) {
  IRT_colors <- color.ls[[IRT_colors]]
} else if (IRT_colors %in% names(color.continuous.ls)) {
  IRT_colors <- color.continuous.ls[[IRT_colors]]
}


color_num <- length(IRT_peptide)
colors <- IRT_colors[(1:color_num)]

# 绘制iRT洗脱时间图
IRT_plot <- ggplot(IRT_data, aes(x = Sample, y = RT, color = Peptide, group = Peptide, fill = Peptide)) +
  geom_point(size = 3, shape = 23) +  # 绘制圆形数据点
  geom_line(linetype = 3, size = 1) +  # 连接不同肽段的点
  labs(title = "iRT Normalized Retention Time Plot",
       x = "Run",
       y = "Normalization RT") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) + # 设置颜色
  theme_minimal() +  # 使用简洁主题
  theme(legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5), 
        axis.line.y = element_line(color = "black", size = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
        
        )

# 样本数目超过一定值
# if (color_num > 20) {
#   IRT_plot <- IRT_plot + scale_x_continuous(breaks = seq(1, color_num, by = 10)) 
# } else {
#   IRT_plot <- IRT_plot + scale_x_discrete() 
# }

ggsave(paste0(output, "./IRT_plot.pdf"), IRT_plot, width = 12, height = 10, dpi = 300)
ggsave(paste0(output, "./IRT_plot.png"), IRT_plot, width = 12, height = 10, dpi = 300)





