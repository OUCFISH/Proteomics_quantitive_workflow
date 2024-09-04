rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(openxlsx)))

# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/2_Data_Summary/Protein_SummaryStat.xlsx"
contrast_file <- "./temp/contrast.xlsx"
output <- "./result/4_Diff_Expressed/4.3_Volcano/"

# 读取数据
data <- read.xlsx(file_path)
contrast_names <- read.xlsx(contrast_file)

# 循环绘制不同对比策略下的volcano plot
for (contrast_item in contrast_names$contrast) {

# print(contrast_item)
if (!dir.exists(paste0(output, contrast_item))) {
  dir.create(paste0(output, contrast_item), recursive = TRUE)
}
# contrast_item <- "M_vs_K" 
# 数据预处理：
data$p_value <- data[[paste0(contrast_item,".pvalue")]]
data$log10Pvalue <- -log10(data[[paste0(contrast_item,".pvalue")]])
data$log2FoldChange <- data[[paste0(contrast_item,".Log2fc")]]

Log10Pvalue <- -log10(pvalue)
Log2FoldChange <- log2(Fold_Change1)


# 添加一个列作为 size_category，根据条件设定为 "Small" 或 "Large"
data$size_category <- ifelse(data$log10Pvalue > Log10Pvalue & abs(data$log2FoldChange) > Log2FoldChange, "Large", "Small")

 # 区分上调和下调基因，设定颜色
data$color <- ifelse(data$log2FoldChange > Log2FoldChange & data$log10Pvalue > Log10Pvalue, "Up regulated", 
                   ifelse(data$log2FoldChange < -Log2FoldChange & data$log10Pvalue > Log10Pvalue, "Down regulated", "Not Significant"))

# 计算每个基因的绝对 Fold Change
data$abs_logFC <- abs(data$log2FoldChange)

# 根据差异标准筛选数据集
df_filtered <- data[data$log10Pvalue > Log10Pvalue & abs(data$log2FoldChange) > Log2FoldChange, ]

# 根据绝对 Fold Change 排序数据集并选择前10个基因
top_Protein_up <- head(df_filtered[df_filtered$log2FoldChange > Log2FoldChange, ][order(df_filtered[df_filtered$log2FoldChange > Log2FoldChange, ]$abs_logFC, decreasing = TRUE), ]$Protein.Accessions, 10)
top_Protein_down <- head(df_filtered[df_filtered$log2FoldChange < -Log2FoldChange, ][order(df_filtered[df_filtered$log2FoldChange < -Log2FoldChange, ]$abs_logFC, decreasing = TRUE), ]$Protein.Accessions, 10)

# 创建带数量信息的图例标签
data$color <- factor(data$color, levels = c("Up regulated", "Down regulated", "Not Significant")) 
legend_labels <- c(paste0("Up regulated (", sum(data$color == "Up regulated"), ")"),
                   paste0("Down regulated (", sum(data$color == "Down regulated"), ")"),
                   paste0("Not Significant (", sum(data$color == "Not Significant"), ")"))
                   
                     
#legend_labels <- factor(legend_labels, levels = c(paste0("Upregulated (", sum(data$color == "Upregulated"), ")"),
#                                                  paste0("Downregulated (", sum(data$color == "Downregulated"), ")"),
#                                                  paste0("Not Significant (", sum(data$color == "Not Significant"), ")")))              

p_value <- pvalue
legend_value <- c((paste0(
                "Fold change 1", " > ", (Fold_Change1), "\n",
                "Fold change 2", " < ", round((Fold_Change2), 2))),
                  paste0("P value < ", p_value))

#print(legend_labels)

# 绘制火山图
# topk
volcano_plot_topk <- ggplot(data, aes(x = log2FoldChange, y = log10Pvalue, color = color, size = size_category)) +
  geom_point(alpha = 0.8, shape = 21, stroke = 0.5, aes(fill = color)) + # 使用映射的颜色作为点的填充色和边框色
  scale_color_manual(name = "Regulation",
                     values = c("Up regulated" = logo_orange, "Down regulated" = logo_blue , "Not Significant" = "#b1b2b4"),
                     labels = legend_labels,
                     #guide = "none"
                    ) + # 设置颜色和图例标签
  scale_fill_manual(name = "Regulation", 
                    values = c("Up regulated" = logo_orange , "Down regulated" = logo_blue, "Not Significant" = "#b1b2b4"),
                    labels = legend_labels,
                     #guide = "none"
                    ) + 
  scale_size_manual(values = c("Small" = 2, "Large" = 2), 
                    name = "Threshold",
                    labels = legend_value,
                    # guide = "none",
                    ) + # 控制点大小，不显示图例
  geom_label_repel(data = subset(data, Protein.Accessions %in% top_Protein_up), # 添加上调前10个基因的标签并防止重叠
                   aes(label = as.character(Protein.Accessions)),
                   size = 3, # 文字大小
                   color = logo_orange,
                   # family = "TT Arial",
                   segment.color = "grey50", # 显示连接线
                   segment.size = 0.2, # 连接线大小
                   box.padding = unit(0.5, "lines"), # 文字标签与边框之间的距离
                   point.padding = unit(0.3, "lines"), # 文字标签与点之间的距离
                   max.overlaps = 20, # 最大重叠数
                   force = 1, # 排斥力，控制标签的分布
                   fill = "white", # 标签背景颜色
                   label.size = 0.5, # 标签边框粗细
                   label.r = unit(0.15, "lines"), # 标签边框圆角半径
                   label.padding = unit(0.25, "lines"), # 标签边框内间距
                   min.segment.length = unit(0.2, "lines"),# 最小连接线长度
                   direction = "both",
                   ) + 
  
  geom_label_repel(data = subset(data, Protein.Accessions %in% top_Protein_down), # 添加下调前五个基因的标签并防止重叠
                   aes(label = as.character(Protein.Accessions)),
                   size = 3, # 文字大小
                   color = logo_blue,
                   # family = "TT Arial",
                   segment.color = "grey50", # 显示连接线
                   segment.size = 0.2, # 连接线大小
                   box.padding = unit(0.5, "lines"), # 文字标签与边框之间的距离
                   point.padding = unit(0.3, "lines"), # 文字标签与点之间的距离
                   max.overlaps = 20, # 最大重叠数
                   force = 1, # 排斥力，控制标签的分布
                   fill = "white", # 标签背景颜色
                   label.size = 0.5, # 标签边框粗细
                   label.r = unit(0.15, "lines"), # 标签边框圆角半径
                   label.padding = unit(0.25, "lines"), # 标签边框内间距
                   min.segment.length = unit(0.2, "lines"), # 最小连接线长度
                   direction = "both",
                   ) +
  guides(color = guide_legend(override.aes = list(size = 2), )) +
  geom_hline(yintercept = Log10Pvalue, linetype = "dotted", color = "darkgrey", size = 0.75) + # 添加p值阈值线
  geom_vline(xintercept = c(-Log2FoldChange, Log2FoldChange), linetype = "dotted", color = "darkgrey", size = 0.75) + # 添加fold change阈值线
  labs(title = "Volcano Plot", 
       x = expression(Log[2] ~ "(Fold change)"), 
       y = expression(-Log[10] ~ "(P value)")) + 
  theme_minimal() + # 使用简约主题
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 14), 
    legend.position = "right", 
    legend.text = element_text(size = 16), 
    legend.title = element_text(size = 16, face = "bold"), 
    legend.key.size = unit(1, "cm"),
    legend.spacing.y = unit(1, "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "grey", size = 1, fill = NA),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.line.x = element_line(color = "grey", linewidth = 0.5),
    axis.line.y = element_line(color = "grey",, linewidth = 0.5),
  )

 #print(volcano_plot_topk)
# 保存图片
ggsave(paste0(output,contrast_item,"/", contrast_item, "_volcano_topk.png"), plot = volcano_plot_topk, width = 10, height = 8, dpi = 300)
ggsave(paste0(output,contrast_item,"/", contrast_item, "_volcano_topk.pdf"), plot = volcano_plot_topk, width = 10, height = 8, dpi = 300)


# no_topk
volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = log10Pvalue, color = color, size = size_category)) +
  geom_point(alpha = 0.8, shape = 21, stroke = 0.5, aes(fill = color)) + # 使用映射的颜色作为点的填充色和边框色
  scale_color_manual(name = "Regulation", 
                     values = c("Up regulated" = logo_orange, "Down regulated" =  logo_blue, "Not Significant" = "#b1b2b4"),
                     labels = legend_labels,
                     #guide = "none"
                    ) + # 设置颜色和图例标签
  scale_fill_manual(name = "Regulation", 
                    values = c("Up regulated" = logo_orange, "Down regulated" = logo_blue, "Not Significant" = "#b1b2b4"),
                    labels = legend_labels,
                    #guide = "none"
                    ) + 
  scale_size_manual(values = c("Small" = 2, "Large" = 2), 
                    name = "Threshold",
                    labels = legend_value,) + 
  guides(color = guide_legend(override.aes = list(size = 2))) +
  geom_hline(yintercept = Log10Pvalue, linetype = "dotted", color = "darkgrey", size = 0.75) + # 添加p值阈值线
  geom_vline(xintercept = c(-Log2FoldChange, Log2FoldChange), linetype = "dotted", color = "darkgrey", size = 0.75) + # 添加fold change阈值线
  labs(title = "Volcano Plot", 
       x = expression(Log[2] ~ "(Fold change)"), 
       y = expression(-Log[10] ~ "(P value)")) +
  theme_minimal() + # 使用简约主题
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 14), 
    legend.position = "right", 
    legend.text = element_text(size = 16), 
    legend.title = element_text(size = 16, face = "bold"), 
    legend.key.size = unit(1, "cm"),
    legend.spacing.y = unit(1, "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "grey", size = 1, fill = NA),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.line.x = element_line(color = "grey", linewidth = 0.5),
    axis.line.y = element_line(color = "grey", linewidth = 0.5),
  )


# print(volcano_plot)
ggsave(paste0(output,contrast_item,"/", contrast_item, "_volcano.png"), plot = volcano_plot, width = 10, height = 8, dpi = 300)
ggsave(paste0(output,contrast_item,"/", contrast_item, "_volcano.pdf"), plot = volcano_plot, width = 10, height = 8, dpi = 300)

}
