rm(list = ls())
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(plotly)))
suppressWarnings(suppressMessages(library(htmlwidgets)))
suppressWarnings(suppressMessages(library(ggplot2)))
output <- "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.8_PCA/"
if (!dir.exists(output)) {
  dir.create(output,recursive = T)
}
protein.data <- read.xlsx("temp/pca_data.xlsx")
pca.data <- as.data.frame(t(protein.data))
colnames(pca.data) <- pca.data[1,]
pca.data <- pca.data[-1,]
pca.data[] <- lapply(pca.data, as.numeric)

pca.data$sample <- rownames(pca.data)
group.data <- read.xlsx("temp/group.xlsx")
# 检查 pca.data$sample 中是否包含 "QC"
if (any(grepl("QC", pca.data$sample))) {
  pca.data$group <- group.data$Group[match(pca.data$sample, group.data$SampleName)]
  pca.data$group[grep(pattern = "QC", pca.data$sample)] <- "QC"
} else {
  pca.data$group <- group.data$Group[match(pca.data$sample, group.data$SampleName)]
}
pca.data <- pca.data[,c((ncol(pca.data)-1):ncol(pca.data), 1:(ncol(pca.data)-2))]
data <- pca.data
# 分离变量和类标签
pca_data <- data[, 3:ncol(data)]  # 假设数据从第三列开始是变量
group <- data[["group"]]
sample_id <- data[["sample"]]

# PCA
pca <- prcomp(pca_data, center = TRUE, scale. = TRUE)
pca_data_transformed <- data.frame(pca$x)
pca_data_transformed$group <- group
pca_data_transformed$ID <- sample_id
pca_variance <- pca$sdev^2 / sum(pca$sdev^2)
pca.data <- pca_data_transformed
pca.data <- pca.data[,c(ncol(pca.data),ncol(pca.data)-1,1:3)]
write.xlsx(x = pca.data,file = paste0(output,"PCA.xlsx"))
# 动态生成颜色调色板
unique_classes <- unique(group)
num_classes <- length(unique_classes)
conf <- "./config/default.conf"
source(conf)
source("./config/set_Theme.R")
if (group_colors %in% names(color.ls)) {
  palette <- color.ls[[group_colors]]
}
colors <- setNames(palette[1:num_classes],unique_classes)
# 设置是否显示标签的逻辑变量
show_labels <- TRUE  # 这里可以根据需要设置为 TRUE 或 FALSE


# 绘制2D PCA图，inherit.aes = TRUE或F控制是否显示分组置信区间置信区间
p2d <- ggplot(pca_data_transformed, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5, color = "grey") + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.5, color = "grey") +
  # stat_ellipse(aes(x = PC1, y = PC2), linetype = 2, size = 0.5, level = 0.95, inherit.aes = TRUE) +
  theme_bw() +
  scale_color_manual(values = colors) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
        
  labs(x = paste0("PC1: ", sprintf("%.1f%%", signif(pca_variance[1] * 100, 3))), 
       y = paste0("PC2: ", sprintf("%.1f%%", signif(pca_variance[2] * 100, 3))), 
       title = paste0("PCA plot")) +
  theme(plot.title = element_text(hjust = 0.5))

# 添加标签显示的条件
if (show_labels) {
  p2d2 <- p2d + geom_text(aes(label = ID), size = 3, hjust = 0.5, vjust = 2, show.legend  = FALSE)
}

# 保存2D PCA图为PNG格式

ggsave(paste0(output, "PCA_2D.png"), plot = p2d, width = 8, height = 6, dpi = 300)
ggsave(paste0(output, "PCA_2D.pdf"), plot = p2d, width = 8, height = 6, dpi = 300)
ggsave(paste0(output, "PCA_2D_label.png"), plot = p2d2, width = 8, height = 6, dpi = 300)
ggsave(paste0(output, "PCA_2D_label.pdf"), plot = p2d2, width = 8, height = 6, dpi = 300)
# 绘制3D PCA图
fig3d <- plotly::plot_ly(pca_data_transformed, 
                 x = ~PC1, y = ~PC2, z = ~PC3, 
                 color = ~group, 
                 colors = colors, 
                 type = "scatter3d",
                 mode = "markers",
                 text = ~text,
                 marker = list(size = 5)) %>%
  layout(scene = list(
    xaxis = list(title = paste0("PC1: ", signif(pca_variance[1] * 100, 3), "%")),
    yaxis = list(title = paste0("PC2: ", signif(pca_variance[2] * 100, 3), "%")),
    zaxis = list(title = paste0("PC3: ", signif(pca_variance[3] * 100, 3), "%"))
  ),
  
  title = "3D PCA Plot")

# 保存3D PCA图为HTML格式
htmlwidgets::saveWidget(fig3d, paste0(output, "3D_PCA_plot.html"))

if (file.exists(paste0(output,"3D_PCA_plot_files"))) {
  # 删除目标文件夹及其内容
  unlink(paste0(output,"3D_PCA_plot_files"), recursive = TRUE)
}























