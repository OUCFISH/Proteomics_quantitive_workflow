rm(list = ls())
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(matrixStats)))
protein.data <- read.xlsx("result/2_Data_Summary/Protein_SummaryStat.xlsx")
group.file <- read.xlsx("./temp/group.xlsx")
protein.intensity <- protein.data[,c(group.file$SampleName)]
group.name <- unique(group.file$Group)
for (i in 1:length(group.name)) {
  rsd.data <- group.name[i]
  group_index <- which(colnames(protein.intensity) %in% group.file$SampleName[group.file$Group == group.name[i]])
  protein.intensity <- protein.intensity %>%
    dplyr::mutate(!!rsd.data := rowSds(as.matrix(dplyr::select(., colnames(protein.intensity)[group_index])), na.rm = TRUE) /
             rowMeans(dplyr::select(., colnames(protein.intensity)[group_index]), na.rm = TRUE))
}
rsd.plot <- protein.intensity[,c(group.name)]

protein_long <- rsd.plot %>%
  pivot_longer(cols = everything(), names_to = "Group", values_to = "RSD")
conf <- "./config/default.conf"
source(conf)
source("./config/set_Theme.R")
if (group_colors %in% names(color.ls)) {
  palette <- color.ls[[group_colors]]
}
colors <- setNames(palette[1:length(unique(group.file$Group))],unique(group.file$Group))

# 绘制小提琴图
violinplot<- ggplot(protein_long, aes(x = Group, y = RSD, fill = Group)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  labs(x = "group",
       y = "Relative Standard Deviation (RSD)") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 18),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        ) + 
  scale_fill_manual(values = colors)

output <- "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.4_RSD/"
if (!dir.exists(output)) {
  dir.create(output,recursive = T)
}
ggsave(paste0(output, "RSD.png"),violinplot,width = 12,height = 7.2,dpi = 300)
ggsave(paste0(output, "RSD.pdf"),violinplot,width = 12,height = 7.2,dpi = 300)





























