# 清除工作区
rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(writexl)))

# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/2_Data_Summary/Protein_SummaryStat.xlsx"
group_file <- "./temp/group.xlsx"
output <- "./result/4_Diff_Expressed/4.1_DiffStats/"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}


# 数据处理
data <- read.xlsx(file_path)
regulation_cols <- names(data)[str_detect(names(data), "Regulation$")]
selected_data <- data[, regulation_cols]
selected_data <- as.data.frame(selected_data)
new_col_names <- sapply(strsplit(regulation_cols, "\\."), `[`, 1)
colnames(selected_data) <- new_col_names

long_data <- pivot_longer(selected_data, cols = everything(), names_to = "Variable", values_to = "Category")

# 对每个变量的每个类别进行计数
category_counts <- long_data %>% dplyr::count(Variable, Category)
data_filtered <- category_counts %>% filter(Category != "Nodiff")

# 绘图参数
fill_label <- factor(data_filtered$Category, levels = c("Up", "Down"))
color_mapping <- c("Up" = logo_orange, "Down" = logo_blue)
x_lable <- data_filtered$Variable
label <- c("Up", "Down" )
 
# 绘制柱状图
diffexp_plot <- ggplot(data_filtered, aes(x = x_lable, y = n, fill = fill_label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5),width = 0.5) +
  scale_fill_manual(values = color_mapping, name = NULL, label = label) +
  labs(title = NULL, x = NULL, y = "Number of Proteins") +
  theme_minimal() +
  theme(legend.position = "right") +
  # theme(axis.text.x = element_blank()) 
  theme(axis.line = element_line(color = "black", linewidth = 1),
      axis.text = element_text( color = "black", size = 18),
      axis.title = element_text(color = "black", size = 24),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.3, "cm"),
      legend.text = element_text(color = "black", size = 18),
      legend.key.size = unit(1, "cm"),
      plot.margin = margin(1, 1, 1, 1, "cm"),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      ) +
  geom_text(aes(label = round(n, 0)),              
          position = position_dodge(width = 0.5),     
          vjust = -0.5,
          size = 6, color = "black",
          parse = FALSE,
          check_overlap = FALSE, na.rm = FALSE) 
  
diffexp_plot <- diffexp_plot + scale_y_continuous(expand = c(0,0),limits = c(0, max(data_filtered$n) * 1.1))

# print(diffexp_plot)

# 保存图片
ggsave(paste0(output, "DEP_stat.png"), plot = diffexp_plot, width = 12, height = 10, dpi = 300)
ggsave(paste0(output, "DEP_stat.pdf"), plot = diffexp_plot, width = 12, height = 10, dpi = 300)


# 输出表格
down_table <- data_filtered %>% filter(Category != "Up")
up_table <- data_filtered %>% filter(Category != "Down")

DEP_stat <-  data.frame(matrix(ncol = 4, nrow = nrow(down_table)))
col_names <- c("Control_vs_Treat","Up Regulation","Down Regulation","Total")
colnames(DEP_stat) = col_names


for (i in (1:nrow(down_table))) {
  DEP_stat[i, "Control_vs_Treat"] <- up_table[i,1]
  DEP_stat[i, "Up Regulation"] <- up_table[i,3]
  DEP_stat[i, "Down Regulation"] <- down_table[i, 3]
  DEP_stat[i, "Total"] <- DEP_stat[i, "Up Regulation"] + DEP_stat[i, "Down Regulation"]
}

write_xlsx(DEP_stat, path = paste0(output,"DEP_Stat_results.xlsx"))


