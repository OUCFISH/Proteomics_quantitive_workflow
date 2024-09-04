# 清除工作区
rm(list = ls())
# 加载库
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(tidyverse)))
# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/2_Data_Summary/Peptide.xlsx"
output <- "./result/3_Quality_Control/3.1_Quality_Control/3.1.1_peptides_length_distribution/"
if (!dir.exists(output)) {
  dir.create(output,recursive = T)
}

# 导入数据
input.data <- read.xlsx(file_path)
if ("Stripped.Sequence"%in%colnames(input.data)) {
  input.data$peptides_length <- nchar(input.data$Stripped.Sequence)
  temp.data <- data.frame(sequence = input.data$Stripped.Sequence,length = input.data$peptides_length)
}else{
  input.data[,1]<-sapply(input.data[,1],function(x) gsub("\\[","",x=x)%>%gsub("\\]","",x=.)%>%gsub("\\.","",x=.))
  input.data$peptides_length <- nchar(input.data$Annotated.Sequence)
  temp.data <- data.frame(sequence = input.data$Annotated.Sequence,length = input.data$peptides_length)
}
# temp.data$length_group <- ifelse(temp.data$length >= 30, ">30", as.character(temp.data$length))
frequency <- table(temp.data$length)
# 将结果转换为数据框
# df <- data.frame(number = names(frequency), frequency = as.numeric(frequency))
df <- data.frame(number = as.numeric(names(frequency)), frequency = as.numeric(frequency))
p <- ggplot(df, aes(x = factor(number),
                               # levels = c(as.character(1:30), ">30")),
                    y = frequency)) +
  geom_bar(stat = "identity", fill = "#0586c0", width = 0.9) +
  # scale_x_discrete(expand = expansion(add = c(2, 2))) +
  labs(x = "Peptides Length", y = "Number of Peptides")+
  theme_bw() +
  theme(
    axis.title = element_text(color = "black", size = 36),
    axis.text = element_text(color = "black", size = 24),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    #panel.border = element_rect(linewidth = 1.5),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1,size = 22),
    ) + 
  geom_text(aes(label = frequency), vjust = -0.7, size = 5.7) + 
  scale_y_continuous(expand = c(0,0),limits = c(0, max(df$frequency) * 1.1))



ggsave(paste0(output, "peptides_distribution.png"), p, width  = 500,height = 400,units = 'mm', dpi = 300)
ggsave(paste0(output, "peptides_distribution.pdf"), p, width  = 500,height = 400,units = 'mm', dpi = 300)




















