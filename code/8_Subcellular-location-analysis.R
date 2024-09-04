rm(list = ls())
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(RColorBrewer)
library(plotly)
library(webshot)

output <- "result/8_Subcellular-location"
if (!dir.exists(output)) {
  dir.create(output)
}
conf <- "config/default.conf"
source(conf)
contrast.name <- read.xlsx(contrast.file)
subdb <- read.csv(Sys.glob(paste0(subcellular_file_path, organism_Scientific_name, "/*csv")))
subdb$UniprotID <- sapply(strsplit(subdb$Protein_ID,split = "\\|"),FUN = function(x)x[2])
sub.db <- dplyr::select(subdb,UniprotID,Localizations)
sub.db$Localizations <- sapply(strsplit(sub.db$Localizations,split = "\\|"),FUN = function(x)x[1])
source("./config/set_Theme.R")
group_colors <- "Summer"
if (group_colors %in% names(color.ls)) {
  palette <- color.ls[[group_colors]]
}
for (contrast in contrast.name$contrast) {
  print(contrast)
  suboutdir <- paste0(output,"/",contrast)
  if (!dir.exists(suboutdir)) {
    dir.create(suboutdir)
  }
  pattern <- paste0("^", contrast, "-DEP_results.xlsx")
  diffExp_path <- dir(indir, pattern = pattern) %>% paste0(indir, .)
  diff.prot <- read.xlsx(diffExp_path) %>% dplyr::rename(Protein=Protein.Accessions,Genes=Gene.Name)
  colnames(diff.prot)[grep(pattern = "FC",x = colnames(diff.prot))] <- "FC"
  diff.prot <- diff.prot[,c("Protein","FC")]
  diff.prot.sub <- left_join(x = diff.prot,y = sub.db,by = c("Protein"="UniprotID"))
  diff.prot.sub <- na.omit(diff.prot.sub)

## pie plot
  pie.df <- as.data.frame(table(diff.prot.sub$Localizations))
  names(pie.df) <- c("Location", "Count")
  pie.df <- pie.df %>%
    mutate(Percent = Count / sum(Count) * 100)
  colors <- setNames(palette[1:length(unique(pie.df$Location))],unique(pie.df$Location))
  
  ## 用Ploty绘制pie图
  pieplot <- plot_ly(data = pie.df, 
                     labels = ~Location, 
                     values = ~Percent,
                     type = 'pie',
                     textinfo = 'percent',  # 显示百分比和标签
                     textposition = 'outside',  # 标签位置
                     size = I(22),
                     marker = list(colors = colors)) %>%  # colors已经定义好
  layout(
    title = list(text = "Subcellular Location Analysis", font = list(size = 24)),
    
    showlegend = TRUE,
    legend = list(
      font = list(size = 24),
      x = 0.9,  # 调整图例位置
      y = 0.5
    ),
    margin = list(l = 100, r = 100, b = 100, t = 100)) # 设置四周的留白
  # pieplot
  
   # 保存为HTML文件
  htmlwidgets::saveWidget(pieplot, file = "temp/sub-pieplot.html")
  webshot("temp/sub-pieplot.html", paste0(suboutdir, "/sub-pieplot.pdf"),
          vwidth = 1500, vheight = 900)
  
  ## barplot
  diff.prot.sub$regulation <- ifelse(test = diff.prot.sub$FC > Fold_Change1,"up","down")
  regulation_count <- diff.prot.sub %>%
    group_by(Localizations, regulation) %>%
    summarise(count = n()) %>%
    spread(regulation, count, fill = 0)
  protein_list <- diff.prot.sub %>%
    group_by(Localizations, regulation) %>%
    summarise(Proteins = paste(Protein, collapse = ",")) %>%
    spread(regulation, Proteins, fill = "")
  colnames(protein_list)[2:3] <- c("down_protein","up_protein")
  combined_data <- left_join(regulation_count, protein_list, by = "Localizations")
  write.xlsx(x = combined_data,file = paste0(suboutdir,"/sub.xlsx"))
  sub.long_data <- tidyr::pivot_longer(combined_data, cols = c(up, down), names_to = "Direction", values_to = "Count")
  sub.long_data$Direction <- factor(sub.long_data$Direction, levels = c("up", "down"))
  barplot <- ggplot(sub.long_data, aes(x = Localizations, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),width = 0.8) +
    geom_text(aes(label = Count), position = position_dodge(width = 0.8), hjust = -0.8,size = 3, color = "black") +
    coord_flip() +  # 旋转坐标轴
    labs(x = "Subcellular location", y = "Number of proteins") +
    scale_fill_manual(values = c("up" = "#f5780c", "down" = "#0586c0")) +
    theme_classic() +
    theme(axis.title = element_text(color = "black", size = 20),
          axis.text.y = element_text(color = "black", size = 10,face = "bold"),
          axis.text.x = element_text(color = "black", size = 15),
          legend.title = element_blank()
    )+
    scale_y_continuous(limits = c(0,max(sub.long_data$Count)+0.5))
  ggsave(paste0(suboutdir, '/sub-barplot.png'), barplot, width  = 350,height = 200,units = 'mm', dpi = 300)
  ggsave(paste0(suboutdir, '/sub-barplot.pdf'), barplot, width  = 350,height = 200,units = 'mm', dpi = 1200)
}