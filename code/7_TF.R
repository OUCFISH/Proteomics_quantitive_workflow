rm(list = ls())
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(ggplot2)))

output <- "result/7_TF"
if (!dir.exists(output)) {
  dir.create(output)
}

conf <- "config/default.conf"
source(conf)
contrast.name <- read.xlsx(contrast.file)

total.protein <- read.xlsx("result/2_Data_Summary/Protein_SummaryStat.xlsx") %>%
  dplyr::select(Protein.Accessions, Gene.Name, TF) %>%
  na.omit()
total.tf <- as.data.frame(table(total.protein$TF))
colnames(total.tf) <- c("TF-Family","Count")
total.tf <- total.tf[order(total.tf$Count,decreasing = T),]
total.tf$`TF-Family` <- factor(total.tf$`TF-Family`,levels = total.tf$`TF-Family`)
colnames(total.protein)[grep(pattern = "TF",x = colnames(total.protein))] <- "TF-Family"
write.xlsx(x = total.protein,paste0(output,"/tf.family.xlsx"))

total.tf.barplot <- ggplot(total.tf, aes(x = `TF-Family`, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = Count), vjust = -0.3) + 
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  labs(x = "TF-Family", y = "Number of Proteins") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 15))

ggsave(paste0(output, "/total_tf-barplot.png"),total.tf.barplot,width = 12,height = 7.2,dpi = 300)
ggsave(paste0(output, "/total_tf-barplot.pdf"),total.tf.barplot,width = 12,height = 7.2,dpi = 300)

for (contrast in contrast.name[[1]]) {
  print(contrast)
  suboutdir <- paste0(output,"/",contrast)
  if (!dir.exists(suboutdir)) {
    dir.create(suboutdir)
  }
  pattern <- paste0("^", contrast, "-DEP_results.xlsx")
  diffExp_path <- dir(indir, pattern = pattern) %>% paste0(indir, .)
  diff.prot <- read.xlsx(diffExp_path)
  colnames(diff.prot)[grep(pattern = "Regulation",x = colnames(diff.prot))] <- "Regulation"
  diff.prot <- diff.prot %>% 
    dplyr::select(Protein.Accessions, Regulation, TF) %>% na.omit()
  diff.prot2 <- diff.prot
  colnames(diff.prot2)[grep(pattern = "TF",x = colnames(diff.prot2))] <- "TF-Family"
  write.xlsx(diff.prot2,file = paste0(suboutdir,"/tf.family.xlsx"))
  if (nrow(diff.prot) != 0) {
    df_summary <- diff.prot %>%
      group_by(TF, Regulation) %>%
      summarise(count = n()) %>%
      ungroup()
    df_summary$Regulation <- factor(df_summary$Regulation,levels = c("Up","Down"))
    diff.barplot <- ggplot(df_summary, aes(x = TF, y = count, fill = Regulation)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8),width = 0.8) +
      geom_text(aes(label = count), position = position_dodge(width = 0.8), hjust = -0.8,size = 3, color = "black") +
      coord_flip() +  # 旋转坐标轴
      labs(x = "TF", y = "Number of proteins") +
      scale_fill_manual(values = c("Up" = "#f5780c", "Down" = "#0586c0")) +
      theme_classic() +
      theme(axis.title = element_text(color = "black", size = 20),
            axis.text.y = element_text(color = "black", size = 10,face = "bold"),
            axis.text.x = element_text(color = "black", size = 15),
            legend.title = element_blank()
      )+
      scale_y_continuous(limits = c(0,max(df_summary$count)+0.5))
    
    ggsave(paste0(suboutdir, "/tf-analysis.png"),diff.barplot,width = 12,height = 7.2,dpi = 300)
    ggsave(paste0(suboutdir, "/tf-analysis.pdf"),diff.barplot,width = 12,height = 7.2,dpi = 300)
  }else{
    writeLines("no_result",paste0(suboutdir, "no_result.txt"))
  }
}





















