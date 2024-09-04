rm(list = ls())
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(tidyverse)))

output <- "result/10_Kinase"
if (!dir.exists(output)) {
  dir.create(output)
}
conf <- "config/default.conf"
source(conf)
contrast.name <- read.xlsx(contrast.file)
kinase.database <- read.xlsx(Sys.glob(paste0(kinase_file_path,organism_category,"/",organism_code,"_",organism_ID,"_Kinase.xlsx")))
for (contrast in contrast.name[[1]]) {
  print(contrast)
  suboutdir <- paste0(output,"/",contrast)
  if (!dir.exists(suboutdir)) {
    dir.create(suboutdir)
  }
  pattern <- paste0("^", contrast, "-DEP_results.xlsx")
  diffExp_path <- dir(indir, pattern = pattern) %>% paste0(indir, .)
  diff.prot <- read.xlsx(diffExp_path)
  colnames(diff.prot)[c(2,3)] <- c("Protein","Genes")
  colnames(diff.prot)[grep(pattern = "Log2fc",x = colnames(diff.prot))] <- "Logfc"
  gene.kinase <- diff.prot[diff.prot$kinase == "TRUE",c("Protein","Genes","Logfc")]
  gene.kinase$kinase.name <- kinase.database$Protein.names[base::match(gene.kinase$Protein,kinase.database$Entry)]
  if (nrow(gene.kinase) == 0) {
    error_plot <- ggplot() +
      theme_void() +
      geom_text(aes(0, 0, label = "The differential experssion protein not be annotated in Kinase database"), size = 10) +
      xlab(NULL)
    ggsave(paste0(suboutdir, "/error_info.png"), error_plot, width = 15,height = 10,dpi = 300, bg = "white")
  }else{
    colnames(gene.kinase) <- c("Protein","GeneName","log2FC","Kinase.Name")
    gene.kinase.up <- gene.kinase[gene.kinase$log2FC > 0,]
    gene.kinase.up.top <- gene.kinase.up[order(gene.kinase.up$log2FC,decreasing = T),]
    gene.kinase.up.top5 <- gene.kinase.up.top[1:min(5,nrow(gene.kinase.up.top)),]
    gene.kinase.down <- gene.kinase[gene.kinase$log2FC < 0,]
    gene.kinase.down.top <- gene.kinase.down[order(gene.kinase.down$log2FC,decreasing = F),]
    gene.kinase.down.top5 <- gene.kinase.down.top[1:min(5,nrow(gene.kinase.down.top)),]
    gene.kinase.top <- rbind(gene.kinase.up.top5,gene.kinase.down.top5)
    
    write.xlsx(gene.kinase,paste0(suboutdir,"/kinase.gene.xlsx"))
    # 绘制双向柱状图
    p <- ggplot(gene.kinase.top, aes(x = log2FC, y = reorder(Protein, log2FC))) +
      geom_bar(stat = "identity", aes(fill = log2FC > 0), show.legend = TRUE) +
      scale_fill_manual(values = c("TRUE" = "#f5780c", "FALSE" = "#0586c0"), labels = c("Down","Up")) +
      labs(title = "",
           x = "log2(Foldchange)",
           y = "Kinase") +
      theme_classic() +
      theme(axis.text.y = element_text(size = 10, hjust = 1,face = "bold"),  # 调整纵坐标标签的字体大小和对齐方式
            axis.text.x = element_text(size = 10, hjust = 1,face = "bold"),  # 调整横坐标标签的字体大小和对齐方式
            axis.title.x = element_text(size = 15),  # 调整横坐标标签的字体大小
            legend.position = "right",  # 将图例放置在右侧
            legend.title = element_blank(),  # 去除图例标题
            panel.grid.major.y = element_blank(),  # 移除横向网格线
            panel.grid.minor.y = element_blank(),  # 移除次级横向网格线
            axis.title = element_text(color = "black", size = 24)) +
      guides(fill=guide_legend(reverse = T))  
    
    ggsave(paste0(suboutdir, "/kinase-barplot.png"),p,width = 8,height = 5,dpi = 300)
    ggsave(paste0(suboutdir, "/kinase-barplot.pdf"),p,width = 8,height = 5,dpi = 1200)
  }
}






















