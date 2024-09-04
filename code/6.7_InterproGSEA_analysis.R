rm(list = ls())
suppressWarnings(suppressMessages(library(tidyr)))# 1.3.0
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(clusterProfiler))) # 4.6.0
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(openxlsx)))
## 获取比较组信息
get_contrast <- function(file_name) {
  data <- read.xlsx(file_name)
  return(data)
}

## 获取输入文件: 5.2结果文件
get_data <- function(file_name) {
  data <- read.xlsx(file_name)
  return(data)
}

## 获取物种对应文件 --> TERM2GENE
get_KEGG_anno_df <- function(pathway2gene) {
  pathwayInfo <- fread(pathway2gene, sep = "\t", header = F)
  # term2gene_pathway <- pathwayInfo[,c(1,5,7)]
  # term2gene_pathway <- as.data.table(cbind(paste0(toupper(pathwayInfo$V5),"(",pathwayInfo$V1,")"),pathwayInfo[,7]))
  term2gene_pathway <- as.data.table(cbind(paste0(toupper(pathwayInfo$V5), "(", pathwayInfo$V1, ")"), pathwayInfo[, 2]))
  colnames(term2gene_pathway) <- c("pathway", "symbol")
  return(term2gene_pathway)
}

## 获取 genelist
get_genelist <- function(diffExp_df, Anno_df, type) {
  diffExp_df <- diffExp_df %>% dplyr::rename(Protein = Protein.Accessions,Genes = Gene.Name)
  colnames(diffExp_df)[grep(pattern = "Log2fc",x = colnames(diffExp_df))] <- "Logfc"
  genelist_df <- diffExp_df[, c("Protein", "Logfc", "Genes")]
  # genelist_df <- genelist_df[genelist_df$log2FoldChange != "-Inf" & genelist_df$log2FoldChange != "Inf", ] # 去除 -Inf 和 Inf 这个需要再次核对是否是合理的
  genelist <- genelist_df$Logfc
  names(genelist) <- genelist_df$Protein
  
  Anno_df <- as.data.frame(Anno_df)
  genelist <- genelist[names(genelist) %in% unique(Anno_df[,3])]
  
  ## log2fc 由高至低排序
  genelist <- sort(genelist, decreasing = T)
  return(genelist)
}

## GSEA 富集分析
GSEA_anno <- function(suboutdir, contrast, genelist, Anno_GO_df, group_info, minGSGOSize, maxGSGOSize, pvalueCutoff, pAdjustMethod, type) {
  library(BiocParallel)
  # 设置并行计算的线程数
  register(MulticoreParam(workers = 4))
  GSEA_result <- GSEA(genelist,
                      TERM2GENE = as.data.frame(Anno_GO_df[, c("interPro.ID","Uniprot.ID")]),
                      # TERM2NAME = as.data.frame(Anno_GO_df[,c("pathway","ID_name")]),
                      minGSSize = minGSGOSize,
                      maxGSSize = maxGSGOSize,
                      pvalueCutoff = pvalueCutoff,
                      pAdjustMethod = pAdjustMethod,
                      verbose = FALSE,
                      eps = 0      # 设置 eps = 0 表示在计算 ES 时，即使窗口内基因得分之和与基因集内基因得分平均值之差非常小（但不为零），也会累积到 ES 上。换句话说，eps = 0 允许 ES 在任何情况下都可能增加，即使是非常微小的基因得分差异。这种设置可能会使 ES 计算更加敏感，能够捕捉到更细微的基因集富集趋势，但也可能导致 ES 计算结果受到噪声或微小波动的影响，从而降低分析的稳健性。
                      # by = 'fgsea'
                      # nPermSimple = 10000 # 设置为10000表示进行10000次随机排列以获得更精确的P值估计。
  )
  
  if (nrow(as.data.table(GSEA_result)) == 0) {
    error_empty_plot(err_info = "The enrichment analysis results are empty.", outdir)
    quit(save = "no")
  }
  
  save(GSEA_result, file = paste0(suboutdir, "/", contrast, "_", type, "_GSEA_result.RData"))
  GSEA_df <- GSEA_result@result
  match.index <- match(GSEA_df$ID, Anno_GO_df$interPro.ID)
  GSEA_df$Description <- Anno_GO_df$interPro.term[match.index]
  colnames(GSEA_df) <- c("ID", "Description", "SIZE", "ES", "NES", "NOM p-val", "FWER p-val", "FDR q-val", "RANK AT MAX", "LEADING EDGE", "CORE ENRICHMENT") ## 修改 GSEA_df 的列名
  GSEA_df$ID <- sub(".*\\((.*)\\).*", "\\1", GSEA_df$ID)
  library(stringr)
  GSEA_df$Description <- str_replace(GSEA_df$Description, "\\(.*?\\)", "")
  GSEA_df$`CORE ENRICHMENT` <- gsub("/",";",GSEA_df$`CORE ENRICHMENT`)
  
  out_GSEA_df <- GSEA_df[order(GSEA_df$"NOM p-val"), ]
  ## 输出作为 11.4 的页面输入
  write.xlsx(as.data.table(out_GSEA_df[, 1]), paste0(suboutdir, "/", contrast, "_", type, "_Pathway.xlsx"))
  
  GSEA_df_first <- GSEA_df[GSEA_df$ES > 0, ]
  write.xlsx(GSEA_df_first, paste0(suboutdir, "/", contrast, "_", type, "_GSEA_enrichment_", group_info[2], ".xlsx")) ## C_vs_A 中 C的结果 如果是 A_vs_C 则是A的结果
  
  GSEA_df_second <- GSEA_df[GSEA_df$ES < 0, ]
  write.xlsx(GSEA_df_second, paste0(suboutdir, "/", contrast, "_", type, "_GSEA_enrichment_", group_info[1], ".xlsx")) ## C_vs_A 中 A的结果 如果是 A_vs_C 则是C的结果
  
  ## 还需要再画两个图
  
  ## Global ES histogram
  plot_freqpoly <- ggplot(GSEA_df, aes(x = ES)) +
    geom_freqpoly(bins = 10, color = "#ff8787", linewidth = 1) +
    scale_x_continuous(limits = c(-1, 1)) +
    theme_classic() +
    theme(
      panel.grid.major = element_line(color = "grey80", linewidth = 0.1, linetype = 2), # 主网格线
      panel.grid.minor = element_line(color = "grey80", linewidth = 0.1, linetype = 2), # 次网格线
      panel.border = element_rect(color = "black", linewidth = 0.1, fill = NA),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Enrichment Score(ES)", y = "# of protein stes", title = paste0(group_info[1], "_versus_", group_info[2]))
  
  # save_plot(plot_freqpoly, outdir = suboutdir, saveName = paste0(contrast, "_", type, "_Global_ES"), width = 8, height = 7, draw = FALSE)
  # ggsave(paste0(outdir,"/",contrast, "_", type, "_Global_ES.svg"),plot=plot_freqpoly,width = 8, height = 7 )
  
  
  ## Plot of p-values vs. NES
  plot_df <- as.data.frame(GSEA_df[, c("NES", "NOM p-val", "FDR q-val")])
  colnames(plot_df) <- c("NES", "NOM_pval", "FDR_qvalue")
  
  plot_vs <- ggplot(plot_df, aes(x = NES)) +
    geom_point(aes(y = NOM_pval, fill = "Nominal P-value")) +
    geom_point(aes(y = rescale(FDR_qvalue, c(0, max(plot_df$NOM_pval))), color = "FDR_qvalue"), shape = 15) +
    ylab("Nominal P-value") +
    scale_y_continuous(
      limits = c(0, max(plot_df$NOM_pval)), # expand = c(0, 0),
      sec.axis = sec_axis(name = "FDR q-value", ~ rescale(., c(0, max(plot_df$FDR_qvalue))))
    ) +
    labs(title = "NES vs.Significance") +
    theme_classic() +
    theme(
      panel.grid.major = element_line(color = "grey80", linewidth = 0.1, linetype = 2), # 主网格线
      panel.grid.minor = element_line(color = "grey80", linewidth = 0.1, linetype = 2), # 次网格线
      panel.border = element_rect(color = "black", linewidth = 0.1, fill = NA),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    ) +
    scale_fill_manual(values = "black") +
    scale_color_manual(values = "red") +
    guides(color = guide_legend(title = " "), fill = guide_legend(title = " "))
  
  # save_plot(plot_vs, outdir = suboutdir, saveName = paste0(contrast, "_", type, "_p_values_vs_NES"), width = 8, height = 7, draw = FALSE)
  # ggsave(paste0(outdir,"/",contrast, "_", type, "_p_values_vs_NES.svg"),plot=plot_vs,width = 8, height = 7 )
  
  return(GSEA_result)
}

error_empty_plot <- function(err_info = "The enrichment analysis results are empty.", outdir) {
  error_plot <- ggplot() +
    theme_void() +
    geom_text(aes(0, 0, label = err_info), size = 20) +
    xlab(NULL)
  
  ggsave(paste0(outdir, "/error_info.pdf", ), p = error_plot)
  return(error_plot)
}

main <- function(conf) {
  showgene <<- FALSE
  gene_label <<- ""
  # source(conf)
  ## 加载颜色方案
  set_theme_module <- paste0(utils_root, "/set_Theme.R")
  source(set_theme_module)
  ## 加载保存方案
  source(paste0(utils_root, "/figure_Save.R"))
  ## 加载字体方案
  AddFont(font_path)
  showtext_auto()
  
  ## 画图
  set_theme_module <- paste0(utils_root, "/set_Theme.R") ## 加载颜色方案
  source(set_theme_module) ## 加载字体方案
  # source(all_plot_R)
  outdir <- "result/6_GSEA/6.4_Domain_GSEA"
  if (!dir.exists(outdir)) {
    dir.create(outdir,recursive = T)
  } ## 创建输出文件夹
  type <- "Domain"
  contrast_df <- get_contrast(contrast.file)
  # Anno_df <- read.xlsx(Annotation)
  # if (type == "GO" && file.exists(GO_RData)) {
  #   load(GO_RData)
  #   term2gene <- term2gene
  #   Anno_sub_df <- term2gene[, c(1,2,4)]
  #   Anno_sub_df <- na.omit(Anno_sub_df)
  #   Anno_sub_df$pathway <- paste0(toupper(Anno_sub_df$Term), "(", Anno_sub_df$GO_ID, ")")
  #   colnames(Anno_sub_df) <- c("GO_ID", "ID_name","Protein", "pathway")
  # }
  if (type == "Domain") {
    interpro_db_file <- Sys.glob(paste0(interpro_file_path,organism_Scientific_name,"/*tsv"))
    interpro_db <- read.delim(interpro_db_file,header = F)
    colnames(interpro_db)<-c("Uniprot.ID","MD5.disest","length"	,"source.database","database.ID","Description","start","stop",	"e-value"	,"status"	,"Date",	"interPro.ID","interPro.term"	,"GO","Pathway")
    interpro_db$Uniprot.ID <- sapply(strsplit(interpro_db$Uniprot.ID,split = "\\|"),FUN = function(x)x[2])
    interpro_db <- interpro_db %>%
      dplyr::select(interPro.ID, interPro.term, Uniprot.ID) %>%
      dplyr::filter(interPro.ID != "-")
    #interpro_db$Uniprot.ID <- sapply(strsplit(interpro_db$Uniprot.ID, split = "\\|"), FUN = function(x) x[2])
    Anno_sub_df <- interpro_db
  }
  ## 对每个比较组的文件进行循环处理
  for (contrast in contrast_df$contrast) {
    print(contrast)
    suboutdir <- paste0(outdir,"/",contrast)
    if (!dir.exists(suboutdir)) {
      dir.create(suboutdir)
    }
    group_info <- unlist(strsplit(contrast, "_vs_"))
    pattern <- paste0("^", contrast, "-total_results.xlsx")
    diffExp_path <- dir(GSEA_indir, pattern = pattern) %>% paste0(GSEA_indir, .)
    ## 判断 5.2 的结果文件是否存在
    if (file.exists(diffExp_path)) {
      diffExp_df <- read.xlsx(diffExp_path) # 文件: 5.2结果文件
      genelist <- get_genelist(diffExp_df, Anno_sub_df, type) ## 获取全部基因list
      GSEA_result <- GSEA_anno(suboutdir, contrast, genelist, Anno_sub_df, group_info, minGSGOSize, maxGSGOSize, pvalueCutoff, pAdjustMethod, type) ## GSEA 富集分析
      GSEA_df <- GSEA_result@result
      
      ## 获取画图数据 --> 按P值排序的前20个
      pathway_names <- head(GSEA_df[order(GSEA_df$pvalue), "ID"], 20)
      # outdir <<- paste0(outdir,"/p_top20")
      # dir.create(paste0(suboutdir, "/p_top20"), showWarnings = F, recursive = T, mode = "0755") ## 创建输出文件夹
      # for (pathway_name in pathway_names) {
      #   # print(pathway_name)
      #   ## 加载7.2画图脚本
      #   # GSEAplot_all <- GSEAplot_fun(Anno_df, GSEA_result, pathway_name, contrast)
      #   GSEAplot_all <- GSEAplot_fun(Anno_df, GSEA_result, pathway_name, title="",base_size = 11, rel_heights = c(1.5, .5, 1), subplots = 1:3,TRUE, contrast)
      #   ## 保存图片和表格
      #   save_plot_fun(GSEAplot_all, paste0(suboutdir, "/p_top20"), contrast, pathway_name, Anno_df, showgene, gene_label)
      # }
    } else {
      ## 容错信息
      err_info <- paste0(diffExp_path, " 文件不存在")
      fwrite(as.data.table(err_info), file = paste0(outdir, "/err.log"), sep = "\t", quote = F, col.names = FALSE)
    }
  }
  ## 压缩 p_top20 文件夹
  # current_dir <- getwd()
  # setwd(outdir)
  # # tar(paste0(outdir, "/p_top20.tar.gz"), files = paste0(outdir, "/p_top20"), compression = "gzip", tar = "tar", extra = "--remove-files")
  # tar("p_top20.tar.gz", files = "p_top20/", compression = "gzip", tar = "tar", extra = "--remove-files")
  # setwd(current_dir)
}

## 开始运行
#debug for GOGSEA_analysis
conf <- "./config/default.conf"
source(conf)
type <- "Domain"
main(conf)












