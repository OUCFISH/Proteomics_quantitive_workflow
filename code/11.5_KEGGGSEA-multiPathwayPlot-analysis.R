rm(list = ls())
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(enrichplot)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(ggraph))) # 2.1.0
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(openxlsx)))
## 加载画图所需文件
get_data <- function(file_name) {
  data <- fread(file_name, header = T, sep = "\t")
  return(data)
}

gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  if (length(object@gene2Symbol) == 0) {
    df$gene <- names(geneList)
  } else {
    df$gene <- object@gene2Symbol[names(geneList)]
  }
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

gseaScores <- getFromNamespace("gseaScores", "DOSE")

GSEAplot_fun <- function(temp_df,GSEA_result, geneSetID, title = "", base_size = 11,
                         rel_heights = c(1.5, .5, 1), subplots = 1:3,
                         pvalue_table = FALSE, contrast = "",
                         first_color, ES_value = TRUE, ES_font_size = 11, ES_font_color = "#000000", ES_font_family = "Arial", ## ES 标签
                         showgene = TRUE, gene_link_class = "line", gene_font_size = 12, gene_font_family = "Arial", gene_font_color = "#000000",
                         gene_link_color = "#000000", gene_link_size = 0.1, ## 上图区 -- 标注基因
                         X_axis_type = "longdash", show_border = TRUE, show_Gridlines = FALSE, ## 图表类型
                         up_y_axis_title = "Enrichment score (ES)", up_y_axis_font_size = 12, up_y_axis_font_color = "#000000", up_y_axis_font_family = "Arial", ## 上图区 y 轴
                         Vertical_link_color = "#000000", Vertical_link_width = 0.1, show_colorbar = TRUE, colorbar_number = 8, colorbar_color = "bl2rd", ## 中图区
                         thrid_color = "bl2rd", show_zero = TRUE, zero_font_size = 4, zero_font_family = "Arial", zero_font_color = "#000000", ## 下图区 标记 0 值
                         show_group = TRUE, group_font_size = 10, group_font_family = "Arial", group_font_color = "#000000", ## 下图区 分组标记
                         down_y_axis_title = "Randed List Metric \n (signal to noise)", down_y_axis_font_size = 12, down_y_axis_font_color = "#000000", down_y_axis_font_family = "Arial", ## 下图区 y 轴
                         main_y_axis_label_size = 12, main_y_axis_label_color = "#000000", main_y_axis_label_family = "Arial", ## y 轴刻度
                         title_font_size = 12, title_font_color = "#000000", title_font_family = "Arial", ## 图 标题
                         main_x_axis_title = "Rank in Ordered Dataset", main_x_axis_font_size = 12, main_x_axis_font_color = "#000000", ## x 轴
                         main_x_axis_font_family = "Arial", main_x_axis_label_size = 12, main_x_axis_label_color = "#000000", main_x_axis_label_family = "Arial", ## x 轴 刻度
                         legend_family,legend_color,legend_size_GSEA
){
  
  plot_data <- as.data.frame(GSEA_result@result)
  rownames(plot_data) <- plot_data[,1]
  geneList <- position <- NULL ## to satisfy codetool
  
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(GSEA_result, geneSetID)
  } else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = GSEA_result))
  }
  gsdata <- left_join(gsdata,temp_df,by=c("Description"="pathway_name"))
  gsdata$Description <- NULL
  gsdata <- gsdata %>% rename(Description=new_name)
  
  ## 统计信息
  if (pvalue_table) {
    pd <- plot_data[geneSetID, c("Description", "NES","pvalue", "p.adjust")]
    rownames(pd) <- pd$Description
    # pd <- pd[, -1]
    pd$NES <- paste0("NES: ",round(pd$NES,2))
    pd$pvalue <- ifelse(as.numeric(pd$pvalue) < 0.001, paste0("Pvalue: < ",0.001),
                        ifelse(as.numeric(pd$pvalue) < 0.01,paste0("Pvalue: < ",0.01),
                               ifelse(as.numeric(pd$pvalue) < 0.05,paste0("Pvalue: < ",0.05), 
                                      ifelse(as.numeric(pd$pvalue) < 0.5,paste0("Pvalue: < ",0.5),
                                             paste0("Pvalue: > ",0.5)))))
    pd$"p.adjust" <- ifelse(as.numeric(pd$"p.adjust") < 0.001, paste0("Ajusted Pvalue: < ",0.001),
                            ifelse(as.numeric(pd$"p.adjust") < 0.01,paste0("Ajusted Pvalue: < ",0.01),
                                   ifelse(as.numeric(pd$"p.adjust") < 0.05,paste0("Ajusted Pvalue: < ",0.05), 
                                          ifelse(as.numeric(pd$"p.adjust") < 0.5,paste0("Ajusted Pvalue: < ",0.5),
                                                 paste0("Ajusted Pvalue: > ",0.5)))))
    
    # pd <- as.data.frame(t(pd))
    pd$label <- paste0(pd$NES,"; ",pd$pvalue,"; ",pd$"p.adjust")
    pd$NES <- pd$pvalue <- pd$"p.adjust" <- NULL
    # print(pd)
    pvalue_df <- full_join(gsdata, pd,by = join_by(Description))
    pvalue_df$label <- paste0(pvalue_df$Description,"\n",pvalue_df$label)
    gsdata <- pvalue_df
    
    # step1 线条颜色、粗细
    es_layer <- geom_line(
      aes(y = runningScore, color = label),
      linewidth = f_line_width, show.legend = TRUE
    )
  }else{
    # step1 线条颜色、粗细
    es_layer <- geom_line(
      aes(y = runningScore, color = Description),
      linewidth = f_line_width, show.legend = TRUE
    )
  }
  p <- ggplot(gsdata, aes(x = x)) +
    xlab(NULL) +
    theme_classic(base_size) +
    # theme(panel.grid.major = element_line(colour = "grey92"),
    #       panel.grid.minor = element_line(colour = "grey92"),
    #       panel.grid.major.y = element_blank(),
    #       panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand=c(0,0))  
  
  link_color_list <- color.continuous.ls[[first_color_GSEA]]
  ################################################### part 1 上绘图区 ###################################################
  
  p.res <- p + es_layer + scale_edge_colour_gradientn(colours = c(link_color_list))
  
  # step2 是否显示 ES 值
  if (ES_value) {
    ES_label_df <- gsdata %>% 
      group_by(Description) %>% 
      summarise(
        runningScore_value = ifelse(sum(runningScore > 0) > sum(runningScore < 0), 
                                    max(runningScore), min(runningScore)),
        x_value = ifelse(sum(runningScore > 0) > sum(runningScore < 0), 
                         x[which.max(runningScore)], x[which.min(runningScore)])
      )
    
    p.res <- p.res + 
      geom_text_repel(data = ES_label_df, 
                      aes(y = runningScore_value, x = x_value, label = paste0("ES:",round(runningScore_value,2))),
                      size = ES_font_size, family = ES_font_family, color = ES_font_color, segment.color = ES_link_color, segment.size = ES_link_size) 
  }
  
  # step3 图表调整 + 标题 + 上y轴
  if (X_axis_type == "longdash") {
    ## 显示X轴的折线（虚线）
    show_line <- geom_hline(yintercept = 0, linetype = "longdash",linewidth = 0.3)
  } else {
    ## 显示X轴的直线 #solid
    show_line <- geom_hline(yintercept = 0,linewidth = 0.3)
  }
  ## 添加边框,不显示网格线
  if (show_border & !show_Gridlines) {
    add_theme <- theme(
      panel.grid.major = element_blank(), # 主网格线
      panel.grid.minor = element_blank(), # 次网格线
      panel.border = element_rect(color = "grey", linewidth = 0.2, fill = NA),
      axis.line.x = element_blank(),
      axis.line.y = element_blank()
    ) # 边框
  }
  ## 添加边框,显示网格线
  if (show_border & show_Gridlines) {
    ## 显示边框，显示网格线
    add_theme <- theme(
      panel.grid.major = element_line(color = "grey", linewidth = 0.2,linetype=2), # 主网格线
      panel.grid.minor = element_line(color = "grey", linewidth = 0.2,linetype=2), # 次网格线
      panel.border = element_rect(color = "grey", linewidth = 0.2, fill = NA),
      axis.line.x = element_blank(),
      axis.line.y = element_blank()
    ) # 边框
  }
  ## 不添加边框,显示网格线
  if (!show_border & show_Gridlines) {
    add_theme <- theme(
      panel.grid.major = element_line(color = "grey", linewidth = 0.2,linetype=2), # 主网格线
      panel.grid.minor = element_line(color = "grey", linewidth = 0.2,linetype=2), # 次网格线
      panel.border = element_blank(),
      axis.line.x = element_line(linewidth = 0.1),
      axis.line.y = element_line(linewidth = 0.1)
    ) # 边框
  }
  ## 不添加边框,不显示网格线
  if (!show_border & !show_Gridlines) {
    add_theme <- theme(
      panel.grid.major = element_blank(), # 主网格线
      panel.grid.minor = element_blank(), # 次网格线
      panel.border = element_blank(),
      axis.line.x = element_line(linewidth = 0.1),
      axis.line.y = element_line(linewidth = 0.1)
    ) # 边框
  }
  
  p.res <- p.res + show_line + add_theme + labs(y = up_y_axis_title) +
    theme(axis.title.y = element_text(size = up_y_axis_font_size, 
                                      colour = up_y_axis_font_color, 
                                      family = up_y_axis_font_family),
          axis.text.y = element_text(size = main_y_axis_label_size, 
                                     colour = main_y_axis_label_color, 
                                     family = main_y_axis_label_family), # y刻度)
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank() ,
          axis.ticks.y= element_line(linewidth = 0.1),
          legend.key.width = unit(1.5,"cm"),
          legend.key.height = unit(1.2,"cm"),
    ) + guides(color = guide_legend(title = "")) +
    theme(legend.text = element_text(family  = legend_family,size = legend_size_GSEA,color = legend_color))
  
  #   ## 整图 标题
  # title <- ifelse(!is.null(title) && !is.na(title) && title != "", title, geneSetID)
  # p.res <- p.res + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, size = title_font_size, colour = title_font_color, family = title_font_family))
  
  
  ################################################### part2 中绘图区 ###################################################
  
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  
  p2 <- ggplot(gsdata, aes(x = x)) +
    geom_linerange(data = subset(gsdata, position == 1), aes(ymin = ymin, ymax = ymax, color=Description), linewidth = Vertical_link_width,show.legend = FALSE) +
    scale_edge_colour_gradientn(colours = c(link_color_list)) + 
    xlab(NULL) +
    ylab(NULL) +
    theme_classic(base_size) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = -.1, b = 0, unit = "cm"),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line.x = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) 
  
  ## 是否显示边框
  if(Vlink_show_border){
    add_theme2 <- theme(
      panel.grid.major = element_blank(), # 主网格线
      panel.grid.minor = element_blank(), # 次网格线
      panel.border = element_rect(color = "grey", linewidth = 0.1, fill = NA),
      axis.line.x = element_blank(),
      axis.line.y = element_blank()
    ) # 边框
    p2 <- p2 + add_theme2
  }
  
  ################################################### part3 下绘图区 ###################################################
  
  df2 <- p$data[p$data$Description==temp_df$new_name[1],]
  # df2 <- p$data[p$data$Description==geneSetID[1],]
  # df2 <- p$data 
  df2$y <- p$data$geneList[df2$x]
  
  # step1 颜色
  if (thrid_color %in% c("BrBG", "PiYG", "RdBu", "bl2rd", "gnbu", "matlab", "Gwr", "Bwr")) {
    ## 渐变色
    thrid_color_list <- color.continuous.ls[[thrid_color]]
    p.pos <- p + geom_segment(data = df2, aes(x = x, xend = x, y = y, yend = 0, color = x), show.legend = F) +
      scale_color_gradientn(colours = thrid_color_list)
  } else {
    ## 纯色
    p.pos <- p + geom_segment(data = df2, aes(x = x, xend = x, y = y, yend = 0), color = thrid_color, show.legend = F)
  }
  
  ## 显示0值标记
  if (show_zero) {
    p.pos <- p.pos + geom_vline(xintercept = mean(df2[df2$y == df2$y[which.min(abs(df2$y))], "x"]), linetype = "dotted", linewidth = 0.3) +
      geom_text(
        data = df2[df2$x == df2[df2$y == df2$y[which.min(abs(df2$y))], "x"][1], ],
        aes(x = mean(x), y = y - 1, label = paste0("Zero cross at ", mean(x))),
        color = zero_font_color, size = zero_font_size, family = zero_font_family
      )
  }
  
  ## 显示分组标记
  if (show_group) {
    # contrast
    group_up <- unlist(strsplit(contrast, "_vs_"))[1]
    group_down <- unlist(strsplit(contrast, "_vs_"))[2]
    
    df2_max_y <- df2[df2$y == max(df2$y), ]
    df2_min_y <- df2[df2$y == min(df2$y), ]
    if (nrow(df2_max_y) > 1 ){df2_max_y <- df2_max_y[df2_max_y$runningScore == max(df2_max_y$runningScore),]}
    if (nrow(df2_min_y) > 1 ){df2_min_y <- df2_min_y[df2_min_y$runningScore == max(df2_min_y$runningScore),]}
    
    group_df <- rbind(df2_max_y, df2_min_y)
    
    
    # group_df <- rbind(df2[df2$y == max(df2$y), ], df2[df2$y == min(df2$y), ])
    group_df$x <- c(group_df$x[2] / 6, group_df$x[2] / 6 * 5)
    group_df$y <- group_df$y / 2
    
    group_df$label <- c(paste0(group_up, " (positively correlated)"), paste0(group_down, " (negtively correlated)"))
    group_text <- geom_text(
      data = group_df, aes(x = x, y = y, label = label),
      color = group_font_color, size = group_font_size, family = group_font_family
    )
  } else {
    group_text <- NULL
  }
  
  p.pos <- p.pos + group_text + show_line
  
  ## 标题设置
  p.pos <- p.pos + labs(y = down_y_axis_title, x = main_x_axis_title) +
    theme(
      axis.title.y = element_text(size = down_y_axis_font_size, colour = down_y_axis_font_color, family = down_y_axis_font_family), # y轴 标题
      axis.text.y = element_text(size = main_y_axis_label_size, colour = main_y_axis_label_color, family = main_y_axis_label_family), # y刻度 
      axis.title.x = element_text(size = main_x_axis_font_size, colour = main_x_axis_font_color, family = main_x_axis_font_family), # x轴 标题
      axis.text.x = element_text(size = main_x_axis_label_size, colour = main_x_axis_label_color, family = main_x_axis_label_family), # x刻度 
      axis.ticks= element_line(linewidth = 0.1),
      axis.line = element_line(colour = "black", linewidth = 0.1, lineend = "square")
    ) + 
    theme(plot.margin = margin(t = -.1, r = .2, b = .2, l = .2, unit = "cm")) + add_theme 
  
  ## 三个图组合
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(
      axis.line.x = element_line(),
      axis.ticks.x = element_line(),
      axis.text.x = element_text()
    )
  
  if (length(subplots) == 1) {
    return(plotlist[[1]] + theme(plot.margin = margin(
      t = .2, r = .2, b = .2,
      l = .2, unit = "cm"
    )))
  }
  
  if (length(rel_heights) > length(subplots)) {
    rel_heights <- rel_heights[subplots]
  }
  
  return(aplot::gglist(gglist = plotlist, ncol = 1, heights = rel_heights))
}

save_plot_fun <- function(GSEAplot_all,suboutdir,contrast,pathway_name){
  GSEAplot <- GSEAplot_all
  ## 保存图片
  GSEAplot1 <- GSEAplot[[1]]
  GSEAplot2 <- GSEAplot[[2]]
  GSEAplot3 <- GSEAplot[[3]]
  combined_plot <- GSEAplot1 / GSEAplot2 / GSEAplot3
  
  
  save_plot(combined_plot, outdir = suboutdir, saveName = paste0(contrast,"_GSEA"), width = 11, height = 8, draw = FALSE)
}

main <- function(conf) {
  # source(conf)
  ## 加载颜色方案
  source(paste0(utils_root,"/set_Theme.R")) 
  ## 加载字体方案
  AddFont(font_path)
  showtext_auto()
  ## 加载保存方案
  source(paste0(utils_root,"/figure_Save.R")) 
  outdir <- "result/11_AdcancedPlot/11.4_GSEA-multiPathwayPlot/KEGG"
  if (!dir.exists(outdir)) {
    dir.create(outdir,recursive = T)
  } ## 创建输出文件夹
  contrast_df <- read.xlsx(contrast.file)
  for(contrast in contrast_df[[1]]){ ## 对每个比较组的文件进行循环处理
    print(contrast)
    suboutdir <- paste0(outdir,"/",contrast,"/")
    if (!dir.exists(suboutdir)) {
      dir.create(suboutdir)
    }
    subname <- ".*.GSEA_result.RData"
    if(grepl("GO",unlist(strsplit(indir,"/"))[length(unlist(strsplit(indir,"/")))])) subname <- "GO_GSEA_result.RData"
    if(grepl("KEGG",unlist(strsplit(indir,"/"))[length(unlist(strsplit(indir,"/")))])) subname <- "KEGG_GSEA_result.RData"
    pattern <- paste0(contrast,"_",subname)
    RData <- dir(paste0(indir,"/",contrast,"/"),pattern = pattern)
    RData_path <- dir(paste0(indir,"/",contrast,"/"),pattern = pattern) %>% paste0(paste0(indir,"/",contrast),"/",.)
    # print(RData_path)
    if (length(RData) != 0) {
      load(RData_path)
      
      if (!exists("pathway_name")) {
        pathway_name <- as.data.table(GSEA_result) %>%
          slice_min(pvalue, n = 5) %>%
          dplyr::select(ID)
        pathway_name <- paste(pathway_name$ID, collapse = ";")     
      } 
      pathway_name <- unlist(strsplit(pathway_name, ";"))
      temp_df <- as.data.table(pathway_name)
      
      temp_df$new_name <- apply(temp_df,1,function(x) {
        namestr <- substr(x, 1, gregexpr("\\(", x)[[1]][1]-1)
        idstr <- substr(x, gregexpr("\\(", x)[[1]][1],nchar(x))
        if(nchar(namestr) > 50) namestr <- paste0(substr(namestr, 1, 50),"...")
        new_pathname <- paste0(namestr,idstr)
        return(new_pathname)
      })
      
      ## 画图
      GSEAplot_all <- GSEAplot_fun(temp_df,GSEA_result, # 画图数据
                                   pathway_name, title,
                                   base_size = 11, rel_heights = c(1.5, .5, 1), subplots = 1:3,
                                   statistic_info, contrast = contrast,
                                   first_color, ES_value, ES_font_size, ES_font_color, ES_font_family, ## ES 标签
                                   showgene, gene_link_class, gene_font_size, gene_font_family, gene_font_color, gene_link_color, gene_link_size, ## 上图区 -- 标注基因
                                   X_axis_type, show_border, show_Gridlines, ## 图表类型
                                   up_y_axis_title, up_y_axis_font_size, up_y_axis_font_color, up_y_axis_font_family, ## 上图区 y 轴
                                   Vertical_link_color, Vertical_link_width, show_colorbar, colorbar_number, colorbar_color, ## 中图区
                                   thrid_color, show_zero, zero_font_size, zero_font_family, zero_font_color, ## 下图区 标记 0 值
                                   show_group, group_font_size, group_font_family, group_font_color, ## 下图区 分组标记
                                   down_y_axis_title, down_y_axis_font_size, down_y_axis_font_color, down_y_axis_font_family, ## 下图区 y 轴
                                   main_y_axis_label_size, main_y_axis_label_color, main_y_axis_label_family, ## y 轴刻度
                                   title_font_size, title_font_color, title_font_family, ## 图 标题
                                   main_x_axis_title, main_x_axis_font_size, main_x_axis_font_color, ## x 轴
                                   main_x_axis_font_family, main_x_axis_label_size, main_x_axis_label_color, main_x_axis_label_family, ## x 轴 刻度
                                   legend_family,legend_color,legend_size_GSEA
      ) 
      ## 保存图片和表格
      save_plot_fun(GSEAplot_all, suboutdir, contrast, pathway_name)
      rm(pathway_name)
    }else{
      writeLines("no_result",paste0(suboutdir, "no_result.txt"))
    }
  }
  
}

## 开始运行
## debug for GSEA_multiKEGGPathway plot
conf <- "config/default.conf"
source(conf)
indir <- "result/6_GSEA/6.2_KEGG_GSEA"
main(conf)


