rm(list = ls())
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(enrichplot)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(ggraph))) # 2.1.0
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(Cairo)))
library(ggplot2,lib.loc = "/home/admin/miniconda3/envs/func/lib/R/library/ggplot2/") # 3.4.4
## 加载画图所需文件
get_data <- function(file_name) {
  data <- read.xlsx(file_name)
  return(data)
}

## 画图
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID)) {
    geneSetID <- object@result[geneSetID, "ID"]
  }
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
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

GSEAplot_fun <- function(Anno_df, GSEA_result, geneSetID, title = "", base_size = 11,
                         rel_heights = c(1.5, .5, 1), subplots = 1:3,
                         pvalue_table = FALSE, contrast = "",
                         first_color_GSEA = "BrBG", ES_value = TRUE, ES_font_size = 4, ES_font_color = "#000000", ES_font_family = "Arial", ## ES 标签
                         showgene = FALSE, gene_link_class = "line", gene_font_size = 2, gene_font_family = "Arial", gene_font_color = "#000000",
                         gene_link_color = "#000000", gene_link_size = 0.1, f_line_width = 0.1, ## 上图区 -- 标注基因
                         X_axis_type = "longdash", show_border = TRUE, show_Gridlines = FALSE, Vlink_show_border = TRUE, ## 图表类型
                         up_y_axis_title = "Enrichment score (ES)", up_y_axis_font_size = 12, up_y_axis_font_color = "#000000", up_y_axis_font_family = "Arial", ## 上图区 y 轴
                         Vertical_link_color = "#000000", Vertical_link_width = 0.1, show_colorbar = TRUE, colorbar_number = 8, colorbar_color = "bl2rd", ## 中图区
                         thrid_color = "bl2rd", show_zero = TRUE, zero_font_size = 4, zero_font_family = "Arial", zero_font_color = "#000000", ## 下图区 标记 0 值
                         show_group = TRUE, group_font_size = 4, group_font_family = "Arial", group_font_color = "#000000", ## 下图区 分组标记
                         down_y_axis_title = "Randed List Metric \n (signal to noise)", down_y_axis_font_size = 12, down_y_axis_font_color = "#000000", down_y_axis_font_family = "Arial", ## 下图区 y 轴
                         main_y_axis_label_size = 12, main_y_axis_label_color = "#000000", main_y_axis_label_family = "Arial", ## y 轴刻度
                         title_font_size = 12, title_font_color = "#000000", title_font_family = "Arial", ## 图 标题
                         main_x_axis_title = "Rank in Ordered Dataset", main_x_axis_font_size = 12, main_x_axis_font_color = "#000000", ## x 轴
                         main_x_axis_font_family = "Arial", main_x_axis_label_size = 10, main_x_axis_label_color = "#000000", main_x_axis_label_family = "Arial", ## x 轴 刻度
                         statistic_info_size = 10, statistic_info_color = "#000000", statistic_info_family = "Arial"
) {
  # plot_data <- as.data.frame(GSEA_result@result)
  plot_data <- as.data.frame(GSEA_result)
  rownames(plot_data) <- plot_data[, 1]
  geneList <- position <- NULL ## to satisfy codetool
  gsdata <- gsInfo(GSEA_result, geneSetID)
  max_or_min <- ifelse(plot_data[plot_data$ID == geneSetID, "enrichmentScore"] > 0, "max", "min")
  group_up <- unlist(strsplit(contrast, "_vs_"))[2]
  group_down <- unlist(strsplit(contrast, "_vs_"))[1]
  p <- ggplot(gsdata, aes(x = x)) +
    xlab(NULL) +
    theme_classic(base_size)
  
  ################################################### part 1 上绘图区 ###################################################
  # step1 线条颜色、粗细
  if (first_color_GSEA %in% c("BrBG", "PiYG", "RdBu", "bl2rd", "gnbu", "matlab", "Gwr", "Bwr")) {
    ## 渐变色
    link_color_list <- color.continuous.ls[[first_color_GSEA]]
    es_layer <- geom_line(
      aes(y = runningScore, color = x),
      linewidth = f_line_width, show.legend = FALSE
    )
    p.res <- p + es_layer + scale_edge_colour_gradientn(colours = c(link_color_list))
  } else {
    ## 纯色
    es_layer <- geom_line(
      aes(y = runningScore),
      color = first_color_GSEA,
      linewidth = f_line_width, show.legend = FALSE
    )
    p.res <- p + es_layer
  }
  # step2 是否显示 ES 值
  if (ES_value) {
    if (max_or_min == "max") {
      ES_x_value <- gsdata[which(gsdata$runningScore == max(gsdata$runningScore)), "x"]
      ES_label <- paste0("ES: ", round(max(gsdata$runningScore), 2))
      p.res <- p.res +
        geom_segment(data = gsdata[which(gsdata$runningScore == max(gsdata$runningScore)), ], aes(x = x, y = runningScore, xend = x, yend = 0), linetype = "dotted", linewidth = 0.3) +
        geom_segment(data = gsdata[which(gsdata$runningScore == max(gsdata$runningScore)), ], aes(x = x, y = runningScore, xend = 0, yend = runningScore), linetype = "dotted", linewidth = 0.3) +
        geom_text(data = gsdata[which(gsdata$runningScore == max(gsdata$runningScore)), ], 
                  aes(y = runningScore, x = x + nrow(gsdata)/20), 
                  label = ES_label, color = ES_font_color, size = ES_font_size, family = ES_font_family) +
        theme(
          legend.position = c(.8, .8), legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent")
        ) + scale_x_continuous(expand = c(0, 0))
    } else {
      ES_x_value <- gsdata[which(gsdata$runningScore == min(gsdata$runningScore)), "x"]
      ES_label <- paste0("ES: ", round(min(gsdata$runningScore), 2))
      p.res <- p.res +
        geom_segment(data = gsdata[which(gsdata$runningScore == min(gsdata$runningScore)), ], aes(x = x, y = runningScore, xend = x, yend = 0), linetype = "dotted", linewidth = 0.3) +
        geom_segment(data = gsdata[which(gsdata$runningScore == min(gsdata$runningScore)), ], aes(x = x, y = runningScore, xend = 0, yend = runningScore), linetype = "dotted", linewidth = 0.3) +
        geom_text(data = gsdata[which(gsdata$runningScore == min(gsdata$runningScore)), ], 
                  aes(y = runningScore, x = x - nrow(gsdata)/20), label = ES_label, color = ES_font_color, size = ES_font_size, family = ES_font_family) +
        theme(
          legend.position = c(.8, .8), legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent")
        ) + scale_x_continuous(expand = c(0, 0))
    }
  }
  # step3 是否显示标记基因
  if (showgene && nchar(gene_label) > 0) {
    ## 显示标记基因为TRUE 并且 gene_label 长度大于0 才进行过滤
    ## 对结果做过滤
    
    # gene_label_list <- unlist(strsplit(gene_label,";"))
    ## 基因标签
    gene_label_list <- as.data.table(unlist(strsplit(gene_label, "/")))
    gene_label_list <- as.data.table(matrix(unlist(strsplit(as.character(gene_label_list$V1), ",")), ncol = 2, byrow = TRUE))
    
    gene_label_df <- as.data.frame(gsdata)
    if(!("SwissProtName" %in% colnames(Anno_df))) Anno_df$SwissProtName <- Anno_df$Name
    gene_label_df <- left_join(gene_label_df, Anno_df[, c("Gene ID", "Name", "SwissProtName"), ], by = c("gene" = "Gene ID"))
    gene_label_df <- gene_label_df[which(gene_label_df$gene %in% gene_label_list$V1), ]
    gene_label_df <- left_join(gene_label_df, gene_label_list, by = c("gene" = "V1"))
    
    
    if (gene_link_class == "arrow") {
      ## 连线为箭头
      arrow_content <- arrow(length = unit(0.01, "npc"))
    } else {
      ## 连线为直线
      arrow_content <- NULL
    }
    plot_show_gene <- geom_text_repel(data = gene_label_df, max.overlaps = 40, aes(y = runningScore, x = x, label = V2, ), size = gene_font_size, family = gene_font_family, color = gene_font_color, segment.color = gene_link_color, segment.size = gene_link_size, arrow = arrow_content, min.segment.length = 0)
  } else {
    plot_show_gene <- NULL
  }
  # step4 图表调整 + 标题 + 上y轴
  if (X_axis_type == "longdash") {
    ## 显示X轴的折线（虚线）
    show_line <- geom_hline(yintercept = 0, linetype = "longdash", linewidth = 0.3)
  } else {
    ## 显示X轴的直线 #solid
    show_line <- geom_hline(yintercept = 0, linewidth = 0.3)
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
      panel.grid.major = element_line(color = "grey", linewidth = 0.2, linetype = 2), # 主网格线
      panel.grid.minor = element_line(color = "grey", linewidth = 0.2, linetype = 2), # 次网格线
      panel.border = element_rect(color = "grey", linewidth = 0.2, fill = NA),
      axis.line.x = element_blank(),
      axis.line.y = element_blank()
    ) # 边框
  }
  ## 不添加边框,显示网格线
  if (!show_border & show_Gridlines) {
    add_theme <- theme(
      panel.grid.major = element_line(color = "grey", linewidth = 0.2, linetype = 2), # 主网格线
      panel.grid.minor = element_line(color = "grey", linewidth = 0.2, linetype = 2), # 次网格线
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
  
  p.res <- p.res + plot_show_gene + show_line + add_theme + labs(y = up_y_axis_title) +
    theme(
      axis.title.y = element_text(size = up_y_axis_font_size, colour = up_y_axis_font_color, family = up_y_axis_font_family),
      axis.text.y = element_text(size = main_y_axis_label_size, colour = main_y_axis_label_color, family = main_y_axis_label_family), # y刻度)
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_line(linewidth = 0.1)
    )
  
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
    geom_linerange(data = subset(gsdata, position == 1), aes(ymin = ymin, ymax = ymax), colour = Vertical_link_color, linewidth = Vertical_link_width) +
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
  
  ## 是否显示色块
  if (show_colorbar) {
    # v <- seq(1, sum(gsdata$position), length.out = (colorbar_number-1))
    v <- seq(1, sum(gsdata$position), length.out = (colorbar_number))
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    # print(unique(inv)) ## 查验具体的色块数量
    if (min(inv) == 0) inv <- inv + 1
    if (max(inv) > 8) {
      inv[inv > 8] <- 8
    }
    # col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    pal <- colorRampPalette(color.continuous.ls[[colorbar_color]])
    bar_col <- pal(colorbar_number)
    
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * .3
    xmin <- which(!duplicated(inv))
    # xmin <- ifelse(xmin!=1,1,xmin)
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(
      ymin = ymin, ymax = yy,
      xmin = xmin,
      xmax = xmax,
      col = bar_col[unique(inv)]
    )
    # print(unique(d$col)) ## 查验具体的颜色
    p2 <- p2 + geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = I(col)
      ),
      data = d,
      # alpha=.9,
      inherit.aes = FALSE
    )
  }
  ## 是否显示边框
  if (Vlink_show_border) {
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
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  
  # step1 颜色
  if (thrid_color %in% c("BrBG", "PiYG", "RdBu", "bl2rd", "gnbu", "matlab", "Gwr", "Bwr")) {
    ## 渐变色
    thrid_color_list <- color.continuous.ls[[thrid_color]]
    p.pos <- p + geom_segment(data = df2, aes(x = x, xend = x, y = y, yend = 0, color = x), show.legend = F) +
      scale_color_gradientn(colours = thrid_color_list) + scale_x_continuous(expand = c(0, 0))
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
    df2_max_y <- df2[df2$y == max(df2$y), ]
    df2_min_y <- df2[df2$y == min(df2$y), ]
    if (nrow(df2_max_y) > 1 ){df2_max_y <- df2_max_y[df2_max_y$runningScore == max(df2_max_y$runningScore),]}
    if (nrow(df2_min_y) > 1 ){df2_min_y <- df2_min_y[df2_min_y$runningScore == max(df2_min_y$runningScore),]}
    
    group_df <- rbind(df2_max_y, df2_min_y)
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
      axis.ticks = element_line(linewidth = 0.1),
      axis.line = element_line(colour = "black", linewidth = 0.1, lineend = "square")
    ) +
    theme(plot.margin = margin(t = -.1, r = .2, b = .2, l = .2, unit = "cm")) + add_theme
  
  ################################################### part4 附件区 ###################################################
  
  ## 统计信息
  if (pvalue_table) {
    pd <- plot_data[geneSetID, c("Description", "NES", "pvalue", "p.adjust")]
    rownames(pd) <- pd$Description
    pd <- pd[, -1]
    pd <- as.data.frame(t(pd))
    pd["NES", ] <- paste0("NES: ", round(pd["NES", ], 2))
    pd["pvalue", ] <- ifelse(as.numeric(pd["pvalue", ]) < 0.001, paste0("Pvalue: < ", 0.001),
                             ifelse(as.numeric(pd["pvalue", ]) < 0.01, paste0("Pvalue: < ", 0.01),
                                    ifelse(as.numeric(pd["pvalue", ]) < 0.05, paste0("Pvalue: < ", 0.05),
                                           ifelse(as.numeric(pd["pvalue", ]) < 0.5, paste0("Pvalue: < ", 0.5),
                                                  paste0("Pvalue: > ", 0.5)
                                           )
                                    )
                             )
    )
    
    pd["p.adjust", ] <- ifelse(as.numeric(pd["p.adjust", ]) < 0.001, paste0("Ajusted Pvalue: < ", 0.001),
                               ifelse(as.numeric(pd["p.adjust", ]) < 0.01, paste0("Ajusted Pvalue: < ", 0.01),
                                      ifelse(as.numeric(pd["p.adjust", ]) < 0.05, paste0("Ajusted Pvalue: < ", 0.05),
                                             ifelse(as.numeric(pd["p.adjust", ]) < 0.5, paste0("Ajusted Pvalue: < ", 0.5),
                                                    paste0("Ajusted Pvalue: > ", 0.5)
                                             )
                                      )
                               )
    )
    
    # colnames(pd) <- " "
    subtitle <- paste(pd[1, 1], pd[2, 1], pd[3, 1], sep = ";  ")
    ## 整图 标题
    title <- ifelse(!is.null(title) && !is.na(title) && title != "", title, geneSetID)
    p.res <- p.res + ggtitle(title, subtitle) + theme(
      plot.title = element_text(hjust = 0.5, size = title_font_size, colour = title_font_color, family = title_font_family),
      plot.subtitle = element_text(hjust = 0.5, size = statistic_info_size, colour = statistic_info_color, family = statistic_info_family)
    )
    
    
    # tp <- tableGrob(pd, rows = NULL ,cols = NULL ,theme = ttheme_default(base_size = 10, base_colour = "black", base_family = "Arial"))
    # p.res <- p.res + theme(legend.position = "none") +
    #   annotation_custom(tp,
    #     xmin = quantile(p.res$data$x, .5),
    #     xmax = quantile(p.res$data$x, 1),
    #     ymin = quantile(p.res$data$runningScore, .75),
    #     ymax = quantile(p.res$data$runningScore, 1)
    #   )
  }
  
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
  return(list(aplot::gglist(gglist = plotlist, nrow = 3, heights = rel_heights), ES_x_value, max_or_min, group_up, group_down))
}

save_plot_fun <- function(GSEAplot_all, outdir, contrast, pathway_name, Anno_df, showgene, gene_label) {
  library(stringr)
  pathway_id <- pathway_name
  # pathway_id <- gsub(":", "_", pathway_id)
  GSEAplot <- GSEAplot_all[[1]]
  ES_x_value <- GSEAplot_all[[2]]
  max_or_min <- GSEAplot_all[[3]]
  group_up <- GSEAplot_all[[4]]
  group_down <- GSEAplot_all[[5]]
  ## 保存图片
  GSEAplot1 <- GSEAplot[[1]]
  GSEAplot2 <- GSEAplot[[2]]
  GSEAplot3 <- GSEAplot[[3]]
  combined_plot <- GSEAplot1 / GSEAplot2 / GSEAplot3
  
  save_group <- ifelse(max_or_min == "max", group_up, group_down)
  save_plot(combined_plot, outdir = outdir, saveName = paste0(contrast, "_(",pathway_id, ")_", save_group), width = 10, height = 8, draw = FALSE)
  # ggsave(paste0(outdir,"/",contrast, pathway_id, "_", save_group,".svg"),plot=combined_plot,width = 10, height = 8 )
  
  
  ## 输出表格
  GSEA_pathway_df <- GSEAplot3$data[, c("x", "gene", "gene", "geneList", "runningScore", "position")]
  GSEAplot3_info <- ggplot_build(GSEAplot3)
  
  GSEA_pathway_df$"RANK METRIC SCORE" <- GSEAplot3_info$data[[1]]$y
  GSEA_pathway_df <- GSEA_pathway_df[GSEA_pathway_df$position == 1, ]
  GSEA_pathway_df$position <- NULL
  GSEA_pathway_df$"CORE ENRICHMENT" <- ifelse(GSEA_pathway_df$x <= ES_x_value, "Yes", "No")
  GSEA_pathway_df$x <- NULL
  colnames(GSEA_pathway_df) <- c("SYMBOL", "TITLE", "RANK IN GENE LIST", "RUNNING ES", "RANK METRIC SCORE", "CORE ENRICHMENT")
  colname_Anno_df <- colnames(Anno_df)
  # if(!("SwissProtName" %in% colnames(Anno_df))) Anno_df$SwissProtName <- Anno_df$symbol
  # GSEA_pathway_df <- left_join(GSEA_pathway_df, Anno_df[, c("pathway", "symbol", "SwissProtName"), ], by = c("SYMBOL" = "symbol"))
  ## 修改顺序
  # colnames(GSEA_pathway_df)[c(7,8)] <- c("Name","SwissProtName")
  # GSEA_pathway_df <- GSEA_pathway_df[, c("Name", "TITLE", "SYMBOL", "SwissProtName", "RANK IN GENE LIST", "RANK METRIC SCORE", "RUNNING ES", "CORE ENRICHMENT")]
  GSEA_pathway_df$TITLE <- NULL
  colnames(GSEA_pathway_df) <- c("Protein", "RANK IN GENE LIST", "RANK METRIC SCORE", "RUNNING ES", "CORE ENRICHMENT")
  # if (showgene & nchar(gene_label) > 0) {
  #   ## 基因标签
  #   gene_label_list <- as.data.table(unlist(strsplit(gene_label, "/")))
  #   gene_label_list <- as.data.table(matrix(unlist(strsplit(as.character(gene_label_list$V1), ",")), ncol = 2, byrow = TRUE))
  #   colnames(gene_label_list) <- c("Gene ID", "ShowLabel")
  #   GSEA_pathway_df <- left_join(GSEA_pathway_df, gene_label_list, by = "Gene ID")
  #   # GSEA_pathway_df <-  GSEA_pathway_df[which(GSEA_pathway_df$gene %in% gene_label_list$V1),]
  # } else {
  #   GSEA_pathway_df <- GSEA_pathway_df %>% mutate(ShowLabel = SYMBOL)
  # }
  # colnames(GSEA_pathway_df)[1] <- "Gene Name"
  GSEA_pathway_df <- GSEA_pathway_df[order(GSEA_pathway_df$"RANK METRIC SCORE", decreasing = TRUE), ]
  
  # if(!("SwissProtName" %in% colname_Anno_df)) GSEA_pathway_df$SwissProtName <- NULL
  
  write.xlsx(GSEA_pathway_df, paste0(outdir, "/", contrast, "_(",pathway_id, ")_", save_group, ".xlsx"))
  
  ## 对应空白图片和空表的生成
  # empty_group <- ifelse(max_or_min == "max", group_down, group_up)
  # empty_plot <- ggplot() +
  #   theme_void() +
  #   geom_text(aes(0, 0, label = "The current group does not contain the\ncurrently pathway name."), size = 10) +
  #   xlab(NULL)
  # save_plot(empty_plot, outdir = outdir, saveName = paste0(contrast, pathway_id, "_", empty_group), width = 10, height = 8, draw = FALSE)
  # 
  # empty_table <- data.table(matrix(ncol = 8, nrow = 0))
  # colnames(empty_table) <- c("Gene Name", "Gene ID", "SwissProtName", "RANK IN GENE LIST", "RANK METRIC SCORE", "RUNNING ES", "CORE ENRICHMENT", "ShowLabel")
  # if(!("SwissProtName" %in% colname_Anno_df)) empty_table$SwissProtName <- NULL
  # fwrite(empty_table, paste0(outdir, "/", contrast, pathway_id, "_", empty_group, ".xls"), sep = "\t")
}

main <- function(conf) {
  # source(conf)
  ## 加载颜色方案
  set_theme_module <- paste0(utils_root, "/set_Theme.R")
  source(set_theme_module)
  ## 加载保存方案
  source(paste0(utils_root, "/figure_Save.R"))
  ## 加载字体方案
  AddFont(font_path)
  showtext_auto()
  outdir <- "result/6_GSEA/6.3_Reactome_GSEA"
  if (!dir.exists(outdir)) {
    dir.create(outdir,recursive = T)
  } ## 创建输出文件夹
  type <- "Reactome"
  contrast_df <- get_data(contrast.file)
  Anno_df <- read.xlsx(Annotation)
  for (contrast in contrast_df[[1]]) { ## 对每个比较组的文件进行循环处理
    print(contrast)
    suboutdir <- paste0(outdir,"/",contrast)
    if (!dir.exists(suboutdir)) {
      dir.create(suboutdir)
    }
    subname <- ".*.GSEA_result.RData"
    # strsplit(indir,"/")
    if(type == "GO") subname <- "GO_GSEA_result.RData"
    if(type == "KEGG") subname <- "KEGG_GSEA_result.RData"
    if(type == "Reactome") subname <- "Reactome_GSEA_result.RData"
    ## 加载文件
    pattern <- paste0(contrast, "_",subname)
    GSEA_result_path <- paste0(suboutdir,"/",pattern)
    # print(paste0("GSEA_result_path:",GSEA_result_path))
    load(GSEA_result_path)
    
    if(!exists("pathway_name") || nchar(pathway_name)==0){
      GSEA_df <- GSEA_result@result
      # 按P值排序的 top1
      pathway_names <- head(GSEA_df[order(GSEA_df$pvalue),"ID"],20)
    }         
    dir.create(paste0(suboutdir, "/p_top20"), showWarnings = F, recursive = T, mode = "0755") ## 创建输出文件夹
    ## 画图
    for (pathway.name in pathway_names) {
      # print(pathway.name)
      GSEAplot_all <- GSEAplot_fun(Anno_df, GSEA_result, pathway.name, title,
                                   base_size = 11, rel_heights = c(1.5, .5, 1), subplots = 1:3,
                                   statistic_info, contrast,
                                   first_color_GSEA, ES_value, ES_font_size, ES_font_color, ES_font_family, ## ES 标签
                                   showgene, gene_link_class, gene_font_size, gene_font_family, gene_font_color, gene_link_color, gene_link_size, f_line_width, ## 上图区 -- 标注基因
                                   X_axis_type, show_border, show_Gridlines, Vlink_show_border, ## 图表类型
                                   up_y_axis_title, up_y_axis_font_size, up_y_axis_font_color, up_y_axis_font_family, ## 上图区 y 轴
                                   Vertical_link_color, Vertical_link_width, show_colorbar, colorbar_number, colorbar_color, ## 中图区
                                   thrid_color, show_zero, zero_font_size, zero_font_family, zero_font_color, ## 下图区 标记 0 值
                                   show_group, group_font_size, group_font_family, group_font_color, ## 下图区 分组标记
                                   down_y_axis_title, down_y_axis_font_size, down_y_axis_font_color, down_y_axis_font_family, ## 下图区 y 轴
                                   main_y_axis_label_size, main_y_axis_label_color, main_y_axis_label_family, ## y 轴刻度
                                   title_font_size, title_font_color, title_font_family, ## 图 标题
                                   main_x_axis_title, main_x_axis_font_size, main_x_axis_font_color, ## x 轴
                                   main_x_axis_font_family, main_x_axis_label_size, main_x_axis_label_color, main_x_axis_label_family, ## x 轴 刻度
                                   statistic_info_size, statistic_info_color, statistic_info_family
      )
      ## 保存图片和表格
      save_plot_fun(GSEAplot_all, paste0(suboutdir, "/p_top20"), contrast, pathway.name, Anno_df, showgene, gene_label)
    }
    
    # rm(pathway_name)
    ## 压缩 p_top20 文件夹
    # current_dir <- getwd()
    # setwd(suboutdir)
    # # tar(paste0(outdir, "/p_top20.tar.gz"), files = paste0(outdir, "/p_top20"), compression = "gzip", tar = "tar", extra = "--remove-files")
    # tar("p_top20.tar.gz", files = "p_top20/", compression = "gzip", tar = "tar", extra = "--remove-files")
    # setwd(current_dir)
  }
}

# 开始运行
# argv <- commandArgs(TRUE)
# if (length(argv) == 1) {
#   source(argv[1])
#   if (!exists("all_plot_R")) {
#     main(argv[1])
#   }
# }

## debug for GOGSEAPlot analysis
conf <- "config/default.conf"
source(conf)
indir <- "result/6_GSEA/6.3_Reactome_GSEA/"
Annotation <- Sys.glob(paste0(uniprot_file_path,organism_category,"/*",organism_code, "_", organism_ID, "_annotation.xlsx"))
type <- "Reactome"
main(conf)



