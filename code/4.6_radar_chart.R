
#Author: zhaowu zhang
#Date: 2024-07-19 05:28:55
#FilePath: /PERSONALBIO1/prot2/test/proteomics_workflow/SP230608621/protein_workflow_v1.0.0/code/4.6_radar_chart.R
#Description: 
#LastEditTime: 2024-07-30 03:14:00
#Copyright (c) 2024 by <zhaowu.zhang@personalbio.cn>, All Rights Reserved. 

rm(list = ls())

# 加载库
suppressWarnings(suppressMessages(library(fmsb)))
suppressWarnings(suppressMessages(library(circlize)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(circlize)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggradar)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(png)))
suppressWarnings(suppressMessages(library(gridBase)))
suppressWarnings(suppressMessages(library(ggraph)))
suppressWarnings(suppressMessages(library(ComplexHeatmap)))




# 文件路径和参数
conf <- "./config/default.conf"
source(conf)
file_path <- "./result/4_Diff_Expressed/4.1_DiffStats/"
group_file_path <- "./temp/group.xlsx"
contrast_file <- "./temp/contrast.xlsx"
output <- "./result/4_Diff_Expressed/4.6_Radar/"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

group_file <- read.xlsx(group_file_path)
group_names <- c(group_file$SampleName)
contrast_names <- read.xlsx(contrast_file)
column_names <- contrast_names[,1]
column_names_num <- nrow(contrast_names)

group_ids <- group_file$Group
groupID_total <- unique(group_ids)





# 加载颜色方案
set_theme_module <-  "./config/set_Theme.R"
source(set_theme_module)

if (group_colors %in% names(color.ls)) {
  group_colors <- color.ls[[group_colors]]
} else if (group_colors %in% names(color.continuous.ls)) {
  group_colors <- color.continuous.ls[[group_colors]]
}

if (sample_colors %in% names(color.ls)) {
  sample_colors <- color.ls[[sample_colors]]
} else if (sample_colors %in% names(color.continuous.ls)) {
  sample_colors <- color.continuous.ls[[sample_colors]]
}

colors <- group_colors[1:length(groupID_total)] 
group_color <- list(type = setNames(colors, groupID_total))



# 循环读取数据
for (contrast_item in contrast_names$contrast) {
  # contrast_item <- "M1_vs_M2_vs_M3"
  if (!dir.exists(paste0(output, contrast_item))) {
    dir.create(paste0(output, contrast_item), recursive = TRUE)
  }
  
  data <- read.xlsx(paste0(file_path,contrast_item,"-DEP_results.xlsx"))
  
  
  # p值排序
  pvalue_names <- paste0(contrast_item,".pvalue")
  pvalue_col <- match(pvalue_names, names(data))
  data_sorted <- data[order(data[,pvalue_col],decreasing = FALSE), ]
  
  group_col <- match(group_names, names(data_sorted))
  group_col <- na.omit(group_col)
  # 多组对比则不绘制图
  if (length(group_col) < 2){
    next
  }
  data_select <- data_sorted[, group_col]
  rownames(data_select) = data_sorted$Protein.Accessions
  
  #datascale
  # 对数缩放
  data_scale <- log10(data_select + 1)
  data_scale <- as.data.frame(data_scale)
  
  # 筛选当前对比策略的分组
  current_sample <- colnames(data_scale)
  current_group <- group_file[group_file$SampleName %in% current_sample, ]
  group_ids <- current_group$Group
  groupID <- unique(group_ids)
  
  # 计算标准化后样本均值
  # 遍历每个组别，计算均值
  for(group in groupID) {
    # 找到属于当前组的列
    # group <- "M1"
    columns_in_group <- which(group_ids == group)
    # 判断 columns_in_group 的长度
    if (length(columns_in_group) < 2) {
      group_mean <- data_scale[, columns_in_group]
    } else {
      group_mean <- rowMeans(data_scale[, columns_in_group], na.rm = TRUE)
    }
    data_scale[paste0(group, "_Mean")] <- group_mean
  }
    

  # 筛选top_k
  total_num  <- nrow(data_scale)
  if (total_num > 40) {
    top_num <- 40
  } else {
    top_num <- total_num
  }
  
  data_scale_mean <- data_scale[,(grep(".*_Mean$", names(data_scale), value = TRUE))]
  # data_mean_k <- head(data_scale_mean,top_num)
  data_mean_k <- data_scale_mean[1:top_num, ]
  
  # radar数据准备
  radar_data <- t(data_mean_k)
  max_value = round(max(apply(data_mean_k, 2, max)),2)
  min_value = round(min(apply(data_mean_k, 2, min)),2)
  radar_data <- rbind(rep(max_value, ncol(radar_data)),rep(0, ncol(radar_data)), radar_data )
  radar_data <- as.data.frame(radar_data)
  radar_data_NULL <- radar_data
  colnames(radar_data_NULL) = vector("character", 0)
  
  
  # 设置分组颜色
  
  # 选颜色
  # colors_in <- group_color$type[groupID[]]
  colors_in <- c(logo_orange, logo_blue)
  alpha_value <- 0.4

  
  
  colors_in_alpha <- sapply(colors_in, function(color) {
    col <- col2rgb(color, alpha = TRUE) / 255  # 转换为0到1范围的RGB值
    rgb(col[1], col[2], col[3], alpha = alpha_value)
  })
  
  
  
  
  # 创建布局并绘制Circos图
  pdf(paste0(output, contrast_item, "/", "Radar.pdf"), width = 12, height = 12)
  
  #绘制圈图
  # 第一层蛋白ID
  protein_ids <- rownames(data_mean_k)
  gap.degree = 0
  circos.par(gap.degree = gap.degree,
             start.degree = (90 + gap.degree - (360/top_num)/2), # 圈图起始角度在90左边，所以减去扇区间隔角度，与雷达图对齐
             clock.wise = FALSE,
             track.height = 0.15,
             canvas.xlim = c(-1.2, 1.2),   # 调整绘图区域的范围
             canvas.ylim = c(-1.2, 1.2),
             cell.padding = c(0, 0, 0, 0), # 单元格内填充
             track.margin = c(0, 0)  # 轨道之间的间距
  )
  
  circos.initialize(factors = protein_ids, xlim = c(0, 1))
  
  # 添加轨道
  circos.track(ylim = c(0, 1),
               track.height = 0.01,
               bg.border = NA,
               track.margin = c(0.08, 0.08),
               panel.fun = function(x, y) {
                 sector.index <- get.cell.meta.data("sector.index")
                 circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.1, sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
               }, bg.col = NA)
  
  
  # 第二层绘制Log2fc数值
  
  log2fc_cols <- grep("Log2fc$", names(data_sorted))
  data_mean_k$logFC <- data_sorted[[log2fc_cols]][1:top_num]
  data_mean_k$logFC <- round(data_mean_k$logFC, 2)
  circos.track(track.index = 2,
               ylim = c(0, 1),
               track.height = 0.01,
               bg.border = NA,
               track.margin = c(0.03, 0.03),
               
               panel.fun = function(x, y) {
                 sector.index <- get.cell.meta.data("sector.index")
                 log_value <- data_mean_k[sector.index, "logFC"]
                 circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.1,
                             format(log_value, digits = 2),
                             facing = "clockwise",
                             niceFacing = TRUE,
                             adj = c(0, 0.5),
                             col = ifelse(log_value > 0, "#FFB4C0", "#87CEEB"))
               }, bg.col = NA)
  
  
  
  # 第三层绘制差异点大小
  
  circos.track(track.index = 3, ylim = c(min(data_mean_k$logFC), max(data_mean_k$logFC)), panel.fun = function(x, y) {
    sector.index <- get.cell.meta.data("sector.index")
    value <- data_mean_k[sector.index, "logFC"]
    # 根据数值大小调整点的大小
    point_size <- sqrt(abs(value)) * 1.5
    
    # 绘制点，颜色可以根据数值正负变化
    circos.points(x = CELL_META$xcenter, y = rep(0.5, length(value)), pch = 16,
                  cex = point_size, col = ifelse(value > 0, "#FFB4C0", "#87CEEB"))
  }, bg.border = NA)
  
  
  # # 图例颜色映射
  sampleA <- colnames(data_mean_k[1])
  sampleB <- colnames(data_mean_k[2])
  legendA <- expression(Log[2] ~ "FC > 0")
  legendB <- expression(Log[2] ~ "FC < 0")
  
  
  par(new = TRUE,mar = c(12, 12, 12, 12))
  
  radar_grob  <- radarchart(radar_data_NULL,
                            axistype = 0,
                            pcol = NA,
                            pfcol = colors_in_alpha,
                            plwd = 1,
                            cglcol = NA,
                            cglty = 2,
                            axislabcol = NA,
                            caxislabels = NA

  )
  
  
  # 绘制图例
  legend_circos <- Legend(
    labels = c(legendA, legendB, sampleA, sampleB),
    graphics = list(
      function(x = 0, y, w, h) {
        r <- h / 3
        grid.circle(x =  w * 0.1  , y = y , r = r, default.units = "npc",
                    gp = gpar(fill = "white", col = "white"))
        grid.circle(x =  w * 0.1  , y = y , r = r, default.units = "npc",
                    gp = gpar(fill = "#FFB4C0", col = "#FFB4C0"))
      },
      function(x = 0, y, w, h) {
        r <-  r <- h / 3
        grid.circle(x =  w * 0.1  , y = y , r = r, default.units = "npc",
                    gp = gpar(fill = "white", col = "white"))
        grid.circle(x = w * 0.1  , y = y , r = r, default.units = "npc",
                    gp = gpar(fill = "#87CEEB", col = "#87CEEB"))
      },
      function(x = 0, y, w, h) { 
        # # 白色边框
        # grid.rect(x = 0, y, 2.1 * w, 1.1 * h,
        #         gp = gpar(fill = NA, col = "white", lwd = 1.5))
        grid.rect(x = 0, y, w*1.5, h/1.5, gp = gpar(fill = "white", col = "white"))
        grid.rect(x = 0, y, w*1.5, h/1.5, gp = gpar(fill = colors_in_alpha[[1]], col = colors_in_alpha[[1]]))
      },
      function(x = 0, y, w, h) {
        # grid.rect(x = 0, y, 2.1 * w, 1.1 * h,
        #         gp = gpar(fill = NA, col = "white", lwd = 1.5))
        grid.rect(x = 0, y, w*1.5, h/1.5, gp = gpar(fill = "white", col = "white"))
        grid.rect(x = 0, y, w*1.5, h/1.5, gp = gpar(fill = colors_in_alpha[[2]], col = colors_in_alpha[[2]]))
      }
    ),
    labels_gp = gpar(fontsize = 10),
    row_gap = unit(2, "mm"),
    gap = unit(15,"mm")
    
  )
  
  center_legend <- packLegend(legend_circos)
  draw(center_legend, just = "center")
  
  # 结束circos图的绘制
  circos.clear()
  
  # #保存图像
  dev.off()
  
}
