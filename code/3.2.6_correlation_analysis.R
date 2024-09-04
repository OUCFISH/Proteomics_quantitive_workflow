#!/usr/bin/env Rscript
rm(list = ls())
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(pheatmap))) # pheatmap:::cluster_mat
suppressWarnings(suppressMessages(library(ape)))
suppressWarnings(suppressMessages(library(igraph)))
suppressWarnings(suppressMessages(library(rjson)))
suppressWarnings(suppressMessages(library(openxlsx)))
scale_log10 <- function(x) {
  x <- log10(x + 1)
  return(x)
}

scale_mat <- function(data, scale) {
  if (!(scale %in% c("none", "log10"))) {
    stop("scale argument shoud take values: 'none' or 'log10'")
  }
  data <- switch(scale,
    none = data,
    log10 = scale_log10(data)
  )
  return(data)
}

make_group <- function(mat, map) {
  mat <- cbind(t(mat), map[, 3])
  mat <- aggregate(mat[, 1:ncol(mat) - 1], by = list(Group = map$Default), FUN = mean)
  rownames(mat) <- mat$Group
  return(t(mat[, 2:ncol(mat)]))
}

# process_mat <- function(mat, map) {
#   mat <- mat[, as.vector(map$SampleID)]
#   colnames(mat) <- map$SampleName
#   return(mat)
#}

path_join <- function(path1, path2) {
  return(paste(path1, path2, sep = ""))
}

build_group_info <- function(info, map) {
  # if(group){
  #    group_info = data.frame(data=info, group=info)
  # } else {
  group_info <- data.frame(data = I(info))
  group_info <- merge(group_info, map,
    by.x = "data",
    by.y = "SampleName", sort = F
  )
  # }
  return(group_info)
}


# clustering
RowColCluster <- function(taxadf, clustering_method, clustering_distance) {
  # row clustering
  if (clustering_method != "none") {
    if (nrow(taxadf) < 2) {
      rowPhylo <- list(messages = "The number of taxon is too low to cluster, please be sure that rowClusterMethod is 'none'!!!")
    } else {
      rowTree <- pheatmap:::cluster_mat(taxadf,
        distance = clustering_distance,
        method = clustering_method
      )
      taxadf <- taxadf[rowTree$order, , drop = FALSE]
      rowPhylo <- as.phylo(rowTree)
    }
  } else {
    rowPhylo <- NULL
  }
  # col clustering
  if (clustering_method != "none") {
    if (ncol(taxadf) < 2) {
      colPhylo <- list(messages = "The number of sample is too low to cluster, please be sure that rowClusterMethod is 'none'!!!")
    } else {
      colTree <- pheatmap:::cluster_mat(t(taxadf),
        distance = clustering_distance,
        method = clustering_method
      )
      taxadf <- taxadf[, colTree$order, drop = FALSE]
      colPhylo <- as.phylo(colTree)
    }
  } else {
    colPhylo <- NULL
  }
  return(list(taxadf = taxadf, rowPhylo = rowPhylo, colPhylo = colPhylo))
}


# tree to json
Phylo2Json <- function(phylo, method, dimension) {
  if (method == "none" || dimension < 2) {
    phyloJson <- toJSON(phylo)
    return(phyloJson)
  } else {
    ftldf <- get.data.frame(ape::as.igraph.phylo(phylo))
    ftldf$length <- phylo$edge.length
    # step1: make all sublist with 2 layers
    fromNodes <- unique(ftldf$from)
    allList <- list()
    for (i in fromNodes) {
      if (i == "Node1") {
        allList[[i]] <- list(name = "Node1", length = 0, children = list())
      } else {
        allList[[i]] <- list(
          name = i,
          length = ftldf[match(i, ftldf$to), 3], children = list()
        )
      }
      pos <- which(ftldf$from %in% i)
      for (j in 1:length(pos)) {
        if (ftldf[pos[j], 2] %in% ftldf$from) {
          allList[[i]]$children[[j]] <- list(
            name = ftldf[pos[j], 2],
            length = as.numeric(ftldf[pos[j], 3]),
            children = list()
          )
        } else {
          allList[[i]]$children[[j]] <- list(
            name = ftldf[pos[j], 2],
            length = as.numeric(ftldf[pos[j], 3])
          )
        }
      }
    }
    # step2: place the lower sublist in the higher sublist
    toNodes <- intersect(ftldf$from, ftldf$to) %>% rev()
    for (i in toNodes) {
      fnodes <- ftldf$from[which(ftldf$to %in% i)]
      for (j in 1:2) {
        if (i == allList[[fnodes]]$children[[j]]$name) {
          allList[[fnodes]]$children[[j]] <- allList[[i]]
        }
      }
    }
    # compensate for the second layer without subnodes
    allList <- allList$Node1
    for (i in 1:2) {
      if (is.null(allList$children[[i]]$children)) {
        allList$children[[i]]$children <- list(allList$children[[i]])
      }
    }
    # list to Json
    phyloJson <- toJSON(allList)
    return(phyloJson)
  }
}


# taxadf to json
data2json <- function(taxadf) {
  taxadf_ori <- taxadf
  taxadf <- t(taxadf) %>%
    data.frame(sample = row.names(.), ., stringsAsFactors = FALSE) %>%
    rbind(c("sample", rownames(taxadf_ori)), .)
  colnames(taxadf) <- NULL
  rownames(taxadf) <- NULL
  # toJSON(taxadf)
  taxaJson <- toJSON(taxadf)
  return(taxaJson)
}



# build group info and convert it to string of json format
GrpInfo2json <- function(taxadf, mapdf) {
  grp <- mapdf[match(colnames(taxadf), mapdf$NewSampleName), 3]

  grpInfo <- data.frame(data = colnames(taxadf), group = grp, stringsAsFactors = F)
  # construct new index of group in json file
  grpIndex <- c(1)
  for (i in 1:length(grp)) {
    ifelse(grp[i + 1] == grp[i],
      grpIndex[i + 1] <- grpIndex[i],
      grpIndex[i + 1] <- grpIndex[i] + 1
    )
  }
  grpInfo$index <- grpIndex
  tmpList <- list()
  for (i in 1:max(grpInfo$index)) {
    tmpList[[i]] <- list(
      group = unique(grpInfo[which(grpInfo$index %in% i), "group"]),
      data = c(grpInfo[which(grpInfo$index %in% i), "data"], "xxxx")
    )
  }
  grpJson <- toJSON(tmpList) %>% gsub(',\"xxxx\"', "", .)
  return(grpJson)
}


group2json <- function(mat, map_d) {
  grpJson <- list()
  grpJson0 <- list()
  index <- 0
  group <- build_group_info(colnames(mat), map_d[, 2:3])
  for (i in unique(group[, 2])) {
    index <- index + 1
    grpJson0[["group"]] <- i
    grpJson0[["data"]] <- group[which(group[, 2] == i), 1]
    grpJson[[index]] <- grpJson0
  }

  grpJson <- toJSON(grpJson)
  return(grpJson)
}



pheatmap_plot <- function(sample_cor, col_points, annotation_col, ann_colors) {
  min_value <- min(apply(sample_cor, 2, min))
  max_value <- max(apply(sample_cor, 2, max))
  # mean_value <- mean(apply(sample_cor, 2, mean))
  pheatmap(sample_cor, # 要绘制热图的矩阵
    color = colorRampPalette(col_points[1:3])(100), # 热图色块颜色是从蓝到红分为100个等级
    border_color = "white", # 热图中每个色块的边框颜色，NA表示无边框
    scale = "none", # 按行进行归一化，"column"表示按列，"none"表示不进行归一化
    cluster_rows = TRUE, # 是否对行进行聚类
    cluster_cols = TRUE, # 是否对列进行聚类
    legend = TRUE, # 是否显示图例
    display_numbers = TRUE, number_color = "black", number_format = "%.4f", # 是否显示数值
    legend_breaks = c(round(min_value, 4)+ 0.0001, round((max_value+min_value)/2, 4), round(max_value), 4), # 设置图例的断点
    annotation_col = annotation_col, annotation_colors = ann_colors, # 分组标签，分组颜色
    # legend_labels = c("low", "", "heigh"), #设置图例断点处的标签
    show_rownames = TRUE, # 是否显示行名
    show_colnames = TRUE, # 是否显示列名
    fontsize = fontsize, # 字体大小，可以通过fontsize_row、fontsize_col参数分别设置行列名的字体大小
    fontfamily = fontfamily, # 字体
    main = "Correlation Heatmap",
    angle_col = 45
  )

}


# build group info and convert it to string of json format
GrpInfo2json <- function(funcdf, mapdf) {
  grp <- mapdf[match(colnames(funcdf), mapdf[, 2]), "Default"]

  grpInfo <- data.frame(data = colnames(funcdf), group = grp, stringsAsFactors = F)
  # construct new index of group in json file
  grpIndex <- c(1)
  for (i in 1:length(grp)) {
    ifelse(grp[i + 1] == grp[i],
      grpIndex[i + 1] <- grpIndex[i],
      grpIndex[i + 1] <- grpIndex[i] + 1
    )
  }
  grpInfo$index <- grpIndex
  tmpList <- list()
  for (i in 1:max(grpInfo$index)) {
    tmpList[[i]] <- list(
      group = unique(grpInfo[which(grpInfo$index %in% i), "group"]),
      data = c(grpInfo[which(grpInfo$index %in% i), "data"], "xxxx")
    )
  }
  grpJson <- toJSON(tmpList) %>% gsub(',\"xxxx\"', "", .)
  return(grpJson)
}



main <- function(conf, offline = "null") {
  if (exists("conf")) {
    source(conf)

    map_d <- read.table(map, skip = 1, sep = "\t", check.names = F)
    primary_data <- read.table(input, header = T, sep = "\t", check.names = F)
  } else {
    primary_data <- read.table(header = T, sep = "\t", text = input, check.names = F)
    map_d <- read.table(skip = 1, sep = "\t", text = map, check.names = F)
  }


  if (!file.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  #primary_data <- process_mat(primary_data, map_d)

  # 过滤表达量为0的列
  primary_data <- primary_data[which(rowSums(primary_data[-1]) > 0), ]

  data <- scale_mat(primary_data, scale)
  data <- as.matrix(data)
  
  numercial_data <- data[,-1]
  coverted_data <- data.frame(apply(numercial_data,2,function(x){as.numeric(x)}))
  sample_cor <- cor(coverted_data, method = cor_method, use = "pairwise.complete.obs")

  mat <- round(sample_cor, 4)
  cor_data <- sample_cor
  # print(cor_data)
  # cor_data <- as.data.frame(sample = rownames(cor_data), cor_data)
  cor_data <- cbind("sample" = rownames(cor_data), cor_data)
  row.names(cor_data) <- NULL
  #write.table(cor_data, file = paste0(outdir, "/correlation.tsv"), row.names = FALSE, append = FALSE, quote = FALSE, sep = "\t")
  cor_data <- as.data.frame(cor_data)
  write.xlsx(cor_data, file = paste0(outdir, "/correlation.xlsx"))


  if (offline == "offline") {
    ## 加载颜色方案
    set_theme_module <- paste0(utils_root, "/set_Theme.R")
    source(set_theme_module)

    if (color_label %in% names(color.ls)) {
      col_label <- color.ls[[color_label]]
    } else if (color_label %in% names(color.continuous.ls)) {
      col_label <- color.continuous.ls[[color_label]]
    }


    if (color_name %in% names(color.ls)) {
      col_points <- color.ls[[color_name]]
    } else if (color_name %in% names(color.continuous.ls)) {
      col_points <- color.continuous.ls[[color_name]]
    }


    ## 加载字体方案
    AddFont(font_path)
    showtext_auto()

    annotation_col <- data.frame(
      Group = map_d[, 3]
    )
    rownames(annotation_col) <- map_d[, 2]


    # 标签颜色
    # ann_colors = list(
    # Group = c("white", "firebrick"),
    # Group = c(A = "#1B9E77", C = "#D95F02")
    # GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
    # )
    q <- 0
    x <- list()
    for (i in unique(map_d[, 3])) {
      q <- q + 1
      x[[i]] <- col_label[q]
    }
    ann_colors <- list(Group = unlist(x))


    # 绘制热图
    pdf(paste0(outdir, "/correlation.pdf"), width = 12, height = 9)
    pheatmap_plot(sample_cor, col_points, annotation_col, ann_colors)
    dev.off()

    png(paste0(outdir, "/correlation.png"), width = 1200, height = 900)
    pheatmap_plot(sample_cor, col_points, annotation_col, ann_colors)
    dev.off()

    svg(paste0(outdir, "/correlation.svg"), width = 12, height = 9)
    pheatmap_plot(sample_cor, col_points, annotation_col, ann_colors)
    dev.off()
  } else {
    taxaDfTreeList <- RowColCluster(mat, clustering_method, clustering_distance)

    taxaDf <- taxaDfTreeList$taxadf
    rowPhylo <- taxaDfTreeList$rowPhylo
    colPhylo <- taxaDfTreeList$colPhylo


    # taxadf to json
    taxaJson <- data2json(taxadf = taxaDf) %>% stringr::str_replace_all(., "\"([0-9.\\-]+)\"", "\\1")
    # row tree to json
    rowJson <- Phylo2Json(rowPhylo, method = clustering_method, dimension = nrow(taxaDf))
    # col tree to json
    colJson <- Phylo2Json(colPhylo, method = clustering_method, dimension = ncol(taxaDf))

    # group to json
    # grpJson <- group2json(mat, map_d)
    grpJson <- GrpInfo2json(funcdf = taxaDf, mapdf = map_d)
    grouporder <- toJSON(levels(factor(map_d[, 3])))

    finJson <- paste0(
      "{\n",
      '"group": ', grpJson, ",\n",
      '"matrix": ', taxaJson, ",\n",
      '"rowJSON": ', rowJson, ",\n",
      '"colJSON": ', colJson, ",\n",
      '"grouporder": ', grouporder,
      "}"
    )

    output <- paste0(outdir, "/data.json")
    cat(finJson, file = output, fill = FALSE, labels = NULL, append = FALSE)
  }
}

# running----
# args <- commandArgs(TRUE)
# if (length(args) == 1) {
#   print("run script with conf online.")
#   conf <- args[1]
#   main(conf, )
# } else if (length(args) == 2) {
#   print("run script with conf offline.")
#   conf <- args[1]
#   offline <- "offline"
#   main(conf, offline)
# } else {
#   print("run script without conf.")
# }

conf <- "config/4.6_correlation_analysis_default.conf"
source(conf)
offline <- "offline"
main(conf, offline)


