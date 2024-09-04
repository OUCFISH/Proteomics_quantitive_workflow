# -*- coding: utf-8 -*-
# @Time    : 2024/7/31 16:02
# @Author  : liangxing.sun
# @FileName: 5.2_KEGG_html.R
# @Software: Rstudio
# @Email: liangxing.sun@personalbio.cn
rm(list = ls())
pacman::p_load(dplyr,magick,base64enc,htmltools,stringr,openxlsx,tidyr,foreach,doParallel)

#' Title 创建元素列表，即只创建HTML中的area，href以及coords
#'
#' @param common_area_df 
#'
#' @return 返回HTML列表
#' @export
#'
#' @examples
create_common_tags <- function(common_area_df) {
  html_common_tags <- list() # 初始化一个空列表来存储HTML tags
  for (i in 1:nrow(common_area_df)) {
    row <- common_area_df[i, ] # 读取每行数据
    href <- paste0("http://www.genome.jp/", row[["V5"]]) # 构建href链接
    if(row[["V3"]] == "poly"){
      coords_temp = as.numeric(unlist(strsplit(row[["coord"]], ",")))
      split_vec <- split(coords_temp, rep(1:ceiling(length(coords_temp)/4), each = 4, length.out = length(coords_temp)))
      temp<-c()
      n=length(split_vec)
      for(n in length(split_vec):1){
        temp<-c(temp,tail(split_vec[[n]], 2), head(split_vec[[n]], 2))
      }
      left<-paste(unlist(split_vec),sep=",",collapse = ",")
      right<-paste(temp,sep=",",collapse = ",")
      coords<-paste(left,right,sep = ",")
    }else{
      coords = row[["coord"]]
    }
    html_common_tags[[i]] <- tags$area(shape = row[["V3"]], coords = coords, href = href) # 使用tags$area函数创造HTML标签并添加到列表
  }
  return(html_common_tags) # 返回HTML tags列表
}
## 构建 单色 area 部分

#' Title 构建 单色 area 映射到本次富集蛋白ID以及pvalue、FC
#'
#' @param single_area_df 单色数据映射
#' @param pathway_up_color 上调颜色，和get_table_color保持一致
#' @param pathway_down_color 下调颜色，和get_table_color保持一致
#' @param pathway_annot_color 注释颜色，和get_table_color保持一致
#'
#' @return 返回一个HTML列表
#' @export
#'
#' @examples
create_single_color_tags <- function(single_area_df,pathway_up_color="red",pathway_down_color="blue",pathway_annot_color="orange") {
  html_single_tags <- list() # 初始化一个空列表来存储HTML tags
  for (i in 1:nrow(single_area_df)) {
    row <- single_area_df[i, ] # 读取每行数据
    href <- paste0("http://www.genome.jp/", row[["V5"]]) # 构建href链接
    inner_li <- tags$li(row[["areaddinfo"]]) # 创建内层的<li>元素
    # 将内层<li>包裹在<ul>内
    inner_ul <- tags$ul(inner_li)
    if(row[["Color"]] == pathway_up_color){
      outer_li <- HTML(paste0("<li style=", "\\\"color: ", row[["Color"]], ";\\\">","Up_Regulation", inner_ul, "</li>"))
    }else if(row[["Color"]] == pathway_down_color){
      outer_li <- HTML(paste0("<li style=", "\\\"color: ", row[["Color"]], ";\\\">","Down_Regulation", inner_ul, "</li>"))
    }else{
      outer_li <- HTML(paste0("<li style=", "\\\"color: ", row[["Color"]], ";\\\">","Annotation", inner_ul, "</li>"))
    }
    onmouseover <- paste0(" onmouseover=", "\'javascript: showInfo(\"", tags$ul(outer_li), "\");'")
    onmouseover <- gsub("\n", "", onmouseover)
    coords<-c()
    if(row[["V3"]] == "poly"){
      coords_temp = as.numeric(unlist(strsplit(row[["coord"]], ",")))
      split_vec <- split(coords_temp, rep(1:ceiling(length(coords_temp)/4), each = 4, length.out = length(coords_temp)))
      temp<-c()
      n=length(split_vec)
      for(n in length(split_vec):1){
        temp<-c(temp,tail(split_vec[[n]], 2), head(split_vec[[n]], 2))
      }
      left<-paste(unlist(split_vec),sep=",",collapse = ",")
      right<-paste(temp,sep=",",collapse = ",")
      coords<-paste(left,right,sep = ",")
    }else{
      coords = row[["coord"]]
    }
    html_single_tags[[i]] <- HTML(paste0("<area shape='", row[["V3"]], "' coords='", coords, "' href='", href, "'", HTML(onmouseover), " />"))
  }
  return(html_single_tags) # 返回HTML tags列表
}
## 构建 双色 area 部分
#' Title构建 双色 area 映射到本次富集蛋白ID以及pvalue、FC
#'
#' @param double_area_df 双色数据映射
#' @param pathway_up_color 上调颜色，和get_table_color保持一致
#' @param pathway_down_color 下调颜色，和get_table_color保持一致
#' @param pathway_annot_color 注释颜色，和get_table_color保持一致
#'
#' @return 返回HTML列表
#' @export
#'
#' @examples
create_double_color_tags <- function(double_area_df,pathway_up_color="red",pathway_down_color="blue",pathway_annot_color="orange") {
  html_double_tags <- list() # 初始化一个空列表来存储HTML tags
  for (i in 1:nrow(double_area_df)) {
    row <- double_area_df[i, ] # 读取每行数据
    href <- paste0("http://www.genome.jp/", row[["V5"]]) # 构建href链接
    # 创建内层的<li>元素
    inner_li1 <- tags$li(row[["up_areaddinfo"]])
    inner_li2 <- tags$li(row[["down_areaddinfo"]])
    # 将内层<li>包裹在<ul>内
    inner_ul1 <- tags$ul(inner_li1)
    inner_ul2 <- tags$ul(inner_li2)
    # 创建外层的<li>元素，将之前创建的内层<ul>作为内容，并添加样式
    outer_li1 <- HTML(paste0("<li style=", "\\\"color: ", row[["up_color"]], ";\\\">","Up_Regulation " ,inner_ul1, "</li>"))
    outer_li2 <- HTML(paste0("<li style=", "\\\"color: ", row[["down_color"]], ";\\\">","Down_Regulation ", inner_ul2, "</li>"))
    onmouseover <- paste0(" onmouseover=", "\'javascript: showInfo(\"", tags$ul(list(outer_li1, outer_li2)), "\");'")
    onmouseover <- gsub("\n", "", onmouseover)
    
    html_double_tags[[i]] <- HTML(paste0("<area shape='", row[["V3"]], "' coords='", row[["coord"]], "' href='", href, "'", HTML(onmouseover), " />"))
  }
  return(html_double_tags) # 返回HTML tags列表
}

#-------------------------- 构建 底图方框颜色 函数 --------------------------#
## up and down color
#' Title png绘制双色框
#'
#' @param row_coords 三个元素，位点，上调颜色，下调颜色
#'
#' @return
#' @export 直接在图上画线
#'
#' @examples
up_down_color_fun <- function(row_coords) {
  coords <- row_coords[1]
  shape <- row_coords[2]
  up_color <- row_coords[3]
  down_color <- row_coords[4]
  coords_vector <- as.numeric(unlist(strsplit(coords, ",")))
  
  if (length(coords_vector) != 4) stop("Each coordinate set should contain exactly four numbers.")
  
  # 提取左上角（x1,y1）和右下角（x2,y2）坐标
  x1 <- coords_vector[1]
  y1 <- coords_vector[2]
  x2 <- coords_vector[3]
  y2 <- coords_vector[4]
  if(shape == "rect"){
    lines(x = c(x1, x1, x2, x2), y = c((y2 - y1) / 2 + y1, y1, y1, (y2 - y1) / 2 + y1), col = up_color, lwd = 4) ## up 的颜色
    lines(x = c(x1, x1, x2, x2), y = c((y2 - y1) / 2 + y1, y2, y2, (y2 - y1) / 2 + y1), col = down_color, lwd = 4) ## down 的颜色
  }else if(shape == "poly"){
    lines(x = c((x2 - x1) / 2 + x1, x1, x1, (x2 - x1) / 2 + x1), y = c(y1,y1,y2,y2), col = up_color, lwd = 4) ## up 的颜色
    lines(x = c((x2 - x1) / 2 + x1, x2, x2, (x2 - x1) / 2 + x1), y = c(y1,y1,y2,y2), col = down_color, lwd = 4) ## down 的颜色
  }else if(shape == "circle"){
    draw.circle_up <- function(center, radius, col) {
      angles <- seq(0, pi, length.out = 100)
      x <- center[1] + radius * cos(angles)
      y <- center[2] + radius * sin(angles)
      polygon(x, y, col = col, border = NA)
    }
    
    draw.circle_down <- function(center, radius, col) {
      angles <- seq( -pi,0, length.out = 100)
      x <- center[1] + radius * cos(angles)
      y <- center[2] + radius * sin(angles)
      polygon(x, y, col = col, border = NA)
    }
    
    center <- c(x, y)
    radius <- 1.5
    
    # 分别绘制圆的上半部分和下半部分
    draw.circle_down(center, radius, down_color)
    draw.circle_up(center, radius, up_color)
    
  }
  
  
}
## single_rect_fun
#' Title png绘制单色框
#'
#' @param row_coords 两元素，位点，颜色
#'
#' @return
#' @export 直接在图上绘制矩形
#'
#' @examples
single_rect_fun <- function(row_coords) {
  coords <- row_coords[1]
  shape <- row_coords[2]
  color <- row_coords[3]
  coords_vector <- as.numeric(unlist(strsplit(coords, ",")))
  
  if (length(coords_vector) != 4) stop("Each coordinate set should contain exactly four numbers.")
  
  # 提取左上角（x1,y1）和右下角（x2,y2）坐标
  x1 <- coords_vector[1]
  y1 <- coords_vector[2]
  x2 <- coords_vector[3]
  y2 <- coords_vector[4]
  
  if(shape == "rect"){
    # 用指定颜色绘制矩形
    graphics::rect(
      xleft = x1, ybottom = y1, xright = x2, ytop = y2,
      border = color, col = NA, lwd = 4
    ) # col=NA 是为了避免填充
  }else if(shape == "poly"){
    lines(x=c(x1,x2),y=c(y1,y2),col = color, lwd = 4)
  }else if(shape == "circle"){
    draw.circle_up <- function(center, radius, col) {
      angles <- seq(0, 2*pi, length.out = 100)
      x <- center[1] + radius * cos(angles)
      y <- center[2] + radius * sin(angles)
      polygon(x, y, col = col, border = NA)
    }
    center <- c(x, y)
    radius <- 1.5
    
    # 分别绘制圆的上半部分和下半部分
    draw.circle(center, radius, color)
  }
  
}
#-------------------------- 构建 颜色和html所需表格 函数 --------------------------#
## 获取 KO 号颜色信息
#' Title
#'
#' @param enrichment_df 富集表
#' @param diff_df 差异表
#' @param gene_type 类型，DEP，up,down,annot
#' @param pathway_up_color 上调颜色："red"
#' @param pathway_down_color 下调颜色："blue"
#' @param pathway_annot_color 注释颜色："orange"
#' @param coord_df 位点表：all_data
#' @param prefix_pathwayid KEGG前缀：rno
#'
#' @return list,list[1]:ID,coords,color,annotation,重复个数；list[2]所相关的信息ID、p、FC、coords等等
#'
#' @examples get_table_color(enrichment_df = enrichment_df,diff_df = diff_df ,gene_type = "DEP",
#'             pathway_up_color = pathway_up_color,
#'             pathway_down_color = pathway_down_color, 
#'             coord_df = all_data, prefix_pathwayid=organism_code)
#'             
get_table_color <- function(enrichment_df,diff_df, gene_type, pathway_up_color="red",
                            pathway_down_color="blue",pathway_annot_color="orange", 
                            coord_df, prefix_pathwayid) {
  if(gene_type != "annot"){
    if(all(enrichment_df$Up_KEGG == "-")){
      gene_type = "down"
    }else if(all(enrichment_df$Down_KEGG == "-")){
      gene_type = "up"
    }else{
      gene_type = "DEP"
    }
  }
  ## 获取 DEP、UP、DOWN 的info
  if (gene_type == "DEP") {
    diff_df<-diff_df%>%dplyr::select(Protein.Accessions,dplyr::contains(c("Log2fc","pvalue","Gene.Name")))
    colnames(diff_df)<-c("Accession","Log2fc","pvalue","Gene.Name")
    kegg_up <- enrichment_df %>%
      dplyr::select(ID, Up_KEGG) %>%
      dplyr::filter(Up_KEGG != "-") %>%
      separate_rows(Up_KEGG, sep = ";") %>%
      separate(Up_KEGG, into = c("Accession", "KEGG"), sep = "\\(")%>%
      left_join(.,diff_df,by=join_by(Accession==Accession))
    # kegg_up$Gene <- paste0(kegg_up$Gene, ")")
    kegg_up$KEGG <- gsub("\\)", "", kegg_up$KEGG)
    kegg_up$Color <- rep(all_of(pathway_up_color))
    kegg_up$Accession<-paste0(kegg_up$Accession,":(Gene Name:",kegg_up$Gene.Name," ,log2fc: ",round(kegg_up$Log2fc,2)," ,pvalue: ",format(kegg_up$pvalue, scientific = T),")") 
    kegg_down <- enrichment_df %>%
      dplyr::select(ID, Down_KEGG) %>%
      dplyr::filter(Down_KEGG != "-") %>%
      separate_rows(Down_KEGG, sep = ";") %>%
      separate(Down_KEGG, into = c("Accession", "KEGG"), sep = "\\(")%>%
      left_join(.,diff_df,by=join_by(Accession==Accession))
    # kegg_down$Gene <- paste0(kegg_down$Gene, ")")
    kegg_down$KEGG <- gsub("\\)", "", kegg_down$KEGG)
    kegg_down$Color <- rep(all_of(pathway_down_color))
    kegg_down$Accession<-paste0(kegg_down$Accession,":(Gene Name:",kegg_down$Gene.Name," ,log2fc: ",round(kegg_down$Log2fc,2)," ,pvalue: ",format(kegg_down$pvalue, scientific = T),")")
    kegg_color_df <- rbind(kegg_up, kegg_down) ## 合并 up 和 down 的颜色
    
    #KEGG_annotation <- read.csv(pathway2gene)
    #Gene_KO <- KEGG_annotation[c("Symbol","KO.ID")]
    #matched_indices <- match(kegg_color_df$Gene, Gene_KO$Symbol)
    #kegg_color_df$KO <- Gene_KO$KO.ID[matched_indices]
  } else if (gene_type == "up") {
    diff_df<-diff_df%>%dplyr::select(Protein.Accessions,dplyr::contains(c("Log2fc","pvalue","Gene.Name")))
    colnames(diff_df)<-c("Accession","Log2fc","pvalue","Gene.Name")
    kegg_color_df <- enrichment_df %>%
      dplyr::select(ID, Up_KEGG) %>%
      dplyr::filter(Up_gene != "-") %>%
      separate_rows(Up_gene, sep = ";") %>%
      separate(Up_gene, into = c("Accession", "KEGG"), sep = "\\(")%>%
      left_join(.,diff_df,by=join_by(Accession==Accession))
    #kegg_color_df$Accession <- paste0(kegg_color_df$Gene, ")")
    kegg_color_df$Accession<-paste0(kegg_color_df$Accession,":(Gene Name:",kegg_color_df$Gene.Name," ,log2fc: ",round(kegg_color_df$Log2fc,2)," ,pvalue: ",format(kegg_color_df$pvalue,scientific = T),")")
    kegg_color_df$KEGG <- gsub("\\)", "", kegg_color_df$KEGG)
    kegg_color_df$Color <- rep(all_of(pathway_up_color))
  } else if(gene_type=="down") {
    diff_df<-diff_df%>%dplyr::select(Protein.Accessions,dplyr::contains(c("Log2fc","pvalue","Gene.Name")))
    colnames(diff_df)<-c("Accession","Log2fc","pvalue","Gene.Name")
    kegg_color_df <- enrichment_df %>%
      dplyr::select(ID, Down_KEGG) %>%
      dplyr::filter(Down_protein != "-") %>%
      separate_rows(Down_protein, sep = ";") %>%
      separate(Down_protein, into = c("Accession", "KEGG"), sep = "\\(")%>%
      left_join(.,diff_df,by=join_by(Accession==Accession))
    #kegg_color_df$Accession <- paste0(kegg_color_df$Gene, ")")
    kegg_color_df$Color <- gsub("\\)", "", kegg_color_df$KO)
    kegg_color_df$Accession<-paste0(kegg_color_df$Accession,":(Gene Name:",kegg_color_df$Gene.Name," ,log2fc: ",round(kegg_color_df$Log2fc,2)," ,pvalue: ",format(kegg_color_df$pvalue,scientific = T),")")
    kegg_color_df$Color <- rep(all_of(pathway_down_color))
  }else if(gene_type=="annot"){
    diff_df<-diff_df%>%dplyr::select(Protein.Accessions,dplyr::contains(c("Log2fc","pvalue","Gene.Name")))
    colnames(diff_df)<-c("Accession","Gene.Name")
    kegg_color_df <- enrichment_df %>%
      dplyr::select(ID, Protein) %>%
      dplyr::filter(Protein != "-") %>%
      separate_rows(Protein, sep = ";") %>%
      separate(Protein, into = c("Accession", "KEGG"), sep = "\\(")%>%
      left_join(.,diff_df,by=join_by(Accession==Accession))
    #kegg_color_df$Gene <- paste0(kegg_color_df$Gene, ")")
    kegg_color_df$KEGG <- gsub("\\)", "", kegg_color_df$KEGG)
    kegg_color_df$Accession<-paste0(kegg_color_df$Accession,":(Gene Name:",kegg_color_df$Gene.Name," )")
    kegg_color_df$Color <- rep(all_of(pathway_annot_color))
  }
  ## ko 替换为 map, 使得ID 的类型统一
  # kegg_color_df$PathwayID <- gsub("[^0-9]+", "map", kegg_color_df$PathwayID)
  #kegg_color_df$PathwayID <- gsub("[^0-9]+", prefix_pathwayid, kegg_color_df$ID)
  kegg_color_df$uniqueid <- paste0(kegg_color_df$ID, "_", kegg_color_df$KEGG)
  
  ## 获取 KO 号坐标信息
  #coord_df <- fread(coord_info, header = F, sep = "\t")
  colnames(coord_df) <- c("PathwayID", "KEGG", "shape", "coord","href")
  # coord_df$PathwayID <- gsub("[^0-9]+", "map", coord_df$PathwayID)
  #coord_df$PathwayID <- gsub("[^0-9]+", prefix_pathwayid, coord_df$PathwayID)
  if(prefix_pathwayid == "ko"){
    coord_df$uniqueid <- paste0(coord_df$PathwayID, "_", coord_df$KEGG)%>%gsub(" ","",x=.)
  }else{
    coord_df$uniqueid <- paste0(coord_df$PathwayID, "_",prefix_pathwayid,":", coord_df$KEGG)%>%gsub(" ","",x=.)
  }
  ## 将坐标信息 和 颜色信息合并为一个表格
  kegg_color_df <- inner_join(kegg_color_df, coord_df[, c("uniqueid", "coord","shape")], by = "uniqueid")
  unique_color_df <- kegg_color_df %>%
    group_by(coord) %>%
    summarise(unique_color = n_distinct(Color))
  kegg_color_df <- inner_join(kegg_color_df, unique_color_df, by = "coord")
  kegg_color_df$unique_areainfo <- paste0(kegg_color_df$Color, "_", kegg_color_df$coord)
  ## 根据unique_areainfo 进行 group_by 拼接 areainfo 的内容
  areainfo_tmpdf <- kegg_color_df %>%
    group_by(unique_areainfo) %>%
    summarise(paste0(Accession, collapse = "<br>"))
  colnames(areainfo_tmpdf) <- c("unique_areainfo", "addinfo")
  areainfo_df <- inner_join(kegg_color_df, areainfo_tmpdf, by = "unique_areainfo") %>%
    mutate(areaddinfo = paste0(KEGG, ":<br> ", addinfo)) %>%
    select(ID, coord, Color, areaddinfo, unique_color)
  
  areainfo_df <- na.omit(areainfo_df)
  combined_info <- paste(areainfo_df$ID, areainfo_df$coord, areainfo_df$Color, sep = "_")
  # 然后用duplicated()检查这个组合向量的重复项
  areainfo_df <- areainfo_df[!duplicated(combined_info), ]
  return(list(areainfo_df, kegg_color_df))
}

#通过map.html可以选取文件夹中的html，html_content是一部分header和javascript
html_content <- "<html>
<head>
<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">
<title>CK-VS-TSG</title>
<style type=\"text/css\">
body {background-color: #fff;}
table {background-color: #000; border-collapse: collapse; border: solid #000 1px; margin: 0 0 50px 0;}
tr {background-color: #fff;}
th, td {border: solid #000 1px;}
</style>
<script type=\"text/javascript\">
<!--
function reSize2() {
	try {
		parent.document.getElementsByTagName(\"iframe\")[0].style.height = document.body.scrollHeight + 10;
		parent.parent.document.getElementsByTagName(\"iframe\")[0].style.height = parent.document.body.scrollHeight;
	} catch(e) {}
}

var preRow = null;
var preColor = null;
function colorRow(trObj) {
	if (preRow != null) {
		preRow.style.backgroundColor = preColor;
	}
	preRow = trObj;
	preColor = trObj.style.backgroundColor;
	trObj.style.backgroundColor = \"#ff0\";
}

function diffColor(tables) {
	var color = [\"#fff\", \"#efefef\"];
	for (var i = 0; i < tables.length; i++) {
		var trObj = tables[i].getElementsByTagName(\"tr\");
		for (var j = 1; j < trObj.length; j++) {
			trObj[j].style.backgroundColor = color[j % color.length];
		}
	}
}

function showPer(tableObj) {
	trObj = tableObj.getElementsByTagName(\"tr\");
	if (trObj.length < 2) {
		return;
	}
	var sum1 = trObj[0].cells[2].innerHTML.replace(/^.*\\((\\d+)\\).*$/, \"$1\");
	var sum2 = 0;
	if (trObj[0].cells.length > 4) {
		sum2 = trObj[0].cells[3].innerHTML.replace(/^.*\\((\\d+)\\).*$/, \"$1\");
	}
	trObj[0].cells[2].innerHTML = \"DEGs with pathway annotation (\" + sum1 + \")\";
	if (trObj[0].cells.length > 4) {
		trObj[0].cells[3].innerHTML = \"All genes with pathway annotation (\" + sum2 + \")\";
	}
	for (var i = 1; i < trObj.length; i++) {
		trObj[i].cells[2].innerHTML += \" (\" + (Math.round(trObj[i].cells[2].innerHTML * 10000 / sum1) / 100) + \"%)\";
		if (trObj[0].cells.length > 4) {
			trObj[i].cells[3].innerHTML += \" (\" + (Math.round(trObj[i].cells[3].innerHTML * 10000 / sum2) / 100) + \"%)\";
		}
	}
}

window.onload = function() {
	setTimeout(\"reSize2()\", 1);
}
//-->
</script>
</head><body>"

#' Title创建引导网页
#'
#' @param enrich_data 富集表
#' @param contrast 比较组
#' @param html_content 引导网页的head 与javascript
#'
#' @return 一个以富集表展示链接到下层目录kegg_plot网页的引导页
#' @export
#'
#' @examples
kegg_index_html <- function(enrich_data, contrast,html_content) {
  html <- paste0("<table><caption style='font-weight: 900;'>",contrast,"</caption>")
  html <- paste0(html_content,html)
  
  # 打开并读取文件
  # enrich_data <- fread(enrich_file,header=T,sep="\t")
  if("URL" %in% colnames(enrich_data)) enrich_data[,"URL"] <- NULL
  html <- paste0(html, "<tr><th>", paste0(colnames(enrich_data), collapse = "</th><th>"), "</th></tr>")
  total_html <- function(lines,contrast) {
    mapid <- lines[1]
    lines[1] <- paste0("<a href='", contrast, "_KEGG_Pathway_plot/", mapid, ".html' title='click to view map'>", lines[1])
    pre_html <- paste0(lines, collapse = "</td><td>")  
    return(pre_html)
  }
  pre_html <- apply(enrich_data, 1, function(x) { total_html(x, contrast) })
  pre_html <- paste0(pre_html,collapse = "</td></tr><tr><td>")
  
  html <- paste0(html, "<tr><td>", pre_html, "</table></body></html>")
  return(html)
}

#' Title 主函数获得KEGG网页与png
#'
#' @param kegg_color_df 从get_table_color函数自己整理的表格[[2]]
#' @param areainfo_df 从get_table_color函数整理的表格[[1]]
#' @param map_path png路径
#' @param area_df org相关信息，
#' @param suboutdir 输出路径
#' @param contrast 比较组
#' @param prefix_pathwayid KEGG物种缩写
#' @param enrichment_df 富集表格为了做一个索引网页
#' @param ncore 并行核数默认4
#'
#' @return
#' @export 创建文件夹其中有html与png图片
#'
#' @examples
get_KEGG_html_file <- function(kegg_color_df, areainfo_df, map_path, area_df, suboutdir,enrichment_df, contrast, prefix_pathwayid,ncore=4) {
  cl<-makeCluster(ncore)
  registerDoParallel(cl)
  area_df$PathwayID <- gsub("[^0-9]+", prefix_pathwayid, area_df$V1)
  foreach (single_pathway = unique(kegg_color_df$ID),
           .errorhandling = "stop",
           .packages =c("dplyr","magick","base64enc","htmltools","stringr","openxlsx","tidyr"),
           .export = c("conf","create_common_tags","create_double_color_tags",
                       "create_single_color_tags","get_KEGG_html_file",
                       "get_table_color","html_content","kegg_index_html",
                       "main","pathway_annot_color","pathway_down_color",
                       "pathway_up_color","single_rect_fun","up_down_color_fun",
                       "organism_code","coord_df","map_path","all_data")
  ) %dopar% {
    except_id = c("01100","01110","01120","01200","01210","01212","01230","01232","01250","01240","01220")
    except_id = paste0(prefix_pathwayid,except_id)
    if(single_pathway %in% except_id){
      next
    }
    #print(single_pathway)
    map_df <- kegg_color_df[kegg_color_df$ID == single_pathway, ] ## 取出对应的坐标和颜色清单，区分单个颜色和多个颜色 --> 画图
    html_df <- areainfo_df[areainfo_df$ID == single_pathway, ] ## 取出对应的坐标 和 附属信息 --> 构建 html 内容
    
    ## 获取pathway 对应表格 --> color 和 coord
    if (2 %in% unique(map_df$unique_color)) {
      ## 画图
      single_color_df <- map_df[map_df$unique_color == 1, c("coord","shape", "Color")] ## 单色
      single_color_df <- single_color_df[!duplicated(single_color_df$coord), ]
      
      ## 双色
      ## 构建 up color 和 down color 根据坐标合并颜色
      double_color_df <- map_df[map_df$unique_color == 2, c("coord","shape", "Color")] %>% mutate(id = paste0(coord, "_", Color))
      double_color_df <- double_color_df[!duplicated(double_color_df$id), ] %>% select(-id)
      up_color_df <- double_color_df[double_color_df$Color == all_of(pathway_up_color), ]
      down_color_df <- double_color_df[double_color_df$Color == all_of(pathway_down_color), ]
      double_color_df <- inner_join(up_color_df, down_color_df, by = c("coord","shape"))
      colnames(double_color_df) <- c("coord","shape", "up_color", "down_color")
      
      ## html
      single_html_df <- html_df[html_df$unique_color == 1, ] ## 单色
      ## 双色
      up_html_df <- html_df[html_df$unique_color == 2 & html_df$Color == all_of(pathway_up_color), ]
      colnames(up_html_df)[3:4] <- c("up_color", "up_areaddinfo")
      down_html_df <- html_df[html_df$unique_color == 2 & html_df$Color == all_of(pathway_down_color), ]
      colnames(down_html_df)[3:4] <- c("down_color", "down_areaddinfo")
      double_html_df <- inner_join(up_html_df[, c("coord", "up_color", "up_areaddinfo")], down_html_df[, c("coord", "down_color", "down_areaddinfo")], by = "coord")
      ## 确保去除所有NA值,将为空的变量删除，方便后续判断
      single_color_df <- na.omit(single_color_df)
      single_html_df <- na.omit(single_html_df)
      double_color_df <- na.omit(double_color_df)
      double_html_df <- na.omit(double_html_df)
      rm(up_html_df)
      rm(down_html_df)
      rm(up_color_df)
      rm(down_color_df)
      if (nrow(single_color_df) == 0) rm(single_color_df)
      if (nrow(single_html_df) == 0) rm(single_html_df)
      if (nrow(double_color_df) == 0) rm(double_color_df)
      if (nrow(double_html_df) == 0) rm(double_html_df)
    } else {
      ## color
      single_color_df <- map_df[map_df$unique_color == 1, c("coord","shape", "Color")]
      single_color_df <- single_color_df[!duplicated(single_color_df$coord), ]
      single_color_df <- na.omit(single_color_df)
      if (nrow(single_color_df) == 0) rm(single_color_df)
      ## html
      single_html_df <- html_df[html_df$unique_color == 1, ]
      single_html_df <- na.omit(single_html_df)
      if (nrow(single_html_df) == 0) rm(single_html_df)
    }
    
    ## 获取底图 --> 未做任何修饰的原图
    image <- image_read(paste0(map_path, "/", single_pathway, ".png")) # 读取 PNG 图片文件
    image <- image_draw(image) # 在绘制图像之前开始图像绘制流程
    
    if (!exists("double_color_df") & !exists("single_color_df")) next ## 若无需要标记的颜色则跳过当前的pathway
    
    ## 根据坐标数据画图
    if (exists("double_color_df")) {
      apply(double_color_df, 1, up_down_color_fun) ## 两个颜色 同时出现一个框里
      rm(double_color_df)
    }
    if (exists("single_color_df")) {
      apply(single_color_df, 1, single_rect_fun) ## 单色框
      rm(single_color_df)
    }
    
    # 完成图像绘制
    dev.off()
    
    ## 获取图片的 base64 编码,将图像保存到一个临时文件
    temp_file <- tempfile(fileext = ".png") # 创建临时文件
    image_write(image, path = temp_file, format = "png") # 将 magick-image 写入临时文件
    base64_string <- base64encode(temp_file) # 读取临时文件内容并转为 base64 编码
    unlink(temp_file) # 清理临时文件
    
    ## 构建 html 的输入文件
    #area_df <- fread(area_info, header = F, sep = "\t")
    # area_df$PathwayID <- gsub("[^0-9]+", "map", area_df$V1)
    
    pathway_area_df <- area_df[area_df$PathwayID == single_pathway, ]
    colnames(pathway_area_df)[4] <- "coord"
    pathway_area_df$V2 <- NULL
    if(any(duplicated(pathway_area_df$coord))){
      pathway_area_df <- pathway_area_df[!duplicated(pathway_area_df$coord),]
    }
    
    
    if (exists("double_html_df") & exists("single_html_df")) {
      common_area_df <- pathway_area_df %>%
        anti_join(single_html_df, by = "coord") %>%
        anti_join(double_html_df, by = "coord")
      double_area_df <- inner_join(double_html_df, pathway_area_df[, c("coord", "V3", "V5")], by = "coord")
      double_area_df <- na.omit(double_area_df)
      single_area_df <- inner_join(single_html_df, pathway_area_df[, c("coord", "V3", "V5")], by = "coord")
      single_area_df <- na.omit(single_area_df)
      if (nrow(double_area_df) > 0) html_double_tags <- create_double_color_tags(double_area_df)
      if (nrow(single_area_df) > 0) html_single_tags <- create_single_color_tags(single_area_df)
      
      rm(double_html_df)
      rm(single_html_df)
      rm(double_area_df)
      rm(single_area_df)
    } else if (!exists("double_html_df") & exists("single_html_df")) {
      common_area_df <- pathway_area_df %>% anti_join(single_html_df, by = "coord")
      single_area_df <- inner_join(single_html_df, pathway_area_df[, c("coord", "V3", "V5")], by = "coord")
      single_area_df <- na.omit(single_area_df)
      if (nrow(single_area_df) > 0) html_single_tags <- create_single_color_tags(single_area_df)
      rm(single_html_df)
      rm(single_area_df)
    } else if (exists("double_html_df") & !exists("single_html_df")) {
      common_area_df <- pathway_area_df %>% anti_join(double_html_df, by = "coord")
      double_area_df <- inner_join(double_html_df, pathway_area_df[, c("coord", "V3", "V5")], by = "coord")
      double_area_df <- na.omit(double_area_df)
      if (nrow(double_area_df) > 0) html_double_tags <- create_double_color_tags(double_area_df)
      rm(double_html_df)
      rm(double_area_df)
    }
    
    html_common_tags <- create_common_tags(common_area_df) ## 构建单个 html 文件
    rm(common_area_df)
    # 构造<head>内容
    # title_id <- gsub("map", "ko", single_pathway)
    title_id <- single_pathway
    head_tags <- tags$header(list(
      tags$meta(http_equiv = "content-type", content = "text/html; charset=utf-8"),
      tags$title(title_id),
      tags$style(HTML("<\narea {cursor: pointer;}\n")),
      tags$script(HTML("\nfunction showInfo(info) {\
            obj = document.getElementById(\"result\");\
            obj.innerHTML = \"<div style='cursor: pointer; position: absolute; right: 5px; color: #000;' onclick='javascript: document.getElementById(\\\"result\\\").style.display = \\\"none\\\";' title='close'>X</div>\" + info;\
            obj.style.top = document.body.scrollTop;\
            obj.style.left = document.body.scrollLeft;\
            obj.style.display = \"\";\
        }\n"))
    ))
    # 构造<body>内容
    if (exists("html_single_tags") & exists("html_double_tags")) {
      body_tags <- tags$body(
        tags$map(name = title_id, html_common_tags, html_single_tags, html_double_tags),
        HTML(paste0("<img src='data:image/png;base64,", base64_string, "' usemap='#", title_id, "' />")),
        HTML(paste0("<div id=", "\'result\'", " style=", "\'position: absolute; width: 50%; border: 1px solid #000; background-color: #fff; filter: alpha(opacity=95); opacity: 0.95; font-size: 12px; padding-right: 20px; display: none;\'", " onmouseover=", "\"javascript: this.style.filter = 'alpha(opacity=100)'; this.style.opacity = 1;\"", " onmouseout=", "\"javascript: this.style.filter = 'alpha(opacity=95)'; this.style.opacity = 0.95;\"></div>"))
      )
      rm(html_double_tags)
      rm(html_single_tags)
    } else if (exists("html_single_tags") & !exists("html_double_tags")) {
      body_tags <- tags$body(
        tags$map(name = title_id, html_common_tags, html_single_tags),
        HTML(paste0("<img src='data:image/png;base64,", base64_string, "' usemap='#", title_id, "' />")),
        HTML(paste0("<div id=", "\'result\'", " style=", "\'position: absolute; width: 50%; border: 1px solid #000; background-color: #fff; filter: alpha(opacity=95); opacity: 0.95; font-size: 12px; padding-right: 20px; display: none;\'", " onmouseover=", "\"javascript: this.style.filter = 'alpha(opacity=100)'; this.style.opacity = 1;\"", " onmouseout=", "\"javascript: this.style.filter = 'alpha(opacity=95)'; this.style.opacity = 0.95;\"></div>"))
      )
      rm(html_single_tags)
    } else if (!exists("html_single_tags") & exists("html_double_tags")) {
      body_tags <- tags$body(
        tags$map(name = title_id, html_common_tags, html_double_tags),
        HTML(paste0("<img src='data:image/png;base64,", base64_string, "' usemap='#", title_id, "' />")),
        HTML(paste0("<div id=", "\'result\'", " style=", "\'position: absolute; width: 50%; border: 1px solid #000; background-color: #fff; filter: alpha(opacity=95); opacity: 0.95; font-size: 12px; padding-right: 20px; display: none;\'", " onmouseover=", "\"javascript: this.style.filter = 'alpha(opacity=100)'; this.style.opacity = 1;\"", " onmouseout=", "\"javascript: this.style.filter = 'alpha(opacity=95)'; this.style.opacity = 0.95;\"></div>"))
      )
      rm(html_double_tags)
    }
    
    html_page <- HTML(as.character(tags$html(list(head_tags, body_tags))))
    rm(map_df)
    rm(html_df)
    rm(html_common_tags)
    rm(body_tags)
    rm(head_tags)
    #dir.create(paste0(suboutdir, "/", contrast, "_KEGG_Pathway_map"), showWarnings = F, recursive = T, mode = "0755")
    #write.table(html_page, file = paste0(suboutdir, "/", contrast, "_KEGG_Pathway_map/", single_pathway, ".html"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    # image_write(image, path = paste0(outdir, "/", contrast, "_KEGG_Pathway_map/", single_pathway, ".png"), format = "png") # 保存 image 为 png 
    
    # 保存为PNG格式
    dir.create(paste0(suboutdir, "/",contrast, "_KEGG_Pathway_plot"), showWarnings = F, recursive = T, mode = "0755")
    write.table(html_page, file = paste0(suboutdir, "/",contrast, "_KEGG_Pathway_plot/", single_pathway, ".html"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    image_write(image, path = paste0(suboutdir, "/",contrast, "_KEGG_Pathway_plot/", single_pathway, ".png"), format = "png") # 保存 image 为 png 
    
  }
  stopImplicitCluster()
  stopCluster(cl)
  html_index_all<-kegg_index_html(enrich_data = enrichment_df,contrast =contrast,html_content =  html_content)
  write.table(html_index_all, file = paste0(suboutdir, "/",contrast, ".KEGG_enrichment_map.html"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#' Title 创建所有HTML、png和index引导网页
#'
#' @param diff_df 差异表
#' @param enrichment_df 富集表
#' @param pathway_up_color 上调颜色
#' @param pathway_down_color 下调颜色
#' @param coord_df 位点表格信息
#' @param organism_code KEGG物种缩写
#' @param contrast 比较组信息
#' @param gene_type 选择DEP、Up、Down、annot
#' @param map_path png图片储存位置
#' @param outdir 输出路径
#'
#' @return
#' @export 在test_KEGG_Pathway_plot文件夹下
#'
#' @examples
main<-function(diff_df,enrichment_df,pathway_up_color = "red",pathway_down_color = "blue",pathway_annot_color="orange",
               contrast,gene_type,map_path,coord_df =all_data,organism_code,outdir){
  kegg_color_areainfo_df <- get_table_color(enrichment_df = enrichment_df,diff_df = diff_df ,gene_type = gene_type,
                                            pathway_up_color = pathway_up_color,
                                            pathway_down_color = pathway_down_color, 
                                            coord_df = all_data, prefix_pathwayid=organism_code)
  
  areainfo_df <- kegg_color_areainfo_df[[1]]
  kegg_color_df <- kegg_color_areainfo_df[[2]]
  
  get_KEGG_html_file(kegg_color_df = kegg_color_df, areainfo_df = areainfo_df, 
                     map_path = map_path, area_df=all_data,enrichment_df=enrichment_df,
                     suboutdir = outdir , contrast = contrast, prefix_pathwayid = organism_code)
  
}

#diff_df = read.xlsx("result/4_Diff_Expressed/4.1_DiffStats/K_vs_M-DEP_results.xlsx")
conf <- "config/default.conf"
source(conf)
#enrichment_df = read.xlsx("rno_KEGG.xlsx")
pathway_up_color = "red"
pathway_down_color = "blue"
pathway_annot_color="orange"
coord_df = paste0("/PERSONALBIO1/prot4/kegg/",organism_code,"/",organism_code,".coord_info")
map_path = paste0("/PERSONALBIO1/prot4/kegg/",organism_code,"/png/")
all_data = read.delim(coord_df,header = F)
contrast_file <- contrast.file
gene_type = "DEP"
contrast<-read.xlsx(contrast_file)
sapply(contrast$contrast,function(x) {
  start=Sys.time()
  print(x)
  kegg_file<-paste0("./result/5_Enrichment/5.2_KEGG/",x,"/",x,"_KEGG_enrichment.xlsx")
  if(file.exists(kegg_file)){
    enrichment_df = read.xlsx(kegg_file)
    diff_file = paste0("result/4_Diff_Expressed/4.1_DiffStats/",x,"-DEP_results.xlsx")
    diff_df<-read.xlsx(diff_file)
    outdir = paste0("./result/5_Enrichment/5.2_KEGG/",x,"/")
    main(diff_df = diff_df, enrichment_df = enrichment_df, pathway_up_color = pathway_up_color, 
         pathway_down_color = pathway_down_color, contrast = x, gene_type = gene_type,
         map_path = map_path, coord_df = all_data,
         organism_code = organism_code, outdir = outdir)
  }else{
    outdir = paste0("./result/5_Enrichment/5.2_KEGG/",x,"/")
    write.table(x="No enrichment results, not analyzed",file = paste0(outdir,"/未富集到结果不进行map分析.txt"),row.names = F,col.names = F,quote = F)
  }
  
  end<-Sys.time()
  print(end-start)
})





