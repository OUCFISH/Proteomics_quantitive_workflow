#--silently library packages----
sLibrary <- function(packages) {
  invisible(suppressMessages(suppressWarnings(
    sapply(packages, function(x){library(x, character.only = TRUE)})
  )))
}
packages <- c("data.table", "magrittr", "tidyverse",  "Biobase", "Mfuzz",
              "e1071", "ComplexHeatmap", "circlize", "showtext", "rsvg",
              "grid", "clusterProfiler")
sLibrary(packages)

#--functions----
# Get data table of DEGs (id, *_vs_*, Name)
getDEG <- function(indir, contrast) {
  DEGgroups <- fread(contrast, header = T) %>% unlist(.)
  
  # function: extract DEG of each sample
  DEGread <- function(DEGgroup, colVec) {
    DEGpattern <- paste0("^",DEGgroup,"-DEP_results.tsv")
    DEGFn <- dir(indir, pattern = DEGpattern) %>% paste0(indir, "/", .) 
    
    if ("Regulation" %in% colVec) {
      DEGs <- fread(DEGFn, select = c("Protein.Accessions", paste0(DEGgroup,".Regulation")), col.names = c("Protein.Accessions", DEGgroup)) %>% 
        .[, (DEGgroup) := gsub(" Regulation","",get(DEGgroup))]
    }else{
      DEGs <- fread(DEGFn, select = colVec)
    }
    return(DEGs)
  }
  
  # construct DEGdf
  DEGdt <- lapply(DEGgroups, DEGread, colVec = c("Protein.Accessions", "Regulation")) %>% 
    purrr::reduce(full_join, by = "Protein.Accessions") %>% 
    .[apply(.[,-1], 1, function(x) {any(x %in% c("Up", "Down"))}), ] 
  DEGdt[is.na(DEGdt)] <- "NoExp"
  
  # add Gene name
  DEGdt1 <- lapply(DEGgroups, DEGread, colVec = c("Protein.Accessions", "Gene.Name")) %>%
    do.call(rbind, .) %>% 
    unique() %>% 
    left_join(DEGdt, ., by = "Protein.Accessions")
  
  # if(length(DEGIds) > 1000 ) DEGIds <- sample(DEGIds, 1000)
  return(DEGdt1)
}

# Get data frame of expression value
getExpDf <- function(abdFn, mapdt, DEGIds, group, scaleMethod) {
  # read abd table
  expdt <- fread(abdFn, select = c("Protein.Accessions", mapdt$SampleID)) %>%
    .[Protein.Accessions %chin% DEGIds] %>% 
    .[!which(rowSums(.[,-1]) == 0)]
  
  # 是否分组处理
  if (group) {
    grps <- unique(mapdt$Default)
    # 分组求平均值
    sapply(grps, function(grp) { subgroup <- unlist(mapdt[mapdt$Default == grp, "SampleID"])
    expdt[, (grp) := rowMeans(.SD), .SDcols = subgroup]
    }
    ) %>% invisible()
    expdt <- expdt[, .SD, .SDcol = c(colnames(expdt)[1], grps)]
  }
  
  # Dt to Df
  expdf <- setDF(expdt, rownames = expdt$Protein.Accessions)
  expdf <- expdf[, -1]
  
  # log 处理
  expdf <- log2(expdf + 1)
  
  # sampleID to sampleName
  if (!group) {
    colnames(expdf) <- mapdt[match(colnames(expdf), SampleID), SampleName]
  }
  
  # scale(none, Z_Score, mean_centering, pareto) for displaying value's difference well
  if (ncol(expdf) < 3) scaleMethod <- "none"
  # scale
  if (scaleMethod == "Z_Score") {
    expdf[,] <- apply(expdf, 1, function(x){scale(x, center = T, scale = T)}) %>% t()
  }else if (scaleMethod == "mean_centering") {
    expdf[,] <- apply(expdf, 1, function(x){scale(x, center = T, scale = F)}) %>% t()
  }else if (scaleMethod == "pareto") {
    expdf[,] <- apply(expdf, 1, function(x){(x - mean(x))/sqrt(sd(x))})
  }
  
  return(expdf)
}

# sort variables for sourcing
packSourceVars <- function(varPattern){
  vars <- ls(name = .GlobalEnv, pattern = varPattern)
  varParas <- sapply(vars, get, simplify = F)
  return(varParas)
}

# get heatmap body
getHeatmapbody <- function(exp,color.continuous.ls, pp, seed){
  # color
  col_range <- c(floor(min(exp)), 0, ceiling(max(exp)))
  col_fun <- colorRamp2(col_range, color.continuous.ls[[pp$hpb_color]])
  
  # col split
  if (pp$hpb_column_split_num > ncol(exp)) {
    print(paste0("The split number of column(", pp$hpb_column_split_num,
                 ") must be lower than the number of sample(", ncol(exp), ")"))
    pp$hpb_column_split_num <- ncol(exp)
  }
  if (pp$hpb_column_split_num == 1) {
    pp$hpb_column_split_num <- NULL
    plotParas$hpb_column_split_num <<- NULL
  }
  # row split
  if (pp$hpb_row_split_num > nrow(exp)) {
    print(paste0("The split number of row(", pp$hpb_row_split_num,
                 ") must be lower than the number of gene(", nrow(exp), ")"))
    pp$hpb_row_split_num <- 1
  }
  if (pp$hpb_row_split_num == 1) {
    pp$hpb_row_split_num <- NULL
    plotParas$hpb_row_split_num <<- NULL
  }
  
  ParasList <- list(matrix = as.matrix(exp),
                    col = col_fun,
                    show_heatmap_legend = FALSE,
                    row_title = NULL,
                    column_title = NULL,
                    width = unit(pp$hpb_width, "cm"),
                    height = unit(pp$hpb_height, "cm"),
                    
                    cluster_rows = pp$hpb_cluster_rows,
                    clustering_distance_rows = pp$hpb_rowDistanceMethod,
                    clustering_method_rows = pp$hpb_rowClusterMethod,
                    row_dend_side = pp$hpb_row_dend_side,
                    row_dend_width = unit(pp$hpb_row_dend_width, "mm"), 
                    show_row_dend = pp$hpb_show_row_dend,
                    
                    cluster_columns = pp$hpb_cluster_columns,
                    clustering_distance_columns = pp$hpb_colDistanceMethod,
                    clustering_method_columns = pp$hpb_colClusterMethod,
                    column_dend_side = pp$hpb_column_dend_side,
                    column_dend_height = unit(pp$hpb_column_dend_height, "mm"),
                    show_column_dend = pp$hpb_show_column_dend,
                    
                    use_raster = TRUE,
                    show_row_names = FALSE,
                    # row_names_gp = gpar(fontsize = 8),
                    show_column_names = pp$hpb_show_column_names,
                    column_names_gp = gpar(fontsize = pp$hpb_colnames_fontsize,
                                           fontfamily = pp$hpb_colnames_fontfamily),
                    column_names_rot = pp$hpb_colnames_rot,
                    
                    row_split = pp$hpb_row_split_num,
                    column_split = pp$hpb_column_split_num
  )
  
  # col cluster
  if (ncol(exp) <= 2) {
    ParasList$cluster_columns <- FALSE
  }
  
  # row cluster
  if (nrow(exp) <= 2) {
    ParasList$row_rows <- FALSE
  }
  
  # k-means or c-means slicing
  hpSubgroup <- list()
  if (!is.null(pp$hpb_row_split_num) && pp$hpb_splitMethod != "cutree") {
    if (pp$hpb_splitMethod == "k-means" && nrow(exp) > 2) {
      ParasList$row_split <- NULL
      ParasList$row_km <- pp$hpb_row_split_num
      
      # get row order
      set.seed(seed)
      kmhp <- do.call(Heatmap, ParasList)
      pdf(NULL)
      kmhp <- draw(kmhp)
      dev.off()
      rowOrderList <- row_order(kmhp)
      rowCluster <- sapply(names(rowOrderList),
                           function(x){ as.numeric(rep(x, length(rowOrderList[[x]]))) }) %>%
                    unlist()
      
      expDt <- exp[unlist(rowOrderList), ] %>% 
               setDT(., keep.rownames = "ProteinID") %>% 
               .[, Cluster := rowCluster] %>% 
               .[order(Cluster)]  %>% # increase
               .[, Cluster := paste0(pp$anno_geneSlice_label, Cluster)]
      
      # would not use parameter of kmeans, so it was set 1
      ParasList$row_km <- 1
    }else if ( pp$hpb_splitMethod == "c-means" && nrow(exp) > 2) {
      myset <- ExpressionSet(assayData = as.matrix(exp)) %>%
               filter.std(., min.std = 0, visu = FALSE)
      
      # cluster step
      m <- mestimate(myset)
      set.seed(seed)
      mfuzzRes <- exprs(myset) %>% cmeans(., centers = pp$hpb_row_split_num, method = "cmeans", m = m)
      
      # prepare data  for plotting
      expDt <- as.data.frame(myset@assayData$exprs) %>%
               setDT(., keep.rownames = "ProteinID") %>%
               .[, Cluster := mfuzzRes$cluster] %>%
               .[order(Cluster)] %>%
               .[, Cluster := paste0(pp$anno_geneSlice_label, Cluster)]
    }
    
    # subgroup
    ParasList$row_split <- expDt$Cluster
    
    # expression data frame with new order
    data <- expDt[, .SD, .SDcols = -c("ProteinID", "Cluster")] %>%
            setDF(., rownames = expDt$ProteinID)
    ParasList$matrix <- data %>% as.matrix(.)
    hpSubgroup$expDf <- data

    # output expDt for heatmap.xls
    hpSubgroup$expDt <- expDt
    hpSubgroup$subgroup <- expDt$Cluster
  }
  
  # heatmap body plot----
  set.seed(seed)
  hp <- do.call(Heatmap, ParasList)
  
  if (pp$hpb_splitMethod == "cutree") {
    # gene order
    pdf(NULL)
    drawedhp <- draw(hp)
    dev.off()
    rowOrderList <- row_order(drawedhp)
    rowCluster <- sapply(1:pp$hpb_row_split_num, function(x){ rep(x, length(rowOrderList[[x]])) }) %>%
      unlist()
    
    expDt <- exp[unlist(rowOrderList), ] %>%
      setDT(., keep.rownames = "ProteinID") %>%
      .[, Cluster := rowCluster] %>%
      .[order(Cluster)] %>%
      .[, Cluster := paste0(pp$anno_geneSlice_label, Cluster)]

    # output expDt for heatmap.xls
    hpSubgroup$expDt <- expDt
    hpSubgroup$subgroup <- expDt$Cluster
  }

  hpSubgroup$hp <- hp
  return(hpSubgroup)
  # return(hp)
}

# Get mark gene
getMarkGene <- function(indir, contrast, exp, top, type, markGeneStr) {
  DEGgroups <- fread(contrast, header = T) %>% unlist(.)
  # function: extract DEG of each sample
  DEGread <- function(DEGgroup, colVec) {
    DEGpattern <- paste0("^",DEGgroup,"-DEP_results.tsv")
    DEGFn <- dir(indir, pattern = DEGpattern) %>% paste0(indir, "/", .) 
    DEGs <- fread(DEGFn, select = c("Protein.Accessions", "Gene.Name", paste0(DEGgroup,".Log2fc"), paste0(DEGgroup,".Q_value")))
    setnames(DEGs,
             old = c(paste0(DEGgroup,".Log2fc"), paste0(DEGgroup,".Q_value")),
             new = c("Logfc", "Q_value"))
    return(DEGs)
  }
  # construct top DEG data table
  DEGdt <- lapply(DEGgroups, DEGread, colVec = c("Protein.Accessions", "Gene.Name", "Logfc", "Q_value")) %>%
    do.call(rbind, .) %>% 
    .[Protein.Accessions %chin% rownames(exp), ] %>% 
    .[order(Q_value, -abs(Logfc))] %>% 
    unique(.,by = c(1))  

  # top 200 table
  if (!exists("top_200")) {
    if (!file.exists(paste0(outdir, "/top200DEPs.tsv"))) {
      top200DEGdt <- DEGdt[1:200, ]
      top200DEGdt$Q_value <- NULL
      top200DEGdt$Logfc <- NULL
      top200DEGdt$ShowLabel <- top200DEGdt$Gene.Name
      # 兼容DEG结果中没有SwissprotName的情况
      if("SwissProtName" %in% colnames(top200DEGdt)){
        colnames(top200DEGdt) <- c("ProteinID","GeneName","SwissProtName","Custom")
      }else{
        colnames(top200DEGdt) <- c("ProteinID","GeneName","Custom")
      }
      fwrite(top200DEGdt, file = paste0(outdir, "/top200DEPs.tsv"), sep = "\t",
             col.names = TRUE, quote = FALSE)
    }
  }
  
  if (!is.null(markGeneStr)) {
    topDEGdt <- str_split(markGeneStr, ";", simplify = TRUE) %>%
      str_split(., ",", simplify = TRUE) %>%
      as.data.frame()
    
    colnames(topDEGdt) <- c("id", "lab")
  }else{
    # construct top DEG data table
    topDEGdt <- DEGdt %>% .[1:top, ]
  }
  
  genePos <- which(rownames(exp) %in% topDEGdt$id)
  markGeneDf <- data.frame(pos = genePos)
  if (type == "ProteinID") markGeneDf$mark <- topDEGdt$id
  if (type == "GeneName") markGeneDf$mark <- topDEGdt$Name
  if (type == "SwissProtName") markGeneDf$mark <- topDEGdt$SwissProtName
  if (type == "Custom") markGeneDf$mark <- topDEGdt$lab
  
  # markGeneDf (position, markNames)
  return(markGeneDf)
}

# Get term list and output enrich_type_enrichment_G-C*.xls
getTermList <- function(enrich_type, term2geneRd, AnnoFile, exp, subgroup, pAdjustMethod, outdir, pp) {
  # data table of background gene set information
  if (enrich_type == "GO") {
    AnnoDt <- fread(AnnoFile, select = c("Entry","Gene Ontology IDs"), col.names = c("id", "GO")) %>%
              separate_rows(., "GO", sep = ";") %>%
              setDT(.) %>%
              .[GO != "-", ] %>%
              .[GO != "", ]
    if (term2geneRd != "") {
      load("temp/term2gene.RData")
    }else{
      text = paste0("Error: ", term2geneRd, " did not exist")
      cat("\n",text,"\n")
      # 写入 空图片
      empty_plot(text = text, outdir = outdir, saveName = "Error")
    }
    term2gene$Term[is.na(term2gene$Term)] <- term2gene$GO_ID[is.na(term2gene$Term)] # process NA
    terms2gene <- term2gene[, c(1,4)]
    term2name <- term2gene[, c(1,2)]
  }else if (enrich_type == "KEGG"){
    AnnoDt <- fread(AnnoFile, header = TRUE, select = c('pathway id','Protein','pathway description'),
                    col.names = c("term", "id", "termName"))

    terms2gene <- AnnoDt[, c("term", "id")]
    term2name <- AnnoDt[, c("term", "termName")]
  }

  # data table of DEG gene id and subgroup
  ProteinIDsubgrp <- data.table(id = rownames(exp), subgroup = subgroup)

  # function of getting term information of each subgroup
  getTermDf <- function(GC){
    geneList <- intersect(ProteinIDsubgrp[subgroup == GC, id], AnnoDt$id)
    enrichmentRes <- enricher(geneList,
                              TERM2GENE = terms2gene,
                              TERM2NAME = term2name,
                              pAdjustMethod = pAdjustMethod,
                              minGSSize = 1,
                              maxGSSize = 10000,
                              qvalueCutoff = 1,
                              pvalueCutoff = 1)

    GOenrichDf <- enrichmentRes@result[1:termNum, "Description", drop = FALSE]
    colnames(GOenrichDf) <- "text"
    GOenrichDf$col <- structure(color.ls[[termcol]][1:termNum], names = GOenrichDf$text)
    GOenrichDf$fontfamily <- termFontfamily
    GOenrichDf$fontsize <- termFontsize

    # output enrich_type_enrichment_G-C*.xls
    fwrite(enrichmentRes@result, file = paste0(outdir,"/",enrich_type, "_enrichment_", GC, ".xls"),
           sep = "\t", col.names = TRUE, quote = FALSE)

    return(GOenrichDf)
  }

  termNum <- pp$anno_box_term_number
  termcol <- pp$anno_box_term_color
  termFontfamily <- pp$anno_box_term_fontfamily
  termFontsize <- pp$anno_box_term_fontsize
  termInfoList <- sapply(unique(ProteinIDsubgrp$subgroup), getTermDf, simplify = FALSE)
  return(termInfoList)
}

# add annotations onto heatmap--
addAnnotation <- function(hpList, markGeneDf, termList, DEGDt, exp, group, mapdt, color.ls, pp){
  hp <- hpList$hp
  subgroup <- hpList$subgroup
  
  pdf(NULL)
  drawedHp <- draw(hp)
  dev.off()
  
  rightAnnoList <- list()
  leftAnnoList <- list()
  # mark gene
  if (pp$anno_mark_gene) {
    if (pp$anno_mark_gene_fontitalic) {
      fontface <-  "italic"
    }else{
      fontface <-  "plain"
    }
    pdf(NULL)
    markGene <- anno_mark(at = markGeneDf$pos,
                          labels = markGeneDf$mark,
                          labels_gp = gpar(fontface = fontface,
                                           fontsize = pp$anno_mark_gene_fontsize,
                                           fontfamily =  pp$anno_mark_gene_fontfamily
                                           ),
                          side = pp$anno_mark_gene_side,
                          which = "row")
    dev.off()
    
    if (pp$anno_mark_gene_side == "right") {
      rightAnnoList$markGene <- markGene
      annotation_width <- unit(((pp$anno_mark_gene_fontsize / 96) * 2.54) * max(nchar(markGeneDf$mark)), "cm")
      rightAnnoList$annotation_width <- annotation_width
      rightAnnoList$annotation_name_gp <- gpar(fontsize = pp$hpb_legend_title_fontsize,
                                               fontfamily =  pp$hpb_legend_title_fontfamily)
      rightAnnoList$annotation_name_side <- pp$anno_geneSlice_name_side
      rightAnnoList$annotation_name_rot <- pp$hpb_colnames_rot
    }else{
      leftAnnoList$markGene <- markGene
      annotation_width <- unit(((pp$anno_mark_gene_fontsize / 96) * 2.54) * max(nchar(markGeneDf$mark)) + 0.2, "cm")
      leftAnnoList$annotation_width <- annotation_width
      leftAnnoList$annotation_name_gp <- gpar(fontsize = pp$hpb_legend_title_fontsize,
                                              fontfamily =  pp$hpb_legend_title_fontfamily)
      leftAnnoList$annotation_name_side <- pp$anno_geneSlice_name_side
      leftAnnoList$annotation_name_rot <- pp$hpb_colnames_rot
    }
  }
  
  # textbox plot annotation
  if (pp$anno_box) {
    textbox <- ComplexHeatmap::anno_textbox(subgroup, termList,
                                           word_wrap = FALSE,
                                           add_new_line = TRUE,
                                           side = pp$anno_box_side,
                                           text_space = unit(8,"pt"),
                                           line_space = unit(pp$anno_box_term_line_space,"cm"),
                                           max_width = unit(pp$anno_box_max_width, "cm"),
                                           background_gp = grid::gpar(fill = "grey95",
                                                                      col = "grey50"),
                                           by = "anno_link")
    
    if (pp$anno_box_side == "right") {
      rightAnnoList$textbox <- textbox
      rightAnnoList$annotation_name_gp <- gpar(fontsize = pp$hpb_legend_title_fontsize,
                                               fontfamily =  pp$hpb_legend_title_fontfamily)
      rightAnnoList$annotation_name_side <- pp$anno_geneSlice_name_side
      rightAnnoList$annotation_name_rot <- pp$hpb_colnames_rot
    }else{
      leftAnnoList$textbox <- textbox
      leftAnnoList$annotation_name_gp <- gpar(fontsize = pp$hpb_legend_title_fontsize,
                                              fontfamily =  pp$hpb_legend_title_fontfamily)
      leftAnnoList$annotation_name_side <- pp$anno_geneSlice_name_side
      leftAnnoList$annotation_name_rot <- pp$hpb_colnames_rot
    }
  }
  
  # Annotation of Regulation
  if (pp$anno_regulation_anno) {
    regDf <- DEGDt[match(rownames(exp), Protein.Accessions)][,`:=`(Protein.Accessions = NULL,Genes = NULL)] %>% setDF()
    colorNames <- unique(unlist(regDf))
    Regulation <- structure(color.ls[[pp$anno_regulation_color]][1:length(colorNames)], names = colorNames)
    colList <- sapply(colnames(regDf), function(x){ colvec <- Regulation[unique(regDf[,x])] }, simplify = F)
    
    # box anno is not compatible with regulation anno, so they must be in different region !!!
    if (pp$anno_box) {
      if ( pp$anno_box_side == "right") {
        pp$anno_regulation_anno_side <- "left"
      }else{
        pp$anno_regulation_anno_side <- "right"
      }
    }
    if (pp$anno_regulation_anno_side == "right") {
      rightAnnoList$df <- regDf
      rightAnnoList$col <- colList
      rightAnnoList$show_legend <- FALSE
      rightAnnoList$show_annotation_name <- TRUE
      rightAnnoList$annotation_label <- colnames(regDf)
      rightAnnoList$annotation_name_gp <- gpar(fontsize = pp$pp$hpb_legend_title_fontsize,
                                               fontfamily =  pp$hpb_legend_title_fontfamily)
      rightAnnoList$annotation_name_side <- pp$anno_regulaiton_annoName_side
      rightAnnoList$annotation_name_rot <- pp$hpb_colnames_rot
    }else{
      leftAnnoList$df <- regDf
      leftAnnoList$col <- colList
      leftAnnoList$show_legend <- FALSE
      leftAnnoList$show_annotation_name <- TRUE
      leftAnnoList$annotation_label <- colnames(regDf)
      leftAnnoList$annotation_name_gp <- gpar(fontsize = pp$hpb_legend_title_fontsize,
                                              fontfamily =  pp$hpb_legend_title_fontfamily)
      leftAnnoList$annotation_name_side <- pp$anno_regulaiton_annoName_side
      leftAnnoList$annotation_name_rot <- pp$hpb_colnames_rot
    }
  }
  
  # Clusters of Gene
  # bug:if pp$anno_regulation_anno_side == pp$anno_geneSlice_anno_side, 
  # then annotation_name_side will be decided by anno_geneSlice_name_side
  row_split <- pp$hpb_row_split_num
  if (!is.null(row_split) && row_split > 1) {
    rowSliceOrder <- 1:row_split
    if (pp$hpb_splitMethod %in% c("k-means", "c-means")) { 
      rowSliceOrder <- gsub(pp$anno_geneSlice_label, "", names(row_order(drawedHp))) %>% as.numeric(.)
    }
    rowSliceColor <- color.ls[[pp$anno_geneSlice_color]][rowSliceOrder]
    if (pp$anno_geneSlice_tag) {
      GeneCluster <- anno_block(gp = gpar(fill = rowSliceColor),
                                which = "row",
                                labels = paste0(pp$anno_geneSlice_label, rowSliceOrder),
                                labels_gp = gpar(fontsize = pp$hpb_legend_label_fontsize,
                                                 fontfamily = pp$hpb_legend_label_fontfamily),
                                show_name = T)
    }else{
      GeneCluster <- anno_block(gp = gpar(fill = rowSliceColor),
                                which = "row",
                                show_name = T)
    }
    
    if (pp$anno_geneSlice_anno_side == "right") {
      rightAnnoList$annotation_label = c(pp$anno_geneSlice_name, rightAnnoList$annotation_label)
      rightAnnoList$annotation_name_gp <- gpar(fontsize = pp$hpb_legend_title_fontsize,
                                               fontfamily =  pp$hpb_legend_title_fontfamily)
      rightAnnoList$annotation_name_side <- pp$anno_geneSlice_name_side
      rightAnnoList$annotation_name_rot <- pp$hpb_colnames_rot
      rightAnnoList$GeneCluster <- GeneCluster
      rightAnnoList <- rev(rightAnnoList)
    }else{
      leftAnnoList$GeneCluster <- GeneCluster
      leftAnnoList$annotation_label = c(pp$anno_geneSlice_name, leftAnnoList$annotation_label)
      leftAnnoList$annotation_name_gp <- gpar(fontsize = pp$hpb_legend_title_fontsize,
                                              fontfamily =  pp$hpb_legend_title_fontfamily)
      leftAnnoList$annotation_name_side <- pp$anno_geneSlice_name_side
      leftAnnoList$annotation_name_rot <- pp$hpb_colnames_rot
    }
  }
  
  topAnnoList <- list()
  bottomAnnoList <- list()
  # Slices of Sample
  column_split <- pp$hpb_column_split_num
  if (!is.null(column_split) && column_split > 1) {
    sliceOrder <- 1:column_split
    sliceColor <- color.ls[[pp$anno_sampleSlice_color]][sliceOrder]
    if(pp$anno_sampleSlice_tag){
      SampleCluster <- anno_block(gp = gpar(fill = sliceColor),
                                  which = "column",
                                  height = unit(6, "mm"),
                                  labels = paste0(pp$anno_sampleSlice_label, sliceOrder),
                                  labels_gp = gpar(fontsize = pp$hpb_legend_label_fontsize,
                                                   fontfamily = pp$hpb_legend_label_fontfamily),
                                  show_name = T)
    }else{
      SampleCluster <- anno_block(gp = gpar(fill = sliceColor),
                                  which = "column",
                                  height = unit(6, "mm"),
                                  show_name = T)
    }
    
    if (pp$anno_sampleSlice_anno_side == "top") {
      topAnnoList$SampleCluster <- SampleCluster
      topAnnoList$annotation_label <- pp$anno_sampleSlice_name
      topAnnoList$annotation_name_gp <- gpar(fontsize = pp$hpb_legend_title_fontsize,
                                             fontfamily =  pp$hpb_legend_title_fontfamily)
      topAnnoList$annotation_name_side <- pp$anno_sampleSlice_name_side
    }else{
      bottomAnnoList$SampleCluster <- SampleCluster
      bottomAnnoList$annotation_label <-  pp$anno_sampleSlice_name
      bottomAnnoList$annotation_name_gp <- gpar(fontsize = pp$hpb_legend_title_fontsize,
                                                fontfamily =  pp$hpb_legend_title_fontfamily)
      bottomAnnoList$annotation_name_side <- pp$anno_sampleSlice_name_side
    }
  }
  
  # Groups of Sample
  # bug:if pp$anno_sampleSlice_anno_side == pp$anno_group_anno_side, 
  # then annotation_name_side will be decided by anno_group_anno_side
  if (pp$anno_group_anno) {
    if (group) {
      grp <- colnames(exp)
      uniqGrp <- colnames(exp)
    }else{
      grp <- mapdt[match(colnames(exp), SampleName), Default]
      uniqGrp <- unique(grp)
    }

    if (pp$anno_group_tag) {
      Group <- anno_simple(grp, pch = grp,
                           which = "column",
                           # pt_size = unit(1, "snpc")*(0.1*pp$anno_group_tag_size),
                           pt_gp = gpar(fontsize = pp$anno_group_tag_size,
                                        fontfamily = pp$anno_group_tag_fontfamily),
                           col = structure(color.ls[[pp$anno_group_color]][1:length(uniqGrp)],names = uniqGrp)
      )
    }else{
      Group <- anno_simple(grp,
                           which = "column",
                           col = structure(color.ls[[pp$anno_group_color]][1:length(uniqGrp)],names = uniqGrp)
      )
    }
    if (pp$anno_group_anno_side == "top") {
      topAnnoList$Group <- Group
      topAnnoList$annotation_label <- c(topAnnoList$annotation_label, pp$anno_group_annoName)
      topAnnoList$annotation_name_gp <- gpar(fontsize = pp$hpb_legend_title_fontsize,
                                             fontfamily =  pp$hpb_legend_title_fontfamily)
      topAnnoList$annotation_name_side <- pp$anno_group_annoName_side
    }else{
      bottomAnnoList$Group <- Group
      bottomAnnoList$annotation_label <- c(bottomAnnoList$annotation_label, pp$anno_group_annoName)
      bottomAnnoList$annotation_name_gp <- gpar(fontsize = pp$hpb_legend_title_fontsize,
                                                fontfamily =  pp$hpb_legend_title_fontfamily)
      bottomAnnoList$annotation_name_side <- pp$anno_group_annoName_side
    }
  }
  
  # plot
  pdf(NULL)
  if (length(rightAnnoList) > 0) {
    rightHa <- do.call(rowAnnotation, rightAnnoList)
    hp <- attach_annotation(hp, rightHa, side =  "right")
  }
  if (length(leftAnnoList) > 0) {
    leftHa <- do.call(rowAnnotation, leftAnnoList)
    hp <- attach_annotation(hp, leftHa, side =  "left")
  }
  if (length(topAnnoList) > 0) {
    topHa <- do.call(HeatmapAnnotation, topAnnoList)
    hp <- attach_annotation(hp, topHa, side =  "top")
  }
  if (length(bottomAnnoList) > 0) {
    bottomHa <- do.call(HeatmapAnnotation, bottomAnnoList)
    hp <- attach_annotation(hp, bottomHa, side =  "bottom")
  }
  dev.off()
  
  return(hp)
}

# add legends onto heatmap--
addLegend <- function(hp, scaleMethod, DEGDt, exp, group, mapdt, color.continuous.ls, color.ls, pp){
  pdf(NULL)
  drawedHp <- draw(hp)
  dev.off()
  lgdList <- list()
  # legend of heatmap
  hp_lgd <- Legend(
    col_fun = colorRamp2(c(floor(min(exp)), 0, ceiling(max(exp))), color.continuous.ls[[pp$hpb_color]]),
    labels_gp = gpar(fontsize = pp$hpb_legend_label_fontsize, 
                     fontfamily = pp$hpb_legend_label_fontfamily),
    title = switch(scaleMethod,
                   Z_Score = "Z_Score",
                   mean_centering = "mean_centering",
                   pareto = "pareto",
                   none = "fpkm"),
    title_gp =  gpar(fontsize = pp$hpb_legend_title_fontsize, 
                     fontfamily = pp$pp$hpb_legend_title_fontfamily)
  )
  
  # legend of Regulation
  regDf <- DEGDt[match(rownames(exp), Protein.Accessions)][,`:=`(Protein.Accessions = NULL,Genes = NULL)] %>% setDF()
  colorNames <- unique(unlist(regDf))
  if (pp$anno_regulation_anno) {
    lgdList$reg_lgd <- Legend(
      title = pp$anno_regulation_legend_title, 
      title_gp = gpar(fontsize = pp$hpb_legend_title_fontsize, 
                      fontfamily = pp$pp$hpb_legend_title_fontfamily),
      labels =  colorNames,
      labels_gp = gpar(fontsize = pp$hpb_legend_label_fontsize, 
                       fontfamily = pp$hpb_legend_label_fontfamily),
      legend_gp = gpar(fill = structure(color.ls[[pp$anno_regulation_color]][1:length(colorNames)], names = colorNames))
    )
  }
  
  # legend of Group
  if (group) {
    grp <- colnames(exp)
    uniqGrp <- colnames(exp)
  }else{
    grp <- mapdt[match(colnames(exp), SampleName), Default]
    uniqGrp <- unique(grp)
  }
  if (pp$anno_group_anno) {
    lgdList$grp_lgd <- Legend(
      title = pp$anno_group_annoName, 
      title_gp = gpar(fontsize = pp$hpb_legend_title_fontsize, 
                      fontfamily = pp$pp$hpb_legend_title_fontfamily),
      labels =  uniqGrp,
      labels_gp = gpar(fontsize = pp$hpb_legend_label_fontsize, 
                       fontfamily = pp$hpb_legend_label_fontfamily),
      legend_gp = gpar(fill = structure(color.ls[[pp$anno_group_color]][1:length(uniqGrp)],names = uniqGrp))
    )
  }
  
  # Slices of Sample
  column_split <- pp$hpb_column_split_num
  if (!is.null(column_split) && column_split > 1) {
    sliceOrder <- 1:column_split
    sliceColor <- color.ls[[pp$anno_sampleSlice_color]][sliceOrder]
    lgdList$SC_lgd <- Legend(title = pp$anno_sampleSlice_name, 
                             title_gp = gpar(fontsize = pp$hpb_legend_title_fontsize, 
                                             fontfamily = pp$pp$hpb_legend_title_fontfamily),
                             labels = paste0(pp$anno_sampleSlice_label, sliceOrder), 
                             labels_gp = gpar(fontsize = pp$hpb_legend_label_fontsize, 
                                              fontfamily = pp$hpb_legend_label_fontfamily),
                             legend_gp = gpar(fill = sliceColor))
  }
  
  # Slices of gene
  row_split <- pp$hpb_row_split_num
  if (!is.null(row_split) && row_split > 1) {
    rowSliceOrder <- 1:row_split
    rowSliceColor <- color.ls[[pp$anno_geneSlice_color]][rowSliceOrder]
    lgdList$GC_lgd <- Legend(title = pp$anno_geneSlice_name, 
                             title_gp = gpar(fontsize = pp$hpb_legend_title_fontsize, 
                                             fontfamily = pp$pp$hpb_legend_title_fontfamily),
                             labels = paste0(pp$anno_geneSlice_label, rowSliceOrder), 
                             labels_gp = gpar(fontsize = pp$hpb_legend_label_fontsize, 
                                              fontfamily = pp$hpb_legend_label_fontfamily),
                             legend_gp = gpar(fill = rowSliceColor))
  }
  
  pdf(NULL)
  hp <- draw(hp, heatmap_legend_list = list(hp_lgd), annotation_legend_list = lgdList)
  dev.off()
  
  return(hp)
}

# calculate size of heatmap--
calc_ht_size <-  function(ht, unit = "inch") {
  pdf(NULL)
  ht <- draw(ht)
  w <- ComplexHeatmap:::width(ht)
  w <- convertX(w, unit, valueOnly = TRUE)
  h <- ComplexHeatmap:::height(ht)
  h <- convertY(h, unit, valueOnly = TRUE)
  dev.off()
  c(w, h)
}

main <- function(conf) {
  if (conf != "") {
    # source("11.3_heatmap_Enrichment.conf")
    source(conf)
  }
  # 创建outdir
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE, mode = "0755")
  
  # 加载通用工具包
  set_theme_module <- paste0(utils_root, "/set_Theme.R")
  source(set_theme_module)
  figure_Save_module <- paste0(utils_root, "/figure_Save_zhaoyu.R")
  source(figure_Save_module)
  
  AddFont(font_path)
  
  showtext_auto()
  
  # 获取分组，map：input/group.tsv，差异对照组；contrast：input/contrast.tsv，样品分组。
  DEGgroups <- fread(contrast, header = T) %>% unlist(.) %>% strsplit(., split = "_vs_") %>% unlist(.) %>% unique(.)
  mapdt <- fread(map, col.names = c("SampleID", "SampleName","Default")) %>% .[Default %chin% DEGgroups, ]
  
  # 从差异表中获取DEG信息
  DEGDt <- getDEG(indir, contrast)
  # abdFn：input/fpkm.tsv；group：FALSE;scaleMethod："Z_Score"
  expDf <- getExpDf(abdFn, mapdt, DEGIds = DEGDt$Protein.Accessions, group, scaleMethod)
  
  plotParas <<- packSourceVars(varPattern = "^hpb_*|^anno_*")
  
  ht_opt(message = FALSE)
  hpSubgroup <- getHeatmapbody(exp = expDf, color.continuous.ls = color.continuous.ls,
                               pp = plotParas, seed = 123)
  if(!is.null(hpSubgroup$expDf)) expDf <- hpSubgroup$expDf
  if(!is.null(hpSubgroup$expDt)) expDt <- hpSubgroup$expDt

  # output heatmap.xls
  if ( hpb_row_split_num > 1 ){
    fwrite(expDt, file = paste0(outdir, "/heatmap.xls"), sep = "\t",
           col.names = TRUE, quote = FALSE)
  }else{
    expDf_out <- data.frame(ProteinID = rownames(expDf), expDf)
    fwrite(expDf_out, file = paste0(outdir, "/heatmap.xls"), sep = "\t",
           col.names = TRUE, quote = FALSE)
  }

  # term List for plotting text box
  if ( enrich_type == "GO") AnnoFile <- GOFile
  if ( enrich_type == "KEGG") AnnoFile <- KEGGFile
  termList <- getTermList(enrich_type = enrich_type, term2geneRd = term2geneRd, AnnoFile = AnnoFile,
                          exp = expDf, subgroup = hpSubgroup$subgroup,
                          pAdjustMethod = pAdjustMethod, outdir = outdir, pp = plotParas)
  
  # get mark gene and top200DEGs.tsv
  if (anno_mark_gene_type != "Custom") anno_mark_gene_str <- NULL
  markGeneDf <- getMarkGene(indir = indir, contrast = contrast, exp = expDf,
                            top = anno_mark_gene_top, type = anno_mark_gene_type,
                            markGeneStr = anno_mark_gene_str)

  # if add mark gene
  if (!anno_mark_gene) markGeneDf <- NULL
  
  hp <- addAnnotation(hpList = hpSubgroup, markGeneDf = markGeneDf,
                      termList = termList, DEGDt = DEGDt,
                      exp = expDf, group = group, mapdt = mapdt, color.ls = color.ls, pp = plotParas)
  
  hp <-  addLegend(hp = hp, scaleMethod = scaleMethod, DEGDt = DEGDt, exp = expDf, group = group, mapdt = mapdt, 
                   color.continuous.ls = color.continuous.ls, color.ls = color.ls,
                   pp = plotParas)
  
  hp_size <- calc_ht_size(ht = hp)
  save_plot(p = hp, outdir = outdir, saveName = "heatmapEnrichment",
            width = hp_size[1] + 1,
            height = hp_size[2] + 1, noggsave = TRUE)
}

# running----
argv <- commandArgs(TRUE)
if (length(argv) == 1) {
    print("run script with conf.")
    main(argv[1])
} else if (length(argv) == 0) {
    print("run script without conf.")
}


#debug 
# conf <- "./config/11.3_heatmap_Enrichment_analysis.conf"
# main(conf)