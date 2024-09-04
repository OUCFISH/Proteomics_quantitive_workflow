# -*- coding: utf-8 -*-
# @Time    : 2024/7/31 16:02
# @Author  : liangxing.sun
# @FileName: 5.2_KEGG_html.R
# @Software: Rstudio
# @Email: liangxing.sun@personalbio.cn

rm(list = ls())
pacman::p_load(clusterProfiler,dplyr,stringr,tidyr,openxlsx)

#' Title GO富集脚本
#'
#' @param DEP_file 差异文件或者注释文件
#' @param GO_annotation_file GO背景文件
#' @param filename 输出富集表名称
#' @param outdir 输出路径
#' @param type 指定富集类型是差异还是注释diff/annot
#'
#' @return
#' @export 指定路径下的富集表
#'
#' @examples
main<-function(DEP_file,GO_annotation_file,filename,outdir,type){
  diff_data = read.xlsx(DEP_file)
  bg_annotation = read.delim(GO_annotation_file)
  colnames(bg_annotation)<-c("ID","Description","Category","ProteinID","Gene.Name")
  bg<-bg_annotation[,c("ID","ProteinID")]
  term<-bg_annotation[,c("ID","Description","Category")]
  term<-unique(term)
  DEP_data<-read.xlsx(DEP_file)
  if(!any(grepl("Log2fc",x = colnames(DEP_data))) && type == "diff"){
    type = "annot"
  }
  if(type == "diff"){
    tg<-select(DEP_data,contains(c("Protein.Accessions","Regulation","Gene.Name")))
    colnames(tg)<-c("Protein.Accessions","Regulation","Gene.Name")
    tg<-tg%>%replace_na(list(Gene.Name = ""))
    enrichment<-enricher(gene =tg$Protein.Accessions,pvalueCutoff = 1,
                         pAdjustMethod = "BH",minGSSize = 1,
                         maxGSSize = 100000,TERM2GENE =bg,qvalueCutoff = 1,
                         TERM2NAME = term)
    if(!is.null(enrichment)){
      Regu<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
          left_join(.,tg,by=join_by(geneID==Protein.Accessions),relationship = "many-to-many")
      if(length(which(Regu$Regulation=="Up")) != 0){
        Up_protein<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
          left_join(.,filter(tg,Regulation=="Up"),by=join_by(geneID==Protein.Accessions),relationship = "many-to-many")%>%
          filter(!is.na(Regulation))%>%
          mutate(Up_protein=paste0(geneID,"(",Gene.Name,")"))%>%
          group_by(ID)%>%
          aggregate(Up_protein~ID,function(x) paste(x,collapse = ";"))
      }else{
        Up_protein <- data.frame("ID" = enrichment@result$ID, "Up_protein" = "-")
      }
      
      if(length(which(Regu$Regulation=="Down")) != 0){
        Down_protein<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
          left_join(.,filter(tg,Regulation=="Down"),by=join_by(geneID==Protein.Accessions),relationship = "many-to-many")%>%
          filter(!is.na(Regulation))%>%
          mutate(Down_protein=paste0(geneID,"(",Gene.Name,")"))%>%
          group_by(ID)%>%
          aggregate(Down_protein~ID,function(x) paste(x,collapse = ";"))
      }else{
        Down_protein <- data.frame("ID" = enrichment@result$ID, "Down_protein" = "-")
      }
      
      enrich_result<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
        left_join(.,tg,by=join_by(geneID==Protein.Accessions),relationship = "many-to-many")%>%
        mutate(Protein=paste0(geneID,"(",Gene.Name,")"))%>%
        group_by(ID)%>%
        aggregate(Protein~ID,function(x) paste(x,collapse = ";"))%>%
        left_join(.,enrichment@result,by=join_by(ID==ID),relationship = "many-to-many")%>%
        left_join(.,select(term,ID,Category),by=join_by(ID==ID))%>%
        left_join(.,Up_protein,by=join_by(ID==ID),relationship = "many-to-many")%>%
        left_join(.,Down_protein,by=join_by(ID==ID),relationship = "many-to-many")%>%
        replace_na(list(Up_protein="-",Down_protein="-"))%>%
        mutate(Up = case_when(Up_protein=="-"~0,T~str_count(Up_protein,pattern=";")+1),
               Down = case_when(Down_protein=="-"~0,T~str_count(Down_protein,pattern=";")+1),
               DEP = case_when(Protein=="-"~0,T~str_count(Protein,pattern=";")+1),
               Total = as.numeric(gsub("(/).*", "", BgRatio)),
               RichFactor = as.numeric(gsub("(/).*", "", GeneRatio))/as.numeric(gsub("(/).*", "", BgRatio)))%>%
        rename(Pvalue=pvalue,adjustPvalue=p.adjust)%>%
        select(ID,Description,Category,Up,Down,DEP,Total,Pvalue,adjustPvalue,Up_protein,Down_protein,GeneRatio,BgRatio,RichFactor,Protein)%>%
        arrange(Pvalue)
    }else{
      enrich_result<-enrichment
    }
    
  }else{
    tg<-select(DEP_data,contains(c("Protein.Accessions","Gene.Name")))
    colnames(tg)<-c("Protein.Accessions","Gene.Name")
    tg<-tg%>%replace_na(list(Gene.Name = ""))
    enrichment<-enricher(gene =tg$Protein.Accessions,pvalueCutoff = 1,
                         pAdjustMethod = "BH",minGSSize = 1,qvalueCutoff = 1,
                         maxGSSize = 100000,TERM2GENE =bg,
                         TERM2NAME = term)
    if(!is.null(enrichment)){
      enrich_result<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
        left_join(.,tg,by=join_by(geneID==Protein.Accessions),relationship = "many-to-many")%>%
        mutate(Protein=paste0(geneID,"(",Gene.Name,")"))%>%
        group_by(ID)%>%
        aggregate(Protein~ID,function(x) paste(x,collapse = ";"))%>%
        left_join(.,enrichment@result,by=join_by(ID==ID),relationship = "many-to-many")%>%
        left_join(.,select(term,ID,Category),by=join_by(ID==ID))%>%
        mutate(DEP = case_when(Protein=="-"~0,T~str_count(Protein,pattern=";")+1),
               Total = as.numeric(gsub("(/).*", "", BgRatio)),
               RichFactor = as.numeric(gsub("(/).*", "", GeneRatio))/as.numeric(gsub("(/).*", "", BgRatio)))%>%
        rename(Pvalue=pvalue,adjustPvalue=p.adjust)%>%
        select(ID,Description,Category,DEP,Total,Pvalue,adjustPvalue,GeneRatio,BgRatio,RichFactor,Protein)%>%
        arrange(Pvalue)
    }else{
      enrich_result<-enrichment
    }
    
  }
  file = paste0(outdir,"/",filename)
  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = T)
  }
  if(is.null(enrich_result)){
    write.table(x="Not enriched for results",file = paste0(outdir,"/未富集到结果.txt"),row.names = F,col.names = F,quote = F)
  }else{
    write.xlsx(enrich_result,file = file)
  }
  
  return(enrichment)
}

GO_DAG <- function(GO_enrichment, Result_enrichment, suboutdir, contrast) {
  prepareMapping <- function(GO_enrichment, GOsep = "/") {
    mat_df <- as.data.frame(GO_enrichment)
    mat_df <- mat_df[mat_df$geneID != "-", ]
    map <- mat_df$geneID
    names(map) <- gsub(" ", "", mat_df$ID)
    map_list <- lapply(map, function(x) gsub(" ", "", strsplit(x, split = GOsep)[[1]]))
    return(map_list)
  }
  dir.create(paste0(suboutdir, "/DAG"), showWarnings = F, recursive = T, mode = "0755")
  
  print("有向无环图")
  
  GO_enrichment@result$ONTOLOGY <- Result_enrichment$Category
  GO2geneID <- prepareMapping(GO_enrichment)
  GO_enrichment <- as.data.frame(GO_enrichment)
  GO_enrichment <- na.omit(GO_enrichment)
  ## 将单个类型的数据分开
  ## 获取类型的list
  type_list <- unique(GO_enrichment$ONTOLOGY)
  for (gotype in type_list) {
    type_df <- GO_enrichment[which(GO_enrichment$ONTOLOGY == gotype), ]
    type_df <- type_df[, c("ONTOLOGY", colnames(type_df)[1:ncol(type_df) - 1])]
    if (nrow(type_df) > 1) {
      geneID_all <- unlist(lapply(type_df$geneID, function(x) unlist(strsplit(x, "/"))))
    } else {
      geneID_all <- unlist(strsplit(type_df$geneID, "/"))
    }
    
    GEO_df <- new("enrichResult",
                  result = type_df, gene = geneID_all, pvalueCutoff = 0.01, pAdjustMethod = "BH",
                  qvalueCutoff = 0.05, ontology = "BP", keytype = "GID", universe = "Unknown", geneSets = GO2geneID, organism = "Unknown", readable = FALSE
    )
    GEO_df@ontology <- gotype
    # suppressMessages(goplot(ego_df, showCategory = 10)) ## igraph 布局方式的有向无环图；如有需要可进行替换
    ## 保存图片 pdf、png、svg格式
    
    save_file <- paste0(suboutdir, "/DAG/", contrast, "_", gotype, "_dag")
    
    pdf(paste0(save_file, ".pdf"))
    suppressMessages(clusterProfiler::plotGOgraph(GEO_df))
    dev.off()
  }
}

conf <- "config/default.conf"
source(conf)
contrast <- contrast.file
outdir <- "result/5_Enrichment/5.1_GO/"
contrast_data<-read.xlsx(contrast)
type = "diff"
if(dir.exists("./background")){
  GO_annotation_file = "./background/GO_db.tsv"
}else{
  GO_annotation_file = Sys.glob(paste0(uniprot_file_path,organism_category,"/*",organism_code, "_", organism_ID, "_godb.tsv"))
  #"/PERSONALBIO1/prote/database/uniprot_annotation/Animals/Rattus_norvegicus_(Norway_rat)_rno_10116_godb.tsv"
}
sapply(contrast_data$contrast,function(x){
  DEP_file = paste0("./result/4_Diff_Expressed/4.1_DiffStats/",x,"-DEP_results.xlsx")
  print(x)
  filename = paste0(x,"_GO_enrichment.xlsx")
  outdir = paste0(outdir,"/",x,"/")
  if(file.exists(GO_annotation_file)){
    GO_enrichment<-main(DEP_file,type = type,
       GO_annotation_file = GO_annotation_file, filename = filename,outdir = outdir)
    if(!is.null(GO_enrichment)){
      data<-read.xlsx(paste0(outdir,"/",filename))
      data<-data%>%rename(pvalue=Pvalue,p.adjust=adjustPvalue)
      GO_DAG(GO_enrichment,data,suboutdir = outdir,contrast = x)
      dag_file<-Sys.glob(paste0(outdir,"/DAG/*"))
      if(length(dag_file)==0){
        write.table(x="Too few enrichment results without graphs",file = paste0(outdir,"/DAG/富集结果过少没有图形.txt"),row.names = F,col.names = F,quote = F)
      }
    }
    
  }else{
    if(!dir.exists("result/5_Enrichment/5.1_GO/")){
      dir.create("result/5_Enrichment/5.1_GO/",recursive = T)
    }
    write.table(x="No background information on this species",file="result/5_Enrichment/5.1_GO/无背景文件不支持此分析.txt",row.names = F,col.names = F,quote = F)
  }
})







