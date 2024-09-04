# -*- coding: utf-8 -*-
# @Time    : 2024/08/06 16:48
# @Author  : liangxing.sun
# @FileName: 5.2_KEGG_html.R
# @Software: Rstudio
# @Email: liangxing.sun@personalbio.cn

rm(list = ls())
pacman::p_load(clusterProfiler,dplyr,stringr,tidyr,openxlsx,parallel)
start = Sys.time()

#' Title 单次读取背景文件（没有并行）
#'
#' @param db_file 文件路径一般是./background/InterPro.tsv
#' @param enrich_type 富集类型c("InterPro","Reactome","DO","WikiPathway")其中一个
#'
#' @return 返回一个数据框内容是背景文件整理之后的
#' @export
#'
#' @examples
read_and_process_data <- function(db_file,enrich_type=enrich_type) {
  
  if(enrich_type == "InterPro"){
    db_df <- read.delim(db_file,header = F)
    colnames(db_df)<-c("Protein_ID","MD5.disest","length"	,"source.database","database.ID","Description_1","start","stop",	"e-value"	,"status"	,"Date",	"ID","Description"	,"GO","Pathway")
    db_df <- db_df %>%
      select(ID, Description,  Protein_ID,source.database) %>%
      filter(ID != "-")
    db_df$Protein_ID <- sapply(strsplit(db_df$Protein_ID, split = "\\|"), FUN = function(x) x[2])
  }else if(enrich_type == "Reactome"){
    db_df <- read.delim(db_file,header = F)
    colnames(db_df)<-c("Protein_ID","ID","URL","Description","GeneName","organism")
    organism_Scientific_name1<-gsub("_"," ",organism_Scientific_name)
    db_df<-db_df%>%filter(organism == organism_Scientific_name1)
    db_df<-db_df[,c("ID","Description","Protein_ID","GeneName")]
  }else if(enrich_type == "DO"){
    db_df<-read.delim(db_file ,skip = 15)
    organism_Scientific_name1<-gsub("_"," ",organism_Scientific_name)
    db_df<-db_df%>%filter(SpeciesName == organism_Scientific_name1,DBobjectType == "gene")%>%
      select(DOID,DOtermName,DBObjectSymbol)
    colnames(db_df)<-c("ID","Description","GeneName")
    db_df <- db_df[,c("ID","Description","GeneName")]
  }else if(enrich_type == "WikiPathway"){
    db_df <- read.delim(db_file,header = T)
    db_df<-db_df%>%select(wpid,name,Entry,gene)%>%
      rename(ID = wpid,Description = name,Protein_ID = Entry,Gene_ID=gene)
  }
  
  return(db_df)
}

# 
#' Title 进行富集分析
#'
#' @param diff_prot 差异表/注释表
#' @param db_df 背景文件
#' @param type diff/annot 表示差异或者注释
#'
#' @return 返回一个数据框，可以直接保存为excel
#' @export
#'
#' @examples
perform_enrichment_analysis <- function(diff_prot, db_df,type,enrich_type) {
  bg <- db_df[, c(1, 3)] %>% unique()
  term <- db_df[,c(1,2)] %>% unique()
  if(!any(grepl("Log2fc",x = colnames(diff_prot))) && type == "diff"){
    type = "annot"
  }
  if(enrich_type == "DO"){
    if(type == "diff"){
      DEP_data<-diff_prot
      tg<-select(DEP_data,contains(c("Protein.Accessions","Regulation","Gene.Name")))
      colnames(tg)<-c("Protein.Accessions","Regulation","Gene.Name")
      tg<-tg%>%replace_na(list(Gene.Name = ""))%>%filter(Gene.Name != "")
      enrichment<-enricher(gene =tg$Gene.Name,pvalueCutoff = 1,
                           pAdjustMethod = "BH",minGSSize = 1,
                           maxGSSize = 100000,
                           TERM2GENE = bg,
                           TERM2NAME = term)
      if(!is.null(enrichment)){
	Regu<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
          left_join(.,tg,by=join_by(geneID==Gene.Name),relationship = "many-to-many")
        if(length(which(Regu$Regulation=="Up")) != 0){
          Up_protein<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
            left_join(.,filter(tg,Regulation=="Up"),by=join_by(geneID==Gene.Name),relationship = "many-to-many")%>%
            filter(!is.na(Regulation))%>%
            mutate(Up_protein=paste0(Protein.Accessions,"(",geneID,")"))%>%
            group_by(ID)%>%
            aggregate(Up_protein~ID,function(x) paste(x,collapse = ";"))
        }else{
          Up_protein <- data.frame("ID" = enrichment@result$ID, "Up_protein" = "-")
        }
        
        if(length(which(Regu$Regulation=="Down")) != 0){
          Down_protein<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
            left_join(.,filter(tg,Regulation=="Down"),by=join_by(geneID==Gene.Name),relationship = "many-to-many")%>%
            filter(!is.na(Regulation))%>%
            mutate(Down_protein=paste0(Protein.Accessions,"(",geneID,")"))%>%
            group_by(ID)%>%
            aggregate(Down_protein~ID,function(x) paste(x,collapse = ";"))
        }else{
          Down_protein <- data.frame("ID" = enrichment@result$ID, "Down_protein" = "-")
        }
        
        enrich_result<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
          left_join(.,tg,by=join_by(geneID==Gene.Name),relationship = "many-to-many")%>%
          mutate(Protein=paste0(Protein.Accessions,"(",geneID,")"))%>%
          group_by(ID)%>%
          aggregate(Protein~ID,function(x) paste(x,collapse = ";"))%>%
          left_join(.,enrichment@result,by=join_by(ID==ID),relationship = "many-to-many")%>%
          left_join(.,Up_protein,by=join_by(ID==ID),relationship = "many-to-many")%>%
          left_join(.,Down_protein,by=join_by(ID==ID),relationship = "many-to-many")%>%
          replace_na(list(Up_protein="-",Down_protein="-"))%>%
          mutate(Up = case_when(Up_protein=="-"~0,T~str_count(Up_protein,pattern=";")+1),
                 Down = case_when(Down_protein=="-"~0,T~str_count(Down_protein,pattern=";")+1),
                 DEP = case_when(Protein=="-"~0,T~str_count(Protein,pattern=";")+1),
                 Total = as.numeric(gsub("(/).*", "", BgRatio)),
                 RichFactor = as.numeric(gsub("(/).*", "", GeneRatio))/as.numeric(gsub("(/).*", "", BgRatio)))%>%
          rename(Pvalue=pvalue,adjustPvalue=p.adjust)%>%
          select(ID,Description,Up,Down,DEP,Total,Pvalue,adjustPvalue,Up_protein,Down_protein,GeneRatio,BgRatio,RichFactor,Protein)%>%
          arrange(Pvalue)
      }else{
        enrich_result<-enrichment
      }
    }else{
      DEP_data<-diff_prot
      tg<-select(DEP_data,contains(c("Protein.Accessions","Gene.Name")))
      colnames(tg)<-c("Protein.Accessions","Gene.Name")
      tg<-tg%>%replace_na(list(Gene.Name = ""))%>%filter(Gene.Name!="")
      enrichment<-enricher(gene =tg$Gene.Name,pvalueCutoff = 1,
                           pAdjustMethod = "BH",minGSSize = 1,
                           maxGSSize = 100000,TERM2GENE =bg,
                           TERM2NAME = term)
      if(!is.null(enrichment)){
        enrich_result<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
          left_join(.,tg,by=join_by(geneID==Gene.Name),relationship = "many-to-many")%>%
          mutate(Protein=paste0(Protein.Accessions,"(",geneID,")"))%>%
          group_by(ID)%>%
          aggregate(Protein~ID,function(x) paste(x,collapse = ";"))%>%
          left_join(.,enrichment@result,by=join_by(ID==ID),relationship = "many-to-many")%>%
          mutate(DEP = case_when(Protein=="-"~0,T~str_count(Protein,pattern=";")+1),
                 Total = as.numeric(gsub("(/).*", "", BgRatio)),
                 RichFactor = as.numeric(gsub("(/).*", "", GeneRatio))/as.numeric(gsub("(/).*", "", BgRatio)))%>%
          rename(Pvalue=pvalue,adjustPvalue=p.adjust)%>%
          select(ID,Description,DEP,Total,Pvalue,adjustPvalue,GeneRatio,BgRatio,RichFactor,Protein)%>%
          arrange(Pvalue)
      }else{
        enrich_result<-enrichment
      }
      
    }
  }else{
    if(type == "diff"){
      DEP_data<-diff_prot
      tg<-select(DEP_data,contains(c("Protein.Accessions","Regulation","Gene.Name")))
      colnames(tg)<-c("Protein.Accessions","Regulation","Gene.Name")
      tg<-tg%>%replace_na(list(Gene.Name = ""))
      enrichment<-enricher(gene =tg$Protein.Accessions,pvalueCutoff = 1,
                           pAdjustMethod = "BH",minGSSize = 1,
                           maxGSSize = 100000,
                           TERM2GENE = bg,
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
          left_join(.,Up_protein,by=join_by(ID==ID),relationship = "many-to-many")%>%
          left_join(.,Down_protein,by=join_by(ID==ID),relationship = "many-to-many")%>%
          replace_na(list(Up_protein="-",Down_protein="-"))%>%
          mutate(Up = case_when(Up_protein=="-"~0,T~str_count(Up_protein,pattern=";")+1),
                 Down = case_when(Down_protein=="-"~0,T~str_count(Down_protein,pattern=";")+1),
                 DEP = case_when(Protein=="-"~0,T~str_count(Protein,pattern=";")+1),
                 Total = as.numeric(gsub("(/).*", "", BgRatio)),
                 RichFactor = as.numeric(gsub("(/).*", "", GeneRatio))/as.numeric(gsub("(/).*", "", BgRatio)))%>%
          rename(Pvalue=pvalue,adjustPvalue=p.adjust)%>%
          select(ID,Description,Up,Down,DEP,Total,Pvalue,adjustPvalue,Up_protein,Down_protein,GeneRatio,BgRatio,RichFactor,Protein)%>%
          arrange(Pvalue)
      }else{
        enrich_result<-enrichment
      }
      
    }else{
      DEP_data<-diff_prot
      tg<-select(DEP_data,contains(c("Protein.Accessions","Gene.Name")))
      colnames(tg)<-c("Protein.Accessions","Gene.Name")
      tg<-tg%>%replace_na(list(Gene.Name = ""))
      enrichment<-enricher(gene =tg$Protein.Accessions,pvalueCutoff = 1,
                           pAdjustMethod = "BH",minGSSize = 1,
                           maxGSSize = 100000,TERM2GENE =bg,
                           TERM2NAME = term)
      if(!is.null(enrichment)){
        enrich_result<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
          left_join(.,tg,by=join_by(geneID==Protein.Accessions),relationship = "many-to-many")%>%
          mutate(Protein=paste0(geneID,"(",Gene.Name,")"))%>%
          group_by(ID)%>%
          aggregate(Protein~ID,function(x) paste(x,collapse = ";"))%>%
          left_join(.,enrichment@result,by=join_by(ID==ID),relationship = "many-to-many")%>%
          mutate(DEP = case_when(Protein=="-"~0,T~str_count(Protein,pattern=";")+1),
                 Total = as.numeric(gsub("(/).*", "", BgRatio)),
                 RichFactor = as.numeric(gsub("(/).*", "", GeneRatio))/as.numeric(gsub("(/).*", "", BgRatio)))%>%
          rename(Pvalue=pvalue,adjustPvalue=p.adjust)%>%
          select(ID,Description,DEP,Total,Pvalue,adjustPvalue,GeneRatio,BgRatio,RichFactor,Protein)%>%
          arrange(Pvalue)
      }else{
        enrich_result<-enrichment
      }
      
    }
  }
  return(enrich_result)
}


#' Title 保存分析结果到Excel文件
#'
#' @param enrichment_final 富集结果数据框
#' @param output_dir 输出路径
#' @param contrast 比较组信息
#' @param enrich_type 富集类型c("InterPro","Reactome","DO","WikiPathway")之一
#'
#' @return
#' @export
#'
#' @examples
save_enrichment_results <- function(enrichment_final, output_dir,contrast,enrich_type) {
  if (!dir.exists(paste0(output_dir,"/",contrast,"/"))) {
    dir.create(paste0(output_dir,"/",contrast,"/"), recursive = TRUE)
  }
  if(is.null(enrichment_final)){
    write.table(x="Not enriched for results",file = paste0(output_dir,"/",contrast,"/未富集到结果.txt"),row.names = F,col.names = F,quote = F)
  }else{
    write.xlsx(enrichment_final, file = paste0(output_dir,"/",contrast,"/",contrast ,"_",enrich_type,"_enrichment.xlsx"))
  }
  
}


#' Title 富集背景文件没有时输出文件内容为error
#'
#' @param enrich_type 富集类型c("InterPro","Reactome","DO","WikiPathway")
#' @param contrast 比较组
#'
#' @return
#' @export 写出同名称表格但内容为没有背景文件信息
#'
#' @examples
creat_error_xlsx<-function(enrich_type,contrast){
  if(enrich_type == "InterPro"){
    outdir<-"./result/5_Enrichment/5.3_domain/"
  }else if(enrich_type == "Reactome"){
    outdir<-"./result/5_Enrichment/5.4_Reactome/"
  }else if(enrich_type == "DO"){
    outdir<-"./result/5_Enrichment/5.5_DO/"
  }else if(enrich_type == "WikiPathway"){
    outdir<-"./result/5_Enrichment/5.6_WikiPathway/"
  }
  if (!dir.exists(paste0(outdir,"/"))) {
    dir.create(paste0(outdir,"/"), recursive = TRUE)
  }
  error_data<-data.frame("error_info"="No background information on this species, analysis failed")
  write.table("No background information on this species",file = paste0(outdir,"/无背景文件不支持此分析.txt"),row.names=F,col.names=F,quote=F)
}

#' Title 单个富集类型进行富集
#'
#' @param diff_prot_file 差异表/注释表
#' @param db_df 背景文件数据框
#' @param enrich_type 富集类型c("InterPro","Reactome","DO","WikiPathway")
#' @param type diff/annot
#' @param contrast 比较组
#'
#' @return
#' @export
#'
#' @examples
all_function<-function(diff_prot_file,db_df,enrich_type,type,contrast){
  library(clusterProfiler)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(openxlsx)
  if(enrich_type == "InterPro"){
    outdir<-"./result/5_Enrichment/5.3_domain/"
  }else if(enrich_type == "Reactome"){
    outdir<-"./result/5_Enrichment/5.4_Reactome/"
  }else if(enrich_type == "DO"){
    outdir<-"./result/5_Enrichment/5.5_DO/"
  }else if(enrich_type == "WikiPathway"){
    outdir<-"./result/5_Enrichment/5.6_WikiPathway/"
  }
  if(length(db_df)==0){
    creat_error_xlsx(enrich_type,contrast)
  }else{
    diff_prot <- read.xlsx(diff_prot_file)
    enrichment_final <- perform_enrichment_analysis(diff_prot = diff_prot, db_df = db_df, type = type,enrich_type=enrich_type)
    save_enrichment_results(enrichment_final, output_dir = outdir, contrast = contrast, enrich_type = enrich_type)
  }
  
}

#' Title 将一个比较组进行并行富集
#'
#' @param contrast 比较组信息
#' @param no_cores 核数
#' @param db_ls 背景文件列表，其中包括每个背景文件c("InterPro","Reactome","DO","WikiPathway")
#'
#' @return
#' @export
#'
#' @examples
all_main<-function(contrast,no_cores,db_ls){
  library(parallel)
  cl<-makeCluster(no_cores,type="FORK")
  diff_prot_file<-paste0("./result/4_Diff_Expressed/4.1_DiffStats/",contrast,"-DEP_results.xlsx")
  parSapply(cl, c("InterPro","Reactome","DO","WikiPathway"),function(x) 
    all_function(diff_prot_file,db_df = db_ls[[x]],enrich_type = x,type,contrast))
  stopCluster(cl)
}

config<-"./config/default.conf"
source(config)
type<-"diff"
no_cores <- 4
contrast_file<-"./temp/contrast.xlsx"
contrast<-read.xlsx(contrast_file)
db_ls<-list()

for (i  in c("InterPro","Reactome","DO","WikiPathway")) {
  if(dir.exists("./background")){
    db_file<-paste0("./background/",i,"_db.tsv")
  }else{
    if(i == "InterPro"){
      db_file = Sys.glob(paste0(interpro_file_path, organism_Scientific_name, "/*tsv"))
    }else if(i == "Reactome"){
      db_file = "/PERSONALBIO1/prote/database/Reactome/UniProt2Reactome_All_Levels.txt"
    }else if(i == "DO"){
      db_file = "/PERSONALBIO1/prote/database/DO_db/DISEASE-ALLIANCE_COMBINED_6.tsv"
    }else if(i== "WikiPathway"){
      db_file = Sys.glob(paste0(WikiPathway_file_path,organism_Scientific_name,".xls"))
    }
  }
  if(!file.exists(db_file)){
    db_ls[[i]]<-NULL
  }else{
    db_ls[[i]]<-read_and_process_data(db_file = db_file,enrich_type = i)
  }
}


sapply(contrast$contrast,function(x) {print(x) 
       all_main(x,no_cores,db_ls)})

end = Sys.time()
print(end - start)




