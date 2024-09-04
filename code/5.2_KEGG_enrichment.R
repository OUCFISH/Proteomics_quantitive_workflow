# -*- coding: utf-8 -*-
# @Time    : 2024/7/31 16:02
# @Author  : liangxing.sun
# @FileName: 5.2_KEGG_html.R
# @Software: Rstudio
# @Email: liangxing.sun@personalbio.cn
rm(list = ls())
pacman::p_load(clusterProfiler,dplyr,stringr,tidyr,openxlsx)
#' Title 单纯富集脚本
#'
#' @param DEP_file 差异表，需要四列"Protein.Accessions","KEGG","Regulation","Gene.Name"注释表需要三列"Protein.Accessions","KEGG","Gene.Name"
#' @param type diff,annot,差异表还是注释表，生成不同格式的富集结果
#' @param bg_annotation_file 背景文件 #"/PERSONALBIO1/prote/database/KEGG_Pathway/Animals/rno_KO_pathway_description.csv"
#' @param filename 生成结果文件
#' @param outdir 生成结果文件路径
#' @param org 物种简写
#'
#' @return
#' @export 富集表
#'
#' @examples main(DEP_file,type = "diff", bg_list = bg_list, bg_annotation_file =bg_annotation_file , filename = filename,outdir = outdir)
main<-function(DEP_file,type = "diff", bg_annotation_file = bg_annotation_file, org = org, filename = filename,outdir = outdir){
  bg_annotation<-read.csv(bg_annotation_file)
  if(org != "ko"){
    bg = bg_annotation[,c("pathway.id","GeneID")]
    bg$KEGGID<-paste(org,bg$GeneID,sep=":")
    bg = bg[,c(1,3)]
  }else{
    bg = bg_annotation[,c("pathway.id","KO.ID")]
  }
  
  term<-bg_annotation[,c("pathway.id","pathway.description","pathway.classA","pathway.classB")]
  colnames(term)<-c("ID","Description","level1","level2")
  term<-unique(term)
  DEP_data<-read.xlsx(DEP_file)
  if(!any(grepl("Log2fc",x = colnames(DEP_data))) && type == "diff"){
    type = "annot"
  }
  if(type == "diff"){
    tg<-select(DEP_data,contains(c("Protein.Accessions","KEGG","Regulation","Gene.Name")))%>%
      separate_rows(KEGG,sep = ";")%>%filter(KEGG!="")
    colnames(tg)<-c("Protein.Accessions","KEGG","Regulation","Gene.Name")
    tg<-tg%>%replace_na(list(Gene.Name = "",KEGG = ""))
    enrichment<-enricher(gene =tg$KEGG,pvalueCutoff = 1,
                         pAdjustMethod = "BH",minGSSize = 1,
                         maxGSSize = 100000,TERM2GENE =bg,
                         TERM2NAME = term[,c(1,2)])
    if(!is.null(enrichment)){
      Regu<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
          left_join(.,tg,by=join_by(geneID==KEGG),relationship = "many-to-many")
      if(length(which(Regu$Regulation=="Up")) != 0){
        Up_KEGG<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
          left_join(.,filter(tg,Regulation=="Up"),by=join_by(geneID==KEGG),relationship = "many-to-many")%>%
          filter(!is.na(Regulation))%>%
          mutate(Up_KEGG=paste0(Protein.Accessions,"(",geneID,")"))%>%
          group_by(ID)%>%
          aggregate(Up_KEGG~ID,function(x) paste(x,collapse = ";"))
        Up_protein<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
          left_join(.,filter(tg,Regulation=="Up"),by=join_by(geneID==KEGG),relationship = "many-to-many")%>%
          filter(!is.na(Regulation))%>%
          mutate(Up_protein=paste0(Protein.Accessions,"(",Gene.Name,")"))%>%
          group_by(ID)%>%
          aggregate(Up_protein~ID,function(x) paste(x,collapse = ";"))
      }else{
        Up_KEGG<-data.frame("ID" = enrichment@result$ID, "Up_KEGG" = "-")
        Up_protein <- data.frame("ID" = enrichment@result$ID, "Up_protein" = "-")
      }
      
      if(length(which(Regu$Regulation=="Down")) != 0){
        Down_KEGG<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
          left_join(.,filter(tg,Regulation=="Down"),by=join_by(geneID==KEGG),relationship = "many-to-many")%>%
          filter(!is.na(Regulation))%>%
          mutate(Down_KEGG=paste0(Protein.Accessions,"(",geneID,")"))%>%
          group_by(ID)%>%
          aggregate(Down_KEGG~ID,function(x) paste(x,collapse = ";"))
        Down_protein<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
          left_join(.,filter(tg,Regulation=="Down"),by=join_by(geneID==KEGG),relationship = "many-to-many")%>%
          filter(!is.na(Regulation))%>%
          mutate(Down_protein=paste0(Protein.Accessions,"(",Gene.Name,")"))%>%
          group_by(ID)%>%
          aggregate(Down_protein~ID,function(x) paste(x,collapse = ";"))
      }else{
        Down_KEGG<-data.frame("ID" = enrichment@result$ID, "Down_KEGG" = "-")
        Down_protein <- data.frame("ID" = enrichment@result$ID, "Down_protein" = "-")
      }
      
      enrich_result<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
        left_join(.,tg,by=join_by(geneID==KEGG),relationship = "many-to-many")%>%
        mutate(Protein=paste0(Protein.Accessions,"(",geneID,")"))%>%
        group_by(ID)%>%
        aggregate(Protein~ID,function(x) paste(x,collapse = ";"))%>%
        left_join(.,enrichment@result,by=join_by(ID==ID),relationship = "many-to-many")%>%
        left_join(.,select(term,ID,level1,level2),by=join_by(ID==ID))%>%
        left_join(.,Up_KEGG,by=join_by(ID==ID),relationship = "many-to-many")%>%
        left_join(.,Down_KEGG,by=join_by(ID==ID),relationship = "many-to-many")%>%
        left_join(.,Up_protein,by=join_by(ID==ID),relationship = "many-to-many")%>%
        left_join(.,Down_protein,by=join_by(ID==ID),relationship = "many-to-many")%>%
        replace_na(list(Up_KEGG="-",Down_KEGG="-",Up_protein="-",Down_protein="-"))%>%
        mutate(Up = case_when(Up_protein=="-"~0,T~str_count(Up_protein,pattern=";")+1),
               Down = case_when(Down_protein=="-"~0,T~str_count(Down_protein,pattern=";")+1),
               DEP = case_when(Protein=="-"~0,T~str_count(Protein,pattern=";")+1),
               Total = as.numeric(gsub("(/).*", "", BgRatio)),
               RichFactor = as.numeric(gsub("(/).*", "", GeneRatio))/as.numeric(gsub("(/).*", "", BgRatio)))%>%
        rename(Pvalue=pvalue,adjustPvalue=p.adjust)%>%
        select(ID,Description,level1,level2,Up,Down,DEP,Total,Pvalue,adjustPvalue,GeneRatio,BgRatio,Up_protein,Down_protein,RichFactor,Protein,Up_KEGG,Down_KEGG)%>%
        arrange(Pvalue)
    }else{
      enrich_result<-enrichment
    }
    
  }else{
    tg<-select(DEP_data,contains(c("Protein.Accessions","KEGG","Gene.Name")))%>%
      separate_rows(KEGG,sep = ";")%>%filter(KEGG!="")
    colnames(tg)<-c("Protein.Accessions","KEGG","Gene.Name")
    tg<-tg%>%replace_na(list(Gene.Name = "",KEGG = ""))
    enrichment<-enricher(gene =tg$KEGG,pvalueCutoff = 1,
                         pAdjustMethod = "BH",minGSSize = 1,
                         maxGSSize = 100000,TERM2GENE =bg,
                         TERM2NAME = term[,c(1,2)])
    if(!is.null(enrichment)){
      enrich_result<-enrichment@result%>%separate_rows(geneID,sep = "/")%>%
        left_join(.,tg,by=join_by(geneID==KEGG),relationship = "many-to-many")%>%
        mutate(Protein=paste0(Protein.Accessions,"(",geneID,")"))%>%
        group_by(ID)%>%
        aggregate(Protein~ID,function(x) paste(x,collapse = ";"))%>%
        left_join(.,enrichment@result,by=join_by(ID==ID),relationship = "many-to-many")%>%
        left_join(.,select(term,ID,level1,level2),by=join_by(ID==ID))%>%
        mutate(DEP = case_when(Protein=="-"~0,T~str_count(Protein,pattern=";")+1),
               Total = as.numeric(gsub("(/).*", "", BgRatio)),
               RichFactor = as.numeric(gsub("(/).*", "", GeneRatio))/as.numeric(gsub("(/).*", "", BgRatio)))%>%
        rename(Pvalue=pvalue,adjustPvalue=p.adjust)%>%
        select(ID,Description,level1,level2,DEP,Total,Pvalue,adjustPvalue,GeneRatio,BgRatio,RichFactor,Protein)%>%
        arrange(Pvalue)
    }else{
      enrich_result<-enrichment
    }
    
  }
  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = T)
  }
  file = paste(outdir,filename,sep="/")
  if(is.null(enrich_result)){
    write.table(x="Not enriched for results",file = paste0(outdir,"/未富集到结果.txt"),row.names = F,col.names = F,quote = F)
  }else{
    write.xlsx(enrich_result,file = file)
  }
}

conf <- "config/default.conf"
source(conf)
contrast <- contrast.file
outdir <- "result/5_Enrichment/5.2_KEGG/"
contrast_data<-read.xlsx(contrast)
type = "diff"
if(dir.exists("./background")){
  bg_annotation_file = "./background/KEGG_db.csv"
}else{
  bg_annotation_file = Sys.glob(paste0(kegg_file_path,organism_category,"/",organism_code, "_KO_pathway_description.csv"))
  #"/PERSONALBIO1/prote/database/KEGG_Pathway/Animals/rno_KO_pathway_description.csv"
}
sapply(contrast_data$contrast,function(x){
	DEP_file = paste0("./result/4_Diff_Expressed/4.1_DiffStats/",x,"-DEP_results.xlsx")
	filename = paste0(x,"_KEGG_enrichment.xlsx")
	outdir = paste0(outdir,"/",x,"/")
	if(file.exists(bg_annotation_file)){
		main(DEP_file,type = type,org = organism_code,
		     bg_annotation_file = bg_annotation_file, filename = filename,outdir = outdir)
	}else{
		if(!dir.exists("result/5_Enrichment/5.2_KEGG/")){
			dir.create("result/5_Enrichment/5.2_KEGG/",recursive = T)
		}
		write.table(x="No background information on this species",file="result/5_Enrichment/5.2_KEGG/无背景文件不支持此分析.txt",row.names = F,col.names = F,quote = F)

	}
})

