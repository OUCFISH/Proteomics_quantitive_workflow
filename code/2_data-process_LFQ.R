# 清除工作区
rm(list = ls())
# 加载必要的库
pacman::p_load(openxlsx,tidyverse,impute,Biostrings,Peptides,SeqKnn,readxl,dplyr)
# 定义文件路径和参数
conf <- "./config/default.conf"
source(conf)
uniprot_file <- Sys.glob(paste0(uniprot_file_path,organism_category,"/*",organism_code, "_", organism_ID, "_annotation.xlsx"))
fasta_file <- list.files(path = "data", pattern = "\\.fasta$", full.names = TRUE)
kegg_file <- Sys.glob(paste0(kegg_file_path,organism_category,"/",organism_code,"_KO_pathway_description.csv"))
interpro_file <- Sys.glob(paste0(interpro_file_path,organism_Scientific_name,"/*",".tsv"))
output1 <- "./result/2_Data_Summary/"
if (!dir.create(output1)) {
  dir.create(output1,recursive = T)
}
output2 <- "./result/4_Diff_Expressed/4.1_DiffStats/"
if (!dir.create(output2)) {
  dir.create(output2,recursive = T)
}
group.file <- read.xlsx("./temp/group.xlsx")
#group_num <- length(unique(group.file$Group))
#sample_num <- length(group.file$Group) / length(unique(group.file$Group))
#nonsamplecol <- 4
group.name <- unique(group.file$Group)
contrast.file <- read.xlsx("./temp/contrast.xlsx")
compare.name <- contrast.file$contrast
foldchange1 <- Fold_Change1
foldchange2 <- Fold_Change2
dif_p <- pvalue
##读入蛋白数据
protein.data.ori <- read_excel(protein_file)
peptide.data.ori <- read_excel(peptide_file)
tag_relate <- read.xlsx("data/Input Files.xlsx")
sample.name <- setNames(tag_relate$File.Name,tag_relate$File.ID)

# 蛋白表保存列名
selected_procol <- c("Protein FDR Confidence: Combined","Accession","Description","Exp. q-value: Combined","Sum PEP Score","Coverage [%]","# Peptides",
                  "# PSMs","# Unique Peptides","# AAs","MW [kDa]","calc. pI","Score Sequest HT: Sequest HT","# Peptides (by Search Engine): Sequest HT",
                  "Gene Symbol","# Protein Pathway Groups","# Razor Peptides")

# 肽段表保存列名
selected_pepcol <- c("Annotated Sequence","Modifications","PEP","q-Value","# Protein Groups","# Proteins","# PSMs","Master Protein Accessions",
                     "Positions in Master Proteins","# Missed Cleavages","Theo. MH+ [Da]","Quan Info","Confidence (by Search Engine): Sequest HT",
                     "q-Value (by Search Engine): Sequest HT","PEP (by Search Engine): Sequest HT","XCorr (by Search Engine): Sequest HT")

process_colnames <- function(data,selected_col,sample.name){
  newdata <- select(data,c(all_of(selected_col),starts_with("Abundances (Scaled):")))
  names(newdata) <- names(newdata) %>% gsub("Abundances \\(Scaled\\): ","",.)%>%gsub(": Sample","",.)
  names(newdata)[match(names(sample.name), names(newdata))] <- sample.name
  return(newdata)
}

protein.data <- process_colnames(data = protein.data.ori,selected_col = selected_procol,sample.name = sample.name)
peptide.data <- process_colnames(data = peptide.data.ori,selected_col = selected_pepcol,sample.name = sample.name)
write.xlsx(protein.data,file = paste0(output1,"ProteinGroup.xlsx"))
write.xlsx(peptide.data,file = paste0(output1,"Peptide.xlsx"))

protein.data <- protein.data %>%  filter(`# Unique Peptides` >= 1, `Score Sequest HT: Sequest HT` > 0) %>% 
  dplyr::rename(Protein.Description = Description,Gene.Name = `Gene Symbol`,Protein.Accessions=Accession,
                `Molecular Weight[KDa]`=`MW [kDa]`,Isoelectric.Point=`calc. pI`,Coverage=`Coverage [%]`)
# 修改列名
# 过滤蛋白质（至少在一个组中有50%以上的样本非零）
keep_rows <- apply(protein.data, 1, function(row) {
  all(sapply(group.name, function(group) {
    # group_cols <- grep(paste0("^",group,"\\d+"), colnames(protein.data), value = TRUE)
    group_cols <- which(colnames(protein.data) %in% group.file$SampleName[group.file$Group == group])
    non_zero_count <- sum(!is.na(row[group_cols]))
    sample_num <- length(group.file[group.file$Group==group,"SampleID"])
    non_zero_count < 0.5 * sample_num
  }))
})
data_filtered <- protein.data[!keep_rows, ]
# 使用mixed imputation strategy填补缺失值
## 最小值的二分之一
nosample_col <- setdiff(names(data_filtered),sample.name)
intensity_matrix <- as.matrix(select(data_filtered, !c(nosample_col))) # 去除前几列非样本列
min_nonzero_value <- min(intensity_matrix[intensity_matrix > 0], na.rm = TRUE)
half_min_value <- min_nonzero_value / 2
# 定义填补缺失值的函数
fill_na <- function(group_data, half_min_value) {
  col_names <- colnames(group_data)
  # 找到全是 NA 的行并用 half_min_value 填充
  all_na_rows <- apply(group_data, 1, function(row) all(is.na(row)))
  for (name in col_names) {
    group_data[all_na_rows, name] <- half_min_value  + runif(1,0,0.001)
  }
  # 使用 KNN 填充剩余的缺失值
  group_data_filled <- as.data.frame(SeqKNN(as.matrix(group_data), k = 10))
  return(group_data_filled)
}
# 对每个组进行填充
protein_data_impute_group <- lapply(group.name, function(group) {
  # group_indices <- grep(pattern = paste0("^", group, "\\d+"), x = colnames(data_filtered))
  group_indices <- which(colnames(data_filtered) %in% group.file$SampleName[group.file$Group == group])
  group_intensity <- data_filtered[, group_indices]
  protein_data_group <- fill_na(group_intensity, half_min_value)
  return(protein_data_group)
})
# 检查colnames(data_filtered)是否包含"QC"
if (any(grepl("QC", colnames(data_filtered)))) {
  # 获取QC列的索引
  qc_indices <- grep(pattern = "QC", x = colnames(data_filtered))
  qc_intensity <- data_filtered[, qc_indices]
  qc_data_impute <- fill_na(qc_intensity, half_min_value)
  
  # 合并填补后的数据
  protein_data_impute_merge <- do.call(cbind, protein_data_impute_group)
  protein_data_impute <- data.frame(data_filtered[, c(1:nonsamplecol)], protein_data_impute_merge)
  protein_data_impute_addqc <- data.frame(data_filtered[, c(1:nonsamplecol)], protein_data_impute_merge, qc_data_impute)
  pca.data <- protein_data_impute_addqc[, c("Protein.Accessions",sample.name)]
} else {
  # 如果不包含"QC"，则仅合并填补后的数据
  protein_data_impute_merge <- do.call(cbind, protein_data_impute_group)
  protein_data_impute <- cbind(data_filtered[, c(nosample_col)],protein_data_impute_merge)
  pca.data <- protein_data_impute[, c("Protein.Accessions",sample.name)]
}
write.xlsx(pca.data, "temp/pca_data.xlsx")

## 按照分组计算 Mean
for (i in 1:length(group.name)) {
  group.col <- paste0("Mean_",group.name[i])
  protein_data_impute <- protein_data_impute %>% 
    mutate(!!group.col := rowMeans(dplyr::select(., which(colnames(protein_data_impute) %in% group.file$SampleName[group.file$Group == group.name[i]])), na.rm = TRUE))
}
# 匹配Uniprot注释并保存结果
annotate_proteins <- function(data, uniprot_file, kinase_file, organism, output_file) {
  uniprot_db <- read.xlsx(uniprot_file)
  uniprot_db_filter <- uniprot_db %>% dplyr::select(Entry, Length,`Gene.Ontology.(GO)`, KEGG)
  data_annotated <- left_join(data, uniprot_db_filter, by = c("Protein.Accessions" = "Entry"))
  data_annotated <- data_annotated %>% select(Protein.Accessions, Gene.Name, Protein.Description,Length,`Molecular Weight[KDa]`,Isoelectric.Point,Coverage, everything())
  data_annotated$Length <- as.numeric(data_annotated$Length)
  # kegg:K-number,pathway_id
  kegg_db <- read.csv(kegg_file)
  kegg_db$GeneID <- paste(organism_code,kegg_db$GeneID,sep=":")
  data_annotated$KO.ID <- sapply(data_annotated$KEGG,function(x){
    paste(kegg_db$KO.ID[which(kegg_db$GeneID %in% unlist(strsplit(x,";")))],collapse = ";")
  })
  data_annotated$pathway.id <- sapply(data_annotated$KEGG,function(x){
    paste(kegg_db$pathway.id[which(kegg_db$GeneID %in% unlist(strsplit(x,";")))],collapse = ";")
  })
  #interpro
  domain_db <- read.delim2(interpro_file)
  colnames(domain_db)<-c("Protein_ID","MD5.disest","length"	,"source.database","database.ID","Description_1","start","stop",	"e-value"	,"status"	,"Date",	"ID","Description"	,"GO","Pathway")
  domain_db$Protein_ID <- sapply(strsplit(domain_db$Protein_ID,split = "\\|"),FUN = function(x)x[2])
  domain_db <- domain_db %>%
    select(ID, Description, Protein_ID,source.database) %>%
    filter(ID != "-")
  data_annotated$InterPro.id <- sapply(data_annotated$Protein.Accessions,function(x){
    paste(domain_db$ID[which(domain_db$Protein_ID %in% x)],collapse = ";")
  })
  # Subcellular database
  subdb <- read.csv(Sys.glob(paste0(subcellular_file_path, organism_Scientific_name, "/*csv")))
  subdb$UniprotID <- sapply(strsplit(subdb$Protein_ID,split = "\\|"),FUN = function(x)x[2])
  sub.db <- dplyr::select(subdb,UniprotID,Localizations)
  sub.db$Localizations <- sapply(strsplit(sub.db$Localizations,split = "\\|"),FUN = function(x)x[1])
  data_annotated$subcellular_location <- sub.db$Localizations[base::match(data_annotated$Protein.Accessions,sub.db$UniprotID)]
  # Kinase database
  kinase.database <- read.xlsx(kinase_file,startRow = 3)
  kinase.db <- read.xlsx(Sys.glob(paste0(kinase_file_path,organism_category,"/",organism_code,"_",organism_ID,"_Kinase.xlsx")))
  if (organism %in% kinase.database$KIN_ORGANISM) {
    kinase.prot <- unique(kinase.database[kinase.database$KIN_ORGANISM == organism,"KIN_ACC_ID"])
    kinase.prot.total <- c(kinase.prot,kinase.db$Entry)
  }else{
    kinase.prot.total <- kinase.db$Entry
  }
  data_annotated$kinase <- ifelse(test = data_annotated$Protein.Accessions %in% kinase.prot.total,yes = "TRUE",no = "FALSE")
  # TF database
  TFdb <- read.table(Sys.glob(paste0(TF_file_path,organism,"/","*_TF.txt")),header = T,sep = "\t")
  TFdb <- TFdb[,c("Symbol","Family")]
  data_annotated$TF <- TFdb$Family[base::match(data_annotated$Gene.Name,TFdb$Symbol)]
  # DO database
  DOdb <- read.delim("/PERSONALBIO1/prote/database/DO_db/DISEASE-ALLIANCE_COMBINED_6.tsv",skip = 15)
  organism_Scientific_name1<-gsub("_"," ",organism_Scientific_name)
  DOdb <- DOdb %>% filter(SpeciesName == organism_Scientific_name1,DBobjectType == "gene")%>%
    select(DOID,DBObjectSymbol)
  data_annotated$DO <- DOdb$DOID[base::match(data_annotated$Gene.Name,DOdb$DBObjectSymbol)]
  # Wikipathway database
  Wikipathway.db <- read.delim(Sys.glob(paste0(WikiPathway_file_path,organism_Scientific_name,".xls")))
  Wikipathway.db <- Wikipathway.db[,c("wpid","Entry")]
  Wikipathway.db.long <- Wikipathway.db %>%
    group_by(Entry) %>%
    summarise(wpid = paste(wpid, collapse = ";"), .groups = 'drop')
  data_annotated <- data_annotated %>%
    left_join(Wikipathway.db.long, by = c("Protein.Accessions" = "Entry"))
  colnames(data_annotated)[which(names(data_annotated) == "wpid")] <- "Wikipathway"
  
  write.xlsx(data_annotated, file = output_file)
}
annotate_proteins_diff <- function(data, uniprot_file, kinase_file, organism, output_file) {
  uniprot_db <- read.xlsx(uniprot_file)
  uniprot_db_filter <- uniprot_db %>% dplyr::select(Entry,`Gene.Ontology.(GO)`, KEGG, Pfam, InterPro, eggNOG)
  data_annotated <- left_join(data, uniprot_db_filter, by = c("Protein.Accessions" = "Entry"))
  new_col_order <- c(
    "Protein.Accessions", "Gene.Name", "Protein.Description", 
    setdiff(colnames(data_annotated), c( "Protein.Accessions", "Gene.Name", "Protein.Description"))
  )
  data_annotated <- data_annotated[, new_col_order]
  # kegg:K-number,pathway_id
  kegg_db <- read.csv(kegg_file)
  kegg_db$GeneID <- paste(organism_code,kegg_db$GeneID,sep=":")
  data_annotated$KO.ID <- sapply(data_annotated$KEGG,function(x){
    paste(kegg_db$KO.ID[which(kegg_db$GeneID %in% unlist(strsplit(x,";")))],collapse = ";")
  })
  data_annotated$pathway.id <- sapply(data_annotated$KEGG,function(x){
    paste(kegg_db$pathway.id[which(kegg_db$GeneID %in% unlist(strsplit(x,";")))],collapse = ";")
  })
  #interpro
  domain_db <- read.delim2(interpro_file)
  colnames(domain_db)<-c("Protein_ID","MD5.disest","length"	,"source.database","database.ID","Description_1","start","stop",	"e-value"	,"status"	,"Date",	"ID","Description"	,"GO","Pathway")
  domain_db$Protein_ID <- sapply(strsplit(domain_db$Protein_ID,split = "\\|"),FUN = function(x)x[2])
  domain_db <- domain_db %>%
    select(ID, Description, Protein_ID,source.database) %>%
    filter(ID != "-")
  data_annotated$InterPro.id <- sapply(data_annotated$Protein.Accessions,function(x){
    paste(domain_db$ID[which(domain_db$Protein_ID %in% x)],collapse = ";")
  })
  # Subcellular database
  subdb <- read.csv(Sys.glob(paste0(subcellular_file_path, organism_Scientific_name, "/*csv")))
  subdb$UniprotID <- sapply(strsplit(subdb$Protein_ID,split = "\\|"),FUN = function(x)x[2])
  sub.db <- dplyr::select(subdb,UniprotID,Localizations)
  sub.db$Localizations <- sapply(strsplit(sub.db$Localizations,split = "\\|"),FUN = function(x)x[1])
  data_annotated$subcellular_location <- sub.db$Localizations[base::match(data_annotated$Protein.Accessions,sub.db$UniprotID)]
  # Kinase database
  kinase.database <- read.xlsx(kinase_file,startRow = 3)
  kinase.db <- read.xlsx(Sys.glob(paste0(kinase_file_path,organism_category,"/",organism_code,"_",organism_ID,"_Kinase.xlsx")))
  if (organism %in% kinase.database$KIN_ORGANISM) {
    kinase.prot <- unique(kinase.database[kinase.database$KIN_ORGANISM == organism,"KIN_ACC_ID"])
    kinase.prot.total <- c(kinase.prot,kinase.db$Entry)
  }else{
    kinase.prot.total <- kinase.db$Entry
  }
  data_annotated$kinase <- ifelse(test = data_annotated$Protein.Accessions %in% kinase.prot.total,yes = "TRUE",no = "FALSE")
  # TF database
  TFdb <- read.table(Sys.glob(paste0(TF_file_path,organism,"/","*_TF.txt")),header = T,sep = "\t")
  TFdb <- TFdb[,c("Symbol","Family")]
  data_annotated$TF <- TFdb$Family[base::match(data_annotated$Gene.Name,TFdb$Symbol)]
  # DO database
  DOdb <- read.delim("/PERSONALBIO1/prote/database/DO_db/DISEASE-ALLIANCE_COMBINED_6.tsv",skip = 15)
  organism_Scientific_name1<-gsub("_"," ",organism_Scientific_name)
  DOdb <- DOdb %>% filter(SpeciesName == organism_Scientific_name1,DBobjectType == "gene")%>%
    select(DOID,DBObjectSymbol)
  data_annotated$DO <- DOdb$DOID[base::match(data_annotated$Gene.Name,DOdb$DBObjectSymbol)]
  # Wikipathway database
  Wikipathway.db <- read.delim(Sys.glob(paste0(WikiPathway_file_path,organism_Scientific_name,".xls")))
  Wikipathway.db <- Wikipathway.db[,c("wpid","Entry")]
  Wikipathway.db.long <- Wikipathway.db %>%
    group_by(Entry) %>%
    summarise(wpid = paste(wpid, collapse = ";"), .groups = 'drop')
  data_annotated <- data_annotated %>%
    left_join(Wikipathway.db.long, by = c("Protein.Accessions" = "Entry"))
  colnames(data_annotated)[which(names(data_annotated) == "wpid")] <- "Wikipathway"
  write.xlsx(data_annotated, file = output_file)
}
## total proteins
## 按照对比策略计算FC log2FC pvalue Qvalue
protein_data_impute.total <- list()
for (j in 1:length(compare.name)) {
  print(compare.name[j])
  compare.group <- unlist(strsplit(x = compare.name[j],split = "_vs_"))
  mean1 <- paste0("Mean_",compare.group[1])
  mean2 <- paste0("Mean_",compare.group[2])
  # 计算 Fold Change (FC) 和 Log2 Fold Change (Log2FC)
  protein_data_impute <- protein_data_impute %>%
    mutate(FC := .[[mean1]] / .[[mean2]],
           Log2fc := log2(.[[mean1]] / .[[mean2]]))
  # subgroup1 <- grep(paste0("^",compare.group[1],"\\d+"), colnames(protein_data_impute), value = TRUE)
  # subgroup2 <- grep(paste0("^",compare.group[2],"\\d+"), colnames(protein_data_impute), value = TRUE)
  subgroup1 <- colnames(protein_data_impute)[which(colnames(protein_data_impute) %in% group.file$SampleName[group.file$Group == compare.group[1]])]
  subgroup2 <- colnames(protein_data_impute)[which(colnames(protein_data_impute) %in% group.file$SampleName[group.file$Group == compare.group[2]])]
  df_a <- protein_data_impute[which(colnames(protein_data_impute) %in% subgroup1)]
  df_b <- protein_data_impute[which(colnames(protein_data_impute) %in% subgroup2)]
  protein_data_impute$p_value <- sapply(1:nrow(df_a), function(k) {
    col_a <- as.numeric(df_a[k,])
    col_b <- as.numeric(df_b[k,])
    if(mean(col_b) == mean(col_a)) {
      return(1)
    } else {
      var_equal <- var.test(col_a, col_b)$p.value > 0.05
      return(t.test(col_a, col_b, var.equal = var_equal)$p.value)
    }
  })
  protein_data_impute$Q_value <- p.adjust(protein_data_impute$p_value, method = "BH", 
                                          n = length(protein_data_impute$p_value))
  protein_data_impute$Regulation <- ifelse(protein_data_impute$p_value > dif_p, "Nodiff",
                                           ifelse(protein_data_impute$FC > foldchange1, "Up",
                                                  ifelse(protein_data_impute$FC < foldchange2, "Down", "Nodiff")))
  colnames(protein_data_impute)[which(colnames(protein_data_impute) == "FC")] <- paste0(compare.group[1],"_vs_",compare.group[2],".FC")
  colnames(protein_data_impute)[which(colnames(protein_data_impute) == "Log2fc")] <- paste0(compare.group[1],"_vs_",compare.group[2],".Log2fc")
  colnames(protein_data_impute)[which(colnames(protein_data_impute) == "p_value")] <- paste0(compare.group[1],"_vs_",compare.group[2],".pvalue")
  colnames(protein_data_impute)[which(colnames(protein_data_impute) == "Q_value")] <- paste0(compare.group[1],"_vs_",compare.group[2],".Q_value")
  colnames(protein_data_impute)[which(colnames(protein_data_impute) == "Regulation")] <- paste0(compare.group[1],"_vs_",compare.group[2],".Regulation")
  protein_data_impute.total[[j]] <- protein_data_impute
  ## diff prot
  protein_data_impute_diff <- select(protein_data_impute,c(all_of(c(nosample_col,mean1,mean2,subgroup1,subgroup2)),starts_with(compare.name[j])))
  protein_data_impute_diff_reg <- protein_data_impute_diff[which(colnames(protein_data_impute_diff) == paste0(compare.group[1],"_vs_",compare.group[2],".Regulation"))]
  protein_data_impute_diff_final <- protein_data_impute_diff[protein_data_impute_diff_reg != "Nodiff",]
  annotate_proteins_diff(protein_data_impute_diff_final, uniprot_file, kinase_file, organism, paste0(output2,compare.name[j],"-DEP_results.xlsx"))
  annotate_proteins_diff(protein_data_impute_diff, uniprot_file, kinase_file, organism, paste0("temp/",compare.name[j],"-total_results.xlsx"))
}
protein_data_impute.total.merge <- do.call(cbind, protein_data_impute.total)
duplicated_cols <- duplicated(colnames(protein_data_impute.total.merge))
protein_data_impute.total.final <- protein_data_impute.total.merge[, !duplicated_cols]
annotate_proteins(protein_data_impute.total.final, uniprot_file, kinase_file, organism, paste0(output1,"Protein_SummaryStat.xlsx"))



