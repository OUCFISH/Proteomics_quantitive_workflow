rm(list = ls())
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(circlize)))
suppressWarnings(suppressMessages(library(ComplexHeatmap)))
# args<-commandArgs(trailingOnly = T)
conf <- "config/default.conf"
source(conf)
source(paste0(utils_root,"/set_Theme.R"))
contrast <- contrast.file
go_dir <- "result/5_Enrichment/5.1_GO/"
gene_dir <- indir
outdir <- "result/5_Enrichment/5.1_GO/"
contrast_data <- read.xlsx(contrast)
GObarplot<-function(go.data,output,contrast){
  go.bp <- go.data %>%
    filter(Category == "BP") %>%
    arrange(Pvalue)
  go.cc <- go.data %>%
    filter(Category == "CC") %>%
    arrange(Pvalue)
  go.mf <- go.data %>%
    filter(Category == "MF") %>%
    arrange(Pvalue)
  go.bp.top <- go.bp[1:min(top_n,nrow(go.bp)), ] %>% arrange(desc(DEP))
  go.cc.top <- go.cc[1:min(top_n,nrow(go.cc)), ] %>% arrange(desc(DEP))
  go.mf.top <- go.mf[1:min(top_n,nrow(go.mf)), ] %>% arrange(desc(DEP))
  go.bp.top$Term <- factor(go.bp.top$Term,levels = go.bp.top$Term)
  go.cc.top$Term <- factor(go.cc.top$Term,levels = go.cc.top$Term)
  go.mf.top$Term <- factor(go.mf.top$Term,levels = go.mf.top$Term)
  enrich.go <- rbind(go.bp.top,go.cc.top,go.mf.top)
  #修整描述的长度
  des_length = sapply(enrich.go$Term,  function(i){nchar(as.vector(i))})
  limit_length = ceiling(unname(quantile(des_length,0.8)))
  enrich.go$Term <- sapply(enrich.go$Term, function(i){
    if(limit_length < nchar(as.vector(i))){
      return(str_wrap(as.vector(i),width=60))
    }else{
      return(as.vector(i))
    }})
  enrich.go$Term <- factor(enrich.go$Term, level = enrich.go$Term)
  p <- ggplot(data = enrich.go, mapping = aes(x = Term, y = DEP, fill = Category)) +
    geom_bar(stat = 'identity') +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.box = "horizontal",
      legend.position = 'right',
      legend.text = element_text(size = 17),
      legend.title = element_text(size = 20),
      axis.title = element_text(color = "black", size = 20),
      axis.text.y = element_text(color = "black", size = 20),
      axis.text.x = element_text(angle = 80, hjust = 1, color = "black", size = 15),
      strip.text = element_text(size = 15)
    ) +
    #scale_fill_manual(
    # name = "GO category",
    # breaks = c("BP","CC","MF"),
    # values = c("BP" = "skyblue", "CC" = "lightgreen", "MF" = "#FD8D62"))+
    labs(x = 'GO Term', y = 'Number of proteins')
    #facet_grid(.~Category, scale = "free")
  if (!dir.exists(output)) {
    dir.create(output,recursive = T)
  }
  output<-paste0(output,"/",contrast)
  ggsave(paste0(output, '_GO_classify.png'), p, width  = 500,height = 300,units = 'mm', dpi = 300,bg = "white")
  ggsave(paste0(output, '_GO_classify.pdf'), p, width  = 500,height = 300,units = 'mm', dpi = 300)
  
}

dodgebarplot<-function(enrichplotdata,output,contrast){
  test<-enrichplotdata%>%group_by(Category)%>%
    arrange(Pvalue,.by_group = T)%>%
    top_n(n = -5,wt = Pvalue)%>%
    select(GO_ID,Category,Term,Up,Down,DEP)%>%
    tidyr::pivot_longer(cols = c(Up,Down),names_to = "variable",values_to = "Count")%>%
    mutate(y_description=str_wrap(Term,width=60))%>%
    arrange(DEP,.by_group = T)
  test$y_description<-factor(test$y_description,levels = unique(test$y_description))
  
  p<-ggplot(test,aes(x = y_description,y = Count,fill = variable))+
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),width = 0.75) +
    geom_text(aes(label = Count), position = position_dodge(width = 0.8), hjust = -0.8,size = 3, color = "black") +
    coord_flip() +  # 旋转坐标轴
    labs(x = paste0("GO Term"), y = "Number of proteins") +
    scale_fill_manual(values = c("Up" = "#f5780c", "Down" = "#0586c0")) +
    theme_classic(base_family = "sans") +
    theme(axis.title = element_text(color = "black", size = 13),
          axis.text.y = element_text(size = 10,color = "black"),
          axis.text.x = element_text(color = "black", size = 10),
          legend.title = element_blank(),
          panel.border = element_rect(color = "black",fill = NA))+
    scale_y_continuous(limits = c(0,max(test$Count)+0.5))+
    facet_grid(Category~.,scales = "free")+
    guides(fill=guide_legend(reverse = T))
  if (!dir.exists(output)) {
    dir.create(output,recursive = T)
  }
  output<-paste0(output,"/",contrast)
  ggsave(paste0(output, '_GO_updownbarplot.pdf'), p, width  = 12,height = 10, dpi = 300)
}

chordfun<-function(enrichfile,dep,enrich_class,outdir,contrast){
  process_data<-function(plotdata,plot_order,string_width,enrich_class){
    
    plotdata<-plotdata%>%arrange(plot_order)
    if(enrich_class=="GO"){
      colnames(plotdata)[1:3]<-c("ID","Category","Description")
    }else if(enrich_class=="KEGG"){
      colnames(plotdata)[1:3]<-c("ID","Description","Category")
    }else if(enrich_class=="InterPro"){
      plotdata<-plotdata%>%rename(DEP=Count,Pvalue=pvalue)
    }else if(enrich_class=="WikiPath"){
      plotdata<-plotdata%>%rename(DEP=Count,Pvalue=pvalue)
    }
    plotdata$Description<-str_wrap(plotdata$Description,width=string_width)
    return(plotdata)
  }
  data<-process_data(plotdata = enrichfile,plot_order = "Pvalue",string_width = 50,enrich_class = enrich_class)
  plotdata<-data[1:min(10,nrow(data)),]
  plotdata<-plotdata%>%select(ID,Description,Pvalue,Up_protein,Down_protein)%>%
    mutate(`-log10(Pvalue)`= -log10(Pvalue))
  
  if(enrich_class=="GO" || enrich_class=="KEGG" || enrich_class=="WikiPath"){
    plotdata$Up_protein<-plotdata$Up_protein%>%gsub("\\(.+?\\)",replacement="",x = .)
    plotdata$Down_protein<-plotdata$Down_protein%>%gsub("\\(.+?\\)",replacement="",x = .)
    plotdata$Down_protein<-plotdata$Down_protein%>%
      gsub("\\(",replacement = "",x = .)%>%gsub("\\)",replacement="",x = .)
    plotdata$Up_protein<-plotdata$Up_protein%>%
      gsub("\\(",replacement = "",x = .)%>%gsub("\\)",replacement="",x = .)
  }else{
    plotdata$Up_protein<-plotdata$Up_protein%>%gsub("\\(",replacement = "",x = .)%>%
      gsub("\\)",replacement = "",x=.)%>%gsub("-",replacement = "",x = .)
    plotdata$Down_protein<-plotdata$Down_protein%>%gsub("\\(",replacement = "",x = .)%>%
      gsub("\\)",replacement = "",x=.)%>%gsub("-",replacement = "",x = .)
  }
  
  plotdata$uniprotID<-paste(plotdata$Down_protein,plotdata$Up_protein,sep = ";")
  plotdata$uniprotID<-gsub("-","",plotdata$uniprotID)
  plotdata$uniprotID<-gsub(",",";",plotdata$uniprotID)
  plotdata<-plotdata%>%separate_rows(uniprotID,sep = ";")
  if("" %in% plotdata$uniprotID){
    plotdata<-plotdata[-which(plotdata$uniprotID==""),]
  }

  df1<-data.frame("from"=plotdata$ID,"to"=plotdata$uniprotID,"value"=1)
  color_fun<-function(log2FC){
    if(min(log2FC)>0){
      col_fun = colorRamp2(c(0, max(log2FC)), c("white", "red"))
      my_color<-col_fun(log2FC)
    }else if(max(log2FC)<0){
      col_fun = colorRamp2(c(min(log2FC),0), c("blue","white"))
      my_color<-col_fun(log2FC)
    }else{
      col_fun = colorRamp2(c(-max(abs(log2FC)),0,max(abs(log2FC))), c("blue","white","red"))
      my_color<-col_fun(log2FC)
    }
    return(my_color)
  }
  #log2FC
  
  dep<-dep%>%select(Protein.Accessions,Gene.Name,contains("Log2FC"))%>%
    arrange(desc(abs(select(.,contains("Log2FC")))))
  
  if(enrich_class=="GO"||enrich_class=="InterPro"||enrich_class=="WikiPath"){
    df2<-left_join(df1,dep,join_by(to==Protein.Accessions),relationship = "many-to-many")%>%unique()
  }else{
    if(any(duplicated(dep$Gene.Name))){
      dep<-dep[-which(duplicated(dep$Gene.Name)),]
    }
    df2<-left_join(df1,dep,join_by(to==Gene.Name),relationship = "many-to-many")%>%unique()
  }
  df2<-df2%>%mutate(ABS=abs(select(df2,contains("Log2FC"))))%>%
    group_by(from)%>%
    arrange(desc(ABS),.by_group = T)%>%
    do(head(., n = 10))%>%unique()
  colnames(df2)[c(2,5)]<-c("GeneName","log2FC")
  df2<-df2[order(df2$log2FC,decreasing = T),]
  gene_color<-color_fun(df2$log2FC)
  #df1<-df1[df1$to %in% df2$GeneName,]
  names(gene_color)<-df2$GeneName
  print(length(unique(df2$GeneName)))
  df3<-df2 %>% left_join(.,select(data,ID,Pvalue),by=join_by(from==ID))%>%arrange(Pvalue)
  col.list<-list("Funny"=c("#7C54FF","#FED840","#2D928F","#69DFDC","#FFA367","#BDB6E4","#5AA4C9","#98D9AD","#5490DE","#C574D5","#FF5C5C","#9777FF","#FEE067","#58A8A6","#86E3E1","#FFB686","#CAC5E9","#7CB6D4","#8BE5A8","#77A7E5","#D190DD","#FF7D7D","#89D7E1","#A780BF","#D2AC7F","#B8C05A","#72A974","#41ABA8","#55C5C2","#9BCBB5","#CDB78E","#E9A990","#D3AFBA","#9CB0DB","#7BAAD2","#6EB5BF","#83C7B6","#81C0BD","#6AA8CD","#7986DB","#9F7DD8","#D86CAC","#EB6484","#DC6592","#B96EC8","#B99ACC","#DBBD99","#C6CD7C","#8FBA91","#67BBB9","#76CFCD","#AED4C2","#D6C5A4","#EDBBA7","#DBC0C8","#B0C0E2","#96BBDB","#81C5C5","#86D5B6","#84D0BC","#7DBBD0","#959FE2","#B397DF","#E089BD","#EF839D","#D79B9E","#B0B9BF"))
  pathway_color<-col.list[["Funny"]][1:length(unique(df3$from))]
  names(pathway_color)<-unique(df3$from)
  grid.col<-c(pathway_color,gene_color)
  
  if(length(unique(df2$GeneName))<50){
    cex=0.9
  }else if(length(unique(df2$GeneName))<70){
    cex=0.8
  }else if(length(unique(df2$GeneName))<90){
    cex=0.7
  }else{
    cex=0.6
  }
  filename=paste0(outdir,"/",contrast,"_GO_chord.pdf")
  #绘图
  pdf(filename,width = 16,height = 10)
  par(mar = c(0, 0, 0, 1))
  circos.par(start.degree = 90,circle.margin=c(0.01,0.8,0.01,0.01),gap.after=1)
  chordDiagram(df3[,c(1,2)], annotationTrack = "grid", order = c(unique(df3$from),rev(df2$GeneName)),grid.col = grid.col,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df3))))))
  # we go back to the first track and customize sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = T, adj = c(0, 0.5),cex =cex)
  }, bg.border = NA) # here set bg.border to NA is important
  circos.clear()
  
  col_fun1<-function(log2FC){
    if(min(log2FC)>0){
      col_fun = colorRamp2(c(0, max(log2FC)), c("white", "red"))
    }else if(max(log2FC)<0){
      col_fun = colorRamp2(c(min(log2FC),0), c("blue","white"))
    }else{
      col_fun = colorRamp2(c(-max(abs(log2FC)),0,max(abs(log2FC))), c("blue","white","red"))
    }
    return(col_fun)
  }
  
  at_break<-function(log2FC){
    if(min(log2FC)>0){
      at<-c(0,round(min(log2FC),1),round(max(log2FC),1))
    }else if(max(log2FC)<0){
      at<-c(round(min(log2FC),1),round(max(log2FC),1),0)
    }else{
      at<-c(round(min(log2FC),1),round(min(log2FC)/2,1),0,round(max(log2FC)/2,1),round(max(log2FC),1))
    }
    return(at)
  }
  
  log2_leged<-Legend(col_fun = col_fun1(df2$log2FC),
                     title = "Log2(Fold change)", direction = "horizontal")
  descipt<-Legend(labels = str_wrap(string = paste(names(pathway_color),data[which(names(pathway_color)%in%data$ID),c("Description")]),width = 50), title = "Description", legend_gp = gpar(fill = pathway_color))
  legendlist<-packLegend(log2_leged, descipt)
  draw(legendlist,x = unit(14, "in"), y = unit(5, "in"))
  
  dev.off()
}



main<-function(go.data,gene.data,output,contrast){
  GObarplot(go.data,output,contrast)
  dodgebarplot(go.data,output,contrast=contrast)
  chordfun(enrichfile=go.data,dep=gene.data,enrich_class="GO",outdir=output,contrast=contrast)
}
sapply(contrast_data[,1],function(x) {
  print(x)
  go_file<-paste0(go_dir,"/",x,"/",x,"_GO_enrichment.xlsx")
  gene_file<-paste0(gene_dir,"/",x,"-DEP_results.xlsx")
  output<-paste0(outdir,"/",x)
  if(!dir.exists(output)){
    dir.create(output,recursive = T)
  }
  go.data<-read.xlsx(go_file)
  go.data$Term<-str_wrap(go.data$Term,width=50)
  gene.data<-read.xlsx(gene_file)
  main(go.data,gene.data,output,contrast=x)
})
