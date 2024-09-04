rm(list = ls())
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggraph)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(tidygraph)))
suppressWarnings(suppressMessages(library(igraph)))
suppressWarnings(suppressMessages(library(showtext)))
source("config/9_PPI.conf")
source("config/default.conf")
source("config/set_Theme.R")
AddFont("config/myFonts")
showtext_auto()
contrast_data <- read.xlsx(contrast.file)
Protein_main<-function(net_file,atrri_file,contrast,outdir){
  set.seed(1234)
  net<-read.xlsx(net_file)
  atrri<-read.xlsx(atrri_file)%>%arrange(desc(Degree))%>%head(100)
  test<-net%>%filter(Score>Score_number)%>%
    right_join(.,select(.data=atrri,Protein,Regulation,Degree),by=join_by(Node1==Protein))%>%
    filter(Node2%in%atrri$Protein)
  graph_tbl<-test[,c(1,2,4)]
  colnames(graph_tbl)<-c("from","to","type")
  graph <- as_tbl_graph(graph_tbl)
  vertex_attr(graph)$type<-left_join(x = data.frame("name"=vertex_attr(graph)$name),
                                     y = graph_tbl[,c(1,3)],by=join_by(name==from))%>%unique()%>%{.[,2]}
  vertex_attr(graph)$Degree<-left_join(x = data.frame("name"=vertex_attr(graph)$name),
                                       y = atrri[,c(4,7)],by=join_by(name==Protein))%>%unique()%>%{.[,2]}
  step<-max(vertex_attr(graph)$Degree)%/%4
  p1<-ggraph(graph = graph,layout = "nicely")+
    geom_edge_link(edge_alpha = 0.2,color="gray")+
    geom_node_point(aes(fill=type,size=Degree),color = "#808080",shape=21)+
    geom_node_text(aes(label=name),repel = T,check_overlap = T,max.overlaps=5)+
    scale_fill_manual(values = c("up"="#f5780c","down"="#0586c0"))+
    scale_size(range = c(3,8),breaks = c(step,step*2,step*3))+
    theme_graph(base_family = "sans")+
    guides(fill=guide_legend(reverse = T,title = "Regulation",order=2,override.aes = list(size=5)),
           size=guide_legend(override.aes = list(fill="black")))
  write.xlsx(x=test,file = paste0(outdir,"/",contrast,"_degree_top_pro.network.xlsx"))
  filename<-paste0(outdir,"/",contrast,".network.pro.")
  ggsave(filename = paste0(filename,"pdf"),plot = p1,width = 10,height = 8)
  ggsave(filename = paste0(filename,"png"),plot = p1,width = 10,height = 8,units = "in",dpi = 300)
}

gene_main<-function(net_file,atrri_file,contrast,outdir){
  set.seed(1234)
  net<-read.xlsx(net_file)
  atrri<-read.xlsx(atrri_file)%>%arrange(desc(Degree))%>%head(100)
  test<-net%>%filter(Score>Score_number)%>%
    right_join(.,select(.data=atrri,GeneName,Regulation,Degree),by=join_by(Node1==GeneName))%>%
    filter(Node2%in%atrri$GeneName)
  graph_tbl<-test[,c(1,2,4)]
  colnames(graph_tbl)<-c("from","to","type")
  graph <- as_tbl_graph(graph_tbl)
  vertex_attr(graph)$type<-left_join(x = data.frame("name"=vertex_attr(graph)$name),
                                     y = graph_tbl[,c(1,3)],by=join_by(name==from))%>%unique()%>%{.[,2]}
  vertex_attr(graph)$Degree<-left_join(x = data.frame("name"=vertex_attr(graph)$name),
                                       y = atrri[,c(1,7)],by=join_by(name==GeneName))%>%unique()%>%{.[,2]}
  step<-max(vertex_attr(graph)$Degree)%/%4
  p1<-ggraph(graph = graph,layout = "nicely")+
    geom_edge_link(edge_alpha = 0.2,color="gray")+
    geom_node_point(aes(fill=type,size=Degree),color = "#808080",shape=21)+
    geom_node_text(aes(label=name),repel = T,check_overlap = T,max.overlaps=5)+
    scale_fill_manual(values = c("up"="#f5780c","down"="#0586c0"))+
    scale_size(range = c(3,8),breaks = c(step,step*2,step*3))+
    theme_graph(base_family = "sans")+
    guides(fill=guide_legend(reverse = T,title = "Regulation",order=2,override.aes = list(size=5)),
           size=guide_legend(override.aes = list(fill="black")))
  write.xlsx(x=test,file = paste0(outdir,"/",contrast,"_degree_top_gene.network.xlsx"))
  filename<-paste0(outdir,"/",contrast,".network.gene.")
  ggsave(filename = paste0(filename,"pdf"),plot = p1,width = 10,height = 8)
  ggsave(filename = paste0(filename,"png"),plot = p1,width = 10,height = 8,units = "in",dpi = 300)
}
sapply(contrast_data$contrast,function(x){
  print(x)
  outdir<-paste0("result/9_PPI/",x)
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  netfile_pro<-paste0("result/9_PPI/",x,"/",x,"_network_pro.xlsx")
  netfile_gene<-paste0("result/9_PPI/",x,"/",x,"_network_gene.xlsx")
  atrri<-paste0("result/9_PPI/",x,"/",x,"_attributes.xlsx")
  Protein_main(net_file=netfile_pro,atrri_file=atrri,contrast=x,outdir=outdir)
  gene_main(net_file=netfile_gene,atrri_file=atrri,contrast=x,outdir=outdir)
})


