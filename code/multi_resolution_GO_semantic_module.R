library(igraph)
library(tidyverse)
library(rvest)
library(rrvgo)
library(org.Hs.eg.db)
library(seriation)
library(viridis)
library(ggridges)
library(UpSetR)
#----------------------

get_GO_ID_fn<-function(x){
  test_term<-read_html(x)
  tbl_l<-test_term%>%html_elements("tr")
  tmp_chr<-as.character(tbl_l[grep("^<tr>\n<th>Exact source</th>\n",test_term%>%html_elements("tr"))])
  if(length(tmp_chr)<1){return(NA)}else {
    return(tmp_chr %>% str_extract(.,"GO:[0-9]+"))
  }}
tbl_in_fn<-function(file){
  tmp_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(tmp_tbl)
}
#----------------------
hm_gene_set<-as_tibble(read.table("~/Documents/multires_bhicect/data/epi_data/Gene_annotation/c5.all.v7.3.entrez.gmt",header = F,sep = "\t",fill=T))

cell_line_files<-grep("HMEC_CAGE",list.files("./data/"),value = T)
hub_files<-grep("hub",cell_line_files,value=T)
cell_line_GO_enrich_set_vec<-unique(unlist(lapply(hub_files,function(file){
  gene_set_enrich_tbl<-get(load(paste0("./data/",file)))
  tmp_obj<-names(mget(load(paste0("./data/",file))))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(gene_set_enrich_tbl %>% filter(FDR<=0.01) %>% distinct(Gene.Set) %>% unlist)
  
})))

cell_line_GO_tbl<-do.call(bind_rows,lapply(cell_line_files,function(file){
  tmp_line<-strsplit(file,split="_")[[1]][2]
  gene_set_enrich_tbl<-get(load(paste0("./data/",file)))
  tmp_obj<-names(mget(load(paste0("./data/",file))))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(gene_set_enrich_tbl %>% filter(Gene.Set %in% cell_line_GO_enrich_set_vec) %>% arrange(FDR) %>% mutate(line=tmp_line))
  
}))

cell_line_GO_tbl %>% ggplot(.,aes(line,Gene.Set,fill=-log10(FDR)))+
  geom_tile()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  scale_fill_viridis_c()


cl_set_combo_tbl<-cell_line_GO_tbl %>% distinct(Gene.Set) %>% 
  left_join(.,hm_gene_set%>%dplyr::select(V1,V2),by=c("Gene.Set" = "V1")) %>% 
  mutate(GO.ID=map_chr(V2,function(x){
    unlist(get_GO_ID_fn(x))
  }))%>% 
  filter(!(is.na(GO.ID)))

simMatrix <- calculateSimMatrix(cl_set_combo_tbl$GO.ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Resnik")
d<-as.dist(1/(simMatrix+1e-3))
order <- seriate(d,method = "OLO")
image(simMatrix[get_order(order),get_order(order)])

main_sub_g<-graph_from_adjacency_matrix(simMatrix,mode = "undirected",weighted = T)
louvain_sample_cluster<-cluster_louvain(main_sub_g)

sample_comm<-unique(louvain_sample_cluster$membership)
comm_edge_tbl<-do.call(bind_rows,lapply(sample_comm,function(i){
  tmp_v<-V(main_sub_g)$name[which(cluster_louvain(main_sub_g)$membership==i)]
  expand_grid(ego=tmp_v,alter=tmp_v) %>% mutate(x=i)
}))
tmp_mat<-matrix(0,nrow = nrow(simMatrix),ncol=ncol(simMatrix),dimnames = dimnames(simMatrix))
tmp_mat[as.matrix(comm_edge_tbl[,1:2])]<-comm_edge_tbl$x
image(tmp_mat[get_order(order),get_order(order)],col=plasma(9))

module_GO_set<-cl_set_combo_tbl %>% 
  inner_join(comm_edge_tbl %>% filter(x==5) %>% 
               summarise(GO.ID=unique(c(ego,alter)))) %>% distinct(Gene.Set) %>% unlist

cell_line_GO_tbl %>% 
  filter(Gene.Set%in%module_GO_set) %>% 
  ggplot(.,aes(line,Gene.Set,fill=log10(OR)))+
  geom_tile()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  scale_fill_viridis_c()

cell_line_GO_tbl %>% 
  filter(Gene.Set%in%module_GO_set) %>%
  filter(line=="5kb") %>% 
  arrange(desc(OR))
