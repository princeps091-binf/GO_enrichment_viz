# Multi-sets semantic module evaluation

library(tidyverse)
library(rvest)
library(rrvgo)
library(org.Hs.eg.db)
library(seriation)
library(viridis)
library(ggridges)
library(UpSetR)
cell_line_files<-grep("hub_5kb_",list.files("./data/"),value = T)

gene_set_enrich_tbl_file<-"./data/hub_5kb_GM12878_CAGE_GOBP_enrich_tbl.Rda"

cell_line_GO_tbl<-do.call(bind_rows,lapply(cell_line_files,function(file){
  tmp_line<-strsplit(file,split="_")[[1]][3]
  gene_set_enrich_tbl<-get(load(paste0("./data/",file)))
  tmp_obj<-names(mget(load(gene_set_enrich_tbl_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(gene_set_enrich_tbl %>% filter(FDR<=0.01) %>% arrange(FDR) %>% mutate(line=tmp_line))
  
}))

GO_line_tbl<-cell_line_GO_tbl %>% group_by(line) %>% 
  summarise(GO.set=list(unique(Gene.Set)))
GO_line_l<-GO_line_tbl$GO.set
names(GO_line_l)<- GO_line_tbl$line

upset(fromList(GO_line_l),order.by = "freq")

inter_set_tbl<-as_tibble(fromList(GO_line_l)) %>% mutate(GO.Set=unique(unlist(GO_line_l)))

cell_line_GO_tbl %>% 
  left_join(.,
            inter_set_tbl %>% 
              filter(HMEC == 0 & H1 == 0 & GM12878 == 1) %>% 
              dplyr::select(GO.Set) %>% mutate(inter.set="ALL"),by=c("Gene.Set"="GO.Set")
  ) %>% 
filter(!(is.na(inter.set))) #%>% 
ggplot(.,aes(-log10(FDR),color=inter.set)) + geom_density()+facet_wrap(inter.set~.,scales="free")


cl_shared_set_combo_tbl<-inter_set_tbl %>% 
  filter(HMEC == 1 & H1 == 1 & GM12878 == 1) %>% 
  dplyr::select(GO.Set)

cl_set_combo_tbl<-cell_line_GO_tbl %>% 
  distinct(Gene.Set)

hm_gene_set<-as_tibble(read.table("~/Documents/multires_bhicect/data/epi_data/Gene_annotation/c5.all.v7.3.entrez.gmt",header = F,sep = "\t",fill=T))

cl_set_combo_tbl<-cl_set_combo_tbl %>%   left_join(.,hm_gene_set%>%dplyr::select(V1,V2),by=c("Gene.Set" = "V1"))



get_GO_ID_fn<-function(x){
  test_term<-read_html(x)
  tbl_l<-test_term%>%html_elements("tr")
  tmp_chr<-as.character(tbl_l[grep("^<tr>\n<th>Exact source</th>\n",test_term%>%html_elements("tr"))])
  if(length(tmp_chr)<1){return(NA)}else {
    return(tmp_chr %>% str_extract(.,"GO:[0-9]+"))
  }}

cl_shared_set_combo_tbl<-cl_shared_set_combo_tbl %>% 
  left_join(.,hm_gene_set%>%dplyr::select(V1,V2),by=c("GO.Set" = "V1")) %>% 
#  dplyr::slice(1:10) %>%  
  mutate(GO.ID=map_chr(V2,function(x){
  unlist(get_GO_ID_fn(x))
}))
cl_shared_set_combo_tbl<-cl_shared_set_combo_tbl %>% 
  filter(is.na(GO.ID)) %>% 
  dplyr::select(V2)


cl_set_combo_tbl<-cl_set_combo_tbl%>%
  mutate(GO.ID=unlist(tmp_tbl_l)) %>% 
  filter(!(is.na(GO.ID))) 

simMatrix <- calculateSimMatrix(cl_set_combo_tbl$GO.ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Resnik")
d<-as.dist(1/(simMatrix+1e-3))
order <- seriate(d,method = "OLO")
image(simMatrix[get_order(order),get_order(order)])

library(igraph)
main_sub_g<-graph_from_adjacency_matrix(simMatrix,mode = "undirected",weighted = T)
louvain_sample_cluster<-cluster_louvain(main_sub_g)

sample_comm<-unique(louvain_sample_cluster$membership)
comm_edge_tbl<-do.call(bind_rows,lapply(sample_comm,function(i){
  tmp_v<-V(main_sub_g)$name[which(cluster_louvain(main_sub_g)$membership==i)]
  expand_grid(ego=tmp_v,alter=tmp_v) %>% mutate(x=i)
}))
tmp_mat<-matrix(0,nrow = nrow(simMatrix),ncol=ncol(simMatrix),dimnames = dimnames(simMatrix))
tmp_mat[as.matrix(comm_edge_tbl[,1:2])]<-comm_edge_tbl$x
image(tmp_mat[get_order(order),get_order(order)],col=plasma(length(sample_comm)+1))
