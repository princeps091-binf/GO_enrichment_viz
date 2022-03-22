library(tidyverse)
library(rvest)
library(rrvgo)
library(org.Hs.eg.db)
library(seriation)
library(viridis)
library(ggridges)
#------------------------------
get_obj_in_fn<-function(file){
  out_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#------------------------------
gene_set_enrich_tbl_file<-"./data/HMEC_TAD_GOBP_enrich_tbl.Rda"
gene_set_enrich_tbl<-get_obj_in_fn(gene_set_enrich_tbl_file)


cl_set_combo_tbl<-gene_set_enrich_tbl %>% 
  filter(FDR<=0.01)

hm_gene_set<-as_tibble(read.table("~/Documents/multires_bhicect/data/epi_data/Gene_annotation/c5.all.v7.3.entrez.gmt",header = F,sep = "\t",fill=T))

cl_set_combo_tbl<-cl_set_combo_tbl %>%   
  left_join(.,hm_gene_set%>%dplyr::select(V1,V2),by=c("Gene.Set" = "V1"))



tmp_tbl_l<-lapply(cl_set_combo_tbl$V2,function(x){
  test_term<-read_html(x)
  tbl_l<-test_term%>%html_elements("tr")
  tmp_chr<-as.character(tbl_l[grep("^<tr>\n<th>Exact source</th>\n",test_term%>%html_elements("tr"))])
  if(length(tmp_chr)<1){return(NA)}else {
    return(tmp_chr %>% str_extract(.,"GO:[0-9]+"))
  }
})
cl_set_combo_tbl<-cl_set_combo_tbl%>%
  mutate(GO.ID=unlist(tmp_tbl_l)) %>% 
  filter(!(is.na(GO.ID))) 

simMatrix <- calculateSimMatrix(cl_set_combo_tbl$GO.ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Resnik")

library(seriation)
d<-as.dist(1/(simMatrix+1e-3))
order <- seriate(d,method = "HC")
image(simMatrix[get_order(order),get_order(order)])

library(igraph)
main_sub_g<-graph_from_adjacency_matrix(simMatrix,mode = "undirected",weighted = T,diag = F)
louvain_sample_cluster<-cluster_louvain(main_sub_g)

sample_comm<-unique(louvain_sample_cluster$membership)

comm_edge_tbl<-do.call(bind_rows,lapply(sample_comm,function(i){
  tmp_v<-V(main_sub_g)$name[which(louvain_sample_cluster$membership==i)]
  expand_grid(ego=tmp_v,alter=tmp_v) %>% mutate(x=i)
}))
tmp_mat<-matrix(0,nrow = nrow(simMatrix),ncol=ncol(simMatrix),dimnames = dimnames(simMatrix))
tmp_mat[as.matrix(comm_edge_tbl[,1:2])]<-comm_edge_tbl$x
image(tmp_mat[get_order(order),get_order(order)],col=plasma(length(sample_comm)+1))

sim_tbl<-do.call(bind_rows,lapply(sample_comm,function(i){
  tmp_w<-E(induced_subgraph(main_sub_g,louvain_sample_cluster$names[which(louvain_sample_cluster$membership==i)]))$weight
  return(tibble(GO.set=i,w=tmp_w))
}))


print(cl_set_combo_tbl %>% 
  inner_join(comm_edge_tbl %>% filter(x==4) %>% 
  summarise(GO.ID=unique(c(ego,alter)))) %>% arrange(FDR),n = 100)

sim_tbl %>% 
  ggplot(.,aes(x=w,y=as.factor(GO.set),fill=as.factor(GO.set)))+
  geom_density_ridges(alpha=0.8)+
  scale_fill_manual(breaks = c(0,sample_comm),values=plasma(9))

cl_set_combo_tbl %>% 
  inner_join(comm_edge_tbl %>% group_by(x) %>% 
               summarise(GO.ID=unique(c(ego,alter)))) %>% 
  dplyr::rename(set=x) %>% 
  ggplot(.,aes(x=-log10(FDR),y=as.factor(set)))+
  geom_density_ridges()

cl_set_combo_tbl %>% 
  inner_join(comm_edge_tbl %>% group_by(x) %>% 
               summarise(GO.ID=unique(c(ego,alter)))) %>% 
  dplyr::rename(set=x) %>% 
  ggplot(.,aes(x=OR,y=as.factor(set)))+
  geom_density_ridges()
#---------------------------------------------------
# Examine overlap in gene content between semantic blocks
semantic_module_content_tbl<-comm_edge_tbl %>% group_by(x) %>% 
  summarise(GO.ID=list(unique(c(ego,alter)))) %>% 
  dplyr::rename(sem.mod=x)

semantic_module_content_tbl<-semantic_module_content_tbl %>% 
  mutate(GO.tbl=map(GO.ID,function(x){
    return(cl_set_combo_tbl %>% filter(GO.ID %in% x) %>% dplyr::select(GO.ID,Gene.Set,FDR,OR))
  })) %>% 
  mutate(entrez.content=map(GO.ID,function(x){
    tmp_names<-cl_set_combo_tbl %>% filter(GO.ID %in% x) %>% dplyr::select(Gene.Set) %>% unlist
    return(hm_gene_set %>% 
      filter(V1 %in% tmp_names) %>% 
      dplyr::select(-V2) %>% 
      pivot_longer(!V1, names_to="gene.id",values_to = "entrez.id") %>% 
      dplyr::select(-gene.id) %>% 
      distinct(entrez.id) %>% unlist %>% as.character)
    
  }))

save(semantic_module_content_tbl,file="./data/semantic_module_tbl/HMEC_semantic_module_tbl.Rda")

sem_mod_gene_l<-lapply(semantic_module_content_tbl$GO.ID,function(i){
  tmp_names<-cl_set_combo_tbl %>% filter(GO.ID %in% i) %>% dplyr::select(Gene.Set) %>% unlist
  hm_gene_set %>% 
    filter(V1 %in% tmp_names) %>% 
    dplyr::select(-V2) %>% 
    pivot_longer(!V1, names_to="gene.id",values_to = "entrez.id") %>% 
    dplyr::select(-gene.id) %>% 
    distinct(entrez.id) %>% unlist %>% as.character
  
  
})

library(UpSetR)
names(sem_mod_gene_l)<-paste0("module.",1:length(sem_mod_gene_l))
upset(fromList(sem_mod_gene_l),order.by = "freq",nsets = length(sample_comm))
