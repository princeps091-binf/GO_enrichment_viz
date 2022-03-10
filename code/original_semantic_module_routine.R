#Explore how to unpack hierchical merges of community detection
library(tidyverse)
library(rvest)
library(rrvgo)
library(org.Hs.eg.db)
library(seriation)
library(viridis)
library(ggridges)
library(data.tree)
#------------------------------
get_obj_in_fn<-function(file){
  out_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#------------------------------
gene_set_enrich_tbl_file<-"./data/H1_5kb_tss_compound_hub_GOBP_enrich_tbl.Rda"
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
sample_greedy_cluster<-cluster_fast_greedy(main_sub_g)

tmp_dendro<-as.dendrogram(sample_greedy_cluster)
GO_semantic_merge<-as.Node(tmp_dendro)
nleaf<-length(GO_semantic_merge$Get('name',filterFun = isLeaf))
GO_semantic_merge$Set(name = 1:(GO_semantic_merge$totalCount-nleaf),filterFun = function(x)!(isLeaf(x)))
node_ancestor<-GO_semantic_merge$Get(function(x){x$Get('name',traversal='ancestor')})
node_set<-names(node_ancestor)
node_set<-node_set[!grepl("GO:",node_set)]
node_children_tbl<-do.call(bind_rows,lapply(node_set,function(x){
  tmp_ch<-names(which(unlist(lapply(node_ancestor,function(y) x %in% y))))
  return(tibble(parent.hub=x,children.hub=tmp_ch))
}))
node_GO_children_tbl<-node_children_tbl %>% 
  filter(grepl("GO:",children.hub))

GO_semantic_cl_tbl<-node_GO_children_tbl %>% 
  group_by(parent.hub) %>% 
  summarise(GO.Set=list(children.hub))

candidate_module_tbl<-GO_semantic_cl_tbl %>% 
  mutate(resnik.stat=map_dbl(GO.Set,function(x){
    sum(apply(simMatrix[x,x],1,function(y)(sum(y>0)/length(y))>0.5))/length(x)
  })) %>% filter(resnik.stat==1)

top_module_set<-candidate_module_tbl$parent.hub[which(unlist(lapply(candidate_module_tbl$parent.hub,function(x){
  
  !(any(candidate_module_tbl$parent.hub %in% node_ancestor[[x]][-1])) 

})))]

top_module_tbl<-candidate_module_tbl %>% 
  filter(parent.hub %in% top_module_set)


sample_comm<-1:nrow(top_module_tbl)

comm_edge_tbl<-do.call(bind_rows,lapply(sample_comm,function(i){
  tmp_v<-top_module_tbl %>% dplyr::slice(i) %>% unnest(GO.Set) %>% dplyr::select(GO.Set) %>% unlist
  expand_grid(ego=tmp_v,alter=tmp_v) %>% mutate(x=i)
}))
tmp_mat<-matrix(0,nrow = nrow(simMatrix),ncol=ncol(simMatrix),dimnames = dimnames(simMatrix))
tmp_mat[as.matrix(comm_edge_tbl[,1:2])]<-comm_edge_tbl$x
image(tmp_mat[get_order(order),get_order(order)],col=plasma(length(sample_comm)+1))

print(cl_set_combo_tbl %>% 
        inner_join(comm_edge_tbl %>% filter(x==1) %>% 
                     summarise(GO.ID=unique(c(ego,alter)))) %>% arrange(FDR),n = 100)
