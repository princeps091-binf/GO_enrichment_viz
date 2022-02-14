# Examine intersection of Gene sets
library(tidyverse)
library(UpSetR)

base::load("./data/GOBP_gene_set_l.Rda")

elements <- unique(unlist(GOBP_set))

data <- do.call(cbind,lapply(GOBP_set, function(x) {
  return(as.numeric(elements %in% x))
}))
data <- data[which(rowSums(data) != 0), ]

gene_set_edge_tbl<-do.call(bind_rows,apply(data,1,function(x){
  
  tmp_set<-names(which(x > 0))
  if(length(tmp_set)<2){
    return(tibble(V1=tmp_set,V2=tmp_set))
  }else{
    return(as_tibble(t(combn(tmp_set,2))))
  }
}))

gene_set_edge_tbl<-gene_set_edge_tbl %>% 
  group_by(V1,V2) %>% 
  summarise(weight=n())

library(igraph)
components(graph_from_data_frame(gene_set_edge_tbl))
