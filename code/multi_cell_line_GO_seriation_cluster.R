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

#----------------------

cell_line_files<-grep("hub_5kb_",list.files("./data/"),value = T)
hm_gene_set<-as_tibble(read.table("~/Documents/multires_bhicect/data/epi_data/Gene_annotation/c5.all.v7.3.entrez.gmt",header = F,sep = "\t",fill=T))


cell_line_GO_enrich_set_vec<-unique(unlist(lapply(cell_line_files,function(file){
  gene_set_enrich_tbl<-get(load(paste0("./data/",file)))
  tmp_obj<-names(mget(load(paste0("./data/",file))))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(gene_set_enrich_tbl %>% filter(FDR<=0.01) %>% distinct(Gene.Set) %>% unlist)
  
})))

cell_line_GO_tbl<-do.call(bind_rows,lapply(cell_line_files,function(file){
  tmp_line<-strsplit(file,split="_")[[1]][3]
  gene_set_enrich_tbl<-get(load(paste0("./data/",file)))
  tmp_obj<-names(mget(load(paste0("./data/",file))))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(gene_set_enrich_tbl %>% filter(Gene.Set %in% cell_line_GO_enrich_set_vec) %>% arrange(FDR) %>% mutate(line=tmp_line))
  
}))

cell_line_GO_tbl %>% ggplot(.,aes(line,Gene.Set,fill=log10(OR)))+
  geom_tile()+
  theme(axis.text.y = element_blank(),
  axis.ticks.y = element_blank())+
  scale_fill_viridis_c()
GO_mat<-cell_line_GO_tbl %>% dplyr::select(Gene.Set,FDR,line) %>% 
  mutate(FDR=-log10(FDR)) %>% 
  pivot_wider(names_from=line,values_from=FDR)
x<-GO_mat %>% dplyr::select(-Gene.Set) %>% as.matrix
d <- dist(x)
o <- seriate(d,method="OLO")
cell_line_GO_tbl %>% 
  left_join(.,tibble(Gene.Set=GO_mat$Gene.Set[get_order(o)],rank=1:nrow(x))) %>% 
  ggplot(.,aes(line,rank,fill=-log10(FDR)))+
  geom_tile()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  scale_fill_viridis_c()
pca_try<-svd(x)
plot(pca_try$u[,c(1,3)])
