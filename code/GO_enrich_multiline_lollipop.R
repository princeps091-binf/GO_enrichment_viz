library(tidyverse)
library(svglite)
library(rjson)
library(simplifyEnrichment)
library(cowplot)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
get_obj_in_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#-------------------------------------------------------------------
cell_line_set<-c("GM12878","HMEC","H1")
result <- fromJSON(file = "./data/MSigDB_tbl/c5.go.bp.v7.5.1.json")
GO_ID_map_vec<-map_chr(result,function(x){
  x$exactSource
})

GO_summary_l<-map(cell_line_set,function(cell_line){
  message(cell_line)
  gene_set_enrich_tbl_file<-paste0("./data/trans_res_hub_GS_tbl/",cell_line,"_trans_res_top_hub_entrez_GOBP_enrich_tbl.Rda")
  gene_set_enrich_tbl<-get_obj_in_fn(gene_set_enrich_tbl_file)
  
  cand_GO<-gene_set_enrich_tbl %>% 
    filter(FDR<=0.01) %>%
    dplyr::select(Gene.Set) %>% 
    unlist
  cand_GO_ID<-GO_ID_map_vec[cand_GO]
  
  mat <- GO_similarity(cand_GO_ID,measure="Resnik",ont = "BP")
  df <- simplifyGO(mat)
  
  semantic_cl<-df %>% 
    as_tibble() %>% 
    group_by(cluster) %>% 
    summarise(n=n()) %>% 
    filter(n>3) %>% 
    arrange(desc(n)) %>% 
    dplyr::select(cluster) %>% 
    unlist
  
  GO_cluster<-df$cluster
  names(GO_cluster)<-df$id
  
  gene_set_enrich_candidate_tbl<-gene_set_enrich_tbl %>% 
    mutate(GO.ID=GO_ID_map_vec[Gene.Set]) %>% 
    filter(FDR<=0.01) %>%
    mutate(GO.cluster=GO_cluster[GO.ID])
  
  out_tbl<-gene_set_enrich_candidate_tbl %>% 
    filter(GO.cluster %in% semantic_cl) %>% 
    group_by(GO.cluster) %>% 
    slice_min(FDR,n = 5) %>% 
    mutate(GO.Set=str_wrap(unlist(lapply(str_split(Gene.Set,"_"),function(x)paste(x[-1],collapse = " "))),30)) %>% 
    mutate(GO.Set=fct_reorder(GO.Set,-log10(FDR)),
           cell.line=cell_line)
  return(out_tbl)
})

do.call(bind_rows,GO_summary_l) %>% 
  ggplot(.,aes(-log10(FDR),GO.Set,size=OR))+
  geom_point()+theme_classic()+
  theme(text=element_text(size=18))+
  facet_wrap(cell.line~.,scales='free')
ggsave("~/Documents/multires_bhicect/weeklies/weekly61/img/GO_tables.png",units = "cm",width = 60, height= 30)
