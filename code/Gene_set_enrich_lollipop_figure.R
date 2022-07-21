renv::install("tidyverse")
renv::install("svglite")

library(tidyverse)
library(svglite)
library(rjson)

gene_set_enrich_tbl_file<-"./data/trans_res_hub_GS_tbl/GM12878_trans_res_top_hub_entrez_GOBP_enrich_tbl.Rda"

gene_set_enrich_tbl<-get(load(gene_set_enrich_tbl_file))
tmp_obj<-names(mget(load(gene_set_enrich_tbl_file)))
rm(list=tmp_obj)
rm(tmp_obj)

gg_pop<-gene_set_enrich_tbl %>% 
  filter(FDR<=0.01) %>% 
  arrange(FDR) %>% 
  slice_head(n=15) %>% 
  mutate(Gene.Set=fct_reorder(Gene.Set,-log10(FDR))) %>% 
  ggplot(.,aes(-log10(FDR),Gene.Set,size=OR))+
  geom_point()+theme_minimal()
gg_pop

ggsave("~/Documents/multires_bhicect/Poster/img/F3/HMEC_hub_tres_GOBP_lollipop.svg",width = 30,height=15,units = "cm")


result <- fromJSON(file = "./data/MSigDB_tbl/c5.go.bp.v7.5.1.json")
GO_ID_map_vec<-map_chr(result,function(x){
  x$exactSource
})

cand_GO<-gene_set_enrich_tbl %>% 
  filter(FDR<=0.01) %>%
  dplyr::select(Gene.Set) %>% 
  unlist
cand_GO_ID<-GO_ID_map_vec[cand_GO]
library(simplifyEnrichment)

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

gene_set_enrich_candidate_tbl %>% 
  filter(GO.cluster %in% semantic_cl) %>% 
  group_by(GO.cluster) %>% 
  slice_min(FDR,n = 5) %>% 
  mutate(GO.Set=str_wrap(unlist(lapply(str_split(Gene.Set,"_"),function(x)paste(x[-1],collapse = " "))),50)) %>% 
  mutate(GO.Set=fct_reorder(GO.Set,-log10(FDR))) %>% 
  ggplot(.,aes(-log10(FDR),GO.Set,size=OR))+
  geom_point()+theme_minimal()

