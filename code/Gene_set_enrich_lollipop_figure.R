renv::install("tidyverse")
renv::install("svglite")

library(tidyverse)
library(svglite)
gene_set_enrich_tbl_file<-"./data/hub_mres_HMEC_CAGE_rich_GOBP_enrich_tbl.Rda"

gene_set_enrich_tbl<-get(load(gene_set_enrich_tbl_file))
tmp_obj<-names(mget(load(gene_set_enrich_tbl_file)))
rm(list=tmp_obj)
rm(tmp_obj)

gg_pop<-gene_set_enrich_tbl %>% 
  filter(FDR<=0.01) %>% 
  arrange(FDR) %>% 
  slice_head(n=30) %>% 
  mutate(Gene.Set=fct_reorder(Gene.Set,-log10(FDR))) %>% 
  ggplot(.,aes(-log10(FDR),Gene.Set,size=OR))+
  geom_point()+theme_minimal()
ggsave("~/Documents/multires_bhicect/weeklies/group_meeting/group_meeting_12_2021/img/GM12878_hub_mres_GOBP_lilopop.svg")
