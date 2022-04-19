renv::install("tidyverse")
renv::install("svglite")

library(tidyverse)
library(svglite)
gene_set_enrich_tbl_file<-"./data/trans_res_hub_GS_tbl/GM12878_trans_res_hub_entrez_GOBP_enrich_tbl.Rda"

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

ggsave("~/Documents/multires_bhicect/weeklies/group_meeting/group_meeting_04_2022/img/GM12878_hub_tres_GOBP_lilopop.svg",width = 30,height=15,units = "cm")
