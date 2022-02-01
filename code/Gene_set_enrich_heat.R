renv::install("tidyverse")
renv::install("rvest")
renv::install("bioc::rrvgo")
renv::install("bioc::org.Hs.eg.db")
renv::install("seriation")

library(tidyverse)
library(rrvgo)
library(rvest)
library(org.Hs.eg.db)
library(seriation)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
tbl_folder<-"./data/"
foi<-c("hub_mres_HMEC_CAGE_rich_GOBP_enrich_tbl.Rda","hub_50kb_HMEC_CAGE_rich_GOBP_enrich_tbl.Rda","hub_5kb_HMEC_CAGE_rich_GOBP_enrich_tbl.Rda","HMEC_CAGE_active_GOBP_enrich_tbl.Rda")

cl_set_combo_tbl<-do.call(bind_rows,lapply (foi, function(cl_file){
  #cl_set<-unlist(lapply(strsplit(cl_file,split="_"),'[',1))
  #cl_set<-paste(unlist(lapply(strsplit(cl_file,split="_"),'[',c(2,3))),collapse = "_")
  
  path_tbl<-get(load(paste0(tbl_folder,cl_file)))
  #load(cl_file)
  return(path_tbl %>%arrange(desc(OR))%>%filter(FDR<=0.01)%>%slice_head(n=15))
  
}))%>%distinct(Gene.Set)
cl_file<-"hub_mres_HMEC_CAGE_rich_GOBP_enrich_tbl.Rda"
cl_set_combo_tbl<-get(load(paste0(tbl_folder,cl_file))) %>% 
                    arrange(desc(OR))%>%filter(FDR<=0.01)%>%slice_head(n=15)

hm_gene_set<-as_tibble(read.table("~/Documents/multires_bhicect/data/epi_data/Gene_annotation/c5.all.v7.3.entrez.gmt",header = F,sep = "\t",fill=T))

cl_set_combo_tbl<-cl_set_combo_tbl %>%   left_join(.,hm_gene_set%>%dplyr::select(V1,V2),by=c("Gene.Set" = "V1"))



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
                                method="Rel")

library(seriation)
d<-as.dist(1/(simMatrix+1))
order <- seriate(d,method = "OLO")
image(simMatrix[get_order(order),get_order(order)])
new_order<-1:nrow(simMatrix)
names(new_order)<-cl_set_combo_tbl$GO.ID[get_order(order)]

id_gene_set_map<-cl_set_combo_tbl$Gene.Set
names(id_gene_set_map)<-cl_set_combo_tbl$GO.ID

clean_name<-unlist(lapply(lapply(strsplit(cl_set_combo_tbl$Gene.Set,split="_"),'[',-1),function(x)paste(x,collapse = " ")))
id_gene_set_map2<-clean_name
names(id_gene_set_map2)<-cl_set_combo_tbl$GO.ID

cl_set_FDR_tbl<-do.call(bind_rows,lapply (foi, function(cl_file){
  #cl_set<-unlist(lapply(strsplit(cl_file,split="_"),'[',1))
  tmp_name<-strsplit(cl_file,split="_")[[1]]
  cl_set<-paste(tmp_name[1:(length(tmp_name)-5)],collapse = "_")
  path_tbl<-get(load(paste0(tbl_folder,cl_file)))
  return(path_tbl %>% filter(Gene.Set %in% cl_set_combo_tbl$Gene.Set)%>%arrange(FDR) %>% mutate(set=cl_set))
  
}))

gg_heat<-cl_set_FDR_tbl %>% 
  mutate(Gene.Set2=unlist(lapply(lapply(strsplit(.$Gene.Set,split="_"),'[',-1),function(x)paste(x,collapse = " "))))%>%
  mutate(Gene.Set=fct_relevel(Gene.Set,id_gene_set_map[names(new_order)]),Gene.Set2=fct_relevel(Gene.Set2,id_gene_set_map2[names(new_order)]))%>%
  #mutate(set=fct_relevel(set,c("hub_50kb","tad_GOBP","inter_GOBP","hub_ex","tad_ex")))%>%
  mutate(OR=ifelse(FDR>0.01,NA,OR))%>%
  ggplot(.,aes(set,Gene.Set2,fill=OR))+geom_tile()+scale_fill_viridis_c()+ ylab("Gene Set")+
  theme(axis.text.y = element_text(size=10))+theme_minimal()
gg_heat
ggsave("~/Documents/multires_bhicect/weeklies/IFI_meeting/img/HMEC_GOBP_heat.svg",gg_heat,width = 30,height=20, units = "cm")
