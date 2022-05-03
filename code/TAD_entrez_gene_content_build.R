library(tidyverse)
library(GenomicRanges)
library(furrr)
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
cl_tbl_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_GM12878_TAD_pval_tbl.Rda"
gene_GRange_file<-"./data/CAGE_GM12878_entrez_gene_GRange.Rda"
out_file<-paste0("./data/GM12878_TAD_CAGE_rich_entrez_tbl.Rda")

cl_tbl<-get_obj_in_fn(cl_tbl_file)

gene_GRange<-get_obj_in_fn(gene_GRange_file)


cl_coord_tbl<-cl_tbl %>% 
  mutate(start=as.numeric(str_split_fixed(ID,"_",3)[,2]),
         end=as.numeric(str_split_fixed(ID,"_",3)[,3])) %>%
  dplyr::select(ID,chr,start,end)
  
TAD_GRange<-GRanges(seqnames=cl_coord_tbl$chr,
                    ranges = IRanges(start=cl_coord_tbl$start,
                                     end=cl_coord_tbl$end
                    ))

tmp_tbl<-findOverlaps(TAD_GRange,gene_GRange) %>% as_tibble %>% 
  mutate(TAD.ID=cl_coord_tbl$ID[queryHits]) %>% 
  group_by(TAD.ID) %>% 
  summarise(entrez.content=list(unique(unlist(mcols(gene_GRange)$entrez[subjectHits]))))

cl_tbl<-cl_tbl %>% 
  left_join(.,tmp_tbl,by=c("ID"="TAD.ID"))

cage_rich_TAD_tbl<-cl_tbl %>% 
  mutate(FDR=p.adjust(emp.pval,method='fdr')) %>% 
  filter(FDR<=0.01)
save(cage_rich_TAD_tbl,file=out_file)
  
