library(tidyverse)
library(GenomicRanges)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
get_obj_in_fn<-function(file){
  out_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

#-------------------------------------------------------------------
cl_tbl_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_HMEC_TAD_pval_tbl.Rda"
gene_GRange_file<-"./data/CAGE_HMEC_gene_GRange.Rda"
out_file<-"./data/TAD_HMEC_CAGE_rich_GRange.Rda"

tad_tbl<-get_obj_in_fn(cl_tbl_file)%>% 
  mutate(FDR=p.adjust(emp.pval,method='fdr')) %>% 
  filter(FDR<=0.01)

tad_GRange<-reduce(do.call("c",unlist(tad_tbl$GRange)))

gene_GRange<-get_obj_in_fn(gene_GRange_file)

TAD_ensg_content_tbl<-tibble(TAD='GM12878',ENSG=list(unique(unlist(mcols(gene_GRange)$ENSG[unique(subjectHits(findOverlaps(tad_GRange,gene_GRange)))]))))

save(TAD_ensg_content_tbl,file=out_file)
