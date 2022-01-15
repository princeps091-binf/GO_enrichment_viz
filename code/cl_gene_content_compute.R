renv::install("tidyverse")
renv::install("bioc::GenomicRanges")
renv::install("furrr")

library(tidyverse)
library(GenomicRanges)
library(furrr)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
cl_GRange_file<-"./data/TAD_H1_CAGE_rich_GRange.Rda"
gene_GRange_file<-"./data/CAGE_H1_gene_GRange.Rda"
out_file<-"./data/TAD_H1_CAGE_rich_ENSG_tbl.Rda"
cl_tbl<-get(load(cl_GRange_file))
tmp_obj<-names(mget(load(cl_GRange_file)))
rm(list=tmp_obj)
rm(tmp_obj)

gene_GRange<-get(load(gene_GRange_file))
tmp_obj<-names(mget(load(gene_GRange_file)))
rm(list=tmp_obj)
rm(tmp_obj)
plan(multisession, workers = 3)


cl_tbl<-cl_tbl %>% mutate(ENSG=future_map(Grange,function(x){
  return(unique(unlist(gene_GRange@elementMetadata$ENSG[unique(subjectHits(findOverlaps(x,gene_GRange)))])))
  
}))

save(cl_tbl,file=out_file)
