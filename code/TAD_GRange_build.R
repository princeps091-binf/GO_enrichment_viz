renv::install("tidyverse")
renv::install("bioc::GenomicRanges")

library(tidyverse)
library(GenomicRanges)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
cl_tbl_file<-"~/Documents/multires_bhicect/data/epi_data/H1/CAGE/TAD_CAGE_pval_peak_shuffle_tbl.Rda"
out_file<-"./data/TAD_H1_CAGE_rich_GRange.Rda"
tad_tbl<-get(load(cl_tbl_file))
tmp_obj<-names(mget(load(cl_tbl_file)))
rm(list=tmp_obj)
rm(tmp_obj)

TAD_GRange_build_fn<-function(chr,start,end){
  nuc_Grange <- GRanges(seqnames=chr,
                        ranges = IRanges(start=start,
                                         end=end,
                                         names=paste(chr,start,sep='_')
                        ))
  return(nuc_Grange) 
  
}

# tad_tbl<-tad_tbl %>% mutate(tad2=tad) %>% 
#   tidyr::separate(tad2,into = c("chr",'start','end'),sep = '_') %>% 
#   mutate(start=as.numeric(start),end=as.numeric(end)) %>% 
#   mutate(Grange=pmap(list(chr,start,end),function(chr,start,end)TAD_GRange_build_fn(chr,start,end))) %>% 
#   dplyr::rename(ID=tad)


tad_tbl<-tad_tbl %>%
  dplyr::rename(chr=chr1,start=x1,end=x2) %>%
  filter(FDR<=0.01) %>% 
  mutate(Grange=pmap(list(chr,start,end),function(chr,start,end)TAD_GRange_build_fn(chr,start,end)))

save(tad_tbl,file=out_file)
