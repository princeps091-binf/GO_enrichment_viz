renv::install("tidyverse")
renv::install("bioc::GenomicRanges")

library(tidyverse)
library(GenomicRanges)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
cage_coord_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/CAGE_tss_coord_HMEC_tbl.Rda"
out_file<-"./data/CAGE_HMEC_gene_GRange.Rda"

cage_tbl<-get(load(cage_coord_file))
tmp_obj<-names(mget(load(cage_coord_file)))
rm(list=tmp_obj)
rm(tmp_obj)

cage_Grange_fn<-function(cage_hmec_a){
  cage_a<-cage_hmec_a%>%filter(!(is.na(start)))
  full_cage_Grange<-   GRanges(seqnames=cage_a$chr,
                               ranges = IRanges(start=cage_a$start,
                                                end=cage_a$end,
                                                names=paste(cage_a$chr,1:nrow(cage_a),sep='_')
                               ))
  mcols(full_cage_Grange)<-tibble(ENSG=cage_a$ENSG)
  return(full_cage_Grange)
}
full_cage_Grange<-cage_Grange_fn(cage_tbl)
save(full_cage_Grange,file=out_file)
