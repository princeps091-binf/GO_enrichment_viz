library(tidyverse)
library(GenomicRanges)
library(furrr)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------
get_obj_in_fn<-function(file){
  tmp_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(tmp_tbl)
}
Build_GRange_fn<-function(bin_set,res,chr,res_num){
  
  return(GRanges(seqnames=chr,
                 ranges = IRanges(start=as.numeric(bin_set),
                                  end=as.numeric(bin_set)+res_num[res]-1
                 )))
  
}
#-------------------------------
compound_mres_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/HMEC_mres_coumpund_hub.Rda"
pval_tbl_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_HMEC_pval_tbl.Rda"
gene_GRange_file<-"./data/CAGE_HMEC_gene_GRange.Rda"

out_file<-paste0("./data/HMEC_compound_hub_ENSG_tbl.Rda")


gene_GRange<-get_obj_in_fn(gene_GRange_file)
compound_mres_hub_tbl<-get_obj_in_fn(compound_mres_hub_file)
pval_tbl<-get_obj_in_fn(pval_tbl_file)


top_hub_tbl<-compound_mres_hub_tbl %>% 
  distinct(chr,top_hub,top_res) %>% 
  filter(!(is.na(top_hub)))

top_hub_tbl<-top_hub_tbl %>% 
  left_join(.,pval_tbl %>% 
              dplyr::select(chr,cl,GRange,emp.pval),by=c("chr"='chr','top_hub'="cl"))

top_hub_GRange<-GenomicRanges::reduce(top_hub_tbl %>% 
                                        filter(top_res %in% c("1Mb","500kb")) %>% dplyr::select(GRange) %>% unlist %>% GRangesList %>% unlist)


ENSG_vec<-unique(unlist(gene_GRange@elementMetadata$ENSG[unique(subjectHits(findOverlaps(top_hub_GRange,gene_GRange)))]))
hub_gene_tbl<-tibble(cl=paste0("compound_hub"),ENSG=list(ENSG_vec))

save(hub_gene_tbl,file=out_file)
