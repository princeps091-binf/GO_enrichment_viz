library(tidyverse)
library(GenomicRanges)
library(furrr)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
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
cl_tbl_file<-"~/Documents/multires_bhicect/BHiCect_poisson_cluster_detect/data/pval_tbl/DAGGER/HMEC_poisson_DAGGER_01.Rda"

cl_spec_res_folder<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"

gene_GRange_file<-"./data/CAGE_HMEC_gene_GRange.Rda"

out_file<-paste0("./data/HMEC_poisson_hub_ENSG_tbl.Rda")

cl_tbl<-get_obj_in_fn(cl_tbl_file)

gene_GRange<-get_obj_in_fn(gene_GRange_file)

#Build BHiCect GRange
chr_set<-unique(cl_tbl$chr)
chr_res_l<-vector('list',length(chr_set))
names(chr_res_l)<-chr_set
for (chromo in chr_set){
  message(chromo)
  load(paste0(cl_spec_res_folder,chromo,"_spec_res.Rda"))
  chr_cl_tbl<-cl_tbl %>% filter(chr==chromo) %>% mutate(bins=chr_spec_res$cl_member[node])
  plan(multisession, workers = 3)
  
  chr_res_l[[chromo]]<-chr_cl_tbl %>% mutate(GRange=future_pmap(list(chr,bins,res),function(chr,bins,res){
    return(GRanges(seqnames=chr,
                   ranges = IRanges(start=as.numeric(bins),
                                    end=as.numeric(bins)+res_num[res]-1
                   )))
    
    
  }))
  plan(sequential)
}

cl_tbl<-do.call(bind_rows,chr_res_l)

plan(multisession, workers = 3)
cl_tbl<-cl_tbl %>% 
  mutate(ENSG.content=future_map(GRange,function(x){
    return(unique(unlist(mcols(gene_GRange)$ENSG[unique(subjectHits(findOverlaps(x,gene_GRange)))])))
  }))
plan(sequential)
cl_tbl
save(cl_tbl,file=out_file)
