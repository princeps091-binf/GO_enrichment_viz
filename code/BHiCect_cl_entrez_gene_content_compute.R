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
cl_tbl_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/GM12878_union_top_trans_res_dagger_tbl.Rda"

cl_spec_res_folder<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"

gene_GRange_file<-"./data/CAGE_GM12878_entrez_gene_GRange.Rda"

out_file<-paste0("./data/GM12878_trans_res_CAGE_top_hub_entrez_tbl.Rda")

cl_tbl<-get_obj_in_fn(cl_tbl_file) %>% 
  mutate(res=str_split_fixed(node,"_",2)[,1])

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
  mutate(entrez.content=future_map(GRange,function(x){
    return(unique(unlist(mcols(gene_GRange)$entrez[unique(subjectHits(findOverlaps(x,gene_GRange)))])))
  }))
plan(sequential)

cl_tbl
save(cl_tbl,file=out_file)
