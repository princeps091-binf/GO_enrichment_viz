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
cl_tbl_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/candidate_compound_hub/HMEC_5kb_tss_compound_hub.Rda"

cl_spec_res_folder<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"

gene_GRange_file<-"./data/CAGE_HMEC_gene_GRange.Rda"

out_file<-paste0("./data/HMEC_5kb_top_compound_hub_ENSG_tbl.Rda")

cl_tbl<-get_obj_in_fn(cl_tbl_file)

top_compound_hub_5kb_tbl<-do.call(bind_rows,map(unique(cl_tbl$chr),function(chromo){
  tmp_tbl<-cl_tbl %>% 
    filter(chr==chromo)
  return(tmp_tbl %>% 
           filter(parent.hub %in% tmp_tbl$parent.hub[which(!(tmp_tbl$parent.hub %in% unique(unlist(tmp_tbl$ch.hub))))])
  )
  
}))

gene_GRange<-get_obj_in_fn(gene_GRange_file)

#Build BHiCect GRange
chr_set<-unique(top_compound_hub_5kb_tbl$chr)
chr_res_l<-vector('list',length(chr_set))
names(chr_res_l)<-chr_set
for (chromo in chr_set){
  message(chromo)
  load(paste0(cl_spec_res_folder,chromo,"_spec_res.Rda"))
  chr_cl_tbl<-top_compound_hub_5kb_tbl %>% filter(chr==chromo) %>% mutate(bins=chr_spec_res$cl_member[parent.hub])
  plan(multisession, workers = 3)
  
  chr_res_l[[chromo]]<-chr_cl_tbl %>% mutate(GRange=future_pmap(list(chr,bins,res),function(chr,bins,res){
    return(GRanges(seqnames=chr,
                   ranges = IRanges(start=as.numeric(bins),
                                    end=as.numeric(bins)+res_num[res]-1
                   )))
    
    
  }))
  
}

top_compound_hub_5kb_tbl<-do.call(bind_rows,chr_res_l)

top_compound_hub_5kb_tbl<-top_compound_hub_5kb_tbl %>% 
  mutate(ENSG.content=future_map(GRange,function(x){
    return(unique(unlist(gene_GRange@elementMetadata$ENSG[unique(subjectHits(findOverlaps(x,gene_GRange)))])))
  }))
top_compound_hub_5kb_tbl
save(top_compound_hub_5kb_tbl,file=out_file)
