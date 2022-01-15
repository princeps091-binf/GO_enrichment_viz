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
cl_tbl_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/CAGE/dagger_mres_fdr_01_multi_cagebin_tbl.Rda"
cl_spec_res_folder<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"

gene_GRange_file<-"./data/CAGE_GM12878_gene_GRange.Rda"
out_file<-"./data/mres_GM12878_hub_ENSG_tbl.Rda"

cl_tbl<-get(load(cl_tbl_file))
tmp_obj<-names(mget(load(cl_tbl_file)))
rm(list=tmp_obj)
rm(tmp_obj)

gene_GRange<-get(load(gene_GRange_file))
tmp_obj<-names(mget(load(gene_GRange_file)))
rm(list=tmp_obj)
rm(tmp_obj)

#Build BHiCect GRange
chr_set<-unique(cl_tbl$chr)
chr_res_l<-vector('list',length(chr_set))
names(chr_res_l)<-chr_set
for (chromo in chr_set){
  print(chromo)
  load(paste0(cl_spec_res_folder,chromo,"_spec_res.Rda"))
  chr_cl_tbl<-cl_tbl %>% filter(chr==chromo) %>% mutate(bins=chr_spec_res$cl_member[node])
  plan(multisession, workers = 3)
  
  chr_cl_tbl<-chr_cl_tbl %>% mutate(GRange=future_pmap(list(chr,bins,res),function(chr,bins,res){
    return(GRanges(seqnames=chr,
                   ranges = IRanges(start=as.numeric(bins),
                                    end=as.numeric(bins)+res_num[res]-1
                   )))
    
    
  }))
  chr_res_l[[chromo]]<-chr_cl_tbl %>% distinct(chr) %>% mutate(Grange=list(reduce(unlist(GRangesList(chr_cl_tbl$GRange)))))
}

grange_tbl<-do.call(bind_rows,chr_res_l)
cl_GRange<-reduce(unlist(GRangesList(grange_tbl$Grange)))
ENSG_vec<-unique(unlist(gene_GRange@elementMetadata$ENSG[unique(subjectHits(findOverlaps(cl_GRange,gene_GRange)))]))
hub_gene_tbl<-tibble(cl="GM12878_hub_mres",ENSG=list(ENSG_vec))
save(hub_gene_tbl,file=out_file)
