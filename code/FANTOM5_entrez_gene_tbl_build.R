library(GenomicRanges)
library(vroom)
library(tidyverse)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
out_file<-"./data/FANTOM5_entrez_gene_tbl.Rda"

FANTOM5_tpm_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/hg19.cage_peak_phase1and2combined_counts_ann.osc.txt",delim="\t",comment = "#",col_select = 1:7)

FANTOM5_entrez_tbl<-FANTOM5_tpm_tbl %>% 
  dplyr::select(`00Annotation`,entrezgene_id) %>% 
  dplyr::rename(peak=`00Annotation`) %>% 
  mutate(entrez.l=str_split(entrezgene_id,",")) %>% 
  mutate(entrez.id= map(entrez.l,function(x){
    str_extract(x,"[0-9]+")
  })) %>% 
  dplyr::select(peak,entrez.id) %>% 
  mutate(chr=str_split_fixed(peak,":|\\.\\.|,",4)[,1],
         start=as.numeric(str_split_fixed(peak,":|\\.\\.|,",4)[,2]),
         end=as.numeric(str_split_fixed(peak,":|\\.\\.|,",4)[,3])) %>% 
  distinct() %>% 
  filter(!(is.na(start)))
#-------------------------------------------------------------------
cage_Grange_fn<-function(cage_hmec_a){
  cage_a<-cage_hmec_a%>%filter(!(is.na(start)))
  full_cage_Grange<-   GRanges(seqnames=cage_a$chr,
                               ranges = IRanges(start=cage_a$start,
                                                end=cage_a$end,
                                                names=paste(cage_a$chr,1:nrow(cage_a),sep='_')
                               ))
  mcols(full_cage_Grange)<-tibble(entrez=cage_a$entrez.id)
  return(full_cage_Grange)
}
#-------------------------------------------------------------------

cage_coord_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/CAGE_tss_coord_GM12878_tbl.Rda"
out_file<-"./data/CAGE_GM12878_entrez_gene_GRange.Rda"

cage_tbl<-get(load(cage_coord_file))
tmp_obj<-names(mget(load(cage_coord_file)))
rm(list=tmp_obj)
rm(tmp_obj)

cage_tbl<-cage_tbl %>% 
  left_join(.,FANTOM5_entrez_tbl,by=c("Id"="peak","chr"="chr","start"="start","end"="end")) %>% 
  dplyr::select(Id,chr,start,end,entrez.id)
full_cage_Grange<-cage_Grange_fn(cage_tbl)

save(full_cage_Grange,file=out_file)
