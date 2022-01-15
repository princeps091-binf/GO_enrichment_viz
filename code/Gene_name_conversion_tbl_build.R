renv::install("tidyverse")
renv::install("bioc::GenomicRanges")
renv::install("furrr")
renv::install("bioc::biomaRt")

library(biomaRt)
library(tidyverse)
library(GenomicRanges)
library(furrr)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
gene_GRange_file<-"./data/CAGE_gene_GRange.Rda"
gene_set1_file<-"./data/GOBP_gene_set_l.Rda"
gene_set2_file<-"./data/Hallmark_gene_set_l.Rda"
gene_set3_file<-"~/Documents/multires_bhicect/data/epi_data/hg19_ENST_to_ENSG.txt"

out_file<-"./data/gene_name_conv_tbl.Rda"

hg19_ENST_to_ENSG <- read_delim(gene_set3_file, 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
hg19_gene<-unique(hg19_ENST_to_ENSG$name2)

gene_GRange<-get(load(gene_GRange_file))
tmp_obj<-names(mget(load(gene_GRange_file)))
rm(list=tmp_obj)
rm(tmp_obj)
gene_set<-unique(unlist(mcols(gene_GRange)$ENSG))

gene_set1_l<-get(load(gene_set1_file))
tmp_obj<-names(mget(load(gene_GRange_file)))
rm(list=tmp_obj)
rm(tmp_obj)
gene_set1<-unique(unlist(gene_set1_l))

gene_set2_l<-get(load(gene_set2_file))
tmp_obj<-names(mget(load(gene_GRange_file)))
rm(list=tmp_obj)
rm(tmp_obj)
gene_set2<-unique(unlist(gene_set2_l))

ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl")

genes_hmec <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=gene_set,
  mart=ensembl)

genes_hg19 <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=hg19_gene,
  mart=ensembl)

genes_hm1 <- getBM(
  filters="entrezgene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=gene_set1,
  mart=ensembl)


genes_hm2 <- getBM(
  filters="entrezgene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=gene_set2,
  mart=ensembl)


tot_gene_conv<-tibble(genes_hmec%>%full_join(.,genes_hm1)%>%full_join(.,genes_hg19) %>% full_join(.,genes_hm2)) 

save(tot_gene_conv,file=out_file)
