library(tidyverse)
library(readr)
library(vroom)
library(parallel)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
# Get current version of gene sets fro MSigDB websites:
# http://www.gsea-msigdb.org/gsea/downloads.jsp

gene_set_file<-"./data/MSigDB_tbl/c8.all.v7.3.entrez.gmt"
out_file<-"./data/Cell_type_gene_set_l.Rda"

# Input and format gene sets of interest
hm_gene_set <- readLines(gene_set_file, warn = FALSE)

#hm_gene_set<-read_delim(gene_set_file,header = F,delim="\t")

cl<-makeCluster(5)
clusterEvalQ(cl, {
  library(dplyr)
})
clusterExport(cl,c('hm_gene_set'))
hallmark_set<-parLapply(cl,1:length(hm_gene_set),function(x){
  tmp<-hm_gene_set[x]
  return(unlist(strsplit(tmp,split = "\\t"))[-c(1,2)])
})
stopCluster(cl)
rm(cl)
names(hallmark_set)<-unlist(lapply(1:length(hm_gene_set),function(x){
  tmp<-hm_gene_set[x]
  return(unlist(strsplit(tmp,split = "\\t"))[1])
  
}))
# Subset the Biological Processes gene sets
save(hallmark_set,file=out_file)
