renv::install("tidyverse")
renv::install("bioc::GenomicRanges")
renv::install("bioc::biomaRt")


library(biomaRt)
library(tidyverse)
library(GenomicRanges)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
# Get current version of gene sets fro MSigDB websites:
# http://www.gsea-msigdb.org/gsea/downloads.jsp

gene_set_file<-"./data/MSigDB_tbl/h.all.v7.5.1.entrez.gmt"
out_file<-"./data/Hallmark_gene_set_l.Rda"

# Input and format gene sets of interest
hm_gene_set<-as_tibble(read.table(gene_set_file,header = F,sep = "\t",fill=T))

cl<-makeCluster(5)
clusterEvalQ(cl, {
  library(dplyr)
})
clusterExport(cl,c('hm_gene_set'))
hallmark_set<-parLapply(cl,1:nrow(hm_gene_set),function(x){
  tmp<-hm_gene_set[x,-c(1,2)]
  return(tmp[!(is.na(tmp))])
})
stopCluster(cl)
rm(cl)
names(hallmark_set)<-hm_gene_set$V1
# Subset the Biological Processes gene sets
GOBP_set<-hallmark_set[grep("GOBP",names(hallmark_set),value=T)]
save(hallmark_set,file=out_file)
