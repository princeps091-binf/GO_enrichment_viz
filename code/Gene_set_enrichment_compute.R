library(tidyverse)
library(GenomicRanges)
#------------------------------
get_obj_in_fn<-function(file){
  out_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

GO_set_enrich_fn<-function(cl_set_gene,cage_active_genes_vec,GOBP_set){
  fn_env<-environment()
  
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(dplyr)
    print('node ready')
  })
  clusterExport(cl,c('cl_set_gene','cage_active_genes_vec'),envir = fn_env)
  go_pval<-parLapply(cl,GOBP_set,function(tmp_set){
    hitInSample<-sum(cl_set_gene %in% tmp_set)
    sampleSize<-length(cl_set_gene)
    hitInPop<-sum(cage_active_genes_vec %in% tmp_set)
    failInPop<-length(cage_active_genes_vec) - hitInPop
    p_val<-phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
    OR_GO<-(hitInSample/sampleSize)/(hitInPop/length(cage_active_genes_vec))
    return(tibble(p.val=p_val,OR=OR_GO,in.gene=hitInSample))
  })
  stopCluster(cl)
  rm(cl)
  path_tbl<-do.call(bind_rows,go_pval)%>%mutate(Gene.Set=names(go_pval),FDR=p.adjust(p.val,method='fdr'))%>%dplyr::select(Gene.Set,FDR,OR,in.gene)
  return(path_tbl)
}

#------------------------------
background_gene_file<-"./data/CAGE_HMEC_gene_GRange.Rda"
foreground_gene_file<-"./data/HMEC_poisson_05_hub_ENSG_tbl.Rda"

out_file<-"./data/HMEC_TAD_GOBP_enrich_tbl.Rda"

gene_conv_tbl_file<-"./data/gene_name_conv_tbl.Rda"

gene_set_file<-"./data/GOBP_gene_set_l.Rda"


background_GRange<-get_obj_in_fn(background_gene_file)

foreground_gene_tbl<-get_obj_in_fn(foreground_gene_file)

gene_conv_tbl<-get_obj_in_fn(gene_conv_tbl_file)

Gene_set_l<-get_obj_in_fn(gene_set_file)


cage_active_gene<-unique(unlist(mcols(background_GRange)$ENSG))

foreground_gene_vec<-foreground_gene_tbl %>% 
  dplyr::select(ENSG.content) %>% 
  unnest(cols=c(ENSG.content)) %>% distinct() %>% 
  inner_join(.,gene_conv_tbl,by=c("ENSG.content"="ensembl_gene_id")) %>% 
  distinct(entrezgene_id) %>% unlist %>% as.character

background_gene_vec<-gene_conv_tbl%>%
  filter(ensembl_gene_id %in% cage_active_gene)%>%
  distinct(entrezgene_id) %>% unlist %>% as.character

#hg19_vec<-gene_conv_tbl %>% distinct(entrezgene_id) %>% unlist


path_tbl<-GO_set_enrich_fn(foreground_gene_vec,background_gene_vec,Gene_set_l)
print(path_tbl %>% 
  filter(FDR<=0.01) %>% 
  arrange(FDR),n=154)
save(path_tbl,file=out_file)

