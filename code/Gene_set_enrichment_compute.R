renv::install("tidyverse")
library(tidyverse)

background_gene_file<-"./data/CAGE_HMEC_gene_GRange.Rda"
foreground_gene_file<-"./data/mres_HMEC_hub_ENSG_tbl.Rda"
out_file<-"./data/hub_mres_HMEC_CAGE_rich_GOBP_enrich_tbl.Rda"

gene_conv_tbl_file<-"./data/gene_name_conv_tbl.Rda"

gene_set_file<-"./data/GOBP_gene_set_l.Rda"


background_GRange<-get(load(background_gene_file))
tmp_obj<-names(mget(load(background_gene_file)))
rm(list=tmp_obj)
rm(tmp_obj)

foreground_gene_tbl<-get(load(foreground_gene_file))
tmp_obj<-names(mget(load(foreground_gene_file)))
rm(list=tmp_obj)
rm(tmp_obj)

gene_conv_tbl<-get(load(gene_conv_tbl_file))
tmp_obj<-names(mget(load(gene_conv_tbl_file)))
rm(list=tmp_obj)
rm(tmp_obj)

Gene_set_l<-get(load(gene_set_file))
tmp_obj<-names(mget(load(gene_set_file)))
rm(list=tmp_obj)
rm(tmp_obj)


cage_active_gene<-unique(unlist(mcols(background_GRange)$ENSG))

foreground_gene_vec<-foreground_gene_tbl %>% dplyr::select(ENSG) %>% 
  unnest(cols=c(ENSG)) %>% distinct() %>% 
  inner_join(.,gene_conv_tbl,by=c("ENSG"="ensembl_gene_id")) %>% 
  distinct(entrezgene_id) %>% unlist %>% as.character

background_gene_vec<-gene_conv_tbl%>%
  filter(ensembl_gene_id %in% cage_active_gene)%>%
  distinct(entrezgene_id) %>% unlist %>% as.character


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


path_tbl<-GO_set_enrich_fn(foreground_gene_vec,background_gene_vec,Gene_set_l)
path_tbl %>% arrange(FDR)
save(path_tbl,file=out_file)

