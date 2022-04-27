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
background_gene_file<-"./data/CAGE_GM12878_entrez_gene_GRange.Rda"
foreground_gene_file<-"./data/GM12878_trans_res_hub_entrez_tbl.Rda"

top_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/GM12878_union_top_trans_res_dagger_tbl.Rda"

out_file<-"./data/trans_res_hub_GS_tbl/H1_trans_res_hub_entrez_GOBP_enrich_tbl.Rda"

gene_set_file<-"./data/Hallmark_gene_set_l.Rda"


background_GRange<-get_obj_in_fn(background_gene_file)

foreground_gene_tbl<-get_obj_in_fn(foreground_gene_file) 

top_hub_tbl<-get_obj_in_fn(top_hub_file)

foreground_gene_tbl<-foreground_gene_tbl %>% 
  inner_join(.,top_hub_tbl)

res_foreground_gene_tbl<-foreground_gene_tbl#  %>% filter(res=="500kb")

Gene_set_l<-get_obj_in_fn(gene_set_file)

foreground_gene_vec<-unique(unlist(res_foreground_gene_tbl$entrez.content))

background_gene_vec<-unique(unlist(mcols(background_GRange)$entrez))


path_tbl<-GO_set_enrich_fn(foreground_gene_vec,background_gene_vec,Gene_set_l)
print(path_tbl %>% 
        filter(FDR<=0.01) %>% 
        arrange(desc(OR))
      #        arrange(FDR)
      ,n=100)
