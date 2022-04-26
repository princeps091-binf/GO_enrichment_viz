library(tidyverse)
library(GenomicRanges)
library(rjson)
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
background_gene_file<-"./data/CAGE_HMEC_entrez_gene_GRange.Rda"
foreground_gene_file<-"./data/HMEC_trans_res_hub_entrez_tbl.Rda"

out_file<-"./data/trans_res_hub_GS_tbl/HMEC_trans_res_hub_entrez_GOBP_enrich_tbl.Rda"

gene_set_file<-"./data/GOBP_gene_set_l.Rda"


result <- fromJSON(file = "./data/MSigDB_tbl/c5.go.bp.v7.5.1.json")
GO_ID_map_vec<-map_chr(result,function(x){
  x$exactSource
})

background_GRange<-get_obj_in_fn(background_gene_file)

foreground_gene_tbl<-get_obj_in_fn(foreground_gene_file) %>% 
#  filter(!(res%in% c("1Mb","500kb","100kb")))
#  filter(!(res %in% c("5kb","10kb")))
  filter(res == "1Mb")
Gene_set_l<-get_obj_in_fn(gene_set_file)

foreground_gene_vec<-unique(unlist(foreground_gene_tbl$entrez.content))

background_gene_vec<-unique(unlist(mcols(background_GRange)$entrez))

#hg19_vec<-gene_conv_tbl %>% distinct(entrezgene_id) %>% unlist


path_tbl<-GO_set_enrich_fn(foreground_gene_vec,background_gene_vec,Gene_set_l)
print(path_tbl %>% 
        filter(FDR<=0.01) %>% 
#        arrange(desc(OR))
              arrange(FDR)
      ,n=100)
#------------------------------------------------
cand_GO<-path_tbl %>% 
  filter(FDR<=0.01) %>%
  dplyr::select(Gene.Set) %>% 
  unlist
cand_GO_ID<-GO_ID_map_vec[cand_GO]
library(simplifyEnrichment)

mat <- GO_similarity(cand_GO_ID,measure="Resnik",ont = "BP")
df <- simplifyGO(mat)
