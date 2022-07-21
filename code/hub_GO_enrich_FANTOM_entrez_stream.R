library(tidyverse)
library(GenomicRanges)
library(furrr)
library(parallel)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

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
active_gene_file<-"./data/CAGE_GM12878_entrez_gene_GRange.Rda"
hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/GM12878_union_trans_res_dagger_tbl.Rda"
spec_res_folder<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"

gene_set_file<-"./data/GOBP_gene_set_l.Rda"

Gene_set_l<-get_obj_in_fn(gene_set_file)

active_peak_GRange<-get_obj_in_fn(active_gene_file)

hub_tbl<-get_obj_in_fn(hub_file)

hub_tbl<-do.call(bind_rows,map(unique(hub_tbl$chr),function(chromo){
  message(chromo)
  base::load(paste0(spec_res_folder,chromo,"_spec_res.Rda"))
  tmp_tbl<-hub_tbl %>% 
    filter(chr==chromo) %>% 
    mutate(bins=chr_spec_res$cl_member[node]) %>% 
    mutate(bins=map(bins,as.numeric)) 
  
})) %>% 
  dplyr::select(chr,node,res,bins)

plan(multisession,workers=4)
hub_tbl<-hub_tbl %>% 
  mutate(GRange=future_pmap(list(chr,bins,res),function(chr,bins,res){
    inter_cl_Grange<-   GRanges(seqnames=chr,
                                ranges = IRanges(start=as.numeric(bins),
                                                 end=as.numeric(bins) + res_num[res]-1
                                ))
    inter_cl_Grange<-GenomicRanges::reduce(inter_cl_Grange)
    return(inter_cl_Grange)
    
  }))
plan(sequential)

tmp_l<-hub_tbl %>% 
  #  filter(res == "100kb") %>% 
  #  filter(res%in% c("1Mb","500kb","100kb")) %>% 
  #  filter(res%in% c("10kb","50kb","5kb")) %>% 
  dplyr::select(GRange) %>% as.list

cl_GRange<-IRanges::reduce(do.call("c",tmp_l$GRange))

peaks_in_hub_entrez_tbl<-findOverlaps(active_peak_GRange,cl_GRange) %>% 
  as_tibble %>% 
  distinct(queryHits) %>% 
  mutate(entrez.gene=mcols(active_peak_GRange)$entrez[queryHits])
in_hub_entrez_vec<-unique(unlist(peaks_in_hub_entrez_tbl$entrez.gene))
out_hub_entrez_vec<-unique(unlist(mcols(active_peak_GRange)$entrez[-peaks_in_hub_entrez_tbl$queryHits]))

full_active_bg<-unique(unlist(mcols(active_peak_GRange)$entrez))

path_tbl<-GO_set_enrich_fn(in_hub_entrez_vec,full_active_bg,Gene_set_l)
print(path_tbl %>% 
        filter(FDR<=0.01) %>% 
        arrange(desc(OR))
      #        arrange(FDR)
      ,n=100)
