#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages({
  library(CellBench)
  library(ggplot2)
  library(scater)
  library(jcolors)
  library(CellMixS)
  library(gridExtra)
  library(purrr)
  library(jcolors)
  library(here)
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(variancePartition)
  library(diffcyt)
  library(ComplexHeatmap)
})
sce <- readRDS(file=args[1])
params<-readRDS(file=args[2])
attach(params)
## Celltype abundance - relative difference in abundances

meta_tib <- as_tibble(colData(sce)) %>% group_by_at(c(batch, celltype)) %>% 
  summarize(n = n()) %>% spread_(batch,"n")
meta_df <- as.data.frame(eval(meta_tib))[,-1]
meta_comb<-combn(meta_df,2,simplify=FALSE)
res<-lapply(meta_comb,function(x){
  cond1<-names(x)[1]
  cond2<-names(x)[2]
  rel_abund_diff<-mapply(function(cond1,cond2) abs(cond1-cond2)/(cond1+cond2), x[,cond1], x[,cond2])
  rel_abund_diff
})

mean_rel_abund_diff<-mean(unlist(res))
min_rel_abund_diff<-min(unlist(res))
max_rel_abund_diff<-max(unlist(res))


rel_diff_abund<-list(mean_rel_abund_diff=mean_rel_abund_diff,min_rel_abund_diff=min_rel_abund_diff,
                     max_rel_abund_diff=max_rel_abund_diff)
saveRDS(rel_diff_abund,file=args[3])

