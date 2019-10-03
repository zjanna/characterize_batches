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
# cms
sce <- cms(sce, group = batch, k = k, cell_min = cell_min, n_dim = n_dim)

sce$cms<-2
saveRDS(sce, file=args[3])
#summarize
mean_cms <- mean(sce$cms)
n_cms_0.01 <- length(which(sce$cms < 0.01))
cluster_mean_cms <- as_tibble(colData(sce)) %>% group_by_at(celltype) %>% summarize(cms_mean = mean(cms))
var_cms <- var(cluster_mean_cms$cms_mean)

cms_summary<-list(mean_cms=mean_cms,n_cms_0.01=n_cms_0.01,
                  cluster_mean_cms=cluster_mean_cms,var_cms=var_cms)
saveRDS(cms_summary,file=args[4])

