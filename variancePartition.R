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
## VariancePartitioning

expr <- as.matrix(assays(sce)$logcounts)
meta_sub <- as.data.frame(colData(sce)[, c(celltype, batch)])
form <- as.formula(paste0("~ (1|", celltype, ") + (1|", batch, ")"))
varPart <- fitExtractVarPartModel(expr, form, meta_sub)
#
# #Sort variables (i.e. columns) by median fraction# of variance explained
vp <- varPart
#
vp_names <- rownames(vp)
vp <-vp  %>% dplyr::mutate(gene= vp_names) %>% dplyr::arrange(-!! rlang::parse_expr(batch))
vp_sub <- vp[1:3] %>% set_rownames(vp$gene)
## Summarize vp
# How many genes have a variance component affected by batch with > 1% ( resp. 10%)
n_batch_gene <- as_tibble(vp) %>%
  dplyr::filter(!! rlang::parse_expr(batch) > 0.01) %>% nrow()/n_genes

n_batch_gene10 <- as_tibble(vp) %>%
  dplyr::filter(!! rlang::parse_expr(batch) > 0.1) %>% nrow()/n_genes

# How many genes have a variance component affected by celltype with > 1%
n_celltype_gene <- as_tibble(vp) %>%
  dplyr::filter(!! rlang::parse_expr(celltype)> 0.01) %>% nrow()/n_genes

# How many genes have a variance component >1% explained by the batch effect relative to the celltype
n_rel <- n_batch_gene/n_celltype_gene


# Mean variance that is explained by the batch effect/celltype
m_batch <- mean(asin(vp[, batch]))
m_celltype <- mean(asin(vp[, celltype]))
m_rel <- m_batch/m_celltype


#Median variance explained by the batch effect/celltype

me_batch <- median(asin(vp[, batch]))
me_celltype <- median(asin(vp[, celltype]))
me_rel <- me_batch/me_celltype


params_summary_vp<-list(vp_sub = vp_sub, n_rel=n_rel,n_batch_gene=n_batch_gene,
                        n_batch_gene10=n_batch_gene10,n_celltype_gene=n_celltype_gene,
                        m_batch=m_batch,m_celltype=m_celltype,
                        me_batch=me_batch,me_celltype=me_celltype)


saveRDS(params_summary_vp , file=args[3])

