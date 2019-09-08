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
params_summary_vp<-readRDS(file=args[1])
attach(params_summary_vp)
rel_diff_abund<-readRDS(file=args[2])
attach(rel_diff_abund)
de_summary<-readRDS(file=args[3])
attach(de_summary)
cms_summary<-readRDS(file=args[4])
attach(cms_summary)
sim<-readRDS(file=args[5])
attach(sim)
sce<-readRDS(file=args[6])
params<-readRDS(file=args[7])
attach(params)

#Size? How much of the variance can be attributed to the batch effect?
size <- data.frame("batch_genes_1per" = n_batch_gene,
                   "batch_genes_10per" = n_batch_gene10,
                   "celltype_gene_1per" = n_celltype_gene,
                   "relative_batch_celltype" = n_rel,
                   "mean_var_batch" = m_batch,
                   "mean_var_celltype" = m_celltype,
                   "median_var_batch" = me_batch,
                   "median_var_celltype" = me_celltype,
                   "n_cells_total" = ncol(sce),
                   "n_genes_total" = nrow(sce))

#Celltype-specificity? How celltype/cluster specific are batch effects? Differences in sample variation between batches?
celltype_spec <- data.frame("mean_cms" = mean_cms,
                            'mean_rel_abund_diff' = mean_rel_abund_diff,
                            'min_rel_abund_diff' = min_rel_abund_diff,
                            'max_rel_abund_diff' = max_rel_abund_diff,
                            "celltype_var_cms" = var_cms,
                            "n_cells_cms_0.01" = n_cms_0.01)

#Gene-specificity? How do they effect genes? Single genes? All genes? Which genes?
gene <- data.frame("mean_mean_n_de_genes" = mean_mean_n_de,
                   "max_mean_n_de_genes" = max_mean_n_de,
                   "min_mean_n_de_genes" = min_mean_n_de,
                   "mean_n_genes_lfc1" = mean_n_genes_lfc1,
                   "min_n_genes_lfc1" = min_n_genes_lfc1,
                   "max_n_genes_lfc1" = max_n_genes_lfc1,
                   "mean_de_overlap" = mean_de_overlap,
                   "min_de_overlap" = min_de_overlap,
                   "max_de_overlap" = max_de_overlap,
                   "mean_rel_cluster_spec"= mean_rel_spec,
                   "min_rel_cluster_spec"= min_rel_spec,
                   "max_rel_cluster_spec"= max_rel_spec
)
# Cell-specificity? How cell-specific are batche effects? Are their differences in within celltype variation between batches?

# Cell-specificity? How cell-specific are batche effects? Are their differences in within celltype variation between batches?

sim <- data.frame("mean_p_be" = mean_p_be,
                  "max_p_be" = max_p_be,
                  "min_p_be" = min_p_be,
                  "sd_p_be" = sd_p_be,
                  "mean_lfc_be" = mean_lfc_be,
                  "min_lfc_be" = min_lfc_be,
                  "max_lfc_be" = max_lfc_be,
                  "mean_lfc_ct" = mean_lfc_ct,
                  "min_lfc_ct" = min_lfc_ct,
                  "max_lfc_ct" = max_lfc_ct,
                  "rel_be" = rel_be,
                  "rel_ct" = rel_ct,
                  "mean_p_ct"= mean_p_ct,
                  "min_p_ct"= min_p_ct,
                  "max_p_ct"= max_p_ct,
                  "sd_p_ct" = sd_p_ct,
                  "mean_lfc_batch" = mean_lfc_batch,
                  "sd_lfc_batch" = sd_lfc_batch,
                  "n_batch_genes" = n_batch_genes,
                  "mean_lfc_type" = mean_lfc_type,
                  "sd_lfc_type" = sd_lfc_type,
                  "n_type" = n_type,
                  "lfc_ct_shape" = lfc_ct_shape,
                  "lfc_ct_scale" = lfc_ct_scale,
                  "lfc_be_shape" = lfc_be_shape,
                  "lfc_be_scale" = lfc_be_scale)
summary <- cbind(size, celltype_spec, gene, sim) %>% set_rownames(dataset_name)
saveRDS(summary, file=args[8])
