#!/usr/bin/env Rscript
args <- (commandArgs(trailingOnly = TRUE))


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
n_genes <- nrow(sce)
n_cells <- ncol(sce)
cols <-c(c(jcolors('pal6'),jcolors('pal8'))[c(1,8,14,5,2:4,6,7,9:13,15:20)],jcolors('pal4'))
names(cols) <- c()
# dataset_name <- substr(args[1] , start = 6 , stop = 17 )
if(grepl('pancreas', args[1])){
  dataset_name<-'pancreas'
  # param
  MultiSample = FALSE
  #variables
  batch <- "tech"
  patient=FALSE
  celltype <- "celltype"
  sample <- NA
  #contrast
  cont<-list("celseq-celseq2" = c(1,-1,0),
             "celseq-smartseq2" = c(1,0,-1),
             "celseq2-smartseq2" = c(0,1,-1))
  k = 120; cell_min = 10; n_dim = 10
  params<-list(n_genes=n_genes, n_cells=n_cells, dataset_name=dataset_name,cols=cols,MultiSample=MultiSample,
               batch=batch,celltype=celltype,sample=sample, patient=patient,cont=cont, k=k, cell_min=cell_min,
               n_dim = n_dim)
  saveRDS(params, file=args[2])
}

if(grepl('cellBench', args[1])){
  dataset_name<-'cell_bench'
  # param
  MultiSample = FALSE
  #variables
  batch <- "protocol"
  patient=FALSE
  celltype <- "cell_line"
  sample <- NA
  #contrast
  cont <-list("sc_10x-sc_celseq" = c(1,-1,0),
              "sc_10x-sc_dropseq" = c(1,0,-1),
              "sc_celseq-sc_dropseq" = c(0,1,-1))
  # cms param
  k = 120; cell_min = 10; n_dim = 10
  params<-list(n_genes=n_genes, n_cells=n_cells, dataset_name=dataset_name,cols=cols,MultiSample=MultiSample,
               batch=batch,celltype=celltype,sample=sample, patient=patient,cont=cont, k=k, cell_min=cell_min,
               n_dim = n_dim)
  saveRDS(params, file=args[2])
}

if(grepl('csf_media', args[1])){
  dataset_name<-'csf_media'
  # param
  #param
  MultiSample = TRUE
  #variables
  batch <- "media"
  patient=FALSE
  celltype <- "cluster.merged"
  sample <- "Sample"
  #contrast
  cont <- list("cryo-fresh" = c(1,-1))
  # cms param
  k = 120; cell_min = 10; n_dim = 10
  params<-list(n_genes=n_genes, n_cells=n_cells, dataset_name=dataset_name,cols=cols,MultiSample=MultiSample,
               batch=batch,celltype=celltype,sample=sample, patient=patient,cont=cont, k=k, cell_min=cell_min,
               n_dim = n_dim)
  saveRDS(params, file=args[2])
}

if(grepl('csf_patient', args[1])){
  dataset_name<-'csf_patient'
  #param
  MultiSample = FALSE
  #variables
  batch <- "patient"
  patient=TRUE
  celltype <- "cluster.merged"
  sample <- "Sample"
  #contrast
  cont <- list("pat2-pat3" = c(-1,1))
  # cms param
  k = 120; cell_min = 10; n_dim = 10
  params<-list(n_genes=n_genes, n_cells=n_cells, dataset_name=dataset_name,cols=cols,MultiSample=MultiSample,
               batch=batch,celltype=celltype,sample=sample, patient=patient,cont=cont, k=k, cell_min=cell_min,
               n_dim = n_dim)
  saveRDS(params, file=args[2])
}

if(grepl('hca', args[1])){
  dataset_name<-'hca'
  #param
  MultiSample = FALSE
  #variables
  batch <- "protocol"
  patient=TRUE
  celltype <- "integrated_snn_res.0.2"
  sample <- NA
  #contrast
  cont<-list("10X2x5Kcell250Kreads-CELseq2" = c(1,-1,0,0),
             "10X2x5Kcell250Kreads-Dropseq" = c(1,0,-1,0),
             "10X2x5Kcell250Kreads-SMARTseq2" = c(1,0,0,-1))
  k = 300; cell_min = 30; n_dim = 10
  params<-list(n_genes=n_genes, n_cells=n_cells, dataset_name=dataset_name,cols=cols,MultiSample=MultiSample,
               batch=batch,celltype=celltype,sample=sample, patient=patient,cont=cont, k=k, cell_min=cell_min,
               n_dim = n_dim)
  saveRDS(params, file=args[2])
}

if(grepl('khang_patient', args[1])){
  dataset_name<-'khang_patient'
  #param
  MultiSample = FALSE
  #variables
  batch <- "patient"
  patient=TRUE
  celltype <- "cell"
  sample <- NA
  #contrast
  cont <- list("1015-1016" = c(0,1,-1,0,0,0,0,0),
               "1015-1244" = c(0,1,0,0,0,-1,0,0),
               "1015-1256" = c(0,1,0,0,0,0,-1,0),
               "1015-1488" = c(0,1,0,0,0,0,0,-1),
               "1016-1244" = c(0,0,1,0,0,-1,0,0),
               "1016-1256" = c(0,0,1,0,0,0,-1,0),
               "1244-1256" = c(0,0,0,0,0,1,-1,0),
               "1244-1488" = c(0,0,0,0,0,1,0,-1),
               "1256-1488" = c(0,0,0,0,0,0,1,-1))
  k = 120; cell_min = 10; n_dim = 10
  params<-list(n_genes=n_genes, n_cells=n_cells, dataset_name=dataset_name,cols=cols,MultiSample=MultiSample,
               batch=batch,celltype=celltype,sample=sample, patient=patient,cont=cont, k=k, cell_min=cell_min,
               n_dim = n_dim)
  saveRDS(params, file=args[2])
}

if(grepl('pbmc_media_storage', args[1])){
  dataset_name<-'pbmc_media_storage'
  #param
  MultiSample = FALSE
  #variables
  batch <- "sample_type"
  patient=FALSE
  celltype <- "cluster.merged"
  sample <- NA
  #contrast
  cont <- list("DMSO-Fresh" = c(0,1,-1,0),
               "PSC-Fresh" = c(0,0,-1,1),
               "CS10-Fresh" = c(1,0,-1,0))
  k = 120; cell_min = 10; n_dim = 10
  params<-list(n_genes=n_genes, n_cells=n_cells, dataset_name=dataset_name,cols=cols,MultiSample=MultiSample,
               batch=batch,celltype=celltype,sample=sample, patient=patient,cont=cont, k=k, cell_min=cell_min,
               n_dim = n_dim)
  saveRDS(params, file=args[2])
}

if(grepl('pbmc2_patient', args[1])){
  dataset_name<-'pbmc2_patient'
  #param
  MultiSample = TRUE
  # variables
  batch <- 'patient'
  patient=TRUE
  celltype <- 'cluster.merged'
  sample <- 'Sample'
  # contrast
  cont <- list('pat1-pat2' = c(1,-1))
  k = 120; cell_min = 10; n_dim = 10
  params<-list(n_genes=n_genes, n_cells=n_cells, dataset_name=dataset_name,cols=cols,MultiSample=MultiSample,
  batch=batch,celltype=celltype,sample=sample, patient=patient,cont=cont, k=k, cell_min=cell_min,
  n_dim = n_dim)
  saveRDS(params, file=args[2])
}

if(grepl('pbmc2_media', args[1])){
  dataset_name<-'pbmc2_media'
  #param
  MultiSample = TRUE
  #variables
  batch <- 'patient'
  patient=TRUE
  celltype <- "cluster.merged"
  sample <- "Sample"
  #contrast
  cont <- list("DMSO-fresh" = c(1,-1,0),
             "MetOH-Fresh" = c(0,-1,1),
             "DMSO-MetOH" = c(1,0,-1))
  k = 120; cell_min = 10; n_dim = 10
  params<-list(n_genes=n_genes, n_cells=n_cells, dataset_name=dataset_name,cols=cols,MultiSample=MultiSample,
               batch=batch,celltype=celltype,sample=sample, patient=patient,cont=cont, k=k, cell_min=cell_min,
               n_dim = n_dim)
  saveRDS(params, file=args[2])
}

