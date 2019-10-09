---
title: "Preprocessing of cellBench dataset"
date: '`r format(Sys.Date(), "%B %d, %Y")`'
output:
  html_document:
  toc: true
---
  
```{r echo = FALSE, warning = FALSE, message = FALSE}
knitr::opts_chunk$set(autodep = TRUE, cache = TRUE)
```

# Libraries

```{r message = FALSE, warning = FALSE}
suppressPackageStartupMessages({
  library(CellBench)
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
  library(scran)
  library(cowplot)
})

```

```{r message = FALSE}
# Parameters
data_path <- here::here("data")
out_path <- here::here("output")
code_path <- here::here("code")
dataset_name <- "cellBench"
seed <- 1234
```

# Load data
```{r message = FALSE}
sc_data <- load_sc_data()
colData(sc_data[[1]])$protocol <- rep(names(sc_data)[1], ncol(sc_data[[1]]))
sce <- sc_data[[1]]

for(i in 2:length(sc_data)){
  colData(sc_data[[i]])$protocol <- rep(names(sc_data)[i], ncol(sc_data[[i]]))
  gene_overlap <- intersect(rownames(sce), rownames(sc_data[[i]]))
  coldata_overlap <- intersect(names(colData(sce)), names(colData(sc_data[[i]])))
  sc_data[[i]] <- sc_data[[i]][gene_overlap,]
  colData(sc_data[[i]]) <- colData(sc_data[[i]])[, coldata_overlap]
  colData(sce) <- colData(sce)[, coldata_overlap]
  sce <- sce[gene_overlap,]
  sce <- cbind(sce, sc_data[[i]])
}

colnames(sce) <- paste0(colnames(sce), "_", sce$protocol)
dim(sce)
## Filter out genes that are not expressed in any cell
sce <- sce[which(rowSums(counts(sce) > 0) > 0), ]
dim(sce)
```


# Normalization
```{r normalization}
# Scater
set.seed(1000)
clusters <- quickCluster(sce, use.ranks=FALSE)
table(clusters)
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters) ##cluster information added
sce <- normalize(sce)   
```


# Dimension reduction
```{r dimred}
sce <- RunPCA(sce, npcs = 20, verbose = FALSE)
sce <- RunTSNE(sce, perplexity = 30,reduction = "pca", dims = seq_len(20),
                  seed.use = seed, do.fast = TRUE, verbose = FALSE)
sce <- RunUMAP(sce, reduction = "pca", dims = seq_len(20),
                  seed.use = seed, verbose = FALSE,n.neighbors = 30, min.dist = 0.5)
```

# Saving data

```{r saving}
# Save data
saveRDS(sce, file = paste0(out_path, "/sce_cellBench.rds"))
saveRDS(seurat, file = paste0(out_path, "/seurat_cellBench.rds"))

```