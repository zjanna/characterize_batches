---
title: "Preprocessing of pancreas dataset"
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
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(umap)
})

```

```{r message = FALSE}
# Parameters
# data_path <- "~/Documents/Zurich/single_cell/Almut/data/pancreas_v3_files/"
#out_path<- "/home/zjanna/preprocessing_codes/output"
#meta_dir <- "20190312 NovaSeqRun Samples Characteristics.xlsx"
seed <- 1234
```

# Load data
```{r message = FALSE}
pancreas.data <- readRDS(file = "~/Documents/Zurich/single_cell/Almut/data/pancreas_v3_files/pancreas_expression_matrix.rds")
sce<-pancreas.data
# load metadata
metadata <- readRDS(file = "~/Documents/Zurich/single_cell/Almut/data/pancreas_v3_files/pancreas_metadata.rds")
```

# Integration
```{r integration}
# create SeuratObject
pancreas <- CreateSeuratObject(pancreas.data, meta.data = metadata)
pancreas.list <- SplitObject(pancreas, split.by = "tech")
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
                                             verbose = FALSE)
}
reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]



# find anchors & integrate
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)


# scale integrated data
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)

```

# Dimension reduction
```{r dimred}
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunTSNE(pancreas.integrated, reduction = "pca", dims = 1:30)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = seq_len(30),
                  seed.use = seed, verbose = FALSE,n.neighbors = 30, min.dist = 0.5)

p1 <- DimPlot(pancreas.integrated, reduction = "tsne", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "tsne", group.by = "celltype", label = TRUE, repel = TRUE) +
  NoLegend()
plot_grid(p1, p2)


```

# Convert seurat to sce

```{r saving}
seurat<-pancreas.integrated
sce <- SingleCellExperiment(
  assays=list(
    counts=seurat@assays$RNA@counts,
    logcounts=seurat@assays$RNA@data
  ),
  colData=seurat@meta.data,
  reducedDims=lapply(seurat@reductions, FUN=function(x) x@cell.embeddings)
)
```



# Saving data
```{r saving data}
# Save data
saveRDS(sce, file = paste0(out_path, "/sce_pancreas.rds"))
saveRDS(seurat, file = paste0(out_path, "/seurat_pancreas.rds"))

```