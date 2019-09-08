######################################## cellBench ########################################
# data 
data_path <- here::here("data")
out_path <- here::here("output")
code_path <- here::here("code")
dataset_name <- "cellBench"

#Mixology dataset
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

sce <- computeSumFactors(sce)
sce <- normalize(sce)
sce <- runPCA(sce, ncomponents = 20)
sce <- runTSNE(sce)
sce <- runUMAP(sce)
saveRDS(sce,file=paste0(out_path,'/cellBench.rds'))
######################################## csf_media_storage ########################################
# data 
data_path <- here::here("data")
out_path <- here::here("output")
code_path <- here::here("code")
dataset_name <- "csf_media_storage"

# loading the data
sce <- readRDS(paste0(data_path, "/sce_csf_Integrated.rds"))
sce <- computeSumFactors(sce)
sce <- normalize(sce)
sce <- runUMAP(sce)
names(colData(sce))[48]<-'cluster.merged'
saveRDS(sce,file=paste0(out_path,'/csf_media.rds'))
######################################## csf_patient ########################################
# data 
data_path <- here::here("data")
out_path <- here::here("output")
code_path <- here::here("code")
dataset_name <- "csf_patient"

# loading the data
sce <- readRDS(paste0(data_path, "/sce_csf_Integrated.rds"))
sce$patient <- ifelse(grepl("E2490", sce$Sample), "pat1", ifelse(grepl("E2528", sce$Sample),"pat2","pat3"))

#Filter patient 2 and 3 (we have only one sample of patient1)
sce <- sce[, sce$patient %in% c("pat2", "pat3")]
sce <- computeSumFactors(sce)
sce <- normalize(sce)
sce <- runUMAP(sce)

names(colData(sce))[which(names(colData(sce)) %in% "integrated_snn_res.0.2")] <- 'cluster.merged'

######################################## hca ########################################
# data
data_path <- here::here("data")
out_path <- here::here("output")
code_path <- here::here("code")
dataset_name <- "hca"


# sce <- readRDS("integrated_protocols_mouse_red.rds")


# 
# sce <- sce[which(rowSums(as.matrix(assays(sce)$counts)>1)>=20),] 
# seurat <- CreateSeuratObject( counts=as.matrix(assays(sce)$counts),
#                               meta.data=as.data.frame(colData(sce)),
#                               project = "10X" )
# 
# seurat.list <- SplitObject(seurat, split.by = "protocol")
# for (i in 1:length(seurat.list)) {
#   seurat.list[[i]] <- NormalizeData(seurat.list[[i]], verbose = FALSE)
#   seurat.list[[i]] <- FindVariableFeatures(seurat.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
# }
# reference.list <- seurat.list[c("10X2x5Kcell250Kreads",'CELseq2', 'Dropseq', 'SMARTseq2' )]
# seurat.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
# 
# seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:30)
# DefaultAssay(seurat.integrated) <- "integrated"
# seurat <- ScaleData(object = seurat.integrated, verbose = FALSE)
# seurat <- RunPCA(object = seurat, npcs = 20, verbose = FALSE)
# seurat <- FindNeighbors(object = seurat, reduction = "pca", dims = 1:20)
# seurat <- FindClusters(object = seurat, resolution = 0.2, random.seed = 1234)
# seurat <- RunTSNE(object = seurat, reduction = "pca", dims = 1:20)
# seurat <- RunUMAP(object = seurat, reduction = "pca", dims = 1:10, n.neighbors = 30, min.dist = 0.5)
# 
# sce <- SingleCellExperiment(
#           assays=list(
#             counts=seurat@assays$RNA@counts,
#             logcounts=seurat@assays$RNA@data
#           ),
#           colData=seurat@meta.data,
#           reducedDims=lapply(seurat@reductions, FUN=function(x) x@cell.embeddings)
#     )
# sce <- runUMAP(sce)
# saveRDS(sce,file='sce_hca_integrated.rds')
sce <-readRDS(paste0(data_path, '/sce_hca_integrated.rds'))
sce <- computeSumFactors(sce)
sce <- normalize(sce)
sce <- runUMAP(sce)
saveRDS(sce, file=paste0(out_path,'/hca.rds'))
######################################## khang_patient ########################################

# data
data_path <- here::here("data")
out_path <- here::here("output")
code_path <- here::here("code")
dataset_name <- "khang_patient"

# loading the data
sce <- readRDS(paste0(data_path, "/sce_khang_patient.rds"))
sce <- computeSumFactors(sce)
sce <- normalize(sce)
sce <- runUMAP(sce)
sce$patient <- as.factor(as.character(colData(sce)[,'ind']))
saveRDS(sce, file=paste0(out_path,'/khang.rds'))
######################################## pancreas ########################################

# data
data_path <- here::here("data")
out_path <- here::here("output")
code_path <- here::here("code")
dataset_name <- "pancreas"

sce <- readRDS(paste0(data_path,"/pancreas.rds"))
sce <- computeSumFactors(sce)
sce <- normalize(sce)
sce <- runUMAP(sce)
saveRDS(sce, file=paste0(out_path,'/pancreas.rds'))
######################################## pbmc_media_storage ########################################

# data 
data_path <- here::here("data")
out_path <- here::here("output")
code_path <- here::here("code")
dataset_name <- "pbmc_media_storage"
# raw data 
#For preprocessing and quality control please check: /code/pbmc_media.rmd
sce <- readRDS(paste0(data_path, "/pbmc_media.rds"))
sce <- computeSumFactors(sce)
sce <- normalize(sce)
sce <- runUMAP(sce)


colData(sce)[,'sample_type'] <- gsub("PBMCs frozen in ","",
                                     colData(sce)[,'sample_type'])
colData(sce)[,'sample_type'] <- gsub(" 7 ?days","",
                                     colData(sce)[,'sample_type'])
colData(sce)[,'sample_type'] <- gsub(" PBMCs","",
                                     colData(sce)[,'sample_type'])

saveRDS(sce, file=paste0(out_path,'/pbmc_media_storage.rds'))
######################################## pbmc2_media_storage ########################################

# data 
data_path <- here::here("data")
out_path <- here::here("output")
code_path <- here::here("code")
dataset_name <- "pbmc2_media_storage"

sce <- readRDS(paste0(data_path, "/pbmc2_media.rds"))
sce <- sce[,!sce$Sample %in% "CR057"]
sce <- computeSumFactors(sce)
sce <- normalize(sce)
sce <- runUMAP(sce)

names(colData(sce))[50]<-'cluster.merged'

######################################## pbmc2_patient ########################################

# data 
data_path <- here::here("data")
out_path <- here::here("output")
code_path <- here::here("code")
dataset_name <- "pbmc2_patient"

sce <- readRDS(paste0(data_path, "/pbmc2_media.rds"))
sce <- computeSumFactors(sce)
sce <- normalize(sce)
sce <- runUMAP(sce)


names(colData(sce))[50]<-'cluster.merged'