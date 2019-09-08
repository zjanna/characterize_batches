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

x<-sce
expr <- as.matrix(assays(sce)$logcounts)
clust <- as.factor(colData(sce)[,celltype])
kids <- levels(clust)
names(kids) <- kids
group <- as.factor(colData(sce)[,batch])
cs <- names(cont)
names(cs) <- cs
es <- expr
ctype <- "contrast"
res_df <- function(k, tt, ct, c) {
  df <- data.frame(
    gene = rownames(tt), cluster_id = k, tt,
    row.names = NULL, stringsAsFactors = FALSE)
  df[[ct]] <- c
  return(df)
}
doDE <- function(x,lfc_cutoff=log2(1.1)){
  res <- lapply(kids, function (k) {
    cat(k, "..", sep = "")
    n <- clust==k
    es_tmp <- es[,n]
    grp <- group[n]
    design <- model.matrix(~0+grp)
    colnames(design)<-levels(group)
    k1 <- rowSums(es_tmp > 0) >= .2*min(table(grp))
    es_tmp <- es_tmp[k1,]
    f <- lmFit(es_tmp, design)
    f <- eBayes(f, trend = TRUE)
    tt <- lapply(cont, function(c) {
      cc<-names(c)
      fc <- contrasts.fit(f, contrasts = c)
      tr <- treat(fc, lfc=lfc_cutoff)
      tt <- topTreat(tr, n=Inf)
      res_df(k, tt, ctype,cc)
    })
    return(list(tt = tt, data = es_tmp))
  })
  # remove empty clusters
  skipped <- vapply(res, is.null, logical(1))
  if (any(skipped))
    message(paste("Cluster(s)", dQuote(kids[skipped]), "skipped due to an",
                  "insufficient number of cells in at least 2 samples per group."))
  res <- res[!skipped]
  kids <- kids[names(res)]
  
  # re-organize by contrast &
  # do global p-value adjustment
  tt <- lapply(res, "[[", "tt")
  tt <- lapply(cs, function(c) map(tt, c))
  
  # return results
  data <- lapply(res, "[[", "data")
  list(table = tt,
       data = data,
       design = design,
       coef = coef)
}
res <- doDE(lfc_cutoff = log2(1.1))

saveRDS(res, file=args[3])

FilterDEGs<-function (degDF=df, filter=c(FDR=5))
{
  rownames(degDF)<-degDF$gene
  pval <- degDF[, grep("adj.P.Val$", colnames(degDF)), drop = FALSE]
  pf <- pval <= filter["FDR"]/100
  pf[is.na(pf)] <- FALSE
  DEGlistUPorDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[, x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
}

result<-list()
m2<-list()

for(jj in 1:length(cs)){
  result[[jj]]<-sapply(res[[1]][[names(cs)[jj]]], function(x) FilterDEGs(x))
  names(result[[jj]])<-kids
  m2[[jj]] = make_comb_mat(result[[jj]], mode = "intersect")
}
names(result)<-names(cs)

# DE genes 
# Number of sig DE genes per celltype and mean of all celltypes.
n_de <- lapply(res[[1]],function(y) vapply(y, function(x) sum(x$adj.P.Val < 0.05), numeric(1)))
n_de_cl <- lapply(res[[1]],function(y) vapply(y, function(x) nrow(x), numeric(1)))
mean_n_de <- lapply(n_de,function(x) mean(x))
mean_mean_n_de <- mean(unlist(mean_n_de))/n_genes
min_mean_n_de <- min(unlist(mean_n_de))/n_genes
max_mean_n_de <- max(unlist(mean_n_de))/n_genes

# Genes with lfc > 1
#Number of genes with a LFC > 1 per celltype and mean of all celltypes.
n_genes_lfc1 <- lapply(res[[1]],function(y) vapply(y, function(x) sum(abs(x$logFC) > 1), numeric(1)))
mean_n_genes_lfc1 <- mean(unlist(n_genes_lfc1))/n_genes
min_n_genes_lfc1 <- min(unlist(n_genes_lfc1))/n_genes
max_n_genes_lfc1 <- max(unlist(n_genes_lfc1))/n_genes

# DE genes overlap between celltypes (celltype specific de genes)
# Genes are "overlapping" if they are present in all clusters with at least 10% of all cells
de_overlap <- lapply(result,function(x){
  result2 <- x[table(colData(sce)[, celltype]) > n_cells * 0.1]
  de_overlap <- length(Reduce(intersect, result2))
  de_overlap
})

mean_de_overlap <- mean(unlist(de_overlap))/n_genes
min_de_overlap <- min(unlist(de_overlap))/n_genes
max_de_overlap <- max(unlist(de_overlap))/n_genes

#Genes unique to single celltypes
unique_genes_matrix <- NULL
unique_genes <- NULL
cb <- length(names(result[[1]]))
unique_genes <- lapply(result,function(x){
  for( i in 1:cb ){
    unique_genes[i] <-as.numeric(length(setdiff(unlist(x[i]),unlist(x[-i]))))
  }
  unique_genes_matrix <- cbind(unique_genes_matrix, unique_genes)
  unique_genes_matrix
})
unique_genes <- Reduce('cbind', unique_genes)
colnames(unique_genes) <- names(result)
rownames(unique_genes) <- names(result[[1]])

# Relative cluster specificity (unique/overlapping)
rel_spec1 <- NULL
for( i in 1:dim(unique_genes)[2] ){
  rel_spec <- unique_genes[,i]/de_overlap[[i]]
  rel_spec1 <- cbind(rel_spec1,rel_spec)
}

mean_rel_spec <- mean(rel_spec1)
min_rel_spec <- min(rel_spec1)
max_rel_spec <- max(rel_spec1)

de_summary<-list(n_de_cl=n_de_cl, de_overlap=de_overlap,result=result,n_de=n_de,mean_n_de=mean_n_de,
                 mean_mean_n_de=mean_mean_n_de,min_mean_n_de=min_mean_n_de,
                 max_mean_n_de=max_mean_n_de,n_genes_lfc1=n_genes_lfc1,
                 mean_n_genes_lfc1=mean_n_genes_lfc1,min_n_genes_lfc1=min_n_genes_lfc1,
                 max_n_genes_lfc1=max_n_genes_lfc1,mean_de_overlap=mean_de_overlap,
                 min_de_overlap=min_de_overlap,max_de_overlap=max_de_overlap,
                 mean_rel_spec=mean_rel_spec,min_rel_spec=min_rel_spec,max_rel_spec=max_rel_spec)
saveRDS(de_summary,file=args[4])
