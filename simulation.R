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
  library(TruncatedDistributions)
  library(scran)
  library(ComplexHeatmap)
})
# Simulation parameter
res<-readRDS(file=args[1])
de_summary<-readRDS(file=args[2])
attach(de_summary)
sce<-readRDS(file=args[3])
params_summary_vp<-readRDS(file=args[4])
attach(params_summary_vp)
params<-readRDS(file=args[5])
attach(params)

# Simulation parameter
# Define parameter to describe the dataset for simulation
block <- colData(sce)[,batch]
kids <- colData(sce)[,celltype]
res_ct <- pairwiseTTests(
  assays(sce)$logcounts,
  clusters = kids, block = block,
  direction = "any", lfc = 0)

tbl <- res_ct$statistics %>% map_depth(1, data.frame)
names(tbl) <- names(tbl) <- res_ct$pairs$first

lfc_ct_all <- tbl %>% map(function(x){
  lfc <- x$logFC
}) %>% unlist()

#variance within celltypes
lfc_by_k <- tbl %>% map(function(x){
  lfc <- x$logFC
}) %>% bind_cols()

# all lfcs merged together
mean_lfc_tib <- tibble("mean_lfc" = colMeans(lfc_by_k), "cluster" = names(tbl)) %>% 
  group_by(cluster) %>% summarize("mean_k" = mean(mean_lfc)) 

lfc_ct_real <- data.frame("lfc" = lfc_ct_all, "type" = rep("real", length(lfc_ct_all)))

#distribution parameter 
mean_lfc_ct <- mean(lfc_ct_all)
min_lfc_ct <- min(lfc_ct_all)
max_lfc_ct <- max(lfc_ct_all)

###simulate overall celltype lfc distribution

#parameter
lfc_ct_shape <- 0.6
lfc_ct_scale <- 0.35

#Gamma distribution that should fit the real lfc distribution 
signs <- sample(c(-1, 1), length(lfc_ct_all), TRUE, prob = c(0.5, 0.5))
tgam_ct <- rtgamma(nrow(lfc_ct_real), lfc_ct_shape, lfc_ct_scale, a=min_lfc_ct, b=max_lfc_ct)*signs
lfc_sim_ct <- data.frame("lfc" = tgam_ct, "type" = rep("sim", nrow(lfc_ct_real)))


##### "type genes" distribution 
type_gene <- which(lfc_ct_all > max(tgam_ct) | lfc_ct_all < min(tgam_ct))

# Number of "type genes"
length(type_gene)

#Distribution of "type genes"
lfc_type <- data.frame("lfc" = lfc_ct_all[type_gene], "type" = rep("type", length(type_gene)))

type_by_ct <- colnames(lfc_by_k) %>% map(function(x){
  v <- lfc_by_k[,x]
  type <- length(which(v > max(tgam_ct) | v < min(tgam_ct)))
}) %>% unlist() %>% set_names(names(tbl))

#Percentage of "type genes" per cluster
type_var <- tibble("type_gene" = type_by_ct, "cluster" = as.numeric(names(type_by_ct))) %>% group_by(cluster) %>% summarise("n_type" = sum(type_gene)) %>% mutate("p" = n_type/(nrow(lfc_by_k) * ncol(lfc_by_k)/nrow(.)))


#Plot real_lfc distribution, fitted gamma distribution and type gene distribution
lfc_comb <- rbind(lfc_sim_ct, lfc_ct_real, lfc_type)
# Extract parameter to simulate "type genes"
mean_lfc_type <- mean(abs(lfc_ct_all[type_gene]))
sd_lfc_type <- sd(abs(lfc_ct_all[type_gene]))
n_type <- length(type_gene)/length(lfc_ct_all)


mean_lfc_cl <- lapply(res[[1]], function(y) vapply(y, function(x){
  de_genes <- which(x$adj.P.Val < 0.05)
  mean_de <- mean(abs(x[, "logFC"]))}
  , numeric(1))) %>% bind_cols()

min_lfc_cl <- lapply(res[[1]],function(y) vapply(y, function(x){
  de_genes <- which(x$adj.P.Val < 0.05)
  min_de <- min(x[, "logFC"])}
  , numeric(1))) %>% bind_cols()

max_lfc_cl <- lapply(res[[1]],function(y) vapply(y, function(x){
  de_genes <- which(x$adj.P.Val < 0.05)
  max_de <- max(x[, "logFC"])}
  , numeric(1))) %>% bind_cols()

mean_lfc_be <- mean(colMeans(mean_lfc_cl, na.rm = TRUE))
min_lfc_be <- min(colMins(as.matrix(min_lfc_cl), na.rm = TRUE))
max_lfc_be <- max(colMaxs(as.matrix(max_lfc_cl), na.rm = TRUE))


# Overall lfc distribution of the batch effect
lfc_dist <- lapply(res[[1]][[1]],function(x){
  mean_de <- x[, "logFC"]}) %>% unlist() %>% as.data.frame() %>% set_colnames("test") %>% mutate("type" = rep("real", nrow(.)))

#Simulation parameter
lfc_be_shape <- 1
lfc_be_scale <- 0.15

#Fit a gamma distribution describimng the overall lfc distribution
signs <- sample(c(-1, 1), nrow(lfc_dist), TRUE, prob = c(0.37, 0.63))
tgam <- rtgamma(nrow(lfc_dist), lfc_be_shape, scale=lfc_be_scale, a=min_lfc_be, b=max_lfc_be)*signs
lfc_sim <- data.frame("test" = tgam, "type" = rep("sim", nrow(lfc_dist)))

#"batch gene" distribution
batch_gene <- which(lfc_dist$test > max(tgam) | lfc_dist$test < min(tgam))
length(batch_gene)
lfc_batch <- data.frame("test" = lfc_dist$test[batch_gene], "type" = rep("batch", length(batch_gene)))


mean_lfc_batch <- mean(abs(lfc_dist$test[batch_gene]))
sd_lfc_batch <- sd(abs(lfc_dist$test[batch_gene]))
n_batch_genes <- length(batch_gene)/nrow(lfc_dist)

#Fit a "batch gene distribution" assuming normality
signs <- sample(c(-1, 1), length(batch_gene), TRUE, prob = c(0.8, 0.2))
tnorm <- rnorm(length(batch_gene), mean = mean_lfc_batch, sd = sd_lfc_batch) * signs
lfc_batch_sim <- data.frame("test" = tnorm, "type" = rep("batch_sim", length(batch_gene)))

# Combine and plot overall and batch gene distribution, and fitted distributions for both
lfc_comb <- rbind(lfc_sim, lfc_dist, lfc_batch, lfc_batch_sim)


#percentage of batch affected genes
cond <- gsub("-.*", "", names(n_de))
cond <- c(cond, unique(gsub(".*-", "", names(n_de))))
cond <- unique(cond)

de_be_tab <- n_de %>% bind_cols()
de_cl_tab <- n_de_cl %>% bind_cols()
de_be <- cond %>% map(function(x){
  de_tab <- de_be_tab[, grep(x, colnames(de_be_tab))]
  de_be <- rowMeans(de_tab)
}) %>% bind_cols() %>% set_colnames(cond)

n_cl <- cond %>% map(function(x){
  cl_tab <- de_cl_tab[, grep(x, colnames(de_cl_tab))]
  de_cl <- rowMeans(cl_tab)
}) %>% bind_cols() %>% set_colnames(cond)


p_be <- de_be/n_cl
mean_p_be <- mean(colMeans(p_be))
min_p_be <- min(colMins(as.matrix(p_be)))
max_p_be <- max(colMaxs(as.matrix(p_be)))
sd_p_be <- mean(colSds(as.matrix(p_be)))
if(is.na(sd_p_be)){ sd_p_be <- 0 }


#### How much does the batch fold change vary between celltypes? "rel_be"
var_lfc_cl <- colSds(as.matrix(mean_lfc_cl), na.rm = TRUE)
rel_be <- mean(var_lfc_cl)

#### How much does the celltype fold change vary between celltypes? "rel_ct"
rel_ct <- sd(mean_lfc_tib$mean_k)

#### Percentage of celltype specific genes "p_ct"
n_de_unique <- lapply(result,function(x){
  de_genes <- unlist(x) %>% unique() %>% length()
  de_genes <- de_genes/length(x)
}) %>% bind_cols()


rel_spec2 <- NULL
for(i in 1:length(de_overlap)){
  rel_spec <- de_overlap[[i]]/mean(n_de[[i]][table(colData(sce)[,celltype])>dim(expr)[2]*0.1])
  rel_spec2 <- cbind(rel_spec2,rel_spec)
}

mean_p_ct <- 1 - mean(rel_spec2)
max_p_ct <- 1 - min(rel_spec2)
min_p_ct <- 1 - max(rel_spec2)
sd_p_ct <- sd(rel_spec2)
if(is.na(sd_p_ct)){ sd_p_ct <- 0 }

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

saveRDS(sim, file=args[6])
