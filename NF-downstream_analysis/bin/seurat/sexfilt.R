#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer",
  'custom_functions', 'm', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set run location
if(length(commandArgs(trailingOnly = TRUE)) == 0){
  cat('No command line arguments provided, user defaults paths are set for running interactively in Rstudio on docker\n')
  opt$runtype = "user"
} else {
  if(is.null(opt$runtype)){
    stop("--runtype must be either 'user' or 'nextflow'")
  }
  if(tolower(opt$runtype) != "user" & tolower(opt$runtype) != "nextflow"){
    stop("--runtype must be either 'user' or 'nextflow'")
  }
  if(tolower(opt$runtype) == "nextflow"){
    if(is.null(opt$custom_functions) | opt$custom_functions == "null"){
      stop("--custom_functions path must be specified in process params config")
    }
  }
}

# Set paths and load data
{
  if (opt$runtype == "user"){
    sapply(list.files('./NF-downstream_analysis/bin/custom_functions/', full.names = T), source)
    plot_path = "./output/NF-downstream_analysis/sexfilt/plots/"
    rds_path = "./output/NF-downstream_analysis/sexfilt/rds_files/"

    data_path = "./output/NF-downstream_analysis/integration_qc/rds_files/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  
  # Load packages - packages are stored within renv in the repository
  reticulate::use_python('/usr/bin/python3.7')
  library(Seurat)
  library(sctransform)
  
  library(future)
  library(dplyr)
  library(cowplot)
  library(clustree)
  library(gridExtra)
  library(grid)
  library(pheatmap)
  library(RColorBrewer)
  library(tidyverse)
}

pre_sexfilt_data <- readRDS(paste0(data_path, 'integration_qc_data.RDS'))

DefaultAssay(pre_sexfilt_data) <- "RNA"

# Log normalize data and find variable features
pre_sexfilt_data <- NormalizeData(pre_sexfilt_data, normalization.method = "LogNormalize", scale.factor = 10000)
pre_sexfilt_data <- FindVariableFeatures(pre_sexfilt_data, selection.method = "vst", nfeatures = 2000)

# Multi-core when running from command line
if(opt$runtype == "nextflow"){
  plan("multiprocess", workers = ncores)
  options(future.globals.maxSize = 32* 1024^3) # 32gb
}

pre_sexfilt_data <- ScaleData(pre_sexfilt_data, features = rownames(pre_sexfilt_data), vars.to.regress = "percent.mt", verbose = FALSE)

# Save RDS after integration
saveRDS(pre_sexfilt_data, paste0(rds_path, "pre_sexfilt_data.RDS"))

# There is a strong sex effect - this plot shows DE genes between clusters 1 and 2 which are predominantly hh4 clusters. Clustering is driven by sex genes
png(paste0(plot_path, 'HM.top15.DE.pre-sexfilt.png'), height = 40, width = 70, units = 'cm', res = 500)
tenx.pheatmap(data = pre_sexfilt_data[,rownames(pre_sexfilt_data@meta.data[pre_sexfilt_data$seurat_clusters == 1 | pre_sexfilt_data$seurat_clusters == 2,])],
              metadata = c("seurat_clusters", "orig.ident"), selected_genes = rownames(FindMarkers(pre_sexfilt_data, ident.1 = 1, ident.2 = 2)),
              hclust_rows = T, gaps_col = "seurat_clusters")
graphics.off()

# plot dimplot for main W gene
#####################################################################################################
#     Heatmap clearly shows clusters segregate by sex - check this and regress out the sex effect   #
#####################################################################################################

# Use W chromosome genes to K-means cluster the cells into male (zz) and female (zw)
W_genes <- as.matrix(pre_sexfilt_data@assays$RNA[grepl("W-", rownames(pre_sexfilt_data@assays$RNA)),])
k_clusters <- kmeans(t(W_genes), 2)
k_clusters <- data.frame(k_clusters$cluster)
pre_sexfilt_data@meta.data$k_clusters <- k_clusters[match(colnames(pre_sexfilt_data@assays$RNA), rownames(k_clusters)),]

# Get rownames for kmeans clusters 1 and 2
k_clus_1 <- rownames(pre_sexfilt_data@meta.data[pre_sexfilt_data@meta.data$k_clusters == 1,])
k_clus_2 <- rownames(pre_sexfilt_data@meta.data[pre_sexfilt_data@meta.data$k_clusters == 2,])

# K clustering identities are stochastic, so I mneed to identify which cluster is male and female
# Sum of W genes is order of magnitude greater in cluster 2 - these are the female cells
sumclus1 <- sum(W_genes[,k_clus_1])
sumclus2 <- sum(W_genes[,k_clus_2])

if(sumclus1 < sumclus2){
  k.male <- k_clus_1
  k.female <- k_clus_2
} else {
  k.female <- k_clus_1
  k.male <- k_clus_2
}

# Add sex data to meta.data
pre_sexfilt_data@meta.data$sex <- unlist(lapply(rownames(pre_sexfilt_data@meta.data), function(x)
  if(x %in% k.male){"male"} else if(x %in% k.female){"female"} else{stop("cell sex is not assigned")}))


###### Next - subset autosomal genes and Z genes - then calculate the average for each gene for both kmeans clustered cells and plot in order to compare whether the two groups show significant differences in their expression.

# This is a test on autosomal genes to try and calculate and compare FC betweeen clusters
# Calculating median is tricky as there are a lot of dropouts in 10x data so you end up with either 0s (when the median  = 0) or 1 (when the median expression in both clusters is the same - probably a result of normalisation resulting in a UMI of 0 or 1 being normalised to a nominal value)

# Make dataframe for mean Z expression in male cells
mean.Z.male <- data.frame(Z.mean = apply(pre_sexfilt_data@assays$RNA[grepl("Z-", rownames(pre_sexfilt_data@assays$RNA)), k.male], 1, mean))
# add 1 before log2 as log2(1) = 0
mean.Z.male <- log2(mean.Z.male + 1)

# Make dataframe for mean Z expression in female cells
mean.Z.female <- data.frame(Z.mean = apply(pre_sexfilt_data@assays$RNA[grepl("Z-", rownames(pre_sexfilt_data@assays$RNA)), k.female], 1, mean))
mean.Z.female <- log2(mean.Z.female + 1)

# Make dataframe for mean autosomal expression in male cells
mean.auto.male <- data.frame(auto.mean = apply(pre_sexfilt_data@assays$RNA[!grepl("Z-", rownames(pre_sexfilt_data@assays$RNA)) & !grepl("W-", rownames(pre_sexfilt_data@assays$RNA)), k.male], 1, mean))
mean.auto.male <- log2(mean.auto.male + 1)

# Make dataframe for mean autosomal expression in male cells
mean.auto.female <- data.frame(auto.mean = apply(pre_sexfilt_data@assays$RNA[!grepl("Z-", rownames(pre_sexfilt_data@assays$RNA)) & !grepl("W-", rownames(pre_sexfilt_data@assays$RNA)), k.female], 1, mean))
mean.auto.female <- log2(mean.auto.female + 1)

# Calculate FC by subtracting log2 expression from each other
FC <- list()
FC$Z <- mean.Z.male - mean.Z.female
FC$auto <-  mean.auto.male - mean.auto.female

# Plot boxplot of Z gene and autosomal expression in male vs female cells
png(paste0(plot_path,"sex_kmeans_log2FC_boxplot.png"), height = 18, width = 18, units = "cm", res = 200)
boxplot(c(FC$Z, FC$auto),  ylab = "male - female log2 FC (mean normalised UMI +1)", names = c("Z chromosome genes", "autosomal genes"))
graphics.off()

# Z genes are upregulated within male genes relative to female genes whereas autosomal genes have a normal distribution of logFCs

#####################################################################################################
#                                     Regress sex effect                                            #
#####################################################################################################

# Init sexscale object 
sexfilt_data <- pre_sexfilt_data

DefaultAssay(sexfilt_data) <- "integrated"

# Log normalize data and find variable features
sexfilt_data <- FindVariableFeatures(sexfilt_data, selection.method = "vst", nfeatures = 2000)

# Multi-core when running from command line
if(opt$runtype == "nextflow"){
  plan("multiprocess", workers = ncores)
  options(future.globals.maxSize = 32* 1024^3) # 32gb
}

sexfilt_data <- ScaleData(sexfilt_data, features = rownames(sexfilt_data), vars.to.regress = c("percent.mt", "sex"), verbose = FALSE)

# PCA
sexfilt_data <- RunPCA(object = sexfilt_data, verbose = FALSE)

png(paste0(plot_path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(sexfilt_data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(sexfilt_data, ndims = 40))
graphics.off()

png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(sexfilt_data, PCA.levels = c(10, 20, 30, 40), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
sexfilt_data <- FindNeighbors(sexfilt_data, dims = 1:30, verbose = FALSE)
sexfilt_data <- RunUMAP(sexfilt_data, dims = 1:30, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = sexfilt_data, by = 0.1, prefix = "integrated_snn_res.")
graphics.off()

# Use clustering resolution = 0.5 to look for contamination clusters
sexfilt_data <- FindClusters(sexfilt_data, resolution = 0.5, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(sexfilt_data)
graphics.off()

# Plot QC for each cluster
png(paste0(plot_path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(sexfilt_data)
graphics.off()


# switch to RNA assay for viewing expression data
DefaultAssay(sexfilt_data) <- "RNA"

# Log normalize data and find variable features
sexfilt_data <- FindVariableFeatures(sexfilt_data, selection.method = "vst", nfeatures = 2000)

# Multi-core when running from command line
if(opt$runtype == "nextflow"){
  plan("multiprocess", workers = ncores)
  options(future.globals.maxSize = 32* 1024^3) # 32gb
}

sexfilt_data <- ScaleData(sexfilt_data, features = rownames(sexfilt_data), vars.to.regress = c("percent.mt", "sex"), verbose = FALSE)

# Save RDS
saveRDS(sexfilt_data, paste0(rds_path, "sexfilt_data.RDS"))

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(sexfilt_data, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster_order = order.cell.stage.clust(seurat_object = sexfilt_data, col.to.sort = seurat_clusters, sort.by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(plot_path, 'HM.top15.DE.post-sexfilt.png'), height = 75, width = 100, units = 'cm', res = 500)
tenx.pheatmap(data = sexfilt_data, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()

