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
    plot_path = "./output/NF-downstream_analysis/seurat_sexfilt/plots/"
    rds_path = "./output/NF-downstream_analysis/seurat_sexfilt/rds_files/"
    data_path = "./output/NF-downstream_analysis/scRNAseq/seurat_integrate/rds_files/"
    
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
  library(STACAS)
  library(tidyverse)
}


integrated_data <- readRDS(paste0(data_path, "integrated_data.RDS"))

# set integrated count data as default
# DefaultAssay(integrated_data) <- "integrated"


integrated_data <- RunPCA(object = integrated_data, verbose = FALSE)
integrated_data <- FindNeighbors(integrated_data, dims = 1:30, verbose = FALSE)
integrated_data <- RunUMAP(integrated_data, dims = 1:30, verbose = FALSE)
integrated_data <- FindClusters(integrated_data, resolution = 0.5, verbose = FALSE)

png(paste0(plot_path, "UMAP_sctransform.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(integrated_data)
graphics.off()




# Change plot path
curr_plot_path <- paste0(plot_path, '0_integrated_data/')
dir.create(curr_plot_path)


#####################################################################################################
#   Seurat's clustering algorithm is based on principle components, so we need to ensure that only the informative PCs are kept! #
#####################################################################################################

# Run PCA analysis
integrated_data <- RunPCA(object = integrated_data, verbose = FALSE)

# Plot heatmap of top variable genes across top principle components
png(paste0(curr_plot_path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(integrated_data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

# another heuristic method is ElbowPlot which ranks PCs based on the % variance explained by each PC
png(paste0(curr_plot_path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(integrated_data, ndims = 40))
graphics.off()

# Run clustering and UMAP at different PCA cutoffs - save this output to compare the optimal number of PCs to be used
png(paste0(curr_plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(integrated_data, PCA.levels = c(5, 10, 20, 40), cluster_res = 0.5)
graphics.off()

# Use PCA=20 as elbow plot is relatively stable across stages
# Use clustering resolution = 0.5 for filtering
integrated_data <- FindNeighbors(integrated_data, dims = 1:20, verbose = FALSE)
integrated_data <- RunUMAP(integrated_data, dims = 1:20, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr_plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = integrated_data, by = 0.1, prefix = 'integrated_snn_res.')
graphics.off()

# Use clustering resolution = 0.5
integrated_data <- FindClusters(integrated_data, resolution = 0.5, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(curr_plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(integrated_data)
graphics.off()

# Plot QC for each cluster
png(paste0(curr_plot_path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(integrated_data)
graphics.off()

# plot dimplot for main W gene

#_####################################################################################################
#     Heatmap clearly shows clusters segregate by sex - check this and regress out the sex effect   #
#####################################################################################################

# Change plot path
curr_plot_path <- paste0(plot_path, '1_sex_filt_integrated/')
dir.create(curr_plot_path)

# There is a strong sex effect - this plot shows DE genes between clusters 1 and 2 which are predominantly hh4 clusters. Clustering is driven by sex genes
png(paste0(curr_plot_path, 'HM.top15.DE.pre-sexfilt.png'), height = 40, width = 70, units = 'cm', res = 500)
tenx.pheatmap(data = integrated_data[,rownames(integrated_data@meta.data[integrated_data$seurat_clusters == 1 | integrated_data$seurat_clusters == 2,])],
              metadata = c("seurat_clusters", "orig.ident"), selected_genes = rownames(FindMarkers(integrated_data, ident.1 = 1, ident.2 = 2)),
              hclust_rows = T, gaps_col = "seurat_clusters")
graphics.off()

# Use W chromosome genes to K-means cluster the cells into male (zz) and female (zw)
W_genes <- as.matrix(integrated_data@assays$RNA[grepl("W-", rownames(integrated_data@assays$RNA)),])
k_clusters <- kmeans(t(W_genes), 2)
k_clusters <- data.frame(k_clusters$cluster)
integrated_data@meta.data$k_clusters <- k_clusters[match(colnames(integrated_data@assays$RNA), rownames(k_clusters)),]

# Get rownames for kmeans clusters 1 and 2
k_clus_1 <- rownames(integrated_data@meta.data[integrated_data@meta.data$k_clusters == 1,])
k_clus_2 <- rownames(integrated_data@meta.data[integrated_data@meta.data$k_clusters == 2,])

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

cell.sex.ID <- list("male.cells" = k.male, "female.cells" = k.female)
saveRDS(cell.sex.ID, paste0(rds.path, "sex_kmeans_integrated.RDS"))

# Add sex data to meta.data
integrated_data@meta.data$sex <- unlist(lapply(rownames(integrated_data@meta.data), function(x)
  if(x %in% k.male){"male"} else if(x %in% k.female){"female"} else{stop("cell sex is not assigned")}))


###### Next - subset autosomal genes and Z genes - then calculate the average for each gene for both kmeans clustered cells and plot in order to compare whether the two groups show significant differences in their expression.

# This is a test on autosomal genes to try and calculate and compare FC betweeen clusters
# Calculating median is tricky as there are a lot of dropouts in 10x data so you end up with either 0s (when the median  = 0) or 1 (when the median expression in both clusters is the same - probably a result of normalisation resulting in a UMI of 0 or 1 being normalised to a nominal value)

# Make dataframe for mean Z expression in male cells
mean.Z.male <- data.frame(Z.mean = apply(integrated_data@assays$RNA[grepl("Z-", rownames(integrated_data@assays$RNA)), k.male], 1, mean))
# add 1 before log2 as log2(1) = 0
mean.Z.male <- log2(mean.Z.male + 1)

# Make dataframe for mean Z expression in female cells
mean.Z.female <- data.frame(Z.mean = apply(integrated_data@assays$RNA[grepl("Z-", rownames(integrated_data@assays$RNA)), k.female], 1, mean))
mean.Z.female <- log2(mean.Z.female + 1)

# Make dataframe for mean autosomal expression in male cells
mean.auto.male <- data.frame(auto.mean = apply(integrated_data@assays$RNA[!grepl("Z-", rownames(integrated_data@assays$RNA)) & !grepl("W-", rownames(integrated_data@assays$RNA)), k.male], 1, mean))
mean.auto.male <- log2(mean.auto.male + 1)

# Make dataframe for mean autosomal expression in male cells
mean.auto.female <- data.frame(auto.mean = apply(integrated_data@assays$RNA[!grepl("Z-", rownames(integrated_data@assays$RNA)) & !grepl("W-", rownames(integrated_data@assays$RNA)), k.female], 1, mean))
mean.auto.female <- log2(mean.auto.female + 1)

# Calculate FC by subtracting log2 expression from each other
FC <- list()
FC$Z <- mean.Z.male - mean.Z.female
FC$auto <-  mean.auto.male - mean.auto.female

# Plot boxplot of Z gene and autosomal expression in male vs female cells
png(paste0(curr_plot_path,"sex_kmeans_log2FC_boxplot.png"), height = 18, width = 18, units = "cm", res = 200)
boxplot(c(FC$Z, FC$auto),  ylab = "male - female log2 FC (mean normalised UMI +1)", names = c("Z chromosome genes", "autosomal genes"))
graphics.off()

# Z genes are upregulated within male genes relative to female genes whereas autosomal genes have a normal distribution of logFCs

#####################################################################################################
#                                     Regress sex effect                                            #
#####################################################################################################

# Init sexscale object 
sexscale_data <- integrated_data

# Re-run findvariablefeatures and scaling
# sexscale_data <- FindVariableFeatures(sexscale_data, selection.method = "vst", nfeatures = 2000)
# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 4000 * 1024^2)

sexscale_data <- SCTransform(sexscale_data, verbose = TRUE, vars.to.regress = c("percent.mt", "sex"))

# sexscale_data <- ScaleData(sexscale_data, features = rownames(sexscale_data), vars.to.regress = c("percent.mt", "sex"))

# Save RDS
saveRDS(sexscale_data, paste0(rds.path, "sexscale_data.RDS"))

# Read in RDS data if needed
# sexscale_data <- readRDS(paste0(rds.path, "sexscale_data.RDS"))

# Set plot path
curr_plot_path <- paste0(plot_path, '1_sex_filt_integrated/')

# PCA
sexscale_data <- RunPCA(object = sexscale_data, verbose = FALSE)

png(paste0(curr_plot_path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(sexscale_data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(curr_plot_path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(sexscale_data, ndims = 40))
graphics.off()

png(paste0(curr_plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(sexscale_data, PCA.levels = c(10, 20, 30, 40), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
sexscale_data <- FindNeighbors(sexscale_data, dims = 1:30, verbose = FALSE)
sexscale_data <- RunUMAP(sexscale_data, dims = 1:30, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr_plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = sexscale_data, by = 0.1, prefix = "integrated_snn_res.")
graphics.off()

# Use clustering resolution = 0.5 to look for contamination clusters
sexscale_data <- FindClusters(sexscale_data, resolution = 0.5, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(curr_plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(sexscale_data)
graphics.off()

# Plot QC for each cluster
png(paste0(curr_plot_path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(sexscale_data)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(sexscale_data, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = sexscale_data, col.to.sort = seurat_clusters, sort.by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr_plot_path, 'HM.top15.DE.post-sexfilt.png'), height = 75, width = 100, units = 'cm', res = 500)
tenx.pheatmap(data = sexscale_data, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = 'integrated')
graphics.off()

