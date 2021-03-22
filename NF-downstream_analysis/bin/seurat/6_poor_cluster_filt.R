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
    plot_path = "./output/NF-downstream_analysis/poor_cluster_filt/plots/"
    rds_path = "./output/NF-downstream_analysis/poor_cluster_filt/rds_files/"

    data_path = "./output/NF-downstream_analysis/contamination_filt/rds_files/"
    
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

contamination_filt_data <- readRDS(paste0(data_path, 'sexfilt_data.RDS'))


############################### Remove poor quality clusters ########################################

norm.data.clustfilt <- rownames(contamination_filt_data@meta.data)[contamination_filt_data@meta.data$seurat_clusters ==  11 |
                                                                  contamination_filt_data@meta.data$seurat_clusters == 15]

norm.data.clustfilt <- subset(contamination_filt_data, cells = norm.data.clustfilt, invert = T)

# Re-run findvariablefeatures and scaling
norm.data.clustfilt <- FindVariableFeatures(norm.data.clustfilt, selection.method = "vst", nfeatures = 2000, assay = 'RNA')

# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

norm.data.clustfilt <- ScaleData(norm.data.clustfilt, features = rownames(norm.data.clustfilt), vars.to.regress = c("percent.mt", "sex"))
saveRDS(norm.data.clustfilt, paste0(rds.path, "norm.data.clustfilt.RDS"))

# Read in RDS data if needed
# norm.data.clustfilt <- readRDS(paste0(rds.path, "norm.data.clustfilt.RDS"))

# Change plot path
curr.plot.path <- paste0(plot.path, "3_cluster_filt/")
dir.create(curr.plot.path)

# PCA
norm.data.clustfilt <- RunPCA(object = norm.data.clustfilt, verbose = FALSE)

png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(norm.data.clustfilt, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(norm.data.clustfilt, ndims = 40))
graphics.off()

png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(norm.data.clustfilt, PCA.levels = c(10, 20, 30, 40), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
norm.data.clustfilt <- FindNeighbors(norm.data.clustfilt, dims = 1:30, verbose = FALSE)
norm.data.clustfilt <- RunUMAP(norm.data.clustfilt, dims = 1:30, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr.plot.path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = norm.data.clustfilt, by = 0.2, prefix = "integrated_snn_res.")
graphics.off()

# Use clustering resolution = 1.2
norm.data.clustfilt <- FindClusters(norm.data.clustfilt, resolution = 1.2)

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(norm.data.clustfilt, stage.col = "stage")
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data.clustfilt, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = norm.data.clustfilt, col.to.sort = seurat_clusters, sort.by = stage)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr.plot.path, 'HM.top15.norm.data.clustfilt.png'), height = 75, width = 100, units = 'cm', res = 500)
tenx.pheatmap(data = norm.data.clustfilt, metadata = c("seurat_clusters", "stage"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = 'integrated')
graphics.off()

saveRDS(norm.data.clustfilt, paste0(rds.path, "seurat_out.RDS"))

# norm.data.clustfilt <- readRDS(paste0(rds.path, "seurat_out.RDS"))

############################### Cluster without HH4 ########################################

norm.data.hh4filt <- rownames(norm.data.clustfilt@meta.data)[grepl('hh4', rownames(norm.data.clustfilt@meta.data))]

norm.data.hh4filt <- subset(norm.data.clustfilt, cells = norm.data.hh4filt, invert = T)

# Re-run findvariablefeatures and scaling
norm.data.hh4filt <- FindVariableFeatures(norm.data.hh4filt, selection.method = "vst", nfeatures = 2000, assay = 'RNA')

# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

norm.data.hh4filt <- ScaleData(norm.data.hh4filt, features = rownames(norm.data.hh4filt), vars.to.regress = c("percent.mt", "sex"))
saveRDS(norm.data.hh4filt, paste0(rds.path, "norm.data.hh4filt.RDS"))

# Read in RDS data if needed
# norm.data.hh4filt <- readRDS(paste0(rds.path, "norm.data.hh4filt.RDS"))

# Change plot path
curr.plot.path <- paste0(plot.path, "4_hh4_filt/")
dir.create(curr.plot.path)

# PCA
norm.data.hh4filt <- RunPCA(object = norm.data.hh4filt, verbose = FALSE)

png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(norm.data.hh4filt, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(norm.data.hh4filt, ndims = 40))
graphics.off()

png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(norm.data.hh4filt, PCA.levels = c(10, 20, 30, 40), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
norm.data.hh4filt <- FindNeighbors(norm.data.hh4filt, dims = 1:30, verbose = FALSE)
norm.data.hh4filt <- RunUMAP(norm.data.hh4filt, dims = 1:30, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr.plot.path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = norm.data.hh4filt, by = 0.2, prefix = "integrated_snn_res.")
graphics.off()

# Use clustering resolution = 1.2
norm.data.hh4filt <- FindClusters(norm.data.hh4filt, resolution = 1.2)

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(norm.data.hh4filt, stage.col = "stage")
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data.hh4filt, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = norm.data.hh4filt, col.to.sort = seurat_clusters, sort.by = stage)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr.plot.path, 'HM.top15.norm.data.hh4filt.png'), height = 75, width = 100, units = 'cm', res = 500)
tenx.pheatmap(data = norm.data.hh4filt, metadata = c("seurat_clusters", "stage"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = 'integrated')
graphics.off()

saveRDS(norm.data.hh4filt, paste0(rds.path, "seurat_out_hh4filt.RDS"))
