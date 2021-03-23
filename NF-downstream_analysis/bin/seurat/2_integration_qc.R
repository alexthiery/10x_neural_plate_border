#!/usr/bin/env Rscript

# Load packages
library(getopt)
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

# Define arguments for Rscript
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
    plot_path = "./output/NF-downstream_analysis/integration_qc/plots/"
    rds_path = "./output/NF-downstream_analysis/integration_qc/rds_files/"

    # switch depending on which integration to work from
    data_path = "./output/NF-downstream_analysis/integration_seurat/rds_files/"
    # data_path = "./output/NF-downstream_analysis/integration_STACAS/rds_files/"
    
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
}

integration_qc_data <- readRDS(paste0(data_path, 'intergration_data.RDS'))

# Set integrated assay as default for clustering
DefaultAssay(integration_qc_data) <- "integrated"

# Run PCA analysis
integration_qc_data <- RunPCA(object = integration_qc_data, verbose = FALSE)

# Plot heatmap of top variable genes across top principle components
png(paste0(plot_path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(integration_qc_data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

# another heuristic method is ElbowPlot which ranks PCs based on the % variance explained by each PC
png(paste0(plot_path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(integration_qc_data, ndims = 40))
graphics.off()

# Run clustering and UMAP at different PCA cutoffs - save this output to compare the optimal number of PCs to be used
png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(integration_qc_data, PCA.levels = c(5, 10, 20, 40), cluster_res = 0.5)
graphics.off()

# Use PCA=20 as elbow plot is relatively stable across stages
# Use clustering resolution = 0.5 for filtering
integration_qc_data <- RunUMAP(integration_qc_data, dims = 1:20, verbose = FALSE)
integration_qc_data <- FindNeighbors(integration_qc_data, dims = 1:20, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = integration_qc_data, by = 0.1, prefix = 'integrated_snn_res.')
graphics.off()

# Use clustering resolution = 0.5
integration_qc_data <- FindClusters(integration_qc_data, resolution = 0.5, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(integration_qc_data)
graphics.off()

# Plot QC for each cluster
png(paste0(plot_path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(integration_qc_data)
graphics.off()

# check whether stages that are resequenced are well integrated
png(paste0(plot_path, "check.integration.png"), width=60, height=20, units = 'cm', res = 200)
check.integration(integration_qc_data)
graphics.off()

# Save RDS
saveRDS(integration_qc_data, paste0(rds_path, "integration_qc_data.RDS"))