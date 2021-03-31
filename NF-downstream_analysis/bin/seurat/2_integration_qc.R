#!/usr/bin/env Rscript

# Load packages
library(getopt)
reticulate::use_python('/usr/bin/python3.7')
library(Seurat)
library(future)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(scHelper)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')

    plot_path = "./output/NF-downstream_analysis/integration_qc/plots/"
    rds_path = "./output/NF-downstream_analysis/integration_qc/rds_files/"
    data_path = "./output/NF-downstream_analysis/integration_STACAS/rds_files/"
    
    ncores = 8

  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')

    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores

    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 32* 1024^3) # 32gb

  } else {
    stop("--runtype must be set to 'nextflow'")
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
PCALevelComparison(integration_qc_data, PCA_levels = c(5, 10, 20, 40), cluster_res = 0.5)
graphics.off()

# Use PCA=20 as elbow plot is relatively stable across stages
# Use clustering resolution = 0.5 for filtering
integration_qc_data <- RunUMAP(integration_qc_data, dims = 1:20, verbose = FALSE)
integration_qc_data <- FindNeighbors(integration_qc_data, dims = 1:20, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = integration_qc_data, prefix = 'integrated_snn_res.')
graphics.off()

############################## Identify poor quality clusters #######################################

# Use higher cluster res in order to filter poor quality clusters
integration_qc_data <- FindClusters(integration_qc_data, resolution = 0.8, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(integration_qc_data)
graphics.off()

# Plot QC for each cluster
png(paste0(plot_path, "QCPlot.png"), width=40, height=14, units = 'cm', res = 200)
QCPlot(integration_qc_data)
graphics.off()

# check whether stages that are resequenced are well integrated
png(paste0(plot_path, "CheckIntegration.png"), width=60, height=20, units = 'cm', res = 200)
CheckIntegration(integration_qc_data)
graphics.off()

# Save RDS
saveRDS(integration_qc_data, paste0(rds_path, "integration_qc_data.RDS"))