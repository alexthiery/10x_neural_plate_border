#!/usr/bin/env Rscript

# Load packages
library(getopt)
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

    plot_path = "./output/NF-downstream_analysis/seurat/3_integration_qc/plots/"
    rds_path = "./output/NF-downstream_analysis/seurat/3_integration_qc/rds_files/"
    data_path = "./output/NF-downstream_analysis/seurat/2_integration/rds_files/"
    
    ncores = 8

  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')

    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores

    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 32gb

  } else {
    stop("--runtype must be set to 'nextflow'")
  }

  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

integration_qc_data <- readRDS(list.files(data_path, full.names = TRUE))


# ################################################################################################
# # Plot batch effect before integration

# # Set RNA assay to check batch effect before integration
# DefaultAssay(integration_qc_data) <- "RNA"

# # Log normalize data and find variable features
# integration_qc_data <- NormalizeData(integration_qc_data, normalization.method = "LogNormalize", scale.factor = 10000)
# integration_qc_data <- FindVariableFeatures(integration_qc_data, selection.method = "vst", nfeatures = 2000)

# # Run PCA analysis
# integration_qc_data <- RunPCA(object = integration_qc_data, verbose = FALSE)
# pc_cutoff <- ElbowCutoff(integration_qc_data)

# integration_qc_data <- RunUMAP(integration_qc_data, dims = 1:pc_cutoff, verbose = FALSE)
# integration_qc_data <- FindNeighbors(integration_qc_data, dims = 1:pc_cutoff, verbose = FALSE)

# # check whether stages that are resequenced are well integrated
# png(paste0(plot_path, "CheckBatchEffect.png"), width=60, height=20, units = 'cm', res = 200)
# CheckIntegration(integration_qc_data)
# graphics.off()

# ################################################################################################


# Set integrated assay as default for clustering
DefaultAssay(integration_qc_data) <- "integrated"

# Run PCA analysis
integration_qc_data <- RunPCA(object = integration_qc_data, verbose = FALSE)

# Plot heatmap of top variable genes across top principle components
png(paste0(plot_path, "dimHM.png"), width=30, height=65, units = 'cm', res = 200)
DimHeatmap(integration_qc_data, dims = 1:40, balanced = TRUE, cells = 500)
graphics.off()

# In order to calculate the PC cutoff for downstream analysis, we unbiasedly calculate the point at which PC elbow starts.
# First we take the larger value of the point where the principal components only contribute 5% of standard deviation and the point where the principal components cumulatively contribute 90% of the standard deviation.
# Next we take the point where the percent change in variation between the consecutive PCs is less than 0.1%.
# The smaller out of these two values is determined at the elbow cutoff

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(integration_qc_data, return = 'plot')
graphics.off()

pc_cutoff <- ElbowCutoff(integration_qc_data)

# Run clustering and UMAP at different PCA cutoffs - save this output to compare the optimal number of PCs to be used
png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCALevelComparison(integration_qc_data, PCA_levels = c(pc_cutoff-5, pc_cutoff, pc_cutoff+5, pc_cutoff+10), cluster_res = 1)
graphics.off()

# Use clustering resolution = 0.5 for filtering
integration_qc_data <- RunUMAP(integration_qc_data, dims = 1:pc_cutoff, verbose = FALSE)
integration_qc_data <- FindNeighbors(integration_qc_data, dims = 1:pc_cutoff, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = integration_qc_data, by = 0.2, prefix = 'integrated_snn_res.')
graphics.off()

############################## Identify poor quality clusters #######################################

# Use higher cluster resolution for filtering poor clusters
integration_qc_data <- FindClusters(integration_qc_data, resolution = 1, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(integration_qc_data)
graphics.off()

# Plot QC for each cluster
png(paste0(plot_path, "QCPlot.png"), width=32, height=28, units = 'cm', res = 200)
QCPlot(integration_qc_data, quantiles = c(0.25, 0.75))
graphics.off()

# Automatically find poor quality clusters
poor_clusters <- IdentifyOutliers(integration_qc_data, group_by = 'seurat_clusters',
                                     metrics = c('nCount_RNA', 'nFeature_RNA'), quantiles = c(0.25, 0.75), intersect_metrics = TRUE)

# Plot UMAP for poor quality clusters
png(paste0(plot_path, "PoorClusters.png"), width=60, height=20, units = 'cm', res = 200)
ClusterDimplot(integration_qc_data, clusters = poor_clusters, plot_title = 'poor quality clusters')
graphics.off()

# check whether stages that are resequenced are well integrated
png(paste0(plot_path, "CheckIntegration.png"), width=60, height=20, units = 'cm', res = 200)
CheckIntegration(integration_qc_data)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(integration_qc_data, only.pos = T, logfc.threshold = 0.25, assay = "RNA")
# get automated cluster order based on percentage of cells in adjacent stages
cluster_order <- OrderCellClusters(seurat_object = integration_qc_data, col_to_sort = 'seurat_clusters', sort_by = 'stage')
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))

png(paste0(plot_path, 'HM.top15.DE.integration_qc_data.png'), height = 75, width = 100, units = 'cm', res = 500)
TenxPheatmap(data = integration_qc_data, metadata = c("seurat_clusters", "stage"), custom_order_column = "seurat_clusters",
              custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = 'RNA')
graphics.off()

# Save RDS
saveRDS(integration_qc_data, paste0(rds_path, "integration_qc_data.RDS"), compress = FALSE)