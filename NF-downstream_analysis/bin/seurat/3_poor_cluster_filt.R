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

    plot_path = "./output/NF-downstream_analysis/3_poor_cluster_filt/plots/"
    rds_path = "./output/NF-downstream_analysis/3_poor_cluster_filt/rds_files/"
    data_path = "./output/NF-downstream_analysis/2_integration_qc/rds_files/"
    
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

poor_cluster_filt_data <- readRDS(paste0(data_path, 'integration_qc_data.RDS'))

############################## Identify and filter poor quality clusters #######################################

# Set RNA to default assay
DefaultAssay(poor_cluster_filt_data) <- "RNA"


# Plot QC for each cluster
png(paste0(plot_path, "QCPlot.png"), width=32, height=28, units = 'cm', res = 200)
QCPlot(poor_cluster_filt_data, plot_quantiles = TRUE)
graphics.off()

# Automatically find poor quality clusters
poor_clusters <- IdentifyOutliers(poor_cluster_filt_data@meta.data, group_by = 'seurat_clusters',
                                     metrics = c('nCount_RNA', 'nFeature_RNA'), intersect_metrics = TRUE)

# Plot UMAP for poor quality clusters
png(paste0(plot_path, "PoorClusters.png"), width=60, height=20, units = 'cm', res = 200)
ClusterDimplot(poor_cluster_filt_data, clusters = poor_clusters, plot_title = 'poor quality clusters')
graphics.off()


# Filter poor quality clusters
poor_quality_cells <- rownames(filter(poor_cluster_filt_data@meta.data, seurat_clusters %in% poor_clusters))
poor_cluster_filt_data <- subset(poor_cluster_filt_data, cells = poor_quality_cells, invert = T)

# Re-run findvariablefeatures and scaling
poor_cluster_filt_data <- FindVariableFeatures(poor_cluster_filt_data, selection.method = "vst", nfeatures = 2000, assay = 'RNA')
poor_cluster_filt_data <- ScaleData(poor_cluster_filt_data, features = rownames(poor_cluster_filt_data), vars.to.regress = "percent.mt")


# Set Integrated to default assay
DefaultAssay(poor_cluster_filt_data) <- "integrated"

# Rescale data on integrated assay
poor_cluster_filt_data <- ScaleData(poor_cluster_filt_data, features = rownames(poor_cluster_filt_data), vars.to.regress = "percent.mt")

# PCA
poor_cluster_filt_data <- RunPCA(object = poor_cluster_filt_data, verbose = FALSE)

png(paste0(plot_path, "dimHM.png"), width=30, height=65, units = 'cm', res = 200)
DimHeatmap(poor_cluster_filt_data, dims = 1:40, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(poor_cluster_filt_data, return = 'plot')
graphics.off()

pc_cutoff <- ElbowCutoff(poor_cluster_filt_data)

png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCALevelComparison(poor_cluster_filt_data, PCA_levels = c(10, 15, 20, 25), cluster_res = 0.5)
graphics.off()

poor_cluster_filt_data <- FindNeighbors(poor_cluster_filt_data, dims = 1:pc_cutoff, verbose = FALSE)
poor_cluster_filt_data <- RunUMAP(poor_cluster_filt_data, dims = 1:pc_cutoff, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = poor_cluster_filt_data, by = 0.2, prefix = "integrated_snn_res.")
graphics.off()

# Use default clustering resolution (0.5)
poor_cluster_filt_data <- FindClusters(poor_cluster_filt_data, resolution = 0.5)

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(poor_cluster_filt_data, stage_col = "stage")
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(poor_cluster_filt_data, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster_order <- OrderCellClusters(seurat_object = poor_cluster_filt_data, col_to_sort = seurat_clusters, sort_by = stage)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))

png(paste0(plot_path, 'HM.top15.DE.poor_cluster_filt_data.png'), height = 75, width = 100, units = 'cm', res = 500)
TenxPheatmap(data = poor_cluster_filt_data, metadata = c("seurat_clusters", "stage"), custom_order_column = "seurat_clusters",
              custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = 'integrated')
graphics.off()

saveRDS(poor_cluster_filt_data, paste0(rds_path, "poor_cluster_filt_data.RDS"))
