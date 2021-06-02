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
  'cores', 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
data_path = "./input/"
stage_cluster_data <- readRDS(list.files(data_path, full.names = TRUE))

# Set stage var based on input
if(length(unique(stage_cluster_data$stage)) == 1){
  stage = unique(stage_cluster_data$stage)
} else {
  stop('Input RDS object contains more than 1 developmental stage')
}

plot_path = paste0("./plots/", stage, "/")
rds_path =paste0("./rds_files/")
dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

# Multi-core when running from command line
ncores = opt$cores
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 4* 1024^3)
cat(paste0("script ran with ", ncores, " cores\n"))

#####################################################################################
#                         Split stage, scale and re-cluster                         #
#####################################################################################
# Set RNA to default assay
DefaultAssay(stage_cluster_data) <- "RNA"

# Re-run findvariablefeatures and scaling
stage_cluster_data <- FindVariableFeatures(stage_cluster_data, selection.method = "vst", nfeatures = 2000, assay = 'RNA')

stage_cluster_data <- ScaleData(stage_cluster_data, features = rownames(stage_cluster_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# Set Integrated to default assay
DefaultAssay(stage_cluster_data) <- "integrated"

# Rescale data on integrated assay
stage_cluster_data <- ScaleData(stage_cluster_data, features = rownames(stage_cluster_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# PCA
stage_cluster_data <- RunPCA(object = stage_cluster_data, verbose = FALSE)

png(paste0(plot_path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(stage_cluster_data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(stage_cluster_data, return = 'plot')
graphics.off()

pc_cutoff <- ElbowCutoff(stage_cluster_data)

png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCALevelComparison(stage_cluster_data, PCA_levels = c(10, 20, 30, 40), cluster_res = 0.5)
graphics.off()

stage_cluster_data <- FindNeighbors(stage_cluster_data, dims = 1:pc_cutoff, verbose = FALSE)
stage_cluster_data <- RunUMAP(stage_cluster_data, dims = 1:pc_cutoff, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = stage_cluster_data, by = 0.1, prefix = "integrated_snn_res.")
graphics.off()

# Use clustering resolution = 0.5 (default)
stage_cluster_data <- FindClusters(stage_cluster_data, resolution = 0.5)

# Plot UMAP for clusters
png(paste0(plot_path, "UMAP.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(stage_cluster_data, group.by = "seurat_clusters") + 
  ggtitle(paste("Clusters")) + theme(plot.title = element_text(hjust = 0.5))
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(stage_cluster_data, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster_order = OrderCellClusters(seurat_object = stage_cluster_data, col_to_sort = seurat_clusters, sort_by = stage)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))

png(paste0(plot_path, 'HM.top15.DE.stage_cluster_data.png'), height = 75, width = 100, units = 'cm', res = 500)
TenxPheatmap(data = stage_cluster_data, metadata = c("seurat_clusters", "stage"), custom_order_column = "seurat_clusters",
             custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = 'integrated')
graphics.off()

saveRDS(stage_cluster_data, paste0(rds_path, stage, "_stage_cluster_data.RDS"), compress = FALSE)