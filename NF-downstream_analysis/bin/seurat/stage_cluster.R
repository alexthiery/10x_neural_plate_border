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
#   'clustres', 'r', 2, "numeric",
#   'npcs', 'p', 2, "integer",
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 32gb
    
    # if(is.null(opt.clustres) && is.null(opt.npcs)){
    #   stop("--clustres and --npcs must be specified")
    # }
    
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

stage_data <- readRDS(list.files(data_path, full.names = TRUE))

#####################################################################################
#                         Split stage, scale and re-cluster                         #
#####################################################################################
# Set RNA to default assay
DefaultAssay(stage_data) <- "RNA"

# Re-run findvariablefeatures and scaling
stage_data <- FindVariableFeatures(stage_data, selection.method = "vst", nfeatures = 2000, assay = 'RNA')

stage_data <- ScaleData(stage_data, features = rownames(stage_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# Set Integrated to default assay
DefaultAssay(stage_data) <- "integrated"

# Rescale data on integrated assay
stage_data <- ScaleData(stage_data, features = rownames(stage_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# PCA
stage_data <- RunPCA(object = stage_data, verbose = FALSE)

png(paste0(plot_path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(stage_data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(stage_data, return = 'plot')
graphics.off()

pc_cutoff <- ElbowCutoff(stage_data)

png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCALevelComparison(stage_data, PCA_levels = c(10, 20, 30, 40), cluster_res = 0.5)
graphics.off()

stage_data <- FindNeighbors(stage_data, dims = 1:pc_cutoff, verbose = FALSE)
stage_data <- RunUMAP(stage_data, dims = 1:pc_cutoff, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = stage_data, by = 0.2, prefix = "integrated_snn_res.")
graphics.off()

# Use clustering resolution = 1 in order to identify any remaining poor quality cells
stage_data <- FindClusters(stage_data, resolution = 1)

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(stage_data, stage_col = "stage")
graphics.off()

# Plot QC for each cluster
png(paste0(plot_path, "QCPlot.png"), width=40, height=28, units = 'cm', res = 200)
QCPlot(stage_data)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(stage_data, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster_order = OrderCellClusters(seurat_object = stage_data, col_to_sort = seurat_clusters, sort_by = stage)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))

png(paste0(plot_path, 'HM.top15.DE.stage_data.png'), height = 75, width = 100, units = 'cm', res = 500)
TenxPheatmap(data = stage_data, metadata = c("seurat_clusters", "stage"), custom_order_column = "seurat_clusters",
             custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = 'integrated')
graphics.off()

saveRDS(stage_data, paste0(rds_path, "stage_data.RDS"), compress = FALSE)