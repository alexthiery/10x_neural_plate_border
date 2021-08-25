#!/usr/bin/env Rscript

# Load packages
library(optparse)
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

# Read in command line opts
option_list <- list(
    make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
    make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
    make_option(c("", "--clustres"), action = "store", type = "double", help = "Clustering resolution. Default is 0.5", default = 0.5),
    make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE)
    )

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

ncores = opt$cores

# Multi-core when running from command line
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 16* 1024^3) # 32gb

# Set paths and load data
plot_path = "./plots/"
rds_path = "./rds_files/"
data_path = "./input/"

dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

seurat_data <- readRDS(list.files(data_path, full.names = TRUE))

# Set RNA to default assay
DefaultAssay(seurat_data) <- "RNA"

# Re-run findvariablefeatures and scaling
seurat_data <- FindVariableFeatures(seurat_data, selection.method = "vst", nfeatures = 2000, assay = 'RNA')

seurat_data <- ScaleData(seurat_data, features = rownames(seurat_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# Set Integrated to default assay
DefaultAssay(seurat_data) <- "integrated"

# Rescale data on integrated assay
seurat_data <- ScaleData(seurat_data, features = rownames(seurat_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# PCA
seurat_data <- RunPCA(object = seurat_data, verbose = FALSE)

png(paste0(plot_path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(seurat_data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(seurat_data, return = 'plot')
graphics.off()

# automatically determine elbow
pc_cutoff <- ElbowCutoff(seurat_data)

png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCALevelComparison(seurat_data, PCA_levels = c(pc_cutoff-5, pc_cutoff, pc_cutoff+5, pc_cutoff+10), cluster_res = opt$clustres)
graphics.off()

seurat_data <- FindNeighbors(seurat_data, dims = 1:pc_cutoff, verbose = FALSE)
seurat_data <- RunUMAP(seurat_data, dims = 1:pc_cutoff, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = seurat_data, by = 0.2, prefix = "integrated_snn_res.")
graphics.off()

# Cluster data
seurat_data <- FindClusters(seurat_data, resolution = opt$clustres)

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(seurat_data, stage_col = "stage")
graphics.off()



# Plot UMAP for developmental stage, clusters and integration (if subset contains more than one batch)
plots <- list()
if(length(unique(seurat_data$stage)) > 1){
    plots$stage_plot <- DimPlot(seurat_data, group.by = "stage") + ggtitle(paste("Developmental stage")) + theme(plot.title = element_text(hjust = 0.5))
}

plots$cluster_plot <- DimPlot(seurat_data, group.by = "seurat_clusters") + ggtitle(paste("Clusters")) + theme(plot.title = element_text(hjust = 0.5))

if(length(unique(seurat_data$run)) > 1){
    plots$integration_plot <- DimPlot(seurat_data, group.by = "run") + ggtitle(paste("Batches")) + theme(plot.title = element_text(hjust = 0.5))
}

png(paste0(plot_path, "UMAP.png"), width=20*length(plots), height=20, units = 'cm', res = 200)
do.call("grid.arrange", c(plots, nrow=1))
graphics.off()


# Plot QC for each cluster
png(paste0(plot_path, "QCPlot.png"), width=40, height=28, units = 'cm', res = 200)
QCPlot(seurat_data)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(seurat_data, only.pos = T, logfc.threshold = 0.25, assay = "RNA")
# get automated cluster order based on percentage of cells in adjacent stages
cluster_order = OrderCellClusters(seurat_object = seurat_data, col_to_sort = seurat_clusters, sort_by = stage)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))

png(paste0(plot_path, 'HM.top15.DE.seurat_data.png'), height = 75, width = 100, units = 'cm', res = 500)
TenxPheatmap(data = seurat_data, metadata = c("seurat_clusters", "stage"), custom_order_column = "seurat_clusters",
             custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = 'RNA')
graphics.off()

saveRDS(seurat_data, paste0(rds_path, "seurat_data.RDS"), compress = FALSE)

