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
    make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE),
    make_option(c("-m", "--meta_col"), action = "store", type = "character", help = "Column name specifying cell type column", default = 'seurat_clusters')
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

# Retrieve seurat object label
label <- sub('_.*', '', list.files(data_path))

# Load seurat data
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

# if pc_cutoff is smaller than 7 then don't run with pc_cutoff-5 as too small to run UMAP
cutoffs = ifelse(pc_cutoff < 7, c(pc_cutoff, pc_cutoff+5, pc_cutoff+10, pc_cutoff+15), c(pc_cutoff-5, pc_cutoff, pc_cutoff+5, pc_cutoff+10))
png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCALevelComparison(seurat_data, PCA_levels = cutoffs, cluster_res = opt$clustres)
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

plots$cluster_plot <- DimPlot(seurat_data, group.by = opt$meta_col) + ggtitle(paste(opt$meta_col)) + theme(plot.title = element_text(hjust = 0.5))

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
cluster_order = OrderCellClusters(seurat_object = seurat_data, col_to_sort = opt$meta_col, sort_by = 'stage')
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))

png(paste0(plot_path, 'HM.top15.DE.seurat_data.png'), height = 75, width = 100, units = 'cm', res = 500)
TenxPheatmap(data = seurat_data, metadata = c(opt$meta_col, "stage"), custom_order_column = opt$meta_col,
             custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = opt$meta_col, assay = 'RNA')
graphics.off()


# Plot feature plots for all variable genes
# Set RNA to default assay
DefaultAssay(seurat_data) <- "RNA"

dir.create(paste0(plot_path, 'feature_plots/'))
for(i in seurat_data@assays$RNA@var.features){
    png(paste0(plot_path, 'feature_plots/', i, '.png'), height = 12, width = 12, units = 'cm', res = 100)
    print(
        FeaturePlot(seurat_data, features = i, pt.size = 1.4) +
            theme_void() +
            theme(plot.title = element_blank(),
                legend.text = element_text(size=16),
                legend.key.size = unit(1, 'cm'))
        )
    graphics.off()
}

system(paste0("zip -rj ", plot_path, "feature_plots.zip ", paste0(plot_path, 'feature_plots/')))
unlink(paste0(plot_path, 'feature_plots/'), recursive=TRUE, force=TRUE)

saveRDS(seurat_data, paste0(rds_path, label, "_clustered_data.RDS"), compress = FALSE)


