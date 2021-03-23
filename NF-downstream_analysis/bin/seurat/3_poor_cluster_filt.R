#!/usr/bin/env Rscript

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
    plot_path = "./output/NF-downstream_analysis/3_poor_cluster_filt/plots/"
    rds_path = "./output/NF-downstream_analysis/3_poor_cluster_filt/rds_files/"

    data_path = "./output/NF-downstream_analysis/2_integration_qc/rds_files/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores

    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 32* 1024^3) # 32gb
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


poor_cluster_filt_data <- readRDS(paste0(data_path, 'integration_qc_data.RDS'))

############################### Remove poor quality clusters ########################################

# Set RNA to default assay
DefaultAssay(poor_cluster_filt_data) <- "RNA"

poor_quality_cells <- rownames(filter(poor_cluster_filt_data@meta.data, seurat_clusters %in% c(11, 15)))

poor_cluster_filt_data <- subset(poor_cluster_filt_data, cells = poor_quality_cells, invert = T)

# Re-run findvariablefeatures and scaling
poor_cluster_filt_data <- FindVariableFeatures(poor_cluster_filt_data, selection.method = "vst", nfeatures = 2000, assay = 'RNA')

poor_cluster_filt_data <- ScaleData(poor_cluster_filt_data, features = rownames(poor_cluster_filt_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

saveRDS(poor_cluster_filt_data, paste0(rds.path, "poor_cluster_filt_data.RDS"))


# Set Integrated to default assay
DefaultAssay(poor_cluster_filt_data) <- "integrated"

# Rescale data on integrated assay
poor_cluster_filt_data <- ScaleData(poor_cluster_filt_data, features = rownames(poor_cluster_filt_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# PCA
poor_cluster_filt_data <- RunPCA(object = poor_cluster_filt_data, verbose = FALSE)

png(paste0(plot_path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(poor_cluster_filt_data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(poor_cluster_filt_data, ndims = 40))
graphics.off()

png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(poor_cluster_filt_data, PCA.levels = c(10, 20, 30, 40), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
poor_cluster_filt_data <- FindNeighbors(poor_cluster_filt_data, dims = 1:30, verbose = FALSE)
poor_cluster_filt_data <- RunUMAP(poor_cluster_filt_data, dims = 1:30, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = poor_cluster_filt_data, by = 0.2, prefix = "integrated_snn_res.")
graphics.off()

# Use clustering resolution = 1.2
poor_cluster_filt_data <- FindClusters(poor_cluster_filt_data, resolution = 1.2)

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(poor_cluster_filt_data, stage.col = "stage")
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(poor_cluster_filt_data, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster_order = order.cell.stage.clust(seurat_object = poor_cluster_filt_data, col.to.sort = seurat_clusters, sort.by = stage)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))

png(paste0(plot_path, 'HM.top15.DE.poor_cluster_filt_data.png'), height = 75, width = 100, units = 'cm', res = 500)
tenx.pheatmap(data = poor_cluster_filt_data, metadata = c("seurat_clusters", "stage"), custom_order_column = "seurat_clusters",
              custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = 'integrated')
graphics.off()

saveRDS(poor_cluster_filt_data, paste0(rds_path, "poor_cluster_filt_data.RDS"))
