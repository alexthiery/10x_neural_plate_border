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
    plot_path = "./output/NF-downstream_analysis/6_contamination_filt/plots/"
    rds_path = "./output/NF-downstream_analysis/6_contamination_filt/rds_files/"

    data_path = "./output/NF-downstream_analysis/5_cell_cycle/rds_files/"
    
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


contamination_filt_data <- readRDS(paste0(data_path, 'cell_cycle_data.RDS'))

# Set RNA to default assay
DefaultAssay(contamination_filt_data) <- "RNA"

#####################################################################################################
#                           Identify and remove contamination (mesoderm and PGCs)                   #
#####################################################################################################

# Identify mesoderm and PGCs
# UMAP plots GOI
genes <- c("EYA2", "SIX1", "TWIST1", "PITX2", "SOX17", "DAZL", "DND1", "CXCR4")

ncol = 4
png(paste0(plot_path, "UMAP_GOI.png"), width = ncol*10, height = 10*ceiling((length(genes)+1)/ncol), units = "cm", res = 200)
multi.feature.plot(seurat.obj = contamination_filt_data, gene.list = genes, plot.clusters = T,
                   plot.stage = T, label = "", cluster.col = "integrated_snn_res.0.5", n.col = ncol)
graphics.off()

# Dotplot for identifying PGCs, Early mesoderm and Late mesoderm
png(paste0(plot_path, "dotplot.GOI.png"), width = 20, height = 12, units = "cm", res = 200)
DotPlot(contamination_filt_data, features = c( "SOX17", "CXCR4","EYA2", "TWIST1", "SIX1",  "PITX2", "DAZL", "CDH5", "TAL1", "HBZ"))
graphics.off()


############################### Remove contaminating cells from clusters ########################################

# Clust 13 = PGC's - expresses dazl very highly

# Clust 11,12 = early mesoderm - expresses sox17, eya2, pitx2, cxcr4
# Clust 8,10 = Late mesoderm - expresses twist1, six1, eya2

filter_cells <- rownames(filter(contamination_filt_data@meta.data, seurat_clusters %in% c(8, 10, 11, 12, 13, 14)))

contamination_filt_data <- subset(contamination_filt_data, cells = filter_cells, invert = T)

# Re-run findvariablefeatures and scaling
contamination_filt_data <- FindVariableFeatures(contamination_filt_data, selection.method = "vst", nfeatures = 2000, assay = 'RNA')

contamination_filt_data <- ScaleData(contamination_filt_data, features = rownames(contamination_filt_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

saveRDS(contamination_filt_data, paste0(rds_path, "contamination_filt_data.RDS")) # remove once script works


# Set Integrated to default assay
DefaultAssay(contamination_filt_data) <- "integrated"

# Rescale data on integrated assay
contamination_filt_data <- ScaleData(contamination_filt_data, features = rownames(contamination_filt_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# PCA
contamination_filt_data <- RunPCA(object = contamination_filt_data, verbose = FALSE)

png(paste0(plot_path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(contamination_filt_data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(contamination_filt_data, ndims = 40))
graphics.off()

png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(contamination_filt_data, PCA.levels = c(10, 20, 30, 40), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
contamination_filt_data <- FindNeighbors(contamination_filt_data, dims = 1:30, verbose = FALSE)
contamination_filt_data <- RunUMAP(contamination_filt_data, dims = 1:30, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = contamination_filt_data, by = 0.2, prefix = "integrated_snn_res.")
graphics.off()

# Use clustering resolution = 1 in order to make lots of clusters and identify any remaining poor quality cells
contamination_filt_data <- FindClusters(contamination_filt_data, resolution = 1)

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(contamination_filt_data, stage.col = "stage")
graphics.off()

# Plot QC for each cluster
png(paste0(plot_path, "cluster.QC.png"), width=60, height=14, units = 'cm', res = 200)
QC.plot(contamination_filt_data)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(contamination_filt_data, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster_order = order.cell.stage.clust(seurat_object = contamination_filt_data, col.to.sort = seurat_clusters, sort.by = stage)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))

png(paste0(plot_path, 'HM.top15.DE.contamination_filt_data.png'), height = 75, width = 100, units = 'cm', res = 500)
tenx.pheatmap(data = contamination_filt_data, metadata = c("seurat_clusters", "stage"), custom_order_column = "seurat_clusters",
              custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = 'integrated')
graphics.off()

saveRDS(contamination_filt_data, paste0(rds_path, "contamination_filt_data.RDS"))