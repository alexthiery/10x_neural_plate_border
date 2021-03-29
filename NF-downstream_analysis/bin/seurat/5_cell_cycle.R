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
    plot_path = "./output/NF-downstream_analysis/5_cell_cycle/plots/"
    rds_path = "./output/NF-downstream_analysis/5_cell_cycle/rds_files/"

    data_path = "./output/NF-downstream_analysis/4_sex_filt/rds_files/"
    
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

pre_cell_cycle_data <- readRDS(paste0(data_path, 'sex_filt_data.RDS'))

####################################################################################
#                            Check for cell cycle effect                           #
####################################################################################

# Calculate cell cycle effect on integrated dataset
DefaultAssay(pre_cell_cycle_data) <- "integrated"

# Calculate cell cycle for each cell
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pre_cell_cycle_data <- CellCycleScoring(pre_cell_cycle_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Scale data and regress cell cycle
cell_cycle_data <- ScaleData(pre_cell_cycle_data, features = rownames(pre_cell_cycle_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# PCA
cell_cycle_data <- RunPCA(object = cell_cycle_data, verbose = FALSE)

png(paste0(plot_path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(cell_cycle_data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(cell_cycle_data, ndims = 40))
graphics.off()

png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(cell_cycle_data, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
cell_cycle_data <- FindNeighbors(cell_cycle_data, dims = 1:15, verbose = FALSE)
cell_cycle_data <- RunUMAP(cell_cycle_data, dims = 1:15, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = cell_cycle_data, by = 0.2, prefix = "integrated_snn_res.")
graphics.off()

# Use clustering resolution = 1.2
cell_cycle_data <- FindClusters(cell_cycle_data, resolution = 1.2)

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(cell_cycle_data)
graphics.off()

# Plot QC for each cluster
png(paste0(plot_path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(cell_cycle_data)
graphics.off()

# UMAP of cell cycle before and after regressing out
png(paste0(plot_path, "cell.cycle.png"), width=40, height=20, units = 'cm', res = 200)
print(gridExtra::grid.arrange(DimPlot(pre_cell_cycle_data, group.by = "Phase", reduction = "umap"),
                              DimPlot(cell_cycle_data, group.by = "Phase", reduction = "umap"),
                              ncol = 2))
graphics.off()


# switch to RNA assay for viewing expression data
DefaultAssay(cell_cycle_data) <- "RNA"

# Find variable features and scale data on RNA assay
cell_cycle_data <- FindVariableFeatures(cell_cycle_data, selection.method = "vst", nfeatures = 2000)
cell_cycle_data <- ScaleData(cell_cycle_data, features = rownames(cell_cycle_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# Find deferentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(cell_cycle_data, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster_order = order.cell.stage.clust(seurat_object = cell_cycle_data, col.to.sort = seurat_clusters, sort.by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))

png(paste0(plot_path, 'HM.top15.DE.png'), height = 75, width = 100, units = 'cm', res = 500)
tenx.pheatmap(data = cell_cycle_data, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()

# Save RDS after regressing cell cycle
saveRDS(cell_cycle_data, paste0(rds_path, "cell_cycle_data.RDS"))