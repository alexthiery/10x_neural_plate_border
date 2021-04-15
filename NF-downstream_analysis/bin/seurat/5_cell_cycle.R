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

    plot_path = "./output/NF-downstream_analysis/5_cell_cycle/plots/"
    rds_path = "./output/NF-downstream_analysis/5_cell_cycle/rds_files/"
    data_path = "./output/NF-downstream_analysis/4_sex_filt/rds_files/"
    
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

png(paste0(plot_path, "dimHM.png"), width=30, height=65, units = 'cm', res = 200)
DimHeatmap(cell_cycle_data, dims = 1:40, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(cell_cycle_data, return = 'plot')
graphics.off()

pc_cutoff <- ElbowCutoff(cell_cycle_data)

png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCALevelComparison(cell_cycle_data, PCA_levels = c(10, 15, 20, 25), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
cell_cycle_data <- FindNeighbors(cell_cycle_data, dims = 1:pc_cutoff, verbose = FALSE)
cell_cycle_data <- RunUMAP(cell_cycle_data, dims = 1:pc_cutoff, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = cell_cycle_data, by = 0.2, prefix = "integrated_snn_res.")
graphics.off()

# Use clustering resolution = 1.2
cell_cycle_data <- FindClusters(cell_cycle_data, resolution = 1.2)

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(cell_cycle_data)
graphics.off()

# Plot QC for each cluster
png(paste0(plot_path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QCPlot(cell_cycle_data)
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
cluster_order <- OrderCellClusters(seurat_object = cell_cycle_data, col_to_sort = seurat_clusters, sort_by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))

png(paste0(plot_path, 'HM.top15.DE.png'), height = 75, width = 100, units = 'cm', res = 500)
TenxPheatmap(data = cell_cycle_data, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()

# Save RDS after regressing cell cycle
saveRDS(cell_cycle_data, paste0(rds_path, "cell_cycle_data.RDS"))