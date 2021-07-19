#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
library(Seurat)
library(scHelper)
library(tidyverse)
library(scatterplot3d)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    plot_path = "./output/NF-downstream_analysis/6_contamination_filt/plots/"
    data_path = "./output/NF-downstream_analysis/5_cell_cycle/rds_files/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/rds_files/"
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }

  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

seurat_data <- readRDS(list.files(data_path, full.names = TRUE))

seurat_data <- RunPCA(object = seurat_data, verbose = FALSE)
# automatically determine elbow
pc_cutoff <- ElbowCutoff(seurat_data)
seurat_data <- FindNeighbors(seurat_data, dims = 1:pc_cutoff, verbose = FALSE)
seurat_data <- RunUMAP(seurat_data, dims = 1:pc_cutoff, verbose = FALSE, n.components = 3L)

# extract dims for plotting
plot_data <- FetchData(object = seurat_data, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters")) %>%
  rownames_to_column(var = "cell_id")

plot_data <- plot_data %>% arrange(seurat_clusters)

plot_colours <- scHelper::ggPlotColours(n = length(unique(plot_data$seurat_clusters)))[
  as.numeric(plot_data$seurat_clusters)]

# plot scatterplots in 3 diff orientations
par(mar = c(0.1, 0.1, 0.1, 0.1))

png(filename = paste0(plot_path, "scatterplot3d_1.png"), width = 8, height = 6, units = "in", res = 300)
scatterplot3d(x = plot_data$UMAP_1, y = plot_data$UMAP_2, z = plot_data$UMAP_3,
              color=plot_colours, cex.symbols = 0.3, pch = 19, angle = -30)
graphics.off()

png(filename = paste0(plot_path, "scatterplot3d_2.png"), width = 8, height = 6, units = "in", res = 300)
scatterplot3d(x = plot_data$UMAP_2, y = plot_data$UMAP_3, z = plot_data$UMAP_1,
              color=plot_colours, cex.symbols = 0.3, pch = 19, angle = -30)
graphics.off()

png(filename = paste0(plot_path, "scatterplot3d_2.png"), width = 8, height = 6, units = "in", res = 300)
scatterplot3d(x = plot_data$UMAP_3, y = plot_data$UMAP_1, z = plot_data$UMAP_2,
              color=plot_colours, cex.symbols = 0.3, pch = 19, angle = -30)
graphics.off()
