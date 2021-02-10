#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
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
    plot_path = "./output/NF-downstream_analysis/2_seurat/plots/"
    rds_path = "./output/NF-downstream_analysis/2_seurat/rds_files/"
    data_path = "./output/NF-downstream_analysis/scRNAseq/seurat_integrate/rds_files/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  
  # Load packages - packages are stored within renv in the repository
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
  library(STACAS)
  library(tidyverse)
}


seurat_integrated_SCTransform <- readRDS(paste0(data_path, "seurat_integrated_SCTransform.RDS"))

# set integrated count data as default
# DefaultAssay(seurat_integrated_SCTransform) <- "integrated"


seurat_integrated_SCTransform <- RunPCA(object = seurat_integrated_SCTransform, verbose = FALSE)
seurat_integrated_SCTransform <- FindNeighbors(seurat_integrated_SCTransform, dims = 1:30, verbose = FALSE)
seurat_integrated_SCTransform <- RunUMAP(seurat_integrated_SCTransform, dims = 1:30, verbose = FALSE)
seurat_integrated_SCTransform <- FindClusters(seurat_integrated_SCTransform, resolution = 0.5, verbose = FALSE)

png(paste0(plot_path, "UMAP_sctransform.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(seurat_integrated_SCTransform)
graphics.off()



seurat_integrated_scale <- readRDS(paste0(data_path, "seurat_integrated_scale.RDS"))

# set integrated count data as default
DefaultAssay(seurat_integrated_scale) <- "integrated"


seurat_integrated_scale <- RunPCA(object = seurat_integrated_scale, verbose = FALSE)
seurat_integrated_scale <- FindNeighbors(seurat_integrated_scale, dims = 1:30, verbose = FALSE)
seurat_integrated_scale <- RunUMAP(seurat_integrated_scale, dims = 1:30, verbose = FALSE)
seurat_integrated_scale <- FindClusters(seurat_integrated_scale, resolution = 0.5, verbose = FALSE)

png(paste0(plot_path, "UMAP_scale.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(seurat_integrated_scale)
graphics.off()

