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
    plot_path = "./output/NF-downstream_analysis/seurat_integrate_new/plots/"
    rds_path = "./output/NF-downstream_analysis/seurat_integrate_new/rds_files/"
    data_path = "./output/NF-scRNAseq_alignment/cellranger/count/filtered_feature_bc_matrix"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/NF-scRNAseq_alignment/cellranger/count/filtered_feature_bc_matrix"
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
  library(tidyverse)
}

# Make dataframe with stage and replicate info extracted from path
input <- list.dirs(data_path, recursive = FALSE, full.names = TRUE)
input <- data.frame(sample = sub('.*/', '', input), run = str_split(sub('.*/', '', input), pattern = "_", simplify = T)[,2], path = input)

# Init list of seurat objects then merge
seurat_list <- apply(input, 1, function(x) CreateSeuratObject(counts= Read10X(data.dir = x[["path"]]), project = x[["sample"]]))
names(seurat_list) <- input$sample
seurat_all <- merge(x = seurat_list[[1]], y=seurat_list[-1], add.cell.ids = names(seurat_list), project = "chick.10x")

# Add metadata col for seq run
seurat_all@meta.data[["run"]] <- gsub(".*_", "", as.character(seurat_all@meta.data$orig.ident))
seurat_all@meta.data[["stage"]] <- gsub("_.*", "", as.character(seurat_all@meta.data$orig.ident))

# Convert metadata character cols to factors
seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)] <- lapply(seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)], as.factor)

# Remove genes expressed in fewer than 5 cells
seurat_all <- DietSeurat(seurat_all, features = names(which(Matrix::rowSums(GetAssayData(seurat_all) > 0) >=5)))

# Store mitochondrial percentage in object meta data
seurat_all <- PercentageFeatureSet(seurat_all, pattern = "^MT-", col.name = "percent.mt")

# Remove data which do not pass filter threshold
seurat_all <- subset(seurat_all, subset = c(nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15))


#####################################################################################################
#                           Integrate data from different 10x runs                                  #
#####################################################################################################


# Split object by run and find integration points
seurat_split <- SplitObject(seurat_all, split.by = "run")

# Multi-core when running from command line
if(opt$runtype == "nextflow"){
  plan("multiprocess", workers = ncores)
  options(future.globals.maxSize = 32* 1024^3) # 32gb
}

# SCTransform replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
seurat_split_SCTransform <- lapply(seurat_split, function(x) SCTransform(x, method = "glmGamPoi", verbose = TRUE, vars.to.regress = "percent.mt"))
# seurat_split_SCTransform <- lapply(seurat_split, function(x) SCTransform(x, verbose = TRUE))


# Save RDS after SCTransform as this step takes time
saveRDS(seurat_split_SCTransform, paste0(rds_path, 'seurat_split_SCTransform.RDS'))


# Find and filter integration anchors
features.sct <- SelectIntegrationFeatures(seurat_split_SCTransform, nfeatures = 2000)
seurat_split_SCTransform <- PrepSCTIntegration(seurat_split_SCTransform, anchor.features = features.sct)
seurat_split_SCTransform <- lapply(seurat_split_SCTransform, FUN = RunPCA, features = features.sct)



seurat_integrated_SCTransform <- FindIntegrationAnchors(object.list = seurat_split_SCTransform, normalization.method = "SCT",
                                  anchor.features = features.sct, dims = 1:30, reduction = "rpca", k.anchor = 5)

seurat_integrated_SCTransform <- IntegrateData(anchorset = seurat_integrated_SCTransform, normalization.method = "SCT", dims = 1:30)


# Change plot path
curr_plot_path <- paste0(plot_path)

# Run PCA analysis
integrated_data <- RunPCA(object = seurat_integrated_SCTransform, verbose = FALSE)

# Plot heatmap of top variable genes across top principle components
png(paste0(curr_plot_path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(integrated_data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

# another heuristic method is ElbowPlot which ranks PCs based on the % variance explained by each PC
png(paste0(curr_plot_path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(integrated_data, ndims = 40))
graphics.off()

# Run clustering and UMAP at different PCA cutoffs - save this output to compare the optimal number of PCs to be used
png(paste0(curr_plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(integrated_data, PCA.levels = c(5, 10, 20, 40), cluster_res = 0.5)
graphics.off()

# Use PCA=20 as elbow plot is relatively stable across stages
# Use clustering resolution = 0.5 for filtering
integrated_data <- FindNeighbors(integrated_data, dims = 1:20, verbose = FALSE)
integrated_data <- RunUMAP(integrated_data, dims = 1:20, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr_plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = integrated_data, by = 0.1, prefix = 'integrated_snn_res.')
graphics.off()

# Use clustering resolution = 0.5
integrated_data <- FindClusters(integrated_data, resolution = 0.5, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(curr_plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(integrated_data)
graphics.off()

# Plot QC for each cluster
png(paste0(curr_plot_path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(integrated_data)
graphics.off()


# Save RDS after integration
saveRDS(integrated_data, paste0(rds_path, "integrated_data.RDS"))


readRDS(paste0(rds_path, "integrated_data.RDS"))





