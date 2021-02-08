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
    plot_path = "./output/NF-downstream_analysis/1_seurat_integrate/plots/"
    rds_path = "./output/NF-downstream_analysis/1_seurat_integrate/rds_files/"
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
  library(STACAS)
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



# #####################################################################################################
# #                           Integrate data from different 10x runs                                  #
# #####################################################################################################

# Split object by run and find integration points
seurat_split <- SplitObject(seurat_all, split.by = "run")

# Multi-core when running from command line
if(opt$runtype == "nextflow"){
  plan("multiprocess", workers = ncores)
  options(future.globals.maxSize = 32* 1024^3) # 32gb
}

# SCTransform replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
seurat_split_SCTransform <- lapply(seurat_split, function(x) SCTransform(x, verbose = TRUE, vars.to.regress = "percent.mt"))

# Save RDS after SCTransform as this step takes time
saveRDS(seurat_split_SCTransform, paste0(rds_path, 'seurat_split_SCTransform.RDS'))


# Find and filter integration anchors
features.sct <- SelectIntegrationFeatures(seurat_split_SCTransform, nfeatures = 500)
seurat_split_SCTransform <- PrepSCTIntegration(seurat_split_SCTransform, anchor.features = features.sct)
seurat_split_SCTransform <- lapply(seurat_split_SCTransform, FUN = RunPCA, features = features.sct)

anchors.sct <- FindAnchors.STACAS(ref.list, anchor.features=features.sct, 
                                  normalization.method = "SCT")

anchors.sct.filtered <- FilterAnchors.STACAS(anchors.sct)

# Run seurat integrate on SCT data
seurat_integrated_SCTransform <- IntegrateData(anchorset=anchors.sct.filtered, dims=1:30, normalization.method = "SCT", preserve.order=T)

# Save RDS after integration
saveRDS(seurat_integrated_SCTransform, paste0(rds_path, "seurat_integrated_SCTransform.RDS"))






########## run integrate on scale data rather than sctransform

# Log normalize data and find variable features
seurat_split_scale <- lapply(seurat_split, function(x) NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000))
seurat_split_scale <- lapply(seurat_split_scale, function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000))

seurat_integrated_scale <- Run.STACAS(seurat_split_scale, dims=1:30, anchor.features=2000, plot.file = paste0(plot_path, 'integrated_scale.png'))
seurat_integrated_scale <- IntegrateData(anchorset = seurat_integrated_scale, dims = 1:30)


# Scale data and regress out MT content
# Enable parallelisation
# Multi-core when running from command line
if(opt$runtype == "nextflow"){
  plan("multiprocess", workers = ncores)
  options(future.globals.maxSize = 32* 1024^3) # 32gb
}

seurat_integrated_scale <- ScaleData(seurat_integrated_scale, features = rownames(seurat_integrated_scale), vars.to.regress = "percent.mt")

# Save RDS after scaling as this step takes time
saveRDS(seurat_integrated_scale, paste0(rds_path, "seurat_integrated_scale.RDS"))


